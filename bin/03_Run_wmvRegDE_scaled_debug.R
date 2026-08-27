#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(Matrix)
  library(Seurat)
  library(Mixscale)
  library(glmGamPoi)
  library(dplyr)
})

print(sessionInfo())

option_list <- list(
  make_option("--obj_rds", type = "character"),
  make_option("--perturb_gene", type = "character"),
  make_option("--target_gene_col", type = "character", default = "target_gene"),
  make_option("--nt_label", type = "character", default = "ONE_INTERGENIC_SITE"),
  make_option("--logfc_threshold", type = "double", default = 0),
  make_option("--min_pct", type = "double", default = 0),
  make_option("--min_cells_group", type = "integer", default = 10),
  make_option("--subsample", type = "character", default = "false",
              help = "Whether to use subsample=TRUE in glm_gp. Default: false"),
  make_option("--save_mixscale_obj", type = "character", default = "false"),
  make_option("--collapsed", type = "character", default = "false",
              help = paste("Fit with the collapsed Gamma-Poisson solver instead of stock",
                           "glm_gp. NT cells carry weight == 0 exactly, so they differ only",
                           "in log_ct = log1p(nCount) and cells sharing a design row share",
                           "mu; every per-cell sum in the fit therefore reduces to a",
                           "per-group sum. The reduction is exact. Stock glm_gp's dense",
                           "n_genes x n_cells Mu becomes n_genes x n_groups, which is what",
                           "makes multi-million-cell fits possible at all. Default: false")),
  make_option("--fc_norm", type = "character", default = "log.norm",
              help = paste("How avg_log2FC is computed: 'log.norm' (default) uses the",
                           "library-size-normalised expression, matching Seurat's",
                           "LogNormalize convention; 'raw' reproduces upstream Mixscale's",
                           "un-normalised counts behaviour. See NOTE in compute_fc_stats().")),
  make_option("--threads", type = "integer", default = 1,
              help = paste("Worker threads for the --collapsed solver only (OpenMP over",
                           "genes in the IRLS/SE passes, forked workers for the",
                           "overdispersion MLE). Genes are independent, so this does not",
                           "change results. Ignored on the stock glm_gp path. Default: 1")),
  make_option("--scale_factor", type = "double", default = 1e4,
              help = "Scale factor used by NormalizeData() in step 02. Default: 1e4")
)

opt <- parse_args(OptionParser(option_list = option_list))
stopifnot(!is.null(opt$obj_rds), !is.null(opt$perturb_gene))

to_bool <- function(x) {
  x <- tolower(as.character(x))
  if (x %in% c("true", "t", "1", "yes", "y")) return(TRUE)
  if (x %in% c("false", "f", "0", "no", "n")) return(FALSE)
  stop("Cannot parse logical value: ", x)
}

subsample_use <- to_bool(opt$subsample)
message("[03] glm_gp subsample: ", subsample_use)

use_collapsed <- tolower(opt$collapsed) %in% c("true", "t", "yes", "1")
if (use_collapsed) {
  .self3 <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
  source(file.path(dirname(.self3), "collapsed", "load_collapsed.R"))
  load_collapsed_glm_gp(file.path(dirname(.self3), "collapsed"))
  if (isTRUE(subsample_use)) {
    stop("[03] --collapsed does not implement glm_gp's subsample path; use --subsample false.")
  }
}
message("[03] collapsed solver: ", use_collapsed)

save_mixscale_obj <- tolower(opt$save_mixscale_obj) %in% c("true", "t", "1", "yes", "y")

`%||%` <- function(a, b) {
  if (!is.null(a)) a else b
}

message("[03] Reading Mixscale object: ", opt$obj_rds)
obj <- readRDS(opt$obj_rds)

# Per-gene log2FC + DE-filter flags, computed in a single pass over the already
# cell-subset sparse counts matrix.
#
# The previous FoldChange_new() called Seurat::FoldChange() (which internally
# subsets object[features, cells.1/2] and computes rowSums / rowSums(x>0)) and
# then re-subset the FULL object matrix two MORE times to recompute min.cell.*.
# For ~60k features that is four full-feature subsets plus two `m > 0` sparse
# allocations per call. Here we instead take two column slices of the small,
# already-subset count matrix and derive everything from rowSums + the stored
# nonzero counts (`@i` via tabulate), so expressing-cell counts are exact
# integers rather than rounded pct * n. Identical avg_log2FC, lower peak RAM.
# NOTE on avg_log2FC (changed deliberately; diverges from upstream Mixscale).
#
# Upstream Mixscale (satijalab/Mixscale, R/scoring_de.R) computes
#     mean.fxn = function(x) log((rowSums(x) + pseudocount.use) / NCOL(x), base)
# on the RAW COUNTS layer. That has two consequences we do not want:
#
#   (a) No library-size normalisation. The DE *test* already controls for depth
#       (log_ct is a covariate in the glm_gp design, hence size_factors = FALSE),
#       but the fold change never sees it, so avg_log2FC tracks sequencing depth.
#   (b) The pseudocount is added to the SUM, so the term is rowMeans + 1/N. With
#       unequal group sizes a gene absent from both groups gets
#       log2(N_nt / N_pert) instead of 0 -- e.g. +5.7 for 319 perturbed vs 16,297
#       NT cells. The bias is largest for lowly expressed genes and for
#       perturbations with few cells, which is exactly where it does most damage.
#
# Upstream's own get_fc() (R/get_fold_change.R) uses the correct convention,
#     log(mean(x) + pseudocount.use)  /  log(mean(expm1(x)) + pseudocount.use)
# and exposes norm.method = 'log.norm'; Run_wmvRegDE() simply never calls it.
#
# We therefore default to fc_norm = "log.norm":
#     avg_log2FC = log2(mean(counts / nCount_RNA * scale.factor) + pseudocount)
# which is identical to Seurat's FoldChange() on the LogNormalize "data" layer,
# because expm1(log1p(y)) == y. Step 02 slims the object to the counts layer, so
# we reconstruct the normalisation here from counts + nCount_RNA rather than
# carrying a second matrix through the pipeline.
#
# Pass --fc_norm raw to reproduce the previous (upstream) behaviour exactly.
# A column subset followed by rowSums is a sparse matrix-vector product. The
# subset only zeroes the columns that were not selected and the Diagonal only
# rescales the ones that were, and a weight vector that is zero outside `cols`
# already says both:
#
#   rowSums(m[, cols, drop = FALSE] %*% Matrix::Diagonal(x = w))  ==  m %*% v
#
# with v[cols] <- w. That allocates one length-ncell weight vector and one
# length-ngene result instead of two matrix-sized copies -- up to 4 x 16 GB at
# 305,110 cells, which is the bulk of step 3's peak. Verified BIT-identical at
# CM scale for both the NT and the perturbed group (max|diff| 0, 0 differing
# elements, and the log2 term identical). A blocked accumulation is NOT
# bit-identical (max|diff| 4.4e-9): it changes the order of the summation,
# whereas the matvec keeps Matrix's own column-by-column order.
row_sums_scaled <- function(m, cols, w) {
  v <- numeric(ncol(m))
  v[cols] <- as.numeric(w)
  as.numeric(m %*% v)
}

# Nonzeros per row over a subset of columns, without subsetting. A dgCMatrix
# column j stores its row indices at @i[(@p[j]+1):@p[j+1]], so the selected
# columns' entries are gathered directly and tabulated. Counting is exact
# integer arithmetic, so how the columns are grouped cannot change the result;
# the grouping only bounds the largest index vector held at once.
n_expr_cols <- function(m, cols, chunk = 2e7) {
  p   <- m@p
  len <- p[cols + 1L] - p[cols]
  out <- integer(nrow(m))
  a <- 1L; n <- length(cols)
  while (a <= n) {
    cs <- cumsum(as.numeric(len[a:n]))
    b  <- a + max(1L, sum(cs <= chunk)) - 1L
    idx <- sequence(len[a:b], from = p[cols[a:b]] + 1L)
    out <- out + tabulate(m@i[idx] + 1L, nrow(m))
    rm(idx)
    a <- b + 1L
  }
  out
}

compute_fc_stats <- function(counts_mat, p_cols, nt_cols,
                             total_counts = NULL,
                             norm.method = "log.norm",
                             scale.factor = 1e4,
                             pseudocount.use = 1, base = 2,
                             logfc.threshold = 0, min.pct = 0,
                             min.cells.group = 10,
                             feats = NULL) {
  if (is.null(feats)) feats <- rownames(counts_mat)
  nfeat <- length(feats)

  # Cells expressing each gene = stored nonzeros per row (counts are > 0 where
  # stored), read off the (p, i) slots so no column subset is ever built.
  n_expr_P <- n_expr_cols(counts_mat, p_cols)
  n_expr_N <- n_expr_cols(counts_mat, nt_cols)

  # n_expr_* are computed above from the RAW counts; normalisation cannot change
  # which entries are non-zero, so the min.pct / min.cells filters are unaffected.
  if (identical(norm.method, "log.norm")) {
    if (is.null(total_counts)) {
      stop("compute_fc_stats(): norm.method='log.norm' requires total_counts.")
    }
    tp <- as.numeric(total_counts[p_cols])
    tn <- as.numeric(total_counts[nt_cols])
    bad_p <- !is.finite(tp) | tp <= 0
    bad_n <- !is.finite(tn) | tn <= 0
    # Only the offending columns are summed; colSums over the full P/N was only
    # ever read at these positions.
    if (any(bad_p))
      tp[bad_p] <- Matrix::colSums(counts_mat[, p_cols[bad_p], drop = FALSE])
    if (any(bad_n))
      tn[bad_n] <- Matrix::colSums(counts_mat[, nt_cols[bad_n], drop = FALSE])
    tp[!is.finite(tp) | tp <= 0] <- 1
    tn[!is.finite(tn) | tn <= 0] <- 1

    # expm1(LogNormalize(counts)) == counts / nCount_RNA * scale.factor, exactly.
    mean_P <- log(row_sums_scaled(counts_mat, p_cols,  scale.factor / tp) /
                    length(p_cols)  + pseudocount.use, base = base)
    mean_N <- log(row_sums_scaled(counts_mat, nt_cols, scale.factor / tn) /
                    length(nt_cols) + pseudocount.use, base = base)
  } else if (identical(norm.method, "raw")) {
    mean_P <- log((row_sums_scaled(counts_mat, p_cols,  rep(1, length(p_cols))) +
                     pseudocount.use) / length(p_cols),  base = base)
    mean_N <- log((row_sums_scaled(counts_mat, nt_cols, rep(1, length(nt_cols))) +
                     pseudocount.use) / length(nt_cols), base = base)
  } else {
    stop("compute_fc_stats(): fc_norm must be 'log.norm' or 'raw', got: ", norm.method)
  }
  avg_log2FC <- mean_P - mean_N

  pass_pct  <- (n_expr_P / length(p_cols)  >= min.pct) |
               (n_expr_N / length(nt_cols) >= min.pct)
  pass_cell <- (n_expr_P >= min.cells.group) | (n_expr_N >= min.cells.group)

  data.frame(
    avg_log2FC = avg_log2FC,
    status     = (abs(avg_log2FC) >= logfc.threshold) & pass_pct & pass_cell,
    status2    = pass_pct & pass_cell,
    row.names  = feats
  )
}

Run_wmvRegDE_scaled_debug <- function(
    object,
    assay = "RNA",
    slot = "counts",
    labels = "gene",
    nt.class.name = "NT",
    verbose = FALSE,
    PRTB_list = NULL,
    split.by = NULL,
    logfc.threshold = 0,
    min.pct = 0.1,
    min.cells.group = 10,
    total_ct_labels = "nCount_RNA",
    fc.norm.method = "log.norm",
    scale.factor = 1e4,
    pseudocount.use = 1,
    base = 2,
    full.results = FALSE,
    subsample = TRUE
) {
  if (!rlang::is_installed("glmGamPoi")) {
    stop("Please install glmGamPoi.")
  }

  if (slot != "counts") {
    stop("The slot must be set to 'counts'.")
  }

  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object) <- assay

  prtb_score <- Tool(object = object, slot = "RunMixscale")

  # ------------------------------------------------------------
  # Important fix:
  # If Mixscale found zero DE genes, RunMixscale can be empty.
  # The old script crashed at prtb_score[[1]] before reaching
  # the binary 1/0 fallback branch.
  # ------------------------------------------------------------
  if (is.null(prtb_score) || length(prtb_score) == 0) {
    message("[03] RunMixscale tool is empty or missing.")
    message("[03] This often happens when Mixscale found 0 DE genes.")
    message("[03] Falling back to standard binary weights: perturbed cells = 1, NT cells = 0.")
    wt_PRTB_list <- character(0)
  } else {
    wt_PRTB_list <- sort(names(prtb_score))
  }

  all_PRTB_list <- sort(unique(object[[labels]][, 1]))
  all_PRTB_list <- all_PRTB_list[all_PRTB_list != nt.class.name]

  wt_PRTB_list <- wt_PRTB_list[wt_PRTB_list %in% all_PRTB_list]

  if (!is.null(PRTB_list)) {
    wt_PRTB_list <- intersect(wt_PRTB_list, PRTB_list)
    all_PRTB_list <- intersect(all_PRTB_list, PRTB_list)

    if (length(all_PRTB_list) == 0) {
      stop("No perturbation left for DE tests. Check PRTB_list, labels, and nt.class.name.")
    }
  }

  if (is.null(split.by)) {
    mat_B <- data.frame(
      cell_label = colnames(object),
      nCount_RNA = object[[total_ct_labels]][, 1],
      cell_type = as.factor(rep("con1", ncol(object))),
      gene = object[[labels]][, 1]
    )
  } else {
    mat_B <- data.frame(
      cell_label = colnames(object),
      nCount_RNA = object[[total_ct_labels]][, 1],
      cell_type = as.factor(object[[split.by]][, 1]),
      gene = object[[labels]][, 1]
    )
  }

  # `gene_idx` selects rows of `dat` to fit. Passing indices rather than
  # `dat[gene_idx, ]` matters at scale: with min_pct = 0 the selection is ~99.7%
  # of rows, so the subset is a near-complete second copy of the counts matrix.
  # The collapsed solver slices its gene chunks out of the transpose anyway, so
  # it can take the index directly and never materialise the subset.
  # dat_genes / dat_cells carry the dimnames that `dat` deliberately does not.
  # They default to dat's own, so the small leave-one-out blocks (which are named)
  # can be passed straight through.
  extract_results <- function(dat, form, meta, gene_idx = NULL,
                              dat_genes = NULL, dat_cells = NULL) {
    if (is.null(dat_genes)) dat_genes <- rownames(dat)
    if (is.null(dat_cells)) dat_cells <- colnames(dat)
    # mat_all is built in the object's column order (see the idx_c block below),
    # so this realignment is normally the identity. Doing it anyway copies a
    # ~100-column data frame -- one column per Mixscale DE gene -- and matches
    # every cell name, once per leave-one-out fit. Check first.
    rownames(meta) <- meta$cell_label
    if (!identical(rownames(meta), dat_cells))
      meta <- meta[dat_cells, , drop = FALSE]

    if (use_collapsed) {
      # Build the model matrix through glm_gp's own handler so the columns,
      # their order and their names are identical to the stock path.
      mm <- glmGamPoi:::handle_design_parameter(form, dat, meta, NULL)$model_matrix
      cf <- collapsed_glm_gp(dat, mm, gene_idx = gene_idx, gnames = dat_genes,
                             verbose = TRUE,
                             threads = max(1L, as.integer(opt$threads)),
                             solver = Sys.getenv("MIXSCALE_SOLVER", "chol"))
      beta_mat <- cf$Beta
      se_mat <- cf$se
      se_mat[!is.finite(se_mat) | se_mat == 0] <- NA_real_
      z2 <- (beta_mat / se_mat)^2
      p_mat <- pchisq(z2, df = 1, lower.tail = FALSE)
      res_df <- cbind(as.data.frame(beta_mat), as.data.frame(p_mat))
      colnames(res_df) <- c(paste0("beta_", colnames(beta_mat)),
                            paste0("p_", colnames(beta_mat)))
      res_df$gene_ID <- rownames(res_df)
      return(res_df)
    }

    if (!is.null(gene_idx)) {
      dat <- dat[gene_idx, , drop = FALSE]
      dimnames(dat) <- list(dat_genes[gene_idx], dat_cells)
    } else if (is.null(dimnames(dat)[[1]])) {
      dimnames(dat) <- list(dat_genes, dat_cells)
    }

    fit <- try(
      glmGamPoi::glm_gp(
        data = dat,
        design = form,
        col_data = meta,
        size_factors = FALSE,
        on_disk = FALSE,
        subsample = subsample
      ),
      silent = FALSE
    )

    if (inherits(fit, "try-error")) {
      message("[03] Sparse glm_gp failed; retrying with dense matrix.")
      fit <- glmGamPoi::glm_gp(
        data = as.matrix(dat),
        design = form,
        col_data = meta,
        size_factors = FALSE,
        on_disk = FALSE,
        subsample = subsample
      )
    }

    beta_mat <- as.matrix(fit$Beta)

    pred <- predict(fit, se.fit = TRUE, newdata = diag(ncol(beta_mat)))
    se_mat <- as.matrix(pred$se.fit)
    se_mat[!is.finite(se_mat) | se_mat == 0] <- NA_real_

    z2 <- (beta_mat / se_mat)^2
    p_mat <- pchisq(z2, df = 1, lower.tail = FALSE)

    res_df <- cbind(as.data.frame(beta_mat), as.data.frame(p_mat))
    colnames(res_df) <- c(
      paste0("beta_", colnames(beta_mat)),
      paste0("p_", colnames(beta_mat))
    )

    res_df$gene_ID <- rownames(res_df)
    res_df
  }

  safe_minmax <- function(x) {
    xr <- range(x, na.rm = TRUE)
    denom <- xr[2] - xr[1]

    if (!is.finite(denom) || denom <= 0) {
      return(x)
    }

    (x - xr[1]) / denom
  }

  has_variation <- function(w) {
    w <- suppressWarnings(as.numeric(w))
    w <- w[is.finite(w)]
    any(w > 0) && any(w == 0)
  }

  all_res <- list()

  for (PRTB in all_PRTB_list) {
    prtb_start <- Sys.time()
    mat_A <- data.frame()

    if (PRTB %in% wt_PRTB_list) {
      message("[03] ", PRTB, ": using weighted Mixscale perturbation scores.")

      celltype_list <- names(prtb_score[[PRTB]])

      for (celltype in celltype_list) {
        tmp <- prtb_score[[PRTB]][[celltype]]

        idx_NT <- which(tmp$gene == nt.class.name)
        idx_gene <- which(tmp$gene == PRTB)

        m_nt <- mean(tmp$pvec[idx_NT], na.rm = TRUE)
        s_nt <- sd(tmp$pvec[idx_NT], na.rm = TRUE)

        std_w <- (tmp$pvec[idx_gene] - m_nt) / s_nt
        std_w[!is.finite(std_w)] <- 0
        std_w[std_w < 0] <- 0

        tmp$weight <- 0
        tmp$weight[idx_gene] <- std_w
        tmp$weight <- safe_minmax(tmp$weight)

        if (ncol(tmp) >= 4) {
          for (idx_col in 3:(ncol(tmp) - 1)) {
            m_loo <- mean(tmp[idx_NT, idx_col], na.rm = TRUE)
            s_loo <- sd(tmp[idx_NT, idx_col], na.rm = TRUE)

            std_loo <- (tmp[idx_gene, idx_col] - m_loo) / s_loo
            std_loo[!is.finite(std_loo)] <- 0
            std_loo[std_loo < 0] <- 0

            wcol <- paste0("weight_", colnames(tmp)[idx_col])
            tmp[[wcol]] <- 0
            tmp[idx_gene, wcol] <- std_loo
            tmp[[wcol]] <- safe_minmax(tmp[[wcol]])
          }
        }

        tmp$cell_label <- rownames(tmp)
        keep_cols <- c("cell_label", "gene", grep("^weight", colnames(tmp), value = TRUE))
        mat_A <- rbind(mat_A, tmp[, keep_cols, drop = FALSE])
      }

      mat_all <- merge(mat_A, mat_B, by = c("cell_label", "gene"))
      DE_FLAG <- "weighted"

    } else {
      message("[03] ", PRTB, ": no valid Mixscale perturbation score found.")
      message("[03] ", PRTB, ": using standard binary weights: perturbed cells = 1, NT cells = 0.")

      celltype_list <- if (is.null(split.by)) {
        "con1"
      } else {
        sort(unique(object[[split.by]][, 1]))
      }

      mat_all <- mat_B[mat_B$gene %in% c(PRTB, nt.class.name), , drop = FALSE]
      mat_all$weight <- ifelse(mat_all$gene == PRTB, 1, 0)
      DE_FLAG <- "standard"
    }

    if (nrow(mat_all) == 0) {
      warning("[03] No cells found for ", PRTB, " versus ", nt.class.name, ". Writing empty result.")

      all_res[[PRTB]] <- data.frame(
        gene_ID = character(),
        log2FC = numeric(),
        beta_weight = numeric(),
        p_weight = numeric(),
        DE_method = character()
      )

      next
    }

    if (!has_variation(mat_all$weight)) {
      warning("[03] Weight has no variation for ", PRTB, ". Writing empty result.")

      all_res[[PRTB]] <- data.frame(
        gene_ID = character(),
        log2FC = numeric(),
        beta_weight = numeric(),
        p_weight = numeric(),
        DE_method = character()
      )

      next
    }

    mat_all$log_ct <- log1p(mat_all$nCount_RNA)
    mat_all$log_ct[!is.finite(mat_all$log_ct)] <- 0

    idx_c <- match(mat_all$cell_label, colnames(object))

    if (any(is.na(idx_c))) {
      warning("[03] Some cells in metadata were not found in object colnames. Removing them.")
      keep <- !is.na(idx_c)
      mat_all <- mat_all[keep, , drop = FALSE]
      idx_c <- idx_c[keep]
    }

    # Reordering the counts matrix to metadata order costs a full extra copy of
    # it -- 17 GB at 300k NT, and it is pure permutation, no data change. The GLM
    # is a sum over cells, so cell order is immaterial: reorder the metadata to
    # the object's column order instead. When the metadata covers every column
    # (the pipeline's case, since step 1 subsets the object to PRTB + NT) the
    # result is the identity and no subsetting happens at all.
    ord <- order(idx_c)
    mat_all <- mat_all[ord, , drop = FALSE]
    idx_c <- idx_c[ord]

    # Seurat v5 stores a layer WITHOUT dimnames and re-attaches them on fetch, so
    # GetAssayData() hands back a copy of the whole counts matrix -- measured at
    # 36.60 GiB peak (nnz = 1.430e9) against 20.44 GiB for taking the raw layer
    # and setting Dimnames on it, exactly one dgCMatrix (12.13 B/nnz) apart. The
    # two results are identical(), verified at CM scale, so the accessor is
    # bypassed and the names are attached directly.
    .a <- object[[assay]]
    count_data_sparse <- if (methods::.hasSlot(.a, "layers")) .a@layers[["counts"]] else NULL
    if (is.null(count_data_sparse))          # v3 Assay, or no counts layer
      count_data_sparse <- GetAssayData(object, assay = assay, layer = "counts")
    rm(.a)
    # Setting the dimnames costs a SECOND full copy: `object` still references the
    # layer, so R duplicates it to write the slot -- exactly the write-back
    # duplication that step 2 hit. Nothing here reads the names off the matrix, so
    # they are carried beside it and handed to the callees that do.
    cd_genes <- rownames(object)
    cd_cells <- colnames(object)
    if (!is.null(dimnames(count_data_sparse)[[1]]))
      cd_genes <- dimnames(count_data_sparse)[[1]]
    if (!is.null(dimnames(count_data_sparse)[[2]]))
      cd_cells <- dimnames(count_data_sparse)[[2]]
    if (!identical(idx_c, seq_len(ncol(count_data_sparse)))) {
      count_data_sparse <- count_data_sparse[, idx_c, drop = FALSE]
      cd_cells <- cd_cells[idx_c]
    }

    # count_data_sparse columns are in mat_all row order, so cell-type / NT vs
    # perturbed selections are plain column indices -- no re-match to the full
    # object and no second copy of the genome-wide counts matrix.
    fc_list <- list()

    for (idx_i in seq_along(celltype_list)) {
      ct <- if (is.null(split.by)) {
        celltype_list[idx_i]
      } else {
        levels(mat_all$cell_type)[idx_i]
      }

      ct_mask <- mat_all$cell_type == ct
      nt_cols <- which(ct_mask & mat_all$gene == nt.class.name)
      p_cols  <- which(ct_mask & mat_all$gene != nt.class.name)

      if (length(nt_cols) == 0 || length(p_cols) == 0) {
        warning("[03] Skipping FoldChange for ", PRTB, ", cell type ", ct,
                ": missing perturbed or NT cells.")
        next
      }

      fc_list[[ct]] <- compute_fc_stats(
        counts_mat      = count_data_sparse,
        feats           = cd_genes,
        p_cols          = p_cols,
        nt_cols         = nt_cols,
        total_counts    = mat_all$nCount_RNA,
        norm.method     = fc.norm.method,
        scale.factor    = scale.factor,
        pseudocount.use = pseudocount.use,
        base            = base,
        logfc.threshold = logfc.threshold,
        min.pct         = min.pct,
        min.cells.group = min.cells.group
      )
    }

    if (length(fc_list) == 0) {
      warning("[03] No valid FoldChange results for ", PRTB, ". Writing empty result.")

      all_res[[PRTB]] <- data.frame(
        gene_ID = character(),
        log2FC = numeric(),
        beta_weight = numeric(),
        p_weight = numeric(),
        DE_method = character()
      )

      next
    }

    for (idx_i in seq_along(fc_list)) {
      ct_name <- names(fc_list)[idx_i]

      if (idx_i == 1) {
        idx_list <- data.frame(ct_name = fc_list[[ct_name]]$status)
        colnames(idx_list) <- ct_name

        idx_list2 <- data.frame(ct_name = fc_list[[ct_name]]$status2)
        colnames(idx_list2) <- ct_name

        fc_mat <- data.frame(ct_name = fc_list[[ct_name]]$avg_log2FC)
        colnames(fc_mat) <- paste0("log2FC_", ct_name)

        rownames(idx_list) <- rownames(fc_list[[ct_name]])
        rownames(idx_list2) <- rownames(fc_list[[ct_name]])
        rownames(fc_mat) <- rownames(fc_list[[ct_name]])

      } else {
        idx_list[[ct_name]] <- fc_list[[ct_name]]$status
        idx_list2[[ct_name]] <- fc_list[[ct_name]]$status2
        fc_mat[[paste0("log2FC_", ct_name)]] <- fc_list[[ct_name]]$avg_log2FC
      }
    }

    # Row-wise OR across cell types, vectorised (apply(., 1, any) is ~60k slow
    # one-row calls; Reduce over the columns is a handful of vector ORs).
    idx_for_DE <- which(Reduce(`|`, idx_list))
    idx_PRTB_id <- which(rownames(object) %in% PRTB)
    idx_for_DE <- unique(c(idx_PRTB_id, idx_for_DE))

    if (length(idx_for_DE) == 0) {
      warning("[03] No genes passed DE filter for ", PRTB, ". Writing empty result.")

      all_res[[PRTB]] <- data.frame(
        gene_ID = character(),
        log2FC = numeric(),
        beta_weight = numeric(),
        p_weight = numeric(),
        DE_method = character()
      )

      next
    }

    fc_mat[!as.matrix(idx_list2)] <- NA

    form_call <- if (length(celltype_list) > 1) {
      ~ 0 + cell_type + weight:cell_type + log_ct
    } else {
      ~ 1 + weight + log_ct
    }

    if (DE_FLAG == "weighted") {
      loo_targets <- gsub("^weight_", "", grep("^weight_", colnames(mat_all), value = TRUE))

      if (PRTB %in% rownames(object)) {
        loo_targets <- unique(c(PRTB, loo_targets))
      }

      keep_loo <- vapply(
        loo_targets,
        function(g) {
          if (g == PRTB) {
            w <- as.numeric(mat_all$gene != nt.class.name)
          } else {
            wcol <- paste0("weight_", g)

            if (!wcol %in% colnames(mat_all)) {
              return(FALSE)
            }

            w <- mat_all[[wcol]]
          }

          has_variation(w)
        },
        logical(1)
      )

      loo_targets <- loo_targets[keep_loo]

      idx_loo_m <- which(rownames(object) %in% loo_targets)
      idx_std_m <- setdiff(idx_for_DE, idx_loo_m)

      # Main genome-wide fit.
      if (length(idx_std_m) > 0) {
        res <- extract_results(count_data_sparse, form_call, mat_all,
                               gene_idx = idx_std_m,
                               dat_genes = cd_genes, dat_cells = cd_cells)
      } else {
        res <- data.frame(gene_ID = character(0))
      }

      # Leave-one-out targets are single genes (cheap) but each uses a different
      # weight column, so they stay serial. Accumulate in a list and rbind once
      # instead of growing `res` quadratically inside the loop.
      loo_res <- vector("list", length(loo_targets))

      # Pull every leave-one-out gene out of the counts ONCE. count_data_sparse
      # is column-major (dgCMatrix), so subsetting a row scans all nonzeros --
      # taking one row costs the same as taking a hundred. Done inside the loop
      # that is ~5 s per fit at 450k cells, ~500 s over 100 fits, and it was the
      # single largest component of step 3. The extracted block is ~100 genes,
      # a few tens of MB, and row-subsetting it afterwards is free.
      loo_mat <- count_data_sparse[idx_loo_m, , drop = FALSE]
      dimnames(loo_mat) <- list(cd_genes[idx_loo_m], cd_cells)

      for (li in seq_along(loo_targets)) {
        g_target <- loo_targets[li]
        idx_target <- which(rownames(loo_mat) == g_target)

        if (length(idx_target) == 0) {
          next
        }

        mat_loo <- mat_all

        if (g_target == PRTB) {
          mat_loo$weight <- 0
          mat_loo$weight[mat_loo$gene != nt.class.name] <- 1
        } else {
          wcol <- paste0("weight_", g_target)
          mat_loo$weight <- suppressWarnings(as.numeric(mat_loo[[wcol]]))
          mat_loo$weight[!is.finite(mat_loo$weight)] <- 0
        }

        mat_loo$log_ct <- suppressWarnings(as.numeric(mat_loo$log_ct))
        mat_loo$log_ct[!is.finite(mat_loo$log_ct)] <- 0
        mat_loo$cell_type <- droplevels(mat_loo$cell_type)

        if (!has_variation(mat_loo$weight)) {
          next
        }

        loo_res[[li]] <- extract_results(
          loo_mat[idx_target, , drop = FALSE],
          form_call,
          mat_loo
        )
      }

      res <- do.call(rbind, c(list(res), loo_res))

    } else {
      res <- extract_results(count_data_sparse, form_call, mat_all,
                             gene_idx = idx_for_DE,
                             dat_genes = cd_genes, dat_cells = cd_cells)
    }

    fc_mat$gene_ID <- rownames(fc_mat)

    res <- merge(res, fc_mat, by = "gene_ID", all.x = TRUE)
    res$DE_method <- DE_FLAG

    if (!full.results) {
      b_cols <- grep("^beta_.*weight", names(res), value = TRUE)
      p_cols <- grep("^p_.*weight", names(res), value = TRUE)
      f_cols <- grep("^log2FC", names(res), value = TRUE)

      keep_cols <- c("gene_ID", f_cols, b_cols, p_cols, "DE_method")
      keep_cols <- keep_cols[keep_cols %in% colnames(res)]

      res <- res[, keep_cols, drop = FALSE]

      if (length(f_cols) == 1 && !"log2FC" %in% colnames(res)) {
        colnames(res)[colnames(res) == f_cols] <- "log2FC"
      }
    }

    all_res[[PRTB]] <- res

    if (verbose) {
      prtb_time <- round(as.numeric(difftime(Sys.time(), prtb_start, units = "mins")), 2)
      message(sprintf(
        "[03] DE for '%s' done in %s min. Remaining = %d",
        PRTB,
        prtb_time,
        length(all_PRTB_list) - match(PRTB, all_PRTB_list)
      ))
    }
  }

  all_res
}

message("[03] avg_log2FC norm.method: ", opt$fc_norm,
        " (scale.factor = ", opt$scale_factor, ")")
message("[03] Running weighted Mixscale DE only for: ", opt$perturb_gene)

de_res <- Run_wmvRegDE_scaled_debug(
  object = obj,
  assay = "RNA",
  slot = "counts",
  total_ct_labels = "nCount_RNA",
  fc.norm.method = opt$fc_norm,
  scale.factor = opt$scale_factor,
  labels = opt$target_gene_col,
  nt.class.name = opt$nt_label,
  PRTB_list = opt$perturb_gene,
  logfc.threshold = opt$logfc_threshold,
  min.pct = opt$min_pct,
  min.cells.group = opt$min_cells_group,
  verbose = TRUE,
  split.by = NULL,
  subsample = subsample_use
)

de_res_df <- de_res[[opt$perturb_gene]]

if (is.null(de_res_df)) {
  message("[03] No DE result returned for ", opt$perturb_gene, "; writing empty result.")
  de_res_df <- data.frame(
    gene_ID = character(),
    log2FC = numeric(),
    beta_weight = numeric(),
    p_weight = numeric(),
    DE_method = character()
  )
}

de_res_df$target_gene <- opt$perturb_gene

# Add BH-adjusted p value using the weighted p-value column produced by
# Run_wmvRegDE_scaled_debug. This avoids relying on a separate `pcols` object.
pval_cols <- grep("^p_.*weight|^p_weight$", colnames(de_res_df), value = TRUE)
message("[03] Weighted p-value columns found: ", paste(pval_cols, collapse = ", "))

if ("p_weight" %in% colnames(de_res_df)) {
  de_res_df$p_adj_bh <- p.adjust(de_res_df$p_weight, method = "BH")
} else if (length(pval_cols) >= 1) {
  de_res_df$p_adj_bh <- p.adjust(de_res_df[[pval_cols[1]]], method = "BH")
} else {
  warning("[03] No weighted p-value column found; setting p_adj_bh = NA.")
  de_res_df$p_adj_bh <- NA_real_
}

prefix <- paste0("pert_", opt$perturb_gene)

write.csv(
  de_res_df,
  paste0(prefix, "_de_res_df.csv"),
  row.names = FALSE,
  quote = FALSE
)

write.csv(
  obj@meta.data,
  paste0(prefix, "_meta_data.csv"),
  row.names = TRUE,
  quote = FALSE
)

if (save_mixscale_obj) {
  saveRDS(obj, paste0(prefix, "_mixscale_obj_after_wmvRegDE.rds"))
}

message("[03] Done.")