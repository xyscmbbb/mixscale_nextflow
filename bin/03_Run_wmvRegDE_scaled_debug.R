#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(Matrix)
  library(Seurat)
  library(Mixscale)
  library(glmGamPoi)
  library(dplyr)
})

option_list <- list(
  make_option("--obj_rds", type = "character"),
  make_option("--perturb_gene", type = "character"),
  make_option("--target_gene_col", type = "character", default = "target_gene"),
  make_option("--nt_label", type = "character", default = "ONE_INTERGENIC_SITE"),
  make_option("--logfc_threshold", type = "double", default = 0),
  make_option("--min_pct", type = "double", default = 0),
  make_option("--min_cells_group", type = "integer", default = 10),
  make_option("--subsample", type = "character", default = "true",
              help = "Whether to use subsample=TRUE in glm_gp. Default: true"),
  make_option("--save_mixscale_obj", type = "character", default = "false")
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

save_mixscale_obj <- tolower(opt$save_mixscale_obj) %in% c("true", "t", "1", "yes", "y")
`%||%` <- function(a, b) if (!is.null(a)) a else b

message("[03] Reading Mixscale object: ", opt$obj_rds)
obj <- readRDS(opt$obj_rds)

FoldChange_new <- function(obj, cells.1, cells.2, mean.fxn, fc.name, features) {
  fc <- Seurat::FoldChange(
    object = obj,
    cells.1 = cells.1,
    cells.2 = cells.2,
    features = features,
    mean.fxn = mean.fxn,
    fc.name = fc.name
  )
  m1 <- obj[features, cells.1, drop = FALSE]
  m2 <- obj[features, cells.2, drop = FALSE]
  fc$min.cell.1 <- Matrix::rowSums(m1 > 0)
  fc$min.cell.2 <- Matrix::rowSums(m2 > 0)
  fc
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
    pseudocount.use = 1,
    base = 2,
    full.results = FALSE,
    subsample = TRUE
) {
    if (!rlang::is_installed("glmGamPoi")) stop("Please install glmGamPoi.")
    if (slot != "counts") stop("The slot must be set to 'counts'.")

    assay <- assay %||% DefaultAssay(object = object)
    DefaultAssay(object) <- assay
    prtb_score <- Tool(object = object, slot = "RunMixscale")

    wt_PRTB_list <- sort(names(prtb_score))
    all_PRTB_list <- sort(unique(object[[labels]][, 1]))
    all_PRTB_list <- all_PRTB_list[all_PRTB_list != nt.class.name]
    wt_PRTB_list <- wt_PRTB_list[wt_PRTB_list %in% all_PRTB_list]

    if (!is.null(PRTB_list)) {
        wt_PRTB_list <- intersect(wt_PRTB_list, PRTB_list)
        all_PRTB_list <- intersect(all_PRTB_list, PRTB_list)
        if (length(all_PRTB_list) == 0) stop("No perturbation left for DE tests. Check PRTB_list.")
    }

    if (is.null(split.by)) {
        mat_B <- data.frame(
            cell_label = colnames(object),
            nCount_RNA = object[[total_ct_labels]][, 1],
            cell_type = as.factor(names(prtb_score[[1]])),
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

    extract_results <- function(dat, form, meta) {
        rownames(meta) <- meta$cell_label
        meta <- meta[colnames(dat), , drop = FALSE]
        fit <- try(glmGamPoi::glm_gp(data = dat, design = form, col_data = meta, size_factors = FALSE, on_disk = FALSE, subsample = subsample), silent = FALSE)
        if (inherits(fit, "try-error")) {
            fit <- glmGamPoi::glm_gp(data = as.matrix(dat), design = form, col_data = meta, size_factors = FALSE, on_disk = FALSE, subsample = subsample)
        }
        beta_mat <- as.matrix(fit$Beta)
        pred <- predict(fit, se.fit = TRUE, newdata = diag(ncol(beta_mat)))
        se_mat <- as.matrix(pred$se.fit)
        se_mat[!is.finite(se_mat) | se_mat == 0] <- NA_real_
        z2 <- (beta_mat / se_mat)^2
        p_mat <- pchisq(z2, df = 1, lower.tail = FALSE)
        res_df <- cbind(as.data.frame(beta_mat), as.data.frame(p_mat))
        colnames(res_df) <- c(paste0("beta_", colnames(beta_mat)), paste0("p_", colnames(beta_mat)))
        res_df$gene_ID <- rownames(res_df)
        res_df
    }

    safe_minmax <- function(x) {
        xr <- range(x, na.rm = TRUE)
        denom <- xr[2] - xr[1]
        if (!is.finite(denom) || denom <= 0) return(x)
        (x - xr[1]) / denom
    }
    has_variation <- function(w) {
        w <- suppressWarnings(as.numeric(w)); w <- w[is.finite(w)]
        any(w > 0) && any(w == 0)
    }

    all_res <- list()
    for (PRTB in all_PRTB_list) {
        prtb_start <- Sys.time()
        mat_A <- data.frame()

        if (PRTB %in% wt_PRTB_list) {
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
            celltype_list <- if (is.null(split.by)) "con1" else sort(unique(object[[split.by]][, 1]))
            mat_all <- mat_B[mat_B$gene %in% c(PRTB, nt.class.name), , drop = FALSE]
            mat_all$weight <- ifelse(mat_all$gene == PRTB, 1, 0)
            DE_FLAG <- "standard"
        }

        mat_all$log_ct <- log1p(mat_all$nCount_RNA)
        mat_all$log_ct[!is.finite(mat_all$log_ct)] <- 0
        idx_c <- match(mat_all$cell_label, colnames(object))
        count_data_sparse <- GetAssayData(object, assay = assay, layer = "counts")[, idx_c, drop = FALSE]

        fc_list <- list()
        for (idx_i in seq_along(celltype_list)) {
            ct <- if (is.null(split.by)) celltype_list[idx_i] else levels(mat_all$cell_type)[idx_i]
            idx_NT_ct <- match(mat_all$cell_label[mat_all$cell_type == ct & mat_all$gene == nt.class.name], colnames(object))
            idx_P_ct <- match(mat_all$cell_label[mat_all$cell_type == ct & mat_all$gene != nt.class.name], colnames(object))
            fc <- FoldChange_new(
                obj = GetAssayData(object, assay = assay, layer = "counts"),
                cells.1 = colnames(object)[idx_P_ct],
                cells.2 = colnames(object)[idx_NT_ct],
                mean.fxn = function(x) log((Matrix::rowSums(x) + pseudocount.use) / NCOL(x), base = base),
                fc.name = "avg_log2FC",
                features = rownames(object)
            )
            fc$status <- abs(fc$avg_log2FC) >= logfc.threshold &
                (fc$pct.1 >= min.pct | fc$pct.2 >= min.pct) &
                (fc$min.cell.1 >= min.cells.group | fc$min.cell.2 >= min.cells.group)
            fc$status2 <- (fc$pct.1 >= min.pct | fc$pct.2 >= min.pct) &
                (fc$min.cell.1 >= min.cells.group | fc$min.cell.2 >= min.cells.group)
            fc_list[[ct]] <- fc
        }

        for (idx_i in seq_along(fc_list)) {
            ct_name <- names(fc_list)[idx_i]
            if (idx_i == 1) {
                idx_list <- data.frame(ct_name = fc_list[[ct_name]]$status); colnames(idx_list) <- ct_name
                idx_list2 <- data.frame(ct_name = fc_list[[ct_name]]$status2); colnames(idx_list2) <- ct_name
                fc_mat <- data.frame(ct_name = fc_list[[ct_name]]$avg_log2FC); colnames(fc_mat) <- paste0("log2FC_", ct_name)
                rownames(idx_list) <- rownames(fc_list[[ct_name]])
                rownames(idx_list2) <- rownames(fc_list[[ct_name]])
                rownames(fc_mat) <- rownames(fc_list[[ct_name]])
            } else {
                idx_list[[ct_name]] <- fc_list[[ct_name]]$status
                idx_list2[[ct_name]] <- fc_list[[ct_name]]$status2
                fc_mat[[paste0("log2FC_", ct_name)]] <- fc_list[[ct_name]]$avg_log2FC
            }
        }

        idx_for_DE <- which(apply(idx_list, 1, any))
        idx_PRTB_id <- which(rownames(object) %in% PRTB)
        idx_for_DE <- unique(c(idx_PRTB_id, idx_for_DE))
        fc_mat[!as.matrix(idx_list2)] <- NA

        form_call <- if (length(celltype_list) > 1) ~ 0 + cell_type + weight:cell_type + log_ct else ~ 1 + weight + log_ct

        if (DE_FLAG == "weighted") {
            loo_targets <- gsub("^weight_", "", grep("^weight_", colnames(mat_all), value = TRUE))
            if (PRTB %in% rownames(object)) loo_targets <- unique(c(PRTB, loo_targets))
            keep_loo <- vapply(loo_targets, function(g) {
                if (g == PRTB) {
                    w <- as.numeric(mat_all$gene != nt.class.name)
                } else {
                    wcol <- paste0("weight_", g)
                    if (!wcol %in% colnames(mat_all)) return(FALSE)
                    w <- mat_all[[wcol]]
                }
                has_variation(w)
            }, logical(1))
            loo_targets <- loo_targets[keep_loo]
            idx_loo_m <- which(rownames(object) %in% loo_targets)
            idx_std_m <- setdiff(idx_for_DE, idx_loo_m)
            if (length(idx_std_m) > 0) {
                res <- extract_results(count_data_sparse[idx_std_m, , drop = FALSE], form_call, mat_all)
            } else {
                res <- data.frame(gene_ID = character(0))
            }
            for (g_target in loo_targets) {
                idx_target <- which(rownames(count_data_sparse) == g_target)
                if (length(idx_target) == 0) next
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
                if (!has_variation(mat_loo$weight)) next
                tmp_res <- extract_results(count_data_sparse[idx_target, , drop = FALSE], form_call, mat_loo)
                res <- rbind(res, tmp_res)
            }
        } else {
            res <- extract_results(count_data_sparse[idx_for_DE, , drop = FALSE], form_call, mat_all)
        }

        fc_mat$gene_ID <- rownames(fc_mat)
        res <- merge(res, fc_mat, by = "gene_ID", all.x = TRUE)
        res$DE_method <- DE_FLAG
        if (!full.results) {
            b_cols <- grep("^beta_.*weight", names(res), value = TRUE)
            p_cols <- grep("^p_.*weight", names(res), value = TRUE)
            f_cols <- grep("^log2FC", names(res), value = TRUE)
            res <- res[, c("gene_ID", f_cols, b_cols, p_cols, "DE_method"), drop = FALSE]
            if (length(f_cols) == 1 && !"log2FC" %in% colnames(res)) colnames(res)[colnames(res) == f_cols] <- "log2FC"
        }
        all_res[[PRTB]] <- res
        if (verbose) {
            prtb_time <- round(as.numeric(difftime(Sys.time(), prtb_start, units = "mins")), 2)
            message(sprintf("[03] DE for '%s' done in %s min. Remaining = %d", PRTB, prtb_time, length(all_PRTB_list) - match(PRTB, all_PRTB_list)))
        }
    }
    all_res
}

message("[03] Running weighted Mixscale DE only for: ", opt$perturb_gene)
de_res <- Run_wmvRegDE_scaled_debug(
  object = obj,
  assay = "RNA",
  slot = "counts",
  total_ct_labels = "nCount_RNA",
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
if (is.null(de_res_df)) de_res_df <- de_res[[1]]
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
write.csv(de_res_df, paste0(prefix, "_de_res_df.csv"), row.names = FALSE, quote = FALSE)
write.csv(obj@meta.data, paste0(prefix, "_meta_data.csv"), row.names = TRUE, quote = FALSE)
message("[03] Done.")
