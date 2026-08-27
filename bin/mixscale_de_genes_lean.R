# Lean re-implementation of the DE-gene selection that RunMixscale performs via
# Seurat:::TopDEGenesMixscape -> FindMarkers(test.use = "wilcox", min.pct = 0.1).
#
# Returns exactly what mixscale_de_genes() returns: a named list of character
# vectors, one per target, in FindMarkers' row order. See de_genes_lean.R for
# the memory argument and the exactness argument.
#
# What is replaced and why (all measured at 305,110 cells):
#   * FoldChange           -- 26.2 GB of transients per call, recomputed per
#                             guide over an NT pool that never changes.
#   * LayerData()          -- Seurat v5 stores layers without dimnames and
#                             attaches them on fetch, duplicating the whole
#                             15.8 GB matrix. The raw layer is used instead and
#                             indexed positionally.
#   * WilcoxDETest         -- re-subsets data.use to c(cells.1, cells.2) twice
#                             more, both times a no-op reorder of the column
#                             order PerformDE already produced. presto is
#                             called directly on PerformDE's matrix.
#
# Everything that decides the returned gene set is reproduced verbatim: the
# min.pct / logfc gates, presto's p-values, the order(p_val, -abs(dpct)) sort,
# and the bonferroni denominator n = nrow(object) (the FULL gene count, not the
# filtered one).

# de_genes_lean.R (the blocked row-stat helpers) must be sourced first; both
# files are loaded side by side by 02_mixscale_preprocess.R.

# Seurat v5 keeps dimnames off the stored layer; fetch it directly when the
# layer covers every cell in the assay, else fall back to the copying accessor.
raw_data_layer <- function(assay, layer = "data") {
  if (methods::.hasSlot(assay, "layers")) {
    m <- assay@layers[[layer]]
    if (!is.null(m) &&
        nrow(m) == length(SeuratObject::Features(assay)) &&
        ncol(m) == length(SeuratObject::Cells(assay))) {
      return(m)
    }
  }
  SeuratObject::LayerData(assay, layer = layer)
}

# presto::wilcoxauc peaks at ~6x the size of the matrix handed to it (measured:
# 96.2 GB peak on a 15.1 GB input at 305,110 cells). The Wilcoxon rank-sum
# statistic is computed independently per gene row, so splitting the feature set
# into row blocks gives bit-identical p-values at a peak proportional to the
# block. Blocks are cut by nonzero count, not gene count, because per-gene
# density spans orders of magnitude.
chunk_by_nnz <- function(nnz_per_row, idx, target_nnz) {
  out <- list(); a <- 1L; n <- length(idx)
  while (a <= n) {
    cs <- cumsum(as.numeric(nnz_per_row[idx[a:n]]))
    take <- max(1L, sum(cs <= target_nnz))
    b <- a + take - 1L
    out[[length(out) + 1L]] <- idx[a:b]
    a <- b + 1L
  }
  out
}

mixscale_de_genes_lean <- function(object,
                                   labels,
                                   nt.class.name,
                                   de.assay = "RNA",
                                   logfc.threshold = 0,
                                   fine.mode = FALSE,
                                   fine.mode.labels = NULL,
                                   pval.cutoff = 5e-2,
                                   min.pct = 0.1,
                                   split.by = NULL,
                                   harmonize = FALSE,
                                   verbose = TRUE,
                                   block = 20000L,
                                   target_nnz = 1.2e8,
                                   outer_nnz = NULL,
                                   threads = 1L) {
  # The outer block is what gets pulled out of d; the inner chunk is what gets
  # handed to presto. Making the outer block bigger buys fewer full scans of d
  # (~5 s each at 305k) at the cost of holding a bigger block, and its column
  # subset, alongside presto's ~6x working set. Scans are already amortised
  # across every guide, so the default trades the seconds for the memory and
  # keeps the two equal.
  outer_nnz <- if (is.null(outer_nnz) || !is.finite(outer_nnz) || outer_nnz <= 0)
    target_nnz else outer_nnz
  # Every gene row is an independent rank-sum, so the inner chunks of one
  # comparison can be fitted concurrently and reassembled in chunk order for a
  # bit-identical result. presto peaks at ~6x the chunk it is handed, so the
  # chunk is divided by the worker count: the aggregate working set stays what
  # the serial path held, and only the wall-clock moves. The workers fork, and
  # everything they read -- d, the outer block, its column subset -- is already
  # built and never written, so copy-on-write keeps them sharing it.
  nthreads  <- max(1L, as.integer(threads))
  inner_nnz <- if (nthreads > 1L) target_nnz / nthreads else target_nnz
  if (!is.null(split.by) || isTRUE(harmonize)) {
    stop("mixscale_de_genes_lean() does not replicate the split.by/harmonize ",
         "path; call RunMixscale directly if you need it.")
  }
  assay <- object[[de.assay]]
  d     <- raw_data_layer(assay, "data")
  genes <- as.character(SeuratObject::Features(assay))
  n_genes_total <- nrow(d)
  n_cells_total <- ncol(d)
  # `d@i + 1L` would allocate a second copy of the index vector -- 5.7 GB at
  # 305k cells -- so the shift is done a slice at a time. tabulate() counts, and
  # integer addition of counts is exact, so the result is identical.
  row_nnz <- integer(n_genes_total)
  .nz <- length(d@i)
  .a  <- 1L
  while (.a <= .nz) {
    .b <- min(.nz, .a + 5e7L - 1L)
    row_nnz <- row_nnz + tabulate(d@i[.a:.b] + 1L, nbins = n_genes_total)
    .a <- .b + 1L
  }
  rm(.nz, .a, .b)

  lab     <- object[[labels]][, 1]
  i2      <- which(lab == nt.class.name)           # NT columns, ascending
  targets <- setdiff(unique(lab), nt.class.name)

  # The NT half of FoldChange: same cells for every guide and every target, so
  # it is computed once here instead of once per FindMarkers call.
  if (verbose) message("[de.genes] NT reference stats over ", length(i2),
                       " cells (", ceiling(length(i2) / block), " blocks)")
  s2 <- nt_row_stats_blocked(d, i2, block = block)

  guide_lab <- if (isTRUE(fine.mode)) object[[fine.mode.labels]][, 1] else NULL

  # ---- pass 1: fold changes and the min.pct/logfc feature mask ---------------
  # Column subsets only, so this never scans the whole matrix.
  sets <- list()
  for (gene in targets) {
    i1 <- which(lab == gene)
    idx_list <- if (isTRUE(fine.mode)) {
      gds <- setdiff(unique(guide_lab[i1]), nt.class.name)
      lapply(gds, function(g) i1[guide_lab[i1] == g])
    } else {
      list(i1)
    }
    for (idx in idx_list) {
      s <- tryCatch({
        s1 <- row_stats_direct(d, idx)
        fc <- fc_from_stats(s1, s2)
        keep <- pmax(fc$pct.1, fc$pct.2) >= min.pct &
                abs(fc$avg_log2FC) >= logfc.threshold
        keep[is.na(keep)] <- FALSE
        fidx <- which(keep)
        cols <- sort(c(idx, i2))
        list(target = gene,
             cols   = cols,
             y      = factor(ifelse(lab[cols] == nt.class.name, "Group2", "Group1"),
                             levels = c("Group1", "Group2")),
             fidx   = fidx,
             pct.1  = fc$pct.1[fidx],
             pct.2  = fc$pct.2[fidx],
             pval   = numeric(length(fidx)),
             failed = FALSE)
      }, error = function(e) NULL)
      if (!is.null(s) && length(s$fidx)) sets[[length(sets) + 1L]] <- s
    }
  }
  if (verbose) message("[de.genes] ", length(sets), " comparison(s) over ",
                       length(targets), " target(s)")

  # ---- pass 2: p-values, blocked over genes ---------------------------------
  # d[rows, ] is a dgCMatrix ROW subset, which scans every nonzero -- so it is
  # done once per outer block over the UNION of every comparison's features,
  # not once per block per comparison. presto is order-invariant (verified
  # bit-identical under both the natural order and a random permutation), so
  # the block is taken in natural column order and each comparison takes the
  # columns it needs from it, which is a cheap column subset.
  if (length(sets)) {
    all_feat <- sort(unique(unlist(lapply(sets, `[[`, "fidx"))))
    outer_blocks <- chunk_by_nnz(row_nnz, all_feat, outer_nnz)
    if (verbose) message("[de.genes] ", length(all_feat), " features pass the gates (",
                         sprintf("%.1f%%", 100 * sum(as.numeric(row_nnz[all_feat])) /
                                 max(1, sum(as.numeric(row_nnz)))),
                         " of nonzeros), ", length(outer_blocks), " block(s)")

    for (ob in outer_blocks) {
      du <- d[ob, , drop = FALSE]
      ob_nnz <- row_nnz[ob]
      for (k in seq_along(sets)) {
        s <- sets[[k]]
        if (isTRUE(s$failed)) next
        sel <- which(s$fidx >= ob[1L] & s$fidx <= ob[length(ob)])
        if (!length(sel)) next
        sub_cols <- length(s$cols) != n_cells_total
        dus <- if (sub_cols) du[, s$cols, drop = FALSE] else du
        ok <- tryCatch({
          rows_local <- match(s$fidx[sel], ob)
          chunks <- chunk_by_nnz(ob_nnz, rows_local, inner_nnz)
          one <- function(ib) {
            m <- dus[ib, , drop = FALSE]
            rownames(m) <- genes[ob[ib]]
            res <- presto::wilcoxauc(X = m, y = s$y)
            res <- res[seq_len(nrow(res) / 2), ]
            stopifnot(identical(as.character(res$feature), genes[ob[ib]]))
            res$pval
          }
          pv <- if (nthreads > 1L && length(chunks) > 1L) {
            parallel::mclapply(chunks, one,
                               mc.cores = min(nthreads, length(chunks)))
          } else {
            lapply(chunks, one)
          }
          # mclapply reports a worker failure as a try-error in the result list
          # rather than raising, so the shape is checked before it is trusted.
          stopifnot(all(vapply(pv, is.numeric, logical(1))),
                    identical(lengths(pv), lengths(chunks)))
          at <- 0L
          for (j in seq_along(chunks)) {
            sets[[k]]$pval[sel[at + seq_along(chunks[[j]])]] <- pv[[j]]
            at <- at + length(chunks[[j]])
          }
          rm(pv)
          TRUE
        }, error = function(e) FALSE)
        if (!ok) sets[[k]]$failed <- TRUE
        # dus is a column subset of du -- 1.4 GB at 305k -- and there is one per
        # comparison per outer block. Without a collect here R holds every one of
        # them until the heap reaches its growth trigger, which is what pinned
        # this stage's peak at the trigger rather than at its live set.
        if (sub_cols) { rm(dus); gc(FALSE) }
      }
      rm(du, ob_nnz); gc(FALSE)
    }
  }

  # ---- the tail of FindMarkers.default --------------------------------------
  out <- setNames(vector("list", length(targets)), targets)
  for (g in targets) out[[g]] <- character(0)
  for (s in sets) {
    if (isTRUE(s$failed)) next
    de <- data.frame(p_val = s$pval, pct.1 = s$pct.1, pct.2 = s$pct.2,
                     row.names = genes[s$fidx])
    de <- de[order(de$p_val, -abs(de$pct.1 - de$pct.2)), ]
    de$p_val_adj <- p.adjust(de$p_val, method = "bonferroni", n = n_genes_total)
    # fine.mode unions across guides; the single-set path is already unique.
    out[[s$target]] <- unique(c(out[[s$target]],
                                rownames(de)[de$p_val_adj < pval.cutoff]))
  }
  out
}
