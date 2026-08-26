# Precompute the DE genes RunMixscale would select, before the PRTB assay exists.
#
# Why this exists: RunMixscale reads the PRTB assay in exactly one place (its body
# line 259) and immediately subsets it to de.genes:
#   dat <- GetAssayData(object[[assay]], layer = "data")[de.genes, all.cells]
# and de.genes comes from TopDEGenesMixscape(de.assay = "RNA"), which never touches
# PRTB. So step 2 builds a dense n_genes x n_cells PRTB matrix and reads at most
# max.de.genes rows of it. Computing the DE genes first lets CalcPerturbSig run on
# just those features, and the list is handed back via RunMixscale(DE.gene = ...).
#
# This mirrors the selection block of RunMixscale exactly (the split.by/harmonize
# path is not replicated because the pipeline never uses it -- see the stop() below).
mixscale_de_genes <- function(object,
                              labels,
                              nt.class.name,
                              de.assay = "RNA",
                              logfc.threshold = 0,
                              fine.mode = FALSE,
                              fine.mode.labels = NULL,
                              pval.cutoff = 5e-2,
                              split.by = NULL,
                              harmonize = FALSE,
                              verbose = TRUE) {
  if (!is.null(split.by) || isTRUE(harmonize)) {
    stop("mixscale_de_genes() does not replicate the split.by/harmonize path; ",
         "call RunMixscale directly if you need it.")
  }
  Idents(object = object) <- labels
  all_cells <- colnames(object)
  genes <- setdiff(unique(object[[labels]][all_cells, 1]), nt.class.name)
  nt.cells <- all_cells[object[[labels]][all_cells, 1] == nt.class.name]

  out <- list()
  for (gene in genes) {
    if (verbose) message("[de.genes] ", gene)
    orig.guide.cells <- all_cells[object[[labels]][all_cells, 1] == gene]

    if (isTRUE(fine.mode)) {
      md <- object[[]][orig.guide.cells, , drop = FALSE]
      guides <- setdiff(unique(md[[fine.mode.labels]]), nt.class.name)
      all.de.genes <- c()
      for (gd in guides) {
        gd.cells <- rownames(md)[which(md[[fine.mode.labels]] == gd)]
        de <- Seurat:::TopDEGenesMixscape(
          object = object, ident.1 = gd.cells, ident.2 = nt.cells,
          de.assay = de.assay, logfc.threshold = logfc.threshold,
          labels = fine.mode.labels, verbose = FALSE, pval.cutoff = pval.cutoff
        )
        all.de.genes <- c(all.de.genes, de)
      }
      all.de.genes <- unique(all.de.genes)
    } else {
      all.de.genes <- Seurat:::TopDEGenesMixscape(
        object = object, ident.1 = orig.guide.cells, ident.2 = nt.cells,
        de.assay = de.assay, logfc.threshold = logfc.threshold,
        labels = labels, verbose = FALSE, pval.cutoff = pval.cutoff
      )
    }
    # Deliberately uncapped: RunMixscale applies max.de.genes itself, and the cap
    # is all.de.genes[1:max.de.genes], so letting it do so keeps the two paths
    # byte-identical rather than relying on the cap being idempotent.
    out[[gene]] <- all.de.genes
  }
  out
}
