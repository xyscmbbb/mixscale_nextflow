# VST.dgCMatrix, verbatim, with Seurat's two Eigen-transposing kernels swapped
# for the transpose-free ones in hvf_kernels.c. Every other line -- the loess
# fit, the 0L initialisation of variance.expected, the ordering, the rank
# assignment -- is copied from Seurat so the returned data.frame is identical
# column by column. See hvf_kernels.c for why the accumulation order matches.
# Three places, in order: the .so the image built (MIXSCALE_HVF_SO), a .so sitting
# next to these sources, and finally a compile from hvf_kernels.c. The compile
# targets a per-process tempdir rather than the source directory, because 5,000
# concurrent jobs share one bin/ and would race on the same output file. It costs
# about a second and exists so a missing or ABI-mismatched .so degrades to slow
# rather than to a failed job.
load_hvf_kernels <- function(dir) {
  if (is.loaded("rowVar2_cm")) return(invisible(TRUE))
  # MIXSCALE_HVF_SO is set by the image, which builds the kernels at image-build
  # time. That is the intended production path: bin/ then carries no binary, and
  # nothing is compiled per job.
  for (so in c(Sys.getenv("MIXSCALE_HVF_SO", ""), file.path(dir, "hvf_kernels.so"))) {
    if (nzchar(so) && file.exists(so) &&
        !inherits(try(dyn.load(so), silent = TRUE), "try-error"))
      return(invisible(TRUE))
  }
  src <- file.path(dir, "hvf_kernels.c")
  if (!file.exists(src))
    stop("neither hvf_kernels.so nor hvf_kernels.c found in ", dir)
  bd <- file.path(tempdir(), "hvf_kernels_build")
  dir.create(bd, showWarnings = FALSE, recursive = TRUE)
  file.copy(src, bd, overwrite = TRUE)
  mv <- file.path(dir, "Makevars")
  if (file.exists(mv)) file.copy(mv, bd, overwrite = TRUE)
  message("[hvf] hvf_kernels.so unusable in ", dir, "; compiling into ", bd)
  # system2() has no working-directory argument, so change into the build dir.
  owd <- setwd(bd); on.exit(setwd(owd), add = TRUE)
  ok <- system2("R", c("CMD", "SHLIB", "hvf_kernels.c"),
                stdout = FALSE, stderr = FALSE,
                env = if (file.exists(mv))
                        paste0("R_MAKEVARS_USER=", file.path(bd, "Makevars"))
                      else character())
  built <- file.path(bd, "hvf_kernels.so")
  if (!identical(ok, 0L) || !file.exists(built))
    stop("could not build hvf_kernels.so -- run: R CMD SHLIB hvf_kernels.c in ", dir)
  dyn.load(built)
  invisible(TRUE)
}

rowVar2_cm <- function(i, p, x, n_genes, n_cells, mu)
  .Call("rowVar2_cm", i, p, x, as.integer(n_genes), as.integer(n_cells), mu)

rowVarStd_cm <- function(i, p, x, n_genes, n_cells, mu, sd, vmax)
  .Call("rowVarStd_cm", i, p, x, as.integer(n_genes), as.integer(n_cells),
        mu, sd, as.double(vmax))

vst_lowmem <- function(data, nselect = 2000L, span = 0.3, clip = NULL,
                       verbose = TRUE) {
  nfeatures <- nrow(x = data)
  hvf.info <- SeuratObject:::EmptyDF(n = nfeatures)
  hvf.info$mean <- Matrix::rowMeans(x = data)
  hvf.info$variance <- rowVar2_cm(data@i, data@p, data@x,
                                  nrow(data), ncol(data), hvf.info$mean)
  hvf.info$variance.expected <- 0L
  not.const <- hvf.info$variance > 0
  fit <- loess(formula = log10(x = variance) ~ log10(x = mean),
               data = hvf.info[not.const, , drop = TRUE], span = span)
  hvf.info$variance.expected[not.const] <- 10^fit$fitted
  hvf.info$variance.standardized <- rowVarStd_cm(
    data@i, data@p, data@x, nrow(data), ncol(data),
    hvf.info$mean, sqrt(x = hvf.info$variance.expected),
    if (is.null(clip)) sqrt(x = ncol(x = data)) else clip)
  hvf.info$variable <- FALSE
  hvf.info$rank <- NA
  vf <- head(x = order(hvf.info$variance.standardized, decreasing = TRUE),
             n = nselect)
  hvf.info$variable[vf] <- TRUE
  hvf.info$rank[vf] <- seq_along(along.with = vf)
  return(hvf.info)
}
