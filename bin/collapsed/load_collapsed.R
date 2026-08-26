# Compile and load the collapsed Gamma-Poisson fitter.
#
# sourceCpp caches by source hash, so the ~20 s build happens once per image and
# every later job in the same container reuses the shared object. Point
# MIXSCALE_CPP_CACHE at a persistent path (e.g. a bind mount) to share the cache
# across containers.
load_collapsed_glm_gp <- function(dir = NULL) {
  if (is.null(dir)) {
    self <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
    dir <- file.path(dirname(self), "collapsed")
  }
  suppressPackageStartupMessages(library(Rcpp))
  cache <- Sys.getenv("MIXSCALE_CPP_CACHE", unset = file.path(tempdir(), "mixscale_cpp"))
  dir.create(cache, showWarnings = FALSE, recursive = TRUE)
  t0 <- Sys.time()
  Rcpp::sourceCpp(file.path(dir, "gp_weighted.cpp"), cacheDir = cache, env = globalenv())
  Rcpp::sourceCpp(file.path(dir, "gp_overdisp_w.cpp"), cacheDir = cache, env = globalenv())
  source(file.path(dir, "collapsed_glm_gp.R"))
  message(sprintf("[03] collapsed glm_gp loaded (%.0f s)",
                  as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  invisible(TRUE)
}
