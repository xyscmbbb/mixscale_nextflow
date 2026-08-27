# Compile and load the collapsed Gamma-Poisson fitter.
#
# Two things make the ~15 s build disappear after the first time:
#
#   MIXSCALE_COLLAPSED_SRC  a directory holding a copy of the .cpp sources that
#                           the image precompiled at build time (see docker/).
#   MIXSCALE_CPP_CACHE      where sourceCpp keeps its build cache.
#
# Image 1.1.10 sets both, so a Terra job loads the solver in well under a second
# instead of compiling it. Neither is required: with both unset this compiles
# into a tempdir exactly as before, which is what happens on 1.1.9 and older.

# Make sure the compiler R's Makeconf names is actually reachable.
#
# This image's R is a conda R: `R CMD config CXX` returns
# x86_64-conda-linux-gnu-c++, which exists only in the conda environment's own
# bin directory -- and that directory is NOT on the image's default PATH. Left
# alone, sourceCpp dies with
#   sh: 1: x86_64-conda-linux-gnu-c++: not found
#   WARNING: The tools required to build C++ code for R were not found.
# after steps 1 and 2 have already run. It never surfaced in benchmarking
# because every bench invocation bind-mounts a prebuilt cache, so the compiler
# was never actually called.
#
# The environment prefix is derived from R.home() rather than hardcoded, so this
# tracks whichever R is running. The directory is APPENDED, never prepended: it
# also contains R, Rscript, make, python and python3, and prepending would
# shadow the /usr/local/bin/Rscript shim the pipeline is invoked through.
ensure_r_toolchain <- function() {
  rbin <- file.path(R.home("bin"), "R")
  cxx <- tryCatch(system2(rbin, c("CMD", "config", "CXX"), stdout = TRUE, stderr = FALSE),
                  error = function(e) character(0))
  if (!length(cxx)) return(invisible(FALSE))
  prog <- strsplit(trimws(cxx[1]), "[[:space:]]+")[[1]][1]
  if (nzchar(Sys.which(prog))) return(invisible(FALSE))

  envdir <- dirname(dirname(R.home()))      # <prefix>/lib/R -> <prefix>
  cand <- file.path(envdir, "bin")
  if (!file.exists(file.path(cand, prog))) {
    stop(sprintf("R's compiler '%s' is not on PATH and is not in %s", prog, cand))
  }
  Sys.setenv(PATH = paste(Sys.getenv("PATH"), cand, sep = .Platform$path.sep))
  message(sprintf("[03] added %s to PATH for the R toolchain (%s)", cand, prog))
  invisible(TRUE)
}

# Decide which directory to compile the .cpp sources from.
#
# sourceCpp's cache key includes the source PATH, not just its contents (checked:
# the same file under two directories builds twice and leaves two sourcecpp_*
# entries in one cache). So a cache baked at image-build time is only reused if
# the job compiles from the very same path -- hence MIXSCALE_COLLAPSED_SRC.
#
# The redirect is gated on an md5 match against the checked-out sources. If this
# branch's .cpp differ from the copy baked into the image, the baked copy is
# ignored and the repo's own sources are compiled. Stale kernels silently
# replacing edited ones is the one failure mode worth spending a hash on.
collapsed_src_dir <- function(dir, files) {
  baked <- Sys.getenv("MIXSCALE_COLLAPSED_SRC", unset = "")
  if (!nzchar(baked) || !dir.exists(baked)) return(dir)
  if (!all(file.exists(file.path(baked, files)))) return(dir)
  a <- unname(tools::md5sum(file.path(dir, files)))
  b <- unname(tools::md5sum(file.path(baked, files)))
  if (anyNA(a) || anyNA(b) || !identical(a, b)) {
    message("[03] ", baked, " does not match the checked-out sources; compiling from ", dir)
    return(dir)
  }
  baked
}

load_collapsed_glm_gp <- function(dir = NULL) {
  if (is.null(dir)) {
    self <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
    dir <- file.path(dirname(self), "collapsed")
  }
  suppressPackageStartupMessages(library(Rcpp))
  ensure_r_toolchain()
  cache <- Sys.getenv("MIXSCALE_CPP_CACHE", unset = file.path(tempdir(), "mixscale_cpp"))
  dir.create(cache, showWarnings = FALSE, recursive = TRUE)

  cpp <- c("gp_weighted.cpp", "gp_overdisp_w.cpp")
  src <- collapsed_src_dir(dir, cpp)

  t0 <- Sys.time()
  for (f in cpp) Rcpp::sourceCpp(file.path(src, f), cacheDir = cache, env = globalenv())
  # The R-level solver always comes from the checked-out tree, never the image.
  source(file.path(dir, "collapsed_glm_gp.R"))
  message(sprintf("[03] collapsed glm_gp loaded (%.1f s%s)",
                  as.numeric(difftime(Sys.time(), t0, units = "secs")),
                  if (identical(src, dir)) "" else ", prebuilt"))
  invisible(TRUE)
}
