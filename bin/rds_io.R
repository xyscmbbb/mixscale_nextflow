# saveRDS() compresses with gzip at level 6 and gives no way to change it. For a
# step-1/step-2 object the payload is a counts matrix -- integer-ish data that
# barely compresses -- so level 6 spends a great deal of CPU for very little.
# Measured on a CM step-1 object (902 MB uncompressed), per GB of object:
#
#   level 6 (saveRDS default)   49.6 s/GB   5.01x smaller
#   level 3                     18.3 s/GB   4.29x
#   level 1                     10.5 s/GB   4.48x
#   no compression               2.0 s/GB   1.00x
#
# At 305,110 cells the object is ~17 GB, so the default is ~840 s -- most of
# step 1's 892 s end-to-end. Level 1 is 4.7x faster for a file 12% larger, which
# is the right trade when the file is an intermediate that gets written once and
# read once. Writing through a gzfile() connection is what exposes the level;
# the result is still an ordinary gzipped RDS, so readRDS() reads it unchanged
# and nothing downstream needs to know.
#
# MIXSCALE_RDS_COMPRESS overrides the level: 0 disables compression entirely
# (2.0 s/GB, but a ~17 GB intermediate to move around), 1-9 pick a gzip level.
save_rds_fast <- function(object, path, level = NULL) {
  if (is.null(level)) {
    lv <- suppressWarnings(as.integer(Sys.getenv("MIXSCALE_RDS_COMPRESS", "1")))
    level <- if (is.na(lv)) 1L else lv
  }
  if (level <= 0L) return(invisible(saveRDS(object, path, compress = FALSE)))
  con <- gzfile(path, "wb", compression = min(9L, as.integer(level)))
  on.exit(close(con), add = TRUE)
  invisible(saveRDS(object, con))
}
