# Adds the transpose-free HVG kernels to the mixscale image.
#
# Seurat's VST calls SparseRowVar2/SparseRowVarStd, which take
# Eigen::SparseMatrix<double> BY VALUE and then transpose in place -- three full
# copies of the counts matrix per call, twice, none of it visible to R's gc().
# At 305k cells that alone set step 2's container peak at 64.6 GB. bin/hvf_kernels.c
# replaces both with a single column-major pass, bit-identical by construction.
#
# Building here rather than in bin/ means no binary is tracked in the repo and no
# job compiles anything at run time. The pipeline finds the result through
# MIXSCALE_HVF_SO; if this variable is unset it falls back to a .so next to the
# scripts, and failing that compiles one into a tempdir.
FROM docker.io/xyscmbbb/r-mixscale:1.1.8

# Makevars sets CC=gcc: this image's R Makeconf points at a conda toolchain that
# is not installed here. The kernels are plain C against R's own API, so the
# system gcc is a safe substitute -- glibc is shared and there is no C++ runtime
# to mismatch, which is exactly why they are not written in C++.
COPY hvf_kernels.c Makevars /opt/mixscale_hvf/

ENV MIXSCALE_HVF_SO=/opt/mixscale_hvf/hvf_kernels.so

RUN cd /opt/mixscale_hvf \
 && R_MAKEVARS_USER=/opt/mixscale_hvf/Makevars R CMD SHLIB hvf_kernels.c \
 && rm -f /opt/mixscale_hvf/hvf_kernels.o \
 && test -f "$MIXSCALE_HVF_SO" \
 && Rscript -e 'dyn.load(Sys.getenv("MIXSCALE_HVF_SO")); stopifnot(is.loaded("rowVar2_cm"), is.loaded("rowVarStd_cm")); cat("hvf kernels loaded OK\n")'
