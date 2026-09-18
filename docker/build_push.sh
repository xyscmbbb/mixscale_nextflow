#!/usr/bin/env bash
# Build and push a mixscale image layer.
#
#   ./build_push.sh              # builds and pushes the default tag below
#   ./build_push.sh 1.1.9        # rebuild the HVG-kernel layer instead
#   PUSH=0 ./build_push.sh       # build and verify only, do not push
#
# Each tag is one layer on the one before it:
#   1.1.9   HVG kernels compiled in, exported as MIXSCALE_HVF_SO
#   1.1.10  collapsed Gamma-Poisson solver precompiled, conda toolchain on PATH
#   1.1.11  gcsfs, so --h5ad can be a gs:// URI read in place (no localisation)
#
# The build context is assembled here from bin/, so the sources in bin/ stay the
# single source of truth and no binary is tracked in the repo. Requires a docker
# login that can push to the xyscmbbb namespace.
set -euo pipefail
TAG=${1:-1.1.11}
IMAGE=xyscmbbb/r-mixscale:$TAG
HERE=$(cd "$(dirname "$0")" && pwd)
BIN=${BIN:-$HERE/../bin}
PUSH=${PUSH:-1}

case "$TAG" in
  1.1.9)
    BASE=docker.io/xyscmbbb/r-mixscale:1.1.8
    DOCKERFILE=$HERE/Dockerfile
    SRC=("$BIN/hvf_kernels.c" "$BIN/Makevars")
    ;;
  1.1.10)
    BASE=docker.io/xyscmbbb/r-mixscale:1.1.9
    DOCKERFILE=$HERE/Dockerfile.1.1.10
    SRC=("$BIN/collapsed/gp_weighted.cpp" "$BIN/collapsed/gp_overdisp_w.cpp")
    ;;
  *)
    # Nothing is compiled from bin/ here, so the build context carries no sources.
    BASE=docker.io/xyscmbbb/r-mixscale:1.1.10
    DOCKERFILE=$HERE/Dockerfile.1.1.11
    SRC=()
    ;;
esac

for f in ${SRC[@]+"${SRC[@]}"}; do
  [ -f "$f" ] || { echo "missing $f"; exit 1; }
done

CTX=$(mktemp -d)
trap 'rm -rf "$CTX"' EXIT
[ ${#SRC[@]} -gt 0 ] && cp "${SRC[@]}" "$CTX/"
cp "$DOCKERFILE" "$CTX/Dockerfile"

echo "=== building $IMAGE from $BASE ==="
docker pull "$BASE"
docker build -t "$IMAGE" "$CTX"

echo "=== verifying $IMAGE ==="
if [ "$TAG" = "1.1.9" ]; then
  # The Dockerfile already fails the build if the .so does not load. This
  # re-checks it in a fresh container and confirms the kernels agree with
  # Seurat's own on a small matrix -- the same identical() gate as
  # validate_hvf_kernels.R, cheap enough to run on every build.
  docker run --rm "$IMAGE" Rscript -e '
    suppressMessages(library(Matrix)); suppressMessages(library(Seurat))
    dyn.load(Sys.getenv("MIXSCALE_HVF_SO"))
    set.seed(1)
    m <- as(matrix(rpois(20000, 0.4), 500, 40), "dgCMatrix")
    mu <- Matrix::rowMeans(m)
    a <- Seurat:::SparseRowVar2(mat = m, mu = mu, display_progress = FALSE)
    b <- .Call("rowVar2_cm", m@i, m@p, m@x, nrow(m), ncol(m), mu)
    sdv <- sqrt(a); vmax <- sqrt(ncol(m))
    c1 <- Seurat:::SparseRowVarStd(mat = m, mu = mu, sd = sdv, vmax = vmax,
                                   display_progress = FALSE)
    c2 <- .Call("rowVarStd_cm", m@i, m@p, m@x, nrow(m), ncol(m), mu, sdv,
                as.double(vmax))
    stopifnot(identical(a, b), identical(c1, c2))
    cat("kernels identical to Seurat: TRUE\n")'
else
  # Two things must hold, and only the second is new. First: the previous
  # layer's HVG kernels still load. Second, and the point of this layer: a
  # fresh container with NO cache mounted and NO PATH help loads the collapsed
  # solver from the baked cache. The elapsed time is the assertion -- anything
  # near 15 s means the cache was missed and every job would pay it again.
  docker run --rm -v "$BIN:/bin_ro:ro" "$IMAGE" Rscript -e '
    dyn.load(Sys.getenv("MIXSCALE_HVF_SO"))
    stopifnot(is.loaded("rowVar2_cm"), is.loaded("rowVarStd_cm"))
    cat("hvf kernels still OK\n")
    source("/bin_ro/collapsed/load_collapsed.R")
    t0 <- Sys.time()
    load_collapsed_glm_gp("/bin_ro/collapsed")
    el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    stopifnot(exists("fit_gp_weighted_cpp"), exists("gp_se_weighted_cpp"))
    cat(sprintf("collapsed solver loaded in %.1f s\n", el))
    if (el > 5) stop("solver was recompiled -- the baked cache was not reused")
    stopifnot(identical(Sys.getenv("RETICULATE_PYTHON"),
                        "/opt/miniconda3/envs/mixscale/bin/python"))
    cat("prebuilt cache reused, env OK\n")'
fi

if [ "$TAG" = "1.1.11" ]; then
  # The layer's whole purpose: h5py must accept a Python file object, because
  # that is how a gs:// h5ad is read without localising it. Checked against a
  # real (local) file object here -- reaching GCS needs credentials this build
  # does not have, but everything except the transport is exercised.
  docker run --rm "$IMAGE" /opt/miniconda3/envs/mixscale/bin/python -c '
import gcsfs, fsspec, h5py, numpy as np, tempfile, os
from anndata.io import read_elem
print("gcsfs", gcsfs.__version__, "fsspec", fsspec.__version__)
p = os.path.join(tempfile.mkdtemp(), "t.h5")
with h5py.File(p, "w") as f:
    f.create_dataset("d", data=np.arange(1000, dtype="int32"))
with h5py.File(open(p, "rb"), "r") as f:
    assert f["d"][10:20].tolist() == list(range(10, 20))
print("h5py over a file object: ranged slice OK")
gcsfs.GCSFileSystem
print("read_elem", read_elem.__name__, "OK")'
fi

if [ "$PUSH" = "1" ]; then
  echo "=== pushing $IMAGE ==="
  docker push "$IMAGE"
  echo "pushed $IMAGE"
else
  echo "PUSH=0, not pushing. Image is local: $IMAGE"
fi

cat <<NOTE

Next: point the pipeline at the new tag --
  nextflow.config   container = '$IMAGE'
  wdl/RunMixscale.wdl   docker_image default
NOTE
