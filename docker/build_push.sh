#!/usr/bin/env bash
# Build and push the mixscale image with the HVG kernels compiled in.
#
#   ./build_push.sh              # builds and pushes xyscmbbb/r-mixscale:1.1.9
#   ./build_push.sh 1.2.0        # different tag
#   PUSH=0 ./build_push.sh       # build and verify only, do not push
#
# The build context is assembled here from bin/, so bin/hvf_kernels.c and
# bin/Makevars stay the single source of truth. Requires a docker login that can
# push to the xyscmbbb namespace.
set -euo pipefail
TAG=${1:-1.1.9}
IMAGE=xyscmbbb/r-mixscale:$TAG
BASE=docker.io/xyscmbbb/r-mixscale:1.1.8
HERE=$(cd "$(dirname "$0")" && pwd)
BIN=${BIN:-$HERE/../bin}
PUSH=${PUSH:-1}

for f in hvf_kernels.c Makevars; do
  [ -f "$BIN/$f" ] || { echo "missing $BIN/$f"; exit 1; }
done

CTX=$(mktemp -d)
trap 'rm -rf "$CTX"' EXIT
cp "$BIN/hvf_kernels.c" "$BIN/Makevars" "$HERE/Dockerfile" "$CTX/"

echo "=== building $IMAGE from $BASE ==="
docker pull "$BASE"
docker build -t "$IMAGE" "$CTX"

# The Dockerfile already fails the build if the .so does not load. This re-checks
# it in a fresh container, and confirms the kernels agree with Seurat's own on a
# small matrix -- the same identical() gate as validate_hvf_kernels.R, cheap
# enough to run on every build.
echo "=== verifying $IMAGE ==="
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
and drop hvf_kernels.so from bin/ (and from promote.sh's file list) once every
runner is on it. The scripts fall back to a .so next to them, then to compiling
one, so an older image keeps working.
NOTE
