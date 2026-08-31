# Scaling the per-perturbation job to a genome-wide screen

Target: ~5,000 perturbations over ~3M cells, one Terra job per perturbation,
each starting from a sharded `.h5ad` and running step 1 -> step 2 -> step 3.

Every number below is measured on this machine (16 vCPU, 125 GB) unless it is
explicitly labelled a projection. Peak memory is the whole process tree's
anonymous RSS sampled at 0.5 s, which is the figure comparable to Cromwell's
`peak_rss` -- NOT cgroup v2 `memory.peak`, which includes page cache and
overstates a job that reads a multi-GB h5ad.

## 1. Why the NT design rows collapse

In a single-cell-line job the design is `~ 1 + weight + log_ct` (p = 3).
Non-targeting cells have `weight == 0` exactly, so an NT cell's design row is
`(1, 0, log_ct)` -- determined solely by `log_ct = log1p(nCount_RNA)`. Cells
sharing a design row share mu, so every per-cell sum in the Gamma-Poisson
likelihood becomes a per-group sum. This is an algebraic identity, not an
approximation: same design, same optimizer, same objective.

The payoff grows with the NT pool because the number of *distinct integer
library sizes* saturates while the cell count does not. Measured on CM's
empirical library-size distribution:

| NT cells | distinct log_ct | NT collapse |
|---|---|---|
| 10,000  |  7,700 |  1.3x |
| 100,000 | 24,400 |  4.1x |
| 300,000 | 28,773 | 10.4x |
| 600,000 | 29,500 | 20.3x |

This is the key scaling property: **the fit cost is O(genes x design rows), and
design rows saturate.** Doubling the NT pool past ~300k adds almost no solver
work.

## 2. Measured: CM shard (16,009 cells = 10,899 NT + 5,110 PDHA1)

At this size the collapse is only 1.29x, so this is close to a worst case for
the method -- and it still wins on both axes because the collapse also removes
the dense `Mu`.

| target | step | stock | collapsed |
|---|---|---|---|
| PDHA1 | step3 | 546.0 s / 24.22 GB | 303.9 s / 8.38 GB |
| RPL9  | step3 | 592.4 s / 25.66 GB | 312.6 s / 8.51 GB |
| MOSMO | step3 | 605.1 s / 26.89 GB | -- |

Steps 1 and 2 at this size: 61.6-69.8 s / 6.95-7.69 GB and 87.9-95.0 s /
8.97-9.81 GB.

### Numerical parity (the reason to trust the above)

Stock vs collapsed on identical input (PDHA1, 22,898 shared genes):

- `avg_log2FC` identical to 0
- median `|dlog10 p|` = 4.6e-14
- among the 3,472 genes significant in either arm: max `|dlog10 p|` = 5.6e-06,
  max `|dbeta|` = 5.1e-08
- **DEG sets identical at every threshold tested** (BH < 0.25/0.1/0.05/0.02/0.01)

The ~50 genes with large beta disagreement are non-identifiable (no expression
in the perturbed group, so no finite MLE); both arms return p ~ 1 for them.

Against the Terra reference (`CM_mixscale_all_genes_log2fc.csv`): log2FC
correlation 0.999384 (MOSMO) and 0.999183 (PDHA1). Collapsed-vs-reference and
stock-vs-reference agree to three decimals, which is what shows the collapse
adds no error of its own.

## 3. Measured: 300k NT (305,110 cells, 31,729 genes)

Step 3, collapsed: **1,084.9 s / 71.82 GB**. Collapse 9.4x (305,110 cells ->
32,424 design rows).

> **Superseded as the step-3 reference.** That figure predates both the LOO
> hoist and the current step-2 object. Re-measured on the same input after
> those landed, step 3 is **512.4 s / 59.59 GB**, and with the zero-copy counts
> fetch (Stage A) **394.8 s / 37.86 GB** at 32,257 design rows. Compare against
> 512.4 s / 59.59 GB, not against 1,084.9 s / 71.82 GB -- the older number
> credits Stage A with savings the LOO hoist had already banked.

Solver breakdown of the main genome-wide fit (297 s of the 1,085 s):

| stage | time |
|---|---|
| sufficient statistics | 32.1 s |
| beta pass 1           | 61.8 s |
| overdispersion MLE    | 132.0 s |
| shrinkage             | 10.1 s |
| beta pass 2 + se      | 60.9 s |

Two things this rules out and one it rules in:

- **Sharing NT statistics across jobs is not worth it.** The sufficient-stats
  pass is the only part that depends on the NT counts, and therefore the only
  part reusable across all 5,000 jobs. It is 32 s of 297 s (14%) -- not worth
  building and distributing a multi-GB shared table for.
- **The solver is not the memory problem.** Its working set is
  `gene_chunk x design_rows` (0.5 GB by default). The peak is the NT counts
  held several times over: the Seurat object's copy, step 3's
  `count_data_sparse` column subset, the fold-change scaled copies, and the
  gene-major transpose, at ~17 GB each. The transpose is real and is at
  `collapsed/collapsed_glm_gp.R:118` (`Yt <- Matrix::t(Y)`), which an earlier
  revision of this document wrongly said did not exist; `Y` is untouched after
  that line, so freeing the caller's chain there is the remaining step-3 lever
  (Stage B, not yet done -- step 3 is already under 40 GB without it).
- **The leave-one-out tail is now the largest runtime component.** Each of the
  ~101 LOO fits does 0.7 s of solver work but re-enters the solver from scratch,
  recomputing the design-row grouping over all 305,110 cells and rebuilding a
  305,110-row model matrix. Across all NT cells the grouping is *identical* in
  every LOO fit -- only the perturbed cells' weight column changes.

## 4. Hard limits

`dgCMatrix` stores row indices as 32-bit signed integers, so a sparse matrix
cannot exceed 2^31-1 = 2.147e9 nonzeros. At CM's density (4,689 nnz/cell) that
caps a single matrix at **~458,000 cells**, regardless of available RAM:

| NT cells | nonzeros | sparse size | fits in dgCMatrix? |
|---|---|---|---|
| 300,000 | 1.41 G | 16.9 GB | yes |
| 450,000 | 2.11 G | 25.3 GB | barely |
| 600,000 | 2.81 G | 33.8 GB | **no** |

Step 1 hits this sooner than step 3 does: it converts `adata.X` to COO and
materializes full-length `i`, `j`, and `x` vectors in R (5.7 + 5.7 + 11.4 GB at
300k) before building the matrix and transposing it.

**The fix that shipped was not block-streaming.** An h5ad's sparse triple
`(indptr, indices, data)` already *is* a dgCMatrix's `(p, i, x)`, so step 1 now
reads the three HDF5 datasets in slices into preallocated R vectors and hands
them to `new("dgCMatrix", ...)` -- no scipy COO, no 1-based copies, no dgT
intermediate, and `new()` shares the vectors rather than copying them. Measured
at 305k: **950.9 s / 112.63 GB -> 270.3 s / 32.61 GB**, DIGEST matching
REFERENCE. A CSR `(cells x genes)` shard is bit-identical to the CSC
`(genes x cells)` Seurat wants, so the transpose is free; writing shards CSR was
expected to halve this again but **measured only +0.09 GB and +3% wall for CSC**,
so it is not worth an upstream change.

Block-streaming remains the answer to the *index cap* rather than to memory: the
counts are touched in exactly one place, `chunk_stats_cpp`, which reduces a gene
chunk to per-group sums. Group membership is per-cell and cell blocks are
disjoint, so `S_r`, `L_r`, group sizes and the count histogram are all
**additive across blocks**. That is what a >458k-cell run would need; nothing at
305k does.

## 4b. Step 2 at 300k NT: seven fixes, and why the first five moved no memory

Measured at 305,110 cells, whole step: **4,734.0 s / 64.62 GB** before any work,
**1,255.6 s / 64.64 GB** after six fixes. Runtime fell 3.8x. **The peak did not
move at all**, and understanding why is the point of this section.

The six fixes, all gated bit-identical on HeLa/RPL9 and CM/PDHA1 (eight checks:
cells, genes, VariableFeatures, assay meta.data, RNA layers, layer counts,
meta.data, `tools$RunMixscale`):

1. **HVG write-back detached.** `obj[["RNA"]]@meta.data <- v` hands the assay out
   through a temporary, so R sees two references and duplicates on assignment --
   and duplicating an `Assay5` deep-copies `@layers`, i.e. the 17 GB counts
   matrix, to store a few columns of statistics. Detaching the assay first makes
   the local the sole reference and R writes in place. Note a function wrapper
   does NOT fix this: an argument bumps the reference count too.
2. **Data-layer write detached**, same reason.
3. **Forked presto** in DE gene selection (each gene row is an independent
   rank-sum; chunks reassemble in chunk order, so the result is bit-identical).
4. **Forked `RANN::nn2`.** Tree build is negligible and query is everything, so
   splitting the query across forks is near-linear.
5. **Raw-layer fetch in CalcPerturbSig.** `GetAssayData`/`LayerData` re-attach
   dimnames on fetch, and Seurat v5 stores layers *without* them -- so the
   accessor copies the whole 17 GB layer merely to label it. Read
   `assay@layers[["..."]]` positionally and name the small result instead.
6. **`scale_pca_hvf`.** ScaleData and RunPCA read only the 2,000 variable
   features, but ScaleData stores its result by writing the assay back -- another
   deep copy of `@layers`. Run both on a variable-features-only temporary object
   and transplant the reduction back; `scale.data` is deliberately not
   transplanted, since the caller drops it two lines later.

### Why none of that lowered the peak: R's heap growth trigger

R does not collect when memory is free; it collects when an allocation would
cross its *gc trigger*, and it grows that trigger to ~1.6x live whenever live
grows. The normalize stage genuinely holds counts and data at once -- 32.9 GB
live, irreducibly (see below) -- which pushed the trigger to **51.3 GiB**. Every
later stage was then allowed to fill to 51.3 GiB with garbage before R would
collect, and two of them did exactly that:

| stage | max used | gc trigger | pinned? |
|---|---|---|---|
| HVG (vst on counts) | 32.66 GB | 38.5 GiB | no -- genuine |
| normalize (blocked) | 43.25 GB | 51.3 GiB | no -- genuine |
| ScaleData + RunPCA  | 51.78 GB | 51.3 GiB | **yes** |
| DE gene selection   | 51.74 GB | 51.3 GiB | **yes** |

A stage whose max-used lands on the trigger is not a stage that needs that
memory. The clinching evidence is fix 6: it provably removes ~9 GB of allocation,
and its own stage peak moved 51.78 -> 51.94 (nothing) while the *following*
stage fell 26.02 -> 16.77 GB. A real saving the stage peak cannot see is what
trigger-pinning looks like.

Setting **`R_GC_MEM_GROW=0`** (a container env var; read at R startup, zero
parity risk) makes the trigger track live instead of 1.6x live:

| stage | default | `R_GC_MEM_GROW=0` |
|---|---|---|
| normalize trigger | 51.3 GiB | 42.8 GiB |
| ScaleData + RunPCA | 51.94 GB | **43.42 GB** |
| ScaleData + RunPCA time | 87.2 s | 89.0 s |
| HVG time | 67.6 s | 96.5 s |

Full stage table at 305,110 cells, `--threads 8`, both arms carrying all six
fixes above:

| stage | default | `R_GC_MEM_GROW=0` |
|---|---|---|
| HVG (vst on counts)         |  67.6 s / 32.66 GB |  68.3 s / 32.57 GB |
| normalize (blocked)         |  51.8 s / 43.25 GB |  53.9 s / 43.21 GB |
| ScaleData + RunPCA          |  87.2 s / 51.94 GB |  89.0 s / **43.42 GB** |
| DietSeurat + drop scale.data|   1.9 s / 16.77 GB |   2.1 s / 16.77 GB |
| DE gene selection (wilcox)  | 342.0 s / 51.74 GB | **279.5 s / 43.26 GB** |
| CalcPerturbSig (kNN + PRTB) | 528.5 s / 17.17 GB | **480.2 s / 17.17 GB** |
| RunMixscale (scores + LOO)  |  67.4 s / 18.77 GB |  66.6 s / 18.77 GB |
| **container peak (cgroup anon)** | **64.64 GB** | **64.52 GB** |

It is close to free on time, and in places it is a win: DE gene selection got
18% *faster* while shedding 8.5 GB, and CalcPerturbSig 9% faster -- less heap
pressure means less copy-on-write traffic through the forked workers.

(A first measurement had HVG at 96.5 s here and it was reported as a +43% cost of
the policy. It was not: a HeLa gate was running on the same box at that moment. A
second run under the identical env var gave 68.3 s against the default arm's
67.6 s. **Never benchmark a memory arm against a shared box** -- the first arm
happened to be the one contended, and the artefact was large enough to look like
a real effect.)

### But it did not lower the container peak, and that is the real lesson

**R's heap fell 51.94 -> 43.42 GB and the cgroup peak did not move (64.64 ->
64.52 GB).** So ~21 GB of the container's peak is not in R's `gc()` accounting at
all.

I first attributed that gap to the forked `mclapply` workers. **That was wrong,
and it is worth recording how it was disproven, because the reasoning was
plausible.** Two measurements kill it:

| arm | forks? | max R-heap stage | container peak |
|---|---|---|---|
| pre-fix `bin` | **no `--threads` at all** | 46.43 GB | 64.62 GB |
| dev, 5 fixes | 8 | 51.78 | 64.62 |
| + `scale_pca_hvf` | 8 | 51.94 | 64.64 |
| + `R_GC_MEM_GROW=0` | 8 | 43.42 | 64.52 |
| + `de_target_nnz 4e7` | 8, 3x smaller chunks | 43.42 | 64.59 |

1. The pre-fix arm ran through `fullrun_ab.sh`, which passes **no `--threads`** to
   step 2. It was fully serial, forked nothing, and peaked at the same 64.62 GB.
2. Cutting the per-worker chunk 3x (`--de_target_nnz` 1.2e8 -> 4e7) moved the
   container peak by 0.07 GB and cost **257 s** (1166.4 -> 1423.7 s), with the DE
   stage going 279.5 -> 555.6 s for an identical 43.26 GB heap.

Five configurations differing in threading, chunk size, gc policy and ~8 GB of
live set span **0.12 GB** of container peak. That is not a fork term; it is a
fixed cost that none of the stage instrumentation observes. **Keep
`de_target_nnz` at 1.2e8** -- shrinking it is a pure loss.

Two consequences that do survive:

- Measure the **cgroup `anon`** figure, never R's `gc()` max-used, when deciding
  a machine size. The two answer different questions and here they disagreed by
  21 GB. This also means `R_GC_MEM_GROW=0` buys ~nothing on the VM size it
  actually provisions -- keep it for the runtime win in the DE stages (18% faster
  DE selection, 9% faster CalcPerturbSig), not as a memory fix.
- A stage timer that reports R's `gc()` peak cannot see the binding constraint.
  To locate the container peak in *time*, sample cgroup `anon` continuously and
  stamp the container start (`dev305_step2prof.sh` does this, with
  `stamp_lines.sh` recording when each log line arrives so the trace can be cut
  into named phases and `analyze_prof.py` reducing the two into a table).

That profiling run is what settled it, and the answer was the first suspect, not
the second: the peak lands **inside** a timed stage, and something allocates
outside R's allocator. See the next section.

### The 48 GB nobody could see: Seurat's VST transposes the counts matrix, twice

`analyze_prof.py` on the 1194.2 s / 64.59 GB arm reported `PEAK 64.59 GB at
t+139.2s`, inside HVG:

| phase | dur | container max |
|---|---|---|
| prologue (R startup + package load) | 47.9 s | 8.11 |
| readRDS | 27.9 | 16.63 |
| **HVG (vst)** | **71.8** | **64.59**  <== PEAK |
| normalize | 54.0 | 43.31 |
| ScaleData + RunPCA | 88.1 | 43.44 |
| DietSeurat | 2.1 | 17.43 |
| DE gene selection | 316.9 | 43.80 |
| CalcPerturbSig | 479.3 | 20.03 |
| re-read counts | 63.2 | 18.89 |
| saveRDS | 38.5 | 18.89 |

`readRDS` is clean: a linear climb 8.11 -> 16.63 GB, exactly one copy. HVG then
makes two identical excursions -- `16.63 -> 32.62 -> 48.60 -> 64.58`, hold 13 s,
back down -- while R's `gc()` reported the whole stage at 32.57 GB.

The cause is in Seurat's own C++ (`src/data_manipulation.cpp`). `VST.dgCMatrix`
calls `SparseRowVar2` and `SparseRowVarStd` once each, and **both take
`Eigen::SparseMatrix<double> mat` by value and open with
`mat = mat.transpose();`** -- they want per-gene statistics out of a cell-major
matrix. That is three full-size copies live at once (Rcpp's copy into Eigen, the
transpose temporary, the assignment result), per call, twice. 3 x 16 GB, and
invisible to `gc()` because none of it is in R's heap.

**Fix (`hvf_kernels.c` + `hvf_lowmem.R`):** a transpose is not needed at all. One
column-major pass over the dgCMatrix slots, read in place, accumulating into a
per-gene vector. Bit-identity is the point, so the accumulation order matches
Seurat's exactly: Seurat walks column *k* of the transposed matrix, i.e. gene
*k*'s nonzeros in ascending original-column order, and visiting original columns
`j = 0..n-1` in order hands each gene its terms in that same ascending-*j*
sequence, so every partial sum is the identical double at every step. (A blocked
or partial-sum accumulation would **not** be bit-identical.) The tail terms and
the divide are applied in Seurat's order too.

Written in plain C, not C++: this image's R `Makeconf` points at a conda
toolchain that is not installed, and building with the system `g++` would link a
different libstdc++ than R itself uses. C against R's own API has no such ABI
surface. `Makevars` sets `CC = gcc`; build with
`R CMD SHLIB hvf_kernels.c` **inside the container**.

Gated on `identical()`, never `all.equal()` -- the latter would hide exactly the
last-bit drift that could reorder the 2,000-gene selection. On CM (38,292 genes x
16,009 cells): `SparseRowVar2` 1.2 s -> 0.1 s, `SparseRowVarStd` 1.2 s -> 0.2 s,
full VST data.frame `identical: TRUE`, top-2000 selection `identical: TRUE`. On
HeLa/RPL9 the full step-2 output object matches on all eight checks
(`STEP2 OBJECTS IDENTICAL`), including the `tools$RunMixscale` slot step 3 reads.

At 305k the HVG excursion disappears: **16.63 -> 32.63 GB and back**, and the
stage runs **71.8 s -> 21.6 s**. What remains above `readRDS` is one counts copy
that R's `gc()` *does* see (heap 32.57 GB), from the assay duplication in
`hvf_detached.R`; it is not worth chasing because normalize and DE selection sit
above it anyway.

Re-profiled end to end, same object, same 8 threads: **1194.2 s / 64.59 GB ->
1158.2 s / 43.94 GB**. The peak moved to DE gene selection, whose window maxes at
43.94 GB, and every other stage now sits below it:

| phase | dur | container max |
|---|---|---|
| prologue | 12.3 s | 0.64 |
| readRDS | 84.1 | 16.66 |
| HVG (vst) | 55.0 | 32.63 |
| normalize | 88.7 | 43.30 |
| ScaleData + RunPCA | 2.1 | 43.37 |
| DietSeurat | 320.3 | 17.43 |
| **DE gene selection** | 486.8 | **43.94**  <== PEAK |
| CalcPerturbSig | 4.2 | 20.03 |
| re-read counts | 63.3 | 18.90 |
| saveRDS | 40.5 | 18.90 |

(The windows are cut at the arrival of the line that *ends* each phase, so each
row's duration belongs to the phase named on the row above it.)

### The layer that nobody needs whole

Independent of where the 21 GB lives, there is one avoidable 17 GB allocation.
Reading what each stage actually consumes:

| stage | reads |
|---|---|
| HVG | **counts** (`vst on counts`) |
| ScaleData + RunPCA | `data`, but only the **2,000 variable features** |
| DE gene selection | `data`, **full width** -- and already blocked over nnz |
| CalcPerturbSig | `data`, but only the **DE genes** (137 on CM, **4** at 305k) |
| RunMixscale | the PRTB assay |
| saveRDS | counts, re-read from disk |

Exactly one consumer needs the full-width `data` layer, and it is the one that is
already chunked. Yet `NormalizeData` materialises all 17.2 GB of it and holds it
live from that point to the end of the run.

The fix is to not materialise it: compute the lognorm per block inside the DE
selection, which already iterates blocks over counts, and materialise a `data`
layer for the HVGs only. It is **exact, not approximate** -- `log1p(x / colSum *
1e4)` per block from counts is the identical arithmetic in the identical order
that `lognorm_blocked` already performs; we simply stop keeping the result. At
305k that is counts 17.2 GB + an HVG-width `data` of ~3 GB + block temporaries,
so a heap floor near 25 GB against today's 43.4 -- and it dissolves the normalize
floor described next, because normalize's 28.5 GB live set stops existing.

### The normalize floor is irreducible, so do not try

It is tempting to rewrite `lognorm_blocked` to overwrite `@x` in place instead of
allocating a second matrix. **Measured: it buys nothing.** Detaching `@x` and
reattaching it are both free, but the first block-subassignment duplicates the
whole vector once regardless (the reference count is saturated and there is no
supported way to clear it), after which the remaining blocks do write in place.
So in-place peaks at `@i + @p + @x + one @x copy` = the same 28.5 GB as today's
version, which already allocates exactly one new `@x` and shares `@i`/`@p`.

That 28.5 GB live is therefore step 2's floor, and via the trigger it sets a
**~43.8 GB R-heap floor for every stage after it**. Going below requires not
holding the counts matrix whole -- streaming it from the RDS in column blocks --
which is a step-1/IO redesign, not a step-2 fix.

### The cgroup peak is not the R heap

R's heap peaked at 51.94 GB but the container peaked at 64.64 GB, and the 13 GB
gap was Eigen memory inside Seurat's VST kernels -- outside R's allocator, so
`gc()` could not see it. Size VMs from the cgroup figure and treat the R heap as
a diagnostic only. Any R package with an `Rcpp`/`RcppEigen` entry point can hide
an arbitrary multiple of its input this way.

## 5. Cost

Terra/Cromwell custom machine type, us-central1, on-demand: $0.033174/vCPU/hr
and $0.004446/GB/hr, plus ~90 s of billed boot and localization per task. RAM is
provisioned at 1.3x measured peak (an OOM kill costs a full retry).

**The extended-memory line is the single most important thing on this page.**
A GCE custom VM may carry 0.9-6.5 GB of RAM per vCPU at the base rate; every
gigabyte above `6.5 x vCPU` bills as *extended memory* at $0.009550/GB/hr, about
2.15x. Two rules follow, and they point in opposite directions to intuition:

- **Cores below the line are nearly free.** At 8 vCPU, 4 / 6 / 8 / 10 / 12 vCPU
  all price the same to four decimals when memory is the binding term.
- **Never buy cores to escape the line.** Widening to 14 vCPU does clear the
  extended tier, and it costs *more* ($0.0780 vs $0.0751 spot), because the
  added vCPU-hours exceed the tier saving. Cut the memory instead.

So the provisioned figure that matters is `1.3 x peak` measured against
`6.5 x vCPU`. At 8 vCPU that line sits at 52 GB provisioned = **40 GB measured
peak**. Landing a step under 40 GB is worth far more than making it faster.

Spot is ~4.7x cheaper on every line item and a sub-hour job carries only ~2-3%
preemption exposure per attempt. On-demand cannot reach the $0.10/perturbation
target by any configuration measured here; spot can.

Cromwell sizes a **separate VM per task**, so steps 1/2/3 are costed against
their own peaks and a single flat `memory` setting wastes the difference on
every job. Keep `cpus_step1..3` / `mem_step1..3` split.

That was written when step 3 alone at 300k NT needed 8 vCPU + 94 GB for 1,175
billed s = **$0.2239/job, $1,120 for 5,000** -- 2.2x over the $500 budget before
steps 1 and 2 were even counted. The 94 GB contributed $0.418/hr against the
8 vCPUs' $0.265/hr, which is the whole reason peak memory, not runtime, drives
this bill.

Where the three steps stand now, all at 8 vCPU, spot, `cost3.py`:

| step-2 peak | provisioned | spot / perturbation | x 5,000 |
|---|---|---|---|
| 64.6 GB (before the VST fix) | 84 G -- extended tier | $0.0867 | $433 |
| 50 GB | 65 G -- extended tier | $0.0724 | $362 |
| **43.94 GB (measured, after the VST fix)** | **58 G -- extended tier** | **$0.0636** | **$318** |
| 40 GB | 52 G -- on the line | $0.0626 | $313 |

The runtime work alone had already cleared the $0.10 target ($0.0986 ->
$0.0867). The VST fix is what brings the *provisioned* figure within sight of
the extended-memory line; the last few GB to actually cross it would have to
come from the `data`-layer change above.

## 6. No per-target shards

Every number above was measured on a `pert_<GENE>.h5ad` shard: that target's cells
plus the whole NT pool, cut upstream and handed to the job. That input does not
survive the screen it was built for.

At 3M cells with 305k NT, a shard is **~99.8% NT cells**, so ~5,000 targets means
writing and storing ~5,000 copies of the same NT pool:

| | per-target shards | whole-file + in-job selection |
|---|---|---|
| artefacts per cell line | ~5,000 | **1** |
| storage | ~16 TB | ~31 GB |
| storage cost | ~$320/month | **~$0.62/month** |
| build time | ~14 months serial, ~27 days at 16x | **~40 min** |
| per job | download a ~3.2 GB shard | ~3.2 GB of ranged reads, nothing on disk |

The build-time figure is the anndata gzip write rate measured here, 12.8 MB/s.
The sharding stage costs more than the analysis it feeds.

### What step 1 does instead

`--h5ad` is now the whole cell-line file. Step 1 computes its own cell set --
`target_gene in {perturb_gene, nt_label}`, optionally AND-ed with a boolean
allowlist column named by `--keep_cell_col` -- and reads only those rows'
nonzeros: one bounded HDF5 slice per contiguous run of wanted cells, straight
into `new("dgCMatrix", i=, p=, x=)`.

**Gated `identical()`.** Reading PDHA1's cells out of the whole 125,439-cell CM
file yields a Seurat object identical to the one built from the pre-made shard --
counts matrix *and* `meta.data`. Measured through Nextflow: **2 runs, 12.7% of the
file's 588M nonzeros, 6.8 s**.

One real bug fell out of that gate and is worth remembering: **anndata drops unused
categorical levels when it subsets**, so a shard arrived carrying only its 2
`target_gene` levels while an in-job subset of the whole file kept all 26. The
values compare equal as character and the counts matrix is untouched, so only
`identical()` on `meta.data` catches it. Step 1 now calls `droplevels()` after
subsetting. At 5,000 targets the un-fixed version would carry 4,999 empty factor
levels into every downstream `model.matrix`.

### End to end, shard input vs whole-file input

The step-1 gate above holds cell ORDER fixed. The repacked file does not: it is
sorted NT-first, so a job reading it gets the same cells in a different order than
the old shard did. That is not free -- floating-point sums are order-dependent --
so the whole pipeline was run twice on the same code, once from `CM_PDHA1.h5ad`
and once from the sorted repack, and the outputs compared:

| quantity | result |
|---|---|
| gene set fitted | identical, 22,898 both |
| `log2FC` | **bit-identical** (`identical()` TRUE) |
| `mixscale_score` | **bit-identical** |
| `p_weight` | median abs diff 1.1e-14, max 1.0e-2 |
| DEG Jaccard, BH<0.05 and BH<0.01, both with and without a 0.25 log2FC gate | **1.000000** |

The p-value drift is the Gamma-Poisson overdispersion MLE reaching a slightly
different stopping point when its per-cell terms are summed in a different order.
It is confined to genes that are not significant: the single largest drift is on a
gene at p = 0.55, and across every gene with p < 0.05 in either run the largest
drift is **1.4e-8**. Nothing crosses a threshold, in either direction, at any gate.

So the input change is output-neutral in every respect that is read downstream. If
you need `p_weight` itself bit-reproducible against an old shard run, repack
without the sort (correctness never depended on it -- only the two-runs-instead-of
-thousands read speed does).

### Two properties of the file, both load-bearing

Neither is optional, and the source files as written satisfy neither:

- **X must be CSR.** A CSR `indptr` is over *cells*, so a cell set is a set of row
  slices. These files are written `csc_matrix`, whose `indptr` is over *genes* --
  there a cell subset touches every nonzero (measured: 15.5 s to scan CM's 588M
  nnz locally, and over GCS it means fetching all 31 GB to use 3.2 GB).
  CSR(cells x genes) is also bit-identical to CSC(genes x cells), which is the
  layout Seurat wants, so this removes step 1's transpose as well (section 4).
- **Rows sorted by `target_gene`, NT first.** Each target is then one contiguous
  range and a job reads exactly two runs. Unsorted, NT+PDHA1 on the real CM file
  is **13,965 scattered runs**.

`tools/repack_h5ad.py` does both, once per cell line, and writes a
`.rowindex.csv` sidecar giving each target's `[start, end)`. It refuses to write
that sidecar if any target ended up split across more than one run, so the file
existing is the assertion that the sort worked.

**It does the transpose in memory and stops above `--max-nnz-in-memory` (2e9).**
At 3M cells (~14e9 nnz) it will not run, and a blocked transpose is the wrong fix:
do the sort and `adata.X = adata.X.tocsr()` in whatever *writes*
`filtered_correct_pairs.h5ad`, where it is nearly free. The PotC producer side is
`lib/repack_h5ad.py` on branch `single-guide-large-scale` of
`PotC_pilot_perturb_analysis`, which replaced that repo's sharding stage.

### Reading it without downloading it

Image **1.1.11** adds `gcsfs`. Step 1 opens a `gs://` URI through it and hands the
Python file object to `h5py`, which issues ranged GETs for just those row runs, so
the 31 GB file never lands on disk. `obs`/`var` come off the same handle via
`anndata.io.read_elem`. Same-region GCS reads are not billed for egress.

Three designs were costed; ranged reads won:

| | transfer | disk_gb | $/job |
|---|---|---|---|
| whole-file `gsutil cp` | ~31 GB | 100 | ~$0.087 |
| NT block + tiny per-target shards | ~3.2 GB | 100 | ~$0.075 |
| **ranged reads via gcsfs** | **~3.2 GB** | **60** | **~$0.071** |

### Cost

Same 8 vCPU / 58 GB spot VM as section 5; what changes is the disk, because the
input no longer has to fit on it. Disk is billed at the full rate with **no spot
discount**, which is why 40 GB of it is worth more than it looks:

| | billed | compute | disk | total | x 5,000 |
|---|---|---|---|---|---|
| localised shard, 100 GB | 32.0 min | $0.0629 | $0.0124 | $0.0754 | $377 |
| **ranged reads, 60 GB** | 32.5 min | $0.0638 | **$0.0076** | **$0.0714** | **$357** |

The billed time goes slightly *up*: ~32 s of `gsutil cp` is replaced by ~60 s of
ranged reads inside step 1, which is a wash against a $0.0048 disk saving. The
shard-build and storage cost, ~16 TB and weeks of wall clock, is what actually
goes away -- and it never appeared on this page because it was somebody else's
line item.

**Still under the $0.10/perturbation target, with the whole sharding stage gone.**

## Status

All three steps are measured at 305,110 cells, 8 vCPU, and every change below is
gated bit-identical against the stage it replaced. End to end: **113 min and a
112.6 GB worst step -> 30 min and a 43.9 GB worst step**, $0.3764 -> **$0.0636**
per perturbation on spot, **$318 for 5,000 targets**:

| step | stock | now | gate |
|---|---|---|---|
| 1 (h5ad -> Seurat) | 950.9 s / 112.63 GB | **270.3 s / 32.61 GB** | DIGEST == REFERENCE |
| 1 input | per-target shard (~16 TB for 5,000) | **one repacked file per cell line** | counts + meta.data `identical()` to the shard build |
| 1 input, end to end | -- | -- | `log2FC` and `mixscale_score` bit-identical, DEG Jaccard 1.000000 |
| 2 (preprocess + Mixscale) | 4734.0 s / 64.62 GB | **1158.2 s / 43.94 GB** | `STEP2 OBJECTS IDENTICAL` on HeLa/RPL9 |
| 3 (weighted DE) | does not fit (dense `Mu` = 77 GB) | **380.1 s / 37.85 GB** | DE CSV byte-identical |

Settled: the collapse identity and its parity, the CM/HeLa/305k/455k numbers,
the dgCMatrix nonzero ceiling, the extended-memory line, and -- for step 2 --
that the container peak was Seurat's VST transposes rather than forking, chunk
size, or R's gc trigger.

Open:

- **`hvf_kernels.so` is hand-built into `bin/`.** It belongs in the image
  (`R CMD SHLIB hvf_kernels.c` in the Dockerfile), otherwise all 5,000 jobs
  either recompile it or depend on a binary tracked in the repo.
- The `data`-layer change in "The layer that nobody needs whole" is designed and
  documented but **not implemented**. It is what would take step 2 across the
  40 GB line.
- Step 3 Stage B (freeing the counts allocation after `Matrix::t(Y)` in
  `collapsed/collapsed_glm_gp.R`) is parked -- 37.85 GB is already under the line.
- The 47.9 s R-startup + package-load prologue is ~4% of step 2's wall clock and
  is pure fixed cost repeated 5,000 times.
- Whether the 1e-8 end-to-end DE difference seen once on CM came from step 2 or
  from 8-thread BLAS reduction order (`determinism_ab.sh` needs re-running).
- The 458k-cell dgCMatrix ceiling is unaddressed; block-streaming (section 4) is
  the answer if a target ever exceeds it.
- **`tools/repack_h5ad.py` will not run at 3M cells** -- its transpose is
  in-memory and it stops above 2e9 nonzeros (section 6). The intended fix is
  upstream (write the source sorted and CSR), not a blocked transpose here.
- **The gcsfs path has not been exercised against real GCS.** Everything except
  the transport is verified in-image -- `gcsfs`/`fsspec` import, `h5py` accepting
  a Python file object and slicing it correctly, `anndata.io.read_elem` on `obs`
  and `var` -- and the whole selection path is gated on a local file. What is
  untested is a live `gs://` open: ranged-GET latency at 3M-cell scale, and
  whether `block_size = 16 MiB` is the right read-ahead.
