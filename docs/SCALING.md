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
  `gene_chunk x design_rows` (0.5 GB by default). The 71.8 GB peak is the NT
  counts held several times over: the Seurat object's copy, step 3's
  `count_data_sparse` column subset, the fold-change scaled copies, and the
  gene-major transpose, at ~17 GB each.
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

The fix follows from the collapse itself. The counts are touched in exactly one
place, `chunk_stats_cpp`, which reduces a gene chunk to per-group sums. Group
membership is per-cell and cell blocks are disjoint, so `S_r`, `L_r`, group
sizes and the count histogram are all **additive across blocks**. Holding the
counts as a list of column blocks and accumulating per block is exact, keeps
every block under the index cap, and -- when blocks are read from the h5ad one
at a time -- makes peak memory O(one block) instead of O(all NT cells).

## 5. Cost

Terra/Cromwell custom machine type, us-central1, on-demand: $0.033174/vCPU/hr
and $0.004446/GB/hr, plus ~90 s of billed boot and localization per task. RAM is
provisioned at 1.3x measured peak (an OOM kill costs a full retry).

Step 3 alone at 300k NT: 8 vCPU + 94 GB for 1,175 billed s = **$0.2239/job,
$1,120 for 5,000** -- already 2.2x over the $500 budget before steps 1 and 2.

Note the split: the 94 GB contributes $0.418/hr against the 8 vCPUs' $0.265/hr.
**Peak memory, not runtime, is what makes this job expensive**, which is why the
block-streaming change in section 4 matters more than any further speedup.

## Status

Measured and settled: the collapse identity, parity, CM-scale numbers, the 300k
step-3 number and its breakdown, the dgCMatrix ceiling.

Still open: step 1 and step 2 at 300k NT (queued), the 450k saturation point
(queued), and the block-streaming and LOO-grouping changes from sections 3-4,
which are prototyped outside the repo and not yet merged here.
