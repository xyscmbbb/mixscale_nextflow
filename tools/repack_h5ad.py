#!/usr/bin/env python
"""
Repack one cell line's h5ad into the layout step 1 can subset without reading it all.

Why this exists
---------------
The pipeline used to be fed a per-target shard: `pert_<GENE>.h5ad` = that target's
cells + ALL NT cells. At genome-wide scale that is ~5,000 shards which are ~99% NT
cells, so the NT block gets written and stored 5,000 times -- tens of TB, and a
shard-generation run that takes longer than the analysis.

Step 1 can now read just the cells one job wants straight out of the whole
filtered_correct_pairs.h5ad, on two conditions:

  1. X must be CSR (cells x genes). A CSR indptr is over CELLS, so a set of cells
     is a set of row slices. The h5ads this pipeline is fed are csc_matrix, whose
     indptr is over GENES -- there a cell subset touches every nonzero in the file,
     which is no better than reading the whole thing.
     CSR(cells x genes) is also bit-identical to CSC(genes x cells), the layout
     Seurat wants, so it removes step 1's transpose as well.

  2. Rows should be SORTED by target_gene, NT first. Then the cells one job wants
     are exactly two contiguous runs -- the NT block and the target's block -- so
     the read is two sequential slices instead of tens of thousands of scattered
     ones. Correctness does not depend on this; speed does.

Output: <dst>.h5ad plus <dst>.h5ad.rowindex.csv (target_gene -> [start, end) row
range), which is what tells you a run really did collapse to two slices.

If you control the code that WRITES filtered_correct_pairs.h5ad, do the sort and
`adata.X = adata.X.tocsr()` there instead and skip this tool entirely.

Usage:
  repack_h5ad.py SRC.h5ad DST.h5ad --nt-label "(INTERGENIC CONTROL)" \
      [--target-col target_gene] [--keep-csv keep_cell_filter.csv] \
      [--keep-col keep_cell] [--cell-col cell] [--max-nnz-in-memory 2e9]
"""
import argparse, os, sys, time

import numpy as np
import h5py


def log(msg):
    print(f"[repack] {msg}", flush=True)


def read_obs_column(f, col):
    """obs column as a str array, categorical or not."""
    g = f["obs"][col]
    if isinstance(g, h5py.Group):  # categorical
        cats = np.asarray([c.decode() if isinstance(c, bytes) else str(c)
                           for c in g["categories"][:]])
        return cats[g["codes"][:]]
    v = g[:]
    if v.dtype.kind == "S":
        return np.asarray([x.decode() for x in v])
    return v.astype(str)


def read_obs_index(f):
    obs = f["obs"]
    name = obs.attrs.get("_index", "_index")
    if isinstance(name, bytes):
        name = name.decode()
    v = obs[name][:]
    return np.asarray([x.decode() if isinstance(x, bytes) else str(x) for x in v])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("src")
    ap.add_argument("dst")
    ap.add_argument("--nt-label", required=True)
    ap.add_argument("--target-col", default="target_gene")
    ap.add_argument("--keep-csv", default=None,
                    help="optional allowlist CSV (the PotC keep_cell_filter.csv)")
    ap.add_argument("--keep-col", default="keep_cell")
    ap.add_argument("--cell-col", default="cell")
    ap.add_argument("--compression", default="gzip",
                    help="gzip (default) or none. gzip roughly thirds the file and "
                         "h5py still decompresses only the chunks a row slice touches.")
    ap.add_argument("--max-nnz-in-memory", type=float, default=2e9,
                    help="above this, transpose in row-block passes instead of at once")
    args = ap.parse_args()

    t0 = time.time()
    import anndata as ad
    import scipy.sparse as sp

    with h5py.File(args.src, "r") as f:
        enc = f["X"].attrs.get("encoding-type", b"")
        enc = enc.decode() if isinstance(enc, bytes) else str(enc)
        shape = tuple(int(x) for x in f["X"].attrs["shape"])
        nnz = int(f["X"]["indices"].shape[0])
        tg = read_obs_column(f, args.target_col)
        cells = read_obs_index(f)
    log(f"src {os.path.basename(args.src)}: {enc}, {shape[0]} cells x {shape[1]} genes, "
        f"nnz = {nnz:,} ({nnz / max(shape[0], 1):.0f}/cell)")

    if args.nt_label not in set(tg):
        sys.exit(f"[repack] nt-label {args.nt_label!r} not found in obs[{args.target_col!r}]")

    keep = np.ones(shape[0], dtype=bool)
    if args.keep_csv:
        import pandas as pd
        kdf = pd.read_csv(args.keep_csv)
        allowed = set(kdf.loc[kdf[args.keep_col].astype(str).str.lower()
                              .isin(["true", "t", "1", "yes", "y"]), args.cell_col].astype(str))
        keep = np.fromiter((c in allowed for c in cells), bool, len(cells))
        log(f"allowlist {os.path.basename(args.keep_csv)}: keeping {keep.sum():,} of {len(keep):,} cells")

    # NT first (sorts before everything as ""), then targets alphabetically.
    # A stable sort keeps each block in the file's original cell order, so a
    # repack never silently reorders cells within a target.
    key = np.where(tg == args.nt_label, "", tg)[keep]
    order = np.flatnonzero(keep)[np.argsort(key, kind="stable")]
    log(f"row order: NT block ({int((tg[order] == args.nt_label).sum()):,} cells) "
        f"then {len(set(tg[order])) - 1} targets")

    if nnz > args.max_nnz_in_memory:
        sys.exit(f"[repack] nnz = {nnz:,} exceeds --max-nnz-in-memory "
                 f"({args.max_nnz_in_memory:,.0f}). The blocked path is not implemented; "
                 f"either raise the limit on a big-RAM box, or -- much better -- write "
                 f"the source h5ad sorted and CSR upstream instead of repacking it here.")

    log("reading ...")
    a = ad.read_h5ad(args.src)
    log(f"read in {time.time() - t0:.1f}s")

    t = time.time(); a = a[order].copy(); log(f"reorder {time.time() - t:.1f}s")
    t = time.time(); a.X = sp.csr_matrix(a.X); log(f"to CSR {time.time() - t:.1f}s")

    comp = None if args.compression.lower() in ("none", "", "0") else args.compression
    t = time.time(); a.write_h5ad(args.dst, compression=comp)
    log(f"write {time.time() - t:.1f}s -> {args.dst} "
        f"({os.path.getsize(args.dst) / 1e9:.2f} GB)")

    # Sidecar row index: proof that each target is one contiguous run, and a
    # cheap way for a caller to see how many cells a job will actually load.
    import pandas as pd
    tg_out = a.obs[args.target_col].astype(str).values
    chg = np.flatnonzero(np.r_[True, tg_out[1:] != tg_out[:-1]])
    ends = np.r_[chg[1:], len(tg_out)]
    idx = pd.DataFrame({"target_gene": tg_out[chg], "start": chg, "end": ends,
                        "n_cells": ends - chg})
    if len(idx) != idx["target_gene"].nunique():
        sys.exit("[repack] a target_gene is split across more than one run -- "
                 "the sort did not group it; refusing to write a misleading row index")
    idx.to_csv(args.dst + ".rowindex.csv", index=False)
    log(f"row index -> {args.dst}.rowindex.csv ({len(idx)} contiguous blocks)")
    log(f"total {time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()
