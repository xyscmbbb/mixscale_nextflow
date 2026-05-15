# mixscale_nextflow_v6

Three-step Nextflow workflow:

1. `01_convert_h5ad_to_obj.R`: h5ad -> Seurat object RDS
2. `02_mixscale_preprocess.R`: Mixscale preprocessing -> Mixscale object RDS
3. `03_Run_wmvRegDE_scaled_debug.R`: weighted Mixscale DE -> final CSV outputs

This v5 pipeline keeps your original `03_Run_wmvRegDE_scaled_debug.R` unchanged. Step 3 runs through `bin/03_wrapper_patch_pcols.R`, which creates a temporary patched copy in the Nextflow work directory to avoid the `object 'pcols' not found` error. The original script file is not modified.

Step 3 uses `stdbuf -oL -eL ... | tee step3_wmvregde.log` and `debug true`, so messages from `glm_gp`/R should stream to the terminal and are also saved in the process work directory as `step3_wmvregde.log`.

## Run example

```bash
nextflow run main.nf -resume \
  --h5ad /home/unix/xuyushan/mixscale/pilot_10/shards/pert_ZNHIT6.h5ad \
  --pair_csv /home/unix/xuyushan/mixscale/pilot_10/VSMC_per_cell_target.csv \
  --perturb_gene ZNHIT6 \
  --target_gene_col target_gene \
  --cell_col cell \
  --guide_col pair_key \
  --nt_label ONE_INTERGENIC_SITE \
  --outdir /home/unix/xuyushan/mixscale/pilot_10/results/ZNHIT6 \
  -ansi-log false \
  -with-trace
```

## Final outputs

Only these two final files are published to `--outdir`:

```text
pert_<GENE>_de_res_df.csv
pert_<GENE>_meta_data.csv
```

Intermediate RDS files stay in Nextflow work directories.


## v6 notes

- Uses `docker.io/xyscmbbb/r-mixscale:1.1.7` by default.
- No wrapper script is used. Step 3 calls `bin/03_Run_wmvRegDE_scaled_debug.R` directly.
- The original notebook is not modified; the R scripts contain copied workflow logic.
- Step 3 uses `stdbuf -oL -eL ... | tee step3_wmvregde.log` so glmGamPoi/glm_gp messages stream during execution.
