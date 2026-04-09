# Mathematica Baseline Generator

This directory contains the reusable baseline generator for the GCA-GCA
intersection construction test data.

Default input:

`/pscratch/sd/h/hyvchen/regrid_bench/GCAGCA/arcs/gca_gca_pairs_seed20251104_N100.csv`

Default output:

`tests/input/gca_gca_pairs_seed20251104_N100_with_baseline.csv`

Perlmutter command:

```bash
module load mathematica
math -script tests/mathematica/generate_gca_gca_baseline.m \
  tests/input/gca_gca_pairs_seed20251104_N100.csv \
  tests/input/gca_gca_pairs_seed20251104_N100_with_baseline.csv
```
If you omit the CSV arguments, the script falls back to the default input and
output paths above.
