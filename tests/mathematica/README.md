# Mathematica Baseline Generator

This directory contains the reusable Mathematica baseline generators for the
construction test data.

## GCA-GCA baseline

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
output paths above. Load the Mathematica module before any `math -script`
invocation on Perlmutter.

## GCA-constant-latitude baseline

Default output:

`tests/input/gca_constlat_cases_with_baseline.csv`

Perlmutter command:

```bash
module load mathematica
math -script tests/mathematica/generate_gca_constlat_baseline.m \
  tests/input/gca_constlat_cases_with_baseline.csv
```

If you omit the CSV argument, the script falls back to the default output path
above. Load the Mathematica module before any `math -script` invocation on
Perlmutter.
