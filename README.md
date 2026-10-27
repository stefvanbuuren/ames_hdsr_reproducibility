# Ames Housing simulations for HDSR

Reproducibility code for **Missing Data for What? Evaluating Imputation Methods by Their Intended Use**.

## Files

- `ames_helpers.R`: preprocessing, reconstruction methods, point and distributional prediction (Goals 1–3), and shared helpers.
- `run_goals123.R`: production run for Table 1 (100 splits; 100 imputations; FCS training m=100 and filling m=1000).
- `ames_goal4.R`: inferential simulation functions (Goal 4), including HC3 standard errors and PMM sensitivity analysis.
- `run_goal4.R`: production run for Table 2 (1,000 samples of 1,000 houses, 20 imputations).
- `run_goals123_parallel.R`: optional parallel equivalent of `run_goals123.R` (Unix/macOS, 12 cores via `parallel::mclapply`); gives identical results.
- `renv.lock`: package lockfile for the production environment.
- `production_outputs/`: raw results, summaries, logs and `sessionInfo` from the verification run, with `SHA256SUMS`.

## Run

Run from the directory containing the four R scripts:

```r
source('run_goals123.R')
source('run_goal4.R')
```

These runs are computationally intensive. For smoke tests, temporarily reduce `nrep` in each runner, and `m_train`/`m_fill` in `run_goals123.R`.

## Dependencies

- R packages: `modeldata`, `mice`, `sandwich`.
- **Pinned `mice` dependency:** All simulations (Goals 1–4) were run using the development branch of `amices/mice` at commit [`df7b02a064ccc3ded3e53ae92c70d91bc938ef77`](https://github.com/amices/mice/commit/df7b02a064ccc3ded3e53ae92c70d91bc938ef77). This commit provides the FCS train/fill interface (`tasks = 'train'` / `'fill'` and `models = ...`). Ordinary CRAN releases may not support this interface.

Install the exact revision before running either simulation:

```r
install.packages('remotes') # if needed
remotes::install_github(
  'amices/mice',
  ref = 'df7b02a064ccc3ded3e53ae92c70d91bc938ef77'
)
```

Verify the installed revision with `packageDescription('mice')[c('Version', 'RemoteSha')]`. Both runners write `sessionInfo` text files.

Alternatively, restore the full package set (66 packages, including `mice` from GitHub at the pinned commit) from the lockfile:

```r
install.packages('renv') # if needed
renv::restore(lockfile = 'renv.lock')
```

### Reported production environment

- **R:** 4.5.2 (2025-10-31), `aarch64-apple-darwin20`
- **OS:** macOS Tahoe 26.5.2 (Apple Silicon)
- **mice:** 3.18.1, development commit pinned above
- **modeldata:** 1.6.0
- **sandwich:** 3.1.3 (confirmed with `packageVersion("sandwich")`)
- **Other key package versions:** Matrix 1.7-6, dplyr 1.2.1, tidyr 1.3.2, MASS 7.3-66, broom 1.0.13, lme4 2.1-0
- **BLAS:** Apple Accelerate; **LAPACK:** 3.12.1
- **Locale:** `en_US.UTF-8`; **time zone:** `Europe/Amsterdam`

The installed `sandwich` version (3.1.3) was confirmed separately using `packageVersion("sandwich")`; it was not loaded when the supplied `sessionInfo()` was captured. See `session_info_production.txt` for the environment details.

## Design notes

- Living area (`Gr_Liv_Area`) is the incomplete predictor X; `log(Sale_Price)` is Y. Auxiliary Z comprises `log(Lot_Area)`, `Lot_Config` (with rare levels merged), and `Bldg_Type`.
- Five houses with living area above 4,000 square feet are excluded (N=2,925).
- Goals 1–3: 100 random 75%/25% training/test splits; one-third MCAR masking independently in each. True test outcomes are withheld from imputation and prediction, and used only for scoring. The predictive specification is Y ~ X + Z; the oracle uses true X.
- The methods with suffix `ZY` may use observed training outcomes for imputation; test outcomes are never available. In the existing single-imputation and MI pipelines, missing test X is generated using the training-fitted X|Z model. The FCS procedure instead fills both X and Y from stored conditional models.
- FCS temporarily masks 10% of training Y to learn/store its conditional model; test Y is entirely masked. Predictive draws are summarized using empirical 5th/95th percentiles.
- Other MI predictive intervals are obtained from the mixture of conditional normal predictive distributions, **not** by averaging individual interval endpoints. Plug-in uncertainty estimates do not fully integrate parameter uncertainty.
- Goal 4: 1,000 independent samples with replacement from the fixed empirical population; one-third MCAR X masking; the inferential estimand is the finite-population linear-projection slope in Y ~ X (without Z). HC3 and Rubin pooling are used for interval construction.
- Goal 1 reconstruction RMSE for MI is the average RMSE across imputed datasets, not the RMSE of the average imputation.
- FCS is a prediction-only procedure; its Goal 1 and slope diagnostics are NA.

## Before publication

1. ~~The `mice` commit and reported R environment are documented; record a package lockfile.~~ See `renv.lock`.
2. ~~Run both production scripts in the release environment and compare summaries against the reference results.~~ Done on 2026-10-26; see *Verification* below. The reference results match every value in Tables 1 and 2 of `HDSR_draft_8.docx` and the Goal 4 numbers quoted in the text (population slope, PMM sensitivity analysis), after rounding.
3. ~~Publish the raw outputs and logs, or provide checksums for archived results.~~ See `production_outputs/` and `production_outputs/SHA256SUMS`.
4. Add a permanent repository URL and release tag/DOI once published.

## Published numerical reference results

The `reference_results/` directory contains the **original production summaries supplied by the authors**:

- `table1_goals123.csv`: 100 train/test replications for Goals 1–3, including observed-X diagnostics.
- `table2_goal4_with_pmm.csv`: 1,000 inferential replications for Goal 4, including the PMM sensitivity analysis.

After running the scripts locally, compare `ames_goals123_FCS_production_summary.csv`
with `reference_results/table1_goals123.csv`, and `ames_goal4_summary.csv`
with `reference_results/table2_goal4_with_pmm.csv`.
Floating-point differences may occur across machines or library versions.

## Verification

On 2026-10-26 the cleaned scripts were rerun at full production settings in an environment that matches the reported production environment: R 4.5.2 on `aarch64-apple-darwin20` (macOS, Apple Silicon), `mice` 3.18.1 at commit `df7b02a`, `modeldata` 1.6.0 and `sandwich` 3.1-3.

- **Goal 4:** `run_goal4.R` was run unmodified (1,000 replications, about 7 minutes).
- **Goals 1–3:** `run_goals123_parallel.R` was used. It calls `ames_one()` for replications 1–100 with the settings in `run_goals123.R` (m=100, m_train=100, m_fill=1000), spread over 12 cores with `parallel::mclapply`, and summarises with `summary_table()` (8 minutes in total, compared with about 1 hour sequentially). Each replication sets its own seed (`seed + rep`), so its result does not depend on run order or parallel execution.

Both regenerated summaries are **byte-identical** to the reference files:

| File | SHA-256 |
|---|---|
| `table1_goals123.csv` / `ames_goals123_FCS_production_summary.csv` | `7385fe8102362301a164875971275cbe87de817ac548ceff22cdc601a71ec8c7` |
| `table2_goal4_with_pmm.csv` / `ames_goal4_summary.csv` | `a5174e55750c3cb34c2aee8427054ab7fc5d9054874d8faa87244fbd44f814b3` |

The raw `.rds` results, summaries, logs and `sessionInfo` files from this run are in `production_outputs/`; verify them with `shasum -a 256 -c SHA256SUMS` from inside that directory. Notes on these files:

- `goals123.log` contains only the elapsed time, because progress messages from the forked workers are not captured.
- `run_goals123_parallel.R` does not write `sessionInfo`, so `sessionInfo_goals123.txt` was captured afterwards (2026-10-27) from the same R installation.
- Both verification sessions ran under the `C.UTF-8` locale rather than the production `en_US.UTF-8`; this did not affect the results.
- `.rds` checksums depend on the R version and serialisation settings, so compare the summary CSVs rather than the `.rds` hashes on other machines.

This confirms that the code reproduces the results in the production environment. It does not cover other platforms, BLAS libraries or package versions. Every value in manuscript Tables 1 and 2 (`HDSR_draft_8.docx`) matches the reference results after rounding.
