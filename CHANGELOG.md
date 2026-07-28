# task_spatial_simulators dev

Bug fixes:
  - `ks_statistic_gene_cell` and `ks_statistic_sc_features`: the metric
    descriptions said Kolmogorov-Smirnov, but both call `ks::kde.test()`, which
    is a kernel density based two-sample test.
  - `file_dataset_sp.yaml`: `logcounts` is a double, not an integer.
  - `downstream`: `clustering_ari` was declared -Inf..+Inf and
    `ctdeconvolute_rmse` 0..+Inf, though both are bounded. Metric labels were
    the ids repeated back.
  - `splatter` and `symsim`: drop the `try()` around the per-cluster loop, which
    let a failed cluster pass silently and only surfaced later as a length
    mismatch.
  - `symsim`: sample gene lengths with replacement, which errored outright on
    datasets holding more genes than `gene_len_pool`.
  - `negative_normal`: round the generated values, so that `counts` holds
    integers as `file_simulated_dataset.yaml` declares.
  - `run_benchmark`: raise `uns_length_cutoff` from 15 to 50, so that
    `extract_uns_metadata` no longer drops the `metric_ids` of components that
    emit more than 15 metrics. All 28 `ks_statistic_gene_cell` metrics were
    reported as 100% missing in `run_2026-07-11_18-00-02` because of this.
  - `ks_statistic_gene_cell`: read `$Tstat` rather than `$tstat`, which
    `ks::kde.test()` does not define. Half of the metric values were being
    dropped on the way out.
  - `ks_statistic_gene_cell`: report NA when `ks::kde.test()` cannot be
    computed, instead of an empirical KS statistic or a penalty of 1. Both
    ranked as very good scores on the unbounded `zstat`/`Tstat` scale, so
    failing was rewarded.
  - `ks_statistic_gene_cell`: retry `ks::kde.test()` with deterministic jitter,
    as `ks_statistic_sc_features` already did, so that ties no longer turn into
    NAs. Inputs without variance report NA rather than surviving the jitter as
    a top score.

# task_spatial_simulators 0.1.0

First release of the spatial simulator benchmark.

Core task documentation and API:
  - Component types: Process Dataset, Method, Metric
  - File formats: Single-Cell Dataset, Spatial Dataset, Solution, Score

Dataset processing components for fetching datasets from the SpatialSimBench figshare:
  - Source: https://figshare.com/articles/dataset/SpatialSimBench_dataset/26054188
  - Transforms the h5ads into standardised components

Simulation methods under `src/methods/`:
  - `scdesign2`
  - `scdesign3_nb`
  - `scdesign3_poisson`
  - `sparsim`
  - `splatter`
  - `srtsim`
  - `symsim`
  - `zinbwave`

Control methods under `src/control_methods/`:
  - `negative_normal`
  - `negative_shuffle`
  - `positive`

Metrics under `src/metrics/`:
  - `downstream`
  - `ks_statistic_gene_cell`
  - `ks_statistic_sc_features`

Documentation:
  - Check the `README.md` and `INSTRUCTIONS.md` for how to use and extend the benchmark.
