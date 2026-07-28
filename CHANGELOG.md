# task_spatial_simulators dev

Bug fixes:
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
  - `downstream` and `ks_statistic_sc_features`: normalise the real and the
    simulated dataset the same way, through a new `compute_logcounts()` helper.
    They used three different transforms between them, so the metrics partly
    measured the difference in normalisation. On the positive control,
    `crosscor_cosine` and `crosscor_mantel` now come out at exactly 1.
  - `downstream`: reuse the `spatial_cluster` computed by
    `process_datasets/generate_sim_spatialcluster` instead of running BayesSpace
    a second time. The two runs disagreed, which put the positive control at
    `clustering_ari` 0.55 while the best simulator scored 0.22.

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
  - `synsim`
  - `zinbwave`

Control methods under `src/control_methods/`:
  - `negative_normal`
  - `negative_shuffle`
  - `positive`

Metrics under `src/metrics/`:
  - `ks_statistic_gene_cell`
  - `ks_statistic_sc_features/`

Documentation:
  - Check the `README.md` and `INSTRUCTIONS.md` for how to use and extend the benchmark.
