# task_spatial_simulators dev

Minor changes:
  - `run_benchmark`: write the commit the workflow ran from into `task_info.yaml`,
    so a published result can be traced back to the code that produced it. Until
    now `task_info.yaml` was a verbatim copy of `_viash.yaml`, whose `version` is
    the revision (`build_main`) rather than a commit.

Bug fixes:
  - `downstream`: normalise the simulated dataset the same way as the real one,
    through a new `compute_logcounts()` helper. The simulated side used
    `log1p(counts)` while the real side used the stored, size-factor normalised
    `logcounts`. On the positive control, `crosscor_cosine` and
    `crosscor_mantel` now come out at exactly 1.
  - `scdesign2`: order the input by `spatial_cluster` before simulating.
    `simulate_count_scDesign2()` returns cells grouped per cell type, so the
    coordinates taken from the unsorted input belonged to different spots than
    the counts they were attached to. Also make `cell_type_sel` and
    `cell_type_prop` agree in order, which they did not for 10 or more clusters.
    The simulated spots are returned in the order they came in.
  - `splatter`, `sparsim`, `symsim` and `zinbwave`: return the simulated spots
    in the order they were given, rather than sorted by `spatial_cluster`.
    `clustering_ari` and `clustering_nmi` compare the real and simulated
    clusterings element-wise, so those four were being scored against
    mismatched spots.
  - `generate_sim_spatialcluster`: check that the simulated dataset still lines
    up with the real one, through a new `check_alignment()` helper.
  - `downstream`: reuse the `spatial_cluster` that
    `process_datasets/generate_sim_spatialcluster` already computed, instead of
    running `BayesSpace::spatialCluster()` a second time. Saves an MCMC run per
    method per dataset and makes the metric report the clustering the pipeline
    stored. This does not make `clustering_ari` reproducible: both before and
    after, a fixed reference is compared against one unseeded run.
  - `generate_cosine()`: filter the Moran's I values rather than the objects
    holding them, so that a single NaN no longer takes the whole metric with
    it. `crosscor_cosine` was NA for 32 of 99 runs.
  - `calculate_precision()`: report 0 rather than NA when a simulation has no
    spatially variable genes at all. Both negative controls were left without a
    score on any dataset, so nothing anchored the bottom of the scale.
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
