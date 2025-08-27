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
