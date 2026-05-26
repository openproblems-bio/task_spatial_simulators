---
description: "Use when writing, fixing, or reviewing method/metric script.py files in src/methods/, src/metrics/, or src/control_methods/. Covers script style, API compliance, and how to verify components."
applyTo: "src/methods/**/script.py,src/metrics/**/script.py,src/control_methods/**/script.py"
---
# Method & Metric Script Guidelines (Python)

## Core Principle

`script.py` should represent **typical bioinformatician usage** of the tool with minimal modifications. Only adapt what is strictly necessary to:
1. Read inputs from the paths provided by `par`
2. Pass the right data structures to the method
3. Convert the method's output back into the expected output structures
4. Write outputs to `par['output']`

Do **not** restructure the method's native API, add abstraction layers, or rewrite the algorithm logic.

## Finding API Specs

Input/output file formats are defined in `src/api/`. Key files:
- `file_dataset_sp.yaml` — spatial dataset input format (contains `layers['counts']`, `obs['row']`, `obs['col']`, `obs['spatial_cluster']`, etc.)
- `file_dataset_sc.yaml` — single-cell dataset input format (metrics only)
- `file_simulated_dataset.yaml` — expected output format for methods
- `file_score.yaml` — expected output format for metrics
- `comp_method.yaml`, `comp_control_method.yaml`, `comp_metric.yaml` — component argument specs

Always check these before deciding what fields to read or write.

## The `## VIASH START` / `## VIASH END` Block

This block is **auto-generated** by viash from the component's `config.vsh.yaml` arguments. It is replaced at build/test time with a real CLI parser. Keep it in the script only as a local development convenience.

- **Do not edit it manually** to add or remove parameters — edit `config.vsh.yaml` instead.
- After adding, removing, or renaming an argument in the config, regenerate the block:
  ```bash
  viash config inject src/methods/<name>/config.vsh.yaml
  ```
- Argument names in the config (`--my_param`) map directly to `par['my_param']` keys.

## Common Patterns

**Method: reading input:**
```python
import anndata as ad
input = ad.read_h5ad(par['input'])
```

**Method: writing simulated dataset output:**
```python
output = ad.AnnData(
    layers={"counts": simulated_counts},  # integer matrix, cells x genes
    obs=input.obs[["row", "col"]],
    var=input.var,
    uns={
        **input.uns,
        "method_id": meta["name"],
    },
)
output.write_h5ad(par['output'], compression="gzip")
```

**Metric: reading inputs:**
```python
import anndata as ad
input_spatial_dataset = ad.read_h5ad(par['input_spatial_dataset'])
input_singlecell_dataset = ad.read_h5ad(par['input_singlecell_dataset'])
input_simulated_dataset = ad.read_h5ad(par['input_simulated_dataset'])
```

**Metric: writing score output:**
```python
output = ad.AnnData(
    uns={
        "dataset_id": input_simulated_dataset.uns["dataset_id"],
        "method_id": input_simulated_dataset.uns["method_id"],
        "metric_ids": ["metric_name_1", "metric_name_2"],
        "metric_values": [score1, score2],
    },
)
output.write_h5ad(par['output'], compression="gzip")
```

## Dependency Fixes

If a library has a dependency conflict (e.g., incompatible with newer `anndata`, `numpy`, etc.), prefer replacing it with an alternative that provides the same model/algorithm natively rather than pinning transitive dependencies.

Update `config.vsh.yaml` to remove the broken package from the `setup` block when replacing it.

## Verification

After any change to a method script or config, verify with:
```bash
viash test src/methods/<name>/config.vsh.yaml
# or
viash test src/metrics/<name>/config.vsh.yaml
```

Both test scripts must succeed (`2 out of 2`).
