## VIASH START
par <- list(
  input_spatial_dataset = "resources_test/datasets/MOBNEW/dataset_sp.h5ad",
  input_singlecell_dataset = "resources_test/datasets/MOBNEW/dataset_sc.h5ad",
  input_simulated_dataset = "resources_test/datasets/MOBNEW/simulated_dataset.h5ad",
  output = "output.h5ad"
)
meta <- list(
  name = "downstream"
  # resources_dir = "target/executable/metrics/ks_statistics_scFeatures"
)

## VIASH END

input_real_sp <- anndata::read_h5ad(par$input_spatial_dataset)
input_sc <- anndata::read_h5ad(par$input_singlecell_dataset)
input_simulated_sp <- anndata::read_h5ad(par$input_simulated_dataset)

cat("spatial variable gene evaluation\n")
real_svg <- generate_svg_sparkx(input_real_sp)
sim_svg <- generate_svg_sparkx(input_simulated_sp)
svg_precision <- calculate_precision(real_svg, sim_svg)
svg_recall <- calculate_recall(real_svg, sim_svg)

cat("cell type deconvolution evaluation\n")
real_ct_prop <- CARD_processing(input_real_sp, input_sc)
sim_ct_prop <- CARD_processing(input_simulated_sp, input_sc)
ctdeconvolute_rmse <- generate_jds(real_ct_prop, sim_ct_prop)
ctdeconcolute_jsd <- generate_rmse(real_ct_prop, sim_ct_prop)

cat("spatial autocorrelation evaluation\n")
counts <- input_simulated_sp$layers[["counts"]]
logcounts <- log1p(counts)
input_simulated_sp$layers[["logcounts"]] <- logcounts
real_moransI <- generate_moransI(input_real_sp)
sim_moransI <- generate_moransI(input_simulated_sp)
crosscor_cosin <- generate_cosine(real_moransI, sim_moransI)
crosscor_mantel <- generate_mantel(real_moransI, sim_moransI)

cat("spatial clustering evaluation\n")
# TODO



cat("Combining metric values\n")
uns_metric_ids <- c(
  "svg_precision",
  "svg_recall",
  "ctdeconvolute_rmse",
  "ctdeconcolute_jsd",
  "crosscor_cosine",
  "crosscor_mantel"
)

uns_metric_values <- c(
  svg_precision,
  svg_recall,
  ctdeconvolute_rmse,
  ctdeconcolute_jsd,
  crosscor_cosine,
  crosscor_mantel
)

cat("Writing output AnnData to file\n")
output <- anndata::AnnData(
  uns = list(
    dataset_id = input_simulated_dataset$uns[["dataset_id"]],
    method_id = input_simulated_dataset$uns[["method_id"]],
    metric_ids = uns_metric_ids,
    metric_values = uns_metric_values
  ),
  shape = c(0L, 0L)
)
output$write_h5ad(par[["output"]], compression = "gzip")