

## VIASH START
par <- list(
  # inputs
  input_sp = "resources_test/datasets/MOBNEW/dataset_sp.h5ad",
  input_sp_sim = "resources_test/datasets/MOBNEW/simulated_dataset.h5ad",
  # outputs
  output_sp = "resources_test/datasets/MOBNEW/simulated_dataset.h5ad"
)
meta <- list(
  resources_dir = "target/executable/process_datasets/generate_sim_sparialCluster"
)
## VIASH END
source(file.path(meta$resources_dir, "utils.R"))

cat("Read input files\n")
input_real_sp <- anndata::read_h5ad(par$input_sp)
input_simulated_sp <- anndata::read_h5ad(par$input_sp_sim)

sim_sce <- scater::logNormCounts(SingleCellExperiment::SingleCellExperiment(
  list(counts = Matrix::t(input_simulated_sp$layers[["counts"]])),
  colData = input_simulated_sp$obs,
  metadata = input_simulated_sp$obsm
))
cat("add spatial cluster in simulated dataset:\n")
sim_cluster <- generate_sim_spatialCluster(input_real_sp, input_simulated_sp)

# need reclassify again
real_cluster <- input_real_sp$obs[,c("spatial_cluster")]
location <- colnames(counts(sim_sce))
sim_new_cluster <- reclassify_simsce(location, real_cluster, sim_cluster)

input_simulated_sp$obs$spatial_cluster <- sim_new_cluster

cat("Writing output to file\n")
input_simulated_sp$write_h5ad(par$output_sp, compression = "gzip")
