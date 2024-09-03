requireNamespace("anndata", quietly = TRUE)
requireNamespace("scFeatures", quietly = TRUE)

## VIASH START
par <- list(
  # inputs
  input_sp = "resources_test/spatialsimbench_mobnew/dataset_sp.h5ad",
  # outputs
  output_sp = "output.h5ad"
)
## VIASH END

cat("Read input files\n")
input_sp <- anndata::read_h5ad(par$input_sp)

cat("Spatial dataset:\n")
print(input_sp)

cat("Fetch data\n")
# get organism
organism_mapping <- c(
  mus_musculus = "Mus musculus",
  homo_sapiens = "Homo sapiens"
)
species <- organism_mapping[[input_sp$uns[["dataset_organism"]]]]

cat("Run scFeatures\n")
scfeatures_result <- scFeatures::scFeatures(
  data = Matrix::t(input_sp$layers[["logcounts"]]),
  sample = rep("sample1", input_sp$n_obs),
  spatialCoords = input_sp$obs[, c("row", "col")],
  feature_types = c("L_stats", "celltype_interaction", "nn_correlation", "morans_I"),
  type = "spatial_t",
  species = species,
  spotProbability =  t(as.matrix(input_sp$obsm[["celltype_proportions"]]))
)

input_sp$uns$L_stats <- scfeatures_result$L_stats
input_sp$uns$celltype_interaction <- scfeatures_result$celltype_interaction
input_sp$uns$nn_correlation <- scfeatures_result$nn_correlation
input_sp$uns$morans_I <- scfeatures_result$morans_I

cat("Output dataset:\n")
print(input_sp)

cat("Write output files\n")
input_sp$write_h5ad(par$output_sp, compression = "gzip")
