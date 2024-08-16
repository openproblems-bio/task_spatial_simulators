requireNamespace("anndata", quietly = TRUE)
requireNamespace("Matrix", quietly = TRUE)
requireNamespace("SPARK", quietly = TRUE)
requireNamespace("spots", quietly = TRUE)

## VIASH START
par <- list(
  # inputs
  input_sp = "resources_test/datasets/MOBNEW/temp_dataset_sp_part2.h5ad",
  # outputs
  output_sp = "resources_test/datasets/MOBNEW/dataset_sp.h5ad"
)
meta <- list(
  resources_dir = "src/helpers"
)
## VIASH END

source(file.path(meta$resources_dir, "utils.R"))

cat("Read input files\n")
input_sp <- anndata::read_h5ad(par$input_sp)

cat("Spatial dataset:\n")
print(input_sp)

cat("Computing SVGs with SPARK\n")
sparkx_out <- generate_svg_sparkx(input_sp)
input_sp$varm$spatial_variable_genes <- sparkx_out$res_mtest

cat("Computing spatial autocorrelation with moransI\n")
moransI_out <- generate_moransI(input_sp)
input_sp$varm$spatial_autocorrelation <- as.data.frame(moransI_out$Morans.I)

cat("Output:\n")
print(input_sp)

cat("Writing output to file\n")
input_sp$write_h5ad(par$output_sp, compression = "gzip")
