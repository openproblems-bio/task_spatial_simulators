## VIASH START
par <- list(
  input = "resources_test/datasets/MOBNEW/dataset_sp.h5ad",
  output = "simulated_dataset.h5ad"
)
meta <- list(
  name = "positive"
)
## VIASH END

cat("Reading input files\n")
input <- anndata::read_h5ad(par$input)

input$uns$method_id <- meta$name

cat("Write output files\n")
input$write_h5ad(par$output, compression = "gzip")