## VIASH START
par <- list(
  input = "resources_test/spatialsimbench_mobnew/dataset_sp.h5ad",
  output = "simulated_dataset.h5ad"
)
meta <- list(
  name = "negative_normal"
)
## VIASH END

cat("Reading input files\n")
input <- anndata::read_h5ad(par$input)

n_rows <- nrow(input)
n_cols <- ncol(input)

shuffled_values <- rnorm(n = n_rows * n_cols, mean = 10, sd = 5)

shuffled_matrix <- matrix(shuffled_values, nrow = n_rows, ncol = n_cols)

cat("Generate outoput file\n")
output <- anndata::AnnData(
  layers = list(
    counts = shuffled_matrix
  ),
  obs = input$obs[c("row", "col")],
  var = input$var,
  uns = c(
    input$uns,
    list(
      method_id = meta$name
    )
  )
)

cat("Write output files\n")
output$write_h5ad(par$output, compression = "gzip")
