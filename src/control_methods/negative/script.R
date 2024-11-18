## VIASH START
par <- list(
  input = "resources_test/datasets/MOBNEW/dataset_sp.h5ad",
  output = "simulated_dataset.h5ad"
)
meta <- list(
  name = "negative"
)
## VIASH END

cat("Reading input files\n")
input <- anndata::read_h5ad(par$input)

count_matrix <- as.matrix(input$layers[['counts']])


n_rows <- nrow(count_matrix)
n_cols <- ncol(count_matrix)

random_matrix_uniform <- matrix(runif(n_rows * n_cols), nrow = n_rows, ncol = n_cols)
rownames(random_matrix_uniform) <- rownames(count_matrix)

cat("Generate outoput file\n")
output <- anndata::AnnData(
  layers = list(
    counts = random_matrix_uniform
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
