## VIASH START
par <- list(
  input = "resources_test/spatialsimbench_mobnew/dataset_sp.h5ad",
  output = "simulated_dataset.h5ad"
)
meta <- list(
  name = "negative"
)
## VIASH END

cat("Reading input files\n")
input <- anndata::read_h5ad(par$input)

count_matrix <- as.matrix(input$layers[['counts']])

shuffled_values <- sample(as.vector(count_matrix))

shuffled_matrix <- matrix(shuffled_values, nrow = nrow(count_matrix), ncol = ncol(count_matrix))


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
