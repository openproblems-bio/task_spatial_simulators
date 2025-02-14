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

# generate random values
n_rows <- nrow(input)
n_cols <- ncol(input)

values <- rnorm(n = n_rows * n_cols, mean = 3, sd = 1)

# make sure all values are positive
values[values < 0] <- abs(values[values < 0])

cat("Generate outoput file\n")
output <- anndata::AnnData(
  layers = list(
    counts = matrix(values, nrow = n_rows, ncol = n_cols)
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
