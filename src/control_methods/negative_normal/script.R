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
input <- anndataR::read_h5ad(par$input)

# generate random values
n_rows <- nrow(input)
n_cols <- ncol(input)

values <- rnorm(n = n_rows * n_cols, mean = 3, sd = 1)

# make sure all values are positive
values[values < 0] <- abs(values[values < 0])

# counts are declared as integer in file_simulated_dataset.yaml, and edgeR,
# SPARK and scran all assume that. Handing them continuous values makes them
# fall over, which loses the run instead of scoring it badly.
values <- round(values)

cat("Generate outoput file\n")
output <- anndataR::AnnData(
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
