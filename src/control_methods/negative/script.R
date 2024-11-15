## VIASH START
par <- list(
  input = "temp_olfactorybulb.negative.ks_statistic_gene_cell/_viash_par/input_spatial_dataset_1/output_sp.h5ad",
  output = "simulated_dataset.h5ad"
)
meta <- list(
  name = "negative"
)
## VIASH END

cat("Reading input files\n")
input <- anndata::read_h5ad(par$input)

count_matrix <- as.matrix(input$layers[['counts']])

# Calculate the total number of elements
total_elements <- length(count_matrix)

# Count the number of zeros
num_zeros <- sum(count_matrix == 0)

# Calculate the proportion of zeros
proportion_zeros <- num_zeros / total_elements

# Print the result
proportion_zeros

n_rows <- nrow(count_matrix)
n_cols <- ncol(count_matrix)

random_matrix_uniform <- matrix(runif(n_rows * n_cols), nrow = n_rows, ncol = n_cols)
rownames(random_matrix_uniform) <- rownames(count_matrix)

# Calculate the number of elements to set to zero
num_elements <- length(random_matrix_uniform)
num_zeros <- ceiling(proportion_zeros * num_elements)

# Randomly select indices to set to zero
indices <- sample(seq_len(num_elements), num_zeros)

# Set the selected indices to zero
random_matrix_uniform[indices] <- 0


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
