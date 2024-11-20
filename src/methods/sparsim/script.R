suppressMessages(library(SingleCellExperiment, quietly = TRUE))
suppressMessages(library(SPARSim, quietly = TRUE))
suppressMessages(library(SummarizedExperiment, quietly = TRUE))

## VIASH START
par <- list(
  input = "resources_test/spatialsimbench_mobnew/dataset_sp.h5ad",
  base = "domain"
)
meta <- list(
  name = "sparsim"
)
## VIASH END

find_cluster_indices <- function(cluster_column) {
  unique_clusters <- sort(unique(cluster_column))
  conditions <- lapply(unique_clusters, function(cluster) {
    indices <- which(cluster_column == cluster)
    seq(min(indices), max(indices))
  })
  names(conditions) <- sprintf("cluster_%s_column_index", LETTERS[seq_along(unique_clusters)])
  return(conditions)
}

cat("Reading input files\n")
input <- anndata::read_h5ad(par$input)

cat("SPARSim simulation start\n")

if (par$base != "domain") {
  stop("Error: Only 'domain' base is supported.")
}

# Order by spatial cluster
ordered_indices <- order(input$obs$spatial_cluster)
input_ordered <- input[ordered_indices]

count_matrix <- as.matrix(Matrix::t(input_ordered$layers[["counts"]]))
sce_scran <- SingleCellExperiment::SingleCellExperiment(
  assays = list(counts = count_matrix),
  colData = input_ordered$obs
)
sce_scran <- scran::computeSumFactors(sce_scran, sizes = seq(20, 100, 5), positive = F) 

# Replace zero or negative size factors with threshold 1e-10
sce_scran$sizeFactor[sce_scran$sizeFactor <= 0] <- 1e-10

# Perform normalization
count_matrix_norm <- scater::normalizeCounts(sce_scran, log = FALSE)

# Find cluster indices for conditions
count_matrix_conditions <- find_cluster_indices(input_ordered$obs[["spatial_cluster"]])

# Estimate SPARSim parameters
SPARSim_sim_param <- SPARSim::SPARSim_estimate_parameter_from_data(
  raw_data = count_matrix,
  norm_data = count_matrix_norm,
  conditions = count_matrix_conditions
)

# Simulate new dataset
sim_result <- SPARSim::SPARSim_simulation(dataset_parameter = SPARSim_sim_param)

# Reorder simulated results
simulated_result_ordered <- sim_result$count_matrix[
  match(rownames(sim_result$count_matrix), rownames(input_ordered$var)),
  match(colnames(sim_result$count_matrix), rownames(input_ordered$obs))
]

cat("Generating output\n")
output <- anndata::AnnData(
  layers = list(counts = t(simulated_result_ordered)),
  obs = input_ordered$obs[c("row", "col")],
  var = input_ordered$var,
  uns = c(
    input$uns,
    list(
      method_id = meta$name
    )
  )
)

cat("Writing output files\n")
output$write_h5ad(par$output, compression = "gzip")