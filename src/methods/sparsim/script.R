suppressMessages(library(SingleCellExperiment, quietly = TRUE))
suppressMessages(library(SPARSim, quietly = TRUE))

## VIASH START
par <- list(
  input = "resources_test/spatialsimbench_mobnew/MOBNEW.rds",
  base = "domain"
)
meta <- list(
  name = "sparsim"
)
## VIASH END

find_cluster_indices <- function(cluster_column) {
  unique_clusters <- sort(unique(cluster_column))
  conditions <- list()

  for (i in seq_along(unique_clusters)) {
    cluster <- unique_clusters[i]
    indices <- which(cluster_column == cluster)
    range_name <- sprintf("cluster_%s_column_index", LETTERS[i])
    assign(range_name, c(min(indices):max(indices)), envir = .GlobalEnv)
    conditions[[range_name]] <- get(range_name)
  }

  return(conditions)
}

cat("Reading input files\n")
input <- anndata::read_h5ad(par$input)

sce <- SingleCellExperiment(
  list(counts = Matrix::t(input$layers[["counts"]])),
  colData = input$obs
)

cat("SPARsim simulation start\n")

ordered_indices <- order(colData(sce)$spatial_cluster)
sce_ordered <- sce[, ordered_indices]

if (par$base != "domain") {
  stop("ONLY domain base")
}

count_matrix <- data.frame(assay(sce_ordered))
sce_scran <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = as.matrix(count_matrix)))
sce_scran <- scran::computeSumFactors(sce_scran, sizes = seq(20, 100, 5), positive = F) 

if (any(sce_scran$sizeFactor <= 0)) {
  threshold <- 1e-10
  sce_scran$sizeFactor[sce_scran$sizeFactor <= 0] <- threshold
}
count_matrix_norm <- scater::normalizeCounts(sce_scran, log = FALSE)
count_matrix_conditions <- find_cluster_indices(sce_ordered@colData$spatial_cluster)
SPARSim_sim_param <- SPARSim_estimate_parameter_from_data(
  raw_data = count_matrix,
  norm_data = count_matrix_norm,
  conditions = count_matrix_conditions
)

sim_result <- SPARSim_simulation(dataset_parameter = SPARSim_sim_param)
colnames(sim_result$count_matrix) <- gsub("\\.", "-", colnames(sim_result$count_matrix))
simulated_result_order <- sce_ordered
assays(simulated_result_order, withDimnames = FALSE) <- list(counts = sim_result$count_matrix)
simulated_result_order <- simulated_result_order[,match(colnames(sce), colnames(simulated_result_order))]
simulated_result_order <- simulated_result_order[match(rownames(sce), rownames(simulated_result_order)),]
new_obs <- as.data.frame(simulated_result_order@colData[c("row", "col")])

cat("Generating output\n")
output <- anndata::AnnData(
  layers = list(
    counts = Matrix::t(counts(simulated_result_order))
  ),
  obs = new_obs,
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
