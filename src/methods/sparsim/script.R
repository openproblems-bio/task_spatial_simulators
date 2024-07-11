suppressMessages(library(SingleCellExperiment, quietly = TRUE))
suppressMessages(library(SPARSim, quietly = TRUE))
suppressMessages(library(Seurat, quietly = TRUE))

## VIASH START
par <- list(
  input = "resources_test/datasets/MOBNEW/dataset_sp.h5ad",
  base = "domain"
)
meta <- list(
  name = "sparsim"
)
## VIASH END

cat("Reading input files\n")
input <- anndata::read_h5ad(par$input)

sce <- SingleCellExperiment(
  list(counts = Matrix::t(input$layers[["counts"]])),
  colData = input$obs
)

cat("SPARsim simulation start\n")

ordered_indices <- order(colData(sce)$spatial_cluster)
sce_ordered <- sce[, ordered_indices]
  
if (base == "domain"){
  message("SPARsim Simulation Start")
}else{
  stop("ONLY domain base")
}
  
count_matrix <- data.frame(assay(sce_ordered)) 
sce_scran <- SingleCellExperiment::SingleCellExperiment(assays=list(counts=as.matrix(count_matrix)))
sce_scran <- scran::computeSumFactors(sce_scran, sizes=seq(20, 100, 5), positive=F) 
  
if (any(sce_scran$sizeFactor <= 0)) {
  threshold <- 1e-10
  sce_scran$sizeFactor[sce_scran$sizeFactor <= 0] <- threshold
}
count_matrix_norm <- scater::normalizeCounts(sce_scran, log = FALSE)
count_matrix_conditions <- find_cluster_indices(sce_ordered@colData$spatial_cluster) # nolint
SPARSim_sim_param <- SPARSim_estimate_parameter_from_data(raw_data = count_matrix,  # nolint
                                                          norm_data = count_matrix_norm,  # nolint
                                                          conditions = count_matrix_conditions) # nolint

sim_result <- SPARSim_simulation(dataset_parameter = SPARSim_sim_param)
colnames(sim_result$count_matrix) <- gsub("\\.", "-", colnames(sim_result$count_matrix))
simulated_result_order <- sce_ordered
assays(simulated_result_order, withDimnames = FALSE) <- list(counts = sim_result$count_matrix)
simulated_result_order <- simulated_result_order[,match(colnames(sce), colnames(simulated_result_order))]
simulated_result_order <- simulated_result_order[match(rownames(sce), rownames(simulated_result_order)),]


if (par$base != "domain") {
  stop("ONLY domain base")
}


cat("Generating output\n")


output <- anndata::AnnData(
  layers = list(
    counts = Matrix::t(counts(simulated_result_order))
  ),
  obs = as.data.frame(simulated_result_order@colData),
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
