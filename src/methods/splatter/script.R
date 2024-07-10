suppressMessages(library(SingleCellExperiment, quietly = TRUE))
suppressMessages(library(splatter, quietly = TRUE))

## VIASH START
par <- list(
  input = "resources_test/datasets/MOBNEW/dataset_sp.h5ad",
  base = "domain"
)
meta <- list(
  name = "splatter"
)
## VIASH END

cat("Reading input files\n")
input <- anndata::read_h5ad(par$input)

sce <- SingleCellExperiment(
  list(counts = Matrix::t(input$layers[["counts"]])),
  colData = input$obs
)

cat("Splatter simulation start\n")

if (par$base != "domain") {
  stop("ONLY domain base")
}

ordered_indices <- order(colData(sce)$spatial_cluster)
sce_ordered <- sce[, ordered_indices]

simulated_result <- NULL
for (spatial_cluster in (unique(sce_ordered$spatial_cluster))) {
  print(spatial_cluster)
  res <- try({
    sce_spatial_cluster <- sce_ordered[, sce_ordered$spatial_cluster == spatial_cluster]
    params <- splatter::splatEstimate(as.matrix(counts(sce_spatial_cluster)))
    sim_spatial_cluster <- splatter::splatSimulate(params)
    sim_spatial_cluster$spatial_cluster <- spatial_cluster
    colnames(sim_spatial_cluster) <- paste0(spatial_cluster, colnames(sim_spatial_cluster))
    names(rowData(sim_spatial_cluster)) <- paste(spatial_cluster, names(rowData(sim_spatial_cluster)))

    # combine the cell types
    if (is.null(simulated_result)) {
      simulated_result <- sim_spatial_cluster
    } else {
      simulated_result <- SingleCellExperiment::cbind(simulated_result, sim_spatial_cluster)
    }
  })
}

colnames(simulated_result) <- colnames(sce_ordered)
rownames(simulated_result) <- rownames(sce_ordered)

cat("Generating output\n")

simulated_result_order <- sce_ordered
counts(simulated_result_order) <- counts(simulated_result)

simulated_result_order <- simulated_result_order[, match(colnames(sce), colnames(simulated_result_order))]
simulated_result_order <- simulated_result_order[match(rownames(sce), rownames(simulated_result_order)), ]

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
