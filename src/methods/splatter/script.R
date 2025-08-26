suppressMessages(library(SingleCellExperiment, quietly = TRUE))
suppressMessages(library(splatter, quietly = TRUE))

## VIASH START
par <- list(
  input = "resources_test/spatialsimbench_mobnew/dataset_sp.h5ad",
  base = "domain"
)
meta <- list(
  name = "splatter"
)
## VIASH END

cat("Reading input files\n")
input <- anndata::read_h5ad(par$input)

cat("Splatter simulation start\n")

if (par$base != "domain") {
  stop("ONLY domain base")
}

ordered_indices <- order(input$obs$spatial_cluster)
input_ordered <- input[ordered_indices]

simulated_result <- NULL
for (spatial_cluster in unique(input_ordered$obs[["spatial_cluster"]])) {
  res <- try({
    input_spatial_cluster <- input_ordered[input_ordered$obs[["spatial_cluster"]] == spatial_cluster]
    params <- splatter::splatEstimate(as.matrix(t(input_spatial_cluster$layers[["counts"]])))
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

colnames(simulated_result) <- rownames(input_ordered$obs)
rownames(simulated_result) <- rownames(input_ordered$var)

simulated_result_ordered <- counts(simulated_result)[
  match(rownames(counts(simulated_result)), rownames(input_ordered$var)),
  match(colnames(counts(simulated_result)), rownames(input_ordered$obs))
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

cat("Write output files\n")
output$write_h5ad(par$output, compression = "gzip")
