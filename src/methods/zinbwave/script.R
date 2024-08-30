suppressMessages(library(SingleCellExperiment, quietly = TRUE))
suppressMessages(library(splatter, quietly = TRUE))
suppressMessages(library(BiocParallel, quietly = TRUE))

## VIASH START
par <- list(
  input = "resources_test/spatialsimbench_mobnew/MOBNEW.rds",
  base = "domain"
)
meta <- list(
  name = "zinbwave"
)
## VIASH END

cat("Reading input files\n")
input <- anndata::read_h5ad(par$input)

sce <- SingleCellExperiment(
  list(counts = Matrix::t(input$layers[["counts"]])),
  colData = input$obs
)

ordered_indices <- order(colData(sce)$spatial_cluster)
sce_ordered <- sce[, ordered_indices]

cat("ZINB-WaVE simulation start\n")

if (par$base != "domain") {
  stop("ONLY domain base")
}

multicoreParam <- MulticoreParam(workers = 8)

X = model.matrix(~spatial_cluster, data=colData(sce_ordered))
params <- splatter::zinbEstimate(as.matrix(counts(sce_ordered)), design.samples = X, BPPARAM = multicoreParam)
simulated_result <- splatter::zinbSimulate(params)

colnames(simulated_result) <- colnames(sce_ordered)
rownames(simulated_result) <- rownames(sce_ordered)

simulated_result_order <- sce_ordered
counts(simulated_result_order) <- counts(simulated_result)
  
simulated_result_order <- simulated_result_order[,match(colnames(sce), colnames(simulated_result_order))]
simulated_result_order <- simulated_result_order[match(rownames(sce), rownames(simulated_result_order)),]
new_obs <- as.data.frame(simulated_result_order@colData[c("row", "col")])

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
