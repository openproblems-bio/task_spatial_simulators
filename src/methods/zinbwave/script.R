suppressMessages(library(SingleCellExperiment, quietly = TRUE))
suppressMessages(library(zinbwave, quietly = TRUE))

## VIASH START
par <- list(
  input = "resources_test/datasets/MOBNEW/dataset_sp.h5ad",
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

ordered_indices <- order(colData(sce)$spatial.cluster)
sce_ordered <- real_sce[, ordered_indices]

cat("ZINB-WaVE simulation start\n")

if (par$base != "domain") {
  stop("ONLY domain base")
}

multicoreParam <- MulticoreParam(workers = 8)

X = model.matrix(~spatial.cluster, data=colData(real_sce_ordered))
params <- zinbEstimate(as.matrix(counts(real_sce_ordered)), design.samples=X, BPPARAM = multicoreParam)
simulated_result <- zinbSimulate(params)

colnames(simulated_result) <- colnames(real_sce_ordered)
rownames(simulated_result) <- rownames(real_sce_ordered)

simulated_result_order <- real_sce_ordered
counts(simulated_result_order) <- counts(simulated_result)
  
simulated_result_order <- simulated_result_order[,match(colnames(real_sce), colnames(simulated_result_order))]
simulated_result_order <- simulated_result_order[match(rownames(real_sce), rownames(simulated_result_order)),]

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
