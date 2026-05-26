suppressMessages(library(SingleCellExperiment, quietly = TRUE))
suppressMessages(library(splatter, quietly = TRUE))
suppressMessages(library(BiocParallel, quietly = TRUE))

## VIASH START
par <- list(
  input = "resources_test/spatialsimbench_mobnew/dataset_sp.h5ad",
  base = "domain"
)
meta <- list(
  name = "zinbwave"
)
## VIASH END

cat("Reading input files\n")
input <- anndata::read_h5ad(par$input)

ordered_indices <- order(input$obs$spatial_cluster)
input_ordered <- input[ordered_indices]

cat("ZINB-WaVE simulation start\n")

if (par$base != "domain") {
  stop("ONLY domain base")
}

cpus <- if (is.null(meta$cpus)) 2L else meta$cpus

multicoreParam <- MulticoreParam(workers = cpus)

X <- model.matrix(~spatial_cluster, data = input_ordered$obs)

params <- splatter::zinbEstimate(as.matrix(t(input_ordered$layers[["counts"]])), design.samples = X, BPPARAM = multicoreParam)
simulated_result <- splatter::zinbSimulate(params)

colnames(simulated_result) <- rownames(input_ordered$obs)
rownames(simulated_result) <- rownames(input_ordered$var)

simulated_result_ordered <- counts(simulated_result)[
  match(rownames(counts(simulated_result)), rownames(input_ordered$var)),
  match(colnames(counts(simulated_result)), rownames(input_ordered$obs))
]

cat("Generating output\n")
output <- anndata::AnnData(
  layers = list(
    counts = Matrix::t(simulated_result_ordered)
  ),
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
