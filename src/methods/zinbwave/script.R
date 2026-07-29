if (!requireNamespace("zinbwave", quietly = TRUE)) {
  stop("The 'zinbwave' package is required but is not installed.")
}
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
input <- anndataR::read_h5ad(par$input)

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

simulated_result_ordered <- counts(simulated_result)

# put the spots back in the order they came in. The metrics compare the
# simulated dataset against the real one spot for spot -- ARI and NMI are
# invariant to a permutation of the labels, but not of the samples.
restore_order <- order(ordered_indices)

cat("Generating output\n")
output <- anndataR::AnnData(
  layers = list(
    counts = Matrix::t(simulated_result_ordered)[restore_order, , drop = FALSE]
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
