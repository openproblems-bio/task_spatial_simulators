suppressMessages(library(SRTsim, quietly = TRUE))

## VIASH START
par <- list(
  input = "resources_test/datasets/MOBNEW/dataset_sp.h5ad",
  output = "simulated_dataset.h5ad",
  base = "domain"
)
meta <- list(
  name = "srtsim"
)
## VIASH END

cat("Reading input files\n")
input <- anndata::read_h5ad(par$input)

cat("SRTsim simulation start\n")

real_count <- Matrix::t(input$layers["counts"])
real_loc <- data.frame(x = input$obs["row"], y = input$obs["col"], region = input$obs["spatial_cluster"])
rownames(real_loc) <- rownames(input$obs)

simSRT <- createSRT(count_in = real_count, loc_in = real_loc)

if (par$base == "domain") {
  simSRT1 <- srtsim_fit(simSRT, sim_schem = "domain")
} else if (par$base == "tissue") {
  simSRT1 <- srtsim_fit(simSRT, sim_schem = "tissue")
} else {
  stop("wrong base parameter")
}

simSRT1 <- srtsim_count(simSRT1)
counts_single <- as.matrix(simSRT1@simCounts)

col_data <- data.frame(
  row = data.frame(simSRT1@simcolData)$x,
  col = data.frame(simSRT1@simcolData)$y,
  row.names = rownames(data.frame(simSRT1@simcolData))
)

cat("Generating output\n")

output <- anndata::AnnData(
  layers = list(
    counts = Matrix::t(counts_single)
  ),
  obs = col_data,
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
