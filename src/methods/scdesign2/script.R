requireNamespace("scDesign2", quietly = TRUE)

## VIASH START
par <- list(
  input = "resources_test/datasets/MOBNEW/dataset_sp.h5ad",
  output = "simulated_dataset.h5ad",
  base = "domain"
)
meta <- list(
  name = "scDesign2",
  cpus = 8L
)
## VIASH END

cat("Reading input files\n")
input <- anndata::read_h5ad(par$input)

if (par$base != "domain") {
  stop("ONLY domain base")
}

cat("scDesign2 simulation start\n")
counts <- as.matrix(Matrix::t(input$layers[["counts"]]))
colnames(counts) <- as.character(input$obs$spatial_cluster)

sim_out <- scDesign2::fit_model_scDesign2(
  data_mat = counts,
  cell_type_sel = unique(colnames(counts)),
  sim_method = "copula",
  ncores = meta$cpus
)
  
sim_out <- simulate_count_scDesign2(copula_result, ncol(counts),sim_method = 'copula',cell_type_prop = prop.table( spatial_cluster_prop))

new_obs <- sce_simu$new_covariate
remap <- c(
  row = "row",
  col = "col"
)
colnames(new_obs) <- remap[colnames(new_obs)]

cat("Generating output\n")

output <- anndata::AnnData(
  layers = list(
    counts = Matrix::t(out)
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
