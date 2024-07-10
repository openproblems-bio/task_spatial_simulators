suppressMessages(library(SingleCellExperiment, quietly = TRUE))
suppressMessages(library(scDesign2, quietly = TRUE))

## VIASH START
par <- list(
  input = "resources_test/datasets/MOBNEW/dataset_sp.h5ad",
  output = "simulated_dataset.h5ad",
  base = "domain"
)
meta <- list(
  name = "scdeisgn2"
)
## VIASH END

cat("Reading input files\n")
input <- anndata::read_h5ad(par$input)

if (par$base != "domain") {
  stop("ONLY domain base")
}

cat("scDesign2 simulation start\n")

traincount <- as.matrix(Matrix::t(input$layers[["counts"]]))
spatial_cluster <- adata$obs$spatial_cluster
spatial_cluster <- as.character(spatial_cluster)
spatial_cluster_sel <- unique(spatial_cluster)
colnames(traincount) <- spatial_cluster
spatial_cluster_prop <- table(spatial_cluster)
copula_result <- fit_model_scDesign2(traincount, spatial_cluster_sel, sim_method = 'copula', ncores = 8)
  
sim_count_copula <- simulate_count_scDesign2(copula_result, ncol(traincount),sim_method = 'copula',cell_type_prop = prop.table( spatial_cluster_prop))

new_obs <- sce_simu$new_covariate
remap <- c(
  row = "row",
  col = "col"
)
colnames(new_obs) <- remap[colnames(new_obs)]

cat("Generating output\n")

output <- anndata::AnnData(
  layers = list(
    counts = Matrix::t(sim_count_copula)
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
