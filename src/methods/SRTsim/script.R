suppressMessages(library(SingleCellExperiment, quietly = TRUE))
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

sce <- SingleCellExperiment(
  list(counts = Matrix::t(input$layers[["counts"]])),
  colData = input$obs
)

real_count <- as.matrix(Matrix::t(input$layers[["counts"]]))
real_loc <- data.frame(x = adata$obs$row,y = adata$obs$col, region = adata$obs$spatial_cluster)
  
simSRT<- createSRT(count_in=real_count,loc_in =real_loc)

cat("SRTsim simulation start\n")

if (base == "domain"){
    simSRT1 <- srtsim_fit(simSRT,sim_schem="domain")
  }else if (base == "tissue"){
    simSRT1 <- srtsim_fit(simSRT,sim_schem="tissue")
  }else{
    stop("wrong base parameter")
  }

simSRT1 <- srtsim_count(simSRT1)
counts_single <- as.matrix(simSRT1@simCounts)

new_obs <- sce_simu$new_covariate
remap <- c(
  cell_type = "spatial_cluster",
  row = "row",
  col = "col"
)
colnames(new_obs) <- remap[colnames(new_obs)]

cat("Generating output\n")

output <- anndata::AnnData(
  layers = list(
    counts = Matrix::t(simSRT1@simCounts)
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
