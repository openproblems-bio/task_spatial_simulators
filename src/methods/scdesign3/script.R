suppressMessages(library(SingleCellExperiment, quietly = TRUE))
suppressMessages(library(scDesign3, quietly = TRUE))

## VIASH START
par <- list(
  input = "resources_test/datasets/MOBNEW/dataset_sp.h5ad",
  base = "domain",
  output = "simulated_dataset.h5ad",
  family = "nb",
  usebam = FALSE
)
meta <- list(
  name = "scdesign3"
)
## VIASH END

cat("Read input files\n")
input <- anndata::read_h5ad(par$input)

sce <- SingleCellExperiment(
  list(counts = Matrix::t(input$layers[["counts"]])),
  colData = input$obs
)

cat("scDesign3 Simulation Start\n")

if (par$base == "tissue") {
  sce@colData$cell_type <- "cell_type"
} else if (par$base == "domain") {
  sce@colData$cell_type <- as.character(sce@colData$spatial_cluster)
} else {
  stop("wrong base parameter")
}

sce_simu <- scdesign3(
  sce = sce,
  assay_use = "counts",
  celltype = "cell_type",
  pseudotime = NULL,
  spatial = c("row", "col"), # spatial location
  other_covariates = NULL,
  mu_formula = "s(row, col, bs = 'gp', k = 100)",
  sigma_formula = "1",
  family_use = par$family, # could be change another distribution
  n_cores = 1,
  usebam = par$usebam,
  corr_formula = "1",
  copula = "gaussian",
  DT = TRUE,
  pseudo_obs = FALSE,
  return_model = FALSE,
  nonzerovar = FALSE,
  parallelization = "mcmapply"
)

cat("Generating output file\n")
new_obs <- sce_simu$new_covariate[c("row", "col")]

output <- anndata::AnnData(
  layers = list(
    counts = Matrix::t(sce_simu$new_count)
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
