requireNamespace("anndata", quietly = TRUE)
requireNamespace("Matrix", quietly = TRUE)
requireNamespace("SPARK", quietly = TRUE)
requireNamespace("spots", quietly = TRUE)

## VIASH START
par <- list(
  # inputs
  input_sp = "resources_test/datasets/MOBNEW/dataset_sp.h5ad",
  # outputs
  output_sp = "resources_test/datasets/MOBNEW/dataset_sp.h5ad"
)
## VIASH END

cat("Read input files\n")
input_sp <- anndata::read_h5ad(par$input_sp)

cat("Spatial dataset:\n")
print(input_sp)

cat("Run sparkx SVG\n")
generate_svg_sparkx <- function(adata) {
  # get count matrix
  sp_count <- Matrix::t(adata$layers[["counts"]])
  
  # format location as a matrix
  location <- as.matrix(adata$obs[, c("col", "row")])
  rownames(location) <- colnames(sp_count)

  # remove mitochondrial genes
  sp_count <- sp_count[!grepl("^(MT|mt)-", rownames(sp_count)), ]

  # run sparkx
  sparkX <- SPARK::sparkx(sp_count, location, numCores = 1, option = "mixture")
  
  return(sparkX)
}

sparkx_out <- generate_svg_sparkx(input_sp)
input_sp$obsm$spatial_variable_genes <- sparkx_out$res_mtest

cat("spatial autocorrelation task:\n")
generate_moransI <- function(adata){
  # get count matrix
  sp_count <- adata$layers[["logcounts"]]

  # get location as a matrix
  loc <- as.matrix(adata$obs[, c("row", "col")])

  # compute inverse distance matrix
  weights <- 1 / as.matrix(dist(loc))
  diag(weights) <- 0

  # run moransI
  res <- spots::BivariateMoransI(X =  sp_count, W = weights)

  return(res)
}

moransI_out <- generate_moransI(input_sp)
input_sp$varm$spatial_autocorrelation <- moransI_out$Morans.I

cat("Output:\n")
print(input_sp)

cat("Writing output to file\n")
input_sp$write_h5ad(par$output_sp, compression = "gzip")
