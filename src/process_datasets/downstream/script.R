suppressMessages(library(SingleCellExperiment, quietly = TRUE))

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

cat("SVG task:\n")
generate_svg_sparkx <- function(adata) {
  sp_count <- Matrix::t(adata$layers[["counts"]])
  
  info <- cbind.data.frame(x=as.numeric(sapply(strsplit(colnames(sp_count),split="x"),"[",1)),
                          y=as.numeric(sapply(strsplit(colnames(sp_count),split="x"),"[",2)))
  
  rownames(info) <- colnames(sp_count)
  
  location <- as.matrix(info)
  mt_idx <- grep("mt-", rownames(sp_count))
  if (length(mt_idx) != 0) {
    sp_count <- sp_count[-mt_idx, ]
  }
  sparkX <- SPARK::sparkx(sp_count, location, numCores = 1, option = "mixture")
  return(sparkX)
}

svg_result <- generate_svg_sparkx(input_sp)$res_mtest

cat("spatial autocorrelation task:\n")
generate_moransI <- function(adata){
  sp_count <- adata$layers[["logcounts"]]
  
  W <- cbind.data.frame(row = as.numeric(sapply(strsplit(rownames(sp_count),split="x"),"[",1)),
                        col = as.numeric(sapply(strsplit(rownames(sp_count),split="x"),"[",2)))
  W <- 1/as.matrix(dist(W))
  diag(W) <- 0
  res <- spots::BivariateMoransI( X =  sp_count , W =  W)
  return(res$Morans.I)
}

moransI_result <- generate_moransI(input_sp)

output_sp <- anndata::AnnData(
  layers = list(
    counts = Matrix::t(assay(input_sp, "counts")),
    logcounts = Matrix::t(assay(input_sp, "logcounts"))
  ),
  obs = data.frame(
    row.names = colnames(input_sp),
    col = colData(input_sp)$col,
    row = colData(input_sp)$row,
    sizeFactor = colData(input_sp)$sizeFactor,
    spatial_cluster = colData(input_sp)$spatial.cluster
  ),
  var = data.frame(
    row.names = rownames(input_sp),
    feature_id = rownames(input_sp),
    feature_name = rownames(input_sp)
  ),
  obsm = list(
    celltype_proportions = celltype_proportions,
    L_stats = scfeatures_result$L_stats,
    celltype_interaction = scfeatures_result$celltype_interaction,
    nn_correlation = scfeatures_result$nn_correlation,
    morans_I = scfeatures_result$morans_I,
    spatial_variable_genes = svg_result,
    spatial_autocorrelation = moransI_result

  ),
  uns = list(
    dataset_id = par$dataset_id,
    dataset_name = par$dataset_name,
    dataset_description = par$dataset_description,
    dataset_url = par$dataset_url,
    dataset_reference = par$dataset_reference,
    dataset_summary = par$dataset_summary,
    dataset_organism = par$dataset_organism
  )
)

cat("Write output files\n")
output_sp$write_h5ad(par$output_sp, compression = "gzip")
