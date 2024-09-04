suppressMessages(library(SingleCellExperiment, quietly = TRUE))

## VIASH START
base <- ".cache/openproblems/figshare/SpatialSimBench_dataset/26054188/PDAC"
par <- list(
  # inputs
  input_sc = paste0(base, "_sc.rds"),
  input_sp = paste0(base, ".rds"),

  # outputs
  output_sc = "resources_test/spatialsimbench_mobnew/MOBNEW_sc.rds",
  output_sp = "resources_test/spatialsimbench_mobnew/MOBNEW.rds",

  # dataset metadata
  dataset_id = "MOBNEW",
  dataset_name = "MOBNEW",
  dataset_description_spatial = "MOBNEW",
  dataset_description_singlecell = "MOBNEW",
  dataset_url_spatial = "...",
  dataset_url_singlecell = "...",
  dataset_reference_singlecell = "...",
  dataset_reference_spatial = "...",
  dataset_summary_singlecell = "...",
  dataset_summary_spatial = "...",
  dataset_organism = "...",
  dataset_assay_spatial = "...",
  dataset_assay_singlecell = "..."
)
## VIASH END

process_matrix <- function(obj, layer_name) {
  if (!layer_name %in% assayNames(obj)) {
    return(NULL)
  }
  assay(obj, layer_name) |>
    as("CsparseMatrix") |>
    Matrix::t()
}

cat("Read input files\n")
input_sc <- readRDS(par$input_sc)
input_sp <- readRDS(par$input_sp)

cat("Single cell dataset:\n")
print(input_sc)

cat("Spatial dataset:\n")
print(input_sp)

cat("Construct uns\n")
uns <- list(
  dataset_id = par$dataset_id,
  dataset_name = par$dataset_name,
  dataset_organism = par$dataset_organism
)

cat("Transforming single cell into AnnData\n")
output_sc <- anndata::AnnData(
  layers = list(
    counts = process_matrix(input_sc, "counts")
  ),
  obs = data.frame(
    row.names = colnames(input_sc),
    cell_type = colData(input_sc)$cellType,
    donor_id = colData(input_sc)$sampleInfo
  ),
  var = data.frame(
    row.names = rownames(input_sc),
    feature_id = rownames(input_sc),
    feature_name = rownames(input_sc)
  ),
  uns = c(
    uns,
    list(
      dataset_summary = par$dataset_summary_singlecell,
      dataset_description = par$dataset_description_singlecell,
      dataset_url = par$dataset_url_singlecell,
      dataset_reference = par$dataset_reference_singlecell,
      dataset_assay = par$dataset_assay_singlecell
    )
  )
)

cat("Transforming spatial into AnnData\n")
celltype_proportions <- as.data.frame(metadata(input_sp)[["celltype_prop"]])

output_sp <- anndata::AnnData(
  layers = list(
    counts = process_matrix(input_sp, "counts"),
    logcounts = process_matrix(input_sp, "logcounts")
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
    celltype_proportions = celltype_proportions
  ),
  uns = c(
    uns,
    list(
      dataset_summary = par$dataset_summary_spatial,
      dataset_description = par$dataset_description_spatial,
      dataset_url = par$dataset_url_spatial,
      dataset_reference = par$dataset_reference_spatial,
      dataset_assay = par$dataset_assay_spatial
    )
  )
)

cat("Write output files\n")
output_sc$write_h5ad(par$output_sc, compression = "gzip")
output_sp$write_h5ad(par$output_sp, compression = "gzip")
