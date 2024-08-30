suppressMessages(library(SingleCellExperiment, quietly = TRUE))

## VIASH START
par <- list(
  # inputs
  input_sc = "resources_test/datasets_raw/MOBNEW/dataset_sc.rds",
  input_sp = "resources_test/datasets_raw/MOBNEW/dataset_sp.rds",

  # outputs
  output_sc = "resources_test/spatialsimbench_mobnew/MOBNEW_sc.rds",
  output_sp = "resources_test/spatialsimbench_mobnew/MOBNEW.rds",

  # dataset metadata
  dataset_id = "MOBNEW",
  dataset_name = "MOBNEW",
  dataset_description = "MOBNEW",
  dataset_url = "...",
  dataset_reference = "...",
  dataset_summary = "...",
  dataset_organism = "...",
  dataset_assay_spatial = "...",
  dataset_assay_singlecell = "..."
)
## VIASH END

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
  dataset_description = par$dataset_description,
  dataset_url = par$dataset_url,
  dataset_reference = par$dataset_reference,
  dataset_summary = par$dataset_summary,
  dataset_organism = par$dataset_organism,
  dataset_assay_spatial = par$dataset_assay_spatial,
  dataset_assay_singlecell = par$dataset_assay_singlecell
)

cat("Transforming single cell into AnnData\n")
output_sc <- anndata::AnnData(
  layers = list(
    counts = Matrix::t(assay(input_sc, "counts")),
    logcounts = Matrix::t(assay(input_sc, "logcounts"))
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
  uns = uns
)

cat("Transforming spatial into AnnData\n")
celltype_proportions <- as.data.frame(metadata(input_sp)[["celltype_prop"]])

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
    celltype_proportions = celltype_proportions
  ),
  uns = uns
)

cat("Write output files\n")
output_sc$write_h5ad(par$output_sc, compression = "gzip")
output_sp$write_h5ad(par$output_sp, compression = "gzip")
