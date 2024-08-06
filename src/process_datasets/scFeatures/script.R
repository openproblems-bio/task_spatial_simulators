suppressMessages(library(SingleCellExperiment, quietly = TRUE))

## VIASH START
par <- list(
  # inputs
  input_sp = "resources_test/datasets/MOBNEW/dataset_sp.h5ad",
  sc_species = "Mus musculus",
  # outputs
  output_sp = "resources_test/datasets/MOBNEW/dataset_sp.h5ad"
)
## VIASH END

cat("Read input files\n")
input_sp <- anndata::read_h5ad(par$input_sp)

cat("Spatial dataset:\n")
print(input_sp)


cat("Run scFeatures\n")
sce <- scater::logNormCounts(SingleCellExperiment(
  list(counts = Matrix::t(input_sp$layers[["counts"]])),
  colData = input_sp$obs,
  metadata = input_sp$obsm
))

feat_types <- c("L_stats","celltype_interaction","nn_correlation","morans_I")
log_count <- assay(sce, "logcounts")

prob_matrix <- sce@metadata$celltype_prop
colnames(prob_matrix) <- paste0("ct", seq_len(ncol(prob_matrix)))
rownames(prob_matrix) <- colnames(log_count)

scfeatures <- scFeatures::scFeatures(log_count,
                                sample = rep("sample1", ncol(log_count)),
                                spatialCoords = list(colData(sce)$row,colData(sce)$col),
                                feature_types = feat_types,
                                type = "spatial_t",
                                species = sc_species,
                                spotProbability =  prob_matrix)
cat("Transforming spatial into AnnData\n")

# scFeatures

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
    L_stats = scFeatures$L_stats,
    celltype_interaction = scFeatures$celltype_interaction,
    nn_correlation = scFeatures$nn_correlation,
    morans_I = scFeatures$morans_I

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
