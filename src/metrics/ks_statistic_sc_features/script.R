requireNamespace("anndata", quietly = TRUE)
requireNamespace("edgeR", quietly = TRUE)
requireNamespace("ks", quietly = TRUE)
requireNamespace("resample", quietly = TRUE)
requireNamespace("reshape2", quietly = TRUE)
requireNamespace("scFeatures", quietly = TRUE)
requireNamespace("CARD", quietly = TRUE)
requireNamespace("SingleCellExperiment", quietly = TRUE)
requireNamespace("SummarizedExperiment", quietly = TRUE)


packages <- c(
  "SingleCellExperiment",
  "SummarizedExperiment",
  "concaveman",
  "sp",
  "Matrix",
  "methods",
  "ggplot2",
  "ggcorrplot",
  "MuSiC",
  "fields",
  "MCMCpack",
  "dplyr",
  "sf",
  "RANN",
  "stats",
  "reshape2",
  "RColorBrewer",
  "scatterpie",
  "grDevices",
  "nnls",
  "pbmcapply",
  "spatstat",
  "gtools",
  "RcppML",
  "NMF"
)

# Load each package
lapply(packages, library, character.only = TRUE)

## VIASH START
par <- list(
  input_spatial_dataset = "temp_brain.zinbwave.ks_statistic_sc_features/_viash_par/input_spatial_dataset_1/output_sp.h5ad",
  input_singlecell_dataset = "temp_brain.zinbwave.ks_statistic_sc_features/_viash_par/input_singlecell_dataset_1/output_sc.h5ad",
  input_simulated_dataset = "temp_brain.zinbwave.ks_statistic_sc_features/_viash_par/input_simulated_dataset_1/brain.zinbwave.generate_sim_spatialcluster.output_sp.h5ad",
  output = "output.h5ad"
)
meta <- list(
  name = "ks_statistic",
  resources_dir = "target/executable/metrics/ks_statistic_sc_features"
)
## VIASH END

source(paste0(meta[["resources_dir"]], "/utils.R"))

input_real_sp <- anndataR::read_h5ad(par$input_spatial_dataset)
input_sc <- anndataR::read_h5ad(par$input_singlecell_dataset)
input_simulated_sp <- anndataR::read_h5ad(par$input_simulated_dataset)

organism_mapping <- c(
  mus_musculus = "Mus musculus",
  homo_sapiens = "Homo sapiens"
)
sc_species <- organism_mapping[[input_real_sp$uns[["dataset_organism"]]]]
if (is.null(sc_species)) {
  sc_species <- input_real_sp$uns[["dataset_organism"]]
}

real_log_count <- t(input_real_sp$layers[["logcounts"]])
real_prob_matrix <- input_real_sp$obsm[["celltype_proportions"]]
colnames(real_prob_matrix) <- paste0("ct", seq_len(ncol(real_prob_matrix)))

sim_sce <- scater::logNormCounts(SingleCellExperiment::SingleCellExperiment(
  list(counts = Matrix::t(input_simulated_sp$layers[["counts"]])),
  colData = input_simulated_sp$obs,
  metadata = input_simulated_sp$obsm
))

sim_log_count <- SummarizedExperiment::assay(sim_sce, "logcounts")

ensure_spot_names <- function(x, spot_names, label) {
  if (is.null(rownames(x))) {
    if (nrow(x) != length(spot_names)) {
      stop(label, " has no row names and cannot be matched to spots.")
    }
    rownames(x) <- spot_names
  }
  x
}

align_spatial_inputs <- function(data, coords, spot_probability, label) {
  coords <- as.data.frame(coords)
  spot_probability <- as.matrix(spot_probability)

  if (is.null(colnames(data))) {
    stop(label, " expression matrix has no column names.")
  }

  coords <- ensure_spot_names(coords, colnames(data), paste(label, "coordinates"))
  spot_probability <- ensure_spot_names(
    spot_probability,
    colnames(data),
    paste(label, "spot probabilities")
  )

  common_spots <- colnames(data)[
    colnames(data) %in% rownames(coords) &
      colnames(data) %in% rownames(spot_probability)
  ]

  if (length(common_spots) < 2) {
    stop(label, " has fewer than two spots shared by expression, coordinates, and spot probabilities.")
  }

  data <- data[, common_spots, drop = FALSE]
  coords <- coords[common_spots, c("row", "col"), drop = FALSE]
  spot_probability <- spot_probability[common_spots, , drop = FALSE]

  keep_spots <- is.finite(as.numeric(coords$row)) &
    is.finite(as.numeric(coords$col)) &
    apply(is.finite(spot_probability), 1, all)

  if (sum(!keep_spots) > 0) {
    warning(label, ": dropping ", sum(!keep_spots), " spots with non-finite coordinates or probabilities.")
  }

  data <- data[, keep_spots, drop = FALSE]
  coords <- coords[keep_spots, , drop = FALSE]
  spot_probability <- spot_probability[keep_spots, , drop = FALSE]

  if (ncol(data) < 2) {
    stop(label, " has fewer than two valid spots after alignment.")
  }

  list(
    data = data,
    sample = rep("sample1", ncol(data)),
    spatialCoords = coords,
    spotProbability = t(spot_probability)
  )
}

as_finite_kde_input <- function(x) {
  if (is.data.frame(x)) {
    x <- as.matrix(x)
  }

  if (is.matrix(x)) {
    x <- matrix(as.numeric(x), nrow = nrow(x), ncol = ncol(x))
    return(x[apply(is.finite(x), 1, all), , drop = FALSE])
  }

  x <- as.numeric(x)
  x[is.finite(x)]
}

add_kde_jitter <- function(x, amount) {
  x_range <- max(x, na.rm = TRUE) - min(x, na.rm = TRUE)
  x_scale <- max(x_range, stats::sd(as.numeric(x), na.rm = TRUE), 1)
  noise <- seq(-amount, amount, length.out = length(x)) * x_scale

  if (is.matrix(x)) {
    matrix(as.numeric(x) + noise, nrow = nrow(x), ncol = ncol(x))
  } else {
    as.numeric(x) + noise
  }
}

try_kde_test <- function(x1, x2) {
  x1 <- as_finite_kde_input(x1)
  x2 <- as_finite_kde_input(x2)

  if (length(x1) == 0 || length(x2) == 0) {
    warning("No finite values available for ks::kde.test; returning NA.")
    return(list(zstat = NA_real_, tstat = NA_real_))
  }

  if (is.matrix(x1) && (nrow(x1) < 2 || nrow(x2) < 2)) {
    warning("Not enough finite rows available for ks::kde.test; returning NA.")
    return(list(zstat = NA_real_, tstat = NA_real_))
  }

  if (!is.matrix(x1) && (length(x1) < 2 || length(x2) < 2)) {
    warning("Not enough finite values available for ks::kde.test; returning NA.")
    return(list(zstat = NA_real_, tstat = NA_real_))
  }

  last_error <- NULL
  for (jitter in c(0, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4)) {
    x1_try <- if (jitter == 0) x1 else add_kde_jitter(x1, jitter)
    x2_try <- if (jitter == 0) x2 else add_kde_jitter(x2, jitter)

    result <- tryCatch(
      ks::kde.test(x1 = x1_try, x2 = x2_try),
      error = function(e) {
        last_error <<- e
        NULL
      }
    )

    if (!is.null(result)) {
      if (jitter > 0) {
        warning(
          "Caught error in ks::kde.test: ",
          last_error$message,
          "\n\nSucceeded after adding deterministic jitter of size ",
          jitter,
          "."
        )
      }
      return(result)
    }
  }

  warning(
    "ks::kde.test failed after deterministic jitter retries: ",
    last_error$message,
    "\n\nReturning NA for this metric."
  )
  list(zstat = NA_real_, tstat = NA_real_)
}

safe_scfeatures <- function(aligned_inputs, label, feat_types) {
  tryCatch(
    scFeatures::scFeatures(
      aligned_inputs$data,
      sample = aligned_inputs$sample,
      spatialCoords = aligned_inputs$spatialCoords,
      feature_types = feat_types,
      type = "spatial_t",
      species = sc_species,
      spotProbability = aligned_inputs$spotProbability
    ),
    error = function(e) {
      warning(label, " scFeatures failed: ", e$message, "\n\nReturning NA for this metric component.")
      list()
    }
  )
}

safe_feature_values <- function(scfeatures_result, feature_name, label) {
  feature <- scfeatures_result[[feature_name]]
  if (is.null(feature)) {
    warning(label, " feature '", feature_name, "' is missing; returning NA.")
    return(NA_real_)
  }

  values <- suppressWarnings(as.numeric(unlist(feature, use.names = FALSE)))
  values <- values[is.finite(values)]
  if (length(values) == 0) {
    warning(label, " feature '", feature_name, "' has no finite values; returning NA.")
    return(NA_real_)
  }

  values
}

cat("Computing cell type proportions in simulated data\n")
sim_prob_matrix <- CARD_processing(input_simulated_sp, input_sc)

feat_types <- c("L_stats", "nn_correlation", "morans_I")

real_inputs <- align_spatial_inputs(
  data = real_log_count,
  coords = input_real_sp$obs[c("row", "col")],
  spot_probability = real_prob_matrix,
  label = "real"
)

sim_inputs <- align_spatial_inputs(
  data = sim_log_count,
  coords = input_simulated_sp$obs[c("row", "col")],
  spot_probability = sim_prob_matrix,
  label = "simulated"
)

real_scfeatures_result <- safe_scfeatures(real_inputs, "real", feat_types)
sim_scfeatures_result <- safe_scfeatures(sim_inputs, "simulated", feat_types)

ks_statistic_L_stats <- try_kde_test(
  x1 = safe_feature_values(real_scfeatures_result, "L_stats", "real"),
  x2 = safe_feature_values(sim_scfeatures_result, "L_stats", "simulated")
)

ks_statistic_nn_correlation <- try_kde_test(
  x1 = safe_feature_values(real_scfeatures_result, "nn_correlation", "real"),
  x2 = safe_feature_values(sim_scfeatures_result, "nn_correlation", "simulated")
)
ks_statistic_morans_I <- try_kde_test(
  x1 = safe_feature_values(real_scfeatures_result, "morans_I", "real"),
  x2 = safe_feature_values(sim_scfeatures_result, "morans_I", "simulated")
)

cat("Combining metric values\n")
uns_metric_ids <- c(
  "ks_statistic_L_stats",
  
  "ks_statistic_nn_correlation",
  "ks_statistic_morans_I"
)

uns_metric_values <- c(
  ks_statistic_L_stats$zstat,
  
  ks_statistic_nn_correlation$zstat,
  ks_statistic_morans_I$zstat
)

cat("Writing output AnnData to file\n")
output <- anndataR::AnnData(
  uns = list(
    dataset_id = input_simulated_sp$uns[["dataset_id"]],
    method_id = input_simulated_sp$uns[["method_id"]],
    metric_ids = uns_metric_ids,
    metric_values = uns_metric_values
  ),
  shape = c(0L, 0L)
)
output$write_h5ad(par[["output"]], compression = "gzip")
