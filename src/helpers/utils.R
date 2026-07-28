# log-normalisation
#
# Both the real and the simulated dataset must be normalised the exact same way
# before they are compared, otherwise the metric picks up the difference in
# normalisation rather than the difference in the data. Do not read the
# `logcounts` layer directly: on the real dataset it was produced by the dataset
# loader, on a simulated dataset it may not exist at all.
compute_logcounts <- function(adata) {
  requireNamespace("scater", quietly = TRUE)
  requireNamespace("SingleCellExperiment", quietly = TRUE)
  requireNamespace("SummarizedExperiment", quietly = TRUE)

  # genes x spots, as expected by SingleCellExperiment
  counts <- Matrix::t(adata$layers[["counts"]])

  sce <- SingleCellExperiment::SingleCellExperiment(list(counts = counts))
  sce <- scater::logNormCounts(sce)

  logcounts <- SummarizedExperiment::assay(sce, "logcounts")
  dimnames(logcounts) <- dimnames(counts)

  # back to spots x genes, matching the anndata layout
  Matrix::t(logcounts)
}

# spatial autocorrelation
generate_moransI <- function(adata) {
  requireNamespace("spots", quietly = TRUE)

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


generate_cosine <- function(real, sim) {
  requireNamespace("lsa", quietly = TRUE)

  real_i <- as.vector(real$Morans.I)
  sim_i <- as.vector(sim$Morans.I)

  # drop the pairs where either side is missing. This used to filter `real` and
  # `sim` themselves, which are lists, so it never dropped anything and a single
  # NaN took the whole metric with it -- 32 of 99 runs came out NA, while
  # generate_mantel() passes na.rm and lost none.
  keep <- is.finite(real_i) & is.finite(sim_i)
  if (sum(keep) < 2) {
    return(NA_real_)
  }

  similarity <- lsa::cosine(lsa::as.textmatrix(cbind(real_i[keep], sim_i[keep])))
  return(mean(similarity))
}

generate_mantel <- function(real, sim) {
  requireNamespace("vegan", quietly = TRUE)
  mantel_test <- vegan::mantel(real$Morans.I, sim$Morans.I, na.rm = TRUE, method = "pearson")
  return(mantel_test$statistic)
}

# Spatial varianble gene
generate_svg_sparkx <- function(adata) {
  requireNamespace("SPARK", quietly = TRUE)

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

calculate_precision <- function(real_svg, sim_svg) {
  filtered_real_data <- dplyr::filter(real_svg$res_mtest, adjustedPval < 0.05)
  filtered_compared_data <- dplyr::filter(sim_svg$res_mtest, adjustedPval < 0.05)
  tp <- length(intersect(row.names(filtered_real_data), row.names(filtered_compared_data)))
  fp <- length(setdiff(row.names(filtered_compared_data), row.names(filtered_real_data)))
  # a simulation that produces no spatially variable genes at all is a bad
  # simulation, not a missing result. Reporting NA meant both negative controls
  # had no score on any dataset, so nothing anchored the bottom of the scale.
  precision <- ifelse((tp + fp) > 0, tp / (tp + fp), 0)
  return(precision)
}

calculate_recall <- function(real_svg, sim_svg) {
  filtered_real_data <- dplyr::filter(real_svg$res_mtest, adjustedPval < 0.05)
  filtered_compared_data <- dplyr::filter(sim_svg$res_mtest, adjustedPval < 0.05)
  tp <- length(intersect(row.names(filtered_real_data), row.names(filtered_compared_data)))
  fn <- length(setdiff(row.names(filtered_real_data), row.names(filtered_compared_data)))
  # unlike precision, an empty denominator here says something about the real
  # dataset rather than about the simulation: there are no spatially variable
  # genes to recover, so no simulator can be scored on it. NA is right.
  recall <- ifelse((tp + fn) > 0, tp / (tp + fn), NA)
  return(recall)
}

# cell type deconvolution
CARD_processing <- function(sp_adata, sc_adata) {
  requireNamespace("MuSiC", quietly = TRUE)
  requireNamespace("CARD", quietly = TRUE)

  spatial_count <- Matrix::t(sp_adata$layers[["counts"]])
  spatial_location <- data.frame(
    x = as.numeric(sp_adata$obs$col),
    y = as.numeric(sp_adata$obs$row)
  )
  sc_count <- Matrix::t(sc_adata$layers[["counts"]])
  sc_meta <- data.frame(
    cellID = rownames(sc_adata),
    cellType = sc_adata$obs$cell_type,
    sampleInfo = sc_adata$obs$donor_id
  )
  rownames(sc_meta) <- sc_meta$cellID
  rownames(spatial_location) <- colnames(spatial_count)

  CARD_obj <- CARD::createCARDObject(
	  sc_count = sc_count,
	  sc_meta = sc_meta,
	  spatial_count = spatial_count,
	  spatial_location = spatial_location,
	  ct.varname = "cellType",
	  ct.select = unique(sc_meta$cellType),
	  sample.varname = "sampleInfo",
	  minCountGene = 100,
	  minCountSpot = 5)
  CARD_obj <- CARD::CARD_deconvolution(CARD_object = CARD_obj)

  Proportion_CARD <- as.matrix(CARD_obj@Proportion_CARD)

  return(Proportion_CARD)

}

generate_jds <- function(real, sim) {
  common_row_names <- intersect(rownames(real), rownames(sim))
  real_common <- real[common_row_names, , drop = FALSE]
  sim_common <- sim[common_row_names, , drop = FALSE]
  jsd_values <- sapply(1:nrow(real_common), function(i) {
    x.count <- rbind(as.vector(real_common[i, ]), as.vector(sim_common[i, ]))
    philentropy::JSD(x.count, est.prob = "empirical")
  })
  average_jsd <- mean(jsd_values)
  return(average_jsd)
}


generate_rmse <- function(real, sim) {
  common_row_names <- intersect(rownames(real), rownames(sim))
  real_common <- real[common_row_names, , drop = FALSE]
  sim_common <- sim[common_row_names, , drop = FALSE]
  rmse_values <- sapply(1:nrow(real_common), function(i) {
    Metrics::rmse(as.vector(real_common[i, ]), as.vector(sim_common[i, ]))
  })
  average_rmse <- mean(rmse_values)
  return(average_rmse)
}


# spatial clustering
reclassify_simsce <- function(location, real_cluster, sim_cluster) {
  flatten_cluster <- function(x, label) {
    if (is.data.frame(x) || is.matrix(x)) {
      if (ncol(x) != 1) {
        stop(label, " must contain exactly one column.")
      }
      x <- x[, 1]
    }
    unname(as.vector(x))
  }

  real_cluster <- flatten_cluster(real_cluster, "Real cluster")
  sim_cluster <- flatten_cluster(sim_cluster, "Simulated cluster")
  location <- as.vector(location)

  if (
    length(location) != length(real_cluster) ||
      length(location) != length(sim_cluster)
  ) {
    stop("Location, real cluster, and simulated cluster vectors must have equal lengths.")
  }

  keep <- !is.na(location) & !is.na(real_cluster) & !is.na(sim_cluster)
  overlap <- table(
    real = as.character(real_cluster[keep]),
    sim = as.character(sim_cluster[keep])
  )

  cluster_map <- stats::setNames(
    rep(NA_character_, ncol(overlap)),
    colnames(overlap)
  )

  for (i in seq_len(min(nrow(overlap), ncol(overlap)))) {
    max_value <- max(overlap)
    if (!is.finite(max_value) || max_value == 0) {
      break
    }

    indices <- which(overlap == max_value, arr.ind = TRUE)[1, ]
    real_index <- indices[1]
    sim_index <- indices[2]
    cluster_map[colnames(overlap)[sim_index]] <- rownames(overlap)[real_index]
    overlap[real_index, ] <- -Inf
    overlap[, sim_index] <- -Inf
  }

  reclassified <- unname(cluster_map[as.character(sim_cluster)])
  if (is.numeric(real_cluster)) {
    reclassified <- as.numeric(reclassified)
  }

  reclassified
}


# generate sparial clustering in simulated data
generate_sim_spatialCluster <- function(real_adata, sim_adata){
  # colnames(sim_adata$obs)[colnames(sim_adata$obs) == "col"] <- "array_col"
  # colnames(sim_adata$obs)[colnames(sim_adata$obs) == "row"] <- "array_row"
  sim_adata$obs["array_col"] <- sim_adata$obs["col"]
  sim_adata$obs["array_row"] <- sim_adata$obs["row"]
  sim_sce <- scater::logNormCounts(SingleCellExperiment::SingleCellExperiment(
    list(counts = Matrix::t(sim_adata$layers[["counts"]])),
    colData = sim_adata$obs,
    metadata = sim_adata$obsm))
  sim_sce <- BayesSpace::spatialPreprocess(sim_sce, n.PCs=7, platform="ST", n.HVGs=2000, log.normalize=FALSE)
  real_clusters <- unique(real_adata$obs[["spatial_cluster"]])
  real_clusters <- real_clusters[!is.na(real_clusters)]
  if (length(real_clusters) < 2) {
    stop("At least two real spatial clusters are required.")
  }

  sim_sce <- BayesSpace::spatialCluster(sim_sce,
    q = length(real_clusters),
    platform = "ST",
    d = 7,
    init.method = "mclust", model = "t", gamma = 2,
    nrep = 1000, burn.in = 100,
    save.chain = TRUE
  )
  sim_cluster <- sim_sce$spatial.cluster
  return(sim_cluster)
}
