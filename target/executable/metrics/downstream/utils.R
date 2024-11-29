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
  real_new <- real[!is.na(real) & !is.na(sim)]
  sim_new <- sim[!is.na(real) & !is.na(sim)]
  similarity <- lsa::cosine(lsa::as.textmatrix(cbind(as.vector(real_new$Morans.I), as.vector(sim_new$Morans.I))))
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
  precision <- ifelse((tp + fp) > 0, tp / (tp + fp), NA)
  return(precision)
}

calculate_recall <- function(real_svg, sim_svg) {
  filtered_real_data <- dplyr::filter(real_svg$res_mtest, adjustedPval < 0.05)
  filtered_compared_data <- dplyr::filter(sim_svg$res_mtest, adjustedPval < 0.05)
  tp <- length(intersect(row.names(filtered_real_data), row.names(filtered_compared_data)))
  fn <- length(setdiff(row.names(filtered_real_data), row.names(filtered_compared_data)))
  recall <- ifelse((tp + fn) > 0, tp / (tp + fn), NA)
  return(recall)
}

# cell type deconvolution
CARD_processing <- function(sp_adata, sc_adata){
  requireNamespace("MuSiC", quietly = TRUE)
  requireNamespace("CARD", quietly = TRUE)
  spatial_count <- Matrix::t(sp_adata$layers[["counts"]])
  spatial_location <- cbind.data.frame(
    x = sp_adata$obs$col,
    y = sp_adata$obs$row
  )
  spatial_location$col <- as.numeric(spatial_location$col)
  spatial_location$row <- as.numeric(spatial_location$row)
  sc_count <- Matrix::t(sc_adata$layers[["counts"]])
  sc_meta <- cbind.data.frame(
    cellID = rownames(sc_adata),
    cellType = input_sc$obs$cell_type,
    sampleInfo = input_sc$obs$donor_id
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
reclassify_simsce <- function(location, real_cluster, sim_cluster){
  test <- data.frame(loc = location, real = real_cluster, sim = sim_cluster)
  matrix_counts <- matrix(0, nrow = 4, ncol = 4,
                          dimnames = list(paste0("real", 1:4), paste0("sim", 1:4)))
  for (r in 1:4) {
    for (s in 1:4) {
      subset_data <- subset(test, real == r & sim == s)
      matrix_counts[r, s] <- length(unique(subset_data$loc))
    }
  }
  reclassification <- numeric(4)
  for (i in 1:4) {
    max_value <- max(matrix_counts)
    if (max_value == 0) break
    indices <- which(matrix_counts == max_value, arr.ind = TRUE)
    real_index <- indices[1, 1]
    sim_index <- indices[1, 2]
    reclassification[sim_index] <- real_index
    matrix_counts[real_index, ] <- -Inf
    matrix_counts[, sim_index] <- -Inf
  }
  cluster_map <- setNames(reclassification, seq_along(reclassification))
  sim_reclassify_cluster <- cluster_map[sim_cluster]
  return(sim_reclassify_cluster)
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
  sim_sce <- BayesSpace::spatialCluster(sim_sce,
    q = max(unique(real_adata$obs[, c("spatial_cluster")])),
    platform = "ST",
    d = 7,
    init.method = "mclust", model = "t", gamma = 2,
    nrep = 1000, burn.in = 100,
    save.chain = TRUE
  )
  sim_cluster <- sim_sce$spatial.cluster
  return(sim_cluster)
}
