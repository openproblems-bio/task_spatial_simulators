# library(dplyr) need this package
requireNamespace("dplyr", quietly = TRUE)
requireNamespace("spots", quietly = TRUE)

# spatial autocorrelation
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

generate_cosine <- function(real, sim){
  real_new <- real[!is.na(real) & !is.na(sim)]
  sim_new <- sim[!is.na(real) & !is.na(sim)]
  similarity <- lsa::cosine(lsa::as.textmatrix(cbind(as.vector(real_new), as.vector(sim_new))))
  return(mean(similarity))
}

generate_mantel <- function(real, sim){
  mantel_test <- vegan::mantel( real , sim, na.rm = T, method="pearson")
  return(mantel_test$statistic)
}

mantel.test <- vegan::mantel( real_moransI , sim_moransI, na.rm = T, method="pearson")

# Spatial varianble gene
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

calculate_precision <- function(real_svg, sim_svg) {
  filtered_real_data <- real_svg$res_mtest %>% 
    dplyr::filter(adjustedPval < 0.05)
  filtered_compared_data <- sim_svg$res_mtest %>% 
    dplyr::filter(adjustedPval < 0.05)
  tp <- length(intersect(row.names(filtered_real_data), row.names(filtered_compared_data)))
  fp <- length(setdiff(row.names(filtered_compared_data), row.names(filtered_real_data)))
  precision <- ifelse((tp + fp) > 0, tp / (tp + fp), NA)
  return(precision)
}

calculate_recall <- function(real_svg, sim_svg) {
  filtered_real_data <- real_svg$res_mtest %>% 
    filter(adjustedPval < 0.05)
  filtered_compared_data <- sim_svg$res_mtest %>% 
    filter(adjustedPval < 0.05)
  tp <- length(intersect(row.names(filtered_real_data), row.names(filtered_compared_data)))
  fn <- length(setdiff(row.names(filtered_real_data), row.names(filtered_compared_data)))
  recall <- ifelse((tp + fn) > 0, tp / (tp + fn), NA)
  return(recall)
}

# cell type deconvolution
CARD_processing <- function(sp_adata, sc_adata){
  # sp_adata <- input_real_sp
  # sc_adata <- input_sc
  spatial_count <- Matrix::t(sp_adata$layers[["counts"]])
  spatial_location <- cbind.data.frame(
    x = as.numeric(sapply(strsplit(colnames(spatial_count),split="x"),"[",1)),
    y = as.numeric(sapply(strsplit(colnames(spatial_count),split="x"),"[",2)))
  # names(spatial_location) <- c("x", "y")
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

