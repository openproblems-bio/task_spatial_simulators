# spatial autocorrelation
generate_moransI <- function(adata){
  sp_count <- adata$layers[["logcounts"]]
  
  W <- cbind.data.frame(row = as.numeric(sapply(strsplit(rownames(sp_count),split="x"),"[",1)),
                        col = as.numeric(sapply(strsplit(rownames(sp_count),split="x"),"[",2)))
  W <- 1/as.matrix(dist(W))
  diag(W) <- 0
  res <- spots::BivariateMoransI( X =  sp_count , W =  W)
  return(res$Morans.I)
}

# Spatial varianble gene
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

# cell type deconvolution
CARD_processing <- function(sp_adata, sc_adata){
  sp_adata <- input_real_sp
  sc_adata <- input_sc
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
  
  metadata(sp_obj)$celltype_prop <- CARD_obj@Proportion_CARD
  
  diff_list <- setdiff(colnames(counts(sp_obj)), rownames(sp_obj@metadata$celltype_prop))
  current_colnames <- colnames(sp_obj)
  columns_to_keep <- !current_colnames %in% diff_list
  sp_obj <- sp_obj[, columns_to_keep]

  return(sp_obj)
}
# spatial clustering