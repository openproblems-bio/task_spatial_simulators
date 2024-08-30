suppressMessages(library(SingleCellExperiment, quietly = TRUE))
suppressMessages(library(SymSim, quietly = TRUE))

## VIASH START
par <- list(
  input = "resources_test/spatialsimbench_mobnew/MOBNEW.rds",
  base = "domain"
)
meta <- list(
  name = "symsim"
)
## VIASH END

cat("Reading input files\n")
input <- anndata::read_h5ad(par$input)

sce <- SingleCellExperiment(
  list(counts = Matrix::t(input$layers[["counts"]])),
  colData = input$obs
)

cat("SymSim simulation start\n")

if (par$base != "domain") {
  stop("ONLY domain base")
}

simulated_result <- NULL
tech <-  "UMI"

ordered_indices <- order(colData(sce)$spatial_cluster)
sce_ordered <- sce[, ordered_indices]

for (thisSpatialCluster in (unique(sce_ordered$spatial_cluster))) {
  res <- try({
    # subset to one cell type
    sce_thiscelltype <- sce_ordered[, sce_ordered$spatial_cluster == thisSpatialCluster]

    # this is because if some genes are 0 , this will cause error in simulation
    keep_feature <- rowSums(counts(sce_thiscelltype) > 0) > 0
    sce_thiscelltype_f <- sce_thiscelltype[keep_feature, ]

    best_matches_UMI <- BestMatchParams(
      tech = "UMI",
      counts = as.matrix(counts(sce_thiscelltype_f)),
      plotfilename = "best_params.umi.qqplot",
      n_optimal = 1
    )

    sim_thiscelltype <- SimulateTrueCounts(
      ncells_total =  dim(sce_thiscelltype)[2],
      ngenes =  dim(sce_thiscelltype)[1],
      evf_type = "one.population",
      randseed = 1,
      Sigma = best_matches_UMI$Sigma[1],
      gene_effects_sd = best_matches_UMI$gene_effects_sd[1],
      scale_s = best_matches_UMI$scale_s[1],
      gene_effect_prob = best_matches_UMI$gene_effect_prob[1],
      prop_hge = best_matches_UMI$prop_hge[1],
      mean_hge = best_matches_UMI$mean_hge[1]
    )

    gene_len <- sample(gene_len_pool, dim(sce_thiscelltype)[1], replace = FALSE)
    sim_thiscelltype <- True2ObservedCounts(
      true_counts = sim_thiscelltype[[1]],
      meta_cell = sim_thiscelltype[[3]],
      protocol = tech,
      alpha_mean = best_matches_UMI$alpha_mean[1],
      alpha_sd = best_matches_UMI$alpha_sd[1],
      gene_len = gene_len,
      depth_mean = best_matches_UMI$depth_mean[1],
      depth_sd = best_matches_UMI$depth_sd[1]
    )

    # tidy up the names
    sim_thiscelltype <- SingleCellExperiment(list(counts = sim_thiscelltype$counts))
    sim_thiscelltype$spatial_cluster <- thisSpatialCluster

    # combine the cell types
    if (is.null(simulated_result)) {
      simulated_result <- sim_thiscelltype
    } else {
      simulated_result <- SingleCellExperiment::cbind(simulated_result, sim_thiscelltype)
    }

  })
}

colnames(simulated_result) <- colnames(sce_ordered)
rownames(simulated_result) <- rownames(sce_ordered)

simulated_result_order <- sce_ordered
counts(simulated_result_order) <- counts(simulated_result)

simulated_result_order <- simulated_result_order[, match(colnames(sce), colnames(simulated_result_order))]
simulated_result_order <- simulated_result_order[match(rownames(sce), rownames(simulated_result_order)), ]
new_obs <- as.data.frame(simulated_result_order@colData[c("row", "col")])

cat("Generating output\n")

output <- anndata::AnnData(
  layers = list(
    counts = Matrix::t(counts(simulated_result_order))
  ),
  obs = new_obs,
  var = input$var,
  uns = c(
    input$uns,
    list(
      method_id = meta$name
    )
  )
)

cat("Write output files\n")
output$write_h5ad(par$output, compression = "gzip")
