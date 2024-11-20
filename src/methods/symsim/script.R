suppressMessages(library(SingleCellExperiment, quietly = TRUE))
suppressMessages(library(SymSim, quietly = TRUE))

## VIASH START
par <- list(
  input = "resources_test/spatialsimbench_mobnew/dataset_sp.h5ad",
  base = "domain"
)
meta <- list(
  name = "symsim"
)
## VIASH END

cat("Reading input files\n")
input <- anndata::read_h5ad(par$input)

cat("SymSim simulation start\n")

if (par$base != "domain") {
  stop("ONLY domain base")
}

simulated_result <- NULL
tech <-  "UMI"

ordered_indices <- order(input$obs$spatial_cluster)
input_ordered <- input[ordered_indices]

for (thisSpatialCluster in unique(input_ordered$obs[["spatial_cluster"]])) {
  res <- try({
    input_thiscelltype <- input_ordered[input_ordered$obs[["spatial_cluster"]] == thisSpatialCluster]

    # this is because if some genes are 0 , this will cause error in simulation
    keep_feature <- colSums(input_thiscelltype$layers[["counts"]] > 0) > 0
    input_thiscelltype_f <- input_thiscelltype[, keep_feature]

    best_matches_UMI <- BestMatchParams(
      tech = "UMI",
      counts = as.matrix(t(input_thiscelltype_f$layers[["counts"]])),
      plotfilename = "best_params.umi.qqplot",
      n_optimal = 1
    )

    sim_thiscelltype <- SimulateTrueCounts(
      ncells_total =  dim(input_thiscelltype)[1],
      ngenes =  dim(input_thiscelltype)[2],
      evf_type = "one.population",
      randseed = 1,
      Sigma = best_matches_UMI$Sigma[1],
      gene_effects_sd = best_matches_UMI$gene_effects_sd[1],
      scale_s = best_matches_UMI$scale_s[1],
      gene_effect_prob = best_matches_UMI$gene_effect_prob[1],
      prop_hge = best_matches_UMI$prop_hge[1],
      mean_hge = best_matches_UMI$mean_hge[1]
    )

    gene_len <- sample(gene_len_pool, dim(input_thiscelltype)[2], replace = FALSE)
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
    sim_thiscelltype <- SingleCellExperiment(list(counts = as.matrix(sim_thiscelltype$counts)))
    sim_thiscelltype$spatial_cluster <- thisSpatialCluster

    # combine the cell types
    if (is.null(simulated_result)) {
      simulated_result <- sim_thiscelltype
    } else {
      simulated_result <- SingleCellExperiment::cbind(simulated_result, sim_thiscelltype)
    }

  })
}

colnames(simulated_result) <- rownames(input_ordered$obs)
rownames(simulated_result) <- rownames(input_ordered$var)

simulated_result_ordered <- counts(simulated_result)[
  match(rownames(counts(simulated_result)), rownames(input_ordered$var)),
  match(colnames(counts(simulated_result)), rownames(input_ordered$obs))
]

cat("Generating output\n")
output <- anndata::AnnData(
  layers = list(
    counts = Matrix::t(simulated_result_ordered)
  ),
  obs = input_ordered$obs[c("row", "col")],
  var = input_ordered$var,
  uns = c(
    input$uns,
    list(
      method_id = meta$name
    )
  )
)

cat("Write output files\n")
output$write_h5ad(par$output, compression = "gzip")
