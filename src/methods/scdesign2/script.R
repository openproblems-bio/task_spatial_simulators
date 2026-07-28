## VIASH START
par <- list(
  input = "resources_test/spatialsimbench_mobnew/dataset_sp.h5ad",
  output = "simulated_dataset.h5ad",
  base = "domain"
)
meta <- list(
  name = "scdesign2",
  cpus = 8L
)
## VIASH END

cat("Reading input files\n")
input <- anndataR::read_h5ad(par$input)

if (par$base != "domain") {
  stop("ONLY domain base")
}

cat("scDesign2 simulation start\n")
# simulate_count_scDesign2() returns the cells grouped per cell type, so order
# the spots by spatial cluster first -- otherwise the coordinates we attach
# below belong to different spots than the counts we attach them to
ordered_indices <- order(input$obs$spatial_cluster)
input_ordered <- input[ordered_indices]

counts <- Matrix::t(input_ordered$layers[["counts"]])
colnames(counts) <- as.character(input_ordered$obs$spatial_cluster)

# cell_type_prop is matched to the fitted models by position, so both have to be
# in the same order -- and that order has to be the one the spots are sorted in,
# which unique() gives us since the columns are already grouped per cluster
cell_types <- unique(colnames(counts))
spatial_cluster_prop <- prop.table(table(colnames(counts))[cell_types])

copula_result <- scDesign2::fit_model_scDesign2(
  data_mat = counts,
  cell_type_sel = cell_types,
  sim_method = "copula",
  ncores = ifelse(!is.null(meta$cpus), meta$cpus, 8L)
)

sim_out <- scDesign2::simulate_count_scDesign2(
  copula_result,
  ncol(counts),
  sim_method = "copula",
  cell_type_prop = spatial_cluster_prop
)

rownames(sim_out) <- input_ordered$var_names
colnames(sim_out) <- input_ordered$obs_names

# put the spots back in the order they came in. The metrics compare the
# simulated dataset against the real one spot for spot -- ARI and NMI are
# invariant to a permutation of the labels, but not of the samples.
restore_order <- order(ordered_indices)

cat("Generating output\n")

output <- anndataR::AnnData(
  layers = list(
    counts = Matrix::t(sim_out)[restore_order, , drop = FALSE]
  ),
  obs = input$obs[c("row", "col")],
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
