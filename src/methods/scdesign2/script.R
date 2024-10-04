# requireNamespace("scDesign2", quietly = TRUE)

## VIASH START
par <- list(
  input = "resources/task_spatial_simulators/datasets/breast/output_sp.h5ad",
  output = "simulated_dataset.h5ad",
  base = "domain"
)
meta <- list(
  name = "scdesign2",
  cpus = 8L
)
## VIASH END

cat("Reading input files\n")
input <- anndata::read_h5ad(par$input)

if (par$base != "domain") {
  stop("ONLY domain base")
}

cat("scDesign2 simulation start\n")
counts <- as.matrix(Matrix::t(input$layers[["counts"]]))
colnames(counts) <- as.character(input$obs$spatial_cluster)
spatial_cluster_prop <- table(input$obs$spatial_cluster)

copula_result <- scDesign2::fit_model_scDesign2(
  data_mat = counts,
  cell_type_sel = unique(colnames(counts)),
  sim_method = "copula",
  ncores = ifelse(!is.null(meta$cpus), meta$cpus, 8L)
)

sim_out <- scDesign2::simulate_count_scDesign2(
  copula_result,
  ncol(counts),
  sim_method = "copula",
  cell_type_prop = prop.table(spatial_cluster_prop)
)

# colnames(sim_out) <- colnames(t(as.matrix(input$layers[["counts"]])))
#rownames(sim_out) <- rownames(t(as.matrix(input$layers[["counts"]])))

rownames(sim_out) <- input$obs_names
colnames(sim_out) <- input$var_names

cat("Generating output\n")

output <- anndata::AnnData(
  layers = list(
    counts = Matrix::t(sim_out)
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
