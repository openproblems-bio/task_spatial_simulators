#!/bin/bash

# download and process all datasets from figshare
# https://figshare.com/articles/dataset/SpatialSimBench_dataset/26054188

cat > /tmp/params.yaml <<'EOF'
param_list:
  - id: breast
    input_sc: resources/datasets_raw/BREAST_sc.rds
    input_sp: resources/datasets_raw/BREAST.rds
    dataset_id: breast
    dataset_name: Breast
    dataset_url: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176078
    dataset_reference: 10.1038/s41588-021-00911-1 
    dataset_summary: A single-cell and spatially resolved atlas of human breast cancers
    dataset_description: "This study presents a single cell and spatially resolved transcriptomics analysis of human breast cancers. We develop a single cell method of intrinsic subtype classification (scSubtype) to reveal recurrent neoplastic cell heterogeneity. Immunophenotyping using CITE-Seq provides high-resolution immune profiles, including novel PD-L1/PD-L2+ macrophage populations associated with clinical outcome. Mesenchymal cells displayed diverse functions and cell surface protein expression through differentiation within 3 major lineages. Stromal-immune niches were spatially organized in tumors, offering insights into anti-tumor immune regulation. Using single cell signatures, we deconvoluted large breast cancer cohorts to stratify them into nine clusters, termed ‘ecotypes’, with unique cellular compositions and clinical outcomes. This study provides a comprehensive transcriptional atlas of the cellular architecture of breast cancer."
    dataset_organism: "homo_sapiens"
    dataset_assay_spatial: Visium # EFO:0010961
    dataset_assay_singlecell: Chromium
  
  - id: osteosarcoma
    input_sc: resources/datasets_raw/HOSTEOSARCOMA_sc.rds
    input_sp: resources/datasets_raw/HOSTEOSARCOMA.rds
    dataset_id: osteosarcoma
    dataset_name: Osteosarcoma
    dataset_url: https://www.pnas.org/doi/full/10.1073/pnas.1912459116
    dataset_reference: 10.1073/pnas.1912459116
    dataset_summary: Spatial profiling of human osteosarcoma cells.
    dataset_description: "Spatial transcriptome profiling by MERFISH reveals subcellular RNA compartmentalization and cell cycle-dependent gene expression."
    dataset_organism: "homo_sapiens"
    # note: Where is the single-cell data coming from? 
    #   The DOI above only seems to contain MERFISH data
    #   -> https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152048 10X chromium BC22

publish_dir: "resources/datasets"
output_sc: "$id/output_sc.h5ad"
output_sp: "$id/output_sp.h5ad"
output_state: "$id/state.yaml"
EOF

nextflow run . \
  -main-script target/nextflow/process_datasets/convert/main.nf \
  -profile docker \
  -params-file /tmp/params.yaml
