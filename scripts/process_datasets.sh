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

  - id: prostate
    input_sc: resources/datasets_raw/HPROSTATE_sc.rds
    input_sp: resources/datasets_raw/HPROSTATE.rds
    dataset_id: prostate
    dataset_name: Prostate
    dataset_url: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE159697
    dataset_reference: 10.1016/j.isci.2021.102640 
    dataset_summary: Spatially resolved gene expression of human protate tissue slices treated with steroid hormones for 8 hours
    dataset_description: "Spatially resolved gene expression was prepard by dissociated hman prostate tissue to single cells, and collected & prepped for RNA-seq using the Visium Spatial Gene Expression kit."
    dataset_organism: "homo_sapiens"
    dataset_assay_spatial: Visium 
    dataset_assay_singlecell: Chromium

  - id: brain
    input_sc: resources/datasets_raw/MBRAIN_sc.rds
    input_sp: resources/datasets_raw/BRAIN.rds
    dataset_id: brain
    dataset_name: Brain
    dataset_url: https://github.com/BayraktarLab/cell2location
    dataset_reference: 10.1038/s41587-021-01139-4
    dataset_summary: 10X Visium spatial RNA-seq from adult mouse brain sections paired to single-nucleus RNA-seq
    dataset_description: "This datasets were generated matched single nucleus (sn, this submission) and Visium spatial RNA-seq (10X Genomics) profiles of adjacent mouse brain sections that contain multiple regions from the telencephalon and diencephalon."
    dataset_organism: "homo_sapiens"
    dataset_assay_spatial: Visium 
    dataset_assay_singlecell: Chromium

  - id: fibrosarcoma
    input_sc: resources/datasets_raw/MCATUMOR_sc.rds
    input_sp: resources/datasets_raw/MCATUMOR.rds
    dataset_id: fibrosarcoma
    dataset_name: Fibrosarcoma
    dataset_url: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE173778
    dataset_reference: 10.1038/s41587-022-01272-8
    dataset_summary: Multi-resolution deconvolution of spatial transcriptomics data reveals continuous patterns of inflammation
    dataset_description: "Murine lymph node and tumor; spatial transcriptomics and scRNA-seq data"
    dataset_organism: "mus_musculus"
    dataset_assay_spatial: Visium 
    dataset_assay_singlecell: Chromium

  - id: cortex
    input_sc: resources/datasets_raw/MCORTEX_sc.rds
    input_sp: resources/datasets_raw/MCORTEX.rds
    dataset_id: cortex
    dataset_name: Cortex
    dataset_url: https://zenodo.org/records/2669683
    dataset_reference: 10.1038/s41586-019-1049-y
    dataset_summary: Scripts and source data for image processing, barcode calling, and cell type annotations in a seqFISH+ experiment.
    dataset_description: "The dataset includes processed image data, cell type annotations with Louvain clusters, gene IDs for transcript locations, and mRNA point locations."
    dataset_organism: "homo_sapiens"
    dataset_assay_spatial: seqFISH+ 
    dataset_assay_singlecell: Chromium

  - id: gastrulation
    input_sc: resources/datasets_raw/MGASTRULA_sc.rds
    input_sp: resources/datasets_raw/MGASTRULA.rds
    dataset_id: gastrulation
    dataset_name: Gastrulation
    dataset_url: https://content.cruk.cam.ac.uk/jmlab/SpatialMouseAtlas2020/
    dataset_reference: 10.1038/s41587-021-01006-2
    dataset_summary: single-cell and spatial transcriptomic molecular map of mouse gastrulation 
    dataset_description: "Single-Cell omics Data across Mouse Gastrulation and Highly multiplexed spatially resolved gene expression profiling of Early Organogenesis."
    dataset_organism: "mus_musculus"
    dataset_assay_spatial: seqFISH 
    dataset_assay_singlecell: Chromium

  - id: olfactorybulb
    input_sc: resources/datasets_raw/MOBNEW_sc.rds
    input_sp: resources/datasets_raw/MOBNEW.rds
    dataset_id: olfactorybulb
    dataset_name: Olfactorybulb
    dataset_url: http://ww1.spatialtranscriptomicsresearch.org/?usid=24&utid=8672855942
    dataset_reference: 10.1126/science.aaf2403
    dataset_summary: Single-cell and spatial transcriptomic of mouse olfactory bulb
    dataset_description: "Single-cell and spatial transcriptomic of mouse olfactory bulb"
    dataset_organism: "mus_musculus"
    dataset_assay_spatial: ST 
    dataset_assay_singlecell: Chromium

  - id: hindlimbmuscle
    input_sc: resources/datasets_raw/MUSCLE_sc.rds
    input_sp: resources/datasets_raw/MUSCLE.rds
    dataset_id: hindlimbmuscle
    dataset_name: Hindlimbmuscle
    dataset_url: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE161318
    dataset_reference: 10.1101/2020.12.01.407460
    dataset_summary: Spatial RNA sequencing of regenerating mouse hindlimb muscle
    dataset_description: "The spatial transcriptomics datasets regenerates mouse muscle tissue generated with the 10x Genomics Visium platform."
    dataset_organism: "mus_musculus"
    dataset_assay_spatial: Visium 
    dataset_assay_singlecell: Chromium
  - id: pdac
    input_sc: resources/datasets_raw/PDAC_sc.rds
    input_sp: resources/datasets_raw/PDAC.rds
    dataset_id: pdac
    dataset_name: pancreatic ductal adenocarcinomas
    dataset_url: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE111672
    dataset_reference: 10.1101/2020.12.01.407460
    dataset_summary: Integrating microarray-based spatial transcriptomics and single-cell RNA-seq reveals tissue architecture in pancreatic ductal adenocarcinomas
    dataset_description: "We developed a multimodal intersection analysis method combining scRNA-seq with spatial transcriptomics to map and characterize the spatial organization and interactions of distinct cell subpopulations in complex tissues, such as primary pancreatic tumors.."
    dataset_organism: "homo_sapiens"
    dataset_assay_spatial: ST 
    dataset_assay_singlecell: Chromium

publish_dir: "resources/datasets"
output_sc: "$id/output_sc.h5ad"
output_sp: "$id/output_sp.h5ad"
output_state: "$id/state.yaml"
EOF

nextflow run . \
  -main-script target/nextflow/process_datasets/convert/main.nf \
  -profile docker \
  -params-file /tmp/params.yaml
