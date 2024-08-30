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
    dataset_url_spatial: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176078
    dataset_summary_spatial: A spatially resolved atlas of human breast cancers
    dataset_description_spatial: "This study presents a spatially resolved transcriptomics analysis of human breast cancers."
    dataset_reference_spatial: 10.1038/s41588-021-00911-1 
    dataset_url_singlecell: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176078
    dataset_summary_singlecell: A single-cell of human breast cancers
    dataset_description_singlecell: "This study presents a single cell analysis of human breast cancers."
    dataset_reference_singlecell: 10.1038/s41588-021-00911-1
    dataset_organism: "homo_sapiens"
    dataset_assay_spatial: Visium # EFO:0010961
    dataset_assay_singlecell: Chromium
  
  - id: osteosarcoma
    input_sc: resources/datasets_raw/HOSTEOSARCOMA_sc.rds
    input_sp: resources/datasets_raw/HOSTEOSARCOMA.rds
    dataset_id: osteosarcoma
    dataset_name: Osteosarcoma
    dataset_url_spatial: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176078 
    dataset_summary_spatial: Spatial profiling of human osteosarcoma cells.
    dataset_description_spatial: "Spatial transcriptome profiling by MERFISH reveals subcellular RNA compartmentalization and cell cycle-dependent gene expression."
    dataset_reference_spatial: 10.1073/pnas.1912459116
    dataset_url_singlecell: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152048
    dataset_summary_singlecell: Single cell analysis of osteosarcoma tissues
    dataset_description_singlecell: "We performed the scRNA-seq anakysis of BC22 osteosarcoma tissues based on the 10X Genomics platform."
    dataset_reference_singlecell: 10.1038/s41467-020-20059-6
    dataset_organism: "homo_sapiens"
    dataset_assay_spatial: Visium 
    dataset_assay_singlecell: Chromium

  - id: prostate
    input_sc: resources/datasets_raw/HPROSTATE_sc.rds
    input_sp: resources/datasets_raw/HPROSTATE.rds
    dataset_id: prostate
    dataset_name: Prostate
    dataset_url_spatial: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE159697
    dataset_summary_spatial: Spatially resolved gene expression of human protate tissue slices treated with steroid hormones for 8 hours
    dataset_description_spatial: "Spatially resolved gene expression was prepard by dissociated hman prostate tissue to single cells, and collected & prepped for RNA-seq using the Visium Spatial Gene Expression kit."
    dataset_reference_spatial: 10.1016/j.isci.2021.102640 
    dataset_url_singlecell: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE142489
    dataset_summary_singlecell: Single cell analysis of primary prostate human epithelial organoids grown in vehicle or 1,25D to day 8
    dataset_description_singlecell: "	Primary prostate epithelial cells derived from a single patient were grown as organoids in vehicle or 10nM 1,25D (vitamin D). At day 8 and day 14 of culture, organoids were dissociated to single cells, and collected & prepped for RNA-seq using the 10x Genomics Chromium Genome Reagent kit v3 Chemistry."
    dataset_reference_singlecell: 10.1016/j.isci.2021.102640
    dataset_organism: "homo_sapiens"
    dataset_assay_spatial: Visium 
    dataset_assay_singlecell: Chromium

  - id: brain
    input_sc: resources/datasets_raw/MBRAIN_sc.rds
    input_sp: resources/datasets_raw/BRAIN.rds
    dataset_id: brain
    dataset_name: Brain
    dataset_url_spatial: https://github.com/BayraktarLab/cell2location
    dataset_reference_spatial: 10.1038/s41587-021-01139-4
    dataset_summary_spatial: 10X Visium spatial RNA-seq from adult mouse brain sections paired to single-nucleus RNA-seq
    dataset_description_spatial: "This datasets were generated matched single nucleus (sn, this submission) and Visium spatial RNA-seq (10X Genomics) profiles of adjacent mouse brain sections that contain multiple regions from the telencephalon and diencephalon."
    dataset_url_singlecell: https://github.com/BayraktarLab/cell2location
    dataset_reference_singlecell: 10.1038/s41587-021-01139-4
    dataset_summary_singlecell: 10X Visium spatial RNA-seq from adult mouse brain sections paired to single-nucleus RNA-seq
    dataset_description_singlecell: "This datasets were generated matched single nucleus (sn, this submission) and Visium spatial RNA-seq (10X Genomics) profiles of adjacent mouse brain sections that contain multiple regions from the telencephalon and diencephalon."
    dataset_organism: "homo_sapiens"
    dataset_assay_spatial: Visium 
    dataset_assay_singlecell: Chromium

  - id: fibrosarcoma
    input_sc: resources/datasets_raw/MCATUMOR_sc.rds
    input_sp: resources/datasets_raw/MCATUMOR.rds
    dataset_id: fibrosarcoma
    dataset_name: Fibrosarcoma
    dataset_url_spatial: https://github.com/romain-lopez/DestVI-reproducibility
    dataset_reference_spatial: 10.1038/s41587-022-01272-8
    dataset_summary_spatial: Multi-resolution deconvolution of spatial transcriptomics data reveals continuous patterns of Tumor A1 of Tissue 1 
    dataset_description_spatial: "Spatial transcriptomics of Tumor A1 of Tissue 1."
    dataset_url_singlecell: https://github.com/romain-lopez/DestVI-reproducibility
    dataset_reference_singlecell: 10.1038/s41587-022-01272-8
    dataset_summary_singlecell: Multi-resolution deconvolution of spatial transcriptomics data reveals continuous patterns of Tumor A1 of Tissue 1 
    dataset_description_singlecell: "Spatial transcriptomics of Tumor A1 of Tissue 1."
    dataset_organism: "mus_musculus"
    dataset_assay_spatial: Visium 
    dataset_assay_singlecell: Chromium

  - id: cortex
    input_sc: resources/datasets_raw/MCORTEX_sc.rds
    input_sp: resources/datasets_raw/MCORTEX.rds
    dataset_id: cortex
    dataset_name: Cortex
    dataset_url_spatial: https://zenodo.org/records/2669683
    dataset_reference_spatial: 10.1038/s41586-019-1049-y
    dataset_summary_spatial: Scripts and source data for image processing, barcode calling, and cell type annotations in a seqFISH+ experiment.
    dataset_description_spatial: "The dataset includes processed image data, cell type annotations with Louvain clusters, gene IDs for transcript locations, and mRNA point locations, with additional data available on Zenodo."
    dataset_url_singlecell: https://portal.brain-map.org/atlases-and-data/rnaseq/mouse-v1-and-alm-smart-seq
    dataset_reference_singlecell: 10.1038/s41586-019-1049-y
    dataset_summary_singlecell: The single cell of mouse primary visual cortex (VISp).
    dataset_description_singlecell: "To investigate the diversity of cell types across the mouse neocortex at the transcriptional level, 15,413 cells were profiled from the primary visual cortex (VISp), and 10,068 cells from the anterior lateral motor cortex (ALM).."
    dataset_organism: "homo_sapiens"
    dataset_assay_spatial: seqFISH+ 
    dataset_assay_singlecell: Smart-seq

  - id: gastrulation
    input_sc: resources/datasets_raw/MGASTRULA_sc.rds
    input_sp: resources/datasets_raw/MGASTRULA.rds
    dataset_id: gastrulation
    dataset_name: Gastrulation
    dataset_url_spatial: https://content.cruk.cam.ac.uk/jmlab/SpatialMouseAtlas2020/
    dataset_reference_spatial: 10.1038/s41587-021-01006-2
    dataset_summary_spatial: single-cell and spatial transcriptomic molecular map of mouse gastrulation 
    dataset_description_spatial: "Single-Cell omics Data across Mouse Gastrulation and Highly multiplexed spatially resolved gene expression profiling of Early Organogenesis."
    dataset_url_singlecell: https://content.cruk.cam.ac.uk/jmlab/SpatialMouseAtlas2020/
    dataset_reference_singlecell: 10.1038/s41587-021-01006-2
    dataset_summary_singlecell: single-cell and spatial transcriptomic molecular map of mouse gastrulation 
    dataset_description_singlecell: "Single-Cell omics Data across Mouse Gastrulation and Highly multiplexed spatially resolved gene expression profiling of Early Organogenesis."
    dataset_organism: "mus_musculus"
    dataset_assay_spatial: seqFISH 
    dataset_assay_singlecell: Chromium

  - id: olfactorybulb
    input_sc: resources/datasets_raw/MOBNEW_sc.rds
    input_sp: resources/datasets_raw/MOBNEW.rds
    dataset_id: olfactorybulb
    dataset_name: Olfactorybulb
    dataset_url_spatial: http://ww1.spatialtranscriptomicsresearch.org/?usid=24&utid=8672855942
    dataset_reference_spatial: 10.1126/science.aaf2403
    dataset_summary_spatial: Single-cell and spatial transcriptomic of mouse olfactory bulb
    dataset_description_spatial: "Single-cell and spatial transcriptomic of mouse olfactory bulb"
    dataset_url_singlecell: http://ww1.spatialtranscriptomicsresearch.org/?usid=24&utid=8672855942
    dataset_reference_singlecell: 10.1126/science.aaf2403
    dataset_summary_singlecell: Single-cell and spatial transcriptomic of mouse olfactory bulb
    dataset_description_singlecell: "Single-cell and spatial transcriptomic of mouse olfactory bulb"
    dataset_organism: "mus_musculus"
    dataset_assay_spatial: ST 
    dataset_assay_singlecell: Chromium

  - id: hindlimbmuscle
    input_sc: resources/datasets_raw/MUSCLE_sc.rds
    input_sp: resources/datasets_raw/MUSCLE.rds
    dataset_id: hindlimbmuscle
    dataset_name: Hindlimbmuscle
    dataset_url_spatial: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE161318
    dataset_reference_spatial: 10.1101/2020.12.01.407460
    dataset_summary_spatial: Spatial RNA sequencing of regenerating mouse hindlimb muscle
    dataset_description_spatial: "The spatial transcriptomics datasets regenerates mouse muscle tissue generated with the 10x Genomics Visium platform."
    dataset_url_singlecell: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE159500
    dataset_reference_singlecell: 10.1038/s42003-021-02810-x
    dataset_summary_singlecell: Single-cell transcriptomic sampling of regenerating mouse muscle tissue
    dataset_description_singlecell: "The overall objective of the study was to generate a high-resolution transcriptomic reference of the distinct cell populations in the mouse muscle at different timepoints after notexin-injury."
    dataset_organism: "mus_musculus"
    dataset_assay_spatial: Visium 
    dataset_assay_singlecell: Chromium

  - id: pdac
    input_sc: resources/datasets_raw/PDAC_sc.rds
    input_sp: resources/datasets_raw/PDAC.rds
    dataset_id: pdac
    dataset_name: pancreatic ductal adenocarcinomas
    dataset_url_spatial: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE111672
    dataset_reference_spatial: 10.1101/2020.12.01.407460
    dataset_summary_spatial: Integrating microarray-based spatial transcriptomics and single-cell RNA-seq reveals tissue architecture in pancreatic ductal adenocarcinomas
    dataset_description_spatial: "We developed a multimodal intersection analysis method combining scRNA-seq with spatial transcriptomics to map and characterize the spatial organization and interactions of distinct cell subpopulations in complex tissues, such as primary pancreatic tumors.."
    dataset_url_singlecell: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE111672
    dataset_reference_singlecell: 10.1101/2020.12.01.407460
    dataset_summary_singlecell: Integrating microarray-based spatial transcriptomics and single-cell RNA-seq reveals tissue architecture in pancreatic ductal adenocarcinomas
    dataset_description_singlecell: "We developed a multimodal intersection analysis method combining scRNA-seq with spatial transcriptomics to map and characterize the spatial organization and interactions of distinct cell subpopulations in complex tissues, such as primary pancreatic tumors.."
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
