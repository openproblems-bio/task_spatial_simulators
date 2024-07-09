#!/bin/bash


cat > /tmp/params.yaml <<'EOF'
param_list:
  - id: BREAST
    input_sc: resources/datasets_raw/BREAST_sc.rds
    input_sp: resources/datasets_raw/BREAST.rds
    dataset_id: BREAST
    dataset_name: "Nicely formatted name"
    dataset_url: "https://example.com"
    dataset_reference: "doi identifier"
    dataset_summary: "A short description of the dataset (1 line)."
    dataset_description: "A longer description of the dataset."
    dataset_organism: "mus_musculus/homo_sapiens"
  - id: HOSTEOSARCOMA
    input_sc: resources/datasets_raw/HOSTEOSARCOMA_sc.rds
    input_sp: resources/datasets_raw/HOSTEOSARCOMA.rds
    dataset_id: HOSTEOSARCOMA
    dataset_name: "Nicely formatted name"
    dataset_url: "https://example.com"
    dataset_reference: "doi identifier"
    dataset_summary: "A short description of the dataset (1 line)."
    dataset_description: "A longer description of the dataset."
    dataset_organism: "mus_musculus/homo_sapiens"
  # todo: add more

publish_dir: "resources/datasets"
output_sc: "$id/output_sc.h5ad"
output_sp: "$id/output_sp.h5ad"
output_state: "$id/state.yaml"
EOF

nextflow run . \
  -main-script target/nextflow/process_datasets/convert/main.nf \
  -profile docker \
  -params-file /tmp/params.yaml
