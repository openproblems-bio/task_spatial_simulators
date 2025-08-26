include { findArgumentSchema } from "${meta.resources_dir}/helper.nf"

workflow auto {
  findStates(params, meta.config)
    | meta.workflow.run(
      auto: [publish: "state"]
    )
}

workflow run_wf {
  take:
  input_ch

  main:

  output_ch = input_ch

    | convert.run(
      fromState: [
        "input_sc",
        "input_sp",
        "dataset_id",
        "dataset_name",
        "dataset_url_spatial",
        "dataset_url_singlecell",
        "dataset_reference_spatial",
        "dataset_reference_singlecell",
        "dataset_summary_spatial",
        "dataset_summary_singlecell",
        "dataset_description_spatial",
        "dataset_description_singlecell",
        "dataset_organism",
        "dataset_assay_spatial",
        "dataset_assay_singlecell"
      ],
      toState: [
        "output_sc",
        "output_sp"
      ]
    )

    | sc_features.run(
      fromState: ["input_sp": "output_sp"],
      toState: ["output_sp"]
    )

    | precompute_downstream.run(
      fromState: ["input_sp": "output_sp"],
      toState: ["output_sp"]
    )

    | setState ([
      "output_sc",
      "output_sp",
      "_meta"
    ])

  emit:
  output_ch
}
