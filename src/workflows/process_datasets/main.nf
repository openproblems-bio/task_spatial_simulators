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
        "dataset_url",
        "dataset_reference",
        "dataset_summary",
        "dataset_description",
        "dataset_organism"
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
