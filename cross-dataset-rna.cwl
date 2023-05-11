#!/usr/bin/env cwl-runner

class: Workflow
cwlVersion: v1.0
label: Pipeline for evaluating differential expression of genes across datasets

inputs:

  enable_manhole:
    label: "Whether to enable remote debugging via 'manhole'"
    type: boolean?

  data_directory:
    label: "List of paths to all processed RNA datasets"
    type: Directory

  uuids_file
    label: "Path to a file containing a list of uuids for the dataset to be indexed"
    type: File


outputs:

  concatenated_file:
    outputSource: annotate-concatenate/concatenated_annotated_file
    type: File
  h5ad_file:
    outputSource: annotate-concatenate/h5ad_file
    type: File
  precompute_file:
    outputSource: annotate-concatenate/precompute_file
    type: File
  pdf_files:
    outputSource: annotate-concatenate/pdf_files
    type: File[]
  hdf5_file:
    outputSource: annotate-concatenate/hdf5_file
    type: File
  precompute_file:
    outputSource: annotate-concatenate/precompute_file
    type: File
  zarr_file:
    outputSource: annotate-concatenate/zarr_store
    type: File

steps:

  - id: annotate-concatenate
    in:
      - id: data_directories
        source: data_directories
      - id: nexus_token
        source: nexus_token
      - id: uuids_file:
        source: uuids_file

    out:
      - hdf5_file
      - h5ad_file
      - precompute_file
      - pdf_files
      - zarr_file

    run: steps/annotate-concatenate.cwl
    label: "Annotates and concatenates h5ad data files in directory"