#!/usr/bin/env cwl-runner

class: Workflow
cwlVersion: v1.0
label: Pipeline for evaluating differential expression of genes across datasets

inputs:

  enable_manhole:
    label: "Whether to enable remote debugging via 'manhole'"
    type: boolean?

  data_directories:
    label: "List of paths to all processed RNA datasets"
    type: Directory[]

  nexus_token:
    label: "Valid nexus token for search-api"
    type: string?

  access_key_id:
    label: "String id containing credentials for writing to s3 bucket"
    type: string

  secret_access_key:
    label: "String key containing credentials for writing to s3 bucket"
    type: string


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
  gene_dictionaries:
    outputSource: annotate-concatenate/gene_dictionaries
    type: File[]

steps:

  - id: annotate-concatenate
    in:
      - id: data_directories
        source: data_directories
      - id: nexus_token
        source: nexus_token
      - id: enable_manhole
        source: enable_manhole
      - id: access_key_id
        source: access_key_id
      - id: secret_access_key
        source: secret_access_key

    out:
      - concatenated_annotated_file
      - gene_dictionaries
      - old_cluster_file
      - hdf5_file
      - h5ad_file
      - precompute_file

    run: steps/annotate-concatenate.cwl
    label: "Annotates and concatenates h5ad data files in directory"