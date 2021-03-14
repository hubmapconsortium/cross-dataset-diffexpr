#!/usr/bin/env cwl-runner

class: Workflow
cwlVersion: v1.0
label: Pipeline for evaluating differential expression of genes across datasets

inputs:

  data_directories:
    label: "List of paths to all processed RNA datasets"
    type: Directory[]

  nexus_token:
    label: "Valid nexus token for search-api"
    type: string?


outputs:

  concatenated_file:
    outputSource: annotate-concatenate/concatenated_annotated_file
    type: File
  h5ad_file:
    outputSource: batch-correct/h5ad_file
    type: File
  pdf_files:
    outputSource: batch-correct/pdf_files
    type: File[]
  hdf5_file:
    outputSource: marker-gene/hdf5_file
    type: File
  csv_file:
    outputSource: marker-gene/csv_file
    type: File
  gene_dictionaries:
    outputSource: annotate-concatenate/gene_dictionaries
    type: File[]
  mini_hdf5_file:
    outputSource: marker-gene/mini_hdf5_file
    type: File
  mini_csv_file:
    outputSource: marker-gene/mini_csv_file
    type: File

steps:

  - id: annotate-concatenate
    in:
      - id: data_directories
        source: data_directories
      - id: nexus_token
        source: nexus_token

    out:
      - concatenated_annotated_file
      - gene_dictionaries
      - old_cluster_file

    run: steps/annotate-concatenate.cwl
    label: "Annotates and concatenates h5ad data files in directory"

  - id: batch-correct
    in:
      - id: concatenated_annotated_file
        source: annotate-concatenate/concatenated_annotated_file
    out:
      - h5ad_file
      - pdf_files

    run: steps/batch-correct.cwl
    label: "Cross dataset secondary analysis via ScanPy, including batch correction, dimensionality reduction and leiden clustering"

  - id: marker-gene
    in:
      - id: batch_corrected_file
        source: batch-correct/h5ad_file
      - id: old_cluster_file
        source: annotate-concatenate/old_cluster_file
    out:
      - hdf5_file
      - csv_file
      - mini_hdf5_file
      - mini_csv_file

    run: steps/marker-gene.cwl
    label: "Cross dataset secondary analysis via ScanPy, including marker gene analysis"
