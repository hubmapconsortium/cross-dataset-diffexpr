#!/usr/bin/env cwl-runner

class: Workflow
cwlVersion: v1.0
label: Pipeline for evaluating differential expression of genes across datasets

inputs:

  data_dir_log:
    label: "Text file containing paths to all processed RNA datasets"
    type: File

  nexus_token:
    label: "Valid nexus token for search-api"
    type: String

outputs:

  concatenated_file:
    outputSource: annotate-concatenate/concatenated_file
    type: File
  h5ad_file:
    outputSource: batch-correct/h5ad_file
    type: File
  pdf_files:
    outputSource: batch-correct/pdf_files
    type: File[]
  csv_files:
    outputSource: marker-gene/csv_files
    type: File[]

steps:

  - id: read-data-dir-log:
    in:
      - id: data_dir_log
        source: data_dir_log
    out:
      - data_directories
    run: steps/read-data-dir-log.cwl
    label: "Reads the log containing processed datasets"

  - id: annotate-concatenate
    in:
      - id: data_directories
        source: read-data-dir-log/data_directories
      - id: nexus_token
        source: nexus_token

    out:
      - concatenated_file

    run: steps/annotate-concatenate.cwl
    label: "Annotates and concatenates h5ad data files in directory"

  - id: batch-correct
    in:
      - id: concatenated_file
        source: annotate-concatenate/concatenated_file
    out:
      - h5ad_file
      - pdf_files

    run: steps/batch-correct.cwl
    label: "Cross dataset secondary analysis via ScanPy, including batch correction, dimensionality reduction and leiden clustering"

  - id: marker-gene
    in:
      - id: batch_corrected_file
        source: batch-correct/h5ad_file
    out:
      - csv_files

    run: steps/marker-gene.cwl
    label: "Cross dataset secondary analysis via ScanPy, including marker gene analysis"
