#!/usr/bin/env cwl-runner

class: Workflow
cwlVersion: v1.0
label: Pipeline for evaluating differential expression of genes across datasets

inputs:
  data_dir:
    label: "Directory containing h5ad data files"
    type: Directory

outputs:
  h5ad_file:
    outputSource: batch-correct/h5ad_file
    type: File
  pdf_files:
    outputSource: batch-correct/pdf_files
    type: File[]
  db_file:
    outputSource: marker-gene/db_file
    type: File

steps:
  - id: annotate-concatenate
    in:
      - id: data_dir
        source: data_dir

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
      - db_file

    run: steps/marker-gene.cwl
    label: "Cross dataset secondary analysis via ScanPy, including marker gene analysis"
