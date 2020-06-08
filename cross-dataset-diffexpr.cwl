#!/usr/bin/env cwl-runner

class: Workflow
cwlVersion: v1.0
label: Pipeline for evaluating differential expression of genes across datasets

inputs:
  data_dir:
    label: "Directory containing h5ad data files"
    type: Directory

outputs:
  h5ad_files:
    outputSource: scanpy_analysis/h5ad_files
    type: File[]
  pdf_files:
    outputSource: scanpy_analysis/pdf_files
    type: File[]
  marker_gene_database:
    outputSource: scanpy_analysis/marker_gene_database
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

  - id: scanpy_analysis
    in:
      - id: concatenated_file
        source: annotate-concatenate/concatenated_file
    out:
      - h5ad_files
      - pdf_files
      - marker_gene_database

    run: steps/scanpy-analysis.cwl
    label: "Cross dataset secondary analysis via ScanPy"
