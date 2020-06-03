#!/usr/bin/env cwl-runner

class: Workflow
cwlVersion: v1.0
label: Pipeline for evaluating differential expression of genes across datasets

inputs:
  data_dir:
    label: "Directory containing h5ad data files"
    type: Directory

outputs:

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
      - qc_results
      - filtered_data
      - umap_pdf
      - marker_gene_plot_t_test
      - marker_gene_plot_logreg
    run: steps/scanpy-analysis.cwl
    label: "Secondary analysis via ScanPy"
