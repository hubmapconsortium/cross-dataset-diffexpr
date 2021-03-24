cwlVersion: v1.0
class: CommandLineTool
label: Batch correction, dimensionality reduction and clustering
hints:
  DockerRequirement:
    dockerPull: docker.pkg.github.com/hubmapconsortium/cross-dataset-diffexpr/cross-dataset-scanpy:latest
baseCommand: /opt/batch_correct_umap_cluster.py

inputs:
  concatenated_annotated_file:
    type: File
    inputBinding:
      position: 1
outputs:
  h5ad_file:
    type: File
    outputBinding:
      glob: "*bc_umap_cluster.h5ad"

  pdf_files:
    type: File[]
    outputBinding:
      glob: "*.pdf"
