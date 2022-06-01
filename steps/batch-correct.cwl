cwlVersion: v1.0
class: CommandLineTool
label: Batch correction, dimensionality reduction and clustering
hints:
  DockerRequirement:
    dockerPull: ghcr.io/hubmapconsortium/cross-dataset-scanpy
baseCommand: /opt/batch_correct_umap_cluster.py

inputs:
  enable_manhole:
    label: "Whether to enable remote debugging via 'manhole'"
    type: boolean?
    inputBinding:
      position: 0
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
