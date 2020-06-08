cwlVersion: v1.0
class: CommandLineTool
label: Dimensionality reduction and clustering
hints:
  DockerRequirement:
    dockerPull: hubmap/cross-dataset-scanpy:latest
baseCommand: /opt/scanpy_entry_point.py

inputs:
  concatenated_file:
    type: File
    inputBinding:
      position: 1
outputs:
  h5ad_files:
    type: File[]
    outputBinding:
      glob: "*.h5ad"

  pdf_files:
    type: File[]
    outputBinding:
      glob: "*.pdf"

  marker_gene_database:
    type: File
    outputBinding:
      glob: marker_genes.db
