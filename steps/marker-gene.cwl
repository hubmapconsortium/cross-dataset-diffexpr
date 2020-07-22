cwlVersion: v1.0
class: CommandLineTool
label: Marker gene analysis
hints:
  DockerRequirement:
    dockerPull: hubmap/cross-dataset-scanpy:latest
baseCommand: /opt/find_marker_genes.py

inputs:
  batch_corrected_file:
    type: File
    inputBinding:
      position: 1
outputs:
  csv_files:
    type: File[]
    outputBinding:
      glob: "*.csv"
