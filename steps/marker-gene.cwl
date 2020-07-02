cwlVersion: v1.0
class: CommandLineTool
label: Marker gene analysis
hints:
  DockerRequirement:
    dockerPull: hubmap/cross-dataset-scanpy:latest
baseCommand: /opt/scanpy_entry_point.py

inputs:
  batch_corrected_file:
    type: File
    inputBinding:
      position: 1
outputs:
  h5ad_files:
    type: File[]
    outputBinding:
      glob: "*.h5ad"
