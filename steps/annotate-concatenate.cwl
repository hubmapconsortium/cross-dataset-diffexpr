cwlVersion: v1.0
class: CommandLineTool
label: Annotates each h5ad file with dataset and tissue type, then concatenates

hints:
  DockerRequirement:
    dockerPull: hubmap/cross-dataset-scanpy:latest
baseCommand: /opt/annotate_concatenate.py

inputs:
  nexus_token:
    type: string
    doc: Valid nexus token for search-api
    inputBinding:
      position: 1

  data_directories:
    type: string[]
    doc: List of paths to processed dataset directories
    inputBinding:
      position: 2

outputs:
  concatenated_file:
    type: File
    outputBinding:
      glob: "concatenated_annotated_data.h5ad"
    doc: Annotated, concatenated dataset in hdf5 format
