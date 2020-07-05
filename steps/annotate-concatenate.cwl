cwlVersion: v1.0
class: CommandLineTool
label: Annotates each h5ad file with dataset and tissue type, then concatenates

hints:
  DockerRequirement:
    dockerPull: hubmap/cross-dataset-scanpy:local-test
baseCommand: /opt/annotate_concatenate.py

inputs:
  data_dir:
    type: Directory
    doc: Base directory to be recursively searched for h5ad files to be annotated and concatenated
    inputBinding:
      position: 1

outputs:
  concatenated_file:
    type: File
    outputBinding:
      glob: "concatenated_annotated_data.h5ad"
    doc: Annotated, concatenated dataset in hdf5 format
