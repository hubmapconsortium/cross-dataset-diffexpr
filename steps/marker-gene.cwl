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

  old_cluster_file:
    type: File
    inputBinding:
      position: 2

outputs:
  hdf5_file:
    type: File
    outputBinding:
      glob: "rna.hdf5"

  csv_file:
    type: File
    outputBinding:
      glob: "rna.csv"

  mini_hdf5_file:
    type: File
    outputBinding:
      glob: "mini_rna.hdf5"

  mini_csv_file:
    type: File
    outputBinding:
      glob: "mini_rna.csv"