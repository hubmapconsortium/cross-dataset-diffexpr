cwlVersion: v1.0
class: CommandLineTool
label: Marker gene analysis
hints:
  DockerRequirement:
    dockerPull: ghcr.io/hubmapconsortium/cross-dataset-scanpy
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
  mini_hdf5_file:
    type: File
    outputBinding:
      glob: "mini_rna.hdf5"

  mini_csv_file:
    type: File
    outputBinding:
      glob: "mini_rna.csv"