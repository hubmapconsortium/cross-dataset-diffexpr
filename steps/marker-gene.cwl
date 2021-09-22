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
  hdf5_file:
    type: File
    outputBinding:
      glob: "rna.hdf5"
    doc: hdf5 file with layers containing dataframes for cell, group, and quant data

  csv_file:
    type: File
    outputBinding:
      glob: "rna.csv"
    doc: csv file containing long narrow expression data

  mini_hdf5_file:
    type: File
    outputBinding:
      glob: "mini_rna.hdf5"

  mini_csv_file:
    type: File
    outputBinding:
      glob: "mini_rna.csv"