cwlVersion: v1.0
class: CommandLineTool
label: Marker gene analysis
hints:
  DockerRequirement:
    dockerPull: ghcr.io/hubmapconsortium/cross-dataset-scanpy
baseCommand: /opt/find_marker_genes.py

inputs:

  enable_manhole:
    label: "Whether to enable remote debugging via 'manhole'"
    type: boolean?
    inputBinding:
      position: 0

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
    doc: hdf5 file with layers containing dataframes for cell, group, and gene data

  had_file:
    type: File
    outputBinding:
      glob: "rna.h5ad"
    doc: h5ad file containing expression data

  mini_hdf5_file:
    type: File
    outputBinding:
      glob: "mini_rna.hdf5"

  mini_csv_file:
    type: File
    outputBinding:
      glob: "mini_rna.csv"