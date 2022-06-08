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

  access_key_id:
    type: str
    inputBinding:
      position: 3

  secret_access_key:
    type: str
    inputBinding:
      position: 4


outputs:
  hdf5_file:
    type: File
    outputBinding:
      glob: "rna.hdf5"
    doc: hdf5 file with layers containing dataframes for cell, group, and gene data

  h5ad_file:
    type: File
    outputBinding:
      glob: "rna.h5ad"
    doc: h5ad file containing expression data

  precompute_file:
    type: File
    outputBinding:
      glob: "rna_precompute.hdf5"
    doc: hdf5 file containing precomputation results for accelerated queries