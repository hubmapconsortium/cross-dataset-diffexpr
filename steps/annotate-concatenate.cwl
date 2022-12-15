cwlVersion: v1.0
class: CommandLineTool
label: Annotates each h5ad file with dataset and tissue type, then concatenates

hints:
  DockerRequirement:
    dockerPull: ghcr.io/hubmapconsortium/cross-dataset-scanpy
baseCommand: /opt/annotate_concatenate.py

inputs:
  enable_manhole:
    label: "Whether to enable remote debugging via 'manhole'"
    type: boolean?
    inputBinding:
      position: 0

  nexus_token:
    type: string?
    doc: Valid nexus token for search-api
    inputBinding:
      position: 1
    default: "None"

  data_directories:
    type: Directory[]
    doc: List of paths to processed dataset directories
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

  precompute_file:
    type: File
    outputBinding:
      glob: "rna_precompute.hdf5"
    doc: hdf5 file with data for accelerated queries

  h5ad_file:
    type: File
    outputBinding:
      glob: "rna.h5ad"
    doc: h5ad file containing expression data

  zarr_file:
    type: File
    outputBinding
      glob: "rna.zarr"
    doc: zarr store containing batch corrected data, metadata, cluster assignments and UMAP coordinates

  pdf_files:
    type: File[]
    outputBinding:
      glob: "*.pdf"