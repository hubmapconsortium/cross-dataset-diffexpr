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

  enable_manhole:
    label: "Whether to enable remote debugging via 'manhole'"
    type: boolean?
    inputBinding:
      position: 0

outputs:
  concatenated_annotated_file:
    type: File
    outputBinding:
      glob: "concatenated_annotated_data.h5ad"
    doc: Annotated, concatenated dataset in hdf5 format

  old_cluster_file:
    type: File
    outputBinding:
      glob: "cluster.hdf5"
    doc: Hdf file containing old cluster p values

  gene_dictionaries:
    type: File[]
    outputBinding:
      glob: '*.json'
    doc: Json files mapping from hugo symbols to versioned ensembl ids