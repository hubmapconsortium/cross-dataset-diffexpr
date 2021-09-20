cwlVersion: v1.0
class: CommandLineTool
label: Finds marker genes associated with different groupings in atac-seq data

hints:
  DockerRequirement:
    dockerPull: hubmap/cross-dataset-scanpy:latest
baseCommand: /opt/precompute_percentages.py

inputs:
  batch_corrected_file:
    type: File
    doc: h5ad file containing batch corrected RNA seq data
    inputBinding:
      position: 1

outputs:
  hdf5_file:
    type: File
    outputBinding:
      glob: "rna_precompute.hdf5"
    doc: hdf5 file with layers containing dataframes for genes and precomputed percentages
