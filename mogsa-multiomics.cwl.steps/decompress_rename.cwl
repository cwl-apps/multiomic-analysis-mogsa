cwlVersion: v1.2
class: CommandLineTool
label: decompress+rename
$namespaces:
  sbg: https://sevenbridges.com

requirements:
- class: ShellCommandRequirement
- class: DockerRequirement
  dockerPull: images.sbgenomics.com/markop/sbg-decompressor:1.0
- class: InlineJavascriptRequirement

inputs:
- id: in_file
  type: File?

outputs:
- id: out_file
  type: File?
  outputBinding:
    glob: '*'

baseCommand: []
arguments:
- prefix: ''
  position: 0
  valueFrom: |
    ${ 
        if (inputs.in_file.nameext.endsWith(".gz")) {

        return "gzip -d " + inputs.in_file.path +  " && mv " + inputs.in_file.path.split(".").slice(0, -1).join(".") + " " + inputs.in_file.metadata.sample_id + ".txt"  
        }
        else {
        return "mv " + inputs.in_file.path + " " + inputs.in_file.metadata.sample_id + ".txt"
        } 
    } 
  shellQuote: false
id: david.roberson/build-multiomic-analysis-mogsa/decompress-rename/1
sbg:appVersion:
- v1.2
sbg:archived: true
sbg:content_hash: a707a4adf324eea0d2ed6e89e7cf18142c5f278f1901c28f203bbf09e9b8399ac
sbg:contributors:
- david.roberson
sbg:createdBy: david.roberson
sbg:createdOn: 1621628321
sbg:id: david.roberson/build-multiomic-analysis-mogsa/decompress-rename/1
sbg:image_url:
sbg:latestRevision: 1
sbg:modifiedBy: david.roberson
sbg:modifiedOn: 1621630774
sbg:project: david.roberson/build-multiomic-analysis-mogsa
sbg:projectName: 'BUILD: Multiomic Analysis MOGSA'
sbg:publisher: sbg
sbg:revision: 1
sbg:revisionNotes: ''
sbg:revisionsInfo:
- sbg:modifiedBy: david.roberson
  sbg:modifiedOn: 1621628321
  sbg:revision: 0
  sbg:revisionNotes: |-
    Copy of david.roberson/copy-of-breast-cancer-multi-omics-analysis/decompress-rename/3
- sbg:modifiedBy: david.roberson
  sbg:modifiedOn: 1621630774
  sbg:revision: 1
  sbg:revisionNotes: ''
sbg:sbgMaintained: false
sbg:validationErrors: []
