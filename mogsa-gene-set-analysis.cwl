cwlVersion: v1.1
class: Workflow
label: MOGSA Gene Set Analysis
$namespaces:
  sbg: https://sevenbridges.com

requirements: []

inputs: []

outputs:
- id: R_workspace
  type: File?
  outputSource:
  - mogsa/R_workspace
  sbg:x: 83.69499969482422
  sbg:y: -326.49700927734375
- id: html
  type: File?
  outputSource:
  - mogsa/html
  sbg:x: 82.69499969482422
  sbg:y: -152.4970245361328

steps:
- id: mogsa
  label: MOGSA
  in: []
  run: mogsa-gene-set-analysis.cwl.steps/mogsa.cwl
  out:
  - id: R_workspace
  - id: html
  sbg:x: -116
  sbg:y: -238
sbg:appVersion:
- v1.1
sbg:content_hash: a047f4621477fcd43893bda4ce2617475830d43ba91cb7117f556c7d6159a04c9
sbg:contributors:
- david.roberson
sbg:createdBy: david.roberson
sbg:createdOn: 1619459415
sbg:id: david.roberson/build-multiomic-analysis-mogsa/mogsa-gene-set-analysis/1
sbg:image_url: |-
  https://cgc.sbgenomics.com/ns/brood/images/david.roberson/build-multiomic-analysis-mogsa/mogsa-gene-set-analysis/1.png
sbg:latestRevision: 1
sbg:modifiedBy: david.roberson
sbg:modifiedOn: 1619459475
sbg:original_source: |-
  https://cgc-api.sbgenomics.com/v2/apps/david.roberson/build-multiomic-analysis-mogsa/mogsa-gene-set-analysis/1/raw/
sbg:project: david.roberson/build-multiomic-analysis-mogsa
sbg:projectName: 'BUILD: Multiomic Analysis MOGSA'
sbg:publisher: sbg
sbg:revision: 1
sbg:revisionNotes: ''
sbg:revisionsInfo:
- sbg:modifiedBy: david.roberson
  sbg:modifiedOn: 1619459415
  sbg:revision: 0
  sbg:revisionNotes: Copy of david.roberson/mogsa-workflow-wrapping/mogsa-gene-set-analysis/2
- sbg:modifiedBy: david.roberson
  sbg:modifiedOn: 1619459475
  sbg:revision: 1
  sbg:revisionNotes: ''
sbg:sbgMaintained: false
sbg:validationErrors: []
