cwlVersion: v1.2
class: Workflow
label: mogsa-multiomics
$namespaces:
  sbg: https://sevenbridges.com

requirements:
- class: ScatterFeatureRequirement
- class: MultipleInputFeatureRequirement
- class: InlineJavascriptRequirement
- class: StepInputExpressionRequirement

inputs:
- id: mogsa_database
  type: File?
  sbg:x: -226
  sbg:y: -254
- id: protein_abundances
  type: File?
  sbg:x: -703
  sbg:y: 53
- id: RNA_ProtGlobal_Non_Phospho
  type: File
  sbg:x: -467
  sbg:y: -324
- id: mrna_files
  label: mRNA files
  type: File[]?
  sbg:x: -900.0237426757812
  sbg:y: -203
- id: custom_input
  label: RUN MOGSA
  doc: RUN untested downstream apps?
  type: boolean
  sbg:x: -559
  sbg:y: 264

outputs:
- id: R_workspace
  type: File?
  outputSource:
  - mogsa/R_workspace
  sbg:x: 409.23529052734375
  sbg:y: -108.63827514648438
- id: html
  type: File?
  outputSource:
  - mogsa/html
  sbg:x: 402.213134765625
  sbg:y: 78.38286590576172

steps:
- id: determine_pcs_mfa
  label: determine-pcs-mfa
  in:
  - id: RNA_ProtGlobal_Non_Phospho
    source: RNA_ProtGlobal_Non_Phospho
  - id: input_files
    source:
    - prepare_protein/cleaned_data
    - prepare_mrna/mRNA_table
  - id: mogsa_database
    source: mogsa_database
  - id: custom_input
    source: custom_input
  run: mogsa-multiomics.cwl.steps/determine_pcs_mfa.cwl
  when: $(inputs.custom_input)
  out:
  - id: output
  sbg:x: -67
  sbg:y: -26
- id: prepare_mrna
  label: prepare-mrna
  in:
  - id: mRNA_matrix
    source:
    - decompress_rename/out_file
    linkMerge: merge_flattened
  run: mogsa-multiomics.cwl.steps/prepare_mrna.cwl
  out:
  - id: mRNA_table
  sbg:x: -494
  sbg:y: -140
- id: prepare_protein
  label: prepare-protein
  in:
  - id: protein_abundances
    source: protein_abundances
  - id: custom_input
    source: custom_input
  run: mogsa-multiomics.cwl.steps/prepare_protein.cwl
  when: $(inputs.custom_input)
  out:
  - id: cleaned_data
  sbg:x: -355
  sbg:y: 34
- id: mogsa
  label: MOGSA
  in:
  - id: mfa_dir
    loadListing: deep_listing
    source: determine_pcs_mfa/output
  - id: mogsa_database
    source: mogsa_database
  - id: custom_input
    source: custom_input
  run: mogsa-multiomics.cwl.steps/mogsa.cwl
  when: $(inputs.custom_input)
  out:
  - id: R_workspace
  - id: html
  sbg:x: 176.1913604736328
  sbg:y: 1
- id: decompress_rename
  label: decompress+rename
  in:
  - id: in_file
    source: mrna_files
  scatter:
  - in_file
  run: mogsa-multiomics.cwl.steps/decompress_rename.cwl
  out:
  - id: out_file
  sbg:x: -655
  sbg:y: -158
sbg:appVersion:
- v1.2
- v1.1
sbg:content_hash: aa027a932ed85ae086e89e7289f4d06ac32f38ff9f8b304a33cf2216bda0bed92
sbg:contributors:
- david.roberson
sbg:createdBy: david.roberson
sbg:createdOn: 1619795588
sbg:id: david.roberson/build-multiomic-analysis-mogsa/mogsa-multiomics/11
sbg:image_url:
sbg:latestRevision: 11
sbg:modifiedBy: david.roberson
sbg:modifiedOn: 1621630817
sbg:original_source: |-
  https://cgc-api.sbgenomics.com/v2/apps/david.roberson/build-multiomic-analysis-mogsa/mogsa-multiomics/11/raw/
sbg:project: david.roberson/build-multiomic-analysis-mogsa
sbg:projectName: 'BUILD: Multiomic Analysis MOGSA'
sbg:publisher: sbg
sbg:revision: 11
sbg:revisionNotes: ''
sbg:revisionsInfo:
- sbg:modifiedBy: david.roberson
  sbg:modifiedOn: 1619795588
  sbg:revision: 0
  sbg:revisionNotes:
- sbg:modifiedBy: david.roberson
  sbg:modifiedOn: 1619795686
  sbg:revision: 1
  sbg:revisionNotes: ''
- sbg:modifiedBy: david.roberson
  sbg:modifiedOn: 1619795781
  sbg:revision: 2
  sbg:revisionNotes: ''
- sbg:modifiedBy: david.roberson
  sbg:modifiedOn: 1619796391
  sbg:revision: 3
  sbg:revisionNotes: ''
- sbg:modifiedBy: david.roberson
  sbg:modifiedOn: 1620238867
  sbg:revision: 4
  sbg:revisionNotes: ''
- sbg:modifiedBy: david.roberson
  sbg:modifiedOn: 1620246730
  sbg:revision: 5
  sbg:revisionNotes: ''
- sbg:modifiedBy: david.roberson
  sbg:modifiedOn: 1621628348
  sbg:revision: 6
  sbg:revisionNotes: ''
- sbg:modifiedBy: david.roberson
  sbg:modifiedOn: 1621628608
  sbg:revision: 7
  sbg:revisionNotes: ''
- sbg:modifiedBy: david.roberson
  sbg:modifiedOn: 1621630206
  sbg:revision: 8
  sbg:revisionNotes: ''
- sbg:modifiedBy: david.roberson
  sbg:modifiedOn: 1621630595
  sbg:revision: 9
  sbg:revisionNotes: '"head"'
- sbg:modifiedBy: david.roberson
  sbg:modifiedOn: 1621630604
  sbg:revision: 10
  sbg:revisionNotes: ''
- sbg:modifiedBy: david.roberson
  sbg:modifiedOn: 1621630817
  sbg:revision: 11
  sbg:revisionNotes: ''
sbg:sbgMaintained: false
sbg:validationErrors: []
