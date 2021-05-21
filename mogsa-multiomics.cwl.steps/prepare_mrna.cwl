cwlVersion: v1.2
class: CommandLineTool
label: prepare-mrna
$namespaces:
  sbg: https://sevenbridges.com

requirements:
- class: DockerRequirement
  dockerPull: cgc-images.sbgenomics.com/david.roberson/mogsa:210521
- class: InitialWorkDirRequirement
  listing:
  - entryname: prepare_mRNA.R
    writable: false
    entry: |
      require(tidyverse)
      library(EnsDb.Hsapiens.v86)


      in_rna = fs::dir_info("./", glob = "*txt") %>% 
        dplyr::select(path) %>% 
        mutate(sample_id = stringr::str_split(.$path, "/")) %>% 
        mutate(sample_id = map_chr(.$sample_id, ~tail(., 1))) %>% 
        mutate(sample_id = stringr::str_split(.$sample_id, ".txt")) %>% 
        mutate(sample_id = map_chr(.$sample_id, ~head(., 1))) %>% 
        group_by(path, sample_id) %>% 
        glimpse() %>% 
        do({
          temp = .
          temp = readr::read_delim(temp$path, delim = "\t", col_names = c("ensembl_gene", "counts"))
        }) %>% 
        ungroup() %>% 
        dplyr::select(ensembl_gene, sample_id, counts) %>% 
        pivot_wider(ensembl_gene, names_from = sample_id, values_from = counts) %>% 
        mutate(gene_symbol = unlist(lapply(.$ensembl_gene, function(x) strsplit(x, '.', fixed = TRUE)[[1]][1]))) %>% 
        mutate(gene_symbol = mapIds(EnsDb.Hsapiens.v86,
                                    keys = .$gene_symbol,
                                    keytype = "GENEID",
                                    column = "SYMBOL")) %>% 
        dplyr::filter(!(is.na(gene_symbol))) %>% 
        dplyr::select(gene_symbol, everything(), -ensembl_gene) %>% 
        write_csv("mRNA_table.csv")
  - $(inputs.mRNA_matrix)
- class: InlineJavascriptRequirement

inputs:
- id: mRNA_matrix
  type: File[]?

outputs:
- id: mRNA_table
  type: File?
  outputBinding:
    glob: '*.csv'

baseCommand:
- Rscript
- prepare_mRNA.R
id: david.roberson/build-multiomic-analysis-mogsa/prepare-mrna/10
sbg:appVersion:
- v1.2
sbg:content_hash: a4dd2cc629d02f62071af53bd68adf419317ac2a8b950ed97b2f835d37f43b9a0
sbg:contributors:
- david.roberson
sbg:createdBy: david.roberson
sbg:createdOn: 1619795187
sbg:id: david.roberson/build-multiomic-analysis-mogsa/prepare-mrna/10
sbg:image_url:
sbg:latestRevision: 10
sbg:modifiedBy: david.roberson
sbg:modifiedOn: 1621630582
sbg:project: david.roberson/build-multiomic-analysis-mogsa
sbg:projectName: 'BUILD: Multiomic Analysis MOGSA'
sbg:publisher: sbg
sbg:revision: 10
sbg:revisionNotes: head
sbg:revisionsInfo:
- sbg:modifiedBy: david.roberson
  sbg:modifiedOn: 1619795187
  sbg:revision: 0
  sbg:revisionNotes:
- sbg:modifiedBy: david.roberson
  sbg:modifiedOn: 1619795218
  sbg:revision: 1
  sbg:revisionNotes: ''
- sbg:modifiedBy: david.roberson
  sbg:modifiedOn: 1620238730
  sbg:revision: 2
  sbg:revisionNotes: added docker registry
- sbg:modifiedBy: david.roberson
  sbg:modifiedOn: 1620238845
  sbg:revision: 3
  sbg:revisionNotes: added in R script
- sbg:modifiedBy: david.roberson
  sbg:modifiedOn: 1621607339
  sbg:revision: 4
  sbg:revisionNotes: ''
- sbg:modifiedBy: david.roberson
  sbg:modifiedOn: 1621607426
  sbg:revision: 5
  sbg:revisionNotes: ''
- sbg:modifiedBy: david.roberson
  sbg:modifiedOn: 1621619654
  sbg:revision: 6
  sbg:revisionNotes: ''
- sbg:modifiedBy: david.roberson
  sbg:modifiedOn: 1621620451
  sbg:revision: 7
  sbg:revisionNotes: ''
- sbg:modifiedBy: david.roberson
  sbg:modifiedOn: 1621627272
  sbg:revision: 8
  sbg:revisionNotes: cgc-images.sbgenomics.com/david.roberson/mogsa:210521
- sbg:modifiedBy: david.roberson
  sbg:modifiedOn: 1621630177
  sbg:revision: 9
  sbg:revisionNotes: ''
- sbg:modifiedBy: david.roberson
  sbg:modifiedOn: 1621630582
  sbg:revision: 10
  sbg:revisionNotes: head
sbg:sbgMaintained: false
sbg:validationErrors: []
