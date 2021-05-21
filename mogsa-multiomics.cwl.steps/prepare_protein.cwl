cwlVersion: v1.1
class: CommandLineTool
label: prepare-protein
$namespaces:
  sbg: https://sevenbridges.com

requirements:
- class: DockerRequirement
  dockerPull: cgc-images.sbgenomics.com/david.roberson/mogsa:1.24.0
- class: InitialWorkDirRequirement
  listing:
  - entryname: prepare-protein.R
    writable: false
    entry: |-
      #Preprocess Protein

      library(rstudioapi)
      library(HGNChelper)

      Select_unique_genes_Anova <- function(dat=dat.dup.gene.list,in.subtypes_root_t= in.subtypes_root_t) {
        
        Pvalue.Prot=NULL;
        for (k in 1:nrow(dat)){
          current.row= dat[k,]
          name.gene= current.row[,1]
          current.row.t = current.row[,-1]
          current.row.t= t(current.row.t)
          current.row.t <- data.frame(current.row.t)
          ## check before binding
          print (identical (row.names(current.row.t) , row.names(in.subtypes_root_t)))
          current.row.t <- cbind(current.row.t, in.subtypes_root_t[,1])
          result <- aov(unlist( current.row.t[,1]) ~ unlist( current.row.t[,2]))
          Pvalue= as.numeric(summary(result)[[1]][["Pr(>F)"]][1])
          Pvalue.Prot= rbind(Pvalue.Prot,data.frame(name.gene, Pvalue,current.row));
        }
        return (Pvalue.Prot);
        
      }


      ### Global protein dataset.
      in.global.protein <-  read.csv(paste0(in.dir,"/APOLLO1_Level3_Gprot_v061418_FINAL.csv"),check.names=F)
      in.global.protein <- in.global.protein [,-(1:2)]
      dim(in.global.protein)







      ## change the format of the samples name
      Sample.name = paste0( sapply(strsplit(colnames(in.global.protein), "[.]"), "[", 1), "-",
                            sapply(strsplit(colnames(in.global.protein), "[.]"), "[", 2))
      colnames(in.global.protein)
      colnames(in.global.protein)[2:101] <- Sample.name[2:101]   
      ###  only need 87 samples
      Index= which (colnames(in.global.protein) %in% row.names(in.subtypes_root_t))
      in.global.protein.samples <- in.global.protein[,Index]
      in.global.protein.samples <- in.global.protein.samples [,row.names(in.subtypes_root_t)]
      identical (row.names(in.subtypes_root_t),colnames(in.global.protein.samples))
      ### cbind back the gene name and samples.
      gb.protein <- cbind(in.global.protein$Gene,in.global.protein.samples)
      ## doublee check 
      identical (in.global.protein$`AP-JK59`,gb.protein$`AP-JK59`)
      ### 
      colnames(gb.protein) [1] <- colnames(in.global.protein)[1];
      ####### check for NA
      Index= which (is.na(gb.protein$Gene))
      ## remove NA
      in.global.protein.filter.na <-  gb.protein[-which (is.na(gb.protein$Gene)==TRUE),]
      ######## check to see there is any duplicate
      dup.genes.list=in.global.protein.filter.na[duplicated(in.global.protein.filter.na[1]) |duplicated(in.global.protein.filter.na[1], fromLast=TRUE),]
      dup.genes.list= dup.genes.list[order(dup.genes.list[1]),]
      dup.genes.list [1:4,1:3]
      dim(dup.genes.list)
      # use anova.
      ## check.
      identical (row.names(in.subtypes_root_t),colnames(dup.genes.list)[2:88])
      ## call function "Select_unique_genes_Anova"
      duplicated.genes <- Select_unique_genes_Anova(dat=dup.genes.list,in.subtypes_root_t = in.subtypes_root_t)
      ## order based on p values from smallest. P val will be column 2
      duplicated.genes=duplicated.genes[order(duplicated.genes[2]),]
      ############## remove duplicated one by keep the first smallest P val.
      uniqe.gene.list= duplicated.genes[!duplicated(duplicated.genes[1]),]
      ######## removed the P values, name of gene
      uniqe.gene.list= uniqe.gene.list[,-(1:2)]
      head(uniqe.gene.list)
      dim(uniqe.gene.list)
      ## 
      dup.rm.Protein= in.global.protein.filter.na[!(duplicated(in.global.protein.filter.na[1]) | duplicated(in.global.protein.filter.na[1], fromLast=TRUE)),]
      ## check colnames.
      identical (colnames(dup.rm.Protein),colnames(uniqe.gene.list))
      ## fix the format of sample name
      Sample.name = paste0( sapply(strsplit(colnames(uniqe.gene.list), "[.]"), "[", 1), "-",
                            sapply(strsplit(colnames(uniqe.gene.list), "[.]"), "[", 2))
      colnames(uniqe.gene.list)[2:88] <- Sample.name[2:88]   
      identical (colnames(dup.rm.Protein),colnames(uniqe.gene.list))

      ## combine if true.
      dat_Prot= rbind(uniqe.gene.list,dup.rm.Protein)
      dat_Prot[1:2,1:3]
      dim(dat_Prot)
      row.names(dat_Prot) <- dat_Prot[,1]
      dat_Prot <- dat_Prot [,-1]
      ########
      identical(row.names(in.subtypes_root_t), colnames(dat_Prot)) 
      identical(colnames(in.rna.filter), colnames(dat_Prot)) 
      dim(dat_Prot)
      ###### ## fix gene symbol if having any
      global.protein.ix.regrex <- grep(pattern = regex, row.names(dat_Prot), ignore.case = TRUE)
      row.names(dat_Prot)[global.protein.ix.regrex] 
      ID <- which(toupper(MonthMappingTable[, "Symbol"]) %in% toupper(row.names(dat_Prot)[global.protein.ix.regrex]))
      Gene_approval <- MonthMappingTable[ID,]
      idx <- c(1,3,11);
      Gene_approval <- Gene_approval[-idx,]
      row.names(Gene_approval) <- Gene_approval[,1]
      Gene_approval <- Gene_approval[toupper(row.names(dat_Prot)[global.protein.ix.regrex]) ,]
      identical (toupper(row.names(dat_Prot)[global.protein.ix.regrex]),Gene_approval[,1] )
      ## check before mapping back.
      row.names(dat_Prot)[global.protein.ix.regrex]

      Gene_approval[,2]

      ## good, mapped back based on global.protein.ix.regrex
      row.names(dat_Prot)[global.protein.ix.regrex] <- Gene_approval[,2]
      row.names(dat_Prot)[global.protein.ix.regrex]
      dim(dat_Prot)
      head(dat_Prot)
      ###### done global.
  - entryname: cwl-input.R
    writable: false

inputs:
- id: protein_abundances
  type: File?

outputs:
- id: cleaned_data
  type: File?
  outputBinding:
    glob: '*.csv'

baseCommand:
- Rscript
- prepare-protein.R
id: david.roberson/build-multiomic-analysis-mogsa/prepare-protein/2
sbg:appVersion:
- v1.1
sbg:content_hash: ace5e6313105330707d9dcdd447e7958e98f82851840d7d58970f1dcbe0751b58
sbg:contributors:
- david.roberson
sbg:createdBy: david.roberson
sbg:createdOn: 1619795092
sbg:id: david.roberson/build-multiomic-analysis-mogsa/prepare-protein/2
sbg:image_url:
sbg:latestRevision: 2
sbg:modifiedBy: david.roberson
sbg:modifiedOn: 1620238957
sbg:project: david.roberson/build-multiomic-analysis-mogsa
sbg:projectName: 'BUILD: Multiomic Analysis MOGSA'
sbg:publisher: sbg
sbg:revision: 2
sbg:revisionNotes: ''
sbg:revisionsInfo:
- sbg:modifiedBy: david.roberson
  sbg:modifiedOn: 1619795092
  sbg:revision: 0
  sbg:revisionNotes:
- sbg:modifiedBy: david.roberson
  sbg:modifiedOn: 1619795158
  sbg:revision: 1
  sbg:revisionNotes: ''
- sbg:modifiedBy: david.roberson
  sbg:modifiedOn: 1620238957
  sbg:revision: 2
  sbg:revisionNotes: ''
sbg:sbgMaintained: false
sbg:validationErrors: []
