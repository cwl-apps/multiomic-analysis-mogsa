cwlVersion: v1.1
class: CommandLineTool
label: MOGSA
$namespaces:
  sbg: https://sevenbridges.com

requirements:
- class: ShellCommandRequirement
- class: LoadListingRequirement
- class: DockerRequirement
  dockerPull: cgc-images.sbgenomics.com/david.roberson/mogsa:1.24.0
- class: InitialWorkDirRequirement
  listing:
  - entryname: Dockerfile
    writable: false
    entry: |-
      FROM bioconductor/bioconductor_docker:RELEASE_3_12
      
      RUN Rscript -e 'BiocManager::install("mogsa")'
      RUN Rscript -e 'install.packages("rmarkdown")'
      RUN Rscript -e 'BiocManager::install("BiocStyle")'
      RUN Rscript -e 'install.packages("magick")'
      RUN Rscript -e 'BiocManager::install("org.Hs.eg.db")'
  - entryname: run_mogsa.R
    writable: false
    entry: |-
      #database_path
      #RNA_ProtGlobal_Non_Phospho
      #subtypes_csv_path
      
      library(rstudioapi)
      library(knitr)
      library(mogsa)
      library(gplots)
      
      source("cwl_input.R")
      source("functions.R")
      
      ################################# do pathway analysis based on the first 3 PCs chosen from previous step.
      
      mutliomics_data = readRDS(multiomics_data_path)
      
      subtypes = read_csv(subtypes_csv_path)
      #/sbgenomics/project-files/MOGSA_Input/apollo1.master.table.prelim.20190606.csv"
      
      
      ###### check that subtypes info matches the dat object
      identical (row.names(in.subtypes_root_t), colnames(dat$RNA));
      identical (row.names(in.subtypes_root_t), colnames(dat$Non_G_MS_RPPA));
      identical (row.names(in.subtypes_root_t), colnames(dat$G.Pro));
      row.names(dat$RNA) <- toupper(row.names(dat$RNA) )
      row.names(dat$G.Pro) <- toupper(row.names(dat$G.Pro) )
      row.names(dat$Non_G_MS_RPPA) <- toupper(row.names(dat$Non_G_MS_RPPA) )
      
      PI.i  <- which(in.subtypes_root_t$subtypePred=="ProxInflam");
      PP.i <- which(in.subtypes_root_t$subtypePred=="ProxProlif");
      TRU.i <- which(in.subtypes_root_t$subtypePred=="TRU");
      len.PI <- length(PI.i);
      len.PP <- length(PP.i);
      len.TRU <- length(TRU.i);
      cancersubType <- c(rep("TRU",len.TRU), rep("ProxProlif",len.PP),rep("ProxInflam",len.PI));
      cancersubType <- factor(cancersubType, levels= c("TRU","ProxProlif","ProxInflam"));
      colcode <- cancersubType
      levels(colcode) <- c("black", "red", "green")
      colcode <- as.character(colcode,stringsAsFactors=FALSE)
      my.palette <- colorpanel(1000,"blue","white","red")
      gene.list <- lapply(dat, rownames)
      
      #pre-made database for pathway database...
      database_path<-"/sbgenomics/project-files/MOGSA_Input/database/";
      file.names <- dir(database_path, pattern <-".gmt")
      file.names <- file.names[-(4:5)]
      
      for ( i in 1: length(file.names)){
        out.dir.sub <- paste0(out.dir,"/Result_",file.names[i]);
        dir.create(out.dir.sub);
        predatabase<- prepMsigDB(file=paste0(path,file.names[i]));
        if (i==3){
          sup.dat <- prepSupMoa(gene.list, geneSets=predatabase, minMatch = 5, maxMatch = 500 );
        }else { sup.dat <- prepSupMoa(gene.list, geneSets=predatabase);}
        mogsa.dat <- mogsa(x = dat, sup=sup.dat , nf=3,
                           proc.row = "center_ssq1", w.data = "lambda1", statis = FALSE)
        ## get the p values for pathway.. but not for p.val corrected.
        p.mat <- getmgsa(mogsa.dat, "p.val");
        ## this is used to check pval and their rows and columns ( logically.) .. rows and column names shoudl be the same in the mATRIX.
        p.val.checked <- mogsa.dat@sup@p.val
        identical (p.val.checked, p.mat)
        ##### same matrix with dimension. 
        p.val.corrected <- mogsa.dat@sup@p.val.corrected
        colnames(p.val.corrected) <- colnames(p.mat)
        row.names(p.val.corrected) <- row.names(p.mat)
        #### sort based on P val. corrected.
        top.gs.based.on.val.corrected <- sort(rowSums(p.val.corrected < 0.01), decreasing = TRUE)
        ## get the name of pathway.
        top.gs.name.corrected <- names(top.gs.based.on.val.corrected)
        ## mapped back in the order.
        p.val.corrected <- p.val.corrected[top.gs.name.corrected, ]
        ## 
        table.count.p <- NULL;
        for (j in 1:nrow(p.val.corrected) ){
          table.count.p <- rbind(table.count.p,length (which (p.val.corrected[j,] <0.01)))
        }
        row.names(table.count.p) <- row.names(p.val.corrected)
        colnames(table.count.p) <- "#ofPTs.having.Pval.corrected.01"
        ## cut off 45/87 percent 50% having P values smaller 0.01
        idx.p.val.45 <- which (table.count.p[,1] >45)
        idx.p.val.45.name <- names(idx.p.val.45)
        scores.mogsa<- getmgsa(mogsa.dat, "score")
        ###### ##  filter based on the significant pathway.
        dat.scores <- scores.mogsa[idx.p.val.45.name, ]
        identical(colnames(dat.scores), row.names(in.subtypes_root_t))
        #### do three glm to obtain the best representative for each subtype.
        TRU.cc.i <- c(rep(1,len.TRU), rep(0,len.PP),rep(0,len.PI));
        
        outfile.name.TRU=paste0(out.dir.sub ,"/glm_out_TRU.txt")
        result.TRU <- glm.function (dat.scores,subtypes.compared=TRU.i, others.subtypes=c(PP.i, PI.i), cc.i= TRU.cc.i,outfile.name=outfile.name.TRU); 
        result.TRU
        ## PP 
        PP.cc.i <- c(rep(0,len.TRU), rep(1,len.PP),rep(0,len.PI));
        outfile.name.PP=paste0(out.dir.sub ,"/glm_out_PP.txt")
        result.PP <- glm.function (dat.scores,subtypes.compared=PP.i, others.subtypes=c(TRU.i,PI.i), cc.i= PP.cc.i,outfile.name=outfile.name.PP); 
        result.PP
        ## PI
        PI.cc.i <- c(rep(0,len.TRU), rep(0,len.PP),rep(1,len.PI));
        outfile.name.PI=paste0(out.dir.sub ,"/glm_out_PI.txt")
        result.PI <- glm.function (dat.scores,subtypes.compared=PI.i, others.subtypes=c(TRU.i,PP.i), cc.i= PI.cc.i,outfile.name=outfile.name.PI); 
        result.PI
        #### 
        #######
        pathways.PI.PP.TRU <- rbind( result.PI, result.PP,  result.TRU)
        ### remove duplicate .
        pathways.PI.PP.TRU.unique <- pathways.PI.PP.TRU[!duplicated( pathways.PI.PP.TRU[1]),];
        #### mapped back to dat.scores.
        dat.scores_filter <- dat.scores[pathways.PI.PP.TRU.unique$Pathway_name,];
        ####### 
        ######
        identical (colnames(dat.scores_filter),row.names(in.subtypes_root_t) );
        pdf(paste0(out.dir.sub ,"/glm.heatmap.pathways.pdf"),15,10); 
        heatmap.2(dat.scores_filter,trace = "n", scale = "row", margins = c(7, 20),ColSideColors=colcode,
                  col=my.palette,offsetRow=-.3,cexCol=0.6,cexRow=0.8,lhei=c(1.5,10), Colv=FALSE)
        legend(0.85,1.03, legend=c("TRU", "PP", "PI"), col= c("#1A1A1A","#E41A1C","#51B25D"),
               cex=0.8, lty=1, lwd=4, xpd=TRUE);
        dev.off();
        
       
        for ( p in 1:length(row.names(dat.scores_filter))){
          ## create directory
          dir.create(paste0(out.dir.sub ,"/",row.names(dat.scores_filter)[p]));
          pdf(paste0(out.dir.sub,"/",row.names(dat.scores_filter)[p],"/",row.names(dat.scores_filter)[p],"GIS.pdf"));
          dat.pathway.GIS= GIS(mogsa.dat, row.names(dat.scores_filter)[p],Fvalue = TRUE, ff = cancersubType,nf = 3);
          dev.off();
          pdf(paste0(out.dir.sub,"/",row.names(dat.scores_filter)[p],"/",row.names(dat.scores_filter)[p],"datawise_decomposition.pdf"))
          ###
          decompose.gs.group(mogsa.dat,row.names(dat.scores_filter)[p],cancersubType, decomp = "data", nf = 3, x.legend = "bottomright", y.legend = NULL, plot = TRUE, main = row.names(dat.scores_filter)[p] )
          dev.off();
        }
      }
      
      ########## end.
  - entryname: functions.R
    writable: false
    entry: |-
      ### functions to be used in the next part.
      ### this function helps to select representative based on each RNA subtype
      
      glm.function<- function(dat.scores,subtypes.compared, others.subtypes, cc.i,outfile.name) {
        PCA.adj <- F; 
        out.s   <- c();
        for (k in 1:nrow(dat.scores)){
          g <- row.names(dat.scores)[k];
          out   <- t.test(dat.scores[k,subtypes.compared], dat.scores[k,others.subtypes]);  
          if (PCA.adj){
            fit   <- summary(glm(cc.i ~ dat.scores[k,] + pca$x[,1] + pca$x[,2], family = binomial));
          } else {
            fit   <- summary(glm(cc.i ~ dat.scores[k,] , family = binomial));    
            a<- try(fit$coef[2,])    
            if (inherits(a, "try-error")) next    
          }
          out.s <- rbind(out.s, c(g, fit$coef[2,3], fit$coef[2,4], out$conf.int[1], out$conf.int[2],out$p.value));
        }
        colnames(out.s) <- c("Pathway_name","glm.t.value","glm.p.value", "T.test.95%CInt.low", "T.test.95%Cint.high","t.test.pval");
        ### choose tthe method: bon, and fdr
        p.adjust.M <- p.adjust.methods[c(4,7)];
        p.adj   <- sapply(p.adjust.M, function(meth) p.adjust(out.s[,3], meth));
        colnames(p.adj) <- c("glm.bon","glm.fdr")
        out.p   <- cbind(out.s, p.adj);
        out.p   <- out.p[order(as.numeric(out.p[,3])),];
        write.table(out.p, outfile.name, sep="\t", row.names=F, quote=F);
        result <- filter.function (filename= outfile.name,out.file.name= paste0(outfile.name,"_filter.top5.lowest5.tval.txt"));
        return (result)
      }
      
      filter.function <- function ( filename,out.file.name) {
        in.file.name <- read.table (filename, stringsAsFactors = F,header=T)
        idx <- which(in.file.name[,3]< 0.05);
        in.file.name.filter <- in.file.name[idx,];
        ## sort based on glm.t.value
        in.file.name.filter <- in.file.name.filter[order(as.numeric(in.file.name.filter[,2]),decreasing=T),]
        dim(in.file.name.filter)
        ## get the last 5 with p val 
        in.file.name.filter <-rbind(head(in.file.name.filter,n=5),tail(in.file.name.filter,n=5));
        write.table(in.file.name.filter, out.file.name, sep="\t", row.names=F, quote=F);
        ####### to check T.test significant.
        return (in.file.name.filter);
      }
- class: InlineJavascriptRequirement

inputs:
- id: mfa_dir
  type: Directory?
  inputBinding:
    position: 0
    shellQuote: false
  loadListing: deep_listing
- id: mogsa_database
  type: File?
  inputBinding:
    position: 0
    shellQuote: false

outputs:
- id: R_workspace
  type: File?
  outputBinding:
    glob: '*.RData'
- id: html
  type: File?
  outputBinding:
    glob: '*.html'
stdout: mogsa.log

baseCommand:
- Rscript
- render_mogsa.R

hints:
- class: sbg:SaveLogs
  value: '*.R'
- class: sbg:SaveLogs
  value: '*.log'
id: david.roberson/build-multiomic-analysis-mogsa/mogsa/2
sbg:appVersion:
- v1.1
sbg:content_hash: aae4a15ae5cce733ecdd6b8a8d6b5066c42870b42543806ac9af844168fb789b4
sbg:contributors:
- david.roberson
sbg:createdBy: david.roberson
sbg:createdOn: 1619459412
sbg:id: david.roberson/build-multiomic-analysis-mogsa/mogsa/2
sbg:image_url:
sbg:latestRevision: 2
sbg:modifiedBy: david.roberson
sbg:modifiedOn: 1620246596
sbg:project: david.roberson/build-multiomic-analysis-mogsa
sbg:projectName: 'BUILD: Multiomic Analysis MOGSA'
sbg:publisher: sbg
sbg:revision: 2
sbg:revisionNotes: ''
sbg:revisionsInfo:
- sbg:modifiedBy: david.roberson
  sbg:modifiedOn: 1619459412
  sbg:revision: 0
  sbg:revisionNotes: Copy of david.roberson/mogsa-workflow-wrapping/mogsa/15
- sbg:modifiedBy: david.roberson
  sbg:modifiedOn: 1619795505
  sbg:revision: 1
  sbg:revisionNotes: ''
- sbg:modifiedBy: david.roberson
  sbg:modifiedOn: 1620246596
  sbg:revision: 2
  sbg:revisionNotes: ''
sbg:sbgMaintained: false
sbg:validationErrors: []
