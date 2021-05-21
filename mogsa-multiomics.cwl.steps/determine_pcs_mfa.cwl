cwlVersion: v1.1
class: CommandLineTool
label: determine-pcs-mfa
$namespaces:
  sbg: https://sevenbridges.com

requirements:
- class: LoadListingRequirement
- class: DockerRequirement
  dockerPull: cgc-images.sbgenomics.com/david.roberson/mogsa:1.24.0
- class: InitialWorkDirRequirement
  listing:
  - entryname: pcs_based_on_mfa.R
    writable: false
    entry: |-
      ### choosing # of PCs based on MFA
      library(rstudioapi)
      library(knitr)
      library(mogsa)
      library(graphite)
      library(VennDiagram);
      library(grid);
      library(ggplot2)
      library(venn)

      source("cwl_input.R")

      dat<- readRDS(file = RNA_ProtGlobal_Non_Phospho_rds_path)



      row.names(dat$RNA) <- toupper(row.names(dat$RNA) )
      dim(dat$RNA)
      row.names(dat$G.Pro) <- toupper(row.names(dat$G.Pro) )
      dim(dat$G.Pro)
      row.names(dat$Non_G_MS_RPPA) <- toupper(row.names(dat$Non_G_MS_RPPA) )
      dim(dat$Non_G_MS_RPPA)
      ## create a directory to store output.

      out.dir <- "MFA_Evaluated_Step2";
      dir.create(out.dir);
      ## first observe distribution of 3 datasets.

      pdf(paste0( out.dir, "/filter.normalized.RNA.hist.pdf"));
      hist(as.numeric(unlist(dat$RNA)),main="RNA", xlab="Log2(TPM)",freq=FALSE,breaks=100)
      dev.off();
      pdf(paste0( out.dir, "/filter.normalized.Global.P.hist.pdf"));
      hist(as.numeric(unlist(dat$G.Pro)), freq=FALSE,main="Global Protein", xlab="Normalized Global Protein",breaks=100)
      dev.off();
      pdf(paste0( out.dir, "/filter.non.Glb_normalized.Phospho.P.hist.pdf"));
      hist(as.numeric(unlist(dat$Non_G_MS_RPPA)), freq=FALSE,main="MS_RPPA", xlab="PHOSPHO_MS_RPPA_combined_Non_global_norm_KNNImpute_K_10_measFrac_0.5_skip_RPPA_multimap",breaks=100)
      dev.off();
      ###### get the gene list.
      gene.list <- lapply(dat, rownames)
      ## Multiple omics data analysis using MFA 
      moa.mfa <- moa(dat, proc.row = "center_ssq1", w.data = "lambda1", statis = FALSE)
      ### plot the corelation of datasets
      pdf(paste0(out.dir,"/mfa.non.Gl.norm.rv.pdf"));
      plot(moa.mfa, value = "RV")
      dev.off();
      ## to help to choose the # of PCs based on bootstrap.
      pdf(paste0(out.dir,"/mfa.method.bootmoaPlot.eig.pdf"),18,10);
      bt <- bootMoa(moa = moa.mfa, proc.row = "center_ssq1", w.data = "lambda1", statis = FALSE, B = 20,plot=TRUE)
      dev.off();
      ## scree plot -- using epigen value
      pdf(paste0(out.dir,"/mfa.method.ScreePlot.eig.pdf"),18,10);

      plot(moa.mfa , value="eig", type=2, n=20);
      dev.off();
      ## the same with eig, but in the eigenvalues are scaled to 1
      pdf(paste0(out.dir,"/mfa.method.ScreePlot.tau.pdf"),18,10);
      plot(moa.mfa , value="tau", type=2, n=20)
      dev.off()
      ######## percentage of variance.
      pdf(paste0(out.dir,"/mfa.method.ScreePlot.eig_percentage.pdf"),18,10);
      moana.mfa.percentage= moa.mfa ;
      moana.mfa.percentage@partial.eig = 100 * (moana.mfa.percentage@partial.eig/sum(moana.mfa.percentage@partial.eig))
      plot(moana.mfa.percentage, value="eig", type=2, n=20);
      dev.off();
      ## to help to choose the # of PCs based on bootstrap.
      pdf(paste0(out.dir,"/mfa.method.bootmoaPlot.percentage.pdf"),18,10);
      bt <- bootMoa(moa = moa.mfa, proc.row = "center_ssq1", w.data = "lambda1", statis = FALSE, B = 20,plot=TRUE)
      dev.off();
      #####
      ######
      ## this one is used to calculate the percentatge of nf=3,4,and 5,6.. 
      ###### 
      from.i  <- 1;
      type.n  <- c(3, 4, 5,6);
      out.s <- NULL;
      for (k in 1:4){
        to.i   <- type.n[k];
        percen.total= sum(moana.mfa.percentage@partial.eig[,from.i:to.i])
        RNA.per=sum(moana.mfa.percentage@partial.eig[1,from.i:to.i])/sum(moana.mfa.percentage@partial.eig[,from.i:to.i]);
        Global.per=sum(moana.mfa.percentage@partial.eig[2,from.i:to.i])/sum(moana.mfa.percentage@partial.eig[,from.i:to.i]);
        phospho.per=sum(moana.mfa.percentage@partial.eig[3,from.i:to.i])/sum(moana.mfa.percentage@partial.eig[,from.i:to.i]);
        
        out.s  <- rbind(out.s, cbind(percen.total,RNA.per,Global.per,phospho.per))
      }
      row.names (out.s) <- c("nf.3","nf.4","nf.5", "nf.6");
      write.csv (out.s, paste0(out.dir, "/nf3.4.5.6.percentatge.csv"))

      ############# loop to get statistics for comparison of nf=3,4,5,6
      ## change to get to the docker.
      path<-"/sbgenomics/project-files/MOGSA_Input/database/"
      ## if local R
      file.names <- dir(path, pattern <- ".gmt")

      for ( i in 1: length(file.names)){
        predatabase<- prepMsigDB(file=paste0(path,file.names[i]));
        length(predatabase)
        sup.dat <- prepSupMoa(gene.list, geneSets=predatabase)
        mogsa.dat.nf.3 <- mogsa(x = dat, sup=sup.dat , nf=3,
                                proc.row = "center_ssq1", w.data = "lambda1", statis = FALSE)
        mogsa.dat.nf.4 <- mogsa(x = dat, sup=sup.dat , nf=4,
                                proc.row = "center_ssq1", w.data = "lambda1", statis = FALSE)
        mogsa.dat.nf.5 <- mogsa(x = dat, sup=sup.dat , nf=5,
                                proc.row = "center_ssq1", w.data = "lambda1", statis = FALSE)
        mogsa.dat.nf.6 <- mogsa(x = dat, sup=sup.dat , nf=6,
                                proc.row = "center_ssq1", w.data = "lambda1", statis = FALSE)
        ##### get P values for the nf=5,nf=10,nf=3
        p.mat.3 <- getmgsa(mogsa.dat.nf.3, "p.val") 
        p.mat.4 <- getmgsa(mogsa.dat.nf.4, "p.val") 
        p.mat.5 <- getmgsa(mogsa.dat.nf.5, "p.val") 
        p.mat.6 <- getmgsa(mogsa.dat.nf.6, "p.val") 
        #### # of pts have p < 0.01
        ######### 
        if (i==6){
          top.gs.3 <- sort(rowSums(p.mat.3 < 0.01), decreasing = TRUE)[1:20]
          top.gs.4 <- sort(rowSums(p.mat.4 < 0.01), decreasing = TRUE)[1:20]
          top.gs.5 <- sort(rowSums(p.mat.5 < 0.01), decreasing = TRUE)[1:20]
          top.gs.6 <- sort(rowSums(p.mat.6 < 0.01), decreasing = TRUE)[1:20]
        }else if (i<3){
          top.gs.3 <- sort(rowSums(p.mat.3 < 0.01), decreasing = TRUE)[1:30]
          top.gs.4 <- sort(rowSums(p.mat.4 < 0.01), decreasing = TRUE)[1:30]
          top.gs.5 <- sort(rowSums(p.mat.5 < 0.01), decreasing = TRUE)[1:30]
          top.gs.6 <- sort(rowSums(p.mat.6 < 0.01), decreasing = TRUE)[1:30]
        }else {
          top.gs.3 <- sort(rowSums(p.mat.3 < 0.01), decreasing = TRUE)[1:50]
          top.gs.4 <- sort(rowSums(p.mat.4 < 0.01), decreasing = TRUE)[1:50]
          top.gs.5 <- sort(rowSums(p.mat.5 < 0.01), decreasing = TRUE)[1:50]
          top.gs.6 <- sort(rowSums(p.mat.6 < 0.01), decreasing = TRUE)[1:50]
        }
        ##  this one for nf=5 top.gs= top.gs, P.mat=p.mat top.gs= top.gs,P.mat=p.mat
        
        top.gs.name.3 <- names(top.gs.3)
        top.gs.name.3
        p.mat.3.filter= p.mat.3[top.gs.name.3,]
        ##  this one for nf=4
        top.gs.name.4 <- names(top.gs.4)
        top.gs.name.4
        p.mat.4.filter= p.mat.4[top.gs.name.4,]
        #############
        ##  this one for nf=5
        top.gs.name.5 <- names(top.gs.5)
        top.gs.name.5
        p.mat.5.filter= p.mat.5[top.gs.name.5,]
        ########
        ### nf=6
        top.gs.name.6 <- names(top.gs.6)
        top.gs.name.6
        p.mat.6.filter= p.mat.6[top.gs.name.6,]
        overlap <- calculate.overlap(
          x = list(
            nf.3 = top.gs.name.3,
            nf.4 = top.gs.name.4,
            nf.5 = top.gs.name.5,
            nf.6=top.gs.name.6)
        );
        overlap$a6
        ##### nf=5, and nf=6,nf=8 filter based on the overlap
        index= which (row.names(p.mat.3.filter) %in% overlap$a6)
        p.mat.3.filter=p.mat.3.filter[index,]
        index= which (row.names(p.mat.4.filter) %in% overlap$a6)
        p.mat.4.filter=p.mat.4.filter[index,]
        index= which (row.names(p.mat.5.filter) %in% overlap$a6)
        p.mat.5.filter=p.mat.5.filter[index,]
        index= which (row.names(p.mat.6.filter) %in% overlap$a6)
        p.mat.6.filter=p.mat.6.filter[index,]
        ### arrange the same order.
        p.mat.6.filter <- p.mat.6.filter[row.names(p.mat.4.filter),]
        p.mat.3.filter <- p.mat.3.filter[row.names(p.mat.4.filter),]
        p.mat.5.filter <- p.mat.5.filter[row.names(p.mat.4.filter),]
        identical(row.names(p.mat.4.filter),row.names(p.mat.3.filter))
        identical(row.names(p.mat.4.filter),row.names(p.mat.6.filter))
        identical(row.names(p.mat.5.filter),row.names(p.mat.6.filter))
        
        ########### get the  overlap for comparison
        p.table.nf.3=NULL
        p.table.nf.4=NULL
        p.table.nf.5=NULL
        p.table.nf.6=NULL
        for (pv in 1:nrow(p.mat.4.filter)){
          p.table.nf.3 <- rbind(p.table.nf.3, length(which(p.mat.3.filter[pv,] <0.01 )))
          p.table.nf.4 <- rbind(p.table.nf.4, length(which(p.mat.4.filter[pv,] <0.01 )))
          p.table.nf.5 <- rbind(p.table.nf.5, length(which(p.mat.5.filter[pv,] <0.01 )))
          p.table.nf.6 <- rbind(p.table.nf.6, length(which(p.mat.6.filter[pv,] <0.01 )))
        }
        p.table=NULL;
        p.table=cbind(p.table.nf.3,p.table.nf.4,p.table.nf.5,p.table.nf.6)
        row.names(p.table) <- row.names(p.mat.4.filter)
        colnames(p.table) <- c("nf3_pval","nf4_pval", "nf5_pval","nf6_pval")
        p.table <- data.frame(p.table)
        ######
        x <- 1:nrow(p.mat.4.filter);
        y1=data.frame(p.table[,1])
        y2= data.frame(p.table[,2])
        y3= data.frame(p.table[,3])
        y4= data.frame(p.table[,4])
        dir.create(paste0(out.dir,"/Result_",file.names[i]))
        pdf(paste0(out.dir,"/Result_",file.names[i],"/compared_nf3_nf4_nf5_nf6.plot.pdf"), 15,10);
        #plot(xdata, y1, type="o", col="blue", pch="o", lty=1, ylim=c(0,110), ylab="y" )
        plot(x, y1$p.table...1., type="o", pch="o", col="blue", ylim=c(0,110), xlab=paste0(length(overlap$a6)," Pathway overlapped between nf.3_nf.4_nf.5_nf.6"), ylab="# P val<0.01 in the pathway across samples")
        points(x, y2$p.table...2., col="red", pch="*")
        lines(x, y2$p.table...2., col="red",lty=2)
        points(x, y3$p.table...3., col="black",pch="+")
        lines(x, y3$p.table...3., col="black", lty=3)
        points(x, y4$p.table...4., col="purple",pch="+")
        lines(x, y4$p.table...4., col="purple", lty=3)
        legend("topright", legend=c("nf3","nf4","nf5", "nf6"),
               col=c("blue", "red","black","purple"), pch=c("o","*","+"),lty=c(1,2,3,4), ncol=1,cex=1.5)
        dev.off();
        ####
        pdf(paste0(out.dir,"/Result_",file.names[i],"/overlap_plot.pdf"));
        my_list <- list(nf.3 = top.gs.name.3,
                        nf.4 = top.gs.name.4,
                        nf.5 = top.gs.name.5,
                        nf.6=top.gs.name.6)
        venn(my_list)
        dev.off()
      }

      ############### review all result in the folder MFA_evaluated_result to choose the best # of PCs.
  - entryname: cwl_input.R
    writable: false
    entry: |2-

      RNA_ProtGlobal_Non_Phospho_rds_path = "$(inputs.RNA_ProtGlobal_Non_Phospho.path)"
- class: InlineJavascriptRequirement

inputs:
- id: RNA_ProtGlobal_Non_Phospho
  type: File
- id: input_files
  type: File[]
- id: mogsa_database
  type: File

outputs:
- id: output
  type: Directory?
  outputBinding:
    glob: MFA_Evaluated_Step2
    loadListing: no_listing

baseCommand:
- Rscript
- pcs_based_on_mfa.R
id: david.roberson/build-multiomic-analysis-mogsa/determine-pcs-mfa/3
sbg:appVersion:
- v1.1
sbg:content_hash: a141aa7ac7b5170ca720d81d2c43ccdb928980773b0b432a6d5d789fd31b55429
sbg:contributors:
- david.roberson
sbg:createdBy: david.roberson
sbg:createdOn: 1619795304
sbg:id: david.roberson/build-multiomic-analysis-mogsa/determine-pcs-mfa/3
sbg:image_url:
sbg:latestRevision: 3
sbg:modifiedBy: david.roberson
sbg:modifiedOn: 1619796359
sbg:project: david.roberson/build-multiomic-analysis-mogsa
sbg:projectName: 'BUILD: Multiomic Analysis MOGSA'
sbg:publisher: sbg
sbg:revision: 3
sbg:revisionNotes: ''
sbg:revisionsInfo:
- sbg:modifiedBy: david.roberson
  sbg:modifiedOn: 1619795304
  sbg:revision: 0
  sbg:revisionNotes:
- sbg:modifiedBy: david.roberson
  sbg:modifiedOn: 1619795446
  sbg:revision: 1
  sbg:revisionNotes: ''
- sbg:modifiedBy: david.roberson
  sbg:modifiedOn: 1619795742
  sbg:revision: 2
  sbg:revisionNotes: ''
- sbg:modifiedBy: david.roberson
  sbg:modifiedOn: 1619796359
  sbg:revision: 3
  sbg:revisionNotes: ''
sbg:sbgMaintained: false
sbg:validationErrors: []
