cwlVersion: v1.1
class: CommandLineTool
label: MOGSA
$namespaces:
  sbg: https://sevenbridges.com

requirements:
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
  - entryname: mogsa_gene_set_analysis.R
    writable: false
    entry: |
      ## ----include=FALSE---------------------------------------------------------
      library(knitr)
      opts_chunk$set(
      concordance=TRUE
      )
      
      ## ----include=FALSE---------------------------------------------------------
      library(knitr)
      opts_chunk$set(
      concordance=TRUE, cache=TRUE, message=FALSE, out.width=".55\\textwidth", echo=TRUE, fig.width=6, fig.height=6, fig.align="center", result="markup", hold=TRUE
      )
      
      ## ----loadLib---------------------------------------------------------------
      # loading gene expression data and supplementary data
      library(mogsa)
      library(gplots) # used for visulizing heatmap
      # loading gene expression data and supplementary data
      data(NCI60_4array_supdata)
      data(NCI60_4arrays)
      
      ## ----dataDim---------------------------------------------------------------
      sapply(NCI60_4arrays, dim) # check dimensions of expression data
      sapply(NCI60_4array_supdata, dim) # check dimensions of supplementary data
      # check if the gene expression data and annotation data are mathced in the same order
      identical(names(NCI60_4arrays), names(NCI60_4array_supdata)) 
      head(rownames(NCI60_4arrays$agilent)) # the type of gene IDs
      
      ## ----dataColMatch----------------------------------------------------------
      dataColNames <- lapply(NCI60_4arrays, colnames)
      supColNames <- lapply(NCI60_4arrays, colnames)
      identical(dataColNames, supColNames)
      
      ## ----defineCancerType------------------------------------------------------
      # define cancer type
      cancerType <- as.factor(substr(colnames(NCI60_4arrays$agilent), 1, 2))
      # define color code to distinguish cancer types
      colcode <- cancerType
      levels(colcode) <- c("black", "red", "green", "blue", 
                           "cyan", "brown", "pink", "gray", "orange")
      colcode <- as.character(colcode)
      
      ## ----mogsaBasicRun---------------------------------------------------------
      mgsa1 <- mogsa(x = NCI60_4arrays, sup=NCI60_4array_supdata, nf=3,
                     proc.row = "center_ssq1", w.data = "inertia", statis = TRUE)
      
      ## ----eigenPlot, fig.cap="The variance of each principal components (PC), the contributions of different data are distinguished by different colors", fig.width=4, fig.height=4----
      eigs <- getmgsa(mgsa1, "partial.eig") # get partial "eigenvalue" for separate data 
      barplot(as.matrix(eigs), legend.text = rownames(eigs))
      
      ## ----scoreMatrix, fig.cap="heatmap showing the gene set score (GSS) matrix"----
      # get the score matrix
      scores <- getmgsa(mgsa1, "score")
      heatmap.2(scores, trace = "n", scale = "r", Colv = NULL, dendrogram = "row", 
                margins = c(6, 10), ColSideColors=colcode)
      
      ## ----subsetScoreMatrix, fig.cap="heatmap showing the gene set score (GSS) matrix for top 20 significant gene sets"----
      p.mat <- getmgsa(mgsa1, "p.val") # get p value matrix
      # select gene sets with most signficant GSS scores.
      top.gs <- sort(rowSums(p.mat < 0.01), decreasing = TRUE)[1:20]
      top.gs.name <- names(top.gs)
      top.gs.name
      heatmap.2(scores[top.gs.name, ], trace = "n", scale = "r", Colv = NULL, dendrogram = "row",
                margins = c(6, 10), ColSideColors=colcode)
      
      ## ----decompGis1_1----------------------------------------------------------
      # gene set score decomposition
      # we explore two gene sets, the first one
      gs1 <- top.gs.name[1] # select the most significant gene set
      gs1
      
      ## ----decompGis1_dc, fig.cap="gene set score (GSS) decomposition. The GSS decomposition are grouped according to the tissue of origin of cell lines. The vertical bar showing the 95\\% of confidence interval of the means."----
      # decompose the gene set score over datasets
      decompose.gs.group(mgsa1, gs1, group = cancerType) 
      
      ## ----decompGis1_gis, fig.cap="The gene influential score (GIS) plot. the GIS are represented as bars and the original data where the gene is from is distingished by different colors."----
      gis1 <- GIS(mgsa1, gs1, barcol = gray.colors(4)) # gene influential score
      head(gis1) # print top 6 influencers
      
      ## ----decompGis2, fig.cap=c("Data-wise decomposed GSS for gene set 'PUJANA ATM PCC NETWORK'", "GIS plot for gene set 'PUJANA ATM PCC NETWORK'")----
      # the section gene set
      gs2 <- "PUJANA_ATM_PCC_NETWORK"
      decompose.gs.group(mgsa1, gs2, group = cancerType, x.legend = "topright")
      gis2 <- GIS(mgsa1, "PUJANA_ATM_PCC_NETWORK", topN = 6, barcol = gray.colors(4))
      gis2
      
      ## ----gsSpace, fig.cap="cell line and gene sets projected on the PC1 and PC2"----
      fs <- getmgsa(mgsa1, "fac.scr") # extract the factor scores for cell lines (cell line space)
      layout(matrix(1:2, 1, 2))
      plot(fs[, 1:2], pch=20, col=colcode, axes = FALSE)
      abline(v=0, h=0)
      legend("topright", col=unique(colcode), pch=20, legend=unique(cancerType), bty = "n")
      plotGS(mgsa1, label.cex = 0.8, center.only = TRUE, topN = 0, label = c(gs1, gs2))
      
      ## ----moa, fig.width=6, fig.cap="cell line and gene sets projected on the PC1 and PC2"----
      # perform multivariate analysis
      ana <- moa(NCI60_4arrays, proc.row = "center_ssq1", w.data = "inertia", statis = TRUE)
      slot(ana, "partial.eig")[, 1:6] # extract the eigenvalue
      # show the eigenvalues in scree plot:
      layout(matrix(1:2, 1, 2)) 
      plot(ana, value="eig", type = 2, n=20, main="variance of PCs") # use '?"moa-class"' to check the help manu
      plot(ana, value="tau", type = 2, n=20, main="Scaled variance of PCs")
      
      ## ----moasup----------------------------------------------------------------
      mgsa2 <- mogsa(x = ana, sup=NCI60_4array_supdata, nf=3)
      identical(mgsa1, mgsa2) # check if the two methods give the same results
      
      ## ----prepGraphite----------------------------------------------------------
      library(graphite)
      keggdb <- prepGraphite(db = pathways("hsapiens", "kegg")[1:50], id = "symbol")
      keggdb <- lapply(keggdb, function(x) sub("SYMBOL:", "", x))
      keggdb[1:2]
      
      ## ----prepMsigDB------------------------------------------------------------
      dir <- system.file(package = "mogsa")
      preGS <- prepMsigDB(file=paste(dir, "/extdata/example_msigdb_data.gmt.gz", sep = ""))
      
      ## ----prepInput-------------------------------------------------------------
      # the prepare
      sup_data1 <- prepSupMoa(NCI60_4arrays, geneSets=keggdb, minMatch = 1)
      mgsa3 <- mogsa(x = NCI60_4arrays, sup=sup_data1, nf=3, proc.row = "center_ssq1", w.data = "inertia", statis = TRUE)
  - entryname: render_mogsa.R
    writable: false
    entry: "\nlibrary(\"rmarkdown\")\n\nrender(\"mogsa_gene_set_analysis.R\", \"html_document\"\
      )"
- class: InlineJavascriptRequirement

inputs: []

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
id: david.roberson/build-multiomic-analysis-mogsa/mogsa/0
sbg:appVersion:
- v1.1
sbg:content_hash: a378b486d6da0ffd8f9512f9f4fe108ec690a7888ca9b59e84c551e109016b32f
sbg:contributors:
- david.roberson
sbg:copyOf: david.roberson/mogsa-workflow-wrapping/mogsa/15
sbg:createdBy: david.roberson
sbg:createdOn: 1619459412
sbg:id: david.roberson/build-multiomic-analysis-mogsa/mogsa/0
sbg:image_url:
sbg:latestRevision: 0
sbg:modifiedBy: david.roberson
sbg:modifiedOn: 1619459412
sbg:project: david.roberson/build-multiomic-analysis-mogsa
sbg:projectName: 'BUILD: Multiomic Analysis MOGSA'
sbg:publisher: sbg
sbg:revision: 0
sbg:revisionNotes: Copy of david.roberson/mogsa-workflow-wrapping/mogsa/15
sbg:revisionsInfo:
- sbg:modifiedBy: david.roberson
  sbg:modifiedOn: 1619459412
  sbg:revision: 0
  sbg:revisionNotes: Copy of david.roberson/mogsa-workflow-wrapping/mogsa/15
sbg:sbgMaintained: false
sbg:validationErrors: []
