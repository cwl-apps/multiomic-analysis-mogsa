#FROM bioconductor/bioconductor_docker:RELEASE_3_12
      
#RUN Rscript -e 'BiocManager::install("mogsa")'
#RUN Rscript -e 'install.packages("rmarkdown")'
#RUN Rscript -e 'BiocManager::install("BiocStyle")'
#RUN Rscript -e 'install.packages("magick")'
#RUN Rscript -e 'BiocManager::install("org.Hs.eg.db")'

FROM cgc-images.sbgenomics.com/david.roberson/mogsa:1.24.0

RUN Rscript -e 'install.packages("tidyverse")'
RUN Rscript -e 'BiocManager::install("EnsDb.Hsapiens.v86", update = F, ask = F)'