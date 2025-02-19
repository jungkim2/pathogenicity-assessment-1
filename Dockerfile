FROM rocker/tidyverse:4.2
MAINTAINER naqvia@chop.edu
WORKDIR /rocker-build/

RUN RSPM="https://packagemanager.rstudio.com/cran/2022-10-07" \
  && echo "options(repos = c(CRAN='$RSPM'), download.file.method = 'libcurl')" >> /usr/local/lib/R/etc/Rprofile.site

COPY scripts/install_bioc.r .
COPY scripts/install_github.r .

## install wget
RUN apt update -y && apt install -y wget bzip2 libbz2-dev

RUN wget https://github.com/samtools/bcftools/releases/download/1.17/bcftools-1.17.tar.bz2  && \
    tar -xvjf bcftools-1.17.tar.bz2 && rm -f bcftools-1.17.tar.bz2 && \
    cd bcftools-1.17 && \
    make && mv /rocker-build/bcftools-1.17/bcftools /bin/.

# install R packages
RUN ./install_bioc.r \
    Biobase \
    BiocManager \
    optparse \
    vcfR \
    vroom

ADD Dockerfile .
