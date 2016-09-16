FROM      ubuntu:14.04

RUN echo "deb http://cran.rstudio.com/bin/linux/ubuntu trusty/" >> /etc/apt/sources.list
RUN gpg --keyserver keyserver.ubuntu.com --recv-key E084DAB9
RUN gpg -a --export E084DAB9 | sudo apt-key add -
RUN apt-get update
RUN apt-get install -y r-base-dev libxml2-dev libssh2-1-dev libcurl4-openssl-dev wget git

WORKDIR /build

RUN wget https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2
RUN tar -xf samtools-1.3.1.tar.bz2
WORKDIR /build/samtools-1.3.1
RUN make
RUN make prefix=/usr/local install

WORKDIR /build

RUN R -e "install.packages(c(\"data.table\",\"gtools\",\"VGAM\",\"devtools\"), \
    repos = \"http://cran.r-project.org\")"

RUN R -e "source(\"http://bioconductor.org/biocLite.R\"); \
    biocLite(\"IRanges\",ask=FALSE); \
    biocLite(\"Rsamtools\",ask=FALSE); \
    biocLite(\"DNAcopy\",ask=FALSE); \
    biocLite(\"Rsubread\",ask=FALSE); \
    biocLite(\"BSgenome.Hsapiens.UCSC.hg19\",ask=FALSE)"

RUN R -e "devtools::install_github(\"mctp/cnatools\")"

WORKDIR /


