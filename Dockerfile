FROM openanalytics/r-base

RUN apt-get update && apt-get install -y \
    sudo \
    pandoc \
    pandoc-citeproc \
    libcurl4-gnutls-dev \
    libcairo2-dev \
    libxt-dev \
    libssl-dev \
    libssh2-1-dev \
    libssl1.0.0

RUN apt-get update -qq && apt-get -y --no-install-recommends install \
  libxml2-dev \
  libsqlite3-dev \
  libmariadbd-dev \
  libmariadb-client-lgpl-dev \
  libpq-dev \
  libssl-dev \
  libcurl4-openssl-dev \
  libssh2-1-dev \
  unixodbc-dev \
  && R -e "source('https://bioconductor.org/biocLite.R')" \
  && install2.r --error \
    --deps TRUE \
    tidyverse \
    dplyr \
    devtools \
    formatR \
    remotes \
    selectr \
    caTools

RUN apt-get install dialog apt-utils -y
RUN apt-get install libssl-dev -y
RUN apt-get install libv8-3.14-dev -y
RUN apt-get install libjpeg-dev -y

RUN R -e "install.packages(c('shiny','rmarkdown'), repos='http://cran.rstudio.com/', dependencies=T)"
RUN R -e "install.packages(c('shinyjs'), repos='http://cran.rstudio.com/', dependencies=T)"
RUN R -e "install.packages(c('openssl'), repos='http://cran.rstudio.com/', dependencies=T)"
RUN R -e "install.packages(c('httr'), repos='http://cran.rstudio.com/', dependencies=T)"

RUN R -e "install.packages(c('ggpubr'), repos='http://cran.rstudio.com/', dependencies=T)"
RUN R -e "source('https://bioconductor.org/biocLite.R'); biocLite('biomaRt')"
RUN R -e "source('https://bioconductor.org/biocLite.R'); biocLite('GenomicFeatures')"
RUN R -e "source('https://bioconductor.org/biocLite.R'); biocLite('Gviz')"
RUN R -e "install.packages(c('DT'), repos='http://cran.rstudio.com/', dependencies=T)"


RUN mkdir /root/GenomeBrowser
COPY GenomeBrowser /root/GenomeBrowser
COPY Rprofile.site /usr/lib/R/etc/

RUN R -e "source('https://bioconductor.org/biocLite.R'); biocLite('regioneR')"
RUN R -e "source('https://bioconductor.org/biocLite.R'); biocLite('GenomicScores')"

EXPOSE 3838

CMD ["R", "-e", "shiny::runApp('/root/GenomeBrowser')"]
 
