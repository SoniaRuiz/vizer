FROM openanalytics/r-base
MAINTAINER Sonia Garcia "s.ruiz@ucl.ac.uk"

RUN apt-get update && apt-get install -y --no-install-recommends \
	libxt-dev \
	libcairo2-dev \
        libxml2-dev \
        libssl-dev \
        libcurl4-openssl-dev \
        libv8-3.14-dev \
        libsqlite3-dev \
        libmariadb-client-lgpl-dev \
		&& rm -rf /var/lib/apt/lists/* \
		&& R -e "install.packages(c('shiny','rmarkdown','ggpubr','shinyjs','openssl','httr','DT', 'RMySQL','tidyverse','RSQLite'), repos='http://cran.rstudio.com/', dependencies=T)" \
		&& R -e "source('https://bioconductor.org/biocLite.R'); biocLite('biomaRt'); biocLite('GenomicFeatures'); biocLite('Gviz'); biocLite('regioneR'); biocLite('GenomicScores')" \
		&& mkdir /root/vizER
COPY vizER /root/vizER
COPY Rprofile.site /usr/lib/R/etc/

RUN mkdir /root/vizER/vizER_data

EXPOSE 3838

CMD ["R", "-e", "shiny::runApp('/root/vizER')"]
