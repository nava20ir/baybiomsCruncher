FROM rocker/shiny:4.2.0


RUN apt-get update \
        && apt-get install -y \
        apt-utils \
        manpages-dev \
        libnetcdf-dev \
        libxml2-dev \
        libglpk-dev \
        libnode-dev  \
        libz-dev


RUN mkdir /home/shiny/app
RUN mkdir /home/shiny/annotDb
RUN mkdir /home/shiny/hashtab
RUN mkdir /home/shiny/projects_baybioms
RUN mkdir /home/shiny/projects_baybioms_mri
RUN mkdir /home/shiny/projects_bioinformatics
RUN mkdir /home/shiny/projects_biotyping
RUN mkdir /home/shiny/projects_kusterlab

COPY ./*.R /home/shiny/app/
COPY ./*.yaml /home/shiny/app/



RUN R -e "install.packages('BiocManager')"
#RUN R -e "BiocManager::install('omicsViewer')"




# installing extra packages
RUN Rscript /home/shiny/app/install_packages.R


#RUN R -e 'devtools::install_github("mengchen18/omicsViewer", dependencies = TRUE)'
RUN R -e "install.packages('shiny')"
RUN R -e "install.packages('shinyBS')"
RUN R -e "install.packages('shinydashboard')"

### DIA GUI installation tested inside the container
RUN R -e 'BiocManager::install("pcaMethods")'
RUN R -e 'BiocManager::install("impute")'
RUN R -e 'install.packages("imputeLCMD")'
RUN R -e 'devtools::install_github("mgerault/DIAgui")'

RUN R -e 'install.packages("remotes")'
RUN R -e 'remotes::install_version("vctrs", "0.5.2", repos = "https://cloud.r-project.org")'
RUN R -e 'remotes::install_version("rlang", version = "1.0.6", repos = "https://cloud.r-project.org")'
RUN R -e 'remotes::install_version("lifecycle", version = "1.0.3", repos = "https://cloud.r-project.org")'
RUN R -e 'remotes::install_version("tidyselect", version = "1.2.0", repos = "https://cloud.r-project.org")'
RUN R -e 'remotes::install_version("dplyr", version = "1.1.0", repos = "https://cloud.r-project.org")'


#COPY  ./omicsViewer /usr/local/lib/R/site-library/omicsViewer

#ADD "https://www.random.org/cgi-bin/randbyte?nbytes=10&format=h" skipcache # to force docker file to ignore cache

CMD ["R", "-e", "shiny::runApp('/home/shiny/app', host = '0.0.0.0', port = 3839)"]
