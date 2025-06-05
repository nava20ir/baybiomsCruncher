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

#COPY  ./omicsViewer /usr/local/lib/R/site-library/omicsViewer

#ADD "https://www.random.org/cgi-bin/randbyte?nbytes=10&format=h" skipcache # to force docker file to ignore cache

CMD ["R", "-e", "shiny::runApp('/home/shiny/app', host = '0.0.0.0', port = 3839)"]
