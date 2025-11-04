FROM nava20ir/cruncherbase
COPY ./*.R /home/shiny/app/
COPY ./*.yaml /home/shiny/app/
CMD ["R", "-e", "shiny::runApp('/home/shiny/app', host = '0.0.0.0', port = 3839)"]
