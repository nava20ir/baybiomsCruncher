FROM mengchen18/mqc_admin:0.3.10
# TODO 
CMD ["R", "-e", "shiny::runApp('/home/shiny/app', host = '0.0.0.0', port = 3840)"]