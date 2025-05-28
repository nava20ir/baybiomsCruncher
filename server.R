#library(MQC)
source('cruncher_server.R')
pars <- yaml::read_yaml("/home/shiny/app/lims.yaml")
server <- function(input, output, session){
     cruncher_server( id = "l0" )
}