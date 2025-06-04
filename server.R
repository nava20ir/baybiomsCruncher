# this the new one
source('cruncher_server.R')
server <- function(input, output, session){
     cruncher_server( id = "l0" )
}