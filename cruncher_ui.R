source('cruncher_functions.R')

cruncher_ui <- function (id, viewOnly = FALSE)
{
    ns <- NS(id)
    dashboardPage(dashboardHeader(title = "Cruncher V2"),
        dashboardSidebar(useShinyjs(), 
        sidebarMenu(id = ns("tabs"),
            menuItem("Data Loading", tabName = ns("loading"),
                icon = icon("file-import"), badgeLabel = "Pending",
                badgeColor = "purple"), 
                
            menuItem("Normalization",
                tabName = ns("normalization"), icon = icon("compress-alt"),
                badgeLabel = "Pending", badgeColor = "purple"),

            menuItem("Stats", tabName = ns("stats"), icon = icon("calculator"),
                startExpanded = TRUE, menuItem("t-test", tabName = ns("ttest"),
                  icon = icon("not-equal"), badgeLabel = "Inactive",
                  badgeColor = "purple"),
                  
            menuItem("correlation",
                  tabName = ns("cortest"), icon = icon("chart-line"),
                  badgeLabel = "Inactive", badgeColor = "purple"),

            menuItem("steps up/down", tabName = ns("steps"),
                  icon = icon("chart-line"), badgeLabel = "Inactive",
                  badgeColor = "purple")),
                  
            menuItem("Functional Annot",
                  tabName = ns("annot"), icon = icon("map-marked-alt"),
                  badgeLabel = "Pending", badgeColor = "purple"),

            menuItem("Miscellaneous", tabName = ns("misc"), icon = icon("play-circle"),
                badgeLabel = "Pending", badgeColor = "purple"),

            menuItem("Results", tabName = ns("viewer"), icon = icon("map-marked-alt"),
                badgeLabel = "Pending", badgeColor = "purple"),
                
            actionButton(ns("reload_btn"), "Reload Page", icon = icon("redo")),
                # and include script binding (use ns to compute id)
                tags$head(tags$script(HTML(sprintf("
                $(document).on('click', '#%s', function(){ location.reload(); });
                ", ns("reload_btn"))))),

            uiOutput(ns("share.ui")))),
            dashboardBody(tabItems(
            
            tabItem(tabName = ns("loading"),
                module_input_ui(id = ns("input"), viewOnly = viewOnly)),

            tabItem(tabName = ns("normalization"),
             module_normalization_ui(ns("body_normalization"),
                viewOnly = viewOnly)), 

            tabItem(tabName = ns("ttest"),
                module_ttest_ui(ns("stats_ttest"), viewOnly = viewOnly)),

            tabItem(tabName = ns("cortest"),
                module_cortest_ui(ns("stats_cortest"),
                tabName = ns("cortest"), viewOnly = viewOnly)),

            tabItem(tabName = ns("steps"),
                module_steps_ui(ns("stats_steps"),
                viewOnly = viewOnly)), tabItem(tabName = ns("annot"),

                module_annot_ui(ns("prot_annot"), annot_dir = normalizePath(pars$annot_dir),
                  viewOnly = viewOnly)), tabItem(tabName = ns("misc"),

                module_submit_ui(ns("misc"), viewOnly = viewOnly)),

            tabItem(tabName = ns("viewer"), tags$style(".container-fluid { margin: 0px; padding: 0px; }"),
                absolutePanel(top = -10, left = 275, tags$h2("Results"),
                  style = "z-index: 1111;"), fluidPage(tags$style(".container-fluid { background-color: white; }"),
                  app_ui(ns("app")))))))
}
