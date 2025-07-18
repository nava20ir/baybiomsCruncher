
module_dian_process <- function (id, obj, config)
{
    moduleServer(id, function(input, output, session) {
        ns <- session$ns
        output$fdataTab <- DT::renderDT({
            req(obj()$fdata)
            formatDTScrollY(obj()$fdata, height = 600)
        })
        observe({
            req(fd <- obj()$fdata)
            ss <- make.names("Majority protein IDs")
            if (!is.null(config$misc$stringDBCol))
                ss <- config$misc$stringDBCol
            if (!ss %in% colnames(fd))
                ss <- "none"
            if (!is.null(config$misc$ptmSeqMotif)) {
                s2 <- unlist(config$misc$ptmSeqMotif)
            }
            else {
                s2 <- grep("sequence.window", colnames(fd), value = TRUE,
                  ignore.case = TRUE)
                if (length(s2) == 0)
                  s2 <- NULL
            }
            updateMultiInput(session, inputId = "fdata_cols",
                choices = c("none", colnames(fd)), selected = ss)
            updateMultiInput(session, inputId = "ptm_cols", choices = colnames(fd),
                selected = s2)
            updateCheckboxInput(session, "outlier", value = config$misc$outlierAnalysis)
        })
        eventReactive(input$save, {
            list(stringDBCol = input$fdata_cols, outlierAnalysis = input$outlier,
                ptmSeqMotif = input$ptm_cols)
        })
    })
}



module_dian_process_ui <- function (id, viewOnly = FALSE)
{
        ns <- NS(id)
        if (viewOnly) {
        tl <- .viewOnlyWidget
        }
        else tl <- absolutePanel(top = 5, right = 20, actionButton(inputId = ns("save"),
        "Save configuration"), style = "z-index: 1111;")
        tagList(absolutePanel(top = -10, left = 275, tags$h2("Miscellaneous and submission"),
        style = "z-index: 1111;"), tl, fluidRow(column(6, wellPanel(style = "background: white; height: 800px",
        tags$h3("Protein information"), DT::dataTableOutput(ns("fdataTab")))),
        column(6, wellPanel(style = "background: white; height: 800px",
        tags$h3("Miscellaneous"), selectInput(ns("fdata_cols"),
        label = "ID column for STRING database query",
        choices = NULL, multiple = FALSE), selectInput(ns("ptm_cols"),
        label = "Sequence window for PTM motif analysis",
        choices = NULL, multiple = TRUE), tags$b("Outlier analysis"),
        checkboxInput(ns("outlier"), label = "Identify outliers features",
        value = FALSE), ))))
}