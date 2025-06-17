source('cruncher_functions.R')

cruncher_server <- function (id, ...)
{
    moduleServer(id, function(input, output, session) {
        ns <- session$ns
        session$onSessionEnded(function() {
            if (tolower(Sys.info()["sysname"]) != "windows")
                return(NULL)
            stopApp()
        })
        config_init <- config <- configNull <- reactiveValues()
        updateBadge <- function(id, v = TRUE) {
            if (v)
                updateTabItemsBadge(session, id, badgeLabel = "Saved",
                  badgeColor = "green")
            else updateTabItemsBadge(session, id, badgeLabel = "Skip",
                badgeColor = "maroon")
        }
        validValue <- function(x) {
            if (inherits(x, "logical"))
                return(x)
            !is.null(x) && length(x) > 0 && !x %in% c("", "none")
        }
        validValues <- function(...) {
            al <- list(...)
            all(sapply(al, validValue))
        }
        writeLock <- reactiveVal(0)
        observe({
            req(config$dataLoading$pathESVProject)
            write_yaml(reactiveValuesToList(config), file = file.path(normalizePath(config$dataLoading$pathESVProject),
                "config.yaml"))
            isolate(if (writeLock() == 1)
                file.create(file.path(normalizePath(config$dataLoading$pathESVProject),
                  ".a.lock")))
            if (is.null(config$dataLoading))
                updateTabItemsBadge(session, "loading", badgeLabel = "Pending",
                  badgeColor = "purple")
            else updateBadge("loading", length(config$dataLoading) >=
                3)
            if (is.null(config$normalization))
                updateTabItemsBadge(session, "normalization",
                  badgeLabel = "Pending", badgeColor = "purple")
            else updateBadge("normalization", length(config$normalization) >=
                3)
            if (is.null(config$ttest))
                updateTabItemsBadge(session, "ttest", badgeLabel = "Pending",
                  badgeColor = "purple")
            else updateBadge("ttest", length(config$ttest) >
                0)
            if (is.null(config$cortest))
                updateTabItemsBadge(session, "cortest", badgeLabel = "Pending",
                  badgeColor = "purple")
            else updateBadge("cortest", length(config$cortest) >
                0)
            if (is.null(config$steps))
                updateTabItemsBadge(session, "steps", badgeLabel = "Pending",
                  badgeColor = "purple")
            else updateBadge("steps", length(config$steps$steps) >
                0)
            if (is.null(config$gsAnnot))
                updateTabItemsBadge(session, "annot", badgeLabel = "Pending",
                  badgeColor = "purple")
            else updateBadge("annot", validValues(config$gsAnnot$gsAnnotFile,
                config$gsAnnot$gsAnnotCol, config$gsAnnot$gsAnnotDbCol))
            if (is.null(config$misc))
                updateTabItemsBadge(session, "misc", badgeLabel = "Pending",
                  badgeColor = "purple")
            else updateBadge("misc", validValue(config$misc$stringDBCol) ||
                validValue(config$misc$outlierAnalysis || validValue(config$misc$ptmSeqMotif)))
            if (file.exists(file.path(normalizePath(config$dataLoading$pathESVProject),
                "results.RDS")))
                updateTabItemsBadge(session, id = "viewer", badgeLabel = "Ready",
                  badgeColor = "green")
            else updateTabItemsBadge(session, "viewer", badgeLabel = "Pending",
                badgeColor = "purple")
            isolate(writeLock(1))
        })
        output$style_backgroundColor <- renderUI({
            if (input$tabs != "viewer")
                return(tags$head(tags$style(HTML(".content-wrapper {background-color:#F5F5F5;}"))))
            tags$head(tags$style(HTML(".content-wrapper {background-color:white;}")))
        })
        path <- list(...)
        if (length(path) == 0) {
            mqparPath <- input_popup(id = "pop0", pars = pars)
            output$share.ui <- renderUI(actionButton(inputId = ns("share"),
                label = "Share"))
            observeEvent(input$share, {
                req(mqparPath())
                tabPath <- file.path(pars$hashtab_dir, "table.tsv")
                r <- read.delim(tabPath, header = FALSE, stringsAsFactors = FALSE)
                colnames(r) <- c("v1", "v2")
                if (!is.na(ip <- match(mqparPath(), r[, 2]))) {
                  pass <- r[ip, 1]
                  status <- "Pre-existing"
                }
                else {
                  pass <- password::password(n = 8)
                  uniqueTest <- TRUE
                  while (uniqueTest) {
                    if (pass %in% r[, 1])
                      pass <- password::password(n = 8)
                    else uniqueTest <- FALSE
                  }
                  status <- "Newly created"
                  df0 <- data.frame(pass, mqparPath())
                  colnames(df0) <- colnames(r)
                  r <- rbind(r, df0)
                  write.table(df0, file = tabPath, append = TRUE,
                    sep = "\t", col.names = FALSE, row.names = FALSE,
                    quote = FALSE)
                }
                info <- sprintf("Please save the passcode and share with collaborator:<br>passcode: <b>%s</b><br>Status: %s<br>Time: %s ",
                  pass, status, Sys.time())
                writeLines(info, con = file.path(dirname(mqparPath()),
                  "share.html"))
                showModal(modalDialog(title = "Passcode generated",
                  HTML(info), footer = modalButton("OK")))
            })
        }
        else {
            mqparPath <- path$path
            updateTabItems(session, inputId = "tabs", selected = ns("viewer"))
        }
        observeEvent(mqparPath(), {
            req(mqparPath())
            cff <- file.path(dirname(mqparPath()), "ESVProject",
                "config.yaml")
            if (file.exists(cff)) {
                cf <- read_yaml(cff)
                for (i in names(cf)) config_init[[i]] <- config[[i]] <- cf[[i]]
            }
            if (is.null(config$dataLoading$pathESVProject))
                return(NULL)
            f1 <- file.exists(file.path(normalizePath(config$dataLoading$pathESVProject),
                ".a.lock"))
            f2 <- file.exists(file.path(normalizePath(config$dataLoading$pathESVProject),
                "results.RDS"))
            if (f1 && f2)
                showModal(modalDialog("The paramters have been updated, but results haven't been recalculated. Please check the parameters and submit to run the analysis.",
                  title = "Warning: unmatched parameters and results",
                  footer = modalButton("Close")))
        })
        obj <- reactiveVal()
        d01 <- module_input(id = "input", dir = mqparPath, config = config_init)
        observeEvent(d01(), {
            req(d01())
            for (i in names(config)) {
                config[[i]] <- NULL
            }
            obj(NULL)
            config$dataLoading <- d01()$config
            if (file.exists(er <- file.path(normalizePath(config$dataLoading$pathESVProject),
                "results.RDS")))
                unlink(er)
            updateTabItems(session, inputId = "tabs", selected = ns("normalization"))
        })
        d02 <- module_normalization(id = "body_normalization",
            object = d01, config = config)
        observeEvent(d02(), {
            obj(d02())
            if (file.exists(er <- file.path(normalizePath(config$dataLoading$pathESVProject),
                "results.RDS")))
                unlink(er)
        })
        observe({
            req(config$dataLoading$pathESVProject)
            f <- file.path(normalizePath(config$dataLoading$pathESVProject),
                "objNorm.RDS")
            if (file.exists(f))
                obj(readRDS(f))
        })
        observeEvent(d02(), {
            req(d02()$pdata)
            config$normalization <- d02()$config
            updateTabItems(session, inputId = "tabs", selected = ns("ttest"))
        })
        d03.1 <- module_ttest(id = "stats_ttest", obj = obj,
            config = config)
        observeEvent(d03.1(), {
            req(d03.1())
            config$ttest <- d03.1()$config
            if (file.exists(er <- file.path(normalizePath(config$dataLoading$pathESVProject),
                "results.RDS")))
                unlink(er)
            updateTabItems(session, inputId = "tabs", selected = ns("cortest"))
        })
        d03.2 <- module_cortest(id = "stats_cortest", obj = obj,
            config = config)
        observeEvent(d03.2(), {
            req(d03.2())
            if (!is.null(d03.2()$config))
                config$cortest <- d03.2()$config
            if (file.exists(er <- file.path(normalizePath(config$dataLoading$pathESVProject),
                "results.RDS")))
                unlink(er)
            updateTabItems(session, inputId = "tabs", selected = ns("steps"))
        })
        d03.3 <- module_steps(id = "stats_steps", obj = obj,
            config = config)
        observeEvent(d03.3(), {
            req(d03.3())
            if (!is.null(d03.3()$config))
                config$steps <- d03.3()$config
            if (file.exists(er <- file.path(normalizePath(config$dataLoading$pathESVProject),
                "results.RDS")))
                unlink(er)
            updateTabItems(session, inputId = "tabs", selected = ns("annot"))
        })
        d04 <- module_annot(id = "prot_annot", obj = obj, annot_dir = normalizePath(pars$annot_dir),
            config = config)
        observeEvent(d04(), {
            req(d04())
            config$gsAnnot <- d04()$config
            if (file.exists(er <- file.path(normalizePath(config$dataLoading$pathESVProject),
                "results.RDS")))
                unlink(er)
            updateTabItems(session, inputId = "tabs", selected = ns("misc"))
        })
        d05 <- module_submit(id = "misc", obj = obj, config = config)
        observeEvent(d05(), {
            req(d05())
            if (length(d05()) > 0)
                config$misc <- d05()
            if (file.exists(er <- file.path(normalizePath(config$dataLoading$pathESVProject),
                "results.RDS")))
                unlink(er)
        })
        observeEvent(d05(), {
            req(config$dataLoading$pathESVProject)
            if (any(sapply(reactiveValuesToList(config), is.null))) {
                showModal(modalDialog("Please check all 'pending' tabs to continue. If you want to skip some anlayses, simply remove the defined \n        analysis and click 'save configuration' button on the top right!",
                  title = "Incomplete analysis configuration",
                  footer = modalButton("Dismiss"), size = "m",
                  easyClose = TRUE, fade = TRUE))
                return(NULL)
            }
            ff <- file.path(normalizePath(config$dataLoading$pathESVProject),
                "config.yaml")
            showModal(modalDialog(procConfigYaml(config), title = "Analysis configuration",
                footer = actionButton(inputId = ns("submit2"),
                  label = "Submit"), size = "l", easyClose = TRUE,
                fade = TRUE))
        })
        res <- reactiveVal(NULL)
        observeEvent(input$submit2, {
            removeModal()
            show_modal_spinner()
            v3 <- mqCrunch(file = file.path(normalizePath(config$dataLoading$pathESVProject),
                "config.yaml"), outputFile = file.path(normalizePath(config$dataLoading$pathESVProject),
                "results.RDS"))
            res(v3)
            print('getting the config')
            print(paste0('this is the imputation method implemented: ',config$normalization$imputationMethod))
            remove_modal_spinner()
            if (!is.null(m <- attr(v3, "message"))) {
                showModal(modalDialog(paste(m, sep = "\n"), title = "Possible problems",
                  footer = modalButton("Dismiss")))
            }
            updateTabItems(session, inputId = "tabs", selected = ns("viewer"))
            updateTabItemsBadge(session, id = "viewer", badgeLabel = "Ready",
                badgeColor = "green")
            if (file.exists(lockFile <- file.path(normalizePath(config$dataLoading$pathESVProject),
                ".a.lock")))
                file.remove(lockFile)
        })
        #callModule(omicsViewer:::app_module, id = "app", .dir = reactive({
        callModule(app_module, id = "app", .dir = reactive({

            print('calling omicsViewer')
            res()
            req(config$dataLoading$pathESVProject)
            req(file.exists(file.path(normalizePath(config$dataLoading$pathESVProject),
                "results.RDS")))
            a <- normalizePath(config$dataLoading$pathESVProject)
            print(a)
            a
        }), additionalTabs = NULL, filePattern = "^results*(.*?).RDS$",
            esetLoader = readESVObj, exprsGetter = exprs,imputationMethod = config$normalization$imputationMethod,
            pDataGetter = pData, fDataGetter = fData, defaultAxisGetter = function(x,
                what = c("sx", "sy", "fx", "fy")[1]) attr(x,
                what), appName = NULL, appVersion = NULL)
    })
}
