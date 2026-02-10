# the way chen did it is that MQCruncher is sitting on the top of omicsViewer. package I first deconvoluted omicsViewer as explained in wiki
source('omicsViewer.R')
source('run_DIA_gui.R') # to calculte maxLFQ using iq and iBAQ with DIA-GUI package
source('map2expr.R') # for mapping raw files to experiment names


pars <- yaml::read_yaml("/home/shiny/app/lims.yaml")
c2n <- function (x)
{
    v <- as.numeric(as.character(x))
    if (is.na(v))
        return(NULL)
    v
}





formatDT <- function (tab, sel = c("single", "multiple")[1], pageLength = 15)
{
    dt <- DT::datatable(na2char(tab), selection = sel, rownames = FALSE,
        filter = "top", class = "table-bordered compact nowrap",
        options = list(scrollX = TRUE, pageLength = pageLength,
            dom = "tip", columnDefs = list(list(targets = unname(which(sapply(tab,
                inherits, c("factor", "character")))) - 1, render = DT::JS("function(data, type, row, meta) {",
                "return type === 'display' && data.length > 30 ?",
                "'<span title=\"' + data + '\">' + data.substr(0, 30) + '...</span>' : data;",
                "}")))))
    DT::formatStyle(dt, columns = (1:ncol(tab)) - 1, fontSize = "90%")
}



getMQParams <- function (x)
{
    mqpar <- fxml_importXMLFlat(x)
    method <- "LF"
    if (length(grep("TMT[0-9]*plex", mqpar$value.)) > 1)
        method <- "TMT"
    ifa <- which(mqpar$elem. == "fastaFilePath")
    if (length(ifa) == 0)
        ifa <- which(mqpar$elem. == "string" & mqpar$level2 ==
            "fastaFiles")
    list(version = mqpar$value.[which(mqpar$elem. == "maxQuantVersion")],
        fasta = mqpar$value.[ifa], enzymes = na.omit(mqpar$value.[which(mqpar$level4 ==
            "enzymes")]), varMod = na.omit(mqpar$value.[which(mqpar$level4 ==
            "variableModifications")]), fixedMod = na.omit(mqpar$value.[which(mqpar$level4 ==
            "fixedModifications")]), mainTol = mqpar$value.[which(mqpar$elem. ==
            "mainSearchTol")], MBR = mqpar$value.[which(mqpar$elem. ==
            "matchBetweenRuns")], peptideFdr = mqpar$value.[which(mqpar$elem. ==
            "peptideFdr")], proteinFdr = mqpar$value.[which(mqpar$elem. ==
            "proteinFdr")], ibaq = mqpar$value.[which(mqpar$elem. ==
            "ibaq")], lfqMode = mqpar$value.[which(mqpar$elem. ==
            "lfqMode")], method = method)
}




getOutliersUp  <- function (x, na.replace = quantile(x, 0.25, na.rm = TRUE), threshold = 0.1)
{
    sec <- na.replace
    f <- apply(x, 1, function(x, sec) {
        x <- sort(x, decreasing = TRUE)
        if (!is.na(x[2]))
            sec <- x[2]
        c(max(x[1] - sec, 0), names(x[1]))
    }, sec = sec)
    df <- data.frame(fold.change.log10 = as.numeric(f[1, ]),
        sample = f[2, ], stringsAsFactors = FALSE)
    df[which(abs(df$fold.change.log10) < threshold), ] <- NA
    df
}


getQCStats <- function (x, method='custom')
{
    print('Running getQCStats with ') # this works fine
    print(method)
    nval <- colSums(!is.na(x))
    nvalCum <- colSums(rowCumsums(apply(!is.na(x), 2, as.integer)) >
        0)
    nvalInt <- colSums(!is.na(rowCumsums(x)))
    expr0 <- x


    if ((method == 'custom') | (method == 'chen-meng')) {
     print('Using the arbitrary method from Chen-meng')
     expr0[is.na(expr0)] <- min(expr0, na.rm = TRUE) - log10(2)
     } else {
     print('using Perseus method')
     expr0 <- impute_perseus(x)
     }
    


    if (ncol(x) <= 2) {
        r1 <- r2 <- NULL
    }
    else {
        iir <- which(!is.na(rowSums(x)))
        if (length(iir) > 2) {
            expra <- x[iir, ]
            r1 <- exprspca(expra, fillNA = FALSE, method = method,
                prefix = "")
        }
        else r1 <- NULL
        r2 <- exprspca(expr0, fillNA = FALSE, method = method, prefix = "")
    }
    list(nval = nval, nvalCum = nvalCum, nvalInt = nvalInt, pcNoImp = r1$samples,
        pcImp = r2$samples)
}


getSteps <- function (fdata, stepList)
{
    sts <- lapply(stepList, function(cc) {
        cc <- unlist(cc)
        pname <- cc[1]
        gname <- cc[-1]
        n <- length(gname)
        m1 <- fdata[[paste("mean", pname, gname[1], sep = "|")]]
        mn <- fdata[[paste("mean", pname, gname[n], sep = "|")]]
        mm <- abs(m1 - mn)
        fillFalse <- function(x) {
            x[is.na(x)] <- FALSE
            x
        }
        s <- lapply(seq_len(length(gname) - 1), function(i) {
            x <- gname[c(i, i + 1)]
            ifdr <- sprintf("ttest|%s_vs_%s|fdr", x[1], x[2])
            imd <- sprintf("ttest|%s_vs_%s|mean.diff", x[1],
                x[2])
            list(signif = cbind(fdata[, ifdr] < 0.01, fdata[,
                ifdr] < 0.05, fdata[, ifdr] < 0.1), down = fillFalse(fdata[,
                imd] > 0), up = fillFalse(fdata[, imd] < 0))
        })
        down <- lapply(s, "[[", "down")
        down <- Reduce("&", down)
        up <- lapply(s, "[[", "up")
        up <- Reduce("&", up)
        sig <- lapply(s, "[[", "signif")
        sig <- Reduce("+", sig)
        sig[!(up | down), ] <- 0
        sig[down, ] <- -sig[down, ]
        colnames(sig) <- paste("x.fdr", c("01", "05", "1"), sep = ".")
        res <- cbind(sig, y.mean.diff = mm)
        colnames(res) <- paste(paste(cc, collapse = ":"), colnames(res),
            sep = "|")
        res
    })
    r <- do.call(cbind, sts)
    colnames(r) <- paste("step", colnames(r), sep = "|")
    r
}


getUPRefProteomeID <- function (domain = c("Eukaryota", "Archaea", "Bacteria", "Viruses")[1])
{
    url <- curl(sprintf("https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/%s/",
        domain))
    on.exit(close(url))
    con <- readLines(url)
    na.omit(stringr::str_match(con, "\\>\\s*(.*?)\\s*/\\<")[,
        2])
}



input_popup <- function (id, pars)
{
    moduleServer(id, function(input, output, session) {
        library(rhandsontable)
        library(shinyFiles)
        ns <- session$ns
        rt <- structure(normalizePath(unlist(pars$project_dir)),
        names = names(pars$project_dir))
        # show modal for the input files DIA


        # Example data for Excel-like component
        table_data <- reactiveVal(data.frame(
        mapping = c("Control", "Treated1", "Treated2"),
        name = c('rawfile_1.raw', 'rawfile2.raw', 'rawfile3.raw')
        ))

        # Render the table inside the modal
        output$excel_table <- renderRHandsontable({
        rhandsontable(table_data(), rowHeaders = NULL)
        })

        # Capture edits from the user
        observe({
        if (!is.null(input$excel_table)) {
            table_data(hot_to_r(input$excel_table))
        }
        })


        showModal(
        modalDialog(
            checkboxInput(ns("checkbox_diann"), "DIA-NN work flow", value = FALSE),
            br(),
            conditionalPanel(
            condition = sprintf("input['%s']", ns("checkbox_diann")),
            radioButtons(
                inputId = ns("method"),
                label = "Choose method",
                choices = c("DIA-NN ProteinGroup file" = "pg", "DIA-NN Report file" = "psm"),
                selected = "pg"
            )
            ),
            br(),
            conditionalPanel(
            condition = sprintf("input['%s']", ns("checkbox_diann")),
            # --- New Excel-like Table ---
            tags$h4("make mapping file between experiment and raw-files"),
            rHandsontableOutput(ns("excel_table")),
            #br(),
            #actionButton(ns("save_table"), "Save Mapping table between raw-file and experiments")
            ),
            #br(),
	        #conditionalPanel(
            #condition = sprintf("input['%s']", ns("checkbox_diann")),
            #shinyFilesButton(
            #id = ns("mappingFile"),
            #label = "Select mapping csv file between RAW and experiments",
            #title = "experiment CSV file",
            #multiple = FALSE
            #)
            #),
	        br(),
            conditionalPanel(
            condition = sprintf("input['%s']", ns("checkbox_diann")),
            checkboxInput(ns("remove_CON"), "Remove contaminats from the output", value = TRUE)
            ),

	        br(),
            conditionalPanel(
            condition = sprintf("input['%s']", ns("checkbox_diann")),
            shinyFilesButton(
            id = ns("diannReportfile"),
            label = "Select diann input",
            title = "DIA-NN input file",
            multiple = FALSE
            )
            ),
            br(),
            conditionalPanel(
            condition = sprintf(
                "input['%s'] && input['%s'] === 'psm'",
                ns("checkbox_diann"),
                ns("method")
            ),
            shinyFilesButton(
            id = ns("fastaFile"),
            label = "Select FASTA file",
            title = "FASTA file",
            multiple = TRUE
            )
            ),
            br(),
            conditionalPanel(
            condition = sprintf("input['%s']", ns("checkbox_diann")),
            actionButton(ns("run_diagui"), "  RUN DIA preparation ")
            ),
	        br(),
  	        conditionalPanel(
            condition = sprintf("input['%s']", ns("checkbox_diann")),
            verbatimTextOutput(ns("log_output"))
            ),
            shinyFilesButton(
            id = ns("mqparOrYaml"),
            label = "Select input file",
            title = "Acceptable file/format: mqpar.xml/.txt/.tsv",
            multiple = FALSE
            ),


            title = "Loading project ...",
            footer = tagList(modalButton("Cancel")),
            size = "m",
            easyClose = FALSE,
            fade = TRUE
            )
            )

        #shinyFileChoose(
        #input = input,
        #id = "mappingFile",
        #roots = rt,
        #defaultRoot = names(rt)[1],
        #session = session,
        #filetypes = c("", "csv"),
        #restrictions = c("AnnotDB", "R-Portable-viewer")
        #)

        shinyFileChoose(
        input = input,
        id = "diannReportfile",
        roots = rt,
        defaultRoot = names(rt)[1],
        session = session,
        filetypes = c("", "tsv"),
        restrictions = c("AnnotDB", "R-Portable-viewer")
        )
        # this is the fasta file selection tool
        shinyFileChoose(
        input = input,
        id = "fastaFile",
        roots = rt,
        defaultRoot = names(rt)[1],
        session = session,
        filetypes = c("", "fasta","fa"),
        restrictions = c("AnnotDB", "R-Portable-viewer")
        )
        shinyFileChoose(
        input = input,
        id = "mqparOrYaml",
        roots = rt,
        defaultRoot = names(rt)[1],
        session = session,
        filetypes = c("", "xml", "txt", "tsv"),
        restrictions = c("AnnotDB", "R-Portable-viewer")
        )

        #observeEvent(input$save_table, {
        #cat("Updated table:\n")
        #print(table_data())
        #})
  
        # Add download handler — triggers file save dialog
        #output$download_table <- downloadHandler(
        #    filename = function() {
        #    paste0("my_table_", Sys.Date(), ".csv")
        #    },
        #    content = function(file) {
        #    write.csv(table_data(), file, row.names = FALSE)
        #    }
        #)

        selected_raw2expr <- reactive({
            req(input$mappingFile)
            req(!inherits(input$mappingFile, "integer"))
            file_info <- parseFilePaths(rt, input$mappingFile)
            req(nrow(file_info) > 0)
            log_content("mapping file uploaded")
            as.character(file_info$datapath[1])
        })

        selected_diann_tsv <- reactive({
            req(input$diannReportfile)
            req(!inherits(input$diannReportfile, "integer"))
            file_info <- parseFilePaths(rt, input$diannReportfile)
            req(nrow(file_info) > 0)
	    log_content("diann report file uploaded")
            as.character(file_info$datapath[1])
        })


        selected_fasta_file <- reactive({
            req(input$fastaFile)
            req(!inherits(input$fastaFile, "integer"))
            file_info <- parseFilePaths(rt, input$fastaFile)
            req(nrow(file_info) > 0)
	    log_content("FASTA file(s) uploaded")
            lapply(file_info$datapath,
            function(x)normalizePath(as.character(x)))
        })

        log_content <- reactiveVal("inactive")
        output$log_output <- renderText({
            log_content()
        })

        # the code to run dia-GUI should be put in here 
        observeEvent(input$run_diagui, {

        req(selected_diann_tsv())  # your diann tsv file
	    #req(selected_raw2expr())
        req(table_data())
        #cat(table_data())
        print(table_data())
	    log_content('wait')
	    log_file = file.path(dirname(selected_diann_tsv()), 'dia_gui_log.txt')
	    con <- file(log_file, open = "a")
	    sink(con, type = 'output')            # redirect output
	    sink(con, type = "message")  # redirect messages

        if (input$method == 'psm') {
        req(selected_fasta_file())  # your FASTA file(s)

	    cat("DIA-NN GUI started at ", Sys.time(), "\n")
	    cat("making FASTA  ", Sys.time(), "\n")
        #mapping_df = read_mapping_raw2expr(selected_raw2expr())
        mapping_df = table_data()

        fasta_combined <- unlist(lapply(selected_fasta_file(), readLines))
        clean_lines <- fasta_combined[nzchar(trimws(fasta_combined))]
	    
	    cat("writing FASTA finsihed  ", Sys.time(), "\n")
        log_content("Writing FASRA...\nProcessing...")

        # Write to a single output file
        writeLines(clean_lines, con =file.path(dirname(selected_diann_tsv()), 'dia_gui.fasta'))
        cat(dirname(selected_diann_tsv()))
	    cat("starting calculation of iBAQ and maxLFQ  ", Sys.time(), "\n")

        res <- baybioms_report_process(
                    selected_diann_tsv(),
                    fasta=file.path(dirname(selected_diann_tsv()), 'dia_gui.fasta'),
                    mapping_dic = mapping_df 
        )
	    cat("making generated tables  ", Sys.time(), "\n")

	    final_ibaq = res$ibaq
	    final_maxlfq = res$max_lfq

        ibaq_folder_path <- file.path(dirname(selected_diann_tsv()), 'dia_ibaq')
        maxlfq_folder_path <- file.path(dirname(selected_diann_tsv()), 'dia_maxlfq')

        if (!dir.exists(ibaq_folder_path)) {
            dir.create(ibaq_folder_path)
        }

        if (!dir.exists(maxlfq_folder_path)) {
            dir.create(maxlfq_folder_path)
        }

        # remove the contaminant if the check box is clicked
        if (input$remove_CON){
            final_ibaq <- final_ibaq[!grepl(final_ibaq[['Protein.Group']],pattern = 'CON_'),]
            final_maxlfq <- final_maxlfq[!grepl(final_maxlfq[['Protein.Group']],pattern = 'CON_'),]
        }

        write.table(final_ibaq, file = file.path(ibaq_folder_path, 'ibaq.tsv'), sep = "\t", row.names = FALSE, quote = FALSE)
        write.table(final_maxlfq, file = file.path(maxlfq_folder_path, 'maxlfq.tsv'), sep = "\t", row.names = FALSE, quote = FALSE)
    	cat("DIA_GUI job finsised  ", Sys.time(), "\n")
	    log_content("Finished DIA-GUI")
        showNotification("File saved successfully!", duration = 8, type = "message")  

        } else { # only mapping for the proteinGroup
        cat('mapping raw file to experiment')

        tryCatch({
        file_path <- file.path(dirname(selected_diann_tsv()), 'mapped_proteinGroup.tsv')
        cat(file_path)
        #final_df <- make_final_pg(selected_diann_tsv(), selected_raw2expr())
        final_df <- make_final_pg(selected_diann_tsv(), table_data())

        # to remove the contaminants if checkbox is clicked
        
        if (input$remove_CON){
            final_df <- final_df[!grepl(final_df[['Protein.Group']],pattern = 'CON_'),]
        }

        write.table(final_df, file = file_path, sep = "\t", row.names = FALSE, quote = FALSE)
        showNotification("File saved successfully!", duration = 8, type = "message")  
        cat("Mapping finished", Sys.time(), "\n")
	    log_content("Finished mapping")
        }, error = function(e) {
        showNotification(paste("Error while saving file:", e$message), type = "error", duration = 8)
        cat(e$message)
        log_content(e$message)
        })
        }
  	    # End logging
  	    sink(type = "message")
 	    sink()
            close(con)
            
        })



        observeEvent(input$mqparOrYaml, { 
        req(!inherits(input$mqparOrYaml, "integer"))       
         removeModal()
            }
        )

        reactive({
            req(input$mqparOrYaml)
            req(!inherits(input$mqparOrYaml, "integer"))
            v <- do.call(file.path, c(rt[[input$mqparOrYaml$root]],
                unlist(input$mqparOrYaml$files, recursive = FALSE)))
                print('showing normalized path')
                cat(normalizePath(v))
            normalizePath(v)
        })
    })

}

landingPage_module <- function (id, codeTable)
{
    moduleServer(id, function(input, output, session) {
        ns <- session$ns
        hashTab <- read.delim(codeTable, header = FALSE, stringsAsFactors = FALSE)
        m0 <- eventReactive(input$go, {
            id <- trimws(input$id)
            i <- which(hashTab[, 1] %in% id)
            if (length(i) != 1)
                return(NULL)
            i
        })
        output$error <- renderUI({
            HTML("There was an error with your passcode, please try again! <br><hr>")
        })
        output$error.ui <- renderUI({
            req(is.null(m0()))
            uiOutput(ns("error"))
        })
        reactive({
            req(m0())
            hashTab[m0(), 2]
        })
    })
}



landingPage_ui  <- function (id)
{
    ns <- NS(id)
    tagList(tags$style("#landingpageid{margin: auto;}"), h1("Welcome to BayBioMS for proteomics!"),
    absolutePanel(id = "landingpageid", class = "myClass",
    fixed = TRUE, draggable = FALSE, left = 0, right = 0,
    top = "35%", width = 500, height = "auto", wellPanel(uiOutput(ns("error.ui")),
    fluidRow(column(10, passwordInput(ns("id"), "Your passcode",
    placeholder = "Enter your passcode here!")),
    column(2, actionButton(ns("go"), "Go!"), style = "padding-top:25px")))))
}



module_annot <- function (id, obj, annot_dir, config)
{
    moduleServer(id, function(input, output, session) {
        ns <- session$ns
        output$fa <- renderText({
            req(config$dataLoading$fastaPath)
            paste(config$dataLoading$fastaPath, sep = "<br>")
        })
        output$fasta_ui <- renderUI({
            tagList(tags$b("Fasta file:"), verbatimTextOutput(ns("fa")))
        })
        observe({
            req(fd <- obj()$fdata)
            ss <- ""
            if (!is.null(config$gsAnnot$gsAnnotCol))
                ss <- config$gsAnnot$gsAnnotCol
            else if (make.names("Majority protein IDs") %in%
                colnames(fd))
                ss <- make.names("Majority protein IDs")
            updateSelectInput(session, inputId = "fdata_cols",
                choices = c("", colnames(fd)), selected = ss)
        })
        observe({
            if (is.null(config$gsAnnot$gsAnnotFile))
                return(NULL)
            updateSelectInput(session, inputId = "annotdb", selected = basename(config$gsAnnot$gsAnnotFile))
        })
        observe({
            if (is.null(config$gsAnnot$gsAnnotColIsoSep) || !config$gsAnnot$gsAnnotColIsoSep %in%
                c("dot (.)", "dash (-)"))
                return(NULL)
            updateSelectInput(session, inputId = "isosep", selected = config$gsAnnot$gsAnnotColIsoSep)
        })
        output$fdata_tab <- DT::renderDT({
            req(obj()$fdata)
            dt <- formatDTScrollY(obj()$fdata[1:100, , drop = FALSE],
                height = 525)
            if (input$fdata_cols %in% colnames(obj()$fdata))
                dt <- formatStyle(dt, input$fdata_cols, backgroundColor = "#85C1E9")
            dt
        }, server = FALSE)
        afa <- reactiveVal(NULL)
        observeEvent(input$annot_new, {
            afa(TRUE)
        })
        observeEvent(input$annotdb, {
            afa(FALSE)
        })
        annotfile <- reactive({
            if (is.null(input$annotdb))
                return(NULL)
            if (input$annotdb == "")
                return(NULL)
            if (!grepl(".annot$", input$annotdb))
                return(NULL)
            file.path(annot_dir, input$annotdb)
        })
        output$afa_text <- renderText("Create a new annotation database using the protein sequences in fasta file! THIS FUNCTION YET TO BE IMPLEMENTED!")
        annotFile <- reactive({
            req(annotfile())
            ret <- read.delim(annotfile(), stringsAsFactors = FALSE,
                nrow = 100)
            i <- colnames(ret) == "source"
            if (!any(i))
                return(ret)
            cls <- rep("NULL", ncol(ret))
            cls[i] <- "character"
            tab <- read.delim(annotfile(), colClasses = cls)
            stats <- table(tab$source)
            lab <- names(stats)
            names(lab) <- paste0(names(stats), " (", stats, ")")
            attr(ret, "source") <- lab
            ret
        })
        annotMatchCol <- reactiveVal()
        observe({
            req(annotFile())
            ss <- colnames(annotFile())
            s1 <- "ACC"
            if (!is.null(config$gsAnnot$gsAnnotDbCol) && nchar(config$gsAnnot$gsAnnotDbCol) >
                0)
                s1 <- config$gsAnnot$gsAnnotDbCol
            updateSelectInput(session = session, inputId = "amatchcol",
                choices = ss, selected = intersect(s1, ss))
        })
        observe({
            annotMatchCol(input$amatchcol)
        })
        observe({
            req(tab <- annotFile())
            lab <- attr(tab, "source")
            if (!is.null(config$gsAnnot$gsAnnotFile)) {
                default <- config$gsAnnot$gsAnnotSource
            }
            else {
                default <- intersect(lab, c("InterPro", "UniPathway",
                  "KEGG", "GO", "ComplexPortal"))
                if (length(default) == 0)
                  default <- lab
            }
            updateMultiInput(session, inputId = "adbs", choices = lab,
                selected = default)
        })
        output$annotdb_tab <- DT::renderDT({
            req(annotFile())
            dt <- formatDTScrollY(annotFile(), height = 300)
            if (!is.null(annotMatchCol()) && nchar(annotMatchCol()) >
                0)
                dt <- formatStyle(dt, annotMatchCol(), backgroundColor = "#85C1E9")
            dt
        }, server = FALSE)
        output$annotSelected <- renderUI({
            req(!is.null(afa()))
            if (afa()) {
                tagAppendAttributes(verbatimTextOutput(ns("afa_text")),
                  style = "white-space:pre-line;")
            }
            else DT::dataTableOutput(ns("annotdb_tab"))
        })
        upids <- eventReactive(input$getupids, {
            req(input$getupids)
            show_modal_spinner(text = "Connecting to UniProt ... this may take a while!")
            l <- getUPRefProteomeID(domain = input$domain)
            remove_modal_spinner()
            l
        })
        output$avaliableupids <- renderUI({
            req(upids())
            selectInput(inputId = ns("refup"), "Select a reference proteome:",
                choices = upids())
        })
        output$orgname <- renderUI({
            req(upids())
            textInput(inputId = ns("orgn"), "Organism name",
                placeholder = "e.g. Homo_sapiens")
        })
        output$upbttn <- renderUI({
            req(upids())
            actionButton(ns("procUPID"), label = "Process")
        })


        observeEvent(input$updateAnnots,{
            
                af <- list.files(annot_dir, pattern = "annot$")
                names(af) <- sub(".annot$", "", af)
                af <- c(none = "", af)
                updateSelectInput(session, inputId = "annotdb", choices = af)
                showNotification('Updated')

        })
   

        observeEvent(input$procUPID, {
            show_modal_spinner(text = "Downloading reference proteome ... ")
            upf <- downloadUPRefProteome(id = input$refup, domain = input$domain,
                destdir = annot_dir)
            show_modal_spinner(text = "Parsing reference proteome ... ")
            fl <- file.path(annot_dir, upf)
            out <- sub(".dat.gz", paste0(make.names(input$orgn),
                ".annot"), fl)
            gsannot <- parseDatTerm(file = fl, outputFile = out)
            remove_modal_spinner()
            af <- list.files(normalizePath(annot_dir), pattern = ".annot$")
            names(af) <- sub(".annot$", "", af)
            updateTabsetPanel(session, inputId = "gsTabset",
                selected = "Select annotation")
            updateSelectInput(session, inputId = "annotdb", choices = af,
                selected = basename(out))
        })
        observeEvent(input$procDatgz, {
            req(input$annotdat)
            fl <- file.path(annot_dir, input$annotdat)
            out <- sub(".dat.gz", paste0(make.names(input$orgnDat),
                ".annot"), fl)
            show_modal_spinner(text = "Parsing reference proteome ... ")
            gsannot <- parseDatTerm(file = fl, outputFile = out)
            remove_modal_spinner()
            af <- list.files(normalizePath(annot_dir), pattern = ".annot$")
            names(af) <- sub(".annot$", "", af)
            updateTabsetPanel(session, inputId = "gsTabset",
                selected = "Select annotation")
            updateSelectInput(session, inputId = "annotdb", choices = af,
                selected = basename(out))
        })
        eventReactive(input$save, {
            config <- list()
            config$gsAnnotFile <- annotfile()
            config$gsAnnotCol <- ""
            if (input$fdata_cols %in% colnames(obj()$fdata))
                config$gsAnnotCol <- input$fdata_cols
            config$gsAnnotColIsoSep <- input$isosep
            config$gsAnnotDbCol <- annotMatchCol()
            config$gsAnnotDir <- annot_dir
            config$gsAnnotUPID <- input$refup
            config$gsAnnotUPDomain <- input$domain
            config$gsAnnotUPName <- input$orgn
            config$gsAnnotSource <- input$adbs
            v <- obj()
            v$config <- config
            v
        })
    })
}



module_annot_ui <- function (id, annot_dir, viewOnly = FALSE)
{
    af <- list.files(annot_dir, pattern = "annot$")
    names(af) <- sub(".annot$", "", af)
    af <- c(none = "", af)
    afdat <- list.files(annot_dir, pattern = ".dat.gz$")
    afdat <- c("", afdat)
    ns <- NS(id)
    if (viewOnly) {
        tl <- .viewOnlyWidget
    }
    else tl <- absolutePanel(top = 5, right = 20, actionButton(inputId = ns("save"),
        "Save configuration"), style = "z-index: 1111;")
        tagList(absolutePanel(top = -10, left = 275, tags$h2("Protein annotation"),
        style = "z-index: 1111;"), tl,
        fluidRow(column(6, wellPanel(style = "background: white; height: 800px",
        tags$h3("Annotation database"), uiOutput(ns("fasta_ui")),
        tabsetPanel(id = ns("gsTabset"), tabPanel("Functional annotation",
            fluidRow(column(8, selectInput(ns("annotdb"), label = "Select annotation database",
                choices = af, multiple = FALSE, selected = "")),
                column(4, selectInput(ns("amatchcol"), label = "ID column",
                  choices = "", multiple = FALSE, selected = ""))),
            fluidRow(column(8,actionButton(ns("updateAnnots"), "Update annotation database"))),
            fluidRow(),
            tags$b("Selected annotation database:"), tabsetPanel(tabPanel("Database",
                uiOutput(ns("annotSelected"))), tabPanel("Sources",
                multiInput(inputId = ns("adbs"), label = "Select annotation source:",
                  choices = "")))), tabPanel("Get from Uniprot reference proteome",
            fluidRow(column(7, radioGroupButtons(ns("domain"),
                label = "Domain", choices = c("Eukaryota", "Archaea",
                  "Bacteria", "Viruses"))), column(5, actionButton(ns("getupids"),
                label = "Get available UP ids", style = "margin-top: 25px"),
                align = "right"), column(12, uiOutput(ns("avaliableupids"))),
                column(12, uiOutput(ns("orgname"))), column(12,
                  uiOutput(ns("upbttn")), align = "right"))),
            tabPanel("Process UP dat file", fluidRow(column(12,
                selectInput(ns("annotdat"), label = "Select annotation database",
                  choices = afdat, multiple = FALSE, selected = "")),
                column(12, textInput(inputId = ns("orgnDat"),
                  "Organism name", placeholder = "e.g. Homo_sapiens")),
                column(12, actionButton(ns("procDatgz"), "Process"),
                  align = "right")))))), column(6, wellPanel(style = "background: white; height: 800px",
        tags$h3("Column to annotate"), fluidRow(column(8, selectInput(ns("fdata_cols"),
            label = "Select ID column", choices = NULL, multiple = FALSE)),
            column(4, selectInput(ns("isosep"), label = "Isoform separator",
                choices = c("none", "dot (.)", "dash (-)"), selected = "none"))),
        DT::dataTableOutput(ns("fdata_tab"))))))
        
}



module_cortest <- function (id, obj, config)
{
    moduleServer(id, function(input, output, session) {
        ns <- session$ns
        ptab <- reactive({
            req(obj()$pdata)
            formatDTScrollY(obj()$pdata, height = "635px")
        })
        output$phenoTab <- DT::renderDT({
            req(dt <- ptab())
            if (!is.null(input$corcol) && any(input$corcol %in%
                colnames(obj()$pdata))) {
                dt <- formatStyle(dt, input$corcol, backgroundColor = "#85C1E9")
            }
            dt
        }, server = FALSE)
        numcols <- reactive({
            req(pd <- obj()$pdata)
            colnames(pd)[sapply(pd, is.numeric)]
        })
        observe({
            req(numcols())
            updateMultiInput(session, inputId = "corcol", choices = numcols(),
                selected = NULL)
        })
        observeEvent(input$add, {
            req(numcols())
            updateMultiInput(session, inputId = "corcol", choices = numcols(),
                selected = numcols()[grep(input$mulcorcol, numcols())])
        })
        observe({
            req(numcols())
            v <- unlist(config$cortest$correlationAnalysisCols)
            updateMultiInput(session, inputId = "corcol", choices = numcols(),
                selected = v)
        })
        eventReactive(input$save, {
            if (length(input$corcol) > 500) {
                showModal(modalDialog("Allowing maximum 500 correlation comparisons",
                  title = "Oops!", footer = modalButton("Dismiss")))
                return(NULL)
            }
            config <- list()
            if (!is.null(input$corcol))
                config$correlationAnalysisCols <- input$corcol
            else config$correlationAnalysisCols <- NULL
            v <- obj()
            v$config <- config
            v
        })
    })
}



module_cortest_ui <- function (id, tabName, viewOnly = FALSE)
{
    ns <- NS(id)
    if (viewOnly) {
    tl <- .viewOnlyWidget
    }
    else tl <- absolutePanel(top = 5, right = 20, actionButton(inputId = ns("save"),
    "Save configuration"), style = "z-index: 1111;")
    tagList(tags$style(HTML(sprintf("\n#shiny-tab-%s .multi-wrapper {\n        height: 642px;\n}\n#shiny-tab-%s .non-selected-wrapper {\n        height: 600px;\n}\n#shiny-tab-%s .selected-wrapper {\n        height: 600px;\n}",
    tabName, tabName, tabName))), absolutePanel(top = -10,
    left = 275, tags$h2("Correlation analysis"), style = "z-index: 1111;"),
    tl, fluidRow(column(6, wellPanel(style = "background: white; height: 800px",
    fluidRow(column(8, textInput(inputId = ns("mulcorcol"),
    label = "Select phenotype data columns for correlateion analysis",
    placeholder = "Select multiple columns using regular expression, such as '^IC50'.",
    width = "100%")), column(4, actionButton(ns("add"),
    label = "Add", style = "margin-top: 25px")),
    column(12, multiInput(inputId = ns("corcol"),
        label = NULL, width = "100%", choices = ""))))),
    column(6, wellPanel(style = "background: white; height: 800px",
    tags$b("Phenotype data"), DT::DTOutput(ns("phenoTab"))))))
}



module_input <- function (id, dir, config)
{
    moduleServer(id, function(input, output, session) {
        ns <- session$ns
        conf <- reactive({
            if (length(config$dataLoading) == 0)
                return(NULL)
            config$dataLoading
        })
        path <- reactiveVal()
        observe({
            req(!is.null(conf()$path) || !is.null(dir()))
            if (!is.null(conf()$path))
                path(conf()$path)
            else path(dirname(dir()))
        })
        output$printPath <- renderUI({
            req(path())
            tags$b(paste("MaxQuant folder selected:", sub("/media/share_baybioms/Projects/002_Proteomics/",
                "", path())))
        })
        datReload <- reactiveVal()
        observe({
            req(path())
            if (!file.exists(objPath <- file.path(path(), "ESVProject/obj.RDS")))
                return(NULL)
            datReload(readRDS(objPath))
        })
        mqparNew <- reactive({
            req(path())
            if (grepl(".xml$", dir(), ignore.case = TRUE)) {
                val <- validMQFolder(path())
                if (!val$valid)
                  return(val)
                c(getMQParams(val$mqpar), val)
            }
            else if (grepl(".txt$|.tsv$", dir(), ignore.case = TRUE)) {
                list(method = "undefined", valid = TRUE, filePath = dir())
            }
        })
        mqpar <- reactiveVal()
        observe({
            if (!is.null(datReload()$mqpar))
                mqpar(datReload()$mqpar)
            else mqpar(mqparNew())
        })
        observe({
            req(mqpar()$valid)
            req(!mqpar()$valid)
            showModal(modalDialog(title = "Invalid MQ serach folder",
                "A 'proteinGroups.txt' file cannot be find in this folder!"))
        })
        dat <- reactiveVal()


        observe({
            req(mqpar()$valid)
            req(mqpar()$method %in% c("LF", "TMT"))
            if (!is.null(datReload()))
                dat(datReload())
            else {
                show_modal_spinner()
                d0 <- try(read.proteinGroups(mqpar()$filePath,
                  quant = mqpar()$method))
                if (class(d0) == "try-error") {
                  showModal(modalDialog(title = "Invalid input file!"))
                  return(NULL)
                }
                else dat(d0)
                remove_modal_spinner()
            }
        })
        observe({
            req(mqpar()$valid)
            req(!mqpar()$method %in% c("LF", "TMT"))
            req(datReload())
            dat(datReload())
        })
        observe({
            req(mqpar()$valid)
            req(!mqpar()$method %in% c("LF", "TMT"))
            if (!is.null(datReload()))
                return(NULL)
            tab <- read.delim(mqpar()$filePath, stringsAsFactors = FALSE,
                nrows = 2100)
            ic <- colnames(tab)[sapply(tab, is.numeric)]
            showModal(modalDialog(fluidRow(column(6, textInput(inputId = ns("mulexprcol"),
                label = "Select columns with pattern:", placeholder = "Using regular expression, such as '^Intensity'.",
                width = "100%")), column(2, actionButton(ns("add"),
                label = "Add", style = "margin-top: 25px")),
                column(2, actionButton(ns("deselect_all"), label = "Deselected all",
                  style = "margin-top: 25px", align = "right")),
                column(2, actionButton(ns("select_all"), label = "Select all",
                  style = "margin-top: 25px", align = "right"))),
                multiInput(ns("exprsCol"), "Columns of expression matrix",
                  choices = ic, width = "100%"), title = "Select columns of expression matrix",
                footer = actionButton(ns("exprsColDone"), label = "Done"),
                size = "l"))
        })
        observeEvent(input$add, {
            if (!is.null(datReload()))
                return(NULL)
            tab <- read.delim(mqpar()$filePath, stringsAsFactors = FALSE,
                nrows = 1)
            choi <- colnames(tab)
            updateMultiInput(session, inputId = "exprsCol", selected = grep(input$mulexprcol,
                choi, value = TRUE))
        })
        observeEvent(input$select_all, {
            if (!is.null(datReload()))
                return(NULL)
            tab <- read.delim(mqpar()$filePath, stringsAsFactors = FALSE,
                nrows = 1)
            choi <- colnames(tab)
            updateMultiInput(session, inputId = "exprsCol", selected = choi)
        })
        observeEvent(input$deselect_all, {
            if (!is.null(datReload()))
                return(NULL)
            updateMultiInput(session, inputId = "exprsCol", selected = "")
        })


  
        observeEvent(input$exprsColDone, {
            if (!is.null(datReload()))
                return(NULL)
            removeModal()
            tab <- read.delim(mqpar()$filePath, stringsAsFactors = FALSE)


            expr <- apply(tab[, input$exprsCol, drop = FALSE],
                2, log10)
            expr[is.infinite(expr)] <- NA
            rownames(expr) <- make.names(rownames(tab))
            d <- list(exprs = expr, annot = tab[, setdiff(colnames(tab),
                input$exprsCol), drop = FALSE])
            attr(d, "label") <- input$exprsCol
            dat(d)
        })

       
        pdata <- reactiveVal()
        observeEvent(dat(), {
            if (!is.null(datReload()))
                pdata(dat()$pdata)
            else {
                aa <- phenoTemplate(attr(dat(), "label"), quant = mqpar()$method)
                pdata(phenoTemplate(attr(dat(), "label"), quant = mqpar()$method))
            }
        })
        output$tab <- renderExcel({
            req(pdata())
            excelTable(data = pdata(), columns = data.frame(type = excelR:::get_col_types(pdata())),
                allowDeleteRow = FALSE, allowInsertRow = FALSE, rowHeaders = NULL,allowCut = FALSE,allowDeleteColumn = TRUE,
                columnSorting = FALSE, colHeaders = colnames(pdata()),
                tableHeight = "700px")
        })
        e2r <- function(x, alt) {
            if (!is.null(x))
                r <- excel_to_R(x)
            else r <- alt
            validNumber <- function(x) !grepl("\\D", gsub("\\.|NA",
                "", x))
            ii <- which(sapply(r, function(x) all(validNumber(x),
                na.rm = TRUE)))
            r[ii] <- lapply(r[ii], function(x) as.numeric(as.character(x)))
            r
        }
        output$saveTemplate <- downloadHandler(filename = function() {
            "phenotypeData.tsv"
        }, content = function(file) {
            write.table(e2r(input$tab, pdata()), file, col.names = TRUE,
                row.names = FALSE, quote = FALSE, sep = "\t")
        })
        observeEvent(input$uploadTemplate$datapath, {
            req(path())
            req(ff <- input$uploadTemplate$datapath)
            t0 <- NULL
            if (grepl(".csv$", ff)) {
                t0 <- read.csv(ff, stringsAsFactors = FALSE,
                  header = TRUE)
            }
            else if (grep(".tsv$|.txt$", ff)) {
                t0 <- read.delim(ff, stringsAsFactors = FALSE,
                  header = TRUE)
            }
            else if (grep(".xlsx$", ff)) {
                t0 <- read.xlsx(ff, sheet = 1)
            }
            else {
                showModal(modalDialog(title = "Unknown file type",
                  easyClose = TRUE))
            }
            req(t0)
            pdr <- e2r(input$tab, pdata())
            if (any(!t0$Label %in% pdr$Label) || any(!pdr$Label %in%
                t0$Label)) {
                showModal(modalDialog(title = "Phenotype data is not valid!",
                  easyClose = TRUE))
                req(NULL)
            }
            pdata(t0[match(pdr$Label, t0$Label), ])
        })
        observe({
            req(dat())
            cn <- colnames(dat()$annot)
            if (is.null(datReload())) {
                sl <- cn[unique(c(1:2, grep("gene|protein.id|protein.name|sequence.window",
                  cn, ignore.case = TRUE)))]
            }
            else {
                sl <- intersect(dat()$fdataHeader, cn)
            }
            updateMultiInput(session, inputId = "proteinHeaders",
                choices = colnames(dat()$annot), selected = sl)
        })
        dt <- reactive({
            req(dat()$annot)
            dt <- formatDT(rbind(head(dat()$annot), tail(dat()$annot)),
                pageLength = 10)
        })
        output$annot <- DT::renderDT({
            dt <- dt()
            if (!is.null(input$proteinHeaders) && any(input$proteinHeaders %in%
                colnames(dat()$annot)))
                dt <- formatStyle(dt, input$proteinHeaders, backgroundColor = "#85C1E9")
            dt
        }, server = FALSE)
        observe({
            sl <- conf()$excludeSamples
            updateMultiInput(session, "exclude", choices = pdata()$Label,
                selected = sl)
        })
        eventReactive(input$save, {
            show_modal_spinner()
            config <- list()
            df <- dat()
            df$path <- path()
            ESVProjectPath <- file.path(df$path, "ESVProject")
            df$fdataHeader <- input$proteinHeaders
            df$fdata <- dat()$annot[, input$proteinHeaders, drop = FALSE]
            df$pdata <- e2r(input$tab, pdata())
            df$time <- c(mTime = file.mtime(mqpar()$filePath),
                processTime = Sys.time())
            df$mqpar <- mqpar()
            if (!file.exists(ESVProjectPath))
                dir.create(ESVProjectPath, recursive = TRUE)
            config$path <- df$path
            config$pathESVProject <- ESVProjectPath
            config$fdataHeader <- df$fdataHeader
            config$excludeSamples <- input$exclude
            config$fastaPath <- mqpar()$fasta
            saveRDS(df, file = file.path(ESVProjectPath, "obj.RDS"))
            df$config <- config
            remove_modal_spinner()
            df
        })
    })
}



module_input_ui <- function (id, viewOnly = FALSE)
{
    ns <- NS(id)
    if (viewOnly) {
    tl <- .viewOnlyWidget
    }
    else tl <- absolutePanel(top = 5, right = 20, actionButton(inputId = ns("save"),
    "Save configuration"), style = "z-index: 1111;")
    tagList(tags$style(HTML("\n.modal-body .multi-wrapper {\n        height: 500px;\n}\n.modal-body .non-selected-wrapper {\n        height: 458px;\n}\n.modal-body .selected-wrapper {\n        height: 458px;\n}")),
    absolutePanel(top = -10, left = 275, tags$h2("Data uploading"),
    style = "z-index: 1111;"), tl, uiOutput(ns("printPath")),
    fluidRow(column(6, wellPanel(style = "background: white; height: 800px",
    tags$h3("Experimental design"), tabsetPanel(tabPanel("Edit sample information",
    wellPanel(style = "background: white; border-color: white",
    fluidRow(column(6, align = "right", style = "padding-left:5px; padding-right:5px; padding-top:0px; padding-bottom:0px;",
    downloadButton(outputId = ns("saveTemplate"),
        label = "Save template")), column(6, style = "padding-left:5px; padding-right:5px; padding-top:0px; padding-bottom:0px;",
    fileInput(inputId = ns("uploadTemplate"),
        label = NULL, placeholder = "Upload self-defined phenotype file",
        multiple = FALSE, accept = c(".csv", ".tsv",
        "xlsx", "txt"))), column(12, div(style = "overflow-x: scroll;",
    tags$h2("** Do not re-order  the samples in the bellow table!  **"),
    tags$h2("** you need a column named cv_group to see the cv plots  **"),
    excelOutput(ns("tab"))))))), tabPanel("Exclude samples",
    wellPanel(style = "background: white; border-color: white;",
    multiInput(inputId = ns("exclude"), label = "Select samples to exclude",
    choices = "", selected = NULL, width = "100%")))))),
    column(6, wellPanel(style = "background: white; height: 800px",
    tags$h3("Protein information"), multiInput(inputId = ns("proteinHeaders"),
    label = "Select columns to be included in further steps",
    choices = "", selected = NULL, width = "100%"),
    tags$b("Examples of available information for proteins"),
    tags$br(), DT::DTOutput(ns("annot"))))))
}



module_normalization <- function (id, object, config)
{
    moduleServer(id, function(input, output, session) {
        ns <- session$ns
        conf <- reactive({
            if (length(config$normalization) == 0)
                return(NULL)
            config$normalization
        })
        obj <- reactiveVal()
        observeEvent(object(), {
            obj(object())
        })
        observe({
            req(config$dataLoading$pathESVProject)
            f <- file.path(config$dataLoading$pathESVProject,
                "obj.RDS")
            if (file.exists(f))
                obj(readRDS(f))
        })

        observe({
            req(obj())
            
            selected = conf()$inputData
            print(paste0('intensity method changed: ',selected))
            cc <- c(`Corrected reporter intensity` = "Reporter.intensity.corrected.log10",
                `LFQ intensity` = "LFQ.intensity", `iBAQ intensity` = "iBAQ",
                `Expression matrix` = "exprs")
            cc <- cc[cc %in% names(obj())]
            updateAwesomeRadio(session = session, inputId = "inputData",
                choices = cc, selected = selected)
        })

        observe({
            req(obj())
            selected = conf()$imputationMethod
            print(paste0('imputation method changed: ',selected))
            cc <- c('none','perseus')
            updateAwesomeRadio(session = session,
                inputId = "imputationMethod",
                choices = cc,
                selected = selected)
        })

        pdata <- reactive({
            req(obj())
            r <- obj()$pdata
            if (!is.null(config$dataLoading$excludeSamples)) {
                i <- match(config$dataLoading$excludeSamples,
                  obj()$pdata$Label)
                r <- r[-i, , drop = FALSE]
                attr(r, "excluded") <- i
            }
            r
        })
        eset <- reactiveVal()
        fDataKeep <- reactiveVal()
        observe({
            req(input$inputData)
            expr <- obj()[[input$inputData]]
            if (!is.null(ric <- attr(pdata(), "excluded")))
                expr <- expr[, -ric, drop = FALSE]
            i <- rowSums(!is.na(expr)) > 0
            req(any(i))
            fDataKeep(obj()$fdata[i, , drop = FALSE])
            eset(list(expr = expr[i, , drop = FALSE], fdata = obj()$fdata[i,
                , drop = FALSE]))
        })
        featureTab_include <- reactiveVal()
        featureTab_exclude <- reactiveVal()
        rv <- reactiveVal()
        observeEvent(obj(), {
            req(pdata())
            t0 <- sapply(pdata(), function(x) max(table(x)))
            t0 <- t0[t0 > 1]
            rv(t0)
            ss <- NULL
            if (!is.null(conf()$filterVar)) {
                ss <- conf()$filterVar
                updateCheckboxInput(session, inputId = "phenoVar",
                  value = TRUE)
            }
            updateSelectInput(session, inputId = "phenoCol",
                choices = names(t0), selected = ss)
            if (!is.null(conf()$filterRowMax)) {
                updateCheckboxInput(session, inputId = "rowMax",
                  value = TRUE)
                updateTextInput(session, inputId = "rowMaxPercent",
                  value = conf()$filterRowMax)
            }
            featureTab_include(eset()$fdata)
            featureTab_exclude(eset()$fdata[numeric(0), ])
        })
        observeEvent(input$phenoCol, {
            req(input$phenoCol)
            req(rv())
            nn <- 2
            if (!is.null(conf()$filterVarMinN))
                nn <- conf()$filterVarMinN
            updateSelectInput(session, inputId = "atLeastN",
                choices = 1:rv()[input$phenoCol], selected = nn)
        })
        output$fTab_include <- DT::renderDT({
            req(tab <- eset()$fdata)
            req(!is.null(featureTab_include()))
            formatDTScrollY(featureTab_include(), height = 285)
        })
        output$fTab_exclude <- DT::renderDT({
            req(tab <- eset()$fdata)
            req(!is.null(featureTab_exclude()))
            formatDTScrollY(featureTab_exclude(), height = 285)
        })
        filterVal <- reactiveVal()
        param <- reactiveValues()
        observeEvent(filterVal(), {
            req(i <- filterVal())
            req(is.logical(i))
            featureTab_include(fDataKeep()[which(i), , drop = FALSE])
            featureTab_exclude(fDataKeep()[which(!i), ])
        })
        observe({


            pdata()
            updateSelectInput(session, inputId = "selectPheno",
                choices = c("None", colnames(pdata())))
            updateSelectInput(session, inputId = "tooltip", choices = colnames(pdata()))
            updateSelectInput(session, inputId = "normCol", selected = conf()$normCol)
            cc <- "None"
            if (refcol <- "Reference" %in% colnames(pdata())) {
                if (sum(pdata()$Reference) > 2)
                  cc <- c(cc, "Reference")
            }
            if ("Batch" %in% colnames(pdata())) {
                cc <- c(cc, "Batch mean")
                if (refcol) {
                  if (all(sapply(split(pdata()$Reference, pdata()$Batch),
                    any)))
                    cc <- c(cc, "Batch reference")
                }
            }
            updateAwesomeRadio(session, inputId = "normRow",
                choices = cc, selected = conf()$normRow)
        })
        expr <- reactiveVal()

        featureData <- reactiveVal()
        observeEvent(input$test, {
            v <- NULL
            p.var <- NULL
            n <- input$atLeastN
            if (nchar(input$rowMaxPercent) > 0 && input$rowMax) {
                v <- c2n(input$rowMaxPercent)
                param$max.quantile <- v
                if (v > 0 && v < 1) {
                  v <- quantile(eset()$expr, na.rm = TRUE, probs = v)
                }
            }
            if (input$phenoCol %in% colnames(pdata()) && input$phenoVar) {
                param$pdata.var <- input$phenoCol
                p.var <- pdata()[, input$phenoCol]
            }
            if (is.null(v) && is.null(p.var))
                filterVal(rep(TRUE, nrow(eset()$expr)))
            else {
                param$pdata.var.n <- as.integer(n)
                bv <- filterRow(eset()$expr, max.value = v, var = p.var,
                  min.rep = as.integer(n))
                filterVal(bv)
            }
            emat <- eset()$expr[filterVal(), ]
            efdata <- eset()$fdata[filterVal(), ]
            d0 <- normalizeData(emat, colWise = input$normCol,
                rowWise = input$normRow, ref = pdata()$Reference,
                batch = pdata()$Batch)
            i <- which(rowSums(!is.na(d0)) > 0)
            expr(d0[i, ])
            featureData(efdata[i, ])
        })
        output$download <- downloadHandler(filename = function() {
            paste("proteinGroups", Sys.Date(), ".xlsx", sep = "")
        }, content = function(file) {
            writeTriplet(expr = expr(), pd = pdata(), fd = eset()$fdata,
                file = file, creator = "BayBioMS")
        })


        

        output$cvscatterplot <- renderPlot({
            req(expr())
            req(pdata())

            df = as.data.frame(expr())
            phenodata = as.data.frame(pdata())



            if ('cv_group' %in% colnames(phenodata)){

            try({
                
            df$ID = rownames(df)
            df = df[,c('ID',colnames(df)[colnames(df) %in% phenodata$Label])]

            make_cv_df <- function(df,phenodata,pheno){
            print(pheno)
            if (pheno == 'global') {
            meta_df = phenodata
            }else{
            meta_df = phenodata[phenodata$cv_group %in% pheno,]
            }
            df = df[,c('ID',colnames(df)[colnames(df) %in% meta_df$Label])]
            df_long <- pivot_longer(df,  cols = -ID,  names_to = "variable",  values_to = "intensity")
            df_long <- df_long[!is.na(df_long$intensity), ]
            df_long = merge(df_long,meta_df,by.x='variable',by.y = 'Label')
            df_long <- df_long %>%  group_by(ID, cv_group) %>%  mutate(count = sum(!is.na(intensity))) %>%  ungroup()
            df_long = as.data.frame(df_long)
            df_long = df_long[df_long['count'] >= 2,]
            df_long <- df_long %>%  mutate(intensity = 10 ^ as.numeric(intensity))
            cv_df <- df_long %>%  group_by(ID, cv_group) %>%  summarise(std  = sd(intensity, na.rm = TRUE), mean = mean(intensity, na.rm = TRUE),.groups = "drop" )
            cv_df$cv = cv_df$std / cv_df$mean
            if (!pheno == 'global') {
            cv_df$pheno = pheno
            }else{
            cv_df$pheno = 'global'
            }
            cv_df
            }



            plot_cv <- function(cv_df, phenodata, legend_title = "Pheno") {
            # Create CV dataframe (your existing function)

            # Remove NA and convert CV to percent
            cv_df <- cv_df[!is.na(cv_df$cv), ]
            cv_df$cv_percent <- cv_df$cv * 100
            
            # Check if 'pheno' column exists
            if(!"pheno" %in% colnames(cv_df)) {
                stop("cv_df must contain a 'pheno' column for grouping")
            }
            
            # ggplot cumulative plot grouped by pheno
            p <- ggplot(cv_df, aes(x = cv_percent, color = pheno)) +
                stat_ecdf(geom = "step", size = 1) +  # step line
                xlim(0, 100) +
                ylim(0, 1) +
                labs(
                x = "% CV",
                y = "Cumulative frequency",
                color = legend_title
                ) +
                theme_minimal() +
                theme(
                legend.position = "bottom",
                legend.title = element_text(size = 12),
                legend.text = element_text(size = 10)
                )
            p
            }

            all_cv_df = lapply(unique(phenodata$cv_group),function(x)make_cv_df(df,phenodata,x))
            final_df = do.call('rbind',all_cv_df)
            global_df = make_cv_df(df,phenodata,'global')
            all_df = rbind(final_df,global_df)
            plot_cv(all_df,phenodata)

            }) # closing try statement
            }  # closing if statement for the column checking            


            })







        stats <- reactive({
            req(expr())
            req(input$imputationMethod)

            print('this is a test of expr')


            imputationMethod = input$imputationMethod

            show_modal_spinner(text = "Calculating ...")
            print('calculating PCA')
            r <- getQCStats(expr(),method = imputationMethod)
            remove_modal_spinner()
            r
        })

        output$barplot <- renderPlotly({

            cutnumorchar <- function(x, n = 60, alt = "") {

                if (is.character(x) || is.factor(x)) {
                  message("too many distinct values, not suitable for color mapping!")
                  v <- alt
                }

                else if (is.numeric(x)) {
                  v <- as.character(cut(x, breaks = n, include.lowest = TRUE,
                    dig.lab = 3))
                }

                else stop("cutnumorchar: x needs to be one of objects: numeric, character, factor")
                v
            }


            req(nrow(pdata()) == length(stats()$nval))
            req(input$selectPheno)

            data <- data.frame(
                x = names(stats()$nval),
                y = stats()$nval,
                dec = stats()$nvalInt, inc = stats()$nvalCum,
                col = "# ID", stringsAsFactors = FALSE
                )

            data$x <- factor(data$x, levels = unique(data$x))
            cc <- "gray"

            if (input$selectPheno %in% colnames(pdata())) {
                data$col <- pdata()[, input$selectPheno]
                cc <- nColors(k = length(unique(data$col)))
                if (length(unique(data$col)) > 60)
                  data$col <- cutnumorchar(data$col, alt = "# ID")
            }

            if (is.numeric(data$col))
                data$col <- as.character(data$col)

            fig <- plot_ly(data)

            fig <- add_trace(fig, x = ~x, y = ~y, type = "bar",
                color = ~col, colors = cc)

            if (!is.null(data$dec))
                fig <- add_trace(fig, x = ~x, y = ~dec, type = "scatter",
                  mode = "lines+markers", name = "Shared")

            if (!is.null(data$inc))
                fig <- add_trace(fig, x = ~x, y = ~inc, type = "scatter",
                  mode = "lines+markers", name = "Cummu.")

            if (input$hideLegend) {

                 layout(
                  fig,
                  xaxis = list(title = ""),
                  yaxis = list(title = "ID"),
                showlegend = FALSE

                  )
            }
            else {
                layout(
                    fig,
                    xaxis = list(title = ""),
                    yaxis = list(title = "ID"),
                    legend = list(orientation = "h", x = 0, y = 1.01,yanchor = "bottom")
                  )
            }



        })
        #

        observe({
            req(stats()$pcImp)
            cn <- colnames(stats()$pcImp)
            names(cn) <- sub("^\\|", "", cn)
            updateSelectInput(session, inputId = "pcx", choices = cn,
                selected = cn[1])
            updateSelectInput(session, inputId = "pcy", choices = cn,
                selected = cn[2])
        })
        callModule(plotly_boxplot_module, id = "boxplotly",
            reactive_param_plotly_boxplot = reactive({
                req(expr())
                list(x = expr(), i = input$fTab_include_rows_selected,
                  ylab = "Intensity")
            }), reactive_checkpoint = reactive(TRUE))
        pcax <- reactive({
            c(match(input$pcx, colnames(stats()$pcImp)), match(input$pcy,
                colnames(stats()$pcImp)))
        })

        
        prepPCScatter <- function(r) {
            req(input$pcx)
            req(input$pcy)
            l <- list()
            l$x <- r[, pcax()[1]]
            l$y <- r[, pcax()[2]]
            l$xlab <- sub("^\\|", "", colnames(r)[pcax()[1]])
            l$ylab <- sub("^\\|", "", colnames(r)[pcax()[2]])
            if (input$selectPheno %in% colnames(pdata()))
                l$color <- pdata()[, input$selectPheno]
            l$tooltips <- pdata()[, input$tooltip]
            l
        }
        stats_pca_noimpute <- reactive({
            req(r <- stats()$pcNoImp)
            prepPCScatter(r)
        })
        v_scatter_noimpute <- callModule(plotly_scatter_module,
            id = "pca_noimpute", reactive_param_plotly_scatter = stats_pca_noimpute)
        stats_pca_impute <- reactive({
            req(r <- stats()$pcImp)
            prepPCScatter(r)
        })
        v_scatter_impute <- callModule(plotly_scatter_module,
            id = "pca_impute", reactive_param_plotly_scatter = stats_pca_impute)
        eventReactive(input$save, {
            if (input$test == 0) {
                showModal(modalDialog("Must test a normaliation method before saving!",
                  title = "Please double check ...", footer = modalButton("Close"),
                  size = c("m", "s", "l")[1]))
                return(NULL)
            }
            conf <- list()
            conf$inputData <- input$inputData
            conf$imputationMethod <- input$imputationMethod
            print("this is imputation method for plotly")
            print(conf$imputationMethod)

            conf$normCol <- input$normCol
            conf$normRow <- input$normRow
            if (input$rowMax)
                conf$filterRowMax <- param$max.quantile
            if (input$phenoVar) {
                conf$filterVar <- param$pdata.var
                conf$filterVarMinN <- param$pdata.var.n
            }

            rl <- list(
                pdata = pdata(),
                fdata = featureData(),
                expr = expr())

            saveRDS(rl, file = file.path(config$dataLoading$pathESVProject,
                "objNorm.RDS"))
            rl$mqpar <- obj()$mqpar
            rl$config <- conf
            rl
        })
    })
}



module_normalization_ui <- function (id, viewOnly = FALSE)
{
    ns <- NS(id)
    if (viewOnly) {
        tl <- .viewOnlyWidget
    }
    else tl <- absolutePanel(top = 5, right = 20, actionButton(inputId = ns("save"),
        "Save configuration"), style = "z-index: 1111;")
    tagList(absolutePanel(top = -10, left = 275, tags$h2("Normalization"),
        style = "z-index: 1111;"), tl, sidebarLayout(sidebarPanel(style = "height: 973px; background-color:white",
        tabsetPanel(tabPanel("Setting",
            awesomeRadio(inputId = ns("inputData"),
            label = "Select input data:", choices = "", selected = "",
            inline = FALSE),
            awesomeRadio(inputId = ns("imputationMethod"),
            label = "Select imputation method - none will be only in t-test:", choices = (c('none','custom','perseus')), selected = "perseus",
            inline = FALSE),
            checkboxInput(ns("rowMax"), label = "Filtering proteins according to row max intensity"),
            conditionalPanel("input.rowMax == true", ns = ns,
                textInputIcon(ns("rowMaxPercent"), label = NULL,
                  icon = list("Quantile:"), placeholder = "e.g. '0.1' means row max intensity lower than 10% intensity quantile will be removed")),
            checkboxInput(ns("phenoVar"), label = "Filtering proteins according to missing value in groups"),
            conditionalPanel("input.phenoVar == true", ns = ns,
                selectInput(ns("phenoCol"), label = "Phenotype variable",
                  choices = NULL), selectInput(ns("atLeastN"),
                  label = "At least N in one of the groups",
                  choices = NULL)), awesomeRadio(inputId = ns("normCol"),
                label = "Column-wise normalization", choices = c("None",
                  "Median centering", "Median centering (shared ID)",
                  "Total sum", "median centering + variance stablization"),
                selected = "None", inline = FALSE), awesomeRadio(inputId = ns("normRow"),
                label = "Row-wise normalization", choices = "None",
                selected = "None", inline = FALSE), actionButton(inputId = ns("test"),
                label = "Test"), downloadButton(outputId = ns("download"),
                label = "Export normalized data")), tabPanel("Filter rows",
            tags$b("Kept features"), DT::DTOutput(ns("fTab_include")),
            tags$hr(), tags$b("Excluded features"), DT::DTOutput(ns("fTab_exclude"))))),
        mainPanel(wellPanel(style = "background: white;", fluidRow(column(3,
            selectInput(ns("selectPheno"), label = "Color to distinguish",
                choices = NULL, selected = NULL, multiple = FALSE,
                selectize = TRUE)), column(3, selectInput(ns("tooltip"),
            label = "Tooltips", choices = NULL, selected = NULL,
            multiple = FALSE, selectize = TRUE)), column(3, selectInput(ns("pcx"),
            label = "PC on x axis", choices = NULL, selected = NULL,
            multiple = FALSE, selectize = TRUE)), column(3, selectInput(ns("pcy"),
            label = "PC on y axis", choices = NULL, selected = NULL,
            multiple = FALSE, selectize = TRUE)), column(width = 6,
            tags$b("PCA no imputation"), plotly_scatter_ui(ns("pca_noimpute"),
                height = "366px")), column(width = 6, tags$b("PCA imputation"),
            plotly_scatter_ui(ns("pca_impute"),
            height = "366px")), column(width = 6, tags$b("Intensity distribution"),
            checkboxInput(
            inputId = ns("hideLegend"),
            label   = "Hide legend in plots",
            value   = FALSE
            ), plotly_boxplot_ui(ns("boxplotly"))),
            column(
                width = 6,
                plotOutput(ns("cvscatterplot"), height = "400px")
            ),
            column(width = 6, tags$b("Protein ID"), 
            plotlyOutput(ns("barplot"))))))))
}



module_steps <- function (id, obj, config)
{
    moduleServer(id, function(input, output, session) {
        ns <- session$ns
        ptab <- reactive({
            req(obj()$pdata)
            formatDTScrollY(obj()$pdata, height = "380px")
        })
        output$phenoTab <- DT::renderDT({
            req(dt <- ptab())
            if (!is.null(tsm()) && nrow(tsm()) > 0 && any(tsm()[,
                1] %in% colnames(obj()$pdata))) {
                ii <- intersect(tsm()[, 1], colnames(obj()$pdata))
                dt <- formatStyle(dt, ii, backgroundColor = "#85C1E9")
            }
            dt
        }, server = FALSE)
        opt <- reactive({
            req(p <- obj()$pdata)
            i <- lapply(p, function(x) {
                tb <- table(x)
                names(tb[tb > 1])
            })
            i[sapply(i, length) > 2]
        })
        observe({
            req(opt())
            updateSelectInput(session, inputId = "column", choices = names(opt()))
        })
        vs <- reactive({
            req(input$column)
            opt()[[input$column]]
        })
        output$bucketOutput <- renderUI({
            req(length(vs()) > 2)
            bucket_list(header = NULL, orientation = "horizontal",
                options = sortable_options(style = "padding:1px"),
                add_rank_list(text = "Available group", labels = vs()),
                add_rank_list(text = "Ordered list", labels = NULL,
                  input_id = ns("orderedList")))
        })
        tsm <- reactiveVal(data.frame(Phenotype = character(0),
            group.1 = character(0), group.2 = character(0), group.3 = character(0)))
        observe({
            req(config$steps$steps)
            req(length(config$steps$steps) > 0)
            m <- t(sapply(config$steps$steps, unlist))
            colnames(m) <- c("Phenotype", paste("Group.", 1:(ncol(m) -
                1)))
            tsm(m)
        })
        observe({
            req(tsm())
            req(nrow(tsm()) > 0)
            req(opt())
            i0 <- apply(tsm(), 1, function(i) {
                any(!i[-1] %in% c("", opt()[[i[1]]]))
            })
            if (any(i0))
                tsm(tsm()[!i0, , drop = FALSE])
        })
        observeEvent(input$clearAll, {
            tsm(tsm()[numeric(0), ])
        })
        frbind <- function(x, y) {
            if (!inherits(x, c("character", "numeric", "matrix",
                "data.frame")) || !inherits(y, c("character",
                "numeric", "matrix", "data.frame")))
                stop("x and y should be one of classes - vector, matrix, data.frame")
            if (is.null(ncol(x)))
                x <- matrix(x, nrow = 1)
            if (is.null(ncol(y)))
                y <- matrix(y, nrow = 1)
            xnc <- ncol(x)
            ync <- ncol(y)
            if (xnc > ync) {
                y <- cbind(y, matrix("", nrow = nrow(y), ncol = xnc -
                  ync))
            }
            else if (xnc < ync) {
                x <- cbind(x, matrix("", nrow = nrow(x), ncol = ync -
                  xnc))
            }
            colnames(x) <- colnames(y) <- c("Phenotype", paste0("Group.",
                1:(ncol(x) - 1)))
            data.frame(rbind(x, y), stringsAsFactors = FALSE)
        }
        observeEvent(input$add, {
            req(input$orderedList)
            if (length(input$orderedList) < 3) {
                showModal(modalDialog(title = "Oops", "The ordered list should contains at least 3 values. Please use 't-test' tab to compare two groups!"))
                return(NULL)
            }
            if (is.null(tsm())) {
                nt <- matrix(c(input$column, input$orderedList),
                  nrow = 1)
            }
            else {
                nt <- frbind(tsm(), c(input$column, input$orderedList))
                nt <- unique(nt)
            }
            if (!is.null(tsm()))
                if (nrow(nt) == nrow(tsm())) {
                  showModal(modalDialog("This comparison has been included!",
                    title = "Oops!", footer = modalButton("Dismiss")))
                  return(NULL)
                }
            tsm(nt)
        })
        output$tabTsm <- renderDT({
            tsm()
            req(tab <- tsm())
            if (nrow(tab) > 0)
                rownames(tab) <- 1:nrow(tab)
            dt <- DT::datatable(tab, selection = "multiple",
                rownames = TRUE, caption = NULL, class = "table-bordered compact nowrap",
                options = list(scrollX = TRUE, scrollY = "150px",
                  dom = "t", ordering = FALSE, pageLength = nrow(tab)))
            DT::formatStyle(dt, columns = 1:ncol(tsm()), fontSize = "90%")
        })
        observeEvent(input$clearSelected, {
            req(i <- input$tabTsm_rows_selected)
            tsm(tsm()[-i, , drop = FALSE])
        })
        eventReactive(input$save, {
            if (nrow(tsm()) > 12) {
                showModal(modalDialog("Allowing maximum 12 steps comparisons.",
                  title = "Oops!", footer = modalButton("Dismiss")))
                return(NULL)
            }
            isolate({
                config <- list()
                if (!is.null(tsm()) && nrow(tsm()) > 0) {
                  config$steps <- split(tsm(), row(tsm()))
                }
                else config$steps <- NULL
                v <- list()
                v$config <- config
                v
            })
        })
    })
}



module_steps_ui <- function (id, viewOnly = FALSE)
{
    ns <- NS(id)
    if (viewOnly) {
    tl <- .viewOnlyWidget
    }
    else tl <- absolutePanel(top = 5, right = 20, actionButton(inputId = ns("save"),
    "Save configuration"), style = "z-index: 1111;")
    tagList(absolutePanel(top = -10, left = 275, tags$h2("Find features step up/down"),
    style = "z-index: 1111;"), tl, fluidRow(column(6, wellPanel(style = "background: white; height: 800px",
    selectInput(inputId = ns("column"), label = "Select phenotype",
    choices = NULL, multiple = FALSE, selectize = TRUE),
    uiOutput(ns("bucketOutput")), column(12, actionButton(inputId = ns("add"),
    label = "Add"), align = "right"))), column(6, wellPanel(style = "background: white; height: 800px",
    tags$b("Phenotype data"), DT::DTOutput(ns("phenoTab")),
    br(), tags$b("Selected comparisons"), DT::DTOutput(ns("tabTsm")),
    fluidRow(column(6, actionButton(ns("clearAll"), "Remove ALL"),
    align = "left"), column(6, actionButton(ns("clearSelected"),
    "Remove selected"), align = "right"))))))
}



module_submit <- function (id, obj, config)
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



module_submit_ui <- function (id, viewOnly = FALSE)
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
        label = "ID column for STRING database querry ",
        choices = NULL, multiple = FALSE), selectInput(ns("ptm_cols"),
        label = "Sequence window for PTM motif analysis",
        choices = NULL, multiple = TRUE), tags$b("Outlier analysis"),
        checkboxInput(ns("outlier"), label = "Identify outliers features",
        value = FALSE), ))))
}



module_ttest <- function (id, obj, config)
{
    moduleServer(id, function(input, output, session) {
        ns <- session$ns
        ptab <- reactive({
            req(obj()$pdata)
            formatDTScrollY(obj()$pdata, height = "635px")
        })
        output$phenoTab <- DT::renderDT({
            req(dt <- ptab())
            if (!is.null(tsm()) && nrow(tsm()) > 0 && any(tsm()[,
                1] %in% colnames(obj()$pdata))) {
                ii <- intersect(tsm()[, 1], colnames(obj()$pdata))
                dt <- formatStyle(dt, ii, backgroundColor = "#85C1E9")
            }
            dt
        }, server = FALSE)
        opt <- reactive({
            req(p <- obj()$pdata)
            i <- lapply(p, function(x) {
                tb <- table(x)
                names(tb[tb > 1])
            })
            i[sapply(i, length) > 1]
        })
        observe({
            req(opt())
            updateSelectInput(session, inputId = "column", choices = names(opt()))
        })
        vs <- reactive({
            req(input$column)
            opt()[[input$column]]
        })
        observe(updateSelectInput(session, inputId = "grp1",
            choices = vs()))
        observe(updateSelectInput(session, inputId = "grp2",
            choices = vs()))
        df0 <- data.frame(Phenotype = character(0), Group.1 = character(0),
            Group.2 = character(0), stringsAsFactors = FALSE)
        tsm <- reactiveVal(df0)
        observe({
            req(config$ttest$ttests)
            m <- t(sapply(config$ttest$ttests, unlist))
            colnames(m) <- c("Phenotype", "Group.1", "Group.2")
            tsm(m)
        })
        observe({
            req(tsm())
            req(nrow(tsm()) > 0)
            req(opt())
            i0 <- apply(tsm(), 1, function(i) {
                any(!i[-1] %in% opt()[[i[1]]])
            })
            if (any(i0))
                tsm(tsm()[!i0, , drop = FALSE])
        })
        observeEvent(input$clearAll, {
            tsm(df0)
        })
        observeEvent(input$clearSelected, {
            req(i <- input$tabTsm_rows_selected)
            tsm(tsm()[-i, , drop = FALSE])
        })
        observeEvent(input$add, {
            if (input$grp1 == input$grp2) {
                showModal(modalDialog("The values of group 1 and group 2 should be different!",
                  title = "Oops!", footer = modalButton("Dismiss")))
                return(NULL)
            }
            nt <- rbind(tsm(), data.frame(Phenotype = input$column,
                Group.1 = input$grp1, Group.2 = input$grp2, stringsAsFactors = FALSE))
            nt <- unique(nt)
            if (!is.null(tsm()))
                if (nrow(nt) == nrow(tsm())) {
                  showModal(modalDialog("This comparison has been included!",
                    title = "Oops!", footer = modalButton("Dismiss")))
                  return(NULL)
                }
            tsm(nt)
        })
        observeEvent(input$allPairs, {
            req(input$column)
            cn <- combn(unique(obj()$pdata[, input$column]),
                m = 2)
            nn <- data.frame(Phenotype = input$column, `Group 1` = cn[1,
                ], `Group 2` = cn[2, ], stringsAsFactors = FALSE)
            tsm(unique(rbind(tsm(), nn)))
        })
        output$tabTsm <- DT::renderDataTable({
            tsm()
            req(tab <- tsm())
            if (nrow(tab) > 0)
                rownames(tab) <- 1:nrow(tab)
            dt <- DT::datatable(tab, extensions = "Scroller",
                rownames = TRUE, filter = "top", class = "table-bordered compact nowrap",
                options = list(scrollX = TRUE, scrollY = 395,
                  dom = "t", ordering = F, scroller = TRUE, pageLength = nrow(tab)))
            DT::formatStyle(dt, columns = 1:ncol(tsm()), fontSize = "90%")
        })
        eventReactive(input$save, {
            if (nrow(tsm()) > 32) {
                showModal(modalDialog("Allowing maximum 32 t-test comparisons.",
                  title = "Oops!", footer = modalButton("Dismiss")))
                return(NULL)
            }
            config <- list()
            if (!is.null(tsm()) && nrow(tsm()) > 0)
                config$ttests <- split(tsm(), row(tsm()))
            else config$ttests <- NULL
            v <- list()
            v$config <- config
            v
        })
    })
}



module_ttest_ui <- function (id, viewOnly = FALSE)
{
    ns <- NS(id)
    if (viewOnly) {
    tl <- .viewOnlyWidget
    }
    else tl <- absolutePanel(top = 5, right = 20, actionButton(inputId = ns("save"),
    "Save configuration"), style = "z-index: 1111;")
    tagList(absolutePanel(top = -10, left = 275, tags$h2("Comparisons of t-test"),
    style = "z-index: 1111;"), tl, fluidRow(column(6, wellPanel(style = "background: white; height: 800px",
    fluidRow(column(8, selectInput(inputId = ns("column"),
    label = "Select phenotype", choices = NULL, multiple = FALSE,
    selectize = TRUE)), column(4, actionButton(ns("allPairs"),
    label = "Add all pairs", width = "100%"), align = "right",
    style = "padding-top:23px"), column(6, selectInput(inputId = ns("grp1"),
    label = "Select group 1", choices = NULL, multiple = FALSE,
    selectize = TRUE)), column(6, selectInput(inputId = ns("grp2"),
    label = "Select group 2", choices = NULL, multiple = FALSE,
    selectize = TRUE)), column(8, br()), column(4, actionButton(inputId = ns("add"),
    label = "Add", width = "100%")), column(12, br()),
    column(12, tags$b("Selected comparisons")), column(4,
        br()), column(4, actionButton(ns("clearAll"),
        "Remove ALL", width = "100%")), column(4, actionButton(ns("clearSelected"),
        "Remove selected", width = "100%")), column(12,
        br()), column(12, DT::dataTableOutput(ns("tabTsm")))))),
    column(6, wellPanel(style = "background: white; height: 800px",
    tags$b("Phenotype data"), DT::DTOutput(ns("phenoTab"))))))
}



mqCrunch <- function (config, file = NULL, outputFile = NULL)
{


    if (missing(config) && !is.null(file))
        config <- yaml::read_yaml(file)
    message <- c()
    path <- file.path(config$dataLoading$path, "ESVProject/objNorm.RDS")
    if (file.exists(path))
        obj <- readRDS(path)
    else {
        obj <- readRDS(file.path(config$dataLoading$path, "ESVProject/obj.RDS"))
        i <- 1:nrow(obj$pdata)
        if (!is.null(config$dataLoading$excludeSamples))
            i <- which(!obj$pdata$Label %in% config$dataLoading$excludeSamples)
        pdata <- obj$pdata[i, ]
        expr <- obj[[config$normalization$inputData]][, i, drop = FALSE]
        i <- rowSums(!is.na(expr)) > 0
        expr <- expr[i, ]


        expr <- normalizeData(expr,
            colWise = config$normalization$normCol,
            rowWise = config$normalization$normRow, ref = pdata$Reference,
            batch = pdata$Batch)
        fdata <- obj$annot[i, config$dataLoading$fdataHeader]
        obj <- list(pdata = pdata, fdata = fdata, expr = expr)
    }
    gs <- NULL
    if (!is.null(config$gsAnnot$gsAnnotFile)) {
        gsannot <- read.delim(config$gsAnnot$gsAnnotFile, stringsAsFactors = FALSE,
            sep = "\t")
        if (length(config$gsAnnot$gsAnnotSource) >= 1)
            gsannot <- gsannot[gsannot$source %in% config$gsAnnot$gsAnnotSource,
                ]
        if (config$gsAnnot$gsAnnotDbCol %in% colnames(gsannot) &&
            config$gsAnnot$gsAnnotCol %in% colnames(obj$fdata)) {
            terms <- data.frame(id = gsannot[[config$gsAnnot$gsAnnotDbCol]],
                term = paste(gsannot$source, trimws(gsannot$term),
                  trimws(substr(gsannot$desc, 1, 50)), sep = "_"),
                stringsAsFactors = FALSE)
            id <- strsplit(obj$fdata[[config$gsAnnot$gsAnnotCol]],
                split = ";")
            if (!is.null(config$gsAnnot$gsAnnotColIsoSep) &&
                config$gsAnnot$gsAnnotColIsoSep %in% c("dot (.)",
                  "dash (-)")) {
                if (config$gsAnnot$gsAnnotColIsoSep == "dot (.)")
                  separator <- "\\."
                else separator <- "-"
                id <- lapply(id, function(x) sapply(strsplit(x,
                  separator), "[", 1))
            }
            gs <- try(gsAnnotIdList(idList = id, gsIdMap = terms),
                silent = TRUE)
            if (inherits(gs, "try-error") || ncol(gs) == 0) {
                message <- c(message, "No gene set annotation mapped to the supplied ID. This may results from a wrong annotation file or annotation column.")
                gs <- NULL
            }
        }
    }
    sdbid <- NULL
    if (config$misc$stringDBCol %in% colnames(obj$fdata))
        sdbid <- sapply(obj$fdata[[config$misc$stringDBCol]],
            "[", 1)
    tss <- t(sapply(config$ttest$ttests, unlist))
    if (length(tss) < 3)
        tss <- NULL
    if (length(steps <- config$steps$steps) > 0) {
        steps <- lapply(steps, setdiff, "")
        tl <- lapply(steps, function(cc) {
            pname <- cc[[1]]
            gname <- unlist(cc[-1])
            n <- length(gname)
            t(sapply(1:(n - 1), function(i) {
                c(pname, gname[i:(i + 1)])
            }))
        })
        tl <- do.call(rbind, tl)
        colnames(tl) <- NULL
        tss <- unique(rbind(tss, tl))
    }
    print('calling the imputation method line 1690')
    #print(config)
    if (!is.null(config$normalization$imputationMethod)) {
    imputationMethod <- config$normalization$imputationMethod
    print(paste0('### this is the imputationMethod before Preomics ### :',imputationMethod))
    } else {
    imputationMethod = 'perseus' # 
    print(' #### sth is wrong with the imputation method so we use perseus method as default ###')
    }
    
    dd <- prepOmicsViewer(expr = obj$expr, pData = obj$pdata,
        fData = obj$fdata, PCA = TRUE, pca.fillNA = TRUE, t.test = tss,
        ttest.fillNA = TRUE, method = imputationMethod, stringDB = sdbid, gs = gs, SummarizedExperiment = FALSE)
    gs <- attr(fData(dd), "GS")
    if (length(steps) > 0) {
        fd <- fData(dd)
        sts <- getSteps(fdata = fd, stepList = steps)
        fd <- cbind(fd, sts)
        fData(dd) <- fd
    }
    mt <- config$misc$ptmSeqMotif
    if (!is.null(mt) && mt != "none" && mt %in% colnames(obj$fdata)) {
        fd <- fData(dd)
        adf <- obj$fdata[, config$misc$ptmSeqMotif, drop = FALSE]
        colnames(adf) <- paste("SeqLogo", "All", colnames(adf),
            sep = "|")
        fd <- cbind(fd, adf)
        fData(dd) <- fd
    }
    if (!is.null(config$misc$outlierAnalysis) && config$misc$outlierAnalysis) {
        fd <- fData(dd)
        o_up <- getOutliersUp(obj$expr)
        colnames(o_up) <- paste("Outlier|Up", colnames(o_up),
            sep = "|")
        o_down <- getOutliersDown(obj$expr)
        colnames(o_down) <- paste("Outlier|Down", colnames(o_down),
            sep = "|")
        fData(dd) <- cbind(fd, o_up, o_down)
    }
    ic <- config$cortest$correlationAnalysisCols
    if (!is.null(ic) & any(ic %in% colnames(obj$pdata))) {
        ic <- intersect(ic, colnames(obj$pdata))
        d0 <- obj$pdata[, ic, drop = FALSE]
        c1 <- correlationAnalysis(obj$expr, d0, min.value = 5)
        c2 <- correlationAnalysis(fillNA(obj$expr, method = config$normaliation$imputationMethod), d0, min.value = 5)
        if (ncol(c2) > 0)
            colnames(c2) <- paste0(colnames(c2), ".impute")
        fd <- fData(dd)
        fData(dd) <- cbind(fd, c1, c2)
    }
    vv <- Biobase::fData(dd)
    attr(vv, "GS") <- gs
    Biobase::fData(dd) <- vv
    tsr <- grep("ttest", colnames(vv), value = TRUE)
    if (length(tsr) >= 5) {
        fx <- tsr[5]
        fy <- tsr[4]
    }
    else if (any(grepl("PCA", colnames(vv)))) {
        fx <- grep("PCA", colnames(vv), value = TRUE)[1]
        fy <- grep("PCA", colnames(vv), value = TRUE)[2]
    }
    else fx <- fy <- NULL
    attr(dd, "fx") <- fx
    attr(dd, "fy") <- fy
    attr(dd, "sx") <- "PCA|All|PC1("
    attr(dd, "sy") <- "PCA|All|PC2("
    attr(dd, "message") <- message
    if (is.null(outputFile))
        outputFile <- file.path(config$dataLoading$pathESVProject,
            "omicsViewerObj.RDS")
    saveRDS(dd, file = outputFile)
    invisible(dd)
}



na2char <- function (x)
{
    cc <- sapply(x, inherits, c("character", "factor"))
    for (i in which(cc)) {
        if (!is.character(x[, i]))
            x[, i] <- as.character(x[, i])
        ir <- is.na(x[, i])
        x[ir, i] <- ""
    }
    x
}



parseDatTerm <- function (file, outputDir = NULL, ...)
{
    message("Reading dat file ...")
    d0 <- readLines(file, ...)
    d0 <- split(d0, cumsum(d0 == "//"))
    org <- trimws(sub("OS", "", grep("^OS", d0[[1]], value = TRUE)))
    org <- make.names(paste(org, collapse = ""))
    while (grepl("\\.\\.", org)) org <- gsub("\\.\\.", ".", org)
    fn <- basename(file)
    vn <- paste0("_", gsub("-", "", Sys.Date()), "_", org, ".annot",
        sep = "")
    if (is.null(outputDir))
        outputDir <- dirname(file)
    outputFile <- file.path(outputDir, sub("(.dat|.dat.gz)$",
        vn, fn))
    message("Processing ...")
    dd <- lapply(d0, function(x) {
        dr <- grep("^DR", x, value = TRUE)
        an <- stringr::str_split_fixed(trimws(sub("^DR", "",
            dr)), ";", 4)
        if (nrow(an) == 0)
            return(NULL)
        ac <- stringr::str_split_fixed(trimws(sub("^AC", "",
            grep("^AC", x, value = TRUE))), ";", n = 2)[1]
        name <- stringr::str_split_fixed(trimws(sub("^ID", "",
            grep("^ID", x, value = TRUE))), " ", 2)[1]
        gn <- trimws(sub("^GN", "", grep("^GN", x, value = TRUE)))
        gn <- grep("Name=", gn, value = TRUE)
        if (length(gn) == 0)
            gn <- NA
        else {
            gn <- gsub("Name=|;$", "", strsplit(gn[1], " ")[[1]][1])
        }
        data.frame(ID = name, ACC = ac, geneName = gn, source = an[,
            1], term = an[, 2], desc = an[, 3], stringsAsFactors = FALSE)
    })
    dd <- do.call(rbind, dd)
    uid <- paste(dd$source, dd$term)
    tb <- table(uid)
    i <- which(uid %in% names(tb[tb >= 5]) & uid %in% names(tb[tb <
        0.25 * length(d0)]))
    dd <- dd[i, ]
    for (i in colnames(dd)) dd[[i]][is.na(dd[[i]])] <- "_NA_"
    if (!is.null(dd)) {
        message("Writing table ...")
        write.table(dd, file = outputFile, col.names = TRUE,
            row.names = FALSE, quote = FALSE, sep = "\t")
    }
    invisible(dd)
}



phenoTemplate  <- function (label, quant = c("LF", "TMT", "undefined")[1])
{
    n <- min(str_count(label, "_")) + 1
    tab <- data.frame(Label = label, stringsAsFactors = FALSE)
    if (quant %in% c("TMT", "undefined")) {
        cn <- str_split_fixed(label, "\\.", 2)
        tab$Batch <- cn[, 2]
        tab$Channel <- str_extract(cn[, 1], "(\\d)+")
    }
    else {
        tab$Batch <- "1"
    }
    tab$Reference <- FALSE
    if (n > 1) {
        m <- str_split_fixed(label, pattern = "_", n = n)
        colnames(m) <- paste0("Var", seq_len(ncol(m)))
        tab <- cbind(tab, m)
    }
    it <- which(vapply(tab, is.factor, logical(1)))
    if (length(it) > 0)
        tab[it] <- lapply(tab[it], as.character)
    tab
}



procConfigYaml <- function (cf)
{
    tagList(c(list(tags$b("Project folder"), tags$p(cf$dataLoading$path)),
    list(tags$b("Input folder"), tags$p(cf$normalization$inputData),
    tags$b("Column-wise normalization"), tags$p(cf$normalization$normCol),
    tags$b("Row-wise normalization"), tags$p(cf$normalization$normRow)),
    list(tags$b("Pairs-wise comparisons using t-test"), tags$ul(lapply(sapply(cf$ttest$ttests,
    function(x) {
        sprintf("%s: %s vs %s", x[1], x[2], x[3])
    }), tags$li))), list(tags$b("Correlation analysis"),
    tags$ul(lapply(cf$cortest$correlationAnalysisCols,
        tags$li))), list(tags$b("Steps up/down"), tags$ul(lapply(sapply(cf$steps$steps,
    function(x) paste(setdiff(x, ""), collapse = " -> ")),
    tags$li))), list(tags$b("Input column for STRING database query"),
    tags$p(cf$misc$stringDBCol)), list(tags$b("Outlier analysis"),
    tags$p(cf$misc$outlierAnalysis)), list(tags$b("PTM sequence window column"),
    tags$p(paste(cf$misc$ptmSeqMotif, collapse = " "))),
    list(tags$b("Functional annotation file"), tags$p(sub(paste0(cf$gsAnnot$gsAnnotDir,
    "/"), "", cf$gsAnnot$gsAnnotFile)), tags$b("Columns mapped to functional annotation file"),
    tags$p(cf$gsAnnot$gsAnnotCol))))
}



read.proteinGroups <- function (x, quant = c("LF", "TMT")[1])
{
    func <- read.proteinGroups.lf
    getName <- function(x) {
        nm <- colnames(x$iBAQ)
        if (is.null(nm))
            nm <- colnames(x[[grep("LFQ", names(x))]])
        if (is.null(nm))
            nm <- colnames(x$Intensity)
        nm
    }
    if (quant == "TMT") {
        getName <- function(x) {
            colnames(v$Reporter.intensity.corrected)
        }
        func <- read.proteinGroups.tmt
    }
    v <- func(x)
    attr(v, "label") <- getName(v)
    id <- str_split_fixed(v$annot$Protein.IDs, pattern = ";",
        2)[, 1]
    gn <- str_split_fixed(v$annot$Gene.names, pattern = ";",
        2)[, 1]
    rn <- make.names(paste(gn, id, sep = "_"))
    ss <- names(which(unlist(vapply(v, nrow, FUN.VALUE = integer(1))) ==
        nrow(v$annot)))
    for (i in ss) rownames(v[[i]]) <- rn
    v
}



read.proteinGroups.lf <- function (file)
{
    pg <- read.delim(file, stringsAsFactors = FALSE)
    df <- data.frame(val = c("iBAQ.", "LFQ.intensity.", "Peptides.",
        "Razor...unique.peptides.", "Unique.peptides.", "Sequence.coverage.",
        "Intensity.", "MS.MS.Count.", "MS.MS.count.", "Identification.type."),
        log = c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE,
            FALSE, FALSE, FALSE), stringsAsFactors = FALSE)
    vi <- vapply(df$val, function(x) length(grep(x, colnames(pg))) >
        0, FUN.VALUE = logical(1))
    df <- df[vi, ]
    if (!"Intensity." %in% df$val)
        stop("The proteinGroup.txt table should have at least Intensity columns!")
    #i <- !(grepl("^REV_", pg$Majority.protein.IDs) | grepl("^CON_",
    #    pg$Majority.protein.IDs) | pg$Only.identified.by.site ==
    #    "+")

    i <- !(pg$Reverse =="+" | pg$Potential.contaminant =="+" | pg$Only.identified.by.site =="+")

    
    annot <- pg[i, -grep(paste(df$val, collapse = "|"), colnames(pg))]
    getExpr <- function(x, type = "iBAQ.", log = TRUE, keep.row = NULL) {
        ic <- grep(type, colnames(x), ignore.case = FALSE, value = TRUE)
        ic <- setdiff(ic, "iBAQ.peptides")
        val <- apply(pg[, ic], 2, as.numeric)
        if (log) {
            val <- log10(val)
            val[is.infinite(val)] <- NA
        }
        if (!is.null(keep.row))
            val <- val[keep.row, ]
        colnames(val) <- gsub(type, "", colnames(val))
        val
    }
    ml <- mapply(function(val, log) getExpr(pg, type = val, log = log,
        keep.row = i), val = df$val, log = df$log)
    names(ml) <- gsub("\\.$", "", df$val)
    ml$annot <- annot
    cnames <- colnames(ml$Intensity)
    for (i in names(ml)) {
        tmp <- ml[[i]]
        if (all(colnames(tmp) %in% cnames) && all(cnames %in%
            colnames(tmp)))
            ml[[i]] <- tmp[, cnames]
    }
    if (!is.null(ml$iBAQ)) {
        i <- which(rowSums(ml$iBAQ, na.rm = TRUE) == 0)
        if (length(i) > 0)
            ml <- lapply(ml, function(x) x[-i, ])
        ml$iBAQ_mc <- sweep(ml$iBAQ, 2, matrixStats::colMedians(ml$iBAQ,
            na.rm = TRUE), "-") + median(ml$iBAQ, na.rm = TRUE)
    }
    ml
}



read.proteinGroups.tmt <- function (file, xref = NULL)
{
    ab <- read.delim(file, stringsAsFactors = FALSE)
    #ir <- c(grep("^CON_", ab$Majority.protein.IDs), grep("^REV_",
    #    ab$Majority.protein.IDs), which(ab$Only.identified.by.site ==
    #    "+"))


    ir <- c(
            which(ab$Only.identified.by.site =="+"),
            which(ab$Reverse =="+"),
            which(ab$Potential.contaminant =="+")
            )    

    eSum <- c("Fraction", "Reporter.intensity.corrected", "Reporter.intensity",
        "Reporter.intensity.count")
    ls <- list()
    for (i in eSum) {
        gb <- grep(paste0(i, ".[0-9]*$"), colnames(ab), value = TRUE)
        ls[[i]] <- apply(ab[-ir, gb, drop = FALSE], 2, as.numeric)
        ab[gb] <- NULL
    }
    lsind <- list()
    eInd <- c("Reporter.intensity.corrected", "Reporter.intensity.count",
        "Reporter.intensity")
    for (i in eInd) {
        gb <- grep(i, colnames(ab), value = TRUE)
        lsind[[i]] <- apply(ab[-ir, gb, drop = FALSE], 2, as.numeric)
        ab[gb] <- NULL
    }
    lsind$Reporter.intensity.corrected.log10 <- log10(lsind$Reporter.intensity.corrected)
    lsind$Reporter.intensity.corrected.log10[is.infinite(lsind$Reporter.intensity.corrected.log10)] <- NA
    lsind$annot <- ab[-ir, ]
    lsind$Summed <- ls
    colnames(lsind$Reporter.intensity.corrected.log10) <- make.names(sub("Reporter.intensity.corrected.",
        "", colnames(lsind$Reporter.intensity.corrected.log10)))
    ec <- intersect(eSum, names(lsind))
    for (i in ec) {
        nn <- sub(i, "", colnames(lsind[[i]]))
        nn <- make.names(sub("^.", "", nn))
        colnames(lsind[[i]]) <- nn
    }
    fn <- make.names(xref$label)
    if (!is.null(xref)) {
        lab <- make.names(paste(xref$channel, xref$mix))
        if (!identical(lab, colnames(lsind[[ec[[1]]]])))
            stop("columne does not match!")
        for (ii in ec) colnames(lsind[[ii]]) <- fn
        if (!is.null(lsind$Reporter.intensity.corrected.log10))
            colnames(lsind$Reporter.intensity.corrected.log10) <- fn
        lsind$xref <- xref
    }
    lsind
}



updateTabItemsBadge <- function (session, id, badgeLabel = "new", badgeColor = "green")
{
    ns <- session$ns
    nid <- sprintf("#shiny-tab-%s", ns(id))
    jscode <- sprintf("var els = document.querySelectorAll(\"a[href='%s']\")[0];\n     var els = els.getElementsByClassName(\"badge\")[0];\n     els.className = \"badge pull-right bg-%s\"",
        nid, badgeColor)
    runjs(jscode)
    html(selector = sprintf("a[href='%s'] > .badge", nid), add = FALSE,
        html = badgeLabel)
}



validMQFolder <- function (dir)
{
    l <- list(valid = FALSE)
    i1 <- file.exists(file.path(dir, "mqpar.xml"))
    i2 <- list.files(dir, pattern = "^proteinGroups.txt$", recursive = TRUE,
        full.names = TRUE)
    if (i1 & length(i2) == 1) {
        l$valid <- TRUE
        l$mqpar <- file.path(dir, "mqpar.xml")
        l$txt <- dirname(i2)
        l$filePath <- file.path(l$txt, "proteinGroups.txt")
    }
    l
}



downloadUPRefProteome <- function (id, domain = c("Eukaryota", "Archaea", "Bacteria",
    "Viruses")[1], destdir = "./")
{
    url <- sprintf("https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/%s/%s",
        domain, id)
    url2 <- curl(url)
    on.exit(close(url2))
    con <- readLines(url2)
    con <- na.omit(stringr::str_match(con, paste0(id, "_s*(.*?)\\s*.dat.gz"))[,
        1])
    con <- con[!grepl("additional", con)]
    durl <- paste(url, con, sep = "/")
    download.file(durl, destfile = file.path(destdir, con))
    con
}


formatDTScrollY <- function (tab, sel = c("single", "multiple")[1], height = 500)
{
    dt <- DT::datatable(na2char(tab), extensions = "Scroller",
        selection = sel, rownames = FALSE, filter = "top", class = "table-bordered compact nowrap",
        options = list(scrollX = TRUE, scrollY = height, dom = "ti",
            scroller = TRUE, pageLength = nrow(tab), columnDefs = list(list(targets = unname(which(sapply(tab,
                inherits, c("factor", "character")))) - 1, render = DT::JS("function(data, type, row, meta) {",
                "return type === 'display' && data.length > 30 ?",
                "'<span title=\"' + data + '\">' + data.substr(0, 30) + '...</span>' : data;",
                "}")))))
    DT::formatStyle(dt, columns = 1:ncol(tab), fontSize = "90%")
}



getOutliersDown <- function (x, threshold = 0.1)
{
    f <- apply(x, 1, function(x) {
        x <- sort(x, decreasing = FALSE)
        if (is.na(x[2]))
            return(c(NA, NA))
        c(min(x[1] - x[2], 0), names(x[1]))
    })
    df <- data.frame(fold.change.log10 = as.numeric(f[1, ]),
        sample = f[2, ], stringsAsFactors = FALSE)
    df[which(abs(df$fold.change.log10) < threshold), ] <- NA
    df
}


validMQFolder <- function (dir)
{
    l <- list(valid = FALSE)
    i1 <- file.exists(file.path(dir, "mqpar.xml"))
    i2 <- list.files(dir, pattern = "^proteinGroups.txt$", recursive = TRUE,
        full.names = TRUE)
    if (i1 & length(i2) == 1) {
        l$valid <- TRUE
        l$mqpar <- file.path(dir, "mqpar.xml")
        l$txt <- dirname(i2)
        l$filePath <- file.path(l$txt, "proteinGroups.txt")
    }
    l
}

writeTriplet <- function (expr, pd, fd, file, creator)
{
    #print(expr)
    td <- function(tab) {
        ic <- which(sapply(tab, is.list))
        if (length(ic) > 0) {
            for (ii in ic) {
                tab[, ii] <- sapply(tab[, ii], paste, collapse = ";")
            }
        }
        tab
    }
    wb <- createWorkbook(creator = creator)
    addWorksheet(wb, sheetName = "Phenotype info")
    addWorksheet(wb, sheetName = "Feature info")
    addWorksheet(wb, sheetName = "Expression")
    id <- paste0("ID", 1:nrow(expr))
    writeData(wb, sheet = "Expression", data.frame(ID = id, expr))
    writeData(wb, sheet = "Feature info", td(cbind(ID = id, fd)))
    writeData(wb, sheet = "Phenotype info", td(pd))
    saveWorkbook(wb, file = file, overwrite = TRUE)
}
