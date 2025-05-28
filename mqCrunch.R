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
        expr <- normalizeData(expr, colWise = config$normalization$normCol,
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
    dd <- prepOmicsViewer(expr = obj$expr, pData = obj$pdata,
        fData = obj$fdata, PCA = TRUE, pca.fillNA = TRUE, t.test = tss,
        ttest.fillNA = TRUE, stringDB = sdbid, gs = gs, SummarizedExperiment = FALSE)
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
        c2 <- correlationAnalysis(fillNA(obj$expr), d0, min.value = 5)
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