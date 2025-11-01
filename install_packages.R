install.packages("V8")
install.packages("SparseM")
# cran
s1 <- c(
  "survminer",
  "survival",
  "fastmatch",
  "reshape2",
  "beeswarm",
  "grDevices",
  "shinycssloaders",
  "shinythemes",
  "networkD3",
  "httr",
  "RColorBrewer",
  "psych",
  "stringr",
  "shiny",
  "rhandsontable",
  "shinydashboard",
  "shinyWidgets",
  "shinybusy",
  "matrixStats",
  "flatxml",
  "excelR",
  "shinyjs",
  "shinyFiles",
  "DT",
  "ROCR",
  "plotly",
  "openxlsx",
  "yaml",
  "curl",
  "sortable",
  "BiocManager",
  "password",
  "ggseqlogo",
  "devtools",
  "RSQLite",
  "readr"
  )

# # BIOC
s2 <- c(
  "Biobase",
  "fgsea",
  "S4Vectors",
  "SummarizedExperiment"
  )

#
lapply(s1, function(x) {
  if (x %in% installed.packages()[, 1])
    return()
  install.packages(x)
})
#
lapply(s2, function(x) {
  if (x %in% installed.packages()[, 1])
    return()
  BiocManager::install(x, update = FALSE)
})
#
a <- installed.packages()[,1 ]
#
xs <- c(s1, s2)
missingPkg <- setdiff(xs, a)

if (length(missingPkg) > 0)
  stop(paste("this packages are missing", paste(missingPkg, collapse = " ")))