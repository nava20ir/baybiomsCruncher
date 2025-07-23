library(DIAgui)


make_final_pg <- function(dian_tsv_path,mapping_csv_path)
{
df <- as.data.frame(data.table::fread(dian_tsv_path, stringsAsFactors = FALSE))
mapping_csv = read.csv(mapping_csv_path,stringsAsFactors = F)

mapcolname2exp <- function(x,mapping_csv){
  if (x %in% mapping_csv$name) {
    return(mapping_csv$mapping[mapping_csv$name == x])
  }else{
    return(x)
  }
}

window2linux <- function(windows_path){
  relative_path <- sub("^[A-Za-z]:", "", windows_path)  # remove 'E:'
  linux_path = gsub("\\\\", "/", relative_path)
  return(basename(linux_path))        # convert backslashes to slashes 
  
}

colnames(df) = sapply(colnames(df),function(x)window2linux(x))
colnames(df) = sapply(colnames(df),function(x)mapcolname2exp(x,mapping_csv))
df = df[,c(setdiff(colnames(df),mapping_csv$mapping),mapping_csv$mapping)]
return(df)
}



mapping_expr_raw_file_2_name <- function(mapping_csv,dian_sourcec){
  
  mapping_csv$name = gsub(mapping_csv$name,pattern = '.raw',replacement ='')
  mapping_dic <- as.character(mapping_csv$mapping)
  names(mapping_dic) <- as.character(mapping_csv$name)
  dian_sourcec$File.Name <- sapply(dian_sourcec$Run, function(x)mapping_dic[x])
  
return(dian_sourcec)
  
}


baybioms_report_process <- function(data, header.id = "Protein.Group", 
                                    sample.id = "File.Name",
                                    quantity.id = "Precursor.Normalised",
                                    secondary.id = "Precursor.Id",
                                    id_to_add = c("Protein.Names", "First.Protein.Description", "Genes"),
                                    qv = 0.01,
                                    pg.qv = 0.01,
                                    p.qv = 1,
                                    gg.qv = 1,
                                    quality = 0.8,
                                    only_proteotypic = FALSE,
                                    get_pep = TRUE,
                                    only_pepall = FALSE,
                                    get_Top3 = FALSE,
                                    get_iBAQ = TRUE,
                                    fasta = NULL,
                                    species = NULL,
                                    peptide_length = c(5,36),
                                    mapping_dic = NULL,
                                    format = c("xlsx", "csv", "txt")){

  
  if (is.null(mapping_dic)) stop('no mapping data') 
  if (!all(c('name','mapping') %in% colnames(mapping_dic))) stop('Two columns: raw and mapping should exist in mapping file')
  report_raw <- diann_load(data)  #load your report file
  
   tryCatch({
     report <- mapping_expr_raw_file_2_name(mapping_dic,report_raw)
   }, error = function(w) {
     stop(conditionMessage(w))
   })

  if(any(is.na(report$File.Name ))) stop('something went wrong with mapping experiment to raw files pleases check the mapping csv')
  
  brut <- report
  if(header.id == "Protein.Group" | header.id == "Genes"){
    if(only_proteotypic){
      report <- report[which(report[["Proteotypic"]] != 0), ]
    }
  }

  report <- report %>% dplyr::filter(Q.Value <= qv & PG.Q.Value <= pg.qv & Protein.Q.Value <= p.qv & GG.Q.Value <= gg.qv)
  if("Quantity.Quality" %in% colnames(report)){
    report <- report[which(report[["Quantity.Quality"]] >= quality),]
  }
    
  #preprocess the data in order to use MaxLFQ (see documentation from iq : https://cran.r-project.org/web/packages/iq/index.html  --> see the vignettes)
  iq_report <- iq::preprocess(report,intensity_col = quantity.id,
                                primary_id = header.id,
                                sample_id  = sample.id,
                                secondary_id = secondary.id,
                                median_normalization = FALSE,
                                pdf_out = file.path(NULL))
   
    iq_report <- iq::fast_MaxLFQ(iq_report) #return a list, check it, call element with '$'
    iq_report <- iq_report$estimate #estimate is the dataset
    iq_report <- as.data.frame(iq_report)
    iq_report = 2^ iq_report # anti log of intensity values
    iq_report <- iq_report[order(rownames(iq_report)),]

    
  ### Add gene names and other informations; reshape data frame
  nc <- ncol(iq_report)
  iq_report[[header.id]] <- rownames(iq_report)
  rownames(iq_report) <- 1:nrow(iq_report)
  report <- report[(report[[header.id]] %in% iq_report[[header.id]]),]
  report <- report[order(report[[header.id]]),]
  col_n <- colnames(report)
  for(i in id_to_add){
    if(i %in% col_n){
      iq_report[[i]] <- unique(report[,c(header.id, i)])[[i]]
    }
  }
  
  iq_report <- iq_report[,c((nc+1):ncol(iq_report), 1:nc)]  # reorder columns
  d_seq <- getallseq(pr_id = iq_report$Protein.Group,
                       fasta_file = TRUE,
                       bank_name = fasta)

    n_cond <- length(unique(report$File.Name))
    n_info <- length((nc+1):ncol(iq_report))
    brut <- diann_matrix(brut, sample.header = sample.id,
                         id.header = "Protein.Group",
                         quantity.header = "Precursor.Quantity",
                         proteotypic.only = TRUE,
                         q = qv, protein.q = p.qv,
                         pg.q = pg.qv, gg.q = gg.qv,
                         method = "sum")
    
    
    brut <- get_iBAQ(brut, proteinDB = d_seq,
                     id_name = "Protein.Group",
                     #ecol = n_info:(n_cond+1),
                     ecol = which(colnames(brut) %in% mapping_dic$mapping),
                     peptideLength = peptide_length,
                     proteaseRegExp = getProtease("trypsin"),keep_original = FALSE,
                     log2_transformed = FALSE)
    
  colnames(brut) <- gsub(colnames(brut),pattern='iBAQ_',replacement='')
  brut <- brut[,c(setdiff(colnames(brut),mapping_dic$mapping),mapping_dic$mapping)]
  iq_report <- iq_report[,c(setdiff(colnames(iq_report),mapping_dic$mapping),mapping_dic$mapping)]
  return(list('ibaq'=brut,'max_lfq'=iq_report))
}







#' Convert a FASTA File to UTF-8 Encoding
#'
#' This function reads a FASTA file that may be encoded in UTF-16LE and converts it
#' to UTF-8 encoding. It automatically detects the file encoding and ensures that
#' the output file is saved in UTF-8, line by line.
#'
#' @param input_path Character. Path to the input FASTA file. The function will guess the encoding.
#' @param output_path Character. Path to the output file that will be written in UTF-8 encoding.
#'
#' @return No return value. The function performs file conversion as a side effect.
#'
#' @details
#' The function uses `readr::guess_encoding()` to detect whether the file is encoded
#' in UTF-16LE. If so, it reads the file using the appropriate encoding. Otherwise, it assumes
#' the file is in the system default encoding. The converted file is written using UTF-8 encoding.
#'
#' @examples
#' \dontrun{
#' convert_fasta_file("input.fasta", "output_utf8.fasta")
#' }
#'
#' @importFrom readr guess_encoding
#' @export
#' 
convert_fasta_file <- function(input_path,output_path){
if (any(grepl(x=readr::guess_encoding(input_path),pattern = '16LE'))){
  print('16LE')
  input <- file(input_path, open = "r", encoding = "UTF-16LE")  
}else{
  input <- file(input_path, open = "r") 
}
output <- file(output_path, open = "w", encoding = "UTF-8")
while(length(line <- readLines(input, n = 1, warn = FALSE)) > 0) {
  writeLines(line, output)
}
close(input)
close(output)
}
