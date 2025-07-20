
library(DIAgui)
#' report_process
#'
#' Process data from DIAnn in a single workflow, save it and also save some quality check plots.
#'
#' @param data Path to data (report from DIAnn)
#' @param header.id Id column name (protein, gene, ...)
#' @param sample.id Id column that contains the fractions names
#' @param quantity.id Quantity column name (Precursor.Normalised)
#' @param secondary.id Id column that contains precursor (usefull when get_pep is TRUE)
#' @param id_to_add Some secondary id you want to keep in your final data
#' @param qv Precursor q-value threshold
#' @param pg.qv Protein group q-value threshold
#' @param p.qv Uniquely identified protein q-value threshold
#' @param gg.qv Gene group q-value threshold
#' @param quality Quantity quality threshold
#' @param only_proteotypic If TRUE, will only keep proteotypic
#' @param get_pep logical; get peptide count ?
#' @param only_pepall logical; should only keep peptide counts all or also peptide counts for each fractions ?
#' @param get_Top3 If TRUE, will calculate Top3 quantification
#'                 (only use when \code{header.id} is set to 'Protein.group' or 'Genes)
#' @param get_iBAQ If TRUE, will calculate iBAQ quantification
#'                 (only use when \code{header.id} is set to 'Protein.group' )
#' @param fasta When iBAQ is set to TRUE, you can add path to fasta files.
#'              If NULL, will get theoretical peptides from swissprot bank (quite long).
#' @param species When fasta is NULL, you need to specify the specie from your data
#' @param peptide_length The minimum and maximum peptide length
#' @param format Format from your saved file (only xlsx, csv and txt are supported)
#'
#' @importFrom utils write.csv write.table
#'
#' @export

baybioms_report_process <- function(data, header.id = "Protein.Group", sample.id = "File.Name",
                           quantity.id = "Precursor.Normalised", secondary.id = "Precursor.Id",
                           id_to_add = c("Protein.Names", "First.Protein.Description", "Genes"),
                           qv = 0.01, pg.qv = 0.01, p.qv = 1, gg.qv = 1, quality = 0.8, only_proteotypic = FALSE,
                           get_pep = TRUE, only_pepall = FALSE, get_Top3 = FALSE, get_iBAQ = TRUE,
                           fasta = NULL, species = NULL, peptide_length = c(5,36),
                           format = c("xlsx", "csv", "txt")){
  
  if(inherits(data, "character")){
    if(file.exists(data)){
      extension <- stringr::str_split(data, "\\.")[[1]]
      if(length(extension) == 1){
        stop("Don't forget to type the format of your file !")
      }
      
      dname <- extension[length(extension)-1]
      dname <- stringr::str_split(dname, "/")[[1]]
      dname <- dname[length(dname)]
      dname <- paste(dname, header.id, sep="_")
      dname <- paste0(format(Sys.time(), "%y%m%d_%H%M_"), dname)
      
      extension <- extension[length(extension)]
    }
    else{
      stop("This file doesn't exists ! Check your file name")
    }
  }
  else {
    stop("Please, provide a file name")
  }
  
  
  if(!file.exists(dname)){
    dir.create(dname)
  }
  
  report <- diann_load(data)  #load your report file
  brut <- report
  
  if(header.id == "Protein.Group" | header.id == "Genes"){
    if(only_proteotypic){
      report <- report[which(report[["Proteotypic"]] != 0), ]
    }
    ### Filtering on q values
    report <- report %>% dplyr::filter(Q.Value <= qv & PG.Q.Value <= pg.qv &
                                         Protein.Q.Value <= p.qv & GG.Q.Value <= gg.qv)
    
    ### Filtering on quality
    if("Quantity.Quality" %in% colnames(report)){
      report <- report[which(report[["Quantity.Quality"]] >= quality),]
    }
    
    ### iq process
    
    #preprocess the data in order to use MaxLFQ (see documentation from iq : https://cran.r-project.org/web/packages/iq/index.html  --> see the vignettes)
    iq_report <- iq::preprocess(report,intensity_col = quantity.id,
                                primary_id = header.id,
                                sample_id  = sample.id,
                                secondary_id = secondary.id,
                                median_normalization = FALSE,
                                pdf_out = file.path(dname, NULL)
    )
  }
  else{
    iq_report <- diann_matrix(report, sample.header = sample.id,
                              id.header = header.id,
                              proteotypic.only = only_proteotypic,
                              q = qv, quality = quality,
                              protein.q = p.qv,
                              pg.q = pg.qv,
                              gg.q = gg.qv,
                              method = "max")
  }
  
  
  if(get_Top3 & (header.id == "Protein.Group" | header.id == "Genes")){
    Top3 <- iq::preprocess(report,
                           intensity_col = "Precursor.Quantity",
                           primary_id = header.id,
                           sample_id  = sample.id,
                           secondary_id = secondary.id,
                           median_normalization = FALSE,
                           pdf_out = NULL)
    Top3 <- Top3 %>% dplyr::group_by(protein_list) %>%
      tidyr::spread(sample_list, quant)
    Top3$id <- NULL
    top3_f <- function(x){
      if(sum(!is.na(x)) < 3){
        x <- NA
      }
      else{
        x <- x[order(x, decreasing = TRUE)][1:3]
        x <- mean(x)
      }
    }
    Top3 <- Top3 %>% dplyr::group_by(protein_list) %>%
      dplyr::summarise(dplyr::across(dplyr::everything(), top3_f))
    
    Top3 <- as.data.frame(Top3)
    rownames(Top3) <- Top3$protein_list
    Top3$protein_list <- NULL
    colnames(Top3) <- paste0("Top3_", colnames(Top3))
    Top3 <- Top3[order(rownames(Top3)),]
  }
  
  if(get_pep & (header.id == "Protein.Group" | header.id == "Genes")){
    ### Get peptides counts from preprocess
    pc <- iq_report %>% dplyr::group_by(protein_list, sample_list) %>%
      dplyr::mutate("countpep" = length(unique(id)))
    pc <- unique(pc[,c("protein_list", "sample_list", "countpep")])
    pc <- tidyr::spread(pc, sample_list, countpep)
    pc[is.na(pc)] <- 0
    pc <- as.data.frame(pc)
    rownames(pc) <- pc$protein_list
    pc$protein_list <- NULL
    pc <- pc[order(rownames(pc)),]
    colnames(pc) <- paste0("pep_count_", colnames(pc))
    pc$peptides_counts_all <- unname(apply(pc, 1, max))
    pc <- pc[,c(ncol(pc), 1:(ncol(pc)-1))]
  }
  
  if(header.id == "Protein.Group" | header.id == "Genes"){
    ### Get the protein group
    ## Perform fast MaxLFQ function from iq
    iq_report <- iq::fast_MaxLFQ(iq_report) #return a list, check it, call element with '$'
    
    iq_report <- iq_report$estimate #estimate is the dataset
    #other element is annotation, depends on your data
    # quick handle of data  --> add column protein group and put it at the beginning, remove rownames
    iq_report <- as.data.frame(iq_report)
    iq_report = 2^ iq_report # anti log of intensity values

    iq_report <- iq_report[order(rownames(iq_report)),]
    
    if(get_Top3){
      iq_report <- cbind(iq_report, Top3)
    }
    ## Add peptide counts
    if(get_pep & only_pepall){ #only keep peptide counts all
      iq_report$peptides_counts_all <- pc$peptides_counts_all
    }
    else if(get_pep){
      iq_report <- cbind(iq_report, pc)
    }
  }
  
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
    else{
      message(paste(i, "is not in your colnames data. Check the id you want to add."))
    }
  }
  iq_report <- iq_report[,c((nc+1):ncol(iq_report), 1:nc)]  # reorder columns
  #max_lfq <<- iq_report
  
  if(get_iBAQ & header.id == "Protein.Group"){
    if(!is.null(fasta)){
      d_seq <- getallseq(pr_id = iq_report$Protein.Group,
                         fasta_file = TRUE,
                         bank_name = fasta)
    }
    else{
      d_seq <- getallseq(pr_id = iq_report$Protein.Group,
                         spec = species)
    }
    n_cond <- length(unique(report$File.Name))
    n_info <- length((nc+1):ncol(iq_report))
    brut <- diann_matrix(brut, sample.header = sample.id,
                         id.header = "Protein.Group",
                         quantity.header = "Precursor.Quantity",
                         proteotypic.only = only_proteotypic,
                         q = qv, protein.q = p.qv,
                         pg.q = pg.qv, gg.q = gg.qv,
                         method = "sum")
    brut$Genes <- NULL
    brut$Protein.Names <- NULL
    brut <- get_iBAQ(brut, proteinDB = d_seq,
                     id_name = "Protein.Group",
                     ecol = n_info:(n_cond+1),
                     peptideLength = peptide_length,
                     proteaseRegExp = getProtease("trypsin"),
                     log2_transformed = FALSE)
    brut <- brut[,-c(n_info:(n_cond+1))]
    
    #iq_report <- dplyr::left_join(iq_report, brut, by = "Protein.Group")
  }
  
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
