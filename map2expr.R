# the ideas is to map raw files to their experiment names in the DIA-NN output

#' Read and Format Mapping File for Expression Data
#'
#' This function reads a CSV mapping file and returns a data frame with
#' relevant columns. If the input file has more than two columns,
#' only the second and third columns are retained (assumed to contain
#' identifiers for mapping).
#'
#' @param mapping_file Character. Path to the mapping CSV file.
#'
#' @return A `data.frame` containing two columns for mapping, or `NULL` if an error occurs during reading.
#'
#' @details
#' The file is read using `read.csv()` with `UTF-8` encoding and `stringsAsFactors = FALSE`.
#' If the input file has more than two columns, only the second and third columns are returned.
#' This is useful when the first column contains row numbers or irrelevant metadata.
#'
#' @examples
#' \dontrun{
#' df <- read_mapping_raw2expr("mapping.csv")
#' head(df)
#' }
#'
#' @export
read_mapping_raw2expr <- function(mapping_file) {
  tryCatch({
    df <- read.csv(mapping_file, stringsAsFactors = FALSE, encoding = "UTF-8")
    if (ncol(df)>2){
    df <- df[,c(2,3)]
    }
    return(df)
  },
  error = function(e) {
    message("❌ Error reading mapping file: ", conditionMessage(e))
    return(NULL)
  })
}


#' Generate Column Mapping Dictionary
#'
#' This function creates a mapping between `raw` file identifiers and experiment labels
#' by matching entries in the `mapping_df` with column names in an intensity data frame.
#'
#' @param mapping_df A `data.frame` with two columns: the first column should contain
#' substrings to match column names in `intensity_df`, and the second should contain experiment labels.
#' @param intensity_df A `data.frame` or `matrix` with column names to match against.
#'
#' @return A `data.frame` with three columns: `raw`, `exp`, and `mapping_col`, where
#' `mapping_col` contains the matched column names in `intensity_df`.
#'
#' @details
#' The function assumes the first column of `mapping_df` contains substrings that can
#' uniquely identify columns in `intensity_df`. The matching is done using `grep`.
#' If no match is found, the corresponding entry in `mapping_col` will be `NA`.
#'
#' @examples
#' \dontrun{
#' mapping_df <- data.frame(raw = c("Sample1", "Sample2"), exp = c("Treated", "Control"))
#' intensity_df <- data.frame(Sample1_A = 1:3, Sample2_B = 4:6, check.names = FALSE)
#' make_dictionary(mapping_df, intensity_df)
#' }
#'
#' @export
make_dictionary <- function(mapping_df,intensity_df){
  colnames(mapping_df) = c('raw','exp')
  mapping_df$mapping_col = apply(mapping_df,1,function(x)colnames(intensity_df)[grep(x= colnames(intensity_df),pattern = x[1])])
return(mapping_df)
}



#' Map Intensity Columns to Experiment Names
#'
#' This function renames columns in an intensity data frame based on a mapping
#' data frame that links raw identifiers to experiment names.
#'
#' @param mapping_df A data frame with at least two columns:
#'   \describe{
#'     \item{`mapping_col`}{Column names in `intensity_df` to match}
#'     \item{`exp`}{Experiment names to assign to matched columns}
#'   }
#' @param intensity_df A data frame whose column names are to be renamed
#' according to the mapping in `mapping_df`.
#'
#' @return A data frame identical to `intensity_df`, but with some columns
#' renamed based on the mapping.
#'
#' @details
#' Columns in `intensity_df` that match `mapping_df$mapping_col` will be renamed
#' using the corresponding `mapping_df$exp` value. If multiple matches are found,
#' only the first match is used and a warning is issued.
#'
#' @examples
#' \dontrun{
#' mapping_df <- data.frame(
#'   raw = c("Raw1", "Raw2"),
#'   exp = c("SampleA", "SampleB"),
#'   mapping_col = c("Intensity_Raw1", "Intensity_Raw2")
#' )
#' intensity_df <- data.frame(
#'   Intensity_Raw1 = 1:3,
#'   Intensity_Raw2 = 4:6
#' )
#' map2experiment(mapping_df, intensity_df)
#' }
#'
#' @export
map2experiment <- function(mapping_df,intensity_df){
  for (i in 1:ncol(intensity_df)){
    col <- colnames(intensity_df)[i]
    if (col %in% mapping_df$mapping_col) colnames(intensity_df)[i] = mapping_df$exp[mapping_df$mapping_col == col]
  }
  return(intensity_df)
}



#' Map iBAQ Columns to Experiment Labels
#'
#' This function maps raw file identifiers in an iBAQ intensity data frame
#' to user-defined experiment labels based on a mapping data frame.
#'
#' @param mapping_df A data frame with at least two columns:
#'   \describe{
#'     \item{`raw`}{Substring to match in column names of `ibaq_df`}
#'     \item{`exp`}{Label to assign to matching columns}
#'   }
#' @param ibaq_df A data frame containing iBAQ intensity values with column names
#' corresponding to raw identifiers.
#'
#' @return A data frame with renamed columns (if matching is successful).
#' If an error occurs, returns the original `ibaq_df`.
#'
#' @examples
#' \dontrun{
#' mapping_df <- data.frame(raw = c("RAW1", "RAW2"), exp = c("Sample1", "Sample2"))
#' ibaq_df <- data.frame(RAW1_int = 1:3, RAW2_int = 4:6)
#' make_final_ibaq(mapping_df, ibaq_df)
#' }
#'
#' @export
make_final_ibaq <- function(mapping_df, ibaq_df) {
  tryCatch({
    mapibaq <- make_dictionary(mapping_df, ibaq_df)
    final_ibaq <- map2experiment(mapibaq, ibaq_df)
    return(final_ibaq)
  }, error = function(e) {
    message("Error in mapping experiment to raw files: ", conditionMessage(e))
    return(ibaq_df)
  })
}

#' Map maxLFQ Columns to Experiment Labels
#'
#' This function processes a maxLFQ intensity data frame by removing count-related columns
#' and renaming intensity columns based on a mapping table of raw file identifiers to experiment names.
#'
#' @param mapping_df A data frame with columns:
#'   \describe{
#'     \item{`raw`}{Substrings to match in `maxlfq_df` column names}
#'     \item{`exp`}{Target experiment names to assign}
#'   }
#' @param maxlfq_df A data frame with maxLFQ intensity values. Column names are expected to include raw identifiers.
#'
#' @return A `data.frame` with renamed columns where applicable. If an error occurs, the original `maxlfq_df` is returned.
#'
#' @details
#' Columns containing the pattern `"_count_"` are removed prior to mapping. Then, column names are matched
#' to the `raw` identifiers and replaced with the corresponding `exp` values.
#'
#' @examples
#' \dontrun{
#' mapping_df <- data.frame(raw = c("RAW1", "RAW2"), exp = c("Sample1", "Sample2"))
#' maxlfq_df <- data.frame(RAW1_A = 1:3, RAW2_B = 4:6, RAW1_count_A = 0:2)
#' make_final_maxlfq(mapping_df, maxlfq_df)
#' }
#'
#' @export
make_final_maxlfq <- function(mapping_df, maxlfq_df) {
  tryCatch({
    # Remove columns matching '_count_'
    maxlfq_df <- maxlfq_df[, !grepl("_count_", colnames(maxlfq_df))]

    # Build mapping dictionary
    maplfq <- make_dictionary(mapping_df, maxlfq_df)

    # Map to experiment
    final_maxlfq <- map2experiment(maplfq, maxlfq_df)

    return(final_maxlfq)
  },
 error = function(e) {
    message("❌ Error in make_final_maxlfq: ", conditionMessage(e))
    return(maxlfq_df)
  })
}


## Test example
# mapping_df = read.csv('/Projects/002_Proteomics/P376_Lorenzo_and_Claassens_DnaAcos/P376_02/amir_dian_test/mapping.csv')
# ibaq_df = read.csv('/Projects/002_Proteomics/P376_Lorenzo_and_Claassens_DnaAcos/P376_02/amir_dian_test/ibaq.tsv',check.names = F,sep = '\t')
# maxlfq_df = read.csv('/Projects/002_Proteomics/P376_Lorenzo_and_Claassens_DnaAcos/P376_02/amir_dian_test/maxlfq.tsv',check.names = F,sep = '\t')
# 
# final_ibaq = make_final_ibaq(mapping_df,ibaq_df)
# final_maxlf = make_final_maxlfq(mapping_df,maxlfq_df)