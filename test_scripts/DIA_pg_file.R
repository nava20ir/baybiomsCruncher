
library(DIAgui)

dian_tsv_path  = '/Projects/002_Proteomics/P376_Lorenzo_and_Claassens_DnaAcos/P376_02/amir_dian_test/results.pg_matrix.tsv'
mapping_csv_path = '/Projects/002_Proteomics/P376_Lorenzo_and_Claassens_DnaAcos/P376_02/amir_dian_test/mapping.csv'

make_final_pg <- function(dian_tsv_path,mapping_csv_path){
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

final_pg <- make_final_pg(dian_tsv_path,mapping_csv_path)