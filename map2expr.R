mapping_df = read.csv('/Projects/002_Proteomics/P376_Lorenzo_and_Claassens_DnaAcos/P376_02/amir_dian_test/mapping.csv')
ibaq_df = read.csv('/Projects/002_Proteomics/P376_Lorenzo_and_Claassens_DnaAcos/P376_02/amir_dian_test/ibaq.tsv',check.names = F,sep = '\t')
maxlfq_df = read.csv('/Projects/002_Proteomics/P376_Lorenzo_and_Claassens_DnaAcos/P376_02/amir_dian_test/maxlfq.tsv',check.names = F,sep = '\t')
mapping_df = mapping_df[,c(2,3)]

colnames(mapping_df) = c('raw','exp')
make_dictionary <- function(mapping_df,intensity_df){
  mapping_df$mapping_col = apply(mapping_df,1,function(x)colnames(intensity_df)[grep(x= colnames(intensity_df),pattern = x[1])])
return(mapping_df)
}

map2experiment <- function(mapping_df,intensity_df){
  for (i in 1:ncol(intensity_df)){
    col <- colnames(intensity_df)[i]
    if (col %in% mapping_df$mapping_col) colnames(intensity_df)[i] = mapping_df$exp[mapping_df$mapping_col == col]
  }
  return(intensity_df)
}

mapibaq = make_dictionary(mapping_df,ibaq_df)
final_ibaq = map2experiment(mapibaq,ibaq_df)

maxlfq_df = maxlfq_df[,!grepl(colnames(maxlfq_df),pattern='_count_')]
maplfq= make_dictionary(mapping_df,maxlfq_df)
final_maxlfq = map2experiment(maplfq,maxlfq_df)
# TODO iBAQ has lesss columns needs to be merged