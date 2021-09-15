
library(tidyverse)
library(readxl)
#Functions

get_coltypes<-function(x,source="SOURCE"){
  summary.default(x) %>% as.data.frame() %>% as_tibble() %>% filter(Var2=="Mode") %>%
    rename(Column_type=Freq,column_name=Var1) %>% select(-Var2)-> coldesc
  coldesc$SOURCE=source
  return(coldesc)
  
}

parse_s_het<-function(infile="../s-het/075523-1.xlsx"){
  dt=readxl::read_xlsx(infile) %>% janitor::clean_names()
  return(dt)
}

parse_gtex<-function(infile="../GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct"){ 
  
}
##decipher
parse_decipher_ddgp<-function(file="../DECIPHER/DDG2P_9_9_2021.csv"){
  
  dt=fread(file) %>% janitor::clean_names()
  
  dt=as_tibble(dt)
  
  print(head(dt))
  return(dt)
}

