

#Make Excel with headers


library(readr)
library(data.table)
library(tidyverse)
library(writexl)
library(DT)

makeHeaderSheet<-function(infile="../dbNSFP4.1_gene.complete.header.tsv"){
  fread(infile) %>% janitor::clean_names()->dbnsfp
  data.frame(colname=names(dbnsfp)) -> ot
  ot$data=unlist(dbnsfp[1,])
  return(ot)
}




dbnsfp=makeHeaderSheet(infile="../dbNSFP4.1_gene.complete.header.tsv")
hpa=makeHeaderSheet(infile="../proteinatlas.header.tsv")
hugo=makeHeaderSheet(infile="../gene_with_protein_product.header.tsv")


write_xlsx(list(
  dbnsfp=dbnsfp,
  HPA=hpa,
  HGNC=hugo
),"../DataFiles_schema.xlsx")