
#Parse OMIM phenotypes

library(tidyverse)
library(data.table)


#OMIM Parser
parse_omim<-function(omimfile="../omim/genemap2.txt"){

fread(omimfile,skip = 3,sep="\t",fill=T) %>% janitor::clean_names() %>% as_tibble() ->omim



omim %>% rowwise() %>% mutate(pheno=get_pheno(phenotypes)$phenotype,
                              pheno_mimid=get_pheno(phenotypes)$mimid,
                              pheno_key=get_pheno(phenotypes)$phenokey,
                              inheritance=get_pheno(phenotypes)$inheritance,
                              pheno_disease=get_pheno(phenotypes)$disease
                              )->omim
omim %>% separate(gene_symbols,c("gene_symbol","gene_symbol2"),sep=",",extra="drop",fill="right") ->omim
saveRDS(omim,"omim.rds")

return(omim)
}




# Gwas cat . < 1e-8 -> high confidence----
parse_gwascat<-function(infile="../GWAS_catalog/gwas_catalog_v1.0.2-associations_e100_r2021-03-25.tsv"){
  dat=fread(infile,quote="")
  janitor::clean_names(dat)->dat
  as.numeric(dat$p_value) ->dat$p_value
  dat$confidence="low"
  dat[p_value <=1e-8,confidence:="high"]
  
  return(dat)
  
}


#Functionn to extract phenotype from omim string
get_pheno<-function(x){
  # {Alkaline phosphatase, plasma level of, QTL 2}, 612367 (2) 
  disease=TRUE
  if (grepl("\\[",x))
    disease=FALSE
  
  x=gsub(";.+","",x) #Take first entry
  
  m1=str_match(pattern="(.+),\\s*(\\d+)\\s+\\((\\d+)\\),*\\s*(.*)",string=x)
  
  phenotype=m1[2]
  
  mimid=m1[3]
  
  phenokey=m1[4]
  
  inheritance=m1[5]
  
  res=list(
    phenotype=phenotype,
    mimid=mimid,
    phenokey=phenokey ,
    inheritance=inheritance,
    disease=disease
  )
  return(res)
}

