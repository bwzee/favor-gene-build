
library(tidyverse)
library(readxl)
#Functions

#' get_coltypes Return column types for a data frame
#'
#' @param x
#' @param source
#'
#' @return
#' @export
#'
#' @examples
get_coltypes<-function(x,source="SOURCE",url="http://someurl",group="Mendelian"){
  summary.default(x) %>% as.data.frame() %>% as_tibble() %>% filter(Var2=="Mode") %>%
    rename(Column_type=Freq,column_name=Var1) %>% select(-Var2)-> coldesc
  coldesc$SOURCE=source
  coldesc$LINK=url
  coldesc$GROUP=group
  return(coldesc)

}

#' Title
#'
#' @param infile
#'
#' @return
#' @export
#'
#' @examples
parse_s_het<-function(infile="../s-het/075523-1.xlsx"){
  dt=readxl::read_xlsx(infile) %>% janitor::clean_names()
  return(dt)
}

#' Title
#'
#' @param infile
#'
#' @return
#' @export
#'
#' @examples
parse_gtex<-function(infile="../GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct"){
   dat=fread(infile)
   return(dat)
}
##decipher
#' Title
#'
#' @param file
#'
#' @return
#' @export
#'
#' @examples
parse_decipher_ddgp<-function(file="../DECIPHER/DDG2P_9_9_2021.csv"){

  dt=fread(file) %>% janitor::clean_names()
  #
  dt=as_tibble(dt)

  #print(head(dt))
  return(dt)
}



## Mouse MGI Mouse knockout model phenotypes ------
#' parse_mouse_mgi_knockout  Mouse MGI Mouse knockout model phenotypes
#'
#' @param infile
#' @param human_mapping_file
#'
#' @return
#' @export
#'
#' @examples
parse_mouse_mgi_knockout<-function(infile="http://www.informatics.jax.org/downloads/reports/MGI_PhenotypicAllele.rpt",  flatten=TRUE,
                                    human_mapping_file="http://www.informatics.jax.org/downloads/reports/HMD_HumanPhenotype.rpt"){
  log_info(glue("Reading {infile}"))
    d=fread(infile,sep="\t",skip = 7,quote="")
  setnames(d,c(
    "MGI Allele Accession ID",
    "Allele Symbol",
    "Allele Name","Allele Type","Allele Attribute","PubMed ID for original reference",
    "MGI Marker Accession ID","Marker Symbol",
     "Marker RefSeq ID","Marker Ensembl ID","
    High-level Mammalian Phenotype ID","Synonyms",
  "Marker Name"))
  d %>% janitor::clean_names()->d

  hmap=fread(human_mapping_file,sep="\t")
  hmap$V6=NULL
  setnames(hmap, c("gene_symbol","human_entrez_id","mouse_marker_symbol","mgi_marker_accession_id","mammalian_phenotype"))
  #merge Human symbols
  d<-left_join(as_tibble(d),as_tibble(hmap))

  #Raw data,collapsed
  #return(d)
    if (flatten==T){
      log_info("Collapsing by gene")
      #Collapse by human gene
    d=as.data.table(d)[,by="gene_symbol",list(
        mgi_allele_accession_id=paste(unique(mgi_allele_accession_id),collapse=","),
      allele_symbol=paste(unique(allele_symbol),collapse=",")
      )]
    print(dim(d))
    return(d)
  }else{
  d=unique(as.data.table(d),by="gene_symbol")
    return(d)
  }
}
