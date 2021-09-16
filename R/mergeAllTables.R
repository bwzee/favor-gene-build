#mergeAll
library(data.table)
library(logger)
library(DT)
library(glue)
library(tidyverse)
library(readxl)
library(writexl)



#' build_main_genedb Main workflow for building FAVOR gene DB
#'
#' @param version
#' @param outdir
#' @param dbnsfp_in
#' @param hpa_in
#' @param hgnc_in
#'
#' @return
#' @export
#'
#' @examples
build_main_genedb<-function(version="v2.2",samplesize=13000,testing=FALSE,
                            outdir="../output/",#Define some input files
                            dbnsfp_in="../dbNSFP4.1_gene.complete",
                            omimfile="../omim/genemap2.txt",
                            shetfile="../s-het/075523-2.xlsx",
                            decipherfile="../DECIPHER/DDG2P_9_9_2021.csv",
                            hpa_in="../proteinatlas.tsv",
                            hgnc_in="../gene_with_protein_product.txt"){
#Output file names
log_info(glue("Setting output files for version {version}"))

if(!dir.exists(outdir))
  dir.create(glue("{outdir}"))

outjson=glue("{outdir}/genedb.{version}.out.json")
tsvout=glue("{outdir}/genedb.{version}.out.tsv")
outrds=glue("{outdir}/genedb.{version}.rds")
coldesc=glue("{outdir}/genedb.column_descriptions.{version}.xlsx")
coldesctsv=glue("{outdir}/genedb.column_descriptions.{version}.tsv")

log_info("Reading in all data")

log_info("Begin Data Reading ..")
#Start with reading DBNSFP
log_info("Reading DBNSFP")
fread(dbnsfp_in) %>% janitor::clean_names()->dbnsfp
if (testing==TRUE){
  dbnsfp<-head(dbnsfp,samplesize)
}

#HPA
log_info("Reading HPA")
fread(hpa_in) %>% janitor::clean_names()-> hpa
as_tibble(hpa) %>% select(-grep("_nx|pathology_prog",names(hpa),value=T)) ->hpa
hpa<-as.data.table(hpa)

log_info("Reading in HGNC")
fread(hgnc_in) %>% janitor::clean_names()->hgnc

log_info("Reading OMIM")
parse_omim(omimfile = omimfile)->omim

log_info("Commencing merge of DBNSFP + HGNC") #gene_name used to merge
merge(dbnsfp,hgnc,by.x="ensembl_gene",by.y="ensembl_gene_id",all.x=T)->mg

log_info("Adding OMIM to merged file")
merge(mg,omim,by.x="ensembl_gene",by.y="ensembl_gene_id",all.x=T)->mg

log_info("Adding S_HET (gene name merge) to merged file")
shet<-parse_s_het()
merge(mg,shet,by.x="gene_symbol",by.y="gene_symbol",all.x=T)->mg

log_info("Reading DECIPHER DD2GP : Gene-phenotyppe")
decipher=parse_decipher_ddgp(file = decipherfile)
#Merge with larger file
merge(mg,decipher,by.x="gene_symbol",by.y="gene_symbol",all.x=T)->mg

log_info("Adding MGI Mouse knockout models (gene_symbol)0")
mgi_mouse=parse_mouse_mgi_knockout()
log_info("..joining tables")
merge(mg,mgi_mouse,by.x="gene_symbol",by.y="gene_symbol",all.x=T)->mg


log_info("Adding HPA")
merge(mg,hpa,by.x="ensembl_gene",by.y="ensembl",all.x=T)-> mg2
setnames(mg2,"gene_name.x","gene_name")
#remove column names that have name.x in them
xnames= grep("\\.x$",names(mg2),value=T)
mg2 %>% select(-xnames) ->mg2


log_info("Saving RDS")
saveRDS(mg2,outrds)

log_info("Writing JSON output")
write(jsonlite::toJSON(mg2,pretty = T),outjson)
log_info("Writing TSV output")
write.table(mg2,tsvout,row.names=F,quote=F,sep="\t")
log_info("Core Database Creation Completed")


### Making Column Assignment Descriptors
log_info("2. Making column definitions file")

log_info("DBNSFP")
#dsn_names is a data frame of the column definitions, their source, group, url(if available)
get_coltypes(dbnsfp,source="DBNSFP_4.1_Gene",url="https://drive.google.com/file/d/1tw5DlroGCX153wV4MqZc8qAMkEP19fWb/view?usp=sharing",group = "CORE_GENE_ANNOTATION")->dsn_names
#HPA
log_info("HPA")
rbind(get_coltypes(hpa,source="HPA_protein_atlas",url="https://www.proteinatlas.org/download/proteinatlas.tsv.zip",group="PROT_EXPRESSION"),dsn_names)-> dsn_names
#HUGO
rbind(get_coltypes(hgnc,source="HGNC_UGO",url="https://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/locus_types/gene_with_protein_product.txt",group="CORE_GENE_ANNOTATION"),dsn_names)-> dsn_names
#OMIM
log_info("OMIM")
rbind(get_coltypes(omim,source="OMIM",url="https://omim.org/static/omim/data/mim2gene.txt",group="OMIM Mendelian Inheritance Disease Phenotypes"),dsn_names)-> dsn_names

log_info("S-HET")
rbind(get_coltypes(shet,source="S-HET",url="075523-1.xlsx",group="POPULATION GENETICS"),dsn_names)-> dsn_names

log_info("DECIPHER DD2GB")
rbind(
get_coltypes(decipher,source="DECIPER_DD2GP",
  url="https://panelapp.genomicsengland.co.uk/panels/484/",group="Human Disease Phenotypes"),
dsn_names)-> dsn_names

log_info("MGI-Mouse")
rbind(get_coltypes(mgi_mouse,source="MGI_MOUSE",
                   url="http://www.informatics.jax.org/downloads/reports/",
                   group="Mouse Knockout Models"),dsn_names)-> dsn_names

dsn_names=as.data.table(dsn_names)
#dsn_names$column_type="character"
#dsn_names[grepl("score|tpm|fpkm|exac_|gnomad_|rvis_|position|location",column),column_type:="numeric"]
writexl::write_xlsx(dsn_names,coldesc)
write.table(dsn_names,coldesctsv,sep="\t",quote=F,row.names = F)
log_info("Created data source column descriptors")

}
