#mergeAll
library(data.table)
library(logger)
library(DT)
source("parse_omim.R")
source("data_readers.R")


outjson="../output/genedb.v2.out.json"
tsvout="../output/genedb.v2.out.tsv"
outrds="../output/genedb.v2.rds"
coldesc="../output/genedb.column_descriptions.v2.xlsx"


log_info("Reading in all data")
dir.create(dirname(outjson),showWarnings = F)

dbnsfp_in="../dbNSFP4.1_gene.complete"
hpa_in="../proteinatlas.tsv"
hgnc_in="../gene_with_protein_product.txt"

fread(dbnsfp_in) %>% janitor::clean_names()->dbnsfp

data.frame(source="DBNSFP_Gene",column=names(dbnsfp)) ->dsn_names


#HPA
log_info("Reading HPA")
fread(hpa_in) %>% janitor::clean_names()-> hpa
as_tibble(hpa) %>% select(-grep("_nx|pathology_prog",names(hpa),value=T)) ->hpa
hpa<-as.data.table(hpa)

rbind(data.frame(source="HPA_protein_atlas",column=names(hpa)), dsn_names)-> dsn_names


log_info("Reading in HGNC")
fread(hgnc_in) %>% janitor::clean_names()->hgnc
rbind(data.frame(source="HGNC_HUGO",column=names(hgnc)), dsn_names)-> dsn_names


log_info("Reading OMIM")
parse_omim()->omim
rbind(data.frame(source="OMIM",column=names(omim)), dsn_names)-> dsn_names


print(dim(dbnsfp))
print(dim(hgnc))
log_info("Commencing merge of DBNSFP + HGNC")
#gene_name
merge(dbnsfp,hgnc,by.x="ensembl_gene",by.y="ensembl_gene_id",all.x=T)->mg

log_info("Adding OMIM")
merge(mg,omim,by.x="ensembl_gene",by.y="ensembl_gene_id",all.x=T)->mg

log_info("Adding S_HET (gene name merge)")
shet<-parse_s_het()
merge(mg,shet,by.x="gene_symbol",by.y="gene_symbol",all.x=T)->mg
rbind(data.frame(source="S-HET",column=names(shet)), dsn_names)-> dsn_names



log_info("Reading DECIPHER DD2GP : Gene-phenotyppe")
decipher=parse_decipher_ddgp()
rbind(data.frame(source="DECIPER_DD2GP",column=names(decipher)), dsn_names)-> dsn_names
merge(mg,decipher,by.x="gene_symbol",by.y="gene_symbol",all.x=T)->mg





log_info("Adding HPA")
merge(mg,hpa,by.x="ensembl_gene",by.y="ensembl",all.x=T)-> mg2

setnames(mg2,"gene_name.x","gene_name")

xnames= grep("\\.x$",names(mg2),value=T)

mg2 %>% select(-xnames) ->mg2



log_info("Saving RDS")
saveRDS(mg2,outrds)

log_info("Writing JSON output")
write(jsonlite::toJSON(mg2,pretty = T),outjson)
write.table(mg2,tsvout,row.names=F,quote=F,sep="\t")

dsn_names=as.data.table(dsn_names)
dsn_names$column_type="character"
dsn_names[grepl("score|tpm|fpkm|exac_|gnomad_|rvis_|position|location",column),column_type:="numeric"]


writexl::write_xlsx(dsn_names,coldesc)

write.table(dsn_names,"../output/genedb.column_descriptions.tsv",sep="\t",quote=F,row.names = F)
