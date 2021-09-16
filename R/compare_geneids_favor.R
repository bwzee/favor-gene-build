#Compare Favor HGNC symbols

compare_geneids_favor<-function(){
setwd("..")

favorids=fread("FAVOR/geneIdDistinct.txt",header=F)$V1
favor_genes=fread("FAVOR/geneNameDistinct.txt",header=F)$V1


dbdat=fread("output/genedb.v1.out.csv",fill=T)

dbgenes=unique(dbdat$gene_name)
dbens=as_tibble(dbdat) %>% filter(ensembl_gene!=".") %>% pull(ensembl_gene) %>% unique()

#Get Ensembl Common annotations (hg38)


annodata=readRDS("ensembl_anno_hg38.rds") %>% janitor::clean_names()
prot_genes=annodata[gene_type=="protein_coding"]
unique(prot_genes,by="gene_stable_id") ->prot_genes



print("Common")
intersect(dbgenes,favor_genes) %>% length()->a
setdiff(dbgenes,favor_genes) %>% length()->b
setdiff(favor_genes,dbgenes) %>% length()->c
setdiff(favor_genes,dbgenes)-> fv_hugo
setdiff(fv_hugo,prot_genes$hgnc_symbol) %>% length() ->d

print("Ensembl Genes")
intersect(dbens,favorids) %>% length() ->ea
setdiff(dbens,favorids) %>% length()->eb
setdiff(favorids,dbens) %>% length()->ec
setdiff(favorids,dbens) -> fv_ens
setdiff(fv_ens,prot_genes$gene_stable_id) %>% length()->ed

mdf=data.frame(Hugo=c(a,b,c,d),Ensembl_Ids=c(ea,eb,ec,ed))

rownames(mdf)=c("Common","Unique_to_DB","Unique_to_FAVOR","FAVOR_not_protein_coding_ensembl")
datatable(mdf)

#How many common genes are favor with Ensembl: 531
intersect(fv_ens,prot_genes$gene_stable_id) %>% length()
#How many favor ensembl IDs are not in protein-coding ensembl: ] 20623
setdiff(fv_ens,prot_genes$gene_stable_id) %>% length()
}

