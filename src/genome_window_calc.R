#Window test

library(reshape) # to rename columns
library(data.table) # to make sliding window dataframe
library(zoo) # to apply rolling function for sliding window
library(ggplot2)

#upload data to dataframe, rename headers, make locus continuous, create subsets
depth <- read.table("sorted.depth", sep="\t", header=F)

depth<-rename(depth,c(V1="Chr", V2="locus", V3="coverageX", V3="coverageY"))
depth$locus <- 1:12157105
Xdepth<-subset(depth, select = c("Chr", "locus","coverageX"))
#genome coverage as sliding window
Xdepth.average<-setDT(Xdepth)[, .(
                window.start = rollapply(locus, width=1000, by=1000, FUN=min, align="left", partial=TRUE),
                window.end = rollapply(locus, width=1000, by=1000, FUN=max, align="left", partial=TRUE),
                coverage = rollapply(coverage, width=1000, by=1000, FUN=mean, align="left", partial=TRUE)
              ), .(Chr)]
              

## tdyvers

as_tibble(gwascat) %>% group_by(chr_id)


as_tibble(gwascat) %>% select(chr_id,chr_pos,p_value)  %>% group_by(chr_id)  %>% do(
  data.frame(
    window.start = rollapply(.$chr_pos, width=wd, by=wd, FUN=min, align="left"),
    window.end = rollapply(.$chr_pos, width=wd, by=wd, FUN=max, align="left"),
    coverage = rollapply(.$chr_pos, width=wd, by=wd, FUN=mean, align="left")
  )
)



setDT(gwascat)
wd=10000
gwascat[,by="chr_id",.(
  window.start = zoo::rollapply(chr_pos, width=wd,by=wd, FUN=min, align="left",partial=T),
  window.end = zoo::rollapply(chr_pos, width=wd,by=wd,  FUN=max, align="left",partial=T),
  coverage = zoo::rollapply(p_value, width=wd,by=wd, FUN=mean, align="left",partial=T)
)]->rr
