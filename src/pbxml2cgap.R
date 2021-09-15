library(xml2)
library(XML)

library(tidyverse)

ifile="~/Downloads/0454_CGAP_UDN537367.pbxml"


xmlobj <- read_xml(ifile)

ls1 <- as_list(xmlobj) #Converts XML to a nested list


person_attributes=ls1$pedigree$people$person %>% names()


xml_find_all(xmlobj,"//pedigree/people") %>%  map_df(~{
  bind_cols(
    relationship= xml_find_all(.x, ".//person/relationships") %>% xml_text() ,
    p = xml_find_all(.x, ".//person/p") %>% xml_text() ,
    sex = xml_find_all(.x, ".//person/sex") %>% xml_text() ,
    age = xml_find_all(.x, ".//person/age") %>% xml_text() ,
    proband = xml_find_all(.x, ".//person/proband") %>% xml_text() ,
    note = xml_find_all(.x, ".//person/note") %>% xml_text() ,
    
    affected1 = xml_find_all(.x, ".//person/affected1") %>% xml_text() ,
    affected2 = xml_find_all(.x, ".//person/affected2") %>% xml_text() 
    
    
    )}) 




# then assign the area name column to the data frame
dat$area_name <- labs
