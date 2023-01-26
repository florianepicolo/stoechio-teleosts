#!/usr/bin/Rscript --vanilla

library(dplyr)
library(tidyr)
library(vroom)
library(purrr)
library(stringr)
library(magrittr)
library(DT)
library(rstatix)
library(ggplot2)
library(ggmosaic)
library(ape)


## FICHIERS
setwd("~/Desktop/Pathways/stoechio-teleosts/")
f.interacts <- "p-min-interactions-ensembl.csv"
f.orthologs <- "g-nbortho-teleosteens-with-NT.csv"
f.sporder <- 'sp-order.csv'

d.interacts <- vroom(f.interacts)
d.orthologs <- vroom(f.orthologs, na = c("", "NA", "NT"))
d.sporder <- vroom(f.sporder, col_names=c(X1="species", X2="order"))

d.orthologs <- d.orthologs %>% pivot_longer(!ensembl_id, names_to='species', values_to='nc')
d.gp <- d.orthologs %>% select(ensembl_id) %>% unique() %>% filter(ensembl_id %in% d.interacts$interact1_ensembl | ensembl_id %in% d.interacts$interact2_ensembl)

# stoichiometrymétrie pour chaque interact et par sp
d.ip.nc <- semi_join(x=d.orthologs, y=d.gp) %>%
  left_join(d.interacts, by=c("ensembl_id"="interact1_ensembl")) %>% 
  filter(! is.na(interact2_ensembl)) %>%
  left_join(semi_join(x=d.orthologs, y=d.gp), by=c("interact2_ensembl"="ensembl_id", "species"="species"), suffix=c(".1", ".2")) %>%
  rename(ip.1=ensembl_id, ip.2=interact2_ensembl)

## sans retirer les 0:0 ni les NA
# do.ip.nc <- d.ip.nc %>% mutate(stoichiometry=sprintf("%s:%s", nc.1, nc.2)) %>% #=if_else(nc.1>nc.2, nc.2, nc.1) // =if_else(nc.1>nc.2, nc.1, nc.2) si on veut trier
#   select(ip.1, ip.2, species, stoichiometry) %>% unique() %>% pivot_wider(., names_from = species, values_from = stoichiometry)

## si on veut retirer tous les 0:0, ni les 0 dans toutes les espèces
do.ip.nc <- d.ip.nc %>% mutate(is.ok = ip.1 != ip.2) %>% filter(is.ok==TRUE) %>% select(-pathway) %>% group_by(ip.1, ip.2) %>% unique() %>%
  mutate(is.ok = !nc.1!=0 & !nc.2!=0 & (!is.na(nc.1) | !is.na(nc.2)) & ip.1 != ip.2, sum = sum(is.ok), sum1 = sum(nc.1!=0), sum2 = sum(nc.2!=0)) %>% filter(sum!=63, sum1!=0, sum2!=0) %>% 
  mutate(stoichiometry=sprintf("%s:%s", nc.1, nc.2)) %>% select(ip.1, ip.2, species, stoichiometry) %>% unique() %>% pivot_wider(., names_from = species, values_from = stoichiometry)

## si on veut garder uniquement les st == partout ou 5% près
do.st.nc <- d.ip.nc %>% mutate(is.ok = ip.1 != ip.2) %>% filter(is.ok==TRUE) %>% select(-pathway) %>% group_by(ip.1, ip.2) %>% unique() %>% 
  mutate(is.ok = !nc.1!=0 & !nc.2!=0 & (!is.na(nc.1) | !is.na(nc.2)) & ip.1 != ip.2, sum = sum(is.ok), sum1 = sum(nc.1!=0), sum2 = sum(nc.2!=0)) %>% 
  filter(sum!=63, sum1!=0, sum2!=0) %>% mutate(is.st = nc.1==nc.2, sumst = sum(is.st)) %>% filter(sumst>=60) %>%  mutate(stoichiometry=sprintf("%s:%s", nc.1, nc.2)) %>%
  select(ip.1, ip.2, species, stoichiometry) %>% unique() %>% pivot_wider(., names_from = species, values_from = stoichiometry)
  
## Supp data 2
# sortie <- write.table(do.ip.nc, file="p-min-interactions-by-teleosts-no00.csv", row.names = FALSE, col.names = TRUE, sep = ";")
# sortie <- write.table(do.st.nc, file="p-min-interactions-by-teleosts-st.csv", row.names = FALSE, col.names = TRUE, sep = ";")

### FIGURE 
# supprime les 0:0 pour toutes les espèces ! 
d.all_rel <- d.ip.nc %>% mutate(is.ok = ip.1 != ip.2) %>% filter(is.ok==TRUE) %>% 
  mutate(is.ok = !nc.1!=0 & !nc.2!=0 & (!is.na(nc.1) | !is.na(nc.2))) %>% group_by(ip.1, ip.2) %>% 
  mutate(sum = sum(is.ok)) %>% filter(sum!=63) %>% group_by(pathway, ip.1, ip.2) %>% 
  mutate(stoichiometry = sprintf("%s:%s", nc.1=if_else(nc.1>nc.2, nc.2, nc.1), nc.2=if_else(nc.1>nc.2, nc.1, nc.2))) %>% left_join(., d.sporder)
d.sp_rel <- d.all_rel %>%
  group_by(species) %>% mutate(nrel=n()) %>%
  group_by(species, nrel, stoichiometry, order) %>% summarize(n=n())

# on récupère les pct les plus élevées (>=5%)
perc.rel <- d.all_rel %>% ungroup %>% select(-pathway) %>% unique() %>% group_by(stoichiometry) %>% summarize(n=n()) %>% arrange(-n) %>% mutate(perc=n/sum(n)*100) 
rel.tokeep <- c('0:0','0:1','0:2','0:3','0:4','1:1','1:2','1:3','1:4','2:2','2:3','2:4','3:3','3:4','4:4')

sp_order <- d.sp_rel %>% filter(stoichiometry %in% rel.tokeep) %>% group_by(species) %>% mutate(freq=n/sum(n)) %>% filter(stoichiometry=="0:1") %>% 
  arrange(desc(order)) %$% species


gg.d.sp_rel <- d.sp_rel %>% 
  filter(stoichiometry %in%  rel.tokeep) %>% group_by(species) %>% mutate(freq=n/sum(n), species=factor(species, levels=c(sp_order))) %>% ggplot(.) +
  geom_mosaic(aes(x=product(species), weight=n, fill=stoichiometry)) + coord_flip() + ggtitle("All pathway")

# ggsave("graph-sp-ratio-order-with-00.png", width = 11, height = 8, plot = gg.d.sp_rel)


### PAREIL MAIS AVEC LES PATHWAYS ! (LA MEME FIGURE)
d.p_rel <- d.all_rel %>% group_by(pathway) %>% mutate(nrel=n()) %>% group_by(pathway, nrel, stoichiometry) %>% summarise(n=n())
gg.d.p_rel <- d.p_rel %>% filter(stoichiometry %in%  rel.tokeep) %>% group_by(pathway) %>% mutate(freq=n/sum(n)) %>% ggplot(.) +
  geom_mosaic(aes(x=product(pathway), weight=n, fill=stoichiometry)) + coord_flip() + ggtitle("All species")

# ggsave("graph-p-ratio-with-00.png", width = 11, height = 8, plot = gg.d.p_rel)




### LES POURCENTAGES DE STOECHIOMETRY
## si on veut garder toutes les rel.tokeep notamment pr les VRAIES STATS ! 
d.all.percbysp <- d.all_rel %>% filter(stoichiometry!="NA:NA") %>% ungroup %>% group_by(species, stoichiometry) %>% summarize(n=n()) %>% 
  arrange(-n) %>% mutate(perc=n/sum(n)*100) %>% select(-n) %>% filter(stoichiometry %in% rel.tokeep) %>% pivot_wider(., names_from = stoichiometry, values_from = perc)

sortie <- write.table(d.all.percbysp, file="p-min-interactions-perc-by-teleosts.csv", row.names = FALSE, col.names = TRUE, sep = ";", dec = ",")



# ###### 3WGD
d.ip.3WGD.nc <- d.ip.nc %>% filter(!species %in% c("Huchen", "Chinook_salmon", "Coho_salmon", "Rainbow_trout", "Atlantic_salmon", "Brown_trout", "Goldfish", "Common_carp", "Golden-line_barbel"))

d3.all_rel <- d.ip.3WGD.nc %>% group_by(pathway, ip.1, ip.2) %>%
  mutate(stoichiometry= sprintf("%s_%s", nc.1=if_else(nc.1>nc.2, nc.2, nc.1), nc.2=if_else(nc.1>nc.2, nc.1, nc.2))) %>%
  left_join(., d.sporder)

d3.sp_rel <- d3.all_rel %>%
  group_by(sp) %>% mutate(nrel=n()) %>%
  group_by(sp, nrel, stoichiometry, order) %>% summarize(n=n())

# on récupère les pct les plus élevées (>=5%)
rel3.tokeep <- d3.all_rel %>% ungroup %>% select(-pathway) %>% unique() %>% group_by(stoichiometry) %>% summarize(n=n()) %>% arrange(-n) %>% mutate(perc=n/sum(n)*100) %>% filter(perc>=5) %$% stoichiometry

# on rajoute les ratios qui ont un haut score également mais flouté par les précédents
########  pas sur de moi car pas de 0_3 ????????
d3.sp_rel %>% mutate(perc=n/sum(n)*100) %>% filter(perc>=5) %>% filter(!stoichiometry %in% rel3.tokeep) %>% arrange(-perc) %>% filter(stoichiometry!='0_0') %>% ungroup %>% select(stoichiometry) %>% unique
rel3.tokeep <- c(rel3.tokeep, '0_4', '0_3', '1_3', '2_3', '2_4', '3_4', '1_4', '4_4') ### j'ai rajouté 0_3 pr que ça suive une logique !!!!
sp3_order <- d3.sp_rel %>% filter(stoichiometry %in% rel.tokeep) %>% group_by(sp) %>% mutate(freq=n/sum(n)) %>% filter(stoichiometry=="0_1") %>% arrange(desc(order)) %$% sp


gg.d3.sp_rel <- d3.sp_rel %>%
  filter(stoichiometry %in%  rel3.tokeep) %>% group_by(sp) %>% mutate(freq=n/sum(n), species=factor(sp, levels=c(sp3_order))) %>% ggplot(.) +
  geom_mosaic(aes(x=product(species), weight=n, fill=stoichiometry)) + coord_flip() + ggtitle("All pathway")

# ###### 4WGD
d.ip.4WGD.nc <- d.ip.nc %>% filter(species %in% c("Huchen", "Chinook_salmon", "Coho_salmon", "Rainbow_trout", "Atlantic_salmon", "Brown_trout", "Goldfish", "Common_carp", "Golden-line_barbel"))

d4.all_rel <- d.ip.4WGD.nc %>% group_by(pathway, ip.1, ip.2) %>%
  mutate(stoichiometry= sprintf("%s_%s", nc.1=if_else(nc.1>nc.2, nc.2, nc.1), nc.2=if_else(nc.1>nc.2, nc.1, nc.2))) %>%
  left_join(., d.sporder)

d4.sp_rel <- d4.all_rel %>%
  group_by(species) %>% mutate(nrel=n()) %>%
  group_by(species, nrel, stoichiometry, order) %>% summarize(n=n())

# on récupère les pct les plus élevées (>=5%)
# rel4.tokeep <- d4.all_rel %>% ungroup %>% select(-pathway) %>% unique() %>% group_by(stoichiometry) %>% summarize(n=n()) %>% arrange(-n) %>% mutate(perc=n/sum(n)*100) %>% filter(perc>=5) %$% stoichiometry
rel4.tokeep <- c('0_1','0_2','0_3','0_4','1_1','1_2','1_3','1_4','2_2','2_3','2_4','3_3','3_4','4_4')

# on rajoute les ratios qui ont un haut score également mais flouté par les précédents
########  pas sur de moi car pas de 0_3 ????????
d4.sp_rel %>% mutate(perc=n/sum(n)*100) %>% filter(perc>=5) %>% filter(!stoichiometry %in% rel3.tokeep) %>% arrange(-perc) %>% filter(stoichiometry!='0_0') %>% ungroup %>% select(stoichiometry) %>% unique
# rel4.tokeep <- c(rel4.tokeep, '0_4', '0_3', '1_3', '2_3', '2_4', '3_4', '1_4', '4_4') ### j'ai rajouté 0_3 pr que ça suive une logique !!!!
sp4_order <- d4.sp_rel %>% filter(stoichiometry %in% rel.tokeep) %>% group_by(species) %>% mutate(freq=n/sum(n)) %>% filter(stoichiometry=="0_1") %>% arrange(desc(order)) %$% species


gg.d4.sp_rel <- d4.sp_rel %>%
  filter(stoichiometry %in%  rel4.tokeep) %>% group_by(species) %>% mutate(freq=n/sum(n), species=factor(species, levels=c(sp4_order))) %>% ggplot(.) +
  geom_mosaic(aes(x=product(species), weight=n, fill=stoichiometry)) + coord_flip() + ggtitle("All pathway")





