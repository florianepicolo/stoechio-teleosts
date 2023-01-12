###!/usr/bin/Rscript --slave


library(dplyr)
library(tidyr)
library(vroom)

setwd("~/Desktop/Pathways/prog_2023")


f.interacts <- "p-interactions-ensembl.csv"
f.orthologs <- "g-nbortho-teleosteens.csv"
d.interacts <- vroom(f.interacts)
d.orthologs <- vroom(f.orthologs)

d.orthologs <- d.orthologs %>% pivot_longer(!ensembl_id, names_to='sp', values_to='nc')

d.interacts
d.orthologs
