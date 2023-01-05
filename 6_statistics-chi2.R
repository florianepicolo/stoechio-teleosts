###!/usr/bin/Rscript --slave


library(dplyr)
library(vroom)

setwd("~/Desktop/Pathways/prog_2023")


f.stats.teleosts <- "p-stats-chi2-teleosteens.csv"
d.stats.teleosts <- vroom(f.stats.teleosts)

d.stats.teleosts <- d.stats.teleosts %>% group_by(species) %>%
  mutate(nb_g_total = sum(nb_g_single, nb_g_dupli, nb_g_tripli),
         nb_p_total = sum(nb_p_single, nb_p_dupli, nb_p_tripli),
         chi2_pval = chisq.test(x=c(nb_p_single, nb_p_dupli + nb_p_tripli), 
                              p=c(nb_g_single/nb_g_total, nb_g_dupli/nb_g_total + nb_g_tripli/nb_g_total))$p.value) %>%
  ungroup() %>% mutate(chi2_pvalBH = p.adjust(chi2_pval, method = "BH", n = length(chi2_pval))) %>%
  group_by(species) %>%
  mutate(hypergeom_pval = phyper((nb_p_tripli + nb_p_dupli - 1), nb_g_tripli + nb_g_dupli,
                            (nb_g_total - nb_g_tripli - nb_g_dupli), nb_p_total, lower.tail = FALSE )) %>%
  ungroup() %>% mutate(hypergeom_pvalBH = p.adjust(hypergeom_pval, method = "BH", n = length(hypergeom_pval))) %>%
  group_by(species) %>%
  mutate(fc_single = (nb_p_single/nb_p_total)/(nb_g_single/nb_g_total),
         fc_dupli = (nb_p_dupli/nb_p_total)/(nb_g_dupli/nb_g_total),
         fc_tripli = (nb_p_tripli/nb_p_total)/(nb_g_tripli/nb_g_total))

do.stats.teleosts <- d.stats.teleosts %>% ungroup 
sortie <- write.table(do.stats.teleosts, file=f.stats.teleosts, row.names = FALSE, col.names = TRUE, sep = ";", dec = ",")


