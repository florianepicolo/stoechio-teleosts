#!/usr/bin/Rscript --vanilla


library(dplyr)
library(vroom)

setwd("~/Desktop/Pathways/stoechio-teleosts")


# f.stats.teleosts <- "p-min-stats-chi2-teleosteens.csv"
f.stats.teleosts <- "p-max-stats-chi2-teleosteens.csv"
d.stats.teleosts <- vroom(f.stats.teleosts)

d.stats.teleosts <- d.stats.teleosts %>% group_by(species) %>%
  mutate(nb_g_total = sum(nb_g_single, nb_g_dupli, nb_g_tripli),
         nb_p_total = sum(nb_p_single, nb_p_dupli, nb_p_tripli),
         d_chi2_pval = chisq.test(x=c(nb_p_single + nb_p_tripli, nb_p_dupli), 
                              p=c(nb_g_single/nb_g_total + nb_g_tripli/nb_g_total, nb_g_dupli/nb_g_total), correct=FALSE)$p.value,
         t_chi2_pval = chisq.test(x=c(nb_p_single + nb_p_dupli,  nb_p_tripli), 
                                  p=c(nb_g_single/nb_g_total + nb_g_dupli/nb_g_total, nb_g_tripli/nb_g_total), correct=FALSE)$p.value) %>%
  ungroup() %>% mutate(d_chi2_pvalBH = p.adjust(d_chi2_pval, method = "BH", n = length(d_chi2_pval)),
                       t_chi2_pvalBH = p.adjust(t_chi2_pval, method = "BH", n = length(t_chi2_pval))) %>%
  group_by(species) %>%
  mutate(d_hypergeom_pval = phyper((nb_p_dupli - 1), nb_g_dupli,
                            (nb_g_total - nb_g_dupli), nb_p_total, lower.tail = FALSE),
         t_hypergeom_pval =  phyper((nb_p_tripli - 1), nb_g_tripli,
                                    (nb_g_total - nb_g_tripli), nb_p_total, lower.tail = FALSE)) %>%
  ungroup() %>% mutate(d_hypergeom_pvalBH = p.adjust(d_hypergeom_pval, method = "BH", n = length(d_hypergeom_pval)),
                       t_hypergeom_pvalBH = p.adjust(t_hypergeom_pval, method = "BH", n = length(t_hypergeom_pval))) %>%
  group_by(species) %>%
  mutate(fc_single = (nb_p_single/nb_p_total)/(nb_g_single/nb_g_total),
         fc_dupli = (nb_p_dupli/nb_p_total)/(nb_g_dupli/nb_g_total),
         fc_tripli = (nb_p_tripli/nb_p_total)/(nb_g_tripli/nb_g_total))

do.stats.teleosts <- d.stats.teleosts %>% ungroup 
sortie <- write.table(do.stats.teleosts, file=f.stats.teleosts, row.names = FALSE, col.names = TRUE, sep = ";", dec = ",")
