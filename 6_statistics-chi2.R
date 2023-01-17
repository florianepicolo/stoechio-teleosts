#!/usr/bin/Rscript --vanilla


library(dplyr)
library(vroom)
library(ggplot2)
library(gridExtra)
library(cowplot)
library(ggsignif)

setwd("~/Desktop/Pathways/stoechio-teleosts")

### Statistics chi2/hypergeom 

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


### Creation du barplot
data3 = read.csv("alamain-data-3WGD.csv", header = T, sep = ",")
data3$mean = rowMeans(data3[,3:56], na.rm = T)
data3$mean_int = round(data3$mean, 0)

str(data3$copy)
data3$copy = as.character(data3$copy)

for (i in 1:6){
  data3$max[i] = max(data3[i,3:56])
  data3$min[i] = min(data3[i,3:56])
  data3$sd[i] = sd(data3[i,3:56], na.rm = T)
  data3$sem[i] = data3$sd[i]/54

}

data4 = read.csv("alamain-data-4WGD.csv", header = T, sep = ",")
data4$mean = rowMeans(data4[,3:11], na.rm = T)
data4$mean_int = round(data4$mean, 0)


str(data4$copy)
data4$copy = as.character(data4$copy)

for (i in 1:6){
  data4$max[i] = max(data4[i,3:11])
  data4$min[i] = min(data4[i,3:11])
  data4$sd[i] = sd(data4[i,3:11], na.rm = T)
  data4$sem[i] = data4$sd[i]/9
  
}

#################
##   BARPLOT   ##
#################

## barplot total

limits = aes(ymax = mean + sem, ymin = mean - sem)

al_bar_3WGD = ggplot(data = data3, aes(x = list, y = mean, fill = copy)) +
  geom_bar(stat = "identity", color = "black", position = position_dodge()) +
  theme_minimal() + scale_fill_manual(values = c("#FFBF00", "#0489B1", "#BDBDBD")) +
  geom_text(aes(label = mean_int), vjust = -1.4, color = "black", position = position_dodge(0.9), size = 3.5) +
  geom_errorbar(limits, colour = 'black', width = 0.15, position = position_dodge(0.9)) +
  ylab("Number of genes") + xlab("") + ggtitle("Teleosts species with 3 WGD", subtitle = "(n=54)") +
  geom_segment(x=1, xend=2, y=5000, yend=5000, col="black") + # trait horizontal
  geom_segment(x=1, xend=1, y=5000, yend=4800, col="black") + # petit trait verticla à gauche 
  geom_segment(x=2, xend=2, y=5000, yend=4800, col="black") + # petit trait vertical à droite
  annotate("text", x = 1.5, y = 5050, label = "***")

al_bar_3WGD

al_bar_4WGD = ggplot(data = data4, aes(x = list, y = mean, fill = copy)) +
  geom_bar(stat = "identity", color = "black", position = position_dodge()) +
  theme_minimal() + scale_fill_manual(values = c("#FFBF00", "#0489B1", "#BDBDBD")) +
  geom_text(aes(label = mean_int), vjust = -1.4, color = "black", position = position_dodge(.9), size = 3.5) +
  geom_errorbar(limits, colour = 'black', width = 0.15, position = position_dodge(0.9)) +
  ylab("Number of genes") + xlab("") + ggtitle("Teleosts species with 4 WGD", subtitle = "(n=9)") +
  geom_segment(x=1.30, xend=2.3, y=5000, yend=5000, col="black") + # trait horizontal
  geom_segment(x=1.30, xend=1.30, y=5000, yend=4900, col="black") + # petit trait verticla à gauche 
  geom_segment(x=2.3, xend=2.3, y=5000, yend=4900, col="black") + # petit trait vertical à droite
  annotate("text", x = 1.8, y = 5050, label = "****")

al_bar_4WGD

plot_grid(al_bar_3WGD, al_bar_4WGD, labels=c("A", "B"), ncol = 2, nrow = 1, scale = c(1,1))


