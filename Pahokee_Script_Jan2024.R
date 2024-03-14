#####PAHOKEE BACTERIA#####

setwd("~/Dropbox (UFL)/Laughinghouse_Lab/MANUSCRIPTS/Pahokee/Jan2024")


######## PACKAGES & SUCH ###########

#if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
#install.packages('ggthemes', dependencies = TRUE)
#install.packages('gg_ordiplot')
#install.packages("remotes")
#remotes::install_github("microbiome/microbiome")
#remotes::install_github("gauravsk/ranacapa")
#devtools::install_github("david-barnett/microViz")
#install.packages("remotes")
#remotes::install_github("jfq3/ggordiplots")
#install.packages("GUniFrac")
#devtools::install_github('microsud/microbiomeutilities')
#install.packages("eulerr")
#devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
#remotes::install_github("microbiota/amplicon", force = TRUE)
#remotes::install_github("ZhonghuiGai/veganEx")

library(veganEx)
library(amplicon)
library(patchwork)
library(pairwiseAdonis)
library(microbiomeutilities)
library(GUniFrac)
library(eulerr)
library(microViz)
library(gclus)
library(ggrepel)
library(ggforce)
library(tibble)
library(BiodiversityR)# also loads vegan
library(ranacapa)
library(microbiome)
library(tidyverse)
library(phyloseq)
library(ggsci)
library(ggthemes)
library(ade4)
library(adegraphics)
library(adespatial)
library(vegan3d)
library(MASS)
library(ellipse)
library(FactoMineR)
library(rrcov)
library(missMDA)
library(Amelia)
library(ggordiplots)
library(indicspecies)








#####Load object#####


Metadata <- read.csv("~/Dropbox (UFL)/Laughinghouse_Lab//MANUSCRIPTS/Pahokee/Jan2024/Pahokee_Metadata.csv") %>%
  mutate(across(c(HAT, Outing, Site), as.factor))

#Reads in bacterial physeq object
physeq = readRDS("~/Dropbox (UFL)/Laughinghouse_Lab/MANUSCRIPTS/Pahokee/Jan2024/PhyloSeq_Objects/Pahokee_Jan19.rds")
x = sample_sums(physeq)
#y = taxa_sums(physeq)
x %>% sum(.)

physeq@sam_data$Outing = as.factor(physeq@sam_data$Outing)
physeq@sam_data$HAT = as.factor(physeq@sam_data$HAT)
physeq@sam_data$Site = as.factor(physeq@sam_data$Site)


physeq_cyano = physeq %>% subset_taxa((Class == "Cyanophyceae"))




#####BAR PLOTS#####
#Phylum
bac_phylum = physeq %>% tax_fix() %>%
  ps_select(Site, Outing) %>%
  phyloseq::merge_samples(group = "Outing") %>%
  comp_barplot(tax_level = "Phylum", n_taxa = 10, nrow = 1, sample_order = "default") +
  labs(x = "Outing",
       y = "Relative Abundance")
bac_phylum

#Class
Class = physeq %>% tax_fix() %>%
  ps_select(Site, HAT) %>%
  phyloseq::merge_samples(group = "HAT") %>%
  comp_barplot(tax_level = "Class", n_taxa = 10, nrow = 1, sample_order = "default") +
  labs(#title = "Bacterial Class",
       #x = "HAT",
       y = "Relative Abundance")
Class

#Order
physeq %>% tax_fix() %>%
  ps_select(Site, Outing) %>%
  phyloseq::merge_samples(group = "Outing") %>%
  comp_barplot(tax_level = "Order", n_taxa = 10, nrow = 1, sample_order = "default") +
  labs(x = "Outing",
       y = "Relative Abundance")

#Family
physeq %>% tax_fix() %>%
  ps_select(Site, HAT) %>%
  phyloseq::merge_samples(group = "HAT") %>%
  comp_barplot(tax_level = "Family", n_taxa = 19, nrow = 1, sample_order = "default") +
  labs(x = "HAT",
       y = "Relative Abundance")

#Family-Cyano
physeq_cyano %>% tax_fix() %>%
  ps_select(Site, Outing) %>%
  phyloseq::merge_samples(group = "Outing") %>%
  comp_barplot(tax_level = "Family", n_taxa = 4, nrow = 1, sample_order = "default") +
  labs(x = "Outing",
       y = "Relative Abundance")


#Genus-Cyano
Genus_C = physeq_cyano %>% tax_fix() %>%
  ps_select(Site, HAT) %>%
  phyloseq::merge_samples(group = "HAT") %>%
  comp_barplot(tax_level = "Genus", n_taxa = 9, nrow = 1, sample_order = "default") +
  labs(#title = "Cyanobacteria",
       #x = "HAT",
       y = "Relative Abundance")

#Genus-Cyano
physeq_cyano %>% tax_fix() %>%
  ps_select(Site, HAT) %>%
  phyloseq::merge_samples(group = "HAT") %>%
  comp_barplot(tax_level = "unique", n_taxa = 10, nrow = 1, sample_order = "default") +
  labs(#x = "Outing",
       y = "Relative Abundance")


(Class / Genus_C / Euk_bar) + plot_annotation(tag_levels = 'A')


#####Alpha#####

p = physeq %>% tax_fix() %>%
    plot_richness(., x="HAT", color="HAT", measures=c("Observed", "Shannon"))

Bacteria_alpha = p + geom_violin(width=1, mapping = aes(fill = HAT)) +
  geom_boxplot(width=0.2, alpha=0.5) +
  geom_point(position = position_jitter(height = 0.2), alpha = 0.5) +
  theme_bw()  +
  labs(#title = "Bacteria"
    )
Bacteria_alpha

alpha_obs = as.data.frame(p$data) %>% dplyr::filter(variable == "Observed")
alpha_shan = as.data.frame(p$data) %>% dplyr::filter(variable == "Shannon")

obs_lm = lm(formula = value ~ HAT, data = alpha_obs) %>%
  emmeans::emmeans(., pairwise ~ HAT, type="response")
obs_lm
#contrast      estimate SE df t.ratio p.value
# HAT0 - HAT4      157.5 87.1 12   1.809  0.3160
#HAT0 - HAT24     173.5 87.1 12   1.993  0.2437
#HAT0 - HAT48     107.2 87.1 12   1.232  0.6198

shan_lm = lm(formula = value ~ HAT, data = alpha_shan) %>%
  emmeans::emmeans(., pairwise ~ HAT, type="response")
shan_lm
#contrast      estimate    SE df t.ratio p.value
#HAT0 - HAT4     0.1077 0.114 12   0.949  0.7799
#HAT0 - HAT24    0.2407 0.114 12   2.120  0.2017
#HAT0 - HAT48    0.2115 0.114 12   1.863  0.2934


#Cyanobacteria
p_C = physeq_cyano %>% tax_fix() %>%
  plot_richness(., x="HAT", color="HAT", measures=c("Observed", "Shannon"))

Cyano_alpha = p_C + geom_violin(width=1, mapping = aes(fill = HAT)) +
  geom_boxplot(width=0.2, alpha=0.5) +
  geom_point(position = position_jitter(height = 0.2), alpha = 0.5) +
  theme_bw()  +
  labs(#title = "Cyanobacteria",
       y = "")
Cyano_alpha


alpha_obs_C = as.data.frame(p_C$data) %>% dplyr::filter(variable == "Observed")
alpha_shan_C = as.data.frame(p_C$data) %>% dplyr::filter(variable == "Shannon")

cyano_obs_lm = lm(formula = value ~ HAT, data = alpha_obs_C) %>%
  emmeans::emmeans(., pairwise ~ HAT, type="response")
cyano_obs_lm
#contrast      estimate SE   df t.ratio p.value
#HAT0 - HAT4      26.75 6.64 12   4.031  0.0078
#HAT0 - HAT24      7.25 6.64 12   1.093  0.7005
#HAT0 - HAT48     -9.75 6.64 12  -1.469  0.4839

cyano_shan_lm = lm(formula = value ~ HAT, data = alpha_shan_C) %>%
  emmeans::emmeans(., pairwise ~ HAT, type="response")
cyano_shan_lm
#contrast      estimate    SE df t.ratio p.value
#HAT0 - HAT4      0.124 0.361 12   0.345  0.9852
#HAT0 - HAT24    -0.397 0.361 12  -1.101  0.6959
#HAT0 - HAT48    -0.541 0.361 12  -1.499  0.4675


Bacteria_alpha + Cyano_alpha  +
  plot_layout(guides = 'collect') + plot_annotation(tag_levels = 'A')

#####BETA#####

physeq_filt = tax_filter(physeq, min_sample_abundance = 100)
physeq_filt #252 taxa

physeq_cyano_filt = physeq_filt %>% subset_taxa((Class == "Cyanophyceae"))


#Bacteria
bacteria_pc = physeq_filt %>% tax_fix() %>%
  tax_transform("identity", rank = "unique") %>%
  dist_calc("aitchison") %>%
  ord_calc("PCoA") %>%
  ord_plot(color = "HAT", shape = "Site", size = 4, auto_caption = NA) +
  scale_colour_brewer(palette = "Dark2") +
  stat_ellipse(aes(colour = HAT), type = "t")



x <- physeq_filt %>%  tax_fix() %>%
  tax_transform("identity", rank = "unique") %>%
  dist_calc("aitchison")

w <- x %>%
  dist_permanova(
    seed = 1234, # for set.seed to ensure reproducibility of random process
    n_processes = 1, n_perms = 99999, # you should use at least 999!
    variables = c("HAT")
  )

# view the permanova results
perm_get(w) %>% as.data.frame()
#         Df SumOfSqs        R2        F Pr(>F)
#HAT       3  886.7498 0.6565475 7.646442  1e-05
#Residual 12  463.8758 0.3434525       NA     NA
#Total    15 1350.6256 1.0000000       NA     NA

b = dist_get(x)
b

metadata4 = Metadata %>% arrange(Name) %>% column_to_rownames(var="Name")


pairwise_site = pairwiseAdonis::pairwise.adonis(b, metadata4[,"HAT"],perm = 99999, p.adjust.m = "holm")
pairwise_site %>% filter(p.adjusted <= 0.05)
#NA

#use ANOSIM
anosim.pairwise(
  b,
  metadata4[,"HAT"],
  perm = 9999
)

#     pairs   anosimR p.value  p.adj
#1  0.vs.24 0.9062500  0.0271 0.1626
#2  0.vs.48 1.0000000  0.0305 0.1830
#3   0.vs.4 0.6354167  0.0275 0.1650


#Cyanobacteria
cyano_pc = physeq_cyano_filt %>% tax_fix() %>%
  tax_transform("identity", rank = "unique") %>% # don't transform!
  dist_calc("aitchison") %>%
  ord_calc("PCoA") %>%
  ord_plot(color = "HAT", shape = "Site", size = 4, auto_caption = NA) +
  scale_colour_brewer(palette = "Dark2") +
  stat_ellipse(aes(colour = HAT), type = "t")


x <- physeq_cyano_filt %>%  tax_fix() %>%
  tax_transform("identity", rank = "unique") %>%
  dist_calc(dist = "aitchison")

w <- x %>%
  dist_permanova(
    seed = 1234, # for set.seed to ensure reproducibility of random process
    n_processes = 1, n_perms = 99999, # you should use at least 999!
    variables = c("HAT")
  )

# view the permanova results
perm_get(w) %>% as.data.frame()
#         Df SumOfSqs        R2        F Pr(>F)
#HAT       3  826.3427 0.2617968 1.418562 0.01438
#Residual 12 2330.0851 0.7382032       NA      NA
#Total    15 3156.4278 1.0000000       NA      NA


b = dist_get(x)
b


pairwise_site = pairwiseAdonis::pairwise.adonis(b, metadata4[,"HAT"],perm = 99999)
pairwise_site %>% filter(p.adjusted <= 0.05)
#     pairs Df SumsOfSqs    F.Model        R2    p.value p.adjusted sig
#1  0 vs 24  1  3.116623  1.3140415 0.1796601 0.25714286  1.0000000
#2  0 vs 48  1  8.114241  3.5406112 0.3711095 0.11428571  0.6857143
#3   0 vs 4  1  1.609696  0.8269668 0.1211324 0.37142857  1.0000000

anosim.pairwise(
  b,
  metadata4[,"HAT"],
  perm = 9999
)

#     pairs   anosimR p.value  p.adj
#1  0.vs.24 0.08333333  0.2205 1.0000
#2  0.vs.48 0.13541667  0.2604 1.0000
#3   0.vs.4 0.38541667  0.0280 0.1680



(bacteria_pc + cyano_pc + Euk_PCoA) +
  plot_layout(guides = 'collect') + plot_annotation(tag_levels = 'A')


#####TOXINS#####

data_sum_mc_tot = plyr::ddply(Metadata, c("HAT"), summarise,
                              N    = length(Total_MC_ppb),
                              mean_MC = mean(Total_MC_ppb),
                              sd   = sd(Total_MC_ppb),
                              se   = sd / sqrt(N)
)

data_sum_mc_aq = plyr::ddply(Metadata, c("HAT"), summarise,
                              N    = length(Aq_MC_ppb),
                              mean_MC = mean(Aq_MC_ppb),
                              sd   = sd(Aq_MC_ppb),
                              se   = sd / sqrt(N)
)

x = ggplot(data_sum_mc_tot, aes(x = HAT, y = mean_MC, fill = HAT, color = HAT))+
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin = mean_MC-se, ymax = mean_MC+se) , width = 0.6, color = "black") +
  theme_bw() +
  theme(legend.position = "none")+
  labs(#title = "Total Microcystins Pahokee Marine",
       #subtitle = "Mean within the marina",
       x = "Hours After Treatment",
       y = "Microcystins ppb"
  )

x2 = ggplot(data_sum_mc_aq, aes(x = HAT, y = mean_MC, fill = HAT, color = HAT))+
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin = mean_MC-se, ymax = mean_MC+se) , width = 0.6, color = "black") +
  theme_bw() +
  theme(legend.position = "none")+
  labs(#title = "Total Microcystins Pahokee Marine",
    #subtitle = "Mean within the marina",
    x = "Hours After Treatment",
    y = ""
  )


tot = ggplot(Metadata, aes(x = HAT, y = Total_MC_ppb, fill = HAT, color = HAT))+
  geom_bar(stat="identity") +
  facet_grid(.~ Site)+
  theme_bw() +
  theme(legend.position = "none")+
  labs(#title = "Total Microcystins",
       #subtitle = "Mean within the marina",
       x = "Hours After Treatment",
       y = "Microcystins ppb"
  )

aq = ggplot(Metadata, aes(x = HAT, y = Aq_MC_ppb, fill = HAT, color = HAT))+
  geom_bar(stat="identity") +
  facet_grid(.~ Site)+
  theme_bw() +
  theme(legend.position = "none")+
  labs(#title = "Aqueous Microcystins",
       #subtitle = "Mean within the marina",
       x = "Hours After Treatment",
       y = " "
  )

(tot + aq) /
(x + x2) + plot_annotation(tag_levels = 'A')

(tot + x) /
(aq + x2) + plot_annotation(tag_levels = 'A')


Toxins_lm = lm(formula = Total_MC_ppb ~ HAT, data = Metadata) %>%
  emmeans::emmeans(., pairwise ~ HAT, type="response")

res_aov <- aov(Metadata$Total_MC_ppb ~ Metadata$HAT)
summary(res_aov)
TukeyHSD(res_aov, conf.level=.95)













#####Abundance#####

chl = ggplot(Metadata, aes(x = HAT, y = ChlaugL, fill = HAT)) +
  geom_bar(stat="identity", fill = "green4") +
  facet_grid(.~ Site)+
  theme_bw() +
  theme(legend.position = "none")+
  labs(#title = "Chlorophyll a",
       #subtitle = "Mean within the marina",
       x = "",
       y = "Chlorophyll a ug/L"
  )


data_sum_mc_chl = plyr::ddply(Metadata, c("HAT"), summarise,
                              N    = length(ChlaugL),
                              mean_MC = mean(ChlaugL),
                              sd   = sd(ChlaugL),
                              se   = sd / sqrt(N)
)

chl_mean = ggplot(data_sum_mc_chl, aes(x = HAT, y = mean_MC, fill = HAT, color = HAT))+
  geom_bar(stat="identity", fill = "green4") +
  geom_errorbar(aes(ymin = mean_MC-se, ymax = mean_MC+se) , width = 0.6, color = "black") +
  theme_bw() +
  theme(legend.position = "none")+
  labs(#title = "Total Microcystins Pahokee Marine",
    #subtitle = "Mean within the marina",
    x = "",
    y = "Chlorophyll a ug/L"
  )


pc = ggplot(Metadata, aes(x = HAT, y = Phycocyanin_ugL, fill = HAT, color = HAT))+
  geom_bar(stat="identity", fill = "steelblue") +
  facet_grid(.~ Site)+
  theme_bw() +
  theme(legend.position = "none")+
  labs(#title = "Phycocyanin",
       #subtitle = "Mean within the marina",
       x = "",
       y = "Phycocyanin mg/L"
  )

data_sum_mc_pc = plyr::ddply(Metadata, c("HAT"), summarise,
                              N    = length(Phycocyanin_ugL),
                              mean_MC = mean(Phycocyanin_ugL),
                              sd   = sd(Phycocyanin_ugL),
                              se   = sd / sqrt(N)
)

PC_mean = ggplot(data_sum_mc_pc, aes(x = HAT, y = mean_MC, fill = HAT, color = HAT))+
  geom_bar(stat="identity", fill = "steelblue") +
  geom_errorbar(aes(ymin = mean_MC-se, ymax = mean_MC+se) , width = 0.6, color = "black") +
  theme_bw() +
  theme(legend.position = "none")+
  labs(#title = "Total Microcystins Pahokee Marine",
    #subtitle = "Mean within the marina",
    x = "",
    y = "Phycocyanin mg/L"
  )

count = ggplot(Metadata, aes(x = HAT, y = Microcysits_counts, fill = HAT, color = HAT))+
  geom_bar(stat="identity") +
  facet_grid(.~ Site)+
  theme_bw() +
  theme(legend.position = "none")+
  labs(#title = "Microcystis counts",
       #subtitle = "Mean within the marina",
       x = "Hours After Treatment",
       y = "Microcystis cells/mL"
  )

data_sum_mc_count = plyr::ddply(Metadata, c("HAT"), summarise,
                             N    = length(Microcysits_counts),
                             mean_MC = mean(Microcysits_counts),
                             sd   = sd(Microcysits_counts),
                             se   = sd / sqrt(N)
)

count_mean = ggplot(data_sum_mc_count, aes(x = HAT, y = mean_MC, fill = HAT, color = HAT))+
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin = mean_MC-se, ymax = mean_MC+se) , width = 0.6, color = "black") +
  theme_bw() +
  theme(legend.position = "none")+
  labs(#title = "Total Microcystins Pahokee Marine",
    #subtitle = "Mean within the marina",
    x = "Hours After Treatment",
    y = "Microcystis cells/mL"
  )


(chl / pc / count) + plot_annotation(tag_levels = 'A')
(chl_mean / PC_mean / count_mean) + plot_annotation(tag_levels = 'A')











#####RDA#####

physeq_filt %>% tax_fix(unknowns = c("Incertae-Sedis")) %>%
  tax_transform("rclr", rank = "Family") %>%
  ord_calc(
    constraints = c("ChlaugL", "Aq_MC_ppb", "Total_MC_ppb"),
    scale_cc = FALSE # doesn't make a difference
  ) %>%
  ord_plot(
    colour = "HAT", size = 2,
    plot_taxa = 1:15, tax_vec_length = 4.5, tax_lab_length = 4.6,
    tax_lab_style = tax_lab_style(
      type = "text", max_angle = 90, fontface = "bold.italic"
    ),
    constraint_vec_style = vec_constraint(1.5, alpha = 0.5),
    constraint_vec_length = 3, constraint_lab_length = 3.3,
    constraint_lab_style = constraint_lab_style(
      alpha = 0.8, size = 3, max_angle = 90, perpendicular = TRUE
    )
  ) +
  stat_ellipse(aes(colour = HAT), type = "t")

#####PCA#####
#Bacteria

PCA = physeq_filt %>% tax_fix() %>%
  tax_transform("identity", rank = "Family") %>%
  #dist_calc("aitchison") %>%
  ord_calc("PCA") %>%
  ord_plot(color = "HAT", shape = "Site", size = 4, auto_caption = NA, plot_taxa = 1:10,
           #tax_vec_length = 4.5, tax_lab_length = 4.6,
           tax_lab_style = tax_lab_style(
             type = "text", max_angle = 90, fontface = "bold"
           ),
           constraint_vec_style = vec_constraint(1.5, alpha = 0.5),
           constraint_vec_length = 3, constraint_lab_length = 3.3,
           constraint_lab_style = constraint_lab_style(
             alpha = 0.8, size = 3, max_angle = 90, perpendicular = TRUE
           )
  ) +
  scale_colour_brewer(palette = "Dark2") +
  stat_ellipse(aes(colour = HAT), type = "t")

(PCA / PCA_EUK) +
  plot_layout(guides = 'collect') + plot_annotation(tag_levels = 'A')
