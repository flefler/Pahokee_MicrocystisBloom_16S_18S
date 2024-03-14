library(phyloseq)
library(tidyverse)

setwd("~/Dropbox (UFL)/Laughinghouse_Lab/MANUSCRIPTS/Pahokee/Jan2024")
#####
#saveRDS(Euks_Data, file="Pahokee_Euks")

Euks_Data = readRDS("E:/Dropbox (UFL)/Laughinghouse_Lab/MANUSCRIPTS/Pahokee/Jan2024/PhyloSeq_Objects/Pahokee_Euks")

#Euks_Data_clean = Euks_Data %>%
#  subset_taxa((Domain == "Eukaryota") | is.na(Domain)) %>%
#  subset_taxa((Division != "Metazoa") | is.na(Division)) %>%
#  subset_taxa((Division != "Fungi") | is.na(Division)) %>%
#  subset_taxa((Supergroup != "NA") | is.na(Supergroup))


#otutable = (Euks_Data_clean@otu_table) %>% as.data.frame()
#taxtable = (Euks_Data_clean@tax_table) %>% as.data.frame()

#vegan::rarecurve(otutable, step=50, cex=0.5)
#?rarefy_even_depth
ps_rarefied <- rarefy_even_depth(Euks_Data,
                                 rngseed = 42069,
                                 sample.size = min(sample_sums(Euks_Data)),
                                 replace = FALSE)
ps_rarefied



summarize_phyloseq(ps_rarefied)

#otutable_rare = (ps_rarefied@otu_table) %>% as.data.frame()

#vegan::rarecurve(otutable_rare, step=50, cex=0.5)


#####


Metadata = as.data.frame(ps_rarefied@sam_data)


library(microViz)
Euk_bar = ps_rarefied %>% tax_fix() %>%
  ps_select(Site, HAT) %>% 
  phyloseq::merge_samples(group = "HAT") %>%
  comp_barplot(tax_level = "Genus", n_taxa = 20, nrow = 1, sample_order = "default") +
  labs(x = "HAT", 
       y = "Relative Abundance") 
Euk_bar

Euk_PCoA = ps_rarefied %>% tax_fix() %>%
  tax_transform("identity", rank = "Genus") %>%
  dist_calc(dist = "bray", binary = FALSE) %>%
  ord_calc("PCoA") %>%
  ord_plot(color = "HAT", size = 4, auto_caption = NA, plot_taxa = 1:10, 
           tax_vec_length = 4.5, tax_lab_length = 4.6,
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


Euk_PCoA


x <- ps_rarefied %>% tax_fix() %>%
  tax_transform("identity", rank = "Genus") %>%
  dist_calc(dist = "bray", binary = FALSE)

w <- x %>%
  dist_permanova(
    seed = 1234, # for set.seed to ensure reproducibility of random process
    n_processes = 1, n_perms = 999, # you should use at least 999!
    variables = c("HAT")
  )

# view the permanova results
perm_get(w) %>% as.data.frame()

#Df SumOfSqs        R2        F  Pr(>F)
#HAT       3 1.223574 0.2146392 1.093201 0.17965
#Residual 12 4.477034 0.7853608       NA      NA
#Total    15 5.700608 1.0000000       NA      NA

b = dist_get(x)
b2 <- sqrt(b)

metadata4 = Metadata %>% arrange(Name) %>% column_to_rownames(var="Name")


pairwise_site = pairwiseAdonis::pairwise.adonis(b2, Metadata$HAT,perm = 999, p.adjust.m = "fdr")

pairwise_site

pairwise_site %>% filter(p.adjusted <= 0.05)
#NA
#remotes::install_github("ZhonghuiGai/veganEx")
library('veganEx')

anosim.pairwise(
  b2,
  Metadata$HAT,
  perm = 999,
  p.adjust.m = "fdr"
)

#?anosim.pairwise

PCA_EUK = ps_rarefied %>% tax_fix() %>%
  tax_transform("identity", rank = "Class") %>%
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

PCA_EUK

ps <- tax_fix(ps_rarefied)

ps%>%
  tax_agg("Genus") %>%
  tax_transform("clr") %>%
  cor_heatmap(
    cor = "spearman", 
    taxa = tax_top(ps, 10, by = max, rank = "Genus"),
    vars = c("Phycocyanin.ug.L", "Total_MC_ppb", "Aq_MC_ppb", "Microcysits_counts")
  )




# compute correlations, with p values, and store in a dataframe
correlations_df <- ps %>% 
  tax_model(
    trans = "clr",
    rank = "Genus", variables =  list("Aq_MC_ppb", "Microcysits_counts"), 
    taxa = tax_top(ps, 10, by = max, rank = "Genus"), type = microViz::cor_test, method = "spearman", 
    return_psx = FALSE, verbose = FALSE, exact=FALSE
  ) %>% 
  tax_models2stats(rank = "Genus")

#?cor


# get nice looking ordering of correlation estimates using hclust
taxa_hclust <- correlations_df %>% 
  dplyr::select(term, taxon, estimate) %>% 
  tidyr::pivot_wider(
    id_cols = taxon, names_from = term, values_from = estimate
  ) %>% 
  tibble::column_to_rownames("taxon") %>% 
  as.matrix() %>% 
  stats::dist(method = "euclidean") %>% 
  hclust(method = "ward.D2") 

taxa_order <- taxa_hclust$labels[taxa_hclust$order]



#library(ggplot2)

correlations_df %>% 
  mutate(p.FDR = p.adjust(p.value, method = "fdr")) %>% 
  ggplot(aes(x = term, y = taxon)) +
  geom_raster(aes(fill = estimate)) +
  geom_point(
    data = function(x) filter(x, p.value < 0.05),
    shape = "asterisk"
  ) +
  geom_point(
    data = function(x) filter(x, p.FDR < 0.05),
    shape = "circle", size = 3
  ) +
  scale_y_discrete(limits = taxa_order) +
  scale_fill_distiller(palette = "BrBG", limits = c(-1, 1)) + 
  theme_minimal() +
  labs(
    x = NULL, y = NULL, fill = "Spearman's\nRank\nCorrelation",
    caption = paste(
      "Asterisk indicates p < 0.05, not FDR adjusted",
      "Filled circle indicates FDR corrected p < 0.05", sep = "\n"
    ))






