library("phyloseq")
library("data.table")
library("plyr")
library("dplyr")
library("tidyr")
library("ggplot2")
library("reshape2")
library("indicspecies")
library("ggnetwork")
library("ape")
library("microbiome")
library("ggthemes")
library("cowplot")
library("ggsignif")

setwd("~/Mote_nutrient_experiment/genotype_7")

# load the rarefied, renamed, relative abundance data table with all species
load(file = "ps_rel_g7.RData") #renamed, rarefied ps object
sample_sums(ps_rel) #should be 1

ps_rel = subset_samples(ps_rel, SampleID != "HP-7-2-T2") #I removed these because they had a ton of RICA and are clearly actually 50
ps_rel = subset_samples(ps_rel, SampleID != "MP-7-3-T2")
mapfile = "~/Mote_nutrient_experiment/data/Mapping-file-full-renamed.txt"
map = import_qiime_sample_data(mapfile)
sample_data(ps_rel) <- map
nsamples(ps_rel)

# melt the data at the Genus level
ps_rel_genus_melt <- ps_rel %>%
  tax_glom(taxrank = "Genus") %>%
  psmelt()
head(ps_rel_genus_melt)

# Get families with mean relative abundance >0.01 across all samples 
gen_sum <- ps_rel_genus_melt %>% group_by(Genus) %>% dplyr::summarise(Aver = mean(Abundance))
gen_sub <- gen_sum[which(gen_sum$Aver > 0.01),]
names <- gen_sub$Genus
names

# Replace genus with <0.01 abundance with "NA"
ps_rel_genus_melt$genus <- ps_rel_genus_melt$Genus

ps_rel_genus_melt$genus[ps_rel_genus_melt$genus != "[Caedibacter] taeniospiralis group" & 
                          ps_rel_genus_melt$genus != "Aquarickettsia" &
                          ps_rel_genus_melt$genus != "Family_Helicobacteraceae" &
                          ps_rel_genus_melt$genus != "Thalassotalea" &
                          ps_rel_genus_melt$genus != "Unclassified_ASV_1" &
                          ps_rel_genus_melt$genus != "Unclassified_ASV_4" &
                          ps_rel_genus_melt$genus != "Unclassified_ASV_5" &
                          ps_rel_genus_melt$genus != "Vibrio"] <- NA

# plot
bar_genus = ggplot(ps_rel_genus_melt, aes(x = reorder(NW_T0_merge, Exposure_weeks), y=Abundance)) + 
  geom_bar(stat="identity", position="fill", aes(fill = reorder(genus, Abundance))) +
  scale_fill_manual(values=c("#660066", "#0072B2", "#CC3300","#FFFF00", "#58B8EC", "#8CD09F", "#AA2384","#56B4E9", "#E27AB7"), 
                    breaks = c("[Caedibacter] taeniospiralis group", "Aquarickettsia", "Family_Helicobacteraceae", "Thalassotalea", "Unclassified_ASV_1", "Unclassified_ASV_4", "Unclassified_ASV_5", "Unclassified_ASV_10", "Vibrio", "Halobacteriovorax"),
                    labels = c("Cysteiniphilum", "Aquarickettsia", "Family Helicobacteraceae", "Thalassotalea", "Unclassified ASV 1", "Unclassified ASV 4", "Unclassified ASV 5", "Unclassified ASV 10", "Vibrio", "Halobacteriovorax"),
                    na.value = "transparent")  +
  facet_grid(~reorder(Nutrient, Nutrient_order), scales = "free_x", space = "free_x") +
  theme_bw() +
  #theme(legend.text = element_text(face = "italic")) +
  ylab("Relative Abundance") +
  xlab("Treatment Weeks") +
  theme(axis.text.x = element_text(angle = 50, hjust = 1)) +
  labs(fill = "Genus", size = "Genus")
bar_genus

ggsave(filename="~/Mote_nutrient_experiment/genotype_7/plots/geno_7_bar_genus.png", plot=bar_genus, device="png",  width = 7, height = 4, dpi=500)

