library(exactRankTests)
library(dplyr)
library(ggplot2)
library(compositions)
library(phyloseq)
library(robustbase)
library(ggpubr)

source("~/Basic_script_formats/ancom_v2.1.R") #sourced from https://github.com/FrederickHuangLin/ANCOM/tree/master/scripts
setwd("~/Mote_nutrient_experiment/genotype_7")

# load unrarefied, renamed data for differential abundance analysis
load(file = "ps_g7.RData")
ps_rel = subset_samples(ps, SampleID != "HP-7-2-T2") #I removed these because they had a ton of RICA and are clearly actually 50
ps_rel = subset_samples(ps, SampleID != "MP-7-3-T2")
sample_sums(ps)

mapfile = "~/Mote_nutrient_experiment/data/Mapping-file-full-renamed.txt"
map = import_qiime_sample_data(mapfile)
sample_data(ps) <- map

# filter to keep only most abundant taxa
ps = filter_taxa(ps, function(x) sum(x > 20) > (0.1*length(x)), TRUE)
# agglomerate to genus
ps <- tax_glom(ps, "Genus")
ntaxa(ps) #35
summary(sample_sums(ps)) #median: 11804

ps <- subset_samples(ps, No_level_no_T0 == "Phosphate")

ps <- subset_samples(ps, No_level_no_T0 != "No Treatment")
ps <- subset_samples(ps, No_level_no_T0 == "No Treatment")
nsamples(ps)
ntaxa(ps)
summary(sample_sums(ps))

OTUdf <- as.data.frame(t(otu_table(ps)))
metadf <- read.delim(file = "~/Mote_nutrient_experiment/data/Mapping-file-full-renamed.txt")
colnames(metadf)[1] <- "Sample.ID"
metadf$Exposure_weeks <- as.factor(metadf$Exposure_weeks)
levels(metadf$Exposure_weeks)

#for reciprocal comparison
metadf$Exposure_weeks <- factor(metadf$Exposure_weeks, levels = c("3", "0", "6")) 
levels(metadf$Exposure_weeks)                    

# Step 1: Data preprocessing
feature_table = OTUdf; meta_data = metadf; sample_var = "Sample.ID"; group_var = NULL
out_cut = 0.05; zero_cut = 0.90; lib_cut = 1000; neg_lb = FALSE
prepro = feature_table_pre_process(feature_table, meta_data, sample_var, group_var, 
                                   out_cut, zero_cut, lib_cut, neg_lb)

#ASV11 is only present in one no treatment sample
#zero_cut: Taxa with proportion of zeroes greater than zero_cut are not included in the analysis.
#out_cut: observations with proportion of mixture distribution less than out_cut (or greater than 1- out_cut) 
#will be detected as outliers

feature_table = prepro$feature_table # Preprocessed feature table
meta_data = prepro$meta_data # Preprocessed metadata
struc_zero = prepro$structure_zeros # Structural zero info

# Step 2: ANCOM
main_var = "Exposure_weeks"; p_adj_method = "fdr"; alpha = 0.05
adj_formula = NULL ; rand_formula = "~ 1 | Tank"

res = ANCOM(feature_table, meta_data, struc_zero, main_var, p_adj_method, 
            alpha, adj_formula, rand_formula)

resdf <- as.data.frame(res$out)

# add taxonomy
tax<-as(tax_table(ps),"matrix")
tax_cols <- c("Kingdom", "Phylum", "Class", "Order","Family","Genus", "Species")
tax<-as.data.frame(tax)
colnames(tax)
tax$taxa_id <- rownames(tax)
rownames(tax) <- NULL
dim(tax)
tax <- tax[,c(8,1:7)] #taxID
keep <- tax[tax$taxa_id %in% rownames(feature_table), ]
tax <- keep
rm(keep)

# Number of taxa except structural zeros
n_taxa = ifelse(is.null(struc_zero), nrow(feature_table), sum(apply(struc_zero, 1, sum) == 0))
# Cutoff values for declaring differentially abundant taxa
cut_off = c(0.9 * (n_taxa -1), 0.8 * (n_taxa -1), 0.7 * (n_taxa -1), 0.6 * (n_taxa -1))
names(cut_off) = c("detected_0.9", "detected_0.8", "detected_0.7", "detected_0.6")

# Annotation data
dat_ann = data.frame(x = min(res$fig$data$x), y = cut_off["detected_0.8"], label = "W[0.8]")

#test figure
fig = res$fig +  
   geom_hline(yintercept = cut_off["detected_0.7"], linetype = "dashed") + 
  geom_text(data = dat_ann, aes(x = x, y = y, label = label), 
              size = 4, vjust = -0.5, hjust = 0, color = "orange", parse = TRUE)
fig  
 
# # specialized plot
figdf <- as.data.frame(fig$data)
figdf <- cbind(figdf,tax, by = "taxa.id")
#figdf <- cbind(figdf,tax, by = "taxa_id")
figdf <- figdf[,-6]
figdf$group <- as.factor(figdf$group)
levels(figdf$group) <- c("3 vs. 0 weeks", "6 vs. 0 weeks")

#reciprocal contrast
levels(figdf$group) <- c("0 vs. 3 weeks", "6 vs. 3 weeks")

# order genus
x = tapply(figdf$y, figdf$Genus, function(x) max(x))
x = sort(x, TRUE)
figdf$Genus = factor(as.character(figdf$Genus), levels=names(x))
figdf$col_genus <- figdf$Genus

#use resdf to identify taxa significant at W statistic of 0.7/0.8

#no treatment
figdf$col_genus[figdf$col_genus != "Ferrimonas"] <- NA

#nitrate
figdf$col_genus[figdf$col_genus != "[Caedibacter] taeniospiralis group" &
                  figdf$col_genus != "Haloferula" &
                  figdf$col_genus != "Candidatus Tenderia"] <- NA

#ammonium
figdf$col_genus[figdf$col_genus != "P3OB-42"] <- NA

#combined
figdf$col_genus[figdf$col_genus != "Unclassified_ASV_5" &
                  figdf$col_genus != "[Caedibacter] taeniospiralis group" &
                  figdf$col_genus != "Ruegeria"] <- NA

#geno 7 nutrient only
figdf$col_genus[figdf$col_genus != "[Caedibacter] taeniospiralis group" &
                  figdf$col_genus != "Unclassified_ASV_5" &
                  figdf$col_genus != "Unclassified_ASV_1" &
                  figdf$col_genus != "Unclassified_ASV_9" &
                  figdf$col_genus != "Aquarickettsia" &
                  figdf$col_genus != "Ruegeria" &
                  figdf$col_genus != "Candidatus Tenderia" &
                  figdf$col_genus != "Haloferula" &
                  figdf$col_genus != "Family_Hyphomonadaceae" &
                  figdf$col_genus != "P3OB-42"] <- NA

#phosphate
figdf$col_genus[figdf$col_genus != "Unclassified_ASV_5" &
                  figdf$col_genus != "Ferrimonas"] <- NA

levels(figdf$col_genus)
#add new factor
figdf$col_genus <- factor(figdf$col_genus, levels = c(levels(figdf$col_genus), "Other"))
# convert NAs to other
figdf$col_genus[is.na(figdf$col_genus)] = "Other"

g7_nitrate_colors = c('#660066', '#51DA9B', '#37A337', '#848482')
g7_phosphate_colors = c('#AA2384',"#ebc650",'#848482')
g7_ammonium_colors = c('#FF5C19', '#848482')
g7_combined_colors = c('#AA2384','#A1CAF1', '#660066', '#848482')

# Taxa above the dashed line are significant with the null‐hypothesis rejected 80% of the time.
ggfig <- ggplot(figdf, aes(x = x, y = y, color = col_genus)) +
  geom_vline(xintercept = 0, color = "grey") +
  geom_point(size = 2) +
  facet_grid(~group) +
  theme_bw() +
  ylab("W statistic") +
  xlim(-3.1, 2) +
  theme(axis.title.x=element_blank()) +
  theme(legend.position = "none") +
  scale_color_manual(name = "Genus", values = g7_phosphate_colors) + #change to color set from above matching contrast of interest
  geom_hline(yintercept = cut_off["detected_0.7"], linetype = "dashed") +
  annotate("text", label = "W = 0.7", size = 3.5, x = -1.5, y = 22.5, color = "blue") #changed y location of text for each contrast
ggfig

ggfig2 <- ggplot(subset(figdf, group == "6 vs. 3 weeks"), aes(x = x, y = y, color = col_genus)) +
   geom_vline(xintercept = 0, color = "grey") +
   geom_point(size = 2) +
   theme_bw() +
   facet_grid(~group) +
   theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
   theme(axis.title.x=element_blank()) +
   xlim(-3, 2) +
   scale_color_manual(name = "Genus", values = g7_phosphate_colors) +
   geom_hline(yintercept = cut_off["detected_0.7"], linetype = "dashed") + 
   annotate("text", label = "W = 0.7", size = 3.5, x = -1.5, y = 22.5, color = "blue")
ggfig2


merged <- ggarrange(ggfig + theme(plot.margin = margin(r = 1)), ggfig2 + theme(plot.margin = margin(l = 1)), nrow = 1, align = "h", widths = c(0.9,1))
merged <- annotate_figure(merged,
                bottom = text_grob("CLR Mean Difference", hjust = 1, color = "black", size = 10))
merged

#panels for figure 5
ggsave(filename="~/Mote_nutrient_experiment/genotype_7/plots/ANCOM_7_nitrate.svg", plot=merged, device="svg",  width = 6.5, height = 3, dpi=500)
ggsave(filename="~/Mote_nutrient_experiment/genotype_7/plots/ANCOM_7_phosphate.svg", plot=merged, device="svg",  width = 5.5, height = 2.95, dpi=500)
ggsave(filename="~/Mote_nutrient_experiment/genotype_7/plots/ANCOM_7_ammonium.png", plot=merged, device="png",  width = 5, height = 3, dpi=500)
ggsave(filename="~/Mote_nutrient_experiment/genotype_7/plots/ANCOM_7_combined.svg", plot=merged, device="svg",  width = 6.5, height = 3, dpi=500)


#no treatment vs nutrient figure, figure 7
g7_ctrl_colors = c("#F3C300", '#848482')

figdf <- (subset(figdf, group == "No Treatment, 6 vs. 0 weeks"))
ggfig <- ggplot(figdf, aes(x = x, y = y, color = col_genus)) +
  geom_vline(xintercept = 0, color = "grey") +
  geom_point(size = 2) +
  facet_grid(~group) +
  theme_bw() +
  ylab("W statistic") +
  xlim(-3, 2) +
  ylim(0,34) +
  theme(axis.title.x=element_blank()) +
  theme(legend.position = "none") +
  scale_color_manual(name = "Genus", values = g7_ctrl_colors) +
  geom_hline(yintercept = cut_off["detected_0.7"], linetype = "dashed") +
  annotate("text", label = "W = 0.7", size = 3.5, x = -1.5, y = 24, color = "blue")
ggfig
ggsave(filename="~/Mote_nutrient_experiment/genotype_7/plots/ANCOM_7_no_trt.svg", plot=ggfig, device="svg",  width = 6.5, height = 3, dpi=500)


g7_nut_colors = c('#AA2384', '#1c9e23',  '#FF5C19', '#A1CAF1','#660066', '#FCA254', '#51DA9B','#0067A5','#58B8EC',"#FB2F2F",  '#848482')
#("#660066", "#0072B2", "#CC3300","#FFFF00", "#58B8EC", "#8CD09F", "#AA2384","#56B4E9", "#E27AB7")
figdf <- (subset(figdf, group == "Nutrient, 6 vs. 0 weeks"))
ggfig2 <- ggplot(figdf, aes(x = x, y = y, color = col_genus)) +
  geom_vline(xintercept = 0, color = "grey") +
  geom_point(size = 2) +
  facet_grid(~group) +
  theme_bw() +
  #ylab("W statistic") +
  xlim(-3, 2) +
  ylim(0,34) +
  theme(axis.title.x=element_blank()) +
  #theme(legend.position = "none") +
  scale_color_manual(name = "Genus", values = g7_nut_colors) +
  geom_hline(yintercept = cut_off["detected_0.7"], linetype = "dashed") +
  annotate("text", label = "W = 0.7", size = 3.5, x = -1.5, y = 25, color = "blue")
ggfig2

merged <- ggarrange(ggfig + theme(plot.margin = margin(r = 1)), ggfig2 + theme(plot.margin = margin(l = 1)), nrow = 1, align = "h", widths = c(0.9,2))
merged <- annotate_figure(merged,
                          bottom = text_grob("CLR Mean Difference", hjust = 1, color = "black", size = 10))
merged

ggsave(filename="~/Mote_nutrient_experiment/genotype_7/plots/ANCOM_7_nutrient_control.svg", plot=merged, device="svg",  width = 6.5, height = 3, dpi=500)

