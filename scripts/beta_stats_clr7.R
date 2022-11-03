library('phyloseq'); packageVersion('phyloseq')
library('vegan'); packageVersion('vegan')
library("usedist")
library("ggplot2")
library("microbiome")
library(multcompView)

rm(list=ls())
setwd("~/Mote_nutrient_experiment/genotype_7")

pairwise.adonis.dm <- function(x,factors,stratum=NULL,p.adjust.m="fdr",perm=999){
  
  library(vegan)
  if(class(x) != "dist") stop("x must be a dissimilarity matrix (dist object)")
  co = as.matrix(combn(unique(factors),2))
  pairs = c()
  F.Model =c()
  R2 = c()
  p.value = c()
  
  
  
  for(elem in 1:ncol(co)){
    sub_inds <- factors %in% c(as.character(co[1,elem]),as.character(co[2,elem]))
    resp <- as.matrix(x)[sub_inds,sub_inds]
    ad = adonis(as.dist(resp) ~
                  
                  factors[sub_inds], strata=stratum[sub_inds], permutations=perm);
    
    pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    F.Model =c(F.Model,ad$aov.tab[1,4]);
    R2 = c(R2,ad$aov.tab[1,5]);
    p.value = c(p.value,ad$aov.tab[1,6])
    
  }
  
  p.adjusted = p.adjust(p.value,method=p.adjust.m)
  pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted)
  return(pairw.res)
}


load(file = "ps_clr_g7.RData") #renamed, clr transformed data
sample_sums(clr)
nsamples(clr) #154

mapfile = "~/Mote_nutrient_experiment/data/Mapping-file-full-renamed.txt" #confirm correct mapping file associated
map = import_qiime_sample_data(mapfile)
sample_data(clr) <- map

clr@sam_data$Nutrient_no_level<- factor(clr@sam_data$Nutrient_no_level, c("T0", "No Treatment", "Nitrate", "Phosphate", "Ammonium", "Combined")) 
clr@sam_data$Exposure_weeks <- as.factor(clr@sam_data$Exposure_weeks)

### Calculate distance matrices
# only calculating euclidean distance for clr transformed data because we have negative values, and with the clr transform we've moved our data into "real space"
clr_euc <- phyloseq::distance(clr, method = "euclidean")

### PERMANOVAs

# PERMANOVA's with Adonis -  Permutational Multivariate Analysis of Variance
# Make a data frame from the sample_data
sampledf <- data.frame(sample_data(clr))
# Adonis 
adonis(clr_euc ~ Exposure_weeks*Nutrient_no_level, data = sampledf) #interaction not significant

pairwise.adonis.dm(clr_euc, sample_data(clr)$Nutrient_no_level, p.adjust.m = "fdr")
#comparisons against T0 significant, nitrate vs ammonium sig, ammonium vs NT
#
pairwise.adonis.dm(clr_euc, sample_data(clr)$Level, p.adjust.m = "fdr")
#H vs L NS

pairwise.adonis.dm(clr_euc, sample_data(clr)$Exposure_weeks, p.adjust.m = "fdr")
# 0 vs 3 and 0 vs 6 sig, not 3 vs 6

adonis(clr_euc ~ No_level_weeks, data = sampledf) #Df = 10
NLweeks <- pairwise.adonis.dm(clr_euc, sample_data(clr)$No_level_weeks, p.adjust.m = "fdr")
write.csv(NLweeks, "NLweeks.csv")

pairwise.adonis.dm(clr_euc, sample_data(clr)$Level, p.adjust.m = "fdr")
#0 vs H and vs L sig, NS H vs L


#PCA ordination with no outliers
ord_clr <- phyloseq::ordinate(clr, "RDA") #RDA without constraints = PCA
phyloseq::plot_scree(ord_clr) + 
  geom_bar(stat="identity", fill = "blue") +
  labs(x = "\nAxis", y = "Proportion of Variance\n") #screen plot plateaus quickly
head(ord_clr$CA$eig)
sapply(ord_clr$CA$eig[1:5], function(x) x / sum(ord_clr$CA$eig))    

clr@sam_data$Exposure_weeks <- as.factor(clr@sam_data$Exposure_weeks)
#clr PCA figure 3A
clr1 <- ord_clr$CA$eig[1] / sum(ord_clr$CA$eig)
clr2 <- ord_clr$CA$eig[2] / sum(ord_clr$CA$eig)
PCA <- phyloseq::plot_ordination(clr, ord_clr, type="samples", color="Nutrient_no_level", shape="Exposure_weeks") + 
  geom_point(size = 2) +
  coord_fixed(((clr2 / clr1))*2) +
  stat_ellipse(aes(group = Nutrient_no_level), linetype = 1) +
  ggtitle("Principal Component Analysis") +
  labs(color = "Treatment") +
  labs(shape = "Exposure Weeks") +
  scale_color_manual(values = c("#BE0032", "#009E73","#F3C300","#848482", "#0067A5", "#000000"))
PCA
ggsave(filename="~/Mote_nutrient_experiment/genotype_7/plots/PCA_7_clr.png", plot=PCA, device="png", dpi=700)

#Dispersion
### PERMDISPs
anova(betadisper(clr_euc, sampledf$Nutrient, bias.adjust = TRUE)) #sig 0.008763
anova(betadisper(clr_euc, sampledf$Nutrient_no_level, bias.adjust = TRUE)) #sig 0.007797
anova(betadisper(clr_euc, sampledf$Exposure_weeks, bias.adjust = TRUE)) #sig 0.004982
anova(betadisper(clr_euc, sampledf$No_level_weeks, bias.adjust = TRUE)) #sig 0.001409

p.adjust(permutest(betadisper(clr_euc, sampledf$Exposure_weeks, 
                              bias.adjust = TRUE), 
                   pairwise=TRUE)$pairwise$permuted, method = 'fdr')
#0-3 NS, 0-6 N, 3-6 S

NLweeks_disp <- p.adjust(permutest(betadisper(clr_euc, sampledf$No_level_weeks_noT0, 
                                              bias.adjust = TRUE), 
                                   pairwise=TRUE)$pairwise$permuted, method = 'fdr')
write.csv(NLweeks_disp, "NLweeks_disp.csv")
#only comparison within group that is significant is A0-A6

NL_disp <- p.adjust(permutest(betadisper(clr_euc, sampledf$Nutrient_no_level, 
                                              bias.adjust = TRUE), 
                                   pairwise=TRUE)$pairwise$permuted, method = 'fdr')
write.csv(NL_disp, "NL_disp.csv")
#ammonium lower dispersion than nitrate or T0

plot(dispr, main = "Ordination Centroids and Dispersion Labeled: Aitchison Distance", sub = "")
boxplot(dispr, main = "", xlab = "") #quick dispersion plot
permutest(dispr)

# extract distance to centroid
clr@sam_data$Exposure_weeks <- as.factor(clr@sam_data$Exposure_weeks)
sampledf <- data.frame(sample_data(clr))
disp <- betadisper(clr_euc, sampledf$No_level_weeks, bias.adjust = TRUE)
dispd <- as.data.frame(disp$distances)
dispd <- cbind(dispd, sample_data(clr))
colnames(dispd)[1] <- "distance"


#Supp Fig. 3
dispersion_stats <- p.adjust(permutest(betadisper(clr_euc, sampledf$No_level_weeks_noT0, 
                                                  bias.adjust = TRUE), 
                                       pairwise=TRUE)$pairwise$permuted, method = 'fdr')
write.table(dispersion_stats, file = "dispersion_stats_clr2.txt", sep = "\t")
#manually turned this into a triangular matrix
dispersion_table <- read.table("dispersion_stats_matrix.txt", header = TRUE, sep = "\t", row.names = 1)
myletters<-multcompLetters(dispersion_table,compare="<=",threshold=0.05,Letters=letters)

dispd$Exposure_weeks <- as.factor(dispd$Exposure_weeks)
levels(dispd$Exposure_weeks) <- c("0","3", "6")
myCols <- c("#0072B2", "#CC3300", "#03B25C")

dispersion_plot <- ggplot(dispd, aes(x=Exposure_weeks, y=distance)) +
  facet_wrap(~No_level_no_T0)+
  geom_boxplot(outlier.shape = NA, color = "gray35") +
  geom_point(aes(color = Exposure_weeks), 
             position = position_jitter(width = .25, height = 0)) +
  ylab("Distance-to-centroid") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = c(0.85, 0.3)) +
  labs(color = "Exposure Weeks") +
  scale_colour_manual(values = myCols) +
  stat_summary(geom = 'text', label = c("a","ab","b", "ab", "a", "ab", "a", "a", "ab", "a", "ab", "ab", "a", "a", "ab"), fun = max, vjust = -0.5, size = 4) +
  ggtitle("Dispersion Over Time, All Treatments") +
  scale_y_continuous(expand = expansion(mult = c(.1)))
dispersion_plot
ggsave(filename="~/Mote_nutrient_experiment/genotype_7/plots/Dispersion_all_clr.png", plot=dispersion_plot, device="png", height = 5, width = 5, dpi=500)

clr_euc <- phyloseq::distance(clr, method = "euclidean")
sampledf <- data.frame(sample_data(clr))
disp <- betadisper(clr_euc, sampledf$Nutrient_no_level, bias.adjust = TRUE)
dispd <- as.data.frame(disp$distances)
dispd <- cbind(dispd, sample_data(clr))
dispd$Nutrient_no_level <- factor(dispd$Nutrient_no_level, c("T0", "No Treatment", "Nitrate", "Phosphate", "Ammonium", "Combined")) 

colnames(dispd)[1] <- "distance"

dispersion_stats <- p.adjust(permutest(betadisper(clr_euc, sampledf$Nutrient_no_level, 
                                                  bias.adjust = TRUE), 
                                       pairwise=TRUE)$pairwise$permuted, method = 'fdr')

#Figure 3B
colors = c("#BE0032", "#009E73","#F3C300","#848482", "#0067A5", "#000000")
dispersion_plot <- ggplot(dispd, aes(x=Nutrient_no_level, y=distance)) +
  #facet_wrap(~Nutrient)+
  geom_boxplot(outlier.shape = NA, color = "gray35") +
  geom_point(aes(color = Nutrient_no_level, shape = Exposure_weeks), 
             position = position_jitter(width = .25, height = 0)) +
  ylab("Distance-to-centroid") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = 'none')  +
  labs(color = "Treatment") +
  scale_colour_manual(values = colors) +
  stat_summary(geom = 'text', label = c("a","ab", "b", "ab", "ab", "b"), fun = min, vjust = 3, size = 4) +
  ggtitle("Dispersion, All Treatments") +
  scale_y_continuous(expand = expansion(mult = c(.1)))
dispersion_plot
ggsave(filename="~/Mote_nutrient_experiment/genotype_7/plots/Dispersion_nutrient.png", plot=dispersion_plot, device="png", height = 7, width = 5.5, dpi=500)


