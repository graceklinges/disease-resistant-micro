library(phyloseq)
library(multcomp)
library(picante)
library(ggplot2 )
library(tidyverse)
library(plyr)
library(dplyr)

rm(list=ls())
setwd("~/Mote_nutrient_experiment/genotype_7")

load(file = "ps_g7.RData") #renamed, filtered, unrarefied ps object

mapfile = "~/Mote_nutrient_experiment/data/Mapping-file-full-renamed.txt"
map = import_qiime_sample_data(mapfile)
sample_data(ps) <- map

#functions
asinTransform <- function(p) { asin(sqrt(p)) }

# Estimate faith's phylogenetic diversity 
estimate_pd <- function(phylo){
  
  # Error if input is not of class phylo
  if(class(phylo) != "phyloseq"){
    stop("Input file is not of class 'phyloseq'.")
  }
  
  # Error if no class phy_tree
  if(!(.hasSlot(phylo, "phy_tree"))){
    stop("Could not find tree slot in phylo object.")
  }
  
  # Transpose if needed
  # Adapted from phyloseq/vegan import
  OTU <- phyloseq::otu_table(phylo)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  
  # Get matrix version of OTU table
  otutable <- as(OTU, "matrix")
  
  # Get phylogenetic tree from pyloseq object
  tree <- phyloseq::phy_tree(phylo)
  
  # Print status message
  message("Calculating Faiths PD-index...")
  
  # If object is greater than 10mb, then print status message
  if(object.size(otutable) > 10000000){
    message("This is a large object, it may take awhile...")
  }
  
  # Calculate Faith's PD-index
  pdtable <- picante::pd(otutable, tree, include.root = F)
  
  # Return data frame of results
  return(pdtable)
}

# Calculate standard error
sderr <- function(x) {sd(x)/sqrt(length(x))}
#+++++++++++++++++++++++++
# Function to calculate the mean and the standard deviation
# for each group
#+++++++++++++++++++++++++
# data : a data frame
# varname : the name of a column containing the variable
#to be summariezed
# groupnames : vector of column names to be used as
# grouping variables
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sderr(x[[col]]), na.rm=TRUE)
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

normality.plots <- function(x) {
  par(mfrow=c(2,2))
  hist(residuals(x), main = "Histogram", xlab = "Values")
  boxplot(residuals(x), main = "Boxplot", ylab = "Values")
  qqnorm(residuals(x))
  qqline(residuals(x))
  plot(density(residuals(x)), main = "Kernel Density Estimate")
}

#Initialize matrices to store alpha diversity estimates
nsamp = nsamples(ps)
observed <- matrix(nrow = nsamp)
row.names(observed) <- sample_names(ps)
simpson <- matrix(nrow =nsamp)
row.names(simpson) <- sample_names(ps)
shannon <- matrix(nrow =nsamp)
row.names(shannon) <- sample_names(ps)

#Calculate statistics
# Options for measures = ("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher")
# Calculate observed
obs <- as.numeric(as.matrix(estimate_richness(ps, measures = "Observed")))
observed[ ,] <- obs
colnames(observed) [1] <- "observed"
# Calculate simpson
simp <- as.numeric(as.matrix(estimate_richness(ps, measures = "Simpson")))
simpson[ ,] <- simp
colnames(simpson) [1] <- "simpson"
# Calculate shannon
shan <- as.numeric(as.matrix(estimate_richness(ps, measures = "Shannon")))
shannon[ ,] <- shan
colnames(shannon) [1] <- "shannon"

#Combine our estimates for observed, simpson, and Shannon indices into one dataframe
alpha <- cbind(observed,simpson,shannon)
head(alpha)
# Add the sample metadata into this dataframe
s <- data.frame(sample_data(ps))
s$Exposure_weeks <- as.factor(s$Exposure_weeks)
levels(s$Exposure_weeks) <- c("0","3", "6")
alphadiv <- cbind(alpha, s)
head(alphadiv)
alphadiv <- alphadiv[,-5]
head(alphadiv)

write.csv(alphadiv, file = "./alphadiv_2000_rare_geno7.csv")

# #Kruskal Wallis -- Shannon's
kruskal.test(shannon ~ Exposure_weeks, data = alphadiv) #sig
kruskal.test(shannon ~ NW_T0_merge, data = alphadiv) #sig
kruskal.test(shannon ~ Nutrient_no_level, data = alphadiv) #sig
kruskal.test(shannon ~ No_level_weeks_noT0, data = alphadiv) #sig
kruskal.test(shannon ~ Nutrient, data = alphadiv) #sig
kruskal.test(shannon ~ No_level_no_T0, data = alphadiv) #NS
kruskal.test(shannon ~ group_weeks, data = alphadiv)

pairwise.wilcox.test(alphadiv$shannon, alphadiv$Exposure_weeks, p.adjust.method = 'fdr')# 0 to 3 NS

pairwise.wilcox.test(alphadiv$shannon, alphadiv$NW_T0_merge, p.adjust.method = 'fdr') #same


#Only ammonium and combined significantly different from T0
KW_shannon_Nutrient_no_level <- pairwise.wilcox.test(alphadiv$shannon, alphadiv$Nutrient_no_level, p.adjust.method = 'fdr')

#most NS except comparisons to T0
KW_shannon_NW_T0_merge <- pairwise.wilcox.test(alphadiv$shannon, alphadiv$NW_T0_merge, p.adjust.method = 'fdr')

#for plot
KW_shannon_No_level_weeks_noT0 <- pairwise.wilcox.test(alphadiv$shannon, alphadiv$No_level_weeks_noT0, p.adjust.method = 'fdr')
write.csv(KW_shannon_No_level_weeks_noT0[["p.value"]], file = "KW_shannon_No_level_weeks_noT0.csv")

#Figure 2
myCols <- c("#0072B2", "#CC3300", "#03B25C")
#goes with stats from KW_shannon_No_level_weeks_noT0, use alpha_div_sig_letters to get letters
A <- ggplot(alphadiv, aes(x=Exposure_weeks, y=shannon)) +
  facet_wrap(~No_level_no_T0)+
  geom_boxplot(outlier.shape = NA, color = "gray35") +
  geom_point(aes(color = Exposure_weeks), 
             position = position_jitter(width = .25, height = 0)) +
  ylab("Shannon Diversity Index") +
  theme_bw() +
  theme(aspect.ratio = 1.8) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = c(0.9,0.3),
        legend.title = element_text(color = "black", size = 12),
        legend.text = element_text(color = "black", size = 12)) +
  scale_colour_manual(values = myCols) +
  stat_summary(geom = 'text', label = c("a","abc","bc","ab", "ab", "c", "a", "ab", "abc", "ab", "abc", "abc","a", "a", "abc"), fun = max, vjust = -1, size = 3.5) + #stats with rick
    ylim(0,6) +
  ggtitle("Shannon Diversity by Nutrient Weeks") +
theme(plot.title = element_text(hjust = 0.5))
A

#Supp figure 3
kw_group_weeks <- pairwise.wilcox.test(alphadiv$shannon, alphadiv$group_weeks, p.adjust.method = 'fdr')
write.csv(kw_group_weeks[["p.value"]], file = "kw_group_weeks.csv")

myCols <- c("#0072B2", "#CC3300", "#03B25C")
#goes with stats from kw_group_weeks, use alpha_div_sig_letters to get letters
A <- ggplot(alphadiv, aes(x=Exposure_weeks, y=shannon)) +
  facet_wrap(~Group)+
  geom_boxplot(outlier.shape = NA, color = "gray35") +
  geom_point(aes(color = Exposure_weeks), 
             position = position_jitter(width = .25, height = 0)) +
  ylab("Shannon Diversity Index") +
  theme_bw() +
  theme(aspect.ratio = 1.8) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = c(0.9,0.9),
        legend.title = element_text(color = "black", size = 12),
        legend.text = element_text(color = "black", size = 12)) +
  scale_colour_manual(values = myCols) +
  stat_summary(geom = 'text', label = c("a", "ab", "ab","a", "a", "b"), fun = max, vjust = -1, size = 3.5) + #letters from alpha div letters script
  ylim(0,6) +
  ggtitle("Shannon Diversity by Nutrient-Treated and No Treatment Groups") +
  theme(plot.title = element_text(hjust = 0.5))
A
ggsave(filename="~/Mote_nutrient_experiment/genotype_7/plots/group_weeks.png", plot=A, device="png", dpi=500,  width = 8, height = 9,)
