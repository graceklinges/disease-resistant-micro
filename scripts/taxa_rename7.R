library("phyloseq")
library("dada2")
library("seqinr")
library("biomformat")
library("data.table")
library("ggplot2")
library("ggthemes")
library("ampvis2")
library("cowplot")
library("phyloseq")

setwd("~/Mote_nutrient_experiment/genotype_7")

load(file = "ps_pruned_geno7.RData")
mapfile = "~/Mote_nutrient_experiment/data/Mapping-file-full-renamed.txt"
map = import_qiime_sample_data(mapfile)
sample_data(ps) <- map
ps
nsamples(ps) #154

#remove negatives
#ps = subset_samples(ps, SampleID != "HP-7-2-T2")
#ps = subset_samples(ps, SampleID != "MP-7-3-T2")

tax.clean <- data.frame(tax_table(ps))
write.table(tax.clean, "taxonomy_unrare_fullseq.txt", sep = "\t")

taxa_names(ps) <- paste0("ASV_", seq(ntaxa(ps)))
tax.clean <- data.frame(tax_table(ps))
write.table(tax.clean, "taxonomy_unrare.txt", sep = "\t")

new_names <- paste0("ASV_", seq(ntaxa(ps)))

for (i in 1:7){ tax.clean[,i] <- as.character(tax.clean[,i])}
tax.clean[is.na(tax.clean)] <- ""

for (i in 1:nrow(tax.clean)){
  
  if (tax.clean[i,6] == ""){
    family <- paste("Family_", tax.clean[i,5], sep = "")
    tax.clean[i, 6:7] <- family
  } else if (tax.clean[i,7] == ""){
    tax.clean$Species[i] <- paste("Genus",tax.clean$Genus[i], sep = "_")
  }
}

for (i in 1:nrow(tax.clean)){
  
  if (tax.clean[i,6] == "Family_"){
    family <- paste("", tax.clean[i,5], sep = "")
    tax.clean[i, 6:7] <- family
  }
  if (tax.clean[i,6] == "MD3-55"){
    tax.clean[i, 6:7] <- "Aquarickettsia"
  }
}

for (i in 1:nrow(tax.clean)){
  
  if (tax.clean[i,6] == ""){
    family <- paste0("Unclassified_",new_names[i], sep = "") 
    tax.clean[i, 5:7] <- family
  }
}


tax_table(ps) <- as.matrix(tax.clean)

save(ps, file = "ps_g7.RData") #unrarefied, rare taxa pruned, renamed

mapfile = "~/Mote_nutrient_experiment/data/Mapping-file-full-renamed.txt"
map = import_qiime_sample_data(mapfile)
sample_data(ps) <- map
ps

#remove negatives
ps = subset_samples(ps, SampleID != "Negative")
ps = subset_samples(ps, SampleID != "HP-7-2-T2")
ps = subset_samples(ps, SampleID != "MP-7-3-T2")
#155 samples

tax.clean <- data.frame(tax_table(ps))
#write.table(tax.clean, "taxonomy_unrare_fullseq.txt", sep = "\t")

taxa_names(ps) <- paste0("ASV_", seq(ntaxa(ps)))
tax.clean <- data.frame(tax_table(ps))
#write.table(tax.clean, "taxonomy_unrare.txt", sep = "\t")

new_names <- paste0("ASV_", seq(ntaxa(ps)))

for (i in 1:7){ tax.clean[,i] <- as.character(tax.clean[,i])}
tax.clean[is.na(tax.clean)] <- ""

for (i in 1:nrow(tax.clean)){
  
  if (tax.clean[i,6] == ""){
    family <- paste("Family_", tax.clean[i,5], sep = "")
    tax.clean[i, 6:7] <- family
  } else if (tax.clean[i,7] == ""){
    tax.clean$Species[i] <- paste("Genus",tax.clean$Genus[i], sep = "_")
  }
}

for (i in 1:nrow(tax.clean)){
  
  if (tax.clean[i,6] == "Family_"){
    family <- paste("", tax.clean[i,5], sep = "")
    tax.clean[i, 6:7] <- family
  }
  if (tax.clean[i,6] == "MD3-55"){
    tax.clean[i, 6:7] <- "Aquarickettsia"
  }
}

for (i in 1:nrow(tax.clean)){
  
  if (tax.clean[i,6] == ""){
    family <- paste0("Unclassified_",new_names[i], sep = "") 
    tax.clean[i, 5:7] <- family
  }
}


tax_table(ps) <- as.matrix(tax.clean)
save(ps, file = "ps_g7.RData") #unrarefied, pruned

pst = fast_melt(ps)
prevdt = pst[, list(Prevalence = sum(count > 0), 
                    TotalCounts = sum(count)),
             by = taxaID]

#write.table(prevdt, file = "genotype7_prevdt_rename.txt", sep = "\t")
#write.table(pst, file = "genotype7_pst_rename.txt", sep = "\t")

ps_rarefied <- rarefy_even_depth(ps, sample.size = 2000, rngseed = 999) #327 removed for pruned

#11 samples removedbecause they contained fewer reads than `sample.size`.
#Up to first five removed samples are: 
  
#LP-7-2-T1MN-7-1-T4MP-7-2-T4HA-7-3-T0LA-7-1-T4

ps_rarefied #2379 for pruned
sum(sample_sums(ps_rarefied)) #288000
sample_sums(ps_rarefied)
rarefied_names <- sample_names(ps_rarefied)
names <- sample_names(ps)
write.table(rarefied_names, file = "rarefied_names.txt", sep = "\t")
write.table(names, file = "full_names.txt", sep = "\t")

#CLR transform using microbiome package for beta div analyses
#CLR transform applies a pseudocount of min(relative abundance)/2 
#to exact zero relative abundance entries in OTU table before taking logs.
clr <- microbiome::transform(ps, 'clr')

#relative abundance transform
ps_rel = filter_taxa(ps, function(x) mean(x) > 0.1, TRUE)
ps_rel = transform_sample_counts(ps_rel, function(x) x / sum(x) )
ps_rel
#543 taxa and 144 samples

rel_abun <- ps_rel@otu_table
rel_abun <- as.matrix(rel_abun)
write.table(rel_abun, file = "rel_abun_new.txt", sep = "\t")

# save ps objects
save(ps_rel, file = "ps_rel_g7.RData") #rel abundance, pruned
save(clr, file = "ps_clr_g7.RData") #clr transformed, pruned

