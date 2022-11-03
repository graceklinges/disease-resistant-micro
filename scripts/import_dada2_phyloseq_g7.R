library("dada2")
library("seqinr")
library("biomformat")
library("data.table")
library("ggplot2")
library("ggthemes")
library("ampvis2")
library("cowplot")
library("phyloseq")

rm(list=ls())
setwd("~/Mote_nutrient_experiment/genotype_7")

# functions
fast_melt = function(physeq){
  # supports "naked" otu_table as `physeq` input.
  otutab = as(otu_table(physeq), "matrix")
  if(!taxa_are_rows(physeq)){otutab <- t(otutab)}
  otudt = data.table(otutab, keep.rownames = TRUE)
  setnames(otudt, "rn", "taxaID")
  # Enforce character taxaID key
  otudt[, taxaIDchar := as.character(taxaID)]
  otudt[, taxaID := NULL]
  setnames(otudt, "taxaIDchar", "taxaID")
  # Melt count table
  mdt = melt.data.table(otudt, 
                        id.vars = "taxaID",
                        variable.name = "SampleID",
                        value.name = "count")
  # Remove zeroes, NAs
  mdt <- mdt[count > 0][!is.na(count)]
  # Calculate relative abundance
  mdt[, RelativeAbundance := count / sum(count), by = SampleID]
  if(!is.null(tax_table(physeq, errorIfNULL = FALSE))){
    # If there is a tax_table, join with it. Otherwise, skip this join.
    taxdt = data.table(as(tax_table(physeq, errorIfNULL = TRUE), "matrix"), keep.rownames = TRUE)
    setnames(taxdt, "rn", "taxaID")
    # Enforce character taxaID key
    taxdt[, taxaIDchar := as.character(taxaID)]
    taxdt[, taxaID := NULL]
    setnames(taxdt, "taxaIDchar", "taxaID")
    # Join with tax table
    setkey(taxdt, "taxaID")
    setkey(mdt, "taxaID")
    mdt <- taxdt[mdt]
  }
  return(mdt)
}

# phyloseq object output from Dada2
load(file = "~/Mote_nutrient_experiment/data/ps_full.RData")
seqtab_7 <- readRDS("~/Mote_nutrient_experiment/genotype_7/seqtab_7.rds")
nsamples(ps)
# sequences object made from Dada2

# import metadata and merge into phyloseq object
mapfile = "~/Mote_nutrient_experiment/data/Mapping-file-full-renamed.txt"
map = import_qiime_sample_data(mapfile)
sample_data(ps) <- map

g7 <- subset_samples(ps, Genotype == "7")

# export taxonomy to import into qiime2
tax<-as(tax_table(ps),"matrix")
tax_cols <- c("Kingdom", "Phylum", "Class", "Order","Family","Genus", "Species")
tax<-as.data.frame(tax)
tax$taxonomy<-do.call(paste, c(tax[tax_cols], sep=";"))
for(co in tax_cols) tax[co]<-NULL
#write.table(tax, "taxonomy_7.txt", quote=FALSE, col.names=FALSE, sep="\t")

# summary of data
ps
#summary(sample_data(ps))
ps_full <- ps

ntaxa(ps) #3941
nsamples(ps) #160
rank_names(ps)
sample_names(ps)[1:5]
sample_variables(ps)

# remove mitochondria and chloroplasts, is.na important becuase if not included
# this command will also remove all Family = NA or Order = NA
ps_with_mito = subset_taxa(ps, (Order!="Chloroplast") | is.na(Order)) #lost 171 taxa
ps_no_mito = subset_taxa(ps_with_mito, (Family!="Mitochondria") | is.na(Family)) #lost 117 taxa
ps_no_Eukaryota = subset_taxa(ps_no_mito, (Kingdom!="Eukaryota") | is.na(Kingdom)) #lost 0 taxa

sum(sample_sums(ps)) - sum(sample_sums(ps_no_Eukaryota)) #1,246,196 reads
ntaxa(ps_no_Eukaryota) #3653
ps = ps_no_Eukaryota #lost 288 taxa total
ntaxa(ps)
ps <- subset_samples(ps, sample_sums(ps) > 1000) 
ntaxa(ps)


summary(sample_sums(ps))
#mean was 13932, min was 268 

nsamples(ps) #6 lost

summary(taxa_sums(ps)) # 6 was first quantile, mean was 610.2
# Filter on prevalence or total counts
pst = fast_melt(ps)
prevdt = pst[, list(Prevalence = sum(count > 0), 
                    TotalCounts = sum(count)),
             by = taxaID]
uniques = prevdt[(Prevalence <2 & TotalCounts >2), taxaID] #2398 taxa present in only one sample
keepTaxa = prevdt[(Prevalence >=0 & TotalCounts >0), taxaID]
ps = prune_taxa(keepTaxa,ps)
ntaxa(ps)
summary(taxa_sums(ps))
save(ps, file = "ps_unpruned_geno7.RData")

keepTaxa2 = prevdt[(Prevalence >=0 & TotalCounts >6), taxaID]
ps_pruned= prune_taxa(keepTaxa2,ps)
summary(taxa_sums(ps_pruned)) #new min is 7, new mean is 855.2
ntaxa(ps_pruned) #2598

write.table(prevdt, file = "genotype7_prevdt.txt", sep = "\t")
write.table(pst, file = "genotype7_pst.txt", sep = "\t")

summary(sample_sums(ps_pruned)) #mean went down to 14427, median 13729, min 1018
sum(sample_sums(ps_pruned)) #2221792
sum(sample_sums(ps)) - sum(sample_sums(ps_pruned)) #3975 reads lost by pruning

ps <- ps_pruned
save(ps, file = "ps_pruned_geno7.RData")

#this summary doesn't like that I removed samples so I'm making two separate tables
summary_tab <- data.frame(init=(as.matrix(sample_sums(ps_full)))[,1],
                           chloros_removed=(as.matrix(sample_sums(ps_with_mito)))[,1],
                           mitos_removed=(as.matrix(sample_sums(ps_no_mito)))[,1])
write.table(summary_tab, "reads_lost_g7_pre_rem.txt", quote=FALSE, col.names=TRUE, sep="\t")

summary_tab <- data.frame(init=(as.matrix(sample_sums(ps)))[,1],
                          pruned=(as.matrix(sample_sums(ps_pruned)))[,1])
write.table(summary_tab, "reads_lost_g7_post_rem.txt", quote=FALSE, col.names=TRUE, sep="\t")

