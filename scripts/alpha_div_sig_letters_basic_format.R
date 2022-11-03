library(multcompView)

tri.to.squ<-function(x)
{
  rn<-row.names(x)
  cn<-colnames(x)
  an<-unique(c(cn,rn))
  myval<-x[!is.na(x)]
  mymat<-matrix(1,nrow=length(an),ncol=length(an),dimnames=list(an,an))
  for(ext in 1:length(cn))
  {
    for(int in 1:length(rn))
    {
      if(is.na(x[row.names(x)==rn[int],colnames(x)==cn[ext]])) next
      mymat[row.names(mymat)==rn[int],colnames(mymat)==cn[ext]]<-x[row.names(x)==rn[int],colnames(x)==cn[ext]]
      mymat[row.names(mymat)==cn[ext],colnames(mymat)==rn[int]]<-x[row.names(x)==rn[int],colnames(x)==cn[ext]]
    }
    
  }
  return(mymat)
}

kw_group_weeks <- pairwise.wilcox.test(alphadiv$shannon, alphadiv$group_weeks, p.adjust.method = 'fdr')
mymat<-tri.to.squ(kw_group_weeks[["p.value"]])
myletters<-multcompLetters(mymat,compare="<=",threshold=0.05,Letters=letters)

KW_simpsons_NW_T0_merge <- pairwise.wilcox.test(alphadiv$simpson, alphadiv$NW_T0_merge, p.adjust.method = 'fdr')
mymat<-tri.to.squ(KW_simpsons_NW_T0_merge[["p.value"]])
myletters<-multcompLetters(mymat,compare="<=",threshold=0.05,Letters=letters)

KW_simpsons_No_level_weeks <- pairwise.wilcox.test(alphadiv$simpson, alphadiv$No_level_weeks, p.adjust.method = 'fdr')
mymat<-tri.to.squ(KW_simpsons_No_level_weeks[["p.value"]])
myletters<-multcompLetters(mymat,compare="<=",threshold=0.05,Letters=letters)

KW_shannon_No_level_weeks_noT0 <- pairwise.wilcox.test(alphadiv$shannon, alphadiv$No_level_weeks_noT0, p.adjust.method = 'fdr')
mymat<-tri.to.squ(KW_shannon_No_level_weeks_noT0[["p.value"]])
myletters<-multcompLetters(mymat,compare="<=",threshold=0.05,Letters=letters)

tle_stats <- read.csv("tle_matrix.csv", header = TRUE, row.names = 1)
myletters<-multcompLetters(tle_stats,compare="<=",threshold=0.05,Letters=letters)
