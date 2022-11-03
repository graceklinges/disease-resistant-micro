setwd("~/Mote_nutrient_experiment/genotype_7")

library(ggplot2)
library(dunn.test)
library(ggpubr)
library(dplyr)
library(agricolae)

tle <- read.csv("tle_7.csv")
tle$Weeks <- factor(tle$Weeks, levels = c("3", "6")) 
tle$trt_weeks <- factor(tle$trt_weeks, levels = c("NT3", "NT6", "LA3" , "LA6", "HA3", "HA6", "LN3", "LN6", "HN3", "HN6", "LP3", "LP6", "HP3", "HP6", "LC3", "LC6", "HC3", "HC6")) 

#tle is not normally distributed
ggdensity(tle$tle)
shapiro.test(tle$tle) #fails shapiro test, low p value
#passes fligner-killeen test of homogeneity of variances for non-normally distributed data, but fails levene test
model <- lm(tle ~ trt_weeks, data = tle)

plot(model, which = 2) #q-q plot
plot(model, which = 3) #variance

#log transform
tle <- mutate(tle, log_tle = log(tle))
model <- lm(log_tle ~ trt_weeks, data = tle)
ggdensity(tle$log_tle) #looks better
shapiro.test(tle$log_tle) #passes
 
#modeling relationship of tle to treatment
model <- lm(log_tle ~ trt_weeks, data = tle)
summary(model)
anova(model)
#omnibus test shows main effect is significant

aov.model <- aov(log_tle ~ trt_weeks, data = tle)
summary(aov.model)
tukey_stats <- TukeyHSD(aov.model)
tukey_groups <- HSD.test(aov.model, "trt_weeks", group=TRUE)
tukey_groups <- tukey_groups$groups
tukey_stats <- as.data.frame(tukey_stats$trt_weeks)
write.table(tukey_stats, file = "Tukey_tle.txt", sep = "\t")
write.table(tukey_groups, file = "Tukey_groups.txt", sep = "\t")

colors = c('#848482', '#FCEDA2', '#FBD92E', '#E26060', '#BE0032', '#74A4D1', '#106EC8', '#70BD64','#1A9906')
colors_full = c('#848482', '#848482', '#FCEDA2','#FCEDA2', '#FBD92E', '#FBD92E', '#E26060', '#E26060', '#BE0032', '#BE0032', '#74A4D1', '#74A4D1', '#106EC8', '#106EC8', '#70BD64', '#70BD64', '#1A9906', '#1A9906')

tle$Treatment <- factor(tle$Treatment, levels = c("NT", "LA", "HA", "LN", "HN", "LP", "HP", "LC", "HC"))
boxplot <- ggplot(tle, aes(x=Treatment, y=tle, fill=trt_weeks)) +
  geom_boxplot() +
  facet_grid(~Weeks) +
  theme_bw() +
  scale_fill_manual(values = colors_full) +
  labs(x="Treatment Weeks", y = "Linear Extension (mm)") +
  stat_summary(geom = 'text', label = c("a", "ab", "a", "a", "a", "a", "a", "a", "a", "c", "c", "c", "abc", "ab", "c", "bc", "c", "bc"), fun = min, vjust = 2, size = 3.5) +
  theme(legend.position = "none")
ggsave(filename="~/Mote_nutrient_experiment/genotype_7/plots/tle_boxplot_full.svg", plot=boxplot, device="svg",  width = 6, height = 5, dpi=500)
boxplot


#treatment_weeks
#tle3 <- tle[ which(tle$Weeks=='3'),]
#tle6 <- tle[ which(tle$Weeks=='6'),]
#
# #3
# aov.model <- aov(log_tle ~ trt_weeks, data = tle3)
# summary(aov.model)
# tukey_stats <- TukeyHSD(aov.model)
# tukey_groups3 <- HSD.test(aov.model, "trt_weeks", group=TRUE)
# tukey_groups3 <- tukey_groups3$groups
# tukey_stats <- as.data.frame(tukey_stats$trt_weeks)
# write.table(tukey_stats, file = "Tukey_tle3.txt", sep = "\t")
# #6 
# aov.model <- aov(log_tle ~ trt_weeks, data = tle6)
# summary(aov.model)
# tukey_stats <- TukeyHSD(aov.model)
# tukey_groups6 <- HSD.test(aov.model, "trt_weeks", group=TRUE) #made this into pairwise matrix manually
# tukey_groups6 <- tukey_groups6$groups
# tukey_stats <- as.data.frame(tukey_stats$trt_weeks)
# write.table(tukey_stats, file = "Tukey_tle6.txt", sep = "\t") 
# 
# tle3$trt_weeks <- factor(tle3$trt_weeks, levels = c("NT3","LA3", "HA3", "LN3", "HN3", "LP3", "HP3", "LC3", "HC3")) 
# tle6$trt_weeks <- factor(tle6$trt_weeks, levels = c("NT6","LA6", "HA6", "LN6", "HN6", "LP6", "HP6", "LC6", "HC6")) 

# boxplot3 <- ggplot(tle3, aes(x=trt_weeks, y=tle, fill=trt_weeks)) + 
#   geom_boxplot() +
#   theme_bw() +
#   scale_fill_manual(values = colors) +
#   labs(x="Treatment Weeks", y = "Linear Extension (mm)") + 
#   ylim(0,14) +
#   stat_summary(geom = 'text', label = c("a","a", "a", "a", "a", "a", "a", "a", "a"), fun = min, vjust = 2, size = 3.5) +
#   theme(legend.position = "none")
# boxplot3
# 
# ggsave(filename="~/Mote_nutrient_experiment/genotype_7/plots/tle_boxplot3.png", plot=boxplot3, device="png",  width = 4, height = 4, dpi=500)
# ggsave(filename="~/Mote_nutrient_experiment/genotype_7/plots/tle_boxplot3.svg", plot=boxplot3, device="svg", dpi=500)
# 
# boxplot6 <- ggplot(tle6, aes(x=trt_weeks, y=tle, fill=trt_weeks)) + 
#   geom_boxplot() +
#   theme_bw() +
#   scale_fill_manual(values = colors) +
#   ylim(0,14) +
#   labs(x="Treatment Weeks", y = "Linear Extension (mm)") + 
#   stat_summary(geom = 'text', label = c("a","a", "ab", "ab", "b", "a", "ab", "ab", "ab"), fun = min, vjust = 2, size = 3.5) +
#   theme(legend.position = "none")
# boxplot6
# 
# ggsave(filename="~/Mote_nutrient_experiment/genotype_7/plots/tle_boxplot6.png", plot=boxplot6, device="png",  width = 4, height = 4, dpi=500)
# ggsave(filename="~/Mote_nutrient_experiment/genotype_7/plots/tle_boxplot6.svg", plot=boxplot6, device="svg", dpi=500)
