library(tidyverse)
library(ggpubr)
library(rstatix)
library(colorspace)

setwd ('D:/jaeeun/experiment/symmetryEnsemble/symmetry_mean/exp1_2AFC/results')

######-outlier
df = read.csv('resultcsvOT.csv')
df$pattern2 <- as.factor(df$pattern2)
df$duration2 <- as.factor(df$duration2)

#summary statistics
ss <- df %>% group_by(pattern2, duration2) %>% get_summary_stats(data2, type = "mean_sd")
ss
ss2 <- df %>% group_by(pattern2) %>% get_summary_stats(data2, type = "mean_sd")
ss2

#Outliers
df %>% group_by(pattern2, duration2) %>% identify_outliers(data2)

#Normality assumption
df %>% group_by(pattern2, duration2) %>% shapiro_test(data2)

ggqqplot(df, "data2", ggtheme = theme_bw()) + facet_grid(pattern2~duration2, labeller = "label_both")

#Computation
ANOVAresOt <- anova_test(df, dv = data2, wid = sub2, within = c(pattern2, duration2))
get_anova_table(ANOVAresOt)
#ano <- df %>% anova_test(data2 ~ pattern2*duration2 + Error(sub2/(pattern2*duration2)))

#Post-hoc tests
#comparisons for pattern variable
df %>% pairwise_t_test(
  data2 ~ pattern2, paired = TRUE,
  p.adjust.method = "bonferroni"
)
#comparisons for duration variable
df %>% pairwise_t_test(
  data2 ~ duration2, paired = TRUE,
  p.adjust.method = "bonferroni"
)

#Visualization
bxp <- ggboxplot(df, x = "duration2", y = "data2", color = "pattern2",  xlab = "duration",
                 ylab = "PSE", title = "PSE of mean size judgments", bxp.errorbar = TRUE,
) + labs(color="pattern") + scale_x_discrete(labels = c("250ms", "500ms")) + 
  scale_color_discrete_qualitative(labels = c("random", "symmetric"), palette = "Set2", order=c(2,1))
bxp


##perceived pattern?
#ANOVAresOt2 <- anova_test(data2, dv = data2, wid = sub2, between = sawPattern1, within = c(pattern2, duration2))
#print(ANOVAresOt2)

## different way to run ANOVA
#aovRe <- aov(data2 ~ pattern2*duration2 + Error(sub2/(pattern2*duration2)), df) #no Greenhouse-Geisser sphericity
#summary(aovRe, effect.size = "ges", detailed = FALSE, observed = NULL)