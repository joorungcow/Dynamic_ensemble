library(rstatix)

setwd ('D:/jaeeun/experiment/symmetryEnsemble/symmetry_mean/exp1_2AFC/results')

data = read.csv('resultcsv.csv')

#ANOVAres <- anova_test(data, dv = data1, wid = sub1, within = c(pattern1, duration1)) 
# get_anova_table(ANOVAres)
#print(ANOVAres) # Print ANOVA table and  Mauchly's Test for Sphericity 

######outlier
data2 = read.csv('resultcsvOT.csv')

ANOVAresOt <- anova_test(data2, dv = data2, wid = sub2, within = c(pattern2, duration2))
print(ANOVAresOt)

ano <- data2 %>% anova_test(data2 ~ pattern2*duration2 + Error(sub2/(pattern2*duration2)))

aovRe <- aov(data2 ~ pattern2*duration2 + Error(sub2/(pattern2*duration2)), data2) #no Greenhouse-Geisser sphericity

data2 %>% group_by(pattern2, duration2) %>% shapiro_test(data2)


###perceived pattern?
#ANOVAres2 <- anova_test(data, dv = data1, wid = sub1, between = sawPattern1, within = c(pattern1, duration1)) 
#print(ANOVAres2)

#ANOVAresOt2 <- anova_test(data2, dv = data2, wid = sub2, between = sawPattern1, within = c(pattern2, duration2))
#print(ANOVAresOt2)


# ## different way to run ANOVA
# res.anova <- aov(value ~ bini*chlabel + Error(ERPset/(bini*chlabel)), data)
# anova_summary(res.anova, effect.size = "ges", detailed = FALSE, observed = NULL)