library(stringi)
library(dplyr)
library(ggplot2)
library(stringr)
library(purrr)
library(ggpubr)
library(reshape2)
library(ggthemes)
library(readxl)
library(rstatix)
library(cowplot)
library(ggsci)
library(DescTools)


mycurrentdirectory = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(mycurrentdirectory)

if (!dir.exists("../Datafiles")) {dir.create("../Datafiles")} else {}
if (!dir.exists("../Datafiles/Statistics")) {dir.create("../Datafiles/Statistics")} else {}

if (!dir.exists("../Figures")) {dir.create("../Figures")} else {}
if (!dir.exists("../Individualplots")) {dir.create("../Individualplots")} else {}

#Valine normalization step
list.files("../Datafiles")

y <- read_excel("../Datafiles//NC_Official_Data.xlsx",sheet = "log2VA")
y <- y %>%
  rename(Compound = Name)
colnames(y)
#now we'll melt the dataframe and create a column called organ min
zmelt<-melt(y,id.vars = "Compound")
zmelt <- zmelt %>%
  mutate(value = as.numeric(value),
         organ = stri_replace_all_regex(variable,"(.*)_(.*)_(.*)","$1"),
         time=stri_replace_all_regex(variable,"(.*)_(.*)_(.*)","$2")) %>%
  mutate(time=ifelse(time == "0min",0,
                     ifelse(time == "30min",30,
                            ifelse(time == "60min",60,
                                   ifelse(time == "120min",120,time))))) %>%
  mutate(time = as.numeric(time))

unique(zmelt$organ)
unique(zmelt$time)

###doing tiny project of only FDR correction by organ
stats.test.onesample.Group.median.perorganFDR<-data.frame()
organtypes<-unique(zmelt$organ)
for (i in 1:length(organtypes)){
  zmeltorgan<-zmelt %>% filter(organ==organtypes[i])
  stats.test.onesample.Group.median.each<-zmeltorgan %>% 
    group_by(Compound,organ) %>% 
    t_test(value~1,mu = 0,alternative = "two.sided") %>% 
    adjust_pvalue(method = "fdr") %>%
    mutate(organ = organtypes[i]) %>%
    add_significance("p.adj")
  stats.test.onesample.Group.median.perorganFDR<-rbind(stats.test.onesample.Group.median.perorganFDR,stats.test.onesample.Group.median.each)
}

zmeltsummary<-zmelt %>%
  group_by(Compound,organ) %>%
  summarise(Medianvalue = median(value,na.rm = T))

stats.test.onesample.Group.median.perorganFDR<-merge(zmeltsummary,stats.test.onesample.Group.median.perorganFDR)

write.csv(stats.test.onesample.Group.median.perorganFDR,file="../Datafiles/Statistics/1_OneSample_PerOrgan_NC.csv",row.names = FALSE)

###############Now we will conduct One-Sample Test for only fed

zmeltfed<-zmelt %>%
  filter(!time == 0)
stats.test.onesample.median<-zmeltfed %>% 
  group_by(Compound,organ) %>% 
  t_test(value~1,mu = 0,alternative = "two.sided") %>% 
  adjust_pvalue(method = "fdr") %>%
  add_significance("p.adj")


write.csv(stats.test.onesample.median,file="../Datafiles/Statistics/1_OneSample_NC_FedOnly.csv",row.names = FALSE)

###doing tiny project of only FDR correction by organ
stats.test.onesample.Group.median.perorganFDR<-data.frame()
organtypes<-unique(zmeltfed$organ)
for (i in 1:length(organtypes)){
  zmeltorgan<-zmeltfed %>% filter(organ==organtypes[i])
  stats.test.onesample.Group.median.each<-zmeltorgan %>% 
    group_by(Compound,organ) %>% 
    t_test(value~1,mu = 0,alternative = "two.sided") %>% 
    adjust_pvalue(method = "fdr") %>%
    mutate(organ = organtypes[i]) %>%
    add_significance("p.adj")
  stats.test.onesample.Group.median.perorganFDR<-rbind(stats.test.onesample.Group.median.perorganFDR,stats.test.onesample.Group.median.each)
}

write.csv(stats.test.onesample.Group.median.perorganFDR,file="../Datafiles/Statistics/1_OneSample_PerOrgan_NC_FedOnly.csv",row.names = FALSE)

###############Now we will conduct One-Sample Test per time
stats.test.onesample.median<-zmelt %>% 
  group_by(Compound,organ,time) %>% 
  t_test(value~1,mu = 0,alternative = "two.sided") %>% 
  adjust_pvalue(method = "fdr") %>%
  add_significance("p.adj")


write.csv(stats.test.onesample.median,file="../Datafiles/Statistics/1.1_OneSample_NC_AllMetabolites.csv",row.names = FALSE)



###doing  FDR correction by organ
stats.test.onesample.Group.median.perorganFDR<-data.frame()
organtypes<-unique(zmelt$organ)
for (i in 1:length(organtypes)){
  zmeltorgan<-zmelt %>% filter(organ==organtypes[i])
  stats.test.onesample.Group.median.each<-zmeltorgan %>% 
    group_by(Compound,organ,time) %>% 
    t_test(value~1,mu = 0,alternative = "two.sided") %>% 
    adjust_pvalue(method = "fdr") %>%
    mutate(organ = organtypes[i]) %>%
    add_significance("p.adj")
  stats.test.onesample.Group.median.perorganFDR<-rbind(stats.test.onesample.Group.median.perorganFDR,stats.test.onesample.Group.median.each)
}

zmelt2median <- zmelt %>%
  group_by(Compound,time,organ) %>%
  summarise(medianvalues = median(value))

zmelt2median <- zmelt2median %>%
  group_by(Compound,organ,time) %>%
  select(Compound,organ,time,#group2,#medianvaluesratioabsolutevalue,
         medianvalues) %>%
  mutate(releaseuptake = ifelse(medianvalues > 0, "release",
                                ifelse(medianvalues < 0, "uptake","none")))

stats.test.onesample.Group.median.perorganFDR2<-merge(stats.test.onesample.Group.median.perorganFDR,zmelt2median)

write.csv(stats.test.onesample.Group.median.perorganFDR2,file="../Datafiles/Statistics/2_OneSample_PerOrganFDR_NC.csv",row.names = FALSE)

###############Now we will conduct ANOVAS per compound, per organ, compared across time

zmelt2<-zmelt %>% group_by(Compound,organ) %>% filter(n()>5)

stats.test.anova.median<-zmelt2 %>% 
  group_by(Compound,organ) %>% 
  anova_test(value ~ time) %>% 
  adjust_pvalue(method = "fdr") %>%
  add_significance("p.adj")

write.csv(stats.test.anova.median,file="../Datafiles/Statistics/3_ANOVA_TimeComparison_NC.csv",row.names = FALSE)


###doing  FDR correction by organ
organtypes<-unique(zmelt2$organ)
stats.test.anova.median.perorganFDR<-data.frame()
for (i in 1:length(organtypes)){
  zmelt2organ<-zmelt2 %>% filter(organ==organtypes[i])
  stats.test.anova.median.each<-zmelt2organ %>% 
    group_by(Compound) %>% 
    anova_test(value ~ time) %>% 
    adjust_pvalue(method = "fdr") %>%
    add_significance("p.adj")
  stats.test.anova.median.perorganFDR<-rbind(stats.test.anova.median.perorganFDR,stats.test.anova.median.each)
}

write.csv(stats.test.anova.median.perorganFDR,file="../Datafiles/Statistics/4_ANOVA_TimeComparison_FDRPerOrgan_NC.csv",row.names = FALSE)

###############Now we will conduct ANOVAS per compound compared across organ
#first create needed columns for ANOVA

stats.test.anova.median<-zmelt %>% 
  group_by(Compound) %>% 
  anova_test(value ~ organ) %>% 
  adjust_pvalue(method = "fdr") %>%
  add_significance("p.adj")

write.csv(stats.test.anova.median,file="../Datafiles/Statistics/3_2_ANOVA_OrganComparison_NC.csv",row.names = FALSE)

stats.test.anova.median<-zmelt %>% 
  group_by(Compound) %>% 
  anova_test(value ~ organ + time) %>% 
  adjust_pvalue(method = "fdr") %>%
  add_significance("p.adj")

write.csv(stats.test.anova.median,file="../Datafiles/Statistics/3_2_ANOVA_OrganPlusTimeComparison_NC.csv",row.names = FALSE)


zmelt3<-zmelt %>% group_by(Compound,organ) %>% filter(!organ == "Aorta") %>%
  mutate(organ_time = paste(organ,time,sep = "_"))

stats.test.anova.median2<-zmelt3 %>% 
  group_by(Compound) %>% 
  anova_test(value ~ organ_time) %>% 
  adjust_pvalue(method = "fdr") %>%
  add_significance("p.adj")

write.csv(stats.test.anova.median2,file="../Datafiles/Statistics/4_2_ANOVA_OrganTime_NC.csv",row.names = FALSE)

#Here I will do t-test 0 min vs Fed minutes
#first create a new column called fed vs not
zmelt2 <- zmelt2 %>%
  mutate(FedvsFast = "Fed") %>%
  mutate(FedvsFast = ifelse(time == 0,"Fast",FedvsFast))

stats.test.t.fastvsfed.median<-zmelt2 %>% 
  group_by(Compound,organ) %>% 
  t_test(value ~ FedvsFast) %>% 
  adjust_pvalue(method = "fdr") %>%
  add_significance("p.adj")

write.csv(stats.test.t.fastvsfed.median,file="../Datafiles/Statistics/5_TTest_FastvsFed_NC.csv",row.names = FALSE)

###doing tiny project of only FDR correction by organ
organtypes<-unique(zmelt2$organ)
stats.test.t.fastvsfed.median.perorganFDR<-data.frame()
for (i in 1:length(organtypes)){
  zmelt2organ<-zmelt2 %>% filter(organ==organtypes[i])
  stats.test.t.fastvsfed.median.each<-zmelt2organ %>% 
    group_by(Compound) %>% 
    t_test(value ~ FedvsFast) %>% 
    adjust_pvalue(method = "fdr") %>%
    add_significance("p.adj")
  stats.test.t.fastvsfed.median.perorganFDR<-rbind(stats.test.t.fastvsfed.median.perorganFDR,stats.test.t.fastvsfed.median.each)
}

write.csv(stats.test.t.fastvsfed.median.perorganFDR,file="../Datafiles/Statistics/6_TTest_FastvsFed_PerOrganFDR_NC.csv",row.names = FALSE)


#Here I will do t-test 0 min vs Fed minutes
#first create a new column called fed vs not
zmelt2 <- zmelt2 %>%
  mutate(FedvsFast = "Fed") %>%
  mutate(FedvsFast = ifelse(time == 0,"Fast",FedvsFast))

stats.test.t.fastvsfed.median<-zmelt2 %>% 
  group_by(Compound,organ) %>% 
  t_test(value ~ time,
         comparisons = list(c(0,30),c(0,60),c(0,120))) %>% 
  adjust_pvalue(method = "fdr") %>%
  add_significance("p.adj")

write.csv(stats.test.t.fastvsfed.median,file="../Datafiles/Statistics/5_2_TTest_0minvs3060120min_NC.csv",row.names = FALSE)

###doing FDR correction by organ
organtypes<-unique(zmelt2$organ)
stats.test.t.fastvsfed.median.perorganFDR<-data.frame()
for (i in 1:length(organtypes)){
  zmelt2organ<-zmelt2 %>% filter(organ==organtypes[i])
  stats.test.t.fastvsfed.median.each<-zmelt2organ %>% 
    group_by(Compound,organ) %>% 
    t_test(value ~ time,
           comparisons = list(c(0,30),c(0,60),c(0,120)
           ),detailed = TRUE) %>% 
    adjust_pvalue(method = "fdr") %>%
    add_significance("p.adj")
  stats.test.t.fastvsfed.median.perorganFDR<-rbind(stats.test.t.fastvsfed.median.perorganFDR,stats.test.t.fastvsfed.median.each)
}

zmelt2median <- zmelt2 %>%
  group_by(Compound,time,organ) %>%
  summarise(medianvalues = median(value),
            meanvalues = mean(value)) 
  

zmelt2median <- zmelt2median %>%
  group_by(Compound,organ) %>%
  mutate(medianvaluesratioabsolutevalue = abs(medianvalues-medianvalues[time==0])) %>%
  mutate(meanvaluesratioabsolutevalue = abs(meanvalues-meanvalues[time==0])) %>%
  filter(!time == 0) %>%
  rename(group2 = time) %>%
  select(Compound,organ,group2,medianvaluesratioabsolutevalue,meanvaluesratioabsolutevalue,medianvalues) %>%
  mutate(releaseuptake = ifelse(medianvalues > 0, "release",
                                ifelse(medianvalues < 0, "uptake","none")))

stats.test.t.fastvsfed.median.perorganFDR_absfoldchange<-merge(stats.test.t.fastvsfed.median.perorganFDR,zmelt2median)

write.csv(stats.test.t.fastvsfed.median.perorganFDR,file="../Datafiles/Statistics/6_2_TTest_0minvs3060120min_PerOrganFDR_NC.csv",row.names = FALSE)
write.csv(stats.test.t.fastvsfed.median.perorganFDR_absfoldchange,file="../Datafiles/Statistics/6_2_TTest_0minvs3060120min_FoldChangeAbsolute_PerOrganFDR_NC.csv",row.names = FALSE)

#One Sample T Test Grouped By Time
stats.test.onesample.GroupTime.median<-zmelt %>% 
  group_by(Compound,organ,time) %>% 
  t_test(value~1,mu = 0,alternative = "two.sided") %>% 
  adjust_pvalue(method = "fdr") %>%
  add_significance("p.adj")

write.csv(stats.test.onesample.GroupTime.median,file="../Datafiles/Statistics/7_OneSample_GroupTime_NC.csv",row.names = FALSE)

###doing FDR correction by organ
organtypes<-unique(zmelt$organ)
stats.test.onesample.GroupTime.median.perorganFDR<-data.frame()
for (i in 1:length(organtypes)){
  zmeltorgan<-zmelt %>% filter(organ==organtypes[i])
  stats.test.onesample.GroupTime.median.each<-zmeltorgan %>% 
    group_by(Compound,organ,time) %>% 
    t_test(value~1,mu = 0,alternative = "two.sided") %>% 
    adjust_pvalue(method = "fdr") %>%
    add_significance("p.adj")
  stats.test.onesample.GroupTime.median.perorganFDR<-rbind(stats.test.onesample.GroupTime.median.perorganFDR,stats.test.onesample.GroupTime.median.each)
}

write.csv(stats.test.onesample.GroupTime.median.perorganFDR,file="../Datafiles/Statistics/8_OneSample_GroupTime_PerOrganFDR_NC2.csv",row.names = FALSE)

