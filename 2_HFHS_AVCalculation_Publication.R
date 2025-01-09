############let's bring everything in together, rbind them all, and then match with values...
library(stringi)
library(dplyr)
library(ggplot2)
library(stringr)
library(purrr)
library(ggpubr)
library(reshape2)
library(ggthemes)
library(readxl)

mycurrentdirectory = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(mycurrentdirectory)

excel_sheets("../Datafiles/HF_Official_Data.xlsx")
xpublish <- read_xlsx("../Datafiles/HF_Official_Data.xlsx",sheet = "Ion Counts")
xpublish<-xpublish %>%
  select(!colnames(xpublish)[grepl("_U_",colnames(xpublish))])

z8 <-melt(xpublish)

organ<-as.character(z8$variable)
unique(organ)
minute<-organ
replicate<-organ
pig<-organ

minute<-stri_replace_all_regex(organ,"(.*)_(.*)_(.*)_(.*)","$2") #Making new column for minute
unique(minute)

organ<-stri_replace_all_regex(organ,"(.*)_(.*)_(.*)_(.*)","$3") #Making new column for organ
unique(organ)

replicate<-stri_replace_all_regex(replicate,"(.*)_(.*)_(.*)_(.*)","$4") #Making new column for minute
unique(replicate)

pig<-stri_replace_all_regex(pig,"(.*)_(.*)_(.*)_(.*)","$1")
pig<-stri_replace_all_regex(pig,"HF",'') 
unique(pig)
z9<-cbind.data.frame(pig,organ,minute,replicate,z8)

z9 <- z9 %>%
  filter(!replicate %in% c(1.2,2.2,3.2))

z9 <- z9 %>%select(!contains(c("qc")))
z9 <- z9[!grepl("qc", z9$pig),]
colnames(z9)
unique(z9$pig)
detach("package:plyr", unload = TRUE)

#Getting AV ratio
z101<-z9 %>% group_by(pig,minute,Name)%>%mutate(AVratio=value/value[organ=="FA" & replicate == "1"])#Normalize by artery, replace "FA" to your artery name
z102<-z9 %>% group_by(pig,minute,Name)%>%mutate(AVratio=value/value[organ=="FA" & replicate == "2"])#Normalize by artery, replace "FA" to your artery name
z103<-z9 %>% group_by(pig,minute,Name)%>%mutate(AVratio=value/value[organ=="FA" & replicate == "3"])#Normalize by artery, replace "FA" to your artery name

z101$combine <- "FA1"
z102$combine <- "FA2"
z103$combine <- "FA3"

zcombine<-rbind(z101,z102,z103)
z10<-zcombine
z11<-z10%>%group_by(pig,minute,organ,Name)%>%summarise(medianvalue=median(AVratio))


#Pig 1 Lung AV
unique(z11$organ)
z12<-z11[!grepl("^H$|LV|RV|A|FA2|FA",z11$organ),]

lung<-z9[grepl("^A$|^RV$",z9$organ) & z9$pig == 1,]
unique(lung$organ)
lung_AV1<-lung %>% group_by(pig,
                            minute,Name)%>%mutate(AVratio=value/value[organ=="RV" & replicate == "1"]) #Normalize for lung
lung_AV2<-lung %>% group_by(pig,
                            minute,Name)%>%mutate(AVratio=value/value[organ=="RV" & replicate == "2"]) #Normalize for lung
lung_AV3<-lung %>% group_by(pig,
                            minute,Name)%>%mutate(AVratio=value/value[organ=="RV" & replicate == "3"]) #Normalize for lung

lung_AV1$combine<-"RV1"
lung_AV2$combine<-"RV2"
lung_AV3$combine<-"RV3"

lungcombine<-rbind(lung_AV1,lung_AV2,lung_AV3)

lungcombine<-lungcombine[!grepl("^RV",lungcombine$organ),]
lungcombine$organ[lungcombine$organ=="A"]<-"lung"
lungcombine_lung1_120min<-lungcombine%>%group_by(pig,minute,organ,Name)%>%summarise(medianvalue=median(AVratio))

#Pig 2-5 Lung AV
lung<-z9[grepl("^LV$|^RV$",z9$organ) & (!z9$pig==1),]
unique(lung$organ)
lung_AV1<-lung %>% group_by(pig,minute,Name)%>%mutate(AVratio=value/value[organ=="RV" & replicate == "1"]) #Normalize for lung
lung_AV2<-lung %>% group_by(pig,minute,Name)%>%mutate(AVratio=value/value[organ=="RV" & replicate == "2"]) #Normalize for lung
lung_AV3<-lung %>% group_by(pig,minute,Name)%>%mutate(AVratio=value/value[organ=="RV" & replicate == "3"]) #Normalize for lung

lung_AV1$combine<-"RV1"
lung_AV2$combine<-"RV2"
lung_AV3$combine<-"RV3"

lungcombine<-rbind(lung_AV1,lung_AV2,lung_AV3)

lungcombine<-lungcombine[!grepl("^RV",lungcombine$organ),]
lungcombine$organ[lungcombine$organ=="LV"]<-"lung"
lungcombine<-lungcombine%>%group_by(pig,minute,organ,Name)%>%summarise(medianvalue=median(AVratio))

lungcombine<-rbind(lungcombine,lungcombine_lung1_120min)

####Liver == 9 replicates ####
liver<-z9[grepl("^H$|^P$|^FA$",z9$organ),]
unique(liver$organ)

liver_AV1<-liver %>% group_by(pig,minute,Name)%>%mutate(AVratio=value/(value[organ=="P" & replicate == "1"]*0.8+value[organ=="FA" & replicate == "1"]*0.2)) #Normalize for liver
liver_AV2<-liver %>% group_by(pig,minute,Name)%>%mutate(AVratio=value/(value[organ=="P" & replicate == "2"]*0.8+value[organ=="FA" & replicate == "2"]*0.2)) #Normalize for liver
liver_AV3<-liver %>% group_by(pig,minute,Name)%>%mutate(AVratio=value/(value[organ=="P" & replicate == "3"]*0.8+value[organ=="FA" & replicate == "3"]*0.2)) #Normalize for liver

liver_AV1$combine<-"PFA1"
liver_AV2$combine<-"PFA2"
liver_AV3$combine<-"PFA3"

livercombine<-rbind(liver_AV1,liver_AV2,liver_AV3)

livercombine<-livercombine[!grepl("^P$|^FA$",livercombine$organ),]
livercombine$organ[livercombine$organ=="H"]<-"liver"
livercombine<-livercombine%>%group_by(pig,minute,organ,Name)%>%summarise(medianvalue=median(AVratio))

z14 <- rbind(z12, lungcombine,livercombine)
z14$AVratio<-log2(z14$medianvalue)

z14<-z14 %>% 
  group_by(organ,minute,Name) %>% 
  filter(length(Name)>=1) 
z14<-z14[!grepl("^FA$",z14$organ),]

#########################
z14$group<-paste(z14$organ,z14$minute,z14$pig,sep="_")
z15 <- dcast(z14, Name~group, value.var="AVratio")
z16<-z15
colnames(z16)<-tolower(colnames(z16))

publisheable<-list(xpublish,z16)
names(publisheable)<-c("Ion Counts","log2VA")
library(writexl)
write_xlsx(publisheable,path = "../Datafiles/HF_Official_Data.xlsx")

