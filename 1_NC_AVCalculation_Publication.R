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

#Valine normalization step
list.files("../Datafiles")
excel_sheets("../Datafiles/NC_Official_Data.xlsx")
xpublish <- read_xlsx("../Datafiles/NC_Official_Data.xlsx",sheet = "Ion Counts")

x <- xpublish %>%
  mutate(metabolite = ifelse(duplicated(tolower(metabolite))==TRUE,paste(metabolite,"duplicate",sep = "_"),metabolite))

colnames(x)
z8 <-melt(x,id.vars = "metabolite")
z8 <- z8 %>%
  mutate(value = as.numeric(value)) %>%
  mutate(variable = as.character(variable)) %>%
  rename(Name = `metabolite`)

####Fixing up the names to remove repeats
variable<-as.character(z8$variable)
repeats<-variable[grepl("2021",variable)]
repeatsunique<-unique(repeats)
repeatsunique<-stri_replace_all_regex(repeatsunique,"_2021[:digit:]+","")
repeatsunique<-paste("^",repeatsunique,"$",sep = "")

####Remove from original dataframe z8 ######
cleaning<-paste(repeatsunique,collapse = "|")
z8<-z8[!grepl(cleaning,z8$variable),]

####fix up 2021
variable<-as.character(z8$variable)
variable<-stri_replace_all_regex(variable,"_2021[:digit:]+","") #Making new column for organ
unique(variable)
z8$variable<-variable

organ<-as.character(z8$variable)
unique(organ)
minute<-organ
replicate<-organ
pig<-organ

minute<-stri_replace_all_regex(organ,"(.*)_(.*)_(.*)","$2") #Making new column for minute
unique(minute)


organ<-stri_replace_all_regex(organ,"(.*)_(.*)_(.*)","$3") #Making new column for organ
unique(organ)

replicate<-stri_replace_all_regex(replicate,"(.*)_(.*)_(.*)","$1") #Making new column for minute
unique(replicate)

pig<-stri_replace_all_regex(pig,"(.*)_(.*)_(.*)","$1")
unique(pig)
z9<-cbind.data.frame(pig,organ,minute,replicate,z8)

zotherorgans<-z9

z9 <- z9 %>%select(!contains(c("qc")))
z9 <- z9[!grepl("qc", z9$pig),]
colnames(z9)
unique(z9$pig)
detach("package:plyr", unload = TRUE)
dput(unique(organ))

#Getting AV ratio
z101<-z9 %>% group_by(pig,minute,Name)%>%mutate(AVratio=value/value[organ=="FA"]) #Normalize by artery, replace "FA" to your artery name

#######################################################

#AV ratio for other organs
lung<-zotherorgans %>%
  filter(organ %in% c("LV","RV"))
unique(zotherorgans$organ)
lung_AV1<-lung %>% group_by(pig,minute,Name)%>%mutate(AVratio=value/value[organ=="RV"]) #Normalize for lung

lungcombine<-lung_AV1 

lungcombine<-lungcombine %>%
  filter(!organ == "RV")
lungcombine$organ[lungcombine$organ=="LV"]<-"Lung"

####Liver ####
liver<-zotherorgans %>% 
  filter(organ %in% c("FA","P","H"))
unique(zotherorgans$organ)
liver_AV1<-liver %>% group_by(pig,minute,Name)%>%mutate(AVratio=value/(value[organ=="P"]*0.8+value[organ=="FA"]*0.2)) #Normalize for liver

livercombine<-liver_AV1 

livercombine<-livercombine %>%
  filter(!organ %in% c("FA","P"))

livercombine$organ[livercombine$organ=="H"]<-"Liver"
unique(z101$organ)
z101 <- z101 %>%
  filter(!organ %in% c("LV","RV","H","FA"))

z10 <- rbind(z101, lungcombine, livercombine)
  
z11<-z10 %>%
  rename(medianvalue = AVratio)

z14 <- z11 
z14$AVratio<-log2(z14$medianvalue)

z14<-z14 %>% 
  group_by(organ,minute,Name) %>% 
  filter(length(Name)>=3) 
z14<-z14[!grepl("FA",z14$organ),]

organtypes<-c("Femoral artery", "Intestine", "Colon", "Liver", 
              "Spleen", "Head (Brain)", "Leg (Muscle)", "Heart", "Kidney", 
              "Ear (Skin)", "Aorta", "Left Ventricle", "Right Ventricle")

organtypetochange<-
  c("FA", "P", "C", "H", 
    "S", "J", "FV", "CS", "R", 
    "E", "A", "LV", "RV")

dforganchange<-data.frame(organtypes = organtypes,
                          organtypetochange = organtypetochange)

for (i in 1:nrow(dforganchange)) {
  z14$organ[z14$organ == dforganchange$organtypetochange[i]] <- dforganchange$organtypes[i]
}

z14$minute<-stri_replace_all_regex(z14$minute,"000m","0min") #Making new column for minute
z14$minute<-stri_replace_all_regex(z14$minute,"030m","30min") #Making new column for minute
z14$minute<-stri_replace_all_regex(z14$minute,"060m","60min") #Making new column for minute
z14$minute<-stri_replace_all_regex(z14$minute,"120m","120min") #Making new column for minute

z14$pig<-stri_replace_all_regex(z14$pig,"FP",'Ct') 

unique(z14$organ)
z14$group<-paste(z14$organ,z14$minute,z14$pig,sep="_")
z15 <- dcast(z14, Name~group, value.var="AVratio")

colnames(z15)

publisheable<-list(x,z15)
names(publisheable)<-c("Ion Counts","log2VA")
library(writexl)
write_xlsx(publisheable,path = "../Datafiles/NC_Official_Data.xlsx")


