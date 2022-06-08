library(meta)
library(metafor)
library(readxl)
library(ggplot2)
library(ggpubr)
library(dmetar)

#########################################
#European stage proportion meta analysis#
#########################################
Data_Europe<-read_xlsx('C:/learning/R/melanoma/manusript/table/stage_dist_meta.xlsx',1)
#Data_after_deletion<-Data[c(-1,-3),]
Europe_meta_stage <- metaprop(III_IV_number, 
                               Total_number, 
                               studlab=Author, 
                               sm="PFT", 
                               data=Data_Europe, 
                               method="Inverse", 
                               method.tau="DL")

summary(Europe_meta_stage)

inf.analysis <- InfluenceAnalysis(x = Europe_meta_stage,random = TRUE)
inf.analysis

Europe_meta_plot<-forest.meta(Europe_meta_stage, 
                              layout="RevMan5", 
                              xlab="Proportion",
                              comb.r=T, 
                              comb.f=F, 
                              xlim = c(0,0.65), 
                              fontsize=10, 
                              digits=3)


######################################
#Asian stage proportion meta analysis#
######################################
Data_Asian<-read_xlsx('C:/learning/R/melanoma/manusript/table/stage_dist_meta.xlsx',2)
#Data_after_deletion<-Data[c(-1,-3),]
Asian_meta_stage <- metaprop(III_IV_number, 
                              Total_number, 
                              studlab=Author, 
                              sm="PFT", 
                              data=Data_Asian, 
                              method="Inverse", 
                              method.tau="DL")

summary(Asian_meta_stage)

inf.analysis <- InfluenceAnalysis(x = Asian_meta_stage,random = TRUE)
inf.analysis

Asian_meta_plot<-forest.meta(Asian_meta_stage, 
                              layout="RevMan5", 
                              xlab="Proportion",
                              comb.r=T, 
                              comb.f=F, 
                              xlim = c(0,0.65), 
                              fontsize=10, 
                              digits=3)

####################################
#Egger's regression Publication bias  
####################################
Data_comb<-rbind(Data_Europe,Data_Asian)
Data_comb$Race<-'Europe'
index<-c(grep('China',Data_comb$Country),
         grep('Korea',Data_comb$Country),
         grep('Japan',Data_comb$Country),
         grep('Taiwan',Data_comb$Country))
Data_comb$Race[index]<-'China'
Data_comb$Race<-as.factor(Data_comb$Race)
All_meta_stage <- metaprop(III_IV_number, 
                              Total_number, 
                              studlab=Author, 
                              sm="PFT", 
                              data=Data_comb, 
                              method="Inverse", 
                              method.tau="DL")

metareg(All_meta_stage, ~ Race)

#############
#funnel plot
#############
funnel.meta(Asian_meta_stage)
funnel.meta(Europe_meta_stage)
metabias(All_meta_stage, method.bias = "peters")
# Sys.setlocale(category = "LC_ALL", locale = "us")
# submeta<-subgroup.analysis.mixed.effects(x = All_meta_stage,subgroups = Data_comb$Race)
# forest(submeta$m.random)
