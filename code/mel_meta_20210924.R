
setwd("C:/learning/R/melanoma/data/data_therapy")
require(readxl)
library(meta)
library(metafor)
library(ggplot2)
library(dplyr)
library(grid)
library(gridExtra)
library(ggrepel)
library(forcats)



#####META 
# function
getlog<-function(upci,lowerci){
  selogHR= (log(upci)-log(lowerci))/(1.96*2)
  return(selogHR)
}

META_random<-function(feature,race){
Data<-read_xlsx(paste0('stage_ulcer_',race,'_meta.xlsx'),feature)
#####get lohHR and selogHR
if (feature==1){
Data$logHR = log(Data$AJCC_HR)
Data$selogHR = getlog(Data$HCI,Data$LCI)
}
if (feature==2){
  Data$logHR = log(Data$Ulceration_HR)
  Data$selogHR = getlog(Data$HCI,Data$LCI)
}
if(race=='asian'&feature==2){Data=Data[-3,]}
if(race=='asian'&feature==1){Data=Data[-4,]}
mg<-metagen(logHR,selogHR,
            studlab = Research,
            data=Data,
            comb.fixed = FALSE,
            comb.random = TRUE,
            hakn = FALSE,
            prediction = TRUE,
            sm='HR')
mgr=summary(mg)
META_logHR=mgr$random[["TE"]]
META_logHRse=mgr$random[["seTE"]]
return(c(META_logHR,META_logHRse))
}


meta_result=META_random(2,"asian")
meta_result_stage=META_random(1,"asian")

##########################################
##1.ACSIAN(CHINA,JAPAN,KOREA) Ulceration##
##########################################
##1.1 data preprocess
Data<-read_xlsx('C:/learning/R/melanoma/table/meta/asian_meta.xlsx',2)
##get lohHR and selogHR
Data$logHR = log(Data$Ulceration_HR)
Data$selogHR = getlog(Data$HCI,Data$LCI)

#1.2 meta analysis
#####DerSimonian-Laird estimator for tau^2
outcome = NULL
for( i in 1:3){
  data = Data[Data$Race_new==i,]
  if (i ==2){#remove extreme value 
    data<-data[-3,]
     }
  meta<-metagen(logHR,selogHR,
          studlab = Research,
          data=data ,
          comb.fixed = FALSE,
          comb.random = TRUE,
          hakn = FALSE,
          prediction = TRUE,
          sm='HR')
  print(meta)
  forest(meta)
  effet_size = c(meta$TE.random,meta$seTE.random)
  outcome<-rbind(outcome,effet_size)
  #grid.text("Asian Ulceration Forest Plot", .5, 0.97, gp=gpar(cex=1.5))
}

#1.3 ulceration seer and meta compare
  tb<-as.data.frame(outcome,row.names =c('China','Japan','Korea') )
  colnames(tb)<-c('logHR','logHRse')
  tb$HR<-exp(tb$logHR)
  tb$lowerci <- exp(tb$logHR-1.96*tb$logHRse)
  tb$upperci <- exp(tb$logHR+1.96*tb$logHRse)
# 1.2.1.2influence analysis
#inf.analysis <- InfluenceAnalysis(x = china_ulcer_meta,random = TRUE)
#inf.analysis

seer_asia=SEERlogHR(dat,"Non-Hispanic Asian or Pacific Islander",'Ulceration')
row = c(seer_asia[1],seer_asia[2],exp(seer_asia[1]),exp(seer_asia[1]-1.96*seer_asia[2]),exp(seer_asia[1]+1.96*seer_asia[2]))
tb<-rbind(tb,row)
row.names(tb)[4]<-'SEER_Ulcer'

# subgroup analysis 
#Sys.setlocale(category = "LC_ALL", locale = "us")
#submeta<-subgroup.analysis.mixed.effects(x = mg,subgroups = Data$Race)
#forest(submeta$m.random)

####################################
##2.ASCIAN(CHINA,JAPAN,KOREA) Stage#
####################################
##1.1 data preprocess
Data<-read_xlsx('C:/learning/R/melanoma/table/meta/asian_meta.xlsx',1)
##get lohHR and selogHR
Data$logHR = log(Data$AJCC_HR)
Data$selogHR = getlog(Data$HCI,Data$LCI)
#2.2 meta analysis
#####DerSimonian-Laird estimator for tau^2
outcome = NULL
l = c('China','Japan','Korea')
for( i in 1:3){
  data = Data[Data$Race==l[i],]
  if (i ==1){#remove extreme value 
    data<-data[-3,]
  }
  meta<-metagen(logHR,selogHR,
                studlab = Research,
                data=data ,
                comb.fixed = FALSE,
                comb.random = TRUE,
                hakn = FALSE,
                prediction = TRUE,
                sm='HR')
  print(meta)
  forest(meta)
  effet_size = c(meta$TE.random,meta$seTE.random)
  outcome<-rbind(outcome,effet_size)
  #grid.text("Asian Ulceration Forest Plot", .5, 0.97, gp=gpar(cex=1.5))
}
#inf.analysis <- InfluenceAnalysis(x = meta,random = TRUE)
#inf.analysis

#2.3 Stage seer and meta compare
tb<-as.data.frame(outcome,row.names =c('China','Japan','Korea') )
colnames(tb)<-c('logHR','logHRse')
tb$HR<-exp(tb$logHR)
tb$lowerci <- exp(tb$logHR-1.96*tb$logHRse)
tb$upperci <- exp(tb$logHR+1.96*tb$logHRse)

seer_asia_stage <-SEERlogHR_stage(dat,"Non-Hispanic Asian or Pacific Islander")
row <- c(seer_asia_stage[1],seer_asia_stage[2],exp(seer_asia_stage[1]),
        exp(seer_asia_stage[1]-1.96*seer_asia_stage[2]),
        exp(seer_asia_stage[1]+1.96*seer_asia_stage[2]))
tb<-rbind(tb,row)
row.names(tb)[4]<-'SEER_stage'

#p_val=SEER_META(seer_result_stage,meta_result)

################################
#3.ASCIAN(CHINA,JAPAN,KOREA) LM#
################################
# 3.1 data preprocess
Data<-read_xlsx('C:/learning/R/melanoma/table/meta/asian_meta.xlsx',4)
#get lohHR and selogHR
Data$logHR = log(Data$LNM)
Data$selogHR = getlog(Data$HCI,Data$LCI)

#3.2 meta analysis
######DerSimonian-Laird estimator for tau^2
outcome = NULL
l = c('China','Japan','Korea')
for( i in 1:3){
  data = Data[Data$Race==l[i],]
  
  meta<-metagen(logHR,selogHR,
                studlab = Author,
                data=data ,
                comb.fixed = FALSE,
                comb.random = TRUE,
                hakn = FALSE,
                prediction = TRUE,
                sm='HR')
  print(meta)
  forest(meta)
  effet_size = c(meta$TE.random,meta$seTE.random)
  outcome<-rbind(outcome,effet_size)
  #grid.text("Asian Ulceration Forest Plot", .5, 0.97, gp=gpar(cex=1.5))
}
#2.3 Stage seer and meta compare
tb<-as.data.frame(outcome,row.names =c('China','Japan','Korea') )
colnames(tb)<-c('logHR','logHRse')
tb$HR<-exp(tb$logHR)
tb$lowerci <- exp(tb$logHR-1.96*tb$logHRse)
tb$upperci <- exp(tb$logHR+1.96*tb$logHRse)

seer_asia_LM <- SEERlogHR(dat,"Non-Hispanic Asian or Pacific Islander",'LN_Metastasis')
row <- c(seer_asia_LM[1],seer_asia_LM[2],exp(seer_asia_LM[1]),
         exp(seer_asia_LM[1]-1.96*seer_asia_LM[2]),
         exp(seer_asia_LM[1]+1.96*seer_asia_LM[2]))
tb<-rbind(tb,row)
row.names(tb)[4]<-'SEER_LM'

######absract meta logHR and selogHR
#mgr=summary(mg)
#META_logHR=mgr$random[["TE"]]
#META_logHRse=mgr$random[["seTE"]]
#meta_result = c(META_logHR,META_logHRse)###1.1362298 0.1432725

###0.3670217 0.3467281 
 
######meta compare seer
#p_val=SEER_META(seer_result,meta_result)###0.02016613
#p_val=SEER_META(seer_result1,meta_result)##<e-8

#######################################
#4.ASCIAN(CHINA,JAPAN,KOREA) THICKNESS#
#######################################
Data<-read_xlsx('C:/learning/R/melanoma/table/meta/asian_meta.xlsx',3)
#####get lohHR and selogHR
# 4.1 data preprocess
Data$logHR = log(Data$Thickness)
Data$selogHR = getlog(Data$HCI,Data$LCI)
# 4.3 meta analysis
######DerSimonian-Laird estimator for tau^2
outcome = NULL
l = c('China','Japan','Korea')
for( i in 1:3){
  data = Data[Data$Race==l[i],]
  
  meta<-metagen(logHR,selogHR,
                studlab = Author,
                data=data ,
                comb.fixed = FALSE,
                comb.random = TRUE,
                hakn = FALSE,
                prediction = TRUE,
                sm='HR')
  print(meta)
  forest(meta)
  effet_size = c(meta$TE.random,meta$seTE.random)
  outcome<-rbind(outcome,effet_size)
  #grid.text("Asian Melanoma Thicknees Forest Plot", .5, .75, gp=gpar(cex=1.5))
}
# 4.3 Stage seer and meta compare
tb<-as.data.frame(outcome,row.names =c('China','Japan','Korea') )
colnames(tb)<-c('logHR','logHRse')
tb$HR<-exp(tb$logHR)
tb$lowerci <- exp(tb$logHR-1.96*tb$logHRse)
tb$upperci <- exp(tb$logHR+1.96*tb$logHRse)

seer_asia_Thick <- SEERlogHR(dat,"Non-Hispanic Asian or Pacific Islander",'Thickness')
row <- c(seer_asia_Thick[1],seer_asia_Thick[2],exp(seer_asia_Thick[1]),
         exp(seer_asia_Thick[1]-1.96*seer_asia_Thick[2]),
         exp(seer_asia_Thick[1]+1.96*seer_asia_Thick[2]))
tb<-rbind(tb,row)
row.names(tb)[4]<-'SEER_Thick'

#p_val=SEER_META(seer_result,meta_result)##0.1719233


#################################
#5.ASCIAN(CHINA,JAPAN,KOREA) AGE#
#################################
Data<-read_xlsx('C:/learning/R/melanoma/table/meta/asian_meta.xlsx',5)
#####get lohHR and selogHR
# 5.1 data preprocess
Data$logHR = log(Data$Age)
Data$selogHR = getlog(Data$HCI,Data$LCI)
# 5.2 meta analysis
######DerSimonian-Laird estimator for tau^2
outcome = NULL
l = c('China','Japan','Korea')
for( i in 1:3){
  data = Data[Data$Race==l[i],]
  
  meta<-metagen(logHR,selogHR,
                studlab = Author,
                data=data ,
                comb.fixed = FALSE,
                comb.random = TRUE,
                hakn = FALSE,
                prediction = TRUE,
                sm='HR')
  print(meta)
  forest(meta)
  effet_size = c(meta$TE.random,meta$seTE.random)
  outcome<-rbind(outcome,effet_size)
  #grid.text("Asian Melanoma Thicknees Forest Plot", .5, .75, gp=gpar(cex=1.5))
}

# 5.3 Stage seer and meta compare
tb<-as.data.frame(outcome,row.names =c('China','Japan','Korea') )
colnames(tb)<-c('logHR','logHRse')
tb$HR<-exp(tb$logHR)
tb$lowerci <- exp(tb$logHR-1.96*tb$logHRse)
tb$upperci <- exp(tb$logHR+1.96*tb$logHRse)

seer_asia_Age <- SEERlogHR(dat,"Non-Hispanic Asian or Pacific Islander",'Age')
row <- c(seer_asia_Age[1],seer_asia_Age[2],exp(seer_asia_Age[1]),
         exp(seer_asia_Age[1]-1.96*seer_asia_Age[2]),
         exp(seer_asia_Age[1]+1.96*seer_asia_Age[2]))
tb<-rbind(tb,row)
row.names(tb)[4]<-'SEER_Age'


#seer_result = SEERlogHR(dat,"Non-Hispanic White")##0.0568906898 0.0004556528
######meta compare seer
#p_val=SEER_META(seer_result,meta_result)##0.161
#p_val=SEER_META(seer_result,meta_result)##0.99










##########################################
##2.EUROPE(CHINA,JAPAN,KOREA) Thickness##
##########################################




seer_result = SEERlogHR(dat,"Non-Hispanic White")####0.0568906898 0.0004556528 
######meta compare seer
p_val=SEER_META(seer_result,meta_result)####1

########################################ulceration


seer_result = SEERlogHR(dat,"Non-Hispanic White")####0.42431110 0.01724422 
######meta compare seer
p_val=SEER_META(seer_result,meta_result)####0.1810782



seer_result = a$coefficients[which(rownames(a$coefficients)=="LNN1"),][c(1,3)]####0.19448486 0.03786802
######meta compare seer
p_val=SEER_META(seer_result,meta_result)####0.0006


seer_result = a$coefficients[which(rownames(a$coefficients)=="StageII"),][c(1,3)]##0.45151609 0.02168693
######meta compare seer
p_val=SEER_META(seer_result,meta_result)####0.0210433









meta_diff<-function(meta1,meta2){
  stat = (meta1[1]-meta2[1])/(sqrt(meta1[2]^2+meta2[2]^2))
  p_val = pnorm(stat)
  if (p_val<0.5){p_val=p_val}
  if (p_val>0.5){p_val=1-p_val}
  return(p_val)
}

####meta e vs a :ulcer
meta1 = c(0.656,0.0955)
meta2 = c(0.503,0.085)
meta_diff(meta1,meta2)

####meta e vs a :age
meta1 = c(0.023, 0.00685)
meta2 = c(0.028,0.0034)
meta_diff(meta1,meta2)
####meta e vs a :thickness
meta1 = c(0.1112,0.0236)
meta2 = c(0.159, 0.0474 )
meta_diff(meta1,meta2)


#####################################################################################
#########seer-meta compare forest plot
###############################################################################
library(forestplot)
library(readxl)

data1<-read_xlsx("C:/learning/R/melanoma/manusript/table/meta/Asian.xlsx")
colnames(data1)[1] =''
data1<- data1[,c(1,4,7)]

data1<-as.matrix(data1)
data1 <- rbind(c(NA,"HR","CI"),data1) 

styles <- fpShapesGp(
  lines = list(
    gpar(col = ""),
    gpar(col = "red"),
    gpar(col = "blue"),
    gpar(col = "yellow"),
    gpar(col = "green"),
    gpar(col = "red"),
    gpar(col = "blue"),
    gpar(col = "yellow"),
    gpar(col = "green"),
    gpar(col = "red"),
    gpar(col = "blue"),
    gpar(col = "yellow"),
    gpar(col = "green"),
    gpar(col = "red"),
    gpar(col = "blue"),
    gpar(col = "yellow"),
    gpar(col = "green"),
    gpar(col = "red"),
    gpar(col = "blue"),
    gpar(col = "yellow"),
    gpar(col = "green")
  ),
  box = list(
    gpar(fill = ""),
    gpar(fill = "red"),
    gpar(fill = "blue"),
    gpar(fill = "yellow"),
    gpar(fill = "green"),
    gpar(fill = "red"),
    gpar(fill = "blue"),
    gpar(fill = "yellow"),
    gpar(fill = "green"),
    gpar(fill = "red"),
    gpar(fill = "blue"),
    gpar(fill = "yellow"),
    gpar(fill = "green"),
    gpar(fill = "red"),
    gpar(fill = "blue"),
    gpar(fill = "yellow"),
    gpar(fill = "green"),
    gpar(fill = "red"),
    gpar(fill = "blue"),
    gpar(fill = "yellow"),
    gpar(fill = "green")
  ) )

forestplot(data1,  
           c(NA,as.numeric(data1[,2][2:21])),
           c(NA,as.numeric(substr(data1[,3][2:21],2,5))), 
           c(NA,as.numeric(substr(data1[,3][2:21],7,10))),
           zero = 1, 
           xlog=FALSE, 
           fn.ci_norm = fpDrawCircleCI,
           boxsize = 0.15, 
           shapes_gp = styles,
           lty.ci = 7,   
           lwd.ci = 3,   
           xticks=seq(1,8,0.2),
           ci.vertices.height = 0.15, 
           txt_gp = fpTxtGp(ticks = gpar(cex = 0.5), xlab = gpar(cex = 0.9), cex = 0.8), 
           lineheight = "auto", 
           xlab="Differences about feature Hazard Ratio between SEER and META" 
)










load("C:/learning/R/melanoma/manusript/data/final_data_survival_20210901.RData")
###############################################################################
###################extract log(HR) and se(log(HR)) of SEER from specific Race
################################################################################
library(survival)
library(survminer)
#########################stage seer HR extract
SEERlogHR_stage<-function(alldata,race){
  subdat<-alldata[alldata$Race==race,]
  subdat$stagecomb<-'Unknown'
  subdat$stagecomb[which(subdat$Stage%in%c('I','II'))]='I/II'
  subdat$stagecomb[which(subdat$Stage%in%c('III','IV'))]='III/IV'
  subdat$stagecomb<-factor(subdat$stagecomb,levels=c('I/II','III/IV','Unknown'))
  model_os<-coxph(Surv(Time,OSevent)~
                    Year+
                    Age+
                    sex+
                    Partner+
                    Site+
                    Histology+
                    Poverty+
                    Education+
                    stagecomb+
                    LN+
                    Ulceration+
                    Thickness+
                    PS_Surg+LN_Surg+Oth_Surg+Radiation+Chemotherapy
                  ,data=subdat)
  report=summary(model_os)
  coef_secoef=report$coefficients[which(rownames(report$coefficients)=="stagecombIII/IV"),][c(1,3)]
  return(coef_secoef)
}

#seer_result_stage=SEERlogHR_stage(dat,"Non-Hispanic Asian or Pacific Islander")#1.358
#seer_result_stage1=SEERlogHR_extr(dat,"Non-Hispanic White")

#subdat<-dat[dat$Race=="Non-Hispanic White",]
#model_os<-coxph(Surv(Time,OSevent)~Year+Age+sex+Partner+Site+Histology+Poverty+
                  #Education+Stage+Ulceration+Thickness
                #,data=subdat)
##############################################################################################
##########################other features seer HR extract
SEERlogHR<-function(alldata,race,para){
  subdat<-alldata[alldata$Race==race,]
  if (para == 'Ulceration')
  { model_os<-coxph(Surv(Time,OSevent)~Year+Age+sex+Partner+Site+
                      Histology+
                      Poverty+
                      Education+
                      Stage+
                      LN+
                      Ulceration+
                      Thickness+
                      PS_Surg+
                      LN_Surg+
                      Oth_Surg+
                      Radiation+
                      Chemotherapy,data=subdat)
    report=summary(model_os)
    coef_secoef=report$coefficients[which(rownames(report$coefficients)=='UlcerationYes'),][c(1,3)]}
  else if (para == 'LN_Metastasis') 
  {model_os<-coxph(Surv(Time,OSevent)~Year+Age+sex+Partner+Site+
                     Histology+
                     Poverty+
                     Education+
                     #Stage+
                     LN_Metastasis+
                     Ulceration+
                     Thickness+
                     PS_Surg+
                     LN_Surg+
                     Oth_Surg+
                     Radiation+
                     Chemotherapy,data=subdat)
    report=summary(model_os)
    coef_secoef=report$coefficients[which(rownames(report$coefficients)=='LN_MetastasisYes'),][c(1,3)]}
  else if (para == 'Thickness')
  { model_os<-coxph(Surv(Time,OSevent)~Year+Age+sex+Partner+Site+
                      Histology+
                      Poverty+
                      Education+
                      #Stage+
                      LN+
                      Ulceration+
                      Thickness+
                      PS_Surg+
                      LN_Surg+
                      Oth_Surg+
                      Radiation+
                      Chemotherapy,data=subdat)
     report=summary(model_os)
     coef_secoef=report$coefficients[which(rownames(report$coefficients)=='Thickness'),][c(1,3)]}
  else if (para == 'Age')
  {model_os<-coxph(Surv(Time,OSevent)~Year+Age+sex+Partner+Site+
                     Histology+
                     Poverty+
                     Education+
                     Stage+
                     LN+
                     Ulceration+
                     Thickness+
                     PS_Surg+
                     LN_Surg+
                     Oth_Surg+
                     Radiation+
                     Chemotherapy,data=subdat)
     report=summary(model_os)
     coef_secoef=report$coefficients[which(rownames(report$coefficients)=='Age'),][c(1,3)]}
  return(coef_secoef)
}

#seer_result=SEERlogHR(dat,"Non-Hispanic Asian or Pacific Islander")##0.3183770 0.1648376 
#seer_result1=SEERlogHR(dat,"Non-Hispanic White")##0.42361387 0.01724056

##########################################################################################
#####compare SEER and META results
##########################################################################################
SEER_META<-function(seer,meta){
  z_statistic = (seer[1]-meta[1])/sqrt(seer[2]^2+meta[2]^2)
  p_val = pnorm(z_statistic)
  if (p_val<0.5){p_val=p_val}
  if (p_val>0.5){p_val=1-p_val}
  return(p_val)
}
SEERlogHR_extr(dat,"Non-Hispanic White")
META_random("stage","asian")

p_value=SEER_META(seer_result,meta_result)###seer asian vs meta asian 0.03816915
p_value1=SEER_META(seer_result1,meta_result)###seer white vs meta asian 0.008328026 
####seer white vs meta europe white 0.433
p_value2=SEER_META(seer_result_stage,meta_result_stage)
p_value3=SEER_META(seer_result_stage1,meta_result_stage)






# test
subdat<-dat[dat$Race=="Non-Hispanic Asian or Pacific Islander",]
model_os<-coxph(Surv(Time,OSevent)~Year+Age+sex+Partner+Site+
                  Histology+
                  Poverty+
                  Education+
                  Stage+
                  LN_Metastasis+
                  Ulceration+
                  Thickness+
                  PS_Surg+
                  LN_Surg+
                  Oth_Surg+
                  Radiation+
                  Chemotherapy,data=subdat)
report=summary(model_os)