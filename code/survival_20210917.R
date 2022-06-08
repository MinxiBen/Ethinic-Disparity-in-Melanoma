library(survival)
library(survminer)

# 1. path setting
path_split<-strsplit(getwd(),split="/")[[1]]
parent_path<-paste(path_split[1:(length(path_split)-1)],collapse="/")
data_dir<-paste(parent_path,"data",sep="/")

data_name<-"final_data_survival_20210901.RData"
data_name1<-"original_data_20210901.RData"
data_path<-paste(data_dir,data_name,sep="/")
data_path1<-paste(data_dir,data_name1,sep="/")
figure_path<-paste(parent_path,"figure",sep="/")

# 2. load data dat
load(data_path)
load(data_path1)

#####################
# 3.Overall Survival#
#####################
surv<-survfit(Surv(Time,OSevent)~Race,data=dat)
# 5-year overall survival for each race
summary(surv)
overall_survplot<-ggsurvplot(surv, data =dat ,
                  risk.table = "nrisk_cumcensor",
                  palette = "jco",
                  tables.height = 0.2, 
                  tables.y.text = FALSE, pval = TRUE,
                  # pval.size = 8,
                  break.time.by=30,
                  xlim=c(0,150),
                  title ="Overall Survival of SEER Melanoma patients from 2004-2014(N=149115)",
                  legend.title="",
                  legend.cex=2,
                  ylab="Overall Survival Probability",
                  xlab="Months",
                  font.tickslab = 15,
                  font.x = 17,
                  font.y = 17,
                  font.legend = 17,
                  font.title = 20,
                  tables.theme = theme_cleantable(),
                  legend.labs=c('White','American African','Asian','American Indian','Hispanics'),
                  ggtheme = theme_classic() + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
)

###############################
# 4.Disease-specific Survival #
###############################
surv_dss<-survfit(Surv(Time,DSevent)~Race,data=dat)
## 5-year disease-specific survival for each race
summary(surv_dss)
disease_specific_survplot<-ggsurvplot(surv_dss, data =dat ,
           risk.table = "nrisk_cumcensor",
           palette = "jco",
           tables.height = 0.2, 
           tables.y.text = FALSE, pval = TRUE,
           break.time.by=30,
           xlim=c(0,150),
           title="Disease-specific Survival of SEER Melanoma patients from 2004-2014(N=149115)",
           legend.title ="",
           legend.cex=2,
           font.tickslab = 15,
           font.x = 17,
           font.y = 17,
           font.legend = 17,
           font.title = 20,
           ylab="Disease-specific Survival Probability",
           xlab="Months",
           tables.theme = theme_cleantable(),
           legend.labs=c('White','American African','Asian','American Indian','Hispanics'),
           ggtheme = theme_classic() + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
)

# arrange the survival plot and print the output
splots<-list()
splots[[1]]<-overall_survplot
splots[[2]]<-disease_specific_survplot
#tiff("Surv_manual.tiff",width=1870,height=900)
setEPS()
postscript("Surv_manual.eps",width=800,height=400)
arrange_ggsurvplots(splots, print = TRUE,
                    ncol = 2, nrow = 1)
dev.off()

#################
# 5. PHCox model#
#################
# check violation of PH 
km<-survfit(Surv(Time,DSevent)~Race,data=dat)
plot(km,fun='cloglog',xlab='time in months using  logarithmic scale',ylab='log-log survival',
     main='log-log curves by Race')

# need remove ulceration unknown
km1<-survfit(Surv(Time,DSevent)~Ulceration,data=dat)
plot(km1,fun='cloglog',xlab='time in months using  logarithmic scale',ylab='log-log survival',
     main='log-log curves by Ulceration')

km2<-survfit(Surv(Time,DSevent)~Stage,data=dat)
plot(km2,fun='cloglog',xlab='time in months using  logarithmic scale',ylab='log-log survival',
     main='log-log curves by Stage')

km3<-survfit(Surv(Time,DSevent)~PS_Surg,data=dat)
plot(km3,fun='cloglog',xlab='time in months using  logarithmic scale',ylab='log-log survival',
     main='log-log curves by PS_Surge')

# Other_Surg
km4<-survfit(Surv(Time,DSevent)~Oth_Surg,data=dat)
plot(km4,fun='cloglog',xlab='time in months using  logarithmic scale',ylab='log-log survival',
     main='log-log curves by Oth_Surge')

# Histoly remove OTHER
km5<-survfit(Surv(Time,DSevent)~Histology,data=dat)
plot(km5,fun='cloglog',xlab='time in months using  logarithmic scale',ylab='log-log survival',
     main='log-log curves by LN_Surge')

km6<-survfit(Surv(Time,DSevent)~Site,data=dat)
plot(km6,fun='cloglog',xlab='time in months using  logarithmic scale',ylab='log-log survival',
     main='log-log curves by LN_Surge')

km7<-survfit(Surv(Time,DSevent)~Partner,data=dat)
plot(km7,fun='cloglog',xlab='time in months using  logarithmic scale',ylab='log-log survival',
     main='log-log curves by Partner')

# Remove LN varaible , violate PH assumption

# Radiation
km9<-survfit(Surv(Time,DSevent)~Radiation,data=dat)
plot(km9,fun='cloglog',xlab='time in months using  logarithmic scale',ylab='log-log survival',
     main='log-log curves by Radiation')

# Chemotherapy
km10<-survfit(Surv(Time,DSevent)~Chemotherapy,data=dat)
plot(km10,fun='cloglog',xlab='time in months using  logarithmic scale',ylab='log-log survival',
     main='log-log curves by Chemotherapy')

# American Indian interact with other races ,try remove American Indian

dat_rm_AIAN<- dat[dat$Race!='Non-Hispanic American Indian/Alaska Native',]
dat_rm_AIAN<-dat_rm_AIAN[dat_rm_AIAN$Ulceration!='Unknown',]
dat_rm_AIAN<-dat_rm_AIAN[dat_rm_AIAN$PS_Surg!='Unknown',]
dat_rm_AIAN<-dat_rm_AIAN[dat_rm_AIAN$Stage!='Unkown',]
dat_rm_AIAN<-dat_rm_AIAN[dat_rm_AIAN$Histology!='Other',]
##############################
#PH COX MODEL Table Function##
##############################
PH_table_structure<-function(model){
CI = summary(model)$conf.int[,-2]
pVal = as.data.frame(summary(model)$coefficients[,5])
table = cbind(CI,pVal)
table[,4] = round(table[,4],2)
table[which(table[,4]<0.001),4]='<0.001'
table[,1:3] = round(table[,1:3],2)
names(table) = c('HR','lower 95%','upper 95%','P value')
print(table)
return(table)
}

# Unadjusted COX model stratified by Race
unadjusted_race<-coxph(Surv(Time,DSevent)~Race,data=dat_rm_AIAN)
summary(unadjusted_race)

# Adjusted COX model 
adjust_coxmodel<-coxph(Surv(Time,DSevent)~Race+sex+Age+Partner+Site+Histology+Poverty+Education+
                 Year+Stage+Ulceration+Thickness+Oth_Surg+PS_Surg+Radiation+Chemotherapy,data=dat_rm_AIAN)
summary(adjust_coxmodel)

# 5.generate COX tables ready to be used
unadjust_model_table<-PH_table_structure(unadjusted_race)
adjust_model_table<-PH_table_structure(adjust_coxmodel)
