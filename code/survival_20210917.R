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

###############################################
# 4.Disease-specific Survival and PH COX Model#
###############################################
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


# Unadjusted COX model stratified by Race
unadjusted_race<-coxph(Surv(Time,DSevent)~Race,data=dat_rm_AIAN)
summary(unadjusted_race)

# Adjusted COX model 
adjust_coxmodel<-coxph(Surv(Time,DSevent)~Race+sex+Age+Partner+Site+Histology+Poverty+Education+
                 Year+Stage+Ulceration+Thickness+Oth_Surg+PS_Surg+Radiation+Chemotherapy,data=dat_rm_AIAN)
summary(adjust_coxmodel)


########
# violation of PH 
# American Indian interact with other races ,try remove American Indian

dat_rm_AIAN<- dat[dat$Race!='Non-Hispanic American Indian/Alaska Native',]
dat_rm_AIAN<-dat_rm_AIAN[dat_rm_AIAN$Ulceration!='Unknown',]
dat_rm_AIAN<-dat_rm_AIAN[dat_rm_AIAN$PS_Surg!='Unknown',]
dat_rm_AIAN<-dat_rm_AIAN[dat_rm_AIAN$Stage!='Unkown',]
dat_rm_AIAN<-dat_rm_AIAN[dat_rm_AIAN$Histology!='Other',]

km<-survfit(Surv(Time,DSevent)~Race,data=dat_rm_AIAN)
plot(km,fun='cloglog',xlab='time in months using  logarithmic scale',ylab='log-log survival',
     main='log-log curves by Race')

# need remove ulceration unknown
km1<-survfit(Surv(Time,DSevent)~Ulceration,data=dat_rm_AIAN)
plot(km1,fun='cloglog',xlab='time in months using  logarithmic scale',ylab='log-log survival',
     main='log-log curves by Ulceration')

km2<-survfit(Surv(Time,DSevent)~Stage,data=dat_rm_AIAN)
plot(km2,fun='cloglog',xlab='time in months using  logarithmic scale',ylab='log-log survival',
     main='log-log curves by Stage')

# 
#dat_rm_stage<- dat[dat$!='Unkown',]
km3<-survfit(Surv(Time,DSevent)~PS_Surg,data=dat_rm_AIAN)
plot(km3,fun='cloglog',xlab='time in months using  logarithmic scale',ylab='log-log survival',
     main='log-log curves by PS_Surge')

# Other_Surg
km4<-survfit(Surv(Time,DSevent)~Oth_Surg,data=dat_rm_AIAN)
plot(km4,fun='cloglog',xlab='time in months using  logarithmic scale',ylab='log-log survival',
     main='log-log curves by Oth_Surge')

# Histoly remove OTHER and Melanoma NOS
km5<-survfit(Surv(Time,DSevent)~Histology,data=dat_rm_AIAN)
plot(km5,fun='cloglog',xlab='time in months using  logarithmic scale',ylab='log-log survival',
     main='log-log curves by LN_Surge')


#dat_rm_Site<- dat[dat$Site!='Unkown',]
km6<-survfit(Surv(Time,DSevent)~Site,data=dat_rm_AIAN)
plot(km6,fun='cloglog',xlab='time in months using  logarithmic scale',ylab='log-log survival',
     main='log-log curves by LN_Surge')


#dat_rm_Site<- dat[dat$Site!='Unkown',]
km7<-survfit(Surv(Time,DSevent)~Partner,data=dat_rm_AIAN)
plot(km7,fun='cloglog',xlab='time in months using  logarithmic scale',ylab='log-log survival',
     main='log-log curves by Partner')

# Remove LN varaible , violate PH assumption

# Radiation
km9<-survfit(Surv(Time,DSevent)~Radiation,data=dat_rm_AIAN)
plot(km9,fun='cloglog',xlab='time in months using  logarithmic scale',ylab='log-log survival',
     main='log-log curves by Radiation')

# Chemotherapy
km10<-survfit(Surv(Time,DSevent)~Chemotherapy,data=dat_rm_AIAN)
plot(km10,fun='cloglog',xlab='time in months using  logarithmic scale',ylab='log-log survival',
     main='log-log curves by Chemotherapy')


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

# 5.generate COX tables ready to be used
unadjust_model_table<-PH_table_structure(unadjusted_race)
adjust_model_table<-PH_table_structure(adjust_coxmodel)



########################
#6. Time-extended model#
########################
# Data Preprocess
dat$id <-sur_data$`Patient ID`
dat[,c("Breslow","N_RN_RM","RLN_Surg","M_stage","LN_Metastasis","Stage_bin","OSevent")]<-NULL

# Convert categorical variables into dummy variables
is_factor<-sapply(dat,is.factor)
factor_varnames<-names(dat)[is_factor]

dummy_all = list()
for (i in 1:14)
  {
  var = factor_varnames[i]
  dummy_all[[i]]<-with(dat, as.data.frame(model.matrix(~get(var))[,-1],col.names=levels(get(var))[-1]))
  }
for(i in 1:14){dat = cbind(dat,dummy_all[[i]])}

colnames(dat)[23:62]=c("Male",
                        "NHB","NHAPI","NHAA","Hispanics",
                        "Unpartner","PartnerUnknown",
                        "Lower","Trunk","Headneck","Overlapping","SiteUnknown",
                        "LM","Other","Desmoplastic","NOS","Amelanotic","Nodular","ALM",
                        "StageII","StageIII","StageIV","StageUnknown",
                        "LNN1","LNN2","LNN3","LNUnknown",
                        "UlcerationUnknown","UlcerationYes",
                        "DMetastasisUnknown","DMetastasisYes",
                        "PS_SurgYes","PS_SurgUnknown",
                        "Oth_SurgPerformed","Oth_SurgUnknown",
                        "N_SurgBiopsy","LN_SurgRemoved","LN_SurgUnknown",
                        "Radiation",
                        "Chemotherapy")

dat[,factor_varnames] = NULL

# time split group by patient id
set.seed(101)
data_split = dat[sample(1:nrow(dat),size=nrow(dat)/2),]
data_split<-data_split[data_split$Time!=0,]
ext_data<-survSplit(data_split,cut = data_split$Time[data_split$DSevent==1],
                                        end="Time",event="DSevent",start="start")
# make time-violation variables interact with log(time)
inter_term<-ext_data[!colnames(ext_data)%in%c("Year","Poverty","Education","NHB","NHAPI","NHAA","Hispanics","start","Time","DSevent","id")]*log(ext_data$Time)
ext_data<-cbind(ext_data,inter_term)
for (i in 50:87){
names(ext_data)[i]<-paste(names(ext_data[i]),"t",sep="")
}

# construct time-extended model 
vars <-names(ext_data)[which(!names(ext_data)%in%c('start','Time','DSevent','id'))]
fomula<-as.formula(paste(paste("Surv(ext_data$start,ext_data$Time,ext_data$DSevent)~", paste(vars, collapse="+")),'+cluster(id)'))
time_extended_model<-coxph(fomula,data=ext_data)
print(summary(time_extended_model))










