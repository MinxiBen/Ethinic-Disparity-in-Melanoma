library("writexl")
library(ggplot2)
library(readxl)

# 1.path setting
path_split <- strsplit(getwd(), split="/")[[1]]
parent_path <- paste(path_split[1:(length(path_split)-1)], collapse = "/")
data_path <- paste(parent_path, "data", sep = "/")
result_path <- paste(parent_path, "table", sep = "/")

load('C:/learning/R/melanoma/manusript/table/rf_repeat/iter_50times.RData')
load("C:/learning/R/melanoma/manusript/table/gbm_50iter.RData")

#2.extract data
############################################################
##mean variable importance plot
##########################################################
imp_comb<-function(outcome){
  imp = NULL
  
  for (i in 1:50){
    #each_imp = summary_rf
    each_imp = outcome[[i]][[3]]
    each_imp=each_imp[each_imp$variable %in%c("Thickness","StageIV","Age","StageIII","UlcerationYes","Year","Education","Poverty","DMetastasisYes"),]
    each_imp= each_imp[order(each_imp$variable),]
    #if (value$variable[7]=='StageIV'){
      #value = value[c(1:6,8,7,9),]
   # }
    importance = each_imp$importance
    imp = cbind(imp,importance)
    imp <-as.data.frame(imp)
  }
  return(imp)
}

rf_imp50<-imp_comb(summary_rf)
gbm_imp50<-imp_comb(summary_gbm)

#3.1 Random forest statistic compute
rf_mean=apply(rf_imp50,1,mean)
rf_var = apply(rf_imp50,1,var)
rf_imp50$mean =rf_mean
rf_se = sqrt(rf_var)/sqrt(50)
rf_imp50$se = rf_se
rf_lowe_ci = rf_mean-qt(0.975,50-1)*rf_se
rf_upper_ci = rf_mean+qt(0.975,50-1)*rf_se
rf_imp50$lower_ci = rf_lowe_ci
rf_imp50$upper_ci = rf_upper_ci

varieble = c('Age','DMetastasisYes','Education','Poverty','StageIII','StageIV','Thickness','UlcerationYes','Year')
rf_imp50 <-as.data.frame(cbind(varieble,rf_imp50))

#3.2 GBM statistic compute
gbm_mean=apply(gbm_imp50,1,mean)
gbm_var = apply(gbm_imp50,1,var)
gbm_imp50$mean =gbm_mean
gbm_se = sqrt(gbm_var)/sqrt(50)
gbm_imp50$se = gbm_se
gbm_lowe_ci = gbm_mean-qt(0.975,50-1)*gbm_se
gbm_upper_ci = gbm_mean+qt(0.975,50-1)*gbm_se
gbm_imp50$lower_ci = gbm_lowe_ci
gbm_imp50$upper_ci = gbm_upper_ci

varieble = c('Age','Education','Poverty','StageIII','StageIV','Thickness','UlcerationYes','Year')
gbm_imp50 <-as.data.frame(cbind(varieble,gbm_imp50))

# 4.export data to rdata and excel
save(rf_imp50,file=paste(result_path,"rf_50imp.RData",sep="/"))
save(gbm_imp50,file=paste(result_path,"gbm_50imp_20210917.RData",sep="/"))
write_xlsx(rf_imp50,paste(result_path,"rf_50imp.xlsx",sep="/") )
write_xlsx(gbm_imp50,paste(result_path,"gbm_50imp_20210917.xlsx",sep="/") )



# 5.bar plot
##############################################################random forest plot genesis
rf_imp<-read_xlsx(paste(result_path,"rf_50imp.xlsx",sep="/"))
colnames(rf_imp)[1]='Variable'
p <- ggplot(rf_imp, aes(x=Variable, y=mean)) + 
  geom_bar(stat="identity", position=position_dodge(),width=0.7, fill="steelblue") +
  scale_x_discrete(limits=c("Thickness","Age","Poverty","Education","Year","StageIV","Metastasis","Ulceration","StageIII"))+
  geom_errorbar(aes(ymin=lower_ci, ymax=upper_ci), width=.2,
                position=position_dodge(.9))

plot1<-p + theme_minimal()+labs(title='Variable Importance From RF Model',x='Features',y='Importance')+
  geom_text(aes(label=round(mean,1)), vjust=1.6, color="white", size=3)+
  theme(axis.text.x = element_text(size = 13,color="black",angle = 40),axis.text.y = element_text(size = 13,color="black"))+
  theme(axis.title.x = element_text(size = 13),axis.title.y = element_text(size = 13))+
  theme(plot.title = element_text(size = 14, face = "bold"))+
  theme(panel.grid=element_blank())

#ggsave("C:/learning/R/melanoma/table/gbm/var_imp_gbm.eps", width = 20, height = 20, units = "cm")
##############################################################GBM plot genesis
gbm_imp<-read_xlsx(paste(result_path,"gbm_50imp_20210917.xlsx",sep="/"))
colnames(gbm_imp)[1]='Variable'
p <- ggplot(gbm_imp, aes(x=Variable, y=mean)) + 
  geom_bar(stat="identity", position=position_dodge(),width=0.7, fill="steelblue") +
  scale_x_discrete(limits=c("Thickness","StageIV","Age","StageIII","Year","Ulceration","Education","Poverty"))+
  geom_errorbar(aes(ymin=lower_ci, ymax=upper_ci), width=.2,
                position=position_dodge(.9))

plot2<-p + theme_minimal()+labs(title='Variable Importance From GBM Model',x='Features',y='Importance')+
  geom_text(aes(label=round(mean,1)), vjust=1.6, color="white", size=3)+
  theme(axis.text.x = element_text(size = 13,color="black",angle=40),axis.text.y = element_text(size = 13,color="black"))+
  theme(axis.title.x = element_text(size = 13),axis.title.y = element_text(size = 13))+
  theme(plot.title = element_text(size = 14, face = "bold"))+
theme(panel.grid=element_blank())

combo<-ggarrange(plot1, plot2, labels = c("A", "B"), ncol = 2, nrow = 1) 
print(combo)





























##############################################################gbm plot genesis
gbm_imp<-read.csv('C:/learning/R/melanoma/table/gbm/gbm_importance_50.csv')
gbm_imp$variable = gbm_imp$X
data = gbm_imp[,52:56]
p <- ggplot(data, aes(x=variable, y=mean)) + 
  geom_bar(stat="identity", position=position_dodge(),width=0.7, fill="steelblue") +
  scale_x_discrete(limits=c("Thickness","StageIV","Age","StageIII","Year","UlcerationYes","Poverty","Education"))+
  geom_errorbar(aes(ymin=lower_ci, ymax=upper_ci), width=.2,
                position=position_dodge(.9))

p + theme_minimal()+labs(title='Mean Feature Importance From Gradient Boosting Trees Model',x='Features',y='Importance')+
  geom_text(aes(label=round(mean,1)), vjust=1.6, color="white", size=3.5)+
  theme(axis.text.x = element_text(size = 13,color="black"),axis.text.y = element_text(size = 13,color="black"))+
  theme(axis.title.x = element_text(size = 13),axis.title.y = element_text(size = 13))+
  theme(plot.title = element_text(size = 16, face = "bold"))
