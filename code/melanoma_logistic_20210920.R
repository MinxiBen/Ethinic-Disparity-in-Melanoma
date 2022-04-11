library(ggplot2)


#1. load data
load("C:/learning/R/melanoma/manusript/data/final_data_ml_20210901.RData")

# 2.logistic interaction regression 
#####################################################
####rebuild logistic model functions 2020/10/12
####################################################
#######Race,Ulceration interaction
outcome <- "DSevent"
x<-c("Year","Age","sex","Partner","Site","Histology","Poverty","Education",
             "Stage","LN","Thickness","DMetastasis","PS_Surg","Oth_Surg",
             "LN_Surg","Radiation","Chemotherapy")
inter_term <-c("Race","Ulceration")

######Race,Thickness interaction
outcome <- "DSevent"
x<-c("Year","Age","sex","Partner","Site","Histology","Poverty","Education",
     "Stage","LN","Ulceration","DMetastasis","PS_Surg","Oth_Surg",
     "LN_Surg","Radiation","Chemotherapy")
inter_term <-c("Race","Breslow")

####Race, Dmetastasis interaction
outcome <- "DSevent"
x<-c("Year","Age","sex","Partner","Site","Histology","Poverty","Education",
     "Stage","LN","Thickness","Ulceration","PS_Surg","Oth_Surg",
     "LN_Surg","Radiation","Chemotherapy")
inter_term <-c("Race","DMetastasis")

##############Stage,Race interaction
dat1$stage = 0
dat1$stage[dat1$Stage%in%c('III','IV')]=1
dat1$Stage = NULL

outcome <- "DSevent"
x<-c("Year","Age","sex","Partner","Site","Histology","Poverty","Education",
     "Ulceration","LN","Thickness","DMetastasis","PS_Surg","Oth_Surg",
     "LN_Surg","Radiation","Chemotherapy")
inter_term <-c("Race","stage")


# 3. logistic regression abstract interaction terms
###########################################################################
##Function(Race_interact): get logistic coefficients and covariances matrix
###parameter:
###inter_y:outcome
###inter_x:variables withour interactions
###inter_term: varaibles which needs to caculate interaction
############################################################################
Race_interact<-function(inter_y,inter_x,interaction_term,inter_data)
  {
 interaction = as.character(paste(interaction_term,collapse = "*"))
 variables<-c(inter_x,interaction)
 f<-as.formula(paste(inter_y,paste(variables,collapse = "+"),sep="~"))
 print(f)
 logistic_model<-glm(f,family=binomial(link='logit'),data=inter_data)
#logistic_model<-glm(DSevent~LN_Surg+PS_Surg+Oth_Surg+Radiation+Chemotherapy+LN
           ## +Thickness+DMetastasis+Age+Site+Partner+Poverty+Education+sex+Histology+Race*Ulceration,family=binomial(link='logit'),data=data)
 summary(logistic_model)
 modvcov = vcov(logistic_model)###variance covariance plot
 raceName = rownames(summary(logistic_model)$coeff)[grep('Race',rownames(summary(logistic_model)$coefficients))]
 #raceName=raceName[c(1,2,3,4,9,10,11,12)]###choose the interaction terms you need
 
 modcvoc_race=modvcov[raceName,raceName]
 raceCoefs = as.data.frame(summary(logistic_model)$coeff[raceName,])
 print(modcvoc_race)
 print(raceCoefs)
 return (list(raceCoefs,modcvoc_race))#####return list of interaction coefficients and covariances
}
####fetch results
result_inter = Race_interact(outcome,x,inter_term,dat1)
result_int_coeff = result_inter[[1]]
result_int_cov =result_inter[[2]]


model1<-glm(DSevent~LN_Surg+PS_Surg+Oth_Surg+Radiation+Chemotherapy+LN+Thickness+DMetastasis+Age+Site+Partner+Poverty+Education+sex+Histology+Race+Ulceration,family=binomial(link='logit'),data=dat1)
model2<-glm(DSevent~LN_Surg+PS_Surg+Oth_Surg+Radiation+Chemotherapy+LN +Thickness+DMetastasis+Age+Site+Partner+Poverty+Education+sex+Histology+Race*Ulceration,family=binomial(link='logit'),data=dat1)
anova(model1,model2,test='LRT')##0.01651 

#4. summary interaction table
###########################################################################
###Function(table of Race interaction):table of race interaction with the feature
###########################################################################
####stage*race
Table_Race_interact<-function(inter_coeff,inter_cov){
  raceinter = data.frame(Estimates = c(inter_coeff$Estimate[1:4],inter_coeff$Estimate[1:4]+inter_coeff$Estimate[5:8]),
                                                                 #inter_coeff$Estimate[1:4]+inter_coeff$Estimate[9:12]),
                                                                 #inter_coeff$Estimate[1:4]+inter_coeff$Estimate[13:16]),
                          stE = sqrt(c(diag(inter_cov[1:4,1:4]),
                                       diag(inter_cov[1:4,1:4])+diag(inter_cov[5:8,5:8])+2*c(diag(inter_cov[1:4,5:8]))
                                       #diag(inter_cov[1:4,1:4])+diag(inter_cov[9:12,9:12])+2*c(diag(inter_cov[1:4,9:12]))
                                       #diag(inter_cov[1:4,1:4])+diag(inter_cov[13:16,13:16])+2*c(diag(inter_cov[1:4,13:16]))
                                       )))###get interCTION term logHR and selogHR
  raceinter$Lb = raceinter$Estimates-qnorm(0.975)*raceinter$stE###lower ci
  raceinter$Ub = raceinter$Estimates+qnorm(0.975)*raceinter$stE###upper ci
  raceinter$Race = rep(c('NHB','NHAPI','NHAA','Hispanics'),2)
  raceinter$Variable  = rep(c('I/II','III/IV'),each=4)
  #raceinter$Variable  = rep(c('No','Yes'),each=4)
  #raceinter$Variable  = rep(c('I','II','III','IV'),each=4)
  race_inter_table<-data.frame(Variable=raceinter$Variable,Race = raceinter$Race,LogOR=raceinter$Estimates,Se=raceinter$stE,
                    LB=raceinter$Lb,UB=raceinter$Ub,OR = exp(raceinter$Estimates),ORLB=exp(raceinter$Lb),ORUB=exp(raceinter$Ub))
  
  print(race_inter_table)
  return(race_inter_table)
  
}

# Ulceration(Dmetastasis)*Race
Table_Race_interact<-function(inter_coeff,inter_cov){
  raceinter = data.frame(Estimates = c(inter_coeff$Estimate[1:4],#inter_coeff$Estimate[1:4]+inter_coeff$Estimate[5:8]),
                         inter_coeff$Estimate[1:4]+inter_coeff$Estimate[9:12]),
                         #inter_coeff$Estimate[1:4]+inter_coeff$Estimate[13:16]),
                         stE = sqrt(c(diag(inter_cov[1:4,1:4]),
                                      #diag(inter_cov[1:4,1:4])+diag(inter_cov[5:8,5:8])+2*c(diag(inter_cov[1:4,5:8]))
                                      diag(inter_cov[1:4,1:4])+diag(inter_cov[9:12,9:12])+2*c(diag(inter_cov[1:4,9:12]))
                                      #diag(inter_cov[1:4,1:4])+diag(inter_cov[13:16,13:16])+2*c(diag(inter_cov[1:4,13:16]))
                         )))###get interCTION term logHR and selogHR
  raceinter$Lb = raceinter$Estimates-qnorm(0.975)*raceinter$stE###lower ci
  raceinter$Ub = raceinter$Estimates+qnorm(0.975)*raceinter$stE###upper ci
  raceinter$Race = rep(c('NHB','NHAPI','NHAA','Hispanics'),2)
  #raceinter$Variable  = rep(c('I/II','III/IV'),each=4)
  raceinter$Variable  = rep(c('No','Yes'),each=4)
  #raceinter$Variable  = rep(c('I','II','III','IV'),each=4)
  race_inter_table<-data.frame(Variable=raceinter$Variable,Race = raceinter$Race,LogOR=raceinter$Estimates,Se=raceinter$stE,
                               LB=raceinter$Lb,UB=raceinter$Ub,OR = exp(raceinter$Estimates),ORLB=exp(raceinter$Lb),ORUB=exp(raceinter$Ub))
  
  print(race_inter_table)
  return(race_inter_table)
  }

####call Table_Race_interact
inter_table<-Table_Race_interact(result_int_coeff,result_int_cov)
##############################################################################
########Function:draw interaction pictures
#######parameter: interaction table, interaction varaible name
#############################################################################
Plot_Inter<-function(inter_tbl,var_name){
  
  inter_plot<-ggplot(inter_tbl,aes(x=Variable,y = round(OR,2),group=Race,color=Race))+
    geom_line(size=2)+
    geom_hline(yintercept = 1,alpha=1,lty=2,size=2)+
    geom_errorbar(aes(ymin=exp(LB),ymax=exp(UB)),width=0.1,alpha=1,lty=2,size=1)+
    geom_text(aes(label=paste0('(',round(OR,2),')')),position = position_dodge(width = 0.5))+
    theme_bw()+
    theme(axis.text.y  = element_text(size = 15),axis.text.x = element_text(size = 15),
          axis.title.x = element_text(size = 15),
          plot.title = element_text(size = 15, face = "bold",hjust = 0.5))+
    labs(x=var_name,y='Odds Ratio',title='OR of other race vs whites at stage I/II & III/V')+
    scale_colour_discrete(labels = c('Hispanics','American Indian','Asian','American African'))+
    scale_color_manual(values=c("#66CC33","#FF0000","#FFCC33","#330000"))
    print(inter_plot)
    return(inter_plot)
  
 }

Plot_Inter(inter_table,"Distant Metastasis")
Plot_Inter(inter_table,"Breslow Thickness")
p1<-Plot_Inter(inter_table,"Ulceration")
p2<-Plot_Inter(inter_table,"Stage")

combo<-ggarrange(p1, p2,common.legend = TRUE,legend = "right", labels = c("A", "B"), ncol = 2, nrow = 1) 
print(combo)

# ggsave("C:/learning/R/melanoma/figure/race_thick_int.eps", width = 30, height = 20, units = "cm")
# ggsave("C:/learning/R/melanoma/figure/race_DM_int.eps", width =30, height = 20, units = "cm")














# depracate!!
inter_plot<-ggplot(inter_table,aes(x=Variable,y = round(OR,2),group=Race,color=Race))+
  geom_line(size=2)+
  geom_hline(yintercept = 1,alpha=1,lty=2,size=2)+
  geom_errorbar(aes(ymin=exp(LB),ymax=exp(UB)),width=0.1,alpha=1,lty=2,size=1)+
  geom_text(aes(label=paste0('(',round(OR,2),')')),position = position_dodge(width = 0.5))+
  theme_bw()+
  theme(axis.text.y  = element_text(size = 15),axis.text.x = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        plot.title = element_text(size = 16, face = "bold"))+
  labs(x='Ulceration',y='Odds Ratio',title='Odds ratio of death of other races relative to NHW')+
  scale_colour_discrete(labels = c('Hispanics','American Indian','Asian','American African'))+
  scale_color_manual(values=c("#66CC33","#FF0000","#FFCC33","#330000"))


print(inter_plot)

inter_plot2<-ggplot(inter_table,aes(x=Variable,y = round(OR,2),group=Race,color=Race))+
  geom_line(size=2)+
  geom_hline(yintercept = 1,alpha=1,lty=2,size=2)+
  geom_errorbar(aes(ymin=exp(LB),ymax=exp(UB)),width=0.1,alpha=1,lty=2,size=1)+
  geom_text(aes(label=paste0('(',round(OR,2),')')),position = position_dodge(width = 0.5))+
  theme_bw()+
  theme(axis.text.y  = element_text(size = 15),axis.text.x = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        plot.title = element_text(size = 16, face = "bold"))+
  labs(x='Distant Metastasis',y='Odds Ratio',title='Odds ratio of death of other races relative to NHW')+
  scale_colour_discrete(labels = c('Hispanics','American Indian','Asian','American African'))+
  scale_color_manual(values=c("#66CC33","#FF0000","#FFCC33","#330000"))


print(inter_plot2)


inter_plot3<-ggplot(inter_table,aes(x=Variable,y = round(OR,2),group=Race,color=Race))+
  geom_line(size=2)+
  geom_hline(yintercept = 1,alpha=1,lty=2,size=2)+
  geom_errorbar(aes(ymin=exp(LB),ymax=exp(UB)),width=0.1,alpha=1,lty=2,size=1)+
  geom_text(aes(label=paste0('(',round(OR,2),')')),position = position_dodge(width = 0.5))+
  theme_bw()+
  theme(axis.text.y  = element_text(size = 15),axis.text.x = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        plot.title = element_text(size = 16, face = "bold"))+
  labs(x='Stage',y='Odds Ratio',title='Odds ratio of death of other races relative to NHW')+
  scale_colour_discrete(labels = c('Hispanics','American Indian','Asian','American African'))+
  scale_color_manual(values=c("#66CC33","#FF0000","#FFCC33","#330000"))


print(inter_plot3)









######data mutate for thickness and race interaction
a<-glm(DSevent ~ Year + Age + sex + Partner + Site + Histology + Poverty + 
         Education + Stage + LN + Thickness + Ulceration + PS_Surg + 
         Oth_Surg + LN_Surg + Radiation + Chemotherapy + Race * DMetastasis,family=binomial(link='logit'),data=dat1)
assignthick<-function(x){
  ret = rep(NA,length(x))
  ret[x<1]='I'
  ret[x>=1&x<2]='II'
  ret[x>=2&x<4]='III'
  ret[x>=4] ='IV'
  return(as.factor(ret))
}
dat1$Breslow<-assignthick(dat1$Thickness)
