#chapter5
library(tableone)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(tidyverse)
library(reshape2)

#1.load data
load("C:/learning/R/melanoma/manusript/data/final_data_survival_20210901.RData")
load("C:/learning/R/melanoma/manusript/data/final_data_ml_20210901.RData")

############################
#2.Create demographic table#
############################
factorVars =c('Histology','Site','Partner','sex','Race','Stage','Ulceration','DMetastasis',
              'PS_Surg','LN_Surg','Oth_Surg','Radiation','Chemotherapy')
vars = c('sex','Age','Partner','Histology','Site','Poverty','Thickness','LN','Ulceration','Stage','DMetastasis',
         'PS_Surg','LN_Surg','Oth_Surg','Radiation','Chemotherapy')
tbl = CreateTableOne(vars=vars, strata = 'Race',factorVars = factorVars,includeNA = F,
                     data=dat1)

table = print(tbl)

#write.csv(table,'C:/learning/R/melanoma/description.csv')


###################################################################
#3.Density Plot Race vs distribution of poverty,education,thickness and age#
###################################################################
figure_plot <- function(x,y){
  ggplot(dat,aes(x=get(x),fill=get(y)))+
    geom_density(alpha=1)+
    #xlim(0,110)+
    #xlim(0,10)+
    xlim(0,50)+
    labs(
      #x='Age at Diagnosis',
      #x='Thicknes(mm)',
      x = 'Poverty Proportion(%)',
      y='Density',
      #title='Distribution of Age at Diagnosis by Race/Ethnicity',
      title = "Distribution of Poverty proportion by Race/Ethnicity",
      #title = "Distribution of Tumor thickness by Race/Ethnicity",
      fill='Race')+
    theme_bw()+
    theme(panel.border = element_blank(),panel.grid.major =element_blank(),
          panel.grid.minor = element_blank(),axis.line = element_line(colour="black"))+
    scale_fill_manual(values=c("#3366FF","#330000","#FFCC33","#66CC33","#FF0000"))
  
}

figure_plot('Age','Race')
figure_plot('Thickness','Race')
figure_plot('Poverty','Race')
#figure_plot('Education','Race')


######################################################################
#4.Box plot of Thickness,poverty with kruskal-wallis and wilcoxn text#
######################################################################
###4.1 thickness
mycomparisons =list(c("Non-Hispanic White","Non-Hispanic Black"),
                    c("Non-Hispanic White","Non-Hispanic Asian or Pacific Islander"),
                    c("Non-Hispanic Black","Non-Hispanic Asian or Pacific Islander"))
#get summary stat mean,se
dat %>%
  group_by(Race)%>%
  get_summary_stats(Thickness,type='common')
#plot without test label
thickness_boxplot<-ggboxplot(dat,x='Race',y='Thickness',color = "Race", palette = "jco")+
  #ylim(0,5)+
  scale_x_discrete(labels=c('NHWs','Blacks','API','AI/AN','Hispanics'))+
  stat_compare_means(comparison=mycomparisons)+
  stat_compare_means()

thickness_boxplot<-ggpar(thickness_boxplot,
          ylab = "Thickness(mm)",
          legend="none",
          title="Tumor thickness by Race/ethnicity")

#4.2 poverty
mycomparisons =list(c("Non-Hispanic White","Non-Hispanic Black"),
                    c("Non-Hispanic White","Non-Hispanic Asian or Pacific Islander"),
                    c("Non-Hispanic Black","Non-Hispanic Asian or Pacific Islander"))
dat %>%
  group_by(Race)%>%
  get_summary_stats(Poverty,type='common')
###plot without test label
poverty_boxplot<-ggboxplot(dat,x='Race',y='Poverty',color = "Race", palette = "jco")+
  scale_x_discrete(labels=c('NHWs','Blacks','API','AI/AN','Hispanics'))+
  stat_compare_means(comparison=mycomparisons)+
  stat_compare_means()

poverty_boxplot<-ggpar(poverty_boxplot,
          ylab = "Proportion below poverty line(%)",
          legend="none",
          title="Poverty proportion by Race/ethnicity")



################################################
##5.Bar plot for Ulceration,Stage versus Race
#################################################
##function of bar plot
race_cat<-function(x,y){
table<-xtabs(~get(x)+get(y),data=dat)
#table<-table[-2,]
table<-table[rownames(table)!='Unkown',]
#table<-table[rownames(table)!='Unknown',]
chi_test<-chisq_test(table)
prop_tb<-prop.table(table,2)
temp<-as.data.frame(t(prop_tb))
temp$Freq<-round(temp$Freq,2)
names(temp)=c('Race','Feat','Proportion')
for (i in 1:5){
  temp[i,3] = 1-temp[i+5,3]-temp[i+10,3]
} 
plot<-ggbarplot(temp,
          x='Race',
          y='Proportion',
          fill='Feat',
          color = 'Feat',
          palette='jco',
          label=TRUE,
          lab.col='Black',lab.pos='out',lab.size=4)+
  labs(subtitle = get_test_label(chi_test,detailed = TRUE),
       title = "Ulceration proportion by Race/Ethnicity",
       #title = "Stage proportion by Race/Ethnicity",
       fill='Ulceration',
       color = 'Ulceration'
       #fill='Stage'
       #fill="Metastasis"
       )+
  scale_x_discrete(labels=c('White','Black','Asian','American Indian','Hispanics'))
 return(plot)
}
plot1<-race_cat('Ulceration','Race')
plot2<-race_cat('Stage','Race')
#race_cat('DMetastasis','Race')


table<-xtabs(~get('Ulceration')+get('Race'),data=dat)
#table<-table[-2,]
table<-table[rownames(table)!='Unkown',]
#table<-table[rownames(table)!='Unknown',]
chi_test<-chisq_test(table)
prop_tb<-prop.table(table,2)
temp<-as.data.frame(t(prop_tb))
temp$Freq<-round(temp$Freq,2)
names(temp)=c('Race','Feat','Proportion')
for (i in 1:5){
  temp[i,3] = 1-temp[i+5,3]-temp[i+10,3]
} 
plot<-ggbarplot(temp,
                x='Race',
                y='Proportion',
                fill='Feat',
                color = 'Feat',
                palette='jco',
                label=TRUE,
                lab.col='Black',lab.pos='out',lab.size=4)+
  labs(subtitle = get_test_label(chi_test,detailed = TRUE),
       title = "Ulceration proportion by Race/Ethnicity",
       #title = "Stage proportion by Race/Ethnicity",
       fill='Ulceration',
       color = 'Ulceration'
       #fill='Stage'
       #fill="Metastasis"
  )+
  scale_x_discrete(labels=c('NHWs','Blacks','API','AI/AN','Hispanics'))



table<-xtabs(~get('Stage')+get('Race'),data=dat)
#table<-table[-2,]
table<-table[rownames(table)!='Unkown',]
#table<-table[rownames(table)!='Unknown',]
chi_test<-chisq_test(table)
prop_tb<-prop.table(table,2)
temp<-as.data.frame(t(prop_tb))
temp$Freq<-round(temp$Freq,2)
names(temp)=c('Race','Feat','Proportion')
for (i in 1:5){
  temp[i,3] = 1-temp[i+5,3]-temp[i+10,3]-temp[i+15,3]
} 
plot1<-ggbarplot(temp,
                x='Race',
                y='Proportion',
                fill='Feat',
                color = 'Feat',
                palette='jco',
                label=TRUE,
                lab.col='Black',lab.pos='out',lab.size=4)+
  labs(subtitle = get_test_label(chi_test,detailed = TRUE),
       #title = "Ulceration proportion by Race/Ethnicity",
       title = "Stage proportion by Race/Ethnicity",
       #fill='Ulceration',
       color = 'Stage',
       fill='Stage'
       #fill="Metastasis"
  )+
  scale_x_discrete(labels=c('NHWs','Blacks','API','AI/AN','Hispanics'))











ggarrange(thickness_boxplot,poverty_boxplot,plot,plot1,common.legend=FALSE, labels = c("A", "B", "C","D"), ncol = 2, nrow = 2) 





