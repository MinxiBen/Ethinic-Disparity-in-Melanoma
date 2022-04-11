library(ggplot2)

# 1.path setting
path_split <- strsplit(getwd(), split="/")[[1]]
parent_path <- paste(path_split[1:(length(path_split)-1)], collapse = "/")
data_path <- paste(parent_path, "data", sep = "/")
result_path <- paste(parent_path, "table", sep = "/")

load("C:/learning/R/melanoma/manusript/table/gbm_50iter.RData")
load('C:/learning/R/melanoma/manusript/table/rf_repeat/iter_50times.RData')

#########################################
#2.Produce mean ROC curve dataframe(rf,gbm) 
#########################################
roc_value<-function(df){

########cut off column for merging
cut_new=rev(seq(0,1,0.001))
cut_new =cut_new[-1001]
dat =NULL

for (j in 1:1000)
  {
  newcut = 1/2*(2*cut_new[j]-0.001)
  dat = c(dat,newcut)
  }

al = as.data.frame(dat)
colnames(al)= 'cutoff'

for (i in 1:50)
{
cut =df[[i]][[4]]$cut
fpr =df[[i]][[4]]$fpr
tpr = df[[i]][[4]]$tpr
roc = data.frame(cut=cut,fpr=fpr,tpr=tpr)
cut_new=rev(seq(0,1,0.001))
cut_new =cut_new[-1001]
newroc = NULL
for (j in 1:1000)
  {
  split = (cut_new[j]-0.001<=cut)&(cut<=cut_new[j])
  newcut = 1/2*(2*cut_new[j]-0.001)
  newfpr = median(roc$fpr[split])
  newtpr = median(roc$tpr[split])
  row = c(newcut,newfpr,newtpr)
  newroc = rbind(newroc,row)
}
a=as.data.frame(newroc)
rownames(a) = 1:1000
colnames(a) = c("cutoff",paste0("fpr",i),paste0("tpr",i))
print(a)
print(nrow(a))
al = merge(al,a,by='cutoff')

}
fpr_mean = apply(al[,seq(2,100,2)],1,mean)
tpr_mean = apply(al[,seq(3,101,2)],1,mean)
sum_al = data.frame(al[,"cutoff"],fpr_mean,tpr_mean)
colnames(sum_al)[1]="cutoff"
return(sum_al)
}

#3.mean roc table of cutoff,tpr and fpr are extracted
mean_roc_gbm<-roc_value(summary_gbm)
mean_roc_rf<-roc_value(summary_rf)
#result=na.omit(result)

save(mean_roc_gbm,file=paste(result_path,"mean_ROC_gbm.RData",sep="/"))
save(mean_roc_rf,file=paste(result_path,"mean_ROC_rf.RData",sep="/"))

#4.plot mean ROC curve of rf and gbm
ggplot(mean_roc_gbm,aes(fpr_mean,tpr_mean))+geom_line()+theme_classic()+
  labs(title = 'Mean ROC Curve of 50 Iterations',
        x ="Mean false positive rate",
        y = "Mean true positive rate")+
  annotate("text", x=0.5, y=0.5, label="Mean AUC")

ggplot(mean_roc_rf,aes(fpr_mean,tpr_mean))+geom_line()+theme_classic()+
  labs(title = 'Mean ROC Curve of 50 Iterations for Rf',
       x ="Mean false positive rate",
       y = "Mean true positive rate")+
  annotate("text", x=0.5, y=0.5, label="Mean AUC")  
#ggsave("C:/learning/R/melanoma/table/gbm/roc_gbm.eps", width = 20, height = 20, units = "cm")



#5.Combine gbm and rf mean ROC curve
mean_roc_rf$group = 0
mean_roc_gbm$group =1
all = rbind(mean_roc_rf,mean_roc_gbm)
#all = rbind(sum_al,sum_al1)

ggplot(all,aes(x=fpr_mean,y=tpr_mean,group=group,col=group))+geom_line()+theme_classic()+
  labs(title = 'RF/GBM Mean ROC Curve of 50 Iterations',
       x ="Mean false positive rate",
       y = "Mean true positive rate")+
  annotate("text", x=0.5, y=0.5, label="RF/GBM Mean AUC:0.967/0.915 ",size=5)+
  guides(color=FALSE)+
  theme(axis.text.x = element_text(size = 14,color="black"),axis.text.y = element_text(size = 13,color="black"))+  ###axis ticks
  theme(axis.title.x = element_text(size = 13),axis.title.y = element_text(size = 13))+  ###axis title
  theme(plot.title = element_text(size = 16, face = "bold")) ###plot title

#ggsave("C:/learning/R/melanoma/table/roc_gbm_rf_comb.eps", width = 20, height = 20, units = "cm")