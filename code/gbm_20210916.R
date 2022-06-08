set.seed(2222)
library(caret)
library(plyr)
library(dplyr)
library(gbm)
library(ROCR)

setwd('C:/learning/R/melanoma/manusript/table')
# 1.path setting
path_split <- strsplit(getwd(), split="/")[[1]]
parent_path <- paste(path_split[1:(length(path_split)-1)], collapse = "/")
data_path <- paste(parent_path, "data", sep = "/")
result_path <- paste(parent_path, "table", sep = "/")

#2. load data
load(paste(data_path,"final_data_ml_20210901.RData",sep = "/"))
dat1$DSevent = as.factor(ifelse(dat1$DSevent==1,'dead','alive'))

#############################################################################
#3.Up-down sampling from 0 and 1 and return train and test data list##
#############################################################################
sampleUpDown<-function(data){
  data$DSevent =as.factor(data$DSevent)
  subset1<-data[data$DSevent=='dead',]
  dead_index<-sample(1:nrow(subset1), size = nrow(subset1)*2, replace = TRUE)
  dead<-subset1[dead_index,]
  subset0<-data[data$DSevent=='alive',]
  live_index<-sample(1:nrow(subset0), size = nrow(subset0)/3, replace =FALSE)
  live<-subset0[live_index,]
  total<-rbind(dead,live)
  total<-total[sample(1:nrow(total), nrow(total)),]###random index
  inTraining<- createDataPartition(y=total$DSevent,p=0.70, list=FALSE)
  train<-total[inTraining,]
  test<-total[-inTraining,]
  alist = list(train,test)
  return(alist)
}

###################################
#4.Tune parameters for train data# 
###################################
tuneGBM<-function(gbm_train){
  
  tuneGrid<-expand.grid(####ranger only 3 tuning parameter
    interaction.depth = c(3,5,7), ###narrow para
    n.trees= c(100,200,300,400,500),
    shrinkage=c(0.1,0.2),
    #####narrow para
    n.minobsinnode=c(10)
  )
  
  fitControl<-trainControl(##5 fold cv
    method="cv",
    number=5,
    summaryFunction = twoClassSummary,##compute AUC,spec and sensi
    classProbs = TRUE, # IMPORTANT!
    verboseIter = FALSE)
  
  gbmfit<-train(
    DSevent~.,
    data = gbm_train,
    tuneGrid = tuneGrid,
    method = "gbm", 
    trControl = fitControl,
    metric= "ROC"
  )
  gbmfit_outcome = gbmfit$results
  whichTwoPct_gbm <- tolerance(gbmfit_outcome, metric = "ROC", 
                               tol = 2, maximize = TRUE)  
  best_gbm_par=gbmfit_outcome[whichTwoPct_gbm,]
  depth= as.numeric(best_gbm_par['interaction.depth'])
  n_trees = as.numeric(best_gbm_par['n.trees'])
  lr = as.numeric(best_gbm_par['shrinkage'])
  min_node_size = as.numeric(best_gbm_par['n.minobsinnode'])
  roc = as.numeric(best_gbm_par['ROC'])
  imp =varImp(gbmfit,scale = FALSE)
  imp =imp$importance
  imp = data.frame(variable=rownames(imp),importance=imp$Overall)
  np = order(imp$importance,decreasing = T)
  imp=imp[np,]
  gbm_rank = imp[1:20,]####variable importance rank table
  update(gbmfit,param = list(interaction.depth= depth,n.trees=n_trees ,shrinkage=lr,n.minobsinnode=min_node_size))
  gbm_results = list(gbmfit,c(depth,n_trees,lr,min_node_size),roc,gbm_rank)
  return(gbm_results)
}


#############
#5.Prediction#
#############
Prediction<-function(model,test){
  Pred=predict(model,newdata=test,type='prob')[,2]
  pred<-prediction(Pred,test$DSevent)
  perf<-performance(pred,"tpr","fpr")
  cutoffs <- data.frame(cut=perf@alpha.values[[1]], fpr=perf@x.values[[1]], 
                        tpr=perf@y.values[[1]])
  ####AUC
  perf1<-performance(pred,"auc")
  true = as.numeric(as.factor(test$DSevent))-1
  predict = ifelse(Pred>0.5,1,0)
  acc_gbm = (nrow(test)-sum(abs(true-predict)))/nrow(test)
  print(paste0("AUC value:",perf1@y.values[[1]]))###94.8%
  print(paste0("Accuracy:",acc_gbm))
  roc_gbm = perf1@y.values[[1]]
  DF = as.data.frame(cbind(true,predict))
  recall_gbm = sum(((DF$true)&(DF$predict)))/sum(DF$true)
  spec = sum(((DF$true==0)&(DF$predict==0)))/sum(DF$true==0)
  print(paste0("Sensitivity:",recall_gbm))
  print(paste0("Specificity:",spec))
  list_pred = c(acc_gbm,roc_gbm,recall_gbm,spec,cutoffs)
  return(list_pred)
}


##########################
#6.iterable function######
##########################
gbm_iter<- function(iterations,data_gbm){###data=dat1
  all_gbm = list()
  for (i in 1:iterations){
    print(paste0("Iteration:",i))
    alist_gbm= sampleUpDown(data_gbm)###train,test data list
    train_gbm = alist_gbm[[1]]
    test_gbm = alist_gbm[[2]]
    result_gbm = tuneGBM(train_gbm)
    print(paste0("Best parameter combination:",result_gbm[[2]]))
    print(paste0("Variable Importance:",result_gbm[[4]]))
    
    pred_list =Prediction(result_gbm[[1]], test_gbm)
    all_gbm[[i]] = list(result_gbm[[2]],result_gbm[[3]],result_gbm[[4]],pred_list)
    gc()
  }
  return(all_gbm)
}


summary_gbm <-gbm_iter(50,dat1)

save(summary_gbm,file = paste(result_path,"gbm_50iter.RData",sep = "/"))


###
#a=summary_gbm
#new2 = c()
#for (i in 1:10){
#  b=a[i]
#  b[[1]][[1]][[1]]=NULL
#  print(b)
#  new2=c(new2,b)
#}

#####new,new1,new2 represent seed500,seed588,seed688
#save(new2,file = "/home/bmx/bmx/seed688.RData")

auc_gbm = NULL
for( i in 1:50){
  auc_gbm<-c(summary_gbm[[i]][[4]][[2]],auc_gbm)
  
}
mean<-mean(auc_gbm)
std <- function(x) sd(x)/sqrt(length(x))
sd<-std(auc_gbm)
CI_lower<-mean-1.96*sd
CI_upper<-mean+1.96*sd
print(paste("(",CI_lower,",",CI_upper,")"))
#( 0.915296385858179 , 0.916145217553098 )

auc_acc = NULL
for( i in 1:50){
  auc_acc<-c(summary_gbm[[i]][[4]][[1]],auc_acc)
  
}
mean<-mean(auc_acc)
std <- function(x) sd(x)/sqrt(length(x))
sd<-std(auc_acc)
CI_lower<-mean-1.96*sd
CI_upper<-mean+1.96*sd
print(paste("(",CI_lower,",",CI_upper,")"))
#( 0.839135748952328 , 0.840344861681117 )
