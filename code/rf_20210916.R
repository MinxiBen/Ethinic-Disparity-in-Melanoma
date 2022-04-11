set.seed(1111)
library(randomForest)
library(ranger)
library(dplyr)
library(ROCR)
library(ggplot2)
library(caret)


#1. load data
load("C:/learning/R/melanoma/manusript/data/final_data_ml_20210901.RData")
dat1$DSevent = as.factor(ifelse(dat1$DSevent==1,'dead','alive'))
#############################################################################
#2.Up-down sampling from 0 and 1 and return train and test data list##
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
#3.Tune parameters for train data# 
###################################
tuneRF<-function(dat){
  
  tuneGrid<-expand.grid(####ranger only 3 tuning parameter
    .mtry = c(5:10), ###narrow para
    .min.node.size = c(1:3),#####narrow para
    .splitrule='gini'
  )
  model <-train(  
    DSevent~.,
    tuneGrid = tuneGrid,
    data = dat,
    method = "ranger",
    importance = "impurity",
    metric = "ROC",#######change from accuracy
    trControl =trainControl(
      method = "cv",
      number = 5,
      classProbs= TRUE,  ####based on ROC
      verboseIter = FALSE,
      summaryFunction = twoClassSummary,#
    )
  )
  rffit = model$results
  whichTwoPct <- tolerance(rffit, metric = "ROC", 
                           tol = 2, maximize = TRUE)  
  best=rffit[whichTwoPct,]
  mtry = as.numeric(best['mtry'])
  nsize = as.numeric(best['min.node.size'])
  avg_roc = as.numeric(best['ROC'])
  imp =varImp(model)$importance
  imp = data.frame(variable=rownames(imp),importance=imp$Overall)
  np = order(imp$importance,decreasing = T)
  imp=imp[np,]
  rank = imp[1:20,]
  update(model, param = list(mtry=mtry,min.node.size=nsize,splitrule='gini'))
  results = list(model,c(mtry,nsize),avg_roc,rank)
  return(results)
}

#############
#4.Prediction##
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
  acc = (nrow(test)-sum(abs(true-predict)))/nrow(test)
  print(paste0("AUC value:",perf1@y.values[[1]]))###94.8%
  print(paste0("Accuracy:",acc))
  roc_rf = perf1@y.values[[1]]
  DF = as.data.frame(cbind(true,predict))
  recall_rf = sum(((DF$true)&(DF$predict)))/sum(DF$true)
  spec = sum(((DF$true==0)&(DF$predict==0)))/sum(DF$true==0)
  print(paste0("Sensitivity:",recall_rf))
  print(paste0("Specificity:",spec))
  list_pred = c(acc,roc_rf,recall_rf ,spec,cutoffs)
  return(list_pred)
}

##########################
#5.iterable function######
##########################
rf_iter<- function(iterations,data_rf){###data= dat1
  all_rf = list()
  for (i in 1:iterations){
    print(paste0("Iteration:",i))
    alist_rf= sampleUpDown(data_rf)###train,test data list
    train_rf = alist_rf[[1]]
    test_rf = alist_rf[[2]]
    result_rf = tuneRF(train_rf)
    print(paste0("Best parameter combination:",result_rf[[2]]))
    print(paste0("Variable Importance:",result_rf[[4]]))
    
    # result_rf[[1]] is fi_model
    pred_list =Prediction(result_rf[[1]], test_rf)
    ##remove model part
    all_rf[[i]] = list(result_rf[[2]], result_rf[[3]], result_rf[[4]], pred_list)
    gc()
  }
  return(all_rf)
}

summary_rf <- rf_iter(50,dat1)

save(summary_rf,file='C:/learning/R/melanoma/manusript/table/rf_repeat/iter_50times.RData') 




#########average AUC , ACC
auc_rf = NULL
for( i in 1:50){
  auc_rf<-c(summary_rf[[i]][[4]][[2]],auc_rf)
  
}
mean<-mean(auc_rf)
std <- function(x) sd(x)/sqrt(length(x))
sd<-std(auc_rf)
CI_lower<-mean-1.96*sd
CI_upper<-mean+1.96*sd
print(paste("(",CI_lower,",",CI_upper,")"))
#( 0.966762393181308 , 0.967441875203871 )

rf_acc = NULL
for( i in 1:50){
  rf_acc<-c(summary_rf[[i]][[4]][[1]],rf_acc)
  
}
mean<-mean(rf_acc)
std <- function(x) sd(x)/sqrt(length(x))
sd<-std(rf_acc)
CI_lower<-mean-1.96*sd
CI_upper<-mean+1.96*sd
print(paste("(",CI_lower,",",CI_upper,")"))
#( 0.839135748952328 , 0.840344861681117 )





