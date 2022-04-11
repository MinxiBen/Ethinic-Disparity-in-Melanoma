library(rcompanion)
library(reshape2)
library(ggplot2)

# 1. path setting
path_split<-strsplit(getwd(),split="/")[[1]]
parent_path<-paste(path_split[1:(length(path_split)-1)],collapse="/")
data_dir<-paste(parent_path,"data",sep="/")

data_name<-"final_data_ml_20210901.RData"
data_path<-paste(data_dir,data_name,sep="/")
figure_path<-paste(parent_path,"figure",sep="/")

# 2. load data
load(data_path)
is_factor<-sapply(dat,is.factor)
factor_varnames<-names(dat)[is_factor]


#######################
# Cramer's v function #
#######################
nominal_associate<-function(dat,name){
  cramerv_mat<-matrix(0,nrow = length(name),ncol= length(name))
  for( i in 1:nrow(cramerv_mat))
    for(j in 1:ncol(cramerv_mat)){
      
      cramerv_mat[i,j] = cramerV(dat[[name[i]]],dat[[name[j]]])
      
    }
  
  colnames(cramerv_mat) = c("Sex","Race","Partner","Site","Histology","Stage","Lymph Node Status",
                            "Ulceration","Distant Metastasis","Primary Site Surgery","Other Surgery",
                            "Lymph node Surgery","Radiation","Chemotherapy")
  rownames(cramerv_mat) =  c("Sex","Race","Partner","Site","Histology","Stage","Lymph Node Status",
                             "Ulceration","Distant Metastasis","Primary Site Surgery","Other Surgery",
                             "Lymph node Surgery","Radiation","Chemotherapy")
  return(cramerv_mat)
}


####################
# Heatmap function #
####################
Heatmap_cramer<-function(cramersummary){
  cramersummary<-as.data.frame(cramersummary)
  cramersummary$ID<-rownames(cramersummary)
  cramersummary$ID<-factor(cramersummary$ID,levels=rownames(cramersummary))
  cramerV_melt<-melt(cramersummary,id.vars = c("ID"))#####long table format
  heatmap_Cramer<-ggplot(cramerV_melt,aes(x=variable,y=ID))+
    geom_tile(aes(fill=value))+
    theme(axis.text.x=element_text(angle=45,hjust=1, vjust=1))+
    scale_fill_gradient(low = "white", high = "red")+
    labs(x="Categorical Variables",y="",title="Categorical Variables Correlation")
  return(heatmap_Cramer)
}

cramer_table<-nominal_associate(dat,factor_varnames)
heatmap<-Heatmap_cramer(cramer_table)

### 3. Save picture ###
png(filename=paste(figure_path,'cramer.png',sep='/'))
heatmap
dev.off()




###################################numeric data correlaton
list =c('Year','Age','Poverty','Education','Thickness') 
dat_numeric = dat1[list]
round(cor(dat_numeric,method="pearson"),2)