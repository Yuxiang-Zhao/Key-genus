

###0. otu表导入
setwd("F:/data7") 
otu <- read.delim("otu.txt", row.names = 1) 
otu <-t(as.matrix(otu)) 
otu.percent <- otu/rowSums(otu)  


###1.‘overall abundant OTUs’:mean relative abundance>0.1%
overall.abundant.OTUs = NULL
for (i in 1:ncol(otu.percent)) {
  if(mean(otu.percent[,i])>=0.001) 
    overall.abundant.OTUs <- rbind(overall.abundant.OTUs, colnames(otu.percent)[i])
} #筛选出mean relative abundance>0.1%的otu
overall.abundant.OTUs <- data.frame(overall.abundant.OTUs) 
names(overall.abundant.OTUs) <- 'overall abundant OTUs' 


###2.‘ubiquitous OTUs’ :occurrence frequency in >80% of all samples。
otu.percent.temp<-otu.percent
otu.percent.temp[otu.percent.temp!=0]<-1 
otu.percent.temp <- otu.percent.temp[,colSums(otu.percent.temp)/nrow(otu.percent.temp)>=0.8]
ubiquitous.OTUs<-colnames(otu.percent.temp) 
ubiquitous.OTUs <- data.frame(ubiquitous.OTUs) 
names(ubiquitous.OTUs) <- 'ubiquitous OTUs' 


###3.‘frequently abundant OTUs’
otu.percent.temp<-t(otu.percent) 
abundant.otu=NULL 
for (i in 1:ncol(otu.percent.temp)) { 
  otu.percent.temp<-otu.percent.temp[order(otu.percent.temp[,i], decreasing=T),]
  count=0 
  for (j in 1:nrow(otu.percent.temp)){ 
    count=count+otu.percent.temp[j,i]  
    abundant.otu<- rbind(abundant.otu, rownames(otu.percent.temp)[j])
    if(count>=0.8) 
      break
  }}
abundant.otu.unique<-unique(abundant.otu) 
frequently.abundant.OTUs <- NULL 
for(i in 1:nrow(abundant.otu.unique)){ 
  count=0 
  for(j in 1:nrow(abundant.otu)){ 
    if(abundant.otu.unique[i]==abundant.otu[j]) 
      count= count+1   
  }
  if(count>=(ncol(otu.percent.temp)*0.5)) 
    frequently.abundant.OTUs<- rbind(frequently.abundant.OTUs, abundant.otu.unique[i])
}
frequently.abundant.OTUs <- data.frame(frequently.abundant.OTUs) 
names(frequently.abundant.OTUs) <- 'frequently abundant OTUs'


###4.判断核心物种，导出到文件
maxlength<-max(nrow(overall.abundant.OTUs),nrow(ubiquitous.OTUs),nrow(frequently.abundant.OTUs))
while(nrow(overall.abundant.OTUs)< maxlength){
  overall.abundant.OTUs<-rbind(overall.abundant.OTUs,NA)
} 
while(nrow(ubiquitous.OTUs)< maxlength){
  ubiquitous.OTUs<-rbind(ubiquitous.OTUs,NA)
}
while(nrow(frequently.abundant.OTUs)< maxlength){
  frequently.abundant.OTUs<-rbind(frequently.abundant.OTUs,NA)
}  #将3个数据框以NA填充到最大行数maxlength，从而使3个数据框长度一致
potential_core_otu=data.frame(overall.abundant.OTUs,ubiquitous.OTUs,frequently.abundant.OTUs)
core_otu=NULL
for(i in 1:nrow(potential_core_otu)){
  if((potential_core_otu$overall.abundant.OTUs[i] %in% potential_core_otu$frequently.abundant.OTUs)& 
     (potential_core_otu$overall.abundant.OTUs[i] %in% potential_core_otu$ubiquitous.OTUs))
    core_otu=rbind(core_otu,as.character(potential_core_otu$overall.abundant.OTUs[i]) ) 
}
colnames(core_otu)="core otu"
write.table(core_otu,'core_otu.txt',
            sep = "\t",quote = FALSE,row.names = F,col.names =T) #导出数据


##########附录#########################附录###################附录##########################
########################### 根据数据分布设定3项判断标准的设定 ##############################

###S0.频率、平均丰度数据获取
otu.percent.temp<-otu.percent 
otu.percent.temp[otu.percent.temp!=0]<-1 
otu.percent.temp<-t(otu.percent.temp) 
result<-data.frame(row.names(otu.percent.temp),
                   rowSums(otu.percent.temp)/ncol(otu.percent.temp),
                   colSums(otu.percent)/nrow(otu.percent))
names(result)<-c('OTU_Id','occurance_frequency','mean_relative_abun') 

#-----------------------------------------------------------------------------------
###S1.按0.1%的标准判断‘overall abundant OTUs’合适吗？
library(ggplot2)
ggplot(result,aes(mean_relative_abun*100, 
                  fill = cut(mean_relative_abun, 50))) + 
  geom_histogram(bins = 50,   
                 show.legend = FALSE)+   
  xlab("Mean relative abundance %")+  
  xlim(0, 1)+  
  theme(axis.title = element_text(size=20,face = 'bold'),     
        axis.text = element_text(size=20))      
ggsave("S1.1 mean_relative_abun的分布.pdf")

#-----------------------------------------------------------------------------------
###S2.按80%的标准判断‘ubiquitous OTUs’合适吗？
frequeny<-as.matrix(result$occurance_frequency) 
frequeny.unique<-unique(frequeny) 
Number_of_OTUs <- NULL 
for(i in 1:nrow(frequeny.unique)){ 
  count=0 
  for(j in 1:nrow(frequeny)){ 
    if(frequeny.unique[i]==frequeny[j]) 
      count= count+1   
  }
  Number_of_OTUs<- rbind(Number_of_OTUs, count)
}
Frequency.distribution<-data.frame(frequeny.unique,Number_of_OTUs)
names(Frequency.distribution)<-c('occurance_frequency','Number_of_OTUs') 
ggplot(Frequency.distribution,aes(x=occurance_frequency,y=Number_of_OTUs))+
  geom_point(size=2) + 
  labs(x="Occurance frequency",              
       y="Number of OTUs")+                  
  theme(axis.title = element_text(size=20,face = 'bold'),     
        axis.text = element_text(size=20))    
ggsave("S2.2出现频率-OTU个数图.pdf")

#-----------------------------------------------------------------------------------
###S3.在>50%的样本中处于top 80%来判断‘frequently abundant OTUs’合理吗？(事实上这是最严格的标准)
otu.percent.temp<-t(otu.percent) 
abundant.otu=NULL 
for (i in 1:ncol(otu.percent.temp)) { 
  otu.percent.temp<-otu.percent.temp[order(otu.percent.temp[,i], decreasing=T),]
  count=0 
  for (j in 1:nrow(otu.percent.temp)){ 
    count=count+otu.percent.temp[j,i] 
    abundant.otu<- rbind(abundant.otu, rownames(otu.percent.temp)[j])
    if(count>=0.8)
      break
  }}
abundant.otu.unique<-unique(abundant.otu)
frequency <- NULL
for(i in 1:nrow(abundant.otu.unique)){ 
  count=0 
  for(j in 1:nrow(abundant.otu)){ 
    if(abundant.otu.unique[i]==abundant.otu[j]) 
      count= count+1   
  }
  frequency<- rbind(frequency, count/nrow(otu.percent))
} 
frequency <- data.frame(frequency) 
library(ggplot2)
ggplot(frequency,aes(frequency)) +  
  geom_histogram(bins = 20, show.legend = FALSE)+  
  xlab("Frequency to be in top 80%")+  
  theme(axis.title = element_text(size=20,face = 'bold'), 
        axis.text = element_text(size=20)) 
ggsave("S3.3 abundant otu进入top 0.8的频率分布.pdf") 
