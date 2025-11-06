# This script describes making plots from STRUCTURE outputs.
# It is heavily based on code written by DrewWham: https://github.com/DrewWham/Genetic-Structure-Tools
# Specifically, the plotSTR.R file can be found here: https://github.com/DrewWham/Genetic-Structure-Tools/blob/6474d60ca34935a7683472a1d943160568a67653/plotSTR.r

# Rversion: 4.3.2

library("dplyr")
source("plotSTR.r")

read.STR<-function(STR.in,STR.out){
  #read in data
  str<-read.table(STR.in,skip=1)	
  str.out<-readLines(STR.out)
  #the next part parses the structure outfile and grabs the q score part
  q.tab.id <- grep("Inferred ancestry of individuals:",str.out,
                   value=FALSE)
  num.ind <- grep("Run parameters:",str.out,
                  value=FALSE)+1
  k.ind<-grep("populations assumed",str.out,
              value=FALSE)
  k<-str.out[k.ind]
  k<-as.numeric(str_extract_all(k,"[0-9.A-Za-z_]+")[[1]][1])
  
  num.ind<-str.out[num.ind]	 				 
  num.ind<-as.numeric(str_extract_all(num.ind,"[0-9.A-Za-z_]+")[[1]][1])
  q.end.id<-1+num.ind+q.tab.id		
  q.tab<-str.out[(q.tab.id+2):q.end.id]
  qs<-t(sapply(str_extract_all(q.tab,"[0-9.A-Za-z_]+"),as.character))
  
  cols_qs<-dim(qs)[2]
  
  qs<-qs[,c(2,(cols_qs-k+1):(cols_qs))]
  qs<-data.frame(qs)
  nclust<- k
  #this part replaces the names with your original names
  qs[,1]<-str[,1]
  write.csv(qs,paste("k",nclust,".csv",sep=""),row.names=F)
}

#### 16S to COI 6kbp amplicon ####
read.STR("Mustelus_16S-COI.txt", "NoAdm_LOC_Ind_1mil_run_2_f")
Amplicon_6kb <- read.csv("k3.csv")
Amplicon_6kb$X2 <- round(Amplicon_6kb$X2)
Amplicon_6kb$X3 <- round(Amplicon_6kb$X3)
Amplicon_6kb$X4 <- round(Amplicon_6kb$X4)
write.csv(Amplicon_6kb, "k3.csv", col.names = T, row.names = F)

Supplementary <- read.csv("Supplementary_data.csv")
Supplementary <- Supplementary[c(1,4)]

#Re-order to custom order
Sorting_order <- c("BoB", "BC", "EC", "IS", "NNS", "SS", "SNS", "WoS")
Supplementary <- Supplementary %>%
  arrange(match(Capture.Region.Code, Sorting_order))

#read data
data<-Amplicon_6kb
plt.order<-Supplementary$Unique_Sample_ID
nclust<-dim(data)[2]-1
#relabel the columns
labs<-c("Specimen",paste("Cluster",seq(1:nclust),sep=""))
og.labs<-paste("X",seq(1:(nclust+1)),sep="")
setnames(data, og.labs, labs)
#reorder, this does not do anything if you dont supply an order file
data <- data[match(Supplementary$Unique_Sample_ID, data$Specimen),]
#I melt here to make the plotting easy
mdata<-melt(data)
names(mdata) <- c("Specimen","Cluster","Probability")
mdata$Specimen <-factor(mdata$Specimen,levels=data$Specimen)
#plot, you can change the palette to something else to explore other color themes
ggplot(mdata,aes(x = Specimen, y = Probability, fill = Cluster)) +
  geom_bar(stat="identity", position="fill") +
  geom_col(width = 0.99) +
  scale_color_viridis_d(option = "D", aesthetics = "fill", alpha = 0.8) +
  theme(axis.text.x = element_text(size = 10, angle = 90, hjust = 1, vjust = -0.00001, colour="black")) +
  theme(axis.title.y = element_blank(), 
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        legend.position="bottom") +
  xlab("Region") +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_discrete(labels = Supplementary$Capture.Region.Code)
#ggsave(paste("k",nclust,"STRplot.pdf",sep=""),width = 20, height = 9)
options(warn = oldw)



 #####   CONTROL REGION ####
setwd("C://Users/berke032/OneDrive - Wageningen University & Research/Documents/Stage Sportvisserij Nederland/STRUCTURE/Mustelus_asterias_20240627_ControlRegion/Manuscript_1mil/Results/")
read.STR("Mustelus_ControlRegion.txt", "Manuscript_1mil_run_2_f")
ControlRegion <- read.csv("../../../Mustelus_asterias_20240627_ControlRegion/Manuscript_1mil/Results/k3.csv")
ControlRegion$X2 <- round(ControlRegion$X2)
ControlRegion$X3 <- round(ControlRegion$X3)
ControlRegion$X4 <- round(ControlRegion$X4)


Sorting_order_CR <- inner_join(Supplementary, ControlRegion, by = c("Unique_Sample_ID" = "X1")) %>%
  arrange(desc(X2), desc(X3), desc(X4)) %>%
  arrange(match(Capture.Region.Code, Sorting_order))

data<-ControlRegion
plt.order<-Sorting_order_CR$Unique_Sample_ID
nclust<-dim(data)[2]-1
#relabel the columns
labs<-c("Specimen",paste("ControlRegionCluster",seq(1:nclust),sep=""))
og.labs<-paste("X",seq(1:(nclust+1)),sep="")
setnames(data, og.labs, labs)
#reorder, this doesnt do anything if you dont supply an order file
data <- data[match(Sorting_order_CR$Unique_Sample_ID, data$Specimen),]
#I melt here to make the plotting easy
mdata<-melt(data)
names(mdata) <- c("Specimen","Cluster","Probability")
mdata$Specimen <-factor(mdata$Specimen,levels=data$Specimen)
#plot, you can change the palette to something else to explore other color themes
threecolours <- pals::cols25(3)
ggplot(mdata,aes(x = Specimen, y = Probability, fill = Cluster)) +
  geom_bar(stat="identity", position="stack") +
  geom_col(width = 0.99) +
  scale_color_viridis_d(option = "D", aesthetics = "fill", alpha = 0.8) +
  #theme_classic() +
  theme(axis.text.x = element_text(size = 10, angle = 90, hjust = 1, vjust = 0.00001, colour="black")) +
  theme(axis.title.y = element_blank(), 
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        legend.position="bottom") +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_discrete(labels = Sorting_order_CR$Capture.Region.Code)
  
#ggsave(paste("k",nclust,"STRplot.pdf",sep=""),width = 20, height = 9)
options(warn = oldw)

 ##### COMPARE CLUSTERS BETWEEEN 16S-COI & ControlRegion ####

ClusterComparison <- inner_join(Amplicon_6kb, ControlRegion)
AA <- sum(subset(ClusterComparison, Cluster1 == 1)$ControlRegionCluster1)
AB <- sum(subset(ClusterComparison, Cluster1 == 1)$ControlRegionCluster2)
AC <- sum(subset(ClusterComparison, Cluster1 == 1)$ControlRegionCluster3)
BA <- sum(subset(ClusterComparison, Cluster2 == 1)$ControlRegionCluster1)
BB <- sum(subset(ClusterComparison, Cluster2 == 1)$ControlRegionCluster2)
BC <- sum(subset(ClusterComparison, Cluster2 == 1)$ControlRegionCluster3)
CA <- sum(subset(ClusterComparison, Cluster3 == 1)$ControlRegionCluster1)
CB <- sum(subset(ClusterComparison, Cluster3 == 1)$ControlRegionCluster2)
CC <- sum(subset(ClusterComparison, Cluster3 == 1)$ControlRegionCluster3)

data.matrix <- matrix(data = c(AA, AB, AC, BA, BB, BC, CA, CB, CC), nrow = 3, ncol = 3, 
                      dimnames = list(c("ControlRegion_Cluster_A", "ControlRegion_Cluster_B", "ControlRegion_Cluster_C"),
                                      c("6kb_Cluster_A", "6kb_Cluster_B", "6kb_Cluster_C")))
mcnemar.test(data.matrix)
library("rcompanion")
install.packages("rcompanion")
