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
colnames(Amplicon_6kb) <- c("Unique_Sample_ID", "Cluster_6kb_A", "Cluster_6kb_B", "Cluster_6kb_C")
Amplicon_6kb <- left_join(Amplicon_6kb, Supplementary)
Sorting_order <- c("BoB", "BC", "EC", "IS", "NNS", "SS", "SNS", "WoS")
Amplicon_6kb <- Amplicon_6kb %>%
  arrange(Cluster_6kb_C, Cluster_6kb_B, Cluster_6kb_A) %>%
  arrange(match(Capture.Region.Code, Sorting_order)) 

Amplicon_6kb$Unique_Sample_ID <- factor(Amplicon_6kb$Unique_Sample_ID, levels = Amplicon_6kb$Unique_Sample_ID)
#read data
plt.order <- Amplicon_6kb$Unique_Sample_ID
#I melt here to make the plotting easy
mdata <- reshape2::melt(Amplicon_6kb)
names(mdata) <- c("Specimen", "Capture.Region.Code", "Cluster", "Probability")
#plot, you can change the palette to something else to explore other color themes
ggplot(mdata, aes(x = Specimen, y = Probability, fill = Cluster)) +
  geom_bar(stat="identity", position="fill") +
  geom_col(width = 1) +
  theme_classic(base_size = 50) +
  theme(axis.text.x = element_text(size = 10, angle = 90, hjust = 1, vjust = -0.00001, colour="black")) +
  theme(axis.title.y = element_blank(), 
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        legend.position = "bottom") +
  xlab("Region") +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_discrete(labels = Amplicon_6kb$Capture.Region.Code) +
  scale_fill_manual(labels = c("Cluster A", "Cluster B", "Cluster C"), values = viridis::viridis(3))
#ggsave(paste("k",nclust,"STRplot.pdf",sep=""),width = 20, height = 9)
options(warn = oldw)



 #####   CONTROL REGION ####
read.STR("Mustelus_ControlRegion.txt", "Manuscript_1mil_run_2_f")
ControlRegion <- read.csv("k3.csv")
ControlRegion$X2 <- round(ControlRegion$X2)
ControlRegion$X3 <- round(ControlRegion$X3)
ControlRegion$X4 <- round(ControlRegion$X4)
# The clusters names are arbitrary. In this case, label them according to the 6kbp fragment, so the colours of the plot are similar.
colnames(ControlRegion) <- c("Unique_Sample_ID", "Cluster_CR_C", "Cluster_CR_B", "Cluster_CR_A")
ControlRegion <- ControlRegion %>% select(Unique_Sample_ID, Cluster_CR_A, Cluster_CR_B, Cluster_CR_C)
ControlRegion <- left_join(ControlRegion, Supplementary)
ControlRegion <- ControlRegion %>%  arrange(Cluster_CR_C, Cluster_CR_B, Cluster_CR_A) %>%
  arrange(match(Capture.Region.Code, Sorting_order)) 

ControlRegion$Unique_Sample_ID <- factor(ControlRegion$Unique_Sample_ID, levels = ControlRegion$Unique_Sample_ID)
#I melt here to make the plotting easy
mdata.CR <- reshape2::melt(ControlRegion)
names(mdata.CR) <- c("Specimen", "Region", "Cluster", "Probability")
#plot, you can change the palette to something else to explore other color themes
ggplot(mdata.CR,aes(x = Specimen, y = Probability, fill = Cluster)) +
  geom_bar(stat="identity", position="stack") +
  geom_col(width = 1) +
  theme_classic(base_size = 50) +
  theme(axis.text.x = element_text(size = 10, angle = 90, hjust = 1, vjust = 0.00001, colour="black")) +
  theme(axis.title.y = element_blank(), 
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        legend.position="bottom") +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_discrete(labels = ControlRegion$Unique_Sample_ID) +
  scale_fill_manual(labels = c("Cluster A", "Cluster B", "Cluster C"), values = viridis::viridis(3))
  
 ##### COMPARE CLUSTERS BETWEEEN 16S-COI & ControlRegion ####
# As mentioned, the cluster names are arbitrary. Rename the Control Region cluster so they overlap with the 6kbp amplicon.
ClusterComparison <- inner_join(Amplicon_6kb, ControlRegion)
colnames(ClusterComparison) <- c("Specimen", "Cluster_6kb_A", "Cluster_6kb_B", "Cluster_6kb_C", "Capture.Region.Code", "Cluster_CR_A", "Cluster_CR_B", "Cluster_CR_C")

AA <- sum(subset(ClusterComparison, Cluster_6kb_A == 1)$Cluster_CR_A)
AB <- sum(subset(ClusterComparison, Cluster_6kb_A == 1)$Cluster_CR_B)
AC <- sum(subset(ClusterComparison, Cluster_6kb_A == 1)$Cluster_CR_C)
BA <- sum(subset(ClusterComparison, Cluster_6kb_B == 1)$Cluster_CR_A)
BB <- sum(subset(ClusterComparison, Cluster_6kb_B == 1)$Cluster_CR_B)
BC <- sum(subset(ClusterComparison, Cluster_6kb_B == 1)$Cluster_CR_C)
CA <- sum(subset(ClusterComparison, Cluster_6kb_C == 1)$Cluster_CR_A)
CB <- sum(subset(ClusterComparison, Cluster_6kb_C == 1)$Cluster_CR_B)
CC <- sum(subset(ClusterComparison, Cluster_6kb_C == 1)$Cluster_CR_C)

data.matrix <- matrix(data = c(AA, AB, AC, BA, BB, BC, CA, CB, CC), nrow = 3, ncol = 3, 
                      dimnames = list(c("Cluster_CR_A", "Cluster_CR_B", "Cluster_CR_C"),
                                      c("Cluster_6kb_A", "Cluster_6kb_B", "Cluster_6kb_C")))
data.matrix
fisher.test(data.matrix, hybrid = T)
