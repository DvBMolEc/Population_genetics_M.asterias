# This script obtains segregating sites in amplicon data from several individuals, stored in a multifasta.
# R version: 4.2.3
# Author: Daniël van Berkel
  # Adopted from: ?
# Date: 2025-11-04

##alignment of sequences
# Untrimmed Sequences ----------------------------------------------------
library("msa")
AllSequences <- readDNAStringSet("Sequences.fasta")

library("DECIPHER")
#Orient .fastas in the same way to ease alignment
AllSequences <- OrientNucleotides(AllSequences)

#align all .fastas with MUSCLE algorithm
AlignedSeqs <- msa(AllSequences, method = "Muscle", type = "dna", order = "input", verbose = TRUE)

library("haplotypes")
AlignedSeqsDNAbin <- as.DNAbin(AlignedSeqs)

#Trim ends
library("ips")
TrimmedSeqs <- trimEnds(AlignedSeqsDNAbin, min.n.seq = nrow(AlignedSeqs))

library("ape")
Alignment <- as.alignment(TrimmedSeqs)
Matrix <- as.matrix.DNAbin(TrimmedSeqs)

#########   EXTRACTION SEQUENCE AND HAPLOTYPE INFORMATION    ###############

nrow(Matrix)#confirmation of number of samples
ncol(Matrix)#confirmation of sequences size

sat2 <- NULL
for (i in 1:nrow(Matrix)) {
  sat2[i] <- paste(Matrix[i, ], collapse="")
}
sat2 <- toupper(sat2) #converts all letters to uppercase
sat3 <- unique(sat2) #gives only unique sequences from all sequences
sat3 #that is, it gives complete sequences of haplotypes
hfreq <- NULL
for (i in 1:length(sat3)) {
  hcount = 0
  s3 <- sat3[i]
  for (j in 1:length(sat2)) {
    s2 <- sat2[j]
    if (s3 == s2) {
      hcount <- (hcount + 1) #counts the number of individuals with the same haplotype sequence. 
      #print(paste(i, "yes", hcount))
    }
    #print(s2)
  }
  hname<-(paste("H",i, sep =""))
  hfreq[i] <- hcount
  #print(paste(hname, hcount, collapse = ""))
}   #haplotype frequency in the all samples

len <- nchar(sat3[1]) #assume all have same length!!!
cnt <- 1
sat4 = list()
for (j in 1:len) {
  same <- TRUE
  first <- substr(sat3[1], j, j)
  for (i in 2:length(sat3)) {
    ch1 <- substr(sat3[i], j, j)
    if (first != ch1) {
      str <- paste(j, first, ch1)
      print(str)
      same <- FALSE
      break
    }
  }
  if (!same) {
    ss <- NULL
    for (i in 1:length(sat3)) {
      ss <- paste(ss, substr(sat3[i], j, j), sep="")
    }
    sat4[cnt] <- ss
    cnt <- cnt + 1
  }
}#it gives the mutation points and the nucleotide substitutions

len <- nchar(sat3[1]) #assume all have same length
cnt <- 1
sat5 = list() 
for (j in 1:len) { #scan all columnns and if all elements are the same do not copy
  same <- TRUE
  first <- substr(sat3[1], j, j)
  scol <- first
  for (i in 2:length(sat3)) {
    ch1 <- substr(sat3[i], j, j)
    scol <- paste(scol, ch1, sep="")
    if (first != ch1) {
      str <- paste(j, first, ch1)
      #print(str)
      same <- FALSE
      #break
    }
  }
  if (!same) {
    scol <- paste("V_", cnt, " ", scol, sep="")
    ss <- NULL
    for (i in 1:length(sat3)) {
      ss <- paste(ss, substr(sat3[i], j, j), sep="")
    } 
    sat5[cnt] <- ss
    cnt <- cnt + 1
  }
}

sat6 <- as.matrix(sat5)
mat6 = matrix(nrow=nrow(sat6), ncol=nchar(sat6[1]))
for (i in 1:nrow(mat6)) {
  s <- as.vector(strsplit(as.character(sat5[i]), ""))
  for (j in 1:ncol(mat6)) {
    mat6[i, j] <- as.character(s[[1]][j])
  }
}
mat7 <- t(mat6) #sequences of haplotypes and variable sites matrix
write.table(mat7,file="Variable_Sites_Mat2024.txt", quote=FALSE, sep="\t")
hname<-paste("H", 1:nrow(mat7), sep = "")
rownames(mat7)=hname
write.table(mat7,file="Variable_Sites_Mat2024.txt", quote=FALSE, sep="\t")

str4 <- NULL
str4[1] <- paste(mat7[1, ], collapse="")
for (i in 2:nrow(mat7)) {
  tmp <- NULL
  for (j in 1:ncol(mat7)) {
    chr = "."
    if(mat7[i, j] != mat7[1, j]) chr = mat7[i, j]
    tmp <- paste(tmp, chr, sep="")
  }
  str4[i] <- paste(tmp, collapse="")
}
nchar(str4[1]) #confirmation of number of variable sites
mstr4<-as.matrix(str4)
rownames(mstr4)<-hname
colnames(mstr4)<-paste("sequences length","(", ncol(mat7), "base pairs", ")")
pct<-round((as.matrix(hfreq)*100/colSums(as.matrix(hfreq))), 2)
colnames(pct)<-c("pct")
cmstr4<-as.data.frame(cbind(mstr4, hfreq, pct))
cmstr4

write.table(cmstr4,file="Variable_Sites_Matrix.txt", quote=FALSE, sep="\t") 
