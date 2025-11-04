# This script describes the generation of a haplotype accumulation curve from amplicon data. The assumed input is .fasta
# R version: 4.2.3
# Author: Daniël van Berkel
# Date: 2025-11-04

library("adegenet")
TrimmedSeqs6kbp <- adegenet::fasta2DNAbin("Sequences_6kbp.fasta")
TrimmedSeqs9kbp <- adegenet::fasta2DNAbin("Sequences_9kbp.fasta")
TrimmedSeqsCR <- adegenet::fasta2DNAbin("Sequences_CR.fasta")
TrimmedSeqsCOI <- adegenet::fasta2DNAbin("Sequences_COI.fasta")

library("spider")
HaploAccum_9kbp <- spider::haploAccum(TrimmedSeqs9kbp, permutations = 1000)
HaploAccum_6kbp <- spider::haploAccum(TrimmedSeqs6kbp, permutations = 1000)
HaploAccum_CR <- spider::haploAccum(TrimmedSeqsCR, permutations = 1000)
HaploAccum_COI <- spider::haploAccum(TrimmedSeqsCOI, permutations = 1000)

plot(x = HaploAccum_COI$sequences, y = HaploAccum_COI$n.haplotypes, 
     type = "l", xlim = c(4.8, 125), ylim = c(2.5, 65), xlab = "Number of sequences", ylab = "Number of haplotypes", axes = T)
# Plot each line one by one
lines(HaploAccum_CR$sequences, HaploAccum_CR$n.haplotypes, 
      type = "l", col = "darkcyan", lwd = 3)
lines(HaploAccum_6kbp$sequences, HaploAccum_6kbp$n.haplotypes, 
      type = "l", col = "darkblue", lwd = 3)
lines(HaploAccum_COI$sequences, HaploAccum_COI$n.haplotypes, 
      type = "l", col = "grey60", lwd = 3)

legend(x = 5, y = 60, legend = c("Control Region", "6 kbp fragment", "COI"), 
       col = c("darkcyan", "darkblue", "grey60"), lty = 1, lwd = 2)
