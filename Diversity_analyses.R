##alignment of sequences
# Untrimmed Sequences
library("msa")
AllSequences <- readDNAStringSet("Sequences.fasta")
#AllSequences <- readDNAStringSet("ControlRegion/20240619_Mustelus_asterias_Control_region_trimmed.fasta")
#Order sequences on "aesthetic haplotype" to make table of segregating sites later
AllSequences <- AllSequences[order(names(AllSequences))]
#align all sequences with MUSCLE algorithm
AlignedSeqs <- msa(AllSequences, method = "Muscle", type = "dna", order = "input", verbose = TRUE)
#Convert to DNAbin object
library("ape")
AlignedSeqsDNAbin <- as.DNAbin(AlignedSeqs)
#Trim ends
library("ips")
TrimmedSeqs <- trimEnds(AlignedSeqsDNAbin, min.n.seq = nrow(AlignedSeqs))
# Convert DNAbin object to matrix
seq_matrix <- as.matrix(as.character(TrimmedSeqs))
# Gaps are not always properly processed. It's safer to convert the gaps to an mismatching nucleotide, thereby creating a SNP.
# For instance,  trimEnds function unfortunately converts all gaps with N's, change this back.
# In this case, none of the gaps/Ns match with a 'G'
# Replace "N" with "G"
seq_matrix[seq_matrix == "n"] <- "G"
# Convert the matrix back to DNAbin object
TrimmedSeqs <- as.DNAbin(seq_matrix)

#Basic statistics 
library("pegas")
hap.div(TrimmedSeqs, T)
nuc.div(TrimmedSeqs, T, pairwise.deletion = F)
nSNPs <- length(seg.sites(TrimmedSeqs, T))
nSamples <- nrow(TrimmedSeqs)
theta.s(nSNPs, nSamples, T)
tajima.test(TrimmedSeqs)
strataG::fusFs(TrimmedSeqs)

#Neighbor Joining Matrix + Tree
Tree <- nj(dist.dna(TrimmedSeqs, model = "K80"))

library("ggtree")
library("ggplot2")
GGT <- ggtree(Tree, cex = 0.8, ladderize = T) +
  scale_color_continuous(high = "lightskyblue1",low = "navyblue") +
  geom_tiplab(align = FALSE, size = 1) +
  geom_treescale(x = .003, y = -5, color = "navyblue", fontsize = 6)
GGT

#Haplotype Distribution
#Composition of haplotypes
library("haplotypes")
HaplotypeFreq <- haplotypes::haplotype(as.dna(TrimmedSeqs), indel = "sic")

#Distance Matrix
TrimmedHaps <- pegas::haplotype(TrimmedSeqs, strict = T)
Dist.Mat <- pegas::dist.hamming(TrimmedHaps)
Dist.Mat <-as.matrix(Dist.Mat)
write.table(Dist.Mat, file = "DistanceMatrixHaplotypes.txt", quote=FALSE, sep="\t")

HaploNames <- paste("H", 1:nrow(TrimmedHaps), sep = "")
rownames(TrimmedHaps) = paste(HaploNames)
#Construct network and plot
TrimmedNet <- haploNet(TrimmedHaps)
plot(TrimmedNet, size = attr(TrimmedNet, "freq"), fast = FALSE)

#Plot with groups in each chart
ind.hap <- with(
  utils::stack(setNames(attr(TrimmedHaps, "index"), rownames(TrimmedHaps))),
  table(hap = ind, individuals = rownames(TrimmedSeqs)[values]))

# transform data to fit haplotype network
mydata <- as.data.frame(ind.hap)
mydata2 <- mydata[mydata$Freq==1,]
subpop <- strsplit(as.character(mydata2$individuals), "-")
subpop <- sapply(subpop, "[[", 3)
subpop <- tm::removeNumbers(subpop)
newhap <- table(mydata2$hap, subpop)

#Eight colours for eight areas
eightcolours <- pals::cols25(8)
#Load exact Coordinates to plot
PlotCoordinates <- readRDS("PlotCoordinates.rds")
#State margins
par(mar = c(0.01, 0.01, 0.01, 0.01))
#plot
plot(TrimmedNet, size = attr(TrimmedNet, "freq"), scale.ratio = 4, pie = newhap, labels = T, shape = "circles", asp = 1,
          show.mutation = 1, lwd = 1, lty = 1, bg = eightcolours, fast = F, threshold = c(1, 1))
pegas::replot(xy = PlotCoordinates)
legend(40, 11.5, c("Bay of Biscay", "Bristol Channel", "English Channel","Irish Sea", "Northern North Sea", "Sicily Strait", "Southern North Sea", "West of Ireland"), 
       col=eightcolours, pch = 20, pt.cex = 3.5, cex = 2)
