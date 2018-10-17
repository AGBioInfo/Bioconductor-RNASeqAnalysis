# What fraction of reads in this file has an A nucleotide in the 5th base of the read?
library(ShortRead) 
library(yeastRNASeq)
fastqFilePath <-  system.file("reads", "wt_1_f.fastq.gz", package = "yeastRNASeq")
fqFile <- FastqFile(fastqFilePath)
reads <- readFastq(fqFile)
reads_set <- sread(reads)
sum(DNAStringSet(reads_set,5,5) == "A") / length(reads_set)

#What is the average numeric quality value of the 5th base of these reads?
qm <- as(quality(reads), "matrix")
mean(qm[,5:5])

#In this interval, how many reads are duplicated by position?
library(leeBamViews)
bamFilePath <- system.file("bam", "isowt5_13e.bam", package="leeBamViews")
bamFile <- BamFile(bamFilePath)
seqinfo(bamFile)
aln <- scanBam(bamFile)
aln <- aln[[1]]  
names(aln)
lapply(aln, function(xx) xx[1])
unique(aln$rname)
gr <- GRanges(seqnames = "Scchr13", ranges = IRanges(start = 800000, end = 801000))
params <- ScanBamParam(which = gr, what = scanBamWhat())
aln <- scanBam(bamFile, param = params)
aln <- aln[[1]]  
aln$pos 
duplicatedValues = unique(aln$pos[duplicated(aln$pos)]) 
sum(aln$pos %in% duplicatedValues)


# What is the absolute value of the log foldchange ( logFC) of the gene with the lowest P.value.
library(limma)
design <- model.matrix(~ normData$group)
fit <- lmFit(normData, design)
fit <- eBayes(fit)
topTable(fit)
abs(topTable(fit, n=1)$logFC) 

#How many genes are differentially expressed between the two groups at an adj.P.value cutoff of 0.05?

topTable(fit, p.value = 0.05)

#What is the mean difference in beta values between the 3 normal samples and the 3 cancer samples, across OpenSea CpGs?
library(minfi)
require(minfiData)
data(RGsetEx)
p <- preprocessFunnorm(RGsetEx)
b <- getBeta(p)
is <- getIslandStatus(p)
pData(p)$status
norm <- b[,c(1,2,5)]
can <- b[,c(3,4,6)]
norm_os <- norm[is == "OpenSea",]
can_os <- can[is == "OpenSea",]
mean(norm_os) - mean(can_os)

#How many of these DNase hypersensitive sites contain one or more CpGs on the 450k array?

library(AnnotationHub)
ah <- AnnotationHub()
qah_h1 <- query(ah, c("Caco2", "AWG"))
h <- qah_h1[["AH22442"]]
sum(countOverlaps(p,h))  
ah_s <- subset(ah, genome == "hg19")
ah_s <- subset(ah, dataprovider == "UCSC")
# write.csv(mcols(ah), "ah.csv")
# g <- ah[["AH5018"]] # assembly
cpg <- ah_s[["AH5086"]] # CpG islands
h_cpg <- subsetByOverlaps(cpg, h)
ov <- subsetByOverlaps(h_cpg, p)

#How many features are differentially expressed between control and treatment (ie. padj <= 0.05)?
library(DESeq2)
library(zebrafishRNASeq)
data(zfGenes)
#exclude spike-in controls
tx <- zfGenes[grep("^ERCC", rownames(zfGenes), invert = T),]
counts_mat <- as.matrix(tx)
colData <- DataFrame(sampleID=colnames(tx), group=as.factor(c("control", "control", "control", "treatment", "treatment", "treatment")))
ddsMat <- DESeqDataSetFromMatrix(counts_mat, colData, design = ~ group)
ddsMat <- DESeq(ddsMat)
res <- results(ddsMat)
res <- res[order(res$padj),]
sigRes <- subset(res, padj <= 0.05)
dim(sigRes)


