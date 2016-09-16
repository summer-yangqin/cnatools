## prepare data for cnatools package

system("mkdir -p inst/data")

## 1) regions to exclude (T-cell receptors, Immunoglobulin loci).  pad by 100 bases
exclude_loci_hg19 <- data.frame(Name=c("TCRA","TCRB","TCRG","IgH","IgL kappa","IgL lambda"),
                   Chr = c(14,7,7,14,2,22),
                   Start = c(22872077,38296693,142447557,106346666,89165115,23180833),
                   End = c(22942224, 38325607, 142469266, 106580733, 89292554, 23229265),
                   stringsAsFactors=F)
exclude_loci_hg19$Start <- exclude_loci_hg19$Start-100
exclude_loci_hg19$End <- exclude_loci_hg19$End+100
devtools::use_data(exclude_loci_hg19,overwrite=TRUE)

## 2) pseudo-autosomal regions
PAR_coords_hg19 <- data.frame(Chr="X",Start=c(60001,154931044),End=c(2699520,155270560),stringsAsFactors=F)
devtools::use_data(PAR_coords_hg19,overwrite=TRUE)

## 3) sample coordinates for genes

G1 <- read.delim("/mctp/mioncoseq/ref/GRCh37/refseq/genes.refseq",stringsAsFactors=F)
G1 <- subset(G1,name2 %in% c("TP53","APC","ATM","BRCA1","BRCA2","PTEN","EGFR","ERBB2"))
t1 <- gsub("chr","",tapply(G1$chrom,G1$name2,unique))
t2 <- tapply(G1$txStart,G1$name2,min)
t3 <- tapply(G1$txEnd,G1$name2,max)
sample_loci <- data.frame(Locus = names(t1), Chr = t1, Start = NA, End = NA, stringsAsFactors = FALSE)
sample_loci$Start <- t2[sample_loci$Locus]
sample_loci$End <- t3[sample_loci$Locus]
rownames(sample_loci) <- NULL
write.table(sample_loci, file = "inst/data/sample_loci.txt", quote = F, sep = "\t", row.names = F)

## 4) sample coverage and zygosity data files for testing. 1-copy loss of chr3, 1-copy gain of chr12
tc <- 0.45
loss_level <- log2(tc*1 + (1-tc)*2)-log2(2)
gain_level <- log2(tc*3 + (1-tc)*2)-log2(2)
loss_shift <- (tc*1 + (1-tc)*1)/(tc*1 + (1-tc)*2)
gain_shift <- (tc*2 + (1-tc)*1)/(tc*3 + (1-tc)*2)
source("R/coverage_functions.R")    
ChrInfo <- get_chrinfo(c(1:22,"X","Y"),numproc=6)
L <- lapply(c(1:22,"X","Y"), function(chrom) {
    data.frame(Chr = chrom,
               Length = ChrInfo$Length_Orig[which(ChrInfo$ChrName==chrom)],
               Start = 200000*(1:(floor(ChrInfo$Length_Orig[which(ChrInfo$ChrName==chrom)]/200000))),
               End = NA, stringsAsFactors = F)
})
Bedfile <- do.call(rbind,L); Bedfile$End <- Bedfile$Start + 150
Bedfile <- subset(Bedfile, End < Length)
Bedfile$Length <- NULL
Bedfile$Name <- paste0("chr",Bedfile$Chr,":",Bedfile$Start+1,"-",Bedfile$End)

CovData <- Bedfile
CovData$Mean <- rep(0,nrow(CovData))
CovData$Mean[which(CovData$Chr=="3")] <- loss_level
CovData$Mean[which(CovData$Chr=="12")] <- gain_level
set.seed(1000)
CovData$MeanCovN <- runif(nrow(CovData),175,225)
set.seed(2000)
CovData$MeanCovT <- 2^(rnorm(nrow(CovData),CovData$Mean,0.2))*CovData$MeanCovN
CovData$Mean <- NULL
CovData$Length <- NULL
CovData <- CovData[,c("Chr","Start","End","MeanCovT","MeanCovN")]

set.seed(2500)
L <- lapply(c(1:22,"X","Y"), function(chrom) {
    data.frame(Chr = chrom, Pos = sort(sample(1:ChrInfo$Length_Orig[which(ChrInfo$ChrName==chrom)], 1000)),
               Ref = "G", Alt = "A", Var.T = NA, Tot.T = NA, Var.N = NA, Tot.N = NA, stringsAsFactors=F)
})
ZygData <- do.call(rbind,L)

set.seed(3000)
ZygData$Tot.T <- sample(175:225,nrow(ZygData),replace=TRUE)
set.seed(4000)
ZygData$Tot.N <- sample(175:225,nrow(ZygData),replace=TRUE)

ZygData$Mean <- 0.5
ZygData$Mean[which(ZygData$Chr=="3")] <- loss_shift
ZygData$Mean[which(ZygData$Chr=="12")] <- gain_shift

set.seed(5000)
ZygData$Var.T <- rbinom(n = nrow(ZygData),size=ZygData$Tot.T, prob = ZygData$Mean)
set.seed(6000)
ZygData$Var.N <- rbinom(n = nrow(ZygData),size=ZygData$Tot.N, prob = 0.5)

## switch ref and alt for some variants
set.seed(7000)
r1 <- rbinom(n = nrow(ZygData), size = 1, prob = 0.5)
ind1 <- which(r1 == 1)
for (i in ind1) {
    ZygData$Var.T[i] <- ZygData$Tot.T[i]-ZygData$Var.T[i]
}
ZygData$Mean <- NULL

write.table(CovData, file = "inst/data/sample_coverage_data.txt",quote=F,sep="\t",row.names=F)
write.table(ZygData, file = "inst/data/sample_zygosity_data.txt",quote=F,sep="\t",row.names=F)
write.table(Bedfile, file = "inst/data/sample_bedfile.txt",quote=F,sep="\t",row.names=F,col.names=F)
 
