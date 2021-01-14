####for github@@@#####
pathToCoordinates<- "/Users/singha30/DataAndAnalysis/MutationRateOfExons/PCAWG/ForPSIbussiness/PCAWGsamplesPSI/CleanCoordinates"
ExonFile <- read.table(paste(pathToCoordinates,"GeneWiseMergedExonsFromGreadsForGRCh37.87.bed", sep = "/"),
										sep = "\t", header = T)
										
										
#######make genes using exon coordinates and return non-overlapping gene coordinats##########
writeGenesfile = TRUE
ExonToGenes <- function(ExonsFile) {
	Coordinates <- aggregate(cbind(starts, ends)~ names + seqnames + strands, data=ExonFile, 
						FUN = function(x) c(smallest = min(x), largest= max(x)))
						
	mat1 <- Coordinates[,4]
	mat2 <- Coordinates[,5]
	mat <- cbind(mat1,mat2)		
	
	Minimum <- apply(mat,1,min)
	Maximum <- apply(mat,1,max)
	
	GenesDf <- data.frame("Chr"  = Coordinates$seqnames,
				 	"Start" = Minimum,
					"End" = Maximum,
					"Gene" = Coordinates$names,
					"Strand" = Coordinates$strands)
	
	return(GenesDf)

}


TriNucleotides <- function(x, n) {
	tri <- substring(x, 1:(nchar(x) - (n - 1)), n:nchar(x))
	counts <- table(tri)
	captured <- possibleTri %in% names(counts)
	missing <- possibleTri[!captured]
	if (length(missing) > 0) {
		vec <- rep(0, length(missing))
		names(vec) <- missing
		counts <- c(counts, vec)
	}
	counts <- counts[order(names(counts))]
	#counts <- data.frame(counts[order(names(counts))])
	return(counts)
} 

PossibleMutations <- function(possibleTri) {
	Output <- vector("list", 192)
	index = 0
	bases <- c('A', 'T', 'G', 'C')
	for (i in 1:length(possibleTri)) {
		vec <- strsplit(possibleTri[i], split = "")
		ref <- bases %in% vec[[1]][2]
		alt <- bases[!ref]
		mutation <- NULL
		for (j in 1:3) {
			index = index + 1
			#mutation[j] <- paste(vec[[1]][1],alt[j],vec[[1]][3], sep = "")
			#mutation <- paste(vec[[1]][1],alt[j],vec[[1]][3], sep = "")
			tri <- paste(vec[[1]][1],vec[[1]][2],vec[[1]][3], sep = "")
			mut <- paste(vec[[1]][1],alt[j],vec[[1]][3], sep = "")
			
			Output[index] <- paste(tri,mut, sep = ":")
		}
		#Output[[i]] <- mutation
		#names(Output)[i] <- possibleTri[i]
	}
	return(Output)
}



GenesDf <- ExonToGenes(ExonsFile)
library(bedr)
GenesDf <- GenesDf[order(GenesDf$Chr, GenesDf$Start), ]
GenesDf$Chr <- paste("chr", GenesDf$Chr, sep = "")

chromosomes <- paste("chr", rep(1:24), sep = "")
chromosomes[23] <- "chrX"
chromosomes[24] <- "chrY"
GenesDf <- GenesDf[GenesDf$Chr %in% chromosomes, ]

if (check.binary("bedtools")) {
    MeregedGenes <- bedr(
        engine = "bedtools", 
        input = list(i = GenesDf), 
        method = "merge", 
        params = "-d 1 -c 4 -o collapse"
        );
    }

indexes <- grep(",", MeregedGenes$V4, invert = T)
NonOverlappingGenesDf <- MeregedGenes[indexes, ]
NonOverlappingGenesDf$V2  <- NonOverlappingGenesDf$V2 + 1 ########convert from bed to normal coordinate system
NonOverlappingGenes <- NonOverlappingGenesDf$V4




################ get exons from non-overlapping genes and merge those################
matches <- ExonFile$names %in% NonOverlappingGenes
filteredExonsDf <- ExonFile[matches, ]	
#filteredExonsDf <- filteredExonsDf[order(filteredExonsDf$Chr, filteredExonsDf$Start), ]
filteredExonsDf$seqnames <- paste("chr", filteredExonsDf$seqnames, sep = "")
MeregedExons <- filteredExonsDf
names(MeregedExons) <- c("Chr", "Start", "End", "Gene", "Score", "Strand")
names(NonOverlappingGenesDf) <- c("Chr", "Start", "End", "Gene")

MeregedIntrons <- bedr(
        input = list(a = NonOverlappingGenesDf, b = MeregedExons), 
        method = "subtract", 
        params = ""
        );											


###################Get the sequences#################

library(GenomicRanges)
library(gtools)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19)
genome <- getBSgenome(BSgenome.Hsapiens.UCSC.hg19)

MeregedExons$Length <- MeregedExons$End - MeregedExons$Start
MeregedIntrons$Length <- MeregedIntrons$End - MeregedIntrons$Start

MeregedExons <- MeregedExons[MeregedExons$Length >=10, ]
MeregedIntrons <- MeregedIntrons[MeregedIntrons$Length >=10, ]

matches <- MeregedExons$Gene %in% MeregedIntrons$Gene  #########Filter Out genes with introns###########
MeregedExons <- MeregedExons[matches, ] #########Filter Out genes with introns###########

matches <- MeregedIntrons$Gene %in% MeregedExons$Gene  #########Filter Out genes with introns###########
MeregedIntrons <- MeregedIntrons[matches, ] #########Filter Out genes with introns###########



MeregedExons$Start <- MeregedExons$Start + 1  #######changeCoordinate from bed to normal
MeregedIntrons$Start <- MeregedIntrons$Start + 1 ###########changeCoordinate from bed to normal

writeFiles <- TRUE

if (writeFiles) {
	temp <- MeregedExons
	temp$Start <- temp$Start - 1
	write.table (temp, "ExonsForReducedMutationBurdenAnalysis.bed", sep = "\t", quote = F, row.names = F)
	temp <- MeregedIntrons
	temp$Start <- temp$Start - 1
	write.table (temp, "IntronsForReducedMutationBurdenAnalysis.bed", sep = "\t", quote = F, row.names = F)
}

exonRanges <- makeGRangesFromDataFrame(MeregedExons, keep.extra.columns = TRUE)
intronRanges <- makeGRangesFromDataFrame(MeregedIntrons, keep.extra.columns = TRUE)

exonSeq <- getSeq(genome, exonRanges)
exonSeq  <- as.character(exonSeq)

intronSeq <- getSeq(genome, intronRanges)
intronSeq  <- as.character(intronSeq)

possibleTri <- permutations(4, 3, c("A", "T", "G", "C"), repeats.allowed = T)
possibleTri <- apply(possibleTri, 1, paste, collapse="")


######get trinucleotides#########
ObservedExontri <- lapply(exonSeq, TriNucleotides, 3)
temp <- grep("N", exonSeq)
ObservedExontri[temp] <- NULL
ObservedExontri <- matrix(unlist(ObservedExontri),ncol = 64, byrow = T)
colnames(ObservedExontri) <- possibleTri
rownames(ObservedExontri) <- names(exonSeq)[-temp]
ObservedExontriDf <- as.data.frame(ObservedExontri)
rm(ObservedExontri)
Genes <- as.character(MeregedExons$Gene)
if (length(temp > 0)) {
	Genes <- as.character(Genes[-temp])
}
ObservedExontriDf$Gene <- Genes

ObservedIntrontri <- lapply(intronSeq, TriNucleotides, 3)
temp <- grep("N", intronSeq)
ObservedIntrontri[temp] <- NULL
ObservedIntrontri <- matrix(unlist(ObservedIntrontri),ncol = 64, byrow = T)
colnames(ObservedIntrontri) <- possibleTri
rownames(ObservedIntrontri) <- names(intronSeq)[-temp]
ObservedIntrontriDf <- as.data.frame(ObservedIntrontri)
rm(ObservedIntrontri)
Genes <- as.character(MeregedIntrons$Gene)
if (length(temp > 0)) {
	Genes <- as.character(Genes[-temp])
}
ObservedIntrontriDf$Gene <- Genes


#PossibleMutations <- function(possibleTri) {
#	Output <- vector("list", 64)
#	bases <- c('A', 'T', 'G', 'C')
#	for (i in 1:length(possibleTri)) {
#		vec <- strsplit(possibleTri[i], split = "")
#		ref <- bases %in% vec[[1]][2]
#		alt <- bases[!ref]
#		mutation <- NULL
#		for (j in 1:3) {
#			mutation[j] <- paste(vec[[1]][1],alt[j],vec[[1]][3], sep = "")
#		}
#		Output[[i]] <- mutation
#		names(Output)[i] <- possibleTri[i]
#	}
#	return(Output)
#}
#Output <- PossibleMutations(possibleTri)


######get possible mutations######

Output <- PossibleMutations(possibleTri)
Output <- unlist(Output)




#CancerTypes <- c('SKCM', 'PRAD', 'COAD', 'BRCA')

CancerTypes <- c("Ovary-AdenoCA", "Liver-HCC", 
				"Panc-Endocrine", "Kidney-RCC", "Prost-AdenoCA", "Lymph-BNHL", 
				"Panc-AdenoCA", "Eso-AdenoCa", "CNS-Medullo", "Lymph-CLL", "Skin-Melanoma", 
				"Stomach-AdenoCA", "Breast-AdenoCa", "Head-SCC", "Lymph-NOS", 
				"Biliary-AdenoCA", "Bone-Osteosarc", "Breast-LobularCa", "Myeloid-MPN", "Bone-Epith")
				
library(stringr)


mutationPath <- "/Users/singha30/DataAndAnalysis/MutationRateOfExons/PCAWG/CancerWiseMutations"

OutputDf <- NULL
library(plyr)
for (mf in 1:length(CancerTypes)) {
	
	paste("Doing For", CancerTypes[mf])
	MutationFile <- read.table(paste(mutationPath, "/PCAWGmutationsIn", CancerTypes[mf], ".txt", sep =""),
							sep = "\t")

	names(MutationFile) <- c("Chr", "Start", "End", "Ref", "Alt", "Name", "Sample")
	MutationFile <- na.omit(MutationFile)	#####because some cosmic muttions had missing coordinats#####
	MutationFile$Chr <- str_replace_all(MutationFile$Chr, "chr23", "chrX")
	MutationFile$Chr <- str_replace_all(MutationFile$Chr, "chr24", "chrY")
	MutationFile$End <- MutationFile$End + 1


	mRanges <- makeGRangesFromDataFrame(MutationFile, keep.extra.columns = TRUE)
	sequences <- getSeq(genome, mRanges)
	sequences <- as.character(sequences)

	sigDf <- array(NA,c(length(sequences), 2))	
		for (i in 1:length(sequences)) {
			f1 <- substr(sequences[i],1,1)
			ref <- substr(sequences[i],2,2)
			f3 <- substr(sequences[i],3,3)
			if(ref == MutationFile[i,4]) {
				alt <- MutationFile[i,5]
			} else {
				alt <- chartr("ATGC","TACG",MutationFile[i,5])
			}
			sigDf[i,1] <- paste(f1,ref,f3, sep ="")
			sigDf[i,2] <- paste(f1,alt,f3, sep ="")
		}
	
	colnames(sigDf) <- c("Ref", "Alt")
	sigDf <- as.data.frame(sigDf)
	
	backGround <- count(sigDf,vars = c("Ref","Alt"))
	backGround <- backGround[backGround$Ref != backGround$Alt, ]
	
	###### if a particular signature did not exist, add it with zero frequecny######
	observedMutations <- paste(backGround$Ref, backGround$Alt, sep = ":")   
	matches <- Output %in% observedMutations   #####Output contains all possible 192 signatures
	missing <- Output[!matches]
	if (length(missing) > 0) {
		MissingBackground <- NULL
		for (k in 1:length(missing)) {
			Ref <- strsplit(missing, ":")[[k]][1]
			Alt <- strsplit(missing, ":")[[k]][2]
			freq <- 0
			missingLine <- data.frame(Ref, Alt, freq)
			MissingBackground <- rbind(MissingBackground, missingLine)
		}
		backGround <- rbind(backGround, MissingBackground)
	}
	
	backGround <- backGround[order(backGround$Ref), ]
	
	GenomeTri <- read.table("trinucleotides.txt", sep = "\t")
	GenomeTri <- GenomeTri[match(backGround$Ref, GenomeTri$V1),2]
	backGround$Prob <- backGround$freq/GenomeTri
	sumProb <- sum(backGround$Prob)
	backGround$Prob <- backGround$Prob/sumProb

	
	#backGround <- read.table("MutationalProbablityFileForTest.txt", sep = "\t", header = T)
	signaturesSum <- aggregate(Prob ~ Ref, data = backGround, sum)
	
	ProductSums <- NULL
	for (i in 1:dim(ObservedExontriDf)[1]) {
		#p = 0
		#for (context in 1:64) {
			#tri <- names(ObservedExontriDf)[context]
			#signatures <- backGround[backGround$Ref == tri, ]
			#tri <- ObservedExontriDf[1, ]
			#signaturesSum <- aggregate(Prob ~ Ref, data = backGround, sum)		
			#p1 <- ObservedExontriDf[i, context] * signatures[1,'Prob']
			#p2 <- ObservedExontriDf[i, context] * signatures[2,'Prob']
			#p3 <- ObservedExontriDf[i, context] * signatures[3,'Prob']
			#p <- p + p1 + p2 + p3
			#p = p + sum(ObservedExontriDf[i, context] * signatures[,'Prob'])
			#}
		#ProductSums[i] <- p
		ProductSums[i] <- sum(ObservedExontriDf[i, 1:64] * signaturesSum[,'Prob'])
	}

	ExonProductSumDf <- data.frame("Gene" = ObservedExontriDf$Gene, 
						"Exon" = rownames(ObservedExontriDf), 
						"ProductSum" <- ProductSums)

					
	rnames <- paste(MeregedExons$Chr, ":", MeregedExons$Start, "-", MeregedExons$End, sep = "")
	ExonProductSumDf$Exon <- rnames
	names(ExonProductSumDf)[3] <- "ProductSums"

	ProductSums <- NULL	
	for (i in 1:dim(ObservedIntrontriDf)[1]) {
		#p = 0
		#for (context in 1:64) {
			#tri <- names(ObservedIntrontriDf)[context]
			#signatures <- backGround[backGround$Ref == tri, ]			
			#p1 <- ObservedIntrontriDf[i, context] * signatures[1,'Prob']
			#p2 <- ObservedIntrontriDf[i, context] * signatures[2,'Prob']
			#p3 <- ObservedIntrontriDf[i, context] * signatures[3,'Prob']
			#p <- p + p1 + p2 + p3
			#p = p + sum(ObservedExontriDf[i, context] * signatures[,'Prob'])
			#}
		ProductSums[i] <- sum(ObservedIntrontriDf[i, 1:64] * signaturesSum[,'Prob'])
	}

	IntronProductSumDf <- data.frame("Gene" = ObservedIntrontriDf$Gene, 
						"Intron" = rownames(ObservedIntrontriDf), 
						"ProductSum" <- ProductSums)
					
	names(IntronProductSumDf)[3] <- "ProductSums"



	###############Map mutations to exons and introns################
	MappedExonicMutations <- bedr(
	        input = list(a = MeregedExons, b = MutationFile), 
	        method = "intersect", 
	        params = "-c"
	        );

	MappedExonicMutations$V8 <- as.numeric(as.character(MappedExonicMutations$V8))

	MappedIntronicMutations <- bedr(
	        input = list(a = MeregedIntrons, b = MutationFile), 
	        method = "intersect", 
	        params = "-c"
	        );

	MappedIntronicMutations$V6 <- as.numeric(as.character(MappedIntronicMutations$V6))


	ExonicMutationsDf <- aggregate(V8 ~ V4, data = MappedExonicMutations, sum)
	IntronicMutationsDf <- aggregate(V6 ~ V4, data = MappedIntronicMutations, sum)
	matches <- ExonicMutationsDf$V4 %in% IntronicMutationsDf$V4
	ExonicMutationsDf <- ExonicMutationsDf[matches, ]

	IntronProductSumDf1 <- aggregate(ProductSums ~ Gene, data = IntronProductSumDf, sum)
	ExonProductSumDf1 <- aggregate(ProductSums ~ Gene, data = ExonProductSumDf, sum)
	matches <- ExonProductSumDf1$Gene %in% IntronProductSumDf1$Gene
	ExonProductSumDf1 <- ExonProductSumDf1[matches, ]
	#########Expected Mutations way 1: each exon different expections############

	GGdata <- data.frame("Gene" = ExonicMutationsDf$V4,
						"ExonicProductSum" = ExonProductSumDf1$ProductSums,
						"IntronicProductSum" = IntronProductSumDf1$ProductSums,
						"ExonicMutations" = ExonicMutationsDf$V8,
						"IntronicMutations" = IntronicMutationsDf$V6)

	GGdata$TotalMutations <- GGdata$ExonicMutations + GGdata$IntronicMutations										
	GGdata$ExonicProbablity <- (GGdata$ExonicProductSum)/(GGdata$ExonicProductSum + GGdata$IntronicProductSum)
	GGdata$IntronicProbablity <- (GGdata$IntronicProductSum)/(GGdata$ExonicProductSum + GGdata$IntronicProductSum)
	GGdata$ExpectedExonic <- GGdata$TotalMutations * GGdata$ExonicProbablity
	GGdata$ExpectedIntronic <- GGdata$TotalMutations * GGdata$IntronicProbablity

	GGsub <- GGdata[GGdata$TotalMutations > 0, ]
	
	TEM <- sum(GGsub$ExonicMutations)
	TIM <- sum(GGsub$IntronicMutations)
	
	EEM <- sum(GGsub$ExpectedExonic)
	EIM <- sum(GGsub$ExpectedIntronic)
	
	lineDf <- data.frame(Cancer = CancerTypes[mf], ObservedExonic = TEM, ObservedIntronic = TIM, ExpectedExonic = EEM, ExpectedIntronic = EIM)
	OutputDf <- rbind(OutputDf, lineDf)
	write.table(lineDf, "OEmutationsRunningFile.txt", sep = "\t", row.names = F, quote = F, append = T)
	paste("Done For", CancerTypes[mf])
	Sys.sleep(5)
}

write.table(OutputDf, "ObservedAndExpectedMutationsAcrossCancerTypes.txt", sep = "\t", row.names = F, quote = F)

##sOutputDf
