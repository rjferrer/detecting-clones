rm(list = ls())

library(data.table)

#load SLIM file to R as data table; skip lines until "CHROM"
#separator is tab
slimdata <- fread("VCFfile PATH", skip = "CHROM", sep = "\t")
str(slimdata) #row 4 = ref, 5 = alt, 10 = first individual
dim(slimdata) #row = 14649; column = 6009

##--------------------------------------------------------------------------SAMPLING THE POPULATIONS-------------------------------------------------------------------
#[1] Subset the Population into 3-----------------------------------------------
pop1 <- slimdata[ , 10:2009]
pop2 <- slimdata[ , 2010:4009]
pop3 <- slimdata[ , 4010:6009]

#[2] Function for sampling individuals------------------------------------------
popFunc <- function(subPop){
  tempPop <- as.list(subPop)
  dupePop <- which(duplicated(tempPop))
  dupeInd <- c(which(duplicated(tempPop, fromLast=TRUE)),dupePop)
  count <- 1:2000
  set.seed(99)
  pick.dupePop <- sample(dupePop, size = 5, replace = FALSE)
  pick.uniPop <- sample(count[!count %in% (c(dupePop,dupeInd))], size = 20, replace = FALSE)
  seqPop <- as.data.frame(c(tempPop[pick.dupePop],tempPop[pick.dupePop],tempPop[pick.uniPop]))
  return(seqPop)
}

#[3] Sample individuals---------------------------------------------------------
subpopData <- cbind(popFunc(pop1), popFunc(pop2), popFunc(pop3))
colnames(subpopData) <- paste(rep(1:3, each = 30),rep("i", each=90),c(1:30,1:30,1:30), sep = "_", collapse = NULL)
dim(subpopData) #14649    90

##------------------------------------------------------------------------CONVERTING VCF TO SEQUENCES------------------------------------------------------------------
#[1] splitting -----------------------------------------------------------------
splitData <- splitstackshape::cSplit(subpopData, 1:length(subpopData),sep = "|", direction = "wide")
dim(splitData)  ##14649   180

splitAlt <- splitstackshape::cSplit(slimdata, 5, sep = ",", direction = "wide")

#check:
alt2Ind <- grep(pattern = ",", x = slimdata$ALT)
checkInd <- which(!is.na(splitAlt$ALT_2))
checkInd
alt2Ind

dim(splitAlt) #14649 6010
snpData <- cbind(slimdata[,1:4], splitAlt$ALT_1, splitAlt$ALT_2, splitData)
colnames(snpData) <- c(colnames(slimdata[,1:4]), "ALT1", "ALT2", colnames(splitData))
View(snpData)

#[2] substitute binary genotype with A, T, C, G --------------------------------
subData <- data.frame(
  apply(snpData, MARGIN = 1, function(x){
    x <- as.character(x) #converts splitdata to list with each column as list element
    ref = x[4] 
    alt1 = x[5]
    alt2 = x[6]
    tempOut <- ifelse(x[7:length(x)]=="0",ref,ifelse(x[7:length(x)]=="1",alt1,alt2))
    return(tempOut)
  }))
dim(subData) #180 14649
subData <-t(subData)
dim(subData) #14649   180
View(subData[1:10, 1:10])

saveRDS(subData, file="PATH/Genotypes.RDS")

#[3] position information-------------------------------------------------------
#locus identifier
locFloor = floor(slimdata$POS/500)
length(locFloor) #14649
length(unique(locFloor)) #4761
#match returns a vector of the positions of (first) matches of its first argument in its second
locSeq = match(locFloor,unique(locFloor)) 
length(locSeq) #14649
tail(sort(locSeq)) #values are 1-4761

#specify base position in locus
locPos = slimdata$POS - (locFloor*500)
length(locPos) #14649
sort(unique(locPos)) #0-499
#replace "position 0" in locPos with 500
locPos <- ifelse(locPos == "0", 500, locPos)
sort(unique(locPos))  #1-500

#[4] generate template sequences------------------------------------------------
#generate bases to fill in positions not in slim output
tempHap <- data.frame(do.call(cbind, replicate(500, "C", simplify = FALSE)))
dim(tempHap) #1 500
#replicate tempHap for each unique locus position in the data set 
#simplify = F outputs list
outSeq <- data.frame(do.call(rbind, replicate(length(unique(locFloor)), tempHap, simplify = FALSE)))
dim(outSeq) #4761  500
#outseq is all loci seqs (4761 loci 500 bp long each) for 1 haplotype
#allseq is all loci seqs for all haplotypes (180)
allSeq <- lapply(1:ncol(subData), function(x) outSeq)
length(allSeq) #180
#check:
dim(allSeq[[180]]) #4761  500

#[5] input genotypes on template sequences--------------------------------------
subData <- readRDS("PATH/Genotypes.RDS")

seqOut <- lapply(1:length(allSeq), function(x,posIn, seqIn, genoIn, tempIn){
  for (i in 1:nrow(genoIn)) {
    tempIn[[x]][seqIn[i],posIn[i]] <- genoIn[i,x]
    message(x)
  } 
  return(tempIn[[x]])
},
posIn = locPos,
seqIn = locSeq,
genoIn = subData,
tempIn = allSeq
)

##--------------------------------------------------------------------------PCR INTRODUCTION OF ERRORS-------------------------------------------------------------------
#[1] Creating function for pcr errors-------------------------------------------
pcrFunc <- function(inSeq){
  numCyc <-  5
  pcrOut <- vector("list", length = numCyc)
  for (cycle in 1:numCyc) {
    if(cycle == 1){
      pcrOut[1] <- inSeq
    }else{
      repSeq <- strsplit(pcrOut[[cycle-1]][1], "")
      set.seed(100+cycle)
      posChange <- which(rbinom(nchar(inSeq), 1, 0.1) == 1)
      repSeq[[1]][posChange] <-  unlist(lapply(repSeq[[1]][posChange], function(x){
        baseList <- c("A", "T", "G", "C")
        return(sample(baseList[which(!(baseList %in% x))], 1))
      }))
      pcrOut[cycle] <- paste(repSeq[[1]], collapse = "")
    }
  }
  return(pcrOut)
}

#[2] Applying pcr function------------------------------------------------------
seqOut2 <- unlist(lapply(1:length(seqOut), function(x) apply(seqOut[[x]], 1, paste, collapse = "")))
pcrOut <- lapply(seqOut2[1:length(seqOut2)], function(x){
  pcrFunc(x)
})

length(pcrOut) #4761*180 = 856980
pcrOut[[1]]# 5 pcr replicates <3

##---------------------------------------------------------------------GENERATION OF FASTQ FILE-------------------------------------------------------------------
#[1] Create df for info in fastq------------------------------------------------
#Create 1 column for all sequences
fastqSeq  <- unlist(pcrOut, recursive = TRUE)
length(fastqSeq) #4284900 = 4761 loci*180 individuals sampled*5 pcr cycles
fastqSeq[[1]] #500 BP seq

fastqTab <- as.data.frame(fastqSeq)
nrow(fastqTab) #4284900
ncol(fastqTab) #1
colnames(fastqTab) <- "sequence"
str(fastqTab)

#Generate column for sequence id
seqID <- as.data.table(rep("@seq", 4284900))
seqID
id <- 1:4284900
seqID <- cbind(seqID, id)
seqID
seqID <- as.data.frame(paste0(seqID$V1, seqID$id, sep=""))
str(seqID)
ncol(seqID) #1
nrow(seqID) #4284900
colnames(seqID) <- "seqID"
str(seqID)

#Generate column for +
plus <- as.data.table(rep("+", 4284900))
ncol(plus) #1
nrow(plus) #4284900
colnames(plus) <- "+"

#generate column for base quality scores
basequal <- paste(rep("~", 500), collapse = "")
basequal
basequalcol <- as.data.table(rep(basequal, 4284900))
ncol(basequalcol) #1
nrow(basequalcol) #4284900
colnames(basequalcol) <- "base_quality_score"

#cbind to 1 df
fastqdf <- cbind(seqID, fastqTab, plus, basequalcol)
View(fastqdf)
nrow(fastqdf) #4284900
ncol(fastqdf) #4
nchar(fastqdf[1,2]) #500

#[2] create fastq file----------------------------------------------------------
fastq <- apply(fastqdf, 1, FUN = function(x) paste0(x)) # 1 means paste is applied row-wise
typeof(fastq) #character vector
write(txt,"PATH", sep="/n")
#separator is enter/return -> each cell is put in the next line
