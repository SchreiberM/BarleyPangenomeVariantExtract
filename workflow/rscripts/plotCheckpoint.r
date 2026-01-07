#!/usr/bin/env Rscript

##############################################################
# Script for size extraction

library(tidyverse)
library(Biostrings)
library(ggplot2)

Args <- commandArgs(trailingOnly=TRUE)

currentDir <- getwd()
outputDir <- unlist(strsplit(Args[3], "/\\s*(?=[^/]+$)", perl=TRUE))[1]

files <- list.files(paste(currentDir, "/", outputDir, "/blastResult/", sep=""))

dir.create(file.path(paste(currentDir, "/", Args[3], sep="")))
dir.create(file.path(paste(currentDir, "/", outputDir, "/report", sep="")))


distancepromoter <- as.integer(Args[2])
geneIdentity <- as.integer(Args[5])
geneLength <- as.integer(Args[7])
geneTotalLength <- as.integer(Args[6])

s <- readDNAStringSet(paste(currentDir, "/", Args[4] , sep=""))
sizeRef <- width(s)

#print(distancepromoter)
furtherChecks <- list()
noInformation <- list()

fileCounter <- length(files)

for(i in 1:length(files)){
    temp <- tryCatch(read.table(paste(currentDir, "/", outputDir, "/blastResult/", files[i], sep=""), header=FALSE), error=function(e) NA)
    if(is.null(dim(temp))){
        name <- unlist(str_split(files[i], "\\."))[1]
        bedfile <- cbind(0, 0, 0, name, 0, ".")
        write.table(bedfile, file=paste(currentDir, "/", outputDir, "/bedFiles/",name,".bed", sep=""), col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
        fileCounter <- fileCounter-1
        noInformation <- append(noInformation, files[i])

    }else if(dim(temp)[1] != 1){
        
        secondFilter <- temp %>% filter(V3>=geneIdentity & V16 >= geneLength)
        if(dim(secondFilter)[1] != 1){
            furtherChecks <- append(furtherChecks, files[i])
        }
    }
}

if(length(-match(c(unlist(furtherChecks), unlist(noInformation)),files))==0){
    newList <- files
}else{
    newList <- files[-match(c(unlist(furtherChecks), unlist(noInformation)),files)]
}

correctMerged <- cbind("chromosome", "start", "end", "genome", "score", "strand", "colour")
colnames(correctMerged) <- c("chromosome", "start", "end", "genome", "score", "strand", "colour")

additionalMerged <- cbind("chromosome", "start", "end", "genome", "score", "strand", "colour")
colnames(additionalMerged) <- c("chromosome", "start", "end", "genome", "score", "strand", "colour")

if(length(newList)>0){
    for(k in 1:length(newList)){
        readBlastResults <- read.table(paste(currentDir, "/", outputDir, "/blastResult/", newList[[k]], sep=""), header=FALSE)
        temp <- readBlastResults %>% filter(V3>=geneIdentity & V16 >= geneTotalLength)
        name <- unlist(str_split(newList[k], "\\."))[1]

        if(dim(temp)[1]==0){
            bedfile <- cbind(0, 0, 0, name, 0, ".")
            write.table(bedfile, file=paste(currentDir, "/", outputDir, "/bedFiles/",name,".bed", sep=""), col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
            next
        }

        if(min(temp[9:10])<distancepromoter){
            minPos <- 0
            maxPos <- max(temp[9:10])+distancepromoter
        }else{
            minPos <- min(temp[9:10])-distancepromoter
            maxPos <- max(temp[9:10])+distancepromoter
        }
        if(temp[,17]=="minus"){
            bedfile <- cbind(temp[2], (minPos), (maxPos), name, 0, "-")
        }else{
            bedfile <- cbind(temp[2], (minPos), (maxPos), name, 0, "+")
        }
        checkpoint <- bedfile   
        toKeep <- (as.integer(checkpoint[,3])-as.integer(checkpoint[,2]))-(2*distancepromoter) > sizeRef/2
        colnames(bedfile) <- c("chromosome", "start", "end", "genome", "score", "strand")
        bedfile$start[bedfile$start<0] <- 0
        write.table(bedfile[toKeep,,drop=FALSE], file=paste(currentDir, "/", outputDir, "/bedFiles/",name,".bed", sep=""), col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
        bedfile$colour <- "#33a02c"
        correctMerged <- rbind(correctMerged, bedfile)
    }
}

correctMerged <- correctMerged[-1,]

############# 

if(length(furtherChecks)>0){
    for(j in 1:length(furtherChecks)){
        readBlastResults <- read.table(paste(currentDir, "/", outputDir, "/blastResult/", furtherChecks[[j]], sep=""), header=FALSE)
        temp <- readBlastResults %>% filter(V3>=geneIdentity & V16 >= geneLength)
        name <- unlist(str_split(furtherChecks[[j]], "\\."))[1]

        if(dim(temp)[1]==0){
            bedfile <- cbind(0, 0, 0, name, 0, ".")
            write.table(bedfile, file=paste(currentDir, "/", outputDir, "/bedFiles/",name,".bed", sep=""), col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
            next
        }

        bedfile <- cbind("chromosome", "start", "end", "genome", "score", "strand")
        colnames(bedfile) <- c("chromosome", "start", "end", "genome", "score", "strand")

        if(all(temp[,16]>geneTotalLength)){
            for(m in 1:dim(temp)[1]){
                if(min(temp[m,9:10])<distancepromoter){
                    minPos <- 0
                    maxPos <- max(temp[m,9:10])+distancepromoter
                }else{
                    minPos <- min(temp[m,9:10])-distancepromoter
                    maxPos <- max(temp[m,9:10])+distancepromoter
                }
                if(temp[m,17]=="minus"){
                    newOut <- cbind(temp[m,2], (minPos), (maxPos), paste(name, m, sep="_"), 0, "-")
                }else{
                    newOut <- cbind(temp[m,2], (minPos), (maxPos), paste(name, m, sep="_"), 0, "+")
                }
                bedfile <- rbind(bedfile, newOut)
            }

            bedfile <- as.data.frame(bedfile[-1,])
            checkpoint <- bedfile
            toKeep <- (as.integer(checkpoint[,3])-as.integer(checkpoint[,2]))-(2*distancepromoter) > sizeRef/2
            bedfile$start[bedfile$start<0] <- 0
            write.table(bedfile[toKeep,,drop=FALSE], file=paste(currentDir, "/", outputDir, "/bedFiles/",name,".bed", sep=""), col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
            checkpoint$colour <- "#1f78b4"
            additionalMerged <- rbind(additionalMerged, checkpoint)
            next
        }

        if(length(table(temp[,2]))==1){
            estimatedSize <- dim(temp)[1]
            df <- data.frame(ID = 1:estimatedSize, scores = temp[,9])
            diffMatrix <- abs(outer(df$scores, df$scores, `-`))
            result <- diffMatrix < (sizeRef + (sizeRef/2))
            for(a in 1:estimatedSize){
                variableQ <- result[a,1:estimatedSize]
                if(sum(variableQ)>1){
                    pos <- which(result[a,]==TRUE)
                    if(min(temp[c(a,pos),9:10])<distancepromoter){
                        minPos <- 0
                        maxPos <- max(temp[c(a,pos),9:10])+distancepromoter
                    }else{
                        minPos <- min(temp[c(a,pos),9:10])-distancepromoter
                        maxPos <- max(temp[c(a,pos),9:10])+distancepromoter
                    } 
                    if(temp[a,17]=="minus"){
                        tempOut <- cbind(temp[a,2], (minPos),(maxPos), paste(name, a, sep="_"), 0, "-")
                    }else{
                        tempOut <- cbind(temp[a,2], (minPos),(maxPos), paste(name, a, sep="_"), 0, "+")
                    }
                    colnames(tempOut) <- c("chromosome", "start", "end", "genome", "score", "strand")
                }else{
                    if(min(temp[a,9:10])<distancepromoter){
                        minPos <- 0
                        maxPos <- max(temp[a,9:10])+distancepromoter
                    }else{
                        minPos <- min(temp[a,9:10])-distancepromoter
                        maxPos <- max(temp[a,9:10])+distancepromoter
                    } 

                    if(temp[a,17]=="minus"){
                        tempOut <- cbind(temp[a,2], (minPos), (maxPos), paste(name, a, sep="_"), 0, "-")
                    }else{
                        tempOut <- cbind(temp[a,2], (minPos), (maxPos), paste(name, a, sep="_"), 0, "+")
                    }
                    colnames(tempOut) <- c("chromosome", "start", "end", "genome", "score", "strand")
                }
                bedfile <- rbind(bedfile, tempOut)
            }

        }else{
            ####### sort multiple chromosome locations next, split matrix into individual chromosome groupings
            numberChromosome <- unique(temp[,2])
            for(b in 1:length(numberChromosome)){
                chromosomeSubset <- temp[temp[,2]==numberChromosome[b],]
                estimatedSize <- dim(chromosomeSubset)[1]
                df <- data.frame(ID = 1:estimatedSize, scores = chromosomeSubset[,9])
                diffMatrix <- abs(outer(df$scores, df$scores, `-`))
                result <- diffMatrix < (sizeRef + (sizeRef/2))
                for(a in 1:estimatedSize){
                    variableQ <- result[a,1:estimatedSize]
                    if(sum(variableQ)>1){
                        pos <- which(result[a,]==TRUE)
                        if(min(chromosomeSubset[c(a,pos),9:10])<distancepromoter){
                            minPos <- 0
                            maxPos <- max(chromosomeSubset[c(a,pos),9:10])+distancepromoter
                        }else{
                            minPos <- min(chromosomeSubset[c(a,pos),9:10])-distancepromoter
                            maxPos <- max(chromosomeSubset[c(a,pos),9:10])+distancepromoter
                        } 
                        if(chromosomeSubset[a,17]=="minus"){
                            tempOut <- cbind(chromosomeSubset[a,2], (minPos),(maxPos), paste(name, numberChromosome[b], a, sep="_"), 0, "-")
                        }else{
                            tempOut <- cbind(chromosomeSubset[a,2], (minPos),(maxPos), paste(name, numberChromosome[b], a, sep="_"), 0, "+")
                        }
                        colnames(tempOut) <- c("chromosome", "start", "end", "genome", "score", "strand")
                    }else{
                        if(min(chromosomeSubset[a,9:10])<distancepromoter){
                            minPos <- 0
                            maxPos <- max(chromosomeSubset[a,9:10])+distancepromoter
                        }else{
                            minPos <- min(chromosomeSubset[a,9:10])-distancepromoter
                            maxPos <- max(chromosomeSubset[a,9:10])+distancepromoter
                        } 
                        if(chromosomeSubset[a,17]=="minus"){
                            tempOut <- cbind(chromosomeSubset[a,2], (minPos), (maxPos), paste(name, numberChromosome[b], a, sep="_"), 0, "-")
                        }else{
                            tempOut <- cbind(chromosomeSubset[a,2], (minPos), (maxPos), paste(name, numberChromosome[b], a, sep="_"), 0, "+")
                        }
                        colnames(tempOut) <- c("chromosome", "start", "end", "genome", "score", "strand")
                    }
                    bedfile <- rbind(bedfile, tempOut)
                }
            }
        }

        bedfile <- as.data.frame(bedfile[-1,])
        extractedUnique <- bedfile[!duplicated(bedfile[,1:3]),]
        extractedUnique$colour <- "#1f78b4"
        toKeep <- (as.integer(extractedUnique[,3])-as.integer(extractedUnique[,2]))-(2*distancepromoter) > sizeRef/2
        extractedUnique$start[extractedUnique$start<0] <- 0
        write.table(extractedUnique[toKeep,,drop=FALSE], file=paste(currentDir, "/", outputDir, "/bedFiles/",name,".bed", sep=""), col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
        additionalMerged <- rbind(additionalMerged, extractedUnique)
    }
}

additionalMerged <- as.data.frame(additionalMerged[-1,])

correctMerged <- rbind(correctMerged, additionalMerged)

input <- correctMerged

colnames(input) <- c("Chromosome", "Start", "End", "GenomeIdentifier", "Score", "Strand", "colour")

input$size <- abs(as.integer(input$End)-as.integer(input$Start))

test <- cbind(input$GenomeIdentifier, input$size, input$colour)
colnames(test) <- c("Genotype", "Size", "colour")

test <- as.data.frame(test)
test$Size <- as.integer(test$Size)

toKeep_p2 <- input$size-(2*distancepromoter) < sizeRef/2

write.table(input[toKeep_p2,], file=paste(currentDir, "/", outputDir, "/report/BlastResults_RemovedSequences.txt", sep=""), col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
write.table(input[!toKeep_p2,], file=paste(currentDir, "/", outputDir, "/report/BlastResults_RetainedSequences.txt", sep=""), col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
test$grouping <- NA
test$grouping[!toKeep_p2] <- "Kept Genome Samples"
test$grouping[toKeep_p2] <- "Removed Genome Samples"

test$grouping <- as.factor(test$grouping)

x <- test[order(test$Size),]

if(dim(x)[1]<=76){
    scaleFactor <- 0.42
}else if(dim(x)[1]<=182){
    toSubstract <- (dim(x)[1]/76)/10
    scaleFactor <- 0.42-toSubstract
}else{
    scaleFactor <- 0.18
}

pdf(Args[1])
dotchart(x$Size, labels=x$Genotype,  cex=scaleFactor, groups = x$grouping, color = x$colour, main=paste("Checkpoint - Size of gene on reference genome.\n", fileCounter," out of ", length(files), " Genomes with blast hit.\n Settings: Percentage Identity ", geneIdentity, " and Percentage Coverage: ", geneLength, sep=""), sub="Green = one blast result, blue = multiple blast results", ylab="Genotype", xlab="Size in bp")
dev.off()



