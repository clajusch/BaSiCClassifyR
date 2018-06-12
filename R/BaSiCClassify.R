
# check following links:
# https://stats.stackexchange.com/questions/122857/how-to-determine-overlap-of-two-empirical-distribution-based-on-quantiles
# https://freakonometrics.hypotheses.org/4199

# well annotated example data: https://storage.googleapis.com/linnarsson-lab-www-blobs/blobs/cortex/expression_mRNA_17-Aug-2014.txt

setwd("C:/Users/scholzcl/Documents/FASTGenomics/deppenClassifyR")

expression <- read.table("Zeisel_genesymbol.tsv",
                         header = TRUE,
                         stringsAsFactors = FALSE,
                         sep = "\t")

metadata <- read.table("Zeisel et al, 2015 cells_metadata.tsv",
                       header = TRUE,
                       stringsAsFactors = FALSE,
                       sep = "\t")

metadata$cellId <- colnames(expression)



summarizeGene <- function(exprsVec, minDetect = 1) {
  cumulated <- sum(exprsVec)
  index <- exprsVec>=minDetect
  N.present <- sum(index)
  N.dropout <- sum(!index)
  if (N.present>0) {
    presentQuantiles <- quantile(exprsVec[index], probs = seq(0, 1, 0.1))
  } else {
    presentQuantiles <- quantile(exprsVec, probs = seq(0, 1, 0.1))
  }
  names(presentQuantiles) <- paste("present.Q", seq(0, 100, 10), sep ="")
  returnVec <- c(cumulated = cumulated,
                 N.present = N.present,
                 N.dropout = N.dropout,
                 presentQuantiles)
  return(returnVec)
}


summarizeCells <- function(exprsMat, minDetect = 1, cellsInColumns = TRUE) {
  shape <- if (cellsInColumns) {
    apply(exprsMat, 1, summarizeGene, minDetect = minDetect)
  } else {
    apply(exprsMat, 2, summarizeGene, minDetect = minDetect)
  }
  return(t(shape))
}


summarizePopulations <- function(exprsMat, population, minDetect = 1, cellsInColumns = TRUE) {
  groups <- unique(population)
  popSummary <- list()
  for (i in groups) {
    tmp <- if (cellsInColumns) {
      exprsMat[, population==i]
    } else {
      exprsMat[population==i,]
    }
    popSummary[[i]] <- data.frame(summarizeCells(tmp, minDetect, cellsInColumns))
  }
  return(popSummary)
}


conditionalDistance <- function(populationSummary, exprsVec, Q=50, weighted = FALSE, ignoreDropouts = TRUE, minDetect = 1) {
  Qname <- paste("present.Q", Q, sep ="")
  rawDist <- abs(populationSummary[, Qname] - exprsVec)
  if (weighted) {
    # the use of weights biases the distance towards populations with higher proportions of dropouts
    weights <- with(populationSummary, N.present/(N.present+N.dropout))
    rawDist <- rawDist * weights
  }
  if (ignoreDropouts) {
    rawDist[exprsVec<minDetect] <- 0
  }
  return(sum(rawDist))
}


distanceFromPopulations <- function(exprsVec, populationSummaries, ...) {
  distances <- sapply(populationSummaries, conditionalDistance, exprsVec = exprsVec, ...)
  return(distances)
}


classifyCells <- function(exprsMat, populationSummaries, ...) {
  distMat <- t(apply(exprsMat, 2, distanceFromPopulations, populationSummaries=populationSummaries))
  classTable <- data.frame(cellId = rownames(distMat),
                           population = apply(distMat, 1, function(x) names(which.min(x))),
                           postProb = apply(distMat, 1, function(x) 1-(min(x)/sum(x))))
}




zeiselLevel1 <- summarizePopulations(expression, population = metadata$level1class)
classLevel1 <- classifyCells(expression, populationSummaries = zeiselLevel1)
classLevel1$truth <- metadata$level1class

as.matrix(ftable(truth~population, data=classLevel1))->m
pheatmap(m, cluster_rows = F, cluster_cols = F)

p <- t(t(m)/colSums(m))
pheatmap(p, cluster_rows = F, cluster_cols = F)


zeiselLevel2 <- summarizePopulations(expression, population = metadata$level2class)
classLevel2 <- classifyCells(expression, populationSummaries = zeiselLevel2)
classLevel2$truth <- metadata$level2class

as.matrix(ftable(truth~population, data=classLevel2))->m
# pheatmap(m, cluster_rows = F, cluster_cols = F)

p <- t(t(m)/colSums(m))
pheatmap(p, cluster_rows = F, cluster_cols = F)
