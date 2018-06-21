
# check following links:
# https://stats.stackexchange.com/questions/122857/how-to-determine-overlap-of-two-empirical-distribution-based-on-quantiles
# https://freakonometrics.hypotheses.org/4199


setwd("C:/Users/scholzcl/Documents/Git/BaSiCClassifyR")

load("data/zeisel2015.RData")


#' Summarize a gene expression vector
#'
#' @param exprsVec A vector of expression values.
#' @param minDetect The minimum expression value of a detected gene. Defaults to 1.
#' @return A vector of summary statistics.
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


#' Generate a gene expression summary profile for a cell poulation.
#'
#' @param exprsMat The gene expression matrix.
#' @param minDetect The minimum expression value of a detected gene. Defaults to 1.
#' @param cellsInColumns Indicates if cells are stored in columns of matrix. Defaults to TRUE.
#' @return A gene expression summary profile.
summarizeCells <- function(exprsMat, minDetect = 1, cellsInColumns = TRUE) {
  shape <- if (cellsInColumns) {
    apply(exprsMat, 1, summarizeGene, minDetect = minDetect)
  } else {
    apply(exprsMat, 2, summarizeGene, minDetect = minDetect)
  }
  return(t(shape))
}


#' Calculate gene expression summary profiles for several cell populations.
#'
#' @param exprsMat The gene expression matrix.
#' @param population A vector of population assignments. Must equal the number of cells in dataset.
#' @param minDetect The minimum expression value of a detected gene. Defaults to 1.
#' @param cellsInColumns Indicates if cells are stored in columns of matrix. Defaults to TRUE.
#' @return A list of gene expression summary profiles.
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


#' Calculate the distance of a drop-out-inflated single-cell expression profile to a gene expression summary profile of a cell population.
#'
#' @param populationSummary
#' @param exprsVec
#' @param Q
#' @param weighted
#' @param ignoreDropouts = TRUE,
#' @param minDetect The minimum expression value of a detected gene. Defaults to 1.
#' @return The drop-out conditional distance.
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


splitTrainAndTest <- function(exprsMat, population, trainProp = 0.8, cellsInColumns = TRUE) {
  splitData <- list(training = list(), test = list())
  cellIDs <- if (cellsInColumns) {
    colnames(exprsMat)
  } else {
    rownames(exprsMat)
  }
  cellIDs <- split(cellIDs, population)
  trainIDs <- unlist(lapply(cellIDs, function(x, p) sample(x, round(length(x)*p)), p=trainProp))
  splitData$training <- if (cellsInColumns) {
    exprsMat[, colnames(exprsMat) %in% trainIDs]
  } else {
    exprsMat[rownames(exprsMat) %in% trainIDs,]
  }
  splitData$test <- if (cellsInColumns) {
    exprsMat[, !colnames(exprsMat) %in% trainIDs]
  } else {
    exprsMat[!rownames(exprsMat) %in% trainIDs,]
  }
  return(splitData)
}


trainClassifyR <- function(exprsMat, population, minDetect = 1, cellsInColumns = TRUE) {

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
