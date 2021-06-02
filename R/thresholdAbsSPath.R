#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

thresholdAbsSPath <- function(expdata, maxBlockSize=50, detectDuplicates=FALSE){

  ## function that take data and provide the clustering output
  ## based on the threshold sample covariance matrix

  ## input expdata (matrix): data to compute the sample covariance matrix S
  ## input maxBockSize (integer): maximum size of a block
  ## input detectDuplicates (bool): detect redundant partitions
  ## output partitionList (list): list of resulting partitions from thresolding S
  ## output lambdaPath (list or numeric): list/vector of threshold parameters

  ## require
  ## library(igraph)

  if(is.matrix(expdata) == FALSE & is.data.frame(expdata) == FALSE)
    stop(paste(sQuote("expdata"), "must be a matrix"))
  detectDuplicates <- detectDuplicates && require(clusterCrit, quietly=TRUE)

  Sabs <- abs(cor(expdata))
  valThres <- Sabs[upper.tri(Sabs)]
  lambdaValue <- sort(unique(valThres), decreasing=TRUE)

  # No need to recompute E at each loop:
  E <- matrix(0, nrow=nrow(Sabs), ncol=ncol(Sabs))
  listPartitions <- list()
  listLambdas <- list()
  prevNbCsize <- 0
  curLastIdx <- 0
  skipped <- 0
  for (lam in seq_along(lambdaValue)) {
    E[Sabs > lambdaValue[lam]] <- 1
    goutput <- graph.adjacency(E, mode="undirected", weighted=NULL)
    cls <- clusters(goutput)
    membership <- as.integer(cls$membership)
    # Stop if a connected component is larger than maxBlockSize:
    if (any(cls$csize > maxBlockSize)) break
    if (detectDuplicates) {
      # Skip if same partition than before (but add lambda anyway):
      if (
        curLastIdx >= 1 && length(cls$csize) == prevNbCsize &&
        clusterCrit::extCriteria(listPartitions[[curLastIdx]]$membership, membership, 'Jaccard') == 1
      ) {
        skipped <- skipped + 1
        listLambdas[[curLastIdx]] <- c(listLambdas[[curLastIdx]], lambdaValue[lam])
        next
      }
    }
    # Else: initialize listLambdas + fill listPartitions
    curLastIdx <- curLastIdx + 1
    listLambdas[[curLastIdx]] <- lambdaValue[lam]
    listPartitions[[curLastIdx]] <- list(membership=membership, csize=as.integer(cls$csize))
    prevNbCsize <- length(cls$csize)
  }

  if (!detectDuplicates || skipped == 0) listLambdas = unlist(listLambdas)
  return(list(partitionList=listPartitions, lambdaPath=listLambdas))
}
