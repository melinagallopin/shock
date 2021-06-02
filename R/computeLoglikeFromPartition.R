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

computeLoglikeFromPartition <- function(labels, expdata, csize=NULL){
  
  if (!is.matrix(expdata) && !is.data.frame(expdata))
    stop(paste(sQuote("expdata"), "must be a matrix"))
  
  if (!is.list(labels) && !is.vector(labels))
    stop(paste(sQuote("labels"), "must be a list"))

  ## function computeLoglikeFromPartition

  ## input labels (list): block labels for each variable
  ## input expdata (matrix): data
  ## input csize (integer): sizes of connected components (recomputed if missing)
  ## output loglike (real): loglikehood of the model with the block diagonal covariance
  ## output df (integer): degree of freedom of the model
  ## output labels (list): labels provided as input
  
  ## require
  ## library(mvtnorm)
  
  covSSS <-matrix(rep(0,ncol(expdata)^2),nrow=ncol(expdata))
  if (is.null(csize)) csize = as.integer(table(labels))
  # NOTE: assuming labels start at 1 and no holes
  for (lab in 1:length(csize)) {
    indices <- which(labels==lab)
    covSSS[indices,indices] <- (1/nrow(expdata)) * t(expdata[,indices]) %*% expdata[,indices]
  }
  list(loglike=sum(dmvnorm(expdata, mean=rep(0,ncol(expdata)), sigma=covSSS, log=TRUE)),
       df=sum(csize*(csize-1)/2),
       labels=labels)
}
