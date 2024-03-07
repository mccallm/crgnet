topology <- function(fits, weighted=FALSE){
  ngene <- nrow(fits[[1]]$parents)
  nexp <- length(fits[[1]]$trajectories)
  grph <- matrix(0, nrow=ngene, ncol=nexp)
  for(k in 1:length(fits)){
    gtmp <- t(fits[[k]]$parents + 1)
    for(j in 1:ncol(gtmp)){
      if(weighted){
        grph[j, gtmp[,j]] <- grph[j, gtmp[,j]] + 1/fits[[k]]$unnormalized_score
      } else{
        grph[j, gtmp[,j]] <- grph[j, gtmp[,j]] + 1
      }
    }
  }
  if(weighted){
    cgrph = grph / sum(sapply(fits, function(x) x$unnormalized_score))
  } else{
    cgrph = grph / length(fits)
  }
  return(cgrph)
}
