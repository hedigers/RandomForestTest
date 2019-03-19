

rcont <- function(n, p, dseq, lambda){
  
  dec <- rbinom(n, size=1, prob=(1-lambda))
  
  baseset <- replicate(p, rnorm(n, 50, 5))
  
  cont_set <- sapply(1:n, function(x){
    if(dec[x]==1){
      baseset[x,] <- replicate(p, rnorm(1, 50, 5))
    } else {
      baseset[x, dseq] <- replicate(length(dseq), rbinom(1, size=100, prob=0.5))
    }
    return(baseset[x,])
  }
  )
  
  return(t(cont_set))
  
}

