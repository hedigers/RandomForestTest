

rtcopula <- function(n,R,df){
  
  p <- ncol(R)
  T <- matrix(NA, n, p)
  
  for(i in 1:n){
    r <- rmvt(1,R,df)
    for(j in 1:p){
      term1 <- pt(r[j], df)
      term2 <- qnorm(term1)
      
      T[i,j] <- term2
    }
  }
  
  return(T)
  
}