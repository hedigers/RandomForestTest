

set.seed(747)
library(ranger)
library(mvtnorm)
library(reticulate)
library(MASS)
library(hypoRF)

n <- 100
p <- 10 
ntrees <- 600
alpha <- 0.05
nodesize <- 4
K <- c(10,20,40,60,80,100,150,200,500,700,1000)

start.time <- Sys.time()

library(parSim)

cores <- 1
rep <- 1
S <- 100

# Run the simulation:
Results <- parSim(
  # Any number of conditions:
  n=n,
  K=K,
  p=p,
  
  # Number of repititions?
  reps = rep, # more is better!
  
  # Parallel?
  nCores = cores,
  
  # Write to file?
  write = F,
  
  # name="test_file",
  
  # Export objects (only needed when nCores > 1):
  export = c("S", "ntrees", "alpha", "nodesize",
             "hypoRF"),
  
  # R expression:
  expression = {
    # Load all R packages in the expression if needed
    
    library(ranger)
    library(mvtnorm)
    library(reticulate)
    library(MASS)
    
    set.seed(747)
    
    # Want to debug? Enter browser() and run the function. Only works with nCores = 1!
    # browser()
    # Enter whatever codes you want. I can use the conditions as objects.
    
    var_permRF <- replicate(S,hypoRF(
      data1=data.frame(replicate(p,rnorm(n))),
      data2=data.frame(replicate(p,rnorm(n))), K=K, num.trees=ntrees,
      min.node.size=nodesize)[[4]])
    
    # Make a data frame with one row to return results (multple rows is
    # also possible but not reccomended):
    data.frame(
      var_permRF
    )
  }
)

end.time <- Sys.time()

(time.taken <- end.time - start.time)


############################# Saving the results and generating plots #######################

# write.table(Results, "var_check_permRF.txt", sep=";")
# 
# Results <- read.table("var_check_permRF.txt", sep=";", header=T)
# 
# library(tikzDevice)
# options(tikzDefaultEngine = 'xetex')
# tikz(file = "plot_varcheck.tex", width=6, height=4)
# 
# plot(Results$K, Results$var_permRF, xlab='$K$',
#      ylab="Estimated Variance", pch=19, cex=0.8)
# abline(h=0.05, lty="dashed")
# 
# dev.off()


