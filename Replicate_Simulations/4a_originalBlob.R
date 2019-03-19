

set.seed(747)
library(ranger)
library(mvtnorm)
library(reticulate)
library(MASS)
library(parSim)
library(hypoRF)

source("func_me.R")
source("func_mmd.R")
source("func_mmdoptimized.R")

blobSampling <- function(n, d, grid.size, Sigma) {

  # centers
  centers <- expand.grid(lapply(1:d, function(i) 1:grid.size))
  # indicators
  ind <- sample(size = n, x = 1:nrow(centers), replace = T)
  s <- sapply(ind, function(x) MASS::mvrnorm(n = 1, mu = as.numeric(centers[x,]), Sigma = Sigma))
  return(list(data=t(s),ind=ind))
}


S_power <- 500
n <- 600
p <- c(2,3) 
sigma <- 0.04
grid.size <- c(2,3)
ntrees <- 600
alpha <- 0.05
nodesize <- 4
K <- 100

start.time <- Sys.time()

cores <- 1
rep <- 1

# Run the simulation:
Results <- parSim(
  # Any number of conditions:
  n=n,
  p=p,
  grid.size=grid.size,

  # Number of repititions?
  reps = rep, # more is better!

  # Parallel?
  nCores = cores,

  # Write to file?
  write = F,

  # name="test_file",

  # Export objects (only needed when nCores > 1):
  export = c("S_power", "ntrees","alpha", "K", "sigma", "nodesize", "me_test",
             "mmd_test", "mmdoptimized_test", "hypoRF", "blobSampling"),

  # R expression:
  expression = {
    # Load all R packages in the expression if needed

    library(ranger)
    library(mvtnorm)
    library(MASS)
    library(reticulate)

    set.seed(747)

    # Want to debug? Enter browser() and run the function. Only works with nCores = 1!
    # browser()
    # Enter whatever codes you want. I can use the conditions as objects.

    n1<-n/2
    n2<-n1

    ################################# DATA GENERATION #################################
    d <- p
    nsim <- p
    R <- 2
    eigenval<-c(0,1)


    Sig <- stats::rWishart(n = 1, df = d,Sigma = diag(d))[,,1]
    R <- Matrix::cov2cor(V = Sig)
    lambda<-0.1
    eigenval<-eigen(R)$value
    i<-0

    while( min(eigenval)/max(eigenval) <= 1-1/sqrt(p) ) {

      (i<-i+1)

      R=(1-lambda)*R + lambda*diag(p)

      eigenval<-eigen(R)$value


    }

    Y <- replicate(S_power, blobSampling(n = n1, d = p, grid.size = grid.size, Sigma = sigma * diag(p))$data)
    X <- replicate(S_power, blobSampling(n = n2, d = p, grid.size = grid.size, Sigma = sigma * R)$data)


    ###################################################################################

    power_me <- mean(sapply(1:S_power, function(x){me_test(X[,,x],Y[,,x])}))
    power_mmd <- mean(sapply(1:S_power, function(x){mmd_test(X[,,x],Y[,,x])}))
    power_mmdopt <- mean(sapply(1:S_power, function(x){mmdoptimized_test(X[,,x],Y[,,x])}))
    power_binom <- mean(sapply(1:S_power, function(x){
      ifelse(hypoRF(data.frame(X[,,x]), data.frame(Y[,,x]), K = 1,
                    num.trees = ntrees, min.node.size = nodesize)$pvalue <= alpha,1,0)}))
    power_perm <- mean(sapply(1:S_power, function(x){
      ifelse(hypoRF(data.frame(X[,,x]), data.frame(Y[,,x]), K = K,
                    num.trees = ntrees, min.node.size = nodesize)$pvalue <= alpha,1,0)}))

    # Make a data frame with one row to return results (multple rows is
    # also possible but not reccomended):
    data.frame(
      power_binom, power_me, power_mmd, power_mmdopt, power_perm
    )
  }
)

end.time <- Sys.time()

(time.taken <- end.time - start.time)


############################# Saving the results #########################################

# write.table(Results, "4b_originalBlob_S500_n600_K100_node4_trees600.txt", sep=";")

# Results <- read.table("4b_originalBlob_S500_n600_K100_node4_trees600.txt", sep=";")


