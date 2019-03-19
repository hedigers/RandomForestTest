

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

threeModeGaussian <- function(n, mean.vec, sd.vec, weights) {
  ind <- sample(x = 0:2, prob = weights, replace = T, size = n)
  return(ifelse(ind == 0, rnorm(n = n, mean = mean.vec[1], sd = sd.vec[1]),
                ifelse(ind == 1, rnorm(n = n, mean = mean.vec[2], sd = sd.vec[2]),
                       rnorm(n = n, mean = mean.vec[3], sd = sd.vec[3]))))
}

multivariateBlob <- function(n, p, weights, mean.vec.h0, mean.vec.h1, sd.vec.h0, sd.vec.h1, ind.h1) {
  l <- lapply(1:p, function(v) {
    if (v %in% ind.h1) {
      return(threeModeGaussian(n = n, mean.vec =  mean.vec.h1, sd.vec = sd.vec.h1, weights = weights))
    } else {
      return(threeModeGaussian(n = n, mean.vec =  mean.vec.h0, sd.vec = sd.vec.h0, weights = weights))
    }
  })
  return(Reduce(x = l, f = cbind))
}

S_power <- 200
n <- 600 
p <- seq(2,32,by=2)
mean.vec.h0 <- c(-5,0,5)
mean.vec.h1 <- c(-5,0,5)
sd.vec.h0 <- c(1,1,1)
sd.vec.h1 <- c(1,2,1)
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

  # Number of repititions?
  reps = rep, # more is better!

  # Parallel?
  nCores = cores,

  # Write to file?
  write = F,

  # name="test_file",

  # Export objects (only needed when nCores > 1):
  export = c("S_power", "ntrees","alpha", "K", "mean.vec.h0",
             "mean.vec.h1", "sd.vec.h0", "sd.vec.h1", "nodesize", "me_test",
             "mmd_test", "mmdoptimized_test", "hypoRF", "threeModeGaussian",
             "multivariateBlob"),
  
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
    Y <- replicate(S_power, multivariateBlob(n = n1, p = p, weights = c(1/3, 1/3, 1/3), mean.vec.h0 = mean.vec.h0,
                                             mean.vec.h1 =mean.vec.h1, sd.vec.h0 =  sd.vec.h0 , sd.vec.h1 = sd.vec.h1,
                                             ind.h1 = c()))
    
    X <- replicate(S_power, multivariateBlob(n = n2, p = p, weights = c(1/3, 1/3, 1/3), mean.vec.h0 = mean.vec.h0,
                                             mean.vec.h1 = mean.vec.h1, sd.vec.h0 =  sd.vec.h0 , sd.vec.h1 = sd.vec.h1,
                                             ind.h1 = c(1:p)))
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


############################# Saving the results and generating plots #######################

# write.table(Results, "4_multivariateblob_S200_n600_K100_node4_trees600_psmall.txt", sep=";")
# 
# Results_s <- read.table("4_multivariateblob_S200_n600_K100_node4_trees600_psmall.txt", sep=";")
# 
# library(tikzDevice)
# options(tikzDefaultEngine = 'xetex')
# tikz(file = "4_plot_Blob_K100_small.tex", width=6, height=4)
# 
# par(oma = c(1, 1, 2, 1))
# plot(Results_s$p, Results_s$power_perm, xlab='$p$',
#      ylab="Power", pch=19, type="b", ylim=c(0,1), col="red4")
# lines(Results_s$p, Results_s$power_binom, col="tomato", pch=17, type="b")
# lines(Results_s$p, Results_s$power_me, col="skyblue3", pch=15, type="b")
# lines(Results_s$p, Results_s$power_mmd, col="snow4", pch=3, type="b")
# lines(Results_s$p, Results_s$power_mmdopt, col="palegreen4", pch=4, type="b")
# abline(h=0.05, lty="dashed")
# 
# leg <- numeric()
# leg[1] <- "HypoRF"
# leg[2] <- "Binomial"
# leg[3] <- "ME-full"
# leg[4] <- "MMDboot"
# leg[5] <- "MMD-full"
# 
# par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
# legend("top", legend=leg, xpd=T, horiz=T, inset=c(0,0),
#        col=c("red4","tomato", "skyblue3", "snow4", "palegreen4"),
#        pch=c(19,17,15,3,4),
#        cex=1, bty="n", bg="transparent", text.width=c(4,4,4,4,4.5))
# 
# dev.off()



