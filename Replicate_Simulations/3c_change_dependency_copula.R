

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
source("FaK_dist.R")

S_power <- 200
n <- 600
p <- 60
ntrees <- 600
alpha <- 0.05
nodesize <- 4 
K <- 100
df <- seq(1,8, by=0.5)
R <- diag(p)

start.time <- Sys.time()

cores <- 1
rep <- 1

# Run the simulation:
Results <- parSim(
  # Any number of conditions:
  n=n,
  p=p,
  df=df,
  
  # Number of repititions?
  reps = rep, # more is better!
  
  # Parallel?
  nCores = cores,
  
  # Write to file?
  write = F,
  
  # name="test_file",
  
  # Export objects (only needed when nCores > 1):
  export = c("S_power", "ntrees","alpha", "K", "rtcopula", "R",
             "nodesize", "me_test", "mmd_test", "mmdoptimized_test", "hypoRF"),
  
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
    X <- replicate(S_power, replicate(p, rnorm(n1)), simplify = F)
    Y <- replicate(S_power, rtcopula(n2,R,df), simplify = F)
    ###################################################################################
    
    power_me <- mean(sapply(1:S_power, function(x){me_test(X[[x]],Y[[x]])}))
    power_mmd <- mean(sapply(1:S_power, function(x){mmd_test(X[[x]],Y[[x]])}))
    power_mmdopt <- mean(sapply(1:S_power, function(x){mmdoptimized_test(X[[x]],Y[[x]])}))
    power_binom <- mean(sapply(1:S_power, function(x){
      ifelse(hypoRF(data.frame(X[[x]]), data.frame(Y[[x]]), K = 1, 
                    num.trees = ntrees, min.node.size = nodesize)$pvalue <= alpha,1,0)}))
    power_perm <- mean(sapply(1:S_power, function(x){
      ifelse(hypoRF(data.frame(X[[x]]), data.frame(Y[[x]]), K = K,
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

# write.table(Results, "3c_copula_S200_p60_trees600_node4_K100.txt", sep=";")
# 
# Results <- read.table("3c_copula_S200_p60_trees600_node4_K100.txt", sep=";")
# 
# library(tikzDevice)
# options(tikzDefaultEngine = 'xetex')
# tikz(file = "3c_plot_copula_K100.tex", width=6, height=4)
# 
# par(oma = c(1, 1, 2, 1))
# plot(Results$df, Results$power_perm, xlab='$v$',
#      ylab="Power", pch=19, type="b", ylim=c(0,1), col="red4")
# lines(Results$df, Results$power_binom, col="tomato", pch=17, type="b")
# lines(Results$df, Results$power_me, col="skyblue3", pch=15, type="b")
# lines(Results$df, Results$power_mmd, col="snow4", pch=3, type="b")
# lines(Results$df, Results$power_mmdopt, col="palegreen4", pch=4, type="b")
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
#        cex=1, bty="n", bg="transparent", text.width=c(1,1,1,0.75,1))
# 
# dev.off()

