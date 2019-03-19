
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
source("contaminated_dist.R")

S_power <- 200
n <- 600
p <- 200
ntrees <- 600
alpha <- 0.05
nodesize <- 4 
K <- 100
lambda <- seq(0.5,1, by=0.05)
dseqind <- c(1,2)
sequences <- list(small=sample(1:p, round(0.1 * p)), all=sample(1:p, round(p)))

start.time <- Sys.time()

cores <- 1
rep <- 1

# Run the simulation:
Results <- parSim(
  # Any number of conditions:
  n = n,
  p = p,
  lambda = lambda,
  dseqind = dseqind,
  
  # Number of repititions?
  reps = rep, # more is better!
  
  # Parallel?
  nCores = cores,
  
  # Write to file?
  write = F,
  
  # name="test_file",
  
  # Export objects (only needed when nCores > 1):
  export = c("S_power", "ntrees","alpha", "K",
             "nodesize", "me_test", "mmd_test", "mmdoptimized_test",
             "hypoRF", "rcont", "sequences"),
  
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
    dseq <- sequences[[dseqind]]
    
    X <- replicate(S_power, rcont(n1, p, dseq, 0), simplify = F)
    Y <- replicate(S_power, rcont(n2, p, dseq, lambda), simplify = F)
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

# write.table(Results, "4_results_cont_S200_trees600_node4_K100.txt", sep=";")
# 
# Results <- read.table("4_results_cont_S200_trees600_node4_K100.txt", sep=";")
# 
# library(tikzDevice)
# options(tikzDefaultEngine = 'xetex')
# tikz(file = "4_plot_contamination_K200.tex", width=6, height=4)
# 
# par(oma = c(1, 1, 2, 1))
# plot(Results$lambda, Results$power_perm, xlab='$\\lambda$',
#     ylab="Power", pch=19, type="b", ylim=c(0,1), col="red4")
# lines(Results$lambda, Results$power_binom, col="tomato", pch=17, type="b")
# lines(Results$lambda, Results$power_me, col="skyblue3", pch=15, type="b")
# lines(Results$lambda, Results$power_mmd, col="snow4", pch=3, type="b")
# lines(Results$lambda, Results$power_mmdopt, col="palegreen4", pch=4, type="b")
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
#        cex=1, bty="n", bg="transparent", text.width=c(0.075,0.075,0.075,0.075,0.08))
# 
# dev.off()






