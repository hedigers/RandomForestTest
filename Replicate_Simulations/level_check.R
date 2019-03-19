

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
source("FaK_dist.R")

n <- 600
p <- 200
ntrees <- 600
alpha <- 0.05
nodesize <- 4
K <- 100


dist_help <- function(n,p,dist_info){
  
  dist <- dist_info[[1]]
  arg.list <- dist_info[[2]]
  multi <- dist_info[[3]]
  arg.list$n <- n
  
  if(multi == "multivariate"){
    sample <- do.call(dist,arg.list)
  } else {
    sample <- replicate(p,do.call(dist,arg.list))
  }
  
  return(sample)
}

rmixture <- function(n, prob, mu1, mu2, Sigma){
  bern <- rbinom(n,1,prob)
  mixture <- sapply(bern, function(x){
    if(x==1){
      mvrnorm(n=1,mu=mu1,Sigma=Sigma)
    } else {
      mvrnorm(n=1,mu=mu2,Sigma=Sigma)
    }})
  
  mixture <- t(mixture)
  return(mixture)
}

cov <- matrix(0.95, nrow=p, ncol=p)
diag(cov) <- 1
dseq <- sample(1:p, 1)

distributions <- list(
  list("rnorm", list(n=n, mean=0, sd=1), "univariate"),
  list("rbinom", list(n=n, size=4, prob=0.7), "univariate"),
  list("rt", list(n=n, df=1), "univariate"),
  list("rpois", list(n=n, lambda=4), "univariate"),
  list("rf", list(n=n, df1=4, df2=12), "univariate"),
  list("runif", list(n=n, min=3, max=10), "univariate"),
  list("rmvt", list(n=n, sigma=diag(p), df=1), "multivariate"),
  list("rexp", list(n=n, rate=4), "univariate"),
  list("rbeta", list(n=n, shape1=2, shape2=3), "univariate"),
  list("mvrnorm", list(n=n, mu=rep(0,p), Sigma=cov), "multivariate"),
  list("rlnorm", list(n=n), "univariate"),
  list("rtcopula", list(n=n, R=diag(p), df=1), "multivariate"),
  list("rweibull", list(n=n, shape=1), "univariate"),
  list("rmixture", list(n=n, prob=0.9, mu1=rep(0,p), mu2=rep(4,p), Sigma=diag(p)),
       "multivariate"),
  list("rcont", list(n=n, p=p, dseq=dseq, lambda=0.6), "multivariate")
)

start.time <- Sys.time()

cores <- 1
rep <- 1
S_power <- 200
ids <- seq(1,length(distributions))

# Run the simulation:
Results <- parSim(
  # Any number of conditions:
  n=n,
  # mean_alt=mean_alt,
  p=p,
  ids=ids,
  
  # Number of repititions?
  reps = rep, # more is better!
  
  # Parallel?
  nCores = cores,
  
  # Write to file?
  write = F,
  
  # name="test_file",
  
  # Export objects (only needed when nCores > 1):
  export = c("S_power", "hypoRF", "dist_help", "rtcopula", "rcont", "K",
             "me_test", "mmd_test", "mmdoptimized_test","alpha", "ntrees",
             "distributions", "nodesize", "rmixture"),
  
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
    distlist <- distributions[[ids]]
    
    X <- replicate(S_power, dist_help(n1,p,distlist), simplify = F)
    Y <- replicate(S_power, dist_help(n2,p,distlist), simplify = F)
    ###################################################################################

    power_binom <- mean(sapply(1:S_power, function(x){
      ifelse(hypoRF(data.frame(X[[x]]), data.frame(Y[[x]]), K = 1, 
                    num.trees = ntrees, min.node.size = nodesize)$pvalue <= alpha,1,0)}))
    power_perm <- mean(sapply(1:S_power, function(x){
      ifelse(hypoRF(data.frame(X[[x]]), data.frame(Y[[x]]), K = K,
                    num.trees = ntrees, min.node.size = nodesize)$pvalue <= alpha,1,0)}))
    
    # Make a data frame with one row to return results (multple rows is
    # also possible but not reccomended):
    # data.frame(
    #   power_me, power_mmd, power_mmdopt, power_binom,  power_perm
    # )
    data.frame(
       power_binom,  power_perm
    )
  }
)

end.time <- Sys.time()

(time.taken <- end.time - start.time)


############################# Saving the results and generating plots #######################

# dists <- sapply(c(1:length(ids)),function(x){distributions[[x]][[1]]})
# 
# Results$dists <- dists
# write.table(Results, "level_check.txt", sep=";")
# 
# Results <- read.table("level_check.txt", sep=";")
# 
# library(tikzDevice)
# options(tikzDefaultEngine = 'xetex')
# tikz(file = "level_check_binom.tex", width=6, height=4)
# 
# par(xaxt="n")
# plot(Results$id, Results$power_binom, xlab='',
#      ylab="Actual Size", pch=19, ylim=c(0,0.1))
# abline(h=0.05, lty="dashed")
# lablist <- Results$dists
# axis(1, at=Results$id, labels = F)
# text(Results$id, par("usr")[3]-0.003, labels = lablist, srt = 45,
#      adj = c(1.1,1.1), xpd = TRUE, cex=0.85)
# title(xlab="Distribution", line=4.2, cex.lab=1)
# dev.off()
# 
# 
# options(tikzDefaultEngine = 'xetex')
# tikz(file = "level_check_perm.tex", width=6, height=4)
# 
# par(xaxt="n")
# plot(Results$id, Results$power_perm, xlab="",
#      ylab="Actual Size", pch=19, ylim=c(0,0.1))
# abline(h=0.05, lty="dashed")
# lablist <- Results$dists
# axis(1, at=Results$id, labels = F)
# text(Results$id, par("usr")[3]-0.003, labels = lablist, srt = 45,
#      adj = c(1.1,1.1), xpd = TRUE, cex=0.85)
# title(xlab="Distribution", line=4.2, cex.lab=1)
# 
# dev.off()


