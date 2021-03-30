#' HypoRF; a Random Forest based Two Sample Test
#'
#' @description Performs a permutation two sample test based on the out-of-bag-error of random forest
#'
#' @param data1 An object of type "data.frame". The first sample.
#' @param data2 An object of type "data.frame". The second sample.
#' @param K A numeric value specifying the number of times the created label is permuted. For K = 1 a binomial test is carried out. The Default is K = 100.
#' @param statistic A character value specifying the statistic for permutation testing. Two options available
#' \itemize{
#' \item\code{PerClassOOB} Sum of OOB per class errors.
#' \item\code{OverallOOB} OOB-error.
#' }. Default is statistic = "PerClassOOB".
#' @param normalapprox A logical value asking for the use of a normal approximation. Default is normalapprox = FALSE.
#' @param seed A numeric value for reproducibility.
#' @param alpha The level of the test. Default is alpha = 0.05.
#' @param ... Arguments to be passed to ranger
#' @return A list with elements
#' \itemize{
#' \item\code{pvalue:} The p-value of the test.
#' \item\code{obs:} The OOB-statistic in case of K>1 or the out-of-sample error in case of K=1 (binomial test).
#' \item\code{val:} The OOB-statistic of the permuted random forests in case of K>1 (otherwise NULL).
#' \item\code{varest:} The estimated variance of the permuted random forest OOB-statistic in case of K>1 (otherwise NULL).
#' \item\code{statistic:} The used OOB-statistic
#' \item\code{importance_ranking:} The variable importance measure, when importance == "impurity".
#' \item\code{cutoff:} The quantile of the importance distribution at level alpha.
#' \item\code{call:} Call to the function.}
#' @seealso \code{\link{ranger}}
#' @examples
#' # Using the default testing procedure (permutation test)
#' x1 <- data.frame(x=rt(100, df=1.5))
#' x2 <- data.frame(x=rnorm(100))
#' hypoRF(x1, x2)
#' # Using the exact binomial test
#' hypoRF(x1, x2, K=1)
#' @export
hypoRF <- function(data1, data2, K = 100,
                   statistic = "PerClassOOB",
                   normalapprox = F, seed = NULL, alpha=0.05, ...) {

  if (!is.null(seed) && (!is.numeric(seed) || !is.finite(seed) || seed < 0)) {
    stop("invalid seed provided")
  }

  if(!is.null(seed)) {
    set.seed(seed)
  }

  add.param <- list(...)

  if (c("probability") %in% names(add.param)) {
    stop("probability parameter not meaningfull for testing")
  }

  if (!is.data.frame(data1) || !is.data.frame(data2)){
    stop("data1 or data2 are not of class data.frame")
  }

  if (ncol(data1)!=ncol(data2)) {
    stop("amount of columns in data1 and data2 not identical")
  }

  if (!identical(colnames(data1),colnames(data2))) {
    stop("colnames of the datasets do not match")
  }

  if (c("y")%in%colnames(data1)) {
    stop("y is a reserved character for internal use")
  }

  if (!is.numeric(K) || !is.finite(K) || K<1) {
    stop("input K not valid")
  }

  if (is.na(statistic) || !(statistic %in% c("PerClassOOB", "OverallOOB"))){
    stop("invalid statistic")
  }

  if (!is.logical(normalapprox)) {
    stop("input normalapprox not logical")
  }

  if ((alpha > 1) || (alpha < 0)) {
    stop("invalid alpha")
  }

  n1 <- nrow(data1) # number of class 0
  n2 <- nrow(data2) # number of class 1
  p  <- ncol(data1)

  if(K > 1){

    y <- factor(c(rep(0, n1) ,rep(1, n2)))

    d <- cbind(rbind(data1, data2), y)

    rf <- ranger::ranger(y~., data=d, probability=T, ...)

    if(statistic == "OverallOOB"){
      obs <- rf$prediction.error
    } else if(statistic == "PerClassOOB") {
      pred_oob <- as.numeric(rf$predictions[,2] > (n2/(n1+n2)))
      tmp <- cbind(d,pred_oob)[,c("y","pred_oob")]

      L0<-sum(tmp[y==0,"pred_oob"])/length(tmp[y==0,"pred_oob"])
      L1<-(length(tmp[y==1,"pred_oob"])-sum(tmp[y==1,"pred_oob"]))/length(tmp[y==1,"pred_oob"])
      obs <-  L0 + L1
    }

    val <- importancedistribution <- numeric()
    for(i in 1:K){
      y <- sample(factor(c(rep(0,n1),rep(1,n2))));
      ds <- cbind(rbind(data1, data2), y);

      rf_i <- ranger::ranger(y~., data=ds, probability = T, ...)
      if(statistic == "OverallOOB"){
        val[i] <- rf_i$prediction.error
      } else if(statistic == "PerClassOOB"){
        pred_oob <- as.numeric(rf_i$predictions[,2] > (n2/(n1+n2)))
        tmpi <- cbind(d,pred_oob)[,c("y","pred_oob")]
        L0i<-sum(tmpi[y==0,"pred_oob"])/length(tmpi[y==0,"pred_oob"])
        L1i<-(length(tmpi[y==1,"pred_oob"])-sum(tmpi[y==1,"pred_oob"]))/length(tmpi[y==1,"pred_oob"])
        val[i] <-  L0i + L1i
      }

      if(!is.null(rf_i$variable.importance)){
        importancedistribution[i] <- max(rf_i$variable.importance)
      }else{
        importancedistribution <- NULL
      }
    }

    if (!normalapprox){
      pvalue = (sum(val < obs)+1) / (K+1)
      varest <- NULL
    }else{
      varest <- stats::var(val)
      mean <- mean(val)
      pvalue <- stats::pnorm(obs, mean=mean, sd=sqrt(varest))
    }

    if(c("importance") %in% names(add.param)){
      if(add.param$importance == "impurity"){
        importance <- as.matrix(ranger::importance(rf))
        importance_ranking <- as.matrix(importance)
        cutoff<-quantile(importancedistribution , 1-alpha)

      }
    } else {
      importance_ranking <- NULL
      cutoff <- NULL
    }

  } else if(K==1) {

    data1 <- as.matrix(data1)
    data2 <- as.matrix(data2)
    n <- n1 + n2
    m1 <- round(0.5*n1) # number of class 0
    m2 <- round(0.5*n2) # number of class 1

    y_train <- as.factor(c(rep(1, n1-m1), rep(0, n2-m2)))
    y_test <- as.factor(c(rep(1, m1), rep(0, m2)))

    idtest1 <- sample(1:n1, m1)
    idtest2 <- sample(1:n2, m2)

    data1train <- data1[-idtest1, , drop=F]
    data2train <- data2[-idtest2, , drop=F]
    data1test <- data1[idtest1, , drop=F]
    data2test <- data2[idtest2, , drop=F]

    rf <- ranger::ranger(Y~., data=data.frame(
      Y=y_train, X=rbind(data1train, data2train)), probability = T,...)

    if(c("importance") %in% names(add.param)){
      if(add.param$importance == "impurity"){
        importance <- as.matrix(ranger::importance(rf))
        importance_ranking <- as.matrix(importance[order(importance[,1], decreasing = T),])
      }
    } else {
      importance_ranking <- NULL
    }

    if(statistic == "OverallOOB"){

      pred_prob<-stats::predict(rf, data=data.frame(Y=y_test, X=rbind(data1test, data2test)))$prediction
      pred_class<-as.numeric(pred_prob[,2] > (n2-m2)/((n1-m1) + (n2-m2)))

      obs <- sum(as.numeric(as.character(y_test) != as.character(pred_class)))

      pvalue <- stats::pbinom(obs, size = (m1+m2), prob=0.5)
    } else if(statistic == "PerClassOOB"){

      pred_prob<-stats::predict(rf, data=data.frame(Y=y_test, X=rbind(data1test, data2test)))$prediction
      pred_class<-as.numeric(pred_prob[,2] > (n2-m2)/((n1-m1) + (n2-m2)))

      tmp<-data.frame(y_test=y_test, pred_class=pred_class)

      obs1<-sum(tmp[y_test==0,"pred_class"])/length(tmp[y_test==0,"pred_class"])
      obs2<-(length(tmp[y_test==1,"pred_class"])-sum(tmp[y_test==1,"pred_class"]))/length(tmp[y_test==1,"pred_class"])

      if(sqrt(obs1*(1-obs1)+obs2*(1-obs2))==0){
        if(obs1+obs2==1){
          pvalue <- 1
        } else {
          pvalue <- 0
        }

      } else {
        obs <- (obs1 + obs2 -1)/sqrt(obs1*(1-obs1)/m1+obs2*(1-obs2)/m2)
        pvalue <- pnorm(obs,mean=0,sd=1)
      }

    }

    val <- NULL
    varest <- NULL
    importancedistribution <- NULL
    cutoff <- NULL

  }

  return(list(pvalue = pvalue,
              obs = obs,
              val = val,
              varest = varest,
              statistic = statistic,
              importance_ranking = importance_ranking,
              importancedistribution=importancedistribution,
              cutoff=cutoff,
              call = match.call()))

}
