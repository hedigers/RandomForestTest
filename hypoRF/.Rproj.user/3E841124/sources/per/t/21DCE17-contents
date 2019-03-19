#' HypoRF; a Random Forest based Two Sample Test
#'
#' @description Performs a permutation two sample test based on the out-of-bag-error of random forest
#'
#' @param data1 An object of type "data.frame". The first sample.
#' @param data2 An object of type "data.frame". The second sample.
#' @param K A numeric value specifying the number of times the created label is permutated. For K = 1 a binomial test is carried out. The Default is K = 100.
#' @param normalapprox A logical value asking for the use of a normal approximation. Default is normalapprox = TRUE.
#' @param seed A numeric value for reproducibility.
#' @param ... Arguments to be passed to ranger
#' @return A list with elements
#' \itemize{
#' \item\code{pvalue:} The p-value of the test.
#' \item\code{obs:} The OOB-error in case of K>1 or the out-of-sample error in case of K=1 (binomial test).
#' \item\code{val:} The OOB-errors of the permutated random forests in case of K>1 (otherwise NULL).
#' \item\code{varest:} The estimated variance of the permutated random forest OOB-errors in case of K>1 (otherwise NULL).
#' \item\code{importance_ranking:} The variable importance measure, when importance == "impurity".
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
hypoRF <- function(data1, data2, K = 100, normalapprox = T, seed = NULL, ...) {

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

  if (nrow(data1)!=nrow(data2)) {
    warning("imbalanced dataset provied to random forest")
  }

  if (!is.numeric(K) || !is.finite(K) || K<1) {
    stop("input K not valid")
  }

  if (!is.logical(normalapprox)) {
    stop("input normalapprox not logical")
  }

  n1 <- nrow(data1)
  n2 <- nrow(data2)
  p  <- ncol(data1)

  if(K > 1){

    y <- factor(c(rep(1, n1) ,rep(0, n2)))

    d <- cbind(rbind(data1, data2), y)

    rf <- ranger::ranger(y~., data=d, ... )
    obs <- rf$prediction.error
    val <- sapply(1:K, function(s) { y<-sample(factor(c(rep(1,n1),rep(0,n2))));
    ds <- cbind(rbind(data1, data2), y);
    return(ranger::ranger(y~., data=ds, ... )$prediction.error)   })

    if (!normalapprox){
      pvalue = (sum(val < obs)+1) / (K+1)
    }else{
      varest <- stats::var(val)
      mean <- mean(val)
      pvalue <- stats::pnorm(obs, mean=mean, sd=sqrt(varest))
    }

    if(c("importance") %in% names(add.param)){
      if(add.param$importance == "impurity"){
        importance <- as.matrix(ranger::importance(rf))
        importance_ranking <- as.matrix(importance[order(importance[,1], decreasing = T),])
      }
    } else {
      importance_ranking <- NULL
    }

  } else if(K==1) {

    data1 <- as.matrix(data1)
    data2 <- as.matrix(data2)
    n <- n1 + n2
    m <- round(0.5 * min(n1, n2)) # maybe set 0.5 as a parameter ???????

    y_train <- as.factor(c(rep(1, n1-m), rep(0, n2-m)))
    y_test <- as.factor(c(rep(1, m), rep(0, m)))

    idtest <- sample(1:min(n1,n2), m, replace=F)

    data1train <- data1[-idtest, , drop=F]
    data2train <- data2[-idtest, , drop=F]
    data1test <- data1[idtest, , drop=F]
    data2test <- data2[idtest, , drop=F]

    rf <- ranger::ranger(Y~., data=data.frame(
      Y=y_train, X=rbind(data1train, data2train)), ...)

    if(c("importance") %in% names(add.param)){
      if(add.param$importance == "impurity"){
        importance <- as.matrix(ranger::importance(rf))
        importance_ranking <- as.matrix(importance[order(importance[,1], decreasing = T),])
      }
    } else {
      importance_ranking <- NULL
    }

    obs <- sum(as.numeric(as.character(y_test) != as.character(
      stats::predict(rf, data=data.frame(Y=y_test, X=rbind(data1test, data2test)))$prediction)))

    pvalue <- stats::pbinom(obs, size = 2*K*m, prob=0.5)

    val <- NULL
    varest <- NULL

  }

  return(list(pvalue = pvalue,
              obs = obs,
              val = val,
              varest = varest,
              importance_ranking = importance_ranking,
              call = match.call()))

}
