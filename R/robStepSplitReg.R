#' 
#' @useDynLib robStepSplitReg
#' @importFrom Rcpp sourceCpp
#' 
#' @importFrom stats coef sd median mad cor
#'
#' @title Robust Stepwise Split Regularized Regression
#' 
#' @description \code{robStepSplitReg} performs robust stepwise split regularized regression.
#' 
#' @param x Design matrix.
#' @param y Response vector.
#' @param n_models Number of models into which the variables are split.
#' @param model_saturation Criterion to determine if a model is saturated. Must be one of "fixed" (default) or "p-value".
#' @param alpha P-value used to determine when the model is saturated
#' @param model_size Size of the models in the ensemble.
#' @param robust Argument to determine if robust measures of location, scale and correlation are used. Default is TRUE.
#' @param compute_coef Argument to determine if coefficients are computed (via adaptive PENSE) for each model. Default is FALSE.
#' @param pense_alpha Elastic net mixing parmeter for model shrinkage in adaptive PENSE. Default is 1/4.
#' @param pense_cv_k Number of folds for the cross-validation procedure in adaptive PENSE. Default is 5.
#' @param pense_cv_repl Number of replications of the cross-validation procedure. Default is 1.
#' @param cl Number of clusters used for the adaptive PENSE fit. If NULL (default), there is no parallelization.
#' 
#' @return An object of class robStepSplitReg.
#' 
#' @export
#' 
#' @author Anthony-Alexander Christidis, \email{anthony.christidis@stat.ubc.ca}
#' 
#' @seealso \code{\link{coef.robStepSplitReg}}, \code{\link{predict.robStepSplitReg}}
#' 
#' @examples 
#' # Required library
#' library(mvnfast)
#' 
#' # Simulation parameters
#' n <- 50
#' p <- 500
#' rho <- 0.8
#' p.active <- 100
#' snr <- 3
#' contamination.prop <- 0.2
#' 
#' # Setting the seed
#' set.seed(0)
#' 
#' # Simulation of beta vector
#' true.beta <- c(runif(p.active, 0, 5)*(-1)^rbinom(p.active, 1, 0.7), rep(0, p - p.active))
#' 
#' # Simulation of uncontaminated data 
#' sigma.mat <- matrix(0, nrow = p, ncol = p)
#' sigma.mat[1:p.active, 1:p.active] <- rho
#' diag(sigma.mat) <- 1
#' x <- mvnfast::rmvn(n, mu = rep(0, p), sigma = sigma.mat)
#' sigma <- as.numeric(sqrt(t(true.beta) %*% sigma.mat %*% true.beta)/sqrt(snr))
#' y <- x %*% true.beta + rnorm(n, 0, sigma)
#' 
#' # Contamination of data 
#' contamination_indices <- 1:floor(n*contamination.prop)
#' k_lev <- 2
#' k_slo <- 100
#' x_train <- x
#' y_train <- y
#' beta_cont <- true.beta
#' beta_cont[true.beta!=0] <- beta_cont[true.beta!=0]*(1 + k_slo)
#' beta_cont[true.beta==0] <- k_slo*max(abs(true.beta))
#' for(cont_id in contamination_indices){
#'   
#'   a <- runif(p, min = -1, max = 1)
#'   a <- a - as.numeric((1/p)*t(a) %*% rep(1, p))
#'   x_train[cont_id,] <- mvnfast::rmvn(1, rep(0, p), 0.1^2*diag(p)) + 
#'     k_lev * a / as.numeric(sqrt(t(a) %*% solve(sigma.mat) %*% a))
#'   y_train[cont_id] <- t(x_train[cont_id,]) %*% beta_cont
#' }
#' 
#' # Ensemble models
#' ensemble_fit <- robStepSplitReg(x_train, y_train,
#'                                 n_models = 2,
#'                                 model_saturation = c("fixed", "p-value")[1],
#'                                 alpha = 0.05, model_size = floor(n/5),
#'                                 robust = TRUE,
#'                                 compute_coef = TRUE,
#'                                 pense_alpha = 1/4, pense_cv_k = 5, pense_cv_repl = 1,
#'                                 cl = NULL)
#' 
#' # Ensemble coefficients
#' ensemble_coefs <- coef(ensemble_fit, group_index = 1:ensemble_fit$n_models)
#' sens_ensemble <- sum(which((ensemble_coefs[-1]!=0)) <= p.active)/p.active
#' spec_ensemble <- sum(which((ensemble_coefs[-1]!=0)) <= p.active)/sum(ensemble_coefs[-1]!=0)
#' 
#' # Simulation of test data
#' m <- 2e3
#' x_test <- mvnfast::rmvn(m, mu = rep(0, p), sigma = sigma.mat)
#' y_test <- x_test %*% true.beta + rnorm(m, 0, sigma)
#' 
#' # Prediction of test samples
#' ensemble_preds <- predict(ensemble_fit, newx = x_test, 
#'                           group_index = 1:ensemble_fit$n_models,
#'                           dynamic = FALSE)
#' mspe_ensemble <- mean((y_test - ensemble_preds)^2)/sigma^2
#' 
robStepSplitReg <- function(x, y, 
                            n_models = 1,
                            model_saturation = c("fixed", "p-value")[1],
                            alpha = 0.05, 
                            model_size = NULL,
                            robust = TRUE,
                            compute_coef = FALSE,
                            pense_alpha = 1/4,
                            pense_cv_k = 5,
                            pense_cv_repl = 1,
                            cl = NULL){
  
  # Data input check
  DataCheck(x, y, 
            n_models,
            model_saturation,
            alpha, 
            model_size,
            robust,
            compute_coef,
            pense_alpha,
            pense_cv_k,
            pense_cv_repl,
            cl)
  
  # Shuffle the data
  n <- nrow(x)
  p <- ncol(x)
  random.permutation <- sample(1:n, n)
  x <- x[random.permutation, ]
  y <- y[random.permutation]
  
  # Model size check
  if(model_saturation == "fixed"){
    
    if(is.null(model_size))
      model_size <- n - 1 else if(model_size >= n)
        stop("The size of the models cannot be equal or exceed n.")
  }
  
  # Standarization predictors and response
  # Computation of correlation for predictors and response
  if(robust){
    
    # Standardization of predictors and response
    x.std <- apply(x, 2, function(x) return((x - median(x))/mad(x)))
    y.std <- (y - median(y))/mad(y)
    
    # Computation of correlations for predictors and response
    xy.std <- cbind(x.std, y.std)
    est_xy <- cellWise::estLocScale(xy.std)
    xy_wrap <- cellWise::wrap(xy.std, est_xy$loc, est_xy$scale)$Xw
    rob.cor <- cor(xy_wrap)
    correlation.predictors <- rob.cor[-nrow(rob.cor),-ncol(rob.cor)]
    correlation.response <- rob.cor[-nrow(rob.cor), ncol(rob.cor)]
    
  } else{
    
    # Standardization of predictors and response
    x.std <- apply(x, 2, function(x) return((x - mean(x))/sd(x)))
    y.std <- (y - mean(y))/sd(y)
    
    # Computation of correlations for predictors and response
    xy.std <- cbind(x.std, y.std)
    CORxy <- cor(xy.std)
    correlation.predictors <- CORxy[-nrow(CORxy), -ncol(CORxy)]
    correlation.response <- CORxy[-nrow(CORxy), ncol(CORxy)]
  }
  
  # Model saturation criterion
  model_saturation_cpp <- ifelse(model_saturation == "fixed", 1, 0)
  
  # Invoking the CPP code for the algorithm
  if(n_models==1)
    selections <- Robust_Stepwise(x.std, y.std,
                                  correlation.predictors, correlation.response,
                                  model_saturation_cpp,
                                  alpha,
                                  model_size) else
                                selections <- Robust_Stepwise_Split(x.std, y.std,
                                                                    correlation.predictors, correlation.response, 
                                                                    model_saturation_cpp,
                                                                    alpha,
                                                                    model_size,
                                                                    n_models)
  
  # Adjusting predictors (incrementing)
  if(n_models == 1)
    selections <- list(selections)
  for(model_id in 1:n_models)
    selections[[model_id]] <- selections[[model_id]] + 1
  
  # Creating the list for the output
  output <- list(x = x, y = y,
    n_models = n_models,
    model_saturation = model_saturation,
    alpha = alpha,
    model_size = model_size,
    robust = robust,
    compute_coef = compute_coef,
    selections = selections,
    intercepts = list(), coefficients = list(),
    DDCx = cellWise::DDC(x, DDCpars = list(fastDDC = TRUE, nbngbrs = p-1)))

  # Computation of final coefficients
  if(compute_coef){
    
    for(model_id in 1:n_models){
      
      adapense_fit <- pense::adapense_cv(x[, output$selections[[model_id]]], y, 
                                         alpha = pense_alpha, cv_k = pense_cv_k, cv_repl = pense_cv_repl,
                                         cl = cl)
      output$intercepts[[model_id]] <- coef(adapense_fit)[1]
      output$coefficients[[model_id]] <- numeric(p)
      output$coefficients[[model_id]][output$selections[[model_id]]] <- coef(adapense_fit)[-1]
    }
  }
  
  # Create the object of class "stepSplitReg"
  class(output) <- append("robStepSplitReg", class(output))
  
  # Returning the output from the stepwise algorithm
  return(output)
}



