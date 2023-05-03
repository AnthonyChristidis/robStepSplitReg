# -------------------------------------------------
# Checking Input Data for robStepSplitReg Function
# -------------------------------------------------
DataCheck <- function(x, y, 
                      n_models,
                      model_saturation,
                      alpha, 
                      model_size,
                      robust,
                      compute_coef,
                      pense_alpha,
                      pense_cv_k,
                      pense_cv_repl,
                      cl){
  
  # Checking the input for the design matrix (x) and the response vector (y)
  if (all(!inherits(x, "matrix"), !inherits(x, "data.frame"))) {
    stop("x should belong to one of the following classes: matrix, data.frame")
  } else if (all(!inherits(y, "matrix"), all(!inherits(y, "numeric")))) {
    stop("y should belong to one of the following classes: matrix, numeric")
  } else if (any(anyNA(x), any(is.nan(x)), any(is.infinite(x)))) {
    stop("x should not have missing, infinite or nan values")
  } else if (any(anyNA(y), any(is.nan(y)), any(is.infinite(y)))) {
    stop("y should not have missing, infinite or nan values")
  } else {
    if(inherits(y, "matrix")) {
      if (ncol(y)>1){
        stop("y should be a vector")
      }
      # Force to vector if input was a matrix
      y <- as.numeric(y)
    }
    len_y <- length(y)
    if (len_y != nrow(x)) {
      stop("y and x should have the same number of rows")
    }
  }
  
  # Checking the input for the number of models
  if(!is.null(n_models)){
    if (!inherits(n_models, "numeric")) {
      stop("n_models should be numeric")
    } else if (any(!n_models == floor(n_models), n_models <= 0)) {
      stop("n_models should be a positive integer greater than 0")
    }
  }
  
  # Checking input for model criterion 
  if(!(model_saturation %in% c("fixed", "p-value")))
    stop("The shrinkage method must be one of: fixed or p-value.")
  
  # Check alpha value
  if(!inherits(alpha, "numeric")) {
    stop("alpha should be numeric.")
  } else if(!all(alpha <= 1, alpha > 0)) {
    stop("alpha should be between 0 and 1.")
  }
  
  # Checking the input for the size of the models
  if(!is.null(model_size)){
    if (!inherits(model_size, "numeric")) {
      stop("model_size should be numeric")
    } else if (any(!(model_size == floor(model_size)), model_size <= 1, model_size > nrow(x))) {
      stop("model_size should be an integer, greater than one and less than the sample size.")
    }
  }
  
  # Check robust parameter
  if(!(robust %in% c(TRUE, FALSE)))
    stop("robust should be TRUE or FALSE.")
  
  # Check parameter to determine whether coefficients are computed
  if(!(compute_coef %in% c(TRUE, FALSE)))
    stop("compute_coef should be TRUE or FALSE.")
  
  # Check PENSE alpha value
  if(!inherits(pense_alpha, "numeric")) {
    stop("pense_alpha should be numeric.")
  } else if(!all(pense_alpha <= 1, pense_alpha > 0)) {
    stop("pense_alpha should be between 0 and 1.")
  }
  
  # Check input for CV folds for PENSE
  if(!inherits(pense_cv_k, "numeric")) {
    stop("pense_cv_k should be numeric.")
  } else if(any(!pense_cv_k == floor(pense_cv_k), pense_cv_k <= 0)) {
    stop("pense_cv_k should be a positive integer.")
  }
  
  # Check input for number of times CV is done for PENSE
  if(!inherits(pense_cv_repl, "numeric")) {
    stop("pense_cv_repl should be numeric.")
  } else if(any(!pense_cv_repl == floor(pense_cv_repl), pense_cv_repl <= 0)) {
    stop("pense_cv_repl should be a positive integer.")
  }
  
  # Check input for number of clusters used for PENSE
  if(!is.null(cl)){
    if(!inherits(cl, "numeric")) {
      stop("cl should be numeric.")
    } else if(any(!cl == floor(cl), cl <= 0)) {
      stop("cl should be a positive integer.")
    }
  }
}