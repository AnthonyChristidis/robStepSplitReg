#' Robust Multi-Model Stepwise Variable Selection
#' 
#' @description Internal function that performs robust stepwise variable selection 
#' for building disjoint ensemble models. Uses robust correlation estimates to 
#' select predictor subsets while ensuring variables are not shared across models.
#' 
#' @param Rx Robust correlation matrix of predictors (p x p matrix).
#' @param Ry Robust correlation vector between response and predictors (p x 1 vector).
#' @param n_models Number of models to build. Default is 1.
#' @param model_saturation Criterion to determine model saturation. Either "p-value" 
#'   or "fixed". Default is "p-value".
#' @param alpha P-value threshold for determining model saturation when 
#'   model_saturation = "p-value". Default is 0.05.
#' @param model_size Maximum number of variables per model when model_saturation = 
#'   "fixed". Default is NULL.
#' @param n Sample size used for computing test statistics. Default is nrow(Rx).
#' 
#' @return If n_models = 1, returns a vector of 0-indexed selected variable indices.
#'   If n_models > 1, returns a list where each element contains the 0-indexed 
#'   selected variable indices for the corresponding model.
#' 
#' @author Anthony-Alexander Christidis, \email{anthony.christidis@stat.ubc.ca}
#' 
#' @keywords internal
#' 
robustStepwiseSplit <- function(Rx, Ry, n_models = 1, model_saturation = "p-value", 
                                alpha = 0.05, model_size = NULL, n = nrow(Rx)) {
    
    p <- length(Ry)
    
    # Pre-allocate all data structures
    model_predictors <- vector("list", n_models)
    model_saturated <- logical(n_models)
    model_rss <- rep(n, n_models)
    candidate_pool <- 1:p
    
    # Use matrices instead of lists for better performance
    P_y <- matrix(n * Ry, nrow = p, ncol = n_models)  # Each column is a model
    P_xx <- array(n * Rx, dim = c(p, p, n_models))    # Each slice is a model
    active_vars <- matrix(TRUE, nrow = p, ncol = n_models)  # Track which vars are available
    
    for(g in 1:n_models) {
        model_predictors[[g]] <- integer(0)
    }
    
    # Pre-allocate working vectors
    test_stats <- numeric(p)
    pvalues <- numeric(n_models)
    rss_decreases <- numeric(n_models)
    best_vars <- integer(n_models)
    
    # Main algorithm loop
    iteration <- 0L
    max_iterations <- 2L * p
    
    while(iteration < max_iterations && any(!model_saturated) && length(candidate_pool) > 0L) {
        iteration <- iteration + 1L
        
        # Reset working vectors
        pvalues[] <- 1.0
        best_vars[] <- 0L
        rss_decreases[] <- 0.0
        
        # Step 1: For each unsaturated model, find optimal next variable
        for(g in 1:n_models) {
            if(model_saturated[g]) next
            
            k <- length(model_predictors[[g]])
            
            # Find available candidates for this model
            available_mask <- active_vars[, g] & (1:p %in% candidate_pool)
            if(!any(available_mask)) {
                model_saturated[g] <- TRUE
                next
            }
            
            if(k == 0L) {
                # First variable: maximize |P_yj|/n
                test_stats[] <- 0.0
                test_stats[available_mask] <- abs(P_y[available_mask, g]) / n
                best_idx <- which.max(test_stats)
                best_vars[g] <- best_idx
                rss_decreases[g] <- (P_y[best_idx, g])^2 / n
            } else {
                # k >= 1: maximize |P_yj / sqrt(P_jj)|
                test_stats[] <- 0.0
                diag_vals <- pmax(P_xx[cbind(1:p, 1:p, g)], 1e-12)
                test_stats[available_mask] <- abs(P_y[available_mask, g]) / sqrt(diag_vals[available_mask])
                best_idx <- which.max(test_stats)
                best_vars[g] <- best_idx
                rss_decreases[g] <- (P_y[best_idx, g])^2 / P_xx[best_idx, best_idx, g]
            }
            
            # Bound RSS decrease
            rss_decreases[g] <- max(0.0, min(rss_decreases[g], model_rss[g] * 0.99))
            
            # Compute F-test p-value inline
            old_rss <- model_rss[g]
            new_rss <- old_rss - rss_decreases[g]
            
            if(k + 1L >= n - 1L || new_rss <= 0.0 || new_rss >= old_rss) {
                pvalues[g] <- 1.0
            } else {
                f_stat <- max(0.0, ((old_rss - new_rss) / 1.0) / (new_rss / (n - k - 2L)))
                pvalues[g] <- pf(f_stat, df1 = 1L, df2 = n - k - 2L, lower.tail = FALSE)
            }
            
            # Check saturation criteria
            if((model_saturation == "p-value" && pvalues[g] >= alpha) ||
               (model_saturation == "fixed" && !is.null(model_size) && k >= model_size) ||
               k >= n - 2L) {
                model_saturated[g] <- TRUE
            }
        }
        
        # Step 2: Find best model
        unsaturated_models <- which(!model_saturated & best_vars > 0L)
        if(length(unsaturated_models) == 0L) break
        
        best_model <- unsaturated_models[which.min(pvalues[unsaturated_models])]
        
        # Step 3: Update if significant
        if(pvalues[best_model] < alpha || model_saturation == "fixed") {
            
            selected_var <- best_vars[best_model]
            if(selected_var == 0L) break
            
            # Add variable to selected model
            model_predictors[[best_model]] <- c(model_predictors[[best_model]], selected_var)
            model_rss[best_model] <- model_rss[best_model] - rss_decreases[best_model]
            
            # Fast P-term updates for selected model using vectorized operations
            if(any(active_vars[, best_model])) {
                P_y_jk <- P_y[selected_var, best_model]
                P_jk_jk <- max(P_xx[selected_var, selected_var, best_model], 1e-12)
                
                # Vectorized update for P_y terms
                remaining_mask <- active_vars[, best_model]
                remaining_mask[selected_var] <- FALSE
                
                if(any(remaining_mask)) {
                    P_j_jk <- P_xx[remaining_mask, selected_var, best_model]
                    P_y[remaining_mask, best_model] <- P_y[remaining_mask, best_model] - (P_j_jk / P_jk_jk) * P_y_jk
                    
                    # Vectorized update for P_xx terms
                    remaining_idx <- which(remaining_mask)
                    P_i_jk <- P_xx[remaining_idx, selected_var, best_model]
                    
                    # Use outer product for efficient matrix update
                    P_xx[remaining_idx, remaining_idx, best_model] <- 
                        P_xx[remaining_idx, remaining_idx, best_model] - outer(P_i_jk, P_i_jk) / P_jk_jk
                }
                
                # Mark variable as inactive for this model
                active_vars[selected_var, best_model] <- FALSE
            }
            
            # Remove from candidate pool
            candidate_pool <- candidate_pool[candidate_pool != selected_var]
            
            # Mark variable as inactive for all other models (no recursive update needed)
            active_vars[selected_var, -best_model] <- FALSE
            
            # Check size-based saturation
            max_size <- if(model_saturation == "fixed" && !is.null(model_size)) model_size else n - 2L
            if(length(model_predictors[[best_model]]) >= max_size) {
                model_saturated[best_model] <- TRUE
            }
        } else {
            break
        }
    }
    
    # Return 0-indexed results
    return(model_predictors)
}