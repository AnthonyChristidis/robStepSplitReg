/*
 * ===========================================================
 * File Type: CPP
 * File Name: robStepReg_Main.cpp
 * Package Name: robStepSplitReg
 *
 * Created by Anthony-A. Christidis.
 * Copyright (c) Anthony-A. Christidis. All rights reserved.
 * ===========================================================
 */

// Header files included
#include "StepModel.hpp"
#include "StepModelFixed.hpp"

// [[Rcpp::export]]
std::vector<arma::uword> Robust_Stepwise(arma::mat& x, arma::vec& y,
                                         arma::mat& correlation_predictors, arma::vec& correlation_response,
                                         arma::uword& model_saturation,
                                         double& sig_level,
                                         arma::uword& model_size) {
  
  // Case with p-value
  if(model_saturation==0){
    
    // Create the stepwise model
    StepModel model(x, y,
                    correlation_predictors, correlation_response, 
                    sig_level);
    
    // Initialize the model through the constructor and add first predictor
    model.Find_First_Predictor(0);
    model.Add_Optimal_Predictor();

    // Variables for model updates
    double p_value = model.Get_P_Value();
    arma::uword n_pred = 0;
    if (p_value < sig_level)
      n_pred++;
    
    // Find new optimal predictor 
    model.Find_Optimal_Predictor();
    p_value = model.Get_P_Value();
    
    // Looping and adding predictors
    while (!model.Get_Full()) {
      
      // Add optimal predictor
      model.Add_Optimal_Predictor();
      
      // Update partial correlations for optimal model
      model.Find_Optimal_Predictor();
    }
    
    // Return model predictors
    return model.Get_Model_Predictors();
  }
  else{ // Case with fixed model size
      
    // Create the stepwise model
    StepModelFixed model(x, y,
                         correlation_predictors, correlation_response, 
                         model_size);
    
    // Initialize the model through the constructor and add first predictor
    model.Find_First_Predictor(0);
    model.Add_Optimal_Predictor();
    
    // Variables for model updates
    double p_value = model.Get_P_Value();
    arma::uword n_pred = 0;
    if (p_value < sig_level)
      n_pred++;
    
    // Find new optimal predictor 
    model.Find_Optimal_Predictor();
    p_value = model.Get_P_Value();
    
    // Looping and adding predictors
    while (!model.Get_Full()) {
      
      // Add optimal predictor
      model.Add_Optimal_Predictor();
      
      // Update partial correlations for optimal model
      model.Find_Optimal_Predictor();
    }
    
    // Return model predictors
    return model.Get_Model_Predictors();
  }
}