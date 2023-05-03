/*
 * ===========================================================
 * File Type: CPP
 * File Name: robStepSplitReg_Main.cpp
 * Package Name: robStepSplitReg
 *
 * Created by Anthony-A. Christidis.
 * Copyright (c) Anthony-A. Christidis. All rights reserved.
 * ===========================================================
 */

// Header files included
#include "StepModel.hpp"
#include "StepModelFixed.hpp"
#include "Generate_Predictors_List.hpp"

// [[Rcpp::export]]
Rcpp::List Robust_Stepwise_Split(arma::mat& x, arma::vec& y,
                                 arma::mat& correlation_predictors, arma::vec& correlation_response,
                                 arma::uword& model_saturation,
                                 double& sig_level,
                                 arma::uword& model_size,
                                 arma::uword& n_models){
  // Case with p-value
  if(model_saturation==0){
    
    // Create the memory for the models (through dynamic allocation)
    std::vector<StepModel*> models;
    
    // Initialize the models through the constructors and add first predictor
    for (arma::uword m = 0; m < n_models; m++) {
      
      models.push_back(new StepModel(x, y, correlation_predictors, correlation_response, sig_level));
      models[m]->Find_First_Predictor(m);
      models[m]->Add_Optimal_Predictor();
    }
    
    // Remove initial predictors already used
    for  (arma::uword m = 0; m < n_models; m++)
      for (arma::uword r = 0; r < n_models; r++){
        if(r != m)
          models[m]->Remove_Available_Predictor(models[r]->Get_Optimal_Predictor());
      }
        
    // Variables for model updates
    arma::vec p_values = arma::ones(n_models);
    arma::uword optimal_model;
    arma::uword n_pred = 0;
    for (arma::uword m = 0; m < n_models; m++) {
      
      p_values(m) = models[m]->Get_P_Value();
      if (!(models[m]->Get_Full()))
        n_pred++;
    }
    
    // Find optimal predictor for unsaturated models
    for (arma::uword m = 0; m < n_models; m++){
      if (!models[m]->Get_Full()) {
        models[m]->Find_Optimal_Predictor();
        p_values[m] = models[m]->Get_P_Value();
      }
    }
    
    // Looping and adding predictors
    while (n_pred < x.n_cols) {
      
      // Find optimal model for update
      optimal_model = p_values.index_min();
      
      // Add optimal predictor
      if (models[optimal_model]->Get_P_Value() < sig_level) {
        models[optimal_model]->Add_Optimal_Predictor();
        n_pred++;
      }
      else
        break; // Update for optimal model is not statistically significant
      
      // Remove optimal predictor for non-optimal models
      for (arma::uword m = 0; m < n_models; m++){
        if ((!models[m]->Get_Full()) && (m != optimal_model))
          models[m]->Remove_Available_Predictor_Update(models[optimal_model]->Get_Optimal_Predictor());
      }
      
      // Update partial correlations for optimal model
      models[optimal_model]->Find_Optimal_Predictor();
      
      // Update p-values
      for (arma::uword m = 0; m < n_models; m++){
        if (!models[m]->Get_Full()) 
          p_values[m] = models[m]->Get_P_Value();
      }
    }
    
    // List with variables in each model
    Rcpp::List final_predictors_list = Generate_Predictors_List(models, n_models);
    
    // Delete the models
    for (arma::uword m = 0; m < n_models; m++)
      delete(models[m]);
    
    return final_predictors_list;
  } 
  else{ // Case with fixed model size
    
    // Create the memory for the models (through dynamic allocation)
    std::vector<StepModelFixed*> models;
    
    // Initialize the models through the constructors and add first predictor
    for (arma::uword m = 0; m < n_models; m++) {
      
      models.push_back(new StepModelFixed(x, y, correlation_predictors, correlation_response, model_size));
      models[m]->Find_First_Predictor(m);
      models[m]->Add_Optimal_Predictor();
    }
    
    // Remove initial predictors already used
    for  (arma::uword m = 0; m < n_models; m++)
      for (arma::uword r = 0; r < n_models; r++)
        if(r != m)
          models[m]->Remove_Available_Predictor(models[r]->Get_Optimal_Predictor());
        
    // Variables for model updates
    arma::vec p_values = arma::ones(n_models);
    arma::uword optimal_model;
    arma::uword full_models = 0; 
    arma::uword n_pred = 0;
    for (arma::uword m = 0; m < n_models; m++) {
      
      if (!(models[m]->Get_Full()))
        n_pred++;
      else
        full_models++;
    }
    
    // Find optimal predictor for unsaturated models
    for (arma::uword m = 0; m < n_models; m++){
      if (!models[m]->Get_Full()) {
        models[m]->Find_Optimal_Predictor();
        p_values[m] = models[m]->Get_P_Value();
      }
    }

    // Looping and adding predictors
    while ((n_pred < x.n_cols) && (full_models < n_models)) {
      
      // Find optimal model for update
      optimal_model = p_values.index_min();
      
      // Add optimal predictor
      if (!(models[optimal_model]->Get_Full())) {
        models[optimal_model]->Add_Optimal_Predictor();
        n_pred++;
        
        // Remove optimal predictor for non-optimal models
        for (arma::uword m = 0; m < n_models; m++){
          if ((!models[m]->Get_Full()) && (m != optimal_model))
            models[m]->Remove_Available_Predictor_Update(models[optimal_model]->Get_Optimal_Predictor());
        }
        
        // Update partial correlations for optimal model 
        models[optimal_model]->Find_Optimal_Predictor();
      } 
      else{
        full_models++;
        p_values[optimal_model] = 2; 
      }
    }
    
    // List with variables in each model
    Rcpp::List final_predictors_list = Generate_Predictors_List_Fixed(models, n_models);
    
    // Delete the models
    for (arma::uword m = 0; m < n_models; m++)
      delete(models[m]);
    
    return final_predictors_list;
  }
}