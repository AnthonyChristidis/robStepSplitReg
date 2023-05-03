/*
 * ===========================================================
 * File Type: HPP
 * File Name: Generate_Predictors_List.hpp
 * Package Name: robStepSplitReg
 *
 * Created by Anthony-A. Christidis.
 * Copyright (c) Anthony-A. Christidis. All rights reserved.
 * ===========================================================
 */

// Header files included
#include "StepModel.hpp"
#include "StepModelFixed.hpp"


// Return a list of vectors with the variables in each model
Rcpp::List Generate_Predictors_List(std::vector<StepModel*> final_models, arma::uword& n_models) {
  
  Rcpp::List final_predictors_list(n_models);
  for (arma::uword m = 0; m < n_models; m++)
    final_predictors_list[m] = final_models[m]->Get_Model_Predictors();
  
  return final_predictors_list;
}

// Return a list of vectors with the variables in each model
Rcpp::List Generate_Predictors_List_Fixed(std::vector<StepModelFixed*> final_models, arma::uword& n_models) {
  
  Rcpp::List final_predictors_list(n_models);
  for (arma::uword m = 0; m < n_models; m++)
    final_predictors_list[m] = final_models[m]->Get_Model_Predictors();
  
  return final_predictors_list;
}