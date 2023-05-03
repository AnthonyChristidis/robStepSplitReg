/*
 * ===========================================================
 * File Type: HPP
 * File Name: StepModelFixed.hpp
 * Package Name: robStepSplitReg
 *
 * Created by Anthony-A. Christidis.
 * Copyright (c) Anthony-A. Christidis. All rights reserved.
 * ===========================================================
 */

#ifndef StepModelFixed_hpp
#define StepModelFixed_hpp

// Libraries included
#include <RcppArmadillo.h>
#include <vector>

class StepModelFixed {
  
private:
  
  // Variables supplied by the user
  arma::mat x;
  arma::vec y;
  arma::mat correlation_predictors;
  arma::vec correlation_response;
  arma::uword model_size;
  
  // Variables created inside class
  arma::uword n;
  arma::uword p;
  std::vector<arma::uword> model_predictors, available_predictors;
  std::vector<arma::uword>::iterator begin_iterator, end_iterator;
  arma::vec partial_correlations;
  arma::uword optimal_predictor;
  arma::mat z_old, z_new;
  double beta_y_optimal;
  arma::vec residuals_old, residuals_new;
  double rss_old, rss_new;
  double F_value;
  double p_value;
  bool model_full;
  
public:
  
  // (+) Model Constructor
  
  StepModelFixed(arma::mat& x, arma::vec& y,
                 arma::mat& correlation_predictors, arma::vec& correlation_response,
                 arma::uword& model_size);
  
  // (+) Functions that update the current state of the model  
  
  // Functions to potentially add a predictor
  void Find_First_Predictor(arma::uword index);
  void Find_Optimal_Predictor();
  void Add_Optimal_Predictor();
  
  // Functions to add or remove a predictor
  void Add_Model_Predictor(arma::uword& predictor);
  void Remove_Available_Predictor(arma::uword predictor);
  void Remove_Available_Predictor_Update(arma::uword predictor);
  
  // Function to update z matrix
  void Update_Z_Matrix();
  
  // Functions to update model status
  void Update_Partial_Correlations();
  void Update_Optimal_Predictor();
  void Update_Beta_Y_Optimal();
  void Update_Residuals();
  void Update_RSS();
  void Update_F_Value();
  void Update_P_Value();
  void Check_Full();
  
  // (+) Functions that return the state of the model
  bool Get_Full();
  double Get_F_Value();
  double Get_P_Value();
  arma::uword Get_Optimal_Predictor();
  std::vector<arma::uword> Get_Model_Predictors();
};

#endif // StepModelFixed_hpp
