/*
 * ===========================================================
 * File Type: CPP
 * File Name: StepModel.cpp
 * Package Name: robStepSplitReg
 *
 * Created by Anthony-A. Christidis.
 * Copyright (c) Anthony-A. Christidis. All rights reserved.
 * ===========================================================
 */

// Header files included
#include "StepModel.hpp"

// (+) Model Constructor

StepModel::StepModel(arma::mat& x, arma::vec& y,
                     arma::mat& correlation_predictors, arma::vec& correlation_response,
                     double& sig_level) :
  x(x), y(y),
  correlation_predictors(correlation_predictors), correlation_response(correlation_response),
  sig_level(sig_level) {
  
  // Initialize dimension of data
  n = x.n_rows;
  p = x.n_cols;
  
  // Initialize available predictors
  for (arma::uword pred_id = 0; pred_id < p; pred_id++)
    available_predictors.push_back(pred_id);
  
  // Initialize partial correlations
  partial_correlations = correlation_response;
  
  // Initialize z matrices
  z_old = z_new = x;
  
  // Initialize residuals
  residuals_old = residuals_new = y;
  rss_old = rss_new = arma::as_scalar(y.t()*y);
  
  // Initialize model saturation
  model_full = false;
}

// (+) Functions that update the current state of the model

// Functions for first predictor
void StepModel::Find_First_Predictor(arma::uword index) {
  
  arma::uvec correlation_index = arma::sort_index(arma::abs(correlation_response), "descend");
  optimal_predictor = correlation_index(index);
  beta_y_optimal = correlation_response(optimal_predictor);
  residuals_new = y - beta_y_optimal*x.col(optimal_predictor);
  Update_RSS();
  Update_F_Value();
  Update_P_Value();
  // Check_Full();
}

// Function for finding optimal predictor (beyond first two predictors)
void StepModel::Find_Optimal_Predictor() {
  
  Update_Z_Matrix();
  Update_Partial_Correlations();
  Update_Optimal_Predictor();
  Update_Beta_Y_Optimal();
  Update_Residuals();
  Update_RSS();
  Update_F_Value();
  Update_P_Value();
  Check_Full();
}

// Function to add optimal predictor to model
void StepModel::Add_Optimal_Predictor() {
  
  if ((!Get_Full()) && (model_predictors.size()<n)) {
    Add_Model_Predictor(optimal_predictor);
    Remove_Available_Predictor(optimal_predictor);
    residuals_old = residuals_new;
    rss_old = rss_new;
    z_old = z_new;
  }
  else
    model_full = true;
}

// Functions to add or remove a predictor
void StepModel::Add_Model_Predictor(arma::uword& predictor) {
  model_predictors.push_back(predictor);
}
void StepModel::Remove_Available_Predictor(arma::uword predictor) {
  
  std::vector<arma::uword>::iterator drop_position = std::find(available_predictors.begin(), available_predictors.end(), predictor);
  if (drop_position != available_predictors.end())
    available_predictors.erase(drop_position);
  partial_correlations(predictor) = 0;
}
void StepModel::Remove_Available_Predictor_Update(arma::uword predictor) {
  
  std::vector<arma::uword>::iterator drop_position = std::find(available_predictors.begin(), available_predictors.end(), predictor);
  if (drop_position != available_predictors.end())
    available_predictors.erase(drop_position);
  partial_correlations(predictor) = 0;
  Update_Optimal_Predictor();
  Update_Beta_Y_Optimal();
  Update_Residuals();
  Update_RSS();
  Update_F_Value();
  Update_P_Value();
  Check_Full();
}

// Function to update z matrix
void StepModel::Update_Z_Matrix() {
  
  if (model_predictors.size() == 1) {
    for (arma::uword pred_id = 0; pred_id < p; pred_id++)
      z_new.col(pred_id) = z_old.col(pred_id) - arma::as_scalar(correlation_predictors(pred_id, optimal_predictor)) * z_old.col(optimal_predictor);
  }
  else {
    for (arma::uword pred_id = 0; pred_id < p; pred_id++)
      z_new.col(pred_id) = z_old.col(pred_id) - (arma::as_scalar(z_old.col(pred_id).t() * z_old.col(optimal_predictor)) / arma::as_scalar(z_old.col(optimal_predictor).t() * z_old.col(optimal_predictor))) * z_old.col(optimal_predictor);
  }
}

// Functions to update model status
void StepModel::Update_Partial_Correlations() {
  
  begin_iterator = available_predictors.begin();
  end_iterator = begin_iterator + available_predictors.size();
  for (auto pred_id = begin_iterator; pred_id != end_iterator; pred_id++) {
    partial_correlations(*pred_id) =
      arma::as_scalar(z_new.col(*pred_id).t() * y) / arma::as_scalar(z_new.col(*pred_id).t() * z_new.col(*pred_id)) / std::sqrt(n);
  }
}
void StepModel::Update_Optimal_Predictor() {
  
  optimal_predictor = arma::abs(partial_correlations).index_max();
}
void StepModel::Update_Beta_Y_Optimal() {
  
  beta_y_optimal = arma::as_scalar((z_new.col(optimal_predictor).t() * y)) / arma::as_scalar((z_new.col(optimal_predictor).t() * z_new.col(optimal_predictor)));
}
void StepModel::Update_Residuals() {
  
  residuals_new = residuals_old - beta_y_optimal * z_new.col(optimal_predictor);
}
void StepModel::Update_RSS() {
  
  rss_new = arma::as_scalar(residuals_new.t() * residuals_new);
}
void StepModel::Update_F_Value() {
  
  F_value = (rss_old - rss_new) / rss_new * (n - model_predictors.size() - 1);
}
void StepModel::Update_P_Value() {
  
  p_value = R::pf(F_value, 1, n - model_predictors.size() - 1, 0, 0);
}

void StepModel::Check_Full() {
  
  if (p_value >= sig_level)
    model_full = true;
}

// (+) Functions that return the state of the model
bool StepModel::Get_Full() {
  return model_full;
}

double StepModel::Get_F_Value() {
  return F_value;
}

double StepModel::Get_P_Value() {
  return p_value;
}

arma::uword StepModel::Get_Optimal_Predictor() {
  return optimal_predictor;
}

std::vector<arma::uword> StepModel::Get_Model_Predictors() {
  return model_predictors;
}
