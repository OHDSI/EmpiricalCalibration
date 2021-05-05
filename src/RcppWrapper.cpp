/*
 * @file RcppWrapper.cpp
 *
 * This file is part of EmpiricalCalibration
 *
 * Copyright 2021 Observational Health Data Sciences and Informatics
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef __RcppWrapper_cpp__
#define __RcppWrapper_cpp__

#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
double gridOneLikelihood(double x, const NumericVector& gridX, const NumericVector& gridY) {
  int n = gridX.size();
  if (x < gridX[0]) {
    double slope = std::max(0.0, (gridY[1] - gridY[0]) / (gridX[1] - gridX[0]));
    return (x - gridX[0]) * slope + gridY[0];
  } 
  if (x >= gridX[n - 1]) {
    double slope = std::min(0.0, (gridY[n - 1] - gridY[n - 2]) / (gridX[n - 1] - gridX[n - 2]));
    return (x - gridX[n - 1]) * slope + gridY[n - 1];
  } 
  for(int i = 0; i < n; ++i) {
    if (gridX[i] >= x) {
      double slope = (gridY[i] - gridY[i - 1]) / (gridX[i] - gridX[i - 1]);
      return (x - gridX[i - 1]) * slope + gridY[i - 1];
    }
  }
  return 0;
}

#endif // __RcppWrapper_cpp__
