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
NumericVector gridLikelihood(NumericVector& x, const NumericVector& row, const NumericVector& gridX) {
  int n = gridX.size();
  NumericVector result(x.size());
  for (int i = 0; i < x.size(); ++i) {
    if (x[i] < gridX[0]) {
      double slope = std::max(0.0, (row[1] - row[0]) / (gridX[1] - gridX[0]));
      result[i] = (x[i] - gridX[0]) * slope + row[0];
    } else if (x[i] >= gridX[n - 1]) {
      double slope = std::min(0.0, (row[n - 1] - row[n - 2]) / (gridX[n - 1] - gridX[n - 2]));
      result[i] = (x[i] - gridX[n - 1]) * slope + row[n - 1];
    } else {
      for (int j = 0; j < n; ++j) {
        if (gridX[j] >= x[i]) {
          double slope = (row[j] - row[j - 1]) / (gridX[j] - gridX[j - 1]);
          result[i] = (x[i] - gridX[j - 1]) * slope + row[j - 1];
          break;
        }
      }
    }
  }
  return result;
}

#endif // __RcppWrapper_cpp__
