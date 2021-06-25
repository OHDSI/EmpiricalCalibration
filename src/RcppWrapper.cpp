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
NumericVector gridLlApproximation(NumericVector& x, const DataFrame& parameters) {
  NumericVector point = parameters["point"];
  NumericVector value = parameters["value"];
  int n = point.size();
  NumericVector result(x.size());
  for (int i = 0; i < x.size(); ++i) {
    if (x[i] < point[0]) {
      double slope = std::max(0.0, (value[1] - value[0]) / (point[1] - point[0]));
      result[i] = (x[i] - point[0]) * slope + value[0];
    } else if (x[i] >= point[n - 1]) {
      double slope = std::min(0.0, (value[n - 1] - value[n - 2]) / (point[n - 1] - point[n - 2]));
      result[i] = (x[i] - point[n - 1]) * slope + value[n - 1];
    } else {
      for (int j = 0; j < n; ++j) {
        if (point[j] >= x[i]) {
          double slope = (value[j] - value[j - 1]) / (point[j] - point[j - 1]);
          result[i] = (x[i] - point[j - 1]) * slope + value[j - 1];
          break;
        }
      }
    }
  }
  return result;
}

#endif // __RcppWrapper_cpp__
