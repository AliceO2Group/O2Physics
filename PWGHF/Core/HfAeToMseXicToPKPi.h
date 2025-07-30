// Copyright 2020-2022 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file HfAeToMseXicToPKPi.h
/// \brief Class to compute the mse for the Autoencoder for Xic+ → p K- π+ analysis selections
/// \author Maria Teresa Camerlingo

#ifndef PWGHF_CORE_HFAETOMSEXICTOPKPI_H_
#define PWGHF_CORE_HFAETOMSEXICTOPKPI_H_

#include "PWGHF/Core/HfMlResponse.h"

#include <vector>
namespace o2::analysis
{
template <typename TypeOutputScore = float>
class HfAeToMseXicToPKPi : public HfMlResponse<TypeOutputScore>
{
 public:
  /// Default constructor
  HfAeToMseXicToPKPi() = default;
  /// Default destructor
  virtual ~HfAeToMseXicToPKPi() = default;

  std::vector<float> yScaled, yOutRescaled;

  void setMinMaxScaling(std::vector<float>& yOut, std::vector<float> yIn, std::vector<float> scaleMin, std::vector<float> scaleMax)
  {
    yOut.clear();                             // initial clear to avoid multiple filling if setMinMax o setScaling are called more than once
    for (size_t j = 0; j < yIn.size(); ++j) { // loop for over the features
      // MinMax scaling of the input features
      LOG(debug) << "--------------> MinMax scaling Debug \t" << scaleMin.at(j) << "\t" << scaleMax.at(j);
      yOut.push_back((yIn.at(j) - scaleMin.at(j)) / (scaleMax.at(j) - scaleMin.at(j)));
      LOG(debug) << "Feature = " << j << " ----> input = " << yIn.at(j) << " scaled feature = " << yOut.at(j);
    }
  }
  // ---- External preprocessing scaling
  void setScaling(bool scaleFlag, int scaleType, /* input features of a candidate */ std::vector<float> yIn, std::vector<float> scaleMin, std::vector<float> scaleMax)
  { // it takes the bool flag and scaling parameters configurables in taskXic
    yScaled.clear();
    if (scaleFlag == false) {
      LOG(debug) << "No external preprocessing transformation will be applied";
      yScaled.assign(yIn.begin(), yIn.end());
    } else {
      if (scaleType == 1) {
        LOG(debug) << "MinMax scaling will be applied";
        setMinMaxScaling(yScaled, yIn, scaleMin, scaleMax);
      } // ... with scaleType > 1 we could add other preprocessing trasformations
    }
  }
  std::vector<float> getPreprocessedFeatures()
  {
    for (size_t j = 0; j < yScaled.size(); ++j)
      LOG(debug) << "Global scaled feature = " << yScaled.at(j);
    return yScaled;
  }
  // Reverse preprocessing - output postprocessing
  void unsetMinMaxScaling(std::vector<float>& yOut, std::vector<float> yIn, std::vector<float> scaleMin, std::vector<float> scaleMax)
  {
    yOut.clear();                             // initial clear to avoid multiple filling if setMinMax o setScaling are called more than once
    for (size_t j = 0; j < yIn.size(); ++j ) { // loop for over the features
      // MinMax scaling of the input features
      LOG(debug) << "--------------> MinMax unscaling Debug \t" << scaleMin.at(j) << "\t" << scaleMax.at(j);
      yOut.push_back(yIn.at(j) * (scaleMax.at(j) - scaleMin.at(j)) + scaleMin.at(j));
      LOG(debug) << "Unscaling output = " << j << " ----> input = " << yIn.at(j) << " rescaled output = " << yOut.at(j);
    }
  }

  void unsetScaling(bool scaleFlag, int scaleType, /*AE output*/ std::vector<float> yIn, std::vector<float> scaleMin, std::vector<float> scaleMax)
  { // it takes the bool flag and scaling parameters configurables in taskXic
    yOutRescaled.clear();
    if (scaleFlag == false) {
      LOG(debug) << "No external preprocessing transformation will be applied";
      yOutRescaled.assign(yIn.begin(), yIn.end());
    } else {
      if (scaleType == 1) {
        LOG(debug) << "MinMax unscaling will be applied";
        unsetMinMaxScaling(yOutRescaled, yIn, scaleMin, scaleMax);
      }  //... with scaleType > 1 we could add other preprocessing trasformations
    }
  }

  std::vector<float> getPostprocessedOutput()
  {
    for (size_t j = 0; j < yOutRescaled.size(); ++j)
      LOG(debug)<<"Global rescaled AE output = "<< yOutRescaled.at(j);
    return yOutRescaled;
  }
  //---- MSE function
  float getMse(std::vector<float> yTrue, std::vector<float> yPred)
  {
    LOG(debug) << "Inside getMse sizes " << yTrue.size() << "\t" << yPred.size();
    float mse= 0.0f;
    float sum = 0.0f;
    for (size_t j = 0; j < yTrue.size(); ++j)
      LOG(debug) << "Local Feature = " << j << " ----> input = " << yTrue.at(j) << " scaled feature = " << yPred.at(j);
    std::vector<float> yTrueScaled = getPreprocessedFeatures(); // to make the input features adimensional
    if (yTrue.size() != yPred.size())
    {
      LOG(debug) << "size of input vector =" << yTrue.size();
      LOG(debug) << "size of AE output vector =" << yPred.size();
      LOG(fatal) << "vectors of input and predictions don't have the same size";
    } else { //MSE
      for (size_t j = 0; j < yPred.size(); ++j) {
        sum += std::pow(((yTrueScaled).at(j) - (yPred).at(j)), 2); //AE model gives adimensional predictions by design choice
        LOG(debug) << "getMse Local feature = " << j << " ----> input = " << yTrueScaled.at(j) << " AE prediction = " << yPred.at(j);
      }
      mse = sum/yPred.size(); // MSE of a candidate
      LOG(debug) << "Local mse " << mse;
    }
    return mse;
  }
}; // end of the class

} // namespace o2::analysis
#endif // PWGHF_CORE_HFAETOMSEXICTOPKPI_H_
