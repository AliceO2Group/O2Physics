// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#ifndef PWGLF_UTILS_STRANGENESSBUILDERHELPER_H_
#define PWGLF_UTILS_STRANGENESSBUILDERHELPER_H_

#include <cstdlib>
#include <cmath>
#include <array>
#include "DCAFitter/DCAFitterN.h"
#include "Framework/AnalysisDataModel.h"

namespace o2
{
namespace pwglf
{
//__________________________________________
// V0 information storage
struct v0candidate{
  // indexing
  int collisionId = -1;
  int negativeTrack = -1;
  int positiveTrack = -1;

  // daughter properties
  std::array<float, 3> positiveMomentum = {0.0f, 0.0f, 0.0f};
  std::array<float, 3> negativeMomentum = {0.0f, 0.0f, 0.0f};
  std::array<float, 3> positivePosition = {0.0f, 0.0f, 0.0f};
  std::array<float, 3> negativePosition = {0.0f, 0.0f, 0.0f};
  float positiveTrackX = 0.0f; 
  float negativeTrackX = 0.0f;
  float positiveDCAxy = 0.0f;
  float negativeDCAxy = 0.0f;
  
  // V0 properties
  std::array<float, 3> position = {0.0f, 0.0f, 0.0f};
  float daughterDCA = 1000.0f;
  float pointingAngle = 0.0f;
  float dcaXY = 0.0f;
  
  // calculated masses for convenience
  float massGamma;
  float massK0Short;
  float massLambda;
  float massAntiLambda;

  // stored for decay chains
  std::array<float, 21> covariance;
};

//__________________________________________
// Cascade information storage
struct cascadeCandidate {
  // indexing
  int collisionId = -1;
  int v0Id = -1;
  int negativeTrack = -1;
  int positiveTrack = -1;
  int bachelorTrack = -1;

  // daughter properties
  std::array<float, 3> positiveMomentum = {0.0f, 0.0f, 0.0f};
  std::array<float, 3> negativeMomentum = {0.0f, 0.0f, 0.0f};
  std::array<float, 3> bachelorMomentum = {0.0f, 0.0f, 0.0f};
  std::array<float, 3> positivePosition = {0.0f, 0.0f, 0.0f};
  std::array<float, 3> negativePosition = {0.0f, 0.0f, 0.0f};
  std::array<float, 3> bachelorPosition = {0.0f, 0.0f, 0.0f};
  float positiveDCAxy = 0.0f;
  float negativeDCAxy = 0.0f;
  float bachelorDCAxy = 0.0f;
  float positiveTrackX = 0.0f; 
  float negativeTrackX = 0.0f;
  float bacheloTrackX = 0.0f;

  // cascade properties
  std::array<float, 3> position = {0.0f, 0.0f, 0.0f};
  int charge = -1; // default: []Minus
  float cascadeDaughterDCA = 1000.0f;
  float v0DaughterDCA = 1000.0f;
  
  float pointingAngle = 0.0f;
  float cascadeDCAxy = 0.0f;
  float cascadeDCAz = 0.0f;
  std::array<float, 3> v0Position = {0.0f, 0.0f, 0.0f};
  std::array<float, 3> v0Momentum = {0.0f, 0.0f, 0.0f};
  std::array<float, 3> cascadePosition = {0.0f, 0.0f, 0.0f};
  std::array<float, 3> cascadeMomentum = {0.0f, 0.0f, 0.0f};

  float bachBaryonPointingAngle = 0.0f;
  float bachBaryonDCAxyToPV = 0.0f;
  float massXi = 0.0f; 
  float massOmega = 0.0f; 

  // stored for decay chains
  std::array<float, 21> covariance;
};

//__________________________________________
// builder helper class
class strangenessBuilderHelper
{
  public:
    strangenessBuilderHelper();
    
    template <typename TTrack>
    bool buildV0Candidate(o2::aod::Collision const& collision,
                          TTrack const& positiveTrack, 
                          TTrack const& negativeTrack, 
                          bool useCollinearFit = false);

    o2::base::MatLayerCylSet* lut; // material LUT for DCA fitter
    o2::vertexing::DCAFitterN<2> fitter; // 2-prong o2 dca fitter

    v0candidate v0;  // storage for V0 candidate properties
    cascadeCandidate cascade; // storage for cascade candidate properties

  private: 
    // internal helper to calculate DCAxy of a straight line to a given PV analytically
    float CalculateDCAStraightToPV(float X, float Y, float Z, float Px, float Py, float Pz, float pvX, float pvY, float pvZ);
};

} // namespace pwglf
} // namespace o2

#endif // PWGLF_UTILS_STRANGENESSBUILDERHELPER_H_
