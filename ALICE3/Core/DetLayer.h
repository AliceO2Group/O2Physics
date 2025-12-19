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

///
/// \file   DetLayer.h
/// \author David Dobrigkeit Chinellato
/// \since  11/03/2021
/// \brief  Basic struct to hold information regarding a detector layer to be used in fast simulation
///

#ifndef ALICE3_CORE_DETLAYER_H_
#define ALICE3_CORE_DETLAYER_H_

#include <TGraph.h>
#include <TString.h>

#include <string>

namespace o2::fastsim
{

struct DetLayer {
 public:
  // Default constructor
  DetLayer() = default;
  // Parametric constructor
  DetLayer(const TString& name, float r, float z, float x0, float xrho,
           float resRPhi = 0.0f, float resZ = 0.0f, float eff = 0.0f, int type = kLayerInert);
  // Copy constructor
  DetLayer(const DetLayer& other);

  // Setters
  void setName(const TString& name) { mName = name; }
  void setRadius(float r) { mR = r; }
  void setZ(float z) { mZ = z; }
  void setRadiationLength(float x0) { mX0 = x0; }
  void setDensity(float xrho) { mXrho = xrho; }
  void setResolutionRPhi(float resRPhi) { mResRPhi = resRPhi; }
  void setResolutionZ(float resZ) { mResZ = resZ; }
  void setEfficiency(float eff) { mEff = eff; }
  void setType(int type) { mType = type; }

  // Dead areas

  /// @brief Add a dead region in phi for this layer
  /// @param phiStart starting angle in radians of the dead region
  /// @param phiEnd ending angle in radians of the dead region
  void addDeadPhiRegion(float phiStart, float phiEnd);

  /// @brief Set the dead regions in phi for this layer with a TGraph containing all regions. The graph should have y=2 for dead regions and y=0 for alive regions.
  /// @param graph graph of the dead regions. Can be nullptr to clear the dead regions.
  void setDeadPhiRegions(TGraph* graph);

  // Getters
  float getRadius() const { return mR; }
  float getZ() const { return mZ; }
  float getRadiationLength() const { return mX0; }
  float getDensity() const { return mXrho; }
  float getResolutionRPhi() const { return mResRPhi; }
  float getResolutionZ() const { return mResZ; }
  float getEfficiency() const { return mEff; }
  int getType() const { return mType; }
  const TString& getName() const { return mName; }
  const TGraph* getDeadPhiRegions() const { return mDeadPhiRegions; }

  // Check layer mType
  bool isInert() const { return mType == kLayerInert; }
  bool isSilicon() const { return mType == kLayerSilicon; }
  bool isGas() const { return mType == kLayerGas; }
  bool isTOF() const { return mType == kLayerTOF; }
  bool isVertex() const { return mType == kLayerVertex; }
  bool isActive() const { return mType != kLayerInert; } // active layers are not inert

  // Utilities
  std::string toString() const;
  friend std::ostream& operator<<(std::ostream& os, const DetLayer& layer)
  {
    os << layer.toString();
    return os;
  }
  /// @brief Check if a given phi angle is in a dead region
  /// @param phi The phi angle to check
  /// @return True if the phi angle is in a dead region, false otherwise
  bool isInDeadPhiRegion(float phi) const
  {
    if (mDeadPhiRegions == nullptr)
      return false;
    return mDeadPhiRegions->Eval(phi) > 1.f;
  };

  static constexpr int kLayerVertex = -1; // vertex layer type (not used in tracking)
  static constexpr int kLayerInert = 0;   // inert/undefined layer
  static constexpr int kLayerSilicon = 1; // silicon layer
  static constexpr int kLayerGas = 2;     // gas/tpc layer
  static constexpr int kLayerTOF = 3;     // TOF layer type (not used in tracking)

 private:
  // TString for holding name
  TString mName;

  // position variables
  float mR; // radius in centimeters
  float mZ; // mZ dimension in centimeters

  // material variables
  float mX0;   // radiation length
  float mXrho; // density

  // resolution variables for active layers
  float mResRPhi; // RPhi resolution in centimeters
  float mResZ;    // Z resolution in centimeters

  // efficiency
  float mEff; // detection efficiency

  // dead regions in phi (in radians)
  TGraph* mDeadPhiRegions = nullptr;

  // layer type
  int mType; // 0: undefined/inert, 1: silicon, 2: gas/tpc
};

} // namespace o2::fastsim

#endif // ALICE3_CORE_DETLAYER_H_
