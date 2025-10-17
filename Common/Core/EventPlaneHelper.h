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
/// \file   EventPlaneHelper.h
/// \author Cindy Mordasini <cindy.mordasini@cern.ch>
/// \author Anna Önnerstad <anna.onnerstad@cern.ch>
///
/// \brief  Helper methods to calculate information needed for the event plane
///         calculations with FIT.
///

#ifndef COMMON_CORE_EVENTPLANEHELPER_H_
#define COMMON_CORE_EVENTPLANEHELPER_H_

#include <FT0Base/Geometry.h>
#include <FV0Base/Geometry.h>

#include <TComplex.h>
#include <TH2.h>

#include <Rtypes.h>

#include <memory>
#include <vector>

class EventPlaneHelper
{
 public:
  EventPlaneHelper() = default;

  // Setters/getters for the data members.
  void SetOffsetFT0A(double offsetX, double offsetY)
  {
    mOffsetFT0AX = offsetX;
    mOffsetFT0AY = offsetY;
  }
  void SetOffsetFT0C(double offsetX, double offsetY)
  {
    mOffsetFT0CX = offsetX;
    mOffsetFT0CY = offsetY;
  }
  void SetOffsetFV0left(double offsetX, double offsetY)
  {
    mOffsetFV0leftX = offsetX;
    mOffsetFV0leftY = offsetY;
  }
  void SetOffsetFV0right(double offsetX, double offsetY)
  {
    mOffsetFV0rightX = offsetX;
    mOffsetFV0rightY = offsetY;
  }

  // Methods to calculate the azimuthal angles for each part of FIT, given the channel number.
  double GetPhiFT0(int chno, o2::ft0::Geometry ft0geom);
  double GetPhiFV0(int chno, o2::fv0::Geometry* fv0geom);

  // Method to get the Q-vector and sum of amplitudes for any channel in FIT, given
  // the detector and amplitude.
  void SumQvectors(int det, int chno, float ampl, int nmod, TComplex& Qvec, float& sum, o2::ft0::Geometry ft0geom, o2::fv0::Geometry* fv0geom);

  // Method to get the bin corresponding to a centrality percentile, according to the
  // centClasses[] array defined in Tasks/qVectorsQA.cxx.
  // Note: Any change in one task should be reflected in the other.
  int GetCentBin(float cent);

  // Method to apply the vector of correction constant passed as a configurable
  // to the calculated Q-vectors.
  void DoCorrections(float& qx, float& qy, const std::vector<float>& corrections);

  void DoRecenter(float& qx, float& qy, float x0, float y0);
  void DoTwist(float& qx, float& qy, float lp, float lm);
  void DoRescale(float& qx, float& qy, float ap, float am);

  // Method to get the recentering correction on the Qx-Qy distribution.
  void GetCorrRecentering(const std::shared_ptr<TH2> histQ, float& meanX, float& meanY);

  // Method to get the std. deviation on the Qx-Qy distribution.
  void GetCorrWidth(const std::shared_ptr<TH2> histQ, float& stdX, float& stdY);

  // Method to get the twist and rescale correction on the Qx-Qy distribution.
  void GetCorrTwistRecale(const std::shared_ptr<TH2> histQ,
                          float& aPlus, float& aMinus,
                          float& lambdaPlus, float& lambdaMinus);

  // Method to calculate the event plane from the provided (Qx, Qy), for n = 2.
  float GetEventPlane(const float qx, const float qy, int nmode = 2);

  // Method to calculate the resolution R2 for the provided profile.
  float GetResolution(const float RefA, const float RefB, int nmode = 2);

 private:
  double mOffsetFT0AX = 0.;     // X-coordinate of the offset of FT0-A.
  double mOffsetFT0AY = 0.;     // Y-coordinate of the offset of FT0-A.
  double mOffsetFT0CX = 0.;     // X-coordinate of the offset of FT0-C.
  double mOffsetFT0CY = 0.;     // Y-coordinate of the offset of FT0-C.
  double mOffsetFV0leftX = 0.;  // X-coordinate of the offset of FV0-A left.
  double mOffsetFV0leftY = 0.;  // Y-coordinate of the offset of FV0-A left.
  double mOffsetFV0rightX = 0.; // X-coordinate of the offset of FV0-A right.
  double mOffsetFV0rightY = 0.; // Y-coordinate of the offset of FV0-A right.

  ClassDefNV(EventPlaneHelper, 2)
};

#endif // COMMON_CORE_EVENTPLANEHELPER_H_
