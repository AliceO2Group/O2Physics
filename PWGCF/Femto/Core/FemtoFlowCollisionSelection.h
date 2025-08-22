// Copyright 2019-2025 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file FemtoFlowCollisionSelection.h
/// \brief FemtoFlowCollisionSelection - event selection within the o2femtoflow framework
/// \author Wenya Wu, TU München, wenya.wu@cern.ch
/// \note The femtoflow borrow and copy the framework from femtodream and femtouniverse

#ifndef PWGCF_FEMTOFLOW_CORE_FEMTOFLOWCOLLISIONSELECTION_H_
#define PWGCF_FEMTOFLOW_CORE_FEMTOFLOWCOLLISIONSELECTION_H_

#include <string>
#include "Common/CCDB/TriggerAliases.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/Logger.h"

#include "Common/Core/EventPlaneHelper.h"
#include "Common/DataModel/Qvectors.h"

using namespace o2;
using namespace o2::framework;

namespace o2::analysis::femto_flow
{

/// \class FemtoFlowCollisionSelection
/// \brief Small selection class to check whether a given collision fulfills the specified selections
class FemtoFlowCollisionSelection
{
 public:
  /// Destructor
  virtual ~FemtoFlowCollisionSelection() = default;

  /// Pass the selection criteria to the class
  /// \param zvtxMax Maximal value of the z-vertex
  /// \param checkTrigger Whether or not to check for the trigger alias
  /// \param trig Requested trigger alias
  /// \param checkOffline Whether or not to check for offline selection criteria
  /// \param checkRun3 To check for the Run3 data
  /// \param centmin Minimum value of centrality selection
  /// \param centmax Maximum value of centrality selection
  void setCuts(float zvtxMax, bool checkTrigger, int trig, bool checkOffline, bool checkRun3, float centmin, float centmax)
  // void setCuts(float zvtxMax, bool checkTrigger, int trig, bool checkOffline, bool checkRun3)
  {
    mCutsSet = true;
    mZvtxMax = zvtxMax;
    mCheckTrigger = checkTrigger;
    mTrigger = static_cast<triggerAliases>(trig);
    mCheckOffline = checkOffline;
    mCheckIsRun3 = checkRun3;
    mCentMin = centmin;
    mCentMax = centmax;
  }

  /// Initializes histograms for the task
  /// \param registry Histogram registry to be passed
  void init(HistogramRegistry* registry)
  {
    if (!mCutsSet) {
      LOGF(error, "Event selection not set - quitting!");
    }
    mHistogramRegistry = registry;
    mHistogramRegistry->add("Event/zvtxhist", "; vtx_{z} (cm); Entries", kTH1F, {{300, -12.5, 12.5}});
    mHistogramRegistry->add("Event/MultV0M", "; vMultV0M; Entries", kTH1F, {{16384, 0, 32768}});
    mHistogramRegistry->add("Event/MultT0M", "; vMultT0M; Entries", kTH1F, {{4096, 0, 8192}});
    mHistogramRegistry->add("Event/MultNTracksPV", "; vMultNTracksPV; Entries", kTH1F, {{120, 0, 120}});
    mHistogramRegistry->add("Event/MultNTracklets", "; vMultNTrackslets; Entries", kTH1F, {{300, 0, 300}});
    mHistogramRegistry->add("Event/MultTPC", "; vMultTPC; Entries", kTH1I, {{600, 0, 600}});
    mHistogramRegistry->add("Event/Sphericity", "; Sphericity; Entries", kTH1I, {{200, 0, 3}});
    mHistogramRegistry->add("Event/qnvector", "; Centrality; qn", kTH2F, {{100, 0, 100}, {1000, 0, 1000}});
    mHistogramRegistry->add("Event/SphrVsqn", "; qn; Sphericity", kTH2F, {{1000, 0, 1000}, {200, 0, 3}});
  }

  /// Print some debug information
  void printCuts()
  {
    LOG(info) << "Debug information for FemtoFlowCollisionSelection";
    LOG(info) << "Max. z-vertex: " << mZvtxMax;
    LOG(info) << "Check trigger: " << mCheckTrigger;
    LOG(info) << "Trigger: " << mTrigger;
    LOG(info) << " Check offline: " << mCheckOffline;
    LOG(info) << " Minimum Centrality: " << mCentMin;
    LOG(info) << " Maximum Centrality: " << mCentMax;
  }

  /// Check whether the collisions fulfills the specified selections
  /// \tparam T type of the collision
  /// \param col Collision
  /// \return whether or not the collisions fulfills the specified selections
  template <typename T>
  bool isSelected(T const& col)
  {
    if (std::abs(col.posZ()) > mZvtxMax) {
      return false;
    }
    if ((col.centFT0C() < mCentMin) || (col.centFT0C() > mCentMax)) {
      return false;
    }
    if (mCheckIsRun3) {
      if (mCheckOffline && !col.sel8()) {
        return false;
      }
    } else {
      if (mCheckTrigger && !col.alias_bit(mTrigger)) {
        return false;
      }
      if (mCheckOffline && !col.sel7()) {
        return false;
      }
    }
    return true;
  }

  /// Check whether the collisions fulfills the specified selections for Run3
  /// \tparam T type of the collision
  /// \param col Collision
  /// \return whether or not the collisions fulfills the specified selections
  template <typename T>
  bool isSelectedRun3(T const& col)
  {
    if (std::abs(col.posZ()) > mZvtxMax) {
      return false;
    }
    if (mCheckOffline && !col.sel8()) {
      return false;
    }
    if ((col.centFT0C() < mCentMin) || (col.centFT0C() > mCentMax)) {
      return false;
    }
    return true;
  }

  /// Some basic QA of the event
  /// \tparam T type of the collision
  /// \param col Collision
  template <typename T>
  void fillQA(T const& col)
  {    

    if (mHistogramRegistry) {
      mHistogramRegistry->fill(HIST("Event/zvtxhist"), col.posZ());
      mHistogramRegistry->fill(HIST("Event/MultT0M"), col.multFT0M());
      mHistogramRegistry->fill(HIST("Event/MultNTracksPV"), col.multNTracksPV());
      mHistogramRegistry->fill(HIST("Event/MultNTracklets"), col.multTracklets());
      mHistogramRegistry->fill(HIST("Event/MultTPC"), col.multTPC());
      if (mCheckIsRun3) {
        mHistogramRegistry->fill(HIST("Event/MultV0M"), col.multFV0M());
      } else {
        mHistogramRegistry->fill(HIST("Event/MultV0M"), 0.5 * (col.multFV0M())); // in AliPhysics, the VOM was defined by (V0A + V0C)/2.
      }
    }
  }

  /// Compute the sphericity of an event
  /// \tparam T1 type of the collision
  /// \tparam T2 type of the tracks
  /// \param col Collision
  /// \param tracks All tracks
  /// \return value of the sphericity of the event
  template <typename T1, typename T2>
  float computeSphericity(T1 const& /*col*/, T2 const& tracks)
  {
    double kS00 = 0;
    double kS11 = 0;
    double kS10 = 0;
    double sumPt = 0;
    int partNumber = 0;
    double spher = 0;

    for (const auto& p : tracks) {
      double phi = p.phi();
      double pT = p.pt();
      double px = pT * std::cos(phi);
      double py = pT * std::sin(phi);

      kS00 = kS00 + px * px / pT;
      kS11 = kS11 + py * py / pT;
      kS10 = kS10 + px * py / pT;
      sumPt = sumPt + pT;
      partNumber++;
    }

    if (sumPt != 0) {
      kS00 = kS00 / sumPt;
      kS11 = kS11 / sumPt;
      kS10 = kS10 / sumPt;

      double lambda1 = (kS00 + kS11 + std::sqrt((kS00 + kS11) * (kS00 + kS11) - 4.0 * (kS00 * kS11 - kS10 * kS10))) / 2.0;
      double lambda2 = (kS00 + kS11 - std::sqrt((kS00 + kS11) * (kS00 + kS11) - 4.0 * (kS00 * kS11 - kS10 * kS10))) / 2.0;

      if ((lambda1 + lambda2) != 0 && partNumber > 2) {
        spher = 2 * lambda2 / (lambda1 + lambda2);
      } else {
        spher = 2;
      }
    } else {
      spher = 2;
    }

    if (mHistogramRegistry) {
      mHistogramRegistry->fill(HIST("Event/Sphericity"), spher);
      mSphericity = spher;
    }
    return spher;
  }

  //Qn-vector calculation
  template <typename T>
  float computeqnVec(T const& col)
  {
    double qn = std::sqrt(col.qvecFT0CReVec()[0] * col.qvecFT0CReVec()[0] + col.qvecFT0CImVec()[0] * col.qvecFT0CImVec()[0]) * std::sqrt(col.sumAmplFT0C());
    if (mHistogramRegistry){
      mHistogramRegistry->fill(HIST("Event/qnvector"), col.centFT0C(), qn);  
      mHistogramRegistry->fill(HIST("Event/SphrVsqn"), qn, mSphericity);  
    }
    return qn;
  }

  //Qn-vector calculation
  template <typename T>
  int myqnBin(T const& col, float centBinLength=10.f)
  { 
    int qnBin = -999;
    float qn = computeqnVec(col);
    int mycentBin = (int)(col.centFT0C() / centBinLength);
    if (mycentBin >= (int)(mCentMax / centBinLength)) return qnBin;

    for (int iqn(0); iqn < static_cast<int>(std::size(mqnBinSeparator[mycentBin]))-1; ++iqn){
      if (qn>mqnBinSeparator[mycentBin][iqn] && qn<=mqnBinSeparator[mycentBin][iqn+1]){
        qnBin = iqn;
        break;
      }else continue;
    }

    return qnBin;
  }

 private:
  HistogramRegistry* mHistogramRegistry = nullptr; ///< For QA output
  bool mCutsSet = false;                           ///< Protection against running without cuts
  bool mCheckTrigger = false;                      ///< Check for trigger
  bool mCheckOffline = false;                      ///< Check for offline criteria (might change)
  bool mCheckIsRun3 = false;                       ///< Check if running on Pilot Beam
  triggerAliases mTrigger = kINT7;                 ///< Trigger to check for
  float mZvtxMax = 999.f;                          ///< Maximal deviation from nominal z-vertex (cm)
  float mCentMin = 0.0;                            ///< Minimum centrality value
  float mCentMax = 100.0;                          ///< Maximum centrality value
  float mSphericity = 2.;
  float mqnBinSeparator [7][11] = {{ 0.0,  63.50,  92.50,  116.50,  139.50,  162.50,  185.50,  212.50,  245.50,  292.50,  877.50},
                                   { 0.0,  57.50,  82.50,  102.50,  121.50,  139.50,  158.50,  178.50,  203.50,  238.50,  616.50},
                                   { 0.0,  49.50,  70.50,  86.50,  102.50,  116.50,  131.50,  148.50,  168.50,  195.50,  483.50},
                                   { 0.0,  38.50,  55.50,  69.50,  82.50,  94.50,  106.50,  120.50,  137.50,  160.50,  375.50},
                                   { 0.0,  29.50,  42.50,  53.50,  63.50,  73.50,  83.50,  95.50,  109.50,  128.50,  322.50},
                                   { 0.0,  21.50,  31.50,  39.50,  47.50,  55.50,  63.50,  72.50,  83.50,  99.50,  266.50},
                                   { 0.0,  15.50,  22.50,  28.50,  33.50,  39.50,  45.50,  52.50,  60.50,  72.50,  232.50}
                                  }; ///< qn bin edge from qn vector distributions, for per 10% centrality, 0-70%
};
} // namespace o2::analysis::femto_flow

#endif // PWGCF_FEMTOFLOW_CORE_FEMTOFLOWCOLLISIONSELECTION_H_
