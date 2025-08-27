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

#ifndef PWGCF_FEMTO_CORE_FEMTOFLOWCOLLISIONSELECTION_H_
#define PWGCF_FEMTO_CORE_FEMTOFLOWCOLLISIONSELECTION_H_

#include "Common/CCDB/TriggerAliases.h"
#include "Common/Core/EventPlaneHelper.h"
#include "Common/DataModel/Qvectors.h"

#include "Framework/HistogramRegistry.h"
#include "Framework/Logger.h"

#include <string>

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
  void init(o2::framework::HistogramRegistry* registry)
  {
    using namespace o2::framework;

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
    using namespace o2::framework;

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

    const double zeroVal = 0;
    const double countPartNumLimit = 2;
    const double defualtSphr = 2;

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

      if ((lambda1 + lambda2) != zeroVal && partNumber > countPartNumLimit) {
        spher = 2 * lambda2 / (lambda1 + lambda2);
      } else {
        spher = defualtSphr;
      }
    } else {
      spher = defualtSphr;
    }

    if (mHistogramRegistry) {
      mHistogramRegistry->fill(HIST("Event/Sphericity"), spher);
      mSphericity = spher;
    }
    return spher;
  }

  // Qn-vector calculation
  template <typename T>
  float computeqnVec(T const& col)
  {
    double qn = std::sqrt(col.qvecFT0CReVec()[0] * col.qvecFT0CReVec()[0] + col.qvecFT0CImVec()[0] * col.qvecFT0CImVec()[0]) * std::sqrt(col.sumAmplFT0C());
    if (mHistogramRegistry) {
      mHistogramRegistry->fill(HIST("Event/qnvector"), col.centFT0C(), qn);
      mHistogramRegistry->fill(HIST("Event/SphrVsqn"), qn, mSphericity);
    }
    return qn;
  }

  // Qn-vector calculation
  template <typename T>
  int myqnBin(T const& col, float centBinLength = 1.f)
  {
    int qnBin = -999;
    float qn = computeqnVec(col);
    int mycentBin = static_cast<int>(col.centFT0C() / centBinLength);
    if (mycentBin >= static_cast<int>(mCentMax / centBinLength))
      return qnBin;

    for (int iqn(0); iqn < static_cast<int>(std::size(mqnBinSeparator[mycentBin])) - 1; ++iqn) {
      if (qn > mqnBinSeparator[mycentBin][iqn] && qn <= mqnBinSeparator[mycentBin][iqn + 1]) {
        qnBin = iqn;
        break;
      } else {
        continue;
      }
    }

    return qnBin;
  }

 private:
  o2::framework::HistogramRegistry* mHistogramRegistry = nullptr; ///< For QA output
  bool mCutsSet = false;                                          ///< Protection against running without cuts
  bool mCheckTrigger = false;                                     ///< Check for trigger
  bool mCheckOffline = false;                                     ///< Check for offline criteria (might change)
  bool mCheckIsRun3 = false;                                      ///< Check if running on Pilot Beam
  triggerAliases mTrigger = kINT7;                                ///< Trigger to check for
  float mZvtxMax = 999.f;                                         ///< Maximal deviation from nominal z-vertex (cm)
  float mCentMin = 0.0;                                           ///< Minimum centrality value
  float mCentMax = 98.0;                                          ///< Maximum centrality value
  float mSphericity = 2.;
  float mqnBinSeparator[98][11] = {
    {0.50, 68.50, 100.50, 126.50, 151.50, 176.50, 203.50, 232.50, 269.50, 322.50, 833.50}, // cent 0 to 1
    {0.50, 66.50, 97.50, 122.50, 147.50, 171.50, 197.50, 226.50, 261.50, 313.50, 821.50},  // cent 1 to 2
    {0.50, 65.50, 95.50, 120.50, 144.50, 168.50, 193.50, 221.50, 256.50, 307.50, 876.50},  // cent 2 to 3
    {0.50, 64.50, 93.50, 118.50, 141.50, 165.50, 190.50, 217.50, 251.50, 300.50, 836.50},  // cent 3 to 4
    {0.50, 63.50, 92.50, 116.50, 139.50, 162.50, 186.50, 214.50, 247.50, 294.50, 732.50},  // cent 4 to 5
    {0.50, 62.50, 91.50, 115.50, 137.50, 160.50, 183.50, 210.50, 242.50, 288.50, 692.50},  // cent 5 to 6
    {0.50, 62.50, 90.50, 114.50, 136.50, 158.50, 181.50, 207.50, 238.50, 283.50, 696.50},  // cent 6 to 7
    {0.50, 61.50, 89.50, 113.50, 134.50, 156.50, 178.50, 203.50, 233.50, 277.50, 681.50},  // cent 7 to 8
    {0.50, 61.50, 88.50, 111.50, 133.50, 154.50, 176.50, 200.50, 229.50, 271.50, 648.50},  // cent 8 to 9
    {0.50, 61.50, 88.50, 110.50, 131.50, 152.50, 173.50, 197.50, 225.50, 265.50, 616.50},  // cent 9 to 10
    {0.50, 60.50, 87.50, 109.50, 129.50, 149.50, 170.50, 193.50, 221.50, 260.50, 615.50},  // cent 10 to 11
    {0.50, 59.50, 86.50, 108.50, 128.50, 147.50, 168.50, 190.50, 217.50, 254.50, 586.50},  // cent 11 to 12
    {0.50, 59.50, 85.50, 106.50, 126.50, 145.50, 165.50, 187.50, 213.50, 249.50, 583.50},  // cent 12 to 13
    {0.50, 58.50, 84.50, 105.50, 124.50, 143.50, 162.50, 183.50, 209.50, 244.50, 542.50},  // cent 13 to 14
    {0.50, 58.50, 83.50, 104.50, 122.50, 141.50, 160.50, 180.50, 205.50, 239.50, 544.50},  // cent 14 to 15
    {0.50, 57.50, 82.50, 102.50, 120.50, 138.50, 157.50, 177.50, 201.50, 234.50, 534.50},  // cent 15 to 16
    {0.50, 56.50, 81.50, 101.50, 119.50, 136.50, 154.50, 174.50, 197.50, 230.50, 530.50},  // cent 16 to 17
    {0.50, 56.50, 80.50, 99.50, 117.50, 134.50, 151.50, 170.50, 193.50, 225.50, 508.50},   // cent 17 to 18
    {0.50, 55.50, 78.50, 98.50, 115.50, 132.50, 149.50, 167.50, 190.50, 221.50, 478.50},   // cent 18 to 19
    {0.50, 54.50, 77.50, 96.50, 113.50, 129.50, 146.50, 164.50, 186.50, 216.50, 508.50},   // cent 19 to 20
    {0.50, 53.50, 76.50, 94.50, 111.50, 127.50, 143.50, 161.50, 183.50, 212.50, 482.50},   // cent 20 to 21
    {0.50, 52.50, 75.50, 93.50, 109.50, 125.50, 141.50, 158.50, 179.50, 208.50, 468.50},   // cent 21 to 22
    {0.50, 51.50, 73.50, 91.50, 107.50, 122.50, 138.50, 155.50, 176.50, 204.50, 436.50},   // cent 22 to 23
    {0.50, 50.50, 72.50, 89.50, 105.50, 120.50, 136.50, 152.50, 172.50, 200.50, 440.50},   // cent 23 to 24
    {0.50, 49.50, 71.50, 88.50, 103.50, 118.50, 133.50, 149.50, 169.50, 196.50, 441.50},   // cent 24 to 25
    {0.50, 48.50, 69.50, 86.50, 101.50, 115.50, 130.50, 146.50, 166.50, 193.50, 412.50},   // cent 25 to 26
    {0.50, 47.50, 68.50, 84.50, 99.50, 113.50, 128.50, 144.50, 162.50, 189.50, 410.50},    // cent 26 to 27
    {0.50, 46.50, 66.50, 82.50, 97.50, 111.50, 125.50, 141.50, 159.50, 185.50, 409.50},    // cent 27 to 28
    {0.50, 46.50, 65.50, 81.50, 95.50, 109.50, 123.50, 138.50, 156.50, 182.50, 405.50},    // cent 28 to 29
    {0.50, 44.50, 63.50, 79.50, 93.50, 106.50, 120.50, 135.50, 153.50, 178.50, 392.50},    // cent 29 to 30
    {0.50, 43.50, 62.50, 77.50, 91.50, 104.50, 118.50, 132.50, 150.50, 174.50, 367.50},    // cent 30 to 31
    {0.50, 42.50, 61.50, 75.50, 89.50, 102.50, 115.50, 130.50, 147.50, 171.50, 368.50},    // cent 31 to 32
    {0.50, 41.50, 59.50, 74.50, 87.50, 100.50, 113.50, 127.50, 144.50, 168.50, 371.50},    // cent 32 to 33
    {0.50, 40.50, 58.50, 72.50, 85.50, 97.50, 110.50, 124.50, 141.50, 164.50, 369.50},     // cent 33 to 34
    {0.50, 39.50, 56.50, 70.50, 83.50, 95.50, 108.50, 121.50, 138.50, 161.50, 363.50},     // cent 34 to 35
    {0.50, 38.50, 55.50, 69.50, 81.50, 93.50, 105.50, 119.50, 135.50, 158.50, 353.50},     // cent 35 to 36
    {0.50, 37.50, 54.50, 67.50, 79.50, 91.50, 103.50, 116.50, 132.50, 154.50, 333.50},     // cent 36 to 37
    {0.50, 36.50, 52.50, 65.50, 77.50, 89.50, 101.50, 114.50, 129.50, 151.50, 333.50},     // cent 37 to 38
    {0.50, 35.50, 51.50, 64.50, 75.50, 87.50, 98.50, 111.50, 126.50, 148.50, 330.50},      // cent 38 to 39
    {0.50, 34.50, 49.50, 62.50, 73.50, 85.50, 96.50, 109.50, 124.50, 145.50, 322.50},      // cent 39 to 40
    {0.50, 33.50, 48.50, 60.50, 71.50, 82.50, 94.50, 106.50, 121.50, 142.50, 319.50},      // cent 40 to 41
    {0.50, 32.50, 47.50, 59.50, 70.50, 80.50, 91.50, 104.50, 118.50, 138.50, 321.50},      // cent 41 to 42
    {0.50, 31.50, 46.50, 57.50, 68.50, 78.50, 89.50, 101.50, 115.50, 135.50, 309.50},      // cent 42 to 43
    {0.50, 31.50, 44.50, 56.50, 66.50, 76.50, 87.50, 99.50, 113.50, 132.50, 311.50},       // cent 43 to 44
    {0.50, 30.50, 43.50, 54.50, 64.50, 74.50, 85.50, 96.50, 110.50, 129.50, 293.50},       // cent 44 to 45
    {0.50, 29.50, 42.50, 53.50, 63.50, 72.50, 82.50, 94.50, 107.50, 126.50, 307.50},       // cent 45 to 46
    {0.50, 28.50, 41.50, 51.50, 61.50, 70.50, 80.50, 91.50, 104.50, 123.50, 290.50},       // cent 46 to 47
    {0.50, 27.50, 39.50, 50.50, 59.50, 68.50, 78.50, 89.50, 102.50, 120.50, 277.50},       // cent 47 to 48
    {0.50, 26.50, 38.50, 48.50, 57.50, 67.50, 76.50, 87.50, 99.50, 117.50, 285.50},        // cent 48 to 49
    {0.50, 25.50, 37.50, 47.50, 56.50, 65.50, 74.50, 84.50, 97.50, 114.50, 264.50},        // cent 49 to 50
    {0.50, 25.50, 36.50, 45.50, 54.50, 63.50, 72.50, 82.50, 94.50, 111.50, 265.50},        // cent 50 to 51
    {0.50, 24.50, 35.50, 44.50, 53.50, 61.50, 70.50, 80.50, 91.50, 108.50, 254.50},        // cent 51 to 52
    {0.50, 23.50, 34.50, 43.50, 51.50, 59.50, 68.50, 77.50, 89.50, 105.50, 256.50},        // cent 52 to 53
    {0.50, 22.50, 33.50, 41.50, 50.50, 58.50, 66.50, 75.50, 86.50, 103.50, 235.50},        // cent 53 to 54
    {0.50, 22.50, 32.50, 40.50, 48.50, 56.50, 64.50, 73.50, 84.50, 100.50, 245.50},        // cent 54 to 55
    {0.50, 21.50, 31.50, 39.50, 47.50, 54.50, 62.50, 71.50, 82.50, 97.50, 239.50},         // cent 55 to 56
    {0.50, 20.50, 30.50, 38.50, 45.50, 52.50, 60.50, 69.50, 79.50, 94.50, 227.50},         // cent 56 to 57
    {0.50, 20.50, 29.50, 36.50, 44.50, 51.50, 58.50, 67.50, 77.50, 91.50, 229.50},         // cent 57 to 58
    {0.50, 19.50, 28.50, 35.50, 42.50, 49.50, 56.50, 65.50, 74.50, 89.50, 227.50},         // cent 58 to 59
    {0.50, 18.50, 27.50, 34.50, 41.50, 48.50, 55.50, 62.50, 72.50, 86.50, 216.50},         // cent 59 to 60
    {0.50, 18.50, 26.50, 33.50, 39.50, 46.50, 53.50, 60.50, 70.50, 83.50, 231.50},         // cent 60 to 61
    {0.50, 17.50, 25.50, 32.50, 38.50, 44.50, 51.50, 58.50, 67.50, 80.50, 194.50},         // cent 61 to 62
    {0.50, 17.50, 24.50, 31.50, 37.50, 43.50, 49.50, 57.50, 65.50, 78.50, 190.50},         // cent 62 to 63
    {0.50, 16.50, 23.50, 30.50, 36.50, 41.50, 48.50, 55.50, 63.50, 75.50, 200.50},         // cent 63 to 64
    {0.50, 15.50, 23.50, 29.50, 34.50, 40.50, 46.50, 53.50, 61.50, 73.50, 183.50},         // cent 64 to 65
    {0.50, 15.50, 22.50, 28.50, 33.50, 39.50, 44.50, 51.50, 59.50, 70.50, 187.50},         // cent 65 to 66
    {0.50, 14.50, 21.50, 27.50, 32.50, 37.50, 43.50, 49.50, 57.50, 68.50, 199.50},         // cent 66 to 67
    {0.50, 14.50, 20.50, 26.50, 31.50, 36.50, 41.50, 47.50, 55.50, 65.50, 171.50},         // cent 67 to 68
    {0.50, 13.50, 19.50, 25.50, 30.50, 34.50, 40.50, 45.50, 53.50, 63.50, 157.50},         // cent 68 to 69
    {0.50, 13.50, 19.50, 24.50, 28.50, 33.50, 38.50, 44.50, 51.50, 60.50, 156.50},         // cent 69 to 70
    {0.50, 12.50, 18.50, 23.50, 27.50, 32.50, 37.50, 42.50, 49.50, 58.50, 157.50},         // cent 70 to 71
    {0.50, 12.50, 17.50, 22.50, 26.50, 31.50, 35.50, 40.50, 47.50, 56.50, 148.50},         // cent 71 to 72
    {0.50, 11.50, 16.50, 21.50, 25.50, 29.50, 34.50, 39.50, 45.50, 54.50, 218.50},         // cent 72 to 73
    {0.50, 11.50, 16.50, 20.50, 24.50, 28.50, 32.50, 37.50, 43.50, 52.50, 201.50},         // cent 73 to 74
    {0.50, 10.50, 15.50, 19.50, 23.50, 27.50, 31.50, 36.50, 41.50, 49.50, 185.50},         // cent 74 to 75
    {0.50, 10.50, 14.50, 18.50, 22.50, 26.50, 30.50, 34.50, 40.50, 47.50, 169.50},         // cent 75 to 76
    {0.50, 9.50, 14.50, 18.50, 21.50, 25.50, 29.50, 33.50, 38.50, 45.50, 156.50},          // cent 76 to 77
    {0.50, 9.50, 13.50, 17.50, 20.50, 24.50, 27.50, 31.50, 36.50, 43.50, 150.50},          // cent 77 to 78
    {0.50, 9.50, 13.50, 16.50, 19.50, 23.50, 26.50, 30.50, 35.50, 42.50, 138.50},          // cent 78 to 79
    {0.50, 8.50, 12.50, 15.50, 19.50, 22.50, 25.50, 29.50, 33.50, 40.50, 127.50},          // cent 79 to 80
    {0.50, 8.50, 12.50, 15.50, 18.50, 21.50, 24.50, 27.50, 32.50, 38.50, 119.50},          // cent 80 to 81
    {0.50, 7.50, 11.50, 14.50, 17.50, 20.50, 23.50, 26.50, 30.50, 36.50, 105.50},          // cent 81 to 82
    {0.50, 7.50, 10.50, 13.50, 16.50, 19.50, 22.50, 25.50, 29.50, 35.50, 107.50},          // cent 82 to 83
    {0.50, 7.50, 10.50, 13.50, 15.50, 18.50, 21.50, 24.50, 28.50, 33.50, 95.50},           // cent 83 to 84
    {0.50, 6.50, 10.50, 12.50, 15.50, 17.50, 20.50, 23.50, 26.50, 31.50, 94.50},           // cent 84 to 85
    {0.50, 6.50, 9.50, 12.50, 14.50, 16.50, 19.50, 22.50, 25.50, 30.50, 89.50},            // cent 85 to 86
    {0.50, 6.50, 9.50, 11.50, 13.50, 15.50, 18.50, 20.50, 24.50, 28.50, 79.50},            // cent 86 to 87
    {0.50, 5.50, 8.50, 10.50, 13.50, 15.50, 17.50, 19.50, 22.50, 27.50, 82.50},            // cent 87 to 88
    {0.50, 5.50, 8.50, 10.50, 12.50, 14.50, 16.50, 18.50, 21.50, 25.50, 91.50},            // cent 88 to 89
    {0.50, 5.50, 7.50, 9.50, 11.50, 13.50, 15.50, 17.50, 20.50, 24.50, 74.50},             // cent 89 to 90
    {0.50, 5.50, 7.50, 9.50, 11.50, 12.50, 14.50, 16.50, 19.50, 23.50, 72.50},             // cent 90 to 91
    {0.50, 4.50, 6.50, 8.50, 10.50, 12.50, 14.50, 16.50, 18.50, 21.50, 63.50},             // cent 91 to 92
    {0.50, 4.50, 6.50, 8.50, 9.50, 11.50, 13.50, 15.50, 17.50, 20.50, 70.50},              // cent 92 to 93
    {0.50, 4.50, 6.50, 7.50, 9.50, 10.50, 12.50, 14.50, 16.50, 19.50, 56.50},              // cent 93 to 94
    {0.50, 4.50, 5.50, 7.50, 8.50, 10.50, 11.50, 13.50, 15.50, 18.50, 62.50},              // cent 94 to 95
    {0.50, 3.50, 5.50, 6.50, 8.50, 9.50, 10.50, 12.50, 14.50, 16.50, 55.50},               // cent 95 to 96
    {0.50, 3.50, 5.50, 6.50, 7.50, 8.50, 10.50, 11.50, 13.50, 15.50, 54.50},               // cent 96 to 97
    {0.50, 3.50, 4.50, 5.50, 6.50, 7.50, 9.50, 10.50, 11.50, 13.50, 44.50}                 // cent 97 to 98
  }; ///< qn bin edge from qn vector distributions, for per 1% centrality, 0-98%
};
} // namespace o2::analysis::femto_flow

#endif // PWGCF_FEMTO_CORE_FEMTOFLOWCOLLISIONSELECTION_H_
