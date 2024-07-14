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

// \brief   Histogram manager for the AC-related analyses.
// \author  Cindy Mordasini (cindy.mordasini@cern.ch)

#ifndef PWGCF_JCORRAN_CORE_FLOWJHISTMANAGER_H_
#define PWGCF_JCORRAN_CORE_FLOWJHISTMANAGER_H_

/* Header files. */
#include <iostream>
#include <vector>
#include <array>
#include <string>
#include <string_view>
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"
#include "TProfile2D.h"

// O2 headers. //
#include "Framework/HistogramRegistry.h"

// O2 Physics headers.

/* Namespaces. */
using namespace o2;
using namespace o2::framework;

// ----------------------------------------------------------------------------
// Histogram manager to fill the general QA common to all flow tasks.
// ----------------------------------------------------------------------------
class FlowJHistManager
{
 public:
  FlowJHistManager() = default;

  // Setters and getters, in the same order as the data members.
  void SetHistRegistryQA(HistogramRegistry* myRegistry)
  {
    mHistRegistryQA = myRegistry;
    LOGF(info, "QA histogram registry successfully set.");
  }
  HistogramRegistry* GetHistRegistryQA() const { return mHistRegistryQA; }

  void SetDebugLog(bool debug)
  {
    mDebugLog = debug;
    LOGF(info, "Debug level: %d", mDebugLog);
  }
  bool GetDebugLog() const { return mDebugLog; }

  void SetObtainNUA(bool nua)
  {
    mObtainNUA = nua;
    LOGF(info, "Obtain 3D Zvtx-eta-phi distribution: %d", mObtainNUA);
  }
  bool GetObtainNUA() const { return mObtainNUA; }

  void SetSaveAllQA(bool saveQA)
  {
    mSaveAllQA = saveQA;
    LOGF(info, "Save the additional QA : %d.", mSaveAllQA);
  }
  bool GetSaveAllQA() const { return mSaveAllQA; }

  void SetSaveQABefore(bool saveQA)
  {
    mSaveQABefore = saveQA;
    LOGF(info, "Save the QA before the selection : %d.", mSaveQABefore);
  }
  bool GetSaveQABefore() const { return mSaveQABefore; }

  void SetUseVariablePtBins(bool myAxis)
  {
    mUseVariablePtBins = myAxis;
    LOGF(info, "Use variable pT binning: %d.", mUseVariablePtBins);
  }
  bool GetUseVariablePtBins() const { return mUseVariablePtBins; }

  /* Methods specific to this class. */
  // The template functions are defined down here to prevent compilation errors.
  void CreateHistQA();
  int GetCentBin(float cValue);

  /// \brief Fill the event QA histograms.
  /// \tparam T Type of collision.
  /// \tparam mode Indicate if Before/ or After/ selection.
  /// \param coll Collision entry of the table.
  /// \param cBin Centrality bin of the collision.
  /// \param cent Centrality percentile of the collision.
  /// \param multi Collision multiplicity at this step.
  template <int mode, typename T>
  void FillEventQA(T const& coll, int cBin, float cent, int multi)
  {
    if (!mHistRegistryQA) {
      LOGF(fatal, "QA histogram registry missing. Quitting...");
      return;
    }

    static constexpr std::string_view subDir[] = {"Before/", "After/"};
    switch (cBin) {
      case 0:
        mHistRegistryQA->fill(HIST(mCentClasses[0]) + HIST(subDir[mode]) + HIST("histCent"), cent);
        mHistRegistryQA->fill(HIST(mCentClasses[0]) + HIST(subDir[mode]) + HIST("histMulti"), multi);
        mHistRegistryQA->fill(HIST(mCentClasses[0]) + HIST(subDir[mode]) + HIST("histZvtx"), coll.posZ());
        break;
      case 1:
        mHistRegistryQA->fill(HIST(mCentClasses[1]) + HIST(subDir[mode]) + HIST("histCent"), cent);
        mHistRegistryQA->fill(HIST(mCentClasses[1]) + HIST(subDir[mode]) + HIST("histMulti"), multi);
        mHistRegistryQA->fill(HIST(mCentClasses[1]) + HIST(subDir[mode]) + HIST("histZvtx"), coll.posZ());
        break;
      case 2:
        mHistRegistryQA->fill(HIST(mCentClasses[2]) + HIST(subDir[mode]) + HIST("histCent"), cent);
        mHistRegistryQA->fill(HIST(mCentClasses[2]) + HIST(subDir[mode]) + HIST("histMulti"), multi);
        mHistRegistryQA->fill(HIST(mCentClasses[2]) + HIST(subDir[mode]) + HIST("histZvtx"), coll.posZ());
        break;
      case 3:
        mHistRegistryQA->fill(HIST(mCentClasses[3]) + HIST(subDir[mode]) + HIST("histCent"), cent);
        mHistRegistryQA->fill(HIST(mCentClasses[3]) + HIST(subDir[mode]) + HIST("histMulti"), multi);
        mHistRegistryQA->fill(HIST(mCentClasses[3]) + HIST(subDir[mode]) + HIST("histZvtx"), coll.posZ());
        break;
      case 4:
        mHistRegistryQA->fill(HIST(mCentClasses[4]) + HIST(subDir[mode]) + HIST("histCent"), cent);
        mHistRegistryQA->fill(HIST(mCentClasses[4]) + HIST(subDir[mode]) + HIST("histMulti"), multi);
        mHistRegistryQA->fill(HIST(mCentClasses[4]) + HIST(subDir[mode]) + HIST("histZvtx"), coll.posZ());
        break;
      case 5:
        mHistRegistryQA->fill(HIST(mCentClasses[5]) + HIST(subDir[mode]) + HIST("histCent"), cent);
        mHistRegistryQA->fill(HIST(mCentClasses[5]) + HIST(subDir[mode]) + HIST("histMulti"), multi);
        mHistRegistryQA->fill(HIST(mCentClasses[5]) + HIST(subDir[mode]) + HIST("histZvtx"), coll.posZ());
        break;
      case 6:
        mHistRegistryQA->fill(HIST(mCentClasses[6]) + HIST(subDir[mode]) + HIST("histCent"), cent);
        mHistRegistryQA->fill(HIST(mCentClasses[6]) + HIST(subDir[mode]) + HIST("histMulti"), multi);
        mHistRegistryQA->fill(HIST(mCentClasses[6]) + HIST(subDir[mode]) + HIST("histZvtx"), coll.posZ());
        break;
      case 7:
        mHistRegistryQA->fill(HIST(mCentClasses[7]) + HIST(subDir[mode]) + HIST("histCent"), cent);
        mHistRegistryQA->fill(HIST(mCentClasses[7]) + HIST(subDir[mode]) + HIST("histMulti"), multi);
        mHistRegistryQA->fill(HIST(mCentClasses[7]) + HIST(subDir[mode]) + HIST("histZvtx"), coll.posZ());
        break;
      case 8:
        mHistRegistryQA->fill(HIST(mCentClasses[8]) + HIST(subDir[mode]) + HIST("histCent"), cent);
        mHistRegistryQA->fill(HIST(mCentClasses[8]) + HIST(subDir[mode]) + HIST("histMulti"), multi);
        mHistRegistryQA->fill(HIST(mCentClasses[8]) + HIST(subDir[mode]) + HIST("histZvtx"), coll.posZ());
        break;
      case 9:
        mHistRegistryQA->fill(HIST(mCentClasses[9]) + HIST(subDir[mode]) + HIST("histCent"), cent);
        mHistRegistryQA->fill(HIST(mCentClasses[9]) + HIST(subDir[mode]) + HIST("histMulti"), multi);
        mHistRegistryQA->fill(HIST(mCentClasses[9]) + HIST(subDir[mode]) + HIST("histZvtx"), coll.posZ());
        break;
    }

    if (mDebugLog) LOGF(info, "The EventQA has been filled.");
  }

  /// \brief Hardcode the cBin for FillThisTrackQA if not constant.
  /// \tparam T Type of track.
  /// \tparam mode Set if we fill Before/ or After/ objects.
  /// \param track Track entry.
  /// \param cBin Centrality bin of the collision.
  /// \param weightNUE Value of the NUE weight to apply to pT.
  /// \param weightNUA Value of the NUA weight to apply to phi.
  template <int mode, typename T>
  void FillTrackQA(T const& track, int cBin,
                   float weightNUE = 1., float weightNUA = 1., float zVtx = 0.)
  {
    if (!mHistRegistryQA) {
      LOGF(fatal, "QA histogram registry missing. Quitting...");
      return;
    }

    switch (cBin) {
      case 0:
        FillThisTrackQA<0, mode>(track, zVtx, weightNUE, weightNUA);
        break;
      case 1:
        FillThisTrackQA<1, mode>(track, zVtx, weightNUE, weightNUA);
        break;
      case 2:
        FillThisTrackQA<2, mode>(track, zVtx, weightNUE, weightNUA);
        break;
      case 3:
        FillThisTrackQA<3, mode>(track, zVtx, weightNUE, weightNUA);
        break;
      case 4:
        FillThisTrackQA<4, mode>(track, zVtx, weightNUE, weightNUA);
        break;
      case 5:
        FillThisTrackQA<5, mode>(track, zVtx, weightNUE, weightNUA);
        break;
      case 6:
        FillThisTrackQA<6, mode>(track, zVtx, weightNUE, weightNUA);
        break;
      case 7:
        FillThisTrackQA<7, mode>(track, zVtx, weightNUE, weightNUA);
        break;
      case 8:
        FillThisTrackQA<8, mode>(track, zVtx, weightNUE, weightNUA);
        break;
      case 9:
        FillThisTrackQA<9, mode>(track, zVtx, weightNUE, weightNUA);
        break;
    }

    if (mDebugLog) {
      LOGF(info, "The TrackQA has been filled for cBin = %d and mode = %d.", cBin, mode);
    }
  }

  /// \brief Fill the track QA histograms in a fixed centrality bin.
  /// \tparam T Type of track.
  /// \tparam cBin Centrality bin of the collision.
  /// \tparam mode Fill the QA before/after the full selection.
  /// \param track Track entry of the table.
  /// \param zVtx Value of the Zvtx of the collision the track belongs to.
  /// \param weightNUE Value of the NUE weight, the inverse is applied to pT.
  /// \param weightNUA Value of the NUA weight, the inverse is applied to phi.
  /// \note This method can be directly used if no switch is previously needed.
  // TODO: Add filling of the weight histograms.
  template <int cBin, int mode, typename T>
  void FillThisTrackQA(T const& track, float zVtx = 0.,
                       float weightNUE = 1., float weightNUA = 1.)
  {
    static constexpr std::string_view subDir[] = {"Before/", "After/"};

    mHistRegistryQA->fill(HIST(mCentClasses[cBin]) + HIST(subDir[mode]) + HIST("histPt"), track.pt());
    mHistRegistryQA->fill(HIST(mCentClasses[cBin]) + HIST(subDir[mode]) + HIST("histEta"), track.eta());
    mHistRegistryQA->fill(HIST(mCentClasses[cBin]) + HIST(subDir[mode]) + HIST("histPhi"), track.phi());
    mHistRegistryQA->fill(HIST(mCentClasses[cBin]) + HIST(subDir[mode]) + HIST("histCharge"), track.sign());

    if (mode == 1) { // 'Weight' distributions are defined only for After/.
      mHistRegistryQA->fill(HIST(mCentClasses[cBin]) + HIST("After/histPtCorrected"),
                            track.pt(), 1. / weightNUE);
      mHistRegistryQA->fill(HIST(mCentClasses[cBin]) + HIST("After/histPhiCorrected"),
                            track.phi(), 1. / weightNUA);

      mHistRegistryQA->fill(HIST(mCentClasses[cBin]) + HIST("After/histNUEWeights"),
                            track.pt(), weightNUE);
      mHistRegistryQA->fill(HIST(mCentClasses[cBin]) + HIST("After/histNUAWeights"),
                            track.phi(), weightNUA);

      // 3D distribution Zvtx-eta-phi.
      if (mObtainNUA) {
        mHistRegistryQA->fill(HIST(mCentClasses[cBin]) + HIST("After/histZvtxEtaPhi"),
                              zVtx, track.eta(), track.phi());
      }
    }

    if (mSaveAllQA) {
      // TPC information.
      mHistRegistryQA->fill(HIST(mCentClasses[cBin]) + HIST(subDir[mode]) + HIST("histTPCNClsFound"),
                            track.tpcNClsFound());
      mHistRegistryQA->fill(HIST(mCentClasses[cBin]) + HIST(subDir[mode]) + HIST("histTPCNClsCrossedRows"),
                            track.tpcNClsCrossedRows());
      mHistRegistryQA->fill(HIST(mCentClasses[cBin]) + HIST(subDir[mode]) + HIST("histTPCCrossedRowsOverFindableCls"),
                            track.tpcCrossedRowsOverFindableCls());
      mHistRegistryQA->fill(HIST(mCentClasses[cBin]) + HIST(subDir[mode]) + HIST("histTPCFoundOverFindableCls"),
                            track.tpcFoundOverFindableCls());
      mHistRegistryQA->fill(HIST(mCentClasses[cBin]) + HIST(subDir[mode]) + HIST("histTPCFractionSharedCls"),
                            track.tpcFractionSharedCls());
      mHistRegistryQA->fill(HIST(mCentClasses[cBin]) + HIST(subDir[mode]) + HIST("histTPCChi2NCl"),
                            track.tpcChi2NCl());

      // ITS information.
      mHistRegistryQA->fill(HIST(mCentClasses[cBin]) + HIST(subDir[mode]) + HIST("histITSNCls"),
                            track.itsNCls());
      mHistRegistryQA->fill(HIST(mCentClasses[cBin]) + HIST(subDir[mode]) + HIST("histITSNClsInnerBarrel"),
                            track.itsNClsInnerBarrel());
      mHistRegistryQA->fill(HIST(mCentClasses[cBin]) + HIST(subDir[mode]) + HIST("histITSChi2NCl"),
                            track.itsChi2NCl());

      // DCA information.
      mHistRegistryQA->fill(HIST(mCentClasses[cBin]) + HIST(subDir[mode]) + HIST("histDCAxy"),
                            track.pt(), track.dcaXY());
      mHistRegistryQA->fill(HIST(mCentClasses[cBin]) + HIST(subDir[mode]) + HIST("histDCAz"),
                            track.dcaZ());
    }

    if (mDebugLog) {
      LOGF(info, "The ThisTrackQA has been filled for cBin = %d and mode = %d.", cBin, mode);
    }
  }

 private:
  HistogramRegistry* mHistRegistryQA = nullptr; ///< For the QA output.

  bool mDebugLog = false;          ///< Enable to print additional log for debug.
  bool mObtainNUA = false;         ///< Enable to get the 3D Zvtx-eta-phi distribution for NUA.
  bool mSaveAllQA = false;         ///< Save the additional QA (true for QA task).
  bool mSaveQABefore = false;      ///< Save the QA output before any selection.
  bool mUseVariablePtBins = false; ///< Enable the use of a variable width pT binning.

  static const int mNcentBins = 10;                    ///< Number of centrality classes.
  static constexpr std::string_view mCentClasses[] = { ///< Centrality classes.
    "Centrality_00-01/", "Centrality_01-02/", "Centrality_02-05/",
    "Centrality_05-10/", "Centrality_10-20/", "Centrality_20-30/",
    "Centrality_30-40/", "Centrality_40-50/", "Centrality_50-60/",
    "Centrality_60-70/"};

  ClassDefNV(FlowJHistManager, 1);
};

#endif // PWGCF_JCORRAN_CORE_FLOWJHISTMANAGER_H_
