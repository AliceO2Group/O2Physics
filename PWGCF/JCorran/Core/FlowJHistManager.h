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

// ----------------------------------------------------------------------------
// Histogram manager to fill the general QA common to all flow tasks.
// ----------------------------------------------------------------------------
class FlowJHistManager
{
 public:
  FlowJHistManager() = default;

  // Setters and getters, in the same order as the data members.
  void setHistRegistryQA(o2::framework::HistogramRegistry* myRegistry)
  {
    mHistRegistryQA = myRegistry;
    LOGF(info, "QA histogram registry successfully set.");
  }
  o2::framework::HistogramRegistry* getHistRegistryQA() const { return mHistRegistryQA; }

  void setDebugLog(bool debug)
  {
    mDebugLog = debug;
    LOGF(info, "Debug level: %d", mDebugLog);
  }
  bool getDebugLog() const { return mDebugLog; }

  void setObtainNUA(bool nua)
  {
    mObtainNUA = nua;
    LOGF(info, "Obtain 3D Zvtx-eta-phi distribution: %d", mObtainNUA);
  }
  bool getObtainNUA() const { return mObtainNUA; }

  void setSaveAllQA(bool saveQA)
  {
    mSaveAllQA = saveQA;
    LOGF(info, "Save the additional QA : %d.", mSaveAllQA);
  }
  bool getSaveAllQA() const { return mSaveAllQA; }

  void setSaveQABefore(bool saveQA)
  {
    mSaveQABefore = saveQA;
    LOGF(info, "Save the QA before the selection : %d.", mSaveQABefore);
  }
  bool getSaveQABefore() const { return mSaveQABefore; }

  void setUseVariablePtBins(bool myAxis)
  {
    mUseVariablePtBins = myAxis;
    LOGF(info, "Use variable pT binning: %d.", mUseVariablePtBins);
  }
  bool getUseVariablePtBins() const { return mUseVariablePtBins; }

  /* Methods specific to this class. */
  // The template functions are defined down here to prevent compilation errors.
  void createHistQA();
  int getCentBin(float cValue);

  /// \brief Fill the event QA histograms.
  /// \tparam T Type of collision.
  /// \tparam mode Indicate if Before/ or After/ selection.
  /// \param coll Collision entry of the table.
  /// \param cBin Centrality bin of the collision.
  /// \param cent Centrality percentile of the collision.
  /// \param multi Collision multiplicity at this step.
  template <int mode, typename T>
  void fillEventQA(T const& coll, int cBin, float cent, int multi)
  {
    if (!mHistRegistryQA) {
      LOGF(fatal, "QA histogram registry missing. Quitting...");
      return;
    }

    static constexpr std::string_view SubDir[] = {"Before/", "After/"};
    switch (cBin) {
      case 0:
        mHistRegistryQA->fill(HIST(MCentClasses[0]) + HIST(SubDir[mode]) + HIST("histCent"), cent);
        mHistRegistryQA->fill(HIST(MCentClasses[0]) + HIST(SubDir[mode]) + HIST("histMulti"), multi);
        mHistRegistryQA->fill(HIST(MCentClasses[0]) + HIST(SubDir[mode]) + HIST("histZvtx"), coll.posZ());
        break;
      case 1:
        mHistRegistryQA->fill(HIST(MCentClasses[1]) + HIST(SubDir[mode]) + HIST("histCent"), cent);
        mHistRegistryQA->fill(HIST(MCentClasses[1]) + HIST(SubDir[mode]) + HIST("histMulti"), multi);
        mHistRegistryQA->fill(HIST(MCentClasses[1]) + HIST(SubDir[mode]) + HIST("histZvtx"), coll.posZ());
        break;
      case 2:
        mHistRegistryQA->fill(HIST(MCentClasses[2]) + HIST(SubDir[mode]) + HIST("histCent"), cent);
        mHistRegistryQA->fill(HIST(MCentClasses[2]) + HIST(SubDir[mode]) + HIST("histMulti"), multi);
        mHistRegistryQA->fill(HIST(MCentClasses[2]) + HIST(SubDir[mode]) + HIST("histZvtx"), coll.posZ());
        break;
      case 3:
        mHistRegistryQA->fill(HIST(MCentClasses[3]) + HIST(SubDir[mode]) + HIST("histCent"), cent);
        mHistRegistryQA->fill(HIST(MCentClasses[3]) + HIST(SubDir[mode]) + HIST("histMulti"), multi);
        mHistRegistryQA->fill(HIST(MCentClasses[3]) + HIST(SubDir[mode]) + HIST("histZvtx"), coll.posZ());
        break;
      case 4:
        mHistRegistryQA->fill(HIST(MCentClasses[4]) + HIST(SubDir[mode]) + HIST("histCent"), cent);
        mHistRegistryQA->fill(HIST(MCentClasses[4]) + HIST(SubDir[mode]) + HIST("histMulti"), multi);
        mHistRegistryQA->fill(HIST(MCentClasses[4]) + HIST(SubDir[mode]) + HIST("histZvtx"), coll.posZ());
        break;
      case 5:
        mHistRegistryQA->fill(HIST(MCentClasses[5]) + HIST(SubDir[mode]) + HIST("histCent"), cent);
        mHistRegistryQA->fill(HIST(MCentClasses[5]) + HIST(SubDir[mode]) + HIST("histMulti"), multi);
        mHistRegistryQA->fill(HIST(MCentClasses[5]) + HIST(SubDir[mode]) + HIST("histZvtx"), coll.posZ());
        break;
      case 6:
        mHistRegistryQA->fill(HIST(MCentClasses[6]) + HIST(SubDir[mode]) + HIST("histCent"), cent);
        mHistRegistryQA->fill(HIST(MCentClasses[6]) + HIST(SubDir[mode]) + HIST("histMulti"), multi);
        mHistRegistryQA->fill(HIST(MCentClasses[6]) + HIST(SubDir[mode]) + HIST("histZvtx"), coll.posZ());
        break;
      case 7:
        mHistRegistryQA->fill(HIST(MCentClasses[7]) + HIST(SubDir[mode]) + HIST("histCent"), cent);
        mHistRegistryQA->fill(HIST(MCentClasses[7]) + HIST(SubDir[mode]) + HIST("histMulti"), multi);
        mHistRegistryQA->fill(HIST(MCentClasses[7]) + HIST(SubDir[mode]) + HIST("histZvtx"), coll.posZ());
        break;
      case 8:
        mHistRegistryQA->fill(HIST(MCentClasses[8]) + HIST(SubDir[mode]) + HIST("histCent"), cent);
        mHistRegistryQA->fill(HIST(MCentClasses[8]) + HIST(SubDir[mode]) + HIST("histMulti"), multi);
        mHistRegistryQA->fill(HIST(MCentClasses[8]) + HIST(SubDir[mode]) + HIST("histZvtx"), coll.posZ());
        break;
      case 9:
        mHistRegistryQA->fill(HIST(MCentClasses[9]) + HIST(SubDir[mode]) + HIST("histCent"), cent);
        mHistRegistryQA->fill(HIST(MCentClasses[9]) + HIST(SubDir[mode]) + HIST("histMulti"), multi);
        mHistRegistryQA->fill(HIST(MCentClasses[9]) + HIST(SubDir[mode]) + HIST("histZvtx"), coll.posZ());
        break;
    }

    if (mDebugLog)
      LOGF(info, "The EventQA has been filled.");
  }

  /// \brief Hardcode the cBin for fillThisTrackQA if not constant.
  /// \tparam T Type of track.
  /// \tparam mode Set if we fill Before/ or After/ objects.
  /// \param track Track entry.
  /// \param cBin Centrality bin of the collision.
  /// \param weightNUE Value of the NUE weight to apply to pT.
  /// \param weightNUA Value of the NUA weight to apply to phi.
  template <int mode, typename T>
  void fillTrackQA(T const& track, int cBin,
                   float weightNUE = 1., float weightNUA = 1., float zVtx = 0.)
  {
    if (!mHistRegistryQA) {
      LOGF(fatal, "QA histogram registry missing. Quitting...");
      return;
    }

    switch (cBin) {
      case 0:
        fillThisTrackQA<0, mode>(track, zVtx, weightNUE, weightNUA);
        break;
      case 1:
        fillThisTrackQA<1, mode>(track, zVtx, weightNUE, weightNUA);
        break;
      case 2:
        fillThisTrackQA<2, mode>(track, zVtx, weightNUE, weightNUA);
        break;
      case 3:
        fillThisTrackQA<3, mode>(track, zVtx, weightNUE, weightNUA);
        break;
      case 4:
        fillThisTrackQA<4, mode>(track, zVtx, weightNUE, weightNUA);
        break;
      case 5:
        fillThisTrackQA<5, mode>(track, zVtx, weightNUE, weightNUA);
        break;
      case 6:
        fillThisTrackQA<6, mode>(track, zVtx, weightNUE, weightNUA);
        break;
      case 7:
        fillThisTrackQA<7, mode>(track, zVtx, weightNUE, weightNUA);
        break;
      case 8:
        fillThisTrackQA<8, mode>(track, zVtx, weightNUE, weightNUA);
        break;
      case 9:
        fillThisTrackQA<9, mode>(track, zVtx, weightNUE, weightNUA);
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
  void fillThisTrackQA(T const& track, float zVtx = 0.,
                       float weightNUE = 1., float weightNUA = 1.)
  {
    static constexpr std::string_view SubDir[] = {"Before/", "After/"};

    mHistRegistryQA->fill(HIST(MCentClasses[cBin]) + HIST(SubDir[mode]) + HIST("histPt"), track.pt());
    mHistRegistryQA->fill(HIST(MCentClasses[cBin]) + HIST(SubDir[mode]) + HIST("histEta"), track.eta());
    mHistRegistryQA->fill(HIST(MCentClasses[cBin]) + HIST(SubDir[mode]) + HIST("histPhi"), track.phi());
    mHistRegistryQA->fill(HIST(MCentClasses[cBin]) + HIST(SubDir[mode]) + HIST("histCharge"), track.sign());

    if (mode == 1) { // 'Weight' distributions are defined only for After/.
      mHistRegistryQA->fill(HIST(MCentClasses[cBin]) + HIST("After/histPtCorrected"),
                            track.pt(), 1. / weightNUE);
      mHistRegistryQA->fill(HIST(MCentClasses[cBin]) + HIST("After/histPhiCorrected"),
                            track.phi(), 1. / weightNUA);

      mHistRegistryQA->fill(HIST(MCentClasses[cBin]) + HIST("After/histNUEWeights"),
                            track.pt(), weightNUE);
      mHistRegistryQA->fill(HIST(MCentClasses[cBin]) + HIST("After/histNUAWeights"),
                            track.phi(), weightNUA);

      // 3D distribution Zvtx-eta-phi.
      if (mObtainNUA) {
        mHistRegistryQA->fill(HIST(MCentClasses[cBin]) + HIST("After/histZvtxEtaPhi"),
                              zVtx, track.eta(), track.phi());
      }
    }

    if (mSaveAllQA) {
      // TPC information.
      /*
      mHistRegistryQA->fill(HIST(MCentClasses[cBin]) + HIST(SubDir[mode]) + HIST("histTPCNClsFound"),
                            track.tpcNClsFound());
      mHistRegistryQA->fill(HIST(MCentClasses[cBin]) + HIST(SubDir[mode]) + HIST("histTPCNClsCrossedRows"),
                            track.tpcNClsCrossedRows());
      mHistRegistryQA->fill(HIST(MCentClasses[cBin]) + HIST(SubDir[mode]) + HIST("histTPCCrossedRowsOverFindableCls"),
                            track.tpcCrossedRowsOverFindableCls());
      mHistRegistryQA->fill(HIST(MCentClasses[cBin]) + HIST(SubDir[mode]) + HIST("histTPCFoundOverFindableCls"),
                            track.tpcFoundOverFindableCls());
      mHistRegistryQA->fill(HIST(MCentClasses[cBin]) + HIST(SubDir[mode]) + HIST("histTPCFractionSharedCls"),
                            track.tpcFractionSharedCls());
      mHistRegistryQA->fill(HIST(MCentClasses[cBin]) + HIST(SubDir[mode]) + HIST("histTPCChi2NCl"),
                            track.tpcChi2NCl());

      // ITS information.
      mHistRegistryQA->fill(HIST(MCentClasses[cBin]) + HIST(SubDir[mode]) + HIST("histITSNCls"),
                            track.itsNCls());
      mHistRegistryQA->fill(HIST(MCentClasses[cBin]) + HIST(SubDir[mode]) + HIST("histITSNClsInnerBarrel"),
                            track.itsNClsInnerBarrel());
      mHistRegistryQA->fill(HIST(MCentClasses[cBin]) + HIST(SubDir[mode]) + HIST("histITSChi2NCl"),
                            track.itsChi2NCl());

      // DCA information.
      mHistRegistryQA->fill(HIST(MCentClasses[cBin]) + HIST(SubDir[mode]) + HIST("histDCAxy"),
                            track.pt(), track.dcaXY());
      mHistRegistryQA->fill(HIST(MCentClasses[cBin]) + HIST(SubDir[mode]) + HIST("histDCAz"),
                            track.dcaZ());
                            */
    }

    if (mDebugLog) {
      LOGF(info, "The ThisTrackQA has been filled for cBin = %d and mode = %d.", cBin, mode);
    }
  }

 private:
  o2::framework::HistogramRegistry* mHistRegistryQA = nullptr; ///< For the QA output.

  bool mDebugLog = false;          ///< Enable to print additional log for debug.
  bool mObtainNUA = false;         ///< Enable to get the 3D Zvtx-eta-phi distribution for NUA.
  bool mSaveAllQA = false;         ///< Save the additional QA (true for QA task).
  bool mSaveQABefore = false;      ///< Save the QA output before any selection.
  bool mUseVariablePtBins = false; ///< Enable the use of a variable width pT binning.

  static const int mNcentBins = 10;                    ///< Number of centrality classes.
  static constexpr std::string_view MCentClasses[] = { ///< Centrality classes.
    "Centrality_00-01/", "Centrality_01-02/", "Centrality_02-05/",
    "Centrality_05-10/", "Centrality_10-20/", "Centrality_20-30/",
    "Centrality_30-40/", "Centrality_40-50/", "Centrality_50-60/",
    "Centrality_60-70/"};

  ClassDefNV(FlowJHistManager, 1);
};

#endif // PWGCF_JCORRAN_CORE_FLOWJHISTMANAGER_H_
