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
/// \brief
/// \author Josué Martínez García, josuem@cern.ch

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "PWGUD/DataModel/UDTables.h"
#include "PWGUD/Core/UPCTauCentralBarrelHelperRL.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct upcPhotonuclearAnalysisJMG {

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // Declare configurables on events/collisions
  Configurable<float> cutMyPosZMin{"cutMyPosZMin", -10., {"My collision cut"}};
  Configurable<float> cutMyPosZMax{"cutMyPosZMax", 10., {"My collision cut"}};
  Configurable<float> cutMyTimeZNA{"cutMyTimeZNA", 2., {"My collision cut"}};
  Configurable<float> cutMyTimeZNC{"cutMyTimeZNC", 2., {"My collision cut"}};
  // Declare configurables on side A gap
  Configurable<float> cutAGapMyEnergyZNA{"cutAGapMyEnergyZNA", 4., {"My collision cut. A Gap"}};
  Configurable<float> cutAGapMyEnergyZNC{"cutAGapMyEnergyZNC", 4., {"My collision cut. A Gap"}};
  // Declare configurables on side C gap
  Configurable<float> cutCGapMyEnergyZNA{"cutCGapMyEnergyZNA", 4., {"My collision cut. C Gap"}};
  Configurable<float> cutCGapMyEnergyZNC{"cutCGapMyEnergyZNC", 4., {"My collision cut. C Gap"}};
  // Declare configurables on both side gap
  Configurable<float> cutBothGapMyEnergyZNA{"cutBothGapMyEnergyZNA", 4., {"My collision cut. Both Gap"}};
  Configurable<float> cutBothGapMyEnergyZNC{"cutBothGapMyEnergyZNC", 4., {"My collision cut. Both Gap"}};
  // Declare configurables on tracks
  Configurable<float> cutMyptMin{"cutMyptMin", 0.15, {"My Track cut"}};
  Configurable<float> cutMyptMax{"cutMyptMax", 10., {"My Track cut"}};
  Configurable<float> cutMyetaMin{"cutMyetaMin", -0.9, {"My Track cut"}};
  Configurable<float> cutMyetaMax{"cutMyetaMax", 0.9, {"My Track cut"}};
  Configurable<float> cutMydcaZmax{"cutMydcaZmax", 2.f, {"My Track cut"}};
  Configurable<float> cutMydcaXYmax{"cutMydcaXYmax", 1e0f, {"My Track cut"}};
  Configurable<bool> cutMydcaXYusePt{"cutMydcaXYusePt", false, {"My Track cut"}};
  Configurable<bool> cutMyHasITS{"cutMyHasITS", true, {"My Track cut"}};
  Configurable<int> cutMyITSNClsMin{"cutMyITSNClsMin", 1, {"My Track cut"}};
  Configurable<float> cutMyITSChi2NClMax{"cutMyITSChi2NClMax", 36.f, {"My Track cut"}};
  Configurable<bool> cutMyHasTPC{"cutMyHasTPC", true, {"MyGlobalTrack cut"}};
  Configurable<int> cutMyTPCNClsMin{"cutMyTPCNClsMin", 0, {"My Track cut"}};
  Configurable<int> cutMyTPCNClsCrossedRowsMin{"cutMyTPCNClsCrossedRowsMin", 50, {"My Track cut"}};
  Configurable<float> cutMyTPCNClsCrossedRowsOverNClsMin{"cutMyTPCNClsCrossedRowsOverNClsMin", 1.f, {"My Track cut"}};
  Configurable<float> cutMyTPCChi2NclMax{"cutMyTPCChi2NclMax", 4.f, {"My Track cut"}};

  using FullSGUDCollision = soa::Join<aod::UDCollisions, aod::UDCollisionsSels, aod::SGCollisions, aod::UDZdcsReduced>::iterator;
  using FullUDTracks = soa::Join<aod::UDTracks, aod::UDTracksExtra, aod::UDTracksDCA, aod::UDTracksFlags>;

  void init(InitContext const&)
  {
    const AxisSpec axisCollision{4, -0.5, 3.5};
    const AxisSpec axisZvtx{40, -20., 20.};
    const AxisSpec axisPt{402, -0.05, 20.05};
    const AxisSpec axisTPCSignal{402, -0.05, 200.05};
    const AxisSpec axisPhi{64, -2 * o2::constants::math::PI, 2 * o2::constants::math::PI};
    const AxisSpec axisEta{50, -1.2, 1.2};
    const AxisSpec axisRapidity{40, -10, 10};
    const AxisSpec axisNch{101, -0.5, 100.5};
    const AxisSpec axisZNEnergy{1002, -0.5, 500.5};
    const AxisSpec axisZNTime{21, -10.5, 10.5};
    const AxisSpec axisNCls{201, -0.5, 200.5};
    const AxisSpec axisChi2NCls{100, 0, 50};
    const AxisSpec axisTPCNClsCrossedRowsMin{100, -0.05, 2.05};

    histos.add("Events/hCountCollisions", "0 total - 1 side A - 2 side C - 3 both side; Number of analysed collision; counts", kTH1F, {axisCollision});

    // histos to selection gap in side A
    histos.add("Tracks/SGsideA/hTrackPt", "#it{p_{T}} distribution; #it{p_{T}}; counts", kTH1F, {axisPt});
    histos.add("Tracks/SGsideA/hTrackPhi", "#it{#phi} distribution; #it{#phi}; counts", kTH1F, {axisPhi});
    histos.add("Tracks/SGsideA/hTrackEta", "#it{#eta} distribution; #it{#eta}; counts", kTH1F, {axisEta});
    histos.add("Tracks/SGsideA/hTrackRapidity", "#it{y} distribution; #it{y}; counts", kTH1F, {axisRapidity});
    histos.add("Tracks/SGsideA/hTrackITSNCls", "#it{N Clusters ITS} distribution; #it{N Clusters ITS}; counts", kTH1F, {axisNCls});
    histos.add("Tracks/SGsideA/hTrackITSChi2NCls", "#it{N Clusters Chi2 ITS} distribution; #it{N Clusters Chi2 ITS}; counts", kTH1F, {axisChi2NCls});
    histos.add("Tracks/SGsideA/hTrackTPCNCls", "#it{N Clusters TPC} distribution; #it{N Clusters TPC}; counts", kTH1F, {axisTPCNClsCrossedRowsMin});
    histos.add("Tracks/SGsideA/hTrackTPCNClsCrossedRows", "#it{Number of crossed TPC Rows} distribution; #it{Number of crossed TPC Rows}; counts", kTH1F, {axisNCls});
    histos.add("Tracks/SGsideA/hTrackTPCNClsFindable", "#it{Findable TPC clusters for this track} distribution; #it{Findable TPC clusters for this track}; counts", kTH1F, {axisNCls});
    histos.add("Tracks/SGsideA/hTrackTPCChi2NCls", "#it{N Clusters Chi2 TPC} distribution; #it{N Clusters Chi2 TPC}; counts", kTH1F, {axisChi2NCls});
    histos.add("Tracks/SGsideA/hTrackITSNClsTPCCls", "#it{ITS Clusters vs TPC Clusters}; #it{TPC Clusters}; #it{ITS Clusters}", kTH2F, {axisNCls, axisNCls});
    histos.add("Tracks/SGsideA/hTrackTPCSignnalP", "#it{TPC dE/dx vs p}; #it{TPC dE/dx}; #it{p}", kTH2F, {axisPt, axisTPCSignal});

    histos.add("Events/SGsideA/hTrackZVtx", "vertex in z; z (cm); counts", kTH1F, {axisZvtx});
    histos.add("Events/SGsideA/hNch", "#it{N_{ch}} distribution; #it{N_{ch}}; counts", kTH1F, {axisNch});
    histos.add("Events/SGsideA/hPtVSNch", "#it{ #LT p_{T} #GT } vs #it{N_{ch}}; #it{N_{ch}}; #it{ #LT p_{T} #GT }", kTH2F, {axisNch, axisPt});
    histos.add("Events/SGsideA/hEnergyZNA", "Energy in side A distribution; Energy in side A; counts", kTH1F, {axisZNEnergy});
    histos.add("Events/SGsideA/hEnergyZNC", "Energy in side C distribution; Energy in side C; counts", kTH1F, {axisZNEnergy});
    histos.add("Events/SGsideA/hEnergyRelationSides", "Energy in side A vs energy in side C; Energy in side A; Energy in side C", kTH2F, {axisZNEnergy, axisZNEnergy});
    histos.add("Events/SGsideA/hTimeZNA", "Time in side A distribution; Time in side A; counts", kTH1F, {axisZNTime});
    histos.add("Events/SGsideA/hTimeZNC", "Time in side C distribution; Time in side C; counts", kTH1F, {axisZNTime});
    histos.add("Events/SGsideA/hTimeRelationSides", "Time in side A vs time in side C; Time in side A; Time in side C", kTH2F, {axisZNTime, axisZNTime});

    // histos to selection gap in side C
    histos.add("Tracks/SGsideC/hTrackPt", "#it{p_{T}} distribution; #it{p_{T}}; counts", kTH1F, {axisPt});
    histos.add("Tracks/SGsideC/hTrackPhi", "#it{#phi} distribution; #it{#phi}; counts", kTH1F, {axisPhi});
    histos.add("Tracks/SGsideC/hTrackEta", "#it{#eta} distribution; #it{#eta}; counts", kTH1F, {axisEta});
    histos.add("Tracks/SGsideC/hTrackRapidity", "#it{y} distribution; #it{y}; counts", kTH1F, {axisRapidity});
    histos.add("Tracks/SGsideC/hTrackITSNClsTPCCls", "#it{ITS Clusters vs TPC Clusters}; #it{TPC Clusters}; #it{ITS Clusters}", kTH2F, {axisNCls, axisNCls});
    histos.add("Tracks/SGsideC/hTrackTPCSignnalP", "#it{TPC dE/dx vs p}; #it{TPC dE/dx}; #it{p}", kTH2F, {axisPt, axisTPCSignal});

    histos.add("Events/SGsideC/hTrackZVtx", "vertex in z; z (cm); counts", kTH1F, {axisZvtx});
    histos.add("Events/SGsideC/hNch", "#it{N_{ch}} distribution; #it{N_{ch}}; counts", kTH1F, {axisNch});
    histos.add("Events/SGsideC/hPtVSNch", "#it{ #LT p_{T} #GT } vs #it{N_{ch}}; #it{N_{ch}}; #it{ #LT p_{T} #GT }", kTH2F, {axisNch, axisPt});
    histos.add("Events/SGsideC/hEnergyZNA", "Energy in side A distribution; Energy in side A; counts", kTH1F, {axisZNEnergy});
    histos.add("Events/SGsideC/hEnergyZNC", "Energy in side C distribution; Energy in side C; counts", kTH1F, {axisZNEnergy});
    histos.add("Events/SGsideC/hEnergyRelationSides", "Energy in side A vs energy in side C; Energy in side A; Energy in side C", kTH2F, {axisZNEnergy, axisZNEnergy});
    histos.add("Events/SGsideC/hTimeZNA", "Time in side A distribution; Time in side A; counts", kTH1F, {axisZNTime});
    histos.add("Events/SGsideC/hTimeZNC", "Time in side C distribution; Time in side C; counts", kTH1F, {axisZNTime});
    histos.add("Events/SGsideC/hTimeRelationSides", "Time in side A vs time in side C; Time in side A; Time in side C", kTH2F, {axisZNTime, axisZNTime});

    // histos to selection gap in both sides
    histos.add("Tracks/SGsideBoth/hTrackPt", "#it{p_{T}} distribution; #it{p_{T}}; counts", kTH1F, {axisPt});
    histos.add("Tracks/SGsideBoth/hTrackPhi", "#it{#phi} distribution; #it{#phi}; counts", kTH1F, {axisPhi});
    histos.add("Tracks/SGsideBoth/hTrackEta", "#it{#eta} distribution; #it{#eta}; counts", kTH1F, {axisEta});
    histos.add("Tracks/SGsideBoth/hTrackRapidity", "#it{y} distribution; #it{y}; counts", kTH1F, {axisRapidity});
    histos.add("Tracks/SGsideBoth/hTrackITSNClsTPCCls", "#it{ITS Clusters vs TPC Clusters}; #it{TPC Clusters}; #it{ITS Clusters}", kTH2F, {axisNCls, axisNCls});
    histos.add("Tracks/SGsideBoth/hTrackTPCSignnalP", "#it{TPC dE/dx vs p}; #it{TPC dE/dx}; #it{p}", kTH2F, {axisPt, axisTPCSignal});

    histos.add("Events/SGsideBoth/hTrackZVtx", "vertex in z; z (cm); counts", kTH1F, {axisZvtx});
    histos.add("Events/SGsideBoth/hNch", "#it{N_{ch}} distribution; #it{N_{ch}}; counts", kTH1F, {axisNch});
    histos.add("Events/SGsideBoth/hPtVSNch", "#it{ #LT p_{T} #GT } vs #it{N_{ch}}; #it{N_{ch}}; #it{ #LT p_{T} #GT }", kTH2F, {axisNch, axisPt});
    histos.add("Events/SGsideBoth/hEnergyZNA", "Energy in side A distribution; Energy in side A; counts", kTH1F, {axisZNEnergy});
    histos.add("Events/SGsideBoth/hEnergyZNC", "Energy in side C distribution; Energy in side C; counts", kTH1F, {axisZNEnergy});
    histos.add("Events/SGsideBoth/hEnergyRelationSides", "Energy in side A vs energy in side C; Energy in side A; Energy in side C", kTH2F, {axisZNEnergy, axisZNEnergy});
    histos.add("Events/SGsideBoth/hTimeZNA", "Time in side A distribution; Time in side A; counts", kTH1F, {axisZNTime});
    histos.add("Events/SGsideBoth/hTimeZNC", "Time in side C distribution; Time in side C; counts", kTH1F, {axisZNTime});
    histos.add("Events/SGsideBoth/hTimeRelationSides", "Time in side A vs time in side C; Time in side A; Time in side C", kTH2F, {axisZNTime, axisZNTime});
  }

  template <typename C>
  bool isGlobalCollisionCut(C const& collision)
  {
    if (collision.posZ() < cutMyPosZMin || cutMyPosZMax < collision.posZ()) {
      return false;
    }
    if ((std::abs(collision.timeZNA()) < cutMyTimeZNA && std::abs(collision.timeZNC()) < cutMyTimeZNC) == false) {
      return false;
    }
    return true;
  }

  template <typename CSG>
  bool isCollisionCutSG(CSG const& collision, int SideGap)
  {
    switch (SideGap) {
      case 0:                                                                                                                  // Gap in A side
        if ((collision.energyCommonZNA() < cutAGapMyEnergyZNA && collision.energyCommonZNC() > cutAGapMyEnergyZNC) == false) { // 0n - A side && Xn - C Side
          return false;
        }
        break;
      case 1:                                                                                                                  // Gap in C side
        if ((collision.energyCommonZNA() > cutCGapMyEnergyZNA && collision.energyCommonZNC() < cutCGapMyEnergyZNC) == false) { // Xn - A side && 0n - C Side
          return false;
        }
        break;
      case 2:                                                                                                                        // Gap in Both Sides
        if ((collision.energyCommonZNA() < cutBothGapMyEnergyZNA && collision.energyCommonZNC() < cutBothGapMyEnergyZNC) == false) { // 0n - A side && 0n - C Side
          return false;
        }
        break;
    }
    return true;
  }

  template <typename T>
  bool isTrackCut(T const& track)
  {
    if (track.pt() < cutMyptMin || track.pt() > cutMyptMax) {
      return false;
    }
    if (eta(track.px(), track.py(), track.pz()) < cutMyetaMin || eta(track.px(), track.py(), track.pz()) > cutMyetaMax) {
      return false;
    }
    if (std::abs(track.dcaZ()) > cutMydcaZmax) {
      return false;
    }
    if (cutMydcaXYusePt) {
      float maxDCA = 0.0105f + 0.0350f / std::pow(track.pt(), 1.1f);
      if (std::abs(track.dcaXY()) > maxDCA) {
        return false;
      }
    } else {
      if (std::abs(track.dcaXY()) > cutMydcaXYmax) {
        return false;
      }
    }
    // Quality Track
    // ITS
    if (cutMyHasITS && !track.hasITS()) {
      return false; // ITS refit
    }
    if (track.itsNCls() < cutMyITSNClsMin) {
      return false;
    }
    if (track.itsChi2NCl() > cutMyITSChi2NClMax) {
      return false;
    }
    // TPC
    if (cutMyHasTPC && !track.hasTPC()) {
      return false; // TPC refit
    }
    if (track.tpcNClsCrossedRows() < cutMyTPCNClsCrossedRowsMin) {
      return false;
    }
    if ((static_cast<float>(track.tpcNClsCrossedRows()) / static_cast<float>(track.tpcNClsFindable())) < cutMyTPCNClsMin) {
      return false; // tpcNClsFound()
    }
    if (track.tpcChi2NCl() > cutMyTPCChi2NclMax) {
      return false; // TPC chi2
    }
    return true;
  }

  void processSG(FullSGUDCollision const& reconstructedCollision, FullUDTracks const& reconstructedTracks)
  {
    histos.fill(HIST("Events/hCountCollisions"), 0);
    int SGside = reconstructedCollision.gapSide();
    int nTracksCharged = 0;
    float sumPt = 0;

    if (isGlobalCollisionCut(reconstructedCollision) == false) {
      return;
    }

    switch (SGside) {
      case 0: // for side A
        if (isCollisionCutSG(reconstructedCollision, 0) == false) {
          return;
        }
        histos.fill(HIST("Events/hCountCollisions"), 1);
        histos.fill(HIST("Events/SGsideA/hEnergyZNA"), reconstructedCollision.energyCommonZNA());
        histos.fill(HIST("Events/SGsideA/hEnergyZNC"), reconstructedCollision.energyCommonZNC());
        histos.fill(HIST("Events/SGsideA/hEnergyRelationSides"), reconstructedCollision.energyCommonZNA(), reconstructedCollision.energyCommonZNC());
        histos.fill(HIST("Events/SGsideA/hTimeZNA"), reconstructedCollision.timeZNA());
        histos.fill(HIST("Events/SGsideA/hTimeZNC"), reconstructedCollision.timeZNC());
        histos.fill(HIST("Events/SGsideA/hTimeRelationSides"), reconstructedCollision.timeZNA(), reconstructedCollision.timeZNC());
        histos.fill(HIST("Events/SGsideA/hTrackZVtx"), reconstructedCollision.posZ());
        for (auto& track : reconstructedTracks) {
          if (track.sign() == 1 || track.sign() == -1) {
            if (isTrackCut(track) == false) {
              continue;
            }
            nTracksCharged++;
            sumPt += track.pt();
            histos.fill(HIST("Tracks/SGsideA/hTrackPt"), track.pt());
            histos.fill(HIST("Tracks/SGsideA/hTrackPhi"), phi(track.px(), track.py()));
            histos.fill(HIST("Tracks/SGsideA/hTrackEta"), eta(track.px(), track.py(), track.pz()));
            histos.fill(HIST("Tracks/SGsideA/hTrackTPCSignnalP"), momentum(track.px(), track.py(), track.pz()), track.tpcSignal());
            // histos.fill(HIST("Tracks/SGsideA/hTrackRapidity"), track.rapidity());

            histos.fill(HIST("Tracks/SGsideA/hTrackITSNCls"), track.itsNCls());
            histos.fill(HIST("Tracks/SGsideA/hTrackITSChi2NCls"), track.itsChi2NCl());
            histos.fill(HIST("Tracks/SGsideA/hTrackTPCNCls"), (static_cast<float>(track.tpcNClsCrossedRows()) / static_cast<float>(track.tpcNClsFindable())));
            histos.fill(HIST("Tracks/SGsideA/hTrackTPCNClsCrossedRows"), track.tpcNClsCrossedRows());
            histos.fill(HIST("Tracks/SGsideA/hTrackTPCNClsFindable"), track.tpcNClsFindable());
            histos.fill(HIST("Tracks/SGsideA/hTrackTPCChi2NCls"), track.tpcChi2NCl());
            histos.fill(HIST("Tracks/SGsideA/hTrackITSNClsTPCCls"), track.tpcNClsFindable(), track.itsNCls());
          }
        }
        histos.fill(HIST("Events/SGsideA/hNch"), nTracksCharged);
        histos.fill(HIST("Events/SGsideA/hPtVSNch"), nTracksCharged, (sumPt / nTracksCharged));
        nTracksCharged = sumPt = 0;
        break;
      case 1: // for side C
        if (isCollisionCutSG(reconstructedCollision, 1) == false) {
          return;
        }
        histos.fill(HIST("Events/hCountCollisions"), 2);
        histos.fill(HIST("Events/SGsideC/hEnergyZNA"), reconstructedCollision.energyCommonZNA());
        histos.fill(HIST("Events/SGsideC/hEnergyZNC"), reconstructedCollision.energyCommonZNC());
        histos.fill(HIST("Events/SGsideC/hEnergyRelationSides"), reconstructedCollision.energyCommonZNA(), reconstructedCollision.energyCommonZNC());
        histos.fill(HIST("Events/SGsideC/hTimeZNA"), reconstructedCollision.timeZNA());
        histos.fill(HIST("Events/SGsideC/hTimeZNC"), reconstructedCollision.timeZNC());
        histos.fill(HIST("Events/SGsideC/hTimeRelationSides"), reconstructedCollision.timeZNA(), reconstructedCollision.timeZNC());
        histos.fill(HIST("Events/SGsideC/hTrackZVtx"), reconstructedCollision.posZ());
        for (auto& track : reconstructedTracks) {
          if (track.sign() == 1 || track.sign() == -1) {
            if (isTrackCut(track) == false) {
              continue;
            }
            nTracksCharged++;
            sumPt += track.pt();
            histos.fill(HIST("Tracks/SGsideC/hTrackPt"), track.pt());
            histos.fill(HIST("Tracks/SGsideC/hTrackPhi"), phi(track.px(), track.py()));
            histos.fill(HIST("Tracks/SGsideC/hTrackEta"), eta(track.px(), track.py(), track.pz()));
            // histos.fill(HIST("Tracks/SGsideC/hTrackRapidity"), track.rapidity());
            histos.fill(HIST("Tracks/SGsideC/hTrackITSNClsTPCCls"), track.tpcNClsFindable(), track.itsNCls());
            histos.fill(HIST("Tracks/SGsideC/hTrackTPCSignnalP"), momentum(track.px(), track.py(), track.pz()), track.tpcSignal());
          }
        }
        histos.fill(HIST("Events/SGsideC/hNch"), nTracksCharged);
        histos.fill(HIST("Events/SGsideC/hPtVSNch"), nTracksCharged, (sumPt / nTracksCharged));
        nTracksCharged = sumPt = 0;
        break;
      case 2: // for both sides
        if (isCollisionCutSG(reconstructedCollision, 2) == false) {
          return;
        }
        histos.fill(HIST("Events/hCountCollisions"), 3);
        histos.fill(HIST("Events/SGsideBoth/hEnergyZNA"), reconstructedCollision.energyCommonZNA());
        histos.fill(HIST("Events/SGsideBoth/hEnergyZNC"), reconstructedCollision.energyCommonZNC());
        histos.fill(HIST("Events/SGsideBoth/hEnergyRelationSides"), reconstructedCollision.energyCommonZNA(), reconstructedCollision.energyCommonZNC());
        histos.fill(HIST("Events/SGsideBoth/hTimeZNA"), reconstructedCollision.timeZNA());
        histos.fill(HIST("Events/SGsideBoth/hTimeZNC"), reconstructedCollision.timeZNC());
        histos.fill(HIST("Events/SGsideBoth/hTimeRelationSides"), reconstructedCollision.timeZNA(), reconstructedCollision.timeZNC());
        histos.fill(HIST("Events/SGsideBoth/hTrackZVtx"), reconstructedCollision.posZ());
        for (auto& track : reconstructedTracks) {
          if (track.sign() == 1 || track.sign() == -1) {
            if (isTrackCut(track) == false) {
              continue;
            }
            nTracksCharged++;
            sumPt += track.pt();
            histos.fill(HIST("Tracks/SGsideBoth/hTrackPt"), track.pt());
            histos.fill(HIST("Tracks/SGsideBoth/hTrackPhi"), phi(track.px(), track.py()));
            histos.fill(HIST("Tracks/SGsideBoth/hTrackEta"), eta(track.px(), track.py(), track.pz()));
            // histos.fill(HIST("Tracks/SGsideBoth/hTrackRapidity"), track.rapidity());
            histos.fill(HIST("Tracks/SGsideBoth/hTrackITSNClsTPCCls"), track.tpcNClsFindable(), track.itsNCls());
            histos.fill(HIST("Tracks/SGsideBoth/hTrackTPCSignnalP"), momentum(track.px(), track.py(), track.pz()), track.tpcSignal());
          }
        }
        histos.fill(HIST("Events/SGsideBoth/hNch"), nTracksCharged);
        histos.fill(HIST("Events/SGsideBoth/hPtVSNch"), nTracksCharged, (sumPt / nTracksCharged));
        nTracksCharged = sumPt = 0;
        break;
      default:
        return;
        break;
    }
  }
  PROCESS_SWITCH(upcPhotonuclearAnalysisJMG, processSG, "Process in UD tables", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<upcPhotonuclearAnalysisJMG>(cfgc, TaskName{"upcphotonuclear"})};
}
