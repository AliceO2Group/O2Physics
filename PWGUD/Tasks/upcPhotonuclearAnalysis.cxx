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

struct upcPhotonuclearAnalysis {

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // Declare configurables on events/collisions
  Configurable<float> cutMyPosZMin{"cutMyPosZMin", -20., {"My collision cut"}};
  Configurable<float> cutMyPosZMax{"cutMyPosZMax", 20., {"My collision cut"}};
  // Declare configurables on tracks
  Configurable<float> cutMyptMin{"cutMyptMin", 0.15, {"My Track cut"}};
  Configurable<float> cutMyptMax{"cutMyptMax", 10., {"My Track cut"}};
  Configurable<float> cutMyetaMin{"cutMyetaMin", -1.0, {"My Track cut"}};
  Configurable<float> cutMyetaMax{"cutMyetaMax", 1.0, {"My Track cut"}};
  Configurable<float> cutMydcaZmax{"cutMydcaZmax", 2.f, {"My Track cut"}};
  Configurable<float> cutMydcaXYmax{"cutMydcaXYmax", 1e0f, {"My Track cut"}};
  Configurable<bool> cutMydcaXYusePt{"cutMydcaXYusePt", false, {"My Track cut"}};

  using FullSGUDCollision = soa::Join<aod::UDCollisions, aod::UDCollisionsSels, aod::SGCollisions, aod::UDZdcsReduced>::iterator;
  using FullUDTracks = soa::Join<aod::UDTracks, aod::UDTracksExtra, aod::UDTracksDCA, aod::UDTracksFlags>;

  void init(InitContext const&)
  {
    const AxisSpec axisCollision{4, -0.5, 3.5};
    const AxisSpec axisZvtx{40, -20., 20.};
    const AxisSpec axisPt{402, -0.05, 20.05};
    const AxisSpec axisPhi{64, -2 * o2::constants::math::PI, 2 * o2::constants::math::PI};
    const AxisSpec axisEta{50, -1.2, 1.2};
    const AxisSpec axisNch{101, -0.5, 100.5};
    const AxisSpec axisZNEnergy{1000, 0, 1000};
    const AxisSpec axisZNTime{40, -20, 20};

    histos.add("Events/hCountCollisions", "0 total - 1 side A - 2 side C - 3 both side; Number of analysed collision; counts", kTH1F, {axisCollision});

    // histos to selection gap in side A
    histos.add("Tracks/SGsideA/hTrackPt", "#it{p_{T}} distribution; #it{p_{T}}; counts", kTH1F, {axisPt});
    histos.add("Tracks/SGsideA/hTrackPhi", "#it{#phi} distribution; #it{#phi}; counts", kTH1F, {axisPhi});
    histos.add("Tracks/SGsideA/hTrackEta", "#it{#eta} distribution; #it{#eta}; counts", kTH1F, {axisEta});

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
  bool isCollisionCut(C const& collision){
    if ( collision.posZ() < || < collision.posZ() ) {
      return false;
    }
    if ( )
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
    return true;
  }

  void processSG(FullSGUDCollision const& reconstructedCollision, FullUDTracks const& reconstructedTracks)
  {
    histos.fill(HIST("Events/hCountCollisions"), 0);
    int SGside = reconstructedCollision.gapSide();
    int nTracksCharged = 0;
    float sumPt = 0;

    switch (SGside) {
      case 0: // for side A
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
          }
        }
        histos.fill(HIST("Events/SGsideA/hNch"), nTracksCharged);
        histos.fill(HIST("Events/SGsideA/hPtVSNch"), nTracksCharged, (sumPt / nTracksCharged));
        nTracksCharged = sumPt = 0;
        break;
      case 1: // for side C
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
          }
        }
        histos.fill(HIST("Events/SGsideC/hNch"), nTracksCharged);
        histos.fill(HIST("Events/SGsideC/hPtVSNch"), nTracksCharged, (sumPt / nTracksCharged));
        nTracksCharged = sumPt = 0;
        break;
      case 2: // for both sides
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
  PROCESS_SWITCH(upcPhotonuclearAnalysis, processSG, "Process in UD tables", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<upcPhotonuclearAnalysis>(cfgc, TaskName{"upcphotonuclear"})};
}
