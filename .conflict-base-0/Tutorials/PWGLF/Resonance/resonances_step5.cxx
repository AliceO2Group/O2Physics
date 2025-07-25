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
/// \brief this is a starting point for the Resonances tutorial
/// \author Hirak Kumar Koley
/// \since 11/10/2024

#include "PWGLF/DataModel/LFResonanceTables.h"

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// STEP 5
// loop over all MC Tracks and produce basic QA histograms
struct resonances_tutorial {
  // Define slice per Resocollision
  SliceCache cache;
  Preslice<aod::ResoTracks> perResoCollision = aod::resodaughter::resoCollisionId;
  Preslice<aod::Tracks> perCollision = aod::track::collisionId;

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // Configurable for number of bins
  Configurable<int> nBins{"nBins", 100, "N bins in all histos"};
  // Configurable for min pT cut
  Configurable<double> cMinPtcut{"cMinPtcut", 0.15, "Track minium pt cut"};

  // Track selection
  // primary track condition
  Configurable<bool> cfgPrimaryTrack{"cfgPrimaryTrack", true, "Primary track selection"};                    // kGoldenChi2 | kDCAxy | kDCAz
  Configurable<bool> cfgGlobalWoDCATrack{"cfgGlobalWoDCATrack", true, "Global track selection without DCA"}; // kQualityTracks (kTrackType | kTPCNCls | kTPCCrossedRows | kTPCCrossedRowsOverNCls | kTPCChi2NDF | kTPCRefit | kITSNCls | kITSChi2NDF | kITSRefit | kITSHits) | kInAcceptanceTracks (kPtRange | kEtaRange)
  Configurable<bool> cfgPVContributor{"cfgPVContributor", true, "PV contributor track selection"};           // PV Contriuibutor

  // DCA Selections
  // DCAr to PV
  Configurable<double> cMaxDCArToPVcut{"cMaxDCArToPVcut", 0.5, "Track DCAr cut to PV Maximum"};
  // DCAz to PV
  Configurable<double> cMaxDCAzToPVcut{"cMaxDCAzToPVcut", 2.0, "Track DCAz cut to PV Maximum"};

  // PID selection
  Configurable<float> nsigmaCutTPC{"nsigmacutTPC", 3.0, "Value of the TPC Nsigma cut"};
  Configurable<float> nsigmacutTOF{"nsigmacutTOF", 3.0, "Value /home/hirak/Desktop/LstarTaskTest/run/dpl-config.jsonof the TOF Nsigma cut"};

  // Initialize the ananlysis task
  void init(o2::framework::InitContext&)
  {
    // register histograms
    histos.add("hVertexZ", "hVertexZ", HistType::kTH1F, {{nBins, -15., 15.}});
    histos.add("hEta", "Eta distribution", kTH1F, {{200, -1.0f, 1.0f}});
    histos.add("hMultiplicityPercent", "Multiplicity Percentile", kTH1F, {{120, 0.0f, 120.0f}});
    histos.add("hDcaxy", "Dcaxy distribution", kTH1F, {{100, -0.5f, 0.5f}});
    histos.add("hDcaz", "Dcaz distribution", kTH1F, {{100, -0.5f, 0.5f}});
    histos.add("hNsigmaKaonTPC", "NsigmaKaon TPC distribution", kTH1F, {{100, -10.0f, 10.0f}});
    histos.add("hNsigmaKaonTOF", "NsigmaKaon TOF distribution", kTH1F, {{100, -10.0f, 10.0f}});

    // Print output histograms statistics
    LOG(info) << "Size of the histograms in resonance tutorial step2:";
    histos.print();
  }

  // Track selection
  template <typename TrackType>
  bool trackCut(const TrackType track)
  {
    // basic track cuts
    if (std::abs(track.pt()) < cMinPtcut)
      return false;
    if (std::abs(track.dcaXY()) > cMaxDCArToPVcut)
      return false;
    if (std::abs(track.dcaZ()) > cMaxDCAzToPVcut)
      return false;
    if (cfgPrimaryTrack && !track.isPrimaryTrack())
      return false;
    if (cfgGlobalWoDCATrack && !track.isGlobalTrackWoDCA())
      return false;
    if (cfgPVContributor && !track.isPVContributor())
      return false;
    return true;
  }

  // PID selection TPC +TOF Veto
  template <typename T>
  bool selectionPID(const T& candidate)
  {
    bool tpcPass = std::abs(candidate.tpcNSigmaKa()) < nsigmaCutTPC;
    bool tofPass = (candidate.hasTOF()) ? std::abs(candidate.tofNSigmaKa()) < nsigmacutTOF : true;
    if (tpcPass && tofPass) {
      return true;
    }
    return false;
  }

  template <typename CollisionType, typename TracksType>
  void fillHistograms(const CollisionType& /* collision */, const TracksType& dTracks)
  {
    for (auto track : dTracks) { // loop over all dTracks
      if (!trackCut(track) || !selectionPID(track)) {
        continue; // track selection and PID selection
      }
      // QA plots
      histos.fill(HIST("hEta"), track.eta());
      histos.fill(HIST("hDcaxy"), track.dcaXY());
      histos.fill(HIST("hDcaz"), track.dcaZ());
      histos.fill(HIST("hNsigmaKaonTPC"), track.tpcNSigmaKa());
      if (track.hasTOF()) {
        histos.fill(HIST("hNsigmaKaonTOF"), track.tofNSigmaKa());
      }
    }
  }

  // Process the data
  void process(aod::ResoCollision& collision,
               soa::Join<aod::ResoTracks, aod::ResoMCTracks> const& resotracks)
  {
    // Fill the event counter
    histos.fill(HIST("hVertexZ"), collision.posZ());
    histos.fill(HIST("hMultiplicityPercent"), collision.cent());

    fillHistograms(collision, resotracks);
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<resonances_tutorial>(cfgc)}; }
