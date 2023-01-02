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
//
// 2-particle correlations for Xi's task
// =============================
//
// Author: Rik Spijkers (rik.spijkers@cern.ch)
//

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/PIDResponse.h"

#include <TFile.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TLorentzVector.h>
#include <Math/Vector4D.h>
#include <TPDGCode.h>
#include <TDatabasePDG.h>
#include <cmath>
#include <array>
#include <cstdlib>
#include "Framework/ASoAHelpers.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;

//use parameters + cov mat non-propagated, aux info + (extension propagated)
using FullTracksExt = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCov, aod::TracksDCA>;
using FullTracksExtIU = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::TracksDCA>;
using FullTracksExtWithPID = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCov, aod::TracksDCA, aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr>;
using FullTracksExtIUWithPID = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::TracksDCA, aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr>;

// Let's try to add a column to the cascdataext table: IsSelected. 
// In the future this will probably also be a way to tell Xi's and Omega's apart
namespace o2::aod
{
namespace cascadeflags
{
DECLARE_SOA_COLUMN(IsSelected, isSelected, int); //~!
} // namespace cascadeflags
DECLARE_SOA_TABLE(CascadeFlags, "AOD", "CASCADEFLAGS", //!
                  cascadeflags::IsSelected);
using CascDataExtSelected = soa::Join<CascDataExt, CascadeFlags>;
}

// Should we make a separate struct with "cascadeSelector" here?
// The process function of this would then produce the CascadeFlags table. 
struct cascadeSelector {
  Produces<aod::CascadeFlags> cascflags;

  // histo's
  HistogramRegistry registry{
    "registry",
    {
      {"hMassXiMinus", "hMassXiMinus", {HistType::kTH1F, {{3000, 0.0f, 3.0f, "Inv. Mass (GeV/c^{2})"}}}},
    },
  };

  void process(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, aod::CascDataExt const& Cascades, aod::V0sLinked const&, aod::V0Datas const&, FullTracksExtIUWithPID const&){
    for (auto& casc : Cascades){
      auto v0 = casc.v0_as<o2::aod::V0sLinked>();
      if (!(v0.has_v0Data())) {
        cascflags(0);
        continue; // reject if no v0data
      }
      auto v0data = v0.v0Data();
      
      // Let's try to do some PID
      // these are the tracks:
      auto bachTrack = casc.bachelor_as<FullTracksExtIUWithPID>();
      auto posTrack = v0data.posTrack_as<FullTracksExtIUWithPID>();
      auto negTrack = v0data.negTrack_as<FullTracksExtIUWithPID>();

      //Bachelor check: Regardless of sign, should be pion (Xi) or kaon (Omega TODO)
      if (TMath::Abs(bachTrack.tpcNSigmaPi()) > 3) {
        cascflags(0);
        continue;
      }

      if (casc.sign() < 0){ // FIXME: only Xi for now, implement Omega's TODO
        //Proton check: 
        if (TMath::Abs(posTrack.tpcNSigmaPr()) > 3) {
          cascflags(0);
          continue;
        }
        //Pion check:
        if (TMath::Abs(negTrack.tpcNSigmaPi()) > 3) {
          cascflags(0);
          continue;
        }
      } // PID checks
      // if we reach here, candidate is good!
      cascflags(1);
    } // cascade loop
  } // process
}; // struct


struct cascadeCorrelations {
  // Basic checks
  HistogramRegistry registry{
    "registry",
    {
      {"hMassXiMinus", "hMassXiMinus", {HistType::kTH1F, {{3000, 0.0f, 3.0f, "Inv. Mass (GeV/c^{2})"}}}},
      {"hMassXiPlus", "hMassXiPlus", {HistType::kTH1F, {{3000, 0.0f, 3.0f, "Inv. Mass (GeV/c^{2})"}}}},
      {"hMassOmegaMinus", "hMassOmegaMinus", {HistType::kTH1F, {{3000, 0.0f, 3.0f, "Inv. Mass (GeV/c^{2})"}}}},
      {"hMassOmegaPlus", "hMassOmegaPlus", {HistType::kTH1F, {{3000, 0.0f, 3.0f, "Inv. Mass (GeV/c^{2})"}}}},

      {"hPhi", "hPhi", {HistType::kTH1F, {{100, 0, 2*PI, "#varphi"}}}},

      {"hDeltaPhiSS", "hDeltaPhiSS", {HistType::kTH1F, {{100, -PI/2, 1.5*PI, "#Delta#varphi"}}}},
      {"hDeltaPhiOS", "hDeltaPhiOS", {HistType::kTH1F, {{100, -PI/2, 1.5*PI, "#Delta#varphi"}}}},
    },
  };

  Filter Selector = aod::cascadeflags::isSelected > 0;

  void process(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, soa::Filtered<aod::CascDataExtSelected> const& Cascades, aod::V0sLinked const&, aod::V0Datas const&, FullTracksExtIU const&)
  {
    // some QA on the cascades
    for (auto& casc : Cascades) {

      auto v0 = casc.v0_as<o2::aod::V0sLinked>();
      if (!(v0.has_v0Data())) {
        continue; // reject if no v0data
      }
      auto v0data = v0.v0Data();
      
      if (casc.sign() < 0) { // FIXME: could be done better...
        // Check if we don't use the same pion twice in single cascade reco
        if (casc.bachelorId() == v0data.negTrackId()) LOGF(info, "autocorrelation in neg cascade! %d %d", casc.bachelorId(), v0data.negTrackId());
        registry.fill(HIST("hMassXiMinus"), casc.mXi());
        registry.fill(HIST("hMassOmegaMinus"), casc.mOmega());
      } else {
        // Check if we don't use the same pion twice in single cascade reco
        if (casc.bachelorId() == v0data.posTrackId()) LOGF(info, "autocorrelation in pos cascade! %d %d", casc.bachelorId(), v0data.negTrackId());
        registry.fill(HIST("hMassXiPlus"), casc.mXi());
        registry.fill(HIST("hMassOmegaPlus"), casc.mOmega());
      }

      registry.fill(HIST("hPhi"), RecoDecay::phi(casc.px(), casc.py()));
    } // casc loop

    for (auto& [c0, c1] : combinations(Cascades, Cascades)) { // combinations automatically applies strictly upper in case of 2 identical tables
      auto lambda0 = c0.v0_as<o2::aod::V0sLinked>();
      auto lambda1 = c1.v0_as<o2::aod::V0sLinked>();
        if (!(lambda0.has_v0Data()) || !(lambda1.has_v0Data())) {
          continue; // reject if no v0data in either of the lambda's
        }
      auto v0data0 = lambda0.v0Data();
      auto v0data1 = lambda1.v0Data();

      LOGF(info, "Found a cascade pair!"); // casc table doesn't have global indices, makes no sense to print them (will just be numbered like 0,1,2,...)
      double phi0 = RecoDecay::phi(c0.px(), c0.py());
      double phi1 = RecoDecay::phi(c1.px(), c1.py());
      double dphi = std::fmod(phi0 - phi1 + 2.5*PI, 2*PI) - 0.5*PI;
      if(c0.sign()*c1.sign() < 0){ // OS
        registry.fill(HIST("hDeltaPhiOS"), dphi);
      } else { // SS
        // Let's see if we have many autocorrelations?
        // Should only be prevalent in SS due to pions, lambda's
        // Let's first check if the lambda's aren't the same:
        if(v0data0.v0Id() == v0data1.v0Id()) {
          LOGF(info, "same v0 in SS correlation! %d %d", v0data0.v0Id(), v0data1.v0Id());
          continue;
        }
        int BachId0 = c0.bachelorId();
        int BachId1 = c1.bachelorId();
        int PosId0 = v0data0.posTrackId();
        int NegId0 = v0data0.negTrackId();
        int PosId1 = v0data1.posTrackId();
        int NegId1 = v0data1.negTrackId();
        if(BachId0 == BachId1){
          LOGF(info, "same bachelor in SS correlation! %d %d", BachId0, BachId1);
          continue;
        }
        // check for same tracks in v0's of cascades
        if (NegId0 == NegId1 || PosId0 == PosId1) {
          LOGF(info, "cascades have a v0-track in common in SS correlation!");
          continue;
        }
        if (c0.sign() < 0){ // min cascade
          if (NegId0 == BachId1 || NegId1 == BachId0) {
            LOGF(info, "bach of casc == v0-pion of other casc in neg SS correlation!");
            continue;
          }
        } else { // pos cascade
          if (PosId0 == BachId1 || PosId1 == BachId0) {
            LOGF(info, "bach of casc == v0-pion of other casc in pos SS correlation!");
            continue;
          }
        }
        registry.fill(HIST("hDeltaPhiSS"), dphi);
      }
    } // correlations
  } // process
}; //struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<cascadeSelector>(cfgc),
    adaptAnalysisTask<cascadeCorrelations>(cfgc)};
}
