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
/// \author
/// \since 23/04/2023

#include <TLorentzVector.h>

#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/runDataProcessing.h"
#include "PWGLF/DataModel/LFResonanceTables.h"
#include "DataFormatsParameters/GRPObject.h"
#include "ReconstructionDataFormats/PID.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// STEP 0
// Starting point: loop over all tracks, produce combinations and fill invariant mass histogram
struct resonances_tutorial {

  // Configurable for number of bins
  Configurable<int> nBins{"nBins", 100, "N bins in all histos"};

  // histogram defined with HistogramRegistry
  HistogramRegistry registry{"registry",
                             {{"hVertexZ", "hVertexZ", {HistType::kTH1F, {{nBins, -15., 15.}}}},
                              {"InputTracks", "InputTracks", {HistType::kTH1F, {{nBins, -5., 5.}}}},
                              {"SignTrk1", "SignTrk1", {HistType::kTH1F, {{2, -2., 2.}}}},
                              {"SignTrk2", "SignTrk2", {HistType::kTH1F, {{2, -2., 2.}}}},
                              {"hTrk1Pt", "hTrk1Pt", {HistType::kTH1F, {{nBins, 0., 5.}}}},
                              {"hTrk2Pt", "hTrk2Pt", {HistType::kTH1F, {{nBins, 0., 5.}}}},
                              {"hMassPhi", "hMassPhi", {HistType::kTH1F, {{200, 0.9f, 1.1f}}}},
                              {"hMassPhiTrue", "hMassPhiTrue", {HistType::kTH1F, {{200, 0.9f, 1.1f}}}},
                              {"hRecoPhiPtTrue", "hRecoPhiPtTrue", {HistType::kTH1F, {{nBins, 0., 5.}}}}}};

  float massKa = o2::track::PID::getMass(o2::track::PID::Kaon);

  // Defining filters for events (event selection)
  // Processed events will be already fulfulling the event selection requirements
  Filter eventFilter = (o2::aod::evsel::sel8 == true);

  // Processed collisions will be already fulfulling the position of the vertex along the z axis
  Filter vtxFilter = (nabs(o2::aod::collision::posZ) < 10.f);

  // Processed tracks will be already fulfulling the quality cuts
  Filter trackFilterEta = (nabs(aod::track::eta) < 0.8f);
  Filter trackFilterITS = (aod::track::itsChi2NCl < 36.f);
  Filter trackFilterTPC = ((aod::track::tpcChi2NCl < 4.f));

  TLorentzVector lDecayDaughter1, lDecayDaughter2, lResonance;
  void process(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& collision,
               soa::Join<aod::ResoTracks, aod::ResoMCTracks> const& resotracks,
               aod::McParticles const& mcParticles)
  {
    // Fill the event counter
    registry.fill(HIST("hVertexZ"), collision.posZ());
    for (const auto& trk : resotracks) {
      registry.fill(HIST("InputTracks"), trk.pt() * trk.sign());
    }

    for (auto& [trk1, trk2] : combinations(o2::soa::CombinationsUpperIndexPolicy(resotracks, resotracks))) {
      registry.fill(HIST("SignTrk1"), trk1.sign());
      registry.fill(HIST("SignTrk2"), trk2.sign());
      // Un-like sign pair only
      if (trk1.sign() * trk2.sign() > 0) {
        continue;
      }

      if (std::abs(trk1.tpcNSigmaKa()) > 2.0) {
        continue;
      }

      if (std::abs(trk2.tpcNSigmaKa()) > 2.0) {
        continue;
      }

      registry.fill(HIST("hTrk1Pt"), trk1.pt());
      registry.fill(HIST("hTrk2Pt"), trk2.pt());

      lDecayDaughter1.SetXYZM(trk1.px(), trk1.py(), trk1.pz(), massKa);
      lDecayDaughter2.SetXYZM(trk2.px(), trk2.py(), trk2.pz(), massKa);
      lResonance = lDecayDaughter1 + lDecayDaughter2;

      if (lResonance.Rapidity() > 0.5 || lResonance.Rapidity() < -0.5) {
        continue;
      }

      registry.fill(HIST("hMassPhi"), lResonance.M());

      // Checking the MC information
      if (abs(trk1.pdgCode()) != kKPlus || abs(trk2.pdgCode()) != kKPlus) // check if the tracks are kaons
      {
        continue;
      }
      auto mother1 = trk1.motherId();
      auto mother2 = trk2.motherId();
      if (mother1 == mother2) {        // Same mother
        if (trk1.motherPDG() == 333) { // Phi
          registry.fill(HIST("hMassPhiTrue"), lResonance.M());
          registry.fill(HIST("hRecoPhiPtTrue"), lResonance.Pt());
        }
      }
    }

    for (auto& part : mcParticles) {             // loop over all MC particles
      if (abs(part.pdgCode()) == 333) {          // Phi
        if (part.y() > 0.5 || part.y() < -0.5) { // rapidity cut
          continue;
        }
        bool isDecaytoKaons = true;
        for (auto& dau : part.daughters_as<aod::McParticles>()) {
          if (abs(dau.pdgCode()) != kKPlus) { // Decay to Kaons
            isDecaytoKaons = false;
            break;
          }
        }
        if (!isDecaytoKaons)
          continue;
        registry.fill(HIST("truephipt"), part.pt());
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<resonances_tutorial>(cfgc)}; }
