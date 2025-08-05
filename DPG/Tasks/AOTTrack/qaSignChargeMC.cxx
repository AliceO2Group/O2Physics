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
/// \file   qaSignChargeMC.cxx
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
/// \brief  Task to analyse the sign and charge of MC particles
/// \since  08/05/2024
///

// O2 includes
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/O2DatabasePDGPlugin.h"

using namespace o2::framework;

struct QaSignChargeMC {
  // Histograms
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  Configurable<int> PDG{"PDG", 2212, "PDG code of the particle of interest"};
  Configurable<bool> absPDG{"absPDG", true, "Check the PDG code in abs. value"};
  Configurable<bool> noFakeHits{"noFakeHits", true, "Check the PDG code in abs. value"};
  Configurable<bool> selPrimaries{"selPrimaries", true, "Select primaries"};
  Configurable<int> trdSel{"trdSel", 0, "TRD selection: -1 = no TRD, 0 = all, 1 = TRD only"};
  Configurable<int> tofSel{"tofSel", 0, "TOF selection: -1 = no TOF, 0 = all, 1 = TOF only"};

  Service<o2::framework::O2DatabasePDG> pdg;

  void init(InitContext&)
  {
    const AxisSpec axisReco{5, -2.5f, 2.5f, "Track sign"};
    const AxisSpec axisCharge{61, -60.5f, 60.5f, "Particle charge"};
    const AxisSpec axisPt{100, 0.f, 5.f, "#it{p}_{T} (GeV/#it{c})"};
    histos.add("sign_charge", "sign_charge", HistType::kTH2F, {axisReco, axisCharge});
    histos.add("sign_charge_pt", "sign_charge_pt", HistType::kTH3F, {axisReco, axisCharge, axisPt});
  }

  template <typename ParticleType>
  double getCharge(ParticleType const& particle)
  {
    auto pdgParticle = pdg->GetParticle(particle.pdgCode());
    if (!pdgParticle) {
      LOG(warning) << "PDG code not found: " << particle.pdgCode() << " returning 0.f as default charge.";
      return 0.f;
    }
    return pdgParticle->Charge();
  }

  void process(o2::soa::Join<o2::aod::TracksIU, o2::aod::TracksExtra, o2::aod::McTrackLabels> const& tracks,
               o2::aod::McParticles const&)
  {
    for (auto const& track : tracks) {
      if (std::abs(track.eta()) > 0.8) {
        continue;
      }
      if (!track.hasITS()) {
        continue;
      }
      if (!track.hasTPC()) {
        continue;
      }
      if (!track.has_mcParticle()) {
        continue;
      }
      if (absPDG) {
        if (std::abs(track.mcParticle().pdgCode()) != PDG) {
          continue;
        }
      } else {
        if (track.mcParticle().pdgCode() != PDG) {
          continue;
        }
      }
      if (selPrimaries && !track.mcParticle().isPhysicalPrimary()) {
        continue;
      }
      switch (trdSel) {
        case 0:
          break;
        case -1:
          if (track.hasTRD()) {
            continue;
          }
          break;
        case 1:
          if (!track.hasTRD()) {
            continue;
          }
          break;
        default:
          LOG(fatal) << "Invalid TRD selection: " << trdSel.value;
          break;
      }
      switch (tofSel) {
        case 0:
          break;
        case -1:
          if (track.hasTOF()) {
            continue;
          }
          break;
        case 1:
          if (!track.hasTOF()) {
            continue;
          }
          break;
        default:
          LOG(fatal) << "Invalid TOF selection: " << tofSel.value;
          break;
      }

      if (noFakeHits) { // Selecting tracks with no fake hits
        bool hasFakeHit = false;
        for (int i = 0; i < 10; i++) { // From ITS to TPC
          if (track.mcMask() & 1 << i) {
            hasFakeHit = true;
            break;
          }
        }
        if (hasFakeHit) {
          continue;
        }
      }
      histos.fill(HIST("sign_charge"), track.sign(), getCharge(track.mcParticle()));
      histos.fill(HIST("sign_charge_pt"), track.sign(), getCharge(track.mcParticle()), track.pt());
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<QaSignChargeMC>(cfgc)}; }
