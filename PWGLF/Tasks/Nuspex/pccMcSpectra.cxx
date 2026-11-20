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

/// \file pccMcSpectra.cxx
/// \brief Task to extract the pt spectra of various particle species as input for a particle-composition correction.
/// \author Maria Calmon Behling <maria.calmon@cern.ch>

#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/runDataProcessing.h>

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <vector>

using namespace o2;
using namespace o2::framework;

struct pccMcSpectra {

  Service<o2::framework::O2DatabasePDG> pdg;

  Configurable<float> etaINEL{"etaINEL", 1.f, "eta range for INEL>0"};
  Configurable<float> etaNch{"etaNch", 0.5f, "eta range for Nch"};
  Configurable<float> yMax{"yMax", 0.5f, "rapidity range for spectra"};
  Configurable<uint32_t> maxMult{"maxMult", 100, "max multiplicity"};

  HistogramRegistry histos{"histos"};

  void init(InitContext const&)
  {
    // Integer multiplicity bins: Nch = 0, 1, ..., 199
    const int nBinsMult = maxMult + 1;
    const AxisSpec multAxis = {nBinsMult, -0.5, nBinsMult - 0.5, "#it{N}_{ch}", "mult"};

    // Transverse momentum: 0 - 10 GeV/c

    std::vector<double> ptBinEdges = {0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75,
                                      0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9,
                                      2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.5, 5.0, 5.5,
                                      6.0, 6.5, 7.0, 8.0, 9.0, 10.0, 20.0, 30.0};
    const AxisSpec ptAxis{ptBinEdges, "#it{p}_{T} (GeV/#it{c})", "pt"};

    histos.add("event/mult",
               "multiplicity distribution;Nch;Events",
               kTH1D, {multAxis});

    histos.add("spectra/pion",
               "Pion spectra;Nch;pT",
               kTH2D, {multAxis, ptAxis});

    histos.add("spectra/kaon",
               "Kaon spectra;Nch;pT",
               kTH2D, {multAxis, ptAxis});

    histos.add("spectra/proton",
               "Proton spectra;Nch;pT",
               kTH2D, {multAxis, ptAxis});

    histos.add("spectra/sigma",
               "Sigma spectra;Nch;pT",
               kTH2D, {multAxis, ptAxis});

    histos.add("spectra/lambda",
               "Lambda spectra;Nch;pT",
               kTH2D, {multAxis, ptAxis});

    histos.add("spectra/xi",
               "Xi spectra;Nch;pT",
               kTH2D, {multAxis, ptAxis});

    histos.add("spectra/omega",
               "Omega spectra;Nch;pT",
               kTH2D, {multAxis, ptAxis});
  }

  void process(aod::McCollisions::iterator const&,
               aod::McParticles const& particles)
  {
    int nch = 0;
    bool isINELgt0 = false;

    // First loop: determine event multiplicity and INEL>0.
    for (const auto& particle : particles) {

      if (!particle.isPhysicalPrimary()) {
        continue;
      }

      auto pdgParticle = pdg->GetParticle(particle.pdgCode());

      if (!pdgParticle || pdgParticle->Charge() == 0.) {
        continue;
      }

      float absEta = std::abs(particle.eta());

      if (absEta < etaINEL) {
        isINELgt0 = true;
      }

      if (absEta < etaNch) {
        ++nch;
      }
    }

    // Reject events that do not satisfy INEL>0.
    if (!isINELgt0) {
      return;
    }

    histos.fill(HIST("event/mult"), nch);

    // Second loop: fill identified-particle spectra.
    for (const auto& particle : particles) {

      if (!particle.isPhysicalPrimary()) {
        continue;
      }

      if (std::abs(particle.y()) >= yMax) {
        continue;
      }

      int absPDG = std::abs(particle.pdgCode());
      float pt = particle.pt();

      switch (absPDG) {

        case 211: // Pion
          histos.fill(HIST("spectra/pion"), nch, pt);
          break;

        case 321: // Kaon
          histos.fill(HIST("spectra/kaon"), nch, pt);
          break;

        case 2212: // Proton
          histos.fill(HIST("spectra/proton"), nch, pt);
          break;

        case 3222: // Sigma+
        case 3112: // Sigma-
          histos.fill(HIST("spectra/sigma"), nch, pt);
          break;

        case 3122: // Lambda
          histos.fill(HIST("spectra/lambda"), nch, pt);
          break;

        case 3312: // Xi
          histos.fill(HIST("spectra/xi"), nch, pt);
          break;

        case 3334: // Omega
          histos.fill(HIST("spectra/omega"), nch, pt);
          break;

        default:
          break;
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<pccMcSpectra>(cfgc)};
}
