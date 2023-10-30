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

/// \file candidateCreatorDstar.cxx
/// \brief Reconstruction of D* decay candidates
///
/// \author Vít Kučera <vit.kucera@cern.ch>, CERN

#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include "PWGHF/DataModel/CandidateReconstructionTables.h"

using namespace o2;
using namespace o2::framework;

/// Reconstruction of D* decay candidates
struct HfCandidateCreatorDstar {
  Configurable<bool> fillHistograms{"fillHistograms", true, "fill histograms"};

  OutputObj<TH1F> hMass{TH1F("hMass", "D* candidates;inv. mass (#pi D^{0}) (GeV/#it{c}^{2});entries", 500, 0., 5.)};
  OutputObj<TH1F> hPtPi{TH1F("hPtPi", "#pi candidates;#it{p}_{T} (GeV/#it{c});entries", 500, 0., 5.)};
  OutputObj<TH1F> hPtD0Prong0{TH1F("hPtD0Prong0", "D^{0} candidates;prong 0 #it{p}_{T} (GeV/#it{c});entries", 500, 0., 5.)};
  OutputObj<TH1F> hPtD0Prong1{TH1F("hPtD0Prong1", "D^{0} candidates;prong 1 #it{p}_{T} (GeV/#it{c});entries", 500, 0., 5.)};
  OutputObj<TH1F> hPtD0{TH1F("hPtD0", "D^{0} candidates;candidate #it{p}_{T} (GeV/#it{c});entries", 500, 0., 5.)};

  void process(aod::Collisions const&,
               aod::HfDstars const& rowsTrackIndexDstar,
               aod::Tracks const&,
               aod::Hf2Prongs const&)
  {
    auto massPi = o2::analysis::pdg::MassPiPlus;
    auto massD0 = o2::analysis::pdg::MassD0;

    // loop over pairs of prong indices
    for (const auto& rowTrackIndexDstar : rowsTrackIndexDstar) {
      auto trackPi = rowTrackIndexDstar.prong0();
      auto prongD0 = rowTrackIndexDstar.prongD0_as<aod::Hf2Prongs>();
      auto trackD0Prong0 = prongD0.prong0();
      auto trackD0Prong1 = prongD0.prong1();
      // auto collisionPiId = trackPi.collisionId();
      // auto collisionD0Id = trackD0Prong0.collisionId();

      // LOGF(info, "Pi collision %ld, D0 collision %ld", collisionPiId, collisionD0Id);

      std::array<float, 3> pVecPi = {trackPi.px(), trackPi.py(), trackPi.pz()};
      std::array<float, 3> pVecD0Prong0 = {trackD0Prong0.px(), trackD0Prong0.py(), trackD0Prong0.pz()};
      std::array<float, 3> pVecD0Prong1 = {trackD0Prong1.px(), trackD0Prong1.py(), trackD0Prong1.pz()};
      auto pVecD0 = RecoDecay::pVec(pVecD0Prong0, pVecD0Prong1);

      // fill histograms
      if (fillHistograms) {
        hPtPi->Fill(RecoDecay::pt(pVecPi));
        hPtD0->Fill(RecoDecay::pt(pVecD0));
        hPtD0Prong0->Fill(RecoDecay::pt(pVecD0Prong0));
        hPtD0Prong1->Fill(RecoDecay::pt(pVecD0Prong1));
        // calculate invariant mass
        auto mass = RecoDecay::m(std::array{pVecPi, pVecD0}, std::array{massPi, massD0});
        hMass->Fill(mass);
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfCandidateCreatorDstar>(cfgc)};
}
