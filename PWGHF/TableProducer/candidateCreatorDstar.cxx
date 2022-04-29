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
// #include "DetectorsVertexing/DCAFitterN.h"
// #include "Common/Core/trackUtilities.h"
// #include "ReconstructionDataFormats/DCA.h"
#include "PWGHF/DataModel/HFSecondaryVertex.h"

#include "Framework/runDataProcessing.h"

using namespace o2;
using namespace o2::framework;
// using namespace o2::aod::hf_cand;
// using namespace o2::aod::hf_cand_prong2;

/// Reconstruction of D* decay candidates
struct HfCandidateCreatorDstar {
  Configurable<bool> fillHistograms{"fillHistograms", true, "fill histograms"};

  OutputObj<TH1F> hMass{TH1F("hMass", "D* candidates;inv. mass (#pi D^{0}) (GeV/#it{c}^{2});entries", 500, 0., 5.)};

  double massPi = RecoDecay::getMassPDG(kPiPlus);
  double massD0 = RecoDecay::getMassPDG(pdg::Code::kD0);

  void process(aod::Collisions const&,
               aod::HfDstar const& rowsTrackIndexDstar,
               aod::BigTracks const&,
               aod::Hf2Prong const&)
  {
    // loop over pairs of prong indices
    for (const auto& rowTrackIndexDstar : rowsTrackIndexDstar) {
      auto trackPi = rowTrackIndexDstar.index0_as<aod::BigTracks>();
      auto prongD0 = rowTrackIndexDstar.indexD0_as<aod::Hf2Prong>();
      auto trackD0Prong0 = prongD0.index0_as<aod::BigTracks>();
      auto trackD0Prong1 = prongD0.index1_as<aod::BigTracks>();
      auto collisionPiId = trackPi.collisionId();
      auto collisionD0Id = trackD0Prong0.collisionId();

      LOGF(info, "Pi collision %ld, D0 collision %ld", collisionPiId, collisionD0Id);

      std::array<float, 3> pVecPi = {trackPi.px(), trackPi.py(), trackPi.pz()};
      std::array<float, 3> pVecD0Prong0 = {trackD0Prong0.px(), trackD0Prong0.py(), trackD0Prong0.pz()};
      std::array<float, 3> pVecD0Prong1 = {trackD0Prong1.px(), trackD0Prong1.py(), trackD0Prong1.pz()};
      auto pVecD0 = RecoDecay::PVec(pVecD0Prong0, pVecD0Prong1);

      // fill histograms
      if (fillHistograms) {
        // calculate invariant mass
        auto mass = RecoDecay::M(std::array{pVecPi, pVecD0}, std::array{massPi, massD0});
        hMass->Fill(mass);
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfCandidateCreatorDstar>(cfgc)
  };
}
