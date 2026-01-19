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
/// \file alice3SeparationPower.cxx
///
/// \brief This task produces the separation power of the ALICE3 detector
///
/// \author Nicolò Jacazio, Universita del Piemonte Orientale (IT)
/// \since  May 13, 2025
///

#include "ALICE3/DataModel/OTFRICH.h"
#include "ALICE3/DataModel/OTFTOF.h"

#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/RunningWorkflowInfo.h>
#include <Framework/runDataProcessing.h>

#include <THashList.h>
#include <TProfile2D.h>

#include <map>
#include <string>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;

std::array<TProfile2D*, 5> separationInnerTOF;
std::array<TProfile2D*, 5> separationOuterTOF;
std::array<TProfile2D*, 5> separationRICH;
struct alice3SeparationPower {

  ConfigurableAxis etaAxis{"etaAxis", {100, -1.f, 1.f}, "Binning in eta"};
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  OutputObj<THashList> listSeparation{"separationPower"};
  void init(o2::framework::InitContext&)
  {
    listSeparation.setObject(new THashList);
    for (int i = 0; i < 5; i++) {
      auto createEfficiency = [&](const char* name, const char* title) {
        TProfile2D* eff = new TProfile2D(Form("%s_%d", name, i),
                                         Form("%s_%d;%s", title, i, "#it{p}_{T} (GeV/#it{c});#it{#eta}"),
                                         100, 0.f, 10.f,
                                         100, 0.f, 10.f);
        listSeparation->Add(eff);
        return eff;
      };
      separationInnerTOF[i] = createEfficiency("separationInnerTOF", "separationInnerTOF");
      separationOuterTOF[i] = createEfficiency("separationOuterTOF", "separationOuterTOF");
      separationRICH[i] = createEfficiency("separationRICH", "separationRICH");
    }
  }

  void process(soa::Join<aod::Collisions, aod::McCollisionLabels>::iterator const& /*collision*/,
               soa::Join<aod::Tracks, aod::TracksCov, aod::McTrackLabels, aod::UpgradeTofMC, aod::UpgradeTof> const& tracks,
               aod::McParticles const&,
               aod::McCollisions const&)
  {
    for (const auto& track : tracks) {
      if (!track.has_mcParticle()) {
        continue;
      }
      //   Check that all the nsigmas are numbers (sanity check)
      for (int i = 0; i < 5; i++) {
        if (std::isnan(track.nSigmaInnerTOF(i)) || std::isnan(track.nSigmaOuterTOF(i))) {
          LOG(warning) << "Unrecognized nsigma for " << i << " " << track.nSigmaInnerTOF(i) << " " << track.nSigmaOuterTOF(i);
        }
      }

      const auto& mcParticle = track.mcParticle();
      // Separation electron pion
      switch (std::abs(mcParticle.pdgCode())) {
        {
          case 211: // electron-pion separation
            separationInnerTOF[0]->Fill(track.pt(), track.eta(), track.nSigmaInnerTOF(0));
            separationOuterTOF[0]->Fill(track.pt(), track.eta(), track.nSigmaOuterTOF(0));
            // separationRICH[0]->Fill(track.pt(), track.eta(), track.nSigmaElectronRich() );
            break;
          case 321: // pion-kaon separation
            separationInnerTOF[1]->Fill(track.pt(), track.eta(), track.nSigmaInnerTOF(1));
            separationOuterTOF[1]->Fill(track.pt(), track.eta(), track.nSigmaInnerTOF(1));
            // separationRICH[1]->Fill(track.pt(), track.eta(), track.nSigmaPionRich() );
            break;
          case 2212: // kaon-proton separation
            separationInnerTOF[2]->Fill(track.pt(), track.eta(), track.nSigmaInnerTOF(2));
            separationOuterTOF[2]->Fill(track.pt(), track.eta(), track.nSigmaInnerTOF(2));
            // separationRICH[2]->Fill((track.nSigmaKaonRich() > 3.f), track.pt(), track.eta());
          default:
            break;
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<alice3SeparationPower>(cfgc)}; }
