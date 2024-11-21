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
/// \author Roberta Ferioli (roberta.ferioli@cern.ch)
/// \since November, 2024

#include <vector>
#include <TMath.h>
#include <TPDGCode.h>
#include <TRandom.h>
#include <TVector2.h>
#include <TVector3.h>
#include <TDatabasePDG.h>
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/DataTypes.h"
#include "ReconstructionDataFormats/Track.h"
#include "ReconstructionDataFormats/PID.h"
#include "ReconstructionDataFormats/DCA.h"
#include "Common/Core/trackUtilities.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/PIDResponse.h"

using namespace std;
using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;
using std::array;

using MCTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TrackSelectionExtension, aod::TracksDCA, aod::McTrackLabels>;

struct he3FromHypertritonMap {
  HistogramRegistry registryMC{
    "registryMC",
    {},
    OutputObjHandlingPolicy::AnalysisObject,
    true,
    true};

  Configurable<int> particle_of_interest{"particle_of_interest", 0, "0=antihelium3, 1=antihypertriton"};

  // Track Parameters
  Configurable<int> min_ITS_nClusters{"min_ITS_nClusters", 7, "minimum number of found ITS clusters"};
  Configurable<int> min_TPC_nClusters{"min_TPC_nClusters", 100, "minimum number of found TPC clusters"};
  Configurable<int> min_TPC_nCrossedRows{"min_TPC_nCrossedRows", 70, "minimum number of TPC crossed pad rows"};
  Configurable<float> max_chi2_TPC{"max_chi2_TPC", 4.0f, "maximum TPC chi^2/Ncls"};
  Configurable<float> min_chi2_TPC{"min_chi2_ITS", 0.5f, "minimum TPC chi^2/Ncls"};
  Configurable<float> min_eta{"min_eta", -0.8f, "minimum_eta"};
  Configurable<float> max_eta{"max_eta", +0.8f, "maximum_eta"};
  Configurable<float> max_dcaxy{"max_dcaxy", 0.05f, "Maximum DCAxy"};
  Configurable<float> max_dcaz{"max_dcaz", 0.05f, "Maximum DCAz"};
  Configurable<float> min_nsigmaTPC{"min_nsigmaTPC", -2.0f, "Minimum nsigma TPC"};
  Configurable<float> max_nsigmaTPC{"max_nsigmaTPC", +2.0f, "Maximum nsigma TPC"};
  Configurable<float> min_pt{"min_pt", 0.0f, "minimum pt of the tracks"};
  Configurable<float> max_pt{"max_pt", 10.0f, "maximum pt of the tracks"};
  Configurable<int> nbin_pt{"nbin_pt", 50, "number of pt bins"};

  int AntihePDG = -1000020030;
  int AntiHypertritonPDG = -1010010030;

  void init(InitContext const&)
  {
    registryMC.add("he3SecPtRec", "he3SecPtRec", HistType::kTH1F, {{nbin_pt, min_pt, max_pt, "p_{T} (GeV/c)"}});
    registryMC.add("hypertritonPtgen", "hypertritonPtGen", HistType::kTH1F, {{nbin_pt, min_pt, max_pt, "p_{T} (GeV/c)"}});
  }

  void processMC(aod::McParticles const& mcParticles, const MCTracks& tracks)
  {
    for (const auto& track : tracks) {
      if (!track.has_mcParticle()) {
        continue;
      }
      auto mcparticle = track.mcParticle();
      if (mcparticle.pdgCode() != AntihePDG || mcparticle.isPhysicalPrimary()) {
        continue;
      }

      for (auto& motherparticle : mcparticle.mothers_as<aod::McParticles>()) {
        if (motherparticle.pdgCode() == AntiHypertritonPDG) {
          if (track.itsNCls() < min_ITS_nClusters ||
              track.tpcNClsFound() < min_TPC_nClusters ||
              track.tpcNClsCrossedRows() < min_TPC_nCrossedRows ||
              track.tpcNClsCrossedRows() < 0.8 * track.tpcNClsFindable() ||
              track.tpcChi2NCl() > 4.f ||
              track.tpcChi2NCl() < min_chi2_TPC ||
              track.eta() < min_eta || track.eta() > max_eta ||
              track.dcaXY() > max_dcaxy || track.dcaXY() < -max_dcaxy ||
              track.dcaZ() > max_dcaz || track.dcaZ() < -max_dcaz ||
              track.itsChi2NCl() > 36.f) {
            continue;
          }
          registryMC.fill(HIST("he3SecPtRec"), track.pt());
          registryMC.fill(HIST("hypertritonPtgen"), motherparticle.pt());
        }
      }
    }
  }
  PROCESS_SWITCH(he3FromHypertritonMap, processMC, "Process MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<he3FromHypertritonMap>(cfgc)};
}
