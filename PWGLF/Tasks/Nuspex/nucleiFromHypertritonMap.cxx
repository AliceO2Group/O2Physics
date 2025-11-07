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

#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/ASoA.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/DataTypes.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/DCA.h"
#include "ReconstructionDataFormats/PID.h"
#include "ReconstructionDataFormats/Track.h"

#include <TDatabasePDG.h>
#include <TMath.h>
#include <TPDGCode.h>
#include <TRandom.h>
#include <TVector2.h>
#include <TVector3.h>

#include <vector>

using namespace std;
using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;
using std::array;

using MCTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TrackSelectionExtension, aod::TracksDCA, aod::McTrackLabels>;

struct nucleiFromHypertritonMap {
  HistogramRegistry registryMC{
    "registryMC",
    {},
    OutputObjHandlingPolicy::AnalysisObject,
    true,
    true};

  // Track Parameters
  Configurable<int> min_ITS_nClusters{"min_ITS_nClusters", 7, "minimum number of found ITS clusters"};
  Configurable<int> min_ITS_InnerBarrel_nClusters{"min_ITS_InnerBarrel_nClusters", 1, "minimum number of found ITS Inner Barrel clusters"};
  Configurable<int> max_ITS_InnerBarrel_nClusters{"max_ITS_InnerBarrel_nClusters", 3, "maximum number of found ITS Inner Barrel clusters"};
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
  Configurable<int> nbin_dca = {"nbin_dca", 50, "number of DCA bins"};
  Configurable<bool> saveHelium{"saveHelium", false, "Save helium candidates"};

  int AntideuteronPDG = -1000010020;
  int AntihePDG = -1000020030;
  int AntiHypertritonPDG = -1010010030;
  int AntiHyperHelium4PDG = -1010020040;

  void init(InitContext const&)
  {
    registryMC.add("hypertritonPtGen", "hypertritonPtGen", HistType::kTH1F, {{nbin_pt, min_pt, max_pt, "p_{T} (GeV/c)"}});
    if (saveHelium) {
      registryMC.add("he3SecPtRec_from_hypertriton", "he3SecPtRec_from_hypertriton", HistType::kTH1F, {{nbin_pt, min_pt, max_pt, "p_{T} (GeV/c)"}});
      registryMC.add("he3SecPtGen_from_hypertriton", "he3SecPtGen_from_hypertriton", HistType::kTH1F, {{nbin_pt, min_pt, max_pt, "p_{T} (GeV/c)"}});
      registryMC.add("hyperHe4PtGen", "hyperHe4PtGen", HistType::kTH1F, {{nbin_pt, min_pt, max_pt, "p_{T} (GeV/c)"}});
      registryMC.add("he3SecPtRec_from_hyperHe4", "he3SecPtRec_from_hyperHe4", HistType::kTH1F, {{nbin_pt, min_pt, max_pt, "p_{T} (GeV/c)"}});
      registryMC.add("he3SecPtGen_from_hyperHe4", "he3SecPtGen_from_hyperHe4", HistType::kTH1F, {{nbin_pt, min_pt, max_pt, "p_{T} (GeV/c)"}});
      registryMC.add("he3PtRec", "he3PtRec", HistType::kTH1F, {{nbin_pt, min_pt, max_pt, "p_{T} (GeV/c)"}});
      registryMC.add("he3PtGen", "he3PtGen", HistType::kTH1F, {{nbin_pt, min_pt, max_pt, "p_{T} (GeV/c)"}});
    } else {
      registryMC.add("deutSecPtRec_from_hypertriton", "deutSecPtRec_from_hypertriton", HistType::kTH1F, {{nbin_pt, min_pt, max_pt, "p_{T} (GeV/c)"}});
      registryMC.add("deutSecPtGen_from_hypertriton", "deutSecPtGen_from_hypertriton", HistType::kTH1F, {{nbin_pt, min_pt, max_pt, "p_{T} (GeV/c)"}});
      registryMC.add("deutPtRec", "deutPtRec", HistType::kTH1F, {{nbin_pt, min_pt, max_pt, "p_{T} (GeV/c)"}});
      registryMC.add("deutPtGen", "deutPtGen", HistType::kTH1F, {{nbin_pt, min_pt, max_pt, "p_{T} (GeV/c)"}});
    }
  }

  void processMC(const aod::McParticles& mcParticles, const MCTracks& tracks)
  {
    int selectedPDG = 0;
    if (saveHelium) {
      selectedPDG = AntihePDG;
    } else {
      selectedPDG = AntideuteronPDG;
    }

    for (const auto& mcparticle : mcParticles) {
      if (((mcparticle.pdgCode() == AntiHypertritonPDG || mcparticle.pdgCode() == AntiHyperHelium4PDG) && mcparticle.has_daughters()) || mcparticle.pdgCode() == selectedPDG) {
        if (mcparticle.pdgCode() == AntiHypertritonPDG) {
          for (auto& daughter : mcparticle.daughters_as<aod::McParticles>()) {
            if (daughter.pdgCode() == selectedPDG) {
              registryMC.fill(HIST("hypertritonPtGen"), mcparticle.pt());
              if (saveHelium) {
                registryMC.fill(HIST("he3SecPtGen_from_hypertriton"), daughter.pt());
              } else {
                registryMC.fill(HIST("deutSecPtGen_from_hypertriton"), daughter.pt());
              }
            }
          }
        }
        if (mcparticle.pdgCode() == AntiHyperHelium4PDG) {
          for (auto& daughter : mcparticle.daughters_as<aod::McParticles>()) {
            if (daughter.pdgCode() == selectedPDG) {
              registryMC.fill(HIST("hyperHe4PtGen"), mcparticle.pt());
              if (saveHelium) {
                registryMC.fill(HIST("he3SecPtGen_from_hyperHe4"), daughter.pt());
              }
            }
          }
        }
        if (mcparticle.pdgCode() == AntihePDG && mcparticle.isPhysicalPrimary()) {
          registryMC.fill(HIST("he3PtGen"), mcparticle.pt());
        }
        if (mcparticle.pdgCode() == AntideuteronPDG && mcparticle.isPhysicalPrimary()) {
          registryMC.fill(HIST("deutPtGen"), mcparticle.pt());
        }
      }
    }

    for (const auto& track : tracks) {
      if (!track.has_mcParticle()) {
        continue;
      }
      auto mcparticle = track.mcParticle();
      if (mcparticle.pdgCode() != selectedPDG) {
        continue;
      }

      if (track.itsNCls() < min_ITS_nClusters ||
          track.itsNClsInnerBarrel() < min_ITS_InnerBarrel_nClusters ||
          track.itsNClsInnerBarrel() > max_ITS_InnerBarrel_nClusters ||
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
      if (mcparticle.pdgCode() == AntideuteronPDG && mcparticle.isPhysicalPrimary()) {
        registryMC.fill(HIST("deutPtRec"), track.pt());
      }
      if (mcparticle.pdgCode() == AntihePDG && mcparticle.isPhysicalPrimary()) {
        registryMC.fill(HIST("he3PtRec"), 2 * track.pt());
      }

      for (auto& motherparticle : mcparticle.mothers_as<aod::McParticles>()) {
        if (motherparticle.pdgCode() == AntiHypertritonPDG || motherparticle.pdgCode() == AntiHyperHelium4PDG) {
          if (motherparticle.pdgCode() == AntiHypertritonPDG) {
            if (mcparticle.pdgCode() == AntihePDG) {
              registryMC.fill(HIST("he3SecPtRec_from_hypertriton"), 2 * track.pt());
            } else {
              registryMC.fill(HIST("deutSecPtRec_from_hypertriton"), track.pt());
            }
          }
          if (motherparticle.pdgCode() == AntiHyperHelium4PDG) {
            registryMC.fill(HIST("he3SecPtRec_from_hyperHe4"), 2 * track.pt());
          }
        }
      }
    }
  }
  PROCESS_SWITCH(nucleiFromHypertritonMap, processMC, "Process MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<nucleiFromHypertritonMap>(cfgc)};
}
