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

/// \file pccQa.cxx
/// \brief Task producing DCA distributions with and without particle-composition correction.
/// \author Mario Krüger <mario.kruger@cern.ch>

#include "PWGLF/DataModel/particleCompositionCorrectionTable.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>
#include <TH2.h>
#include <TMCProcess.h>
#include <TProfile.h>

#include <fmt/format.h>

#include <cstdint>
#include <cstdlib>
#include <memory>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using aod::track::TrackSelectionFlags;

int getMotherCategory(int pdg)
{
  switch (std::abs(pdg)) {
    case 310:
      return 1; // K0S
    case 130:
      return 2; // K0L
    case 3122:
      return 3; // Lambda

    case 3222: // Sigma+
    case 3112:
      return 4; // Sigma-

    case 3312:
      return 5; // Xi-

    case 3334:
      return 6; // Omega

    case 211:
      return 7; // Pion
    case 321:
      return 8; // Kaon
    case 2112:
      return 9; // Neutron
    case 13:
      return 10; // Muon

    default:
      return 11; // Other
  }
}

static std::string getMotherName(int category)
{
  static const std::vector<std::string> names = {
    "K0S",
    "K0L",
    "Lambda",
    "SigmaCharged",
    "XiCharged",
    "Omega",
    "Pi",
    "K",
    "n",
    "Muon",
    "Other"};

  return names.at(category - 1);
}

struct PccQa {
  HistogramRegistry histos;
  Service<o2::framework::O2DatabasePDG> pdg;

  std::vector<std::shared_ptr<TH2>> motherDCA;
  std::vector<std::shared_ptr<TH2>> motherDCAWeighted;
  std::vector<std::shared_ptr<TProfile>> motherWeights;

  static constexpr float MaxVtxZ = 10.f;

  Configurable<uint32_t> maxMult{"maxMult", 100, "max multiplicity"};

  void init(InitContext const&);

  template <bool IS_MC, typename C, typename T>
  void processMeas(const C& collision, const T& tracks);

  using CollisionTableData = soa::Join<aod::Collisions, aod::EvSels>;
  using TrackTableData = soa::Join<aod::FullTracks, aod::TracksDCA, aod::TrackSelection>;
  void processData(CollisionTableData::iterator const& collision, TrackTableData const& tracks);
  PROCESS_SWITCH(PccQa, processData, "process data", false);

  using CollisionTableMCTrue = aod::McCollisions;
  using CollisionTableMC = soa::SmallGroups<soa::Join<aod::McCollisionLabels, aod::Collisions, aod::EvSels>>;
  using TrackTableMC = soa::Join<aod::FullTracks, aod::TracksDCA, aod::TrackSelection, aod::McTrackLabels>;
  using ParticleTableMC = soa::Join<aod::McParticles, aod::ParticleCompositionCorrection>;
  Preslice<TrackTableMC> perCollision = aod::track::collisionId;
  void processMC(CollisionTableMCTrue::iterator const& mcCollision, TrackTableMC const& tracks, CollisionTableMC const& collisions, ParticleTableMC const&);
  PROCESS_SWITCH(PccQa, processMC, "process mc", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<PccQa>(cfgc)};
}

void PccQa::init(InitContext const&)
{

  // const AxisSpec dcaAxis = {1000, -1., 1., "#it{DCA}_{xy}", "dca"};
  const AxisSpec dcaAxis = {99, -1., 1., "#it{DCA}_{xy}", "dca"};
  std::vector<double> ptBinEdges = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0};
  const AxisSpec ptAxis{ptBinEdges, "#it{p}_{T} (GeV/#it{c})", "pt"};

  const int nBinsMult = 25;
  const AxisSpec nchAxis = {nBinsMult, 1., static_cast<double>(maxMult), "#it{N}_{ch}", "nch"};

  histos.add("multDist", "", kTH1D, {nchAxis});

  histos.add("DCAxyVsPtVsMult", "", kTH3D, {nchAxis, ptAxis, dcaAxis});

  if (doprocessMC) {
    histos.add("DCAxyVsPtVsMult_weighted", "", kTH3D, {nchAxis, ptAxis, dcaAxis});

    histos.add("prim/DCAxyVsPtVsMult", "", kTH3D, {nchAxis, ptAxis, dcaAxis});
    histos.add("prim/DCAxyVsPtVsMult_weighted", "", kTH3D, {nchAxis, ptAxis, dcaAxis});

    histos.add("sec/DCAxyVsPtVsMult", "", kTH3D, {nchAxis, ptAxis, dcaAxis});
    histos.add("sec/DCAxyVsPtVsMult_weighted", "", kTH3D, {nchAxis, ptAxis, dcaAxis});

    histos.add("sec/dec/DCAxyVsPtVsMult", "", kTH3D, {nchAxis, ptAxis, dcaAxis});
    histos.add("sec/dec/DCAxyVsPtVsMult_weighted", "", kTH3D, {nchAxis, ptAxis, dcaAxis});

    histos.add("sec/mat/DCAxyVsPtVsMult", "", kTH3D, {nchAxis, ptAxis, dcaAxis});
    histos.add("sec/mat/DCAxyVsPtVsMult_weighted", "", kTH3D, {nchAxis, ptAxis, dcaAxis});

    histos.add("prim/weights", "", kTProfile, {ptAxis});
    histos.add("sec/weights", "", kTProfile, {ptAxis});
    histos.add("sec/dec/weights", "", kTProfile, {ptAxis});
    histos.add("sec/mat/weights", "", kTProfile, {ptAxis});

    std::vector<std::string> labels = {
      "K^{0}_{S}", "K^{0}_{L}", "#Lambda",
      "#Sigma^{#pm}",
      "#Xi^{#pm}",
      "#Omega^{#pm}",
      "#pi^{#pm}", "K^{#pm}", "n",
      "#mu^{#pm}",
      "Other"};

    const int nMotherCategories = static_cast<int>(labels.size());

    histos.add(
      "sec/dec/motherPDG",
      "Mother particles of decay secondaries",
      kTH1D,
      {{nMotherCategories, 0.5, nMotherCategories + 0.5}});

    auto hMother = histos.get<TH1>(HIST("sec/dec/motherPDG"));

    for (int i = 0; i < nMotherCategories; ++i) {
      hMother->GetXaxis()->SetBinLabel(i + 1, labels[i].c_str());
    }

    for (int category = 1; category <= nMotherCategories; ++category) {
      const auto mother = getMotherName(category);

      motherDCA.push_back(histos.add<TH2>(
        fmt::format("sec/dec/{}/DCAxyVsPt", mother),
        "", kTH2D, {ptAxis, dcaAxis}));

      motherDCAWeighted.push_back(histos.add<TH2>(
        fmt::format("sec/dec/{}/DCAxyVsPt_weighted", mother),
        "", kTH2D, {ptAxis, dcaAxis}));

      motherWeights.push_back(histos.add<TProfile>(
        fmt::format("sec/dec/{}/weights", mother),
        "", kTProfile, {ptAxis}));
    }
  }
}

void PccQa::processData(CollisionTableData::iterator const& collision, TrackTableData const& tracks)
{
  processMeas<false>(collision, tracks);
}
void PccQa::processMC(CollisionTableMCTrue::iterator const&, TrackTableMC const& tracks, CollisionTableMC const& collisions, ParticleTableMC const&)
{

  for (const auto& collision : collisions) {
    auto curTracks = tracks.sliceBy(perCollision, collision.globalIndex());
    processMeas<true>(collision, curTracks);
    break;
  }
}

template <bool IS_MC, typename C, typename T>
void PccQa::processMeas(const C& collision, const T& tracks)
{

  if ((std::abs(collision.posZ()) > MaxVtxZ) || !collision.sel8()) {
    return;
  }

  auto MultCounterMeas = 0;

  for (const auto& track : tracks) {

    if (!TrackSelectionFlags::checkFlag(track.trackCutFlag(), TrackSelectionFlags::kGlobalTrackWoDCA)) {
      continue;
    }

    MultCounterMeas++;
  }

  /* LOGP(info, "MCB:      --------------------------");
  LOGP(info, "MCB:      ++ New event with Nch = {} ",MultCounterMeas);
  LOGP(info, "MCB:      --------------------------"); */

  histos.fill(HIST("multDist"), MultCounterMeas);

  for (const auto& track : tracks) {
    if (!TrackSelectionFlags::checkFlag(track.trackCutFlag(), TrackSelectionFlags::kGlobalTrackWoDCA)) {
      continue;
    }

    histos.fill(HIST("DCAxyVsPtVsMult"), MultCounterMeas, track.pt(), track.dcaXY());

    if constexpr (IS_MC) {
      if (!track.has_mcParticle()) {
        continue;
      }

      const auto& particle = track.template mcParticle_as<ParticleTableMC>();

      histos.fill(HIST("DCAxyVsPtVsMult_weighted"), MultCounterMeas, track.pt(), track.dcaXY(), particle.pccWeight());

      if (particle.isPhysicalPrimary()) {

        histos.fill(HIST("prim/weights"), particle.pt(), particle.pccWeight());
        histos.fill(HIST("prim/DCAxyVsPtVsMult"), MultCounterMeas, track.pt(), track.dcaXY());
        histos.fill(HIST("prim/DCAxyVsPtVsMult_weighted"), MultCounterMeas, track.pt(), track.dcaXY(), particle.pccWeight());

      } else {

        if (!particle.has_mothers()) {
          continue;
        }

        histos.fill(HIST("sec/weights"), particle.pt(), particle.pccWeight());

        histos.fill(HIST("sec/DCAxyVsPtVsMult"), MultCounterMeas, track.pt(), track.dcaXY());
        histos.fill(HIST("sec/DCAxyVsPtVsMult_weighted"), MultCounterMeas, track.pt(), track.dcaXY(), particle.pccWeight());

        if (particle.getProcess() == TMCProcess::kPDecay || particle.getProcess() == TMCProcess::kPPrimary) {

          histos.fill(HIST("sec/dec/DCAxyVsPtVsMult"),
                      MultCounterMeas, track.pt(), track.dcaXY());

          histos.fill(HIST("sec/dec/DCAxyVsPtVsMult_weighted"),
                      MultCounterMeas, track.pt(), track.dcaXY(),
                      particle.pccWeight());

          histos.fill(HIST("sec/dec/weights"),
                      particle.pt(), particle.pccWeight());

          /* LOGP(info,
            " -> Particle={}, PDG={}, process={}, "
            "nMoms={}, 1. mom ID={}",
            particle.globalIndex(),
            particle.pdgCode(),
            particle.getProcess(),
            particle.mothersIds().size(),
            particle.mothersIds()[0]); */

          // Identify the immediate mother
          // auto mother = particle.template mothers_first_as<aod::McParticles>();
          auto mother = particle.template mothers_first_as<ParticleTableMC>();

          const int motherCategory = getMotherCategory(mother.pdgCode());
          // const std::string motherName = getMotherName(motherCategory);

          histos.fill(
            HIST("sec/dec/motherPDG"),
            motherCategory);

          const int i = motherCategory - 1;

          /* LOGP(info,
              "  --> Mother PDG={}",
              mother.pdgCode()
              ); */

          motherDCA[i]->Fill(track.pt(), track.dcaXY());

          motherDCAWeighted[i]->Fill(
            track.pt(), track.dcaXY(), particle.pccWeight());

          motherWeights[i]->Fill(
            particle.pt(), particle.pccWeight());

        } else {

          histos.fill(HIST("sec/mat/weights"), particle.pt(), particle.pccWeight());
          histos.fill(HIST("sec/mat/DCAxyVsPtVsMult"), MultCounterMeas, track.pt(), track.dcaXY());
          histos.fill(HIST("sec/mat/DCAxyVsPtVsMult_weighted"), MultCounterMeas, track.pt(), track.dcaXY(), particle.pccWeight());
        }
      }
    }
  }
}
