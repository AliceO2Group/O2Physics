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

/// \file AntiNucleiTask.cxx
/// \brief A task to analyse Anti-nuclei
/// \author Arkaprabha Saha <arkaprabha.saha@cern.ch>

#include "Common/Core/PID/TPCPIDResponse.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "DataFormatsTPC/BetheBlochAleph.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include <TMCProcess.h>
#include <TParameter.h>

#include <cmath>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;

using CollisionWithEvSel = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs>;
using TotalTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::pidTOFDe>;
using CollisionMC = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs, aod::McCollisionLabels>;
using TracksMC = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::pidTOFDe, aod::McTrackLabels>;

namespace
{
static const int kPdgCodeHe3 = 1000020030;
static const float kMaxEta = 0.8f;
static const float kMaxGenRapidity = 0.5f;
static const float kMaxVertexZ = 10.f;
static const int kMinTpcCrossedRows = 70;
static const double kBetheBlochDefault[6]{-1.e32, -1.e32, -1.e32, -1.e32, -1.e32, -1.e32};
static const std::vector<std::string> kParamNamesBB{"p0", "p1", "p2", "p3", "p4", "resolution"};
static const std::vector<std::string> kParticleNamesBB{"He3"};
} // namespace

struct AntiNucleiTask {
  HistogramRegistry histo{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<float> mMinClsTPC{"minClsTPC", 70.0f, ""};
  Configurable<float> mMinClsITS{"minClsITS", 4.0f, ""};
  Configurable<float> mMaxChi2TPC{"maxChi2TPC", 4.0f, ""};
  Configurable<float> mMinChi2TPC{"minChi2TPC", 0.0f, ""};
  Configurable<float> mMaxChi2ITS{"maxChi2ITS", 36.0f, ""};
  Configurable<float> mMaxDCAZ{"maxDCAZ", 0.02f, ""};
  Configurable<float> mMaxDCAXY{"maxDCAXY", 0.02f, ""};
  Configurable<float> mMaxNsigmaTPC{"maxNsigmaTPC", 3.0f, ""};
  Configurable<LabeledArray<double>> mBetheBlochParams{"betheBlochParams", {kBetheBlochDefault, 1, 6, kParticleNamesBB, kParamNamesBB}, ""};

  void init(InitContext const&)
  {
    ConfigurableAxis axisEta{"eta", {16, -0.8, +0.8}, "#eta"};
    ConfigurableAxis axisPhi{"phi", {70, 0.f, 7.f}, "#phi (rad)"};
    ConfigurableAxis axisZVtx{"zVtx", {100, -20.f, 20.f}, "Vertex z (cm)"};
    ConfigurableAxis axisNSigma{"nSigma", {50, -5.f, 5.f}, "n#sigma"};
    ConfigurableAxis axisPt{"pt", {200, -10.0f, 10.0f}, "#it{p}_{T} (GeV/#it{c})"};
    ConfigurableAxis axisMomTPC{"momTPC", {5.e2, 0.f, 5.f}, "#it{p}_{TPC} (GeV/#it{c})"};
    ConfigurableAxis axisSignalTPC{"signalTPC", {4.e2, 0.f, 4.e3f}, "d#it{E}/d#it{x}_{TPC} (a.u.)"};
    ConfigurableAxis axisCent{"centAxis", {100, 0, 100.0f}, "Centrality"};

    histo.add("ZVtx_Raw", "Raw Z-vertex", kTH1F, {{axisZVtx}});
    histo.add("ZVtx", "Z-vertex", kTH1F, {{axisZVtx}});
    histo.add("Eta_Raw", "Raw Eta", kTH1F, {{axisEta}});
    histo.add("Eta", "Eta", kTH1F, {{axisEta}});
    histo.add("Phi_Raw", "Raw Phi", kTH1F, {{axisPhi}});
    histo.add("Phi", "Phi", kTH1F, {{axisPhi}});
    histo.add("Pt_Raw", "Raw Pt", kTH1F, {{axisPt}});
    histo.add("Pt", "Pt", kTH1F, {{axisPt}});
    histo.add("TPCSignal", "TPC dEdx vs Momentum", kTH2F, {{axisMomTPC}, {axisSignalTPC}});
    histo.add("NSigmaTPC_Raw", "Raw TPC nSigma", kTH3F, {{axisCent}, {axisPt}, {axisNSigma}});
    histo.add("NSigmaTOF_Raw", "Raw TOF nSigma", kTH3F, {{axisCent}, {axisPt}, {axisNSigma}});
    histo.add("NSigmaTPC", "TPC nSigma", kTH3F, {{axisCent}, {axisPt}, {axisNSigma}});
    histo.add("NSigmaTOF", "TOF nSigma", kTH3F, {{axisCent}, {axisPt}, {axisNSigma}});
    histo.add("GenPt", "Generated He3 (|y|<0.5)", kTH1F, {{axisPt}});
    histo.add("RecPt", "Reconstructed He3 (|#eta|<0.8)", kTH1F, {{axisPt}});
  }

  template <typename T>
  bool passTrackSelection(const T& track)
  {
    if (std::abs(track.eta()) > kMaxEta)
      return false;
    if (track.tpcNClsFound() < mMinClsTPC)
      return false;
    if (track.tpcNClsCrossedRows() < kMinTpcCrossedRows)
      return false;
    if (track.itsNCls() < mMinClsITS)
      return false;
    if (track.tpcChi2NCl() > mMaxChi2TPC || track.tpcChi2NCl() < mMinChi2TPC)
      return false;
    if (track.itsChi2NCl() > mMaxChi2ITS)
      return false;
    if (std::abs(track.dcaXY()) > mMaxDCAXY || std::abs(track.dcaZ()) > mMaxDCAZ)
      return false;
    return true;
  }

  static float getSignedPtMC(auto const& particle)
  {
    return (particle.pdgCode() > 0) ? particle.pt() : -particle.pt();
  }

  static bool isPrimaryHe3(auto const& particle)
  {
    return std::abs(particle.pdgCode()) == kPdgCodeHe3 && particle.isPhysicalPrimary();
  }

  void process(CollisionWithEvSel::iterator const& collision, TotalTracks const& tracks)
  {
    if (std::abs(collision.posZ()) > kMaxVertexZ)
      return;
    histo.fill(HIST("ZVtx_Raw"), collision.posZ());

    if (!collision.sel8())
      return;
    histo.fill(HIST("ZVtx"), collision.posZ());

    for (const auto& track : tracks) {
      double expectedSignal = tpc::BetheBlochAleph(
        static_cast<double>(track.tpcInnerParam()),
        mBetheBlochParams->get("p0"), mBetheBlochParams->get("p1"),
        mBetheBlochParams->get("p2"), mBetheBlochParams->get("p3"),
        mBetheBlochParams->get("p4"));

      float nSigmaTPC = static_cast<float>((track.tpcSignal() - expectedSignal) / (expectedSignal * mBetheBlochParams->get("resolution")));
      float signedPt = (track.sign() > 0) ? 2 * track.pt() : -2 * track.pt();

      histo.fill(HIST("Eta_Raw"), track.eta());
      histo.fill(HIST("Phi_Raw"), track.phi());
      histo.fill(HIST("Pt_Raw"), signedPt);
      histo.fill(HIST("NSigmaTPC_Raw"), collision.centFT0C(), signedPt, nSigmaTPC);
      histo.fill(HIST("NSigmaTOF_Raw"), collision.centFT0C(), signedPt, track.tofNSigmaDe());

      if (passTrackSelection(track)) {
        histo.fill(HIST("Eta"), track.eta());
        histo.fill(HIST("Phi"), track.phi());
        histo.fill(HIST("Pt"), signedPt);
        histo.fill(HIST("NSigmaTPC"), collision.centFT0C(), signedPt, nSigmaTPC);
        histo.fill(HIST("TPCSignal"), track.tpcInnerParam(), track.tpcSignal());

        if (std::abs(nSigmaTPC) < mMaxNsigmaTPC) {
          histo.fill(HIST("NSigmaTOF"), collision.centFT0C(), signedPt, track.tofNSigmaDe());
        }
      }
    }
  }

  void processMC(CollisionMC const& collisions, aod::McCollisions const&, TracksMC const& tracks, aod::McParticles const& particlesMC)
  {
    for (const auto& collision : collisions) {
      if (std::abs(collision.posZ()) > kMaxVertexZ || !collision.sel8())
        continue;

      for (const auto& particle : particlesMC) {
        if (particle.mcCollisionId() == collision.mcCollisionId() && isPrimaryHe3(particle)) {
          if (std::abs(particle.y()) < kMaxGenRapidity) {
            histo.fill(HIST("GenPt"), getSignedPtMC(particle));
          }
        }
      }

      for (const auto& track : tracks) {
        if (track.collisionId() != collision.index() || !passTrackSelection(track) || track.mcParticleId() == -1)
          continue;

        auto particle = particlesMC.iteratorAt(track.mcParticleId());
        if (isPrimaryHe3(particle)) {
          if (std::abs(track.rapidity(MassHelium3)) < kMaxGenRapidity) {
            histo.fill(HIST("RecPt"), getSignedPtMC(particle));
          }
        }
      }
    }
  }

  PROCESS_SWITCH(AntiNucleiTask, processMC, "processMC", true);
  PROCESS_SWITCH(AntiNucleiTask, process, "process", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<AntiNucleiTask>(cfgc)};
}
