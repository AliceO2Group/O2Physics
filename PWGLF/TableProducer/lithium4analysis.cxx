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
// Analysis task for anti-lithium4 analysis

#include <TH1F.h>
#include <TDirectory.h>
#include <THn.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TFile.h>
#include <TH2F.h>
#include <TLorentzVector.h>
#include <TPDGCode.h>
#include <TDatabasePDG.h>

#include <cmath>
#include <array>
#include <cstdlib>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/StepTHn.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Framework/ASoAHelpers.h"
#include "DataFormatsTPC/BetheBlochAleph.h"
#include "CCDB/BasicCCDBManager.h"

#include "PWGLF/DataModel/LFLithium4Tables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;

namespace
{
constexpr double betheBlochDefault[1][6]{{-1.e32, -1.e32, -1.e32, -1.e32, -1.e32, -1.e32}};
static const std::vector<std::string> betheBlochParNames{"p0", "p1", "p2", "p3", "p4", "resolution"};
static const std::vector<std::string> particleNames{"He3"};

constexpr float he3Mass = o2::constants::physics::MassHelium3;
constexpr float protonMass = o2::constants::physics::MassProton;
constexpr int lithium4PDG = 1000030040;
constexpr int protonPDG = 2212;
constexpr int he3PDG = 1000020030;

} // namespace

struct lithium4Candidate {

  float recoPtHe3() const { return std::hypot(momHe3[0], momHe3[1]); }
  float recoPhiHe3() const { return std::atan2(momHe3[1], momHe3[0]); }
  float recoEtaHe3() const { return std::asinh(momHe3[2] / recoPtHe3()); }
  float recoPtPr() const { return std::hypot(momPr[0], momPr[1]); }
  float recoPhiPr() const { return std::atan2(momPr[1], momPr[0]); }
  float recoEtaPr() const { return std::asinh(momPr[2] / recoPtPr()); }

  std::array<float, 3> momHe3;
  std::array<float, 3> momPr;

  float nSigmaHe3 = -10;
  float he3DCAXY = -10;
  float he3DCAZ = -10;
  float protonDCAXY = -10;
  float protonDCAZ = -10;
  uint16_t tpcSignalHe3 = 0u;
  float momHe3TPC = -10.f;
  uint8_t nTPCClustersHe3 = 0u;

  bool isBkgUS = false;
  bool isBkgEM = false;
  bool isMatter = false;

  float l4PtMC = -10.f;
  float l4MassMC = -10.f;
};

struct lithium4analysis {

  Produces<aod::Lithium4Table> outputDataTable;
  Produces<aod::Lithium4TableMC> outputMCTable;

  std::vector<lithium4Candidate> l4Candidates;

  SliceCache cache;
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  // events
  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};
  // track
  Configurable<float> cfgCutPT{"cfgCutPT", 0.2, "PT cut on daughter track"};
  Configurable<float> cfgCutMaxPrPT{"cfgCutMaxPrPT", 1.8, "Max PT cut on proton"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8, "Eta cut on daughter track"};
  Configurable<float> cfgCutDCAxy{"cfgCutDCAxy", 2.0f, "DCAxy range for tracks"};
  Configurable<float> cfgCutDCAz{"cfgCutDCAz", 2.0f, "DCAz range for tracks"};
  Configurable<float> nsigmaCutTPC{"nsigmacutTPC", 3.0, "Value of the TPC Nsigma cut"};
  Configurable<float> nsigmaCutTOF{"nsigmaCutTOF", 3.0, "Value of the TOF Nsigma cut"};
  Configurable<int> cfgNoMixedEvents{"cfgNoMixedEvents", 5, "Number of mixed events per event"};
  Configurable<bool> cfgEnableBkgUS{"cfgEnableBkgUS", false, "Enable US background"};

  // bethe bloch parameters
  std::array<float, 6> mBBparamsHe;
  Configurable<LabeledArray<double>> cfgBetheBlochParams{"cfgBetheBlochParams", {betheBlochDefault[0], 1, 6, particleNames, betheBlochParNames}, "TPC Bethe-Bloch parameterisation for He3"};
  // MC
  Configurable<bool> isMC{"isMC", false, "Run MC"};
  void init(o2::framework::InitContext&)
  {

    histos.add("hCentrality", "Centrality distribution", kTH1F, {{2001, -0.5, 2000.5}});
    histos.add("hVtxZ", "Vertex distribution in Z;Z (cm)", kTH1F, {{400, -20.0, 20.0}});
    histos.add("hNcontributor", "Number of primary vertex contributor", kTH1F, {{2000, 0.0f, 2000.0f}});
    histos.add("hHe3Dcaxy", "Dcaxy distribution", kTH1F, {{200, -1.0f, 1.0f}});
    histos.add("hHe3Dcaz", "Dcaz distribution", kTH1F, {{200, -1.0f, 1.0f}});
    histos.add("hLitInvMass", "Invariant mass distribution", kTH1F, {{50, 3.74f, 3.85f}});
    histos.add("hHe3Pt", "Pt distribution", kTH1F, {{200, 0.0f, 6.0f}});
    histos.add("hProtonPt", "Pt distribution", kTH1F, {{200, 0.0f, 3.0f}});
    histos.add("h2NsigmaHe3TPC", "NsigmaHe3 TPC distribution", kTH2F, {{20, -5.0f, 5.0f}, {200, -5.0f, 5.0f}});
    histos.add("h2NsigmaHe3TOF", "NsigmaHe3 TOF distribution", kTH2F, {{20, -5.0f, 5.0f}, {200, -10.0f, 10.0f}});
    histos.add("h2NsigmaProtonTPC", "NsigmaProton TPC distribution", kTH2F, {{20, -3.0f, 3.0f}, {200, -5.0f, 5.0f}});
    histos.add("h2NsigmaProtonTOF", "NsigmaProton TOF distribution", kTH2F, {{20, -3.0f, 3.0f}, {200, -10.0f, 10.0f}});

    for (int i = 0; i < 5; i++) {
      mBBparamsHe[i] = cfgBetheBlochParams->get("He3", Form("p%i", i));
    }
    mBBparamsHe[5] = cfgBetheBlochParams->get("He3", "resolution");
  }

  template <typename T>
  bool selectionTrack(const T& candidate)
  {
    if (!(candidate.isGlobalTrack() || candidate.isPVContributor())) {
      return false;
    } else {
      if (!candidate.isGlobalTrackWoDCA() || std::abs(candidate.dcaXY()) > cfgCutDCAxy || std::abs(candidate.dcaZ()) > cfgCutDCAz || !candidate.isPVContributor()) {
        return false;
      }
    }
    return true;
  }

  template <typename T>
  bool selectionPIDProton(const T& candidate)
  {
    if (candidate.hasTOF()) {
      if (std::abs(candidate.tofNSigmaPr()) < nsigmaCutTOF && std::abs(candidate.tpcNSigmaPr()) < nsigmaCutTPC) {
        histos.fill(HIST("h2NsigmaProtonTPC"), candidate.tpcInnerParam(), candidate.tpcNSigmaPr());
        histos.fill(HIST("h2NsigmaProtonTOF"), candidate.p(), candidate.tofNSigmaPr());
        return true;
      }
    } else if (std::abs(candidate.tpcNSigmaPr()) < nsigmaCutTPC) {
      histos.fill(HIST("h2NsigmaProtonTPC"), candidate.tpcInnerParam(), candidate.tpcNSigmaPr());
      return true;
    }
    return false;
  }

  template <typename T>
  float computeNSigmaHe3(const T& candidate)
  {
    float expTPCSignal = o2::tpc::BetheBlochAleph(static_cast<float>(candidate.tpcInnerParam() * 2 / constants::physics::MassHelium3), mBBparamsHe[0], mBBparamsHe[1], mBBparamsHe[2], mBBparamsHe[3], mBBparamsHe[4]);
    double resoTPC{expTPCSignal * mBBparamsHe[5]};
    return static_cast<float>((candidate.tpcSignal() - expTPCSignal) / resoTPC);
  }

  template <typename T>
  bool selectionPIDHe3(const T& candidate)
  {
    auto nSigmaHe3 = computeNSigmaHe3(candidate);
    if (std::abs(nSigmaHe3) < nsigmaCutTPC) {
      histos.fill(HIST("h2NsigmaHe3TPC"), candidate.tpcInnerParam() * candidate.sign(), nSigmaHe3);
      return true;
    }
    return false;
  }

  template <typename T1, typename T2>
  void FillCandidateInfo(const T1& candidateHe3, const T2& candidatePr, bool mix)
  {
    lithium4Candidate l4Cand;

    l4Cand.momHe3 = array{2 * candidateHe3.px(), 2 * candidateHe3.py(), 2 * candidateHe3.pz()};
    l4Cand.momPr = array{candidatePr.px(), candidatePr.py(), candidatePr.pz()};

    float invMass = RecoDecay::m(array{l4Cand.momHe3, l4Cand.momPr}, array{he3Mass, protonMass});

    if (invMass < 3.74 || invMass > 3.85 || candidatePr.pt() > cfgCutMaxPrPT) {
      return;
    }

    histos.fill(HIST("hHe3Pt"), l4Cand.recoPtHe3());
    histos.fill(HIST("hProtonPt"), l4Cand.recoPtPr());

    int he3Sign = candidateHe3.sign();
    int protonSign = candidatePr.sign();

    l4Cand.isBkgUS = he3Sign * protonSign < 0;
    l4Cand.isBkgEM = mix;
    l4Cand.isMatter = he3Sign > 0;
    l4Cand.he3DCAXY = candidateHe3.dcaXY();
    l4Cand.he3DCAZ = candidateHe3.dcaZ();
    l4Cand.protonDCAXY = candidatePr.dcaXY();
    l4Cand.protonDCAZ = candidatePr.dcaZ();
    l4Cand.tpcSignalHe3 = candidateHe3.tpcSignal();
    l4Cand.momHe3TPC = candidateHe3.tpcInnerParam();
    l4Cand.nTPCClustersHe3 = candidateHe3.tpcNClsFound();
    l4Cand.nSigmaHe3 = computeNSigmaHe3(candidateHe3);
    l4Candidates.push_back(l4Cand);
  }

  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter acceptanceFilter = (nabs(aod::track::eta) < cfgCutEta && nabs(aod::track::pt) > cfgCutPT);
  Filter DCAcutFilter = (nabs(aod::track::dcaXY) < cfgCutDCAxy) && (nabs(aod::track::dcaZ) < cfgCutDCAz);

  using EventCandidates = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>;
  using TrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullPr, aod::pidTOFFullPr>>;
  using TrackCandidatesMC = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullPr, aod::pidTOFFullPr, aod::McTrackLabels>>;

  Preslice<TrackCandidates> perCol = aod::track::collisionId;
  Preslice<TrackCandidatesMC> perColMC = aod::track::collisionId;

  // binning for EM background
  ConfigurableAxis axisVertex{"axisVertex", {20, -10, 10}, "vertex axis for bin"};
  using BinningType = ColumnBinningPolicy<aod::collision::PosZ>;
  BinningType binningOnPositions{{axisVertex}, true};
  SameKindPair<EventCandidates, TrackCandidates, BinningType> pair{binningOnPositions, cfgNoMixedEvents, -1, &cache};

  void processSameEvent(soa::Join<aod::Collisions, aod::EvSels> const& collisions, TrackCandidates const& tracks, aod::BCs const&)
  {
    l4Candidates.clear();

    for (auto& collision : collisions) {
      if (!collision.sel8() || std::abs(collision.posZ()) > cfgCutVertex) {
        continue;
      }
      histos.fill(HIST("hNcontributor"), collision.numContrib());
      histos.fill(HIST("hVtxZ"), collision.posZ());

      const uint64_t collIdx = collision.globalIndex();
      auto TrackTable_thisCollision = tracks.sliceBy(perCol, collIdx);
      TrackTable_thisCollision.bindExternalIndices(&tracks);

      for (auto track1 : TrackTable_thisCollision) {
        if (!selectionTrack(track1) || !selectionPIDHe3(track1)) {
          continue;
        }
        histos.fill(HIST("hHe3Dcaxy"), track1.dcaXY());
        histos.fill(HIST("hHe3Dcaz"), track1.dcaZ());

        for (auto track2 : TrackTable_thisCollision) {

          if (!cfgEnableBkgUS) {
            if (track1.sign() * track2.sign() < 0) {
              continue;
            }
          }

          if (!selectionTrack(track2) || !selectionPIDProton(track2)) {
            continue;
          }

          FillCandidateInfo(track1, track2, false);
        }
      }
    }

    for (auto& l4Cand : l4Candidates) {
      outputDataTable(l4Cand.isMatter, l4Cand.recoPtHe3(), l4Cand.recoEtaHe3(), l4Cand.recoPhiHe3(),
                      l4Cand.recoPtPr(), l4Cand.recoEtaPr(), l4Cand.recoPhiPr(),
                      l4Cand.he3DCAXY, l4Cand.he3DCAZ, l4Cand.protonDCAXY, l4Cand.protonDCAZ,
                      l4Cand.tpcSignalHe3, l4Cand.momHe3TPC, l4Cand.nTPCClustersHe3, l4Cand.nSigmaHe3,
                      l4Cand.isBkgUS, l4Cand.isBkgEM);
    }
  }
  PROCESS_SWITCH(lithium4analysis, processSameEvent, "Process Same event", false);

  void processMixedEvent(soa::Join<aod::Collisions, aod::EvSels> const& collisions, TrackCandidates const& tracks)
  {
    l4Candidates.clear();

    for (auto& [c1, tracks1, c2, tracks2] : pair) {
      if (!c1.sel8()) {
        continue;
      }
      if (!c2.sel8()) {
        continue;
      }

      for (auto& [t1, t2] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {
        if (!selectionTrack(t1)) {
          continue;
        }
        if (!selectionTrack(t2)) {
          continue;
        }
        if ((selectionPIDHe3(t1) && selectionPIDProton(t2)) || (selectionPIDHe3(t2) && selectionPIDProton(t1))) {
          FillCandidateInfo(t1, t2, true);
        }
      }
    }

    for (auto& l4Cand : l4Candidates) {
      outputDataTable(l4Cand.isMatter, l4Cand.recoPtHe3(), l4Cand.recoEtaHe3(), l4Cand.recoPhiHe3(),
                      l4Cand.recoPtPr(), l4Cand.recoEtaPr(), l4Cand.recoPhiPr(),
                      l4Cand.he3DCAXY, l4Cand.he3DCAZ, l4Cand.protonDCAXY, l4Cand.protonDCAZ,
                      l4Cand.tpcSignalHe3, l4Cand.momHe3TPC, l4Cand.nTPCClustersHe3, l4Cand.nSigmaHe3,
                      l4Cand.isBkgUS, l4Cand.isBkgEM);
    }
  }
  PROCESS_SWITCH(lithium4analysis, processMixedEvent, "Process Mixed event", false);

  void processMC(soa::Join<aod::Collisions, aod::EvSels> const& collisions, aod::BCs const&, TrackCandidatesMC const& tracks, aod::McParticles const& mcParticles)
  {
    std::vector<unsigned int> filledMothers;
    l4Candidates.clear();

    for (auto& collision : collisions) {

      if (!collision.sel8() || std::abs(collision.posZ()) > cfgCutVertex) {
        continue;
      }

      histos.fill(HIST("hCentrality"), 1);
      histos.fill(HIST("hNcontributor"), collision.numContrib());
      histos.fill(HIST("hVtxZ"), collision.posZ());

      const uint64_t collIdx = collision.globalIndex();
      auto TrackTable_thisCollision = tracks.sliceBy(perColMC, collIdx);
      TrackTable_thisCollision.bindExternalIndices(&tracks);

      for (auto track1 : TrackTable_thisCollision) {

        if (!track1.has_mcParticle()) {
          continue;
        }
        if (!selectionTrack(track1) || !selectionPIDHe3(track1)) {
          continue;
        }

        for (auto track2 : TrackTable_thisCollision) {
          if (!track2.has_mcParticle()) {
            continue;
          }
          if (!selectionTrack(track2) || !selectionPIDProton(track2)) {
            continue;
          }

          if (track1.sign() * track2.sign() < 0) {
            continue;
          }

          const auto mctrackHe3 = track1.mcParticle();
          const auto mctrackPr = track2.mcParticle();

          if (std::abs(mctrackHe3.pdgCode()) != he3PDG || std::abs(mctrackPr.pdgCode()) != protonPDG) {
            continue;
          }

          for (auto& mothertrack : mctrackHe3.mothers_as<aod::McParticles>()) {
            for (auto& mothertrackPr : mctrackPr.mothers_as<aod::McParticles>()) {

              if (mothertrack != mothertrackPr || mothertrack.pdgCode() != lithium4PDG) {
                continue;
              }

              if (std::abs(mothertrack.y()) > 1) {
                continue;
              }

              FillCandidateInfo(track1, track2, false);
              auto& l4Candidate = l4Candidates.back();
              l4Candidate.l4PtMC = mothertrack.pt();
              l4Candidate.l4MassMC = std::sqrt(mothertrack.e() * mothertrack.e() - mothertrack.p() * mothertrack.p());
              filledMothers.push_back(mothertrack.globalIndex());
            }
          }
        }
      }
    }

    for (auto& mcParticle : mcParticles) {

      if (std::abs(mcParticle.pdgCode()) != lithium4PDG) {
        continue;
      }

      if (std::abs(mcParticle.y()) > 1) {
        continue;
      }

      if (std::find(filledMothers.begin(), filledMothers.end(), mcParticle.globalIndex()) != filledMothers.end()) {
        continue;
      }

      auto kDaughters = mcParticle.daughters_as<aod::McParticles>();
      auto daughtHe3 = false;
      auto daughtPr = false;
      for (auto kCurrentDaughter : kDaughters) {
        if (std::abs(kCurrentDaughter.pdgCode()) == he3PDG) {
          daughtHe3 = true;
        } else if (std::abs(kCurrentDaughter.pdgCode()) == protonPDG) {
          daughtPr = true;
        }
      }
      if (daughtHe3 && daughtPr) {
        lithium4Candidate l4Candidate;
        int sign = mcParticle.pdgCode() > 0 ? 1 : -1;
        l4Candidate.l4PtMC = mcParticle.pt() * sign;
        l4Candidate.l4MassMC = std::sqrt(mcParticle.e() * mcParticle.e() - mcParticle.p() * mcParticle.p());
        l4Candidates.push_back(l4Candidate);
      }
    }

    for (auto& l4Cand : l4Candidates) {
      outputMCTable(l4Cand.isMatter, l4Cand.recoPtHe3(), l4Cand.recoEtaHe3(), l4Cand.recoPhiHe3(),
                    l4Cand.recoPtPr(), l4Cand.recoEtaPr(), l4Cand.recoPhiPr(),
                    l4Cand.he3DCAXY, l4Cand.he3DCAZ, l4Cand.protonDCAXY, l4Cand.protonDCAZ,
                    l4Cand.tpcSignalHe3, l4Cand.momHe3TPC, l4Cand.nTPCClustersHe3, l4Cand.nSigmaHe3,
                    l4Cand.isBkgUS, l4Cand.isBkgEM,
                    l4Cand.l4PtMC, l4Cand.l4MassMC);
    }
  }
  PROCESS_SWITCH(lithium4analysis, processMC, "Process MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<lithium4analysis>(cfgc, TaskName{"lithium4analysis"})};
}
