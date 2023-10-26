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

#include "Common/Core/PID/PIDTOF.h"
#include "Common/TableProducer/PID/pidTOFBase.h"

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

constexpr float he3Mass = o2::constants::physics::MassHelium3;
constexpr float protonMass = o2::constants::physics::MassProton;
constexpr int lithium4PDG = 1000030040;
constexpr int protonPDG = 2212;
constexpr int he3PDG = 1000020030;

enum Selections {
  kNoCuts = 0,
  kGlobalTrack,
  kTrackCuts,
  kPID,
  kAll
};

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
  float nSigmaPr = -10;
  float massTOFHe3 = -10;
  float massTOFPr = -10;

  float he3DCAXY = -10;
  float he3DCAZ = -10;
  float protonDCAXY = -10;
  float protonDCAZ = -10;
  uint16_t tpcSignalHe3 = 0u;
  float momHe3TPC = -10.f;
  float momPrTPC = -10.f;
  float invMass = -10.f;
  uint8_t nTPCClustersHe3 = 0u;

  uint8_t sharedClustersHe3 = 0u;
  uint8_t sharedClustersPr = 0u;

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
  Configurable<LabeledArray<double>> cfgBetheBlochParams{"cfgBetheBlochParams", {betheBlochDefault[0], 1, 6, {"He3"}, betheBlochParNames}, "TPC Bethe-Bloch parameterisation for He3"};
  // MC
  Configurable<bool> isMC{"isMC", false, "Run MC"};
  void init(o2::framework::InitContext&)
  {

    histos.add("hCentrality", "Centrality distribution", kTH1F, {{2001, -0.5, 2000.5}});
    histos.add("hVtxZ", "Vertex distribution in Z;Z (cm)", kTH1F, {{400, -20.0, 20.0}});
    histos.add("hNcontributor", "Number of primary vertex contributor", kTH1F, {{2000, 0.0f, 2000.0f}});
    histos.add("hHe3Dcaxy", ";DCA_{xy} (cm)", kTH1F, {{200, -1.0f, 1.0f}});
    histos.add("hHe3Dcaz", ";DCA_{z} (cm)", kTH1F, {{200, -1.0f, 1.0f}});
    histos.add("hLitInvMass", "; M(^{3}He + p) (GeV/#it{c}^{2})", kTH1F, {{50, 3.74f, 3.85f}});
    histos.add("hHe3Pt", "#it{p}_{T} distribution; #it{p}_{T} (GeV/#it{c})", kTH1F, {{200, 0.0f, 6.0f}});
    histos.add("hProtonPt", "Pt distribution; #it{p}_{T} (GeV/#it{c})", kTH1F, {{200, 0.0f, 3.0f}});
    histos.add("h2NsigmaHe3TPC", "NsigmaHe3 TPC distribution; Signed #it{p}/#it{z} (GeV/#it{c}); n#sigma_{TPC}({}^{3}He)", kTH2F, {{20, -5.0f, 5.0f}, {200, -5.0f, 5.0f}});
    histos.add("h2NsigmaProtonTPC", "NsigmaProton TPC distribution; Signed #it{p}/#it{z} (GeV/#it{c}); n#sigma_{TPC}(p)", kTH2F, {{20, -5.0f, 5.0f}, {200, -5.0f, 5.0f}});
    histos.add("h2NsigmaProtonTOF", "NsigmaProton TOF distribution; #it{p}_{T} (GeV/#it{c}); n#sigma_{TOF}(p)", kTH2F, {{20, -5.0f, 5.0f}, {200, -5.0f, 5.0f}});
    histos.add("hTrackSel", "Accepted tracks", kTH1F, {{Selections::kAll, -0.5, Selections::kAll - 0.5}});

    for (int i = 0; i < 5; i++) {
      mBBparamsHe[i] = cfgBetheBlochParams->get("He3", Form("p%i", i));
    }
    mBBparamsHe[5] = cfgBetheBlochParams->get("He3", "resolution");

    std::vector<std::string> labels = {"All", "Global track", "Track selection", "PID {}^{3}He"};
    for (int i = 0; i < Selections::kAll; i++) {
      histos.get<TH1>(HIST("hTrackSel"))->GetXaxis()->SetBinLabel(i + 1, labels[i].c_str());
    }
  }

  template <typename T>
  bool selectionTrack(const T& candidate)
  {

    if (candidate.itsNCls() < 5 ||
        candidate.tpcNClsFound() < 70 ||
        candidate.tpcNClsCrossedRows() < 70 ||
        candidate.tpcNClsCrossedRows() < 0.8 * candidate.tpcNClsFindable() ||
        candidate.tpcChi2NCl() > 4.f ||
        candidate.itsChi2NCl() > 36.f) {
      return false;
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
      return true;
    }
    return false;
  }

  template <typename T1, typename T2>
  bool FillCandidateInfo(const T1& candidateHe3, const T2& candidatePr, bool mix, bool isMC = false)
  {
    lithium4Candidate l4Cand;

    l4Cand.momHe3 = array{2 * candidateHe3.px(), 2 * candidateHe3.py(), 2 * candidateHe3.pz()};
    l4Cand.momPr = array{candidatePr.px(), candidatePr.py(), candidatePr.pz()};

    float invMass = RecoDecay::m(array{l4Cand.momHe3, l4Cand.momPr}, array{he3Mass, protonMass});

    if (invMass < 3.74 || invMass > 3.85 || candidatePr.pt() > cfgCutMaxPrPT) {
      return false;
    }

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
    l4Cand.momPrTPC = candidatePr.tpcInnerParam();
    l4Cand.invMass = invMass;

    l4Cand.nTPCClustersHe3 = candidateHe3.tpcNClsFound();
    l4Cand.nSigmaHe3 = computeNSigmaHe3(candidateHe3);
    l4Cand.nSigmaPr = candidatePr.tpcNSigmaPr();

    l4Cand.sharedClustersHe3 = candidateHe3.tpcNClsShared();
    l4Cand.sharedClustersPr = candidatePr.tpcNClsShared();

    l4Candidates.push_back(l4Cand);
    return true;
  }

  template <typename T>
  void fillHistograms(const T& l4cand)
  {
    int candSign = l4cand.isMatter ? 1 : -1;
    histos.fill(HIST("hHe3Pt"), l4cand.recoPtHe3());
    histos.fill(HIST("hProtonPt"), l4cand.recoPtPr());
    histos.fill(HIST("hLitInvMass"), l4cand.invMass);
    histos.fill(HIST("hHe3Dcaxy"), l4cand.he3DCAXY);
    histos.fill(HIST("hHe3Dcaz"), l4cand.he3DCAZ);
    histos.fill(HIST("h2NsigmaHe3TPC"), candSign * l4cand.momHe3TPC, l4cand.nSigmaHe3);
    histos.fill(HIST("h2NsigmaProtonTPC"), candSign * l4cand.momPrTPC, l4cand.nSigmaPr);
    histos.fill(HIST("h2NsigmaProtonTOF"), l4cand.recoPtPr(), l4cand.nSigmaPr);
  }

  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter acceptanceFilter = (nabs(aod::track::eta) < cfgCutEta && nabs(aod::track::pt) > cfgCutPT);
  Filter DCAcutFilter = (nabs(aod::track::dcaXY) < cfgCutDCAxy) && (nabs(aod::track::dcaZ) < cfgCutDCAz);

  using EventCandidates = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>;
  using TrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullPr, aod::pidTOFFullPr, aod::TOFSignal, aod::TOFEvTime>>;
  using TrackCandidatesMC = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullPr, aod::pidTOFFullPr, aod::TOFSignal, aod::TOFEvTime, aod::McTrackLabels>>;
  o2::pid::tof::Beta<TrackCandidates::iterator> responseBeta;
  o2::pid::tof::Beta<TrackCandidatesMC::iterator> responseBetaMC;

  Preslice<TrackCandidates> perCol = aod::track::collisionId;
  Preslice<TrackCandidatesMC> perColMC = aod::track::collisionId;

  // binning for EM background
  ConfigurableAxis axisVertex{"axisVertex", {30, -10, 10}, "vertex axis for bin"};
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

        histos.fill(HIST("hTrackSel"), Selections::kNoCuts);

        if (!track1.isGlobalTrackWoDCA()) {
          continue;
        }
        histos.fill(HIST("hTrackSel"), Selections::kGlobalTrack);

        if (!selectionTrack(track1)) {
          continue;
        }
        histos.fill(HIST("hTrackSel"), Selections::kTrackCuts);

        if (!selectionPIDHe3(track1)) {
          continue;
        }
        histos.fill(HIST("hTrackSel"), Selections::kPID);

        for (auto track2 : TrackTable_thisCollision) {

          if (track1 == track2) {
            continue;
          }

          if (!cfgEnableBkgUS) {
            if (track1.sign() * track2.sign() < 0) {
              continue;
            }
          }

          if (!track2.isGlobalTrackWoDCA()) {
            continue;
          }

          if (!selectionTrack(track2)) {
            continue;
          }

          if (!selectionPIDProton(track2)) {
            continue;
          }

          if (!FillCandidateInfo(track1, track2, false)) {
            continue;
          }
          // fill TOF info outside to avoide responseBeta crash
          auto& cand = l4Candidates.back();
          if (track1.hasTOF()) {
            float beta = responseBeta.GetBeta(track1);
            beta = std::min(1.f - 1.e-6f, std::max(1.e-4f, beta)); /// sometimes beta > 1 or < 0, to be checked
            cand.massTOFHe3 = track1.tpcInnerParam() * 2 * std::sqrt(1.f / (beta * beta) - 1.f);
          }
          if (track2.hasTOF()) {
            float beta = responseBeta.GetBeta(track2);
            beta = std::min(1.f - 1.e-6f, std::max(1.e-4f, beta)); /// sometimes beta > 1 or < 0, to be checked
            cand.massTOFPr = track2.tpcInnerParam() * std::sqrt(1.f / (beta * beta) - 1.f);
          }
          fillHistograms(cand);
        }
      }
    }

    for (auto& l4Cand : l4Candidates) {
      outputDataTable(l4Cand.isMatter, l4Cand.recoPtHe3(), l4Cand.recoEtaHe3(), l4Cand.recoPhiHe3(),
                      l4Cand.recoPtPr(), l4Cand.recoEtaPr(), l4Cand.recoPhiPr(),
                      l4Cand.he3DCAXY, l4Cand.he3DCAZ, l4Cand.protonDCAXY, l4Cand.protonDCAZ,
                      l4Cand.tpcSignalHe3, l4Cand.momHe3TPC, l4Cand.nTPCClustersHe3, l4Cand.nSigmaHe3,
                      l4Cand.nSigmaPr, l4Cand.massTOFHe3, l4Cand.massTOFPr, l4Cand.sharedClustersHe3,
                      l4Cand.sharedClustersPr, l4Cand.isBkgUS, l4Cand.isBkgEM);
    }
  }
  PROCESS_SWITCH(lithium4analysis, processSameEvent, "Process Same event", false);

  void processMixedEvent(EventCandidates& collisions, TrackCandidates const& tracks)
  {
    l4Candidates.clear();
    for (auto& [c1, tracks1, c2, tracks2] : pair) {
      if (!c1.sel8()) {
        continue;
      }
      if (!c2.sel8()) {
        continue;
      }
      histos.fill(HIST("hNcontributor"), c1.numContrib());
      histos.fill(HIST("hVtxZ"), c1.posZ());

      for (auto& [t1, t2] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {

        if (!t1.isGlobalTrackWoDCA()) {
          continue;
        }

        if (!selectionTrack(t1)) {
          continue;
        }

        if (!t2.isGlobalTrackWoDCA()) {
          continue;
        }

        if (!selectionTrack(t2)) {
          continue;
        }

        TrackCandidates::iterator he3Cand, protonCand;
        bool passPID = false;
        if (selectionPIDHe3(t1) && selectionPIDProton(t2)) {
          he3Cand = t1, protonCand = t2;
          passPID = true;
        }
        if (selectionPIDHe3(t2) && selectionPIDProton(t1)) {
          he3Cand = t2, protonCand = t1;
          passPID = true;
        }
        if (!passPID) {
          continue;
        }

        if (!FillCandidateInfo(he3Cand, protonCand, true)) {
          continue;
        }
        // fill TOF info outside to avoide responseBeta crash
        auto& cand = l4Candidates.back();
        if (he3Cand.hasTOF()) {
          float beta = responseBeta.GetBeta(he3Cand);
          beta = std::min(1.f - 1.e-6f, std::max(1.e-4f, beta)); /// sometimes beta > 1 or < 0, to be checked
          cand.massTOFHe3 = he3Cand.tpcInnerParam() * 2 * std::sqrt(1.f / (beta * beta) - 1.f);
        }
        if (protonCand.hasTOF()) {
          float beta = responseBeta.GetBeta(protonCand);
          beta = std::min(1.f - 1.e-6f, std::max(1.e-4f, beta)); /// sometimes beta > 1 or < 0, to be checked
          cand.massTOFPr = protonCand.tpcInnerParam() * std::sqrt(1.f / (beta * beta) - 1.f);
        }
        fillHistograms(cand);
      }
    }

    for (auto& l4Cand : l4Candidates) {
      outputDataTable(l4Cand.isMatter, l4Cand.recoPtHe3(), l4Cand.recoEtaHe3(), l4Cand.recoPhiHe3(),
                      l4Cand.recoPtPr(), l4Cand.recoEtaPr(), l4Cand.recoPhiPr(),
                      l4Cand.he3DCAXY, l4Cand.he3DCAZ, l4Cand.protonDCAXY, l4Cand.protonDCAZ,
                      l4Cand.tpcSignalHe3, l4Cand.momHe3TPC, l4Cand.nTPCClustersHe3, l4Cand.nSigmaHe3,
                      l4Cand.nSigmaPr, l4Cand.massTOFHe3, l4Cand.massTOFPr, l4Cand.sharedClustersHe3,
                      l4Cand.sharedClustersPr, l4Cand.isBkgUS, l4Cand.isBkgEM);
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

      histos.fill(HIST("hNcontributor"), collision.numContrib());
      histos.fill(HIST("hVtxZ"), collision.posZ());

      const uint64_t collIdx = collision.globalIndex();
      auto TrackTable_thisCollision = tracks.sliceBy(perColMC, collIdx);
      TrackTable_thisCollision.bindExternalIndices(&tracks);

      for (auto track1 : TrackTable_thisCollision) {

        if (!track1.has_mcParticle()) {
          continue;
        }

        histos.fill(HIST("hTrackSel"), Selections::kNoCuts);

        if (!track1.isGlobalTrackWoDCA()) {
          continue;
        }
        histos.fill(HIST("hTrackSel"), Selections::kGlobalTrack);

        if (!selectionTrack(track1)) {
          continue;
        }
        histos.fill(HIST("hTrackSel"), Selections::kTrackCuts);

        if (!selectionPIDHe3(track1)) {
          continue;
        }
        histos.fill(HIST("hTrackSel"), Selections::kPID);

        for (auto track2 : TrackTable_thisCollision) {
          if (!track2.has_mcParticle()) {
            continue;
          }

          if (!track2.isGlobalTrackWoDCA()) {
            continue;
          }

          if (!selectionTrack(track2)) {
            continue;
          }

          if (!selectionPIDProton(track2)) {
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

              if (mothertrack != mothertrackPr || std::abs(mothertrack.pdgCode()) != lithium4PDG) {
                continue;
              }

              if (std::abs(mothertrack.y()) > 1) {
                continue;
              }

              if (!FillCandidateInfo(track1, track2, false, true)) {
                continue;
              }

              // fill TOF info outside to avoide responseBeta crash
              auto& cand = l4Candidates.back();
              if (track1.hasTOF()) {
                float beta = responseBetaMC.GetBeta(track1);
                beta = std::min(1.f - 1.e-6f, std::max(1.e-4f, beta)); /// sometimes beta > 1 or < 0, to be checked
                cand.massTOFHe3 = track1.tpcInnerParam() * 2 * std::sqrt(1.f / (beta * beta) - 1.f);
              }
              if (track2.hasTOF()) {
                float beta = responseBetaMC.GetBeta(track2);
                beta = std::min(1.f - 1.e-6f, std::max(1.e-4f, beta)); /// sometimes beta > 1 or < 0, to be checked
                cand.massTOFPr = track2.tpcInnerParam() * std::sqrt(1.f / (beta * beta) - 1.f);
              }
              cand.l4PtMC = mothertrack.pt();
              double eLit = mctrackHe3.e() + mctrackPr.e();
              cand.l4MassMC = std::sqrt(eLit * eLit - mothertrack.p() * mothertrack.p());
              filledMothers.push_back(mothertrack.globalIndex());
              fillHistograms(cand);
            }
          }
        }
      }
    }

    for (auto& mcParticle : mcParticles) {

      if (std::abs(mcParticle.pdgCode()) != lithium4PDG) {
        continue;
      }

      if (std::abs(mcParticle.y()) > 1 || mcParticle.isPhysicalPrimary() == false) {
        continue;
      }

      if (std::find(filledMothers.begin(), filledMothers.end(), mcParticle.globalIndex()) != filledMothers.end()) {
        continue;
      }

      auto kDaughters = mcParticle.daughters_as<aod::McParticles>();
      auto daughtHe3 = false;
      auto daughtPr = false;
      double eLit = 0;
      for (auto kCurrentDaughter : kDaughters) {
        if (std::abs(kCurrentDaughter.pdgCode()) == he3PDG) {
          daughtHe3 = true;
          eLit += kCurrentDaughter.e();
        } else if (std::abs(kCurrentDaughter.pdgCode()) == protonPDG) {
          daughtPr = true;
          eLit += kCurrentDaughter.e();
        }
      }
      if (daughtHe3 && daughtPr) {
        lithium4Candidate l4Candidate;
        int sign = mcParticle.pdgCode() > 0 ? 1 : -1;
        l4Candidate.l4PtMC = mcParticle.pt() * sign;
        l4Candidate.l4MassMC = std::sqrt(eLit * eLit - mcParticle.p() * mcParticle.p());
        l4Candidates.push_back(l4Candidate);
      }
    }

    for (auto& l4Cand : l4Candidates) {
      outputMCTable(l4Cand.isMatter, l4Cand.recoPtHe3(), l4Cand.recoEtaHe3(), l4Cand.recoPhiHe3(),
                    l4Cand.recoPtPr(), l4Cand.recoEtaPr(), l4Cand.recoPhiPr(),
                    l4Cand.he3DCAXY, l4Cand.he3DCAZ, l4Cand.protonDCAXY, l4Cand.protonDCAZ,
                    l4Cand.tpcSignalHe3, l4Cand.momHe3TPC, l4Cand.nTPCClustersHe3, l4Cand.nSigmaHe3,
                    l4Cand.nSigmaPr, l4Cand.massTOFHe3, l4Cand.massTOFPr, l4Cand.sharedClustersHe3,
                    l4Cand.sharedClustersPr, l4Cand.isBkgUS, l4Cand.isBkgEM,
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
