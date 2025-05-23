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
//
/// \file hyperhelium4sigmaAnalysis.cxx
/// \brief QA and analysis task for hyper-helium4sigma (He4S)
/// \author Yuanzhe Wang <yuanzhe.wang@cern.ch>

#include <vector>
#include <memory>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponse.h"
#include "CommonConstants/PhysicsConstants.h"
#include "PWGLF/DataModel/LFKinkDecayTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using CollisionsFull = soa::Join<aod::Collisions, aod::EvSels>;
using FullTracksExtIU = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::pidTPCFullPr, aod::pidTPCFullAl, aod::pidTPCFullTr, aod::pidTPCFullPi>;
using MCLabeledTracksIU = soa::Join<FullTracksExtIU, aod::McTrackLabels>;

namespace
{
constexpr int kITSLayers = 7;
constexpr int kITSInnerBarrelLayers = 3;
// constexpr int kITSOuterBarrelLayers = 4;
std::shared_ptr<TH1> hMotherCounter;
std::shared_ptr<TH1> hDauAlphaCounter;
std::shared_ptr<TH1> hDauTritonCounter;
std::shared_ptr<TH1> hDauProtonCounter;
std::shared_ptr<TH1> hDauPionCounter;
} // namespace

//-------------------------------Check the decay channel of H4S-------------------------------
enum Channel {
  k2body = 0, // helium4, pion0
  k3body_p,   // triton, proton, pion0
  k3body_n,   // triton, neutron, pion+
  kNDecayChannel
};

template <class TMCTrackTo, typename TMCParticle>
Channel getDecayChannelH4S(TMCParticle const& particle, std::vector<int>& list)
{
  if (std::abs(particle.pdgCode()) != o2::constants::physics::Pdg::kHyperHelium4Sigma) {
    return kNDecayChannel;
  }

  // list: charged, charged or empty, neutral
  list.clear();
  list.resize(3, -1);

  bool haveAlpha = false, haveTriton = false, haveProton = false, haveNeuteron = false;
  bool haveAntiAlpha = false, haveAntiTriton = false, haveAntiProton = false, haveAntiNeuteron = false;
  bool havePionPlus = false, havePionMinus = false, havePion0 = false;
  auto daughters = particle.template daughters_as<TMCTrackTo>();
  for (const auto& mcDaughter : particle.template daughters_as<TMCTrackTo>()) {
    if (mcDaughter.pdgCode() == o2::constants::physics::Pdg::kAlpha) {
      haveAlpha = true;
      list[0] = mcDaughter.globalIndex();
    }
    if (mcDaughter.pdgCode() == -o2::constants::physics::Pdg::kAlpha) {
      haveAntiAlpha = true;
      list[0] = mcDaughter.globalIndex();
    }
    if (mcDaughter.pdgCode() == o2::constants::physics::Pdg::kTriton) {
      haveTriton = true;
      list[0] = mcDaughter.globalIndex();
    }
    if (mcDaughter.pdgCode() == -o2::constants::physics::Pdg::kTriton) {
      haveAntiTriton = true;
      list[0] = mcDaughter.globalIndex();
    }
    if (mcDaughter.pdgCode() == PDG_t::kProton) {
      haveProton = true;
      list[1] = mcDaughter.globalIndex();
    }
    if (mcDaughter.pdgCode() == -PDG_t::kProton) {
      haveAntiProton = true;
      list[1] = mcDaughter.globalIndex();
    }
    if (mcDaughter.pdgCode() == PDG_t::kNeutron) {
      haveNeuteron = true;
      list[2] = mcDaughter.globalIndex();
    }
    if (mcDaughter.pdgCode() == -PDG_t::kNeutron) {
      haveAntiNeuteron = true;
      list[2] = mcDaughter.globalIndex();
    }
    if (mcDaughter.pdgCode() == PDG_t::kPiPlus) {
      havePionPlus = true;
      list[1] = mcDaughter.globalIndex();
    }
    if (mcDaughter.pdgCode() == -PDG_t::kPiPlus) {
      havePionMinus = true;
      list[1] = mcDaughter.globalIndex();
    }
    if (mcDaughter.pdgCode() == PDG_t::kPi0) {
      havePion0 = true;
      list[2] = mcDaughter.globalIndex();
    }
  }

  if ((haveAlpha && havePion0) || (haveAntiAlpha && havePion0)) {
    return k2body;
  } else if ((haveTriton && haveProton && havePion0) || (haveAntiTriton && haveAntiProton && havePion0)) {
    return k3body_p;
  } else if ((haveTriton && haveNeuteron && havePionPlus) || (haveAntiTriton && haveAntiNeuteron && havePionMinus)) {
    return k3body_n;
  }

  return kNDecayChannel;
}

//--------------------------------------------------------------
// check if the mcparticle is daughter of hyperhelium4sigma
template <typename TMCParticle>
bool isDaughterTrack(TMCParticle const& mcparticle, int pdgcode = o2::constants::physics::Pdg::kHyperHelium4Sigma)
{
  if (!mcparticle.has_mothers()) {
    return false;
  }
  for (const auto& particleMother : mcparticle.template mothers_as<aod::McParticles>()) {
    if (std::abs(particleMother.pdgCode()) == pdgcode) {
      return true;
    }
  }
  return false;
}

//--------------------------------------------------------------
// analysis task for hyperhelium4sigma 2-body decay
struct Hyperhelium4sigmaAnalysis {
  // Histograms are defined with HistogramRegistry
  HistogramRegistry registry{"registry", {}};

  // Configurable for event selection
  Configurable<bool> doEventCut{"doEventCut", true, "Apply event selection"};
  Configurable<float> cutzvertex{"cutzvertex", 10.0f, "Accepted z-vertex range (cm)"};
  Configurable<float> cutNSigmaAl{"cutNSigmaAl", 5, "NSigmaTPCAlpha"};

  void init(InitContext const&)
  {
    // Axes
    const AxisSpec vertexZAxis{100, -15., 15., "vrtx_{Z} [cm]"};
    const AxisSpec ptAxis{50, -10, 10, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec nSigmaAxis{120, -6.f, 6.f, "n#sigma_{#alpha}"};
    const AxisSpec massAxis{100, 3.85, 4.25, "m (GeV/#it{c}^{2})"};

    registry.add<TH1>("hEventCounter", "hEventCounter", HistType::kTH1F, {{2, 0, 2}});
    registry.add<TH1>("hVertexZCollision", "hVertexZCollision", HistType::kTH1F, {vertexZAxis});

    if (doprocessData == true) {
      registry.add<TH1>("hCandidateCounter", "hCandidateCounter", HistType::kTH1F, {{3, 0, 3}});
    }

    if (doprocessMC == true) {
      registry.add<TH1>("hCandidateCounter", "hCandidateCounter", HistType::kTH1F, {{6, 0, 6}});
      registry.add<TH1>("hDiffSVx", ";;#Delta x (cm)", HistType::kTH1F, {{200, -10, 10}});
      registry.add<TH1>("hDiffSVy", ";;#Delta y (cm)", HistType::kTH1F, {{200, -10, 10}});
      registry.add<TH1>("hDiffSVz", ";;#Delta z (cm)", HistType::kTH1F, {{200, -10, 10}});
      registry.add<TH1>("hDiffDauPx", ";#Delta p_{x} (GeV/#it{c}); ", HistType::kTH1D, {{200, -10, 10}});
      registry.add<TH1>("hDiffDauPy", ";#Delta p_{y} (GeV/#it{c}); ", HistType::kTH1D, {{200, -10, 10}});
      registry.add<TH1>("hDiffDauPz", ";#Delta p_{z} (GeV/#it{c}); ", HistType::kTH1D, {{200, -10, 10}});
    }

    registry.add<TH2>("h2MassHyperhelium4sigmaPt", "h2MassHyperhelium4sigmaPt", HistType::kTH2F, {{ptAxis, massAxis}});
    registry.add<TH2>("h2NSigmaAlPt", "h2NSigmaAlPt", HistType::kTH2F, {{ptAxis, nSigmaAxis}});
  }

  void processData(CollisionsFull const& collisions, aod::KinkCands const& KinkCands, FullTracksExtIU const&)
  {
    for (const auto& collision : collisions) {
      registry.fill(HIST("hEventCounter"), 0);
      if (doEventCut && (std::abs(collision.posZ()) > cutzvertex || !collision.sel8())) {
        continue;
      }
      registry.fill(HIST("hEventCounter"), 1);
      registry.fill(HIST("hVertexZCollision"), collision.posZ());
    }

    for (const auto& kinkCand : KinkCands) {
      registry.fill(HIST("hCandidateCounter"), 0);
      auto collision = kinkCand.collision_as<CollisionsFull>();
      if (doEventCut && (std::abs(collision.posZ()) > cutzvertex || !collision.sel8())) {
        continue;
      }
      registry.fill(HIST("hCandidateCounter"), 1);

      auto dauTrack = kinkCand.trackDaug_as<FullTracksExtIU>();
      if (std::abs(dauTrack.tpcNSigmaAl()) > cutNSigmaAl) {
        continue;
      }
      float invMass = RecoDecay::m(std::array{std::array{kinkCand.pxDaug(), kinkCand.pyDaug(), kinkCand.pzDaug()}, std::array{kinkCand.pxDaugNeut(), kinkCand.pyDaugNeut(), kinkCand.pzDaugNeut()}}, std::array{o2::constants::physics::MassAlpha, o2::constants::physics::MassPi0});
      registry.fill(HIST("hCandidateCounter"), 2);
      registry.fill(HIST("h2MassHyperhelium4sigmaPt"), kinkCand.mothSign() * kinkCand.ptMoth(), invMass);
      registry.fill(HIST("h2NSigmaAlPt"), kinkCand.mothSign() * kinkCand.ptDaug(), dauTrack.tpcNSigmaAl());
    }
  }
  PROCESS_SWITCH(Hyperhelium4sigmaAnalysis, processData, "process data", true);

  void processMC(CollisionsFull const& collisions, aod::KinkCands const& KinkCands, MCLabeledTracksIU const&, aod::McParticles const&)
  {
    std::vector<int> dauIDList;

    for (const auto& collision : collisions) {
      registry.fill(HIST("hEventCounter"), 0);
      if (doEventCut && (!collision.selection_bit(aod::evsel::kIsTriggerTVX) || !collision.selection_bit(aod::evsel::kNoTimeFrameBorder) || std::abs(collision.posZ()) > cutzvertex)) {
        continue;
      }
      registry.fill(HIST("hEventCounter"), 1);
      registry.fill(HIST("hVertexZCollision"), collision.posZ());
    }

    for (const auto& kinkCand : KinkCands) {
      registry.fill(HIST("hCandidateCounter"), 0);
      auto collision = kinkCand.collision_as<CollisionsFull>();
      if (doEventCut && (!collision.selection_bit(aod::evsel::kIsTriggerTVX) || !collision.selection_bit(aod::evsel::kNoTimeFrameBorder) || std::abs(collision.posZ()) > cutzvertex)) {
        continue;
      }
      registry.fill(HIST("hCandidateCounter"), 1);

      auto motherTrack = kinkCand.trackMoth_as<MCLabeledTracksIU>();
      auto dauTrack = kinkCand.trackDaug_as<MCLabeledTracksIU>();
      if (!motherTrack.has_mcParticle() || !dauTrack.has_mcParticle()) {
        continue;
      }
      registry.fill(HIST("hCandidateCounter"), 2);

      auto mcMotherTrack = motherTrack.mcParticle_as<aod::McParticles>();
      auto mcDauTrack = dauTrack.mcParticle_as<aod::McParticles>();
      auto dChannel = getDecayChannelH4S<aod::McParticles>(mcMotherTrack, dauIDList);
      if (dChannel != k2body) {
        continue;
      }
      registry.fill(HIST("hCandidateCounter"), 3);

      if (dauIDList[0] != mcDauTrack.globalIndex()) {
        continue;
      }
      registry.fill(HIST("hCandidateCounter"), 4);

      if (std::abs(dauTrack.tpcNSigmaAl()) > cutNSigmaAl) {
        continue;
      }
      registry.fill(HIST("hCandidateCounter"), 5);

      float invMass = RecoDecay::m(std::array{std::array{kinkCand.pxDaug(), kinkCand.pyDaug(), kinkCand.pzDaug()}, std::array{kinkCand.pxDaugNeut(), kinkCand.pyDaugNeut(), kinkCand.pzDaugNeut()}}, std::array{o2::constants::physics::MassAlpha, o2::constants::physics::MassPi0});
      registry.fill(HIST("hDiffSVx"), kinkCand.xDecVtx() - mcDauTrack.vx());
      registry.fill(HIST("hDiffSVy"), kinkCand.yDecVtx() - mcDauTrack.vy());
      registry.fill(HIST("hDiffSVz"), kinkCand.zDecVtx() - mcDauTrack.vz());
      registry.fill(HIST("hDiffDauPx"), kinkCand.pxDaug() - mcDauTrack.px());
      registry.fill(HIST("hDiffDauPy"), kinkCand.pyDaug() - mcDauTrack.py());
      registry.fill(HIST("hDiffDauPz"), kinkCand.pzDaug() - mcDauTrack.pz());
      registry.fill(HIST("h2MassHyperhelium4sigmaPt"), kinkCand.mothSign() * kinkCand.ptMoth(), invMass);
      registry.fill(HIST("h2NSigmaAlPt"), kinkCand.mothSign() * kinkCand.ptDaug(), dauTrack.tpcNSigmaAl());
    }
  }
  PROCESS_SWITCH(Hyperhelium4sigmaAnalysis, processMC, "process MC", false);
};

//--------------------------------------------------------------
// check the performance of mcparticle
struct Hyperhelium4sigmaQa {

  HistogramRegistry genQAHist{"genQAHist", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry recoQAHist{"recoQAHist", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  ConfigurableAxis ptBins{"ptBins", {200, 0.f, 10.f}, "Binning for #it{p}_{T} (GeV/#it{c})"};
  ConfigurableAxis ctBins{"ctBins", {100, 0.f, 25.f}, "Binning for c#it{t} (cm)"};
  ConfigurableAxis rigidityBins{"rigidityBins", {200, -10.f, 10.f}, "Binning for #it{p}_{T} (GeV/#it{c})"};
  ConfigurableAxis nsigmaBins{"nsigmaBins", {120, -6.f, 6.f}, "Binning for n sigma"};
  ConfigurableAxis invMassBins{"invMassBins", {100, 3.85f, 4.15f}, "Binning for invariant mass (GeV/#it{c}^{2})"};
  ConfigurableAxis radiusBins{"radiusBins", {40, 0.f, 40.f}, "Binning for radius in xy plane (cm)"};

  void init(InitContext&)
  {
    if (doprocessMC == true) {
      const AxisSpec ptAxis{ptBins, "#it{p}_{T} (GeV/#it{c})"};
      const AxisSpec ctAxis{ctBins, "c#it{t} (cm)"};
      const AxisSpec rigidityAxis{rigidityBins, "p/z (GeV/#it{c})"};
      const AxisSpec nsigmaAxis{nsigmaBins, "TPC n#sigma"};
      const AxisSpec invMassAxis{invMassBins, "Inv Mass (GeV/#it{c}^{2})"};
      const AxisSpec diffPtAxis{200, -10.f, 10.f, "#Delta p_{T} (GeV/#it{c})"};
      const AxisSpec radiusAxis{radiusBins, "R (cm)"};

      auto hCollCounter = genQAHist.add<TH1>("hCollCounter", "hCollCounter", HistType::kTH1F, {{2, 0.0f, 2.0f}});
      hCollCounter->GetXaxis()->SetBinLabel(1, "Reconstructed Collisions");
      hCollCounter->GetXaxis()->SetBinLabel(2, "Selected");
      auto hMcCollCounter = genQAHist.add<TH1>("hMcCollCounter", "hMcCollCounter", HistType::kTH1F, {{2, 0.0f, 2.0f}});
      hMcCollCounter->GetXaxis()->SetBinLabel(1, "MC Collisions");
      hMcCollCounter->GetXaxis()->SetBinLabel(2, "Reconstructed");

      auto hGenHyperHelium4SigmaCounter = genQAHist.add<TH1>("hGenHyperHelium4SigmaCounter", "", HistType::kTH1F, {{10, 0.f, 10.f}});
      hGenHyperHelium4SigmaCounter->GetXaxis()->SetBinLabel(1, "H4S All");
      hGenHyperHelium4SigmaCounter->GetXaxis()->SetBinLabel(2, "Matter");
      hGenHyperHelium4SigmaCounter->GetXaxis()->SetBinLabel(3, "AntiMatter");
      hGenHyperHelium4SigmaCounter->GetXaxis()->SetBinLabel(4, "#alpha + #pi^{0}");
      hGenHyperHelium4SigmaCounter->GetXaxis()->SetBinLabel(5, "#bar{#alpha} + #pi^{0}");
      hGenHyperHelium4SigmaCounter->GetXaxis()->SetBinLabel(6, "t + p + #pi^{0}");
      hGenHyperHelium4SigmaCounter->GetXaxis()->SetBinLabel(7, "#bar{t} + #bar{p} + #pi^{0}");
      hGenHyperHelium4SigmaCounter->GetXaxis()->SetBinLabel(8, "t + n + #pi^{+}");
      hGenHyperHelium4SigmaCounter->GetXaxis()->SetBinLabel(9, "#bar{t} + #bar{n} + #pi^{+}");
      hGenHyperHelium4SigmaCounter->GetXaxis()->SetBinLabel(10, "Unexpected");

      auto hEvtSelectedHyperHelium4SigmaCounter = genQAHist.add<TH1>("hEvtSelectedHyperHelium4SigmaCounter", "", HistType::kTH1F, {{2, 0.f, 2.f}});
      hEvtSelectedHyperHelium4SigmaCounter->GetXaxis()->SetBinLabel(1, "Generated");
      hEvtSelectedHyperHelium4SigmaCounter->GetXaxis()->SetBinLabel(2, "Survived");

      genQAHist.add<TH1>("hGenHyperHelium4SigmaP", "", HistType::kTH1F, {ptAxis});
      genQAHist.add<TH1>("hGenHyperHelium4SigmaPt", "", HistType::kTH1F, {ptAxis});
      genQAHist.add<TH1>("hGenHyperHelium4SigmaCt", "", HistType::kTH1F, {ctAxis});
      genQAHist.add<TH1>("hMcRecoInvMass", "", HistType::kTH1F, {invMassAxis});

      // efficiency/criteria studies for tracks which are true candidates
      hMotherCounter = recoQAHist.add<TH1>("hMotherCounter", "", HistType::kTH1F, {{9, 0.f, 9.f}});
      hMotherCounter->GetXaxis()->SetBinLabel(1, "Generated");
      hMotherCounter->GetXaxis()->SetBinLabel(2, "Reconstructed");
      hMotherCounter->GetXaxis()->SetBinLabel(3, "eta");
      hMotherCounter->GetXaxis()->SetBinLabel(4, "has collision");
      hMotherCounter->GetXaxis()->SetBinLabel(5, "ITSonly");
      hMotherCounter->GetXaxis()->SetBinLabel(6, "ITS hits");
      hMotherCounter->GetXaxis()->SetBinLabel(7, "ITS IR");
      hMotherCounter->GetXaxis()->SetBinLabel(8, "ITS chi2");
      hMotherCounter->GetXaxis()->SetBinLabel(9, "pt");
      recoQAHist.add<TH2>("hTrueMotherRVsDiffPt", ";#Delta p_{T} (GeV/#it{c});R (cm);", HistType::kTH2F, {diffPtAxis, radiusAxis});
      recoQAHist.add<TH2>("hTrueMotherRVsDiffPz", ";#Delta p_{z} (GeV/#it{c});R (cm);", HistType::kTH2F, {diffPtAxis, radiusAxis});
      recoQAHist.add<TH2>("hGoodMotherRVsDiffPt", ";#Delta p_{T} (GeV/#it{c});R (cm);", HistType::kTH2F, {diffPtAxis, radiusAxis});
      recoQAHist.add<TH2>("hGoodMotherRVsDiffPz", ";#Delta p_{z} (GeV/#it{c});R (cm);", HistType::kTH2F, {diffPtAxis, radiusAxis});

      hDauAlphaCounter = recoQAHist.add<TH1>("hDauAlphaCounter", "", HistType::kTH1F, {{7, 0.f, 7.f}});
      hDauTritonCounter = recoQAHist.add<TH1>("hDauTritonCounter", "", HistType::kTH1F, {{7, 0.f, 7.f}});
      hDauProtonCounter = recoQAHist.add<TH1>("hDauProtonCounter", "", HistType::kTH1F, {{7, 0.f, 7.f}});
      hDauPionCounter = recoQAHist.add<TH1>("hDauPionCounter", "", HistType::kTH1F, {{7, 0.f, 7.f}});

      recoQAHist.add<TH2>("hDauAlphaTPCNSigma", "", HistType::kTH2F, {rigidityAxis, nsigmaAxis});
      recoQAHist.add<TH2>("hDauTritonTPCNSigma", "", HistType::kTH2F, {rigidityAxis, nsigmaAxis});
      recoQAHist.add<TH2>("hDauProtonTPCNSigma", "", HistType::kTH2F, {rigidityAxis, nsigmaAxis});
      recoQAHist.add<TH2>("hDauPionTPCNSigma", "", HistType::kTH2F, {rigidityAxis, nsigmaAxis});
    }
  }

  Configurable<bool> skipRejectedEvents{"skipRejectedEvents", false, "Flag to skip events that fail event selection cuts"};
  Configurable<bool> mcEventCut{"mcEventCut", true, "flag to enable mc event selection: kIsTriggerTVX and kNoTimeFrameBorder"};
  Configurable<bool> eventPosZCut{"eventPosZCut", true, "flag to enable event posZ selection"};
  Configurable<float> maxPosZ{"maxPosZ", 10.f, "max pv posZ for event selection"};

  Configurable<float> etaMax{"etaMax", 1., "eta cut for tracks"};
  Configurable<float> minPtMoth{"minPtMoth", 0.5, "Minimum pT of the hypercandidate"};
  Configurable<float> tpcPidNsigmaCut{"tpcPidNsigmaCut", 5, "tpcPidNsigmaCut"};
  Configurable<float> nTPCClusMinDaug{"nTPCClusMinDaug", 80, "daug NTPC clusters cut"};
  Configurable<float> itsMaxChi2{"itsMaxChi2", 36, "max chi2 for ITS"};
  Configurable<float> minRatioTPCNCls{"minRatioTPCNCls", 0.8, "min ratio of TPC crossed rows to findable clusters"};

  Preslice<aod::McParticles> permcCollision = o2::aod::mcparticle::mcCollisionId;

  // construct index array from mcParticle to track
  std::vector<int64_t> mcPartIndices;

  template <typename TTrackTable>
  void setTrackIDForMC(aod::McParticles const& particlesMC, TTrackTable const& tracks)
  {
    mcPartIndices.clear();
    mcPartIndices.resize(particlesMC.size());
    std::fill(mcPartIndices.begin(), mcPartIndices.end(), -1);
    for (const auto& track : tracks) {
      if (track.has_mcParticle()) {
        auto mcparticle = track.template mcParticle_as<aod::McParticles>();
        if (mcPartIndices[mcparticle.globalIndex()] == -1) {
          mcPartIndices[mcparticle.globalIndex()] = track.globalIndex();
        } else {
          auto candTrack = tracks.rawIteratorAt(mcPartIndices[mcparticle.globalIndex()]);
          // Use the track which has innest information (also best quality?
          if (track.x() < candTrack.x()) {
            mcPartIndices[mcparticle.globalIndex()] = track.globalIndex();
          }
        }
      }
    }
  }

  // qa for mother track selection
  template <typename TTrack>
  bool motherTrackCheck(const TTrack& track, const std::shared_ptr<TH1> hist)
  {
    hist->Fill(1);

    if (std::abs(track.eta()) > etaMax) {
      return false;
    }
    hist->Fill(2);

    if (!track.has_collision()) {
      return false;
    }
    hist->Fill(3);

    if (!track.hasITS() || track.hasTPC() || track.hasTOF()) {
      return false;
    }
    hist->Fill(4);

    if (track.itsNCls() >= kITSLayers - 1) {
      return false;
    }
    hist->Fill(5);

    if (track.itsNClsInnerBarrel() != kITSInnerBarrelLayers) {
      return false;
    }
    hist->Fill(6);

    if (track.itsChi2NCl() >= itsMaxChi2) {
      return false;
    }
    hist->Fill(7);

    if (track.pt() <= minPtMoth) {
      return false;
    }
    hist->Fill(8);

    return true;
  }

  // qa for daughter track selection
  template <typename TTrack>
  void daughterTrackCheck(const TTrack& track, const std::shared_ptr<TH1> hist, float tpcNSigma)
  {
    hist->Fill(1);

    if (std::abs(track.eta()) > etaMax) {
      return;
    }
    hist->Fill(2);

    if (!track.hasITS() || !track.hasTPC()) {
      return;
    }
    hist->Fill(3);

    if (track.itsNClsInnerBarrel() != 0 && track.itsNCls() > kITSInnerBarrelLayers && track.tpcNClsCrossedRows() <= minRatioTPCNCls * track.tpcNClsFindable() && track.tpcNClsFound() <= nTPCClusMinDaug) {
      return;
    }
    hist->Fill(4);

    if (std::abs(tpcNSigma) > tpcPidNsigmaCut) {
      return;
    }
    hist->Fill(5);

    if (track.hasTOF()) {
      return;
    }
    hist->Fill(6);
  }

  void processData(o2::aod::Collisions const&)
  {
    // dummy process function;
  }
  PROCESS_SWITCH(Hyperhelium4sigmaQa, processData, "process data", true);

  void processMC(aod::McCollisions const& mcCollisions, aod::McParticles const& particlesMC, o2::soa::Join<o2::aod::Collisions, o2::aod::McCollisionLabels, o2::aod::EvSels> const& collisions, MCLabeledTracksIU const& tracks)
  {

    // check mcparticles
    setTrackIDForMC(particlesMC, tracks);
    std::vector<int64_t> selectedEvents(collisions.size());
    std::vector<int> dauIDList;
    int nevts = 0;
    for (const auto& collision : collisions) {
      genQAHist.fill(HIST("hCollCounter"), 0.5);
      if (mcEventCut && (!collision.selection_bit(aod::evsel::kIsTriggerTVX) || !collision.selection_bit(aod::evsel::kNoTimeFrameBorder))) {
        continue;
      }
      if (eventPosZCut && std::abs(collision.posZ()) > maxPosZ) { // 10cm
        continue;
      }
      genQAHist.fill(HIST("hCollCounter"), 1.5);
      selectedEvents[nevts++] = collision.mcCollision_as<aod::McCollisions>().globalIndex();
    }
    selectedEvents.resize(nevts);

    for (const auto& mcCollision : mcCollisions) {
      genQAHist.fill(HIST("hMcCollCounter"), 0.5);
      const auto evtReconstructedAndSelected = std::find(selectedEvents.begin(), selectedEvents.end(), mcCollision.globalIndex()) != selectedEvents.end();
      if (evtReconstructedAndSelected) { // Check that the event is reconstructed and that the reconstructed events pass the selection
        genQAHist.fill(HIST("hMcCollCounter"), 1.5);
      } else {
        if (skipRejectedEvents) {
          continue;
        }
      }

      const auto& dparticlesMC = particlesMC.sliceBy(permcCollision, mcCollision.globalIndex());

      for (const auto& mcparticle : dparticlesMC) {

        bool isMatter;
        if (mcparticle.pdgCode() == o2::constants::physics::Pdg::kHyperHelium4Sigma) {
          genQAHist.fill(HIST("hGenHyperHelium4SigmaCounter"), 1.5);
          isMatter = true;
        } else if (mcparticle.pdgCode() == -o2::constants::physics::Pdg::kHyperHelium4Sigma) {
          genQAHist.fill(HIST("hGenHyperHelium4SigmaCounter"), 2.5);
          isMatter = false;
        } else {
          continue;
        }

        genQAHist.fill(HIST("hGenHyperHelium4SigmaCounter"), 0.5);
        genQAHist.fill(HIST("hEvtSelectedHyperHelium4SigmaCounter"), 0.5);
        if (evtReconstructedAndSelected) {
          genQAHist.fill(HIST("hEvtSelectedHyperHelium4SigmaCounter"), 1.5);
        }

        auto dChannel = getDecayChannelH4S<aod::McParticles>(mcparticle, dauIDList);
        if (dChannel == kNDecayChannel) {
          genQAHist.fill(HIST("hGenHyperHelium4SigmaCounter"), 9.5);
          continue;
        }

        // bool isAllDauReconstructed = mcPartIndices[dauIDList[0]] != -1 && (dChannel == k2body ? true : mcPartIndices[dauIDList[1]] != -1);

        // check for mother tracks
        recoQAHist.fill(HIST("hMotherCounter"), 0);

        double svPos[3] = {-999, -999, -999};
        double dauAlphaMom[3] = {-999, -999, -999};
        double dauTritonMom[3] = {-999, -999, -999};
        double dauProtonMom[3] = {-999, -999, -999};
        double dauNeuteronMom[3] = {-999, -999, -999};
        double dauChargedPionMom[3] = {-999, -999, -999};
        double dauPion0Mom[3] = {-999, -999, -999};
        for (const auto& mcparticleDaughter : mcparticle.daughters_as<aod::McParticles>()) {
          if (std::abs(mcparticleDaughter.pdgCode()) == o2::constants::physics::Pdg::kAlpha) {
            dauAlphaMom[0] = mcparticleDaughter.px();
            dauAlphaMom[1] = mcparticleDaughter.py();
            dauAlphaMom[2] = mcparticleDaughter.pz();

            // get SV position for 2body decay
            svPos[0] = mcparticleDaughter.vx();
            svPos[1] = mcparticleDaughter.vy();
            svPos[2] = mcparticleDaughter.vz();

            recoQAHist.fill(HIST("hDauAlphaCounter"), 0);
            if (mcPartIndices[mcparticleDaughter.globalIndex()] != -1) {
              auto track = tracks.rawIteratorAt(mcPartIndices[mcparticleDaughter.globalIndex()]);
              recoQAHist.fill(HIST("hDauAlphaTPCNSigma"), track.p() * track.sign(), track.tpcNSigmaAl());
              daughterTrackCheck(track, hDauAlphaCounter, track.tpcNSigmaAl());
            }
          } else if (std::abs(mcparticleDaughter.pdgCode()) == o2::constants::physics::Pdg::kTriton) {
            dauTritonMom[0] = mcparticleDaughter.px();
            dauTritonMom[1] = mcparticleDaughter.py();
            dauTritonMom[2] = mcparticleDaughter.pz();

            // get SV position for 3body decay
            svPos[0] = mcparticleDaughter.vx();
            svPos[1] = mcparticleDaughter.vy();
            svPos[2] = mcparticleDaughter.vz();

            recoQAHist.fill(HIST("hDauTritonCounter"), 0);
            if (mcPartIndices[mcparticleDaughter.globalIndex()] != -1) {
              auto track = tracks.rawIteratorAt(mcPartIndices[mcparticleDaughter.globalIndex()]);
              recoQAHist.fill(HIST("hDauTritonTPCNSigma"), track.p() * track.sign(), track.tpcNSigmaTr());
              daughterTrackCheck(track, hDauTritonCounter, track.tpcNSigmaTr());
            }
          } else if (std::abs(mcparticleDaughter.pdgCode()) == PDG_t::kProton) {
            dauProtonMom[0] = mcparticleDaughter.px();
            dauProtonMom[1] = mcparticleDaughter.py();
            dauProtonMom[2] = mcparticleDaughter.pz();

            recoQAHist.fill(HIST("hDauProtonCounter"), 0);
            if (mcPartIndices[mcparticleDaughter.globalIndex()] != -1) {
              auto track = tracks.rawIteratorAt(mcPartIndices[mcparticleDaughter.globalIndex()]);
              recoQAHist.fill(HIST("hDauProtonTPCNSigma"), track.p() * track.sign(), track.tpcNSigmaPr());
              daughterTrackCheck(track, hDauProtonCounter, track.tpcNSigmaPr());
            }
          } else if (std::abs(mcparticleDaughter.pdgCode()) == PDG_t::kNeutron) {
            dauNeuteronMom[0] = mcparticleDaughter.px();
            dauNeuteronMom[1] = mcparticleDaughter.py();
            dauNeuteronMom[2] = mcparticleDaughter.pz();
          } else if (std::abs(mcparticleDaughter.pdgCode()) == PDG_t::kPiPlus) {
            dauChargedPionMom[0] = mcparticleDaughter.px();
            dauChargedPionMom[1] = mcparticleDaughter.py();
            dauChargedPionMom[2] = mcparticleDaughter.pz();

            recoQAHist.fill(HIST("hDauPionCounter"), 0);
            if (mcPartIndices[mcparticleDaughter.globalIndex()] != -1) {
              auto track = tracks.rawIteratorAt(mcPartIndices[mcparticleDaughter.globalIndex()]);
              recoQAHist.fill(HIST("hDauPionTPCNSigma"), track.p() * track.sign(), track.tpcNSigmaPi());
              daughterTrackCheck(track, hDauPionCounter, track.tpcNSigmaPi());
            }
          } else if (mcparticleDaughter.pdgCode() == PDG_t::kPi0) {
            dauPion0Mom[0] = mcparticleDaughter.px();
            dauPion0Mom[1] = mcparticleDaughter.py();
            dauPion0Mom[2] = mcparticleDaughter.pz();
          }
        }

        genQAHist.fill(HIST("hGenHyperHelium4SigmaP"), mcparticle.p());
        genQAHist.fill(HIST("hGenHyperHelium4SigmaPt"), mcparticle.pt());
        double ct = RecoDecay::sqrtSumOfSquares(svPos[0] - mcparticle.vx(), svPos[1] - mcparticle.vy(), svPos[2] - mcparticle.vz()) * o2::constants::physics::MassHyperHelium4Sigma / mcparticle.p();
        genQAHist.fill(HIST("hGenHyperHelium4SigmaCt"), ct);

        if (mcPartIndices[mcparticle.globalIndex()] != -1) {
          auto motherTrack = tracks.rawIteratorAt(mcPartIndices[mcparticle.globalIndex()]);
          bool isGoodMother = motherTrackCheck(motherTrack, hMotherCounter);
          double svR = RecoDecay::sqrtSumOfSquares(svPos[0], svPos[1]);
          recoQAHist.fill(HIST("hTrueMotherRVsDiffPt"), mcparticle.pt() - 2 * motherTrack.pt(), svR);
          recoQAHist.fill(HIST("hTrueMotherRVsDiffPz"), mcparticle.pz() - 2 * motherTrack.pz(), svR);
          if (isGoodMother) {
            recoQAHist.fill(HIST("hGoodMotherRVsDiffPt"), mcparticle.pt() - 2 * motherTrack.pt(), svR);
            recoQAHist.fill(HIST("hGoodMotherRVsDiffPz"), mcparticle.pz() - 2 * motherTrack.pz(), svR);
          }
        }

        if (dChannel == k2body) {
          if (isMatter) {
            genQAHist.fill(HIST("hGenHyperHelium4SigmaCounter"), 3.5);
          } else {
            genQAHist.fill(HIST("hGenHyperHelium4SigmaCounter"), 4.5);
          }
          double hyperHelium4SigmaMCMass = RecoDecay::m(std::array{std::array{dauAlphaMom[0], dauAlphaMom[1], dauAlphaMom[2]}, std::array{dauPion0Mom[0], dauPion0Mom[1], dauPion0Mom[2]}}, std::array{o2::constants::physics::MassAlpha, o2::constants::physics::MassPi0});
          genQAHist.fill(HIST("hMcRecoInvMass"), hyperHelium4SigmaMCMass);
        } else if (dChannel == k3body_p) {
          if (isMatter) {
            genQAHist.fill(HIST("hGenHyperHelium4SigmaCounter"), 5.5);
          } else {
            genQAHist.fill(HIST("hGenHyperHelium4SigmaCounter"), 6.5);
          }
          double hyperHelium4SigmaMCMass = RecoDecay::m(std::array{std::array{dauTritonMom[0], dauTritonMom[1], dauTritonMom[2]}, std::array{dauProtonMom[0], dauProtonMom[1], dauProtonMom[2]}, std::array{dauPion0Mom[0], dauPion0Mom[1], dauPion0Mom[2]}}, std::array{o2::constants::physics::MassTriton, o2::constants::physics::MassProton, o2::constants::physics::MassPi0});
          genQAHist.fill(HIST("hMcRecoInvMass"), hyperHelium4SigmaMCMass);
        } else if (dChannel == k3body_n) {
          if (isMatter) {
            genQAHist.fill(HIST("hGenHyperHelium4SigmaCounter"), 7.5);
          } else {
            genQAHist.fill(HIST("hGenHyperHelium4SigmaCounter"), 8.5);
          }
          double hyperHelium4SigmaMCMass = RecoDecay::m(std::array{std::array{dauTritonMom[0], dauTritonMom[1], dauTritonMom[2]}, std::array{dauNeuteronMom[0], dauNeuteronMom[1], dauNeuteronMom[2]}, std::array{dauChargedPionMom[0], dauChargedPionMom[1], dauChargedPionMom[2]}}, std::array{o2::constants::physics::MassTriton, o2::constants::physics::MassNeutron, o2::constants::physics::MassPionCharged});
          genQAHist.fill(HIST("hMcRecoInvMass"), hyperHelium4SigmaMCMass);
        }
      }
    }
  }
  PROCESS_SWITCH(Hyperhelium4sigmaQa, processMC, "do QA for MC prodcutions", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<Hyperhelium4sigmaAnalysis>(cfgc),
    adaptAnalysisTask<Hyperhelium4sigmaQa>(cfgc),
  };
}
