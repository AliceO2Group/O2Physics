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

/// \author Dushmanta Sahu (dushmanta.sahu@cern.ch)
/// \file multiplicityPt.cxx
/// \brief Analysis to do PID with MC - Full correction factors for pions, kaons, protons

#include "PWGLF/DataModel/mcCentrality.h"
#include "PWGLF/Utils/inelGt.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/MathConstants.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/EndOfStreamContext.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/runDataProcessing.h>

#include <TF1.h>
#include <TPDGCode.h>
#include <TString.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::math;
using namespace constants::physics;

using BCsRun3 = soa::Join<aod::BCs, aod::Timestamps, aod::BcSels,
                          aod::Run3MatchedToBCSparse>;

struct MultiplicityPt {

  Service<o2::framework::O2DatabasePDG> pdg;
  Service<ccdb::BasicCCDBManager> ccdb;
  static constexpr int CentBinMax = 100;
  static constexpr int MultBinMax = 200;
  static constexpr int RecMultBinMax = 100;

  enum INELCutSelection : int {
    INEL = 0,
    INELgt0 = 1,
    INELgt1 = 2
  };

  Configurable<bool> isRun3{"isRun3", true, "is Run3 dataset"};
  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};
  Configurable<int> cfgINELCut{"cfgINELCut", 0, "INEL event selection: 0 no sel, 1 INEL>0, 2 INEL>1"};
  Configurable<float> cfgCutEtaMax{"cfgCutEtaMax", 0.8f, "Max eta range for tracks"};
  Configurable<float> cfgCutEtaMin{"cfgCutEtaMin", -0.8f, "Min eta range for tracks"};
  Configurable<float> cfgCutNsigma{"cfgCutNsigma", 3.0f, "nsigma cut range for tracks"};

  Configurable<bool> enablePIDHistograms{"enablePIDHistograms", true, "Enable PID histograms"};
  Configurable<bool> useCustomTrackCuts{"useCustomTrackCuts", true, "Flag to use custom track cuts"};
  Configurable<int> itsPattern{"itsPattern", 0, "0 = Run3ITSibAny, 1 = Run3ITSallAny, 2 = Run3ITSall7Layers, 3 = Run3ITSibTwo"};
  Configurable<bool> requireITS{"requireITS", true, "Additional cut on the ITS requirement"};
  Configurable<bool> requireTPC{"requireTPC", true, "Additional cut on the TPC requirement"};
  Configurable<bool> requireGoldenChi2{"requireGoldenChi2", true, "Additional cut on the GoldenChi2"};
  Configurable<float> minNCrossedRowsTPC{"minNCrossedRowsTPC", 70.f, "Additional cut on the minimum number of crossed rows in the TPC"};
  Configurable<float> minNCrossedRowsOverFindableClustersTPC{"minNCrossedRowsOverFindableClustersTPC", 0.8f, "Additional cut on the minimum value of the ratio between crossed rows and findable clusters in the TPC"};
  Configurable<float> maxChi2PerClusterTPC{"maxChi2PerClusterTPC", 4.f, "Additional cut on the maximum value of the chi2 per cluster in the TPC"};
  Configurable<float> minChi2PerClusterTPC{"minChi2PerClusterTPC", 0.5f, "Additional cut on the minimum value of the chi2 per cluster in the TPC"};
  Configurable<float> maxChi2PerClusterITS{"maxChi2PerClusterITS", 36.f, "Additional cut on the maximum value of the chi2 per cluster in the ITS"};
  Configurable<float> maxDcaXYFactor{"maxDcaXYFactor", 1.f, "Additional cut on the maximum value of the DCA xy (multiplicative factor)"};
  Configurable<float> maxDcaZ{"maxDcaZ", 0.1f, "Additional cut on the maximum value of the DCA z"};
  Configurable<float> minTPCNClsFound{"minTPCNClsFound", 70.0f, "min number of found TPC clusters"};
  Configurable<float> minTPCNClsPID{"minTPCNClsPID", 130.0f, "min number of PID TPC clusters"};
  Configurable<bool> nClTPCFoundCut{"nClTPCFoundCut", false, "Apply TPC found clusters cut"};
  Configurable<bool> nClTPCPIDCut{"nClTPCPIDCut", true, "Apply TPC clusters for PID cut"};

  Configurable<bool> applyPhiCut{"applyPhiCut", false, "Apply phi sector cut to remove problematic TPC regions"};
  Configurable<float> pTthresholdPhiCut{"pTthresholdPhiCut", 2.0f, "pT threshold above which to apply phi cut"};
  Configurable<double> phiCutLowParam1{"phiCutLowParam1", 0.119297, "First parameter for low phi cut"};
  Configurable<double> phiCutLowParam2{"phiCutLowParam2", 0.000379693, "Second parameter for low phi cut"};
  Configurable<double> phiCutHighParam1{"phiCutHighParam1", 0.16685, "First parameter for high phi cut"};
  Configurable<double> phiCutHighParam2{"phiCutHighParam2", 0.00981942, "Second parameter for high phi cut"};

  Configurable<float> cfgTrkLowPtCut{"cfgTrkLowPtCut", 0.15f, "Minimum constituent pT"};

  Configurable<float> cfgCutNsigmaPi{"cfgCutNsigmaPi", 3.0f, "nsigma cut for pions"};
  Configurable<float> cfgCutNsigmaKa{"cfgCutNsigmaKa", 2.5f, "nsigma cut for kaons"};
  Configurable<float> cfgCutNsigmaPr{"cfgCutNsigmaPr", 2.5f, "nsigma cut for protons"};

  TrackSelection customTrackCuts;
  TF1* fphiCutLow = nullptr;
  TF1* fphiCutHigh = nullptr;

  HistogramRegistry ue;

  using CollisionTableData = soa::Join<aod::Collisions, aod::EvSels, aod::McCentFT0Ms>;
  using TrackTableData = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection,
                                   aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr>;
  using TrackTableMC = soa::Join<aod::Tracks, aod::McTrackLabels, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection,
                                 aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr>;
  using ParticlesMC = aod::McParticles;
  using McCollisions = aod::McCollisions;
  using RecoCollisions = aod::Collisions;

  Preslice<aod::McParticles> perMCCol = aod::mcparticle::mcCollisionId;

  enum ParticleSpecies : int {
    PartPion = 0,
    PartKaon = 1,
    PartProton = 2,
  };

  int getMagneticField(uint64_t timestamp)
  {
    static o2::parameters::GRPMagField* grpo = nullptr;
    if (grpo == nullptr) {
      grpo = ccdb->getForTimeStamp<o2::parameters::GRPMagField>("GLO/Config/GRPMagField", timestamp);
      if (grpo == nullptr) {
        LOGF(fatal, "GRP object not found for timestamp %llu", timestamp);
        return 0;
      }
      LOGF(info, "Retrieved GRP for timestamp %llu with magnetic field of %d kG", timestamp, grpo->getNominalL3Field());
    }
    return grpo->getNominalL3Field();
  }

  float getTransformedPhi(const float phi, const int charge, const float magField) const
  {
    float transformedPhi = phi;
    if (magField < 0)
      transformedPhi = o2::constants::math::TwoPI - transformedPhi;
    if (charge < 0)
      transformedPhi = o2::constants::math::TwoPI - transformedPhi;
    transformedPhi += o2::constants::math::PI / 18.0f;
    transformedPhi = std::fmod(transformedPhi, o2::constants::math::PI / 9.0f);
    return transformedPhi;
  }

  template <typename TrackType>
  bool passedPhiCut(const TrackType& track, float magField) const
  {
    if (!applyPhiCut.value)
      return true;
    if (track.pt() < pTthresholdPhiCut.value)
      return true;

    float phi = track.phi();
    int charge = track.sign();

    if (magField < 0)
      phi = o2::constants::math::TwoPI - phi;
    if (charge < 0)
      phi = o2::constants::math::TwoPI - phi;
    phi += o2::constants::math::PI / 18.0f;
    phi = std::fmod(phi, o2::constants::math::PI / 9.0f);

    if (phi < fphiCutHigh->Eval(track.pt()) && phi > fphiCutLow->Eval(track.pt()))
      return false;
    return true;
  }

  template <typename ParticleContainer>
  int countGeneratedChargedPrimaries(const ParticleContainer& particles, float etaMax, float ptMin) const
  {
    int count = 0;
    for (const auto& particle : particles) {
      auto pdgParticle = pdg->GetParticle(particle.pdgCode());
      if (!pdgParticle || pdgParticle->Charge() == 0.)
        continue;
      if (!particle.isPhysicalPrimary())
        continue;
      if (std::abs(particle.eta()) > etaMax)
        continue;
      if (particle.pt() < ptMin)
        continue;
      count++;
    }
    return count;
  }

  template <typename T>
  bool passedNClTPCFoundCut(const T& trk) const
  {
    return !nClTPCFoundCut.value || trk.tpcNClsFound() >= minTPCNClsFound.value;
  }

  template <typename T>
  bool passedNClTPCPIDCut(const T& trk) const
  {
    return !nClTPCPIDCut.value || trk.tpcNClsPID() >= minTPCNClsPID.value;
  }

  template <typename TrackType>
  bool passesCutWoDCA(TrackType const& track) const
  {
    if (useCustomTrackCuts.value) {
      for (int i = 0; i < static_cast<int>(TrackSelection::TrackCuts::kNCuts); i++) {
        if (i == static_cast<int>(TrackSelection::TrackCuts::kDCAxy) ||
            i == static_cast<int>(TrackSelection::TrackCuts::kDCAz))
          continue;
        if (!customTrackCuts.IsSelected(track, static_cast<TrackSelection::TrackCuts>(i)))
          return false;
      }
      return true;
    }
    return track.isGlobalTrackWoDCA();
  }

  template <typename TrackType>
  bool passesDCAxyCut(TrackType const& track) const
  {
    if (useCustomTrackCuts.value) {
      if (!passesCutWoDCA(track))
        return false;
      constexpr float DcaXYConst = 0.0105f;
      constexpr float DcaXYPtScale = 0.0350f;
      constexpr float DcaXYPtPower = 1.1f;
      const float maxDcaXY = maxDcaXYFactor.value * (DcaXYConst + DcaXYPtScale / std::pow(track.pt(), DcaXYPtPower));
      return std::abs(track.dcaXY()) <= maxDcaXY;
    }
    return track.isGlobalTrack();
  }

  template <typename TrackType>
  bool passesTrackSelection(TrackType const& track, float magField = 0) const
  {
    if (track.eta() < cfgCutEtaMin.value || track.eta() > cfgCutEtaMax.value)
      return false;
    if (track.tpcChi2NCl() < minChi2PerClusterTPC.value || track.tpcChi2NCl() > maxChi2PerClusterTPC.value)
      return false;
    if (!passesCutWoDCA(track))
      return false;
    if (!passesDCAxyCut(track))
      return false;
    if (!passedNClTPCFoundCut(track))
      return false;
    if (!passedNClTPCPIDCut(track))
      return false;
    if (!passedPhiCut(track, magField))
      return false;
    return true;
  }

  template <typename TrackType>
  bool passesPIDSelection(const TrackType& track, int species) const
  {
    float nsigmaTPC = 0.f;
    float cutValue = 0.f;

    if (species == PartPion) {
      nsigmaTPC = track.tpcNSigmaPi();
      cutValue = cfgCutNsigmaPi.value;
    } else if (species == PartKaon) {
      nsigmaTPC = track.tpcNSigmaKa();
      cutValue = cfgCutNsigmaKa.value;
    } else if (species == PartProton) {
      nsigmaTPC = track.tpcNSigmaPr();
      cutValue = cfgCutNsigmaPr.value;
    }

    return std::abs(nsigmaTPC) < cutValue;
  }

  template <typename ParticleType>
  bool isGoodPrimary(ParticleType const& particle) const
  {
    auto pdgParticle = pdg->GetParticle(particle.pdgCode());
    if (!pdgParticle || pdgParticle->Charge() == 0.)
      return false;
    if (!particle.isPhysicalPrimary())
      return false;
    if (std::abs(particle.eta()) >= cfgCutEtaMax.value)
      return false;
    if (particle.pt() < cfgTrkLowPtCut.value)
      return false;
    return true;
  }

  void processData(CollisionTableData::iterator const& collision, TrackTableData const& tracks, BCsRun3 const& bcs);
  PROCESS_SWITCH(MultiplicityPt, processData, "process data", false);

  void processMC(TrackTableMC const& tracks, aod::McParticles const& particles, aod::McCollisions const& mcCollisions,
                 RecoCollisions const& collisions, aod::McCollisionLabels const& labels, aod::McCentFT0Ms const& centTable, BCsRun3 const& bcs);
  PROCESS_SWITCH(MultiplicityPt, processMC, "process MC", true);

  void init(InitContext const&);
  void endOfStream(EndOfStreamContext& /*eos*/) {}
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<MultiplicityPt>(cfgc)};
}

void MultiplicityPt::init(InitContext const&)
{
  LOG(info) << "==================================================";
  LOG(info) << "Initializing MultiplicityPt task - FULL CORRECTION FACTORS";
  LOG(info) << "==================================================";

  if (applyPhiCut.value) {
    fphiCutLow = new TF1("StandardPhiCutLow", Form("%f/x/x+pi/18.0-%f", phiCutLowParam1.value, phiCutLowParam2.value), 0, 50);
    fphiCutHigh = new TF1("StandardPhiCutHigh", Form("%f/x+pi/18.0+%f", phiCutHighParam1.value, phiCutHighParam2.value), 0, 50);
  }

  if (useCustomTrackCuts.value) {
    customTrackCuts = getGlobalTrackSelectionRun3ITSMatch(itsPattern.value);
    customTrackCuts.SetRequireITSRefit(requireITS.value);
    customTrackCuts.SetRequireTPCRefit(requireTPC.value);
    customTrackCuts.SetRequireGoldenChi2(requireGoldenChi2.value);
    customTrackCuts.SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC.value);
    customTrackCuts.SetMaxChi2PerClusterITS(maxChi2PerClusterITS.value);
    customTrackCuts.SetMinNCrossedRowsTPC(minNCrossedRowsTPC.value);
    customTrackCuts.SetMinNClustersTPC(minTPCNClsFound.value);
    customTrackCuts.SetMinNCrossedRowsOverFindableClustersTPC(minNCrossedRowsOverFindableClustersTPC.value);
    customTrackCuts.SetMaxDcaXYPtDep([](float /*pt*/) { return 10000.f; });
    customTrackCuts.SetMaxDcaZ(maxDcaZ.value);
    customTrackCuts.print();
  }

  // Axis definitions
  ConfigurableAxis ptBinning{"ptBinning", {VARIABLE_WIDTH, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 18.0, 20.0}, "pT bin limits"};
  AxisSpec ptAxis = {ptBinning, "#it{p}_{T} (GeV/#it{c})"};

  std::vector<double> centBinningStd = {0., 1., 5., 10., 15., 20., 30., 40., 50., 60., 70., 80., 90., 100.};
  AxisSpec centAxis = {centBinningStd, "FT0M Centrality (%)"};

  std::vector<double> centBinningFine;
  for (int i = 0; i <= CentBinMax; i++)
    centBinningFine.push_back(static_cast<double>(i));
  AxisSpec centFineAxis = {centBinningFine, "FT0M Centrality (%)"};

  std::vector<double> multBins;
  for (int i = 0; i <= MultBinMax; i++)
    multBins.push_back(static_cast<double>(i));
  AxisSpec multAxis = {multBins, "N_{ch}^{gen} (|#eta|<0.8)"};

  std::vector<double> recoMultBins;
  for (int i = 0; i <= RecMultBinMax; i++)
    recoMultBins.push_back(static_cast<double>(i));
  AxisSpec recoMultAxis = {recoMultBins, "N_{ch}^{reco}"};

  ue.add("Centrality/hCentRaw", "Raw FT0M Centrality;Centrality (%);Counts", HistType::kTH1D, {centFineAxis});
  ue.add("Centrality/hCentAfterVtx", "Centrality after vertex cut;Centrality (%);Counts", HistType::kTH1D, {centFineAxis});
  ue.add("Centrality/hCentAfterAll", "Centrality after all cuts;Centrality (%);Counts", HistType::kTH1D, {centFineAxis});
  ue.add("Centrality/hCentVsMult", "Centrality vs Generated Multiplicity;Centrality (%);N_{ch}^{gen}", HistType::kTH2D, {centFineAxis, multAxis});
  ue.add("Centrality/hMultVsCent", "Generated Multiplicity vs Centrality;N_{ch}^{gen};Centrality (%)", HistType::kTH2D, {multAxis, centFineAxis});
  ue.add("Centrality/hRecoMultVsCent", "Reconstructed Track Multiplicity vs Centrality;Centrality (%);N_{tracks}^{reco}", HistType::kTH2D, {centFineAxis, recoMultAxis});

  ue.add("CutFlow/hCutStats", "Cut Statistics;Cut Stage;Counts", HistType::kTH1D, {{5, 0.5, 5.5}});
  auto hCut = ue.get<TH1>(HIST("CutFlow/hCutStats"));
  hCut->GetXaxis()->SetBinLabel(1, "All collisions");
  hCut->GetXaxis()->SetBinLabel(2, "Has MC match");
  hCut->GetXaxis()->SetBinLabel(3, "Has centrality");
  hCut->GetXaxis()->SetBinLabel(4, "Pass vertex");
  hCut->GetXaxis()->SetBinLabel(5, "Selected");

  ue.add("hEventLossBreakdown", "Event loss breakdown", HistType::kTH1D, {{3, 0.5, 3.5}});
  auto hLoss = ue.get<TH1>(HIST("hEventLossBreakdown"));
  hLoss->GetXaxis()->SetBinLabel(1, "Physics selected");
  hLoss->GetXaxis()->SetBinLabel(2, "Reconstructed");
  hLoss->GetXaxis()->SetBinLabel(3, "Selected");

  ue.add("MC/GenRecoCollisions", "Generated and Reconstructed MC Collisions", HistType::kTH1D, {{5, 0.5, 5.5}});
  auto hColl = ue.get<TH1>(HIST("MC/GenRecoCollisions"));
  hColl->GetXaxis()->SetBinLabel(1, "Collisions generated");
  hColl->GetXaxis()->SetBinLabel(2, "INEL>0");
  hColl->GetXaxis()->SetBinLabel(3, "INEL>1");
  hColl->GetXaxis()->SetBinLabel(4, "Reconstructed");
  hColl->GetXaxis()->SetBinLabel(5, "Selected");

  ue.add("MC/EventLoss/GenMultVsCent", "Generated charged particles vs FT0M centrality;FT0M Centrality (%);N_{ch}^{gen}", HistType::kTH2D, {centAxis, multAxis});
  ue.add("MC/EventLoss/GenMultVsCent_Selected", "Generated vs FT0M centrality (selected events);FT0M Centrality (%);N_{ch}^{gen}", HistType::kTH2D, {centAxis, multAxis});

  ue.add("Inclusive/hPtGen", "Generated primaries (physics selected);#it{p}_{T} (GeV/#it{c});Counts", HistType::kTH1D, {ptAxis});
  ue.add("Inclusive/hPtReco", "All reconstructed tracks (track selection only);#it{p}_{T} (GeV/#it{c});Counts", HistType::kTH1D, {ptAxis});
  ue.add("Inclusive/hPtPrimReco", "Reconstructed primaries (MC matched);#it{p}_{T} (GeV/#it{c});Counts", HistType::kTH1D, {ptAxis});
  ue.add("Inclusive/hPtSecReco", "Reconstructed secondaries (MC matched);#it{p}_{T} (GeV/#it{c});Counts", HistType::kTH1D, {ptAxis});

  ue.add("Pion/hPtPrimGenAll", "All generated #pi^{#pm} (no cuts);#it{p}_{T} (GeV/#it{c});Counts", HistType::kTH1D, {ptAxis});
  ue.add("Pion/hPtGenINEL", "Generated #pi^{#pm} in INEL>0 events (|#eta|<0.8, p_{T}>0.15);#it{p}_{T} (GeV/#it{c});Counts", HistType::kTH1D, {ptAxis});
  ue.add("Pion/hPtGenRecoEvent", "Generated #pi^{#pm} in reconstructed events (|#eta|<0.8, p_{T}>0.15);#it{p}_{T} (GeV/#it{c});Counts", HistType::kTH1D, {ptAxis});
  ue.add("Pion/hPtGen", "Generated #pi^{#pm} (physics selected);#it{p}_{T} (GeV/#it{c});Counts", HistType::kTH1D, {ptAxis});
  ue.add("Pion/hPtReco", "Reconstructed #pi^{#pm} (MC matched, any status);#it{p}_{T} (GeV/#it{c});Counts", HistType::kTH1D, {ptAxis});
  ue.add("Pion/hPtPrimReco", "Reconstructed primary #pi^{#pm} (MC matched);#it{p}_{T} (GeV/#it{c});Counts", HistType::kTH1D, {ptAxis});
  ue.add("Pion/hPtMeasured", "Measured #pi^{#pm} (PID selected);#it{p}_{T} (GeV/#it{c});Counts", HistType::kTH1D, {ptAxis});

  ue.add("Pion/hPtGenRecoEvent_Mult", "Generated #pi^{#pm} in reconstructed events vs multiplicity;#it{p}_{T} (GeV/#it{c});N_{ch}^{gen}", HistType::kTH2D, {ptAxis, multAxis});
  ue.add("Pion/hPtGenINEL_Mult", "Generated #pi^{#pm} in INEL>0 events vs multiplicity;#it{p}_{T} (GeV/#it{c});N_{ch}^{gen}", HistType::kTH2D, {ptAxis, multAxis});

  if (enablePIDHistograms) {
    ue.add("Pion/hNsigmaTPC", "TPC n#sigma #pi^{#pm};#it{p}_{T} (GeV/#it{c});n#sigma_{TPC}", HistType::kTH2D, {ptAxis, {200, -10, 10}});
  }

  ue.add("Kaon/hPtPrimGenAll", "All generated K^{#pm} (no cuts);#it{p}_{T} (GeV/#it{c});Counts", HistType::kTH1D, {ptAxis});
  ue.add("Kaon/hPtGenINEL", "Generated K^{#pm} in INEL>0 events (|#eta|<0.8, p_{T}>0.15);#it{p}_{T} (GeV/#it{c});Counts", HistType::kTH1D, {ptAxis});
  ue.add("Kaon/hPtGenRecoEvent", "Generated K^{#pm} in reconstructed events (|#eta|<0.8, p_{T}>0.15);#it{p}_{T} (GeV/#it{c});Counts", HistType::kTH1D, {ptAxis});
  ue.add("Kaon/hPtGen", "Generated K^{#pm} (physics selected);#it{p}_{T} (GeV/#it{c});Counts", HistType::kTH1D, {ptAxis});
  ue.add("Kaon/hPtReco", "Reconstructed K^{#pm} (MC matched, any status);#it{p}_{T} (GeV/#it{c});Counts", HistType::kTH1D, {ptAxis});
  ue.add("Kaon/hPtPrimReco", "Reconstructed primary K^{#pm} (MC matched);#it{p}_{T} (GeV/#it{c});Counts", HistType::kTH1D, {ptAxis});
  ue.add("Kaon/hPtMeasured", "Measured K^{#pm} (PID selected);#it{p}_{T} (GeV/#it{c});Counts", HistType::kTH1D, {ptAxis});

  ue.add("Kaon/hPtGenRecoEvent_Mult", "Generated K^{#pm} in reconstructed events vs multiplicity;#it{p}_{T} (GeV/#it{c});N_{ch}^{gen}", HistType::kTH2D, {ptAxis, multAxis});
  ue.add("Kaon/hPtGenINEL_Mult", "Generated K^{#pm} in INEL>0 events vs multiplicity;#it{p}_{T} (GeV/#it{c});N_{ch}^{gen}", HistType::kTH2D, {ptAxis, multAxis});

  if (enablePIDHistograms) {
    ue.add("Kaon/hNsigmaTPC", "TPC n#sigma K^{#pm};#it{p}_{T} (GeV/#it{c});n#sigma_{TPC}", HistType::kTH2D, {ptAxis, {200, -10, 10}});
  }

  ue.add("Proton/hPtPrimGenAll", "All generated p+#bar{p} (no cuts);#it{p}_{T} (GeV/#it{c});Counts", HistType::kTH1D, {ptAxis});
  ue.add("Proton/hPtGenINEL", "Generated p+#bar{p} in INEL>0 events (|#eta|<0.8, p_{T}>0.15);#it{p}_{T} (GeV/#it{c});Counts", HistType::kTH1D, {ptAxis});
  ue.add("Proton/hPtGenRecoEvent", "Generated p+#bar{p} in reconstructed events (|#eta|<0.8, p_{T}>0.15);#it{p}_{T} (GeV/#it{c});Counts", HistType::kTH1D, {ptAxis});
  ue.add("Proton/hPtGen", "Generated p+#bar{p} (physics selected);#it{p}_{T} (GeV/#it{c});Counts", HistType::kTH1D, {ptAxis});
  ue.add("Proton/hPtReco", "Reconstructed p+#bar{p} (MC matched, any status);#it{p}_{T} (GeV/#it{c});Counts", HistType::kTH1D, {ptAxis});
  ue.add("Proton/hPtPrimReco", "Reconstructed primary p+#bar{p} (MC matched);#it{p}_{T} (GeV/#it{c});Counts", HistType::kTH1D, {ptAxis});
  ue.add("Proton/hPtMeasured", "Measured p+#bar{p} (PID selected);#it{p}_{T} (GeV/#it{c});Counts", HistType::kTH1D, {ptAxis});

  ue.add("Proton/hPtGenRecoEvent_Mult", "Generated p+#bar{p} in reconstructed events vs multiplicity;#it{p}_{T} (GeV/#it{c});N_{ch}^{gen}", HistType::kTH2D, {ptAxis, multAxis});
  ue.add("Proton/hPtGenINEL_Mult", "Generated p+#bar{p} in INEL>0 events vs multiplicity;#it{p}_{T} (GeV/#it{c});N_{ch}^{gen}", HistType::kTH2D, {ptAxis, multAxis});

  if (enablePIDHistograms) {
    ue.add("Proton/hNsigmaTPC", "TPC n#sigma p+#bar{p};#it{p}_{T} (GeV/#it{c});n#sigma_{TPC}", HistType::kTH2D, {ptAxis, {200, -10, 10}});
  }

  ue.add("hEventsReco_Cent", "Reconstructed events vs centrality;FT0M Centrality (%);Counts", HistType::kTH1D, {centAxis});
  ue.add("hEventsINEL_Cent", "INEL>0 events vs centrality;FT0M Centrality (%);Counts", HistType::kTH1D, {centAxis});

  ue.add("hNclFoundTPC", "Number of TPC found clusters", HistType::kTH1D, {{200, 0, 200}});
  ue.add("hNclPIDTPC", "Number of TPC PID clusters", HistType::kTH1D, {{200, 0, 200}});
  ue.add("hEta", "Track eta;#eta;Counts", HistType::kTH1D, {{20, -0.8, 0.8}});
  ue.add("hPhi", "Track phi;#varphi (rad);Counts", HistType::kTH1D, {{64, 0, TwoPI}});
  ue.add("hvtxZ", "Vertex Z (data);Vertex Z (cm);Events", HistType::kTH1F, {{40, -20.0, 20.0}});

  LOG(info) << "=== Initialization complete ===";
}

void MultiplicityPt::processMC(TrackTableMC const& tracks,
                               aod::McParticles const& particles,
                               aod::McCollisions const& mcCollisions,
                               RecoCollisions const& collisions,
                               aod::McCollisionLabels const& labels,
                               aod::McCentFT0Ms const& centTable,
                               BCsRun3 const& /*bcs*/)
{

  std::map<int64_t, int> mcCollisionToNch;
  std::set<int64_t> physicsSelectedMCCollisions;
  std::map<int64_t, int> mcCollisionToINELClass;
  std::set<int64_t> inel0MCCollisions;
  std::map<int64_t, float> mcCollToCentFromReco; // MC collision ID -> centrality from reco

  ue.fill(HIST("MC/GenRecoCollisions"), 1.f, mcCollisions.size());

  for (const auto& mcCollision : mcCollisions) {
    int64_t mcCollId = mcCollision.globalIndex();
    auto particlesInCollision = particles.sliceBy(perMCCol, mcCollId);

    // Check if event has at least 1 particle in acceptance (INEL>0)
    bool hasParticleInAcceptance = false;
    for (const auto& particle : particlesInCollision) {
      auto pdgParticle = pdg->GetParticle(particle.pdgCode());
      if (!pdgParticle || pdgParticle->Charge() == 0.)
        continue;
      if (!particle.isPhysicalPrimary())
        continue;
      if (std::abs(particle.eta()) >= cfgCutEtaMax.value)
        continue;
      if (particle.pt() < cfgTrkLowPtCut.value)
        continue;
      hasParticleInAcceptance = true;
      break;
    }

    if (hasParticleInAcceptance) {
      inel0MCCollisions.insert(mcCollId);
    }

    int nGenCharged = countGeneratedChargedPrimaries(particlesInCollision, cfgCutEtaMax.value, cfgTrkLowPtCut.value);
    mcCollisionToNch[mcCollId] = nGenCharged;

    bool inel0 = o2::pwglf::isINELgt0mc(particlesInCollision, pdg);
    bool inel1 = o2::pwglf::isINELgt1mc(particlesInCollision, pdg);
    int inelClass = inel1 ? 2 : (inel0 ? 1 : 0);
    mcCollisionToINELClass[mcCollId] = inelClass;

    bool physicsSelected = (std::abs(mcCollision.posZ()) <= cfgCutVertex.value);
    if (cfgINELCut.value == INELgt0 && !inel0)
      physicsSelected = false;
    if (cfgINELCut.value == INELgt1 && !inel1)
      physicsSelected = false;

    if (physicsSelected) {
      physicsSelectedMCCollisions.insert(mcCollId);
      if (inel0)
        ue.fill(HIST("MC/GenRecoCollisions"), 2.f);
      if (inel1)
        ue.fill(HIST("MC/GenRecoCollisions"), 3.f);
    }

    // Fill generated particle spectra
    for (const auto& particle : particlesInCollision) {
      auto pdgParticle = pdg->GetParticle(particle.pdgCode());
      if (!pdgParticle || pdgParticle->Charge() == 0.)
        continue;
      if (!particle.isPhysicalPrimary())
        continue;

      int pdgCode = std::abs(particle.pdgCode());
      float pt = particle.pt();

      // Fill hPtPrimGenAll for ALL generated particles (NO CUTS)
      if (pdgCode == PDG_t::kPiPlus) {
        ue.fill(HIST("Pion/hPtPrimGenAll"), pt);
      } else if (pdgCode == PDG_t::kKPlus) {
        ue.fill(HIST("Kaon/hPtPrimGenAll"), pt);
      } else if (pdgCode == PDG_t::kProton) {
        ue.fill(HIST("Proton/hPtPrimGenAll"), pt);
      }

      // Fill hPtGenINEL for particles in INEL>0 events (with acceptance cuts)
      if (hasParticleInAcceptance) {
        if (std::abs(particle.eta()) >= cfgCutEtaMax.value)
          continue;
        if (particle.pt() < cfgTrkLowPtCut.value)
          continue;

        if (pdgCode == PDG_t::kPiPlus) {
          ue.fill(HIST("Pion/hPtGenINEL"), pt);
          ue.fill(HIST("Pion/hPtGenINEL_Mult"), pt, nGenCharged);
        } else if (pdgCode == PDG_t::kKPlus) {
          ue.fill(HIST("Kaon/hPtGenINEL"), pt);
          ue.fill(HIST("Kaon/hPtGenINEL_Mult"), pt, nGenCharged);
        } else if (pdgCode == PDG_t::kProton) {
          ue.fill(HIST("Proton/hPtGenINEL"), pt);
          ue.fill(HIST("Proton/hPtGenINEL_Mult"), pt, nGenCharged);
        }
      }

      // Apply acceptance cuts for physics-selected spectra
      if (std::abs(particle.eta()) >= cfgCutEtaMax.value)
        continue;
      if (particle.pt() < cfgTrkLowPtCut.value)
        continue;

      // Fill generated spectra (physics selected only)
      if (physicsSelected) {
        ue.fill(HIST("Inclusive/hPtGen"), pt);

        if (pdgCode == PDG_t::kPiPlus) {
          ue.fill(HIST("Pion/hPtGen"), pt);
        } else if (pdgCode == PDG_t::kKPlus) {
          ue.fill(HIST("Kaon/hPtGen"), pt);
        } else if (pdgCode == PDG_t::kProton) {
          ue.fill(HIST("Proton/hPtGen"), pt);
        }
      }
    }
  }

  std::map<int64_t, int64_t> recoToMcMap;
  std::map<int64_t, float> recoToCentMap;

  size_t nPairs = std::min(labels.size(), collisions.size());
  for (size_t i = 0; i < nPairs; ++i) {
    const auto& collision = collisions.iteratorAt(i);
    const auto& label = labels.iteratorAt(i);
    recoToMcMap[collision.globalIndex()] = label.mcCollisionId();
  }

  size_t nCentPairs = std::min(centTable.size(), collisions.size());
  for (size_t i = 0; i < nCentPairs; ++i) {
    const auto& collision = collisions.iteratorAt(i);
    const auto& cent = centTable.iteratorAt(i);
    float centValue = cent.centFT0M();
    recoToCentMap[collision.globalIndex()] = centValue;
    ue.fill(HIST("Centrality/hCentRaw"), centValue);
  }

  std::set<int64_t> reconstructedMCCollisions;
  std::set<int64_t> selectedMCCollisions;

  for (const auto& collision : collisions) {
    ue.fill(HIST("CutFlow/hCutStats"), 1);

    int64_t collId = collision.globalIndex();

    // MC matching
    auto mcIt = recoToMcMap.find(collId);
    if (mcIt == recoToMcMap.end())
      continue;
    ue.fill(HIST("CutFlow/hCutStats"), 2);

    int64_t mcCollId = mcIt->second;
    auto nchIt = mcCollisionToNch.find(mcCollId);
    if (nchIt == mcCollisionToNch.end())
      continue;
    int nGenCharged = nchIt->second;

    auto inelIt = mcCollisionToINELClass.find(mcCollId);
    int inelClass = (inelIt != mcCollisionToINELClass.end()) ? inelIt->second : 0;

    // Centrality
    auto centIt = recoToCentMap.find(collId);
    if (centIt == recoToCentMap.end())
      continue;
    ue.fill(HIST("CutFlow/hCutStats"), 3);

    float cent = centIt->second;
    if (cent < 0 || cent > CentBinMax)
      continue;

    // Store centrality for this MC collision (used later for MC truth plots)
    mcCollToCentFromReco[mcCollId] = cent;

    // Vertex cut
    bool passVertex = std::abs(collision.posZ()) <= cfgCutVertex.value;
    if (!passVertex)
      continue;
    ue.fill(HIST("CutFlow/hCutStats"), 4);
    ue.fill(HIST("Centrality/hCentAfterVtx"), cent);

    // INEL cut
    bool passINEL = true;
    if (cfgINELCut.value == INELgt0 && inelClass < INELgt0)
      passINEL = false;
    if (cfgINELCut.value == INELgt1 && inelClass < INELgt1)
      passINEL = false;
    if (!passINEL)
      continue;

    // Event passed all cuts
    ue.fill(HIST("CutFlow/hCutStats"), 5);
    ue.fill(HIST("Centrality/hCentAfterAll"), cent);

    ue.fill(HIST("MC/EventLoss/GenMultVsCent_Selected"), cent, nGenCharged);
    ue.fill(HIST("Centrality/hCentVsMult"), cent, nGenCharged);
    ue.fill(HIST("Centrality/hMultVsCent"), nGenCharged, cent);
    ue.fill(HIST("hvtxZ"), collision.posZ());

    selectedMCCollisions.insert(mcCollId);
    reconstructedMCCollisions.insert(mcCollId);

    // Get magnetic field for phi cut
    float magField = 0;
    if (applyPhiCut.value) {
      const auto& bc = collision.bc_as<BCsRun3>();
      magField = getMagneticField(bc.timestamp());
    }

    int nTracksInEvent = 0;

    for (const auto& track : tracks) {
      if (!track.has_collision())
        continue;
      if (track.collisionId() != collId)
        continue;

      if (!passesTrackSelection(track, magField))
        continue;

      nTracksInEvent++;

      ue.fill(HIST("hNclFoundTPC"), track.tpcNClsFound());
      ue.fill(HIST("hNclPIDTPC"), track.tpcNClsPID());
      ue.fill(HIST("hEta"), track.eta());
      ue.fill(HIST("hPhi"), track.phi());
      ue.fill(HIST("Inclusive/hPtReco"), track.pt());

      if (track.has_mcParticle()) {
        const auto& particle = track.mcParticle();
        int pdgCode = std::abs(particle.pdgCode());

        if (pdgCode == PDG_t::kPiPlus) {
          ue.fill(HIST("Pion/hPtReco"), track.pt());
          if (particle.isPhysicalPrimary()) {
            ue.fill(HIST("Pion/hPtPrimReco"), track.pt());
            ue.fill(HIST("Inclusive/hPtPrimReco"), track.pt());
          } else {
            ue.fill(HIST("Inclusive/hPtSecReco"), track.pt());
          }
        } else if (pdgCode == PDG_t::kKPlus) {
          ue.fill(HIST("Kaon/hPtReco"), track.pt());
          if (particle.isPhysicalPrimary()) {
            ue.fill(HIST("Kaon/hPtPrimReco"), track.pt());
            ue.fill(HIST("Inclusive/hPtPrimReco"), track.pt());
          } else {
            ue.fill(HIST("Inclusive/hPtSecReco"), track.pt());
          }
        } else if (pdgCode == PDG_t::kProton) {
          ue.fill(HIST("Proton/hPtReco"), track.pt());
          if (particle.isPhysicalPrimary()) {
            ue.fill(HIST("Proton/hPtPrimReco"), track.pt());
            ue.fill(HIST("Inclusive/hPtPrimReco"), track.pt());
          } else {
            ue.fill(HIST("Inclusive/hPtSecReco"), track.pt());
          }
        }
      }

      if (passesPIDSelection(track, PartPion)) {
        ue.fill(HIST("Pion/hPtMeasured"), track.pt());
        if (enablePIDHistograms) {
          ue.fill(HIST("Pion/hNsigmaTPC"), track.pt(), track.tpcNSigmaPi());
        }
      }

      if (passesPIDSelection(track, PartKaon)) {
        ue.fill(HIST("Kaon/hPtMeasured"), track.pt());
        if (enablePIDHistograms) {
          ue.fill(HIST("Kaon/hNsigmaTPC"), track.pt(), track.tpcNSigmaKa());
        }
      }

      if (passesPIDSelection(track, PartProton)) {
        ue.fill(HIST("Proton/hPtMeasured"), track.pt());
        if (enablePIDHistograms) {
          ue.fill(HIST("Proton/hNsigmaTPC"), track.pt(), track.tpcNSigmaPr());
        }
      }
    }

    ue.fill(HIST("Centrality/hRecoMultVsCent"), cent, nTracksInEvent);
  }

  for (const auto& mcCollision : mcCollisions) {
    int64_t mcCollId = mcCollision.globalIndex();

    if (reconstructedMCCollisions.find(mcCollId) == reconstructedMCCollisions.end())
      continue;

    auto centIt = mcCollToCentFromReco.find(mcCollId);
    if (centIt == mcCollToCentFromReco.end())
      continue;

    auto nchIt = mcCollisionToNch.find(mcCollId);
    int nGenCharged = (nchIt != mcCollisionToNch.end()) ? nchIt->second : 0;

    auto particlesInCollision = particles.sliceBy(perMCCol, mcCollId);

    for (const auto& particle : particlesInCollision) {
      auto pdgParticle = pdg->GetParticle(particle.pdgCode());
      if (!pdgParticle || pdgParticle->Charge() == 0.)
        continue;
      if (!particle.isPhysicalPrimary())
        continue;

      if (std::abs(particle.eta()) >= cfgCutEtaMax.value)
        continue;
      if (particle.pt() < cfgTrkLowPtCut.value)
        continue;

      int pdgCode = std::abs(particle.pdgCode());
      float pt = particle.pt();

      if (pdgCode == PDG_t::kPiPlus) {
        ue.fill(HIST("Pion/hPtGenRecoEvent"), pt);
        ue.fill(HIST("Pion/hPtGenRecoEvent_Mult"), pt, nGenCharged);
      } else if (pdgCode == PDG_t::kKPlus) {
        ue.fill(HIST("Kaon/hPtGenRecoEvent"), pt);
        ue.fill(HIST("Kaon/hPtGenRecoEvent_Mult"), pt, nGenCharged);
      } else if (pdgCode == PDG_t::kProton) {
        ue.fill(HIST("Proton/hPtGenRecoEvent"), pt);
        ue.fill(HIST("Proton/hPtGenRecoEvent_Mult"), pt, nGenCharged);
      }
    }
  }

  for (const auto& mcCollId : inel0MCCollisions) {
    auto centIt = mcCollToCentFromReco.find(mcCollId);
    if (centIt != mcCollToCentFromReco.end()) {
      float cent = centIt->second;
      int nGenCharged = mcCollisionToNch[mcCollId];
      ue.fill(HIST("MC/EventLoss/GenMultVsCent"), cent, nGenCharged);
      ue.fill(HIST("hEventsINEL_Cent"), cent);
    }
  }

  for (const auto& mcCollId : reconstructedMCCollisions) {
    auto centIt = mcCollToCentFromReco.find(mcCollId);
    if (centIt != mcCollToCentFromReco.end()) {
      ue.fill(HIST("hEventsReco_Cent"), centIt->second);
    }
  }

  ue.fill(HIST("hEventLossBreakdown"), 1.f, physicsSelectedMCCollisions.size());
  ue.fill(HIST("hEventLossBreakdown"), 2.f, reconstructedMCCollisions.size());
  ue.fill(HIST("hEventLossBreakdown"), 3.f, selectedMCCollisions.size());

  ue.fill(HIST("MC/GenRecoCollisions"), 4.f, reconstructedMCCollisions.size());
  ue.fill(HIST("MC/GenRecoCollisions"), 5.f, selectedMCCollisions.size());

  LOG(info) << "\n=== EVENT STATISTICS ===";
  LOG(info) << "INEL>0 MC collisions: " << inel0MCCollisions.size();
  LOG(info) << "Physics selected MC collisions: " << physicsSelectedMCCollisions.size();
  LOG(info) << "Reconstructed collisions: " << reconstructedMCCollisions.size();
  LOG(info) << "Selected collisions: " << selectedMCCollisions.size();

  if (physicsSelectedMCCollisions.size() > 0) {
    float efficiency = 100.f * selectedMCCollisions.size() / physicsSelectedMCCollisions.size();
    LOG(info) << "Final efficiency: " << efficiency << "%";
  }
}

void MultiplicityPt::processData(CollisionTableData::iterator const& /*collision*/,
                                 TrackTableData const& /*tracks*/,
                                 BCsRun3 const& /*bcs*/)
{
}
