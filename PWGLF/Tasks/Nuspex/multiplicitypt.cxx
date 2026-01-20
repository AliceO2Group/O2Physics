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

#include "PWGLF/Utils/inelGt.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/Logger.h" #include "ReconstructionDataFormats/Track.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/StaticFor.h"
#include "Framework/runDataProcessing.h"

#include "TDatabasePDG.h"
#include <TF1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TRandom.h>

#include <cmath>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;
using BCsRun3 = soa::Join<aod::BCs, aod::Timestamps, aod::BcSels,
                          aod::Run3MatchedToBCSparse>;

struct MultiplicityPt {

  // Service
  Service<o2::framework::O2DatabasePDG> pdg;

  Configurable<bool> isRun3{"isRun3", true, "is Run3 dataset"};
  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};
  Configurable<int> cfgINELCut{"cfgINELCut", 0, "INEL event selection: 0 no sel, 1 INEL>0, 2 INEL>1"};
  Configurable<bool> askForCustomTVX{"askForCustomTVX", false, "Ask for custom TVX rather than sel8"};
  Configurable<bool> removeITSROFrameBorder{"removeITSROFrameBorder", false, "Remove ITS Read-Out Frame border"};
  Configurable<bool> removeNoSameBunchPileup{"removeNoSameBunchPileup", false, "Remove no same bunch pileup"};
  Configurable<bool> requireIsGoodZvtxFT0vsPV{"requireIsGoodZvtxFT0vsPV", false, "Require good Z vertex FT0 vs PV"};
  Configurable<bool> requireIsVertexITSTPC{"requireIsVertexITSTPC", false, "Require vertex ITSTPC"};
  Configurable<bool> removeNoTimeFrameBorder{"removeNoTimeFrameBorder", false, "Remove no time frame border"};
  Configurable<float> cfgCutEtaMax{"cfgCutEtaMax", 0.8f, "Max eta range for tracks"};
  Configurable<float> cfgCutEtaMin{"cfgCutEtaMin", -0.8f, "Min eta range for tracks"};
  Configurable<float> cfgCutY{"cfgCutY", 0.5f, "Y range for tracks"};
  Configurable<float> cfgCutNsigma{"cfgCutNsigma", 3.0f, "nsigma cut range for tracks"};
  Configurable<int> lastRequiredTrdCluster{"lastRequiredTrdCluster", -1, "Last cluster to require in TRD"};
  Configurable<bool> requireTrdOnly{"requireTrdOnly", false, "Require only tracks from TRD"};
  Configurable<bool> requireNoTrd{"requireNoTrd", false, "Require tracks without TRD"};
  Configurable<int> multiplicityEstimator{"multiplicityEstimator", 6,
                                          "Multiplicity estimator: 0=NoMult, 1=MultFV0M, 2=MultFT0M, 3=MultFDDM, 4=MultTracklets, 5=MultTPC, 6=MultNTracksPV, 7=MultNTracksPVeta1, 8=CentFT0C, 9=CentFT0M, 10=CentFV0A"};

  // Analysis switches
  Configurable<bool> enableDCAHistograms{"enableDCAHistograms", false, "Enable DCA histograms"};
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
  Configurable<float> minTPCNClsFound{"minTPCNClsFound", 100.f, "Additional cut on the minimum value of the number of found clusters in the TPC"};
  Configurable<int> min_ITS_nClusters{"min_ITS_nClusters", 5, "minimum number of found ITS clusters"};

  // Basic track cuts
  Configurable<float> cfgTrkEtaCut{"cfgTrkEtaCut", 0.8f, "Eta range for tracks"};
  Configurable<float> cfgTrkLowPtCut{"cfgTrkLowPtCut", 0.15f, "Minimum constituent pT"};

  // Custom track cuts matching spectraTOF
  TrackSelection customTrackCuts;

  // Histogram Registry
  HistogramRegistry ue;

  // Table definitions - EXACT spectraTOF approach
  using CollisionTableData = soa::Join<aod::Collisions, aod::EvSels, aod::TPCMults, aod::PVMults, aod::MultZeqs>;
  using CollisionTableMC = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::TPCMults, aod::PVMults, aod::MultZeqs>;

  // Track tables - TPC PID only
  using TrackTableData = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection,
                                   aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr>;
  using TrackTableMC = soa::Join<aod::Tracks, aod::McTrackLabels, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection,
                                 aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr>;

  // MC tables - EXACT spectraTOF approach
  using CollisionTableMCTrue = aod::McCollisions;
  using ParticleTableMC = aod::McParticles;

  // Preslice for MC particles (like spectraTOF)
  Preslice<aod::McParticles> perMCCol = aod::mcparticle::mcCollisionId;

  // Multiplicity estimator enum (like spectraTOF)
  enum MultCodes : int {
    kNoMultiplicity = 0,
    kMultFV0M = 1,
    kMultFT0M = 2,
    kMultFDDM = 3,
    kMultTracklets = 4,
    kMultTPC = 5,
    kMultNTracksPV = 6,
    kMultNTracksPVeta1 = 7,
    kCentralityFT0C = 8,
    kCentralityFT0M = 9,
    kCentralityFV0A = 10
  };

  // Particle species enum (from spectraTOF)
  enum ParticleSpecies : int {
    kPion = 0,
    kKaon = 1,
    kProton = 2,
    kNSpecies = 3
  };

  // PDG codes
  static constexpr int PDGPion = Pdg::kPiPlus;
  static constexpr int PDGKaon = Pdg::kKPlus;
  static constexpr int PDGProton = Pdg::kProton;

  // ========================================================================
  // PROCESS FUNCTION DECLARATIONS - SPECTRATOF STYLE
  // ========================================================================

  // Data processing
  void processData(CollisionTableData::iterator const& collision,
                   TrackTableData const& tracks);
  PROCESS_SWITCH(multiplicitypt, processData, "process data", false);

  // MC processing - EXACT spectraTOF approach
  void processMC(TrackTableMC const& tracks,
                 aod::McParticles const& particles,
                 CollisionTableMCTrue const& mcCollisions,
                 CollisionTableMC const& collisions);
  PROCESS_SWITCH(multiplicitypt, processMC, "process MC", true);

  // True MC processing - EXACT spectraTOF approach
  void processTrue(CollisionTableMCTrue const& mcCollisions,
                   ParticleTableMC const& particles);
  PROCESS_SWITCH(multiplicitypt, processTrue, "process true MC", true);

  // ========================================================================
  // TRACK SELECTION FUNCTIONS - MATCHING spectraTOF
  // ========================================================================

  template <typename TrackType>
  bool passesCutWoDCA(TrackType const& track) const
  {
    if (useCustomTrackCuts.value) {
      for (int i = 0; i < static_cast<int>(TrackSelection::TrackCuts::kNCuts); i++) {
        if (i == static_cast<int>(TrackSelection::TrackCuts::kDCAxy) ||
            i == static_cast<int>(TrackSelection::TrackCuts::kDCAz)) {
          continue;
        }
        if (!customTrackCuts.IsSelected(track, static_cast<TrackSelection::TrackCuts>(i))) {
          return false;
        }
      }
      return true;
    }
    return track.isGlobalTrackWoDCA();
  }

  template <typename TrackType>
  bool passesDCAxyCut(TrackType const& track) const
  {
    if (useCustomTrackCuts.value) {
      if (!passesCutWoDCA(track)) {
        return false;
      }
      const float maxDcaXY = maxDcaXYFactor.value * (0.0105f + 0.0350f / std::pow(track.pt(), 1.1f));
      if (std::abs(track.dcaXY()) > maxDcaXY) {
        return false;
      }
      return true;
    }
    return track.isGlobalTrack();
  }

  template <typename TrackType>
  bool passesTrackSelection(TrackType const& track) const
  {
    if (track.eta() < cfgCutEtaMin.value || track.eta() > cfgCutEtaMax.value)
      return false;

    if (track.tpcChi2NCl() < minChi2PerClusterTPC.value || track.tpcChi2NCl() > maxChi2PerClusterTPC.value)
      return false;

    if (!passesCutWoDCA(track))
      return false;

    return passesDCAxyCut(track);
  }

  // ========================================================================
  // PID SELECTION FUNCTIONS - TPC ONLY (OLD NON-EXCLUSIVE METHOD)
  // ========================================================================

  template <int species, typename TrackType>
  bool passesPIDSelection(TrackType const& track) const
  {
    float nsigmaTPC = 0.f;

    if constexpr (species == kPion) {
      nsigmaTPC = track.tpcNSigmaPi();
    } else if constexpr (species == kKaon) {
      nsigmaTPC = track.tpcNSigmaKa();
    } else if constexpr (species == kProton) {
      nsigmaTPC = track.tpcNSigmaPr();
    }

    // TPC-only PID (works for all pT, but better at low pT < 1 GeV/c)
    return (std::abs(nsigmaTPC) < cfgCutNsigma.value);
  }

  // ========================================================================
  // EXCLUSIVE PID SELECTION - Returns best hypothesis for a track
  // ========================================================================

  template <typename TrackType>
  int getBestPIDHypothesis(TrackType const& track) const
  {
    // Return values: -1 = no ID, 0 = pion, 1 = kaon, 2 = proton

    float nsigmaPi = std::abs(track.tpcNSigmaPi());
    float nsigmaKa = std::abs(track.tpcNSigmaKa());
    float nsigmaPr = std::abs(track.tpcNSigmaPr());

    // Find the hypothesis with smallest |nσ| that passes the cut
    float minNSigma = 999.0f;
    int bestSpecies = -1;

    if (nsigmaPi < cfgCutNsigma.value && nsigmaPi < minNSigma) {
      minNSigma = nsigmaPi;
      bestSpecies = kPion;
    }
    if (nsigmaKa < cfgCutNsigma.value && nsigmaKa < minNSigma) {
      minNSigma = nsigmaKa;
      bestSpecies = kKaon;
    }
    if (nsigmaPr < cfgCutNsigma.value && nsigmaPr < minNSigma) {
      minNSigma = nsigmaPr;
      bestSpecies = kProton;
    }

    return bestSpecies;
  }

  // ========================================================================
  // EVENT SELECTION FUNCTION - EXACT spectraTOF
  // ========================================================================

  template <bool fillHistograms = false, typename CollisionType>
  bool isEventSelected(CollisionType const& collision)
  {
    if constexpr (fillHistograms) {
      ue.fill(HIST("evsel"), 1.f);
      if (collision.isInelGt0())
        ue.fill(HIST("evsel"), 2.f);
      if (collision.isInelGt1())
        ue.fill(HIST("evsel"), 3.f);
    }

    if (askForCustomTVX.value) {
      if (!collision.selection_bit(aod::evsel::kIsTriggerTVX))
        return false;
    } else {
      if (!collision.sel8())
        return false;
    }

    if constexpr (fillHistograms)
      ue.fill(HIST("evsel"), 4.f);

    if (removeITSROFrameBorder.value && !collision.selection_bit(aod::evsel::kNoITSROFrameBorder))
      return false;
    if constexpr (fillHistograms)
      ue.fill(HIST("evsel"), 5.f);

    if (removeNoSameBunchPileup.value && !collision.selection_bit(aod::evsel::kNoSameBunchPileup))
      return false;
    if constexpr (fillHistograms)
      ue.fill(HIST("evsel"), 6.f);

    if (requireIsGoodZvtxFT0vsPV.value && !collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV))
      return false;
    if constexpr (fillHistograms)
      ue.fill(HIST("evsel"), 7.f);

    if (requireIsVertexITSTPC.value && !collision.selection_bit(aod::evsel::kIsVertexITSTPC))
      return false;
    if constexpr (fillHistograms)
      ue.fill(HIST("evsel"), 8.f);

    if (removeNoTimeFrameBorder.value && !collision.selection_bit(aod::evsel::kNoTimeFrameBorder))
      return false;
    if constexpr (fillHistograms)
      ue.fill(HIST("evsel"), 9.f);

    if (std::abs(collision.posZ()) > cfgCutVertex.value)
      return false;

    if constexpr (fillHistograms) {
      ue.fill(HIST("evsel"), 13.f);
      if (collision.isInelGt0())
        ue.fill(HIST("evsel"), 14.f);
      if (collision.isInelGt1())
        ue.fill(HIST("evsel"), 15.f);
    }

    if (cfgINELCut.value == 1 && !collision.isInelGt0())
      return false;
    if (cfgINELCut.value == 2 && !collision.isInelGt1())
      return false;

    return true;
  }

  // ========================================================================
  // PRIMARY SELECTION - MATCHING spectraTOF
  // ========================================================================

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

    if (std::abs(particle.y()) > cfgCutY.value)
      return false;

    return true;
  }

  // Particle-specific primary selection
  template <int species, typename ParticleType>
  bool isGoodPrimarySpecies(ParticleType const& particle) const
  {
    int pdgCode = std::abs(particle.pdgCode());
    int expectedPDG = 0;

    if constexpr (species == kPion)
      expectedPDG = PDGPion;
    else if constexpr (species == kKaon)
      expectedPDG = PDGKaon;
    else if constexpr (species == kProton)
      expectedPDG = PDGProton;

    if (pdgCode != expectedPDG)
      return false;

    return isGoodPrimary(particle);
  }

  void init(InitContext const&);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<MultiplicityPt>(cfgc)};
}

void MultiplicityPt::init(InitContext const&)
{
  // ========================================================================
  // CUSTOM TRACK CUTS INITIALIZATION - MATCHING spectraTOF
  // ========================================================================

  if (useCustomTrackCuts.value) {
    LOG(info) << "Using custom track cuts matching spectraTOF approach";
    customTrackCuts = getGlobalTrackSelectionRun3ITSMatch(itsPattern.value);

    customTrackCuts.SetRequireITSRefit(requireITS.value);
    customTrackCuts.SetRequireTPCRefit(requireTPC.value);
    customTrackCuts.SetMinNClustersITS(min_ITS_nClusters.value);
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

  // ========================================================================
  // AXIS DEFINITIONS
  // ========================================================================

  ConfigurableAxis ptBinning{
    "ptBinning",
    {0.0, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45,
     0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4,
     1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8,
     3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0,
     12.0, 14.0, 16.0, 18.0, 20.0, 25.0, 30.0, 40.0, 50.0},
    "pT bin limits"};
  AxisSpec ptAxis = {ptBinning, "#it{p}_{T} (GeV/#it{c})"};

  // ========================================================================
  // HISTOGRAM REGISTRY - INCLUSIVE + PARTICLE-SPECIFIC
  // ========================================================================

  // Event counting - EXACT spectraTOF approach
  ue.add("MC/GenRecoCollisions", "Generated and Reconstructed MC Collisions", HistType::kTH1D, {{10, 0.5, 10.5}});
  auto hColl = ue.get<TH1>(HIST("MC/GenRecoCollisions"));
  hColl->GetXaxis()->SetBinLabel(1, "Collisions generated");
  hColl->GetXaxis()->SetBinLabel(2, "Collisions reconstructed");

  // CRITICAL: Complete event counting system
  ue.add("hEventsAllGen", "All generated events", HistType::kTH1F, {{1, 0.5, 1.5}});
  ue.add("hEventsPassPhysicsSelection", "Events passing physics selection", HistType::kTH1F, {{1, 0.5, 1.5}});
  ue.add("hEventsReconstructable", "Physics-selected events with reconstruction", HistType::kTH1F, {{1, 0.5, 1.5}});
  ue.add("hEventsSelectedReco", "Selected reconstructed events", HistType::kTH1F, {{1, 0.5, 1.5}});

  // Event loss breakdown histogram
  ue.add("hEventLossBreakdown", "Event loss breakdown", HistType::kTH1D, {{4, 0.5, 4.5}});
  auto hLoss = ue.get<TH1>(HIST("hEventLossBreakdown"));
  hLoss->GetXaxis()->SetBinLabel(1, "Physics selected");
  hLoss->GetXaxis()->SetBinLabel(2, "Reconstructed");
  hLoss->GetXaxis()->SetBinLabel(3, "Selected");
  hLoss->GetXaxis()->SetBinLabel(4, "Final efficiency");

  // ========================================================================
  // INCLUSIVE CHARGED PARTICLE HISTOGRAMS
  // ========================================================================

  // ALL generated primaries (before any physics selection)
  ue.add("Inclusive/hPtPrimGenAll", "All generated primaries (no cuts);#it{p}_{T} (GeV/#it{c});Counts",
         HistType::kTH1D, {ptAxis});

  // Generated primaries AFTER physics selection
  ue.add("Inclusive/hPtPrimGen", "Generated primaries (after physics selection);#it{p}_{T} (GeV/#it{c});Counts",
         HistType::kTH1D, {ptAxis});

  // Tracking Efficiency
  ue.add("Inclusive/hPtNumEff", "Tracking efficiency numerator;#it{p}_{T} (GeV/#it{c});Counts",
         HistType::kTH1D, {ptAxis});
  ue.add("Inclusive/hPtDenEff", "Tracking efficiency denominator;#it{p}_{T} (GeV/#it{c});Counts",
         HistType::kTH1D, {ptAxis});

  // Primary Fraction
  ue.add("Inclusive/hPtAllReco", "All reconstructed tracks;#it{p}_{T} (GeV/#it{c});Counts",
         HistType::kTH1D, {ptAxis});
  ue.add("Inclusive/hPtPrimReco", "Reconstructed primaries;#it{p}_{T} (GeV/#it{c});Counts",
         HistType::kTH1D, {ptAxis});
  ue.add("Inclusive/hPtSecReco", "Reconstructed secondaries;#it{p}_{T} (GeV/#it{c});Counts",
         HistType::kTH1D, {ptAxis});

  // Measured spectra
  ue.add("Inclusive/hPtMeasured", "All measured tracks;#it{p}_{T} (GeV/#it{c});Counts",
         HistType::kTH1D, {ptAxis});

  // ========================================================================
  // PARTICLE-SPECIFIC HISTOGRAMS (Pions, Kaons, Protons)
  // ========================================================================

  const std::array<std::string, kNSpecies> particleNames = {"Pion", "Kaon", "Proton"};
  const std::array<std::string, kNSpecies> particleSymbols = {"#pi^{#pm}", "K^{#pm}", "p+#bar{p}"};

  for (int iSpecies = 0; iSpecies < kNSpecies; ++iSpecies) {
    const auto& name = particleNames[iSpecies];
    const auto& symbol = particleSymbols[iSpecies];

    // Generated histograms
    ue.add(Form("%s/hPtPrimGenAll", name.c_str()),
           Form("All generated %s (no cuts);#it{p}_{T} (GeV/#it{c});Counts", symbol.c_str()),
           HistType::kTH1D, {ptAxis});

    ue.add(Form("%s/hPtPrimGen", name.c_str()),
           Form("Generated %s (after physics selection);#it{p}_{T} (GeV/#it{c});Counts", symbol.c_str()),
           HistType::kTH1D, {ptAxis});

    // Tracking efficiency
    ue.add(Form("%s/hPtNumEff", name.c_str()),
           Form("%s tracking efficiency numerator;#it{p}_{T} (GeV/#it{c});Counts", symbol.c_str()),
           HistType::kTH1D, {ptAxis});

    ue.add(Form("%s/hPtDenEff", name.c_str()),
           Form("%s tracking efficiency denominator;#it{p}_{T} (GeV/#it{c});Counts", symbol.c_str()),
           HistType::kTH1D, {ptAxis});

    // Reconstructed histograms
    ue.add(Form("%s/hPtAllReco", name.c_str()),
           Form("All reconstructed %s;#it{p}_{T} (GeV/#it{c});Counts", symbol.c_str()),
           HistType::kTH1D, {ptAxis});

    ue.add(Form("%s/hPtPrimReco", name.c_str()),
           Form("Reconstructed primary %s;#it{p}_{T} (GeV/#it{c});Counts", symbol.c_str()),
           HistType::kTH1D, {ptAxis});

    ue.add(Form("%s/hPtSecReco", name.c_str()),
           Form("Reconstructed secondary %s;#it{p}_{T} (GeV/#it{c});Counts", symbol.c_str()),
           HistType::kTH1D, {ptAxis});

    // Measured spectra
    ue.add(Form("%s/hPtMeasured", name.c_str()),
           Form("Measured %s;#it{p}_{T} (GeV/#it{c});Counts", symbol.c_str()),
           HistType::kTH1D, {ptAxis});

    // PID quality histograms - TPC ONLY
    if (enablePIDHistograms) {
      ue.add(Form("%s/hNsigmaTPC", name.c_str()),
             Form("TPC n#sigma %s;#it{p}_{T} (GeV/#it{c});n#sigma_{TPC}", symbol.c_str()),
             HistType::kTH2D, {ptAxis, {200, -10, 10}});
    }
  }

  // ========================================================================
  // MONITORING HISTOGRAMS
  // ========================================================================
  ue.add("evsel", "Event selection", HistType::kTH1D, {{20, 0.5, 20.5}});
  auto h = ue.get<TH1>(HIST("evsel"));
  h->GetXaxis()->SetBinLabel(1, "Events read");
  h->GetXaxis()->SetBinLabel(2, "INEL>0");
  h->GetXaxis()->SetBinLabel(3, "INEL>1");
  h->GetXaxis()->SetBinLabel(4, "Trigger passed");
  h->GetXaxis()->SetBinLabel(5, "NoITSROFrameBorder");
  h->GetXaxis()->SetBinLabel(6, "NoSameBunchPileup");
  h->GetXaxis()->SetBinLabel(7, "IsGoodZvtxFT0vsPV");
  h->GetXaxis()->SetBinLabel(8, "IsVertexITSTPC");
  h->GetXaxis()->SetBinLabel(9, "NoTimeFrameBorder");
  h->GetXaxis()->SetBinLabel(13, "posZ passed");
  h->GetXaxis()->SetBinLabel(14, "INEL>0 (final)");
  h->GetXaxis()->SetBinLabel(15, "INEL>1 (final)");

  ue.add("hEta", "Track eta;#eta;Counts", HistType::kTH1D, {{20, -0.8, 0.8}});
  ue.add("hPhi", "Track phi;#varphi (rad);Counts", HistType::kTH1D, {{64, 0, 2.0 * o2::constants::math::PI}});
  ue.add("hvtxZ", "Vertex Z (data);Vertex Z (cm);Events", HistType::kTH1F, {{40, -20.0, 20.0}});
  ue.add("hvtxZmc", "MC vertex Z;Vertex Z (cm);Events", HistType::kTH1F, {{40, -20.0, 20.0}});

  LOG(info) << "Initialized multiplicitypt task with EXCLUSIVE PID for INCLUSIVE + PARTICLE-SPECIFIC (Pi, K, p) analysis";
}

// ========================================================================
// DATA PROCESSING - WITH EXCLUSIVE PID
// ========================================================================
void MultiplicityPt::processData(CollisionTableData::iterator const& collision, TrackTableData const& tracks)
{
  if (!isEventSelected<true>(collision)) {
    return;
  }
  ue.fill(HIST("hvtxZ"), collision.posZ());

  for (const auto& track : tracks) {
    if (!passesTrackSelection(track)) {
      continue;
    }

    // Inclusive charged particle (always filled)
    ue.fill(HIST("Inclusive/hPtMeasured"), track.pt());
    ue.fill(HIST("hEta"), track.eta());
    ue.fill(HIST("hPhi"), track.phi());

    // Exclusive particle identification
    int bestSpecies = getBestPIDHypothesis(track);

    if (bestSpecies == kPion) {
      ue.fill(HIST("Pion/hPtMeasured"), track.pt());
      if (enablePIDHistograms) {
        ue.fill(HIST("Pion/hNsigmaTPC"), track.pt(), track.tpcNSigmaPi());
      }
    } else if (bestSpecies == kKaon) {
      ue.fill(HIST("Kaon/hPtMeasured"), track.pt());
      if (enablePIDHistograms) {
        ue.fill(HIST("Kaon/hNsigmaTPC"), track.pt(), track.tpcNSigmaKa());
      }
    } else if (bestSpecies == kProton) {
      ue.fill(HIST("Proton/hPtMeasured"), track.pt());
      if (enablePIDHistograms) {
        ue.fill(HIST("Proton/hNsigmaTPC"), track.pt(), track.tpcNSigmaPr());
      }
    }
  }
}

// ========================================================================
// MC PROCESSING - WITH FIXED PRIMARY FRACTION CALCULATION
// ========================================================================
void MultiplicityPt::processMC(TrackTableMC const& tracks,
                               aod::McParticles const& particles,
                               CollisionTableMCTrue const& mcCollisions,
                               CollisionTableMC const& collisions)
{
  LOG(info) << "=== DEBUG processMC START ===";
  LOG(info) << "MC collisions: " << mcCollisions.size();
  LOG(info) << "Reconstructed collisions: " << collisions.size();

  // ========================================================================
  // STEP 1: Identify which MC collisions are reconstructable
  // ========================================================================

  std::set<int64_t> reconstructableMCCollisions;

  for (const auto& mcCollision : mcCollisions) {
    auto particlesInCollision = particles.sliceBy(perMCCol, mcCollision.globalIndex());

    if (std::abs(mcCollision.posZ()) > cfgCutVertex.value) {
      continue;
    }
    if (cfgINELCut.value == 1 && !o2::pwglf::isINELgt0mc(particlesInCollision, pdg)) {
      continue;
    }
    if (cfgINELCut.value == 2 && !o2::pwglf::isINELgt1mc(particlesInCollision, pdg)) {
      continue;
    }

    reconstructableMCCollisions.insert(mcCollision.globalIndex());
  }

  LOG(info) << "DEBUG: Physics-selected MC collisions: " << reconstructableMCCollisions.size();

  // ========================================================================
  // STEP 2: Track reconstruction outcomes
  // ========================================================================

  std::set<int64_t> reconstructedMCCollisions;
  std::set<int64_t> selectedMCCollisions;
  std::set<int64_t> selectedCollisionIndices;

  for (const auto& collision : collisions) {
    if (!collision.has_mcCollision()) {
      continue;
    }

    const auto& mcCollision = collision.mcCollision_as<CollisionTableMCTrue>();
    int64_t mcCollId = mcCollision.globalIndex();

    if (reconstructableMCCollisions.find(mcCollId) == reconstructableMCCollisions.end()) {
      continue;
    }

    reconstructedMCCollisions.insert(mcCollId);

    if (isEventSelected<false>(collision)) {
      selectedMCCollisions.insert(mcCollId);
      selectedCollisionIndices.insert(collision.globalIndex());
      ue.fill(HIST("hvtxZ"), collision.posZ());
    }
  }

  auto hEventsReconstructable = ue.get<TH1>(HIST("hEventsReconstructable"));
  auto hEventsSelectedReco = ue.get<TH1>(HIST("hEventsSelectedReco"));

  hEventsReconstructable->SetBinContent(1, reconstructedMCCollisions.size());
  hEventsSelectedReco->SetBinContent(1, selectedMCCollisions.size());

  int nReconstructableTotal = reconstructableMCCollisions.size();
  int nReconstructableWithReco = reconstructedMCCollisions.size();
  int nSelectedReco = selectedMCCollisions.size();

  LOG(info) << "DEBUG: Reconstructed MC collisions: " << nReconstructableWithReco;
  LOG(info) << "DEBUG: Selected MC collisions: " << nSelectedReco;

  if (nReconstructableTotal > 0) {
    ue.fill(HIST("hEventLossBreakdown"), 1, nReconstructableTotal);
    ue.fill(HIST("hEventLossBreakdown"), 2, nReconstructableWithReco);
    ue.fill(HIST("hEventLossBreakdown"), 3, nSelectedReco);
    ue.fill(HIST("hEventLossBreakdown"), 4, (nSelectedReco * 100.0 / nReconstructableTotal));
  }

  // ========================================================================
  // STEP 3: Process tracks with EXCLUSIVE PID - FIXED PRIMARY FRACTION
  // ========================================================================

  int totalTracksProcessed = 0;
  int tracksFromSelectedEvents = 0;
  int tracksPassingSelection = 0;

  std::array<int, kNSpecies> particleTracksIdentified = {0};
  std::array<int, kNSpecies> particleTracksPrimary = {0};
  std::array<int, kNSpecies> particleTracksSecondary = {0};

  for (const auto& track : tracks) {
    totalTracksProcessed++;

    if (!track.has_collision())
      continue;

    const auto& collision = track.collision_as<CollisionTableMC>();

    if (selectedCollisionIndices.find(collision.globalIndex()) == selectedCollisionIndices.end()) {
      continue;
    }
    tracksFromSelectedEvents++;

    if (!passesTrackSelection(track))
      continue;
    tracksPassingSelection++;

    // ========================================================================
    // INCLUSIVE CHARGED PARTICLE ANALYSIS
    // ========================================================================

    ue.fill(HIST("Inclusive/hPtMeasured"), track.pt());
    ue.fill(HIST("Inclusive/hPtAllReco"), track.pt());
    ue.fill(HIST("hEta"), track.eta());
    ue.fill(HIST("hPhi"), track.phi());

    // ========================================================================
    // EFFICIENCY NUMERATOR: Fill based on TRUE particle type
    // ========================================================================

    if (track.has_mcParticle()) {
      const auto& particle = track.mcParticle();
      int pdgCode = std::abs(particle.pdgCode());

      if (particle.isPhysicalPrimary()) {
        ue.fill(HIST("Inclusive/hPtNumEff"), particle.pt());
        ue.fill(HIST("Inclusive/hPtPrimReco"), track.pt());

        // Fill particle-specific efficiency numerator based on TRUE type
        if (pdgCode == PDGPion) {
          ue.fill(HIST("Pion/hPtNumEff"), particle.pt());
        }
        if (pdgCode == PDGKaon) {
          ue.fill(HIST("Kaon/hPtNumEff"), particle.pt());
        }
        if (pdgCode == PDGProton) {
          ue.fill(HIST("Proton/hPtNumEff"), particle.pt());
        }
      } else {
        ue.fill(HIST("Inclusive/hPtSecReco"), track.pt());
      }
    }

    // ========================================================================
    // EXCLUSIVE PID - FIXED PRIMARY FRACTION LOGIC
    // ========================================================================

    int bestSpecies = getBestPIDHypothesis(track);

    // ========================================================================
    // PION CHANNEL
    // ========================================================================
    if (bestSpecies == kPion) {
      ue.fill(HIST("Pion/hPtMeasured"), track.pt());
      ue.fill(HIST("Pion/hPtAllReco"), track.pt());
      particleTracksIdentified[kPion]++;

      if (enablePIDHistograms) {
        ue.fill(HIST("Pion/hNsigmaTPC"), track.pt(), track.tpcNSigmaPi());
      }

      if (track.has_mcParticle()) {
        const auto& particle = track.mcParticle();

        // KEY FIX: Primary fraction = fraction of identified pions that are primary
        // This includes correctly identified pions AND misidentified kaons/protons
        // that happen to be primary particles
        if (particle.isPhysicalPrimary()) {
          ue.fill(HIST("Pion/hPtPrimReco"), track.pt());
          particleTracksPrimary[kPion]++;
        } else {
          ue.fill(HIST("Pion/hPtSecReco"), track.pt());
          particleTracksSecondary[kPion]++;
        }
      }
    }

    // ========================================================================
    // KAON CHANNEL
    // ========================================================================
    else if (bestSpecies == kKaon) {
      ue.fill(HIST("Kaon/hPtMeasured"), track.pt());
      ue.fill(HIST("Kaon/hPtAllReco"), track.pt());
      particleTracksIdentified[kKaon]++;

      if (enablePIDHistograms) {
        ue.fill(HIST("Kaon/hNsigmaTPC"), track.pt(), track.tpcNSigmaKa());
      }

      if (track.has_mcParticle()) {
        const auto& particle = track.mcParticle();

        // KEY FIX: Primary fraction of identified kaons
        // A misidentified pion that is primary still contributes to primary fraction
        if (particle.isPhysicalPrimary()) {
          ue.fill(HIST("Kaon/hPtPrimReco"), track.pt());
          particleTracksPrimary[kKaon]++;
        } else {
          ue.fill(HIST("Kaon/hPtSecReco"), track.pt());
          particleTracksSecondary[kKaon]++;
        }
      }
    }

    // ========================================================================
    // PROTON CHANNEL
    // ========================================================================
    else if (bestSpecies == kProton) {
      ue.fill(HIST("Proton/hPtMeasured"), track.pt());
      ue.fill(HIST("Proton/hPtAllReco"), track.pt());
      particleTracksIdentified[kProton]++;

      if (enablePIDHistograms) {
        ue.fill(HIST("Proton/hNsigmaTPC"), track.pt(), track.tpcNSigmaPr());
      }

      if (track.has_mcParticle()) {
        const auto& particle = track.mcParticle();

        // KEY FIX: Primary fraction of identified protons
        if (particle.isPhysicalPrimary()) {
          ue.fill(HIST("Proton/hPtPrimReco"), track.pt());
          particleTracksPrimary[kProton]++;
        } else {
          ue.fill(HIST("Proton/hPtSecReco"), track.pt());
          particleTracksSecondary[kProton]++;
        }
      }
    }
  }

  LOG(info) << "=== DEBUG TRACK COUNTING ===";
  LOG(info) << "Total tracks processed: " << totalTracksProcessed;
  LOG(info) << "Tracks from selected events: " << tracksFromSelectedEvents;
  LOG(info) << "Tracks passing selection: " << tracksPassingSelection;
  LOG(info) << "Pions identified: " << particleTracksIdentified[kPion]
            << ", primary: " << particleTracksPrimary[kPion]
            << ", secondary: " << particleTracksSecondary[kPion];
  LOG(info) << "Kaons identified: " << particleTracksIdentified[kKaon]
            << ", primary: " << particleTracksPrimary[kKaon]
            << ", secondary: " << particleTracksSecondary[kKaon];
  LOG(info) << "Protons identified: " << particleTracksIdentified[kProton]
            << ", primary: " << particleTracksPrimary[kProton]
            << ", secondary: " << particleTracksSecondary[kProton];

  // Calculate and log primary fractions
  if (particleTracksIdentified[kPion] > 0) {
    float pionPrimFrac = (float)particleTracksPrimary[kPion] / particleTracksIdentified[kPion];
    LOG(info) << "Pion primary fraction: " << pionPrimFrac * 100.0 << "%";
  }
  if (particleTracksIdentified[kKaon] > 0) {
    float kaonPrimFrac = (float)particleTracksPrimary[kKaon] / particleTracksIdentified[kKaon];
    LOG(info) << "Kaon primary fraction: " << kaonPrimFrac * 100.0 << "%";
  }
  if (particleTracksIdentified[kProton] > 0) {
    float protonPrimFrac = (float)particleTracksPrimary[kProton] / particleTracksIdentified[kProton];
    LOG(info) << "Proton primary fraction: " << protonPrimFrac * 100.0 << "%";
  }

  LOG(info) << "=== DEBUG processMC END ===";
}

// ========================================================================
// TRUE MC PROCESSING - WITH PARTICLE-SPECIFIC SIGNAL LOSS
// ========================================================================
void MultiplicityPt::processTrue(CollisionTableMCTrue const& mcCollisions,
                                 ParticleTableMC const& particles)
{
  LOG(info) << "=== DEBUG processTrue START ===";
  LOG(info) << "Number of MC collisions: " << mcCollisions.size();

  int nPassPhysicsSelection = 0;
  int nParticlesFilledAll = 0;
  int nParticlesFilledAfterPS = 0;

  std::array<int, kNSpecies> particleCountAll = {0};
  std::array<int, kNSpecies> particleCountAfterPS = {0};

  for (const auto& mcCollision : mcCollisions) {
    // Count EVERY generated event
    ue.fill(HIST("hEventsAllGen"), 1.0);

    ue.fill(HIST("hvtxZmc"), mcCollision.posZ());
    auto particlesInCollision = particles.sliceBy(perMCCol, mcCollision.globalIndex());

    // ========================================================================
    // Fill ALL generated primaries BEFORE physics selection
    // ========================================================================
    for (const auto& particle : particlesInCollision) {
      if (isGoodPrimary(particle)) {
        ue.fill(HIST("Inclusive/hPtPrimGenAll"), particle.pt());
        nParticlesFilledAll++;
      }

      if (isGoodPrimarySpecies<kPion>(particle)) {
        ue.fill(HIST("Pion/hPtPrimGenAll"), particle.pt());
        particleCountAll[kPion]++;
      }

      if (isGoodPrimarySpecies<kKaon>(particle)) {
        ue.fill(HIST("Kaon/hPtPrimGenAll"), particle.pt());
        particleCountAll[kKaon]++;
      }

      if (isGoodPrimarySpecies<kProton>(particle)) {
        ue.fill(HIST("Proton/hPtPrimGenAll"), particle.pt());
        particleCountAll[kProton]++;
      }
    }

    // ========================================================================
    // Apply physics selection
    // ========================================================================
    if (std::abs(mcCollision.posZ()) > cfgCutVertex.value)
      continue;

    if (cfgINELCut.value == 1 && !o2::pwglf::isINELgt0mc(particlesInCollision, pdg))
      continue;
    if (cfgINELCut.value == 2 && !o2::pwglf::isINELgt1mc(particlesInCollision, pdg))
      continue;

    // Count physics-selected events
    ue.fill(HIST("hEventsPassPhysicsSelection"), 1.0);
    nPassPhysicsSelection++;

    // Fill primaries AFTER physics selection
    for (const auto& particle : particlesInCollision) {
      if (isGoodPrimary(particle)) {
        ue.fill(HIST("Inclusive/hPtDenEff"), particle.pt());
        ue.fill(HIST("Inclusive/hPtPrimGen"), particle.pt());
        nParticlesFilledAfterPS++;
      }

      if (isGoodPrimarySpecies<kPion>(particle)) {
        ue.fill(HIST("Pion/hPtDenEff"), particle.pt());
        ue.fill(HIST("Pion/hPtPrimGen"), particle.pt());
        particleCountAfterPS[kPion]++;
      }

      if (isGoodPrimarySpecies<kKaon>(particle)) {
        ue.fill(HIST("Kaon/hPtDenEff"), particle.pt());
        ue.fill(HIST("Kaon/hPtPrimGen"), particle.pt());
        particleCountAfterPS[kKaon]++;
      }

      if (isGoodPrimarySpecies<kProton>(particle)) {
        ue.fill(HIST("Proton/hPtDenEff"), particle.pt());
        ue.fill(HIST("Proton/hPtPrimGen"), particle.pt());
        particleCountAfterPS[kProton]++;
      }
    }
  }

  LOG(info) << "=== DEBUG processTrue END ===";
  LOG(info) << "All generated events: " << mcCollisions.size();
  LOG(info) << "Passing physics selection: " << nPassPhysicsSelection;
  LOG(info) << "Total primaries (before PS): " << nParticlesFilledAll;
  LOG(info) << "Total primaries (after PS): " << nParticlesFilledAfterPS;

  LOG(info) << "=== PARTICLE-SPECIFIC STATISTICS ===";
  LOG(info) << "Pions - All: " << particleCountAll[kPion]
            << ", After PS: " << particleCountAfterPS[kPion];
  LOG(info) << "Kaons - All: " << particleCountAll[kKaon]
            << ", After PS: " << particleCountAfterPS[kKaon];
  LOG(info) << "Protons - All: " << particleCountAll[kProton]
            << ", After PS: " << particleCountAfterPS[kProton];
}
