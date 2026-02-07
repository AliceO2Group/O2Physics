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

#include "PWGLF/DataModel/LFParticleIdentification.h"
#include "PWGLF/DataModel/spectraTOF.h"
#include "PWGLF/Utils/inelGt.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/McCollisionExtra.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/Logger.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/StaticFor.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

#include "TPDGCode.h"
#include <TF1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TRandom.h>

#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using BCsRun3 = soa::Join<aod::BCs, aod::Timestamps, aod::BcSels,
                          aod::Run3MatchedToBCSparse>;

struct MultiplicityPt {

  // Service
  Service<o2::framework::O2DatabasePDG> pdg;

  // Add CCDB service for magnetic field
  Service<ccdb::BasicCCDBManager> ccdb;

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
  Configurable<float> minTPCNClsFound{"minTPCNClsFound", 70.f, "Additional cut on the minimum value of the number of found clusters in the TPC"};
  Configurable<int> min_ITS_nClusters{"min_ITS_nClusters", 5, "minimum number of found ITS clusters"};

  // Phi cut parameters
  Configurable<bool> applyPhiCut{"applyPhiCut", true, "Apply phi sector cut to remove problematic TPC regions"};
  Configurable<float> pTthresholdPhiCut{"pTthresholdPhiCut", 2.0f, "pT threshold above which to apply phi cut"};
  Configurable<double> phiCutLowParam1{"phiCutLowParam1", 0.119297, "First parameter for low phi cut"};
  Configurable<double> phiCutLowParam2{"phiCutLowParam2", 0.000379693, "Second parameter for low phi cut"};
  Configurable<double> phiCutHighParam1{"phiCutHighParam1", 0.16685, "First parameter for high phi cut"};
  Configurable<double> phiCutHighParam2{"phiCutHighParam2", 0.00981942, "Second parameter for high phi cut"};

  // Basic track cuts
  Configurable<float> cfgTrkEtaCut{"cfgTrkEtaCut", 0.8f, "Eta range for tracks"};
  Configurable<float> cfgTrkLowPtCut{"cfgTrkLowPtCut", 0.15f, "Minimum constituent pT"};

  // Custom track cuts matching spectraTOF
  TrackSelection customTrackCuts;

  // TF1 pointers for phi cuts
  TF1* fphiCutLow = nullptr;
  TF1* fphiCutHigh = nullptr;

  // Histogram Registry
  HistogramRegistry ue;

  // ========================================================================
  // CENTRALITY/MULTIPLICITY CLASSES - Using same bins as before for consistency
  // ========================================================================
  static constexpr int kCentralityClasses = 10;
  static constexpr double CentClasses[kCentralityClasses + 1] = {0.0, 1.0, 5.0, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 70.0, 100.0};

  // Multiplicity percentile boundaries (computed on first pass)
  std::vector<float> multPercentileboundaries;
  bool percentilesComputed = false;

  // Storage for multiplicity distribution (for percentile calculation)
  std::vector<float> multiplicityValues;

  // Table definitions - NO McCentFT0Ms dependency
  using CollisionTableData = soa::Join<aod::Collisions, aod::EvSels, aod::TPCMults, aod::PVMults, aod::MultZeqs>;
  using CollisionTableMC = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::TPCMults, aod::PVMults, aod::MultZeqs>;

  // Track tables - TPC PID only
  using TrackTableData = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection,
                                   aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr>;
  using TrackTableMC = soa::Join<aod::Tracks, aod::McTrackLabels, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection,
                                 aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr>;

  // MC tables - NO McCentFT0Ms
  using CollisionTableMCTrue = aod::McCollisions;
  using ParticleTableMC = aod::McParticles;

  // Preslice for MC particles
  Preslice<aod::McParticles> perMCCol = aod::mcparticle::mcCollisionId;

  // Multiplicity estimator enum
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

  // Particle species enum
  enum ParticleSpecies : int {
    kPion = 0,
    kKaon = 1,
    kProton = 2,
    kNSpecies = 3
  };

  // PDG codes
  static constexpr int PDGPion = 211;
  static constexpr int PDGKaon = 321;
  static constexpr int PDGProton = 2212;

  void processData(CollisionTableData::iterator const& collision,
                   TrackTableData const& tracks,
                   BCsRun3 const& bcs);
  PROCESS_SWITCH(MultiplicityPt, processData, "process data", false);

  // MC processing - First pass to build percentiles
  void processPercentileCalibration(CollisionTableMCTrue const& mcCollisions,
                                    ParticleTableMC const& particles);
  PROCESS_SWITCH(MultiplicityPt, processPercentileCalibration, "Build multiplicity percentile calibration (run first)", false);

  // MC processing - Main analysis
  void processMC(TrackTableMC const& tracks,
                 aod::McParticles const& particles,
                 CollisionTableMCTrue const& mcCollisions,
                 CollisionTableMC const& collisions,
                 BCsRun3 const& bcs);
  PROCESS_SWITCH(MultiplicityPt, processMC, "process MC", true);

  // True MC processing
  void processTrue(CollisionTableMCTrue const& mcCollisions,
                   ParticleTableMC const& particles);
  PROCESS_SWITCH(MultiplicityPt, processTrue, "process true MC", true);

  // ========================================================================
  // MULTIPLICITY GETTER FUNCTIONS - Using raw charged particle count
  // ========================================================================

  // Count charged primaries in |eta| < 1.0
  template <typename MCCollisionType>
  int countChargedPrimaries(const MCCollisionType& mcCollision, const ParticleTableMC& particles) const
  {
    int nCharged = 0;
    auto particlesInColl = particles.sliceBy(perMCCol, mcCollision.globalIndex());
    for (const auto& p : particlesInColl) {
      if (!p.isPhysicalPrimary())
        continue;
      auto pdgParticle = pdg->GetParticle(p.pdgCode());
      if (!pdgParticle || pdgParticle->Charge() == 0.)
        continue;
      if (std::abs(p.eta()) < 1.0)
        nCharged++;
    }
    return nCharged;
  }

  // For reconstructed collisions
  template <typename CollisionType>
  float getMultiplicity(const CollisionType& collision) const
  {
    switch (multiplicityEstimator.value) {
      case kNoMultiplicity:
        return 50.f;
      case kMultFV0M:
        return collision.multZeqFV0A();
      case kMultFT0M:
        return collision.multZeqFT0A() + collision.multZeqFT0C();
      case kMultFDDM:
        return collision.multZeqFDDA() + collision.multZeqFDDC();
      case kMultTracklets:
        return 0.f;
      case kMultTPC:
        return collision.multTPC();
      case kMultNTracksPV:
        return collision.multZeqNTracksPV();
      case kMultNTracksPVeta1:
        return collision.multNTracksPVeta1();
      case kCentralityFT0C:
      case kCentralityFT0M:
      case kCentralityFV0A:
        return collision.multZeqNTracksPV();
      default:
        return 0.f;
    }
  }

  // For MC collisions - returns RAW multiplicity
  template <typename MCCollisionType>
  float getMultiplicityMC(const MCCollisionType& mcCollision, const ParticleTableMC& particles) const
  {
    return static_cast<float>(countChargedPrimaries(mcCollision, particles));
  }

  // Convert raw multiplicity to percentile
  float multiplicityToPercentile(float rawMult) const
  {
    if (!percentilesComputed || multPercentileboundaries.empty()) {
      // If percentiles not computed, return raw multiplicity
      return rawMult;
    }

    // Find which percentile bin this multiplicity falls into
    for (size_t i = 0; i < multPercentileboundaries.size() - 1; ++i) {
      if (rawMult >= multPercentileboundaries[i] && rawMult < multPercentileboundaries[i + 1]) {
        // Return the CENTER of the percentile bin
        return CentClasses[i] + (CentClasses[i + 1] - CentClasses[i]) / 2.0;
      }
    }

    // Handle edge cases
    if (rawMult < multPercentileboundaries[0]) {
      return CentClasses[0];
    }
    return CentClasses[kCentralityClasses];
  }

  // Get centrality class index from raw multiplicity
  int getCentralityClass(float rawMult) const
  {
    if (!percentilesComputed || multPercentileboundaries.empty()) {
      // Fallback: divide into equal bins
      float maxMult = 150.0f; // Assumed maximum
      int bin = static_cast<int>((rawMult / maxMult) * kCentralityClasses);
      return std::min(bin, kCentralityClasses - 1);
    }

    // Use computed percentiles
    for (int i = 0; i < kCentralityClasses; ++i) {
      if (rawMult >= multPercentileboundaries[i] && rawMult < multPercentileboundaries[i + 1]) {
        return i;
      }
    }

    // Outside range
    if (rawMult < multPercentileboundaries[0])
      return 0;
    return kCentralityClasses - 1;
  }

  // ========================================================================
  // COMPUTE PERCENTILE BOUNDARIES
  // ========================================================================
  void computePercentileBoundaries()
  {
    if (multiplicityValues.empty()) {
      LOG(warning) << "No multiplicity values to compute percentiles from!";
      return;
    }

    // Sort multiplicity values
    std::sort(multiplicityValues.begin(), multiplicityValues.end());

    LOG(info) << "Computing percentile boundaries from " << multiplicityValues.size() << " events";

    // Compute percentile boundaries
    multPercentileboundaries.clear();
    multPercentileboundaries.reserve(kCentralityClasses + 1);

    for (int i = 0; i <= kCentralityClasses; ++i) {
      float percentile = CentClasses[i];
      size_t index = static_cast<size_t>(percentile / 100.0 * multiplicityValues.size());
      if (index >= multiplicityValues.size()) {
        index = multiplicityValues.size() - 1;
      }
      float boundary = multiplicityValues[index];
      multPercentileboundaries.push_back(boundary);
      LOG(info) << "Percentile " << percentile << "% -> Multiplicity >= " << boundary;
    }

    percentilesComputed = true;

    LOG(info) << "=== Percentile Boundaries Computed ===";
    for (int i = 0; i < kCentralityClasses; ++i) {
      LOG(info) << "Class " << i << ": [" << CentClasses[i] << "%-" << CentClasses[i + 1]
                << "%] = Mult [" << multPercentileboundaries[i] << "-" << multPercentileboundaries[i + 1] << ")";
    }
  }

  // ========================================================================
  // MAGNETIC FIELD FUNCTION
  // ========================================================================
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

  // ========================================================================
  // PHI CUT FUNCTION
  // ========================================================================
  template <typename TrackType>
  bool passedPhiCut(const TrackType& track, float magField) const
  {
    if (!applyPhiCut.value) {
      return true;
    }

    if (track.pt() < pTthresholdPhiCut.value) {
      return true;
    }

    float pt = track.pt();
    float phi = track.phi();
    int charge = track.sign();

    if (magField < 0) {
      phi = o2::constants::math::TwoPI - phi;
    }
    if (charge < 0) {
      phi = o2::constants::math::TwoPI - phi;
    }

    phi += o2::constants::math::PI / 18.0f;
    phi = std::fmod(phi, o2::constants::math::PI / 9.0f);

    if (phi < fphiCutHigh->Eval(pt) && phi > fphiCutLow->Eval(pt)) {
      return false;
    }

    return true;
  }

  float getTransformedPhi(const float phi, const int charge, const float magField) const
  {
    float transformedPhi = phi;
    if (magField < 0) {
      transformedPhi = o2::constants::math::TwoPI - transformedPhi;
    }
    if (charge < 0) {
      transformedPhi = o2::constants::math::TwoPI - transformedPhi;
    }
    transformedPhi += o2::constants::math::PI / 18.0f;
    transformedPhi = std::fmod(transformedPhi, o2::constants::math::PI / 9.0f);
    return transformedPhi;
  }

  // ========================================================================
  // TRACK SELECTION FUNCTIONS
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
      constexpr float dcaXYConst = 0.0105f;
      constexpr float dcaXYPtScale = 0.0350f;
      constexpr float dcaXYPtPower = 1.1f;
      const float maxDcaXY = maxDcaXYFactor.value * (dcaXYConst + dcaXYPtScale / std::pow(track.pt(), dcaXYPtPower));
      if (std::abs(track.dcaXY()) > maxDcaXY) {
        return false;
      }
      return true;
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

    if (applyPhiCut.value && !passedPhiCut(track, magField))
      return false;

    return passesDCAxyCut(track);
  }

  // ========================================================================
  // PID SELECTION FUNCTIONS
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

    return (std::abs(nsigmaTPC) < cfgCutNsigma.value);
  }

  template <typename TrackType>
  int getBestPIDHypothesis(TrackType const& track) const
  {
    float nsigmaPi = std::abs(track.tpcNSigmaPi());
    float nsigmaKa = std::abs(track.tpcNSigmaKa());
    float nsigmaPr = std::abs(track.tpcNSigmaPr());

    constexpr float largeNSigmaValue = 999.0f;
    float minNSigma = largeNSigmaValue;
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
  // EVENT SELECTION FUNCTION
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
  // PRIMARY SELECTION
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
  // CUSTOM TRACK CUTS INITIALIZATION
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
  // PHI CUT INITIALIZATION
  // ========================================================================

  if (applyPhiCut.value) {
    fphiCutLow = new TF1("StandardPhiCutLow",
                         Form("%f/x/x+pi/18.0-%f",
                              phiCutLowParam1.value, phiCutLowParam2.value),
                         0, 50);
    fphiCutHigh = new TF1("StandardPhiCutHigh",
                          Form("%f/x+pi/18.0+%f",
                               phiCutHighParam1.value, phiCutHighParam2.value),
                          0, 50);

    LOGF(info, "=== Phi Cut Parameters ===");
    LOGF(info, "Low cut: %.6f/x² + pi/18 - %.6f",
         phiCutLowParam1.value, phiCutLowParam2.value);
    LOGF(info, "High cut: %.6f/x + pi/18 + %.6f",
         phiCutHighParam1.value, phiCutHighParam2.value);
    LOGF(info, "Applied for pT > %.1f GeV/c", pTthresholdPhiCut.value);
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

  // Multiplicity axis - initially raw multiplicity, will represent percentiles after calibration
  std::vector<double> centBins(CentClasses, CentClasses + kCentralityClasses + 1);
  AxisSpec multAxis = {centBins, "Centrality/Multiplicity Class (%)"};

  // Raw multiplicity axis for calibration
  AxisSpec rawMultAxis = {150, 0, 150, "N_{ch} (|#eta| < 1.0)"};

  // ========================================================================
  // HISTOGRAM REGISTRY
  // ========================================================================

  // Multiplicity distribution for percentile calibration
  ue.add("Calibration/hRawMultiplicity", "Raw multiplicity distribution;N_{ch};Events",
         HistType::kTH1D, {rawMultAxis});

  // Event counting
  ue.add("MC/GenRecoCollisions", "Generated and Reconstructed MC Collisions", HistType::kTH1D, {{10, 0.5, 10.5}});
  auto hColl = ue.get<TH1>(HIST("MC/GenRecoCollisions"));
  hColl->GetXaxis()->SetBinLabel(1, "Collisions generated");
  hColl->GetXaxis()->SetBinLabel(2, "Collisions reconstructed");

  // Event loss histograms
  ue.add("MC/EventLoss/MultGenerated", "Generated events vs multiplicity",
         HistType::kTH1D, {multAxis});
  ue.add("MC/EventLoss/MultBadVertex", "Events with bad vertex vs multiplicity",
         HistType::kTH1D, {multAxis});
  ue.add("MC/EventLoss/MultPhysicsSelected", "Physics-selected events vs multiplicity",
         HistType::kTH1D, {multAxis});
  ue.add("MC/EventLoss/MultReconstructed", "Reconstructed events vs multiplicity",
         HistType::kTH1D, {multAxis});
  ue.add("MC/EventLoss/MultRecoSelected", "Reconstructed+selected events vs multiplicity",
         HistType::kTH1D, {multAxis});

  ue.add("hEventLossBreakdown", "Event loss breakdown", HistType::kTH1D, {{4, 0.5, 4.5}});
  auto hLoss = ue.get<TH1>(HIST("hEventLossBreakdown"));
  hLoss->GetXaxis()->SetBinLabel(1, "Physics selected");
  hLoss->GetXaxis()->SetBinLabel(2, "Reconstructed");
  hLoss->GetXaxis()->SetBinLabel(3, "Selected");
  hLoss->GetXaxis()->SetBinLabel(4, "Final efficiency");

  // ========================================================================
  // INCLUSIVE CHARGED PARTICLE HISTOGRAMS
  // ========================================================================

  ue.add("Inclusive/hPtPrimGenAll", "All generated primaries (no cuts);#it{p}_{T} (GeV/#it{c});Counts",
         HistType::kTH1D, {ptAxis});
  ue.add("Inclusive/hPtPrimGenAllVsMult", "All generated primaries vs mult;#it{p}_{T} (GeV/#it{c});Mult Class (%)",
         HistType::kTH2D, {ptAxis, multAxis});

  ue.add("Inclusive/hPtPrimBadVertex", "Generated primaries (bad vertex);#it{p}_{T} (GeV/#it{c});Counts",
         HistType::kTH1D, {ptAxis});
  ue.add("Inclusive/hPtPrimBadVertexVsMult", "Generated primaries (bad vertex) vs mult;#it{p}_{T} (GeV/#it{c});Mult Class (%)",
         HistType::kTH2D, {ptAxis, multAxis});

  ue.add("Inclusive/hPtPrimGen", "Generated primaries (after physics selection);#it{p}_{T} (GeV/#it{c});Counts",
         HistType::kTH1D, {ptAxis});
  ue.add("Inclusive/hPtPrimGenVsMult", "Generated primaries (after phys sel) vs mult;#it{p}_{T} (GeV/#it{c});Mult Class (%)",
         HistType::kTH2D, {ptAxis, multAxis});

  ue.add("Inclusive/hPtPrimRecoEv", "Generated primaries (reco events);#it{p}_{T} (GeV/#it{c});Counts",
         HistType::kTH1D, {ptAxis});
  ue.add("Inclusive/hPtPrimRecoEvVsMult", "Generated primaries (reco events) vs mult;#it{p}_{T} (GeV/#it{c});Mult Class (%)",
         HistType::kTH2D, {ptAxis, multAxis});

  ue.add("Inclusive/hPtPrimGoodEv", "Generated primaries (good events);#it{p}_{T} (GeV/#it{c});Counts",
         HistType::kTH1D, {ptAxis});
  ue.add("Inclusive/hPtPrimGoodEvVsMult", "Generated primaries (good events) vs mult;#it{p}_{T} (GeV/#it{c});Mult Class (%)",
         HistType::kTH2D, {ptAxis, multAxis});

  ue.add("Inclusive/hPtNumEff", "Tracking efficiency numerator;#it{p}_{T} (GeV/#it{c});Counts",
         HistType::kTH1D, {ptAxis});
  ue.add("Inclusive/hPtNumEffVsMult", "Tracking efficiency numerator vs mult;#it{p}_{T} (GeV/#it{c});Mult Class (%)",
         HistType::kTH2D, {ptAxis, multAxis});

  ue.add("Inclusive/hPtDenEff", "Tracking efficiency denominator;#it{p}_{T} (GeV/#it{c});Counts",
         HistType::kTH1D, {ptAxis});
  ue.add("Inclusive/hPtDenEffVsMult", "Tracking efficiency denominator vs mult;#it{p}_{T} (GeV/#it{c});Mult Class (%)",
         HistType::kTH2D, {ptAxis, multAxis});

  ue.add("Inclusive/hPtAllReco", "All reconstructed tracks;#it{p}_{T} (GeV/#it{c});Counts",
         HistType::kTH1D, {ptAxis});
  ue.add("Inclusive/hPtAllRecoVsMult", "All reconstructed tracks vs mult;#it{p}_{T} (GeV/#it{c});Mult Class (%)",
         HistType::kTH2D, {ptAxis, multAxis});

  ue.add("Inclusive/hPtPrimReco", "Reconstructed primaries;#it{p}_{T} (GeV/#it{c});Counts",
         HistType::kTH1D, {ptAxis});
  ue.add("Inclusive/hPtPrimRecoVsMult", "Reconstructed primaries vs mult;#it{p}_{T} (GeV/#it{c});Mult Class (%)",
         HistType::kTH2D, {ptAxis, multAxis});

  ue.add("Inclusive/hPtSecReco", "Reconstructed secondaries;#it{p}_{T} (GeV/#it{c});Counts",
         HistType::kTH1D, {ptAxis});
  ue.add("Inclusive/hPtSecRecoVsMult", "Reconstructed secondaries vs mult;#it{p}_{T} (GeV/#it{c});Mult Class (%)",
         HistType::kTH2D, {ptAxis, multAxis});

  ue.add("Inclusive/hPtMeasured", "All measured tracks;#it{p}_{T} (GeV/#it{c});Counts",
         HistType::kTH1D, {ptAxis});
  ue.add("Inclusive/hPtMeasuredVsMult", "All measured tracks vs mult;#it{p}_{T} (GeV/#it{c});Mult Class (%)",
         HistType::kTH2D, {ptAxis, multAxis});

  // ========================================================================
  // PARTICLE-SPECIFIC HISTOGRAMS
  // ========================================================================

  const std::array<std::string, kNSpecies> particleNames = {"Pion", "Kaon", "Proton"};
  const std::array<std::string, kNSpecies> particleSymbols = {"#pi^{#pm}", "K^{#pm}", "p+#bar{p}"};

  for (int iSpecies = 0; iSpecies < kNSpecies; ++iSpecies) {
    const auto& name = particleNames[iSpecies];
    const auto& symbol = particleSymbols[iSpecies];

    // 1D versions
    ue.add(Form("%s/hPtPrimGenAll", name.c_str()),
           Form("All generated %s (no cuts);#it{p}_{T} (GeV/#it{c});Counts", symbol.c_str()),
           HistType::kTH1D, {ptAxis});

    ue.add(Form("%s/hPtPrimBadVertex", name.c_str()),
           Form("Generated %s (bad vertex);#it{p}_{T} (GeV/#it{c});Counts", symbol.c_str()),
           HistType::kTH1D, {ptAxis});

    ue.add(Form("%s/hPtPrimGen", name.c_str()),
           Form("Generated %s (after physics selection);#it{p}_{T} (GeV/#it{c});Counts", symbol.c_str()),
           HistType::kTH1D, {ptAxis});

    ue.add(Form("%s/hPtPrimRecoEv", name.c_str()),
           Form("Generated %s (reco events);#it{p}_{T} (GeV/#it{c});Counts", symbol.c_str()),
           HistType::kTH1D, {ptAxis});

    ue.add(Form("%s/hPtPrimGoodEv", name.c_str()),
           Form("Generated %s (good events);#it{p}_{T} (GeV/#it{c});Counts", symbol.c_str()),
           HistType::kTH1D, {ptAxis});

    // 2D versions (vs multiplicity class)
    ue.add(Form("%s/hPtPrimGenAllVsMult", name.c_str()),
           Form("All generated %s vs mult;#it{p}_{T} (GeV/#it{c});Mult Class (%%)", symbol.c_str()),
           HistType::kTH2D, {ptAxis, multAxis});

    ue.add(Form("%s/hPtPrimBadVertexVsMult", name.c_str()),
           Form("Generated %s (bad vertex) vs mult;#it{p}_{T} (GeV/#it{c});Mult Class (%%)", symbol.c_str()),
           HistType::kTH2D, {ptAxis, multAxis});

    ue.add(Form("%s/hPtPrimGenVsMult", name.c_str()),
           Form("Generated %s (after phys sel) vs mult;#it{p}_{T} (GeV/#it{c});Mult Class (%%)", symbol.c_str()),
           HistType::kTH2D, {ptAxis, multAxis});

    ue.add(Form("%s/hPtPrimRecoEvVsMult", name.c_str()),
           Form("Generated %s (reco events) vs mult;#it{p}_{T} (GeV/#it{c});Mult Class (%%)", symbol.c_str()),
           HistType::kTH2D, {ptAxis, multAxis});

    ue.add(Form("%s/hPtPrimGoodEvVsMult", name.c_str()),
           Form("Generated %s (good events) vs mult;#it{p}_{T} (GeV/#it{c});Mult Class (%%)", symbol.c_str()),
           HistType::kTH2D, {ptAxis, multAxis});

    // Tracking efficiency
    ue.add(Form("%s/hPtNumEff", name.c_str()),
           Form("%s tracking efficiency numerator;#it{p}_{T} (GeV/#it{c});Counts", symbol.c_str()),
           HistType::kTH1D, {ptAxis});
    ue.add(Form("%s/hPtNumEffVsMult", name.c_str()),
           Form("%s tracking eff numerator vs mult;#it{p}_{T} (GeV/#it{c});Mult Class (%%)", symbol.c_str()),
           HistType::kTH2D, {ptAxis, multAxis});

    ue.add(Form("%s/hPtDenEff", name.c_str()),
           Form("%s tracking efficiency denominator;#it{p}_{T} (GeV/#it{c});Counts", symbol.c_str()),
           HistType::kTH1D, {ptAxis});
    ue.add(Form("%s/hPtDenEffVsMult", name.c_str()),
           Form("%s tracking eff denominator vs mult;#it{p}_{T} (GeV/#it{c});Mult Class (%%)", symbol.c_str()),
           HistType::kTH2D, {ptAxis, multAxis});

    // Primary fraction
    ue.add(Form("%s/hPtAllReco", name.c_str()),
           Form("All reconstructed %s;#it{p}_{T} (GeV/#it{c});Counts", symbol.c_str()),
           HistType::kTH1D, {ptAxis});
    ue.add(Form("%s/hPtAllRecoVsMult", name.c_str()),
           Form("All reconstructed %s vs mult;#it{p}_{T} (GeV/#it{c});Mult Class (%%)", symbol.c_str()),
           HistType::kTH2D, {ptAxis, multAxis});

    ue.add(Form("%s/hPtPrimReco", name.c_str()),
           Form("Reconstructed primary %s;#it{p}_{T} (GeV/#it{c});Counts", symbol.c_str()),
           HistType::kTH1D, {ptAxis});
    ue.add(Form("%s/hPtPrimRecoVsMult", name.c_str()),
           Form("Reconstructed primary %s vs mult;#it{p}_{T} (GeV/#it{c});Mult Class (%%)", symbol.c_str()),
           HistType::kTH2D, {ptAxis, multAxis});

    ue.add(Form("%s/hPtSecReco", name.c_str()),
           Form("Reconstructed secondary %s;#it{p}_{T} (GeV/#it{c});Counts", symbol.c_str()),
           HistType::kTH1D, {ptAxis});
    ue.add(Form("%s/hPtSecRecoVsMult", name.c_str()),
           Form("Reconstructed secondary %s vs mult;#it{p}_{T} (GeV/#it{c});Mult Class (%%)", symbol.c_str()),
           HistType::kTH2D, {ptAxis, multAxis});

    // Measured spectra
    ue.add(Form("%s/hPtMeasured", name.c_str()),
           Form("Measured %s;#it{p}_{T} (GeV/#it{c});Counts", symbol.c_str()),
           HistType::kTH1D, {ptAxis});
    ue.add(Form("%s/hPtMeasuredVsMult", name.c_str()),
           Form("Measured %s vs mult;#it{p}_{T} (GeV/#it{c});Mult Class (%%)", symbol.c_str()),
           HistType::kTH2D, {ptAxis, multAxis});

    // PID quality
    if (enablePIDHistograms) {
      ue.add(Form("%s/hNsigmaTPC", name.c_str()),
             Form("TPC n#sigma %s;#it{p}_{T} (GeV/#it{c});n#sigma_{TPC}", symbol.c_str()),
             HistType::kTH2D, {ptAxis, {200, -10, 10}});
    }
  }

  // ========================================================================
  // PHI CUT MONITORING
  // ========================================================================

  if (applyPhiCut.value) {
    ue.add("PhiCut/hPtVsPhiPrimeBefore", "pT vs φ' before cut;p_{T} (GeV/c);φ'",
           HistType::kTH2F, {{100, 0, 10}, {100, 0, 0.4}});
    ue.add("PhiCut/hPtVsPhiPrimeAfter", "pT vs φ' after cut;p_{T} (GeV/c);φ'",
           HistType::kTH2F, {{100, 0, 10}, {100, 0, 0.4}});
    ue.add("PhiCut/hRejectionRate", "Track rejection rate by phi cut;p_{T} (GeV/c);Rejection Rate",
           HistType::kTProfile, {{100, 0, 10}});
  }

  // ========================================================================
  // EVENT SELECTION HISTOGRAM
  // ========================================================================

  constexpr int nEvSelBins = 20;
  constexpr float evSelMin = 0.5f;
  constexpr float evSelMax = 20.5f;
  ue.add("evsel", "Event selection", HistType::kTH1D, {{nEvSelBins, evSelMin, evSelMax}});
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
  ue.add("hPhi", "Track phi;#varphi (rad);Counts", HistType::kTH1D, {{64, 0, 2.0 * M_PI}});
  ue.add("hvtxZ", "Vertex Z (data);Vertex Z (cm);Events", HistType::kTH1F, {{40, -20.0, 20.0}});
  ue.add("hvtxZmc", "MC vertex Z;Vertex Z (cm);Events", HistType::kTH1F, {{40, -20.0, 20.0}});

  LOG(info) << "=== Initialized MultiplicityPt task with ON-THE-FLY PERCENTILE COMPUTATION ===";
  LOG(info) << "Centrality classes: " << kCentralityClasses;
  LOG(info) << "Multiplicity estimator: " << multiplicityEstimator.value;
  LOG(info) << "IMPORTANT: Run processPercentileCalibration FIRST to build percentile boundaries!";
  if (applyPhiCut.value) {
    LOG(info) << "Phi cut ENABLED for pT > " << pTthresholdPhiCut.value << " GeV/c";
  }
}

// ========================================================================
// PERCENTILE CALIBRATION PASS
// ========================================================================
void MultiplicityPt::processPercentileCalibration(CollisionTableMCTrue const& mcCollisions,
                                                  ParticleTableMC const& particles)
{
  LOG(info) << "=== PERCENTILE CALIBRATION PASS ===";
  LOG(info) << "Processing " << mcCollisions.size() << " MC collisions";

  multiplicityValues.clear();
  multiplicityValues.reserve(mcCollisions.size());

  for (const auto& mcCollision : mcCollisions) {
    // Apply basic cuts
    if (std::abs(mcCollision.posZ()) > cfgCutVertex.value)
      continue;

    auto particlesInCollision = particles.sliceBy(perMCCol, mcCollision.globalIndex());

    // Apply INEL cuts
    if (cfgINELCut.value == 1 && !o2::pwglf::isINELgt0mc(particlesInCollision, pdg))
      continue;
    if (cfgINELCut.value == 2 && !o2::pwglf::isINELgt1mc(particlesInCollision, pdg))
      continue;

    // Calculate multiplicity
    float mcMult = getMultiplicityMC(mcCollision, particles);
    multiplicityValues.push_back(mcMult);

    ue.fill(HIST("Calibration/hRawMultiplicity"), mcMult);
  }

  // Compute percentile boundaries
  computePercentileBoundaries();

  LOG(info) << "=== PERCENTILE CALIBRATION COMPLETE ===";
  LOG(info) << "Processed " << multiplicityValues.size() << " events";
  LOG(info) << "Now run processMC and processTrue with these percentiles";
}

// ========================================================================
// DATA PROCESSING
// ========================================================================
void MultiplicityPt::processData(CollisionTableData::iterator const& collision,
                                 TrackTableData const& tracks,
                                 BCsRun3 const& /*bcs*/)
{
  if (!isEventSelected<true>(collision)) {
    return;
  }
  ue.fill(HIST("hvtxZ"), collision.posZ());

  float magField = 0;
  if (applyPhiCut.value) {
    const auto& bc = collision.bc_as<BCsRun3>();
    magField = getMagneticField(bc.timestamp());
  }

  for (const auto& track : tracks) {
    if (applyPhiCut.value && track.pt() >= pTthresholdPhiCut.value) {
      float phiPrime = getTransformedPhi(track.phi(), track.sign(), magField);
      ue.fill(HIST("PhiCut/hPtVsPhiPrimeBefore"), track.pt(), phiPrime);
    }

    if (!passesTrackSelection(track, magField)) {
      continue;
    }

    if (applyPhiCut.value && track.pt() >= pTthresholdPhiCut.value) {
      float phiPrime = getTransformedPhi(track.phi(), track.sign(), magField);
      ue.fill(HIST("PhiCut/hPtVsPhiPrimeAfter"), track.pt(), phiPrime);
    }

    ue.fill(HIST("Inclusive/hPtMeasured"), track.pt());
    ue.fill(HIST("hEta"), track.eta());
    ue.fill(HIST("hPhi"), track.phi());

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
// MC PROCESSING - Using computed percentiles
// ========================================================================
void MultiplicityPt::processMC(TrackTableMC const& tracks,
                               aod::McParticles const& particles,
                               CollisionTableMCTrue const& mcCollisions,
                               CollisionTableMC const& collisions,
                               BCsRun3 const& /*bcs*/)
{
  if (!percentilesComputed) {
    LOG(warning) << "Percentiles not computed yet! Run processPercentileCalibration first!";
    LOG(warning) << "Using fallback linear binning for now...";
  }

  LOG(info) << "=== DEBUG processMC START ===";
  LOG(info) << "MC collisions: " << mcCollisions.size();
  LOG(info) << "Reconstructed collisions: " << collisions.size();

  ue.fill(HIST("MC/GenRecoCollisions"), 1.f, mcCollisions.size());
  ue.fill(HIST("MC/GenRecoCollisions"), 2.f, collisions.size());

  std::set<int64_t> physicsSelectedMCCollisions;
  std::set<int64_t> reconstructedMCCollisions;
  std::set<int64_t> selectedMCCollisions;

  std::map<int64_t, float> mcCollisionMultiplicity;
  std::map<int64_t, float> mcCollisionPercentile;

  // First pass: classify MC collisions
  for (const auto& mcCollision : mcCollisions) {
    int64_t mcCollId = mcCollision.globalIndex();

    float mcMult = getMultiplicityMC(mcCollision, particles);
    mcCollisionMultiplicity[mcCollId] = mcMult;

    // Convert to percentile
    float percentile = multiplicityToPercentile(mcMult);
    mcCollisionPercentile[mcCollId] = percentile;

    ue.fill(HIST("MC/EventLoss/MultGenerated"), percentile);

    auto particlesInCollision = particles.sliceBy(perMCCol, mcCollId);

    if (std::abs(mcCollision.posZ()) > cfgCutVertex.value) {
      ue.fill(HIST("MC/EventLoss/MultBadVertex"), percentile);
      continue;
    }

    if (cfgINELCut.value == 1 && !o2::pwglf::isINELgt0mc(particlesInCollision, pdg)) {
      continue;
    }
    if (cfgINELCut.value == 2 && !o2::pwglf::isINELgt1mc(particlesInCollision, pdg)) {
      continue;
    }

    physicsSelectedMCCollisions.insert(mcCollId);
    ue.fill(HIST("MC/EventLoss/MultPhysicsSelected"), percentile);
  }

  LOG(info) << "Physics-selected MC collisions: " << physicsSelectedMCCollisions.size();

  // Second pass: track reconstructed events
  std::set<int64_t> selectedCollisionIndices;

  for (const auto& collision : collisions) {
    if (!collision.has_mcCollision()) {
      continue;
    }

    const auto& mcCollision = collision.mcCollision_as<CollisionTableMCTrue>();
    int64_t mcCollId = mcCollision.globalIndex();

    if (physicsSelectedMCCollisions.find(mcCollId) == physicsSelectedMCCollisions.end()) {
      continue;
    }

    float percentile = mcCollisionPercentile[mcCollId];

    if (reconstructedMCCollisions.find(mcCollId) == reconstructedMCCollisions.end()) {
      reconstructedMCCollisions.insert(mcCollId);
      ue.fill(HIST("MC/EventLoss/MultReconstructed"), percentile);
    }

    if (isEventSelected<false>(collision)) {
      if (selectedMCCollisions.find(mcCollId) == selectedMCCollisions.end()) {
        selectedMCCollisions.insert(mcCollId);
        ue.fill(HIST("MC/EventLoss/MultRecoSelected"), percentile);
      }
      selectedCollisionIndices.insert(collision.globalIndex());
      ue.fill(HIST("hvtxZ"), collision.posZ());
    }
  }

  LOG(info) << "Reconstructed MC collisions: " << reconstructedMCCollisions.size();
  LOG(info) << "Selected MC collisions: " << selectedMCCollisions.size();

  int nPhysicsSelected = physicsSelectedMCCollisions.size();
  int nReconstructed = reconstructedMCCollisions.size();
  int nSelected = selectedMCCollisions.size();

  if (nPhysicsSelected > 0) {
    ue.fill(HIST("hEventLossBreakdown"), 1, nPhysicsSelected);
    ue.fill(HIST("hEventLossBreakdown"), 2, nReconstructed);
    ue.fill(HIST("hEventLossBreakdown"), 3, nSelected);
    ue.fill(HIST("hEventLossBreakdown"), 4, (nSelected * 100.0 / nPhysicsSelected));
  }

  // Process tracks
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

    if (!collision.has_mcCollision())
      continue;

    const auto& mcCollision = collision.mcCollision_as<CollisionTableMCTrue>();
    float percentile = mcCollisionPercentile[mcCollision.globalIndex()];

    float magField = 0;
    if (applyPhiCut.value) {
      const auto& bc = collision.bc_as<BCsRun3>();
      magField = getMagneticField(bc.timestamp());
    }

    if (!passesTrackSelection(track, magField)) {
      continue;
    }
    tracksPassingSelection++;

    // Inclusive charged particle
    ue.fill(HIST("Inclusive/hPtMeasured"), track.pt());
    ue.fill(HIST("Inclusive/hPtMeasuredVsMult"), track.pt(), percentile);
    ue.fill(HIST("Inclusive/hPtAllReco"), track.pt());
    ue.fill(HIST("Inclusive/hPtAllRecoVsMult"), track.pt(), percentile);
    ue.fill(HIST("hEta"), track.eta());
    ue.fill(HIST("hPhi"), track.phi());

    // Efficiency numerator
    if (track.has_mcParticle()) {
      const auto& particle = track.mcParticle();
      int pdgCode = std::abs(particle.pdgCode());

      if (particle.isPhysicalPrimary()) {
        ue.fill(HIST("Inclusive/hPtNumEff"), particle.pt());
        ue.fill(HIST("Inclusive/hPtNumEffVsMult"), particle.pt(), percentile);
        ue.fill(HIST("Inclusive/hPtPrimReco"), track.pt());
        ue.fill(HIST("Inclusive/hPtPrimRecoVsMult"), track.pt(), percentile);

        if (pdgCode == PDGPion) {
          ue.fill(HIST("Pion/hPtNumEff"), particle.pt());
          ue.fill(HIST("Pion/hPtNumEffVsMult"), particle.pt(), percentile);
        }
        if (pdgCode == PDGKaon) {
          ue.fill(HIST("Kaon/hPtNumEff"), particle.pt());
          ue.fill(HIST("Kaon/hPtNumEffVsMult"), particle.pt(), percentile);
        }
        if (pdgCode == PDGProton) {
          ue.fill(HIST("Proton/hPtNumEff"), particle.pt());
          ue.fill(HIST("Proton/hPtNumEffVsMult"), particle.pt(), percentile);
        }
      } else {
        ue.fill(HIST("Inclusive/hPtSecReco"), track.pt());
        ue.fill(HIST("Inclusive/hPtSecRecoVsMult"), track.pt(), percentile);
      }
    }

    // Identified particle analysis
    int bestSpecies = getBestPIDHypothesis(track);

    if (bestSpecies == kPion) {
      ue.fill(HIST("Pion/hPtMeasured"), track.pt());
      ue.fill(HIST("Pion/hPtMeasuredVsMult"), track.pt(), percentile);
      ue.fill(HIST("Pion/hPtAllReco"), track.pt());
      ue.fill(HIST("Pion/hPtAllRecoVsMult"), track.pt(), percentile);
      particleTracksIdentified[kPion]++;

      if (enablePIDHistograms) {
        ue.fill(HIST("Pion/hNsigmaTPC"), track.pt(), track.tpcNSigmaPi());
      }

      if (track.has_mcParticle()) {
        const auto& particle = track.mcParticle();
        if (particle.isPhysicalPrimary()) {
          ue.fill(HIST("Pion/hPtPrimReco"), track.pt());
          ue.fill(HIST("Pion/hPtPrimRecoVsMult"), track.pt(), percentile);
          particleTracksPrimary[kPion]++;
        } else {
          ue.fill(HIST("Pion/hPtSecReco"), track.pt());
          ue.fill(HIST("Pion/hPtSecRecoVsMult"), track.pt(), percentile);
          particleTracksSecondary[kPion]++;
        }
      }

    } else if (bestSpecies == kKaon) {
      ue.fill(HIST("Kaon/hPtMeasured"), track.pt());
      ue.fill(HIST("Kaon/hPtMeasuredVsMult"), track.pt(), percentile);
      ue.fill(HIST("Kaon/hPtAllReco"), track.pt());
      ue.fill(HIST("Kaon/hPtAllRecoVsMult"), track.pt(), percentile);
      particleTracksIdentified[kKaon]++;

      if (enablePIDHistograms) {
        ue.fill(HIST("Kaon/hNsigmaTPC"), track.pt(), track.tpcNSigmaKa());
      }

      if (track.has_mcParticle()) {
        const auto& particle = track.mcParticle();
        if (particle.isPhysicalPrimary()) {
          ue.fill(HIST("Kaon/hPtPrimReco"), track.pt());
          ue.fill(HIST("Kaon/hPtPrimRecoVsMult"), track.pt(), percentile);
          particleTracksPrimary[kKaon]++;
        } else {
          ue.fill(HIST("Kaon/hPtSecReco"), track.pt());
          ue.fill(HIST("Kaon/hPtSecRecoVsMult"), track.pt(), percentile);
          particleTracksSecondary[kKaon]++;
        }
      }

    } else if (bestSpecies == kProton) {
      ue.fill(HIST("Proton/hPtMeasured"), track.pt());
      ue.fill(HIST("Proton/hPtMeasuredVsMult"), track.pt(), percentile);
      ue.fill(HIST("Proton/hPtAllReco"), track.pt());
      ue.fill(HIST("Proton/hPtAllRecoVsMult"), track.pt(), percentile);
      particleTracksIdentified[kProton]++;

      if (enablePIDHistograms) {
        ue.fill(HIST("Proton/hNsigmaTPC"), track.pt(), track.tpcNSigmaPr());
      }

      if (track.has_mcParticle()) {
        const auto& particle = track.mcParticle();
        if (particle.isPhysicalPrimary()) {
          ue.fill(HIST("Proton/hPtPrimReco"), track.pt());
          ue.fill(HIST("Proton/hPtPrimRecoVsMult"), track.pt(), percentile);
          particleTracksPrimary[kProton]++;
        } else {
          ue.fill(HIST("Proton/hPtSecReco"), track.pt());
          ue.fill(HIST("Proton/hPtSecRecoVsMult"), track.pt(), percentile);
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

  LOG(info) << "=== DEBUG processMC END ===";
}

// ========================================================================
// TRUE MC PROCESSING - Using computed percentiles
// ========================================================================
void MultiplicityPt::processTrue(CollisionTableMCTrue const& mcCollisions,
                                 ParticleTableMC const& particles)
{
  if (!percentilesComputed) {
    LOG(warning) << "Percentiles not computed yet! Run processPercentileCalibration first!";
  }

  LOG(info) << "=== DEBUG processTrue START ===";
  LOG(info) << "Number of MC collisions: " << mcCollisions.size();

  int nAllGenerated = 0;
  int nBadVertex = 0;
  int nPhysicsSelected = 0;

  std::array<int, kNSpecies> particleCountAll = {0};
  std::array<int, kNSpecies> particleCountBadVertex = {0};
  std::array<int, kNSpecies> particleCountAfterPS = {0};

  for (const auto& mcCollision : mcCollisions) {
    nAllGenerated++;

    float mcMult = getMultiplicityMC(mcCollision, particles);
    float percentile = multiplicityToPercentile(mcMult);

    ue.fill(HIST("hvtxZmc"), mcCollision.posZ());
    auto particlesInCollision = particles.sliceBy(perMCCol, mcCollision.globalIndex());

    // Fill ALL generated primaries BEFORE any cuts
    for (const auto& particle : particlesInCollision) {
      if (isGoodPrimary(particle)) {
        ue.fill(HIST("Inclusive/hPtPrimGenAll"), particle.pt());
        ue.fill(HIST("Inclusive/hPtPrimGenAllVsMult"), particle.pt(), percentile);
      }

      if (isGoodPrimarySpecies<kPion>(particle)) {
        ue.fill(HIST("Pion/hPtPrimGenAll"), particle.pt());
        ue.fill(HIST("Pion/hPtPrimGenAllVsMult"), particle.pt(), percentile);
        particleCountAll[kPion]++;
      }

      if (isGoodPrimarySpecies<kKaon>(particle)) {
        ue.fill(HIST("Kaon/hPtPrimGenAll"), particle.pt());
        ue.fill(HIST("Kaon/hPtPrimGenAllVsMult"), particle.pt(), percentile);
        particleCountAll[kKaon]++;
      }

      if (isGoodPrimarySpecies<kProton>(particle)) {
        ue.fill(HIST("Proton/hPtPrimGenAll"), particle.pt());
        ue.fill(HIST("Proton/hPtPrimGenAllVsMult"), particle.pt(), percentile);
        particleCountAll[kProton]++;
      }
    }

    // Apply vertex cut
    if (std::abs(mcCollision.posZ()) > cfgCutVertex.value) {
      nBadVertex++;

      for (const auto& particle : particlesInCollision) {
        if (isGoodPrimary(particle)) {
          ue.fill(HIST("Inclusive/hPtPrimBadVertex"), particle.pt());
          ue.fill(HIST("Inclusive/hPtPrimBadVertexVsMult"), particle.pt(), percentile);
        }

        if (isGoodPrimarySpecies<kPion>(particle)) {
          ue.fill(HIST("Pion/hPtPrimBadVertex"), particle.pt());
          ue.fill(HIST("Pion/hPtPrimBadVertexVsMult"), particle.pt(), percentile);
          particleCountBadVertex[kPion]++;
        }

        if (isGoodPrimarySpecies<kKaon>(particle)) {
          ue.fill(HIST("Kaon/hPtPrimBadVertex"), particle.pt());
          ue.fill(HIST("Kaon/hPtPrimBadVertexVsMult"), particle.pt(), percentile);
          particleCountBadVertex[kKaon]++;
        }

        if (isGoodPrimarySpecies<kProton>(particle)) {
          ue.fill(HIST("Proton/hPtPrimBadVertex"), particle.pt());
          ue.fill(HIST("Proton/hPtPrimBadVertexVsMult"), particle.pt(), percentile);
          particleCountBadVertex[kProton]++;
        }
      }
      continue;
    }

    // Apply INEL cuts
    if (cfgINELCut.value == 1 && !o2::pwglf::isINELgt0mc(particlesInCollision, pdg))
      continue;
    if (cfgINELCut.value == 2 && !o2::pwglf::isINELgt1mc(particlesInCollision, pdg))
      continue;

    nPhysicsSelected++;

    // Fill primaries AFTER physics selection (denominator for efficiency)
    for (const auto& particle : particlesInCollision) {
      if (isGoodPrimary(particle)) {
        ue.fill(HIST("Inclusive/hPtDenEff"), particle.pt());
        ue.fill(HIST("Inclusive/hPtDenEffVsMult"), particle.pt(), percentile);
        ue.fill(HIST("Inclusive/hPtPrimGen"), particle.pt());
        ue.fill(HIST("Inclusive/hPtPrimGenVsMult"), particle.pt(), percentile);
      }

      if (isGoodPrimarySpecies<kPion>(particle)) {
        ue.fill(HIST("Pion/hPtDenEff"), particle.pt());
        ue.fill(HIST("Pion/hPtDenEffVsMult"), particle.pt(), percentile);
        ue.fill(HIST("Pion/hPtPrimGen"), particle.pt());
        ue.fill(HIST("Pion/hPtPrimGenVsMult"), particle.pt(), percentile);
        particleCountAfterPS[kPion]++;
      }

      if (isGoodPrimarySpecies<kKaon>(particle)) {
        ue.fill(HIST("Kaon/hPtDenEff"), particle.pt());
        ue.fill(HIST("Kaon/hPtDenEffVsMult"), particle.pt(), percentile);
        ue.fill(HIST("Kaon/hPtPrimGen"), particle.pt());
        ue.fill(HIST("Kaon/hPtPrimGenVsMult"), particle.pt(), percentile);
        particleCountAfterPS[kKaon]++;
      }

      if (isGoodPrimarySpecies<kProton>(particle)) {
        ue.fill(HIST("Proton/hPtDenEff"), particle.pt());
        ue.fill(HIST("Proton/hPtDenEffVsMult"), particle.pt(), percentile);
        ue.fill(HIST("Proton/hPtPrimGen"), particle.pt());
        ue.fill(HIST("Proton/hPtPrimGenVsMult"), particle.pt(), percentile);
        particleCountAfterPS[kProton]++;
      }
    }
  }

  LOG(info) << "=== DEBUG processTrue END ===";
  LOG(info) << "All generated events: " << nAllGenerated;
  LOG(info) << "Events with bad vertex: " << nBadVertex;
  LOG(info) << "Passing physics selection: " << nPhysicsSelected;

  LOG(info) << "=== PARTICLE-SPECIFIC STATISTICS ===";
  LOG(info) << "Pions - All: " << particleCountAll[kPion]
            << ", Bad vertex: " << particleCountBadVertex[kPion]
            << ", After PS: " << particleCountAfterPS[kPion];
  LOG(info) << "Kaons - All: " << particleCountAll[kKaon]
            << ", Bad vertex: " << particleCountBadVertex[kKaon]
            << ", After PS: " << particleCountAfterPS[kKaon];
  LOG(info) << "Protons - All: " << particleCountAll[kProton]
            << ", Bad vertex: " << particleCountBadVertex[kProton]
            << ", After PS: " << particleCountAfterPS[kProton];
}
