// Copyright 2019-2022 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file femtoUniverseEfficiencyBase.cxx
/// \brief Tasks that reads the track tables used for the pairing and builds pairs of two tracks
/// \author Zuzanna Chochulska, WUT Warsaw & CTU Prague, zchochul@cern.ch
/// \author Alicja Płachta, WUT Warsaw, alicja.plachta@cern.ch

#include <vector>
#include <random>
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/StepTHn.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "TDatabasePDG.h"

#include "PWGCF/FemtoUniverse/DataModel/FemtoDerived.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseParticleHisto.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseEventHisto.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniversePairCleaner.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseContainer.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseDetaDphiStar.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUtils.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseTrackSelection.h"

using namespace o2;
using namespace o2::analysis::femtoUniverse;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

namespace
{
static constexpr int nPart = 2;
static constexpr int nCuts = 5;
static const std::vector<std::string> partNames{"PartOne", "PartTwo"};
static const std::vector<std::string> cutNames{"MaxPt", "PIDthr", "nSigmaTPC", "nSigmaTPCTOF", "MaxP"};
static const float cutsTable[nPart][nCuts]{
  {4.05f, 1.f, 3.f, 3.f, 100.f},
  {4.05f, 1.f, 3.f, 3.f, 100.f}};
} // namespace

struct femtoUniverseEfficiencyBase {
  SliceCache cache;
  using FemtoFullParticles = soa::Join<aod::FDParticles, aod::FDExtParticles, aod::FDMCLabels>;
  Preslice<FemtoFullParticles> perCol = aod::femtouniverseparticle::fdCollisionId;

  /// Particle selection part
  /// Configurables for both particles
  Configurable<int> ConfNspecies{"ConfNspecies", 2, "Number of particle spieces with PID info"};
  Configurable<LabeledArray<float>> ConfCutTable{"ConfCutTable", {cutsTable[0], nPart, nCuts, partNames, cutNames}, "Particle selections"};
  Configurable<bool> ConfUse3D{"ConfUse3D", false, "Enable three dimensional histogramms (to be used only for analysis with high statistics): k* vs mT vs multiplicity"};
  Configurable<float> ConfEtaMax{"ConfEtaMax", 0.8f, "Higher limit for |Eta| (the same for both particles)"};

  /// Particle 1
  Configurable<int32_t> ConfPDGCodePartOne{"ConfPDGCodePartOne", 2212, "Particle 1 - PDG code"};
  Configurable<uint8_t> ConfParticleTypePartOne{"ConfParticleTypePartOne", aod::femtouniverseparticle::ParticleType::kTrack, "Particle 1 - particle type"};
  Configurable<bool> ConfNoPDGPartOne{"ConfNoPDGPartOne", false, "0: selecting part by PDG, 1: no PID selection"};
  Configurable<float> ConfPtLowPart1{"ConfPtLowPart1", 0.2, "Lower limit for Pt for the first particle"};
  Configurable<float> ConfPtHighPart1{"ConfPtHighPart1", 2.5, "Higher limit for Pt for the first particle"};

  /// Partition for particle 1
  Partition<FemtoFullParticles> partsOneMCGen = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kMCTruthTrack)) && aod::femtouniverseparticle::pt < ConfPtHighPart1 && aod::femtouniverseparticle::pt > ConfPtLowPart1&& nabs(aod::femtouniverseparticle::eta) < ConfEtaMax;

  Partition<FemtoFullParticles> partsTrackOneMCReco = aod::femtouniverseparticle::pt < ConfPtHighPart1 && aod::femtouniverseparticle::pt > ConfPtLowPart1&& nabs(aod::femtouniverseparticle::eta) < ConfEtaMax;

  /// Histogramming for particle 1
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kMCTruthTrack, 1> trackHistoPartOneGen;
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kTrack, 1> trackHistoPartOneRec;
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kV0, 1> trackHistoV0OneRec;
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kV0Child, 3> trackHistoV0OneChildPosRec;
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kV0Child, 4> trackHistoV0OneChildNegRec;

  /// Particle 2
  Configurable<bool> ConfIsSame{"ConfIsSame", false, "Pairs of the same particle"};
  Configurable<int32_t> ConfPDGCodePartTwo{"ConfPDGCodePartTwo", 333, "Particle 2 - PDG code"};
  Configurable<uint8_t> ConfParticleTypePartTwo{"ConfParticleTypePartTwo", aod::femtouniverseparticle::ParticleType::kTrack, "Particle 2 - particle type"};
  Configurable<bool> ConfNoPDGPartTwo{"ConfNoPDGPartTwo", false, "0: selecting part by PDG, 1: no PID selection"};
  Configurable<float> ConfPtLowPart2{"ConfPtLowPart2", 0.2, "Lower limit for Pt for the second particle"};
  Configurable<float> ConfPtHighPart2{"ConfPtHighPart2", 2.5, "Higher limit for Pt for the second particle"};

  /// Partition for particle 2
  Partition<FemtoFullParticles> partsTwoGen = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kMCTruthTrack)) && aod::femtouniverseparticle::pt < ConfPtHighPart2 && aod::femtouniverseparticle::pt > ConfPtLowPart2&& nabs(aod::femtouniverseparticle::eta) < ConfEtaMax;

  Partition<FemtoFullParticles> partsTrackTwoMCReco = aod::femtouniverseparticle::pt < ConfPtHighPart2 && aod::femtouniverseparticle::pt > ConfPtLowPart2&& nabs(aod::femtouniverseparticle::eta) < ConfEtaMax;

  /// Histogramming for particle 2
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kMCTruthTrack, 2> trackHistoPartTwoGen;
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kTrack, 2> trackHistoPartTwoRec;
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kV0, 2> trackHistoV0TwoRec;
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kV0Child, 3> trackHistoV0TwoChildPosRec;
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kV0Child, 4> trackHistoV0TwoChildNegRec;

  /// Histogramming for Event
  FemtoUniverseEventHisto eventHisto;

  /// The configurables need to be passed to an std::vector
  int vPIDPartOne, vPIDPartTwo;
  std::vector<float> kNsigma;

  /// particle part
  ConfigurableAxis ConfTempFitVarpTBins{"ConfTempFitVarpTBins", {20, 0.5, 4.05}, "pT binning of the pT vs. TempFitVar plot"};
  ConfigurableAxis ConfTempFitVarPDGBins{"ConfDTempFitVarInvMassBins", {6000, -2300, 2300}, "binning of the TempFitVar in the pT vs. TempFitVar plot"};

  /// Correlation part
  ConfigurableAxis ConfMultBins{"ConfMultBins", {VARIABLE_WIDTH, 0.0f, 4.0f, 8.0f, 12.0f, 16.0f, 20.0f, 24.0f, 28.0f, 32.0f, 36.0f, 40.0f, 44.0f, 48.0f, 52.0f, 56.0f, 60.0f, 64.0f, 68.0f, 72.0f, 76.0f, 80.0f, 84.0f, 88.0f, 92.0f, 96.0f, 100.0f, 200.0f, 99999.f}, "Mixing bins - multiplicity"}; // \todo to be obtained from the hash task
  ConfigurableAxis ConfVtxBins{"ConfVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};

  ConfigurableAxis ConfmTBins3D{"ConfmTBins3D", {VARIABLE_WIDTH, 1.02f, 1.14f, 1.20f, 1.26f, 1.38f, 1.56f, 1.86f, 4.50f}, "mT Binning for the 3Dimensional plot: k* vs multiplicity vs mT (set <<ConfUse3D>> to true in order to use)"};
  ConfigurableAxis ConfmultBins3D{"ConfmultBins3D", {VARIABLE_WIDTH, 0.0f, 20.0f, 30.0f, 40.0f, 99999.0f}, "multiplicity Binning for the 3Dimensional plot: k* vs multiplicity vs mT (set <<ConfUse3D>> to true in order to use)"};

  ColumnBinningPolicy<aod::collision::PosZ, aod::femtouniversecollision::MultNtr> colBinning{{ConfVtxBins, ConfMultBins}, true};

  ConfigurableAxis ConfkstarBins{"ConfkstarBins", {1500, 0., 6.}, "binning kstar"};
  ConfigurableAxis ConfkTBins{"ConfkTBins", {150, 0., 9.}, "binning kT"};
  ConfigurableAxis ConfmTBins{"ConfmTBins", {225, 0., 7.5}, "binning mT"};
  Configurable<int> ConfNEventsMix{"ConfNEventsMix", 5, "Number of events for mixing"};
  Configurable<bool> ConfIsCPR{"ConfIsCPR", true, "Close Pair Rejection"};
  Configurable<bool> ConfCPRPlotPerRadii{"ConfCPRPlotPerRadii", false, "Plot CPR per radii"};
  Configurable<float> ConfCPRdeltaPhiCutMax{"ConfCPRdeltaPhiCutMax", 0.0, "Delta Phi max cut for Close Pair Rejection"};
  Configurable<float> ConfCPRdeltaPhiCutMin{"ConfCPRdeltaPhiCutMin", 0.0, "Delta Phi min cut for Close Pair Rejection"};
  Configurable<float> ConfCPRdeltaEtaCutMax{"ConfCPRdeltaEtaCutMax", 0.0, "Delta Eta max cut for Close Pair Rejection"};
  Configurable<float> ConfCPRdeltaEtaCutMin{"ConfCPRdeltaEtaCutMin", 0.0, "Delta Eta min cut for Close Pair Rejection"};
  Configurable<float> ConfCPRChosenRadii{"ConfCPRChosenRadii", 0.80, "Delta Eta cut for Close Pair Rejection"};
  Configurable<int> ConfPhiBins{"ConfPhiBins", 29, "Number of phi bins in deta dphi"};
  Configurable<int> ConfEtaBins{"ConfEtaBins", 29, "Number of eta bins in deta dphi"};

  struct : o2::framework::ConfigurableGroup {
    Configurable<float> ConfNsigmaCombinedProton{"ConfNsigmaCombinedProton", 3.0, "TPC and TOF Proton Sigma (combined) for momentum > 0.5"};
    Configurable<float> ConfNsigmaTPCProton{"ConfNsigmaTPCProton", 3.0, "TPC Proton Sigma for momentum < 0.5"};
    Configurable<float> ConfNsigmaCombinedPion{"ConfNsigmaCombinedPion", 3.0, "TPC and TOF Pion Sigma (combined) for momentum > 0.5"};
    Configurable<float> ConfNsigmaTPCPion{"ConfNsigmaTPCPion", 3.0, "TPC Pion Sigma for momentum < 0.5"};
  } ConfBothTracks;

  // Lambda cuts
  Configurable<float> ConfV0InvMassLowLimit{"ConfV0InvV0MassLowLimit", 1.10, "Lower limit of the V0 invariant mass"};
  Configurable<float> ConfV0InvMassUpLimit{"ConfV0InvV0MassUpLimit", 1.13, "Upper limit of the V0 invariant mass"};

  FemtoUniverseContainer<femtoUniverseContainer::EventType::same, femtoUniverseContainer::Observable::kstar> sameEventCont;
  FemtoUniverseTrackSelection trackCuts;
  /// Histogram output
  HistogramRegistry qaRegistry{"TrackQA", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry registryPDG{"PDGHistos", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};
  HistogramRegistry resultRegistry{"Correlations", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry MixQaRegistry{"MixQaRegistry", {}, OutputObjHandlingPolicy::AnalysisObject};

  FemtoUniversePairCleaner<aod::femtouniverseparticle::ParticleType::kV0, aod::femtouniverseparticle::ParticleType::kV0> pairCleaner;
  FemtoUniverseDetaDphiStar<aod::femtouniverseparticle::ParticleType::kV0, aod::femtouniverseparticle::ParticleType::kV0> pairCloseRejection;

  /// @brief Counter for particle swapping
  int fNeventsProcessed = 0;

  void init(InitContext&)
  {

    eventHisto.init(&qaRegistry);
    trackHistoPartOneGen.init(&qaRegistry, ConfTempFitVarpTBins, ConfTempFitVarPDGBins, 0, ConfPDGCodePartOne, false);
    trackHistoPartOneRec.init(&qaRegistry, ConfTempFitVarpTBins, ConfTempFitVarPDGBins, 0, ConfPDGCodePartOne, false);
    trackHistoV0OneRec.init(&qaRegistry, ConfTempFitVarpTBins, ConfTempFitVarPDGBins, 0, ConfPDGCodePartOne, true);
    trackHistoV0OneChildPosRec.init(&qaRegistry, ConfTempFitVarpTBins, ConfTempFitVarPDGBins, 0, 0, true, "posChildV0_1");
    trackHistoV0OneChildNegRec.init(&qaRegistry, ConfTempFitVarpTBins, ConfTempFitVarPDGBins, 0, 0, true, "negChildV0_1");

    registryPDG.add("part1/PDGvspT", "PDG;#it{p}_{T} (GeV/c); PDG", {HistType::kTH2F, {{500, 0, 5}, {16001, -8000.5, 8000.5}}});
    registryPDG.add("part1/dpositive/PDGvspT", "PDG;#it{p}_{T} (GeV/c); PDG", {HistType::kTH2F, {{500, 0, 5}, {16001, -8000.5, 8000.5}}});
    registryPDG.add("part1/dnegative/PDGvspT", "PDG;#it{p}_{T} (GeV/c); PDG", {HistType::kTH2F, {{500, 0, 5}, {16001, -8000.5, 8000.5}}});

    if (!ConfIsSame) {
      trackHistoPartTwoGen.init(&qaRegistry, ConfTempFitVarpTBins, ConfTempFitVarPDGBins, 0, ConfPDGCodePartTwo, false);
      trackHistoPartTwoRec.init(&qaRegistry, ConfTempFitVarpTBins, ConfTempFitVarPDGBins, 0, ConfPDGCodePartTwo, false);
      trackHistoV0TwoRec.init(&qaRegistry, ConfTempFitVarpTBins, ConfTempFitVarPDGBins, 0, ConfPDGCodePartTwo, true);
      trackHistoV0TwoChildPosRec.init(&qaRegistry, ConfTempFitVarpTBins, ConfTempFitVarPDGBins, 0, 0, true, "posChildV0_2");
      trackHistoV0TwoChildNegRec.init(&qaRegistry, ConfTempFitVarpTBins, ConfTempFitVarPDGBins, 0, 0, true, "negChildV0_2");
      registryPDG.addClone("part1/", "part2/");
    }

    MixQaRegistry.add("MixingQA/hSECollisionBins", ";bin;Entries", kTH1F, {{120, -0.5, 119.5}});

    sameEventCont.init(&resultRegistry, ConfkstarBins, ConfMultBins, ConfkTBins, ConfmTBins, ConfmultBins3D, ConfmTBins3D, ConfEtaBins, ConfPhiBins, 0, ConfUse3D);
    sameEventCont.setPDGCodes(ConfPDGCodePartOne, ConfPDGCodePartTwo);
    pairCleaner.init(&qaRegistry);
    if (ConfIsCPR.value) {
      pairCloseRejection.init(&resultRegistry, &qaRegistry, ConfCPRdeltaPhiCutMin.value, ConfCPRdeltaPhiCutMax.value, ConfCPRdeltaEtaCutMin.value, ConfCPRdeltaEtaCutMax.value, ConfCPRChosenRadii.value, ConfCPRPlotPerRadii.value);
    }
  }

  bool IsProtonNSigma(float mom, float nsigmaTPCPr, float nsigmaTOFPr) // previous version from: https://github.com/alisw/AliPhysics/blob/master/PWGCF/FEMTOSCOPY/AliFemtoUser/AliFemtoMJTrackCut.cxx
  {
    if (mom < 0.5) {
      if (TMath::Abs(nsigmaTPCPr) < ConfBothTracks.ConfNsigmaTPCProton) {
        return true;
      } else {
        return false;
      }
    } else if (mom > 0.4) {
      if (TMath::Hypot(nsigmaTOFPr, nsigmaTPCPr) < ConfBothTracks.ConfNsigmaCombinedProton) {
        return true;
      } else {
        return false;
      }
    }
    return false;
  }

  bool IsKaonNSigma(float mom, float nsigmaTPCK, float nsigmaTOFK)
  {
    if (mom < 0.3) { // 0.0-0.3
      if (TMath::Abs(nsigmaTPCK) < 3.0) {
        return true;
      } else {
        return false;
      }
    } else if (mom < 0.45) { // 0.30 - 0.45
      if (TMath::Abs(nsigmaTPCK) < 2.0) {
        return true;
      } else {
        return false;
      }
    } else if (mom < 0.55) { // 0.45-0.55
      if (TMath::Abs(nsigmaTPCK) < 1.0) {
        return true;
      } else {
        return false;
      }
    } else if (mom < 1.5) { // 0.55-1.5 (now we use TPC and TOF)
      if ((TMath::Abs(nsigmaTOFK) < 3.0) && (TMath::Abs(nsigmaTPCK) < 3.0)) {
        {
          return true;
        }
      } else {
        return false;
      }
    } else if (mom > 1.5) { // 1.5 -
      if ((TMath::Abs(nsigmaTOFK) < 2.0) && (TMath::Abs(nsigmaTPCK) < 3.0)) {
        return true;
      } else {
        return false;
      }
    } else {
      return false;
    }
  }

  bool IsPionNSigma(float mom, float nsigmaTPCPi, float nsigmaTOFPi)
  {
    if (true) {
      if (mom < 0.5) {
        if (TMath::Abs(nsigmaTPCPi) < ConfBothTracks.ConfNsigmaTPCPion) {
          return true;
        } else {
          return false;
        }
      } else if (mom > 0.5) {
        if (TMath::Hypot(nsigmaTOFPi, nsigmaTPCPi) < ConfBothTracks.ConfNsigmaCombinedPion) {
          return true;
        } else {
          return false;
        }
      }
    }
    return false;
  }

  bool IsParticleNSigma(float mom, float nsigmaTPCPr, float nsigmaTOFPr, float nsigmaTPCPi, float nsigmaTOFPi, float nsigmaTPCK, float nsigmaTOFK)
  {
    switch (ConfPDGCodePartOne) {
      case 2212:  // Proton
      case -2212: // anty Proton
        return IsProtonNSigma(mom, nsigmaTPCPr, nsigmaTOFPr);
        break;
      case 211:  // Pion
      case -211: // Pion-
        return IsPionNSigma(mom, nsigmaTPCPi, nsigmaTOFPi);
        break;
      case 321:  // Kaon+
      case -321: // Kaon-
        return IsKaonNSigma(mom, nsigmaTPCK, nsigmaTOFK);
        break;
      default:
        return false;
    }
  }

  bool IsParticleTPC(int pdg, float nsigmaTPC)
  {
    switch (pdg) {
      case 2212:  // Proton
      case -2212: // anty Proton
        if (TMath::Abs(nsigmaTPC) < ConfBothTracks.ConfNsigmaTPCProton) {
          return true;
        } else {
          return false;
        }
        break;
      case 211:  // Pion
      case -211: // Pion-
        if (TMath::Abs(nsigmaTPC) < ConfBothTracks.ConfNsigmaTPCPion) {
          return true;
        } else {
          return false;
        }
        break;
      default:
        return false;
    }
  }

  bool invMLambda(float invMassLambda, float invMassAntiLambda)
  {
    if ((invMassLambda < ConfV0InvMassLowLimit || invMassLambda > ConfV0InvMassUpLimit) && (invMassAntiLambda < ConfV0InvMassLowLimit || invMassAntiLambda > ConfV0InvMassUpLimit)) {
      return false;
    }
    return true;
  }

  template <typename CollisionType>
  void fillCollision(CollisionType col)
  {
    MixQaRegistry.fill(HIST("MixingQA/hSECollisionBins"), colBinning.getBin({col.posZ(), col.multNtr()}));
    eventHisto.fillQA(col);
  }

  /// This function processes the same event and takes care of all the histogramming
  /// \todo the trivial loops over the tracks should be factored out since they will be common to all combinations of T-T, T-V0, V0-V0, ...
  /// @tparam PartitionType
  /// @tparam PartType
  /// @tparam isMC: enables Monte Carlo truth specific histograms
  /// @param grouppartsOneMCGen partition for the first particle passed by the process function
  /// @param grouppartsTwoGen partition for the second particle passed by the process function
  /// @param parts femtoUniverseParticles table (in case of Monte Carlo joined with FemtoUniverseMCLabels)
  /// @param magFieldTesla magnetic field of the collision
  /// @param multCol multiplicity of the collision
  template <bool isMC, typename PartitionType>
  void doMCGen(PartitionType grouppartsOneMCGen, PartitionType grouppartsTwoGen, float /*magFieldTesla*/, int multCol)
  {
    bool swpart = fNeventsProcessed % 2;
    fNeventsProcessed++;

    /// Histogramming same event
    for (auto& part : grouppartsOneMCGen) {
      if (!ConfNoPDGPartOne && part.pidcut() != ConfPDGCodePartOne) {
        continue;
      }
      trackHistoPartOneGen.fillQA<isMC, false>(part);
    }

    if (!ConfIsSame) {
      for (auto& part : grouppartsTwoGen) {
        if (!ConfNoPDGPartOne && part.pidcut() != ConfPDGCodePartTwo) {
          continue;
        }
        trackHistoPartTwoGen.fillQA<isMC, false>(part);
      }
    }
    /// Now build the combinations
    for (auto& [p1, p2] : combinations(CombinationsStrictlyUpperIndexPolicy(grouppartsOneMCGen, grouppartsTwoGen))) {

      if (!ConfNoPDGPartOne && (p1.pidcut() != ConfPDGCodePartOne || p2.pidcut() != ConfPDGCodePartTwo)) {
        continue;
      }

      if (swpart)
        sameEventCont.setPair<isMC>(p1, p2, multCol, ConfUse3D);
      else
        sameEventCont.setPair<isMC>(p2, p1, multCol, ConfUse3D);

      swpart = !swpart;
    }
  }

  /// This function processes the same event and takes care of all the histogramming
  /// \todo the trivial loops over the tracks should be factored out since they will be common to all combinations of T-T, T-V0, V0-V0, ...
  /// @tparam PartitionType
  /// @tparam PartType
  /// @tparam isMC: enables Monte Carlo truth specific histograms
  /// @param grouppartsOneMCGen partition for the first particle passed by the process function
  /// @param grouppartsTwoGen partition for the second particle passed by the process function
  /// @param parts femtoUniverseParticles table (in case of Monte Carlo joined with FemtoUniverseMCLabels)
  /// @param magFieldTesla magnetic field of the collision
  /// @param multCol multiplicity of the collision
  template <bool isMC, typename PartitionType>
  void doMCRecTrackTrack(PartitionType grouppartsOneMCGen, PartitionType grouppartsTwoGen, float /*magFieldTesla*/, int multCol)
  {
    bool swpart = fNeventsProcessed % 2;
    fNeventsProcessed++;

    /// Histogramming same event
    for (auto& part : grouppartsOneMCGen) {
      if (!IsParticleNSigma(part.p(), trackCuts.getNsigmaTPC(part, o2::track::PID::Proton), trackCuts.getNsigmaTOF(part, o2::track::PID::Proton), trackCuts.getNsigmaTPC(part, o2::track::PID::Pion), trackCuts.getNsigmaTOF(part, o2::track::PID::Pion), trackCuts.getNsigmaTPC(part, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(part, o2::track::PID::Kaon))) {
        continue;
      }
      trackHistoPartOneRec.fillQA<isMC, false>(part);
    }

    if (!ConfIsSame) {
      for (auto& part : grouppartsTwoGen) {
        if (!IsParticleNSigma(part.p(), trackCuts.getNsigmaTPC(part, o2::track::PID::Proton), trackCuts.getNsigmaTOF(part, o2::track::PID::Proton), trackCuts.getNsigmaTPC(part, o2::track::PID::Pion), trackCuts.getNsigmaTOF(part, o2::track::PID::Pion), trackCuts.getNsigmaTPC(part, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(part, o2::track::PID::Kaon))) {
          continue;
        }

        trackHistoPartTwoRec.fillQA<isMC, false>(part);
      }
    }
    /// Now build the combinations
    for (auto& [p1, p2] : combinations(CombinationsStrictlyUpperIndexPolicy(grouppartsOneMCGen, grouppartsTwoGen))) {
      if (!IsParticleNSigma(p1.p(), trackCuts.getNsigmaTPC(p1, o2::track::PID::Proton), trackCuts.getNsigmaTOF(p1, o2::track::PID::Proton), trackCuts.getNsigmaTPC(p1, o2::track::PID::Pion), trackCuts.getNsigmaTOF(p1, o2::track::PID::Pion), trackCuts.getNsigmaTPC(p1, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(p1, o2::track::PID::Kaon))) {
        continue;
      }
      if (!IsParticleNSigma(p2.p(), trackCuts.getNsigmaTPC(p2, o2::track::PID::Proton), trackCuts.getNsigmaTOF(p2, o2::track::PID::Proton), trackCuts.getNsigmaTPC(p2, o2::track::PID::Pion), trackCuts.getNsigmaTOF(p2, o2::track::PID::Pion), trackCuts.getNsigmaTPC(p2, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(p2, o2::track::PID::Kaon))) {
        continue;
      }
      if (swpart)
        sameEventCont.setPair<isMC>(p1, p2, multCol, ConfUse3D);
      else
        sameEventCont.setPair<isMC>(p2, p1, multCol, ConfUse3D);

      swpart = !swpart;
    }
  }

  /// This function processes the same event and takes care of all the histogramming
  /// \todo the trivial loops over the tracks should be factored out since they will be common to all combinations of T-T, T-V0, V0-V0, ...
  /// @tparam PartitionType
  /// @tparam PartType
  /// @tparam isMC: enables Monte Carlo truth specific histograms
  /// @param grouppartsOneMCGen partition for the first particle passed by the process function
  /// @param grouppartsTwoGen partition for the second particle passed by the process function
  /// @param parts femtoUniverseParticles table (in case of Monte Carlo joined with FemtoUniverseMCLabels)
  /// @param magFieldTesla magnetic field of the collision
  /// @param multCol multiplicity of the collision
  template <bool isMC, typename PartitionType>
  void doMCRecTrackPhi(PartitionType grouppartsOneMCGen, PartitionType grouppartsTwoGen, float /*magFieldTesla*/, int multCol)
  { // part1 is track and part2 is Phi
    bool swpart = fNeventsProcessed % 2;
    fNeventsProcessed++;

    /// Histogramming same event
    for (auto& part : grouppartsOneMCGen) {
      if (!IsParticleNSigma(part.p(), trackCuts.getNsigmaTPC(part, o2::track::PID::Proton), trackCuts.getNsigmaTOF(part, o2::track::PID::Proton), trackCuts.getNsigmaTPC(part, o2::track::PID::Pion), trackCuts.getNsigmaTOF(part, o2::track::PID::Pion), trackCuts.getNsigmaTPC(part, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(part, o2::track::PID::Kaon))) {
        continue;
      }
      trackHistoPartOneRec.fillQA<isMC, false>(part);
    }

    if (!ConfIsSame) {
      for (auto& part : grouppartsTwoGen) {

        trackHistoPartTwoRec.fillQA<isMC, false>(part);
      }
    }
    /// Now build the combinations
    for (auto& [p1, p2] : combinations(CombinationsStrictlyUpperIndexPolicy(grouppartsOneMCGen, grouppartsTwoGen))) {
      if (!IsParticleNSigma(p1.p(), trackCuts.getNsigmaTPC(p1, o2::track::PID::Proton), trackCuts.getNsigmaTOF(p1, o2::track::PID::Proton), trackCuts.getNsigmaTPC(p1, o2::track::PID::Pion), trackCuts.getNsigmaTOF(p1, o2::track::PID::Pion), trackCuts.getNsigmaTPC(p1, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(p1, o2::track::PID::Kaon))) {
        continue;
      }

      if (swpart)
        sameEventCont.setPair<isMC>(p1, p2, multCol, ConfUse3D);
      else
        sameEventCont.setPair<isMC>(p2, p1, multCol, ConfUse3D);

      swpart = !swpart;
    }
  }

  /// This function processes the same event and takes care of all the histogramming
  /// \todo the trivial loops over the tracks should be factored out since they will be common to all combinations of T-T, T-V0, V0-V0, ...
  /// @tparam PartitionType
  /// @tparam PartType
  /// @tparam isMC: enables Monte Carlo truth specific histograms
  /// @param grouppartsOneMCGen partition for the first particle passed by the process function
  /// @param grouppartsTwoGen partition for the second particle passed by the process function
  /// @param parts femtoUniverseParticles table (in case of Monte Carlo joined with FemtoUniverseMCLabels)
  /// @param magFieldTesla magnetic field of the collision
  /// @param multCol multiplicity of the collision
  template <bool isMC, typename PartitionType, typename ParticlesType>
  void doMCRecV0antiV0(PartitionType grouppartsOneMCGen, PartitionType grouppartsTwoGen, float magFieldTesla, int multCol, ParticlesType parts)
  {

    bool swpart = fNeventsProcessed % 2;
    fNeventsProcessed++;

    /// Histogramming same event
    for (auto& part : grouppartsOneMCGen) {

      if (part.partType() != uint8_t(ConfParticleTypePartOne) || !invMLambda(part.mLambda(), part.mAntiLambda())) {
        continue;
      }

      const auto& posChild = parts.iteratorAt(part.index() - 2);
      const auto& negChild = parts.iteratorAt(part.index() - 1);

      /// Daughters that do not pass this condition are not selected
      if (!IsParticleTPC(2212, trackCuts.getNsigmaTPC(posChild, o2::track::PID::Proton)) || !IsParticleTPC(-211, trackCuts.getNsigmaTPC(negChild, o2::track::PID::Pion)))
        continue;

      trackHistoV0OneRec.fillQA<isMC, true>(part);
      trackHistoV0OneChildPosRec.fillQABase<isMC, true>(posChild, HIST("posChildV0_1"));
      trackHistoV0OneChildNegRec.fillQABase<isMC, true>(negChild, HIST("negChildV0_1"));

      if (!posChild.has_fdMCParticle() || !negChild.has_fdMCParticle() || !part.has_fdMCParticle()) {
        continue;
      }
      const auto mcParticle = part.fdMCParticle();
      const auto mcPosChild = posChild.fdMCParticle();
      const auto mcNegChild = negChild.fdMCParticle();

      registryPDG.fill(HIST("part1/PDGvspT"), part.pt(), mcParticle.pdgMCTruth());
      registryPDG.fill(HIST("part1/dpositive/PDGvspT"), part.pt(), mcPosChild.pdgMCTruth());
      registryPDG.fill(HIST("part1/dnegative/PDGvspT"), part.pt(), mcNegChild.pdgMCTruth());
    }

    if (!ConfIsSame) {
      for (auto& part : grouppartsTwoGen) {

        if (part.partType() != uint8_t(ConfParticleTypePartTwo) || !invMLambda(part.mLambda(), part.mAntiLambda())) {
          continue;
        }

        const auto& posChild = parts.iteratorAt(part.index() - 2);
        const auto& negChild = parts.iteratorAt(part.index() - 1);

        /// Daughters that do not pass this condition are not selected
        if (!IsParticleTPC(211, trackCuts.getNsigmaTPC(posChild, o2::track::PID::Pion)) || !IsParticleTPC(-2212, trackCuts.getNsigmaTPC(negChild, o2::track::PID::Proton)))
          continue;

        trackHistoV0TwoRec.fillQA<isMC, true>(part);
        trackHistoV0TwoChildPosRec.fillQABase<isMC, true>(posChild, HIST("posChildV0_2"));
        trackHistoV0TwoChildNegRec.fillQABase<isMC, true>(negChild, HIST("negChildV0_2"));

        if (!posChild.has_fdMCParticle() || !negChild.has_fdMCParticle() || !part.has_fdMCParticle()) {
          continue;
        }
        const auto mcParticle = part.fdMCParticle();
        const auto mcPosChild = posChild.fdMCParticle();
        const auto mcNegChild = negChild.fdMCParticle();

        registryPDG.fill(HIST("part2/PDGvspT"), part.pt(), mcParticle.pdgMCTruth());
        registryPDG.fill(HIST("part2/dpositive/PDGvspT"), part.pt(), mcPosChild.pdgMCTruth());
        registryPDG.fill(HIST("part2/dnegative/PDGvspT"), part.pt(), mcNegChild.pdgMCTruth());
      }
    }

    /// Now build the combinations
    for (auto& [p1, p2] : combinations(CombinationsStrictlyUpperIndexPolicy(grouppartsOneMCGen, grouppartsTwoGen))) {

      if (p1.partType() != uint8_t(ConfParticleTypePartOne) || p2.partType() != uint8_t(ConfParticleTypePartTwo)) {
        continue;
      }

      // Lambda invariant mass cut for p1
      if (!invMLambda(p1.mLambda(), p1.mAntiLambda()))
        continue;
      // Lambda invariant mass cut for p2
      if (!invMLambda(p2.mLambda(), p2.mAntiLambda()))
        continue;

      const auto& posChild1 = parts.iteratorAt(p1.index() - 2);
      const auto& negChild1 = parts.iteratorAt(p1.index() - 1);
      /// Daughters that do not pass this condition are not selected
      if (!IsParticleTPC(2212, trackCuts.getNsigmaTPC(posChild1, o2::track::PID::Proton)) || !IsParticleTPC(-211, trackCuts.getNsigmaTPC(negChild1, o2::track::PID::Pion))) // Lambda
        continue;

      const auto& posChild2 = parts.iteratorAt(p2.index() - 2);
      const auto& negChild2 = parts.iteratorAt(p2.index() - 1);
      /// Daughters that do not pass this condition are not selected
      if (!IsParticleTPC(211, trackCuts.getNsigmaTPC(posChild2, o2::track::PID::Pion)) || !IsParticleTPC(-2212, trackCuts.getNsigmaTPC(negChild2, o2::track::PID::Proton))) // Antilambda
        continue;

      // track cleaning
      if (!pairCleaner.isCleanPair(p1, p2, parts)) {
        continue;
      }
      if (ConfIsCPR.value) {
        if (pairCloseRejection.isClosePair(p1, p2, parts, magFieldTesla, femtoUniverseContainer::EventType::same)) {
          continue;
        }
      }

      if (swpart)
        sameEventCont.setPair<isMC>(p1, p2, multCol, ConfUse3D);
      else
        sameEventCont.setPair<isMC>(p2, p1, multCol, ConfUse3D);

      swpart = !swpart;
    }
  }

  /// process function for to call doMCPlots with Data
  /// \param col subscribe to the collision table (Data)
  /// \param parts subscribe to the femtoUniverseParticleTable
  void processTrackTrack(o2::aod::FDCollision& col,
                         FemtoFullParticles&)
  {
    fillCollision(col);
    // MCGen
    auto thegrouppartsOneMCGen = partsOneMCGen->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    auto thegrouppartsTwoGen = partsTwoGen->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    doMCGen<false>(thegrouppartsOneMCGen, thegrouppartsTwoGen, col.magField(), col.multNtr());
    // MCRec
    auto thegroupPartsTrackOneRec = partsTrackOneMCReco->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    auto thegroupPartsTrackTwoReco = partsTrackTwoMCReco->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    doMCRecTrackTrack<false>(thegroupPartsTrackOneRec, thegroupPartsTrackTwoReco, col.magField(), col.multNtr());
  }
  PROCESS_SWITCH(femtoUniverseEfficiencyBase, processTrackTrack, "Enable processing track-track efficiency task", false);

  /// process function for to call doMCPlots with Data
  /// \param col subscribe to the collision table (Data)
  /// \param parts subscribe to the femtoUniverseParticleTable
  void processTrackPhi(o2::aod::FDCollision& col,
                       FemtoFullParticles&)
  {
    fillCollision(col);
    // MCGen
    auto thegrouppartsOneMCGen = partsOneMCGen->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    auto thegrouppartsTwoGen = partsTwoGen->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    doMCGen<false>(thegrouppartsOneMCGen, thegrouppartsTwoGen, col.magField(), col.multNtr());
    // MCRec
    auto thegroupPartsTrackOneRec = partsTrackOneMCReco->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    auto thegroupPartsTrackTwoReco = partsTrackTwoMCReco->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    doMCRecTrackPhi<false>(thegroupPartsTrackOneRec, thegroupPartsTrackTwoReco, col.magField(), col.multNtr());
  }
  PROCESS_SWITCH(femtoUniverseEfficiencyBase, processTrackPhi, "Enable processing track-phi efficiency task", true);

  /// process function for to call doMCPlots with Data
  /// \param col subscribe to the collision table (Data)
  /// \param parts subscribe to the femtoUniverseParticleTable
  void processV0antiV0(o2::aod::FDCollision& col,
                       FemtoFullParticles& parts, aod::FDMCParticles const& mcparticles)
  {
    fillCollision(col);
    // MCGen
    auto thegrouppartsOneMCGen = partsOneMCGen->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    auto thegrouppartsTwoGen = partsTwoGen->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    doMCGen<false>(thegrouppartsOneMCGen, thegrouppartsTwoGen, col.magField(), col.multNtr());

    // MCRec
    auto thegroupPartsTrackOneRec = partsTrackOneMCReco->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    auto thegroupPartsTrackTwoReco = partsTrackTwoMCReco->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    doMCRecV0antiV0<false>(thegroupPartsTrackOneRec, thegroupPartsTrackTwoReco, col.magField(), col.multNtr(), parts);
  }
  PROCESS_SWITCH(femtoUniverseEfficiencyBase, processV0antiV0, "Enable processing V0-antiV0 efficiency task", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<femtoUniverseEfficiencyBase>(cfgc),
  };
  return workflow;
}
