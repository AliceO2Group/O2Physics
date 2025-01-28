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
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/RunningWorkflowInfo.h"

#include "PWGCF/FemtoUniverse/Core/FemtoUniverseParticleHisto.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseEventHisto.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseTrackSelection.h"

using namespace o2;
using namespace o2::analysis::femto_universe;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

struct femtoUniverseEfficiencyBase {
  SliceCache cache;
  using FemtoFullParticles = soa::Join<aod::FDParticles, aod::FDExtParticles, aod::FDMCLabels>;
  Preslice<FemtoFullParticles> perCol = aod::femtouniverseparticle::fdCollisionId;

  Configurable<bool> ConfIsDebug{"ConfIsDebug", true, "Enable debug histograms"};

  // Collisions
  Configurable<float> ConfZVertex{"ConfZVertex", 10.f, "Event sel: Maximum z-Vertex (cm)"};

  Filter collisionFilter = (nabs(aod::collision::posZ) < ConfZVertex);
  using FilteredFDCollisions = soa::Filtered<o2::aod::FdCollisions>;
  using FilteredFDCollision = FilteredFDCollisions::iterator;

  /// Particle selection part
  /// Configurables for both particles
  ConfigurableAxis ConfTempFitVarpTBins{"ConfTempFitVarpTBins", {20, 0.5, 4.05}, "Binning of the pT in the pT vs. TempFitVar plot"};
  ConfigurableAxis ConfTempFitVarPDGBins{"ConfTempFitVarPDGBins", {6000, -2300, 2300}, "Binning of the PDG code in the pT vs. TempFitVar plot"};
  ConfigurableAxis ConfTempFitVarCPABins{"ConfTempFitVarCPABins", {1000, 0.9, 1}, "Binning of the pointing angle cosinus in the pT vs. TempFitVar plot"};
  ConfigurableAxis ConfTempFitVarDCABins{"ConfTempFitVarDCABins", {1000, -5, 5}, "Binning of the PDG code in the pT vs. TempFitVar plot"};

  struct : o2::framework::ConfigurableGroup {
    Configurable<float> ConfEtaMax{"ConfEtaMax", 0.8f, "Higher limit for |Eta| (the same for both particles)"};
    Configurable<float> ConfMomProton{"ConfMomProton", 0.75, "Momentum threshold for proton identification using TOF"};
    Configurable<float> ConfMomPion{"ConfMomPion", 0.75, "Momentum threshold for pion identification using TOF"};
    Configurable<float> ConfNsigmaCombinedProton{"ConfNsigmaCombinedProton", 3.0, "TPC and TOF Proton Sigma (combined) for momentum > ConfMomProton"};
    Configurable<float> ConfNsigmaTPCProton{"ConfNsigmaTPCProton", 3.0, "TPC Proton Sigma for momentum < ConfMomProton"};
    Configurable<float> ConfNsigmaCombinedPion{"ConfNsigmaCombinedPion", 3.0, "TPC and TOF Pion Sigma (combined) for momentum > ConfMomPion"};
    Configurable<float> ConfNsigmaTPCPion{"ConfNsigmaTPCPion", 3.0, "TPC Pion Sigma for momentum < ConfMomPion"};
  } ConfBothTracks;

  /// Lambda cuts
  Configurable<float> ConfV0InvMassLowLimit{"ConfV0InvV0MassLowLimit", 1.10, "Lower limit of the V0 invariant mass"};
  Configurable<float> ConfV0InvMassUpLimit{"ConfV0InvV0MassUpLimit", 1.13, "Upper limit of the V0 invariant mass"};

  /// Kaon configurable
  Configurable<bool> IsKaonRun2{"IsKaonRun2", false, "Enable kaon selection used in Run2"}; // to check consistency with Run2 results

  /// Deuteron configurables
  struct : o2::framework::ConfigurableGroup {
    Configurable<float> ConfNsigmaTPCDe{"ConfNsigmaTPCDe", 2.0f, "TPC Deuteron Sigma for momentum < ConfTOFpMinDe"};
    Configurable<float> ConfNsigmaTOFDe{"ConfNsigmaTOFDe", 2.0f, "TOF Deuteron Sigma"};
    Configurable<float> ConfTOFpMinDe{"ConfTOFpMinDe", 0.5f, "Min. momentum for deuterons for which TOF is required for PID"};
    Configurable<float> ConfPLowDe{"ConfPLowDe", 0.8f, "Lower limit for momentum for deuterons"};
    Configurable<float> ConfPHighDe{"ConfPHighDe", 1.8f, "Higher limit for momentum for deuterons"};
  } deuteronconfigs;

  /// Particle 1
  Configurable<int32_t> ConfPDGCodePartOne{"ConfPDGCodePartOne", 2212, "Particle 1 - PDG code"};
  Configurable<uint> ConfParticleTypePartOne{"ConfParticleTypePartOne", aod::femtouniverseparticle::ParticleType::kTrack, "Particle 1 - particle type: 0 - track, 2 - V0, 6 - phi"};
  Configurable<bool> ConfNoPDGPartOne{"ConfNoPDGPartOne", false, "0: selecting part one by PDG, 1: no PID selection"};
  Configurable<float> ConfPtLowPart1{"ConfPtLowPart1", 0.2, "Lower limit for Pt for the first particle"};
  Configurable<float> ConfPtHighPart1{"ConfPtHighPart1", 2.5, "Higher limit for Pt for the first particle"};
  Configurable<int> ConfChargePart1{"ConfChargePart1", 1, "Charge of the first particle"};

  /// Partition for particle 1
  Partition<FemtoFullParticles> partsOneMCGen = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kMCTruthTrack)) && aod::femtouniverseparticle::pt < ConfPtHighPart1 && aod::femtouniverseparticle::pt > ConfPtLowPart1&& nabs(aod::femtouniverseparticle::eta) < ConfBothTracks.ConfEtaMax;

  Partition<FemtoFullParticles> partsTrackOneMCReco = aod::femtouniverseparticle::pt < ConfPtHighPart1 && aod::femtouniverseparticle::pt > ConfPtLowPart1&& nabs(aod::femtouniverseparticle::eta) < ConfBothTracks.ConfEtaMax;

  /// Histogramming for particle 1
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kMCTruthTrack, 1> trackHistoPartOneGen;
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kTrack, 1> trackHistoPartOneRec;
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kV0, 1> trackHistoV0OneRec;
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kV0Child, 3> trackHistoV0OneChildPosRec;
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kV0Child, 4> trackHistoV0OneChildNegRec;

  /// Particle 2
  Configurable<bool> ConfIsSame{"ConfIsSame", false, "Pairs of the same particle"};
  Configurable<int32_t> ConfPDGCodePartTwo{"ConfPDGCodePartTwo", 333, "Particle 2 - PDG code"};
  Configurable<uint> ConfParticleTypePartTwo{"ConfParticleTypePartTwo", aod::femtouniverseparticle::ParticleType::kTrack, "Particle 2 - particle type:  0 - track, 2 - V0, 6 - phi"};
  Configurable<bool> ConfNoPDGPartTwo{"ConfNoPDGPartTwo", false, "0: selecting part two by PDG, 1: no PID selection"};
  Configurable<float> ConfPtLowPart2{"ConfPtLowPart2", 0.2, "Lower limit for Pt for the second particle"};
  Configurable<float> ConfPtHighPart2{"ConfPtHighPart2", 2.5, "Higher limit for Pt for the second particle"};
  Configurable<int> ConfChargePart2{"ConfChargePart2", 1, "Charge of the second particle"};

  /// Partition for particle 2
  Partition<FemtoFullParticles> partsTwoGen = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kMCTruthTrack)) && aod::femtouniverseparticle::pt < ConfPtHighPart2 && aod::femtouniverseparticle::pt > ConfPtLowPart2&& nabs(aod::femtouniverseparticle::eta) < ConfBothTracks.ConfEtaMax;

  Partition<FemtoFullParticles> partsTrackTwoMCReco = aod::femtouniverseparticle::pt < ConfPtHighPart2 && aod::femtouniverseparticle::pt > ConfPtLowPart2&& nabs(aod::femtouniverseparticle::eta) < ConfBothTracks.ConfEtaMax;

  /// Histogramming for particle 2
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kMCTruthTrack, 2> trackHistoPartTwoGen;
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kTrack, 2> trackHistoPartTwoRec;
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kV0, 2> trackHistoV0TwoRec;
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kV0Child, 3> trackHistoV0TwoChildPosRec;
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kV0Child, 4> trackHistoV0TwoChildNegRec;

  /// Histogramming for Event
  FemtoUniverseEventHisto eventHisto;

  FemtoUniverseTrackSelection trackCuts;

  /// Histogram output
  HistogramRegistry qaRegistry{"TrackQA", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry registryPDG{"PDGHistos", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};
  HistogramRegistry registryMCOrigin{"MCOriginHistos", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};

  void init(InitContext&)
  {

    eventHisto.init(&qaRegistry);
    trackHistoPartOneGen.init(&qaRegistry, ConfTempFitVarpTBins, ConfTempFitVarPDGBins, 0, ConfPDGCodePartOne, false);
    trackHistoPartOneRec.init(&qaRegistry, ConfTempFitVarpTBins, ConfTempFitVarDCABins, 0, ConfPDGCodePartOne, ConfIsDebug);
    registryMCOrigin.add("part1/hPt", " ;#it{p}_{T} (GeV/c); Entries", {HistType::kTH1F, {{240, 0, 6}}});
    registryPDG.add("part1/PDGvspT", "PDG;#it{p}_{T} (GeV/c); PDG", {HistType::kTH2F, {{500, 0, 5}, {16001, -8000.5, 8000.5}}});
    if (ConfParticleTypePartOne == uint8_t(aod::femtouniverseparticle::ParticleType::kV0)) {
      trackHistoV0OneRec.init(&qaRegistry, ConfTempFitVarpTBins, ConfTempFitVarCPABins, 0, ConfPDGCodePartOne, ConfIsDebug);
      trackHistoV0OneChildPosRec.init(&qaRegistry, ConfTempFitVarpTBins, ConfTempFitVarDCABins, 0, 0, ConfIsDebug, "posChildV0_1");
      trackHistoV0OneChildNegRec.init(&qaRegistry, ConfTempFitVarpTBins, ConfTempFitVarDCABins, 0, 0, ConfIsDebug, "negChildV0_1");
      registryPDG.add("part1/dpositive/PDGvspT", "PDG;#it{p}_{T} (GeV/c); PDG", {HistType::kTH2F, {{500, 0, 5}, {16001, -8000.5, 8000.5}}});
      registryPDG.add("part1/dnegative/PDGvspT", "PDG;#it{p}_{T} (GeV/c); PDG", {HistType::kTH2F, {{500, 0, 5}, {16001, -8000.5, 8000.5}}});
    }

    registryPDG.add("part2/PDGvspT", "PDG;#it{p}_{T} (GeV/c); PDG", {HistType::kTH2F, {{500, 0, 5}, {16001, -8000.5, 8000.5}}});
    if (!ConfIsSame) {
      trackHistoPartTwoGen.init(&qaRegistry, ConfTempFitVarpTBins, ConfTempFitVarPDGBins, 0, ConfPDGCodePartTwo, false);
      trackHistoPartTwoRec.init(&qaRegistry, ConfTempFitVarpTBins, ConfTempFitVarDCABins, 0, ConfPDGCodePartTwo, ConfIsDebug);
      registryMCOrigin.add("part2/hPt", " ;#it{p}_{T} (GeV/c); Entries", {HistType::kTH1F, {{240, 0, 6}}});
      if (ConfParticleTypePartTwo == uint8_t(aod::femtouniverseparticle::ParticleType::kV0)) {
        trackHistoV0TwoRec.init(&qaRegistry, ConfTempFitVarpTBins, ConfTempFitVarCPABins, 0, ConfPDGCodePartTwo, ConfIsDebug);
        trackHistoV0TwoChildPosRec.init(&qaRegistry, ConfTempFitVarpTBins, ConfTempFitVarDCABins, 0, 0, ConfIsDebug, "posChildV0_2");
        trackHistoV0TwoChildNegRec.init(&qaRegistry, ConfTempFitVarpTBins, ConfTempFitVarDCABins, 0, 0, ConfIsDebug, "negChildV0_2");
        registryPDG.add("part2/dpositive/PDGvspT", "PDG;#it{p}_{T} (GeV/c); PDG", {HistType::kTH2F, {{500, 0, 5}, {16001, -8000.5, 8000.5}}});
        registryPDG.add("part2/dnegative/PDGvspT", "PDG;#it{p}_{T} (GeV/c); PDG", {HistType::kTH2F, {{500, 0, 5}, {16001, -8000.5, 8000.5}}});
      }
    }
  }

  bool IsProtonNSigma(float mom, float nsigmaTPCPr, float nsigmaTOFPr) // previous version from: https://github.com/alisw/AliPhysics/blob/master/PWGCF/FEMTOSCOPY/AliFemtoUser/AliFemtoMJTrackCut.cxx
  {
    if (mom < ConfBothTracks.ConfMomProton) {
      if (TMath::Abs(nsigmaTPCPr) < ConfBothTracks.ConfNsigmaTPCProton) {
        return true;
      } else {
        return false;
      }
    } else {
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
    if (IsKaonRun2 == true) {
      if (mom < 0.4) {
        return TMath::Abs(nsigmaTPCK) < 2;
      } else if (mom > 0.4 && mom < 0.45) {
        return TMath::Abs(nsigmaTPCK) < 1;
      } else if (mom > 0.45 && mom < 0.8) {
        return (TMath::Abs(nsigmaTPCK) < 3 && TMath::Abs(nsigmaTOFK) < 2);
      } else if (mom > 0.8 && mom < 1.5) {
        return (TMath::Abs(nsigmaTPCK) < 3 && TMath::Abs(nsigmaTOFK) < 1.5);
      } else {
        return false;
      }
    } else {
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
  }

  bool IsPionNSigma(float mom, float nsigmaTPCPi, float nsigmaTOFPi)
  {
    if (mom < ConfBothTracks.ConfMomPion) {
      if (TMath::Abs(nsigmaTPCPi) < ConfBothTracks.ConfNsigmaTPCPion) {
        return true;
      } else {
        return false;
      }
    } else {
      if (TMath::Hypot(nsigmaTOFPi, nsigmaTPCPi) < ConfBothTracks.ConfNsigmaCombinedPion) {
        return true;
      } else {
        return false;
      }
    }
    return false;
  }

  bool IsDeuteronNSigma(float mom, float nsigmaTPCDe, float nsigmaTOFDe)
  {
    if (mom > deuteronconfigs.ConfPLowDe && mom < deuteronconfigs.ConfPHighDe) {
      if (mom < deuteronconfigs.ConfTOFpMinDe) {
        return (TMath::Abs(nsigmaTPCDe) < deuteronconfigs.ConfNsigmaTPCDe);
      } else {
        return (TMath::Abs(nsigmaTOFDe) < deuteronconfigs.ConfNsigmaTOFDe && (TMath::Abs(nsigmaTPCDe) < deuteronconfigs.ConfNsigmaTPCDe));
      }
    } else {
      return false;
    }
  }

  bool IsParticleNSigma(int pdgCode, float mom, float nsigmaTPCPr, float nsigmaTOFPr, float nsigmaTPCPi, float nsigmaTOFPi, float nsigmaTPCK, float nsigmaTOFK, float nsigmaTPCDe, float nsigmaTOFDe)
  {
    switch (pdgCode) {
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
      case 1000010020:  // Deuteron
      case -1000010020: // Antideuteron
        return IsDeuteronNSigma(mom, nsigmaTPCDe, nsigmaTOFDe);
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
    eventHisto.fillQA(col);
  }

  /// This function processes the same event and takes care of all the histogramming
  /// @tparam PartitionType
  /// @tparam isMC: enables Monte Carlo truth specific histograms
  /// @param grouppartsOneMCGen partition for the first particle passed by the process function
  /// @param grouppartsTwoMCGen partition for the second particle passed by the process function
  template <bool isMC, typename PartitionType>
  void doMCGen(PartitionType grouppartsOneMCGen, PartitionType grouppartsTwoMCGen)
  {
    /// Histogramming same event
    for (auto& part : grouppartsOneMCGen) {
      if (!ConfNoPDGPartOne && part.tempFitVar() != ConfPDGCodePartOne) {
        continue;
      }
      trackHistoPartOneGen.fillQA<isMC, false>(part);
    }

    if (!ConfIsSame) {
      for (auto& part : grouppartsTwoMCGen) {
        if (!ConfNoPDGPartTwo && part.tempFitVar() != ConfPDGCodePartTwo) {
          continue;
        }
        trackHistoPartTwoGen.fillQA<isMC, false>(part);
      }
    }
  }

  /// This function processes the same event and takes care of all the histogramming
  /// @tparam PartitionType
  /// @tparam isMC: enables Monte Carlo truth specific histograms
  /// @tparam isDebug: enables debug histograms
  /// @param grouppartsOneMCRec partition for the first particle passed by the process function
  /// @param grouppartsTwoMCRec partition for the second particle passed by the process function
  template <bool isMC, bool isDebug, typename PartitionType>
  void doMCRecTrackTrack(PartitionType grouppartsOneMCRec, PartitionType grouppartsTwoMCRec)
  {
    /// Histogramming same event
    for (auto& part : grouppartsOneMCRec) {
      if (part.partType() != ConfParticleTypePartOne || part.sign() != ConfChargePart1 || !IsParticleNSigma(ConfPDGCodePartOne, part.p(), trackCuts.getNsigmaTPC(part, o2::track::PID::Proton), trackCuts.getNsigmaTOF(part, o2::track::PID::Proton), trackCuts.getNsigmaTPC(part, o2::track::PID::Pion), trackCuts.getNsigmaTOF(part, o2::track::PID::Pion), trackCuts.getNsigmaTPC(part, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(part, o2::track::PID::Kaon), trackCuts.getNsigmaTPC(part, o2::track::PID::Deuteron), trackCuts.getNsigmaTOF(part, o2::track::PID::Deuteron))) {
        continue;
      }
      trackHistoPartOneRec.fillQA<isMC, isDebug>(part);

      if (!part.has_fdMCParticle()) {
        continue;
      }
      const auto mcParticle = part.fdMCParticle();

      registryPDG.fill(HIST("part1/PDGvspT"), part.pt(), mcParticle.pdgMCTruth());
      registryMCOrigin.fill(HIST("part1/hPt"), mcParticle.pt());
    }

    if (!ConfIsSame) {
      for (auto& part : grouppartsTwoMCRec) {
        if (part.partType() != ConfParticleTypePartTwo || part.sign() != ConfChargePart2 || !IsParticleNSigma(ConfPDGCodePartTwo, part.p(), trackCuts.getNsigmaTPC(part, o2::track::PID::Proton), trackCuts.getNsigmaTOF(part, o2::track::PID::Proton), trackCuts.getNsigmaTPC(part, o2::track::PID::Pion), trackCuts.getNsigmaTOF(part, o2::track::PID::Pion), trackCuts.getNsigmaTPC(part, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(part, o2::track::PID::Kaon), trackCuts.getNsigmaTPC(part, o2::track::PID::Deuteron), trackCuts.getNsigmaTOF(part, o2::track::PID::Deuteron))) {
          continue;
        }

        trackHistoPartTwoRec.fillQA<isMC, isDebug>(part);

        if (!part.has_fdMCParticle()) {
          continue;
        }
        const auto mcParticle = part.fdMCParticle();

        registryPDG.fill(HIST("part2/PDGvspT"), part.pt(), mcParticle.pdgMCTruth());
        registryMCOrigin.fill(HIST("part2/hPt"), mcParticle.pt());
      }
    }
  }

  /// This function processes the same event and takes care of all the histogramming
  /// @tparam PartitionType
  /// @tparam isMC: enables Monte Carlo truth specific histograms
  /// @tparam isDebug: enables debug histograms
  /// @param grouppartsOneMCRec partition for the first particle passed by the process function
  /// @param grouppartsTwoMCRec partition for the second particle passed by the process function
  template <bool isMC, bool isDebug, typename PartitionType>
  void doMCRecTrackPhi(PartitionType grouppartsOneMCRec, PartitionType grouppartsTwoMCRec)
  { // part1 is track and part2 is Phi

    for (auto& part : grouppartsOneMCRec) {
      if (part.partType() != ConfParticleTypePartOne || part.sign() != ConfChargePart1 || !IsParticleNSigma(ConfPDGCodePartOne, part.p(), trackCuts.getNsigmaTPC(part, o2::track::PID::Proton), trackCuts.getNsigmaTOF(part, o2::track::PID::Proton), trackCuts.getNsigmaTPC(part, o2::track::PID::Pion), trackCuts.getNsigmaTOF(part, o2::track::PID::Pion), trackCuts.getNsigmaTPC(part, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(part, o2::track::PID::Kaon), trackCuts.getNsigmaTPC(part, o2::track::PID::Deuteron), trackCuts.getNsigmaTOF(part, o2::track::PID::Deuteron))) {
        continue;
      }
      trackHistoPartOneRec.fillQA<isMC, isDebug>(part);

      if (!part.has_fdMCParticle()) {
        continue;
      }
      const auto mcParticle = part.fdMCParticle();

      registryPDG.fill(HIST("part1/PDGvspT"), part.pt(), mcParticle.pdgMCTruth());
      registryMCOrigin.fill(HIST("part1/hPt"), mcParticle.pt());
    }

    if (!ConfIsSame) {
      for (auto& part : grouppartsTwoMCRec) {
        if (part.partType() != ConfParticleTypePartTwo || part.sign() != ConfChargePart2) {
          continue;
        }

        trackHistoPartTwoRec.fillQA<isMC, isDebug>(part);

        if (!part.has_fdMCParticle()) {
          continue;
        }
        const auto mcParticle = part.fdMCParticle();

        registryPDG.fill(HIST("part2/PDGvspT"), part.pt(), mcParticle.pdgMCTruth());
        registryMCOrigin.fill(HIST("part2/hPt"), mcParticle.pt());
      }
    }
  }

  /// This function processes the same event and takes care of all the histogramming
  /// @tparam PartitionType
  /// @tparam ParticlesType
  /// @tparam isMC: enables Monte Carlo truth specific histograms
  /// @tparam isDebug: enables debug histograms
  /// @param grouppartsOneMCRec partition for the first particle passed by the process function
  /// @param grouppartsTwoMCRec partition for the second particle passed by the process function
  /// @param parts all tracks
  template <bool isMC, bool isDebug, typename PartitionType, typename ParticlesType>
  void doMCRecV0V0(PartitionType grouppartsOneMCRec, PartitionType grouppartsTwoMCRec, ParticlesType parts)
  {
    /// Histogramming same event
    for (auto& part : grouppartsOneMCRec) {

      if (part.partType() != ConfParticleTypePartOne || !invMLambda(part.mLambda(), part.mAntiLambda())) {
        continue;
      }

      const auto& posChild = parts.iteratorAt(part.index() - 2);
      const auto& negChild = parts.iteratorAt(part.index() - 1);

      if (ConfPDGCodePartOne > 0 && (!IsProtonNSigma(0, trackCuts.getNsigmaTPC(posChild, o2::track::PID::Proton), trackCuts.getNsigmaTOF(posChild, o2::track::PID::Proton)) || !IsPionNSigma(0, trackCuts.getNsigmaTPC(negChild, o2::track::PID::Pion), trackCuts.getNsigmaTOF(negChild, o2::track::PID::Pion)))) { // give momentum as 0 to only check TPC nSigma, not combined with TOF
        continue;
      }
      if (ConfPDGCodePartOne < 0 && (!IsProtonNSigma(0, trackCuts.getNsigmaTPC(negChild, o2::track::PID::Proton), trackCuts.getNsigmaTOF(negChild, o2::track::PID::Proton)) || !IsPionNSigma(0, trackCuts.getNsigmaTPC(posChild, o2::track::PID::Pion), trackCuts.getNsigmaTOF(posChild, o2::track::PID::Pion)))) { // give momentum as 0 to only check TPC nSigma, not combined with TOF
        continue;
      }

      trackHistoV0OneRec.fillQA<isMC, isDebug>(part);
      trackHistoV0OneChildPosRec.fillQABase<isMC, isDebug>(posChild, HIST("posChildV0_1"));
      trackHistoV0OneChildNegRec.fillQABase<isMC, isDebug>(negChild, HIST("negChildV0_1"));

      if (!posChild.has_fdMCParticle() || !negChild.has_fdMCParticle() || !part.has_fdMCParticle()) {
        continue;
      }
      const auto mcParticle = part.fdMCParticle();
      const auto mcPosChild = posChild.fdMCParticle();
      const auto mcNegChild = negChild.fdMCParticle();

      registryPDG.fill(HIST("part1/PDGvspT"), part.pt(), mcParticle.pdgMCTruth());
      registryPDG.fill(HIST("part1/dpositive/PDGvspT"), part.pt(), mcPosChild.pdgMCTruth());
      registryPDG.fill(HIST("part1/dnegative/PDGvspT"), part.pt(), mcNegChild.pdgMCTruth());
      registryMCOrigin.fill(HIST("part1/hPt"), mcParticle.pt());
    }

    if (!ConfIsSame) {
      for (auto& part : grouppartsTwoMCRec) {

        if (part.partType() != ConfParticleTypePartTwo || !invMLambda(part.mLambda(), part.mAntiLambda())) {
          continue;
        }

        const auto& posChild = parts.iteratorAt(part.index() - 2);
        const auto& negChild = parts.iteratorAt(part.index() - 1);

        if (ConfPDGCodePartTwo > 0 && (!IsProtonNSigma(0, trackCuts.getNsigmaTPC(posChild, o2::track::PID::Proton), trackCuts.getNsigmaTOF(posChild, o2::track::PID::Proton)) || !IsPionNSigma(0, trackCuts.getNsigmaTPC(negChild, o2::track::PID::Pion), trackCuts.getNsigmaTOF(negChild, o2::track::PID::Pion)))) { // give momentum as 0 to only check TPC nSigma, not combined with TOF
          continue;
        }
        if (ConfPDGCodePartTwo < 0 && (!IsProtonNSigma(0, trackCuts.getNsigmaTPC(negChild, o2::track::PID::Proton), trackCuts.getNsigmaTOF(negChild, o2::track::PID::Proton)) || !IsPionNSigma(0, trackCuts.getNsigmaTPC(posChild, o2::track::PID::Pion), trackCuts.getNsigmaTOF(posChild, o2::track::PID::Pion)))) { // give momentum as 0 to only check TPC nSigma, not combined with TOF
          continue;
        }

        trackHistoV0TwoRec.fillQA<isMC, isDebug>(part);
        trackHistoV0TwoChildPosRec.fillQABase<isMC, isDebug>(posChild, HIST("posChildV0_2"));
        trackHistoV0TwoChildNegRec.fillQABase<isMC, isDebug>(negChild, HIST("negChildV0_2"));

        if (!posChild.has_fdMCParticle() || !negChild.has_fdMCParticle() || !part.has_fdMCParticle()) {
          continue;
        }
        const auto mcParticle = part.fdMCParticle();
        const auto mcPosChild = posChild.fdMCParticle();
        const auto mcNegChild = negChild.fdMCParticle();

        registryPDG.fill(HIST("part2/PDGvspT"), part.pt(), mcParticle.pdgMCTruth());
        registryPDG.fill(HIST("part2/dpositive/PDGvspT"), part.pt(), mcPosChild.pdgMCTruth());
        registryPDG.fill(HIST("part2/dnegative/PDGvspT"), part.pt(), mcNegChild.pdgMCTruth());
        registryMCOrigin.fill(HIST("part2/hPt"), mcParticle.pt());
      }
    }
  }

  /// This function processes the same event and takes care of all the histogramming
  /// @tparam PartitionType
  /// @tparam ParticlesType
  /// @tparam isMC: enables Monte Carlo truth specific histograms
  /// @tparam isDebug: enables debug histograms
  /// @param grouppartsOneMCRec partition for the first particle passed by the process function
  /// @param grouppartsTwoMCRec partition for the second particle passed by the process function
  /// @param parts all tracks
  template <bool isMC, bool isDebug, typename PartitionType, typename ParticlesType>
  void doMCRecTrackV0(PartitionType grouppartsOneMCRec, PartitionType grouppartsTwoMCRec, ParticlesType const& parts)
  { // part1 is track and part2 is V0

    /// Histogramming same event
    for (auto& part : grouppartsOneMCRec) {
      if (part.partType() != ConfParticleTypePartOne || part.sign() != ConfChargePart1 || !IsParticleNSigma(ConfPDGCodePartOne, part.p(), trackCuts.getNsigmaTPC(part, o2::track::PID::Proton), trackCuts.getNsigmaTOF(part, o2::track::PID::Proton), trackCuts.getNsigmaTPC(part, o2::track::PID::Pion), trackCuts.getNsigmaTOF(part, o2::track::PID::Pion), trackCuts.getNsigmaTPC(part, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(part, o2::track::PID::Kaon), trackCuts.getNsigmaTPC(part, o2::track::PID::Deuteron), trackCuts.getNsigmaTOF(part, o2::track::PID::Deuteron))) {
        continue;
      }

      trackHistoPartOneRec.fillQA<isMC, isDebug>(part);
      if (!part.has_fdMCParticle())
        continue;

      const auto mcParticle = part.fdMCParticle();
      registryPDG.fill(HIST("part1/PDGvspT"), part.pt(), mcParticle.pdgMCTruth());
      registryMCOrigin.fill(HIST("part1/hPt"), mcParticle.pt());
    }

    if (!ConfIsSame) {
      for (auto& part : grouppartsTwoMCRec) {

        if (part.partType() != ConfParticleTypePartTwo || !invMLambda(part.mLambda(), part.mAntiLambda())) {
          continue;
        }
        const auto& posChild = parts.iteratorAt(part.index() - 2);
        const auto& negChild = parts.iteratorAt(part.index() - 1);

        if (ConfPDGCodePartTwo > 0 && (!IsProtonNSigma(0, trackCuts.getNsigmaTPC(posChild, o2::track::PID::Proton), trackCuts.getNsigmaTOF(posChild, o2::track::PID::Proton)) || !IsPionNSigma(0, trackCuts.getNsigmaTPC(negChild, o2::track::PID::Pion), trackCuts.getNsigmaTOF(negChild, o2::track::PID::Pion)))) { // give momentum as 0 to only check TPC nSigma, not combined with TOF
          continue;
        }
        if (ConfPDGCodePartTwo < 0 && (!IsProtonNSigma(0, trackCuts.getNsigmaTPC(negChild, o2::track::PID::Proton), trackCuts.getNsigmaTOF(negChild, o2::track::PID::Proton)) || !IsPionNSigma(0, trackCuts.getNsigmaTPC(posChild, o2::track::PID::Pion), trackCuts.getNsigmaTOF(posChild, o2::track::PID::Pion)))) { // give momentum as 0 to only check TPC nSigma, not combined with TOF
          continue;
        }

        trackHistoV0TwoRec.fillQA<isMC, isDebug>(part);
        trackHistoV0TwoChildPosRec.fillQABase<isMC, isDebug>(posChild, HIST("posChildV0_2"));
        trackHistoV0TwoChildNegRec.fillQABase<isMC, isDebug>(negChild, HIST("negChildV0_2"));

        if (!posChild.has_fdMCParticle() || !negChild.has_fdMCParticle() || !part.has_fdMCParticle()) {
          continue;
        }
        const auto mcParticle = part.fdMCParticle();
        const auto mcPosChild = posChild.fdMCParticle();
        const auto mcNegChild = negChild.fdMCParticle();

        registryPDG.fill(HIST("part2/PDGvspT"), part.pt(), mcParticle.pdgMCTruth());
        registryPDG.fill(HIST("part2/dpositive/PDGvspT"), part.pt(), mcPosChild.pdgMCTruth());
        registryPDG.fill(HIST("part2/dnegative/PDGvspT"), part.pt(), mcNegChild.pdgMCTruth());
        registryMCOrigin.fill(HIST("part2/hPt"), mcParticle.pt());
      }
    }
  }

  /// process function for to call doMCRecTrackTrack with Data
  /// \param col subscribe to the collision table (Data)
  void processTrackTrack(FilteredFDCollision& col,
                         FemtoFullParticles&, aod::FdMCParticles const&)
  {
    fillCollision(col);
    // MCGen
    auto thegrouppartsOneMCGen = partsOneMCGen->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    auto thegrouppartsTwoMCGen = partsTwoGen->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    doMCGen<false>(thegrouppartsOneMCGen, thegrouppartsTwoMCGen);
    // MCRec
    auto thegroupPartsTrackOneRec = partsTrackOneMCReco->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    auto thegroupPartsTrackTwoRec = partsTrackTwoMCReco->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    if (ConfIsDebug) {
      doMCRecTrackTrack<false, true>(thegroupPartsTrackOneRec, thegroupPartsTrackTwoRec);
    } else {
      doMCRecTrackTrack<false, false>(thegroupPartsTrackOneRec, thegroupPartsTrackTwoRec);
    }
  }
  PROCESS_SWITCH(femtoUniverseEfficiencyBase, processTrackTrack, "Enable processing track-track efficiency task", true);

  /// process function for to call doMCRecTrackPhi with Data
  /// \param col subscribe to the collision table (Data)
  void processTrackPhi(FilteredFDCollision& col,
                       FemtoFullParticles&, aod::FdMCParticles const&)
  {
    fillCollision(col);
    // MCGen
    auto thegrouppartsOneMCGen = partsOneMCGen->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    auto thegrouppartsTwoMCGen = partsTwoGen->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    doMCGen<false>(thegrouppartsOneMCGen, thegrouppartsTwoMCGen);
    // MCRec
    auto thegroupPartsTrackOneRec = partsTrackOneMCReco->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    auto thegroupPartsTrackTwoRec = partsTrackTwoMCReco->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    if (ConfIsDebug) {
      doMCRecTrackPhi<false, true>(thegroupPartsTrackOneRec, thegroupPartsTrackTwoRec);
    } else {
      doMCRecTrackPhi<false, false>(thegroupPartsTrackOneRec, thegroupPartsTrackTwoRec);
    }
  }
  PROCESS_SWITCH(femtoUniverseEfficiencyBase, processTrackPhi, "Enable processing track-phi efficiency task", false);

  /// process function for to call doMCRecV0V0 with Data
  /// \param col subscribe to the collision table (Data)
  /// \param parts subscribe to the femtoUniverseParticleTable
  void processV0V0(FilteredFDCollision& col,
                   FemtoFullParticles& parts, aod::FdMCParticles const&)
  {
    fillCollision(col);
    // MCGen
    auto thegrouppartsOneMCGen = partsOneMCGen->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    auto thegrouppartsTwoMCGen = partsTwoGen->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    doMCGen<false>(thegrouppartsOneMCGen, thegrouppartsTwoMCGen);

    // MCRec
    auto thegroupPartsTrackOneRec = partsTrackOneMCReco->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    auto thegroupPartsTrackTwoRec = partsTrackTwoMCReco->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    if (ConfIsDebug) {
      doMCRecV0V0<false, true>(thegroupPartsTrackOneRec, thegroupPartsTrackTwoRec, parts);
    } else {
      doMCRecV0V0<false, false>(thegroupPartsTrackOneRec, thegroupPartsTrackTwoRec, parts);
    }
  }
  PROCESS_SWITCH(femtoUniverseEfficiencyBase, processV0V0, "Enable processing V0-V0 efficiency task", false);

  /// process function for to call doMCRecTrackV0 with Data
  /// \param col subscribe to the collision table (Data)
  /// \param parts subscribe to the femtoUniverseParticleTable
  void processTrackV0(FilteredFDCollision& col,
                      FemtoFullParticles& parts, aod::FdMCParticles const&)
  {
    fillCollision(col);
    // MCGen
    auto thegrouppartsOneMCGen = partsOneMCGen->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    auto thegrouppartsTwoMCGen = partsTwoGen->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    doMCGen<false>(thegrouppartsOneMCGen, thegrouppartsTwoMCGen);
    // MCRec
    auto thegroupPartsTrackOneRec = partsTrackOneMCReco->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    auto thegroupPartsTrackTwoRec = partsTrackTwoMCReco->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    if (ConfIsDebug) {
      doMCRecTrackV0<false, true>(thegroupPartsTrackOneRec, thegroupPartsTrackTwoRec, parts);
    } else {
      doMCRecTrackV0<false, false>(thegroupPartsTrackOneRec, thegroupPartsTrackTwoRec, parts);
    }
  }
  PROCESS_SWITCH(femtoUniverseEfficiencyBase, processTrackV0, "Enable processing track-V0 efficiency task", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<femtoUniverseEfficiencyBase>(cfgc),
  };
  return workflow;
}
