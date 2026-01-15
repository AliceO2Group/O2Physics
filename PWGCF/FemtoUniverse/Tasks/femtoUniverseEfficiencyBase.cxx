// Copyright 2019-2025 CERN and copyright holders of ALICE O2.
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

#include "PWGCF/FemtoUniverse/Core/FemtoUniverseEventHisto.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseParticleHisto.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseTrackSelection.h"

#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/runDataProcessing.h"

#include <vector>

using namespace o2;
using namespace o2::analysis::femto_universe;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

struct FemtoUniverseEfficiencyBase {
  SliceCache cache;
  using FemtoFullParticles = soa::Join<aod::FDParticles, aod::FDExtParticles, aod::FDMCLabels>;
  Preslice<FemtoFullParticles> perCol = aod::femtouniverseparticle::fdCollisionId;

  Configurable<bool> confIsDebug{"confIsDebug", true, "Enable debug histograms"};
  Configurable<bool> confIsMCGen{"confIsMCGen", false, "Enable QA histograms for MC Gen"};
  Configurable<bool> confIsMCReco{"confIsMCReco", false, "Enable QA histograms for MC Reco"};

  // Collisions
  Configurable<float> confZVertex{"confZVertex", 10.f, "Event sel: Maximum z-Vertex (cm)"};

  Filter collisionFilter = (nabs(aod::collision::posZ) < confZVertex);
  using FilteredFDCollisions = soa::Filtered<o2::aod::FdCollisions>;
  using FilteredFDCollision = FilteredFDCollisions::iterator;

  /// Particle selection part
  /// Configurables for both particles
  ConfigurableAxis confTempFitVarpTBins{"confTempFitVarpTBins", {20, 0.5, 4.05}, "Binning of the pT in the pT vs. TempFitVar plot"};
  ConfigurableAxis confTempFitVarPDGBins{"confTempFitVarPDGBins", {6000, -2300, 2300}, "Binning of the PDG code in the pT vs. TempFitVar plot"};
  ConfigurableAxis confTempFitVarCPABins{"confTempFitVarCPABins", {1000, 0.9, 1}, "Binning of the pointing angle cosinus in the pT vs. TempFitVar plot"};
  ConfigurableAxis confTempFitVarDCABins{"confTempFitVarDCABins", {1000, -5, 5}, "Binning of the PDG code in the pT vs. TempFitVar plot"};

  struct : o2::framework::ConfigurableGroup {
    Configurable<float> confEtaMax{"confEtaMax", 0.8f, "Higher limit for |Eta| (the same for both particles)"};
    Configurable<float> confMomProton{"confMomProton", 0.75, "Momentum threshold for proton identification using TOF"};
    Configurable<float> confMomPion{"confMomPion", 0.75, "Momentum threshold for pion identification using TOF"};
    Configurable<float> confNsigmaCombinedProton{"confNsigmaCombinedProton", 3.0, "TPC and TOF Proton Sigma (combined) for momentum > confMomProton"};
    Configurable<float> confNsigmaTPCProton{"confNsigmaTPCProton", 3.0, "TPC Proton Sigma for momentum < confMomProton"};
    Configurable<float> confNsigmaCombinedPion{"confNsigmaCombinedPion", 3.0, "TPC and TOF Pion Sigma (combined) for momentum > confMomPion"};
    Configurable<float> confNsigmaTPCPion{"confNsigmaTPCPion", 3.0, "TPC Pion Sigma for momentum < confMomPion"};
    Configurable<bool> confPDGCheckMCReco{"confPDGCheckMCReco", true, "Check PDG code of MC reco paricles"};

  } ConfBothTracks;

  /// Lambda cuts
  Configurable<float> confV0InvMassLowLimit{"confV0InvMassLowLimit", 1.10, "Lower limit of the V0 invariant mass"};
  Configurable<float> confV0InvMassUpLimit{"confV0InvMassUpLimit", 1.13, "Upper limit of the V0 invariant mass"};

  /// Kaon configurable
  Configurable<bool> isKaonRun2{"isKaonRun2", false, "Enable kaon selection used in Run2"}; // to check consistency with Run2 results
  Configurable<bool> isKaonLF{"isKaonLF", false, "Enable kaon selection used in LF group"}; // select kaons according to the selection in LF group
  struct : o2::framework::ConfigurableGroup {
    // Momentum thresholds for Run2 and Run3
    Configurable<float> confMomKaonRun2{"confMomKaonRun2", 0.4, "Momentum threshold for kaon identification using ToF (Run2)"};
    Configurable<float> confMomKaonRun3{"confMomKaonRun3", 0.3, "Momentum threshold for kaon identification using ToF (Run3)"};
    Configurable<float> confMomKaon045{"confMomKaon045", 0.45, "Momentum threshold for kaon identification pT = 0.45 GeV/c"};
    Configurable<float> confMomKaon055{"confMomKaon055", 0.55, "Momentum threshold for kaon identification pT = 0.55 GeV/c"};
    Configurable<float> confMomKaon08{"confMomKaon08", 0.8, "Momentum threshold for kaon identification pT = 0.8 GeV/c"};
    Configurable<float> confMomKaon15{"confMomKaon15", 1.5, "Momentum threshold for kaon identification pT = 1.5 GeV/c"};
    // n sigma cuts for Run 2
    Configurable<float> confKaonNsigmaTPCbelow04Run2{"confKaonNsigmaTPCbelow04Run2", 2.0, "Reject kaons with pT below 0.4 if TPC n sigma is above this value."};
    Configurable<float> confKaonNsigmaTPCfrom04to045Run2{"confKaonNsigmaTPCfrom04to045Run2", 1.0, "Reject kaons within pT from 0.4 to 0.45 if TPC n sigma is above this value."};
    Configurable<float> confKaonNsigmaTPCfrom045to08Run2{"confKaonNsigmaTPCfrom045to08Run2", 3.0, "Reject kaons within pT from 0.45 to 0.8 if TPC n sigma is above this value."};
    Configurable<float> confKaonNsigmaTOFfrom045to08Run2{"confKaonNsigmaTOFfrom045to08Run2", 2.0, "Reject kaons within pT from 0.45 to 0.8 if ToF n sigma is above this value."};
    Configurable<float> confKaonNsigmaTPCfrom08to15Run2{"confKaonNsigmaTPCfrom08to15Run2", 3.0, "Reject kaons within pT from 0.8 to 1.5 if TPC n sigma is above this value."};
    Configurable<float> confKaonNsigmaTOFfrom08to15Run2{"confKaonNsigmaTOFfrom08to15Run2", 1.5, "Reject kaons within pT from 0.8 to 1.5 if ToF n sigma is above this value."};
    // n sigma cuts for Run 3
    Configurable<float> confKaonNsigmaTPCfrom0to03{"confKaonNsigmaTPCfrom0to03", 3.0, "Reject kaons within pT from 0.0 to 0.3 if TPC n sigma is above this value."};
    Configurable<float> confKaonNsigmaTPCfrom03to045{"confKaonNsigmaTPCfrom03to045", 2.0, "Reject kaons within pT from 0.3 to 0.45 if TPC n sigma is above this value."};
    Configurable<float> confKaonNsigmaTPCfrom045to055{"confKaonNsigmaTPCfrom045to055", 1.0, "Reject kaons within pT from 0.45 to 0.55 if TPC n sigma is above this value."};
    Configurable<float> confKaonNsigmaTPCfrom055to15{"confKaonNsigmaTPCfrom055to15", 3.0, "Reject kaons within pT from 0.55 to 1.5 if TPC n sigma is above this value."};
    Configurable<float> confKaonNsigmaTOFfrom055to15{"confKaonNsigmaTOFfrom055to15", 3.0, "Reject kaons within pT from 0.55 to 1.5 if ToF n sigma is above this value."};
    Configurable<float> confKaonNsigmaTPCfrom15{"confKaonNsigmaTPCfrom15", 3.0, "Reject kaons with pT above 1.5 if TPC n sigma is above this value."};
    Configurable<float> confKaonNsigmaTOFfrom15{"confKaonNsigmaTOFfrom15", 2.0, "Reject kaons with pT above 1.5 if ToF n sigma is above this value.."};
    // n sigma cuts as in LF
    Configurable<float> confMomKaonLF{"confMomKaonLF", 0.5, "Momentum threshold for kaon identification using TOF (LF selection)"};
    Configurable<float> confNSigmaTPCKaonLF{"confNSigmaTPCKaonLF", 3.0, "TPC Kaon Sigma as in LF"};
    Configurable<float> confNSigmaCombKaonLF{"confNSigmaCombKaonLF", 3.0, "TPC and TOF Kaon Sigma (combined) as in LF"};
  } ConfKaonSelection;

  /// Deuteron configurables
  struct : o2::framework::ConfigurableGroup {
    Configurable<float> confNsigmaTPCDe{"confNsigmaTPCDe", 2.0f, "TPC Deuteron Sigma for momentum < confTOFpMinDe"};
    Configurable<float> confNsigmaTOFDe{"confNsigmaTOFDe", 2.0f, "TOF Deuteron Sigma"};
    Configurable<float> confTOFpMinDe{"confTOFpMinDe", 0.5f, "Min. momentum for deuterons for which TOF is required for PID"};
    Configurable<float> confPLowDe{"confPLowDe", 0.8f, "Lower limit for momentum for deuterons"};
    Configurable<float> confPHighDe{"confPHighDe", 1.8f, "Higher limit for momentum for deuterons"};
  } deuteronconfigs;

  /// Particle 1
  Configurable<int32_t> confPDGCodePartOne{"confPDGCodePartOne", 2212, "Particle 1 - PDG code"};
  Configurable<uint> confParticleTypePartOne{"confParticleTypePartOne", aod::femtouniverseparticle::ParticleType::kTrack, "Particle 1 - particle type: 0 - track, 2 - V0, 6 - phi"};
  Configurable<bool> confNoPDGPartOne{"confNoPDGPartOne", false, "0: selecting part one by PDG, 1: no PID selection"};
  Configurable<float> confPtLowPart1{"confPtLowPart1", 0.2, "Lower limit for Pt for the first particle"};
  Configurable<float> confPtHighPart1{"confPtHighPart1", 2.5, "Higher limit for Pt for the first particle"};
  Configurable<int> confChargePart1{"confChargePart1", 1, "Charge of the first particle"};

  /// Partition for particle 1
  Partition<FemtoFullParticles> partsOneMCGen = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kMCTruthTrack)) && (aod::femtouniverseparticle::pt < confPtHighPart1) && (aod::femtouniverseparticle::pt > confPtLowPart1) && (nabs(aod::femtouniverseparticle::eta) < ConfBothTracks.confEtaMax);

  Partition<FemtoFullParticles> partsTrackOneMCReco = (aod::femtouniverseparticle::pt < confPtHighPart1) && (aod::femtouniverseparticle::pt > confPtLowPart1) && (nabs(aod::femtouniverseparticle::eta) < ConfBothTracks.confEtaMax);

  /// Histogramming for particle 1
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kMCTruthTrack, 1> trackHistoPartOneGen;
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kTrack, 1> trackHistoPartOneRec;
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kV0, 1> trackHistoV0OneRec;
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kV0Child, 3> trackHistoV0OneChildPosRec;
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kV0Child, 4> trackHistoV0OneChildNegRec;

  /// Particle 2
  Configurable<bool> confIsSame{"confIsSame", false, "Pairs of the same particle"};
  Configurable<int32_t> confPDGCodePartTwo{"confPDGCodePartTwo", 333, "Particle 2 - PDG code"};
  Configurable<uint> confParticleTypePartTwo{"confParticleTypePartTwo", aod::femtouniverseparticle::ParticleType::kTrack, "Particle 2 - particle type:  0 - track, 2 - V0, 6 - phi"};
  Configurable<bool> confNoPDGPartTwo{"confNoPDGPartTwo", false, "0: selecting part two by PDG, 1: no PID selection"};
  Configurable<float> confPtLowPart2{"confPtLowPart2", 0.2, "Lower limit for Pt for the second particle"};
  Configurable<float> confPtHighPart2{"confPtHighPart2", 2.5, "Higher limit for Pt for the second particle"};
  Configurable<int> confChargePart2{"confChargePart2", 1, "Charge of the second particle"};

  /// Partition for particle 2
  Partition<FemtoFullParticles> partsTwoMCGen = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kMCTruthTrack)) && (aod::femtouniverseparticle::pt < confPtHighPart2) && (aod::femtouniverseparticle::pt > confPtLowPart2) && (nabs(aod::femtouniverseparticle::eta) < ConfBothTracks.confEtaMax);

  Partition<FemtoFullParticles> partsTrackTwoMCReco = (aod::femtouniverseparticle::pt < confPtHighPart2) && (aod::femtouniverseparticle::pt > confPtLowPart2) && (nabs(aod::femtouniverseparticle::eta) < ConfBothTracks.confEtaMax);

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
  HistogramRegistry registryCuts{"CutsPtHistos", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};
  HistogramRegistry registryMCOrigin{"MCOriginHistos", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};

  void init(InitContext&)
  {

    eventHisto.init(&qaRegistry);
    registryCuts.add("part1/cutsVspT", ";#it{p}_{T} (GeV/c) ;Cut no.", {HistType::kTH2F, {{500, 0, 5}, {7, 0, 7}}});
    trackHistoPartOneGen.init(&qaRegistry, confTempFitVarpTBins, confTempFitVarPDGBins, confIsMCGen, confPDGCodePartOne, false);
    trackHistoPartOneRec.init(&qaRegistry, confTempFitVarpTBins, confTempFitVarDCABins, confIsMCReco, confPDGCodePartOne, confIsDebug);
    registryMCOrigin.add("part1/hPt", " ;#it{p}_{T} (GeV/c); Entries", {HistType::kTH1F, {{240, 0, 6}}});
    registryPDG.add("part1/PDGvspT", "PDG;#it{p}_{T} (GeV/c); PDG", {HistType::kTH2F, {{500, 0, 5}, {16001, -8000.5, 8000.5}}});
    registryPDG.add("part1/PDGvspTall", "PDG;#it{p}_{T} (GeV/c); PDG", {HistType::kTH2F, {{500, 0, 5}, {16001, -8000.5, 8000.5}}});
    if (confParticleTypePartOne == uint8_t(aod::femtouniverseparticle::ParticleType::kV0)) {
      trackHistoV0OneRec.init(&qaRegistry, confTempFitVarpTBins, confTempFitVarCPABins, 0, confPDGCodePartOne, confIsDebug);
      trackHistoV0OneChildPosRec.init(&qaRegistry, confTempFitVarpTBins, confTempFitVarDCABins, 0, 0, confIsDebug, "posChildV0_1");
      trackHistoV0OneChildNegRec.init(&qaRegistry, confTempFitVarpTBins, confTempFitVarDCABins, 0, 0, confIsDebug, "negChildV0_1");
      registryPDG.add("part1/dpositive/PDGvspT", "PDG;#it{p}_{T} (GeV/c); PDG", {HistType::kTH2F, {{500, 0, 5}, {16001, -8000.5, 8000.5}}});
      registryPDG.add("part1/dnegative/PDGvspT", "PDG;#it{p}_{T} (GeV/c); PDG", {HistType::kTH2F, {{500, 0, 5}, {16001, -8000.5, 8000.5}}});
    }

    registryPDG.add("part2/PDGvspT", "PDG;#it{p}_{T} (GeV/c); PDG", {HistType::kTH2F, {{500, 0, 5}, {16001, -8000.5, 8000.5}}});
    registryPDG.add("part2/PDGvspTall", "PDG;#it{p}_{T} (GeV/c); PDG", {HistType::kTH2F, {{500, 0, 5}, {16001, -8000.5, 8000.5}}});
    if (!confIsSame) {
      registryCuts.add("part2/cutsVspT", ";#it{p}_{T} (GeV/c) ;Cut no.", {HistType::kTH2F, {{500, 0, 5}, {7, 0, 7}}});
      trackHistoPartTwoGen.init(&qaRegistry, confTempFitVarpTBins, confTempFitVarPDGBins, confIsMCGen, confPDGCodePartTwo, false);
      trackHistoPartTwoRec.init(&qaRegistry, confTempFitVarpTBins, confTempFitVarDCABins, confIsMCReco, confPDGCodePartTwo, confIsDebug);
      registryMCOrigin.add("part2/hPt", " ;#it{p}_{T} (GeV/c); Entries", {HistType::kTH1F, {{240, 0, 6}}});
      if (confParticleTypePartTwo == uint8_t(aod::femtouniverseparticle::ParticleType::kV0)) {
        trackHistoV0TwoRec.init(&qaRegistry, confTempFitVarpTBins, confTempFitVarCPABins, 0, confPDGCodePartTwo, confIsDebug);
        trackHistoV0TwoChildPosRec.init(&qaRegistry, confTempFitVarpTBins, confTempFitVarDCABins, 0, 0, confIsDebug, "posChildV0_2");
        trackHistoV0TwoChildNegRec.init(&qaRegistry, confTempFitVarpTBins, confTempFitVarDCABins, 0, 0, confIsDebug, "negChildV0_2");
        registryPDG.add("part2/dpositive/PDGvspT", "PDG;#it{p}_{T} (GeV/c); PDG", {HistType::kTH2F, {{500, 0, 5}, {16001, -8000.5, 8000.5}}});
        registryPDG.add("part2/dnegative/PDGvspT", "PDG;#it{p}_{T} (GeV/c); PDG", {HistType::kTH2F, {{500, 0, 5}, {16001, -8000.5, 8000.5}}});
      }
    }
  }

  bool isProtonNSigma(float mom, float nsigmaTPCPr, float nsigmaTOFPr) // previous version from: https://github.com/alisw/AliPhysics/blob/master/PWGCF/FEMTOSCOPY/AliFemtoUser/AliFemtoMJTrackCut.cxx
  {
    if (mom < ConfBothTracks.confMomProton) {
      if (std::abs(nsigmaTPCPr) < ConfBothTracks.confNsigmaTPCProton) {
        return true;
      } else {
        return false;
      }
    } else {
      if (std::hypot(nsigmaTOFPr, nsigmaTPCPr) < ConfBothTracks.confNsigmaCombinedProton) {
        return true;
      } else {
        return false;
      }
    }
    return false;
  }

  bool isKaonNSigma(float mom, float nsigmaTPCK, float nsigmaTOFK)
  {
    if (isKaonRun2 == true) {
      if (mom < ConfKaonSelection.confMomKaonRun2) { // < 0.4 GeV/c
        return std::abs(nsigmaTPCK) < ConfKaonSelection.confKaonNsigmaTPCbelow04Run2;
      } else if (mom > ConfKaonSelection.confMomKaonRun2 && mom < ConfKaonSelection.confMomKaon045) { // 0.4 - 0.45
        return std::abs(nsigmaTPCK) < ConfKaonSelection.confKaonNsigmaTPCfrom04to045Run2;
      } else if (mom > ConfKaonSelection.confMomKaon045 && mom < ConfKaonSelection.confMomKaon08) { // 0.45 - 0.8
        return (std::abs(nsigmaTPCK) < ConfKaonSelection.confKaonNsigmaTPCfrom045to08Run2 && std::abs(nsigmaTOFK) < ConfKaonSelection.confKaonNsigmaTOFfrom045to08Run2);
      } else if (mom > ConfKaonSelection.confMomKaon08 && mom < ConfKaonSelection.confMomKaon15) { // 0.8 - 1.5
        return (std::abs(nsigmaTPCK) < ConfKaonSelection.confKaonNsigmaTPCfrom08to15Run2 && std::abs(nsigmaTOFK) < ConfKaonSelection.confKaonNsigmaTOFfrom08to15Run2);
      } else {
        return false;
      }
    } else {
      if (mom < ConfKaonSelection.confMomKaonRun3) { // 0.0-0.3
        if (std::abs(nsigmaTPCK) < ConfKaonSelection.confKaonNsigmaTPCfrom0to03) {
          return true;
        } else {
          return false;
        }
      } else if (mom < ConfKaonSelection.confMomKaon045) { // 0.30 - 0.45
        if (std::abs(nsigmaTPCK) < ConfKaonSelection.confKaonNsigmaTPCfrom03to045) {
          return true;
        } else {
          return false;
        }
      } else if (mom < ConfKaonSelection.confMomKaon055) { // 0.45-0.55
        if (std::abs(nsigmaTPCK) < ConfKaonSelection.confKaonNsigmaTPCfrom045to055) {
          return true;
        } else {
          return false;
        }
      } else if (mom < ConfKaonSelection.confMomKaon15) { // 0.55-1.5 (now we use TPC and TOF)
        if ((std::abs(nsigmaTOFK) < ConfKaonSelection.confKaonNsigmaTOFfrom055to15) && (std::abs(nsigmaTPCK) < ConfKaonSelection.confKaonNsigmaTPCfrom055to15)) {
          {
            return true;
          }
        } else {
          return false;
        }
      } else if (mom > ConfKaonSelection.confMomKaon15) { // > 1.5 GeV/c
        if ((std::abs(nsigmaTOFK) < ConfKaonSelection.confKaonNsigmaTOFfrom15) && (std::abs(nsigmaTPCK) < ConfKaonSelection.confKaonNsigmaTPCfrom15)) {
          return true;
        } else {
          return false;
        }
      } else {
        return false;
      }
    }
  }

  bool isKaonNSigmaLF(float mom, float nsigmaTPCK, float nsigmaTOFK)
  {
    if (mom < ConfKaonSelection.confMomKaonLF) {
      if (std::abs(nsigmaTPCK) < ConfKaonSelection.confNSigmaTPCKaonLF) {
        return true;
      } else {
        return false;
      }
    } else if (mom >= ConfKaonSelection.confMomKaonLF) {
      if (std::sqrt(nsigmaTPCK * nsigmaTPCK + nsigmaTOFK * nsigmaTOFK) < ConfKaonSelection.confNSigmaCombKaonLF) {
        return true;
      } else {
        return false;
      }
    } else {
      return false;
    }
  }

  bool isPionNSigma(float mom, float nsigmaTPCPi, float nsigmaTOFPi)
  {
    if (mom < ConfBothTracks.confMomPion) {
      if (std::abs(nsigmaTPCPi) < ConfBothTracks.confNsigmaTPCPion) {
        return true;
      } else {
        return false;
      }
    } else {
      if (std::hypot(nsigmaTOFPi, nsigmaTPCPi) < ConfBothTracks.confNsigmaCombinedPion) {
        return true;
      } else {
        return false;
      }
    }
    return false;
  }

  bool isDeuteronNSigma(float mom, float nsigmaTPCDe, float nsigmaTOFDe)
  {
    if (mom > deuteronconfigs.confPLowDe && mom < deuteronconfigs.confPHighDe) {
      if (mom < deuteronconfigs.confTOFpMinDe) {
        return (std::abs(nsigmaTPCDe) < deuteronconfigs.confNsigmaTPCDe);
      } else {
        return (std::abs(nsigmaTOFDe) < deuteronconfigs.confNsigmaTOFDe && (std::abs(nsigmaTPCDe) < deuteronconfigs.confNsigmaTPCDe));
      }
    } else {
      return false;
    }
  }

  bool isParticleNSigma(int pdgCode, float mom, float nsigmaTPCPr, float nsigmaTOFPr, float nsigmaTPCPi, float nsigmaTOFPi, float nsigmaTPCK, float nsigmaTOFK, float nsigmaTPCDe, float nsigmaTOFDe)
  {
    switch (pdgCode) {
      case 2212:  // Proton
      case -2212: // anty Proton
        return isProtonNSigma(mom, nsigmaTPCPr, nsigmaTOFPr);
        break;
      case 211:  // Pion
      case -211: // Pion-
        return isPionNSigma(mom, nsigmaTPCPi, nsigmaTOFPi);
        break;
      case 321:  // Kaon+
      case -321: // Kaon-
        if (isKaonLF) {
          return isKaonNSigmaLF(mom, nsigmaTPCK, nsigmaTOFK);
        } else {
          return isKaonNSigma(mom, nsigmaTPCK, nsigmaTOFK);
        }
        break;
      case 1000010020:  // Deuteron
      case -1000010020: // Antideuteron
        return isDeuteronNSigma(mom, nsigmaTPCDe, nsigmaTOFDe);
        break;
      default:
        return false;
    }
  }

  bool invMLambda(float invMassLambda, float invMassAntiLambda)
  {
    if ((invMassLambda < confV0InvMassLowLimit || invMassLambda > confV0InvMassUpLimit) && (invMassAntiLambda < confV0InvMassLowLimit || invMassAntiLambda > confV0InvMassUpLimit)) {
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
    for (const auto& part : grouppartsOneMCGen) {
      if (!confNoPDGPartOne && part.tempFitVar() != confPDGCodePartOne) {
        continue;
      }
      trackHistoPartOneGen.fillQA<isMC, false>(part);
    }

    if (!confIsSame) {
      for (const auto& part : grouppartsTwoMCGen) {
        if (!confNoPDGPartTwo && part.tempFitVar() != confPDGCodePartTwo) {
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
    for (const auto& part : grouppartsOneMCRec) {

      // only partition
      registryCuts.fill(HIST("part1/cutsVspT"), part.pt(), 0);

      if (part.partType() != confParticleTypePartOne) {
        continue;
      }
      registryCuts.fill(HIST("part1/cutsVspT"), part.pt(), 1);

      if (part.sign() != confChargePart1) {
        continue;
      }
      registryCuts.fill(HIST("part1/cutsVspT"), part.pt(), 2);

      if (!isParticleNSigma(confPDGCodePartOne, part.p(), trackCuts.getNsigmaTPC(part, o2::track::PID::Proton), trackCuts.getNsigmaTOF(part, o2::track::PID::Proton), trackCuts.getNsigmaTPC(part, o2::track::PID::Pion), trackCuts.getNsigmaTOF(part, o2::track::PID::Pion), trackCuts.getNsigmaTPC(part, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(part, o2::track::PID::Kaon), trackCuts.getNsigmaTPC(part, o2::track::PID::Deuteron), trackCuts.getNsigmaTOF(part, o2::track::PID::Deuteron))) {
        continue;
      }
      registryCuts.fill(HIST("part1/cutsVspT"), part.pt(), 3);

      if (!part.has_fdMCParticle()) {
        continue;
      }
      registryCuts.fill(HIST("part1/cutsVspT"), part.pt(), 4);

      const auto mcParticle = part.fdMCParticle();

      registryPDG.fill(HIST("part1/PDGvspTall"), part.pt(), mcParticle.pdgMCTruth());
      trackHistoPartOneRec.fillQA<isMC, isDebug>(part);

      if (!(mcParticle.partOriginMCTruth() == aod::femtouniverse_mc_particle::ParticleOriginMCTruth::kPrimary)) {
        continue;
      }
      registryCuts.fill(HIST("part1/cutsVspT"), part.pt(), 5);

      if (ConfBothTracks.confPDGCheckMCReco && !(std::abs(mcParticle.pdgMCTruth()) == std::abs(confPDGCodePartOne))) {
        continue;
      }
      registryCuts.fill(HIST("part1/cutsVspT"), part.pt(), 6);

      registryPDG.fill(HIST("part1/PDGvspT"), part.pt(), mcParticle.pdgMCTruth());
      registryMCOrigin.fill(HIST("part1/hPt"), mcParticle.pt());
    }

    if (!confIsSame) {
      for (const auto& part : grouppartsTwoMCRec) {
        registryCuts.fill(HIST("part2/cutsVspT"), part.pt(), 0);

        if (part.partType() != confParticleTypePartOne) {
          continue;
        }
        registryCuts.fill(HIST("part2/cutsVspT"), part.pt(), 1);

        if (part.sign() != confChargePart1) {
          continue;
        }
        registryCuts.fill(HIST("part2/cutsVspT"), part.pt(), 2);

        if (!isParticleNSigma(confPDGCodePartOne, part.p(), trackCuts.getNsigmaTPC(part, o2::track::PID::Proton), trackCuts.getNsigmaTOF(part, o2::track::PID::Proton), trackCuts.getNsigmaTPC(part, o2::track::PID::Pion), trackCuts.getNsigmaTOF(part, o2::track::PID::Pion), trackCuts.getNsigmaTPC(part, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(part, o2::track::PID::Kaon), trackCuts.getNsigmaTPC(part, o2::track::PID::Deuteron), trackCuts.getNsigmaTOF(part, o2::track::PID::Deuteron))) {
          continue;
        }
        registryCuts.fill(HIST("part2/cutsVspT"), part.pt(), 3);

        if (!part.has_fdMCParticle()) {
          continue;
        }
        registryCuts.fill(HIST("part2/cutsVspT"), part.pt(), 4);
        const auto mcParticle = part.fdMCParticle();

        registryPDG.fill(HIST("part2/PDGvspTall"), part.pt(), mcParticle.pdgMCTruth());
        trackHistoPartTwoRec.fillQA<isMC, isDebug>(part);

        if (!(mcParticle.partOriginMCTruth() == aod::femtouniverse_mc_particle::ParticleOriginMCTruth::kPrimary)) {
          continue;
        }
        registryCuts.fill(HIST("part2/cutsVspT"), part.pt(), 5);

        if (ConfBothTracks.confPDGCheckMCReco && !(std::abs(mcParticle.pdgMCTruth()) == std::abs(confPDGCodePartTwo))) {
          continue;
        }
        registryCuts.fill(HIST("part2/cutsVspT"), part.pt(), 6);

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

    for (const auto& part : grouppartsOneMCRec) {
      if (part.partType() != confParticleTypePartOne || part.sign() != confChargePart1 || !isParticleNSigma(confPDGCodePartOne, part.p(), trackCuts.getNsigmaTPC(part, o2::track::PID::Proton), trackCuts.getNsigmaTOF(part, o2::track::PID::Proton), trackCuts.getNsigmaTPC(part, o2::track::PID::Pion), trackCuts.getNsigmaTOF(part, o2::track::PID::Pion), trackCuts.getNsigmaTPC(part, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(part, o2::track::PID::Kaon), trackCuts.getNsigmaTPC(part, o2::track::PID::Deuteron), trackCuts.getNsigmaTOF(part, o2::track::PID::Deuteron))) {
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

    if (!confIsSame) {
      for (const auto& part : grouppartsTwoMCRec) {
        if (part.partType() != confParticleTypePartTwo || part.sign() != confChargePart2) {
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
    for (const auto& part : grouppartsOneMCRec) {

      if (part.partType() != confParticleTypePartOne || !invMLambda(part.mLambda(), part.mAntiLambda())) {
        continue;
      }

      const auto& posChild = parts.iteratorAt(part.index() - 2);
      const auto& negChild = parts.iteratorAt(part.index() - 1);

      if (confPDGCodePartOne > 0 && (!isProtonNSigma(0, trackCuts.getNsigmaTPC(posChild, o2::track::PID::Proton), trackCuts.getNsigmaTOF(posChild, o2::track::PID::Proton)) || !isPionNSigma(0, trackCuts.getNsigmaTPC(negChild, o2::track::PID::Pion), trackCuts.getNsigmaTOF(negChild, o2::track::PID::Pion)))) { // give momentum as 0 to only check TPC nSigma, not combined with TOF
        continue;
      }
      if (confPDGCodePartOne < 0 && (!isProtonNSigma(0, trackCuts.getNsigmaTPC(negChild, o2::track::PID::Proton), trackCuts.getNsigmaTOF(negChild, o2::track::PID::Proton)) || !isPionNSigma(0, trackCuts.getNsigmaTPC(posChild, o2::track::PID::Pion), trackCuts.getNsigmaTOF(posChild, o2::track::PID::Pion)))) { // give momentum as 0 to only check TPC nSigma, not combined with TOF
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

    if (!confIsSame) {
      for (const auto& part : grouppartsTwoMCRec) {

        if (part.partType() != confParticleTypePartTwo || !invMLambda(part.mLambda(), part.mAntiLambda())) {
          continue;
        }

        const auto& posChild = parts.iteratorAt(part.index() - 2);
        const auto& negChild = parts.iteratorAt(part.index() - 1);

        if (confPDGCodePartTwo > 0 && (!isProtonNSigma(0, trackCuts.getNsigmaTPC(posChild, o2::track::PID::Proton), trackCuts.getNsigmaTOF(posChild, o2::track::PID::Proton)) || !isPionNSigma(0, trackCuts.getNsigmaTPC(negChild, o2::track::PID::Pion), trackCuts.getNsigmaTOF(negChild, o2::track::PID::Pion)))) { // give momentum as 0 to only check TPC nSigma, not combined with TOF
          continue;
        }
        if (confPDGCodePartTwo < 0 && (!isProtonNSigma(0, trackCuts.getNsigmaTPC(negChild, o2::track::PID::Proton), trackCuts.getNsigmaTOF(negChild, o2::track::PID::Proton)) || !isPionNSigma(0, trackCuts.getNsigmaTPC(posChild, o2::track::PID::Pion), trackCuts.getNsigmaTOF(posChild, o2::track::PID::Pion)))) { // give momentum as 0 to only check TPC nSigma, not combined with TOF
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
    for (const auto& part : grouppartsOneMCRec) {
      if (part.partType() != confParticleTypePartOne || part.sign() != confChargePart1 || !isParticleNSigma(confPDGCodePartOne, part.p(), trackCuts.getNsigmaTPC(part, o2::track::PID::Proton), trackCuts.getNsigmaTOF(part, o2::track::PID::Proton), trackCuts.getNsigmaTPC(part, o2::track::PID::Pion), trackCuts.getNsigmaTOF(part, o2::track::PID::Pion), trackCuts.getNsigmaTPC(part, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(part, o2::track::PID::Kaon), trackCuts.getNsigmaTPC(part, o2::track::PID::Deuteron), trackCuts.getNsigmaTOF(part, o2::track::PID::Deuteron))) {
        continue;
      }

      trackHistoPartOneRec.fillQA<isMC, isDebug>(part);
      if (!part.has_fdMCParticle())
        continue;

      const auto mcParticle = part.fdMCParticle();
      registryPDG.fill(HIST("part1/PDGvspT"), part.pt(), mcParticle.pdgMCTruth());
      registryMCOrigin.fill(HIST("part1/hPt"), mcParticle.pt());
    }

    if (!confIsSame) {
      for (const auto& part : grouppartsTwoMCRec) {

        if (part.partType() != confParticleTypePartTwo || !invMLambda(part.mLambda(), part.mAntiLambda())) {
          continue;
        }
        const auto& posChild = parts.iteratorAt(part.index() - 2);
        const auto& negChild = parts.iteratorAt(part.index() - 1);

        if (confPDGCodePartTwo > 0 && (!isProtonNSigma(0, trackCuts.getNsigmaTPC(posChild, o2::track::PID::Proton), trackCuts.getNsigmaTOF(posChild, o2::track::PID::Proton)) || !isPionNSigma(0, trackCuts.getNsigmaTPC(negChild, o2::track::PID::Pion), trackCuts.getNsigmaTOF(negChild, o2::track::PID::Pion)))) { // give momentum as 0 to only check TPC nSigma, not combined with TOF
          continue;
        }
        if (confPDGCodePartTwo < 0 && (!isProtonNSigma(0, trackCuts.getNsigmaTPC(negChild, o2::track::PID::Proton), trackCuts.getNsigmaTOF(negChild, o2::track::PID::Proton)) || !isPionNSigma(0, trackCuts.getNsigmaTPC(posChild, o2::track::PID::Pion), trackCuts.getNsigmaTOF(posChild, o2::track::PID::Pion)))) { // give momentum as 0 to only check TPC nSigma, not combined with TOF
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
  void processTrackTrack(FilteredFDCollision const& col,
                         FemtoFullParticles const&, aod::FdMCParticles const&)
  {
    fillCollision(col);
    // MCGen
    auto thegrouppartsOneMCGen = partsOneMCGen->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    auto thegrouppartsTwoMCGen = partsTwoMCGen->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    if (confIsMCGen) {
      doMCGen<true>(thegrouppartsOneMCGen, thegrouppartsTwoMCGen);
    } else {
      doMCGen<false>(thegrouppartsOneMCGen, thegrouppartsTwoMCGen);
    }
    // MCRec
    auto thegroupPartsTrackOneRec = partsTrackOneMCReco->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    auto thegroupPartsTrackTwoRec = partsTrackTwoMCReco->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    if (confIsDebug) {
      if (confIsMCGen) {
        doMCRecTrackTrack<true, true>(thegroupPartsTrackOneRec, thegroupPartsTrackTwoRec);
      } else {
        doMCRecTrackTrack<false, true>(thegroupPartsTrackOneRec, thegroupPartsTrackTwoRec);
      }
    } else {
      doMCRecTrackTrack<false, false>(thegroupPartsTrackOneRec, thegroupPartsTrackTwoRec);
    }
  }
  PROCESS_SWITCH(FemtoUniverseEfficiencyBase, processTrackTrack, "Enable processing track-track efficiency task", true);

  /// process function for to call doMCRecTrackPhi with Data
  /// \param col subscribe to the collision table (Data)
  void processTrackPhi(FilteredFDCollision const& col,
                       FemtoFullParticles const&, aod::FdMCParticles const&)
  {
    fillCollision(col);
    // MCGen
    auto thegrouppartsOneMCGen = partsOneMCGen->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    auto thegrouppartsTwoMCGen = partsTwoMCGen->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    doMCGen<false>(thegrouppartsOneMCGen, thegrouppartsTwoMCGen);
    // MCRec
    auto thegroupPartsTrackOneRec = partsTrackOneMCReco->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    auto thegroupPartsTrackTwoRec = partsTrackTwoMCReco->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    if (confIsDebug) {
      doMCRecTrackPhi<false, true>(thegroupPartsTrackOneRec, thegroupPartsTrackTwoRec);
    } else {
      doMCRecTrackPhi<false, false>(thegroupPartsTrackOneRec, thegroupPartsTrackTwoRec);
    }
  }
  PROCESS_SWITCH(FemtoUniverseEfficiencyBase, processTrackPhi, "Enable processing track-phi efficiency task", false);

  /// process function for to call doMCRecV0V0 with Data
  /// \param col subscribe to the collision table (Data)
  /// \param parts subscribe to the femtoUniverseParticleTable
  void processV0V0(FilteredFDCollision const& col,
                   FemtoFullParticles const& parts, aod::FdMCParticles const&)
  {
    fillCollision(col);
    // MCGen
    auto thegrouppartsOneMCGen = partsOneMCGen->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    auto thegrouppartsTwoMCGen = partsTwoMCGen->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    doMCGen<false>(thegrouppartsOneMCGen, thegrouppartsTwoMCGen);

    // MCRec
    auto thegroupPartsTrackOneRec = partsTrackOneMCReco->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    auto thegroupPartsTrackTwoRec = partsTrackTwoMCReco->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    if (confIsDebug) {
      doMCRecV0V0<false, true>(thegroupPartsTrackOneRec, thegroupPartsTrackTwoRec, parts);
    } else {
      doMCRecV0V0<false, false>(thegroupPartsTrackOneRec, thegroupPartsTrackTwoRec, parts);
    }
  }
  PROCESS_SWITCH(FemtoUniverseEfficiencyBase, processV0V0, "Enable processing V0-V0 efficiency task", false);

  /// process function for to call doMCRecTrackV0 with Data
  /// \param col subscribe to the collision table (Data)
  /// \param parts subscribe to the femtoUniverseParticleTable
  void processTrackV0(FilteredFDCollision const& col,
                      FemtoFullParticles const& parts, aod::FdMCParticles const&)
  {
    fillCollision(col);
    // MCGen
    auto thegrouppartsOneMCGen = partsOneMCGen->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    auto thegrouppartsTwoMCGen = partsTwoMCGen->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    doMCGen<false>(thegrouppartsOneMCGen, thegrouppartsTwoMCGen);
    // MCRec
    auto thegroupPartsTrackOneRec = partsTrackOneMCReco->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    auto thegroupPartsTrackTwoRec = partsTrackTwoMCReco->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    if (confIsDebug) {
      doMCRecTrackV0<false, true>(thegroupPartsTrackOneRec, thegroupPartsTrackTwoRec, parts);
    } else {
      doMCRecTrackV0<false, false>(thegroupPartsTrackOneRec, thegroupPartsTrackTwoRec, parts);
    }
  }
  PROCESS_SWITCH(FemtoUniverseEfficiencyBase, processTrackV0, "Enable processing track-V0 efficiency task", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<FemtoUniverseEfficiencyBase>(cfgc),
  };
  return workflow;
}
