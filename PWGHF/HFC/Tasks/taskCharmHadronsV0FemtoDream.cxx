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

/// \file taskCharmHadronsTrackFemtoDream.cxx
/// \brief Tasks that reads the V0 and CharmHadrons tables used for the pairing and builds pairs
/// \author Biao Zhang, Heidelberg University, biao.zhang@cern.ch

#include "PWGCF/DataModel/FemtoDerived.h"
#include "PWGCF/FemtoDream/Core/femtoDreamContainer.h"
#include "PWGCF/FemtoDream/Core/femtoDreamDetaDphiStar.h"
#include "PWGCF/FemtoDream/Core/femtoDreamEventHisto.h"
#include "PWGCF/FemtoDream/Core/femtoDreamMath.h"
#include "PWGCF/FemtoDream/Core/femtoDreamPairCleaner.h"
#include "PWGCF/FemtoDream/Core/femtoDreamParticleHisto.h"
#include "PWGCF/FemtoDream/Core/femtoDreamUtils.h"
#include "PWGHF/Core/DecayChannels.h"

#include "Common/Core/RecoDecay.h"

#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoAHelpers.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/BinningPolicy.h>
#include <Framework/Configurable.h>
#include <Framework/Expressions.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <TPDGCode.h>

#include <array>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::aod;
using namespace o2::soa;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::analysis::femtoDream;
using namespace o2::hf_decay;
using namespace o2::hf_decay::hf_cand_3prong;
using namespace o2::constants::physics;

struct HfTaskCharmHadronsV0FemtoDream {

  enum TrackCharge {
    PositiveCharge = 1,
    NegativeCharge = -1
  };

  enum PairSign {
    PairNotDefined = 0,
    LikeSignPair = 1,
    UnLikeSignPair = 2
  };
  // decay channels
  enum DecayChannel { DplusToPiKPi = 0,
                      LcToPKPi,
                      D0ToPiK,
                      DstarToD0Pi
  };

  constexpr static int OriginRecPrompt = 1;
  constexpr static int OriginRecFD = 2;
  constexpr static int CutBitChargePositive = 2;

  Produces<o2::aod::FDHfCharmTrkPairs> rowFemtoResultPairs;
  Produces<o2::aod::FDHfCharm3Prong> rowFemtoResultCharm3Prong;
  Produces<o2::aod::FDHfCharm2Prong> rowFemtoResultCharm2Prong;
  Produces<o2::aod::FDHfCharmDstar> rowFemtoResultCharmDstar;
  Produces<o2::aod::FDHfV0> rowFemtoResultV0;
  Produces<o2::aod::FDHfColl> rowFemtoResultColl;

  Configurable<int> confTempFitVarMomentum{"confTempFitVarMomentum", 0, "Momentum used for binning: 0 -> pt; 1 -> preco; 2 -> ptpc"};

  /// Particle 2 (Charm Hadrons)
  struct : ConfigurableGroup {
    Configurable<float> charmHadBkgBDTmax{"charmHadBkgBDTmax", 1., "Maximum background bdt score for Charm Hadron (particle 2)"};
    Configurable<int> charmHadCandSel{"charmHadCandSel", 1, "candidate selection for charm hadron"};
    Configurable<int> charmHadMcSel{"charmHadMcSel", DecayChannelMain::LcToPKPi, "charm hadron selection for mc, DplusToPiKPi = 1, LcToPKPi = 17"};
    Configurable<float> charmHadFdBDTmin{"charmHadFdBDTmin", 0., "Minimum feed-down bdt score Charm Hadron (particle 2)"};
    Configurable<float> charmHadFdBDTmax{"charmHadFdBDTmax", 1., "Maximum feed-down bdt score Charm Hadron (particle 2)"};
    Configurable<float> charmHadMaxInvMass{"charmHadMaxInvMass", 2.45, "Maximum invariant mass of Charm Hadron (particle 2)"};
    Configurable<float> charmHadMinInvMass{"charmHadMinInvMass", 2.15, "Minimum invariant mass of Charm Hadron (particle 2)"};
    Configurable<float> charmHadMinPt{"charmHadMinPt", 0., "Minimum pT of Charm Hadron (particle 2)"};
    Configurable<float> charmHadMaxPt{"charmHadMaxPt", 999., "Maximum pT of Charm Hadron (particle 2)"};
    Configurable<int> charmHadPDGCode{"charmHadPDGCode", 4122, "PDG code of particle 2 Charm Hadron"};
    Configurable<float> charmHadPromptBDTmin{"charmHadPromptBDTmin", 0., "Minimum prompt bdt score Charm Hadron (particle 2)"};
    Configurable<float> charmHadPromptBDTmax{"charmHadPromptBDTmax", 1., "Maximum prompt bdt score Charm Hadron (particle 2)"};
  } charmSel;
  /// General options
  Configurable<float> cprDeltaEtaMax{"cprDeltaEtaMax", 0.01, "Max. Delta Eta for Close Pair Rejection"};
  Configurable<float> cprDeltaPhiMax{"cprDeltaPhiMax", 0.01, "Max. Delta Phi for Close Pair Rejection"};
  Configurable<bool> cprPlotPerRadii{"cprPlotPerRadii", false, "Plot CPR per radii"};
  Configurable<bool> extendedPlots{"extendedPlots", false, "Enable additional three dimensional histogramms. High memory consumption. Use for debugging"};
  Configurable<float> highkstarCut{"highkstarCut", 100000., "Set a cut for high k*, above which the pairs are rejected"};
  Configurable<bool> isMc{"isMc", false, "Set true in the case of a MonteCarlo Run"};
  Configurable<bool> smearingByOrigin{"smearingByOrigin", false, "Obtain the smearing matrix differential in the MC origin of particle 1 and particle 2. High memory consumption. Use with care!"};
  Configurable<bool> use4D{"use4D", false, "Enable four dimensional histogramms (to be used only for analysis with high statistics): k* vs multiplicity vs multiplicity percentil vs mT"};
  Configurable<bool> useCPR{"useCPR", false, "Close Pair Rejection"};
  Configurable<bool> fillTableWithCharm{"fillTableWithCharm", true, "Write charm/tracks/collision table only if >=1 charm hadron in this collision"};

  // Mixing configurables
  struct : ConfigurableGroup {
    Configurable<bool> doMixEvent{"doMixEvent", false, "choose do mix-event online"};
    Configurable<int> mixingBinPolicy{"mixingBinPolicy", 0, "Binning policy for mixing - 0: multiplicity, 1: multipliciy percentile, 2: both"};
    Configurable<int> mixingDepth{"mixingDepth", 5, "Number of events for mixing"};
  } mixSetting;
  /// Event selection
  struct : ConfigurableGroup {
    std::string prefix = "eventSel";
    Configurable<int> multMin{"multMin", 0, "Minimum Multiplicity (MultNtr)"};
    Configurable<int> multMax{"multMax", 99999, "Maximum Multiplicity (MultNtr)"};
    Configurable<float> multPercentileMin{"multPercentileMin", 0, "Maximum Multiplicity Percentile"};
    Configurable<float> multPercentileMax{"multPercentileMax", 100, "Minimum Multiplicity Percentile"};
  } eventSel;

  /// particle 1 (V0), Λ or Ks0
  struct : ConfigurableGroup {
    std::string prefix = std::string("v0Sel");
    Configurable<int> pdgCodeV0{"pdgCodeV0", 310, "PDG code of V0 (310: K0S, 3122: Lambda)"};
    Configurable<femtodreamparticle::cutContainerType> cutBit{"cutBit", 7518, "Selection bit for particle 1 (V0)"};
    Configurable<femtodreamparticle::cutContainerType> childPosCutBit{"childPosCutBit", 210, "Selection bit for positive child of V01"};
    Configurable<femtodreamparticle::cutContainerType> childPosTPCBit{"childPosTPCBit", 64, "PID TPC bit for positive child of V01"};
    Configurable<femtodreamparticle::cutContainerType> childNegCutBit{"childNegCutBit", 209, "Selection bit for negative child of V01"};
    Configurable<femtodreamparticle::cutContainerType> childNegTPCBit{"childNegTPCBit", 256, "PID TPC bit for negative child of V01"};
    Configurable<float> invMassV0Min{"invMassV0Min", 0.45, "Minimum invariant mass of Partricle 1 (particle) (V0)"};
    Configurable<float> invMassV0Max{"invMassV0Max", 0.55, "Maximum invariant mass of Partricle 1 (particle) (V0)"};
    Configurable<float> invMassAntiV0Min{"invMassAntiV0Min", 0.45, "Minimum invariant mass of Partricle 1 (antiparticle) (V0)"};
    Configurable<float> invMassAntiV0Max{"invMassAntiV0Max", 0.55, "Maximum invariant mass of Partricle 1 (antiparticle) (V0)"};
    Configurable<float> ptV0Min{"ptV0Min", 0., "Minimum pT of Partricle 1 (V0)"};
    Configurable<float> ptV0Max{"ptV0Max", 999., "Maximum pT of Partricle 1 (V0)"};
    Configurable<float> etaV0Min{"etaV0Min", -10., "Minimum eta of Partricle 1 (V0)"};
    Configurable<float> etaV0Max{"etaV0Max", 10., "Maximum eta of Partricle 1 (V0)"};
  } v0Sel;

  SliceCache cache;

  using FilteredCharmCand3Prongs = soa::Filtered<aod::FDHfCand3Prong>;
  using FilteredCharmCand3Prong = FilteredCharmCand3Prongs::iterator;

  using FilteredCharmCand2Prongs = soa::Filtered<aod::FDHfCand2Prong>;
  using FilteredCharmCand2Prong = FilteredCharmCand2Prongs::iterator;

  using FilteredCharmCandDstars = soa::Filtered<aod::FDHfCandDstar>;
  using FilteredCharmCandDstar = FilteredCharmCandDstars::iterator;

  using FilteredCharmMcCand3Prongs = soa::Filtered<soa::Join<aod::FDHfCand3Prong, aod::FDHfCandMC>>;
  using FilteredCharmMcCand3Prong = FilteredCharmMcCand3Prongs::iterator;

  using FilteredCharmMcCand2Prongs = soa::Filtered<soa::Join<aod::FDHfCand2Prong, aod::FDHfCandMC>>;
  using FilteredCharmMcCand2Prong = FilteredCharmMcCand2Prongs::iterator;

  using FilteredCharmMcCandDstars = soa::Filtered<soa::Join<aod::FDHfCandDstar, aod::FDHfCandMC>>;
  using FilteredCharmMcCandDstar = FilteredCharmMcCandDstars::iterator;

  using FilteredCollisions = soa::Filtered<soa::Join<FDCollisions, FDColMasks>>;
  using FilteredCollision = FilteredCollisions::iterator;

  using FilteredMcColisions = soa::Filtered<soa::Join<aod::FDCollisions, FDColMasks, aod::FDMCCollLabels>>;
  using FilteredMcColision = FilteredMcColisions::iterator;

  using FilteredFDMcParts = soa::Filtered<soa::Join<aod::FDParticles, aod::FDParticlesIndex, aod::FDExtParticles, aod::FDMCLabels, aod::FDExtMCLabels>>;
  using FilteredFDMcPart = FilteredFDMcParts::iterator;

  using FilteredFDParticles = soa::Filtered<soa::Join<aod::FDParticles, aod::FDExtParticles, aod::FDParticlesIndex, aod::FDTrkTimeStamp>>;
  using FilteredFDParticle = FilteredFDParticles::iterator;

  Filter eventMultiplicity = aod::femtodreamcollision::multNtr >= eventSel.multMin && aod::femtodreamcollision::multNtr <= eventSel.multMax;
  Filter eventMultiplicityPercentile = aod::femtodreamcollision::multV0M >= eventSel.multPercentileMin && aod::femtodreamcollision::multV0M <= eventSel.multPercentileMax;
  Filter hfCandSelFilter = aod::fdhf::candidateSelFlag >= charmSel.charmHadCandSel;
  Filter hfMcSelFilter = (nabs(aod::fdhf::flagMc) == charmSel.charmHadMcSel);

  Preslice<FilteredFDParticles> perCol = aod::femtodreamparticle::fdCollisionId;
  Preslice<FilteredCharmCand3Prongs> perHf3ProngByCol = aod::femtodreamparticle::fdCollisionId;
  Preslice<FilteredCharmCand2Prongs> perHf2ProngByCol = aod::femtodreamparticle::fdCollisionId;
  Preslice<FilteredCharmCandDstars> perHfDstarByCol = aod::femtodreamparticle::fdCollisionId;

  /// Partition for particle Lambda
  Partition<FilteredFDParticles> partitionLambda = (aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kV0)) &&
                                                   ((aod::femtodreamparticle::cut & v0Sel.cutBit) == v0Sel.cutBit) &&
                                                   (aod::femtodreamparticle::pt > v0Sel.ptV0Min) &&
                                                   (aod::femtodreamparticle::pt < v0Sel.ptV0Max) &&
                                                   (aod::femtodreamparticle::eta > v0Sel.etaV0Min) &&
                                                   (aod::femtodreamparticle::eta < v0Sel.etaV0Max) &&
                                                   (aod::femtodreamparticle::mLambda > v0Sel.invMassV0Min) &&
                                                   (aod::femtodreamparticle::mLambda < v0Sel.invMassV0Max) &&
                                                   (aod::femtodreamparticle::mAntiLambda > v0Sel.invMassAntiV0Min) &&
                                                   (aod::femtodreamparticle::mAntiLambda < v0Sel.invMassAntiV0Max);

  /// Partition for particle K0Short
  Partition<FilteredFDParticles> partitionK0Short = (aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kV0K0Short)) &&
                                                    ((aod::femtodreamparticle::cut & v0Sel.cutBit) == v0Sel.cutBit) &&
                                                    (aod::femtodreamparticle::pt > v0Sel.ptV0Min) &&
                                                    (aod::femtodreamparticle::pt < v0Sel.ptV0Max) &&
                                                    (aod::femtodreamparticle::eta > v0Sel.etaV0Min) &&
                                                    (aod::femtodreamparticle::eta < v0Sel.etaV0Max) &&
                                                    (aod::femtodreamparticle::mLambda > v0Sel.invMassV0Min) &&
                                                    (aod::femtodreamparticle::mLambda < v0Sel.invMassV0Max) &&
                                                    (aod::femtodreamparticle::mAntiLambda > v0Sel.invMassAntiV0Min) &&
                                                    (aod::femtodreamparticle::mAntiLambda < v0Sel.invMassAntiV0Max);

  /// Partition for particle 2
  Partition<FilteredCharmCand3Prongs> partitionCharmHadron3Prong = aod::fdhf::bdtBkg < charmSel.charmHadBkgBDTmax && aod::fdhf::bdtFD < charmSel.charmHadFdBDTmax && aod::fdhf::bdtFD > charmSel.charmHadFdBDTmin&& aod::fdhf::bdtPrompt<charmSel.charmHadPromptBDTmax && aod::fdhf::bdtPrompt> charmSel.charmHadPromptBDTmin;
  Partition<FilteredCharmCand2Prongs> partitionCharmHadron2Prong = aod::fdhf::bdtBkg < charmSel.charmHadBkgBDTmax && aod::fdhf::bdtFD < charmSel.charmHadFdBDTmax && aod::fdhf::bdtFD > charmSel.charmHadFdBDTmin&& aod::fdhf::bdtPrompt<charmSel.charmHadPromptBDTmax && aod::fdhf::bdtPrompt> charmSel.charmHadPromptBDTmin;
  Partition<FilteredCharmCandDstars> partitionCharmHadronDstar = aod::fdhf::bdtBkg < charmSel.charmHadBkgBDTmax && aod::fdhf::bdtFD < charmSel.charmHadFdBDTmax && aod::fdhf::bdtFD > charmSel.charmHadFdBDTmin&& aod::fdhf::bdtPrompt<charmSel.charmHadPromptBDTmax && aod::fdhf::bdtPrompt> charmSel.charmHadPromptBDTmin;

  Partition<FilteredCharmMcCand3Prongs> partitionMcCharmHadron3Prong = aod::fdhf::originMcRec == OriginRecPrompt || aod::fdhf::originMcRec == OriginRecFD;
  Partition<FilteredCharmMcCand2Prongs> partitionMcCharmHadron2Prong = aod::fdhf::originMcRec == OriginRecPrompt || aod::fdhf::originMcRec == OriginRecFD;
  Partition<FilteredCharmMcCandDstars> partitionMcCharmHadronDstar = aod::fdhf::originMcRec == OriginRecPrompt || aod::fdhf::originMcRec == OriginRecFD;
  struct : ConfigurableGroup {
    /// Axis configurables
    ConfigurableAxis dummy{"dummy", {1, 0, 1}, "dummy axis"};
    /// Binning configurables
    ConfigurableAxis bin4Dkstar{"bin4Dkstar", {1500, 0., 6.}, "binning kstar for the 4Dimensional plot: k* vs multiplicity vs multiplicity percentile vs mT (set <<Confuse4D>> to true in order to use)"};
    ConfigurableAxis bin4DMult{"bin4DMult", {VARIABLE_WIDTH, 0.0f, 4.0f, 8.0f, 12.0f, 16.0f, 20.0f, 24.0f, 28.0f, 32.0f, 36.0f, 40.0f, 44.0f, 48.0f, 52.0f, 56.0f, 60.0f, 64.0f, 68.0f, 72.0f, 76.0f, 80.0f, 84.0f, 88.0f, 92.0f, 96.0f, 100.0f, 200.0f}, "multiplicity Binning for the 4Dimensional plot: k* vs multiplicity vs multiplicity percentile vs mT (set <<Confuse4D>> to true in order to use)"};
    ConfigurableAxis bin4DmT{"bin4DmT", {VARIABLE_WIDTH, 1.02f, 1.14f, 1.20f, 1.26f, 1.38f, 1.56f, 1.86f, 4.50f}, "mT Binning for the 4Dimensional plot: k* vs multiplicity vs multiplicity percentile vs mT (set <<Confuse4D>> to true in order to use)"};
    ConfigurableAxis bin4DmultPercentile{"bin4DmultPercentile", {10, 0.0f, 100.0f}, "multiplicity percentile Binning for the 4Dimensional plot: k* vs multiplicity vs multiplicity percentile vs mT (set <<Confuse4D>> to true in order to use)"};
    ConfigurableAxis binInvMass{"binInvMass", {400, 2.10, 2.50}, "InvMass binning"};
    ConfigurableAxis binpTCharm{"binpTCharm", {360, 0, 36}, "pT binning of charm hadron"};
    ConfigurableAxis binTempFitVarV0{"binTempFitVarV0", {300, 0.9, 1}, "binning of the TempFitVar in the pT vs. TempFitVar plot (charm))"};
    ConfigurableAxis binTempFitVarV0Child{"binTempFitVarV0Child", {300, -0.15, 0.15}, "binning of the TempFitVar in the pT vs. TempFitVar plot (V0 child)"};
    ConfigurableAxis binmT{"binmT", {225, 0., 7.5}, "binning mT"};
    ConfigurableAxis binmultTempFit{"binmultTempFit", {1, 0, 1}, "multiplicity Binning for the TempFitVar plot"};
    ConfigurableAxis binMulPercentile{"binMulPercentile", {10, 0.0f, 100.0f}, "multiplicity percentile Binning"};
    ConfigurableAxis binpTV0{"binpTV0", {50, 0.5, 10.05}, "pT binning of the pT vs. TempFitVar plot (Track)"};
    ConfigurableAxis binpTV0Child{"binpTV0Child", {20, 0.5, 4.05}, "pT binning of the pT vs. TempFitVar plot (V0 Child)"};
    ConfigurableAxis binEta{"binEta", {{200, -1.5, 1.5}}, "eta binning"};
    ConfigurableAxis binPhi{"binPhi", {{200, 0, o2::constants::math::TwoPI}}, "phi binning"};
    ConfigurableAxis binkT{"binkT", {150, 0., 9.}, "binning kT"};
    ConfigurableAxis binkstar{"binkstar", {1500, 0., 6.}, "binning kstar"};
    ConfigurableAxis binNSigmaTPC{"binNSigmaTPC", {1600, -8, 8}, "Binning of Nsigma TPC plot"};
    ConfigurableAxis binNSigmaTOF{"binNSigmaTOF", {3000, -15, 15}, "Binning of the Nsigma TOF plot"};
    ConfigurableAxis binNSigmaTPCTOF{"binNSigmaTPCTOF", {3000, -15, 15}, "Binning of the Nsigma TPC+TOF plot"};
    ConfigurableAxis binTPCClusters{"binTPCClusters", {163, -0.5, 162.5}, "Binning of TPC found clusters plot"};
    ConfigurableAxis mixingBinMult{"mixingBinMult", {VARIABLE_WIDTH, 0.0f, 20.0f, 60.0f, 200.0f}, "Mixing bins - multiplicity"};
    ConfigurableAxis mixingBinMultPercentile{"mixingBinMultPercentile", {VARIABLE_WIDTH, 0.0f, 100.f}, "Mixing bins - multiplicity percentile"};
    ConfigurableAxis mixingBinVztx{"mixingBinVztx", {VARIABLE_WIDTH, -10.0f, -4.f, 0.f, 4.f, 10.f}, "Mixing bins - z-vertex"};
  } AxisBinning;

  ColumnBinningPolicy<aod::collision::PosZ, aod::femtodreamcollision::MultNtr> colBinningMult{{AxisBinning.mixingBinVztx, AxisBinning.mixingBinMult}, true};
  ColumnBinningPolicy<aod::collision::PosZ, aod::femtodreamcollision::MultV0M> colBinningMultPercentile{{AxisBinning.mixingBinVztx, AxisBinning.mixingBinMultPercentile}, true};
  ColumnBinningPolicy<aod::collision::PosZ, aod::femtodreamcollision::MultNtr, aod::femtodreamcollision::MultV0M> colBinningMultMultPercentile{{AxisBinning.mixingBinVztx, AxisBinning.mixingBinMult, AxisBinning.mixingBinMultPercentile}, true};

  FemtoDreamContainer<femtoDreamContainer::EventType::same, femtoDreamContainer::Observable::kstar> sameEventCont;
  FemtoDreamContainer<femtoDreamContainer::EventType::mixed, femtoDreamContainer::Observable::kstar> mixedEventCont;
  FemtoDreamPairCleaner<aod::femtodreamparticle::ParticleType::kV0, aod::femtodreamparticle::ParticleType::kCharmHadron3Prong> pairCleaner3Prong;
  FemtoDreamDetaDphiStar<aod::femtodreamparticle::ParticleType::kV0, aod::femtodreamparticle::ParticleType::kCharmHadron3Prong> pairCloseRejectionSE3Prong;
  FemtoDreamDetaDphiStar<aod::femtodreamparticle::ParticleType::kV0, aod::femtodreamparticle::ParticleType::kCharmHadron3Prong> pairCloseRejectionME3Prong;

  FemtoDreamPairCleaner<aod::femtodreamparticle::ParticleType::kV0, aod::femtodreamparticle::ParticleType::kCharmHadron2Prong> pairCleaner2Prong;
  FemtoDreamDetaDphiStar<aod::femtodreamparticle::ParticleType::kV0, aod::femtodreamparticle::ParticleType::kCharmHadron2Prong> pairCloseRejectionSE2Prong;
  FemtoDreamDetaDphiStar<aod::femtodreamparticle::ParticleType::kV0, aod::femtodreamparticle::ParticleType::kCharmHadron2Prong> pairCloseRejectionME2Prong;

  FemtoDreamPairCleaner<aod::femtodreamparticle::ParticleType::kV0, aod::femtodreamparticle::ParticleType::kCharmHadronDstar> pairCleanerDstar;
  FemtoDreamDetaDphiStar<aod::femtodreamparticle::ParticleType::kV0, aod::femtodreamparticle::ParticleType::kCharmHadronDstar> pairCloseRejectionSEDstar;
  FemtoDreamDetaDphiStar<aod::femtodreamparticle::ParticleType::kV0, aod::femtodreamparticle::ParticleType::kCharmHadronDstar> pairCloseRejectionMEDstar;

  femtodreamcollision::BitMaskType bitMask = 1 << 0;

  /// Histogramming for particle 1
  FemtoDreamParticleHisto<aod::femtodreamparticle::ParticleType::kV0, 1> v0HistoPartOne;
  FemtoDreamParticleHisto<aod::femtodreamparticle::ParticleType::kV0, 5> v0HistoPartOneSelected;
  FemtoDreamParticleHisto<aod::femtodreamparticle::ParticleType::kV0Child, 3> posChildHistos;
  FemtoDreamParticleHisto<aod::femtodreamparticle::ParticleType::kV0Child, 4> negChildHistos;
  /// Histogramming for Event
  FemtoDreamEventHisto eventHisto;
  /// Histogram output
  HistogramRegistry registry{"CorrelationsAndQA", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry registryMixQa{"registryMixQa"};
  HistogramRegistry registryCharmHadronQa{"registryCharmHadronQa"};

  float massOne = o2::analysis::femtoDream::getMass(v0Sel.pdgCodeV0);
  float massTwo = o2::analysis::femtoDream::getMass(charmSel.charmHadPDGCode);
  int8_t partSign = 0;
  int64_t processType = 0;

  void init(InitContext& /*context*/)
  {
    std::array<bool, 4> processes = {doprocessDataLcK0s, doprocessDataDplusK0s, doprocessDataD0K0s, doprocessDataDstarK0s};
    if (std::accumulate(processes.begin(), processes.end(), 0) != 1) {
      LOGP(fatal, "One and only one process function must be enabled at a time.");
    }
    bool process3Prong = doprocessDataLcK0s || doprocessDataDplusK0s;
    bool process2Prong = doprocessDataD0K0s;
    bool processDstar = doprocessDataDstarK0s;

    // setup columnpolicy for binning
    colBinningMult = {{AxisBinning.mixingBinVztx, AxisBinning.mixingBinMult}, true};
    colBinningMultPercentile = {{AxisBinning.mixingBinVztx, AxisBinning.mixingBinMultPercentile}, true};
    colBinningMultMultPercentile = {{AxisBinning.mixingBinVztx, AxisBinning.mixingBinMult, AxisBinning.mixingBinMultPercentile}, true};
    eventHisto.init(&registry);
    v0HistoPartOne.init(&registry, AxisBinning.binmultTempFit, AxisBinning.dummy, AxisBinning.binpTV0, AxisBinning.dummy, AxisBinning.dummy, AxisBinning.binTempFitVarV0, AxisBinning.dummy, AxisBinning.dummy, AxisBinning.dummy, AxisBinning.dummy, AxisBinning.dummy, AxisBinning.dummy, isMc, v0Sel.pdgCodeV0);
    v0HistoPartOneSelected.init(&registry, AxisBinning.binmultTempFit, AxisBinning.dummy, AxisBinning.binpTV0, AxisBinning.dummy, AxisBinning.dummy, AxisBinning.binTempFitVarV0, AxisBinning.dummy, AxisBinning.dummy, AxisBinning.dummy, AxisBinning.dummy, AxisBinning.dummy, AxisBinning.dummy, isMc, v0Sel.pdgCodeV0);
    posChildHistos.init(&registry, AxisBinning.binmultTempFit, AxisBinning.dummy, AxisBinning.binpTV0Child, AxisBinning.dummy, AxisBinning.dummy, AxisBinning.binTempFitVarV0Child, AxisBinning.dummy, AxisBinning.dummy, AxisBinning.dummy, AxisBinning.dummy, AxisBinning.dummy, AxisBinning.dummy, false, 0);
    negChildHistos.init(&registry, AxisBinning.binmultTempFit, AxisBinning.dummy, AxisBinning.binpTV0Child, AxisBinning.dummy, AxisBinning.dummy, AxisBinning.binTempFitVarV0Child, AxisBinning.dummy, AxisBinning.dummy, AxisBinning.dummy, AxisBinning.dummy, AxisBinning.dummy, AxisBinning.dummy, false, 0);
    sameEventCont.init(&registry,
                       AxisBinning.binkstar, AxisBinning.binpTV0, AxisBinning.binkT, AxisBinning.binmT, AxisBinning.mixingBinMult, AxisBinning.mixingBinMultPercentile,
                       AxisBinning.bin4Dkstar, AxisBinning.bin4DmT, AxisBinning.bin4DMult, AxisBinning.bin4DmultPercentile,
                       isMc, use4D, extendedPlots,
                       highkstarCut,
                       smearingByOrigin, AxisBinning.binInvMass);

    sameEventCont.setPDGCodes(v0Sel.pdgCodeV0, charmSel.charmHadPDGCode);
    mixedEventCont.init(&registry,
                        AxisBinning.binkstar, AxisBinning.binpTV0, AxisBinning.binkT, AxisBinning.binmT, AxisBinning.mixingBinMult, AxisBinning.mixingBinMultPercentile,
                        AxisBinning.bin4Dkstar, AxisBinning.bin4DmT, AxisBinning.bin4DMult, AxisBinning.bin4DmultPercentile,
                        isMc, use4D, extendedPlots,
                        highkstarCut,
                        smearingByOrigin, AxisBinning.binInvMass);

    mixedEventCont.setPDGCodes(v0Sel.pdgCodeV0, charmSel.charmHadPDGCode);
    registryMixQa.add("MixingQA/hSECollisionBins", "; bin; Entries", kTH1F, {{120, -0.5, 119.5}});
    registryMixQa.add("MixingQA/hSECollisionPool", "; Vz (cm); Mul", kTH2F, {{100, -10, 10}, {200, 0, 200}});
    registryMixQa.add("MixingQA/hMECollisionBins", "; bin; Entries", kTH1F, {{120, -0.5, 119.5}});
    registryCharmHadronQa.add("CharmHadronQA/hPtVsMass", "; #it{p}_{T} (GeV/#it{c}); inv. mass (GeV/#it{c}^{2})", kTH2F, {AxisBinning.binpTCharm, AxisBinning.binInvMass});

    if (useCPR.value && process3Prong) {
      pairCleaner3Prong.init(&registry);
      pairCloseRejectionSE3Prong.init(&registry, &registry, cprDeltaPhiMax.value, cprDeltaEtaMax.value, cprPlotPerRadii.value, 1);
      pairCloseRejectionME3Prong.init(&registry, &registry, cprDeltaPhiMax.value, cprDeltaEtaMax.value, cprPlotPerRadii.value, 2);
    } else if (useCPR.value && process2Prong) {
      pairCleaner2Prong.init(&registry);
      pairCloseRejectionSE2Prong.init(&registry, &registry, cprDeltaPhiMax.value, cprDeltaEtaMax.value, cprPlotPerRadii.value, 1);
      pairCloseRejectionME2Prong.init(&registry, &registry, cprDeltaPhiMax.value, cprDeltaEtaMax.value, cprPlotPerRadii.value, 2);
    } else if (useCPR.value && processDstar) {
      pairCleanerDstar.init(&registry);
      pairCloseRejectionSEDstar.init(&registry, &registry, cprDeltaPhiMax.value, cprDeltaEtaMax.value, cprPlotPerRadii.value, 1);
      pairCloseRejectionMEDstar.init(&registry, &registry, cprDeltaPhiMax.value, cprDeltaEtaMax.value, cprPlotPerRadii.value, 2);
    }
  }

  template <typename CollisionType>
  void fillCollision(CollisionType const& col)
  {
    registryMixQa.fill(HIST("MixingQA/hSECollisionBins"), colBinningMult.getBin({col.posZ(), col.multNtr()}));
    registryMixQa.fill(HIST("MixingQA/hSECollisionPool"), col.posZ(), col.multNtr());
  }

  /// Compute the charm hadron candidates mass with the daughter masses
  /// assumes the candidate is either a D+ or Λc+ or D0 or Dstar
  template <DecayChannel Channel, typename Candidate>
  float getCharmHadronMass(const Candidate& cand, bool ReturnDaughMass = false)
  {
    float invMass = 0.0f;
    if constexpr (Channel == DecayChannel::LcToPKPi) {
      if (cand.candidateSelFlag() == 1) {
        invMass = cand.m(std::array{MassProton, MassKPlus, MassPiPlus});
        return invMass;
      }
      invMass = cand.m(std::array{MassPiPlus, MassKPlus, MassProton});
      return invMass;
    } else if constexpr (Channel == DecayChannel::DplusToPiKPi) { // D+ → π K π (PDG: 411)
      invMass = cand.m(std::array{MassPiPlus, MassKPlus, MassPiPlus});
      return invMass;
    } else if constexpr (Channel == DecayChannel::D0ToPiK) { // D0 → π K  (PDG: 421)
      if (cand.candidateSelFlag() == 1) {
        invMass = cand.m(std::array{MassPiPlus, MassKPlus});
        return invMass;
      } else {
        invMass = cand.m(std::array{MassKPlus, MassPiPlus});
        return invMass;
      }
    } else if constexpr (Channel == DecayChannel::DstarToD0Pi) { // D* → D0π (PDG: 413)
      float mDstar = 0.f;
      float mD0 = 0.f;
      if (cand.charge() > 0.f) {
        mDstar = cand.m(std::array{MassPiPlus, MassKPlus, MassPiPlus});
        mD0 = cand.mDaughD0(std::array{MassPiPlus, MassKPlus});
      } else {
        mDstar = cand.m(std::array{MassKPlus, MassPiPlus, MassPiPlus});
        mD0 = cand.mDaughD0(std::array{MassKPlus, MassPiPlus});
      }
      if (ReturnDaughMass) {
        return mD0;
      } else {
        return mDstar - mD0;
      }
    }
    // Add more channels as needed
    return 0.f;
  }

  /// This function processes the same event and takes care of all the histogramming
  template <bool IsMc, DecayChannel Channel, typename CandType, typename PartitionType, typename FDParticles, typename Collision>
  void doSameEvent(CandType& sliceCharmHad, PartitionType& sliceV01, FDParticles const& femtoParts, Collision const& col)
  {
    fillCollision(col);
    processType = 1; // for same event
    for (auto const& [p1, p2] : combinations(CombinationsFullIndexPolicy(sliceV01, sliceCharmHad))) {

      if constexpr (Channel == DecayChannel::D0ToPiK) {
        if (p1.childrenIds()[0] == p2.prong0Id() || p1.childrenIds()[0] == p2.prong1Id() || p1.childrenIds()[1] == p2.prong0Id() || p1.childrenIds()[1] == p2.prong1Id())
          continue;

        if (useCPR.value) {
          if (pairCloseRejectionSE2Prong.isClosePair(p1, p2, femtoParts, col.magField())) {
            continue;
          }
        }
        if (!pairCleaner2Prong.isCleanPair(p1, p2, femtoParts)) {
          continue;
        }
      } else if constexpr (Channel == DecayChannel::LcToPKPi || Channel == DecayChannel::DplusToPiKPi) {
        if (p1.childrenIds()[0] == p2.prong0Id() || p1.childrenIds()[0] == p2.prong1Id() || p1.childrenIds()[0] == p2.prong2Id() || p1.childrenIds()[1] == p2.prong0Id() || p1.childrenIds()[1] == p2.prong1Id() || p1.childrenIds()[1] == p2.prong2Id())
          continue;
        if (useCPR.value) {
          if (pairCloseRejectionSE3Prong.isClosePair(p1, p2, femtoParts, col.magField())) {
            continue;
          }
        }
        if (!pairCleaner3Prong.isCleanPair(p1, p2, femtoParts)) {
          continue;
        }
      } else if constexpr (Channel == DecayChannel::DstarToD0Pi) {
        if (p1.childrenIds()[0] == p2.prong0Id() || p1.childrenIds()[0] == p2.prong1Id() || p1.childrenIds()[0] == p2.prong2Id() || p1.childrenIds()[1] == p2.prong0Id() || p1.childrenIds()[1] == p2.prong1Id() || p1.childrenIds()[1] == p2.prong2Id())
          continue;
        if (useCPR.value) {
          if (pairCloseRejectionSEDstar.isClosePair(p1, p2, femtoParts, col.magField())) {
            continue;
          }
        }

        if (!pairCleanerDstar.isCleanPair(p1, p2, femtoParts)) {
          continue;
        }
      }
      // v0 daughters selection
      const auto& posChild = femtoParts.iteratorAt(p1.index() - 2);
      const auto& negChild = femtoParts.iteratorAt(p1.index() - 1);
      if (((posChild.cut() & v0Sel.childPosCutBit) == v0Sel.childPosCutBit) &&
          ((posChild.pidcut() & v0Sel.childPosTPCBit) == v0Sel.childPosTPCBit) &&
          ((negChild.cut() & v0Sel.childNegCutBit) == v0Sel.childNegCutBit) &&
          ((negChild.pidcut() & v0Sel.childNegTPCBit) == v0Sel.childNegTPCBit)) {

        v0HistoPartOneSelected.fillQA<IsMc, true>(p1, static_cast<aod::femtodreamparticle::MomentumType>(confTempFitVarMomentum.value), col.multNtr(), col.multV0M());
        posChildHistos.fillQA<false, false>(posChild, static_cast<aod::femtodreamparticle::MomentumType>(confTempFitVarMomentum.value), col.multNtr(), col.multV0M());
        negChildHistos.fillQA<false, false>(negChild, static_cast<aod::femtodreamparticle::MomentumType>(confTempFitVarMomentum.value), col.multNtr(), col.multV0M());
      } else {
        continue;
      }

      float kstar = FemtoDreamMath::getkstar(p1, massOne, p2, massTwo);
      if (kstar > highkstarCut) {
        continue;
      }
      float invMass = getCharmHadronMass<Channel>(p2);

      if (invMass < charmSel.charmHadMinInvMass || invMass > charmSel.charmHadMaxInvMass) {
        continue;
      }

      if (p2.pt() < charmSel.charmHadMinPt || p2.pt() > charmSel.charmHadMaxPt) {
        continue;
      }

      float invMassV0 = 0.f;
      if (p1.sign() > 0)
        invMassV0 = p1.mLambda();
      else
        invMassV0 = p1.mAntiLambda();

      float chargeV0 = 0.;
      if ((p1.cut() & p1.sign()) == CutBitChargePositive) {
        chargeV0 = PositiveCharge;
      } else {
        chargeV0 = NegativeCharge;
      }
      int pairSign = 0;
      if (chargeV0 == p2.charge()) {
        pairSign = LikeSignPair;
      } else {
        pairSign = UnLikeSignPair;
      }

      int charmHadMc = 0;
      int originType = 0;
      if constexpr (IsMc) {
        charmHadMc = p2.flagMc();
        originType = p2.originMcRec();
      }

      rowFemtoResultPairs(
        invMass,
        p2.pt(),
        p1.pt(),
        p2.bdtBkg(),
        p2.bdtPrompt(),
        p2.bdtFD(),
        kstar,
        FemtoDreamMath::getkT(p1, massOne, p2, massTwo),
        FemtoDreamMath::getmT(p1, massOne, p2, massTwo),
        col.multNtr(),
        col.multV0M(),
        p2.charge(),
        pairSign,
        invMassV0,
        processType,
        charmHadMc,
        originType);

      sameEventCont.setPair<IsMc, true>(p1, p2, col.multNtr(), col.multV0M(), use4D, extendedPlots, smearingByOrigin);
    }
  }

  template <bool IsMc, DecayChannel Channel, typename CollisionType, typename PartitionType1, typename PartitionType2, typename FDParticles, typename BinningType>
  void doMixedEvent(CollisionType const& cols, PartitionType1& charms, PartitionType2& v0s, FDParticles const& femtoParts, BinningType policy)
  {
    processType = 2; // for mixed event
    // Mixed events that contain the pair of interest
    Partition<CollisionType> partitionMaskedCol1 = (aod::femtodreamcollision::bitmaskTrackOne & bitMask) == bitMask;
    partitionMaskedCol1.bindTable(cols);

    Partition<CollisionType> partitionMaskedCol2 = (aod::femtodreamcollision::bitmaskTrackTwo & bitMask) == bitMask;
    partitionMaskedCol2.bindTable(cols);

    for (auto const& [collision1, collision2] : combinations(soa::CombinationsBlockFullIndexPolicy(policy, mixSetting.mixingDepth, -1, *partitionMaskedCol1.mFiltered, *partitionMaskedCol2.mFiltered))) {
      // make sure that tracks in the same events are not mixed
      if (collision1.globalIndex() == collision2.globalIndex()) {
        continue;
      }

      const int multiplicityCol = collision1.multNtr();
      registryMixQa.fill(HIST("MixingQA/hMECollisionBins"), colBinningMult.getBin({collision1.posZ(), multiplicityCol}));

      auto sliceV01 = v0s->sliceByCached(aod::femtodreamparticle::fdCollisionId, collision1.globalIndex(), cache);
      auto sliceCharmHad = charms->sliceByCached(aod::femtodreamparticle::fdCollisionId, collision2.globalIndex(), cache);
      for (const auto& [p1, p2] : combinations(CombinationsFullIndexPolicy(sliceV01, sliceCharmHad))) {

        if constexpr (Channel == DecayChannel::D0ToPiK) {

          if (useCPR.value) {
            if (pairCloseRejectionME2Prong.isClosePair(p1, p2, femtoParts, collision1.magField())) {
              continue;
            }
          }

          if (!pairCleaner2Prong.isCleanPair(p1, p2, femtoParts)) {
            continue;
          }
        }

        if constexpr (Channel == DecayChannel::DplusToPiKPi || Channel == DecayChannel::LcToPKPi) {

          if (useCPR.value) {
            if (pairCloseRejectionME3Prong.isClosePair(p1, p2, femtoParts, collision1.magField())) {
              continue;
            }
          }

          if (!pairCleaner3Prong.isCleanPair(p1, p2, femtoParts)) {
            continue;
          }
        }

        if constexpr (Channel == DecayChannel::DstarToD0Pi) {

          if (useCPR.value) {
            if (pairCloseRejectionME3Prong.isClosePair(p1, p2, femtoParts, collision1.magField())) {
              continue;
            }
          }

          if (!pairCleanerDstar.isCleanPair(p1, p2, femtoParts)) {
            continue;
          }
        }

        float kstar = FemtoDreamMath::getkstar(p1, massOne, p2, massTwo);
        if (kstar > highkstarCut) {
          continue;
        }

        float invMass = getCharmHadronMass<Channel>(p2);

        if (invMass < charmSel.charmHadMinInvMass || invMass > charmSel.charmHadMaxInvMass) {
          continue;
        }

        if (p2.pt() < charmSel.charmHadMinPt || p2.pt() > charmSel.charmHadMaxPt) {
          continue;
        }

        float invMassV0 = 0.f;
        if (p1.sign() > 0)
          invMassV0 = p1.mLambda();
        else
          invMassV0 = p1.mAntiLambda();

        float chargeV0 = 0.;
        if ((p1.cut() & p1.sign()) == CutBitChargePositive) {
          chargeV0 = PositiveCharge;
        } else {
          chargeV0 = NegativeCharge;
        }
        int pairSign = 0;
        if (chargeV0 == p2.charge()) {
          pairSign = LikeSignPair;
        } else {
          pairSign = UnLikeSignPair;
        }

        int charmHadMc = 0;
        int originType = 0;
        if constexpr (IsMc) {
          charmHadMc = p2.flagMc();
          originType = p2.originMcRec();
        }

        rowFemtoResultPairs(
          invMass,
          p2.pt(),
          p1.pt(),
          p2.bdtBkg(),
          p2.bdtPrompt(),
          p2.bdtFD(),
          kstar,
          FemtoDreamMath::getkT(p1, massOne, p2, massTwo),
          FemtoDreamMath::getmT(p1, massOne, p2, massTwo),
          collision1.multNtr(),
          collision1.multV0M(),
          p2.charge(),
          pairSign,
          invMassV0,
          processType,
          charmHadMc,
          originType);

        mixedEventCont.setPair<IsMc, true>(p1, p2, collision1.multNtr(), collision1.multV0M(), use4D, extendedPlots, smearingByOrigin);
      }
    }
  }
  template <bool IsMc, DecayChannel Channel, typename V0Type, typename CandType, typename CollType, typename FDParticles>
  void fillTables(const CollType& col,
                  const V0Type& sliceV01,
                  const CandType& sliceCharmHad,
                  const FDParticles& femtoParts)
  {
    int64_t timeStamp = -999;

    // ---- Fill Charm-Hadron Table ----
    for (auto const& part : sliceCharmHad) {
      float invMass = getCharmHadronMass<Channel>(part);
      registryCharmHadronQa.fill(HIST("CharmHadronQA/hPtVsMass"), part.pt(), invMass);
      timeStamp = part.timeStamp();

      if constexpr (Channel == DecayChannel::DplusToPiKPi || Channel == DecayChannel::LcToPKPi) {

        rowFemtoResultCharm3Prong(
          col.globalIndex(),
          timeStamp,
          invMass,
          part.pt(),
          part.eta(),
          part.phi(),
          part.prong0Id(),
          part.prong1Id(),
          part.prong2Id(),
          part.charge(),
          part.bdtBkg(),
          part.bdtPrompt(),
          part.bdtFD());
      } else if constexpr (Channel == DecayChannel::D0ToPiK) {
        rowFemtoResultCharm2Prong(
          col.globalIndex(),
          timeStamp,
          invMass,
          part.pt(),
          part.eta(),
          part.phi(),
          part.prong0Id(),
          part.prong1Id(),
          part.charge(),
          part.bdtBkg(),
          part.bdtPrompt(),
          part.bdtFD());
      } else if constexpr (Channel == DecayChannel::DstarToD0Pi) {
        float invMassD0 = getCharmHadronMass<Channel>(part, true);
        rowFemtoResultCharmDstar(
          col.globalIndex(),
          timeStamp,
          invMass,
          invMassD0,
          part.pt(),
          part.eta(),
          part.phi(),
          part.prong0Id(),
          part.prong1Id(),
          part.prong2Id(),
          part.charge(),
          part.bdtBkg(),
          part.bdtPrompt(),
          part.bdtFD());
      }
    }

    // ---- Fill V0 Table ----
    for (auto const& part : sliceV01) {

      // v0 daughters selection
      const auto& posChild = femtoParts.iteratorAt(part.index() - 2);
      const auto& negChild = femtoParts.iteratorAt(part.index() - 1);
      if (((posChild.cut() & v0Sel.childPosCutBit) == v0Sel.childPosCutBit) &&
          ((posChild.pidcut() & v0Sel.childPosTPCBit) == v0Sel.childPosTPCBit) &&
          ((negChild.cut() & v0Sel.childNegCutBit) == v0Sel.childNegCutBit) &&
          ((negChild.pidcut() & v0Sel.childNegTPCBit) == v0Sel.childNegTPCBit)) {
        v0HistoPartOne.fillQA<IsMc, true>(part, static_cast<aod::femtodreamparticle::MomentumType>(confTempFitVarMomentum.value), col.multNtr(), col.multV0M());
      } else {
        continue;
      }

      std::vector<int> childIds = {part.childrenIds()[0], part.childrenIds()[1]};
      timeStamp = part.timeStamp();
      rowFemtoResultV0(
        col.globalIndex(),
        timeStamp,
        part.pt(),
        part.eta(),
        part.phi(),
        childIds,
        part.sign(),
        part.mLambda(),
        part.mAntiLambda(),
        part.tpcNClsFound(),
        part.tpcNClsFindable(),
        part.tpcNClsCrossedRows(),
        part.tpcNSigmaPr(),
        part.tofNSigmaPr());
    }

    // ---- Fill Collision Table ----
    if (sliceCharmHad.size() > 0 || sliceV01.size() > 0) {
      rowFemtoResultColl(
        col.globalIndex(),
        timeStamp,
        col.posZ(),
        col.multNtr());
    }
  }

  void processDataLcK0s(FilteredCollisions const& cols,
                        FilteredFDParticles const& parts,
                        FilteredCharmCand3Prongs const&)
  {
    for (const auto& col : cols) {
      eventHisto.fillQA(col);
      auto sliceCharmHad = partitionCharmHadron3Prong->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
      if (fillTableWithCharm.value && sliceCharmHad.size() == 0) {
        continue;
      }
      if (v0Sel.pdgCodeV0 == kLambda0) {
        auto sliceV0 = partitionLambda->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
        fillTables<false, DecayChannel::LcToPKPi>(col, sliceV0, sliceCharmHad, parts);
        if (sliceCharmHad.size() > 0 && sliceV0.size() > 0) {
          doSameEvent<false, DecayChannel::LcToPKPi, FilteredCharmCand3Prongs>(sliceCharmHad, sliceV0, parts, col);
        }

      } else if (v0Sel.pdgCodeV0 == kK0Short) {
        auto sliceV0 = partitionK0Short->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
        fillTables<false, DecayChannel::LcToPKPi>(col, sliceV0, sliceCharmHad, parts);

        if (sliceCharmHad.size() > 0 && sliceV0.size() > 0) {
          doSameEvent<false, DecayChannel::LcToPKPi, FilteredCharmCand3Prongs>(sliceCharmHad, sliceV0, parts, col);
        }
      } else {
        LOG(fatal) << "Unsupported V0 PDG: " << v0Sel.pdgCodeV0 << " (allowed: 3122, 310)";
      }
    }

    if (mixSetting.doMixEvent) {
      switch (mixSetting.mixingBinPolicy) {
        case femtodreamcollision::kMult:
          doMixedEvent<false, DecayChannel::LcToPKPi, FilteredCollisions>(cols, partitionCharmHadron3Prong, partitionK0Short, parts, colBinningMult);
          break;
        case femtodreamcollision::kMultPercentile:
          doMixedEvent<false, DecayChannel::LcToPKPi, FilteredCollisions>(cols, partitionCharmHadron3Prong, partitionK0Short, parts, colBinningMultPercentile);
          break;
        case femtodreamcollision::kMultMultPercentile:
          doMixedEvent<false, DecayChannel::LcToPKPi, FilteredCollisions>(cols, partitionCharmHadron3Prong, partitionK0Short, parts, colBinningMultMultPercentile);
          break;
        default:
          LOG(fatal) << "Invalid binning policiy specifed. Breaking...";
      }
    }
  }
  PROCESS_SWITCH(HfTaskCharmHadronsV0FemtoDream, processDataLcK0s, "Enable processing LcToPKPi and V0 correlation", false);

  void processDataDplusK0s(FilteredCollisions const& cols,
                           FilteredFDParticles const& parts,
                           FilteredCharmCand3Prongs const&)
  {
    for (const auto& col : cols) {
      eventHisto.fillQA(col);
      auto sliceCharmHad = partitionCharmHadron3Prong->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
      if (fillTableWithCharm.value && sliceCharmHad.size() == 0) {
        continue;
      }
      if (v0Sel.pdgCodeV0 == kLambda0) {
        auto sliceV0 = partitionLambda->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
        fillTables<false, DecayChannel::DplusToPiKPi>(col, sliceV0, sliceCharmHad, parts);
        if (sliceCharmHad.size() > 0 && sliceV0.size() > 0) {
          doSameEvent<false, DecayChannel::DplusToPiKPi, FilteredCharmCand3Prongs>(sliceCharmHad, sliceV0, parts, col);
        }

      } else if (v0Sel.pdgCodeV0 == kK0Short) {
        auto sliceV0 = partitionK0Short->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
        fillTables<false, DecayChannel::DplusToPiKPi>(col, sliceV0, sliceCharmHad, parts);

        if (sliceCharmHad.size() > 0 && sliceV0.size() > 0) {
          doSameEvent<false, DecayChannel::DplusToPiKPi, FilteredCharmCand3Prongs>(sliceCharmHad, sliceV0, parts, col);
        }
      } else {
        LOG(fatal) << "Unsupported V0 PDG: " << v0Sel.pdgCodeV0 << " (allowed: 3122, 310)";
      }
    }

    if (mixSetting.doMixEvent) {
      switch (mixSetting.mixingBinPolicy) {
        case femtodreamcollision::kMult:
          doMixedEvent<false, DecayChannel::DplusToPiKPi, FilteredCollisions>(cols, partitionCharmHadron3Prong, partitionK0Short, parts, colBinningMult);
          break;
        case femtodreamcollision::kMultPercentile:
          doMixedEvent<false, DecayChannel::DplusToPiKPi, FilteredCollisions>(cols, partitionCharmHadron3Prong, partitionK0Short, parts, colBinningMultPercentile);
          break;
        case femtodreamcollision::kMultMultPercentile:
          doMixedEvent<false, DecayChannel::DplusToPiKPi, FilteredCollisions>(cols, partitionCharmHadron3Prong, partitionK0Short, parts, colBinningMultMultPercentile);
          break;
        default:
          LOG(fatal) << "Invalid binning policiy specifed. Breaking...";
      }
    }
  }
  PROCESS_SWITCH(HfTaskCharmHadronsV0FemtoDream, processDataDplusK0s, "Enable processing DplusToPiKPi and V0 correlation", false);

  void processDataD0K0s(FilteredCollisions const& cols,
                        FilteredFDParticles const& parts,
                        FilteredCharmCand2Prongs const&)
  {
    for (const auto& col : cols) {
      eventHisto.fillQA(col);
      auto sliceCharmHad = partitionCharmHadron2Prong->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
      if (fillTableWithCharm.value && sliceCharmHad.size() == 0) {
        continue;
      }
      if (v0Sel.pdgCodeV0 == kLambda0) {
        auto sliceV0 = partitionLambda->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
        fillTables<false, DecayChannel::D0ToPiK>(col, sliceV0, sliceCharmHad, parts);
        if (sliceCharmHad.size() > 0 && sliceV0.size() > 0) {
          doSameEvent<false, DecayChannel::D0ToPiK, FilteredCharmCand2Prongs>(sliceCharmHad, sliceV0, parts, col);
        }

      } else if (v0Sel.pdgCodeV0 == kK0Short) {
        auto sliceV0 = partitionK0Short->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
        fillTables<false, DecayChannel::D0ToPiK>(col, sliceV0, sliceCharmHad, parts);

        if (sliceCharmHad.size() > 0 && sliceV0.size() > 0) {
          doSameEvent<false, DecayChannel::D0ToPiK, FilteredCharmCand2Prongs>(sliceCharmHad, sliceV0, parts, col);
        }
      } else {
        LOG(fatal) << "Unsupported V0 PDG: " << v0Sel.pdgCodeV0 << " (allowed: 3122, 310)";
      }
    }

    if (mixSetting.doMixEvent) {
      switch (mixSetting.mixingBinPolicy) {
        case femtodreamcollision::kMult:
          doMixedEvent<false, DecayChannel::D0ToPiK, FilteredCollisions>(cols, partitionCharmHadron2Prong, partitionK0Short, parts, colBinningMult);
          break;
        case femtodreamcollision::kMultPercentile:
          doMixedEvent<false, DecayChannel::D0ToPiK, FilteredCollisions>(cols, partitionCharmHadron2Prong, partitionK0Short, parts, colBinningMultPercentile);
          break;
        case femtodreamcollision::kMultMultPercentile:
          doMixedEvent<false, DecayChannel::D0ToPiK, FilteredCollisions>(cols, partitionCharmHadron2Prong, partitionK0Short, parts, colBinningMultMultPercentile);
          break;
        default:
          LOG(fatal) << "Invalid binning policiy specifed. Breaking...";
      }
    }
  }
  PROCESS_SWITCH(HfTaskCharmHadronsV0FemtoDream, processDataD0K0s, "Enable processing D0ToPiK and V0 correlation", false);

  void processDataDstarK0s(FilteredCollisions const& cols,
                           FilteredFDParticles const& parts,
                           FilteredCharmCandDstars const&)
  {
    for (const auto& col : cols) {
      eventHisto.fillQA(col);
      auto sliceCharmHad = partitionCharmHadronDstar->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
      if (fillTableWithCharm.value && sliceCharmHad.size() == 0) {
        continue;
      }
      if (v0Sel.pdgCodeV0 == kLambda0) {
        auto sliceV0 = partitionLambda->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
        fillTables<false, DecayChannel::DstarToD0Pi>(col, sliceV0, sliceCharmHad, parts);
        if (sliceCharmHad.size() > 0 && sliceV0.size() > 0) {
          doSameEvent<false, DecayChannel::DstarToD0Pi, FilteredCharmCandDstars>(sliceCharmHad, sliceV0, parts, col);
        }

      } else if (v0Sel.pdgCodeV0 == kK0Short) {
        auto sliceV0 = partitionK0Short->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
        fillTables<false, DecayChannel::DstarToD0Pi>(col, sliceV0, sliceCharmHad, parts);

        if (sliceCharmHad.size() > 0 && sliceV0.size() > 0) {
          doSameEvent<false, DecayChannel::DstarToD0Pi, FilteredCharmCandDstars>(sliceCharmHad, sliceV0, parts, col);
        }
      } else {
        LOG(fatal) << "Unsupported V0 PDG: " << v0Sel.pdgCodeV0 << " (allowed: 3122, 310)";
      }
    }

    if (mixSetting.doMixEvent) {
      switch (mixSetting.mixingBinPolicy) {
        case femtodreamcollision::kMult:
          doMixedEvent<false, DecayChannel::DstarToD0Pi, FilteredCollisions>(cols, partitionCharmHadronDstar, partitionK0Short, parts, colBinningMult);
          break;
        case femtodreamcollision::kMultPercentile:
          doMixedEvent<false, DecayChannel::DstarToD0Pi, FilteredCollisions>(cols, partitionCharmHadronDstar, partitionK0Short, parts, colBinningMultPercentile);
          break;
        case femtodreamcollision::kMultMultPercentile:
          doMixedEvent<false, DecayChannel::DstarToD0Pi, FilteredCollisions>(cols, partitionCharmHadronDstar, partitionK0Short, parts, colBinningMultMultPercentile);
          break;
        default:
          LOG(fatal) << "Invalid binning policiy specifed. Breaking...";
      }
    }
  }
  PROCESS_SWITCH(HfTaskCharmHadronsV0FemtoDream, processDataDstarK0s, "Enable processing DstarToD0Pi and V0 correlation", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskCharmHadronsV0FemtoDream>(cfgc)};
}
