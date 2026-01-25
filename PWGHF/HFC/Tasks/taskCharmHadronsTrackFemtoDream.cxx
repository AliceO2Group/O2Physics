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

/// \file taskCharmHadronsFemtoDream.cxx
/// \brief Tasks that reads the track tables used for the pairing and builds pairs of two tracks
/// \author Ravindra SIngh, GSI, ravindra.singh@cern.ch
/// \author Biao Zhang, Heidelberg University, biao.zhang@cern.ch
/// \author Yunfan Liu, Central China Normal University, yunfan.l@cern.ch

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

using namespace o2;
using namespace o2::aod;
using namespace o2::soa;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::analysis::femtoDream;
using namespace o2::hf_decay;
using namespace o2::hf_decay::hf_cand_3prong;
using namespace o2::constants::physics;

inline o2::framework::expressions::Node coshEta(o2::framework::expressions::Node&& eta)
{
  auto e1 = std::move(eta);
  auto e2 = e1;

  return (nexp(std::move(e1)) + nexp(std::move(e2) * (-1.0f))) * 0.5f;
}

struct HfTaskCharmHadronsTrackFemtoDream {

  enum TrackCharge {
    PositiveCharge = 1,
    NegativeCharge = -1
  };

  enum PairSign {
    PairNotDefined = 0,
    LikeSignPair = 1,
    UnLikeSignPair = 2,
    ReflectedPair = 3
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
  Produces<o2::aod::FDHfTrk> rowFemtoResultTrk;
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

  /// Particle 1 (track)
  struct : ConfigurableGroup {
    Configurable<femtodreamparticle::cutContainerType> cutBitTrack1{"cutBitTrack1", 8188, "Particle 1 (Track) - Selection bit from cutCulator"};
    Configurable<int> pdgCodeTrack1{"pdgCodeTrack1", 2212, "PDG code of Particle 1 (Track)"};
    Configurable<float> pidThresTrack1{"pidThresTrack1", 0.75, "Momentum threshold for PID selection for particle 1 (Track)"};
    Configurable<femtodreamparticle::cutContainerType> tpcBitTrack1{"tpcBitTrack1", 4, "PID TPC bit from cutCulator for particle 1 (Track)"};
    Configurable<femtodreamparticle::cutContainerType> tpcTofBitTrack1{"tpcTofBitTrack1", 2, "PID TPCTOF bit from cutCulator for particle 1 (Track)"};
    Configurable<float> etaTrack1Max{"etaTrack1Max", 10., "Maximum eta of partricle 1 (Track)"};
    Configurable<float> ptTrack1Max{"ptTrack1Max", 999., "Maximum pT of partricle 1 (Track)"};
    Configurable<float> etaTrack1Min{"etaTrack1Min", -10., "Minimum eta of partricle 1 (Track)"};
    Configurable<float> ptTrack1Min{"ptTrack1Min", 0., "Minimum pT of partricle 1 (Track)"};
  } trackSel;
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
  Filter trackEtaFilterLow = ifnode(aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kTrack), aod::femtodreamparticle::eta < trackSel.etaTrack1Max, true);
  Filter trackEtaFilterUp = ifnode(aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kTrack), aod::femtodreamparticle::eta > trackSel.etaTrack1Min, true);
  Filter trackPtFilterLow = ifnode(aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kTrack), aod::femtodreamparticle::pt < trackSel.ptTrack1Max, true);
  Filter trackPtFilterUp = ifnode(aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kTrack), aod::femtodreamparticle::pt > trackSel.ptTrack1Min, true);

  Preslice<FilteredFDParticles> perCol = aod::femtodreamparticle::fdCollisionId;
  Preslice<FilteredCharmCand3Prongs> perHf3ProngByCol = aod::femtodreamparticle::fdCollisionId;
  Preslice<FilteredCharmCand2Prongs> perHf2ProngByCol = aod::femtodreamparticle::fdCollisionId;
  Preslice<FilteredCharmCandDstars> perHfDstarByCol = aod::femtodreamparticle::fdCollisionId;

  /// Partition for particle 1
  Partition<FilteredFDParticles> partitionTrk1 = (aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kTrack)) && (ncheckbit(aod::femtodreamparticle::cut, trackSel.cutBitTrack1)) && ifnode(aod::femtodreamparticle::pt * coshEta(aod::femtodreamparticle::eta) <= trackSel.pidThresTrack1, ncheckbit(aod::femtodreamparticle::pidcut, trackSel.tpcBitTrack1), ncheckbit(aod::femtodreamparticle::pidcut, trackSel.tpcTofBitTrack1));

  Partition<FilteredFDMcParts> partitionMcTrk1 = (aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kTrack)) &&
                                                 (ncheckbit(aod::femtodreamparticle::cut, trackSel.cutBitTrack1)) &&
                                                 ifnode(aod::femtodreamparticle::pt * coshEta(aod::femtodreamparticle::eta) <= trackSel.pidThresTrack1, ncheckbit(aod::femtodreamparticle::pidcut, trackSel.tpcBitTrack1), ncheckbit(aod::femtodreamparticle::pidcut, trackSel.tpcTofBitTrack1));

  /// Partition for particle 2
  Partition<FilteredCharmCand3Prongs> partitionCharmHadron3Prong = aod::fdhf::bdtBkg < charmSel.charmHadBkgBDTmax && aod::fdhf::bdtFD < charmSel.charmHadFdBDTmax && aod::fdhf::bdtFD > charmSel.charmHadFdBDTmin&& aod::fdhf::bdtPrompt<charmSel.charmHadPromptBDTmax && aod::fdhf::bdtPrompt> charmSel.charmHadPromptBDTmin;
  Partition<FilteredCharmCand2Prongs> partitionCharmHadron2Prong = aod::fdhf::bdtBkg < charmSel.charmHadBkgBDTmax && aod::fdhf::bdtFD < charmSel.charmHadFdBDTmax && aod::fdhf::bdtFD > charmSel.charmHadFdBDTmin&& aod::fdhf::bdtPrompt<charmSel.charmHadPromptBDTmax && aod::fdhf::bdtPrompt> charmSel.charmHadPromptBDTmin;
  Partition<FilteredCharmCandDstars> partitionCharmHadronDstar = aod::fdhf::bdtBkg < charmSel.charmHadBkgBDTmax && aod::fdhf::bdtFD < charmSel.charmHadFdBDTmax && aod::fdhf::bdtFD > charmSel.charmHadFdBDTmin&& aod::fdhf::bdtPrompt<charmSel.charmHadPromptBDTmax && aod::fdhf::bdtPrompt> charmSel.charmHadPromptBDTmin;

  Partition<FilteredCharmMcCand3Prongs> partitionMcCharmHadron3Prong = aod::fdhf::originMcRec == OriginRecPrompt || aod::fdhf::originMcRec == OriginRecFD;
  Partition<FilteredCharmMcCand2Prongs> partitionMcCharmHadron2Prong = aod::fdhf::originMcRec == OriginRecPrompt || aod::fdhf::originMcRec == OriginRecFD;
  Partition<FilteredCharmMcCandDstars> partitionMcCharmHadronDstar = aod::fdhf::originMcRec == OriginRecPrompt || aod::fdhf::originMcRec == OriginRecFD;

  /// Axis configurables
  ConfigurableAxis dummy{"dummy", {1, 0, 1}, "dummy axis"};
  /// Binning configurables
  ConfigurableAxis bin4Dkstar{"bin4Dkstar", {1500, 0., 6.}, "binning kstar for the 4Dimensional plot: k* vs multiplicity vs multiplicity percentile vs mT (set <<Confuse4D>> to true in order to use)"};
  ConfigurableAxis bin4DMult{"bin4DMult", {VARIABLE_WIDTH, 0.0f, 4.0f, 8.0f, 12.0f, 16.0f, 20.0f, 24.0f, 28.0f, 32.0f, 36.0f, 40.0f, 44.0f, 48.0f, 52.0f, 56.0f, 60.0f, 64.0f, 68.0f, 72.0f, 76.0f, 80.0f, 84.0f, 88.0f, 92.0f, 96.0f, 100.0f, 200.0f}, "multiplicity Binning for the 4Dimensional plot: k* vs multiplicity vs multiplicity percentile vs mT (set <<Confuse4D>> to true in order to use)"};
  ConfigurableAxis bin4DmT{"bin4DmT", {VARIABLE_WIDTH, 1.02f, 1.14f, 1.20f, 1.26f, 1.38f, 1.56f, 1.86f, 4.50f}, "mT Binning for the 4Dimensional plot: k* vs multiplicity vs multiplicity percentile vs mT (set <<Confuse4D>> to true in order to use)"};
  ConfigurableAxis bin4DmultPercentile{"bin4DmultPercentile", {10, 0.0f, 100.0f}, "multiplicity percentile Binning for the 4Dimensional plot: k* vs multiplicity vs multiplicity percentile vs mT (set <<Confuse4D>> to true in order to use)"};
  ConfigurableAxis binInvMass{"binInvMass", {400, 2.10, 2.50}, "InvMass binning"};
  ConfigurableAxis binpTCharm{"binpTCharm", {360, 0, 36}, "pT binning of charm hadron"};
  ConfigurableAxis binTempFitVarTrack{"binTempFitVarTrack", {300, -0.15, 0.15}, "binning of the TempFitVar in the pT vs. TempFitVar plot (Track)"};
  ConfigurableAxis binmT{"binmT", {225, 0., 7.5}, "binning mT"};
  ConfigurableAxis binmultTempFit{"binmultTempFit", {1, 0, 1}, "multiplicity Binning for the TempFitVar plot"};
  ConfigurableAxis binMulPercentile{"binMulPercentile", {10, 0.0f, 100.0f}, "multiplicity percentile Binning"};
  ConfigurableAxis binpTTrack{"binpTTrack", {50, 0.5, 10.05}, "pT binning of the pT vs. TempFitVar plot (Track)"};
  ConfigurableAxis binEta{"binEta", {{200, -1.5, 1.5}}, "eta binning"};
  ConfigurableAxis binPhi{"binPhi", {{200, 0, o2::constants::math::TwoPI}}, "phi binning"};
  ConfigurableAxis binkT{"binkT", {150, 0., 9.}, "binning kT"};
  ConfigurableAxis binkstar{"binkstar", {1500, 0., 6.}, "binning kstar"};
  ConfigurableAxis binNSigmaTPC{"binNSigmaTPC", {1600, -8, 8}, "Binning of Nsigma TPC plot"};
  ConfigurableAxis binNSigmaTOF{"binNSigmaTOF", {3000, -15, 15}, "Binning of the Nsigma TOF plot"};
  ConfigurableAxis binNSigmaTPCTOF{"binNSigmaTPCTOF", {3000, -15, 15}, "Binning of the Nsigma TPC+TOF plot"};
  ConfigurableAxis binTPCClusters{"binTPCClusters", {163, -0.5, 162.5}, "Binning of TPC found clusters plot"};
  // Mixing axis configurables
  ConfigurableAxis mixingBinMult{"mixingBinMult", {VARIABLE_WIDTH, 0.0f, 20.0f, 60.0f, 200.0f}, "Mixing bins - multiplicity"};
  ConfigurableAxis mixingBinMultPercentile{"mixingBinMultPercentile", {VARIABLE_WIDTH, 0.0f, 100.f}, "Mixing bins - multiplicity percentile"};
  ConfigurableAxis mixingBinVztx{"mixingBinVztx", {VARIABLE_WIDTH, -10.0f, -4.f, 0.f, 4.f, 10.f}, "Mixing bins - z-vertex"};

  ColumnBinningPolicy<aod::collision::PosZ, aod::femtodreamcollision::MultNtr> colBinningMult{{mixingBinVztx, mixingBinMult}, true};
  ColumnBinningPolicy<aod::collision::PosZ, aod::femtodreamcollision::MultV0M> colBinningMultPercentile{{mixingBinVztx, mixingBinMultPercentile}, true};
  ColumnBinningPolicy<aod::collision::PosZ, aod::femtodreamcollision::MultNtr, aod::femtodreamcollision::MultV0M> colBinningMultMultPercentile{{mixingBinVztx, mixingBinMult, mixingBinMultPercentile}, true};

  FemtoDreamContainer<femtoDreamContainer::EventType::same, femtoDreamContainer::Observable::kstar> sameEventCont;
  FemtoDreamContainer<femtoDreamContainer::EventType::mixed, femtoDreamContainer::Observable::kstar> mixedEventCont;
  FemtoDreamPairCleaner<aod::femtodreamparticle::ParticleType::kTrack, aod::femtodreamparticle::ParticleType::kCharmHadron3Prong> pairCleaner3Prong;
  FemtoDreamDetaDphiStar<aod::femtodreamparticle::ParticleType::kTrack, aod::femtodreamparticle::ParticleType::kCharmHadron3Prong> pairCloseRejectionSE3Prong;
  FemtoDreamDetaDphiStar<aod::femtodreamparticle::ParticleType::kTrack, aod::femtodreamparticle::ParticleType::kCharmHadron3Prong> pairCloseRejectionME3Prong;

  FemtoDreamPairCleaner<aod::femtodreamparticle::ParticleType::kTrack, aod::femtodreamparticle::ParticleType::kCharmHadron2Prong> pairCleaner2Prong;
  FemtoDreamDetaDphiStar<aod::femtodreamparticle::ParticleType::kTrack, aod::femtodreamparticle::ParticleType::kCharmHadron2Prong> pairCloseRejectionSE2Prong;
  FemtoDreamDetaDphiStar<aod::femtodreamparticle::ParticleType::kTrack, aod::femtodreamparticle::ParticleType::kCharmHadron2Prong> pairCloseRejectionME2Prong;

  FemtoDreamPairCleaner<aod::femtodreamparticle::ParticleType::kTrack, aod::femtodreamparticle::ParticleType::kCharmHadronDstar> pairCleanerDstar;
  FemtoDreamDetaDphiStar<aod::femtodreamparticle::ParticleType::kTrack, aod::femtodreamparticle::ParticleType::kCharmHadronDstar> pairCloseRejectionSEDstar;
  FemtoDreamDetaDphiStar<aod::femtodreamparticle::ParticleType::kTrack, aod::femtodreamparticle::ParticleType::kCharmHadronDstar> pairCloseRejectionMEDstar;

  femtodreamcollision::BitMaskType bitMask = 1 << 0;

  /// Histogramming for particle 1
  FemtoDreamParticleHisto<aod::femtodreamparticle::ParticleType::kTrack, 1> allTrackHisto;
  FemtoDreamParticleHisto<aod::femtodreamparticle::ParticleType::kTrack, 5> selectedTrackHisto;

  /// Histogramming for Event
  FemtoDreamEventHisto eventHisto;
  /// Histogram output
  HistogramRegistry registry{"CorrelationsAndQA", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry registryMixQa{"registryMixQa"};
  HistogramRegistry registryCharmHadronQa{"registryCharmHadronQa"};

  float massOne = o2::analysis::femtoDream::getMass(trackSel.pdgCodeTrack1);
  float massTwo = o2::analysis::femtoDream::getMass(charmSel.charmHadPDGCode);
  int8_t partSign = 0;
  int64_t processType = 0;

  void init(InitContext& /*context*/)
  {
    std::array<bool, 8> processes = {doprocessDataLcTrk, doprocessDataDplusTrk, doprocessDataD0Trk, doprocessDataDstarTrk, doprocessMcLcTrk, doprocessMcDplusTrk, doprocessMcD0Trk, doprocessMcDstarTrk};
    if (std::accumulate(processes.begin(), processes.end(), 0) != 1) {
      LOGP(fatal, "One and only one process function must be enabled at a time.");
    }
    bool process3Prong = doprocessDataLcTrk || doprocessDataDplusTrk || doprocessMcLcTrk || doprocessMcDplusTrk;
    bool process2Prong = doprocessDataD0Trk || doprocessMcD0Trk;
    bool processDstar = doprocessDataDstarTrk || doprocessMcDstarTrk;

    // setup columnpolicy for binning
    colBinningMult = {{mixingBinVztx, mixingBinMult}, true};
    colBinningMultPercentile = {{mixingBinVztx, mixingBinMultPercentile}, true};
    colBinningMultMultPercentile = {{mixingBinVztx, mixingBinMult, mixingBinMultPercentile}, true};
    eventHisto.init(&registry);
    allTrackHisto.init(&registry, binmultTempFit, binMulPercentile, binpTTrack, binEta, binPhi, binTempFitVarTrack, binNSigmaTPC, binNSigmaTOF, binNSigmaTPCTOF, binTPCClusters, dummy, dummy, isMc, trackSel.pdgCodeTrack1, true);
    selectedTrackHisto.init(&registry, binmultTempFit, binMulPercentile, binpTTrack, binEta, binPhi, binTempFitVarTrack, binNSigmaTPC, binNSigmaTOF, binNSigmaTPCTOF, binTPCClusters, dummy, dummy, isMc, trackSel.pdgCodeTrack1, true);

    sameEventCont.init(&registry,
                       binkstar, binpTTrack, binkT, binmT, mixingBinMult, mixingBinMultPercentile,
                       bin4Dkstar, bin4DmT, bin4DMult, bin4DmultPercentile,
                       isMc, use4D, extendedPlots,
                       highkstarCut,
                       smearingByOrigin, binInvMass);

    sameEventCont.setPDGCodes(trackSel.pdgCodeTrack1, charmSel.charmHadPDGCode);
    mixedEventCont.init(&registry,
                        binkstar, binpTTrack, binkT, binmT, mixingBinMult, mixingBinMultPercentile,
                        bin4Dkstar, bin4DmT, bin4DMult, bin4DmultPercentile,
                        isMc, use4D, extendedPlots,
                        highkstarCut,
                        smearingByOrigin, binInvMass);

    mixedEventCont.setPDGCodes(trackSel.pdgCodeTrack1, charmSel.charmHadPDGCode);
    registryMixQa.add("MixingQA/hSECollisionBins", "; bin; Entries", kTH1F, {{120, -0.5, 119.5}});
    registryMixQa.add("MixingQA/hSECollisionPool", "; Vz (cm); Mul", kTH2F, {{100, -10, 10}, {200, 0, 200}});
    registryMixQa.add("MixingQA/hMECollisionBins", "; bin; Entries", kTH1F, {{120, -0.5, 119.5}});
    registryCharmHadronQa.add("CharmHadronQA/hPtVsMass", "; #it{p}_{T} (GeV/#it{c}); inv. mass (GeV/#it{c}^{2})", kTH2F, {binpTCharm, binInvMass});

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

  template <DecayChannel Channel, typename Candidate, typename Track>
  float getCharmHadronTrackMass(const Candidate& cand,
                                const Track& trk,
                                int trkPDGCode)
  {
    auto pVecProng0 = RecoDecayPtEtaPhi::pVector(cand.prong0Pt(), cand.prong0Eta(), cand.prong0Phi());
    auto pVecProng1 = RecoDecayPtEtaPhi::pVector(cand.prong1Pt(), cand.prong1Eta(), cand.prong1Phi());
    auto pVecTrack = RecoDecayPtEtaPhi::pVector(trk.pt(), trk.eta(), trk.phi());

    double trackMassHyp = 0.;
    switch (trkPDGCode) {
      case kProton:
        trackMassHyp = MassProton;
        break;
      case kPiPlus:
        trackMassHyp = MassPiPlus;
        break;
      case kKPlus:
        trackMassHyp = MassKPlus;
        break;
      case kDeuteron:
        trackMassHyp = MassDeuteron;
        break;
      default:
        LOG(fatal) << "Invalid PDG code for track mass hypothesis: " << trkPDGCode;
    }

    // D0 → K π + track (2-prong)
    if constexpr (Channel == DecayChannel::D0ToPiK) {
      const auto pVecCharmTrk = std::array{pVecProng0, pVecProng1, pVecTrack};
      std::array<double, 3> massCharmTrk{};

      if (cand.candidateSelFlag() == 1) {
        massCharmTrk = {MassPiPlus, MassKPlus, trackMassHyp};
      } else {
        massCharmTrk = {MassKPlus, MassPiPlus, trackMassHyp};
      }

      return static_cast<float>(RecoDecay::m(pVecCharmTrk, massCharmTrk));
    }

    // 3-prong：Λc → p K π, D+ → π K π + track
    if constexpr (Channel == DecayChannel::LcToPKPi || Channel == DecayChannel::DplusToPiKPi || Channel == DecayChannel::DstarToD0Pi) {
      auto pVecProng2 = RecoDecayPtEtaPhi::pVector(cand.prong2Pt(), cand.prong2Eta(), cand.prong2Phi());
      const auto pVecCharmTrk = std::array{pVecProng0, pVecProng1, pVecProng2, pVecTrack};
      std::array<double, 4> massCharmTrk{};

      if constexpr (Channel == DecayChannel::LcToPKPi) {
        // Λc⁺ → p K π
        if (cand.candidateSelFlag() == 1) {
          massCharmTrk = {MassProton, MassKPlus, MassPiPlus, trackMassHyp};
        } else {
          massCharmTrk = {MassPiPlus, MassKPlus, MassProton, trackMassHyp};
        }
      } else if constexpr (Channel == DecayChannel::DplusToPiKPi) {
        // D⁺ → π K π
        massCharmTrk = {MassPiPlus, MassKPlus, MassPiPlus, trackMassHyp};
      } else if constexpr (Channel == DecayChannel::DstarToD0Pi) {
        // D* → D0π
        if (cand.candidateSelFlag() == 1) {
          massCharmTrk = {MassPiPlus, MassKPlus, MassPiPlus, trackMassHyp};
        } else {
          massCharmTrk = {MassKPlus, MassPiPlus, MassPiPlus, trackMassHyp};
        }
      }

      return static_cast<float>(RecoDecay::m(pVecCharmTrk, massCharmTrk));
    }
    return -1.f;
  }

  /// This function processes the same event and takes care of all the histogramming
  template <bool IsMc, DecayChannel Channel, typename CandType, typename PartitionType, typename TableTracks, typename Collision>
  void doSameEvent(CandType& sliceCharmHad, PartitionType& sliceTrk1, TableTracks const& parts, Collision const& col)
  {
    fillCollision(col);
    processType = 1; // for same event
    for (auto const& [p1, p2] : combinations(CombinationsFullIndexPolicy(sliceTrk1, sliceCharmHad))) {

      if constexpr (Channel == DecayChannel::D0ToPiK) {
        if (p1.trackId() == p2.prong0Id() || p1.trackId() == p2.prong1Id())
          continue;

        if (useCPR.value) {
          if (pairCloseRejectionSE2Prong.isClosePair(p1, p2, parts, col.magField())) {
            continue;
          }
        }

        if (!pairCleaner2Prong.isCleanPair(p1, p2, parts)) {
          continue;
        }
      }

      if constexpr (Channel == DecayChannel::LcToPKPi || Channel == DecayChannel::DplusToPiKPi) {
        if (p1.trackId() == p2.prong0Id() || p1.trackId() == p2.prong1Id() || p1.trackId() == p2.prong2Id())
          continue;
        if (useCPR.value) {
          if (pairCloseRejectionSE3Prong.isClosePair(p1, p2, parts, col.magField())) {
            continue;
          }
        }

        if (!pairCleaner3Prong.isCleanPair(p1, p2, parts)) {
          continue;
        }
      }

      if constexpr (Channel == DecayChannel::DstarToD0Pi) {
        if (p1.trackId() == p2.prong0Id() || p1.trackId() == p2.prong1Id() || p1.trackId() == p2.prong2Id())
          continue;
        if (useCPR.value) {
          if (pairCloseRejectionSEDstar.isClosePair(p1, p2, parts, col.magField())) {
            continue;
          }
        }

        if (!pairCleanerDstar.isCleanPair(p1, p2, parts)) {
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

      float deltaInvMassPair = getCharmHadronTrackMass<Channel>(p2, p1, trackSel.pdgCodeTrack1.value) - invMass;

      // proton track charge
      float chargeTrack = 0.;
      if ((p1.cut() & CutBitChargePositive) == CutBitChargePositive) {
        chargeTrack = PositiveCharge;
      } else {
        chargeTrack = NegativeCharge;
      }
      int pairSign = 0;
      if (chargeTrack == p2.charge()) {
        pairSign = LikeSignPair;
      } else if (chargeTrack == -p2.charge()) {
        pairSign = UnLikeSignPair;
      } else {
        pairSign = ReflectedPair;
      }

      /// Filling QA histograms of the selected tracks
      selectedTrackHisto.fillQA<IsMc, true>(p1, static_cast<aod::femtodreamparticle::MomentumType>(confTempFitVarMomentum.value), col.multNtr(), col.multV0M());

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
        deltaInvMassPair,
        processType,
        charmHadMc,
        originType);

      sameEventCont.setPair<IsMc, true>(p1, p2, col.multNtr(), col.multV0M(), use4D, extendedPlots, smearingByOrigin);
    }
  }

  template <bool IsMc, DecayChannel Channel, typename CollisionType, typename PartitionType1, typename PartitionType2, typename TableTracks, typename BinningType>
  void doMixedEvent(CollisionType const& cols, PartitionType1& charms, PartitionType2& trks, TableTracks const& parts, BinningType policy)
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

      auto sliceTrk1 = trks->sliceByCached(aod::femtodreamparticle::fdCollisionId, collision1.globalIndex(), cache);
      auto sliceCharmHad = charms->sliceByCached(aod::femtodreamparticle::fdCollisionId, collision2.globalIndex(), cache);
      for (const auto& [p1, p2] : combinations(CombinationsFullIndexPolicy(sliceTrk1, sliceCharmHad))) {

        if constexpr (Channel == DecayChannel::D0ToPiK) {

          if (useCPR.value) {
            if (pairCloseRejectionME2Prong.isClosePair(p1, p2, parts, collision1.magField())) {
              continue;
            }
          }

          if (!pairCleaner2Prong.isCleanPair(p1, p2, parts)) {
            continue;
          }
        }

        if constexpr (Channel == DecayChannel::DplusToPiKPi || Channel == DecayChannel::LcToPKPi) {

          if (useCPR.value) {
            if (pairCloseRejectionME3Prong.isClosePair(p1, p2, parts, collision1.magField())) {
              continue;
            }
          }

          if (!pairCleaner3Prong.isCleanPair(p1, p2, parts)) {
            continue;
          }
        }

        if constexpr (Channel == DecayChannel::DstarToD0Pi) {

          if (useCPR.value) {
            if (pairCloseRejectionME3Prong.isClosePair(p1, p2, parts, collision1.magField())) {
              continue;
            }
          }

          if (!pairCleanerDstar.isCleanPair(p1, p2, parts)) {
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

        float deltaInvMassPair = getCharmHadronTrackMass<Channel>(p2, p1, trackSel.pdgCodeTrack1.value) - invMass;

        // proton track charge
        float chargeTrack = 0.;
        if ((p1.cut() & CutBitChargePositive) == CutBitChargePositive) {
          chargeTrack = PositiveCharge;
        } else {
          chargeTrack = NegativeCharge;
        }
        int pairSign = 0;
        if (chargeTrack == p2.charge()) {
          pairSign = LikeSignPair;
        } else if (chargeTrack == -p2.charge()) {
          pairSign = UnLikeSignPair;
        } else {
          pairSign = ReflectedPair;
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
          deltaInvMassPair,
          processType,
          charmHadMc,
          originType);

        mixedEventCont.setPair<IsMc, true>(p1, p2, collision1.multNtr(), collision1.multV0M(), use4D, extendedPlots, smearingByOrigin);
      }
    }
  }
  template <bool IsMc, DecayChannel Channel, typename TrackType, typename CandType, typename CollType>
  void fillTables(const CollType& col,
                  const TrackType& sliceTrk1,
                  const CandType& sliceCharmHad)
  {
    int64_t timeStamp = -999;

    // ---- Fill Charm-Hadron Table ----
    for (auto const& part : sliceCharmHad) {
      float invMass = getCharmHadronMass<Channel>(part);
      registryCharmHadronQa.fill(
        HIST("CharmHadronQA/hPtVsMass"),
        part.pt(), invMass);

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
    // ---- Fill Track Table ----
    for (auto const& part : sliceTrk1) {
      allTrackHisto.fillQA<IsMc, true>(
        part,
        static_cast<aod::femtodreamparticle::MomentumType>(confTempFitVarMomentum.value),
        col.multNtr(),
        col.multV0M());

      float chargeTrack = ((part.cut() & CutBitChargePositive) == CutBitChargePositive)
                            ? PositiveCharge
                            : NegativeCharge;

      timeStamp = part.timeStamp();
      float tpcNSigma = 999.f;
      float tofNSigma = 999.f;
      switch (trackSel.pdgCodeTrack1.value) {
        case kProton:
          tpcNSigma = part.tpcNSigmaPr();
          tofNSigma = part.tofNSigmaPr();
          break;
        case kPiPlus:
          tpcNSigma = part.tpcNSigmaPi();
          tofNSigma = part.tofNSigmaPi();
          break;
        case kKPlus:
          tpcNSigma = part.tpcNSigmaKa();
          tofNSigma = part.tofNSigmaKa();
          break;
        case kDeuteron:
          tpcNSigma = part.tpcNSigmaDe();
          tofNSigma = part.tofNSigmaDe();
          break;
        default:
          LOG(fatal) << "Unhandled PDG code in PID switch: "
                     << trackSel.pdgCodeTrack1.value;
          break;
      }
      rowFemtoResultTrk(
        col.globalIndex(),
        timeStamp,
        part.pt(),
        part.eta(),
        part.phi(),
        part.trackId(),
        chargeTrack,
        part.tpcNClsFound(),
        part.tpcNClsFindable(),
        part.tpcNClsCrossedRows(),
        tpcNSigma,
        tofNSigma);
    }

    // ---- Fill Collision Table ----
    if (sliceCharmHad.size() > 0 || sliceTrk1.size() > 0) {
      rowFemtoResultColl(
        col.globalIndex(),
        timeStamp,
        col.posZ(),
        col.multNtr());
    }
  }

  void processDataLcTrk(FilteredCollisions const& cols,
                        FilteredFDParticles const& parts,
                        FilteredCharmCand3Prongs const&)
  {
    for (const auto& col : cols) {
      eventHisto.fillQA(col);
      auto sliceTrk1 = partitionTrk1->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
      auto sliceCharmHad = partitionCharmHadron3Prong->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
      if (fillTableWithCharm.value && sliceCharmHad.size() == 0) {
        continue;
      } else {
        fillTables<false, DecayChannel::LcToPKPi>(col, sliceTrk1, sliceCharmHad);
      }
      if (sliceCharmHad.size() > 0 && sliceTrk1.size() > 0) {
        doSameEvent<false, DecayChannel::LcToPKPi, FilteredCharmCand3Prongs>(sliceCharmHad, sliceTrk1, parts, col);
      }
    }
    if (mixSetting.doMixEvent) {
      switch (mixSetting.mixingBinPolicy) {
        case femtodreamcollision::kMult:
          doMixedEvent<false, DecayChannel::LcToPKPi, FilteredCollisions>(cols, partitionCharmHadron3Prong, partitionTrk1, parts, colBinningMult);
          break;
        case femtodreamcollision::kMultPercentile:
          doMixedEvent<false, DecayChannel::LcToPKPi, FilteredCollisions>(cols, partitionCharmHadron3Prong, partitionTrk1, parts, colBinningMultPercentile);
          break;
        case femtodreamcollision::kMultMultPercentile:
          doMixedEvent<false, DecayChannel::LcToPKPi, FilteredCollisions>(cols, partitionCharmHadron3Prong, partitionTrk1, parts, colBinningMultMultPercentile);
          break;
        default:
          LOG(fatal) << "Invalid binning policiy specifed. Breaking...";
      }
    }
  }
  PROCESS_SWITCH(HfTaskCharmHadronsTrackFemtoDream, processDataLcTrk, "Enable processing LcToPKPi and Tracks correlation", false);

  void processDataDplusTrk(FilteredCollisions const& cols,
                           FilteredFDParticles const& parts,
                           FilteredCharmCand3Prongs const&)
  {
    for (const auto& col : cols) {
      eventHisto.fillQA(col);
      auto sliceTrk1 = partitionTrk1->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
      auto sliceCharmHad = partitionCharmHadron3Prong->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);

      if (fillTableWithCharm.value && sliceCharmHad.size() == 0) {
        continue;
      } else {
        fillTables<false, DecayChannel::DplusToPiKPi>(col, sliceTrk1, sliceCharmHad);
      }
      if (sliceCharmHad.size() > 0 && sliceTrk1.size() > 0) {
        doSameEvent<false, DecayChannel::DplusToPiKPi, FilteredCharmCand3Prongs>(sliceCharmHad, sliceTrk1, parts, col);
      }
    }
    if (mixSetting.doMixEvent) {
      switch (mixSetting.mixingBinPolicy) {
        case femtodreamcollision::kMult:
          doMixedEvent<false, DecayChannel::DplusToPiKPi, FilteredCollisions>(cols, partitionCharmHadron3Prong, partitionTrk1, parts, colBinningMult);
          break;
        case femtodreamcollision::kMultPercentile:
          doMixedEvent<false, DecayChannel::DplusToPiKPi, FilteredCollisions>(cols, partitionCharmHadron3Prong, partitionTrk1, parts, colBinningMultPercentile);
          break;
        case femtodreamcollision::kMultMultPercentile:
          doMixedEvent<false, DecayChannel::DplusToPiKPi, FilteredCollisions>(cols, partitionCharmHadron3Prong, partitionTrk1, parts, colBinningMultMultPercentile);
          break;
        default:
          LOG(fatal) << "Invalid binning policiy specifed. Breaking...";
      }
    }
  }
  PROCESS_SWITCH(HfTaskCharmHadronsTrackFemtoDream, processDataDplusTrk, "Enable processing DplusToPiKPi and Tracks correlation", false);

  void processDataD0Trk(FilteredCollisions const& cols,
                        FilteredFDParticles const& parts,
                        FilteredCharmCand2Prongs const&)
  {
    for (const auto& col : cols) {
      eventHisto.fillQA(col);
      auto sliceTrk1 = partitionTrk1->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
      auto sliceCharmHad = partitionCharmHadron2Prong->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
      if (fillTableWithCharm.value && sliceCharmHad.size() == 0) {
        continue;
      } else {
        fillTables<false, DecayChannel::D0ToPiK>(col, sliceTrk1, sliceCharmHad);
      }
      if (sliceCharmHad.size() > 0 && sliceTrk1.size() > 0) {
        doSameEvent<false, DecayChannel::D0ToPiK, FilteredCharmCand2Prongs>(sliceCharmHad, sliceTrk1, parts, col);
      }
    }
    if (mixSetting.doMixEvent) {
      switch (mixSetting.mixingBinPolicy) {
        case femtodreamcollision::kMult:
          doMixedEvent<false, DecayChannel::D0ToPiK, FilteredCollisions>(cols, partitionCharmHadron2Prong, partitionTrk1, parts, colBinningMult);
          break;
        case femtodreamcollision::kMultPercentile:
          doMixedEvent<false, DecayChannel::D0ToPiK, FilteredCollisions>(cols, partitionCharmHadron2Prong, partitionTrk1, parts, colBinningMultPercentile);
          break;
        case femtodreamcollision::kMultMultPercentile:
          doMixedEvent<false, DecayChannel::D0ToPiK, FilteredCollisions>(cols, partitionCharmHadron2Prong, partitionTrk1, parts, colBinningMultMultPercentile);
          break;
        default:
          LOG(fatal) << "Invalid binning policiy specifed. Breaking...";
      }
    }
  }
  PROCESS_SWITCH(HfTaskCharmHadronsTrackFemtoDream, processDataD0Trk, "Enable processing D0ToPiK and Tracks correlation", false);

  void processDataDstarTrk(FilteredCollisions const& cols,
                           FilteredFDParticles const& parts,
                           FilteredCharmCandDstars const&)
  {
    for (const auto& col : cols) {
      eventHisto.fillQA(col);
      auto sliceTrk1 = partitionTrk1->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
      auto sliceCharmHad = partitionCharmHadronDstar->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
      if (fillTableWithCharm.value && sliceCharmHad.size() == 0) {
        continue;
      } else {
        fillTables<false, DecayChannel::DstarToD0Pi>(col, sliceTrk1, sliceCharmHad);
      }
      if (sliceCharmHad.size() > 0 && sliceTrk1.size() > 0) {
        doSameEvent<false, DecayChannel::DstarToD0Pi, FilteredCharmCandDstars>(sliceCharmHad, sliceTrk1, parts, col);
      }
    }
    if (mixSetting.doMixEvent) {
      switch (mixSetting.mixingBinPolicy) {
        case femtodreamcollision::kMult:
          doMixedEvent<false, DecayChannel::DstarToD0Pi, FilteredCollisions>(cols, partitionCharmHadronDstar, partitionTrk1, parts, colBinningMult);
          break;
        case femtodreamcollision::kMultPercentile:
          doMixedEvent<false, DecayChannel::DstarToD0Pi, FilteredCollisions>(cols, partitionCharmHadronDstar, partitionTrk1, parts, colBinningMultPercentile);
          break;
        case femtodreamcollision::kMultMultPercentile:
          doMixedEvent<false, DecayChannel::DstarToD0Pi, FilteredCollisions>(cols, partitionCharmHadronDstar, partitionTrk1, parts, colBinningMultMultPercentile);
          break;
        default:
          LOG(fatal) << "Invalid binning policiy specifed. Breaking...";
      }
    }
  }
  PROCESS_SWITCH(HfTaskCharmHadronsTrackFemtoDream, processDataDstarTrk, "Enable processing DstarToD0Pi and Tracks correlation", false);

  void processMcLcTrk(FilteredMcColisions const& cols,
                      FilteredFDMcParts const& parts,
                      o2::aod::FDMCParticles const&,
                      o2::aod::FDExtMCParticles const&,
                      FilteredCharmMcCand3Prongs const&)
  {
    for (const auto& col : cols) {
      eventHisto.fillQA(col);
      auto sliceMcTrk1 = partitionMcTrk1->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
      auto sliceMcCharmHad = partitionMcCharmHadron3Prong->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
      if ((col.bitmaskTrackOne() & bitMask) != bitMask || (col.bitmaskTrackTwo() & bitMask) != bitMask) {
        continue;
      }
      doSameEvent<true, DecayChannel::LcToPKPi, FilteredCharmMcCand3Prongs>(sliceMcCharmHad, sliceMcTrk1, parts, col);
    }
    switch (mixSetting.mixingBinPolicy) {
      case femtodreamcollision::kMult:
        doMixedEvent<true, DecayChannel::LcToPKPi, FilteredMcColisions>(cols, partitionMcCharmHadron3Prong, partitionMcTrk1, parts, colBinningMult);
        break;
      case femtodreamcollision::kMultPercentile:
        doMixedEvent<true, DecayChannel::LcToPKPi, FilteredMcColisions>(cols, partitionMcCharmHadron3Prong, partitionMcTrk1, parts, colBinningMultPercentile);
        break;
      case femtodreamcollision::kMultMultPercentile:
        doMixedEvent<true, DecayChannel::LcToPKPi, FilteredMcColisions>(cols, partitionMcCharmHadron3Prong, partitionMcTrk1, parts, colBinningMultMultPercentile);
        break;
      default:
        LOG(fatal) << "Invalid binning policiy specifed. Breaking...";
    }
  }
  PROCESS_SWITCH(HfTaskCharmHadronsTrackFemtoDream, processMcLcTrk, "Enable processing LcToPKPi and Tracks correlation for Monte Carlo", false);

  void processMcDplusTrk(FilteredMcColisions const& cols,
                         FilteredFDMcParts const& parts,
                         o2::aod::FDMCParticles const&,
                         o2::aod::FDExtMCParticles const&,
                         FilteredCharmMcCand3Prongs const&)
  {
    for (const auto& col : cols) {
      eventHisto.fillQA(col);
      auto sliceMcTrk1 = partitionMcTrk1->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
      auto sliceMcCharmHad = partitionMcCharmHadron3Prong->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
      if ((col.bitmaskTrackOne() & bitMask) != bitMask || (col.bitmaskTrackTwo() & bitMask) != bitMask) {
        continue;
      }
      doSameEvent<true, DecayChannel::DplusToPiKPi, FilteredCharmMcCand3Prongs>(sliceMcCharmHad, sliceMcTrk1, parts, col);
    }
    switch (mixSetting.mixingBinPolicy) {
      case femtodreamcollision::kMult:
        doMixedEvent<true, DecayChannel::DplusToPiKPi, FilteredMcColisions>(cols, partitionMcCharmHadron3Prong, partitionMcTrk1, parts, colBinningMult);
        break;
      case femtodreamcollision::kMultPercentile:
        doMixedEvent<true, DecayChannel::DplusToPiKPi, FilteredMcColisions>(cols, partitionMcCharmHadron3Prong, partitionMcTrk1, parts, colBinningMultPercentile);
        break;
      case femtodreamcollision::kMultMultPercentile:
        doMixedEvent<true, DecayChannel::DplusToPiKPi, FilteredMcColisions>(cols, partitionMcCharmHadron3Prong, partitionMcTrk1, parts, colBinningMultMultPercentile);
        break;
      default:
        LOG(fatal) << "Invalid binning policiy specifed. Breaking...";
    }
  }
  PROCESS_SWITCH(HfTaskCharmHadronsTrackFemtoDream, processMcDplusTrk, "Enable processing DplusToPiKPi and Tracks correlation for Monte Carlo", false);

  void processMcD0Trk(FilteredMcColisions const& cols,
                      FilteredFDMcParts const& parts,
                      o2::aod::FDMCParticles const&,
                      o2::aod::FDExtMCParticles const&,
                      FilteredCharmMcCand2Prongs const&)
  {
    for (const auto& col : cols) {
      eventHisto.fillQA(col);
      auto sliceMcTrk1 = partitionMcTrk1->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
      auto sliceMcCharmHad = partitionMcCharmHadron2Prong->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
      if ((col.bitmaskTrackOne() & bitMask) != bitMask || (col.bitmaskTrackTwo() & bitMask) != bitMask) {
        continue;
      }
      doSameEvent<true, DecayChannel::D0ToPiK, FilteredCharmMcCand2Prongs>(sliceMcCharmHad, sliceMcTrk1, parts, col);
    }
    switch (mixSetting.mixingBinPolicy) {
      case femtodreamcollision::kMult:
        doMixedEvent<true, DecayChannel::D0ToPiK, FilteredMcColisions>(cols, partitionMcCharmHadron2Prong, partitionMcTrk1, parts, colBinningMult);
        break;
      case femtodreamcollision::kMultPercentile:
        doMixedEvent<true, DecayChannel::D0ToPiK, FilteredMcColisions>(cols, partitionMcCharmHadron2Prong, partitionMcTrk1, parts, colBinningMultPercentile);
        break;
      case femtodreamcollision::kMultMultPercentile:
        doMixedEvent<true, DecayChannel::D0ToPiK, FilteredMcColisions>(cols, partitionMcCharmHadron2Prong, partitionMcTrk1, parts, colBinningMultMultPercentile);
        break;
      default:
        LOG(fatal) << "Invalid binning policiy specifed. Breaking...";
    }
  }
  PROCESS_SWITCH(HfTaskCharmHadronsTrackFemtoDream, processMcD0Trk, "Enable processing D0ToPiK and Tracks correlation for Monte Carlo", false);

  void processMcDstarTrk(FilteredMcColisions const& cols,
                         FilteredFDMcParts const& parts,
                         o2::aod::FDMCParticles const&,
                         o2::aod::FDExtMCParticles const&,
                         FilteredCharmMcCandDstars const&)
  {
    for (const auto& col : cols) {
      eventHisto.fillQA(col);
      auto sliceMcTrk1 = partitionMcTrk1->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
      auto sliceMcCharmHad = partitionMcCharmHadronDstar->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
      if ((col.bitmaskTrackOne() & bitMask) != bitMask || (col.bitmaskTrackTwo() & bitMask) != bitMask) {
        continue;
      }
      doSameEvent<true, DecayChannel::DstarToD0Pi, FilteredCharmMcCandDstars>(sliceMcCharmHad, sliceMcTrk1, parts, col);
    }
    switch (mixSetting.mixingBinPolicy) {
      case femtodreamcollision::kMult:
        doMixedEvent<true, DecayChannel::DstarToD0Pi, FilteredMcColisions>(cols, partitionMcCharmHadronDstar, partitionMcTrk1, parts, colBinningMult);
        break;
      case femtodreamcollision::kMultPercentile:
        doMixedEvent<true, DecayChannel::DstarToD0Pi, FilteredMcColisions>(cols, partitionMcCharmHadronDstar, partitionMcTrk1, parts, colBinningMultPercentile);
        break;
      case femtodreamcollision::kMultMultPercentile:
        doMixedEvent<true, DecayChannel::DstarToD0Pi, FilteredMcColisions>(cols, partitionMcCharmHadronDstar, partitionMcTrk1, parts, colBinningMultMultPercentile);
        break;
      default:
        LOG(fatal) << "Invalid binning policiy specifed. Breaking...";
    }
  }
  PROCESS_SWITCH(HfTaskCharmHadronsTrackFemtoDream, processMcDstarTrk, "Enable processing DstarToD0Pi and Tracks correlation for Monte Carlo", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskCharmHadronsTrackFemtoDream>(cfgc)};
}
