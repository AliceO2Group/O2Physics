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

inline o2::framework::expressions::Node coshEta(o2::framework::expressions::Node&& eta)
{
  return (nexp(std::move(eta)) + nexp(0.0f - std::move(eta))) * 0.5f;
}

struct HfTaskCharmHadronsFemtoDream {

  enum TrackCharge {
    PositiveCharge = 1,
    NegativeCharge = -1
  };

  enum PairSign {
    PairNotDefined = 0,
    LikeSignPair = 1,
    UnLikeSignPair = 2
  };

  constexpr static int OriginRecPrompt = 1;
  constexpr static int OriginRecFD = 2;
  constexpr static int CutBitChargePositive = 2;

  Produces<o2::aod::FDHfPairs> rowFemtoResultPairs;
  Produces<o2::aod::FDHfCharm> rowFemtoResultCharm;
  Produces<o2::aod::FDHfTrk> rowFemtoResultTrk;
  Produces<o2::aod::FDHfColl> rowFemtoResultColl;

  Configurable<int> confTempFitVarMomentum{"confTempFitVarMomentum", 0, "Momentum used for binning: 0 -> pt; 1 -> preco; 2 -> ptpc"};

  /// Particle 2 (Charm Hadrons)
  Configurable<float> charmHadBkgBDTmax{"charmHadBkgBDTmax", 1., "Maximum background bdt score for Charm Hadron (particle 2)"};
  Configurable<int> charmHadCandSel{"charmHadCandSel", 1, "candidate selection for charm hadron"};
  Configurable<int> charmHadMcSel{"charmHadMcSel", DecayChannelMain::LcToPKPi, "charm hadron selection for mc, DplusToPiKPi = 1, LcToPKPi = 17"};
  Configurable<float> charmHadFdBDTmin{"charmHadFdBDTmin", 0., "Minimum feed-down bdt score Charm Hadron (particle 2)"};
  Configurable<float> charmHadFdBDTmax{"charmHadFdBDTmax", 1., "Maximum feed-down bdt score Charm Hadron (particle 2)"};
  Configurable<float> charmHadMaxInvMass{"charmHadMaxInvMass", 2.45, "Maximum invariant mass of Charm Hadron (particle 2)"};
  Configurable<float> charmHadMaxPt{"charmHadMaxPt", 999., "Maximum pT of Charm Hadron (particle 2)"};
  Configurable<float> charmHadMinInvMass{"charmHadMinInvMass", 2.15, "Minimum invariant mass of Charm Hadron (particle 2)"};
  Configurable<float> charmHadMinPt{"charmHadMinPt", 0., "Minimum pT of Charm Hadron (particle 2)"};
  Configurable<int> charmHadPDGCode{"charmHadPDGCode", 4122, "PDG code of particle 2 Charm Hadron"};
  Configurable<float> charmHadPromptBDTmin{"charmHadPromptBDTmin", 0., "Minimum prompt bdt score Charm Hadron (particle 2)"};
  Configurable<float> charmHadPromptBDTmax{"charmHadPromptBDTmax", 1., "Maximum prompt bdt score Charm Hadron (particle 2)"};

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
  Configurable<femtodreamparticle::cutContainerType> cutBitTrack1{"cutBitTrack1", 8188, "Particle 1 (Track) - Selection bit from cutCulator"};
  Configurable<int> pdgCodeTrack1{"pdgCodeTrack1", 2212, "PDG code of Particle 1 (Track)"};
  Configurable<float> pidThresTrack1{"pidThresTrack1", 0.75, "Momentum threshold for PID selection for particle 1 (Track)"};
  Configurable<femtodreamparticle::cutContainerType> tpcBitTrack1{"tpcBitTrack1", 4, "PID TPC bit from cutCulator for particle 1 (Track)"};
  Configurable<femtodreamparticle::cutContainerType> tpcTofBitTrack1{"tpcTofBitTrack1", 2, "PID TPCTOF bit from cutCulator for particle 1 (Track)"};
  Configurable<float> etaTrack1Max{"etaTrack1Max", 10., "Maximum eta of partricle 1 (Track)"};
  Configurable<float> ptTrack1Max{"ptTrack1Max", 999., "Maximum pT of partricle 1 (Track)"};
  Configurable<float> etaTrack1Min{"etaTrack1Min", -10., "Minimum eta of partricle 1 (Track)"};
  Configurable<float> ptTrack1Min{"ptTrack1Min", 0., "Minimum pT of partricle 1 (Track)"};

  SliceCache cache;

  using FilteredCharmCands = soa::Filtered<aod::FDHfCand>;
  using FilteredCharmCand = FilteredCharmCands::iterator;

  using FilteredCharmMcCands = soa::Filtered<soa::Join<aod::FDHfCand, aod::FDHfCandMC>>;
  using FilteredCharmMcCand = FilteredCharmMcCands::iterator;

  using FilteredColisions = soa::Filtered<soa::Join<FDCollisions, FDColMasks>>;
  using FilteredColision = FilteredColisions::iterator;

  using FilteredMcColisions = soa::Filtered<soa::Join<aod::FDCollisions, FDColMasks, aod::FDMCCollLabels>>;
  using FilteredMcColision = FilteredMcColisions::iterator;

  using FilteredFDMcParts = soa::Filtered<soa::Join<aod::FDParticles, aod::FDParticlesIndex, aod::FDExtParticles, aod::FDMCLabels, aod::FDExtMCLabels>>;
  using FilteredFDMcPart = FilteredFDMcParts::iterator;

  using FilteredFDParticles = soa::Filtered<soa::Join<aod::FDParticles, aod::FDExtParticles, aod::FDParticlesIndex, aod::FDTrkTimeStamp>>;
  using FilteredFDParticle = FilteredFDParticles::iterator;

  Filter eventMultiplicity = aod::femtodreamcollision::multNtr >= eventSel.multMin && aod::femtodreamcollision::multNtr <= eventSel.multMax;
  Filter eventMultiplicityPercentile = aod::femtodreamcollision::multV0M >= eventSel.multPercentileMin && aod::femtodreamcollision::multV0M <= eventSel.multPercentileMax;
  Filter hfCandSelFilter = aod::fdhf::candidateSelFlag >= charmHadCandSel;
  Filter hfMcSelFilter = (nabs(aod::fdhf::flagMc) == charmHadMcSel);
  Filter trackEtaFilterLow = ifnode(aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kTrack), aod::femtodreamparticle::eta < etaTrack1Max, true);
  Filter trackEtaFilterUp = ifnode(aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kTrack), aod::femtodreamparticle::eta > etaTrack1Min, true);
  Filter trackPtFilterLow = ifnode(aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kTrack), aod::femtodreamparticle::pt < ptTrack1Max, true);
  Filter trackPtFilterUp = ifnode(aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kTrack), aod::femtodreamparticle::pt > ptTrack1Min, true);

  Preslice<FilteredFDParticles> perCol = aod::femtodreamparticle::fdCollisionId;
  Preslice<FilteredCharmCands> perHfByCol = aod::femtodreamparticle::fdCollisionId;

  /// Partition for particle 1
  Partition<FilteredFDParticles> partitionTrk1 = (aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kTrack)) && (ncheckbit(aod::femtodreamparticle::cut, cutBitTrack1)) && ifnode(aod::femtodreamparticle::pt * coshEta(aod::femtodreamparticle::eta) <= pidThresTrack1, ncheckbit(aod::femtodreamparticle::pidcut, tpcBitTrack1), ncheckbit(aod::femtodreamparticle::pidcut, tpcTofBitTrack1));

  Partition<FilteredFDMcParts> partitionMcTrk1 = (aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kTrack)) &&
                                                 (ncheckbit(aod::femtodreamparticle::cut, cutBitTrack1)) &&
                                                 ifnode(aod::femtodreamparticle::pt * coshEta(aod::femtodreamparticle::eta) <= pidThresTrack1, ncheckbit(aod::femtodreamparticle::pidcut, tpcBitTrack1), ncheckbit(aod::femtodreamparticle::pidcut, tpcTofBitTrack1));

  /// Partition for particle 2
  Partition<FilteredCharmCands> partitionCharmHadron = aod::fdhf::bdtBkg < charmHadBkgBDTmax && aod::fdhf::bdtFD < charmHadFdBDTmax && aod::fdhf::bdtFD > charmHadFdBDTmin&& aod::fdhf::bdtPrompt<charmHadPromptBDTmax && aod::fdhf::bdtPrompt> charmHadPromptBDTmin;
  Partition<FilteredCharmMcCands> partitionMcCharmHadron = aod::fdhf::originMcRec == OriginRecPrompt || aod::fdhf::originMcRec == OriginRecFD;

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
  FemtoDreamPairCleaner<aod::femtodreamparticle::ParticleType::kTrack, aod::femtodreamparticle::ParticleType::kCharmHadron> pairCleaner;
  FemtoDreamDetaDphiStar<aod::femtodreamparticle::ParticleType::kTrack, aod::femtodreamparticle::ParticleType::kCharmHadron> pairCloseRejectionSE;
  FemtoDreamDetaDphiStar<aod::femtodreamparticle::ParticleType::kTrack, aod::femtodreamparticle::ParticleType::kCharmHadron> pairCloseRejectionME;

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

  float massOne = o2::analysis::femtoDream::getMass(pdgCodeTrack1);
  float massTwo = o2::analysis::femtoDream::getMass(charmHadPDGCode);
  int8_t partSign = 0;
  int64_t processType = 0;

  void init(InitContext& /*context*/)
  {
    // setup columnpolicy for binning
    colBinningMult = {{mixingBinVztx, mixingBinMult}, true};
    colBinningMultPercentile = {{mixingBinVztx, mixingBinMultPercentile}, true};
    colBinningMultMultPercentile = {{mixingBinVztx, mixingBinMult, mixingBinMultPercentile}, true};
    eventHisto.init(&registry);
    allTrackHisto.init(&registry, binmultTempFit, binMulPercentile, binpTTrack, binEta, binPhi, binTempFitVarTrack, binNSigmaTPC, binNSigmaTOF, binNSigmaTPCTOF, binTPCClusters, dummy, dummy, isMc, pdgCodeTrack1, true);
    selectedTrackHisto.init(&registry, binmultTempFit, binMulPercentile, binpTTrack, binEta, binPhi, binTempFitVarTrack, binNSigmaTPC, binNSigmaTOF, binNSigmaTPCTOF, binTPCClusters, dummy, dummy, isMc, pdgCodeTrack1, true);

    sameEventCont.init<true>(&registry,
                             binkstar, binpTTrack, binkT, binmT, mixingBinMult, mixingBinMultPercentile,
                             bin4Dkstar, bin4DmT, bin4DMult, bin4DmultPercentile,
                             isMc, use4D, extendedPlots,
                             highkstarCut,
                             smearingByOrigin, binInvMass);

    sameEventCont.setPDGCodes(pdgCodeTrack1, charmHadPDGCode);
    mixedEventCont.init<true>(&registry,
                              binkstar, binpTTrack, binkT, binmT, mixingBinMult, mixingBinMultPercentile,
                              bin4Dkstar, bin4DmT, bin4DMult, bin4DmultPercentile,
                              isMc, use4D, extendedPlots,
                              highkstarCut,
                              smearingByOrigin, binInvMass);

    mixedEventCont.setPDGCodes(pdgCodeTrack1, charmHadPDGCode);
    registryMixQa.add("MixingQA/hSECollisionBins", "; bin; Entries", kTH1F, {{120, -0.5, 119.5}});
    registryMixQa.add("MixingQA/hSECollisionPool", "; Vz (cm); Mul", kTH2F, {{100, -10, 10}, {200, 0, 200}});
    registryMixQa.add("MixingQA/hMECollisionBins", "; bin; Entries", kTH1F, {{120, -0.5, 119.5}});
    registryCharmHadronQa.add("CharmHadronQA/hPtVsMass", "; #it{p}_{T} (GeV/#it{c}); inv. mass (GeV/#it{c}^{2})", kTH2F, {binpTCharm, binInvMass});

    pairCleaner.init(&registry);
    if (useCPR.value) {
      pairCloseRejectionSE.init(&registry, &registry, cprDeltaPhiMax.value, cprDeltaEtaMax.value, cprPlotPerRadii.value, 1);
      pairCloseRejectionME.init(&registry, &registry, cprDeltaPhiMax.value, cprDeltaEtaMax.value, cprPlotPerRadii.value, 2);
    }
  }

  template <typename CollisionType>
  void fillCollision(CollisionType const& col)
  {
    registryMixQa.fill(HIST("MixingQA/hSECollisionBins"), colBinningMult.getBin({col.posZ(), col.multNtr()}));
    registryMixQa.fill(HIST("MixingQA/hSECollisionPool"), col.posZ(), col.multNtr());
  }

  /// Compute the charm hadron candidates mass with the daughter masses
  /// assumes the candidate is either a D+ or Λc+
  template <typename Candidate>
  float getCharmHadronMass(const Candidate& cand)
  {
    float invMass = 0.0f;
    if (charmHadPDGCode == o2::constants::physics::Pdg::kLambdaCPlus) {
      if (cand.candidateSelFlag() == 1) {
        invMass = cand.m(std::array{o2::constants::physics::MassProton, o2::constants::physics::MassKPlus, o2::constants::physics::MassPiPlus});
        return invMass;
      }
      invMass = cand.m(std::array{o2::constants::physics::MassPiPlus, o2::constants::physics::MassKPlus, o2::constants::physics::MassProton});
      return invMass;
    }
    // D+ → π K π (PDG: 411)
    if (charmHadPDGCode == o2::constants::physics::Pdg::kDPlus) {
      invMass = cand.m(std::array{o2::constants::physics::MassPiPlus, o2::constants::physics::MassKPlus, o2::constants::physics::MassPiPlus});
      return invMass;
    }
    // Add more channels as needed
    return invMass;
  }

  template <typename Candidate, typename Track>
  float getCharmHadronTrackMass(const Candidate& cand,
                                const Track& trk,
                                double trackMassHyp = o2::constants::physics::MassProton)
  {

    auto pVecProng0 = RecoDecayPtEtaPhi::pVector(cand.prong0Pt(), cand.prong0Eta(), cand.prong0Phi());
    auto pVecProng1 = RecoDecayPtEtaPhi::pVector(cand.prong1Pt(), cand.prong1Eta(), cand.prong1Phi());
    auto pVecProng2 = RecoDecayPtEtaPhi::pVector(cand.prong2Pt(), cand.prong2Eta(), cand.prong2Phi());
    auto pVecTrack = RecoDecayPtEtaPhi::pVector(trk.pt(), trk.eta(), trk.phi());
    const auto pVecCharmTrk = std::array{pVecProng0, pVecProng1, pVecProng2, pVecTrack};

    std::array<double, 4> massCharmTrk{};

    if (charmHadPDGCode == o2::constants::physics::Pdg::kLambdaCPlus) {
      // Λc⁺ → p K π
      if (cand.candidateSelFlag() == 1) {
        massCharmTrk = {
          o2::constants::physics::MassProton,
          o2::constants::physics::MassKPlus,
          o2::constants::physics::MassPiPlus,
          trackMassHyp};
      } else {
        // prong0=π, prong1=K, prong2=p
        massCharmTrk = {
          o2::constants::physics::MassPiPlus,
          o2::constants::physics::MassKPlus,
          o2::constants::physics::MassProton,
          trackMassHyp};
      }
    } else if (charmHadPDGCode == o2::constants::physics::Pdg::kDPlus) {
      // D⁺ → π K π
      massCharmTrk = {
        o2::constants::physics::MassPiPlus,
        o2::constants::physics::MassKPlus,
        o2::constants::physics::MassPiPlus,
        trackMassHyp};
    } else {
      return -1.f;
    }
    return static_cast<float>(RecoDecay::m(pVecCharmTrk, massCharmTrk));
  }

  /// This function processes the same event and takes care of all the histogramming
  template <bool IsMc, typename PartitionType, typename CandType, typename TableTracks, typename Collision>
  void doSameEvent(PartitionType& sliceTrk1, CandType& sliceCharmHad, TableTracks const& parts, Collision const& col)
  {
    fillCollision(col);
    processType = 1; // for same event
    for (auto const& [p1, p2] : combinations(CombinationsFullIndexPolicy(sliceTrk1, sliceCharmHad))) {

      if (p1.trackId() == p2.prong0Id() || p1.trackId() == p2.prong1Id() || p1.trackId() == p2.prong2Id()) {
        continue;
      }

      if (useCPR.value) {
        if (pairCloseRejectionSE.isClosePair(p1, p2, parts, col.magField())) {
          continue;
        }
      }

      if (!pairCleaner.isCleanPair(p1, p2, parts)) {
        continue;
      }

      float kstar = FemtoDreamMath::getkstar(p1, massOne, p2, massTwo);
      if (kstar > highkstarCut) {
        continue;
      }
      float invMass = getCharmHadronMass(p2);

      if (invMass < charmHadMinInvMass || invMass > charmHadMaxInvMass) {
        continue;
      }

      if (p2.pt() < charmHadMinPt || p2.pt() > charmHadMaxPt) {
        continue;
      }

      float deltaInvMassPair = getCharmHadronTrackMass(p2, p1, o2::constants::physics::MassProton) - invMass;

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
      } else {
        pairSign = UnLikeSignPair;
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

  template <bool IsMc, typename CollisionType, typename PartType, typename PartitionType1, typename PartitionType2, typename BinningType>
  void doMixedEvent(CollisionType const& cols, PartType const& parts, PartitionType1& part1, PartitionType2& part2, BinningType policy)
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

      auto sliceTrk1 = part1->sliceByCached(aod::femtodreamparticle::fdCollisionId, collision1.globalIndex(), cache);
      auto sliceCharmHad = part2->sliceByCached(aod::femtodreamparticle::fdCollisionId, collision2.globalIndex(), cache);
      for (const auto& [p1, p2] : combinations(CombinationsFullIndexPolicy(sliceTrk1, sliceCharmHad))) {

        if (useCPR.value) {
          if (pairCloseRejectionME.isClosePair(p1, p2, parts, collision1.magField())) {
            continue;
          }
        }
        if (!pairCleaner.isCleanPair(p1, p2, parts)) {
          continue;
        }

        float kstar = FemtoDreamMath::getkstar(p1, massOne, p2, massTwo);
        if (kstar > highkstarCut) {
          continue;
        }

        float invMass = getCharmHadronMass(p2);

        if (invMass < charmHadMinInvMass || invMass > charmHadMaxInvMass) {
          continue;
        }

        if (p2.pt() < charmHadMinPt || p2.pt() > charmHadMaxPt) {
          continue;
        }

        float deltaInvMassPair = getCharmHadronTrackMass(p2, p1, o2::constants::physics::MassProton) - invMass;

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
          deltaInvMassPair,
          processType,
          charmHadMc,
          originType);

        mixedEventCont.setPair<IsMc, true>(p1, p2, collision1.multNtr(), collision1.multV0M(), use4D, extendedPlots, smearingByOrigin);
      }
    }
  }

  void processSameEvent(FilteredColision const& col,
                        FilteredFDParticles const& parts,
                        FilteredCharmCands const&)
  {
    eventHisto.fillQA(col);
    auto sliceTrk1 = partitionTrk1->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
    auto sliceCharmHad = partitionCharmHadron->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
    if (fillTableWithCharm.value && sliceCharmHad.size() == 0) {
      return;
    }

    int64_t timeStamp = -999;

    for (auto const& part : sliceCharmHad) {
      float invMass = getCharmHadronMass(part);
      registryCharmHadronQa.fill(HIST("CharmHadronQA/hPtVsMass"), part.pt(), invMass);
      timeStamp = part.timeStamp();
      rowFemtoResultCharm(
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
    }

    for (auto const& part : sliceTrk1) {
      allTrackHisto.fillQA<false, true>(part,
                                        static_cast<aod::femtodreamparticle::MomentumType>(confTempFitVarMomentum.value),
                                        col.multNtr(), col.multV0M());

      float chargeTrack = ((part.cut() & CutBitChargePositive) == CutBitChargePositive)
                            ? PositiveCharge
                            : NegativeCharge;
      timeStamp = part.timeStamp();
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
        part.tpcNSigmaPr(),
        part.tofNSigmaPr());
    }
    if (sliceCharmHad.size() > 0 || sliceTrk1.size() > 0) {
      rowFemtoResultColl(
        col.globalIndex(),
        timeStamp,
        col.posZ(),
        col.multNtr());
    } else {
      return;
    }

    doSameEvent<false>(sliceTrk1, sliceCharmHad, parts, col);
  }
  PROCESS_SWITCH(HfTaskCharmHadronsFemtoDream, processSameEvent, "Enable processing same event", false);

  void processMixedEvent(FilteredColisions const& cols,
                         FilteredFDParticles const& parts,
                         FilteredCharmCands const&)
  {
    switch (mixSetting.mixingBinPolicy) {
      case femtodreamcollision::kMult:
        doMixedEvent<false>(cols, parts, partitionTrk1, partitionCharmHadron, colBinningMult);
        break;
      case femtodreamcollision::kMultPercentile:
        doMixedEvent<false>(cols, parts, partitionTrk1, partitionCharmHadron, colBinningMultPercentile);
        break;
      case femtodreamcollision::kMultMultPercentile:
        doMixedEvent<false>(cols, parts, partitionTrk1, partitionCharmHadron, colBinningMultMultPercentile);
        break;
      default:
        LOG(fatal) << "Invalid binning policiy specifed. Breaking...";
    }
  }
  PROCESS_SWITCH(HfTaskCharmHadronsFemtoDream, processMixedEvent, "Enable processing mixed events", false);

  /// process function for to call doSameEvent with Monte Carlo
  /// \param col subscribe to the collision table (Monte Carlo Reconstructed reconstructed)
  /// \param parts subscribe to joined table FemtoDreamParticles and FemtoDreamMCLables to access Monte Carlo truth
  /// \param FemtoDreamMCParticles subscribe to the Monte Carlo truth table
  void processSameEventMc(FilteredMcColision const& col,
                          FilteredFDMcParts const& parts,
                          o2::aod::FDMCParticles const&,
                          o2::aod::FDExtMCParticles const&,
                          FilteredCharmMcCands const&)
  {
    eventHisto.fillQA(col);

    auto sliceMcTrk1 = partitionMcTrk1->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
    auto sliceMcCharmHad = partitionMcCharmHadron->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);

    if (sliceMcCharmHad.size() > 0) {
      for (auto const& part : sliceMcCharmHad) {
        registryCharmHadronQa.fill(HIST("CharmHadronQA/hPtVsMass"), part.pt(), getCharmHadronMass(part));
      }
    }
    if (sliceMcTrk1.size() == 0 && sliceMcCharmHad.size() == 0) {
      return;
    }
    /// Filling QA histograms of the all mc tracks before pairing
    for (auto const& part : sliceMcTrk1) {
      allTrackHisto.fillQA<true, true>(part, static_cast<aod::femtodreamparticle::MomentumType>(confTempFitVarMomentum.value), col.multNtr(), col.multV0M());
    }

    if ((col.bitmaskTrackOne() & bitMask) != bitMask || (col.bitmaskTrackTwo() & bitMask) != bitMask) {
      return;
    }
    doSameEvent<true>(sliceMcTrk1, sliceMcCharmHad, parts, col);
  }
  PROCESS_SWITCH(HfTaskCharmHadronsFemtoDream, processSameEventMc, "Enable processing same event for Monte Carlo", false);

  /// brief process function for to call doMixedEvent with Monte Carlo
  /// @param cols subscribe to the collisions table (Monte Carlo Reconstructed reconstructed)
  /// @param parts subscribe to joined table FemtoDreamParticles and FemtoDreamMCLables to access Monte Carlo truth
  /// @param FemtoDreamMCParticles subscribe to the Monte Carlo truth table
  void processMixedEventMc(FilteredMcColisions const& cols,
                           FilteredFDMcParts const& parts,
                           o2::aod::FDMCParticles const&,
                           o2::aod::FDExtMCParticles const&,
                           FilteredCharmMcCands const&)
  {
    switch (mixSetting.mixingBinPolicy) {
      case femtodreamcollision::kMult:
        doMixedEvent<true>(cols, parts, partitionMcTrk1, partitionMcCharmHadron, colBinningMult);
        break;
      case femtodreamcollision::kMultPercentile:
        doMixedEvent<true>(cols, parts, partitionMcTrk1, partitionMcCharmHadron, colBinningMultPercentile);
        break;
      case femtodreamcollision::kMultMultPercentile:
        doMixedEvent<true>(cols, parts, partitionMcTrk1, partitionMcCharmHadron, colBinningMultMultPercentile);
        break;
      default:
        LOG(fatal) << "Invalid binning policiy specifed. Breaking...";
    }
  }
  PROCESS_SWITCH(HfTaskCharmHadronsFemtoDream, processMixedEventMc, "Enable processing mixed events MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskCharmHadronsFemtoDream>(cfgc)};
}
