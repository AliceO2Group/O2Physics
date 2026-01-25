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

/// \file femtoDreamPairTaskV0Reso.cxx
/// \brief Tasks that reads the track tables used for the pairing and builds pairs of two tracks
/// \author Christopher Klumm, TU München, christopher.klumm@cern.ch
/// \author Anton Riedel, TU München, anton.riedel@cern.ch
/// \author Georgios Mantzaridis, TU München, georgios.mantzaridis@cern.ch

#include "PWGCF/DataModel/FemtoDerived.h"
#include "PWGCF/FemtoDream/Core/femtoDreamContainer.h"
#include "PWGCF/FemtoDream/Core/femtoDreamDetaDphiStar.h"
#include "PWGCF/FemtoDream/Core/femtoDreamEventHisto.h"
#include "PWGCF/FemtoDream/Core/femtoDreamPairCleaner.h"
#include "PWGCF/FemtoDream/Core/femtoDreamParticleHisto.h"
#include "PWGCF/FemtoDream/Core/femtoDreamUtils.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/Configurable.h"
#include "Framework/Expressions.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/StepTHn.h"
#include "Framework/runDataProcessing.h"

#include "TRandom3.h"

#include <bitset>
#include <cstdint>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::aod;
using namespace o2::soa;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::analysis::femtoDream;

struct FemtoDreamPairTaskV0Reso {
  SliceCache cache;
  Preslice<aod::FDParticles> perCol = aod::femtodreamparticle::fdCollisionId;

  /// Histogramming for Event
  FemtoDreamEventHisto eventHisto;

  FemtoDreamContainer<femtoDreamContainer::EventType::same, femtoDreamContainer::Observable::kstar> sameEventCont;
  FemtoDreamContainer<femtoDreamContainer::EventType::mixed, femtoDreamContainer::Observable::kstar> mixedEventCont;
  // FemtoDreamPairCleaner<aod::femtodreamparticle::ParticleType::kV0, aod::femtodreamparticle::ParticleType::kReso> pairCleaner;
  FemtoDreamDetaDphiStar<aod::femtodreamparticle::ParticleType::kV0, aod::femtodreamparticle::ParticleType::kReso> pairCloseRejectionSE;
  FemtoDreamDetaDphiStar<aod::femtodreamparticle::ParticleType::kV0, aod::femtodreamparticle::ParticleType::kReso> pairCloseRejectionME;

  using FilteredMCCollisions = soa::Filtered<soa::Join<FDCollisions, aod::FDMCCollLabels>>;
  using FilteredMCCollision = FilteredMCCollisions::iterator;

  using FDMCParts = soa::Join<aod::FDParticles, aod::FDMCLabels>;
  using FDMCPart = FDMCParts::iterator;

  /// General options
  struct : ConfigurableGroup {
    std::string prefix = std::string("Option");
    Configurable<bool> isMC{"isMC", false, "Enable additional Histogramms in the case of a MonteCarlo Run"};
    Configurable<bool> use4D{"use4D", false, "Enable four dimensional histogramms (to be used only for analysis with high statistics): k* vs multiplicity vs multiplicity percentil vs mT"};
    Configurable<bool> extendedPlots{"extendedPlots", true, "Enable additional three dimensional histogramms. High memory consumption. Use for debugging"};
    Configurable<float> highkstarCut{"highkstarCut", -1., "Set a cut for high k*, above which the pairs are rejected. Set it to -1 to deactivate it"};
    Configurable<bool> cPROn{"cPROn", false, "Close Pair Rejection"};
    Configurable<bool> cPROld{"cPROld", false, "Set to FALSE to use fixed version of CPR (for testing now, will be default soon)"};
    Configurable<bool> cPRPlotPerRadii{"cPRPlotPerRadii", false, "Plot CPR per radii"};
    Configurable<float> cPRdeltaPhiMax{"cPRdeltaPhiMax", 0.01, "Max. Delta Phi for Close Pair Rejection"};
    Configurable<float> cPRdeltaEtaMax{"cPRdeltaEtaMax", 0.01, "Max. Delta Eta for Close Pair Rejection"};
    Configurable<bool> dCACutPtDep{"dCACutPtDep", false, "Use pt dependent dca cut"};
    Configurable<bool> mixEventWithPairs{"mixEventWithPairs", true, "Only use events that contain particle 1 and partile 2 for the event mixing"};
    Configurable<bool> smearingByOrigin{"smearingByOrigin", false, "Obtain the smearing matrix differential in the MC origin of particle 1 and particle 2. High memory consumption. Use with care!"};
    ConfigurableAxis dummy{"dummy", {1, 0, 1}, "dummy axis"};
  } Option;

  /// Event selection
  struct : ConfigurableGroup {
    std::string prefix = std::string("EventSel");
    Configurable<int> multMin{"multMin", 0, "Minimum Multiplicity (MultNtr)"};
    Configurable<int> multMax{"multMax", 99999, "Maximum Multiplicity (MultNtr)"};
    Configurable<float> multPercentileMin{"multPercentileMin", 0, "Minimum Multiplicity Percentile"};
    Configurable<float> multPercentileMax{"multPercentileMax", 100, "Maximum Multiplicity Percentile"};
  } EventSel;

  /// Binning configurables
  struct : ConfigurableGroup {
    std::string prefix = std::string("Binning");
    ConfigurableAxis tempFitVarReso{"tempFitVarReso", {300, 0.9, 1}, "binning of the TempFitVar in the pT vs. TempFitVar plot (V0)"};
    ConfigurableAxis tempFitVarV0{"tempFitVarV0", {300, 0.9, 1}, "binning of the TempFitVar in the pT vs. TempFitVar plot (Reso))"};
    ConfigurableAxis tempFitVarV0Child{"tempFitVarV0Child", {300, -0.15, 0.15}, "binning of the TempFitVar in the pT vs. TempFitVar plot (V0 child)"};
    ConfigurableAxis tempFitVarResoChild{"tempFitVarResoChild", {300, -0.15, 0.15}, "binning of the TempFitVar in the pT vs. TempFitVar plot (Reso child)"};
    ConfigurableAxis invMass{"invMass", {1500, 0.9, 1.13}, "invMass binning"};
    ConfigurableAxis pTTrack{"pTTrack", {20, 0.5, 4.05}, "pT binning of the pT vs. TempFitVar plot (Track)"};
    ConfigurableAxis pTV0{"pTV0", {20, 0.5, 4.05}, "pT binning of the pT vs. TempFitVar plot (V0)"};
    ConfigurableAxis pTReso{"pTReso", {20, 0.5, 4.05}, "pT binning of the pT vs. TempFitVar plot (Reso)"};
    ConfigurableAxis pTV0Child{"pTV0Child", {20, 0.5, 4.05}, "pT binning of the pT vs. TempFitVar plot (V0 Child)"};
    ConfigurableAxis pTResoChild{"pTResoChild", {20, 0.5, 4.05}, "pT binning of the pT vs. TempFitVar plot (Reso Child)"};
    ConfigurableAxis pT{"pT", {20, 0.5, 4.05}, "pT binning"};
    ConfigurableAxis kstar{"kstar", {1500, 0., 6.}, "binning kstar"};
    ConfigurableAxis kT{"kT", {150, 0., 9.}, "binning kT"};
    ConfigurableAxis mT{"mT", {225, 0., 7.5}, "binning mT"};
    ConfigurableAxis multTempFit{"multTempFit", {1, 0, 1}, "multiplicity for the TempFitVar plot"};
  } Binning;

  struct : ConfigurableGroup {
    std::string prefix = std::string("Binning4D");
    ConfigurableAxis kstar{"kstar", {1500, 0., 6.}, "binning kstar for the 4Dimensional plot: k* vs multiplicity vs multiplicity percentile vs mT (set <<Confuse4D>> to true in order to use)"};
    ConfigurableAxis mT{"mT", {VARIABLE_WIDTH, 1.02f, 1.14f, 1.20f, 1.26f, 1.38f, 1.56f, 1.86f, 4.50f}, "mT Binning for the 4Dimensional plot: k* vs multiplicity vs multiplicity percentile vs mT (set <<Confuse4D>> to true in order to use)"};
    ConfigurableAxis mult{"mult", {VARIABLE_WIDTH, 0.0f, 4.0f, 8.0f, 12.0f, 16.0f, 20.0f, 24.0f, 28.0f, 32.0f, 36.0f, 40.0f, 44.0f, 48.0f, 52.0f, 56.0f, 60.0f, 64.0f, 68.0f, 72.0f, 76.0f, 80.0f, 84.0f, 88.0f, 92.0f, 96.0f, 100.0f, 200.0f}, "multiplicity Binning for the 4Dimensional plot: k* vs multiplicity vs multiplicity percentile vs mT (set <<Confuse4D>> to true in order to use)"};
    ConfigurableAxis multPercentile{"multPercentile", {10, 0.0f, 100.0f}, "multiplicity percentile Binning for the 4Dimensional plot: k* vs multiplicity vs multiplicity percentile vs mT (set <<Confuse4D>> to true in order to use)"};
  } Binning4D;

  // Mixing configurables
  struct : ConfigurableGroup {
    std::string prefix = std::string("Mixing");
    ConfigurableAxis binMult{"binMult", {VARIABLE_WIDTH, 0.0f, 4.0f, 8.0f, 12.0f, 16.0f, 20.0f, 24.0f, 28.0f, 32.0f, 36.0f, 40.0f, 44.0f, 48.0f, 52.0f, 56.0f, 60.0f, 64.0f, 68.0f, 72.0f, 76.0f, 80.0f, 84.0f, 88.0f, 92.0f, 96.0f, 100.0f, 200.0f}, "bins - multiplicity"};
    ConfigurableAxis binMultPercentile{"binMultPercentile", {VARIABLE_WIDTH, 0.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f, 100.f}, "bins - multiplicity percentile"};
    ConfigurableAxis binVztx{"binVztx", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "bins - z-vertex"};
    Configurable<int> depth{"depth", 5, "Number of events for mixing"};
    Configurable<int> policy{"policy", 0, "Binning policy for mixing - 0: multiplicity, 1: multipliciy percentile, 2: both"};
  } Mixing;

  /// particle 1 (V01), Λ
  struct : ConfigurableGroup {
    std::string prefix = std::string("V01");
    Configurable<int> pdgCode{"pdgCode", 3122, "PDG code of particle 1 (V0)"};
    Configurable<femtodreamparticle::cutContainerType> cutBit{"cutBit", 7518, "Selection bit for particle 1 (V0)"};
    Configurable<femtodreamparticle::cutContainerType> childPosCutBit{"childPosCutBit", 210, "Selection bit for positive child of V01"};
    Configurable<femtodreamparticle::cutContainerType> childPosTPCBit{"childPosTPCBit", 64, "PID TPC bit for positive child of V01"};
    Configurable<femtodreamparticle::cutContainerType> childNegCutBit{"childNegCutBit", 209, "Selection bit for negative child of V01"};
    Configurable<femtodreamparticle::cutContainerType> childNegTPCBit{"childNegTPCBit", 256, "PID TPC bit for negative child of V01"};

    Configurable<float> invMassMin{"invMassMin", 1.08, "Minimum invariant mass of Partricle 1 (particle) (V0)"};
    Configurable<float> invMassMax{"invMassMax", 1.3, "Maximum invariant mass of Partricle 1 (particle) (V0)"};
    Configurable<float> invMassAntiMin{"invMassAntiMin", 0., "Minimum invariant mass of Partricle 1 (antiparticle) (V0)"}; // should be the same as for Lambda...
    Configurable<float> invMassAntiMax{"invMassAntiMax", 999., "Maximum invariant mass of Partricle 1 (antiparticle) (V0)"};

    Configurable<float> ptMin{"ptMin", 0., "Minimum pT of Partricle 1 (V0)"};
    Configurable<float> ptMax{"ptMax", 999., "Maximum pT of Partricle 1 (V0)"};
    Configurable<float> etaMin{"etaMin", -10., "Minimum eta of Partricle 1 (V0)"};
    Configurable<float> etaMax{"etaMax", 10., "Maximum eta of Partricle 1 (V0)"};
  } V01;

  /// particle 2, (Resonance)   (needs implementation phi in cut bit )
  struct : ConfigurableGroup {
    std::string prefix = std::string("Reso2");
    Configurable<int> pdgCode{"pdgCode", 333, "PDG code of particle 2 (V0)"};

    Configurable<float> invMassMin{"invMassMin", 0.86, "Minimum invariant mass of Partricle 2 (particle) (V0)"}; // phi values  for inv mass
    Configurable<float> invMassMax{"invMassMax", 1.3, "Maximum invariant mass of Partricle 2 (particle) (V0)"};
    Configurable<float> ptMin{"ptMin", 0., "Minimum pT of Partricle 2 (V0)"};
    Configurable<float> ptMax{"ptMax", 999., "Maximum pT of Partricle 2 (V0)"};
    Configurable<float> etaMin{"etaMin", -10., "Minimum eta of Partricle 2 (V0)"}; // change values
    Configurable<float> etaMax{"etaMax", 10., "Maximum eta of Partricle 2 (V0)"};  // change values

    Configurable<femtodreamparticle::cutContainerType> daughPosCutBit{"daughPosCutBit", 2401446, "Selection bit for positive child of V02"}; // K+
    Configurable<femtodreamparticle::cutContainerType> daughPosTPCBit{"daughPosTPCBit", 4096, "PID TPC bit for positive child of V02"};
    Configurable<femtodreamparticle::cutContainerType> daughPosTPCTOFBit{"daughPosTPCTOFBit", 2048, "PID TOF bit for positive child of V02"};
    Configurable<femtodreamparticle::cutContainerType> daughNegCutBit{"daughNegCutBit", 2401445, "Selection bit for negative child of V02"}; // K-
    Configurable<femtodreamparticle::cutContainerType> daughNegMergedTPCBit{"daughNegMergedTPCBit", 16386, "PID TPC bit for negative child of V02"};
    Configurable<femtodreamparticle::cutContainerType> daughNegMergedTPCTOFBit{"daughNegMergedTPCTOFBit", 8194, "PID TOF bit for negative child of V02"};

    Configurable<float> dcaXYPar0{"dcaXYPar0", 0.004, "first parameter for pt dependent dcaXY cut"};
    Configurable<float> dcaXYPar1{"dcaXYPar1", 0.013, "second parameter for pt dependent dcaXY cut"};
  } Reso2;

  /// Partition for particle Lambda
  Partition<aod::FDParticles> partitionV01 = (aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kV0)) &&
                                             ((aod::femtodreamparticle::cut & V01.cutBit) == V01.cutBit) &&
                                             (aod::femtodreamparticle::pt > V01.ptMin) &&
                                             (aod::femtodreamparticle::pt < V01.ptMax) &&
                                             (aod::femtodreamparticle::eta > V01.etaMin) &&
                                             (aod::femtodreamparticle::eta < V01.etaMax) &&
                                             (aod::femtodreamparticle::mLambda > V01.invMassMin) &&
                                             (aod::femtodreamparticle::mLambda < V01.invMassMax) &&
                                             (aod::femtodreamparticle::mAntiLambda > V01.invMassAntiMin) &&
                                             (aod::femtodreamparticle::mAntiLambda < V01.invMassAntiMax);

  /// Partition for particle Phi
  Partition<aod::FDParticles> partitionReso2 = (ifnode(aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kResoPosdaughTPC_NegdaughTPC), ncheckbit(aod::femtodreamparticle::pidcut, Reso2.daughPosTPCBit) && ncheckbit(aod::femtodreamparticle::cut, Reso2.daughNegMergedTPCBit), false) ||
                                                ifnode(aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kResoPosdaughTOF_NegdaughTOF), ncheckbit(aod::femtodreamparticle::pidcut, Reso2.daughPosTPCTOFBit) && ncheckbit(aod::femtodreamparticle::cut, Reso2.daughNegMergedTPCTOFBit), false) ||
                                                ifnode(aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kResoPosdaughTOF_NegdaughTPC), ncheckbit(aod::femtodreamparticle::pidcut, Reso2.daughPosTPCTOFBit) && ncheckbit(aod::femtodreamparticle::cut, Reso2.daughNegMergedTPCBit), false) ||
                                                ifnode(aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kResoPosdaughTPC_NegdaughTOF), ncheckbit(aod::femtodreamparticle::pidcut, Reso2.daughPosTPCBit) && ncheckbit(aod::femtodreamparticle::cut, Reso2.daughNegMergedTPCTOFBit), false)) &&
                                               (aod::femtodreamparticle::pt < Reso2.ptMax) &&
                                               (aod::femtodreamparticle::eta > Reso2.etaMin) &&
                                               (aod::femtodreamparticle::eta < Reso2.etaMax) &&
                                               (aod::femtodreamparticle::mLambda > Reso2.invMassMin) &&
                                               (aod::femtodreamparticle::mLambda < Reso2.invMassMax);

  /// Partitions for K0Short and KStar

  /// Partition for particle K0Short
  Partition<aod::FDParticles> partitionK0Short1 = (aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kV0K0Short)) &&
                                                  ((aod::femtodreamparticle::cut & V01.cutBit) == V01.cutBit) &&
                                                  (aod::femtodreamparticle::pt > V01.ptMin) &&
                                                  (aod::femtodreamparticle::pt < V01.ptMax) &&
                                                  (aod::femtodreamparticle::eta > V01.etaMin) &&
                                                  (aod::femtodreamparticle::eta < V01.etaMax) &&
                                                  (aod::femtodreamparticle::mLambda > V01.invMassMin) &&
                                                  (aod::femtodreamparticle::mLambda < V01.invMassMax) &&
                                                  (aod::femtodreamparticle::mAntiLambda > V01.invMassAntiMin) &&
                                                  (aod::femtodreamparticle::mAntiLambda < V01.invMassAntiMax);

  /// Partition for particle KStar
  Partition<aod::FDParticles> partitionKStar2 = (ifnode(aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kResoKStarPosdaughTPC_NegdaughTPC), ncheckbit(aod::femtodreamparticle::pidcut, Reso2.daughPosTPCBit) && ncheckbit(aod::femtodreamparticle::cut, Reso2.daughNegMergedTPCBit), false) ||
                                                 ifnode(aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kResoKStarPosdaughTOF_NegdaughTOF), ncheckbit(aod::femtodreamparticle::pidcut, Reso2.daughPosTPCTOFBit) && ncheckbit(aod::femtodreamparticle::cut, Reso2.daughNegMergedTPCTOFBit), false) ||
                                                 ifnode(aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kResoKStarPosdaughTOF_NegdaughTPC), ncheckbit(aod::femtodreamparticle::pidcut, Reso2.daughPosTPCTOFBit) && ncheckbit(aod::femtodreamparticle::cut, Reso2.daughNegMergedTPCBit), false) ||
                                                 ifnode(aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kResoKStarPosdaughTPC_NegdaughTOF), ncheckbit(aod::femtodreamparticle::pidcut, Reso2.daughPosTPCBit) && ncheckbit(aod::femtodreamparticle::cut, Reso2.daughNegMergedTPCTOFBit), false)) &&
                                                (aod::femtodreamparticle::pt < Reso2.ptMax) &&
                                                (aod::femtodreamparticle::eta > Reso2.etaMin) &&
                                                (aod::femtodreamparticle::eta < Reso2.etaMax) &&
                                                (aod::femtodreamparticle::mLambda > Reso2.invMassMin) &&
                                                (aod::femtodreamparticle::mLambda < Reso2.invMassMax);

  ColumnBinningPolicy<aod::collision::PosZ, aod::femtodreamcollision::MultNtr> colBinningMult{{Mixing.binVztx, Mixing.binMult}, true};
  ColumnBinningPolicy<aod::collision::PosZ, aod::femtodreamcollision::MultV0M> colBinningMultPercentile{{Mixing.binVztx, Mixing.binMultPercentile}, true};
  ColumnBinningPolicy<aod::collision::PosZ, aod::femtodreamcollision::MultNtr, aod::femtodreamcollision::MultV0M> colBinningMultMultPercentile{{Mixing.binVztx, Mixing.binMult, Mixing.binMultPercentile}, true};

  Filter eventMultiplicity = aod::femtodreamcollision::multNtr >= EventSel.multMin && aod::femtodreamcollision::multNtr <= EventSel.multMax;
  Filter eventMultiplicityPercentile = aod::femtodreamcollision::multV0M >= EventSel.multPercentileMin && aod::femtodreamcollision::multV0M <= EventSel.multPercentileMax;

  using FilteredCollisions = soa::Filtered<FDCollisions>;
  using FilteredCollision = FilteredCollisions::iterator;

  /// Histogramming for particle 1
  FemtoDreamParticleHisto<aod::femtodreamparticle::ParticleType::kV0, 1> v0HistoPartOne;
  FemtoDreamParticleHisto<aod::femtodreamparticle::ParticleType::kV0Child, 3> posChildHistos;
  FemtoDreamParticleHisto<aod::femtodreamparticle::ParticleType::kV0Child, 4> negChildHistos;

  /// Histogramming for particle 2
  FemtoDreamParticleHisto<aod::femtodreamparticle::ParticleType::kReso, 2> resoHistoPartTwo;
  FemtoDreamParticleHisto<aod::femtodreamparticle::ParticleType::kResoChild, 3> resoposChildHistos;
  FemtoDreamParticleHisto<aod::femtodreamparticle::ParticleType::kResoChild, 4> resonegChildHistos;

  /// Histogram output
  HistogramRegistry registry{"Output", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry resoRegistry{"ResodcaXY", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext&) // InitContext& context
  {

    // setup binnnig policy for mixing
    colBinningMult = {{Mixing.binVztx, Mixing.binMult}, true};
    colBinningMultPercentile = {{Mixing.binVztx, Mixing.binMultPercentile}, true};
    colBinningMultMultPercentile = {{Mixing.binVztx, Mixing.binMult, Mixing.binMultPercentile}, true};

    eventHisto.init(&registry, Option.isMC);

    v0HistoPartOne.init(&registry, Binning.multTempFit, Option.dummy, Binning.pTTrack, Option.dummy, Option.dummy, Binning.tempFitVarV0, Option.dummy, Option.dummy, Option.dummy, Option.dummy, Binning.invMass, Option.dummy, Option.isMC, V01.pdgCode);
    posChildHistos.init(&registry, Binning.multTempFit, Option.dummy, Binning.pTV0Child, Option.dummy, Option.dummy, Binning.tempFitVarV0Child, Option.dummy, Option.dummy, Option.dummy, Option.dummy, Option.dummy, Option.dummy, false, 0);
    negChildHistos.init(&registry, Binning.multTempFit, Option.dummy, Binning.pTV0Child, Option.dummy, Option.dummy, Binning.tempFitVarV0Child, Option.dummy, Option.dummy, Option.dummy, Option.dummy, Option.dummy, Option.dummy, false, 0);

    resoHistoPartTwo.init(&registry, Binning.multTempFit, Option.dummy, Binning.pTReso, Option.dummy, Option.dummy, Binning.tempFitVarReso, Option.dummy, Option.dummy, Option.dummy, Option.dummy, Binning.invMass, Option.dummy, Option.isMC, Reso2.pdgCode);
    resoposChildHistos.init(&registry, Binning.multTempFit, Option.dummy, Binning.pTResoChild, Option.dummy, Option.dummy, Binning.tempFitVarResoChild, Option.dummy, Option.dummy, Option.dummy, Option.dummy, Option.dummy, Option.dummy, false, 0);
    resonegChildHistos.init(&registry, Binning.multTempFit, Option.dummy, Binning.pTResoChild, Option.dummy, Option.dummy, Binning.tempFitVarResoChild, Option.dummy, Option.dummy, Option.dummy, Option.dummy, Option.dummy, Option.dummy, false, 0);

    sameEventCont.init(&registry,
                       Binning.kstar, Binning.pT, Binning.kT, Binning.mT, Mixing.binMult, Mixing.binMultPercentile,
                       Binning4D.kstar, Binning4D.mT, Binning4D.mult, Binning4D.multPercentile,
                       Option.isMC, Option.use4D, Option.extendedPlots,
                       Option.highkstarCut,
                       Option.smearingByOrigin, Binning.invMass);

    sameEventCont.setPDGCodes(V01.pdgCode, Reso2.pdgCode);
    mixedEventCont.init(&registry,
                        Binning.kstar, Binning.pT, Binning.kT, Binning.mT, Mixing.binMult, Mixing.binMultPercentile,
                        Binning4D.kstar, Binning4D.mT, Binning4D.mult, Binning4D.multPercentile,
                        Option.isMC, Option.use4D, Option.extendedPlots,
                        Option.highkstarCut,
                        Option.smearingByOrigin, Binning.invMass);

    mixedEventCont.setPDGCodes(V01.pdgCode, Reso2.pdgCode);

    resoRegistry.add("Before/DCAxyPt", "DCAxyvsPt", HistType::kTH2F, {{100, -0.8, 0.8}, {100, 0.0, 4}});
    resoRegistry.add("After/DCAxyPt", "DCAxyvsPt", HistType::kTH2F, {{100, -0.8, 0.8}, {100, 0.0, 4}});

    // pairCleaner.init(&registry);
    if (Option.cPROn.value) {
      pairCloseRejectionSE.init(&registry, &registry, Option.cPRdeltaPhiMax.value, Option.cPRdeltaEtaMax.value, Option.cPRPlotPerRadii.value, 1, Option.cPROld.value);
      pairCloseRejectionME.init(&registry, &registry, Option.cPRdeltaPhiMax.value, Option.cPRdeltaEtaMax.value, Option.cPRPlotPerRadii.value, 2, Option.cPROld.value, 99, true);
    }
  }

  template <bool isMC, typename sliceType1, typename sliceType2, typename TableTracks, typename Collision>
  void doSameEvent(sliceType1& sliceV01, sliceType2& sliceReso2, TableTracks const& parts, Collision const& col)
  {
    /// Histogramming for same event missing

    for (const auto& v0 : sliceV01) {
      const auto& posChild = parts.iteratorAt(v0.index() - 2);
      const auto& negChild = parts.iteratorAt(v0.index() - 1);

      if (((posChild.cut() & V01.childPosCutBit) == V01.childPosCutBit) &&
          ((posChild.pidcut() & V01.childPosTPCBit) == V01.childPosTPCBit) &&
          ((negChild.cut() & V01.childNegCutBit) == V01.childNegCutBit) &&
          ((negChild.pidcut() & V01.childNegTPCBit) == V01.childNegTPCBit)) {
        v0HistoPartOne.fillQA<isMC, false>(v0, aod::femtodreamparticle::kPt, col.multNtr(), col.multV0M());
        posChildHistos.fillQA<false, false>(posChild, aod::femtodreamparticle::kPt, col.multNtr(), col.multV0M());
        negChildHistos.fillQA<false, false>(negChild, aod::femtodreamparticle::kPt, col.multNtr(), col.multV0M());
      }
    }

    for (const auto& reso : sliceReso2) {
      const auto& posresoChild = parts.iteratorAt(reso.index() - 2);
      const auto& negresoChild = parts.iteratorAt(reso.index() - 1);

      resoRegistry.fill(HIST("Before/DCAxyPt"), posresoChild.tempFitVar(), posresoChild.pt());
      resoRegistry.fill(HIST("Before/DCAxyPt"), negresoChild.tempFitVar(), negresoChild.pt());

      if ((std::abs(posresoChild.tempFitVar()) > Reso2.dcaXYPar0 + Reso2.dcaXYPar1 * std::pow(posresoChild.pt(), -1)) || (std::abs(negresoChild.tempFitVar()) > Reso2.dcaXYPar0 + Reso2.dcaXYPar1 * std::pow(negresoChild.pt(), -1))) {
        continue;
      }

      resoRegistry.fill(HIST("After/DCAxyPt"), posresoChild.tempFitVar(), posresoChild.pt());
      resoRegistry.fill(HIST("After/DCAxyPt"), negresoChild.tempFitVar(), negresoChild.pt());

      if (ncheckbit(posresoChild.cut(), Reso2.daughPosCutBit)) {
        resoposChildHistos.fillQA<false, false>(posresoChild, aod::femtodreamparticle::kPt, col.multNtr(), col.multV0M());
      }

      if (ncheckbit(negresoChild.cut(), Reso2.daughNegCutBit)) {
        resonegChildHistos.fillQA<false, false>(negresoChild, aod::femtodreamparticle::kPt, col.multNtr(), col.multV0M());
      }

      if (((posresoChild.cut() & Reso2.daughPosCutBit) == Reso2.daughPosCutBit) &&
          ((negresoChild.cut() & Reso2.daughNegCutBit) == Reso2.daughNegCutBit)) {
        resoHistoPartTwo.fillQA<isMC, false>(reso, aod::femtodreamparticle::kPt, col.multNtr(), col.multV0M()); // improve
      }
    }

    /// build particle combinations
    for (auto const& [p1, p2] : combinations(CombinationsFullIndexPolicy(sliceV01, sliceReso2))) {

      const auto& posChild = parts.iteratorAt(p1.index() - 2);
      const auto& negChild = parts.iteratorAt(p1.index() - 1);

      const auto& posresoChild = parts.iteratorAt(p2.index() - 2);
      const auto& negresoChild = parts.iteratorAt(p2.index() - 1);

      if (((posChild.cut() & V01.childPosCutBit) == V01.childPosCutBit) &&
          ((posChild.pidcut() & V01.childPosTPCBit) == V01.childPosTPCBit) &&
          ((negChild.cut() & V01.childNegCutBit) == V01.childNegCutBit) &&
          ((negChild.pidcut() & V01.childNegTPCBit) == V01.childNegTPCBit) &&

          ((posresoChild.cut() & Reso2.daughPosCutBit) == Reso2.daughPosCutBit) &&
          ((negresoChild.cut() & Reso2.daughNegCutBit) == Reso2.daughNegCutBit)) {

        if (Option.cPROn.value) {
          if (pairCloseRejectionSE.isClosePair(p1, p2, parts, col.magField())) {
            continue;
          }
        }
        sameEventCont.setPair<isMC>(p1, p2, col.multNtr(), col.multV0M(), Option.use4D, Option.extendedPlots, Option.smearingByOrigin);
      }
    }
  }

  template <bool isMC, typename CollisionType, typename PartType, typename PartitionType, typename BinningType>
  void doMixedEvent(CollisionType const& cols, PartType const& parts, PartitionType& part1, PartitionType& part2, BinningType policy)
  {

    if (Option.mixEventWithPairs.value) {
      for (auto const& [collision1, collision2] : soa::selfCombinations(policy, Mixing.depth.value, -1, cols, cols)) {
        // make sure that tracks in same events are not mixed
        if (collision1.globalIndex() == collision2.globalIndex()) {
          continue;
        }

        auto sliceV01 = part1.sliceByCached(aod::femtodreamparticle::fdCollisionId, collision1.globalIndex(), cache); // maybe use .
        auto sliceReso2 = part2.sliceByCached(aod::femtodreamparticle::fdCollisionId, collision2.globalIndex(), cache);

        if (sliceV01.size() == 0 || sliceReso2.size() == 0) {
          continue;
        }

        for (const auto& [p1, p2] : combinations(CombinationsFullIndexPolicy(sliceV01, sliceReso2))) {

          const auto& posChild = parts.iteratorAt(p1.globalIndex() - 2);
          const auto& negChild = parts.iteratorAt(p1.globalIndex() - 1);

          const auto& posresoChild = parts.iteratorAt(p2.globalIndex() - 2);
          const auto& negresoChild = parts.iteratorAt(p2.globalIndex() - 1);

          if ((std::abs(posresoChild.tempFitVar()) > Reso2.dcaXYPar0 + Reso2.dcaXYPar1 * std::pow(posresoChild.pt(), -1)) || (std::abs(negresoChild.tempFitVar()) > Reso2.dcaXYPar0 + Reso2.dcaXYPar1 * std::pow(negresoChild.pt(), -1))) {
            continue;
          }

          // why pass if fullfilled??
          if ((((posChild.cut() & V01.childPosCutBit) == V01.childPosCutBit) &&
               ((posChild.pidcut() & V01.childPosTPCBit) == V01.childPosTPCBit) &&
               ((negChild.cut() & V01.childNegCutBit) == V01.childNegCutBit) &&
               ((negChild.pidcut() & V01.childNegTPCBit) == V01.childNegTPCBit) &&

               ((posresoChild.cut() & Reso2.daughPosCutBit) == Reso2.daughPosCutBit) &&
               ((negresoChild.cut() & Reso2.daughNegCutBit) == Reso2.daughNegCutBit))) {

            if (Option.cPROn.value) {
              if (pairCloseRejectionME.isClosePair(p1, p2, parts, collision1.magField())) {
                continue;
              }
            }
            mixedEventCont.setPair<isMC>(p1, p2, collision1.multNtr(), collision1.multV0M(), Option.use4D, Option.extendedPlots, Option.smearingByOrigin);
          }
        }
      }
    }
  }

  void processSameEvent(const FilteredCollision& col, const FDParticles& parts)
  {
    // fillCollision<false>(col);
    auto sliceV01 = partitionV01.sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
    auto sliceReso2 = partitionReso2.sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);

    // if (sliceV01.size() == 0 && sliceReso2.size() == 0) {
    //   return;
    // }
    eventHisto.fillQA<false>(col);
    doSameEvent<false>(sliceV01, sliceReso2, parts, col);
  }
  PROCESS_SWITCH(FemtoDreamPairTaskV0Reso, processSameEvent, "Enable processing same event", true);

  void processMixedEvent(const FilteredCollisions& cols, const FDParticles& parts)
  {
    switch (Mixing.policy.value) {
      case femtodreamcollision::kMult:
        doMixedEvent<false>(cols, parts, partitionV01, partitionReso2, colBinningMult);
        break;
      case femtodreamcollision::kMultPercentile:
        doMixedEvent<false>(cols, parts, partitionV01, partitionReso2, colBinningMultPercentile);
        break;
      case femtodreamcollision::kMultMultPercentile:
        doMixedEvent<false>(cols, parts, partitionV01, partitionReso2, colBinningMultMultPercentile);
        break;
      default:
        LOG(fatal) << "Invalid binning policiy specifed. Breaking...";
    }
  }
  PROCESS_SWITCH(FemtoDreamPairTaskV0Reso, processMixedEvent, "Enable processing mixed event", true);

  //////////////////////////////////////
  /// procees functions for K0Short-KStar

  void processSameEventK0ShortKStar(const FilteredCollision& col, const FDParticles& parts)
  {
    auto sliceV01 = partitionK0Short1.sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
    auto sliceReso2 = partitionKStar2.sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);

    eventHisto.fillQA<false>(col);
    doSameEvent<false>(sliceV01, sliceReso2, parts, col);
  }
  PROCESS_SWITCH(FemtoDreamPairTaskV0Reso, processSameEventK0ShortKStar, "Enable processing same event K0Short-KStar", false);

  void processMixedEventK0ShortKStar(const FilteredCollisions& cols, const FDParticles& parts)
  {
    switch (Mixing.policy.value) {
      case femtodreamcollision::kMult:
        doMixedEvent<false>(cols, parts, partitionK0Short1, partitionKStar2, colBinningMult);
        break;
      case femtodreamcollision::kMultPercentile:
        doMixedEvent<false>(cols, parts, partitionK0Short1, partitionKStar2, colBinningMultPercentile);
        break;
      case femtodreamcollision::kMultMultPercentile:
        doMixedEvent<false>(cols, parts, partitionK0Short1, partitionKStar2, colBinningMultMultPercentile);
        break;
      default:
        LOG(fatal) << "Invalid binning policiy specifed. Breaking...";
    }
  }
  PROCESS_SWITCH(FemtoDreamPairTaskV0Reso, processMixedEventK0ShortKStar, "Enable processing mixed event K0Short-KStar", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<FemtoDreamPairTaskV0Reso>(cfgc),
  };
  return workflow;
}
