
// comments

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

struct femtoDreamPairTaskV0Reso {
  SliceCache cache;
  Preslice<aod::FDParticles> perCol = aod::femtodreamparticle::fdCollisionId;

  Filter EventMultiplicity = aod::femtodreamcollision::multNtr >= EventSel.MultMin && aod::femtodreamcollision::multNtr <= EventSel.MultMax;
  Filter EventMultiplicityPercentile = aod::femtodreamcollision::multV0M >= EventSel.MultPercentileMin && aod::femtodreamcollision::multV0M <= EventSel.MultPercentileMax;

  using FilteredCollisions = soa::Filtered<FDCollisions>;
  using FilteredCollision = FilteredCollisions::iterator;
  // no masked yet

  /// Histogramming for Event
  FemtoDreamEventHisto eventHisto; // are my cases included check & add!

  FemtoDreamContainer<femtoDreamContainer::EventType::same, femtoDreamContainer::Observable::kstar> sameEventCont;
  FemtoDreamContainer<femtoDreamContainer::EventType::mixed, femtoDreamContainer::Observable::kstar> mixedEventCont;
  // FemtoDreamPairCleaner<aod::femtodreamparticle::ParticleType::kV0, aod::femtodreamparticle::ParticleType::kReso> pairCleaner;
  FemtoDreamDetaDphiStar<aod::femtodreamparticle::ParticleType::kV0, aod::femtodreamparticle::ParticleType::kReso> pairCloseRejectionSE;
  FemtoDreamDetaDphiStar<aod::femtodreamparticle::ParticleType::kV0, aod::femtodreamparticle::ParticleType::kReso> pairCloseRejectionME;

  /// General options
  struct : ConfigurableGroup {
    std::string prefix = std::string("Option");
    Configurable<bool> IsMC{"IsMC", false, "Enable additional Histogramms in the case of a MonteCarlo Run"};
    Configurable<bool> Use4D{"Use4D", false, "Enable four dimensional histogramms (to be used only for analysis with high statistics): k* vs multiplicity vs multiplicity percentil vs mT"};
    Configurable<bool> ExtendedPlots{"ExtendedPlots", false, "Enable additional three dimensional histogramms. High memory consumption. Use for debugging"};
    Configurable<float> HighkstarCut{"HighkstarCut", -1., "Set a cut for high k*, above which the pairs are rejected. Set it to -1 to deactivate it"};
    Configurable<bool> CPROn{"CPROn", true, "Close Pair Rejection"};
    Configurable<bool> CPROld{"CPROld", false, "Set to FALSE to use fixed version of CPR (for testing now, will be default soon)"};
    Configurable<bool> CPRPlotPerRadii{"CPRPlotPerRadii", false, "Plot CPR per radii"};
    Configurable<float> CPRdeltaPhiMax{"CPRdeltaPhiMax", 0.01, "Max. Delta Phi for Close Pair Rejection"};
    Configurable<float> CPRdeltaEtaMax{"CPRdeltaEtaMax", 0.01, "Max. Delta Eta for Close Pair Rejection"};
    Configurable<bool> DCACutPtDep{"DCACutPtDep", false, "Use pt dependent dca cut"};
    Configurable<bool> MixEventWithPairs{"MixEventWithPairs", true, "Only use events that contain particle 1 and partile 2 for the event mixing"};
    Configurable<bool> smearingByOrigin{"smearingByOrigin", false, "Obtain the smearing matrix differential in the MC origin of particle 1 and particle 2. High memory consumption. Use with care!"};
    ConfigurableAxis Dummy{"Dummy", {1, 0, 1}, "Dummy axis"};
  } Option;

  /// Event selection
  struct : ConfigurableGroup {
    std::string prefix = std::string("EventSel");
    Configurable<int> MultMin{"MultMin", 0, "Minimum Multiplicity (MultNtr)"};
    Configurable<int> MultMax{"MultMax", 99999, "Maximum Multiplicity (MultNtr)"};
    Configurable<float> MultPercentileMin{"MultPercentileMin", 0, "Minimum Multiplicity Percentile"};
    Configurable<float> MultPercentileMax{"MultPercentileMax", 100, "Maximum Multiplicity Percentile"};
  } EventSel;

  /// Binning configurables
  struct : ConfigurableGroup {
    std::string prefix = std::string("Binning");
    ConfigurableAxis TempFitVarReso{"TempFitVarReso", {300, 0.9, 1}, "binning of the TempFitVar in the pT vs. TempFitVar plot (V0)"};
    ConfigurableAxis TempFitVarV0{"TempFitVarV0", {300, 0.9, 1}, "binning of the TempFitVar in the pT vs. TempFitVar plot (Reso))"};
    ConfigurableAxis TempFitVarV0Child{"TempFitVarV0Child", {300, -0.15, 0.15}, "binning of the TempFitVar in the pT vs. TempFitVar plot (V0 child)"};
    ConfigurableAxis TempFitVarResoChild{"TempFitVarResoChild", {300, -0.15, 0.15}, "binning of the TempFitVar in the pT vs. TempFitVar plot (Reso child)"};
    ConfigurableAxis InvMass{"InvMass", {1500, 0.9, 1.13}, "InvMass binning"};
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
    ConfigurableAxis kstar{"kstar", {1500, 0., 6.}, "binning kstar for the 4Dimensional plot: k* vs multiplicity vs multiplicity percentile vs mT (set <<ConfUse4D>> to true in order to use)"};
    ConfigurableAxis mT{"mT", {VARIABLE_WIDTH, 1.02f, 1.14f, 1.20f, 1.26f, 1.38f, 1.56f, 1.86f, 4.50f}, "mT Binning for the 4Dimensional plot: k* vs multiplicity vs multiplicity percentile vs mT (set <<ConfUse4D>> to true in order to use)"};
    ConfigurableAxis Mult{"mult", {VARIABLE_WIDTH, 0.0f, 4.0f, 8.0f, 12.0f, 16.0f, 20.0f, 24.0f, 28.0f, 32.0f, 36.0f, 40.0f, 44.0f, 48.0f, 52.0f, 56.0f, 60.0f, 64.0f, 68.0f, 72.0f, 76.0f, 80.0f, 84.0f, 88.0f, 92.0f, 96.0f, 100.0f, 200.0f}, "multiplicity Binning for the 4Dimensional plot: k* vs multiplicity vs multiplicity percentile vs mT (set <<ConfUse4D>> to true in order to use)"};
    ConfigurableAxis multPercentile{"multPercentile", {10, 0.0f, 100.0f}, "multiplicity percentile Binning for the 4Dimensional plot: k* vs multiplicity vs multiplicity percentile vs mT (set <<ConfUse4D>> to true in order to use)"};
  } Binning4D;

  // Mixing configurables
  struct : ConfigurableGroup {
    std::string prefix = std::string("Mixing");
    ConfigurableAxis BinMult{"BinMult", {VARIABLE_WIDTH, 0.0f, 4.0f, 8.0f, 12.0f, 16.0f, 20.0f, 24.0f, 28.0f, 32.0f, 36.0f, 40.0f, 44.0f, 48.0f, 52.0f, 56.0f, 60.0f, 64.0f, 68.0f, 72.0f, 76.0f, 80.0f, 84.0f, 88.0f, 92.0f, 96.0f, 100.0f, 200.0f}, "bins - multiplicity"};
    ConfigurableAxis BinMultPercentile{"BinMultPercentile", {VARIABLE_WIDTH, 0.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f, 100.f}, "bins - multiplicity percentile"};
    ConfigurableAxis BinVztx{"BinVztx", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "bins - z-vertex"};
    Configurable<int> Depth{"Depth", 5, "Number of events for mixing"};
    Configurable<int> Policy{"BinPolicy", 0, "Binning policy for mixing - 0: multiplicity, 1: multipliciy percentile, 2: both"};
  } Mixing;

  /// particle 1 (V01), Λ
  struct : ConfigurableGroup {
    std::string prefix = std::string("V01");
    Configurable<int> PDGCode{"PDGCode", 3122, "PDG code of particle 1 (V0)"};
    Configurable<femtodreamparticle::cutContainerType> CutBit{"CutBit", 7518, "Selection bit for particle 1 (V0)"};
    Configurable<femtodreamparticle::cutContainerType> ChildPos_CutBit{"ChildPos_CutBit", 210, "Selection bit for positive child of V01"};
    Configurable<femtodreamparticle::cutContainerType> ChildPos_TPCBit{"ChildPos_TPCBit", 64, "PID TPC bit for positive child of V01"};
    Configurable<femtodreamparticle::cutContainerType> ChildNeg_CutBit{"ChildNeg_CutBit", 209, "Selection bit for negative child of V01"};
    Configurable<femtodreamparticle::cutContainerType> ChildNeg_TPCBit{"ChildNeg_TPCBit", 256, "PID TPC bit for negative child of V01"};

    Configurable<float> InvMassMin{"InvMassMin", 1.08, "Minimum invariant mass of Partricle 1 (particle) (V0)"};
    Configurable<float> InvMassMax{"InvMassMax", 1.15, "Maximum invariant mass of Partricle 1 (particle) (V0)"};
    Configurable<float> InvMassAntiMin{"InvMassAntiMin", 0., "Minimum invariant mass of Partricle 1 (antiparticle) (V0)"}; // should be the same as for Lambda...
    Configurable<float> InvMassAntiMax{"InvMassAntiMax", 999., "Maximum invariant mass of Partricle 1 (antiparticle) (V0)"};

    Configurable<float> PtMin{"PtMin", 0., "Minimum pT of Partricle 1 (V0)"};
    Configurable<float> PtMax{"PtMax", 999., "Maximum pT of Partricle 1 (V0)"};
    Configurable<float> EtaMin{"EtaMin", -10., "Minimum eta of Partricle 1 (V0)"};
    Configurable<float> EtaMax{"EtaMax", 10., "Maximum eta of Partricle 1 (V0)"};
  } V01; // hier evtl noch weiter Configurables einfügen...

  /// particle 2, (Resonance)   (needs implementation phi in cut bit )
  struct : ConfigurableGroup {
    std::string prefix = std::string("Reso2");
    Configurable<int> PDGCode{"PDGCode", 333, "PDG code of particle 2 (V0)"};
    Configurable<femtodreamparticle::cutContainerType> mask_TPC_TPC{"mask_TPC_TPC", 136, "bitmask for TPC and TPC selection for the reconstructd particle"}; // selection masks for the 4 types
    Configurable<femtodreamparticle::cutContainerType> mask_TOF_TOF{"mask_TOF_TOF", 528, "bitmask for TOF and TOF selection for the reconstructd particle"};
    Configurable<femtodreamparticle::cutContainerType> mask_TOF_TPC{"mask_TOF_TPC", 144, "bitmask for TOF and TPC selection for the reconstructd particle"};
    Configurable<femtodreamparticle::cutContainerType> mask_TPC_TOF{"mask_TPC_TOF", 520, "bitmask for TPC and TOF selection for the reconstructd particle"};

    Configurable<float> InvMassMin{"InvMassMin", 1.017, "Minimum invariant mass of Partricle 2 (particle) (V0)"}; // phi values  ofr inv mass
    Configurable<float> InvMassMax{"InvMassMax", 1.027, "Maximum invariant mass of Partricle 2 (particle) (V0)"};
    Configurable<float> PtMin{"PtMin", 0., "Minimum pT of Partricle 2 (V0)"};
    Configurable<float> PtMax{"PtMax", 999., "Maximum pT of Partricle 2 (V0)"};
    Configurable<float> EtaMin{"EtaMin", -10., "Minimum eta of Partricle 2 (V0)"}; // change values
    Configurable<float> EtaMax{"EtaMax", 10., "Maximum eta of Partricle 2 (V0)"};  // change values

    Configurable<femtodreamparticle::cutContainerType> DaughPos_CutBit{"DaughPos_CutBit", 4860458, "Selection bit for positive child of V02"}; // K+
    Configurable<femtodreamparticle::cutContainerType> DaughPos_TPCBit{"DaughPos_TPCBit", 16, "PID TPC bit for positive child of V02"};        // NSigma_TPC = 2.5
    Configurable<femtodreamparticle::cutContainerType> DaughPos_TPCTOFBit{"DaughPos_TOFBit", 8, "PID TOF bit for positive child of V02"};      // NSigma_TOF = 2.5
    Configurable<femtodreamparticle::cutContainerType> DaughNeg_CutBit{"DaughNeg_CutBit", 4860457, "Selection bit for negative child of V02"}; // K-
    Configurable<femtodreamparticle::cutContainerType> DaughNeg_TPCBit{"DaughNeg_TPCBit", 16, "PID TPC bit for negative child of V02"};        // NSigma_TPC = 2.5
    Configurable<femtodreamparticle::cutContainerType> DaughNeg_TPCTOFBit{"DaughNeg_TOFBit", 8, "PID TOF bit for negative child of V02"};      // NSigma_TOF = 2.5
  } Reso2;

  /// Partition for particle 1
  Partition<aod::FDParticles> PartitionV01 = (aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kV0)) &&
                                             ((aod::femtodreamparticle::cut & V01.CutBit) == V01.CutBit) &&
                                             (aod::femtodreamparticle::pt > V01.PtMin) &&
                                             (aod::femtodreamparticle::pt < V01.PtMax) &&
                                             (aod::femtodreamparticle::eta > V01.EtaMin) &&
                                             (aod::femtodreamparticle::eta < V01.EtaMax) &&
                                             (aod::femtodreamparticle::mLambda > V01.InvMassMin) &&
                                             (aod::femtodreamparticle::mLambda < V01.InvMassMax) &&
                                             (aod::femtodreamparticle::mAntiLambda > V01.InvMassAntiMin) &&
                                             (aod::femtodreamparticle::mAntiLambda < V01.InvMassAntiMax);

  /// Partition for particle 2
  Partition<aod::FDParticles> PartitionReso2 = ((ifnode(aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kPhiPosdaughTPC_NegdaughTPC), ncheckbit(aod::femtodreamparticle::pidcut, Reso2.mask_TPC_TPC), false)) ||
                                                (ifnode(aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kPhiPosdaughTOF_NegdaughTOF), ncheckbit(aod::femtodreamparticle::pidcut, Reso2.mask_TOF_TOF), false)) ||
                                                (ifnode(aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kPhiPosdaughTOF_NegdaughTPC), ncheckbit(aod::femtodreamparticle::pidcut, Reso2.mask_TOF_TPC), false)) ||
                                                (ifnode(aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kPhiPosdaughTPC_NegdaughTOF), ncheckbit(aod::femtodreamparticle::pidcut, Reso2.mask_TPC_TOF), false))) &&
                                               (aod::femtodreamparticle::pt < Reso2.PtMax) &&
                                               (aod::femtodreamparticle::eta > Reso2.EtaMin) &&
                                               (aod::femtodreamparticle::eta < Reso2.EtaMax) &&
                                               (aod::femtodreamparticle::mLambda > Reso2.InvMassMin) &&
                                               (aod::femtodreamparticle::mLambda < Reso2.InvMassMax);

  ColumnBinningPolicy<aod::collision::PosZ, aod::femtodreamcollision::MultNtr> colBinningMult{{Mixing.BinVztx, Mixing.BinMult}, true};
  ColumnBinningPolicy<aod::collision::PosZ, aod::femtodreamcollision::MultV0M> colBinningMultPercentile{{Mixing.BinVztx, Mixing.BinMultPercentile}, true};
  ColumnBinningPolicy<aod::collision::PosZ, aod::femtodreamcollision::MultNtr, aod::femtodreamcollision::MultV0M> colBinningMultMultPercentile{{Mixing.BinVztx, Mixing.BinMult, Mixing.BinMultPercentile}, true};

  /// Histogramming for particle 1
  FemtoDreamParticleHisto<aod::femtodreamparticle::ParticleType::kV0, 1> V0HistoPartOne;
  FemtoDreamParticleHisto<aod::femtodreamparticle::ParticleType::kV0Child, 3> posChildHistos;
  FemtoDreamParticleHisto<aod::femtodreamparticle::ParticleType::kV0Child, 4> negChildHistos;

  /// Histogramming for particle 2
  /// prob need to add cases in fillQA, fillDebug in femtoDreamParticleHisto
  FemtoDreamParticleHisto<aod::femtodreamparticle::ParticleType::kReso, 2> ResoHistoPartTwo;
  FemtoDreamParticleHisto<aod::femtodreamparticle::ParticleType::kResoChild, 3> ResoposChildHistos;
  FemtoDreamParticleHisto<aod::femtodreamparticle::ParticleType::kResoChild, 4> ResonegChildHistos;

  /// Histogram output
  HistogramRegistry Registry{"Output", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext&) // InitContext& context
  {

    // setup binnnig policy for mixing
    colBinningMult = {{Mixing.BinVztx, Mixing.BinMult}, true};
    colBinningMultPercentile = {{Mixing.BinVztx, Mixing.BinMultPercentile}, true};
    colBinningMultMultPercentile = {{Mixing.BinVztx, Mixing.BinMult, Mixing.BinMultPercentile}, true};

    eventHisto.init(&Registry, Option.IsMC);
    // change them !!
    V0HistoPartOne.init(&Registry, Binning.multTempFit, Option.Dummy, Binning.pTTrack, Option.Dummy, Option.Dummy, Binning.TempFitVarV0, Option.Dummy, Option.Dummy, Option.Dummy, Option.Dummy, Option.Dummy, Option.Dummy, Option.IsMC, V01.PDGCode);
    posChildHistos.init(&Registry, Binning.multTempFit, Option.Dummy, Binning.pTV0Child, Option.Dummy, Option.Dummy, Binning.TempFitVarV0Child, Option.Dummy, Option.Dummy, Option.Dummy, Option.Dummy, Option.Dummy, Option.Dummy, false, 0);
    negChildHistos.init(&Registry, Binning.multTempFit, Option.Dummy, Binning.pTV0Child, Option.Dummy, Option.Dummy, Binning.TempFitVarV0Child, Option.Dummy, Option.Dummy, Option.Dummy, Option.Dummy, Option.Dummy, Option.Dummy, false, 0);

    ResoHistoPartTwo.init(&Registry, Binning.multTempFit, Option.Dummy, Binning.pTReso, Option.Dummy, Option.Dummy, Binning.TempFitVarReso, Option.Dummy, Option.Dummy, Option.Dummy, Option.Dummy, Binning.InvMass, Option.Dummy, Option.IsMC, Reso2.PDGCode);
    ResoposChildHistos.init(&Registry, Binning.multTempFit, Option.Dummy, Binning.pTResoChild, Option.Dummy, Option.Dummy, Binning.TempFitVarResoChild, Option.Dummy, Option.Dummy, Option.Dummy, Option.Dummy, Option.Dummy, Option.Dummy, false, 0);
    ResonegChildHistos.init(&Registry, Binning.multTempFit, Option.Dummy, Binning.pTResoChild, Option.Dummy, Option.Dummy, Binning.TempFitVarResoChild, Option.Dummy, Option.Dummy, Option.Dummy, Option.Dummy, Option.Dummy, Option.Dummy, false, 0);

    sameEventCont.init(&Registry,
                       Binning.kstar, Binning.pT, Binning.kT, Binning.mT, Mixing.BinMult, Mixing.BinMultPercentile,
                       Binning4D.kstar, Binning4D.mT, Binning4D.Mult, Binning4D.multPercentile,
                       Option.IsMC, Option.Use4D, Option.ExtendedPlots,
                       Option.HighkstarCut,
                       Option.smearingByOrigin, Binning.InvMass);

    sameEventCont.setPDGCodes(V01.PDGCode, Reso2.PDGCode);
    mixedEventCont.init(&Registry,
                        Binning.kstar, Binning.pT, Binning.kT, Binning.mT, Mixing.BinMult, Mixing.BinMultPercentile,
                        Binning4D.kstar, Binning4D.mT, Binning4D.Mult, Binning4D.multPercentile,
                        Option.IsMC, Option.Use4D, Option.ExtendedPlots,
                        Option.HighkstarCut,
                        Option.smearingByOrigin, Binning.InvMass);

    mixedEventCont.setPDGCodes(V01.PDGCode, Reso2.PDGCode);
    // pairCleaner.init(&Registry);
    if (Option.CPROn.value) {
      pairCloseRejectionSE.init(&Registry, &Registry, Option.CPRdeltaPhiMax.value, Option.CPRdeltaEtaMax.value, Option.CPRPlotPerRadii.value, 1, Option.CPROld.value);
      pairCloseRejectionME.init(&Registry, &Registry, Option.CPRdeltaPhiMax.value, Option.CPRdeltaEtaMax.value, Option.CPRPlotPerRadii.value, 2, Option.CPROld.value, 99, true);
    }
  }

  template <typename PartitionType, typename TableTracks, typename Collision>
  void doSameEvent(PartitionType& SliceV01, PartitionType& SliceReso2, TableTracks const& parts, Collision const& col)
  {
    /// Histogramming for same event missing

    for (auto& v0 : SliceV01) {
      const auto& posChild = parts.iteratorAt(v0.index() - 2);
      const auto& negChild = parts.iteratorAt(v0.index() - 1);

      if (((posChild.cut() & V01.ChildPos_CutBit) == V01.ChildPos_CutBit) &&
          ((posChild.pidcut() & V01.ChildPos_TPCBit) == V01.ChildPos_TPCBit) &&
          ((negChild.cut() & V01.ChildNeg_CutBit) == V01.ChildNeg_CutBit) &&
          ((negChild.pidcut() & V01.ChildNeg_TPCBit) == V01.ChildNeg_TPCBit)) {
        V0HistoPartOne.fillQA<false, false>(v0, aod::femtodreamparticle::kPt, col.multNtr(), col.multV0M()); // fillQA<Option.IsMC, false>, here IsDebug == true, false??
        posChildHistos.fillQA<false, false>(posChild, aod::femtodreamparticle::kPt, col.multNtr(), col.multV0M());
        negChildHistos.fillQA<false, false>(negChild, aod::femtodreamparticle::kPt, col.multNtr(), col.multV0M());
      }
    }

    for (auto& reso : SliceReso2) {
      const auto& posresoChild = parts.iteratorAt(reso.index() - 2);
      const auto& negresoChild = parts.iteratorAt(reso.index() - 1);

      if (ncheckbit(posresoChild.cut(), Reso2.DaughPos_CutBit)) {
        ResoposChildHistos.fillQA<false, false>(negresoChild, aod::femtodreamparticle::kPt, col.multNtr(), col.multV0M());
      }

      if (ncheckbit(negresoChild.cut(), Reso2.DaughNeg_CutBit)) {
        ResonegChildHistos.fillQA<false, false>(posresoChild, aod::femtodreamparticle::kPt, col.multNtr(), col.multV0M());
      }

      if (((posresoChild.cut() & Reso2.DaughPos_CutBit) == Reso2.DaughPos_CutBit) &&
          ((negresoChild.cut() & Reso2.DaughNeg_CutBit) == Reso2.DaughNeg_CutBit)) {
        ResoHistoPartTwo.fillQA<false, false>(reso, aod::femtodreamparticle::kPt, col.multNtr(), col.multV0M());
      }
    }

    /// Now build particle combinations vorerst nur not Samespecies!!!
    for (auto const& [p1, p2] : combinations(CombinationsFullIndexPolicy(SliceV01, SliceReso2))) {

      const auto& posChild = parts.iteratorAt(p1.index() - 2);
      const auto& negChild = parts.iteratorAt(p1.index() - 1);

      const auto& posresoChild = parts.iteratorAt(p2.index() - 2);
      const auto& negresoChild = parts.iteratorAt(p2.index() - 1);

      // cuts on V0 children still need to be applied
      if (((posChild.cut() & V01.ChildPos_CutBit) == V01.ChildPos_CutBit) &&
          ((posChild.pidcut() & V01.ChildPos_TPCBit) == V01.ChildPos_TPCBit) &&
          ((negChild.cut() & V01.ChildNeg_CutBit) == V01.ChildNeg_CutBit) &&
          ((negChild.pidcut() & V01.ChildNeg_TPCBit) == V01.ChildNeg_TPCBit) &&

          ((posresoChild.cut() & Reso2.DaughPos_CutBit) == Reso2.DaughPos_CutBit) &&
          ((negresoChild.cut() & Reso2.DaughNeg_CutBit) == Reso2.DaughNeg_CutBit) // TPC & TOF checked in partition...
      ) {
        if (Option.CPROn.value) {
          if (pairCloseRejectionSE.isClosePair(p1, p2, parts, col.magField())) {
            continue;
          }
        }
        sameEventCont.setPair<false>(p1, p2, col.multNtr(), col.multV0M(), Option.Use4D, Option.ExtendedPlots, Option.smearingByOrigin);
      }
    }
  }

  template <bool isMC, typename CollisionType, typename PartType, typename PartitionType, typename BinningType>
  void doMixedEvent(CollisionType const& cols, PartType const& parts, PartitionType& part1, PartitionType& part2, BinningType policy)
  {

    if (Option.MixEventWithPairs.value) {
      for (auto const& [collision1, collision2] : soa::selfCombinations(policy, Mixing.Depth.value, -1, cols, cols)) {
        // make sure that tracks in same events are not mixed
        if (collision1.globalIndex() == collision2.globalIndex()) {
          continue;
        }

        auto SliceV01 = part1->sliceByCached(aod::femtodreamparticle::fdCollisionId, collision1.globalIndex(), cache);
        auto SliceReso2 = part2->sliceByCached(aod::femtodreamparticle::fdCollisionId, collision2.globalIndex(), cache);

        if (SliceV01.size() == 0 || SliceReso2.size() == 0) {
          continue;
        }

        for (auto& [p1, p2] : combinations(CombinationsFullIndexPolicy(SliceV01, SliceReso2))) {

          const auto& posChild = parts.iteratorAt(p1.index() - 2);
          const auto& negChild = parts.iteratorAt(p1.index() - 1);

          const auto& posresoChild = parts.iteratorAt(p2.index() - 2);
          const auto& negresoChild = parts.iteratorAt(p2.index() - 1);

          // why pass if fullfilled??
          if (((posChild.cut() & V01.ChildPos_CutBit) == V01.ChildPos_CutBit) &&
              ((posChild.pidcut() & V01.ChildPos_TPCBit) == V01.ChildPos_TPCBit) &&
              ((negChild.cut() & V01.ChildNeg_CutBit) == V01.ChildNeg_CutBit) &&
              ((negChild.pidcut() & V01.ChildNeg_TPCBit) == V01.ChildNeg_TPCBit) &&

              ((posresoChild.cut() & Reso2.DaughPos_CutBit) == Reso2.DaughPos_CutBit) &&
              ((posresoChild.pidcut() & Reso2.DaughPos_TPCBit) == Reso2.DaughPos_TPCBit) &&
              ((posresoChild.pidcut() & Reso2.DaughPos_TPCTOFBit) == Reso2.DaughPos_TPCTOFBit) && // not really needed already sleceted in partition..
              ((negresoChild.cut() & Reso2.DaughNeg_CutBit) == Reso2.DaughNeg_CutBit) &&
              ((negresoChild.pidcut() & Reso2.DaughNeg_TPCBit) == Reso2.DaughNeg_TPCBit) &&
              ((negresoChild.pidcut() & Reso2.DaughNeg_TPCTOFBit) == Reso2.DaughNeg_TPCTOFBit)) // not really needed already sleceted in partition..
          {
            continue;
          }
          if (Option.CPROn.value) {
            if (pairCloseRejectionME.isClosePair(p1, p2, parts, collision1.magField())) {
              continue;
            }
          }

          mixedEventCont.setPair<false>(p1, p2, collision1.multNtr(), collision1.multV0M(), Option.Use4D, Option.ExtendedPlots, Option.smearingByOrigin);
        }
      }
    }
  }

  void processSameEvent(FilteredCollision& col, FDParticles& parts) // try this.
  {
    // fillCollision<false>(col);
    auto SliceV01 = PartitionV01->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
    auto SliceReso2 = PartitionReso2->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);

    if (SliceV01.size() == 0 && SliceReso2.size() == 0) {
      return;
    }
    eventHisto.fillQA<false>(col);
    doSameEvent(SliceV01, SliceReso2, parts, col);
  }
  PROCESS_SWITCH(femtoDreamPairTaskV0Reso, processSameEvent, "Enable processing same event", true);

  void processMixedEvent(FilteredCollisions& cols, FDParticles const& parts)
  {

    switch (Mixing.Policy.value) {
      case femtodreamcollision::kMult:
        doMixedEvent<false>(cols, parts, PartitionV01, PartitionReso2, colBinningMult);
        break;
      case femtodreamcollision::kMultPercentile:
        doMixedEvent<false>(cols, parts, PartitionV01, PartitionReso2, colBinningMultPercentile);
        break;
      case femtodreamcollision::kMultMultPercentile:
        doMixedEvent<false>(cols, parts, PartitionV01, PartitionReso2, colBinningMultMultPercentile);
        break;
      default:
        LOG(fatal) << "Invalid binning policiy specifed. Breaking...";
    }
  }
  PROCESS_SWITCH(femtoDreamPairTaskV0Reso, processMixedEvent, "Enable processing mixed event", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<femtoDreamPairTaskV0Reso>(cfgc),
  };
  return workflow;
}
