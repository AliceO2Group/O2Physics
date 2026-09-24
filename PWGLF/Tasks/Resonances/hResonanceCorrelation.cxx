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
//
/// \file hResonanceCorrelation.cxx
/// \brief This task serves to do hadron-resonance correlation studies.
///  The yield will be calculated using the two-particle correlation method.
///  Trigger particle : Hadrons
///  Associated Particles : Phi, K*0
///  this task requires the hResonanceCorrelationFilter to have been run before.
///
/// \author Hirak Kumar Koley (hirak.koley@cern.ch)

#include "PWGLF/DataModel/LFHResonanceCorrelationTables.h"
#include "PWGLF/DataModel/mcCentrality.h"
#include "PWGLF/Utils/inelGt.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/Zorro.h"
#include "Common/Core/ZorroSummary.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <Framework/ASoA.h>
#include <Framework/ASoAHelpers.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/BinningPolicy.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/StaticFor.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <THn.h>
#include <TList.h>
#include <TPDGCode.h>
#include <TString.h>

#include <fmt/format.h>

#include <Rtypes.h>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <string_view>
#include <variant>
#include <vector>

using namespace o2;
using namespace o2::constants::math;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;

using TracksComplete = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA>;
using AssocPhis = soa::Join<aod::AssocPhis>;
using AssocKstars = soa::Join<aod::AssocKstars>;

struct HResonanceCorrelation {
  // for efficiency corrections if requested
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  Service<o2::framework::O2DatabasePDG> pdgDB;
  o2::pwglf::ParticleCounter<o2::framework::O2DatabasePDG> mCounter;

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // event filtering
  Configurable<std::string> zorroMask{"zorroMask", "", "zorro trigger class to select on (empty: none)"};

  Zorro zorro;
  OutputObj<ZorroSummary> zorroSummary{"zorroSummary"};

  struct : ConfigurableGroup {
    std::string prefix = "masterConfigurations";
    Configurable<bool> doPPAnalysis{"doPPAnalysis", true, "if in pp, set to true"};
    Configurable<int> collisionHasTriggOrAssoc{"collisionHasTriggOrAssoc", 0, "require the collisions containing (0:no requirement 1:trig 2:assoc 3:trig or assoc 4:trig and assoc"};
    Configurable<bool> doFullCorrelationStudy{"doFullCorrelationStudy", true, "if true, do full correlation study by creating all THnSparse histograms for the correlation function"};
    Configurable<bool> doCorrelationHadron{"doCorrelationHadron", false, "do Hadron correlation"};
    Configurable<bool> doCorrelationPhi{"doCorrelationPhi", false, "do Phi correlation"};
    Configurable<bool> doCorrelationKstar{"doCorrelationKstar", false, "do K*0 correlation"};
    Configurable<bool> doCorrelationPion{"doCorrelationPion", false, "do Pion correlation"};
    Configurable<bool> doGenEventSelection{"doGenEventSelection", true, "use event selections when performing closure test for the gen events"};
    Configurable<bool> selectINELgtZERO{"selectINELgtZERO", true, "select INEL>0 events"};
    Configurable<float> zVertexCut{"zVertexCut", 10, "Cut on PV position"};
    Configurable<bool> requireAllGoodITSLayers{"requireAllGoodITSLayers", false, " require that in the event all ITS are good"};
    Configurable<bool> requireGoodTriggerTVX{"requireGoodTriggerTVX", false, " require acceptable FT0C-FT0A time difference"};
    Configurable<bool> requireGoodZvtxFT0vsPV{"requireGoodZvtxFT0vsPV", false, " require small difference between z-vertex from PV and from FT0"};
    Configurable<bool> skipUnderOverflowInTHn{"skipUnderOverflowInTHn", false, "skip under/overflow in THns"};
    Configurable<int> mixingParameter{"mixingParameter", 10, "how many events are mixed"};
    Configurable<bool> doMCassociation{"doMCassociation", false, "fill everything only for MC associated"};
    Configurable<bool> doTriggPhysicalPrimary{"doTriggPhysicalPrimary", false, "require physical primary for trigger particles"};
    Configurable<bool> applyNewMCSelection{"applyNewMCSelection", false, "apply new MC Generated selection"};
    Configurable<bool> doSeparateFT0Prediction{"doSeparateFT0Prediction", false, "separate FT0M to FT0A and FT0C in prediction process"};
    Configurable<bool> useCentralityinPrediction{"useCentralityinPrediction", false, "if true, use centrality instead of multiplisity"};
    Configurable<bool> doMirroringInDelataEta{"doMirroringInDelataEta", false, "if true, fill only positive delta eta and mirror the negative side in post processing, Adjust the delta axis!"};
    Configurable<bool> fillCorrelationHistWithMass{"fillCorrelationHistWithMass", false, "if true, fill correlation histograms with particle mass"};
  } masterConfigurations;

  // master analysis switches
  Configurable<bool> doAssocPhysicalPrimary{"doAssocPhysicalPrimary", false, "require physical primary for associated particles"};
  Configurable<bool> doAssocPhysicalPrimaryInGen{"doAssocPhysicalPrimaryInGen", false, "require physical primary for associated particles in Generated Partilces"};
  Configurable<bool> doAutocorrelationRejection{"doAutocorrelationRejection", true, "reject pairs where trigger Id is the same as daughter particle Id"};
  Configurable<bool> doMixingQAandEventQA{"doMixingQAandEventQA", true, "if true, add EvnetQA and MixingQA hist to histos"};
  Configurable<bool> doITSClustersQA{"doITSClustersQA", true, "if true, add ITSCluster hist to histos"};
  Configurable<bool> doDeltaPhiStarCheck{"doDeltaPhiStarCheck", false, "if true, create and fill delta phi star histograms"};

  Configurable<int> triggerBinToSelect{"triggerBinToSelect", 0, "trigger bin to select on if processSelectEventWithTrigger enabled"};
  Configurable<int> triggerParticleCharge{"triggerParticleCharge", 0, "For checks, if 0 all charged tracks, if -1 only neg., if 1 only positive"};
  Configurable<float> etaSel{"etaSel", 0.8, "Selection in eta for trigger and associated particles"};
  Configurable<float> ySel{"ySel", 0.5, "Selection in rapidity for consistency checks"};

  Configurable<bool> useTheLeadingParticleAsTrigger{"useTheLeadingParticleAsTrigger", false, "if true, use the leading particle in the event as trigger particle"};
  // used for event selections in Pb-Pb
  Configurable<int> cfgCutOccupancyHigh{"cfgCutOccupancyHigh", 3000, "High cut on TPC occupancy"};
  Configurable<int> cfgCutOccupancyLow{"cfgCutOccupancyLow", 0, "Low cut on TPC occupancy"};

  // Axes - configurable for smaller sizes
  struct : ConfigurableGroup {
    std::string prefix = "axesConfigurations";
    ConfigurableAxis axisMult{"axisMult", {VARIABLE_WIDTH, 0.0f, 0.01f, 1.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 70.0f, 100.0f}, "Mixing bins - multiplicity"};
    ConfigurableAxis axisVtxZ{"axisVtxZ", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
    ConfigurableAxis axisPhi{"axisPhi", {72, 0, TwoPI}, "#phi"};
    ConfigurableAxis axisEta{"axisEta", {80, -0.8, +0.8}, "#eta"};
    ConfigurableAxis axisDeltaPhi{"axisDeltaPhi", {72, -PIHalf, PIHalf * 3}, "delta #varphi axis for histograms"};
    ConfigurableAxis axisDeltaEta{"axisDeltaEta", {50, -1.6, 1.6}, "delta eta axis for histograms"};
    ConfigurableAxis axisPtAssoc{"axisPtAssoc", {VARIABLE_WIDTH, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 10.0}, "pt associated axis for histograms"};
    ConfigurableAxis axisPtTrigger{"axisPtTrigger", {VARIABLE_WIDTH, 0.0, 1.0, 2.0, 3.0, 100}, "pt associated axis for histograms"};
    ConfigurableAxis axisPtQA{"axisPtQA", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 17.0f, 19.0f, 21.0f, 23.0f, 25.0f, 30.0f, 35.0f, 40.0f, 50.0f}, "pt axis for QA histograms"};
    ConfigurableAxis axisMassNSigma{"axisMassNSigma", {40, -2, 2}, "Axis for mass Nsigma"};
    ConfigurableAxis axisPhiMass{"axisPhiMass", {300, 1.01f, 1.31f}, "Inv. Mass (GeV/c^{2})"};
    ConfigurableAxis axisKstarMass{"axisKstarMass", {300, 0.596f, 1.196f}, "Inv. Mass (GeV/c^{2})"};
    ConfigurableAxis axisMultiplicity{"axisMultiplicity", {VARIABLE_WIDTH, 0, 20, 40, 60, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300}, "Binning of the Multiplicity axis in model prediction process"};
    ConfigurableAxis axisMidrapidityMultiplicity{"axisMidrapidityMultiplicity", {VARIABLE_WIDTH, 0, 20, 40, 60, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300}, "Binning of the Midrapidity Multiplicity axis in model prediction process"};
  } axesConfigurations;

  // for topo var QA
  // Per-species, peak-relative mass windows: a single shared
  // `massWindowConfigurations` comparing raw assocCandidate.mass() directly
  // against these bounds would break isLeftBg structurally (it requires a
  // NEGATIVE mass, which raw mass never is), and sharing one window between
  // Phi and K*0 leaves one species' tuned range (Phi's, ~1.0-1.04 GeV)
  // entirely outside the other's candidate mass range (K*0's filter-level
  // cut is 0.796-0.996 GeV), making K*0's RightBg structurally unreachable.
  //
  // Instead, classify using delta = mass - peakMass (a real signed quantity,
  // so isLeftBg is reachable), with each species' own peak mass and its own
  // effective sigma, in its own struct. maxPeakNSigma/minBgNSigma/maxBgNSigma
  // are genuinely "number of sigma" (2, 3, 6 by default), multiplied by
  // `sigma` (in GeV) to get the actual mass-offset bounds.
  //
  // `sigma` defaults here are PDG Gamma/2.355 (treating the natural
  // Breit-Wigner width as if it were a Gaussian FWHM) -- a reasonable
  // placeholder for K*0, where the ~47 MeV natural width likely dominates
  // over detector resolution, but LIKELY AN UNDERESTIMATE for Phi, whose
  // ~4.2 MeV natural width is probably smaller than your actual track
  // momentum resolution. Refit both from your own reconstructed invariant
  // mass peak (Voigt / BW (x) Gaussian) once you have statistics, and update
  // `sigma` here -- everything else (the Nsigma multipliers) stays valid.
  struct : ConfigurableGroup {
    std::string prefix = "massWindowConfigurationsPhi";
    Configurable<float> peakMass{"peakMass", 1.019455f, "Phi(1020) PDG mass (GeV), used as delta = mass - peakMass"};
    Configurable<float> sigma{"sigma", 0.0018f, "effective width (GeV) for the Nsigma windows below -- PDG Gamma/2.355 placeholder, refit from your own peak"};
    Configurable<float> maxPeakNSigma{"maxPeakNSigma", 2, "Signal region half-width, in units of sigma"};
    Configurable<float> minBgNSigma{"minBgNSigma", 3, "Bg region edge closest to peak, in units of sigma"};
    Configurable<float> maxBgNSigma{"maxBgNSigma", 6, "Bg region edge furthest from peak, in units of sigma"};
  } massWindowConfigurationsPhi; // allows for gap between peak and bg in case someone wants to

  struct : ConfigurableGroup {
    std::string prefix = "massWindowConfigurationsKstar";
    Configurable<float> peakMass{"peakMass", 0.89555f, "K*0(892) PDG mass (GeV), used as delta = mass - peakMass"};
    Configurable<float> sigma{"sigma", 0.0201f, "effective width (GeV) for the Nsigma windows below -- PDG Gamma/2.355 placeholder, refit from your own peak"};
    Configurable<float> maxPeakNSigma{"maxPeakNSigma", 2, "Signal region half-width, in units of sigma"};
    Configurable<float> minBgNSigma{"minBgNSigma", 3, "Bg region edge closest to peak, in units of sigma"};
    Configurable<float> maxBgNSigma{"maxBgNSigma", 6, "Bg region edge furthest from peak, in units of sigma"};
  } massWindowConfigurationsKstar; // allows for gap between peak and bg in case someone wants to

  // Implementation of on-the-spot efficiency correction
  struct : ConfigurableGroup {
    std::string prefix = "efficiencyFlags";
    Configurable<bool> applyEfficiencyCorrection{"applyEfficiencyCorrection", false, "apply efficiency correction"};
    Configurable<bool> applyEfficiencyForTrigger{"applyEfficiencyForTrigger", false, "apply efficiency correction for the trigger particle"};
    Configurable<bool> applyEfficiencyPropagation{"applyEfficiencyPropagation", false, "propagate also the efficiency uncertainty"};
    Configurable<bool> applyPurityHadron{"applyPurityHadron", false, "apply the purity correction for associated hadrons"};
    Configurable<bool> applyPurityTrigger{"applyPurityTrigger", false, "apply the purity correction for trigger particle"};
    Configurable<bool> applyEffAsFunctionOfMult{"applyEffAsFunctionOfMult", false, "apply efficiency as a function of multiplicity as well"};
    Configurable<bool> applyEffAsFunctionOfMultAndPhi{"applyEffAsFunctionOfMultAndPhi", false, "apply efficiency as a function of multiplicity and phi"};
  } efficiencyFlags;
  Configurable<std::string> ccdburl{"ccdburl", "http://alice-ccdb.cern.ch", "url of the ccdb repository to use"};
  Configurable<std::string> efficiencyCCDBPath{"efficiencyCCDBPath", "GLO/Config/GeometryAligned", "Path of the efficiency corrections"};

  // Configurables for doing subwagon systematics
  struct : ConfigurableGroup {
    std::string prefix = "trackSelection";
    // --- Track quality variations (single track, both trigger and assoc daughters)
    Configurable<int> minTPCNCrossedRowsTrigger{"minTPCNCrossedRowsTrigger", 70, "Minimum TPC crossed rows (trigger)"};
    Configurable<int> minTPCNCrossedRowsAssociated{"minTPCNCrossedRowsAssociated", 70, "Minimum TPC crossed rows (associated)"};
    Configurable<bool> triggerRequireITS{"triggerRequireITS", true, "require ITS signal in trigger tracks"};
    Configurable<bool> assocRequireITS{"assocRequireITS", true, "require ITS signal in associated primary tracks"};
    Configurable<int> triggerMaxTPCSharedClusters{"triggerMaxTPCSharedClusters", 200, "maximum number of shared TPC clusters (inclusive)"};
    Configurable<int> assocMaxTPCSharedClusters{"assocMaxTPCSharedClusters", 200, "maximum number of shared TPC clusters (inclusive) for assoc primary tracks"};
    Configurable<bool> triggerRequireL0{"triggerRequireL0", false, "require ITS L0 cluster for trigger"};
    Configurable<bool> assocRequireL0{"assocRequireL0", true, "require ITS L0 cluster for assoc primary track"};
    Configurable<bool> requireDCAzCut{"requireDCAzCut", false, "require DCAz cut for trigger and associated primary tracks"};

    // --- Trigger: DCA variation from basic formula: |DCAxy| <  0.004f + (0.013f / pt)
    Configurable<float> dcaXYconstant{"dcaXYconstant", 0.004, "[0] in |DCAxy| < [0]+[1]/pT"};
    Configurable<float> dcaXYpTdep{"dcaXYpTdep", 0.013, "[1] in |DCAxy| < [0]+[1]/pT"};
    // --- Assoc track: DCA variation from basic formula: |DCAxy| <  0.004f + (0.013f / pt)
    Configurable<float> dcaXYconstantAssoc{"dcaXYconstantAssoc", 0.004, "[0] in |DCAxy| < [0]+[1]/pT"};
    Configurable<float> dcaXYpTdepAssoc{"dcaXYpTdepAssoc", 0.013, "[1] in |DCAxy| < [0]+[1]/pT"};

    Configurable<float> dcaZconstant{"dcaZconstant", 0.004, "[0] in |DCAz| < [0]+[1]/pT"};
    Configurable<float> dcaZpTdep{"dcaZpTdep", 0.013, "[1] in |DCAz| < [0]+[1]/pT"};
    Configurable<float> dcaZconstantAssoc{"dcaZconstantAssoc", 0.004, "[0] in |DCAz| < [0]+[1]/pT"};
    Configurable<float> dcaZpTdepAssoc{"dcaZpTdepAssoc", 0.013, "[1] in |DCAz| < [0]+[1]/pT"};
  } trackSelection;

  struct : ConfigurableGroup {
    std::string prefix = "checks";

    // on the fly correction instead of mixingParameter
    Configurable<bool> doOnTheFlyFlattening{"doOnTheFlyFlattening", 0, "enable an on-the-fly correction instead of using mixing"};
  } checks;

  struct ValidCollision {
    struct ValidParticle {
      float eta;
      float phi;
      float pt;
      int region;
      float efficiency;
      float efficiencyError;
      // no species tag stored here: the `type` argument to addValidParticle
      // below (-1 = trigger, else Phi/K*0) only decides which vector this
      // particle goes into; once stored, nothing reads it back.
    };
    float pvz = 0.f;
    float mult = 0.f;
    std::vector<ValidParticle> trigParticles;
    std::vector<ValidParticle> assocParticles;
    void addValidParticle(float eta, float phi, float pt, int region, float efficiency, float efficiencyError, int type)
    {
      ValidParticle particle{eta, phi, pt, region, efficiency, efficiencyError};

      if (type == -1) {
        trigParticles.push_back(particle);
      } else {
        assocParticles.push_back(particle);
      }
    }
  };

  using ValidCollisions = std::vector<std::vector<ValidCollision>>;
  ValidCollisions validCollisions;

  // objects to use for efficiency corrections
  TH2F* hEfficiencyTrigger = nullptr;
  TH3F* hEfficiencyTriggerMult = nullptr;
  THnF* hEfficiencyTriggerMultVsPhi = nullptr;
  TH2F* hEfficiencyPion = nullptr;
  TH2F* hEfficiencyPhi = nullptr;
  THnF* hEfficiencyPhiMultVsPhi = nullptr;
  TH2F* hEfficiencyKstar = nullptr;
  THnF* hEfficiencyKstarMultVsPhi = nullptr;
  TH2F* hEfficiencyHadron = nullptr;
  TH3F* hEfficiencyHadronMult = nullptr;
  TH1F* hPurityHadron = nullptr;
  TH2F* hPurityHadronMult = nullptr;

  // objects to propagate the efficiency uncertainty
  TH2F* hEfficiencyUncertaintyTrigger = nullptr;
  TH3F* hEfficiencyUncertaintyTriggerMult = nullptr;
  TH2F* hEfficiencyUncertaintyPion = nullptr;
  TH2F* hEfficiencyUncertaintyPhi = nullptr;
  TH2F* hEfficiencyUncertaintyKstar = nullptr;
  TH2F* hEfficiencyUncertaintyHadron = nullptr;
  TH3F* hEfficiencyUncertaintyHadronMult = nullptr;
  TH1F* hPurityUncertaintyHadron = nullptr;
  TH2F* hPurityUncertaintyHadronMult = nullptr;

  using BinningTypePP = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0M>;
  using BinningTypePbPb = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0C>;

  // collision slicing for mixed events
  Preslice<aod::TriggerTracks> collisionSliceTracks = aod::triggerTracks::collisionId;
  // Preslice<aod::AssocPions> collisionSlicePions = aod::assocHadrons::collisionId;
  Preslice<aod::AssocHadrons> collisionSliceHadrons = aod::assocHadrons::collisionId;
  Preslice<aod::McParticles> perCollision = aod::mcparticle::mcCollisionId;
  Preslice<aod::AssocPhis> collisionSlicePhis = aod::assocPhis::collisionId;
  Preslice<aod::AssocKstars> collisionSliceKstars = aod::assocKstars::collisionId;

  static constexpr std::string_view Particlenames[] = {"Phi", "Kstar0", "Pion", "Hadron"};
  static constexpr int PdgCodes[] = {333, 313, 211, 0}; // Hadron has no single PDG code; 0 is a harmless placeholder

  static constexpr int IndexPhi = 0;
  static constexpr int IndexKstar = 1;
  static constexpr int IndexPion = 2;

  uint16_t doCorrelation = 0;
  int mRunNumber = 0;
  int mRunNumberZorro = 0;

  std::vector<std::vector<float>> axisRanges;

  static constexpr float MinRadiusTPC = 0.8;
  static constexpr float MaxRadiusTPC = 2.5;

  static constexpr float Neutral = 0.0;

  static constexpr int AssocParticleTypes = 4;         // Phi, Kstar0, Pion, Hadron
  static constexpr int AssocParticleTypesNoHadron = 3; // Phi, Kstar0, Pion

  /// Function to aid in calculating delta-phi
  /// \param phi1 first phi value
  /// \param phi2 second phi value
  double computeDeltaPhi(double phi1, double phi2)
  {
    double deltaPhi = phi1 - phi2;
    double shiftedDeltaPhi = RecoDecay::constrainAngle(deltaPhi, -PIHalf);
    return shiftedDeltaPhi;
  }

  /// Function to load zorro
  /// \param bc provided such that the run number + timestamp can be used
  void initZorro(aod::BCsWithTimestamps::iterator const& bc)
  {
    if (mRunNumberZorro == bc.runNumber()) {
      return;
    }

    zorro.initCCDB(ccdb.service, bc.runNumber(), bc.timestamp(), zorroMask.value);
    zorro.populateHistRegistry(histos, bc.runNumber());

    mRunNumberZorro = bc.runNumber();
  }

  /// Function to load efficiencies to memory from CCDB
  /// \param bc provided such that the run number can be used
  void initEfficiencyFromCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }
    mRunNumber = bc.runNumber();
    LOG(info) << "Loading efficiencies from CCDB for run " << mRunNumber << " now...";
    auto timeStamp = bc.timestamp();

    TList* listEfficiencies = ccdb->getForTimeStamp<TList>(efficiencyCCDBPath, timeStamp);

    if (!listEfficiencies) {
      LOG(fatal) << "Problem getting TList object with efficiencies!";
    }

    hEfficiencyTrigger = static_cast<TH2F*>(listEfficiencies->FindObject("hEfficiencyTrigger"));
    hEfficiencyTriggerMult = static_cast<TH3F*>(listEfficiencies->FindObject("hEfficiencyTriggerMult"));
    hEfficiencyTriggerMultVsPhi = static_cast<THnF*>(listEfficiencies->FindObject("hEfficiencyTriggerMultVsPhi"));
    hEfficiencyPhi = static_cast<TH2F*>(listEfficiencies->FindObject("hEfficiencyPhi"));
    hEfficiencyPhiMultVsPhi = static_cast<THnF*>(listEfficiencies->FindObject("hEfficiencyPhiMultVsPhi"));
    hEfficiencyKstar = static_cast<TH2F*>(listEfficiencies->FindObject("hEfficiencyKstar"));
    hEfficiencyKstarMultVsPhi = static_cast<THnF*>(listEfficiencies->FindObject("hEfficiencyKstarMultVsPhi"));
    hEfficiencyHadron = static_cast<TH2F*>(listEfficiencies->FindObject("hEfficiencyHadron"));
    hEfficiencyHadronMult = static_cast<TH3F*>(listEfficiencies->FindObject("hEfficiencyHadronMult"));
    hEfficiencyPion = static_cast<TH2F*>(listEfficiencies->FindObject("hEfficiencyPion"));
    hPurityHadron = static_cast<TH1F*>(listEfficiencies->FindObject("hPurityHadron"));
    hPurityHadronMult = static_cast<TH2F*>(listEfficiencies->FindObject("hPurityHadronMult"));
    hEfficiencyUncertaintyTrigger = static_cast<TH2F*>(listEfficiencies->FindObject("hEfficiencyUncertaintyTrigger"));
    hEfficiencyUncertaintyTriggerMult = static_cast<TH3F*>(listEfficiencies->FindObject("hEfficiencyUncertaintyTriggerMult"));
    hEfficiencyUncertaintyPhi = static_cast<TH2F*>(listEfficiencies->FindObject("hEfficiencyUncertaintyPhi"));
    hEfficiencyUncertaintyKstar = static_cast<TH2F*>(listEfficiencies->FindObject("hEfficiencyUncertaintyKstar"));
    hEfficiencyUncertaintyPion = static_cast<TH2F*>(listEfficiencies->FindObject("hEfficiencyUncertaintyPion"));
    hEfficiencyUncertaintyHadron = static_cast<TH2F*>(listEfficiencies->FindObject("hEfficiencyUncertaintyHadron"));
    hEfficiencyUncertaintyHadronMult = static_cast<TH3F*>(listEfficiencies->FindObject("hEfficiencyUncertaintyHadronMult"));
    hPurityUncertaintyHadron = static_cast<TH1F*>(listEfficiencies->FindObject("hPurityUncertaintyHadron"));
    hPurityUncertaintyHadronMult = static_cast<TH2F*>(listEfficiencies->FindObject("hPurityUncertaintyHadronMult"));
    if (efficiencyFlags.applyEfficiencyPropagation && !efficiencyFlags.applyEffAsFunctionOfMultAndPhi && !hEfficiencyUncertaintyTrigger)
      LOG(fatal) << "Problem getting hEfficiencyUncertaintyTrigger!";
    if (efficiencyFlags.applyEffAsFunctionOfMult && !hEfficiencyTriggerMult)
      LOG(fatal) << "Problem getting hEfficiencyTriggerMult!";
    LOG(info) << "Efficiencies now loaded for " << mRunNumber;
  }

  void init(InitContext const&)
  {
    zorroSummary.setObject(zorro.getZorroSummary());
    mRunNumber = 0;
    mRunNumberZorro = 0;
    hEfficiencyPion = 0x0;
    hEfficiencyPhi = 0x0;
    hEfficiencyKstar = 0x0;
    hEfficiencyUncertaintyTrigger = 0x0;
    hEfficiencyUncertaintyPion = 0x0;
    hEfficiencyUncertaintyPhi = 0x0;
    hEfficiencyUncertaintyKstar = 0x0;

    hEfficiencyHadron = 0x0;
    hPurityHadron = 0x0;
    hPurityUncertaintyHadron = 0x0;
    hEfficiencyUncertaintyHadron = 0x0;

    // set bitmap for convenience
    doCorrelation = 0;
    if (masterConfigurations.doCorrelationPhi)
      SETBIT(doCorrelation, IndexPhi);
    if (masterConfigurations.doCorrelationKstar)
      SETBIT(doCorrelation, IndexKstar);
    if (masterConfigurations.doCorrelationPion)
      SETBIT(doCorrelation, IndexPion);
    if (masterConfigurations.doCorrelationHadron)
      SETBIT(doCorrelation, 3);

    // Store axis ranges to prevent spurious filling
    // axis status:
    // --- Delta-phi is safe -> math forbids insanity
    // --- Delta-eta depends on pre-filter -> check
    // --- pT assoc depends on binning -> check
    // --- vertex Z is safe -> skipped at evsel level
    // --- multiplicity -> check

    // grab axis edge from ConfigurableAxes
    const AxisSpec preAxisDeltaPhi{axesConfigurations.axisDeltaPhi, "#Delta#varphi"};
    const AxisSpec preAxisDeltaEta{axesConfigurations.axisDeltaEta, "#Delta#eta"};
    const AxisSpec preAxisPtAssoc{axesConfigurations.axisPtAssoc, "#it{p}_{T}^{assoc} (GeV/c)"};
    const AxisSpec preAxisPtTrigger{axesConfigurations.axisPtTrigger, "#it{p}_{T}^{trigger} (GeV/c)"};
    const AxisSpec preAxisVtxZ{axesConfigurations.axisVtxZ, "vertex Z (cm)"};
    const AxisSpec preAxisMult{axesConfigurations.axisMult, "mult percentile"};
    const AxisSpec axisPtPhi{axesConfigurations.axisPtAssoc, "#it{p}_{T}^{#phi} (GeV/c)"};
    const AxisSpec preAxisMultiplicity{axesConfigurations.axisMultiplicity, "multiplicity"};

    // store the original axes in specific TH1Cs for completeness
    histos.add("axes/hDeltaPhiAxis", "", kTH1C, {preAxisDeltaPhi});
    histos.add("axes/hDeltaEtaAxis", "", kTH1C, {preAxisDeltaEta});
    histos.add("axes/hPtAssocAxis", "", kTH1C, {preAxisPtAssoc});
    histos.add("axes/hPtTriggerAxis", "", kTH1C, {preAxisPtTrigger});
    histos.add("axes/hVertexZAxis", "", kTH1C, {preAxisVtxZ});
    histos.add("axes/hMultAxis", "", kTH1C, {preAxisMult});
    histos.add("axes/hMultiplicityAxis", "", kTH1C, {preAxisMultiplicity});

    std::vector<double> edgesDeltaPhiOrig = preAxisDeltaPhi.binEdges;
    std::vector<double> edgesDeltaEtaOrig = preAxisDeltaEta.binEdges;
    std::vector<double> edgesPtAssocOrig = preAxisPtAssoc.binEdges;
    std::vector<double> edgesPtTriggerOrig = preAxisPtTrigger.binEdges;
    std::vector<double> edgesVtxZOrig = preAxisVtxZ.binEdges;
    std::vector<double> edgesMultOrig = preAxisMult.binEdges;
    std::vector<double> edgesMultiplicityOrig = preAxisMultiplicity.binEdges;

    std::vector<float> rangesDeltaPhi = {static_cast<float>(edgesDeltaPhiOrig[0]), static_cast<float>(edgesDeltaPhiOrig[edgesDeltaPhiOrig.size() - 1])};
    std::vector<float> rangesDeltaEta = {static_cast<float>(edgesDeltaEtaOrig[0]), static_cast<float>(edgesDeltaEtaOrig[edgesDeltaEtaOrig.size() - 1])};
    std::vector<float> rangesPtAssoc = {static_cast<float>(edgesPtAssocOrig[0]), static_cast<float>(edgesPtAssocOrig[edgesPtAssocOrig.size() - 1])};
    std::vector<float> rangesPtTrigger = {static_cast<float>(edgesPtTriggerOrig[0]), static_cast<float>(edgesPtTriggerOrig[edgesPtTriggerOrig.size() - 1])};
    std::vector<float> rangesVtxZ = {static_cast<float>(edgesVtxZOrig[0]), static_cast<float>(edgesVtxZOrig[edgesVtxZOrig.size() - 1])};
    std::vector<float> rangesMult = {static_cast<float>(edgesMultOrig[0]), static_cast<float>(edgesMultOrig[edgesMultOrig.size() - 1])};
    std::vector<float> rangesMultiplicity = {static_cast<float>(edgesMultiplicityOrig[0]), static_cast<float>(edgesMultiplicityOrig[edgesMultiplicityOrig.size() - 1])};

    axisRanges.emplace_back(rangesDeltaPhi);
    axisRanges.emplace_back(rangesDeltaEta);
    axisRanges.emplace_back(rangesPtAssoc);
    axisRanges.emplace_back(rangesPtTrigger);
    axisRanges.emplace_back(rangesVtxZ);
    axisRanges.emplace_back(rangesMult);
    axisRanges.emplace_back(rangesMultiplicity);

    std::vector<double> edgesDeltaPhi;
    std::vector<double> edgesDeltaEta;
    std::vector<double> edgesPtAssoc;
    std::vector<double> edgesPtTrigger;
    std::vector<double> edgesVtxZ;
    std::vector<double> edgesMult;
    std::vector<double> edgesMultiplicity;

    // v--- skipUnderOverflowInTHn ---v
    //
    // if enabled, this will change the axes such that they will solely cover the interval from
    // edge[1] to edge[n-1]; this will mean that the bin 1 and bin N will be stored in
    // under / overflow bins and will have to be manually unpacked. Do not forget to do the manual
    // unpacking a posteriori!
    //
    // this feature is meant to save memory conveniently.
    // it should actually be implemented centrally in ROOT but ok, this will do it for now.

    int offset = masterConfigurations.skipUnderOverflowInTHn ? 1 : 0;
    // ===] delta-phi [===
    if (!preAxisDeltaPhi.nBins.has_value()) {
      // variable binning, use bins provided
      for (int i = offset; i < static_cast<int>(edgesDeltaPhiOrig.size()) - offset; i++)
        edgesDeltaPhi.emplace_back(edgesDeltaPhiOrig[i]);
    } else {
      // fixed binning, generate the bin edges on-the-spot
      double min = edgesDeltaPhiOrig[0];
      double delta = (edgesDeltaPhiOrig[1] - edgesDeltaPhiOrig[0]) / preAxisDeltaPhi.nBins.value();
      for (int i = offset; i < preAxisDeltaPhi.nBins.value() + 1 - offset; i++)
        edgesDeltaPhi.emplace_back(min + static_cast<double>(i) * delta);
    }
    // ===] delta-eta [===
    if (!preAxisDeltaEta.nBins.has_value()) {
      // variable binning, use bins provided
      for (int i = offset; i < static_cast<int>(edgesDeltaEtaOrig.size()) - offset; i++)
        edgesDeltaEta.emplace_back(edgesDeltaEtaOrig[i]);
    } else {
      // fixed binning, generate the bin edges on-the-spot
      double min = edgesDeltaEtaOrig[0];
      double delta = (edgesDeltaEtaOrig[1] - edgesDeltaEtaOrig[0]) / preAxisDeltaEta.nBins.value();
      for (int i = offset; i < preAxisDeltaEta.nBins.value() + 1 - offset; i++)
        edgesDeltaEta.emplace_back(min + static_cast<double>(i) * delta);
    }
    // ===] pt assoc [===
    if (!preAxisPtAssoc.nBins.has_value()) {
      // variable binning, use bins provided
      for (int i = offset; i < static_cast<int>(edgesPtAssocOrig.size()) - offset; i++)
        edgesPtAssoc.emplace_back(edgesPtAssocOrig[i]);
    } else {
      // fixed binning, generate the bin edges on-the-spot
      double min = edgesPtAssocOrig[0];
      double delta = (edgesPtAssocOrig[1] - edgesPtAssocOrig[0]) / preAxisPtAssoc.nBins.value();
      for (int i = offset; i < preAxisPtAssoc.nBins.value() + 1 - offset; i++)
        edgesPtAssoc.emplace_back(min + static_cast<double>(i) * delta);
    }
    // ===] pt trigger [===
    if (!preAxisPtTrigger.nBins.has_value()) {
      // variable binning, use bins provided
      for (int i = offset; i < static_cast<int>(edgesPtTriggerOrig.size()) - offset; i++)
        edgesPtTrigger.emplace_back(edgesPtTriggerOrig[i]);
    } else {
      // fixed binning, generate the bin edges on-the-spot
      double min = edgesPtTriggerOrig[0];
      double delta = (edgesPtTriggerOrig[1] - edgesPtTriggerOrig[0]) / preAxisPtTrigger.nBins.value();
      for (int i = offset; i < preAxisPtTrigger.nBins.value() + 1 - offset; i++)
        edgesPtTrigger.emplace_back(min + static_cast<double>(i) * delta);
    }
    // ===] vtx Z [===
    if (!preAxisVtxZ.nBins.has_value()) {
      // variable binning, use bins provided
      for (int i = offset; i < static_cast<int>(edgesVtxZOrig.size()) - offset; i++)
        edgesVtxZ.emplace_back(edgesVtxZOrig[i]);
    } else {
      // fixed binning, generate the bin edges on-the-spot
      double min = edgesVtxZOrig[0];
      double delta = (edgesVtxZOrig[1] - edgesVtxZOrig[0]) / preAxisVtxZ.nBins.value();
      for (int i = offset; i < preAxisVtxZ.nBins.value() + 1 - offset; i++)
        edgesVtxZ.emplace_back(min + static_cast<double>(i) * delta);
    }
    // ===] mult percentile [===
    if (!preAxisMult.nBins.has_value()) {
      // variable binning, use bins provided
      for (int i = offset; i < static_cast<int>(edgesMultOrig.size()) - offset; i++)
        edgesMult.emplace_back(edgesMultOrig[i]);
    } else {
      // fixed binning, generate the bin edges on-the-spot
      double min = edgesMultOrig[0];
      double delta = (edgesMultOrig[1] - edgesMultOrig[0]) / preAxisMult.nBins.value();
      for (int i = offset; i < preAxisMult.nBins.value() + 1 - offset; i++)
        edgesMult.emplace_back(min + static_cast<double>(i) * delta);
    }
    // ===] multiplicity count [===
    if (!preAxisMultiplicity.nBins.has_value()) {
      // variable binning, use bins provided
      for (int i = offset; i < static_cast<int>(edgesMultiplicityOrig.size()) - offset; i++)
        edgesMultiplicity.emplace_back(edgesMultiplicityOrig[i]);
    } else {
      // fixed binning, generate the bin edges on-the-spot
      double min = edgesMultiplicityOrig[0];
      double delta = (edgesMultiplicityOrig[1] - edgesMultiplicityOrig[0]) / preAxisMultiplicity.nBins.value();
      for (int i = offset; i < preAxisMultiplicity.nBins.value() + 1 - offset; i++)
        edgesMultiplicity.emplace_back(min + static_cast<double>(i) * delta);
    }

    LOGF(info, "Initialized THnF axis delta-phi with %i bins.", edgesDeltaPhi.size() - 1);
    LOGF(info, "Initialized THnF axis delta-eta with %i bins.", edgesDeltaEta.size() - 1);
    LOGF(info, "Initialized THnF axis pTassoc with %i bins.", edgesPtAssoc.size() - 1);
    LOGF(info, "Initialized THnF axis pTtrigger with %i bins.", edgesPtTrigger.size() - 1);
    LOGF(info, "Initialized THnF axis vertex-Z with %i bins.", edgesVtxZ.size() - 1);
    LOGF(info, "Initialized THnF axis mult cent with %i bins.", edgesMult.size() - 1);
    LOGF(info, "Initialized THnF axis multiplicity with %i bins.", edgesMultiplicity.size() - 1);

    const AxisSpec axisDeltaPhiNDim{edgesDeltaPhi, "#Delta#varphi"};
    const AxisSpec axisDeltaEtaNDim{edgesDeltaEta, "#Delta#eta"};
    const AxisSpec axisPtAssocNDim{edgesPtAssoc, "#it{p}_{T}^{assoc} (GeV/c)"};
    const AxisSpec axisPtTriggerNDim{edgesPtTrigger, "#it{p}_{T}^{trigger} (GeV/c)"};
    const AxisSpec axisVtxZNDim{edgesVtxZ, "vertex Z (cm)"};
    const AxisSpec axisMultNDim{edgesMult, "mult percentile"};
    const AxisSpec axisMultiplicityNDim{edgesMultiplicity, "Multiplicity"};

    if (doprocessMixedEventHPhisInBuffer || doprocessMixedEventHKstarsInBuffer) {
      validCollisions.resize(histos.get<TH1>(HIST("axes/hMultAxis"))->GetNbinsX() * histos.get<TH1>(HIST("axes/hVertexZAxis"))->GetNbinsX());
      for (size_t i = 0; i < validCollisions.size(); ++i) {
        validCollisions[i].reserve(masterConfigurations.mixingParameter);
      }
    }
    if (!masterConfigurations.doPPAnalysis) {
      // event selections in Pb-Pb
      histos.add("hEventSelection", "hEventSelection", kTH1F, {{10, 0, 10}});
      TString eventSelLabel[] = {"all", "sel8", "kIsTriggerTVX", "PV_{z}", "kIsGoodITSLayersAll", "kIsGoodZvtxFT0vsPV", "OccupCut", "kNoTimeFrameBorder", "kNoITSROFrameBorder", "kNoSameBunchPileup "};
      for (int i = 1; i <= histos.get<TH1>(HIST("hEventSelection"))->GetNbinsX(); i++) {
        histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(i, eventSelLabel[i - 1]);
      }
    }

    // ========================================================================
    // SHARED QA HISTOGRAMS
    // ========================================================================
    bool anySameEvent = doprocessSameEventHPhis || doprocessSameEventHKstars || doprocessSameEventHPions || doprocessSameEventHHadrons;
    bool anyMixedEvent = doprocessMixedEventHPhis || doprocessMixedEventHPhisInBuffer || doprocessMixedEventHKstars || doprocessMixedEventHKstarsInBuffer || doprocessMixedEventHPions || doprocessMixedEventHHadrons;

    if (doMixingQAandEventQA && (anySameEvent || anyMixedEvent)) {
      if (anySameEvent)
        histos.add("MixingQA/hSECollisionBins", ";bin;Entries", kTH1F, {{140, -0.5, 139.5}});
      if (anyMixedEvent) {
        histos.add("MixingQA/hMECollisionBins", ";bin;Entries", kTH1F, {{140, -0.5, 139.5}});
        histos.add("MixingQA/hMEpvz1", ";pvz;Entries", kTH1F, {{30, -15, 15}});
        histos.add("MixingQA/hMEpvz2", ";pvz;Entries", kTH1F, {{30, -15, 15}});
        histos.add("MixingQA/hMixingQA", "mixing QA", kTH1F, {{2, -0.5, 1.5}});
      }
      histos.add("EventQA/hMult", "Multiplicity", kTH1F, {axesConfigurations.axisMult});
      histos.add("EventQA/hPvz", ";pvz;Entries", kTH1F, {{30, -15, 15}});
    }

    if (anySameEvent) {
      histos.add("hTriggerAllSelectedEtaVsPt", "hTriggerAllSelectedEtaVsPt", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisMult});
      histos.add("hTriggerPtResolution", ";p_{T}^{reconstructed} (GeV/c); p_{T}^{generated} (GeV/c)", kTH2F, {axesConfigurations.axisPtQA, axesConfigurations.axisPtQA});
      histos.add("hTrackEtaVsPtVsPhi", "hTrackEtaVsPtVsPhi", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisPhi});
      if (!efficiencyFlags.applyEffAsFunctionOfMultAndPhi) {
        histos.add("hTriggerPrimaryEtaVsPt", "hTriggerPrimaryEtaVsPt", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisMult});
      } else {
        histos.add("hTriggerPrimaryEtaVsPt", "hTriggerPrimaryEtaVsPt", kTHnF, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisPhi, axesConfigurations.axisMult});
      }
    }

    // ========================================================================
    // MAIN CORRELATION LOOP (SE & ME)
    // ========================================================================
    bool needsSignalBgCloneSE = false;
    bool needsSignalBgCloneME = false;
    for (int i = 0; i < AssocParticleTypes; i++) {
      if (TESTBIT(doCorrelation, i)) {

        bool doSE = false;
        bool doME = false;
        if (i == IndexPhi) { // Phi
          doSE = doprocessSameEventHPhis;
          doME = doprocessMixedEventHPhis || doprocessMixedEventHPhisInBuffer;
        } else if (i == IndexKstar) { // K*0
          doSE = doprocessSameEventHKstars;
          doME = doprocessMixedEventHKstars || doprocessMixedEventHKstarsInBuffer;
        } else if (i == IndexPion) { // Pion
          doSE = doprocessSameEventHPions;
          doME = doprocessMixedEventHPions;
        } else { // Hadron
          doSE = doprocessSameEventHHadrons;
          doME = doprocessMixedEventHHadrons;
        }

        // Kinematic QA (Only for Phi)
        if (i < IndexPion && (doSE || doME || doprocessPrediction || doprocessClosureTest)) {
          if (!efficiencyFlags.applyEffAsFunctionOfMultAndPhi) {
            histos.add(fmt::format("h{}EtaVsPtVsPhi", Particlenames[i]).c_str(), "", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisPhi});
            histos.add(fmt::format("h{}EtaVsPtVsPhiBg", Particlenames[i]).c_str(), "", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisPhi});
          } else {
            histos.add(fmt::format("h{}EtaVsPtVsPhiVsCent", Particlenames[i]).c_str(), "", kTHnF, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisPhi, axesConfigurations.axisMult});
            histos.add(fmt::format("h{}EtaVsPtVsPhiVsCentBg", Particlenames[i]).c_str(), "", kTHnF, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisPhi, axesConfigurations.axisMult});
          }
          histos.add(fmt::format("h3d{}Spectrum", Particlenames[i]).c_str(), fmt::format("h3d{}Spectrum", Particlenames[i]).c_str(), kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisMult, axesConfigurations.axisMassNSigma});
          histos.add(fmt::format("h3d{}SpectrumY", Particlenames[i]).c_str(), fmt::format("h3d{}SpectrumY", Particlenames[i]).c_str(), kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisMult, axesConfigurations.axisMassNSigma});

          // REMOVED (2026-09-24): hITSClusters{Phi,Kstar0}{Positive,Negative}Daughter{Toward,Transverse}.
          // These were registered but never filled anywhere in this file (confirmed by
          // exhaustive grep for any fill using these names). The V0/cascade analog in
          // hStrangeCorrelation.cxx fills the equivalent histograms with the DAUGHTER's
          // ITS cluster count vs the V0's DECAY RADIUS (assoc.v0radius()), split into
          // "Toward"/"Transverse" via a delta-phi cut using checks.towardDeltaEtaRange /
          // checks.transwerseDeltaEtaRangeMin/Max -- but neither of those two things has
          // an equivalent here: Phi/K*0 candidates are reconstructed at the primary
          // vertex from track pairs, not displaced V0 decays, so there is no "decay
          // radius" analog on AssocPhis/AssocKstars, and the towardDeltaEtaRange /
          // transwerseDeltaEtaRangeMin/Max configurables that drive the region split
          // don't exist in this task's `checks` struct at all (only doOnTheFlyFlattening
          // does). Implementing this would mean inventing both a physical observable for
          // the third axis and the region-cut configurables, so it's removed rather than
          // guessed. If you want it back, tell me what the third axis should represent
          // for a non-displaced daughter pair (e.g. a DCA quantity, or drop to a TH2),
          // and I'll add the config fields and fill it properly.
        }

        // Core Correlation THnFs
        if (masterConfigurations.doFullCorrelationStudy) {
          if (doSE) {
            if (masterConfigurations.fillCorrelationHistWithMass && (i == IndexPhi || i == IndexKstar)) {
              histos.add(fmt::format("sameEvent/Signal/{}", Particlenames[i]).c_str(), "", kTHnF, {axisDeltaPhiNDim, (i == IndexPhi ? axesConfigurations.axisPhiMass : axesConfigurations.axisKstarMass), axisPtAssocNDim, axisPtTriggerNDim, axisVtxZNDim, axisMultNDim});
            } else {
              histos.add(fmt::format("sameEvent/Signal/{}", Particlenames[i]).c_str(), "", kTHnF, {axisDeltaPhiNDim, axisDeltaEtaNDim, axisPtAssocNDim, axisPtTriggerNDim, axisVtxZNDim, axisMultNDim});
            }
            if (doDeltaPhiStarCheck)
              histos.add(fmt::format("sameEvent/Signal/{}DeltaPhiStar", Particlenames[i]).c_str(), "", kTH3F, {{100, -0.3, 0.3}, {50, -0.05, 0.05}, {2, -1, 1}});
            if ((i == IndexPhi || i == IndexKstar) && !masterConfigurations.fillCorrelationHistWithMass) {
              needsSignalBgCloneSE = true;
            }
          }

          if (doME) {
            if (masterConfigurations.fillCorrelationHistWithMass && (i == IndexPhi || i == IndexKstar)) {
              histos.add(fmt::format("mixedEvent/Signal/{}", Particlenames[i]).c_str(), "", kTHnF, {axisDeltaPhiNDim, (i == IndexPhi ? axesConfigurations.axisPhiMass : axesConfigurations.axisKstarMass), axisPtAssocNDim, axisPtTriggerNDim, axisVtxZNDim, axisMultNDim});
            } else {
              histos.add(fmt::format("mixedEvent/Signal/{}", Particlenames[i]).c_str(), "", kTHnF, {axisDeltaPhiNDim, axisDeltaEtaNDim, axisPtAssocNDim, axisPtTriggerNDim, axisVtxZNDim, axisMultNDim});
            }
            if (doDeltaPhiStarCheck)
              histos.add(fmt::format("mixedEvent/Signal/{}DeltaPhiStar", Particlenames[i]).c_str(), "", kTH3F, {{100, -0.3, 0.3}, {50, -0.05, 0.05}, {2, -1, 1}});
            if ((i == IndexPhi || i == IndexKstar) && !masterConfigurations.fillCorrelationHistWithMass) {
              needsSignalBgCloneME = true;
            }
          }
        }
      }
    }
    // Clone Signal/ into LeftBg/ and RightBg/ once, after all species for this
    // loop have booked their Signal/ histograms (Phi and/or K*0 use sideband
    // subtraction; avoids cloning the same prefix twice, which the histogram
    // registry does not allow).
    if (needsSignalBgCloneSE) {
      histos.addClone("sameEvent/Signal/", "sameEvent/LeftBg/");
      histos.addClone("sameEvent/Signal/", "sameEvent/RightBg/");
    }
    if (needsSignalBgCloneME) {
      histos.addClone("mixedEvent/Signal/", "mixedEvent/LeftBg/");
      histos.addClone("mixedEvent/Signal/", "mixedEvent/RightBg/");
    }

    // ========================================================================
    // TRIGGER AND ASSOC SPECIFIC QA
    // ========================================================================
    if (TESTBIT(doCorrelation, IndexPhi) && doprocessSameEventHPhis && masterConfigurations.doFullCorrelationStudy) {
      histos.add("sameEvent/TriggerParticlesPhi", "TriggersPhi", kTH2F, {axesConfigurations.axisPtQA, axesConfigurations.axisMult});
      histos.add("hNumberOfRejectedPairsPhi", "hNumberOfRejectedPairsPhi", kTH1F, {{1, 0, 1}});
    }

    if (TESTBIT(doCorrelation, IndexKstar) && doprocessSameEventHKstars && masterConfigurations.doFullCorrelationStudy) {
      histos.add("sameEvent/TriggerParticlesKstar0", "TriggersKstar0", kTH2F, {axesConfigurations.axisPtQA, axesConfigurations.axisMult});
      histos.add("hNumberOfRejectedPairsKstar0", "hNumberOfRejectedPairsKstar0", kTH1F, {{1, 0, 1}});
    }

    if (TESTBIT(doCorrelation, IndexPion) && doprocessSameEventHPions) {
      if (masterConfigurations.doFullCorrelationStudy) {
        histos.add("sameEvent/TriggerParticlesPion", "TriggersPion", kTH2F, {axesConfigurations.axisPtQA, axesConfigurations.axisMult});
      }
      histos.add("hNumberOfRejectedPairsPion", "hNumberOfRejectedPairsPion", kTH1F, {{1, 0, 1}});
      histos.add("hPionEtaVsPtAllSelected", "", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisMult});
      histos.add("hPionEtaVsPt", "", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisMult});
      histos.add("hPositivePionEtaVsPt", "", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisMult});
      histos.add("hNegativePionEtaVsPt", "", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisMult});
    }

    if (TESTBIT(doCorrelation, 3) && doprocessSameEventHHadrons) {
      if (masterConfigurations.doFullCorrelationStudy) {
        histos.add("sameEvent/TriggerParticlesHadron", "TriggersHadron", kTH2F, {axesConfigurations.axisPtQA, axesConfigurations.axisMult});
      }
      histos.add("hNumberOfRejectedPairsHadron", "hNumberOfRejectedPairsHadron", kTH1F, {{1, 0, 1}});
      histos.add("hDCAzTriggerHadron", "hDCAzTriggerHadron", kTH2F, {{200, -0.5, 0.5}, axesConfigurations.axisPtQA});
      histos.add("hDCAxyTriggerHadron", "hDCAxyTriggerHadron", kTH2F, {{200, -0.5, 0.5}, axesConfigurations.axisPtQA});
      histos.add("hDCAzAssociatedHadron", "hDCAzAssociatedHadron", kTH2F, {{200, -0.5, 0.5}, axesConfigurations.axisPtQA});
      histos.add("hDCAxyAssociatedHadron", "hDCAxyAssociatedHadron", kTH2F, {{200, -0.5, 0.5}, axesConfigurations.axisPtQA});
      histos.add("hAsssocTrackEtaVsPtVsPhi", "", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisPhi});
      histos.add("hAssocPrimaryEtaVsPt", "", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisMult});
      histos.add("hAssocHadronsAllSelectedEtaVsPt", "", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisMult});
      histos.add("hAssocPtResolution", ";p_{T}^{reconstructed} (GeV/c); p_{T}^{generated} (GeV/c)", kTH2F, {axesConfigurations.axisPtQA, axesConfigurations.axisPtQA});
    }

    LOGF(info, "Init THnFs done");

    // ========================================================================
    // MC GENERATED, CLOSURE TEST, AND PREDICTION
    // ========================================================================
    if (doprocessMCGenerated) {
      histos.add("hGeneratedQAPtTrigger", "hGeneratedQAPtTrigger", kTH2F, {axesConfigurations.axisPtQA, {5, -0.5f, 4.5f}});
      histos.add("hGeneratedQAPtAssociatedPhi", "hGeneratedQAPtAssociatedPhi", kTH2F, {axesConfigurations.axisPtQA, {5, -0.5f, 4.5f}});
      histos.add("hClosureTestEventCounter", "hClosureTestEventCounter", kTH1F, {{10, 0, 10}});

      if (!efficiencyFlags.applyEffAsFunctionOfMultAndPhi) {
        histos.add("Generated/hTrigger", "", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisMult});
      } else {
        histos.add("Generated/hTrigger", "", kTHnF, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisPhi, axesConfigurations.axisMult});
      }
      histos.add("Generated/hPositiveTrigger", "", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisMult});
      histos.add("Generated/hNegativeTrigger", "", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisMult});

      for (int i = 0; i < AssocParticleTypes; i++) {
        if (TESTBIT(doCorrelation, i)) { // Fixed missing check here
          if (!efficiencyFlags.applyEffAsFunctionOfMultAndPhi) {
            histos.add(fmt::format("Generated/h{}", Particlenames[i]).c_str(), "", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisMult});
          } else {
            histos.add(fmt::format("Generated/h{}", Particlenames[i]).c_str(), "", kTHnF, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisPhi, axesConfigurations.axisMult});
          }
          if (i == IndexPion) {
            histos.add(fmt::format("Generated/hPositive{}", Particlenames[i]).c_str(), "", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisMult});
            histos.add(fmt::format("Generated/hNegative{}", Particlenames[i]).c_str(), "", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisMult});
          }
        }
      }
      histos.addClone("Generated/", "GeneratedWithPV/");

      for (int i = 0; i < AssocParticleTypesNoHadron; i++) {
        if (TESTBIT(doCorrelation, i)) {
          histos.add(fmt::format("GeneratedWithPV/h{}_MidYVsMult", Particlenames[i]).c_str(), "", kTH2F, {axesConfigurations.axisPtQA, axesConfigurations.axisMult});
          histos.add(fmt::format("GeneratedWithPV/h{}_MidYVsMult_TwoPVsOrMore", Particlenames[i]).c_str(), "", kTH2F, {axesConfigurations.axisPtQA, axesConfigurations.axisMult});
        }
      }
    }

    if (doprocessClosureTest) {
      if (!doprocessMCGenerated)
        histos.add("hClosureTestEventCounter", "hClosureTestEventCounter", kTH1F, {{10, 0, 10}});
      histos.add("hClosureQAPtTrigger", "hClosureQAPtTrigger", kTH2F, {axesConfigurations.axisPtQA, {5, -0.5f, 4.5f}});
      histos.add("hClosureQAPtAssociatedPhi", "hClosureQAPtAssociatedPhi", kTH2F, {axesConfigurations.axisPtQA, {5, -0.5f, 4.5f}});

      for (int i = 0; i < AssocParticleTypes; i++) {
        if (TESTBIT(doCorrelation, i)) {
          histos.add(fmt::format("ClosureTest/sameEvent/{}", Particlenames[i]).c_str(), "", kTHnF, {axisDeltaPhiNDim, axisDeltaEtaNDim, axisPtAssocNDim, axisPtTriggerNDim, axisVtxZNDim, axisMultNDim});
          histos.add(fmt::format("ClosureTest/h{}", Particlenames[i]).c_str(), "", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisPhi});
        }
      }
      histos.add("ClosureTest/hTrigger", "Trigger Tracks", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisMult});
    }

    if (doprocessPrediction) {
      mCounter.mPdgDatabase = pdgDB.service;
      mCounter.mSelectPrimaries = doAssocPhysicalPrimary.value;
      histos.add("Prediction/hEventSelection", "hEventSelection", kTH1F, {{3, 0, 3}});
      TString eventSelLabel[] = {"Read", "INELgt0", "|Z|<10"};
      for (int i = 1; i <= histos.get<TH1>(HIST("Prediction/hEventSelection"))->GetNbinsX(); i++) {
        histos.get<TH1>(HIST("Prediction/hEventSelection"))->GetXaxis()->SetBinLabel(i, eventSelLabel[i - 1]);
      }
      if (masterConfigurations.useCentralityinPrediction) {
        if (masterConfigurations.doSeparateFT0Prediction) {
          histos.add("Prediction/hTriggerFT0A", "Trigger Tracks FT0A", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisMult});
          histos.add("Prediction/hTriggerFT0C", "Trigger Tracks FT0C", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisMult});
        }
        histos.add("Prediction/hTrigger", "Trigger Tracks", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisMult});
      } else {
        if (masterConfigurations.doSeparateFT0Prediction) {
          histos.add("Prediction/hTriggerFT0A", "Trigger Tracks FT0A", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisMultiplicity});
          histos.add("Prediction/hFT0AvsNchEta08", "Nch in 0.8 vs FT0A multiplicity", kTH2F, {axesConfigurations.axisMultiplicity, axesConfigurations.axisMidrapidityMultiplicity});
          histos.add("Prediction/hFT0AvsNchEta05", "Nch in 0.5 vs FT0A multiplicity", kTH2F, {axesConfigurations.axisMultiplicity, axesConfigurations.axisMidrapidityMultiplicity});
          histos.add("Prediction/hTriggerFT0C", "Trigger Tracks FT0C", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisMultiplicity});
          histos.add("Prediction/hFT0CvsNchEta08", "Nch in 0.8 vs FT0C multiplicity", kTH2F, {axesConfigurations.axisMultiplicity, axesConfigurations.axisMidrapidityMultiplicity});
          histos.add("Prediction/hFT0CvsNchEta05", "Nch in 0.5 vs FT0C multiplicity", kTH2F, {axesConfigurations.axisMultiplicity, axesConfigurations.axisMidrapidityMultiplicity});
        }
        histos.add("Prediction/hTrigger", "Trigger Tracks", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisMultiplicity});
        histos.add("Prediction/hFT0MvsNchEta08", "Nch in 0.8 vs FT0M multiplicity", kTH2F, {axesConfigurations.axisMultiplicity, axesConfigurations.axisMidrapidityMultiplicity});
        histos.add("Prediction/hFT0MvsNchEta05", "Nch in 0.5 vs FT0M multiplicity", kTH2F, {axesConfigurations.axisMultiplicity, axesConfigurations.axisMidrapidityMultiplicity});
      }
      for (int i = 0; i < AssocParticleTypes; i++) {
        if (TESTBIT(doCorrelation, i)) {
          histos.add(fmt::format("Prediction/h{}", Particlenames[i]).c_str(), "", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisPhi});
          if (masterConfigurations.useCentralityinPrediction) {
            histos.add(fmt::format("Prediction/sameEvent/{}", Particlenames[i]).c_str(), "", kTHnF, {axisDeltaPhiNDim, axisDeltaEtaNDim, axisPtAssocNDim, axisPtTriggerNDim, axisVtxZNDim, axisMultNDim});
          } else {
            histos.add(fmt::format("Prediction/sameEvent/{}", Particlenames[i]).c_str(), "", kTHnF, {axisDeltaPhiNDim, axisDeltaEtaNDim, axisPtAssocNDim, axisPtTriggerNDim, axisVtxZNDim, axisMultiplicityNDim});
          }
        }
      }
      if (masterConfigurations.doSeparateFT0Prediction) {
        histos.addClone("Prediction/sameEvent/", "Prediction/sameEventFT0A/");
        histos.addClone("Prediction/sameEvent/", "Prediction/sameEventFT0C/");
      }
    }

    // visual inspection of sizes
    histos.print();

    // initialize CCDB *only* if efficiency correction requested
    // skip if not requested, saves a bit of time
    if (efficiencyFlags.applyEfficiencyCorrection) {
      ccdb->setURL(ccdburl);
      ccdb->setCaching(true);
      ccdb->setLocalObjectValidityChecking();
      ccdb->setFatalWhenNull(false);
    }
  }

  // this function allows for all event selections to be done in a modular way
  template <typename TCollision>
  bool isisCollisionSelect(TCollision const& collision)
  {
    // ________________________________________________
    // Perform basic event selection
    if (!collision.sel8()) {
      return false;
    }
    if (std::abs(collision.posZ()) > masterConfigurations.zVertexCut) {
      return false;
    }
    if (collision.centFT0M() > axisRanges[5][1] || collision.centFT0M() < axisRanges[5][0]) {
      return false;
    }
    if (!collision.isInelGt0() && masterConfigurations.selectINELgtZERO) {
      return false;
    }
    if (!collision.selection_bit(aod::evsel::kIsGoodITSLayersAll) && masterConfigurations.requireAllGoodITSLayers) {
      return false;
    }
    if (zorroMask.value != "") {
      auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
      initZorro(bc);
      bool zorroSelected = zorro.isSelected(collision.template bc_as<aod::BCsWithTimestamps>().globalBC()); /// Just let Zorro do the accounting
      if (!zorroSelected) {
        return false;
      }
    }
    return true;
  }

  // event selections in Pb-Pb
  template <typename TCollision>
  bool isisCollisionSelectPbPb(TCollision collision, bool fillHists)
  {
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 0.5 /* all collisions */);

    // Perform basic event selection
    if (!collision.sel8()) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 1.5 /* collisions  after sel8*/);

    if (!collision.selection_bit(aod::evsel::kIsTriggerTVX) && masterConfigurations.requireGoodTriggerTVX) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 2.5 /* FT0 vertex (acceptable FT0C-FT0A time difference) collisions */);

    if (std::abs(collision.posZ()) > masterConfigurations.zVertexCut) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 3.5 /* collisions  after sel pvz sel*/);

    if (!collision.selection_bit(aod::evsel::kIsGoodITSLayersAll) && masterConfigurations.requireAllGoodITSLayers) {
      // cut time intervals with dead ITS staves
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 4.5 /* collisions  after cut time intervals with dead ITS staves*/);

    if (!collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV) && masterConfigurations.requireGoodZvtxFT0vsPV) {
      // removes collisions with large differences between z of PV by tracks and z of PV from FT0 A-C time difference
      // use this cut at low multiplicities with caution
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 5.5 /* removes collisions with large differences between z of PV by tracks and z of PV from FT0 A-C time difference*/);

    auto occupancy = collision.trackOccupancyInTimeRange();
    if (occupancy < cfgCutOccupancyLow || occupancy > cfgCutOccupancyHigh)
      return false;
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 6.5 /* Below min occupancy and Above max occupancy*/);

    /*
    if (collision.alias_bit(kTVXinTRD)) {
      // TRD triggered
      return false;
    }
    */

    if (!collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
      // reject collisions close to Time Frame borders
      // O2-4623
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 7.5 /* reject collisions close to Time Frame borders*/);

    if (!collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
      // reject events affected by the ITS ROF border
      // O2-4309
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 8.5 /* reject events affected by the ITS ROF border*/);

    if (!collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      // rejects collisions which are associated with the same "found-by-T0" bunch crossing
      // https://indico.cern.ch/event/1396220/#1-event-selection-with-its-rof
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 9.5 /* rejects collisions which are associated with the same "found-by-T0" bunch crossing*/);
    return true;
  }

  template <class TTrack>
  bool isValidTrigger(TTrack track, bool isLeading)
  {
    if (track.tpcNClsCrossedRows() < trackSelection.minTPCNCrossedRowsTrigger) {
      return false; // crossed rows
    }
    if (!track.hasITS() && trackSelection.triggerRequireITS) {
      return false; // skip, doesn't have ITS signal (skips lots of TPC-only!)
    }
    if (track.tpcNClsShared() > trackSelection.triggerMaxTPCSharedClusters) {
      return false; // skip, has shared clusters
    }
    if (!(TESTBIT(track.itsClusterMap(), 0)) && trackSelection.triggerRequireL0) {
      return false; // skip, doesn't have cluster in ITS L0
    }
    // systematic variations: trigger DCAxy
    if (std::abs(track.dcaXY()) > trackSelection.dcaXYconstant + trackSelection.dcaXYpTdep * std::abs(track.signed1Pt())) {
      return false;
    }
    // systematic variations: trigger DCAz
    if (trackSelection.requireDCAzCut && std::abs(track.dcaZ()) > trackSelection.dcaZconstant + trackSelection.dcaZpTdep * std::abs(track.signed1Pt())) {
      return false;
    }
    if (track.pt() > axisRanges[3][1] || track.pt() < axisRanges[3][0]) {
      return false;
    }
    if (triggerParticleCharge > 0 && track.sign() < 0) {
      return false;
    }
    if (triggerParticleCharge < 0 && track.sign() > 0) {
      return false;
    }
    if (useTheLeadingParticleAsTrigger && !isLeading) {
      return false;
    }
    return true;
  }
  template <class TTrack>
  bool isValidAssocHadron(TTrack track)
  {
    if (track.tpcNClsCrossedRows() < trackSelection.minTPCNCrossedRowsAssociated) {
      return false; // crossed rows
    }
    if (!track.hasITS() && trackSelection.assocRequireITS) {
      return false; // skip, doesn't have ITS signal (skips lots of TPC-only!)
    }
    if (track.tpcNClsShared() > trackSelection.assocMaxTPCSharedClusters) {
      return false; // skip, has shared clusters
    }
    if (!(TESTBIT(track.itsClusterMap(), 0)) && trackSelection.assocRequireL0) {
      return false; // skip, doesn't have cluster in ITS L0
    }
    // systematic variations: trigger DCAxy
    if (std::abs(track.dcaXY()) > trackSelection.dcaXYconstantAssoc + trackSelection.dcaXYpTdepAssoc * std::abs(track.signed1Pt())) {
      return false;
    }
    // systematic variations: trigger DCAz
    if (trackSelection.requireDCAzCut && std::abs(track.dcaZ()) > trackSelection.dcaZconstantAssoc + trackSelection.dcaZpTdepAssoc * std::abs(track.signed1Pt())) {
      return false;
    }
    if (track.pt() > axisRanges[2][1] || track.pt() < axisRanges[2][0]) {
      return false;
    }
    return true;
  }

  double calculateAverageDeltaPhiStar(const double* trigg, const double* assoc, double B)
  {
    double dPhiStar = 0;
    double dPhiStarMean = 0;

    double dPhi = assoc[0] - trigg[0];
    double phaseProton = (-0.3 * B * assoc[2]) / (2 * assoc[1]);
    double phaseTrack = (-0.3 * B * trigg[2]) / (2 * trigg[1]);

    for (double r = MinRadiusTPC; r <= MaxRadiusTPC; r += 0.05) {
      dPhiStar = dPhi + std::asin(phaseProton * r) - std::asin(phaseTrack * r);
      dPhiStarMean += (dPhiStar / 34);
    }

    return dPhiStarMean;
  }
  void fillTriggerHistogram(std::shared_ptr<TH2> hist, double pt, double mult, float eff, float effUncert, float purity, float purityErr)
  {
    int binx = hist->GetXaxis()->FindBin(pt);
    int biny = hist->GetYaxis()->FindBin(mult);
    float previousContent = hist->GetBinContent(binx, biny);
    float previousUncert = hist->GetBinError(binx, biny);
    float newContent = previousContent + purity / eff;
    float newUncert = std::sqrt(previousUncert * previousUncert + std::pow(purity / eff, 2) + std::pow(purityErr / eff, 2) + std::pow(effUncert, 2) / std::pow(eff, 4));
    hist->SetBinContent(binx, biny, newContent);
    hist->SetBinError(binx, biny, newUncert);
  }
  void fillCorrelationHistogram(std::shared_ptr<THn> hist, double binFillThn[], float etaWeight, float efficiency, float totalEffUncert, float purity, float totalPurityUncert)
  {
    float previousContent, previousError2, currentContent, currentError2;
    int bin = hist->GetBin(binFillThn);
    previousContent = hist->GetBinContent(bin);
    previousError2 = hist->GetBinError2(bin);
    currentContent = previousContent + etaWeight * purity / (efficiency);
    currentError2 = previousError2 + std::pow(etaWeight * purity / (efficiency), 2) + std::pow(etaWeight * totalPurityUncert / (efficiency), 2) + std::pow(totalEffUncert * purity * etaWeight, 2) / std::pow(efficiency, 4);
    hist->SetBinContent(bin, currentContent);
    hist->SetBinError2(bin, currentError2);
  }

  template <typename TTriggers, typename THadrons>
  void fillCorrelationsHadron(TTriggers const& triggers, THadrons const& assocs, bool mixing, float pvz, float mult, double bField)
  {

    for (auto const& triggerTrack : triggers) {
      if (masterConfigurations.doTriggPhysicalPrimary && !triggerTrack.mcPhysicalPrimary())
        continue;
      auto trigg = triggerTrack.template track_as<TracksComplete>();
      if (!isValidTrigger(trigg, triggerTrack.isLeading()))
        continue;

      float efficiencyTrigger = 1.0f;
      float efficiencyTriggerError = 0.0f;
      float purityTrigger = 1.0f;
      float purityTriggerError = 0.0f;
      if (efficiencyFlags.applyEfficiencyForTrigger) {
        if (efficiencyFlags.applyEffAsFunctionOfMult)
          efficiencyTrigger = hEfficiencyTriggerMult->Interpolate(trigg.pt(), trigg.eta(), mult);
        else
          efficiencyTrigger = hEfficiencyTrigger->Interpolate(trigg.pt(), trigg.eta());
        if (efficiencyFlags.applyPurityTrigger) {
          if (efficiencyFlags.applyEffAsFunctionOfMult)
            purityTrigger = hPurityHadronMult->Interpolate(trigg.pt(), mult);
          else
            purityTrigger = hPurityHadron->Interpolate(trigg.pt());
        }
        if (efficiencyFlags.applyEfficiencyPropagation) {
          if (efficiencyFlags.applyEffAsFunctionOfMult)
            efficiencyTriggerError = hEfficiencyUncertaintyTriggerMult->Interpolate(trigg.pt(), trigg.eta(), mult);
          else
            efficiencyTriggerError = hEfficiencyUncertaintyTrigger->Interpolate(trigg.pt(), trigg.eta());
          if (efficiencyFlags.applyPurityTrigger) {
            if (efficiencyFlags.applyEffAsFunctionOfMult)
              purityTriggerError = hPurityUncertaintyHadronMult->Interpolate(trigg.pt(), mult);
            else
              purityTriggerError = hPurityUncertaintyHadron->Interpolate(trigg.pt());
          }
        }
        if (efficiencyTrigger == 0) { // check for zero efficiency, do not apply if the case
          efficiencyTrigger = 1;
          efficiencyTriggerError = 0;
        }
      }
      if (!mixing) {
        if constexpr (requires { triggerTrack.extra(); })
          fillTriggerHistogram(histos.get<TH2>(HIST("sameEvent/TriggerParticlesPion")), trigg.pt(), mult, efficiencyTrigger, efficiencyTriggerError, purityTrigger, purityTriggerError);
        else
          fillTriggerHistogram(histos.get<TH2>(HIST("sameEvent/TriggerParticlesHadron")), trigg.pt(), mult, efficiencyTrigger, efficiencyTriggerError, purityTrigger, purityTriggerError);
      }
      double triggSign = trigg.sign();
      const double triggForDeltaPhiStar[] = {trigg.phi(), trigg.pt(), triggSign};
      for (auto const& assocTrack : assocs) {
        auto assoc = assocTrack.template track_as<TracksComplete>();

        //---] removing autocorrelations [---
        if (doAutocorrelationRejection) {
          if (trigg.globalIndex() == assoc.globalIndex()) {
            if constexpr (requires { assocTrack.nSigmaTPCPi(); })
              histos.fill(HIST("hNumberOfRejectedPairsPion"), 0.5);
            else
              histos.fill(HIST("hNumberOfRejectedPairsHadron"), 0.5);
            continue;
          }
        }
        //---] track quality check [---
        if (!isValidAssocHadron(assoc))
          continue;
        if (doAssocPhysicalPrimary && !assocTrack.mcPhysicalPrimary()) {
          continue;
        }
        float deltaphi = computeDeltaPhi(trigg.phi(), assoc.phi());
        float deltaeta = trigg.eta() - assoc.eta();
        if (masterConfigurations.doMirroringInDelataEta) {
          deltaeta = std::abs(deltaeta);
        }
        float ptassoc = assoc.pt();
        float pttrigger = trigg.pt();

        double assocSign = assoc.sign();
        const double assocForDeltaPhiStar[] = {assoc.phi(), assoc.pt(), assocSign};

        float etaWeight = 1.;
        if (checks.doOnTheFlyFlattening) {
          float preWeight = 1 - std::abs(deltaeta) / 1.6;
          etaWeight = preWeight != 0 ? 1.0f / preWeight : 1.0f;
        }

        // skip if basic ranges not met
        if (deltaphi < axisRanges[0][0] || deltaphi > axisRanges[0][1])
          continue;
        if (deltaeta < axisRanges[1][0] || deltaeta > axisRanges[1][1])
          continue;
        if (ptassoc < axisRanges[2][0] || ptassoc > axisRanges[2][1])
          continue;

        float efficiency = 1;
        float purity = 1.0f;
        float purityUncertainty = 0.0f;
        float totalEffUncert = 0.0;
        float efficiencyUncertainty = 0.0f;
        float totalPurityUncert = 0.0;
        if (efficiencyFlags.applyEfficiencyCorrection) {
          if constexpr (requires { assocTrack.nSigmaTPCPi(); }) {
            efficiency = hEfficiencyPion->Interpolate(ptassoc, assoc.eta());
            if (efficiencyFlags.applyEfficiencyPropagation)
              efficiencyUncertainty = hEfficiencyUncertaintyPion->Interpolate(ptassoc, assoc.eta());
          } else {
            if (efficiencyFlags.applyEffAsFunctionOfMult)
              efficiency = hEfficiencyHadronMult->Interpolate(ptassoc, assoc.eta(), mult);
            else
              efficiency = hEfficiencyHadron->Interpolate(ptassoc, assoc.eta());
            if (efficiencyFlags.applyPurityHadron) {
              if (efficiencyFlags.applyEffAsFunctionOfMult)
                purity = hPurityHadronMult->Interpolate(ptassoc, mult);
              else
                purity = hPurityHadron->Interpolate(ptassoc);
            }
            if (efficiencyFlags.applyEfficiencyPropagation) {
              if (efficiencyFlags.applyEffAsFunctionOfMult)
                efficiencyUncertainty = hEfficiencyUncertaintyHadronMult->Interpolate(ptassoc, assoc.eta(), mult);
              else
                efficiencyUncertainty = hEfficiencyUncertaintyHadron->Interpolate(ptassoc, assoc.eta());
              if (efficiencyFlags.applyPurityHadron) {
                if (efficiencyFlags.applyEffAsFunctionOfMult)
                  purityUncertainty = hPurityUncertaintyHadronMult->Interpolate(ptassoc, mult);
                else
                  purityUncertainty = hPurityUncertaintyHadron->Interpolate(ptassoc);
              }
            }
          }
        }
        if (efficiency == 0) { // check for zero efficiency, do not apply if the case
          efficiency = 1;
          efficiencyUncertainty = 0.0;
        }
        if (efficiencyFlags.applyEfficiencyPropagation) {
          totalEffUncert = std::sqrt(std::pow(efficiencyTrigger * efficiencyUncertainty, 2) + std::pow(efficiencyTriggerError * efficiency, 2));
          totalPurityUncert = std::sqrt(std::pow(purityTrigger * purityUncertainty, 2) + std::pow(purity * purityTriggerError, 2));
        }

        double binFillThn[6] = {deltaphi, deltaeta, ptassoc, pttrigger, pvz, mult};
        double deltaPhiStar = calculateAverageDeltaPhiStar(triggForDeltaPhiStar, assocForDeltaPhiStar, bField);
        if (!mixing) {
          if constexpr (requires { assocTrack.nSigmaTPCPi(); }) {
            fillCorrelationHistogram(histos.get<THn>(HIST("sameEvent/Signal/Pion")), binFillThn, etaWeight, efficiency * efficiencyTrigger, totalEffUncert, purity * purityTrigger, totalPurityUncert);
            if (triggSign == assocSign && doDeltaPhiStarCheck) {
              histos.fill(HIST("sameEvent/Signal/Pion") + HIST("DeltaPhiStar"), deltaPhiStar, trigg.eta() - assoc.eta(), 0.5);
            } else if (doDeltaPhiStarCheck) {
              histos.fill(HIST("sameEvent/Signal/Pion") + HIST("DeltaPhiStar"), deltaPhiStar, trigg.eta() - assoc.eta(), -0.5);
            }
          } else {
            if (triggSign == assocSign && doDeltaPhiStarCheck) {
              histos.fill(HIST("sameEvent/Signal/Hadron") + HIST("DeltaPhiStar"), deltaPhiStar, trigg.eta() - assoc.eta(), 0.5);
            } else if (doDeltaPhiStarCheck) {
              histos.fill(HIST("sameEvent/Signal/Hadron") + HIST("DeltaPhiStar"), deltaPhiStar, trigg.eta() - assoc.eta(), -0.5);
            }
            fillCorrelationHistogram(histos.get<THn>(HIST("sameEvent/Signal/Hadron")), binFillThn, etaWeight, efficiency * efficiencyTrigger, totalEffUncert, purity * purityTrigger, totalPurityUncert);
          }
        } else {
          if constexpr (requires { assocTrack.nSigmaTPCPi(); }) {
            fillCorrelationHistogram(histos.get<THn>(HIST("mixedEvent/Signal/Pion")), binFillThn, 1, efficiency * efficiencyTrigger, totalEffUncert, purity * purityTrigger, totalPurityUncert);
          } else {
            fillCorrelationHistogram(histos.get<THn>(HIST("mixedEvent/Signal/Hadron")), binFillThn, 1, efficiency * efficiencyTrigger, totalEffUncert, purity * purityTrigger, totalPurityUncert);
          }
        }
      }
    }
  }

  // ---- OPTIMIZED (2026-09-24): daughter-track resolution, efficiency lookup and
  // mass-window classification hoisted out of the trigger loop; see comment block
  // originally delivered in fillCorrelationsPhi_Kstar_optimized.cxx for rationale ----
  void fillCorrelationsPhi(aod::TriggerTracks const& triggers, aod::AssocPhis const& assocs, bool mixing, bool mixingInBf, float pvz, float mult, double bField)
  {
    ValidCollision currentCollision;
    int binMult = 0;
    int nBinsMult = 0;
    int nBinsVtxZ = 0;
    int binVtxZ = 0;
    currentCollision.pvz = pvz;
    currentCollision.mult = mult;
    if (mixingInBf) {
      nBinsMult = histos.get<TH1>(HIST("axes/hMultAxis"))->GetNbinsX();
      binMult = histos.get<TH1>(HIST("axes/hMultAxis"))->GetXaxis()->FindBin(mult) - 1;
      nBinsVtxZ = histos.get<TH1>(HIST("axes/hVertexZAxis"))->GetNbinsX();
      binVtxZ = histos.get<TH1>(HIST("axes/hVertexZAxis"))->GetXaxis()->FindBin(pvz) - 1;
    }

    // ---------------------------------------------------------------------
    // Pre-resolve everything that does NOT depend on the trigger, once per
    // assoc candidate instead of once per (trigger, assoc) pair.
    // ---------------------------------------------------------------------
    struct CachedPhi {
      float phi, eta, pt, mass;
      float phiKplus, etaKplus, ptKplus, signKplus;
      float phiKminus, etaKminus, ptKminus, signKminus;
      int32_t globalIndexPos, globalIndexNeg;
      float efficiency, efficiencyError;
      bool isLeftBg, isSignal, isRightBg;
      bool passesMcSelection; // (!doMCassociation || mcTruePhi()) && (!doAssocPhysicalPrimary || mcPhysicalPrimary())
    };
    std::vector<CachedPhi> cache;
    cache.reserve(assocs.size());
    for (auto const& assocCandidate : assocs) {
      auto postrack = assocCandidate.posTrackDaughter_as<TracksComplete>();
      auto negtrack = assocCandidate.negTrackDaughter_as<TracksComplete>();

      CachedPhi c;
      c.phi = assocCandidate.phi();
      c.eta = assocCandidate.eta();
      c.pt = assocCandidate.pt();
      c.mass = assocCandidate.mass();

      c.phiKplus = postrack.phi();
      c.etaKplus = postrack.eta();
      c.ptKplus = postrack.pt();
      c.signKplus = postrack.sign();
      c.phiKminus = negtrack.phi();
      c.etaKminus = negtrack.eta();
      c.ptKminus = negtrack.pt();
      c.signKminus = negtrack.sign();
      c.globalIndexPos = postrack.globalIndex();
      c.globalIndexNeg = negtrack.globalIndex();

      c.efficiency = 1.0f;
      c.efficiencyError = 0.0f;
      if (efficiencyFlags.applyEfficiencyCorrection) {
        if (!efficiencyFlags.applyEffAsFunctionOfMultAndPhi) {
          c.efficiency = hEfficiencyPhi->Interpolate(c.pt, c.eta);
          if (efficiencyFlags.applyEfficiencyPropagation)
            c.efficiencyError = hEfficiencyUncertaintyPhi->Interpolate(c.pt, c.eta);
        } else {
          double bin[4] = {c.pt, c.eta, c.phi, mult};
          c.efficiency = hEfficiencyPhiMultVsPhi->GetBinContent(hEfficiencyPhiMultVsPhi->GetBin(bin));
          if (efficiencyFlags.applyEfficiencyPropagation)
            c.efficiencyError = hEfficiencyPhiMultVsPhi->GetBinError(hEfficiencyPhiMultVsPhi->GetBin(bin));
        }
      }

      // Peak-relative classification using Phi's own mass window.
      float delta = c.mass - massWindowConfigurationsPhi.peakMass;
      float sig = massWindowConfigurationsPhi.sigma;
      c.isLeftBg = (-massWindowConfigurationsPhi.maxBgNSigma * sig < delta && delta < -massWindowConfigurationsPhi.minBgNSigma * sig);
      c.isSignal = (-massWindowConfigurationsPhi.maxPeakNSigma * sig < delta && delta < massWindowConfigurationsPhi.maxPeakNSigma * sig);
      c.isRightBg = (massWindowConfigurationsPhi.minBgNSigma * sig < delta && delta < massWindowConfigurationsPhi.maxBgNSigma * sig);

      c.passesMcSelection = (!masterConfigurations.doMCassociation || assocCandidate.mcTruePhi()) &&
                            (!doAssocPhysicalPrimary || assocCandidate.mcPhysicalPrimary());

      // Candidate-level kinematic/spectrum QA, ported from the equivalent
      // block in hStrangeCorrelation.cxx (h3d<Species>Spectrum[Y],
      // h<Species>EtaVsPtVsPhi[Bg]). These are registered (histos.add(...)
      // above) but only filled here. Filled once per candidate per
      // SAME EVENT only (never during mixing -- these are candidate QA, not
      // correlation-pair QA, so filling them again per mixed pair would
      // overcount each candidate by the mixing-pool depth).
      // "SpectrumY" mirrors the V0 analog's rapidity cut (std::abs(rapidity) <
      // ySel); AssocPhis has no exposed rapidity() accessor in this file, so
      // pseudorapidity (c.eta) is used as a stand-in -- swap in a real
      // rapidity() call here if/when the table exposes one.
      if (!mixing) {
        float nSigmaVal = delta / sig;
        histos.fill(HIST("h3dPhiSpectrum"), c.pt, mult, nSigmaVal);
        if (std::abs(c.eta) < ySel) {
          histos.fill(HIST("h3dPhiSpectrumY"), c.pt, mult, nSigmaVal);
        }
        // matches the two registration branches above (TH3F vs THnF-with-mult,
        // selected by efficiencyFlags.applyEffAsFunctionOfMultAndPhi)
        if (!efficiencyFlags.applyEffAsFunctionOfMultAndPhi) {
          if (c.isSignal) {
            histos.fill(HIST("hPhiEtaVsPtVsPhi"), c.pt, c.eta, c.phi);
          }
          if (c.isLeftBg || c.isRightBg) {
            histos.fill(HIST("hPhiEtaVsPtVsPhiBg"), c.pt, c.eta, c.phi);
          }
        } else {
          if (c.isSignal) {
            histos.fill(HIST("hPhiEtaVsPtVsPhiVsCent"), c.pt, c.eta, c.phi, mult);
          }
          if (c.isLeftBg || c.isRightBg) {
            histos.fill(HIST("hPhiEtaVsPtVsPhiVsCentBg"), c.pt, c.eta, c.phi, mult);
          }
        }
      }

      cache.push_back(c);
    }
    // ---------------------------------------------------------------------

    bool firstLoop = false;
    for (auto const& triggerTrack : triggers) {
      if (masterConfigurations.doTriggPhysicalPrimary && !triggerTrack.mcPhysicalPrimary())
        continue;
      auto trigg = triggerTrack.track_as<TracksComplete>();
      if (!isValidTrigger(trigg, triggerTrack.isLeading()))
        continue;

      float efficiencyTrigg = 1.0f;
      float efficiencyTriggError = 0.0f;
      float purityTrigg = 1.0f;
      float purityTriggErr = 0.0;
      double bintrig[4] = {trigg.pt(), trigg.eta(), trigg.phi(), mult};

      if (efficiencyFlags.applyEfficiencyForTrigger) {
        if (!efficiencyFlags.applyEffAsFunctionOfMultAndPhi) {
          efficiencyTrigg = hEfficiencyTrigger->Interpolate(trigg.pt(), trigg.eta());
          if (efficiencyFlags.applyPurityTrigger)
            purityTrigg = hPurityHadron->Interpolate(trigg.pt());
          if (efficiencyFlags.applyEfficiencyPropagation) {
            efficiencyTriggError = hEfficiencyUncertaintyTrigger->Interpolate(trigg.pt(), trigg.eta());
            if (efficiencyFlags.applyPurityTrigger)
              purityTriggErr = hPurityHadron->Interpolate(trigg.pt());
          }
        } else {
          efficiencyTrigg = hEfficiencyTriggerMultVsPhi->GetBinContent(hEfficiencyTriggerMultVsPhi->GetBin(bintrig));
          if (efficiencyFlags.applyEfficiencyPropagation) {
            efficiencyTriggError = hEfficiencyTriggerMultVsPhi->GetBinError(hEfficiencyTriggerMultVsPhi->GetBin(bintrig));
          }
        }
        if (efficiencyTrigg == 0) {
          efficiencyTrigg = 1;
          efficiencyTriggError = 0;
        }
      }

      if (!mixing) {
        fillTriggerHistogram(histos.get<TH2>(HIST("sameEvent/TriggerParticlesPhi")), trigg.pt(), mult, efficiencyTrigg, efficiencyTriggError, purityTrigg, purityTriggErr);
      }

      double triggSign = trigg.sign();
      const double triggForDeltaPhiStar[] = {trigg.phi(), trigg.pt(), triggSign};

      if (mixingInBf) {
        currentCollision.addValidParticle(trigg.eta(), trigg.phi(), trigg.pt(), -1, efficiencyTrigg, efficiencyTriggError, -1);
        if (firstLoop)
          continue;
      }

      for (auto const& c : cache) {
        firstLoop = true;

        //---] removing autocorrelations [---
        if (doAutocorrelationRejection) {
          if (trigg.globalIndex() == c.globalIndexPos || trigg.globalIndex() == c.globalIndexNeg) {
            histos.fill(HIST("hNumberOfRejectedPairsPhi"), 0.5);
            continue;
          }
        }

        float deltaphi = computeDeltaPhi(trigg.phi(), c.phi);
        float deltaeta = trigg.eta() - c.eta;
        if (masterConfigurations.doMirroringInDelataEta) {
          deltaeta = std::abs(deltaeta);
        }
        float ptassoc = c.pt;
        float pttrigger = trigg.pt();
        float massassoc = c.mass;

        // skip if basic ranges not met
        if (deltaphi < axisRanges[0][0] || deltaphi > axisRanges[0][1])
          continue;
        if (deltaeta < axisRanges[1][0] || deltaeta > axisRanges[1][1])
          continue;
        if (ptassoc < axisRanges[2][0] || ptassoc > axisRanges[2][1])
          continue;

        float etaWeight = 1;
        if (checks.doOnTheFlyFlattening) {
          float preWeight = 1 - std::abs(deltaeta) / 1.6;
          etaWeight = preWeight != 0 ? 1.0f / preWeight : 1.0f;
        }

        float totalEffUncert = 0.0f;
        if (efficiencyFlags.applyEfficiencyPropagation) {
          totalEffUncert = std::sqrt(std::pow(efficiencyTrigg * c.efficiencyError, 2) + std::pow(efficiencyTriggError * c.efficiency, 2));
        }

        double binFillThn[6] = {deltaphi, deltaeta, ptassoc, pttrigger, pvz, mult};

        const double assocForDeltaPhiStarPlus[] = {c.phiKplus, c.ptKplus, c.signKplus};
        const double assocForDeltaPhiStarMinus[] = {c.phiKminus, c.ptKminus, c.signKminus};

        if (!c.passesMcSelection)
          continue;

        if (!mixing && c.isLeftBg && !masterConfigurations.fillCorrelationHistWithMass) {
          fillCorrelationHistogram(histos.get<THn>(HIST("sameEvent/LeftBg/Phi")), binFillThn, etaWeight, c.efficiency * efficiencyTrigg, totalEffUncert, purityTrigg, purityTriggErr);
          if (doDeltaPhiStarCheck) {
            double deltaPhiStarPlus = calculateAverageDeltaPhiStar(triggForDeltaPhiStar, assocForDeltaPhiStarPlus, bField);
            double deltaPhiStarMinus = calculateAverageDeltaPhiStar(triggForDeltaPhiStar, assocForDeltaPhiStarMinus, bField);
            if (triggSign > 0) {
              histos.fill(HIST("sameEvent/Signal/PhiDeltaPhiStar"), deltaPhiStarPlus, trigg.eta() - c.etaKplus, +0.5);
              histos.fill(HIST("sameEvent/Signal/PhiDeltaPhiStar"), deltaPhiStarMinus, trigg.eta() - c.etaKminus, -0.5);
            } else {
              histos.fill(HIST("sameEvent/Signal/PhiDeltaPhiStar"), deltaPhiStarPlus, trigg.eta() - c.etaKplus, -0.5);
              histos.fill(HIST("sameEvent/Signal/PhiDeltaPhiStar"), deltaPhiStarMinus, trigg.eta() - c.etaKminus, +0.5);
            }
          }
        }
        if (!mixing && (masterConfigurations.fillCorrelationHistWithMass || c.isSignal)) {
          if (masterConfigurations.fillCorrelationHistWithMass) {
            binFillThn[1] = massassoc;
          }
          fillCorrelationHistogram(histos.get<THn>(HIST("sameEvent/Signal/Phi")), binFillThn, etaWeight, c.efficiency * efficiencyTrigg, totalEffUncert, purityTrigg, purityTriggErr);
          if (doDeltaPhiStarCheck) {
            double deltaPhiStarPlus = calculateAverageDeltaPhiStar(triggForDeltaPhiStar, assocForDeltaPhiStarPlus, bField);
            double deltaPhiStarMinus = calculateAverageDeltaPhiStar(triggForDeltaPhiStar, assocForDeltaPhiStarMinus, bField);
            if (triggSign > 0) {
              histos.fill(HIST("sameEvent/Signal/PhiDeltaPhiStar"), deltaPhiStarPlus, trigg.eta() - c.etaKplus, +0.5);
              histos.fill(HIST("sameEvent/Signal/PhiDeltaPhiStar"), deltaPhiStarMinus, trigg.eta() - c.etaKminus, -0.5);
            } else {
              histos.fill(HIST("sameEvent/Signal/PhiDeltaPhiStar"), deltaPhiStarPlus, trigg.eta() - c.etaKplus, -0.5);
              histos.fill(HIST("sameEvent/Signal/PhiDeltaPhiStar"), deltaPhiStarMinus, trigg.eta() - c.etaKminus, +0.5);
            }
          }
        }
        if (!mixing && c.isRightBg && !masterConfigurations.fillCorrelationHistWithMass) {
          fillCorrelationHistogram(histos.get<THn>(HIST("sameEvent/RightBg/Phi")), binFillThn, etaWeight, c.efficiency * efficiencyTrigg, totalEffUncert, purityTrigg, purityTriggErr);
          if (doDeltaPhiStarCheck) {
            double deltaPhiStarPlus = calculateAverageDeltaPhiStar(triggForDeltaPhiStar, assocForDeltaPhiStarPlus, bField);
            double deltaPhiStarMinus = calculateAverageDeltaPhiStar(triggForDeltaPhiStar, assocForDeltaPhiStarMinus, bField);
            if (triggSign > 0) {
              histos.fill(HIST("sameEvent/Signal/PhiDeltaPhiStar"), deltaPhiStarPlus, trigg.eta() - c.etaKplus, +0.5);
              histos.fill(HIST("sameEvent/Signal/PhiDeltaPhiStar"), deltaPhiStarMinus, trigg.eta() - c.etaKminus, -0.5);
            } else {
              histos.fill(HIST("sameEvent/Signal/PhiDeltaPhiStar"), deltaPhiStarPlus, trigg.eta() - c.etaKplus, -0.5);
              histos.fill(HIST("sameEvent/Signal/PhiDeltaPhiStar"), deltaPhiStarMinus, trigg.eta() - c.etaKminus, +0.5);
            }
          }
        }

        // Mixed Event Logic buffered in ValidCollision
        if (mixing && c.isLeftBg) {
          if (mixingInBf) {
            currentCollision.addValidParticle(c.eta, c.phi, c.pt, 0, c.efficiency, c.efficiencyError, 0 /* 0 = Phi type */);
          } else {
            fillCorrelationHistogram(histos.get<THn>(HIST("mixedEvent/LeftBg/Phi")), binFillThn, 1, c.efficiency * efficiencyTrigg, totalEffUncert, purityTrigg, purityTriggErr);
          }
        }
        if (mixing && c.isSignal) {
          if (mixingInBf) {
            currentCollision.addValidParticle(c.eta, c.phi, c.pt, 1, c.efficiency, c.efficiencyError, 0);
          } else {
            fillCorrelationHistogram(histos.get<THn>(HIST("mixedEvent/Signal/Phi")), binFillThn, 1, c.efficiency * efficiencyTrigg, totalEffUncert, purityTrigg, purityTriggErr);
          }
        }
        if (mixing && c.isRightBg) {
          if (mixingInBf) {
            currentCollision.addValidParticle(c.eta, c.phi, c.pt, 2, c.efficiency, c.efficiencyError, 0);
          } else {
            fillCorrelationHistogram(histos.get<THn>(HIST("mixedEvent/RightBg/Phi")), binFillThn, 1, c.efficiency * efficiencyTrigg, totalEffUncert, purityTrigg, purityTriggErr);
          }
        }
      }
    }

    // Process Buffer for Mixed Events -- unchanged from original.
    if (!mixingInBf || binVtxZ < 0 || binVtxZ > nBinsVtxZ - 1 || binMult < 0 || binMult > nBinsMult - 1)
      return;

    int binnumb = binMult * nBinsVtxZ + binVtxZ;
    int hastirgorassoc = masterConfigurations.collisionHasTriggOrAssoc;
    if ((hastirgorassoc == 1 && currentCollision.trigParticles.empty()) ||
        (hastirgorassoc == 2 && currentCollision.assocParticles.empty()) ||
        (hastirgorassoc == 3 && currentCollision.trigParticles.empty() && currentCollision.assocParticles.empty()) ||
        (hastirgorassoc == 4 && (currentCollision.trigParticles.empty() || currentCollision.assocParticles.empty())))
      return;

    for (const auto& collision : validCollisions[binnumb]) {
      BinningTypePP colBinning{{axesConfigurations.axisVtxZ, axesConfigurations.axisMult}, true};
      histos.fill(HIST("MixingQA/hMECollisionBins"), colBinning.getBin({collision.pvz, collision.mult}));

      for (const auto& trigger : collision.trigParticles) {
        for (const auto& assoc : currentCollision.assocParticles) {
          float deltaeta = trigger.eta - assoc.eta;
          float deltaphi = computeDeltaPhi(trigger.phi, assoc.phi);
          float efficiencyTrigg = trigger.efficiency;
          float efficiencyAssoc = assoc.efficiency;
          float efficiencyTriggError = trigger.efficiencyError;
          float efficiencyAssocError = assoc.efficiencyError;
          float totalEffUncert = 0.0;
          float ptassoc = assoc.pt;
          float pttrigger = trigger.pt;

          if (efficiencyFlags.applyEfficiencyPropagation) {
            totalEffUncert = std::sqrt(std::pow(efficiencyTrigg * efficiencyAssocError, 2) + std::pow(efficiencyTriggError * efficiencyAssoc, 2));
          }
          if (masterConfigurations.doMirroringInDelataEta) {
            deltaeta = std::abs(deltaeta);
          }
          double binFillThn[6] = {deltaphi, deltaeta, ptassoc, pttrigger, pvz, mult};

          if (assoc.region == 0) {
            fillCorrelationHistogram(histos.get<THn>(HIST("mixedEvent/LeftBg/Phi")), binFillThn, 1, efficiencyTrigg * efficiencyAssoc, totalEffUncert, 1., 0.);
          }
          if (assoc.region == 1) {
            fillCorrelationHistogram(histos.get<THn>(HIST("mixedEvent/Signal/Phi")), binFillThn, 1, efficiencyTrigg * efficiencyAssoc, totalEffUncert, 1., 0.);
          }
          if (assoc.region == 2) {
            fillCorrelationHistogram(histos.get<THn>(HIST("mixedEvent/RightBg/Phi")), binFillThn, 1, efficiencyTrigg * efficiencyAssoc, totalEffUncert, 1., 0.);
          }
        }
      }
    }

    if (validCollisions[binnumb].size() >= static_cast<size_t>(masterConfigurations.mixingParameter)) {
      validCollisions[binnumb].erase(validCollisions[binnumb].begin());
    }
    if (!currentCollision.trigParticles.empty())
      validCollisions[binnumb].push_back(currentCollision);
  }

  // ============================================================================
  // fillCorrelationsKstar -- identical restructuring, mirrors fillCorrelationsPhi
  // exactly (K*0 has no separate LeftBg/Signal/RightBg physics difference from
  // Phi in this file -- same three-region sideband pattern, same daughter-pair
  // caching opportunity). mcTruePhi() -> mcTrueKstar() is the only selection
  // difference, and hEfficiencyPhi* -> hEfficiencyKstar* for the efficiency
  // histograms.
  // ============================================================================

  void fillCorrelationsKstar(aod::TriggerTracks const& triggers, aod::AssocKstars const& assocs, bool mixing, bool mixingInBf, float pvz, float mult, double bField)
  {
    ValidCollision currentCollision;
    int binMult = 0;
    int nBinsMult = 0;
    int nBinsVtxZ = 0;
    int binVtxZ = 0;
    currentCollision.pvz = pvz;
    currentCollision.mult = mult;
    if (mixingInBf) {
      nBinsMult = histos.get<TH1>(HIST("axes/hMultAxis"))->GetNbinsX();
      binMult = histos.get<TH1>(HIST("axes/hMultAxis"))->GetXaxis()->FindBin(mult) - 1;
      nBinsVtxZ = histos.get<TH1>(HIST("axes/hVertexZAxis"))->GetNbinsX();
      binVtxZ = histos.get<TH1>(HIST("axes/hVertexZAxis"))->GetXaxis()->FindBin(pvz) - 1;
    }

    // ---------------------------------------------------------------------
    // Same pre-pass as fillCorrelationsPhi -- this matters most for K*0,
    // since assocKstars typically has far more candidates per event than
    // assocPhis, so the N_trig-fold redundant daughter-track resolution is
    // multiplied by a much bigger candidate count here than for Phi.
    // ---------------------------------------------------------------------
    struct CachedKstar {
      float phi, eta, pt, mass;
      float phiPos, etaPos, ptPos, signPos;
      float phiNeg, etaNeg, ptNeg, signNeg;
      int32_t globalIndexPos, globalIndexNeg;
      float efficiency, efficiencyError;
      bool isLeftBg, isSignal, isRightBg;
      bool passesMcSelection; // (!doMCassociation || mcTrueKstar()) && (!doAssocPhysicalPrimary || mcPhysicalPrimary())
    };
    std::vector<CachedKstar> cache;
    cache.reserve(assocs.size());
    for (auto const& assocCandidate : assocs) {
      auto postrack = assocCandidate.posTrackDaughter_as<TracksComplete>();
      auto negtrack = assocCandidate.negTrackDaughter_as<TracksComplete>();

      CachedKstar c;
      c.phi = assocCandidate.phi();
      c.eta = assocCandidate.eta();
      c.pt = assocCandidate.pt();
      c.mass = assocCandidate.mass();

      c.phiPos = postrack.phi();
      c.etaPos = postrack.eta();
      c.ptPos = postrack.pt();
      c.signPos = postrack.sign();
      c.phiNeg = negtrack.phi();
      c.etaNeg = negtrack.eta();
      c.ptNeg = negtrack.pt();
      c.signNeg = negtrack.sign();
      c.globalIndexPos = postrack.globalIndex();
      c.globalIndexNeg = negtrack.globalIndex();

      c.efficiency = 1.0f;
      c.efficiencyError = 0.0f;
      if (efficiencyFlags.applyEfficiencyCorrection) {
        if (!efficiencyFlags.applyEffAsFunctionOfMultAndPhi) {
          c.efficiency = hEfficiencyKstar->Interpolate(c.pt, c.eta);
          if (efficiencyFlags.applyEfficiencyPropagation)
            c.efficiencyError = hEfficiencyUncertaintyKstar->Interpolate(c.pt, c.eta);
        } else {
          double bin[4] = {c.pt, c.eta, c.phi, mult};
          c.efficiency = hEfficiencyKstarMultVsPhi->GetBinContent(hEfficiencyKstarMultVsPhi->GetBin(bin));
          if (efficiencyFlags.applyEfficiencyPropagation)
            c.efficiencyError = hEfficiencyKstarMultVsPhi->GetBinError(hEfficiencyKstarMultVsPhi->GetBin(bin));
        }
      }

      // Peak-relative classification using K*0's own mass window.
      float delta = c.mass - massWindowConfigurationsKstar.peakMass;
      float sig = massWindowConfigurationsKstar.sigma;
      c.isLeftBg = (-massWindowConfigurationsKstar.maxBgNSigma * sig < delta && delta < -massWindowConfigurationsKstar.minBgNSigma * sig);
      c.isSignal = (-massWindowConfigurationsKstar.maxPeakNSigma * sig < delta && delta < massWindowConfigurationsKstar.maxPeakNSigma * sig);
      c.isRightBg = (massWindowConfigurationsKstar.minBgNSigma * sig < delta && delta < massWindowConfigurationsKstar.maxBgNSigma * sig);

      c.passesMcSelection = (!masterConfigurations.doMCassociation || assocCandidate.mcTrueKstar()) &&
                            (!doAssocPhysicalPrimary || assocCandidate.mcPhysicalPrimary());

      // Same candidate-level Spectrum/EtaVsPtVsPhi QA as fillCorrelationsPhi
      // above -- see the comment there for rationale.
      if (!mixing) {
        float nSigmaVal = delta / sig;
        histos.fill(HIST("h3dKstar0Spectrum"), c.pt, mult, nSigmaVal);
        if (std::abs(c.eta) < ySel) {
          histos.fill(HIST("h3dKstar0SpectrumY"), c.pt, mult, nSigmaVal);
        }
        if (!efficiencyFlags.applyEffAsFunctionOfMultAndPhi) {
          if (c.isSignal) {
            histos.fill(HIST("hKstar0EtaVsPtVsPhi"), c.pt, c.eta, c.phi);
          }
          if (c.isLeftBg || c.isRightBg) {
            histos.fill(HIST("hKstar0EtaVsPtVsPhiBg"), c.pt, c.eta, c.phi);
          }
        } else {
          if (c.isSignal) {
            histos.fill(HIST("hKstar0EtaVsPtVsPhiVsCent"), c.pt, c.eta, c.phi, mult);
          }
          if (c.isLeftBg || c.isRightBg) {
            histos.fill(HIST("hKstar0EtaVsPtVsPhiVsCentBg"), c.pt, c.eta, c.phi, mult);
          }
        }
      }

      cache.push_back(c);
    }
    // ---------------------------------------------------------------------

    bool firstLoop = false;
    for (auto const& triggerTrack : triggers) {
      if (masterConfigurations.doTriggPhysicalPrimary && !triggerTrack.mcPhysicalPrimary())
        continue;
      auto trigg = triggerTrack.track_as<TracksComplete>();
      if (!isValidTrigger(trigg, triggerTrack.isLeading()))
        continue;

      float efficiencyTrigg = 1.0f;
      float efficiencyTriggError = 0.0f;
      float purityTrigg = 1.0f;
      float purityTriggErr = 0.0;
      double bintrig[4] = {trigg.pt(), trigg.eta(), trigg.phi(), mult};

      if (efficiencyFlags.applyEfficiencyForTrigger) {
        if (!efficiencyFlags.applyEffAsFunctionOfMultAndPhi) {
          efficiencyTrigg = hEfficiencyTrigger->Interpolate(trigg.pt(), trigg.eta());
          if (efficiencyFlags.applyPurityTrigger)
            purityTrigg = hPurityHadron->Interpolate(trigg.pt());
          if (efficiencyFlags.applyEfficiencyPropagation) {
            efficiencyTriggError = hEfficiencyUncertaintyTrigger->Interpolate(trigg.pt(), trigg.eta());
            if (efficiencyFlags.applyPurityTrigger)
              purityTriggErr = hPurityHadron->Interpolate(trigg.pt());
          }
        } else {
          efficiencyTrigg = hEfficiencyTriggerMultVsPhi->GetBinContent(hEfficiencyTriggerMultVsPhi->GetBin(bintrig));
          if (efficiencyFlags.applyEfficiencyPropagation) {
            efficiencyTriggError = hEfficiencyTriggerMultVsPhi->GetBinError(hEfficiencyTriggerMultVsPhi->GetBin(bintrig));
          }
        }
        if (efficiencyTrigg == 0) {
          efficiencyTrigg = 1;
          efficiencyTriggError = 0;
        }
      }

      if (!mixing) {
        fillTriggerHistogram(histos.get<TH2>(HIST("sameEvent/TriggerParticlesKstar0")), trigg.pt(), mult, efficiencyTrigg, efficiencyTriggError, purityTrigg, purityTriggErr);
      }

      double triggSign = trigg.sign();
      const double triggForDeltaPhiStar[] = {trigg.phi(), trigg.pt(), triggSign};

      if (mixingInBf) {
        currentCollision.addValidParticle(trigg.eta(), trigg.phi(), trigg.pt(), -1, efficiencyTrigg, efficiencyTriggError, -1);
        if (firstLoop)
          continue;
      }

      for (auto const& c : cache) {
        firstLoop = true;

        //---] removing autocorrelations [---
        if (doAutocorrelationRejection) {
          if (trigg.globalIndex() == c.globalIndexPos || trigg.globalIndex() == c.globalIndexNeg) {
            histos.fill(HIST("hNumberOfRejectedPairsKstar0"), 0.5);
            continue;
          }
        }

        float deltaphi = computeDeltaPhi(trigg.phi(), c.phi);
        float deltaeta = trigg.eta() - c.eta;
        if (masterConfigurations.doMirroringInDelataEta) {
          deltaeta = std::abs(deltaeta);
        }
        float ptassoc = c.pt;
        float pttrigger = trigg.pt();
        float massassoc = c.mass;

        // skip if basic ranges not met
        if (deltaphi < axisRanges[0][0] || deltaphi > axisRanges[0][1])
          continue;
        if (deltaeta < axisRanges[1][0] || deltaeta > axisRanges[1][1])
          continue;
        if (ptassoc < axisRanges[2][0] || ptassoc > axisRanges[2][1])
          continue;

        float etaWeight = 1;
        if (checks.doOnTheFlyFlattening) {
          float preWeight = 1 - std::abs(deltaeta) / 1.6;
          etaWeight = preWeight != 0 ? 1.0f / preWeight : 1.0f;
        }

        float totalEffUncert = 0.0f;
        if (efficiencyFlags.applyEfficiencyPropagation) {
          totalEffUncert = std::sqrt(std::pow(efficiencyTrigg * c.efficiencyError, 2) + std::pow(efficiencyTriggError * c.efficiency, 2));
        }

        double binFillThn[6] = {deltaphi, deltaeta, ptassoc, pttrigger, pvz, mult};

        const double assocForDeltaPhiStarPlus[] = {c.phiPos, c.ptPos, c.signPos};
        const double assocForDeltaPhiStarMinus[] = {c.phiNeg, c.ptNeg, c.signNeg};

        if (!c.passesMcSelection)
          continue;

        if (!mixing && c.isLeftBg && !masterConfigurations.fillCorrelationHistWithMass) {
          fillCorrelationHistogram(histos.get<THn>(HIST("sameEvent/LeftBg/Kstar0")), binFillThn, etaWeight, c.efficiency * efficiencyTrigg, totalEffUncert, purityTrigg, purityTriggErr);
          if (doDeltaPhiStarCheck) {
            double deltaPhiStarPlus = calculateAverageDeltaPhiStar(triggForDeltaPhiStar, assocForDeltaPhiStarPlus, bField);
            double deltaPhiStarMinus = calculateAverageDeltaPhiStar(triggForDeltaPhiStar, assocForDeltaPhiStarMinus, bField);
            if (triggSign > 0) {
              histos.fill(HIST("sameEvent/Signal/Kstar0DeltaPhiStar"), deltaPhiStarPlus, trigg.eta() - c.etaPos, +0.5);
              histos.fill(HIST("sameEvent/Signal/Kstar0DeltaPhiStar"), deltaPhiStarMinus, trigg.eta() - c.etaNeg, -0.5);
            } else {
              histos.fill(HIST("sameEvent/Signal/Kstar0DeltaPhiStar"), deltaPhiStarPlus, trigg.eta() - c.etaPos, -0.5);
              histos.fill(HIST("sameEvent/Signal/Kstar0DeltaPhiStar"), deltaPhiStarMinus, trigg.eta() - c.etaNeg, +0.5);
            }
          }
        }
        if (!mixing && (masterConfigurations.fillCorrelationHistWithMass || c.isSignal)) {
          if (masterConfigurations.fillCorrelationHistWithMass) {
            binFillThn[1] = massassoc;
          }
          fillCorrelationHistogram(histos.get<THn>(HIST("sameEvent/Signal/Kstar0")), binFillThn, etaWeight, c.efficiency * efficiencyTrigg, totalEffUncert, purityTrigg, purityTriggErr);
          if (doDeltaPhiStarCheck) {
            double deltaPhiStarPlus = calculateAverageDeltaPhiStar(triggForDeltaPhiStar, assocForDeltaPhiStarPlus, bField);
            double deltaPhiStarMinus = calculateAverageDeltaPhiStar(triggForDeltaPhiStar, assocForDeltaPhiStarMinus, bField);
            if (triggSign > 0) {
              histos.fill(HIST("sameEvent/Signal/Kstar0DeltaPhiStar"), deltaPhiStarPlus, trigg.eta() - c.etaPos, +0.5);
              histos.fill(HIST("sameEvent/Signal/Kstar0DeltaPhiStar"), deltaPhiStarMinus, trigg.eta() - c.etaNeg, -0.5);
            } else {
              histos.fill(HIST("sameEvent/Signal/Kstar0DeltaPhiStar"), deltaPhiStarPlus, trigg.eta() - c.etaPos, -0.5);
              histos.fill(HIST("sameEvent/Signal/Kstar0DeltaPhiStar"), deltaPhiStarMinus, trigg.eta() - c.etaNeg, +0.5);
            }
          }
        }
        if (!mixing && c.isRightBg && !masterConfigurations.fillCorrelationHistWithMass) {
          fillCorrelationHistogram(histos.get<THn>(HIST("sameEvent/RightBg/Kstar0")), binFillThn, etaWeight, c.efficiency * efficiencyTrigg, totalEffUncert, purityTrigg, purityTriggErr);
          if (doDeltaPhiStarCheck) {
            double deltaPhiStarPlus = calculateAverageDeltaPhiStar(triggForDeltaPhiStar, assocForDeltaPhiStarPlus, bField);
            double deltaPhiStarMinus = calculateAverageDeltaPhiStar(triggForDeltaPhiStar, assocForDeltaPhiStarMinus, bField);
            if (triggSign > 0) {
              histos.fill(HIST("sameEvent/Signal/Kstar0DeltaPhiStar"), deltaPhiStarPlus, trigg.eta() - c.etaPos, +0.5);
              histos.fill(HIST("sameEvent/Signal/Kstar0DeltaPhiStar"), deltaPhiStarMinus, trigg.eta() - c.etaNeg, -0.5);
            } else {
              histos.fill(HIST("sameEvent/Signal/Kstar0DeltaPhiStar"), deltaPhiStarPlus, trigg.eta() - c.etaPos, -0.5);
              histos.fill(HIST("sameEvent/Signal/Kstar0DeltaPhiStar"), deltaPhiStarMinus, trigg.eta() - c.etaNeg, +0.5);
            }
          }
        }

        // Mixed-event logic: buffered into ValidCollision when mixing via the
        // collision pool (processMixedEventHKstarsInBuffer), filled directly
        // otherwise (processMixedEventHKstars).
        if (mixing && c.isLeftBg) {
          if (mixingInBf) {
            currentCollision.addValidParticle(c.eta, c.phi, c.pt, 0, c.efficiency, c.efficiencyError, 1 /* 1 = Kstar0 type */);
          } else {
            fillCorrelationHistogram(histos.get<THn>(HIST("mixedEvent/LeftBg/Kstar0")), binFillThn, 1, c.efficiency * efficiencyTrigg, totalEffUncert, purityTrigg, purityTriggErr);
          }
        }
        if (mixing && c.isSignal) {
          if (mixingInBf) {
            currentCollision.addValidParticle(c.eta, c.phi, c.pt, 1, c.efficiency, c.efficiencyError, 1);
          } else {
            fillCorrelationHistogram(histos.get<THn>(HIST("mixedEvent/Signal/Kstar0")), binFillThn, 1, c.efficiency * efficiencyTrigg, totalEffUncert, purityTrigg, purityTriggErr);
          }
        }
        if (mixing && c.isRightBg) {
          if (mixingInBf) {
            currentCollision.addValidParticle(c.eta, c.phi, c.pt, 2, c.efficiency, c.efficiencyError, 1);
          } else {
            fillCorrelationHistogram(histos.get<THn>(HIST("mixedEvent/RightBg/Kstar0")), binFillThn, 1, c.efficiency * efficiencyTrigg, totalEffUncert, purityTrigg, purityTriggErr);
          }
        }
      }
    }

    // Process buffer for mixed events (mirrors fillCorrelationsPhi; reached
    // when called from processMixedEventHKstarsInBuffer with mixingInBf=true)
    if (!mixingInBf || binVtxZ < 0 || binVtxZ > nBinsVtxZ - 1 || binMult < 0 || binMult > nBinsMult - 1)
      return;

    int binnumb = binMult * nBinsVtxZ + binVtxZ;
    int hastirgorassoc = masterConfigurations.collisionHasTriggOrAssoc;
    if ((hastirgorassoc == 1 && currentCollision.trigParticles.empty()) ||
        (hastirgorassoc == 2 && currentCollision.assocParticles.empty()) ||
        (hastirgorassoc == 3 && currentCollision.trigParticles.empty() && currentCollision.assocParticles.empty()) ||
        (hastirgorassoc == 4 && (currentCollision.trigParticles.empty() || currentCollision.assocParticles.empty())))
      return;

    for (const auto& collision : validCollisions[binnumb]) {
      BinningTypePP colBinning{{axesConfigurations.axisVtxZ, axesConfigurations.axisMult}, true};
      histos.fill(HIST("MixingQA/hMECollisionBins"), colBinning.getBin({collision.pvz, collision.mult}));

      for (const auto& trigger : collision.trigParticles) {
        for (const auto& assoc : currentCollision.assocParticles) {
          float deltaeta = trigger.eta - assoc.eta;
          float deltaphi = computeDeltaPhi(trigger.phi, assoc.phi);
          float efficiencyTrigg = trigger.efficiency;
          float efficiencyAssoc = assoc.efficiency;
          float efficiencyTriggError = trigger.efficiencyError;
          float efficiencyAssocError = assoc.efficiencyError;
          float totalEffUncert = 0.0;
          float ptassoc = assoc.pt;
          float pttrigger = trigger.pt;

          if (efficiencyFlags.applyEfficiencyPropagation) {
            totalEffUncert = std::sqrt(std::pow(efficiencyTrigg * efficiencyAssocError, 2) + std::pow(efficiencyTriggError * efficiencyAssoc, 2));
          }
          if (masterConfigurations.doMirroringInDelataEta) {
            deltaeta = std::abs(deltaeta);
          }
          double binFillThn[6] = {deltaphi, deltaeta, ptassoc, pttrigger, pvz, mult};

          if (assoc.region == 0) {
            fillCorrelationHistogram(histos.get<THn>(HIST("mixedEvent/LeftBg/Kstar0")), binFillThn, 1, efficiencyTrigg * efficiencyAssoc, totalEffUncert, 1., 0.);
          }
          if (assoc.region == 1) {
            fillCorrelationHistogram(histos.get<THn>(HIST("mixedEvent/Signal/Kstar0")), binFillThn, 1, efficiencyTrigg * efficiencyAssoc, totalEffUncert, 1., 0.);
          }
          if (assoc.region == 2) {
            fillCorrelationHistogram(histos.get<THn>(HIST("mixedEvent/RightBg/Kstar0")), binFillThn, 1, efficiencyTrigg * efficiencyAssoc, totalEffUncert, 1., 0.);
          }
        }
      }
    }

    if (validCollisions[binnumb].size() >= static_cast<size_t>(masterConfigurations.mixingParameter)) {
      validCollisions[binnumb].erase(validCollisions[binnumb].begin());
    }
    if (!currentCollision.trigParticles.empty())
      validCollisions[binnumb].push_back(currentCollision);
  }

  double getMagneticField(uint64_t timestamp)
  {
    static parameters::GRPMagField* grpo = nullptr;
    if (grpo == nullptr) {
      grpo = ccdb->getForTimeStamp<parameters::GRPMagField>("GLO/Config/GRPMagField", timestamp);
      if (grpo == nullptr) {
        LOGF(fatal, "GRP object not found for timestamp %llu", timestamp);
        return 0;
      }
      LOGF(info, "Retrieved GRP for timestamp %llu with magnetic field of %d kG", timestamp, grpo->getNominalL3Field());
    }

    return 0.1 * (grpo->getNominalL3Field()); // 1 T = 10 kG
  }
  // if this process function is enabled, it will be such that only events with trigger particles within a given
  // trigger pt bin are taken for the entire processing. This allows for the calculation of e.g. efficiencies
  // within an event class that has a trigger (which may differ with respect to other cases, to be checked)

  // for map determining which trigger bins are present and which aren't
  std::vector<uint32_t> triggerPresenceMap;

  void processSelectEventWithTrigger(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::PVMults> const& collisions, aod::TriggerTracks const& triggerTracks, TracksComplete const&, aod::BCsWithTimestamps const&)
  {
    // setup
    triggerPresenceMap.clear();
    triggerPresenceMap.resize(collisions.size(), 0);

    for (auto const& collision : collisions) {
      // ________________________________________________
      // Perform basic event selection
      if (!isisCollisionSelect(collision)) {
        continue;
      }

      // do not forget to re-group ...
      auto slicedTriggerTracks = triggerTracks.sliceBy(collisionSliceTracks, collision.globalIndex());

      for (auto const& triggerTrack : slicedTriggerTracks) {
        auto track = triggerTrack.track_as<TracksComplete>();
        if (!isValidTrigger(track, triggerTrack.isLeading())) {
          continue;
        }
        auto binNumber = histos.get<TH1>(HIST("axes/hPtTriggerAxis"))->FindFixBin(track.pt()) - 1;
        SETBIT(triggerPresenceMap[collision.globalIndex()], binNumber);
      }
    }
  }

  void processSameEventHHadrons(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::PVMults>::iterator const& collision, aod::AssocHadrons const& assocHadrons, aod::TriggerTracks const& triggerTracks, TracksComplete const&, aod::BCsWithTimestamps const&)
  {
    LOGF(info, "SameEventHadron: collisions=%d triggers=%zu assocHadrons=%zu", collision.globalIndex(), triggerTracks.size(), assocHadrons.size());

    BinningTypePP colBinning{{axesConfigurations.axisVtxZ, axesConfigurations.axisMult}, true}; // true is for 'ignore overflows' (true by default). Underflows and overflows will have bin -1.

    // ________________________________________________
    // skip if desired trigger not found
    if (triggerPresenceMap.size() > 0 && !TESTBIT(triggerPresenceMap[collision.globalIndex()], triggerBinToSelect)) {
      return;
    }

    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    auto bField = getMagneticField(bc.timestamp());

    if (efficiencyFlags.applyEfficiencyCorrection) {
      initEfficiencyFromCCDB(bc);
    }

    // ________________________________________________
    // Perform basic event selection
    if (!isisCollisionSelect(collision)) {
      return;
    }
    // ________________________________________________
    if (!doprocessSameEventHPhis && !doprocessSameEventHKstars && !doprocessSameEventHPions && doMixingQAandEventQA) {
      histos.fill(HIST("MixingQA/hSECollisionBins"), colBinning.getBin({collision.posZ(), collision.centFT0M()}));
      histos.fill(HIST("EventQA/hMult"), collision.centFT0M());
      histos.fill(HIST("EventQA/hPvz"), collision.posZ());
    }

    // Do basic QA
    if (!doprocessSameEventHPhis && !doprocessSameEventHKstars && !doprocessSameEventHPions) {
      for (auto const& triggerTrack : triggerTracks) {
        auto track = triggerTrack.track_as<TracksComplete>();
        if (!isValidTrigger(track, triggerTrack.isLeading()))
          continue;
        histos.fill(HIST("hDCAzTriggerHadron"), track.dcaZ(), track.pt());
        histos.fill(HIST("hDCAxyTriggerHadron"), track.dcaXY(), track.pt());
        float efficiency = 1.0f;
        if (efficiencyFlags.applyEfficiencyCorrection) {
          efficiency = hEfficiencyTrigger->Interpolate(track.pt(), track.eta());
        }
        if (efficiency == 0) { // check for zero efficiency, do not apply if the case
          efficiency = 1;
        }
        float weight = efficiencyFlags.applyEfficiencyCorrection ? 1. / efficiency : 1.0f;
        histos.fill(HIST("hTriggerAllSelectedEtaVsPt"), track.pt(), track.eta(), collision.centFT0M());
        histos.fill(HIST("hTriggerPtResolution"), track.pt(), triggerTrack.mcOriginalPt());
        if (masterConfigurations.doTriggPhysicalPrimary && !triggerTrack.mcPhysicalPrimary())
          continue;
        histos.fill(HIST("hTriggerPrimaryEtaVsPt"), track.pt(), track.eta(), collision.centFT0M());
        histos.fill(HIST("hTrackEtaVsPtVsPhi"), track.pt(), track.eta(), track.phi(), weight);
      }
    }
    for (auto const& assocTrack : assocHadrons) {
      auto assoc = assocTrack.track_as<TracksComplete>();
      if (!isValidAssocHadron(assoc))
        continue;
      float efficiency = 1.0f;
      float purity = 1.0f;
      histos.fill(HIST("hDCAzAssociatedHadron"), assoc.dcaZ(), assoc.pt());
      histos.fill(HIST("hDCAxyAssociatedHadron"), assoc.dcaXY(), assoc.pt());
      if (efficiencyFlags.applyEfficiencyCorrection) {
        efficiency = hEfficiencyHadron->Interpolate(assoc.pt(), assoc.eta());
        if (efficiencyFlags.applyPurityHadron)
          purity = hPurityHadron->Interpolate(assoc.pt());
      }
      if (efficiency == 0) { // check for zero efficiency, do not apply if the case
        efficiency = 1;
      }
      float weight = efficiencyFlags.applyEfficiencyCorrection ? purity / efficiency : 1.0f;
      histos.fill(HIST("hAssocHadronsAllSelectedEtaVsPt"), assoc.pt(), assoc.eta(), collision.centFT0M(), weight);
      histos.fill(HIST("hAssocPtResolution"), assoc.pt(), assocTrack.mcOriginalPt());
      if (doAssocPhysicalPrimary && !assocTrack.mcPhysicalPrimary())
        continue;
      histos.fill(HIST("hAssocPrimaryEtaVsPt"), assoc.pt(), assoc.eta(), collision.centFT0M());
      histos.fill(HIST("hAsssocTrackEtaVsPtVsPhi"), assoc.pt(), assoc.eta(), assoc.phi(), weight);
    }

    // ________________________________________________
    // Do hadron - hadron correlations
    if (masterConfigurations.doFullCorrelationStudy)
      fillCorrelationsHadron(triggerTracks, assocHadrons, false, collision.posZ(), collision.centFT0M(), bField);
  }

  void processSameEventHPions(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::PVMults>::iterator const& collision, soa::Join<aod::AssocHadrons, aod::AssocPID> const& associatedPions, soa::Join<aod::TriggerTracks, aod::TriggerTrackExtras> const& triggerTracks, TracksComplete const&, aod::BCsWithTimestamps const&)
  {
    LOGF(info, "SameEventPion: collisions=%d triggers=%zu assocPions=%zu", collision.globalIndex(), triggerTracks.size(), associatedPions.size());

    BinningTypePP colBinning{{axesConfigurations.axisVtxZ, axesConfigurations.axisMult}, true};
    // ________________________________________________
    // skip if desired trigger not found
    if (triggerPresenceMap.size() > 0 && !TESTBIT(triggerPresenceMap[collision.globalIndex()], triggerBinToSelect)) {
      return;
    }
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    auto bField = getMagneticField(bc.timestamp());

    if (efficiencyFlags.applyEfficiencyCorrection) {
      initEfficiencyFromCCDB(bc);
    }

    // ________________________________________________
    // Perform basic event selection
    if (!isisCollisionSelect(collision)) {
      return;
    }
    // ________________________________________________
    if (!doprocessSameEventHPhis && !doprocessSameEventHKstars && doMixingQAandEventQA) {
      histos.fill(HIST("MixingQA/hSECollisionBins"), colBinning.getBin({collision.posZ(), collision.centFT0M()}));
      histos.fill(HIST("EventQA/hMult"), collision.centFT0M());
      histos.fill(HIST("EventQA/hPvz"), collision.posZ());
    }
    // Do basic QA
    for (auto const& pion : associatedPions) {
      auto pionTrack = pion.track_as<TracksComplete>();
      if (!isValidAssocHadron(pionTrack))
        continue;

      histos.fill(HIST("hPionEtaVsPtAllSelected"), pionTrack.pt(), pionTrack.eta(), collision.centFT0M());
      if (doAssocPhysicalPrimary && !pion.mcPhysicalPrimary())
        continue;
      if (masterConfigurations.doMCassociation && std::abs(pion.pdgCode()) != PdgCodes[IndexPion])
        continue;
      histos.fill(HIST("hPionEtaVsPt"), pionTrack.pt(), pionTrack.eta(), collision.centFT0M());
      if (pionTrack.sign() > 0)
        histos.fill(HIST("hPositivePionEtaVsPt"), pionTrack.pt(), pionTrack.eta(), collision.centFT0M());
      else
        histos.fill(HIST("hNegativePionEtaVsPt"), pionTrack.pt(), pionTrack.eta(), collision.centFT0M());
    }
    if (!doprocessSameEventHPhis && !doprocessSameEventHKstars) {
      for (auto const& triggerTrack : triggerTracks) {
        auto track = triggerTrack.track_as<TracksComplete>();
        if (!isValidTrigger(track, triggerTrack.isLeading()))
          continue;
        histos.fill(HIST("hTriggerAllSelectedEtaVsPt"), track.pt(), track.eta(), collision.centFT0M());
        histos.fill(HIST("hTriggerPtResolution"), track.pt(), triggerTrack.mcOriginalPt());
        if (masterConfigurations.doTriggPhysicalPrimary && !triggerTrack.mcPhysicalPrimary())
          continue;
        histos.fill(HIST("hTriggerPrimaryEtaVsPt"), track.pt(), track.eta(), collision.centFT0M());
        histos.fill(HIST("hTrackEtaVsPtVsPhi"), track.pt(), track.eta(), track.phi());
      }
    }

    // ________________________________________________
    // Do hadron - Pion correlations
    if (masterConfigurations.doFullCorrelationStudy)
      fillCorrelationsHadron(triggerTracks, associatedPions, false, collision.posZ(), collision.centFT0M(), bField);
  }

  void processSameEventHPhis(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::PVMults>::iterator const& collision, aod::AssocPhis const& associatedPhis, aod::TriggerTracks const& triggerTracks, TracksComplete const&, aod::BCsWithTimestamps const&)
  {
    LOGF(info, "SameEventPhi: collisions=%d triggers=%zu assocPhi=%zu", collision.globalIndex(), triggerTracks.size(), associatedPhis.size());

    BinningTypePP colBinning{{axesConfigurations.axisVtxZ, axesConfigurations.axisMult}, true};

    // Skip if desired trigger not found
    if (triggerPresenceMap.size() > 0 && !TESTBIT(triggerPresenceMap[collision.globalIndex()], triggerBinToSelect)) {
      return;
    }

    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    auto bField = getMagneticField(bc.timestamp());

    if (efficiencyFlags.applyEfficiencyCorrection) {
      initEfficiencyFromCCDB(bc);
    }

    if (!isisCollisionSelect(collision)) {
      return;
    }

    if (doMixingQAandEventQA) {
      histos.fill(HIST("MixingQA/hSECollisionBins"), colBinning.getBin({collision.posZ(), collision.centFT0M()}));
      histos.fill(HIST("EventQA/hMult"), collision.centFT0M());
      histos.fill(HIST("EventQA/hPvz"), collision.posZ());
    }

    // -----------------------------
    // Trigger QA
    // -----------------------------
    for (auto const& triggerTrack : triggerTracks) {
      auto track = triggerTrack.track_as<TracksComplete>();

      if (!isValidTrigger(track, triggerTrack.isLeading())) {
        continue;
      }

      histos.fill(HIST("hTriggerAllSelectedEtaVsPt"), track.pt(), track.eta(), collision.centFT0M());

      histos.fill(HIST("hTriggerPtResolution"), track.pt(), triggerTrack.mcOriginalPt());

      if (masterConfigurations.doTriggPhysicalPrimary && !triggerTrack.mcPhysicalPrimary()) {
        continue;
      }

      histos.fill(HIST("hTriggerPrimaryEtaVsPt"), track.pt(), track.eta(), collision.centFT0M());

      histos.fill(HIST("hTrackEtaVsPtVsPhi"), track.pt(), track.eta(), track.phi());
    }

    // -----------------------------
    // h-phi correlation
    // -----------------------------
    fillCorrelationsPhi(triggerTracks, associatedPhis, false, false, collision.posZ(), collision.centFT0M(), bField);
  }

  void processSameEventHKstars(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::PVMults>::iterator const& collision, aod::AssocKstars const& associatedKstars, aod::TriggerTracks const& triggerTracks, TracksComplete const&, aod::BCsWithTimestamps const&)
  {
    LOGF(info, "SameEventKstar: collisions=%d triggers=%zu assocKstar=%zu", collision.globalIndex(), triggerTracks.size(), associatedKstars.size());

    BinningTypePP colBinning{{axesConfigurations.axisVtxZ, axesConfigurations.axisMult}, true};

    // Skip if desired trigger not found
    if (triggerPresenceMap.size() > 0 && !TESTBIT(triggerPresenceMap[collision.globalIndex()], triggerBinToSelect)) {
      return;
    }

    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    auto bField = getMagneticField(bc.timestamp());

    if (efficiencyFlags.applyEfficiencyCorrection) {
      initEfficiencyFromCCDB(bc);
    }

    if (!isisCollisionSelect(collision)) {
      return;
    }

    if (doMixingQAandEventQA) {
      histos.fill(HIST("MixingQA/hSECollisionBins"), colBinning.getBin({collision.posZ(), collision.centFT0M()}));
      histos.fill(HIST("EventQA/hMult"), collision.centFT0M());
      histos.fill(HIST("EventQA/hPvz"), collision.posZ());
    }

    // -----------------------------
    // Trigger QA
    // -----------------------------
    for (auto const& triggerTrack : triggerTracks) {
      auto track = triggerTrack.track_as<TracksComplete>();

      if (!isValidTrigger(track, triggerTrack.isLeading())) {
        continue;
      }

      histos.fill(HIST("hTriggerAllSelectedEtaVsPt"), track.pt(), track.eta(), collision.centFT0M());

      histos.fill(HIST("hTriggerPtResolution"), track.pt(), triggerTrack.mcOriginalPt());

      if (masterConfigurations.doTriggPhysicalPrimary && !triggerTrack.mcPhysicalPrimary()) {
        continue;
      }

      histos.fill(HIST("hTriggerPrimaryEtaVsPt"), track.pt(), track.eta(), collision.centFT0M());

      histos.fill(HIST("hTrackEtaVsPtVsPhi"), track.pt(), track.eta(), track.phi());
    }

    // -----------------------------
    // h-K*0 correlation
    // -----------------------------
    fillCorrelationsKstar(triggerTracks, associatedKstars, false, false, collision.posZ(), collision.centFT0M(), bField);
  }

  void
    processMixedEventHHadrons(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::PVMults> const& collisions, aod::AssocHadrons const& assocHadrons, aod::TriggerTracks const& triggerTracks, TracksComplete const&, aod::BCsWithTimestamps const&)
  {
    BinningTypePP colBinning{{axesConfigurations.axisVtxZ, axesConfigurations.axisMult}, true};
    for (auto const& [collision1, collision2] : soa::selfCombinations(colBinning, masterConfigurations.mixingParameter, -1, collisions, collisions)) {
      auto bc = collision1.bc_as<aod::BCsWithTimestamps>();
      auto bField = getMagneticField(bc.timestamp());
      // ________________________________________________
      if (efficiencyFlags.applyEfficiencyCorrection) {
        initEfficiencyFromCCDB(bc);
      }
      // ________________________________________________
      // skip if desired trigger not found
      if (triggerPresenceMap.size() > 0 && (!TESTBIT(triggerPresenceMap[collision1.globalIndex()], triggerBinToSelect) || !TESTBIT(triggerPresenceMap[collision2.globalIndex()], triggerBinToSelect))) {
        return;
      }

      // ________________________________________________
      // Perform basic event selection on both collisions
      if (!isisCollisionSelect(collision1) || !isisCollisionSelect(collision2)) {
        continue;
      }
      if (collision1.centFT0M() > axisRanges[5][1] || collision1.centFT0M() < axisRanges[5][0])
        continue;
      if (collision2.centFT0M() > axisRanges[5][1] || collision2.centFT0M() < axisRanges[5][0])
        continue;
      if (doMixingQAandEventQA) {
        if (collision1.globalIndex() == collision2.globalIndex()) {
          histos.fill(HIST("MixingQA/hMixingQA"), 0.0f); // same-collision pair counting
        }
        histos.fill(HIST("MixingQA/hMEpvz1"), collision1.posZ());
        histos.fill(HIST("MixingQA/hMEpvz2"), collision2.posZ());
        histos.fill(HIST("MixingQA/hMECollisionBins"), colBinning.getBin({collision1.posZ(), collision1.centFT0M()}));
      }
      // ________________________________________________
      // Do slicing
      auto slicedTriggerTracks = triggerTracks.sliceBy(collisionSliceTracks, collision1.globalIndex());
      auto slicedAssocHadrons = assocHadrons.sliceBy(collisionSliceHadrons, collision2.globalIndex());
      // ________________________________________________
      // Do hadron - hadron correlations
      if (masterConfigurations.doFullCorrelationStudy)
        fillCorrelationsHadron(slicedTriggerTracks, slicedAssocHadrons, true, collision1.posZ(), collision1.centFT0M(), bField);
    }
  }

  void processMixedEventHPions(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::PVMults> const& collisions, soa::Join<aod::AssocHadrons, aod::AssocPID> const& assocPions, soa::Join<aod::TriggerTracks, aod::TriggerTrackExtras> const& triggerTracks, TracksComplete const&, aod::BCsWithTimestamps const&)
  {
    BinningTypePP colBinning{{axesConfigurations.axisVtxZ, axesConfigurations.axisMult}, true};
    for (auto const& [collision1, collision2] : soa::selfCombinations(colBinning, masterConfigurations.mixingParameter, -1, collisions, collisions)) {
      auto bc = collision1.bc_as<aod::BCsWithTimestamps>();
      auto bField = getMagneticField(bc.timestamp());
      // ________________________________________________
      if (efficiencyFlags.applyEfficiencyCorrection) {
        initEfficiencyFromCCDB(bc);
      }
      // ________________________________________________
      // skip if desired trigger not found
      if (triggerPresenceMap.size() > 0 && (!TESTBIT(triggerPresenceMap[collision1.globalIndex()], triggerBinToSelect) || !TESTBIT(triggerPresenceMap[collision2.globalIndex()], triggerBinToSelect))) {
        continue;
      }

      // ________________________________________________
      // Perform basic event selection on both collisions
      if (!isisCollisionSelect(collision1) || !isisCollisionSelect(collision2)) {
        continue;
      }
      if (collision1.centFT0M() > axisRanges[5][1] || collision1.centFT0M() < axisRanges[5][0])
        continue;
      if (collision2.centFT0M() > axisRanges[5][1] || collision2.centFT0M() < axisRanges[5][0])
        continue;
      if (doMixingQAandEventQA) {
        if (collision1.globalIndex() == collision2.globalIndex()) {
          histos.fill(HIST("MixingQA/hMixingQA"), 0.0f); // same-collision pair counting
        }
        histos.fill(HIST("MixingQA/hMEpvz1"), collision1.posZ());
        histos.fill(HIST("MixingQA/hMEpvz2"), collision2.posZ());
        histos.fill(HIST("MixingQA/hMECollisionBins"), colBinning.getBin({collision1.posZ(), collision1.centFT0M()}));
      }
      // ________________________________________________
      // Do slicing
      auto slicedTriggerTracks = triggerTracks.sliceBy(collisionSliceTracks, collision1.globalIndex());
      auto slicedAssocPions = assocPions.sliceBy(collisionSliceHadrons, collision2.globalIndex());
      // ________________________________________________
      // Do hadron - phi correlations
      if (masterConfigurations.doFullCorrelationStudy)
        fillCorrelationsHadron(slicedTriggerTracks, slicedAssocPions, true, collision1.posZ(), collision1.centFT0M(), bField);
    }
  }

  void processMixedEventHPhis(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::PVMults> const& collisions, aod::AssocPhis const& assocPhis, aod::TriggerTracks const& triggerTracks, TracksComplete const&, aod::BCsWithTimestamps const&)
  {
    BinningTypePP colBinning{{axesConfigurations.axisVtxZ, axesConfigurations.axisMult}, true};

    for (auto const& [collision1, collision2] : soa::selfCombinations(colBinning, masterConfigurations.mixingParameter, -1, collisions, collisions)) {

      auto bc = collision1.bc_as<aod::BCsWithTimestamps>();
      auto bField = getMagneticField(bc.timestamp());

      if (efficiencyFlags.applyEfficiencyCorrection) {
        initEfficiencyFromCCDB(bc);
      }

      if (triggerPresenceMap.size() > 0 && (!TESTBIT(triggerPresenceMap[collision1.globalIndex()], triggerBinToSelect) || !TESTBIT(triggerPresenceMap[collision2.globalIndex()], triggerBinToSelect))) {
        continue;
      }

      if (!isisCollisionSelect(collision1) || !isisCollisionSelect(collision2)) {
        continue;
      }

      if (collision1.centFT0M() > axisRanges[5][1] || collision1.centFT0M() < axisRanges[5][0]) {
        continue;
      }

      if (collision2.centFT0M() > axisRanges[5][1] || collision2.centFT0M() < axisRanges[5][0]) {
        continue;
      }

      if (doMixingQAandEventQA) {
        histos.fill(HIST("MixingQA/hMEpvz1"), collision1.posZ());
        histos.fill(HIST("MixingQA/hMEpvz2"), collision2.posZ());
        histos.fill(HIST("MixingQA/hMECollisionBins"), colBinning.getBin({collision1.posZ(), collision1.centFT0M()}));
      }

      auto slicedTriggerTracks = triggerTracks.sliceBy(collisionSliceTracks, collision1.globalIndex());

      auto slicedAssocPhis = assocPhis.sliceBy(collisionSlicePhis, collision2.globalIndex());

      if (masterConfigurations.doFullCorrelationStudy) {
        // Use the per-collision slices computed above, not the full unsliced
        // triggerTracks/assocPhis tables -- passing the whole dataframe here
        // would mix every trigger with every candidate in the TF regardless
        // of which two collisions are actually being paired.
        fillCorrelationsPhi(slicedTriggerTracks, slicedAssocPhis, true, false, collision1.posZ(), collision1.centFT0M(), bField);
      }
    }
  }

  void processMixedEventHKstars(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::PVMults> const& collisions, aod::AssocKstars const& assocKstars, aod::TriggerTracks const& triggerTracks, TracksComplete const&, aod::BCsWithTimestamps const&)
  {
    BinningTypePP colBinning{{axesConfigurations.axisVtxZ, axesConfigurations.axisMult}, true};

    for (auto const& [collision1, collision2] : soa::selfCombinations(colBinning, masterConfigurations.mixingParameter, -1, collisions, collisions)) {

      auto bc = collision1.bc_as<aod::BCsWithTimestamps>();
      auto bField = getMagneticField(bc.timestamp());

      if (efficiencyFlags.applyEfficiencyCorrection) {
        initEfficiencyFromCCDB(bc);
      }

      if (triggerPresenceMap.size() > 0 && (!TESTBIT(triggerPresenceMap[collision1.globalIndex()], triggerBinToSelect) || !TESTBIT(triggerPresenceMap[collision2.globalIndex()], triggerBinToSelect))) {
        continue;
      }

      if (!isisCollisionSelect(collision1) || !isisCollisionSelect(collision2)) {
        continue;
      }

      if (collision1.centFT0M() > axisRanges[5][1] || collision1.centFT0M() < axisRanges[5][0]) {
        continue;
      }

      if (collision2.centFT0M() > axisRanges[5][1] || collision2.centFT0M() < axisRanges[5][0]) {
        continue;
      }

      if (doMixingQAandEventQA) {
        histos.fill(HIST("MixingQA/hMEpvz1"), collision1.posZ());
        histos.fill(HIST("MixingQA/hMEpvz2"), collision2.posZ());
        histos.fill(HIST("MixingQA/hMECollisionBins"), colBinning.getBin({collision1.posZ(), collision1.centFT0M()}));
      }

      auto slicedTriggerTracks = triggerTracks.sliceBy(collisionSliceTracks, collision1.globalIndex());

      auto slicedAssocKstars = assocKstars.sliceBy(collisionSliceKstars, collision2.globalIndex());

      if (masterConfigurations.doFullCorrelationStudy) {
        // Same reasoning as processMixedEventHPhis above -- use this collision
        // pair's own slices, not the full unsliced triggerTracks/assocKstars.
        fillCorrelationsKstar(slicedTriggerTracks, slicedAssocKstars, true, false, collision1.posZ(), collision1.centFT0M(), bField);
      }
    }
  }

  void processMCGenerated(aod::McCollision const& /*mcCollision*/, soa::SmallGroups<soa::Join<aod::McCollisionLabels, aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::CentFT0Cs, aod::PVMults>> const& collisions, aod::McParticles const& mcParticles)
  {
    histos.fill(HIST("hClosureTestEventCounter"), 2.5f);

    //---------------------------------------------
    // Generated particles before event selection
    //---------------------------------------------

    for (auto const& mcParticle : mcParticles) {

      if (std::abs(mcParticle.eta()) > etaSel)
        continue;

      float pt = mcParticle.pt();

      // trigger particles
      if (std::abs(mcParticle.pdgCode()) == PDG_t::kPiPlus || std::abs(mcParticle.pdgCode()) == PDG_t::kKPlus || std::abs(mcParticle.pdgCode()) == PDG_t::kProton || std::abs(mcParticle.pdgCode()) == PDG_t::kElectron || std::abs(mcParticle.pdgCode()) == PDG_t::kMuonMinus) {
        if (!masterConfigurations.doTriggPhysicalPrimary || mcParticle.isPhysicalPrimary()) {
          histos.fill(HIST("hGeneratedQAPtTrigger"), pt, 0.0f);
        }
      }

      // generated phi

      if (std::abs(mcParticle.pdgCode()) == Pdg::kPhi && masterConfigurations.doCorrelationPhi) {

        // if (!doAssocPhysicalPrimaryInGen ||
        //   mcParticle.isPhysicalPrimary()) {

        histos.fill(HIST("hGeneratedQAPtAssociatedPhi"), pt, 0.0f);

        //}
      }
    }

    for (auto const& mcParticle : mcParticles) {

      // phi
      if (std::abs(mcParticle.pdgCode()) == Pdg::kPhi) {
        if (efficiencyFlags.applyEffAsFunctionOfMultAndPhi) {
          histos.fill(HIST("Generated/hPhi"), mcParticle.pt(), mcParticle.eta(), mcParticle.phi(), 1.0);
        } else {
          histos.fill(HIST("Generated/hPhi"), mcParticle.pt(), mcParticle.eta(), 1.0);
        }
      }

      if (doAssocPhysicalPrimaryInGen && !mcParticle.isPhysicalPrimary())
        continue;
      // NOTE: intentionally restricted to Pion only (IndexPion..IndexPion); Kstar0 is
      // not extended to processMCGenerated/ClosureTest/Prediction in this version.
      static_for<IndexPion, IndexPion>([&](auto i) {
        constexpr int Index = i.value;
        if (i == IndexPion && mcParticle.pdgCode() > Neutral) {
          histos.fill(HIST("Generated/hPositive") + HIST(Particlenames[Index]), mcParticle.pt(), mcParticle.eta(), 1);
        } else if (i == IndexPion && mcParticle.pdgCode() < Neutral) {
          histos.fill(HIST("Generated/hNegative") + HIST(Particlenames[Index]), mcParticle.pt(), mcParticle.eta(), 1);
        } else if (mcParticle.pdgCode() == PdgCodes[i]) {
          if (efficiencyFlags.applyEffAsFunctionOfMultAndPhi) {
            histos.fill(HIST("Generated/h") + HIST(Particlenames[Index]), mcParticle.pt(), mcParticle.eta(), mcParticle.phi(), 1);
          } else {
            histos.fill(HIST("Generated/h") + HIST(Particlenames[Index]), mcParticle.pt(), mcParticle.eta(), 1);
          }
        }
      });
    }

    //---------------------------------------------
    // Require reconstructed collision
    //---------------------------------------------

    if (collisions.size() < 1)
      return;

    int biggestNContribs = -1;
    float bestCollisionVtxZ = 0.f;
    float bestCollisionFT0Mpercentile = -1.f;
    float bestCollisionFT0Cpercentile = -1.f;
    bool bestCollisionSel8 = false;
    bool bestCollisionINELgtZERO = false;
    bool isCollisionSelect = false;
    uint32_t bestCollisionTriggerPresenceMap = 0;

    for (auto const& collision : collisions) {

      if (collision.numContrib() <= biggestNContribs)
        continue;

      biggestNContribs = collision.numContrib();

      bestCollisionFT0Mpercentile = collision.centFT0M();
      bestCollisionFT0Cpercentile = collision.centFT0C();

      if (masterConfigurations.applyNewMCSelection) {
        isCollisionSelect = (masterConfigurations.doPPAnalysis && isisCollisionSelect(collision)) || (!masterConfigurations.doPPAnalysis && isisCollisionSelectPbPb(collision, false));
      } else {
        bestCollisionSel8 = collision.sel8();
        bestCollisionVtxZ = collision.posZ();
        bestCollisionINELgtZERO = collision.isInelGt0();
      }

      if (triggerPresenceMap.size() > 0)
        bestCollisionTriggerPresenceMap = triggerPresenceMap[collision.globalIndex()];
    }

    if (collisions.size() > 1) {
      for (auto const& mcParticle : mcParticles) {

        if (std::abs(mcParticle.y()) > ySel)
          continue;
        if (std::abs(mcParticle.pdgCode()) == Pdg::kPhi) {
          histos.fill(HIST("GeneratedWithPV/hPhi") + HIST("_MidYVsMult_TwoPVsOrMore"), mcParticle.pt(), bestCollisionFT0Mpercentile);
        }

        if (doAssocPhysicalPrimaryInGen && !mcParticle.isPhysicalPrimary())
          continue;

        // NOTE: Pion only, see note above; Kstar0 not extended here.
        static_for<IndexPion, IndexPion>([&](auto i) {
          constexpr int Index = i.value;

          if (std::abs(mcParticle.pdgCode()) == std::abs(PdgCodes[i]))
            histos.fill(HIST("GeneratedWithPV/h") + HIST(Particlenames[Index]) + HIST("_MidYVsMult_TwoPVsOrMore"), mcParticle.pt(), bestCollisionFT0Mpercentile);
        });
      }
    }

    //---------------------------------------------
    // Collision selection
    //---------------------------------------------

    if (triggerPresenceMap.size() > 0 && !TESTBIT(bestCollisionTriggerPresenceMap, triggerBinToSelect))
      return;

    if (masterConfigurations.applyNewMCSelection) {

      if (!isCollisionSelect)
        return;

    } else {

      if (!bestCollisionSel8)
        return;

      if (std::abs(bestCollisionVtxZ) > masterConfigurations.zVertexCut)
        return;

      if (!bestCollisionINELgtZERO)
        return;
    }

    histos.fill(HIST("hClosureTestEventCounter"), 3.5f);

    //---------------------------------------------
    // Generated particles after event selection
    //---------------------------------------------

    for (auto const& mcParticle : mcParticles) {

      if (std::abs(mcParticle.eta()) > etaSel)
        continue;

      float pt = mcParticle.pt();

      // trigger QA

      if (std::abs(mcParticle.pdgCode()) == PDG_t::kPiPlus || std::abs(mcParticle.pdgCode()) == PDG_t::kKPlus || std::abs(mcParticle.pdgCode()) == PDG_t::kProton || std::abs(mcParticle.pdgCode()) == PDG_t::kElectron || std::abs(mcParticle.pdgCode()) == PDG_t::kMuonMinus) {

        if (!masterConfigurations.doTriggPhysicalPrimary || mcParticle.isPhysicalPrimary()) {
          histos.fill(HIST("hGeneratedQAPtTrigger"), pt, 1.0f);
        }
      }

      // generated phi

      if (std::abs(mcParticle.pdgCode()) != Pdg::kPhi)
        continue;

      // if (doAssocPhysicalPrimaryInGen &&
      //   !mcParticle.isPhysicalPrimary())
      // continue;

      histos.fill(HIST("hGeneratedQAPtAssociatedPhi"), pt, 1.0f);
    }

    for (auto const& mcParticle : mcParticles) {

      double eta = mcParticle.eta();
      double pt = mcParticle.pt();

      // generated phi

      if (std::abs(mcParticle.pdgCode()) == Pdg::kPhi) {

        if (efficiencyFlags.applyEffAsFunctionOfMultAndPhi) {

          histos.fill(HIST("GeneratedWithPV/hPhi"), pt, mcParticle.eta(), mcParticle.phi(), bestCollisionFT0Cpercentile);

        } else {
          histos.fill(HIST("GeneratedWithPV/hPhi"), pt, mcParticle.eta(), bestCollisionFT0Mpercentile);
        }
        if (std::abs(mcParticle.y()) < ySel)
          histos.fill(HIST("GeneratedWithPV/hPhi") + HIST("_MidYVsMult"), pt, bestCollisionFT0Mpercentile);
      }

      if (doAssocPhysicalPrimaryInGen && !mcParticle.isPhysicalPrimary()) {
        continue;
      }
      if (std::abs(mcParticle.pdgCode()) == PDG_t::kPiPlus || std::abs(mcParticle.pdgCode()) == PDG_t::kKPlus || std::abs(mcParticle.pdgCode()) == PDG_t::kProton || std::abs(mcParticle.pdgCode()) == PDG_t::kElectron || std::abs(mcParticle.pdgCode()) == PDG_t::kMuonMinus) {
        if (efficiencyFlags.applyEffAsFunctionOfMultAndPhi) {
          histos.fill(HIST("GeneratedWithPV/hTrigger"), pt, eta, mcParticle.phi(), bestCollisionFT0Cpercentile);
        } else {
          histos.fill(HIST("GeneratedWithPV/hTrigger"), pt, eta, bestCollisionFT0Mpercentile);
        }
        if (mcParticle.pdgCode() > 0)
          histos.fill(HIST("GeneratedWithPV/hPositiveTrigger"), pt, eta, bestCollisionFT0Mpercentile);
        else
          histos.fill(HIST("GeneratedWithPV/hNegativeTrigger"), pt, eta, bestCollisionFT0Mpercentile);
      }

      // NOTE: Pion only, see note above; Kstar0 not extended here.
      static_for<IndexPion, IndexPion>([&](auto i) {
        constexpr int Index = i.value;
        if (i == IndexPion && std::abs(mcParticle.pdgCode()) == std::abs(PdgCodes[i]) && mcParticle.pdgCode() > Neutral) {
          histos.fill(HIST("GeneratedWithPV/hPositive") + HIST(Particlenames[Index]), mcParticle.pt(), mcParticle.eta(), bestCollisionFT0Mpercentile);
        } else if (i == IndexPion && std::abs(mcParticle.pdgCode()) == std::abs(PdgCodes[i]) && mcParticle.pdgCode() < Neutral) {
          histos.fill(HIST("GeneratedWithPV/hNegative") + HIST(Particlenames[Index]), mcParticle.pt(), mcParticle.eta(), bestCollisionFT0Mpercentile);
        }

        if (std::abs(mcParticle.pdgCode()) == std::abs(PdgCodes[i])) {
          if (efficiencyFlags.applyEffAsFunctionOfMultAndPhi) {
            histos.fill(HIST("GeneratedWithPV/h") + HIST(Particlenames[Index]), pt, eta, mcParticle.phi(), bestCollisionFT0Cpercentile);
          } else {
            histos.fill(HIST("GeneratedWithPV/h") + HIST(Particlenames[Index]), pt, eta, bestCollisionFT0Mpercentile);
          }
          if (std::abs(mcParticle.y()) < ySel)
            histos.fill(HIST("GeneratedWithPV/h") + HIST(Particlenames[Index]) + HIST("_MidYVsMult"), pt, bestCollisionFT0Mpercentile);
        }
      });
    }
  }

  void processClosureTest(aod::McCollision const& /*mcCollision*/, soa::SmallGroups<soa::Join<aod::McCollisionLabels, aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::PVMults>> const& recCollisions, aod::McParticles const& mcParticles)
  {

    std::vector<uint32_t> triggerIndices;
    std::vector<std::vector<uint32_t>> associatedIndices;
    std::vector<uint32_t> assocHadronIndices;
    std::vector<uint32_t> piIndices;
    std::vector<uint32_t> phiIndices;
    std::vector<uint32_t> kstarIndices; // intentionally left unpopulated: Kstar0 is not
                                        // extended to this process function in this version

    for (auto const& mcParticle : mcParticles) {
      double geta = mcParticle.eta();
      if (std::abs(geta) > etaSel) {
        continue;
      }
      double gpt = mcParticle.pt();
      if (std::abs(mcParticle.pdgCode()) == PDG_t::kPiPlus || std::abs(mcParticle.pdgCode()) == PDG_t::kKPlus || std::abs(mcParticle.pdgCode()) == PDG_t::kProton || std::abs(mcParticle.pdgCode()) == PDG_t::kElectron || std::abs(mcParticle.pdgCode()) == PDG_t::kMuonMinus) {
        if (!masterConfigurations.doTriggPhysicalPrimary || mcParticle.isPhysicalPrimary()) {
          histos.fill(HIST("hClosureQAPtTrigger"), gpt, 0.0f); // step 1: no event selection whatsoever
        }
      }

      // if (!doAssocPhysicalPrimary || mcParticle.isPhysicalPrimary()) {
      if (std::abs(mcParticle.pdgCode()) == Pdg::kPhi && masterConfigurations.doCorrelationPhi) {
        histos.fill(HIST("hClosureQAPtAssociatedPhi"), gpt, 0.0f); // step 1: no event selection whatsoever
      }
      // }
    }

    histos.fill(HIST("hClosureTestEventCounter"), 0.5f);

    int bestCollisionFT0Mpercentile = -1;
    float bestCollisionVtxZ = 0.0f;
    bool bestCollisionSel8 = false;
    bool bestCollisionINELgtZERO = false;
    int biggestNContribs = -1;
    uint32_t bestCollisionTriggerPresenceMap = 0;

    for (auto const& recCollision : recCollisions) {
      if (biggestNContribs < recCollision.numContrib()) {
        biggestNContribs = recCollision.numContrib();
        bestCollisionFT0Mpercentile = recCollision.centFT0M();
        bestCollisionSel8 = recCollision.sel8();
        bestCollisionVtxZ = recCollision.posZ();
        bestCollisionINELgtZERO = recCollision.isInelGt0();
        if (triggerPresenceMap.size() > 0)
          bestCollisionTriggerPresenceMap = triggerPresenceMap[recCollision.globalIndex()];
      }
    }
    // ________________________________________________
    // skip if desired trigger not found
    if (triggerPresenceMap.size() > 0 && !TESTBIT(bestCollisionTriggerPresenceMap, triggerBinToSelect)) {
      return;
    }

    if (masterConfigurations.doGenEventSelection) {
      if (!bestCollisionSel8)
        return;
      if (std::abs(bestCollisionVtxZ) > masterConfigurations.zVertexCut)
        return;
      if (!bestCollisionINELgtZERO)
        return;
      if (bestCollisionFT0Mpercentile > axisRanges[5][1] || bestCollisionFT0Mpercentile < axisRanges[5][0]) {
        return;
      }
    }

    histos.fill(HIST("hClosureTestEventCounter"), 1.5f);

    for (auto const& mcParticle : mcParticles) {
      double geta = mcParticle.eta();
      if (std::abs(geta) > etaSel) {
        continue;
      }
      double gpt = mcParticle.pt();
      if (std::abs(mcParticle.pdgCode()) == PDG_t::kPiPlus || std::abs(mcParticle.pdgCode()) == PDG_t::kKPlus || std::abs(mcParticle.pdgCode()) == PDG_t::kProton || std::abs(mcParticle.pdgCode()) == PDG_t::kElectron || std::abs(mcParticle.pdgCode()) == PDG_t::kMuonMinus) {
        if (!masterConfigurations.doTriggPhysicalPrimary || mcParticle.isPhysicalPrimary()) {
          histos.fill(HIST("hClosureQAPtTrigger"), gpt, 1.0f); // step 2: after event selection
        }
      }

      // if (!doAssocPhysicalPrimary || mcParticle.isPhysicalPrimary()) {
      if (std::abs(mcParticle.pdgCode()) == Pdg::kPhi && masterConfigurations.doCorrelationPhi) {
        histos.fill(HIST("hClosureQAPtAssociatedPhi"), gpt, 1.0f); // step 2: after event selection
      }
      // }
    }

    int iteratorNum = -1;
    for (auto const& mcParticle : mcParticles) {
      iteratorNum = iteratorNum + 1;
      double geta = mcParticle.eta();
      double gpt = mcParticle.pt();
      double gphi = mcParticle.phi();
      if (std::abs(geta) > etaSel) {
        continue;
      }
      if (std::abs(mcParticle.pdgCode()) == PDG_t::kPiPlus || std::abs(mcParticle.pdgCode()) == PDG_t::kKPlus || std::abs(mcParticle.pdgCode()) == PDG_t::kProton || std::abs(mcParticle.pdgCode()) == PDG_t::kElectron || std::abs(mcParticle.pdgCode()) == PDG_t::kMuonMinus) {
        if (!masterConfigurations.doTriggPhysicalPrimary || mcParticle.isPhysicalPrimary()) {
          triggerIndices.emplace_back(iteratorNum);
          histos.fill(HIST("ClosureTest/hTrigger"), gpt, geta, bestCollisionFT0Mpercentile);
        }
        if (masterConfigurations.doCorrelationHadron) {
          if (!doAssocPhysicalPrimary || mcParticle.isPhysicalPrimary()) {
            assocHadronIndices.emplace_back(iteratorNum);
            histos.fill(HIST("ClosureTest/hHadron"), gpt, geta, gphi);
          }
        }
      }
      if (!doAssocPhysicalPrimary || mcParticle.isPhysicalPrimary()) {
        if (std::abs(mcParticle.pdgCode()) == PDG_t::kPiPlus && masterConfigurations.doCorrelationPion) {
          piIndices.emplace_back(iteratorNum);
          histos.fill(HIST("ClosureTest/hPion"), gpt, geta, gphi);
        }
      }
      if (std::abs(mcParticle.pdgCode()) == Pdg::kPhi && masterConfigurations.doCorrelationPhi) {
        phiIndices.emplace_back(iteratorNum);
        histos.fill(HIST("ClosureTest/hPhi"), gpt, geta, gphi);
      }
    }

    associatedIndices.emplace_back(phiIndices);         // IndexPhi   = 0
    associatedIndices.emplace_back(kstarIndices);       // IndexKstar = 1 (always empty here)
    associatedIndices.emplace_back(piIndices);          // IndexPion  = 2
    associatedIndices.emplace_back(assocHadronIndices); // Hadron = 3

    for (std::size_t iTrigger = 0; iTrigger < triggerIndices.size(); iTrigger++) {
      auto triggerParticle = mcParticles.iteratorAt(triggerIndices[iTrigger]);
      // check range of trigger particle
      if (triggerParticle.pt() > axisRanges[3][1] || triggerParticle.pt() < axisRanges[3][0]) {
        continue;
      }
      double getatrigger = triggerParticle.eta();
      double gphitrigger = triggerParticle.phi();
      double pttrigger = triggerParticle.pt();

      auto mothers = triggerParticle.mothers_as<aod::McParticles>();
      auto globalIndex = (mothers.size() > 0) ? triggerParticle.mothers_first_as<aod::McParticles>().globalIndex() : -1;

      static_for<0, 3>([&](auto i) { // associated loop (Kstar0 slot is always empty, see above)
        constexpr int Index = i.value;
        for (std::size_t iassoc = 0; iassoc < associatedIndices[Index].size(); iassoc++) {
          auto assocParticle = mcParticles.iteratorAt(associatedIndices[Index][iassoc]);
          if (triggerIndices[iTrigger] != associatedIndices[Index][iassoc] && globalIndex != assocParticle.globalIndex()) { // avoid self
            double getaassoc = assocParticle.eta();
            double gphiassoc = assocParticle.phi();
            double ptassoc = assocParticle.pt();
            double deltaphi = computeDeltaPhi(gphitrigger, gphiassoc);
            double deltaeta = getatrigger - getaassoc;

            // skip if basic ranges not met
            if (deltaphi < axisRanges[0][0] || deltaphi > axisRanges[0][1])
              continue;
            if (deltaeta < axisRanges[1][0] || deltaeta > axisRanges[1][1])
              continue;
            if (ptassoc < axisRanges[2][0] || ptassoc > axisRanges[2][1])
              continue;
            if (TESTBIT(doCorrelation, i))
              histos.fill(HIST("ClosureTest/sameEvent/") + HIST(Particlenames[Index]), computeDeltaPhi(gphitrigger, gphiassoc), deltaeta, ptassoc, pttrigger, bestCollisionVtxZ, bestCollisionFT0Mpercentile);
          }
        }
      });
    }
  }

  void processMixedEventHPhisInBuffer(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::PVMults> const& collisions, aod::TriggerTracks const& triggerTracks, aod::AssocPhis const& assocPhis, TracksComplete const&, aod::BCsWithTimestamps const&)
  {
    for (auto const& collision : collisions) {

      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      auto dBz = getMagneticField(bc.timestamp());

      // Slice down to this collision's own triggers/candidates before
      // buffering into ValidCollision -- without this, every collision would
      // buffer the entire dataframe's triggers/candidates instead of just
      // its own, defeating the point of per-collision buffering.
      auto slicedTriggerTracks = triggerTracks.sliceBy(collisionSliceTracks, collision.globalIndex());
      auto slicedAssocPhis = assocPhis.sliceBy(collisionSlicePhis, collision.globalIndex());

      fillCorrelationsPhi(slicedTriggerTracks, slicedAssocPhis, true, true, collision.posZ(), collision.centFT0M(), dBz);
    }
  }

  void processMixedEventHKstarsInBuffer(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::PVMults> const& collisions, aod::TriggerTracks const& triggerTracks, aod::AssocKstars const& assocKstars, TracksComplete const&, aod::BCsWithTimestamps const&)
  {
    for (auto const& collision : collisions) {

      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      auto dBz = getMagneticField(bc.timestamp());

      auto slicedTriggerTracks = triggerTracks.sliceBy(collisionSliceTracks, collision.globalIndex());
      auto slicedAssocKstars = assocKstars.sliceBy(collisionSliceKstars, collision.globalIndex());

      fillCorrelationsKstar(slicedTriggerTracks, slicedAssocKstars, true, true, collision.posZ(), collision.centFT0M(), dBz);
    }
  }

  void processPrediction(soa::Join<aod::McCollisions, aod::McCentFT0Ms, aod::McCentFT0Cs, aod::McCentFT0As>::iterator const& mcCollision, aod::McParticles const& mcParticles)
  {
    std::vector<uint32_t> triggerIndices;
    std::vector<std::vector<uint32_t>> associatedIndices;
    std::vector<uint32_t> assocHadronIndices;
    std::vector<uint32_t> piIndices;
    std::vector<uint32_t> phiIndices;
    std::vector<uint32_t> kstarIndices; // intentionally left unpopulated: Kstar0 is not
                                        // extended to this process function in this version
    float centMultFT0M = -1;
    float centMultFT0A = -1;
    float centMultFT0C = -1;
    float multFT0M = -1;
    float multFT0A = -1;
    float multFT0C = -1;
    float multEta08 = -1;
    float multEta05 = -1;
    histos.fill(HIST("Prediction/hEventSelection"), 0.5);
    if (masterConfigurations.selectINELgtZERO && !o2::pwglf::isINELgt0mc(mcParticles, pdgDB)) {
      return;
    }
    histos.fill(HIST("Prediction/hEventSelection"), 1.5);
    if (std::abs(mcCollision.posZ()) > masterConfigurations.zVertexCut) {
      return;
    }
    histos.fill(HIST("Prediction/hEventSelection"), 2.5);
    if (masterConfigurations.useCentralityinPrediction) {
      centMultFT0M = mcCollision.centFT0M();
      centMultFT0A = mcCollision.centFT0A();
      centMultFT0C = mcCollision.centFT0C();
    } else {
      multFT0M = mCounter.countFT0A(mcParticles) + mCounter.countFT0C(mcParticles);
      multFT0A = mCounter.countFT0A(mcParticles);
      multFT0C = mCounter.countFT0C(mcParticles);
      multEta08 = mCounter.countEta08(mcParticles);
      multEta05 = mCounter.countEta05(mcParticles);
    }
    if (!masterConfigurations.useCentralityinPrediction) {
      if (masterConfigurations.doSeparateFT0Prediction) {
        histos.fill(HIST("Prediction/hFT0AvsNchEta08"), multFT0A, multEta08);
        histos.fill(HIST("Prediction/hFT0AvsNchEta05"), multFT0A, multEta05);
        histos.fill(HIST("Prediction/hFT0CvsNchEta08"), multFT0C, multEta08);
        histos.fill(HIST("Prediction/hFT0CvsNchEta05"), multFT0C, multEta05);
      }
      histos.fill(HIST("Prediction/hFT0MvsNchEta08"), multFT0M, multEta08);
      histos.fill(HIST("Prediction/hFT0MvsNchEta05"), multFT0M, multEta05);
    }
    int iteratorNum = -1;
    for (auto const& mcParticle : mcParticles) {
      iteratorNum = iteratorNum + 1;
      double geta = mcParticle.eta();
      double gpt = mcParticle.pt();
      double gphi = mcParticle.phi();
      if (std::abs(geta) > etaSel) {
        continue;
      }
      if (std::abs(mcParticle.pdgCode()) == PDG_t::kPiPlus || std::abs(mcParticle.pdgCode()) == PDG_t::kKPlus || std::abs(mcParticle.pdgCode()) == PDG_t::kProton || std::abs(mcParticle.pdgCode()) == PDG_t::kElectron || std::abs(mcParticle.pdgCode()) == PDG_t::kMuonMinus) {
        if (!masterConfigurations.doTriggPhysicalPrimary || mcParticle.isPhysicalPrimary()) {
          triggerIndices.emplace_back(iteratorNum);
          if (masterConfigurations.useCentralityinPrediction) {
            if (masterConfigurations.doSeparateFT0Prediction) {
              histos.fill(HIST("Prediction/hTriggerFT0A"), gpt, geta, centMultFT0A);
              histos.fill(HIST("Prediction/hTriggerFT0C"), gpt, geta, centMultFT0C);
            }
            histos.fill(HIST("Prediction/hTrigger"), gpt, geta, centMultFT0M);
          } else {
            if (masterConfigurations.doSeparateFT0Prediction) {
              histos.fill(HIST("Prediction/hTriggerFT0A"), gpt, geta, multFT0A);
              histos.fill(HIST("Prediction/hTriggerFT0C"), gpt, geta, multFT0C);
            }
            histos.fill(HIST("Prediction/hTrigger"), gpt, geta, multFT0M);
          }
        }
        if (masterConfigurations.doCorrelationHadron) {
          if (!doAssocPhysicalPrimary || mcParticle.isPhysicalPrimary()) {
            assocHadronIndices.emplace_back(iteratorNum);
            histos.fill(HIST("Prediction/hHadron"), gpt, geta, gphi);
          }
        }
      }
      if (!doAssocPhysicalPrimary || mcParticle.isPhysicalPrimary()) {
        if (std::abs(mcParticle.pdgCode()) == PDG_t::kPiPlus && masterConfigurations.doCorrelationPion) {
          piIndices.emplace_back(iteratorNum);
          histos.fill(HIST("Prediction/hPion"), gpt, geta, gphi);
        }
      }
      if (std::abs(mcParticle.pdgCode()) == Pdg::kPhi && masterConfigurations.doCorrelationPhi) {
        phiIndices.emplace_back(iteratorNum);
        histos.fill(HIST("Prediction/hPhi"), gpt, geta, gphi);
      }
    }

    associatedIndices.emplace_back(phiIndices);         // IndexPhi   = 0
    associatedIndices.emplace_back(kstarIndices);       // IndexKstar = 1 (always empty here)
    associatedIndices.emplace_back(piIndices);          // IndexPion  = 2
    associatedIndices.emplace_back(assocHadronIndices); // Hadron = 3
    for (std::size_t iTrigger = 0; iTrigger < triggerIndices.size(); iTrigger++) {
      auto triggerParticle = mcParticles.iteratorAt(triggerIndices[iTrigger]);
      // check range of trigger particle
      if (triggerParticle.pt() > axisRanges[3][1] || triggerParticle.pt() < axisRanges[3][0]) {
        continue;
      }
      double getatrigger = triggerParticle.eta();
      double gphitrigger = triggerParticle.phi();
      double pttrigger = triggerParticle.pt();
      auto mothers = triggerParticle.mothers_as<aod::McParticles>();
      auto globalIndex = (mothers.size() > 0) ? triggerParticle.mothers_first_as<aod::McParticles>().globalIndex() : -1;
      static_for<0, 3>([&](auto i) { // associated loop (Kstar0 slot is always empty, see above)
        constexpr int Index = i.value;
        for (std::size_t iassoc = 0; iassoc < associatedIndices[Index].size(); iassoc++) {
          auto assocParticle = mcParticles.iteratorAt(associatedIndices[Index][iassoc]);
          if (triggerIndices[iTrigger] != associatedIndices[Index][iassoc] && globalIndex != assocParticle.globalIndex()) { // avoid self
            double getaassoc = assocParticle.eta();
            double gphiassoc = assocParticle.phi();
            double ptassoc = assocParticle.pt();
            double deltaphi = computeDeltaPhi(gphitrigger, gphiassoc);
            double deltaeta = getatrigger - getaassoc;

            // skip if basic ranges not met
            if (deltaphi < axisRanges[0][0] || deltaphi > axisRanges[0][1])
              continue;
            if (deltaeta < axisRanges[1][0] || deltaeta > axisRanges[1][1])
              continue;
            if (ptassoc < axisRanges[2][0] || ptassoc > axisRanges[2][1])
              continue;
            if (TESTBIT(doCorrelation, i)) {
              if (masterConfigurations.useCentralityinPrediction) {
                histos.fill(HIST("Prediction/sameEvent/") + HIST(Particlenames[Index]), computeDeltaPhi(gphitrigger, gphiassoc), deltaeta, ptassoc, pttrigger, mcCollision.posZ(), centMultFT0M);
                if (masterConfigurations.doSeparateFT0Prediction) {
                  histos.fill(HIST("Prediction/sameEventFT0A/") + HIST(Particlenames[Index]), computeDeltaPhi(gphitrigger, gphiassoc), deltaeta, ptassoc, pttrigger, mcCollision.posZ(), centMultFT0A);
                  histos.fill(HIST("Prediction/sameEventFT0C/") + HIST(Particlenames[Index]), computeDeltaPhi(gphitrigger, gphiassoc), deltaeta, ptassoc, pttrigger, mcCollision.posZ(), centMultFT0C);
                }
              } else {
                histos.fill(HIST("Prediction/sameEvent/") + HIST(Particlenames[Index]), computeDeltaPhi(gphitrigger, gphiassoc), deltaeta, ptassoc, pttrigger, mcCollision.posZ(), multFT0M);
                if (masterConfigurations.doSeparateFT0Prediction) {
                  histos.fill(HIST("Prediction/sameEventFT0A/") + HIST(Particlenames[Index]), computeDeltaPhi(gphitrigger, gphiassoc), deltaeta, ptassoc, pttrigger, mcCollision.posZ(), multFT0A);
                  histos.fill(HIST("Prediction/sameEventFT0C/") + HIST(Particlenames[Index]), computeDeltaPhi(gphitrigger, gphiassoc), deltaeta, ptassoc, pttrigger, mcCollision.posZ(), multFT0C);
                }
              }
            }
          }
        }
      });
    }
  }

  PROCESS_SWITCH(HResonanceCorrelation, processSelectEventWithTrigger, "Select events with trigger only", true);

  PROCESS_SWITCH(HResonanceCorrelation, processSameEventHPions, "Process same events, h-Pion", true);
  PROCESS_SWITCH(HResonanceCorrelation, processSameEventHHadrons, "Process same events, h-h", true);
  PROCESS_SWITCH(HResonanceCorrelation, processSameEventHPhis, "Process same events, h-phi", true);
  PROCESS_SWITCH(HResonanceCorrelation, processSameEventHKstars, "Process same events, h-K*0", true);

  PROCESS_SWITCH(HResonanceCorrelation, processMixedEventHPions, "Process mixed events, h-Pion", true);
  PROCESS_SWITCH(HResonanceCorrelation, processMixedEventHHadrons, "Process mixed events, h-h", true);
  PROCESS_SWITCH(HResonanceCorrelation, processMixedEventHPhis, "Process mixed events, h-phi", true);
  PROCESS_SWITCH(HResonanceCorrelation, processMixedEventHPhisInBuffer, "Process mixed-event, h-phi using collision buffer", false);
  PROCESS_SWITCH(HResonanceCorrelation, processMixedEventHKstars, "Process mixed events, h-K*0", true);
  PROCESS_SWITCH(HResonanceCorrelation, processMixedEventHKstarsInBuffer, "Process mixed-event, h-K*0 using collision buffer", false);

  PROCESS_SWITCH(HResonanceCorrelation, processMCGenerated, "Process MC generated", false);
  PROCESS_SWITCH(HResonanceCorrelation, processClosureTest, "Process Closure Test", false);
  PROCESS_SWITCH(HResonanceCorrelation, processPrediction, "process model prediction", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HResonanceCorrelation>(cfgc)};
}
