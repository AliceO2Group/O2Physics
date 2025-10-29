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

/// \file strangenessBuilderModule.h
/// \brief strangeness builder module
/// \author ALICE

#ifndef PWGLF_UTILS_STRANGENESSBUILDERMODULE_H_
#define PWGLF_UTILS_STRANGENESSBUILDERMODULE_H_

// simple checkers, but ensure 8 bit integers
#define BITSET(var, nbit) ((var) |= (static_cast<uint8_t>(1) << static_cast<uint8_t>(nbit)))

#include "TableHelper.h"

#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/Utils/strangenessBuilderHelper.h"

#include "Common/Core/TPCVDriftManager.h"

#include "DataFormatsCalibration/MeanVertexObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisHelpers.h"
#include "Framework/Configurable.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/HistogramSpec.h"

#include <array>
#include <cmath>
#include <cstdlib>
#include <memory>
#include <string>
#include <vector>

//__________________________________________
// strangeness builder module

namespace o2
{
namespace pwglf
{
namespace strangenessbuilder // avoid polluting other namespaces
{

// statics necessary for the configurables in this namespace
static constexpr int nParameters = 1;
static const std::vector<std::string> tableNames{
  "V0Indices",          //.0 (standard analysis: V0Cores)
  "V0CoresBase",        //.1 (standard analyses: main table)
  "V0Covs",             //.2 (joinable with V0Cores)
  "CascIndices",        //.3 (standard analyses: CascData)
  "KFCascIndices",      //.4 (standard analyses: KFCascData)
  "TraCascIndices",     //.5 (standard analyses: TraCascData)
  "StoredCascCores",    //.6 (standard analyses: CascData, main table)
  "StoredKFCascCores",  //.7 (standard analyses: KFCascData, main table)
  "StoredTraCascCores", //.8 (standard analyses: TraCascData, main table)
  "CascCovs",           //.9 (joinable with CascData)
  "KFCascCovs",         //.10 (joinable with KFCascData)
  "TraCascCovs",        //.11 (joinable with TraCascData)
  "V0TrackXs",          //.12 (joinable with V0Data)
  "CascTrackXs",        //.13 (joinable with CascData)
  "CascBBs",            //.14 (standard, bachelor-baryon vars)
  "V0DauCovs",          //.15 (requested: tracking studies)
  "V0DauCovIUs",        //.16 (requested: tracking studies)
  "V0TraPosAtDCAs",     //.17 (requested: tracking studies)
  "V0TraPosAtIUs",      //.18 (requested: tracking studies)
  "V0Ivanovs",          //.19 (requested: tracking studies)
  "McV0Labels",         //.20 (MC/standard analysis)
  "V0MCCores",          //.21 (MC, all generated desired V0s)
  "V0CoreMCLabels",     //.22 (MC, refs V0Cores to V0MCCores)
  "V0MCCollRefs",       //.23 (MC, refs V0MCCores to McCollisions)
  "McCascLabels",       //.24 (MC/standard analysis)
  "McKFCascLabels",     //.25 (MC, refs KFCascCores to CascMCCores)
  "McTraCascLabels",    //.26 (MC, refs TraCascCores to CascMCCores)
  "McCascBBTags",       //.27 (MC, joinable with CascCores, tags reco-ed)
  "CascMCCores",        //.28 (MC, all generated desired cascades)
  "CascCoreMCLabels",   //.29 (MC, refs CascCores to CascMCCores)
  "CascMCCollRefs",     // 30 (MC, refs CascMCCores to McCollisions)
  "CascToTraRefs",      //.31 (interlink CascCores -> TraCascCores)
  "CascToKFRefs",       //.32 (interlink CascCores -> KFCascCores)
  "TraToCascRefs",      //.33 (interlink TraCascCores -> CascCores)
  "KFToCascRefs",       //.34 (interlink KFCascCores -> CascCores)
  "V0FoundTags",        //.35 (tags found vs findable V0s in findable mode)
  "CascFoundTags"       //.36 (tags found vs findable Cascades in findable mode)
};

static constexpr int nTablesConst = 37;

static const std::vector<std::string> parameterNames{"enable"};
static const int defaultParameters[nTablesConst][nParameters]{
  {-1},
  {-1},
  {-1},
  {-1},
  {-1},
  {-1},
  {-1},
  {-1},
  {-1},
  {-1}, // index 9
  {-1},
  {-1},
  {-1},
  {-1},
  {-1},
  {-1},
  {-1},
  {-1},
  {-1},
  {-1}, // index 19
  {-1},
  {-1},
  {-1},
  {-1},
  {-1},
  {-1},
  {-1},
  {-1},
  {-1},
  {-1}, // index 29
  {-1},
  {-1},
  {-1},
  {-1},
  {-1},
  {-1},
  {-1}};

// table index : match order above
enum tableIndex { kV0Indices = 0,
                  kV0CoresBase,
                  kV0Covs,
                  kCascIndices,
                  kKFCascIndices,
                  kTraCascIndices,
                  kStoredCascCores,
                  kStoredKFCascCores,
                  kStoredTraCascCores,
                  kCascCovs,
                  kKFCascCovs,
                  kTraCascCovs,
                  kV0TrackXs,
                  kCascTrackXs,
                  kCascBBs,
                  kV0DauCovs,
                  kV0DauCovIUs,
                  kV0TraPosAtDCAs,
                  kV0TraPosAtIUs,
                  kV0Ivanovs,
                  kMcV0Labels,
                  kV0MCCores,
                  kV0CoreMCLabels,
                  kV0MCCollRefs,
                  kMcCascLabels,
                  kMcKFCascLabels,
                  kMcTraCascLabels,
                  kMcCascBBTags,
                  kCascMCCores,
                  kCascCoreMCLabels,
                  kCascMCCollRefs,
                  kCascToTraRefs,
                  kCascToKFRefs,
                  kTraToCascRefs,
                  kKFToCascRefs,
                  kV0FoundTags,
                  kCascFoundTags,
                  nTables };

enum V0PreSelection : uint8_t { selGamma = 0,
                                selK0Short,
                                selLambda,
                                selAntiLambda };

enum CascPreSelection : uint8_t { selXiMinus = 0,
                                  selXiPlus,
                                  selOmegaMinus,
                                  selOmegaPlus };

static constexpr float defaultK0MassWindowParameters[1][4] = {{2.81882e-03, 1.14057e-03, 1.72138e-03, 5.00262e-01}};
static constexpr float defaultLambdaWindowParameters[1][4] = {{1.17518e-03, 1.24099e-04, 5.47937e-03, 3.08009e-01}};
static constexpr float defaultXiMassWindowParameters[1][4] = {{1.43210e-03, 2.03561e-04, 2.43187e-03, 7.99668e-01}};
static constexpr float defaultOmMassWindowParameters[1][4] = {{1.43210e-03, 2.03561e-04, 2.43187e-03, 7.99668e-01}};

static constexpr float defaultLifetimeCuts[1][4] = {{20, 60, 40, 20}};

struct products : o2::framework::ProducesGroup {
  //__________________________________________________
  // V0 tables
  o2::framework::Produces<aod::V0Indices> v0indices; // standard part of V0Datas
  o2::framework::Produces<aod::V0CoresBase> v0cores; // standard part of V0Datas
  o2::framework::Produces<aod::V0Covs> v0covs;       // for decay chain reco

  //__________________________________________________
  // cascade tables
  o2::framework::Produces<aod::CascIndices> cascidx;            // standard part of CascDatas
  o2::framework::Produces<aod::KFCascIndices> kfcascidx;        // standard part of KFCascDatas
  o2::framework::Produces<aod::TraCascIndices> tracascidx;      // standard part of TraCascDatas
  o2::framework::Produces<aod::StoredCascCores> cascdata;       // standard part of CascDatas
  o2::framework::Produces<aod::StoredKFCascCores> kfcascdata;   // standard part of KFCascDatas
  o2::framework::Produces<aod::StoredTraCascCores> tracascdata; // standard part of TraCascDatas
  o2::framework::Produces<aod::CascCovs> casccovs;              // for decay chain reco
  o2::framework::Produces<aod::KFCascCovs> kfcasccovs;          // for decay chain reco
  o2::framework::Produces<aod::TraCascCovs> tracasccovs;        // for decay chain reco

  //__________________________________________________
  // interlink tables
  o2::framework::Produces<aod::V0DataLink> v0dataLink;           // de-refs V0s -> V0Data
  o2::framework::Produces<aod::CascDataLink> cascdataLink;       // de-refs Cascades -> CascData
  o2::framework::Produces<aod::KFCascDataLink> kfcascdataLink;   // de-refs Cascades -> KFCascData
  o2::framework::Produces<aod::TraCascDataLink> tracascdataLink; // de-refs Cascades -> TraCascData

  //__________________________________________________
  // secondary auxiliary tables
  o2::framework::Produces<aod::V0TrackXs> v0trackXs;     // for decay chain reco
  o2::framework::Produces<aod::CascTrackXs> cascTrackXs; // for decay chain reco

  //__________________________________________________
  // further auxiliary / optional if desired
  o2::framework::Produces<aod::CascBBs> cascbb;
  o2::framework::Produces<aod::V0DauCovs> v0daucovs;            // covariances of daughter tracks
  o2::framework::Produces<aod::V0DauCovIUs> v0daucovIUs;        // covariances of daughter tracks
  o2::framework::Produces<aod::V0TraPosAtDCAs> v0dauPositions;  // auxiliary debug information
  o2::framework::Produces<aod::V0TraPosAtIUs> v0dauPositionsIU; // auxiliary debug information
  o2::framework::Produces<aod::V0Ivanovs> v0ivanovs;            // information for Marian's tests

  //__________________________________________________
  // MC information: V0
  o2::framework::Produces<aod::McV0Labels> v0labels;           // MC labels for V0s
  o2::framework::Produces<aod::V0MCCores> v0mccores;           // mc info storage
  o2::framework::Produces<aod::V0CoreMCLabels> v0CoreMCLabels; // interlink V0Cores -> V0MCCores
  o2::framework::Produces<aod::V0MCCollRefs> v0mccollref;      // references collisions from V0MCCores

  // MC information: Cascades
  o2::framework::Produces<aod::McCascLabels> casclabels;           // MC labels for cascades
  o2::framework::Produces<aod::McKFCascLabels> kfcasclabels;       // MC labels for KF cascades
  o2::framework::Produces<aod::McTraCascLabels> tracasclabels;     // MC labels for tracked cascades
  o2::framework::Produces<aod::McCascBBTags> bbtags;               // bb tags (inv structure tagging in mc)
  o2::framework::Produces<aod::CascMCCores> cascmccores;           // mc info storage
  o2::framework::Produces<aod::CascCoreMCLabels> cascCoreMClabels; // interlink CascCores -> CascMCCores
  o2::framework::Produces<aod::CascMCCollRefs> cascmccollrefs;     // references MC collisions from MC cascades

  //__________________________________________________
  // cascade interlinks
  // FIXME: commented out until strangederivedbuilder adjusted accordingly
  // o2::framework::Produces<aod::CascToTraRefs> cascToTraRefs; // cascades -> tracked
  // o2::framework::Produces<aod::CascToKFRefs> cascToKFRefs;   // cascades -> KF
  // o2::framework::Produces<aod::TraToCascRefs> traToCascRefs; // tracked -> cascades
  // o2::framework::Produces<aod::KFToCascRefs> kfToCascRefs;   // KF -> cascades

  //__________________________________________________
  // Findable tags
  o2::framework::Produces<aod::V0FoundTags> v0FoundTag;
  o2::framework::Produces<aod::CascFoundTags> cascFoundTag;
};

// strangenessBuilder: 1st-order configurables
struct coreConfigurables : o2::framework::ConfigurableGroup {
  o2::framework::Configurable<o2::framework::LabeledArray<int>> enabledTables{"enabledTables",
                                                                              {defaultParameters[0], nTables, nParameters, tableNames, parameterNames},
                                                                              "Produce this table: -1 for autodetect; otherwise, 0/1 is false/true"};
  std::vector<int> mEnabledTables; // Vector of enabled tables

  // first order deduplication implementation
  // more algorithms to be added as necessary
  o2::framework::Configurable<int> deduplicationAlgorithm{"deduplicationAlgorithm", 1, "0: disabled; 1: best pointing angle wins; 2: best DCA daughters wins; 3: best pointing and best DCA wins"};

  // V0 buffer for V0s used in cascades: master switch
  // exchanges CPU (generate V0s again) with memory (save pre-generated V0s)
  o2::framework::Configurable<bool> useV0BufferForCascades{"useV0BufferForCascades", false, "store array of V0s for cascades or not. False (default): save RAM, use more CPU; true: save CPU, use more RAM"};

  o2::framework::Configurable<int> mc_findableMode{"mc_findableMode", 0, "0: disabled; 1: add findable-but-not-found to existing V0s from AO2D; 2: reset V0s and generate only findable-but-not-found"};

  // test the possibility of refitting with material corrections (DCA Fitter option)
  o2::framework::Configurable<bool> refitWithMaterialCorrection{"refitWithMaterialCorrection", false, "do refit after material corrections were applied"};
};

// strangenessBuilder: V0 building options
struct v0Configurables : o2::framework::ConfigurableGroup {
  std::string prefix = "v0BuilderOpts";
  o2::framework::Configurable<bool> generatePhotonCandidates{"generatePhotonCandidates", false, "generate gamma conversion candidates (V0s using TPC-only tracks)"};
  o2::framework::Configurable<bool> moveTPCOnlyTracks{"moveTPCOnlyTracks", true, "if dealing with TPC-only tracks, move them according to TPC drift / time info"};

  // baseline conditionals of V0 building
  o2::framework::Configurable<int> minCrossedRows{"minCrossedRows", 50, "minimum TPC crossed rows for daughter tracks"};
  o2::framework::Configurable<float> dcanegtopv{"dcanegtopv", .1, "DCA Neg To PV"};
  o2::framework::Configurable<float> dcapostopv{"dcapostopv", .1, "DCA Pos To PV"};
  o2::framework::Configurable<double> v0cospa{"v0cospa", 0.95, "V0 CosPA"}; // double -> N.B. dcos(x)/dx = 0 at x=0)
  o2::framework::Configurable<float> dcav0dau{"dcav0dau", 1.0, "DCA V0 Daughters"};
  o2::framework::Configurable<float> v0radius{"v0radius", 0.9, "v0radius"};
  o2::framework::Configurable<float> maxDaughterEta{"maxDaughterEta", 5.0, "Maximum daughter eta (in abs value)"};

  // MC builder options
  o2::framework::Configurable<bool> mc_populateV0MCCoresSymmetric{"mc_populateV0MCCoresSymmetric", false, "populate V0MCCores table for derived data analysis, keep V0MCCores joinable with V0Cores"};
  o2::framework::Configurable<bool> mc_populateV0MCCoresAsymmetric{"mc_populateV0MCCoresAsymmetric", true, "populate V0MCCores table for derived data analysis, create V0Cores -> V0MCCores interlink. Saves only labeled V0s."};
  o2::framework::Configurable<bool> mc_treatPiToMuDecays{"mc_treatPiToMuDecays", true, "if true, will correctly capture pi -> mu and V0 label will still point to originating V0 decay in those cases. Nota bene: prong info will still be for the muon!"};
  o2::framework::Configurable<float> mc_rapidityWindow{"mc_rapidityWindow", 0.5, "rapidity window to save non-recoed candidates"};
  o2::framework::Configurable<bool> mc_keepOnlyPhysicalPrimary{"mc_keepOnlyPhysicalPrimary", false, "Keep only physical primary generated V0s if not recoed"};
  o2::framework::Configurable<bool> mc_addGeneratedK0Short{"mc_addGeneratedK0Short", true, "add V0MCCore entry for generated, not-recoed K0Short"};
  o2::framework::Configurable<bool> mc_addGeneratedLambda{"mc_addGeneratedLambda", true, "add V0MCCore entry for generated, not-recoed Lambda"};
  o2::framework::Configurable<bool> mc_addGeneratedAntiLambda{"mc_addGeneratedAntiLambda", true, "add V0MCCore entry for generated, not-recoed AntiLambda"};
  o2::framework::Configurable<bool> mc_addGeneratedGamma{"mc_addGeneratedGamma", false, "add V0MCCore entry for generated, not-recoed Gamma"};
  o2::framework::Configurable<bool> mc_addGeneratedGammaMakeCollinear{"mc_addGeneratedGammaMakeCollinear", true, "when adding findable gammas, mark them as collinear"};
  o2::framework::Configurable<bool> mc_findableDetachedV0{"mc_findableDetachedV0", false, "if true, generate findable V0s that have collisionId -1. Caution advised."};
};

// strangenessBuilder: cascade building options
struct cascadeConfigurables : o2::framework::ConfigurableGroup {
  std::string prefix = "cascadeBuilderOpts";
  o2::framework::Configurable<bool> useCascadeMomentumAtPrimVtx{"useCascadeMomentumAtPrimVtx", false, "use cascade momentum at PV"};

  // conditionals
  o2::framework::Configurable<int> minCrossedRows{"minCrossedRows", 50, "minimum TPC crossed rows for daughter tracks"};
  o2::framework::Configurable<float> dcabachtopv{"dcabachtopv", .05, "DCA Bach To PV"};
  o2::framework::Configurable<float> cascradius{"cascradius", 0.9, "cascradius"};
  o2::framework::Configurable<float> casccospa{"casccospa", 0.95, "casccospa"};
  o2::framework::Configurable<float> dcacascdau{"dcacascdau", 1.0, "DCA cascade Daughters"};
  o2::framework::Configurable<float> lambdaMassWindow{"lambdaMassWindow", .010, "Distance from Lambda mass (does not apply to KF path)"};
  o2::framework::Configurable<float> maxDaughterEta{"maxDaughterEta", 5.0, "Maximum daughter eta (in abs value)"};

  // KF building specific
  o2::framework::Configurable<bool> kfTuneForOmega{"kfTuneForOmega", false, "if enabled, take main cascade properties from Omega fit instead of Xi fit (= default)"};
  o2::framework::Configurable<int> kfConstructMethod{"kfConstructMethod", 2, "KF Construct Method"};
  o2::framework::Configurable<bool> kfUseV0MassConstraint{"kfUseV0MassConstraint", true, "KF: use Lambda mass constraint"};
  o2::framework::Configurable<bool> kfUseCascadeMassConstraint{"kfUseCascadeMassConstraint", false, "KF: use Cascade mass constraint - WARNING: not adequate for inv mass analysis of Xi"};
  o2::framework::Configurable<bool> kfDoDCAFitterPreMinimV0{"kfDoDCAFitterPreMinimV0", true, "KF: do DCAFitter pre-optimization before KF fit to include material corrections for V0"};
  o2::framework::Configurable<bool> kfDoDCAFitterPreMinimCasc{"kfDoDCAFitterPreMinimCasc", true, "KF: do DCAFitter pre-optimization before KF fit to include material corrections for Xi"};

  // MC builder options
  o2::framework::Configurable<bool> mc_populateCascMCCoresSymmetric{"mc_populateCascMCCoresSymmetric", false, "populate CascMCCores table for derived data analysis, keep CascMCCores joinable with CascCores"};
  o2::framework::Configurable<bool> mc_populateCascMCCoresAsymmetric{"mc_populateCascMCCoresAsymmetric", true, "populate CascMCCores table for derived data analysis, create CascCores -> CascMCCores interlink. Saves only labeled Cascades."};
  o2::framework::Configurable<bool> mc_addGeneratedXiMinus{"mc_addGeneratedXiMinus", true, "add CascMCCore entry for generated, not-recoed XiMinus"};
  o2::framework::Configurable<bool> mc_addGeneratedXiPlus{"mc_addGeneratedXiPlus", true, "add CascMCCore entry for generated, not-recoed XiPlus"};
  o2::framework::Configurable<bool> mc_addGeneratedOmegaMinus{"mc_addGeneratedOmegaMinus", true, "add CascMCCore entry for generated, not-recoed OmegaMinus"};
  o2::framework::Configurable<bool> mc_addGeneratedOmegaPlus{"mc_addGeneratedOmegaPlus", true, "add CascMCCore entry for generated, not-recoed OmegaPlus"};
  o2::framework::Configurable<bool> mc_treatPiToMuDecays{"mc_treatPiToMuDecays", true, "if true, will correctly capture pi -> mu and V0 label will still point to originating V0 decay in those cases. Nota bene: prong info will still be for the muon!"};
  o2::framework::Configurable<float> mc_rapidityWindow{"mc_rapidityWindow", 0.5, "rapidity window to save non-recoed candidates"};
  o2::framework::Configurable<bool> mc_keepOnlyPhysicalPrimary{"mc_keepOnlyPhysicalPrimary", false, "Keep only physical primary generated cascades if not recoed"};
  o2::framework::Configurable<bool> mc_findableDetachedCascade{"mc_findableDetachedCascade", false, "if true, generate findable cascades that have collisionId -1. Caution advised."};
};

// preselection options
struct preSelectOpts : o2::framework::ConfigurableGroup {
  std::string prefix = "preSelectOpts";
  o2::framework::Configurable<bool> preselectOnlyDesiredV0s{"preselectOnlyDesiredV0s", false, "preselect only V0s with compatible TPC PID and mass info"};
  o2::framework::Configurable<bool> preselectOnlyDesiredCascades{"preselectOnlyDesiredCascades", false, "preselect only Cascades with compatible TPC PID and mass info"};

  // lifetime preselection options
  // apply lifetime cuts to V0 and cascade candidates
  // unit of measurement: centimeters
  // lifetime of K0Short ~2.6844 cm, no feeddown and plenty to cut
  // lifetime of Lambda ~7.9 cm but keep in mind feeddown from cascades
  // lifetime of Xi ~4.91 cm
  // lifetime of Omega ~2.461 cm
  o2::framework::Configurable<o2::framework::LabeledArray<float>> lifetimeCut{"lifetimeCut", {defaultLifetimeCuts[0], 4, {"lifetimeCutK0S", "lifetimeCutLambda", "lifetimeCutXi", "lifetimeCutOmega"}}, "Lifetime cut for V0s and cascades (cm)"};

  // mass preselection options
  o2::framework::Configurable<float> massCutPhoton{"massCutPhoton", 0.3, "Photon max mass"};
  o2::framework::Configurable<o2::framework::LabeledArray<float>> massCutK0{"massCutK0", {defaultK0MassWindowParameters[0], 4, {"constant", "linear", "expoConstant", "expoRelax"}}, "mass parameters for K0"};
  o2::framework::Configurable<o2::framework::LabeledArray<float>> massCutLambda{"massCutLambda", {defaultLambdaWindowParameters[0], 4, {"constant", "linear", "expoConstant", "expoRelax"}}, "mass parameters for Lambda"};
  o2::framework::Configurable<o2::framework::LabeledArray<float>> massCutXi{"massCutXi", {defaultXiMassWindowParameters[0], 4, {"constant", "linear", "expoConstant", "expoRelax"}}, "mass parameters for Xi"};
  o2::framework::Configurable<o2::framework::LabeledArray<float>> massCutOm{"massCutOm", {defaultOmMassWindowParameters[0], 4, {"constant", "linear", "expoConstant", "expoRelax"}}, "mass parameters for Omega"};
  o2::framework::Configurable<float> massWindownumberOfSigmas{"massWindownumberOfSigmas", 20, "number of sigmas around mass peaks to keep"};
  o2::framework::Configurable<float> massWindowSafetyMargin{"massWindowSafetyMargin", 0.001, "Extra mass window safety margin (in GeV/c2)"};

  // TPC PID preselection options
  o2::framework::Configurable<float> maxTPCpidNsigma{"maxTPCpidNsigma", 5.0, "Maximum TPC PID N sigma (in abs value)"};
};

class BuilderModule
{
 public:
  BuilderModule()
  {
    // constructor
  }

  // mass windows
  float getMassSigmaK0Short(float pt)
  {
    return preSelectOpts.massCutK0->get("constant") + pt * preSelectOpts.massCutK0->get("linear") + preSelectOpts.massCutK0->get("expoConstant") * TMath::Exp(-pt / preSelectOpts.massCutK0->get("expoRelax"));
  }
  float getMassSigmaLambda(float pt)
  {
    return preSelectOpts.massCutLambda->get("constant") + pt * preSelectOpts.massCutLambda->get("linear") + preSelectOpts.massCutLambda->get("expoConstant") * TMath::Exp(-pt / preSelectOpts.massCutLambda->get("expoRelax"));
  }
  float getMassSigmaXi(float pt)
  {
    return preSelectOpts.massCutXi->get("constant") + pt * preSelectOpts.massCutXi->get("linear") + preSelectOpts.massCutXi->get("expoConstant") * TMath::Exp(-pt / preSelectOpts.massCutXi->get("expoRelax"));
  }
  float getMassSigmaOmega(float pt)
  {
    return preSelectOpts.massCutOm->get("constant") + pt * preSelectOpts.massCutOm->get("linear") + preSelectOpts.massCutOm->get("expoConstant") * TMath::Exp(-pt / preSelectOpts.massCutOm->get("expoRelax"));
  }

  int nEnabledTables = 0;

  // helper object
  o2::pwglf::strangenessBuilderHelper straHelper;

  // for handling TPC-only tracks (photons)
  int mRunNumber;
  o2::aod::common::TPCVDriftManager mVDriftMgr;

  // for establishing CascData/KFData/TraCascData interlinks
  struct {
    std::vector<int> cascCoreToCascades;
    std::vector<int> kfCascCoreToCascades;
    std::vector<int> traCascCoreToCascades;
    std::vector<int> cascadeToCascCores;
    std::vector<int> cascadeToKFCascCores;
    std::vector<int> cascadeToTraCascCores;
  } interlinks;

  //*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
  // struct to add abstraction layer between V0s, Cascades and build indices
  // desirable for adding extra (findable, etc) V0s, Cascades to built list
  struct trackEntry {
    int globalId = -1;
    int originId = -1;
    int mcCollisionId = -1;
    int pdgCode = -1;
  };
  struct v0Entry {
    int globalId = -1;
    int collisionId = -1;
    int posTrackId = -1;
    int negTrackId = -1;
    int v0Type = 0;
    int pdgCode = 0;     // undefined if not MC - useful for faster finding
    int particleId = -1; // de-reference the V0 particle if necessary
    bool isCollinearV0 = false;
    bool found = false;
  };
  struct cascadeEntry {
    int globalId = -1;
    int collisionId = -1;
    int v0Id = -1;
    int posTrackId = -1;
    int negTrackId = -1;
    int bachTrackId = -1;
    bool found = false;
  };

  //*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
  // Helper struct to contain V0MCCore information prior to filling
  struct mcV0info {
    int label = -1;
    int motherLabel = -1;
    int pdgCode = 0;
    int pdgCodeMother = 0;
    int pdgCodePositive = 0;
    int pdgCodeNegative = 0;
    int mcCollision = -1;
    bool isPhysicalPrimary = false;
    int processPositive = -1;
    int processNegative = -1;
    std::array<float, 3> xyz;
    std::array<float, 3> posP;
    std::array<float, 3> negP;
    std::array<float, 3> momentum;
  };
  mcV0info thisInfo;
  //*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
  // Helper struct to contain CascMCCore information prior to filling
  struct mcCascinfo {
    int label;
    int motherLabel;
    int mcCollision;
    int pdgCode;
    int pdgCodeMother;
    int pdgCodeV0;
    int pdgCodePositive;
    int pdgCodeNegative;
    int pdgCodeBachelor;
    bool isPhysicalPrimary;
    int processPositive = -1;
    int processNegative = -1;
    int processBachelor = -1;
    std::array<float, 3> xyz;
    std::array<float, 3> lxyz;
    std::array<float, 3> posP;
    std::array<float, 3> negP;
    std::array<float, 3> bachP;
    std::array<float, 3> momentum;
    int mcParticlePositive;
    int mcParticleNegative;
    int mcParticleBachelor;
  };
  mcCascinfo thisCascInfo;
  //*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*

  std::vector<v0Entry> v0List;
  std::vector<cascadeEntry> cascadeList;
  std::vector<std::size_t> sorted_v0;
  std::vector<std::size_t> sorted_cascade;

  // for tagging V0s used in cascades
  std::vector<o2::pwglf::v0candidate> v0sFromCascades; // Vector of v0 candidates used in cascades
  std::vector<int> ao2dV0toV0List;                     // index to relate v0s -> v0List
  std::vector<int> v0Map;                              // index to relate v0List -> v0sFromCascades

  // declaration of structs here
  // (N.B.: will be invisible to the outside, create your own copies)
  o2::pwglf::strangenessbuilder::coreConfigurables baseOpts;
  o2::pwglf::strangenessbuilder::v0Configurables v0BuilderOpts;
  o2::pwglf::strangenessbuilder::cascadeConfigurables cascadeBuilderOpts;
  o2::pwglf::strangenessbuilder::preSelectOpts preSelectOpts;

  template <typename TBaseConfigurables, typename TV0Configurables, typename TCascadeConfigurables, typename TPreSelOpts, typename THistoRegistry, typename TInitContext>
  void init(TBaseConfigurables const& inputBaseOpts, TV0Configurables const& inputV0BuilderOpts, TCascadeConfigurables const& inputCascadeBuilderOpts, TPreSelOpts const& inputPreSelectOpts, THistoRegistry& histos, TInitContext& context)
  {
    // read in configurations from the task where it's used
    // could be grouped even further, but should work
    baseOpts = inputBaseOpts;
    v0BuilderOpts = inputV0BuilderOpts;
    cascadeBuilderOpts = inputCascadeBuilderOpts;
    preSelectOpts = inputPreSelectOpts;

    baseOpts.mEnabledTables.resize(nTables, 0);

    LOGF(info, "Checking if strangeness building is required");
    auto& workflows = context.services().template get<o2::framework::RunningWorkflowInfo const>();

    nEnabledTables = 0;

    constexpr int kTablesConst = nTables; // silence warning
    TString listOfRequestors[kTablesConst];
    for (int i = 0; i < nTables; i++) {
      int f = baseOpts.enabledTables->get(tableNames[i].c_str(), "enable");
      if (f == 1) {
        baseOpts.mEnabledTables[i] = 1;
        listOfRequestors[i] = "manual enabling";
      }
      if (f == -1) {
        // autodetect this table in other devices
        for (o2::framework::DeviceSpec const& device : workflows.devices) {
          // Step 1: check if this device subscribed to the V0data table
          for (auto const& input : device.inputs) {
            if (o2::framework::DataSpecUtils::partialMatch(input.matcher, o2::header::DataOrigin("AOD"))) {
              auto&& [origin, description, version] = o2::framework::DataSpecUtils::asConcreteDataMatcher(input.matcher);
              std::string tableNameWithVersion = tableNames[i];
              if (version > 0) {
                tableNameWithVersion += Form("_%03d", version);
              }
              if (input.matcher.binding == tableNameWithVersion) {
                LOGF(info, "Device %s has subscribed to %s (version %i)", device.name, tableNames[i], version);
                listOfRequestors[i].Append(Form("%s ", device.name.c_str()));
                baseOpts.mEnabledTables[i] = 1;
                nEnabledTables++;
              }
            }
          }
        }
      }
    }

    if (nEnabledTables == 0) {
      LOGF(info, "Strangeness building not required. Will suppress all functionality, including logs, from this point forward.");
      return;
    }

    // setup bookkeeping histogram
    auto h = histos.template add<TH1>("hTableBuildingStatistics", "hTableBuildingStatistics", o2::framework::kTH1D, {{nTablesConst, -0.5f, static_cast<float>(nTablesConst)}});
    auto h2 = histos.template add<TH1>("hInputStatistics", "hInputStatistics", o2::framework::kTH1D, {{nTablesConst, -0.5f, static_cast<float>(nTablesConst)}});
    h2->SetTitle("Input table sizes");

    if (v0BuilderOpts.generatePhotonCandidates.value == true) {
      auto hDeduplicationStatistics = histos.template add<TH1>("hDeduplicationStatistics", "hDeduplicationStatistics", o2::framework::kTH1D, {{2, -0.5f, 1.5f}});
      hDeduplicationStatistics->GetXaxis()->SetBinLabel(1, "AO2D V0s");
      hDeduplicationStatistics->GetXaxis()->SetBinLabel(2, "Deduplicated V0s");
    }

    if (preSelectOpts.preselectOnlyDesiredV0s.value == true) {
      auto hPreselectionV0s = histos.template add<TH1>("hPreselectionV0s", "hPreselectionV0s", o2::framework::kTH1D, {{16, -0.5f, 15.5f}});
      hPreselectionV0s->GetXaxis()->SetBinLabel(1, "Not preselected");
      hPreselectionV0s->GetXaxis()->SetBinLabel(2, "#gamma");
      hPreselectionV0s->GetXaxis()->SetBinLabel(3, "K^{0}_{S}");
      hPreselectionV0s->GetXaxis()->SetBinLabel(4, "#gamma, K^{0}_{S}");
      hPreselectionV0s->GetXaxis()->SetBinLabel(5, "#Lambda");
      hPreselectionV0s->GetXaxis()->SetBinLabel(6, "#gamma, #Lambda");
      hPreselectionV0s->GetXaxis()->SetBinLabel(7, "K^{0}_{S}, #Lambda");
      hPreselectionV0s->GetXaxis()->SetBinLabel(8, "#gamma, K^{0}_{S}, #Lambda");
      hPreselectionV0s->GetXaxis()->SetBinLabel(9, "#bar{#Lambda}");
      hPreselectionV0s->GetXaxis()->SetBinLabel(10, "#gamma, #bar{#Lambda}");
      hPreselectionV0s->GetXaxis()->SetBinLabel(11, "K^{0}_{S}, #bar{#Lambda}");
      hPreselectionV0s->GetXaxis()->SetBinLabel(12, "#gamma, K^{0}_{S}, #bar{#Lambda}");
      hPreselectionV0s->GetXaxis()->SetBinLabel(13, "#Lambda, #bar{#Lambda}");
      hPreselectionV0s->GetXaxis()->SetBinLabel(14, "#gamma, #Lambda, #bar{#Lambda}");
      hPreselectionV0s->GetXaxis()->SetBinLabel(15, "K^{0}_{S}, #Lambda, #bar{#Lambda}");
      hPreselectionV0s->GetXaxis()->SetBinLabel(16, "#gamma, K^{0}_{S}, #Lambda, #bar{#Lambda}");
    }

    if (preSelectOpts.preselectOnlyDesiredCascades.value == true) {
      auto hPreselectionCascades = histos.template add<TH1>("hPreselectionCascades", "hPreselectionCascades", o2::framework::kTH1D, {{16, -0.5f, 15.5f}});
      hPreselectionCascades->GetXaxis()->SetBinLabel(1, "Not preselected");
      hPreselectionCascades->GetXaxis()->SetBinLabel(2, "#Xi^{-}");
      hPreselectionCascades->GetXaxis()->SetBinLabel(3, "#Xi^{+}");
      hPreselectionCascades->GetXaxis()->SetBinLabel(4, "#Xi^{-}, #Xi^{+}");
      hPreselectionCascades->GetXaxis()->SetBinLabel(5, "#Omega^{-}");
      hPreselectionCascades->GetXaxis()->SetBinLabel(6, "#Xi^{-}, #Omega^{-}");
      hPreselectionCascades->GetXaxis()->SetBinLabel(7, "#Xi^{+}, #Omega^{-}");
      hPreselectionCascades->GetXaxis()->SetBinLabel(8, "#Xi^{-}, #Xi^{+}, #Omega^{-}");
      hPreselectionCascades->GetXaxis()->SetBinLabel(9, "#Omega^{+}");
      hPreselectionCascades->GetXaxis()->SetBinLabel(10, "#Xi^{-}, #Omega^{+}");
      hPreselectionCascades->GetXaxis()->SetBinLabel(11, "#Xi^{+}, #Omega^{+}");
      hPreselectionCascades->GetXaxis()->SetBinLabel(12, "#Xi^{-}, #Xi^{+}, #Omega^{+}");
      hPreselectionCascades->GetXaxis()->SetBinLabel(13, "#Omega^{-}, #Omega^{+}");
      hPreselectionCascades->GetXaxis()->SetBinLabel(14, "#Xi^{-}, #Omega^{-}, #Omega^{+}");
      hPreselectionCascades->GetXaxis()->SetBinLabel(15, "#Xi^{+}, #Omega^{-}, #Omega^{+}");
      hPreselectionCascades->GetXaxis()->SetBinLabel(16, "#Xi^{-}, #Xi^{+}, #Omega^{-}, #Omega^{+}");
    }

    if (baseOpts.mc_findableMode.value > 0) {
      // save statistics of findable candidate processing
      auto hFindable = histos.template add<TH1>("hFindableStatistics", "hFindableStatistics", o2::framework::kTH1D, {{6, -0.5f, 5.5f}});
      hFindable->SetTitle(Form("Findable mode: %i", static_cast<int>(baseOpts.mc_findableMode.value)));
      hFindable->GetXaxis()->SetBinLabel(1, "AO2D V0s");
      hFindable->GetXaxis()->SetBinLabel(2, "V0s to be built");
      hFindable->GetXaxis()->SetBinLabel(3, "V0s with collId -1");
      hFindable->GetXaxis()->SetBinLabel(4, "AO2D Cascades");
      hFindable->GetXaxis()->SetBinLabel(5, "Cascades to be built");
      hFindable->GetXaxis()->SetBinLabel(6, "Cascades with collId -1");
    }

    auto hPrimaryV0s = histos.template add<TH1>("hPrimaryV0s", "hPrimaryV0s", o2::framework::kTH1D, {{2, -0.5f, 1.5f}});
    hPrimaryV0s->GetXaxis()->SetBinLabel(1, "All V0s");
    hPrimaryV0s->GetXaxis()->SetBinLabel(2, "Primary V0s");

    mRunNumber = 0;

    for (int i = 0; i < nTables; i++) {
      // adjust bookkeeping histogram
      h->GetXaxis()->SetBinLabel(i + 1, tableNames[i].c_str());
      h2->GetXaxis()->SetBinLabel(i + 1, tableNames[i].c_str());
      if (baseOpts.mEnabledTables[i] == false) {
        h->SetBinContent(i + 1, -1); // mark disabled tables, distinguish from zero counts
      }
    }

    LOGF(info, "*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*");
    LOGF(info, " Strangeness builder: basic configuration listing");
    LOGF(info, "*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*");

    if (baseOpts.mc_findableMode.value == 1) {
      LOGF(info, " ===> findable mode 1 is enabled: complement reco-ed with non-found findable");
    }
    if (baseOpts.mc_findableMode.value == 2) {
      LOGF(info, " ===> findable mode 2 is enabled: re-generate all findable from scratch");
    }

    // list enabled tables

    for (int i = 0; i < nTables; i++) {
      // printout to be improved in the future
      if (baseOpts.mEnabledTables[i]) {
        LOGF(info, " -~> Table enabled: %s, requested by %s", tableNames[i], listOfRequestors[i].Data());
        h->SetBinContent(i + 1, 0); // mark enabled
      }
    }
    LOGF(info, "*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*");
    // print base cuts
    LOGF(info, "-~> V0 | min crossed rows ..............: %i", v0BuilderOpts.minCrossedRows.value);
    LOGF(info, "-~> V0 | DCA pos track to PV ...........: %f", v0BuilderOpts.dcapostopv.value);
    LOGF(info, "-~> V0 | DCA neg track to PV ...........: %f", v0BuilderOpts.dcanegtopv.value);
    LOGF(info, "-~> V0 | V0 cosine of PA ...............: %f", v0BuilderOpts.v0cospa.value);
    LOGF(info, "-~> V0 | DCA between V0 daughters ......: %f", v0BuilderOpts.dcav0dau.value);
    LOGF(info, "-~> V0 | V0 2D decay radius ............: %f", v0BuilderOpts.v0radius.value);
    LOGF(info, "-~> V0 | Maximum daughter eta ..........: %f", v0BuilderOpts.maxDaughterEta.value);

    LOGF(info, "-~> Cascade | min crossed rows .........: %i", cascadeBuilderOpts.minCrossedRows.value);
    LOGF(info, "-~> Cascade | DCA bach track to PV .....: %f", cascadeBuilderOpts.dcabachtopv.value);
    LOGF(info, "-~> Cascade | Cascade cosine of PA .....: %f", cascadeBuilderOpts.casccospa.value);
    LOGF(info, "-~> Cascade | Cascade daughter DCA .....: %f", cascadeBuilderOpts.dcacascdau.value);
    LOGF(info, "-~> Cascade | Cascade radius ...........: %f", cascadeBuilderOpts.cascradius.value);
    LOGF(info, "-~> Cascade | Lambda mass window .......: %f", cascadeBuilderOpts.lambdaMassWindow.value);
    LOGF(info, "-~> Cascade | Maximum daughter eta .....: %f", cascadeBuilderOpts.maxDaughterEta.value);
    LOGF(info, "*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*");

    // set V0 parameters in the helper
    straHelper.v0selections.minCrossedRows = v0BuilderOpts.minCrossedRows;
    straHelper.v0selections.dcanegtopv = v0BuilderOpts.dcanegtopv;
    straHelper.v0selections.dcapostopv = v0BuilderOpts.dcapostopv;
    straHelper.v0selections.v0cospa = v0BuilderOpts.v0cospa;
    straHelper.v0selections.dcav0dau = v0BuilderOpts.dcav0dau;
    straHelper.v0selections.v0radius = v0BuilderOpts.v0radius;
    straHelper.v0selections.maxDaughterEta = v0BuilderOpts.maxDaughterEta;

    // set cascade parameters in the helper
    straHelper.cascadeselections.minCrossedRows = cascadeBuilderOpts.minCrossedRows;
    straHelper.cascadeselections.dcabachtopv = cascadeBuilderOpts.dcabachtopv;
    straHelper.cascadeselections.cascradius = cascadeBuilderOpts.cascradius;
    straHelper.cascadeselections.casccospa = cascadeBuilderOpts.casccospa;
    straHelper.cascadeselections.dcacascdau = cascadeBuilderOpts.dcacascdau;
    straHelper.cascadeselections.lambdaMassWindow = cascadeBuilderOpts.lambdaMassWindow;
    straHelper.cascadeselections.maxDaughterEta = cascadeBuilderOpts.maxDaughterEta;

    // Set option to refit with material corrections
    straHelper.fitter.setRefitWithMatCorr(baseOpts.refitWithMaterialCorrection.value);
  }

  // for sorting
  template <typename T>
  std::vector<std::size_t> sort_indices(const std::vector<T>& v, bool doSorting = false)
  {
    std::vector<std::size_t> idx(v.size());
    std::iota(idx.begin(), idx.end(), 0);
    if (doSorting) {
      // do sorting only if requested (not always necessary)
      std::stable_sort(idx.begin(), idx.end(),
                       [&v](std::size_t i1, std::size_t i2) { return v[i1].collisionId < v[i2].collisionId; });
    }
    return idx;
  }

  template <typename TCollisions, typename TCCDB>
  bool initCCDB(TCCDB& ccdb, aod::BCsWithTimestamps const& bcs, TCollisions const& collisions)
  {
    auto bc = collisions.size() ? collisions.begin().template bc_as<aod::BCsWithTimestamps>() : bcs.begin();
    if (!bcs.size()) {
      LOGF(warn, "No BC found, skipping this DF.");
      return false; // signal to skip this DF
    }

    if (mRunNumber == bc.runNumber()) {
      return true;
    }

    auto timestamp = bc.timestamp();

    // Fetch magnetic field from ccdb for current collision
    auto magneticField = o2::base::Propagator::Instance()->getNominalBz();
    LOG(info) << "Configuring for timestamp " << timestamp << " with magnetic field of " << magneticField << " kG";

    // Set magnetic field value once known
    straHelper.fitter.setBz(magneticField);

    LOG(info) << "Fully configured for run: " << bc.runNumber();
    // mmark this run as configured
    mRunNumber = bc.runNumber();

    if (v0BuilderOpts.generatePhotonCandidates.value && v0BuilderOpts.moveTPCOnlyTracks.value) {
      // initialize only if needed, avoid unnecessary CCDB calls
      mVDriftMgr.init(&ccdb->instance());
      mVDriftMgr.update(timestamp);
    }

    return true;
  }

  //__________________________________________________
  void resetInterlinks()
  {
    interlinks.cascCoreToCascades.clear();
    interlinks.kfCascCoreToCascades.clear();
    interlinks.traCascCoreToCascades.clear();
    interlinks.cascadeToCascCores.clear();
    interlinks.cascadeToKFCascCores.clear();
    interlinks.cascadeToTraCascCores.clear();
  }

  //__________________________________________________
  void populateCascadeInterlinks()
  {
    // if (mEnabledTables[kCascToKFRefs]) {
    //   for (const auto& cascCore : interlinks.cascCoreToCascades) {
    //     cascToKFRefs(interlinks.cascadeToKFCascCores[cascCore]);
    //     histos.fill(HIST("hTableBuildingStatistics"), kCascToKFRefs);
    //   }
    // }
    // if (mEnabledTables[kCascToTraRefs]) {
    //   for (const auto& cascCore : interlinks.cascCoreToCascades) {
    //     cascToTraRefs(interlinks.cascadeToTraCascCores[cascCore]);
    //     histos.fill(HIST("hTableBuildingStatistics"), kCascToTraRefs);
    //   }
    // }
    // if (mEnabledTables[kKFToCascRefs]) {
    //   for (const auto& kfCascCore : interlinks.kfCascCoreToCascades) {
    //     kfToCascRefs(interlinks.cascadeToCascCores[kfCascCore]);
    //     histos.fill(HIST("hTableBuildingStatistics"), kKFToCascRefs);
    //   }
    // }
    // if (mEnabledTables[kTraToCascRefs]) {
    //   for (const auto& traCascCore : interlinks.traCascCoreToCascades) {
    //     traToCascRefs(interlinks.cascadeToCascCores[traCascCore]);
    //     histos.fill(HIST("hTableBuildingStatistics"), kTraToCascRefs);
    //   }
    // }
  }

  //__________________________________________________
  template <class TBCs, typename THistoRegistry, typename TCollisions, typename TMCCollisions, typename TV0s, typename TCascades, typename TTracks, typename TMCParticles>
  void prepareBuildingLists(THistoRegistry& histos, TCollisions const& collisions, TMCCollisions const& mcCollisions, TV0s const& v0s, TCascades const& cascades, TTracks const& tracks, TMCParticles const& mcParticles)
  {
    // this function prepares the v0List and cascadeList depending on
    // how the task has been set up. Standard operation simply uses
    // the existing V0s and Cascades from AO2D, while findable MC
    // operation either complements with all findable-but-not-found
    // or resets and fills with all findable.
    //
    // Whenever using findable candidates, they will be appropriately
    // marked for posterior analysis using 'tag' variables.
    //
    // findable mode legend:
    // 0: simple passthrough of V0s, Cascades in AO2Ds
    //    (in data, this is the only mode possible!)
    // 1: add extra findable that haven't been found
    // 2: generate only findable (no background)

    // redo lists from scratch
    v0List.clear();
    cascadeList.clear();
    sorted_v0.clear();
    sorted_cascade.clear();
    ao2dV0toV0List.clear();

    trackEntry currentTrackEntry;
    v0Entry currentV0Entry;
    cascadeEntry currentCascadeEntry;

    std::vector<int> bestCollisionArray;          // stores McCollision -> Collision map
    std::vector<int> bestCollisionNContribsArray; // stores Ncontribs for biggest coll assoc to mccoll

    int collisionLessV0s = 0;
    int collisionLessCascades = 0;

    if (baseOpts.mc_findableMode.value > 0) {
      if constexpr (soa::is_table<TMCCollisions>) {
        // if mcCollisions exist, assemble mcColl -> bestRecoColl map here
        bestCollisionArray.clear();
        bestCollisionNContribsArray.clear();
        bestCollisionArray.resize(mcCollisions.size(), -1);          // marks not reconstructed
        bestCollisionNContribsArray.resize(mcCollisions.size(), -1); // marks not reconstructed

        // single loop over double loop at a small cost in memory for extra array
        for (const auto& collision : collisions) {
          if (collision.has_mcCollision()) {
            if (collision.numContrib() > bestCollisionNContribsArray[collision.mcCollisionId()]) {
              bestCollisionArray[collision.mcCollisionId()] = collision.globalIndex();
              bestCollisionNContribsArray[collision.mcCollisionId()] = collision.numContrib();
            }
          }
        } // end collision loop
      } // end is_table<TMCCollisions>
    } // end findable mode check

    if (baseOpts.mc_findableMode.value < 2) {
      // keep all unless de-duplication active
      ao2dV0toV0List.resize(v0s.size(), -1); // -1 means keep, -2 means do not keep

      if (baseOpts.deduplicationAlgorithm > 0 && v0BuilderOpts.generatePhotonCandidates) {
        // handle duplicates explicitly: group V0s according to (p,n) indices
        // will provide a list of collisionIds (in V0group), allowing for
        // easy de-duplication when passing to the v0List
        std::vector<o2::pwglf::V0group> v0tableGrouped = o2::pwglf::groupDuplicates(v0s);
        histos.fill(HIST("hDeduplicationStatistics"), 0.0, v0s.size());
        histos.fill(HIST("hDeduplicationStatistics"), 1.0, v0tableGrouped.size());

        // process grouped duplicates, remove 'bad' ones
        for (size_t iV0 = 0; iV0 < v0tableGrouped.size(); iV0++) {
          auto pTrack = tracks.rawIteratorAt(v0tableGrouped[iV0].posTrackId);
          auto nTrack = tracks.rawIteratorAt(v0tableGrouped[iV0].negTrackId);

          bool isPosTPCOnly = (pTrack.hasTPC() && !pTrack.hasITS() && !pTrack.hasTRD() && !pTrack.hasTOF());
          bool isNegTPCOnly = (nTrack.hasTPC() && !nTrack.hasITS() && !nTrack.hasTRD() && !nTrack.hasTOF());

          // skip single copy V0s
          if (v0tableGrouped[iV0].collisionIds.size() == 1) {
            continue;
          }

          // don't try to de-duplicate if no track is TPC only
          if (!isPosTPCOnly && !isNegTPCOnly) {
            continue;
          }

          // fitness criteria defined here
          float bestPointingAngle = 10; // a nonsense angle, anything's better
          size_t bestPointingAngleIndex = -1;

          float bestDCADaughters = 1e+3; // an excessively large DCA
          size_t bestDCADaughtersIndex = -1;

          for (size_t ic = 0; ic < v0tableGrouped[iV0].collisionIds.size(); ic++) {
            // get track parametrizations, collisions
            auto posTrackPar = getTrackParCov(pTrack);
            auto negTrackPar = getTrackParCov(nTrack);
            auto const& collision = collisions.rawIteratorAt(v0tableGrouped[iV0].collisionIds[ic]);

            // handle TPC-only tracks properly (photon conversions)
            if (v0BuilderOpts.moveTPCOnlyTracks) {
              if (isPosTPCOnly) {
                // Nota bene: positive is TPC-only -> this entire V0 merits treatment as photon candidate
                posTrackPar.setPID(o2::track::PID::Electron);
                negTrackPar.setPID(o2::track::PID::Electron);
                if (!mVDriftMgr.moveTPCTrack<TBCs, TCollisions>(collision, pTrack, posTrackPar)) {
                  continue;
                }
              }
              if (isNegTPCOnly) {
                // Nota bene: negative is TPC-only -> this entire V0 merits treatment as photon candidate
                posTrackPar.setPID(o2::track::PID::Electron);
                negTrackPar.setPID(o2::track::PID::Electron);
                if (!mVDriftMgr.moveTPCTrack<TBCs, TCollisions>(collision, nTrack, negTrackPar)) {
                  continue;
                }
              }
            } // end TPC drift treatment

            // process candidate with helper, generate properties for consulting
            // first 'false' : do not apply selections: do as much as possible to preserve
            // second 'false': do not calculate prong DCA to PV, unnecessary, costly if XIU = 83.1f
            // candidate at this level and do not select with topo selections
            if (straHelper.buildV0Candidate<false, false>(v0tableGrouped[iV0].collisionIds[ic], collision.posX(), collision.posY(), collision.posZ(), pTrack, nTrack, posTrackPar, negTrackPar, true, false, true)) {
              // candidate built, check pointing angle
              if (straHelper.v0.pointingAngle < bestPointingAngle) {
                bestPointingAngle = straHelper.v0.pointingAngle;
                bestPointingAngleIndex = ic;
              }
              if (straHelper.v0.daughterDCA < bestDCADaughters) {
                bestDCADaughters = straHelper.v0.daughterDCA;
                bestDCADaughtersIndex = ic;
              }
            } // end build V0
          } // end candidate loop

          // mark de-duplicated candidates
          for (size_t ic = 0; ic < v0tableGrouped[iV0].collisionIds.size(); ic++) {
            ao2dV0toV0List[v0tableGrouped[iV0].V0Ids[ic]] = -2;
            // algorithm 1: best pointing angle
            if (bestPointingAngleIndex == ic && baseOpts.deduplicationAlgorithm.value == 1) {
              ao2dV0toV0List[v0tableGrouped[iV0].V0Ids[ic]] = -1; // keep best only
            }
            if (bestDCADaughtersIndex == ic && baseOpts.deduplicationAlgorithm.value == 2) {
              ao2dV0toV0List[v0tableGrouped[iV0].V0Ids[ic]] = -1; // keep best only
            }
            if (bestDCADaughtersIndex == ic && bestPointingAngleIndex == ic && baseOpts.deduplicationAlgorithm.value == 3) {
              ao2dV0toV0List[v0tableGrouped[iV0].V0Ids[ic]] = -1; // keep best only
            }
          }
        } // end V0 loop
      } // end de-duplication process

      for (const auto& v0 : v0s) {
        if (ao2dV0toV0List[v0.globalIndex()] == -1) {       // keep only de-duplicated
          ao2dV0toV0List[v0.globalIndex()] = v0List.size(); // maps V0s to the corresponding v0List entry
          currentV0Entry.globalId = v0.globalIndex();
          currentV0Entry.collisionId = v0.collisionId();
          currentV0Entry.posTrackId = v0.posTrackId();
          currentV0Entry.negTrackId = v0.negTrackId();
          currentV0Entry.v0Type = v0.v0Type();
          currentV0Entry.pdgCode = 0;
          currentV0Entry.particleId = -1;
          currentV0Entry.isCollinearV0 = v0.isCollinearV0();
          currentV0Entry.found = true;
          v0List.push_back(currentV0Entry);
        }
      }
    }
    // any mode other than 0 will require mcParticles
    if constexpr (soa::is_table<TMCCollisions>) {
      if (baseOpts.mc_findableMode.value > 0) {
        // for search if existing or not
        int v0ListReconstructedSize = v0List.size();

        // find extra candidates, step 1: find subset of tracks that interest
        std::vector<trackEntry> positiveTrackArray;
        std::vector<trackEntry> negativeTrackArray;
        // vector elements: track index, origin index [, mc collision id, pdg code]
        int dummy = -1; // unnecessary in this path
        for (const auto& track : tracks) {
          if (!track.has_mcParticle()) {
            continue; // skip this, it's trouble
          }
          auto particle = track.template mcParticle_as<TMCParticles>();
          int originParticleIndex = getOriginatingParticle(particle, dummy, v0BuilderOpts.mc_treatPiToMuDecays);
          if (originParticleIndex < 0) {
            continue; // skip this, it's trouble (2)
          }
          auto originParticle = mcParticles.rawIteratorAt(originParticleIndex);

          bool trackIsInteresting = false;
          if (
            (originParticle.pdgCode() == PDG_t::kK0Short && v0BuilderOpts.mc_addGeneratedK0Short.value > 0) ||
            (originParticle.pdgCode() == PDG_t::kLambda0 && v0BuilderOpts.mc_addGeneratedLambda.value > 0) ||
            (originParticle.pdgCode() == PDG_t::kLambda0Bar && v0BuilderOpts.mc_addGeneratedAntiLambda.value > 0) ||
            (originParticle.pdgCode() == PDG_t::kGamma && v0BuilderOpts.mc_addGeneratedGamma.value > 0)) {
            trackIsInteresting = true;
          }
          if (!trackIsInteresting) {
            continue; // skip this, it's uninteresting
          }

          currentTrackEntry.globalId = static_cast<int>(track.globalIndex());
          currentTrackEntry.originId = originParticleIndex;
          currentTrackEntry.mcCollisionId = originParticle.mcCollisionId();
          currentTrackEntry.pdgCode = originParticle.pdgCode();

          // now separate according to particle species
          if (track.sign() < 0) {
            negativeTrackArray.push_back(currentTrackEntry);
          } else {
            positiveTrackArray.push_back(currentTrackEntry);
          }
        }

        // Nested loop only with valuable tracks
        for (const auto& positiveTrackIndex : positiveTrackArray) {
          for (const auto& negativeTrackIndex : negativeTrackArray) {
            if (positiveTrackIndex.originId != negativeTrackIndex.originId) {
              continue; // not the same originating particle
            }
            // findable mode 1: add non-reconstructed as v0Type 8
            if (baseOpts.mc_findableMode.value == 1) {
              bool detected = false;
              for (int ii = 0; ii < v0ListReconstructedSize; ii++) {
                // check if this particular combination already exists in v0List
                if (v0List[ii].posTrackId == positiveTrackIndex.globalId &&
                    v0List[ii].negTrackId == negativeTrackIndex.globalId) {
                  detected = true;
                  // override pdg code with something useful for cascade findable math
                  v0List[ii].pdgCode = positiveTrackIndex.pdgCode;
                  break;
                }
              }
              if (detected == false) {
                // collision index: from best-version-of-this-mcCollision
                // nota bene: this could be negative, caution advised
                currentV0Entry.globalId = -1;
                currentV0Entry.collisionId = bestCollisionArray[positiveTrackIndex.mcCollisionId];
                currentV0Entry.posTrackId = positiveTrackIndex.globalId;
                currentV0Entry.negTrackId = negativeTrackIndex.globalId;
                currentV0Entry.v0Type = 1;
                currentV0Entry.pdgCode = positiveTrackIndex.pdgCode;
                currentV0Entry.particleId = positiveTrackIndex.originId;
                currentV0Entry.isCollinearV0 = false;
                if (v0BuilderOpts.mc_addGeneratedGammaMakeCollinear.value && currentV0Entry.pdgCode == PDG_t::kGamma) {
                  currentV0Entry.isCollinearV0 = true;
                }
                currentV0Entry.found = false;
                if (bestCollisionArray[positiveTrackIndex.mcCollisionId] < 0) {
                  collisionLessV0s++;
                }
                if (v0BuilderOpts.mc_findableDetachedV0.value || currentV0Entry.collisionId >= 0) {
                  v0List.push_back(currentV0Entry);
                }
              }
            }
            // findable mode 2
            if (baseOpts.mc_findableMode.value == 2) {
              currentV0Entry.globalId = -1;
              currentV0Entry.collisionId = bestCollisionArray[positiveTrackIndex.mcCollisionId];
              currentV0Entry.posTrackId = positiveTrackIndex.globalId;
              currentV0Entry.negTrackId = negativeTrackIndex.globalId;
              currentV0Entry.v0Type = 1;
              currentV0Entry.pdgCode = positiveTrackIndex.pdgCode;
              currentV0Entry.particleId = positiveTrackIndex.originId;
              currentV0Entry.isCollinearV0 = false;
              if (v0BuilderOpts.mc_addGeneratedGammaMakeCollinear.value && currentV0Entry.pdgCode == PDG_t::kGamma) {
                currentV0Entry.isCollinearV0 = true;
              }
              currentV0Entry.found = false;
              for (const auto& v0 : v0s) {
                if (v0.posTrackId() == positiveTrackIndex.globalId &&
                    v0.negTrackId() == negativeTrackIndex.globalId) {
                  // this will override type, but not collision index
                  // N.B.: collision index checks still desirable!
                  currentV0Entry.globalId = v0.globalIndex();
                  currentV0Entry.v0Type = v0.v0Type();
                  currentV0Entry.isCollinearV0 = v0.isCollinearV0();
                  currentV0Entry.found = true;
                  break;
                }
              }
              if (v0BuilderOpts.mc_findableDetachedV0.value || currentV0Entry.collisionId >= 0) {
                v0List.push_back(currentV0Entry);
              }
            }
          }
        } // end positive / negative track loops

        // fill findable statistics table
        histos.fill(HIST("hFindableStatistics"), 0.0, v0s.size());
        histos.fill(HIST("hFindableStatistics"), 1.0, v0List.size());
        histos.fill(HIST("hFindableStatistics"), 2.0, collisionLessV0s);

      } // end findableMode > 0 check
    } // end soa::is_table<TMCCollisions>

    // determine properly collision-id-sorted index array for later use
    // N.B.: necessary also before cascade part
    sorted_v0.clear();
    sorted_v0 = sort_indices(v0List, (baseOpts.mc_findableMode.value > 0));

    // Cascade part if cores are requested, skip otherwise
    if (baseOpts.mEnabledTables[kStoredCascCores] || baseOpts.mEnabledTables[kStoredKFCascCores]) {
      if (baseOpts.mc_findableMode.value < 2) {
        // simple passthrough: copy existing cascades to build list
        for (const auto& cascade : cascades) {
          auto const& v0 = cascade.v0();
          currentCascadeEntry.globalId = cascade.globalIndex();
          currentCascadeEntry.collisionId = cascade.collisionId();
          currentCascadeEntry.v0Id = ao2dV0toV0List[v0.globalIndex()];
          currentCascadeEntry.posTrackId = v0.posTrackId();
          currentCascadeEntry.negTrackId = v0.negTrackId();
          currentCascadeEntry.bachTrackId = cascade.bachelorId();
          currentCascadeEntry.found = true;
          cascadeList.push_back(currentCascadeEntry);
        }
      }

      // any mode other than 0 will require mcParticles
      if constexpr (soa::is_table<TMCCollisions>) {
        if (baseOpts.mc_findableMode.value > 0) {
          // for search if existing or not
          size_t cascadeListReconstructedSize = cascadeList.size();

          // determine which tracks are of interest
          std::vector<trackEntry> bachelorTrackArray;
          // vector elements: track index, origin index, mc collision id, pdg code]
          int dummy = -1; // unnecessary in this path
          for (const auto& track : tracks) {
            if (!track.has_mcParticle()) {
              continue; // skip this, it's trouble
            }
            auto particle = track.template mcParticle_as<TMCParticles>();
            int originParticleIndex = getOriginatingParticle(particle, dummy, cascadeBuilderOpts.mc_treatPiToMuDecays);
            if (originParticleIndex < 0) {
              continue; // skip this, it's trouble (2)
            }
            auto originParticle = mcParticles.rawIteratorAt(originParticleIndex);

            bool trackIsInteresting = false;
            if (
              (originParticle.pdgCode() == PDG_t::kXiMinus && cascadeBuilderOpts.mc_addGeneratedXiMinus.value > 0) ||
              (originParticle.pdgCode() == PDG_t::kXiPlusBar && cascadeBuilderOpts.mc_addGeneratedXiPlus.value > 0) ||
              (originParticle.pdgCode() == PDG_t::kOmegaMinus && cascadeBuilderOpts.mc_addGeneratedOmegaMinus.value > 0) ||
              (originParticle.pdgCode() == PDG_t::kOmegaPlusBar && cascadeBuilderOpts.mc_addGeneratedOmegaPlus.value > 0)) {
              trackIsInteresting = true;
            }
            if (!trackIsInteresting) {
              continue; // skip this, it's uninteresting
            }

            currentTrackEntry.globalId = static_cast<int>(track.globalIndex());
            currentTrackEntry.originId = originParticleIndex;
            currentTrackEntry.mcCollisionId = originParticle.mcCollisionId();
            currentTrackEntry.pdgCode = originParticle.pdgCode();

            // populate list of bachelor tracks to pair
            bachelorTrackArray.push_back(currentTrackEntry);
          }

          // determine which V0s are of interest to pair and do pairing
          for (size_t v0i = 0; v0i < v0List.size(); v0i++) {
            auto v0 = v0List[sorted_v0[v0i]];

            if (std::abs(v0.pdgCode) != PDG_t::kLambda0) {
              continue; // this V0 isn't a lambda, can't come from a cascade: skip
            }
            if (v0.particleId < 0) {
              continue; // no de-referencing possible (e.g. background, ...)
            }
            auto v0Particle = mcParticles.rawIteratorAt(v0.particleId);

            int v0OriginParticleIndex = -1;
            if (v0Particle.has_mothers()) {
              auto const& motherList = v0Particle.template mothers_as<TMCParticles>();
              if (motherList.size() == 1) {
                for (const auto& mother : motherList) {
                  v0OriginParticleIndex = mother.globalIndex();
                }
              }
            }
            if (v0OriginParticleIndex < 0) {
              continue;
            }
            auto v0OriginParticle = mcParticles.rawIteratorAt(v0OriginParticleIndex);

            if (std::abs(v0OriginParticle.pdgCode()) != PDG_t::kXiMinus && std::abs(v0OriginParticle.pdgCode()) != PDG_t::kOmegaMinus) {
              continue; // this V0 does not come from any particle of interest, don't try
            }
            for (const auto& bachelorTrackIndex : bachelorTrackArray) {
              if (bachelorTrackIndex.originId != v0OriginParticle.globalIndex()) {
                continue;
              }
              // if we are here: v0 origin is 3312 or 3334, bachelor origin matches V0 origin
              // findable mode 1: add non-reconstructed as cascadeType 1
              if (baseOpts.mc_findableMode.value == 1) {
                bool detected = false;
                for (size_t ii = 0; ii < cascadeListReconstructedSize; ii++) {
                  // check if this particular combination already exists in cascadeList
                  // caution: use track indices (immutable) but not V0 indices (re-indexing)
                  if (cascadeList[ii].posTrackId == v0.posTrackId &&
                      cascadeList[ii].negTrackId == v0.negTrackId &&
                      cascadeList[ii].bachTrackId == bachelorTrackIndex.globalId) {
                    detected = true;
                    break;
                  }
                }
                if (detected == false) {
                  // collision index: from best-version-of-this-mcCollision
                  // nota bene: this could be negative, caution advised
                  currentCascadeEntry.globalId = -1;
                  currentCascadeEntry.collisionId = bestCollisionArray[bachelorTrackIndex.mcCollisionId];
                  currentCascadeEntry.v0Id = v0i; // correct information here
                  currentCascadeEntry.posTrackId = v0.posTrackId;
                  currentCascadeEntry.negTrackId = v0.negTrackId;
                  currentCascadeEntry.bachTrackId = bachelorTrackIndex.globalId;
                  currentCascadeEntry.found = false;
                  cascadeList.push_back(currentCascadeEntry);
                  if (bestCollisionArray[bachelorTrackIndex.mcCollisionId] < 0) {
                    collisionLessCascades++;
                  }
                  if (cascadeBuilderOpts.mc_findableDetachedCascade.value || currentCascadeEntry.collisionId >= 0) {
                    cascadeList.push_back(currentCascadeEntry);
                  }
                }
              }

              // findable mode 2: determine type based on cascade table,
              // with type 1 being reserved to findable-but-not-found
              if (baseOpts.mc_findableMode.value == 2) {
                currentCascadeEntry.globalId = -1;
                currentCascadeEntry.collisionId = bestCollisionArray[bachelorTrackIndex.mcCollisionId];
                currentCascadeEntry.v0Id = v0i; // fill this in one go later
                currentCascadeEntry.posTrackId = v0.posTrackId;
                currentCascadeEntry.negTrackId = v0.negTrackId;
                currentCascadeEntry.bachTrackId = bachelorTrackIndex.globalId;
                currentCascadeEntry.found = false;
                if (bestCollisionArray[bachelorTrackIndex.mcCollisionId] < 0) {
                  collisionLessCascades++;
                }
                for (const auto& cascade : cascades) {
                  auto const& v0fromAOD = cascade.v0();
                  if (v0fromAOD.posTrackId() == v0.posTrackId &&
                      v0fromAOD.negTrackId() == v0.negTrackId &&
                      cascade.bachelorId() == bachelorTrackIndex.globalId) {
                    // this will override type, but not collision index
                    // N.B.: collision index checks still desirable!
                    currentCascadeEntry.found = true;
                    currentCascadeEntry.globalId = cascade.globalIndex();
                    break;
                  }
                }
                if (cascadeBuilderOpts.mc_findableDetachedCascade.value || currentCascadeEntry.collisionId >= 0) {
                  cascadeList.push_back(currentCascadeEntry);
                }
              }
            } // end bachelorTrackArray loop
          } // end v0List loop

          // at this stage, cascadeList is alright, but the v0 indices are still not
          // correct. We'll have to loop over all V0s and find the appropriate matches
          // ---> but only in mode 1, and only for AO2D-native V0s
          if (baseOpts.mc_findableMode.value == 1) {
            for (size_t casci = 0; casci < cascadeListReconstructedSize; casci++) {
              // loop over v0List to find corresponding v0 index, but do it in sorted way
              for (size_t v0i = 0; v0i < v0List.size(); v0i++) {
                auto v0 = v0List[sorted_v0[v0i]];
                if (cascadeList[casci].posTrackId == v0.posTrackId &&
                    cascadeList[casci].negTrackId == v0.negTrackId) {
                  cascadeList[casci].v0Id = v0i; // fix, point to correct V0 index
                  break;
                }
              }
            }
          }
          // we should now be done! collect statistics
          histos.fill(HIST("hFindableStatistics"), 3.0, cascades.size());
          histos.fill(HIST("hFindableStatistics"), 4.0, cascadeList.size());
          histos.fill(HIST("hFindableStatistics"), 5.0, collisionLessCascades);

        } // end findable mode check
      } // end soa::is_table<TMCCollisions>

      // we need to allow for sorted use of cascadeList
      sorted_cascade.clear();
      sorted_cascade = sort_indices(cascadeList, (baseOpts.mc_findableMode.value > 0));
    }

    LOGF(debug, "AO2D input: %i V0s, %i cascades. Building list sizes: %i V0s, %i cascades", v0s.size(), cascades.size(), v0List.size(), cascadeList.size());
  }

  //__________________________________________________
  template <typename TV0s, typename TCascades, typename TTrackedCascades>
  void markV0sUsedInCascades(TV0s const& v0s, TCascades const& cascades, TTrackedCascades const& trackedCascades)
  {
    int v0sUsedInCascades = 0;
    v0sFromCascades.clear();
    v0Map.clear();
    v0Map.resize(v0List.size(), -2); // marks not used
    if (baseOpts.useV0BufferForCascades.value == false) {
      return; // don't attempt to mark needlessly
    }
    if (baseOpts.mEnabledTables[kStoredCascCores]) {
      for (const auto& cascade : cascadeList) {
        if (cascade.v0Id < 0)
          continue;
        if (v0Map[cascade.v0Id] == -2) {
          v0sUsedInCascades++;
        }
        v0Map[cascade.v0Id] = -1; // marks used (but isn't the index of a properly built V0, which would be >= 0)
      }
    }
    int trackedCascadeCount = 0;
    if constexpr (soa::is_table<TTrackedCascades>) {
      // tracked only created outside of findable mode
      if (baseOpts.mEnabledTables[kStoredTraCascCores] && baseOpts.mc_findableMode.value == 0) {
        trackedCascadeCount = trackedCascades.size();
        for (const auto& trackedCascade : trackedCascades) {
          auto const& cascade = trackedCascade.cascade();
          if (v0Map[ao2dV0toV0List[cascade.v0Id()]] == -2) {
            v0sUsedInCascades++;
          }
          v0Map[ao2dV0toV0List[cascade.v0Id()]] = -1; // marks used (but isn't the index of a built V0, which would be >= 0)
        }
      }
    }
    LOGF(debug, "V0 total %i, Cascade total %i, Tracked cascade total %i, V0s flagged used in cascades: %i", v0s.size(), cascades.size(), trackedCascadeCount, v0sUsedInCascades);
  }

  //__________________________________________________
  template <class TBCs, typename THistoRegistry, typename TCollisions, typename TTracks, typename TV0s, typename TMCParticles, typename TProducts>
  void buildV0s(THistoRegistry& histos, TCollisions const& collisions, TV0s const& v0s, TTracks const& tracks, TMCParticles const& mcParticles, TProducts& products)
  {
    // prepare MC containers (not necessarily used)
    std::vector<mcV0info> mcV0infos; // V0MCCore information
    std::vector<bool> mcParticleIsReco;

    if constexpr (soa::is_table<TMCParticles>) {
      // do this if provided with a mcParticle table as well
      mcParticleIsReco.resize(mcParticles.size(), false);
    }

    int nV0s = 0;
    // Loops over all V0s in the time frame
    histos.fill(HIST("hInputStatistics"), kV0CoresBase, v0s.size());
    for (size_t iv0 = 0; iv0 < v0List.size(); iv0++) {
      const auto& v0 = v0List[sorted_v0[iv0]];

      if (!v0BuilderOpts.generatePhotonCandidates.value && v0.v0Type > 1) {
        // skip photons if not requested
        products.v0dataLink(-1, -1);
        continue;
      }

      if (!baseOpts.mEnabledTables[kV0CoresBase] && v0Map[iv0] == -2) {
        // this v0 hasn't been used by cascades and we're not generating V0s, so skip it
        products.v0dataLink(-1, -1);
        continue;
      }

      // Get tracks and generate candidate
      // if collisionId positive: get vertex, negative: origin
      // could be replaced by mean vertex (but without much benefit...)
      float pvX = 0.0f, pvY = 0.0f, pvZ = 0.0f;
      if (v0.collisionId >= 0) {
        auto const& collision = collisions.rawIteratorAt(v0.collisionId);
        pvX = collision.posX();
        pvY = collision.posY();
        pvZ = collision.posZ();
      }
      auto const& posTrack = tracks.rawIteratorAt(v0.posTrackId);
      auto const& negTrack = tracks.rawIteratorAt(v0.negTrackId);

      auto posTrackPar = getTrackParCov(posTrack);
      auto negTrackPar = getTrackParCov(negTrack);

      // handle TPC-only tracks properly (photon conversions)
      if (v0BuilderOpts.moveTPCOnlyTracks) {
        bool isPosTPCOnly = (posTrack.hasTPC() && !posTrack.hasITS() && !posTrack.hasTRD() && !posTrack.hasTOF());
        if (isPosTPCOnly) {
          // Nota bene: positive is TPC-only -> this entire V0 merits treatment as photon candidate
          posTrackPar.setPID(o2::track::PID::Electron);
          negTrackPar.setPID(o2::track::PID::Electron);

          auto const& collision = collisions.rawIteratorAt(v0.collisionId);
          if (!mVDriftMgr.moveTPCTrack<TBCs, TCollisions>(collision, posTrack, posTrackPar)) {
            products.v0dataLink(-1, -1);
            continue;
          }
        }

        bool isNegTPCOnly = (negTrack.hasTPC() && !negTrack.hasITS() && !negTrack.hasTRD() && !negTrack.hasTOF());
        if (isNegTPCOnly) {
          // Nota bene: negative is TPC-only -> this entire V0 merits treatment as photon candidate
          posTrackPar.setPID(o2::track::PID::Electron);
          negTrackPar.setPID(o2::track::PID::Electron);

          auto const& collision = collisions.rawIteratorAt(v0.collisionId);
          if (!mVDriftMgr.moveTPCTrack<TBCs, TCollisions>(collision, negTrack, negTrackPar)) {
            products.v0dataLink(-1, -1);
            continue;
          }
        }
      }

      if (!straHelper.buildV0Candidate(v0.collisionId, pvX, pvY, pvZ, posTrack, negTrack, posTrackPar, negTrackPar, v0.isCollinearV0, baseOpts.mEnabledTables[kV0Covs], v0BuilderOpts.generatePhotonCandidates)) {
        products.v0dataLink(-1, -1);
        continue;
      }
      if constexpr (requires { posTrack.tpcNSigmaEl(); }) {
        if (preSelectOpts.preselectOnlyDesiredV0s) {
          float lPt = RecoDecay::sqrtSumOfSquares(
            straHelper.v0.positiveMomentum[0] + straHelper.v0.negativeMomentum[0],
            straHelper.v0.positiveMomentum[1] + straHelper.v0.negativeMomentum[1]);

          float lPtot = RecoDecay::sqrtSumOfSquares(
            straHelper.v0.positiveMomentum[0] + straHelper.v0.negativeMomentum[0],
            straHelper.v0.positiveMomentum[1] + straHelper.v0.negativeMomentum[1],
            straHelper.v0.positiveMomentum[2] + straHelper.v0.negativeMomentum[2]);

          float lLengthTraveled = RecoDecay::sqrtSumOfSquares(
            straHelper.v0.position[0] - pvX,
            straHelper.v0.position[1] - pvY,
            straHelper.v0.position[2] - pvZ);

          uint8_t maskV0Preselection = 0;

          if ( // photon PID, mass, lifetime selection
            std::abs(posTrack.tpcNSigmaEl()) < preSelectOpts.maxTPCpidNsigma &&
            std::abs(negTrack.tpcNSigmaEl()) < preSelectOpts.maxTPCpidNsigma &&
            std::abs(straHelper.v0.massGamma) < preSelectOpts.massCutPhoton) {
            BITSET(maskV0Preselection, selGamma);
          }

          if ( // K0Short PID, mass, lifetime selection
            std::abs(posTrack.tpcNSigmaPi()) < preSelectOpts.maxTPCpidNsigma &&
            std::abs(negTrack.tpcNSigmaPi()) < preSelectOpts.maxTPCpidNsigma &&
            o2::constants::physics::MassKaonNeutral * lLengthTraveled / (lPtot + 1e-13) < preSelectOpts.lifetimeCut->get("lifetimeCutK0S") &&
            std::abs(straHelper.v0.massK0Short - o2::constants::physics::MassKaonNeutral) < preSelectOpts.massWindownumberOfSigmas * getMassSigmaK0Short(lPt) + preSelectOpts.massWindowSafetyMargin) {
            BITSET(maskV0Preselection, selK0Short);
          }

          if ( // Lambda PID, mass, lifetime selection
            std::abs(posTrack.tpcNSigmaPr()) < preSelectOpts.maxTPCpidNsigma &&
            std::abs(negTrack.tpcNSigmaPi()) < preSelectOpts.maxTPCpidNsigma &&
            o2::constants::physics::MassLambda * lLengthTraveled / (lPtot + 1e-13) < preSelectOpts.lifetimeCut->get("lifetimeCutLambda") &&
            std::abs(straHelper.v0.massLambda - o2::constants::physics::MassLambda) < preSelectOpts.massWindownumberOfSigmas * getMassSigmaLambda(lPt) + preSelectOpts.massWindowSafetyMargin) {
            BITSET(maskV0Preselection, selLambda);
          }

          if ( // antiLambda PID, mass, lifetime selection
            std::abs(posTrack.tpcNSigmaPi()) < preSelectOpts.maxTPCpidNsigma &&
            std::abs(negTrack.tpcNSigmaPr()) < preSelectOpts.maxTPCpidNsigma &&
            o2::constants::physics::MassLambda * lLengthTraveled / (lPtot + 1e-13) < preSelectOpts.lifetimeCut->get("lifetimeCutLambda") &&
            std::abs(straHelper.v0.massAntiLambda - o2::constants::physics::MassLambda) < preSelectOpts.massWindownumberOfSigmas * getMassSigmaLambda(lPt) + preSelectOpts.massWindowSafetyMargin) {
            BITSET(maskV0Preselection, selAntiLambda);
          }

          histos.fill(HIST("hPreselectionV0s"), maskV0Preselection);

          if (maskV0Preselection == 0) {
            products.v0dataLink(-1, -1);
            continue;
          }
        }
      }
      if (v0Map[iv0] == -1 && baseOpts.useV0BufferForCascades) {
        v0Map[iv0] = v0sFromCascades.size(); // provide actual valid index in buffer
        v0sFromCascades.push_back(straHelper.v0);
      }
      // fill requested cursors only if type is not 0
      if (v0.v0Type == 1 || (v0.v0Type > 1 && v0BuilderOpts.generatePhotonCandidates)) {
        nV0s++;
        if (baseOpts.mEnabledTables[kV0Indices]) {
          // for referencing (especially - but not only - when using derived data)
          products.v0indices(v0.posTrackId, v0.negTrackId,
                             v0.collisionId, iv0);
          histos.fill(HIST("hTableBuildingStatistics"), kV0Indices);
        }
        if (baseOpts.mEnabledTables[kV0TrackXs]) {
          // further decay chains may need this
          products.v0trackXs(straHelper.v0.positiveTrackX, straHelper.v0.negativeTrackX);
          histos.fill(HIST("hTableBuildingStatistics"), kV0TrackXs);
        }
        if (baseOpts.mEnabledTables[kV0CoresBase]) {
          // standard analysis
          products.v0cores(straHelper.v0.position[0], straHelper.v0.position[1], straHelper.v0.position[2],
                           straHelper.v0.positiveMomentum[0], straHelper.v0.positiveMomentum[1], straHelper.v0.positiveMomentum[2],
                           straHelper.v0.negativeMomentum[0], straHelper.v0.negativeMomentum[1], straHelper.v0.negativeMomentum[2],
                           straHelper.v0.daughterDCA,
                           straHelper.v0.positiveDCAxy,
                           straHelper.v0.negativeDCAxy,
                           TMath::Cos(straHelper.v0.pointingAngle),
                           straHelper.v0.dcaToPV,
                           v0.v0Type);
          products.v0dataLink(products.v0cores.lastIndex(), -1);
          histos.fill(HIST("hTableBuildingStatistics"), kV0CoresBase);
        }
        if (baseOpts.mEnabledTables[kV0TraPosAtDCAs]) {
          // for tracking studies
          products.v0dauPositions(straHelper.v0.positivePosition[0], straHelper.v0.positivePosition[1], straHelper.v0.positivePosition[2],
                                  straHelper.v0.negativePosition[0], straHelper.v0.negativePosition[1], straHelper.v0.negativePosition[2]);
          histos.fill(HIST("hTableBuildingStatistics"), kV0TraPosAtDCAs);
        }
        if (baseOpts.mEnabledTables[kV0TraPosAtIUs]) {
          // for tracking studies
          std::array<float, 3> positivePositionIU;
          std::array<float, 3> negativePositionIU;
          o2::track::TrackPar positiveTrackParam = getTrackPar(posTrack);
          o2::track::TrackPar negativeTrackParam = getTrackPar(negTrack);
          positiveTrackParam.getXYZGlo(positivePositionIU);
          negativeTrackParam.getXYZGlo(negativePositionIU);
          products.v0dauPositionsIU(positivePositionIU[0], positivePositionIU[1], positivePositionIU[2],
                                    negativePositionIU[0], negativePositionIU[1], negativePositionIU[2]);
          histos.fill(HIST("hTableBuildingStatistics"), kV0TraPosAtIUs);
        }
        if (baseOpts.mEnabledTables[kV0Covs]) {
          products.v0covs(straHelper.v0.positionCovariance, straHelper.v0.momentumCovariance);
          histos.fill(HIST("hTableBuildingStatistics"), kV0Covs);
        }

        //_________________________________________________________
        // MC handling part
        if constexpr (soa::is_table<TMCParticles>) {
          // only worry about this if someone else worried about this
          if ((baseOpts.mEnabledTables[kV0MCCores] || baseOpts.mEnabledTables[kMcV0Labels] || baseOpts.mEnabledTables[kV0MCCollRefs])) {
            thisInfo.label = -1;
            thisInfo.motherLabel = -1;
            thisInfo.pdgCode = 0;
            thisInfo.pdgCodeMother = 0;
            thisInfo.pdgCodePositive = 0;
            thisInfo.pdgCodeNegative = 0;
            thisInfo.mcCollision = -1;
            thisInfo.xyz[0] = thisInfo.xyz[1] = thisInfo.xyz[2] = 0.0f;
            thisInfo.posP[0] = thisInfo.posP[1] = thisInfo.posP[2] = 0.0f;
            thisInfo.negP[0] = thisInfo.negP[1] = thisInfo.negP[2] = 0.0f;
            thisInfo.momentum[0] = thisInfo.momentum[1] = thisInfo.momentum[2] = 0.0f;

            // Association check
            // There might be smarter ways of doing this in the future
            if (negTrack.has_mcParticle() && posTrack.has_mcParticle()) {
              auto lMCNegTrack = negTrack.template mcParticle_as<aod::McParticles>();
              auto lMCPosTrack = posTrack.template mcParticle_as<aod::McParticles>();

              thisInfo.pdgCodePositive = lMCPosTrack.pdgCode();
              thisInfo.pdgCodeNegative = lMCNegTrack.pdgCode();
              thisInfo.processPositive = lMCPosTrack.getProcess();
              thisInfo.processNegative = lMCNegTrack.getProcess();
              thisInfo.posP[0] = lMCPosTrack.px();
              thisInfo.posP[1] = lMCPosTrack.py();
              thisInfo.posP[2] = lMCPosTrack.pz();
              thisInfo.negP[0] = lMCNegTrack.px();
              thisInfo.negP[1] = lMCNegTrack.py();
              thisInfo.negP[2] = lMCNegTrack.pz();

              // check for pi -> mu + antineutrino decay
              // if present, de-reference original V0 correctly and provide label to original object
              // NOTA BENE: the prong info will still correspond to a muon, treat carefully!
              int negOriginating = -1, posOriginating = -1, particleForDecayPositionIdx = -1;
              negOriginating = getOriginatingParticle(lMCNegTrack, particleForDecayPositionIdx, v0BuilderOpts.mc_treatPiToMuDecays);
              posOriginating = getOriginatingParticle(lMCPosTrack, particleForDecayPositionIdx, v0BuilderOpts.mc_treatPiToMuDecays);

              if (negOriginating > -1 && negOriginating == posOriginating) {
                auto originatingV0 = mcParticles.rawIteratorAt(negOriginating);
                auto particleForDecayPosition = mcParticles.rawIteratorAt(particleForDecayPositionIdx);

                thisInfo.label = originatingV0.globalIndex();
                thisInfo.xyz[0] = particleForDecayPosition.vx();
                thisInfo.xyz[1] = particleForDecayPosition.vy();
                thisInfo.xyz[2] = particleForDecayPosition.vz();

                if (originatingV0.has_mcCollision()) {
                  thisInfo.mcCollision = originatingV0.mcCollisionId(); // save this reference, please
                }

                // acquire information
                thisInfo.pdgCode = originatingV0.pdgCode();
                thisInfo.isPhysicalPrimary = originatingV0.isPhysicalPrimary();
                thisInfo.momentum[0] = originatingV0.px();
                thisInfo.momentum[1] = originatingV0.py();
                thisInfo.momentum[2] = originatingV0.pz();

                if (originatingV0.has_mothers()) {
                  for (const auto& lV0Mother : originatingV0.template mothers_as<aod::McParticles>()) {
                    thisInfo.pdgCodeMother = lV0Mother.pdgCode();
                    thisInfo.motherLabel = lV0Mother.globalIndex();
                  }
                }
              }

            } // end association check
            // Construct label table (note: this will be joinable with V0Datas!)
            if (baseOpts.mEnabledTables[kMcV0Labels]) {
              products.v0labels(thisInfo.label, thisInfo.motherLabel);
              histos.fill(HIST("hTableBuildingStatistics"), kMcV0Labels);
            }

            // Construct found tag
            if (baseOpts.mEnabledTables[kV0FoundTags]) {
              products.v0FoundTag(v0.found);
              histos.fill(HIST("hTableBuildingStatistics"), kV0FoundTags);
            }

            // Mark mcParticle as recoed (no searching necessary afterwards)
            if (thisInfo.label > -1) {
              mcParticleIsReco[thisInfo.label] = true;
            }

            // ---] Symmetric populate [---
            // in this approach, V0Cores will be joinable with V0MCCores.
            // this is the most pedagogical approach, but it is also more limited
            // and it might use more disk space unnecessarily.
            if (v0BuilderOpts.mc_populateV0MCCoresSymmetric) {
              if (baseOpts.mEnabledTables[kV0MCCores]) {
                products.v0mccores(
                  thisInfo.label, thisInfo.pdgCode,
                  thisInfo.pdgCodeMother, thisInfo.pdgCodePositive, thisInfo.pdgCodeNegative,
                  thisInfo.isPhysicalPrimary, thisInfo.xyz[0], thisInfo.xyz[1], thisInfo.xyz[2],
                  thisInfo.posP[0], thisInfo.posP[1], thisInfo.posP[2],
                  thisInfo.negP[0], thisInfo.negP[1], thisInfo.negP[2],
                  thisInfo.momentum[0], thisInfo.momentum[1], thisInfo.momentum[2]);
                histos.fill(HIST("hTableBuildingStatistics"), kV0MCCores);
                histos.fill(HIST("hPrimaryV0s"), 0);
                if (thisInfo.isPhysicalPrimary)
                  histos.fill(HIST("hPrimaryV0s"), 1);
              }
              if (baseOpts.mEnabledTables[kV0MCCollRefs]) {
                products.v0mccollref(thisInfo.mcCollision);
                histos.fill(HIST("hTableBuildingStatistics"), kV0MCCollRefs);
              }

              // n.b. placing the interlink index here allows for the writing of
              //      code that is agnostic with respect to the joinability of
              //      V0Cores and V0MCCores (always dereference -> safe)
              if (baseOpts.mEnabledTables[kV0CoreMCLabels]) {
                products.v0CoreMCLabels(iv0); // interlink index
                histos.fill(HIST("hTableBuildingStatistics"), kV0CoreMCLabels);
              }
            }
            // ---] Asymmetric populate [---
            // in this approach, V0Cores will NOT be joinable with V0MCCores.
            // an additional reference to V0MCCore that IS joinable with V0Cores
            // will be provided to the user.
            if (v0BuilderOpts.mc_populateV0MCCoresAsymmetric) {
              int thisV0MCCoreIndex = -1;
              // step 1: check if this element is already provided in the table
              //         using the packedIndices variable calculated above
              for (uint32_t ii = 0; ii < mcV0infos.size(); ii++) {
                if (thisInfo.label == mcV0infos[ii].label && mcV0infos[ii].label > -1) {
                  thisV0MCCoreIndex = ii;
                  break; // this exists already in list
                }
              }
              if (thisV0MCCoreIndex < 0 && thisInfo.label > -1) {
                // this V0MCCore does not exist yet. Create it and reference it
                thisV0MCCoreIndex = mcV0infos.size();
                mcV0infos.push_back(thisInfo);
              }
              if (baseOpts.mEnabledTables[kV0CoreMCLabels]) {
                products.v0CoreMCLabels(thisV0MCCoreIndex); // interlink index
                histos.fill(HIST("hTableBuildingStatistics"), kV0CoreMCLabels);
              }
            }
          } // enabled tables check
        } // constexpr requires check
      } else {
        products.v0dataLink(-1, -1);
      }
    }

    // finish populating V0MCCores if in asymmetric mode
    if constexpr (soa::is_table<TMCParticles>) {
      if (v0BuilderOpts.mc_populateV0MCCoresAsymmetric && (baseOpts.mEnabledTables[kV0MCCores] || baseOpts.mEnabledTables[kV0MCCollRefs])) {
        // first step: add any un-recoed v0mmcores that were requested
        for (const auto& mcParticle : mcParticles) {
          thisInfo.label = -1;
          thisInfo.motherLabel = -1;
          thisInfo.pdgCode = 0;
          thisInfo.pdgCodeMother = -1;
          thisInfo.pdgCodePositive = -1;
          thisInfo.pdgCodeNegative = -1;
          thisInfo.mcCollision = -1;
          thisInfo.xyz[0] = thisInfo.xyz[1] = thisInfo.xyz[2] = 0.0f;
          thisInfo.posP[0] = thisInfo.posP[1] = thisInfo.posP[2] = 0.0f;
          thisInfo.negP[0] = thisInfo.negP[1] = thisInfo.negP[2] = 0.0f;
          thisInfo.momentum[0] = thisInfo.momentum[1] = thisInfo.momentum[2] = 0.0f;

          if (mcParticleIsReco[mcParticle.globalIndex()] == true)
            continue; // skip if already created in list

          if (std::fabs(mcParticle.y()) > v0BuilderOpts.mc_rapidityWindow)
            continue; // skip outside midrapidity

          if (v0BuilderOpts.mc_keepOnlyPhysicalPrimary && !mcParticle.isPhysicalPrimary())
            continue; // skip secondary MC V0s

          if (
            (v0BuilderOpts.mc_addGeneratedK0Short && mcParticle.pdgCode() == PDG_t::kK0Short) ||
            (v0BuilderOpts.mc_addGeneratedLambda && mcParticle.pdgCode() == PDG_t::kLambda0) ||
            (v0BuilderOpts.mc_addGeneratedAntiLambda && mcParticle.pdgCode() == PDG_t::kLambda0Bar) ||
            (v0BuilderOpts.mc_addGeneratedGamma && mcParticle.pdgCode() == PDG_t::kGamma)) {
            thisInfo.pdgCode = mcParticle.pdgCode();
            thisInfo.isPhysicalPrimary = mcParticle.isPhysicalPrimary();
            thisInfo.label = mcParticle.globalIndex();

            if (mcParticle.has_mcCollision()) {
              thisInfo.mcCollision = mcParticle.mcCollisionId(); // save this reference, please
            }

            //
            thisInfo.momentum[0] = mcParticle.px();
            thisInfo.momentum[1] = mcParticle.py();
            thisInfo.momentum[2] = mcParticle.pz();

            if (mcParticle.has_mothers()) {
              auto const& mother = mcParticle.template mothers_first_as<aod::McParticles>();
              thisInfo.pdgCodeMother = mother.pdgCode();
              thisInfo.motherLabel = mother.globalIndex();
            }
            if (mcParticle.has_daughters()) {
              auto const& daughters = mcParticle.template daughters_as<aod::McParticles>();

              for (const auto& dau : daughters) {
                if (dau.getProcess() != TMCProcess::kPDecay)
                  continue;

                if (dau.pdgCode() > 0) {
                  thisInfo.pdgCodePositive = dau.pdgCode();
                  thisInfo.processPositive = dau.getProcess();
                  thisInfo.posP[0] = dau.px();
                  thisInfo.posP[1] = dau.py();
                  thisInfo.posP[2] = dau.pz();
                  thisInfo.xyz[0] = dau.vx();
                  thisInfo.xyz[1] = dau.vy();
                  thisInfo.xyz[2] = dau.vz();
                }
                if (dau.pdgCode() < 0) {
                  thisInfo.pdgCodeNegative = dau.pdgCode();
                  thisInfo.processNegative = dau.getProcess();
                  thisInfo.negP[0] = dau.px();
                  thisInfo.negP[1] = dau.py();
                  thisInfo.negP[2] = dau.pz();
                }
              }
            }

            // if I got here, it means this MC particle was not recoed and is of interest. Add it please
            mcV0infos.push_back(thisInfo);
          }
        }

        for (const auto& info : mcV0infos) {
          if (baseOpts.mEnabledTables[kV0MCCores]) {
            products.v0mccores(
              info.label, info.pdgCode,
              info.pdgCodeMother, info.pdgCodePositive, info.pdgCodeNegative,
              info.isPhysicalPrimary, info.xyz[0], info.xyz[1], info.xyz[2],
              info.posP[0], info.posP[1], info.posP[2],
              info.negP[0], info.negP[1], info.negP[2],
              info.momentum[0], info.momentum[1], info.momentum[2]);
            histos.fill(HIST("hTableBuildingStatistics"), kV0MCCores);
            histos.fill(HIST("hPrimaryV0s"), 0);
            if (info.isPhysicalPrimary)
              histos.fill(HIST("hPrimaryV0s"), 1);
          }
          if (baseOpts.mEnabledTables[kV0MCCollRefs]) {
            products.v0mccollref(info.mcCollision);
            histos.fill(HIST("hTableBuildingStatistics"), kV0MCCollRefs);
          }
        }
      } // end V0MCCores filling in case of MC
    } // end constexpr requires mcParticles

    LOGF(debug, "V0s in DF: %i, V0s built: %i, V0s built and buffered for cascades: %i.", v0s.size(), nV0s, v0sFromCascades.size());
  }

  //__________________________________________________
  template <typename TTrack, typename TMCParticles>
  void extractMonteCarloProperties(TTrack const& posTrack, TTrack const& negTrack, TTrack const& bachTrack, TMCParticles const& mcParticles)
  {
    // encapsulates acquisition of MC properties from MC
    thisCascInfo.pdgCode = -1, thisCascInfo.pdgCodeMother = -1;
    thisCascInfo.pdgCodePositive = -1, thisCascInfo.pdgCodeNegative = -1;
    thisCascInfo.pdgCodeBachelor = -1, thisCascInfo.pdgCodeV0 = -1;
    thisCascInfo.isPhysicalPrimary = false;
    thisCascInfo.xyz[0] = -999.0f, thisCascInfo.xyz[1] = -999.0f, thisCascInfo.xyz[2] = -999.0f;
    thisCascInfo.lxyz[0] = -999.0f, thisCascInfo.lxyz[1] = -999.0f, thisCascInfo.lxyz[2] = -999.0f;
    thisCascInfo.posP[0] = -999.0f, thisCascInfo.posP[1] = -999.0f, thisCascInfo.posP[2] = -999.0f;
    thisCascInfo.negP[0] = -999.0f, thisCascInfo.negP[1] = -999.0f, thisCascInfo.negP[2] = -999.0f;
    thisCascInfo.bachP[0] = -999.0f, thisCascInfo.bachP[1] = -999.0f, thisCascInfo.bachP[2] = -999.0f;
    thisCascInfo.momentum[0] = -999.0f, thisCascInfo.momentum[1] = -999.0f, thisCascInfo.momentum[2] = -999.0f;
    thisCascInfo.label = -1, thisCascInfo.motherLabel = -1;
    thisCascInfo.mcParticlePositive = -1;
    thisCascInfo.mcParticleNegative = -1;
    thisCascInfo.mcParticleBachelor = -1;

    // Association check
    // There might be smarter ways of doing this in the future
    if (negTrack.has_mcParticle() && posTrack.has_mcParticle() && bachTrack.has_mcParticle()) {
      auto lMCBachTrack = bachTrack.template mcParticle_as<aod::McParticles>();
      auto lMCNegTrack = negTrack.template mcParticle_as<aod::McParticles>();
      auto lMCPosTrack = posTrack.template mcParticle_as<aod::McParticles>();

      thisCascInfo.mcParticlePositive = lMCPosTrack.globalIndex();
      thisCascInfo.mcParticleNegative = lMCNegTrack.globalIndex();
      thisCascInfo.mcParticleBachelor = lMCBachTrack.globalIndex();
      thisCascInfo.pdgCodePositive = lMCPosTrack.pdgCode();
      thisCascInfo.pdgCodeNegative = lMCNegTrack.pdgCode();
      thisCascInfo.pdgCodeBachelor = lMCBachTrack.pdgCode();
      thisCascInfo.posP[0] = lMCPosTrack.px();
      thisCascInfo.posP[1] = lMCPosTrack.py();
      thisCascInfo.posP[2] = lMCPosTrack.pz();
      thisCascInfo.negP[0] = lMCNegTrack.px();
      thisCascInfo.negP[1] = lMCNegTrack.py();
      thisCascInfo.negP[2] = lMCNegTrack.pz();
      thisCascInfo.bachP[0] = lMCBachTrack.px();
      thisCascInfo.bachP[1] = lMCBachTrack.py();
      thisCascInfo.bachP[2] = lMCBachTrack.pz();
      thisCascInfo.processPositive = lMCPosTrack.getProcess();
      thisCascInfo.processNegative = lMCNegTrack.getProcess();
      thisCascInfo.processBachelor = lMCBachTrack.getProcess();

      // Step 0: treat pi -> mu + antineutrino
      // if present, de-reference original V0 correctly and provide label to original object
      // NOTA BENE: the prong info will still correspond to a muon, treat carefully!
      int negOriginating = -1, posOriginating = -1, bachOriginating = -1;
      int particleForLambdaDecayPositionIdx = -1, particleForCascadeDecayPositionIdx = -1;
      negOriginating = getOriginatingParticle(lMCNegTrack, particleForLambdaDecayPositionIdx, cascadeBuilderOpts.mc_treatPiToMuDecays);
      posOriginating = getOriginatingParticle(lMCPosTrack, particleForLambdaDecayPositionIdx, cascadeBuilderOpts.mc_treatPiToMuDecays);
      bachOriginating = getOriginatingParticle(lMCBachTrack, particleForCascadeDecayPositionIdx, cascadeBuilderOpts.mc_treatPiToMuDecays);

      if (negOriginating > -1 && negOriginating == posOriginating) {
        auto originatingV0 = mcParticles.rawIteratorAt(negOriginating);
        auto particleForLambdaDecayPosition = mcParticles.rawIteratorAt(particleForLambdaDecayPositionIdx);

        thisCascInfo.label = originatingV0.globalIndex();
        thisCascInfo.lxyz[0] = particleForLambdaDecayPosition.vx();
        thisCascInfo.lxyz[1] = particleForLambdaDecayPosition.vy();
        thisCascInfo.lxyz[2] = particleForLambdaDecayPosition.vz();
        thisCascInfo.pdgCodeV0 = originatingV0.pdgCode();

        if (originatingV0.has_mothers()) {
          for (const auto& lV0Mother : originatingV0.template mothers_as<aod::McParticles>()) {
            if (lV0Mother.globalIndex() == bachOriginating) { // found mother particle
              thisCascInfo.label = lV0Mother.globalIndex();

              if (lV0Mother.has_mcCollision()) {
                thisCascInfo.mcCollision = lV0Mother.mcCollisionId(); // save this reference, please
              }

              thisCascInfo.pdgCode = lV0Mother.pdgCode();
              thisCascInfo.isPhysicalPrimary = lV0Mother.isPhysicalPrimary();
              thisCascInfo.xyz[0] = originatingV0.vx();
              thisCascInfo.xyz[1] = originatingV0.vy();
              thisCascInfo.xyz[2] = originatingV0.vz();
              thisCascInfo.momentum[0] = lV0Mother.px();
              thisCascInfo.momentum[1] = lV0Mother.py();
              thisCascInfo.momentum[2] = lV0Mother.pz();
              if (lV0Mother.has_mothers()) {
                for (const auto& lV0GrandMother : lV0Mother.template mothers_as<aod::McParticles>()) {
                  thisCascInfo.pdgCodeMother = lV0GrandMother.pdgCode();
                  thisCascInfo.motherLabel = lV0GrandMother.globalIndex();
                }
              }
            }
          } // end v0 mother loop
        } // end has_mothers check for V0
      } // end conditional of pos/neg originating being the same
    } // end association check
  }

  //__________________________________________________
  template <typename THistoRegistry, typename TCollisions, typename TTracks, typename TCascades, typename TMCParticles, typename TProducts>
  void buildCascades(THistoRegistry& histos, TCollisions const& collisions, TCascades const& cascades, TTracks const& tracks, TMCParticles const& mcParticles, TProducts& products)
  {
    // prepare MC containers (not necessarily used)
    std::vector<mcCascinfo> mcCascinfos; // V0MCCore information
    std::vector<bool> mcParticleIsReco;

    if constexpr (soa::is_table<TMCParticles>) {
      // do this if provided with a mcParticle table as well
      mcParticleIsReco.resize(mcParticles.size(), false);
    }

    if (!baseOpts.mEnabledTables[kStoredCascCores]) {
      return; // don't do if no request for cascades in place
    }
    int nCascades = 0;
    // Loops over all cascades in the time frame
    histos.fill(HIST("hInputStatistics"), kStoredCascCores, cascades.size());
    for (size_t icascade = 0; icascade < cascades.size(); icascade++) {
      // Get tracks and generate candidate
      auto const& cascade = cascades[sorted_cascade[icascade]];
      // if collisionId positive: get vertex, negative: origin
      // could be replaced by mean vertex (but without much benefit...)
      float pvX = 0.0f, pvY = 0.0f, pvZ = 0.0f;
      if (cascade.collisionId >= 0) {
        auto const& collision = collisions.rawIteratorAt(cascade.collisionId);
        pvX = collision.posX();
        pvY = collision.posY();
        pvZ = collision.posZ();
      }
      auto const& posTrack = tracks.rawIteratorAt(cascade.posTrackId);
      auto const& negTrack = tracks.rawIteratorAt(cascade.negTrackId);
      auto const& bachTrack = tracks.rawIteratorAt(cascade.bachTrackId);
      if (baseOpts.useV0BufferForCascades) {
        // this processing path uses a buffer of V0s so that no
        // additional minimization step is redone. It consumes less
        // CPU at the cost of more memory. Since memory is a more
        // limited commodity, this isn't the default option.

        // check if cached - if not, skip
        if (cascade.v0Id < 0 || v0Map[cascade.v0Id] < 0) {
          // this V0 hasn't been stored / cached
          products.cascdataLink(-1);
          interlinks.cascadeToCascCores.push_back(-1);
          continue; // didn't work out, skip
        }

        if (!straHelper.buildCascadeCandidate(cascade.collisionId, pvX, pvY, pvZ,
                                              v0sFromCascades[v0Map[cascade.v0Id]],
                                              posTrack,
                                              negTrack,
                                              bachTrack,
                                              baseOpts.mEnabledTables[kCascBBs],
                                              cascadeBuilderOpts.useCascadeMomentumAtPrimVtx,
                                              baseOpts.mEnabledTables[kCascCovs])) {
          products.cascdataLink(-1);
          interlinks.cascadeToCascCores.push_back(-1);
          continue; // didn't work out, skip
        }
      } else {
        // this processing path generates the entire cascade
        // from tracks, without any need to have V0s generated.
        if (!straHelper.buildCascadeCandidate(cascade.collisionId, pvX, pvY, pvZ,
                                              posTrack,
                                              negTrack,
                                              bachTrack,
                                              baseOpts.mEnabledTables[kCascBBs],
                                              cascadeBuilderOpts.useCascadeMomentumAtPrimVtx,
                                              baseOpts.mEnabledTables[kCascCovs])) {
          products.cascdataLink(-1);
          interlinks.cascadeToCascCores.push_back(-1);
          continue; // didn't work out, skip
        }
      }
      nCascades++;

      if constexpr (requires { posTrack.tpcNSigmaEl(); }) {
        if (preSelectOpts.preselectOnlyDesiredCascades) {
          float lPt = RecoDecay::sqrtSumOfSquares(
            straHelper.cascade.bachelorMomentum[0] + straHelper.cascade.positiveMomentum[0] + straHelper.cascade.negativeMomentum[0],
            straHelper.cascade.bachelorMomentum[1] + straHelper.cascade.positiveMomentum[1] + straHelper.cascade.negativeMomentum[1]);

          float lPtot = RecoDecay::sqrtSumOfSquares(
            straHelper.cascade.bachelorMomentum[0] + straHelper.cascade.positiveMomentum[0] + straHelper.cascade.negativeMomentum[0],
            straHelper.cascade.bachelorMomentum[1] + straHelper.cascade.positiveMomentum[1] + straHelper.cascade.negativeMomentum[1],
            straHelper.cascade.bachelorMomentum[2] + straHelper.cascade.positiveMomentum[2] + straHelper.cascade.negativeMomentum[2]);

          float lV0Ptot = RecoDecay::sqrtSumOfSquares(
            straHelper.cascade.positiveMomentum[0] + straHelper.cascade.negativeMomentum[0],
            straHelper.cascade.positiveMomentum[1] + straHelper.cascade.negativeMomentum[1],
            straHelper.cascade.positiveMomentum[2] + straHelper.cascade.negativeMomentum[2]);

          float lLengthTraveled = RecoDecay::sqrtSumOfSquares(
            straHelper.cascade.cascadePosition[0] - pvX,
            straHelper.cascade.cascadePosition[1] - pvY,
            straHelper.cascade.cascadePosition[2] - pvZ);

          float lV0LengthTraveled = RecoDecay::sqrtSumOfSquares(
            straHelper.cascade.v0Position[0] - straHelper.cascade.cascadePosition[0],
            straHelper.cascade.v0Position[1] - straHelper.cascade.cascadePosition[1],
            straHelper.cascade.v0Position[2] - straHelper.cascade.cascadePosition[2]);

          uint8_t maskCascadePreselection = 0;

          if ( // XiMinus PID and mass selection
            straHelper.cascade.charge < 0 &&
            std::abs(posTrack.tpcNSigmaPr()) < preSelectOpts.maxTPCpidNsigma &&
            std::abs(negTrack.tpcNSigmaPi()) < preSelectOpts.maxTPCpidNsigma &&
            std::abs(bachTrack.tpcNSigmaPi()) < preSelectOpts.maxTPCpidNsigma &&
            o2::constants::physics::MassLambda * lV0LengthTraveled / (lV0Ptot + 1e-13) < preSelectOpts.lifetimeCut->get("lifetimeCutLambda") &&
            o2::constants::physics::MassXiMinus * lLengthTraveled / (lPtot + 1e-13) < preSelectOpts.lifetimeCut->get("lifetimeCutXi") &&
            std::abs(straHelper.cascade.massXi - o2::constants::physics::MassXiMinus) < preSelectOpts.massWindownumberOfSigmas * getMassSigmaXi(lPt) + preSelectOpts.massWindowSafetyMargin) {
            BITSET(maskCascadePreselection, selXiMinus);
          }

          if ( // XiPlus PID and mass selection
            straHelper.cascade.charge > 0 &&
            std::abs(posTrack.tpcNSigmaPi()) < preSelectOpts.maxTPCpidNsigma &&
            std::abs(negTrack.tpcNSigmaPr()) < preSelectOpts.maxTPCpidNsigma &&
            std::abs(bachTrack.tpcNSigmaPi()) < preSelectOpts.maxTPCpidNsigma &&
            o2::constants::physics::MassLambda * lV0LengthTraveled / (lV0Ptot + 1e-13) < preSelectOpts.lifetimeCut->get("lifetimeCutLambda") &&
            o2::constants::physics::MassXiMinus * lLengthTraveled / (lPtot + 1e-13) < preSelectOpts.lifetimeCut->get("lifetimeCutXi") &&
            std::abs(straHelper.cascade.massXi - o2::constants::physics::MassXiMinus) < preSelectOpts.massWindownumberOfSigmas * getMassSigmaXi(lPt) + preSelectOpts.massWindowSafetyMargin) {
            BITSET(maskCascadePreselection, selXiPlus);
          }

          if ( // OmegaMinus PID and mass selection
            straHelper.cascade.charge < 0 &&
            std::abs(posTrack.tpcNSigmaPr()) < preSelectOpts.maxTPCpidNsigma &&
            std::abs(negTrack.tpcNSigmaPi()) < preSelectOpts.maxTPCpidNsigma &&
            std::abs(bachTrack.tpcNSigmaKa()) < preSelectOpts.maxTPCpidNsigma &&
            o2::constants::physics::MassLambda * lV0LengthTraveled / (lV0Ptot + 1e-13) < preSelectOpts.lifetimeCut->get("lifetimeCutLambda") &&
            o2::constants::physics::MassOmegaMinus * lLengthTraveled / (lPtot + 1e-13) < preSelectOpts.lifetimeCut->get("lifetimeCutOmega") &&
            std::abs(straHelper.cascade.massOmega - o2::constants::physics::MassOmegaMinus) < preSelectOpts.massWindownumberOfSigmas * getMassSigmaOmega(lPt) + preSelectOpts.massWindowSafetyMargin) {
            BITSET(maskCascadePreselection, selOmegaMinus);
          }

          if ( // OmegaPlus PID and mass selection
            straHelper.cascade.charge > 0 &&
            std::abs(posTrack.tpcNSigmaPi()) < preSelectOpts.maxTPCpidNsigma &&
            std::abs(negTrack.tpcNSigmaPr()) < preSelectOpts.maxTPCpidNsigma &&
            std::abs(bachTrack.tpcNSigmaKa()) < preSelectOpts.maxTPCpidNsigma &&
            o2::constants::physics::MassLambda * lV0LengthTraveled / (lV0Ptot + 1e-13) < preSelectOpts.lifetimeCut->get("lifetimeCutLambda") &&
            o2::constants::physics::MassOmegaMinus * lLengthTraveled / (lPtot + 1e-13) < preSelectOpts.lifetimeCut->get("lifetimeCutOmega") &&
            std::abs(straHelper.cascade.massOmega - o2::constants::physics::MassOmegaMinus) < preSelectOpts.massWindownumberOfSigmas * getMassSigmaOmega(lPt) + preSelectOpts.massWindowSafetyMargin) {
            BITSET(maskCascadePreselection, selOmegaPlus);
          }

          histos.fill(HIST("hPreselectionCascades"), maskCascadePreselection);

          if (maskCascadePreselection == 0) {
            products.cascdataLink(-1);
            interlinks.cascadeToCascCores.push_back(-1);
            continue;
          }
        }
      }

      // generate analysis tables as required
      if (baseOpts.mEnabledTables[kCascIndices]) {
        products.cascidx(cascade.globalId,
                         straHelper.cascade.positiveTrack, straHelper.cascade.negativeTrack,
                         straHelper.cascade.bachelorTrack, straHelper.cascade.collisionId);
        histos.fill(HIST("hTableBuildingStatistics"), kCascIndices);
      }
      if (baseOpts.mEnabledTables[kStoredCascCores]) {
        products.cascdata(straHelper.cascade.charge, straHelper.cascade.massXi, straHelper.cascade.massOmega,
                          straHelper.cascade.cascadePosition[0], straHelper.cascade.cascadePosition[1], straHelper.cascade.cascadePosition[2],
                          straHelper.cascade.v0Position[0], straHelper.cascade.v0Position[1], straHelper.cascade.v0Position[2],
                          straHelper.cascade.positiveMomentum[0], straHelper.cascade.positiveMomentum[1], straHelper.cascade.positiveMomentum[2],
                          straHelper.cascade.negativeMomentum[0], straHelper.cascade.negativeMomentum[1], straHelper.cascade.negativeMomentum[2],
                          straHelper.cascade.bachelorMomentum[0], straHelper.cascade.bachelorMomentum[1], straHelper.cascade.bachelorMomentum[2],
                          straHelper.cascade.cascadeMomentum[0], straHelper.cascade.cascadeMomentum[1], straHelper.cascade.cascadeMomentum[2],
                          straHelper.cascade.v0DaughterDCA, straHelper.cascade.cascadeDaughterDCA,
                          straHelper.cascade.positiveDCAxy, straHelper.cascade.negativeDCAxy,
                          straHelper.cascade.bachelorDCAxy, straHelper.cascade.cascadeDCAxy, straHelper.cascade.cascadeDCAz);
        histos.fill(HIST("hTableBuildingStatistics"), kStoredCascCores);

        // interlink always produced if cascades generated
        products.cascdataLink(products.cascdata.lastIndex());
        interlinks.cascCoreToCascades.push_back(cascade.globalId);
        interlinks.cascadeToCascCores.push_back(products.cascdata.lastIndex());
      }

      if (baseOpts.mEnabledTables[kCascTrackXs]) {
        products.cascTrackXs(straHelper.cascade.positiveTrackX, straHelper.cascade.negativeTrackX, straHelper.cascade.bachelorTrackX);
        histos.fill(HIST("hTableBuildingStatistics"), kCascTrackXs);
      }
      if (baseOpts.mEnabledTables[kCascBBs]) {
        products.cascbb(straHelper.cascade.bachBaryonCosPA, straHelper.cascade.bachBaryonDCAxyToPV);
        histos.fill(HIST("hTableBuildingStatistics"), kCascBBs);
      }
      if (baseOpts.mEnabledTables[kCascCovs]) {
        products.casccovs(straHelper.cascade.covariance);
        histos.fill(HIST("hTableBuildingStatistics"), kCascCovs);
      }

      //_________________________________________________________
      // MC handling part
      if constexpr (soa::is_table<TMCParticles>) {
        // only worry about this if someone else worried about this
        if ((baseOpts.mEnabledTables[kCascMCCores] || baseOpts.mEnabledTables[kMcCascLabels] || baseOpts.mEnabledTables[kCascMCCollRefs])) {
          extractMonteCarloProperties(posTrack, negTrack, bachTrack, mcParticles);

          // Construct label table (note: this will be joinable with CascDatas)
          if (baseOpts.mEnabledTables[kMcCascLabels]) {
            products.casclabels(
              thisCascInfo.label, thisCascInfo.motherLabel);
            histos.fill(HIST("hTableBuildingStatistics"), kMcCascLabels);
          }

          // Construct found tag
          if (baseOpts.mEnabledTables[kCascFoundTags]) {
            products.cascFoundTag(cascade.found);
            histos.fill(HIST("hTableBuildingStatistics"), kCascFoundTags);
          }

          // Mark mcParticle as recoed (no searching necessary afterwards)
          if (thisCascInfo.label > -1) {
            mcParticleIsReco[thisCascInfo.label] = true;
          }

          if (cascadeBuilderOpts.mc_populateCascMCCoresSymmetric) {
            if (baseOpts.mEnabledTables[kCascMCCores]) {
              products.cascmccores(
                thisCascInfo.pdgCode, thisCascInfo.pdgCodeMother, thisCascInfo.pdgCodeV0, thisCascInfo.isPhysicalPrimary,
                thisCascInfo.pdgCodePositive, thisCascInfo.pdgCodeNegative, thisCascInfo.pdgCodeBachelor,
                thisCascInfo.xyz[0], thisCascInfo.xyz[1], thisCascInfo.xyz[2],
                thisCascInfo.lxyz[0], thisCascInfo.lxyz[1], thisCascInfo.lxyz[2],
                thisCascInfo.posP[0], thisCascInfo.posP[1], thisCascInfo.posP[2],
                thisCascInfo.negP[0], thisCascInfo.negP[1], thisCascInfo.negP[2],
                thisCascInfo.bachP[0], thisCascInfo.bachP[1], thisCascInfo.bachP[2],
                thisCascInfo.momentum[0], thisCascInfo.momentum[1], thisCascInfo.momentum[2]);
              histos.fill(HIST("hTableBuildingStatistics"), kCascMCCores);
            }
            if (baseOpts.mEnabledTables[kCascMCCollRefs]) {
              products.cascmccollrefs(thisCascInfo.mcCollision);
              histos.fill(HIST("hTableBuildingStatistics"), kCascMCCollRefs);
            }
          }

          if (cascadeBuilderOpts.mc_populateCascMCCoresAsymmetric) {
            int thisCascMCCoreIndex = -1;
            // step 1: check if this element is already provided in the table
            //         using the packedIndices variable calculated above
            for (uint32_t ii = 0; ii < mcCascinfos.size(); ii++) {
              if (thisCascInfo.label == mcCascinfos[ii].label && mcCascinfos[ii].label > -1) {
                thisCascMCCoreIndex = ii;
                break; // this exists already in list
              }
            }
            if (thisCascMCCoreIndex < 0) {
              // this CascMCCore does not exist yet. Create it and reference it
              thisCascMCCoreIndex = mcCascinfos.size();
              mcCascinfos.push_back(thisCascInfo);
            }
            if (baseOpts.mEnabledTables[kCascCoreMCLabels]) {
              products.cascCoreMClabels(thisCascMCCoreIndex); // interlink: reconstructed -> MC index
              histos.fill(HIST("hTableBuildingStatistics"), kCascCoreMCLabels);
            }
          }

        } // enabled tables check

        // if BB tags requested, generate them now
        if (baseOpts.mEnabledTables[kMcCascBBTags]) {
          bool bbTag = false;
          if (bachTrack.has_mcParticle()) {
            auto bachelorParticle = bachTrack.template mcParticle_as<aod::McParticles>();
            if (bachelorParticle.pdgCode() == PDG_t::kPiPlus) { // pi+, look for antiproton in negative prong
              if (negTrack.has_mcParticle()) {
                auto baryonParticle = negTrack.template mcParticle_as<aod::McParticles>();
                if (baryonParticle.has_mothers() && bachelorParticle.has_mothers() && baryonParticle.pdgCode() == PDG_t::kProtonBar) {
                  for (const auto& baryonMother : baryonParticle.template mothers_as<aod::McParticles>()) {
                    for (const auto& pionMother : bachelorParticle.template mothers_as<aod::McParticles>()) {
                      if (baryonMother.globalIndex() == pionMother.globalIndex() && baryonMother.pdgCode() == PDG_t::kLambda0Bar) {
                        bbTag = true;
                      }
                    }
                  }
                }
              }
            } // end if-pion
            if (bachelorParticle.pdgCode() == PDG_t::kPiMinus) { // pi-, look for proton in positive prong
              if (posTrack.has_mcParticle()) {
                auto baryonParticle = posTrack.template mcParticle_as<aod::McParticles>();
                if (baryonParticle.has_mothers() && bachelorParticle.has_mothers() && baryonParticle.pdgCode() == PDG_t::kProton) {
                  for (const auto& baryonMother : baryonParticle.template mothers_as<aod::McParticles>()) {
                    for (const auto& pionMother : bachelorParticle.template mothers_as<aod::McParticles>()) {
                      if (baryonMother.globalIndex() == pionMother.globalIndex() && baryonMother.pdgCode() == PDG_t::kLambda0) {
                        bbTag = true;
                      }
                    }
                  }
                }
              }
            } // end if-pion
          } // end bachelor has mcparticle
          // Construct label table (note: this will be joinable with CascDatas)
          products.bbtags(bbTag);
          histos.fill(HIST("hTableBuildingStatistics"), kMcCascBBTags);
        } // end BB tag table enabled check

      } // constexpr requires mcParticles check
    } // cascades loop

    //_________________________________________________________
    // MC handling part
    if constexpr (soa::is_table<TMCParticles>) {
      if ((baseOpts.mEnabledTables[kCascMCCores] || baseOpts.mEnabledTables[kMcCascLabels] || baseOpts.mEnabledTables[kCascMCCollRefs])) {
        // now populate V0MCCores if in asymmetric mode
        if (cascadeBuilderOpts.mc_populateCascMCCoresAsymmetric) {
          // first step: add any un-recoed v0mmcores that were requested
          for (const auto& mcParticle : mcParticles) {
            thisCascInfo.pdgCode = -1, thisCascInfo.pdgCodeMother = -1;
            thisCascInfo.pdgCodePositive = -1, thisCascInfo.pdgCodeNegative = -1;
            thisCascInfo.pdgCodeBachelor = -1, thisCascInfo.pdgCodeV0 = -1;
            thisCascInfo.isPhysicalPrimary = false;
            thisCascInfo.xyz[0] = 0.0f, thisCascInfo.xyz[1] = 0.0f, thisCascInfo.xyz[2] = 0.0f;
            thisCascInfo.lxyz[0] = 0.0f, thisCascInfo.lxyz[1] = 0.0f, thisCascInfo.lxyz[2] = 0.0f;
            thisCascInfo.posP[0] = 0.0f, thisCascInfo.posP[1] = 0.0f, thisCascInfo.posP[2] = 0.0f;
            thisCascInfo.negP[0] = 0.0f, thisCascInfo.negP[1] = 0.0f, thisCascInfo.negP[2] = 0.0f;
            thisCascInfo.bachP[0] = 0.0f, thisCascInfo.bachP[1] = 0.0f, thisCascInfo.bachP[2] = 0.0f;
            thisCascInfo.momentum[0] = 0.0f, thisCascInfo.momentum[1] = 0.0f, thisCascInfo.momentum[2] = 0.0f;
            thisCascInfo.label = -1, thisCascInfo.motherLabel = -1;
            thisCascInfo.mcParticlePositive = -1;
            thisCascInfo.mcParticleNegative = -1;
            thisCascInfo.mcParticleBachelor = -1;

            if (mcParticleIsReco[mcParticle.globalIndex()] == true)
              continue; // skip if already created in list

            if (std::fabs(mcParticle.y()) > cascadeBuilderOpts.mc_rapidityWindow)
              continue; // skip outside midrapidity

            if (cascadeBuilderOpts.mc_keepOnlyPhysicalPrimary && !mcParticle.isPhysicalPrimary())
              continue; // skip secondary MC cascades

            if (
              (cascadeBuilderOpts.mc_addGeneratedXiMinus && mcParticle.pdgCode() == PDG_t::kXiMinus) ||
              (cascadeBuilderOpts.mc_addGeneratedXiPlus && mcParticle.pdgCode() == PDG_t::kXiPlusBar) ||
              (cascadeBuilderOpts.mc_addGeneratedOmegaMinus && mcParticle.pdgCode() == PDG_t::kOmegaMinus) ||
              (cascadeBuilderOpts.mc_addGeneratedOmegaPlus && mcParticle.pdgCode() == PDG_t::kOmegaPlusBar)) {
              thisCascInfo.pdgCode = mcParticle.pdgCode();
              thisCascInfo.isPhysicalPrimary = mcParticle.isPhysicalPrimary();

              if (mcParticle.has_mcCollision()) {
                thisCascInfo.mcCollision = mcParticle.mcCollisionId(); // save this reference, please
              }
              thisCascInfo.momentum[0] = mcParticle.px();
              thisCascInfo.momentum[1] = mcParticle.py();
              thisCascInfo.momentum[2] = mcParticle.pz();
              thisCascInfo.label = mcParticle.globalIndex();

              if (mcParticle.has_daughters()) {
                auto const& daughters = mcParticle.template daughters_as<aod::McParticles>();
                for (const auto& dau : daughters) {
                  if (dau.getProcess() != TMCProcess::kPDecay) // check whether the daughter comes from a decay
                    continue;

                  if (std::abs(dau.pdgCode()) == PDG_t::kPiPlus || std::abs(dau.pdgCode()) == PDG_t::kKPlus) {
                    thisCascInfo.pdgCodeBachelor = dau.pdgCode();
                    thisCascInfo.bachP[0] = dau.px();
                    thisCascInfo.bachP[1] = dau.py();
                    thisCascInfo.bachP[2] = dau.pz();
                    thisCascInfo.xyz[0] = dau.vx();
                    thisCascInfo.xyz[1] = dau.vy();
                    thisCascInfo.xyz[2] = dau.vz();
                    thisCascInfo.mcParticleBachelor = dau.globalIndex();
                  }
                  if (std::abs(dau.pdgCode()) == PDG_t::kProton) {
                    thisCascInfo.pdgCodeV0 = dau.pdgCode();

                    for (const auto& v0Dau : dau.template daughters_as<aod::McParticles>()) {
                      if (v0Dau.getProcess() != TMCProcess::kPDecay)
                        continue;

                      if (v0Dau.pdgCode() > 0) {
                        thisCascInfo.pdgCodePositive = v0Dau.pdgCode();
                        thisCascInfo.processPositive = v0Dau.getProcess();
                        thisCascInfo.posP[0] = v0Dau.px();
                        thisCascInfo.posP[1] = v0Dau.py();
                        thisCascInfo.posP[2] = v0Dau.pz();
                        thisCascInfo.lxyz[0] = v0Dau.vx();
                        thisCascInfo.lxyz[1] = v0Dau.vy();
                        thisCascInfo.lxyz[2] = v0Dau.vz();
                        thisCascInfo.mcParticlePositive = v0Dau.globalIndex();
                      }
                      if (v0Dau.pdgCode() < 0) {
                        thisCascInfo.pdgCodeNegative = v0Dau.pdgCode();
                        thisCascInfo.processNegative = v0Dau.getProcess();
                        thisCascInfo.negP[0] = v0Dau.px();
                        thisCascInfo.negP[1] = v0Dau.py();
                        thisCascInfo.negP[2] = v0Dau.pz();
                        thisCascInfo.mcParticleNegative = v0Dau.globalIndex();
                      }
                    }
                  }
                }
              }

              // if I got here, it means this MC particle was not recoed and is of interest. Add it please
              mcCascinfos.push_back(thisCascInfo);
            }
          }

          for (const auto& thisInfoToFill : mcCascinfos) {
            if (baseOpts.mEnabledTables[kCascMCCores]) {
              products.cascmccores( // a lot of the info below will be compressed in case of not-recoed MC (good!)
                thisInfoToFill.pdgCode, thisInfoToFill.pdgCodeMother, thisInfoToFill.pdgCodeV0, thisInfoToFill.isPhysicalPrimary,
                thisInfoToFill.pdgCodePositive, thisInfoToFill.pdgCodeNegative, thisInfoToFill.pdgCodeBachelor,
                thisInfoToFill.xyz[0], thisInfoToFill.xyz[1], thisInfoToFill.xyz[2],
                thisInfoToFill.lxyz[0], thisInfoToFill.lxyz[1], thisInfoToFill.lxyz[2],
                thisInfoToFill.posP[0], thisInfoToFill.posP[1], thisInfoToFill.posP[2],
                thisInfoToFill.negP[0], thisInfoToFill.negP[1], thisInfoToFill.negP[2],
                thisInfoToFill.bachP[0], thisInfoToFill.bachP[1], thisInfoToFill.bachP[2],
                thisInfoToFill.momentum[0], thisInfoToFill.momentum[1], thisInfoToFill.momentum[2]);
              histos.fill(HIST("hTableBuildingStatistics"), kCascMCCores);
            }
            if (baseOpts.mEnabledTables[kCascMCCollRefs]) {
              products.cascmccollrefs(thisInfoToFill.mcCollision);
              histos.fill(HIST("hTableBuildingStatistics"), kCascMCCollRefs);
            }
          }
        }
      } // enabled tables check
    } // constexpr requires mcParticles check

    LOGF(debug, "Cascades in DF: %i, cascades built: %i", cascades.size(), nCascades);
  }

  //__________________________________________________
  template <typename THistoRegistry, typename TCollisions, typename TCascades, typename TTracks, typename TMCParticles, typename TProducts>
  void buildKFCascades(THistoRegistry& histos, TCollisions const& collisions, TCascades const& cascades, TTracks const& tracks, TMCParticles const& mcParticles, TProducts& products)
  {
    if (!baseOpts.mEnabledTables[kStoredKFCascCores]) {
      return; // don't do if no request for cascades in place
    }
    int nCascades = 0;
    // Loops over all cascades in the time frame
    histos.fill(HIST("hInputStatistics"), kStoredKFCascCores, cascades.size());
    for (size_t icascade = 0; icascade < cascades.size(); icascade++) {
      // Get tracks and generate candidate
      auto const& cascade = cascades[sorted_cascade[icascade]];
      // if collisionId positive: get vertex, negative: origin
      // could be replaced by mean vertex (but without much benefit...)
      float pvX = 0.0f, pvY = 0.0f, pvZ = 0.0f;
      if (cascade.collisionId >= 0) {
        auto const& collision = collisions.rawIteratorAt(cascade.collisionId);
        pvX = collision.posX();
        pvY = collision.posY();
        pvZ = collision.posZ();
      }
      auto const& posTrack = tracks.rawIteratorAt(cascade.posTrackId);
      auto const& negTrack = tracks.rawIteratorAt(cascade.negTrackId);
      auto const& bachTrack = tracks.rawIteratorAt(cascade.bachTrackId);
      if (!straHelper.buildCascadeCandidateWithKF(cascade.collisionId, pvX, pvY, pvZ,
                                                  posTrack,
                                                  negTrack,
                                                  bachTrack,
                                                  baseOpts.mEnabledTables[kCascBBs],
                                                  cascadeBuilderOpts.kfConstructMethod,
                                                  cascadeBuilderOpts.kfTuneForOmega,
                                                  cascadeBuilderOpts.kfUseV0MassConstraint,
                                                  cascadeBuilderOpts.kfUseCascadeMassConstraint,
                                                  cascadeBuilderOpts.kfDoDCAFitterPreMinimV0,
                                                  cascadeBuilderOpts.kfDoDCAFitterPreMinimCasc)) {
        products.kfcascdataLink(-1);
        interlinks.cascadeToKFCascCores.push_back(-1);
        continue; // didn't work out, skip
      }
      nCascades++;

      // generate analysis tables as required
      if (baseOpts.mEnabledTables[kKFCascIndices]) {
        products.kfcascidx(cascade.globalId,
                           straHelper.cascade.positiveTrack, straHelper.cascade.negativeTrack,
                           straHelper.cascade.bachelorTrack, straHelper.cascade.collisionId);
        histos.fill(HIST("hTableBuildingStatistics"), kKFCascIndices);
      }
      if (baseOpts.mEnabledTables[kStoredKFCascCores]) {
        products.kfcascdata(straHelper.cascade.charge, straHelper.cascade.massXi, straHelper.cascade.massOmega,
                            straHelper.cascade.cascadePosition[0], straHelper.cascade.cascadePosition[1], straHelper.cascade.cascadePosition[2],
                            straHelper.cascade.v0Position[0], straHelper.cascade.v0Position[1], straHelper.cascade.v0Position[2],
                            straHelper.cascade.positivePosition[0], straHelper.cascade.positivePosition[1], straHelper.cascade.positivePosition[2],
                            straHelper.cascade.negativePosition[0], straHelper.cascade.negativePosition[1], straHelper.cascade.negativePosition[2],
                            straHelper.cascade.positiveMomentum[0], straHelper.cascade.positiveMomentum[1], straHelper.cascade.positiveMomentum[2],
                            straHelper.cascade.negativeMomentum[0], straHelper.cascade.negativeMomentum[1], straHelper.cascade.negativeMomentum[2],
                            straHelper.cascade.bachelorMomentum[0], straHelper.cascade.bachelorMomentum[1], straHelper.cascade.bachelorMomentum[2],
                            straHelper.cascade.v0Momentum[0], straHelper.cascade.v0Momentum[1], straHelper.cascade.v0Momentum[2],
                            straHelper.cascade.cascadeMomentum[0], straHelper.cascade.cascadeMomentum[1], straHelper.cascade.cascadeMomentum[2],
                            straHelper.cascade.v0DaughterDCA, straHelper.cascade.cascadeDaughterDCA,
                            straHelper.cascade.positiveDCAxy, straHelper.cascade.negativeDCAxy,
                            straHelper.cascade.bachelorDCAxy, straHelper.cascade.cascadeDCAxy, straHelper.cascade.cascadeDCAz,
                            straHelper.cascade.kfMLambda, straHelper.cascade.kfV0Chi2, straHelper.cascade.kfCascadeChi2);
        histos.fill(HIST("hTableBuildingStatistics"), kStoredKFCascCores);

        // interlink always produced if cascades generated
        products.kfcascdataLink(products.kfcascdata.lastIndex());
        interlinks.kfCascCoreToCascades.push_back(cascade.globalId);
        interlinks.cascadeToKFCascCores.push_back(products.kfcascdata.lastIndex());
      }
      if (baseOpts.mEnabledTables[kKFCascCovs]) {
        products.kfcasccovs(straHelper.cascade.covariance, straHelper.cascade.kfTrackCovarianceV0, straHelper.cascade.kfTrackCovariancePos, straHelper.cascade.kfTrackCovarianceNeg);
        histos.fill(HIST("hTableBuildingStatistics"), kKFCascCovs);
      }

      //_________________________________________________________
      // MC handling part (labels only)
      if constexpr (soa::is_table<TMCParticles>) {
        // only worry about this if someone else worried about this
        if ((baseOpts.mEnabledTables[kMcKFCascLabels])) {
          extractMonteCarloProperties(posTrack, negTrack, bachTrack, mcParticles);

          // Construct label table (note: this will be joinable with KFCascDatas)
          products.kfcasclabels(thisCascInfo.label);
          histos.fill(HIST("hTableBuildingStatistics"), kMcKFCascLabels);
        } // enabled tables check
      } // constexpr requires mcParticles check
    } // end loop over cascades

    LOGF(debug, "KF Cascades in DF: %i, KF cascades built: %i", cascades.size(), nCascades);
  }

  //__________________________________________________
  template <class TTracks, typename THistoRegistry, typename TCollisions, typename TStrangeTracks, typename TMCParticles, typename TProducts>
  void buildTrackedCascades(THistoRegistry& histos, TCollisions const& collisions, TStrangeTracks const& cascadeTracks, TMCParticles const& mcParticles, TProducts& products)
  {
    if (!baseOpts.mEnabledTables[kStoredTraCascCores] || baseOpts.mc_findableMode.value != 0) {
      return; // don't do if no request for cascades in place or findable mode used
    }
    int nCascades = 0;
    std::vector<int> traCascIndices(cascadeList.size(), -1);
    // Loops over all V0s in the time frame
    histos.fill(HIST("hInputStatistics"), kStoredTraCascCores, cascadeTracks.size());
    for (const auto& cascadeTrack : cascadeTracks) {
      // Get tracks and generate candidate
      if (!cascadeTrack.has_track()) {
        continue; // safety (should be fine but depends on future stratrack dev)
      }

      auto const& strangeTrack = cascadeTrack.template track_as<TTracks>();

      // if collisionId positive: get vertex, negative: origin
      // could be replaced by mean vertex (but without much benefit...)
      float pvX = 0.0f, pvY = 0.0f, pvZ = 0.0f;
      if (strangeTrack.has_collision()) {
        auto const& collision = collisions.rawIteratorAt(strangeTrack.collisionId());
        pvX = collision.posX();
        pvY = collision.posY();
        pvZ = collision.posZ();
      }
      auto const& cascade = cascadeTrack.cascade();
      auto const& v0 = cascade.v0();
      auto const& posTrack = v0.template posTrack_as<TTracks>();
      auto const& negTrack = v0.template negTrack_as<TTracks>();
      auto const& bachTrack = cascade.template bachelor_as<TTracks>();
      if (!straHelper.buildCascadeCandidate(strangeTrack.collisionId(), pvX, pvY, pvZ,
                                            posTrack,
                                            negTrack,
                                            bachTrack,
                                            baseOpts.mEnabledTables[kCascBBs],
                                            cascadeBuilderOpts.useCascadeMomentumAtPrimVtx,
                                            baseOpts.mEnabledTables[kCascCovs])) {
        continue; // didn't work out, skip
      }

      // recalculate DCAxy, DCAz with strange track
      auto strangeTrackParCov = getTrackParCov(strangeTrack);
      std::array<float, 2> dcaInfo;
      strangeTrackParCov.setPID(o2::track::PID::XiMinus); // FIXME: not OK for omegas
      o2::base::Propagator::Instance()->propagateToDCABxByBz({pvX, pvY, pvZ}, strangeTrackParCov, 2.f, straHelper.fitter.getMatCorrType(), &dcaInfo);
      straHelper.cascade.cascadeDCAxy = dcaInfo[0];
      straHelper.cascade.cascadeDCAz = dcaInfo[1];

      // get momentum from strange track (should not be very different)
      strangeTrackParCov.getPxPyPzGlo(straHelper.cascade.cascadeMomentum);

      // accounting
      nCascades++;

      // generate analysis tables as required
      if (baseOpts.mEnabledTables[kTraCascIndices]) {
        products.tracascidx(cascade.globalIndex(),
                            straHelper.cascade.positiveTrack, straHelper.cascade.negativeTrack,
                            straHelper.cascade.bachelorTrack, cascadeTrack.trackId(), straHelper.cascade.collisionId);
        histos.fill(HIST("hTableBuildingStatistics"), kTraCascIndices);
      }
      if (baseOpts.mEnabledTables[kStoredTraCascCores]) {
        products.tracascdata(straHelper.cascade.charge, cascadeTrack.xiMass(), cascadeTrack.omegaMass(),
                             cascadeTrack.decayX(), cascadeTrack.decayY(), cascadeTrack.decayZ(),
                             straHelper.cascade.v0Position[0], straHelper.cascade.v0Position[1], straHelper.cascade.v0Position[2],
                             straHelper.cascade.positiveMomentum[0], straHelper.cascade.positiveMomentum[1], straHelper.cascade.positiveMomentum[2],
                             straHelper.cascade.negativeMomentum[0], straHelper.cascade.negativeMomentum[1], straHelper.cascade.negativeMomentum[2],
                             straHelper.cascade.bachelorMomentum[0], straHelper.cascade.bachelorMomentum[1], straHelper.cascade.bachelorMomentum[2],
                             straHelper.cascade.cascadeMomentum[0], straHelper.cascade.cascadeMomentum[1], straHelper.cascade.cascadeMomentum[2],
                             straHelper.cascade.v0DaughterDCA, straHelper.cascade.cascadeDaughterDCA,
                             straHelper.cascade.positiveDCAxy, straHelper.cascade.negativeDCAxy,
                             straHelper.cascade.bachelorDCAxy, straHelper.cascade.cascadeDCAxy, straHelper.cascade.cascadeDCAz,
                             cascadeTrack.matchingChi2(), cascadeTrack.topologyChi2(), cascadeTrack.itsClsSize());
        histos.fill(HIST("hTableBuildingStatistics"), kStoredTraCascCores);

        // interlink always produced if base core table generated
        traCascIndices[cascade.globalIndex()] = products.tracascdata.lastIndex();
      }
      if (baseOpts.mEnabledTables[kCascCovs]) {
        std::array<float, o2::track::kLabCovMatSize> traCovMat = {0.};
        strangeTrackParCov.getCovXYZPxPyPzGlo(traCovMat);
        float traCovMatArray[o2::track::kLabCovMatSize];
        for (int ii = 0; ii < o2::track::kLabCovMatSize; ii++) {
          traCovMatArray[ii] = traCovMat[ii];
        }
        products.tracasccovs(traCovMatArray);
        histos.fill(HIST("hTableBuildingStatistics"), kCascCovs);
      }

      //_________________________________________________________
      // MC handling part (labels only)
      if constexpr (soa::is_table<TMCParticles>) {
        // only worry about this if someone else worried about this
        if ((baseOpts.mEnabledTables[kMcTraCascLabels])) {
          extractMonteCarloProperties(posTrack, negTrack, bachTrack, mcParticles);

          // Construct label table (note: this will be joinable with KFCascDatas)
          products.tracasclabels(thisCascInfo.label);
          histos.fill(HIST("hTableBuildingStatistics"), kMcTraCascLabels);
        } // enabled tables check
      } // constexpr requires mcParticles check
    } // end loop over cascades

    for (std::size_t icascade = 0; icascade < cascadeList.size(); icascade++) {
      products.tracascdataLink(traCascIndices[icascade]);
    }
    LOGF(debug, "Tracked cascades in DF: %i, tracked cascades built: %i", cascadeTracks.size(), nCascades);
  }

  //__________________________________________________
  // MC kink handling
  template <typename mcpart>
  int getOriginatingParticle(mcpart const& part, int& indexForPositionOfDecay, bool treatPiToMuDecays)
  {
    int returnValue = -1;
    if (part.has_mothers()) {
      auto const& motherList = part.template mothers_as<aod::McParticles>();
      if (motherList.size() == 1) {
        for (const auto& mother : motherList) {
          if (std::abs(part.pdgCode()) == PDG_t::kMuonMinus && treatPiToMuDecays) {
            // muon decay, de-ref mother twice
            if (mother.has_mothers()) {
              auto grandMotherList = mother.template mothers_as<aod::McParticles>();
              if (grandMotherList.size() == 1) {
                for (const auto& grandMother : grandMotherList) {
                  returnValue = grandMother.globalIndex();
                  indexForPositionOfDecay = mother.globalIndex(); // for V0 decay position: grab muon
                }
              }
            }
          } else {
            returnValue = mother.globalIndex();
            indexForPositionOfDecay = part.globalIndex();
          }
        }
      }
    }
    return returnValue;
  }

  //__________________________________________________
  template <typename TCCDB, typename THistoRegistry, typename TCollisions, typename TMCCollisions, typename TV0s, typename TCascades, typename TTrackedCascades, typename TTracks, typename TBCs, typename TMCParticles, typename TProducts>
  void dataProcess(TCCDB& ccdb, THistoRegistry& histos, TCollisions const& collisions, TMCCollisions const& mccollisions, TV0s const& v0s, TCascades const& cascades, TTrackedCascades const& trackedCascades, TTracks const& tracks, TBCs const& bcs, TMCParticles const& mcParticles, TProducts& products)
  {
    if (nEnabledTables == 0) {
      return; // fully suppressed
    }

    if (!initCCDB(ccdb, bcs, collisions))
      return;

    // reset vectors for cascade interlinks
    resetInterlinks();

    // prepare v0List, cascadeList
    prepareBuildingLists<TBCs>(histos, collisions, mccollisions, v0s, cascades, tracks, mcParticles);

    // mark V0s that will be buffered for the cascade building
    markV0sUsedInCascades(v0List, cascadeList, trackedCascades);

    // build V0s
    buildV0s<TBCs>(histos, collisions, v0List, tracks, mcParticles, products);

    // build cascades
    buildCascades(histos, collisions, cascadeList, tracks, mcParticles, products);
    buildKFCascades(histos, collisions, cascadeList, tracks, mcParticles, products);

    // build tracked cascades only if subscription is Run 3 like (doesn't exist in Run 2)
    if constexpr (soa::is_table<TTrackedCascades>) {
      buildTrackedCascades<TTracks>(histos, collisions, trackedCascades, mcParticles, products);
    }

    populateCascadeInterlinks();
  }
}; // end BuilderModule

} // namespace strangenessbuilder
} // namespace pwglf
} // namespace o2

#endif // PWGLF_UTILS_STRANGENESSBUILDERMODULE_H_
