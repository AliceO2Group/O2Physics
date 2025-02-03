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
// Strangeness builder task
// ========================
//
// This task produces all tables that may be necessary for
// strangeness analyses. A single device is provided to
// ensure better computing resource (memory) management.
//
//  process functions:
//
//  -- processRealData[Run2] .........: use this OR processMonteCarlo but NOT both
//  -- processMonteCarlo[Run2] .......: use this OR processRealData but NOT both
//

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Common/DataModel/PIDResponse.h"
#include "TableHelper.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/LFStrangenessPIDTables.h"
#include "PWGLF/Utils/strangenessBuilderHelper.h"
#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"

using namespace o2;
using namespace o2::framework;

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
  "KFToCascRefs"        //.34 (interlink KFCascCores -> CascCores)
};

static constexpr int nTablesConst = 35;

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
  {-1}};

// use parameters + cov mat non-propagated, aux info + (extension propagated)
using FullTracksExt = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCov>;
using FullTracksExtIU = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU>;
using FullTracksExtLabeled = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCov, aod::McTrackLabels>;
using FullTracksExtLabeledIU = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::McTrackLabels>;
using TracksWithExtra = soa::Join<aod::Tracks, aod::TracksExtra>;

// For dE/dx association in pre-selection
using TracksExtraWithPID = soa::Join<aod::TracksExtra, aod::pidTPCFullEl, aod::pidTPCFullPi, aod::pidTPCFullPr, aod::pidTPCFullHe>;

struct StrangenessBuilder {
  // helper object
  o2::pwglf::strangenessBuilderHelper straHelper;

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
                    nTables };

  //__________________________________________________
  // V0 tables
  Produces<aod::V0Indices> v0indices; // standard part of V0Datas
  Produces<aod::V0CoresBase> v0cores; // standard part of V0Datas
  Produces<aod::V0Covs> v0covs;       // for decay chain reco

  //__________________________________________________
  // cascade tables
  Produces<aod::CascIndices> cascidx;            // standard part of CascDatas
  Produces<aod::KFCascIndices> kfcascidx;        // standard part of KFCascDatas
  Produces<aod::TraCascIndices> tracascidx;      // standard part of TraCascDatas
  Produces<aod::StoredCascCores> cascdata;       // standard part of CascDatas
  Produces<aod::StoredKFCascCores> kfcascdata;   // standard part of KFCascDatas
  Produces<aod::StoredTraCascCores> tracascdata; // standard part of TraCascDatas
  Produces<aod::CascCovs> casccovs;              // for decay chain reco
  Produces<aod::KFCascCovs> kfcasccovs;          // for decay chain reco
  Produces<aod::TraCascCovs> tracasccovs;        // for decay chain reco

  //__________________________________________________
  // interlink tables
  Produces<aod::V0DataLink> v0dataLink;           // de-refs V0s -> V0Data
  Produces<aod::CascDataLink> cascdataLink;       // de-refs Cascades -> CascData
  Produces<aod::KFCascDataLink> kfcascdataLink;   // de-refs Cascades -> KFCascData
  Produces<aod::TraCascDataLink> tracascdataLink; // de-refs Cascades -> TraCascData

  //__________________________________________________
  // secondary auxiliary tables
  Produces<aod::V0TrackXs> v0trackXs;     // for decay chain reco
  Produces<aod::CascTrackXs> cascTrackXs; // for decay chain reco

  //__________________________________________________
  // further auxiliary / optional if desired
  Produces<aod::CascBBs> cascbb;
  Produces<aod::V0DauCovs> v0daucovs;            // covariances of daughter tracks
  Produces<aod::V0DauCovIUs> v0daucovIUs;        // covariances of daughter tracks
  Produces<aod::V0TraPosAtDCAs> v0dauPositions;  // auxiliary debug information
  Produces<aod::V0TraPosAtIUs> v0dauPositionsIU; // auxiliary debug information
  Produces<aod::V0Ivanovs> v0ivanovs;            // information for Marian's tests

  //__________________________________________________
  // MC information: V0
  Produces<aod::McV0Labels> v0labels;           // MC labels for V0s
  Produces<aod::V0MCCores> v0mccores;           // mc info storage
  Produces<aod::V0CoreMCLabels> v0CoreMCLabels; // interlink V0Cores -> V0MCCores
  Produces<aod::V0MCCollRefs> v0mccollref;      // references collisions from V0MCCores

  // MC information: Cascades
  Produces<aod::McCascLabels> casclabels;           // MC labels for cascades
  Produces<aod::McKFCascLabels> kfcasclabels;       // MC labels for KF cascades
  Produces<aod::McTraCascLabels> tracasclabels;     // MC labels for tracked cascades
  Produces<aod::McCascBBTags> bbtags;               // bb tags (inv structure tagging in mc)
  Produces<aod::CascMCCores> cascmccores;           // mc info storage
  Produces<aod::CascCoreMCLabels> cascCoreMClabels; // interlink CascCores -> CascMCCores
  Produces<aod::CascMCCollRefs> cascmccollrefs;     // references MC collisions from MC cascades

  //__________________________________________________
  // cascade interlinks
  Produces<aod::CascToTraRefs> cascToTraRefs; // cascades -> tracked
  Produces<aod::CascToKFRefs> cascToKFRefs;   // cascades -> KF
  Produces<aod::TraToCascRefs> traToCascRefs; // tracked -> cascades
  Produces<aod::KFToCascRefs> kfToCascRefs;   // KF -> cascades

  Configurable<LabeledArray<int>> enabledTables{"enabledTables",
                                                {defaultParameters[0], nTables, nParameters, tableNames, parameterNames},
                                                "Produce this table: -1 for autodetect; otherwise, 0/1 is false/true"};
  std::vector<int> mEnabledTables; // Vector of enabled tables

  // CCDB options
  struct : ConfigurableGroup {
    std::string prefix = "ccdb";
    Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
    Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
    Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
    Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
    Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
  } ccdbConfigurations;

  // V0 building options
  struct : ConfigurableGroup {
    std::string prefix = "v0BuilderOpts";
    Configurable<bool> generatePhotonCandidates{"generatePhotonCandidates", false, "generate gamma conversion candidates (V0s using TPC-only tracks)"};

    // baseline conditionals of V0 building
    Configurable<int> minCrossedRows{"minCrossedRows", 50, "minimum TPC crossed rows for daughter tracks"};
    Configurable<float> dcanegtopv{"dcanegtopv", .1, "DCA Neg To PV"};
    Configurable<float> dcapostopv{"dcapostopv", .1, "DCA Pos To PV"};
    Configurable<double> v0cospa{"v0cospa", 0.95, "V0 CosPA"}; // double -> N.B. dcos(x)/dx = 0 at x=0)
    Configurable<float> dcav0dau{"dcav0dau", 1.0, "DCA V0 Daughters"};
    Configurable<float> v0radius{"v0radius", 0.9, "v0radius"};
    Configurable<float> maxDaughterEta{"maxDaughterEta", 5.0, "Maximum daughter eta (in abs value)"};

    // MC builder options
    Configurable<bool> mc_populateV0MCCoresSymmetric{"mc_populateV0MCCoresSymmetric", false, "populate V0MCCores table for derived data analysis, keep V0MCCores joinable with V0Cores"};
    Configurable<bool> mc_populateV0MCCoresAsymmetric{"mc_populateV0MCCoresAsymmetric", true, "populate V0MCCores table for derived data analysis, create V0Cores -> V0MCCores interlink. Saves only labeled V0s."};
    Configurable<bool> mc_treatPiToMuDecays{"mc_treatPiToMuDecays", true, "if true, will correctly capture pi -> mu and V0 label will still point to originating V0 decay in those cases. Nota bene: prong info will still be for the muon!"};
    Configurable<float> mc_rapidityWindow{"mc_rapidityWindow", 0.5, "rapidity window to save non-recoed candidates"};
    Configurable<bool> mc_addGeneratedK0Short{"mc_addGeneratedK0Short", false, "add V0MCCore entry for generated, not-recoed K0Short"};
    Configurable<bool> mc_addGeneratedLambda{"mc_addGeneratedLambda", false, "add V0MCCore entry for generated, not-recoed Lambda"};
    Configurable<bool> mc_addGeneratedAntiLambda{"mc_addGeneratedAntiLambda", false, "add V0MCCore entry for generated, not-recoed AntiLambda"};
    Configurable<bool> mc_addGeneratedGamma{"mc_addGeneratedGamma", false, "add V0MCCore entry for generated, not-recoed Gamma"};
  } v0BuilderOpts;

  // cascade building options
  struct : ConfigurableGroup {
    std::string prefix = "cascadeBuilderOpts";
    Configurable<bool> useCascadeMomentumAtPrimVtx{"useCascadeMomentumAtPrimVtx", false, "use cascade momentum at PV"};

    // conditionals
    Configurable<int> minCrossedRows{"minCrossedRows", 50, "minimum TPC crossed rows for daughter tracks"};
    Configurable<float> dcabachtopv{"dcabachtopv", .05, "DCA Bach To PV"};
    Configurable<float> cascradius{"cascradius", 0.9, "cascradius"};
    Configurable<float> casccospa{"casccospa", 0.95, "casccospa"};
    Configurable<float> dcacascdau{"dcacascdau", 1.0, "DCA cascade Daughters"};
    Configurable<float> lambdaMassWindow{"lambdaMassWindow", .010, "Distance from Lambda mass (does not apply to KF path)"};
    Configurable<float> maxDaughterEta{"maxDaughterEta", 5.0, "Maximum daughter eta (in abs value)"};

    // KF building specific
    Configurable<bool> kfTuneForOmega{"kfTuneForOmega", false, "if enabled, take main cascade properties from Omega fit instead of Xi fit (= default)"};
    Configurable<int> kfConstructMethod{"kfConstructMethod", 2, "KF Construct Method"};
    Configurable<bool> kfUseV0MassConstraint{"kfUseV0MassConstraint", true, "KF: use Lambda mass constraint"};
    Configurable<bool> kfUseCascadeMassConstraint{"kfUseCascadeMassConstraint", false, "KF: use Cascade mass constraint - WARNING: not adequate for inv mass analysis of Xi"};
    Configurable<bool> kfDoDCAFitterPreMinimV0{"kfDoDCAFitterPreMinimV0", true, "KF: do DCAFitter pre-optimization before KF fit to include material corrections for V0"};
    Configurable<bool> kfDoDCAFitterPreMinimCasc{"kfDoDCAFitterPreMinimCasc", true, "KF: do DCAFitter pre-optimization before KF fit to include material corrections for Xi"};

    // MC builder options
    Configurable<bool> mc_populateCascMCCoresSymmetric{"mc_populateCascMCCoresSymmetric", false, "populate CascMCCores table for derived data analysis, keep CascMCCores joinable with CascCores"};
    Configurable<bool> mc_populateCascMCCoresAsymmetric{"mc_populateCascMCCoresAsymmetric", true, "populate CascMCCores table for derived data analysis, create CascCores -> CascMCCores interlink. Saves only labeled Cascades."};
    Configurable<bool> mc_addGeneratedXiMinus{"mc_addGeneratedXiMinus", false, "add CascMCCore entry for generated, not-recoed XiMinus"};
    Configurable<bool> mc_addGeneratedXiPlus{"mc_addGeneratedXiPlus", false, "add CascMCCore entry for generated, not-recoed XiPlus"};
    Configurable<bool> mc_addGeneratedOmegaMinus{"mc_addGeneratedOmegaMinus", false, "add CascMCCore entry for generated, not-recoed OmegaMinus"};
    Configurable<bool> mc_addGeneratedOmegaPlus{"mc_addGeneratedOmegaPlus", false, "add CascMCCore entry for generated, not-recoed OmegaPlus"};
    Configurable<bool> mc_treatPiToMuDecays{"mc_treatPiToMuDecays", true, "if true, will correctly capture pi -> mu and V0 label will still point to originating V0 decay in those cases. Nota bene: prong info will still be for the muon!"};
    Configurable<float> mc_rapidityWindow{"mc_rapidityWindow", 0.5, "rapidity window to save non-recoed candidates"};
  } cascadeBuilderOpts;

  o2::ccdb::CcdbApi ccdbApi;
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  int mRunNumber;
  o2::base::MatLayerCylSet* lut = nullptr;

  // for tagging V0s used in cascades
  std::vector<o2::pwglf::v0candidate> v0sFromCascades; // Vector of v0 candidates used in cascades
  std::vector<int> v0Map;                              // index to relate V0s -> v0sFromCascades

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

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext& context)
  {
    // setup bookkeeping histogram
    auto h = histos.add<TH1>("hTableBuildingStatistics", "hTableBuildingStatistics", kTH1D, {{nTablesConst, -0.5f, static_cast<float>(nTablesConst)}});
    auto h2 = histos.add<TH1>("hInputStatistics", "hInputStatistics", kTH1D, {{nTablesConst, -0.5f, static_cast<float>(nTablesConst)}});
    h2->SetTitle("Input table sizes");

    mRunNumber = 0;

    mEnabledTables.resize(nTables, 0);

    LOGF(info, "Configuring tables to generate");
    auto& workflows = context.services().get<RunningWorkflowInfo const>();

    for (int i = 0; i < nTables; i++) {
      // adjust bookkeeping histogram
      h->GetXaxis()->SetBinLabel(i + 1, tableNames[i].c_str());
      h2->GetXaxis()->SetBinLabel(i + 1, tableNames[i].c_str());
      h->SetBinContent(i + 1, -1); // mark all as disabled to start

      int f = enabledTables->get(tableNames[i].c_str(), "enable");
      if (f == 1) {
        mEnabledTables[i] = 1;
      }
      if (f == -1) {
        // autodetect this table in other devices
        for (DeviceSpec const& device : workflows.devices) {
          // Step 1: check if this device subscribed to the V0data table
          for (auto const& input : device.inputs) {
            if (device.name.compare("strangenessbuilder-initializer") == 0)
              continue; // don't listen to the initializer
            if (input.matcher.binding == tableNames[i]) {
              LOGF(info, "Device %s has subscribed to %s", device.name, tableNames[i]);
              mEnabledTables[i] = 1;
            }
          }
        }
      }
    }

    // list enabled tables
    LOGF(info, "*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*");
    LOGF(info, " Strangeness builder: enabled table listing");
    LOGF(info, "*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*");
    for (int i = 0; i < nTables; i++) {
      // printout to be improved in the future
      if (mEnabledTables[i]) {
        LOGF(info, " -~> Table enabled: %s", tableNames[i]);
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

    ccdb->setURL(ccdbConfigurations.ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);

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
  }

  bool initCCDB(aod::BCsWithTimestamps const& bcs, aod::Collisions const& collisions)
  {
    auto bc = collisions.size() ? collisions.begin().bc_as<aod::BCsWithTimestamps>() : bcs.begin();
    if (!bcs.size()) {
      LOGF(warn, "No BC found, skipping this DF.");
      return false; // signal to skip this DF
    }

    if (mRunNumber == bc.runNumber()) {
      return true;
    }

    auto timestamp = bc.timestamp();
    o2::parameters::GRPObject* grpo = 0x0;
    o2::parameters::GRPMagField* grpmag = 0x0;
    if (doprocessRealDataRun2) {
      grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>(ccdbConfigurations.grpPath, timestamp);
      if (!grpo) {
        LOG(fatal) << "Got nullptr from CCDB for path " << ccdbConfigurations.grpPath << " of object GRPObject for timestamp " << timestamp;
      }
      o2::base::Propagator::initFieldFromGRP(grpo);
    } else {
      grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(ccdbConfigurations.grpmagPath, timestamp);
      if (!grpmag) {
        LOG(fatal) << "Got nullptr from CCDB for path " << ccdbConfigurations.grpmagPath << " of object GRPMagField for timestamp " << timestamp;
      }
      o2::base::Propagator::initFieldFromGRP(grpmag);
    }
    // Fetch magnetic field from ccdb for current collision
    auto magneticField = o2::base::Propagator::Instance()->getNominalBz();
    LOG(info) << "Retrieved GRP for timestamp " << timestamp << " with magnetic field of " << magneticField << " kG";

    // Set magnetic field value once known
    straHelper.fitter.setBz(magneticField);

    // acquire LUT for this timestamp
    LOG(info) << "Loading material look-up table for timestamp: " << timestamp;
    lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->getForTimeStamp<o2::base::MatLayerCylSet>(ccdbConfigurations.lutPath, timestamp));
    o2::base::Propagator::Instance()->setMatLUT(lut);
    straHelper.lut = lut;

    LOG(info) << "Fully configured for run: " << bc.runNumber();
    // mmark this run as configured
    mRunNumber = bc.runNumber();

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
    if (mEnabledTables[kCascToKFRefs]) {
      for (auto& cascCore : interlinks.cascCoreToCascades) {
        cascToKFRefs(interlinks.cascadeToKFCascCores[cascCore]);
        histos.fill(HIST("hTableBuildingStatistics"), kCascToKFRefs);
      }
    }
    if (mEnabledTables[kCascToTraRefs]) {
      for (auto& cascCore : interlinks.cascCoreToCascades) {
        cascToTraRefs(interlinks.cascadeToTraCascCores[cascCore]);
        histos.fill(HIST("hTableBuildingStatistics"), kCascToTraRefs);
      }
    }
    if (mEnabledTables[kKFToCascRefs]) {
      for (auto& kfCascCore : interlinks.kfCascCoreToCascades) {
        kfToCascRefs(interlinks.cascadeToCascCores[kfCascCore]);
        histos.fill(HIST("hTableBuildingStatistics"), kKFToCascRefs);
      }
    }
    if (mEnabledTables[kTraToCascRefs]) {
      for (auto& traCascCore : interlinks.traCascCoreToCascades) {
        traToCascRefs(interlinks.cascadeToCascCores[traCascCore]);
        histos.fill(HIST("hTableBuildingStatistics"), kTraToCascRefs);
      }
    }
  }

  //__________________________________________________
  template <typename TV0s, typename TCascades, typename TTrackedCascades>
  void markV0sUsedInCascades(TV0s const& v0s, TCascades const& cascades, TTrackedCascades const& trackedCascades)
  {
    int v0sUsedInCascades = 0;
    v0sFromCascades.clear();
    v0Map.clear();
    v0Map.resize(v0s.size(), -2); // marks not used
    if (mEnabledTables[kStoredCascCores]) {
      for (auto& cascade : cascades) {
        if (v0Map[cascade.v0Id()] == -2) {
          v0sUsedInCascades++;
        }
        v0Map[cascade.v0Id()] = -1; // marks used (but isn't the index of a properly built V0, which would be >= 0)
      }
    }
    int trackedCascadeCount = 0;
    if constexpr (soa::is_table<TTrackedCascades>) {
      if (mEnabledTables[kStoredTraCascCores]) {
        trackedCascadeCount = trackedCascades.size();
        for (auto& trackedCascade : trackedCascades) {
          auto const& cascade = trackedCascade.cascade();
          if (v0Map[cascade.v0Id()] == -2) {
            v0sUsedInCascades++;
          }
          v0Map[cascade.v0Id()] = -1; // marks used (but isn't the index of a properly built V0, which would be >= 0)
        }
      }
    }
    LOGF(debug, "V0 total %i, Cascade total %i, Tracked cascade total %i, V0s flagged used in cascades: %i", v0s.size(), cascades.size(), trackedCascadeCount, v0sUsedInCascades);
  }

  //__________________________________________________
  template <class TTracks, typename TV0s, typename TMCParticles>
  void buildV0s(TV0s const& v0s, TMCParticles const& mcParticles)
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
    for (auto& v0 : v0s) {
      if (!mEnabledTables[kV0CoresBase] && v0Map[v0.globalIndex()] == -2) {
        // this v0 hasn't been used by cascades and we're not generating V0s, so skip it
        v0dataLink(-1, -1);
        continue;
      }

      // Get tracks and generate candidate
      auto const& collision = v0.collision();
      auto const& posTrack = v0.template posTrack_as<TTracks>();
      auto const& negTrack = v0.template negTrack_as<TTracks>();
      if (!straHelper.buildV0Candidate(collision, posTrack, negTrack, v0.isCollinearV0(), mEnabledTables[kV0Covs])) {
        v0dataLink(-1, -1);
        continue;
      }
      if (v0Map[v0.globalIndex()] == -1) {
        v0Map[v0.globalIndex()] = v0sFromCascades.size(); // provide actual valid index in buffer
        v0sFromCascades.push_back(straHelper.v0);
      }
      // fill requested cursors only if type is not 0
      if (v0.v0Type() == 1 || (v0.v0Type() > 1 && v0BuilderOpts.generatePhotonCandidates)) {
        nV0s++;
        if (mEnabledTables[kV0Indices]) {
          // for referencing (especially - but not only - when using derived data)
          v0indices(v0.posTrackId(), v0.negTrackId(),
                    v0.collisionId(), v0.globalIndex());
          histos.fill(HIST("hTableBuildingStatistics"), kV0Indices);
        }
        if (mEnabledTables[kV0TrackXs]) {
          // further decay chains may need this
          v0trackXs(straHelper.v0.positiveTrackX, straHelper.v0.negativeTrackX);
          histos.fill(HIST("hTableBuildingStatistics"), kV0TrackXs);
        }
        if (mEnabledTables[kV0CoresBase]) {
          // standard analysis
          v0cores(straHelper.v0.position[0], straHelper.v0.position[1], straHelper.v0.position[2],
                  straHelper.v0.positiveMomentum[0], straHelper.v0.positiveMomentum[1], straHelper.v0.positiveMomentum[2],
                  straHelper.v0.negativeMomentum[0], straHelper.v0.negativeMomentum[1], straHelper.v0.negativeMomentum[2],
                  straHelper.v0.daughterDCA,
                  straHelper.v0.positiveDCAxy,
                  straHelper.v0.negativeDCAxy,
                  TMath::Cos(straHelper.v0.pointingAngle),
                  straHelper.v0.dcaXY,
                  v0.v0Type());
          v0dataLink(v0cores.lastIndex(), -1);
          histos.fill(HIST("hTableBuildingStatistics"), kV0CoresBase);
        }
        if (mEnabledTables[kV0TraPosAtDCAs]) {
          // for tracking studies
          v0dauPositions(straHelper.v0.positivePosition[0], straHelper.v0.positivePosition[1], straHelper.v0.positivePosition[2],
                         straHelper.v0.negativePosition[0], straHelper.v0.negativePosition[1], straHelper.v0.negativePosition[2]);
          histos.fill(HIST("hTableBuildingStatistics"), kV0TraPosAtDCAs);
        }
        if (mEnabledTables[kV0TraPosAtIUs]) {
          // for tracking studies
          std::array<float, 3> positivePositionIU;
          std::array<float, 3> negativePositionIU;
          o2::track::TrackPar positiveTrackParam = getTrackPar(posTrack);
          o2::track::TrackPar negativeTrackParam = getTrackPar(negTrack);
          positiveTrackParam.getXYZGlo(positivePositionIU);
          negativeTrackParam.getXYZGlo(negativePositionIU);
          v0dauPositionsIU(positivePositionIU[0], positivePositionIU[1], positivePositionIU[2],
                           negativePositionIU[0], negativePositionIU[1], negativePositionIU[2]);
          histos.fill(HIST("hTableBuildingStatistics"), kV0TraPosAtIUs);
        }
        if (mEnabledTables[kV0Covs]) {
          v0covs(straHelper.v0.positionCovariance, straHelper.v0.momentumCovariance);
          histos.fill(HIST("hTableBuildingStatistics"), kV0Covs);
        }

        //_________________________________________________________
        // MC handling part
        if constexpr (soa::is_table<TMCParticles>) {
          // only worry about this if someone else worried about this
          if ((mEnabledTables[kV0MCCores] || mEnabledTables[kMcV0Labels] || mEnabledTables[kV0MCCollRefs])) {
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
                  for (auto& lV0Mother : originatingV0.template mothers_as<aod::McParticles>()) {
                    thisInfo.pdgCodeMother = lV0Mother.pdgCode();
                    thisInfo.motherLabel = lV0Mother.globalIndex();
                  }
                }
              }

            } // end association check
            // Construct label table (note: this will be joinable with V0Datas!)
            if (mEnabledTables[kMcV0Labels]) {
              v0labels(thisInfo.label, thisInfo.motherLabel);
              histos.fill(HIST("hTableBuildingStatistics"), kMcV0Labels);
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
              if (mEnabledTables[kV0MCCores]) {
                v0mccores(
                  thisInfo.label, thisInfo.pdgCode,
                  thisInfo.pdgCodeMother, thisInfo.pdgCodePositive, thisInfo.pdgCodeNegative,
                  thisInfo.isPhysicalPrimary, thisInfo.xyz[0], thisInfo.xyz[1], thisInfo.xyz[2],
                  thisInfo.posP[0], thisInfo.posP[1], thisInfo.posP[2],
                  thisInfo.negP[0], thisInfo.negP[1], thisInfo.negP[2],
                  thisInfo.momentum[0], thisInfo.momentum[1], thisInfo.momentum[2]);
                histos.fill(HIST("hTableBuildingStatistics"), kV0MCCores);
              }
              if (mEnabledTables[kV0MCCollRefs]) {
                v0mccollref(thisInfo.mcCollision);
                histos.fill(HIST("hTableBuildingStatistics"), kV0MCCollRefs);
              }

              // n.b. placing the interlink index here allows for the writing of
              //      code that is agnostic with respect to the joinability of
              //      V0Cores and V0MCCores (always dereference -> safe)
              if (mEnabledTables[kV0CoreMCLabels]) {
                v0CoreMCLabels(v0.globalIndex()); // interlink index
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
              if (mEnabledTables[kV0CoreMCLabels]) {
                v0CoreMCLabels(thisV0MCCoreIndex); // interlink index
                histos.fill(HIST("hTableBuildingStatistics"), kV0CoreMCLabels);
              }
            }
          } // enabled tables check
        } // constexpr requires check
      }
    }

    // finish populating V0MCCores if in asymmetric mode
    if constexpr (soa::is_table<TMCParticles>) {
      if (v0BuilderOpts.mc_populateV0MCCoresAsymmetric && (mEnabledTables[kV0MCCores] || mEnabledTables[kV0MCCollRefs])) {
        // first step: add any un-recoed v0mmcores that were requested
        for (auto& mcParticle : mcParticles) {
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

          if (TMath::Abs(mcParticle.y()) > v0BuilderOpts.mc_rapidityWindow)
            continue; // skip outside midrapidity

          if (
            (v0BuilderOpts.mc_addGeneratedK0Short && mcParticle.pdgCode() == 310) ||
            (v0BuilderOpts.mc_addGeneratedLambda && mcParticle.pdgCode() == 3122) ||
            (v0BuilderOpts.mc_addGeneratedAntiLambda && mcParticle.pdgCode() == -3122) ||
            (v0BuilderOpts.mc_addGeneratedGamma && mcParticle.pdgCode() == 22)) {
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

              for (auto& dau : daughters) {
                if (dau.getProcess() != 4)
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

        for (auto info : mcV0infos) {
          if (mEnabledTables[kV0MCCores]) {
            v0mccores(
              info.label, info.pdgCode,
              info.pdgCodeMother, info.pdgCodePositive, info.pdgCodeNegative,
              info.isPhysicalPrimary, info.xyz[0], info.xyz[1], info.xyz[2],
              info.posP[0], info.posP[1], info.posP[2],
              info.negP[0], info.negP[1], info.negP[2],
              info.momentum[0], info.momentum[1], info.momentum[2]);
            histos.fill(HIST("hTableBuildingStatistics"), kV0MCCores);
          }
          if (mEnabledTables[kV0MCCollRefs]) {
            v0mccollref(info.mcCollision);
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
          for (auto& lV0Mother : originatingV0.template mothers_as<aod::McParticles>()) {
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
                for (auto& lV0GrandMother : lV0Mother.template mothers_as<aod::McParticles>()) {
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
  template <class TTracks, typename TCascades, typename TMCParticles>
  void buildCascades(TCascades const& cascades, TMCParticles const& mcParticles)
  {
    // prepare MC containers (not necessarily used)
    std::vector<mcCascinfo> mcCascinfos; // V0MCCore information
    std::vector<bool> mcParticleIsReco;

    if constexpr (soa::is_table<TMCParticles>) {
      // do this if provided with a mcParticle table as well
      mcParticleIsReco.resize(mcParticles.size(), false);
    }

    if (!mEnabledTables[kStoredCascCores]) {
      return; // don't do if no request for cascades in place
    }
    int nCascades = 0;
    // Loops over all cascades in the time frame
    histos.fill(HIST("hInputStatistics"), kStoredCascCores, cascades.size());
    for (auto& cascade : cascades) {
      // Get tracks and generate candidate
      auto const& collision = cascade.collision();
      auto const& v0 = cascade.v0();
      auto const& posTrack = v0.template posTrack_as<TTracks>();
      auto const& negTrack = v0.template negTrack_as<TTracks>();
      auto const& bachTrack = cascade.template bachelor_as<TTracks>();
      if (v0Map[v0.globalIndex()] < 0) {
        // this V0 hasn't been stored / cached
        cascdataLink(-1);
        interlinks.cascadeToCascCores.push_back(-1);
        continue; // didn't work out, skip
      }
      if (!straHelper.buildCascadeCandidate(collision,
                                            v0sFromCascades[v0Map[v0.globalIndex()]],
                                            posTrack,
                                            negTrack,
                                            bachTrack,
                                            mEnabledTables[kCascBBs],
                                            cascadeBuilderOpts.useCascadeMomentumAtPrimVtx,
                                            mEnabledTables[kCascCovs])) {
        cascdataLink(-1);
        interlinks.cascadeToCascCores.push_back(-1);
        continue; // didn't work out, skip
      }
      nCascades++;

      // generate analysis tables as required
      if (mEnabledTables[kCascIndices]) {
        cascidx(cascade.globalIndex(),
                straHelper.cascade.positiveTrack, straHelper.cascade.negativeTrack,
                straHelper.cascade.bachelorTrack, straHelper.cascade.collisionId);
        histos.fill(HIST("hTableBuildingStatistics"), kCascIndices);
      }
      if (mEnabledTables[kStoredCascCores]) {
        cascdata(straHelper.cascade.charge, straHelper.cascade.massXi, straHelper.cascade.massOmega,
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
        cascdataLink(cascdata.lastIndex());
        interlinks.cascCoreToCascades.push_back(cascade.globalIndex());
        interlinks.cascadeToCascCores.push_back(cascdata.lastIndex());
      }

      if (mEnabledTables[kCascTrackXs]) {
        cascTrackXs(straHelper.cascade.positiveTrackX, straHelper.cascade.negativeTrackX, straHelper.cascade.bachelorTrackX);
        histos.fill(HIST("hTableBuildingStatistics"), kCascTrackXs);
      }
      if (mEnabledTables[kCascBBs]) {
        cascbb(straHelper.cascade.bachBaryonCosPA, straHelper.cascade.bachBaryonDCAxyToPV);
        histos.fill(HIST("hTableBuildingStatistics"), kCascBBs);
      }
      if (mEnabledTables[kCascCovs]) {
        casccovs(straHelper.cascade.covariance);
        histos.fill(HIST("hTableBuildingStatistics"), kCascCovs);
      }

      //_________________________________________________________
      // MC handling part
      if constexpr (soa::is_table<TMCParticles>) {
        // only worry about this if someone else worried about this
        if ((mEnabledTables[kCascMCCores] || mEnabledTables[kMcCascLabels] || mEnabledTables[kCascMCCollRefs])) {
          extractMonteCarloProperties(posTrack, negTrack, bachTrack, mcParticles);

          // Construct label table (note: this will be joinable with CascDatas)
          if (mEnabledTables[kMcCascLabels]) {
            casclabels(
              thisCascInfo.label, thisCascInfo.motherLabel);
            histos.fill(HIST("hTableBuildingStatistics"), kMcCascLabels);
          }

          // Mark mcParticle as recoed (no searching necessary afterwards)
          if (thisCascInfo.label > -1) {
            mcParticleIsReco[thisCascInfo.label] = true;
          }

          if (cascadeBuilderOpts.mc_populateCascMCCoresSymmetric) {
            if (mEnabledTables[kCascMCCores]) {
              cascmccores(
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
            if (mEnabledTables[kCascMCCollRefs]) {
              cascmccollrefs(thisCascInfo.mcCollision);
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
            if (mEnabledTables[kCascCoreMCLabels]) {
              cascCoreMClabels(thisCascMCCoreIndex); // interlink: reconstructed -> MC index
              histos.fill(HIST("hTableBuildingStatistics"), kCascCoreMCLabels);
            }
          }

        } // enabled tables check

        // if BB tags requested, generate them now
        if (mEnabledTables[kMcCascBBTags]) {
          bool bbTag = false;
          if (bachTrack.has_mcParticle()) {
            auto bachelorParticle = bachTrack.template mcParticle_as<aod::McParticles>();
            if (bachelorParticle.pdgCode() == 211) { // pi+, look for antiproton in negative prong
              if (negTrack.has_mcParticle()) {
                auto baryonParticle = negTrack.template mcParticle_as<aod::McParticles>();
                if (baryonParticle.has_mothers() && bachelorParticle.has_mothers() && baryonParticle.pdgCode() == -2212) {
                  for (auto& baryonMother : baryonParticle.template mothers_as<aod::McParticles>()) {
                    for (auto& pionMother : bachelorParticle.template mothers_as<aod::McParticles>()) {
                      if (baryonMother.globalIndex() == pionMother.globalIndex() && baryonMother.pdgCode() == -3122) {
                        bbTag = true;
                      }
                    }
                  }
                }
              }
            } // end if-pion
            if (bachelorParticle.pdgCode() == -211) { // pi-, look for proton in positive prong
              if (posTrack.has_mcParticle()) {
                auto baryonParticle = posTrack.template mcParticle_as<aod::McParticles>();
                if (baryonParticle.has_mothers() && bachelorParticle.has_mothers() && baryonParticle.pdgCode() == 2212) {
                  for (auto& baryonMother : baryonParticle.template mothers_as<aod::McParticles>()) {
                    for (auto& pionMother : bachelorParticle.template mothers_as<aod::McParticles>()) {
                      if (baryonMother.globalIndex() == pionMother.globalIndex() && baryonMother.pdgCode() == 3122) {
                        bbTag = true;
                      }
                    }
                  }
                }
              }
            } // end if-pion
          } // end bachelor has mcparticle
          // Construct label table (note: this will be joinable with CascDatas)
          bbtags(bbTag);
          histos.fill(HIST("hTableBuildingStatistics"), kMcCascBBTags);
        } // end BB tag table enabled check

      } // constexpr requires mcParticles check
    } // cascades loop

    //_________________________________________________________
    // MC handling part
    if constexpr (soa::is_table<TMCParticles>) {
      if ((mEnabledTables[kCascMCCores] || mEnabledTables[kMcCascLabels] || mEnabledTables[kCascMCCollRefs])) {
        // now populate V0MCCores if in asymmetric mode
        if (cascadeBuilderOpts.mc_populateCascMCCoresAsymmetric) {
          // first step: add any un-recoed v0mmcores that were requested
          for (auto& mcParticle : mcParticles) {
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

            if (TMath::Abs(mcParticle.y()) > cascadeBuilderOpts.mc_rapidityWindow)
              continue; // skip outside midrapidity

            if (
              (cascadeBuilderOpts.mc_addGeneratedXiMinus && mcParticle.pdgCode() == 3312) ||
              (cascadeBuilderOpts.mc_addGeneratedXiPlus && mcParticle.pdgCode() == -3312) ||
              (cascadeBuilderOpts.mc_addGeneratedOmegaMinus && mcParticle.pdgCode() == 3334) ||
              (cascadeBuilderOpts.mc_addGeneratedOmegaPlus && mcParticle.pdgCode() == -3334)) {
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
                for (auto& dau : daughters) {
                  if (dau.getProcess() != 4) // check whether the daughter comes from a decay
                    continue;

                  if (TMath::Abs(dau.pdgCode()) == 211 || TMath::Abs(dau.pdgCode()) == 321) {
                    thisCascInfo.pdgCodeBachelor = dau.pdgCode();
                    thisCascInfo.bachP[0] = dau.px();
                    thisCascInfo.bachP[1] = dau.py();
                    thisCascInfo.bachP[2] = dau.pz();
                    thisCascInfo.xyz[0] = dau.vx();
                    thisCascInfo.xyz[1] = dau.vy();
                    thisCascInfo.xyz[2] = dau.vz();
                    thisCascInfo.mcParticleBachelor = dau.globalIndex();
                  }
                  if (TMath::Abs(dau.pdgCode()) == 2212) {
                    thisCascInfo.pdgCodeV0 = dau.pdgCode();

                    for (auto& v0Dau : dau.template daughters_as<aod::McParticles>()) {
                      if (v0Dau.getProcess() != 4)
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

          for (auto thisInfoToFill : mcCascinfos) {
            if (mEnabledTables[kCascMCCores]) {
              cascmccores( // a lot of the info below will be compressed in case of not-recoed MC (good!)
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
            if (mEnabledTables[kCascMCCollRefs]) {
              cascmccollrefs(thisInfoToFill.mcCollision);
              histos.fill(HIST("hTableBuildingStatistics"), kCascMCCollRefs);
            }
          }
        }
      } // enabled tables check
    } // constexpr requires mcParticles check

    LOGF(debug, "Cascades in DF: %i, cascades built: %i", cascades.size(), nCascades);
  }

  //__________________________________________________
  template <class TTracks, typename TCascades, typename TMCParticles>
  void buildKFCascades(TCascades const& cascades, TMCParticles const& mcParticles)
  {
    if (!mEnabledTables[kStoredKFCascCores]) {
      return; // don't do if no request for cascades in place
    }
    int nCascades = 0;
    // Loops over all cascades in the time frame
    histos.fill(HIST("hInputStatistics"), kStoredKFCascCores, cascades.size());
    for (auto& cascade : cascades) {
      // Get tracks and generate candidate
      auto const& collision = cascade.collision();
      auto const& v0 = cascade.v0();
      auto const& posTrack = v0.template posTrack_as<TTracks>();
      auto const& negTrack = v0.template negTrack_as<TTracks>();
      auto const& bachTrack = cascade.template bachelor_as<TTracks>();
      if (!straHelper.buildCascadeCandidateWithKF(collision,
                                                  posTrack,
                                                  negTrack,
                                                  bachTrack,
                                                  mEnabledTables[kCascBBs],
                                                  cascadeBuilderOpts.kfConstructMethod,
                                                  cascadeBuilderOpts.kfTuneForOmega,
                                                  cascadeBuilderOpts.kfUseV0MassConstraint,
                                                  cascadeBuilderOpts.kfUseCascadeMassConstraint,
                                                  cascadeBuilderOpts.kfDoDCAFitterPreMinimV0,
                                                  cascadeBuilderOpts.kfDoDCAFitterPreMinimCasc)) {
        kfcascdataLink(-1);
        interlinks.cascadeToKFCascCores.push_back(-1);
        continue; // didn't work out, skip
      }
      nCascades++;

      // generate analysis tables as required
      if (mEnabledTables[kKFCascIndices]) {
        kfcascidx(cascade.globalIndex(),
                  straHelper.cascade.positiveTrack, straHelper.cascade.negativeTrack,
                  straHelper.cascade.bachelorTrack, straHelper.cascade.collisionId);
        histos.fill(HIST("hTableBuildingStatistics"), kKFCascIndices);
      }
      if (mEnabledTables[kStoredKFCascCores]) {
        kfcascdata(straHelper.cascade.charge, straHelper.cascade.massXi, straHelper.cascade.massOmega,
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
        kfcascdataLink(kfcascdata.lastIndex());
        interlinks.kfCascCoreToCascades.push_back(cascade.globalIndex());
        interlinks.cascadeToKFCascCores.push_back(kfcascdata.lastIndex());
      }
      if (mEnabledTables[kKFCascCovs]) {
        kfcasccovs(straHelper.cascade.covariance, straHelper.cascade.kfTrackCovarianceV0, straHelper.cascade.kfTrackCovariancePos, straHelper.cascade.kfTrackCovarianceNeg);
        histos.fill(HIST("hTableBuildingStatistics"), kKFCascCovs);
      }

      //_________________________________________________________
      // MC handling part (labels only)
      if constexpr (soa::is_table<TMCParticles>) {
        // only worry about this if someone else worried about this
        if ((mEnabledTables[kMcKFCascLabels])) {
          extractMonteCarloProperties(posTrack, negTrack, bachTrack, mcParticles);

          // Construct label table (note: this will be joinable with KFCascDatas)
          kfcasclabels(thisCascInfo.label);
          histos.fill(HIST("hTableBuildingStatistics"), kMcKFCascLabels);
        } // enabled tables check
      } // constexpr requires mcParticles check
    } // end loop over cascades

    LOGF(debug, "KF Cascades in DF: %i, KF cascades built: %i", cascades.size(), nCascades);
  }

  //__________________________________________________
  template <class TTracks, typename TStrangeTracks, typename TMCParticles>
  void buildTrackedCascades(TStrangeTracks const& cascadeTracks, TMCParticles const& mcParticles)
  {
    if (!mEnabledTables[kStoredTraCascCores]) {
      return; // don't do if no request for cascades in place
    }
    int nCascades = 0;
    // Loops over all V0s in the time frame
    histos.fill(HIST("hInputStatistics"), kStoredTraCascCores, cascadeTracks.size());
    for (auto& cascadeTrack : cascadeTracks) {
      // Get tracks and generate candidate
      if (!cascadeTrack.has_track())
        continue; // safety (should be fine but depends on future stratrack dev)

      auto const& strangeTrack = cascadeTrack.template track_as<TTracks>();
      auto const& collision = strangeTrack.collision();
      auto const& cascade = cascadeTrack.cascade();
      auto const& v0 = cascade.v0();
      auto const& posTrack = v0.template posTrack_as<TTracks>();
      auto const& negTrack = v0.template negTrack_as<TTracks>();
      auto const& bachTrack = cascade.template bachelor_as<TTracks>();
      if (!straHelper.buildCascadeCandidate(collision,
                                            posTrack,
                                            negTrack,
                                            bachTrack,
                                            mEnabledTables[kCascBBs],
                                            cascadeBuilderOpts.useCascadeMomentumAtPrimVtx,
                                            mEnabledTables[kCascCovs])) {
        tracascdataLink(-1);
        interlinks.cascadeToTraCascCores.push_back(-1);
        continue; // didn't work out, skip
      }

      // recalculate DCAxy, DCAz with strange track
      auto strangeTrackParCov = getTrackParCov(strangeTrack);
      gpu::gpustd::array<float, 2> dcaInfo;
      strangeTrackParCov.setPID(o2::track::PID::XiMinus); // FIXME: not OK for omegas
      o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, strangeTrackParCov, 2.f, straHelper.fitter.getMatCorrType(), &dcaInfo);
      straHelper.cascade.cascadeDCAxy = dcaInfo[0];
      straHelper.cascade.cascadeDCAz = dcaInfo[1];

      // get momentum from strange track (should not be very different)
      strangeTrackParCov.getPxPyPzGlo(straHelper.cascade.cascadeMomentum);

      // accounting
      nCascades++;

      // generate analysis tables as required
      if (mEnabledTables[kTraCascIndices]) {
        tracascidx(cascade.globalIndex(),
                   straHelper.cascade.positiveTrack, straHelper.cascade.negativeTrack,
                   straHelper.cascade.bachelorTrack, cascadeTrack.trackId(), straHelper.cascade.collisionId);
        histos.fill(HIST("hTableBuildingStatistics"), kTraCascIndices);
      }
      if (mEnabledTables[kStoredTraCascCores]) {
        tracascdata(straHelper.cascade.charge, cascadeTrack.xiMass(), cascadeTrack.omegaMass(),
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
        tracascdataLink(tracascdata.lastIndex());
        interlinks.traCascCoreToCascades.push_back(cascade.globalIndex());
        interlinks.cascadeToTraCascCores.push_back(tracascdata.lastIndex());
      }
      if (mEnabledTables[kCascCovs]) {
        std::array<float, 21> traCovMat = {0.};
        strangeTrackParCov.getCovXYZPxPyPzGlo(traCovMat);
        float traCovMatArray[21];
        for (int ii = 0; ii < 21; ii++) {
          traCovMatArray[ii] = traCovMat[ii];
        }
        tracasccovs(traCovMatArray);
        histos.fill(HIST("hTableBuildingStatistics"), kCascCovs);
      }

      //_________________________________________________________
      // MC handling part (labels only)
      if constexpr (soa::is_table<TMCParticles>) {
        // only worry about this if someone else worried about this
        if ((mEnabledTables[kMcTraCascLabels])) {
          extractMonteCarloProperties(posTrack, negTrack, bachTrack, mcParticles);

          // Construct label table (note: this will be joinable with KFCascDatas)
          tracasclabels(thisCascInfo.label);
          histos.fill(HIST("hTableBuildingStatistics"), kMcTraCascLabels);
        } // enabled tables check
      } // constexpr requires mcParticles check
    } // end loop over cascades
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
          if (std::abs(part.pdgCode()) == 13 && treatPiToMuDecays) {
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
  template <typename TCollisions, typename TV0s, typename TCascades, typename TTrackedCascades, typename TTracks, typename TBCs, typename TMCParticles>
  void dataProcess(TCollisions const& collisions, TV0s const& v0s, TCascades const& cascades, TTrackedCascades const& trackedCascades, TTracks const&, TBCs const& bcs, TMCParticles const& mcParticles)
  {
    if (!initCCDB(bcs, collisions))
      return;

    // reset vectors for cascade interlinks
    resetInterlinks();

    // mark V0s that will be buffered for the cascade building
    markV0sUsedInCascades(v0s, cascades, trackedCascades);

    // build V0s
    buildV0s<TTracks>(v0s, mcParticles);

    // build cascades
    buildCascades<TTracks>(cascades, mcParticles);
    buildKFCascades<TTracks>(cascades, mcParticles);

    // build tracked cascades only if subscription is Run 3 like (doesn't exist in Run 2)
    if constexpr (soa::is_table<TTrackedCascades>) {
      buildTrackedCascades<TTracks>(trackedCascades, mcParticles);
    }

    populateCascadeInterlinks();
  }

  void processRealData(aod::Collisions const& collisions, aod::V0s const& v0s, aod::Cascades const& cascades, aod::TrackedCascades const& trackedCascades, FullTracksExtIU const& tracks, aod::BCsWithTimestamps const& bcs)
  {
    dataProcess(collisions, v0s, cascades, trackedCascades, tracks, bcs, (TObject*)nullptr);
  }

  void processRealDataRun2(aod::Collisions const& collisions, aod::V0s const& v0s, aod::Cascades const& cascades, FullTracksExt const& tracks, aod::BCsWithTimestamps const& bcs)
  {
    dataProcess(collisions, v0s, cascades, (TObject*)nullptr, tracks, bcs, (TObject*)nullptr);
  }

  void processMonteCarlo(aod::Collisions const& collisions, aod::V0s const& v0s, aod::Cascades const& cascades, aod::TrackedCascades const& trackedCascades, FullTracksExtLabeledIU const& tracks, aod::BCsWithTimestamps const& bcs, aod::McParticles const& mcParticles)
  {
    dataProcess(collisions, v0s, cascades, trackedCascades, tracks, bcs, mcParticles);
  }

  void processMonteCarloRun2(aod::Collisions const& collisions, aod::V0s const& v0s, aod::Cascades const& cascades, FullTracksExtLabeled const& tracks, aod::BCsWithTimestamps const& bcs, aod::McParticles const& mcParticles)
  {
    dataProcess(collisions, v0s, cascades, (TObject*)nullptr, tracks, bcs, mcParticles);
  }

  void processSimulationFindable(aod::Collisions const& collisions, aod::FindableV0s const& v0s, aod::Cascades const& cascades, FullTracksExtIU const& tracks, aod::BCsWithTimestamps const& bcs)
  {
    dataProcess(collisions, v0s, cascades, (TObject*)nullptr, tracks, bcs, (TObject*)nullptr);
  }

  PROCESS_SWITCH(StrangenessBuilder, processRealData, "process real data", true);
  PROCESS_SWITCH(StrangenessBuilder, processRealDataRun2, "process real data (Run 2)", false);
  PROCESS_SWITCH(StrangenessBuilder, processMonteCarlo, "process monte carlo", false);
  PROCESS_SWITCH(StrangenessBuilder, processMonteCarloRun2, "process monte carlo (Run 2)", false);
  PROCESS_SWITCH(StrangenessBuilder, processSimulationFindable, "process simulation findable (requires lambdakzeromcfinder)", false);
};

// Extends the v0data table with expression columns
struct strangenessbuilderInitializer {
  Spawns<aod::V0Cores> v0cores;
  Spawns<aod::CascCores> cascdataext;
  Spawns<aod::KFCascCores> kfcascdataext;
  Spawns<aod::TraCascCores> tracascdataext;
  void init(InitContext const&) {}
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<StrangenessBuilder>(cfgc),
    adaptAnalysisTask<strangenessbuilderInitializer>(cfgc)};
}
