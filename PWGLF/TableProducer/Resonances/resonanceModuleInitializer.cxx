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
/// \file resonanceModuleInitializer.cxx
/// \brief Initializes variables for the resonance candidate producers
///
/// \author Bong-Hwi Lim <bong-hwi.lim@cern.ch>, Minjae Kim <minjae.kim@cern.ch>
/// \since  Aug.31 2026

#include "PWGLF/DataModel/LFResonanceTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/mcCentrality.h"
#include "PWGLF/Utils/collisionCuts.h"

#include "Common/CCDB/RCTSelectionFlags.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <DataFormatsParameters/GRPObject.h>
#include <DetectorsBase/Propagator.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>
#include <THnSparse.h>
#include <TPDGCode.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::constants::physics;
using namespace o2::constants::math;
using namespace o2::aod::rctsel;

/**
 * @brief Initializer for the event pool for resonance study
 *
 * This struct is responsible for initializing and processing collision data
 * for resonance studies. It handles event selection, centrality estimation,
 * and QA histogram filling.
 */
// Framework Service members are populated by the analysis framework.
// NOLINTNEXTLINE(cppcoreguidelines-pro-type-member-init)
struct ResonanceModuleInitializer {
  static constexpr double BzOverrideThreshold = -990.;
  static constexpr double MinimumNonzeroBz = 1.e-5;
  static constexpr int CentralityFT0M = 0;
  static constexpr int CentralityFT0C = 1;
  static constexpr int CentralityFT0A = 2;
  static constexpr int CentralityFV0A = 3;
  static constexpr int MultiplicityNTracksPV = 0;
  static constexpr int MultiplicityNTracksPVeta1 = 1;
  static constexpr int MultiplicityNTracksPVetaHalf = 2;
  static constexpr int MultiplicityFT0M = 3;
  static constexpr int MultiplicityFT0A = 4;
  static constexpr int MultiplicityFT0C = 5;
  static constexpr int MultiplicityFV0A = 6;
  static constexpr int DetailedQARCTStage = o2::analysis::CollisonCuts::kAllpassed + 1;
  static constexpr int DetailedQAStages = DetailedQARCTStage + 1;
  static constexpr float MCVertexZMax = 10.f;
  // PDG codes used by the persistent resonance-parent selection. Named O2/ROOT
  // values are preferred where available; the remaining resonances are kept as
  // local named constants because neither PDG_t nor PhysicsConstants defines them.
  static constexpr int PdgKStar0 = o2::constants::physics::Pdg::kK0Star892;
  static constexpr int PdgKStarCharged = o2::constants::physics::Pdg::kKPlusStar892;
  static constexpr int PdgPhi = o2::constants::physics::Pdg::kPhi;
  static constexpr int F0Code980 = 9010221;
  static constexpr int F0Code1370 = 10221;
  static constexpr int F0Code1500 = 9030221;
  static constexpr int F0Code1710 = 10331;
  static constexpr int F1Code1285 = 20223;
  static constexpr int F1Code1420 = 20333;
  static constexpr int F2PrimeCode1525 = 335;
  static constexpr int PdgRho0 = PDG_t::kRho770_0;
  static constexpr int PdgRhoCharged = PDG_t::kRho770Plus;
  static constexpr int SigmaStarPlusCode = 3224;
  static constexpr int PdgLambda1520 = o2::constants::physics::Pdg::kLambda1520_Py;
  static constexpr int Xi1530Code = 3324;
  static constexpr int PdgK1Plus1270 = o2::constants::physics::Pdg::kK1_1270Plus;
  static constexpr int Xi1820NeutralCode = 123314;
  static constexpr int Xi1820MinusCode = 123324;
  static constexpr int Omega2012MinusCode = 123334;
  static constexpr int PdgProton = PDG_t::kProton;
  static constexpr int PdgLambda0 = PDG_t::kLambda0;
  static constexpr int PdgXiMinus = PDG_t::kXiMinus;
  static constexpr int PdgXi0 = o2::constants::physics::Pdg::kXi0;
  static constexpr int PdgOmegaMinus = PDG_t::kOmegaMinus;

  int mRunNumber = 0;                       ///< Run number for the current data
  int multEstimator = CentralityFT0M;       ///< Centrality estimator type
  float dBz = 0.f;                          ///< Magnetic field value
  float centrality = 0.f;                   ///< Centrality value for the event
  Service<o2::ccdb::BasicCCDBManager> ccdb; ///< CCDB manager service

  Produces<aod::ResoCollisions_001> resoCollisions;              ///< Output table for resonance collisions
  Produces<aod::ResoCollisionColls> resoCollisionColls;          ///< Optional source collision soft links
  Produces<aod::ResoCollisionGroups_001> resoCollisionGroups001; ///< Scalar original-collision grouping keys
  Produces<aod::ResoMCCollisions_001> resoMCCollisions001;       ///< Generator-only MC collision extension
  Produces<aod::ResoMCCollisionIds> resoMCCollisionIds;          ///< Optional source generator-collision soft links
  Produces<aod::ResoMCParents_001> reso2mcparents;               ///< Generated parents with scalar source-particle IDs

  // CCDB options
  struct : ConfigurableGroup {
    Configurable<std::string> ccdbURL{"ccdbURL", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
    Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
    Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
    Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
    Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
    Configurable<bool> cfgFatalWhenNull{"cfgFatalWhenNull", true, "Fatal when null on ccdb access"};
    Configurable<bool> cfgBypassCollIndexFill{"cfgBypassCollIndexFill", false, "Deprecated compatibility option; collision mapping tables are always written"};
  } CCDB;

  // General event options
  struct : ConfigurableGroup {
    Configurable<double> dBzInput{"dBzInput", -999, "bz field, -999 is automatic"};
    Configurable<bool> cfgFillQA{"cfgFillQA", true, "Fill QA histograms"};
    Configurable<bool> cfgFillDetailedQA{"cfgFillDetailedQA", true, "Fill the Run 3 event-selection stage vs vertex-z vs centrality vs multiplicity THnSparse"};
    Configurable<bool> cfgBypassCCDB{"cfgBypassCCDB", true, "Bypass loading CCDB part to save CPU time and memory"}; // will be affected to b_z value.
    Configurable<std::string> cfgMultName{"cfgMultName", "FT0M", "Centrality estimator: FT0M, FT0C, FT0A, or FV0A"};
    Configurable<int> cfgMultiplicityEstimator{"cfgMultiplicityEstimator", 3,
                                               "Stored multiplicity (NOT percentile): 0 -> NTracksPV, 1 -> NTracksPVeta1, 2 -> NTracksPVetaHalf, 3 -> FT0M, 4 -> FT0A, 5 -> FT0C, 6 -> FV0A"};
    ConfigurableAxis binsCent{"binsCent", {VARIABLE_WIDTH, 0., 0.01, 0.1, 1., 5., 10., 15., 20., 30., 40., 50., 60., 70., 80., 90., 100., 105.}, "Binning of the centrality axis"};
    ConfigurableAxis binsMultiplicity{"binsMultiplicity", {500, 0.f, 5000.f}, "Binning of the reconstructed multiplicity axis for detailed collision QA"};
    ConfigurableAxis cfgVtxBins{"cfgVtxBins", {400, -20.f, 20.f}, "Binning of the collision vertex-z axis for detailed QA"};
  } EventConfig;

  /// Event cuts
  o2::analysis::CollisonCuts colCuts;
  struct : ConfigurableGroup {
    Configurable<float> cfgEvtZvtx{"cfgEvtZvtx", 10.f, "Evt sel: Max. z-Vertex (cm)"};
    Configurable<int> cfgEvtOccupancyInTimeRange{"cfgEvtOccupancyInTimeRange", -1, "Evt sel: maximum track occupancy"};
    Configurable<bool> cfgEvtTriggerCheck{"cfgEvtTriggerCheck", false, "Evt sel: check for trigger"};
    Configurable<bool> cfgEvtOfflineCheck{"cfgEvtOfflineCheck", false, "Evt sel: check for offline selection (sel8)"};
    Configurable<bool> cfgEvtTriggerTVXSel{"cfgEvtTriggerTVXSel", true, "Evt sel: triggerTVX selection (MB)"};
    Configurable<bool> cfgEvtTFBorderCut{"cfgEvtTFBorderCut", true, "Evt sel: apply TF border cut"};
    Configurable<bool> cfgEvtUseITSTPCvertex{"cfgEvtUseITSTPCvertex", false, "Evt sel: use at lease on ITS-TPC track for vertexing"};
    Configurable<bool> cfgEvtCollInTimeRangeNarrow{"cfgEvtCollInTimeRangeNarrow", false, "Evt sel: apply NoCollInTimeRangeNarrow"};
    Configurable<bool> cfgEvtZvertexTimedifference{"cfgEvtZvertexTimedifference", false, "Evt sel: apply Z-vertex time difference"};
    Configurable<bool> cfgEvtPileupRejection{"cfgEvtPileupRejection", false, "Evt sel: apply pileup rejection"};
    Configurable<bool> cfgEvtNoITSROBorderCut{"cfgEvtNoITSROBorderCut", false, "Evt sel: apply NoITSRO border cut"};
    Configurable<bool> cfgEvtRun2AliEventCuts{"cfgEvtRun2AliEventCuts", false, "Evt sel: apply Run2 AliEventCuts"};
    Configurable<bool> cfgEvtRun2INELgtZERO{"cfgEvtRun2INELgtZERO", false, "Evt sel: apply Run2 INELgtZERO"};
    Configurable<bool> cfgEvtUseRCTFlagChecker{"cfgEvtUseRCTFlagChecker", true, "Evt sel: use RCT flag checker"};
    Configurable<bool> cfgEvtBCRCT{"cfgEvtBCRCT", false, "Evt sel: check RCT on the nominal associated BCSEL instead of the collision EVSEL"};
    Configurable<std::string> cfgEvtRCTFlagCheckerLabel{"cfgEvtRCTFlagCheckerLabel", "CBT_hadronPID", "Evt sel: RCT flag checker label"};
    Configurable<bool> cfgEvtRCTFlagCheckerZDCCheck{"cfgEvtRCTFlagCheckerZDCCheck", false, "Evt sel: RCT flag checker ZDC check"};
    Configurable<bool> cfgEvtRCTFlagCheckerLimitAcceptAsBad{"cfgEvtRCTFlagCheckerLimitAcceptAsBad", false, "Evt sel: RCT flag checker treat Limited Acceptance As Bad"};
    Configurable<bool> cfgEvtRCTCheckTableValidity{"cfgEvtRCTCheckTableValidity", false, "Evt sel: reject collisions when the RCT CCDB payload is unavailable"};
  } EventCuts;
  RCTFlagsChecker recoRCTChecker;

  // Generator-level event and resonance QA
  struct : ConfigurableGroup {
    Configurable<bool> cfgGenBCRCT{"cfgGenBCRCT", false, "GenEvent: apply the RCT flag checker to the associated BC"};
    Configurable<bool> cfgGenRCTCheckTableValidity{"cfgGenRCTCheckTableValidity", false, "GenEvent: reject MC collisions when the RCT CCDB payload is unavailable"};
    Configurable<bool> cfgGenMult05{"cfgGenMult05", true, "GenEvent: multiplicity in |eta| < 0.5"};
    Configurable<bool> cfgGenMult10{"cfgGenMult10", false, "GenEvent: multiplicity in |eta| < 1.0"};
    Configurable<bool> cfgGenMultFT0M{"cfgGenMultFT0M", false, "GenEvent: generated charged-particle multiplicity in the FT0A + FT0C acceptance"};
    Configurable<bool> cfgGenMultFT0C{"cfgGenMultFT0C", false, "GenEvent: generated charged-particle multiplicity in the FT0C acceptance"};
    Configurable<bool> cfgGenMultFV0A{"cfgGenMultFV0A", false, "GenEvent: generated charged-particle multiplicity in the FV0A acceptance"};
    Configurable<bool> cfgGenMultPercentile{"cfgGenMultPercentile", true, "Use the configured FT0M, FT0C, or FV0A percentile from the MC centrality wagon"};
    Configurable<bool> cfgFillMCCollisionSoftLink{"cfgFillMCCollisionSoftLink", false,
                                                  "Write the source MC-collision soft link; enable only with the original AO2D as linked parent"};
    Configurable<bool> isZvtxcutGen{"isZvtxcutGen", true, "Apply the generator-collision z-vertex cut"};
    Configurable<float> cutzvertexGen{"cutzvertexGen", 10.f, "Maximum absolute generator-collision z vertex (cm)"};
    Configurable<bool> checkIsTrueINELgt0{"checkIsTrueINELgt0", true, "Classify true INEL>0 generator collisions"};
    ConfigurableAxis binsCentGen{"binsCentGen",
                                 {VARIABLE_WIDTH, 0., 0.01, 0.1, 1., 5., 10., 15., 20., 30., 40., 50., 60., 70., 80., 90., 100., 105.},
                                 "Generator centrality axis"};
    ConfigurableAxis ptAxisGen{"ptAxisGen", {400, 0.f, 20.f}, "#it{p}_{T} (GeV/#it{c})"};
    ConfigurableAxis multNTracksAxis{"multNTracksAxis", {500, 0.f, 5000.f}, "Number of charged particles"};
    ConfigurableAxis impactParameterAxis{"impactParameterAxis", {500, 0.f, 50.f}, "Impact parameter (fm)"};
    Configurable<bool> isDaughterCheck{"isDaughterCheck", true, "Require the configured two-body decay"};
    Configurable<float> cfgRapidityCutMinGen{"cfgRapidityCutMinGen", -0.5f, "Minimum generated-particle rapidity"};
    Configurable<float> cfgRapidityCutMaxGen{"cfgRapidityCutMaxGen", 0.5f, "Maximum generated-particle rapidity"};
    Configurable<int> pdgTruthMother{"pdgTruthMother", static_cast<int>(Xi1530Code), "Absolute PDG code of the generated mother"};
    Configurable<int> pdgTruthDaughter1{"pdgTruthDaughter1", static_cast<int>(PdgXiMinus), "Absolute PDG code of the first daughter"};
    Configurable<int> pdgTruthDaughter2{"pdgTruthDaughter2", PDG_t::kPiPlus, "Absolute PDG code of the second daughter"};
    Configurable<bool> cfgDoSignalLoss{"cfgDoSignalLoss", false, "Save reference particles for mT-scaling signal-loss studies"};
  } GenCuts;
  RCTFlagsChecker genRCTChecker;

  // Keep the established ResoMCParents content compatible with the legacy
  // initializer. The additional stable-particle species are written only for
  // signal-loss studies and are filtered in fillMCParents.
  Partition<aod::McParticles> selectedMCParticles = (nabs(aod::mcparticle::pdgCode) == PdgKStar0)             // K*(892)0
                                                    || (nabs(aod::mcparticle::pdgCode) == PdgKStarCharged)    // K*(892)+
                                                    || (nabs(aod::mcparticle::pdgCode) == PdgPhi)             // phi(1020)
                                                    || (nabs(aod::mcparticle::pdgCode) == F0Code980)          // f0(980)
                                                    || (nabs(aod::mcparticle::pdgCode) == F0Code1370)         // f0(1370)
                                                    || (nabs(aod::mcparticle::pdgCode) == F0Code1500)         // f0(1500)
                                                    || (nabs(aod::mcparticle::pdgCode) == F0Code1710)         // f0(1710)
                                                    || (nabs(aod::mcparticle::pdgCode) == F1Code1285)         // f1(1285)
                                                    || (nabs(aod::mcparticle::pdgCode) == F1Code1420)         // f1(1420)
                                                    || (nabs(aod::mcparticle::pdgCode) == F2PrimeCode1525)    // f2'(1525)
                                                    || (nabs(aod::mcparticle::pdgCode) == PdgRho0)            // rho(770)0
                                                    || (nabs(aod::mcparticle::pdgCode) == PdgRhoCharged)      // rho(770)+
                                                    || (nabs(aod::mcparticle::pdgCode) == SigmaStarPlusCode)  // Sigma(1385)+
                                                    || (nabs(aod::mcparticle::pdgCode) == PdgLambda1520)      // Lambda(1520)
                                                    || (nabs(aod::mcparticle::pdgCode) == Xi1530Code)         // Xi(1530)0
                                                    || (nabs(aod::mcparticle::pdgCode) == PdgK1Plus1270)      // K1(1270)+
                                                    || (nabs(aod::mcparticle::pdgCode) == Xi1820NeutralCode)  // Xi(1820)0
                                                    || (nabs(aod::mcparticle::pdgCode) == Xi1820MinusCode)    // Xi(1820)-
                                                    || (nabs(aod::mcparticle::pdgCode) == Omega2012MinusCode) // Omega(2012)-
                                                    || (nabs(aod::mcparticle::pdgCode) == PdgProton)          // proton
                                                    || (nabs(aod::mcparticle::pdgCode) == PdgLambda0)         // Lambda0
                                                    || (nabs(aod::mcparticle::pdgCode) == PdgXiMinus)         // Xi-
                                                    || (nabs(aod::mcparticle::pdgCode) == PdgXi0)             // Xi0
                                                    || (nabs(aod::mcparticle::pdgCode) == PdgOmegaMinus);     // Omega-
  Preslice<aod::McParticles> mcParticlesPerMcCollision = aod::mcparticle::mcCollisionId;

  HistogramRegistry qaRegistry{"QAHistos", {}, OutputObjHandlingPolicy::AnalysisObject};

  /**
   * @brief Initializes the task
   *
   * @param context Initialization context
   */
  void init(InitContext&)
  {
    mRunNumber = 0;
    dBz = 0;
    centrality = 0;
    // Determine the centrality estimator based on the configuration.
    if (EventConfig.cfgMultName.value == "FT0M") {
      multEstimator = CentralityFT0M;
    } else if (EventConfig.cfgMultName.value == "FT0C") {
      multEstimator = CentralityFT0C;
    } else if (EventConfig.cfgMultName.value == "FT0A") {
      multEstimator = CentralityFT0A;
    } else if (EventConfig.cfgMultName.value == "FV0A") {
      multEstimator = CentralityFV0A;
    } else {
      LOGF(fatal, "Unsupported cfgMultName '%s'; choose FT0M, FT0C, FT0A, or FV0A", EventConfig.cfgMultName.value.c_str());
    }
    LOGF(info, "Centrality estimator: %d, %s", multEstimator, EventConfig.cfgMultName.value.c_str());
    if (EventConfig.cfgMultiplicityEstimator.value < MultiplicityNTracksPV ||
        EventConfig.cfgMultiplicityEstimator.value > MultiplicityFV0A) {
      LOG(fatal) << "cfgMultiplicityEstimator must be in the range [0, 6]";
    }
    LOGF(info, "Stored collision multiplicity estimator: %d", EventConfig.cfgMultiplicityEstimator.value);

    // Run 2 and Run 3 callbacks require different event-selection semantics.
    const bool anyRun2Process = doprocessRun2 || doprocessRun2MC;
    const bool anyRun3Process = doprocessRun3 || doprocessRun3MC || doprocessMCgen;
    if (anyRun2Process && anyRun3Process) {
      LOG(fatal) << "Run 2 and Run 3 processes cannot be enabled in the same ResonanceModuleInitializer";
    }
    if (doprocessRun2 && doprocessRun2MC) {
      LOG(fatal) << "processRun2MC writes both ResoCollisions and ResoMCCollisions_001; do not enable processRun2 with it";
    }
    if (doprocessRun3 && doprocessRun3MC) {
      LOG(fatal) << "processRun3MC writes both ResoCollisions and ResoMCCollisions_001; do not enable processRun3 with it";
    }
    const int enabledGenMultiplicityEstimators = static_cast<int>(GenCuts.cfgGenMult05.value) +
                                                 static_cast<int>(GenCuts.cfgGenMult10.value) +
                                                 static_cast<int>(GenCuts.cfgGenMultFT0M.value) +
                                                 static_cast<int>(GenCuts.cfgGenMultFT0C.value) +
                                                 static_cast<int>(GenCuts.cfgGenMultFV0A.value);
    if ((doprocessMCgen || doprocessRun2MC || doprocessRun3MC) && enabledGenMultiplicityEstimators > 1) {
      LOG(fatal) << "Only one generator multiplicity estimator can be enabled: cfgGenMult05, cfgGenMult10, cfgGenMultFT0M, cfgGenMultFT0C, or cfgGenMultFV0A";
    }
    if (doprocessMCgen) {
      if (GenCuts.cfgGenMultPercentile && multEstimator != CentralityFT0M &&
          multEstimator != CentralityFT0C && multEstimator != CentralityFV0A) {
        LOGF(fatal, "cfgGenMultPercentile supports cfgMultName=FT0M, FT0C, or FV0A");
      }
      if (GenCuts.isZvtxcutGen &&
          (!std::isfinite(GenCuts.cutzvertexGen.value) || GenCuts.cutzvertexGen.value <= 0.f)) {
        LOG(fatal) << "cutzvertexGen must be finite and positive when the generator vertex cut is enabled";
      }
      if (!std::isfinite(GenCuts.cfgRapidityCutMinGen.value) ||
          !std::isfinite(GenCuts.cfgRapidityCutMaxGen.value) ||
          GenCuts.cfgRapidityCutMinGen.value >= GenCuts.cfgRapidityCutMaxGen.value) {
        LOG(fatal) << "Generator rapidity limits must be finite and satisfy cfgRapidityCutMinGen < cfgRapidityCutMaxGen";
      }
    }

    // Initialize event selection cuts based on the process type
    if (anyRun2Process) {
      colCuts.setCuts(EventCuts.cfgEvtZvtx, EventCuts.cfgEvtTriggerCheck, EventCuts.cfgEvtOfflineCheck, false);
    } else if (anyRun3Process) {
      colCuts.setCuts(EventCuts.cfgEvtZvtx, EventCuts.cfgEvtTriggerCheck, EventCuts.cfgEvtOfflineCheck, true, false, EventCuts.cfgEvtOccupancyInTimeRange);
    }
    colCuts.init(&qaRegistry);
    colCuts.setTriggerTVX(EventCuts.cfgEvtTriggerTVXSel);
    colCuts.setApplyTFBorderCut(EventCuts.cfgEvtTFBorderCut);
    colCuts.setApplyITSTPCvertex(EventCuts.cfgEvtUseITSTPCvertex);
    colCuts.setApplyCollInTimeRangeNarrow(EventCuts.cfgEvtCollInTimeRangeNarrow);
    colCuts.setApplyZvertexTimedifference(EventCuts.cfgEvtZvertexTimedifference);
    colCuts.setApplyPileupRejection(EventCuts.cfgEvtPileupRejection);
    colCuts.setApplyNoITSROBorderCut(EventCuts.cfgEvtNoITSROBorderCut);
    colCuts.setApplyRun2AliEventCuts(EventCuts.cfgEvtRun2AliEventCuts);
    colCuts.setApplyRun2INELgtZERO(EventCuts.cfgEvtRun2INELgtZERO);

    if (EventConfig.cfgFillDetailedQA && (doprocessRun3 || doprocessRun3MC)) {
      AxisSpec selectionStageAxis{DetailedQAStages, -0.5f, static_cast<float>(DetailedQAStages) - 0.5f, "Passed event-selection stage"};
      AxisSpec vertexAxis{EventConfig.cfgVtxBins, "Collision vertex z (cm)"};
      AxisSpec centralityAxis{EventConfig.binsCent, "Centrality (%)"};
      AxisSpec multiplicityAxis{EventConfig.binsMultiplicity, "Multiplicity"};
      qaRegistry.add("Event/h4EventSelectionDetail", "Event-selection cut flow", kTHnSparseD,
                     {selectionStageAxis, vertexAxis, centralityAxis, multiplicityAxis});

      auto detailedQA = qaRegistry.get<THnSparse>(HIST("Event/h4EventSelectionDetail"));
      auto cutCounts = qaRegistry.get<TH1>(HIST("CollCutCounts"));
      for (int stage = o2::analysis::CollisonCuts::kAllEvent;
           stage <= o2::analysis::CollisonCuts::kAllpassed; ++stage) {
        detailedQA->GetAxis(0)->SetBinLabel(stage + 1, cutCounts->GetXaxis()->GetBinLabel(colCuts.binLabel(stage)));
      }
      detailedQA->GetAxis(0)->SetBinLabel(DetailedQARCTStage + 1, "RCT");
    }

    recoRCTChecker.init(EventCuts.cfgEvtRCTFlagCheckerLabel,
                        EventCuts.cfgEvtRCTFlagCheckerZDCCheck,
                        EventCuts.cfgEvtRCTFlagCheckerLimitAcceptAsBad,
                        EventCuts.cfgEvtRCTCheckTableValidity.value);
    genRCTChecker.init(EventCuts.cfgEvtRCTFlagCheckerLabel,
                       EventCuts.cfgEvtRCTFlagCheckerZDCCheck,
                       EventCuts.cfgEvtRCTFlagCheckerLimitAcceptAsBad,
                       GenCuts.cfgGenRCTCheckTableValidity.value);

    // Configure CCDB access if not bypassed
    if (!EventConfig.cfgBypassCCDB) {
      ccdb->setURL(CCDB.ccdbURL.value);
      ccdb->setCaching(true);
      ccdb->setLocalObjectValidityChecking();
      ccdb->setFatalWhenNull(CCDB.cfgFatalWhenNull);
      uint64_t now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
      ccdb->setCreatedNotAfter(now); // TODO must become global parameter from the train creation time
    }

    if (doprocessMCgen) {
      constexpr std::array<char const*, 5> MCEventLabels{"All", "z vertex", "BC RCT", "INEL", "INEL>0"};
      AxisSpec centAxisGen = {GenCuts.binsCentGen, "Centrality (%)"};
      AxisSpec eventTypeAxis = {2, 0.f, 2.f, "Event type"};
      qaRegistry.add("EventGen/hNEventsMC", "Generator event selection", kTH1D, {{5, 0.f, 5.f}});
      auto eventCounter = qaRegistry.get<TH1>(HIST("EventGen/hNEventsMC"));
      for (std::size_t i = 0; i < MCEventLabels.size(); ++i) {
        eventCounter->GetXaxis()->SetBinLabel(static_cast<int>(i + 1), MCEventLabels[i]);
      }
      qaRegistry.add("EventGen/h5ResonanceTruth", "Generated resonance", kTHnSparseD,
                     {eventTypeAxis, GenCuts.ptAxisGen, centAxisGen, GenCuts.multNTracksAxis, GenCuts.impactParameterAxis});
      qaRegistry.add("EventGen/h5ResonanceTruthAnti", "Generated anti-resonance", kTHnSparseD,
                     {eventTypeAxis, GenCuts.ptAxisGen, centAxisGen, GenCuts.multNTracksAxis, GenCuts.impactParameterAxis});
      qaRegistry.add("EventGen/hZCollisionGen", "Generator collision z vertex", kTH1D, {{100, -20.f, 20.f}});
      qaRegistry.add("EventGen/h4MultCent_genMC", "Generator-event multiplicity and centrality", kTHnSparseD,
                     {eventTypeAxis, centAxisGen, GenCuts.multNTracksAxis, GenCuts.impactParameterAxis});
      qaRegistry.add("EventGen/h4MultCent_recMC", "Reconstructed-event multiplicity and centrality", kTHnSparseD,
                     {eventTypeAxis, centAxisGen, GenCuts.multNTracksAxis, GenCuts.impactParameterAxis});
      qaRegistry.add("EventGen/h2CentralityVsMultMC", "Representative reconstructed centrality vs generator multiplicity", kTH2D,
                     {centAxisGen, GenCuts.multNTracksAxis});
    }
  }

  /**
   * @brief Initializes CCDB for a given BC
   *
   * @param bc BC iterator
   */
  template <typename BCType>
  void initCCDB(BCType const& bc) // Simple copy from LambdaKzeroFinder.cxx
  {
    if (EventConfig.cfgBypassCCDB) {
      return;
    }
    if (mRunNumber == bc.runNumber()) {
      return;
    }

    // In case override, don't proceed, please - no CCDB access required
    if (EventConfig.dBzInput > BzOverrideThreshold) {
      dBz = EventConfig.dBzInput;
      o2::parameters::GRPMagField grpmag;
      if (std::fabs(dBz) > MinimumNonzeroBz) {
        grpmag.setL3Current(30000.f / (dBz / 5.0f));
      }
      o2::base::Propagator::initFieldFromGRP(&grpmag);
      mRunNumber = bc.runNumber();
      return;
    }

    auto run3grpTimestamp = bc.timestamp();
    auto* grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>(CCDB.grpPath, run3grpTimestamp);
    o2::parameters::GRPMagField* grpmag = nullptr;
    if (grpo) {
      o2::base::Propagator::initFieldFromGRP(grpo);
      // Fetch magnetic field from ccdb for current collision
      dBz = grpo->getNominalL3Field();
      LOG(info) << "Retrieved GRP for timestamp " << run3grpTimestamp << " with magnetic field of " << dBz << " kZG";
    } else {
      grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(CCDB.grpmagPath, run3grpTimestamp);
      if (!grpmag) {
        LOG(fatal) << "Got nullptr from CCDB for path " << CCDB.grpmagPath << " of object GRPMagField and " << CCDB.grpPath << " of object GRPObject for timestamp " << run3grpTimestamp;
      }
      o2::base::Propagator::initFieldFromGRP(grpmag);
      // Fetch magnetic field from ccdb for current collision
      dBz = std::lround(5.f * grpmag->getL3Current() / 30000.f);
      LOG(info) << "Retrieved GRP for timestamp " << run3grpTimestamp << " with magnetic field of " << dBz << " kZG";
    }
    mRunNumber = bc.runNumber();
    // Set magnetic field value once known
    LOGF(info, "Bz set to %f for run: ", dBz, mRunNumber);
  }

  /**
   * @brief Centrality estimator selection
   *
   * @tparam ResoColl Type of resonance collision
   * @param resoEvents Resonance events
   * @return Centrality value
   */
  template <typename ResoColl>
  float centEst(ResoColl const& resoEvents)
  {
    switch (multEstimator) {
      case CentralityFT0M:
        return resoEvents.centFT0M();
      case CentralityFT0C:
        return resoEvents.centFT0C();
      case CentralityFT0A:
        return resoEvents.centFT0A();
      case CentralityFV0A:
        return resoEvents.centFV0A();
      default:
        return -999.f;
    }
  }

  /**
   * @brief Returns the configured reconstructed multiplicity estimator
   */
  template <typename CollisionType>
  float collisionMultiplicity(CollisionType const& collision)
  {
    switch (EventConfig.cfgMultiplicityEstimator.value) {
      case MultiplicityNTracksPV:
        return collision.multNTracksPV();
      case MultiplicityNTracksPVeta1:
        return collision.multNTracksPVeta1();
      case MultiplicityNTracksPVetaHalf:
        return collision.multNTracksPVetaHalf();
      case MultiplicityFT0M:
        return collision.multFT0M();
      case MultiplicityFT0A:
        return collision.multFT0A();
      case MultiplicityFT0C:
        return collision.multFT0C();
      case MultiplicityFV0A:
        return collision.multFV0A();
      default:
        return -1.f;
    }
  }

  /// Fill the detailed Run 3 collision QA at one passed selection stage.
  template <typename CollisionType>
  void fillDetailedCollisionQA(CollisionType const& collision, int selectionStage)
  {
    qaRegistry.fill(HIST("Event/h4EventSelectionDetail"),
                    selectionStage,
                    collision.posZ(),
                    centEst(collision),
                    collisionMultiplicity(collision));
  }

  /// Reproduce the configured Run 3 collision cut flow for detailed QA only.
  /// The authoritative event decision remains CollisonCuts::isSelected().
  template <typename CollisionType>
  void fillDetailedRun3SelectionQA(CollisionType const& collision)
  {
    fillDetailedCollisionQA(collision, o2::analysis::CollisonCuts::kAllEvent);
    if (std::abs(collision.posZ()) > EventCuts.cfgEvtZvtx.value) {
      return;
    }
    fillDetailedCollisionQA(collision, o2::analysis::CollisonCuts::kFlagZvertex);

// Keep this QA-only mapping synchronized with CollisonCuts without exposing
// its internal selection registry.
// NOLINTNEXTLINE(cppcoreguidelines-macro-usage) -- X-macro interface of EventSelectionFlagsMapping.def
#define EVSEL_FLAG(enumVal, member, defaultVal, evtSelEnum, setter, getter, label, desc) \
  if (colCuts.getSelection(o2::analysis::CollisonCuts::evtSelEnum)) {                    \
    if (!collision.selection_bit(o2::aod::evsel::enumVal)) {                             \
      return;                                                                            \
    }                                                                                    \
    fillDetailedCollisionQA(collision, o2::analysis::CollisonCuts::evtSelEnum);          \
  }
#include "PWGLF/Utils/EventSelectionFlagsMapping.def" // NOLINT(build/include)
#undef EVSEL_FLAG

    if (EventCuts.cfgEvtOfflineCheck.value && !collision.sel8()) {
      return;
    }
    fillDetailedCollisionQA(collision, o2::analysis::CollisonCuts::kFlagSel8);

    if (EventCuts.cfgEvtOccupancyInTimeRange.value > 0 &&
        collision.trackOccupancyInTimeRange() > EventCuts.cfgEvtOccupancyInTimeRange.value) {
      return;
    }
    fillDetailedCollisionQA(collision, o2::analysis::CollisonCuts::kFlagOccupancy);
    fillDetailedCollisionQA(collision, o2::analysis::CollisonCuts::kAllpassed);
  }

  using GenMCCollisions = soa::Join<aod::McCollisions, aod::McCentFT0Ms, aod::McCentFT0Cs, aod::McCentFV0As, aod::MultsExtraMC>;
  using Run3MCCollisions = soa::Join<aod::McCollisions, aod::MultsExtraMC>;
  using Run2MCCollisions = soa::Join<aod::McCollisions, aod::MultsExtraMC>;
  using GenRecoCollisions = soa::Join<aod::ResoCollisionCandidates, aod::MultsExtra, aod::PVMults, aod::McCollisionLabels>;
  using BCsWithRCT = soa::Join<aod::BCsWithTimestamps, aod::BcSels>;

  /**
   * @brief Applies the configured RCT selection to a reconstructed collision
   *
   * The collision-level mode reads the RCT value stored in EVSEL. The BC-level
   * mode reads the value directly from the nominal BC referenced by collision.
   */
  template <typename CollisionType>
  bool isRecoRCTSelected(CollisionType const& collision)
  {
    if (!EventCuts.cfgEvtUseRCTFlagChecker.value) {
      return true;
    }
    if (!EventCuts.cfgEvtBCRCT.value) {
      return recoRCTChecker(collision);
    }
    const auto bc = collision.template bc_as<BCsWithRCT>();
    return recoRCTChecker(bc);
  }

  /// Apply the reconstructed Run 3 event and RCT selections. Detailed cut-flow
  /// QA is evaluated independently and never controls the event decision.
  template <typename CollisionType>
  bool isRun3CollisionSelected(CollisionType const& collision)
  {
    if (EventConfig.cfgFillDetailedQA) {
      fillDetailedRun3SelectionQA(collision);
    }
    if (!colCuts.isSelected(collision, EventConfig.cfgFillQA) || !isRecoRCTSelected(collision)) {
      return false;
    }
    if (EventConfig.cfgFillDetailedQA) {
      fillDetailedCollisionQA(collision, DetailedQARCTStage);
    }
    return true;
  }

  /**
   * @brief Fills generator-level resonance QA
   *
   * @tparam MCParticlesType Type of MC-particle group
   * @param mcParticles MC particles grouped by generator collision
   * @param generatorCentrality Generator or representative reconstructed centrality
   * @param multiplicity Generator-level charged-particle multiplicity
   * @param impactParameter Generator collision impact parameter
   * @param eventType INEL/INEL>0 category
   */
  template <typename MCParticlesType>
  void fillMCGenParticles(MCParticlesType const& mcParticles,
                          float generatorCentrality,
                          float multiplicity,
                          float impactParameter,
                          int eventType)
  {
    for (auto const& mcPart : mcParticles) {
      if (std::abs(mcPart.pdgCode()) != std::abs(GenCuts.pdgTruthMother.value)) {
        continue;
      }
      if (mcPart.y() <= GenCuts.cfgRapidityCutMinGen || mcPart.y() >= GenCuts.cfgRapidityCutMaxGen) {
        continue;
      }

      std::array<int, 2> daughterPDGs{-1, -1};
      const bool hasDaughters = mcPart.has_daughters();
      const auto daughterIds = mcPart.daughtersIds();
      if (hasDaughters) {
        auto daughter01 = mcParticles.rawIteratorAt(daughterIds[0] - mcParticles.offset());
        auto daughter02 = mcParticles.rawIteratorAt(daughterIds[1] - mcParticles.offset());
        daughterPDGs = {daughter01.pdgCode(), daughter02.pdgCode()};
      }
      if (GenCuts.isDaughterCheck) {
        const int daughter1 = std::abs(GenCuts.pdgTruthDaughter1.value);
        const int daughter2 = std::abs(GenCuts.pdgTruthDaughter2.value);
        const int firstDaughterPDG = std::abs(daughterPDGs[0]);
        const int secondDaughterPDG = std::abs(daughterPDGs[1]);
        // Match the configured decay to two distinct daughter slots. This is
        // important when both configured absolute PDGs are equal (e.g. K+K-):
        // one matching daughter must not satisfy both requirements.
        const bool matchesConfiguredDecay =
          (firstDaughterPDG == daughter1 && secondDaughterPDG == daughter2) ||
          (firstDaughterPDG == daughter2 && secondDaughterPDG == daughter1);
        if (!matchesConfiguredDecay) {
          continue;
        }
      }

      if (mcPart.pdgCode() > 0) {
        qaRegistry.fill(HIST("EventGen/h5ResonanceTruth"), eventType, mcPart.pt(), generatorCentrality, multiplicity, impactParameter);
      } else {
        qaRegistry.fill(HIST("EventGen/h5ResonanceTruthAnti"), eventType, mcPart.pt(), generatorCentrality, multiplicity, impactParameter);
      }
    }
  }

  /**
   * @brief Fills generated resonance parents for one reduced collision
   *
   * @tparam SelectedMCParticlesType Type of the selected MC-particle slice
   * @tparam MCParticlesType Type of the complete MC-particle table
   * @param reducedCollisionId Reduced collision referenced by the output rows
   * @param selectedParents Selected parent particles in the associated MC collision
   * @param mcParticles Complete MC-particle table used to resolve daughter indices
   *
   * The source MC-particle global index is persisted as a scalar in
   * ResoMCParents_001; it is not a relation requiring McParticles at merge time.
   */
  template <typename SelectedMCParticlesType, typename MCParticlesType>
  void fillMCParents(int64_t reducedCollisionId,
                     SelectedMCParticlesType const& selectedParents,
                     MCParticlesType const& mcParticles)
  {
    for (auto const& mcPart : selectedParents) {
      if (!GenCuts.cfgDoSignalLoss) {
        const int absPdg = std::abs(mcPart.pdgCode());
        if (absPdg == PdgProton || absPdg == PdgLambda0 || absPdg == PdgXiMinus || absPdg == PdgXi0 || absPdg == PdgOmegaMinus) {
          continue;
        }
      }

      std::array<int, 2> daughterPDGs{-1, -1};
      if (mcPart.has_daughters()) {
        const auto daughter1 = mcParticles.rawIteratorAt(mcPart.daughtersIds()[0] - mcParticles.offset());
        const auto daughter2 = mcParticles.rawIteratorAt(mcPart.daughtersIds()[1] - mcParticles.offset());
        daughterPDGs = {daughter1.pdgCode(), daughter2.pdgCode()};
      }
      reso2mcparents(reducedCollisionId,
                     mcPart.globalIndex(),
                     mcPart.pdgCode(),
                     daughterPDGs[0],
                     daughterPDGs[1],
                     mcPart.isPhysicalPrimary(),
                     mcPart.producedByGenerator(),
                     mcPart.pt(),
                     mcPart.px(),
                     mcPart.py(),
                     mcPart.pz(),
                     mcPart.y(),
                     mcPart.e(),
                     mcPart.statusCode());
    }
  }

  /**
   * @brief Returns the configured generator-level multiplicity estimator
   */
  template <typename MCCollision>
  float getMCMultiplicity(MCCollision const& mcCollision)
  {
    if (GenCuts.cfgGenMult05) {
      return mcCollision.multMCNParticlesEta05();
    }
    if (GenCuts.cfgGenMult10) {
      return mcCollision.multMCNParticlesEta10();
    }
    if (GenCuts.cfgGenMultFT0M) {
      return mcCollision.multMCFT0A() + mcCollision.multMCFT0C();
    }
    if (GenCuts.cfgGenMultFT0C) {
      return mcCollision.multMCFT0C();
    }
    if (GenCuts.cfgGenMultFV0A) {
      return mcCollision.multMCFV0A();
    }
    return -1.f;
  }

  /**
   * @brief Returns the configured generator-level centrality percentile
   */
  template <typename MCCollision>
  float getMCCentrality(MCCollision const& mcCollision) const
  {
    if (!GenCuts.cfgGenMultPercentile.value) {
      return 100.5f;
    }
    switch (multEstimator) {
      case CentralityFT0M:
        return mcCollision.centFT0M();
      case CentralityFT0C:
        return mcCollision.centFT0C();
      case CentralityFV0A:
        return mcCollision.centFV0A();
      default:
        return 100.5f;
    }
  }

  /**
   * @brief Fills the generator-only MC collision extension and optional source link
   *
   * Accepting only the generator collision here prevents reconstructed event
   * selection state from accidentally entering ResoMCCollisions_001.
   */
  template <typename MCCollision>
  void fillMCCollision001(MCCollision const& mcCollision)
  {
    const float mcMultiplicity = getMCMultiplicity(mcCollision);
    const bool inTrueVtx10 = std::abs(mcCollision.posZ()) < MCVertexZMax;
    // Keep the generator-level definition identical to MultsExtraMC:
    // at least one physical-primary charged particle within |eta| < 1.
    const bool isTrueINELgt0 = mcCollision.isInelGt0();
    resoMCCollisions001(inTrueVtx10,
                        isTrueINELgt0,
                        mcCollision.impactParameter(),
                        mcMultiplicity);
    if (GenCuts.cfgFillMCCollisionSoftLink) {
      resoMCCollisionIds(mcCollision.globalIndex());
    }
  }

  /// Write the Run 3 base collision, optional soft link, and scalar grouping key.
  template <typename CollisionType>
  void fillRun3Collision(CollisionType const& collision)
  {
    const bool isRecINELgt0 = collision.isInelGt0();
    centrality = centEst(collision);
    resoCollisions(collisionMultiplicity(collision),
                   collision.posX(),
                   collision.posY(),
                   collision.posZ(),
                   centrality,
                   dBz,
                   isRecINELgt0);
    resoCollisionColls(collision.globalIndex());
    resoCollisionGroups001(collision.globalIndex());
  }

  /// Write the Run 2 base collision, optional soft link, and scalar grouping key.
  template <typename CollisionType>
  void fillRun2Collision(CollisionType const& collision)
  {
    centrality = collision.centRun2V0M();
    // The configured Run 3 estimators are not available in the Run 2 input
    // schema. Preserve the previous zero-filled Run 2 behaviour.
    resoCollisions(0.f,
                   collision.posX(),
                   collision.posY(),
                   collision.posZ(),
                   centrality,
                   dBz,
                   0);
    resoCollisionColls(collision.globalIndex());
    resoCollisionGroups001(collision.globalIndex());
  }

  /**
   * @brief Processes Dummy
   *
   * @param collision Collision data
   */
  void processDummy(aod::Collisions const&)
  {
  }
  PROCESS_SWITCH(ResonanceModuleInitializer, processDummy, "process Dummy", true);

  /**
   * @brief Processes Run3 data
   *
   * @param collision Collision data
   * @param bc BC data
   */
  void processRun3(aod::ResoCollisionCandidates::iterator const& collision,
                   BCsWithRCT const&)
  {
    auto bc = collision.bc_as<BCsWithRCT>();
    initCCDB(bc);
    // Default event selection
    if (!isRun3CollisionSelected(collision)) {
      return;
    }
    if (EventConfig.cfgFillQA) {
      colCuts.fillQA(collision);
    }
    fillRun3Collision(collision);
  }
  PROCESS_SWITCH(ResonanceModuleInitializer, processRun3, "Default process for RUN3", false);

  /**
   * @brief Processes Run2 data
   *
   * @param collision Collision data
   * @param bc BC data
   */
  void processRun2(aod::ResoRun2CollisionCandidates::iterator const& collision,
                   aod::BCsWithRun2Info const&)
  {
    // auto bc = collision.bc_as<aod::BCsWithRun2Info>();
    // Default event selection
    if (!colCuts.isSelected(collision, EventConfig.cfgFillQA)) {
      return;
    }
    if (EventConfig.cfgFillQA) {
      colCuts.fillQARun2(collision);
    }
    fillRun2Collision(collision);
  }
  PROCESS_SWITCH(ResonanceModuleInitializer, processRun2, "process for RUN2", false);

  /**
   * @brief Processes generator-level MC event and resonance QA
   *
   * The original MC collision is the grouping key. MC particles are therefore
   * grouped automatically, while reconstructed collisions arrive as a 0..N
   * SmallGroup through their McCollisionLabels relation. This auxiliary
   * callback intentionally fills generator-level QA only and does not write
   * reduced AOD tables. RCT quality is evaluated through the generator
   * collision's associated BC because it is a run-condition property.
   */
  void processMCgen(GenMCCollisions::iterator const& mcCollision,
                    aod::McParticles const& mcParticles,
                    soa::SmallGroups<GenRecoCollisions> const& collisions,
                    BCsWithRCT const&)
  {
    auto bc = mcCollision.bc_as<BCsWithRCT>();
    initCCDB(bc);

    const auto getReconstructedCentrality = [&](auto const& collision) {
      return centEst(collision);
    };

    const float generatorCentrality = getMCCentrality(mcCollision);
    const float impactParameter = mcCollision.impactParameter();
    const float multiplicity = getMCMultiplicity(mcCollision);

    qaRegistry.fill(HIST("EventGen/hNEventsMC"), 0.5);
    if (GenCuts.isZvtxcutGen && std::abs(mcCollision.posZ()) > GenCuts.cutzvertexGen) {
      return;
    }
    qaRegistry.fill(HIST("EventGen/hNEventsMC"), 1.5);
    if (GenCuts.cfgGenBCRCT && !genRCTChecker(bc)) {
      return;
    }
    qaRegistry.fill(HIST("EventGen/hNEventsMC"), 2.5);
    qaRegistry.fill(HIST("EventGen/hZCollisionGen"), mcCollision.posZ());

    int eventType = 0;
    qaRegistry.fill(HIST("EventGen/hNEventsMC"), 3.5);
    if (GenCuts.checkIsTrueINELgt0 && mcCollision.isInelGt0()) {
      eventType = 1;
      qaRegistry.fill(HIST("EventGen/hNEventsMC"), 4.5);
    }

    bool hasSelectedRecoCollision = false;
    int largestNContributors = -1;
    float reconstructedCentrality = 100.5f;
    for (auto const& collision : collisions) {
      if (!isRecoRCTSelected(collision)) {
        continue;
      }
      if (!colCuts.isSelected(collision, false)) {
        continue;
      }
      const int nContributors = static_cast<int>(collision.multPVTotalContributors());
      if (nContributors > largestNContributors) {
        largestNContributors = nContributors;
        reconstructedCentrality = getReconstructedCentrality(collision);
      }
      hasSelectedRecoCollision = true;
    }

    if (GenCuts.cfgGenMultPercentile) {
      fillMCGenParticles(mcParticles, generatorCentrality, multiplicity, impactParameter, eventType);
      qaRegistry.fill(HIST("EventGen/h4MultCent_genMC"), eventType, generatorCentrality, multiplicity, impactParameter);
    } else {
      fillMCGenParticles(mcParticles, reconstructedCentrality, multiplicity, impactParameter, eventType);
      qaRegistry.fill(HIST("EventGen/h4MultCent_genMC"), eventType, reconstructedCentrality, multiplicity, impactParameter);
      qaRegistry.fill(HIST("EventGen/h2CentralityVsMultMC"), reconstructedCentrality, multiplicity);
    }
    if (hasSelectedRecoCollision) {
      qaRegistry.fill(HIST("EventGen/h4MultCent_recMC"), eventType, reconstructedCentrality, multiplicity, impactParameter);
    }
  }
  PROCESS_SWITCH(ResonanceModuleInitializer, processMCgen, "Process generator-level MC QA", false);

  /**
   * @brief Processes Run3 MC data
   *
   * Reconstructed collisions pass the same event selection and QA sequence as
   * Run 3 data before the reduced collision and its MC extension are written.
   *
   * @param collision Collision data
   * @param mcCollisions MC collisions
   * @param mcParticles MC particles used to fill reduced resonance parents
   */
  void processRun3MC(aod::ResoCollisionCandidatesMC::iterator const& collision,
                     Run3MCCollisions const&,
                     aod::McParticles const& mcParticles,
                     BCsWithRCT const&)
  {
    auto bc = collision.bc_as<BCsWithRCT>();
    initCCDB(bc);
    if (!isRun3CollisionSelected(collision)) {
      return;
    }
    if (EventConfig.cfgFillQA) {
      colCuts.fillQA(collision);
    }
    if (!collision.has_mcCollision()) {
      return;
    }
    const auto& mcCollision = collision.mcCollision_as<Run3MCCollisions>();

    // ResoMCCollisions_001 is a positional extension of ResoCollisions. When
    // enabled, ResoMCCollisionIds is written in the same callback and is 1:1.
    fillRun3Collision(collision);
    const int64_t reducedCollisionId = resoCollisions.lastIndex();
    fillMCCollision001(mcCollision);

    // A generator collision can be associated with more than one reconstructed
    // collision. Write one parent set for each reduced collision, as in the
    // legacy collision-wise producer, so ResoCollisionId remains unambiguous.
    auto selectedParents = selectedMCParticles->sliceBy(mcParticlesPerMcCollision, mcCollision.globalIndex());
    fillMCParents(reducedCollisionId, selectedParents, mcParticles);
  }
  PROCESS_SWITCH(ResonanceModuleInitializer, processRun3MC, "process MC for RUN3", false);

  /**
   * @brief Processes Run2 MC data
   *
   * Reconstructed collisions pass the same event selection and QA sequence as
   * Run 2 data before the reduced collision and its MC extension are written.
   *
   * @param collision Collision data
   */
  void processRun2MC(aod::ResoRun2CollisionCandidatesMC::iterator const& collision,
                     Run2MCCollisions const&)
  {
    if (!colCuts.isSelected(collision, EventConfig.cfgFillQA)) {
      return;
    }
    if (EventConfig.cfgFillQA) {
      colCuts.fillQARun2(collision);
    }
    if (!collision.has_mcCollision()) {
      return;
    }
    const auto& mcCollision = collision.mcCollision_as<Run2MCCollisions>();

    // Keep the base and MC extension one-to-one; the optional source link is
    // written by fillMCCollision001 in the same order.
    fillRun2Collision(collision);
    fillMCCollision001(mcCollision);
  }
  PROCESS_SWITCH(ResonanceModuleInitializer, processRun2MC, "process MC for RUN2", false);
};

/**
 * @brief Initializer for the resonance daughters producer
 *
 * This struct initializes and processes daughters for resonance studies.
 * It applies daughter selection criteria and fills QA histograms for daughter properties.
 */
struct ResonanceDaughterInitializer {
  enum class UltraMicroPidSpecies : uint8_t {
    Pion,
    Kaon,
    Proton
  };

  static constexpr int TrackSelectionNone = 0;
  static constexpr int TrackSelectionGlobal = 1;
  static constexpr int TrackSelectionGlobalWoPtEta = 2;
  static constexpr int TrackSelectionGlobalWoDCA = 3;
  static constexpr int TrackSelectionQuality = 4;
  static constexpr int TrackSelectionInAcceptance = 5;
  static constexpr int PairGateModeConfigured = 0;
  static constexpr int PairGateModeEither = 1;
  static constexpr float MomentumQuantizationScale = 1000.f;
  static constexpr std::size_t StoredMCRelationCount = 2;
  static constexpr std::size_t MaxCandidateDaughters = 3;

  /// Selected-candidate state and the optional global daughter-ID veto set.
  /// By default only candidate existence is recorded and daughter reuse is
  /// rejected later for each concrete pair through the stored trackId.
  struct SelectedCandidateDaughters {
    std::vector<int64_t> allDaughterIds;
    bool hasSelectedCandidate = false;
    bool useGlobalDaughterVeto = false;

    explicit SelectedCandidateDaughters(bool globalDaughterVeto = false)
      : useGlobalDaughterVeto(globalDaughterVeto)
    {
    }

    void addCandidate()
    {
      hasSelectedCandidate = true;
    }

    template <std::size_t DaughterCount>
    void addCandidate(std::array<int64_t, DaughterCount> const& daughterIds)
    {
      static_assert(DaughterCount <= MaxCandidateDaughters);
      allDaughterIds.insert(allDaughterIds.end(), daughterIds.begin(), daughterIds.end());
      hasSelectedCandidate = true;
    }

    void finalize()
    {
      if (!useGlobalDaughterVeto) {
        return;
      }
      std::sort(allDaughterIds.begin(), allDaughterIds.end());
      allDaughterIds.erase(std::unique(allDaughterIds.begin(), allDaughterIds.end()), allDaughterIds.end());
    }

    [[nodiscard]] bool accepts(int64_t trackId) const
    {
      if (!hasSelectedCandidate) {
        return false;
      }
      if (useGlobalDaughterVeto) {
        return !std::binary_search(allDaughterIds.begin(), allDaughterIds.end(), trackId);
      }
      return true;
    }
  };

  /// Producer-side track retention for the candidate types enabled by a callback.
  struct PairTrackSelection {
    SelectedCandidateDaughters v0Candidates;
    SelectedCandidateDaughters cascadeCandidates;
    bool useV0Candidates = false;
    bool useCascadeCandidates = false;
    bool useGlobalDaughterVeto = false;

    explicit PairTrackSelection(bool globalDaughterVeto = false)
      : v0Candidates(globalDaughterVeto),
        cascadeCandidates(globalDaughterVeto),
        useGlobalDaughterVeto(globalDaughterVeto)
    {
    }

    template <typename TrackType>
    bool operator()(TrackType const& track) const
    {
      if (!useGlobalDaughterVeto) {
        return (useV0Candidates && v0Candidates.hasSelectedCandidate) ||
               (useCascadeCandidates && cascadeCandidates.hasSelectedCandidate);
      }

      const auto trackId = static_cast<int64_t>(track.globalIndex());
      bool hasCandidate = false;
      if (useV0Candidates && v0Candidates.hasSelectedCandidate) {
        hasCandidate = true;
        if (!v0Candidates.accepts(trackId)) {
          return false;
        }
      }
      if (useCascadeCandidates && cascadeCandidates.hasSelectedCandidate) {
        hasCandidate = true;
        if (!cascadeCandidates.accepts(trackId)) {
          return false;
        }
      }
      return hasCandidate;
    }
  };

  struct KeepAllTracks {
    template <typename TrackType>
    bool operator()(TrackType const&) const
    {
      return true;
    }
  };

  UltraMicroPidSpecies ultraMicroPidSpecies = UltraMicroPidSpecies::Pion;
  bool warnedUltraMicroMomentumRange = false;

  Preslice<soa::Filtered<aod::ResoTrackCandidates>> tracksPerCollision = aod::track::collisionId;
  Preslice<soa::Filtered<aod::ResoTrackCandidatesMC>> tracksMCPerCollision = aod::track::collisionId;
  Preslice<aod::ResoV0Candidates> v0sPerCollision = aod::v0data::collisionId;
  Preslice<aod::ResoV0CandidatesMC> v0sMCPerCollision = aod::v0data::collisionId;
  Preslice<aod::ResoCascadesCandidates> cascadesPerCollision = aod::cascdata::collisionId;
  Preslice<aod::ResoCascadesCandidatesMC> cascadesMCPerCollision = aod::cascdata::collisionId;
  Produces<aod::ResoTracks> reso2trks;                                ///< Output table for resonance tracks
  Produces<aod::ResoTrackTracks> resoTrackTracks;                     ///< Output table for original track row IDs
  Produces<aod::ResoMicroTracks_001> reso2microtrks;                  ///< Output table for resonance microtracks
  Produces<aod::ResoMicroTrackTracks> resoMicroTrackTracks;           ///< Positional original-track soft links for microtracks
  Produces<aod::ResoMCMicroTracks_001> reso2mcmicrotrks;              ///< Positional MC extension for resonance microtracks
  Produces<aod::ResoUltraMicroTracks> reso2ultramicrotrks;            ///< Output table for resonance ultra-microtracks
  Produces<aod::ResoUltraMicroTrackTracks> resoUltraMicroTrackTracks; ///< Output table for original ultra-microtrack row IDs
  Produces<aod::ResoMCTracks> reso2mctracks;                          ///< Output table for MC resonance tracks
  Produces<aod::ResoV0s> reso2v0s;                                    ///< Output table for resonance V0s
  Produces<aod::ResoV0V0s> resoV0V0s;                                 ///< Output table for original V0 row IDs
  Produces<aod::ResoMCV0s> reso2mcv0s;                                ///< Output table for MC resonance V0s
  Produces<aod::ResoCascades> reso2cascades;                          ///< Output table for resonance cascades
  Produces<aod::ResoCascadeCascades> resoCascadeCascades;             ///< Output table for original cascade row IDs
  Produces<aod::ResoMCCascades> reso2mccascades;                      ///< Output table for MC resonance cascades

  // General daughter output options
  Configurable<bool> cfgFillQA{"cfgFillQA", false, "Fill QA histograms"};
  Configurable<bool> cfgDetailTrackQA{"cfgDetailTrackQA", false, "Fill detailed QA histograms for enabled track output tables"};

  // Track pre-selection and DCA cuts
  struct : ConfigurableGroup {
    Configurable<float> cfgCutEta{"cfgCutEta", 0.8f, "Eta range for tracks"};
    Configurable<float> cfgCutMinPt{"cfgCutMinPt", 0.1f, "Minimum pT for tracks (GeV/c)"};
    Configurable<float> cfgCutMaxPt{"cfgCutMaxPt", 999.0f, "Maximum pT for tracks (GeV/c)"};
    Configurable<float> pidnSigmaPreSelectionCut{"pidnSigmaPreSelectionCut", 5.0f, "TPC PID half-width around the configured species mean (loose preselection)"};
    Configurable<float> pidnSigmaPreSelectionCutTOF{"pidnSigmaPreSelectionCutTOF", 5.0f, "TOF PID half-width around the configured species mean (loose preselection)"};
    Configurable<float> pidnSigmaPreSelectionMeanPion{"pidnSigmaPreSelectionMeanPion", 0.000f, "Offset for TPC PID mean for pions"};
    Configurable<float> pidnSigmaPreSelectionMeanKaon{"pidnSigmaPreSelectionMeanKaon", 0.000f, "Offset for TPC PID mean for kaons"};
    Configurable<float> pidnSigmaPreSelectionMeanProton{"pidnSigmaPreSelectionMeanProton", 0.000f, "Offset for TPC PID mean for protons"};
    Configurable<float> pidnSigmaPreSelectionMeanTOFPion{"pidnSigmaPreSelectionMeanTOFPion", 0.000f, "Offset for TOF PID mean for pions"};
    Configurable<float> pidnSigmaPreSelectionMeanTOFKaon{"pidnSigmaPreSelectionMeanTOFKaon", 0.000f, "Offset for TOF PID mean for kaons"};
    Configurable<float> pidnSigmaPreSelectionMeanTOFProton{"pidnSigmaPreSelectionMeanTOFProton", 0.000f, "Offset for TOF PID mean for protons"};
    Configurable<bool> cfgUseTOFPIDPreSelection{"cfgUseTOFPIDPreSelection", false,
                                                "Apply the TOF PID cut to tracks with TOF; false ignores TOF PID, and tracks without TOF use TPC only"};
    Configurable<int> trackSelection{"trackSelection", 3, "Track selection: 0 -> No Cut, 1 -> kGlobalTrack, 2 -> kGlobalTrackWoPtEta, 3 -> kGlobalTrackWoDCA, 4 -> kQualityTracks, 5 -> kInAcceptanceTracks"};
    Configurable<double> cMaxDCArToPVcut{"cMaxDCArToPVcut", 0.5f, "Track DCAr cut to PV Maximum"};
    Configurable<double> cMaxDCAzToPVcut{"cMaxDCAzToPVcut", 1.0f, "Track DCAz cut to PV Maximum"};
    Configurable<double> cMinDCAzToPVcut{"cMinDCAzToPVcut", 0.0f, "Track DCAz cut to PV Minimum"};
    Configurable<bool> cfgApplyTightDCAPtDepSelection{"cfgApplyTightDCAPtDepSelection", true, "Apply the pT-dependent tight DCA selection"};
    Configurable<float> cfgTightDCAOffset{"cfgTightDCAOffset", 0.004f, "Constant term of the tight DCA threshold (cm)"};
    Configurable<float> cfgTightDCAPtCoefficient{"cfgTightDCAPtCoefficient", 0.013f, "Coefficient of the pT-dependent tight DCA threshold"};
    Configurable<float> cfgTightDCAPtPower{"cfgTightDCAPtPower", 1.f, "Power in tight DCA = offset + coefficient / pT^power"};
  } TrackCuts;

  // V0 and V0-daughter cuts : based on loose cuts in V0s production analysis in pp 13.6 TeV
  struct : ConfigurableGroup {
    Configurable<int> mincrossedrowsV0s{"mincrossedrowsV0s", 70, "Minimum crossed rows for V0 daughter tracks"};
    Configurable<double> cMinV0PosDCArToPVcut{"cMinV0PosDCArToPVcut", 0.05f, "V0 Positive Track DCAr cut to PV Minimum"};
    Configurable<double> cMinV0NegDCArToPVcut{"cMinV0NegDCArToPVcut", 0.05f, "V0 Negative Track DCAr cut to PV Minimum"};
    Configurable<double> cMinV0Radius{"cMinV0Radius", 0.9f, "Minimum V0 radius from PV"};
    Configurable<double> cMaxV0Radius{"cMaxV0Radius", 200.0f, "Maximum V0 radius from PV"};
    Configurable<double> cMinV0CosPA{"cMinV0CosPA", 0.95f, "Minimum V0 CosPA to PV"};
  } V0Cuts;

  // Additional K0s and Lambda0 selections not covered by V0Cuts: based on loose cuts in V0s production analysis in pp 13.6 TeV
  struct : ConfigurableGroup {
    Configurable<bool> cfgSecondaryRequire{"cfgSecondaryRequire", true, "Secondary cuts on/off"};
    Configurable<bool> cfgSecondaryArmenterosCut{"cfgSecondaryArmenterosCut", false, "cut on Armenteros-Podolanski graph"};
    Configurable<bool> cfgSecondaryCrossMassHypothesisCut{"cfgSecondaryCrossMassHypothesisCut", false, "Apply cut based on the lambda mass hypothesis"};
    Configurable<bool> cfgByPassDauPIDSelection{"cfgByPassDauPIDSelection", false, "Bypass TPC PID preselection for V0 daughters"};
    Configurable<float> cfgSecondaryDauDCAMax{"cfgSecondaryDauDCAMax", 1.0f, "Maximum DCA between V0 daughters"};
    Configurable<float> cfgSecondaryPtMin{"cfgSecondaryPtMin", 0.0f, "Minimum transverse momentum of Secondary"};
    Configurable<float> cfgSecondaryRapidityMax{"cfgSecondaryRapidityMax", 0.5f, "Maximum rapidity of Secondary"};
    Configurable<float> cfgSecondaryDCAtoPVMax{"cfgSecondaryDCAtoPVMax", 0.4f, "Maximum DCA Secondary to PV"};
    Configurable<float> cfgSecondaryProperLifetimeMax{"cfgSecondaryProperLifetimeMax", 40.f, "Maximum Secondary Lifetime"};
    Configurable<float> cfgSecondaryparamArmenterosCut{"cfgSecondaryparamArmenterosCut", 0.2f, "parameter for Armenteros Cut"};
    Configurable<float> cfgSecondaryMassWindow{"cfgSecondaryMassWindow", 0.03f, "Secondary inv mass selection window (GeV/c^2)"};
    Configurable<float> cfgSecondaryCrossMassCutWindow{"cfgSecondaryCrossMassCutWindow", 0.02f, "Secondary inv mass selection window with (anti)lambda hypothesis (GeV/c^2)"};
  } SecondaryCuts;

  // Cascade and cascade-daughter cuts: based on loose cuts in cascades production analysis in pp 13.6 TeV
  struct : ConfigurableGroup {
    Configurable<int> cfgMinCrossedRowsCascBach{"cfgMinCrossedRowsCascBach", 50, "min crossed rows for bachelor track from cascade"};
    Configurable<float> cMinCascBachDCArToPVcut{"cMinCascBachDCArToPVcut", 0.05f, "Cascade Bachelor Track DCAr cut to PV Minimum"};
    Configurable<float> cMaxCascBachDCArToPVcut{"cMaxCascBachDCArToPVcut", 200.0f, "Cascade Bachelor Track DCAr cut to PV Maximum"};
    Configurable<float> cMaxCascDCAV0Daughters{"cMaxCascDCAV0Daughters", 0.5f, "Cascade DCA between V0 daughters Maximum"};
    Configurable<float> cMaxCascDCACascDaughters{"cMaxCascDCACascDaughters", 1.2f, "Cascade DCA between Casc daughters Maximum"};
    Configurable<float> cMinCascV0CosPA{"cMinCascV0CosPA", 0.98f, "Minimum Cascade V0 CosPA to PV"};
    Configurable<float> cMaxCascV0Radius{"cMaxCascV0Radius", 200.0f, "Maximum Cascade V0 radius from PV"};
    Configurable<float> cMinCascV0Radius{"cMinCascV0Radius", 0.4f, "Minimum Cascade V0 radius from PV"};
    Configurable<float> cMinCascRadius{"cMinCascRadius", 0.4f, "Minimum Cascade radius from PV"};
    Configurable<float> cMaxCascRadius{"cMaxCascRadius", 200.0f, "Maximum Cascade radius from PV"};
    Configurable<float> cMinCascCosPA{"cMinCascCosPA", 0.99, "Minimum Cascade CosPA to PV"};
    Configurable<float> cMaxXiMassWindow{"cMaxXiMassWindow", 0.02f, "Xi mass Window (GeV/c^2)"};
  } CascadeCuts;

  // Derived dataset selections
  struct : ConfigurableGroup {
    Configurable<bool> cfgFillPionTracks{"cfgFillPionTracks", false, "Apply the pion PID filter to every enabled track table"};
    Configurable<bool> cfgFillKaonTracks{"cfgFillKaonTracks", false, "Apply the kaon PID filter to every enabled track table"};
    Configurable<bool> cfgFillProtonTracks{"cfgFillProtonTracks", false, "Apply the proton PID filter to every enabled track table"};
    Configurable<bool> cfgFillK0s{"cfgFillK0s", true, "Fill K0s"};
    Configurable<bool> cfgFillLambda0{"cfgFillLambda0", true, "Fill Lambda0"};
    Configurable<int> cfgPairGateMode{"cfgPairGateMode", 0,
                                      "Combined *WithPairGate mode: 0 -> require every enabled candidate type, "
                                      "1 -> require a selected V0 or cascade"};
    Configurable<bool> cfgBypassNoPairV0s{"cfgBypassNoPairV0s", false,
                                          "Require a selected V0 and at least one selected collision track in pair-gate mode 0"};
    Configurable<bool> cfgBypassNoPairCascades{"cfgBypassNoPairCascades", true,
                                               "Require a selected cascade and at least one selected collision track in pair-gate mode 0"};
    Configurable<bool> cfgGlobalDaughterVeto{"cfgGlobalDaughterVeto", false,
                                             "Pair-gate callbacks only: reject a track if its original ID is a daughter "
                                             "of any selected V0 or cascade; false stores all selected collision tracks "
                                             "for mandatory candidate-local rejection through trackId"};
    Configurable<bool> cfgFillMicroTracks{"cfgFillMicroTracks", false, "Fill micro tracks"};
    Configurable<bool> cfgFillUltraMicroTracks{"cfgFillUltraMicroTracks", false, "Fill ultra micro tracks"};
    Configurable<bool> cfgBypassTrackFill{"cfgBypassTrackFill", false, "Bypass the full ResoTracks table fill"};
    Configurable<bool> cfgBypassTrackIndexFill{"cfgBypassTrackIndexFill", false,
                                               "Bypass optional source-object soft-link side tables; "
                                               "ResoMicroTracks_001 keeps its inline scalar trackId, but Full/Ultra "
                                               "outputs then cannot perform candidate-local daughter rejection"};
  } FilterForDerivedTables;

  // Track selection filter based on configuration
  Filter trackFilter = (TrackCuts.trackSelection.node() == TrackSelectionNone) ||
                       ((TrackCuts.trackSelection.node() == TrackSelectionGlobal) && requireGlobalTrackInFilter()) ||                                          // kGlobalTrack = kQualityTracks | kPrimaryTracks | kInAcceptanceTracks
                       ((TrackCuts.trackSelection.node() == TrackSelectionGlobalWoPtEta) && requireGlobalTrackWoPtEtaInFilter()) ||                            // kGlobalTrackWoPtEta = kQualityTracks | kPrimaryTracks
                       ((TrackCuts.trackSelection.node() == TrackSelectionGlobalWoDCA) && requireGlobalTrackWoDCAInFilter()) ||                                // kGlobalTrackWoDCA = kQualityTracks | kInAcceptanceTracks
                       ((TrackCuts.trackSelection.node() == TrackSelectionQuality) && requireQualityTracksInFilter()) ||                                       // kQualityTracks = kQualityTracksITS | kQualityTracksTPC
                       ((TrackCuts.trackSelection.node() == TrackSelectionInAcceptance) && requireTrackCutInFilter(TrackSelectionFlags::kInAcceptanceTracks)); // kInAcceptanceTracks = kPtRange | kEtaRange
  Filter trackKinematicsFilter = nabs(aod::track::eta) < TrackCuts.cfgCutEta &&
                                 aod::track::pt >= TrackCuts.cfgCutMinPt &&
                                 aod::track::pt <= TrackCuts.cfgCutMaxPt;
  HistogramRegistry qaRegistry{"QAHistos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // The daughter task treats fIndexResoCollisions as a scalar v001 row ID. The
  // data-model relation keeps its legacy v000 default target; v001 consumers
  // must use resoCollision_as<T>() with the exact bound v001 table type for
  // dereferencing (aod::ResoCollisions_001 when the parent is not joined).
  // Original-collision grouping is provided by the scalar-only version-1
  // ResoCollisionGroups_001 table, avoiding a hard source-AO2D relation.
  // Collision mappings are always written; cfgBypassCollIndexFill is retained
  // only so existing configuration files remain accepted.
  using ResoCollisionWithIndex = soa::Join<aod::ResoCollisions_001, aod::ResoCollisionColls>;
  using SelectedResoCollisions = soa::Join<aod::ResoCollisions_001, aod::ResoCollisionGroups_001>;
  PresliceUnsorted<SelectedResoCollisions> reducedCollisionsPerOriginalCollision = aod::resocollisiongroup001::originalCollisionId;

  /**
   * @brief Initializes the task
   *
   * @param context Initialization context
   */
  void init(InitContext&)
  {
    if (FilterForDerivedTables.cfgPairGateMode.value < PairGateModeConfigured ||
        FilterForDerivedTables.cfgPairGateMode.value > PairGateModeEither) {
      LOGF(fatal, "cfgPairGateMode must be 0 (configured gates) or 1 (V0 or cascade)");
    }
    const bool useEitherPairGate = FilterForDerivedTables.cfgPairGateMode.value == PairGateModeEither;
    const bool processTrackDataEnabled = doprocessData || doprocessDataHybrid || doprocessDataWithPairGate ||
                                         doprocessDataWithV0PairGate || doprocessDataWithCascPairGate;
    const bool processTrackMCEnabled = doprocessMC || doprocessMCWithPairGate ||
                                       doprocessMCWithV0PairGate || doprocessMCWithCascPairGate;
    const bool pairTrackProcessEnabled = doprocessDataWithPairGate || doprocessDataWithV0PairGate ||
                                         doprocessDataWithCascPairGate || doprocessMCWithPairGate ||
                                         doprocessMCWithV0PairGate || doprocessMCWithCascPairGate;
    const bool processV0DataEnabled = doprocessV0Data;
    const bool processCascDataEnabled = doprocessCascData;
    const bool anyDataProcessEnabled = processTrackDataEnabled || processV0DataEnabled || processCascDataEnabled;
    const bool anyMCProcessEnabled = processTrackMCEnabled || doprocessV0MC || doprocessCascMC;
    const int enabledTrackProcesses = static_cast<int>(doprocessData) +
                                      static_cast<int>(doprocessDataHybrid) +
                                      static_cast<int>(doprocessDataWithPairGate) +
                                      static_cast<int>(doprocessDataWithV0PairGate) +
                                      static_cast<int>(doprocessDataWithCascPairGate) +
                                      static_cast<int>(doprocessMC) +
                                      static_cast<int>(doprocessMCWithPairGate) +
                                      static_cast<int>(doprocessMCWithV0PairGate) +
                                      static_cast<int>(doprocessMCWithCascPairGate);

    if (enabledTrackProcesses > 1) {
      LOGF(fatal, "Only one track process can be enabled in ResonanceDaughterInitializer");
    }
    if (pairTrackProcessEnabled &&
        FilterForDerivedTables.cfgBypassTrackFill &&
        !FilterForDerivedTables.cfgFillMicroTracks &&
        !FilterForDerivedTables.cfgFillUltraMicroTracks) {
      LOGF(fatal, "A pair-gate process requires at least one enabled Full, Micro, or UltraMicro track output");
    }
    if (pairTrackProcessEnabled && !FilterForDerivedTables.cfgGlobalDaughterVeto &&
        FilterForDerivedTables.cfgBypassTrackIndexFill &&
        (!FilterForDerivedTables.cfgBypassTrackFill || FilterForDerivedTables.cfgFillUltraMicroTracks)) {
      LOGF(warn,
           "Default pair-gate mode defers daughter reuse rejection to analysis trackId comparisons, "
           "but cfgBypassTrackIndexFill removes that ID from Full/Ultra track outputs");
    }
    if (useEitherPairGate &&
        (doprocessDataWithV0PairGate || doprocessDataWithCascPairGate ||
         doprocessMCWithV0PairGate || doprocessMCWithCascPairGate)) {
      LOGF(fatal, "cfgPairGateMode 1 requires the combined processDataWithPairGate or processMCWithPairGate callback");
    }
    if (static_cast<int>(doprocessV0Data) + static_cast<int>(doprocessV0MC) > 1) {
      LOGF(fatal, "Only one V0 process can be enabled in ResonanceDaughterInitializer");
    }
    if (static_cast<int>(doprocessCascData) + static_cast<int>(doprocessCascMC) > 1) {
      LOGF(fatal, "Only one cascade process can be enabled in ResonanceDaughterInitializer");
    }
    if ((doprocessData || doprocessDataHybrid || doprocessMC) &&
        (useEitherPairGate || FilterForDerivedTables.cfgBypassNoPairV0s || FilterForDerivedTables.cfgBypassNoPairCascades ||
         FilterForDerivedTables.cfgGlobalDaughterVeto)) {
      LOGF(warn, "Pair-gate options are ignored by processData/processDataHybrid/processMC; enable the matching *WithPairGate process to apply them");
    }
    const auto validatePairGateOutputs = [&](bool pairProcessEnabled,
                                             bool v0OutputEnabled,
                                             bool cascadeOutputEnabled,
                                             char const* processName) {
      if (!pairProcessEnabled) {
        return;
      }
      if (useEitherPairGate) {
        if (!v0OutputEnabled || !cascadeOutputEnabled) {
          LOGF(fatal, "%s with pair-gate mode 1 requires both V0 and cascade output processes", processName);
        }
      } else {
        if (!FilterForDerivedTables.cfgBypassNoPairV0s && !FilterForDerivedTables.cfgBypassNoPairCascades) {
          LOGF(fatal, "%s with pair-gate mode 0 requires at least one enabled V0/cascade gate", processName);
        }
        if (FilterForDerivedTables.cfgBypassNoPairV0s && !v0OutputEnabled) {
          LOGF(fatal, "%s requires a V0 output process when cfgBypassNoPairV0s is enabled", processName);
        }
        if (FilterForDerivedTables.cfgBypassNoPairCascades && !cascadeOutputEnabled) {
          LOGF(fatal, "%s requires a cascade output process when cfgBypassNoPairCascades is enabled", processName);
        }
      }
    };
    validatePairGateOutputs(doprocessDataWithPairGate,
                            processV0DataEnabled,
                            processCascDataEnabled,
                            "processDataWithPairGate");
    validatePairGateOutputs(doprocessMCWithPairGate,
                            doprocessV0MC,
                            doprocessCascMC,
                            "processMCWithPairGate");
    if (doprocessDataWithV0PairGate && !processV0DataEnabled) {
      LOGF(fatal, "processDataWithV0PairGate requires processV0Data");
    }
    if (doprocessDataWithCascPairGate && !processCascDataEnabled) {
      LOGF(fatal, "processDataWithCascPairGate requires processCascData");
    }
    if (doprocessMCWithV0PairGate && !doprocessV0MC) {
      LOGF(fatal, "processMCWithV0PairGate requires processV0MC");
    }
    if (doprocessMCWithCascPairGate && !doprocessCascMC) {
      LOGF(fatal, "processMCWithCascPairGate requires processCascMC");
    }
    if (!std::isfinite(TrackCuts.cfgCutMinPt.value) ||
        !std::isfinite(TrackCuts.cfgCutMaxPt.value) ||
        TrackCuts.cfgCutMinPt.value < 0.f ||
        TrackCuts.cfgCutMaxPt.value <= TrackCuts.cfgCutMinPt.value) {
      LOGF(fatal, "Track pT limits must be finite and satisfy 0 <= cfgCutMinPt < cfgCutMaxPt");
    }
    if (!std::isfinite(TrackCuts.cfgCutEta.value) || TrackCuts.cfgCutEta.value <= 0.f) {
      LOGF(fatal, "cfgCutEta must be finite and positive");
    }
    if (TrackCuts.trackSelection.value < TrackSelectionNone || TrackCuts.trackSelection.value > TrackSelectionInAcceptance) {
      LOGF(fatal, "trackSelection must be in [0, 5]");
    }
    if (!std::isfinite(TrackCuts.cMaxDCArToPVcut.value) ||
        !std::isfinite(TrackCuts.cMaxDCAzToPVcut.value) ||
        !std::isfinite(TrackCuts.cMinDCAzToPVcut.value) ||
        TrackCuts.cMaxDCArToPVcut.value < 0. ||
        TrackCuts.cMinDCAzToPVcut.value < 0. ||
        TrackCuts.cMaxDCAzToPVcut.value < TrackCuts.cMinDCAzToPVcut.value) {
      LOGF(fatal, "Track DCA limits must be finite and satisfy 0 <= cMinDCAzToPVcut <= cMaxDCAzToPVcut and cMaxDCArToPVcut >= 0");
    }
    if (!std::isfinite(TrackCuts.pidnSigmaPreSelectionCut.value) ||
        TrackCuts.pidnSigmaPreSelectionCut.value < 0.f) {
      LOGF(fatal, "pidnSigmaPreSelectionCut must be finite and non-negative");
    }
    const std::array tpcPidMeans{TrackCuts.pidnSigmaPreSelectionMeanPion.value,
                                 TrackCuts.pidnSigmaPreSelectionMeanKaon.value,
                                 TrackCuts.pidnSigmaPreSelectionMeanProton.value};
    if (!std::all_of(tpcPidMeans.begin(), tpcPidMeans.end(), [](float mean) {
          return std::isfinite(mean);
        })) {
      LOGF(fatal, "All TPC PID preselection means must be finite");
    }
    if (TrackCuts.cfgUseTOFPIDPreSelection.value) {
      if (!std::isfinite(TrackCuts.pidnSigmaPreSelectionCutTOF.value) ||
          TrackCuts.pidnSigmaPreSelectionCutTOF.value < 0.f) {
        LOGF(fatal, "pidnSigmaPreSelectionCutTOF must be finite and non-negative when TOF PID selection is enabled");
      }
      const std::array tofPidMeans{TrackCuts.pidnSigmaPreSelectionMeanTOFPion.value,
                                   TrackCuts.pidnSigmaPreSelectionMeanTOFKaon.value,
                                   TrackCuts.pidnSigmaPreSelectionMeanTOFProton.value};
      if (!std::all_of(tofPidMeans.begin(), tofPidMeans.end(), [](float mean) {
            return std::isfinite(mean);
          })) {
        LOGF(fatal, "All TOF PID preselection means must be finite when TOF PID selection is enabled");
      }
    }
    if (TrackCuts.cfgApplyTightDCAPtDepSelection.value &&
        (!std::isfinite(TrackCuts.cfgTightDCAOffset.value) ||
         !std::isfinite(TrackCuts.cfgTightDCAPtCoefficient.value) ||
         !std::isfinite(TrackCuts.cfgTightDCAPtPower.value) ||
         TrackCuts.cfgTightDCAOffset.value < 0.f ||
         TrackCuts.cfgTightDCAPtCoefficient.value < 0.f ||
         TrackCuts.cfgTightDCAPtPower.value < 0.f)) {
      LOGF(fatal, "Tight-DCA offset, pT coefficient, and power must be finite and non-negative");
    }
    if (!std::isfinite(V0Cuts.cMinV0CosPA.value) ||
        V0Cuts.cMinV0CosPA.value < -1. || V0Cuts.cMinV0CosPA.value >= 1.) {
      LOGF(fatal, "cMinV0CosPA must be finite and satisfy -1 <= cMinV0CosPA < 1");
    }
    if (!std::isfinite(CascadeCuts.cMinCascCosPA.value) ||
        CascadeCuts.cMinCascCosPA.value < -1. || CascadeCuts.cMinCascCosPA.value >= 1.) {
      LOGF(fatal, "cMinCascCosPA must be finite and satisfy -1 <= cMinCascCosPA < 1");
    }
    if (!std::isfinite(CascadeCuts.cMinCascRadius.value) ||
        !std::isfinite(CascadeCuts.cMaxCascRadius.value) ||
        CascadeCuts.cMinCascRadius.value < 0. ||
        CascadeCuts.cMaxCascRadius.value <= CascadeCuts.cMinCascRadius.value) {
      LOGF(fatal, "Cascade radius limits must be finite and satisfy 0 <= cMinCascRadius < cMaxCascRadius");
    }
    if (FilterForDerivedTables.cfgFillUltraMicroTracks) {
      int enabledUltraMicroSpecies = 0;
      enabledUltraMicroSpecies += FilterForDerivedTables.cfgFillPionTracks ? 1 : 0;
      enabledUltraMicroSpecies += FilterForDerivedTables.cfgFillKaonTracks ? 1 : 0;
      enabledUltraMicroSpecies += FilterForDerivedTables.cfgFillProtonTracks ? 1 : 0;
      if (enabledUltraMicroSpecies != 1) {
        LOGF(fatal, "Exactly one of cfgFillPionTracks, cfgFillKaonTracks, or cfgFillProtonTracks must be enabled when filling ultra-micro tracks");
      }
      if (TrackCuts.pidnSigmaPreSelectionCut.value > o2::aod::resoultramicrodaughter::PidNSigma::MaxNSigma) {
        LOGF(fatal, "Ultra-micro TPC PID encoding requires pidnSigmaPreSelectionCut <= 5");
      }
      if (TrackCuts.cfgUseTOFPIDPreSelection.value &&
          TrackCuts.pidnSigmaPreSelectionCutTOF.value > o2::aod::resoultramicrodaughter::PidNSigma::MaxNSigma) {
        LOGF(fatal, "Ultra-micro TOF PID encoding requires pidnSigmaPreSelectionCutTOF <= 5");
      }
      if (FilterForDerivedTables.cfgFillKaonTracks) {
        ultraMicroPidSpecies = UltraMicroPidSpecies::Kaon;
      } else if (FilterForDerivedTables.cfgFillProtonTracks) {
        ultraMicroPidSpecies = UltraMicroPidSpecies::Proton;
      } else {
        ultraMicroPidSpecies = UltraMicroPidSpecies::Pion;
      }
    }

    if (cfgFillQA) {
      AxisSpec idxAxis = {8, 0.0, 8.0, "Cumulative selection stage"};
      AxisSpec ptAxis = {100, 0.0f, 10.0f, "#it{p}_{T} (GeV/#it{c})"};
      AxisSpec dcaPtAxis = {300, 0.f, 30.f, "#it{p}_{T} (GeV/#it{c})"};
      AxisSpec etaAxis = {100, -1.0f, 1.0f, "#eta"};
      // Keep the azimuthal QA bins aligned with the 18 TPC sectors (10 bins per sector).
      AxisSpec phiAxis = {180, 0.0f, TwoPI, "#phi"};
      AxisSpec dcaXYAxis = {1000, -0.5, 0.5, "DCA_{xy} (cm)"};
      AxisSpec dcaZAxis = {1000, -0.5, 0.5, "DCA_{z} (cm)"};
      // Tiny PID has 255 discrete values from -6.35 to 6.35 in 0.05 steps.
      // Extend the histogram limits by half a step so every decoded value is
      // located at a bin centre instead of on a floating-point bin boundary.
      constexpr int NTinyTPCPidBins = aod::pidtpc_tiny::binning::nbins + 1;
      constexpr float HalfTinyTPCPidBinWidth = 0.5f * aod::pidtpc_tiny::binning::bin_width;
      constexpr int NTinyTOFPidBins = aod::pidtof_tiny::binning::nbins + 1;
      constexpr float HalfTinyTOFPidBinWidth = 0.5f * aod::pidtof_tiny::binning::bin_width;
      AxisSpec nSigmaTPCAxis = {NTinyTPCPidBins,
                                aod::pidtpc_tiny::binning::binned_min - HalfTinyTPCPidBinWidth,
                                aod::pidtpc_tiny::binning::binned_max + HalfTinyTPCPidBinWidth,
                                "TPC N#sigma"};
      AxisSpec nSigmaTOFAxis = {NTinyTOFPidBins,
                                aod::pidtof_tiny::binning::binned_min - HalfTinyTOFPidBinWidth,
                                aod::pidtof_tiny::binning::binned_max + HalfTinyTOFPidBinWidth,
                                "TOF N#sigma"};

      if (processTrackDataEnabled || processTrackMCEnabled) {
        qaRegistry.add("QA/hGoodTrackIndices", "hGoodTrackIndices", kTH1D, {idxAxis});
        auto trackSelection = qaRegistry.get<TH1>(HIST("QA/hGoodTrackIndices"));
        trackSelection->GetXaxis()->SetBinLabel(1, "Before DCA cuts");
        trackSelection->GetXaxis()->SetBinLabel(2, "Finite DCA, |DCA_{xy}| #leq max");
        trackSelection->GetXaxis()->SetBinLabel(3, "min #leq |DCA_{z}| #leq max");
        trackSelection->GetXaxis()->SetBinLabel(8, "Pass DCA selection");
        if (processTrackMCEnabled) {
          qaRegistry.add("QA/hGoodMCTrackIndices", "hGoodMCTrackIndices", kTH1D, {idxAxis});
          qaRegistry.get<TH1>(HIST("QA/hGoodMCTrackIndices"))->GetXaxis()->SetBinLabel(1, "MC path: before DCA cuts");
        }
        if (cfgDetailTrackQA) {
          if (!FilterForDerivedTables.cfgBypassTrackFill) {
            qaRegistry.add("QA/h4TrackPtEtaPhi", "ResoTracks pT, eta, phi", kTHnSparseD, {ptAxis, etaAxis, phiAxis});
            qaRegistry.add("QA/h2TrackDCAxyVsPt", "ResoTracks DCAxy vs pT", kTH2D, {dcaPtAxis, dcaXYAxis});
            qaRegistry.add("QA/h2TrackDCAzVsPt", "ResoTracks DCAz vs pT", kTH2D, {dcaPtAxis, dcaZAxis});
            qaRegistry.add("QA/h4TrackTPCnSigma", "ResoTracks TPC nSigma Pi, Ka, Pr as pT", kTHnSparseD, {ptAxis, nSigmaTPCAxis, nSigmaTPCAxis, nSigmaTPCAxis});
            qaRegistry.add("QA/h4TrackTOFnSigma", "ResoTracks TOF nSigma Pi, Ka, Pr as pT", kTHnSparseD, {ptAxis, nSigmaTOFAxis, nSigmaTOFAxis, nSigmaTOFAxis});
          }
          if (FilterForDerivedTables.cfgFillMicroTracks) {
            qaRegistry.add("QA/h4MicroTrackPtEtaPhi", "ResoMicroTracks pT, eta, phi", kTHnSparseD, {ptAxis, etaAxis, phiAxis});
            qaRegistry.add("QA/h2MicroTrackDCAxyVsPt", "ResoMicroTracks DCAxy vs pT", kTH2D, {dcaPtAxis, dcaXYAxis});
            qaRegistry.add("QA/h2MicroTrackDCAzVsPt", "ResoMicroTracks DCAz vs pT", kTH2D, {dcaPtAxis, dcaZAxis});
            qaRegistry.add("QA/h4MicroTrackTPCnSigma", "ResoMicroTracks TPC nSigma Pi, Ka, Pr as pT", kTHnSparseD, {ptAxis, nSigmaTPCAxis, nSigmaTPCAxis, nSigmaTPCAxis});
            qaRegistry.add("QA/h4MicroTrackTOFnSigma", "ResoMicroTracks TOF nSigma Pi, Ka, Pr as pT", kTHnSparseD, {ptAxis, nSigmaTOFAxis, nSigmaTOFAxis, nSigmaTOFAxis});
          }
          if (FilterForDerivedTables.cfgFillUltraMicroTracks) {
            qaRegistry.add("QA/h4UltraMicroTrackPtEtaPhi", "ResoUltraMicroTracks pT, eta, phi", kTHnSparseD, {ptAxis, etaAxis, phiAxis});
            qaRegistry.add("QA/h2UltraMicroTrackDCAxyVsPt", "ResoUltraMicroTracks DCAxy vs pT", kTH2D, {dcaPtAxis, dcaXYAxis});
            qaRegistry.add("QA/h2UltraMicroTrackDCAzVsPt", "ResoUltraMicroTracks DCAz vs pT", kTH2D, {dcaPtAxis, dcaZAxis});
            qaRegistry.add("QA/h4UltraMicroTrackTPCnSigma", "ResoUltraMicroTracks TPC nSigma Pi, Ka, Pr as pT", kTHnSparseD, {ptAxis, nSigmaTPCAxis, nSigmaTPCAxis, nSigmaTPCAxis});
            qaRegistry.add("QA/h4UltraMicroTrackTOFnSigma", "ResoUltraMicroTracks TOF nSigma Pi, Ka, Pr as pT", kTHnSparseD, {ptAxis, nSigmaTOFAxis, nSigmaTOFAxis, nSigmaTOFAxis});
          }
        }
      }

      if (processV0DataEnabled || doprocessV0MC) {
        qaRegistry.add("QA/hGoodV0Indices", "hGoodV0Indices", kTH1D, {idxAxis});
        auto v0Selection = qaRegistry.get<TH1>(HIST("QA/hGoodV0Indices"));
        v0Selection->GetXaxis()->SetBinLabel(1, "Before V0 cuts");
        v0Selection->GetXaxis()->SetBinLabel(2, "Daughter TPC rows");
        v0Selection->GetXaxis()->SetBinLabel(3, "Daughter |DCA_{xy}| to PV");
        v0Selection->GetXaxis()->SetBinLabel(4, "V0 radius window");
        v0Selection->GetXaxis()->SetBinLabel(5, "V0 cosPA cut passed");
        if (doprocessV0MC) {
          qaRegistry.add("QA/hGoodMCV0Indices", "hGoodMCV0Indices", kTH1D, {idxAxis});
          qaRegistry.get<TH1>(HIST("QA/hGoodMCV0Indices"))->GetXaxis()->SetBinLabel(1, "Reco V0 cuts passed (MC)");
        }
      }

      if (processCascDataEnabled || doprocessCascMC) {
        qaRegistry.add("QA/hGoodCascIndices", "hGoodCascIndices", kTH1D, {idxAxis});
        auto cascadeSelection = qaRegistry.get<TH1>(HIST("QA/hGoodCascIndices"));
        cascadeSelection->GetXaxis()->SetBinLabel(1, "Before cascade cuts");
        cascadeSelection->GetXaxis()->SetBinLabel(2, "Bachelor TPC rows");
        cascadeSelection->GetXaxis()->SetBinLabel(3, "Bachelor DCA_{xy} window");
        cascadeSelection->GetXaxis()->SetBinLabel(4, "V0/cascade daughter DCA");
        cascadeSelection->GetXaxis()->SetBinLabel(5, "Cascade/V0 cosPA");
        cascadeSelection->GetXaxis()->SetBinLabel(6, "V0 radius window");
        cascadeSelection->GetXaxis()->SetBinLabel(7, "Cascade radius window");
        cascadeSelection->GetXaxis()->SetBinLabel(8, "#Xi mass window");
        if (doprocessCascMC) {
          qaRegistry.add("QA/hGoodMCCascIndices", "hGoodMCCascIndices", kTH1D, {idxAxis});
          qaRegistry.get<TH1>(HIST("QA/hGoodMCCascIndices"))->GetXaxis()->SetBinLabel(1, "Reco cascade cuts passed (MC)");
        }
      }
    }
    if (processTrackDataEnabled || processTrackMCEnabled) {
      LOGF(info, "ResonanceDaughterInitializer initialized with tracks");
    }
    if (processV0DataEnabled || doprocessV0MC) {
      LOGF(info, "ResonanceDaughterInitializer initialized with V0s");
    }
    if (processCascDataEnabled || doprocessCascMC) {
      LOGF(info, "ResonanceDaughterInitializer initialized with cascades");
    }

    // Check if none of the processes are enabled
    if (!doprocessDummy && !anyDataProcessEnabled && !anyMCProcessEnabled) {
      LOGF(fatal, "ResonanceDaughterInitializer not initialized, enable at least one process");
    }
  }
  static bool passesCenteredPID(float nSigma, float mean, float cut)
  {
    return std::isfinite(nSigma) && std::abs(nSigma - mean) < cut;
  }

  bool passesPIDPreSelection(float tpcNSigma,
                             float tofNSigma,
                             bool hasTOF,
                             float tpcMean,
                             float tofMean,
                             float tpcCut,
                             float tofCut) const
  {
    if (!passesCenteredPID(tpcNSigma, tpcMean, tpcCut)) {
      return false;
    }
    return !TrackCuts.cfgUseTOFPIDPreSelection.value ||
           !hasTOF || passesCenteredPID(tofNSigma, tofMean, tofCut);
  }

  template <typename T>
  bool filterTrack(T const& track)
  {
    // if no selection is requested, return true
    if (!FilterForDerivedTables.cfgFillPionTracks && !FilterForDerivedTables.cfgFillKaonTracks && !FilterForDerivedTables.cfgFillProtonTracks) {
      return true;
    }
    const float tpcCut = TrackCuts.pidnSigmaPreSelectionCut.value;
    const float tofCut = TrackCuts.pidnSigmaPreSelectionCutTOF.value;
    if (FilterForDerivedTables.cfgFillPionTracks &&
        passesPIDPreSelection(track.tpcNSigmaPi(), track.tofNSigmaPi(), track.hasTOF(),
                              TrackCuts.pidnSigmaPreSelectionMeanPion.value,
                              TrackCuts.pidnSigmaPreSelectionMeanTOFPion.value, tpcCut, tofCut)) {
      return true;
    }
    if (FilterForDerivedTables.cfgFillKaonTracks &&
        passesPIDPreSelection(track.tpcNSigmaKa(), track.tofNSigmaKa(), track.hasTOF(),
                              TrackCuts.pidnSigmaPreSelectionMeanKaon.value,
                              TrackCuts.pidnSigmaPreSelectionMeanTOFKaon.value, tpcCut, tofCut)) {
      return true;
    }
    if (FilterForDerivedTables.cfgFillProtonTracks &&
        passesPIDPreSelection(track.tpcNSigmaPr(), track.tofNSigmaPr(), track.hasTOF(),
                              TrackCuts.pidnSigmaPreSelectionMeanProton.value,
                              TrackCuts.pidnSigmaPreSelectionMeanTOFProton.value, tpcCut, tofCut)) {
      return true;
    }
    return false;
  }

  template <typename CollisionType, typename V0Type, typename TrackType>
  bool filterV0(CollisionType const& collision, V0Type const& v0, TrackType const&)
  {
    if (!FilterForDerivedTables.cfgFillK0s && !FilterForDerivedTables.cfgFillLambda0) {
      return true;
    }
    if (!SecondaryCuts.cfgSecondaryRequire) {
      return true;
    }
    if (v0.dcaV0daughters() > SecondaryCuts.cfgSecondaryDauDCAMax ||
        v0.pt() < SecondaryCuts.cfgSecondaryPtMin ||
        v0.dcav0topv() > SecondaryCuts.cfgSecondaryDCAtoPVMax) {
      return false;
    }
    if (SecondaryCuts.cfgSecondaryArmenterosCut &&
        v0.qtarm() < SecondaryCuts.cfgSecondaryparamArmenterosCut * std::abs(v0.alpha())) {
      return false;
    }

    const float decayLengthOverMomentum = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ());
    const auto posTrack = v0.template posTrack_as<TrackType>();
    const auto negTrack = v0.template negTrack_as<TrackType>();
    const bool bypassDaughterPID = SecondaryCuts.cfgByPassDauPIDSelection;
    const float daughterPIDCut = TrackCuts.pidnSigmaPreSelectionCut.value;
    bool selected = false;
    if (FilterForDerivedTables.cfgFillK0s) {
      const bool passesK0DaughterPID = bypassDaughterPID ||
                                       (passesCenteredPID(posTrack.tpcNSigmaPi(), TrackCuts.pidnSigmaPreSelectionMeanPion.value, daughterPIDCut) &&
                                        passesCenteredPID(negTrack.tpcNSigmaPi(), TrackCuts.pidnSigmaPreSelectionMeanPion.value, daughterPIDCut));
      const bool passesK0 = std::fabs(v0.yK0Short()) <= SecondaryCuts.cfgSecondaryRapidityMax &&
                            decayLengthOverMomentum * MassK0Short <= SecondaryCuts.cfgSecondaryProperLifetimeMax &&
                            std::fabs(v0.mK0Short() - MassK0Short) <= SecondaryCuts.cfgSecondaryMassWindow &&
                            passesK0DaughterPID &&
                            (!SecondaryCuts.cfgSecondaryCrossMassHypothesisCut ||
                             (std::fabs(v0.mLambda() - MassLambda0) >= SecondaryCuts.cfgSecondaryCrossMassCutWindow &&
                              std::fabs(v0.mAntiLambda() - MassLambda0Bar) >= SecondaryCuts.cfgSecondaryCrossMassCutWindow));
      selected = passesK0;
    }
    if (FilterForDerivedTables.cfgFillLambda0) {
      const bool passesLambdaPID = bypassDaughterPID ||
                                   (passesCenteredPID(posTrack.tpcNSigmaPr(), TrackCuts.pidnSigmaPreSelectionMeanProton.value, daughterPIDCut) &&
                                    passesCenteredPID(negTrack.tpcNSigmaPi(), TrackCuts.pidnSigmaPreSelectionMeanPion.value, daughterPIDCut));
      const bool passesAntiLambdaPID = bypassDaughterPID ||
                                       (passesCenteredPID(posTrack.tpcNSigmaPi(), TrackCuts.pidnSigmaPreSelectionMeanPion.value, daughterPIDCut) &&
                                        passesCenteredPID(negTrack.tpcNSigmaPr(), TrackCuts.pidnSigmaPreSelectionMeanProton.value, daughterPIDCut));
      const bool passesLambdaMassAndPID =
        (std::fabs(v0.mLambda() - MassLambda0) <= SecondaryCuts.cfgSecondaryMassWindow && passesLambdaPID) ||
        (std::fabs(v0.mAntiLambda() - MassLambda0Bar) <= SecondaryCuts.cfgSecondaryMassWindow && passesAntiLambdaPID);
      const bool passesLambda = std::fabs(v0.yLambda()) <= SecondaryCuts.cfgSecondaryRapidityMax &&
                                decayLengthOverMomentum * MassLambda0 <= SecondaryCuts.cfgSecondaryProperLifetimeMax &&
                                passesLambdaMassAndPID &&
                                (!SecondaryCuts.cfgSecondaryCrossMassHypothesisCut ||
                                 std::fabs(v0.mK0Short() - MassK0Short) >= SecondaryCuts.cfgSecondaryCrossMassCutWindow);
      selected = selected || passesLambda;
    }
    return selected;
  }

  template <bool isMC, typename CollisionType, typename TrackType>
  bool isMicroTrackSelected(CollisionType const&, TrackType const& track)
  {
    if (!std::isfinite(track.dcaXY()) || !std::isfinite(track.dcaZ())) {
      return false;
    }
    if (std::fabs(track.dcaXY()) > TrackCuts.cMaxDCArToPVcut) {
      return false;
    }
    if (std::fabs(track.dcaZ()) > TrackCuts.cMaxDCAzToPVcut || std::fabs(track.dcaZ()) < TrackCuts.cMinDCAzToPVcut) {
      return false;
    }

    return true;
  }

  template <bool isMC, typename CollisionType, typename TrackType>
  bool isTrackSelected(CollisionType const&, TrackType const& track, bool fillSelectionQA = true)
  {
    if (cfgFillQA && fillSelectionQA) {
      qaRegistry.fill(HIST("QA/hGoodTrackIndices"), 0.5);
      if constexpr (isMC) {
        qaRegistry.fill(HIST("QA/hGoodMCTrackIndices"), 0.5);
      }
    }
    if (!std::isfinite(track.dcaXY()) || !std::isfinite(track.dcaZ())) {
      return false;
    }
    if (std::fabs(track.dcaXY()) > TrackCuts.cMaxDCArToPVcut) {
      return false;
    }
    if (cfgFillQA && fillSelectionQA) {
      qaRegistry.fill(HIST("QA/hGoodTrackIndices"), 1.5);
    }
    if (std::fabs(track.dcaZ()) > TrackCuts.cMaxDCAzToPVcut || std::fabs(track.dcaZ()) < TrackCuts.cMinDCAzToPVcut) {
      return false;
    }
    if (cfgFillQA && fillSelectionQA) {
      qaRegistry.fill(HIST("QA/hGoodTrackIndices"), 2.5);
      qaRegistry.fill(HIST("QA/hGoodTrackIndices"), 7.5);
    }
    return true;
  }

  template <bool isMC, typename CollisionType, typename V0Type, typename TrackType>
  bool isV0Selected(CollisionType const&, V0Type const& v0, TrackType const&, bool fillSelectionQA = true)
  {
    if (cfgFillQA && fillSelectionQA) {
      qaRegistry.fill(HIST("QA/hGoodV0Indices"), 0.5);
    }

    auto posTrack = v0.template posTrack_as<TrackType>();
    auto negTrack = v0.template negTrack_as<TrackType>();
    if (posTrack.tpcNClsCrossedRows() < V0Cuts.mincrossedrowsV0s || negTrack.tpcNClsCrossedRows() < V0Cuts.mincrossedrowsV0s) {
      return false;
    }
    if (cfgFillQA && fillSelectionQA) {
      qaRegistry.fill(HIST("QA/hGoodV0Indices"), 1.5);
    }
    if (std::fabs(posTrack.dcaXY()) < V0Cuts.cMinV0PosDCArToPVcut ||
        std::fabs(negTrack.dcaXY()) < V0Cuts.cMinV0NegDCArToPVcut) {
      return false;
    }
    if (cfgFillQA && fillSelectionQA) {
      qaRegistry.fill(HIST("QA/hGoodV0Indices"), 2.5);
    }
    if (v0.v0radius() > V0Cuts.cMaxV0Radius || v0.v0radius() < V0Cuts.cMinV0Radius) {
      return false;
    }
    if (cfgFillQA && fillSelectionQA) {
      qaRegistry.fill(HIST("QA/hGoodV0Indices"), 3.5);
    }
    if (v0.v0cosPA() < V0Cuts.cMinV0CosPA) {
      return false;
    }
    if (cfgFillQA && fillSelectionQA) {
      qaRegistry.fill(HIST("QA/hGoodV0Indices"), 4.5);
      if constexpr (isMC) {
        qaRegistry.fill(HIST("QA/hGoodMCV0Indices"), 0.5);
      }
    }
    return true;
  }

  template <bool isMC, typename CollisionType, typename CascType, typename TrackType>
  bool isCascSelected(CollisionType const& collision, CascType const& casc, TrackType const&, bool fillSelectionQA = true)
  {
    if (cfgFillQA && fillSelectionQA) {
      qaRegistry.fill(HIST("QA/hGoodCascIndices"), 0.5);
    }

    auto bachelor = casc.template bachelor_as<TrackType>();
    if (bachelor.tpcNClsCrossedRows() < CascadeCuts.cfgMinCrossedRowsCascBach) {
      return false;
    }
    if (cfgFillQA && fillSelectionQA) {
      qaRegistry.fill(HIST("QA/hGoodCascIndices"), 1.5);
    }
    if (std::fabs(bachelor.dcaXY()) < CascadeCuts.cMinCascBachDCArToPVcut ||
        std::fabs(bachelor.dcaXY()) > CascadeCuts.cMaxCascBachDCArToPVcut) {
      return false;
    }
    if (cfgFillQA && fillSelectionQA) {
      qaRegistry.fill(HIST("QA/hGoodCascIndices"), 2.5);
    }
    if (casc.dcaV0daughters() > CascadeCuts.cMaxCascDCAV0Daughters ||
        casc.dcacascdaughters() > CascadeCuts.cMaxCascDCACascDaughters) {
      return false;
    }
    if (cfgFillQA && fillSelectionQA) {
      qaRegistry.fill(HIST("QA/hGoodCascIndices"), 3.5);
    }
    if (casc.casccosPA(collision.posX(), collision.posY(), collision.posZ()) < CascadeCuts.cMinCascCosPA ||
        casc.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) < CascadeCuts.cMinCascV0CosPA) {
      return false;
    }
    if (cfgFillQA && fillSelectionQA) {
      qaRegistry.fill(HIST("QA/hGoodCascIndices"), 4.5);
    }
    auto v0Radius = casc.v0radius();
    if (v0Radius > CascadeCuts.cMaxCascV0Radius || v0Radius < CascadeCuts.cMinCascV0Radius) {
      return false;
    }
    if (cfgFillQA && fillSelectionQA) {
      qaRegistry.fill(HIST("QA/hGoodCascIndices"), 5.5);
    }
    auto cascRadius = casc.cascradius();
    if (cascRadius > CascadeCuts.cMaxCascRadius || cascRadius < CascadeCuts.cMinCascRadius) {
      return false;
    }
    if (cfgFillQA && fillSelectionQA) {
      qaRegistry.fill(HIST("QA/hGoodCascIndices"), 6.5);
    }
    if (std::abs(casc.mXi() - MassXiMinus) > CascadeCuts.cMaxXiMassWindow) {
      return false;
    }
    if (cfgFillQA && fillSelectionQA) {
      qaRegistry.fill(HIST("QA/hGoodCascIndices"), 7.5);
      if constexpr (isMC) {
        qaRegistry.fill(HIST("QA/hGoodMCCascIndices"), 0.5);
      }
    }
    return true;
  }

  /// @brief Find a selected V0 and collect daughter IDs only for the optional global veto
  template <bool isMC, typename CollisionType, typename V0Type, typename TrackType>
  SelectedCandidateDaughters collectSelectedV0Daughters(CollisionType const& collision,
                                                        V0Type const& v0s,
                                                        TrackType const& tracks,
                                                        bool useGlobalDaughterVeto)
  {
    SelectedCandidateDaughters selectedCandidates{useGlobalDaughterVeto};
    for (auto const& v0 : v0s) {
      if (!isV0Selected<isMC>(collision, v0, tracks, false) ||
          !filterV0(collision, v0, tracks)) {
        continue;
      }
      if (!useGlobalDaughterVeto) {
        selectedCandidates.addCandidate();
        break;
      }
      selectedCandidates.addCandidate(
        std::array<int64_t, 2>{static_cast<int64_t>(v0.posTrackId()),
                               static_cast<int64_t>(v0.negTrackId())});
    }
    selectedCandidates.finalize();
    return selectedCandidates;
  }

  /// @brief Find a selected cascade and collect daughter IDs only for the optional global veto
  template <bool isMC, typename CollisionType, typename CascType, typename TrackType>
  SelectedCandidateDaughters collectSelectedCascadeDaughters(CollisionType const& collision,
                                                             CascType const& cascades,
                                                             TrackType const& tracks,
                                                             bool useGlobalDaughterVeto)
  {
    SelectedCandidateDaughters selectedCandidates{useGlobalDaughterVeto};
    for (auto const& casc : cascades) {
      if (!isCascSelected<isMC>(collision, casc, tracks, false)) {
        continue;
      }
      if (!useGlobalDaughterVeto) {
        selectedCandidates.addCandidate();
        break;
      }
      selectedCandidates.addCandidate(
        std::array<int64_t, 3>{static_cast<int64_t>(casc.posTrackId()),
                               static_cast<int64_t>(casc.negTrackId()),
                               static_cast<int64_t>(casc.bachelorId())});
    }
    selectedCandidates.finalize();
    return selectedCandidates;
  }

  /// @brief Build the per-track selection while preserving the configured collision gate
  template <bool isMC,
            typename CollisionType,
            typename TrackType,
            typename V0Type,
            typename CascType,
            typename V0PresliceType,
            typename CascPresliceType>
  bool preparePairTrackSelection(CollisionType const& collision,
                                 TrackType const& tracks,
                                 V0Type const& v0s,
                                 CascType const& cascades,
                                 V0PresliceType const& v0Preslice,
                                 CascPresliceType const& cascPreslice,
                                 PairTrackSelection& selection)
  {
    const bool useEitherPairGate = FilterForDerivedTables.cfgPairGateMode.value == PairGateModeEither;
    selection.useV0Candidates = useEitherPairGate || FilterForDerivedTables.cfgBypassNoPairV0s;
    selection.useCascadeCandidates = useEitherPairGate || FilterForDerivedTables.cfgBypassNoPairCascades;

    if (selection.useV0Candidates) {
      auto v0sThisCollision = v0s.sliceBy(v0Preslice, collision.collisionId());
      selection.v0Candidates =
        collectSelectedV0Daughters<isMC>(collision, v0sThisCollision, tracks, selection.useGlobalDaughterVeto);
    }
    if (selection.useCascadeCandidates &&
        (!useEitherPairGate || selection.useGlobalDaughterVeto || !selection.v0Candidates.hasSelectedCandidate)) {
      auto cascadesThisCollision = cascades.sliceBy(cascPreslice, collision.collisionId());
      selection.cascadeCandidates =
        collectSelectedCascadeDaughters<isMC>(collision, cascadesThisCollision, tracks, selection.useGlobalDaughterVeto);
    }

    if (useEitherPairGate) {
      return selection.v0Candidates.hasSelectedCandidate ||
             selection.cascadeCandidates.hasSelectedCandidate;
    }

    if (FilterForDerivedTables.cfgBypassNoPairV0s && !selection.v0Candidates.hasSelectedCandidate) {
      return false;
    }
    if (FilterForDerivedTables.cfgBypassNoPairCascades && !selection.cascadeCandidates.hasSelectedCandidate) {
      return false;
    }
    return true;
  }

  /// @brief Round and saturate a floating-point value into a persistent integer column
  template <typename Integer>
  static Integer quantizeSaturated(float value, double scale)
  {
    constexpr Integer LowerLimit = std::numeric_limits<Integer>::lowest();
    constexpr Integer UpperLimit = std::numeric_limits<Integer>::max();
    if (std::isnan(value)) {
      return UpperLimit;
    }
    if (!std::isfinite(value)) {
      return std::signbit(value) ? LowerLimit : UpperLimit;
    }
    const double rounded = std::round(static_cast<double>(value) * scale);
    if (rounded <= static_cast<double>(LowerLimit)) {
      return LowerLimit;
    }
    if (rounded >= static_cast<double>(UpperLimit)) {
      return UpperLimit;
    }
    return static_cast<Integer>(rounded);
  }

  static bool quantizeP(float p, int16_t& quantized)
  {
    if (!std::isfinite(p)) {
      return false;
    }
    const double rounded = std::round(static_cast<double>(p) * MomentumQuantizationScale);
    if (rounded < static_cast<double>(std::numeric_limits<int16_t>::min()) ||
        rounded > static_cast<double>(std::numeric_limits<int16_t>::max())) {
      return false;
    }
    quantized = static_cast<int16_t>(rounded);
    return true;
  }

  float tightDCAThreshold(float pt) const
  {
    if (!std::isfinite(pt) || pt <= 0.f) {
      return -1.f;
    }
    if (TrackCuts.cfgTightDCAPtCoefficient.value == 0.f) {
      return TrackCuts.cfgTightDCAOffset.value;
    }
    const float threshold = TrackCuts.cfgTightDCAOffset.value +
                            TrackCuts.cfgTightDCAPtCoefficient.value /
                              std::pow(pt, TrackCuts.cfgTightDCAPtPower.value);
    return std::isfinite(threshold) ? threshold : -1.f;
  }

  template <typename TrackType>
  bool evaluatePtDependentDCA(TrackType const& track,
                              bool& passedPtDependentDCAxy,
                              bool& passedPtDependentDCAz) const
  {
    passedPtDependentDCAxy = false;
    passedPtDependentDCAz = false;
    if (!TrackCuts.cfgApplyTightDCAPtDepSelection.value) {
      return true;
    }
    const float dcaThreshold = tightDCAThreshold(track.pt());
    if (dcaThreshold < 0.f) {
      return false;
    }
    passedPtDependentDCAxy = std::isfinite(track.dcaXY()) && std::abs(track.dcaXY()) < dcaThreshold;
    passedPtDependentDCAz = std::isfinite(track.dcaZ()) && std::abs(track.dcaZ()) < dcaThreshold;
    return passedPtDependentDCAxy && passedPtDependentDCAz;
  }

  template <bool isMC, typename CollisionType, typename TrackType>
  bool isFullTrackOutputSelected(CollisionType const& collision,
                                 TrackType const& track,
                                 bool fillSelectionQA = true)
  {
    return isTrackSelected<isMC>(collision, track, fillSelectionQA) && filterTrack(track);
  }

  template <bool isMC, typename CollisionType, typename TrackType>
  bool isMicroTrackOutputSelected(CollisionType const& collision,
                                  TrackType const& track,
                                  bool& passedPtDependentDCAxy,
                                  bool& passedPtDependentDCAz)
  {
    return isMicroTrackSelected<isMC>(collision, track) &&
           filterTrack(track) &&
           evaluatePtDependentDCA(track, passedPtDependentDCAxy, passedPtDependentDCAz);
  }

  template <bool isMC, typename CollisionType, typename TrackType>
  bool isUltraMicroTrackOutputSelected(CollisionType const& collision, TrackType const& track)
  {
    return isMicroTrackSelected<isMC>(collision, track) &&
           filterTrack(track) &&
           o2::aod::resoultramicrodaughter::DCAEncoding::isValid(track.dcaXY()) &&
           o2::aod::resoultramicrodaughter::DCAEncoding::isValid(track.dcaZ());
  }

  template <typename TrackType>
  static bool quantizeUltraMicroMomentum(TrackType const& track,
                                         int16_t& px1000,
                                         int16_t& py1000,
                                         int16_t& pz1000)
  {
    return quantizeP(track.px(), px1000) &&
           quantizeP(track.py(), py1000) &&
           quantizeP(track.pz(), pz1000);
  }

  /// Check that every enabled track output will receive at least one row.
  template <bool isMC, typename CollisionType, typename TrackTableType, typename TrackPredicate>
  bool hasTracksForEnabledOutputs(CollisionType const& collision,
                                  TrackTableType const& tracks,
                                  TrackPredicate const& keepTrack)
  {
    constexpr uint8_t FullTrackOutput = 1u << 0;
    constexpr uint8_t MicroTrackOutput = 1u << 1;
    constexpr uint8_t UltraMicroTrackOutput = 1u << 2;
    uint8_t requiredOutputs = 0;
    if (!FilterForDerivedTables.cfgBypassTrackFill.value) {
      requiredOutputs |= FullTrackOutput;
    }
    if (FilterForDerivedTables.cfgFillMicroTracks.value) {
      requiredOutputs |= MicroTrackOutput;
    }
    if (FilterForDerivedTables.cfgFillUltraMicroTracks.value) {
      requiredOutputs |= UltraMicroTrackOutput;
    }

    uint8_t availableOutputs = 0;
    for (auto const& track : tracks) {
      if (!keepTrack(track)) {
        continue;
      }
      if ((requiredOutputs & FullTrackOutput) != 0 &&
          (availableOutputs & FullTrackOutput) == 0 &&
          isFullTrackOutputSelected<isMC>(collision, track, false)) {
        availableOutputs |= FullTrackOutput;
      }
      if ((requiredOutputs & MicroTrackOutput) != 0 &&
          (availableOutputs & MicroTrackOutput) == 0) {
        bool passedPtDependentDCAxy = false;
        bool passedPtDependentDCAz = false;
        if (isMicroTrackOutputSelected<isMC>(collision, track,
                                             passedPtDependentDCAxy,
                                             passedPtDependentDCAz)) {
          availableOutputs |= MicroTrackOutput;
        }
      }
      if ((requiredOutputs & UltraMicroTrackOutput) != 0 &&
          (availableOutputs & UltraMicroTrackOutput) == 0 &&
          isUltraMicroTrackOutputSelected<isMC>(collision, track)) {
        int16_t px1000 = 0;
        int16_t py1000 = 0;
        int16_t pz1000 = 0;
        if (quantizeUltraMicroMomentum(track, px1000, py1000, pz1000)) {
          availableOutputs |= UltraMicroTrackOutput;
        }
      }
      if (availableOutputs == requiredOutputs) {
        return true;
      }
    }
    return false;
  }

  template <bool isMC, typename TrackType, typename CollisionType, typename TrackPredicate>
  void fillUltraMicroTracks(CollisionType const& collision, TrackType const& tracks, TrackPredicate const& keepTrack)
  {
    // Loop over tracks
    for (auto const& track : tracks) {
      if (!keepTrack(track)) {
        continue;
      }
      if (!isUltraMicroTrackOutputSelected<isMC>(collision, track)) {
        continue;
      }
      o2::aod::resoultramicrodaughter::DCAEncoding dcaEncoding(track.dcaXY(), track.dcaZ());
      int16_t px1000 = 0;
      int16_t py1000 = 0;
      int16_t pz1000 = 0;
      if (!quantizeUltraMicroMomentum(track, px1000, py1000, pz1000)) {
        if (!warnedUltraMicroMomentumRange) {
          LOGF(warn, "Skipping ultra-micro tracks with non-finite or out-of-range momentum components");
          warnedUltraMicroMomentumRange = true;
        }
        continue;
      }
      uint8_t trackFlags = (track.passedITSRefit() << 0) |
                           (track.passedTPCRefit() << 1) |
                           (track.isGlobalTrackWoDCA() << 2) |
                           (track.isGlobalTrack() << 3) |
                           (track.isPrimaryTrack() << 4) |
                           (track.isPVContributor() << 5) |
                           (track.hasTOF() << 6) |
                           ((track.sign() > 0) << 7); // sign +1: 1, -1: 0
      uint8_t pidFlag = 0;
      switch (ultraMicroPidSpecies) {
        case UltraMicroPidSpecies::Pion:
          pidFlag = static_cast<uint8_t>(o2::aod::resoultramicrodaughter::PidNSigma(
            track.tpcNSigmaPi(), track.tofNSigmaPi(), track.hasTOF()));
          break;
        case UltraMicroPidSpecies::Kaon:
          pidFlag = static_cast<uint8_t>(o2::aod::resoultramicrodaughter::PidNSigma(
            track.tpcNSigmaKa(), track.tofNSigmaKa(), track.hasTOF()));
          break;
        case UltraMicroPidSpecies::Proton:
          pidFlag = static_cast<uint8_t>(o2::aod::resoultramicrodaughter::PidNSigma(
            track.tpcNSigmaPr(), track.tofNSigmaPr(), track.hasTOF()));
          break;
      }
      if (cfgFillQA && cfgDetailTrackQA) {
        qaRegistry.fill(HIST("QA/h4UltraMicroTrackPtEtaPhi"), track.pt(), track.eta(), track.phi());
        qaRegistry.fill(HIST("QA/h2UltraMicroTrackDCAxyVsPt"), track.pt(), track.dcaXY());
        qaRegistry.fill(HIST("QA/h2UltraMicroTrackDCAzVsPt"), track.pt(), track.dcaZ());
        qaRegistry.fill(HIST("QA/h4UltraMicroTrackTPCnSigma"), track.pt(), track.tpcNSigmaPi(), track.tpcNSigmaKa(), track.tpcNSigmaPr());
        if (track.hasTOF()) {
          qaRegistry.fill(HIST("QA/h4UltraMicroTrackTOFnSigma"), track.pt(), track.tofNSigmaPi(), track.tofNSigmaKa(), track.tofNSigmaPr());
        }
      }
      reso2ultramicrotrks(collision.globalIndex(),
                          px1000,
                          py1000,
                          pz1000,
                          pidFlag,
                          static_cast<uint8_t>(dcaEncoding),
                          trackFlags);
      if (!FilterForDerivedTables.cfgBypassTrackIndexFill) {
        resoUltraMicroTrackTracks(track.globalIndex());
      }
    }
  }

  /**
   * @brief Fills track data
   *
   * @tparam isMC Boolean indicating if it's MC
   * @tparam TrackType Type of track
   * @tparam CollisionType Type of collision
   * @param collision Collision data
   * @param tracks Track data
   * @note ResoMicroTracks_001 is intended for fixed producer-side selections.
   *       Exact downstream cuts on its decoded values require zero-centred
   *       strict PID windows |nSigma| < C with
   *       C in {2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5}, and strict DCA windows
   *       |DCA| < C with C = N * 0.025 cm (N = 1, ..., 6). The listed PID
   *       guarantee assumes a zero mean; shifted PID windows require both
   *       interval edges to align with encoding boundaries and separate
   *       validation. PID cuts below 2 sigma and off-grid DCA/PID cuts cannot
   *       be reconstructed exactly. In particular, pT-dependent DCA cuts must
   *       be applied to the unquantised values here via
   *       cfgApplyTightDCAPtDepSelection; consumers should retain that producer
   *       decision rather than retune the decoded lower-edge DCA values.
   */
  template <bool isMC, typename TrackType, typename CollisionType, typename TrackPredicate>
  void fillMicroTracks(CollisionType const& collision, TrackType const& tracks, TrackPredicate const& keepTrack)
  {
    // Loop over tracks
    for (auto const& track : tracks) {
      if (!keepTrack(track)) {
        continue;
      }
      bool passedPtDependentDCAxy = false;
      bool passedPtDependentDCAz = false;
      if (!isMicroTrackOutputSelected<isMC>(collision, track,
                                            passedPtDependentDCAxy,
                                            passedPtDependentDCAz)) {
        continue;
      }
      const o2::aod::resomicrodaughter001::DCAEncoding trackSelFlag(
        track.dcaXY(), track.dcaZ(), passedPtDependentDCAxy, passedPtDependentDCAz);
      uint8_t trackFlags = (track.passedITSRefit() << 0) |
                           (track.passedTPCRefit() << 1) |
                           (track.isGlobalTrackWoDCA() << 2) |
                           (track.isGlobalTrack() << 3) |
                           (track.isPrimaryTrack() << 4) |
                           (track.isPVContributor() << 5) |
                           (track.hasTOF() << 6) |
                           ((track.sign() > 0) << 7); // sign +1: 1, -1: 0
      if (cfgFillQA && cfgDetailTrackQA) {
        qaRegistry.fill(HIST("QA/h4MicroTrackPtEtaPhi"), track.pt(), track.eta(), track.phi());
        qaRegistry.fill(HIST("QA/h2MicroTrackDCAxyVsPt"), track.pt(), track.dcaXY());
        qaRegistry.fill(HIST("QA/h2MicroTrackDCAzVsPt"), track.pt(), track.dcaZ());
        qaRegistry.fill(HIST("QA/h4MicroTrackTPCnSigma"), track.pt(), track.tpcNSigmaPi(), track.tpcNSigmaKa(), track.tpcNSigmaPr());
        if (track.hasTOF()) {
          qaRegistry.fill(HIST("QA/h4MicroTrackTOFnSigma"), track.pt(), track.tofNSigmaPi(), track.tofNSigmaKa(), track.tofNSigmaPr());
        }
      }
      reso2microtrks(collision.globalIndex(),
                     track.globalIndex(),
                     track.px(),
                     track.py(),
                     track.pz(),
                     static_cast<uint8_t>(o2::aod::resomicrodaughter001::PidNSigma(track.tpcNSigmaPi(), track.tofNSigmaPi(), track.hasTOF())),
                     static_cast<uint8_t>(o2::aod::resomicrodaughter001::PidNSigma(track.tpcNSigmaKa(), track.tofNSigmaKa(), track.hasTOF())),
                     static_cast<uint8_t>(o2::aod::resomicrodaughter001::PidNSigma(track.tpcNSigmaPr(), track.tofNSigmaPr(), track.hasTOF())),
                     static_cast<uint8_t>(trackSelFlag),
                     trackFlags);
      if (!FilterForDerivedTables.cfgBypassTrackIndexFill) {
        resoMicroTrackTracks(track.globalIndex());
      }
      if constexpr (isMC) {
        fillMCTrack(track, reso2mcmicrotrks);
      }
    }
  }

  /**
   * @brief Fills track data
   *
   * @tparam isMC Boolean indicating if it's MC
   * @tparam TrackType Type of track
   * @tparam CollisionType Type of collision
   * @param collision Collision data
   * @param tracks Track data
   */
  template <bool isMC, typename TrackType, typename CollisionType, typename TrackPredicate>
  void fillTracks(CollisionType const& collision, TrackType const& tracks, TrackPredicate const& keepTrack)
  {
    if (FilterForDerivedTables.cfgBypassTrackFill) {
      return;
    }
    // Loop over tracks
    for (auto const& track : tracks) {
      if (!keepTrack(track)) {
        continue;
      }
      if (!isFullTrackOutputSelected<isMC>(collision, track)) {
        continue;
      }
      uint8_t trackFlags = (track.passedITSRefit() << 0) |
                           (track.passedTPCRefit() << 1) |
                           (track.isGlobalTrackWoDCA() << 2) |
                           (track.isGlobalTrack() << 3) |
                           (track.isPrimaryTrack() << 4) |
                           (track.isPVContributor() << 5) |
                           (track.hasTOF() << 6) |
                           ((track.sign() > 0) << 7); // sign +1: 1, -1: 0
      if (cfgFillQA && cfgDetailTrackQA) {
        qaRegistry.fill(HIST("QA/h4TrackPtEtaPhi"), track.pt(), track.eta(), track.phi());
        qaRegistry.fill(HIST("QA/h2TrackDCAxyVsPt"), track.pt(), track.dcaXY());
        qaRegistry.fill(HIST("QA/h2TrackDCAzVsPt"), track.pt(), track.dcaZ());
        qaRegistry.fill(HIST("QA/h4TrackTPCnSigma"), track.pt(), track.tpcNSigmaPi(), track.tpcNSigmaKa(), track.tpcNSigmaPr());
        if (track.hasTOF()) {
          qaRegistry.fill(HIST("QA/h4TrackTOFnSigma"), track.pt(), track.tofNSigmaPi(), track.tofNSigmaKa(), track.tofNSigmaPr());
        }
      }
      reso2trks(collision.globalIndex(),
                track.pt(),
                track.px(),
                track.py(),
                track.pz(),
                static_cast<uint8_t>(track.tpcNClsCrossedRows()),
                static_cast<uint8_t>(track.tpcNClsFound()),
                quantizeSaturated<int16_t>(track.dcaXY(), 10000.),
                quantizeSaturated<int16_t>(track.dcaZ(), 10000.),
                quantizeSaturated<int8_t>(track.tpcNSigmaPi(), 10.),
                quantizeSaturated<int8_t>(track.tpcNSigmaKa(), 10.),
                quantizeSaturated<int8_t>(track.tpcNSigmaPr(), 10.),
                quantizeSaturated<int8_t>(track.tofNSigmaPi(), 10.),
                quantizeSaturated<int8_t>(track.tofNSigmaKa(), 10.),
                quantizeSaturated<int8_t>(track.tofNSigmaPr(), 10.),
                quantizeSaturated<int16_t>(track.tpcSignal(), 100.),
                trackFlags);
      if (!FilterForDerivedTables.cfgBypassTrackIndexFill) {
        resoTrackTracks(track.globalIndex());
      }
      if constexpr (isMC) {
        fillMCTrack(track, reso2mctracks);
      }
    }
  }

  /**
   * @brief Fills MC track data
   *
   * @tparam TrackType Type of track
   * @tparam MCOutputTable Type of positional MC extension to fill
   * @param track Track data
   * @param output Positional MC extension writer
   */
  template <typename TrackType, typename MCOutputTable>
  void fillMCTrack(TrackType const& track, Produces<MCOutputTable>& output)
  {
    // ------ Temporal lambda function to prevent error in build
    auto getMothersIndeces = [&](auto const& theMcParticle) {
      std::vector<int> lMothersIndeces{};
      for (auto const& lMother : theMcParticle.template mothers_as<aod::McParticles>()) {
        LOGF(debug, "   mother index lMother: %d", lMother.globalIndex());
        lMothersIndeces.push_back(lMother.globalIndex());
      }
      return lMothersIndeces;
    };
    auto getMothersPDGCodes = [&](auto const& theMcParticle) {
      std::vector<int> lMothersPDGs{};
      for (auto const& lMother : theMcParticle.template mothers_as<aod::McParticles>()) {
        LOGF(debug, "   mother pdgcode lMother: %d", lMother.pdgCode());
        lMothersPDGs.push_back(lMother.pdgCode());
      }
      return lMothersPDGs;
    };
    auto getSiblingsIndeces = [&](auto const& theMcParticle) {
      std::vector<int> lSiblingsIndeces{};
      for (auto const& lMother : theMcParticle.template mothers_as<aod::McParticles>()) {
        LOGF(debug, "   mother index lMother: %d", lMother.globalIndex());
        for (auto const& lDaughter : lMother.template daughters_as<aod::McParticles>()) {
          LOGF(debug, "   daughter index lDaughter: %d", lDaughter.globalIndex());
          if (lDaughter.globalIndex() != theMcParticle.globalIndex()) {
            lSiblingsIndeces.push_back(lDaughter.globalIndex());
          }
        }
      }
      return lSiblingsIndeces;
    };
    // ------
    std::vector<int> mothers = {-1, -1};
    std::vector<int> motherPDGs = {-1, -1};
    std::array<int, StoredMCRelationCount> siblings{-1, -1};
    std::vector<int> siblingsTemp{};
    if (track.has_mcParticle()) {
      // Get the MC particle
      const auto& particle = track.mcParticle();
      if (particle.has_mothers()) {
        mothers = getMothersIndeces(particle);
        motherPDGs = getMothersPDGCodes(particle);
        siblingsTemp = getSiblingsIndeces(particle);
      }
      mothers.resize(StoredMCRelationCount, -1);
      motherPDGs.resize(StoredMCRelationCount, -1);
      if (!siblingsTemp.empty()) {
        siblings[0] = siblingsTemp[0];
      }
      if (siblingsTemp.size() > 1) {
        siblings[1] = siblingsTemp[1];
      }
      output(particle.pdgCode(),
             mothers[0],
             motherPDGs[0],
             siblings.data(),
             particle.isPhysicalPrimary(),
             particle.producedByGenerator());
    } else {
      // No MC particle associated
      output(0,
             mothers[0],
             motherPDGs[0],
             siblings.data(),
             0,
             0);
    }
  }

  /**
   * @brief Fills V0 data
   *
   * @tparam isMC Boolean indicating if it's MC
   * @tparam CollisionType Type of collision
   * @tparam V0Type Type of V0
   * @tparam TrackType Type of track
   * @param collision Collision data
   * @param v0s V0 data
   * @param tracks Track data
   */
  template <bool isMC, typename CollisionType, typename V0Type, typename TrackType>
  void fillV0s(CollisionType const& collision, V0Type const& v0s, TrackType const& tracks)
  {
    for (auto const& v0 : v0s) {
      if (!isV0Selected<isMC>(collision, v0, tracks)) {
        continue;
      }
      if (!filterV0(collision, v0, tracks)) {
        continue;
      }
      const auto posTrack = v0.template posTrack_as<TrackType>();
      const auto negTrack = v0.template negTrack_as<TrackType>();
      const std::array<int, 2> childIDs{v0.posTrackId(), v0.negTrackId()}; // Original track IDs for downstream pair-level shared-daughter rejection
      reso2v0s(collision.globalIndex(),
               v0.pt(),
               v0.px(),
               v0.py(),
               v0.pz(),
               childIDs.data(),
               quantizeSaturated<int8_t>(posTrack.tpcNSigmaPi(), 10.),
               quantizeSaturated<int8_t>(posTrack.tpcNSigmaKa(), 10.),
               quantizeSaturated<int8_t>(posTrack.tpcNSigmaPr(), 10.),
               quantizeSaturated<int8_t>(negTrack.tpcNSigmaPi(), 10.),
               quantizeSaturated<int8_t>(negTrack.tpcNSigmaKa(), 10.),
               quantizeSaturated<int8_t>(negTrack.tpcNSigmaPr(), 10.),
               quantizeSaturated<int8_t>(posTrack.tofNSigmaPi(), 10.),
               quantizeSaturated<int8_t>(posTrack.tofNSigmaKa(), 10.),
               quantizeSaturated<int8_t>(posTrack.tofNSigmaPr(), 10.),
               quantizeSaturated<int8_t>(negTrack.tofNSigmaPi(), 10.),
               quantizeSaturated<int8_t>(negTrack.tofNSigmaKa(), 10.),
               quantizeSaturated<int8_t>(negTrack.tofNSigmaPr(), 10.),
               v0.v0cosPA(),
               v0.dcaV0daughters(),
               v0.dcapostopv(),
               v0.dcanegtopv(),
               v0.dcav0topv(),
               static_cast<uint8_t>(posTrack.tpcNClsCrossedRows()),
               static_cast<uint8_t>(negTrack.tpcNClsCrossedRows()),
               v0.mLambda(),
               v0.mAntiLambda(),
               v0.mK0Short(),
               v0.v0radius(), v0.x(), v0.y(), v0.z(),
               v0.alpha(), v0.qtarm());
      if (!FilterForDerivedTables.cfgBypassTrackIndexFill) {
        resoV0V0s(v0.v0Id());
      }
      if constexpr (isMC) {
        fillMCV0(v0);
      }
    }
  }

  /**
   * @brief Fills MC V0 data
   *
   * @tparam V0Type Type of V0
   * @param v0 V0 data
   */
  template <typename V0Type>
  void fillMCV0(V0Type const& v0)
  {
    // ------ Temporal lambda function to prevent error in build
    auto getMothersIndeces = [&](auto const& theMcParticle) {
      std::vector<int> lMothersIndeces{};
      for (auto const& lMother : theMcParticle.template mothers_as<aod::McParticles>()) {
        LOGF(debug, "   mother index lMother: %d", lMother.globalIndex());
        lMothersIndeces.push_back(lMother.globalIndex());
      }
      return lMothersIndeces;
    };
    auto getMothersPDGCodes = [&](auto const& theMcParticle) {
      std::vector<int> lMothersPDGs{};
      for (auto const& lMother : theMcParticle.template mothers_as<aod::McParticles>()) {
        LOGF(debug, "   mother pdgcode lMother: %d", lMother.pdgCode());
        lMothersPDGs.push_back(lMother.pdgCode());
      }
      return lMothersPDGs;
    };
    auto getMothersPt = [&](auto const& theMcParticle) {
      std::vector<float> lMothersPts{};
      for (auto const& lMother : theMcParticle.template mothers_as<aod::McParticles>()) {
        LOGF(debug, "   mother pt lMother: %f", lMother.pt());
        lMothersPts.push_back(lMother.pt());
      }
      return lMothersPts;
    };
    auto getMothersRap = [&](auto const& theMcParticle) {
      std::vector<float> lMothersRaps{};
      for (auto const& lMother : theMcParticle.template mothers_as<aod::McParticles>()) {
        LOGF(debug, "   mother rap lMother: %f", lMother.y());
        lMothersRaps.push_back(lMother.y());
      }
      return lMothersRaps;
    };
    auto getDaughtersIndeces = [&](auto const& theMcParticle) {
      std::vector<int> lDaughtersIndeces{};
      for (auto const& lDaughter : theMcParticle.template daughters_as<aod::McParticles>()) {
        LOGF(debug, "   daughter index lDaughter: %d", lDaughter.globalIndex());
        lDaughtersIndeces.push_back(lDaughter.globalIndex());
      }
      return lDaughtersIndeces;
    };
    auto getDaughtersPDGCodes = [&](auto const& theMcParticle) {
      std::vector<int> lDaughtersPDGs{};
      for (auto const& lDaughter : theMcParticle.template daughters_as<aod::McParticles>()) {
        LOGF(debug, "   daughter pdgcode lDaughter: %d", lDaughter.pdgCode());
        lDaughtersPDGs.push_back(lDaughter.pdgCode());
      }
      return lDaughtersPDGs;
    };
    // ------
    std::vector<int> mothers = {-1, -1};
    std::vector<int> motherPDGs = {-1, -1};
    std::vector<float> mothersPts = {-1.0f, -1.0f};
    std::vector<float> mothersRaps = {-1.0f, -1.0f};
    std::vector<int> daughters = {-1, -1};
    std::vector<int> daughterPDGs = {-1, -1};
    if (v0.has_mcParticle()) {
      auto v0mc = v0.mcParticle();
      if (v0mc.has_mothers()) {
        mothers = getMothersIndeces(v0mc);
        motherPDGs = getMothersPDGCodes(v0mc);
        mothersPts = getMothersPt(v0mc);
        mothersRaps = getMothersRap(v0mc);
      }
      mothers.resize(StoredMCRelationCount, -1);
      motherPDGs.resize(StoredMCRelationCount, -1);
      mothersPts.resize(StoredMCRelationCount, -1.0f);
      mothersRaps.resize(StoredMCRelationCount, -1.0f);
      if (v0mc.has_daughters()) {
        daughters = getDaughtersIndeces(v0mc);
        daughterPDGs = getDaughtersPDGCodes(v0mc);
      }
      // if (daughters.size() > StoredMCRelationCount) {
      //   LOGF(info, "daughters.size() is larger than 2");
      // }
      daughters.resize(StoredMCRelationCount, -1);
      daughterPDGs.resize(StoredMCRelationCount, -1);
      reso2mcv0s(v0mc.pdgCode(),
                 mothers[0],
                 motherPDGs[0],
                 mothersPts[0],
                 mothersRaps[0],
                 daughters[0],
                 daughters[1],
                 daughterPDGs[0],
                 daughterPDGs[1],
                 v0mc.isPhysicalPrimary(),
                 v0mc.producedByGenerator());
    } else {
      reso2mcv0s(0,
                 mothers[0],
                 motherPDGs[0],
                 mothersPts[0],
                 mothersRaps[0],
                 daughters[0],
                 daughters[1],
                 daughterPDGs[0],
                 daughterPDGs[1],
                 0,
                 0);
    }
  }

  /**
   * @brief Fills cascade data
   *
   * @tparam isMC Boolean indicating if it's MC
   * @tparam CollisionType Type of collision
   * @tparam CascType Type of cascade
   * @tparam TrackType Type of track
   * @param collision Collision data
   * @param cascades Cascade data
   * @param tracks Track data
   */
  template <bool isMC, typename CollisionType, typename CascType, typename TrackType>
  void fillCascades(CollisionType const& collision, CascType const& cascades, TrackType const& tracks)
  {
    for (auto const& casc : cascades) {
      if (!isCascSelected<isMC>(collision, casc, tracks)) {
        continue;
      }
      const auto posTrack = casc.template posTrack_as<TrackType>();
      const auto negTrack = casc.template negTrack_as<TrackType>();
      const auto bachelor = casc.template bachelor_as<TrackType>();
      const std::array<int, 3> childIDs{casc.posTrackId(), casc.negTrackId(), casc.bachelorId()}; // Original track IDs for downstream pair-level shared-daughter rejection
      reso2cascades(collision.globalIndex(),
                    casc.pt(),
                    casc.px(),
                    casc.py(),
                    casc.pz(),
                    childIDs.data(),
                    quantizeSaturated<int8_t>(posTrack.tpcNSigmaPi(), 10.),
                    quantizeSaturated<int8_t>(posTrack.tpcNSigmaKa(), 10.),
                    quantizeSaturated<int8_t>(posTrack.tpcNSigmaPr(), 10.),
                    quantizeSaturated<int8_t>(negTrack.tpcNSigmaPi(), 10.),
                    quantizeSaturated<int8_t>(negTrack.tpcNSigmaKa(), 10.),
                    quantizeSaturated<int8_t>(negTrack.tpcNSigmaPr(), 10.),
                    quantizeSaturated<int8_t>(bachelor.tpcNSigmaPi(), 10.),
                    quantizeSaturated<int8_t>(bachelor.tpcNSigmaKa(), 10.),
                    quantizeSaturated<int8_t>(bachelor.tpcNSigmaPr(), 10.),
                    quantizeSaturated<int8_t>(posTrack.tofNSigmaPi(), 10.),
                    quantizeSaturated<int8_t>(posTrack.tofNSigmaKa(), 10.),
                    quantizeSaturated<int8_t>(posTrack.tofNSigmaPr(), 10.),
                    quantizeSaturated<int8_t>(negTrack.tofNSigmaPi(), 10.),
                    quantizeSaturated<int8_t>(negTrack.tofNSigmaKa(), 10.),
                    quantizeSaturated<int8_t>(negTrack.tofNSigmaPr(), 10.),
                    quantizeSaturated<int8_t>(bachelor.tofNSigmaPi(), 10.),
                    quantizeSaturated<int8_t>(bachelor.tofNSigmaKa(), 10.),
                    quantizeSaturated<int8_t>(bachelor.tofNSigmaPr(), 10.),
                    casc.v0cosPA(collision.posX(), collision.posY(), collision.posZ()),
                    casc.casccosPA(collision.posX(), collision.posY(), collision.posZ()),
                    casc.dcaV0daughters(),
                    casc.dcacascdaughters(),
                    casc.dcapostopv(),
                    casc.dcanegtopv(),
                    casc.dcabachtopv(),
                    casc.dcav0topv(collision.posX(), collision.posY(), collision.posZ()),
                    casc.dcaXYCascToPV(),
                    casc.dcaZCascToPV(),
                    casc.sign(),
                    static_cast<uint8_t>(posTrack.tpcNClsCrossedRows()),
                    static_cast<uint8_t>(negTrack.tpcNClsCrossedRows()),
                    static_cast<uint8_t>(bachelor.tpcNClsCrossedRows()),
                    casc.mLambda(),
                    casc.mXi(),
                    casc.v0radius(), casc.cascradius(), casc.x(), casc.y(), casc.z());
      if (!FilterForDerivedTables.cfgBypassTrackIndexFill) {
        resoCascadeCascades(casc.cascadeId());
      }
      if constexpr (isMC) {
        fillMCCascade(casc);
      }
    }
  }

  /**
   * @brief Fills MC cascade data
   *
   * @tparam CascType Type of cascade
   * @param casc Cascade data
   */
  template <typename CascType>
  void fillMCCascade(CascType const& casc)
  {
    // ------ Temporal lambda function to prevent error in build
    auto getMothersIndeces = [&](auto const& theMcParticle) {
      std::vector<int> lMothersIndeces{};
      for (auto const& lMother : theMcParticle.template mothers_as<aod::McParticles>()) {
        LOGF(debug, "   mother index lMother: %d", lMother.globalIndex());
        lMothersIndeces.push_back(lMother.globalIndex());
      }
      return lMothersIndeces;
    };
    auto getMothersPDGCodes = [&](auto const& theMcParticle) {
      std::vector<int> lMothersPDGs{};
      for (auto const& lMother : theMcParticle.template mothers_as<aod::McParticles>()) {
        LOGF(debug, "   mother pdgcode lMother: %d", lMother.pdgCode());
        lMothersPDGs.push_back(lMother.pdgCode());
      }
      return lMothersPDGs;
    };
    auto getMothersPt = [&](auto const& theMcParticle) {
      std::vector<float> lMothersPts{};
      for (auto const& lMother : theMcParticle.template mothers_as<aod::McParticles>()) {
        LOGF(debug, "   mother pt lMother: %f", lMother.pt());
        lMothersPts.push_back(lMother.pt());
      }
      return lMothersPts;
    };
    auto getMothersRap = [&](auto const& theMcParticle) {
      std::vector<float> lMothersRaps{};
      for (auto const& lMother : theMcParticle.template mothers_as<aod::McParticles>()) {
        LOGF(debug, "   mother rap lMother: %f", lMother.y());
        lMothersRaps.push_back(lMother.y());
      }
      return lMothersRaps;
    };
    auto getDaughtersIndeces = [&](auto const& theMcParticle) {
      std::vector<int> lDaughtersIndeces{};
      for (auto const& lDaughter : theMcParticle.template daughters_as<aod::McParticles>()) {
        LOGF(debug, "   daughter index lDaughter: %d", lDaughter.globalIndex());
        lDaughtersIndeces.push_back(lDaughter.globalIndex());
      }
      return lDaughtersIndeces;
    };
    auto getDaughtersPDGCodes = [&](auto const& theMcParticle) {
      std::vector<int> lDaughtersPDGs{};
      for (auto const& lDaughter : theMcParticle.template daughters_as<aod::McParticles>()) {
        LOGF(debug, "   daughter pdgcode lDaughter: %d", lDaughter.pdgCode());
        lDaughtersPDGs.push_back(lDaughter.pdgCode());
      }
      return lDaughtersPDGs;
    };
    // ------
    std::vector<int> mothers = {-1, -1};
    std::vector<int> motherPDGs = {-1, -1};
    std::vector<int> daughters = {-1, -1};
    std::vector<int> daughterPDGs = {-1, -1};
    std::vector<float> mothersPts = {-1.0f, -1.0f};
    std::vector<float> mothersRaps = {-1.0f, -1.0f};
    if (casc.has_mcParticle()) {
      auto cascmc = casc.mcParticle();
      if (cascmc.has_mothers()) {
        mothers = getMothersIndeces(cascmc);
        motherPDGs = getMothersPDGCodes(cascmc);
        mothersPts = getMothersPt(cascmc);
        mothersRaps = getMothersRap(cascmc);
      }
      mothers.resize(StoredMCRelationCount, -1);
      motherPDGs.resize(StoredMCRelationCount, -1);
      mothersPts.resize(StoredMCRelationCount, -1.0f);
      mothersRaps.resize(StoredMCRelationCount, -1.0f);
      if (cascmc.has_daughters()) {
        daughters = getDaughtersIndeces(cascmc);
        daughterPDGs = getDaughtersPDGCodes(cascmc);
      }
      // if (daughters.size() > StoredMCRelationCount) {
      //   LOGF(info, "daughters.size() is larger than 2");
      // }
      daughters.resize(StoredMCRelationCount, -1);
      daughterPDGs.resize(StoredMCRelationCount, -1);
      reso2mccascades(cascmc.pdgCode(),
                      mothers[0],
                      motherPDGs[0],
                      mothersPts[0],
                      mothersRaps[0],
                      daughters[0],
                      daughters[1],
                      daughterPDGs[0],
                      daughterPDGs[1],
                      cascmc.isPhysicalPrimary(),
                      cascmc.producedByGenerator());
    } else {
      reso2mccascades(0,
                      mothers[0],
                      motherPDGs[0],
                      mothersPts[0],
                      mothersRaps[0],
                      daughters[0],
                      daughters[1],
                      daughterPDGs[0],
                      daughterPDGs[1],
                      0,
                      0);
    }
  }

  /**
   * @brief Processes dummy
   *
   * @param collision Collision data
   */
  void processDummy(aod::ResoCollisions_001::iterator const&)
  {
  }
  PROCESS_SWITCH(ResonanceDaughterInitializer, processDummy, "Process dummy", true);

  /**
   * @brief Fills all enabled track tables from an already grouped track slice
   *
   * @tparam isMC Boolean indicating if it's MC
   * @param collision Reduced collision used as the output foreign key
   * @param tracks Tracks belonging to the corresponding original collision
   */
  template <bool isMC, typename CollisionType, typename TrackTableType, typename TrackPredicate>
  void fillTrackTables(CollisionType const& collision, TrackTableType const& tracks, TrackPredicate const& keepTrack)
  {
    fillTracks<isMC>(collision, tracks, keepTrack);
    if (FilterForDerivedTables.cfgFillMicroTracks) {
      fillMicroTracks<isMC>(collision, tracks, keepTrack);
    }
    if (FilterForDerivedTables.cfgFillUltraMicroTracks) {
      fillUltraMicroTracks<isMC>(collision, tracks, keepTrack);
    }
  }

  template <bool isMC, typename CollisionType, typename TrackTableType>
  void fillTrackTables(CollisionType const& collision, TrackTableType const& tracks)
  {
    fillTrackTables<isMC>(collision, tracks, KeepAllTracks{});
  }

  /**
   * @brief Fills track tables for one original collision
   *
   * @tparam isMC Boolean indicating if it's MC
   * @tparam CollisionType Type of reduced collision
   * @tparam TrackTableType Type of input track table
   * @tparam PresliceType Type of the original-collision preslice
   * @param collision Reduced collision with the original collision index
   * @param tracks Input track table
   */
  template <bool isMC, typename CollisionType, typename TrackTableType, typename PresliceType>
  void fillTrackTablesForCollision(CollisionType const& collision, TrackTableType const& tracks, PresliceType const& perCollision)
  {
    auto tracksThisCollision = tracks.sliceBy(perCollision, collision.collisionId());
    fillTrackTables<isMC>(collision, tracksThisCollision);
  }

  /**
   * @brief Processes data tracks
   *
   * @param collision Collision data
   * @param tracks Track data
   */
  void processData(ResoCollisionWithIndex::iterator const& collision,
                   soa::Filtered<aod::ResoTrackCandidates> const& tracks)
  {
    fillTrackTablesForCollision<false>(collision, tracks, tracksPerCollision);
  }
  PROCESS_SWITCH(ResonanceDaughterInitializer, processData, "Process tracks for data", false);

  /**
   * @brief Processes data tracks using the two-stage hybrid grouping
   *
   * GroupSlicer associates tracks automatically to the original
   * aod::Collision. Reduced collisions retain a scalar original-collision row
   * number and are explicitly sliced from the much smaller mapping table. The
   * tracks argument is already the selected slice and must not be sliced again.
   */
  void processDataHybrid(aod::Collision const& originalCollision,
                         SelectedResoCollisions const& reducedCollisions,
                         soa::Filtered<aod::ResoTrackCandidates> const& tracks)
  {
    auto reducedCollisionsThisCollision = reducedCollisions.sliceBy(reducedCollisionsPerOriginalCollision, originalCollision.globalIndex());
    if (reducedCollisionsThisCollision.size() == 0) {
      return;
    }
    if (reducedCollisionsThisCollision.size() > 1) {
      LOGF(error, "Found %zu reduced collisions for one original collision; skipping the ambiguous association", reducedCollisionsThisCollision.size());
      return;
    }
    auto reducedCollision = reducedCollisionsThisCollision.begin();
    fillTrackTables<false>(reducedCollision, tracks);
  }
  PROCESS_SWITCH(ResonanceDaughterInitializer, processDataHybrid, "Process data tracks with the two-stage hybrid grouping", false);

  /**
   * @brief Processes data tracks with configurable selected V0 and cascade gates
   */
  void processDataWithPairGate(ResoCollisionWithIndex::iterator const& collision,
                               soa::Filtered<aod::ResoTrackCandidates> const& tracks,
                               aod::ResoV0Candidates const& v0s,
                               aod::ResoCascadesCandidates const& cascades)
  {
    auto tracksThisCollision = tracks.sliceBy(tracksPerCollision, collision.collisionId());
    PairTrackSelection pairSelection{FilterForDerivedTables.cfgGlobalDaughterVeto.value};
    if (!preparePairTrackSelection<false>(collision, tracks, v0s, cascades,
                                          v0sPerCollision, cascadesPerCollision, pairSelection)) {
      return;
    }
    if (!hasTracksForEnabledOutputs<false>(collision, tracksThisCollision, pairSelection)) {
      return;
    }
    fillTrackTables<false>(collision, tracksThisCollision, pairSelection);
  }
  PROCESS_SWITCH(ResonanceDaughterInitializer, processDataWithPairGate, "Process data tracks with the configured pair-gate mode", false);

  /**
   * @brief Processes data tracks when a selectex0 and collision track exist
   *
   * This dedicated callback deliberately has no cascade input, so enabling
   * the V0 pair gate does not activate the cascade upstream dependency.
   */
  void processDataWithV0PairGate(ResoCollisionWithIndex::iterator const& collision,
                                 soa::Filtered<aod::ResoTrackCandidates> const& tracks,
                                 aod::ResoV0Candidates const& v0s)
  {
    auto tracksThisCollision = tracks.sliceBy(tracksPerCollision, collision.collisionId());
    auto v0sThisCollision = v0s.sliceBy(v0sPerCollision, collision.collisionId());
    PairTrackSelection pairSelection{FilterForDerivedTables.cfgGlobalDaughterVeto.value};
    pairSelection.useV0Candidates = true;
    pairSelection.v0Candidates =
      collectSelectedV0Daughters<false>(collision, v0sThisCollision, tracks, pairSelection.useGlobalDaughterVeto);
    if (!hasTracksForEnabledOutputs<false>(collision, tracksThisCollision, pairSelection)) {
      return;
    }
    fillTrackTables<false>(collision, tracksThisCollision, pairSelection);
  }
  PROCESS_SWITCH(ResonanceDaughterInitializer, processDataWithV0PairGate, "Process data tracks requiring a selected V0", false);

  /**
   * @brief Processes data tracks when a selected cascade and collision track exist
   *
   * This dedicated callback deliberately has no V0 input, so enabling the
   * cascade pair gate does not activate the V0 upstream dependency.
   */
  void processDataWithCascPairGate(ResoCollisionWithIndex::iterator const& collision,
                                   soa::Filtered<aod::ResoTrackCandidates> const& tracks,
                                   aod::ResoCascadesCandidates const& cascades)
  {
    auto tracksThisCollision = tracks.sliceBy(tracksPerCollision, collision.collisionId());
    auto cascadesThisCollision = cascades.sliceBy(cascadesPerCollision, collision.collisionId());
    PairTrackSelection pairSelection{FilterForDerivedTables.cfgGlobalDaughterVeto.value};
    pairSelection.useCascadeCandidates = true;
    pairSelection.cascadeCandidates =
      collectSelectedCascadeDaughters<false>(collision, cascadesThisCollision, tracks, pairSelection.useGlobalDaughterVeto);
    if (!hasTracksForEnabledOutputs<false>(collision, tracksThisCollision, pairSelection)) {
      return;
    }
    fillTrackTables<false>(collision, tracksThisCollision, pairSelection);
  }
  PROCESS_SWITCH(ResonanceDaughterInitializer, processDataWithCascPairGate, "Process data tracks requiring a selected cascade", false);

  /**
   * @brief Processes MC tracks
   *
   * @param collision Collision data
   * @param tracks Track data
   * @param mcParticles MC particles
   */
  void processMC(ResoCollisionWithIndex::iterator const& collision,
                 soa::Filtered<aod::ResoTrackCandidatesMC> const& tracks,
                 aod::McParticles const&)
  {
    fillTrackTablesForCollision<true>(collision, tracks, tracksMCPerCollision);
  }
  PROCESS_SWITCH(ResonanceDaughterInitializer, processMC, "Process tracks for MC", false);

  /**
   * @brief Processes MC tracks with configurable selected V0 and cascade gates
   */
  void processMCWithPairGate(ResoCollisionWithIndex::iterator const& collision,
                             soa::Filtered<aod::ResoTrackCandidatesMC> const& tracks,
                             aod::ResoV0CandidatesMC const& v0s,
                             aod::ResoCascadesCandidatesMC const& cascades,
                             aod::McParticles const&)
  {
    auto tracksThisCollision = tracks.sliceBy(tracksMCPerCollision, collision.collisionId());
    PairTrackSelection pairSelection{FilterForDerivedTables.cfgGlobalDaughterVeto.value};
    if (!preparePairTrackSelection<true>(collision, tracks, v0s, cascades,
                                         v0sMCPerCollision, cascadesMCPerCollision, pairSelection)) {
      return;
    }
    if (!hasTracksForEnabledOutputs<true>(collision, tracksThisCollision, pairSelection)) {
      return;
    }
    fillTrackTables<true>(collision, tracksThisCollision, pairSelection);
  }
  PROCESS_SWITCH(ResonanceDaughterInitializer, processMCWithPairGate, "Process MC tracks with the configured pair-gate mode", false);

  /**
   * @brief Processes MC tracks when a selected V0 and collision track exist
   */
  void processMCWithV0PairGate(ResoCollisionWithIndex::iterator const& collision,
                               soa::Filtered<aod::ResoTrackCandidatesMC> const& tracks,
                               aod::ResoV0CandidatesMC const& v0s,
                               aod::McParticles const&)
  {
    auto tracksThisCollision = tracks.sliceBy(tracksMCPerCollision, collision.collisionId());
    auto v0sThisCollision = v0s.sliceBy(v0sMCPerCollision, collision.collisionId());
    PairTrackSelection pairSelection{FilterForDerivedTables.cfgGlobalDaughterVeto.value};
    pairSelection.useV0Candidates = true;
    pairSelection.v0Candidates =
      collectSelectedV0Daughters<true>(collision, v0sThisCollision, tracks, pairSelection.useGlobalDaughterVeto);
    if (!hasTracksForEnabledOutputs<true>(collision, tracksThisCollision, pairSelection)) {
      return;
    }
    fillTrackTables<true>(collision, tracksThisCollision, pairSelection);
  }
  PROCESS_SWITCH(ResonanceDaughterInitializer, processMCWithV0PairGate, "Process MC tracks requiring a selected V0", false);

  /**
   * @brief Processes MC tracks when a selected cascade and collision track exist
   */
  void processMCWithCascPairGate(ResoCollisionWithIndex::iterator const& collision,
                                 soa::Filtered<aod::ResoTrackCandidatesMC> const& tracks,
                                 aod::ResoCascadesCandidatesMC const& cascades,
                                 aod::McParticles const&)
  {
    auto tracksThisCollision = tracks.sliceBy(tracksMCPerCollision, collision.collisionId());
    auto cascadesThisCollision = cascades.sliceBy(cascadesMCPerCollision, collision.collisionId());
    PairTrackSelection pairSelection{FilterForDerivedTables.cfgGlobalDaughterVeto.value};
    pairSelection.useCascadeCandidates = true;
    pairSelection.cascadeCandidates =
      collectSelectedCascadeDaughters<true>(collision, cascadesThisCollision, tracks, pairSelection.useGlobalDaughterVeto);
    if (!hasTracksForEnabledOutputs<true>(collision, tracksThisCollision, pairSelection)) {
      return;
    }
    fillTrackTables<true>(collision, tracksThisCollision, pairSelection);
  }
  PROCESS_SWITCH(ResonanceDaughterInitializer, processMCWithCascPairGate, "Process MC tracks requiring a selected cascade", false);

  /**
   * @brief Processes V0 data
   *
   * @param collision Collision data
   * @param v0s V0 data
   * @param tracks Track data
   */
  void processV0Data(ResoCollisionWithIndex::iterator const& collision, aod::ResoV0Candidates const& v0s, aod::ResoTrackCandidates const& tracks)
  {
    auto v0sThisCollision = v0s.sliceBy(v0sPerCollision, collision.collisionId());
    fillV0s<false>(collision, v0sThisCollision, tracks);
  }
  PROCESS_SWITCH(ResonanceDaughterInitializer, processV0Data, "Process V0s for data", false);

  /**
   * @brief Processes MC V0 data
   *
   * @param collision Collision data
   * @param v0s V0 data
   * @param tracks Track data
   */
  void processV0MC(ResoCollisionWithIndex::iterator const& collision, aod::ResoV0CandidatesMC const& v0s, aod::ResoTrackCandidatesMC const& tracks, aod::McParticles const&)
  {
    auto v0sThisCollision = v0s.sliceBy(v0sMCPerCollision, collision.collisionId());
    fillV0s<true>(collision, v0sThisCollision, tracks);
  }
  PROCESS_SWITCH(ResonanceDaughterInitializer, processV0MC, "Process V0s for MC", false);

  /**
   * @brief Processes cascade data
   *
   * @param collision Collision data
   * @param cascades Cascade data
   * @param tracks Track data
   */
  void processCascData(ResoCollisionWithIndex::iterator const& collision, aod::ResoCascadesCandidates const& cascades, aod::ResoTrackCandidates const& tracks)
  {
    auto cascadesThisCollision = cascades.sliceBy(cascadesPerCollision, collision.collisionId());
    fillCascades<false>(collision, cascadesThisCollision, tracks);
  }
  PROCESS_SWITCH(ResonanceDaughterInitializer, processCascData, "Process Cascades for data", false);

  /**
   * @brief Processes MC cascade data
   *
   * @param collision Collision data
   * @param cascades Cascade data
   * @param tracks Track data
   */
  void processCascMC(ResoCollisionWithIndex::iterator const& collision, aod::ResoCascadesCandidatesMC const& cascades, aod::ResoTrackCandidatesMC const& tracks, aod::McParticles const&)
  {
    auto cascadesThisCollision = cascades.sliceBy(cascadesMCPerCollision, collision.collisionId());
    fillCascades<true>(collision, cascadesThisCollision, tracks);
  }
  PROCESS_SWITCH(ResonanceDaughterInitializer, processCascMC, "Process Cascades for MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& context)
{
  return WorkflowSpec{
    adaptAnalysisTask<ResonanceModuleInitializer>(context),
    adaptAnalysisTask<ResonanceDaughterInitializer>(context)};
}
