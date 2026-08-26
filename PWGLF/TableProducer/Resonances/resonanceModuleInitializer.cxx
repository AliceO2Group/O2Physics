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
/// \since  Aug.18 2026

#include "PWGLF/DataModel/LFResonanceTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/mcCentrality.h"
#include "PWGLF/Utils/collisionCuts.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/CCDB/RCTSelectionFlags.h"
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
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

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
  static constexpr double MinimumChargedParticleCharge = 3.; // ROOT particle charge is stored in units of e/3
  static constexpr int MCCentralityRecoEstimator = 0;
  static constexpr int MCCentralityGeneratorEstimator = 1;
  static constexpr int MCCentralityImpactParameterEstimator = 2;
  static constexpr float MCVertexZMax = 10.f;

  int mRunNumber = 0;                        ///< Run number for the current data
  int multEstimator = 0;                     ///< Multiplicity estimator type
  float dBz = 0.f;                           ///< Magnetic field value
  float centrality = 0.f;                    ///< Centrality value for the event
  Service<o2::ccdb::BasicCCDBManager> ccdb;  ///< CCDB manager service
  Service<o2::framework::O2DatabasePDG> pdg; ///< PDG database service

  Produces<aod::ResoCollisions> resoCollisions;           ///< Output table for resonance collisions
  Produces<aod::ResoCollisionColls> resoCollisionColls;   ///< Output table for collision references
  Produces<aod::ResoCollisionGroups> resoCollisionGroups; ///< Canonical original-collision grouping references
  Produces<aod::ResoMCCollisions> resoMCCollisions;       ///< Output table for MC resonance collisions

  // CCDB options
  struct : ConfigurableGroup {
    Configurable<std::string> ccdbURL{"ccdbURL", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
    Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
    Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
    Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
    Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
    Configurable<bool> cfgFatalWhenNull{"cfgFatalWhenNull", true, "Fatal when null on ccdb access"};
    Configurable<bool> cfgBypassCollIndexFill{"cfgBypassCollIndexFill", false, "Unsupported in the modular workflow; must remain false"};
  } CCDB;

  // General event options
  struct : ConfigurableGroup {
    Configurable<double> dBzInput{"dBzInput", -999, "bz field, -999 is automatic"};
    Configurable<bool> cfgFillQA{"cfgFillQA", true, "Fill QA histograms"};
    Configurable<bool> cfgBypassCCDB{"cfgBypassCCDB", true, "Bypass loading CCDB part to save CPU time and memory"}; // will be affected to b_z value.
    Configurable<std::string> cfgMultName{"cfgMultName", "FT0M", "The name of multiplicity estimator"};
    Configurable<int> cfgCentralityMC{"cfgCentralityMC", 0, "Centrality estimator for MC (0: Reco, 1: MC, 2: impact parameter)"};
    ConfigurableAxis binsCent{"binsCent", {VARIABLE_WIDTH, 0., 0.01, 0.1, 1.0, 5.0, 10., 15., 20., 30., 40., 50., 70., 100.0, 105.}, "Binning of the centrality axis"};
    ConfigurableAxis cfgVtxBins{"cfgVtxBins", {VARIABLE_WIDTH, -20, -15, -10, -7, -5, -3, -2, -1, 0, 1, 2, 3, 5, 7, 10, 15, 20}, "Mixing bins - z-vertex"};
  } EventConfig;

  /// Event cuts
  o2::analysis::CollisonCuts colCuts;
  struct : ConfigurableGroup {
    Configurable<float> cfgEvtZvtx{"cfgEvtZvtx", 10.f, "Evt sel: Max. z-Vertex (cm)"};
    Configurable<int> cfgEvtOccupancyInTimeRange{"cfgEvtOccupancyInTimeRange", -1, "Evt sel: maximum track occupancy"};
    Configurable<bool> cfgEvtTriggerCheck{"cfgEvtTriggerCheck", false, "Evt sel: check for trigger"};
    Configurable<bool> cfgEvtOfflineCheck{"cfgEvtOfflineCheck", true, "Evt sel: check for offline selection"};
    Configurable<bool> cfgEvtTriggerTVXSel{"cfgEvtTriggerTVXSel", false, "Evt sel: triggerTVX selection (MB)"};
    Configurable<bool> cfgEvtTFBorderCut{"cfgEvtTFBorderCut", false, "Evt sel: apply TF border cut"};
    Configurable<bool> cfgEvtUseITSTPCvertex{"cfgEvtUseITSTPCvertex", false, "Evt sel: use at lease on ITS-TPC track for vertexing"};
    Configurable<bool> cfgEvtCollInTimeRangeNarrow{"cfgEvtCollInTimeRangeNarrow", false, "Evt sel: apply NoCollInTimeRangeNarrow"};
    Configurable<bool> cfgEvtZvertexTimedifference{"cfgEvtZvertexTimedifference", false, "Evt sel: apply Z-vertex time difference"};
    Configurable<bool> cfgEvtPileupRejection{"cfgEvtPileupRejection", false, "Evt sel: apply pileup rejection"};
    Configurable<bool> cfgEvtNoITSROBorderCut{"cfgEvtNoITSROBorderCut", false, "Evt sel: apply NoITSRO border cut"};
    Configurable<bool> cfgEvtRun2AliEventCuts{"cfgEvtRun2AliEventCuts", true, "Evt sel: apply Run2 AliEventCuts"};
    Configurable<bool> cfgEvtRun2INELgtZERO{"cfgEvtRun2INELgtZERO", false, "Evt sel: apply Run2 INELgtZERO"};
    Configurable<bool> cfgEvtUseRCTFlagChecker{"cfgEvtUseRCTFlagChecker", false, "Evt sel: use RCT flag checker"};
    Configurable<std::string> cfgEvtRCTFlagCheckerLabel{"cfgEvtRCTFlagCheckerLabel", "CBT_hadronPID", "Evt sel: RCT flag checker label"};
    Configurable<bool> cfgEvtRCTFlagCheckerZDCCheck{"cfgEvtRCTFlagCheckerZDCCheck", false, "Evt sel: RCT flag checker ZDC check"};
    Configurable<bool> cfgEvtRCTFlagCheckerLimitAcceptAsBad{"cfgEvtRCTFlagCheckerLimitAcceptAsBad", false, "Evt sel: RCT flag checker treat Limited Acceptance As Bad"};
  } EventCuts;
  RCTFlagsChecker rctChecker;

  HistogramRegistry qaRegistry{"QAHistos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Filter collisionFilter = nabs(aod::collision::posZ) < EventCuts.cfgEvtZvtx;

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
    // Determine the multiplicity estimator based on the configuration
    multEstimator = 0;
    if (EventConfig.cfgMultName.value == "FT0M") {
      multEstimator = 0;
    } else if (EventConfig.cfgMultName.value == "FT0C") {
      multEstimator = 1;
    } else if (EventConfig.cfgMultName.value == "FT0A") {
      multEstimator = 2;
    }
    LOGF(info, "Mult estimator: %d, %s", multEstimator, EventConfig.cfgMultName.value.c_str());

    // Ensure that only one process type is active at a time
    if (doprocessRun3 && doprocessRun2) {
      LOG(fatal) << "You cannot run both Run2 and Run3 processes at the same time";
    }
    if (doprocessRun2MC && doprocessRun3MC) {
      LOG(fatal) << "You cannot run both Run2 and Run3 MC processes at the same time";
    }
    if (CCDB.cfgBypassCollIndexFill) {
      LOG(fatal) << "cfgBypassCollIndexFill is incompatible with ResonanceDaughterInitializer";
    }

    // Initialize event selection cuts based on the process type
    if (doprocessRun2) {
      colCuts.setCuts(EventCuts.cfgEvtZvtx, EventCuts.cfgEvtTriggerCheck, EventCuts.cfgEvtOfflineCheck, false);
    } else if (doprocessRun3) {
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

    rctChecker.init(EventCuts.cfgEvtRCTFlagCheckerLabel, EventCuts.cfgEvtRCTFlagCheckerZDCCheck, EventCuts.cfgEvtRCTFlagCheckerLimitAcceptAsBad);

    // Configure CCDB access if not bypassed
    if (!EventConfig.cfgBypassCCDB) {
      ccdb->setURL(CCDB.ccdbURL.value);
      ccdb->setCaching(true);
      ccdb->setLocalObjectValidityChecking();
      ccdb->setFatalWhenNull(CCDB.cfgFatalWhenNull);
      uint64_t now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
      ccdb->setCreatedNotAfter(now); // TODO must become global parameter from the train creation time
    }

    // Initialize QA histograms if required
    if (EventConfig.cfgFillQA && (doprocessRun3MC || doprocessRun2MC)) {
      AxisSpec centAxis = {EventConfig.binsCent, "Centrality (%)"};
      AxisSpec idxMCAxis = {26, -0.5, 25.5, "Index"};
      qaRegistry.add("Event/hMCEventIndices", "hMCEventIndices", kTH2D, {centAxis, idxMCAxis});
    }
  }

  /**
   * @brief Initializes CCDB for a given BC
   *
   * @param bc BC iterator
   */
  void initCCDB(aod::BCsWithTimestamps::iterator const& bc) // Simple copy from LambdaKzeroFinder.cxx
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
   * @brief Checks if the collision is INEL>0
   *
   * @tparam MCPart Type of MC particles
   * @param mcparts MC particles
   * @return true if INEL>0, false otherwise
   */
  template <typename MCPart>
  bool isTrueINEL0(MCPart const& mcparts)
  {
    for (auto const& mcparticle : mcparts) {
      if (!mcparticle.isPhysicalPrimary()) {
        continue;
      }
      auto p = pdg->GetParticle(mcparticle.pdgCode());
      if (p != nullptr) {
        if (std::abs(p->Charge()) >= MinimumChargedParticleCharge) {
          if (std::abs(mcparticle.eta()) < 1) {
            return true;
          }
        }
      }
    }
    return false;
  }

  /**
   * @brief Centrality estimator selection
   *
   * @tparam ResoColl Type of resonance collision
   * @tparam isMC Boolean indicating if it's MC
   * @param resoEvents Resonance events
   * @return Centrality value
   */
  template <typename ResoColl, bool isMC = false>
  float centEst(ResoColl const& resoEvents)
  {
    float returnValue = -999.0;
    switch (multEstimator) {
      case 0:
        returnValue = resoEvents.centFT0M();
        break;
      case 1:
        if constexpr (isMC) {
          LOG(fatal) << "CentFT0C is not available for MC";
          return returnValue;
        } else {
          returnValue = resoEvents.centFT0C();
          break;
        }
      case 2:
        if constexpr (isMC) {
          LOG(fatal) << "CentFT0A is not available for MC";
          return returnValue;
        } else {
          returnValue = resoEvents.centFT0A();
          break;
        }
      default:
        returnValue = resoEvents.centFT0M();
        break;
    }
    return returnValue;
  }
  using GenMCCollisions = soa::Join<aod::McCollisions, aod::McCentFT0Ms, aod::MultsExtraMC>;
  float centEstMC(const GenMCCollisions::iterator& collision) { return centEst<GenMCCollisions::iterator, true>(collision); }

  /**
   * @brief Fills MC particles
   *
   * @tparam CollisionType Type of collision
   * @tparam SelectedMCPartType Type of selected MC particles
   * @tparam TotalMCParts Type of total MC particles
   * @param collision Collision data
   * @param mcParts Selected MC particles
   * @param mcParticles Total MC particles
   */
  template <typename CollisionType, typename SelectedMCPartType, typename TotalMCParts>
  void fillMCParticles(CollisionType collision, SelectedMCPartType const& mcParts, TotalMCParts const& mcParticles)
  {
    for (auto const& mcPart : mcParts) {
      std::vector<int> daughterPDGs;
      if (mcPart.has_daughters()) {
        auto daughter01 = mcParticles.rawIteratorAt(mcPart.daughtersIds()[0] - mcParticles.offset());
        auto daughter02 = mcParticles.rawIteratorAt(mcPart.daughtersIds()[1] - mcParticles.offset());
        daughterPDGs = {daughter01.pdgCode(), daughter02.pdgCode()};
      } else {
        daughterPDGs = {-1, -1};
      }
      reso2mcparents(collision.globalIndex(),
                     mcPart.globalIndex(),
                     mcPart.pdgCode(),
                     daughterPDGs[0], daughterPDGs[1],
                     mcPart.isPhysicalPrimary(),
                     mcPart.producedByGenerator(),
                     mcPart.pt(),
                     mcPart.px(),
                     mcPart.py(),
                     mcPart.pz(),
                     mcPart.eta(),
                     mcPart.phi(),
                     mcPart.y());
      daughterPDGs.clear();
    }
  }

  /**
   * @brief Fills MC collision data
   *
   * @tparam isRun2 Boolean indicating if it's Run2
   * @tparam MCCol Type of MC collision
   * @tparam MCPart Type of MC particles
   * @param mccol MC collision data
   * @param mcparts MC particles
   */
  template <bool isRun2, typename MCCol, typename MCPart>
  void fillMCCollision(MCCol const& mccol, MCPart const& mcparts)
  {
    const auto& mcColg = mccol.template mcCollision_as<GenMCCollisions>();
    float mcCent = 999.0;
    if constexpr (isRun2) {
      if (EventConfig.cfgCentralityMC == MCCentralityRecoEstimator) {
        mcCent = mccol.centRun2V0M();
      } else {
        mcCent = mcColg.impactParameter();
      }
    } else {
      if (EventConfig.cfgCentralityMC == MCCentralityRecoEstimator) {
        mcCent = centEst(mccol);
      } else if (EventConfig.cfgCentralityMC == MCCentralityGeneratorEstimator) {
        mcCent = centEstMC(mcColg);
      } else if (EventConfig.cfgCentralityMC == MCCentralityImpactParameterEstimator) {
        mcCent = mcColg.impactParameter();
      }
    }
    const bool inVtx10 = !(std::abs(mcColg.posZ()) > MCVertexZMax);
    bool isTrueINELgt0 = isTrueINEL0(mcparts);
    bool isTriggerTVX = mccol.selection_bit(aod::evsel::kIsTriggerTVX);
    bool isSel8 = mccol.sel8();
    bool isSelected = colCuts.isSelected(mccol, EventConfig.cfgFillQA);
    resoMCCollisions(inVtx10, isTrueINELgt0, isTriggerTVX, isSel8, isSelected, mcCent, -1.0f);

    if (EventConfig.cfgFillQA) {
      // QA for trigger efficiency
      qaRegistry.fill(HIST("Event/hMCEventIndices"), mcCent, aod::resocollision::kINEL);
      if (inVtx10) {
        qaRegistry.fill(HIST("Event/hMCEventIndices"), mcCent, aod::resocollision::kINEL10);
      }
      if (isTrueINELgt0) {
        qaRegistry.fill(HIST("Event/hMCEventIndices"), mcCent, aod::resocollision::kINELg0);
      }
      if (inVtx10 && isTrueINELgt0) {
        qaRegistry.fill(HIST("Event/hMCEventIndices"), mcCent, aod::resocollision::kINELg010);
      }

      // TVX MB trigger
      if (isTriggerTVX) {
        qaRegistry.fill(HIST("Event/hMCEventIndices"), mcCent, aod::resocollision::kTrig);
      }
      if (isTriggerTVX && inVtx10) {
        qaRegistry.fill(HIST("Event/hMCEventIndices"), mcCent, aod::resocollision::kTrig10);
      }
      if (isTriggerTVX && isTrueINELgt0) {
        qaRegistry.fill(HIST("Event/hMCEventIndices"), mcCent, aod::resocollision::kTrigINELg0);
      }
      if (isTriggerTVX && isTrueINELgt0 && inVtx10) {
        qaRegistry.fill(HIST("Event/hMCEventIndices"), mcCent, aod::resocollision::kTrigINELg010);
      }

      // Sel8 event selection
      if (isSel8) {
        qaRegistry.fill(HIST("Event/hMCEventIndices"), mcCent, aod::resocollision::kSel8);
      }
      if (isSel8 && inVtx10) {
        qaRegistry.fill(HIST("Event/hMCEventIndices"), mcCent, aod::resocollision::kSel810);
      }
      if (isSel8 && isTrueINELgt0) {
        qaRegistry.fill(HIST("Event/hMCEventIndices"), mcCent, aod::resocollision::kSel8INELg0);
      }
      if (isSel8 && isTrueINELgt0 && inVtx10) {
        qaRegistry.fill(HIST("Event/hMCEventIndices"), mcCent, aod::resocollision::kSel8INELg010);
      }

      // CollisionCuts selection
      if (isSelected) {
        qaRegistry.fill(HIST("Event/hMCEventIndices"), mcCent, aod::resocollision::kAllCuts);
      }
      if (isSelected && inVtx10) {
        qaRegistry.fill(HIST("Event/hMCEventIndices"), mcCent, aod::resocollision::kAllCuts10);
      }
      if (isSelected && isTrueINELgt0) {
        qaRegistry.fill(HIST("Event/hMCEventIndices"), mcCent, aod::resocollision::kAllCutsINELg0);
      }
      if (isSelected && isTrueINELgt0 && inVtx10) {
        qaRegistry.fill(HIST("Event/hMCEventIndices"), mcCent, aod::resocollision::kAllCutsINELg010);
      }
    }
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
  void processRun3(soa::Filtered<aod::ResoCollisionCandidates>::iterator const& collision,
                   aod::BCsWithTimestamps const&)
  {
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    initCCDB(bc);
    // Default event selection
    if (!colCuts.isSelected(collision, EventConfig.cfgFillQA)) {
      return;
    }
    if (EventCuts.cfgEvtUseRCTFlagChecker && !rctChecker(collision)) {
      return;
    }
    if (EventConfig.cfgFillQA) {
      colCuts.fillQA(collision);
    }
    const bool isRecINELgt0 = collision.isInelGt0();
    centrality = centEst(collision);

    resoCollisions(collision.multNTracksPV(), collision.multNTracksPVeta1(), collision.multNTracksPVetaHalf(), collision.posX(), collision.posY(), collision.posZ(), centEst(collision), dBz, isRecINELgt0);
    resoCollisionColls(collision.globalIndex());
    resoCollisionGroups(collision.globalIndex());
  }
  PROCESS_SWITCH(ResonanceModuleInitializer, processRun3, "Default process for RUN3", false);

  /**
   * @brief Processes Run2 data
   *
   * @param collision Collision data
   * @param bc BC data
   */
  void processRun2(soa::Filtered<aod::ResoRun2CollisionCandidates>::iterator const& collision,
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
    centrality = collision.centRun2V0M();

    resoCollisions(0, 0, 0, collision.posX(), collision.posY(), collision.posZ(), centrality, dBz, 0);
    resoCollisionColls(collision.globalIndex());
    resoCollisionGroups(collision.globalIndex());
  }
  PROCESS_SWITCH(ResonanceModuleInitializer, processRun2, "process for RUN2", false);

  /**
   * @brief Processes Run3 MC data
   *
   * @param collision Collision data
   * @param mcParticles MC particles
   * @param mcCollisions MC collisions
   */
  void processRun3MC(soa::Filtered<aod::ResoCollisionCandidatesMC>::iterator const& collision,
                     aod::McParticles const& mcParticles, GenMCCollisions const&)
  {
    if (EventCuts.cfgEvtUseRCTFlagChecker && !rctChecker(collision)) {
      return;
    }
    fillMCCollision<false>(collision, mcParticles);
  }
  PROCESS_SWITCH(ResonanceModuleInitializer, processRun3MC, "process MC for RUN3", false);

  /**
   * @brief Processes Run2 MC data
   *
   * @param collision Collision data
   * @param mcParticles MC particles
   */
  void processRun2MC(soa::Filtered<aod::ResoRun2CollisionCandidatesMC>::iterator const& collision,
                     aod::McParticles const& mcParticles)
  {
    fillMCCollision<true>(collision, mcParticles);
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
  static constexpr float MomentumQuantizationScale = 1000.f;
  static constexpr std::size_t StoredMCRelationCount = 2;

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
  Produces<aod::ResoMicroTracks> reso2microtrks;                      ///< Output table for resonance microtracks
  Produces<aod::ResoMicroTrackTracks> resoMicroTrackTracks;           ///< Output table for original microtrack row IDs
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
    Configurable<float> cfgCutMaxPt{"cfgCutMaxPt", 32.767f, "Maximum pT for tracks (GeV/c)"};
    Configurable<float> pidnSigmaPreSelectionCut{"pidnSigmaPreSelectionCut", 5.0f, "TPC PID cut (loose, improve performance)"};
    Configurable<int> mincrossedrows{"mincrossedrows", 70, "Minimum crossed rows for V0 daughter tracks"};
    Configurable<int> trackSelection{"trackSelection", 3, "Track selection: 0 -> No Cut, 1 -> kGlobalTrack, 2 -> kGlobalTrackWoPtEta, 3 -> kGlobalTrackWoDCA, 4 -> kQualityTracks, 5 -> kInAcceptanceTracks"};
    Configurable<double> cMaxDCArToPVcut{"cMaxDCArToPVcut", 2.0, "Track DCAr cut to PV Maximum"};
    Configurable<double> cMaxDCAzToPVcut{"cMaxDCAzToPVcut", 2.0, "Track DCAz cut to PV Maximum"};
    Configurable<double> cMinDCAzToPVcut{"cMinDCAzToPVcut", 0.0, "Track DCAz cut to PV Minimum"};
    Configurable<bool> cfgApplyTightDCAPtDepSelection{"cfgApplyTightDCAPtDepSelection", true, "Apply the pT-dependent tight DCA selection"};
    Configurable<float> cfgTightDCAOffset{"cfgTightDCAOffset", 0.004f, "Constant term of the tight DCA threshold (cm)"};
    Configurable<float> cfgTightDCAPtCoefficient{"cfgTightDCAPtCoefficient", 0.013f, "Coefficient of the pT-dependent tight DCA threshold"};
    Configurable<float> cfgTightDCAPtPower{"cfgTightDCAPtPower", 1.f, "Power in tight DCA = offset + coefficient / pT^power"};
  } TrackCuts;

  // V0 and V0-daughter cuts
  struct : ConfigurableGroup {
    Configurable<double> cMinV0PosDCArToPVcut{"cMinV0PosDCArToPVcut", 0.05f, "V0 Positive Track DCAr cut to PV Minimum"};
    Configurable<double> cMinV0NegDCArToPVcut{"cMinV0NegDCArToPVcut", 0.05f, "V0 Negative Track DCAr cut to PV Minimum"};
    Configurable<double> cMinV0Radius{"cMinV0Radius", 0.0, "Minimum V0 radius from PV"};
    Configurable<double> cMaxV0Radius{"cMaxV0Radius", 200.0, "Maximum V0 radius from PV"};
    Configurable<double> cMinV0CosPA{"cMinV0CosPA", 0.995, "Minimum V0 CosPA to PV"};
  } V0Cuts;

  // Cascade and cascade-daughter cuts
  struct : ConfigurableGroup {
    Configurable<int> cfgMinCrossedRowsCascBach{"cfgMinCrossedRowsCascBach", 70, "min crossed rows for bachelor track from cascade"};
    Configurable<double> cMinCascBachDCArToPVcut{"cMinCascBachDCArToPVcut", 0.05f, "Cascade Bachelor Track DCAr cut to PV Minimum"};
    Configurable<double> cMaxCascBachDCArToPVcut{"cMaxCascBachDCArToPVcut", 999.0f, "Cascade Bachelor Track DCAr cut to PV Maximum"};
    Configurable<double> cMaxCascDCAV0Daughters{"cMaxCascDCAV0Daughters", 1.6, "Cascade DCA between V0 daughters Maximum"};
    Configurable<double> cMaxCascDCACascDaughters{"cMaxCascDCACascDaughters", 1.6, "Cascade DCA between Casc daughters Maximum"};
    Configurable<double> cMinCascV0CosPA{"cMinCascV0CosPA", 0.97, "Minimum Cascade V0 CosPA to PV"};
    Configurable<double> cMaxCascV0Radius{"cMaxCascV0Radius", 200.0, "Maximum Cascade V0 radius from PV"};
    Configurable<double> cMinCascV0Radius{"cMinCascV0Radius", 0.0, "Minimum Cascade V0 radius from PV"};
    Configurable<double> cMinCascRadius{"cMinCascRadius", 0.0, "Minimum Cascade radius from PV"};
    Configurable<double> cMaxCascRadius{"cMaxCascRadius", 200.0, "Maximum Cascade radius from PV"};
    Configurable<double> cMinCascCosPA{"cMinCascCosPA", 0.97, "Minimum Cascade CosPA to PV"};
    Configurable<double> cCascMassResol{"cCascMassResol", 999, "Cascade mass resolution"};
  } CascadeCuts;

  // Derived dataset selections
  struct : ConfigurableGroup {
    Configurable<bool> cfgFillPionTracks{"cfgFillPionTracks", false, "Fill pion tracks"};
    Configurable<bool> cfgFillKaonTracks{"cfgFillKaonTracks", false, "Fill kaon tracks"};
    Configurable<bool> cfgFillProtonTracks{"cfgFillProtonTracks", false, "Fill proton tracks"};
    Configurable<bool> cfgFillPionMicroTracks{"cfgFillPionMicroTracks", false, "Fill pion micro tracks"};
    Configurable<bool> cfgFillKaonMicroTracks{"cfgFillKaonMicroTracks", false, "Fill kaon micro tracks"};
    Configurable<bool> cfgFillProtonMicroTracks{"cfgFillProtonMicroTracks", false, "Fill proton micro tracks"};
    Configurable<bool> cfgFillPionUltraMicroTracks{"cfgFillPionUltraMicroTracks", true, "Fill pion ultra micro tracks"};
    Configurable<bool> cfgFillKaonUltraMicroTracks{"cfgFillKaonUltraMicroTracks", false, "Fill kaon ultra micro tracks"};
    Configurable<bool> cfgFillProtonUltraMicroTracks{"cfgFillProtonUltraMicroTracks", false, "Fill proton ultra micro tracks"};
    Configurable<bool> cfgFillK0s{"cfgFillK0s", false, "Fill K0s"};
    Configurable<bool> cfgFillLambda0{"cfgFillLambda0", false, "Fill Lambda0"};
    Configurable<bool> cfgBypassNoPairV0s{"cfgBypassNoPairV0s", false, "In a *WithPairGate process, bypass track fill if no V0 passes the configured selections"};
    Configurable<bool> cfgBypassNoPairCascades{"cfgBypassNoPairCascades", true, "In a *WithPairGate process, bypass track fill if no cascade passes the configured selections"};
    Configurable<bool> cfgFillMicroTracks{"cfgFillMicroTracks", false, "Fill micro tracks"};
    Configurable<bool> cfgFillUltraMicroTracks{"cfgFillUltraMicroTracks", false, "Fill ultra micro tracks"};
    Configurable<bool> cfgBypassTrackFill{"cfgBypassTrackFill", false, "Bypass the full ResoTracks table fill"};
    Configurable<bool> cfgBypassTrackIndexFill{"cfgBypassTrackIndexFill", false, "Bypass original daughter ID table fill"};
  } FilterForDerivedTables;

  // Secondary selections for K0s and Lambda0
  struct : ConfigurableGroup {
    Configurable<bool> cfgSecondaryRequire{"cfgSecondaryRequire", false, "Secondary cuts on/off"};
    Configurable<bool> cfgSecondaryArmenterosCut{"cfgSecondaryArmenterosCut", false, "cut on Armenteros-Podolanski graph"};
    Configurable<bool> cfgSecondaryCrossMassHypothesisCut{"cfgSecondaryCrossMassHypothesisCut", false, "Apply cut based on the lambda mass hypothesis"};
    Configurable<bool> cfgByPassDauPIDSelection{"cfgByPassDauPIDSelection", true, "Bypass TPC PID preselection for V0 daughters"};
    Configurable<float> cfgSecondaryDauDCAMax{"cfgSecondaryDauDCAMax", 0.2, "Maximum DCA Secondary daughters to PV"};
    Configurable<float> cfgSecondaryDauPosDCAtoPVMin{"cfgSecondaryDauPosDCAtoPVMin", 0.0, "Minimum DCA Secondary positive daughters to PV"};
    Configurable<float> cfgSecondaryDauNegDCAtoPVMin{"cfgSecondaryDauNegDCAtoPVMin", 0.0, "Minimum DCA Secondary negative daughters to PV"};
    Configurable<float> cfgSecondaryPtMin{"cfgSecondaryPtMin", 0.f, "Minimum transverse momentum of Secondary"};
    Configurable<float> cfgSecondaryRapidityMax{"cfgSecondaryRapidityMax", 0.5, "Maximum rapidity of Secondary"};
    Configurable<float> cfgSecondaryRadiusMin{"cfgSecondaryRadiusMin", 0.0, "Minimum transverse radius of Secondary"};
    Configurable<float> cfgSecondaryRadiusMax{"cfgSecondaryRadiusMax", 999.9, "Maximum transverse radius of Secondary"};
    Configurable<float> cfgSecondaryCosPAMin{"cfgSecondaryCosPAMin", 0.998, "Mininum cosine pointing angle of Secondary"};
    Configurable<float> cfgSecondaryDCAtoPVMax{"cfgSecondaryDCAtoPVMax", 0.4, "Maximum DCA Secondary to PV"};
    Configurable<float> cfgSecondaryProperLifetimeMax{"cfgSecondaryProperLifetimeMax", 20., "Maximum Secondary Lifetime"};
    Configurable<float> cfgSecondaryparamArmenterosCut{"cfgSecondaryparamArmenterosCut", 0.2, "parameter for Armenteros Cut"};
    Configurable<float> cfgSecondaryMassWindow{"cfgSecondaryMassWindow", 0.03, "Secondary inv mass selection window"};
    Configurable<float> cfgSecondaryCrossMassCutWindow{"cfgSecondaryCrossMassCutWindow", 0.05, "Secondary inv mass selection window with (anti)lambda hypothesis"};
  } SecondaryCuts;

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

  // The daughter task needs this row-wise mapping back to the original collision.
  // Keep ResonanceModuleInitializer::cfgBypassCollIndexFill disabled and enable
  // the matching Run 2/Run 3 base event process for MC workflows.
  using ResoCollisionWithIndex = soa::Join<aod::ResoCollisions, aod::ResoCollisionColls>;
  using SelectedResoCollisions = soa::Join<aod::ResoCollisions, aod::ResoCollisionGroups>;

  /**
   * @brief Initializes the task
   *
   * @param context Initialization context
   */
  void init(InitContext&)
  {
    const bool processTrackDataEnabled = doprocessData || doprocessDataHybrid || doprocessDataWithPairGate;
    const bool processTrackMCEnabled = doprocessMC || doprocessMCWithPairGate;
    const bool processV0DataEnabled = doprocessV0Data || doprocessV0DataHybrid;
    const bool processCascDataEnabled = doprocessCascData || doprocessCascDataHybrid;
    const int enabledTrackProcesses = static_cast<int>(doprocessData) +
                                      static_cast<int>(doprocessDataHybrid) +
                                      static_cast<int>(doprocessDataWithPairGate) +
                                      static_cast<int>(doprocessMC) +
                                      static_cast<int>(doprocessMCWithPairGate);

    if (enabledTrackProcesses > 1) {
      LOGF(fatal, "Only one track process can be enabled in ResonanceDaughterInitializer");
    }
    if (static_cast<int>(doprocessV0Data) + static_cast<int>(doprocessV0DataHybrid) + static_cast<int>(doprocessV0MC) > 1) {
      LOGF(fatal, "Only one V0 process can be enabled in ResonanceDaughterInitializer");
    }
    if (static_cast<int>(doprocessCascData) + static_cast<int>(doprocessCascDataHybrid) + static_cast<int>(doprocessCascMC) > 1) {
      LOGF(fatal, "Only one cascade process can be enabled in ResonanceDaughterInitializer");
    }
    if ((doprocessData || doprocessDataHybrid || doprocessMC) &&
        (FilterForDerivedTables.cfgBypassNoPairV0s || FilterForDerivedTables.cfgBypassNoPairCascades)) {
      LOGF(warn, "Pair-gate options are ignored by processData/processDataHybrid/processMC; enable the matching *WithPairGate process to apply them");
    }
    if (doprocessDataWithPairGate && FilterForDerivedTables.cfgBypassNoPairV0s && !processV0DataEnabled) {
      LOGF(fatal, "cfgBypassNoPairV0s requires processV0Data or processV0DataHybrid so an accepted V0 is written for every retained collision");
    }
    if (doprocessDataWithPairGate && FilterForDerivedTables.cfgBypassNoPairCascades && !processCascDataEnabled) {
      LOGF(fatal, "cfgBypassNoPairCascades requires processCascData or processCascDataHybrid so an accepted cascade is written for every retained collision");
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

    if (TrackCuts.cfgApplyTightDCAPtDepSelection.value &&
        (!std::isfinite(TrackCuts.cfgTightDCAOffset.value) ||
         !std::isfinite(TrackCuts.cfgTightDCAPtCoefficient.value) ||
         !std::isfinite(TrackCuts.cfgTightDCAPtPower.value) ||
         TrackCuts.cfgTightDCAOffset.value < 0.f ||
         TrackCuts.cfgTightDCAPtCoefficient.value < 0.f ||
         TrackCuts.cfgTightDCAPtPower.value < 0.f)) {
      LOGF(fatal, "Tight-DCA offset, pT coefficient, and power must be finite and non-negative");
    }

    if (FilterForDerivedTables.cfgFillUltraMicroTracks) {
      constexpr float MaxQuantizedMomentum = static_cast<float>(std::numeric_limits<int16_t>::max()) / MomentumQuantizationScale;
      const float maxLongitudinalMomentum = TrackCuts.cfgCutMaxPt.value * std::sinh(TrackCuts.cfgCutEta.value);
      if (TrackCuts.cfgCutMaxPt.value > MaxQuantizedMomentum ||
          !std::isfinite(maxLongitudinalMomentum) ||
          maxLongitudinalMomentum > MaxQuantizedMomentum) {
        LOGF(fatal, "Ultra-micro cfgCutMaxPt/cfgCutEta allow momentum components beyond the int16_t quantization range");
      }
      int enabledUltraMicroSpecies = 0;
      enabledUltraMicroSpecies += FilterForDerivedTables.cfgFillPionUltraMicroTracks ? 1 : 0;
      enabledUltraMicroSpecies += FilterForDerivedTables.cfgFillKaonUltraMicroTracks ? 1 : 0;
      enabledUltraMicroSpecies += FilterForDerivedTables.cfgFillProtonUltraMicroTracks ? 1 : 0;
      if (enabledUltraMicroSpecies != 1) {
        LOGF(fatal, "Exactly one pion/kaon/proton PID species must be enabled when filling ultra-micro tracks");
      }
      if (TrackCuts.pidnSigmaPreSelectionCut.value > o2::aod::resoultramicrodaughter::PidNSigma::MaxNSigma) {
        LOGF(fatal, "Ultra-micro PID encoding requires pidnSigmaPreSelectionCut <= 5");
      }
      if (FilterForDerivedTables.cfgFillKaonUltraMicroTracks) {
        ultraMicroPidSpecies = UltraMicroPidSpecies::Kaon;
      } else if (FilterForDerivedTables.cfgFillProtonUltraMicroTracks) {
        ultraMicroPidSpecies = UltraMicroPidSpecies::Proton;
      } else {
        ultraMicroPidSpecies = UltraMicroPidSpecies::Pion;
      }
    }

    if (cfgFillQA) {
      AxisSpec idxAxis = {8, 0.0, 8.0, "Index"};
      AxisSpec ptAxis = {100, 0.0f, 10.0f, "#it{p}_{T} (GeV/#it{c})"};
      // The DCA maps cover the full configured pT range in 0.1 GeV/c steps.
      constexpr float DcaPtBinWidth = 0.1f;
      const int maxDcaPtBin = static_cast<int>(std::ceil(TrackCuts.cfgCutMaxPt.value / DcaPtBinWidth));
      const int nDcaPtBins = maxDcaPtBin + 1;
      const float dcaPtAxisHalfBin = 0.5f * DcaPtBinWidth;
      AxisSpec dcaPtAxis = {nDcaPtBins,
                            -dcaPtAxisHalfBin,
                            maxDcaPtBin * DcaPtBinWidth + dcaPtAxisHalfBin,
                            "#it{p}_{T} (GeV/#it{c})"};
      AxisSpec etaAxis = {100, -1.0f, 1.0f, "#eta"};
      AxisSpec phiAxis = {100, 0.0f, TwoPI, "#phi"};
      // Keep the configured signed DCA limits at bin centres so tracks exactly
      // on an accepted selection boundary are not moved to overflow.
      constexpr int NDcaBins = 201;
      constexpr float DcaAxisPaddingFraction = 1.f / 200.f;
      const auto configuredDcaXYMax = static_cast<float>(TrackCuts.cMaxDCArToPVcut.value);
      const auto configuredDcaZMax = static_cast<float>(TrackCuts.cMaxDCAzToPVcut.value);
      const float dcaXYAxisHalfRange = configuredDcaXYMax > 0.f ? configuredDcaXYMax * (1.f + DcaAxisPaddingFraction) : 1.e-4f;
      const float dcaZAxisHalfRange = configuredDcaZMax > 0.f ? configuredDcaZMax * (1.f + DcaAxisPaddingFraction) : 1.e-4f;
      AxisSpec dcaXYAxis = {NDcaBins, -dcaXYAxisHalfRange, dcaXYAxisHalfRange, "DCA_{xy} (cm)"};
      AxisSpec dcaZAxis = {NDcaBins, -dcaZAxisHalfRange, dcaZAxisHalfRange, "DCA_{z} (cm)"};
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
        if (processTrackMCEnabled) {
          qaRegistry.add("QA/hGoodMCTrackIndices", "hGoodMCTrackIndices", kTH1D, {idxAxis});
        }
        if (cfgDetailTrackQA) {
          if (!FilterForDerivedTables.cfgBypassTrackFill) {
            qaRegistry.add("QA/h4TrackPtEtaPhi", "ResoTracks pT, eta, phi", kTH3F, {ptAxis, etaAxis, phiAxis});
            qaRegistry.add("QA/h2TrackDCAxyVsPt", "ResoTracks DCAxy vs pT", kTH2F, {dcaPtAxis, dcaXYAxis});
            qaRegistry.add("QA/h2TrackDCAzVsPt", "ResoTracks DCAz vs pT", kTH2F, {dcaPtAxis, dcaZAxis});
            qaRegistry.add("QA/h4TrackTPCnSigma", "ResoTracks TPC nSigma Pi, Ka, Pr as pT", kTHnSparseD, {ptAxis, nSigmaTPCAxis, nSigmaTPCAxis, nSigmaTPCAxis});
            qaRegistry.add("QA/h4TrackTOFnSigma", "ResoTracks TOF nSigma Pi, Ka, Pr as pT", kTHnSparseD, {ptAxis, nSigmaTOFAxis, nSigmaTOFAxis, nSigmaTOFAxis});
          }
          if (FilterForDerivedTables.cfgFillMicroTracks) {
            qaRegistry.add("QA/h4MicroTrackPtEtaPhi", "ResoMicroTracks pT, eta, phi", kTH3F, {ptAxis, etaAxis, phiAxis});
            qaRegistry.add("QA/h2MicroTrackDCAxyVsPt", "ResoMicroTracks DCAxy vs pT", kTH2F, {dcaPtAxis, dcaXYAxis});
            qaRegistry.add("QA/h2MicroTrackDCAzVsPt", "ResoMicroTracks DCAz vs pT", kTH2F, {dcaPtAxis, dcaZAxis});
            qaRegistry.add("QA/h4MicroTrackTPCnSigma", "ResoMicroTracks TPC nSigma Pi, Ka, Pr as pT", kTHnSparseD, {ptAxis, nSigmaTPCAxis, nSigmaTPCAxis, nSigmaTPCAxis});
            qaRegistry.add("QA/h4MicroTrackTOFnSigma", "ResoMicroTracks TOF nSigma Pi, Ka, Pr as pT", kTHnSparseD, {ptAxis, nSigmaTOFAxis, nSigmaTOFAxis, nSigmaTOFAxis});
          }
          if (FilterForDerivedTables.cfgFillUltraMicroTracks) {
            qaRegistry.add("QA/h4UltraMicroTrackPtEtaPhi", "ResoUltraMicroTracks pT, eta, phi", kTH3F, {ptAxis, etaAxis, phiAxis});
            qaRegistry.add("QA/h2UltraMicroTrackDCAxyVsPt", "ResoUltraMicroTracks DCAxy vs pT", kTH2F, {dcaPtAxis, dcaXYAxis});
            qaRegistry.add("QA/h2UltraMicroTrackDCAzVsPt", "ResoUltraMicroTracks DCAz vs pT", kTH2F, {dcaPtAxis, dcaZAxis});
            qaRegistry.add("QA/h4UltraMicroTrackTPCnSigma", "ResoUltraMicroTracks TPC nSigma Pi, Ka, Pr as pT", kTHnSparseD, {ptAxis, nSigmaTPCAxis, nSigmaTPCAxis, nSigmaTPCAxis});
            qaRegistry.add("QA/h4UltraMicroTrackTOFnSigma", "ResoUltraMicroTracks TOF nSigma Pi, Ka, Pr as pT", kTHnSparseD, {ptAxis, nSigmaTOFAxis, nSigmaTOFAxis, nSigmaTOFAxis});
          }
        }
      }

      if (processV0DataEnabled || doprocessV0MC) {
        qaRegistry.add("QA/hGoodV0Indices", "hGoodV0Indices", kTH1D, {idxAxis});
        if (doprocessV0MC) {
          qaRegistry.add("QA/hGoodMCV0Indices", "hGoodMCV0Indices", kTH1D, {idxAxis});
        }
        AxisSpec radiusAxis = {100, 0.0, 200.0, "V0 Radius"};
        AxisSpec cosPAAxis = {100, 0.995, 1.0, "V0 CosPA"};
        qaRegistry.add("QA/hV0Radius", "V0 Radius", kTH1F, {radiusAxis});
        qaRegistry.add("QA/hV0CosPA", "V0 CosPA", kTH1F, {cosPAAxis});
      }

      if (processCascDataEnabled || doprocessCascMC) {
        AxisSpec radiusAxis = {100, 0.0, 200.0, "Cascade Radius"};
        AxisSpec cosPAAxis = {100, 0.97, 1.0, "Cascade CosPA"};
        qaRegistry.add("QA/hGoodCascIndices", "hGoodCascIndices", kTH1D, {idxAxis});
        if (doprocessCascMC) {
          qaRegistry.add("QA/hGoodMCCascIndices", "hGoodMCCascIndices", kTH1D, {idxAxis});
        }
        qaRegistry.add("QA/hCascRadius", "Cascade Radius", kTH1F, {radiusAxis});
        qaRegistry.add("QA/hCascCosPA", "Cascade CosPA", kTH1F, {cosPAAxis});
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

    // Check if the module is initialized with both data and MC
    if ((processTrackDataEnabled && processTrackMCEnabled) || (processV0DataEnabled && doprocessV0MC) || (processCascDataEnabled && doprocessCascMC)) {
      LOGF(fatal, "ResonanceDaughterInitializer initialized with both data and MC");
    }
    // Check if none of the processes are enabled
    if (!doprocessDummy && !processTrackDataEnabled && !processTrackMCEnabled && !processV0DataEnabled && !doprocessV0MC && !processCascDataEnabled && !doprocessCascMC) {
      LOGF(fatal, "ResonanceDaughterInitializer not initialized, enable at least one process");
    }
  }
  template <typename T>
  bool filterMicroTrack(T const& track)
  {
    // if no selection is requested, return true
    if (!FilterForDerivedTables.cfgFillPionMicroTracks && !FilterForDerivedTables.cfgFillKaonMicroTracks && !FilterForDerivedTables.cfgFillProtonMicroTracks) {
      return true;
    }
    if (FilterForDerivedTables.cfgFillPionMicroTracks) {
      if (std::abs(track.tpcNSigmaPi()) < TrackCuts.pidnSigmaPreSelectionCut) {
        return true;
      }
    }
    if (FilterForDerivedTables.cfgFillKaonMicroTracks) {
      if (std::abs(track.tpcNSigmaKa()) < TrackCuts.pidnSigmaPreSelectionCut) {
        return true;
      }
    }
    if (FilterForDerivedTables.cfgFillProtonMicroTracks) {
      if (std::abs(track.tpcNSigmaPr()) < TrackCuts.pidnSigmaPreSelectionCut) {
        return true;
      }
    }
    return false;
  }

  template <typename T>
  bool filterUltraMicroTrack(T const& track)
  {
    switch (ultraMicroPidSpecies) {
      case UltraMicroPidSpecies::Pion:
        return std::abs(track.tpcNSigmaPi()) < TrackCuts.pidnSigmaPreSelectionCut;
      case UltraMicroPidSpecies::Kaon:
        return std::abs(track.tpcNSigmaKa()) < TrackCuts.pidnSigmaPreSelectionCut;
      case UltraMicroPidSpecies::Proton:
        return std::abs(track.tpcNSigmaPr()) < TrackCuts.pidnSigmaPreSelectionCut;
    }
    return false;
  }

  template <typename T>
  bool filterTrack(T const& track)
  {
    // if no selection is requested, return true
    if (!FilterForDerivedTables.cfgFillPionTracks && !FilterForDerivedTables.cfgFillKaonTracks && !FilterForDerivedTables.cfgFillProtonTracks) {
      return true;
    }
    if (FilterForDerivedTables.cfgFillPionTracks) {
      if (std::abs(track.tpcNSigmaPi()) < TrackCuts.pidnSigmaPreSelectionCut) {
        return true;
      }
    }
    if (FilterForDerivedTables.cfgFillKaonTracks) {
      if (std::abs(track.tpcNSigmaKa()) < TrackCuts.pidnSigmaPreSelectionCut) {
        return true;
      }
    }
    if (FilterForDerivedTables.cfgFillProtonTracks) {
      if (std::abs(track.tpcNSigmaPr()) < TrackCuts.pidnSigmaPreSelectionCut) {
        return true;
      }
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
        std::abs(v0.dcapostopv()) < SecondaryCuts.cfgSecondaryDauPosDCAtoPVMin ||
        std::abs(v0.dcanegtopv()) < SecondaryCuts.cfgSecondaryDauNegDCAtoPVMin ||
        v0.pt() < SecondaryCuts.cfgSecondaryPtMin ||
        v0.v0radius() < SecondaryCuts.cfgSecondaryRadiusMin ||
        v0.v0radius() > SecondaryCuts.cfgSecondaryRadiusMax ||
        v0.dcav0topv() > SecondaryCuts.cfgSecondaryDCAtoPVMax ||
        v0.v0cosPA() < SecondaryCuts.cfgSecondaryCosPAMin) {
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
    bool selected = false;
    if (FilterForDerivedTables.cfgFillK0s) {
      const bool passesK0DaughterPID = bypassDaughterPID ||
                                       (std::abs(posTrack.tpcNSigmaPi()) < TrackCuts.pidnSigmaPreSelectionCut &&
                                        std::abs(negTrack.tpcNSigmaPi()) < TrackCuts.pidnSigmaPreSelectionCut);
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
                                   (std::abs(posTrack.tpcNSigmaPr()) < TrackCuts.pidnSigmaPreSelectionCut &&
                                    std::abs(negTrack.tpcNSigmaPi()) < TrackCuts.pidnSigmaPreSelectionCut);
      const bool passesAntiLambdaPID = bypassDaughterPID ||
                                       (std::abs(posTrack.tpcNSigmaPi()) < TrackCuts.pidnSigmaPreSelectionCut &&
                                        std::abs(negTrack.tpcNSigmaPr()) < TrackCuts.pidnSigmaPreSelectionCut);
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
  bool isTrackSelected(CollisionType const&, TrackType const& track)
  {
    if (cfgFillQA) {
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
    if (cfgFillQA) {
      qaRegistry.fill(HIST("QA/hGoodTrackIndices"), 1.5);
    }
    if (std::fabs(track.dcaZ()) > TrackCuts.cMaxDCAzToPVcut || std::fabs(track.dcaZ()) < TrackCuts.cMinDCAzToPVcut) {
      return false;
    }
    if (cfgFillQA) {
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
    if (posTrack.tpcNClsCrossedRows() < TrackCuts.mincrossedrows || negTrack.tpcNClsCrossedRows() < TrackCuts.mincrossedrows) {
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
    if (std::abs(casc.mXi() - MassXiMinus) > CascadeCuts.cCascMassResol) {
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

  /// @brief Check whether a collision has at least one V0 that would be written
  template <bool isMC, typename CollisionType, typename V0Type, typename TrackType>
  bool hasSelectedV0(CollisionType const& collision, V0Type const& v0s, TrackType const& tracks)
  {
    for (auto const& v0 : v0s) {
      if (isV0Selected<isMC>(collision, v0, tracks, false) && filterV0(collision, v0, tracks)) {
        return true;
      }
    }
    return false;
  }

  /// @brief Check whether a collision has at least one cascade that would be written
  template <bool isMC, typename CollisionType, typename CascType, typename TrackType>
  bool hasSelectedCascade(CollisionType const& collision, CascType const& cascades, TrackType const& tracks)
  {
    for (auto const& casc : cascades) {
      if (isCascSelected<isMC>(collision, casc, tracks, false)) {
        return true;
      }
    }
    return false;
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

  template <bool isMC, typename TrackType, typename CollisionType>
  void fillUltraMicroTracks(CollisionType const& collision, TrackType const& tracks)
  {
    // Loop over tracks
    for (auto const& track : tracks) {
      if (!isMicroTrackSelected<isMC>(collision, track)) {
        continue;
      }
      if (!filterUltraMicroTrack(track)) {
        continue;
      }
      if (!o2::aod::resoultramicrodaughter::DCAEncoding::isValid(track.dcaXY()) ||
          !o2::aod::resoultramicrodaughter::DCAEncoding::isValid(track.dcaZ())) {
        continue;
      }
      o2::aod::resoultramicrodaughter::DCAEncoding dcaEncoding(track.dcaXY(), track.dcaZ());
      int16_t px1000 = 0;
      int16_t py1000 = 0;
      int16_t pz1000 = 0;
      if (!quantizeP(track.px(), px1000) ||
          !quantizeP(track.py(), py1000) ||
          !quantizeP(track.pz(), pz1000)) {
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
   */
  template <bool isMC, typename TrackType, typename CollisionType>
  void fillMicroTracks(CollisionType const& collision, TrackType const& tracks)
  {
    // Loop over tracks
    for (auto const& track : tracks) {
      if (!isMicroTrackSelected<isMC>(collision, track)) {
        continue;
      }
      if (!filterMicroTrack(track)) {
        continue;
      }
      o2::aod::resomicrodaughter::ResoMicroTrackSelFlag trackSelFlag(track.dcaXY(), track.dcaZ());
      if (TrackCuts.cfgApplyTightDCAPtDepSelection) {
        const float dcaThreshold = tightDCAThreshold(track.pt());
        if (dcaThreshold >= 0.f && std::abs(track.dcaXY()) < dcaThreshold) {
          trackSelFlag.setDCAxy0();
        }
        if (dcaThreshold >= 0.f && std::abs(track.dcaZ()) < dcaThreshold) {
          trackSelFlag.setDCAz0();
        }
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
        qaRegistry.fill(HIST("QA/h4MicroTrackPtEtaPhi"), track.pt(), track.eta(), track.phi());
        qaRegistry.fill(HIST("QA/h2MicroTrackDCAxyVsPt"), track.pt(), track.dcaXY());
        qaRegistry.fill(HIST("QA/h2MicroTrackDCAzVsPt"), track.pt(), track.dcaZ());
        qaRegistry.fill(HIST("QA/h4MicroTrackTPCnSigma"), track.pt(), track.tpcNSigmaPi(), track.tpcNSigmaKa(), track.tpcNSigmaPr());
        if (track.hasTOF()) {
          qaRegistry.fill(HIST("QA/h4MicroTrackTOFnSigma"), track.pt(), track.tofNSigmaPi(), track.tofNSigmaKa(), track.tofNSigmaPr());
        }
      }
      reso2microtrks(collision.globalIndex(),
                     track.px(),
                     track.py(),
                     track.pz(),
                     static_cast<uint8_t>(o2::aod::resomicrodaughter::PidNSigma(std::abs(track.tpcNSigmaPi()), std::abs(track.tofNSigmaPi()), track.hasTOF())),
                     static_cast<uint8_t>(o2::aod::resomicrodaughter::PidNSigma(std::abs(track.tpcNSigmaKa()), std::abs(track.tofNSigmaKa()), track.hasTOF())),
                     static_cast<uint8_t>(o2::aod::resomicrodaughter::PidNSigma(std::abs(track.tpcNSigmaPr()), std::abs(track.tofNSigmaPr()), track.hasTOF())),
                     static_cast<uint8_t>(trackSelFlag),
                     trackFlags);
      if (!FilterForDerivedTables.cfgBypassTrackIndexFill) {
        resoMicroTrackTracks(track.globalIndex());
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
  template <bool isMC, typename TrackType, typename CollisionType>
  void fillTracks(CollisionType const& collision, TrackType const& tracks)
  {
    if (FilterForDerivedTables.cfgBypassTrackFill) {
      return;
    }
    // Loop over tracks
    for (auto const& track : tracks) {
      if (!isTrackSelected<isMC>(collision, track)) {
        continue;
      }
      if (!filterTrack(track)) {
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
                static_cast<int16_t>(std::round(track.dcaXY() * 10000)),
                static_cast<int16_t>(std::round(track.dcaZ() * 10000)),
                static_cast<int8_t>(std::round(track.tpcNSigmaPi() * 10)),
                static_cast<int8_t>(std::round(track.tpcNSigmaKa() * 10)),
                static_cast<int8_t>(std::round(track.tpcNSigmaPr() * 10)),
                static_cast<int8_t>(std::round(track.tofNSigmaPi() * 10)),
                static_cast<int8_t>(std::round(track.tofNSigmaKa() * 10)),
                static_cast<int8_t>(std::round(track.tofNSigmaPr() * 10)),
                static_cast<int16_t>(std::round(track.tpcSignal() * 100)),
                trackFlags);
      if (!FilterForDerivedTables.cfgBypassTrackIndexFill) {
        resoTrackTracks(track.globalIndex());
      }
      if constexpr (isMC) {
        fillMCTrack(track);
      }
    }
  }

  /**
   * @brief Fills MC track data
   *
   * @tparam TrackType Type of track
   * @param track Track data
   */
  template <typename TrackType>
  void fillMCTrack(TrackType const& track)
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
      reso2mctracks(particle.pdgCode(),
                    mothers[0],
                    motherPDGs[0],
                    siblings.data(),
                    particle.isPhysicalPrimary(),
                    particle.producedByGenerator());
    } else {
      // No MC particle associated
      reso2mctracks(0,
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
      if (cfgFillQA) {
        qaRegistry.fill(HIST("QA/hV0Radius"), v0.v0radius());
        qaRegistry.fill(HIST("QA/hV0CosPA"), v0.v0cosPA());
      }
      const std::array<int, 2> childIDs{v0.posTrackId(), v0.negTrackId()}; // Original track IDs for downstream pair-level shared-daughter rejection
      reso2v0s(collision.globalIndex(),
               v0.pt(),
               v0.px(),
               v0.py(),
               v0.pz(),
               childIDs.data(),
               (int8_t)(v0.template posTrack_as<TrackType>().tpcNSigmaPi() * 10),
               (int8_t)(v0.template posTrack_as<TrackType>().tpcNSigmaKa() * 10),
               (int8_t)(v0.template posTrack_as<TrackType>().tpcNSigmaPr() * 10),
               (int8_t)(v0.template negTrack_as<TrackType>().tpcNSigmaPi() * 10),
               (int8_t)(v0.template negTrack_as<TrackType>().tpcNSigmaKa() * 10),
               (int8_t)(v0.template negTrack_as<TrackType>().tpcNSigmaPr() * 10),
               (int8_t)(v0.template posTrack_as<TrackType>().tofNSigmaPi() * 10),
               (int8_t)(v0.template posTrack_as<TrackType>().tofNSigmaKa() * 10),
               (int8_t)(v0.template posTrack_as<TrackType>().tofNSigmaPr() * 10),
               (int8_t)(v0.template negTrack_as<TrackType>().tofNSigmaPi() * 10),
               (int8_t)(v0.template negTrack_as<TrackType>().tofNSigmaKa() * 10),
               (int8_t)(v0.template negTrack_as<TrackType>().tofNSigmaPr() * 10),
               v0.v0cosPA(),
               v0.dcaV0daughters(),
               v0.dcapostopv(),
               v0.dcanegtopv(),
               v0.dcav0topv(),
               static_cast<uint8_t>(v0.template posTrack_as<TrackType>().tpcNClsCrossedRows()),
               static_cast<uint8_t>(v0.template negTrack_as<TrackType>().tpcNClsCrossedRows()),
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
      if (daughters.size() > StoredMCRelationCount) {
        LOGF(info, "daughters.size() is larger than 2");
      }
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
      if (cfgFillQA) {
        qaRegistry.fill(HIST("QA/hCascRadius"), casc.cascradius());
        qaRegistry.fill(HIST("QA/hCascCosPA"), casc.casccosPA(collision.posX(), collision.posY(), collision.posZ()));
      }
      const std::array<int, 3> childIDs{casc.posTrackId(), casc.negTrackId(), casc.bachelorId()}; // Original track IDs for downstream pair-level shared-daughter rejection
      reso2cascades(collision.globalIndex(),
                    casc.pt(),
                    casc.px(),
                    casc.py(),
                    casc.pz(),
                    childIDs.data(),
                    (int8_t)(casc.template posTrack_as<TrackType>().tpcNSigmaPi() * 10),
                    (int8_t)(casc.template posTrack_as<TrackType>().tpcNSigmaKa() * 10),
                    (int8_t)(casc.template posTrack_as<TrackType>().tpcNSigmaPr() * 10),
                    (int8_t)(casc.template negTrack_as<TrackType>().tpcNSigmaPi() * 10),
                    (int8_t)(casc.template negTrack_as<TrackType>().tpcNSigmaKa() * 10),
                    (int8_t)(casc.template negTrack_as<TrackType>().tpcNSigmaPr() * 10),
                    (int8_t)(casc.template bachelor_as<TrackType>().tpcNSigmaPi() * 10),
                    (int8_t)(casc.template bachelor_as<TrackType>().tpcNSigmaKa() * 10),
                    (int8_t)(casc.template bachelor_as<TrackType>().tpcNSigmaPr() * 10),
                    (int8_t)(casc.template posTrack_as<TrackType>().tofNSigmaPi() * 10),
                    (int8_t)(casc.template posTrack_as<TrackType>().tofNSigmaKa() * 10),
                    (int8_t)(casc.template posTrack_as<TrackType>().tofNSigmaPr() * 10),
                    (int8_t)(casc.template negTrack_as<TrackType>().tofNSigmaPi() * 10),
                    (int8_t)(casc.template negTrack_as<TrackType>().tofNSigmaKa() * 10),
                    (int8_t)(casc.template negTrack_as<TrackType>().tofNSigmaPr() * 10),
                    (int8_t)(casc.template bachelor_as<TrackType>().tofNSigmaPi() * 10),
                    (int8_t)(casc.template bachelor_as<TrackType>().tofNSigmaKa() * 10),
                    (int8_t)(casc.template bachelor_as<TrackType>().tofNSigmaPr() * 10),
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
                    static_cast<uint8_t>(casc.template posTrack_as<TrackType>().tpcNClsCrossedRows()),
                    static_cast<uint8_t>(casc.template negTrack_as<TrackType>().tpcNClsCrossedRows()),
                    static_cast<uint8_t>(casc.template bachelor_as<TrackType>().tpcNClsCrossedRows()),
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
      if (daughters.size() > StoredMCRelationCount) {
        LOGF(info, "daughters.size() is larger than 2");
      }
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
  void processDummy(aod::ResoCollision const&)
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
  template <bool isMC, typename CollisionType, typename TrackTableType>
  void fillTrackTables(CollisionType const& collision, TrackTableType const& tracks)
  {
    fillTracks<isMC>(collision, tracks);
    if (FilterForDerivedTables.cfgFillMicroTracks) {
      fillMicroTracks<isMC>(collision, tracks);
    }
    if (FilterForDerivedTables.cfgFillUltraMicroTracks) {
      fillUltraMicroTracks<isMC>(collision, tracks);
    }
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
   * The canonical fIndexCollisions column in ResoCollisionGroups lets
   * GroupSlicer associate both reduced collisions and tracks to the same
   * original aod::Collision. The tracks argument is therefore already the
   * selected slice for this collision and must not be sliced again.
   */
  void processDataHybrid(aod::Collision const&,
                         soa::SmallGroups<SelectedResoCollisions> const& reducedCollisions,
                         soa::Filtered<aod::ResoTrackCandidates> const& tracks)
  {
    if (reducedCollisions.size() == 0) {
      return;
    }
    if (reducedCollisions.size() != 1) {
      LOGF(fatal, "Expected exactly one reduced collision for an original collision, found %zu", reducedCollisions.size());
    }
    auto reducedCollision = reducedCollisions.begin();
    fillTrackTables<false>(reducedCollision, tracks);
  }
  PROCESS_SWITCH(ResonanceDaughterInitializer, processDataHybrid, "Process data tracks with the two-stage hybrid grouping", false);

  /**
   * @brief Processes data tracks with configurable selected-V0 and selected-cascade gates
   */
  void processDataWithPairGate(ResoCollisionWithIndex::iterator const& collision,
                               soa::Filtered<aod::ResoTrackCandidates> const& tracks,
                               aod::ResoV0Candidates const& v0s,
                               aod::ResoCascadesCandidates const& cascades)
  {
    if (FilterForDerivedTables.cfgBypassNoPairV0s) {
      auto v0sThisCollision = v0s.sliceBy(v0sPerCollision, collision.collisionId());
      if (!hasSelectedV0<false>(collision, v0sThisCollision, tracks)) {
        return;
      }
    }
    if (FilterForDerivedTables.cfgBypassNoPairCascades) {
      auto cascadesThisCollision = cascades.sliceBy(cascadesPerCollision, collision.collisionId());
      if (!hasSelectedCascade<false>(collision, cascadesThisCollision, tracks)) {
        return;
      }
    }
    fillTrackTablesForCollision<false>(collision, tracks, tracksPerCollision);
  }
  PROCESS_SWITCH(ResonanceDaughterInitializer, processDataWithPairGate, "Process data tracks with configurable pair gates", false);

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
   * @brief Processes MC tracks with configurable V0 and cascade candidate gates
   */
  void processMCWithPairGate(ResoCollisionWithIndex::iterator const& collision,
                             soa::Filtered<aod::ResoTrackCandidatesMC> const& tracks,
                             aod::ResoV0CandidatesMC const& v0s,
                             aod::ResoCascadesCandidatesMC const& cascades,
                             aod::McParticles const&)
  {
    auto v0sThisCollision = v0s.sliceBy(v0sMCPerCollision, collision.collisionId());
    if (FilterForDerivedTables.cfgBypassNoPairV0s && v0sThisCollision.size() < 1) {
      return;
    }
    auto cascadesThisCollision = cascades.sliceBy(cascadesMCPerCollision, collision.collisionId());
    if (FilterForDerivedTables.cfgBypassNoPairCascades && cascadesThisCollision.size() < 1) {
      return;
    }
    fillTrackTablesForCollision<true>(collision, tracks, tracksMCPerCollision);
  }
  PROCESS_SWITCH(ResonanceDaughterInitializer, processMCWithPairGate, "Process MC tracks with configurable pair gates", false);

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
   * @brief Processes data V0s grouped automatically by their original collision
   *
   * Both V0s and tracks are already restricted to the current original
   * aod::Collision by GroupSlicer. The unfiltered track table is required for
   * resolving the positive and negative daughter indices.
   */
  void processV0DataHybrid(aod::Collision const&,
                           soa::SmallGroups<SelectedResoCollisions> const& reducedCollisions,
                           aod::ResoV0Candidates const& v0s,
                           aod::ResoTrackCandidates const& tracks)
  {
    if (reducedCollisions.size() == 0) {
      return;
    }
    if (reducedCollisions.size() != 1) {
      LOGF(fatal, "Expected exactly one reduced collision for an original collision, found %zu", reducedCollisions.size());
    }
    auto reducedCollision = reducedCollisions.begin();
    fillV0s<false>(reducedCollision, v0s, tracks);
  }
  PROCESS_SWITCH(ResonanceDaughterInitializer, processV0DataHybrid, "Process data V0s with the two-stage hybrid grouping", false);

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
   * @brief Processes data cascades grouped automatically by their original collision
   *
   * Cascades and tracks arrive as original-collision groups. Keeping this as a
   * separate callback from V0 processing avoids enabling either upstream input
   * dependency unless its process switch is selected.
   */
  void processCascDataHybrid(aod::Collision const&,
                             soa::SmallGroups<SelectedResoCollisions> const& reducedCollisions,
                             aod::ResoCascadesCandidates const& cascades,
                             aod::ResoTrackCandidates const& tracks)
  {
    if (reducedCollisions.size() == 0) {
      return;
    }
    if (reducedCollisions.size() != 1) {
      LOGF(fatal, "Expected exactly one reduced collision for an original collision, found %zu", reducedCollisions.size());
    }
    auto reducedCollision = reducedCollisions.begin();
    fillCascades<false>(reducedCollision, cascades, tracks);
  }
  PROCESS_SWITCH(ResonanceDaughterInitializer, processCascDataHybrid, "Process data cascades with the two-stage hybrid grouping", false);

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
