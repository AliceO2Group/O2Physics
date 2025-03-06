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
/// \author Bong-Hwi Lim <bong-hwi.lim@cern.ch>

#include <string>
#include <vector>
#include "CCDB/BasicCCDBManager.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/Qvectors.h"
#include "Common/Core/EventPlaneHelper.h"
#include "CommonConstants/PhysicsConstants.h"
#include "CommonConstants/MathConstants.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DetectorsBase/Propagator.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/LFResonanceTables.h"
#include "PWGLF/DataModel/mcCentrality.h"
#include "PWGLF/Utils/collisionCuts.h"
#include "ReconstructionDataFormats/Track.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::constants::physics;
using namespace o2::constants::math;

/**
 * @brief Initializer for the event pool for resonance study
 *
 * This struct is responsible for initializing and processing collision data
 * for resonance studies. It handles event selection, centrality estimation,
 * and QA histogram filling.
 */
struct ResonanceModuleInitializer {
  int mRunNumber;                            ///< Run number for the current data
  int multEstimator;                         ///< Multiplicity estimator type
  float dBz;                                 ///< Magnetic field value
  float centrality;                          ///< Centrality value for the event
  Service<o2::ccdb::BasicCCDBManager> ccdb;  ///< CCDB manager service
  Service<o2::framework::O2DatabasePDG> pdg; ///< PDG database service
  EventPlaneHelper helperEP;                 ///< Helper for event plane calculations

  Produces<aod::ResoCollisions> resoCollisions;             ///< Output table for resonance collisions
  Produces<aod::ResoMCCollisions> resoMCCollisions;         ///< Output table for MC resonance collisions
  Produces<aod::ResoSpheroCollisions> resoSpheroCollisions; ///< Output table for spherocity
  Produces<aod::ResoEvtPlCollisions> resoEvtPlCollisions;   ///< Output table for event plane

  // CCDB options
  Configurable<std::string> ccdbURL{"ccdbURL", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};

  Configurable<bool> cfgFatalWhenNull{"cfgFatalWhenNull", true, "Fatal when null on ccdb access"};

  // Configurables
  Configurable<double> dBzInput{"dBzInput", -999, "bz field, -999 is automatic"};
  Configurable<bool> cfgFillQA{"cfgFillQA", false, "Fill QA histograms"};
  Configurable<bool> cfgBypassCCDB{"cfgBypassCCDB", true, "Bypass loading CCDB part to save CPU time and memory"}; // will be affected to b_z value.
  Configurable<std::string> cfgMultName{"cfgMultName", "FT0M", "The name of multiplicity estimator"};
  Configurable<int> cfgCentralityMC{"cfgCentralityMC", 0, "Centrality estimator for MC (0: Reco, 1: MC, 2: impact parameter)"};

  // EventCorrection for MC
  ConfigurableAxis binsCent{"binsCent", {VARIABLE_WIDTH, 0., 0.01, 0.1, 1.0, 5.0, 10., 15., 20., 30., 40., 50., 70., 100.0, 105.}, "Binning of the centrality axis"};
  ConfigurableAxis cfgVtxBins{"cfgVtxBins", {VARIABLE_WIDTH, -20, -15, -10, -7, -5, -3, -2, -1, 0, 1, 2, 3, 5, 7, 10, 15, 20}, "Mixing bins - z-vertex"};

  /// Event cuts
  o2::analysis::CollisonCuts colCuts;
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

  // Spherocity configuration
  Configurable<int> cfgTrackSphMin{"cfgTrackSphMin", 10, "Number of tracks for Spherocity Calculation"};
  Configurable<int> cfgTrackSphDef{"cfgTrackSphDef", 0, "Spherocity Definition: |pT| = 1 -> 0, otherwise -> 1"};

  // Qvector configuration
  Configurable<int> cfgEvtPl{"cfgEvtPl", 40500, "Configuration of three subsystems for the event plane and its resolution, 10000*RefA + 100*RefB + S, where FT0C:0, FT0A:1, FT0M:2, FV0A:3, BPos:5, BNeg:6"};

  int evtPlRefAId = static_cast<int>(cfgEvtPl / 10000);
  int evtPlRefBId = static_cast<int>((cfgEvtPl - evtPlRefAId * 10000) / 100);
  int evtPlDetId = cfgEvtPl - evtPlRefAId * 10000 - evtPlRefBId * 100;

  HistogramRegistry qaRegistry{"QAHistos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Filter collisionFilter = nabs(aod::collision::posZ) < cfgEvtZvtx;

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
    if (cfgMultName.value == "FT0M") {
      multEstimator = 0;
    } else if (cfgMultName.value == "FT0C") {
      multEstimator = 1;
    } else if (cfgMultName.value == "FT0A") {
      multEstimator = 2;
    }
    LOGF(info, "Mult estimator: %d, %s", multEstimator, cfgMultName.value.c_str());

    // Ensure that only one process type is active at a time
    if (doprocessRun3 && doprocessRun2) {
      LOG(fatal) << "You cannot run both Run2 and Run3 processes at the same time";
    }
    if (doprocessRun2MC && doprocessRun3MC) {
      LOG(fatal) << "You cannot run both Run2 and Run3 MC processes at the same time";
    }

    // Initialize event selection cuts based on the process type
    if (doprocessRun2) {
      colCuts.setCuts(cfgEvtZvtx, cfgEvtTriggerCheck, cfgEvtOfflineCheck, false);
    } else if (doprocessRun3) {
      colCuts.setCuts(cfgEvtZvtx, cfgEvtTriggerCheck, cfgEvtOfflineCheck, true, false, cfgEvtOccupancyInTimeRange);
    }
    colCuts.init(&qaRegistry);
    colCuts.setTriggerTVX(cfgEvtTriggerTVXSel);
    colCuts.setApplyTFBorderCut(cfgEvtTFBorderCut);
    colCuts.setApplyITSTPCvertex(cfgEvtUseITSTPCvertex);
    colCuts.setApplyCollInTimeRangeNarrow(cfgEvtCollInTimeRangeNarrow);
    colCuts.setApplyZvertexTimedifference(cfgEvtZvertexTimedifference);
    colCuts.setApplyPileupRejection(cfgEvtPileupRejection);
    colCuts.setApplyNoITSROBorderCut(cfgEvtNoITSROBorderCut);
    colCuts.setApplyRun2AliEventCuts(cfgEvtRun2AliEventCuts);
    colCuts.setApplyRun2INELgtZERO(cfgEvtRun2INELgtZERO);

    // Configure CCDB access if not bypassed
    if (!cfgBypassCCDB) {
      ccdb->setURL(ccdbURL.value);
      ccdb->setCaching(true);
      ccdb->setLocalObjectValidityChecking();
      ccdb->setFatalWhenNull(cfgFatalWhenNull);
      uint64_t now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
      ccdb->setCreatedNotAfter(now); // TODO must become global parameter from the train creation time
    }

    // Initialize QA histograms if required
    if (doprocessRun3MC || doprocessRun2MC) {
      AxisSpec centAxis = {binsCent, "Centrality (%)"};
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
    if (cfgBypassCCDB)
      return;
    if (mRunNumber == bc.runNumber()) {
      return;
    }

    // In case override, don't proceed, please - no CCDB access required
    if (dBzInput > -990) {
      dBz = dBzInput;
      ;
      o2::parameters::GRPMagField grpmag;
      if (std::fabs(dBz) > 1e-5) {
        grpmag.setL3Current(30000.f / (dBz / 5.0f));
      }
      o2::base::Propagator::initFieldFromGRP(&grpmag);
      mRunNumber = bc.runNumber();
      return;
    }

    auto run3grpTimestamp = bc.timestamp();
    o2::parameters::GRPObject* grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>(grpPath, run3grpTimestamp);
    o2::parameters::GRPMagField* grpmag = 0x0;
    if (grpo) {
      o2::base::Propagator::initFieldFromGRP(grpo);
      // Fetch magnetic field from ccdb for current collision
      dBz = grpo->getNominalL3Field();
      LOG(info) << "Retrieved GRP for timestamp " << run3grpTimestamp << " with magnetic field of " << dBz << " kZG";
    } else {
      grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, run3grpTimestamp);
      if (!grpmag) {
        LOG(fatal) << "Got nullptr from CCDB for path " << grpmagPath << " of object GRPMagField and " << grpPath << " of object GRPObject for timestamp " << run3grpTimestamp;
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
      if (!mcparticle.isPhysicalPrimary())
        continue;
      auto p = pdg->GetParticle(mcparticle.pdgCode());
      if (p != nullptr) {
        if (std::abs(p->Charge()) >= 3) {
          if (std::abs(mcparticle.eta()) < 1)
            return true;
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
   * @param ResoEvents Resonance events
   * @return Centrality value
   */
  template <typename ResoColl, bool isMC = false>
  float centEst(ResoColl ResoEvents)
  {
    float returnValue = -999.0;
    switch (multEstimator) {
      case 0:
        returnValue = ResoEvents.centFT0M();
        break;
      case 1:
        if constexpr (isMC) {
          LOG(fatal) << "CentFT0C is not available for MC";
          return returnValue;
        } else {
          returnValue = ResoEvents.centFT0C();
          break;
        }
      case 2:
        if constexpr (isMC) {
          LOG(fatal) << "CentFT0A is not available for MC";
          return returnValue;
        } else {
          returnValue = ResoEvents.centFT0A();
          break;
        }
      default:
        returnValue = ResoEvents.centFT0M();
        break;
    }
    return returnValue;
  }
  using GenMCCollisions = soa::Join<aod::McCollisions, aod::McCentFT0Ms, aod::MultsExtraMC>;
  float centEstMC(const GenMCCollisions::iterator& collision) { return centEst<GenMCCollisions::iterator, true>(collision); }

  /**
   * @brief Computes the spherocity of an event
   *
   * @tparam T Type of the tracks
   * @param tracks All tracks
   * @param nTracksMin Minimum number of tracks
   * @param spdef Spherocity definition
   * @return Spherocity value
   */
  template <typename T>
  float computeSpherocity(T const& tracks, int nTracksMin, int spdef)
  {
    // if number of tracks is not enough for spherocity estimation.
    int ntrks = tracks.size();
    if (ntrks < nTracksMin)
      return -99.;

    // start computing spherocity

    float ptSum = 0.;
    for (auto const& track : tracks) {
      if (cfgFillQA) {
        qaRegistry.fill(HIST("Phi"), track.phi());
      }
      if (spdef == 0) {
        ptSum += 1.;
      } else {
        ptSum += track.pt();
      }
    }

    float tempSph = 1.;
    for (int i = 0; i < 360 / 0.1; ++i) {
      float sum = 0., pt = 0.;
      float phiparm = (PI * i * 0.1) / 180.;
      float nx = std::cos(phiparm);
      float ny = std::sin(phiparm);
      for (auto const& trk : tracks) {
        pt = trk.pt();
        if (spdef == 0) {
          pt = 1.;
        }
        float phi = trk.phi();
        float px = pt * std::cos(phi);
        float py = pt * std::sin(phi);
        // sum += pt * abs(sin(phiparm - phi));
        sum += std::abs(px * ny - py * nx);
      }
      float sph = std::pow((sum / ptSum), 2);
      if (sph < tempSph)
        tempSph = sph;
    }

    return std::pow(PIHalf, 2) * tempSph;
  }

  /**
   * @brief Gets the event plane
   *
   * @tparam ResoColl Type of resonance collision
   * @param ResoEvents Resonance events
   * @return Event plane value
   */
  template <typename ResoColl>
  float getEvtPl(ResoColl ResoEvents)
  {
    float returnValue = -999.0;
    if (ResoEvents.qvecAmp()[evtPlDetId] > 1e-8)
      returnValue = helperEP.GetEventPlane(ResoEvents.qvecRe()[evtPlDetId * 4 + 3], ResoEvents.qvecIm()[evtPlDetId * 4 + 3], 2);
    return returnValue;
  }

  /**
   * @brief Gets the event plane resolution
   *
   * @tparam ResoColl Type of resonance collision
   * @param ResoEvents Resonance events
   * @param a First index
   * @param b Second index
   * @return Event plane resolution
   */
  template <typename ResoColl>
  float getEvtPlRes(ResoColl ResoEvents, int a, int b)
  {
    float returnValue = -999.0;
    if (ResoEvents.qvecAmp()[a] < 1e-8 || ResoEvents.qvecAmp()[b] < 1e-8)
      return returnValue;
    returnValue = helperEP.GetResolution(helperEP.GetEventPlane(ResoEvents.qvecRe()[a * 4 + 3], ResoEvents.qvecIm()[a * 4 + 3], 2), helperEP.GetEventPlane(ResoEvents.qvecRe()[b * 4 + 3], ResoEvents.qvecIm()[b * 4 + 3], 2), 2);
    return returnValue;
  }

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
      if (cfgCentralityMC == 0) {
        mcCent = mccol.centRun2V0M();
      } else {
        mcCent = mcColg.impactParameter();
      }
    } else {
      if (cfgCentralityMC == 0) {
        mcCent = centEst(mccol);
      } else if (cfgCentralityMC == 1) {
        mcCent = centEstMC(mcColg);
      } else if (cfgCentralityMC == 2) {
        mcCent = mcColg.impactParameter();
      }
    }
    bool inVtx10 = (std::abs(mcColg.posZ()) > 10.) ? false : true;
    bool isTrueINELgt0 = isTrueINEL0(mcparts);
    bool isTriggerTVX = mccol.selection_bit(aod::evsel::kIsTriggerTVX);
    bool isSel8 = mccol.sel8();
    bool isSelected = colCuts.isSelected(mccol);
    resoMCCollisions(inVtx10, isTrueINELgt0, isTriggerTVX, isSel8, isSelected, mcCent);

    // QA for Trigger efficiency
    qaRegistry.fill(HIST("Event/hMCEventIndices"), mcCent, aod::resocollision::kINEL);
    if (inVtx10)
      qaRegistry.fill(HIST("Event/hMCEventIndices"), mcCent, aod::resocollision::kINEL10);
    if (isTrueINELgt0)
      qaRegistry.fill(HIST("Event/hMCEventIndices"), mcCent, aod::resocollision::kINELg0);
    if (inVtx10 && isTrueINELgt0)
      qaRegistry.fill(HIST("Event/hMCEventIndices"), mcCent, aod::resocollision::kINELg010);

    // TVX MB trigger
    if (isTriggerTVX)
      qaRegistry.fill(HIST("Event/hMCEventIndices"), mcCent, aod::resocollision::kTrig);
    if (isTriggerTVX && inVtx10)
      qaRegistry.fill(HIST("Event/hMCEventIndices"), mcCent, aod::resocollision::kTrig10);
    if (isTriggerTVX && isTrueINELgt0)
      qaRegistry.fill(HIST("Event/hMCEventIndices"), mcCent, aod::resocollision::kTrigINELg0);
    if (isTriggerTVX && isTrueINELgt0 && inVtx10)
      qaRegistry.fill(HIST("Event/hMCEventIndices"), mcCent, aod::resocollision::kTrigINELg010);

    // Sel8 event selection
    if (isSel8)
      qaRegistry.fill(HIST("Event/hMCEventIndices"), mcCent, aod::resocollision::kSel8);
    if (isSel8 && inVtx10)
      qaRegistry.fill(HIST("Event/hMCEventIndices"), mcCent, aod::resocollision::kSel810);
    if (isSel8 && isTrueINELgt0)
      qaRegistry.fill(HIST("Event/hMCEventIndices"), mcCent, aod::resocollision::kSel8INELg0);
    if (isSel8 && isTrueINELgt0 && inVtx10)
      qaRegistry.fill(HIST("Event/hMCEventIndices"), mcCent, aod::resocollision::kSel8INELg010);

    // CollisionCuts selection
    if (isSelected)
      qaRegistry.fill(HIST("Event/hMCEventIndices"), mcCent, aod::resocollision::kAllCuts);
    if (isSelected && inVtx10)
      qaRegistry.fill(HIST("Event/hMCEventIndices"), mcCent, aod::resocollision::kAllCuts10);
    if (isSelected && isTrueINELgt0)
      qaRegistry.fill(HIST("Event/hMCEventIndices"), mcCent, aod::resocollision::kAllCutsINELg0);
    if (isSelected && isTrueINELgt0 && inVtx10)
      qaRegistry.fill(HIST("Event/hMCEventIndices"), mcCent, aod::resocollision::kAllCutsINELg010);
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
    if (!colCuts.isSelected(collision))
      return;
    colCuts.fillQA(collision);
    centrality = centEst(collision);

    resoCollisions(collision.globalIndex(), 0, collision.posX(), collision.posY(), collision.posZ(), centrality, dBz);
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
    if (!colCuts.isSelected(collision))
      return;
    colCuts.fillQARun2(collision);
    centrality = collision.centRun2V0M();

    resoCollisions(collision.globalIndex(), 0, collision.posX(), collision.posY(), collision.posZ(), centrality, dBz);
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

  /**
   * @brief Processes Spherocity
   *
   * @param collision Collision data
   * @param tracks Track data
   */
  void processSpherocity(soa::Filtered<aod::ResoCollisionCandidates>::iterator const& collision, aod::ResoTrackCandidates const& tracks)
  {
    float spherocity = computeSpherocity(tracks, cfgTrackSphMin, cfgTrackSphDef);
    resoSpheroCollisions(collision.globalIndex(), spherocity);
  }
  PROCESS_SWITCH(ResonanceModuleInitializer, processSpherocity, "process Spherocity", false);

  /**
   * @brief Processes Event Plane
   *
   * @param collision Collision data with Qvectors
   * @param tracks Track data
   */
  void processEventPlane(soa::Filtered<soa::Join<aod::ResoCollisionCandidates, aod::Qvectors>>::iterator const& collision)
  {
    resoEvtPlCollisions(collision.globalIndex(), getEvtPl(collision), getEvtPlRes(collision, evtPlDetId, evtPlRefAId), getEvtPlRes(collision, evtPlDetId, evtPlRefBId), getEvtPlRes(collision, evtPlRefAId, evtPlRefBId));
  }
  PROCESS_SWITCH(ResonanceModuleInitializer, processEventPlane, "process Event Plane", false);
};

/**
 * @brief Initializer for the resonance daughters producer
 *
 * This struct initializes and processes daughters for resonance studies.
 * It applies daughter selection criteria and fills QA histograms for daughter properties.
 */
struct ResonanceDaughterInitializer {
  SliceCache cache;
  Produces<aod::ResoTracks> reso2trks;           ///< Output table for resonance tracks
  Produces<aod::ResoMCTracks> reso2mctracks;     ///< Output table for MC resonance tracks
  Produces<aod::ResoV0s> reso2v0s;               ///< Output table for resonance V0s
  Produces<aod::ResoMCV0s> reso2mcv0s;           ///< Output table for MC resonance V0s
  Produces<aod::ResoCascades> reso2cascades;     ///< Output table for resonance cascades
  Produces<aod::ResoMCCascades> reso2mccascades; ///< Output table for MC resonance cascades

  // Configurables
  Configurable<bool> cfgFillQA{"cfgFillQA", false, "Fill QA histograms"};

  // Configurables for tracks
  Configurable<float> cMaxDCArToPVcut{"cMaxDCArToPVcut", 2.0, "Track DCAr cut to PV Maximum"};
  Configurable<float> cMinDCArToPVcut{"cMinDCArToPVcut", 0.0, "Track DCAr cut to PV Minimum"};
  Configurable<float> cMaxDCAzToPVcut{"cMaxDCAzToPVcut", 2.0, "Track DCAz cut to PV Maximum"};
  Configurable<float> cMinDCAzToPVcut{"cMinDCAzToPVcut", 0.0, "Track DCAz cut to PV Minimum"};
  Configurable<int> trackSelection{"trackSelection", 1, "Track selection: 0 -> No Cut, 1 -> kGlobalTrack, 2 -> kGlobalTrackWoPtEta, 3 -> kGlobalTrackWoDCA, 4 -> kQualityTracks, 5 -> kInAcceptanceTracks"};

  // Configurables for V0s
  Configurable<double> cMinV0Radius{"cMinV0Radius", 0.0, "Minimum V0 radius from PV"};
  Configurable<double> cMaxV0Radius{"cMaxV0Radius", 200.0, "Maximum V0 radius from PV"};
  Configurable<double> cMinV0CosPA{"cMinV0CosPA", 0.995, "Minimum V0 CosPA to PV"};

  // Configurables for cascades
  Configurable<double> cMinCascRadius{"cMinCascRadius", 0.0, "Minimum Cascade radius from PV"};
  Configurable<double> cMaxCascRadius{"cMaxCascRadius", 200.0, "Maximum Cascade radius from PV"};
  Configurable<double> cMinCascCosPA{"cMinCascCosPA", 0.97, "Minimum Cascade CosPA to PV"};

  // Filters
  Filter dcaXYFilter = nabs(aod::track::dcaXY) < cMaxDCArToPVcut && nabs(aod::track::dcaXY) > cMinDCArToPVcut;
  Filter dcaZFilter = nabs(aod::track::dcaZ) < cMaxDCAzToPVcut && nabs(aod::track::dcaZ) > cMinDCAzToPVcut;
  Preslice<aod::McParticles> perMcCollision = aod::mcparticle::mcCollisionId;

  // Track selection filter based on configuration
  Filter trackFilter = (trackSelection.node() == 0) ||
                       ((trackSelection.node() == 1) && requireGlobalTrackInFilter()) ||                                    // kGlobalTrack = kQualityTracks | kPrimaryTracks | kInAcceptanceTracks
                       ((trackSelection.node() == 2) && requireGlobalTrackWoPtEtaInFilter()) ||                             // kGlobalTrackWoPtEta = kQualityTracks | kPrimaryTracks
                       ((trackSelection.node() == 3) && requireGlobalTrackWoDCAInFilter()) ||                               // kGlobalTrackWoDCA = kQualityTracks | kInAcceptanceTracks
                       ((trackSelection.node() == 4) && requireQualityTracksInFilter()) ||                                  // kQualityTracks = kQualityTracksITS | kQualityTracksTPC
                       ((trackSelection.node() == 5) && requireTrackCutInFilter(TrackSelectionFlags::kInAcceptanceTracks)); // kInAcceptanceTracks = kPtRange | kEtaRange
  HistogramRegistry qaRegistry{"QAHistos", {}, OutputObjHandlingPolicy::AnalysisObject};

  /**
   * @brief Initializes the task
   *
   * @param context Initialization context
   */
  void init(InitContext&)
  {
    if (cfgFillQA) {
      AxisSpec idxAxis = {8, -0.5, 7.5, "Index"};
      AxisSpec ptAxis = {100, 0.0f, 10.0f, "#it{p}_{T} (GeV/#it{c})"};
      AxisSpec etaAxis = {100, -1.0f, 1.0f, "#eta"};
      AxisSpec phiAxis = {100, 0.0f, TwoPI, "#phi"};

      qaRegistry.add("QA/hGoodTrackIndices", "hGoodTrackIndices", kTH1D, {idxAxis});
      if (doprocessMC) {
        qaRegistry.add("QA/hGoodMCTrackIndices", "hGoodMCTrackIndices", kTH1D, {idxAxis});
      }
      qaRegistry.add("QA/hTrackPt", "Track pT", kTH1F, {ptAxis});
      qaRegistry.add("QA/hTrackEta", "Track eta", kTH1F, {etaAxis});
      qaRegistry.add("QA/hTrackPhi", "Track phi", kTH1F, {phiAxis});

      if (doprocessV0Data || doprocessV0MC) {
        qaRegistry.add("QA/hGoodV0Indices", "hGoodV0Indices", kTH1D, {idxAxis});
        if (doprocessMC) {
          qaRegistry.add("QA/hGoodMCV0Indices", "hGoodMCV0Indices", kTH1D, {idxAxis});
        }
        AxisSpec radiusAxis = {100, 0.0, 200.0, "V0 Radius"};
        AxisSpec cosPAAxis = {100, 0.995, 1.0, "V0 CosPA"};
        qaRegistry.add("QA/hV0Radius", "V0 Radius", kTH1F, {radiusAxis});
        qaRegistry.add("QA/hV0CosPA", "V0 CosPA", kTH1F, {cosPAAxis});
      }

      if (doprocessCascData || doprocessCascMC) {
        AxisSpec radiusAxis = {100, 0.0, 200.0, "Cascade Radius"};
        AxisSpec cosPAAxis = {100, 0.97, 1.0, "Cascade CosPA"};
        qaRegistry.add("QA/hGoodCascIndices", "hGoodCascIndices", kTH1D, {idxAxis});
        if (doprocessMC) {
          qaRegistry.add("QA/hGoodMCCascIndices", "hGoodMCCascIndices", kTH1D, {idxAxis});
        }
        qaRegistry.add("QA/hCascRadius", "Cascade Radius", kTH1F, {radiusAxis});
        qaRegistry.add("QA/hCascCosPA", "Cascade CosPA", kTH1F, {cosPAAxis});
      }
    }
    if (doprocessData || doprocessMC) {
      LOGF(info, "ResonanceDaughterInitializer initialized with tracks");
    }
    if (doprocessV0Data || doprocessV0MC) {
      LOGF(info, "ResonanceDaughterInitializer initialized with V0s");
    }
    if (doprocessCascData || doprocessCascMC) {
      LOGF(info, "ResonanceDaughterInitializer initialized with cascades");
    }

    // Check if the module is initialized with both data and MC
    if ((doprocessData && doprocessMC) || (doprocessV0Data && doprocessV0MC) || (doprocessCascData && doprocessCascMC)) {
      LOGF(fatal, "ResonanceDaughterInitializer initialized with both data and MC");
    }
    // Check if none of the processes are enabled
    if (!doprocessDummy && !doprocessData && !doprocessMC && !doprocessV0Data && !doprocessV0MC && !doprocessCascData && !doprocessCascMC) {
      LOGF(fatal, "ResonanceDaughterInitializer not initialized, enable at least one process");
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
    // Loop over tracks
    for (auto const& track : tracks) {
      if (cfgFillQA) {
        qaRegistry.fill(HIST("QA/hGoodTrackIndices"), 0);
        qaRegistry.fill(HIST("QA/hTrackPt"), track.pt());
        qaRegistry.fill(HIST("QA/hTrackEta"), track.eta());
        qaRegistry.fill(HIST("QA/hTrackPhi"), track.phi());
      }
      uint8_t trackFlags = (track.passedITSRefit() << 0) |
                           (track.passedTPCRefit() << 1) |
                           (track.isGlobalTrackWoDCA() << 2) |
                           (track.isGlobalTrack() << 3) |
                           (track.isPrimaryTrack() << 4) |
                           (track.isPVContributor() << 5) |
                           (track.hasTOF() << 6) |
                           ((track.sign() > 0) << 7); // sign +1: 1, -1: 0
      reso2trks(collision.globalIndex(),
                track.globalIndex(),
                track.pt(),
                track.px(),
                track.py(),
                track.pz(),
                (uint8_t)track.tpcNClsCrossedRows(),
                (uint8_t)track.tpcNClsFound(),
                static_cast<int16_t>(track.dcaXY() * 10000),
                static_cast<int16_t>(track.dcaZ() * 10000),
                (int8_t)(track.tpcNSigmaPi() * 10),
                (int8_t)(track.tpcNSigmaKa() * 10),
                (int8_t)(track.tpcNSigmaPr() * 10),
                (int8_t)(track.tofNSigmaPi() * 10),
                (int8_t)(track.tofNSigmaKa() * 10),
                (int8_t)(track.tofNSigmaPr() * 10),
                (int8_t)(track.tpcSignal() * 10),
                trackFlags);
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
          if (lDaughter.globalIndex() != 0 && lDaughter.globalIndex() != theMcParticle.globalIndex()) {
            lSiblingsIndeces.push_back(lDaughter.globalIndex());
          }
        }
      }
      return lSiblingsIndeces;
    };
    // ------
    std::vector<int> mothers = {-1, -1};
    std::vector<int> motherPDGs = {-1, -1};
    int siblings[2] = {0, 0};
    std::vector<int> siblingsTemp = {-1, -1};
    if (track.has_mcParticle()) {
      if (cfgFillQA) {
        qaRegistry.fill(HIST("QA/hGoodMCTrackIndices"), 0);
      }
      // Get the MC particle
      const auto& particle = track.mcParticle();
      if (particle.has_mothers()) {
        mothers = getMothersIndeces(particle);
        motherPDGs = getMothersPDGCodes(particle);
        siblingsTemp = getSiblingsIndeces(particle);
      }
      while (mothers.size() > 2) {
        mothers.pop_back();
        motherPDGs.pop_back();
      }
      if (siblingsTemp.size() > 0)
        siblings[0] = siblingsTemp[0];
      if (siblingsTemp.size() > 1)
        siblings[1] = siblingsTemp[1];
      reso2mctracks(particle.pdgCode(),
                    mothers[0],
                    motherPDGs[0],
                    siblings,
                    particle.isPhysicalPrimary(),
                    particle.producedByGenerator());
    } else {
      // No MC particle associated
      reso2mctracks(0,
                    mothers[0],
                    motherPDGs[0],
                    siblings,
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
  void fillV0s(CollisionType const& collision, V0Type const& v0s, TrackType const&)
  {
    int childIDs[2] = {0, 0}; // these IDs are necessary to keep track of the children
    for (auto const& v0 : v0s) {
      if (cfgFillQA) {
        qaRegistry.fill(HIST("QA/hGoodV0Indices"), 0);
        qaRegistry.fill(HIST("QA/hV0Radius"), v0.v0radius());
        qaRegistry.fill(HIST("QA/hV0CosPA"), v0.v0cosPA());
      }
      childIDs[0] = v0.posTrackId();
      childIDs[1] = v0.negTrackId();
      reso2v0s(collision.globalIndex(),
               v0.globalIndex(),
               v0.pt(),
               v0.px(),
               v0.py(),
               v0.pz(),
               childIDs,
               (int8_t)(v0.template posTrack_as<TrackType>().tpcNSigmaPi() * 10),
               (int8_t)(v0.template posTrack_as<TrackType>().tpcNSigmaKa() * 10),
               (int8_t)(v0.template posTrack_as<TrackType>().tpcNSigmaPr() * 10),
               (int8_t)(v0.template negTrack_as<TrackType>().tpcNSigmaPi() * 10),
               (int8_t)(v0.template negTrack_as<TrackType>().tpcNSigmaKa() * 10),
               (int8_t)(v0.template negTrack_as<TrackType>().tpcNSigmaPr() * 10),
               (int8_t)(v0.template negTrack_as<TrackType>().tofNSigmaPi() * 10),
               (int8_t)(v0.template negTrack_as<TrackType>().tofNSigmaKa() * 10),
               (int8_t)(v0.template negTrack_as<TrackType>().tofNSigmaPr() * 10),
               (int8_t)(v0.template posTrack_as<TrackType>().tofNSigmaPi() * 10),
               (int8_t)(v0.template posTrack_as<TrackType>().tofNSigmaKa() * 10),
               (int8_t)(v0.template posTrack_as<TrackType>().tofNSigmaPr() * 10),
               v0.v0cosPA(),
               v0.dcaV0daughters(),
               v0.dcapostopv(),
               v0.dcanegtopv(),
               v0.dcav0topv(),
               v0.mLambda(),
               v0.mAntiLambda(),
               v0.mK0Short(),
               v0.v0radius(), v0.x(), v0.y(), v0.z());
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
    auto getDaughtersIndeces = [&](auto const& theMcParticle) {
      std::vector<int> lDaughtersIndeces{};
      for (auto const& lDaughter : theMcParticle.template daughters_as<aod::McParticles>()) {
        LOGF(debug, "   daughter index lDaughter: %d", lDaughter.globalIndex());
        if (lDaughter.globalIndex() != 0) {
          lDaughtersIndeces.push_back(lDaughter.globalIndex());
        }
      }
      return lDaughtersIndeces;
    };
    auto getDaughtersPDGCodes = [&](auto const& theMcParticle) {
      std::vector<int> lDaughtersPDGs{};
      for (auto const& lDaughter : theMcParticle.template daughters_as<aod::McParticles>()) {
        LOGF(debug, "   daughter pdgcode lDaughter: %d", lDaughter.pdgCode());
        if (lDaughter.globalIndex() != 0) {
          lDaughtersPDGs.push_back(lDaughter.pdgCode());
        }
      }
      return lDaughtersPDGs;
    };
    // ------
    std::vector<int> mothers = {-1, -1};
    std::vector<int> motherPDGs = {-1, -1};
    std::vector<int> daughters = {-1, -1};
    std::vector<int> daughterPDGs = {-1, -1};
    if (v0.has_mcParticle()) {
      if (cfgFillQA) {
        qaRegistry.fill(HIST("QA/hGoodMCV0Indices"), 0);
      }
      auto v0mc = v0.mcParticle();
      if (v0mc.has_mothers()) {
        mothers = getMothersIndeces(v0mc);
        motherPDGs = getMothersPDGCodes(v0mc);
      }
      while (mothers.size() > 2) {
        mothers.pop_back();
        motherPDGs.pop_back();
      }
      if (v0mc.has_daughters()) {
        daughters = getDaughtersIndeces(v0mc);
        daughterPDGs = getDaughtersPDGCodes(v0mc);
      }
      while (daughters.size() > 2) {
        LOGF(info, "daughters.size() is larger than 2");
        daughters.pop_back();
        daughterPDGs.pop_back();
      }
      reso2mcv0s(v0mc.pdgCode(),
                 mothers[0],
                 motherPDGs[0],
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
  void fillCascades(CollisionType const& collision, CascType const& cascades, TrackType const&)
  {
    int childIDs[3] = {0, 0, 0}; // these IDs are necessary to keep track of the children
    for (auto const& casc : cascades) {
      if (cfgFillQA) {
        qaRegistry.fill(HIST("QA/hGoodCascIndices"), 0);
        qaRegistry.fill(HIST("QA/hCascRadius"), casc.cascradius());
        qaRegistry.fill(HIST("QA/hCascCosPA"), casc.casccosPA(collision.posX(), collision.posY(), collision.posZ()));
      }
      childIDs[0] = casc.posTrackId();
      childIDs[1] = casc.negTrackId();
      childIDs[2] = casc.bachelorId();
      reso2cascades(collision.globalIndex(),
                    casc.globalIndex(),
                    casc.pt(),
                    casc.px(),
                    casc.py(),
                    casc.pz(),
                    childIDs,
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
                    casc.mLambda(),
                    casc.mXi(),
                    casc.v0radius(), casc.cascradius(), casc.x(), casc.y(), casc.z());
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
    auto getDaughtersIndeces = [&](auto const& theMcParticle) {
      std::vector<int> lDaughtersIndeces{};
      for (auto const& lDaughter : theMcParticle.template daughters_as<aod::McParticles>()) {
        LOGF(debug, "   daughter index lDaughter: %d", lDaughter.globalIndex());
        if (lDaughter.globalIndex() != 0) {
          lDaughtersIndeces.push_back(lDaughter.globalIndex());
        }
      }
      return lDaughtersIndeces;
    };
    auto getDaughtersPDGCodes = [&](auto const& theMcParticle) {
      std::vector<int> lDaughtersPDGs{};
      for (auto const& lDaughter : theMcParticle.template daughters_as<aod::McParticles>()) {
        LOGF(debug, "   daughter pdgcode lDaughter: %d", lDaughter.pdgCode());
        if (lDaughter.globalIndex() != 0) {
          lDaughtersPDGs.push_back(lDaughter.pdgCode());
        }
      }
      return lDaughtersPDGs;
    };
    // ------
    std::vector<int> mothers = {-1, -1};
    std::vector<int> motherPDGs = {-1, -1};
    std::vector<int> daughters = {-1, -1};
    std::vector<int> daughterPDGs = {-1, -1};
    if (casc.has_mcParticle()) {
      if (cfgFillQA) {
        qaRegistry.fill(HIST("QA/hGoodMCCascIndices"), 0);
      }
      auto cascmc = casc.mcParticle();
      if (cascmc.has_mothers()) {
        mothers = getMothersIndeces(cascmc);
        motherPDGs = getMothersPDGCodes(cascmc);
      }
      while (mothers.size() > 2) {
        mothers.pop_back();
        motherPDGs.pop_back();
      }
      if (cascmc.has_daughters()) {
        daughters = getDaughtersIndeces(cascmc);
        daughterPDGs = getDaughtersPDGCodes(cascmc);
      }
      while (daughters.size() > 2) {
        LOGF(info, "daughters.size() is larger than 2");
        daughters.pop_back();
        daughterPDGs.pop_back();
      }
      reso2mccascades(cascmc.pdgCode(),
                      mothers[0],
                      motherPDGs[0],
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
   * @brief Processes data tracks
   *
   * @param collision Collision data
   * @param tracks Track data
   */
  void processData(aod::ResoCollision const& collision,
                   soa::Filtered<aod::ResoTrackCandidates> const& tracks)
  {
    fillTracks<false>(collision, tracks);
  }
  PROCESS_SWITCH(ResonanceDaughterInitializer, processData, "Process tracks for data", false);

  /**
   * @brief Processes MC tracks
   *
   * @param collision Collision data
   * @param tracks Track data
   * @param mcParticles MC particles
   */
  void processMC(aod::ResoCollision const& collision,
                 soa::Filtered<aod::ResoTrackCandidatesMC> const& tracks,
                 aod::McParticles const&)
  {
    fillTracks<true>(collision, tracks);
  }
  PROCESS_SWITCH(ResonanceDaughterInitializer, processMC, "Process tracks for MC", false);

  /**
   * @brief Processes V0 data
   *
   * @param collision Collision data
   * @param v0s V0 data
   * @param tracks Track data
   */
  void processV0Data(aod::ResoCollision const& collision, aod::ResoV0Candidates const& v0s, aod::ResoTrackCandidates const& tracks)
  {
    fillV0s<false>(collision, v0s, tracks);
  }
  PROCESS_SWITCH(ResonanceDaughterInitializer, processV0Data, "Process V0s for data", false);

  /**
   * @brief Processes MC V0 data
   *
   * @param collision Collision data
   * @param v0s V0 data
   * @param tracks Track data
   */
  void processV0MC(aod::ResoCollision const& collision, aod::ResoV0CandidatesMC const& v0s, aod::ResoTrackCandidatesMC const& tracks)
  {
    fillV0s<true>(collision, v0s, tracks);
  }
  PROCESS_SWITCH(ResonanceDaughterInitializer, processV0MC, "Process V0s for MC", false);

  /**
   * @brief Processes cascade data
   *
   * @param collision Collision data
   * @param cascades Cascade data
   * @param tracks Track data
   */
  void processCascData(aod::ResoCollision const& collision, aod::ResoCascadesCandidates const& cascades, aod::ResoTrackCandidates const& tracks)
  {
    fillCascades<false>(collision, cascades, tracks);
  }
  PROCESS_SWITCH(ResonanceDaughterInitializer, processCascData, "Process Cascades for data", false);

  /**
   * @brief Processes MC cascade data
   *
   * @param collision Collision data
   * @param cascades Cascade data
   * @param tracks Track data
   */
  void processCascMC(aod::ResoCollision const& collision, aod::ResoCascadesCandidatesMC const& cascades, aod::ResoTrackCandidatesMC const& tracks)
  {
    fillCascades<true>(collision, cascades, tracks);
  }
  PROCESS_SWITCH(ResonanceDaughterInitializer, processCascMC, "Process Cascades for MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<ResonanceModuleInitializer>(cfgc),
    adaptAnalysisTask<ResonanceDaughterInitializer>(cfgc)};
}
