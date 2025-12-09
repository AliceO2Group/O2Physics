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

/// \file dataCreatorCharmResoToDplusReduced.cxx
/// \brief Creation of D+ V0 and D+ track pairs
///
/// \author Luca Aglietta <luca.aglietta@cern.ch>, UniTO Turin
/// \author Fabrizio Grosa <fabrizio.grosa@cern.ch>, CERN

#include "PWGHF/Core/CentralityEstimation.h"
#include "PWGHF/D2H/Core/DataCreationCharmReso.h"
#include "PWGHF/D2H/DataModel/ReducedDataModel.h"
#include "PWGHF/D2H/Utils/utilsRedDataFormat.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/Utils/utilsBfieldCCDB.h"
#include "PWGHF/Utils/utilsEvSelHf.h"

#include "Common/Core/ZorroSummary.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CCDB/CcdbApi.h>
#include <DCAFitter/DCAFitterN.h>
#include <DetectorsBase/MatLayerCylSet.h>
#include <DetectorsBase/Propagator.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/DeviceSpec.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/RunningWorkflowInfo.h>
#include <Framework/runDataProcessing.h>

#include <chrono>
#include <cmath>
#include <string>

using namespace o2;
using namespace o2::analysis;
using namespace o2::aod;
using namespace o2::analysis::hf_charm_reso;
using namespace o2::constants::physics;
using namespace o2::framework;
using namespace o2::framework::expressions;

/// Creation of D-V0 pairs
struct HfDataCreatorCharmResoToDplusReduced {

  // Produces AOD tables to store collision information
  Produces<aod::HfRedCollisions> hfReducedCollision; // Defined in PWGHF/D2H/DataModel/ReducedDataModel.h
  Produces<aod::HfOrigColCounts> hfCollisionCounter; // Defined in PWGHF/D2H/DataModel/ReducedDataModel.h
  // tracks, V0 and D candidates reduced tables
  Produces<aod::HfRedVzeros> hfCandV0;            // Defined in PWGHF/D2H/DataModel/ReducedDataModel.h
  Produces<aod::HfRedTrkNoParams> hfTrackNoParam; // Defined in PWGHF/D2H/DataModel/ReducedDataModel.h
  Produces<aod::HfRed3PrNoTrks> hfCandD3Pr;       // Defined in PWGHF/D2H/DataModel/ReducedDataModel.h
  // ML optional Tables
  Produces<aod::HfRed3ProngsMl> hfCandD3PrMl; // Defined in PWGHF/D2H/DataModel/ReducedDataModel.h
  // MC Tables
  Produces<aod::HfMcGenRedResos> rowHfResoMcGenReduced;
  Produces<aod::Hf3PrV0McRec> rowHf3PrV0McRecReduced;
  Produces<aod::Hf3PrTrkMcRec> rowHf3PrTrkMcRecReduced;

  // selection D
  struct : ConfigurableGroup {
    std::string prefix = "dmesons";
    Configurable<int> selectionFlagDplus{"selectionFlagDplus", 7, "Selection Flag for D"};
  } cfgDmesCuts;

  // selection V0
  HfResoConfigV0Cuts cfgV0Cuts;
  // selection single tracks
  HfResoConfigSingleTrackCuts cfgSingleTrackCuts;
  // QA histograms
  HfResoConfigQaPlots cfgQaPlots;
  // other configurables
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> ccdbPathGrpMag{"ccdbPathGrpMag", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object (Run 3)"};
  Configurable<bool> doMcRecQa{"doMcRecQa", true, "Fill QA histograms for Mc matching"};
  Configurable<bool> rejectPairsWithCommonDaughter{"rejectPairsWithCommonDaughter", true, "flag to reject already at this stage the pairs that share a daughter track"};
  Configurable<bool> rejectCollisionsWithBadEvSel{"rejectCollisionsWithBadEvSel", true, "flag to reject collisions with bad event selection"};

  o2::hf_evsel::HfEventSelection hfEvSel;
  o2::hf_evsel::HfEventSelectionMc hfEvSelMc;

  // CCDB service
  o2::ccdb::CcdbApi ccdbApi;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  double bz{0.};
  int runNumber{0}; // needed to detect if the run changed and trigger update of calibrations etc.

  // material correction for track propagation
  o2::base::MatLayerCylSet* lut{};
  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;

  // O2DatabasePDG service
  Service<o2::framework::O2DatabasePDG> pdg;

  // vertex fitter
  o2::vertexing::DCAFitterN<2> fitter;

  // Dplus
  using CandsDplusFiltered = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDplusToPiKPi>>;
  using CandsDplusFilteredWithMl = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDplusToPiKPi, aod::HfMlDplusToPiKPi>>;
  using CandsDplusFilteredWithMc = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDplusToPiKPi, aod::HfCand3ProngMcRec>>;
  using CandsDplusFilteredWithMlAndMc = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDplusToPiKPi, aod::HfMlDplusToPiKPi, aod::HfCand3ProngMcRec>>;
  // Tracks
  using TracksWithPID = soa::Join<aod::Tracks, aod::TracksCov, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullKa, aod::pidTOFFullKa, aod::pidTPCFullPr, aod::pidTOFFullPr>;
  using TracksWithPIDAndMC = soa::Join<TracksWithPID, aod::McTrackLabels>;
  using TracksIUWithPID = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::pidTPCPi, aod::pidTOFFullPi, aod::pidTPCPr, aod::pidTOFFullPr>;
  using TracksIUWithPIDAndMC = soa::Join<TracksIUWithPID, aod::McTrackLabels>;
  // Collisions MC
  using BCsInfo = soa::Join<aod::BCs, aod::Timestamps, aod::BcSels>;
  using McCollisionsNoCents = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>;

  Filter filterSelectDplus = (aod::hf_sel_candidate_dplus::isSelDplusToPiKPi >= cfgDmesCuts.selectionFlagDplus);

  Preslice<CandsDplusFiltered> candsDplusPerCollision = aod::hf_cand::collisionId;
  Preslice<CandsDplusFilteredWithMl> candsDplusPerCollisionWithMl = aod::hf_cand::collisionId;
  Preslice<aod::V0s> candsV0PerCollision = aod::v0::collisionId;
  Preslice<aod::TrackAssoc> trackIndicesPerCollision = aod::track_association::collisionId;
  Preslice<aod::McParticles> mcParticlesPerMcCollision = aod::mcparticle::mcCollisionId;
  PresliceUnsorted<McCollisionsNoCents> colPerMcCollision = aod::mccollisionlabel::mcCollisionId;

  HistogramRegistry registry{"registry"};
  OutputObj<ZorroSummary> zorroSummary{"zorroSummary"};

  void init(InitContext& initContext)
  {
    // histograms
    if (doprocessDplusV0MC || doprocessDplusTrackMC || doprocessDplusV0AndTrackMC || doprocessDplusV0MCWithMl || doprocessDplusTrackMCWithMl || doprocessDplusV0AndTrackMCWithMl) {
      addHistograms<true, DMesonType::Dplus>(registry);
    } else {
      addHistograms<false, DMesonType::Dplus>(registry);
    }

    // Configure CCDB access
    ccdb->setURL(ccdbUrl.value);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setCreatedNotAfter(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
    ccdbApi.init(ccdbUrl);
    runNumber = 0;
    lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>("GLO/Param/MatLUT"));

    // Configure DCA fitter
    fitter.setPropagateToPCA(true);
    fitter.setMaxR(200.);
    fitter.setMinParamChange(1e-3);
    fitter.setMinRelChi2Change(0.9);
    fitter.setMaxDZIni(1e9);
    fitter.setMaxDXYIni(4);
    fitter.setMaxChi2(1e9);
    fitter.setUseAbsDCA(true);
    fitter.setWeightedFinalPCA(false);

    // init HF event selection helper
    hfEvSel.init(registry, zorroSummary);

    const auto& workflows = initContext.services().get<RunningWorkflowInfo const>();
    for (const DeviceSpec& device : workflows.devices) {
      if (device.name == "hf-data-creator-charm-reso-to-dplus-reduced") {
        // init HF event selection helper
        hfEvSelMc.init(device, registry);
        break;
      }
    }
  }

  // Process functions
  // No ML
  // Data
  void processDplusV0(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                      CandsDplusFiltered const& candsDplus,
                      aod::V0s const& v0s,
                      TracksIUWithPID const& tracksIU,
                      aod::BCsWithTimestamps const&)
  {
    int zvtxColl{0};
    int sel8Coll{0};
    int zvtxAndSel8Coll{0};
    int zvtxAndSel8CollAndSoftTrig{0};
    int allSelColl{0};
    int dummyTable{0};
    for (const auto& collision : collisions) {
      const auto hfRejMap = o2::hf_evsel::getEvSel<true, o2::hf_centrality::CentralityEstimator::None, aod::BCsWithTimestamps>(collision, hfEvSel, zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl, ccdb, registry);
      if (rejectCollisionsWithBadEvSel && hfRejMap != 0) {
        continue;
      }
      auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
      if (runNumber != bc.runNumber()) {
        LOG(info) << ">>>>>>>>>>>> Current run number: " << runNumber;
        initCCDB(bc, runNumber, ccdb, ccdbPathGrpMag, lut, false);
        bz = o2::base::Propagator::Instance()->getNominalBz();
        LOG(info) << ">>>>>>>>>>>> Magnetic field: " << bz;
      }
      fitter.setBz(bz);
      auto thisCollId = collision.globalIndex();
      auto candsDThisColl = candsDplus.sliceBy(candsDplusPerCollision, thisCollId);
      auto v0sThisColl = v0s.sliceBy(candsV0PerCollision, thisCollId);
      runDataCreation<false, false, DMesonType::Dplus, PairingType::V0Only>(collision, candsDThisColl, v0sThisColl, nullptr, tracksIU, tracksIU, nullptr, hfRejMap, bz, pdg, registry, matCorr, fitter, cfgDmesCuts, cfgSingleTrackCuts, cfgV0Cuts, cfgQaPlots, rejectPairsWithCommonDaughter, hfReducedCollision, hfCandD3Pr, hfCandV0, dummyTable, dummyTable, dummyTable, dummyTable);
    }
    // handle normalization by the right number of collisions
    hfCollisionCounter(collisions.tableSize(), zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
  }
  PROCESS_SWITCH(HfDataCreatorCharmResoToDplusReduced, processDplusV0, "Process Dplus candidates paired with V0s", false);

  void processDplusTrack(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                         CandsDplusFiltered const& candsDplus,
                         aod::TrackAssoc const& trackIndices,
                         TracksWithPID const& tracks,
                         aod::BCsWithTimestamps const&)
  {
    int zvtxColl{0};
    int sel8Coll{0};
    int zvtxAndSel8Coll{0};
    int zvtxAndSel8CollAndSoftTrig{0};
    int allSelColl{0};
    int dummyTable{0};
    for (const auto& collision : collisions) {
      const auto hfRejMap = o2::hf_evsel::getEvSel<true, o2::hf_centrality::CentralityEstimator::None, aod::BCsWithTimestamps>(collision, hfEvSel, zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl, ccdb, registry);
      if (rejectCollisionsWithBadEvSel && hfRejMap != 0) {
        continue;
      }
      auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
      if (runNumber != bc.runNumber()) {
        LOG(info) << ">>>>>>>>>>>> Current run number: " << runNumber;
        initCCDB(bc, runNumber, ccdb, ccdbPathGrpMag, lut, false);
        bz = o2::base::Propagator::Instance()->getNominalBz();
        LOG(info) << ">>>>>>>>>>>> Magnetic field: " << bz;
      }
      fitter.setBz(bz);
      auto thisCollId = collision.globalIndex();
      auto candsDThisColl = candsDplus.sliceBy(candsDplusPerCollision, thisCollId);
      auto trackIdsThisColl = trackIndices.sliceBy(trackIndicesPerCollision, thisCollId);
      runDataCreation<false, false, DMesonType::Dplus, PairingType::TrackOnly>(collision, candsDThisColl, nullptr, trackIdsThisColl, tracks, tracks, nullptr, hfRejMap, bz, pdg, registry, matCorr, fitter, cfgDmesCuts, cfgSingleTrackCuts, cfgV0Cuts, cfgQaPlots, rejectPairsWithCommonDaughter, hfReducedCollision, hfCandD3Pr, dummyTable, hfTrackNoParam, dummyTable, dummyTable, dummyTable);
    }
    // handle normalization by the right number of collisions
    hfCollisionCounter(collisions.tableSize(), zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
  }
  PROCESS_SWITCH(HfDataCreatorCharmResoToDplusReduced, processDplusTrack, "Process Dplus candidates paired with Tracks", false);

  void processDplusV0AndTrack(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                              CandsDplusFiltered const& candsDplus,
                              aod::V0s const& v0s,
                              aod::TrackAssoc const& trackIndices,
                              TracksWithPID const& tracks,
                              TracksIUWithPID const& tracksIU,
                              aod::BCsWithTimestamps const&)
  {
    int zvtxColl{0};
    int sel8Coll{0};
    int zvtxAndSel8Coll{0};
    int zvtxAndSel8CollAndSoftTrig{0};
    int allSelColl{0};
    int dummyTable{0};
    for (const auto& collision : collisions) {
      const auto hfRejMap = o2::hf_evsel::getEvSel<true, o2::hf_centrality::CentralityEstimator::None, aod::BCsWithTimestamps>(collision, hfEvSel, zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl, ccdb, registry);
      if (rejectCollisionsWithBadEvSel && hfRejMap != 0) {
        continue;
      }
      auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
      if (runNumber != bc.runNumber()) {
        LOG(info) << ">>>>>>>>>>>> Current run number: " << runNumber;
        initCCDB(bc, runNumber, ccdb, ccdbPathGrpMag, lut, false);
        bz = o2::base::Propagator::Instance()->getNominalBz();
        LOG(info) << ">>>>>>>>>>>> Magnetic field: " << bz;
      }
      fitter.setBz(bz);
      auto thisCollId = collision.globalIndex();
      auto candsDThisColl = candsDplus.sliceBy(candsDplusPerCollision, thisCollId);
      auto v0sThisColl = v0s.sliceBy(candsV0PerCollision, thisCollId);
      auto trackIdsThisColl = trackIndices.sliceBy(trackIndicesPerCollision, thisCollId);
      runDataCreation<false, false, DMesonType::Dplus, PairingType::V0AndTrack>(collision, candsDThisColl, v0sThisColl, trackIdsThisColl, tracks, tracksIU, nullptr, hfRejMap, bz, pdg, registry, matCorr, fitter, cfgDmesCuts, cfgSingleTrackCuts, cfgV0Cuts, cfgQaPlots, rejectPairsWithCommonDaughter, hfReducedCollision, hfCandD3Pr, hfCandV0, hfTrackNoParam, dummyTable, dummyTable, dummyTable);
    }
    // handle normalization by the right number of collisions
    hfCollisionCounter(collisions.tableSize(), zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
  }
  PROCESS_SWITCH(HfDataCreatorCharmResoToDplusReduced, processDplusV0AndTrack, "Process Dplus candidates paired with V0s and Tracks", false);

  // ML
  // Data
  void processDplusV0WithMl(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                            CandsDplusFilteredWithMl const& candsDplus,
                            aod::V0s const& v0s,
                            TracksIUWithPID const& tracksIU,
                            aod::BCsWithTimestamps const&)
  {
    int zvtxColl{0};
    int sel8Coll{0};
    int zvtxAndSel8Coll{0};
    int zvtxAndSel8CollAndSoftTrig{0};
    int allSelColl{0};
    int dummyTable{0};
    for (const auto& collision : collisions) {
      const auto hfRejMap = o2::hf_evsel::getEvSel<true, o2::hf_centrality::CentralityEstimator::None, aod::BCsWithTimestamps>(collision, hfEvSel, zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl, ccdb, registry);
      if (rejectCollisionsWithBadEvSel && hfRejMap != 0) {
        continue;
      }
      auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
      if (runNumber != bc.runNumber()) {
        LOG(info) << ">>>>>>>>>>>> Current run number: " << runNumber;
        initCCDB(bc, runNumber, ccdb, ccdbPathGrpMag, lut, false);
        bz = o2::base::Propagator::Instance()->getNominalBz();
        LOG(info) << ">>>>>>>>>>>> Magnetic field: " << bz;
      }
      fitter.setBz(bz);
      auto thisCollId = collision.globalIndex();
      auto candsDThisColl = candsDplus.sliceBy(candsDplusPerCollisionWithMl, thisCollId);
      auto v0sThisColl = v0s.sliceBy(candsV0PerCollision, thisCollId);
      runDataCreation<true, false, DMesonType::Dplus, PairingType::V0Only>(collision, candsDThisColl, v0sThisColl, nullptr, tracksIU, tracksIU, nullptr, hfRejMap, bz, pdg, registry, matCorr, fitter, cfgDmesCuts, cfgSingleTrackCuts, cfgV0Cuts, cfgQaPlots, rejectPairsWithCommonDaughter, hfReducedCollision, hfCandD3Pr, hfCandV0, dummyTable, dummyTable, dummyTable, hfCandD3PrMl);
    }
    // handle normalization by the right number of collisions
    hfCollisionCounter(collisions.tableSize(), zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
  }
  PROCESS_SWITCH(HfDataCreatorCharmResoToDplusReduced, processDplusV0WithMl, "Process Dplus candidates paired with V0s with ML info", false);

  void processDplusTrackWithMl(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                               CandsDplusFilteredWithMl const& candsDplus,
                               aod::TrackAssoc const& trackIndices,
                               TracksWithPID const& tracks,
                               aod::BCsWithTimestamps const&)
  {
    int zvtxColl{0};
    int sel8Coll{0};
    int zvtxAndSel8Coll{0};
    int zvtxAndSel8CollAndSoftTrig{0};
    int allSelColl{0};
    int dummyTable{0};
    for (const auto& collision : collisions) {
      const auto hfRejMap = o2::hf_evsel::getEvSel<true, o2::hf_centrality::CentralityEstimator::None, aod::BCsWithTimestamps>(collision, hfEvSel, zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl, ccdb, registry);
      if (rejectCollisionsWithBadEvSel && hfRejMap != 0) {
        continue;
      }
      auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
      if (runNumber != bc.runNumber()) {
        LOG(info) << ">>>>>>>>>>>> Current run number: " << runNumber;
        initCCDB(bc, runNumber, ccdb, ccdbPathGrpMag, lut, false);
        bz = o2::base::Propagator::Instance()->getNominalBz();
        LOG(info) << ">>>>>>>>>>>> Magnetic field: " << bz;
      }
      fitter.setBz(bz);
      auto thisCollId = collision.globalIndex();
      auto candsDThisColl = candsDplus.sliceBy(candsDplusPerCollisionWithMl, thisCollId);
      auto trackIdsThisColl = trackIndices.sliceBy(trackIndicesPerCollision, thisCollId);
      runDataCreation<true, false, DMesonType::Dplus, PairingType::TrackOnly>(collision, candsDThisColl, nullptr, trackIdsThisColl, tracks, tracks, nullptr, hfRejMap, bz, pdg, registry, matCorr, fitter, cfgDmesCuts, cfgSingleTrackCuts, cfgV0Cuts, cfgQaPlots, rejectPairsWithCommonDaughter, hfReducedCollision, hfCandD3Pr, dummyTable, hfTrackNoParam, dummyTable, dummyTable, hfCandD3PrMl);
    }
    // handle normalization by the right number of collisions
    hfCollisionCounter(collisions.tableSize(), zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
  }
  PROCESS_SWITCH(HfDataCreatorCharmResoToDplusReduced, processDplusTrackWithMl, "Process Dplus candidates paired with Tracks with ML info", false);

  void processDplusV0AndTrackWithMl(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                                    CandsDplusFilteredWithMl const& candsDplus,
                                    aod::V0s const& v0s,
                                    aod::TrackAssoc const& trackIndices,
                                    TracksWithPID const& tracks,
                                    TracksIUWithPID const& tracksIU,
                                    aod::BCsWithTimestamps const&)
  {
    int zvtxColl{0};
    int sel8Coll{0};
    int zvtxAndSel8Coll{0};
    int zvtxAndSel8CollAndSoftTrig{0};
    int allSelColl{0};
    int dummyTable{0};
    for (const auto& collision : collisions) {
      const auto hfRejMap = o2::hf_evsel::getEvSel<true, o2::hf_centrality::CentralityEstimator::None, aod::BCsWithTimestamps>(collision, hfEvSel, zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl, ccdb, registry);
      if (rejectCollisionsWithBadEvSel && hfRejMap != 0) {
        continue;
      }
      auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
      if (runNumber != bc.runNumber()) {
        LOG(info) << ">>>>>>>>>>>> Current run number: " << runNumber;
        initCCDB(bc, runNumber, ccdb, ccdbPathGrpMag, lut, false);
        bz = o2::base::Propagator::Instance()->getNominalBz();
        LOG(info) << ">>>>>>>>>>>> Magnetic field: " << bz;
      }
      fitter.setBz(bz);
      auto thisCollId = collision.globalIndex();
      auto candsDThisColl = candsDplus.sliceBy(candsDplusPerCollisionWithMl, thisCollId);
      auto v0sThisColl = v0s.sliceBy(candsV0PerCollision, thisCollId);
      auto trackIdsThisColl = trackIndices.sliceBy(trackIndicesPerCollision, thisCollId);
      runDataCreation<true, false, DMesonType::Dplus, PairingType::V0AndTrack>(collision, candsDThisColl, v0sThisColl, trackIdsThisColl, tracks, tracksIU, nullptr, hfRejMap, bz, pdg, registry, matCorr, fitter, cfgDmesCuts, cfgSingleTrackCuts, cfgV0Cuts, cfgQaPlots, rejectPairsWithCommonDaughter, hfReducedCollision, hfCandD3Pr, hfCandV0, hfTrackNoParam, dummyTable, dummyTable, hfCandD3PrMl);
    }
    // handle normalization by the right number of collisions
    hfCollisionCounter(collisions.tableSize(), zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
  }
  PROCESS_SWITCH(HfDataCreatorCharmResoToDplusReduced, processDplusV0AndTrackWithMl, "Process Dplus candidates paired with V0s and Tracks with ML info", false);

  // Dplus
  void processDplusV0MC(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                        CandsDplusFilteredWithMc const& candsDplus,
                        aod::V0s const& v0s,
                        TracksIUWithPIDAndMC const& tracksIU,
                        aod::McParticles const& particlesMc,
                        BCsInfo const& bcs,
                        McCollisionsNoCents const& collInfos,
                        aod::McCollisions const& mcCollisions)
  {
    int zvtxColl{0};
    int sel8Coll{0};
    int zvtxAndSel8Coll{0};
    int zvtxAndSel8CollAndSoftTrig{0};
    int allSelColl{0};
    int dummyTable{0};
    for (const auto& collision : collisions) {
      const auto hfRejMap = o2::hf_evsel::getEvSel<true, o2::hf_centrality::CentralityEstimator::None, BCsInfo>(collision, hfEvSel, zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl, ccdb, registry);
      if (rejectCollisionsWithBadEvSel && hfRejMap != 0) {
        continue;
      }
      auto bc = collision.template bc_as<BCsInfo>();
      if (runNumber != bc.runNumber()) {
        LOG(info) << ">>>>>>>>>>>> Current run number: " << runNumber;
        initCCDB(bc, runNumber, ccdb, ccdbPathGrpMag, lut, false);
        bz = o2::base::Propagator::Instance()->getNominalBz();
        LOG(info) << ">>>>>>>>>>>> Magnetic field: " << bz;
      }
      fitter.setBz(bz);
      auto thisCollId = collision.globalIndex();
      auto candsDThisColl = candsDplus.sliceBy(candsDplusPerCollision, thisCollId);
      auto v0sThisColl = v0s.sliceBy(candsV0PerCollision, thisCollId);
      runDataCreation<false, true, DMesonType::Dplus, PairingType::V0Only>(collision, candsDThisColl, v0sThisColl, nullptr, tracksIU, tracksIU, particlesMc, hfRejMap, bz, pdg, registry, matCorr, fitter, cfgDmesCuts, cfgSingleTrackCuts, cfgV0Cuts, cfgQaPlots, rejectPairsWithCommonDaughter, hfReducedCollision, hfCandD3Pr, hfCandV0, dummyTable, rowHf3PrV0McRecReduced, dummyTable, dummyTable);
    }
    runMcGen<DMesonType::Dplus, PairingType::V0Only>(particlesMc, mcParticlesPerMcCollision, collInfos, colPerMcCollision, mcCollisions, hfEvSelMc, rejectCollisionsWithBadEvSel, registry, pdg, rowHfResoMcGenReduced, bcs);
    // handle normalization by the right number of collisions
    hfCollisionCounter(collisions.tableSize(), zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
  }
  PROCESS_SWITCH(HfDataCreatorCharmResoToDplusReduced, processDplusV0MC, "Process Dplus candidates paired with V0s with MC matching", false);

  void processDplusTrackMC(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                           CandsDplusFilteredWithMc const& candsDplus,
                           TracksWithPIDAndMC const& tracks,
                           aod::TrackAssoc const& trackIndices,
                           aod::McParticles const& particlesMc,
                           BCsInfo const& bcs,
                           McCollisionsNoCents const& collInfos,
                           aod::McCollisions const& mcCollisions)
  {
    int zvtxColl{0};
    int sel8Coll{0};
    int zvtxAndSel8Coll{0};
    int zvtxAndSel8CollAndSoftTrig{0};
    int allSelColl{0};
    int dummyTable{0};
    for (const auto& collision : collisions) {
      const auto hfRejMap = o2::hf_evsel::getEvSel<true, o2::hf_centrality::CentralityEstimator::None, BCsInfo>(collision, hfEvSel, zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl, ccdb, registry);
      if (rejectCollisionsWithBadEvSel && hfRejMap != 0) {
        continue;
      }
      auto bc = collision.template bc_as<BCsInfo>();
      if (runNumber != bc.runNumber()) {
        LOG(info) << ">>>>>>>>>>>> Current run number: " << runNumber;
        initCCDB(bc, runNumber, ccdb, ccdbPathGrpMag, lut, false);
        bz = o2::base::Propagator::Instance()->getNominalBz();
        LOG(info) << ">>>>>>>>>>>> Magnetic field: " << bz;
      }
      fitter.setBz(bz);
      auto thisCollId = collision.globalIndex();
      auto candsDThisColl = candsDplus.sliceBy(candsDplusPerCollision, thisCollId);
      auto trackIdsThisColl = trackIndices.sliceBy(trackIndicesPerCollision, thisCollId);
      runDataCreation<false, true, DMesonType::Dplus, PairingType::TrackOnly>(collision, candsDThisColl, nullptr, trackIdsThisColl, tracks, tracks, particlesMc, hfRejMap, bz, pdg, registry, matCorr, fitter, cfgDmesCuts, cfgSingleTrackCuts, cfgV0Cuts, cfgQaPlots, rejectPairsWithCommonDaughter, hfReducedCollision, hfCandD3Pr, dummyTable, hfTrackNoParam, dummyTable, rowHf3PrTrkMcRecReduced, dummyTable);
    }
    runMcGen<DMesonType::Dplus, PairingType::TrackOnly>(particlesMc, mcParticlesPerMcCollision, collInfos, colPerMcCollision, mcCollisions, hfEvSelMc, rejectCollisionsWithBadEvSel, registry, pdg, rowHfResoMcGenReduced, bcs);
    // handle normalization by the right number of collisions
    hfCollisionCounter(collisions.tableSize(), zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
  }
  PROCESS_SWITCH(HfDataCreatorCharmResoToDplusReduced, processDplusTrackMC, "Process Dplus candidates paired with tracks with MC matching", false);

  void processDplusV0AndTrackMC(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                                CandsDplusFilteredWithMc const& candsDplus,
                                aod::V0s const& v0s,
                                aod::TrackAssoc const& trackIndices,
                                TracksWithPIDAndMC const& tracks,
                                TracksIUWithPIDAndMC const& tracksIU,
                                aod::McParticles const& particlesMc,
                                BCsInfo const& bcs,
                                McCollisionsNoCents const& collInfos,
                                aod::McCollisions const& mcCollisions)
  {
    int zvtxColl{0};
    int sel8Coll{0};
    int zvtxAndSel8Coll{0};
    int zvtxAndSel8CollAndSoftTrig{0};
    int allSelColl{0};
    int dummyTable{0};
    for (const auto& collision : collisions) {
      const auto hfRejMap = o2::hf_evsel::getEvSel<true, o2::hf_centrality::CentralityEstimator::None, BCsInfo>(collision, hfEvSel, zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl, ccdb, registry);
      if (rejectCollisionsWithBadEvSel && hfRejMap != 0) {
        continue;
      }
      auto bc = collision.template bc_as<BCsInfo>();
      if (runNumber != bc.runNumber()) {
        LOG(info) << ">>>>>>>>>>>> Current run number: " << runNumber;
        initCCDB(bc, runNumber, ccdb, ccdbPathGrpMag, lut, false);
        bz = o2::base::Propagator::Instance()->getNominalBz();
        LOG(info) << ">>>>>>>>>>>> Magnetic field: " << bz;
      }
      fitter.setBz(bz);
      auto thisCollId = collision.globalIndex();
      auto candsDThisColl = candsDplus.sliceBy(candsDplusPerCollision, thisCollId);
      auto v0sThisColl = v0s.sliceBy(candsV0PerCollision, thisCollId);
      auto trackIdsThisColl = trackIndices.sliceBy(trackIndicesPerCollision, thisCollId);
      runDataCreation<false, true, DMesonType::Dplus, PairingType::TrackOnly>(collision, candsDThisColl, v0sThisColl, trackIdsThisColl, tracks, tracksIU, particlesMc, hfRejMap, bz, pdg, registry, matCorr, fitter, cfgDmesCuts, cfgSingleTrackCuts, cfgV0Cuts, cfgQaPlots, rejectPairsWithCommonDaughter, hfReducedCollision, hfCandD3Pr, hfCandV0, hfTrackNoParam, rowHf3PrV0McRecReduced, rowHf3PrTrkMcRecReduced, dummyTable);
    }
    runMcGen<DMesonType::Dplus, PairingType::V0AndTrack>(particlesMc, mcParticlesPerMcCollision, collInfos, colPerMcCollision, mcCollisions, hfEvSelMc, rejectCollisionsWithBadEvSel, registry, pdg, rowHfResoMcGenReduced, bcs);
    // handle normalization by the right number of collisions
    hfCollisionCounter(collisions.tableSize(), zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
  }
  PROCESS_SWITCH(HfDataCreatorCharmResoToDplusReduced, processDplusV0AndTrackMC, "Process Dplus candidates paired with V0s and tracks with MC matching", false);

  // ML
  void processDplusV0MCWithMl(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                              CandsDplusFilteredWithMlAndMc const& candsDplus,
                              aod::V0s const& v0s,
                              TracksIUWithPIDAndMC const& tracksIU,
                              aod::McParticles const& particlesMc,
                              BCsInfo const& bcs,
                              McCollisionsNoCents const& collInfos,
                              aod::McCollisions const& mcCollisions)
  {
    int zvtxColl{0};
    int sel8Coll{0};
    int zvtxAndSel8Coll{0};
    int zvtxAndSel8CollAndSoftTrig{0};
    int allSelColl{0};
    int dummyTable{0};
    for (const auto& collision : collisions) {
      const auto hfRejMap = o2::hf_evsel::getEvSel<true, o2::hf_centrality::CentralityEstimator::None, BCsInfo>(collision, hfEvSel, zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl, ccdb, registry);
      if (rejectCollisionsWithBadEvSel && hfRejMap != 0) {
        continue;
      }
      auto bc = collision.template bc_as<BCsInfo>();
      if (runNumber != bc.runNumber()) {
        LOG(info) << ">>>>>>>>>>>> Current run number: " << runNumber;
        initCCDB(bc, runNumber, ccdb, ccdbPathGrpMag, lut, false);
        bz = o2::base::Propagator::Instance()->getNominalBz();
        LOG(info) << ">>>>>>>>>>>> Magnetic field: " << bz;
      }
      fitter.setBz(bz);
      auto thisCollId = collision.globalIndex();
      auto candsDThisColl = candsDplus.sliceBy(candsDplusPerCollisionWithMl, thisCollId);
      auto v0sThisColl = v0s.sliceBy(candsV0PerCollision, thisCollId);
      runDataCreation<true, true, DMesonType::Dplus, PairingType::V0Only>(collision, candsDThisColl, v0sThisColl, nullptr, tracksIU, tracksIU, particlesMc, hfRejMap, bz, pdg, registry, matCorr, fitter, cfgDmesCuts, cfgSingleTrackCuts, cfgV0Cuts, cfgQaPlots, rejectPairsWithCommonDaughter, hfReducedCollision, hfCandD3Pr, hfCandV0, dummyTable, rowHf3PrV0McRecReduced, dummyTable, hfCandD3PrMl);
    }
    runMcGen<DMesonType::Dplus, PairingType::V0Only>(particlesMc, mcParticlesPerMcCollision, collInfos, colPerMcCollision, mcCollisions, hfEvSelMc, rejectCollisionsWithBadEvSel, registry, pdg, rowHfResoMcGenReduced, bcs);
    // handle normalization by the right number of collisions
    hfCollisionCounter(collisions.tableSize(), zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
  }
  PROCESS_SWITCH(HfDataCreatorCharmResoToDplusReduced, processDplusV0MCWithMl, "Process Dplus candidates paired with V0s with MC matching and with ML info", false);

  void processDplusTrackMCWithMl(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                                 CandsDplusFilteredWithMlAndMc const& candsDplus,
                                 TracksWithPIDAndMC const& tracks,
                                 aod::TrackAssoc const& trackIndices,
                                 aod::McParticles const& particlesMc,
                                 BCsInfo const& bcs,
                                 McCollisionsNoCents const& collInfos,
                                 aod::McCollisions const& mcCollisions)
  {
    int zvtxColl{0};
    int sel8Coll{0};
    int zvtxAndSel8Coll{0};
    int zvtxAndSel8CollAndSoftTrig{0};
    int allSelColl{0};
    int dummyTable{0};
    for (const auto& collision : collisions) {
      const auto hfRejMap = o2::hf_evsel::getEvSel<true, o2::hf_centrality::CentralityEstimator::None, BCsInfo>(collision, hfEvSel, zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl, ccdb, registry);
      if (rejectCollisionsWithBadEvSel && hfRejMap != 0) {
        continue;
      }
      auto bc = collision.template bc_as<BCsInfo>();
      if (runNumber != bc.runNumber()) {
        LOG(info) << ">>>>>>>>>>>> Current run number: " << runNumber;
        initCCDB(bc, runNumber, ccdb, ccdbPathGrpMag, lut, false);
        bz = o2::base::Propagator::Instance()->getNominalBz();
        LOG(info) << ">>>>>>>>>>>> Magnetic field: " << bz;
      }
      fitter.setBz(bz);
      auto thisCollId = collision.globalIndex();
      auto candsDThisColl = candsDplus.sliceBy(candsDplusPerCollisionWithMl, thisCollId);
      auto trackIdsThisColl = trackIndices.sliceBy(trackIndicesPerCollision, thisCollId);
      runDataCreation<true, true, DMesonType::Dplus, PairingType::TrackOnly>(collision, candsDThisColl, nullptr, trackIdsThisColl, tracks, tracks, particlesMc, hfRejMap, bz, pdg, registry, matCorr, fitter, cfgDmesCuts, cfgSingleTrackCuts, cfgV0Cuts, cfgQaPlots, rejectPairsWithCommonDaughter, hfReducedCollision, hfCandD3Pr, dummyTable, hfTrackNoParam, dummyTable, rowHf3PrTrkMcRecReduced, hfCandD3PrMl);
    }
    runMcGen<DMesonType::Dplus, PairingType::TrackOnly>(particlesMc, mcParticlesPerMcCollision, collInfos, colPerMcCollision, mcCollisions, hfEvSelMc, rejectCollisionsWithBadEvSel, registry, pdg, rowHfResoMcGenReduced, bcs);
    // handle normalization by the right number of collisions
    hfCollisionCounter(collisions.tableSize(), zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
  }
  PROCESS_SWITCH(HfDataCreatorCharmResoToDplusReduced, processDplusTrackMCWithMl, "Process Dplus candidates paired with tracks with MC matching and with ML info", false);

  void processDplusV0AndTrackMCWithMl(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                                      CandsDplusFilteredWithMlAndMc const& candsDplus,
                                      aod::V0s const& v0s,
                                      aod::TrackAssoc const& trackIndices,
                                      TracksWithPIDAndMC const& tracks,
                                      TracksIUWithPIDAndMC const& tracksIU,
                                      aod::McParticles const& particlesMc,
                                      BCsInfo const& bcs,
                                      McCollisionsNoCents const& collInfos,
                                      aod::McCollisions const& mcCollisions)
  {
    int zvtxColl{0};
    int sel8Coll{0};
    int zvtxAndSel8Coll{0};
    int zvtxAndSel8CollAndSoftTrig{0};
    int allSelColl{0};
    for (const auto& collision : collisions) {
      const auto hfRejMap = o2::hf_evsel::getEvSel<true, o2::hf_centrality::CentralityEstimator::None, BCsInfo>(collision, hfEvSel, zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl, ccdb, registry);
      if (rejectCollisionsWithBadEvSel && hfRejMap != 0) {
        continue;
      }
      auto bc = collision.template bc_as<BCsInfo>();
      if (runNumber != bc.runNumber()) {
        LOG(info) << ">>>>>>>>>>>> Current run number: " << runNumber;
        initCCDB(bc, runNumber, ccdb, ccdbPathGrpMag, lut, false);
        bz = o2::base::Propagator::Instance()->getNominalBz();
        LOG(info) << ">>>>>>>>>>>> Magnetic field: " << bz;
      }
      fitter.setBz(bz);
      auto thisCollId = collision.globalIndex();
      auto candsDThisColl = candsDplus.sliceBy(candsDplusPerCollisionWithMl, thisCollId);
      auto v0sThisColl = v0s.sliceBy(candsV0PerCollision, thisCollId);
      auto trackIdsThisColl = trackIndices.sliceBy(trackIndicesPerCollision, thisCollId);
      runDataCreation<true, true, DMesonType::Dplus, PairingType::V0AndTrack>(collision, candsDThisColl, v0sThisColl, trackIdsThisColl, tracks, tracksIU, particlesMc, hfRejMap, bz, pdg, registry, matCorr, fitter, cfgDmesCuts, cfgSingleTrackCuts, cfgV0Cuts, cfgQaPlots, rejectPairsWithCommonDaughter, hfReducedCollision, hfCandD3Pr, hfCandV0, hfTrackNoParam, rowHf3PrV0McRecReduced, rowHf3PrTrkMcRecReduced, hfCandD3PrMl);
    }
    runMcGen<DMesonType::Dplus, PairingType::V0AndTrack>(particlesMc, mcParticlesPerMcCollision, collInfos, colPerMcCollision, mcCollisions, hfEvSelMc, rejectCollisionsWithBadEvSel, registry, pdg, rowHfResoMcGenReduced, bcs);
    // handle normalization by the right number of collisions
    hfCollisionCounter(collisions.tableSize(), zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
  }
  PROCESS_SWITCH(HfDataCreatorCharmResoToDplusReduced, processDplusV0AndTrackMCWithMl, "Process Dplus candidates paired with V0s and tracks with MC matching and with ML info", false);
}; // struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfDataCreatorCharmResoToDplusReduced>(cfgc)};
}
