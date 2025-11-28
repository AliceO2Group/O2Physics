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

/// \file dataCreatorCharmResoToDstarReduced.cxx
/// \brief Creation of D* V0 and D* track pairs
///
/// \author Luca Aglietta <luca.aglietta@cern.ch>, UniTO Turin
/// \author Fabrizio Grosa <fabrizio.grosa@cern.ch>, CERN

#include "PWGHF/Core/CentralityEstimation.h"
#include "PWGHF/Core/DecayChannels.h"
#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/D2H/Core/DataCreationCharmReso.h"
#include "PWGHF/D2H/DataModel/ReducedDataModel.h"
#include "PWGHF/D2H/Utils/utilsRedDataFormat.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/Utils/utilsBfieldCCDB.h"
#include "PWGHF/Utils/utilsEvSelHf.h"
#include "PWGHF/Utils/utilsMcMatching.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CCDB/CcdbApi.h>
#include <CommonConstants/PhysicsConstants.h>
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
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/RunningWorkflowInfo.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>
#include <TPDGCode.h>

#include <Rtypes.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <map>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <vector>

using namespace o2;
using namespace o2::analysis;
using namespace o2::aod;
using namespace o2::analysis::hf_charm_reso;
using namespace o2::constants::physics;
using namespace o2::framework;
using namespace o2::framework::expressions;


/// Creation of D-V0 pairs
struct HfDataCreatorCharmResoToDstarReduced {

  // Produces AOD tables to store collision information
  Produces<aod::HfRedCollisions> hfReducedCollision; // Defined in PWGHF/D2H/DataModel/ReducedDataModel.h
  Produces<aod::HfOrigColCounts> hfCollisionCounter; // Defined in PWGHF/D2H/DataModel/ReducedDataModel.h
  // tracks, V0 and D candidates reduced tables
  Produces<aod::HfRedVzeros> hfCandV0;            // Defined in PWGHF/D2H/DataModel/ReducedDataModel.h
  Produces<aod::HfRedTrkNoParams> hfTrackNoParam; // Defined in PWGHF/D2H/DataModel/ReducedDataModel.h
  Produces<aod::HfRedDstarNoTrks> hfCandDstar;    // Defined in PWGHF/D2H/DataModel/ReducedDataModel.h
  // ML optional Tables
  Produces<aod::HfRed3ProngsMl> hfCandD3PrMl; // Defined in PWGHF/D2H/DataModel/ReducedDataModel.h
  // MC Tables
  Produces<aod::HfMcGenRedResos> rowHfResoMcGenReduced;
  Produces<aod::HfDstarV0McRec> rowHfDstarV0McRecReduced;
  Produces<aod::HfDstarTrkMcRec> rowHfDstarTrkMcRecReduced;

  // selection D
  struct : ConfigurableGroup {
    std::string prefix = "dmesons";
    Configurable<bool> selectionFlagDstarToD0Pi{"selectionFlagDstarToD0Pi", true, "Selection Flag for D* decay to D0 & Pi"};
  } cfgDmesCuts;

  // selection V0
  struct : ConfigurableGroup {
    std::string prefix = "v0s";
    Configurable<float> deltaMassK0s{"deltaMassK0s", 0.02, "delta mass cut for K0S"};
    Configurable<float> deltaMassLambda{"deltaMassLambda", 0.01, "delta mass cut for Lambda"};
    Configurable<float> etaMax{"etaMax", 0.8f, "maximum eta"};
    Configurable<float> etaMaxDau{"etaMaxDau", 5.f, "maximum eta V0 daughters"};
    Configurable<float> trackNclusItsCut{"trackNclusItsCut", 0, "Minimum number of ITS clusters for V0 daughter"};
    Configurable<int> trackNCrossedRowsTpc{"trackNCrossedRowsTpc", 50, "Minimum TPC crossed rows"};
    Configurable<float> trackNsharedClusTpc{"trackNsharedClusTpc", 1000, "Maximum number of shared TPC clusters for V0 daughter"};
    Configurable<float> trackFracMaxindableTpcCls{"trackFracMaxindableTpcCls", 0.8f, "Maximum fraction of findable TPC clusters for V0 daughter"};
    Configurable<float> dcaDau{"dcaDau", 1.f, "DCA V0 daughters"};
    Configurable<float> dcaMaxDauToPv{"dcaMaxDauToPv", 0.1f, "Maximum daughter's DCA to PV"};
    Configurable<float> dcaPv{"dcaPv", 1.f, "DCA V0 to PV"};
    Configurable<double> cosPa{"cosPa", 0.99f, "V0 CosPA"};
    Configurable<float> radiusMin{"radiusMin", 0.9f, "Minimum v0 radius accepted"};
    Configurable<float> nSigmaTpc{"nSigmaTpc", 4.f, "Nsigmatpc"};
    Configurable<float> nSigmaTofPr{"nSigmaTofPr", 4.f, "N sigma TOF for protons only"};
    Configurable<bool> propagateV0toPV{"propagateV0toPV", false, "Enable or disable V0 propagation to V0"};
  } cfgV0Cuts;

  // selection single tracks
  struct : ConfigurableGroup {
    std::string prefix = "singleTracks";
    Configurable<int> setTrackSelections{"setTrackSelections", 2, "flag to apply track selections: 0=none; 1=global track w/o DCA selection; 2=global track; 3=only ITS quality"};
    Configurable<float> maxEta{"maxEta", 0.8, "maximum pseudorapidity for single tracks to be paired with D mesons"};
    Configurable<float> minPt{"minPt", 0.1, "minimum pT for single tracks to be paired with D mesons"};
    Configurable<float> maxNsigmaTpcPi{"maxNsigmaTpcPi", -1., "maximum pion NSigma in TPC for single tracks to be paired with D mesons; set negative to reject"};
    Configurable<float> maxNsigmaTpcKa{"maxNsigmaTpcKa", -1., "maximum kaon NSigma in TPC for single tracks to be paired with D mesons; set negative to reject"};
    Configurable<float> maxNsigmaTpcPr{"maxNsigmaTpcPr", 3., "maximum proton NSigma in TPC for single tracks to be paired with D mesons; set negative to reject"};
  } cfgSingleTrackCuts;

  // QA histograms
  struct : ConfigurableGroup {
    std::string prefix = "qaPlots";
    Configurable<bool> applyCutsForQaHistograms{"applyCutsForQaHistograms", true, "flag to apply cuts to QA histograms"};
    Configurable<float> cutMassDstarMin{"cutMassDstarMin", 0.143, "minimum mass for Dstar candidates"};
    Configurable<float> cutMassDstarMax{"cutMassDstarMax", 0.155, "maximum mass for Dstar candidates"};
    Configurable<float> cutMassK0sMin{"cutMassK0sMin", 0.485, "minimum mass for K0s candidates"};
    Configurable<float> cutMassK0sMax{"cutMassK0sMax", 0.509, "maximum mass for K0s candidates"};
    Configurable<float> cutMassLambdaMin{"cutMassLambdaMin", 1.11, "minimum mass for Lambda candidates"};
    Configurable<float> cutMassLambdaMax{"cutMassLambdaMax", 1.12, "maximum mass for Lambda candidates"};
  } cfgQaPlots;
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

  // Dstar
  using CandsDstarFiltered = soa::Filtered<soa::Join<aod::HfD0FromDstar, aod::HfCandDstars, aod::HfSelDstarToD0Pi>>;
  using CandsDstarFilteredWithMl = soa::Filtered<soa::Join<aod::HfD0FromDstar, aod::HfCandDstars, aod::HfSelDstarToD0Pi, aod::HfMlDstarToD0Pi>>;
  using CandsDstarFilteredWithMc = soa::Filtered<soa::Join<aod::HfD0FromDstar, aod::HfCandDstars, aod::HfSelDstarToD0Pi, aod::HfCandDstarMcRec>>;
  using CandsDstarFilteredWithMlAndMc = soa::Filtered<soa::Join<aod::HfD0FromDstar, aod::HfCandDstars, aod::HfSelDstarToD0Pi, aod::HfMlDstarToD0Pi, aod::HfCandDstarMcRec>>;
  // Tracks
  using TracksWithPID = soa::Join<aod::Tracks, aod::TracksCov, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullKa, aod::pidTOFFullKa, aod::pidTPCFullPr, aod::pidTOFFullPr>;
  using TracksWithPIDAndMC = soa::Join<TracksWithPID, aod::McTrackLabels>;
  using TracksIUWithPID = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::pidTPCPi, aod::pidTOFFullPi, aod::pidTPCPr, aod::pidTOFFullPr>;
  using TracksIUWithPIDAndMC = soa::Join<TracksIUWithPID, aod::McTrackLabels>;
  // Collisions MC
  using BCsInfo = soa::Join<aod::BCs, aod::Timestamps, aod::BcSels>;
  using McCollisionsNoCents = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>;

  Filter filterSelectedCandDstar = (aod::hf_sel_candidate_dstar::isSelDstarToD0Pi == cfgDmesCuts.selectionFlagDstarToD0Pi);

  Preslice<CandsDstarFiltered> candsDstarPerCollision = aod::hf_cand::collisionId;
  Preslice<CandsDstarFilteredWithMl> candsDstarPerCollisionWithMl = aod::hf_cand::collisionId;
  Preslice<aod::V0s> candsV0PerCollision = aod::v0::collisionId;
  Preslice<aod::TrackAssoc> trackIndicesPerCollision = aod::track_association::collisionId;
  Preslice<aod::McParticles> mcParticlesPerMcCollision = aod::mcparticle::mcCollisionId;
  PresliceUnsorted<McCollisionsNoCents> colPerMcCollision = aod::mccollisionlabel::mcCollisionId;

  HistogramRegistry registry{"registry"};

  void init(InitContext& initContext)
  {
    // histograms
    constexpr int kNBinsEvents = kNEvent;
    std::string labels[kNBinsEvents];
    labels[Event::Processed] = "processed";
    labels[Event::NoDV0Selected] = "without DV0 pairs";
    labels[Event::DV0Selected] = "with DV0 pairs";
    static const AxisSpec axisEvents = {kNBinsEvents, 0.5, kNBinsEvents + 0.5, ""};
    registry.add("hEvents", "Events;;entries", HistType::kTH1D, {axisEvents});
    for (int iBin = 0; iBin < kNBinsEvents; iBin++) {
      registry.get<TH1>(HIST("hEvents"))->GetXaxis()->SetBinLabel(iBin + 1, labels[iBin].data());
    }

    const AxisSpec axisPt{50, 0.f, 50.f, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec axisP{100, 0.f, 10.f, "#it{p} (GeV/#it{c})"};
    const AxisSpec axisDeDx{500, 0.f, 1000.f, ""};
    const AxisSpec axisMassDstar{200, 0.139f, 0.179f, "delta inv. mass (GeV/#it{c}^{2})"}; // o2-linter: disable=pdg/explicit-mass (false positive)
    const AxisSpec axisMassLambda{100, 1.05f, 1.35f, "inv. mass (GeV/#it{c}^{2})"};
    const AxisSpec axisMassKzero{100, 0.35f, 0.65f, "inv. mass (GeV/#it{c}^{2})"};
    const AxisSpec axisDeltaMassToK{500, 0.49, 1.49, "inv. mass (GeV/#it{c}^{2})"};
    const AxisSpec axisDeltaMassToPi{500, 0.13, 1.13, "inv. mass (GeV/#it{c}^{2})"};
    const AxisSpec axisDeltaMassToPr{500, 0.93, 1.93, "inv. mass (GeV/#it{c}^{2})"};
    const AxisSpec axisDeltaMassToLambda{500, 1.05, 2.05, "inv. mass (GeV/#it{c}^{2})"};
    const AxisSpec axisMassDsj{400, 0.49f, 0.89f, ""}; // Ds1 and Ds2Star legacy
    registry.add("hMassVsPtK0s", "K0^{s} candidates;#it{p}_{T} (GeV/#it{c});inv. mass (#pi^{#plus}#pi^{#minus}) (GeV/#it{c}^{2});entries", {HistType::kTH2D, {axisPt, axisMassKzero}});
    registry.add("hMassVsPtLambda", "Lambda candidates;#it{p}_{T} (GeV/#it{c});inv. mass (p #pi^{#minus}) (GeV/#it{c}^{2});entries", {HistType::kTH2D, {axisPt, axisMassLambda}});
    registry.add("hdEdxVsP", "Tracks;#it{p} (GeV/#it{c});d#it{E}/d#it{x};entries", {HistType::kTH2D, {axisP, axisDeDx}});
    registry.add("hMassVsPtDstarAll", "Dstar candidates (all, regardless the pairing with V0s);#it{p}_{T} (GeV/#it{c});inv. mass (GeV/#it{c}^{2});entries", {HistType::kTH2D, {axisPt, axisMassDstar}});
    registry.add("hMassVsPtDstarPaired", "Dstar candidates (paired with V0s);#it{p}_{T} (GeV/#it{c});inv. mass (GeV/#it{c}^{2});entries", {HistType::kTH2D, {axisPt, axisMassDstar}});
    registry.add("hMassDstarPi", "DstarPi candidates; m_{D^{*+}#pi^{-}} (GeV/#it{c}^{2});entries", {HistType::kTH2D, {axisPt, axisDeltaMassToPi}});
    registry.add("hMassDstarK", "DstarK candidates; m_{D^{*+}#pi^{-}} (GeV/#it{c}^{2});entries", {HistType::kTH2D, {axisPt, axisDeltaMassToK}});
    registry.add("hMassDstarProton", "DstarProton candidates; m_{D^{*}p} (GeV/#it{c}^{2});entries", {HistType::kTH2D, {axisPt, axisDeltaMassToPr}});
    registry.add("hMassDstarK0s", "DstarK0s candidates; m_{D^{*}K^{0}_{S}} (GeV/#it{c}^{2});entries", {HistType::kTH2D, {axisPt, axisDeltaMassToK}});
    registry.add("hMassDstarLambda", "DstarLambda candidates; m_{D^{*}#Lambda} (GeV/#it{c}^{2});entries", {HistType::kTH2D, {axisPt, axisDeltaMassToLambda}});
    if (doprocessDstarV0MC || doprocessDstarTrackMC || doprocessDstarV0AndTrackMC || doprocessDstarV0MCWithMl || doprocessDstarTrackMCWithMl || doprocessDstarV0AndTrackMCWithMl) {
      // MC Rec
      int const nChannels = hf_decay::hf_cand_reso::DecayChannelMain::NChannelsMain;
      registry.add("hMCRecCounter", "Number of Reconstructed MC Matched candidates per channel", {HistType::kTH1D, {{2 * nChannels + 1, -(nChannels + 0.5), nChannels + 0.5}}});
      registry.add("hMCRecDebug", "Debug of MC Reco", {HistType::kTH1D, {{551, -0.5, 550.5}}});
      registry.add("hMCRecOrigin", "Origin of Matched particles", {HistType::kTH1D, {{3, -0.5, 2.5}}});
      registry.add("hMCRecMassGen", "Generated inv. mass of resoncances", {HistType::kTH1D, {{2000, 1.8, 3.8}}});
      registry.add("hMCRecCharmDau", "Charm daughter flag", {HistType::kTH1D, {{57, -28.5, 28.5}}});
      // MC Gen
      registry.add("hMCGenCounter", "Number of Generated particles; Decay Channel Flag; pT [GeV/c]", {HistType::kTH2D, {{17, -8.5, 8.5}, {100, 0, 50}}});
      registry.add("hMCGenOrigin", "Origin of Generated particles", {HistType::kTH1D, {{3, -0.5, 2.5}}});
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
    hfEvSel.init(registry);

    const auto& workflows = initContext.services().get<RunningWorkflowInfo const>();
    for (const DeviceSpec& device : workflows.devices) {
      if (device.name == "hf-data-creator-charm-reso-to-dstar-reduced") {
        // init HF event selection helper
        hfEvSelMc.init(device, registry);
        break;
      }
    }
  }

  // Process functions
  // No ML
  // Data
  void processDstarV0(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                      CandsDstarFiltered const& candsDstar,
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
      o2::hf_evsel::checkEvSel<true, o2::hf_centrality::CentralityEstimator::None, aod::BCsWithTimestamps>(collision, hfEvSel, zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl, ccdb, registry);
      float centrality = -1.f;
      const auto hfRejMap = hfEvSel.getHfCollisionRejectionMask<true, o2::hf_centrality::CentralityEstimator::None, aod::BCsWithTimestamps>(collision, centrality, ccdb, registry);
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
      auto candsDThisColl = candsDstar.sliceBy(candsDstarPerCollision, thisCollId);
      auto v0sThisColl = v0s.sliceBy(candsV0PerCollision, thisCollId);

      runDataCreation<false, false, DType::Dstar, PairingType::V0Only>(collision, candsDThisColl, v0sThisColl, nullptr, tracksIU, tracksIU, nullptr, hfRejMap, bz, pdg, registry, matCorr, fitter, cfgDmesCuts, cfgSingleTrackCuts, cfgV0Cuts, cfgQaPlots, rejectPairsWithCommonDaughter, hfReducedCollision, hfCandDstar, hfCandV0, dummyTable, dummyTable, dummyTable, dummyTable);
    }
    // handle normalization by the right number of collisions
    hfCollisionCounter(collisions.tableSize(), zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
  }
  PROCESS_SWITCH(HfDataCreatorCharmResoToDstarReduced, processDstarV0, "Process Dstar candidates paired with V0s", true);

  void processDstarTrack(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                         CandsDstarFiltered const& candsDstar,
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
      o2::hf_evsel::checkEvSel<true, o2::hf_centrality::CentralityEstimator::None, aod::BCsWithTimestamps>(collision, hfEvSel, zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl, ccdb, registry);
      float centrality = -1.f;
      const auto hfRejMap = hfEvSel.getHfCollisionRejectionMask<true, o2::hf_centrality::CentralityEstimator::None, aod::BCsWithTimestamps>(collision, centrality, ccdb, registry);
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
      auto candsDThisColl = candsDstar.sliceBy(candsDstarPerCollision, thisCollId);
      auto trackIdsThisColl = trackIndices.sliceBy(trackIndicesPerCollision, thisCollId);

      runDataCreation<false, false, DType::Dstar, PairingType::TrackOnly>(collision, candsDThisColl, nullptr, trackIdsThisColl, tracks, tracks, nullptr, hfRejMap, bz, pdg, registry, matCorr, fitter, cfgDmesCuts, cfgSingleTrackCuts, cfgV0Cuts, cfgQaPlots, rejectPairsWithCommonDaughter, hfReducedCollision, hfCandDstar, dummyTable, hfTrackNoParam, dummyTable, dummyTable, dummyTable);
    }
    // handle normalization by the right number of collisions
    hfCollisionCounter(collisions.tableSize(), zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
  }
  PROCESS_SWITCH(HfDataCreatorCharmResoToDstarReduced, processDstarTrack, "Process Dstar candidates paired with Tracks", false);

  void processDstarV0AndTrack(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                              CandsDstarFiltered const& candsDstar,
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
      o2::hf_evsel::checkEvSel<true, o2::hf_centrality::CentralityEstimator::None, aod::BCsWithTimestamps>(collision, hfEvSel, zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl, ccdb, registry);
      float centrality = -1.f;
      const auto hfRejMap = hfEvSel.getHfCollisionRejectionMask<true, o2::hf_centrality::CentralityEstimator::None, aod::BCsWithTimestamps>(collision, centrality, ccdb, registry);
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
      auto candsDThisColl = candsDstar.sliceBy(candsDstarPerCollision, thisCollId);
      auto v0sThisColl = v0s.sliceBy(candsV0PerCollision, thisCollId);
      auto trackIdsThisColl = trackIndices.sliceBy(trackIndicesPerCollision, thisCollId);
      runDataCreation<false, false, DType::Dstar, PairingType::V0AndTrack>(collision, candsDThisColl, v0sThisColl, trackIdsThisColl, tracks, tracksIU, nullptr, hfRejMap, bz, pdg, registry, matCorr, fitter, cfgDmesCuts, cfgSingleTrackCuts, cfgV0Cuts, cfgQaPlots, rejectPairsWithCommonDaughter, hfReducedCollision, hfCandDstar, hfCandV0, hfTrackNoParam, dummyTable, dummyTable, dummyTable);
    }
    // handle normalization by the right number of collisions
    hfCollisionCounter(collisions.tableSize(), zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
  }
  PROCESS_SWITCH(HfDataCreatorCharmResoToDstarReduced, processDstarV0AndTrack, "Process Dstar candidates paired with V0s and Tracks", false);

  // ML
  // Data
  void processDstarV0WithMl(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                            CandsDstarFilteredWithMl const& candsDstar,
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
      o2::hf_evsel::checkEvSel<true, o2::hf_centrality::CentralityEstimator::None, aod::BCsWithTimestamps>(collision, hfEvSel, zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl, ccdb, registry);
      float centrality = -1.f;
      const auto hfRejMap = hfEvSel.getHfCollisionRejectionMask<true, o2::hf_centrality::CentralityEstimator::None, aod::BCsWithTimestamps>(collision, centrality, ccdb, registry);
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
      auto candsDThisColl = candsDstar.sliceBy(candsDstarPerCollisionWithMl, thisCollId);
      auto v0sThisColl = v0s.sliceBy(candsV0PerCollision, thisCollId);
      runDataCreation<true, false, DType::Dstar, PairingType::V0Only>(collision, candsDThisColl, v0sThisColl, nullptr, tracksIU, tracksIU, nullptr, hfRejMap, bz, pdg, registry, matCorr, fitter, cfgDmesCuts, cfgSingleTrackCuts, cfgV0Cuts, cfgQaPlots, rejectPairsWithCommonDaughter, hfReducedCollision, hfCandDstar, hfCandV0, dummyTable, dummyTable, dummyTable, hfCandD3PrMl);
    }
    // handle normalization by the right number of collisions
    hfCollisionCounter(collisions.tableSize(), zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
  }
  PROCESS_SWITCH(HfDataCreatorCharmResoToDstarReduced, processDstarV0WithMl, "Process Dstar candidates paired with V0s with ML info", false);

  void processDstarTrackWithMl(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                               CandsDstarFilteredWithMl const& candsDstar,
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
      o2::hf_evsel::checkEvSel<true, o2::hf_centrality::CentralityEstimator::None, aod::BCsWithTimestamps>(collision, hfEvSel, zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl, ccdb, registry);
      float centrality = -1.f;
      const auto hfRejMap = hfEvSel.getHfCollisionRejectionMask<true, o2::hf_centrality::CentralityEstimator::None, aod::BCsWithTimestamps>(collision, centrality, ccdb, registry);
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
      auto candsDThisColl = candsDstar.sliceBy(candsDstarPerCollisionWithMl, thisCollId);
      auto trackIdsThisColl = trackIndices.sliceBy(trackIndicesPerCollision, thisCollId);
      runDataCreation<true, false, DType::Dstar, PairingType::TrackOnly>(collision, candsDThisColl, nullptr, trackIdsThisColl, tracks, tracks, nullptr, hfRejMap, bz, pdg, registry, matCorr, fitter, cfgDmesCuts, cfgSingleTrackCuts, cfgV0Cuts, cfgQaPlots, rejectPairsWithCommonDaughter, hfReducedCollision, hfCandDstar, dummyTable, hfTrackNoParam, dummyTable, dummyTable, hfCandD3PrMl);
    }
    // handle normalization by the right number of collisions
    hfCollisionCounter(collisions.tableSize(), zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
  }
  PROCESS_SWITCH(HfDataCreatorCharmResoToDstarReduced, processDstarTrackWithMl, "Process Dstar candidates paired with Tracks with ML info", false);

  void processDstarV0AndTrackWithMl(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                                    CandsDstarFilteredWithMl const& candsDstar,
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
      o2::hf_evsel::checkEvSel<true, o2::hf_centrality::CentralityEstimator::None, aod::BCsWithTimestamps>(collision, hfEvSel, zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl, ccdb, registry);
      float centrality = -1.f;
      const auto hfRejMap = hfEvSel.getHfCollisionRejectionMask<true, o2::hf_centrality::CentralityEstimator::None, aod::BCsWithTimestamps>(collision, centrality, ccdb, registry);
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
      auto candsDThisColl = candsDstar.sliceBy(candsDstarPerCollisionWithMl, thisCollId);
      auto v0sThisColl = v0s.sliceBy(candsV0PerCollision, thisCollId);
      auto trackIdsThisColl = trackIndices.sliceBy(trackIndicesPerCollision, thisCollId);
      runDataCreation<true, false, DType::Dstar, PairingType::V0AndTrack>(collision, candsDThisColl, v0sThisColl, trackIdsThisColl, tracks, tracksIU, nullptr, hfRejMap, bz, pdg, registry, matCorr, fitter, cfgDmesCuts, cfgSingleTrackCuts, cfgV0Cuts, cfgQaPlots, rejectPairsWithCommonDaughter, hfReducedCollision, hfCandDstar, hfCandV0, hfTrackNoParam, dummyTable, dummyTable, hfCandD3PrMl);
    }
    // handle normalization by the right number of collisions
    hfCollisionCounter(collisions.tableSize(), zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
  }
  PROCESS_SWITCH(HfDataCreatorCharmResoToDstarReduced, processDstarV0AndTrackWithMl, "Process Dstar candidates paired with V0s and Tracks with ML info", false);

  // MC
  // No ML
  void processDstarV0MC(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                        CandsDstarFilteredWithMc const& candsDstar,
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
      o2::hf_evsel::checkEvSel<true, o2::hf_centrality::CentralityEstimator::None, BCsInfo>(collision, hfEvSel, zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl, ccdb, registry);
      float centrality = -1.f;
      const auto hfRejMap = hfEvSel.getHfCollisionRejectionMask<true, o2::hf_centrality::CentralityEstimator::None, BCsInfo>(collision, centrality, ccdb, registry);
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
      auto candsDThisColl = candsDstar.sliceBy(candsDstarPerCollision, thisCollId);
      auto v0sThisColl = v0s.sliceBy(candsV0PerCollision, thisCollId);
      runDataCreation<false, true, DType::Dstar, PairingType::V0Only>(collision, candsDThisColl, v0sThisColl, nullptr, tracksIU, tracksIU, particlesMc, hfRejMap, bz, pdg, registry, matCorr, fitter, cfgDmesCuts, cfgSingleTrackCuts, cfgV0Cuts, cfgQaPlots, rejectPairsWithCommonDaughter, hfReducedCollision, hfCandDstar, hfCandV0, dummyTable, rowHfDstarV0McRecReduced, dummyTable, dummyTable);
    }
    runMcGen<DType::Dstar, PairingType::V0Only>(particlesMc, mcParticlesPerMcCollision, collInfos, colPerMcCollision, mcCollisions, hfEvSelMc, rejectCollisionsWithBadEvSel, registry, pdg, rowHfResoMcGenReduced, bcs);
    // handle normalization by the right number of collisions
    hfCollisionCounter(collisions.tableSize(), zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
  }
  PROCESS_SWITCH(HfDataCreatorCharmResoToDstarReduced, processDstarV0MC, "Process Dstar candidates paired with V0s with MC matching", false);

  void processDstarTrackMC(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                           CandsDstarFilteredWithMc const& candsDstar,
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
      o2::hf_evsel::checkEvSel<true, o2::hf_centrality::CentralityEstimator::None, BCsInfo>(collision, hfEvSel, zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl, ccdb, registry);
      float centrality = -1.f;
      const auto hfRejMap = hfEvSel.getHfCollisionRejectionMask<true, o2::hf_centrality::CentralityEstimator::None, BCsInfo>(collision, centrality, ccdb, registry);
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
      auto candsDThisColl = candsDstar.sliceBy(candsDstarPerCollision, thisCollId);
      auto trackIdsThisColl = trackIndices.sliceBy(trackIndicesPerCollision, thisCollId);
      runDataCreation<false, true, DType::Dstar, PairingType::TrackOnly>(collision, candsDThisColl, nullptr, trackIdsThisColl, tracks, tracks, particlesMc, hfRejMap, bz, pdg, registry, matCorr, fitter, cfgDmesCuts, cfgSingleTrackCuts, cfgV0Cuts, cfgQaPlots, rejectPairsWithCommonDaughter, hfReducedCollision, hfCandDstar, dummyTable, hfTrackNoParam, dummyTable, rowHfDstarTrkMcRecReduced, dummyTable);
    }
    runMcGen<DType::Dstar, PairingType::TrackOnly>(particlesMc, mcParticlesPerMcCollision, collInfos, colPerMcCollision, mcCollisions, hfEvSelMc, rejectCollisionsWithBadEvSel, registry, pdg, rowHfResoMcGenReduced, bcs);
    // handle normalization by the right number of collisions
    hfCollisionCounter(collisions.tableSize(), zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
  }
  PROCESS_SWITCH(HfDataCreatorCharmResoToDstarReduced, processDstarTrackMC, "Process Dstar candidates paired with tracks with MC matching", false);

  void processDstarV0AndTrackMC(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                                CandsDstarFilteredWithMc const& candsDstar,
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
      o2::hf_evsel::checkEvSel<true, o2::hf_centrality::CentralityEstimator::None, BCsInfo>(collision, hfEvSel, zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl, ccdb, registry);
      float centrality = -1.f;
      const auto hfRejMap = hfEvSel.getHfCollisionRejectionMask<true, o2::hf_centrality::CentralityEstimator::None, BCsInfo>(collision, centrality, ccdb, registry);
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
      auto candsDThisColl = candsDstar.sliceBy(candsDstarPerCollision, thisCollId);
      auto v0sThisColl = v0s.sliceBy(candsV0PerCollision, thisCollId);
      auto trackIdsThisColl = trackIndices.sliceBy(trackIndicesPerCollision, thisCollId);
      runDataCreation<false, true, DType::Dstar, PairingType::V0AndTrack>(collision, candsDThisColl, v0sThisColl, trackIdsThisColl, tracks, tracksIU, particlesMc, hfRejMap, bz, pdg, registry, matCorr, fitter, cfgDmesCuts, cfgSingleTrackCuts, cfgV0Cuts, cfgQaPlots, rejectPairsWithCommonDaughter, hfReducedCollision, hfCandDstar, hfCandV0, hfTrackNoParam, rowHfDstarV0McRecReduced, rowHfDstarTrkMcRecReduced, dummyTable);
    }
    runMcGen<DType::Dstar, PairingType::V0AndTrack>(particlesMc, mcParticlesPerMcCollision, collInfos, colPerMcCollision, mcCollisions, hfEvSelMc, rejectCollisionsWithBadEvSel, registry, pdg, rowHfResoMcGenReduced, bcs);
    // handle normalization by the right number of collisions
    hfCollisionCounter(collisions.tableSize(), zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
  }
  PROCESS_SWITCH(HfDataCreatorCharmResoToDstarReduced, processDstarV0AndTrackMC, "Process Dstar candidates paired with V0s and tracks with MC matching", false);

  // ML
  void processDstarV0MCWithMl(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                              CandsDstarFilteredWithMlAndMc const& candsDstar,
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
      o2::hf_evsel::checkEvSel<true, o2::hf_centrality::CentralityEstimator::None, BCsInfo>(collision, hfEvSel, zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl, ccdb, registry);
      float centrality = -1.f;
      const auto hfRejMap = hfEvSel.getHfCollisionRejectionMask<true, o2::hf_centrality::CentralityEstimator::None, BCsInfo>(collision, centrality, ccdb, registry);
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
      auto candsDThisColl = candsDstar.sliceBy(candsDstarPerCollisionWithMl, thisCollId);
      auto v0sThisColl = v0s.sliceBy(candsV0PerCollision, thisCollId);
      runDataCreation<true, true, DType::Dstar, PairingType::V0Only>(collision, candsDThisColl, v0sThisColl, nullptr, tracksIU, tracksIU, particlesMc, hfRejMap, bz, pdg, registry, matCorr, fitter, cfgDmesCuts, cfgSingleTrackCuts, cfgV0Cuts, cfgQaPlots, rejectPairsWithCommonDaughter, hfReducedCollision, hfCandDstar, hfCandV0, dummyTable, rowHfDstarV0McRecReduced, dummyTable, hfCandD3PrMl);
    }
    runMcGen<DType::Dstar, PairingType::V0Only>(particlesMc, mcParticlesPerMcCollision, collInfos, colPerMcCollision, mcCollisions, hfEvSelMc, rejectCollisionsWithBadEvSel, registry, pdg, rowHfResoMcGenReduced, bcs);
    // handle normalization by the right number of collisions
    hfCollisionCounter(collisions.tableSize(), zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
  }
  PROCESS_SWITCH(HfDataCreatorCharmResoToDstarReduced, processDstarV0MCWithMl, "Process Dstar candidates paired with V0s with MC matching and with ML info", false);

  void processDstarTrackMCWithMl(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                                 CandsDstarFilteredWithMlAndMc const& candsDstar,
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
      o2::hf_evsel::checkEvSel<true, o2::hf_centrality::CentralityEstimator::None, BCsInfo>(collision, hfEvSel, zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl, ccdb, registry);
      float centrality = -1.f;
      const auto hfRejMap = hfEvSel.getHfCollisionRejectionMask<true, o2::hf_centrality::CentralityEstimator::None, BCsInfo>(collision, centrality, ccdb, registry);
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
      auto candsDThisColl = candsDstar.sliceBy(candsDstarPerCollisionWithMl, thisCollId);
      auto trackIdsThisColl = trackIndices.sliceBy(trackIndicesPerCollision, thisCollId);
      runDataCreation<true, true, DType::Dstar, PairingType::TrackOnly>(collision, candsDThisColl, nullptr, trackIdsThisColl, tracks, tracks, particlesMc, hfRejMap, bz, pdg, registry, matCorr, fitter, cfgDmesCuts, cfgSingleTrackCuts, cfgV0Cuts, cfgQaPlots, rejectPairsWithCommonDaughter, hfReducedCollision, hfCandDstar, dummyTable, hfTrackNoParam, dummyTable, rowHfDstarTrkMcRecReduced, hfCandD3PrMl);
    }
    runMcGen<DType::Dstar, PairingType::TrackOnly>(particlesMc, mcParticlesPerMcCollision, collInfos, colPerMcCollision, mcCollisions, hfEvSelMc, rejectCollisionsWithBadEvSel, registry, pdg, rowHfResoMcGenReduced, bcs);
    // handle normalization by the right number of collisions
    hfCollisionCounter(collisions.tableSize(), zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
  }
  PROCESS_SWITCH(HfDataCreatorCharmResoToDstarReduced, processDstarTrackMCWithMl, "Process Dstar candidates paired with tracks with MC matching and with ML info", false);

  void processDstarV0AndTrackMCWithMl(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                                      CandsDstarFilteredWithMlAndMc const& candsDstar,
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
      o2::hf_evsel::checkEvSel<true, o2::hf_centrality::CentralityEstimator::None, BCsInfo>(collision, hfEvSel, zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl, ccdb, registry);
      float centrality = -1.f;
      const auto hfRejMap = hfEvSel.getHfCollisionRejectionMask<true, o2::hf_centrality::CentralityEstimator::None, BCsInfo>(collision, centrality, ccdb, registry);
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
      auto candsDThisColl = candsDstar.sliceBy(candsDstarPerCollisionWithMl, thisCollId);
      auto v0sThisColl = v0s.sliceBy(candsV0PerCollision, thisCollId);
      auto trackIdsThisColl = trackIndices.sliceBy(trackIndicesPerCollision, thisCollId);
      runDataCreation<true, true, DType::Dstar, PairingType::V0AndTrack>(collision, candsDThisColl, v0sThisColl, trackIdsThisColl, tracks, tracksIU, particlesMc, hfRejMap, bz, pdg, registry, matCorr, fitter, cfgDmesCuts, cfgSingleTrackCuts, cfgV0Cuts, cfgQaPlots, rejectPairsWithCommonDaughter, hfReducedCollision, hfCandDstar, hfCandV0, hfTrackNoParam, rowHfDstarV0McRecReduced, rowHfDstarTrkMcRecReduced, hfCandD3PrMl);
    }
    runMcGen<DType::Dstar, PairingType::V0AndTrack>(particlesMc, mcParticlesPerMcCollision, collInfos, colPerMcCollision, mcCollisions, hfEvSelMc, rejectCollisionsWithBadEvSel, registry, pdg, rowHfResoMcGenReduced, bcs);
    // handle normalization by the right number of collisions
    hfCollisionCounter(collisions.tableSize(), zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
  }
  PROCESS_SWITCH(HfDataCreatorCharmResoToDstarReduced, processDstarV0AndTrackMCWithMl, "Process Dstar candidates paired with V0s and tracks with MC matching and with ML info", false);
}; // struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfDataCreatorCharmResoToDstarReduced>(cfgc)};
}
