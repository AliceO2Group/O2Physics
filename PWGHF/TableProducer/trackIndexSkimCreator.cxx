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

/// \file trackIndexSkimCreator.cxx
/// \brief Pre-selection of 2-prong, 3-prong, D*, and cascade secondary vertices of heavy-flavour decay candidates
///
/// \author Gian Michele Innocenti <gian.michele.innocenti@cern.ch>, CERN
/// \author Vít Kučera <vit.kucera@cern.ch>, CERN
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>, CERN
/// \author Mattia Faggin <mfaggin@cern.ch>, University and INFN Padova
/// \author Jinjoo Seo <jseo@cern.ch>, Inha University
/// \author Fabrizio Grosa <fgrosa@cern.ch>, CERN
/// \author Federica Zanone <federica.zanone@cern.ch>, Heidelberg University
/// \author Ruiqi Yin <ruiqi.yin@cern.ch>, Fudan University

#include "PWGHF/Core/CentralityEstimation.h"
#include "PWGHF/Core/SelectorCuts.h"
#include "PWGHF/DataModel/AliasTables.h"
#include "PWGHF/DataModel/TrackIndexSkimmingTables.h"
#include "PWGHF/Utils/utilsAnalysis.h"
#include "PWGHF/Utils/utilsBfieldCCDB.h"
#include "PWGHF/Utils/utilsEvSelHf.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/CCDB/TriggerAliases.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelectorPID.h"
#include "Common/Core/ZorroSummary.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Tools/ML/MlResponse.h"

#include <CCDB/BasicCCDBManager.h> // for PV refit
#include <CCDB/CcdbApi.h>
#include <CommonConstants/PhysicsConstants.h>
#include <CommonUtils/ConfigurableParam.h>
#include <DCAFitter/DCAFitterN.h>
#include <DetectorsBase/MatLayerCylSet.h>
#include <DetectorsBase/Propagator.h>     // for PV refit
#include <DetectorsVertexing/PVertexer.h> // for PV refit
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Array2D.h>
#include <Framework/Configurable.h>
#include <Framework/DataTypes.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/runDataProcessing.h>
#include <ReconstructionDataFormats/Track.h>
#include <ReconstructionDataFormats/Vertex.h> // for PV refit

#include <TH1.h>
#include <TString.h>

#include <sys/types.h>

#include <Rtypes.h>

#include <algorithm> // std::find
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iterator> // std::distance
#include <numeric>
#include <string>  // std::string
#include <utility> // std::forward
#include <vector>  // std::vector

using namespace o2;
using namespace o2::analysis;
using namespace o2::hf_evsel;
using namespace o2::aod;
using namespace o2::hf_centrality;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;

// enum for candidate type
enum CandidateType {
  Cand2Prong = 0,
  Cand3Prong,
  CandV0bachelor,
  CandDstar,
  CandCascadeBachelor,
  NCandidateTypes
};

// enum for proton PID strategy (only proton for baryons)
enum ProtonPidStrategy {
  NoPid = 0,
  PidTpcOnly,
  PidTofOnly,
  PidTpcOrTof,
  PidTpcAndTof,
  NPidStatus
};

// enum for proton PID channels
enum ChannelsProtonPid {
  LcToPKPi = 0,
  XicToPKPi,
  LcToPK0S,
  NChannelsProtonPid
};
// kaon PID (opposite-sign track in 3-prong decays)
constexpr int ChannelKaonPid = ChannelsProtonPid::NChannelsProtonPid;
constexpr int ChannelsDeuteronPid = ChannelsProtonPid::NChannelsProtonPid + 1;

/// Event selection
struct HfTrackIndexSkimCreatorTagSelCollisions {
  Produces<aod::HfSelCollision> rowSelectedCollision;

  Configurable<bool> fillHistograms{"fillHistograms", true, "fill histograms"};
  Configurable<std::string> triggerClassName{"triggerClassName", "kINT7", "Run 2 trigger class, only for Run 2 converted data"};
  HfEventSelection hfEvSel;                   // event selection and monitoring
  Service<o2::ccdb::BasicCCDBManager> ccdb{}; // needed for evSelection

  // QA histos
  HistogramRegistry registry{"registry"};
  OutputObj<ZorroSummary> zorroSummary{"zorroSummary"};

  void init(InitContext const&)
  {
    const std::array<int, 7> doProcess = {doprocessTrigAndCentFT0ASel, doprocessTrigAndCentFT0CSel, doprocessTrigAndCentFT0MSel, doprocessTrigAndCentFV0ASel, doprocessTrigSel, doprocessNoTrigSel, doprocessUpcSel};
    if (std::accumulate(doProcess.begin(), doProcess.end(), 0) != 1) {
      LOGP(fatal, "One and only one process function for collision selection can be enabled at a time!");
    }

    // set numerical value of the Run 2 trigger class
    auto* const triggerAlias = std::find(aliasLabels, aliasLabels + kNaliases, triggerClassName.value.data());
    if (triggerAlias != aliasLabels + kNaliases) {
      hfEvSel.triggerClass.value = std::distance(aliasLabels, triggerAlias);
    }

    hfEvSel.init(registry, &zorroSummary); // collision monitoring
    if (fillHistograms) {
      if (doprocessTrigAndCentFT0ASel || doprocessTrigAndCentFT0CSel || doprocessTrigAndCentFT0MSel || doprocessTrigAndCentFV0ASel) {
        const AxisSpec axisCentrality{200, 0., 100., "centrality percentile"};
        registry.add("hCentralitySelected", "Centrality percentile of selected events in the centrality interval; centrality percentile;entries", {HistType::kTH1D, {axisCentrality}});
        registry.add("hCentralityRejected", "Centrality percentile of selected events outside the centrality interval; centrality percentile;entries", {HistType::kTH1D, {axisCentrality}});
      }
    }
  }

  /// Collision selection
  /// \param collision  collision table with
  template <bool ApplyTrigSel, bool ApplyUpcSel, o2::hf_centrality::CentralityEstimator CentEstimator, typename Col, typename BCsType>
  void selectCollision(const Col& collision,
                       const BCsType& bcs)
  {
    float centrality{-1.f};
    o2::hf_evsel::HfCollisionRejectionMask rejectionMask{};

    if constexpr (ApplyUpcSel) {
      rejectionMask = hfEvSel.getHfCollisionRejectionMaskWithUpc<ApplyTrigSel, CentEstimator>(
        collision, centrality, ccdb, registry, bcs);
    } else {
      rejectionMask = hfEvSel.getHfCollisionRejectionMask<ApplyTrigSel, CentEstimator, BCsType>(
        collision, centrality, ccdb, registry);
    }

    if (fillHistograms) {
      hfEvSel.fillHistograms(collision, rejectionMask, centrality);
      // additional centrality histos
      if constexpr (CentEstimator != o2::hf_centrality::None) {
        if (rejectionMask == 0) {
          registry.fill(HIST("hCentralitySelected"), centrality);
        } else if (rejectionMask == BIT(EventRejection::Centrality)) { // rejected by centrality only
          registry.fill(HIST("hCentralityRejected"), centrality);
        }
      }
    }

    // fill table row
    rowSelectedCollision(rejectionMask);
  }

  /// Event selection with trigger and FT0A centrality selection
  void processTrigAndCentFT0ASel(soa::Join<aod::Collisions,
                                           aod::EvSels, aod::CentFT0As>::iterator const& collision,
                                 aod::BcFullInfos const& bcs)
  {
    selectCollision<true, false, CentralityEstimator::FT0A>(collision, bcs);
  }
  PROCESS_SWITCH(HfTrackIndexSkimCreatorTagSelCollisions, processTrigAndCentFT0ASel, "Use trigger and centrality selection with FT0A", false);

  /// Event selection with trigger and FT0C centrality selection
  void processTrigAndCentFT0CSel(soa::Join<aod::Collisions,
                                           aod::EvSels, aod::CentFT0Cs>::iterator const& collision,
                                 aod::BcFullInfos const& bcs)
  {
    selectCollision<true, false, CentralityEstimator::FT0C>(collision, bcs);
  }
  PROCESS_SWITCH(HfTrackIndexSkimCreatorTagSelCollisions, processTrigAndCentFT0CSel, "Use trigger and centrality selection with FT0C", false);

  /// Event selection with trigger and FT0M centrality selection
  void processTrigAndCentFT0MSel(soa::Join<aod::Collisions,
                                           aod::EvSels, aod::CentFT0Ms>::iterator const& collision,
                                 aod::BcFullInfos const& bcs)
  {
    selectCollision<true, false, CentralityEstimator::FT0M>(collision, bcs);
  }
  PROCESS_SWITCH(HfTrackIndexSkimCreatorTagSelCollisions, processTrigAndCentFT0MSel, "Use trigger and centrality selection with FT0M", false);

  /// Event selection with trigger and FV0A centrality selection
  void processTrigAndCentFV0ASel(soa::Join<aod::Collisions,
                                           aod::EvSels, aod::CentFV0As>::iterator const& collision,
                                 aod::BcFullInfos const& bcs)
  {
    selectCollision<true, false, CentralityEstimator::FV0A>(collision, bcs);
  }
  PROCESS_SWITCH(HfTrackIndexSkimCreatorTagSelCollisions, processTrigAndCentFV0ASel, "Use trigger and centrality selection with FV0A", false);

  /// Event selection with trigger selection
  void processTrigSel(soa::Join<aod::Collisions,
                                aod::EvSels>::iterator const& collision,
                      aod::BcFullInfos const& bcs)
  {
    selectCollision<true, false, CentralityEstimator::None>(collision, bcs);
  }
  PROCESS_SWITCH(HfTrackIndexSkimCreatorTagSelCollisions, processTrigSel, "Use trigger selection", false);

  /// Event selection without trigger selection
  void processNoTrigSel(aod::Collision const& collision,
                        aod::BcFullInfos const& bcs)
  {
    selectCollision<false, false, CentralityEstimator::None>(collision, bcs);
  }
  PROCESS_SWITCH(HfTrackIndexSkimCreatorTagSelCollisions, processNoTrigSel, "Do not use trigger selection", true);

  /// Event selection with UPC
  void processUpcSel(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision,
                     aod::BcFullInfos const& bcs,
                     aod::FT0s const& /*ft0s*/,
                     aod::FV0As const& /*fv0as*/,
                     aod::FDDs const& /*fdds*/,
                     aod::Zdcs const& /*zdcs*/)
  {
    selectCollision<true, true, CentralityEstimator::None>(collision, bcs);
  }
  PROCESS_SWITCH(HfTrackIndexSkimCreatorTagSelCollisions, processUpcSel, "Use UPC event selection", false);
};

/// Track selection
struct HfTrackIndexSkimCreatorTagSelTracks {
  Produces<aod::HfSelTrack> rowSelectedTrack;
  Produces<aod::HfPvRefitTrack> tabPvRefitTrack;

  struct : ConfigurableGroup {
    double etaMinDefault{-99999.};
    Configurable<bool> isRun2{"isRun2", false, "enable Run 2 or Run 3 GRP objects for magnetic field"};
    Configurable<bool> doPvRefit{"doPvRefit", false, "do PV refit excluding the considered track"};
    Configurable<bool> fillHistograms{"fillHistograms", true, "fill histograms"};
    Configurable<bool> debugPvRefit{"debugPvRefit", false, "debug lines for primary vertex refit"};
    // Configurable<double> bz{"bz", 5., "bz field"};
    // quality cut
    Configurable<bool> doCutQuality{"doCutQuality", true, "apply quality cuts"};
    Configurable<bool> useIsGlobalTrack{"useIsGlobalTrack", false, "check isGlobalTrack status for tracks, for Run3 studies"};
    Configurable<bool> useIsGlobalTrackWoDCA{"useIsGlobalTrackWoDCA", false, "check isGlobalTrackWoDCA status for tracks, for Run3 studies"};
    Configurable<int> tpcNClsFoundMin{"tpcNClsFoundMin", 70, "min. number of found TPC clusters"};
    // pT bins for single-track cuts
    Configurable<std::vector<double>> binsPtTrack{"binsPtTrack", std::vector<double>{hf_cuts_single_track::vecBinsPtTrack}, "track pT bin limits for 2-prong DCA XY pT-dependent cut"};
    // 2-prong cuts
    Configurable<double> ptMinTrack2Prong{"ptMinTrack2Prong", -1., "min. track pT for 2 prong candidate"};
    Configurable<LabeledArray<double>> cutsTrack2Prong{"cutsTrack2Prong", {hf_cuts_single_track::CutsTrack[0], hf_cuts_single_track::NBinsPtTrack, hf_cuts_single_track::NCutVarsTrack, hf_cuts_single_track::labelsPtTrack, hf_cuts_single_track::labelsCutVarTrack}, "Single-track selections per pT bin for 2-prong candidates"};
    Configurable<double> etaMinTrack2Prong{"etaMinTrack2Prong", std::forward<double>(etaMinDefault), "min. pseudorapidity for 2 prong candidate"};
    Configurable<double> etaMaxTrack2Prong{"etaMaxTrack2Prong", 4., "max. pseudorapidity for 2 prong candidate"};
    // 3-prong cuts
    Configurable<double> ptMinTrack3Prong{"ptMinTrack3Prong", -1., "min. track pT for 3 prong candidate"};
    Configurable<LabeledArray<double>> cutsTrack3Prong{"cutsTrack3Prong", {hf_cuts_single_track::CutsTrack[0], hf_cuts_single_track::NBinsPtTrack, hf_cuts_single_track::NCutVarsTrack, hf_cuts_single_track::labelsPtTrack, hf_cuts_single_track::labelsCutVarTrack}, "Single-track selections per pT bin for 3-prong candidates"};
    Configurable<double> etaMinTrack3Prong{"etaMinTrack3Prong", std::forward<double>(etaMinDefault), "min. pseudorapidity for 3 prong candidate"};
    Configurable<double> etaMaxTrack3Prong{"etaMaxTrack3Prong", 4., "max. pseudorapidity for 3 prong candidate"};
    // bachelor cuts (V0 + bachelor decays)
    Configurable<double> ptMinTrackBach{"ptMinTrackBach", 0.3, "min. track pT for bachelor in cascade candidate"}; // 0.5 for PbPb 2015?
    Configurable<LabeledArray<double>> cutsTrackBach{"cutsTrackBach", {hf_cuts_single_track::CutsTrack[0], hf_cuts_single_track::NBinsPtTrack, hf_cuts_single_track::NCutVarsTrack, hf_cuts_single_track::labelsPtTrack, hf_cuts_single_track::labelsCutVarTrack}, "Single-track selections per pT bin for the bachelor of V0-bachelor candidates"};
    Configurable<double> etaMinTrackBach{"etaMinTrackBach", std::forward<double>(etaMinDefault), "min. pseudorapidity for bachelor in cascade candidate"};
    Configurable<double> etaMaxTrackBach{"etaMaxTrackBach", 0.8, "max. pseudorapidity for bachelor in cascade candidate"};
    // bachelor cuts (cascade + bachelor decays)
    Configurable<double> ptMinTrackBachLfCasc{"ptMinTrackBachLfCasc", 0.1, "min. track pT for bachelor in cascade + bachelor decays"}; // 0.5 for PbPb 2015?
    Configurable<LabeledArray<double>> cutsTrackBachLfCasc{"cutsTrackBachLfCasc", {hf_cuts_single_track::CutsTrack[0], hf_cuts_single_track::NBinsPtTrack, hf_cuts_single_track::NCutVarsTrack, hf_cuts_single_track::labelsPtTrack, hf_cuts_single_track::labelsCutVarTrack}, "Single-track selections per pT bin for the bachelor in cascade + bachelor decays"};
    Configurable<double> etaMinTrackBachLfCasc{"etaMinTrackBachLfCasc", std::forward<double>(etaMinDefault), "min. pseudorapidity for bachelor in cascade + bachelor decays"};
    Configurable<double> etaMaxTrackBachLfCasc{"etaMaxTrackBachLfCasc", 1.1, "max. pseudorapidity for bachelor in cascade + bachelor decays"};
    Configurable<bool> useIsGlobalTrackForBachLfCasc{"useIsGlobalTrackForBachLfCasc", false, "check isGlobalTrack status for bachelor in cascade + bachelor decays"};
    Configurable<bool> useIsGlobalTrackWoDCAForBachLfCasc{"useIsGlobalTrackWoDCAForBachLfCasc", false, "check isGlobalTrackWoDCA status for bachelor in cascade + bachelor decays"};
    Configurable<bool> useIsQualityTrackITSForBachLfCasc{"useIsQualityTrackITSForBachLfCasc", true, "check isQualityTrackITS status for bachelor in cascade + bachelor decays"};
    // soft pion cuts for D*
    Configurable<double> ptMinSoftPionForDstar{"ptMinSoftPionForDstar", 0.05, "min. track pT for soft pion in D* candidate"};
    Configurable<double> etaMinSoftPionForDstar{"etaMinSoftPionForDstar", std::forward<double>(etaMinDefault), "min. pseudorapidity for soft pion in D* candidate"};
    Configurable<double> etaMaxSoftPionForDstar{"etaMaxSoftPionForDstar", 0.8, "max. pseudorapidity for soft pion in D* candidate"};
    Configurable<LabeledArray<double>> cutsTrackDstar{"cutsTrackDstar", {hf_cuts_single_track::CutsTrackPrimary[0], hf_cuts_single_track::NBinsPtTrack, hf_cuts_single_track::NCutVarsTrack, hf_cuts_single_track::labelsPtTrack, hf_cuts_single_track::labelsCutVarTrack}, "Single-track selections per pT bin for the soft pion of D* candidates"};
    Configurable<bool> useIsGlobalTrackForSoftPion{"useIsGlobalTrackForSoftPion", false, "check isGlobalTrack status for soft pion tracks"};
    Configurable<bool> useIsGlobalTrackWoDCAForSoftPion{"useIsGlobalTrackWoDCAForSoftPion", false, "check isGlobalTrackWoDCA status for soft pion tracks"};
    Configurable<bool> useIsQualityTrackITSForSoftPion{"useIsQualityTrackITSForSoftPion", true, "check qualityTracksITS status for soft pion tracks"};
    // proton PID, applied only if corresponding process function enabled
    Configurable<LabeledArray<float>> selectionsPid{"selectionsPid", {hf_presel_pid::CutsPid[0], 5, 6, hf_presel_pid::labelsRowsPid, hf_presel_pid::labelsCutsPid}, "PID selections for proton / kaon / deuteron applied if proper process function enabled"};
    // CCDB
    Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
    Configurable<std::string> ccdbPathLut{"ccdbPathLut", "GLO/Param/MatLUT", "Path for LUT parametrization"};
    Configurable<std::string> ccdbPathGrp{"ccdbPathGrp", "GLO/GRP/GRP", "Path of the grp file (Run 2)"};
    Configurable<std::string> ccdbPathGrpMag{"ccdbPathGrpMag", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object (Run 3)"};
  } config;

  SliceCache cache;

  // Needed for PV refitting
  Service<o2::ccdb::BasicCCDBManager> ccdb{};
  o2::base::MatLayerCylSet* lut{};
  o2::base::Propagator::MatCorrType noMatCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;
  int runNumber{};

  using TracksWithSelAndDca = soa::Join<aod::TracksWCovDcaExtra, aod::TrackSelection>;
  using TracksWithSelAndDcaAndPidTpc = soa::Join<aod::TracksWCovDcaExtra, aod::TrackSelection, aod::pidTPCFullPr, aod::pidTPCFullKa, aod::pidTPCFullDe>;
  using TracksWithSelAndDcaAndPidTof = soa::Join<aod::TracksWCovDcaExtra, aod::TrackSelection, aod::pidTOFFullPr, aod::pidTOFFullKa, aod::pidTOFFullDe>;
  using TracksWithSelAndDcaAndPidTpcTof = soa::Join<aod::TracksWCovDcaExtra, aod::TrackSelection, aod::pidTPCFullPr, aod::pidTOFFullPr, aod::pidTPCFullKa, aod::pidTOFFullKa, aod::pidTPCFullDe, aod::pidTOFFullDe>;

  Preslice<aod::Tracks> perCol = aod::track::collisionId;
  Preslice<TrackAssoc> trackIndicesPerCollision = aod::track_association::collisionId;

  // single-track cuts
  static const int nCuts = 4;
  // array of 2-prong and 3-prong cuts
  std::array<LabeledArray<double>, CandidateType::NCandidateTypes> cutsSingleTrack{};
  // proton PID, if enabled
  std::array<TrackSelectorPr, ChannelsProtonPid::NChannelsProtonPid> selectorProton{};
  TrackSelectorKa selectorKaon;
  TrackSelectorDe selectorDeuteron;

  Partition<TracksWithSelAndDca> pvContributors = ((aod::track::flags & static_cast<uint32_t>(aod::track::PVContributor)) == static_cast<uint32_t>(aod::track::PVContributor));
  Partition<TracksWithSelAndDcaAndPidTpc> pvContributorsWithPidTpc = ((aod::track::flags & static_cast<uint32_t>(aod::track::PVContributor)) == static_cast<uint32_t>(aod::track::PVContributor));
  Partition<TracksWithSelAndDcaAndPidTof> pvContributorsWithPidTof = ((aod::track::flags & static_cast<uint32_t>(aod::track::PVContributor)) == static_cast<uint32_t>(aod::track::PVContributor));
  Partition<TracksWithSelAndDcaAndPidTpcTof> pvContributorsWithPidTpcTof = ((aod::track::flags & static_cast<uint32_t>(aod::track::PVContributor)) == static_cast<uint32_t>(aod::track::PVContributor));

  // QA of PV refit
  ConfigurableAxis axisPvRefitDeltaX{"axisPvRefitDeltaX", {1000, -0.5f, 0.5f}, "DeltaX binning PV refit"};
  ConfigurableAxis axisPvRefitDeltaY{"axisPvRefitDeltaY", {1000, -0.5f, 0.5f}, "DeltaY binning PV refit"};
  ConfigurableAxis axisPvRefitDeltaZ{"axisPvRefitDeltaZ", {1000, -0.5f, 0.5f}, "DeltaZ binning PV refit"};

  HistogramRegistry registry{"registry"};

  void init(InitContext const&)
  {
    const std::array<int, 5> doProcess = {doprocessNoPid, doprocessProtonPidTpc, doprocessProtonPidTof, doprocessProtonPidTpcOrTof, doprocessProtonPidTpcAndTof};
    if (std::accumulate(doProcess.begin(), doProcess.end(), 0) != 1) {
      LOGP(fatal, "One and only one process function for the different PID selection strategies can be enabled at a time!");
    }

    cutsSingleTrack = {config.cutsTrack2Prong, config.cutsTrack3Prong, config.cutsTrackBach, config.cutsTrackDstar, config.cutsTrackBachLfCasc};

    if (config.etaMinTrack2Prong == config.etaMinDefault) {
      config.etaMinTrack2Prong.value = -config.etaMaxTrack2Prong;
    }
    if (config.etaMinTrack3Prong == config.etaMinDefault) {
      config.etaMinTrack3Prong.value = -config.etaMaxTrack3Prong;
    }
    if (config.etaMinTrackBach == config.etaMinDefault) {
      config.etaMinTrackBach.value = -config.etaMaxTrackBach;
    }
    if (config.etaMinSoftPionForDstar == config.etaMinDefault) {
      config.etaMinSoftPionForDstar.value = -config.etaMaxSoftPionForDstar;
    }
    if (config.etaMinTrackBachLfCasc == config.etaMinDefault) {
      config.etaMinTrackBachLfCasc.value = -config.etaMaxTrackBachLfCasc;
    }

    if (config.fillHistograms) {
      // general tracks
      registry.add("hRejTracks", "Tracks;;entries", {HistType::kTH1D, {{25, 0.5, 25.5}}});
      registry.add("hPtNoCuts", "all tracks;#it{p}_{T}^{track} (GeV/#it{c});entries", {HistType::kTH1D, {{360, 0., 36.}}});

      // 2-prong histograms
      registry.add("hPtCuts2Prong", "tracks selected for 2-prong vertexing;#it{p}_{T}^{track} (GeV/#it{c});entries", {HistType::kTH1D, {{360, 0., 36.}}});
      registry.add("hDCAToPrimXYVsPtCuts2Prong", "tracks selected for 2-prong vertexing;#it{p}_{T}^{track} (GeV/#it{c});DCAxy to prim. vtx. (cm);entries", {HistType::kTH2D, {{360, 0., 36.}, {400, -2., 2.}}});
      registry.add("hEtaCuts2Prong", "tracks selected for 2-prong vertexing;#it{#eta};entries", {HistType::kTH1D, {{static_cast<int>(0.6 * (config.etaMaxTrack2Prong - config.etaMinTrack2Prong) * 100), -1.2 * config.etaMinTrack2Prong, 1.2 * config.etaMaxTrack2Prong}}});
      // 3-prong histograms
      registry.add("hPtCuts3Prong", "tracks selected for 3-prong vertexing;#it{p}_{T}^{track} (GeV/#it{c});entries", {HistType::kTH1D, {{360, 0., 36.}}});
      registry.add("hDCAToPrimXYVsPtCuts3Prong", "tracks selected for 3-prong vertexing;#it{p}_{T}^{track} (GeV/#it{c});DCAxy to prim. vtx. (cm);entries", {HistType::kTH2D, {{360, 0., 36.}, {400, -2., 2.}}});
      registry.add("hEtaCuts3Prong", "tracks selected for 3-prong vertexing;#it{#eta};entries", {HistType::kTH1D, {{static_cast<int>(0.6 * (config.etaMaxTrack3Prong - config.etaMinTrack3Prong) * 100), -1.2 * config.etaMinTrack3Prong, 1.2 * config.etaMaxTrack3Prong}}});
      // bachelor (for V0 + bachelor decays) histograms
      registry.add("hPtCutsV0bachelor", "tracks selected for V0-bachelor vertexing;#it{p}_{T}^{track} (GeV/#it{c});entries", {HistType::kTH1D, {{360, 0., 36.}}});
      registry.add("hDCAToPrimXYVsPtCutsV0bachelor", "tracks selected for V0-bachelor vertexing;#it{p}_{T}^{track} (GeV/#it{c});DCAxy to prim. vtx. (cm);entries", {HistType::kTH2D, {{360, 0., 36.}, {400, -2., 2.}}});
      registry.add("hEtaCutsV0bachelor", "tracks selected for V0-bachelor vertexing;#it{#eta};entries", {HistType::kTH1D, {{static_cast<int>(0.6 * (config.etaMaxTrackBach - config.etaMinTrackBach) * 100), -1.2 * config.etaMinTrackBach, 1.2 * config.etaMaxTrackBach}}});
      // soft pion (for D*) histograms
      registry.add("hPtCutsSoftPionForDstar", "tracks selected for D* soft pion;#it{p}_{T}^{track} (GeV/#it{c});entries", {HistType::kTH1D, {{360, 0., 36.}}});
      registry.add("hDCAToPrimXYVsPtCutsSoftPionForDstar", "tracks selected for D* soft pion;#it{p}_{T}^{track} (GeV/#it{c});DCAxy to prim. vtx. (cm);entries", {HistType::kTH2D, {{360, 0., 36.}, {400, -2., 2.}}});
      registry.add("hEtaCutsSoftPionForDstar", "tracks selected for D* soft pion;#it{#eta};entries", {HistType::kTH1D, {{static_cast<int>(0.6 * (config.etaMaxSoftPionForDstar - config.etaMinSoftPionForDstar) * 100), -1.2 * config.etaMinSoftPionForDstar, 1.2 * config.etaMaxSoftPionForDstar}}});
      // bachelor (for cascade + bachelor decays) histograms
      registry.add("hPtCutsCascadeBachelor", "tracks selected for cascade-bachelor vertexing;#it{p}_{T}^{track} (GeV/#it{c});entries", {HistType::kTH1D, {{360, 0., 36.}}});
      registry.add("hDCAToPrimXYVsPtCutsCascadeBachelor", "tracks selected for cascade-bachelor vertexing;#it{p}_{T}^{track} (GeV/#it{c});DCAxy to prim. vtx. (cm);entries", {HistType::kTH2D, {{360, 0., 36.}, {400, -2., 2.}}});
      registry.add("hEtaCutsCascadeBachelor", "tracks selected for cascade-bachelor vertexing;#it{#eta};entries", {HistType::kTH1D, {{static_cast<int>(0.6 * (config.etaMaxTrackBachLfCasc - config.etaMinTrackBachLfCasc) * 100), -1.2 * config.etaMinTrackBachLfCasc, 1.2 * config.etaMaxTrackBachLfCasc}}});

      const std::string cutNames[nCuts + 1] = {"selected", "rej pT", "rej eta", "rej track quality", "rej dca"};
      const std::string candNames[CandidateType::NCandidateTypes] = {"2-prong", "3-prong", "bachelor", "dstar", "lfCascBachelor"};
      for (int iCandType = 0; iCandType < CandidateType::NCandidateTypes; iCandType++) {
        for (int iCut = 0; iCut < nCuts + 1; iCut++) {
          registry.get<TH1>(HIST("hRejTracks"))->GetXaxis()->SetBinLabel((nCuts + 1) * iCandType + iCut + 1, Form("%s %s", candNames[iCandType].data(), cutNames[iCut].data()));
        }
      }
    }

    // Needed for PV refitting
    if (config.doPvRefit) {
      if (config.fillHistograms) {
        const AxisSpec axisCollisionX{100, -20.f, 20.f, "X (cm)"};
        const AxisSpec axisCollisionY{100, -20.f, 20.f, "Y (cm)"};
        const AxisSpec axisCollisionZ{100, -20.f, 20.f, "Z (cm)"};
        const AxisSpec axisCollisionXOriginal{100, -2.f, 2.f, "X original PV (cm)"};
        const AxisSpec axisCollisionYOriginal{100, -2.f, 2.f, "Y original PV (cm)"};
        const AxisSpec axisCollisionZOriginal{100, -2.f, 2.f, "Z original PV (cm)"};
        const AxisSpec axisCollisionNContrib{1000, 0, 1000, "Number of contributors"};
        const AxisSpec axisCollisionDeltaX{axisPvRefitDeltaX, "#Delta x_{PV} (cm)"};
        const AxisSpec axisCollisionDeltaY{axisPvRefitDeltaY, "#Delta y_{PV} (cm)"};
        const AxisSpec axisCollisionDeltaZ{axisPvRefitDeltaZ, "#Delta z_{PV} (cm)"};

        registry.add("PvRefit/hVerticesPerTrack", "", kTH1D, {{3, 0.5f, 3.5f, ""}});
        registry.get<TH1>(HIST("PvRefit/hVerticesPerTrack"))->GetXaxis()->SetBinLabel(1, "All PV");
        registry.get<TH1>(HIST("PvRefit/hVerticesPerTrack"))->GetXaxis()->SetBinLabel(2, "PV refit doable");
        registry.get<TH1>(HIST("PvRefit/hVerticesPerTrack"))->GetXaxis()->SetBinLabel(3, "PV refit #chi^{2}!=-1");
        registry.add("PvRefit/hPvDeltaXvsNContrib", "", kTH2D, {axisCollisionNContrib, axisCollisionDeltaX});
        registry.add("PvRefit/hPvDeltaYvsNContrib", "", kTH2D, {axisCollisionNContrib, axisCollisionDeltaY});
        registry.add("PvRefit/hPvDeltaZvsNContrib", "", kTH2D, {axisCollisionNContrib, axisCollisionDeltaZ});
        registry.add("PvRefit/hChi2vsNContrib", "", kTH2D, {axisCollisionNContrib, {102, -1.5, 100.5, "#chi^{2} PV refit"}});
        registry.add("PvRefit/hPvRefitXChi2Minus1", "PV refit with #it{#chi}^{2}==#minus1", kTH2D, {axisCollisionX, axisCollisionXOriginal});
        registry.add("PvRefit/hPvRefitYChi2Minus1", "PV refit with #it{#chi}^{2}==#minus1", kTH2D, {axisCollisionY, axisCollisionYOriginal});
        registry.add("PvRefit/hPvRefitZChi2Minus1", "PV refit with #it{#chi}^{2}==#minus1", kTH2D, {axisCollisionZ, axisCollisionZOriginal});
        registry.add("PvRefit/hNContribPvRefitNotDoable", "N. contributors for PV refit not doable", kTH1D, {axisCollisionNContrib});
        registry.add("PvRefit/hNContribPvRefitChi2Minus1", "N. contributors original PV for PV refit #it{#chi}^{2}==#minus1", kTH1D, {axisCollisionNContrib});
      }

      ccdb->setURL(config.ccdbUrl);
      ccdb->setCaching(true);
      ccdb->setLocalObjectValidityChecking();

      lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>(config.ccdbPathLut));
      runNumber = 0;
    }

    // configure proton PID
    for (auto iChannel{0u}; iChannel < ChannelsProtonPid::NChannelsProtonPid; ++iChannel) {
      selectorProton[iChannel].setRangePtTpc(config.selectionsPid->get(iChannel, 0u), config.selectionsPid->get(iChannel, 1u));      // 0u == "minPtTpc", 1u == "maxPtTpc"
      selectorProton[iChannel].setRangePtTof(config.selectionsPid->get(iChannel, 3u), config.selectionsPid->get(iChannel, 4u));      // 3u == "minPtTof, 4u == "maxPtTof"
      selectorProton[iChannel].setRangeNSigmaTpc(-config.selectionsPid->get(iChannel, 2u), config.selectionsPid->get(iChannel, 2u)); // 2u == "nSigmaMaxTpc"
      selectorProton[iChannel].setRangeNSigmaTof(-config.selectionsPid->get(iChannel, 5u), config.selectionsPid->get(iChannel, 5u)); // 5u == "nSigmaMaxTof"
    }
    // after the proton PID we have the kaon PID
    selectorKaon.setRangePtTpc(config.selectionsPid->get(ChannelKaonPid, 0u), config.selectionsPid->get(ChannelKaonPid, 1u));      // 0u == "minPtTpc", 1u == "maxPtTpc"
    selectorKaon.setRangePtTof(config.selectionsPid->get(ChannelKaonPid, 3u), config.selectionsPid->get(ChannelKaonPid, 4u));      // 3u == "minPtTof, 4u == "maxPtTof"
    selectorKaon.setRangeNSigmaTpc(-config.selectionsPid->get(ChannelKaonPid, 2u), config.selectionsPid->get(ChannelKaonPid, 2u)); // 2u == "nSigmaMaxTpc"
    selectorKaon.setRangeNSigmaTof(-config.selectionsPid->get(ChannelKaonPid, 5u), config.selectionsPid->get(ChannelKaonPid, 5u)); // 5u == "nSigmaMaxTof"

    selectorDeuteron.setRangePtTpc(config.selectionsPid->get(ChannelsDeuteronPid, 0u), config.selectionsPid->get(ChannelsDeuteronPid, 1u));      // 0u == "minPtTpc", 1u == "maxPtTpc"
    selectorDeuteron.setRangePtTof(config.selectionsPid->get(ChannelsDeuteronPid, 3u), config.selectionsPid->get(ChannelsDeuteronPid, 4u));      // 3u == "minPtTof, 4u == "maxPtTof"
    selectorDeuteron.setRangeNSigmaTpc(-config.selectionsPid->get(ChannelsDeuteronPid, 2u), config.selectionsPid->get(ChannelsDeuteronPid, 2u)); // 2u == "nSigmaMaxTpc"
    selectorDeuteron.setRangeNSigmaTof(-config.selectionsPid->get(ChannelsDeuteronPid, 5u), config.selectionsPid->get(ChannelsDeuteronPid, 5u)); // 5u == "nSigmaMaxTof"
  }

  /// PID track cuts (for proton only)
  /// \param hfTrack is a track
  /// \return true if the track is compatible with a proton hypothesis
  template <int PidStrategy, typename T>
  uint8_t isSelectedPid(const T& hfTrack)
  {
    std::array statusPid{TrackSelectorPID::Accepted, TrackSelectorPID::Accepted, TrackSelectorPID::Accepted, TrackSelectorPID::Accepted, TrackSelectorPID::Accepted};
    if constexpr (PidStrategy == ProtonPidStrategy::PidTofOnly) {
      if (hfTrack.hasTOF()) {
        for (auto iChannel{0u}; iChannel < ChannelsProtonPid::NChannelsProtonPid; ++iChannel) {
          statusPid[iChannel] = selectorProton[iChannel].statusTof(hfTrack);
        }
        statusPid[ChannelKaonPid] = selectorKaon.statusTof(hfTrack);
        statusPid[ChannelsDeuteronPid] = selectorDeuteron.statusTof(hfTrack);
      }
    }
    if constexpr (PidStrategy == ProtonPidStrategy::PidTpcOnly) {
      if (hfTrack.hasTPC()) {
        for (auto iChannel{0u}; iChannel < ChannelsProtonPid::NChannelsProtonPid; ++iChannel) {
          statusPid[iChannel] = selectorProton[iChannel].statusTpc(hfTrack);
        }
        statusPid[ChannelKaonPid] = selectorKaon.statusTpc(hfTrack);
        statusPid[ChannelsDeuteronPid] = selectorDeuteron.statusTpc(hfTrack);
      }
    }
    if constexpr (PidStrategy == ProtonPidStrategy::PidTpcOrTof) {
      for (auto iChannel{0u}; iChannel < ChannelsProtonPid::NChannelsProtonPid; ++iChannel) {
        statusPid[iChannel] = selectorProton[iChannel].statusTpcOrTof(hfTrack);
      }
      statusPid[ChannelKaonPid] = selectorKaon.statusTpcOrTof(hfTrack);
      statusPid[ChannelsDeuteronPid] = selectorDeuteron.statusTpcOrTof(hfTrack);
    }
    if constexpr (PidStrategy == ProtonPidStrategy::PidTpcAndTof) {
      for (auto iChannel{0u}; iChannel < ChannelsProtonPid::NChannelsProtonPid; ++iChannel) {
        statusPid[iChannel] = selectorProton[iChannel].statusTpcAndTof(hfTrack);
      }
      statusPid[ChannelKaonPid] = selectorKaon.statusTpcAndTof(hfTrack);
      statusPid[ChannelsDeuteronPid] = selectorDeuteron.statusTpcAndTof(hfTrack);
    }

    int8_t flag = BIT(ChannelsProtonPid::NChannelsProtonPid + 2) - 1; // all bits on (including the kaon one)
    for (auto iChannel{0u}; iChannel < ChannelsProtonPid::NChannelsProtonPid + 2; ++iChannel) {
      if (statusPid[iChannel] == TrackSelectorPID::Rejected) {
        CLRBIT(flag, iChannel);
      }
    }

    return flag;
  }

  /// Single-track cuts for 2-prongs or 3-prongs
  /// \param hfTrack is a track
  /// \param trackPt is the track pt
  /// \param trackEta is the track eta
  /// \param dca is a 2-element array with dca in transverse and longitudinal directions
  /// \param statusProng is the selection flag
  template <typename T>
  void isSelectedTrack(const T& hfTrack, const float trackPt, const float trackEta, const std::array<float, 2>& dca, int& statusProng)
  {
    if (config.fillHistograms) {
      registry.fill(HIST("hPtNoCuts"), trackPt);
    }

    int iCut{2};
    // pT cut
    if (trackPt < config.ptMinTrack2Prong) {
      CLRBIT(statusProng, CandidateType::Cand2Prong); // set the nth bit to 0
      if (config.fillHistograms) {
        registry.fill(HIST("hRejTracks"), (nCuts + 1) * CandidateType::Cand2Prong + iCut);
      }
    }
    if (trackPt < config.ptMinTrack3Prong) {
      CLRBIT(statusProng, CandidateType::Cand3Prong);
      if (config.fillHistograms) {
        registry.fill(HIST("hRejTracks"), (nCuts + 1) * CandidateType::Cand3Prong + iCut);
      }
    }

    if (trackPt < config.ptMinTrackBach) {
      CLRBIT(statusProng, CandidateType::CandV0bachelor);
      if (config.fillHistograms) {
        registry.fill(HIST("hRejTracks"), (nCuts + 1) * CandidateType::CandV0bachelor + iCut);
      }
    }
    if (trackPt < config.ptMinSoftPionForDstar) {
      CLRBIT(statusProng, CandidateType::CandDstar);
      if (config.fillHistograms) {
        registry.fill(HIST("hRejTracks"), (nCuts + 1) * CandidateType::CandDstar + iCut);
      }
    }
    if (trackPt < config.ptMinTrackBachLfCasc) {
      CLRBIT(statusProng, CandidateType::CandCascadeBachelor);
      if (config.fillHistograms) {
        registry.fill(HIST("hRejTracks"), (nCuts + 1) * CandidateType::CandCascadeBachelor + iCut);
      }
    }

    iCut = 3;
    // eta cut
    if (TESTBIT(statusProng, CandidateType::Cand2Prong) && (trackEta > config.etaMaxTrack2Prong || trackEta < config.etaMinTrack2Prong)) {
      CLRBIT(statusProng, CandidateType::Cand2Prong);
      if (config.fillHistograms) {
        registry.fill(HIST("hRejTracks"), (nCuts + 1) * CandidateType::Cand2Prong + iCut);
      }
    }
    if (TESTBIT(statusProng, CandidateType::Cand3Prong) && (trackEta > config.etaMaxTrack3Prong || trackEta < config.etaMinTrack3Prong)) {
      CLRBIT(statusProng, CandidateType::Cand3Prong);
      if (config.fillHistograms) {
        registry.fill(HIST("hRejTracks"), (nCuts + 1) * CandidateType::Cand3Prong + iCut);
      }
    }

    if (TESTBIT(statusProng, CandidateType::CandV0bachelor) && (trackEta > config.etaMaxTrackBach || trackEta < config.etaMinTrackBach)) {
      CLRBIT(statusProng, CandidateType::CandV0bachelor);
      if (config.fillHistograms) {
        registry.fill(HIST("hRejTracks"), (nCuts + 1) * CandidateType::CandV0bachelor + iCut);
      }
    }

    if (TESTBIT(statusProng, CandidateType::CandDstar) && (trackEta > config.etaMaxSoftPionForDstar || trackEta < config.etaMinSoftPionForDstar)) {
      CLRBIT(statusProng, CandidateType::CandDstar);
      if (config.fillHistograms) {
        registry.fill(HIST("hRejTracks"), (nCuts + 1) * CandidateType::CandDstar + iCut);
      }
    }

    if (TESTBIT(statusProng, CandidateType::CandCascadeBachelor) && (trackEta > config.etaMaxTrackBachLfCasc || trackEta < config.etaMinTrackBachLfCasc)) {
      CLRBIT(statusProng, CandidateType::CandCascadeBachelor);
      if (config.fillHistograms) {
        registry.fill(HIST("hRejTracks"), (nCuts + 1) * CandidateType::CandCascadeBachelor + iCut);
      }
    }

    // quality cut
    iCut = 4;
    bool hasGoodQuality = true;
    if (config.doCutQuality.value && statusProng > 0) { // FIXME to make a more complete selection e.g track.flags() & o2::aod::track::TPCrefit && track.flags() & o2::aod::track::GoldenChi2 &&
      if (config.useIsGlobalTrack) {
        if (!hfTrack.isGlobalTrack()) {
          hasGoodQuality = false;
        }
      } else if (config.useIsGlobalTrackWoDCA) {
        if (!hfTrack.isGlobalTrackWoDCA()) {
          hasGoodQuality = false;
        }
      } else {
        const auto clustermap = hfTrack.itsClusterMap();
        if (!(hfTrack.tpcNClsFound() >= config.tpcNClsFoundMin.value && // is this the number of TPC clusters? It should not be used
              TESTBIT(hfTrack.flags(), o2::aod::track::ITSrefit) &&
              (TESTBIT(clustermap, 0) || TESTBIT(clustermap, 1)))) {
          hasGoodQuality = false;
        }
      }
      if (!hasGoodQuality) {
        for (int iCandType = 0; iCandType < CandidateType::NCandidateTypes; iCandType++) {
          if (iCandType == CandidateType::CandDstar) { // different quality criteria for D* soft pions
            continue;
          }
          CLRBIT(statusProng, iCandType);
          if (config.fillHistograms) {
            registry.fill(HIST("hRejTracks"), (nCuts + 1) * iCandType + iCut);
          }
        }
      }
    }

    // quality cut for soft pion
    hasGoodQuality = true;
    if (config.doCutQuality.value && TESTBIT(statusProng, CandidateType::CandDstar)) {
      if (config.useIsGlobalTrackForSoftPion) {
        if (!hfTrack.isGlobalTrack()) {
          hasGoodQuality = false;
        }
      } else if (config.useIsGlobalTrackWoDCAForSoftPion) {
        if (!hfTrack.isGlobalTrackWoDCA()) {
          hasGoodQuality = false;
        }
      } else if (config.useIsQualityTrackITSForSoftPion) {
        if (!hfTrack.isQualityTrackITS()) {
          hasGoodQuality = false;
        }
      } else { // selections for Run2 converted data
        const auto clustermap = hfTrack.itsClusterMap();
        if (!(TESTBIT(hfTrack.flags(), o2::aod::track::ITSrefit) && (TESTBIT(clustermap, 0) || TESTBIT(clustermap, 1)))) {
          hasGoodQuality = false;
        }
      }
      if (!hasGoodQuality) {
        CLRBIT(statusProng, CandidateType::CandDstar);
        if (config.fillHistograms) {
          registry.fill(HIST("hRejTracks"), (nCuts + 1) * CandidateType::CandDstar + iCut);
        }
      }
    }

    // quality cut for bachelor in cascade + bachelor decays
    hasGoodQuality = true;
    if (config.doCutQuality.value && TESTBIT(statusProng, CandidateType::CandCascadeBachelor)) {
      if (config.useIsGlobalTrackForBachLfCasc) {
        if (!hfTrack.isGlobalTrack()) {
          hasGoodQuality = false;
        }
      } else if (config.useIsGlobalTrackWoDCAForBachLfCasc) {
        if (!hfTrack.isGlobalTrackWoDCA()) {
          hasGoodQuality = false;
        }
      } else if (config.useIsQualityTrackITSForBachLfCasc) {
        if (!hfTrack.isQualityTrackITS()) {
          hasGoodQuality = false;
        }
      } else { // selections for Run2 converted data
        const auto clustermap = hfTrack.itsClusterMap();
        if (!(TESTBIT(hfTrack.flags(), o2::aod::track::ITSrefit) && (TESTBIT(clustermap, 0) || TESTBIT(clustermap, 1)))) {
          hasGoodQuality = false;
        }
      }
      if (!hasGoodQuality) {
        CLRBIT(statusProng, CandidateType::CandCascadeBachelor);
        if (config.fillHistograms) {
          registry.fill(HIST("hRejTracks"), (nCuts + 1) * CandidateType::CandCascadeBachelor + iCut);
        }
      }
    }

    // DCA cut
    iCut = 5;
    if (statusProng > 0) {
      for (int iCandType = 0; iCandType < CandidateType::NCandidateTypes; ++iCandType) {
        if (TESTBIT(statusProng, iCandType) && !isSelectedTrackDca(config.binsPtTrack, &cutsSingleTrack[iCandType], trackPt, dca[0], dca[1])) {
          CLRBIT(statusProng, iCandType);
          if (config.fillHistograms) {
            registry.fill(HIST("hRejTracks"), (nCuts + 1) * iCandType + iCut);
          }
        }
      }
    }

    // fill histograms
    if (config.fillHistograms) {
      iCut = 1;
      if (TESTBIT(statusProng, CandidateType::Cand2Prong)) {
        registry.fill(HIST("hPtCuts2Prong"), trackPt);
        registry.fill(HIST("hEtaCuts2Prong"), trackEta);
        registry.fill(HIST("hDCAToPrimXYVsPtCuts2Prong"), trackPt, dca[0]);
        registry.fill(HIST("hRejTracks"), (nCuts + 1) * CandidateType::Cand2Prong + iCut);
      }
      if (TESTBIT(statusProng, CandidateType::Cand3Prong)) {
        registry.fill(HIST("hPtCuts3Prong"), trackPt);
        registry.fill(HIST("hEtaCuts3Prong"), trackEta);
        registry.fill(HIST("hDCAToPrimXYVsPtCuts3Prong"), trackPt, dca[0]);
        registry.fill(HIST("hRejTracks"), (nCuts + 1) * CandidateType::Cand3Prong + iCut);
      }
      if (TESTBIT(statusProng, CandidateType::CandV0bachelor)) {
        registry.fill(HIST("hPtCutsV0bachelor"), trackPt);
        registry.fill(HIST("hEtaCutsV0bachelor"), trackEta);
        registry.fill(HIST("hDCAToPrimXYVsPtCutsV0bachelor"), trackPt, dca[0]);
        registry.fill(HIST("hRejTracks"), (nCuts + 1) * CandidateType::CandV0bachelor + iCut);
      }
      if (TESTBIT(statusProng, CandidateType::CandDstar)) {
        registry.fill(HIST("hPtCutsSoftPionForDstar"), trackPt);
        registry.fill(HIST("hEtaCutsSoftPionForDstar"), trackEta);
        registry.fill(HIST("hDCAToPrimXYVsPtCutsSoftPionForDstar"), trackPt, dca[0]);
        registry.fill(HIST("hRejTracks"), (nCuts + 1) * CandidateType::CandDstar + iCut);
      }
      if (TESTBIT(statusProng, CandidateType::CandCascadeBachelor)) {
        registry.fill(HIST("hPtCutsCascadeBachelor"), trackPt);
        registry.fill(HIST("hEtaCutsCascadeBachelor"), trackEta);
        registry.fill(HIST("hDCAToPrimXYVsPtCutsCascadeBachelor"), trackPt, dca[0]);
        registry.fill(HIST("hRejTracks"), (nCuts + 1) * CandidateType::CandCascadeBachelor + iCut);
      }
    }
  }

  /// Method for the PV refit and DCA recalculation for tracks with a collision assigned
  /// \param collision is a collision
  /// \param bcWithTimeStamps is a table of bunch crossing joined with timestamps used to query the CCDB for B and material budget
  /// \param vecPvContributorGlobId is a vector containing the global ID of PV contributors for the current collision
  /// \param vecPvContributorTrackParCov is a vector containing the TrackParCov of PV contributors for the current collision
  /// \param trackToRemove is the track to be removed, if contributor, from the PV refit
  /// \param pvCoord is an array containing the coordinates of the refitted PV
  /// \param pvCovMatrix is an array containing the covariance matrix values of the refitted PV
  /// \param dcaXYdcaZ is an array containing the dcaXY and dcaZ of trackToRemove with respect to the refitted PV
  template <typename TTrack>
  void performPvRefitTrack(aod::Collision const& collision,
                           aod::BCsWithTimestamps const&,
                           std::vector<int64_t> const& vecPvContributorGlobId,
                           std::vector<o2::track::TrackParCov> const& vecPvContributorTrackParCov,
                           TTrack const& trackToRemove,
                           std::array<float, 3>& pvCoord,
                           std::array<float, 6>& pvCovMatrix,
                           std::array<float, 2>& dcaXYdcaZ)
  {
    std::vector<bool> vecPvRefitContributorUsed(vecPvContributorGlobId.size(), true);

    /// Prepare the vertex refitting
    // set the magnetic field from CCDB
    const auto bc = collision.bc_as<o2::aod::BCsWithTimestamps>();
    initCCDB(bc, runNumber, ccdb, config.isRun2 ? config.ccdbPathGrp : config.ccdbPathGrpMag, lut, config.isRun2);
    /*if (runNumber != bc.runNumber()) {

      if (isRun2) { // Run 2 GRP object
        o2::parameters::GRPObject* grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>(ccdbPathGrp, bc.timestamp());
        if (grpo != nullptr) {
          o2::base::Propagator::initFieldFromGRP(grpo);
          o2::base::Propagator::Instance()->setMatLUT(lut);
          LOGF(info, "Setting magnetic field to %d kG for run %d from its GRP CCDB object (type o2::parameters::GRPObject)", grpo->getNominalL3Field(), bc.runNumber());
        } else {
          LOGF(fatal, "Run 2 GRP object (type o2::parameters::GRPObject) is not available in CCDB for run=%d at timestamp=%llu", bc.runNumber(), bc.timestamp());
        }
      } else { // Run 3 GRP object
        o2::parameters::GRPMagField* grpo = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(ccdbPathGrpMag, bc.timestamp());
        if (grpo != nullptr) {
          o2::base::Propagator::initFieldFromGRP(grpo);
          o2::base::Propagator::Instance()->setMatLUT(lut);
          LOG(info) << "Setting magnetic field to current" << grpo->getL3Current() << " A for run" << bc.runNumber() << " from its GRP CCDB object (type o2::parameters::GRPMagField)";
        } else {
          LOGF(fatal, "Run 3 GRP object (type o2::parameters::GRPMagField) is not available in CCDB for run=%d at timestamp=%llu", bc.runNumber(), bc.timestamp());
        }
      }

      runNumber = bc.runNumber();
    }*/

    // build the VertexBase to initialize the vertexer
    o2::dataformats::VertexBase primVtx;
    primVtx.setX(collision.posX());
    primVtx.setY(collision.posY());
    primVtx.setZ(collision.posZ());
    primVtx.setCov(collision.covXX(), collision.covXY(), collision.covYY(), collision.covXZ(), collision.covYZ(), collision.covZZ());
    // configure PVertexer
    o2::vertexing::PVertexer vertexer;
    o2::conf::ConfigurableParam::updateFromString("pvertexer.useMeanVertexConstraint=false"); /// remove diamond constraint (let's keep it at the moment...)
    vertexer.init();
    const bool pvRefitDoable = vertexer.prepareVertexRefit(vecPvContributorTrackParCov, primVtx);
    if (!pvRefitDoable) {
      LOG(info) << "Not enough tracks accepted for the refit";
      if (config.doPvRefit && config.fillHistograms) {
        registry.fill(HIST("PvRefit/hNContribPvRefitNotDoable"), collision.numContrib());
      }
    }
    if (config.debugPvRefit) {
      LOG(info) << "prepareVertexRefit = " << pvRefitDoable << " Ncontrib= " << vecPvContributorTrackParCov.size() << " Ntracks= " << collision.numContrib() << " Vtx= " << primVtx.asString();
    }

    if (config.fillHistograms) {
      registry.fill(HIST("PvRefit/hVerticesPerTrack"), 1);
      if (pvRefitDoable) {
        registry.fill(HIST("PvRefit/hVerticesPerTrack"), 2);
      }
    }
    /// PV refitting, if the tracks contributed to this at the beginning
    o2::dataformats::VertexBase primVtxBaseRecalc;
    bool recalcImpPar = false;
    if (config.doPvRefit && pvRefitDoable) {
      recalcImpPar = true;
      const auto trackIterator = std::find(vecPvContributorGlobId.begin(), vecPvContributorGlobId.end(), trackToRemove.globalIndex()); /// track global index
      if (trackIterator != vecPvContributorGlobId.end()) {

        /// this track contributed to the PV fit: let's do the refit without it
        const int entry = std::distance(vecPvContributorGlobId.begin(), trackIterator);

        vecPvRefitContributorUsed[entry] = false; /// remove the track from the PV refitting

        const auto primVtxRefitted = vertexer.refitVertex(vecPvRefitContributorUsed, primVtx); // vertex refit
        // LOG(info) << "refit " << cnt << "/" << ntr << " result = " << primVtxRefitted.asString();
        if (config.debugPvRefit) {
          LOG(info) << "refit for track with global index " << static_cast<int>(trackToRemove.globalIndex()) << " " << primVtxRefitted.asString();
        }
        if (primVtxRefitted.getChi2() < 0) {
          if (config.debugPvRefit) {
            LOG(info) << "---> Refitted vertex has bad chi2 = " << primVtxRefitted.getChi2();
          }
          if (config.fillHistograms) {
            registry.fill(HIST("PvRefit/hPvRefitXChi2Minus1"), primVtxRefitted.getX(), collision.posX());
            registry.fill(HIST("PvRefit/hPvRefitYChi2Minus1"), primVtxRefitted.getY(), collision.posY());
            registry.fill(HIST("PvRefit/hPvRefitZChi2Minus1"), primVtxRefitted.getZ(), collision.posZ());
            registry.fill(HIST("PvRefit/hNContribPvRefitChi2Minus1"), collision.numContrib());
          }
          recalcImpPar = false;
        } else if (config.fillHistograms) {
          registry.fill(HIST("PvRefit/hVerticesPerTrack"), 3);
        }
        if (config.fillHistograms) {
          registry.fill(HIST("PvRefit/hChi2vsNContrib"), primVtxRefitted.getNContributors(), primVtxRefitted.getChi2());
        }

        vecPvRefitContributorUsed[entry] = true; /// restore the track for the next PV refitting (probably not necessary here)

        if (recalcImpPar) {
          // fill the histograms for refitted PV with good Chi2
          const double deltaX = primVtx.getX() - primVtxRefitted.getX();
          const double deltaY = primVtx.getY() - primVtxRefitted.getY();
          const double deltaZ = primVtx.getZ() - primVtxRefitted.getZ();
          if (config.fillHistograms) {
            registry.fill(HIST("PvRefit/hPvDeltaXvsNContrib"), primVtxRefitted.getNContributors(), deltaX);
            registry.fill(HIST("PvRefit/hPvDeltaYvsNContrib"), primVtxRefitted.getNContributors(), deltaY);
            registry.fill(HIST("PvRefit/hPvDeltaZvsNContrib"), primVtxRefitted.getNContributors(), deltaZ);
          }

          // fill the newly calculated PV
          primVtxBaseRecalc.setX(primVtxRefitted.getX());
          primVtxBaseRecalc.setY(primVtxRefitted.getY());
          primVtxBaseRecalc.setZ(primVtxRefitted.getZ());
          primVtxBaseRecalc.setCov(primVtxRefitted.getSigmaX2(), primVtxRefitted.getSigmaXY(), primVtxRefitted.getSigmaY2(), primVtxRefitted.getSigmaXZ(), primVtxRefitted.getSigmaYZ(), primVtxRefitted.getSigmaZ2());
        }

        // cnt++;
      }
    } /// end 'if (doPvRefit && pvRefitDoable)'

    // updated value after PV recalculation
    if (recalcImpPar) {

      /// Track propagation to the PV refit considering also the material budget
      /// Mandatory for tracks updated at most only to the innermost ITS layer
      auto trackPar = getTrackPar(trackToRemove);
      std::array dcaInfo{-999.f, -999.f};
      if (o2::base::Propagator::Instance()->propagateToDCABxByBz({primVtxBaseRecalc.getX(), primVtxBaseRecalc.getY(), primVtxBaseRecalc.getZ()}, trackPar, 2.f, noMatCorr, &dcaInfo)) {
        pvCoord[0] = primVtxBaseRecalc.getX();
        pvCoord[1] = primVtxBaseRecalc.getY();
        pvCoord[2] = primVtxBaseRecalc.getZ();
        pvCovMatrix[0] = primVtxBaseRecalc.getSigmaX2();
        pvCovMatrix[1] = primVtxBaseRecalc.getSigmaXY();
        pvCovMatrix[2] = primVtxBaseRecalc.getSigmaY2();
        pvCovMatrix[3] = primVtxBaseRecalc.getSigmaXZ();
        pvCovMatrix[4] = primVtxBaseRecalc.getSigmaYZ();
        pvCovMatrix[5] = primVtxBaseRecalc.getSigmaZ2();
        dcaXYdcaZ[0] = dcaInfo[0]; // [cm]
        dcaXYdcaZ[1] = dcaInfo[1]; // [cm]
        // TODO: add DCAxy and DCAz uncertainties?
      }

      /// Track propagation to the PV refit done only ia geometrical way
      /// Correct only if no further material budget is crossed, namely for tracks already propagated to the original PV
      // o2::dataformats::DCA impactParameter;
      // if (getTrackParCov(myTrack).propagateToDCA(primVtxBaseRecalc, o2::base::Propagator::Instance()->getNominalBz(), &impactParameter)) {
      //   if (debug) {
      //     LOG(info) << "===> nominal Bz: " << o2::base::Propagator::Instance()->getNominalBz();
      //   }
      //   pvCoord[0] = primVtxBaseRecalc.getX();
      //   pvCoord[1] = primVtxBaseRecalc.getY();
      //   pvCoord[2] = primVtxBaseRecalc.getZ();
      //   pvCovMatrix[0] = primVtxBaseRecalc.getSigmaX2();
      //   pvCovMatrix[1] = primVtxBaseRecalc.getSigmaXY();
      //   pvCovMatrix[2] = primVtxBaseRecalc.getSigmaY2();
      //   pvCovMatrix[3] = primVtxBaseRecalc.getSigmaXZ();
      //   pvCovMatrix[4] = primVtxBaseRecalc.getSigmaYZ();
      //   pvCovMatrix[5] = primVtxBaseRecalc.getSigmaZ2();
      //   dcaXYdcaZ[0] = impactParameter.getY(); // [cm]
      //   dcaXYdcaZ[1] = impactParameter.getZ();    // [cm]
      //   // TODO: add DCAxy and DCAz uncertainties?
      // }
    }
  } /// end of performPvRefitTrack function

  /// Selection tag for tracks
  /// \tparam TTracks is the type of the track table
  /// \param collision is the collision iterator
  /// \param trackIndicesCollision are the track indices associated to this collision (from track-to-collision-associator)
  /// \param pvContrCollision are the PV contributors of this collision
  /// \param bcWithTimeStamps is the bc with timestamp for PVrefit
  /// \param pvRefitDcaPerTrack is a vector to be filled with track dcas after PV refit
  /// \param pvRefitPvCoordPerTrack is a vector to be filled with PV coordinates after PV refit
  /// \param pvRefitPvCovMatrixPerTrack is a vector to be filled with PV coordinate covariances after PV refit
  template <int PidStrategy, typename TTracks, typename GroupedTrackIndices, typename GroupedPvContributors>
  void runTagSelTracks(aod::Collision const& collision,
                       TTracks const&,
                       GroupedTrackIndices const& trackIndicesCollision,
                       GroupedPvContributors const& pvContrCollision,
                       aod::BCsWithTimestamps const& bcWithTimeStamps,
                       std::vector<std::array<float, 2>>& pvRefitDcaPerTrack,
                       std::vector<std::array<float, 3>>& pvRefitPvCoordPerTrack,
                       std::vector<std::array<float, 6>>& pvRefitPvCovMatrixPerTrack)
  {
    const auto thisCollId = collision.globalIndex();
    for (const auto& trackId : trackIndicesCollision) {
      int statusProng = BIT(CandidateType::NCandidateTypes) - 1; // all bits on
      const auto track = trackId.template track_as<TTracks>();
      float trackPt = track.pt();
      float trackEta = track.eta();

      std::array pvRefitDcaXYDcaZ{track.dcaXY(), track.dcaZ()};
      std::array pvRefitPvCoord{0.f, 0.f, 0.f};
      std::array pvRefitPvCovMatrix{1e10f, 1e10f, 1e10f, 1e10f, 1e10f, 1e10f};

      // PV refit and DCA recalculation only for tracks with an assigned collision
      if (config.doPvRefit && track.has_collision() && track.collisionId() == thisCollId && track.isPVContributor()) {
        pvRefitPvCoord = {collision.posX(), collision.posY(), collision.posZ()};
        pvRefitPvCovMatrix = {collision.covXX(), collision.covXY(), collision.covYY(), collision.covXZ(), collision.covYZ(), collision.covZZ()};

        /// retrieve PV contributors for the current collision
        std::vector<int64_t> vecPvContributorGlobId{};
        std::vector<o2::track::TrackParCov> vecPvContributorTrackParCov{};

        for (const auto& contributor : pvContrCollision) {
          vecPvContributorGlobId.push_back(contributor.globalIndex());
          vecPvContributorTrackParCov.push_back(getTrackParCov(contributor));
        }
        if (config.debugPvRefit) {
          LOG(info) << "### vecPvContributorGlobId.size()=" << vecPvContributorGlobId.size() << ", vecPvContributorTrackParCov.size()=" << vecPvContributorTrackParCov.size() << ", N. original contributors=" << collision.numContrib();
          /// Perform the PV refit only for tracks with an assigned collision
          LOG(info) << "[BEFORE performPvRefitTrack] track.collision().globalIndex(): " << collision.globalIndex();
        }
        performPvRefitTrack(collision, bcWithTimeStamps, vecPvContributorGlobId, vecPvContributorTrackParCov, track, pvRefitPvCoord, pvRefitPvCovMatrix, pvRefitDcaXYDcaZ);
        // we subtract the offset since trackIdx is the global index referred to the total track table
        const auto trackIdx = track.globalIndex();
        pvRefitDcaPerTrack[trackIdx] = pvRefitDcaXYDcaZ;
        pvRefitPvCoordPerTrack[trackIdx] = pvRefitPvCoord;
        pvRefitPvCovMatrixPerTrack[trackIdx] = pvRefitPvCovMatrix;
      } else if (track.collisionId() != thisCollId) {
        const auto bc = collision.bc_as<o2::aod::BCsWithTimestamps>();
        initCCDB(bc, runNumber, ccdb, config.isRun2 ? config.ccdbPathGrp : config.ccdbPathGrpMag, lut, config.isRun2);
        auto trackPar = getTrackPar(track);
        std::array dcaInfo{-999.f, -999.f};
        o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackPar, 2.f, noMatCorr, &dcaInfo);
        trackPt = trackPar.getPt();
        trackEta = trackPar.getEta();
        pvRefitDcaXYDcaZ[0] = dcaInfo[0];
        pvRefitDcaXYDcaZ[1] = dcaInfo[1];
      }

      // bool cutStatus[CandidateType::NCandidateTypes][nCuts];
      // if (debug) {
      //   for (int iCandType = 0; iCandType < CandidateType::NCandidateTypes; iCandType++) {
      //     for (int iCut = 0; iCut < nCuts; iCut++) {
      //       cutStatus[iCandType][iCut] = true;
      //     }
      //   }
      // }

      isSelectedTrack(track, trackPt, trackEta, pvRefitDcaXYDcaZ, statusProng);
      const int8_t isIdentifiedPid = isSelectedPid<PidStrategy>(track);
      const bool isPositive = track.sign() > 0;
      rowSelectedTrack(statusProng, isIdentifiedPid, isPositive);
    }
  }

  /// Helper function to fill PVrefit table
  /// \param pvRefitDcaPerTrack is a vector filled with track dcas after PV refit
  /// \param pvRefitPvCoordPerTrack is a vector filled with PV coordinates after PV refit
  /// \param pvRefitPvCovMatrixPerTrack is a vector filled with PV coordinate covariances after PV refit
  void fillPvRefitTable(std::vector<std::array<float, 2>> const& pvRefitDcaPerTrack,
                        std::vector<std::array<float, 3>> const& pvRefitPvCoordPerTrack,
                        std::vector<std::array<float, 6>> const& pvRefitPvCovMatrixPerTrack)
  {
    for (auto iTrack{0u}; iTrack < pvRefitDcaPerTrack.size(); ++iTrack) {
      tabPvRefitTrack(pvRefitPvCoordPerTrack[iTrack][0], pvRefitPvCoordPerTrack[iTrack][1], pvRefitPvCoordPerTrack[iTrack][2],
                      pvRefitPvCovMatrixPerTrack[iTrack][0], pvRefitPvCovMatrixPerTrack[iTrack][1], pvRefitPvCovMatrixPerTrack[iTrack][2], pvRefitPvCovMatrixPerTrack[iTrack][3], pvRefitPvCovMatrixPerTrack[iTrack][4], pvRefitPvCovMatrixPerTrack[iTrack][5],
                      pvRefitDcaPerTrack[iTrack][0], pvRefitDcaPerTrack[iTrack][1]);
    }
  }

  void processNoPid(aod::Collisions const& collisions,
                    TrackAssoc const& trackIndices,
                    TracksWithSelAndDca const& tracks,
                    aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    rowSelectedTrack.reserve(tracks.size());
    // prepare vectors to cache quantities needed for PV refit
    std::vector<std::array<float, 2>> pvRefitDcaPerTrack{};
    std::vector<std::array<float, 3>> pvRefitPvCoordPerTrack{};
    std::vector<std::array<float, 6>> pvRefitPvCovMatrixPerTrack{};
    if (config.doPvRefit) {
      const auto numTracks = tracks.size();
      pvRefitDcaPerTrack.resize(numTracks);
      pvRefitPvCoordPerTrack.resize(numTracks);
      pvRefitPvCovMatrixPerTrack.resize(numTracks);
      tabPvRefitTrack.reserve(numTracks);
    }

    for (const auto& collision : collisions) {
      const auto thisCollId = collision.globalIndex();
      const auto groupedTrackIndices = trackIndices.sliceBy(trackIndicesPerCollision, thisCollId);
      const auto pvContrCollision = pvContributors->sliceByCached(aod::track::collisionId, thisCollId, cache);
      runTagSelTracks<NoPid>(collision, tracks, groupedTrackIndices, pvContrCollision, bcWithTimeStamps, pvRefitDcaPerTrack, pvRefitPvCoordPerTrack, pvRefitPvCovMatrixPerTrack);
    }

    if (config.doPvRefit) { /// fill table with PV refit info (it has to be filled per track and not track index)
      fillPvRefitTable(pvRefitDcaPerTrack, pvRefitPvCoordPerTrack, pvRefitPvCovMatrixPerTrack);
    }
  }
  PROCESS_SWITCH(HfTrackIndexSkimCreatorTagSelTracks, processNoPid, "Process without PID selections", true);

  void processProtonPidTpc(aod::Collisions const& collisions,
                           TrackAssoc const& trackIndices,
                           TracksWithSelAndDcaAndPidTpc const& tracks,
                           aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    rowSelectedTrack.reserve(tracks.size());
    // prepare vectors to cache quantities needed for PV refit
    std::vector<std::array<float, 2>> pvRefitDcaPerTrack{};
    std::vector<std::array<float, 3>> pvRefitPvCoordPerTrack{};
    std::vector<std::array<float, 6>> pvRefitPvCovMatrixPerTrack{};
    if (config.doPvRefit) {
      const auto numTracks = tracks.size();
      pvRefitDcaPerTrack.resize(numTracks);
      pvRefitPvCoordPerTrack.resize(numTracks);
      pvRefitPvCovMatrixPerTrack.resize(numTracks);
      tabPvRefitTrack.reserve(numTracks);
    }

    for (const auto& collision : collisions) {
      const auto thisCollId = collision.globalIndex();
      const auto groupedTrackIndices = trackIndices.sliceBy(trackIndicesPerCollision, thisCollId);
      const auto pvContrCollision = pvContributorsWithPidTpc->sliceByCached(aod::track::collisionId, thisCollId, cache);
      runTagSelTracks<PidTpcOnly>(collision, tracks, groupedTrackIndices, pvContrCollision, bcWithTimeStamps, pvRefitDcaPerTrack, pvRefitPvCoordPerTrack, pvRefitPvCovMatrixPerTrack);
    }

    if (config.doPvRefit) { /// fill table with PV refit info (it has to be filled per track and not track index)
      fillPvRefitTable(pvRefitDcaPerTrack, pvRefitPvCoordPerTrack, pvRefitPvCovMatrixPerTrack);
    }
  }
  PROCESS_SWITCH(HfTrackIndexSkimCreatorTagSelTracks, processProtonPidTpc, "Process with proton TPC PID selection", false);

  void processProtonPidTof(aod::Collisions const& collisions,
                           TrackAssoc const& trackIndices,
                           TracksWithSelAndDcaAndPidTof const& tracks,
                           aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    rowSelectedTrack.reserve(tracks.size());
    // prepare vectors to cache quantities needed for PV refit
    std::vector<std::array<float, 2>> pvRefitDcaPerTrack{};
    std::vector<std::array<float, 3>> pvRefitPvCoordPerTrack{};
    std::vector<std::array<float, 6>> pvRefitPvCovMatrixPerTrack{};
    if (config.doPvRefit) {
      const auto numTracks = tracks.size();
      pvRefitDcaPerTrack.resize(numTracks);
      pvRefitPvCoordPerTrack.resize(numTracks);
      pvRefitPvCovMatrixPerTrack.resize(numTracks);
      tabPvRefitTrack.reserve(numTracks);
    }

    for (const auto& collision : collisions) {
      const auto thisCollId = collision.globalIndex();
      const auto groupedTrackIndices = trackIndices.sliceBy(trackIndicesPerCollision, thisCollId);
      const auto pvContrCollision = pvContributorsWithPidTof->sliceByCached(aod::track::collisionId, thisCollId, cache);
      runTagSelTracks<PidTofOnly>(collision, tracks, groupedTrackIndices, pvContrCollision, bcWithTimeStamps, pvRefitDcaPerTrack, pvRefitPvCoordPerTrack, pvRefitPvCovMatrixPerTrack);
    }

    if (config.doPvRefit) { /// fill table with PV refit info (it has to be filled per track and not track index)
      fillPvRefitTable(pvRefitDcaPerTrack, pvRefitPvCoordPerTrack, pvRefitPvCovMatrixPerTrack);
    }
  }
  PROCESS_SWITCH(HfTrackIndexSkimCreatorTagSelTracks, processProtonPidTof, "Process with proton TOF PID selection", false);

  void processProtonPidTpcOrTof(aod::Collisions const& collisions,
                                TrackAssoc const& trackIndices,
                                TracksWithSelAndDcaAndPidTpcTof const& tracks,
                                aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    rowSelectedTrack.reserve(tracks.size());
    // prepare vectors to cache quantities needed for PV refit
    std::vector<std::array<float, 2>> pvRefitDcaPerTrack{};
    std::vector<std::array<float, 3>> pvRefitPvCoordPerTrack{};
    std::vector<std::array<float, 6>> pvRefitPvCovMatrixPerTrack{};
    if (config.doPvRefit) {
      const auto numTracks = tracks.size();
      pvRefitDcaPerTrack.resize(numTracks);
      pvRefitPvCoordPerTrack.resize(numTracks);
      pvRefitPvCovMatrixPerTrack.resize(numTracks);
      tabPvRefitTrack.reserve(numTracks);
    }

    for (const auto& collision : collisions) {
      const auto thisCollId = collision.globalIndex();
      const auto groupedTrackIndices = trackIndices.sliceBy(trackIndicesPerCollision, thisCollId);
      const auto pvContrCollision = pvContributorsWithPidTpcTof->sliceByCached(aod::track::collisionId, thisCollId, cache);
      runTagSelTracks<PidTpcOrTof>(collision, tracks, groupedTrackIndices, pvContrCollision, bcWithTimeStamps, pvRefitDcaPerTrack, pvRefitPvCoordPerTrack, pvRefitPvCovMatrixPerTrack);
    }

    if (config.doPvRefit) { /// fill table with PV refit info (it has to be filled per track and not track index)
      fillPvRefitTable(pvRefitDcaPerTrack, pvRefitPvCoordPerTrack, pvRefitPvCovMatrixPerTrack);
    }
  }
  PROCESS_SWITCH(HfTrackIndexSkimCreatorTagSelTracks, processProtonPidTpcOrTof, "Process with proton PID selection (TPC or TOF logic)", false);

  void processProtonPidTpcAndTof(aod::Collisions const& collisions,
                                 TrackAssoc const& trackIndices,
                                 TracksWithSelAndDcaAndPidTpcTof const& tracks,
                                 aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    rowSelectedTrack.reserve(tracks.size());
    // prepare vectors to cache quantities needed for PV refit
    std::vector<std::array<float, 2>> pvRefitDcaPerTrack{};
    std::vector<std::array<float, 3>> pvRefitPvCoordPerTrack{};
    std::vector<std::array<float, 6>> pvRefitPvCovMatrixPerTrack{};
    if (config.doPvRefit) {
      const auto numTracks = tracks.size();
      pvRefitDcaPerTrack.resize(numTracks);
      pvRefitPvCoordPerTrack.resize(numTracks);
      pvRefitPvCovMatrixPerTrack.resize(numTracks);
      tabPvRefitTrack.reserve(numTracks);
    }

    for (const auto& collision : collisions) {
      const auto thisCollId = collision.globalIndex();
      const auto groupedTrackIndices = trackIndices.sliceBy(trackIndicesPerCollision, thisCollId);
      const auto pvContrCollision = pvContributorsWithPidTpcTof->sliceByCached(aod::track::collisionId, thisCollId, cache);
      runTagSelTracks<PidTpcAndTof>(collision, tracks, groupedTrackIndices, pvContrCollision, bcWithTimeStamps, pvRefitDcaPerTrack, pvRefitPvCoordPerTrack, pvRefitPvCovMatrixPerTrack);
    }

    if (config.doPvRefit) { /// fill table with PV refit info (it has to be filled per track and not track index)
      fillPvRefitTable(pvRefitDcaPerTrack, pvRefitPvCoordPerTrack, pvRefitPvCovMatrixPerTrack);
    }
  }
  PROCESS_SWITCH(HfTrackIndexSkimCreatorTagSelTracks, processProtonPidTpcAndTof, "Process with proton PID selection (TPC and TOF logic)", false);
};

//____________________________________________________________________________________________________________________________________________

/// Pre-selection of 2-prong and 3-prong secondary vertices
struct HfTrackIndexSkimCreator {
  Produces<aod::Hf2Prongs> rowTrackIndexProng2;
  Produces<aod::HfCutStatus2Prong> rowProng2CutStatus;
  Produces<aod::HfPvRefit2Prong> rowProng2PVrefit;
  Produces<aod::Hf3Prongs> rowTrackIndexProng3;
  Produces<aod::HfCutStatus3Prong> rowProng3CutStatus;
  Produces<aod::HfPvRefit3Prong> rowProng3PVrefit;
  Produces<aod::HfDstars> rowTrackIndexDstar;
  Produces<aod::HfCutStatusDstar> rowDstarCutStatus;
  Produces<aod::HfPvRefitDstar> rowDstarPVrefit;
  // Tables with ML scores for HF Filters
  Produces<aod::Hf2ProngMlProbs> rowTrackIndexMlScoreProng2;
  Produces<aod::Hf3ProngMlProbs> rowTrackIndexMlScoreProng3;

  struct : ConfigurableGroup {
    Configurable<bool> isRun2{"isRun2", false, "enable Run 2 or Run 3 GRP objects for magnetic field"};
    Configurable<bool> do3Prong{"do3Prong", false, "do 3 prong"};
    Configurable<bool> doDstar{"doDstar", false, "do D* candidates"};
    Configurable<bool> debug{"debug", false, "debug mode"};
    Configurable<bool> debugPvRefit{"debugPvRefit", false, "debug lines for primary vertex refit"};
    Configurable<bool> fillHistograms{"fillHistograms", true, "fill histograms"};
    // Configurable<int> nCollsMax{"nCollsMax", -1, "Max collisions per file"}; //can be added to run over limited collisions per file - for tesing purposes
    // preselection
    Configurable<double> ptTolerance{"ptTolerance", 0.1, "pT tolerance in GeV/c for applying preselections before vertex reconstruction"};
    // preselection of 3-prongs using the decay length computed only with the first two tracks
    Configurable<double> minTwoTrackDecayLengthFor3Prongs{"minTwoTrackDecayLengthFor3Prongs", 0., "Minimum decay length computed with 2 tracks for 3-prongs to speedup combinatorial"};
    Configurable<double> maxTwoTrackChi2PcaFor3Prongs{"maxTwoTrackChi2PcaFor3Prongs", 1.e10, "Maximum chi2 pca computed with 2 tracks for 3-prongs to speedup combinatorial"};
    // vertexing
    // Configurable<double> bz{"bz", 5., "magnetic field kG"};
    Configurable<bool> propagateToPCA{"propagateToPCA", true, "create tracks version propagated to PCA"};
    Configurable<bool> useAbsDCA{"useAbsDCA", false, "Minimise abs. distance rather than chi2"};
    Configurable<bool> useWeightedFinalPCA{"useWeightedFinalPCA", false, "Recalculate vertex position using track covariances, effective only if useAbsDCA is true"};
    Configurable<double> maxR{"maxR", 200., "reject PCA's above this radius"};
    Configurable<double> maxDZIni{"maxDZIni", 4., "reject (if>0) PCA candidate if tracks DZ exceeds threshold"};
    Configurable<double> minParamChange{"minParamChange", 1.e-3, "stop iterations if largest change of any X is smaller than this"};
    Configurable<double> minRelChi2Change{"minRelChi2Change", 0.9, "stop iterations if chi2/chi2old > this"};
    // CCDB
    Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
    Configurable<std::string> ccdbPathLut{"ccdbPathLut", "GLO/Param/MatLUT", "Path for LUT parametrization"};
    Configurable<std::string> ccdbPathGrp{"ccdbPathGrp", "GLO/GRP/GRP", "Path of the grp file (Run 2)"};
    Configurable<std::string> ccdbPathGrpMag{"ccdbPathGrpMag", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object (Run 3)"};

    // D0 cuts
    Configurable<std::vector<double>> binsPtD0ToPiK{"binsPtD0ToPiK", std::vector<double>{hf_cuts_presel_2prong::vecBinsPt}, "pT bin limits for D0->piK pT-dependent cuts"};
    Configurable<LabeledArray<double>> cutsD0ToPiK{"cutsD0ToPiK", {hf_cuts_presel_2prong::Cuts[0], hf_cuts_presel_2prong::NBinsPt, hf_cuts_presel_2prong::NCutVars, hf_cuts_presel_2prong::labelsPt, hf_cuts_presel_2prong::labelsCutVar}, "D0->piK selections per pT bin"};
    // Jpsi -> ee cuts
    Configurable<std::vector<double>> binsPtJpsiToEE{"binsPtJpsiToEE", std::vector<double>{hf_cuts_presel_2prong::vecBinsPt}, "pT bin limits for Jpsi->ee pT-dependent cuts"};
    Configurable<LabeledArray<double>> cutsJpsiToEE{"cutsJpsiToEE", {hf_cuts_presel_2prong::Cuts[0], hf_cuts_presel_2prong::NBinsPt, hf_cuts_presel_2prong::NCutVars, hf_cuts_presel_2prong::labelsPt, hf_cuts_presel_2prong::labelsCutVar}, "Jpsi->ee selections per pT bin"};
    // Jpsi -> mumu cuts
    Configurable<std::vector<double>> binsPtJpsiToMuMu{"binsPtJpsiToMuMu", std::vector<double>{hf_cuts_presel_2prong::vecBinsPt}, "pT bin limits for Jpsi->mumu pT-dependent cuts"};
    Configurable<LabeledArray<double>> cutsJpsiToMuMu{"cutsJpsiToMuMu", {hf_cuts_presel_2prong::Cuts[0], hf_cuts_presel_2prong::NBinsPt, hf_cuts_presel_2prong::NCutVars, hf_cuts_presel_2prong::labelsPt, hf_cuts_presel_2prong::labelsCutVar}, "Jpsi->mumu selections per pT bin"};
    // D+ cuts
    Configurable<std::vector<double>> binsPtDplusToPiKPi{"binsPtDplusToPiKPi", std::vector<double>{hf_cuts_presel_3prong::vecBinsPt}, "pT bin limits for D+->piKpi pT-dependent cuts"};
    Configurable<LabeledArray<double>> cutsDplusToPiKPi{"cutsDplusToPiKPi", {hf_cuts_presel_3prong::Cuts[0], hf_cuts_presel_3prong::NBinsPt, hf_cuts_presel_3prong::NCutVars, hf_cuts_presel_3prong::labelsPt, hf_cuts_presel_3prong::labelsCutVar}, "D+->piKpi selections per pT bin"};
    // Ds+ cuts
    Configurable<std::vector<double>> binsPtDsToKKPi{"binsPtDsToKKPi", std::vector<double>{hf_cuts_presel_ds::vecBinsPt}, "pT bin limits for Ds+->KKPi pT-dependent cuts"};
    Configurable<LabeledArray<double>> cutsDsToKKPi{"cutsDsToKKPi", {hf_cuts_presel_ds::Cuts[0], hf_cuts_presel_ds::NBinsPt, hf_cuts_presel_ds::NCutVars, hf_cuts_presel_ds::labelsPt, hf_cuts_presel_ds::labelsCutVar}, "Ds+->KKPi selections per pT bin"};
    // Lc+ cuts
    Configurable<std::vector<double>> binsPtLcToPKPi{"binsPtLcToPKPi", std::vector<double>{hf_cuts_presel_3prong::vecBinsPt}, "pT bin limits for Lc->pKpi pT-dependent cuts"};
    Configurable<LabeledArray<double>> cutsLcToPKPi{"cutsLcToPKPi", {hf_cuts_presel_3prong::Cuts[0], hf_cuts_presel_3prong::NBinsPt, hf_cuts_presel_3prong::NCutVars, hf_cuts_presel_3prong::labelsPt, hf_cuts_presel_3prong::labelsCutVar}, "Lc->pKpi selections per pT bin"};
    // Xic+ cuts
    Configurable<std::vector<double>> binsPtXicToPKPi{"binsPtXicToPKPi", std::vector<double>{hf_cuts_presel_3prong::vecBinsPt}, "pT bin limits for Xic->pKpi pT-dependent cuts"};
    Configurable<LabeledArray<double>> cutsXicToPKPi{"cutsXicToPKPi", {hf_cuts_presel_3prong::Cuts[0], hf_cuts_presel_3prong::NBinsPt, hf_cuts_presel_3prong::NCutVars, hf_cuts_presel_3prong::labelsPt, hf_cuts_presel_3prong::labelsCutVar}, "Xic->pKpi selections per pT bin"};
    // Cd cuts
    Configurable<std::vector<double>> binsPtCdToDeKPi{"binsPtCdToDeKPi", std::vector<double>{hf_cuts_presel_3prong::vecBinsPt}, "pT bin limits for Cd->DeKpi pT-dependent cuts"};
    Configurable<LabeledArray<double>> cutsCdToDeKPi{"cutsCdToDeKPi", {hf_cuts_presel_3prong::Cuts[0], hf_cuts_presel_3prong::NBinsPt, hf_cuts_presel_3prong::NCutVars, hf_cuts_presel_3prong::labelsPt, hf_cuts_presel_3prong::labelsCutVar}, "Cd->deKpi selections per pT bin"};

    // D*+ cuts
    Configurable<std::vector<double>> binsPtDstarToD0Pi{"binsPtDstarToD0Pi", std::vector<double>{hf_cuts_presel_dstar::vecBinsPt}, "pT bin limits for D*+->D0pi pT-dependent cuts"};
    Configurable<LabeledArray<double>> cutsDstarToD0Pi{"cutsDstarToD0Pi", {hf_cuts_presel_dstar::Cuts[0], hf_cuts_presel_dstar::NBinsPt, hf_cuts_presel_dstar::NCutVars, hf_cuts_presel_dstar::labelsPt, hf_cuts_presel_dstar::labelsCutVar}, "D*+->D0pi selections per pT bin"};

    // proton PID selections for Lc and Xic
    Configurable<bool> applyProtonPidForLcToPKPi{"applyProtonPidForLcToPKPi", false, "Apply proton PID for Lc->pKpi"};
    Configurable<bool> applyProtonPidForXicToPKPi{"applyProtonPidForXicToPKPi", false, "Apply proton PID for Xic->pKpi"};
    Configurable<bool> applyKaonPidIn3Prongs{"applyKaonPidIn3Prongs", false, "Apply kaon PID for opposite-sign track in 3-prong and D* decays"};
    Configurable<bool> applyDeuteronPidForCdToDeKPi{"applyDeuteronPidForCdToDeKPi", false, "Require deuteron PID for Cd->deKpi"};
    // ML models for triggers
    Configurable<bool> applyMlForHfFilters{"applyMlForHfFilters", false, "Flag to enable ML application for HF Filters"};
    Configurable<std::string> mlModelPathCCDB{"mlModelPathCCDB", "EventFiltering/PWGHF/BDTSmeared", "Path on CCDB of ML models for HF Filters"};
    Configurable<int64_t> timestampCcdbForHfFilters{"timestampCcdbForHfFilters", 1657032422771, "timestamp of the ONNX file for ML model used to query in CCDB"};
    Configurable<bool> loadMlModelsFromCCDB{"loadMlModelsFromCCDB", true, "Flag to enable or disable the loading of ML models from CCDB"};

    Configurable<LabeledArray<std::string>> onnxFileNames{"onnxFileNames", {hf_cuts_bdt_multiclass::onnxFileNameSpecies[0], 5, 1, hf_cuts_bdt_multiclass::labelsSpecies, hf_cuts_bdt_multiclass::labelsModels}, "ONNX file names for ML models"};

    Configurable<LabeledArray<double>> thresholdMlScoreD0ToKPi{"thresholdMlScoreD0ToKPi", {hf_cuts_bdt_multiclass::Cuts[0], hf_cuts_bdt_multiclass::NBinsPt, hf_cuts_bdt_multiclass::NCutBdtScores, hf_cuts_bdt_multiclass::labelsPt, hf_cuts_bdt_multiclass::labelsCutBdt}, "Threshold values for Ml output scores of D0 candidates"};
    Configurable<LabeledArray<double>> thresholdMlScoreDplusToPiKPi{"thresholdMlScoreDplusToPiKPi", {hf_cuts_bdt_multiclass::Cuts[0], hf_cuts_bdt_multiclass::NBinsPt, hf_cuts_bdt_multiclass::NCutBdtScores, hf_cuts_bdt_multiclass::labelsPt, hf_cuts_bdt_multiclass::labelsCutBdt}, "Threshold values for Ml output scores of D+ candidates"};
    Configurable<LabeledArray<double>> thresholdMlScoreDsToPiKK{"thresholdMlScoreDsToPiKK", {hf_cuts_bdt_multiclass::Cuts[0], hf_cuts_bdt_multiclass::NBinsPt, hf_cuts_bdt_multiclass::NCutBdtScores, hf_cuts_bdt_multiclass::labelsPt, hf_cuts_bdt_multiclass::labelsCutBdt}, "Threshold values for Ml output scores of Ds+ candidates"};
    Configurable<LabeledArray<double>> thresholdMlScoreLcToPiKP{"thresholdMlScoreLcToPiKP", {hf_cuts_bdt_multiclass::Cuts[0], hf_cuts_bdt_multiclass::NBinsPt, hf_cuts_bdt_multiclass::NCutBdtScores, hf_cuts_bdt_multiclass::labelsPt, hf_cuts_bdt_multiclass::labelsCutBdt}, "Threshold values for Ml output scores of Lc+ candidates"};
    Configurable<LabeledArray<double>> thresholdMlScoreXicToPiKP{"thresholdMlScoreXicToPiKP", {hf_cuts_bdt_multiclass::Cuts[0], hf_cuts_bdt_multiclass::NBinsPt, hf_cuts_bdt_multiclass::NCutBdtScores, hf_cuts_bdt_multiclass::labelsPt, hf_cuts_bdt_multiclass::labelsCutBdt}, "Threshold values for Ml output scores of Xic+ candidates"};
  } config;

  SliceCache cache;
  o2::vertexing::DCAFitterN<2> df2; // 2-prong vertex fitter
  o2::vertexing::DCAFitterN<3> df3; // 3-prong vertex fitter
  // Needed for PV refitting
  Service<o2::ccdb::BasicCCDBManager> ccdb{};
  o2::base::MatLayerCylSet* lut{};
  o2::base::Propagator::MatCorrType noMatCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;
  int runNumber{};

  // int nColls{0}; //can be added to run over limited collisions per file - for tesing purposes

  static constexpr int kN2ProngDecays = hf_cand_2prong::DecayType::N2ProngDecays;                                                                                                                                                    // number of 2-prong hadron types
  static constexpr int kN3ProngDecays = hf_cand_3prong::DecayType::N3ProngDecays;                                                                                                                                                    // number of 3-prong hadron types
  static constexpr int kNCuts2Prong[kN2ProngDecays] = {hf_cuts_presel_2prong::NCutVars, hf_cuts_presel_2prong::NCutVars, hf_cuts_presel_2prong::NCutVars};                                                                           // how many different selections are made on 2-prongs
  static constexpr int kNCuts3Prong[kN3ProngDecays] = {hf_cuts_presel_3prong::NCutVars, hf_cuts_presel_3prong::NCutVars + 1, hf_cuts_presel_ds::NCutVars, hf_cuts_presel_3prong::NCutVars + 1, hf_cuts_presel_3prong::NCutVars + 1}; // how many different selections are made on 3-prongs (Lc and Xic have also PID potentially)
  static constexpr int kNCutsDstar = 3;                                                                                                                                                                                              // how many different selections are made on Dstars
  std::array<std::array<std::array<double, 2>, 2>, kN2ProngDecays> arrMass2Prong{};
  std::array<std::array<std::array<double, 3>, 2>, kN3ProngDecays> arrMass3Prong{};
  // arrays of 2-prong and 3-prong cuts
  std::array<LabeledArray<double>, kN2ProngDecays> cut2Prong{};
  std::array<std::vector<double>, kN2ProngDecays> binsPt2Prong{};
  std::array<LabeledArray<double>, kN3ProngDecays> cut3Prong{};
  std::array<std::vector<double>, kN3ProngDecays> binsPt3Prong{};

  // ML response
  o2::analysis::MlResponse<float> hfMlResponse2Prongs;                               // only D0
  std::array<o2::analysis::MlResponse<float>, kN3ProngDecays> hfMlResponse3Prongs{}; // D+, Lc, Ds, Xic
  std::array<bool, kN3ProngDecays> hasMlModel3Prong{false};
  o2::ccdb::CcdbApi ccdbApi;

  using SelectedCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::HfSelCollision>>;
  using TracksWithPVRefitAndDCA = soa::Join<aod::TracksWCovDcaExtra, aod::HfPvRefitTrack>;
  using FilteredTrackAssocSel = soa::Filtered<soa::Join<aod::TrackAssoc, aod::HfSelTrack>>;

  // filter collisions
  Filter filterSelectCollisions = (aod::hf_sel_collision::whyRejectColl == static_cast<o2::hf_evsel::HfCollisionRejectionMask>(0));
  // filter track indices
  Filter filterSelectTrackIds = ((aod::hf_sel_track::isSelProng & static_cast<uint32_t>(BIT(CandidateType::Cand2Prong))) != 0u) || ((aod::hf_sel_track::isSelProng & static_cast<uint32_t>(BIT(CandidateType::Cand3Prong))) != 0u) || ((aod::hf_sel_track::isSelProng & static_cast<uint32_t>(BIT(CandidateType::CandDstar))) != 0u);

  Preslice<FilteredTrackAssocSel> trackIndicesPerCollision = aod::track_association::collisionId;
  // define slice of track indices per collisions
  Preslice<TracksWithPVRefitAndDCA> tracksPerCollision = aod::track::collisionId; // needed for PV refit

  // define partitions
  Partition<FilteredTrackAssocSel> positiveFor2And3Prongs = aod::hf_sel_track::isPositive == true && (((aod::hf_sel_track::isSelProng & static_cast<uint32_t>(BIT(CandidateType::Cand2Prong))) != 0u) || ((aod::hf_sel_track::isSelProng & static_cast<uint32_t>(BIT(CandidateType::Cand3Prong))) != 0u));
  Partition<FilteredTrackAssocSel> negativeFor2And3Prongs = aod::hf_sel_track::isPositive == false && (((aod::hf_sel_track::isSelProng & static_cast<uint32_t>(BIT(CandidateType::Cand2Prong))) != 0u) || ((aod::hf_sel_track::isSelProng & static_cast<uint32_t>(BIT(CandidateType::Cand3Prong))) != 0u));
  Partition<FilteredTrackAssocSel> positiveSoftPions = aod::hf_sel_track::isPositive == true && ((aod::hf_sel_track::isSelProng & static_cast<uint32_t>(BIT(CandidateType::CandDstar))) != 0u);
  Partition<FilteredTrackAssocSel> negativeSoftPions = aod::hf_sel_track::isPositive == false && ((aod::hf_sel_track::isSelProng & static_cast<uint32_t>(BIT(CandidateType::CandDstar))) != 0u);

  // QA of PV refit
  ConfigurableAxis axisPvRefitDeltaX{"axisPvRefitDeltaX", {1000, -0.5f, 0.5f}, "DeltaX binning PV refit"};
  ConfigurableAxis axisPvRefitDeltaY{"axisPvRefitDeltaY", {1000, -0.5f, 0.5f}, "DeltaY binning PV refit"};
  ConfigurableAxis axisPvRefitDeltaZ{"axisPvRefitDeltaZ", {1000, -0.5f, 0.5f}, "DeltaZ binning PV refit"};

  HistogramRegistry registry{"registry"};

  void init(InitContext const&)
  {
    if (!doprocess2And3ProngsWithPvRefit && !doprocess2And3ProngsNoPvRefit && !doprocess2And3ProngsWithPvRefitWithPidForHfFiltersBdt && !doprocess2And3ProngsNoPvRefitWithPidForHfFiltersBdt) {
      return;
    }

    arrMass2Prong[hf_cand_2prong::DecayType::D0ToPiK] = std::array{std::array{MassPiPlus, MassKPlus},
                                                                   std::array{MassKPlus, MassPiPlus}};

    arrMass2Prong[hf_cand_2prong::DecayType::JpsiToEE] = std::array{std::array{MassElectron, MassElectron},
                                                                    std::array{MassElectron, MassElectron}};

    arrMass2Prong[hf_cand_2prong::DecayType::JpsiToMuMu] = std::array{std::array{MassMuonPlus, MassMuonPlus},
                                                                      std::array{MassMuonPlus, MassMuonPlus}};

    arrMass3Prong[hf_cand_3prong::DecayType::DplusToPiKPi] = std::array{std::array{MassPiPlus, MassKPlus, MassPiPlus},
                                                                        std::array{MassPiPlus, MassKPlus, MassPiPlus}};

    arrMass3Prong[hf_cand_3prong::DecayType::LcToPKPi] = std::array{std::array{MassProton, MassKPlus, MassPiPlus},
                                                                    std::array{MassPiPlus, MassKPlus, MassProton}};

    arrMass3Prong[hf_cand_3prong::DecayType::DsToKKPi] = std::array{std::array{MassKPlus, MassKPlus, MassPiPlus},
                                                                    std::array{MassPiPlus, MassKPlus, MassKPlus}};

    arrMass3Prong[hf_cand_3prong::DecayType::XicToPKPi] = std::array{std::array{MassProton, MassKPlus, MassPiPlus},
                                                                     std::array{MassPiPlus, MassKPlus, MassProton}};

    arrMass3Prong[hf_cand_3prong::DecayType::CdToDeKPi] = std::array{std::array{MassDeuteron, MassKPlus, MassPiPlus},
                                                                     std::array{MassPiPlus, MassKPlus, MassDeuteron}};

    // cuts for 2-prong decays retrieved by json. the order must be then one in hf_cand_2prong::DecayType
    cut2Prong = {config.cutsD0ToPiK, config.cutsJpsiToEE, config.cutsJpsiToMuMu};
    binsPt2Prong = {config.binsPtD0ToPiK, config.binsPtJpsiToEE, config.binsPtJpsiToMuMu};
    // cuts for 3-prong decays retrieved by json. the order must be then one in hf_cand_3prong::DecayType
    cut3Prong = {config.cutsDplusToPiKPi, config.cutsLcToPKPi, config.cutsDsToKKPi, config.cutsXicToPKPi, config.cutsCdToDeKPi};
    binsPt3Prong = {config.binsPtDplusToPiKPi, config.binsPtLcToPKPi, config.binsPtDsToKKPi, config.binsPtXicToPKPi, config.binsPtCdToDeKPi};

    df2.setPropagateToPCA(config.propagateToPCA);
    df2.setMaxR(config.maxR);
    df2.setMaxDZIni(config.maxDZIni);
    df2.setMinParamChange(config.minParamChange);
    df2.setMinRelChi2Change(config.minRelChi2Change);
    df2.setUseAbsDCA(config.useAbsDCA);
    df2.setWeightedFinalPCA(config.useWeightedFinalPCA);

    df3.setPropagateToPCA(config.propagateToPCA);
    df3.setMaxR(config.maxR);
    df3.setMaxDZIni(config.maxDZIni);
    df3.setMinParamChange(config.minParamChange);
    df3.setMinRelChi2Change(config.minRelChi2Change);
    df3.setUseAbsDCA(config.useAbsDCA);
    df3.setWeightedFinalPCA(config.useWeightedFinalPCA);

    ccdb->setURL(config.ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>(config.ccdbPathLut));
    runNumber = 0;

    if (config.fillHistograms) {
      const AxisSpec axisNumTracks{500, -0.5f, 499.5f, "Number of tracks"};
      const AxisSpec axisNumCands{1000, -0.5f, 999.5f, "Number of candidates"};
      registry.add("hNTracks", "Number of selected tracks;# of selected tracks;entries", {HistType::kTH1D, {axisNumTracks}});
      // 2-prong histograms
      registry.add("hVtx2ProngX", "2-prong candidates;#it{x}_{sec. vtx.} (cm);entries", {HistType::kTH1D, {{1000, -2., 2.}}});
      registry.add("hVtx2ProngY", "2-prong candidates;#it{y}_{sec. vtx.} (cm);entries", {HistType::kTH1D, {{1000, -2., 2.}}});
      registry.add("hVtx2ProngZ", "2-prong candidates;#it{z}_{sec. vtx.} (cm);entries", {HistType::kTH1D, {{1000, -20., 20.}}});
      registry.add("hNCand2Prong", "2-prong candidates preselected;# of candidates;entries", {HistType::kTH1D, {axisNumCands}});
      registry.add("hNCand2ProngVsNTracks", "2-prong candidates preselected;# of selected tracks;# of candidates;entries", {HistType::kTH2D, {axisNumTracks, axisNumCands}});
      registry.add("hMassD0ToPiK", "D^{0} candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH1D, {{500, 0., 5.}}});
      registry.add("hMassJpsiToEE", "J/#psi candidates;inv. mass (e^{#plus} e^{#minus}) (GeV/#it{c}^{2});entries", {HistType::kTH1D, {{500, 0., 5.}}});
      registry.add("hMassJpsiToMuMu", "J/#psi candidates;inv. mass (#mu^{#plus} #mu^{#minus}) (GeV/#it{c}^{2});entries", {HistType::kTH1D, {{500, 0., 5.}}});
      // 3-prong histograms
      registry.add("hVtx3ProngX", "3-prong candidates;#it{x}_{sec. vtx.} (cm);entries", {HistType::kTH1D, {{1000, -2., 2.}}});
      registry.add("hVtx3ProngY", "3-prong candidates;#it{y}_{sec. vtx.} (cm);entries", {HistType::kTH1D, {{1000, -2., 2.}}});
      registry.add("hVtx3ProngZ", "3-prong candidates;#it{z}_{sec. vtx.} (cm);entries", {HistType::kTH1D, {{1000, -20., 20.}}});
      registry.add("hNCand3Prong", "3-prong candidates preselected;# of candidates;entries", {HistType::kTH1D, {axisNumCands}});
      registry.add("hNCand3ProngVsNTracks", "3-prong candidates preselected;# of selected tracks;# of candidates;entries", {HistType::kTH2D, {axisNumTracks, axisNumCands}});
      registry.add("hMassDPlusToPiKPi", "D^{#plus} candidates;inv. mass (#pi K #pi) (GeV/#it{c}^{2});entries", {HistType::kTH1D, {{500, 0., 5.}}});
      registry.add("hMassLcToPKPi", "#Lambda_{c}^{#plus} candidates;inv. mass (p K #pi) (GeV/#it{c}^{2});entries", {HistType::kTH1D, {{500, 0., 5.}}});
      registry.add("hMassDsToKKPi", "D_{s}^{#plus} candidates;inv. mass (K K #pi) (GeV/#it{c}^{2});entries", {HistType::kTH1D, {{500, 0., 5.}}});
      registry.add("hMassXicToPKPi", "#Xi_{c}^{#plus} candidates;inv. mass (p K #pi) (GeV/#it{c}^{2});entries", {HistType::kTH1D, {{500, 0., 5.}}});
      registry.add("hMassDstarToD0Pi", "D^{*#plus} candidates;inv. mass (K #pi #pi) - mass (K #pi) (GeV/#it{c}^{2});entries", {HistType::kTH1D, {{500, 0.135, 0.185}}});
      registry.add("hMassCdToDeKPi", "C Deuteron candidates;inv. mass (De K #pi) (GeV/#it{c}^{2});entries", {HistType::kTH1D, {{500, 0., 5.}}});

      // needed for PV refitting
      if (doprocess2And3ProngsWithPvRefit || doprocess2And3ProngsWithPvRefitWithPidForHfFiltersBdt) {
        const AxisSpec axisCollisionX{100, -20.f, 20.f, "X (cm)"};
        const AxisSpec axisCollisionY{100, -20.f, 20.f, "Y (cm)"};
        const AxisSpec axisCollisionZ{100, -20.f, 20.f, "Z (cm)"};
        const AxisSpec axisCollisionXOriginal{1000, -20.f, 20.f, "X original PV (cm)"};
        const AxisSpec axisCollisionYOriginal{1000, -20.f, 20.f, "Y original PV (cm)"};
        const AxisSpec axisCollisionZOriginal{1000, -20.f, 20.f, "Z original PV (cm)"};
        const AxisSpec axisCollisionNContrib{1000, 0, 1000, "Number of contributors"};
        const AxisSpec axisCollisionDeltaX{axisPvRefitDeltaX, "#Delta x_{PV} (cm)"};
        const AxisSpec axisCollisionDeltaY{axisPvRefitDeltaY, "#Delta y_{PV} (cm)"};
        const AxisSpec axisCollisionDeltaZ{axisPvRefitDeltaZ, "#Delta z_{PV} (cm)"};
        registry.add("PvRefit/verticesPerCandidate", "", kTH1D, {{6, 0.5f, 6.5f, ""}});
        registry.get<TH1>(HIST("PvRefit/verticesPerCandidate"))->GetXaxis()->SetBinLabel(1, "All PV");
        registry.get<TH1>(HIST("PvRefit/verticesPerCandidate"))->GetXaxis()->SetBinLabel(2, "PV refit doable");
        registry.get<TH1>(HIST("PvRefit/verticesPerCandidate"))->GetXaxis()->SetBinLabel(3, "PV refit #chi^{2}!=-1");
        registry.get<TH1>(HIST("PvRefit/verticesPerCandidate"))->GetXaxis()->SetBinLabel(4, "PV refit #it{#chi}^{2}==#minus1");
        registry.get<TH1>(HIST("PvRefit/verticesPerCandidate"))->GetXaxis()->SetBinLabel(5, "1 daughter contr.");
        registry.get<TH1>(HIST("PvRefit/verticesPerCandidate"))->GetXaxis()->SetBinLabel(6, "no PV refit");
        registry.add("PvRefit/hPvDeltaXvsNContrib", "", kTH2D, {axisCollisionNContrib, axisCollisionDeltaX});
        registry.add("PvRefit/hPvDeltaYvsNContrib", "", kTH2D, {axisCollisionNContrib, axisCollisionDeltaY});
        registry.add("PvRefit/hPvDeltaZvsNContrib", "", kTH2D, {axisCollisionNContrib, axisCollisionDeltaZ});
        registry.add("PvRefit/hChi2vsNContrib", "", kTH2D, {axisCollisionNContrib, {102, -1.5, 100.5, "#chi^{2} PV refit"}});
        registry.add("PvRefit/hPvRefitXChi2Minus1", "PV refit with #it{#chi}^{2}==#minus1", kTH2D, {axisCollisionX, axisCollisionXOriginal});
        registry.add("PvRefit/hPvRefitYChi2Minus1", "PV refit with #it{#chi}^{2}==#minus1", kTH2D, {axisCollisionY, axisCollisionYOriginal});
        registry.add("PvRefit/hPvRefitZChi2Minus1", "PV refit with #it{#chi}^{2}==#minus1", kTH2D, {axisCollisionZ, axisCollisionZOriginal});
        registry.add("PvRefit/hNContribPvRefitNotDoable", "N. contributors for PV refit not doable", kTH1D, {axisCollisionNContrib});
        registry.add("PvRefit/hNContribPvRefitChi2Minus1", "N. contributors original PV for PV refit #it{#chi}^{2}==#minus1", kTH1D, {axisCollisionNContrib});
      }

      if (config.applyMlForHfFilters) {
        const AxisSpec axisBdtScore{100, 0.f, 1.f};
        registry.add("ML/hMlScoreBkgD0", "Bkg ML score for D^{0} candidates;Bkg ML score;entries", kTH1D, {axisBdtScore});
        registry.add("ML/hMlScorePromptD0", "Prompt ML score for D^{0} candidates;Prompt ML score;entries", kTH1D, {axisBdtScore});
        registry.add("ML/hMlScoreNonpromptD0", "Non-prompt ML score for D^{0} candidates;Non-prompt ML score;entries", kTH1D, {axisBdtScore});
        registry.add("ML/hMlScoreBkgDplus", "Bkg ML score for D^{#plus} candidates;Bkg ML score;entries", kTH1D, {axisBdtScore});
        registry.add("ML/hMlScorePromptDplus", "Prompt ML score for D^{#plus} candidates;Prompt ML score;entries", kTH1D, {axisBdtScore});
        registry.add("ML/hMlScoreNonpromptDplus", "Non-prompt ML score for D^{#plus} candidates;Non-prompt ML score;entries", kTH1D, {axisBdtScore});
        registry.add("ML/hMlScoreBkgDs", "Bkg ML score for D_{s}^{#plus} candidates;Bkg ML score;entries", kTH1D, {axisBdtScore});
        registry.add("ML/hMlScorePromptDs", "Prompt ML score for D_{s}^{#plus} candidates;Prompt ML score;entries", kTH1D, {axisBdtScore});
        registry.add("ML/hMlScoreNonpromptDs", "Non-prompt ML score for D_{s}^{#plus} candidates;Non-prompt ML score;entries", kTH1D, {axisBdtScore});
        registry.add("ML/hMlScoreBkgLc", "Bkg ML score for #Lambda_{c}^{#plus} candidates;Bkg ML score;entries", kTH1D, {axisBdtScore});
        registry.add("ML/hMlScorePromptLc", "Prompt ML score for #Lambda_{c}^{#plus} candidates;Prompt ML score;entries", kTH1D, {axisBdtScore});
        registry.add("ML/hMlScoreNonpromptLc", "Non-prompt ML score for #Lambda_{c}^{#plus} candidates;Non-prompt ML score;entries", kTH1D, {axisBdtScore});
        registry.add("ML/hMlScoreBkgXic", "Bkg ML score for #Xi_{c}^{#plus} candidates;Bkg ML score;entries", kTH1D, {axisBdtScore});
        registry.add("ML/hMlScorePromptXic", "Prompt ML score for #Xi_{c}^{#plus} candidates;Prompt ML score;entries", kTH1D, {axisBdtScore});
        registry.add("ML/hMlScoreNonpromptXic", "Non-prompt ML score for #Xi_{c}^{#plus} candidates;Non-prompt ML score;entries", kTH1D, {axisBdtScore});
      }
    }

    if (config.applyMlForHfFilters) {
      const std::vector<std::string> onnxFileNames2Prongs{config.onnxFileNames->get(0u, 0u)};
      // Exclude Cd from the 3-prong list, as it is not included in the pp trigger program
      const std::array<std::vector<std::string>, kN3ProngDecays - 1> onnxFileNames3Prongs{std::vector<std::string>{config.onnxFileNames->get(1u, 0u)}, std::vector<std::string>{config.onnxFileNames->get(2u, 0u)}, std::vector<std::string>{config.onnxFileNames->get(3u, 0u)}, std::vector<std::string>{config.onnxFileNames->get(4u, 0u)}};
      const std::vector<std::string> mlModelPathCcdb2Prongs{config.mlModelPathCCDB.value + "D0"};
      const std::array<std::vector<std::string>, kN3ProngDecays - 1> mlModelPathCcdb3Prongs{std::vector<std::string>{config.mlModelPathCCDB.value + "Dplus"}, std::vector<std::string>{config.mlModelPathCCDB.value + "Lc"}, std::vector<std::string>{config.mlModelPathCCDB.value + "Ds"}, std::vector<std::string>{config.mlModelPathCCDB.value + "Xic"}};
      const std::vector<double> ptBinsMl{0., 1.e10};
      const std::vector<int> cutDirMl{o2::cuts_ml::CutDirection::CutGreater, o2::cuts_ml::CutDirection::CutSmaller, o2::cuts_ml::CutDirection::CutSmaller};
      const std::array<LabeledArray<double>, kN3ProngDecays - 1> thresholdMlScore3Prongs{config.thresholdMlScoreDplusToPiKPi, config.thresholdMlScoreLcToPiKP, config.thresholdMlScoreDsToPiKK, config.thresholdMlScoreXicToPiKP};

      // initialise 2-prong ML response
      hfMlResponse2Prongs.configure(ptBinsMl, config.thresholdMlScoreD0ToKPi, cutDirMl, 3);
      if (config.loadMlModelsFromCCDB) {
        ccdbApi.init(config.ccdbUrl);
        hfMlResponse2Prongs.setModelPathsCCDB(onnxFileNames2Prongs, ccdbApi, mlModelPathCcdb2Prongs, config.timestampCcdbForHfFilters);
      } else {
        hfMlResponse2Prongs.setModelPathsLocal(onnxFileNames2Prongs);
      }
      hfMlResponse2Prongs.init();

      // initialise 3-prong ML responses
      for (int iDecay3P{0}; iDecay3P < kN3ProngDecays - 1; ++iDecay3P) {
        if (onnxFileNames3Prongs[iDecay3P][0].empty()) { // 3-prong species to be skipped
          continue;
        }
        hasMlModel3Prong[iDecay3P] = true;
        hfMlResponse3Prongs[iDecay3P].configure(ptBinsMl, thresholdMlScore3Prongs[iDecay3P], cutDirMl, 3);
        if (config.loadMlModelsFromCCDB) {
          ccdbApi.init(config.ccdbUrl);
          hfMlResponse3Prongs[iDecay3P].setModelPathsCCDB(onnxFileNames3Prongs[iDecay3P], ccdbApi, mlModelPathCcdb3Prongs[iDecay3P], config.timestampCcdbForHfFilters);
        } else {
          hfMlResponse3Prongs[iDecay3P].setModelPathsLocal(onnxFileNames3Prongs[iDecay3P]);
        }
        hfMlResponse3Prongs[iDecay3P].init();
      }
    }
  }

  /// Method to perform selections for 2-prong candidates before vertex reconstruction
  /// \param pVecTrack0 is the momentum array of the first daughter track
  /// \param pVecTrack1 is the momentum array of the second daughter track
  /// \param dcaTrack0 is the dcaXY of the first daughter track
  /// \param dcaTrack1 is the dcaXY of the second daughter track
  /// \param cutStatus is a 2D array with outcome of each selection (filled only in debug mode)
  /// \param whichHypo information of the mass hypoteses that were selected
  /// \param isSelected is a bitmap with selection outcome
  /// \param pt2Prong is the pt of the 2-prong candidate
  template <typename T1, typename T2, typename T3, typename T4>
  void applyPreselection2Prong(T1 const& pVecTrack0, T1 const& pVecTrack1, T2 const& dcaTrack0, T2 const& dcaTrack1, T3& cutStatus, T4& whichHypo, auto& isSelected, float& pt2Prong)
  {
    whichHypo[kN2ProngDecays] = 0; // D0 for D*

    pt2Prong = RecoDecay::pt(pVecTrack0, pVecTrack1);
    const auto pt = pt2Prong + config.ptTolerance; // add tolerance because of no reco decay vertex

    for (int iDecay2P = 0; iDecay2P < kN2ProngDecays; iDecay2P++) {

      // pT
      const auto binPt = findBin(&binsPt2Prong[iDecay2P], pt);
      // return immediately if it is outside the defined pT bins
      if (binPt == -1) {
        CLRBIT(isSelected, iDecay2P);
        if (config.debug) {
          cutStatus[iDecay2P][0] = false;
        }
        continue;
      }

      // invariant mass
      double massHypos[2] = {0., 0.};
      whichHypo[iDecay2P] = 3; // 2 bits on

      if (config.debug || TESTBIT(isSelected, iDecay2P)) {
        const double minMass = cut2Prong[iDecay2P].get(binPt, 0u);
        const double maxMass = cut2Prong[iDecay2P].get(binPt, 1u);
        if (minMass >= 0. && maxMass > 0.) {
          const std::array arrMom{pVecTrack0, pVecTrack1};
          massHypos[0] = RecoDecay::m2(arrMom, arrMass2Prong[iDecay2P][0]);
          massHypos[1] = (iDecay2P == hf_cand_2prong::DecayType::D0ToPiK) ? RecoDecay::m2(arrMom, arrMass2Prong[iDecay2P][1]) : massHypos[0];
          const double min2 = minMass * minMass;
          const double max2 = maxMass * maxMass;
          if (massHypos[0] < min2 || massHypos[0] >= max2) {
            CLRBIT(whichHypo[iDecay2P], 0);
          }
          if (massHypos[1] < min2 || massHypos[1] >= max2) {
            CLRBIT(whichHypo[iDecay2P], 1);
          }
          if (whichHypo[iDecay2P] == 0) {
            CLRBIT(isSelected, iDecay2P);
            if (config.debug) {
              cutStatus[iDecay2P][1] = false;
            }
          }
        }
      }

      // imp. par. product cut
      if (config.debug || TESTBIT(isSelected, iDecay2P)) {
        const auto impParProduct = dcaTrack0 * dcaTrack1;
        if (impParProduct > cut2Prong[iDecay2P].get(binPt, 3u)) {
          CLRBIT(isSelected, iDecay2P);
          if (config.debug) {
            cutStatus[iDecay2P][2] = false;
          }
        }
      }

      // additional check for D0 to be used in D* finding
      if (iDecay2P == hf_cand_2prong::DecayType::D0ToPiK && config.doDstar && TESTBIT(isSelected, iDecay2P)) {
        const auto binPtDstar = findBin(config.binsPtDstarToD0Pi, pt * 1.2); // assuming the D* pT about 20% higher than the one of the D0 to be safe
        if (binPtDstar >= 0) {
          whichHypo[kN2ProngDecays] = whichHypo[hf_cand_2prong::DecayType::D0ToPiK];
          const double deltaMass = config.cutsDstarToD0Pi->get(binPtDstar, 1u);

          if (TESTBIT(whichHypo[iDecay2P], 0) && (massHypos[0] > (MassD0 + deltaMass) * (MassD0 + deltaMass) || massHypos[0] < (MassD0 - deltaMass) * (MassD0 - deltaMass))) {
            CLRBIT(whichHypo[kN2ProngDecays], 0);
          }
          if (TESTBIT(whichHypo[iDecay2P], 1) && (massHypos[1] > (MassD0 + deltaMass) * (MassD0 + deltaMass) || massHypos[1] < (MassD0 - deltaMass) * (MassD0 - deltaMass))) {
            CLRBIT(whichHypo[kN2ProngDecays], 1);
          }
        }
      }
    }
  }

  /// Method to perform selections on difference from nominal mass for phi decay
  /// \param binPt pt bin for the cuts
  /// \param pVecTrack0 is the momentum array of the first daughter track
  /// \param pVecTrack1 is the momentum array of the second daughter track
  /// \param pVecTrack2 is the momentum array of the third daughter track
  /// \param cutStatus is a 2D array with outcome of each selection (filled only in debug mode)
  /// \param whichHypo information of the mass hypoteses that were selected
  /// \param isSelected is a bitmap with selection outcome
  template <typename T1, typename T2, typename T3>
  void applyPreselectionPhiDecay(const int binPt, T1 const& pVecTrack0, T1 const& pVecTrack1, T1 const& pVecTrack2, T2& cutStatus, T3& whichHypo, auto& isSelected)
  {
    const double deltaMassMax = cut3Prong[hf_cand_3prong::DecayType::DsToKKPi].get(binPt, 4u);
    if (TESTBIT(whichHypo[hf_cand_3prong::DecayType::DsToKKPi], 0)) {
      const double mass2PhiKKPi = RecoDecay::m2(std::array{pVecTrack0, pVecTrack1}, std::array{arrMass3Prong[hf_cand_3prong::DecayType::DsToKKPi][0][0], arrMass3Prong[hf_cand_3prong::DecayType::DsToKKPi][0][1]});
      if (mass2PhiKKPi > (MassPhi + deltaMassMax) * (MassPhi + deltaMassMax) || mass2PhiKKPi < (MassPhi - deltaMassMax) * (MassPhi - deltaMassMax)) {
        CLRBIT(whichHypo[hf_cand_3prong::DecayType::DsToKKPi], 0);
      }
    }
    if (TESTBIT(whichHypo[hf_cand_3prong::DecayType::DsToKKPi], 1)) {
      const double mass2PhiPiKK = RecoDecay::m2(std::array{pVecTrack1, pVecTrack2}, std::array{arrMass3Prong[hf_cand_3prong::DecayType::DsToKKPi][1][1], arrMass3Prong[hf_cand_3prong::DecayType::DsToKKPi][1][2]});
      if (mass2PhiPiKK > (MassPhi + deltaMassMax) * (MassPhi + deltaMassMax) || mass2PhiPiKK < (MassPhi - deltaMassMax) * (MassPhi - deltaMassMax)) {
        CLRBIT(whichHypo[hf_cand_3prong::DecayType::DsToKKPi], 1);
      }
    }
    if (whichHypo[hf_cand_3prong::DecayType::DsToKKPi] == 0) {
      CLRBIT(isSelected, hf_cand_3prong::DecayType::DsToKKPi);
      if (config.debug) {
        cutStatus[hf_cand_3prong::DecayType::DsToKKPi][4] = false;
      }
    }
  }

  /// Method to perform selections for 3-prong candidates before vertex reconstruction
  /// \param pVecTrack0 is the momentum array of the first daughter track
  /// \param pVecTrack1 is the momentum array of the second daughter track
  /// \param pVecTrack2 is the momentum array of the third daughter track
  /// \param isIdentifiedPidTrack0 is the flag that tells if the track 0 has been tagged as a proton
  /// \param isIdentifiedPidTrack2 is the flag that tells if the track 2 has been tagged as a proton
  /// \param cutStatus is a 2D array with outcome of each selection (filled only in debug mode)
  /// \param whichHypo information of the mass hypoteses that were selected
  /// \param isSelected is a bitmap with selection outcome
  template <typename T1, typename T2, typename T3>
  void applyPreselection3Prong(T1 const& pVecTrack0, T1 const& pVecTrack1, T1 const& pVecTrack2, const auto isIdentifiedPidTrack0, const auto isIdentifiedPidTrack2, T2& cutStatus, T3& whichHypo, auto& isSelected)
  {
    const auto pt = RecoDecay::pt(pVecTrack0, pVecTrack1, pVecTrack2) + config.ptTolerance; // add tolerance because of no reco decay vertex

    for (int iDecay3P = 0; iDecay3P < kN3ProngDecays; iDecay3P++) {

      // check proton PID for Lc and Xic
      whichHypo[iDecay3P] = 3; // 2 bits on

      if ((iDecay3P == hf_cand_3prong::DecayType::LcToPKPi && config.applyProtonPidForLcToPKPi) || (iDecay3P == hf_cand_3prong::DecayType::XicToPKPi && config.applyProtonPidForXicToPKPi) || (iDecay3P == hf_cand_3prong::DecayType::CdToDeKPi && config.applyDeuteronPidForCdToDeKPi)) {

        if ((iDecay3P == hf_cand_3prong::DecayType::LcToPKPi && !TESTBIT(isIdentifiedPidTrack0, ChannelsProtonPid::LcToPKPi)) || (iDecay3P == hf_cand_3prong::DecayType::XicToPKPi && !TESTBIT(isIdentifiedPidTrack0, ChannelsProtonPid::XicToPKPi)) || (iDecay3P == hf_cand_3prong::DecayType::CdToDeKPi && !TESTBIT(isIdentifiedPidTrack0, ChannelsDeuteronPid))) {
          CLRBIT(whichHypo[iDecay3P], 0);
        }
        if ((iDecay3P == hf_cand_3prong::DecayType::LcToPKPi && !TESTBIT(isIdentifiedPidTrack2, ChannelsProtonPid::LcToPKPi)) || (iDecay3P == hf_cand_3prong::DecayType::XicToPKPi && !TESTBIT(isIdentifiedPidTrack2, ChannelsProtonPid::XicToPKPi)) || (iDecay3P == hf_cand_3prong::DecayType::CdToDeKPi && !TESTBIT(isIdentifiedPidTrack2, ChannelsDeuteronPid))) {
          CLRBIT(whichHypo[iDecay3P], 1);
        }
        if (whichHypo[iDecay3P] == 0) {
          CLRBIT(isSelected, iDecay3P);
          if (config.debug) {
            cutStatus[iDecay3P][hf_cuts_presel_3prong::NCutVars] = false; // PID
          }
          continue; // no need to check further for this particle hypothesis
        }
      }
      // pT
      const auto binPt = findBin(&binsPt3Prong[iDecay3P], pt);
      // return immediately if it is outside the defined pT bins
      if (binPt == -1) {
        CLRBIT(isSelected, iDecay3P);
        whichHypo[iDecay3P] = 0;
        if (config.debug) {
          cutStatus[iDecay3P][0] = false;
        }
        continue;
      }

      // invariant mass
      if ((config.debug || TESTBIT(isSelected, iDecay3P))) {
        const double minMass = cut3Prong[iDecay3P].get(binPt, 0u);
        const double maxMass = cut3Prong[iDecay3P].get(binPt, 1u);
        if (minMass >= 0. && maxMass > 0.) { // no need to check isSelected but to avoid mistakes
          double massHypos[2] = {0., 0.};
          const std::array arrMom{pVecTrack0, pVecTrack1, pVecTrack2};
          const double min2 = minMass * minMass;
          const double max2 = maxMass * maxMass;
          massHypos[0] = RecoDecay::m2(arrMom, arrMass3Prong[iDecay3P][0]);
          massHypos[1] = (iDecay3P != hf_cand_3prong::DecayType::DplusToPiKPi) ? RecoDecay::m2(arrMom, arrMass3Prong[iDecay3P][1]) : massHypos[0];
          if (massHypos[0] < min2 || massHypos[0] >= max2) {
            CLRBIT(whichHypo[iDecay3P], 0);
          }
          if (massHypos[1] < min2 || massHypos[1] >= max2) {
            CLRBIT(whichHypo[iDecay3P], 1);
          }
          if (whichHypo[iDecay3P] == 0) {
            CLRBIT(isSelected, iDecay3P);
            if (config.debug) {
              cutStatus[iDecay3P][1] = false;
            }
          }
        }
      }

      if ((config.debug || TESTBIT(isSelected, iDecay3P)) && iDecay3P == hf_cand_3prong::DecayType::DsToKKPi) {
        applyPreselectionPhiDecay(binPt, pVecTrack0, pVecTrack1, pVecTrack2, cutStatus, whichHypo, isSelected);
      }
    }
  }

  /// Method to perform selections for 2-prong candidates after vertex reconstruction
  /// \param pVecCand is the array for the candidate momentum after reconstruction of secondary vertex
  /// \param secVtx is the secondary vertex
  /// \param primVtx is the primary vertex
  /// \param cutStatus is a 2D array with outcome of each selection (filled only in debug mode)
  /// \param isSelected ia s bitmap with selection outcome
  template <typename T1, typename T2, typename T3, typename T4>
  void applySelection2Prong(const T1& pVecCand, const T2& secVtx, const T3& primVtx, T4& cutStatus, auto& isSelected)
  {
    if (config.debug || isSelected > 0) {

      for (int iDecay2P = 0; iDecay2P < kN2ProngDecays; iDecay2P++) {

        // pT
        const auto binPt = findBin(&binsPt2Prong[iDecay2P], RecoDecay::pt(pVecCand));
        if (binPt == -1) { // cut if it is outside the defined pT bins
          CLRBIT(isSelected, iDecay2P);
          if (config.debug) {
            cutStatus[iDecay2P][0] = false;
          }
          continue;
        }

        // cos of pointing angle
        if (config.debug || TESTBIT(isSelected, iDecay2P)) {
          const auto cpa = RecoDecay::cpa(primVtx, secVtx, pVecCand);
          if (cpa < cut2Prong[iDecay2P].get(binPt, 2u)) { // 2u == "cospIndex[iDecay2P]"
            CLRBIT(isSelected, iDecay2P);
            if (config.debug) {
              cutStatus[iDecay2P][3] = false;
            }
          }
        }
      }
    }
  }

  /// Method to perform ML selections for 2-prong candidates after the rectangular selections
  /// \param featuresCand is the vector with the candidate features
  /// \param outputScores is the vector with the output scores to be filled
  /// \param isSelected ia s bitmap with selection outcome
  void applyMlSelectionForHfFilters2Prong(std::vector<float> featuresCand, std::vector<float>& outputScores, auto& isSelected)
  {
    if (!TESTBIT(isSelected, hf_cand_2prong::DecayType::D0ToPiK)) {
      return;
    }
    const float ptDummy = 1.; // dummy pT value (only one pT bin)
    const bool isSelMl = hfMlResponse2Prongs.isSelectedMl(featuresCand, ptDummy, outputScores);
    if (config.fillHistograms) {
      registry.fill(HIST("ML/hMlScoreBkgD0"), outputScores[0]);
      registry.fill(HIST("ML/hMlScorePromptD0"), outputScores[1]);
      registry.fill(HIST("ML/hMlScoreNonpromptD0"), outputScores[2]);
    }
    if (!isSelMl) {
      CLRBIT(isSelected, hf_cand_2prong::DecayType::D0ToPiK);
    }
  }

  /// Method to perform selections for 2-prong candidates after vertex reconstruction
  /// \param secVtx is the secondary vertex
  /// \param primVtx is the primary vertex
  /// \param dcaFitter is the DCAFitter used for the 2-track vertex
  /// \returns true if the candidate is selected
  template <typename T1, typename T2, typename T3>
  bool isTwoTrackVertexSelectedFor3Prongs(const T1& secVtx, const T2& primVtx, const T3& dcaFitter)
  {
    if (dcaFitter.getChi2AtPCACandidate() > config.maxTwoTrackChi2PcaFor3Prongs) {
      return false;
    }
    const auto decLen = RecoDecay::distance(primVtx, secVtx);
    return static_cast<bool>(decLen >= config.minTwoTrackDecayLengthFor3Prongs);
  }

  /// Method to perform selections for 3-prong candidates after vertex reconstruction
  /// \param pVecCand is the array for the candidate momentum after reconstruction of secondary vertex
  /// \param secVtx is the secondary vertex
  /// \param primVtx is the primary vertex
  /// \param cutStatus is a 2D array with outcome of each selection (filled only in debug mode)
  /// \param isSelected ia s bitmap with selection outcome
  template <typename T1, typename T2, typename T3, typename T4>
  void applySelection3Prong(const T1& pVecCand, const T2& secVtx, const T3& primVtx, T4& cutStatus, auto& isSelected)
  {
    if (config.debug || isSelected > 0) {

      const auto pt = RecoDecay::pt(pVecCand);
      for (int iDecay3P = 0; iDecay3P < kN3ProngDecays; iDecay3P++) {

        // pT
        const auto binPt = findBin(&binsPt3Prong[iDecay3P], pt);
        if (binPt == -1) { // cut if it is outside the defined pT bins
          CLRBIT(isSelected, iDecay3P);
          if (config.debug) {
            cutStatus[iDecay3P][0] = false;
          }
          continue;
        }

        // cos of pointing angle
        if (config.debug || TESTBIT(isSelected, iDecay3P)) {
          const auto cpa = RecoDecay::cpa(primVtx, secVtx, pVecCand);
          if (cpa < cut3Prong[iDecay3P].get(binPt, 2u)) { // 2u == cospIndex[iDecay3P]
            CLRBIT(isSelected, iDecay3P);
            if (config.debug) {
              cutStatus[iDecay3P][2] = false;
            }
          }
        }

        // decay length
        if ((config.debug || TESTBIT(isSelected, iDecay3P))) {
          const auto decayLength = RecoDecay::distance(primVtx, secVtx);
          if (decayLength < cut3Prong[iDecay3P].get(binPt, 3u)) { // 3u == decLenIndex[iDecay3P]
            CLRBIT(isSelected, iDecay3P);
            if (config.debug) {
              cutStatus[iDecay3P][3] = false;
            }
          }
        }
      }
    }
  }

  /// Method to perform ML selections for 2-prong candidates after the rectangular selections
  /// \tparam usePidForHfFiltersBdt is the flag to determine whether to use also the PID features for the Lc BDT
  /// \param featuresCand is the vector with the candidate features
  /// \param featuresCandPid is the vector with the candidate PID features
  /// \param outputScores is the array of vectors with the output scores to be filled
  /// \param isSelected ia s bitmap with selection outcome
  template <bool UsePidForHfFiltersBdt>
  void applyMlSelectionForHfFilters3Prong(std::vector<float> featuresCand, std::vector<float> featuresCandPid, std::array<std::vector<float>, kN3ProngDecays - 1>& outputScores, auto& isSelected)
  {
    if (isSelected == 0) {
      return;
    }

    const float ptDummy = 1.f; // dummy pT value (only one pT bin)
    for (int iDecay3P{0}; iDecay3P < kN3ProngDecays - 1; ++iDecay3P) {
      if (TESTBIT(isSelected, iDecay3P) && hasMlModel3Prong[iDecay3P]) {
        bool isMlSel = false;
        if constexpr (UsePidForHfFiltersBdt) {
          if (iDecay3P != hf_cand_3prong::DecayType::LcToPKPi && iDecay3P != hf_cand_3prong::DecayType::XicToPKPi) {
            isMlSel = hfMlResponse3Prongs[iDecay3P].isSelectedMl(featuresCand, ptDummy, outputScores[iDecay3P]);
          } else {
            std::vector<float> featuresCandWithPid{featuresCand};
            featuresCandWithPid.insert(featuresCandWithPid.end(), featuresCandPid.begin(), featuresCandPid.end());
            isMlSel = hfMlResponse3Prongs[iDecay3P].isSelectedMl(featuresCandWithPid, ptDummy, outputScores[iDecay3P]);
          }
        } else {
          isMlSel = hfMlResponse3Prongs[iDecay3P].isSelectedMl(featuresCand, ptDummy, outputScores[iDecay3P]);
        }
        if (config.fillHistograms) {
          switch (iDecay3P) {
            case hf_cand_3prong::DecayType::DplusToPiKPi: {
              registry.fill(HIST("ML/hMlScoreBkgDplus"), outputScores[iDecay3P][0]);
              registry.fill(HIST("ML/hMlScorePromptDplus"), outputScores[iDecay3P][1]);
              registry.fill(HIST("ML/hMlScoreNonpromptDplus"), outputScores[iDecay3P][2]);
              break;
            }
            case hf_cand_3prong::DecayType::LcToPKPi: {
              registry.fill(HIST("ML/hMlScoreBkgLc"), outputScores[iDecay3P][0]);
              registry.fill(HIST("ML/hMlScorePromptLc"), outputScores[iDecay3P][1]);
              registry.fill(HIST("ML/hMlScoreNonpromptLc"), outputScores[iDecay3P][2]);
              break;
            }
            case hf_cand_3prong::DecayType::DsToKKPi: {
              registry.fill(HIST("ML/hMlScoreBkgDs"), outputScores[iDecay3P][0]);
              registry.fill(HIST("ML/hMlScorePromptDs"), outputScores[iDecay3P][1]);
              registry.fill(HIST("ML/hMlScoreNonpromptDs"), outputScores[iDecay3P][2]);
              break;
            }
            case hf_cand_3prong::DecayType::XicToPKPi: {
              registry.fill(HIST("ML/hMlScoreBkgXic"), outputScores[iDecay3P][0]);
              registry.fill(HIST("ML/hMlScorePromptXic"), outputScores[iDecay3P][1]);
              registry.fill(HIST("ML/hMlScoreNonpromptXic"), outputScores[iDecay3P][2]);
              break;
            }
          }
        }
        if (!isMlSel) {
          CLRBIT(isSelected, iDecay3P);
        }
      }
    }
  }

  /// Method to perform selections for D* candidates before vertex reconstruction
  /// \param pVecTrack0 is the momentum array of the first daughter track (same charge)
  /// \param pVecTrack1 is the momentum array of the second daughter track (opposite charge)
  /// \param pVecTrack2 is the momentum array of the third daughter track (same charge)
  /// \param cutStatus is the cut status (filled only in case of debug)
  /// \param deltaMass is the M(Kpipi) - M(Kpi) value to be filled in the control histogram
  /// \return a bitmap with selection outcome
  template <typename T1, typename T2>
  uint8_t applySelectionDstar(T1 const& pVecTrack0, T1 const& pVecTrack1, T1 const& pVecTrack2, uint8_t& cutStatus, T2& deltaMass)
  {
    uint8_t isSelected{1};
    const std::array arrMom{pVecTrack0, pVecTrack1, pVecTrack2};
    const std::array arrMomD0{pVecTrack0, pVecTrack1};
    const auto pt = RecoDecay::pt(pVecTrack0, pVecTrack1, pVecTrack2) + config.ptTolerance; // add tolerance because of no reco decay vertex

    // pT
    const auto binPt = findBin(config.binsPtDstarToD0Pi, pt);
    // return immediately if it is outside the defined pT bins
    if (binPt == -1) {
      isSelected = 0;
      if (config.debug) {
        CLRBIT(cutStatus, 0);
      }
      return isSelected;
    }

    // D0 mass
    const double deltaMassD0 = config.cutsDstarToD0Pi->get(binPt, 1u); // 1u == deltaMassD0Index
    const double invMassD0 = RecoDecay::m(arrMomD0, std::array{MassPiPlus, MassKPlus});
    if (std::abs(invMassD0 - MassD0) > deltaMassD0) {
      isSelected = 0;
      if (config.debug) {
        CLRBIT(cutStatus, 1);
      }
      return isSelected;
    }

    // D*+ mass
    const double maxDeltaMass = config.cutsDstarToD0Pi->get(binPt, 0u); // 0u == deltaMassIndex
    const double invMassDstar = RecoDecay::m(arrMom, std::array{MassPiPlus, MassKPlus, MassPiPlus});
    deltaMass = invMassDstar - invMassD0;
    if (deltaMass > maxDeltaMass) {
      isSelected = 0;
      if (config.debug) {
        CLRBIT(cutStatus, 1);
      }
      return isSelected;
    }

    return isSelected;
  }

  /// Method for the PV refit excluding the candidate daughters
  /// \param collision is a collision
  /// \param bcWithTimeStamps is a table of bunch crossing joined with timestamps used to query the CCDB for B and material budget
  /// \param vecPvContributorGlobId is a vector containing the global ID of PV contributors for the current collision
  /// \param vecPvContributorTrackParCov is a vector containing the TrackParCov of PV contributors for the current collision
  /// \param vecCandPvContributorGlobId is a vector containing the global indices of daughter tracks that contributed to the original PV refit
  /// \param pvCoord is a vector where to store X, Y and Z values of refitted PV
  /// \param pvCovMatrix is a vector where to store the covariance matrix values of refitted PV
  void performPvRefitCandProngs(SelectedCollisions::iterator const& collision,
                                aod::BCsWithTimestamps const&,
                                std::vector<int64_t> const& vecPvContributorGlobId,
                                std::vector<o2::track::TrackParCov> const& vecPvContributorTrackParCov,
                                std::vector<int64_t> const& vecCandPvContributorGlobId,
                                std::array<float, 3>& pvCoord,
                                std::array<float, 6>& pvCovMatrix)
  {
    std::vector<bool> vecPvRefitContributorUsed(vecPvContributorGlobId.size(), true);

    /// Prepare the vertex refitting
    // set the magnetic field from CCDB
    const auto bc = collision.bc_as<o2::aod::BCsWithTimestamps>();
    initCCDB(bc, runNumber, ccdb, config.isRun2 ? config.ccdbPathGrp : config.ccdbPathGrpMag, lut, config.isRun2);

    // build the VertexBase to initialize the vertexer
    o2::dataformats::VertexBase primVtx;
    primVtx.setX(collision.posX());
    primVtx.setY(collision.posY());
    primVtx.setZ(collision.posZ());
    primVtx.setCov(collision.covXX(), collision.covXY(), collision.covYY(), collision.covXZ(), collision.covYZ(), collision.covZZ());
    // configure PVertexer
    o2::vertexing::PVertexer vertexer;
    o2::conf::ConfigurableParam::updateFromString("pvertexer.useMeanVertexConstraint=false"); /// remove diamond constraint (let's keep it at the moment...)
    vertexer.init();
    const bool pvRefitDoable = vertexer.prepareVertexRefit(vecPvContributorTrackParCov, primVtx);
    if (!pvRefitDoable) {
      LOG(info) << "Not enough tracks accepted for the refit";
      if ((doprocess2And3ProngsWithPvRefit || doprocess2And3ProngsWithPvRefitWithPidForHfFiltersBdt) && config.fillHistograms) {
        registry.fill(HIST("PvRefit/hNContribPvRefitNotDoable"), collision.numContrib());
      }
    }
    if (config.debugPvRefit) {
      LOG(info) << "prepareVertexRefit = " << pvRefitDoable << " Ncontrib= " << vecPvContributorTrackParCov.size() << " Ntracks= " << collision.numContrib() << " Vtx= " << primVtx.asString();
    }

    /// PV refitting, if the tracks contributed to this at the beginning
    o2::dataformats::VertexBase primVtxBaseRecalc;
    if ((doprocess2And3ProngsWithPvRefit || doprocess2And3ProngsWithPvRefitWithPidForHfFiltersBdt) && pvRefitDoable) {
      if (config.fillHistograms) {
        registry.fill(HIST("PvRefit/verticesPerCandidate"), 2);
      }
      bool recalcPvRefit = true;
      int nCandContr = 0;
      for (const uint64_t myGlobalID : vecCandPvContributorGlobId) {                                              // o2-linter: disable=const-ref-in-for-loop (small type)
        auto trackIterator = std::find(vecPvContributorGlobId.begin(), vecPvContributorGlobId.end(), myGlobalID); /// track global index
        if (trackIterator != vecPvContributorGlobId.end()) {
          /// this is a contributor, let's remove it for the PV refit
          const int entry = std::distance(vecPvContributorGlobId.begin(), trackIterator);
          vecPvRefitContributorUsed[entry] = false; /// remove the track from the PV refitting
          nCandContr++;
        }
      }

      /// do the PV refit excluding the candidate daughters that originally contributed to fit it
      if (config.debugPvRefit) {
        LOG(info) << "### PV refit after removing " << nCandContr << " tracks";
      }
      auto primVtxRefitted = vertexer.refitVertex(vecPvRefitContributorUsed, primVtx); // vertex refit
      // LOG(info) << "refit " << cnt << "/" << ntr << " result = " << primVtxRefitted.asString();
      // LOG(info) << "refit for track with global index " << static_cast<int>(myTrack.globalIndex()) << " " << primVtxRefitted.asString();
      if (primVtxRefitted.getChi2() < 0) {
        if (config.debugPvRefit) {
          LOG(info) << "---> Refitted vertex has bad chi2 = " << primVtxRefitted.getChi2();
        }
        if (config.fillHistograms) {
          registry.fill(HIST("PvRefit/verticesPerCandidate"), 4);
          registry.fill(HIST("PvRefit/hPvRefitXChi2Minus1"), primVtxRefitted.getX(), collision.posX());
          registry.fill(HIST("PvRefit/hPvRefitYChi2Minus1"), primVtxRefitted.getY(), collision.posY());
          registry.fill(HIST("PvRefit/hPvRefitZChi2Minus1"), primVtxRefitted.getZ(), collision.posZ());
          registry.fill(HIST("PvRefit/hNContribPvRefitChi2Minus1"), collision.numContrib());
        }
        recalcPvRefit = false;
      } else if (config.fillHistograms) {
        registry.fill(HIST("PvRefit/verticesPerCandidate"), 3);
      }
      if (config.fillHistograms) {
        registry.fill(HIST("PvRefit/hChi2vsNContrib"), primVtxRefitted.getNContributors(), primVtxRefitted.getChi2());
      }

      if (recalcPvRefit) {
        // fill the histograms for refitted PV with good Chi2
        const double deltaX = primVtx.getX() - primVtxRefitted.getX();
        const double deltaY = primVtx.getY() - primVtxRefitted.getY();
        const double deltaZ = primVtx.getZ() - primVtxRefitted.getZ();
        if (config.fillHistograms) {
          registry.fill(HIST("PvRefit/hPvDeltaXvsNContrib"), primVtxRefitted.getNContributors(), deltaX);
          registry.fill(HIST("PvRefit/hPvDeltaYvsNContrib"), primVtxRefitted.getNContributors(), deltaY);
          registry.fill(HIST("PvRefit/hPvDeltaZvsNContrib"), primVtxRefitted.getNContributors(), deltaZ);
        }

        // fill the newly calculated PV
        primVtxBaseRecalc.setX(primVtxRefitted.getX());
        primVtxBaseRecalc.setY(primVtxRefitted.getY());
        primVtxBaseRecalc.setZ(primVtxRefitted.getZ());
        primVtxBaseRecalc.setCov(primVtxRefitted.getSigmaX2(), primVtxRefitted.getSigmaXY(), primVtxRefitted.getSigmaY2(), primVtxRefitted.getSigmaXZ(), primVtxRefitted.getSigmaYZ(), primVtxRefitted.getSigmaZ2());

      } else {
        /// copy the original collision PV
        primVtxBaseRecalc.setX(primVtx.getX());
        primVtxBaseRecalc.setY(primVtx.getY());
        primVtxBaseRecalc.setZ(primVtx.getZ());
        primVtxBaseRecalc.setCov(primVtx.getSigmaX2(), primVtx.getSigmaXY(), primVtx.getSigmaY2(), primVtx.getSigmaXZ(), primVtx.getSigmaYZ(), primVtx.getSigmaZ2());
      }

      // fill the output
      pvCoord[0] = primVtxBaseRecalc.getX();
      pvCoord[1] = primVtxBaseRecalc.getY();
      pvCoord[2] = primVtxBaseRecalc.getZ();
      pvCovMatrix[0] = primVtxBaseRecalc.getSigmaX2();
      pvCovMatrix[1] = primVtxBaseRecalc.getSigmaXY();
      pvCovMatrix[2] = primVtxBaseRecalc.getSigmaY2();
      pvCovMatrix[3] = primVtxBaseRecalc.getSigmaXZ();
      pvCovMatrix[4] = primVtxBaseRecalc.getSigmaYZ();
      pvCovMatrix[5] = primVtxBaseRecalc.getSigmaZ2();

      // cnt++;

    } /// end 'if (doprocess2And3ProngsWithPvRefit && pvRefitDoable)'

  } /// end of performPvRefitCandProngs function

  template <bool DoPvRefit, bool UsePidForHfFiltersBdt, typename TTracks>
  void run2And3Prongs(SelectedCollisions const& collisions,
                      aod::BCsWithTimestamps const& bcWithTimeStamps,
                      FilteredTrackAssocSel const&,
                      TTracks const& tracks)
  {

    // can be added to run over limited collisions per file - for tesing purposes
    /*
    if (nCollsMax > -1){
      if (nColls == nCollMax){
        return;
        //can be added to run over limited collisions per file - for tesing purposes
      }
      nColls++;
    }
    */

    for (const auto& collision : collisions) {

      /// retrieve PV contributors for the current collision
      std::vector<int64_t> vecPvContributorGlobId{};
      std::vector<o2::track::TrackParCov> vecPvContributorTrackParCov{};
      std::vector<bool> vecPvRefitContributorUsed{};
      if constexpr (DoPvRefit) {
        auto groupedTracksUnfiltered = tracks.sliceBy(tracksPerCollision, collision.globalIndex());
        const int nTrk = groupedTracksUnfiltered.size();
        int nContrib = 0;
        int nNonContrib = 0;
        for (const auto& trackUnfiltered : groupedTracksUnfiltered) {
          if (!trackUnfiltered.isPVContributor()) {
            /// the track did not contribute to fit the primary vertex
            nNonContrib++;
            continue;
          }
          vecPvContributorGlobId.push_back(trackUnfiltered.globalIndex());
          vecPvContributorTrackParCov.push_back(getTrackParCov(trackUnfiltered));
          nContrib++;
          if (config.debugPvRefit) {
            LOG(info) << "---> a contributor! stuff saved";
            LOG(info) << "vec_contrib size: " << vecPvContributorTrackParCov.size() << ", nContrib: " << nContrib;
          }
        }
        if (config.debugPvRefit) {
          LOG(info) << "===> nTrk: " << nTrk << ",   nContrib: " << nContrib << ",   nNonContrib: " << nNonContrib;
          if (static_cast<uint16_t>(vecPvContributorTrackParCov.size()) != collision.numContrib() || static_cast<uint16_t>(nContrib != collision.numContrib())) {
            LOG(info) << "!!! Some problem here !!! vecPvContributorTrackParCov.size()= " << vecPvContributorTrackParCov.size() << ", nContrib=" << nContrib << ", collision.numContrib()" << collision.numContrib();
          }
        }
        vecPvRefitContributorUsed = std::vector<bool>(vecPvContributorGlobId.size(), true);
      }

      // auto centrality = collision.centV0M(); //FIXME add centrality when option for variations to the process function appears

      const auto n2ProngBit = BIT(kN2ProngDecays) - 1; // bit value for 2-prong candidates where each candidate is one bit and they are all set to 1
      const auto n3ProngBit = BIT(kN3ProngDecays) - 1; // bit value for 3-prong candidates where each candidate is one bit and they are all set to 1

      std::array<std::vector<bool>, kN2ProngDecays> cutStatus2Prong{};
      std::array<std::vector<bool>, kN3ProngDecays> cutStatus3Prong{};
      uint8_t nCutStatus2ProngBit[kN2ProngDecays]; // bit value for selection status for each 2-prong candidate where each selection is one bit and they are all set to 1
      uint8_t nCutStatus3ProngBit[kN3ProngDecays]; // bit value for selection status for each 3-prong candidate where each selection is one bit and they are all set to 1

      for (int iDecay2P = 0; iDecay2P < kN2ProngDecays; iDecay2P++) {
        nCutStatus2ProngBit[iDecay2P] = BIT(kNCuts2Prong[iDecay2P]) - 1;
        cutStatus2Prong[iDecay2P] = std::vector<bool>(kNCuts2Prong[iDecay2P], true);
      }
      for (int iDecay3P = 0; iDecay3P < kN3ProngDecays; iDecay3P++) {
        nCutStatus3ProngBit[iDecay3P] = BIT(kNCuts3Prong[iDecay3P]) - 1;
        cutStatus3Prong[iDecay3P] = std::vector<bool>(kNCuts3Prong[iDecay3P], true);
      }

      int whichHypo2Prong[kN2ProngDecays + 1]; // we also put D0 for D* in the last slot
      int whichHypo3Prong[kN3ProngDecays];

      // set the magnetic field from CCDB
      const auto bc = collision.bc_as<o2::aod::BCsWithTimestamps>();
      initCCDB(bc, runNumber, ccdb, config.isRun2 ? config.ccdbPathGrp : config.ccdbPathGrpMag, lut, config.isRun2);
      df2.setBz(o2::base::Propagator::Instance()->getNominalBz());
      df3.setBz(o2::base::Propagator::Instance()->getNominalBz());

      // used to calculate number of candidiates per event
      auto nCand2 = rowTrackIndexProng2.lastIndex();
      auto nCand3 = rowTrackIndexProng3.lastIndex();

      // if there isn't at least a positive and a negative track, continue immediately
      // if (tracksPos.size() < 1 || tracksNeg.size() < 1) {
      //  return;
      //}

      const auto thisCollId = collision.globalIndex();

      // first loop over positive tracks
      const auto groupedTrackIndicesPos1 = positiveFor2And3Prongs->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
      int lastFilledD0 = -1; // index to be filled in table for D* mesons
      for (auto trackIndexPos1 = groupedTrackIndicesPos1.begin(); trackIndexPos1 != groupedTrackIndicesPos1.end(); ++trackIndexPos1) {
        const auto trackPos1 = trackIndexPos1.template track_as<TTracks>();

        // retrieve the selection flag that corresponds to this collision
        const auto isSelProngPos1 = trackIndexPos1.isSelProng();
        const bool sel2ProngStatusPos = TESTBIT(isSelProngPos1, CandidateType::Cand2Prong);
        const bool sel3ProngStatusPos1 = TESTBIT(isSelProngPos1, CandidateType::Cand3Prong);

        auto trackParVarPos1 = getTrackParCov(trackPos1);
        std::array pVecTrackPos1{trackPos1.pVector()};
        std::array dcaInfoPos1{trackPos1.dcaXY(), trackPos1.dcaZ()};
        if (thisCollId != trackPos1.collisionId()) { // this is not the "default" collision for this track, we have to re-propagate it
          o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParVarPos1, 2.f, noMatCorr, &dcaInfoPos1);
          getPxPyPz(trackParVarPos1, pVecTrackPos1);
        }

        // first loop over negative tracks
        const auto groupedTrackIndicesNeg1 = negativeFor2And3Prongs->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
        for (auto trackIndexNeg1 = groupedTrackIndicesNeg1.begin(); trackIndexNeg1 != groupedTrackIndicesNeg1.end(); ++trackIndexNeg1) {
          const auto trackNeg1 = trackIndexNeg1.template track_as<TTracks>();

          // retrieve the selection flag that corresponds to this collision
          const auto isSelProngNeg1 = trackIndexNeg1.isSelProng();
          const bool sel2ProngStatusNeg = TESTBIT(isSelProngNeg1, CandidateType::Cand2Prong);
          const bool sel3ProngStatusNeg1 = TESTBIT(isSelProngNeg1, CandidateType::Cand3Prong);

          auto trackParVarNeg1 = getTrackParCov(trackNeg1);
          std::array pVecTrackNeg1{trackNeg1.pVector()};
          std::array dcaInfoNeg1{trackNeg1.dcaXY(), trackNeg1.dcaZ()};
          if (thisCollId != trackNeg1.collisionId()) { // this is not the "default" collision for this track, we have to re-propagate it
            o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParVarNeg1, 2.f, noMatCorr, &dcaInfoNeg1);
            getPxPyPz(trackParVarNeg1, pVecTrackNeg1);
          }

          uint isSelected2ProngCand = n2ProngBit; // bitmap for checking status of two-prong candidates (1 is true, 0 is rejected)

          if (config.debug) {
            for (int iDecay2P = 0; iDecay2P < kN2ProngDecays; iDecay2P++) {
              for (int iCut = 0; iCut < kNCuts2Prong[iDecay2P]; iCut++) {
                cutStatus2Prong[iDecay2P][iCut] = true;
              }
            }
          }

          // initialise PV refit coordinates and cov matrix for 2-prongs already here for D*
          std::array pvRefitCoord2Prong = {collision.posX(), collision.posY(), collision.posZ()}; /// initialize to the original PV
          std::array pvRefitCovMatrix2Prong = getPrimaryVertex(collision).getCov();               /// initialize to the original PV

          // 2-prong vertex reconstruction
          float pt2Prong{-1.};
          bool is2ProngCandidateGoodFor3Prong{sel3ProngStatusPos1 && sel3ProngStatusNeg1};
          int nVtxFrom2ProngFitter = 0;
          if (sel2ProngStatusPos && sel2ProngStatusNeg) {

            // 2-prong preselections
            // TODO: in case of PV refit, the single-track DCA is calculated wrt two different PV vertices (only 1 track excluded)
            applyPreselection2Prong(pVecTrackPos1, pVecTrackNeg1, dcaInfoPos1[0], dcaInfoNeg1[0], cutStatus2Prong, whichHypo2Prong, isSelected2ProngCand, pt2Prong);

            if (isSelected2ProngCand > 0) {
              // secondary vertex reconstruction and further 2-prong selections
              try {
                nVtxFrom2ProngFitter = df2.process(trackParVarPos1, trackParVarNeg1);
              } catch (...) {
              }

              if (nVtxFrom2ProngFitter > 0) { // should it be this or > 0 or are they equivalent
                // get secondary vertex
                const auto& secondaryVertex2 = df2.getPCACandidate();
                // get track momenta
                std::array<float, 3> pvec0{};
                std::array<float, 3> pvec1{};
                df2.getTrack(0).getPxPyPzGlo(pvec0);
                df2.getTrack(1).getPxPyPzGlo(pvec1);

                /// PV refit excluding the candidate daughters, if contributors
                if constexpr (DoPvRefit) {
                  if (config.fillHistograms) {
                    registry.fill(HIST("PvRefit/verticesPerCandidate"), 1);
                  }
                  int nCandContr = 2;
                  auto trackFirstIt = std::find(vecPvContributorGlobId.begin(), vecPvContributorGlobId.end(), trackPos1.globalIndex());
                  auto trackSecondIt = std::find(vecPvContributorGlobId.begin(), vecPvContributorGlobId.end(), trackNeg1.globalIndex());
                  bool isTrackFirstContr = true;
                  bool isTrackSecondContr = true;
                  if (trackFirstIt == vecPvContributorGlobId.end()) {
                    /// This track did not contribute to the original PV refit
                    if (config.debugPvRefit) {
                      LOG(info) << "--- [2 Prong] trackPos1 with globalIndex " << trackPos1.globalIndex() << " was not a PV contributor";
                    }
                    nCandContr--;
                    isTrackFirstContr = false;
                  }
                  if (trackSecondIt == vecPvContributorGlobId.end()) {
                    /// This track did not contribute to the original PV refit
                    if (config.debugPvRefit) {
                      LOG(info) << "--- [2 Prong] trackNeg1 with globalIndex " << trackNeg1.globalIndex() << " was not a PV contributor";
                    }
                    nCandContr--;
                    isTrackSecondContr = false;
                  }
                  if (nCandContr == 2) { // o2-linter: disable="magic-number" (see comment below)
                    /// Both the daughter tracks were used for the original PV refit, let's refit it after excluding them
                    if (config.debugPvRefit) {
                      LOG(info) << "### [2 Prong] Calling performPvRefitCandProngs for HF 2 prong candidate";
                    }
                    performPvRefitCandProngs(collision, bcWithTimeStamps, vecPvContributorGlobId, vecPvContributorTrackParCov, {trackPos1.globalIndex(), trackNeg1.globalIndex()}, pvRefitCoord2Prong, pvRefitCovMatrix2Prong);
                  } else if (nCandContr == 1) {
                    /// Only one daughter was a contributor, let's use then the PV recalculated by excluding only it
                    if (config.debugPvRefit) {
                      LOG(info) << "####### [2 Prong] nCandContr==" << nCandContr << " ---> just 1 contributor!";
                    }
                    if (config.fillHistograms) {
                      registry.fill(HIST("PvRefit/verticesPerCandidate"), 5);
                    }
                    if (isTrackFirstContr && !isTrackSecondContr) {
                      /// the first daughter is contributor, the second is not
                      pvRefitCoord2Prong = {trackPos1.pvRefitX(), trackPos1.pvRefitY(), trackPos1.pvRefitZ()};
                      pvRefitCovMatrix2Prong = {trackPos1.pvRefitSigmaX2(), trackPos1.pvRefitSigmaXY(), trackPos1.pvRefitSigmaY2(), trackPos1.pvRefitSigmaXZ(), trackPos1.pvRefitSigmaYZ(), trackPos1.pvRefitSigmaZ2()};
                    } else if (!isTrackFirstContr && isTrackSecondContr) {
                      ///  the second daughter is contributor, the first is not
                      pvRefitCoord2Prong = {trackNeg1.pvRefitX(), trackNeg1.pvRefitY(), trackNeg1.pvRefitZ()};
                      pvRefitCovMatrix2Prong = {trackNeg1.pvRefitSigmaX2(), trackNeg1.pvRefitSigmaXY(), trackNeg1.pvRefitSigmaY2(), trackNeg1.pvRefitSigmaXZ(), trackNeg1.pvRefitSigmaYZ(), trackNeg1.pvRefitSigmaZ2()};
                    }
                  } else {
                    /// 0 contributors among the HF candidate daughters
                    if (config.fillHistograms) {
                      registry.fill(HIST("PvRefit/verticesPerCandidate"), 6);
                    }
                    if (config.debugPvRefit) {
                      LOG(info) << "####### [2 Prong] nCandContr==" << nCandContr << " ---> some of the candidate daughters did not contribute to the original PV fit, PV refit not redone";
                    }
                  }
                }

                const auto pVecCandProng2 = RecoDecay::pVec(pvec0, pvec1);
                // 2-prong selections after secondary vertex
                std::array pvCoord2Prong = {collision.posX(), collision.posY(), collision.posZ()};
                if constexpr (DoPvRefit) {
                  pvCoord2Prong[0] = pvRefitCoord2Prong[0];
                  pvCoord2Prong[1] = pvRefitCoord2Prong[1];
                  pvCoord2Prong[2] = pvRefitCoord2Prong[2];
                }
                applySelection2Prong(pVecCandProng2, secondaryVertex2, pvCoord2Prong, cutStatus2Prong, isSelected2ProngCand);
                if (is2ProngCandidateGoodFor3Prong && config.do3Prong) {
                  is2ProngCandidateGoodFor3Prong = isTwoTrackVertexSelectedFor3Prongs(secondaryVertex2, pvCoord2Prong, df2);
                }

                std::vector<float> mlScoresD0{};
                if (config.applyMlForHfFilters) {
                  const auto trackParVarPcaPos1 = df2.getTrack(0);
                  const auto trackParVarPcaNeg1 = df2.getTrack(1);
                  const std::vector<float> inputFeatures{trackParVarPcaPos1.getPt(), dcaInfoPos1[0], dcaInfoPos1[1], trackParVarPcaNeg1.getPt(), dcaInfoNeg1[0], dcaInfoNeg1[1]};
                  applyMlSelectionForHfFilters2Prong(inputFeatures, mlScoresD0, isSelected2ProngCand);
                }

                if (isSelected2ProngCand > 0) {
                  // fill table row
                  rowTrackIndexProng2(thisCollId, trackPos1.globalIndex(), trackNeg1.globalIndex(), isSelected2ProngCand);
                  if (config.applyMlForHfFilters) {
                    rowTrackIndexMlScoreProng2(mlScoresD0);
                  }
                  if (TESTBIT(isSelected2ProngCand, hf_cand_2prong::DecayType::D0ToPiK)) {
                    lastFilledD0 = rowTrackIndexProng2.lastIndex();
                  }

                  if constexpr (DoPvRefit) {
                    // fill table row with coordinates of PV refit
                    rowProng2PVrefit(pvRefitCoord2Prong[0], pvRefitCoord2Prong[1], pvRefitCoord2Prong[2],
                                     pvRefitCovMatrix2Prong[0], pvRefitCovMatrix2Prong[1], pvRefitCovMatrix2Prong[2], pvRefitCovMatrix2Prong[3], pvRefitCovMatrix2Prong[4], pvRefitCovMatrix2Prong[5]);
                  }

                  if (config.debug) {
                    uint8_t prong2CutStatus[kN2ProngDecays];
                    for (int iDecay2P = 0; iDecay2P < kN2ProngDecays; iDecay2P++) {
                      prong2CutStatus[iDecay2P] = nCutStatus2ProngBit[iDecay2P];
                      for (int iCut = 0; iCut < kNCuts2Prong[iDecay2P]; iCut++) {
                        if (!cutStatus2Prong[iDecay2P][iCut]) {
                          CLRBIT(prong2CutStatus[iDecay2P], iCut);
                        }
                      }
                    }
                    rowProng2CutStatus(prong2CutStatus[0], prong2CutStatus[1], prong2CutStatus[2]); // FIXME when we can do this by looping over kN2ProngDecays
                  }

                  // fill histograms
                  if (config.fillHistograms) {
                    registry.fill(HIST("hVtx2ProngX"), secondaryVertex2[0]);
                    registry.fill(HIST("hVtx2ProngY"), secondaryVertex2[1]);
                    registry.fill(HIST("hVtx2ProngZ"), secondaryVertex2[2]);
                    const std::array arrMom{pvec0, pvec1};
                    for (int iDecay2P = 0; iDecay2P < kN2ProngDecays; iDecay2P++) {
                      if (TESTBIT(isSelected2ProngCand, iDecay2P)) {
                        if (TESTBIT(whichHypo2Prong[iDecay2P], 0)) {
                          const auto mass2Prong = RecoDecay::m(arrMom, arrMass2Prong[iDecay2P][0]);
                          switch (iDecay2P) {
                            case hf_cand_2prong::DecayType::D0ToPiK:
                              registry.fill(HIST("hMassD0ToPiK"), mass2Prong);
                              break;
                            case hf_cand_2prong::DecayType::JpsiToEE:
                              registry.fill(HIST("hMassJpsiToEE"), mass2Prong);
                              break;
                            case hf_cand_2prong::DecayType::JpsiToMuMu:
                              registry.fill(HIST("hMassJpsiToMuMu"), mass2Prong);
                              break;
                          }
                        }
                        if (TESTBIT(whichHypo2Prong[iDecay2P], 1)) {
                          const auto mass2Prong = RecoDecay::m(arrMom, arrMass2Prong[iDecay2P][1]);
                          if (iDecay2P == hf_cand_2prong::DecayType::D0ToPiK) {
                            registry.fill(HIST("hMassD0ToPiK"), mass2Prong);
                          }
                        }
                      }
                    }
                  }
                }
              } else {
                isSelected2ProngCand = 0; // reset to 0 not to use the D0 to build a D* meson
              }
            } else {
              isSelected2ProngCand = 0; // reset to 0 not to use the D0 to build a D* meson
            }
          }

          // if the cut on the decay length of 3-prongs computed with the first two tracks is enabled and the vertex was not computed for the D0, we compute it now
          if (config.do3Prong && is2ProngCandidateGoodFor3Prong && (config.minTwoTrackDecayLengthFor3Prongs > 0.f || config.maxTwoTrackChi2PcaFor3Prongs < 1.e9f) && nVtxFrom2ProngFitter == 0) { // o2-linter: disable="magic-number" (default maxTwoTrackChi2PcaFor3Prongs is 1.e10)
            try {
              nVtxFrom2ProngFitter = df2.process(trackParVarPos1, trackParVarNeg1);
            } catch (...) {
            }
            if (nVtxFrom2ProngFitter > 0) {
              const auto& secondaryVertex2 = df2.getPCACandidate();
              const std::array pvCoord2Prong{collision.posX(), collision.posY(), collision.posZ()};
              is2ProngCandidateGoodFor3Prong = isTwoTrackVertexSelectedFor3Prongs(secondaryVertex2, pvCoord2Prong, df2);
            } else {
              is2ProngCandidateGoodFor3Prong = false;
            }
          }

          if (config.do3Prong && is2ProngCandidateGoodFor3Prong) { // if 3 prongs are enabled and the first 2 tracks are selected for the 3-prong channels
            // second loop over positive tracks
            for (auto trackIndexPos2 = trackIndexPos1 + 1; trackIndexPos2 != groupedTrackIndicesPos1.end(); ++trackIndexPos2) {

              uint isSelected3ProngCand = n3ProngBit;
              if (!TESTBIT(trackIndexPos2.isSelProng(), CandidateType::Cand3Prong)) { // continue immediately
                if (!config.debug) {
                  continue;
                }
                isSelected3ProngCand = 0;
              }

              if (config.applyKaonPidIn3Prongs && !TESTBIT(trackIndexNeg1.isIdentifiedPid(), ChannelKaonPid)) { // continue immediately if kaon PID enabled and opposite-sign track not a kaon
                if (!config.debug) {
                  continue;
                }
                isSelected3ProngCand = 0;
              }

              const auto trackPos2 = trackIndexPos2.template track_as<TTracks>();
              auto trackParVarPos2 = getTrackParCov(trackPos2);
              std::array dcaInfoPos2{trackPos2.dcaXY(), trackPos2.dcaZ()};

              // preselection of 3-prong candidates
              if (isSelected3ProngCand) {
                std::array pVecTrackPos2{trackPos2.pVector()};
                if (thisCollId != trackPos2.collisionId()) { // this is not the "default" collision for this track and we still did not re-propagate it, we have to re-propagate it
                  o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParVarPos2, 2.f, noMatCorr, &dcaInfoPos2);
                  getPxPyPz(trackParVarPos2, pVecTrackPos2);
                }

                if (config.debug) {
                  for (int iDecay3P = 0; iDecay3P < kN3ProngDecays; iDecay3P++) {
                    for (int iCut = 0; iCut < kNCuts3Prong[iDecay3P]; iCut++) {
                      cutStatus3Prong[iDecay3P][iCut] = true;
                    }
                  }
                }

                // 3-prong preselections
                const auto isIdentifiedPidTrackPos1 = trackIndexPos1.isIdentifiedPid();
                const auto isIdentifiedPidTrackPos2 = trackIndexPos2.isIdentifiedPid();
                applyPreselection3Prong(pVecTrackPos1, pVecTrackNeg1, pVecTrackPos2, isIdentifiedPidTrackPos1, isIdentifiedPidTrackPos2, cutStatus3Prong, whichHypo3Prong, isSelected3ProngCand);
                if (!config.debug && isSelected3ProngCand == 0) {
                  continue;
                }
              }

              /// PV refit excluding the candidate daughters, if contributors
              std::array pvRefitCoord3Prong2Pos1Neg{collision.posX(), collision.posY(), collision.posZ()}; /// initialize to the original PV
              std::array pvRefitCovMatrix3Prong2Pos1Neg{getPrimaryVertex(collision).getCov()};             /// initialize to the original PV
              if constexpr (DoPvRefit) {
                if (config.fillHistograms) {
                  registry.fill(HIST("PvRefit/verticesPerCandidate"), 1);
                }
                int nCandContr = 3;
                auto trackFirstIt = std::find(vecPvContributorGlobId.begin(), vecPvContributorGlobId.end(), trackPos1.globalIndex());
                auto trackSecondIt = std::find(vecPvContributorGlobId.begin(), vecPvContributorGlobId.end(), trackNeg1.globalIndex());
                auto trackThirdIt = std::find(vecPvContributorGlobId.begin(), vecPvContributorGlobId.end(), trackPos2.globalIndex());
                bool isTrackFirstContr = true;
                bool isTrackSecondContr = true;
                bool isTrackThirdContr = true;
                if (trackFirstIt == vecPvContributorGlobId.end()) {
                  /// This track did not contribute to the original PV refit
                  if (config.debugPvRefit) {
                    LOG(info) << "--- [3 prong] trackPos1 with globalIndex " << trackPos1.globalIndex() << " was not a PV contributor";
                  }
                  nCandContr--;
                  isTrackFirstContr = false;
                }
                if (trackSecondIt == vecPvContributorGlobId.end()) {
                  /// This track did not contribute to the original PV refit
                  if (config.debugPvRefit) {
                    LOG(info) << "--- [3 prong] trackNeg1 with globalIndex " << trackNeg1.globalIndex() << " was not a PV contributor";
                  }
                  nCandContr--;
                  isTrackSecondContr = false;
                }
                if (trackThirdIt == vecPvContributorGlobId.end()) {
                  /// This track did not contribute to the original PV refit
                  if (config.debugPvRefit) {
                    LOG(info) << "--- [3 prong] trackPos2 with globalIndex " << trackPos2.globalIndex() << " was not a PV contributor";
                  }
                  nCandContr--;
                  isTrackThirdContr = false;
                }

                // Fill a vector with global ID of candidate daughters that are contributors
                std::vector<int64_t> vecCandPvContributorGlobId = {};
                if (isTrackFirstContr) {
                  vecCandPvContributorGlobId.push_back(trackPos1.globalIndex());
                }
                if (isTrackSecondContr) {
                  vecCandPvContributorGlobId.push_back(trackNeg1.globalIndex());
                }
                if (isTrackThirdContr) {
                  vecCandPvContributorGlobId.push_back(trackPos2.globalIndex());
                }

                if (nCandContr == 3 || nCandContr == 2) { // o2-linter: disable="magic-number" (see comment below)
                  /// At least two of the daughter tracks were used for the original PV refit, let's refit it after excluding them
                  if (config.debugPvRefit) {
                    LOG(info) << "### [3 prong] Calling performPvRefitCandProngs for HF 3 prong candidate, removing " << nCandContr << " daughters";
                  }
                  performPvRefitCandProngs(collision, bcWithTimeStamps, vecPvContributorGlobId, vecPvContributorTrackParCov, vecCandPvContributorGlobId, pvRefitCoord3Prong2Pos1Neg, pvRefitCovMatrix3Prong2Pos1Neg);
                } else if (nCandContr == 1) {
                  /// Only one daughter was a contributor, let's use then the PV recalculated by excluding only it
                  if (config.debugPvRefit) {
                    LOG(info) << "####### [3 Prong] nCandContr==" << nCandContr << " ---> just 1 contributor!";
                  }
                  if (config.fillHistograms) {
                    registry.fill(HIST("PvRefit/verticesPerCandidate"), 5);
                  }
                  if (isTrackFirstContr && !isTrackSecondContr && !isTrackThirdContr) {
                    /// the first daughter is contributor, the second and the third are not
                    pvRefitCoord3Prong2Pos1Neg = {trackPos1.pvRefitX(), trackPos1.pvRefitY(), trackPos1.pvRefitZ()};
                    pvRefitCovMatrix3Prong2Pos1Neg = {trackPos1.pvRefitSigmaX2(), trackPos1.pvRefitSigmaXY(), trackPos1.pvRefitSigmaY2(), trackPos1.pvRefitSigmaXZ(), trackPos1.pvRefitSigmaYZ(), trackPos1.pvRefitSigmaZ2()};
                  } else if (!isTrackFirstContr && isTrackSecondContr && !isTrackThirdContr) {
                    /// the second daughter is contributor, the first and the third are not
                    pvRefitCoord3Prong2Pos1Neg = {trackNeg1.pvRefitX(), trackNeg1.pvRefitY(), trackNeg1.pvRefitZ()};
                    pvRefitCovMatrix3Prong2Pos1Neg = {trackNeg1.pvRefitSigmaX2(), trackNeg1.pvRefitSigmaXY(), trackNeg1.pvRefitSigmaY2(), trackNeg1.pvRefitSigmaXZ(), trackNeg1.pvRefitSigmaYZ(), trackNeg1.pvRefitSigmaZ2()};
                  } else if (!isTrackFirstContr && !isTrackSecondContr && isTrackThirdContr) {
                    /// the third daughter is contributor, the first and the second are not
                    pvRefitCoord3Prong2Pos1Neg = {trackPos2.pvRefitX(), trackPos2.pvRefitY(), trackPos2.pvRefitZ()};
                    pvRefitCovMatrix3Prong2Pos1Neg = {trackPos2.pvRefitSigmaX2(), trackPos2.pvRefitSigmaXY(), trackPos2.pvRefitSigmaY2(), trackPos2.pvRefitSigmaXZ(), trackPos2.pvRefitSigmaYZ(), trackPos2.pvRefitSigmaZ2()};
                  }
                } else {
                  /// 0 contributors among the HF candidate daughters
                  if (config.fillHistograms) {
                    registry.fill(HIST("PvRefit/verticesPerCandidate"), 6);
                  }
                  if (config.debugPvRefit) {
                    LOG(info) << "####### [3 prong] nCandContr==" << nCandContr << " ---> some of the candidate daughters did not contribute to the original PV fit, PV refit not redone";
                  }
                }
              }

              // reconstruct the 3-prong secondary vertex
              int nVtxFrom3ProngFitter = 0;
              try {
                nVtxFrom3ProngFitter = df3.process(trackParVarPos1, trackParVarNeg1, trackParVarPos2);
              } catch (...) {
                continue;
              }

              if (nVtxFrom3ProngFitter == 0) {
                continue;
              }
              // get secondary vertex
              const auto& secondaryVertex3 = df3.getPCACandidate();
              // get track momenta
              std::array<float, 3> pvec0{};
              std::array<float, 3> pvec1{};
              std::array<float, 3> pvec2{};
              const auto trackParVarPcaPos1 = df3.getTrack(0);
              const auto trackParVarPcaNeg1 = df3.getTrack(1);
              const auto trackParVarPcaPos2 = df3.getTrack(2);
              trackParVarPcaPos1.getPxPyPzGlo(pvec0);
              trackParVarPcaNeg1.getPxPyPzGlo(pvec1);
              trackParVarPcaPos2.getPxPyPzGlo(pvec2);
              const auto pVecCandProng3Pos = RecoDecay::pVec(pvec0, pvec1, pvec2);

              // 3-prong selections after secondary vertex
              applySelection3Prong(pVecCandProng3Pos, secondaryVertex3, pvRefitCoord3Prong2Pos1Neg, cutStatus3Prong, isSelected3ProngCand);

              std::array<std::vector<float>, kN3ProngDecays - 1> mlScores3Prongs;
              if (config.applyMlForHfFilters) {
                const std::vector<float> inputFeatures{trackParVarPcaPos1.getPt(), dcaInfoPos1[0], dcaInfoPos1[1], trackParVarPcaNeg1.getPt(), dcaInfoNeg1[0], dcaInfoNeg1[1], trackParVarPcaPos2.getPt(), dcaInfoPos2[0], dcaInfoPos2[1]};
                std::vector<float> inputFeaturesLcPid{};
                if constexpr (UsePidForHfFiltersBdt) {
                  inputFeaturesLcPid.push_back(trackPos1.tpcNSigmaPr());
                  inputFeaturesLcPid.push_back(trackPos2.tpcNSigmaPr());
                  inputFeaturesLcPid.push_back(trackPos1.tpcNSigmaPi());
                  inputFeaturesLcPid.push_back(trackPos2.tpcNSigmaPi());
                  inputFeaturesLcPid.push_back(trackNeg1.tpcNSigmaKa());
                }
                applyMlSelectionForHfFilters3Prong<UsePidForHfFiltersBdt>(inputFeatures, inputFeaturesLcPid, mlScores3Prongs, isSelected3ProngCand);
              }

              if (!config.debug && isSelected3ProngCand == 0) {
                continue;
              }

              // fill table row
              rowTrackIndexProng3(thisCollId, trackPos1.globalIndex(), trackNeg1.globalIndex(), trackPos2.globalIndex(), isSelected3ProngCand);
              if (config.applyMlForHfFilters) {
                rowTrackIndexMlScoreProng3(mlScores3Prongs[0], mlScores3Prongs[1], mlScores3Prongs[2], mlScores3Prongs[3]);
              }
              if constexpr (DoPvRefit) {
                // fill table row of coordinates of PV refit
                rowProng3PVrefit(pvRefitCoord3Prong2Pos1Neg[0], pvRefitCoord3Prong2Pos1Neg[1], pvRefitCoord3Prong2Pos1Neg[2],
                                 pvRefitCovMatrix3Prong2Pos1Neg[0], pvRefitCovMatrix3Prong2Pos1Neg[1], pvRefitCovMatrix3Prong2Pos1Neg[2], pvRefitCovMatrix3Prong2Pos1Neg[3], pvRefitCovMatrix3Prong2Pos1Neg[4], pvRefitCovMatrix3Prong2Pos1Neg[5]);
              }

              if (config.debug) {
                uint8_t prong3CutStatus[kN3ProngDecays];
                for (int iDecay3P = 0; iDecay3P < kN3ProngDecays; iDecay3P++) {
                  prong3CutStatus[iDecay3P] = nCutStatus3ProngBit[iDecay3P];
                  for (int iCut = 0; iCut < kNCuts3Prong[iDecay3P]; iCut++) {
                    if (!cutStatus3Prong[iDecay3P][iCut]) {
                      CLRBIT(prong3CutStatus[iDecay3P], iCut);
                    }
                  }
                }
                rowProng3CutStatus(prong3CutStatus[0], prong3CutStatus[1], prong3CutStatus[2], prong3CutStatus[3]); // FIXME when we can do this by looping over kN3ProngDecays
              }

              // fill histograms
              if (config.fillHistograms) {
                registry.fill(HIST("hVtx3ProngX"), secondaryVertex3[0]);
                registry.fill(HIST("hVtx3ProngY"), secondaryVertex3[1]);
                registry.fill(HIST("hVtx3ProngZ"), secondaryVertex3[2]);
                const std::array arr3Mom{pvec0, pvec1, pvec2};
                for (int iDecay3P = 0; iDecay3P < kN3ProngDecays; iDecay3P++) {
                  if (TESTBIT(isSelected3ProngCand, iDecay3P)) {
                    if (TESTBIT(whichHypo3Prong[iDecay3P], 0)) {
                      const auto mass3Prong = RecoDecay::m(arr3Mom, arrMass3Prong[iDecay3P][0]);
                      switch (iDecay3P) {
                        case hf_cand_3prong::DecayType::DplusToPiKPi:
                          registry.fill(HIST("hMassDPlusToPiKPi"), mass3Prong);
                          break;
                        case hf_cand_3prong::DecayType::DsToKKPi:
                          registry.fill(HIST("hMassDsToKKPi"), mass3Prong);
                          break;
                        case hf_cand_3prong::DecayType::LcToPKPi:
                          registry.fill(HIST("hMassLcToPKPi"), mass3Prong);
                          break;
                        case hf_cand_3prong::DecayType::XicToPKPi:
                          registry.fill(HIST("hMassXicToPKPi"), mass3Prong);
                          break;
                        case hf_cand_3prong::DecayType::CdToDeKPi:
                          registry.fill(HIST("hMassCdToDeKPi"), mass3Prong);
                          break;
                      }
                    }
                    if (TESTBIT(whichHypo3Prong[iDecay3P], 1)) {
                      const auto mass3Prong = RecoDecay::m(arr3Mom, arrMass3Prong[iDecay3P][1]);
                      switch (iDecay3P) {
                        case hf_cand_3prong::DecayType::DsToKKPi:
                          registry.fill(HIST("hMassDsToKKPi"), mass3Prong);
                          break;
                        case hf_cand_3prong::DecayType::LcToPKPi:
                          registry.fill(HIST("hMassLcToPKPi"), mass3Prong);
                          break;
                        case hf_cand_3prong::DecayType::XicToPKPi:
                          registry.fill(HIST("hMassXicToPKPi"), mass3Prong);
                          break;
                        case hf_cand_3prong::DecayType::CdToDeKPi:
                          registry.fill(HIST("hMassCdToDeKPi"), mass3Prong);
                          break;
                      }
                    }
                  }
                }
              }
            }

            // second loop over negative tracks
            for (auto trackIndexNeg2 = trackIndexNeg1 + 1; trackIndexNeg2 != groupedTrackIndicesNeg1.end(); ++trackIndexNeg2) {

              int isSelected3ProngCand = n3ProngBit;
              if (!TESTBIT(trackIndexNeg2.isSelProng(), CandidateType::Cand3Prong)) { // continue immediately
                if (!config.debug) {
                  continue;
                }
                isSelected3ProngCand = 0;
              }

              if (config.applyKaonPidIn3Prongs && !TESTBIT(trackIndexPos1.isIdentifiedPid(), ChannelKaonPid)) { // continue immediately if kaon PID enabled and opposite-sign track not a kaon
                if (!config.debug) {
                  continue;
                }
                isSelected3ProngCand = 0;
              }

              auto trackNeg2 = trackIndexNeg2.template track_as<TTracks>();
              auto trackParVarNeg2 = getTrackParCov(trackNeg2);
              std::array dcaInfoNeg2{trackNeg2.dcaXY(), trackNeg2.dcaZ()};

              // preselection of 3-prong candidates
              if (isSelected3ProngCand) {
                std::array pVecTrackNeg2{trackNeg2.pVector()};
                if (thisCollId != trackNeg2.collisionId()) { // this is not the "default" collision for this track and we still did not re-propagate it, we have to re-propagate it
                  o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParVarNeg2, 2.f, noMatCorr, &dcaInfoNeg2);
                  getPxPyPz(trackParVarNeg2, pVecTrackNeg2);
                }

                if (config.debug) {
                  for (int iDecay3P = 0; iDecay3P < kN3ProngDecays; iDecay3P++) {
                    for (int iCut = 0; iCut < kNCuts3Prong[iDecay3P]; iCut++) {
                      cutStatus3Prong[iDecay3P][iCut] = true;
                    }
                  }
                }

                // 3-prong preselections
                int8_t const isIdentifiedPidTrackNeg1 = trackIndexNeg1.isIdentifiedPid();
                int8_t const isIdentifiedPidTrackNeg2 = trackIndexNeg2.isIdentifiedPid();
                applyPreselection3Prong(pVecTrackNeg1, pVecTrackPos1, pVecTrackNeg2, isIdentifiedPidTrackNeg1, isIdentifiedPidTrackNeg2, cutStatus3Prong, whichHypo3Prong, isSelected3ProngCand);
                if (!config.debug && isSelected3ProngCand == 0) {
                  continue;
                }
              }

              /// PV refit excluding the candidate daughters, if contributors
              std::array pvRefitCoord3Prong1Pos2Neg{collision.posX(), collision.posY(), collision.posZ()}; /// initialize to the original PV
              std::array pvRefitCovMatrix3Prong1Pos2Neg{getPrimaryVertex(collision).getCov()};             /// initialize to the original PV
              if constexpr (DoPvRefit) {
                if (config.fillHistograms) {
                  registry.fill(HIST("PvRefit/verticesPerCandidate"), 1);
                }
                int nCandContr = 3;
                auto trackFirstIt = std::find(vecPvContributorGlobId.begin(), vecPvContributorGlobId.end(), trackPos1.globalIndex());
                auto trackSecondIt = std::find(vecPvContributorGlobId.begin(), vecPvContributorGlobId.end(), trackNeg1.globalIndex());
                auto trackThirdIt = std::find(vecPvContributorGlobId.begin(), vecPvContributorGlobId.end(), trackNeg2.globalIndex());
                bool isTrackFirstContr = true;
                bool isTrackSecondContr = true;
                bool isTrackThirdContr = true;
                if (trackFirstIt == vecPvContributorGlobId.end()) {
                  /// This track did not contribute to the original PV refit
                  if (config.debugPvRefit) {
                    LOG(info) << "--- [3 prong] trackPos1 with globalIndex " << trackPos1.globalIndex() << " was not a PV contributor";
                  }
                  nCandContr--;
                  isTrackFirstContr = false;
                }
                if (trackSecondIt == vecPvContributorGlobId.end()) {
                  /// This track did not contribute to the original PV refit
                  if (config.debugPvRefit) {
                    LOG(info) << "--- [3 prong] trackNeg1 with globalIndex " << trackNeg1.globalIndex() << " was not a PV contributor";
                  }
                  nCandContr--;
                  isTrackSecondContr = false;
                }
                if (trackThirdIt == vecPvContributorGlobId.end()) {
                  /// This track did not contribute to the original PV refit
                  if (config.debugPvRefit) {
                    LOG(info) << "--- [3 prong] trackNeg2 with globalIndex " << trackNeg2.globalIndex() << " was not a PV contributor";
                  }
                  nCandContr--;
                  isTrackThirdContr = false;
                }

                // Fill a vector with global ID of candidate daughters that are contributors
                std::vector<int64_t> vecCandPvContributorGlobId = {};
                if (isTrackFirstContr) {
                  vecCandPvContributorGlobId.push_back(trackPos1.globalIndex());
                }
                if (isTrackSecondContr) {
                  vecCandPvContributorGlobId.push_back(trackNeg1.globalIndex());
                }
                if (isTrackThirdContr) {
                  vecCandPvContributorGlobId.push_back(trackNeg2.globalIndex());
                }

                if (nCandContr == 3 || nCandContr == 2) { // o2-linter: disable="magic-number" (see comment below)
                  /// At least two of the daughter tracks were used for the original PV refit, let's refit it after excluding them
                  if (config.debugPvRefit) {
                    LOG(info) << "### [3 prong] Calling performPvRefitCandProngs for HF 3 prong candidate, removing " << nCandContr << " daughters";
                  }
                  performPvRefitCandProngs(collision, bcWithTimeStamps, vecPvContributorGlobId, vecPvContributorTrackParCov, vecCandPvContributorGlobId, pvRefitCoord3Prong1Pos2Neg, pvRefitCovMatrix3Prong1Pos2Neg);
                } else if (nCandContr == 1) {
                  /// Only one daughter was a contributor, let's use then the PV recalculated by excluding only it
                  if (config.debugPvRefit) {
                    LOG(info) << "####### [3 Prong] nCandContr==" << nCandContr << " ---> just 1 contributor!";
                  }
                  if (config.fillHistograms) {
                    registry.fill(HIST("PvRefit/verticesPerCandidate"), 5);
                  }
                  if (isTrackFirstContr && !isTrackSecondContr && !isTrackThirdContr) {
                    /// the first daughter is contributor, the second and the third are not
                    pvRefitCoord3Prong1Pos2Neg = {trackPos1.pvRefitX(), trackPos1.pvRefitY(), trackPos1.pvRefitZ()};
                    pvRefitCovMatrix3Prong1Pos2Neg = {trackPos1.pvRefitSigmaX2(), trackPos1.pvRefitSigmaXY(), trackPos1.pvRefitSigmaY2(), trackPos1.pvRefitSigmaXZ(), trackPos1.pvRefitSigmaYZ(), trackPos1.pvRefitSigmaZ2()};
                  } else if (!isTrackFirstContr && isTrackSecondContr && !isTrackThirdContr) {
                    /// the second daughter is contributor, the first and the third are not
                    pvRefitCoord3Prong1Pos2Neg = {trackNeg1.pvRefitX(), trackNeg1.pvRefitY(), trackNeg1.pvRefitZ()};
                    pvRefitCovMatrix3Prong1Pos2Neg = {trackNeg1.pvRefitSigmaX2(), trackNeg1.pvRefitSigmaXY(), trackNeg1.pvRefitSigmaY2(), trackNeg1.pvRefitSigmaXZ(), trackNeg1.pvRefitSigmaYZ(), trackNeg1.pvRefitSigmaZ2()};
                  } else if (!isTrackFirstContr && !isTrackSecondContr && isTrackThirdContr) {
                    /// the third daughter is contributor, the first and the second are not
                    pvRefitCoord3Prong1Pos2Neg = {trackNeg2.pvRefitX(), trackNeg2.pvRefitY(), trackNeg2.pvRefitZ()};
                    pvRefitCovMatrix3Prong1Pos2Neg = {trackNeg2.pvRefitSigmaX2(), trackNeg2.pvRefitSigmaXY(), trackNeg2.pvRefitSigmaY2(), trackNeg2.pvRefitSigmaXZ(), trackNeg2.pvRefitSigmaYZ(), trackNeg2.pvRefitSigmaZ2()};
                  }
                } else {
                  /// 0 contributors among the HF candidate daughters
                  if (config.fillHistograms) {
                    registry.fill(HIST("PvRefit/verticesPerCandidate"), 6);
                  }
                  if (config.debugPvRefit) {
                    LOG(info) << "####### [3 prong] nCandContr==" << nCandContr << " ---> some of the candidate daughters did not contribute to the original PV fit, PV refit not redone";
                  }
                }
              }

              // reconstruct the 3-prong secondary vertex
              int nVtxFrom3ProngFitterSecondLoop = 0;
              try {
                nVtxFrom3ProngFitterSecondLoop = df3.process(trackParVarNeg1, trackParVarPos1, trackParVarNeg2);
              } catch (...) {
                continue;
              }

              if (nVtxFrom3ProngFitterSecondLoop == 0) {
                continue;
              }
              // get secondary vertex
              const auto& secondaryVertex3 = df3.getPCACandidate();
              // get track momenta
              std::array<float, 3> pvec0{};
              std::array<float, 3> pvec1{};
              std::array<float, 3> pvec2{};
              const auto trackParVarPcaNeg1 = df3.getTrack(0);
              const auto trackParVarPcaPos1 = df3.getTrack(1);
              const auto trackParVarPcaNeg2 = df3.getTrack(2);
              trackParVarPcaNeg1.getPxPyPzGlo(pvec0);
              trackParVarPcaPos1.getPxPyPzGlo(pvec1);
              trackParVarPcaNeg2.getPxPyPzGlo(pvec2);

              const auto pVecCandProng3Neg = RecoDecay::pVec(pvec0, pvec1, pvec2);

              // 3-prong selections after secondary vertex
              applySelection3Prong(pVecCandProng3Neg, secondaryVertex3, pvRefitCoord3Prong1Pos2Neg, cutStatus3Prong, isSelected3ProngCand);

              std::array<std::vector<float>, kN3ProngDecays - 1> mlScores3Prongs{};
              if (config.applyMlForHfFilters) {
                const std::vector<float> inputFeatures{trackParVarPcaNeg1.getPt(), dcaInfoNeg1[0], dcaInfoNeg1[1], trackParVarPcaPos1.getPt(), dcaInfoPos1[0], dcaInfoPos1[1], trackParVarPcaNeg2.getPt(), dcaInfoNeg2[0], dcaInfoNeg2[1]};
                std::vector<float> inputFeaturesLcPid{};
                if constexpr (UsePidForHfFiltersBdt) {
                  inputFeaturesLcPid.push_back(trackNeg1.tpcNSigmaPr());
                  inputFeaturesLcPid.push_back(trackNeg2.tpcNSigmaPr());
                  inputFeaturesLcPid.push_back(trackNeg1.tpcNSigmaPi());
                  inputFeaturesLcPid.push_back(trackNeg2.tpcNSigmaPi());
                  inputFeaturesLcPid.push_back(trackPos1.tpcNSigmaKa());
                }
                applyMlSelectionForHfFilters3Prong<UsePidForHfFiltersBdt>(inputFeatures, inputFeaturesLcPid, mlScores3Prongs, isSelected3ProngCand);
              }

              if (!config.debug && isSelected3ProngCand == 0) {
                continue;
              }

              // fill table row
              rowTrackIndexProng3(thisCollId, trackNeg1.globalIndex(), trackPos1.globalIndex(), trackNeg2.globalIndex(), isSelected3ProngCand);
              if (config.applyMlForHfFilters) {
                rowTrackIndexMlScoreProng3(mlScores3Prongs[0], mlScores3Prongs[1], mlScores3Prongs[2], mlScores3Prongs[3]);
              }
              if constexpr (DoPvRefit) {
                // fill table row of coordinates of PV refit
                rowProng3PVrefit(pvRefitCoord3Prong1Pos2Neg[0], pvRefitCoord3Prong1Pos2Neg[1], pvRefitCoord3Prong1Pos2Neg[2],
                                 pvRefitCovMatrix3Prong1Pos2Neg[0], pvRefitCovMatrix3Prong1Pos2Neg[1], pvRefitCovMatrix3Prong1Pos2Neg[2], pvRefitCovMatrix3Prong1Pos2Neg[3], pvRefitCovMatrix3Prong1Pos2Neg[4], pvRefitCovMatrix3Prong1Pos2Neg[5]);
              }

              if (config.debug) {
                int prong3CutStatus[kN3ProngDecays];
                for (int iDecay3P = 0; iDecay3P < kN3ProngDecays; iDecay3P++) {
                  prong3CutStatus[iDecay3P] = nCutStatus3ProngBit[iDecay3P];
                  for (int iCut = 0; iCut < kNCuts3Prong[iDecay3P]; iCut++) {
                    if (!cutStatus3Prong[iDecay3P][iCut]) {
                      CLRBIT(prong3CutStatus[iDecay3P], iCut);
                    }
                  }
                }
                rowProng3CutStatus(prong3CutStatus[0], prong3CutStatus[1], prong3CutStatus[2], prong3CutStatus[3]); // FIXME when we can do this by looping over kN3ProngDecays
              }

              // fill histograms
              if (config.fillHistograms) {
                registry.fill(HIST("hVtx3ProngX"), secondaryVertex3[0]);
                registry.fill(HIST("hVtx3ProngY"), secondaryVertex3[1]);
                registry.fill(HIST("hVtx3ProngZ"), secondaryVertex3[2]);
                const std::array arr3Mom{pvec0, pvec1, pvec2};
                for (int iDecay3P = 0; iDecay3P < kN3ProngDecays; iDecay3P++) {
                  if (TESTBIT(isSelected3ProngCand, iDecay3P)) {
                    if (TESTBIT(whichHypo3Prong[iDecay3P], 0)) {
                      const auto mass3Prong = RecoDecay::m(arr3Mom, arrMass3Prong[iDecay3P][0]);
                      switch (iDecay3P) {
                        case hf_cand_3prong::DecayType::DplusToPiKPi:
                          registry.fill(HIST("hMassDPlusToPiKPi"), mass3Prong);
                          break;
                        case hf_cand_3prong::DecayType::DsToKKPi:
                          registry.fill(HIST("hMassDsToKKPi"), mass3Prong);
                          break;
                        case hf_cand_3prong::DecayType::LcToPKPi:
                          registry.fill(HIST("hMassLcToPKPi"), mass3Prong);
                          break;
                        case hf_cand_3prong::DecayType::XicToPKPi:
                          registry.fill(HIST("hMassXicToPKPi"), mass3Prong);
                          break;
                        case hf_cand_3prong::DecayType::CdToDeKPi:
                          registry.fill(HIST("hMassCdToDeKPi"), mass3Prong);
                          break;
                      }
                    }
                    if (TESTBIT(whichHypo3Prong[iDecay3P], 1)) {
                      const auto mass3Prong = RecoDecay::m(arr3Mom, arrMass3Prong[iDecay3P][1]);
                      switch (iDecay3P) {
                        case hf_cand_3prong::DecayType::DsToKKPi:
                          registry.fill(HIST("hMassDsToKKPi"), mass3Prong);
                          break;
                        case hf_cand_3prong::DecayType::LcToPKPi:
                          registry.fill(HIST("hMassLcToPKPi"), mass3Prong);
                          break;
                        case hf_cand_3prong::DecayType::XicToPKPi:
                          registry.fill(HIST("hMassXicToPKPi"), mass3Prong);
                          break;
                        case hf_cand_3prong::DecayType::CdToDeKPi:
                          registry.fill(HIST("hMassCdToDeKPi"), mass3Prong);
                          break;
                      }
                    }
                  }
                }
              }
            }
          }

          if (config.doDstar && TESTBIT(isSelected2ProngCand, hf_cand_2prong::DecayType::D0ToPiK) && (pt2Prong + config.ptTolerance) * 1.2 > config.binsPtDstarToD0Pi->at(0) && whichHypo2Prong[kN2ProngDecays] != 0) { // o2-linter: disable="magic-number" (see comment below)
                                                                                                                                                                                                                        // if D* enabled and pt of the D0 is larger than the minimum of the D* one within 20% (D* and D0 momenta are very similar, always within 20% according to PYTHIA8)
            // second loop over positive tracks
            if (TESTBIT(whichHypo2Prong[kN2ProngDecays], 0) && (!config.applyKaonPidIn3Prongs || TESTBIT(trackIndexNeg1.isIdentifiedPid(), ChannelKaonPid))) { // only for D0 candidates; moreover if kaon PID enabled, apply to the negative track
              auto groupedTrackIndicesSoftPionsPos = positiveSoftPions->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
              for (auto trackIndexPos2 = groupedTrackIndicesSoftPionsPos.begin(); trackIndexPos2 != groupedTrackIndicesSoftPionsPos.end(); ++trackIndexPos2) {
                if (trackIndexPos2 == trackIndexPos1) {
                  continue;
                }
                auto trackPos2 = trackIndexPos2.template track_as<TTracks>();
                std::array pVecTrackPos2{trackPos2.pVector()};
                if (thisCollId != trackPos2.collisionId()) { // this is not the "default" collision for this track, we have to re-propagate it
                  auto trackParVarPos2 = getTrackParCov(trackPos2);
                  std::array dcaInfoPos2{trackPos2.dcaXY(), trackPos2.dcaZ()};
                  o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParVarPos2, 2.f, noMatCorr, &dcaInfoPos2);
                  getPxPyPz(trackParVarPos2, pVecTrackPos2);
                }

                uint8_t isSelectedDstar{0};
                uint8_t cutStatus{BIT(kNCutsDstar) - 1};
                float deltaMass{-1.};
                isSelectedDstar = applySelectionDstar(pVecTrackPos1, pVecTrackNeg1, pVecTrackPos2, cutStatus, deltaMass); // we do not compute the D* decay vertex at this stage because we are not interested in applying topological selections
                if (isSelectedDstar) {
                  rowTrackIndexDstar(thisCollId, trackPos2.globalIndex(), lastFilledD0);
                  if (config.fillHistograms) {
                    registry.fill(HIST("hMassDstarToD0Pi"), deltaMass);
                  }
                  if constexpr (DoPvRefit) {
                    // fill table row with coordinates of PV refit (same as 2-prong because we do not remove the soft pion)
                    rowDstarPVrefit(pvRefitCoord2Prong[0], pvRefitCoord2Prong[1], pvRefitCoord2Prong[2],
                                    pvRefitCovMatrix2Prong[0], pvRefitCovMatrix2Prong[1], pvRefitCovMatrix2Prong[2], pvRefitCovMatrix2Prong[3], pvRefitCovMatrix2Prong[4], pvRefitCovMatrix2Prong[5]);
                  }
                }
                if (config.debug) {
                  rowDstarCutStatus(cutStatus);
                }
              }
            }

            // second loop over negative tracks
            if (TESTBIT(whichHypo2Prong[kN2ProngDecays], 1) && (!config.applyKaonPidIn3Prongs || TESTBIT(trackIndexPos1.isIdentifiedPid(), ChannelKaonPid))) { // only for D0bar candidates; moreover if kaon PID enabled, apply to the positive track
              auto groupedTrackIndicesSoftPionsNeg = negativeSoftPions->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
              for (auto trackIndexNeg2 = groupedTrackIndicesSoftPionsNeg.begin(); trackIndexNeg2 != groupedTrackIndicesSoftPionsNeg.end(); ++trackIndexNeg2) {
                if (trackIndexNeg1 == trackIndexNeg2) {
                  continue;
                }
                auto trackNeg2 = trackIndexNeg2.template track_as<TTracks>();
                std::array pVecTrackNeg2{trackNeg2.pVector()};
                if (thisCollId != trackNeg2.collisionId()) { // this is not the "default" collision for this track, we have to re-propagate it
                  auto trackParVarNeg2 = getTrackParCov(trackNeg2);
                  std::array dcaInfoNeg2{trackNeg2.dcaXY(), trackNeg2.dcaZ()};
                  o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParVarNeg2, 2.f, noMatCorr, &dcaInfoNeg2);
                  getPxPyPz(trackParVarNeg2, pVecTrackNeg2);
                }

                uint8_t isSelectedDstar{0};
                uint8_t cutStatus{BIT(kNCutsDstar) - 1};
                float deltaMass{-1.};
                isSelectedDstar = applySelectionDstar(pVecTrackNeg1, pVecTrackPos1, pVecTrackNeg2, cutStatus, deltaMass); // we do not compute the D* decay vertex at this stage because we are not interested in applying topological selections
                if (isSelectedDstar) {
                  rowTrackIndexDstar(thisCollId, trackNeg2.globalIndex(), lastFilledD0);
                  if (config.fillHistograms) {
                    registry.fill(HIST("hMassDstarToD0Pi"), deltaMass);
                  }
                  if constexpr (DoPvRefit) {
                    // fill table row with coordinates of PV refit (same as 2-prong because we do not remove the soft pion)
                    rowDstarPVrefit(pvRefitCoord2Prong[0], pvRefitCoord2Prong[1], pvRefitCoord2Prong[2],
                                    pvRefitCovMatrix2Prong[0], pvRefitCovMatrix2Prong[1], pvRefitCovMatrix2Prong[2], pvRefitCovMatrix2Prong[3], pvRefitCovMatrix2Prong[4], pvRefitCovMatrix2Prong[5]);
                  }
                }
                if (config.debug) {
                  rowDstarCutStatus(cutStatus);
                }
              }
            }
          } // end of D*
        }
      }

      const int nTracks = 0;
      // auto nTracks = trackIndicesPerCollision.lastIndex() - trackIndicesPerCollision.firstIndex(); // number of tracks passing 2 and 3 prong selection in this collision
      nCand2 = rowTrackIndexProng2.lastIndex() - nCand2; // number of 2-prong candidates in this collision
      nCand3 = rowTrackIndexProng3.lastIndex() - nCand3; // number of 3-prong candidates in this collision

      if (config.fillHistograms) {
        registry.fill(HIST("hNTracks"), nTracks);
        registry.fill(HIST("hNCand2Prong"), nCand2);
        registry.fill(HIST("hNCand3Prong"), nCand3);
        registry.fill(HIST("hNCand2ProngVsNTracks"), nTracks, nCand2);
        registry.fill(HIST("hNCand3ProngVsNTracks"), nTracks, nCand3);
      }
    }
  } /// end of run2And3Prongs function

  void processNo2And3Prongs(SelectedCollisions const&)
  {
    // dummy
  }
  PROCESS_SWITCH(HfTrackIndexSkimCreator, processNo2And3Prongs, "Do not process 2-prongs and 3-prongs", false);

  void process2And3ProngsWithPvRefit( // soa::Join<aod::Collisions, aod::CentV0Ms>::iterator const& collision, //FIXME add centrality when option for variations to the process function appears
    SelectedCollisions const& collisions,
    aod::BCsWithTimestamps const& bcWithTimeStamps,
    FilteredTrackAssocSel const& trackIndices,
    TracksWithPVRefitAndDCA const& tracks)
  {
    run2And3Prongs<true, false>(collisions, bcWithTimeStamps, trackIndices, tracks);
  }
  PROCESS_SWITCH(HfTrackIndexSkimCreator, process2And3ProngsWithPvRefit, "Process 2-prong and 3-prong skim with PV refit", false);

  void process2And3ProngsNoPvRefit( // soa::Join<aod::Collisions, aod::CentV0Ms>::iterator const& collision, //FIXME add centrality when option for variations to the process function appears
    SelectedCollisions const& collisions,
    aod::BCsWithTimestamps const& bcWithTimeStamps,
    FilteredTrackAssocSel const& trackIndices,
    aod::TracksWCovDcaExtra const& tracks)
  {
    run2And3Prongs<false, false>(collisions, bcWithTimeStamps, trackIndices, tracks);
  }
  PROCESS_SWITCH(HfTrackIndexSkimCreator, process2And3ProngsNoPvRefit, "Process 2-prong and 3-prong skim without PV refit", true);

  void process2And3ProngsWithPvRefitWithPidForHfFiltersBdt( // soa::Join<aod::Collisions, aod::CentV0Ms>::iterator const& collision, //FIXME add centrality when option for variations to the process function appears
    SelectedCollisions const& collisions,
    aod::BCsWithTimestamps const& bcWithTimeStamps,
    FilteredTrackAssocSel const& trackIndices,
    soa::Join<TracksWithPVRefitAndDCA, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr> const& tracks)
  {
    run2And3Prongs<true, true>(collisions, bcWithTimeStamps, trackIndices, tracks);
  }
  PROCESS_SWITCH(HfTrackIndexSkimCreator, process2And3ProngsWithPvRefitWithPidForHfFiltersBdt, "Process 2-prong and 3-prong skim with PV refit and PID for software trigger BDTs (Lc and Xic only)", false);

  void process2And3ProngsNoPvRefitWithPidForHfFiltersBdt( // soa::Join<aod::Collisions, aod::CentV0Ms>::iterator const& collision, //FIXME add centrality when option for variations to the process function appears
    SelectedCollisions const& collisions,
    aod::BCsWithTimestamps const& bcWithTimeStamps,
    FilteredTrackAssocSel const& trackIndices,
    soa::Join<aod::TracksWCovDcaExtra, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr> const& tracks)
  {
    run2And3Prongs<false, true>(collisions, bcWithTimeStamps, trackIndices, tracks);
  }
  PROCESS_SWITCH(HfTrackIndexSkimCreator, process2And3ProngsNoPvRefitWithPidForHfFiltersBdt, "Process 2-prong and 3-prong skim without PV refit and PID for software trigger BDTs (Lc and Xic only)", false);
};

//________________________________________________________________________________________________________________________

/// Pre-selection of cascade secondary vertices
/// It will produce in any case a Hf2Prongs object, but mixing a V0
/// with a track, instead of 2 tracks

/// to run: o2-analysis-weak-decay-indices --aod-file AO2D.root -b | o2-analysis-lambdakzerobuilder -b |
///         o2-analysis-trackextension -b | o2-analysis-hf-track-index-skim-creator -b

struct HfTrackIndexSkimCreatorCascades {
  Produces<aod::HfCascades> rowTrackIndexCasc;

  struct : ConfigurableGroup {
    double etaMinDefault{-99999.};
    Configurable<bool> isRun2{"isRun2", false, "enable Run 2 or Run 3 GRP objects for magnetic field"};
    Configurable<bool> fillHistograms{"fillHistograms", true, "fill histograms"};
    // vertexing
    Configurable<bool> useDCAFitter{"useDCAFitter", true, "flag to optionally turn on/off the vertex reconstruction with the DCAFitter"};
    Configurable<bool> propagateToPCA{"propagateToPCA", true, "create tracks version propagated to PCA"};
    Configurable<double> maxR{"maxR", 200., "reject PCA's above this radius"};
    Configurable<double> maxDZIni{"maxDZIni", 4., "reject (if>0) PCA candidate if tracks DZ exceeds threshold"};
    Configurable<double> minParamChange{"minParamChange", 1.e-3, "stop iterations if largest change of any X is smaller than this"};
    Configurable<double> minRelChi2Change{"minRelChi2Change", 0.9, "stop iterations if chi2/chi2old > this"};
    Configurable<bool> useAbsDCA{"useAbsDCA", true, "Minimise abs. distance rather than chi2"};
    Configurable<bool> useWeightedFinalPCA{"useWeightedFinalPCA", true, "Recalculate vertex position using track covariances, effective only if useAbsDCA is true"};
    // track cuts for V0 daughters
    Configurable<double> etaMinV0Daugh{"etaMinV0Daugh", std::forward<double>(etaMinDefault), "min. pseudorapidity V0 daughters"};
    Configurable<double> etaMaxV0Daugh{"etaMaxV0Daugh", 1.1, "max. pseudorapidity V0 daughters"};
    Configurable<double> ptMinV0Daugh{"ptMinV0Daugh", 0.05, "min. pT V0 daughters"};
    // v0 cuts
    Configurable<double> cpaV0Min{"cpaV0Min", 0.95, "min. cos PA V0"}; // as in the task that create the V0s
    Configurable<double> cutInvMassV0{"cutInvMassV0", 0.05, "V0 candidate invariant mass difference wrt PDG"};
    // cascade cuts
    Configurable<double> ptCascCandMin{"ptCascCandMin", -1., "min. pT of the cascade candidate"};                    // PbPb 2018: use 1
    Configurable<double> cutInvMassCascLc{"cutInvMassCascLc", 1., "Lc candidate invariant mass difference wrt PDG"}; // for PbPb 2018: use 0.2
    // Configurable<double> cutCascDCADaughters{"cutCascDCADaughters", .1, "DCA between V0 and bachelor in cascade"};
    // proton PID
    Configurable<bool> applyProtonPid{"applyProtonPid", false, "Apply proton PID for Lc->pK0S"};

    // CCDB
    Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
    Configurable<std::string> ccdbPathLut{"ccdbPathLut", "GLO/Param/MatLUT", "Path for LUT parametrization"};
    Configurable<std::string> ccdbPathGrp{"ccdbPathGrp", "GLO/GRP/GRP", "Path of the grp file (Run 2)"};
    Configurable<std::string> ccdbPathGrpMag{"ccdbPathGrpMag", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object (Run 3)"};
  } config;

  o2::vertexing::DCAFitterN<2> df2; // 2-prong vertex fitter
  // Needed for PV refitting
  Service<o2::ccdb::BasicCCDBManager> ccdb{};
  o2::base::MatLayerCylSet* lut{};
  o2::base::Propagator::MatCorrType noMatCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;
  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;
  int runNumber{0};

  using SelectedCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::HfSelCollision>>;
  using FilteredTrackAssocSel = soa::Filtered<soa::Join<aod::TrackAssoc, aod::HfSelTrack>>;

  Filter filterSelectCollisions = (aod::hf_sel_collision::whyRejectColl == static_cast<o2::hf_evsel::HfCollisionRejectionMask>(0));
  Filter filterSelectTrackIds = (aod::hf_sel_track::isSelProng & static_cast<uint32_t>(BIT(CandidateType::CandV0bachelor))) != 0u && (!config.applyProtonPid || (aod::hf_sel_track::isIdentifiedPid & static_cast<uint32_t>(BIT(ChannelsProtonPid::LcToPK0S))) != 0u);

  Preslice<FilteredTrackAssocSel> trackIndicesPerCollision = aod::track_association::collisionId;
  Preslice<aod::V0Datas> v0sPerCollision = aod::v0data::collisionId;

  // histograms
  HistogramRegistry registry{"registry"};

  void init(InitContext const&)
  {
    if (!doprocessCascades) {
      return;
    }

    if (config.etaMinV0Daugh == config.etaMinDefault) {
      config.etaMinV0Daugh.value = -config.etaMaxV0Daugh;
    }

    if (config.useDCAFitter) {
      df2.setPropagateToPCA(config.propagateToPCA);
      df2.setMaxR(config.maxR);
      df2.setMinParamChange(config.minParamChange);
      df2.setMinRelChi2Change(config.minRelChi2Change);
      df2.setMaxDZIni(config.maxDZIni);
      df2.setUseAbsDCA(config.useAbsDCA);
      df2.setWeightedFinalPCA(config.useWeightedFinalPCA);
    }

    ccdb->setURL(config.ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>(config.ccdbPathLut));

    if (config.fillHistograms) {
      registry.add("hVtx2ProngX", "2-prong candidates;#it{x}_{sec. vtx.} (cm);entries", {HistType::kTH1D, {{1000, -2., 2.}}});
      registry.add("hVtx2ProngY", "2-prong candidates;#it{y}_{sec. vtx.} (cm);entries", {HistType::kTH1D, {{1000, -2., 2.}}});
      registry.add("hVtx2ProngZ", "2-prong candidates;#it{z}_{sec. vtx.} (cm);entries", {HistType::kTH1D, {{1000, -2., 2.}}});
      registry.add("hMassLcToPK0S", "#Lambda_{c}^ candidates;inv. mass (p K_{S}^{0}) (GeV/#it{c}^{2});entries", {HistType::kTH1D, {{500, 0., 5.}}});
    }
  }

  void processNoCascades(SelectedCollisions const&)
  {
    // dummy
  }
  PROCESS_SWITCH(HfTrackIndexSkimCreatorCascades, processNoCascades, "Do not skim HF -> V0 cascades", true);

  void processCascades(SelectedCollisions const& collisions,
                       soa::Join<aod::V0Datas, aod::V0Covs> const& v0s,
                       FilteredTrackAssocSel const& trackIndices,
                       aod::TracksWCovDcaExtra const&,
                       aod::BCsWithTimestamps const&)
  {
    // set the magnetic field from CCDB
    for (const auto& collision : collisions) {
      const auto bc = collision.bc_as<o2::aod::BCsWithTimestamps>();
      initCCDB(bc, runNumber, ccdb, config.isRun2 ? config.ccdbPathGrp : config.ccdbPathGrpMag, lut, config.isRun2);
      if (config.useDCAFitter) {
        df2.setBz(o2::base::Propagator::Instance()->getNominalBz());
        df2.setMatCorrType(matCorr);
      }

      const auto thisCollId = collision.globalIndex();
      const auto groupedBachTrackIndices = trackIndices.sliceBy(trackIndicesPerCollision, thisCollId);
      const auto groupedV0s = v0s.sliceBy(v0sPerCollision, thisCollId);

      // fist we loop over the bachelor candidate
      for (const auto& bachIdx : groupedBachTrackIndices) {

        const auto bach = bachIdx.track_as<aod::TracksWCovDcaExtra>();
        std::array pVecBach{bach.pVector()};
        auto trackBach = getTrackParCov(bach);
        if (thisCollId != bach.collisionId()) { // this is not the "default" collision for this track, we have to re-propagate it
          std::array dcaInfoBach{bach.dcaXY(), bach.dcaZ()};
          o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackBach, 2.f, noMatCorr, &dcaInfoBach);
          getPxPyPz(trackBach, pVecBach);
        }

        // now we loop over the V0s
        for (const auto& v0 : groupedV0s) {
          // selections on the V0 daughters
          const auto& trackV0DaughPos = v0.posTrack_as<aod::TracksWCovDcaExtra>(); // only used for indices and track cuts (TPC clusters, TPC refit)
          const auto& trackV0DaughNeg = v0.negTrack_as<aod::TracksWCovDcaExtra>(); // only used for indices and track cuts (TPC clusters, TPC refit)

          // check not to take the same track twice (as bachelor and V0 daughter)
          if (trackV0DaughPos.globalIndex() == bach.globalIndex() || trackV0DaughNeg.globalIndex() == bach.globalIndex()) {
            continue;
          }

          const std::array pVecPos{v0.pxpos(), v0.pypos(), v0.pzpos()};
          const std::array pVecNeg{v0.pxneg(), v0.pyneg(), v0.pzneg()};

          const float ptPos = RecoDecay::pt(pVecPos);
          const float ptNeg = RecoDecay::pt(pVecNeg);
          if (ptPos < config.ptMinV0Daugh || // to the filters? I can't for now, it is not in the tables
              ptNeg < config.ptMinV0Daugh) {
            continue;
          }

          const float etaPos = RecoDecay::eta(pVecPos);
          const float etaNeg = RecoDecay::eta(pVecNeg);
          if ((etaPos > config.etaMaxV0Daugh || etaPos < config.etaMinV0Daugh) || // to the filters? I can't for now, it is not in the tables
              (etaNeg > config.etaMaxV0Daugh || etaNeg < config.etaMinV0Daugh)) {
            continue;
          }

          // V0 invariant mass selection
          if (std::abs(v0.mK0Short() - MassK0Short) > config.cutInvMassV0) {
            continue; // should go to the filter, but since it is a dynamic column, I cannot use it there
          }

          // V0 cosPointingAngle selection
          if (v0.v0cosPA() < config.cpaV0Min) {
            continue;
          }

          std::array pVecV0{v0.px(), v0.py(), v0.pz()};

          // invariant-mass cut: we do it here, before updating the momenta of bach and V0 during the fitting to save CPU
          // TODO: but one should better check that the value here and after the fitter do not change significantly!!!
          double mass2K0sP = RecoDecay::m(std::array{pVecBach, pVecV0}, std::array{MassProton, MassK0Short});
          if ((config.cutInvMassCascLc >= 0.) && (std::abs(mass2K0sP - MassLambdaCPlus) > config.cutInvMassCascLc)) {
            continue;
          }

          // now we find the DCA between the V0 and the bachelor, for the cascade
          if (config.useDCAFitter) {

            const std::array vertexV0{v0.x(), v0.y(), v0.z()};
            // we build the neutral track to then build the cascade
            std::array<float, 21> covV{};
            constexpr std::size_t NIndicesMom{6u};
            constexpr std::size_t MomInd[NIndicesMom] = {9, 13, 14, 18, 19, 20}; // cov matrix elements for momentum component
            for (std::size_t i = 0; i < NIndicesMom; i++) {
              covV[MomInd[i]] = v0.momentumCovMat()[i];
              covV[i] = v0.positionCovMat()[i];
            }
            auto trackV0 = o2::track::TrackParCov(vertexV0, pVecV0, covV, 0, true);
            trackV0.setAbsCharge(0);
            trackV0.setPID(o2::track::PID::K0);

            int nCand2 = 0;
            try {
              nCand2 = df2.process(trackV0, trackBach);
            } catch (...) {
              continue;
            }

            if (nCand2 == 0) {
              continue;
            }
            df2.propagateTracksToVertex();        // propagate the bach and V0 to the Lc vertex
            df2.getTrack(0).getPxPyPzGlo(pVecV0); // take the momentum at the Lc vertex
            df2.getTrack(1).getPxPyPzGlo(pVecBach);
          }

          // cascade candidate pT cut
          const auto ptCascCand = RecoDecay::pt(pVecBach, pVecV0);
          if (ptCascCand < config.ptCascCandMin) {
            continue;
          }

          // invariant mass
          // re-calculate invariant masses with updated momenta, to fill the histogram
          mass2K0sP = RecoDecay::m(std::array{pVecBach, pVecV0}, std::array{MassProton, MassK0Short});
          std::array posCasc{0., 0., 0.};
          if (config.useDCAFitter) {
            const auto& cascVtx = df2.getPCACandidate();
            for (int iCoord{0}; iCoord < 3; ++iCoord) { // o2-linter: disable="magic-number" ({x, y, z} coordinates})
              posCasc[iCoord] = cascVtx[iCoord];
            }
          }

          // fill table row
          rowTrackIndexCasc(thisCollId, bach.globalIndex(), v0.v0Id());
          // fill histograms
          if (config.fillHistograms) {
            registry.fill(HIST("hVtx2ProngX"), posCasc[0]);
            registry.fill(HIST("hVtx2ProngY"), posCasc[1]);
            registry.fill(HIST("hVtx2ProngZ"), posCasc[2]);
            registry.fill(HIST("hMassLcToPK0S"), mass2K0sP);
          }
        } // loop over V0s
      } // loop over tracks
    } // loop over collisions
  } // processCascades
  PROCESS_SWITCH(HfTrackIndexSkimCreatorCascades, processCascades, "Skim HF -> V0 cascades", false);
};

struct HfTrackIndexSkimCreatorLfCascades {
  Produces<aod::HfCascLf2Prongs> rowTrackIndexCasc2Prong;
  Produces<aod::HfCascLf3Prongs> rowTrackIndexCasc3Prong;

  struct : ConfigurableGroup {
    // whether to do or not validation plots
    Configurable<bool> fillHistograms{"fillHistograms", true, "fill histograms"};

    Configurable<bool> do3Prong{"do3Prong", false, "do 3-prong cascade"};
    Configurable<bool> rejDiffCollTrack{"rejDiffCollTrack", false, "Reject tracks coming from different collisions"};
    Configurable<double> ptTolerance{"ptTolerance", 0.1, "pT tolerance in GeV/c for applying preselections before vertex reconstruction"};

    // charm baryon invariant mass spectra limits
    Configurable<double> massXiPiMin{"massXiPiMin", 2.1, "Invariant mass lower limit for xi pi decay channel"};
    Configurable<double> massXiPiMax{"massXiPiMax", 3., "Invariant mass upper limit for xi pi decay channel"};
    Configurable<double> massOmegaCharmBachelorMin{"massOmegaCharmBachelorMin", 2.4, "Invariant mass lower limit for omega pi and omega k decay channel"};
    Configurable<double> massOmegaCharmBachelorMax{"massOmegaCharmBachelorMax", 3., "Invariant mass upper limit for omega pi and omega k decay channel"};
    Configurable<double> massXiPiPiMin{"massXiPiPiMin", 2.1, "Invariant mass lower limit for xi pi pi decay channel"};
    Configurable<double> massXiPiPiMax{"massXiPiPiMax", 2.8, "Invariant mass upper limit for xi pi pi decay channel"};

    // DCAFitter settings
    Configurable<bool> propagateToPCA{"propagateToPCA", false, "create tracks version propagated to PCA"};
    Configurable<bool> useAbsDCA{"useAbsDCA", true, "Minimise abs. distance rather than chi2"};
    Configurable<bool> useWeightedFinalPCA{"useWeightedFinalPCA", true, "Recalculate vertex position using track covariances, effective only if useAbsDCA is true"};
    Configurable<double> maxR{"maxR", 200., "reject PCA's above this radius"};
    Configurable<double> maxDZIni{"maxDZIni", 4., "reject (if>0) PCA candidate if tracks DZ exceeds threshold"};
    Configurable<double> maxDXYIni{"maxDXYIni", 4., "reject (if>0) PCA candidate if tracks DXY exceeds threshold"};
    Configurable<double> minParamChange{"minParamChange", 1.e-3, "stop iterations if largest change of any X is smaller than this"};
    Configurable<double> minRelChi2Change{"minRelChi2Change", 0.9, "stop iterations is chi2/chi2old > this"};
    Configurable<double> maxChi2{"maxChi2", 100., "discard vertices with chi2/Nprongs > this (or sum{DCAi^2}/Nprongs for abs. distance minimization)"};
    Configurable<bool> refitWithMatCorr{"refitWithMatCorr", true, "when doing propagateTracksToVertex, propagate tracks to vtx with material corrections and rerun minimization"};

    // Selection criteria
    // selections have been set to run2 lambda dedicated cuts
    // selections for cascade have been set to the loosest value between xi and omega
    // a tolerance has been added to be more conservative  ptMinOmegaczeroToOmegaKaLfCasc ptMinXicZeroOmegacZeroToXiPiLfCasc
    Configurable<float> ptMinOmegacZeroToOmegaPiLfCasc{"ptMinOmegacZeroToOmegaPiLfCasc", 0.f, "min. pT for Omegaczero in Omega + Pi decays"};
    Configurable<float> ptMinOmegaczeroToOmegaKaLfCasc{"ptMinOmegaczeroToOmegaKaLfCasc", 0.f, "min. pT for Omegaczero in Omega + Ka decays"};
    Configurable<float> ptMinXicZeroOmegacZeroToXiPiLfCasc{"ptMinXicZeroOmegacZeroToXiPiLfCasc", 0.f, "min. pT for XicZeroOmegacZero in Xi + Pi decays"};
    Configurable<float> ptMinXicplusLfCasc{"ptMinXicplusLfCasc", 0.f, "min. pT for Xicplus in Xi + Pi + Pi decays"};
    Configurable<float> v0TransvRadius{"v0TransvRadius", 1.f, "V0 radius in xy plane"};           // 1.2 (xi) and 1.1 (omega) in run2
    Configurable<float> cascTransvRadius{"cascTransvRadius", 0.4f, "Cascade radius in xy plane"}; // 0.5 cm (xi) and 0.6 (omega) in run2
    Configurable<float> decayLengthXicMin{"decayLengthXicMin", -1.f, "Min. decay length of Xic"}; // ...
    Configurable<float> dcaBachToPv{"dcaBachToPv", 0.03f, "DCA Bach To PV"};                      // 0.04 in run2
    Configurable<float> dcaV0ToPv{"dcaV0ToPv", 0.02f, "DCA V0 To PV"};                            // 0.03 in run2
    Configurable<double> v0CosPA{"v0CosPA", 0.95, "V0 CosPA"};                                    // 0.97 in run2 - KEEP LOSE to re-cut after PVRefit! - double -> N.B. dcos(x)/dx = 0 at x=0)
    Configurable<double> cascCosPA{"cascCosPA", 0.95, "Casc CosPA"};                              // 0.97 in run2 - KEEP LOSE to re-cut after PVRefit! - double -> N.B. dcos(x)/dx = 0 at x=0)
    Configurable<double> xicCosPA{"xicCosPA", 0.95, "Xic CosPA"};                                 // ...
    Configurable<float> dcaV0Dau{"dcaV0Dau", 2.f, "DCA V0 Daughters"};                            // conservative, a cut ar 1.0 should also be fine
    Configurable<float> dcaCascDau{"dcaCascDau", 2.f, "DCA Casc Daughters"};                      // conservative, a cut ar 1.0 should also be fine
    Configurable<float> dcaNegToPv{"dcaNegToPv", 0.05f, "DCA Neg To PV"};                         // 0.06 in run2
    Configurable<float> dcaPosToPv{"dcaPosToPv", 0.05f, "DCA Pos To PV"};                         // 0.06 in run2
    Configurable<float> v0MassWindow{"v0MassWindow", 0.01f, "V0 mass window"};                    // 0.008 in run2
    Configurable<float> cascadeMassWindow{"cascadeMassWindow", 0.01f, "Cascade mass window"};

    // magnetic field setting from CCDB
    Configurable<bool> isRun2{"isRun2", false, "enable Run 2 or Run 3 GRP objects for magnetic field"};
    Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
    Configurable<std::string> ccdbPathLut{"ccdbPathLut", "GLO/Param/MatLUT", "Path for LUT parametrization"};
    Configurable<std::string> ccdbPathGrp{"ccdbPathGrp", "GLO/GRP/GRP", "Path of the grp file (Run 2)"};
    Configurable<std::string> ccdbPathGrpMag{"ccdbPathGrpMag", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object (Run 3)"};
  } config;

  o2::vertexing::DCAFitterN<2> df2; // 2-prong vertex fitter
  Service<o2::ccdb::BasicCCDBManager> ccdb{};
  o2::base::MatLayerCylSet* lut{};
  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;
  int runNumber{};

  // array of PDG masses of possible charm baryon daughters
  static constexpr int kN2ProngDecays = hf_cand_casc_lf::DecayType2Prong::N2ProngDecays; // number of 2-prong hadron types
  static constexpr int kN3ProngDecays = hf_cand_casc_lf::DecayType3Prong::N3ProngDecays; // number of 3-prong hadron types
  std::array<std::array<double, 2>, kN2ProngDecays> arrMass2Prong{};
  std::array<std::array<double, 3>, kN3ProngDecays> arrMass3Prong{};

  using SelectedCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::HfSelCollision>>;
  using SelectedHfTrackAssoc = soa::Filtered<soa::Join<aod::TrackAssoc, aod::HfSelTrack>>;
  using CascFull = soa::Join<aod::CascDatas, aod::CascCovs>;
  using V0Full = soa::Join<aod::V0Datas, aod::V0Covs>;

  Filter filterSelectCollisions = (aod::hf_sel_collision::whyRejectColl == static_cast<o2::hf_evsel::HfCollisionRejectionMask>(0));
  Filter filterSelectTrackIds = (aod::hf_sel_track::isSelProng & static_cast<uint32_t>(BIT(CandidateType::CandCascadeBachelor))) != 0u;

  Preslice<aod::TracksWCovDca> tracksPerCollision = aod::track::collisionId;                     // needed for PV refit
  Preslice<SelectedHfTrackAssoc> trackIndicesPerCollision = aod::track_association::collisionId; // aod::hf_track_association::collisionId
  Preslice<CascFull> cascadesPerCollision = aod::cascdata::collisionId;

  // histograms
  HistogramRegistry registry{"registry"};

  void init(InitContext const&)
  {
    if (!doprocessLfCascades) {
      return;
    }

    arrMass2Prong[hf_cand_casc_lf::DecayType2Prong::XiczeroOmegaczeroToXiPi] = std::array{MassXiMinus, MassPiPlus};
    arrMass2Prong[hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaPi] = std::array{MassOmegaMinus, MassPiPlus};
    arrMass2Prong[hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaK] = std::array{MassOmegaMinus, MassKPlus};
    arrMass3Prong[hf_cand_casc_lf::DecayType3Prong::XicplusToXiPiPi] = std::array{MassXiMinus, MassPiPlus, MassPiPlus};

    df2.setPropagateToPCA(config.propagateToPCA);
    df2.setMaxR(config.maxR);
    df2.setMaxDZIni(config.maxDZIni);
    df2.setMinParamChange(config.minParamChange);
    df2.setMinRelChi2Change(config.minRelChi2Change);
    df2.setUseAbsDCA(config.useAbsDCA);
    df2.setWeightedFinalPCA(config.useWeightedFinalPCA);

    ccdb->setURL(config.ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>(config.ccdbPathLut));
    runNumber = 0;

    if (config.fillHistograms) {
      const AxisSpec ptAxis = {400, 0.0f, 20.0f, "it{p}_{T} (GeV/c)"};
      const AxisSpec massAxisXi = {200, 1.222f, 1.422f, "Inv. Mass (GeV/c^{2})"};
      const AxisSpec massAxisOmega = {200, 1.572f, 1.772f, "Inv. Mass (GeV/c^{2})"};

      registry.add("hCandidateCounter", "hCandidateCounter", {HistType::kTH1D, {{10, 0.0f, 10.0f}}});

      // Cascade mass spectra
      registry.add("hMassXiMinus", "hMassXiMinus", {HistType::kTH1D, {{400, 1.122f, 1.522f, "Inv. Mass (GeV/c^{2})"}}});
      registry.add("hMassXiPlus", "hMassXiPlus", {HistType::kTH1D, {{400, 1.122f, 1.522f, "Inv. Mass (GeV/c^{2}²)"}}});
      registry.add("hMassOmegaMinus", "hMassOmegaMinus", {HistType::kTH1D, {{400, 1.472f, 1.872f, "Inv. Mass (GeV/c^{2})"}}});
      registry.add("hMassOmegaPlus", "hMassOmegaPlus", {HistType::kTH1D, {{400, 1.472f, 1.872f, "Inv. Mass (GeV/c^{2})"}}});

      // Cascade topology
      registry.add("hV0Radius", "hV0Radius", {HistType::kTH1D, {{500, 0.0, 100.0, "cm"}}});
      registry.add("hCascRadius", "hCascRadius", {HistType::kTH1D, {{500, 0.0, 100.0, "cm"}}});
      registry.add("hV0CosPA", "hV0CosPA", {HistType::kTH1D, {{100, 0.9f, 1.0f}}});
      registry.add("hCascCosPA", "hCascCosPA", {HistType::kTH1D, {{100, 0.9f, 1.0f}}});
      registry.add("hDCAPosToPV", "hDCAPosToPV", {HistType::kTH1D, {{1000, -10.0f, 10.0f, "cm"}}});
      registry.add("hDCANegToPV", "hDCANegToPV", {HistType::kTH1D, {{1000, -10.0f, 10.0f, "cm"}}});
      registry.add("hDCABachToPV", "hDCABachToPV", {HistType::kTH1D, {{1000, -10.0f, 10.0f, "cm"}}});
      registry.add("hDCAV0ToPV", "hDCAV0ToPV", {HistType::kTH1D, {{1000, -10.0f, 10.0f, "cm"}}});
      registry.add("hDCAV0Dau", "hDCAV0Dau", {HistType::kTH1D, {{500, 0.0f, 5.0f, "cm^{2}"}}});
      registry.add("hDCACascDau", "hDCACascDau", {HistType::kTH1D, {{500, 0.0f, 5.0f, "cm^{2}"}}});
      registry.add("hLambdaMass", "hLambdaMass", {HistType::kTH1D, {{400, 0.916f, 1.316f, "Inv. Mass (GeV/c^{2})"}}});

      // pT rej
      registry.add("hPtCutsXicZeroOmegacZeroToXiPi", "Omegac/Xic to Xi Pi tracks selected by pT;#it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1D, {{500, 0., 50.}}});
      registry.add("hPtCutsOmegacZeroToOmegaPi", "Omegac to Omega Pi tracks selected by pT;#it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1D, {{500, 0., 50.}}});
      registry.add("hPtCutsOmegacZeroToOmegaKa", "Omegac to Omega Ka tracks selected by pT;#it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1D, {{500, 0., 50.}}});
      registry.add("hPtCutsXicPlusToXiPiPi", "Xicplus to Xi Pi Pi tracks selected by pT;#it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1D, {{500, 0., 50.}}});
      registry.add("hRejpTStatusXicZeroOmegacZeroToXiPi", "XicZeroOmegacZeroToXiPi rejected by pT status;status;entries", {HistType::kTH1D, {{2, -0.5, 1.5}}}); // pass dcafitter --> 0, pT>pTmin --> 1
      registry.add("hRejpTStatusOmegacZeroToOmegaPi", "OmegacZeroToOmegaPi rejected by pT status;status;entries", {HistType::kTH1D, {{2, -0.5, 1.5}}});
      registry.add("hRejpTStatusOmegacZeroToOmegaKa", "OmegacZeroToOmegaKa rejected by pT status;status;entries", {HistType::kTH1D, {{2, -0.5, 1.5}}});
      registry.add("hRejpTStatusXicPlusToXiPiPi", "XicPlusToXiPiPi rejected by pT status;status;entries", {HistType::kTH1D, {{2, -0.5, 1.5}}});

      // mass spectra
      registry.add("hMassXicZeroOmegacZeroToXiPi", "2-prong candidates;inv. mass (#Xi #pi) (GeV/#it{c}^{2});entries", {HistType::kTH1D, {{500, 2., 3.}}});
      registry.add("hMassOmegacZeroToOmegaPi", "2-prong candidates;inv. mass (#Omega #pi) (GeV/#it{c}^{2});entries", {HistType::kTH1D, {{500, 2., 3.}}});
      registry.add("hMassOmegacZeroToOmegaK", "2-prong candidates;inv. mass (#Omega K) (GeV/#it{c}^{2});entries", {HistType::kTH1D, {{500, 2., 3.}}});
      registry.add("hMassXicPlusToXiPiPi", "3-prong candidates;inv. mass (#Xi #pi #pi) (GeV/#it{c}^{2});entries", {HistType::kTH1D, {{500, 2., 3.}}});

      // dcaFitter exception counter
      registry.add("hFitterStatusXi2Prong", "Charm DCAFitter status (xi hyp. - 2prong);status;entries", {HistType::kTH1D, {{2, -0.5, 1.5}}});       // 0 --> successful call of DCAFitter 1 --> exception found by DCAFitter
      registry.add("hFitterStatusOmega2Prong", "Charm DCAFitter status (omega hyp. - 2prong);status;entries", {HistType::kTH1D, {{2, -0.5, 1.5}}}); // 0 --> successful call of DCAFitter 1 --> exception found by DCAFitter
      registry.add("hFitterStatusXi3Prong", "Charm DCAFitter status (xi hyp. - 3prong);status;entries", {HistType::kTH1D, {{2, -0.5, 1.5}}});       // 0 --> successful call of DCAFitter 1 --> exception found by DCAFitter
    }
  }

  /// Single-cascade cuts
  template <typename TCascade>
  bool isPreselectedCascade(const TCascade& casc, const float pvx, const float pvy, const float pvz)
  {
    registry.fill(HIST("hCandidateCounter"), 2.5);

    if (casc.v0cosPA(pvx, pvy, pvz) > config.v0CosPA &&
        casc.casccosPA(pvx, pvy, pvz) > config.cascCosPA &&
        casc.dcacascdaughters() < config.dcaCascDau &&
        casc.dcaV0daughters() < config.dcaV0Dau &&
        std::abs(casc.dcanegtopv()) > config.dcaNegToPv &&
        std::abs(casc.dcapostopv()) > config.dcaPosToPv &&
        std::abs(casc.dcabachtopv()) > config.dcaBachToPv &&
        std::abs(casc.dcav0topv(pvx, pvy, pvz)) > config.dcaV0ToPv &&
        casc.v0radius() > config.v0TransvRadius &&
        casc.cascradius() > config.cascTransvRadius &&
        std::abs(casc.mLambda() - MassLambda0) < config.v0MassWindow) {

      registry.fill(HIST("hCandidateCounter"), 3.5); // pass cascade selections

      if (config.fillHistograms) {
        // The basic eleven!
        registry.fill(HIST("hV0Radius"), casc.v0radius());
        registry.fill(HIST("hCascRadius"), casc.cascradius());
        registry.fill(HIST("hV0CosPA"), casc.v0cosPA(pvx, pvy, pvz));
        registry.fill(HIST("hCascCosPA"), casc.casccosPA(pvx, pvy, pvz));
        registry.fill(HIST("hDCAPosToPV"), casc.dcapostopv());
        registry.fill(HIST("hDCANegToPV"), casc.dcanegtopv());
        registry.fill(HIST("hDCABachToPV"), casc.dcabachtopv());
        registry.fill(HIST("hDCAV0ToPV"), casc.dcav0topv(pvx, pvy, pvz));
        registry.fill(HIST("hDCAV0Dau"), casc.dcaV0daughters());
        registry.fill(HIST("hDCACascDau"), casc.dcacascdaughters());
        registry.fill(HIST("hLambdaMass"), casc.mLambda());
        if (casc.sign() < 0) {
          registry.fill(HIST("hMassXiMinus"), casc.mXi());
          registry.fill(HIST("hMassOmegaMinus"), casc.mOmega());
        } else {
          registry.fill(HIST("hMassXiPlus"), casc.mXi());
          registry.fill(HIST("hMassOmegaPlus"), casc.mOmega());
        }
      }
      return true;
    }
    return false;
  }

  /// Method to perform selections for Xic 3-prong candidates before vertex reconstruction
  /// \param pVecXi is the momentum array of the Xi daughter track
  /// \param pVecPi1 is the momentum array of the first pion daughter track
  /// \param pVecPi2 is the momentum array of the second pion daughter track
  /// \return selection outcome
  template <typename T1>
  bool isPreselectedCandidateXic(T1 const& pVecXi, T1 const& pVecPi1, T1 const& pVecPi2)
  {
    // pT
    if (config.ptMinXicplusLfCasc > 0.f) {
      const auto pt = RecoDecay::pt(pVecXi, pVecPi1, pVecPi2) + config.ptTolerance; // add tolerance because of no reco decay vertex
      if (pt < config.ptMinXicplusLfCasc) {
        return false;
      }
    }

    // invariant mass
    if (config.massXiPiPiMin >= 0.f && config.massXiPiPiMax > 0.f) {
      const double invMassMin = config.massXiPiPiMin;
      const double invMassMax = config.massXiPiPiMax;
      const std::array arrMom{pVecXi, pVecPi1, pVecPi2};
      const auto invMass2 = RecoDecay::m2(arrMom, arrMass3Prong[hf_cand_casc_lf::DecayType3Prong::XicplusToXiPiPi]);
      if (invMass2 < invMassMin * invMassMin || invMass2 >= invMassMax * invMassMax) {
        return false;
      }
    }

    return true;
  }

  /// Method to perform selections for Xic 3-prong candidates after vertex reconstruction
  /// \param pVecCand is the momentum array of the candidate after the reconstruction of the secondary vertex
  /// \param secVtx is the secondary vertex
  /// \param primVtx is the primary vertex
  /// \return selection outcome
  template <typename T1, typename T2, typename T3>
  bool isSelectedCandidateXic(const T1& pVecCand, const T2& secVtx, const T3& primVtx)
  {
    // pT
    if (config.ptMinXicplusLfCasc > 0.f) {
      const auto pt = RecoDecay::pt(pVecCand);
      if (pt < config.ptMinXicplusLfCasc) {
        return false;
      }
    }

    // CPA
    if (config.xicCosPA > -1.f) {
      const auto cpa = RecoDecay::cpa(primVtx, secVtx, pVecCand);
      if (cpa < config.xicCosPA) {
        return false;
      }
    }

    // decay length
    if (config.decayLengthXicMin > 0.f) {
      const auto decayLength = RecoDecay::distance(primVtx, secVtx);
      if (decayLength < config.decayLengthXicMin) {
        return false;
      }
    }

    return true;
  }

  void processNoLfCascades(SelectedCollisions const&)
  {
    // dummy
  }
  PROCESS_SWITCH(HfTrackIndexSkimCreatorLfCascades, processNoLfCascades, "Do not skim LF cascades", true);

  void processLfCascades(SelectedCollisions const& collisions,
                         CascFull const& cascades,
                         SelectedHfTrackAssoc const& trackIndices,
                         aod::TracksWCovDca const&,
                         aod::BCsWithTimestamps const&,
                         V0Full const&)
  {
    uint8_t hfFlag = 0;
    bool isGoogForXi2Prong = true;
    bool isGoogForOmega2Prong = true;

    for (const auto& collision : collisions) {

      // set the magnetic field from CCDB
      const auto bc = collision.bc_as<o2::aod::BCsWithTimestamps>();
      initCCDB(bc, runNumber, ccdb, config.isRun2 ? config.ccdbPathGrp : config.ccdbPathGrpMag, lut, config.isRun2);
      const auto magneticField = o2::base::Propagator::Instance()->getNominalBz(); // z component

      df2.setBz(magneticField);
      df2.setRefitWithMatCorr(config.refitWithMatCorr);

      // cascade loop
      const auto thisCollId = collision.globalIndex();
      const auto groupedCascades = cascades.sliceBy(cascadesPerCollision, thisCollId);

      for (const auto& casc : groupedCascades) {

        registry.fill(HIST("hCandidateCounter"), 0.5); // all cascade candidates

        //----------------accessing particles in the decay chain-------------
        // cascade daughter - charged particle
        const auto trackCascDauCharged = casc.bachelor_as<aod::TracksWCovDca>(); // meson <- xi track
        // cascade daughter - V0
        const auto trackV0PosDau = casc.posTrack_as<aod::TracksWCovDca>(); // p <- V0 track (positive track) 0
        // V0 negative daughter
        const auto trackV0NegDau = casc.negTrack_as<aod::TracksWCovDca>(); // pion <- V0 track (negative track) 1

        // check that particles come from the same collision
        if (config.rejDiffCollTrack) {
          if (trackV0PosDau.collisionId() != trackV0NegDau.collisionId()) {
            continue;
          }
          if (trackCascDauCharged.collisionId() != trackV0PosDau.collisionId()) {
            continue;
          }
        }

        if (trackV0PosDau.globalIndex() == trackV0NegDau.globalIndex() || trackV0PosDau.globalIndex() == trackCascDauCharged.globalIndex() || trackV0NegDau.globalIndex() == trackCascDauCharged.globalIndex()) {
          continue;
        }

        if (!(isPreselectedCascade(casc, collision.posX(), collision.posY(), collision.posZ()))) {
          continue;
        }

        const std::array vertexCasc{casc.x(), casc.y(), casc.z()};
        const std::array pVecCasc{casc.px(), casc.py(), casc.pz()};
        std::array<float, 21> covCasc{};
        constexpr std::size_t NIndicesMom{6u};
        constexpr std::size_t MomInd[NIndicesMom] = {9, 13, 14, 18, 19, 20}; // cov matrix elements for momentum component
        for (std::size_t i = 0; i < NIndicesMom; i++) {
          covCasc[MomInd[i]] = casc.momentumCovMat()[i];
          covCasc[i] = casc.positionCovMat()[i];
        }
        // create cascade track
        o2::track::TrackParCov trackParCovCascXi;
        if (trackCascDauCharged.sign() > 0) {
          trackParCovCascXi = o2::track::TrackParCov(vertexCasc, pVecCasc, covCasc, 1, true);
        } else if (trackCascDauCharged.sign() < 0) {
          trackParCovCascXi = o2::track::TrackParCov(vertexCasc, pVecCasc, covCasc, -1, true);
        } else {
          continue;
        }
        trackParCovCascXi.setAbsCharge(1);

        auto trackParCovCascOmega = trackParCovCascXi;

        trackParCovCascXi.setPID(o2::track::PID::XiMinus);
        trackParCovCascOmega.setPID(o2::track::PID::OmegaMinus);

        //--------------combining cascade and pion tracks--------------
        const auto groupedBachTrackIndices = trackIndices.sliceBy(trackIndicesPerCollision, thisCollId);
        for (auto trackIdCharmBachelor1 = groupedBachTrackIndices.begin(); trackIdCharmBachelor1 != groupedBachTrackIndices.end(); ++trackIdCharmBachelor1) {

          hfFlag = 0;
          isGoogForXi2Prong = true;
          isGoogForOmega2Prong = true;

          const auto trackCharmBachelor1 = trackIdCharmBachelor1.track_as<aod::TracksWCovDca>();

          if ((config.rejDiffCollTrack) && (trackCascDauCharged.collisionId() != trackCharmBachelor1.collisionId())) {
            continue;
          }

          // ask for opposite sign daughters
          if (trackCharmBachelor1.sign() * trackCascDauCharged.sign() >= 0) {
            continue;
          }

          // check not to take the same particle twice in the decay chain
          if (trackCharmBachelor1.globalIndex() == trackCascDauCharged.globalIndex() || trackCharmBachelor1.globalIndex() == trackV0PosDau.globalIndex() || trackCharmBachelor1.globalIndex() == trackV0NegDau.globalIndex()) {
            continue;
          }

          // primary pion track to be processed with DCAFitter
          const auto trackParCovCharmBachelor1 = getTrackParCov(trackCharmBachelor1);

          // find charm baryon decay using xi PID hypothesis (xi pi channel)
          int nVtxFrom2ProngFitterXiHyp = 0;
          try {
            nVtxFrom2ProngFitterXiHyp = df2.process(trackParCovCascXi, trackParCovCharmBachelor1);
          } catch (...) {
            if (config.fillHistograms) {
              registry.fill(HIST("hFitterStatusXi2Prong"), 1);
            }
            isGoogForXi2Prong = false;
          }
          if (isGoogForXi2Prong && config.fillHistograms) {
            registry.fill(HIST("hFitterStatusXi2Prong"), 0);
          }

          if (nVtxFrom2ProngFitterXiHyp > 0) {

            df2.propagateTracksToVertex();

            if (df2.isPropagateTracksToVertexDone()) {
              std::array<float, 3> pVecXi{};
              std::array<float, 3> pVecPion1XiHyp{};
              df2.getTrack(0).getPxPyPzGlo(pVecXi);
              df2.getTrack(1).getPxPyPzGlo(pVecPion1XiHyp);
              const float ptXic = RecoDecay::pt(pVecXi, pVecPion1XiHyp);

              const std::array arrMomToXi{pVecXi, pVecPion1XiHyp};
              const auto mass2ProngXiHyp = RecoDecay::m(arrMomToXi, arrMass2Prong[hf_cand_casc_lf::DecayType2Prong::XiczeroOmegaczeroToXiPi]);

              if ((std::abs(casc.mXi() - MassXiMinus) < config.cascadeMassWindow) && (mass2ProngXiHyp >= config.massXiPiMin) && (mass2ProngXiHyp <= config.massXiPiMax)) {
                registry.fill(HIST("hRejpTStatusXicZeroOmegacZeroToXiPi"), 0);
                if (ptXic >= config.ptMinXicZeroOmegacZeroToXiPiLfCasc) {
                  SETBIT(hfFlag, aod::hf_cand_casc_lf::DecayType2Prong::XiczeroOmegaczeroToXiPi);
                  registry.fill(HIST("hRejpTStatusXicZeroOmegacZeroToXiPi"), 1);
                }
              }

              // fill histograms
              if (config.fillHistograms && (TESTBIT(hfFlag, aod::hf_cand_casc_lf::DecayType2Prong::XiczeroOmegaczeroToXiPi))) {
                registry.fill(HIST("hMassXicZeroOmegacZeroToXiPi"), mass2ProngXiHyp);
                registry.fill(HIST("hPtCutsXicZeroOmegacZeroToXiPi"), ptXic);
              }
            } else if (df2.isPropagationFailure()) {
              LOGF(info, "Exception caught: failed to propagate tracks (2prong - xi) to charm baryon decay vtx");
            }
          }

          // find charm baryon decay using omega PID hypothesis to be combined with the charm bachelor (either pion or kaon)
          int nVtxFrom2ProngFitterOmegaHyp = 0;
          try {
            nVtxFrom2ProngFitterOmegaHyp = df2.process(trackParCovCascOmega, trackParCovCharmBachelor1);
          } catch (...) {
            if (config.fillHistograms) {
              registry.fill(HIST("hFitterStatusOmega2Prong"), 1);
            }
            isGoogForOmega2Prong = false;
          }
          if (isGoogForOmega2Prong && config.fillHistograms) {
            registry.fill(HIST("hFitterStatusOmega2Prong"), 0);
          }

          if (nVtxFrom2ProngFitterOmegaHyp > 0) {

            df2.propagateTracksToVertex();

            if (df2.isPropagateTracksToVertexDone()) {

              std::array<float, 3> pVecOmega{};
              std::array<float, 3> pVecCharmBachelor1OmegaHyp{};
              df2.getTrack(0).getPxPyPzGlo(pVecOmega);
              df2.getTrack(1).getPxPyPzGlo(pVecCharmBachelor1OmegaHyp);
              const float ptOmegac = RecoDecay::pt(pVecOmega, pVecCharmBachelor1OmegaHyp);

              const std::array arrMomToOmega{pVecOmega, pVecCharmBachelor1OmegaHyp};
              const auto mass2ProngOmegaPiHyp = RecoDecay::m(arrMomToOmega, arrMass2Prong[hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaPi]);
              const auto mass2ProngOmegaKHyp = RecoDecay::m(arrMomToOmega, arrMass2Prong[hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaK]);

              if (std::abs(casc.mOmega() - MassOmegaMinus) < config.cascadeMassWindow) {
                if ((mass2ProngOmegaPiHyp >= config.massOmegaCharmBachelorMin) && (mass2ProngOmegaPiHyp <= config.massOmegaCharmBachelorMax)) {
                  registry.fill(HIST("hRejpTStatusOmegacZeroToOmegaPi"), 0);
                  if (ptOmegac >= config.ptMinOmegacZeroToOmegaPiLfCasc) {
                    SETBIT(hfFlag, aod::hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaPi);
                    registry.fill(HIST("hRejpTStatusOmegacZeroToOmegaPi"), 1);
                  }
                }
                if ((mass2ProngOmegaKHyp >= config.massOmegaCharmBachelorMin) && (mass2ProngOmegaKHyp <= config.massOmegaCharmBachelorMax)) {
                  registry.fill(HIST("hRejpTStatusOmegacZeroToOmegaKa"), 0);
                  if (ptOmegac >= config.ptMinOmegaczeroToOmegaKaLfCasc) {
                    SETBIT(hfFlag, aod::hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaK);
                    registry.fill(HIST("hRejpTStatusOmegacZeroToOmegaKa"), 1);
                  }
                }
              }

              // fill histograms
              if (config.fillHistograms) {
                if (TESTBIT(hfFlag, aod::hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaPi)) {
                  registry.fill(HIST("hMassOmegacZeroToOmegaPi"), mass2ProngOmegaPiHyp);
                  registry.fill(HIST("hPtCutsOmegacZeroToOmegaPi"), ptOmegac);
                }
                if (TESTBIT(hfFlag, aod::hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaK)) {
                  registry.fill(HIST("hMassOmegacZeroToOmegaK"), mass2ProngOmegaKHyp);
                  registry.fill(HIST("hPtCutsOmegacZeroToOmegaKa"), ptOmegac);
                }
              }
            } else if (df2.isPropagationFailure()) {
              LOGF(info, "Exception caught: failed to propagate tracks (2prong - omega) to charm baryon decay vtx");
            }
          }

          // fill table row only if a vertex was found
          if (hfFlag != 0) {
            rowTrackIndexCasc2Prong(thisCollId,
                                    casc.cascadeId(),
                                    trackCharmBachelor1.globalIndex(),
                                    hfFlag);
          }

          // Xic -> Xi pi pi
          if (config.do3Prong) {
            // Xi mass cut
            if (std::abs(casc.mXi() - MassXiMinus) > config.cascadeMassWindow) {
              continue;
            }

            // second loop over tracks
            for (auto trackIdCharmBachelor2 = trackIdCharmBachelor1 + 1; trackIdCharmBachelor2 != groupedBachTrackIndices.end(); ++trackIdCharmBachelor2) {

              if (!TESTBIT(trackIdCharmBachelor2.isSelProng(), CandidateType::CandCascadeBachelor)) {
                continue;
              }

              const auto trackCharmBachelor2 = trackIdCharmBachelor2.track_as<aod::TracksWCovDca>();

              if ((config.rejDiffCollTrack) && (trackCascDauCharged.collisionId() != trackCharmBachelor2.collisionId())) {
                continue;
              }

              // ask for same sign daughters
              if (trackCharmBachelor2.sign() * trackCharmBachelor1.sign() <= 0) {
                continue;
              }

              // check not to take the same particle twice in the decay chain
              if (trackCharmBachelor2.globalIndex() == trackCharmBachelor1.globalIndex() || trackCharmBachelor2.globalIndex() == trackCascDauCharged.globalIndex() || trackCharmBachelor2.globalIndex() == trackV0PosDau.globalIndex() || trackCharmBachelor2.globalIndex() == trackV0NegDau.globalIndex()) {
                continue;
              }

              if (!isPreselectedCandidateXic(pVecCasc, trackCharmBachelor1.pVector(), trackCharmBachelor2.pVector())) {
                continue;
              }

              // reconstruct Xic with DCAFitter
              // Use only bachelor tracks for vertex reconstruction because the Xi track has large uncertainties.
              int nVtxFrom3ProngFitterXiHyp = 0;
              try {
                const auto trackParCovCharmBachelor2 = getTrackParCov(trackCharmBachelor2);
                nVtxFrom3ProngFitterXiHyp = df2.process(trackParCovCharmBachelor1, trackParCovCharmBachelor2);
              } catch (...) {
                if (config.fillHistograms) {
                  registry.fill(HIST("hFitterStatusXi3Prong"), 1);
                }
                continue;
              }
              if (config.fillHistograms) {
                registry.fill(HIST("hFitterStatusXi3Prong"), 0);
              }

              if (nVtxFrom3ProngFitterXiHyp > 0) {
                df2.propagateTracksToVertex();
                if (df2.isPropagateTracksToVertexDone()) {
                  std::array<float, 3> pVecPi1{};
                  std::array<float, 3> pVecPi2{};
                  // get bachelor momenta at the Xic vertex
                  df2.getTrack(0).getPxPyPzGlo(pVecPi1);
                  df2.getTrack(1).getPxPyPzGlo(pVecPi2);
                  const auto pVecCand = RecoDecay::pVec(pVecCasc, pVecPi1, pVecPi2);
                  const auto ptCand = RecoDecay::pt(pVecCand);
                  const std::array primaryVertex{collision.posX(), collision.posY(), collision.posZ()}; // primary vertex
                  const auto& secondaryVertex = df2.getPCACandidate();                                  // secondary vertex

                  registry.fill(HIST("hRejpTStatusXicPlusToXiPiPi"), 0);
                  if (ptCand >= config.ptMinXicplusLfCasc) {
                    registry.fill(HIST("hRejpTStatusXicPlusToXiPiPi"), 1);
                  }

                  if (!isSelectedCandidateXic(pVecCand, secondaryVertex, primaryVertex)) {
                    continue;
                  }

                  // fill histograms
                  if (config.fillHistograms) {
                    const std::array arr3Mom{pVecCasc, pVecPi1, pVecPi2};
                    const auto mass3Prong = RecoDecay::m(arr3Mom, arrMass3Prong[hf_cand_casc_lf::DecayType3Prong::XicplusToXiPiPi]);
                    registry.fill(HIST("hMassXicPlusToXiPiPi"), mass3Prong);
                    registry.fill(HIST("hPtCutsXicPlusToXiPiPi"), ptCand);
                  }

                  // fill table row if a vertex was found
                  rowTrackIndexCasc3Prong(thisCollId,
                                          casc.cascadeId(),
                                          trackCharmBachelor1.globalIndex(),
                                          trackCharmBachelor2.globalIndex());
                } else if (df2.isPropagationFailure()) {
                  LOGF(info, "Exception caught: failed to propagate tracks (3prong) to charm baryon decay vtx");
                }
              }
            } // end 3prong loop
          } // end 3prong condition

        } // loop over pion
      } // loop over cascade
    } // loop over collisions
  } // processLfCascades

  PROCESS_SWITCH(HfTrackIndexSkimCreatorLfCascades, processLfCascades, "Skim HF -> LF cascade + bachelor", false);
};

//________________________________________________________________________________________________________________________
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{};
  workflow.push_back(adaptAnalysisTask<HfTrackIndexSkimCreatorTagSelCollisions>(cfgc));
  workflow.push_back(adaptAnalysisTask<HfTrackIndexSkimCreatorTagSelTracks>(cfgc));
  workflow.push_back(adaptAnalysisTask<HfTrackIndexSkimCreator>(cfgc));
  workflow.push_back(adaptAnalysisTask<HfTrackIndexSkimCreatorCascades>(cfgc));
  workflow.push_back(adaptAnalysisTask<HfTrackIndexSkimCreatorLfCascades>(cfgc));
  return workflow;
}
