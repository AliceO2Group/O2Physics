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

/// \file dataCreatorHiddenCharmReduced.cxx
/// \brief Reduced data creator for hidden-charm analyses at midrapidity
///
/// \author A. Palasciano, <antonio.palasciano@cern.ch>, INFN Bari
/// \author S. Politanò <stefano.politano@cern.ch>, CERN

#include "PWGHF/Core/CentralityEstimation.h"
#include "PWGHF/Core/SelectorCuts.h"
#include "PWGHF/D2H/DataModel/ReducedDataModel.h"
#include "PWGHF/D2H/Utils/utilsRedDataFormat.h"
#include "PWGHF/Utils/utilsAnalysis.h"
#include "PWGHF/Utils/utilsEvSelHf.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/TPCVDriftManager.h"
#include "Common/Core/ZorroSummary.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CCDB/CcdbApi.h>
#include <CommonConstants/PhysicsConstants.h>
#include <DetectorsBase/MatLayerCylSet.h>
#include <DetectorsBase/Propagator.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Array2D.h>
#include <Framework/Configurable.h>
#include <Framework/DeviceSpec.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/RunningWorkflowInfo.h>
#include <Framework/runDataProcessing.h>

#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::analysis;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::hf_centrality;

enum TrackType : uint8_t {
  Pion = 0,
  Kaon,
  Proton
};

struct HfDataCreatorHiddenCharmReduced {

  Produces<aod::HfRedCollisions> hfReducedCollision;
  Produces<aod::HfOrigColCounts> hfCollisionCounter;
  Produces<aod::HcSelTracks> hfTrackLite;

  struct : ConfigurableGroup {
    // track quality
    Configurable<bool> fillHistograms{"fillHistograms", true, "Fill proton QA histograms"};
    Configurable<bool> selectProtons{"selectProtons", true, "Select protons"};
    Configurable<int> itsNClsMin{"itsNClsMin", 5, "Minimum number of ITS clusters"};
    Configurable<int> tpcNClsFoundMin{"tpcNClsFoundMin", 50, "Minimum number of found TPC clusters"};
    Configurable<int> tpcNClsCrossedRowsMin{"tpcNClsCrossedRowsMin", 80, "Minimum number of crossed TPC rows"};
    Configurable<double> ptMinTrack{"ptMinTrack", 0.5, "Minimum proton-track pT"};
    Configurable<double> etaMaxTrack{"etaMaxTrack", 0.8, "Maximum proton-track |eta|"};
    Configurable<float> trackChi2Cut{"trackChi2Cut", 4.f, "Maximum chi2/ncls in TPC"};
    Configurable<float> trackMinChi2Cut{"trackMinChi2Cut", 0.f, "Minimum chi2/ncls in TPC"};
    Configurable<float> trackMaxChi2ITS{"trackMaxChi2ITS", 36.f, "Maximum chi2/ncls in ITS"};
    Configurable<std::vector<double>> binsPtTrack{"binsPtTrack", std::vector<double>{hf_cuts_single_track::vecBinsPtTrack}, "Track pT bin limits for DCA cuts"};
    Configurable<LabeledArray<double>> cutsTrack{"cutsTrack", {hf_cuts_single_track::CutsTrack[0], hf_cuts_single_track::NBinsPtTrack, hf_cuts_single_track::NCutVarsTrack, hf_cuts_single_track::labelsPtTrack, hf_cuts_single_track::labelsCutVarTrack}, "Single-track DCA selections per pT bin"};
    // DCA
    Configurable<std::vector<float>> paramsDCAxyPtDep{"paramsDCAxyPtDep", std::vector<float>{0.0010, 0.0080, 0.73}, "Parameters for pT-dependent DCAxy cut: p0, p1, p2 for cut = p0 + p1/pt^p2"};
    Configurable<std::vector<float>> paramsDCAzPtDep{"paramsDCAzPtDep", std::vector<float>{-0.0044, 0.0152, 0.47}, "Parameters for pT-dependent DCAz cut: p0, p1, p2 for cut = p0 + p1/pt^p2"};
    // PID
    Configurable<float> momForCombinedPid{"momForCombinedPid", 0.75f, "Momentum threshold above which combined TPC+TOF proton PID is used"};
    Configurable<float> maxNsigmaTofPi{"maxNsigmaTofPi", 2.f, "Maximum pion n-sigma in TOF for proton rejection"};
    Configurable<float> maxNsigmaTofKa{"maxNsigmaTofKa", 2.f, "Maximum kaon n-sigma in TOF for proton rejection"};
    Configurable<float> maxNsigmaTofPr{"maxNsigmaTofPr", 3.f, "Maximum proton n-sigma in TOF"};
    Configurable<float> maxNsigmaCombinedPr{"maxNsigmaCombinedPr", 3.f, "Maximum combined proton n-sigma from TPC and TOF"};
    Configurable<float> maxNsigmaTpcPi{"maxNsigmaTpcPi", 2.f, "Maximum pion n-sigma in TPC for proton rejection"};
    Configurable<float> maxNsigmaTpcKa{"maxNsigmaTpcKa", 2.f, "Maximum kaon n-sigma in TPC for proton rejection"};
    Configurable<float> maxNsigmaTpcPr{"maxNsigmaTpcPr", 3.f, "Maximum proton n-sigma in TPC"};
    Configurable<bool> forceTOF{"forceTOF", false, "Require TOF PID information for proton selection"};
  } config;

  // Configurable (others)
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> ccdbPathGrpMag{"ccdbPathGrpMag", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object (Run 3)"};
  Configurable<bool> doMcRecQa{"doMcRecQa", true, "Fill QA histograms for Mc matching"};
  Configurable<bool> rejectPairsWithCommonDaughter{"rejectPairsWithCommonDaughter", true, "flag to reject already at this stage the pairs that share a daughter track"};
  Configurable<bool> rejectCollisionsWithBadEvSel{"rejectCollisionsWithBadEvSel", true, "flag to reject collisions with bad event selection"};

  o2::hf_evsel::HfEventSelection hfEvSel;
  o2::hf_evsel::HfEventSelectionMc hfEvSelMc;

  double bz{0.};
  int runNumber{0}; // needed to detect if the run changed and trigger update of calibrations etc.

  // material correction for track propagation
  o2::base::MatLayerCylSet* lut{};
  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;
  o2::aod::common::TPCVDriftManager vDriftMgr;
  o2::ccdb::CcdbApi ccdbApi;
  Service<o2::ccdb::BasicCCDBManager> ccdb{};

  using BCsInfo = soa::Join<aod::BCs, aod::Timestamps, aod::BcSels>;
  using TracksWithPID = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullKa, aod::pidTOFFullKa, aod::pidTPCFullPr, aod::pidTOFFullPr>;

  Preslice<aod::TrackAssoc> trackIndicesPerCollision = aod::track_association::collisionId;

  HistogramRegistry registry{"registry"};
  OutputObj<ZorroSummary> zorroSummary{"zorroSummary"};

  void init(InitContext& initContext)
  {
    // Configure CCDB access
    ccdb->setURL(ccdbUrl.value);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setCreatedNotAfter(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
    ccdbApi.init(ccdbUrl);
    runNumber = 0;
    lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>("GLO/Param/MatLUT"));

    if (config.fillHistograms) {
      const AxisSpec axisPt{360, 0., 36., "#it{p}_{T}^{proton} (GeV/#it{c})"};
      const AxisSpec axisEta{100, -1., 1., "#eta"};
      const AxisSpec axisDca{400, -2., 2., "DCA_{xy} to primary vertex (cm)"};
      const AxisSpec axisNSigma{100, -5., 5., "n#sigma"};

      registry.add("hPzVtx", "Z position of primary vertex for selected tracks;z_{vtx} (cm);entries", {HistType::kTH1D, {AxisSpec{200, -20., 20., "z_{vtx} (cm)"}}});
      registry.add("hPtNoCuts", "All associated tracks;#it{p}_{T}^{track} (GeV/#it{c});entries", {HistType::kTH1D, {axisPt}});
      registry.add("hPtCutsProton", "Selected proton tracks;#it{p}_{T}^{track} (GeV/#it{c});entries", {HistType::kTH1D, {axisPt}});
      registry.add("hEtaCutsProton", "Selected proton tracks;#eta;entries", {HistType::kTH1D, {axisEta}});
      registry.add("hDCAToPrimXYVsPtCutsProton", "Selected proton tracks;#it{p}_{T}^{track} (GeV/#it{c});DCA_{xy} to primary vertex (cm)", {HistType::kTH2D, {axisPt, axisDca}});
      registry.add("hNSigmaTPCProton", "Selected proton tracks;#it{p}_{T}^{track} (GeV/#it{c});n#sigma_{TPC}", {HistType::kTH2D, {axisPt, axisNSigma}});
      registry.add("hNSigmaTOFProton", "Selected proton tracks;#it{p}_{T}^{track} (GeV/#it{c});n#sigma_{TOF}", {HistType::kTH2D, {axisPt, axisNSigma}});
      registry.add("hInvMass", "Invariant mass of selected proton with all other tracks in the event;#it{p}_{T}^{proton} (GeV/#it{c});invariant mass with other tracks (GeV/#it{c}^{2})", {HistType::kTH2D, {axisPt, AxisSpec{100, 2.85, 3.25, "invariant mass with other tracks (GeV/#it{c}^{2})"}}});
      registry.add("hDeDxTPCProton", "Selected proton tracks;#it{p}_{T}^{track} (GeV/#it{c});TPC dE/dx (a.u.)", {HistType::kTH2D, {axisPt, AxisSpec{100, 0., 200., "TPC dE/dx (a.u.)"}}});
    }

    // init HF event selection helper
    hfEvSel.init(registry, &zorroSummary);

    const auto& workflows = initContext.services().get<RunningWorkflowInfo const>();
    for (const DeviceSpec& device : workflows.devices) {
      if (device.name == "hf-data-creator-charm-reso-to-d0-reduced") {
        // init HF event selection helper
        hfEvSelMc.init(device, registry);
        break;
      }
    }
  }

  static float dcaSigma(float const& pt, float const& p0, float const& p1, float const& p2)
  {
    return p0 + p1 / std::pow(std::abs(pt), p2);
  }

  template <typename TTrack>
  bool isSelectedPid(TTrack const& track) const
  {
    const float momForCombinedPid = config.momForCombinedPid.value;
    const float maxNsigmaTpcPr = config.maxNsigmaTpcPr.value;
    const float maxNsigmaTpcPi = config.maxNsigmaTpcPi.value;
    const float maxNsigmaTpcKa = config.maxNsigmaTpcKa.value;
    const float maxNsigmaCombinedPr = config.maxNsigmaCombinedPr.value;
    const float maxNsigmaTofPi = config.maxNsigmaTofPi.value;
    const float maxNsigmaTofKa = config.maxNsigmaTofKa.value;
    const bool forceTOF = config.forceTOF.value;

    const float mom = std::hypot(track.px(), track.py(), track.pz());
    const float nSigmaTPCPr = track.tpcNSigmaPr();
    const float nSigmaTPCPi = track.tpcNSigmaPi();
    const float nSigmaTPCKa = track.tpcNSigmaKa();
    const bool hasTOF = track.hasTOF();
    if (!hasTOF && forceTOF) {
      return false;
    }

    bool isProton = false;
    bool rejectAsPion = false;
    bool rejectAsKaon = false;

    if (mom < momForCombinedPid) {
      isProton = std::abs(nSigmaTPCPr) < maxNsigmaTpcPr;
      rejectAsPion = std::abs(nSigmaTPCPi) < maxNsigmaTpcPi;
      rejectAsKaon = std::abs(nSigmaTPCKa) < maxNsigmaTpcKa;

      if (hasTOF) {
        rejectAsPion = rejectAsPion || std::abs(track.tofNSigmaPi()) < maxNsigmaTofPi;
        rejectAsKaon = rejectAsKaon || std::abs(track.tofNSigmaKa()) < maxNsigmaTofKa;
        isProton = isProton || std::abs(track.tofNSigmaPr()) < maxNsigmaCombinedPr;
      }

    } else {
      const float nSigmaTOFPr = track.tofNSigmaPr();
      const float nSigmaTOFPi = track.tofNSigmaPi();
      const float nSigmaTOFKa = track.tofNSigmaKa();
      isProton = std::hypot(nSigmaTPCPr, nSigmaTOFPr) < maxNsigmaCombinedPr;
      rejectAsPion = std::hypot(nSigmaTPCPi, nSigmaTOFPi) < maxNsigmaTofPi;
      rejectAsKaon = std::hypot(nSigmaTPCKa, nSigmaTOFKa) < maxNsigmaTofKa;
    }
    return isProton && !rejectAsPion && !rejectAsKaon;
  }

  template <typename TTrack>
  bool isSelectedTrack(TTrack const& track) const
  {
    const int tpcNClsFoundMin = config.tpcNClsFoundMin.value;
    const int tpcNClsCrossedRowsMin = config.tpcNClsCrossedRowsMin.value;
    const int itsNClsMin = config.itsNClsMin.value;
    const double etaMaxTrack = config.etaMaxTrack.value;
    const double ptMinTrack = config.ptMinTrack.value;
    const float trackChi2Cut = config.trackChi2Cut.value;
    const float trackMinChi2Cut = config.trackMinChi2Cut.value;
    const float trackMaxChi2ITS = config.trackMaxChi2ITS.value;
    const float dcaXY = track.dcaXY();
    const float dcaZ = track.dcaZ();

    if (!track.isGlobalTrackWoDCA()) {
      return false;
    }
    if (track.tpcNClsFound() < tpcNClsFoundMin) {
      return false;
    }
    if (track.tpcNClsCrossedRows() < tpcNClsCrossedRowsMin) {
      return false;
    }
    if (track.itsNCls() < itsNClsMin) {
      return false;
    }
    if (track.tpcChi2NCl() > trackChi2Cut || track.tpcChi2NCl() < trackMinChi2Cut) {
      return false;
    }
    if (track.itsChi2NCl() > trackMaxChi2ITS) {
      return false;
    }
    if (std::abs(track.eta()) > etaMaxTrack) {
      return false;
    }
    if (track.pt() < ptMinTrack) {
      return false;
    }
    if (!isSelectedTrackDca(config.binsPtTrack, config.cutsTrack, track.pt(), track.dcaXY(), track.dcaZ())) {
      return false;
    }
    if (dcaSigma(track.pt(), config.paramsDCAxyPtDep.value[0], config.paramsDCAxyPtDep.value[1], config.paramsDCAxyPtDep.value[2]) > std::abs(dcaXY) || dcaSigma(track.pt(), config.paramsDCAzPtDep.value[0], config.paramsDCAzPtDep.value[1], config.paramsDCAzPtDep.value[2]) > std::abs(dcaZ)) {
      return false;
    }
    return isSelectedPid(track);
  }

  void processEtaCTrack(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                        aod::TrackAssoc const& trackIndices,
                        TracksWithPID const& tracks,
                        aod::BCsWithTimestamps const&)
  {
    hfTrackLite.reserve(tracks.size());

    int zvtxColl{0};
    int sel8Coll{0};
    int zvtxAndSel8Coll{0};
    int zvtxAndSel8CollAndSoftTrig{0};
    int allSelColl{0};
    constexpr int NDaughtersCharmMeson = 2;
    for (const auto& collision : collisions) {
      const auto hfRejMap = o2::hf_evsel::getEvSel<true, o2::hf_centrality::CentralityEstimator::None, aod::BCsWithTimestamps>(collision, hfEvSel, zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl, ccdb, registry);
      if (rejectCollisionsWithBadEvSel && hfRejMap != 0) {
        continue;
      }
      const auto thisCollId = collision.globalIndex();
      const auto trackIds = trackIndices.sliceBy(trackIndicesPerCollision, thisCollId);
      if (config.fillHistograms) {
        registry.fill(HIST("hPzVtx"), collision.posZ());
      }
      std::vector<int> selectedTrackIds;
      for (const auto& trkId : trackIds) {
        auto trk = trkId.track_as<TracksWithPID>();
        if (config.fillHistograms) {
          registry.fill(HIST("hPtNoCuts"), trk.pt());
        }
        if (!isSelectedTrack(trk)) {
          continue;
        }
        std::array pVecProton{trk.pVector()};
        hfTrackLite(trk.globalIndex(), collision.globalIndex(), pVecProton[0], pVecProton[1], pVecProton[2], trk.sign(), static_cast<uint8_t>(TrackType::Proton));
        selectedTrackIds.push_back(trk.globalIndex());
        if (config.fillHistograms) {
          registry.fill(HIST("hPtCutsProton"), trk.pt());
          registry.fill(HIST("hEtaCutsProton"), trk.eta());
          registry.fill(HIST("hDCAToPrimXYVsPtCutsProton"), trk.pt(), trk.dcaXY());
          registry.fill(HIST("hNSigmaTPCProton"), trk.pt(), trk.tpcNSigmaPr());
          registry.fill(HIST("hDeDxTPCProton"), trk.pt(), trk.tpcSignal());
          if (trk.hasTOF()) {
            registry.fill(HIST("hNSigmaTOFProton"), trk.pt(), trk.tofNSigmaPr());
          }
        }
      }
      if (selectedTrackIds.size() < NDaughtersCharmMeson) {
        continue;
      }
      hfReducedCollision(collision.posX(), collision.posY(), collision.posZ(), collision.numContrib(), hfRejMap, bz);
      for (size_t i = 0; i < selectedTrackIds.size(); ++i) {
        auto t1 = tracks.rawIteratorAt(selectedTrackIds[i]);
        std::array pVec1{t1.pVector()};
        for (size_t j = i + 1; j < selectedTrackIds.size(); ++j) {
          auto t2 = tracks.rawIteratorAt(selectedTrackIds[j]);
          if (t1.sign() * t2.sign() > 0) {
            continue;
          }
          std::array pVec2{t2.pVector()};
          float invMass = RecoDecay::m(std::array{pVec1, pVec2}, std::array{o2::constants::physics::MassProton, o2::constants::physics::MassProton});
          float ptEtac = RecoDecay::pt(RecoDecay::sumOfVec(pVec1, pVec2));
          if (config.fillHistograms) {
            registry.fill(HIST("hInvMass"), ptEtac, invMass);
          }
        }
      }
    }
    hfCollisionCounter(collisions.tableSize(), zvtxColl, sel8Coll, zvtxAndSel8Coll, zvtxAndSel8CollAndSoftTrig, allSelColl);
  }

  PROCESS_SWITCH(HfDataCreatorHiddenCharmReduced, processEtaCTrack, "EtaC -> p pbar reconstruction", true);
};

// struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{};
  workflow.push_back(adaptAnalysisTask<HfDataCreatorHiddenCharmReduced>(cfgc));
  return workflow;
}
