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
// ========================
//
// This code creates parameters used in track tuner.
//    Please write to: daiki.sekihata@cern.ch

#include "PWGEM/Dilepton/Utils/MCUtilities.h"
#include "PWGEM/Dilepton/Utils/PairUtilities.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/CCDB/RCTSelectionFlags.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/Zorro.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/Qvectors.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <DataFormatsCalibration/MeanVertexObject.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <DataFormatsParameters/GRPObject.h>
#include <DetectorsBase/MatLayerCylSet.h>
#include <DetectorsBase/Propagator.h>
#include <Framework/ASoAHelpers.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/DataTypes.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>
#include <ReconstructionDataFormats/DCA.h>
#include <ReconstructionDataFormats/PID.h>

#include <Math/Vector4D.h> // IWYU pragma: keep (do not replace with Math/Vector4Dfwd.h)
#include <Math/Vector4Dfwd.h>
#include <TH1.h>
#include <TString.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iterator>
#include <optional>
#include <string>
#include <string_view>
#include <type_traits>
#include <vector>

#include <math.h>

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

struct createTTP {
  // Configurables
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
  Configurable<int> cfgEventGeneratorType{"cfgEventGeneratorType", -1, "if positive, select event generator type. i.e. gap or signal"};

  ConfigurableAxis ConfPtBins{"ConfPtBins", {VARIABLE_WIDTH, 0.00, 0.10, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.5, 6.0, 7.5, 7.0, 7.5, 8.0, 9.0, 10.0, 15, 20}, "pT bins for output histograms"};
  Configurable<int> cfgNbinsEta{"cfgNbinsEta", 20, "number of eta bins for output histograms"};
  Configurable<int> cfgNbinsPhi{"cfgNbinsPhi", 36, "number of phi bins for output histograms"};
  ConfigurableAxis ConfDCABins{"ConfDCABins", {200, -1000, +1000}, "DCA bins for output histograms"};
  ConfigurableAxis ConfDCASigmaBins{"ConfDCASigmaBins", {200, -10, +10}, "DCA in sigma bins for output histograms"};

  struct : ConfigurableGroup {
    std::string prefix = "eventCut";
    Configurable<float> cfgZvtxMin{"cfgZvtxMin", -10, "min. Zvtx"};
    Configurable<float> cfgZvtxMax{"cfgZvtxMax", +10, "max. Zvtx"};
    Configurable<bool> cfgRequireFT0AND{"cfgRequireFT0AND", true, "require FT0AND in event cut"};
    Configurable<bool> cfgRequireNoTFB{"cfgRequireNoTFB", true, "require No time frame border in event cut"};
    Configurable<bool> cfgRequireNoITSROFB{"cfgRequireNoITSROFB", false, "require no ITS readout frame border in event cut"};
    Configurable<bool> cfgRequireVertexITSTPC{"cfgRequireVertexITSTPC", false, "require Vertex ITSTPC in event cut"};             // ITS-TPC matched track contributes PV.
    Configurable<bool> cfgRequireVertexTOFmatched{"cfgRequireVertexTOFmatched", false, "require Vertex TOFmatched in event cut"}; // ITS-TPC-TOF matched track contributes PV.
    Configurable<bool> cfgRequireNoSameBunchPileup{"cfgRequireNoSameBunchPileup", false, "require no same bunch pileup in event cut"};
    Configurable<bool> cfgRequireGoodZvtxFT0vsPV{"cfgRequireGoodZvtxFT0vsPV", false, "require good Zvtx between FT0 vs. PV in event cut"};
    Configurable<int> cfgTrackOccupancyMin{"cfgTrackOccupancyMin", -2, "min. track occupancy"};
    Configurable<int> cfgTrackOccupancyMax{"cfgTrackOccupancyMax", 1000000000, "max. track occupancy"};
    Configurable<float> cfgFT0COccupancyMin{"cfgFT0COccupancyMin", -2, "min. FT0C occupancy"};
    Configurable<float> cfgFT0COccupancyMax{"cfgFT0COccupancyMax", 1000000000, "max. FT0C occupancy"};

    Configurable<int> cfgCentEstimator{"cfgCentEstimator", 2, "FT0M:0, FT0A:1, FT0C:2"};
    Configurable<float> cfgCentMin{"cfgCentMin", -1.f, "min. centrality"};
    Configurable<float> cfgCentMax{"cfgCentMax", 999.f, "max. centrality"};
    Configurable<uint16_t> cfgNumContribMin{"cfgNumContribMin", 0, "min. numContrib"};
    Configurable<uint16_t> cfgNumContribMax{"cfgNumContribMax", 65000, "max. numContrib"};

    // for RCT
    Configurable<bool> cfgRequireGoodRCT{"cfgRequireGoodRCT", true, "require good detector flag in run condtion table"};
    Configurable<std::string> cfgRCTLabel{"cfgRCTLabel", "CBT_hadronPID", "select 1 [CBT, CBT_hadronPID, CBT_muon_glo] see O2Physics/Common/CCDB/RCTSelectionFlags.h"};
    Configurable<bool> cfgCheckZDC{"cfgCheckZDC", false, "set ZDC flag for PbPb"};
    Configurable<bool> cfgTreatLimitedAcceptanceAsBad{"cfgTreatLimitedAcceptanceAsBad", false, "reject all events where the detectors relevant for the specified Runlist are flagged as LimitedAcceptance"};
  } eventCut;

  struct : ConfigurableGroup {
    std::string prefix = "tagCut";
    Configurable<float> min_pt_track{"min_pt_track", 0.4, "min pT for single track"};
    Configurable<float> max_pt_track{"max_pt_track", 1e+10, "max pT for single track"};
    Configurable<float> min_eta_track{"min_eta_track", -0.8, "min eta for single track"};
    Configurable<float> max_eta_track{"max_eta_track", +0.8, "max eta for single track"};
    Configurable<int> min_ncluster_tpc{"min_ncluster_tpc", 0, "min ncluster tpc"};
    Configurable<int> min_ncrossedrows_tpc{"min_ncrossedrows_tpc", 120, "min. crossed rows"};
    Configurable<int> min_ncluster_its{"min_ncluster_its", 5, "min ncluster its"};
    Configurable<int> min_ncluster_itsib{"min_ncluster_itsib", 3, "min ncluster itsib"};
    Configurable<float> min_tpc_cr_findable_ratio{"min_tpc_cr_findable_ratio", 0.8, "min. TPC Ncr/Nf ratio"};
    Configurable<float> max_frac_shared_clusters_tpc{"max_frac_shared_clusters_tpc", 0.7, "max fraction of shared clusters in TPC"};
    Configurable<float> max_chi2tpc{"max_chi2tpc", 2.5, "max. chi2/NclsTPC"};
    Configurable<float> max_chi2its{"max_chi2its", 36.0, "max. chi2/NclsITS"};
    Configurable<float> max_dcaxy{"max_dcaxy", 0.2, "max DCAxy in cm"};
    Configurable<float> max_dcaz{"max_dcaz", 0.2, "max DCAz in cm"};
    Configurable<float> max_dca_in_sigma{"max_dca_in_sigma", 1e+10, "max dca in sigma for a single track"};

    Configurable<float> min_TPCNsigmaEl{"min_TPCNsigmaEl", -2, "min. TPC n sigma for electron inclusion"};
    Configurable<float> max_TPCNsigmaEl{"max_TPCNsigmaEl", +3, "max. TPC n sigma for electron inclusion"};
    Configurable<float> min_TOFNsigmaEl{"min_TOFNsigmaEl", -3, "min. TOF n sigma for electron inclusion"};
    Configurable<float> max_TOFNsigmaEl{"max_TOFNsigmaEl", +3, "max. TOF n sigma for electron inclusion"};
  } tagCut;

  struct : ConfigurableGroup {
    std::string prefix = "probeCut";
    Configurable<float> min_pt_track{"min_pt_track", 0.05, "min pT for single track"};
    Configurable<float> max_pt_track{"max_pt_track", 1e+10, "max pT for single track"};
    Configurable<float> min_eta_track{"min_eta_track", -0.8, "min eta for single track"};
    Configurable<float> max_eta_track{"max_eta_track", +0.8, "max eta for single track"};
    Configurable<int> min_ncrossedrows_tpc{"min_ncrossedrows_tpc", 80, "min ncrossed rows"};
    Configurable<int> min_ncluster_tpc{"min_ncluster_tpc", 0, "min ncluster tpc"};
    Configurable<int> min_ncluster_its{"min_ncluster_its", 5, "min ncluster its"};
    Configurable<int> min_ncluster_itsib{"min_ncluster_itsib", 1, "min ncluster itsib"};
    Configurable<float> min_cr2findable_ratio_tpc{"min_cr2findable_ratio_tpc", 0.8, "min. TPC Ncr/Nf ratio"};
    Configurable<float> max_frac_shared_clusters_tpc{"max_frac_shared_clusters_tpc", 0.7f, "max fraction of shared clusters in TPC"};
    Configurable<float> max_chi2tpc{"max_chi2tpc", 4.0, "max chi2/NclsTPC"};
    Configurable<float> max_chi2its{"max_chi2its", 36.0, "max chi2/NclsITS"};
    Configurable<float> max_dcaxy{"max_dcaxy", 1, "max dca XY for single track in cm"};
    Configurable<float> max_dcaz{"max_dcaz", 1, "max dca Z for single track in cm"};

    Configurable<float> min_TPCNsigmaEl{"min_TPCNsigmaEl", -2, "min n sigma e in TPC"};
    Configurable<float> max_TPCNsigmaEl{"max_TPCNsigmaEl", +3, "max n sigma e in TPC"};
    Configurable<float> min_TOFNsigmaEl{"min_TOFNsigmaEl", -3, "min n sigma e in TOF"};
    Configurable<float> max_TOFNsigmaEl{"max_TOFNsigmaEl", +3, "max n sigma e in TOF"};
  } probeCut;

  struct : ConfigurableGroup {
    std::string prefix = "pairCut";
    Configurable<float> minMee{"minMee", 0.00, "min mee for pi0 -> ee"};
    Configurable<float> maxMee{"maxMee", 0.01, "max mee for pi0 -> ee"};
    Configurable<float> minPhiV{"minPhiV", 0.f, "min phiv for pi0 -> ee"};
    Configurable<float> maxPhiV{"maxPhiV", M_PI / 2, "max phiv for pi0 -> ee"};
  } pairCut;

  o2::aod::rctsel::RCTFlagsChecker rctChecker;

  int mRunNumber{0};
  float d_bz{0};
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;
  o2::dataformats::VertexBase mVtx;
  o2::base::MatLayerCylSet* lut = nullptr;

  HistogramRegistry fRegistry{"output", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};

  void init(InitContext&)
  {
    mRunNumber = 0;
    d_bz = 0;
    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);
    rctChecker.init(eventCut.cfgRCTLabel.value, eventCut.cfgCheckZDC.value, eventCut.cfgTreatLimitedAcceptanceAsBad.value);

    addhistograms();
  }

  ~createTTP() {}

  template <typename TBC>
  void initCCDB(TBC const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }

    // load matLUT for this timestamp
    if (!lut) {
      LOG(info) << "Loading material look-up table for timestamp: " << bc.timestamp();
      lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->getForTimeStamp<o2::base::MatLayerCylSet>(lutPath, bc.timestamp()));
    } else {
      LOG(info) << "Material look-up table already in place. Not reloading.";
    }

    auto run3grp_timestamp = bc.timestamp();
    o2::parameters::GRPObject* grpo = 0x0;
    o2::parameters::GRPMagField* grpmag = 0x0;
    if (grpo) {
      o2::base::Propagator::initFieldFromGRP(grpo);
      o2::base::Propagator::Instance()->setMatLUT(lut);
      // Fetch magnetic field from ccdb for current collision
      d_bz = grpo->getNominalL3Field();
      LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << d_bz << " kZG";
    } else {
      grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, run3grp_timestamp);
      if (!grpmag) {
        LOG(fatal) << "Got nullptr from CCDB for path " << grpmagPath << " of object GRPMagField and " << grpPath << " of object GRPObject for timestamp " << run3grp_timestamp;
      }
      o2::base::Propagator::initFieldFromGRP(grpmag);
      o2::base::Propagator::Instance()->setMatLUT(lut);
      // Fetch magnetic field from ccdb for current collision
      d_bz = std::lround(5.f * grpmag->getL3Current() / 30000.f);
      LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << d_bz << " kZG";
    }

    mRunNumber = bc.runNumber();
  }

  void addhistograms()
  {
    // event info
    auto hCollisionCounter = fRegistry.add<TH1>("Event/hCollisionCounter", "collision counter;;Number of events", kTH1D, {{2, -0.5, 1.5}}, false);
    hCollisionCounter->GetXaxis()->SetBinLabel(1, "all");
    hCollisionCounter->GetXaxis()->SetBinLabel(2, "accepted");

    const AxisSpec axis_pt{ConfPtBins, "p_{T} (GeV/c)"};
    const AxisSpec axis_eta{cfgNbinsEta, -1, +1, "#eta"};
    const AxisSpec axis_phi{cfgNbinsPhi, 0.0, 2 * M_PI, "#varphi (rad.)"};
    const AxisSpec axis_sign{3, -1.5, +1.5, "sign"};
    const AxisSpec axis_mass{400, 0, 4, "m_{ee} (GeV/c^{2})"};
    const AxisSpec axis_ptee{100, 0, 10, "p_{T,ee} (GeV/c)"};

    const AxisSpec axis_mean_dcaXY{ConfDCABins, "DCA_{xy} (#mum)"};
    const AxisSpec axis_mean_dcaZ{ConfDCABins, "DCA_{z} (#mum)"};
    const AxisSpec axis_pull_dcaXY{ConfDCASigmaBins, "DCA_{xy}/#sigma_{xy}^{DCA}"};
    const AxisSpec axis_pull_dcaZ{ConfDCASigmaBins, "DCA_{z}/#sigma_{z}^{DCA}"};

    fRegistry.add("Track/hs", "electron", kTHnSparseD, {axis_pt, axis_eta, axis_phi, axis_sign, axis_mean_dcaXY, axis_mean_dcaZ, axis_pull_dcaXY, axis_pull_dcaZ}, false);
    fRegistry.add("Track/hTPCNsigmaEl", "TPC n sigma el;p_{in} (GeV/c);n #sigma_{e}^{TPC}", kTH2F, {{1000, 0, 10}, {100, -5.f, +5.f}}, false);
    fRegistry.add("Pair/hMvsPhiV", "m_{ee} vs. #varphi_{V} ULS;#varphi_{V} (rad.);m_{ee} (GeV/c^{2})", kTH2F, {{90, 0.f, M_PI}, {100, 0, 0.1}});

    if (doprocessMC) {
      fRegistry.add("Pair/hMvsPt_omega", "#omega->ee", kTH2D, {axis_mass, axis_ptee}, false);
      fRegistry.add("Pair/hMvsPt_phi", "#phi->ee", kTH2D, {axis_mass, axis_ptee}, false);
      fRegistry.add("Pair/hMvsPt_jpsi", "J/#psi->ee", kTH2D, {axis_mass, axis_ptee}, false);
    }
  }

  template <typename TTrack>
  bool isTag(TTrack const& track, const float pt, const float eta, const float dcaXY, const float dcaZ, const float cYY, const float cZZ, const float cZY)
  {
    if (!track.hasITS() || !track.hasTPC()) {
      return false;
    }

    if (!isTagElectron(track)) {
      return false;
    }

    if (pt < tagCut.min_pt_track || tagCut.max_pt_track < pt) {
      return false;
    }

    if (eta < tagCut.min_eta_track || tagCut.max_eta_track < eta) {
      return false;
    }

    if (std::fabs(dcaXY) > tagCut.max_dcaxy || std::fabs(dcaZ) > tagCut.max_dcaz) {
      return false;
    }

    if (dca3DinSigmaOTF(dcaXY, dcaZ, cYY, cZZ, cZY) > tagCut.max_dca_in_sigma) {
      return false;
    }

    if (track.itsChi2NCl() < 0.f || tagCut.max_chi2its < track.itsChi2NCl()) {
      return false;
    }
    if (track.itsNCls() < tagCut.min_ncluster_its) {
      return false;
    }
    if (track.itsNClsInnerBarrel() < tagCut.min_ncluster_itsib) {
      return false;
    }

    if (track.tpcChi2NCl() < 0.f || tagCut.max_chi2tpc < track.tpcChi2NCl()) {
      return false;
    }

    if (track.tpcNClsFound() < tagCut.min_ncluster_tpc) {
      return false;
    }

    if (track.tpcNClsCrossedRows() < tagCut.min_ncrossedrows_tpc) {
      return false;
    }

    if (track.tpcCrossedRowsOverFindableCls() < tagCut.min_tpc_cr_findable_ratio) {
      return false;
    }

    if (track.tpcFractionSharedCls() > tagCut.max_frac_shared_clusters_tpc) {
      return false;
    }

    return true;
  }

  template <typename TTrack>
  bool isProbe(TTrack const& track, const float pt, const float eta, const float dcaXY, const float dcaZ)
  {
    if (!track.hasITS() || !track.hasTPC()) {
      return false;
    }

    if (!isProbeElectron(track)) {
      return false;
    }

    if (pt < probeCut.min_pt_track || probeCut.max_pt_track < pt) {
      return false;
    }

    if (eta < probeCut.min_eta_track || probeCut.max_eta_track < eta) {
      return false;
    }

    if (std::fabs(dcaXY) > probeCut.max_dcaxy) {
      return false;
    }

    if (std::fabs(dcaZ) > probeCut.max_dcaz) {
      return false;
    }

    if (track.itsChi2NCl() > probeCut.max_chi2its) {
      return false;
    }

    if (track.itsNCls() < probeCut.min_ncluster_its) {
      return false;
    }

    if (track.itsNClsInnerBarrel() < probeCut.min_ncluster_itsib) {
      return false;
    }

    if (track.tpcChi2NCl() > probeCut.max_chi2tpc) {
      return false;
    }

    if (track.tpcNClsFound() < probeCut.min_ncluster_tpc) {
      return false;
    }

    if (track.tpcNClsCrossedRows() < probeCut.min_ncrossedrows_tpc) {
      return false;
    }

    if (track.tpcCrossedRowsOverFindableCls() < probeCut.min_cr2findable_ratio_tpc) {
      return false;
    }

    if (track.tpcFractionSharedCls() > probeCut.max_frac_shared_clusters_tpc) {
      return false;
    }

    return true;
  }

  template <typename TTrack>
  bool isTagElectron(TTrack const& track)
  {
    bool is_El_TPC = tagCut.min_TPCNsigmaEl < track.tpcNSigmaEl() && track.tpcNSigmaEl() < tagCut.max_TPCNsigmaEl;
    bool is_El_TOF = tagCut.min_TOFNsigmaEl < track.tofNSigmaEl() && track.tofNSigmaEl() < tagCut.max_TOFNsigmaEl;
    return is_El_TPC && is_El_TOF;
  }

  template <typename TTrack>
  bool isProbeElectron(TTrack const& track)
  {
    bool is_El_TPC = probeCut.min_TPCNsigmaEl < track.tpcNSigmaEl() && track.tpcNSigmaEl() < probeCut.max_TPCNsigmaEl;
    bool is_El_TOF = track.hasTOF() ? probeCut.min_TOFNsigmaEl < track.tofNSigmaEl() && track.tofNSigmaEl() < probeCut.max_TOFNsigmaEl : true;
    return is_El_TPC && is_El_TOF;
  }

  template <typename TCollision>
  bool isSelectedCollision(TCollision const& collision)
  {
    if (collision.posZ() < eventCut.cfgZvtxMin || eventCut.cfgZvtxMax < collision.posZ()) {
      return false;
    }

    if (eventCut.cfgRequireFT0AND && !collision.selection_bit(o2::aod::evsel::kIsTriggerTVX)) {
      return false;
    }

    if (eventCut.cfgRequireNoTFB && !collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
      return false;
    }

    if (eventCut.cfgRequireNoITSROFB && !collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
      return false;
    }

    if (eventCut.cfgRequireVertexITSTPC && !collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) {
      return false;
    }

    if (eventCut.cfgRequireVertexTOFmatched && !collision.selection_bit(o2::aod::evsel::kIsVertexTOFmatched)) {
      return false;
    }

    if (eventCut.cfgRequireNoSameBunchPileup && !collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      return false;
    }

    if (eventCut.cfgRequireGoodZvtxFT0vsPV && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return false;
    }

    if (!(eventCut.cfgTrackOccupancyMin <= collision.trackOccupancyInTimeRange() && collision.trackOccupancyInTimeRange() < eventCut.cfgTrackOccupancyMax)) {
      return false;
    }

    if (!(eventCut.cfgFT0COccupancyMin <= collision.ft0cOccupancyInTimeRange() && collision.ft0cOccupancyInTimeRange() < eventCut.cfgFT0COccupancyMax)) {
      return false;
    }

    if (eventCut.cfgRequireGoodRCT && !rctChecker.checkTable(collision)) {
      return false;
    }

    return true;
  }

  float dca3DinSigmaOTF(const float dcaXY, const float dcaZ, const float cYY, const float cZZ, const float cZY)
  {
    float det = cYY * cZZ - cZY * cZY; // determinant
    if (det < 0) {
      return 999.f;
    } else {
      return std::sqrt(std::fabs((dcaXY * dcaXY * cZZ + dcaZ * dcaZ * cYY - 2. * dcaXY * dcaZ * cZY) / det / 2.)); // dca 3d in sigma
    }
  }

  using MyBCs = soa::Join<aod::BCsWithTimestamps, aod::BcSels>;

  using MyCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs>;

  using MyTracks = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU,
                             aod::pidTPCFullEl, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr,
                             aod::pidTOFFullEl, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFbeta>;

  Filter collisionFilter_zvtx = eventCut.cfgZvtxMin < o2::aod::collision::posZ && o2::aod::collision::posZ < eventCut.cfgZvtxMax;
  Filter collisionFilter_centrality = (eventCut.cfgCentMin < o2::aod::cent::centFT0M && o2::aod::cent::centFT0M < eventCut.cfgCentMax) || (eventCut.cfgCentMin < o2::aod::cent::centFT0A && o2::aod::cent::centFT0A < eventCut.cfgCentMax) || (eventCut.cfgCentMin < o2::aod::cent::centFT0C && o2::aod::cent::centFT0C < eventCut.cfgCentMax);
  Filter collisionFilter_numContrib = eventCut.cfgNumContribMin <= o2::aod::collision::numContrib && o2::aod::collision::numContrib < eventCut.cfgNumContribMax;
  Filter collisionFilter_track_occupancy = eventCut.cfgTrackOccupancyMin <= o2::aod::evsel::trackOccupancyInTimeRange && o2::aod::evsel::trackOccupancyInTimeRange < eventCut.cfgTrackOccupancyMax;
  Filter collisionFilter_ft0c_occupancy = eventCut.cfgFT0COccupancyMin <= o2::aod::evsel::ft0cOccupancyInTimeRange && o2::aod::evsel::ft0cOccupancyInTimeRange < eventCut.cfgFT0COccupancyMax;

  using FilteredMyCollisions = soa::Filtered<MyCollisions>;
  using FilteredMyCollision = FilteredMyCollisions::iterator;

  Filter trackFilter_itstpc = ncheckbit(aod::track::v001::detectorMap, (uint8_t)o2::aod::track::ITS) == true && ncheckbit(aod::track::v001::detectorMap, (uint8_t)o2::aod::track::TPC) == true;
  Filter trackFilter_tpcpid = probeCut.min_TPCNsigmaEl < o2::aod::pidtpc::tpcNSigmaEl && o2::aod::pidtpc::tpcNSigmaEl < probeCut.max_TPCNsigmaEl;
  using FilteredMyTracks = soa::Filtered<MyTracks>;

  SliceCache cache;
  Preslice<aod::TracksIU> perCol = o2::aod::track::collisionId;
  Partition<FilteredMyTracks> posTracks = o2::aod::track::signed1Pt > 0.f;
  Partition<FilteredMyTracks> negTracks = o2::aod::track::signed1Pt < 0.f;

  template <typename TBCs, typename TCollisions, typename TTracks>
  void runTAP(TBCs const&, TCollisions const& collisions, TTracks const&)
  {
    for (const auto& collision : collisions) {
      auto bc = collision.template bc_as<TBCs>();
      initCCDB(bc);

      float centralities[3] = {collision.centFT0M(), collision.centFT0A(), collision.centFT0C()};
      if (centralities[eventCut.cfgCentEstimator] < eventCut.cfgCentMin || eventCut.cfgCentMax < centralities[eventCut.cfgCentEstimator]) {
        continue;
      }

      fRegistry.fill(HIST("Event/hCollisionCounter"), 0);
      if (!isSelectedCollision(collision)) {
        continue;
      }
      fRegistry.fill(HIST("Event/hCollisionCounter"), 1);

      mVtx.setPos({collision.posX(), collision.posY(), collision.posZ()});
      mVtx.setCov(collision.covXX(), collision.covXY(), collision.covYY(), collision.covXZ(), collision.covYZ(), collision.covZZ());

      auto posTracks_per_coll = posTracks->sliceByCached(o2::aod::track::collisionId, collision.globalIndex(), cache);
      auto negTracks_per_coll = negTracks->sliceByCached(o2::aod::track::collisionId, collision.globalIndex(), cache);

      o2::dataformats::DCA mDcaInfoCov;
      for (const auto& tag : posTracks_per_coll) { // positron is tag
        mDcaInfoCov.set(999, 999, 999, 999, 999);
        auto tagParCov = getTrackParCov(tag);
        tagParCov.setPID(o2::track::PID::Electron);
        bool isPropOK = o2::base::Propagator::Instance()->propagateToDCABxByBz(mVtx, tagParCov, 2.f, matCorr, &mDcaInfoCov);
        if (!isPropOK) {
          continue;
        }
        if (!isTag(tag, tagParCov.getPt(), tagParCov.getEta(), mDcaInfoCov.getY(), mDcaInfoCov.getZ(), tagParCov.getSigmaY2(), tagParCov.getSigmaZ2(), tagParCov.getSigmaZY())) {
          continue;
        }

        for (const auto& probe : negTracks_per_coll) { // electron is probe.
          mDcaInfoCov.set(999, 999, 999, 999, 999);
          auto probeParCov = getTrackParCov(probe);
          probeParCov.setPID(o2::track::PID::Electron);
          bool isPropOK = o2::base::Propagator::Instance()->propagateToDCABxByBz(mVtx, probeParCov, 2.f, matCorr, &mDcaInfoCov);
          if (!isPropOK) {
            continue;
          }
          if (!isProbe(probe, probeParCov.getPt(), probeParCov.getEta(), mDcaInfoCov.getY(), mDcaInfoCov.getZ())) {
            continue;
          }

          ROOT::Math::PtEtaPhiMVector v1(tagParCov.getPt(), tagParCov.getEta(), RecoDecay::constrainAngle(tagParCov.getPhi(), 0, 1U), o2::constants::physics::MassElectron);
          ROOT::Math::PtEtaPhiMVector v2(probeParCov.getPt(), probeParCov.getEta(), RecoDecay::constrainAngle(probeParCov.getPhi(), 0, 1U), o2::constants::physics::MassElectron);
          ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
          float mee = v12.M();
          float phiv = o2::aod::pwgem::dilepton::utils::pairutil::getPhivPair(v1.Px(), v1.Py(), v1.Pz(), v2.Px(), v2.Py(), v2.Pz(), tag.sign(), probe.sign(), d_bz);
          fRegistry.fill(HIST("Pair/hMvsPhiV"), phiv, mee);

          if ((pairCut.minMee < mee && mee < pairCut.maxMee) && (pairCut.minPhiV < phiv && phiv < pairCut.maxPhiV)) {
            fRegistry.fill(HIST("Track/hs"), probeParCov.getPt(), probeParCov.getEta(), RecoDecay::constrainAngle(probeParCov.getPhi(), 0, 1U), probe.sign(), mDcaInfoCov.getY() * 1e+4, mDcaInfoCov.getZ() * 1e+4, mDcaInfoCov.getY() / std::sqrt(probeParCov.getSigmaY2()), mDcaInfoCov.getZ() / std::sqrt(probeParCov.getSigmaZ2()));
            fRegistry.fill(HIST("Track/hTPCNsigmaEl"), probe.tpcInnerParam(), probe.tpcNSigmaEl());
          }

        } // end of electron loop
      } // end of positron loop

      for (const auto& tag : negTracks_per_coll) { // electron is tag
        mDcaInfoCov.set(999, 999, 999, 999, 999);
        auto tagParCov = getTrackParCov(tag);
        tagParCov.setPID(o2::track::PID::Electron);
        bool isPropOK = o2::base::Propagator::Instance()->propagateToDCABxByBz(mVtx, tagParCov, 2.f, matCorr, &mDcaInfoCov);
        if (!isPropOK) {
          continue;
        }
        if (!isTag(tag, tagParCov.getPt(), tagParCov.getEta(), mDcaInfoCov.getY(), mDcaInfoCov.getZ(), tagParCov.getSigmaY2(), tagParCov.getSigmaZ2(), tagParCov.getSigmaZY())) {
          continue;
        }

        for (const auto& probe : posTracks_per_coll) { // positron is probe.
          mDcaInfoCov.set(999, 999, 999, 999, 999);
          auto probeParCov = getTrackParCov(probe);
          probeParCov.setPID(o2::track::PID::Electron);
          bool isPropOK = o2::base::Propagator::Instance()->propagateToDCABxByBz(mVtx, probeParCov, 2.f, matCorr, &mDcaInfoCov);
          if (!isPropOK) {
            continue;
          }
          if (!isProbe(probe, probeParCov.getPt(), probeParCov.getEta(), mDcaInfoCov.getY(), mDcaInfoCov.getZ())) {
            continue;
          }

          ROOT::Math::PtEtaPhiMVector v1(tagParCov.getPt(), tagParCov.getEta(), RecoDecay::constrainAngle(tagParCov.getPhi(), 0, 1U), o2::constants::physics::MassElectron);
          ROOT::Math::PtEtaPhiMVector v2(probeParCov.getPt(), probeParCov.getEta(), RecoDecay::constrainAngle(probeParCov.getPhi(), 0, 1U), o2::constants::physics::MassElectron);
          ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
          float mee = v12.M();
          float phiv = o2::aod::pwgem::dilepton::utils::pairutil::getPhivPair(v1.Px(), v1.Py(), v1.Pz(), v2.Px(), v2.Py(), v2.Pz(), tag.sign(), probe.sign(), d_bz);
          fRegistry.fill(HIST("Pair/hMvsPhiV"), phiv, mee);

          if ((pairCut.minMee < mee && mee < pairCut.maxMee) && (pairCut.minPhiV < phiv && phiv < pairCut.maxPhiV)) {
            fRegistry.fill(HIST("Track/hs"), probeParCov.getPt(), probeParCov.getEta(), RecoDecay::constrainAngle(probeParCov.getPhi(), 0, 1U), probe.sign(), mDcaInfoCov.getY() * 1e+4, mDcaInfoCov.getZ() * 1e+4, mDcaInfoCov.getY() / std::sqrt(probeParCov.getSigmaY2()), mDcaInfoCov.getZ() / std::sqrt(probeParCov.getSigmaZ2()));
            fRegistry.fill(HIST("Track/hTPCNsigmaEl"), probe.tpcInnerParam(), probe.tpcNSigmaEl());
          }

        } // end of positron loop
      } // end of electron loop

    } // end of collision loop
  } // end of runTAP

  void processData(MyBCs const& bcs, FilteredMyCollisions const& collisions, FilteredMyTracks const& tracks)
  {
    runTAP(bcs, collisions, tracks);
  }
  PROCESS_SWITCH(createTTP, processData, "process data", true);

  struct lepton {
    float pt{0};
    float eta{0};
    float phi{0};
    int mcParticleId{-1};
  };

  using MyCollisionsMC = soa::Join<MyCollisions, aod::McCollisionLabels>;
  using MyTracksMC = soa::Join<MyTracks, aod::McTrackLabels>;

  using FilteredMyCollisionsMC = soa::Filtered<MyCollisionsMC>;
  using FilteredMyTracksMC = soa::Filtered<MyTracksMC>;

  void processMC(MyBCs const&, FilteredMyCollisionsMC const& collisions, FilteredMyTracksMC const& tracks, aod::McCollisions const&, aod::McParticles const& mcParticles)
  {
    for (const auto& collision : collisions) {
      if (!collision.has_mcCollision()) {
        continue;
      }

      auto bc = collision.template bc_as<MyBCs>();
      initCCDB(bc);

      float centralities[3] = {collision.centFT0M(), collision.centFT0A(), collision.centFT0C()};
      if (centralities[eventCut.cfgCentEstimator] < eventCut.cfgCentMin || eventCut.cfgCentMax < centralities[eventCut.cfgCentEstimator]) {
        continue;
      }

      fRegistry.fill(HIST("Event/hCollisionCounter"), 0);
      if (!isSelectedCollision(collision)) {
        continue;
      }
      fRegistry.fill(HIST("Event/hCollisionCounter"), 1);

      mVtx.setPos({collision.posX(), collision.posY(), collision.posZ()});
      mVtx.setCov(collision.covXX(), collision.covXY(), collision.covYY(), collision.covXZ(), collision.covYZ(), collision.covZZ());

      std::vector<lepton> posLeptons;
      std::vector<lepton> negLeptons;

      o2::dataformats::DCA mDcaInfoCov;
      auto tracks_per_coll = tracks.sliceBy(perCol, collision.globalIndex());
      for (const auto& track : tracks_per_coll) {
        if (!track.has_mcParticle()) {
          continue;
        }
        auto mcParticle = track.template mcParticle_as<aod::McParticles>();
        auto mcCollision = mcParticle.template mcCollision_as<aod::McCollisions>();
        if (std::abs(mcParticle.pdgCode()) != 11) {
          continue;
        }
        if (!(mcParticle.isPhysicalPrimary() || mcParticle.producedByGenerator())) {
          continue;
        }
        if (!mcParticle.has_mothers()) {
          continue;
        }

        if (cfgEventGeneratorType >= 0 && mcCollision.getSubGeneratorId() != cfgEventGeneratorType) {
          continue;
        }

        mDcaInfoCov.set(999, 999, 999, 999, 999);
        auto trackParCov = getTrackParCov(track);
        trackParCov.setPID(o2::track::PID::Electron);
        bool isPropOK = o2::base::Propagator::Instance()->propagateToDCABxByBz(mVtx, trackParCov, 2.f, matCorr, &mDcaInfoCov);
        if (!isPropOK) {
          continue;
        }
        if (!isProbe(track, trackParCov.getPt(), trackParCov.getEta(), mDcaInfoCov.getY(), mDcaInfoCov.getZ())) {
          continue;
        }

        fRegistry.fill(HIST("Track/hs"), trackParCov.getPt(), trackParCov.getEta(), RecoDecay::constrainAngle(trackParCov.getPhi(), 0, 1U), track.sign(), mDcaInfoCov.getY() * 1e+4, mDcaInfoCov.getZ() * 1e+4, mDcaInfoCov.getY() / std::sqrt(trackParCov.getSigmaY2()), mDcaInfoCov.getZ() / std::sqrt(trackParCov.getSigmaZ2()));

        lepton tmp;
        tmp.pt = trackParCov.getPt();
        tmp.eta = trackParCov.getEta();
        tmp.phi = RecoDecay::constrainAngle(trackParCov.getPhi(), 0, 1U);
        tmp.mcParticleId = mcParticle.globalIndex();

        if (track.sign() > 0) {
          posLeptons.emplace_back(tmp);
        } else {
          negLeptons.emplace_back(tmp);
        }
      } // end of track loop

      for (const auto& pos : posLeptons) {
        for (const auto& neg : negLeptons) {
          ROOT::Math::PtEtaPhiMVector v1(pos.pt, pos.eta, pos.phi, o2::constants::physics::MassElectron);
          ROOT::Math::PtEtaPhiMVector v2(neg.pt, neg.eta, neg.phi, o2::constants::physics::MassElectron);
          ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;

          auto posmc = mcParticles.rawIteratorAt(pos.mcParticleId);
          auto negmc = mcParticles.rawIteratorAt(neg.mcParticleId);

          int omegaId = o2::aod::pwgem::dilepton::utils::mcutil::FindCommonMotherFrom2Prongs(posmc, negmc, -11, 11, 223, mcParticles);
          int phiId = o2::aod::pwgem::dilepton::utils::mcutil::FindCommonMotherFrom2Prongs(posmc, negmc, -11, 11, 333, mcParticles);
          int jpsiId = o2::aod::pwgem::dilepton::utils::mcutil::FindCommonMotherFrom2Prongs(posmc, negmc, -11, 11, 443, mcParticles);
          if (omegaId > 0) {
            auto mcMother = mcParticles.rawIteratorAt(omegaId);
            int ndau = mcMother.daughtersIds()[1] - mcMother.daughtersIds()[0] + 1;
            if (ndau == 2) {
              fRegistry.fill(HIST("Pair/hMvsPt_omega"), v12.M(), v12.Pt());
            }
          } else if (phiId > 0) {
            auto mcMother = mcParticles.rawIteratorAt(phiId);
            int ndau = mcMother.daughtersIds()[1] - mcMother.daughtersIds()[0] + 1;
            if (ndau == 2) {
              fRegistry.fill(HIST("Pair/hMvsPt_phi"), v12.M(), v12.Pt());
            }
          } else if (jpsiId > 0) {
            auto mcMother = mcParticles.rawIteratorAt(jpsiId);
            int ndau = mcMother.daughtersIds()[1] - mcMother.daughtersIds()[0] + 1;
            if (ndau == 2) {
              fRegistry.fill(HIST("Pair/hMvsPt_jpsi"), v12.M(), v12.Pt());
            }
          }
        } // end of electron loop
      } // end of positron loop

      posLeptons.clear();
      posLeptons.shrink_to_fit();
      negLeptons.clear();
      negLeptons.shrink_to_fit();

    } // end of collision loop
  }
  PROCESS_SWITCH(createTTP, processMC, "process MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<createTTP>(cfgc, TaskName{"create-ttp"})};
}
