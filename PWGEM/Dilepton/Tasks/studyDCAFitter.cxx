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

/// \file studyDCAFitter.cxx
/// \brief a task to study tagging e from charm hadron decays in MC
/// \author daiki.sekihata@cern.ch

#include "PWGEM/Dilepton/DataModel/lmeeMLTables.h"
#include "PWGEM/Dilepton/Utils/MCUtilities.h"
#include "PWGEM/Dilepton/Utils/PairUtilities.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/PhysicsConstants.h>
#include <DCAFitter/DCAFitterN.h>
// #include <DataFormatsCalibration/MeanVertexObject.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <DataFormatsParameters/GRPObject.h>
#include <DetectorsBase/MatLayerCylSet.h>
#include <DetectorsBase/Propagator.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>
#include <ReconstructionDataFormats/DCA.h>
#include <ReconstructionDataFormats/PID.h>

#include <Math/Vector4D.h> // IWYU pragma: keep (do not replace with Math/Vector4Dfwd.h)
#include <Math/Vector4Dfwd.h>
#include <TH1.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <random>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::soa;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;
using namespace o2::aod::pwgem::dilepton::utils::mcutil;
using namespace o2::aod::pwgem::dilepton::utils::pairutil;

struct studyDCAFitter {
  using MyCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::McCollisionLabels>;
  using MyTracks = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::McTrackLabels>;
  using MyBCs = soa::Join<aod::BCsWithTimestamps, aod::BcSels>;

  Produces<aod::EMMLEvents> eventTable;
  Produces<aod::EMMLDielectronsAtSV> dileptonTable;

  // Configurables
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
  // Configurable<std::string> mVtxPath{"mVtxPath", "GLO/Calib/MeanVertex", "Path of the mean vertex file"};
  Configurable<float> cfgDownSampling{"cfgDownSampling", 1.1, "down sampling for wrongly found SV"};

  struct : ConfigurableGroup {
    std::string prefix = "electronCut";
    Configurable<float> cfg_min_pt_track{"cfg_min_pt_track", 0.4, "min pT for single track"};
    Configurable<float> cfg_max_pt_track{"cfg_max_pt_track", 1e+10, "max pT for single track"};
    Configurable<float> cfg_min_eta_track{"cfg_min_eta_track", -0.8, "min eta for single track"};
    Configurable<float> cfg_max_eta_track{"cfg_max_eta_track", +0.8, "max eta for single track"};
    Configurable<float> cfg_min_cr2findable_ratio_tpc{"cfg_min_cr2findable_ratio_tpc", 0.8, "min. TPC Ncr/Nf ratio"};
    Configurable<float> cfg_max_frac_shared_clusters_tpc{"cfg_max_frac_shared_clusters_tpc", 0.7, "max fraction of shared clusters in TPC"};
    Configurable<int> cfg_min_ncrossedrows_tpc{"cfg_min_ncrossedrows_tpc", 80, "min ncrossed rows"};
    Configurable<int> cfg_min_ncluster_tpc{"cfg_min_ncluster_tpc", 0, "min ncluster tpc"};
    Configurable<int> cfg_min_ncluster_its{"cfg_min_ncluster_its", 5, "min ncluster its"};
    Configurable<int> cfg_min_ncluster_itsib{"cfg_min_ncluster_itsib", 1, "min ncluster itsib"};
    Configurable<float> cfg_max_chi2tpc{"cfg_max_chi2tpc", 4.0, "max chi2/NclsTPC"};
    Configurable<float> cfg_max_chi2its{"cfg_max_chi2its", 5.0, "max chi2/NclsITS"};
    Configurable<float> cfg_max_dcaxy{"cfg_max_dcaxy", 1.0, "max dca XY for single track in cm"};
    Configurable<float> cfg_max_dcaz{"cfg_max_dcaz", 1.0, "max dca Z for single track in cm"};
  } electronCut;

  struct : ConfigurableGroup {
    std::string prefix = "svCut";
    Configurable<bool> d_UseAbsDCA{"d_UseAbsDCA", true, "Use Abs DCAs"};
    Configurable<bool> d_UseWeightedPCA{"d_UseWeightedPCA", false, "Vertices use cov matrices"};
    Configurable<float> cfg_max_chi2PCA{"cfg_max_chi2PCA", 1.0, "max chi2PCA"};
  } svCut;

  struct : ConfigurableGroup {
    std::string prefix = "eventCut";
    Configurable<int> cfgEventGeneratorType{"cfgEventGeneratorType", -1, "if positive, select event generator type. i.e. gap or signal"};
    Configurable<int> cfgRejectEventGenerator{"cfgRejectEventGenerator", 999, "reject event generator. e.g. reject tracks from gap events"};
    Configurable<int> cfgCentEstimator{"cfgCentEstimator", 2, "FT0M:0, FT0A:1, FT0C:2"};
    Configurable<float> cfgCentMin{"cfgCentMin", -1.f, "min. centrality"};
    Configurable<float> cfgCentMax{"cfgCentMax", 999.f, "max. centrality"};
    Configurable<float> cfgZvtxMin{"cfgZvtxMin", -10.f, "min. Zvtx"};
    Configurable<float> cfgZvtxMax{"cfgZvtxMax", 10.f, "max. Zvtx"};
    Configurable<bool> cfgRequireFT0AND{"cfgRequireFT0AND", true, "require FT0AND"};
    Configurable<bool> cfgRequireNoTFB{"cfgRequireNoTFB", true, "require No time frame border"};
    Configurable<bool> cfgRequireNoITSROFB{"cfgRequireNoITSROFB", false, "require no ITS readout frame border"};
    Configurable<bool> cfgRequireNoSameBunchPileup{"cfgRequireNoSameBunchPileup", false, "require no same bunch pileup in event cut"};
    Configurable<bool> cfgRequireGoodZvtxFT0vsPV{"cfgRequireGoodZvtxFT0vsPV", false, "require good Zvtx between FT0 vs. PV in event cut"};
    // for RCT
    o2::framework::Configurable<bool> cfgRequireGoodRCT{"cfgRequireGoodRCT", true, "require good detector flag in run condtion table"};
    o2::framework::Configurable<std::string> cfgRCTLabel{"cfgRCTLabel", "CBT_hadronPID", "select 1 [CBT, CBT_hadronPID] see O2Physics/Common/CCDB/RCTSelectionFlags.h"};
    o2::framework::Configurable<bool> cfgCheckZDC{"cfgCheckZDC", false, "set ZDC flag for AA"};
    o2::framework::Configurable<bool> cfgTreatLimitedAcceptanceAsBad{"cfgTreatLimitedAcceptanceAsBad", false, "reject all events where the detectors relevant for the specified Runlist are flagged as LimitedAcceptance"};
  } eventCut;

  o2::framework::ConfigurableAxis ConfMllBins{"ConfMllBins", {o2::framework::VARIABLE_WIDTH, 0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.30, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.40, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.50, 0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59, 0.60, 0.61, 0.62, 0.63, 0.64, 0.65, 0.66, 0.67, 0.68, 0.69, 0.70, 0.71, 0.72, 0.73, 0.74, 0.75, 0.76, 0.77, 0.78, 0.79, 0.80, 0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 0.87, 0.88, 0.89, 0.90, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 1.00, 1.01, 1.02, 1.03, 1.04, 1.05, 1.06, 1.07, 1.08, 1.09, 1.10, 1.11, 1.12, 1.13, 1.14, 1.15, 1.16, 1.17, 1.18, 1.19, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00, 2.10, 2.20, 2.30, 2.40, 2.50, 2.60, 2.70, 2.75, 2.80, 2.85, 2.90, 2.95, 3.00, 3.05, 3.10, 3.15, 3.20, 3.25, 3.30, 3.35, 3.40, 3.45, 3.50, 3.55, 3.60, 3.65, 3.70, 3.75, 3.80, 3.85, 3.90, 3.95, 4.00, 4.50, 5.00, 5.50, 6.00, 6.50, 7.00, 7.50, 8.00}, "mll bins for output histograms"};
  o2::framework::ConfigurableAxis ConfPtllBins{"ConfPtllBins", {o2::framework::VARIABLE_WIDTH, 0.00, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00, 2.50, 3.00, 3.50, 4.00, 4.50, 5.00, 6.00, 7.00, 8.00, 9.00, 10.00}, "pTll bins for output histograms"};
  o2::framework::ConfigurableAxis ConfDCAllBins{"ConfDCAllBins", {o2::framework::VARIABLE_WIDTH, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0}, "DCAll bins for output histograms"};

  HistogramRegistry fRegistry{"fRegistry"};

  void init(o2::framework::InitContext&)
  {
    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);

    fitter.setPropagateToPCA(true);
    fitter.setMaxR(200.f);
    fitter.setMinParamChange(1e-3);
    fitter.setMinRelChi2Change(0.9);
    fitter.setMaxDZIni(1e9);
    fitter.setMaxChi2(1e9);
    fitter.setUseAbsDCA(svCut.d_UseAbsDCA);
    fitter.setWeightedFinalPCA(svCut.d_UseWeightedPCA);
    fitter.setMatCorrType(matCorr);

    mRunNumber = 0;
    d_bz = 0;

    std::random_device seed_gen;
    engine = std::mt19937(seed_gen());
    dist01 = std::uniform_real_distribution<float>(0.0f, 1.0f);

    addHistograms();
  }

  int mRunNumber{0};
  float d_bz{0};
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  // o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;
  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;
  // const o2::dataformats::MeanVertexObject* mMeanVtx = nullptr;
  o2::base::MatLayerCylSet* lut = nullptr;
  o2::vertexing::DCAFitterN<2> fitter;
  o2::dataformats::DCA mDcaInfoCov;
  o2::dataformats::VertexBase mVtx;

  std::mt19937 engine;
  std::uniform_real_distribution<float> dist01;

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
      // mMeanVtx = ccdb->getForTimeStamp<o2::dataformats::MeanVertexObject>(mVtxPath, bc.timestamp());
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
      // mMeanVtx = ccdb->getForTimeStamp<o2::dataformats::MeanVertexObject>(mVtxPath, bc.timestamp());

      // Fetch magnetic field from ccdb for current collision
      d_bz = std::lround(5.f * grpmag->getL3Current() / 30000.f);
      LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << d_bz << " kZG";
    }
    mRunNumber = bc.runNumber();
    fitter.setBz(d_bz);
  }

  void addHistograms()
  {
    auto hCollisionCounter = fRegistry.add<TH1>("Event/hCollisionCounter", "collision counter", kTH1D, {{5, -0.5f, 4.5f}}, false);
    hCollisionCounter->GetXaxis()->SetBinLabel(1, "all");
    hCollisionCounter->GetXaxis()->SetBinLabel(2, "accepted");

    fRegistry.add("Event/hZvtx", "vertex z; Z_{vtx} (cm)", kTH1F, {{100, -50, +50}}, false);
    fRegistry.add("Event/hMultNTracksPV", "hMultNTracksPV; N_{track} to PV", kTH1F, {{6001, -0.5, 6000.5}}, false);
    fRegistry.add("Event/hMultNTracksPVeta1", "hMultNTracksPVeta1; N_{track} to PV", kTH1F, {{6001, -0.5, 6000.5}}, false);
    fRegistry.add("Event/hMultFT0", "hMultFT0;mult. FT0A;mult. FT0C", kTH2F, {{200, 0, 200000}, {60, 0, 60000}}, false);
    fRegistry.add("Event/hCentFT0A", "hCentFT0A;centrality FT0A (%)", kTH1F, {{110, 0, 110}}, false);
    fRegistry.add("Event/hCentFT0C", "hCentFT0C;centrality FT0C (%)", kTH1F, {{110, 0, 110}}, false);
    fRegistry.add("Event/hCentFT0M", "hCentFT0M;centrality FT0M (%)", kTH1F, {{110, 0, 110}}, false);
    fRegistry.add("Event/hCentFT0CvsMultNTracksPV", "hCentFT0CvsMultNTracksPV;centrality FT0C (%);N_{track} to PV", kTH2F, {{110, 0, 110}, {600, 0, 6000}}, false);
    fRegistry.add("Event/hMultFT0CvsMultNTracksPV", "hMultFT0CvsMultNTracksPV;mult. FT0C;N_{track} to PV", kTH2F, {{60, 0, 60000}, {600, 0, 6000}}, false);

    const o2::framework::AxisSpec axis_mass{ConfMllBins, "m_{ll} (GeV/c^{2})"};
    const o2::framework::AxisSpec axis_pt{ConfPtllBins, "p_{T,ee} (GeV/c)"};
    const o2::framework::AxisSpec axis_dca{ConfDCAllBins, "DCA_{ee}^{3D} (#sigma)"};

    // for pairs
    fRegistry.add("Pair/PV/Zboson/uls/hs", "hs;m_{ee} (GeV/c^{2});p_{T,ee} (GeV/c);DCA_{ee}^{3D} (#sigma);", kTHnSparseF, {axis_mass, axis_pt, axis_dca}, false);
    fRegistry.addClone("Pair/PV/Zboson/uls/", "Pair/PV/Zboson/lspp/");
    fRegistry.addClone("Pair/PV/Zboson/uls/", "Pair/PV/Zboson/lsmm/");
    fRegistry.addClone("Pair/PV/Zboson/", "Pair/PV/PromptJpsi/");
    fRegistry.addClone("Pair/PV/Zboson/", "Pair/PV/NonPromptJpsi/");
    fRegistry.addClone("Pair/PV/Zboson/", "Pair/PV/c2e_c2e/");
    fRegistry.addClone("Pair/PV/Zboson/", "Pair/PV/b2e_b2e/");
    fRegistry.addClone("Pair/PV/Zboson/", "Pair/PV/b2c2e_b2c2e/");
    fRegistry.addClone("Pair/PV/Zboson/", "Pair/PV/b2c2e_b2e_sameb/");
    fRegistry.addClone("Pair/PV/Zboson/", "Pair/PV/b2c2e_b2e_diffb/");

    fRegistry.add("Pair/SV/Zboson/uls/hs", "hs;m_{ee} (GeV/c^{2});p_{T,ee} (GeV/c);DCA_{ee}^{3D} (#sigma);", kTHnSparseF, {axis_mass, axis_pt}, false);
    fRegistry.add("Pair/SV/Zboson/uls/hCosPA", "cosPA;cosPA;", kTH1F, {{200, -1, 1}}, false);
    fRegistry.add("Pair/SV/Zboson/uls/hChi2PCA", "chi2 at PCA;log_{10}(#chi^{2}_{PCA});", kTH1F, {{1000, -10, 0}}, false);
    fRegistry.addClone("Pair/SV/Zboson/uls/", "Pair/SV/Zboson/lspp/");
    fRegistry.addClone("Pair/SV/Zboson/uls/", "Pair/SV/Zboson/lsmm/");
    fRegistry.addClone("Pair/SV/Zboson/", "Pair/SV/PromptJpsi/");
    fRegistry.addClone("Pair/SV/Zboson/", "Pair/SV/NonPromptJpsi/");
    fRegistry.addClone("Pair/SV/Zboson/", "Pair/SV/c2e_c2e/");
    fRegistry.addClone("Pair/SV/Zboson/", "Pair/SV/b2e_b2e/");
    fRegistry.addClone("Pair/SV/Zboson/", "Pair/SV/b2c2e_b2c2e/");
    fRegistry.addClone("Pair/SV/Zboson/", "Pair/SV/b2c2e_b2e_sameb/");
    fRegistry.addClone("Pair/SV/Zboson/", "Pair/SV/b2c2e_b2e_diffb/");
  }

  template <typename TCollision>
  bool isSelectedCollision(TCollision const& collision)
  {
    if (!(eventCut.cfgZvtxMin < collision.posZ() && collision.posZ() < eventCut.cfgZvtxMax)) {
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
    if (eventCut.cfgRequireNoSameBunchPileup && !collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      return false;
    }
    if (eventCut.cfgRequireGoodZvtxFT0vsPV && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return false;
    }
    return true;
  }

  template <typename TTrack, typename TTrackParCov>
  bool isSelectedTrack(TTrack const& track, TTrackParCov const& trackParCov, const float dcaXY, const float dcaZ)
  {

    if (!track.hasITS() || !track.hasTPC()) {
      return false;
    }

    if (trackParCov.getPt() < electronCut.cfg_min_pt_track || electronCut.cfg_max_pt_track < trackParCov.getPt()) {
      return false;
    }

    if (trackParCov.getEta() < electronCut.cfg_min_eta_track || electronCut.cfg_max_eta_track < trackParCov.getEta()) {
      return false;
    }

    if (std::fabs(dcaXY) > electronCut.cfg_max_dcaxy) {
      return false;
    }

    if (std::fabs(dcaZ) > electronCut.cfg_max_dcaz) {
      return false;
    }

    if (track.itsChi2NCl() > electronCut.cfg_max_chi2its) {
      return false;
    }

    if (track.itsNCls() < electronCut.cfg_min_ncluster_its) {
      return false;
    }

    if (track.itsNClsInnerBarrel() < electronCut.cfg_min_ncluster_itsib) {
      return false;
    }

    if (track.tpcChi2NCl() > electronCut.cfg_max_chi2tpc) {
      return false;
    }

    if (track.tpcNClsFound() < electronCut.cfg_min_ncluster_tpc) {
      return false;
    }

    if (track.tpcNClsCrossedRows() < electronCut.cfg_min_ncrossedrows_tpc) {
      return false;
    }

    if (track.tpcCrossedRowsOverFindableCls() < electronCut.cfg_min_cr2findable_ratio_tpc) {
      return false;
    }

    if (track.tpcFractionSharedCls() > electronCut.cfg_max_frac_shared_clusters_tpc) {
      return false;
    }

    return true;
  }

  template <typename TCollision>
  void fillEventHistograms(TCollision const& collision)
  {
    fRegistry.fill(HIST("Event/hZvtx"), collision.posZ());
    fRegistry.fill(HIST("Event/hMultNTracksPV"), collision.multNTracksPV());
    fRegistry.fill(HIST("Event/hMultNTracksPVeta1"), collision.multNTracksPVeta1());
    fRegistry.fill(HIST("Event/hMultFT0"), collision.multFT0A(), collision.multFT0C());
    fRegistry.fill(HIST("Event/hCentFT0A"), collision.centFT0A());
    fRegistry.fill(HIST("Event/hCentFT0C"), collision.centFT0C());
    fRegistry.fill(HIST("Event/hCentFT0M"), collision.centFT0M());
    fRegistry.fill(HIST("Event/hCentFT0CvsMultNTracksPV"), collision.centFT0C(), collision.multNTracksPV());
    fRegistry.fill(HIST("Event/hMultFT0CvsMultNTracksPV"), collision.multFT0C(), collision.multNTracksPV());
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

  template <typename TTrack, typename TMCParticles>
  int FindCommonMother(TTrack const& posmc, TTrack const& negmc, TMCParticles const& mcparticles)
  {
    int arr[] = {
      FindCommonMotherFrom2Prongs(posmc, negmc, -11, 11, 23, mcparticles), // Z/gamma*
      FindCommonMotherFrom2Prongs(posmc, negmc, -11, 11, 111, mcparticles),
      FindCommonMotherFrom2Prongs(posmc, negmc, -11, 11, 221, mcparticles),
      FindCommonMotherFrom2Prongs(posmc, negmc, -11, 11, 331, mcparticles),
      FindCommonMotherFrom2Prongs(posmc, negmc, -11, 11, 113, mcparticles),
      FindCommonMotherFrom2Prongs(posmc, negmc, -11, 11, 223, mcparticles),
      FindCommonMotherFrom2Prongs(posmc, negmc, -11, 11, 333, mcparticles),
      FindCommonMotherFrom2Prongs(posmc, negmc, -11, 11, 443, mcparticles),
      FindCommonMotherFrom2Prongs(posmc, negmc, -11, 11, 100443, mcparticles),
      FindCommonMotherFrom2Prongs(posmc, negmc, -11, 11, 553, mcparticles),
      FindCommonMotherFrom2Prongs(posmc, negmc, -11, 11, 100553, mcparticles),
      FindCommonMotherFrom2Prongs(posmc, negmc, -11, 11, 200553, mcparticles),
      FindCommonMotherFrom2Prongs(posmc, negmc, -11, 11, 300553, mcparticles)};
    int size = sizeof(arr) / sizeof(*arr);
    int max = *std::max_element(arr, arr + size);
    return max;
  }

  template <uint8_t signType, typename TTrack, typename TMCParticles>
  void runPairingAtPV(TTrack const& t1, TTrack const& t2, TMCParticles const& mcParticles)
  {
    mDcaInfoCov.set(999, 999, 999, 999, 999);
    auto trackParCov1 = getTrackParCov(t1);
    trackParCov1.setPID(o2::track::PID::Electron);
    o2::base::Propagator::Instance()->propagateToDCABxByBz(mVtx, trackParCov1, 2.f, matCorr, &mDcaInfoCov);
    float dcaXY1 = mDcaInfoCov.getY();
    float dcaZ1 = mDcaInfoCov.getZ();

    mDcaInfoCov.set(999, 999, 999, 999, 999);
    auto trackParCov2 = getTrackParCov(t2);
    trackParCov2.setPID(o2::track::PID::Electron);
    o2::base::Propagator::Instance()->propagateToDCABxByBz(mVtx, trackParCov2, 2.f, matCorr, &mDcaInfoCov);
    float dcaXY2 = mDcaInfoCov.getY();
    float dcaZ2 = mDcaInfoCov.getZ();

    ROOT::Math::PtEtaPhiMVector v1(trackParCov1.getPt(), trackParCov1.getEta(), trackParCov1.getPhi(), o2::constants::physics::MassElectron);
    ROOT::Math::PtEtaPhiMVector v2(trackParCov2.getPt(), trackParCov2.getEta(), trackParCov2.getPhi(), o2::constants::physics::MassElectron);
    ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
    float dca3DinSigma1 = dca3DinSigmaOTF(dcaXY1, dcaZ1, trackParCov1.getSigmaY2(), trackParCov1.getSigmaZ2(), trackParCov1.getSigmaZY());
    float dca3DinSigma2 = dca3DinSigmaOTF(dcaXY2, dcaZ2, trackParCov2.getSigmaY2(), trackParCov2.getSigmaZ2(), trackParCov2.getSigmaZY());
    float pair_dca = pairDCAQuadSum(dca3DinSigma1, dca3DinSigma2);

    auto t1mc = t1.template mcParticle_as<aod::McParticles>();
    auto t2mc = t2.template mcParticle_as<aod::McParticles>();

    int mcCommonMotherId = FindCommonMother(t1mc, t2mc, mcParticles);
    int hfee_type = IsHF(t1mc, t2mc, mcParticles);

    if (mcCommonMotherId > -1) {
      auto cmp = mcParticles.rawIteratorAt(mcCommonMotherId);
      switch (std::abs(cmp.pdgCode())) {
        case 23:
          if constexpr (signType == 0) { // ULS
            fRegistry.fill(HIST("Pair/PV/Zboson/uls/hs"), v12.M(), v12.Pt(), pair_dca);
          } else if constexpr (signType == 1) { // LS++
            fRegistry.fill(HIST("Pair/PV/Zboson/lspp/hs"), v12.M(), v12.Pt(), pair_dca);
          } else if constexpr (signType == 2) { // LS--
            fRegistry.fill(HIST("Pair/PV/Zboson/lsmm/hs"), v12.M(), v12.Pt(), pair_dca);
          }
          break;
        case 443:
          if (IsFromCharm(cmp, mcParticles) < 0 && IsFromBeauty(cmp, mcParticles) < 0) { // prompt
            if constexpr (signType == 0) {                                               // ULS
              fRegistry.fill(HIST("Pair/PV/PromptJpsi/uls/hs"), v12.M(), v12.Pt(), pair_dca);
            } else if constexpr (signType == 1) { // LS++
              fRegistry.fill(HIST("Pair/PV/PromptJpsi/lspp/hs"), v12.M(), v12.Pt(), pair_dca);
            } else if constexpr (signType == 2) { // LS--
              fRegistry.fill(HIST("Pair/PV/PromptJpsi/lsmm/hs"), v12.M(), v12.Pt(), pair_dca);
            }
          } else {                         // nonprompt
            if constexpr (signType == 0) { // ULS
              fRegistry.fill(HIST("Pair/PV/NonPromptJpsi/uls/hs"), v12.M(), v12.Pt(), pair_dca);
            } else if constexpr (signType == 1) { // LS++
              fRegistry.fill(HIST("Pair/PV/NonPromptJpsi/lspp/hs"), v12.M(), v12.Pt(), pair_dca);
            } else if constexpr (signType == 2) { // LS--
              fRegistry.fill(HIST("Pair/PV/NonPromptJpsi/lsmm/hs"), v12.M(), v12.Pt(), pair_dca);
            }
          }
          break;
        default:
          break;
      } // end of switch for LF
    } else if (hfee_type > -1) {
      switch (hfee_type) {
        case static_cast<int>(EM_HFeeType::kCe_Ce):
          if constexpr (signType == 0) { // ULS
            fRegistry.fill(HIST("Pair/PV/c2e_c2e/uls/hs"), v12.M(), v12.Pt(), pair_dca);
          } else if constexpr (signType == 1) { // LS++
            fRegistry.fill(HIST("Pair/PV/c2e_c2e/lspp/hs"), v12.M(), v12.Pt(), pair_dca);
          } else if constexpr (signType == 2) { // LS--
            fRegistry.fill(HIST("Pair/PV/c2e_c2e/lsmm/hs"), v12.M(), v12.Pt(), pair_dca);
          }
          break;
        case static_cast<int>(EM_HFeeType::kBe_Be):
          if constexpr (signType == 0) { // ULS
            fRegistry.fill(HIST("Pair/PV/b2e_b2e/uls/hs"), v12.M(), v12.Pt(), pair_dca);
          } else if constexpr (signType == 1) { // LS++
            fRegistry.fill(HIST("Pair/PV/b2e_b2e/lspp/hs"), v12.M(), v12.Pt(), pair_dca);
          } else if constexpr (signType == 2) { // LS--
            fRegistry.fill(HIST("Pair/PV/b2e_b2e/lsmm/hs"), v12.M(), v12.Pt(), pair_dca);
          }
          break;
        case static_cast<int>(EM_HFeeType::kBCe_BCe):
          if constexpr (signType == 0) { // ULS
            fRegistry.fill(HIST("Pair/PV/b2c2e_b2c2e/uls/hs"), v12.M(), v12.Pt(), pair_dca);
          } else if constexpr (signType == 1) { // LS++
            fRegistry.fill(HIST("Pair/PV/b2c2e_b2c2e/lspp/hs"), v12.M(), v12.Pt(), pair_dca);
          } else if constexpr (signType == 2) { // LS--
            fRegistry.fill(HIST("Pair/PV/b2c2e_b2c2e/lsmm/hs"), v12.M(), v12.Pt(), pair_dca);
          }
          break;
        case static_cast<int>(EM_HFeeType::kBCe_Be_SameB):
          if constexpr (signType == 0) { // ULS
            fRegistry.fill(HIST("Pair/PV/b2c2e_b2e_sameb/uls/hs"), v12.M(), v12.Pt(), pair_dca);
          } else if constexpr (signType == 1) { // LS++
            fRegistry.fill(HIST("Pair/PV/b2c2e_b2e_sameb/lspp/hs"), v12.M(), v12.Pt(), pair_dca);
          } else if constexpr (signType == 2) { // LS--
            fRegistry.fill(HIST("Pair/PV/b2c2e_b2e_sameb/lsmm/hs"), v12.M(), v12.Pt(), pair_dca);
          }
          break;
        case static_cast<int>(EM_HFeeType::kBCe_Be_DiffB):
          if constexpr (signType == 0) { // ULS
            fRegistry.fill(HIST("Pair/PV/b2c2e_b2e_diffb/uls/hs"), v12.M(), v12.Pt(), pair_dca);
          } else if constexpr (signType == 1) { // LS+diff
            fRegistry.fill(HIST("Pair/PV/b2c2e_b2e_diffb/lspp/hs"), v12.M(), v12.Pt(), pair_dca);
          } else if constexpr (signType == 2) { // LS-diff
            fRegistry.fill(HIST("Pair/PV/b2c2e_b2e_diffb/lsmm/hs"), v12.M(), v12.Pt(), pair_dca);
          }
          break;
        default:
          break;
      } // end of switch for HFee
    } // end of HFee
  }

  template <uint8_t signType, typename TCollision, typename TTrack, typename TMCParticles>
  void runSVFinder(TCollision const& collision, TTrack const& t1, TTrack const& t2, TMCParticles const& mcParticles)
  {
    auto trackParCov1 = getTrackParCov(t1);
    trackParCov1.setPID(o2::track::PID::Electron);
    mDcaInfoCov.set(999, 999, 999, 999, 999);
    trackParCov1.setPID(o2::track::PID::Electron);
    o2::base::Propagator::Instance()->propagateToDCABxByBz(mVtx, trackParCov1, 2.f, matCorr, &mDcaInfoCov);
    float dcaXY1 = mDcaInfoCov.getY();
    float dcaZ1 = mDcaInfoCov.getZ();
    float CYY1 = trackParCov1.getSigmaY2();
    float CZY1 = trackParCov1.getSigmaZY();
    float CZZ1 = trackParCov1.getSigmaZ2();
    float signed1Pt1 = trackParCov1.getQ2Pt(); // at PV
    float eta1 = trackParCov1.getEta();        // at PV

    auto trackParCov2 = getTrackParCov(t2);
    trackParCov2.setPID(o2::track::PID::Electron);
    mDcaInfoCov.set(999, 999, 999, 999, 999);
    trackParCov2.setPID(o2::track::PID::Electron);
    o2::base::Propagator::Instance()->propagateToDCABxByBz(mVtx, trackParCov2, 2.f, matCorr, &mDcaInfoCov);
    float dcaXY2 = mDcaInfoCov.getY();
    float dcaZ2 = mDcaInfoCov.getZ();
    float CYY2 = trackParCov2.getSigmaY2();
    float CZY2 = trackParCov2.getSigmaZY();
    float CZZ2 = trackParCov2.getSigmaZ2();
    float signed1Pt2 = trackParCov2.getQ2Pt(); // at PV
    float eta2 = trackParCov2.getEta();        // at PV

    std::array<float, 3> pVtx = {collision.posX(), collision.posY(), collision.posZ()};
    std::array<float, 3> svpos = {0.}; // secondary vertex position
    std::array<float, 3> pvec0 = {0.};
    std::array<float, 3> pvec1 = {0.};

    int nCand = 0;
    try {
      nCand = fitter.process(trackParCov1, trackParCov2);
    } catch (...) {
      LOG(error) << "Exception caught in DCA fitter process call!";
      return;
    }
    if (nCand == 0) {
      return;
    }

    fitter.propagateTracksToVertex(); // propagate ee to SV
    if (!fitter.isPropagateTracksToVertexDone()) {
      return;
    }
    const auto& vtx = fitter.getPCACandidate();
    for (int i = 0; i < 3; i++) {
      svpos[i] = vtx[i];
    }
    fitter.getTrack(0).getPxPyPzGlo(pvec0);
    fitter.getTrack(1).getPxPyPzGlo(pvec1);
    std::array<float, 3> pvecSum = {pvec0[0] + pvec1[0], pvec0[1] + pvec1[1], pvec0[2] + pvec1[2]};

    float cpa = RecoDecay::cpa(pVtx, svpos, pvecSum);
    float cpaXY = RecoDecay::cpaXY(pVtx, svpos, pvecSum);
    float cpaRZ = RecoDecay::cpaRZ(pVtx, svpos, pvecSum);
    float chi2PCA = fitter.getChi2AtPCACandidate();
    if (svCut.cfg_max_chi2PCA < chi2PCA) {
      return;
    }

    float lxy = std::sqrt(std::pow(svpos[0] - collision.posX(), 2) + std::pow(svpos[1] - collision.posY(), 2)); // in cm
    float lz = svpos[2] - collision.posZ();                                                                     // in cm
    float lxyz = std::sqrt(std::pow(lxy, 2) + std::pow(lz, 2));

    auto primaryVertex = getPrimaryVertex(collision);
    std::array<float, 6> covVtx = fitter.calcPCACovMatrixFlat();
    double phi{}, theta{};
    getPointDirection(std::array{primaryVertex.getX(), primaryVertex.getY(), primaryVertex.getZ()}, svpos, phi, theta);
    float lxyErr = std::sqrt(getRotatedCovMatrixXX(primaryVertex.getCov(), phi, 0.) + getRotatedCovMatrixXX(covVtx, phi, 0.));
    float lzErr = std::sqrt(getRotatedCovMatrixXX(primaryVertex.getCov(), 0, theta) + getRotatedCovMatrixXX(covVtx, 0, theta));
    float lxyzErr = std::sqrt(getRotatedCovMatrixXX(primaryVertex.getCov(), phi, theta) + getRotatedCovMatrixXX(covVtx, phi, theta));

    float meeAtSV = RecoDecay::m(std::array{pvec0, pvec1}, std::array{o2::constants::physics::MassElectron, o2::constants::physics::MassElectron});
    float pteeAtSV = RecoDecay::sqrtSumOfSquares(pvecSum[0], pvecSum[1]);
    float yeeAtSV = RecoDecay::y(pvecSum, meeAtSV);

    auto t1mc = t1.template mcParticle_as<aod::McParticles>(); // true lepton
    auto t2mc = t2.template mcParticle_as<aod::McParticles>(); // true lepton
    bool isCorrectCollision1 = t1mc.mcCollisionId() == collision.mcCollisionId();
    bool isCorrectCollision2 = t2mc.mcCollisionId() == collision.mcCollisionId();
    int pdgCodeMother1 = 0;
    int pdgCodeMother2 = 0;

    int mcCommonMotherId = FindCommonMother(t1mc, t2mc, mcParticles);
    int hfee_type = IsHF(t1mc, t2mc, mcParticles);
    bool keepSignal = false;
    uint8_t dileptonType = 0;

    if (mcCommonMotherId > -1) {
      auto cmp = mcParticles.rawIteratorAt(mcCommonMotherId);
      pdgCodeMother1 = cmp.pdgCode();
      pdgCodeMother2 = cmp.pdgCode();
      switch (std::abs(cmp.pdgCode())) {
        case 23:
          keepSignal = true;
          dileptonType = 1;
          if constexpr (signType == 0) { // ULS
            fRegistry.fill(HIST("Pair/SV/Zboson/uls/hs"), meeAtSV, pteeAtSV);
            fRegistry.fill(HIST("Pair/SV/Zboson/uls/hCosPA"), cpa);
            fRegistry.fill(HIST("Pair/SV/Zboson/uls/hChi2PCA"), std::log10(chi2PCA));
          } else if constexpr (signType == 1) { // LS++
            fRegistry.fill(HIST("Pair/SV/Zboson/lspp/hs"), meeAtSV, pteeAtSV);
            fRegistry.fill(HIST("Pair/SV/Zboson/lspp/hCosPA"), cpa);
            fRegistry.fill(HIST("Pair/SV/Zboson/lspp/hChi2PCA"), std::log10(chi2PCA));
          } else if constexpr (signType == 2) { // LS--
            fRegistry.fill(HIST("Pair/SV/Zboson/lsmm/hs"), meeAtSV, pteeAtSV);
            fRegistry.fill(HIST("Pair/SV/Zboson/lsmm/hCosPA"), cpa);
            fRegistry.fill(HIST("Pair/SV/Zboson/lsmm/hChi2PCA"), std::log10(chi2PCA));
          }
          break;
        case 113: // rho0
          keepSignal = true;
          if (IsFromCharm(cmp, mcParticles) < 0 && IsFromBeauty(cmp, mcParticles) < 0) { // prompt
            dileptonType = 1;
          } else { // nonprompt
            dileptonType = 2;
          }
          break;
        case 223: // omega
          keepSignal = true;
          if (IsFromCharm(cmp, mcParticles) < 0 && IsFromBeauty(cmp, mcParticles) < 0) { // prompt
            dileptonType = 1;
          } else { // nonprompt
            dileptonType = 2;
          }
          break;
        case 333: // phi
          keepSignal = true;
          if (IsFromCharm(cmp, mcParticles) < 0 && IsFromBeauty(cmp, mcParticles) < 0) { // prompt
            dileptonType = 1;
          } else { // nonprompt
            dileptonType = 2;
          }
          break;
        case 443:
          keepSignal = true;
          if (IsFromCharm(cmp, mcParticles) < 0 && IsFromBeauty(cmp, mcParticles) < 0) { // prompt
            dileptonType = 1;
            if constexpr (signType == 0) { // ULS
              fRegistry.fill(HIST("Pair/SV/PromptJpsi/uls/hs"), meeAtSV, pteeAtSV);
              fRegistry.fill(HIST("Pair/SV/PromptJpsi/uls/hCosPA"), cpa);
              fRegistry.fill(HIST("Pair/SV/PromptJpsi/uls/hChi2PCA"), std::log10(chi2PCA));
            } else if constexpr (signType == 1) { // LS++
              fRegistry.fill(HIST("Pair/SV/PromptJpsi/lspp/hs"), meeAtSV, pteeAtSV);
              fRegistry.fill(HIST("Pair/SV/PromptJpsi/lspp/hCosPA"), cpa);
              fRegistry.fill(HIST("Pair/SV/PromptJpsi/lspp/hChi2PCA"), std::log10(chi2PCA));
            } else if constexpr (signType == 2) { // LS--
              fRegistry.fill(HIST("Pair/SV/PromptJpsi/lsmm/hs"), meeAtSV, pteeAtSV);
              fRegistry.fill(HIST("Pair/SV/PromptJpsi/lsmm/hCosPA"), cpa);
              fRegistry.fill(HIST("Pair/SV/PromptJpsi/lsmm/hChi2PCA"), std::log10(chi2PCA));
            }
          } else { // nonprompt
            dileptonType = 2;
            if constexpr (signType == 0) { // ULS
              fRegistry.fill(HIST("Pair/SV/NonPromptJpsi/uls/hs"), meeAtSV, pteeAtSV);
              fRegistry.fill(HIST("Pair/SV/NonPromptJpsi/uls/hCosPA"), cpa);
              fRegistry.fill(HIST("Pair/SV/NonPromptJpsi/uls/hChi2PCA"), std::log10(chi2PCA));
            } else if constexpr (signType == 1) { // LS++
              fRegistry.fill(HIST("Pair/SV/NonPromptJpsi/lspp/hs"), meeAtSV, pteeAtSV);
              fRegistry.fill(HIST("Pair/SV/NonPromptJpsi/lspp/hCosPA"), cpa);
              fRegistry.fill(HIST("Pair/SV/NonPromptJpsi/lspp/hChi2PCA"), std::log10(chi2PCA));
            } else if constexpr (signType == 2) { // LS--
              fRegistry.fill(HIST("Pair/SV/NonPromptJpsi/lsmm/hs"), meeAtSV, pteeAtSV);
              fRegistry.fill(HIST("Pair/SV/NonPromptJpsi/lsmm/hCosPA"), cpa);
              fRegistry.fill(HIST("Pair/SV/NonPromptJpsi/lsmm/hChi2PCA"), std::log10(chi2PCA));
            }
          }
          break;
        default:
          keepSignal = false;
          break;
      } // end of switch for LF
    } else if (hfee_type > -1) {
      keepSignal = true;
      auto t1mcMother = t1mc.template mothers_first_as<aod::McParticles>();
      auto t2mcMother = t2mc.template mothers_first_as<aod::McParticles>();
      pdgCodeMother1 = t1mcMother.pdgCode();
      pdgCodeMother2 = t2mcMother.pdgCode();
      switch (hfee_type) {
        case static_cast<int>(EM_HFeeType::kCe_Ce):
          dileptonType = 3;
          if constexpr (signType == 0) { // ULS
            fRegistry.fill(HIST("Pair/SV/c2e_c2e/uls/hs"), meeAtSV, pteeAtSV);
            fRegistry.fill(HIST("Pair/SV/c2e_c2e/uls/hCosPA"), cpa);
            fRegistry.fill(HIST("Pair/SV/c2e_c2e/uls/hChi2PCA"), std::log10(chi2PCA));
          } else if constexpr (signType == 1) { // LS++
            fRegistry.fill(HIST("Pair/SV/c2e_c2e/lspp/hs"), meeAtSV, pteeAtSV);
            fRegistry.fill(HIST("Pair/SV/c2e_c2e/lspp/hCosPA"), cpa);
            fRegistry.fill(HIST("Pair/SV/c2e_c2e/lspp/hChi2PCA"), std::log10(chi2PCA));
          } else if constexpr (signType == 2) { // LS--
            fRegistry.fill(HIST("Pair/SV/c2e_c2e/lsmm/hs"), meeAtSV, pteeAtSV);
            fRegistry.fill(HIST("Pair/SV/c2e_c2e/lsmm/hCosPA"), cpa);
            fRegistry.fill(HIST("Pair/SV/c2e_c2e/lsmm/hChi2PCA"), std::log10(chi2PCA));
          }
          break;
        case static_cast<int>(EM_HFeeType::kBe_Be):
          dileptonType = 4;
          if constexpr (signType == 0) { // ULS
            fRegistry.fill(HIST("Pair/SV/b2e_b2e/uls/hs"), meeAtSV, pteeAtSV);
            fRegistry.fill(HIST("Pair/SV/b2e_b2e/uls/hCosPA"), cpa);
            fRegistry.fill(HIST("Pair/SV/b2e_b2e/uls/hChi2PCA"), std::log10(chi2PCA));
          } else if constexpr (signType == 1) { // LS++
            fRegistry.fill(HIST("Pair/SV/b2e_b2e/lspp/hs"), meeAtSV, pteeAtSV);
            fRegistry.fill(HIST("Pair/SV/b2e_b2e/lspp/hCosPA"), cpa);
            fRegistry.fill(HIST("Pair/SV/b2e_b2e/lspp/hChi2PCA"), std::log10(chi2PCA));
          } else if constexpr (signType == 2) { // LS--
            fRegistry.fill(HIST("Pair/SV/b2e_b2e/lsmm/hs"), meeAtSV, pteeAtSV);
            fRegistry.fill(HIST("Pair/SV/b2e_b2e/lsmm/hCosPA"), cpa);
            fRegistry.fill(HIST("Pair/SV/b2e_b2e/lsmm/hChi2PCA"), std::log10(chi2PCA));
          }
          break;
        case static_cast<int>(EM_HFeeType::kBCe_BCe):
          dileptonType = 5;
          if constexpr (signType == 0) { // ULS
            fRegistry.fill(HIST("Pair/SV/b2c2e_b2c2e/uls/hs"), meeAtSV, pteeAtSV);
            fRegistry.fill(HIST("Pair/SV/b2c2e_b2c2e/uls/hCosPA"), cpa);
            fRegistry.fill(HIST("Pair/SV/b2c2e_b2c2e/uls/hChi2PCA"), std::log10(chi2PCA));
          } else if constexpr (signType == 1) { // LS++
            fRegistry.fill(HIST("Pair/SV/b2c2e_b2c2e/lspp/hs"), meeAtSV, pteeAtSV);
            fRegistry.fill(HIST("Pair/SV/b2c2e_b2c2e/lspp/hCosPA"), cpa);
            fRegistry.fill(HIST("Pair/SV/b2c2e_b2c2e/lspp/hChi2PCA"), std::log10(chi2PCA));
          } else if constexpr (signType == 2) { // LS--
            fRegistry.fill(HIST("Pair/SV/b2c2e_b2c2e/lsmm/hs"), meeAtSV, pteeAtSV);
            fRegistry.fill(HIST("Pair/SV/b2c2e_b2c2e/lsmm/hCosPA"), cpa);
            fRegistry.fill(HIST("Pair/SV/b2c2e_b2c2e/lsmm/hChi2PCA"), std::log10(chi2PCA));
          }
          break;
        case static_cast<int>(EM_HFeeType::kBCe_Be_SameB):
          dileptonType = 6;
          if constexpr (signType == 0) { // ULS
            fRegistry.fill(HIST("Pair/SV/b2c2e_b2e_sameb/uls/hs"), meeAtSV, pteeAtSV);
            fRegistry.fill(HIST("Pair/SV/b2c2e_b2e_sameb/uls/hCosPA"), cpa);
            fRegistry.fill(HIST("Pair/SV/b2c2e_b2e_sameb/uls/hChi2PCA"), std::log10(chi2PCA));
          } else if constexpr (signType == 1) { // LS++
            fRegistry.fill(HIST("Pair/SV/b2c2e_b2e_sameb/lspp/hs"), meeAtSV, pteeAtSV);
            fRegistry.fill(HIST("Pair/SV/b2c2e_b2e_sameb/lspp/hCosPA"), cpa);
            fRegistry.fill(HIST("Pair/SV/b2c2e_b2e_sameb/lspp/hChi2PCA"), std::log10(chi2PCA));
          } else if constexpr (signType == 2) { // LS--
            fRegistry.fill(HIST("Pair/SV/b2c2e_b2e_sameb/lsmm/hs"), meeAtSV, pteeAtSV);
            fRegistry.fill(HIST("Pair/SV/b2c2e_b2e_sameb/lsmm/hCosPA"), cpa);
            fRegistry.fill(HIST("Pair/SV/b2c2e_b2e_sameb/lsmm/hChi2PCA"), std::log10(chi2PCA));
          }
          break;
        case static_cast<int>(EM_HFeeType::kBCe_Be_DiffB):
          dileptonType = 7;
          if constexpr (signType == 0) { // ULS
            fRegistry.fill(HIST("Pair/SV/b2c2e_b2e_diffb/uls/hs"), meeAtSV, pteeAtSV);
            fRegistry.fill(HIST("Pair/SV/b2c2e_b2e_diffb/uls/hCosPA"), cpa);
            fRegistry.fill(HIST("Pair/SV/b2c2e_b2e_diffb/uls/hChi2PCA"), std::log10(chi2PCA));
          } else if constexpr (signType == 1) { // LS+diff
            fRegistry.fill(HIST("Pair/SV/b2c2e_b2e_diffb/lspp/hs"), meeAtSV, pteeAtSV);
            fRegistry.fill(HIST("Pair/SV/b2c2e_b2e_diffb/lspp/hCosPA"), cpa);
            fRegistry.fill(HIST("Pair/SV/b2c2e_b2e_diffb/lspp/hChi2PCA"), std::log10(chi2PCA));
          } else if constexpr (signType == 2) { // LS-diff
            fRegistry.fill(HIST("Pair/SV/b2c2e_b2e_diffb/lsmm/hs"), meeAtSV, pteeAtSV);
            fRegistry.fill(HIST("Pair/SV/b2c2e_b2e_diffb/lsmm/hCosPA"), cpa);
            fRegistry.fill(HIST("Pair/SV/b2c2e_b2e_diffb/lsmm/hChi2PCA"), std::log10(chi2PCA));
          }
          break;
        default:
          break;
      } // end of switch for HFee
    } else {
      keepSignal = true;
      dileptonType = 0; // bkg
      auto t1mcMother = t1mc.template mothers_first_as<aod::McParticles>();
      auto t2mcMother = t2mc.template mothers_first_as<aod::McParticles>();
      pdgCodeMother1 = t1mcMother.pdgCode();
      pdgCodeMother2 = t2mcMother.pdgCode();
    }

    if (keepSignal) {
      if (dileptonType == 0 && dist01(engine) > cfgDownSampling) { // random sampling, if necessary
        return;
      }

      dileptonTable(eventTable.lastIndex() + 1,
                    signed1Pt1, eta1, dcaXY1, dcaZ1, CYY1, CZY1, CZZ1, isCorrectCollision1, pdgCodeMother1,
                    signed1Pt2, eta2, dcaXY2, dcaZ2, CYY2, CZY2, CZZ2, isCorrectCollision2, pdgCodeMother2,
                    meeAtSV, pteeAtSV, yeeAtSV,
                    chi2PCA,
                    cpa, cpaXY, cpaRZ,
                    lxy, lz, lxyz,
                    lxyErr, lzErr, lxyzErr,
                    dileptonType);
    }
  }

  template <typename TBCs, typename TCollisions, typename TTracks, typename TTrackAssoc, typename TMCCollisions, typename TMCParticles>
  void run(TBCs const&, TCollisions const& collisions, TTracks const& tracks, TTrackAssoc const& trackIndices, TMCCollisions const&, TMCParticles const& mcParticles)
  {
    for (const auto& collision : collisions) {
      auto bc = collision.template bc_as<TBCs>();
      initCCDB(bc);
      fRegistry.fill(HIST("Event/hCollisionCounter"), 0);

      if (!isSelectedCollision(collision)) {
        continue;
      }

      float centralities[3] = {collision.centFT0M(), collision.centFT0A(), collision.centFT0C()};
      if (centralities[eventCut.cfgCentEstimator] < eventCut.cfgCentMin || eventCut.cfgCentMax < centralities[eventCut.cfgCentEstimator]) {
        continue;
      }
      fRegistry.fill(HIST("Event/hCollisionCounter"), 1);
      fillEventHistograms(collision);

      auto trackIdsThisCollision = trackIndices.sliceBy(trackIndicesPerCollision, collision.globalIndex());
      electronIds.reserve(trackIdsThisCollision.size());
      positronIds.reserve(trackIdsThisCollision.size());
      mVtx.setPos({collision.posX(), collision.posY(), collision.posZ()});
      mVtx.setCov(collision.covXX(), collision.covXY(), collision.covYY(), collision.covXZ(), collision.covYZ(), collision.covZZ());

      for (const auto& trackId : trackIdsThisCollision) {
        auto track = trackId.template track_as<MyTracks>();
        if (!track.hasITS() || !track.hasTPC()) {
          continue;
        }

        if (!track.has_mcParticle()) {
          continue;
        }
        auto mctrack = track.template mcParticle_as<aod::McParticles>();
        auto mcCollision = mctrack.template mcCollision_as<aod::McCollisions>();
        if (eventCut.cfgEventGeneratorType >= 0 && mcCollision.getSubGeneratorId() != eventCut.cfgEventGeneratorType) {
          continue;
        }
        if (!mctrack.has_mothers()) {
          continue;
        }
        if (std::abs(mctrack.pdgCode()) != 11) {
          continue;
        }

        auto mcMother = mctrack.template mothers_first_as<aod::McParticles>();
        if (std::abs(mcMother.pdgCode()) > 1e+9) {
          continue;
        }
        if (!(mctrack.isPhysicalPrimary() || mctrack.producedByGenerator())) {
          continue;
        }

        mDcaInfoCov.set(999, 999, 999, 999, 999);
        auto trackParCov = getTrackParCov(track);
        trackParCov.setPID(o2::track::PID::Electron);
        o2::base::Propagator::Instance()->propagateToDCABxByBz(mVtx, trackParCov, 2.f, matCorr, &mDcaInfoCov);
        float dcaXY = mDcaInfoCov.getY();
        float dcaZ = mDcaInfoCov.getZ();

        if (isSelectedTrack(track, trackParCov, dcaXY, dcaZ)) {
          if (track.sign() > 0) { // positron
            positronIds.emplace_back(track.globalIndex());
          } else { // electron
            electronIds.emplace_back(track.globalIndex());
          }
        }
      } // end of track loop for electron selection

      for (const auto& posId : positronIds) {
        auto pos = tracks.rawIteratorAt(posId);
        for (const auto& eleId : electronIds) {
          auto ele = tracks.rawIteratorAt(eleId);
          runSVFinder<0>(collision, pos, ele, mcParticles);
          runPairingAtPV<0>(pos, ele, mcParticles);
        } // end of electron loop
      } // end of positron loop

      for (const auto& posId1 : positronIds) {
        auto pos1 = tracks.rawIteratorAt(posId1);
        for (const auto& posId2 : positronIds) {
          auto pos2 = tracks.rawIteratorAt(posId2);
          if (pos1.globalIndex() == pos2.globalIndex()) {
            continue;
          }
          runSVFinder<1>(collision, pos1, pos2, mcParticles);
          runPairingAtPV<1>(pos1, pos2, mcParticles);
        } // end of positron loop
      } // end of positron loop

      for (const auto& eleId1 : electronIds) {
        auto ele1 = tracks.rawIteratorAt(eleId1);
        for (const auto& eleId2 : electronIds) {
          auto ele2 = tracks.rawIteratorAt(eleId2);
          if (ele1.globalIndex() == ele2.globalIndex()) {
            continue;
          }
          runSVFinder<2>(collision, ele1, ele2, mcParticles);
          runPairingAtPV<2>(ele1, ele2, mcParticles);
        } // end of electron loop
      } // end of electron loop

      if (positronIds.size() + electronIds.size() > 1) { // fill event table if at least 1 dilepton exists.
        eventTable(collision.numContrib(), collision.trackOccupancyInTimeRange(), collision.ft0cOccupancyInTimeRange(), 0);
      }

      electronIds.clear();
      electronIds.shrink_to_fit();
      positronIds.clear();
      positronIds.shrink_to_fit();

    } // end of collision loop
  }

  SliceCache cache;
  Preslice<aod::TracksIU> perCol = o2::aod::track::collisionId;

  Filter collisionFilter_evsel = eventCut.cfgZvtxMin < o2::aod::collision::posZ && o2::aod::collision::posZ < eventCut.cfgZvtxMax;
  Filter collisionFilter_centrality = (eventCut.cfgCentMin < o2::aod::cent::centFT0M && o2::aod::cent::centFT0M < eventCut.cfgCentMax) || (eventCut.cfgCentMin < o2::aod::cent::centFT0A && o2::aod::cent::centFT0A < eventCut.cfgCentMax) || (eventCut.cfgCentMin < o2::aod::cent::centFT0C && o2::aod::cent::centFT0C < eventCut.cfgCentMax);
  using FilteredMyCollisions = soa::Filtered<MyCollisions>;

  Preslice<aod::TrackAssoc> trackIndicesPerCollision = aod::track_association::collisionId;
  // Partition<MyFilteredTracks> posTracks = o2::aod::track::signed1Pt > 0.f;
  // Partition<MyFilteredTracks> negTracks = o2::aod::track::signed1Pt < 0.f;

  std::vector<int> electronIds;
  std::vector<int> positronIds;

  void processMC(FilteredMyCollisions const& collisions, MyBCs const& bcs, MyTracks const& tracks, aod::TrackAssoc const& trackIndices, aod::McCollisions const& mcCollisions, aod::McParticles const& mcParticles)
  {
    run(bcs, collisions, tracks, trackIndices, mcCollisions, mcParticles);
  }
  PROCESS_SWITCH(studyDCAFitter, processMC, "processMC", true);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<studyDCAFitter>(cfgc, TaskName{"study-dcafitter"})};
}
