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

/// \brief write relevant information about primary electrons.
/// \author daiki.sekihata@cern.ch

#include <unordered_map>
#include <string>
#include <vector>
#include <utility>

#include "Math/Vector4D.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "CCDB/BasicCCDBManager.h"
#include "Common/Core/trackUtilities.h"
#include "CommonConstants/PhysicsConstants.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/Core/TableHelper.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/PIDResponseITS.h"

#include "PWGEM/Dilepton/DataModel/dileptonTables.h"
#include "PWGEM/Dilepton/Utils/PairUtilities.h"

using namespace o2;
using namespace o2::soa;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;

using MyCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::EMEvSels>;
using MyCollisionsWithSWT = soa::Join<MyCollisions, aod::EMSWTriggerInfosTMP>;

using MyTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TracksCov,
                           aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr,
                           aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFbeta>;
using MyTrack = MyTracks::iterator;
using MyTracksMC = soa::Join<MyTracks, aod::McTrackLabels>;
using MyTrackMC = MyTracksMC::iterator;

struct skimmerPrimaryElectron {
  SliceCache cache;
  Preslice<aod::Tracks> perCol = o2::aod::track::collisionId;
  Produces<aod::EMPrimaryElectrons> emprimaryelectrons;
  Produces<aod::EMPrimaryElectronsCov> emprimaryelectronscov;

  // Configurables
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<bool> skipGRPOquery{"skipGRPOquery", true, "skip grpo query"};

  // Operation and minimisation criteria
  Configurable<bool> fillQAHistogram{"fillQAHistogram", false, "flag to fill QA histograms"};
  Configurable<float> d_bz_input{"d_bz_input", -999, "bz field in kG, -999 is automatic"};
  Configurable<int> min_ncluster_tpc{"min_ncluster_tpc", 10, "min ncluster tpc"};
  Configurable<int> mincrossedrows{"mincrossedrows", 70, "min. crossed rows"};
  Configurable<float> min_tpc_cr_findable_ratio{"min_tpc_cr_findable_ratio", 0.8, "min. TPC Ncr/Nf ratio"};
  Configurable<int> min_ncluster_its{"min_ncluster_its", 4, "min ncluster its"};
  Configurable<int> min_ncluster_itsib{"min_ncluster_itsib", 1, "min ncluster itsib"};
  Configurable<float> maxchi2tpc{"maxchi2tpc", 5.0, "max. chi2/NclsTPC"};
  Configurable<float> maxchi2its{"maxchi2its", 6.0, "max. chi2/NclsITS"};
  Configurable<float> minpt{"minpt", 0.15, "min pt for track"};
  Configurable<float> maxeta{"maxeta", 0.9, "eta acceptance"};
  Configurable<float> dca_xy_max{"dca_xy_max", 1.0f, "max DCAxy in cm"};
  Configurable<float> dca_z_max{"dca_z_max", 1.0f, "max DCAz in cm"};
  Configurable<float> dca_3d_sigma_max{"dca_3d_sigma_max", 1e+10, "max DCA 3D in sigma"};
  Configurable<float> minTPCNsigmaEl{"minTPCNsigmaEl", -2.5, "min. TPC n sigma for electron inclusion"};
  Configurable<float> maxTPCNsigmaEl{"maxTPCNsigmaEl", 3.5, "max. TPC n sigma for electron inclusion"};
  Configurable<float> maxTOFNsigmaEl{"maxTOFNsigmaEl", 3.5, "max. TOF n sigma for electron inclusion"};
  Configurable<float> maxTPCNsigmaPi{"maxTPCNsigmaPi", 2.5, "max. TPC n sigma for pion exclusion"};
  Configurable<float> minTPCNsigmaPi{"minTPCNsigmaPi", -1e+10, "min. TPC n sigma for pion exclusion"}; // set to -2 for lowB, -1e+10 for nominalB
  Configurable<float> maxTPCNsigmaKa{"maxTPCNsigmaKa", 2.5, "max. TPC n sigma for kaon exclusion"};
  Configurable<float> minTPCNsigmaKa{"minTPCNsigmaKa", -2.5, "min. TPC n sigma for kaon exclusion"};
  Configurable<float> maxTPCNsigmaPr{"maxTPCNsigmaPr", 2.5, "max. TPC n sigma for proton exclusion"};
  Configurable<float> minTPCNsigmaPr{"minTPCNsigmaPr", -2.5, "min. TPC n sigma for proton exclusion"};
  Configurable<bool> requireTOF{"requireTOF", false, "require TOF hit"};
  Configurable<float> max_pin_for_pion_rejection{"max_pin_for_pion_rejection", 1e+10, "pion rejection is applied below this pin"};
  Configurable<float> max_frac_shared_clusters_tpc{"max_frac_shared_clusters_tpc", 999.f, "max fraction of shared clusters in TPC"};

  HistogramRegistry fRegistry{"output", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};

  int mRunNumber;
  float d_bz;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;

  void init(InitContext&)
  {
    mRunNumber = 0;
    d_bz = 0;

    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);

    if (fillQAHistogram) {
      fRegistry.add("Track/hPt", "pT;p_{T} (GeV/c)", kTH1F, {{1000, 0.0f, 10}}, false);
      fRegistry.add("Track/hQoverPt", "q/pT;q/p_{T} (GeV/c)^{-1}", kTH1F, {{400, -20, 20}}, false);
      fRegistry.add("Track/hEtaPhi", "#eta vs. #varphi;#varphi (rad.);#eta", kTH2F, {{180, 0, 2 * M_PI}, {20, -1.0f, 1.0f}}, false);
      fRegistry.add("Track/hDCAxyz", "DCA xy vs. z;DCA_{xy} (cm);DCA_{z} (cm)", kTH2F, {{200, -1.0f, 1.0f}, {200, -1.0f, 1.0f}}, false);
      fRegistry.add("Track/hDCAxyzSigma", "DCA xy vs. z;DCA_{xy} (#sigma);DCA_{z} (#sigma)", kTH2F, {{200, -10.0f, 10.0f}, {200, -10.0f, 10.0f}}, false);
      fRegistry.add("Track/hDCAxyRes_Pt", "DCA_{xy} resolution vs. pT;p_{T} (GeV/c);DCA_{xy} resolution (#mum)", kTH2F, {{1000, 0, 10}, {500, 0., 500}}, false);
      fRegistry.add("Track/hDCAzRes_Pt", "DCA_{z} resolution vs. pT;p_{T} (GeV/c);DCA_{z} resolution (#mum)", kTH2F, {{1000, 0, 10}, {500, 0., 500}}, false);
      fRegistry.add("Track/hNclsTPC", "number of TPC clusters", kTH1F, {{161, -0.5, 160.5}}, false);
      fRegistry.add("Track/hNcrTPC", "number of TPC crossed rows", kTH1F, {{161, -0.5, 160.5}}, false);
      fRegistry.add("Track/hChi2TPC", "chi2/number of TPC clusters", kTH1F, {{100, 0, 10}}, false);
      fRegistry.add("Track/hChi2TOF", "chi2 of TOF", kTH1F, {{100, 0, 10}}, false);
      fRegistry.add("Track/hTPCNcr2Nf", "TPC Ncr/Nfindable", kTH1F, {{200, 0, 2}}, false);
      fRegistry.add("Track/hTPCNcls2Nf", "TPC Ncls/Nfindable", kTH1F, {{200, 0, 2}}, false);
      fRegistry.add("Track/hTPCNclsShared", "TPC Ncls shared/Ncls;p_{T} (GeV/c);N_{cls}^{shared}/N_{cls} in TPC", kTH2F, {{1000, 0, 10}, {100, 0, 1}}, false);
      fRegistry.add("Track/hNclsITS", "number of ITS clusters", kTH1F, {{8, -0.5, 7.5}}, false);
      fRegistry.add("Track/hChi2ITS", "chi2/number of ITS clusters", kTH1F, {{100, 0, 10}}, false);
      fRegistry.add("Track/hITSClusterMap", "ITS cluster map", kTH1F, {{128, -0.5, 127.5}}, false);
      fRegistry.add("Track/hTPCdEdx", "TPC dE/dx;p_{in} (GeV/c);TPC dE/dx (a.u.)", kTH2F, {{1000, 0, 10}, {200, 0, 200}}, false);
      fRegistry.add("Track/hTPCNsigmaEl", "TPC n sigma el;p_{in} (GeV/c);n #sigma_{e}^{TPC}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
      fRegistry.add("Track/hTPCNsigmaMu", "TPC n sigma mu;p_{in} (GeV/c);n #sigma_{#mu}^{TPC}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
      fRegistry.add("Track/hTPCNsigmaPi", "TPC n sigma pi;p_{in} (GeV/c);n #sigma_{#pi}^{TPC}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
      fRegistry.add("Track/hTPCNsigmaKa", "TPC n sigma ka;p_{in} (GeV/c);n #sigma_{K}^{TPC}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
      fRegistry.add("Track/hTPCNsigmaPr", "TPC n sigma pr;p_{in} (GeV/c);n #sigma_{p}^{TPC}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
      fRegistry.add("Track/hTOFbeta", "TOF beta;p_{pv} (GeV/c);#beta", kTH2F, {{1000, 0, 10}, {240, 0, 1.2}}, false);
      fRegistry.add("Track/hTOFNsigmaEl", "TOF n sigma el;p_{in} (GeV/c);n #sigma_{e}^{TOF}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
      fRegistry.add("Track/hTOFNsigmaMu", "TOF n sigma mu;p_{in} (GeV/c);n #sigma_{#mu}^{TOF}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
      fRegistry.add("Track/hTOFNsigmaPi", "TOF n sigma pi;p_{in} (GeV/c);n #sigma_{#pi}^{TOF}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
      fRegistry.add("Track/hTOFNsigmaKa", "TOF n sigma ka;p_{in} (GeV/c);n #sigma_{K}^{TOF}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
      fRegistry.add("Track/hTOFNsigmaPr", "TOF n sigma pr;p_{in} (GeV/c);n #sigma_{p}^{TOF}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
      fRegistry.add("Track/hMeanClusterSizeITS", "mean cluster size ITS;p_{pv} (GeV/c);<cluster size> on ITS #times cos(#lambda)", kTH2F, {{1000, 0, 10}, {150, 0, 15}}, false);
      fRegistry.add("Track/hITSNsigmaEl", "ITS n sigma el;p_{pv} (GeV/c);n #sigma_{e}^{ITS}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
      fRegistry.add("Track/hITSNsigmaMu", "ITS n sigma mu;p_{pv} (GeV/c);n #sigma_{#mu}^{ITS}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
      fRegistry.add("Track/hITSNsigmaPi", "ITS n sigma pi;p_{pv} (GeV/c);n #sigma_{#pi}^{ITS}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
      fRegistry.add("Track/hITSNsigmaKa", "ITS n sigma ka;p_{pv} (GeV/c);n #sigma_{K}^{ITS}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
      fRegistry.add("Track/hITSNsigmaPr", "ITS n sigma pr;p_{pv} (GeV/c);n #sigma_{p}^{ITS}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    }
  }

  void initCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }

    // In case override, don't proceed, please - no CCDB access required
    if (d_bz_input > -990) {
      d_bz = d_bz_input;
      o2::parameters::GRPMagField grpmag;
      if (std::fabs(d_bz) > 1e-5) {
        grpmag.setL3Current(30000.f / (d_bz / 5.0f));
      }
      o2::base::Propagator::initFieldFromGRP(&grpmag);
      mRunNumber = bc.runNumber();
      return;
    }

    auto run3grp_timestamp = bc.timestamp();
    o2::parameters::GRPObject* grpo = 0x0;
    o2::parameters::GRPMagField* grpmag = 0x0;
    if (!skipGRPOquery)
      grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>(grpPath, run3grp_timestamp);
    if (grpo) {
      o2::base::Propagator::initFieldFromGRP(grpo);
      // Fetch magnetic field from ccdb for current collision
      d_bz = grpo->getNominalL3Field();
      LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << d_bz << " kZG";
    } else {
      grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, run3grp_timestamp);
      if (!grpmag) {
        LOG(fatal) << "Got nullptr from CCDB for path " << grpmagPath << " of object GRPMagField and " << grpPath << " of object GRPObject for timestamp " << run3grp_timestamp;
      }
      o2::base::Propagator::initFieldFromGRP(grpmag);
      // Fetch magnetic field from ccdb for current collision
      d_bz = std::lround(5.f * grpmag->getL3Current() / 30000.f);
      LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << d_bz << " kZG";
    }
    mRunNumber = bc.runNumber();
  }

  template <bool isMC, typename TCollision, typename TTrack>
  bool checkTrack(TCollision const& collision, TTrack const& track)
  {
    if constexpr (isMC) {
      if (!track.has_mcParticle()) {
        return false;
      }
    }

    if (requireTOF && !(track.hasTOF() && std::fabs(track.tofNSigmaEl()) < maxTOFNsigmaEl)) {
      return false;
    }

    if (track.tpcChi2NCl() > maxchi2tpc) {
      return false;
    }

    if (track.itsChi2NCl() > maxchi2its) {
      return false;
    }

    if (!track.hasITS() || !track.hasTPC()) {
      return false;
    }
    if (track.itsNCls() < min_ncluster_its) {
      return false;
    }
    if (track.itsNClsInnerBarrel() < min_ncluster_itsib) {
      return false;
    }

    if (track.tpcNClsFound() < min_ncluster_tpc) {
      return false;
    }

    if (track.tpcNClsCrossedRows() < mincrossedrows) {
      return false;
    }

    if (track.tpcCrossedRowsOverFindableCls() < min_tpc_cr_findable_ratio) {
      return false;
    }

    if (track.tpcFractionSharedCls() > max_frac_shared_clusters_tpc) {
      return false;
    }

    if (track.hasTOF() && (maxTOFNsigmaEl < std::fabs(track.tofNSigmaEl()))) {
      return false;
    }

    gpu::gpustd::array<float, 2> dcaInfo;
    auto track_par_cov_recalc = getTrackParCov(track);
    track_par_cov_recalc.setPID(o2::track::PID::Electron);
    // std::array<float, 3> pVec_recalc = {0, 0, 0}; // px, py, pz
    o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, track_par_cov_recalc, 2.f, matCorr, &dcaInfo);
    // getPxPyPz(track_par_cov_recalc, pVec_recalc);
    float dcaXY = dcaInfo[0];
    float dcaZ = dcaInfo[1];

    if (std::fabs(dcaXY) > dca_xy_max || std::fabs(dcaZ) > dca_z_max) {
      return false;
    }

    if (track_par_cov_recalc.getPt() < minpt || std::fabs(track_par_cov_recalc.getEta()) > maxeta) {
      return false;
    }

    float dca_3d = 999.f;
    float det = track_par_cov_recalc.getSigmaY2() * track_par_cov_recalc.getSigmaZ2() - track_par_cov_recalc.getSigmaZY() * track_par_cov_recalc.getSigmaZY();
    if (det < 0) {
      dca_3d = 999.f;
    } else {
      float chi2 = (dcaXY * dcaXY * track_par_cov_recalc.getSigmaZ2() + dcaZ * dcaZ * track_par_cov_recalc.getSigmaY2() - 2. * dcaXY * dcaZ * track_par_cov_recalc.getSigmaZY()) / det;
      dca_3d = std::sqrt(std::fabs(chi2) / 2.);
    }
    if (dca_3d > dca_3d_sigma_max) {
      return false;
    }

    return true;
  }

  template <typename TTrack>
  bool isElectron(TTrack const& track)
  {
    return isElectron_TPChadrej(track) || isElectron_TOFreq(track);
  }

  template <typename TTrack>
  bool isElectron_TPChadrej(TTrack const& track)
  {
    if (track.tpcNSigmaEl() < minTPCNsigmaEl || maxTPCNsigmaEl < track.tpcNSigmaEl()) {
      return false;
    }
    if (minTPCNsigmaPi < track.tpcNSigmaPi() && track.tpcNSigmaPi() < maxTPCNsigmaPi && track.tpcInnerParam() < max_pin_for_pion_rejection) {
      return false;
    }
    if (minTPCNsigmaKa < track.tpcNSigmaKa() && track.tpcNSigmaKa() < maxTPCNsigmaKa) {
      return false;
    }
    if (minTPCNsigmaPr < track.tpcNSigmaPr() && track.tpcNSigmaPr() < maxTPCNsigmaPr) {
      return false;
    }
    if (track.hasTOF() && (maxTOFNsigmaEl < std::fabs(track.tofNSigmaEl()))) {
      return false;
    }
    return true;
  }

  template <typename TTrack>
  bool isElectron_TOFreq(TTrack const& track)
  {
    if (minTPCNsigmaPi < track.tpcNSigmaPi() && track.tpcNSigmaPi() < maxTPCNsigmaPi && track.tpcInnerParam() < max_pin_for_pion_rejection) {
      return false;
    }
    return minTPCNsigmaEl < track.tpcNSigmaEl() && track.tpcNSigmaEl() < maxTPCNsigmaEl && std::fabs(track.tofNSigmaEl()) < maxTOFNsigmaEl;
  }

  template <typename TCollision, typename TTrack>
  void fillTrackTable(TCollision const& collision, TTrack const& track)
  {
    if (std::find(stored_trackIds.begin(), stored_trackIds.end(), std::pair<int, int>{collision.globalIndex(), track.globalIndex()}) == stored_trackIds.end()) {
      gpu::gpustd::array<float, 2> dcaInfo;
      auto track_par_cov_recalc = getTrackParCov(track);
      track_par_cov_recalc.setPID(o2::track::PID::Electron);
      // std::array<float, 3> pVec_recalc = {0, 0, 0}; // px, py, pz
      o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, track_par_cov_recalc, 2.f, matCorr, &dcaInfo);
      // getPxPyPz(track_par_cov_recalc, pVec_recalc);
      float dcaXY = dcaInfo[0];
      float dcaZ = dcaInfo[1];

      float pt_recalc = track_par_cov_recalc.getPt();
      float eta_recalc = track_par_cov_recalc.getEta();
      float phi_recalc = track_par_cov_recalc.getPhi();

      bool isAssociatedToMPC = collision.globalIndex() == track.collisionId();

      emprimaryelectrons(collision.globalIndex(), track.globalIndex(), track.sign(),
                         pt_recalc, eta_recalc, phi_recalc, dcaXY, dcaZ,
                         track.tpcNClsFindable(), track.tpcNClsFindableMinusFound(), track.tpcNClsFindableMinusCrossedRows(), track.tpcNClsShared(),
                         track.tpcChi2NCl(), track.tpcInnerParam(),
                         track.tpcSignal(), track.tpcNSigmaEl(), track.tpcNSigmaMu(), track.tpcNSigmaPi(), track.tpcNSigmaKa(), track.tpcNSigmaPr(),
                         track.beta(), track.tofNSigmaEl(), track.tofNSigmaMu(), track.tofNSigmaPi(), track.tofNSigmaKa(), track.tofNSigmaPr(),
                         track.itsClusterSizes(), track.itsNSigmaEl(), track.itsNSigmaMu(), track.itsNSigmaPi(), track.itsNSigmaKa(), track.itsNSigmaPr(),
                         track.itsChi2NCl(), track.tofChi2(), track.detectorMap(),
                         track_par_cov_recalc.getX(), track_par_cov_recalc.getAlpha(), track_par_cov_recalc.getY(), track_par_cov_recalc.getZ(), track_par_cov_recalc.getSnp(), track_par_cov_recalc.getTgl(), isAssociatedToMPC);

      emprimaryelectronscov(
        track_par_cov_recalc.getSigmaY2(),
        track_par_cov_recalc.getSigmaZY(),
        track_par_cov_recalc.getSigmaZ2(),
        track_par_cov_recalc.getSigmaSnpY(),
        track_par_cov_recalc.getSigmaSnpZ(),
        track_par_cov_recalc.getSigmaSnp2(),
        track_par_cov_recalc.getSigmaTglY(),
        track_par_cov_recalc.getSigmaTglZ(),
        track_par_cov_recalc.getSigmaTglSnp(),
        track_par_cov_recalc.getSigmaTgl2(),
        track_par_cov_recalc.getSigma1PtY(),
        track_par_cov_recalc.getSigma1PtZ(),
        track_par_cov_recalc.getSigma1PtSnp(),
        track_par_cov_recalc.getSigma1PtTgl(),
        track_par_cov_recalc.getSigma1Pt2());

      stored_trackIds.emplace_back(std::pair<int, int>{collision.globalIndex(), track.globalIndex()});

      if (fillQAHistogram) {
        uint32_t itsClusterSizes = track.itsClusterSizes();
        int total_cluster_size = 0, nl = 0;
        for (unsigned int layer = 0; layer < 7; layer++) {
          int cluster_size_per_layer = (itsClusterSizes >> (layer * 4)) & 0xf;
          if (cluster_size_per_layer > 0) {
            nl++;
          }
          total_cluster_size += cluster_size_per_layer;
        }

        fRegistry.fill(HIST("Track/hPt"), pt_recalc);
        fRegistry.fill(HIST("Track/hQoverPt"), track.sign() / pt_recalc);
        fRegistry.fill(HIST("Track/hEtaPhi"), phi_recalc, eta_recalc);
        fRegistry.fill(HIST("Track/hDCAxyz"), dcaXY, dcaZ);
        fRegistry.fill(HIST("Track/hDCAxyzSigma"), dcaXY / sqrt(track_par_cov_recalc.getSigmaY2()), dcaZ / sqrt(track_par_cov_recalc.getSigmaZ2()));
        fRegistry.fill(HIST("Track/hDCAxyRes_Pt"), pt_recalc, sqrt(track_par_cov_recalc.getSigmaY2()) * 1e+4); // convert cm to um
        fRegistry.fill(HIST("Track/hDCAzRes_Pt"), pt_recalc, sqrt(track_par_cov_recalc.getSigmaZ2()) * 1e+4);  // convert cm to um
        fRegistry.fill(HIST("Track/hNclsITS"), track.itsNCls());
        fRegistry.fill(HIST("Track/hNclsTPC"), track.tpcNClsFound());
        fRegistry.fill(HIST("Track/hNcrTPC"), track.tpcNClsCrossedRows());
        fRegistry.fill(HIST("Track/hTPCNcr2Nf"), track.tpcCrossedRowsOverFindableCls());
        fRegistry.fill(HIST("Track/hTPCNcls2Nf"), track.tpcFoundOverFindableCls());
        fRegistry.fill(HIST("Track/hTPCNclsShared"), track.pt(), track.tpcFractionSharedCls());
        fRegistry.fill(HIST("Track/hChi2TPC"), track.tpcChi2NCl());
        fRegistry.fill(HIST("Track/hChi2ITS"), track.itsChi2NCl());
        fRegistry.fill(HIST("Track/hChi2TOF"), track.tofChi2());
        fRegistry.fill(HIST("Track/hITSClusterMap"), track.itsClusterMap());
        fRegistry.fill(HIST("Track/hTPCdEdx"), track.tpcInnerParam(), track.tpcSignal());
        fRegistry.fill(HIST("Track/hTPCNsigmaEl"), track.tpcInnerParam(), track.tpcNSigmaEl());
        fRegistry.fill(HIST("Track/hTPCNsigmaMu"), track.tpcInnerParam(), track.tpcNSigmaMu());
        fRegistry.fill(HIST("Track/hTPCNsigmaPi"), track.tpcInnerParam(), track.tpcNSigmaPi());
        fRegistry.fill(HIST("Track/hTPCNsigmaKa"), track.tpcInnerParam(), track.tpcNSigmaKa());
        fRegistry.fill(HIST("Track/hTPCNsigmaPr"), track.tpcInnerParam(), track.tpcNSigmaPr());
        fRegistry.fill(HIST("Track/hTOFbeta"), track.p(), track.beta());
        fRegistry.fill(HIST("Track/hTOFNsigmaEl"), track.tpcInnerParam(), track.tofNSigmaEl());
        fRegistry.fill(HIST("Track/hTOFNsigmaMu"), track.tpcInnerParam(), track.tofNSigmaMu());
        fRegistry.fill(HIST("Track/hTOFNsigmaPi"), track.tpcInnerParam(), track.tofNSigmaPi());
        fRegistry.fill(HIST("Track/hTOFNsigmaKa"), track.tpcInnerParam(), track.tofNSigmaKa());
        fRegistry.fill(HIST("Track/hTOFNsigmaPr"), track.tpcInnerParam(), track.tofNSigmaPr());
        fRegistry.fill(HIST("Track/hMeanClusterSizeITS"), track.p(), static_cast<float>(total_cluster_size) / static_cast<float>(nl) * std::cos(std::atan(track.tgl())));
        fRegistry.fill(HIST("Track/hITSNsigmaEl"), track.p(), track.itsNSigmaEl());
        fRegistry.fill(HIST("Track/hITSNsigmaMu"), track.p(), track.itsNSigmaMu());
        fRegistry.fill(HIST("Track/hITSNsigmaPi"), track.p(), track.itsNSigmaPi());
        fRegistry.fill(HIST("Track/hITSNsigmaKa"), track.p(), track.itsNSigmaKa());
        fRegistry.fill(HIST("Track/hITSNsigmaPr"), track.p(), track.itsNSigmaPr());
      }
    }
  }

  Preslice<aod::TrackAssoc> trackIndicesPerCollision = aod::track_association::collisionId;
  std::vector<std::pair<int, int>> stored_trackIds;
  Filter trackFilter = o2::aod::track::pt > minpt&& nabs(o2::aod::track::eta) < maxeta&& o2::aod::track::tpcChi2NCl < maxchi2tpc&& o2::aod::track::itsChi2NCl < maxchi2its&& ncheckbit(aod::track::v001::detectorMap, (uint8_t)o2::aod::track::ITS) == true && ncheckbit(aod::track::v001::detectorMap, (uint8_t)o2::aod::track::TPC) == true;
  Filter pidFilter = minTPCNsigmaEl < o2::aod::pidtpc::tpcNSigmaEl && o2::aod::pidtpc::tpcNSigmaEl < maxTPCNsigmaEl;
  using MyFilteredTracks = soa::Filtered<MyTracks>;

  Partition<MyFilteredTracks> posTracks = o2::aod::track::signed1Pt > 0.f;
  Partition<MyFilteredTracks> negTracks = o2::aod::track::signed1Pt < 0.f;

  // ---------- for data ----------

  void processRec_SA(MyCollisions const& collisions, aod::BCsWithTimestamps const&, MyFilteredTracks const& tracks)
  {
    auto tracksWithITSPid = soa::Attach<MyFilteredTracks, aod::pidits::ITSNSigmaEl, aod::pidits::ITSNSigmaMu, aod::pidits::ITSNSigmaPi, aod::pidits::ITSNSigmaKa, aod::pidits::ITSNSigmaPr>(tracks);
    stored_trackIds.reserve(tracks.size());

    for (auto& collision : collisions) {
      auto bc = collision.template foundBC_as<aod::BCsWithTimestamps>();
      initCCDB(bc);

      if (!collision.isSelected()) {
        continue;
      }

      auto tracks_per_coll = tracksWithITSPid.sliceBy(perCol, collision.globalIndex());
      for (auto& track : tracks_per_coll) {
        if (!checkTrack<false>(collision, track) || !isElectron(track)) {
          continue;
        }
        fillTrackTable(collision, track);
      }

    } // end of collision loop

    stored_trackIds.clear();
    stored_trackIds.shrink_to_fit();
  }
  PROCESS_SWITCH(skimmerPrimaryElectron, processRec_SA, "process reconstructed info only", true); // standalone

  void processRec_TTCA(MyCollisions const& collisions, aod::BCsWithTimestamps const&, MyTracks const& tracks, aod::TrackAssoc const& trackIndices)
  {
    auto tracksWithITSPid = soa::Attach<MyTracks, aod::pidits::ITSNSigmaEl, aod::pidits::ITSNSigmaMu, aod::pidits::ITSNSigmaPi, aod::pidits::ITSNSigmaKa, aod::pidits::ITSNSigmaPr>(tracks);
    stored_trackIds.reserve(tracks.size() * 2);

    for (auto& collision : collisions) {
      auto bc = collision.template foundBC_as<aod::BCsWithTimestamps>();
      initCCDB(bc);

      if (!collision.isSelected()) {
        continue;
      }

      auto trackIdsThisCollision = trackIndices.sliceBy(trackIndicesPerCollision, collision.globalIndex());

      for (auto& trackId : trackIdsThisCollision) {
        // auto track = trackId.template track_as<MyTracks>();
        auto track = tracksWithITSPid.rawIteratorAt(trackId.trackId());
        if (!checkTrack<false>(collision, track) || !isElectron(track)) {
          continue;
        }
        fillTrackTable(collision, track);
      }
    } // end of collision loop

    stored_trackIds.clear();
    stored_trackIds.shrink_to_fit();
  }
  PROCESS_SWITCH(skimmerPrimaryElectron, processRec_TTCA, "process reconstructed info only", false); // with TTCA

  void processRec_SA_SWT(MyCollisionsWithSWT const& collisions, aod::BCsWithTimestamps const&, MyFilteredTracks const& tracks)
  {
    auto tracksWithITSPid = soa::Attach<MyFilteredTracks, aod::pidits::ITSNSigmaEl, aod::pidits::ITSNSigmaMu, aod::pidits::ITSNSigmaPi, aod::pidits::ITSNSigmaKa, aod::pidits::ITSNSigmaPr>(tracks);
    stored_trackIds.reserve(tracks.size());

    for (auto& collision : collisions) {
      auto bc = collision.template foundBC_as<aod::BCsWithTimestamps>();
      initCCDB(bc);

      if (!collision.isSelected()) {
        continue;
      }

      if (collision.swtaliastmp_raw() == 0) {
        continue;
      }

      auto tracks_per_coll = tracksWithITSPid.sliceBy(perCol, collision.globalIndex());
      for (auto& track : tracks_per_coll) {
        if (!checkTrack<false>(collision, track) || !isElectron(track)) {
          continue;
        }
        fillTrackTable(collision, track);
      }

    } // end of collision loop

    stored_trackIds.clear();
    stored_trackIds.shrink_to_fit();
  }
  PROCESS_SWITCH(skimmerPrimaryElectron, processRec_SA_SWT, "process reconstructed info only", false); // standalone with swt

  void processRec_TTCA_SWT(MyCollisionsWithSWT const& collisions, aod::BCsWithTimestamps const&, MyTracks const& tracks, aod::TrackAssoc const& trackIndices)
  {
    auto tracksWithITSPid = soa::Attach<MyTracks, aod::pidits::ITSNSigmaEl, aod::pidits::ITSNSigmaMu, aod::pidits::ITSNSigmaPi, aod::pidits::ITSNSigmaKa, aod::pidits::ITSNSigmaPr>(tracks);
    stored_trackIds.reserve(tracks.size() * 2);

    for (auto& collision : collisions) {
      auto bc = collision.template foundBC_as<aod::BCsWithTimestamps>();
      initCCDB(bc);

      if (!collision.isSelected()) {
        continue;
      }
      if (collision.swtaliastmp_raw() == 0) {
        continue;
      }

      auto trackIdsThisCollision = trackIndices.sliceBy(trackIndicesPerCollision, collision.globalIndex());

      for (auto& trackId : trackIdsThisCollision) {
        // auto track = trackId.template track_as<MyTracks>();
        auto track = tracksWithITSPid.rawIteratorAt(trackId.trackId());
        if (!checkTrack<false>(collision, track) || !isElectron(track)) {
          continue;
        }
        fillTrackTable(collision, track);
      }
    } // end of collision loop

    stored_trackIds.clear();
    stored_trackIds.shrink_to_fit();
  }
  PROCESS_SWITCH(skimmerPrimaryElectron, processRec_TTCA_SWT, "process reconstructed info only", false); // with TTCA with swt

  // ---------- for MC ----------

  using MyFilteredTracksMC = soa::Filtered<MyTracksMC>;
  Partition<MyFilteredTracksMC> posTracksMC = o2::aod::track::signed1Pt > 0.f;
  Partition<MyFilteredTracksMC> negTracksMC = o2::aod::track::signed1Pt < 0.f;
  void processMC_SA(soa::Join<MyCollisions, aod::McCollisionLabels> const& collisions, aod::McCollisions const&, aod::BCsWithTimestamps const&, MyFilteredTracksMC const& tracks)
  {
    auto tracksWithITSPid = soa::Attach<MyFilteredTracksMC, aod::pidits::ITSNSigmaEl, aod::pidits::ITSNSigmaMu, aod::pidits::ITSNSigmaPi, aod::pidits::ITSNSigmaKa, aod::pidits::ITSNSigmaPr>(tracks);
    stored_trackIds.reserve(tracks.size());

    for (auto& collision : collisions) {
      if (!collision.has_mcCollision()) {
        continue;
      }
      auto bc = collision.template foundBC_as<aod::BCsWithTimestamps>();
      initCCDB(bc);

      if (!collision.isSelected()) {
        continue;
      }

      auto tracks_per_coll = tracksWithITSPid.sliceBy(perCol, collision.globalIndex());
      for (auto& track : tracks_per_coll) {
        if (!checkTrack<true>(collision, track) || !isElectron(track)) {
          continue;
        }
        fillTrackTable(collision, track);
      }
    } // end of collision loop

    stored_trackIds.clear();
    stored_trackIds.shrink_to_fit();
  }
  PROCESS_SWITCH(skimmerPrimaryElectron, processMC_SA, "process reconstructed and MC info ", false);

  void processMC_TTCA(soa::Join<MyCollisions, aod::McCollisionLabels> const& collisions, aod::McCollisions const&, aod::BCsWithTimestamps const&, MyTracksMC const& tracks, aod::TrackAssoc const& trackIndices)
  {
    auto tracksWithITSPid = soa::Attach<MyTracksMC, aod::pidits::ITSNSigmaEl, aod::pidits::ITSNSigmaMu, aod::pidits::ITSNSigmaPi, aod::pidits::ITSNSigmaKa, aod::pidits::ITSNSigmaPr>(tracks);
    stored_trackIds.reserve(tracks.size() * 2);

    for (auto& collision : collisions) {
      if (!collision.has_mcCollision()) {
        continue;
      }
      auto bc = collision.template foundBC_as<aod::BCsWithTimestamps>();
      initCCDB(bc);

      if (!collision.isSelected()) {
        continue;
      }

      auto trackIdsThisCollision = trackIndices.sliceBy(trackIndicesPerCollision, collision.globalIndex());

      for (auto& trackId : trackIdsThisCollision) {
        // auto track = trackId.template track_as<MyTracksMC>();
        auto track = tracksWithITSPid.rawIteratorAt(trackId.trackId());
        if (!checkTrack<true>(collision, track) || !isElectron(track)) {
          continue;
        }
        fillTrackTable(collision, track);
      }
    } // end of collision loop

    stored_trackIds.clear();
    stored_trackIds.shrink_to_fit();
  }
  PROCESS_SWITCH(skimmerPrimaryElectron, processMC_TTCA, "process reconstructed info only", false); // with TTCA
};

struct prefilterPrimaryElectron {
  Produces<aod::EMPrimaryElectronsPrefilterBit> ele_pfb;

  SliceCache cache;
  Preslice<aod::Tracks> perCol_track = o2::aod::track::collisionId;
  PresliceUnsorted<aod::EMPrimaryElectrons> perCol_ele = o2::aod::emprimaryelectron::collisionId;

  // CCDB options
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<bool> skipGRPOquery{"skipGRPOquery", true, "skip grpo query"};

  // Operation and minimisation criteria
  Configurable<double> d_bz_input{"d_bz", -999, "bz field, -999 is automatic"};

  Configurable<float> max_dcaxy{"max_dcaxy", 0.3, "DCAxy To PV for loose track sample"};
  Configurable<float> max_dcaz{"max_dcaz", 0.3, "DCAz To PV for loose track sample"};
  Configurable<float> minpt{"minpt", 0.1, "min pt for track for loose track sample"};
  Configurable<float> maxeta{"maxeta", 0.9, "eta acceptance for loose track sample"};
  Configurable<int> min_ncluster_tpc{"min_ncluster_tpc", 0, "min ncluster tpc"};
  Configurable<int> mincrossedrows{"mincrossedrows", 70, "min crossed rows"};
  Configurable<float> max_frac_shared_clusters_tpc{"max_frac_shared_clusters_tpc", 999.f, "max fraction of shared clusters in TPC"};
  Configurable<float> min_tpc_cr_findable_ratio{"min_tpc_cr_findable_ratio", 0.8, "min. TPC Ncr/Nf ratio"};
  Configurable<float> maxchi2tpc{"maxchi2tpc", 5.0, "max chi2/NclsTPC"};
  Configurable<float> maxchi2its{"maxchi2its", 6.0, "max chi2/NclsITS"};
  Configurable<int> min_ncluster_its{"min_ncluster_its", 4, "min ncluster its"};
  Configurable<int> min_ncluster_itsib{"min_ncluster_itsib", 1, "min ncluster itsib"};
  Configurable<float> minTPCNsigmaEl{"minTPCNsigmaEl", -3.0, "min. TPC n sigma for electron inclusion"};
  Configurable<float> maxTPCNsigmaEl{"maxTPCNsigmaEl", 3.0, "max. TPC n sigma for electron inclusion"};
  Configurable<float> slope{"slope", 0.0185, "slope for m vs. phiv"};
  Configurable<float> intercept{"intercept", -0.0280, "intercept for m vs. phiv"};

  Configurable<std::vector<float>> max_mee_vec{"max_mee_vec", std::vector<float>{0.08, 0.10, 0.12}, "vector fo max mee for prefilter in ULS. Please sort this by increasing order."}; // currently, 3 thoresholds are allowed.

  HistogramRegistry fRegistry{"output", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};

  int mRunNumber;
  float d_bz;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;

  void init(InitContext&)
  {
    mRunNumber = 0;
    d_bz = 0;

    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);

    if (!doprocessDummy) {
      addHistograms();
    }
  }

  void addHistograms()
  {
    fRegistry.add("Track/hPt", "pT;p_{T} (GeV/c)", kTH1F, {{1000, 0.0f, 10}}, false);
    fRegistry.add("Track/hEtaPhi", "#eta vs. #varphi;#varphi (rad.);#eta", kTH2F, {{180, 0, 2 * M_PI}, {40, -1.0f, 1.0f}}, false);
    fRegistry.add("Track/hTPCNsigmaEl", "loose track TPC PID", kTH2F, {{1000, 0.f, 10}, {100, -5, +5}});
    fRegistry.add("Pair/before/uls/hMvsPt", "mass vs. pT;m_{ee} (GeV/c^{2});p_{T,ee} (GeV/c)", kTH2F, {{400, 0, 4}, {100, 0, 10}});
    fRegistry.add("Pair/before/uls/hMvsPhiV", "mass vs. phiv;#varphi_{V} (rad.);m_{ee} (GeV/c^{2})", kTH2F, {{90, 0.f, M_PI}, {100, 0, 1.f}});
    fRegistry.addClone("Pair/before/uls/", "Pair/before/lspp/");
    fRegistry.addClone("Pair/before/uls/", "Pair/before/lsmm/");
    fRegistry.addClone("Pair/before/", "Pair/after/");
  }

  void initCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }

    // In case override, don't proceed, please - no CCDB access required
    if (d_bz_input > -990) {
      d_bz = d_bz_input;
      o2::parameters::GRPMagField grpmag;
      if (std::fabs(d_bz) > 1e-5) {
        grpmag.setL3Current(30000.f / (d_bz / 5.0f));
      }
      o2::base::Propagator::initFieldFromGRP(&grpmag);
      mRunNumber = bc.runNumber();
      return;
    }

    auto run3grp_timestamp = bc.timestamp();
    o2::parameters::GRPObject* grpo = 0x0;
    o2::parameters::GRPMagField* grpmag = 0x0;
    if (!skipGRPOquery)
      grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>(grpPath, run3grp_timestamp);
    if (grpo) {
      o2::base::Propagator::initFieldFromGRP(grpo);
      // Fetch magnetic field from ccdb for current collision
      d_bz = grpo->getNominalL3Field();
      LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << d_bz << " kZG";
    } else {
      grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, run3grp_timestamp);
      if (!grpmag) {
        LOG(fatal) << "Got nullptr from CCDB for path " << grpmagPath << " of object GRPMagField and " << grpPath << " of object GRPObject for timestamp " << run3grp_timestamp;
      }
      o2::base::Propagator::initFieldFromGRP(grpmag);
      // Fetch magnetic field from ccdb for current collision
      d_bz = std::lround(5.f * grpmag->getL3Current() / 30000.f);
      LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << d_bz << " kZG";
    }
    mRunNumber = bc.runNumber();
  }

  o2::base::Propagator::MatCorrType noMatCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;

  template <typename TCollision, typename TTrack>
  bool checkTrack(TCollision const& collision, TTrack const& track)
  {
    if (!track.hasITS()) {
      return false;
    }
    if (track.itsChi2NCl() > maxchi2its) {
      return false;
    }
    if (track.itsNCls() < min_ncluster_its) {
      return false;
    }
    if (track.itsNClsInnerBarrel() < min_ncluster_itsib) {
      return false;
    }

    if (!track.hasTPC()) {
      return false;
    }
    if (track.tpcNSigmaEl() < minTPCNsigmaEl || maxTPCNsigmaEl < track.tpcNSigmaEl()) {
      return false;
    }
    if (track.tpcNClsFound() < min_ncluster_tpc) {
      return false;
    }
    if (track.tpcNClsCrossedRows() < mincrossedrows) {
      return false;
    }
    if (track.tpcCrossedRowsOverFindableCls() < min_tpc_cr_findable_ratio) {
      return false;
    }
    if (track.tpcFractionSharedCls() > max_frac_shared_clusters_tpc) {
      return false;
    }
    if (track.tpcChi2NCl() > maxchi2its) {
      return false;
    }

    gpu::gpustd::array<float, 2> dcaInfo;
    auto track_par_cov_recalc = getTrackParCov(track);
    track_par_cov_recalc.setPID(o2::track::PID::Electron);
    // std::array<float, 3> pVec_recalc = {0, 0, 0}; // px, py, pz
    o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, track_par_cov_recalc, 2.f, matCorr, &dcaInfo);
    // getPxPyPz(track_par_cov_recalc, pVec_recalc);

    if (std::fabs(dcaInfo[0]) > max_dcaxy || std::fabs(dcaInfo[1]) > max_dcaz) {
      return false;
    }

    if (track_par_cov_recalc.getPt() < minpt || std::fabs(track_par_cov_recalc.getEta()) > maxeta) {
      return false;
    }

    return true;
  }

  template <int loose_track_sign, typename TCollision, typename TTrack1, typename TTrack2>
  bool reconstructPC(TCollision const& collision, TTrack1 const& ele, TTrack2 const& pos)
  {
    float mee = 0, phiv = 0;
    gpu::gpustd::array<float, 2> dcaInfo;
    std::array<float, 3> pVec_recalc = {0, 0, 0}; // px, py, pz

    if constexpr (loose_track_sign > 0) { // positive track is loose track
      auto track_par_cov_recalc = getTrackParCov(pos);
      o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, track_par_cov_recalc, 2.f, matCorr, &dcaInfo);
      getPxPyPz(track_par_cov_recalc, pVec_recalc);

      ROOT::Math::PtEtaPhiMVector v1(ele.pt(), ele.eta(), ele.phi(), o2::constants::physics::MassElectron);
      ROOT::Math::PtEtaPhiMVector v2(track_par_cov_recalc.getPt(), track_par_cov_recalc.getEta(), track_par_cov_recalc.getPhi(), o2::constants::physics::MassElectron);
      ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
      mee = v12.M();
      phiv = o2::aod::pwgem::dilepton::utils::pairutil::getPhivPair(pVec_recalc[0], pVec_recalc[1], pVec_recalc[2], ele.px(), ele.py(), ele.pz(), pos.sign(), ele.sign(), d_bz);
    } else {
      auto track_par_cov_recalc = getTrackParCov(ele);
      o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, track_par_cov_recalc, 2.f, matCorr, &dcaInfo);
      getPxPyPz(track_par_cov_recalc, pVec_recalc);

      ROOT::Math::PtEtaPhiMVector v1(track_par_cov_recalc.getPt(), track_par_cov_recalc.getEta(), track_par_cov_recalc.getPhi(), o2::constants::physics::MassElectron);
      ROOT::Math::PtEtaPhiMVector v2(pos.pt(), pos.eta(), pos.phi(), o2::constants::physics::MassElectron);
      ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
      mee = v12.M();
      phiv = o2::aod::pwgem::dilepton::utils::pairutil::getPhivPair(pos.px(), pos.py(), pos.pz(), pVec_recalc[0], pVec_recalc[1], pVec_recalc[2], pos.sign(), ele.sign(), d_bz);
    }

    if (mee < slope * phiv + intercept) {
      return true;
    } else {
      return false;
    }
  }

  Preslice<aod::TrackAssoc> trackIndicesPerCollision = aod::track_association::collisionId;

  Filter trackFilter = o2::aod::track::pt > minpt&& nabs(o2::aod::track::eta) < maxeta&& ncheckbit(aod::track::v001::detectorMap, (uint8_t)o2::aod::track::ITS) == true && ncheckbit(aod::track::v001::detectorMap, (uint8_t)o2::aod::track::TPC) == true;
  using MyFilteredTracks = soa::Filtered<MyTracks>;
  Partition<MyFilteredTracks> posTracks = o2::aod::track::signed1Pt > 0.f;
  Partition<MyFilteredTracks> negTracks = o2::aod::track::signed1Pt < 0.f;

  Partition<aod::EMPrimaryElectrons> positrons = o2::aod::emprimaryelectron::sign > int8_t(0);
  Partition<aod::EMPrimaryElectrons> electrons = o2::aod::emprimaryelectron::sign < int8_t(0);

  void processPrefilter_TTCA(MyCollisions const& collisions, aod::BCsWithTimestamps const&, MyTracks const&, aod::EMPrimaryElectrons const& primaryelectrons, aod::TrackAssoc const& trackIndices)
  {
    std::unordered_map<int, uint8_t> pfb_map; // map track.globalIndex -> prefilter bit

    for (auto& collision : collisions) {
      auto bc = collision.template foundBC_as<aod::BCsWithTimestamps>();
      initCCDB(bc);
      if (!collision.isSelected()) {
        continue;
      }

      auto positrons_per_coll = positrons->sliceByCachedUnsorted(o2::aod::emprimaryelectron::collisionId, collision.globalIndex(), cache); // signal sample
      auto electrons_per_coll = electrons->sliceByCachedUnsorted(o2::aod::emprimaryelectron::collisionId, collision.globalIndex(), cache); // signal sample

      auto trackIdsThisCollision = trackIndices.sliceBy(trackIndicesPerCollision, collision.globalIndex());
      std::vector<MyTrack> posTracks_per_coll;
      std::vector<MyTrack> negTracks_per_coll;
      posTracks_per_coll.reserve(trackIdsThisCollision.size());
      negTracks_per_coll.reserve(trackIdsThisCollision.size());

      for (auto& trackId : trackIdsThisCollision) {
        auto track = trackId.template track_as<MyTracks>();
        if (!checkTrack(collision, track)) {
          continue;
        }
        fRegistry.fill(HIST("Track/hPt"), track.pt());
        fRegistry.fill(HIST("Track/hEtaPhi"), track.phi(), track.eta());
        if (track.sign() > 0) {
          posTracks_per_coll.emplace_back(track);
        } else {
          negTracks_per_coll.emplace_back(track);
        }
      }

      for (auto& ele : negTracks_per_coll) {
        if (!checkTrack(collision, ele)) {
          continue;
        }
        gpu::gpustd::array<float, 2> dcaInfo;
        std::array<float, 3> pVec_recalc = {0, 0, 0}; // px, py, pz
        auto track_par_cov_recalc = getTrackParCov(ele);
        track_par_cov_recalc.setPID(o2::track::PID::Electron);
        o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, track_par_cov_recalc, 2.f, matCorr, &dcaInfo);
        getPxPyPz(track_par_cov_recalc, pVec_recalc);

        for (auto& empos : positrons_per_coll) {
          if (empos.trackId() == ele.globalIndex()) {
            continue;
          }

          ROOT::Math::PtEtaPhiMVector v1(track_par_cov_recalc.getPt(), track_par_cov_recalc.getEta(), track_par_cov_recalc.getPhi(), o2::constants::physics::MassElectron); // loose track
          ROOT::Math::PtEtaPhiMVector v2(empos.pt(), empos.eta(), empos.phi(), o2::constants::physics::MassElectron);                                                       // signal track
          ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
          float phiv = o2::aod::pwgem::dilepton::utils::pairutil::getPhivPair(empos.px(), empos.py(), empos.pz(), pVec_recalc[0], pVec_recalc[1], pVec_recalc[2], empos.sign(), ele.sign(), d_bz);
          fRegistry.fill(HIST("Pair/before/uls/hMvsPhiV"), phiv, v12.M());
          fRegistry.fill(HIST("Pair/before/uls/hMvsPt"), v12.M(), v12.Pt());
          if (v12.M() < max_mee_vec->at(static_cast<int>(max_mee_vec->size()) - 1)) {
            fRegistry.fill(HIST("Track/hTPCNsigmaEl"), ele.tpcInnerParam(), ele.tpcNSigmaEl());
          }
          for (int i = 0; i < static_cast<int>(max_mee_vec->size()); i++) {
            if (v12.M() < max_mee_vec->at(i)) {
              pfb_map[empos.globalIndex()] |= (uint8_t(1) << (static_cast<int>(o2::aod::pwgem::dilepton::utils::pairutil::DileptonPrefilterBit::kElFromPi0_1) + i));
            }
          }

          if (v12.M() < slope * phiv + intercept) {
            pfb_map[empos.globalIndex()] |= (uint8_t(1) << static_cast<int>(o2::aod::pwgem::dilepton::utils::pairutil::DileptonPrefilterBit::kElFromPC));
          }

        } // end of signal positon loop
      } // end of loose electron loop

      for (auto& pos : posTracks_per_coll) {
        if (!checkTrack(collision, pos)) { // track cut is applied to loose sample
          continue;
        }
        gpu::gpustd::array<float, 2> dcaInfo;
        std::array<float, 3> pVec_recalc = {0, 0, 0}; // px, py, pz
        auto track_par_cov_recalc = getTrackParCov(pos);
        track_par_cov_recalc.setPID(o2::track::PID::Electron);
        o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, track_par_cov_recalc, 2.f, matCorr, &dcaInfo);
        getPxPyPz(track_par_cov_recalc, pVec_recalc);
        for (auto& emele : electrons_per_coll) {
          if (emele.trackId() == pos.globalIndex()) {
            continue;
          }

          ROOT::Math::PtEtaPhiMVector v1(emele.pt(), emele.eta(), emele.phi(), o2::constants::physics::MassElectron);                                                       // signal track
          ROOT::Math::PtEtaPhiMVector v2(track_par_cov_recalc.getPt(), track_par_cov_recalc.getEta(), track_par_cov_recalc.getPhi(), o2::constants::physics::MassElectron); // loose track
          ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
          float phiv = o2::aod::pwgem::dilepton::utils::pairutil::getPhivPair(pVec_recalc[0], pVec_recalc[1], pVec_recalc[2], emele.px(), emele.py(), emele.pz(), pos.sign(), emele.sign(), d_bz);
          fRegistry.fill(HIST("Pair/before/uls/hMvsPhiV"), phiv, v12.M());
          fRegistry.fill(HIST("Pair/before/uls/hMvsPt"), v12.M(), v12.Pt());
          if (v12.M() < max_mee_vec->at(static_cast<int>(max_mee_vec->size()) - 1)) {
            fRegistry.fill(HIST("Track/hTPCNsigmaEl"), pos.tpcInnerParam(), pos.tpcNSigmaEl());
          }
          for (int i = 0; i < static_cast<int>(max_mee_vec->size()); i++) {
            if (v12.M() < max_mee_vec->at(i)) {
              pfb_map[emele.globalIndex()] |= (uint8_t(1) << (static_cast<int>(o2::aod::pwgem::dilepton::utils::pairutil::DileptonPrefilterBit::kElFromPi0_1) + i));
            }
          }

          if (v12.M() < slope * phiv + intercept) {
            pfb_map[emele.globalIndex()] |= (uint8_t(1) << static_cast<int>(o2::aod::pwgem::dilepton::utils::pairutil::DileptonPrefilterBit::kElFromPC));
          }
        } // end of signal electron loop
      } // end of loose positon loop

      for (auto& pos : posTracks_per_coll) {
        if (!checkTrack(collision, pos)) { // track cut is applied to loose sample
          continue;
        }
        gpu::gpustd::array<float, 2> dcaInfo;
        std::array<float, 3> pVec_recalc = {0, 0, 0}; // px, py, pz
        auto track_par_cov_recalc = getTrackParCov(pos);
        track_par_cov_recalc.setPID(o2::track::PID::Electron);
        o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, track_par_cov_recalc, 2.f, matCorr, &dcaInfo);
        getPxPyPz(track_par_cov_recalc, pVec_recalc);
        for (auto& empos : positrons_per_coll) {
          if (empos.trackId() == pos.globalIndex()) {
            continue;
          }

          ROOT::Math::PtEtaPhiMVector v1(empos.pt(), empos.eta(), empos.phi(), o2::constants::physics::MassElectron);                                                       // signal track
          ROOT::Math::PtEtaPhiMVector v2(track_par_cov_recalc.getPt(), track_par_cov_recalc.getEta(), track_par_cov_recalc.getPhi(), o2::constants::physics::MassElectron); // loose track
          ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
          float phiv = o2::aod::pwgem::dilepton::utils::pairutil::getPhivPair(pVec_recalc[0], pVec_recalc[1], pVec_recalc[2], empos.px(), empos.py(), empos.pz(), pos.sign(), empos.sign(), d_bz);
          fRegistry.fill(HIST("Pair/before/lspp/hMvsPhiV"), phiv, v12.M());
          fRegistry.fill(HIST("Pair/before/lspp/hMvsPt"), v12.M(), v12.Pt());
        } // end of signal positron loop
      } // end of loose positon loop

      for (auto& ele : negTracks_per_coll) {
        if (!checkTrack(collision, ele)) {
          continue;
        }
        gpu::gpustd::array<float, 2> dcaInfo;
        std::array<float, 3> pVec_recalc = {0, 0, 0}; // px, py, pz
        auto track_par_cov_recalc = getTrackParCov(ele);
        track_par_cov_recalc.setPID(o2::track::PID::Electron);
        o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, track_par_cov_recalc, 2.f, matCorr, &dcaInfo);
        getPxPyPz(track_par_cov_recalc, pVec_recalc);

        for (auto& emele : electrons_per_coll) {
          if (emele.trackId() == ele.globalIndex()) {
            continue;
          }

          ROOT::Math::PtEtaPhiMVector v1(track_par_cov_recalc.getPt(), track_par_cov_recalc.getEta(), track_par_cov_recalc.getPhi(), o2::constants::physics::MassElectron); // loose track
          ROOT::Math::PtEtaPhiMVector v2(emele.pt(), emele.eta(), emele.phi(), o2::constants::physics::MassElectron);                                                       // signal track
          ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
          float phiv = o2::aod::pwgem::dilepton::utils::pairutil::getPhivPair(emele.px(), emele.py(), emele.pz(), pVec_recalc[0], pVec_recalc[1], pVec_recalc[2], emele.sign(), ele.sign(), d_bz);
          fRegistry.fill(HIST("Pair/before/lsmm/hMvsPhiV"), phiv, v12.M());
          fRegistry.fill(HIST("Pair/before/lsmm/hMvsPt"), v12.M(), v12.Pt());

        } // end of signal electron loop
      } // end of loose electron loop

      posTracks_per_coll.clear();
      negTracks_per_coll.clear();
      posTracks_per_coll.shrink_to_fit();
      negTracks_per_coll.shrink_to_fit();
    } // end of collision loop

    for (auto& ele : primaryelectrons) {
      ele_pfb(pfb_map[ele.globalIndex()]);
    }

    // check prefilter
    for (auto& collision : collisions) {
      auto positrons_per_coll = positrons->sliceByCachedUnsorted(o2::aod::emprimaryelectron::collisionId, collision.globalIndex(), cache); // signal sample
      auto electrons_per_coll = electrons->sliceByCachedUnsorted(o2::aod::emprimaryelectron::collisionId, collision.globalIndex(), cache); // signal sample

      for (auto& [ele, pos] : combinations(CombinationsFullIndexPolicy(electrons_per_coll, positrons_per_coll))) {
        if (pfb_map[ele.globalIndex()] != 0 || pfb_map[pos.globalIndex()] != 0) {
          continue;
        }

        ROOT::Math::PtEtaPhiMVector v1(ele.pt(), ele.eta(), ele.phi(), o2::constants::physics::MassElectron);
        ROOT::Math::PtEtaPhiMVector v2(pos.pt(), pos.eta(), pos.phi(), o2::constants::physics::MassElectron);
        ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
        float phiv = o2::aod::pwgem::dilepton::utils::pairutil::getPhivPair(pos.px(), pos.py(), pos.pz(), ele.px(), ele.py(), ele.pz(), pos.sign(), ele.sign(), d_bz);
        fRegistry.fill(HIST("Pair/after/uls/hMvsPhiV"), phiv, v12.M());
        fRegistry.fill(HIST("Pair/after/uls/hMvsPt"), v12.M(), v12.Pt());

      } // end of ULS pairing
    } // end of collision loop

    pfb_map.clear();
  }
  PROCESS_SWITCH(prefilterPrimaryElectron, processPrefilter_TTCA, "process prefilter with TTCA", false);

  void processPrefilter_SA(MyCollisions const& collisions, aod::BCsWithTimestamps const&, MyFilteredTracks const&, aod::EMPrimaryElectrons const& primaryelectrons)
  {
    std::unordered_map<int, uint8_t> pfb_map; // map track.globalIndex -> prefilter bit

    for (auto& collision : collisions) {
      auto bc = collision.template foundBC_as<aod::BCsWithTimestamps>();
      initCCDB(bc);
      if (!collision.isSelected()) {
        continue;
      }

      auto posTracks_per_coll = posTracks->sliceByCached(o2::aod::track::collisionId, collision.globalIndex(), cache); // loose track sample
      auto negTracks_per_coll = negTracks->sliceByCached(o2::aod::track::collisionId, collision.globalIndex(), cache); // loose track sample

      auto positrons_per_coll = positrons->sliceByCachedUnsorted(o2::aod::emprimaryelectron::collisionId, collision.globalIndex(), cache); // signal sample
      auto electrons_per_coll = electrons->sliceByCachedUnsorted(o2::aod::emprimaryelectron::collisionId, collision.globalIndex(), cache); // signal sample

      for (auto& pos : posTracks_per_coll) {
        if (!checkTrack(collision, pos)) { // track cut is applied to loose sample
          continue;
        }
        fRegistry.fill(HIST("Track/hPt"), pos.pt());
        fRegistry.fill(HIST("Track/hEtaPhi"), pos.phi(), pos.eta());
      }
      for (auto& neg : negTracks_per_coll) {
        if (!checkTrack(collision, neg)) { // track cut is applied to loose sample
          continue;
        }
        fRegistry.fill(HIST("Track/hPt"), neg.pt());
        fRegistry.fill(HIST("Track/hEtaPhi"), neg.phi(), neg.eta());
      }

      for (auto& [ele, empos] : combinations(CombinationsFullIndexPolicy(negTracks_per_coll, positrons_per_coll))) {
        // auto pos = tracks.rawIteratorAt(empos.trackId()); // use rawIterator, if the table is filtered.
        if (!checkTrack(collision, ele)) { // track cut is applied to loose sample
          continue;
        }
        if (empos.trackId() == ele.globalIndex()) {
          continue;
        }

        ROOT::Math::PtEtaPhiMVector v1(ele.pt(), ele.eta(), ele.phi(), o2::constants::physics::MassElectron);       // loose track
        ROOT::Math::PtEtaPhiMVector v2(empos.pt(), empos.eta(), empos.phi(), o2::constants::physics::MassElectron); // signal track
        ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
        float phiv = o2::aod::pwgem::dilepton::utils::pairutil::getPhivPair(empos.px(), empos.py(), empos.pz(), ele.px(), ele.py(), ele.pz(), empos.sign(), ele.sign(), d_bz);
        fRegistry.fill(HIST("Pair/before/uls/hMvsPhiV"), phiv, v12.M());
        fRegistry.fill(HIST("Pair/before/uls/hMvsPt"), v12.M(), v12.Pt());
        if (v12.M() < max_mee_vec->at(static_cast<int>(max_mee_vec->size()) - 1)) {
          fRegistry.fill(HIST("Track/hTPCNsigmaEl"), ele.tpcInnerParam(), ele.tpcNSigmaEl());
        }
        for (int i = 0; i < static_cast<int>(max_mee_vec->size()); i++) {
          if (v12.M() < max_mee_vec->at(i)) {
            pfb_map[empos.globalIndex()] |= (uint8_t(1) << (static_cast<int>(o2::aod::pwgem::dilepton::utils::pairutil::DileptonPrefilterBit::kElFromPi0_1) + i));
          }
        }

        if (v12.M() < slope * phiv + intercept) {
          pfb_map[empos.globalIndex()] |= (uint8_t(1) << static_cast<int>(o2::aod::pwgem::dilepton::utils::pairutil::DileptonPrefilterBit::kElFromPC));
        }

      } // end of ULS pairing

      for (auto& [pos, emele] : combinations(CombinationsFullIndexPolicy(posTracks_per_coll, electrons_per_coll))) {
        // auto ele = tracks.rawIteratorAt(emele.trackId()); // use rawIterator, if the table is filtered.
        if (!checkTrack(collision, pos)) { // track cut is applied to loose sample
          continue;
        }
        if (emele.trackId() == pos.globalIndex()) {
          continue;
        }

        ROOT::Math::PtEtaPhiMVector v1(emele.pt(), emele.eta(), emele.phi(), o2::constants::physics::MassElectron); // signal track
        ROOT::Math::PtEtaPhiMVector v2(pos.pt(), pos.eta(), pos.phi(), o2::constants::physics::MassElectron);       // loose track
        ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
        float phiv = o2::aod::pwgem::dilepton::utils::pairutil::getPhivPair(pos.px(), pos.py(), pos.pz(), emele.px(), emele.py(), emele.pz(), pos.sign(), emele.sign(), d_bz);
        fRegistry.fill(HIST("Pair/before/uls/hMvsPhiV"), phiv, v12.M());
        fRegistry.fill(HIST("Pair/before/uls/hMvsPt"), v12.M(), v12.Pt());
        if (v12.M() < max_mee_vec->at(static_cast<int>(max_mee_vec->size()) - 1)) {
          fRegistry.fill(HIST("Track/hTPCNsigmaEl"), pos.tpcInnerParam(), pos.tpcNSigmaEl());
        }
        for (int i = 0; i < static_cast<int>(max_mee_vec->size()); i++) {
          if (v12.M() < max_mee_vec->at(i)) {
            pfb_map[emele.globalIndex()] |= (uint8_t(1) << (static_cast<int>(o2::aod::pwgem::dilepton::utils::pairutil::DileptonPrefilterBit::kElFromPi0_1) + i));
          }
        }

        if (v12.M() < slope * phiv + intercept) {
          pfb_map[emele.globalIndex()] |= (uint8_t(1) << static_cast<int>(o2::aod::pwgem::dilepton::utils::pairutil::DileptonPrefilterBit::kElFromPC));
        }

      } // end of ULS pairing

      for (auto& [pos, empos] : combinations(CombinationsFullIndexPolicy(posTracks_per_coll, positrons_per_coll))) {
        // auto pos = tracks.rawIteratorAt(empos.trackId()); // use rawIterator, if the table is filtered.
        if (!checkTrack(collision, pos)) { // track cut is applied to loose sample
          continue;
        }
        if (empos.trackId() == pos.globalIndex()) {
          continue;
        }

        ROOT::Math::PtEtaPhiMVector v1(pos.pt(), pos.eta(), pos.phi(), o2::constants::physics::MassElectron);       // loose track
        ROOT::Math::PtEtaPhiMVector v2(empos.pt(), empos.eta(), empos.phi(), o2::constants::physics::MassElectron); // signal track
        ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
        float phiv = o2::aod::pwgem::dilepton::utils::pairutil::getPhivPair(empos.px(), empos.py(), empos.pz(), pos.px(), pos.py(), pos.pz(), empos.sign(), pos.sign(), d_bz);
        fRegistry.fill(HIST("Pair/before/lspp/hMvsPhiV"), phiv, v12.M());
        fRegistry.fill(HIST("Pair/before/lspp/hMvsPt"), v12.M(), v12.Pt());
      } // end of LS++ pairing

      for (auto& [ele, emele] : combinations(CombinationsFullIndexPolicy(negTracks_per_coll, electrons_per_coll))) {
        // auto ele = tracks.rawIteratorAt(emele.trackId()); // use rawIterator, if the table is filtered.
        if (!checkTrack(collision, ele)) { // track cut is applied to loose sample
          continue;
        }
        if (emele.trackId() == ele.globalIndex()) {
          continue;
        }

        ROOT::Math::PtEtaPhiMVector v1(ele.pt(), ele.eta(), ele.phi(), o2::constants::physics::MassElectron);       // loose track
        ROOT::Math::PtEtaPhiMVector v2(emele.pt(), emele.eta(), emele.phi(), o2::constants::physics::MassElectron); // signal track
        ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
        float phiv = o2::aod::pwgem::dilepton::utils::pairutil::getPhivPair(emele.px(), emele.py(), emele.pz(), ele.px(), ele.py(), ele.pz(), emele.sign(), ele.sign(), d_bz);
        fRegistry.fill(HIST("Pair/before/lsmm/hMvsPhiV"), phiv, v12.M());
        fRegistry.fill(HIST("Pair/before/lsmm/hMvsPt"), v12.M(), v12.Pt());
      } // end of LS-- pairing

    } // end of collision loop

    for (auto& ele : primaryelectrons) {
      ele_pfb(pfb_map[ele.globalIndex()]);
    }

    // check prefilter
    for (auto& collision : collisions) {
      auto positrons_per_coll = positrons->sliceByCachedUnsorted(o2::aod::emprimaryelectron::collisionId, collision.globalIndex(), cache); // signal sample
      auto electrons_per_coll = electrons->sliceByCachedUnsorted(o2::aod::emprimaryelectron::collisionId, collision.globalIndex(), cache); // signal sample

      for (auto& [ele, pos] : combinations(CombinationsFullIndexPolicy(electrons_per_coll, positrons_per_coll))) {
        if (pfb_map[ele.globalIndex()] != 0 || pfb_map[pos.globalIndex()] != 0) {
          continue;
        }

        ROOT::Math::PtEtaPhiMVector v1(ele.pt(), ele.eta(), ele.phi(), o2::constants::physics::MassElectron);
        ROOT::Math::PtEtaPhiMVector v2(pos.pt(), pos.eta(), pos.phi(), o2::constants::physics::MassElectron);
        ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
        float phiv = o2::aod::pwgem::dilepton::utils::pairutil::getPhivPair(pos.px(), pos.py(), pos.pz(), ele.px(), ele.py(), ele.pz(), pos.sign(), ele.sign(), d_bz);
        fRegistry.fill(HIST("Pair/after/uls/hMvsPhiV"), phiv, v12.M());
        fRegistry.fill(HIST("Pair/after/uls/hMvsPt"), v12.M(), v12.Pt());
      } // end of ULS pairing
    } // end of collision loop

    pfb_map.clear();
  }
  PROCESS_SWITCH(prefilterPrimaryElectron, processPrefilter_SA, "process prefilter standalone", false);

  void processDummy(aod::EMPrimaryElectrons const& primaryelectrons)
  {
    for (int i = 0; i < primaryelectrons.size(); i++) {
      ele_pfb(0);
    }
  }
  PROCESS_SWITCH(prefilterPrimaryElectron, processDummy, "process dummy", true);
};

struct associateAmbiguousElectron {
  Produces<aod::EMAmbiguousElectronSelfIds> em_amb_ele_ids;

  SliceCache cache;
  PresliceUnsorted<aod::EMPrimaryElectrons> perTrack = o2::aod::emprimaryelectron::trackId;
  std::vector<int> ambele_self_Ids;

  void process(aod::EMPrimaryElectrons const& electrons)
  {
    for (auto& electron : electrons) {
      auto electrons_with_same_trackId = electrons.sliceBy(perTrack, electron.trackId());
      ambele_self_Ids.reserve(electrons_with_same_trackId.size());
      for (auto& amb_ele : electrons_with_same_trackId) {
        if (amb_ele.globalIndex() == electron.globalIndex()) { // don't store myself.
          continue;
        }
        ambele_self_Ids.emplace_back(amb_ele.globalIndex());
      }
      em_amb_ele_ids(ambele_self_Ids);
      ambele_self_Ids.clear();
      ambele_self_Ids.shrink_to_fit();
    }
  }
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<skimmerPrimaryElectron>(cfgc, TaskName{"skimmer-primary-electron"}),
    adaptAnalysisTask<prefilterPrimaryElectron>(cfgc, TaskName{"prefilter-primary-electron"}),
    adaptAnalysisTask<associateAmbiguousElectron>(cfgc, TaskName{"associate-ambiguous-electron"})};
}
