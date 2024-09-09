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

#include "PWGEM/Dilepton/DataModel/dileptonTables.h"
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGEM/Dilepton/Utils/PairUtilities.h"
#include "PWGEM/Dilepton/Utils/EMTrackUtilities.h"

using namespace o2;
using namespace o2::soa;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;

using namespace o2::aod::pwgem::dilepton::utils::emtrackutil;

using MyTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TracksCov,
                           aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr,
                           aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFbeta>;
using MyTrack = MyTracks::iterator;
using MyTracksMC = soa::Join<MyTracks, aod::McTrackLabels>;
using MyTrackMC = MyTracksMC::iterator;

struct filterDielectronEvent {
  using MyCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::EMEvSels>;
  using MyCollisionsWithSWT = soa::Join<MyCollisions, aod::EMSWTriggerInfosTMP>;

  SliceCache cache;
  Preslice<aod::Tracks> perCol = o2::aod::track::collisionId;
  Produces<aod::EMPrimaryElectrons> emprimaryelectrons;
  Produces<aod::EMPrimaryElectronsCov> emprimaryelectronscov;
  Produces<aod::EMEventsNee> filter;

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
  Configurable<float> max_mean_its_cluster_size{"max_mean_its_cluster_size", 16.f, "max. <ITS cluster size> x cos(lambda)"};
  Configurable<float> max_p_for_its_cluster_size{"max_p_for_its_cluster_size", 0.2, "its cluster size cut is applied below this p"};
  Configurable<float> max_pin_for_pion_rejection{"max_pin_for_pion_rejection", -1, "pion rejection is applied below this pin"};
  Configurable<int> minitsncls{"minitsncls", 4, "min. number of ITS clusters"};
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
  Configurable<float> minTPCNsigmaPi{"minTPCNsigmaPi", -1e+10, "min. TPC n sigma for pion exclusion"};
  Configurable<float> maxTPCNsigmaPi{"maxTPCNsigmaPi", 2.0, "max. TPC n sigma for pion exclusion"};
  Configurable<float> maxMee{"maxMee", 0.02, "max mee for virtual photon selection"};

  Configurable<bool> apply_phiv{"apply_phiv", true, "flag to apply phiv cut"};
  Configurable<float> slope{"slope", 0.0181, "slope for mee vs. phiv"};
  Configurable<float> intercept{"intercept", -0.0370, "intercept for mee vs. phiv"};

  HistogramRegistry fRegistry{"output", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};

  std::pair<int8_t, std::set<uint8_t>> itsRequirement = {1, {0, 1, 2}}; // any hits on 3 ITS ib layers.

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
      fRegistry.add("Track/hRelSigma1Pt", "relative p_{T} resolution;p_{T} (GeV/c);#sigma_{1/p_{T}} #times p_{T}", kTH2F, {{1000, 0, 10}, {100, 0, 0.1}}, false);
      fRegistry.add("Track/hEtaPhi", "#eta vs. #varphi;#varphi (rad.);#eta", kTH2F, {{180, 0, 2 * M_PI}, {20, -1.0f, 1.0f}}, false);
      fRegistry.add("Track/hDCAxyz", "DCA xy vs. z;DCA_{xy} (cm);DCA_{z} (cm)", kTH2F, {{200, -1.0f, 1.0f}, {200, -1.0f, 1.0f}}, false);
      fRegistry.add("Track/hDCAxyzSigma", "DCA xy vs. z;DCA_{xy} (#sigma);DCA_{z} (#sigma)", kTH2F, {{200, -10.0f, 10.0f}, {200, -10.0f, 10.0f}}, false);
      fRegistry.add("Track/hDCAxyRes_Pt", "DCA_{xy} resolution vs. pT;p_{T} (GeV/c);DCA_{xy} resolution (#mum)", kTH2F, {{1000, 0, 10}, {500, 0., 500}}, false);
      fRegistry.add("Track/hDCAzRes_Pt", "DCA_{z} resolution vs. pT;p_{T} (GeV/c);DCA_{z} resolution (#mum)", kTH2F, {{1000, 0, 10}, {500, 0., 500}}, false);
      fRegistry.add("Track/hNclsTPC", "number of TPC clusters", kTH1F, {{161, -0.5, 160.5}}, false);
      fRegistry.add("Track/hNcrTPC", "number of TPC crossed rows", kTH1F, {{161, -0.5, 160.5}}, false);
      fRegistry.add("Track/hChi2TPC", "chi2/number of TPC clusters", kTH1F, {{100, 0, 10}}, false);
      fRegistry.add("Track/hTPCdEdx", "TPC dE/dx;p_{in} (GeV/c);TPC dE/dx (a.u.)", kTH2F, {{1000, 0, 10}, {200, 0, 200}}, false);
      fRegistry.add("Track/hTPCNsigmaEl", "TPC n sigma el;p_{in} (GeV/c);n #sigma_{e}^{TPC}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
      fRegistry.add("Track/hTPCNsigmaMu", "TPC n sigma mu;p_{in} (GeV/c);n #sigma_{#mu}^{TPC}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
      fRegistry.add("Track/hTPCNsigmaPi", "TPC n sigma pi;p_{in} (GeV/c);n #sigma_{#pi}^{TPC}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
      fRegistry.add("Track/hTPCNsigmaKa", "TPC n sigma ka;p_{in} (GeV/c);n #sigma_{K}^{TPC}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
      fRegistry.add("Track/hTPCNsigmaPr", "TPC n sigma pr;p_{in} (GeV/c);n #sigma_{p}^{TPC}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
      fRegistry.add("Track/hTOFbeta", "TOF beta;p_{pv} (GeV/c);#beta", kTH2F, {{1000, 0, 10}, {240, 0, 1.2}}, false);
      fRegistry.add("Track/hTOFNsigmaEl", "TOF n sigma el;p_{pv} (GeV/c);n #sigma_{e}^{TOF}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
      fRegistry.add("Track/hTOFNsigmaMu", "TOF n sigma mu;p_{pv} (GeV/c);n #sigma_{#mu}^{TOF}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
      fRegistry.add("Track/hTOFNsigmaPi", "TOF n sigma pi;p_{pv} (GeV/c);n #sigma_{#pi}^{TOF}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
      fRegistry.add("Track/hTOFNsigmaKa", "TOF n sigma ka;p_{pv} (GeV/c);n #sigma_{K}^{TOF}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
      fRegistry.add("Track/hTOFNsigmaPr", "TOF n sigma pr;p_{pv} (GeV/c);n #sigma_{p}^{TOF}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
      fRegistry.add("Track/hTPCNcr2Nf", "TPC Ncr/Nfindable", kTH1F, {{200, 0, 2}}, false);
      fRegistry.add("Track/hTPCNcls2Nf", "TPC Ncls/Nfindable", kTH1F, {{200, 0, 2}}, false);
      fRegistry.add("Track/hNclsITS", "number of ITS clusters", kTH1F, {{8, -0.5, 7.5}}, false);
      fRegistry.add("Track/hChi2ITS", "chi2/number of ITS clusters", kTH1F, {{100, 0, 10}}, false);
      fRegistry.add("Track/hITSClusterMap", "ITS cluster map", kTH1F, {{128, -0.5, 127.5}}, false);
      fRegistry.add("Track/hMeanClusterSizeITS", "mean cluster size ITS;p_{pv} (GeV/c);<cluster size> on ITS #times cos(#lambda)", kTH2F, {{1000, 0, 10}, {160, 0, 16}}, false);
      fRegistry.add("Pair/before/hMvsPt", "m_{ee} vs. p_{T,ee};m_{ee} (GeV/c^{2});p_{T,ee} (GeV/c)", kTH2F, {{100, 0, 0.1}, {100, 0, 1}}, false);
      fRegistry.add("Pair/before/hMvsPhiV", "mee vs. phiv;#varphi_{V} (rad.);m_{ee} (GeV/c^{2})", kTH2F, {{90, 0, M_PI}, {100, 0, 0.1}}, false);
      fRegistry.addClone("Pair/before/", "Pair/after/");
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
      if (fabs(d_bz) > 1e-5) {
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

    if (track.tpcChi2NCl() > maxchi2tpc) {
      return false;
    }

    if (track.itsChi2NCl() > maxchi2its) {
      return false;
    }

    if (!track.hasITS() || !track.hasTPC()) {
      return false;
    }
    if (track.itsNCls() < minitsncls) {
      return false;
    }

    auto hits = std::count_if(itsRequirement.second.begin(), itsRequirement.second.end(), [&](auto&& requiredLayer) { return track.itsClusterMap() & (1 << requiredLayer); });
    if (hits < itsRequirement.first) {
      return false;
    }

    uint32_t itsClusterSizes = track.itsClusterSizes();
    int total_cluster_size = 0, nl = 0;
    for (unsigned int layer = 0; layer < 7; layer++) {
      int cluster_size_per_layer = (itsClusterSizes >> (layer * 4)) & 0xf;
      if (cluster_size_per_layer > 0) {
        nl++;
      }
      total_cluster_size += cluster_size_per_layer;
    }
    if (static_cast<float>(total_cluster_size) / static_cast<float>(nl) * std::cos(std::atan(track.tgl())) > max_mean_its_cluster_size && track.p() < max_p_for_its_cluster_size) {
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

    gpu::gpustd::array<float, 2> dcaInfo;
    auto track_par_cov_recalc = getTrackParCov(track);
    track_par_cov_recalc.setPID(o2::track::PID::Electron);
    // std::array<float, 3> pVec_recalc = {0, 0, 0}; // px, py, pz
    o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, track_par_cov_recalc, 2.f, matCorr, &dcaInfo);
    // getPxPyPz(track_par_cov_recalc, pVec_recalc);
    float dcaXY = dcaInfo[0];
    float dcaZ = dcaInfo[1];

    if (fabs(dcaXY) > dca_xy_max || fabs(dcaZ) > dca_z_max) {
      return false;
    }

    if (track_par_cov_recalc.getPt() < minpt || fabs(track_par_cov_recalc.getEta()) > maxeta) {
      return false;
    }

    float dca_3d = 999.f;
    float det = track_par_cov_recalc.getSigmaY2() * track_par_cov_recalc.getSigmaZ2() - track_par_cov_recalc.getSigmaZY() * track_par_cov_recalc.getSigmaZY();
    if (det < 0) {
      dca_3d = 999.f;
    } else {
      float chi2 = (dcaXY * dcaXY * track_par_cov_recalc.getSigmaZ2() + dcaZ * dcaZ * track_par_cov_recalc.getSigmaY2() - 2. * dcaXY * dcaZ * track_par_cov_recalc.getSigmaZY()) / det;
      dca_3d = std::sqrt(fabs(chi2) / 2.);
    }
    if (dca_3d > dca_3d_sigma_max) {
      return false;
    }

    return true;
  }

  template <typename TTrack>
  bool isElectron(TTrack const& track)
  {
    if (track.tpcInnerParam() < max_pin_for_pion_rejection && (minTPCNsigmaPi < track.tpcNSigmaPi() && track.tpcNSigmaPi() < maxTPCNsigmaPi)) {
      return false;
    }
    if (track.hasTOF()) {
      return minTPCNsigmaEl < track.tpcNSigmaEl() && track.tpcNSigmaEl() < maxTPCNsigmaEl && fabs(track.tofNSigmaEl()) < maxTOFNsigmaEl;
    } else {
      return minTPCNsigmaEl < track.tpcNSigmaEl() && track.tpcNSigmaEl() < maxTPCNsigmaEl;
    }
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
                         track.tpcNClsFindable(), track.tpcNClsFindableMinusFound(), track.tpcNClsFindableMinusCrossedRows(),
                         track.tpcChi2NCl(), track.tpcInnerParam(),
                         track.tpcSignal(), track.tpcNSigmaEl(), track.tpcNSigmaMu(), track.tpcNSigmaPi(), track.tpcNSigmaKa(), track.tpcNSigmaPr(),
                         track.beta(), track.tofNSigmaEl(), track.tofNSigmaMu(), track.tofNSigmaPi(), track.tofNSigmaKa(), track.tofNSigmaPr(),
                         track.itsClusterSizes(), track.itsChi2NCl(), track.detectorMap(),
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
        fRegistry.fill(HIST("Track/hRelSigma1Pt"), pt_recalc, std::sqrt(track_par_cov_recalc.getSigma1Pt2()) * pt_recalc);
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
        fRegistry.fill(HIST("Track/hChi2TPC"), track.tpcChi2NCl());
        fRegistry.fill(HIST("Track/hChi2ITS"), track.itsChi2NCl());
        fRegistry.fill(HIST("Track/hITSClusterMap"), track.itsClusterMap());
        fRegistry.fill(HIST("Track/hMeanClusterSizeITS"), track.p(), static_cast<float>(total_cluster_size) / static_cast<float>(nl) * std::cos(std::atan(track.tgl())));
        fRegistry.fill(HIST("Track/hTPCdEdx"), track.tpcInnerParam(), track.tpcSignal());
        fRegistry.fill(HIST("Track/hTPCNsigmaEl"), track.tpcInnerParam(), track.tpcNSigmaEl());
        fRegistry.fill(HIST("Track/hTPCNsigmaMu"), track.tpcInnerParam(), track.tpcNSigmaMu());
        fRegistry.fill(HIST("Track/hTPCNsigmaPi"), track.tpcInnerParam(), track.tpcNSigmaPi());
        fRegistry.fill(HIST("Track/hTPCNsigmaKa"), track.tpcInnerParam(), track.tpcNSigmaKa());
        fRegistry.fill(HIST("Track/hTPCNsigmaPr"), track.tpcInnerParam(), track.tpcNSigmaPr());
        fRegistry.fill(HIST("Track/hTOFbeta"), track.p(), track.beta());
        fRegistry.fill(HIST("Track/hTOFNsigmaEl"), track.p(), track.tofNSigmaEl());
        fRegistry.fill(HIST("Track/hTOFNsigmaMu"), track.p(), track.tofNSigmaMu());
        fRegistry.fill(HIST("Track/hTOFNsigmaPi"), track.p(), track.tofNSigmaPi());
        fRegistry.fill(HIST("Track/hTOFNsigmaKa"), track.p(), track.tofNSigmaKa());
        fRegistry.fill(HIST("Track/hTOFNsigmaPr"), track.p(), track.tofNSigmaPr());
      }
    }
  }

  template <typename TCollision, typename TTrack>
  o2::track::TrackParCov propagateTrack(TCollision const& collision, TTrack const& track)
  {
    gpu::gpustd::array<float, 2> dcaInfo;
    auto track_par_cov_recalc = getTrackParCov(track);
    track_par_cov_recalc.setPID(o2::track::PID::Electron);
    // std::array<float, 3> pVec_recalc = {0, 0, 0}; // px, py, pz
    o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, track_par_cov_recalc, 2.f, matCorr, &dcaInfo);
    // getPxPyPz(track_par_cov_recalc, pVec_recalc);
    return track_par_cov_recalc;
  }

  std::vector<std::pair<int, int>> stored_trackIds;
  std::vector<std::pair<int, int>> stored_pairIds;
  Filter trackFilter = o2::aod::track::pt > minpt&& nabs(o2::aod::track::eta) < maxeta&& o2::aod::track::tpcChi2NCl < maxchi2tpc&& o2::aod::track::itsChi2NCl < maxchi2its&& ncheckbit(aod::track::v001::detectorMap, (uint8_t)o2::aod::track::ITS) == true && ncheckbit(aod::track::v001::detectorMap, (uint8_t)o2::aod::track::TPC) == true;
  Filter pidFilter = minTPCNsigmaEl < o2::aod::pidtpc::tpcNSigmaEl && o2::aod::pidtpc::tpcNSigmaEl < maxTPCNsigmaEl;
  using MyFilteredTracks = soa::Filtered<MyTracks>;

  Partition<MyFilteredTracks> posTracks = o2::aod::track::signed1Pt > 0.f;
  Partition<MyFilteredTracks> negTracks = o2::aod::track::signed1Pt < 0.f;

  // ---------- for data ----------

  void processRec_SA(MyCollisions const& collisions, aod::BCsWithTimestamps const&, MyFilteredTracks const&)
  {
    stored_trackIds.reserve(posTracks.size() + negTracks.size());

    for (auto& collision : collisions) {
      auto bc = collision.template foundBC_as<aod::BCsWithTimestamps>();
      initCCDB(bc);

      if (!collision.isSelected()) {
        filter(0, 0, 0);
        continue;
      }

      int nee_uls = 0;
      auto posTracks_per_coll = posTracks->sliceByCached(o2::aod::track::collisionId, collision.globalIndex(), cache);
      auto negTracks_per_coll = negTracks->sliceByCached(o2::aod::track::collisionId, collision.globalIndex(), cache);

      for (auto& [pos, ele] : combinations(CombinationsFullIndexPolicy(posTracks_per_coll, negTracks_per_coll))) {
        if (!checkTrack<false>(collision, pos) || !checkTrack<false>(collision, ele)) {
          continue;
        }
        if (!isElectron(pos) || !isElectron(ele)) {
          continue;
        }

        ROOT::Math::PtEtaPhiMVector v1(pos.pt(), pos.eta(), pos.phi(), o2::constants::physics::MassElectron);
        ROOT::Math::PtEtaPhiMVector v2(ele.pt(), ele.eta(), ele.phi(), o2::constants::physics::MassElectron);
        ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
        float phiv = o2::aod::pwgem::dilepton::utils::pairutil::getPhivPair(pos.px(), pos.py(), pos.pz(), ele.px(), ele.py(), ele.pz(), pos.sign(), ele.sign(), d_bz);

        if (fillQAHistogram) {
          fRegistry.fill(HIST("Pair/before/hMvsPt"), v12.M(), v12.Pt());
          fRegistry.fill(HIST("Pair/before/hMvsPhiV"), phiv, v12.M());
        }
        if (apply_phiv ? (v12.M() < maxMee && slope * phiv + intercept < v12.M()) : (v12.M() < maxMee)) {
          fillTrackTable(collision, pos);
          fillTrackTable(collision, ele);
          if (fillQAHistogram) {
            fRegistry.fill(HIST("Pair/after/hMvsPt"), v12.M(), v12.Pt());
            fRegistry.fill(HIST("Pair/after/hMvsPhiV"), phiv, v12.M());
          }
          nee_uls++;
        }

      } // end of pairing loop

      if (nee_uls < 1) {
        filter(nee_uls, 0, 0);
        continue;
      }
      filter(nee_uls, 0, 0);

    } // end of collision loop

    stored_trackIds.clear();
    stored_trackIds.shrink_to_fit();
    stored_pairIds.clear();
    stored_pairIds.shrink_to_fit();
  }
  PROCESS_SWITCH(filterDielectronEvent, processRec_SA, "process reconstructed info only", true); // standalone

  Preslice<aod::TrackAssoc> trackIndicesPerCollision = aod::track_association::collisionId;
  void processRec_TTCA(MyCollisions const& collisions, aod::BCsWithTimestamps const&, MyTracks const& tracks, aod::TrackAssoc const& trackIndices)
  {
    stored_trackIds.reserve(tracks.size() * 2);

    for (auto& collision : collisions) {
      auto bc = collision.template foundBC_as<aod::BCsWithTimestamps>();
      initCCDB(bc);

      if (!collision.isSelected()) {
        filter(0, 0, 0);
        continue;
      }

      int nee_uls = 0;
      auto trackIdsThisCollision = trackIndices.sliceBy(trackIndicesPerCollision, collision.globalIndex());
      std::vector<MyTrack> posTracks_per_coll;
      std::vector<MyTrack> negTracks_per_coll;
      posTracks_per_coll.reserve(trackIdsThisCollision.size());
      negTracks_per_coll.reserve(trackIdsThisCollision.size());

      for (auto& trackId : trackIdsThisCollision) {
        auto track = trackId.template track_as<MyTracks>();
        if (!checkTrack<false>(collision, track) || !isElectron(track)) {
          continue;
        }

        if (track.sign() > 0) {
          posTracks_per_coll.emplace_back(track);
        } else {
          negTracks_per_coll.emplace_back(track);
        }
      } // end of track loop

      for (auto& pos : posTracks_per_coll) {
        for (auto& ele : negTracks_per_coll) {

          auto pos_prop = propagateTrack(collision, pos);
          auto ele_prop = propagateTrack(collision, ele);

          std::array<float, 3> pVec_pos = {0, 0, 0}; // px, py, pz
          getPxPyPz(pos_prop, pVec_pos);
          std::array<float, 3> pVec_ele = {0, 0, 0}; // px, py, pz
          getPxPyPz(ele_prop, pVec_ele);

          ROOT::Math::PtEtaPhiMVector v1(pos_prop.getPt(), pos_prop.getEta(), pos_prop.getPhi(), o2::constants::physics::MassElectron);
          ROOT::Math::PtEtaPhiMVector v2(ele_prop.getPt(), ele_prop.getEta(), ele_prop.getPhi(), o2::constants::physics::MassElectron);
          ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
          float phiv = o2::aod::pwgem::dilepton::utils::pairutil::getPhivPair(pVec_pos[0], pVec_pos[1], pVec_pos[2], pVec_ele[0], pVec_ele[1], pVec_ele[2], pos.sign(), ele.sign(), d_bz);

          if (fillQAHistogram) {
            fRegistry.fill(HIST("Pair/before/hMvsPt"), v12.M(), v12.Pt());
            fRegistry.fill(HIST("Pair/before/hMvsPhiV"), phiv, v12.M());
          }
          if (apply_phiv ? (v12.M() < maxMee && slope * phiv + intercept < v12.M()) : (v12.M() < maxMee)) {
            fillTrackTable(collision, pos);
            fillTrackTable(collision, ele);
            if (fillQAHistogram) {
              fRegistry.fill(HIST("Pair/after/hMvsPt"), v12.M(), v12.Pt());
              fRegistry.fill(HIST("Pair/after/hMvsPhiV"), phiv, v12.M());
            }
            nee_uls++;
          }

        } // end of negative track loop
      } // end of postive track loop

      if (nee_uls < 1) {
        filter(nee_uls, 0, 0);
        continue;
      }

      filter(nee_uls, 0, 0);

      posTracks_per_coll.clear();
      negTracks_per_coll.clear();
      posTracks_per_coll.shrink_to_fit();
      negTracks_per_coll.shrink_to_fit();

    } // end of collision loop

    stored_trackIds.clear();
    stored_trackIds.shrink_to_fit();
    stored_pairIds.clear();
    stored_pairIds.shrink_to_fit();
  }
  PROCESS_SWITCH(filterDielectronEvent, processRec_TTCA, "process reconstructed info only", false); // with TTCA

  // ---------- for data ----------

  void processRec_SA_SWT(MyCollisionsWithSWT const& collisions, aod::BCsWithTimestamps const&, MyFilteredTracks const&)
  {
    stored_trackIds.reserve(posTracks.size() + negTracks.size());

    for (auto& collision : collisions) {
      auto bc = collision.template foundBC_as<aod::BCsWithTimestamps>();
      initCCDB(bc);

      if (!collision.isSelected()) {
        filter(0, 0, 0);
        continue;
      }
      if (collision.swtaliastmp_raw() == 0) {
        filter(0, 0, 0);
        continue;
      }

      int nee_uls = 0;
      auto posTracks_per_coll = posTracks->sliceByCached(o2::aod::track::collisionId, collision.globalIndex(), cache);
      auto negTracks_per_coll = negTracks->sliceByCached(o2::aod::track::collisionId, collision.globalIndex(), cache);

      for (auto& [pos, ele] : combinations(CombinationsFullIndexPolicy(posTracks_per_coll, negTracks_per_coll))) {
        if (!checkTrack<false>(collision, pos) || !checkTrack<false>(collision, ele)) {
          continue;
        }
        if (!isElectron(pos) || !isElectron(ele)) {
          continue;
        }

        ROOT::Math::PtEtaPhiMVector v1(pos.pt(), pos.eta(), pos.phi(), o2::constants::physics::MassElectron);
        ROOT::Math::PtEtaPhiMVector v2(ele.pt(), ele.eta(), ele.phi(), o2::constants::physics::MassElectron);
        ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
        float phiv = o2::aod::pwgem::dilepton::utils::pairutil::getPhivPair(pos.px(), pos.py(), pos.pz(), ele.px(), ele.py(), ele.pz(), pos.sign(), ele.sign(), d_bz);

        if (fillQAHistogram) {
          fRegistry.fill(HIST("Pair/before/hMvsPt"), v12.M(), v12.Pt());
          fRegistry.fill(HIST("Pair/before/hMvsPhiV"), phiv, v12.M());
        }
        if (apply_phiv ? (v12.M() < maxMee && slope * phiv + intercept < v12.M()) : (v12.M() < maxMee)) {
          fillTrackTable(collision, pos);
          fillTrackTable(collision, ele);
          if (fillQAHistogram) {
            fRegistry.fill(HIST("Pair/after/hMvsPt"), v12.M(), v12.Pt());
            fRegistry.fill(HIST("Pair/after/hMvsPhiV"), phiv, v12.M());
          }
          nee_uls++;
        }

      } // end of pairing loop

      if (nee_uls < 1) {
        filter(nee_uls, 0, 0);
        continue;
      }
      filter(nee_uls, 0, 0);

    } // end of collision loop

    stored_trackIds.clear();
    stored_trackIds.shrink_to_fit();
    stored_pairIds.clear();
    stored_pairIds.shrink_to_fit();
  }
  PROCESS_SWITCH(filterDielectronEvent, processRec_SA_SWT, "process reconstructed info only", false); // standalone

  void processRec_TTCA_SWT(MyCollisionsWithSWT const& collisions, aod::BCsWithTimestamps const&, MyTracks const& tracks, aod::TrackAssoc const& trackIndices)
  {
    stored_trackIds.reserve(tracks.size() * 2);

    for (auto& collision : collisions) {
      auto bc = collision.template foundBC_as<aod::BCsWithTimestamps>();
      initCCDB(bc);

      if (!collision.isSelected()) {
        filter(0, 0, 0);
        continue;
      }
      if (collision.swtaliastmp_raw() == 0) {
        filter(0, 0, 0);
        continue;
      }

      int nee_uls = 0;
      auto trackIdsThisCollision = trackIndices.sliceBy(trackIndicesPerCollision, collision.globalIndex());
      std::vector<MyTrack> posTracks_per_coll;
      std::vector<MyTrack> negTracks_per_coll;
      posTracks_per_coll.reserve(trackIdsThisCollision.size());
      negTracks_per_coll.reserve(trackIdsThisCollision.size());

      for (auto& trackId : trackIdsThisCollision) {
        auto track = trackId.template track_as<MyTracks>();
        if (!checkTrack<false>(collision, track) || !isElectron(track)) {
          continue;
        }

        if (track.sign() > 0) {
          posTracks_per_coll.emplace_back(track);
        } else {
          negTracks_per_coll.emplace_back(track);
        }
      } // end of track loop

      for (auto& pos : posTracks_per_coll) {
        for (auto& ele : negTracks_per_coll) {

          auto pos_prop = propagateTrack(collision, pos);
          auto ele_prop = propagateTrack(collision, ele);

          std::array<float, 3> pVec_pos = {0, 0, 0}; // px, py, pz
          getPxPyPz(pos_prop, pVec_pos);
          std::array<float, 3> pVec_ele = {0, 0, 0}; // px, py, pz
          getPxPyPz(ele_prop, pVec_ele);

          ROOT::Math::PtEtaPhiMVector v1(pos_prop.getPt(), pos_prop.getEta(), pos_prop.getPhi(), o2::constants::physics::MassElectron);
          ROOT::Math::PtEtaPhiMVector v2(ele_prop.getPt(), ele_prop.getEta(), ele_prop.getPhi(), o2::constants::physics::MassElectron);
          ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
          float phiv = o2::aod::pwgem::dilepton::utils::pairutil::getPhivPair(pVec_pos[0], pVec_pos[1], pVec_pos[2], pVec_ele[0], pVec_ele[1], pVec_ele[2], pos.sign(), ele.sign(), d_bz);

          if (fillQAHistogram) {
            fRegistry.fill(HIST("Pair/before/hMvsPt"), v12.M(), v12.Pt());
            fRegistry.fill(HIST("Pair/before/hMvsPhiV"), phiv, v12.M());
          }
          if (apply_phiv ? (v12.M() < maxMee && slope * phiv + intercept < v12.M()) : (v12.M() < maxMee)) {
            fillTrackTable(collision, pos);
            fillTrackTable(collision, ele);
            if (fillQAHistogram) {
              fRegistry.fill(HIST("Pair/after/hMvsPt"), v12.M(), v12.Pt());
              fRegistry.fill(HIST("Pair/after/hMvsPhiV"), phiv, v12.M());
            }
            nee_uls++;
          }

        } // end of negative track loop
      } // end of postive track loop

      if (nee_uls < 1) {
        filter(nee_uls, 0, 0);
        continue;
      }

      filter(nee_uls, 0, 0);

      posTracks_per_coll.clear();
      negTracks_per_coll.clear();
      posTracks_per_coll.shrink_to_fit();
      negTracks_per_coll.shrink_to_fit();

    } // end of collision loop

    stored_trackIds.clear();
    stored_trackIds.shrink_to_fit();
    stored_pairIds.clear();
    stored_pairIds.shrink_to_fit();
  }
  PROCESS_SWITCH(filterDielectronEvent, processRec_TTCA_SWT, "process reconstructed info only", false); // with TTCA

  // ---------- for MC ----------

  using MyFilteredTracksMC = soa::Filtered<MyTracksMC>;
  Partition<MyFilteredTracksMC> posTracksMC = o2::aod::track::signed1Pt > 0.f;
  Partition<MyFilteredTracksMC> negTracksMC = o2::aod::track::signed1Pt < 0.f;
  void processMC_SA(soa::Join<MyCollisions, aod::McCollisionLabels> const& collisions, aod::McCollisions const&, aod::BCsWithTimestamps const&, MyFilteredTracksMC const& tracks)
  {
    stored_trackIds.reserve(tracks.size());

    for (auto& collision : collisions) {
      if (!collision.has_mcCollision()) {
        continue;
      }
      auto bc = collision.template foundBC_as<aod::BCsWithTimestamps>();
      initCCDB(bc);

      if (!collision.isSelected()) {
        filter(0, 0, 0);
        continue;
      }

      int nee_uls = 0;
      auto posTracks_per_coll = posTracksMC->sliceByCached(o2::aod::track::collisionId, collision.globalIndex(), cache);
      auto negTracks_per_coll = negTracksMC->sliceByCached(o2::aod::track::collisionId, collision.globalIndex(), cache);

      for (auto& [pos, ele] : combinations(CombinationsFullIndexPolicy(posTracks_per_coll, negTracks_per_coll))) {
        if (!checkTrack<true>(collision, pos) || !checkTrack<true>(collision, ele)) {
          continue;
        }
        if (!isElectron(pos) || !isElectron(ele)) {
          continue;
        }

        ROOT::Math::PtEtaPhiMVector v1(pos.pt(), pos.eta(), pos.phi(), o2::constants::physics::MassElectron);
        ROOT::Math::PtEtaPhiMVector v2(ele.pt(), ele.eta(), ele.phi(), o2::constants::physics::MassElectron);
        ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
        float phiv = o2::aod::pwgem::dilepton::utils::pairutil::getPhivPair(pos.px(), pos.py(), pos.pz(), ele.px(), ele.py(), ele.pz(), pos.sign(), ele.sign(), d_bz);
        if (fillQAHistogram) {
          fRegistry.fill(HIST("Pair/before/hMvsPt"), v12.M(), v12.Pt());
          fRegistry.fill(HIST("Pair/before/hMvsPhiV"), phiv, v12.M());
        }
        if (apply_phiv ? (v12.M() < maxMee && slope * phiv + intercept < v12.M()) : (v12.M() < maxMee)) {
          fillTrackTable(collision, pos);
          fillTrackTable(collision, ele);
          if (fillQAHistogram) {
            fRegistry.fill(HIST("Pair/after/hMvsPt"), v12.M(), v12.Pt());
            fRegistry.fill(HIST("Pair/after/hMvsPhiV"), phiv, v12.M());
          }
          nee_uls++;
        }

      } // end of pairing loop

      if (nee_uls < 1) {
        filter(nee_uls, 0, 0);
        continue;
      }
      filter(nee_uls, 0, 0);
    } // end of collision loop

    stored_trackIds.clear();
    stored_trackIds.shrink_to_fit();
    stored_pairIds.clear();
    stored_pairIds.shrink_to_fit();
  }
  PROCESS_SWITCH(filterDielectronEvent, processMC_SA, "process reconstructed and MC info ", false);

  void processMC_TTCA(soa::Join<MyCollisions, aod::McCollisionLabels> const& collisions, aod::McCollisions const&, aod::BCsWithTimestamps const&, MyTracksMC const& tracks, aod::TrackAssoc const& trackIndices)
  {
    stored_trackIds.reserve(tracks.size() * 2);

    for (auto& collision : collisions) {
      if (!collision.has_mcCollision()) {
        continue;
      }
      auto bc = collision.template foundBC_as<aod::BCsWithTimestamps>();
      initCCDB(bc);

      if (!collision.isSelected()) {
        filter(0, 0, 0);
        continue;
      }

      int nee_uls = 0;
      auto trackIdsThisCollision = trackIndices.sliceBy(trackIndicesPerCollision, collision.globalIndex());
      std::vector<MyTrackMC> posTracks_per_coll;
      std::vector<MyTrackMC> negTracks_per_coll;
      posTracks_per_coll.reserve(trackIdsThisCollision.size());
      negTracks_per_coll.reserve(trackIdsThisCollision.size());

      for (auto& trackId : trackIdsThisCollision) {
        auto track = trackId.template track_as<MyTracksMC>();
        if (!checkTrack<true>(collision, track) || !isElectron(track)) {
          continue;
        }

        if (track.sign() > 0) {
          posTracks_per_coll.emplace_back(track);
        } else {
          negTracks_per_coll.emplace_back(track);
        }
      } // end of track loop

      for (auto& pos : posTracks_per_coll) {
        for (auto& ele : negTracks_per_coll) {
          auto pos_prop = propagateTrack(collision, pos);
          auto ele_prop = propagateTrack(collision, ele);
          std::array<float, 3> pVec_pos = {0, 0, 0}; // px, py, pz
          getPxPyPz(pos_prop, pVec_pos);
          std::array<float, 3> pVec_ele = {0, 0, 0}; // px, py, pz
          getPxPyPz(ele_prop, pVec_ele);

          ROOT::Math::PtEtaPhiMVector v1(pos_prop.getPt(), pos_prop.getEta(), pos_prop.getPhi(), o2::constants::physics::MassElectron);
          ROOT::Math::PtEtaPhiMVector v2(ele_prop.getPt(), ele_prop.getEta(), ele_prop.getPhi(), o2::constants::physics::MassElectron);
          ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
          float phiv = o2::aod::pwgem::dilepton::utils::pairutil::getPhivPair(pVec_pos[0], pVec_pos[1], pVec_pos[2], pVec_ele[0], pVec_ele[1], pVec_ele[2], pos.sign(), ele.sign(), d_bz);
          if (fillQAHistogram) {
            fRegistry.fill(HIST("Pair/before/hMvsPt"), v12.M(), v12.Pt());
            fRegistry.fill(HIST("Pair/before/hMvsPhiV"), phiv, v12.M());
          }
          if (apply_phiv ? (v12.M() < maxMee && slope * phiv + intercept < v12.M()) : (v12.M() < maxMee)) {
            fillTrackTable(collision, pos);
            fillTrackTable(collision, ele);
            if (fillQAHistogram) {
              fRegistry.fill(HIST("Pair/after/hMvsPt"), v12.M(), v12.Pt());
              fRegistry.fill(HIST("Pair/after/hMvsPhiV"), phiv, v12.M());
            }
            nee_uls++;
          }

        } // end of negative track loop
      } // end of postive track loop

      if (nee_uls < 1) {
        filter(nee_uls, 0, 0);
        continue;
      }
      filter(nee_uls, 0, 0);

      posTracks_per_coll.clear();
      negTracks_per_coll.clear();
      posTracks_per_coll.shrink_to_fit();
      negTracks_per_coll.shrink_to_fit();

    } // end of collision loop

    stored_trackIds.clear();
    stored_trackIds.shrink_to_fit();
    stored_pairIds.clear();
    stored_pairIds.shrink_to_fit();
  }
  PROCESS_SWITCH(filterDielectronEvent, processMC_TTCA, "process reconstructed info only", false); // with TTCA
};
struct prefilterPrimaryElectron {
  Produces<aod::EMPrimaryElectronsPrefilterBit> ele_pfb;
  void process(aod::EMPrimaryElectrons const& primaryelectrons)
  {
    for (int i = 0; i < primaryelectrons.size(); i++) {
      ele_pfb(0);
    }
  }
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
struct createEMEvent2VP {
  using MyBCs = soa::Join<aod::BCsWithTimestamps, aod::BcSels>;
  using MyQvectors = soa::Join<aod::QvectorFT0CVecs, aod::QvectorFT0AVecs, aod::QvectorFT0MVecs, aod::QvectorBPosVecs, aod::QvectorBNegVecs, aod::QvectorBTotVecs>;

  using MyCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::EMEventsNee, aod::EMEventsNgPCM, aod::EMEvSels>;
  using MyCollisions_Cent = soa::Join<MyCollisions, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::CentNTPVs>; // centrality table has dependency on multiplicity table.
  using MyCollisions_Cent_Qvec = soa::Join<MyCollisions_Cent, MyQvectors>;

  using MyCollisionsWithSWT = soa::Join<MyCollisions, aod::EMSWTriggerInfosTMP>;
  using MyCollisionsWithSWT_Cent = soa::Join<MyCollisionsWithSWT, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::CentNTPVs>; // centrality table has dependency on multiplicity table.
  using MyCollisionsWithSWT_Cent_Qvec = soa::Join<MyCollisionsWithSWT_Cent, MyQvectors>;

  using MyCollisionsMC = soa::Join<MyCollisions, aod::McCollisionLabels>;
  using MyCollisionsMC_Cent = soa::Join<MyCollisionsMC, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::CentNTPVs>; // centrality table has dependency on multiplicity table.
  using MyCollisionsMC_Cent_Qvec = soa::Join<MyCollisionsMC_Cent, MyQvectors>;

  Produces<o2::aod::EMEvents> event;
  // Produces<o2::aod::EMEventsCov> eventcov;
  Produces<o2::aod::EMEventsMult> event_mult;
  Produces<o2::aod::EMEventsCent> event_cent;
  Produces<o2::aod::EMEventsQvec> event_qvec;
  Produces<o2::aod::EMSWTriggerInfos> emswtbit;

  enum class EMEventType : int {
    kEvent = 0,
    kEvent_Cent = 1,
    kEvent_Cent_Qvec = 2,
  };

  HistogramRegistry registry{"registry"};
  void init(o2::framework::InitContext&)
  {
    auto hEventCounter = registry.add<TH1>("hEventCounter", "hEventCounter", kTH1I, {{7, 0.5f, 7.5f}});
    hEventCounter->GetXaxis()->SetBinLabel(1, "all");
    hEventCounter->GetXaxis()->SetBinLabel(2, "sel8");

    registry.add("hNInspectedTVX", "N inspected TVX;run number;N_{TVX}", kTProfile, {{80000, 520000.5, 600000.5}}, true);
  }

  ~createEMEvent2VP()
  {
    swt_names.clear();
    swt_names.shrink_to_fit();
  }

  std::vector<int> mTOIidx;
  std::vector<std::string> swt_names;
  uint64_t mNinspectedTVX{0};

  int mRunNumber;

  template <bool isMC, bool isTriggerAnalysis, EMEventType eventype, typename TCollisions, typename TBCs>
  void skimEvent(TCollisions const& collisions, TBCs const&)
  {
    for (auto& collision : collisions) {
      if constexpr (isMC) {
        if (!collision.has_mcCollision()) {
          continue;
        }
      }

      if constexpr (isTriggerAnalysis) {
        if (collision.swtaliastmp_raw() == 0) {
          continue;
        }
      }

      auto bc = collision.template foundBC_as<TBCs>();

      if (!collision.isSelected()) {
        continue;
      }

      if (!(collision.neeuls() >= 1 || collision.neeuls() + collision.ngpcm() >= 2)) {
        continue;
      }

      if constexpr (isTriggerAnalysis) {
        emswtbit(collision.swtaliastmp_raw(), collision.nInspectedTVX());
      }

      // LOGF(info, "collision.neeuls() = %d, collision.ngpcm() = %d", collision.neeuls(), collision.ngpcm());
      // LOGF(info, "collision.multNTracksPV() = %d, collision.multFT0A() = %f, collision.multFT0C() = %f", collision.multNTracksPV(), collision.multFT0A(), collision.multFT0C());

      registry.fill(HIST("hEventCounter"), 1);

      event(collision.globalIndex(), bc.runNumber(), bc.globalBC(), collision.alias_raw(), collision.selection_raw(), bc.timestamp(),
            collision.posX(), collision.posY(), collision.posZ(),
            collision.numContrib(), collision.trackOccupancyInTimeRange());

      // eventcov(collision.covXX(), collision.covXY(), collision.covXZ(), collision.covYY(), collision.covYZ(), collision.covZZ(), collision.chi2());

      event_mult(collision.multFT0A(), collision.multFT0C(), collision.multTPC(), collision.multNTracksPV(), collision.multNTracksPVeta1(), collision.multNTracksPVetaHalf());

      if constexpr (eventype == EMEventType::kEvent) {
        event_cent(105.f, 105.f, 105.f, 105.f);
        event_qvec(
          999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f,
          999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f);
      } else if constexpr (eventype == EMEventType::kEvent_Cent) {
        event_cent(collision.centFT0M(), collision.centFT0A(), collision.centFT0C(), collision.centNTPV());
        event_qvec(
          999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f,
          999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f);
      } else if constexpr (eventype == EMEventType::kEvent_Cent_Qvec) {
        event_cent(collision.centFT0M(), collision.centFT0A(), collision.centFT0C(), collision.centNTPV());
        float q2xft0m = 999.f, q2yft0m = 999.f, q2xft0a = 999.f, q2yft0a = 999.f, q2xft0c = 999.f, q2yft0c = 999.f, q2xbpos = 999.f, q2ybpos = 999.f, q2xbneg = 999.f, q2ybneg = 999.f, q2xbtot = 999.f, q2ybtot = 999.f;
        float q3xft0m = 999.f, q3yft0m = 999.f, q3xft0a = 999.f, q3yft0a = 999.f, q3xft0c = 999.f, q3yft0c = 999.f, q3xbpos = 999.f, q3ybpos = 999.f, q3xbneg = 999.f, q3ybneg = 999.f, q3xbtot = 999.f, q3ybtot = 999.f;

        if (collision.qvecFT0CReVec().size() >= 2) { // harmonics 2,3
          q2xft0m = collision.qvecFT0MReVec()[0], q2xft0a = collision.qvecFT0AReVec()[0], q2xft0c = collision.qvecFT0CReVec()[0], q2xbpos = collision.qvecBPosReVec()[0], q2xbneg = collision.qvecBNegReVec()[0], q2xbtot = collision.qvecBTotReVec()[0];
          q2yft0m = collision.qvecFT0MImVec()[0], q2yft0a = collision.qvecFT0AImVec()[0], q2yft0c = collision.qvecFT0CImVec()[0], q2ybpos = collision.qvecBPosImVec()[0], q2ybneg = collision.qvecBNegImVec()[0], q2ybtot = collision.qvecBTotImVec()[0];
          q3xft0m = collision.qvecFT0MReVec()[1], q3xft0a = collision.qvecFT0AReVec()[1], q3xft0c = collision.qvecFT0CReVec()[1], q3xbpos = collision.qvecBPosReVec()[1], q3xbneg = collision.qvecBNegReVec()[1], q3xbtot = collision.qvecBTotReVec()[1];
          q3yft0m = collision.qvecFT0MImVec()[1], q3yft0a = collision.qvecFT0AImVec()[1], q3yft0c = collision.qvecFT0CImVec()[1], q3ybpos = collision.qvecBPosImVec()[1], q3ybneg = collision.qvecBNegImVec()[1], q3ybtot = collision.qvecBTotImVec()[1];
        } else if (collision.qvecFT0CReVec().size() >= 1) { // harmonics 2
          q2xft0m = collision.qvecFT0MReVec()[0], q2xft0a = collision.qvecFT0AReVec()[0], q2xft0c = collision.qvecFT0CReVec()[0], q2xbpos = collision.qvecBPosReVec()[0], q2xbneg = collision.qvecBNegReVec()[0], q2xbtot = collision.qvecBTotReVec()[0];
          q2yft0m = collision.qvecFT0MImVec()[0], q2yft0a = collision.qvecFT0AImVec()[0], q2yft0c = collision.qvecFT0CImVec()[0], q2ybpos = collision.qvecBPosImVec()[0], q2ybneg = collision.qvecBNegImVec()[0], q2ybtot = collision.qvecBTotImVec()[0];
        }
        event_qvec(
          q2xft0m, q2yft0m, q2xft0a, q2yft0a, q2xft0c, q2yft0c, q2xbpos, q2ybpos, q2xbneg, q2ybneg, q2xbtot, q2ybtot,
          q3xft0m, q3yft0m, q3xft0a, q3yft0a, q3xft0c, q3yft0c, q3xbpos, q3ybpos, q3xbneg, q3ybneg, q3xbtot, q3ybtot);
      } else {
        event_cent(105.f, 105.f, 105.f, 105.f);
        event_qvec(
          999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f,
          999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f);
      }
    } // end of collision loop
  } // end of skimEvent

  void processEvent(MyCollisions const& collisions, MyBCs const& bcs)
  {
    skimEvent<false, false, EMEventType::kEvent>(collisions, bcs);
  }
  PROCESS_SWITCH(createEMEvent2VP, processEvent, "process event info", false);

  void processEvent_Cent(MyCollisions_Cent const& collisions, MyBCs const& bcs)
  {
    skimEvent<false, false, EMEventType::kEvent_Cent>(collisions, bcs);
  }
  PROCESS_SWITCH(createEMEvent2VP, processEvent_Cent, "process event info", false);

  void processEvent_Cent_Qvec(MyCollisions_Cent_Qvec const& collisions, MyBCs const& bcs)
  {
    skimEvent<false, false, EMEventType::kEvent_Cent_Qvec>(collisions, bcs);
  }
  PROCESS_SWITCH(createEMEvent2VP, processEvent_Cent_Qvec, "process event info", false);

  void processEvent_SWT(MyCollisionsWithSWT const& collisions, MyBCs const& bcs)
  {
    skimEvent<false, true, EMEventType::kEvent>(collisions, bcs);
  }
  PROCESS_SWITCH(createEMEvent2VP, processEvent_SWT, "process event info", false);

  void processEvent_SWT_Cent(MyCollisionsWithSWT_Cent const& collisions, MyBCs const& bcs)
  {
    skimEvent<false, true, EMEventType::kEvent_Cent>(collisions, bcs);
  }
  PROCESS_SWITCH(createEMEvent2VP, processEvent_SWT_Cent, "process event info", false);

  void processEvent_SWT_Cent_Qvec(MyCollisionsWithSWT_Cent_Qvec const& collisions, MyBCs const& bcs)
  {
    skimEvent<false, true, EMEventType::kEvent_Cent_Qvec>(collisions, bcs);
  }
  PROCESS_SWITCH(createEMEvent2VP, processEvent_SWT_Cent_Qvec, "process event info", false);

  void processEventMC(MyCollisionsMC const& collisions, MyBCs const& bcs)
  {
    skimEvent<true, false, EMEventType::kEvent>(collisions, bcs);
  }
  PROCESS_SWITCH(createEMEvent2VP, processEventMC, "process event info", false);

  void processEventMC_Cent(MyCollisionsMC_Cent const& collisions, MyBCs const& bcs)
  {
    skimEvent<true, false, EMEventType::kEvent_Cent>(collisions, bcs);
  }
  PROCESS_SWITCH(createEMEvent2VP, processEventMC_Cent, "process event info", false);

  void processEventMC_Cent_Qvec(MyCollisionsMC_Cent_Qvec const& collisions, MyBCs const& bcs)
  {
    skimEvent<true, false, EMEventType::kEvent_Cent_Qvec>(collisions, bcs);
  }
  PROCESS_SWITCH(createEMEvent2VP, processEventMC_Cent_Qvec, "process event info", false);

  void processDummy(aod::Collisions const&) {}
  PROCESS_SWITCH(createEMEvent2VP, processDummy, "processDummy", true);
};
struct AssociateDileptonToEMEvent2VP {
  Produces<o2::aod::V0KFEMEventIds> v0kfeventid;
  Produces<o2::aod::EMPrimaryElectronEMEventIds> prmeleventid;

  Preslice<aod::V0PhotonsKF> perCollision_pcm = aod::v0photonkf::collisionId;
  PresliceUnsorted<aod::EMPrimaryElectrons> perCollision_el = aod::emprimaryelectron::collisionId;

  void init(o2::framework::InitContext&) {}

  template <typename TCollisions, typename TLeptons, typename TEventIds, typename TPreslice>
  void fillEventId(TCollisions const& collisions, TLeptons const& leptons, TEventIds& eventIds, TPreslice const& perCollision)
  {
    for (auto& collision : collisions) {
      auto leptons_coll = leptons.sliceBy(perCollision, collision.collisionId());
      int nl = leptons_coll.size();
      // LOGF(info, "collision.collisionId() = %d , nl = %d", collision.collisionId(), nl);
      for (int il = 0; il < nl; il++) {
        eventIds(collision.globalIndex());
      } // end of photon loop
    } // end of collision loop
  }

  // This struct is for both data and MC.
  // Note that reconstructed collisions without mc collisions are already rejected in CreateEMEventDilepton in MC.

  void processPCM(aod::EMEvents const& collisions, aod::V0PhotonsKF const& photons)
  {
    fillEventId(collisions, photons, v0kfeventid, perCollision_pcm);
  }

  void processElectron(aod::EMEvents const& collisions, aod::EMPrimaryElectrons const& tracks)
  {
    fillEventId(collisions, tracks, prmeleventid, perCollision_el);
  }

  void processDummy(aod::EMEvents const&) {}

  PROCESS_SWITCH(AssociateDileptonToEMEvent2VP, processPCM, "process pcm-event indexing", false);
  PROCESS_SWITCH(AssociateDileptonToEMEvent2VP, processElectron, "process dalitzee-event indexing", false);
  PROCESS_SWITCH(AssociateDileptonToEMEvent2VP, processDummy, "process dummy", true);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<filterDielectronEvent>(cfgc, TaskName{"filter-dielectron-event"}),
    adaptAnalysisTask<prefilterPrimaryElectron>(cfgc, TaskName{"prefilter-primary-electron"}),
    adaptAnalysisTask<associateAmbiguousElectron>(cfgc, TaskName{"associate-ambiguous-electron"}),
    adaptAnalysisTask<createEMEvent2VP>(cfgc, TaskName{"create-emevent-2vp"}),
    adaptAnalysisTask<AssociateDileptonToEMEvent2VP>(cfgc, TaskName{"associate-dilepton-to-emevent2VP"}),
  };
}
