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
//
// Analysis task to produce resolution mapfor electrons/muons in dilepton analysis
//    Please write to: daiki.sekihata@cern.ch

#include <map>
#include <string>
#include <utility>
#include <set>
#include <vector>
#include <array>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/ASoA.h"
#include "Framework/DataTypes.h"
#include "Framework/HistogramRegistry.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"

#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "TGeoGlobalMagField.h"
#include "Field/MagneticField.h"

#include "DetectorsBase/Propagator.h"
#include "GlobalTracking/MatchGlobalFwd.h"
#include "MCHTracking/TrackExtrap.h"
#include "MCHTracking/TrackParam.h"
#include "ReconstructionDataFormats/TrackFwd.h"

#include "PWGEM/Dilepton/DataModel/dileptonTables.h"
#include "PWGEM/Dilepton/Utils/MCUtilities.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod;
using namespace o2::soa;
using namespace o2::aod::pwgem::dilepton::utils::mcutil;

using MyCollisions = Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>;
using MyCollision = MyCollisions::iterator;

using MyCollisionsCent = soa::Join<MyCollisions, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs>;
using MyCollisionCent = MyCollisionsCent::iterator;

using MyMCTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TracksCov, aod::McTrackLabels>;
using MyMCTrack = MyMCTracks::iterator;

using MyMCFwdTracks = soa::Join<aod::FwdTracks, aod::FwdTracksCov, aod::McFwdTrackLabels>;
using MyMCFwdTrack = MyMCFwdTracks::iterator;

struct CreateResolutionMap {
  // Index used to set different options for Muon propagation
  enum class MuonExtrapolation : int {
    kToVertex = 0, // propagtion to vertex by default
    kToDCA = 1,
    kToRabs = 2,
  };
  using SMatrix55 = ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepSym<double, 5>>;
  using SMatrix5 = ROOT::Math::SVector<double, 5>;

  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
  Configurable<int> cfgEventGeneratorType{"cfgEventGeneratorType", -1, "if positive, select event generator type. i.e. gap or signal"};
  Configurable<int> cfgCentEstimator{"cfgCentEstimator", 2, "FT0M:0, FT0A:1, FT0C:2"};
  Configurable<bool> cfg_require_true_mc_collision_association{"cfg_require_true_mc_collision_association", true, "flag to require true mc collision association"};
  Configurable<bool> cfg_reject_fake_match_its_tpc{"cfg_reject_fake_match_its_tpc", false, "flag to reject fake match between ITS-TPC"};
  // Configurable<bool> cfg_reject_fake_match_its_tpc_tof{"cfg_reject_fake_match_its_tpc_tof", false, "flag to reject fake match between ITS-TPC-TOF"};
  Configurable<bool> cfg_reject_fake_match_mft_mch{"cfg_reject_fake_match_mft_mch", false, "flag to reject fake match between MFT-MCH"};

  ConfigurableAxis ConfPtGenBins{"ConfPtGenBins", {VARIABLE_WIDTH, 0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.30, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.40, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.50, 0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59, 0.60, 0.61, 0.62, 0.63, 0.64, 0.65, 0.66, 0.67, 0.68, 0.69, 0.70, 0.71, 0.72, 0.73, 0.74, 0.75, 0.76, 0.77, 0.78, 0.79, 0.80, 0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 0.87, 0.88, 0.89, 0.90, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 1.00, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00, 2.10, 2.20, 2.30, 2.40, 2.50, 2.60, 2.70, 2.80, 2.90, 3.00, 3.10, 3.20, 3.30, 3.40, 3.50, 3.60, 3.70, 3.80, 3.90, 4.00, 4.10, 4.20, 4.30, 4.40, 4.50, 4.60, 4.70, 4.80, 4.90, 5.00, 5.50, 6.00, 6.50, 7.00, 7.50, 8.00, 8.50, 9.00, 9.50, 10.00, 11.00, 12.00, 13.00, 14.00, 15.00, 16.00, 17.00, 18.00, 19.00, 20.00}, "gen. pT bins for output histograms"};
  ConfigurableAxis ConfCentBins{"ConfCentBins", {VARIABLE_WIDTH, 0, 10, 30, 50, 110}, "centrality (%) bins for output histograms"};

  ConfigurableAxis ConfEtaCBGenBins{"ConfEtaCBGenBins", {30, -1.5, +1.5}, "gen. eta bins at midrapidity for output histograms"};
  ConfigurableAxis ConfEtaFWDGenBins{"ConfEtaFWDGenBins", {40, -5.5, -1.5}, "gen. eta bins at forward rapidity for output histograms"};
  ConfigurableAxis ConfPhiGenBins{"ConfPhiGenBins", {72, 0, 2.f * M_PI}, "gen. eta bins at forward rapidity for output histograms"};

  ConfigurableAxis ConfRelDeltaPtBins{"ConfRelDeltaPtBins", {200, -1.f, +1.f}, "rel. dpt for output histograms"};
  ConfigurableAxis ConfDeltaEtaBins{"ConfDeltaEtaBins", {100, -0.1f, +0.1f}, "deta bins for output histograms"};
  ConfigurableAxis ConfDeltaPhiBins{"ConfDeltaPhiBins", {100, -0.1f, +0.1f}, "dphi bins for output histograms"};

  struct : ConfigurableGroup {
    std::string prefix = "electroncut_group";
    Configurable<float> cfg_min_pt_track{"cfg_min_pt_track", 0.01, "min pT for single track"};
    Configurable<float> cfg_min_eta_track{"cfg_min_eta_track", -1.5, "min eta for single track"};
    Configurable<float> cfg_max_eta_track{"cfg_max_eta_track", +1.5, "max eta for single track"};
    Configurable<int> cfg_min_ncluster_tpc{"cfg_min_ncluster_tpc", 0, "min ncluster tpc"};
    Configurable<int> cfg_min_ncluster_its{"cfg_min_ncluster_its", 4, "min ncluster its"};
    Configurable<int> cfg_min_ncluster_itsib{"cfg_min_ncluster_itsib", 1, "min ncluster itsib"};
    Configurable<int> cfg_min_ncrossedrows{"cfg_min_ncrossedrows", 80, "min ncrossed rows"};
    Configurable<float> cfg_max_chi2tpc{"cfg_max_chi2tpc", 4.0, "max chi2/NclsTPC"};
    Configurable<float> cfg_max_chi2its{"cfg_max_chi2its", 5.0, "max chi2/NclsITS"};
    Configurable<float> cfg_min_tpc_cr_findable_ratio{"cfg_min_tpc_cr_findable_ratio", 0.8, "min. TPC Ncr/Nf ratio"};
    Configurable<float> cfg_max_frac_shared_clusters_tpc{"cfg_max_frac_shared_clusters_tpc", 999.f, "max fraction of shared clusters in TPC"};
    Configurable<float> cfg_max_dcaxy{"cfg_max_dcaxy", 0.3, "max dca XY for single track in cm"};
    Configurable<float> cfg_max_dcaz{"cfg_max_dcaz", 0.3, "max dca Z for single track in cm"};
    Configurable<bool> cfg_require_itsib_1st{"cfg_require_itsib_1st", false, "flag to require ITS ib 1st hit"};
    Configurable<bool> cfg_enableTTCA{"enableTTCA", true, "Flag to enable or disable TTCA"};
  } electroncuts;

  struct : ConfigurableGroup {
    std::string prefix = "muoncut_group";
    Configurable<uint> cfg_track_type{"cfg_track_type", 3, "muon track type [0: MFT-MCH-MID, 3: MCH-MID]"};
    Configurable<float> cfg_min_pt_track{"cfg_min_pt_track", 0.01, "min pT for single track"};
    Configurable<float> cfg_min_eta_track{"cfg_min_eta_track", -5.5, "min eta for single track"};
    Configurable<float> cfg_max_eta_track{"cfg_max_eta_track", -1.5, "max eta for single track"};
    Configurable<int> cfg_min_ncluster_mft{"cfg_min_ncluster_mft", 5, "min ncluster MFT"};
    Configurable<int> cfg_min_ncluster_mch{"cfg_min_ncluster_mch", 5, "min ncluster MCH"};
    Configurable<float> cfg_max_chi2{"cfg_max_chi2", 1e+10, "max chi2/NclsTPC"};
    Configurable<float> cfg_max_matching_chi2_mftmch{"cfg_max_matching_chi2_mftmch", 1e+10, "max chi2 for MFT-MCH matching"};
    Configurable<float> cfg_max_matching_chi2_mchmid{"cfg_max_matching_chi2_mchmid", 1e+10, "max chi2 for MCH-MID matching"};
    Configurable<float> cfg_max_dcaxy{"cfg_max_dcaxy", 1e+10, "max dca XY for single track in cm"};
    Configurable<float> cfg_min_rabs{"cfg_min_rabs", 17.6, "min Radius at the absorber end"};
    Configurable<float> cfg_max_rabs{"cfg_max_rabs", 89.5, "max Radius at the absorber end"};
    // Configurable<bool> cfg_enableTTCA{"enableTTCA", true, "Flag to enable or disable TTCA"};
  } muoncuts;

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};

  o2::ccdb::CcdbApi ccdbApi;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::globaltracking::MatchGlobalFwd mMatching;
  int mRunNumber = 0;

  void init(o2::framework::InitContext&)
  {
    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);
    ccdbApi.init(ccdburl);

    const AxisSpec axis_cent{ConfCentBins, "centrality (%)"};
    const AxisSpec axis_pt_gen{ConfPtGenBins, "p_{T,l}^{gen} (GeV/c)"};
    const AxisSpec axis_eta_cb_gen{ConfEtaCBGenBins, "#eta_{l}^{gen}"};
    const AxisSpec axis_eta_fwd_gen{ConfEtaFWDGenBins, "#eta_{l}^{gen}"};
    const AxisSpec axis_phi_gen{ConfPhiGenBins, "#varphi_{l}^{gen} (rad.)"};
    const AxisSpec axis_dpt{ConfRelDeltaPtBins, "(p_{T,l}^{gen} - p_{T,l}^{rec})/p_{T,l}^{gen}"};
    const AxisSpec axis_deta{ConfDeltaEtaBins, "#eta_{l}^{gen} - #eta_{l}^{rec}"};
    const AxisSpec axis_dphi{ConfDeltaPhiBins, "#varphi_{l}^{gen} - #varphi_{l}^{rec} (rad.)"};
    const AxisSpec axis_charge_gen{3, -1.5, +1.5, "true sign"};

    registry.add("Event/hImpPar_Centrality", "true imapact parameter vs. estimated centrality;impact parameter (fm);centrality (%)", kTH2F, {{200, 0, 20}, {110, 0, 110}}, true);
    registry.add("Electron/Ptgen_RelDeltaPt", "resolution", kTH2F, {{axis_pt_gen}, {axis_dpt}}, true);
    registry.add("Electron/Ptgen_DeltaEta", "resolution", kTH2F, {{axis_pt_gen}, {axis_deta}}, true);
    registry.add("Electron/Ptgen_DeltaPhi_Pos", "resolution", kTH2F, {{axis_pt_gen}, {axis_dphi}}, true);
    registry.add("Electron/Ptgen_DeltaPhi_Neg", "resolution", kTH2F, {{axis_pt_gen}, {axis_dphi}}, true);
    registry.addClone("Electron/", "StandaloneMuon/");
    registry.addClone("Electron/", "GlobalMuon/");

    registry.add("Electron/hs_reso", "8D resolution positive", kTHnSparseF, {axis_cent, axis_pt_gen, axis_eta_cb_gen, axis_phi_gen, axis_charge_gen, axis_dpt, axis_deta, axis_dphi}, true);
    registry.add("StandaloneMuon/hs_reso", "8D resolution positive", kTHnSparseF, {axis_cent, axis_pt_gen, axis_eta_fwd_gen, axis_phi_gen, axis_charge_gen, axis_dpt, axis_deta, axis_dphi}, true);
    registry.add("GlobalMuon/hs_reso", "8D resolution positive", kTHnSparseF, {axis_cent, axis_pt_gen, axis_eta_fwd_gen, axis_phi_gen, axis_charge_gen, axis_dpt, axis_deta, axis_dphi}, true);
  }

  void initCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }

    mRunNumber = bc.runNumber();
    std::map<string, string> metadata;
    auto soreor = o2::ccdb::BasicCCDBManager::getRunDuration(ccdbApi, mRunNumber);
    auto ts = soreor.first;
    auto grpmag = ccdbApi.retrieveFromTFileAny<o2::parameters::GRPMagField>(grpmagPath, metadata, ts);
    o2::base::Propagator::initFieldFromGRP(grpmag);
    if (!o2::base::GeometryManager::isGeometryLoaded()) {
      ccdb->get<TGeoManager>(geoPath);
    }
    o2::mch::TrackExtrap::setField();
  }

  std::pair<int8_t, std::set<uint8_t>> itsRequirement_ibany = {1, {0, 1, 2}}; // any hits on 3 ITS ib layers.
  std::pair<int8_t, std::set<uint8_t>> itsRequirement_ib1st = {1, {0}};       // first hit on ITS ib layers.

  template <typename TTrack>
  bool checkTrack(TTrack const& track)
  {
    if (track.tpcChi2NCl() > electroncuts.cfg_max_chi2tpc) {
      return false;
    }

    if (track.itsChi2NCl() > electroncuts.cfg_max_chi2its) {
      return false;
    }

    if (track.itsNCls() < electroncuts.cfg_min_ncluster_its) {
      return false;
    }
    if (track.itsNClsInnerBarrel() < electroncuts.cfg_min_ncluster_itsib) {
      return false;
    }

    auto hits = std::count_if(itsRequirement_ibany.second.begin(), itsRequirement_ibany.second.end(), [&](auto&& requiredLayer) { return track.itsClusterMap() & (1 << requiredLayer); });
    if (hits < itsRequirement_ibany.first) {
      return false;
    }
    if (electroncuts.cfg_require_itsib_1st) {
      auto hit_ib1st = std::count_if(itsRequirement_ib1st.second.begin(), itsRequirement_ib1st.second.end(), [&](auto&& requiredLayer) { return track.itsClusterMap() & (1 << requiredLayer); });
      if (hit_ib1st < itsRequirement_ib1st.first) {
        return false;
      }
    }

    if (track.tpcNClsFound() < electroncuts.cfg_min_ncluster_tpc) {
      return false;
    }

    if (track.tpcNClsCrossedRows() < electroncuts.cfg_min_ncrossedrows) {
      return false;
    }

    if (track.tpcCrossedRowsOverFindableCls() < electroncuts.cfg_min_tpc_cr_findable_ratio) {
      return false;
    }

    if (track.tpcFractionSharedCls() > electroncuts.cfg_max_frac_shared_clusters_tpc) {
      return false;
    }

    return true;
  }

  template <typename T, typename C>
  o2::dataformats::GlobalFwdTrack PropagateMuon(T const& muon, C const& collision, const CreateResolutionMap::MuonExtrapolation endPoint)
  {
    double chi2 = muon.chi2();
    SMatrix5 tpars(muon.x(), muon.y(), muon.phi(), muon.tgl(), muon.signed1Pt());
    std::vector<float> v1{muon.cXX(), muon.cXY(), muon.cYY(), muon.cPhiX(), muon.cPhiY(),
                          muon.cPhiPhi(), muon.cTglX(), muon.cTglY(), muon.cTglPhi(), muon.cTglTgl(),
                          muon.c1PtX(), muon.c1PtY(), muon.c1PtPhi(), muon.c1PtTgl(), muon.c1Pt21Pt2()};
    SMatrix55 tcovs(v1.begin(), v1.end());
    o2::track::TrackParCovFwd fwdtrack{muon.z(), tpars, tcovs, chi2};
    o2::dataformats::GlobalFwdTrack propmuon;

    if (static_cast<int>(muon.trackType()) > 2) { // MCH-MID or MCH standalone
      o2::dataformats::GlobalFwdTrack track;
      track.setParameters(tpars);
      track.setZ(fwdtrack.getZ());
      track.setCovariances(tcovs);
      auto mchTrack = mMatching.FwdtoMCH(track);

      if (endPoint == CreateResolutionMap::MuonExtrapolation::kToVertex) {
        o2::mch::TrackExtrap::extrapToVertex(mchTrack, collision.posX(), collision.posY(), collision.posZ(), collision.covXX(), collision.covYY());
      }
      if (endPoint == CreateResolutionMap::MuonExtrapolation::kToDCA) {
        o2::mch::TrackExtrap::extrapToVertexWithoutBranson(mchTrack, collision.posZ());
      }
      if (endPoint == CreateResolutionMap::MuonExtrapolation::kToRabs) {
        o2::mch::TrackExtrap::extrapToZ(mchTrack, -505.);
      }

      auto proptrack = mMatching.MCHtoFwd(mchTrack);
      propmuon.setParameters(proptrack.getParameters());
      propmuon.setZ(proptrack.getZ());
      propmuon.setCovariances(proptrack.getCovariances());
    } else if (static_cast<int>(muon.trackType()) < 2) { // MFT-MCH-MID
      double centerMFT[3] = {0, 0, -61.4};
      o2::field::MagneticField* field = static_cast<o2::field::MagneticField*>(TGeoGlobalMagField::Instance()->GetField());
      auto Bz = field->getBz(centerMFT); // Get field at centre of MFT
      auto geoMan = o2::base::GeometryManager::meanMaterialBudget(muon.x(), muon.y(), muon.z(), collision.posX(), collision.posY(), collision.posZ());
      auto x2x0 = static_cast<float>(geoMan.meanX2X0);
      fwdtrack.propagateToVtxhelixWithMCS(collision.posZ(), {collision.posX(), collision.posY()}, {collision.covXX(), collision.covYY()}, Bz, x2x0);
      propmuon.setParameters(fwdtrack.getParameters());
      propmuon.setZ(fwdtrack.getZ());
      propmuon.setCovariances(fwdtrack.getCovariances());
    }

    v1.clear();
    v1.shrink_to_fit();

    return propmuon;
  }

  template <typename TMuon, typename TCollision>
  bool checkFwdTrack(TMuon const& muon, TCollision const& collision, const float centrality)
  {
    o2::dataformats::GlobalFwdTrack propmuonAtPV = PropagateMuon(muon, collision, CreateResolutionMap::MuonExtrapolation::kToVertex);
    float pt = propmuonAtPV.getPt();
    float eta = propmuonAtPV.getEta();
    float phi = propmuonAtPV.getPhi();

    if (pt < muoncuts.cfg_min_pt_track) {
      return false;
    }

    if (eta < muoncuts.cfg_min_eta_track || muoncuts.cfg_max_eta_track < eta) {
      return false;
    }

    o2::math_utils::bringTo02Pi(phi);
    if (phi < 0.f || 2.f * M_PI < phi) {
      return false;
    }

    float rAtAbsorberEnd = muon.rAtAbsorberEnd();
    if (muon.trackType() == static_cast<uint8_t>(o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack)) {
      o2::dataformats::GlobalFwdTrack propmuonAtRabs = PropagateMuon(muon, collision, CreateResolutionMap::MuonExtrapolation::kToRabs);
      float xAbs = propmuonAtRabs.getX();
      float yAbs = propmuonAtRabs.getY();
      rAtAbsorberEnd = std::sqrt(xAbs * xAbs + yAbs * yAbs); // Redo propagation only for muon tracks // propagation of MFT tracks alredy done in reconstruction
    }

    if (rAtAbsorberEnd < muoncuts.cfg_min_rabs || muoncuts.cfg_max_rabs < rAtAbsorberEnd) {
      return false;
    }

    if (rAtAbsorberEnd < 26.5) {
      if (muon.pDca() > 594.f) {
        return false;
      }
    } else {
      if (muon.pDca() > 324.f) {
        return false;
      }
    }

    if (muon.trackType() == static_cast<uint8_t>(o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack) && muon.chi2MatchMCHMFT() > muoncuts.cfg_max_matching_chi2_mftmch) {
      return false;
    }

    if (cfg_reject_fake_match_mft_mch && muon.trackType() == static_cast<uint8_t>(o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack) && o2::aod::pwgem::dilepton::utils::mcutil::hasFakeMatchMFTMCH(muon)) {
      return false;
    }

    auto mctrack = muon.template mcParticle_as<aod::McParticles>();
    if (muon.trackType() == static_cast<uint8_t>(o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack)) {
      registry.fill(HIST("StandaloneMuon/hs_reso"), centrality, mctrack.pt(), mctrack.eta(), mctrack.phi(), -mctrack.pdgCode() / 13, (mctrack.pt() - pt) / mctrack.pt(), mctrack.eta() - eta, mctrack.phi() - phi);
      registry.fill(HIST("StandaloneMuon/Ptgen_RelDeltaPt"), mctrack.pt(), (mctrack.pt() - pt) / mctrack.pt());
      registry.fill(HIST("StandaloneMuon/Ptgen_DeltaEta"), mctrack.pt(), mctrack.eta() - eta);
      if (mctrack.pdgCode() == -13) { // positive muon
        registry.fill(HIST("StandaloneMuon/Ptgen_DeltaPhi_Pos"), mctrack.pt(), mctrack.phi() - phi);
      } else if (mctrack.pdgCode() == 13) { // negative muon
        registry.fill(HIST("StandaloneMuon/Ptgen_DeltaPhi_Neg"), mctrack.pt(), mctrack.phi() - phi);
      }
    } else if (muon.trackType() == static_cast<uint8_t>(o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack)) {
      registry.fill(HIST("GlobalMuon/hs_reso"), centrality, mctrack.pt(), mctrack.eta(), mctrack.phi(), -mctrack.pdgCode() / 13, (mctrack.pt() - pt) / mctrack.pt(), mctrack.eta() - eta, mctrack.phi() - phi);
      registry.fill(HIST("GlobalMuon/Ptgen_RelDeltaPt"), mctrack.pt(), (mctrack.pt() - pt) / mctrack.pt());
      registry.fill(HIST("GlobalMuon/Ptgen_DeltaEta"), mctrack.pt(), mctrack.eta() - eta);
      if (mctrack.pdgCode() == -13) { // positive muon
        registry.fill(HIST("GlobalMuon/Ptgen_DeltaPhi_Pos"), mctrack.pt(), mctrack.phi() - phi);
      } else if (mctrack.pdgCode() == 13) { // negative muon
        registry.fill(HIST("GlobalMuon/Ptgen_DeltaPhi_Neg"), mctrack.pt(), mctrack.phi() - phi);
      }
    }
    return true;
  }

  SliceCache cache;
  Preslice<aod::Tracks> perCollision_mid = o2::aod::track::collisionId;
  Preslice<aod::FwdTracks> perCollision_fwd = o2::aod::fwdtrack::collisionId;

  Filter collisionFilter = o2::aod::evsel::sel8 == true && nabs(o2::aod::collision::posZ) < 10.f;
  using MyFilteredCollisions = soa::Filtered<MyCollisions>;
  using MyFilteredCollisionsCent = soa::Filtered<MyCollisionsCent>;

  Filter trackFilter_mid = o2::aod::track::pt > electroncuts.cfg_min_pt_track&& electroncuts.cfg_min_eta_track < o2::aod::track::eta&& o2::aod::track::eta < electroncuts.cfg_max_eta_track&& o2::aod::track::tpcChi2NCl < electroncuts.cfg_max_chi2tpc&& o2::aod::track::itsChi2NCl < electroncuts.cfg_max_chi2its&& ncheckbit(aod::track::v001::detectorMap, (uint8_t)o2::aod::track::ITS) == true && ncheckbit(aod::track::v001::detectorMap, (uint8_t)o2::aod::track::TPC) == true && nabs(o2::aod::track::dcaXY) < electroncuts.cfg_max_dcaxy&& nabs(o2::aod::track::dcaZ) < electroncuts.cfg_max_dcaz;
  Filter ttcaFilter_electron = ifnode(electroncuts.cfg_enableTTCA.node(), o2::aod::emprimaryelectron::isAssociatedToMPC == true || o2::aod::emprimaryelectron::isAssociatedToMPC == false, o2::aod::emprimaryelectron::isAssociatedToMPC == true);
  using MyFilteredMCTracks = soa::Filtered<MyMCTracks>;

  Partition<MyMCFwdTracks> sa_muons = o2::aod::fwdtrack::trackType == uint8_t(o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack); // MCH-MID
  Partition<MyMCFwdTracks> global_muons = o2::aod::fwdtrack::trackType == uint8_t(o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack); // MFT-MCH-MID

  template <typename TCollisions>
  void process(TCollisions const& collisions, aod::BCsWithTimestamps const&, MyFilteredMCTracks const& tracks, MyMCFwdTracks const&, aod::McCollisions const&, aod::McParticles const&)
  {
    for (auto& collision : collisions) {
      auto bc = collision.template foundBC_as<aod::BCsWithTimestamps>();
      initCCDB(bc);

      if (!collision.has_mcCollision()) {
        continue;
      }

      auto mccollision = collision.template mcCollision_as<aod::McCollisions>();
      if (cfgEventGeneratorType >= 0 && mccollision.getSubGeneratorId() != cfgEventGeneratorType) {
        continue;
      }

      float centrality = 105.f;
      if constexpr (std::is_same_v<std::decay_t<TCollisions>, MyFilteredCollisionsCent>) {
        centrality = std::array{collision.centFT0M(), collision.centFT0A(), collision.centFT0C()}[cfgCentEstimator];
      }

      registry.fill(HIST("Event/hImpPar_Centrality"), mccollision.impactParameter(), centrality);

      auto tracks_per_coll = tracks.sliceBy(perCollision_mid, collision.globalIndex());
      for (auto& track : tracks_per_coll) {
        if (!track.has_mcParticle()) {
          continue;
        }
        auto mctrack = track.template mcParticle_as<aod::McParticles>();
        if (cfg_require_true_mc_collision_association && mctrack.mcCollisionId() != collision.mcCollisionId()) {
          continue;
        }
        if (abs(mctrack.pdgCode()) != 11 || !(mctrack.isPhysicalPrimary() || mctrack.producedByGenerator())) {
          continue;
        }
        if (!checkTrack(track)) {
          continue;
        }

        if (cfg_reject_fake_match_its_tpc && o2::aod::pwgem::dilepton::utils::mcutil::hasFakeMatchITSTPC(track)) {
          continue;
        }

        registry.fill(HIST("Electron/hs_reso"), centrality, mctrack.pt(), mctrack.eta(), mctrack.phi(), -mctrack.pdgCode() / 11, (mctrack.pt() - track.pt()) / mctrack.pt(), mctrack.eta() - track.eta(), mctrack.phi() - track.phi());
        registry.fill(HIST("Electron/Ptgen_RelDeltaPt"), mctrack.pt(), (mctrack.pt() - track.pt()) / mctrack.pt());
        registry.fill(HIST("Electron/Ptgen_DeltaEta"), mctrack.pt(), mctrack.eta() - track.eta());
        if (mctrack.pdgCode() == -11) { // positron
          registry.fill(HIST("Electron/Ptgen_DeltaPhi_Pos"), mctrack.pt(), mctrack.phi() - track.phi());
        } else if (mctrack.pdgCode() == 11) { // electron
          registry.fill(HIST("Electron/Ptgen_DeltaPhi_Neg"), mctrack.pt(), mctrack.phi() - track.phi());
        }

      } // end of track loop

      auto sa_muons_per_coll = sa_muons->sliceByCached(o2::aod::fwdtrack::collisionId, collision.globalIndex(), cache);
      auto global_muons_per_coll = global_muons->sliceByCached(o2::aod::fwdtrack::collisionId, collision.globalIndex(), cache);

      for (auto& muon : sa_muons_per_coll) {
        if (!muon.has_mcParticle()) {
          continue;
        }
        auto mctrack = muon.template mcParticle_as<aod::McParticles>();
        if (cfg_require_true_mc_collision_association && mctrack.mcCollisionId() != collision.mcCollisionId()) {
          continue;
        }
        if (abs(mctrack.pdgCode()) != 13 || !(mctrack.isPhysicalPrimary() || mctrack.producedByGenerator())) {
          continue;
        }
        if (!checkFwdTrack(muon, collision, centrality)) {
          continue;
        }
      } // end of standalone muon loop

      for (auto& muon : global_muons_per_coll) {
        if (!muon.has_mcParticle()) {
          continue;
        }
        auto mctrack = muon.template mcParticle_as<aod::McParticles>();
        if (cfg_require_true_mc_collision_association && mctrack.mcCollisionId() != collision.mcCollisionId()) {
          continue;
        }
        if (abs(mctrack.pdgCode()) != 13 || !(mctrack.isPhysicalPrimary() || mctrack.producedByGenerator())) {
          continue;
        }
        if (!checkFwdTrack(muon, collision, centrality)) {
          continue;
        }

      } // end of global muon loop

    } // end of collision loop
  }
  PROCESS_SWITCH_FULL(CreateResolutionMap, process<MyFilteredCollisionsCent>, processWithCent, "create resolution map wit centrality", true);
  PROCESS_SWITCH_FULL(CreateResolutionMap, process<MyFilteredCollisions>, processWithoutCent, "create resolution map without centrality", false);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<CreateResolutionMap>(cfgc, TaskName{"create-resolution-map"})};
}
