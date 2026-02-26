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
// This code runs loop over dalitz ee table for dalitz QC.
//    Please write to: daiki.sekihata@cern.ch

#include "PWGEM/Dilepton/Utils/PairUtilities.h"
#include "PWGEM/PhotonMeson/Core/DalitzEECut.h"
#include "PWGEM/PhotonMeson/Core/EMPhotonEventCut.h"
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGEM/PhotonMeson/Utils/EventHistograms.h"

#include "Common/Core/RecoDecay.h"

#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DetectorsBase/GeometryManager.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include "Math/Vector4D.h"
#include "TString.h"

#include <string>
#include <vector>

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::aod::pwgem::dilepton::utils::emtrackutil;

using MyCollisions = soa::Join<aod::EMEvents_004, aod::EMEventsAlias, aod::EMEventsMult_000, aod::EMEventsCent_000>;
using MyCollision = MyCollisions::iterator;

using MyTracks = soa::Join<aod::EMPrimaryElectronsFromDalitz, aod::EMPrimaryElectronDaEMEventIds>;
using MyTrack = MyTracks::iterator;

struct DalitzEEQC {
  // Configurables
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<bool> skipGRPOquery{"skipGRPOquery", true, "skip grpo query"};
  Configurable<float> d_bz_input{"d_bz_input", -999, "bz field in kG, -999 is automatic"};

  Configurable<int> cfgCentEstimator{"cfgCentEstimator", 2, "FT0M:0, FT0A:1, FT0C:2"};
  Configurable<float> cfgCentMin{"cfgCentMin", 0, "min. centrality"};
  Configurable<float> cfgCentMax{"cfgCentMax", 999.f, "max. centrality"};
  Configurable<float> maxY{"maxY", 0.9, "maximum rapidity for reconstructed particles"};

  EMPhotonEventCut fEMEventCut;
  struct : ConfigurableGroup {
    std::string prefix = "eventcut_group";
    Configurable<float> cfgZvtxMax{"cfgZvtxMax", 10.f, "max. Zvtx"};
    Configurable<bool> cfgRequireSel8{"cfgRequireSel8", true, "require sel8 in event cut"};
    Configurable<bool> cfgRequireFT0AND{"cfgRequireFT0AND", true, "require FT0AND in event cut"};
    Configurable<bool> cfgRequireNoTFB{"cfgRequireNoTFB", false, "require No time frame border in event cut"};
    Configurable<bool> cfgRequireNoITSROFB{"cfgRequireNoITSROFB", false, "require no ITS readout frame border in event cut"};
    Configurable<bool> cfgRequireNoSameBunchPileup{"cfgRequireNoSameBunchPileup", false, "require no same bunch pileup in event cut"};
    Configurable<bool> cfgRequireVertexITSTPC{"cfgRequireVertexITSTPC", false, "require Vertex ITSTPC in event cut"}; // ITS-TPC matched track contributes PV.
    Configurable<bool> cfgRequireGoodZvtxFT0vsPV{"cfgRequireGoodZvtxFT0vsPV", false, "require good Zvtx between FT0 vs. PV in event cut"};
    Configurable<int> cfgOccupancyMin{"cfgOccupancyMin", -1, "min. occupancy"};
    Configurable<int> cfgOccupancyMax{"cfgOccupancyMax", 1000000000, "max. occupancy"};
  } eventcuts;

  DalitzEECut fDileptonCut;
  struct : ConfigurableGroup {
    std::string prefix = "dileptoncut_group";
    Configurable<float> cfg_min_mass{"cfg_min_mass", 0.0, "min mass"};
    Configurable<float> cfg_max_mass{"cfg_max_mass", 0.02, "max mass"};
    Configurable<bool> cfg_apply_phiv{"cfg_apply_phiv", true, "flag to apply phiv cut"};
    Configurable<bool> cfg_require_itsib_any{"cfg_require_itsib_any", true, "flag to require ITS ib any hits"};
    Configurable<bool> cfg_require_itsib_1st{"cfg_require_itsib_1st", false, "flag to require ITS ib 1st hit"};
    Configurable<float> cfg_phiv_slope{"cfg_phiv_slope", 0.0185, "slope for m vs. phiv"};
    Configurable<float> cfg_phiv_intercept{"cfg_phiv_intercept", -0.0280, "intercept for m vs. phiv"};

    Configurable<float> cfg_min_pt_track{"cfg_min_pt_track", 0.1, "min pT for single track"};
    Configurable<float> cfg_max_eta_track{"cfg_max_eta_track", 0.9, "max eta for single track"};
    Configurable<int> cfg_min_ncluster_tpc{"cfg_min_ncluster_tpc", 0, "min ncluster tpc"};
    Configurable<int> cfg_min_ncluster_its{"cfg_min_ncluster_its", 5, "min ncluster its"};
    Configurable<int> cfg_min_ncrossedrows{"cfg_min_ncrossedrows", 70, "min ncrossed rows"};
    Configurable<float> cfg_max_frac_shared_clusters_tpc{"cfg_max_frac_shared_clusters_tpc", 999.f, "max fraction of shared clusters in TPC"};
    Configurable<float> cfg_max_chi2tpc{"cfg_max_chi2tpc", 4.0, "max chi2/NclsTPC"};
    Configurable<float> cfg_max_chi2its{"cfg_max_chi2its", 5.0, "max chi2/NclsITS"};
    Configurable<float> cfg_max_dcaxy{"cfg_max_dcaxy", 1.0, "max dca XY for single track in cm"};
    Configurable<float> cfg_max_dcaz{"cfg_max_dcaz", 1.0, "max dca Z for single track in cm"};
    Configurable<float> cfg_max_dca3dsigma_track{"cfg_max_dca3dsigma_track", 1.5, "max DCA 3D in sigma"};
    Configurable<bool> includeITSsa{"includeITSsa", false, "Flag to enable ITSsa tracks"};
    Configurable<float> cfg_max_pt_track_ITSsa{"cfg_max_pt_track_ITSsa", 0.15, "max pt for ITSsa tracks"};

    Configurable<int> cfg_pid_scheme{"cfg_pid_scheme", static_cast<int>(DalitzEECut::PIDSchemes::kTOFif), "pid scheme [kTOFif : 0, kTPConly : 1]"};
    Configurable<float> cfg_min_TPCNsigmaEl{"cfg_min_TPCNsigmaEl", -2.0, "min. TPC n sigma for electron inclusion"};
    Configurable<float> cfg_max_TPCNsigmaEl{"cfg_max_TPCNsigmaEl", +3.0, "max. TPC n sigma for electron inclusion"};
    Configurable<float> cfg_min_TPCNsigmaPi{"cfg_min_TPCNsigmaPi", -0.0, "min. TPC n sigma for pion exclusion"};
    Configurable<float> cfg_max_TPCNsigmaPi{"cfg_max_TPCNsigmaPi", +0.0, "max. TPC n sigma for pion exclusion"};
    Configurable<float> cfg_min_TOFNsigmaEl{"cfg_min_TOFNsigmaEl", -3.0, "min. TOF n sigma for electron inclusion"};
    Configurable<float> cfg_max_TOFNsigmaEl{"cfg_max_TOFNsigmaEl", +3.0, "max. TOF n sigma for electron inclusion"};
  } dileptoncuts;

  o2::ccdb::CcdbApi ccdbApi;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  int mRunNumber;
  float d_bz;

  HistogramRegistry fRegistry{"output", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};
  static constexpr std::string_view event_cut_types[2] = {"before/", "after/"};

  void init(InitContext& /*context*/)
  {
    DefineEMEventCut();
    DefineDileptonCut();
    addhistograms();

    mRunNumber = 0;
    d_bz = 0;

    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);
  }

  template <typename TCollision>
  void initCCDB(TCollision const& collision)
  {
    if (mRunNumber == collision.runNumber()) {
      return;
    }

    // In case override, don't proceed, please - no CCDB access required
    if (d_bz_input > -990) {
      d_bz = d_bz_input;
      o2::parameters::GRPMagField grpmag;
      if (std::fabs(d_bz) > 1e-5) {
        grpmag.setL3Current(30000.f / (d_bz / 5.0f));
      }
      mRunNumber = collision.runNumber();
      return;
    }

    auto run3grp_timestamp = collision.timestamp();
    o2::parameters::GRPObject* grpo = 0x0;
    o2::parameters::GRPMagField* grpmag = 0x0;
    if (!skipGRPOquery)
      grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>(grpPath, run3grp_timestamp);
    if (grpo) {
      // Fetch magnetic field from ccdb for current collision
      d_bz = grpo->getNominalL3Field();
      LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << d_bz << " kZG";
    } else {
      grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, run3grp_timestamp);
      if (!grpmag) {
        LOG(fatal) << "Got nullptr from CCDB for path " << grpmagPath << " of object GRPMagField and " << grpPath << " of object GRPObject for timestamp " << run3grp_timestamp;
      }
      // Fetch magnetic field from ccdb for current collision
      d_bz = std::lround(5.f * grpmag->getL3Current() / 30000.f);
      LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << d_bz << " kZG";
    }
    mRunNumber = collision.runNumber();
  }

  ~DalitzEEQC() {}

  void addhistograms()
  {
    // event info
    o2::aod::pwgem::photonmeson::utils::eventhistogram::addEventHistograms(&fRegistry);

    // pair info
    const AxisSpec axis_mass{200, 0, 0.2, "m_{ee} (GeV/c^{2})"};
    const AxisSpec axis_pt{200, 0, 2, "p_{T,ee} (GeV/c)"};

    fRegistry.add("Pair/same/hMvsPt", "m_{ee} vs. p_{T,ee};m_{ee} (GeV/c^{2});p_{T,ee} (GeV/c)", kTH2F, {axis_mass, axis_pt}, true);
    fRegistry.add("Pair/same/hMvsPhiV", "m_{ee} vs. #varphi_{V};#varphi (rad.);m_{ee} (GeV/c^{2})", kTH2F, {{90, 0, M_PI}, {100, 0.0f, 0.1f}}, true);

    fRegistry.add("Track/hPt", "pT;p_{T} (GeV/c)", kTH1F, {{1000, 0.0f, 10}}, false);
    fRegistry.add("Track/hQoverPt", "q/pT;q/p_{T} (GeV/c)^{-1}", kTH1F, {{400, -20, 20}}, false);
    fRegistry.add("Track/hEtaPhi", "#eta vs. #varphi;#varphi (rad.);#eta", kTH2F, {{180, 0, 2 * M_PI}, {40, -2.0f, 2.0f}}, false);
    fRegistry.add("Track/hDCAxyz", "DCA xy vs. z;DCA_{xy} (cm);DCA_{z} (cm)", kTH2F, {{200, -0.1f, 0.1f}, {200, -0.1f, 0.1f}}, false);
    fRegistry.add("Track/hDCAxyzSigma", "DCA xy vs. z;DCA_{xy} (#sigma);DCA_{z} (#sigma)", kTH2F, {{200, -10.0f, 10.0f}, {200, -10.0f, 10.0f}}, false);
    fRegistry.add("Track/hDCAxyRes_Pt", "DCA_{xy} resolution vs. pT;p_{T} (GeV/c);DCA_{xy} resolution (#mum)", kTH2F, {{200, 0, 10}, {500, 0., 500}}, false);
    fRegistry.add("Track/hDCAzRes_Pt", "DCA_{z} resolution vs. pT;p_{T} (GeV/c);DCA_{z} resolution (#mum)", kTH2F, {{200, 0, 10}, {500, 0., 500}}, false);
    fRegistry.add("Track/hNclsTPC", "number of TPC clusters", kTH1F, {{161, -0.5, 160.5}}, false);
    fRegistry.add("Track/hNcrTPC", "number of TPC crossed rows", kTH1F, {{161, -0.5, 160.5}}, false);
    fRegistry.add("Track/hChi2TPC", "chi2/number of TPC clusters", kTH1F, {{100, 0, 10}}, false);
    fRegistry.add("Track/hTPCNcr2Nf", "TPC Ncr/Nfindable", kTH1F, {{200, 0, 2}}, false);
    fRegistry.add("Track/hTPCNcls2Nf", "TPC Ncls/Nfindable", kTH1F, {{200, 0, 2}}, false);
    fRegistry.add("Track/hTPCNclsShared", "TPC Ncls shared/Ncls;p_{T} (GeV/c);N_{cls}^{shared}/N_{cls} in TPC", kTH2F, {{1000, 0, 10}, {100, 0, 1}}, false);

    fRegistry.add("Track/hNclsITS", "number of ITS clusters", kTH1F, {{8, -0.5, 7.5}}, false);
    fRegistry.add("Track/hChi2ITS", "chi2/number of ITS clusters", kTH1F, {{100, 0, 10}}, false);
    fRegistry.add("Track/hITSClusterMap", "ITS cluster map", kTH1F, {{128, -0.5, 127.5}}, false);
    fRegistry.add("Track/hMeanClusterSizeITS", "mean cluster size ITS;p_{pv} (GeV/c);<cluster size> on ITS #times cos(#lambda)", kTH2F, {{1000, 0, 10}, {150, 0, 15}}, false);

    fRegistry.add("Track/hTPCdEdx", "TPC dE/dx;p_{in} (GeV/c);TPC dE/dx (a.u.)", kTH2F, {{1000, 0, 10}, {200, 0, 200}}, false);
    fRegistry.add("Track/hTPCNsigmaEl", "TPC n sigma el;p_{in} (GeV/c);n #sigma_{e}^{TPC}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    fRegistry.add("Track/hTPCNsigmaPi", "TPC n sigma pi;p_{in} (GeV/c);n #sigma_{#pi}^{TPC}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);

    fRegistry.add("Track/hChi2TOF", "chi2 of TOF", kTH1F, {{100, 0, 10}}, false);
    fRegistry.add("Track/hTOFbeta", "TOF beta;p_{pv} (GeV/c);#beta", kTH2F, {{1000, 0, 10}, {240, 0, 1.2}}, false);
    fRegistry.add("Track/hTOFNsigmaEl", "TOF n sigma el;p_{pv} (GeV/c);n #sigma_{e}^{TOF}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
  }

  void DefineEMEventCut()
  {
    fEMEventCut = EMPhotonEventCut("fEMEventCut", "fEMEventCut");
    fEMEventCut.SetRequireSel8(eventcuts.cfgRequireSel8);
    fEMEventCut.SetRequireFT0AND(eventcuts.cfgRequireFT0AND);
    fEMEventCut.SetZvtxRange(-eventcuts.cfgZvtxMax, +eventcuts.cfgZvtxMax);
    fEMEventCut.SetRequireNoTFB(eventcuts.cfgRequireNoTFB);
    fEMEventCut.SetRequireNoITSROFB(eventcuts.cfgRequireNoITSROFB);
    fEMEventCut.SetRequireNoSameBunchPileup(eventcuts.cfgRequireNoSameBunchPileup);
    fEMEventCut.SetRequireVertexITSTPC(eventcuts.cfgRequireVertexITSTPC);
    fEMEventCut.SetRequireGoodZvtxFT0vsPV(eventcuts.cfgRequireGoodZvtxFT0vsPV);
  }

  void DefineDileptonCut()
  {
    fDileptonCut = DalitzEECut("fDileptonCut", "fDileptonCut");

    // for pair
    fDileptonCut.SetMeeRange(dileptoncuts.cfg_min_mass, dileptoncuts.cfg_max_mass);
    fDileptonCut.SetMaxPhivPairMeeDep([&](float mll) { return (mll - dileptoncuts.cfg_phiv_intercept) / dileptoncuts.cfg_phiv_slope; });
    fDileptonCut.ApplyPhiV(dileptoncuts.cfg_apply_phiv);
    fDileptonCut.RequireITSibAny(dileptoncuts.cfg_require_itsib_any);
    fDileptonCut.RequireITSib1st(dileptoncuts.cfg_require_itsib_1st);

    // for track
    fDileptonCut.SetTrackPtRange(dileptoncuts.cfg_min_pt_track, 1e+10f);
    fDileptonCut.SetTrackEtaRange(-dileptoncuts.cfg_max_eta_track, +dileptoncuts.cfg_max_eta_track);
    fDileptonCut.SetMinNClustersTPC(dileptoncuts.cfg_min_ncluster_tpc);
    fDileptonCut.SetMinNCrossedRowsTPC(dileptoncuts.cfg_min_ncrossedrows);
    fDileptonCut.SetMinNCrossedRowsOverFindableClustersTPC(0.8);
    fDileptonCut.SetMaxFracSharedClustersTPC(dileptoncuts.cfg_max_frac_shared_clusters_tpc);
    fDileptonCut.SetChi2PerClusterTPC(0.0, dileptoncuts.cfg_max_chi2tpc);
    fDileptonCut.SetChi2PerClusterITS(0.0, dileptoncuts.cfg_max_chi2its);
    fDileptonCut.SetNClustersITS(dileptoncuts.cfg_min_ncluster_its, 7);
    fDileptonCut.SetMaxDcaXY(dileptoncuts.cfg_max_dcaxy);
    fDileptonCut.SetMaxDcaZ(dileptoncuts.cfg_max_dcaz);
    fDileptonCut.SetTrackDca3DRange(0.f, dileptoncuts.cfg_max_dca3dsigma_track); // in sigma
    fDileptonCut.IncludeITSsa(dileptoncuts.includeITSsa, dileptoncuts.cfg_max_pt_track_ITSsa);

    // for eID
    fDileptonCut.SetPIDScheme(dileptoncuts.cfg_pid_scheme);
    fDileptonCut.SetTPCNsigmaElRange(dileptoncuts.cfg_min_TPCNsigmaEl, dileptoncuts.cfg_max_TPCNsigmaEl);
    fDileptonCut.SetTPCNsigmaPiRange(dileptoncuts.cfg_min_TPCNsigmaPi, dileptoncuts.cfg_max_TPCNsigmaPi);
    fDileptonCut.SetTOFNsigmaElRange(dileptoncuts.cfg_min_TOFNsigmaEl, dileptoncuts.cfg_max_TOFNsigmaEl);
  }

  template <typename TCollision, typename TTrack1, typename TTrack2>
  bool fillPairInfo(TCollision const&, TTrack1 const& t1, TTrack2 const& t2)
  {
    if (t1.trackId() == t2.trackId()) { // this is protection against pairing identical 2 tracks. This happens, when TTCA is used. TTCA can assign a track to several possible collisions.
      return false;
    }

    if (!fDileptonCut.IsSelectedTrack(t1) || !fDileptonCut.IsSelectedTrack(t2)) {
      return false;
    }

    if (!fDileptonCut.IsSelectedPair(t1, t2, d_bz)) {
      return false;
    }

    ROOT::Math::PtEtaPhiMVector v1(t1.pt(), t1.eta(), t1.phi(), o2::constants::physics::MassElectron);
    ROOT::Math::PtEtaPhiMVector v2(t2.pt(), t2.eta(), t2.phi(), o2::constants::physics::MassElectron);
    ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
    if (std::fabs(v12.Rapidity()) > maxY) {
      return false;
    }
    float phiv = o2::aod::pwgem::dilepton::utils::pairutil::getPhivPair(t1.px(), t1.py(), t1.pz(), t2.px(), t2.py(), t2.pz(), t1.sign(), t2.sign(), d_bz);

    if (t1.sign() * t2.sign() < 0) { // ULS
      fRegistry.fill(HIST("Pair/same/hMvsPt"), v12.M(), v12.Pt());
      fRegistry.fill(HIST("Pair/same/hMvsPhiV"), phiv, v12.M());
    }

    if (std::find(used_trackIds.begin(), used_trackIds.end(), t1.globalIndex()) == used_trackIds.end()) {
      used_trackIds.emplace_back(t1.globalIndex());
      fillTrackInfo(t1);
    }
    if (std::find(used_trackIds.begin(), used_trackIds.end(), t2.globalIndex()) == used_trackIds.end()) {
      used_trackIds.emplace_back(t2.globalIndex());
      fillTrackInfo(t2);
    }
    return true;
  }

  template <typename TTrack>
  void fillTrackInfo(TTrack const& track)
  {
    fRegistry.fill(HIST("Track/hPt"), track.pt());
    fRegistry.fill(HIST("Track/hQoverPt"), track.sign() / track.pt());
    fRegistry.fill(HIST("Track/hEtaPhi"), track.phi(), track.eta());
    fRegistry.fill(HIST("Track/hDCAxyz"), track.dcaXY(), track.dcaZ());
    fRegistry.fill(HIST("Track/hDCAxyzSigma"), track.dcaXY() / std::sqrt(track.cYY()), track.dcaZ() / std::sqrt(track.cZZ()));
    fRegistry.fill(HIST("Track/hDCAxyRes_Pt"), track.pt(), std::sqrt(track.cYY()) * 1e+4); // convert cm to um
    fRegistry.fill(HIST("Track/hDCAzRes_Pt"), track.pt(), std::sqrt(track.cZZ()) * 1e+4);  // convert cm to um
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
    fRegistry.fill(HIST("Track/hTPCNsigmaPi"), track.tpcInnerParam(), track.tpcNSigmaPi());
    fRegistry.fill(HIST("Track/hTOFNsigmaEl"), track.p(), track.tofNSigmaEl());
    fRegistry.fill(HIST("Track/hTOFbeta"), track.p(), track.beta());
    fRegistry.fill(HIST("Track/hMeanClusterSizeITS"), track.p(), track.meanClusterSizeITS() * std::cos(std::atan(track.tgl())));
  }

  Filter collisionFilter_centrality = (cfgCentMin < o2::aod::cent::centFT0M && o2::aod::cent::centFT0M < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0A && o2::aod::cent::centFT0A < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0C && o2::aod::cent::centFT0C < cfgCentMax);
  using FilteredMyCollisions = soa::Filtered<MyCollisions>;

  SliceCache cache;
  Preslice<MyTracks> perCollision_track = aod::emprimaryelectronda::emphotoneventId;
  Filter trackFilter = dileptoncuts.cfg_min_pt_track < o2::aod::track::pt && nabs(o2::aod::track::eta) < dileptoncuts.cfg_max_eta_track && o2::aod::track::tpcChi2NCl < dileptoncuts.cfg_max_chi2tpc && o2::aod::track::itsChi2NCl < dileptoncuts.cfg_max_chi2its && nabs(o2::aod::track::dcaXY) < dileptoncuts.cfg_max_dcaxy && nabs(o2::aod::track::dcaZ) < dileptoncuts.cfg_max_dcaz;
  Filter pidFilter = dileptoncuts.cfg_min_TPCNsigmaEl < o2::aod::pidtpc::tpcNSigmaEl && o2::aod::pidtpc::tpcNSigmaEl < dileptoncuts.cfg_max_TPCNsigmaEl && (o2::aod::pidtpc::tpcNSigmaPi < dileptoncuts.cfg_min_TPCNsigmaPi || dileptoncuts.cfg_max_TPCNsigmaPi < o2::aod::pidtpc::tpcNSigmaPi);
  using FilteredMyTracks = soa::Filtered<MyTracks>;
  Partition<FilteredMyTracks> posTracks = o2::aod::emprimaryelectron::sign > int8_t(0);
  Partition<FilteredMyTracks> negTracks = o2::aod::emprimaryelectron::sign < int8_t(0);

  std::vector<int> used_trackIds;
  void processQC(FilteredMyCollisions const& collisions, FilteredMyTracks const& /*tracks*/)
  {
    for (auto& collision : collisions) {
      initCCDB(collision);
      const float centralities[3] = {collision.centFT0M(), collision.centFT0A(), collision.centFT0C()};
      if (centralities[cfgCentEstimator] < cfgCentMin || cfgCentMax < centralities[cfgCentEstimator]) {
        continue;
      }

      o2::aod::pwgem::photonmeson::utils::eventhistogram::fillEventInfo<0>(&fRegistry, collision);
      if (!fEMEventCut.IsSelected(collision)) {
        continue;
      }
      if (!(eventcuts.cfgOccupancyMin <= collision.trackOccupancyInTimeRange() && collision.trackOccupancyInTimeRange() < eventcuts.cfgOccupancyMax)) {
        continue;
      }
      o2::aod::pwgem::photonmeson::utils::eventhistogram::fillEventInfo<1>(&fRegistry, collision);
      fRegistry.fill(HIST("Event/before/hCollisionCounter"), 12.0); // accepted
      fRegistry.fill(HIST("Event/after/hCollisionCounter"), 12.0);  // accepted

      auto posTracks_per_coll = posTracks->sliceByCached(o2::aod::emprimaryelectronda::emphotoneventId, collision.globalIndex(), cache);
      auto negTracks_per_coll = negTracks->sliceByCached(o2::aod::emprimaryelectronda::emphotoneventId, collision.globalIndex(), cache);
      // LOGF(info, "centrality = %f , posTracks_per_coll.size() = %d, negTracks_per_coll.size() = %d", centralities[cfgCentEstimator], posTracks_per_coll.size(), negTracks_per_coll.size());

      for (auto& [pos, ele] : combinations(CombinationsFullIndexPolicy(posTracks_per_coll, negTracks_per_coll))) { // ULS
        fillPairInfo(collision, pos, ele);
      }
    } // end of collision loop

    used_trackIds.clear();
    used_trackIds.shrink_to_fit();

  } // end of process
  PROCESS_SWITCH(DalitzEEQC, processQC, "run Dalitz QC", true);

  void processDummy(MyCollisions const&) {}
  PROCESS_SWITCH(DalitzEEQC, processDummy, "Dummy function", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<DalitzEEQC>(cfgc, TaskName{"dalitz-ee-qc"})};
}
