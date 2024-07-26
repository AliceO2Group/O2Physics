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
// Analysis task for dimuon in MC.
//    Please write to: daiki.sekihata@cern.ch

#include "TString.h"
#include "Math/Vector4D.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"

#include "CCDB/BasicCCDBManager.h"
#include "Common/Core/RecoDecay.h"

#include "PWGEM/Dilepton/DataModel/dileptonTables.h"
#include "PWGEM/Dilepton/Core/DimuonCut.h"
#include "PWGEM/Dilepton/Core/EMEventCut.h"
#include "PWGEM/Dilepton/Utils/MCUtilities.h"
#include "PWGEM/Dilepton/Utils/EventHistograms.h"
#include "PWGEM/Dilepton/Utils/EMTrackUtilities.h"
#include "PWGEM/Dilepton/Utils/PairUtilities.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::aod::pwgem::dilepton::utils::mcutil;
using namespace o2::aod::pwgem::dilepton::utils::emtrackutil;

using MyCollisions = soa::Join<aod::EMEvents, aod::EMEventsMult, aod::EMEventsCent, aod::EMMCEventLabels>;
using MyCollision = MyCollisions::iterator;

using MyMCTracks = soa::Join<aod::EMPrimaryMuons, aod::EMPrimaryMuonsCov, aod::EMPrimaryMuonEMEventIds, aod::EMPrimaryMuonMCLabels>;
using MyMCTrack = MyMCTracks::iterator;

struct singleMuonQCMC {

  // Configurables
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};

  Configurable<int> cfgEventGeneratorType{"cfgEventGeneratorType", -1, "if positive, select event generator type. i.e. gap or signal"};
  Configurable<int> cfgCentEstimator{"cfgCentEstimator", 2, "FT0M:0, FT0A:1, FT0C:2"};
  Configurable<float> cfgCentMin{"cfgCentMin", 0, "min. centrality"};
  Configurable<float> cfgCentMax{"cfgCentMax", 999.f, "max. centrality"};

  ConfigurableAxis ConfPtmuBins{"ConfPtmuBins", {VARIABLE_WIDTH, 0.00, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00, 2.10, 2.20, 2.30, 2.40, 2.50, 2.60, 2.70, 2.80, 2.90, 3.00, 3.10, 3.20, 3.30, 3.40, 3.50, 3.60, 3.70, 3.80, 3.90, 4.00, 4.10, 4.20, 4.30, 4.40, 4.50, 4.60, 4.70, 4.80, 4.90, 5.00, 5.50, 6.00, 6.50, 7.00, 7.50, 8.00, 8.50, 9.00, 9.50, 10.00}, "pTmu bins for output histograms"};

  EMEventCut fEMEventCut;
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

  DimuonCut fDimuonCut;
  struct : ConfigurableGroup {
    std::string prefix = "dimuoncut_group";
    Configurable<float> cfg_min_mass{"cfg_min_mass", 0.0, "min mass"};
    Configurable<float> cfg_max_mass{"cfg_max_mass", 1e+10, "max mass"};
    Configurable<float> cfg_min_pair_pt{"cfg_min_pair_pt", 0.0, "min pair pt"};
    Configurable<float> cfg_max_pair_pt{"cfg_max_pair_pt", 1e+10, "max pair pt"};
    Configurable<float> cfg_min_pair_dcaxy{"cfg_min_pair_dcaxy", 0.0, "min pair dca3d in sigma"};
    Configurable<float> cfg_max_pair_dcaxy{"cfg_max_pair_dcaxy", 1e+10, "max pair dca3d in sigma"};

    Configurable<uint8_t> cfg_track_type{"cfg_track_type", 3, "muon track type [0: MFT-MCH-MID, 3: MCH-MID]"};
    Configurable<float> cfg_min_pt_track{"cfg_min_pt_track", 0.1, "min pT for single track"};
    Configurable<float> cfg_min_eta_track{"cfg_min_eta_track", -4.0, "min eta for single track"};
    Configurable<float> cfg_max_eta_track{"cfg_max_eta_track", -2.5, "max eta for single track"};
    Configurable<int> cfg_min_ncluster_mft{"cfg_min_ncluster_mft", 5, "min ncluster MFT"};
    Configurable<int> cfg_min_ncluster_mch{"cfg_min_ncluster_mch", 5, "min ncluster MCH"};
    Configurable<float> cfg_max_chi2{"cfg_max_chi2", 1e+10, "max chi2/NclsTPC"};
    Configurable<float> cfg_max_matching_chi2_mftmch{"cfg_max_matching_chi2_mftmch", 1e+10, "max chi2 for MFT-MCH matching"};
    Configurable<float> cfg_max_matching_chi2_mchmid{"cfg_max_matching_chi2_mchmid", 1e+10, "max chi2 for MCH-MID matching"};
    Configurable<float> cfg_max_dcaxy{"cfg_max_dcaxy", 1e+10, "max dca XY for single track in cm"};
    Configurable<float> cfg_min_rabs{"cfg_min_rabs", 17.6, "min Radius at the absorber end"};
    Configurable<float> cfg_max_rabs{"cfg_max_rabs", 89.5, "max Radius at the absorber end"};
    Configurable<bool> enableTTCA{"enableTTCA", true, "Flag to enable or disable TTCA"};
  } dimuoncuts;

  o2::ccdb::CcdbApi ccdbApi;
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  struct : ConfigurableGroup {
    std::string prefix = "mctrackcut_group";
    Configurable<float> min_mcPt{"min_mcPt", 0.05, "min. MC pT"};
    Configurable<float> max_mcPt{"max_mcPt", 1e+10, "max. MC pT"};
    Configurable<float> min_mcEta{"min_mcEta", -4.0, "min. MC eta"};
    Configurable<float> max_mcEta{"max_mcEta", -2.5, "max. MC eta"};
  } mctrackcuts;

  HistogramRegistry fRegistry{"output", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};
  static constexpr std::string_view event_cut_types[2] = {"before/", "after/"};
  static constexpr std::string_view muon_source_types[9] = {"lf/", "Photon/", "PromptJPsi/", "NonPromptJPsi/", "PromptPsi2S/", "NonPromptPsi2S/", "c2mu/", "b2mu/", "b2c2mu/"};

  ~singleMuonQCMC() {}

  void addhistograms()
  {
    // event info
    o2::aod::pwgem::dilepton::utils::eventhistogram::addEventHistograms(&fRegistry);

    const AxisSpec axis_pt{ConfPtmuBins, "p_{T,#mu} (GeV/c)"};
    const AxisSpec axis_eta{25, -4.5, -2.0, "#eta_{#mu}"};
    const AxisSpec axis_phi{36, 0.0, 2 * M_PI, "#varphi_{#mu} (rad.)"};
    const AxisSpec axis_charge_rec{3, -1.5, +1.5, "charge"};
    const AxisSpec axis_charge_gen{3, -1.5, +1.5, "true charge"};

    // generated info
    fRegistry.add("Generated/lf/hs", "gen. single muon", kTHnSparseF, {axis_pt, axis_eta, axis_phi, axis_charge_gen}, true);
    fRegistry.addClone("Generated/lf/", "Generated/PromptJPsi/");
    fRegistry.addClone("Generated/lf/", "Generated/NonPromptJPsi/");
    fRegistry.addClone("Generated/lf/", "Generated/PromptPsi2S/");
    fRegistry.addClone("Generated/lf/", "Generated/NonPromptPsi2S/");
    fRegistry.addClone("Generated/lf/", "Generated/c2mu/");
    fRegistry.addClone("Generated/lf/", "Generated/b2mu/");
    fRegistry.addClone("Generated/lf/", "Generated/b2c2mu/");

    // track info
    fRegistry.add("Track/lf/hs", "rec. single muon", kTHnSparseF, {axis_pt, axis_eta, axis_phi, axis_charge_rec, axis_charge_gen}, true);
    fRegistry.add("Track/lf/hQoverPt", "q/pT;q/p_{T} (GeV/c)^{-1}", kTH1F, {{400, -20, 20}}, false);
    fRegistry.add("Track/lf/hTrackType", "track type", kTH1F, {{6, -0.5f, 5.5}}, false);
    fRegistry.add("Track/lf/hDCAxy", "DCA x vs. y;DCA_{x} (cm);DCA_{y} (cm)", kTH2F, {{200, -1.0f, 1.0f}, {200, -1.0f, 1.0f}}, false);
    fRegistry.add("Track/lf/hDCAxySigma", "DCA x vs. y;DCA_{x} (#sigma);DCA_{y} (#sigma)", kTH2F, {{200, -10.0f, 10.0f}, {200, -10.0f, 10.0f}}, false);
    fRegistry.add("Track/lf/hDCA2DSigma", "DCA xy;DCA_{xy} (#sigma)", kTH1F, {{100, 0.0f, 10.0f}}, false);
    fRegistry.add("Track/lf/hDCAxRes_Pt", "DCA_{x} resolution vs. pT;p_{T} (GeV/c);DCA_{x} resolution (#mum)", kTH2F, {{1000, 0, 10}, {200, 0., 200}}, false);
    fRegistry.add("Track/lf/hDCAyRes_Pt", "DCA_{y} resolution vs. pT;p_{T} (GeV/c);DCA_{y} resolution (#mum)", kTH2F, {{1000, 0, 10}, {200, 0., 200}}, false);
    fRegistry.add("Track/lf/hNclsMCH", "number of MCH clusters", kTH1F, {{21, -0.5, 20.5}}, false);
    fRegistry.add("Track/lf/hNclsMFT", "number of MFT clusters", kTH1F, {{11, -0.5, 10.5}}, false);
    fRegistry.add("Track/lf/hPDCA", "pDCA;p_{T} at PV (GeV/c);p #times DCA (GeV/c #upoint cm)", kTH2F, {{100, 0, 10}, {100, 0.0f, 1000}}, false);
    fRegistry.add("Track/lf/hChi2", "chi2;chi2", kTH1F, {{100, 0.0f, 100}}, false);
    fRegistry.add("Track/lf/hChi2MatchMCHMID", "chi2 match MCH-MID;chi2", kTH1F, {{100, 0.0f, 100}}, false);
    fRegistry.add("Track/lf/hChi2MatchMCHMFT", "chi2 match MCH-MFT;chi2", kTH1F, {{100, 0.0f, 100}}, false);
    fRegistry.add("Track/lf/hMFTClusterMap", "MFT cluster map", kTH1F, {{1024, -0.5, 1023.5}}, false);
    fRegistry.add("Track/lf/hPtGen_DeltaPtOverPtGen", "muon p_{T} resolution;p_{T}^{gen} (GeV/c);(p_{T}^{rec} - p_{T}^{gen})/p_{T}^{gen}", kTH2F, {{1000, 0, 10}, {400, -1.0f, 1.0f}}, true);
    fRegistry.add("Track/lf/hPtGen_DeltaEta", "muon #eta resolution;p_{T}^{gen} (GeV/c);#eta^{rec} - #eta^{gen}", kTH2F, {{1000, 0, 10}, {400, -1.0f, 1.0f}}, true);
    fRegistry.add("Track/lf/hPtGen_DeltaPhi", "muon #varphi resolution;p_{T}^{gen} (GeV/c);#varphi^{rec} - #varphi^{gen} (rad.)", kTH2F, {{1000, 0, 10}, {400, -1.0f, 1.0f}}, true);
    fRegistry.addClone("Track/lf/", "Track/Photon/"); // this is not for efficiency! only for contamination. We don't store generated photon conversions.
    fRegistry.addClone("Track/lf/", "Track/PromptJPsi/");
    fRegistry.addClone("Track/lf/", "Track/NonPromptJPsi/");
    fRegistry.addClone("Track/lf/", "Track/PromptPsi2S/");
    fRegistry.addClone("Track/lf/", "Track/NonPromptPsi2S/");
    fRegistry.addClone("Track/lf/", "Track/c2mu/");
    fRegistry.addClone("Track/lf/", "Track/b2mu/");
    fRegistry.addClone("Track/lf/", "Track/b2c2mu/");
  }

  void init(InitContext&)
  {
    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);

    DefineEMEventCut();
    DefineDimuonCut();
    addhistograms();
  }

  void DefineEMEventCut()
  {
    fEMEventCut = EMEventCut("fEMEventCut", "fEMEventCut");
    fEMEventCut.SetRequireSel8(eventcuts.cfgRequireSel8);
    fEMEventCut.SetRequireFT0AND(eventcuts.cfgRequireFT0AND);
    fEMEventCut.SetZvtxRange(-eventcuts.cfgZvtxMax, +eventcuts.cfgZvtxMax);
    fEMEventCut.SetRequireNoTFB(eventcuts.cfgRequireNoTFB);
    fEMEventCut.SetRequireNoITSROFB(eventcuts.cfgRequireNoITSROFB);
    fEMEventCut.SetRequireNoSameBunchPileup(eventcuts.cfgRequireNoSameBunchPileup);
    fEMEventCut.SetRequireVertexITSTPC(eventcuts.cfgRequireVertexITSTPC);
    fEMEventCut.SetRequireGoodZvtxFT0vsPV(eventcuts.cfgRequireGoodZvtxFT0vsPV);
    fEMEventCut.SetOccupancyRange(eventcuts.cfgOccupancyMin, eventcuts.cfgOccupancyMax);
  }

  void DefineDimuonCut()
  {
    fDimuonCut = DimuonCut("fDimuonCut", "fDimuonCut");

    // for pair
    fDimuonCut.SetMassRange(dimuoncuts.cfg_min_mass, dimuoncuts.cfg_max_mass);
    fDimuonCut.SetPairPtRange(dimuoncuts.cfg_min_pair_pt, dimuoncuts.cfg_max_pair_pt);
    fDimuonCut.SetPairDCAxyRange(dimuoncuts.cfg_min_pair_dcaxy, dimuoncuts.cfg_max_pair_dcaxy); // DCAxy in cm

    // for track
    fDimuonCut.SetTrackType(dimuoncuts.cfg_track_type);
    fDimuonCut.SetTrackPtRange(dimuoncuts.cfg_min_pt_track, 1e10f);
    fDimuonCut.SetTrackEtaRange(dimuoncuts.cfg_min_eta_track, dimuoncuts.cfg_max_eta_track);
    fDimuonCut.SetNClustersMFT(dimuoncuts.cfg_min_ncluster_mft, 10);
    fDimuonCut.SetNClustersMCHMID(dimuoncuts.cfg_min_ncluster_mch, 16);
    fDimuonCut.SetChi2(0.f, dimuoncuts.cfg_max_chi2);
    fDimuonCut.SetMatchingChi2MCHMFT(0.f, dimuoncuts.cfg_max_matching_chi2_mftmch);
    fDimuonCut.SetMatchingChi2MCHMID(0.f, dimuoncuts.cfg_max_matching_chi2_mchmid);
    fDimuonCut.SetDCAxy(0.f, dimuoncuts.cfg_max_dcaxy);
    fDimuonCut.SetRabs(dimuoncuts.cfg_min_rabs, dimuoncuts.cfg_max_rabs);
    fDimuonCut.SetMaxPDCARabsDep([&](float rabs) { return (rabs < 26.5 ? 594.f : 324.f); });
  }

  template <typename T>
  bool isInAcceptance(T const& t1)
  {
    if ((mctrackcuts.min_mcPt < t1.pt() && t1.pt() < mctrackcuts.max_mcPt) && (mctrackcuts.min_mcEta < t1.eta() && t1.eta() < mctrackcuts.max_mcEta)) {
      return true;
    } else {
      return false;
    }
  }

  template <int mu_source_id, typename TMCParticles, typename TTrack>
  void fillTrackInfo(TTrack const& track)
  {
    auto mctrack = track.template emmcparticle_as<TMCParticles>();
    float dca_xy = fwdDcaXYinSigma(track);
    fRegistry.fill(HIST("Track/") + HIST(muon_source_types[mu_source_id]) + HIST("hs"), track.pt(), track.eta(), track.phi(), track.sign(), -mctrack.pdgCode() / 11);
    fRegistry.fill(HIST("Track/") + HIST(muon_source_types[mu_source_id]) + HIST("hQoverPt"), track.sign() / track.pt());
    fRegistry.fill(HIST("Track/") + HIST(muon_source_types[mu_source_id]) + HIST("hTrackType"), track.trackType());
    fRegistry.fill(HIST("Track/") + HIST(muon_source_types[mu_source_id]) + HIST("hDCAxy"), track.fwdDcaX(), track.fwdDcaY());
    fRegistry.fill(HIST("Track/") + HIST(muon_source_types[mu_source_id]) + HIST("hDCAxySigma"), track.fwdDcaX() / std::sqrt(track.cXX()), track.fwdDcaY() / std::sqrt(track.cYY()));
    fRegistry.fill(HIST("Track/") + HIST(muon_source_types[mu_source_id]) + HIST("hDCA2DSigma"), dca_xy);
    fRegistry.fill(HIST("Track/") + HIST(muon_source_types[mu_source_id]) + HIST("hDCAxRes_Pt"), track.pt(), std::sqrt(track.cXX()) * 1e+4);
    fRegistry.fill(HIST("Track/") + HIST(muon_source_types[mu_source_id]) + HIST("hDCAyRes_Pt"), track.pt(), std::sqrt(track.cYY()) * 1e+4);
    fRegistry.fill(HIST("Track/") + HIST(muon_source_types[mu_source_id]) + HIST("hNclsMCH"), track.nClusters());
    fRegistry.fill(HIST("Track/") + HIST(muon_source_types[mu_source_id]) + HIST("hNclsMFT"), track.nClustersMFT());
    fRegistry.fill(HIST("Track/") + HIST(muon_source_types[mu_source_id]) + HIST("hPDCA"), track.pt(), track.pDca());
    fRegistry.fill(HIST("Track/") + HIST(muon_source_types[mu_source_id]) + HIST("hChi2"), track.chi2());
    fRegistry.fill(HIST("Track/") + HIST(muon_source_types[mu_source_id]) + HIST("hChi2MatchMCHMID"), track.chi2MatchMCHMID());
    fRegistry.fill(HIST("Track/") + HIST(muon_source_types[mu_source_id]) + HIST("hChi2MatchMCHMFT"), track.chi2MatchMCHMFT());
    fRegistry.fill(HIST("Track/") + HIST(muon_source_types[mu_source_id]) + HIST("hMFTClusterMap"), track.mftClusterMap());
    fRegistry.fill(HIST("Track/") + HIST(muon_source_types[mu_source_id]) + HIST("hPtGen_DeltaPtOverPtGen"), mctrack.pt(), (track.pt() - mctrack.pt()) / mctrack.pt());
    fRegistry.fill(HIST("Track/") + HIST(muon_source_types[mu_source_id]) + HIST("hPtGen_DeltaEta"), mctrack.pt(), track.eta() - mctrack.eta());
    fRegistry.fill(HIST("Track/") + HIST(muon_source_types[mu_source_id]) + HIST("hPtGen_DeltaPhi"), mctrack.pt(), track.phi() - mctrack.phi());
  }

  SliceCache cache;
  Preslice<MyMCTracks> perCollision_track = aod::emprimarymuon::emeventId;
  Filter trackFilter = o2::aod::fwdtrack::trackType == dimuoncuts.cfg_track_type && dimuoncuts.cfg_min_pt_track < o2::aod::fwdtrack::pt && dimuoncuts.cfg_min_eta_track < o2::aod::fwdtrack::eta && o2::aod::fwdtrack::eta < dimuoncuts.cfg_max_eta_track;
  Filter ttcaFilter = ifnode(dimuoncuts.enableTTCA.node(), o2::aod::emprimarymuon::isAssociatedToMPC == true || o2::aod::emprimarymuon::isAssociatedToMPC == false, o2::aod::emprimarymuon::isAssociatedToMPC == true);
  using FilteredMyMCTracks = soa::Filtered<MyMCTracks>;

  Filter collisionFilter_centrality = (cfgCentMin < o2::aod::cent::centFT0M && o2::aod::cent::centFT0M < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0A && o2::aod::cent::centFT0A < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0C && o2::aod::cent::centFT0C < cfgCentMax);
  using FilteredMyCollisions = soa::Filtered<MyCollisions>;

  void processQCMC(FilteredMyCollisions const& collisions, FilteredMyMCTracks const& tracks, aod::EMMCParticles const& mcparticles, aod::EMMCEvents const&)
  {
    for (auto& collision : collisions) {
      float centralities[3] = {collision.centFT0M(), collision.centFT0A(), collision.centFT0C()};
      if (centralities[cfgCentEstimator] < cfgCentMin || cfgCentMax < centralities[cfgCentEstimator]) {
        continue;
      }

      o2::aod::pwgem::dilepton::utils::eventhistogram::fillEventInfo<0, -1>(&fRegistry, collision);
      if (!fEMEventCut.IsSelected(collision)) {
        continue;
      }
      o2::aod::pwgem::dilepton::utils::eventhistogram::fillEventInfo<1, -1>(&fRegistry, collision);
      fRegistry.fill(HIST("Event/before/hCollisionCounter"), o2::aod::pwgem::dilepton::utils::eventhistogram::nbin_ev); // accepted
      fRegistry.fill(HIST("Event/after/hCollisionCounter"), o2::aod::pwgem::dilepton::utils::eventhistogram::nbin_ev);  // accepted

      auto mccollision = collision.emmcevent_as<aod::EMMCEvents>();
      if (cfgEventGeneratorType >= 0 && mccollision.getSubGeneratorId() != cfgEventGeneratorType) {
        continue;
      }

      auto tracks_per_coll = tracks.sliceBy(perCollision_track, collision.globalIndex());

      for (auto& track : tracks_per_coll) {
        auto mctrack = track.template emmcparticle_as<aod::EMMCParticles>();
        if (abs(mctrack.pdgCode()) != 13) {
          continue;
        }
        if (!fDimuonCut.IsSelectedTrack(track)) {
          continue;
        }

        if (!mctrack.has_mothers()) {
          continue;
        }
        auto mcmother = mcparticles.iteratorAt(mctrack.mothersIds()[0]);
        int pdg_mother = abs(mcmother.pdgCode());

        if (mctrack.isPhysicalPrimary() || mctrack.producedByGenerator()) {
          if (pdg_mother == 111 || pdg_mother == 221 || pdg_mother == 331 || pdg_mother == 113 || pdg_mother == 223 || pdg_mother == 333) {
            fillTrackInfo<0, aod::EMMCParticles>(track);
          } else if (pdg_mother == 443) {
            if (IsFromBeauty(mcmother, mcparticles) > 0) { // b is found in full decay chain.
              fillTrackInfo<3, aod::EMMCParticles>(track);
            } else {
              fillTrackInfo<2, aod::EMMCParticles>(track);
            }
          } else if (pdg_mother == 100443) {
            if (IsFromBeauty(mcmother, mcparticles) > 0) { // b is found in full decay chain.
              fillTrackInfo<5, aod::EMMCParticles>(track);
            } else {
              fillTrackInfo<4, aod::EMMCParticles>(track);
            }
          } else if (IsFromBeauty(mctrack, mcparticles) > 0) { // b is found in full decay chain.
            if (IsFromCharm(mctrack, mcparticles) > 0) {       // c is found in full decay chain.
              fillTrackInfo<8, aod::EMMCParticles>(track);
            } else {
              fillTrackInfo<7, aod::EMMCParticles>(track);
            }
          } else if (IsFromCharm(mctrack, mcparticles) > 0) { // c is found in full decay chain. Not from b.
            fillTrackInfo<6, aod::EMMCParticles>(track);
          }
        } else {
          fillTrackInfo<1, aod::EMMCParticles>(track);
        }

      } // end of track loop

    } // end of collision loop

  } // end of process
  PROCESS_SWITCH(singleMuonQCMC, processQCMC, "run QC MC", true);

  Partition<aod::EMMCParticles> muonsMC = nabs(o2::aod::mcparticle::pdgCode) == 13; // mu+, mu-
  PresliceUnsorted<aod::EMMCParticles> perMcCollision = aod::emmcparticle::emmceventId;
  void processGen(MyCollisions const& collisions, aod::EMMCEvents const&, aod::EMMCParticles const& mcparticles)
  {
    // loop over mc stack and fill histograms for pure MC truth signals
    // all MC tracks which belong to the MC event corresponding to the current reconstructed event

    for (auto& collision : collisions) {
      float centralities[3] = {collision.centFT0M(), collision.centFT0A(), collision.centFT0C()};
      if (centralities[cfgCentEstimator] < cfgCentMin || cfgCentMax < centralities[cfgCentEstimator]) {
        continue;
      }

      if (!fEMEventCut.IsSelected(collision)) {
        continue;
      }

      auto mccollision = collision.emmcevent_as<aod::EMMCEvents>();
      // LOGF(info, "mccollision.getGeneratorId() = %d", mccollision.getGeneratorId());
      // LOGF(info, "mccollision.getSubGeneratorId() = %d", mccollision.getSubGeneratorId());
      // LOGF(info, "mccollision.getSourceId() = %d", mccollision.getSourceId());
      if (cfgEventGeneratorType >= 0 && mccollision.getSubGeneratorId() != cfgEventGeneratorType) {
        continue;
      }

      auto muonsMC_per_coll = muonsMC->sliceByCachedUnsorted(o2::aod::emmcparticle::emmceventId, mccollision.globalIndex(), cache);

      for (auto& muon : muonsMC_per_coll) {
        if (!(muon.isPhysicalPrimary() || muon.producedByGenerator())) {
          continue;
        }
        if (!isInAcceptance(muon)) {
          continue;
        }
        if (!muon.has_mothers()) {
          continue;
        }
        auto mcmother = mcparticles.iteratorAt(muon.mothersIds()[0]);
        int pdg_mother = abs(mcmother.pdgCode());
        if (pdg_mother == 111 || pdg_mother == 221 || pdg_mother == 331 || pdg_mother == 113 || pdg_mother == 223 || pdg_mother == 333) {
          fRegistry.fill(HIST("Generated/lf/hs"), muon.pt(), muon.eta(), muon.phi(), -muon.pdgCode() / 13);
        } else if (pdg_mother == 443) {
          if (IsFromBeauty(mcmother, mcparticles) > 0) { // b is found in full decay chain.
            fRegistry.fill(HIST("Generated/NonPromptJPsi/hs"), muon.pt(), muon.eta(), muon.phi(), -muon.pdgCode() / 13);
          } else {
            fRegistry.fill(HIST("Generated/PromptJPsi/hs"), muon.pt(), muon.eta(), muon.phi(), -muon.pdgCode() / 13);
          }
        } else if (pdg_mother == 100443) {
          if (IsFromBeauty(mcmother, mcparticles) > 0) { // b is found in full decay chain.
            fRegistry.fill(HIST("Generated/NonPromptPsi2S/hs"), muon.pt(), muon.eta(), muon.phi(), -muon.pdgCode() / 13);
          } else {
            fRegistry.fill(HIST("Generated/PromptPsi2S/hs"), muon.pt(), muon.eta(), muon.phi(), -muon.pdgCode() / 13);
          }
        } else if (IsFromBeauty(muon, mcparticles) > 0) { // b is found in full decay chain.
          if (IsFromCharm(muon, mcparticles) > 0) {       // c is found in full decay chain.
            fRegistry.fill(HIST("Generated/b2c2mu/hs"), muon.pt(), muon.eta(), muon.phi(), -muon.pdgCode() / 13);
          } else {
            fRegistry.fill(HIST("Generated/b2mu/hs"), muon.pt(), muon.eta(), muon.phi(), -muon.pdgCode() / 13);
          }
        } else if (IsFromCharm(muon, mcparticles) > 0) { // c is found in full decay chain. Not from b.
          fRegistry.fill(HIST("Generated/c2mu/hs"), muon.pt(), muon.eta(), muon.phi(), -muon.pdgCode() / 13);
        }
      }

    } // end of collision loop
  }
  PROCESS_SWITCH(singleMuonQCMC, processGen, "run genrated info", true);

  void processDummy(MyCollisions const&) {}
  PROCESS_SWITCH(singleMuonQCMC, processDummy, "Dummy function", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<singleMuonQCMC>(cfgc, TaskName{"single-muon-qc-mc"})};
}
