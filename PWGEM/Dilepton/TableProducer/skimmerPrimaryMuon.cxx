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

/// \brief write relevant information for muons.
/// \author daiki.sekihata@cern.ch

#include <string>
#include <map>
#include <utility>
#include <vector>

#include "Math/Vector4D.h"
#include "Math/SMatrix.h"

#include "Framework/DataTypes.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "CommonConstants/PhysicsConstants.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/Core/TableHelper.h"
#include "Common/DataModel/PropagatedFwdTrackTables.h"

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

using namespace o2;
using namespace o2::soa;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;

struct skimmerPrimaryMuon {
  using MyCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::EMEvSels>;
  using MyCollisionsWithSWT = soa::Join<MyCollisions, aod::EMSWTriggerInfosTMP>;

  using MyPropagatedFwdTracks = soa::Join<aod::PropagatedFwdTracks, aod::PropagatedFwdTracksCov>; // muon tracks are repeated. i.e. not exclusive.
  using MyPropagatedFwdTrack = MyPropagatedFwdTracks::iterator;

  using MyFwdTracks = soa::Join<aod::FwdTracks, aod::FwdTracksCov>; // muon tracks are repeated. i.e. not exclusive.
  using MyFwdTrack = MyFwdTracks::iterator;

  using MyFwdTracksMC = soa::Join<MyFwdTracks, aod::McFwdTrackLabels>;
  using MyFwdTrackMC = MyFwdTracksMC::iterator;

  using MFTTracksMC = soa::Join<o2::aod::MFTTracks, aod::McMFTTrackLabels>;
  using MFTTrackMC = MFTTracksMC::iterator;

  Produces<aod::EMPrimaryMuons> emprimarymuons;
  Produces<aod::EMPrimaryMuonsCov> emprimarymuonscov;

  SliceCache cache;
  Preslice<aod::PropagatedFwdTracks> perCollision = o2::aod::fwdtrack::collisionId;
  Preslice<aod::MFTTracks> perCollision_mft = o2::aod::fwdtrack::collisionId;

  // Configurables
  Configurable<bool> fillQAHistogram{"fillQAHistogram", false, "flag to fill QA histograms"};
  // Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  // Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  // Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  // Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
  Configurable<float> minpt{"minpt", 0.2, "min pt for muon"};
  Configurable<float> mineta{"mineta", -4.0, "eta acceptance"};
  Configurable<float> maxeta{"maxeta", -2.5, "eta acceptance"};
  Configurable<float> mineta_mft{"mineta_mft", -3.6, "eta acceptance"};
  Configurable<float> maxeta_mft{"maxeta_mft", -2.5, "eta acceptance"};
  Configurable<float> minRabs{"minRabs", 17.6, "min. R at absorber end"};
  Configurable<float> maxRabs{"maxRabs", 89.5, "max. R at absorber end"};

  // o2::ccdb::CcdbApi ccdbApi;
  // Service<o2::ccdb::BasicCCDBManager> ccdb;
  int mRunNumber;

  // o2::globaltracking::MatchGlobalFwd mMatching;
  HistogramRegistry fRegistry{"output", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};
  static constexpr std::string_view muon_types[5] = {"MFTMCHMID/", "MFTMCHMIDOtherMatch/", "MFTMCH/", "MCHMID/", "MCH/"};

  void init(InitContext&)
  {
    // ccdb->setURL(ccdburl);
    // ccdb->setCaching(true);
    // ccdb->setLocalObjectValidityChecking();
    // ccdb->setFatalWhenNull(false);
    // ccdbApi.init(ccdburl);

    if (fillQAHistogram) {
      addHistograms();
    }
    mRunNumber = 0;
  }

  void initCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }
    mRunNumber = bc.runNumber();

    // std::map<string, string> metadata;
    // auto soreor = o2::ccdb::BasicCCDBManager::getRunDuration(ccdbApi, mRunNumber);
    // auto ts = soreor.first;
    // auto grpmag = ccdbApi.retrieveFromTFileAny<o2::parameters::GRPMagField>(grpmagPath, metadata, ts);
    // o2::base::Propagator::initFieldFromGRP(grpmag);
    // if (!o2::base::GeometryManager::isGeometryLoaded()) {
    //   ccdb->get<TGeoManager>(geoPath);
    // }
    // o2::mch::TrackExtrap::setField();
  }

  void addHistograms()
  {
    // for track
    auto hMuonType = fRegistry.add<TH1>("Track/hMuonType", "muon type", kTH1F, {{5, -0.5f, 4.5f}});
    hMuonType->GetXaxis()->SetBinLabel(1, "MFT-MCH-MID (global muon)");
    hMuonType->GetXaxis()->SetBinLabel(2, "MFT-MCH-MID (global muon other match)");
    hMuonType->GetXaxis()->SetBinLabel(3, "MFT-MCH");
    hMuonType->GetXaxis()->SetBinLabel(4, "MCH-MID");
    hMuonType->GetXaxis()->SetBinLabel(5, "MCH standalone");

    fRegistry.add("Track/MFTMCHMID/hPt", "pT;p_{T} (GeV/c)", kTH1F, {{1000, 0.0f, 10}}, false);
    fRegistry.add("Track/MFTMCHMID/hEtaPhi", "#eta vs. #varphi;#varphi (rad.);#eta", kTH2F, {{360, 0, 2 * M_PI}, {30, -5.0f, -2.0f}}, false);
    fRegistry.add("Track/MFTMCHMID/hNclusters", "Nclusters;Nclusters", kTH1F, {{21, -0.5f, 20.5}}, false);
    fRegistry.add("Track/MFTMCHMID/hNclustersMFT", "NclustersMFT;Nclusters MFT", kTH1F, {{11, -0.5f, 10.5}}, false);
    fRegistry.add("Track/MFTMCHMID/hRatAbsorberEnd", "R at absorber end;R at absorber end (cm)", kTH1F, {{100, 0.0f, 100}}, false);
    fRegistry.add("Track/MFTMCHMID/hPDCA", "pDCA;r at absorber end (cm);p #times DCA (GeV/c #upoint cm)", kTH2F, {{100, 0, 100}, {100, 0.0f, 1000}}, false);
    fRegistry.add("Track/MFTMCHMID/hChi2", "chi2;chi2", kTH1F, {{100, 0.0f, 100}}, false);
    fRegistry.add("Track/MFTMCHMID/hChi2MatchMCHMID", "chi2 match MCH-MID;chi2", kTH1F, {{100, 0.0f, 100}}, false);
    fRegistry.add("Track/MFTMCHMID/hChi2MatchMCHMFT", "chi2 match MCH-MFT;chi2", kTH1F, {{100, 0.0f, 100}}, false);
    fRegistry.add("Track/MFTMCHMID/hDCAxy2D", "DCA xy;DCA_{x} (cm);DCA_{y} (cm)", kTH2F, {{200, -10, 10}, {200, -10, +10}}, false);
    fRegistry.add("Track/MFTMCHMID/hDCAxy2DinSigma", "DCA xy;DCA_{x} (#sigma);DCA_{y} (#sigma)", kTH2F, {{200, -10, 10}, {200, -10, +10}}, false);
    fRegistry.add("Track/MFTMCHMID/hDCAxResolutionvsPt", "DCA_{x} vs. p_{T,#mu};p_{T,#mu} (GeV/c);DCA_{x} resolution (#mum);", kTH2F, {{100, 0, 10.f}, {100, 0, 100}}, false);
    fRegistry.add("Track/MFTMCHMID/hDCAyResolutionvsPt", "DCA_{y} vs. p_{T,#mu};p_{T,#mu} (GeV/c);DCA_{y} resolution (#mum);", kTH2F, {{100, 0, 10.f}, {100, 0, 100}}, false);
    fRegistry.addClone("Track/MFTMCHMID/", "Track/MCHMID/");
  }

  template <int mu_id, typename TTrack>
  void fillTrackHistogram(TTrack const& track)
  {
    fRegistry.fill(HIST("Track/hMuonType"), track.trackType());
    fRegistry.fill(HIST("Track/") + HIST(muon_types[mu_id]) + HIST("hPt"), track.pt());
    fRegistry.fill(HIST("Track/") + HIST(muon_types[mu_id]) + HIST("hEtaPhi"), track.phi(), track.eta());
    fRegistry.fill(HIST("Track/") + HIST(muon_types[mu_id]) + HIST("hNclusters"), track.nClusters());
    fRegistry.fill(HIST("Track/") + HIST(muon_types[mu_id]) + HIST("hPDCA"), track.rAtAbsorberEnd(), track.pDca());
    fRegistry.fill(HIST("Track/") + HIST(muon_types[mu_id]) + HIST("hRatAbsorberEnd"), track.rAtAbsorberEnd());
    fRegistry.fill(HIST("Track/") + HIST(muon_types[mu_id]) + HIST("hChi2"), track.chi2());
    fRegistry.fill(HIST("Track/") + HIST(muon_types[mu_id]) + HIST("hChi2MatchMCHMID"), track.chi2MatchMCHMID());
    fRegistry.fill(HIST("Track/") + HIST(muon_types[mu_id]) + HIST("hChi2MatchMCHMFT"), track.chi2MatchMCHMFT());
    fRegistry.fill(HIST("Track/") + HIST(muon_types[mu_id]) + HIST("hDCAxy2D"), track.fwdDcaX(), track.fwdDcaY());
    fRegistry.fill(HIST("Track/") + HIST(muon_types[mu_id]) + HIST("hDCAxy2DinSigma"), track.fwdDcaX() / std::sqrt(track.cXX()), track.fwdDcaY() / std::sqrt(track.cYY()));
    fRegistry.fill(HIST("Track/") + HIST(muon_types[mu_id]) + HIST("hDCAxResolutionvsPt"), track.pt(), std::sqrt(track.cXXatDCA()) * 1e+4); // convert cm to um
    fRegistry.fill(HIST("Track/") + HIST(muon_types[mu_id]) + HIST("hDCAyResolutionvsPt"), track.pt(), std::sqrt(track.cYYatDCA()) * 1e+4); // convert cm to um
  }

  template <bool isMC, bool isTriggerAnalysis, typename TCollisions, typename TBCs, typename TSAMuons, typename TGlobalMuons, typename TFwdTracks, typename TMFTTracks>
  void run(TCollisions const& collisions, TBCs const&, TSAMuons const& saMuons, TGlobalMuons const& glMuons, TFwdTracks const&, TMFTTracks const&)
  {
    for (const auto& collision : collisions) {
      auto bc = collision.template foundBC_as<aod::BCsWithTimestamps>();
      initCCDB(bc);
      if (!collision.isSelected()) {
        continue;
      }

      if constexpr (isTriggerAnalysis) {
        if (collision.swtaliastmp_raw() == 0) {
          continue;
        }
      }

      if constexpr (isMC) {
        if (!collision.has_mcCollision()) {
          continue;
        }
      }

      auto sa_muons_per_coll = saMuons.sliceByCached(o2::aod::fwdtrack::collisionId, collision.globalIndex(), cache);
      auto global_muons_per_coll = glMuons.sliceByCached(o2::aod::fwdtrack::collisionId, collision.globalIndex(), cache);

      for (const auto& muon : sa_muons_per_coll) {
        if (fillQAHistogram) {
          fillTrackHistogram<3>(muon);
        }

        const auto& fwdtrack = muon.template fwdtrack_as<TFwdTracks>();

        if constexpr (isMC) {
          if (!fwdtrack.has_mcParticle()) {
            continue;
          }
        }
        emprimarymuons(collision.globalIndex(), fwdtrack.globalIndex(), -1, -1, muon.trackType(),
                       muon.pt(), muon.eta(), muon.phi(), muon.sign(), muon.fwdDcaX(), muon.fwdDcaY(), muon.cXXatDCA(), muon.cYYatDCA(), muon.cXYatDCA(), muon.etaMatchedMCHMID(), muon.phiMatchedMCHMID(),
                       // muon.x(), muon.y(), muon.z(), muon.tgl(),
                       muon.nClusters(), muon.pDca(), muon.rAtAbsorberEnd(), muon.chi2(), muon.chi2MatchMCHMID(), muon.chi2MatchMCHMFT(),
                       muon.mchBitMap(), muon.midBitMap(), muon.midBoards(), 0, 999999.f, muon.isAssociatedToMPC(), muon.isAmbiguous());

        emprimarymuonscov(
          muon.cXX(),
          muon.cXY(),
          muon.cYY(),
          muon.cPhiX(),
          muon.cPhiY(),
          muon.cPhiPhi(),
          muon.cTglX(),
          muon.cTglY(),
          muon.cTglPhi(),
          muon.cTglTgl(),
          muon.c1PtX(),
          muon.c1PtY(),
          muon.c1PtPhi(),
          muon.c1PtTgl(),
          muon.c1Pt21Pt2());

      } // end of standalone muon loop
      for (const auto& muon : global_muons_per_coll) {
        if (fillQAHistogram) {
          fillTrackHistogram<0>(muon);
        }

        const auto& fwdtrack = muon.template fwdtrack_as<TFwdTracks>();
        const auto& mfttrack = muon.template matchMFTTrack_as<TMFTTracks>();
        const auto& mchtrack = muon.template matchMCHTrack_as<TFwdTracks>();

        if constexpr (isMC) {
          if (!fwdtrack.has_mcParticle()) {
            continue;
          }
        }

        emprimarymuons(collision.globalIndex(), fwdtrack.globalIndex(), mfttrack.globalIndex(), mchtrack.globalIndex(), muon.trackType(),
                       muon.pt(), muon.eta(), muon.phi(), muon.sign(), muon.fwdDcaX(), muon.fwdDcaY(), muon.cXXatDCA(), muon.cYYatDCA(), muon.cXYatDCA(), muon.etaMatchedMCHMID(), muon.phiMatchedMCHMID(),
                       // muon.x(), muon.y(), muon.z(), muon.tgl(),
                       muon.nClusters(), muon.pDca(), muon.rAtAbsorberEnd(), muon.chi2(), muon.chi2MatchMCHMID(), muon.chi2MatchMCHMFT(),
                       muon.mchBitMap(), muon.midBitMap(), muon.midBoards(), mfttrack.mftClusterSizesAndTrackFlags(), mfttrack.chi2(), muon.isAssociatedToMPC(), muon.isAmbiguous());

        emprimarymuonscov(
          muon.cXX(),
          muon.cXY(),
          muon.cYY(),
          muon.cPhiX(),
          muon.cPhiY(),
          muon.cPhiPhi(),
          muon.cTglX(),
          muon.cTglY(),
          muon.cTglPhi(),
          muon.cTglTgl(),
          muon.c1PtX(),
          muon.c1PtY(),
          muon.c1PtPhi(),
          muon.c1PtTgl(),
          muon.c1Pt21Pt2());

      } // end of global muon loop
    } // end of collision loop
  }

  std::map<std::pair<int, int>, int> map_new_sa_muon_index; // new standalone muon index

  Partition<MyPropagatedFwdTracks> global_muons = o2::aod::fwdtrack::trackType == uint8_t(o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack); // MFT-MCH-MID
  Partition<MyPropagatedFwdTracks> sa_muons = o2::aod::fwdtrack::trackType == uint8_t(o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack); // MCH-MID

  void processRec(MyCollisions const& collisions, aod::BCsWithTimestamps const& bcs, MyPropagatedFwdTracks const&, MyFwdTracks const& fwdtracks, aod::MFTTracks const& mfttracks)
  {
    run<false, false>(collisions, bcs, sa_muons, global_muons, fwdtracks, mfttracks);
  }
  PROCESS_SWITCH(skimmerPrimaryMuon, processRec, "process reconstructed info", true);

  void processRec_SWT(MyCollisionsWithSWT const& collisions, aod::BCsWithTimestamps const& bcs, MyPropagatedFwdTracks const&, MyFwdTracks const& fwdtracks, aod::MFTTracks const& mfttracks)
  {
    run<false, true>(collisions, bcs, sa_muons, global_muons, fwdtracks, mfttracks);
  }
  PROCESS_SWITCH(skimmerPrimaryMuon, processRec_SWT, "process reconstructed info only with standalone", false);

  void processMC(soa::Join<MyCollisions, aod::McCollisionLabels> const& collisions, aod::BCsWithTimestamps const& bcs, MyPropagatedFwdTracks const&, MyFwdTracksMC const& fwdtracks, MFTTracksMC const& mfttracks)
  {
    run<true, false>(collisions, bcs, sa_muons, global_muons, fwdtracks, mfttracks);
  }
  PROCESS_SWITCH(skimmerPrimaryMuon, processMC, "process reconstructed and MC info", false);
};
struct associateAmbiguousMuon {
  Produces<aod::EMAmbiguousMuonSelfIds> em_amb_muon_ids;

  SliceCache cache;
  PresliceUnsorted<aod::EMPrimaryMuons> perTrack = o2::aod::emprimarymuon::fwdtrackId;
  std::vector<int> ambmuon_self_Ids;

  void process(aod::EMPrimaryMuons const& muons)
  {
    for (const auto& muon : muons) {
      auto muons_with_same_trackId = muons.sliceBy(perTrack, muon.fwdtrackId());
      ambmuon_self_Ids.reserve(muons_with_same_trackId.size());
      for (const auto& amb_muon : muons_with_same_trackId) {
        if (amb_muon.globalIndex() == muon.globalIndex()) { // don't store myself.
          continue;
        }
        ambmuon_self_Ids.emplace_back(amb_muon.globalIndex());
      }
      em_amb_muon_ids(ambmuon_self_Ids);
      ambmuon_self_Ids.clear();
      ambmuon_self_Ids.shrink_to_fit();
    }
  }
};
struct associateSameMFT {
  Produces<aod::EMGlobalMuonSelfIds> em_same_mft_ids;

  SliceCache cache;
  PresliceUnsorted<aod::EMPrimaryMuons> perMFTTrack = o2::aod::emprimarymuon::mfttrackId;
  std::vector<int> self_Ids;

  void process(aod::EMPrimaryMuons const& muons)
  {
    for (const auto& muon : muons) {
      if (muon.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack) {
        auto muons_with_same_mfttrackId = muons.sliceBy(perMFTTrack, muon.mfttrackId());
        self_Ids.reserve(muons_with_same_mfttrackId.size());
        for (const auto& global_muon : muons_with_same_mfttrackId) {
          if (global_muon.globalIndex() == muon.globalIndex()) { // don't store myself.
            continue;
          }
          if (global_muon.collisionId() == muon.collisionId()) {
            self_Ids.emplace_back(global_muon.globalIndex());
          }
        }
        em_same_mft_ids(self_Ids);
        self_Ids.clear();
        self_Ids.shrink_to_fit();
      } else {                               // for standalone muons
        em_same_mft_ids(std::vector<int>{}); // empty
      }
    } // end of muon loop
  }
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<skimmerPrimaryMuon>(cfgc, TaskName{"skimmer-primary-muon"}),
    adaptAnalysisTask<associateAmbiguousMuon>(cfgc, TaskName{"associate-ambiguous-muon"}),
    adaptAnalysisTask<associateSameMFT>(cfgc, TaskName{"associate-same-mft"})};
}
