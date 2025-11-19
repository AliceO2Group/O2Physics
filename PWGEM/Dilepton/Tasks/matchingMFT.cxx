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

/// \file matchingMFT.cxx
/// \brief a task to study matching MFT-[MCH-MID] in MC
/// \author daiki.sekihata@cern.ch

#include "TableHelper.h"

#include "Common/CCDB/RCTSelectionFlags.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/fwdtrackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"

#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DetectorsBase/Propagator.h"
#include "Field/MagneticField.h"
#include "Framework/AnalysisTask.h"
#include "Framework/DataTypes.h"
#include "Framework/runDataProcessing.h"
#include "GlobalTracking/MatchGlobalFwd.h"
#include "MCHTracking/TrackExtrap.h"
#include "MCHTracking/TrackParam.h"
#include "ReconstructionDataFormats/TrackFwd.h"

#include "TGeoGlobalMagField.h"

#include <format>
#include <map>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

using namespace o2;
using namespace o2::soa;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;
using namespace o2::aod::fwdtrackutils;

struct matchingMFT {
  using MyCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::McCollisionLabels>;
  using MyFwdTracks = soa::Join<aod::FwdTracks, aod::FwdTracksCov, aod::McFwdTrackLabels>;
  using MyMFTTracks = soa::Join<o2::aod::MFTTracks, aod::McMFTTrackLabels>;

  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
  Configurable<float> minPt{"minPt", 0.01, "min pt for muon"};
  Configurable<float> maxPt{"maxPt", 1e+10, "max pt for muon"};
  Configurable<float> minEtaSA{"minEtaSA", -4.0, "min. eta acceptance for MCH-MID"};
  Configurable<float> maxEtaSA{"maxEtaSA", -2.5, "max. eta acceptance for MCH-MID"};
  Configurable<float> minEtaGL{"minEtaGL", -3.6, "min. eta acceptance for MFT-MCH-MID"};
  Configurable<float> maxEtaGL{"maxEtaGL", -2.5, "max. eta acceptance for MFT-MCH-MID"};
  Configurable<float> minRabs{"minRabs", 17.6, "min. R at absorber end"};
  Configurable<float> midRabs{"midRabs", 26.5, "middle R at absorber end for pDCA cut"};
  Configurable<float> maxRabs{"maxRabs", 89.5, "max. R at absorber end"};
  Configurable<float> maxDCAxy{"maxDCAxy", 1e+10, "max. DCAxy for global muons"};
  Configurable<float> maxPDCAforLargeR{"maxPDCAforLargeR", 324.f, "max. pDCA for large R at absorber end"};
  Configurable<float> maxPDCAforSmallR{"maxPDCAforSmallR", 594.f, "max. pDCA for small R at absorber end"};
  Configurable<float> maxMatchingChi2MCHMFT{"maxMatchingChi2MCHMFT", 1e+10, "max. chi2 for MCH-MFT matching"};
  Configurable<float> maxChi2SA{"maxChi2SA", 1e+6f, "max. chi2 for standalone muon"};
  Configurable<float> maxChi2GL{"maxChi2GL", 1e+6f, "max. chi2 for global muon"};
  Configurable<float> maxChi2MFT{"maxChi2MFT", 1e+6f, "max. chi2/ndf for MFT track in global muon"};
  Configurable<int> minNclustersMFT{"minNclustersMFT", 5, "min nclusters MFT"};
  Configurable<bool> refitGlobalMuon{"refitGlobalMuon", true, "flag to refit global muon"};
  Configurable<bool> propagateToDCAhelix{"propagateToDCAhelix", false, "flag to use propagateToDCAhelix"};

  Configurable<bool> requireTrueAssociation{"requireTrueAssociation", false, "flag to require true mc collision association"};
  Configurable<float> maxRelDPt{"maxRelDPt", 1e+10f, "max. relative dpt between MFT-MCH-MID and MCH-MID"};
  Configurable<float> maxDEta{"maxDEta", 1e+10f, "max. deta between MFT-MCH-MID and MCH-MID"};
  Configurable<float> maxDPhi{"maxDPhi", 1e+10f, "max. dphi between MFT-MCH-MID and MCH-MID"};
  Configurable<float> maxDEtaMP{"maxDEtaMP", 1e+10f, "max. deta between MFT and MCH-MID at matching plane Z"};
  Configurable<float> maxDPhiMP{"maxDPhiMP", 1e+10f, "max. dphi between MFT and MCH-MID at matching plane Z"};
  Configurable<bool> requireMFTHitMap{"requireMFTHitMap", false, "flag to require MFT hit map"};
  Configurable<std::vector<int>> requiredMFTDisks{"requiredMFTDisks", std::vector<int>{0}, "hit map on MFT disks [0,1,2,3,4]. logical-OR of each double-sided disk"};
  Configurable<float> matchingZ{"matchingZ", -77.5, "z position where matching is performed"};
  Configurable<int> cfgBestMatchFinder{"cfgBestMatchFinder", 0, "matching chi2:0, dr:1"};

  Configurable<int> cfgCentEstimator{"cfgCentEstimator", 2, "FT0M:0, FT0A:1, FT0C:2"};
  Configurable<float> cfgCentMin{"cfgCentMin", -1.f, "min. centrality"};
  Configurable<float> cfgCentMax{"cfgCentMax", 999.f, "max. centrality"};
  Configurable<float> cfgZvtxMin{"cfgZvtxMin", -10.f, "min. Zvtx"};
  Configurable<float> cfgZvtxMax{"cfgZvtxMax", 10.f, "max. Zvtx"};

  // for RCT
  Configurable<bool> cfgRequireGoodRCT{"cfgRequireGoodRCT", false, "require good detector flag in run condtion table"};
  Configurable<std::string> cfgRCTLabel{"cfgRCTLabel", "CBT_muon_glo", "select 1 [CBT_muon, CBT_muon_glo] see O2Physics/Common/CCDB/RCTSelectionFlags.h"};
  Configurable<bool> cfgCheckZDC{"cfgCheckZDC", false, "set ZDC flag for PbPb"};
  Configurable<bool> cfgTreatLimitedAcceptanceAsBad{"cfgTreatLimitedAcceptanceAsBad", false, "reject all events where the detectors relevant for the specified Runlist are flagged as LimitedAcceptance"};

  o2::aod::rctsel::RCTFlagsChecker rctChecker;

  HistogramRegistry fRegistry{"fRegistry"};
  static constexpr std::string_view muon_types[5] = {"MFTMCHMID/", "MFTMCHMIDOtherMatch/", "MFTMCH/", "MCHMID/", "MCH/"};

  void init(o2::framework::InitContext&)
  {
    if (doprocessWithoutFTTCA && doprocessWithFTTCA) {
      LOGF(fatal, "Cannot enable doprocessWithoutFTTCA and doprocessWithFTTCA at the same time. Please choose one.");
    }

    if (refitGlobalMuon && propagateToDCAhelix) {
      LOGF(fatal, "Cannot enable doprocessWithoutFTTCA and doprocessWithFTTCA at the same time. Please choose one.");
    }

    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);
    ccdbApi.init(ccdburl);
    rctChecker.init(cfgRCTLabel.value, cfgCheckZDC.value, cfgTreatLimitedAcceptanceAsBad.value);

    addHistograms();
  }

  o2::ccdb::CcdbApi ccdbApi;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  int mRunNumber = -1;
  float mBz = 0;

  template <typename TBC>
  void initCCDB(TBC const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }
    mRunNumber = bc.runNumber();
    LOGF(info, "mRunNumber = %d", mRunNumber);
    std::map<std::string, std::string> metadata;
    auto soreor = o2::ccdb::BasicCCDBManager::getRunDuration(ccdbApi, mRunNumber);
    auto ts = soreor.first;
    auto grpmag = ccdbApi.retrieveFromTFileAny<o2::parameters::GRPMagField>(grpmagPath, metadata, ts);
    o2::base::Propagator::initFieldFromGRP(grpmag);
    if (!o2::base::GeometryManager::isGeometryLoaded()) {
      ccdb->get<TGeoManager>(geoPath);
    }
    o2::mch::TrackExtrap::setField();
    const double centerMFT[3] = {0, 0, -61.4};
    o2::field::MagneticField* field = static_cast<o2::field::MagneticField*>(TGeoGlobalMagField::Instance()->GetField());
    mBz = field->getBz(centerMFT); // Get field at centre of MFT
    LOGF(info, "Bz at center of MFT = %f kZG", mBz);
  }

  void addHistograms()
  {
    auto hCollisionCounter = fRegistry.add<TH1>("Event/hCollisionCounter", "collision counter", kTH1F, {{5, -0.5f, 4.5f}}, false);
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

    auto hMuonType = fRegistry.add<TH1>("hMuonType", "muon type", kTH1F, {{5, -0.5f, 4.5f}}, false);
    hMuonType->GetXaxis()->SetBinLabel(1, "MFT-MCH-MID (global muon)");
    hMuonType->GetXaxis()->SetBinLabel(2, "MFT-MCH-MID (global muon other match)");
    hMuonType->GetXaxis()->SetBinLabel(3, "MFT-MCH");
    hMuonType->GetXaxis()->SetBinLabel(4, "MCH-MID");
    hMuonType->GetXaxis()->SetBinLabel(5, "MCH standalone");

    const AxisSpec axis_pt{{0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10}, "p_{T}^{gl} (GeV/c)"};

    fRegistry.add("MFTMCHMID/primary/correct/hPt", "pT;p_{T} (GeV/c)", kTH1F, {{100, 0.0f, 10}}, false);
    fRegistry.add("MFTMCHMID/primary/correct/hEtaPhi", "#eta vs. #varphi;#varphi (rad.);#eta", kTH2F, {{180, 0, 2 * M_PI}, {80, -5.f, -1.f}}, false);
    fRegistry.add("MFTMCHMID/primary/correct/hEtaPhi_MatchedMCHMID", "#eta vs. #varphi;#varphi (rad.);#eta", kTH2F, {{180, 0, 2 * M_PI}, {80, -5.f, -1.f}}, false);
    fRegistry.add("MFTMCHMID/primary/correct/hsDelta", "diff. between GL and associated SA;p_{T}^{gl} (GeV/c);(p_{T}^{sa} - p_{T}^{gl})/p_{T}^{gl};#eta^{sa} - #eta^{gl};#varphi^{sa} - #varphi^{gl} (rad.);", kTHnSparseF, {axis_pt, {200, -0.5, +0.5}, {200, -1, +1}, {90, -M_PI / 4, M_PI / 4}}, false);
    fRegistry.add("MFTMCHMID/primary/correct/hsDeltaAtMP", "diff. between MFT and MCH-MID;p_{T}^{gl} (GeV/c);#varphi^{MCH-MID} - #varphi^{MFT} (rad.);#eta^{MCH-MID} - #eta^{MFT};", kTHnSparseF, {axis_pt, {90, -M_PI / 4, M_PI / 4}, {200, -1, +1}}, false);
    fRegistry.add("MFTMCHMID/primary/correct/hDiffCollId", "difference in collision index;collisionId_{TTCA} - collisionId_{MP}", kTH1F, {{41, -20.5, +20.5}}, false);
    fRegistry.add("MFTMCHMID/primary/correct/hSign", "sign;sign", kTH1F, {{3, -1.5, +1.5}}, false);
    fRegistry.add("MFTMCHMID/primary/correct/hNclusters", "Nclusters;Nclusters", kTH1F, {{21, -0.5f, 20.5}}, false);
    fRegistry.add("MFTMCHMID/primary/correct/hNclustersMFT", "NclustersMFT;Nclusters MFT", kTH1F, {{11, -0.5f, 10.5}}, false);
    fRegistry.add("MFTMCHMID/primary/correct/hMFTClusterMap", "MFT cluster map", kTH1F, {{1024, -0.5, 1023.5}}, false);
    fRegistry.add("MFTMCHMID/primary/correct/hRatAbsorberEnd", "R at absorber end;R at absorber end (cm)", kTH1F, {{100, 0.0f, 100}}, false);
    fRegistry.add("MFTMCHMID/primary/correct/hPDCA_Rabs", "pDCA vs. Rabs;R at absorber end (cm);p #times DCA (GeV/c #upoint cm)", kTH2F, {{100, 0, 100}, {100, 0.0f, 1000}}, false);
    fRegistry.add("MFTMCHMID/primary/correct/hChi2", "chi2;chi2/ndf", kTH1F, {{100, 0.0f, 10}}, false);
    fRegistry.add("MFTMCHMID/primary/correct/hChi2MFT", "chi2 MFT/ndf;chi2 MFT/ndf", kTH1F, {{100, 0.0f, 10}}, false);
    fRegistry.add("MFTMCHMID/primary/correct/hChi2MatchMCHMID", "chi2 match MCH-MID;chi2", kTH1F, {{100, 0.0f, 100}}, false);
    fRegistry.add("MFTMCHMID/primary/correct/hChi2MatchMCHMFT", "chi2 match MCH-MFT;chi2", kTH1F, {{100, 0.0f, 100}}, false);
    fRegistry.add("MFTMCHMID/primary/correct/hDCAxy2D", "DCA x vs. y;DCA_{x} (cm);DCA_{y} (cm)", kTH2F, {{200, -0.5, 0.5}, {200, -0.5, +0.5}}, false);
    fRegistry.add("MFTMCHMID/primary/correct/hDCAxy2DinSigma", "DCA x vs. y in sigma;DCA_{x} (#sigma);DCA_{y} (#sigma)", kTH2F, {{200, -10, 10}, {200, -10, +10}}, false);
    fRegistry.add("MFTMCHMID/primary/correct/hDCAxy", "DCAxy;DCA_{xy} (cm);", kTH1F, {{100, 0, 1}}, false);
    fRegistry.add("MFTMCHMID/primary/correct/hDCAxyinSigma", "DCAxy in sigma;DCA_{xy} (#sigma);", kTH1F, {{100, 0, 10}}, false);
    fRegistry.add("MFTMCHMID/primary/correct/hDCAxResolutionvsPt", "DCA_{x} resolution vs. p_{T};p_{T} (GeV/c);DCA_{x} resolution (#mum);", kTH2F, {{100, 0, 10.f}, {500, 0, 500}}, false);
    fRegistry.add("MFTMCHMID/primary/correct/hDCAyResolutionvsPt", "DCA_{y} resolution vs. p_{T};p_{T} (GeV/c);DCA_{y} resolution (#mum);", kTH2F, {{100, 0, 10.f}, {500, 0, 500}}, false);
    fRegistry.add("MFTMCHMID/primary/correct/hDCAxyResolutionvsPt", "DCA_{xy} resolution vs. p_{T};p_{T} (GeV/c);DCA_{xy} resolution (#mum);", kTH2F, {{100, 0, 10.f}, {500, 0, 500}}, false);
    fRegistry.add("MFTMCHMID/primary/correct/hDCAz", "DCAz;DCA_{z} (cm);", kTH1F, {{200, -0.1, 0.1}}, false);
    fRegistry.add("MFTMCHMID/primary/correct/hDCAxyz", "DCA xy vs. z;DCA_{xy} (cm);DCA_{z} (cm);", kTH2F, {{100, 0, 1}, {200, -0.1, 0.1}}, false);
    fRegistry.add("MFTMCHMID/primary/correct/hMCHBitMap", "MCH bit map;MCH bit map", kTH1F, {{1024, -0.5, 1023.5}}, false);
    fRegistry.add("MFTMCHMID/primary/correct/hMIDBitMap", "MID bit map;MID bit map", kTH1F, {{256, -0.5, 255.5}}, false);

    fRegistry.add("MFTMCHMID/primary/correct/hProdVtxZ", "prod. vtx Z of muon;V_{z} (cm)", kTH1F, {{200, -100, 100}}, false);
    fRegistry.add("MFTMCHMID/primary/correct/hRelDeltaPt", "pT resolution;p_{T}^{gen} (GeV/c);(p_{T}^{rec} - p_{T}^{gen})/p_{T}^{gen}", kTH2F, {{100, 0, 10}, {200, -1, +1}}, false);
    fRegistry.add("MFTMCHMID/primary/correct/hDeltaEta_Pos", "#eta resolution;p_{T}^{gen} (GeV/c);#eta^{rec} - #eta^{gen}", kTH2F, {{100, 0, 10}, {400, -0.2, +0.2}}, false);
    fRegistry.add("MFTMCHMID/primary/correct/hDeltaEta_Neg", "#eta resolution;p_{T}^{gen} (GeV/c);#eta^{rec} - #eta^{gen}", kTH2F, {{100, 0, 10}, {400, -0.2, +0.2}}, false);
    fRegistry.add("MFTMCHMID/primary/correct/hDeltaPhi_Pos", "#varphi resolution;p_{T}^{gen} (GeV/c);#varphi^{rec} - #varphi^{gen} (rad.)", kTH2F, {{100, 0, 10}, {400, -0.2, +0.2}}, false);
    fRegistry.add("MFTMCHMID/primary/correct/hDeltaPhi_Neg", "#varphi resolution;p_{T}^{gen} (GeV/c);#varphi^{rec} - #varphi^{gen} (rad.)", kTH2F, {{100, 0, 10}, {400, -0.2, +0.2}}, false);
    fRegistry.addClone("MFTMCHMID/primary/correct/", "MFTMCHMID/primary/wrong/");
    fRegistry.addClone("MFTMCHMID/primary/", "MFTMCHMID/secondary/");
    // fRegistry.addClone("MFTMCHMID/", "MCHMID/");
  }

  bool isSelected(const float pt, const float eta, const float rAtAbsorberEnd, const float pDCA, const float chi2_per_ndf, const uint8_t trackType, const float dcaXY)
  {
    if (pt < minPt || maxPt < pt) {
      return false;
    }
    if (rAtAbsorberEnd < minRabs || maxRabs < rAtAbsorberEnd) {
      return false;
    }
    if (rAtAbsorberEnd < midRabs ? pDCA > maxPDCAforSmallR : pDCA > maxPDCAforLargeR) {
      return false;
    }

    if (trackType == static_cast<uint8_t>(o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack)) {
      if (eta < minEtaGL || maxEtaGL < eta) {
        return false;
      }
      if (maxDCAxy < dcaXY) {
        return false;
      }
      if (chi2_per_ndf < 0.f || maxChi2GL < chi2_per_ndf) {
        return false;
      }
    } else if (trackType == static_cast<uint8_t>(o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack)) {
      if (eta < minEtaSA || maxEtaSA < eta) {
        return false;
      }
      if (chi2_per_ndf < 0.f || maxChi2SA < chi2_per_ndf) {
        return false;
      }
    } else {
      return false;
    }

    return true;
  }

  template <typename T>
  uint16_t mftClusterMap(T const& track)
  {
    uint64_t mftClusterSizesAndTrackFlags = track.mftClusterSizesAndTrackFlags();
    uint16_t clmap = 0;
    for (unsigned int layer = 0; layer < 10; layer++) {
      if ((mftClusterSizesAndTrackFlags >> (layer * 6)) & 0x3f) {
        clmap |= (1 << layer);
      }
    }
    return clmap;
  }

  template <int begin = 0, int end = 9, typename T>
  bool hasMFT(T const& track)
  {
    // logical-OR
    uint64_t mftClusterSizesAndTrackFlags = track.mftClusterSizesAndTrackFlags();
    uint16_t clmap = 0;
    for (unsigned int layer = begin; layer <= end; layer++) {
      if ((mftClusterSizesAndTrackFlags >> (layer * 6)) & 0x3f) {
        clmap |= (1 << layer);
      }
    }
    return (clmap > 0);
  }

  // template <typename T>
  // float meanClusterSizeMFT(T const& track)
  // {
  //   uint64_t mftClusterSizesAndTrackFlags = track.mftClusterSizesAndTrackFlags();
  //   uint16_t clsSize = 0;
  //   uint16_t n = 0;
  //   for (unsigned int layer = 0; layer < 10; layer++) {
  //     uint16_t size_per_layer = (mftClusterSizesAndTrackFlags >> (layer * 6)) & 0x3f;
  //     clsSize += size_per_layer;
  //     if (size_per_layer > 0) {
  //       n++;
  //     }
  //     // LOGF(info, "track.globalIndex() = %d, layer = %d, size_per_layer = %d", track.globalIndex(), layer, size_per_layer);
  //   }

  //   if (n > 0) {
  //     return static_cast<float>(clsSize) / static_cast<float>(n) * std::fabs(std::sin(std::atan(track.tgl())));
  //   } else {
  //     return 0.f;
  //   }
  // }

  template <typename TFwdTracks, typename TMFTTracks, typename TCollision, typename TFwdTrack, typename TMFTrackCov>
  void getDeltaEtaDeltaPhiAtMatchingPlane(TCollision const& collision, TFwdTrack const& fwdtrack, TMFTrackCov const& mftCovs, float& deta, float& dphi)
  {
    if (fwdtrack.trackType() != o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack) {
      deta = 999.f;
      dphi = 999.f;
      return; // do nothing
    }

    auto mchtrack = fwdtrack.template matchMCHTrack_as<TFwdTracks>(); // MCH-MID
    auto mfttrack = fwdtrack.template matchMFTTrack_as<TMFTTracks>();

    if (mchtrack.trackType() != o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack) {
      deta = 999.f;
      dphi = 999.f;
      return; // do nothing
    }

    o2::dataformats::GlobalFwdTrack propmuonAtPV = propagateMuon(mchtrack, mchtrack, collision, propagationPoint::kToMatchingPlane, matchingZ, mBz);
    auto mfttrackcov = mftCovs.rawIteratorAt(map_mfttrackcovs[mfttrack.globalIndex()]);

    auto muonAtMP = propagateMuon(mchtrack, mchtrack, collision, propagationPoint::kToMatchingPlane, matchingZ, mBz); // propagated to matching plane
    o2::track::TrackParCovFwd mftsaAtMP = getTrackParCovFwd(mfttrack, mfttrackcov);                                   // values at innermost update
    mftsaAtMP.propagateToZhelix(matchingZ, mBz);                                                                      // propagated to matching plane
    deta = muonAtMP.getEta() - mftsaAtMP.getEta();
    dphi = muonAtMP.getPhi() - mftsaAtMP.getPhi();
    o2::math_utils::bringToPMPi(dphi);
    // reldpt = (muonAtMP.getPt() - mftsaAtMP.getPt()) / muonAtMP.getPt();

    // o2::track::TrackParCovFwd mftsa = getTrackParCovFwd(mfttrack, mfttrackcov); // values at innermost update
    // o2::dataformats::GlobalFwdTrack globalMuonRefit = o2::aod::fwdtrackutils::refitGlobalMuonCov(propmuonAtPV, mftsa); // this is track at IU.
    // auto globalMuonAtDCA = o2::aod::fwdtrackutils::propagateTrackParCovFwd(globalMuonRefit, fwdtrack.trackType(), collision, propagationPoint::kToMatchingPlane, matchingZ);
    // deta = propmuonAtPV.getEta() - globalMuonAtDCA.getEta();
    // dphi = propmuonAtPV.getPhi() - globalMuonAtDCA.getPhi();
    // o2::math_utils::bringToPMPi(dphi);
    // reldpt = (globalMuonAtDCA.getPt() - propmuonAtPV.getPt()) / propmuonAtPV.getPt();
  }

  template <bool withMFTCov = false, typename TCollision, typename TFwdTrack, typename TFwdTracks, typename TMFTTracks, typename TMFTTracksCov>
  void fillHistograms(TCollision const& collision, TFwdTrack fwdtrack, TFwdTracks const&, TMFTTracks const&, TMFTTracksCov const& mftCovs)
  {
    if (fwdtrack.trackType() != o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack) { // only for protection
      return;
    }

    auto mchtrack = fwdtrack.template matchMCHTrack_as<TFwdTracks>(); // MCH-MID
    auto mfttrack = fwdtrack.template matchMFTTrack_as<TMFTTracks>();

    if (!fwdtrack.has_mcParticle() || !mchtrack.has_mcParticle() || !mfttrack.has_mcParticle()) {
      return;
    }

    auto mcParticle_MFTMCHMID = fwdtrack.template mcParticle_as<aod::McParticles>(); // this is identical to mcParticle_MCHMID
    auto mcParticle_MCHMID = mchtrack.template mcParticle_as<aod::McParticles>();    // this is identical to mcParticle_MFTMCHMID
    auto mcParticle_MFT = mfttrack.template mcParticle_as<aod::McParticles>();
    // LOGF(info, "mcParticle_MFTMCHMID.pdgCode() = %d, mcParticle_MCHMID.pdgCode() = %d, mcParticle_MFT.pdgCode() = %d", mcParticle_MFTMCHMID.pdgCode(), mcParticle_MCHMID.pdgCode(), mcParticle_MFT.pdgCode());
    // LOGF(info, "mcParticle_MFTMCHMID.globalIndex() = %d, mcParticle_MCHMID.globalIndex() = %d, mcParticle_MFT.globalIndex() = %d", mcParticle_MFTMCHMID.globalIndex(), mcParticle_MCHMID.globalIndex(), mcParticle_MFT.globalIndex());

    int nClustersMFT = mfttrack.nClusters();
    float chi2mft = mfttrack.chi2() / (2.f * nClustersMFT - 5.f);
    if (chi2mft < 0.f || maxChi2MFT < chi2mft) {
      return;
    }

    if (fwdtrack.chi2MatchMCHMFT() > maxMatchingChi2MCHMFT) {
      return;
    }

    if (fwdtrack.chi2() < 0.f || maxChi2GL < fwdtrack.chi2() / (2.f * (mchtrack.nClusters() + nClustersMFT) - 5.f)) {
      return;
    }

    if (fwdtrack.rAtAbsorberEnd() < minRabs || maxRabs < fwdtrack.rAtAbsorberEnd()) {
      return;
    }

    if (nClustersMFT < minNclustersMFT) {
      return;
    }

    if (std::abs(mcParticle_MCHMID.pdgCode()) != 13) { // select true muon
      return;
    }

    if (requireTrueAssociation && (mcParticle_MCHMID.mcCollisionId() != collision.mcCollisionId())) {
      return;
    }

    bool isPrimary = mcParticle_MCHMID.isPhysicalPrimary() || mcParticle_MCHMID.producedByGenerator();
    bool isMatched = (mcParticle_MFT.globalIndex() == mcParticle_MCHMID.globalIndex()) && (mcParticle_MFT.mcCollisionId() == mcParticle_MCHMID.mcCollisionId());

    o2::dataformats::GlobalFwdTrack propmuonAtPV = propagateMuon(fwdtrack, fwdtrack, collision, propagationPoint::kToVertex, matchingZ, mBz);
    o2::dataformats::GlobalFwdTrack propmuonAtDCA = propagateMuon(fwdtrack, fwdtrack, collision, propagationPoint::kToDCA, matchingZ, mBz);
    o2::dataformats::GlobalFwdTrack propmuonAtDCA_Matched = propagateMuon(mchtrack, mchtrack, collision, propagationPoint::kToDCA, matchingZ, mBz);
    o2::dataformats::GlobalFwdTrack propmuonAtPV_Matched = propagateMuon(mchtrack, mchtrack, collision, propagationPoint::kToVertex, matchingZ, mBz);

    float pt = propmuonAtPV.getPt();
    float eta = propmuonAtPV.getEta();
    float phi = propmuonAtPV.getPhi();
    o2::math_utils::bringTo02Pi(phi);

    float cXXatDCA = propmuonAtDCA.getSigma2X();
    float cYYatDCA = propmuonAtDCA.getSigma2Y();
    float cXYatDCA = propmuonAtDCA.getSigmaXY();

    float dcaX = propmuonAtDCA.getX() - collision.posX();
    float dcaY = propmuonAtDCA.getY() - collision.posY();
    float dcaZ = propmuonAtDCA.getZ() - collision.posZ(); // 0 at this point.
    float rAtAbsorberEnd = fwdtrack.rAtAbsorberEnd();     // this works only for GlobalMuonTrack
    float dcaXY = std::sqrt(dcaX * dcaX + dcaY * dcaY);
    float det = cXXatDCA * cYYatDCA - cXYatDCA * cXYatDCA; // determinanat

    float dcaXYinSigma = 999.f;
    if (det < 0) {
      dcaXYinSigma = 999.f;
    } else {
      dcaXYinSigma = std::sqrt(std::fabs((dcaX * dcaX * cYYatDCA + dcaY * dcaY * cXXatDCA - 2. * dcaX * dcaY * cXYatDCA) / det / 2.)); // dca xy in sigma
    }
    float sigma_dcaXY = dcaXY / dcaXYinSigma;

    float dcaX_Matched = propmuonAtDCA_Matched.getX() - collision.posX();
    float dcaY_Matched = propmuonAtDCA_Matched.getY() - collision.posY();
    float dcaXY_Matched = std::sqrt(dcaX_Matched * dcaX_Matched + dcaY_Matched * dcaY_Matched);
    float pDCA = mchtrack.p() * dcaXY_Matched;
    float dphiMP = 999.f, detaMP = 999.f;

    if (refitGlobalMuon) {
      // eta = mfttrack.eta();
      // phi = mfttrack.phi();
      // o2::math_utils::bringTo02Pi(phi);
      eta = propmuonAtDCA.getEta();
      phi = propmuonAtDCA.getPhi();
      o2::math_utils::bringTo02Pi(phi);
      pt = propmuonAtPV_Matched.getP() * std::sin(2.f * std::atan(std::exp(-eta)));

      if constexpr (withMFTCov) {
        // auto mfttrackcov = mftCovs.rawIteratorAt(map_mfttrackcovs[mfttrack.globalIndex()]);
        // o2::track::TrackParCovFwd mftsa = getTrackParCovFwd(mfttrack, mfttrackcov);                                                // values at innermost update
        // o2::dataformats::GlobalFwdTrack globalMuonRefit = o2::aod::fwdtrackutils::refitGlobalMuonCov(propmuonAtPV_Matched, mftsa); // this is track at IU.
        // auto globalMuonAtDCA = o2::aod::fwdtrackutils::propagateTrackParCovFwd(globalMuonRefit, fwdtrack.trackType(), collision, propagationPoint::kToDCA, matchingZ, mBz);

        // eta = globalMuonAtDCA.getEta();
        // phi = globalMuonAtDCA.getPhi();
        // o2::math_utils::bringTo02Pi(phi);
        // pt = globalMuonAtDCA.getPt();
        // p = globalMuonAtDCA.getP();

        // eta = globalMuonRefit.getEta();
        // phi = globalMuonRefit.getPhi();
        // o2::math_utils::bringTo02Pi(phi);
        // pt = globalMuonRefit.getPt();
        // p = globalMuonRefit.getP();

        // cXXatDCA = globalMuonAtDCA.getSigma2X();
        // cYYatDCA = globalMuonAtDCA.getSigma2Y();
        // cXYatDCA = globalMuonAtDCA.getSigmaXY();
        // dcaX = globalMuonAtDCA.getX() - collision.posX();
        // dcaY = globalMuonAtDCA.getY() - collision.posY();
        // dcaXY = std::sqrt(dcaX * dcaX + dcaY * dcaY);
        // det = cXXatDCA * cYYatDCA - cXYatDCA * cXYatDCA; // determinanat
        // if (det < 0) {
        //   dcaXYinSigma = 999.f;
        // } else {
        //   dcaXYinSigma = std::sqrt(std::fabs((dcaX * dcaX * cYYatDCA + dcaY * dcaY * cXXatDCA - 2. * dcaX * dcaY * cXYatDCA) / det / 2.)); // dca xy in sigma
        // }
        // sigma_dcaXY = dcaXY / dcaXYinSigma;

        // cXXatDCA = mftsaAtDCA.getSigma2X();
        // cYYatDCA = mftsaAtDCA.getSigma2Y();
        // cXYatDCA = mftsaAtDCA.getSigmaXY();
        // dcaX = mftsaAtDCA.getX() - collision.posX();
        // dcaY = mftsaAtDCA.getY() - collision.posY();
        // dcaXY = std::sqrt(dcaX * dcaX + dcaY * dcaY);
        // det = cXXatDCA * cYYatDCA - cXYatDCA * cXYatDCA; // determinanat
        // if (det < 0) {
        //   dcaXYinSigma = 999.f;
        // } else {
        //   dcaXYinSigma = std::sqrt(std::fabs((dcaX * dcaX * cYYatDCA + dcaY * dcaY * cXXatDCA - 2. * dcaX * dcaY * cXYatDCA) / det / 2.)); // dca xy in sigma
        // }
        // sigma_dcaXY = dcaXY / dcaXYinSigma;

        getDeltaEtaDeltaPhiAtMatchingPlane<TFwdTracks, TMFTTracks>(collision, fwdtrack, mftCovs, detaMP, dphiMP);
        o2::math_utils::bringToPMPi(dphiMP);
      }
    } else if (propagateToDCAhelix) {
      auto fwdTrackParCov = getTrackParCovFwd(fwdtrack, fwdtrack); // values at innermost update
      std::array<double, 3> dcaInfOrig{999.f, 999.f, 999.f};
      fwdTrackParCov.propagateToDCAhelix(mBz, {collision.posX(), collision.posY(), collision.posZ()}, dcaInfOrig);

      eta = fwdTrackParCov.getEta();
      phi = fwdTrackParCov.getPhi();
      o2::math_utils::bringTo02Pi(phi);
      // pt = fwdTrackParCov.getPt();
      pt = propmuonAtPV_Matched.getP() * std::sin(2.f * std::atan(std::exp(-eta)));

      cXXatDCA = fwdTrackParCov.getSigma2X();
      cYYatDCA = fwdTrackParCov.getSigma2Y();
      cXYatDCA = fwdTrackParCov.getSigmaXY();

      dcaX = fwdTrackParCov.getX() - collision.posX();
      dcaY = fwdTrackParCov.getY() - collision.posY();
      dcaZ = fwdTrackParCov.getZ() - collision.posZ();
      dcaXY = std::sqrt(dcaX * dcaX + dcaY * dcaY);
      det = cXXatDCA * cYYatDCA - cXYatDCA * cXYatDCA; // determinanat

      dcaXYinSigma = 999.f;
      if (det < 0) {
        dcaXYinSigma = 999.f;
      } else {
        dcaXYinSigma = std::sqrt(std::fabs((dcaX * dcaX * cYYatDCA + dcaY * dcaY * cXXatDCA - 2. * dcaX * dcaY * cXYatDCA) / det / 2.)); // dca xy in sigma
      }
      sigma_dcaXY = dcaXY / dcaXYinSigma;
    }

    float ptMatchedMCHMID = propmuonAtPV_Matched.getPt();
    float etaMatchedMCHMID = propmuonAtPV_Matched.getEta();
    float phiMatchedMCHMID = propmuonAtPV_Matched.getPhi();
    o2::math_utils::bringTo02Pi(phiMatchedMCHMID);
    float dpt = (ptMatchedMCHMID - pt) / pt;
    float deta = etaMatchedMCHMID - eta;
    float dphi = phiMatchedMCHMID - phi;
    o2::math_utils::bringToPMPi(dphi);
    if (std::sqrt(std::pow(deta / maxDEta, 2) + std::pow(dphi / maxDPhi, 2)) > 1.f) {
      return;
    }
    if (std::sqrt(std::pow(detaMP / maxDEtaMP, 2) + std::pow(dphiMP / maxDPhiMP, 2)) > 1.f) {
      return;
    }
    if (std::fabs(dpt) > maxRelDPt) {
      return;
    }

    if (!isSelected(pt, eta, rAtAbsorberEnd, pDCA, fwdtrack.chi2() / (2.f * (mchtrack.nClusters() + nClustersMFT) - 5.f), fwdtrack.trackType(), dcaXY)) {
      return;
    }

    if (requireMFTHitMap) {
      std::vector<bool> hasMFTs{hasMFT<0, 1>(mfttrack), hasMFT<2, 3>(mfttrack), hasMFT<4, 5>(mfttrack), hasMFT<6, 7>(mfttrack), hasMFT<8, 9>(mfttrack)};
      for (int i = 0; i < static_cast<int>(requiredMFTDisks->size()); i++) {
        if (!hasMFTs[requiredMFTDisks->at(i)]) {
          return;
        }
      }
    }

    fRegistry.fill(HIST("hMuonType"), fwdtrack.trackType());
    if (isPrimary) {
      if (isMatched) {
        fRegistry.fill(HIST("MFTMCHMID/primary/correct/hPt"), pt);
        fRegistry.fill(HIST("MFTMCHMID/primary/correct/hEtaPhi"), phi, eta);
        fRegistry.fill(HIST("MFTMCHMID/primary/correct/hEtaPhi_MatchedMCHMID"), phiMatchedMCHMID, etaMatchedMCHMID);
        fRegistry.fill(HIST("MFTMCHMID/primary/correct/hsDelta"), pt, dpt, deta, dphi);
        fRegistry.fill(HIST("MFTMCHMID/primary/correct/hsDeltaAtMP"), pt, dphiMP, detaMP);
        fRegistry.fill(HIST("MFTMCHMID/primary/correct/hDiffCollId"), collision.globalIndex() - fwdtrack.collisionId());
        fRegistry.fill(HIST("MFTMCHMID/primary/correct/hSign"), fwdtrack.sign());
        fRegistry.fill(HIST("MFTMCHMID/primary/correct/hNclusters"), fwdtrack.nClusters());
        fRegistry.fill(HIST("MFTMCHMID/primary/correct/hNclustersMFT"), nClustersMFT);
        fRegistry.fill(HIST("MFTMCHMID/primary/correct/hMFTClusterMap"), mftClusterMap(mfttrack));
        fRegistry.fill(HIST("MFTMCHMID/primary/correct/hPDCA_Rabs"), rAtAbsorberEnd, pDCA);
        fRegistry.fill(HIST("MFTMCHMID/primary/correct/hRatAbsorberEnd"), rAtAbsorberEnd);
        fRegistry.fill(HIST("MFTMCHMID/primary/correct/hChi2"), fwdtrack.chi2() / (2.f * (fwdtrack.nClusters() + nClustersMFT) - 5.f));
        fRegistry.fill(HIST("MFTMCHMID/primary/correct/hChi2MFT"), chi2mft);
        fRegistry.fill(HIST("MFTMCHMID/primary/correct/hChi2MatchMCHMID"), fwdtrack.chi2MatchMCHMID());
        fRegistry.fill(HIST("MFTMCHMID/primary/correct/hChi2MatchMCHMFT"), fwdtrack.chi2MatchMCHMFT());
        fRegistry.fill(HIST("MFTMCHMID/primary/correct/hDCAxy2D"), dcaX, dcaY);
        fRegistry.fill(HIST("MFTMCHMID/primary/correct/hDCAxy2DinSigma"), dcaX / std::sqrt(cXXatDCA), dcaY / std::sqrt(cYYatDCA));
        fRegistry.fill(HIST("MFTMCHMID/primary/correct/hDCAxy"), dcaXY);
        fRegistry.fill(HIST("MFTMCHMID/primary/correct/hDCAz"), dcaZ);
        fRegistry.fill(HIST("MFTMCHMID/primary/correct/hDCAxyz"), dcaXY, dcaZ);
        fRegistry.fill(HIST("MFTMCHMID/primary/correct/hDCAxyinSigma"), dcaXYinSigma);
        fRegistry.fill(HIST("MFTMCHMID/primary/correct/hMCHBitMap"), fwdtrack.mchBitMap());
        fRegistry.fill(HIST("MFTMCHMID/primary/correct/hMIDBitMap"), fwdtrack.midBitMap());
        fRegistry.fill(HIST("MFTMCHMID/primary/correct/hDCAxResolutionvsPt"), pt, std::sqrt(cXXatDCA) * 1e+4); // convert cm to um
        fRegistry.fill(HIST("MFTMCHMID/primary/correct/hDCAyResolutionvsPt"), pt, std::sqrt(cYYatDCA) * 1e+4); // convert cm to um
        fRegistry.fill(HIST("MFTMCHMID/primary/correct/hDCAxyResolutionvsPt"), pt, sigma_dcaXY * 1e+4);        // convert cm to um
        fRegistry.fill(HIST("MFTMCHMID/primary/correct/hProdVtxZ"), mcParticle_MFTMCHMID.vz());

        fRegistry.fill(HIST("MFTMCHMID/primary/correct/hRelDeltaPt"), mcParticle_MFTMCHMID.pt(), (pt - mcParticle_MFTMCHMID.pt()) / mcParticle_MFTMCHMID.pt());
        if (mcParticle_MFTMCHMID.pdgCode() > 0) {
          fRegistry.fill(HIST("MFTMCHMID/primary/correct/hDeltaEta_Neg"), mcParticle_MFTMCHMID.pt(), eta - mcParticle_MFTMCHMID.eta());
          fRegistry.fill(HIST("MFTMCHMID/primary/correct/hDeltaPhi_Neg"), mcParticle_MFTMCHMID.pt(), phi - mcParticle_MFTMCHMID.phi());
        } else {
          fRegistry.fill(HIST("MFTMCHMID/primary/correct/hDeltaEta_Pos"), mcParticle_MFTMCHMID.pt(), eta - mcParticle_MFTMCHMID.eta());
          fRegistry.fill(HIST("MFTMCHMID/primary/correct/hDeltaPhi_Pos"), mcParticle_MFTMCHMID.pt(), phi - mcParticle_MFTMCHMID.phi());
        }
      } else {
        fRegistry.fill(HIST("MFTMCHMID/primary/wrong/hPt"), pt);
        fRegistry.fill(HIST("MFTMCHMID/primary/wrong/hEtaPhi"), phi, eta);
        fRegistry.fill(HIST("MFTMCHMID/primary/wrong/hEtaPhi_MatchedMCHMID"), phiMatchedMCHMID, etaMatchedMCHMID);
        fRegistry.fill(HIST("MFTMCHMID/primary/wrong/hsDelta"), pt, dpt, deta, dphi);
        fRegistry.fill(HIST("MFTMCHMID/primary/wrong/hsDeltaAtMP"), pt, dphiMP, detaMP);
        fRegistry.fill(HIST("MFTMCHMID/primary/wrong/hDiffCollId"), collision.globalIndex() - fwdtrack.collisionId());
        fRegistry.fill(HIST("MFTMCHMID/primary/wrong/hSign"), fwdtrack.sign());
        fRegistry.fill(HIST("MFTMCHMID/primary/wrong/hNclusters"), fwdtrack.nClusters());
        fRegistry.fill(HIST("MFTMCHMID/primary/wrong/hNclustersMFT"), nClustersMFT);
        fRegistry.fill(HIST("MFTMCHMID/primary/wrong/hMFTClusterMap"), mftClusterMap(mfttrack));
        fRegistry.fill(HIST("MFTMCHMID/primary/wrong/hPDCA_Rabs"), rAtAbsorberEnd, pDCA);
        fRegistry.fill(HIST("MFTMCHMID/primary/wrong/hRatAbsorberEnd"), rAtAbsorberEnd);
        fRegistry.fill(HIST("MFTMCHMID/primary/wrong/hChi2"), fwdtrack.chi2() / (2.f * (fwdtrack.nClusters() + nClustersMFT) - 5.f));
        fRegistry.fill(HIST("MFTMCHMID/primary/wrong/hChi2MFT"), chi2mft);
        fRegistry.fill(HIST("MFTMCHMID/primary/wrong/hChi2MatchMCHMID"), fwdtrack.chi2MatchMCHMID());
        fRegistry.fill(HIST("MFTMCHMID/primary/wrong/hChi2MatchMCHMFT"), fwdtrack.chi2MatchMCHMFT());
        fRegistry.fill(HIST("MFTMCHMID/primary/wrong/hDCAxy2D"), dcaX, dcaY);
        fRegistry.fill(HIST("MFTMCHMID/primary/wrong/hDCAxy2DinSigma"), dcaX / std::sqrt(cXXatDCA), dcaY / std::sqrt(cYYatDCA));
        fRegistry.fill(HIST("MFTMCHMID/primary/wrong/hDCAxy"), dcaXY);
        fRegistry.fill(HIST("MFTMCHMID/primary/wrong/hDCAz"), dcaZ);
        fRegistry.fill(HIST("MFTMCHMID/primary/wrong/hDCAxyz"), dcaXY, dcaZ);
        fRegistry.fill(HIST("MFTMCHMID/primary/wrong/hDCAxyinSigma"), dcaXYinSigma);
        fRegistry.fill(HIST("MFTMCHMID/primary/wrong/hMCHBitMap"), fwdtrack.mchBitMap());
        fRegistry.fill(HIST("MFTMCHMID/primary/wrong/hMIDBitMap"), fwdtrack.midBitMap());
        fRegistry.fill(HIST("MFTMCHMID/primary/wrong/hDCAxResolutionvsPt"), pt, std::sqrt(cXXatDCA) * 1e+4); // convert cm to um
        fRegistry.fill(HIST("MFTMCHMID/primary/wrong/hDCAyResolutionvsPt"), pt, std::sqrt(cYYatDCA) * 1e+4); // convert cm to um
        fRegistry.fill(HIST("MFTMCHMID/primary/wrong/hDCAxyResolutionvsPt"), pt, sigma_dcaXY * 1e+4);        // convert cm to um
        fRegistry.fill(HIST("MFTMCHMID/primary/wrong/hProdVtxZ"), mcParticle_MFTMCHMID.vz());

        fRegistry.fill(HIST("MFTMCHMID/primary/wrong/hRelDeltaPt"), mcParticle_MFTMCHMID.pt(), (pt - mcParticle_MFTMCHMID.pt()) / mcParticle_MFTMCHMID.pt());
        if (mcParticle_MFTMCHMID.pdgCode() > 0) {
          fRegistry.fill(HIST("MFTMCHMID/primary/wrong/hDeltaEta_Neg"), mcParticle_MFTMCHMID.pt(), eta - mcParticle_MFTMCHMID.eta());
          fRegistry.fill(HIST("MFTMCHMID/primary/wrong/hDeltaPhi_Neg"), mcParticle_MFTMCHMID.pt(), phi - mcParticle_MFTMCHMID.phi());
        } else {
          fRegistry.fill(HIST("MFTMCHMID/primary/wrong/hDeltaEta_Pos"), mcParticle_MFTMCHMID.pt(), eta - mcParticle_MFTMCHMID.eta());
          fRegistry.fill(HIST("MFTMCHMID/primary/wrong/hDeltaPhi_Pos"), mcParticle_MFTMCHMID.pt(), phi - mcParticle_MFTMCHMID.phi());
        }
      }
    } else { // secondary
      if (isMatched) {
        fRegistry.fill(HIST("MFTMCHMID/secondary/correct/hPt"), pt);
        fRegistry.fill(HIST("MFTMCHMID/secondary/correct/hEtaPhi"), phi, eta);
        fRegistry.fill(HIST("MFTMCHMID/secondary/correct/hEtaPhi_MatchedMCHMID"), phiMatchedMCHMID, etaMatchedMCHMID);
        fRegistry.fill(HIST("MFTMCHMID/secondary/correct/hsDelta"), pt, dpt, deta, dphi);
        fRegistry.fill(HIST("MFTMCHMID/secondary/correct/hsDeltaAtMP"), pt, dphiMP, detaMP);
        fRegistry.fill(HIST("MFTMCHMID/secondary/correct/hDiffCollId"), collision.globalIndex() - fwdtrack.collisionId());
        fRegistry.fill(HIST("MFTMCHMID/secondary/correct/hSign"), fwdtrack.sign());
        fRegistry.fill(HIST("MFTMCHMID/secondary/correct/hNclusters"), fwdtrack.nClusters());
        fRegistry.fill(HIST("MFTMCHMID/secondary/correct/hNclustersMFT"), nClustersMFT);
        fRegistry.fill(HIST("MFTMCHMID/secondary/correct/hMFTClusterMap"), mftClusterMap(mfttrack));
        fRegistry.fill(HIST("MFTMCHMID/secondary/correct/hPDCA_Rabs"), rAtAbsorberEnd, pDCA);
        fRegistry.fill(HIST("MFTMCHMID/secondary/correct/hRatAbsorberEnd"), rAtAbsorberEnd);
        fRegistry.fill(HIST("MFTMCHMID/secondary/correct/hChi2"), fwdtrack.chi2() / (2.f * (fwdtrack.nClusters() + nClustersMFT) - 5.f));
        fRegistry.fill(HIST("MFTMCHMID/secondary/correct/hChi2MFT"), chi2mft);
        fRegistry.fill(HIST("MFTMCHMID/secondary/correct/hChi2MatchMCHMID"), fwdtrack.chi2MatchMCHMID());
        fRegistry.fill(HIST("MFTMCHMID/secondary/correct/hChi2MatchMCHMFT"), fwdtrack.chi2MatchMCHMFT());
        fRegistry.fill(HIST("MFTMCHMID/secondary/correct/hDCAxy2D"), dcaX, dcaY);
        fRegistry.fill(HIST("MFTMCHMID/secondary/correct/hDCAxy2DinSigma"), dcaX / std::sqrt(cXXatDCA), dcaY / std::sqrt(cYYatDCA));
        fRegistry.fill(HIST("MFTMCHMID/secondary/correct/hDCAxy"), dcaXY);
        fRegistry.fill(HIST("MFTMCHMID/secondary/correct/hDCAz"), dcaZ);
        fRegistry.fill(HIST("MFTMCHMID/secondary/correct/hDCAxyz"), dcaXY, dcaZ);
        fRegistry.fill(HIST("MFTMCHMID/secondary/correct/hDCAxyinSigma"), dcaXYinSigma);
        fRegistry.fill(HIST("MFTMCHMID/secondary/correct/hMCHBitMap"), fwdtrack.mchBitMap());
        fRegistry.fill(HIST("MFTMCHMID/secondary/correct/hMIDBitMap"), fwdtrack.midBitMap());
        fRegistry.fill(HIST("MFTMCHMID/secondary/correct/hDCAxResolutionvsPt"), pt, std::sqrt(cXXatDCA) * 1e+4); // convert cm to um
        fRegistry.fill(HIST("MFTMCHMID/secondary/correct/hDCAyResolutionvsPt"), pt, std::sqrt(cYYatDCA) * 1e+4); // convert cm to um
        fRegistry.fill(HIST("MFTMCHMID/secondary/correct/hDCAxyResolutionvsPt"), pt, sigma_dcaXY * 1e+4);        // convert cm to um
        fRegistry.fill(HIST("MFTMCHMID/secondary/correct/hProdVtxZ"), mcParticle_MFTMCHMID.vz());

        fRegistry.fill(HIST("MFTMCHMID/secondary/correct/hRelDeltaPt"), mcParticle_MFTMCHMID.pt(), (pt - mcParticle_MFTMCHMID.pt()) / mcParticle_MFTMCHMID.pt());
        if (mcParticle_MFTMCHMID.pdgCode() > 0) {
          fRegistry.fill(HIST("MFTMCHMID/secondary/correct/hDeltaEta_Neg"), mcParticle_MFTMCHMID.pt(), eta - mcParticle_MFTMCHMID.eta());
          fRegistry.fill(HIST("MFTMCHMID/secondary/correct/hDeltaPhi_Neg"), mcParticle_MFTMCHMID.pt(), phi - mcParticle_MFTMCHMID.phi());
        } else {
          fRegistry.fill(HIST("MFTMCHMID/secondary/correct/hDeltaEta_Pos"), mcParticle_MFTMCHMID.pt(), eta - mcParticle_MFTMCHMID.eta());
          fRegistry.fill(HIST("MFTMCHMID/secondary/correct/hDeltaPhi_Pos"), mcParticle_MFTMCHMID.pt(), phi - mcParticle_MFTMCHMID.phi());
        }
      } else {
        fRegistry.fill(HIST("MFTMCHMID/secondary/wrong/hPt"), pt);
        fRegistry.fill(HIST("MFTMCHMID/secondary/wrong/hEtaPhi"), phi, eta);
        fRegistry.fill(HIST("MFTMCHMID/secondary/wrong/hEtaPhi_MatchedMCHMID"), phiMatchedMCHMID, etaMatchedMCHMID);
        fRegistry.fill(HIST("MFTMCHMID/secondary/wrong/hsDelta"), pt, dpt, deta, dphi);
        fRegistry.fill(HIST("MFTMCHMID/secondary/wrong/hsDeltaAtMP"), pt, dphiMP, detaMP);
        fRegistry.fill(HIST("MFTMCHMID/secondary/wrong/hDiffCollId"), collision.globalIndex() - fwdtrack.collisionId());
        fRegistry.fill(HIST("MFTMCHMID/secondary/wrong/hSign"), fwdtrack.sign());
        fRegistry.fill(HIST("MFTMCHMID/secondary/wrong/hNclusters"), fwdtrack.nClusters());
        fRegistry.fill(HIST("MFTMCHMID/secondary/wrong/hNclustersMFT"), nClustersMFT);
        fRegistry.fill(HIST("MFTMCHMID/secondary/wrong/hMFTClusterMap"), mftClusterMap(mfttrack));
        fRegistry.fill(HIST("MFTMCHMID/secondary/wrong/hPDCA_Rabs"), rAtAbsorberEnd, pDCA);
        fRegistry.fill(HIST("MFTMCHMID/secondary/wrong/hRatAbsorberEnd"), rAtAbsorberEnd);
        fRegistry.fill(HIST("MFTMCHMID/secondary/wrong/hChi2"), fwdtrack.chi2() / (2.f * (fwdtrack.nClusters() + nClustersMFT) - 5.f));
        fRegistry.fill(HIST("MFTMCHMID/secondary/wrong/hChi2MFT"), chi2mft);
        fRegistry.fill(HIST("MFTMCHMID/secondary/wrong/hChi2MatchMCHMID"), fwdtrack.chi2MatchMCHMID());
        fRegistry.fill(HIST("MFTMCHMID/secondary/wrong/hChi2MatchMCHMFT"), fwdtrack.chi2MatchMCHMFT());
        fRegistry.fill(HIST("MFTMCHMID/secondary/wrong/hDCAxy2D"), dcaX, dcaY);
        fRegistry.fill(HIST("MFTMCHMID/secondary/wrong/hDCAxy2DinSigma"), dcaX / std::sqrt(cXXatDCA), dcaY / std::sqrt(cYYatDCA));
        fRegistry.fill(HIST("MFTMCHMID/secondary/wrong/hDCAxy"), dcaXY);
        fRegistry.fill(HIST("MFTMCHMID/secondary/wrong/hDCAz"), dcaZ);
        fRegistry.fill(HIST("MFTMCHMID/secondary/wrong/hDCAxyz"), dcaXY, dcaZ);
        fRegistry.fill(HIST("MFTMCHMID/secondary/wrong/hDCAxyinSigma"), dcaXYinSigma);
        fRegistry.fill(HIST("MFTMCHMID/secondary/wrong/hMCHBitMap"), fwdtrack.mchBitMap());
        fRegistry.fill(HIST("MFTMCHMID/secondary/wrong/hMIDBitMap"), fwdtrack.midBitMap());
        fRegistry.fill(HIST("MFTMCHMID/secondary/wrong/hDCAxResolutionvsPt"), pt, std::sqrt(cXXatDCA) * 1e+4); // convert cm to um
        fRegistry.fill(HIST("MFTMCHMID/secondary/wrong/hDCAyResolutionvsPt"), pt, std::sqrt(cYYatDCA) * 1e+4); // convert cm to um
        fRegistry.fill(HIST("MFTMCHMID/secondary/wrong/hDCAxyResolutionvsPt"), pt, sigma_dcaXY * 1e+4);        // convert cm to um
        fRegistry.fill(HIST("MFTMCHMID/secondary/wrong/hProdVtxZ"), mcParticle_MFTMCHMID.vz());

        fRegistry.fill(HIST("MFTMCHMID/secondary/wrong/hRelDeltaPt"), mcParticle_MFTMCHMID.pt(), (pt - mcParticle_MFTMCHMID.pt()) / mcParticle_MFTMCHMID.pt());
        if (mcParticle_MFTMCHMID.pdgCode() > 0) {
          fRegistry.fill(HIST("MFTMCHMID/secondary/wrong/hDeltaEta_Neg"), mcParticle_MFTMCHMID.pt(), eta - mcParticle_MFTMCHMID.eta());
          fRegistry.fill(HIST("MFTMCHMID/secondary/wrong/hDeltaPhi_Neg"), mcParticle_MFTMCHMID.pt(), phi - mcParticle_MFTMCHMID.phi());
        } else {
          fRegistry.fill(HIST("MFTMCHMID/secondary/wrong/hDeltaEta_Pos"), mcParticle_MFTMCHMID.pt(), eta - mcParticle_MFTMCHMID.eta());
          fRegistry.fill(HIST("MFTMCHMID/secondary/wrong/hDeltaPhi_Pos"), mcParticle_MFTMCHMID.pt(), phi - mcParticle_MFTMCHMID.phi());
        }
      }
    }
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

  std::vector<std::tuple<int, int, int>> vec_min_chi2MatchMCHMFT; // std::pair<globalIndex of global muon, globalIndex of matched MCH-MID, globalIndex of MFT> -> chi2MatchMCHMFT;
  std::vector<std::tuple<int, int, int>> vec_min_dr;              // std::pair<globalIndex of global muon, globalIndex of matched MCH-MID, globalIndex of MFT> -> deta + dphi;
  std::map<std::tuple<int, int, int>, bool> mapCorrectMatch;

  template <bool withMFTCov, typename TCollision, typename TFwdTrack, typename TFwdTracks, typename TMFTTracks, typename TMFTTracksCov>
  void findBestMatchPerMCHMID(TCollision const& collision, TFwdTrack const& fwdtrack, TFwdTracks const& fwdtracks, TMFTTracks const&, TMFTTracksCov const& mftCovs)
  {
    if (fwdtrack.trackType() != o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack) {
      return;
    }
    if (!fwdtrack.has_mcParticle()) {
      return;
    }

    std::tuple<int, int, int> tupleIds_at_min_chi2mftmch;
    std::tuple<int, int, int> tupleIds_at_min_dr;
    float min_chi2MatchMCHMFT = 1e+10;
    float min_dr = 1e+10;
    auto muons_per_MCHMID = fwdtracks.sliceBy(fwdtracksPerMCHTrack, fwdtrack.globalIndex());
    // LOGF(info, "muons_per_MCHMID.size() = %d", muons_per_MCHMID.size());

    for (const auto& muon_tmp : muons_per_MCHMID) {
      if (muon_tmp.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack) {

        auto tupleId = std::make_tuple(muon_tmp.globalIndex(), muon_tmp.matchMCHTrackId(), muon_tmp.matchMFTTrackId());
        auto mchtrack = muon_tmp.template matchMCHTrack_as<MyFwdTracks>(); // MCH-MID
        auto mfttrack = muon_tmp.template matchMFTTrack_as<MyMFTTracks>();

        if (!muon_tmp.has_mcParticle() || !mchtrack.has_mcParticle() || !mfttrack.has_mcParticle()) {
          continue;
        }

        float deta = 999.f, dphi = 999.f;
        getDeltaEtaDeltaPhiAtMatchingPlane<TFwdTracks, TMFTTracks>(collision, muon_tmp, mftCovs, deta, dphi);
        float dr = std::sqrt(deta * deta + dphi * dphi);

        // auto mcParticle_MFTMCHMID = muon_tmp.template mcParticle_as<aod::McParticles>(); // this is identical to mcParticle_MCHMID
        auto mcParticle_MCHMID = mchtrack.template mcParticle_as<aod::McParticles>(); // this is identical to mcParticle_MFTMCHMID
        auto mcParticle_MFT = mfttrack.template mcParticle_as<aod::McParticles>();
        // float chi2ndf = muon_tmp.chi2() / (2.f * (mchtrack.nClusters() + mfttrack.nClusters()) - 5.f);

        if (mcParticle_MFT.globalIndex() == mcParticle_MCHMID.globalIndex()) {
          mapCorrectMatch[tupleId] = true;
        } else {
          mapCorrectMatch[tupleId] = false;
        }

        // if (std::abs(mcParticle_MCHMID.pdgCode()) == 13 && mcParticle_MCHMID.isPhysicalPrimary()) {
        //   if (mcParticle_MFT.globalIndex() == mcParticle_MCHMID.globalIndex()) {
        //     LOGF(info, "This is correct match between MFT and MCH-MID: muon_tmp.globalIndex() = %d, chi2/ndf = %f, matching chi2/ndf = %f, mcParticle.pt() = %f, mcParticle.eta() = %f, mcParticle.phi() = %f, reldpt = %f, deta = %f, dphi = %f, dr = %f", muon_tmp.globalIndex(), chi2ndf, muon_tmp.chi2MatchMCHMFT(),  mcParticle_MCHMID.pt(), mcParticle_MCHMID.eta(), mcParticle_MCHMID.phi(), reldpt, deta, dphi, dr);
        //   } else {
        //     LOGF(info, "This is wrong match between MFT and MCH-MID: muon_tmp.globalIndex() = %d, chi2/ndf = %f, matching chi2/ndf = %f , mcParticle.pt() = %f, mcParticle.eta() = %f, mcParticle.phi() = %f, reldpt = %f, deta = %f, dphi = %f, dr = %f" , muon_tmp.globalIndex(), chi2ndf, muon_tmp.chi2MatchMCHMFT(), mcParticle_MCHMID.pt(), mcParticle_MCHMID.eta(), mcParticle_MCHMID.phi(), reldpt, deta, dphi, dr);
        //   }
        // }

        if (0.f < muon_tmp.chi2MatchMCHMFT() && muon_tmp.chi2MatchMCHMFT() < min_chi2MatchMCHMFT) {
          min_chi2MatchMCHMFT = muon_tmp.chi2MatchMCHMFT();
          tupleIds_at_min_chi2mftmch = tupleId;
        }

        if (dr < min_dr) {
          min_dr = dr;
          tupleIds_at_min_dr = tupleId;
        }

      } // end of if global muon
    } // end of candidates loop

    vec_min_chi2MatchMCHMFT.emplace_back(tupleIds_at_min_chi2mftmch);
    vec_min_dr.emplace_back(tupleIds_at_min_dr);

    // auto mcParticleTMP = fwdtrack.template mcParticle_as<aod::McParticles>(); // this is identical to mcParticle_MFTMCHMID
    // if (std::abs(mcParticleTMP.pdgCode()) == 13 && mcParticleTMP.isPhysicalPrimary()) {
    //   LOGF(info, "min chi2: muon_tmp.globalIndex() = %d, muon_tmp.matchMCHTrackId() = %d, muon_tmp.matchMFTTrackId() = %d, muon_tmp.chi2MatchMCHMFT() = %f, correct match = %d", std::get<0>(tupleIds_at_min_chi2mftmch), std::get<1>(tupleIds_at_min_chi2mftmch), std::get<2>(tupleIds_at_min_chi2mftmch), min_chi2MatchMCHMFT, mapCorrectMatch[tupleIds_at_min_chi2mftmch]);
    //   LOGF(info, "min dr: muon_tmp.globalIndex() = %d, muon_tmp.matchMCHTrackId() = %d, muon_tmp.matchMFTTrackId() = %d, dr = %f, correct match = %d", std::get<0>(tupleIds_at_min_dr), std::get<1>(tupleIds_at_min_dr), std::get<2>(tupleIds_at_min_dr), min_dr, mapCorrectMatch[tupleIds_at_min_dr]);
    // }
  }

  SliceCache cache;
  PresliceUnsorted<aod::FwdTracks> fwdtracksPerMCHTrack = aod::fwdtrack::matchMCHTrackId;
  PresliceUnsorted<aod::FwdTracks> perMFTTrack = o2::aod::fwdtrack::matchMFTTrackId;
  Preslice<aod::FwdTracks> perCollision = o2::aod::fwdtrack::collisionId;
  Preslice<aod::FwdTrackAssoc> fwdtrackIndicesPerCollision = aod::track_association::collisionId;
  PresliceUnsorted<aod::FwdTrackAssoc> fwdtrackIndicesPerFwdTrack = aod::track_association::fwdtrackId;

  Filter collisionFilter_evsel = o2::aod::evsel::sel8 == true && (cfgZvtxMin < o2::aod::collision::posZ && o2::aod::collision::posZ < cfgZvtxMax);
  Filter collisionFilter_centrality = (cfgCentMin < o2::aod::cent::centFT0M && o2::aod::cent::centFT0M < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0A && o2::aod::cent::centFT0A < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0C && o2::aod::cent::centFT0C < cfgCentMax);
  using FilteredMyCollisions = soa::Filtered<MyCollisions>;

  void processWithoutFTTCA(FilteredMyCollisions const& collisions, MyFwdTracks const& fwdtracks, MyMFTTracks const& mfttracks, aod::BCsWithTimestamps const&, aod::McParticles const&)
  {
    for (const auto& collision : collisions) {
      auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);
      fRegistry.fill(HIST("Event/hCollisionCounter"), 0);
      if (!collision.has_mcCollision()) {
        continue;
      }
      if (cfgRequireGoodRCT && !rctChecker.checkTable(collision)) {
        continue;
      }
      float centralities[3] = {collision.centFT0M(), collision.centFT0A(), collision.centFT0C()};
      if (centralities[cfgCentEstimator] < cfgCentMin || cfgCentMax < centralities[cfgCentEstimator]) {
        continue;
      }
      fRegistry.fill(HIST("Event/hCollisionCounter"), 1);
      fillEventHistograms(collision);

      auto fwdtracks_per_coll = fwdtracks.sliceBy(perCollision, collision.globalIndex());

      for (const auto& fwdtrack : fwdtracks_per_coll) {
        if (fwdtrack.trackType() != o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack) {
          continue;
        }
        fillHistograms<false>(collision, fwdtrack, fwdtracks, mfttracks, nullptr);
      } // end of fwdtrack loop
    } // end of collision loop

    vec_min_chi2MatchMCHMFT.clear();
    vec_min_chi2MatchMCHMFT.shrink_to_fit();
    vec_min_dr.clear();
    vec_min_dr.shrink_to_fit();
  }
  PROCESS_SWITCH(matchingMFT, processWithoutFTTCA, "process without FTTCA", false);

  void processWithFTTCA(FilteredMyCollisions const& collisions, MyFwdTracks const& fwdtracks, MyMFTTracks const& mfttracks, aod::BCsWithTimestamps const&, aod::FwdTrackAssoc const& fwdtrackIndices, aod::McParticles const&)
  {
    for (const auto& collision : collisions) {
      auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);
      fRegistry.fill(HIST("Event/hCollisionCounter"), 0);
      if (!collision.has_mcCollision()) {
        continue;
      }
      if (cfgRequireGoodRCT && !rctChecker.checkTable(collision)) {
        continue;
      }
      float centralities[3] = {collision.centFT0M(), collision.centFT0A(), collision.centFT0C()};
      if (centralities[cfgCentEstimator] < cfgCentMin || cfgCentMax < centralities[cfgCentEstimator]) {
        continue;
      }
      fRegistry.fill(HIST("Event/hCollisionCounter"), 1);
      fillEventHistograms(collision);

      auto fwdtrackIdsThisCollision = fwdtrackIndices.sliceBy(fwdtrackIndicesPerCollision, collision.globalIndex());
      for (const auto& fwdtrackId : fwdtrackIdsThisCollision) {
        auto fwdtrack = fwdtrackId.template fwdtrack_as<MyFwdTracks>();
        if (fwdtrack.trackType() != o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack) {
          continue;
        }
        fillHistograms<false>(collision, fwdtrack, fwdtracks, mfttracks, nullptr);
      } // end of fwdtrack loop
    } // end of collision loop

    vec_min_chi2MatchMCHMFT.clear();
    vec_min_chi2MatchMCHMFT.shrink_to_fit();
    vec_min_dr.clear();
    vec_min_dr.shrink_to_fit();
  }
  PROCESS_SWITCH(matchingMFT, processWithFTTCA, "process with FTTCA", true);

  std::unordered_map<int, int> map_mfttrackcovs;

  void processWithFTTCA_withMFTCov(FilteredMyCollisions const& collisions, MyFwdTracks const& fwdtracks, MyMFTTracks const& mfttracks, aod::MFTTracksCov const& mftCovs, aod::BCsWithTimestamps const&, aod::FwdTrackAssoc const& fwdtrackIndices, aod::McParticles const&)
  {
    for (const auto& mfttrackConv : mftCovs) {
      map_mfttrackcovs[mfttrackConv.matchMFTTrackId()] = mfttrackConv.globalIndex();
    }

    vec_min_chi2MatchMCHMFT.reserve(fwdtracks.size());
    vec_min_dr.reserve(fwdtracks.size());

    for (const auto& collision : collisions) {
      auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);
      auto fwdtrackIdsThisCollision = fwdtrackIndices.sliceBy(fwdtrackIndicesPerCollision, collision.globalIndex());
      for (const auto& fwdtrackId : fwdtrackIdsThisCollision) {
        auto fwdtrack = fwdtrackId.template fwdtrack_as<MyFwdTracks>();
        findBestMatchPerMCHMID<true>(collision, fwdtrack, fwdtracks, mfttracks, mftCovs);
      } // end of fwdtrack loop
    } // end of collision loop

    for (const auto& collision : collisions) {
      auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);
      fRegistry.fill(HIST("Event/hCollisionCounter"), 0);
      if (!collision.has_mcCollision()) {
        continue;
      }
      if (cfgRequireGoodRCT && !rctChecker.checkTable(collision)) {
        continue;
      }
      float centralities[3] = {collision.centFT0M(), collision.centFT0A(), collision.centFT0C()};
      if (centralities[cfgCentEstimator] < cfgCentMin || cfgCentMax < centralities[cfgCentEstimator]) {
        continue;
      }
      fRegistry.fill(HIST("Event/hCollisionCounter"), 1);
      fillEventHistograms(collision);

      auto fwdtrackIdsThisCollision = fwdtrackIndices.sliceBy(fwdtrackIndicesPerCollision, collision.globalIndex());
      for (const auto& fwdtrackId : fwdtrackIdsThisCollision) {
        auto fwdtrack = fwdtrackId.template fwdtrack_as<MyFwdTracks>();

        if (cfgBestMatchFinder == 0) { // chi2
          if (fwdtrack.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack && std::find(vec_min_chi2MatchMCHMFT.begin(), vec_min_chi2MatchMCHMFT.end(), std::make_tuple(fwdtrack.globalIndex(), fwdtrack.matchMCHTrackId(), fwdtrack.matchMFTTrackId())) == vec_min_chi2MatchMCHMFT.end()) {
            continue;
          }
        } else if (cfgBestMatchFinder == 1) { // dr
          if (fwdtrack.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack && std::find(vec_min_dr.begin(), vec_min_dr.end(), std::make_tuple(fwdtrack.globalIndex(), fwdtrack.matchMCHTrackId(), fwdtrack.matchMFTTrackId())) == vec_min_dr.end()) {
            continue;
          }
        } else { // best match is not selected. Histograms are filled with all global muons.
          if (fwdtrack.trackType() != o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack) {
            continue;
          }
        }
        fillHistograms<true>(collision, fwdtrack, fwdtracks, mfttracks, mftCovs);
      } // end of fwdtrack loop
    } // end of collision loop

    map_mfttrackcovs.clear();
    vec_min_chi2MatchMCHMFT.clear();
    vec_min_chi2MatchMCHMFT.shrink_to_fit();
    vec_min_dr.clear();
    vec_min_dr.shrink_to_fit();
  }
  PROCESS_SWITCH(matchingMFT, processWithFTTCA_withMFTCov, "process with FTTCA with MFTCov", false);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<matchingMFT>(cfgc, TaskName{"matching-mft"})};
}
