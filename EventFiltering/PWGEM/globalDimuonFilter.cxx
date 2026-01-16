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

// \brief software trigger for global dimuons
// \author daiki.sekihata@cern.ch

#include "Common/Core/TableHelper.h"
#include "Common/Core/fwdtrackUtilities.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/DataModel/EventSelection.h"
#include "EventFiltering/filterTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "CommonConstants/PhysicsConstants.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DetectorsBase/Propagator.h"
#include "Field/MagneticField.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/DataTypes.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "GlobalTracking/MatchGlobalFwd.h"
#include "MCHTracking/TrackExtrap.h"
#include "MCHTracking/TrackParam.h"
#include "ReconstructionDataFormats/TrackFwd.h"

#include "Math/Vector4D.h"

#include <algorithm>
#include <map>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

using namespace o2;
using namespace o2::soa;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;
using namespace o2::aod::fwdtrackutils;

struct globalDimuonFilter {
  Produces<aod::GlobalDimuonFilters> tags;

  // Configurables
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
  Configurable<bool> fillQAHistograms{"fillQAHistograms", true, "flag to fill QA histograms"};
  Configurable<int> minNmuon{"minNmuon", 2, "min number of global muon candidates per collision"};

  struct : ConfigurableGroup {
    std::string prefix = "eventCutGroup";
    Configurable<float> minZvtx{"minZvtx", -10.f, "min. Zvtx of collision"};
    Configurable<float> maxZvtx{"maxZvtx", +10.f, "max. Zvtx of collision"};
  } eventCutGroup;

  struct : ConfigurableGroup {
    std::string prefix = "glMuonCutGroup";
    Configurable<float> minPt{"minPt", 0.01, "min pt for muon"};
    Configurable<float> maxPt{"maxPt", 1e+10, "max pt for muon"};
    Configurable<float> minEta{"minEta", -3.6, "min. eta acceptance for MFT-MCH-MID"};
    Configurable<float> maxEta{"maxEta", -2.5, "max. eta acceptance for MFT-MCH-MID"};
    Configurable<float> minRabs{"minRabs", 17.6, "min. R at absorber end for global muon (min. eta = -3.6)"}; // std::tan(2.f * std::atan(std::exp(- -3.6)) ) * -505. = 27.6
    // Configurable<float> midRabs{"midRabs", 26.5, "middle R at absorber end for pDCA cut"};
    Configurable<float> maxRabs{"maxRabs", 89.5, "max. R at absorber end"};
    Configurable<float> maxDCAxy{"maxDCAxy", 1.0, "max. DCAxy for global muons"};
    // Configurable<float> maxPDCAforLargeR{"maxPDCAforLargeR", 324.f, "max. pDCA for large R at absorber end"};
    // Configurable<float> maxPDCAforSmallR{"maxPDCAforSmallR", 594.f, "max. pDCA for small R at absorber end"};
    Configurable<float> maxMatchingChi2MCHMFT{"maxMatchingChi2MCHMFT", 100.f, "max. chi2 for MCH-MFT matching"};
    Configurable<float> maxChi2{"maxChi2", 10.f, "max. chi2/ndf for global muon"};
    Configurable<float> maxChi2MFT{"maxChi2MFT", 1e+10, "max. chi2/ndf for MFTsa"};
    Configurable<bool> refitGlobalMuon{"refitGlobalMuon", true, "flag to refit global muon"};
    Configurable<float> matchingZ{"matchingZ", -77.5, "z position where matching is performed"};
    Configurable<float> maxDEta{"maxDEta", 0.2, "max. deta between MFT-MCH-MID and MCH-MID"};
    Configurable<float> maxDPhi{"maxDPhi", 0.2, "max. dphi between MFT-MCH-MID and MCH-MID"};
  } glMuonCutGroup;

  HistogramRegistry fRegistry{"output", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};
  o2::ccdb::CcdbApi ccdbApi;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  int mRunNumber = 0;
  float mBz = 0;

  void init(o2::framework::InitContext&)
  {
    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);
    ccdbApi.init(ccdburl);
    mRunNumber = 0;
    mBz = 0;

    addHistograms();
  }

  template <typename TBC>
  void initCCDB(TBC const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }
    mRunNumber = bc.runNumber();

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
    auto hCollisionCounter{std::get<std::shared_ptr<TH1>>(fRegistry.add("hCollisionCounter", "Number of collisions", HistType::kTH1D, {{10, 0.5, 10.5}}))};
    hCollisionCounter->GetXaxis()->SetBinLabel(1, "all");
    hCollisionCounter->GetXaxis()->SetBinLabel(2, "TVX");
    hCollisionCounter->GetXaxis()->SetBinLabel(3, "No TFB");
    hCollisionCounter->GetXaxis()->SetBinLabel(4, "No ITS ROFB");
    hCollisionCounter->GetXaxis()->SetBinLabel(5, "|Z_{vtx}| < 10 cm");
    hCollisionCounter->GetXaxis()->SetBinLabel(6, "sel8");
    hCollisionCounter->GetXaxis()->SetBinLabel(7, "sel8 && |Z_{vtx}| < 10 cm");
    hCollisionCounter->GetXaxis()->SetBinLabel(8, "sel8 && |Z_{vtx}| < 10 cm && global dimuon");
    hCollisionCounter->GetXaxis()->SetBinLabel(9, "sel8 && |Z_{vtx}| < 10 cm && global dimuon for QA");

    fRegistry.add("hNmu", "#mu multiplicity;N_{#mu} per collision", kTH1D, {{11, -0.5, 10.5}}, false);

    if (fillQAHistograms) {
      fRegistry.add("MFTMCHMID/hPt", "pT;p_{T} (GeV/c)", kTH1D, {{200, 0.0f, 10}}, false);
      fRegistry.add("MFTMCHMID/hEtaPhi", "#eta vs. #varphi;#varphi (rad.);#eta", kTH2D, {{180, 0, 2 * M_PI}, {80, -4.f, -2.f}}, false);
      fRegistry.add("MFTMCHMID/hEtaPhi_MatchedMCHMID", "#eta vs. #varphi;#varphi (rad.);#eta", kTH2D, {{180, 0, 2 * M_PI}, {80, -4.f, -2.f}}, false);
      fRegistry.add("MFTMCHMID/hDEtaDPhi", "#Delta#eta vs. #Delta#varphi;#Delta#varphi (rad.);#Delta#eta", kTH2D, {{90, -M_PI / 4, M_PI / 4}, {100, -0.5, +0.5}}, false);
      fRegistry.add("MFTMCHMID/hSign", "sign;sign", kTH1D, {{3, -1.5, +1.5}}, false);
      fRegistry.add("MFTMCHMID/hNclusters", "Nclusters;Nclusters", kTH1D, {{21, -0.5f, 20.5}}, false);
      fRegistry.add("MFTMCHMID/hNclustersMFT", "NclustersMFT;Nclusters MFT", kTH1D, {{11, -0.5f, 10.5}}, false);
      fRegistry.add("MFTMCHMID/hRatAbsorberEnd", "R at absorber end;R at absorber end (cm)", kTH1D, {{200, 0, 100}}, false);
      fRegistry.add("MFTMCHMID/hPDCA_Rabs", "pDCA vs. Rabs;R at absorber end (cm);p #times DCA (GeV/c #upoint cm)", kTH2D, {{200, 0, 100}, {100, 0, 1000}}, false);
      fRegistry.add("MFTMCHMID/hChi2", "chi2;#chi^{2}/ndf", kTH1D, {{100, 0.0f, 10}}, false);
      fRegistry.add("MFTMCHMID/hChi2MFT", "chi2 MFT;#chi^{2} MFT/ndf", kTH1D, {{100, 0.0f, 10}}, false);
      fRegistry.add("MFTMCHMID/hChi2MatchMCHMID", "chi2 match MCH-MID;#chi^{2}", kTH1D, {{200, 0.0f, 20}}, false);
      fRegistry.add("MFTMCHMID/hChi2MatchMCHMFT", "chi2 match MCH-MFT;#chi^{2}", kTH1D, {{100, 0.0f, 100}}, false);
      fRegistry.add("MFTMCHMID/hDCAxy2D", "DCA x vs. y;DCA_{x} (cm);DCA_{y} (cm)", kTH2D, {{400, -1, 1}, {400, -1, +1}}, false);
      fRegistry.add("MFTMCHMID/hDCAxy2DinSigma", "DCA x vs. y in sigma;DCA_{x} (#sigma);DCA_{y} (#sigma)", kTH2D, {{200, -10, 10}, {200, -10, +10}}, false);
      fRegistry.add("MFTMCHMID/hDCAxy", "DCAxy;DCA_{xy} (cm);", kTH1D, {{100, 0, 1}}, false);
      fRegistry.add("MFTMCHMID/hDCAxyz", "DCA xy vs. z;DCA_{xy} (cm);DCA_{z} (cm)", kTH2D, {{100, 0, 1}, {200, -0.1, 0.1}}, false);
      fRegistry.add("MFTMCHMID/hDCAxyinSigma", "DCAxy in sigma;DCA_{xy} (#sigma);", kTH1D, {{100, 0, 10}}, false);
      fRegistry.add("MFTMCHMID/hDCAxResolutionvsPt", "DCA_{x} resolution vs. p_{T};p_{T} (GeV/c);DCA_{x} resolution (#mum);", kTH2D, {{100, 0, 10.f}, {500, 0, 500}}, false);
      fRegistry.add("MFTMCHMID/hDCAyResolutionvsPt", "DCA_{y} resolution vs. p_{T};p_{T} (GeV/c);DCA_{y} resolution (#mum);", kTH2D, {{100, 0, 10.f}, {500, 0, 500}}, false);
      fRegistry.add("MFTMCHMID/hDCAxyResolutionvsPt", "DCA_{xy} resolution vs. p_{T};p_{T} (GeV/c);DCA_{xy} resolution (#mum);", kTH2D, {{100, 0, 10.f}, {500, 0, 500}}, false);

      // const AxisSpec axisMass{{0.20, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.30, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.40, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.50, 0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59, 0.60, 0.61, 0.62, 0.63, 0.64, 0.65, 0.66, 0.67, 0.68, 0.69, 0.70, 0.71, 0.72, 0.73, 0.74, 0.75, 0.76, 0.77, 0.78, 0.79, 0.80, 0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 0.87, 0.88, 0.89, 0.90, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 1.00, 1.01, 1.02, 1.03, 1.04, 1.05, 1.06, 1.07, 1.08, 1.09, 1.10, 1.11, 1.12, 1.13, 1.14, 1.15, 1.16, 1.17, 1.18, 1.19, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00, 2.05, 2.10, 2.15, 2.20, 2.25, 2.30, 2.35, 2.40, 2.45, 2.50, 2.55, 2.60, 2.65, 2.70, 2.75, 2.80, 2.85, 2.90, 2.95, 3.00, 3.05, 3.10, 3.15, 3.20, 3.25, 3.30, 3.35, 3.40, 3.45, 3.50, 3.55, 3.60, 3.65, 3.70, 3.75, 3.80, 3.85, 3.90, 3.95, 4.00}, "m_{#mu#mu} (GeV/c^{2})"};
      const AxisSpec axisPt{{0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10}, "p_{T,#mu#mu} (GeV/c)"};
      // const AxisSpec axisY{30, -4.0, -2.5, "#y_{#mu#mu}"};

      fRegistry.add("Pair/same/uls/hs", "dimuon;m_{#mu#mu} (GeV/c^{2});p_{T,#mu#mu} (GeV/c);y_{#mu#mu};", kTHnSparseD, {{380, 0.2, 4}, axisPt, {30, -4, -2.5}}, false);
      fRegistry.add("Pair/same/lspp/hs", "dimuon;m_{#mu#mu} (GeV/c^{2});p_{T,#mu#mu} (GeV/c);y_{#mu#mu};", kTHnSparseD, {{380, 0.2, 4}, axisPt, {30, -4, -2.5}}, false);
      fRegistry.add("Pair/same/lsmm/hs", "dimuon;m_{#mu#mu} (GeV/c^{2});p_{T,#mu#mu} (GeV/c);y_{#mu#mu};", kTHnSparseD, {{380, 0.2, 4}, axisPt, {30, -4, -2.5}}, false);
    }
  }

  std::vector<std::tuple<int, int, int>> vec_min_chi2MatchMCHMFT; // std::pair<globalIndex of global muon, globalIndex of matched MCH-MID, globalIndex of MFT> -> chi2MatchMCHMFT;
  template <typename TMuons>
  void findBestMatchPerMCHMID(TMuons const& muons)
  {
    vec_min_chi2MatchMCHMFT.reserve(muons.size());
    for (const auto& muon : muons) {
      if (muon.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack) {
        const auto& muons_per_MCHMID = muons.sliceBy(fwdtracksPerMCHTrack, muon.globalIndex());
        // LOGF(info, "stanadalone: muon.globalIndex() = %d, muon.chi2MatchMCHMFT() = %f", muon.globalIndex(), muon.chi2MatchMCHMFT());
        // LOGF(info, "muons_per_MCHMID.size() = %d", muons_per_MCHMID.size());

        float min_chi2MatchMCHMFT = 1e+10;
        std::tuple<int, int, int> tupleIds_at_min;
        for (const auto& muon_tmp : muons_per_MCHMID) {
          if (muon_tmp.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack) {
            // LOGF(info, "muon_tmp.globalIndex() = %d, muon_tmp.matchMCHTrackId() = %d, muon_tmp.matchMFTTrackId() = %d, muon_tmp.chi2MatchMCHMFT() = %f", muon_tmp.globalIndex(), muon_tmp.matchMCHTrackId(), muon_tmp.matchMFTTrackId(), muon_tmp.chi2MatchMCHMFT());
            if (0.f < muon_tmp.chi2MatchMCHMFT() && muon_tmp.chi2MatchMCHMFT() < min_chi2MatchMCHMFT) {
              min_chi2MatchMCHMFT = muon_tmp.chi2MatchMCHMFT();
              tupleIds_at_min = std::make_tuple(muon_tmp.globalIndex(), muon_tmp.matchMCHTrackId(), muon_tmp.matchMFTTrackId());
            }
          }
        }
        vec_min_chi2MatchMCHMFT.emplace_back(tupleIds_at_min);
        // LOGF(info, "min: muon_tmp.globalIndex() = %d, muon_tmp.matchMCHTrackId() = %d, muon_tmp.matchMFTTrackId() = %d, muon_tmp.chi2MatchMCHMFT() = %f", std::get<0>(tupleIds_at_min), std::get<1>(tupleIds_at_min), std::get<2>(tupleIds_at_min), min_chi2MatchMCHMFT);
      }
    } // end of muon loop
  }

  PresliceUnsorted<aod::FwdTracks> perMFTTrack = o2::aod::fwdtrack::matchMFTTrackId;
  template <typename TCollision, typename TFwdTrack, typename TFwdTracks, typename TMFTTracks, typename TMFTCovs>
  bool isBestMatch(TCollision const& collision, TFwdTrack const& fwdtrack, TFwdTracks const& fwdtracks, TMFTTracks const& mfttracks, TMFTCovs const& mftCovs)
  {
    if (fwdtrack.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack) {
      std::map<int, float> map_chi2MCHMFT;
      map_chi2MCHMFT[fwdtrack.globalIndex()] = fwdtrack.chi2MatchMCHMFT(); // add myself
      // LOGF(info, "add myself: fwdtrack.globalIndex() = %d, fwdtrack.chi2MatchMCHMFT() = %f", fwdtrack.globalIndex(), fwdtrack.chi2MatchMCHMFT());

      auto candidates = fwdtracks.sliceBy(perMFTTrack, fwdtrack.matchMFTTrackId()); // global muon candidates including this fwdtrack.
      for (const auto& candidate : candidates) {
        if (candidate.globalIndex() == fwdtrack.globalIndex()) {
          continue; // don't add fwdtrack.globalIndex() again.
        }

        float pt = 999.f, eta = 999.f, phi = 999.f;
        if (isSelectedGlobalMuon(collision, candidate, fwdtracks, mfttracks, mftCovs, pt, eta, phi)) {
          map_chi2MCHMFT[candidate.globalIndex()] = candidate.chi2MatchMCHMFT();
          // LOGF(info, "same MFT found: candidate.globalIndex() = %d, candidate.chi2MatchMCHMFT() = %f", candidate.globalIndex(), candidate.chi2MatchMCHMFT());
        }
      }

      auto it = std::min_element(map_chi2MCHMFT.begin(), map_chi2MCHMFT.end(), [](decltype(map_chi2MCHMFT)::value_type& l, decltype(map_chi2MCHMFT)::value_type& r) -> bool { return l.second < r.second; });

      // LOGF(info, "min: globalIndex = %d, chi2 = %f", it->first, it->second);
      // LOGF(info, "bool = %d", it->first == fwdtrack.globalIndex());

      if (it->first == fwdtrack.globalIndex()) { // search for minimum matching-chi2
        map_chi2MCHMFT.clear();
        return true;
      } else {
        map_chi2MCHMFT.clear();
        return false;
      }
    } else {
      return true;
    }
  }

  template <typename TCollision, typename TFwdTrack, typename TFwdTracks, typename TMFTTracks, typename TMFTCovs>
  bool isSelectedGlobalMuon(TCollision const& collision, TFwdTrack const& fwdtrack, TFwdTracks const&, TMFTTracks const&, TMFTCovs const& mftCovs, float& pt, float& eta, float& phi)
  {
    if (fwdtrack.trackType() != o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack) {
      return false;
    }

    if (std::find(vec_min_chi2MatchMCHMFT.begin(), vec_min_chi2MatchMCHMFT.end(), std::make_tuple(fwdtrack.globalIndex(), fwdtrack.matchMCHTrackId(), fwdtrack.matchMFTTrackId())) == vec_min_chi2MatchMCHMFT.end()) {
      return false;
    }

    if (fwdtrack.chi2MatchMCHMFT() > glMuonCutGroup.maxMatchingChi2MCHMFT) {
      return false;
    }

    if (fwdtrack.chi2MatchMCHMID() < 0.f) { // this should never happen. only for protection.
      return false;
    }

    if (fwdtrack.rAtAbsorberEnd() < glMuonCutGroup.minRabs || glMuonCutGroup.maxRabs < fwdtrack.rAtAbsorberEnd()) {
      return false;
    }

    auto mchtrack = fwdtrack.template matchMCHTrack_as<TFwdTracks>(); // MCH-MID
    auto mfttrack = fwdtrack.template matchMFTTrack_as<TMFTTracks>(); // MFTsa

    float rAtAbsorberEnd = fwdtrack.rAtAbsorberEnd(); // this works only for GlobalMuonTrack
    int nClustersMFT = mfttrack.nClusters();
    int ndf_mchmft = 2.f * (mchtrack.nClusters() + nClustersMFT) - 5.f;
    float chi2 = fwdtrack.chi2() / ndf_mchmft;
    if (glMuonCutGroup.maxChi2 < chi2) {
      return false;
    }

    int ndf_mft = 2.f * nClustersMFT - 5.f;
    float chi2mft = mfttrack.chi2() / ndf_mft;
    if (glMuonCutGroup.maxChi2MFT < chi2mft) {
      return false;
    }

    o2::dataformats::GlobalFwdTrack propmuonAtPV_Matched = propagateMuon(mchtrack, mchtrack, collision, propagationPoint::kToVertex, glMuonCutGroup.matchingZ, mBz);
    // float ptMatchedMCHMID = propmuonAtPV_Matched.getPt();
    float etaMatchedMCHMID = propmuonAtPV_Matched.getEta();
    float phiMatchedMCHMID = propmuonAtPV_Matched.getPhi();
    o2::math_utils::bringTo02Pi(phiMatchedMCHMID);

    auto mfttrackcov = mftCovs.rawIteratorAt(map_mfttrackcovs[mfttrack.globalIndex()]);
    o2::track::TrackParCovFwd mftsa = getTrackParCovFwd(mfttrack, mfttrackcov);                                                // values at innermost update
    o2::dataformats::GlobalFwdTrack globalMuonRefit = o2::aod::fwdtrackutils::refitGlobalMuonCov(propmuonAtPV_Matched, mftsa); // this is track at IU.
    auto globalMuon = o2::aod::fwdtrackutils::propagateTrackParCovFwd(globalMuonRefit, fwdtrack.trackType(), collision, propagationPoint::kToVertex, glMuonCutGroup.matchingZ, mBz);
    pt = globalMuon.getPt();
    eta = globalMuon.getEta();
    phi = globalMuon.getPhi();
    o2::math_utils::bringTo02Pi(phi);
    if (glMuonCutGroup.refitGlobalMuon) {
      pt = propmuonAtPV_Matched.getP() * std::sin(2.f * std::atan(std::exp(-eta)));
    }

    float deta = etaMatchedMCHMID - eta;
    float dphi = phiMatchedMCHMID - phi;
    o2::math_utils::bringToPMPi(dphi);

    if (std::sqrt(std::pow(deta / glMuonCutGroup.maxDEta, 2) + std::pow(dphi / glMuonCutGroup.maxDPhi, 2)) > 1.f) {
      return false;
    }

    float dcaX = globalMuon.getX() - collision.posX();
    float dcaY = globalMuon.getY() - collision.posY();
    float dcaZ = globalMuon.getZ() - collision.posZ();
    float dcaXY = std::sqrt(dcaX * dcaX + dcaY * dcaY);
    float cXX = globalMuon.getSigma2X();
    float cYY = globalMuon.getSigma2Y();
    float cXY = globalMuon.getSigmaXY();
    if (glMuonCutGroup.maxDCAxy < dcaXY) {
      return false;
    }

    float det = cXX * cYY - cXY * cXY; // determinanat
    float dcaXYinSigma = 999.f;
    if (det < 0) {
      dcaXYinSigma = 999.f;
    } else {
      dcaXYinSigma = std::sqrt(std::fabs((dcaX * dcaX * cYY + dcaY * dcaY * cXX - 2.f * dcaX * dcaY * cXY) / det / 2.f)); // dca xy in sigma
    }
    float sigma_dcaXY = dcaXY / dcaXYinSigma;

    o2::dataformats::GlobalFwdTrack propmuonAtDCA_Matched = propagateMuon(mchtrack, mchtrack, collision, propagationPoint::kToDCA, glMuonCutGroup.matchingZ, mBz);
    float dcaX_Matched = propmuonAtDCA_Matched.getX() - collision.posX();
    float dcaY_Matched = propmuonAtDCA_Matched.getY() - collision.posY();
    float dcaXY_Matched = std::sqrt(dcaX_Matched * dcaX_Matched + dcaY_Matched * dcaY_Matched);
    float pDCA = mchtrack.p() * dcaXY_Matched;

    if (pt < glMuonCutGroup.minPt || glMuonCutGroup.maxPt < pt) {
      return false;
    }

    if (eta < glMuonCutGroup.minEta || glMuonCutGroup.maxEta < eta) {
      return false;
    }

    if (fillQAHistograms) {
      fRegistry.fill(HIST("MFTMCHMID/hPt"), pt);
      fRegistry.fill(HIST("MFTMCHMID/hEtaPhi"), phi, eta);
      fRegistry.fill(HIST("MFTMCHMID/hEtaPhi_MatchedMCHMID"), phiMatchedMCHMID, etaMatchedMCHMID);
      fRegistry.fill(HIST("MFTMCHMID/hDEtaDPhi"), dphi, deta);
      fRegistry.fill(HIST("MFTMCHMID/hSign"), fwdtrack.sign());
      fRegistry.fill(HIST("MFTMCHMID/hNclusters"), fwdtrack.nClusters());
      fRegistry.fill(HIST("MFTMCHMID/hNclustersMFT"), nClustersMFT);
      fRegistry.fill(HIST("MFTMCHMID/hPDCA_Rabs"), rAtAbsorberEnd, pDCA);
      fRegistry.fill(HIST("MFTMCHMID/hRatAbsorberEnd"), rAtAbsorberEnd);
      fRegistry.fill(HIST("MFTMCHMID/hChi2"), chi2);
      fRegistry.fill(HIST("MFTMCHMID/hChi2MFT"), chi2mft);
      fRegistry.fill(HIST("MFTMCHMID/hChi2MatchMCHMID"), fwdtrack.chi2MatchMCHMID());
      fRegistry.fill(HIST("MFTMCHMID/hChi2MatchMCHMFT"), fwdtrack.chi2MatchMCHMFT());
      fRegistry.fill(HIST("MFTMCHMID/hDCAxy2D"), dcaX, dcaY);
      fRegistry.fill(HIST("MFTMCHMID/hDCAxy2DinSigma"), dcaX / std::sqrt(cXX), dcaY / std::sqrt(cYY));
      fRegistry.fill(HIST("MFTMCHMID/hDCAxy"), dcaXY);
      fRegistry.fill(HIST("MFTMCHMID/hDCAxyz"), dcaXY, dcaZ);
      fRegistry.fill(HIST("MFTMCHMID/hDCAxyinSigma"), dcaXYinSigma);
      fRegistry.fill(HIST("MFTMCHMID/hDCAxResolutionvsPt"), pt, std::sqrt(cXX) * 1e+4); // convert cm to um
      fRegistry.fill(HIST("MFTMCHMID/hDCAyResolutionvsPt"), pt, std::sqrt(cYY) * 1e+4); // convert cm to um
      fRegistry.fill(HIST("MFTMCHMID/hDCAxyResolutionvsPt"), pt, sigma_dcaXY * 1e+4);   // convert cm to um
    }

    return true;
  }

  template <typename TCollision, typename TMuons, typename TFwdTracks, typename TMFTTracks, typename TMFTCovs>
  void runPairing(TCollision const& collision, TMuons const& posMuons, TMuons const& negMuons, TFwdTracks const& fwdtracks, TMFTTracks const& mfttracks, TMFTCovs const& mftCovs)
  {
    for (const auto& pos : posMuons) {
      auto v1 = std::get<1>(pos);
      auto fwdtrack1 = fwdtracks.rawIteratorAt(std::get<0>(pos));
      if (!isBestMatch(collision, fwdtrack1, fwdtracks, mfttracks, mftCovs)) {
        continue;
      }

      for (const auto& neg : negMuons) {
        auto v2 = std::get<1>(neg);
        auto fwdtrack2 = fwdtracks.rawIteratorAt(std::get<0>(neg));
        if (!isBestMatch(collision, fwdtrack2, fwdtracks, mfttracks, mftCovs)) {
          continue;
        }

        auto v12 = v1 + v2;
        fRegistry.fill(HIST("Pair/same/uls/hs"), v12.M(), v12.Pt(), v12.Rapidity());

      } // end end of negative muon loop
    } // end end of positive muon loop

    for (int i1 = 0; i1 < static_cast<int>(posMuons.size()) - 1; i1++) {
      auto pos1 = posMuons[i1];
      auto v1 = std::get<1>(pos1);
      auto fwdtrack1 = fwdtracks.rawIteratorAt(std::get<0>(pos1));
      if (!isBestMatch(collision, fwdtrack1, fwdtracks, mfttracks, mftCovs)) {
        continue;
      }

      for (int i2 = i1 + 1; i2 < static_cast<int>(posMuons.size()); i2++) {
        auto pos2 = posMuons[i2];
        auto v2 = std::get<1>(pos2);
        auto fwdtrack2 = fwdtracks.rawIteratorAt(std::get<0>(pos2));
        // if (std::get<0>(pos1) == std::get<0>(pos2)) {
        //   continue;
        // }

        if (!isBestMatch(collision, fwdtrack2, fwdtracks, mfttracks, mftCovs)) {
          continue;
        }

        auto v12 = v1 + v2;
        fRegistry.fill(HIST("Pair/same/lspp/hs"), v12.M(), v12.Pt(), v12.Rapidity());
      } // end end of positive muon loop
    } // end end of positive muon loop

    for (int i1 = 0; i1 < static_cast<int>(negMuons.size()) - 1; i1++) {
      auto neg1 = negMuons[i1];
      auto v1 = std::get<1>(neg1);
      auto fwdtrack1 = fwdtracks.rawIteratorAt(std::get<0>(neg1));
      if (!isBestMatch(collision, fwdtrack1, fwdtracks, mfttracks, mftCovs)) {
        continue;
      }

      for (int i2 = i1 + 1; i2 < static_cast<int>(negMuons.size()); i2++) {
        auto neg2 = negMuons[i2];
        auto v2 = std::get<1>(neg2);
        auto fwdtrack2 = fwdtracks.rawIteratorAt(std::get<0>(neg2));
        // if (std::get<0>(neg1) == std::get<0>(neg2)) {
        //   continue;
        // }

        if (!isBestMatch(collision, fwdtrack2, fwdtracks, mfttracks, mftCovs)) {
          continue;
        }

        auto v12 = v1 + v2;
        fRegistry.fill(HIST("Pair/same/lsmm/hs"), v12.M(), v12.Pt(), v12.Rapidity());
      } // end end of negative muon loop
    } // end end of negative muon loop
  }

  using MyCollisions = soa::Join<aod::Collisions, aod::EvSels>;
  using MyCollision = MyCollisions::iterator;

  using MyBCs = soa::Join<aod::BCsWithTimestamps, aod::BcSels>;
  using MyBC = MyBCs::iterator;

  using MyFwdTracks = soa::Join<aod::FwdTracks, aod::FwdTracksCov>; // muon tracks are repeated. i.e. not exclusive.
  using MyFwdTrack = MyFwdTracks::iterator;

  SliceCache cache;
  Preslice<aod::FwdTracks> perCollision = o2::aod::fwdtrack::collisionId;
  Preslice<aod::FwdTrackAssoc> fwdtrackIndicesPerCollision = aod::track_association::collisionId;
  PresliceUnsorted<aod::FwdTrackAssoc> fwdtrackIndicesPerFwdTrack = aod::track_association::fwdtrackId;
  PresliceUnsorted<aod::FwdTracks> fwdtracksPerMCHTrack = aod::fwdtrack::matchMCHTrackId;

  std::unordered_map<int, int> map_mfttrackcovs;
  std::unordered_multimap<int, int> map_sa2gl; // sa muon index -> gl muon index
  std::unordered_map<int, bool> mapSelectedCollisions;

  void processSA(MyCollisions const& collisions, MyBCs const&, MyFwdTracks const& fwdtracks, aod::MFTTracks const& mfttracks, aod::MFTTracksCov const& mftCovs)
  {
    for (const auto& mfttrackConv : mftCovs) {
      map_mfttrackcovs[mfttrackConv.matchMFTTrackId()] = mfttrackConv.globalIndex();
    }
    findBestMatchPerMCHMID(fwdtracks);

    for (const auto& fwdtrack : fwdtracks) {
      // LOGF(info, "fwdtrack.globalIndex() = %d, fwdtrack.trackType() = %d, fwdtrack.matchMCHTrackId() = %d, fwdtrack.matchMFTTrackId() = %d", fwdtrack.globalIndex(), fwdtrack.trackType(), fwdtrack.matchMCHTrackId(), fwdtrack.matchMFTTrackId());
      map_sa2gl.insert({fwdtrack.matchMCHTrackId(), fwdtrack.globalIndex()});
    }

    for (const auto& collision : collisions) {
      auto bc = collision.template bc_as<MyBCs>();
      initCCDB(bc);
      fRegistry.fill(HIST("hCollisionCounter"), 1);

      if (collision.selection_bit(o2::aod::evsel::kIsTriggerTVX)) {
        fRegistry.fill(HIST("hCollisionCounter"), 2);
      }
      if (collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
        fRegistry.fill(HIST("hCollisionCounter"), 3);
      }
      if (collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
        fRegistry.fill(HIST("hCollisionCounter"), 4);
      }
      if (eventCutGroup.minZvtx < collision.posZ() && collision.posZ() < eventCutGroup.maxZvtx) {
        fRegistry.fill(HIST("hCollisionCounter"), 5);
      }
      if (collision.sel8()) {
        fRegistry.fill(HIST("hCollisionCounter"), 6);
      }
      if (collision.sel8() && eventCutGroup.minZvtx < collision.posZ() && collision.posZ() < eventCutGroup.maxZvtx) {
        fRegistry.fill(HIST("hCollisionCounter"), 7);
      }

      if (collision.sel8() && eventCutGroup.minZvtx < collision.posZ() && collision.posZ() < eventCutGroup.maxZvtx) {
        int nGlobalMuon = 0;
        // LOGF(info, "collision.globalIndex() = %d, fwdtrackIdsThisCollision = %d", collision.globalIndex(), fwdtrackIdsThisCollision.size());
        auto fwdtracks_per_coll = fwdtracks.sliceBy(perCollision, collision.globalIndex());
        for (const auto& fwdtrack : fwdtracks_per_coll) {
          float pt = 999.f, eta = 999.f, phi = 999.f;
          if (isSelectedGlobalMuon(collision, fwdtrack, fwdtracks, mfttracks, mftCovs, pt, eta, phi)) {
            nGlobalMuon++;
          }
        }

        // LOGF(info, "nGlobalMuon = %d", nGlobalMuon);
        fRegistry.fill(HIST("hNmu"), nGlobalMuon);

        if (nGlobalMuon >= minNmuon) {
          fRegistry.fill(HIST("hCollisionCounter"), 8);
          mapSelectedCollisions[collision.globalIndex()] = true;
          tags(true);
        } else {
          mapSelectedCollisions[collision.globalIndex()] = false;
          tags(false);
        }
      } else {
        mapSelectedCollisions[collision.globalIndex()] = false;
        tags(false);
      }
    } // end of collision loop

    // trigger QA
    for (const auto& collision : collisions) {
      auto bc = collision.template bc_as<MyBCs>();
      initCCDB(bc);

      if (!mapSelectedCollisions[collision.globalIndex()]) {
        continue;
      }
      fRegistry.fill(HIST("hCollisionCounter"), 9);

      auto fwdtracks_per_coll = fwdtracks.sliceBy(perCollision, collision.globalIndex());
      std::vector<std::tuple<int, ROOT::Math::PtEtaPhiMVector>> posMuons;
      std::vector<std::tuple<int, ROOT::Math::PtEtaPhiMVector>> negMuons;
      posMuons.reserve(fwdtracks_per_coll.size());
      negMuons.reserve(fwdtracks_per_coll.size());

      for (const auto& fwdtrack : fwdtracks_per_coll) {
        float pt = 999.f, eta = 999.f, phi = 999.f;
        if (isSelectedGlobalMuon(collision, fwdtrack, fwdtracks, mfttracks, mftCovs, pt, eta, phi)) {
          if (fwdtrack.sign() > 0) {
            posMuons.emplace_back(std::make_tuple(fwdtrack.globalIndex(), ROOT::Math::PtEtaPhiMVector(pt, eta, phi, o2::constants::physics::MassMuon)));
          } else {
            negMuons.emplace_back(std::make_tuple(fwdtrack.globalIndex(), ROOT::Math::PtEtaPhiMVector(pt, eta, phi, o2::constants::physics::MassMuon)));
          }
        }
      } // end of fwdtrack loop

      // make pairs
      runPairing(collision, posMuons, negMuons, fwdtracks, mfttracks, mftCovs);
    } // end of collision loop

    map_mfttrackcovs.clear();
    map_sa2gl.clear();
    mapSelectedCollisions.clear();
    vec_min_chi2MatchMCHMFT.clear();
    vec_min_chi2MatchMCHMFT.shrink_to_fit();
  }
  PROCESS_SWITCH(globalDimuonFilter, processSA, "process without FTTCA", false);

  void processFTTCA(MyCollisions const& collisions, MyBCs const&, MyFwdTracks const& fwdtracks, aod::MFTTracks const& mfttracks, aod::MFTTracksCov const& mftCovs, aod::FwdTrackAssoc const& fwdtrackIndices)
  {
    for (const auto& mfttrackConv : mftCovs) {
      map_mfttrackcovs[mfttrackConv.matchMFTTrackId()] = mfttrackConv.globalIndex();
    }
    findBestMatchPerMCHMID(fwdtracks);

    for (const auto& fwdtrack : fwdtracks) {
      // LOGF(info, "fwdtrack.globalIndex() = %d, fwdtrack.trackType() = %d, fwdtrack.matchMCHTrackId() = %d, fwdtrack.matchMFTTrackId() = %d", fwdtrack.globalIndex(), fwdtrack.trackType(), fwdtrack.matchMCHTrackId(), fwdtrack.matchMFTTrackId());
      map_sa2gl.insert({fwdtrack.matchMCHTrackId(), fwdtrack.globalIndex()});
    }

    for (const auto& collision : collisions) {
      auto bc = collision.template bc_as<MyBCs>();
      initCCDB(bc);

      fRegistry.fill(HIST("hCollisionCounter"), 1);

      if (collision.selection_bit(o2::aod::evsel::kIsTriggerTVX)) {
        fRegistry.fill(HIST("hCollisionCounter"), 2);
      }
      if (collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
        fRegistry.fill(HIST("hCollisionCounter"), 3);
      }
      if (collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
        fRegistry.fill(HIST("hCollisionCounter"), 4);
      }
      if (eventCutGroup.minZvtx < collision.posZ() && collision.posZ() < eventCutGroup.maxZvtx) {
        fRegistry.fill(HIST("hCollisionCounter"), 5);
      }
      if (collision.sel8()) {
        fRegistry.fill(HIST("hCollisionCounter"), 6);
      }

      if (collision.sel8() && eventCutGroup.minZvtx < collision.posZ() && collision.posZ() < eventCutGroup.maxZvtx) {
        fRegistry.fill(HIST("hCollisionCounter"), 7);

        int nGlobalMuon = 0;
        auto fwdtrackIdsThisCollision = fwdtrackIndices.sliceBy(fwdtrackIndicesPerCollision, collision.globalIndex());
        // LOGF(info, "collision.globalIndex() = %d, fwdtrackIdsThisCollision = %d", collision.globalIndex(), fwdtrackIdsThisCollision.size());
        for (const auto& fwdtrackId : fwdtrackIdsThisCollision) {
          auto fwdtrack = fwdtrackId.template fwdtrack_as<MyFwdTracks>();
          float pt = 999.f, eta = 999.f, phi = 999.f;
          if (isSelectedGlobalMuon(collision, fwdtrack, fwdtracks, mfttracks, mftCovs, pt, eta, phi)) {
            nGlobalMuon++;
          }
        }

        // LOGF(info, "nGlobalMuon = %d", nGlobalMuon);
        fRegistry.fill(HIST("hNmu"), nGlobalMuon);

        if (nGlobalMuon >= minNmuon) {
          fRegistry.fill(HIST("hCollisionCounter"), 8);
          mapSelectedCollisions[collision.globalIndex()] = true;
          tags(true);
        } else {
          mapSelectedCollisions[collision.globalIndex()] = false;
          tags(false);
        }
      } else {
        mapSelectedCollisions[collision.globalIndex()] = false;
        tags(false);
      }

    } // end of collision loop

    // trigger QA
    for (const auto& collision : collisions) {
      auto bc = collision.template bc_as<MyBCs>();
      initCCDB(bc);

      if (!mapSelectedCollisions[collision.globalIndex()]) {
        continue;
      }
      fRegistry.fill(HIST("hCollisionCounter"), 9);

      auto fwdtrackIdsThisCollision = fwdtrackIndices.sliceBy(fwdtrackIndicesPerCollision, collision.globalIndex());
      std::vector<std::tuple<int, ROOT::Math::PtEtaPhiMVector>> posMuons;
      std::vector<std::tuple<int, ROOT::Math::PtEtaPhiMVector>> negMuons;
      posMuons.reserve(fwdtrackIdsThisCollision.size());
      negMuons.reserve(fwdtrackIdsThisCollision.size());

      for (const auto& fwdtrackId : fwdtrackIdsThisCollision) {
        auto fwdtrack = fwdtrackId.template fwdtrack_as<MyFwdTracks>();
        float pt = 999.f, eta = 999.f, phi = 999.f;
        if (isSelectedGlobalMuon(collision, fwdtrack, fwdtracks, mfttracks, mftCovs, pt, eta, phi)) {
          if (fwdtrack.sign() > 0) {
            posMuons.emplace_back(std::make_tuple(fwdtrack.globalIndex(), ROOT::Math::PtEtaPhiMVector(pt, eta, phi, o2::constants::physics::MassMuon)));
          } else {
            negMuons.emplace_back(std::make_tuple(fwdtrack.globalIndex(), ROOT::Math::PtEtaPhiMVector(pt, eta, phi, o2::constants::physics::MassMuon)));
          }
        }
      } // end of fwdtrack loop

      // make pairs
      runPairing(collision, posMuons, negMuons, fwdtracks, mfttracks, mftCovs);

      posMuons.clear();
      posMuons.shrink_to_fit();
      negMuons.clear();
      negMuons.shrink_to_fit();
    } // end of collision loop

    map_mfttrackcovs.clear();
    map_sa2gl.clear();
    mapSelectedCollisions.clear();
    vec_min_chi2MatchMCHMFT.clear();
    vec_min_chi2MatchMCHMFT.shrink_to_fit();
  }
  PROCESS_SWITCH(globalDimuonFilter, processFTTCA, "process with FTTCA", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfg)
{
  return WorkflowSpec{
    adaptAnalysisTask<globalDimuonFilter>(cfg, TaskName{"em-global-dimuon-filter"})};
}
