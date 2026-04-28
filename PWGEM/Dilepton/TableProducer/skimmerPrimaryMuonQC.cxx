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

#include "PWGEM/Dilepton/DataModel/dileptonTables.h"

#include "Common/Core/TableHelper.h"
#include "Common/Core/fwdtrackUtilities.h"
#include "Common/DataModel/CollisionAssociationTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "CommonConstants/PhysicsConstants.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DetectorsBase/Propagator.h"
#include "Field/MagneticField.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/DataTypes.h"
#include "Framework/runDataProcessing.h"
#include "GlobalTracking/MatchGlobalFwd.h"
#include "MCHTracking/TrackExtrap.h"
#include "MCHTracking/TrackParam.h"
#include "ReconstructionDataFormats/TrackFwd.h"

#include "Math/SMatrix.h"
#include "Math/Vector4D.h"
#include "TGeoGlobalMagField.h"

#include <algorithm>
#include <map>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::soa;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;
using namespace o2::aod::fwdtrackutils;

struct skimmerPrimaryMuonQC {
  using MyCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::EMEvSels>;
  using MyCollisionsWithSWT = soa::Join<MyCollisions, aod::EMSWTriggerBitsTMP>;

  using MyFwdTracks = soa::Join<aod::FwdTracks, aod::FwdTracksCov>; // muon tracks are repeated. i.e. not exclusive.
  using MyFwdTrack = MyFwdTracks::iterator;

  using MyFwdTracksMC = soa::Join<MyFwdTracks, aod::McFwdTrackLabels>;
  using MyFwdTrackMC = MyFwdTracksMC::iterator;

  using MFTTracksMC = soa::Join<o2::aod::MFTTracks, aod::McMFTTrackLabels>;
  using MFTTrackMC = MFTTracksMC::iterator;

  Produces<aod::EMPrimaryMuons> emprimarymuons;
  Produces<aod::EMPrimaryMuonsCov> emprimarymuonscov;

  // Configurables
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
  Configurable<bool> fillQAHistograms{"fillQAHistograms", false, "flag to fill QA histograms"};

  // for z shift for propagation
  Configurable<bool> cfgApplyZShiftFromCCDB{"cfgApplyZShiftFromCCDB", false, "flag to apply z shift"};
  Configurable<std::string> cfgZShiftPath{"cfgZShiftPath", "Users/m/mcoquet/ZShift", "CCDB path for z shift to apply to forward tracks"};
  Configurable<float> cfgManualZShift{"cfgManualZShift", 0, "manual z-shift for propagation of global muon to PV"};
  Configurable<float> matchingZ{"matchingZ", -77.5, "z position where matching is performed"};
  Configurable<bool> refitGlobalMuon{"refitGlobalMuon", true, "flag to refit global muon"};

  struct : ConfigurableGroup { // tight cut
    std::string prefix = "tagMuonCut";
    Configurable<float> minPt{"minPt", 0.8, "min pt for muon"};
    Configurable<float> maxPt{"maxPt", 1e+10, "max pt for muon"};
    Configurable<float> minEta{"minEta", -3.6, "min. eta acceptance for MFT-MCH-MID"};
    Configurable<float> maxEta{"maxEta", -2.5, "max. eta acceptance for MFT-MCH-MID"};
    Configurable<float> minRabs{"minRabs", 27.6, "min. R at absorber end for global muon (min. eta = -3.6)"}; // std::tan(2.f * std::atan(std::exp(- -3.6)) ) * -505. = 27.6
    Configurable<float> maxRabs{"maxRabs", 89.5, "max. R at absorber end"};
    Configurable<float> maxDCAxy{"maxDCAxy", 0.06, "max. DCAxy for global muons"};
    Configurable<float> maxMatchingChi2MCHMFT{"maxMatchingChi2MCHMFT", 40.f, "max. chi2 for MCH-MFT matching"};
    Configurable<float> maxChi2{"maxChi2", 4.f, "max. chi2/ndf for global muon"};
    // Configurable<int> minNclsMFT{"minNclsMFT", 5, "min ncluster of MFT"};
    // Configurable<int> minNclsMCH{"minNclsMCH", 5, "min ncluster of MCH"};
    Configurable<float> maxDEta{"maxDEta", 0.08, "max. deta between MFT-MCH-MID and MCH-MID"};
    Configurable<float> maxDPhi{"maxDPhi", 0.08, "max. dphi between MFT-MCH-MID and MCH-MID"};
  } tagMuonCut;

  struct : ConfigurableGroup { // loose cut
    std::string prefix = "probeMuonCut";
    Configurable<float> minPt{"minPt", 0.01, "min pt for muon"};
    Configurable<float> maxPt{"maxPt", 1e+10, "max pt for muon"};
    Configurable<float> minEtaSA{"minEtaSA", -4.0, "min. eta acceptance for MFT-MCH"};
    Configurable<float> maxEtaSA{"maxEtaSA", -2.5, "max. eta acceptance for MFT-MCH"};
    Configurable<float> minEtaGL{"minEtaGL", -3.6, "min. eta acceptance for MFT-MCH-MID"};
    Configurable<float> maxEtaGL{"maxEtaGL", -2.5, "max. eta acceptance for MFT-MCH-MID"};
    Configurable<float> minRabs{"minRabs", 17.6, "min. R at absorber end for global muon (min. eta = -3.6)"}; // std::tan(2.f * std::atan(std::exp(- -3.6)) ) * -505. = 27.6
    Configurable<float> midRabs{"midRabs", 26.5, "middle R at absorber end for pDCA cut"};
    Configurable<float> maxRabs{"maxRabs", 89.5, "max. R at absorber end"};
    Configurable<float> maxDCAxy{"maxDCAxy", 1.f, "max. DCAxy for global muons"};
    Configurable<float> maxPDCAforLargeR{"maxPDCAforLargeR", 324.f, "max. pDCA for large R at absorber end"};
    Configurable<float> maxPDCAforSmallR{"maxPDCAforSmallR", 594.f, "max. pDCA for small R at absorber end"};
    Configurable<float> maxMatchingChi2MCHMFT{"maxMatchingChi2MCHMFT", 100, "max. chi2 for MCH-MFT matching"};
    Configurable<float> maxChi2{"maxChi2", 1e+10, "max. chi2/ndf for global muon"};
    // Configurable<int> minNclsMFT{"minNclsMFT", 5, "min ncluster of MFT"};
    // Configurable<int> minNclsMCH{"minNclsMCH", 5, "min ncluster of MCH"};
    Configurable<float> maxDEta{"maxDEta", 1e+10, "max. deta between MFT-MCH-MID and MCH-MID"};
    Configurable<float> maxDPhi{"maxDPhi", 1e+10, "max. dphi between MFT-MCH-MID and MCH-MID"};
  } probeMuonCut;

  struct : ConfigurableGroup {
    std::string prefix = "pairCuts";
    Configurable<float> minMass{"minMass", 0.21, "min mass"};
    Configurable<float> maxMass{"maxMass", 0.30, "max mass"};
  } pairCuts;

  o2::ccdb::CcdbApi ccdbApi;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  int mRunNumber = 0;
  float mBz = 0;
  float mZShift = 0;

  HistogramRegistry fRegistry{"output", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};
  // static constexpr std::string_view muon_types[5] = {"MFTMCHMID/", "MFTMCHMIDOtherMatch/", "MFTMCH/", "MCHMID/", "MCH/"};

  void init(InitContext&)
  {
    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);
    ccdbApi.init(ccdburl);

    if (fillQAHistograms) {
      addHistograms();
    }
    mRunNumber = 0;
    mBz = 0;
    mZShift = 0;
  }

  void initCCDB(aod::BCsWithTimestamps::iterator const& bc)
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

    if (cfgApplyZShiftFromCCDB) {
      auto* zShift = ccdb->getForTimeStamp<std::vector<float>>(cfgZShiftPath, bc.timestamp());
      if (zShift != nullptr && !zShift->empty()) {
        LOGF(info, "reading z shift %f from %s", (*zShift)[0], cfgZShiftPath.value);
        mZShift = (*zShift)[0];
      } else {
        LOGF(info, "z shift is not found in ccdb path %s. set to 0 cm", cfgZShiftPath.value);
        mZShift = 0;
      }
    } else {
      LOGF(info, "z shift is manually set to %f cm", cfgManualZShift.value);
      mZShift = cfgManualZShift;
    }
  }

  void addHistograms()
  {
    // auto hMuonType = fRegistry.add<TH1>("hMuonType", "muon type", kTH1F, {{5, -0.5f, 4.5f}}, false);
    // hMuonType->GetXaxis()->SetBinLabel(1, "MFT-MCH-MID (global muon)");
    // hMuonType->GetXaxis()->SetBinLabel(2, "MFT-MCH-MID (global muon other match)");
    // hMuonType->GetXaxis()->SetBinLabel(3, "MFT-MCH");
    // hMuonType->GetXaxis()->SetBinLabel(4, "MCH-MID");
    // hMuonType->GetXaxis()->SetBinLabel(5, "MCH standalone");

    // fRegistry.add("MFTMCHMID/hPt", "pT;p_{T} (GeV/c)", kTH1F, {{200, 0.0f, 10}}, false);
    // fRegistry.add("MFTMCHMID/hEtaPhi", "#eta vs. #varphi;#varphi (rad.);#eta", kTH2F, {{180, 0, 2 * M_PI}, {80, -4.f, -2.f}}, false);
    // fRegistry.add("MFTMCHMID/hEtaPhi_MatchedMCHMID", "#eta vs. #varphi;#varphi (rad.);#eta", kTH2F, {{180, 0, 2 * M_PI}, {80, -4.f, -2.f}}, false);
    // fRegistry.add("MFTMCHMID/hDeltaPt_Pt", "#Deltap_{T}/p_{T} vs. p_{T};p_{T}^{gl} (GeV/c);(p_{T}^{sa} - p_{T}^{gl})/p_{T}^{gl}", kTH2F, {{100, 0, 10}, {200, -0.5, +0.5}}, false);
    // fRegistry.add("MFTMCHMID/hDeltaEta_Pt", "#Delta#eta vs. p_{T};p_{T}^{gl} (GeV/c);#Delta#eta", kTH2F, {{100, 0, 10}, {200, -0.5, +0.5}}, false);
    // fRegistry.add("MFTMCHMID/hDeltaPhi_Pt", "#Delta#varphi vs. p_{T};p_{T}^{gl} (GeV/c);#Delta#varphi (rad.)", kTH2F, {{100, 0, 10}, {200, -0.5, +0.5}}, false);
    // fRegistry.add("MFTMCHMID/hSign", "sign;sign", kTH1F, {{3, -1.5, +1.5}}, false);
    // fRegistry.add("MFTMCHMID/hNclusters", "Nclusters;Nclusters", kTH1F, {{21, -0.5f, 20.5}}, false);
    // fRegistry.add("MFTMCHMID/hNclustersMFT", "NclustersMFT;Nclusters MFT", kTH1F, {{11, -0.5f, 10.5}}, false);
    // fRegistry.add("MFTMCHMID/hRatAbsorberEnd", "R at absorber end;R at absorber end (cm)", kTH1F, {{100, 0.0f, 100}}, false);
    // fRegistry.add("MFTMCHMID/hPDCA_Rabs", "pDCA vs. Rabs;R at absorber end (cm);p #times DCA (GeV/c #upoint cm)", kTH2F, {{100, 0, 100}, {100, 0.0f, 1000}}, false);
    // fRegistry.add("MFTMCHMID/hChi2", "chi2;chi2/ndf", kTH1F, {{200, 0.0f, 20}}, false);
    // fRegistry.add("MFTMCHMID/hChi2MFT", "chi2 MFT;chi2 MFT/ndf", kTH1F, {{200, 0.0f, 20}}, false);
    // fRegistry.add("MFTMCHMID/hChi2MatchMCHMID", "chi2 match MCH-MID;chi2", kTH1F, {{200, 0.0f, 20}}, false);
    // fRegistry.add("MFTMCHMID/hChi2MatchMCHMFT", "chi2 match MCH-MFT;chi2", kTH1F, {{200, 0.0f, 100}}, false);
    // fRegistry.add("MFTMCHMID/hDCAxy2D", "DCA x vs. y;DCA_{x} (cm);DCA_{y} (cm)", kTH2F, {{200, -1, 1}, {200, -1, +1}}, false);
    // fRegistry.add("MFTMCHMID/hDCAxy2DinSigma", "DCA x vs. y in sigma;DCA_{x} (#sigma);DCA_{y} (#sigma)", kTH2F, {{200, -10, 10}, {200, -10, +10}}, false);
    // fRegistry.add("MFTMCHMID/hDCAxy", "DCAxy;DCA_{xy} (cm);", kTH1F, {{100, 0, 1}}, false);
    // fRegistry.add("MFTMCHMID/hDCAxyz", "DCA xy vs. z;DCA_{xy} (cm);DCA_{z} (cm)", kTH2F, {{100, 0, 1}, {200, -0.1, 0.1}}, false);
    // fRegistry.add("MFTMCHMID/hDCAxyinSigma", "DCAxy in sigma;DCA_{xy} (#sigma);", kTH1F, {{100, 0, 10}}, false);
    // fRegistry.add("MFTMCHMID/hDCAx_PosZ", "DCAx vs. posZ;Z_{vtx} (cm);DCA_{x} (cm)", kTH2F, {{200, -10, +10}, {400, -0.2, +0.2}}, false);
    // fRegistry.add("MFTMCHMID/hDCAy_PosZ", "DCAy vs. posZ;Z_{vtx} (cm);DCA_{y} (cm)", kTH2F, {{200, -10, +10}, {400, -0.2, +0.2}}, false);
    // fRegistry.add("MFTMCHMID/hDCAx_Phi", "DCAx vs. #varphi;#varphi (rad.);DCA_{x} (cm)", kTH2F, {{90, 0, 2 * M_PI}, {400, -0.2, +0.2}}, false);
    // fRegistry.add("MFTMCHMID/hDCAy_Phi", "DCAy vs. #varphi;#varphi (rad.);DCA_{y} (cm)", kTH2F, {{90, 0, 2 * M_PI}, {400, -0.2, +0.2}}, false);
    // fRegistry.add("MFTMCHMID/hNmu", "#mu multiplicity;N_{#mu} per collision", kTH1F, {{21, -0.5, 20.5}}, false);

    // fRegistry.addClone("MFTMCHMID/", "MCHMID/");
    // fRegistry.add("MFTMCHMID/hDCAxResolutionvsPt", "DCA_{x} vs. p_{T};p_{T} (GeV/c);DCA_{x} resolution (#mum);", kTH2F, {{100, 0, 10.f}, {500, 0, 500}}, false);
    // fRegistry.add("MFTMCHMID/hDCAyResolutionvsPt", "DCA_{y} vs. p_{T};p_{T} (GeV/c);DCA_{y} resolution (#mum);", kTH2F, {{100, 0, 10.f}, {500, 0, 500}}, false);
    // fRegistry.add("MFTMCHMID/hDCAxyResolutionvsPt", "DCA_{xy} vs. p_{T};p_{T} (GeV/c);DCA_{y} resolution (#mum);", kTH2F, {{100, 0, 10.f}, {500, 0, 500}}, false);
    // fRegistry.add("MCHMID/hDCAxResolutionvsPt", "DCA_{x} vs. p_{T};p_{T} (GeV/c);DCA_{x} resolution (#mum);", kTH2F, {{100, 0, 10.f}, {500, 0, 5e+5}}, false);
    // fRegistry.add("MCHMID/hDCAyResolutionvsPt", "DCA_{y} vs. p_{T};p_{T} (GeV/c);DCA_{y} resolution (#mum);", kTH2F, {{100, 0, 10.f}, {500, 0, 5e+5}}, false);
    // fRegistry.add("MCHMID/hDCAxyResolutionvsPt", "DCA_{xy} vs. p_{T};p_{T} (GeV/c);DCA_{y} resolution (#mum);", kTH2F, {{100, 0, 10.f}, {500, 0, 5e+5}}, false);

    fRegistry.add("Pair/uls/gl_gl/hMvsPt", "dimuon;m_{#mu#mu} (GeV/c^{2});p_{T,#mu} (GeV/c);", kTH2F, {{380, 0.2, 4.f}, {100, 0, 10}}, false);
    fRegistry.add("Pair/uls/gl_sa/hMvsPt", "dimuon;m_{#mu#mu} (GeV/c^{2});p_{T,#mu} (GeV/c);", kTH2F, {{380, 0.2, 4.f}, {100, 0, 10}}, false);
  }

  struct Muon {
    int globalIndex = -1;
    int collisionId = -1;
    int matchMCHTrackId = -1;
    int matchMFTTrackId = -1;
    uint8_t trackType = 99;
    int8_t sign = 0;
    float pt = 0;
    float eta = 0;
    float phi = 0;
    float dcaX = 0;  // in cm
    float dcaY = 0;  // in cm
    float dcaXY = 0; // in cm
    float cXX = 0;
    float cYY = 0;
    float cXY = 0;
    float rAtAbsorberEnd = 0;
    float pDCA = 0;
    float chi2ndf = 0;
    float chi2MatchMCHMID = 0;

    // only for global muons
    float ptMatchedMCHMID = 0;
    float etaMatchedMCHMID = 0;
    float phiMatchedMCHMID = 0;
    float chi2MatchMCHMFT = 0;
    float chi2mft = -999.f;
    uint64_t mftClusterSizesAndTrackFlags = 0;
  };

  bool isSelected(Muon const& muon)
  {
    if (muon.pt < probeMuonCut.minPt || probeMuonCut.maxPt < muon.pt) {
      return false;
    }

    if (muon.rAtAbsorberEnd < probeMuonCut.minRabs || probeMuonCut.maxRabs < muon.rAtAbsorberEnd) {
      return false;
    }

    if (muon.trackType == static_cast<uint8_t>(o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack)) {
      if (muon.eta < probeMuonCut.minEtaGL || probeMuonCut.maxEtaGL < muon.eta) {
        return false;
      }
      if (probeMuonCut.maxDCAxy < muon.dcaXY) {
        return false;
      }
      if (probeMuonCut.maxMatchingChi2MCHMFT < muon.chi2MatchMCHMFT) {
        return false;
      }
    } else if (muon.trackType == static_cast<uint8_t>(o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack)) {
      if (muon.eta < probeMuonCut.minEtaSA || probeMuonCut.maxEtaSA < muon.eta) {
        return false;
      }
      if (muon.rAtAbsorberEnd < probeMuonCut.midRabs ? muon.pDCA > probeMuonCut.maxPDCAforSmallR : muon.pDCA > probeMuonCut.maxPDCAforLargeR) {
        return false;
      }
    } else {
      return false;
    }

    return true;
  }

  bool isSelectedTight(Muon const& muon)
  {
    if (muon.trackType != static_cast<uint8_t>(o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack)) { // tag muon should be tight.
      return false;
    }

    if (muon.pt < tagMuonCut.minPt || tagMuonCut.maxPt < muon.pt) {
      return false;
    }

    if (muon.rAtAbsorberEnd < tagMuonCut.minRabs || tagMuonCut.maxRabs < muon.rAtAbsorberEnd) {
      return false;
    }

    if (muon.chi2ndf < 0.f || tagMuonCut.maxChi2 < muon.chi2ndf) {
      return false;
    }

    if (muon.eta < tagMuonCut.minEta || tagMuonCut.maxEta < muon.eta) {
      return false;
    }

    if (tagMuonCut.maxDCAxy < muon.dcaXY) {
      return false;
    }

    if (tagMuonCut.maxMatchingChi2MCHMFT < muon.chi2MatchMCHMFT) {
      return false;
    }

    float deta = muon.etaMatchedMCHMID - muon.eta;
    float dphi = muon.phiMatchedMCHMID - muon.phi;
    // LOGF(info, "muon.trackType = %d, deta = %f, dphi = %f", muon.trackType, deta, dphi);
    if (std::sqrt(std::pow(deta / tagMuonCut.maxDEta, 2) + std::pow(dphi / tagMuonCut.maxDPhi, 2)) > 1.f) {
      return false;
    }

    return true;
  }

  template <bool isMC, typename TFwdTracks, typename TMFTTracks, typename TCollision, typename TFwdTrack>
  bool fillMuonInfo(TCollision const& collision, TFwdTrack fwdtrack)
  {
    if (fwdtrack.trackType() != static_cast<uint8_t>(o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack) && fwdtrack.trackType() != static_cast<uint8_t>(o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack)) {
      return false;
    }

    if (fwdtrack.chi2MatchMCHMID() < 0.f) { // this should never happen. only for protection.
      return false;
    }

    if (fwdtrack.chi2() < 0.f) { // this should never happen. only for protection.
      return false;
    }

    o2::dataformats::GlobalFwdTrack propmuonAtPV = propagateMuon(fwdtrack, fwdtrack, collision, propagationPoint::kToVertex, matchingZ, mBz, mZShift);
    float pt = propmuonAtPV.getPt();
    float eta = propmuonAtPV.getEta();
    float phi = propmuonAtPV.getPhi();
    o2::math_utils::bringTo02Pi(phi);

    float dcaX = propmuonAtPV.getX() - collision.posX();
    float dcaY = propmuonAtPV.getY() - collision.posY();
    float dcaXY = std::sqrt(dcaX * dcaX + dcaY * dcaY);
    float rAtAbsorberEnd = fwdtrack.rAtAbsorberEnd(); // this works only for GlobalMuonTrack
    float cXX = propmuonAtPV.getSigma2X();
    float cYY = propmuonAtPV.getSigma2Y();
    float cXY = propmuonAtPV.getSigmaXY();

    float pDCA = propmuonAtPV.getP() * dcaXY;
    int nClustersMFT = 0;
    float ptMatchedMCHMID = propmuonAtPV.getPt();
    float etaMatchedMCHMID = propmuonAtPV.getEta();
    float phiMatchedMCHMID = propmuonAtPV.getPhi();
    o2::math_utils::bringTo02Pi(phiMatchedMCHMID);
    float chi2mft = -999.f;
    uint64_t mftClusterSizesAndTrackFlags = 0;
    int ndf_mchmft = 1;

    if (fwdtrack.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack) {
      if (fwdtrack.chi2MatchMCHMFT() < 0.f) {
        return false;
      } // Users have to decide the best match between MFT and MCH-MID at analysis level. The same global muon is repeatedly stored.

      auto mchtrack = fwdtrack.template matchMCHTrack_as<TFwdTracks>(); // MCH-MID
      auto mfttrack = fwdtrack.template matchMFTTrack_as<TMFTTracks>(); // MFTsa
      if (mfttrack.chi2() < 0.f) {
        return false;
      }

      if constexpr (isMC) {
        if (!mfttrack.has_mcParticle() || !mchtrack.has_mcParticle() || !fwdtrack.has_mcParticle()) {
          return false;
        }
        // auto mcParticle_MFTMCHMID = fwdtrack.template mcParticle_as<aod::McParticles>(); // this is identical to mcParticle_MCHMID
        auto mcParticle_MCHMID = mchtrack.template mcParticle_as<aod::McParticles>(); // this is identical to mcParticle_MFTMCHMID
        auto mcParticle_MFT = mfttrack.template mcParticle_as<aod::McParticles>();
      }

      nClustersMFT = mfttrack.nClusters();
      mftClusterSizesAndTrackFlags = mfttrack.mftClusterSizesAndTrackFlags();
      ndf_mchmft = 2.f * (mchtrack.nClusters() + nClustersMFT) - 5.f;
      chi2mft = mfttrack.chi2();

      o2::dataformats::GlobalFwdTrack propmuonAtPV_Matched = propagateMuon(mchtrack, mchtrack, collision, propagationPoint::kToVertex, matchingZ, mBz, mZShift);
      ptMatchedMCHMID = propmuonAtPV_Matched.getPt();
      etaMatchedMCHMID = propmuonAtPV_Matched.getEta();
      phiMatchedMCHMID = propmuonAtPV_Matched.getPhi();
      o2::math_utils::bringTo02Pi(phiMatchedMCHMID);

      o2::dataformats::GlobalFwdTrack propmuonAtDCA_Matched = propagateMuon(mchtrack, mchtrack, collision, propagationPoint::kToDCA, matchingZ, mBz, mZShift);
      float dcaX_Matched = propmuonAtDCA_Matched.getX() - collision.posX();
      float dcaY_Matched = propmuonAtDCA_Matched.getY() - collision.posY();
      float dcaXY_Matched = std::sqrt(dcaX_Matched * dcaX_Matched + dcaY_Matched * dcaY_Matched);
      pDCA = mchtrack.p() * dcaXY_Matched;

      if (refitGlobalMuon) {
        pt = propmuonAtPV_Matched.getP() * std::sin(2.f * std::atan(std::exp(-eta)));
      }
    } else if (fwdtrack.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack) {
      o2::dataformats::GlobalFwdTrack propmuonAtRabs = propagateMuon(fwdtrack, fwdtrack, collision, propagationPoint::kToRabs, matchingZ, mBz, mZShift); // this is necessary only for MuonStandaloneTrack
      float xAbs = propmuonAtRabs.getX();
      float yAbs = propmuonAtRabs.getY();
      rAtAbsorberEnd = std::sqrt(xAbs * xAbs + yAbs * yAbs); // Redo propagation only for muon tracks // propagation of MFT tracks alredy done in reconstruction
      ndf_mchmft = 1;                                        // chi2 is already normalized by ndf for MCH-MID tracks.

      o2::dataformats::GlobalFwdTrack propmuonAtDCA = propagateMuon(fwdtrack, fwdtrack, collision, propagationPoint::kToDCA, matchingZ, mBz, mZShift);
      cXX = propmuonAtDCA.getSigma2X();
      cYY = propmuonAtDCA.getSigma2Y();
      cXY = propmuonAtDCA.getSigmaXY();
      dcaX = propmuonAtDCA.getX() - collision.posX();
      dcaY = propmuonAtDCA.getY() - collision.posY();
      dcaXY = std::sqrt(dcaX * dcaX + dcaY * dcaY);
      pDCA = fwdtrack.p() * dcaXY;
    } else {
      return false;
    }

    Muon muon;
    muon.globalIndex = fwdtrack.globalIndex();
    muon.collisionId = collision.globalIndex();
    muon.trackType = fwdtrack.trackType();
    muon.sign = fwdtrack.sign();
    muon.pt = pt;
    muon.eta = eta;
    muon.phi = phi;
    muon.dcaX = dcaX;
    muon.dcaY = dcaY;
    muon.dcaXY = dcaXY;
    muon.cXX = cXX;
    muon.cYY = cYY;
    muon.cXY = cXY;
    muon.rAtAbsorberEnd = rAtAbsorberEnd;
    muon.pDCA = pDCA;
    muon.chi2ndf = fwdtrack.chi2() / ndf_mchmft;
    muon.ptMatchedMCHMID = ptMatchedMCHMID;
    muon.etaMatchedMCHMID = etaMatchedMCHMID;
    muon.phiMatchedMCHMID = phiMatchedMCHMID;
    muon.chi2mft = chi2mft;
    muon.matchMCHTrackId = fwdtrack.matchMCHTrackId();
    muon.matchMFTTrackId = fwdtrack.matchMFTTrackId();
    muon.mftClusterSizesAndTrackFlags = mftClusterSizesAndTrackFlags;

    vecMuons.emplace_back(muon);
    return true;
  }

  template <typename TCollision, typename TMuon, typename TFwdTrack>
  bool fillFwdTrackTable(TCollision const& collision, TMuon const& muon, TFwdTrack const& fwdtrack)
  {
    emprimarymuons(collision.globalIndex(), fwdtrack.globalIndex(), fwdtrack.matchMFTTrackId(), fwdtrack.matchMCHTrackId(), fwdtrack.trackType(),
                   muon.pt, muon.eta, muon.phi, fwdtrack.sign(), muon.dcaX, muon.dcaY, muon.cXX, muon.cYY, muon.cXY, muon.ptMatchedMCHMID, muon.etaMatchedMCHMID, muon.phiMatchedMCHMID,
                   fwdtrack.nClusters(), muon.pDCA, muon.rAtAbsorberEnd, fwdtrack.chi2(), fwdtrack.chi2MatchMCHMID(), fwdtrack.chi2MatchMCHMFT(),
                   fwdtrack.mchBitMap(), fwdtrack.midBitMap(), fwdtrack.midBoards(), muon.mftClusterSizesAndTrackFlags, muon.chi2mft, true, false);

    // if (fillQAHistograms) {
    //   fRegistry.fill(HIST("hMuonType"), fwdtrack.trackType());
    //   if (fwdtrack.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack) {
    //     fRegistry.fill(HIST("MFTMCHMID/hPt"), pt);
    //     fRegistry.fill(HIST("MFTMCHMID/hEtaPhi"), phi, eta);
    //     fRegistry.fill(HIST("MFTMCHMID/hEtaPhi_MatchedMCHMID"), phiMatchedMCHMID, etaMatchedMCHMID);
    //     fRegistry.fill(HIST("MFTMCHMID/hDeltaPt_Pt"), pt, dpt);
    //     fRegistry.fill(HIST("MFTMCHMID/hDeltaEta_Pt"), pt, deta);
    //     fRegistry.fill(HIST("MFTMCHMID/hDeltaPhi_Pt"), pt, dphi);
    //     fRegistry.fill(HIST("MFTMCHMID/hSign"), fwdtrack.sign());
    //     fRegistry.fill(HIST("MFTMCHMID/hNclusters"), fwdtrack.nClusters());
    //     fRegistry.fill(HIST("MFTMCHMID/hNclustersMFT"), nClustersMFT);
    //     fRegistry.fill(HIST("MFTMCHMID/hPDCA_Rabs"), rAtAbsorberEnd, pDCA);
    //     fRegistry.fill(HIST("MFTMCHMID/hRatAbsorberEnd"), rAtAbsorberEnd);
    //     fRegistry.fill(HIST("MFTMCHMID/hChi2"), fwdtrack.chi2() / ndf_mchmft);
    //     fRegistry.fill(HIST("MFTMCHMID/hChi2MFT"), chi2mft / ndf_mft);
    //     fRegistry.fill(HIST("MFTMCHMID/hChi2MatchMCHMID"), fwdtrack.chi2MatchMCHMID());
    //     fRegistry.fill(HIST("MFTMCHMID/hChi2MatchMCHMFT"), fwdtrack.chi2MatchMCHMFT());
    //     fRegistry.fill(HIST("MFTMCHMID/hDCAxy2D"), dcaX, dcaY);
    //     fRegistry.fill(HIST("MFTMCHMID/hDCAxy2DinSigma"), dcaX / std::sqrt(cXX), dcaY / std::sqrt(cYY));
    //     fRegistry.fill(HIST("MFTMCHMID/hDCAxy"), dcaXY);
    //     fRegistry.fill(HIST("MFTMCHMID/hDCAxyz"), dcaXY, dcaZ);
    //     fRegistry.fill(HIST("MFTMCHMID/hDCAxyinSigma"), dcaXYinSigma);
    //     fRegistry.fill(HIST("MFTMCHMID/hDCAxResolutionvsPt"), pt, std::sqrt(cXX) * 1e+4); // convert cm to um
    //     fRegistry.fill(HIST("MFTMCHMID/hDCAyResolutionvsPt"), pt, std::sqrt(cYY) * 1e+4); // convert cm to um
    //     fRegistry.fill(HIST("MFTMCHMID/hDCAxyResolutionvsPt"), pt, sigma_dcaXY * 1e+4);   // convert cm to um
    //     fRegistry.fill(HIST("MFTMCHMID/hDCAx_PosZ"), collision.posZ(), dcaX);
    //     fRegistry.fill(HIST("MFTMCHMID/hDCAy_PosZ"), collision.posZ(), dcaY);
    //     fRegistry.fill(HIST("MFTMCHMID/hDCAx_Phi"), phi, dcaX);
    //     fRegistry.fill(HIST("MFTMCHMID/hDCAy_Phi"), phi, dcaY);
    //   } else if (fwdtrack.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack) {
    //     fRegistry.fill(HIST("MCHMID/hPt"), pt);
    //     fRegistry.fill(HIST("MCHMID/hEtaPhi"), phi, eta);
    //     fRegistry.fill(HIST("MCHMID/hEtaPhi_MatchedMCHMID"), phiMatchedMCHMID, etaMatchedMCHMID);
    //     fRegistry.fill(HIST("MCHMID/hDeltaPt_Pt"), pt, dpt);
    //     fRegistry.fill(HIST("MCHMID/hDeltaEta_Pt"), pt, deta);
    //     fRegistry.fill(HIST("MCHMID/hDeltaPhi_Pt"), pt, dphi);
    //     fRegistry.fill(HIST("MCHMID/hSign"), fwdtrack.sign());
    //     fRegistry.fill(HIST("MCHMID/hNclusters"), fwdtrack.nClusters());
    //     fRegistry.fill(HIST("MCHMID/hNclustersMFT"), nClustersMFT);
    //     fRegistry.fill(HIST("MCHMID/hPDCA_Rabs"), rAtAbsorberEnd, pDCA);
    //     fRegistry.fill(HIST("MCHMID/hRatAbsorberEnd"), rAtAbsorberEnd);
    //     fRegistry.fill(HIST("MCHMID/hChi2"), fwdtrack.chi2());
    //     fRegistry.fill(HIST("MCHMID/hChi2MFT"), chi2mft / ndf_mft);
    //     fRegistry.fill(HIST("MCHMID/hChi2MatchMCHMID"), fwdtrack.chi2MatchMCHMID());
    //     fRegistry.fill(HIST("MCHMID/hChi2MatchMCHMFT"), fwdtrack.chi2MatchMCHMFT());
    //     fRegistry.fill(HIST("MCHMID/hDCAxy2D"), dcaX, dcaY);
    //     fRegistry.fill(HIST("MCHMID/hDCAxy2DinSigma"), dcaX / std::sqrt(cXX), dcaY / std::sqrt(cYY));
    //     fRegistry.fill(HIST("MCHMID/hDCAxy"), dcaXY);
    //     fRegistry.fill(HIST("MCHMID/hDCAxyz"), dcaXY, dcaZ);
    //     fRegistry.fill(HIST("MCHMID/hDCAxyinSigma"), dcaXYinSigma);
    //     fRegistry.fill(HIST("MCHMID/hDCAxResolutionvsPt"), pt, std::sqrt(cXX) * 1e+4); // convert cm to um
    //     fRegistry.fill(HIST("MCHMID/hDCAyResolutionvsPt"), pt, std::sqrt(cYY) * 1e+4); // convert cm to um
    //     fRegistry.fill(HIST("MCHMID/hDCAxyResolutionvsPt"), pt, sigma_dcaXY * 1e+4);   // convert cm to um
    //   }
    // }
    return true;
  }

  SliceCache cache;
  Preslice<aod::FwdTracks> perCollision = o2::aod::fwdtrack::collisionId;
  Preslice<aod::FwdTrackAssoc> fwdtrackIndicesPerCollision = aod::track_association::collisionId;
  PresliceUnsorted<aod::FwdTrackAssoc> fwdtrackIndicesPerFwdTrack = aod::track_association::fwdtrackId;
  PresliceUnsorted<aod::FwdTracks> fwdtracksPerMCHTrack = aod::fwdtrack::matchMCHTrackId;

  // Filter trackFilter = o2::aod::fwdtrack::trackType == static_cast<uint8_t>(o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack) || o2::aod::fwdtrack::trackType == static_cast<uint8_t>(o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack);
  // using filteredMyFwdTracks = soa::Filtered<MyFwdTracks>;
  // Partition<MyFwdTracks> posTracks = o2::aod::fwdtrack::signed1Pt > 0.f;
  // Partition<MyFwdTracks> negTracks = o2::aod::fwdtrack::signed1Pt < 0.f;

  std::vector<Muon> vecMuons;

  void processRec(MyCollisions const& collisions, MyFwdTracks const& fwdtracks, aod::MFTTracks const&, aod::BCsWithTimestamps const&)
  {
    vecMuons.reserve(fwdtracks.size());

    for (const auto& collision : collisions) {
      auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);

      if (!collision.isSelected()) {
        continue;
      }

      auto fwdtracks_per_coll = fwdtracks.sliceBy(perCollision, collision.globalIndex());
      for (const auto& fwdtrack : fwdtracks_per_coll) {
        fillMuonInfo<false, MyFwdTracks, aod::MFTTracks>(collision, fwdtrack);
      }

      auto pos_muons_per_col = std::views::filter(vecMuons, [](Muon muon) { return muon.sign > 0; });
      auto neg_muons_per_col = std::views::filter(vecMuons, [](Muon muon) { return muon.sign < 0; });

      // ULS
      for (const auto& pos : pos_muons_per_col) {
        if (!isSelectedTight(pos)) { // pos is tag, neg is probe
          continue;
        }
        for (const auto& neg : neg_muons_per_col) {
          if (!isSelected(neg)) {
            continue;
          }
          ROOT::Math::PtEtaPhiMVector v1(pos.pt, pos.eta, pos.phi, o2::constants::physics::MassMuon); // tag
          ROOT::Math::PtEtaPhiMVector v2(neg.pt, neg.eta, neg.phi, o2::constants::physics::MassMuon); // probe
          ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
          if (neg.trackType == static_cast<uint8_t>(o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack)) {
            if (pos.matchMCHTrackId == neg.matchMCHTrackId || pos.matchMFTTrackId == neg.matchMFTTrackId) { // this should not happen in ULS. only for protection.
              continue;
            }
            fRegistry.fill(HIST("Pair/uls/gl_gl/hMvsPt"), v12.M(), v2.Pt());
          } else if (neg.trackType == static_cast<uint8_t>(o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack)) {
            if (pos.matchMCHTrackId == neg.globalIndex) { // this should not happen in ULS. only for protection.
              continue;
            }
            fRegistry.fill(HIST("Pair/uls/gl_sa/hMvsPt"), v12.M(), v2.Pt());
          }
          if (pairCuts.minMass < v12.M() && v12.M() < pairCuts.maxMass) {
            fillFwdTrackTable(collision, neg, fwdtracks.rawIteratorAt(neg.globalIndex));
          }
        } // end of neg
      } // end of pos

      // ULS
      for (const auto& neg : neg_muons_per_col) {
        if (!isSelectedTight(neg)) { // neg is tag, pos is probe
          continue;
        }
        for (const auto& pos : pos_muons_per_col) {
          if (!isSelected(pos)) {
            continue;
          }
          ROOT::Math::PtEtaPhiMVector v1(neg.pt, neg.eta, neg.phi, o2::constants::physics::MassMuon); // tag
          ROOT::Math::PtEtaPhiMVector v2(pos.pt, pos.eta, pos.phi, o2::constants::physics::MassMuon); // probe
          ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
          if (pos.trackType == static_cast<uint8_t>(o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack)) {
            if (pos.matchMCHTrackId == neg.matchMCHTrackId || pos.matchMFTTrackId == neg.matchMFTTrackId) { // this should not happen in ULS. only for protection.
              continue;
            }
            fRegistry.fill(HIST("Pair/uls/gl_gl/hMvsPt"), v12.M(), v2.Pt());
          } else if (pos.trackType == static_cast<uint8_t>(o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack)) {
            if (neg.matchMCHTrackId == pos.globalIndex) { // this should not happen in ULS. only for protection.
              continue;
            }
            fRegistry.fill(HIST("Pair/uls/gl_sa/hMvsPt"), v12.M(), v2.Pt());
          }
          if (pairCuts.minMass < v12.M() && v12.M() < pairCuts.maxMass) {
            fillFwdTrackTable(collision, pos, fwdtracks.rawIteratorAt(pos.globalIndex));
          }
        } // end of pos
      } // end of neg

    } // end of collision loop

    vecMuons.clear();
    vecMuons.shrink_to_fit();
  }
  PROCESS_SWITCH(skimmerPrimaryMuonQC, processRec, "process reconstructed info", false);

  void processRec_SWT(MyCollisionsWithSWT const& collisions, MyFwdTracks const& fwdtracks, aod::MFTTracks const&, aod::BCsWithTimestamps const&)
  {
    vecMuons.reserve(fwdtracks.size());

    for (const auto& collision : collisions) {
      const auto& bc = collision.template bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);

      if (!collision.isSelected()) {
        continue;
      }

      if (collision.swtaliastmp_raw() == 0) {
        continue;
      }

      // const auto& fwdtracks_per_coll = fwdtracks.sliceBy(perCollision, collision.globalIndex());
      // for (const auto& fwdtrack : fwdtracks_per_coll) {
      //   if (fwdtrack.trackType() != o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack && fwdtrack.trackType() != o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack) {
      //     continue;
      //   }

      //   if (!fillFwdTrackTable<false, filteredMyFwdTracks, aod::MFTTracks, false>(collision, fwdtrack, false)) {
      //     continue;
      //   }

      // } // end of fwdtrack loop
    } // end of collision loop

    vecMuons.clear();
    vecMuons.shrink_to_fit();
  }
  PROCESS_SWITCH(skimmerPrimaryMuonQC, processRec_SWT, "process reconstructed info only with standalone", false);

  using filteredMyFwdTracksMC = soa::Filtered<MyFwdTracksMC>;
  void processMC(soa::Join<MyCollisions, aod::McCollisionLabels> const& collisions, MyFwdTracksMC const& fwdtracks, MFTTracksMC const&, aod::BCsWithTimestamps const&, aod::McParticles const&)
  {
    vecMuons.reserve(fwdtracks.size());

    for (const auto& collision : collisions) {
      auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);
      if (!collision.isSelected()) {
        continue;
      }
      if (!collision.has_mcCollision()) {
        continue;
      }

      // auto fwdtracks_per_coll = fwdtracks.sliceBy(perCollision, collision.globalIndex());
      // for (const auto& fwdtrack : fwdtracks_per_coll) {
      //   if (!fwdtrack.has_mcParticle()) {
      //     continue;
      //   }
      //   if (fwdtrack.trackType() != o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack && fwdtrack.trackType() != o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack) {
      //     continue;
      //   }

      //   if (!fillFwdTrackTable<true, filteredMyFwdTracksMC, MFTTracksMC, false>(collision, fwdtrack, false)) {
      //     continue;
      //   }

      // } // end of fwdtrack loop
    } // end of collision loop

    vecMuons.clear();
    vecMuons.shrink_to_fit();
  }
  PROCESS_SWITCH(skimmerPrimaryMuonQC, processMC, "process reconstructed and MC info", false);

  void processDummy(aod::Collisions const&) {}
  PROCESS_SWITCH(skimmerPrimaryMuonQC, processDummy, "process dummy", true);
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

struct associateSameMuonElement {
  Produces<aod::EMGlobalMuonSelfIds> glmuon_same_ids;

  SliceCache cache;
  PresliceUnsorted<aod::EMPrimaryMuons> perMFTTrack = o2::aod::emprimarymuon::mfttrackId;
  PresliceUnsorted<aod::EMPrimaryMuons> perMCHTrack = o2::aod::emprimarymuon::mchtrackId;
  std::vector<int> selfIds_per_MFT;
  std::vector<int> selfIds_per_MCHMID;

  // Multiple MCH-MID tracks can match with the same MFTsa. This function is to reject such global muons.
  void process(aod::EMPrimaryMuons const& muons)
  {
    for (const auto& muon : muons) {
      if (muon.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack) {
        auto muons_with_same_mfttrackId = muons.sliceBy(perMFTTrack, muon.mfttrackId());
        auto muons_with_same_mchtrackId = muons.sliceBy(perMCHTrack, muon.mchtrackId());
        selfIds_per_MFT.reserve(muons_with_same_mfttrackId.size());
        selfIds_per_MCHMID.reserve(muons_with_same_mchtrackId.size());
        // LOGF(info, "muons_with_same_mchtrackId.size() = %d, muons_with_same_mfttrackId.size() = %d", muons_with_same_mchtrackId.size(), muons_with_same_mfttrackId.size());

        for (const auto& global_muon : muons_with_same_mfttrackId) {
          // LOGF(info, "same MFT: global_muon.globalIndex() = %d, global_muon.mchtrackId() = %d, global_muon.mfttrackId() = %d, global_muon.collisionId() = %d", global_muon.globalIndex(), global_muon.mchtrackId(), global_muon.mfttrackId(), global_muon.collisionId());
          if (global_muon.globalIndex() == muon.globalIndex()) { // don't store myself.
            continue;
          }
          if (global_muon.collisionId() == muon.collisionId()) { // the same global muon is repeatedly stored and associated to different collisions if FTTCA is used.
            selfIds_per_MFT.emplace_back(global_muon.globalIndex());
          }
        }

        for (const auto& global_muon : muons_with_same_mchtrackId) {
          // LOGF(info, "same MCH: global_muon.globalIndex() = %d, global_muon.mchtrackId() = %d, global_muon.mfttrackId() = %d, global_muon.collisionId() = %d", global_muon.globalIndex(), global_muon.mchtrackId(), global_muon.mfttrackId(), global_muon.collisionId());
          if (global_muon.globalIndex() == muon.globalIndex()) { // don't store myself.
            continue;
          }
          if (global_muon.collisionId() == muon.collisionId()) { // the same global muon is repeatedly stored and associated to different collisions if FTTCA is used.
            selfIds_per_MCHMID.emplace_back(global_muon.globalIndex());
          }
        }

        glmuon_same_ids(selfIds_per_MCHMID, selfIds_per_MFT);
        selfIds_per_MFT.clear();
        selfIds_per_MFT.shrink_to_fit();
        selfIds_per_MCHMID.clear();
        selfIds_per_MCHMID.shrink_to_fit();
      } else {
        glmuon_same_ids(std::vector<int>{}, std::vector<int>{}); // empty for standalone muons
        selfIds_per_MFT.clear();
        selfIds_per_MFT.shrink_to_fit();
        selfIds_per_MCHMID.clear();
        selfIds_per_MCHMID.shrink_to_fit();
      }
    } // end of muon loop
  }
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<skimmerPrimaryMuonQC>(cfgc, TaskName{"skimmer-primary-muon-qc"}),
    adaptAnalysisTask<associateAmbiguousMuon>(cfgc, TaskName{"associate-ambiguous-muon"}),
    adaptAnalysisTask<associateSameMuonElement>(cfgc, TaskName{"associate-same-muon-element"})};
}
