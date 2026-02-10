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

#include <map>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::soa;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;
using namespace o2::aod::fwdtrackutils;

struct skimmerPrimaryMuon {
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
  Produces<aod::EMPrimaryMuonsMatchMC> emprimarymuonsmatchmc;

  // Configurables
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
  Configurable<bool> fillQAHistograms{"fillQAHistograms", false, "flag to fill QA histograms"};
  Configurable<float> minPt{"minPt", 0.01, "min pt for muon"};
  Configurable<float> maxPt{"maxPt", 1e+10, "max pt for muon"};
  Configurable<float> minEtaSA{"minEtaSA", -4.0, "min. eta acceptance for MCH-MID"};
  Configurable<float> maxEtaSA{"maxEtaSA", -2.5, "max. eta acceptance for MCH-MID"};
  Configurable<float> minEtaGL{"minEtaGL", -3.6, "min. eta acceptance for MFT-MCH-MID"};
  Configurable<float> maxEtaGL{"maxEtaGL", -2.5, "max. eta acceptance for MFT-MCH-MID"};
  Configurable<float> minRabsGL{"minRabsGL", 17.6, "min. R at absorber end for global muon (min. eta = -3.6)"}; // std::tan(2.f * std::atan(std::exp(- -3.6)) ) * -505.
  Configurable<float> minRabs{"minRabs", 17.6, "min. R at absorber end"};
  Configurable<float> midRabs{"midRabs", 26.5, "middle R at absorber end for pDCA cut"};
  Configurable<float> maxRabs{"maxRabs", 89.5, "max. R at absorber end"};
  Configurable<float> maxDCAxy{"maxDCAxy", 1, "max. DCAxy for global muons"};
  Configurable<float> maxPDCAforLargeR{"maxPDCAforLargeR", 324.f, "max. pDCA for large R at absorber end"};
  Configurable<float> maxPDCAforSmallR{"maxPDCAforSmallR", 594.f, "max. pDCA for small R at absorber end"};
  Configurable<float> maxMatchingChi2MCHMFT{"maxMatchingChi2MCHMFT", 50.f, "max. chi2 for MCH-MFT matching"};
  Configurable<float> maxChi2SA{"maxChi2SA", 1e+6, "max. chi2 for standalone muon"};
  Configurable<float> maxChi2GL{"maxChi2GL", 10, "max. chi2 for global muon"};
  Configurable<bool> refitGlobalMuon{"refitGlobalMuon", true, "flag to refit global muon"};
  Configurable<float> matchingZ{"matchingZ", -77.5, "z position where matching is performed"};
  Configurable<int> minNmuon{"minNmuon", 0, "min number of muon candidates per collision"};
  Configurable<float> maxDEta{"maxDEta", 1e+10f, "max. deta between MFT-MCH-MID and MCH-MID"};
  Configurable<float> maxDPhi{"maxDPhi", 1e+10f, "max. dphi between MFT-MCH-MID and MCH-MID"};
  Configurable<bool> cfgApplyPreselectionInBestMatch{"cfgApplyPreselectionInBestMatch", false, "flag to apply preselection in find best match function"};
  Configurable<float> cfgSlope_dr_chi2MatchMFTMCH{"cfgSlope_dr_chi2MatchMFTMCH", -0.15 / 30, "slope of chiMatchMCHMFT vs. dR"};
  Configurable<float> cfgIntercept_dr_chi2MatchMFTMCH{"cfgIntercept_dr_chi2MatchMFTMCH", 1e+10f, "intercept of chiMatchMCHMFT vs. dR"};
  Configurable<float> cfgPeakPosition_chi2MatchMFTMCH{"cfgPeakPosition_chi2MatchMFTMCH", 0.f, "peak position of chiMatchMCHMFT distribution"}; // 2
  Configurable<float> cfgPeakPosition_dr{"cfgPeakPosition_dr", 0.f, "peak position of dr distribution"};                                       // 0.01

  // for z shift for propagation
  Configurable<bool> cfgApplyZShiftFromCCDB{"cfgApplyZShiftFromCCDB", false, "flag to apply z shift"};
  Configurable<std::string> cfgZShiftPath{"cfgZShiftPath", "Users/m/mcoquet/ZShift", "CCDB path for z shift to apply to forward tracks"};
  Configurable<float> cfgManualZShift{"cfgManualZShift", 0, "manual z-shift for propagation of global muon to PV"};

  o2::ccdb::CcdbApi ccdbApi;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  int mRunNumber = 0;
  float mBz = 0;
  float mZShift = 0;

  HistogramRegistry fRegistry{"output", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};
  static constexpr std::string_view muon_types[5] = {"MFTMCHMID/", "MFTMCHMIDOtherMatch/", "MFTMCH/", "MCHMID/", "MCH/"};

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
    auto hMuonType = fRegistry.add<TH1>("hMuonType", "muon type", kTH1F, {{5, -0.5f, 4.5f}}, false);
    hMuonType->GetXaxis()->SetBinLabel(1, "MFT-MCH-MID (global muon)");
    hMuonType->GetXaxis()->SetBinLabel(2, "MFT-MCH-MID (global muon other match)");
    hMuonType->GetXaxis()->SetBinLabel(3, "MFT-MCH");
    hMuonType->GetXaxis()->SetBinLabel(4, "MCH-MID");
    hMuonType->GetXaxis()->SetBinLabel(5, "MCH standalone");

    fRegistry.add("MFTMCHMID/hPt", "pT;p_{T} (GeV/c)", kTH1F, {{100, 0.0f, 10}}, false);
    fRegistry.add("MFTMCHMID/hEtaPhi", "#eta vs. #varphi;#varphi (rad.);#eta", kTH2F, {{180, 0, 2 * M_PI}, {60, -5.f, -2.f}}, false);
    fRegistry.add("MFTMCHMID/hEtaPhi_MatchedMCHMID", "#eta vs. #varphi;#varphi (rad.);#eta", kTH2F, {{180, 0, 2 * M_PI}, {60, -5.f, -2.f}}, false);
    fRegistry.add("MFTMCHMID/hDeltaPt_Pt", "#Deltap_{T}/p_{T} vs. p_{T};p_{T}^{gl} (GeV/c);(p_{T}^{sa} - p_{T}^{gl})/p_{T}^{gl}", kTH2F, {{100, 0, 10}, {200, -0.5, +0.5}}, false);
    fRegistry.add("MFTMCHMID/hDeltaEta_Pt", "#Delta#eta vs. p_{T};p_{T}^{gl} (GeV/c);#Delta#eta", kTH2F, {{100, 0, 10}, {200, -0.5, +0.5}}, false);
    fRegistry.add("MFTMCHMID/hDeltaPhi_Pt", "#Delta#varphi vs. p_{T};p_{T}^{gl} (GeV/c);#Delta#varphi (rad.)", kTH2F, {{100, 0, 10}, {200, -0.5, +0.5}}, false);
    fRegistry.add("MFTMCHMID/hDeltaEtaAtMP_Pt", "#Delta#eta vs. p_{T} at MP;p_{T}^{gl} (GeV/c);#Delta#eta", kTH2F, {{100, 0, 10}, {200, -0.5, +0.5}}, false);
    fRegistry.add("MFTMCHMID/hDeltaPhiAtMP_Pt", "#Delta#varphi vs. p_{T} at MP;p_{T}^{gl} (GeV/c);#Delta#varphi (rad.)", kTH2F, {{100, 0, 10}, {200, -0.5, +0.5}}, false);
    fRegistry.add("MFTMCHMID/hSign", "sign;sign", kTH1F, {{3, -1.5, +1.5}}, false);
    fRegistry.add("MFTMCHMID/hNclusters", "Nclusters;Nclusters", kTH1F, {{21, -0.5f, 20.5}}, false);
    fRegistry.add("MFTMCHMID/hNclustersMFT", "NclustersMFT;Nclusters MFT", kTH1F, {{11, -0.5f, 10.5}}, false);
    fRegistry.add("MFTMCHMID/hRatAbsorberEnd", "R at absorber end;R at absorber end (cm)", kTH1F, {{100, 0.0f, 100}}, false);
    fRegistry.add("MFTMCHMID/hPDCA_Rabs", "pDCA vs. Rabs;R at absorber end (cm);p #times DCA (GeV/c #upoint cm)", kTH2F, {{100, 0, 100}, {100, 0.0f, 1000}}, false);
    fRegistry.add("MFTMCHMID/hChi2", "chi2;chi2/ndf", kTH1F, {{100, 0.0f, 10}}, false);
    fRegistry.add("MFTMCHMID/hChi2MFT", "chi2 MFT;chi2 MFT/ndf", kTH1F, {{100, 0.0f, 10}}, false);
    fRegistry.add("MFTMCHMID/hChi2MatchMCHMID", "chi2 match MCH-MID;chi2", kTH1F, {{100, 0.0f, 100}}, false);
    fRegistry.add("MFTMCHMID/hChi2MatchMCHMFT", "chi2 match MCH-MFT;chi2", kTH1F, {{100, 0.0f, 100}}, false);
    fRegistry.add("MFTMCHMID/hDCAxy2D", "DCA x vs. y;DCA_{x} (cm);DCA_{y} (cm)", kTH2F, {{200, -1, 1}, {200, -1, +1}}, false);
    fRegistry.add("MFTMCHMID/hDCAxy2DinSigma", "DCA x vs. y in sigma;DCA_{x} (#sigma);DCA_{y} (#sigma)", kTH2F, {{200, -10, 10}, {200, -10, +10}}, false);
    fRegistry.add("MFTMCHMID/hDCAxy", "DCAxy;DCA_{xy} (cm);", kTH1F, {{100, 0, 1}}, false);
    fRegistry.add("MFTMCHMID/hDCAxyz", "DCA xy vs. z;DCA_{xy} (cm);DCA_{z} (cm)", kTH2F, {{100, 0, 1}, {200, -0.1, 0.1}}, false);
    fRegistry.add("MFTMCHMID/hDCAxyinSigma", "DCAxy in sigma;DCA_{xy} (#sigma);", kTH1F, {{100, 0, 10}}, false);
    fRegistry.add("MFTMCHMID/hDCAx_PosZ", "DCAx vs. posZ;Z_{vtx} (cm);DCA_{x} (cm)", kTH2F, {{200, -10, +10}, {400, -0.2, +0.2}}, false);
    fRegistry.add("MFTMCHMID/hDCAy_PosZ", "DCAy vs. posZ;Z_{vtx} (cm);DCA_{y} (cm)", kTH2F, {{200, -10, +10}, {400, -0.2, +0.2}}, false);
    fRegistry.add("MFTMCHMID/hDCAx_Phi", "DCAx vs. #varphi;#varphi (rad.);DCA_{x} (cm)", kTH2F, {{90, 0, 2 * M_PI}, {400, -0.2, +0.2}}, false);
    fRegistry.add("MFTMCHMID/hDCAy_Phi", "DCAy vs. #varphi;#varphi (rad.);DCA_{y} (cm)", kTH2F, {{90, 0, 2 * M_PI}, {400, -0.2, +0.2}}, false);
    fRegistry.add("MFTMCHMID/hNmu", "#mu multiplicity;N_{#mu} per collision", kTH1F, {{21, -0.5, 20.5}}, false);
    fRegistry.add("MFTMCHMID/hdR_Chi2MatchMCHMFT", "dr vs. matching chi2 MCH-MFT;chi2 match MCH-MFT;#DeltaR", kTH2F, {{200, 0, 50}, {200, 0, 0.5}}, false);
    fRegistry.add("MFTMCHMID/hdR_Chi2", "dr vs. chi2;global chi2/ndf;#DeltaR", kTH2F, {{100, 0, 10}, {200, 0, 0.5}}, false);
    fRegistry.add("MFTMCHMID/hChi2_Chi2MatchMCHMFT", "chi2 vs. matching chi2 MCH-MFT;chi2 match MCH-MFT;global chi2/ndf", kTH2F, {{200, 0, 50}, {100, 0, 10}}, false);

    fRegistry.addClone("MFTMCHMID/", "MCHMID/");
    fRegistry.add("MFTMCHMID/hDCAxResolutionvsPt", "DCA_{x} vs. p_{T};p_{T} (GeV/c);DCA_{x} resolution (#mum);", kTH2F, {{100, 0, 10.f}, {500, 0, 500}}, false);
    fRegistry.add("MFTMCHMID/hDCAyResolutionvsPt", "DCA_{y} vs. p_{T};p_{T} (GeV/c);DCA_{y} resolution (#mum);", kTH2F, {{100, 0, 10.f}, {500, 0, 500}}, false);
    fRegistry.add("MFTMCHMID/hDCAxyResolutionvsPt", "DCA_{xy} vs. p_{T};p_{T} (GeV/c);DCA_{y} resolution (#mum);", kTH2F, {{100, 0, 10.f}, {500, 0, 500}}, false);
    fRegistry.add("MCHMID/hDCAxResolutionvsPt", "DCA_{x} vs. p_{T};p_{T} (GeV/c);DCA_{x} resolution (#mum);", kTH2F, {{100, 0, 10.f}, {500, 0, 5e+5}}, false);
    fRegistry.add("MCHMID/hDCAyResolutionvsPt", "DCA_{y} vs. p_{T};p_{T} (GeV/c);DCA_{y} resolution (#mum);", kTH2F, {{100, 0, 10.f}, {500, 0, 5e+5}}, false);
    fRegistry.add("MCHMID/hDCAxyResolutionvsPt", "DCA_{xy} vs. p_{T};p_{T} (GeV/c);DCA_{y} resolution (#mum);", kTH2F, {{100, 0, 10.f}, {500, 0, 5e+5}}, false);
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
      if (maxChi2GL < chi2_per_ndf) {
        return false;
      }
      if (rAtAbsorberEnd < minRabsGL || maxRabs < rAtAbsorberEnd) {
        return false;
      }
    } else if (trackType == static_cast<uint8_t>(o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack)) {
      if (eta < minEtaSA || maxEtaSA < eta) {
        return false;
      }
      if (maxChi2SA < chi2_per_ndf) {
        return false;
      }
    } else {
      return false;
    }

    return true;
  }

  template <bool isMC, bool withMFTCov, typename TFwdTracks, typename TMFTTracks, bool fillTable, typename TCollision, typename TFwdTrack, typename TMFTTracksCov>
  bool fillFwdTrackTable(TCollision const& collision, TFwdTrack fwdtrack, TMFTTracksCov const& mftCovs, const bool isAmbiguous)
  {
    if (fwdtrack.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack && fwdtrack.chi2MatchMCHMFT() > maxMatchingChi2MCHMFT) {
      return false;
    } // Users have to decide the best match between MFT and MCH-MID at analysis level. The same global muon is repeatedly stored.

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
    float dcaZ = propmuonAtPV.getZ() - collision.posZ();
    float dcaXY = std::sqrt(dcaX * dcaX + dcaY * dcaY);
    float rAtAbsorberEnd = fwdtrack.rAtAbsorberEnd(); // this works only for GlobalMuonTrack
    float cXX = propmuonAtPV.getSigma2X();
    float cYY = propmuonAtPV.getSigma2Y();
    float cXY = propmuonAtPV.getSigmaXY();

    float det = cXX * cYY - cXY * cXY; // determinanat
    float dcaXYinSigma = 999.f;
    if (det < 0) {
      dcaXYinSigma = 999.f;
    } else {
      dcaXYinSigma = std::sqrt(std::fabs((dcaX * dcaX * cYY + dcaY * dcaY * cXX - 2.f * dcaX * dcaY * cXY) / det / 2.f)); // dca xy in sigma
    }
    float sigma_dcaXY = dcaXY / dcaXYinSigma;

    float pDCA = propmuonAtPV.getP() * dcaXY;
    int nClustersMFT = 0;
    float ptMatchedMCHMID = propmuonAtPV.getPt();
    float etaMatchedMCHMID = propmuonAtPV.getEta();
    float phiMatchedMCHMID = propmuonAtPV.getPhi();
    o2::math_utils::bringTo02Pi(phiMatchedMCHMID);
    // float x = fwdtrack.x();
    // float y = fwdtrack.y();
    // float z = fwdtrack.z();
    // float tgl = fwdtrack.tgl();
    float chi2mft = 0.f;
    uint64_t mftClusterSizesAndTrackFlags = 0;
    int ndf_mchmft = 1;
    int ndf_mft = 1;

    float etaMatchedMCHMIDatMP = 999.f;
    float phiMatchedMCHMIDatMP = 999.f;
    float etaMatchedMFTatMP = 999.f;
    float phiMatchedMFTatMP = 999.f;

    float deta = 999.f;
    float dphi = 999.f;
    bool isCorrectMatchMFTMCH = true; // by default, it is true. it is evaluated for global muons in MC.

    if (fwdtrack.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack) {
      // apply r-absorber cut here to minimize the number of calling propagateMuon.
      if (fwdtrack.rAtAbsorberEnd() < minRabsGL || maxRabs < fwdtrack.rAtAbsorberEnd()) {
        return false;
      }

      // apply dca cut here to minimize the number of calling propagateMuon.
      if (maxDCAxy < dcaXY) {
        return false;
      }

      auto mchtrack = fwdtrack.template matchMCHTrack_as<TFwdTracks>(); // MCH-MID
      auto mfttrack = fwdtrack.template matchMFTTrack_as<TMFTTracks>(); // MFTsa

      if constexpr (isMC) {
        if (!mfttrack.has_mcParticle() || !mchtrack.has_mcParticle() || !fwdtrack.has_mcParticle()) {
          return false;
        }
        // auto mcParticle_MFTMCHMID = fwdtrack.template mcParticle_as<aod::McParticles>(); // this is identical to mcParticle_MCHMID
        auto mcParticle_MCHMID = mchtrack.template mcParticle_as<aod::McParticles>(); // this is identical to mcParticle_MFTMCHMID
        auto mcParticle_MFT = mfttrack.template mcParticle_as<aod::McParticles>();
        isCorrectMatchMFTMCH = static_cast<bool>(mcParticle_MCHMID.globalIndex() == mcParticle_MFT.globalIndex());
      }

      nClustersMFT = mfttrack.nClusters();
      mftClusterSizesAndTrackFlags = mfttrack.mftClusterSizesAndTrackFlags();
      ndf_mchmft = 2.f * (mchtrack.nClusters() + nClustersMFT) - 5.f;
      ndf_mft = 2.f * nClustersMFT - 5.f;
      chi2mft = mfttrack.chi2();
      // chi2mft = mfttrack.chi2() / (2.f * nClustersMFT - 5.f);

      // apply chi2/ndf cut here to minimize the number of calling propagateMuon.
      if (maxChi2GL < fwdtrack.chi2() / ndf_mchmft) {
        return false;
      }

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

      if constexpr (withMFTCov) {
        auto mfttrackcov = mftCovs.rawIteratorAt(map_mfttrackcovs[mfttrack.globalIndex()]);
        auto muonAtMP = propagateMuon(mchtrack, mchtrack, collision, propagationPoint::kToMatchingPlane, matchingZ, mBz, mZShift); // propagated to matching plane
        o2::track::TrackParCovFwd mftsaAtMP = getTrackParCovFwdShift(mfttrack, mZShift, mfttrackcov);                              // values at innermost update
        mftsaAtMP.propagateToZhelix(matchingZ, mBz);                                                                               // propagated to matching plane
        etaMatchedMFTatMP = mftsaAtMP.getEta();
        phiMatchedMFTatMP = mftsaAtMP.getPhi();
        etaMatchedMCHMIDatMP = muonAtMP.getEta();
        phiMatchedMCHMIDatMP = muonAtMP.getPhi();
        o2::math_utils::bringTo02Pi(phiMatchedMCHMIDatMP);
        o2::math_utils::bringTo02Pi(phiMatchedMFTatMP);

        o2::track::TrackParCovFwd mftsa = getTrackParCovFwdShift(mfttrack, mZShift, mfttrackcov);                                  // values at innermost update
        o2::dataformats::GlobalFwdTrack globalMuonRefit = o2::aod::fwdtrackutils::refitGlobalMuonCov(propmuonAtPV_Matched, mftsa); // this is track at IU.
        auto globalMuon = o2::aod::fwdtrackutils::propagateTrackParCovFwd(globalMuonRefit, fwdtrack.trackType(), collision, propagationPoint::kToVertex, matchingZ, mBz);
        pt = globalMuon.getPt();
        eta = globalMuon.getEta();
        phi = globalMuon.getPhi();
        o2::math_utils::bringTo02Pi(phi);

        cXX = globalMuon.getSigma2X();
        cYY = globalMuon.getSigma2Y();
        cXY = globalMuon.getSigmaXY();
        dcaX = globalMuon.getX() - collision.posX();
        dcaY = globalMuon.getY() - collision.posY();
        dcaZ = globalMuon.getZ() - collision.posZ();
        dcaXY = std::sqrt(dcaX * dcaX + dcaY * dcaY);
        det = cXX * cYY - cXY * cXY; // determinanat
        dcaXYinSigma = 999.f;
        if (det < 0) {
          dcaXYinSigma = 999.f;
        } else {
          dcaXYinSigma = std::sqrt(std::fabs((dcaX * dcaX * cYY + dcaY * dcaY * cXX - 2.f * dcaX * dcaY * cXY) / det / 2.f)); // dca xy in sigma
        }
        sigma_dcaXY = dcaXY / dcaXYinSigma;
      }

      deta = etaMatchedMCHMID - eta;
      dphi = phiMatchedMCHMID - phi;
      o2::math_utils::bringToPMPi(dphi);

      if (std::sqrt(std::pow(deta / maxDEta, 2) + std::pow(dphi / maxDPhi, 2)) > 1.f) {
        return false;
      }

      float dr = std::sqrt(deta * deta + dphi * dphi);
      if (cfgSlope_dr_chi2MatchMFTMCH * fwdtrack.chi2MatchMCHMFT() + cfgIntercept_dr_chi2MatchMFTMCH < dr) {
        return false;
      }

      if (refitGlobalMuon) {
        pt = propmuonAtPV_Matched.getP() * std::sin(2.f * std::atan(std::exp(-eta)));
      }
    } else if (fwdtrack.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack) {
      o2::dataformats::GlobalFwdTrack propmuonAtRabs = propagateMuon(fwdtrack, fwdtrack, collision, propagationPoint::kToRabs, matchingZ, mBz, mZShift); // this is necessary only for MuonStandaloneTrack
      float xAbs = propmuonAtRabs.getX();
      float yAbs = propmuonAtRabs.getY();
      rAtAbsorberEnd = std::sqrt(xAbs * xAbs + yAbs * yAbs); // Redo propagation only for muon tracks // propagation of MFT tracks alredy done in reconstruction

      o2::dataformats::GlobalFwdTrack propmuonAtDCA = propagateMuon(fwdtrack, fwdtrack, collision, propagationPoint::kToDCA, matchingZ, mBz, mZShift);
      cXX = propmuonAtDCA.getSigma2X();
      cYY = propmuonAtDCA.getSigma2Y();
      cXY = propmuonAtDCA.getSigmaXY();
      dcaX = propmuonAtDCA.getX() - collision.posX();
      dcaY = propmuonAtDCA.getY() - collision.posY();
      dcaZ = propmuonAtDCA.getZ() - collision.posZ();
      dcaXY = std::sqrt(dcaX * dcaX + dcaY * dcaY);
      pDCA = fwdtrack.p() * dcaXY;

      det = cXX * cYY - cXY * cXY; // determinanat
      dcaXYinSigma = 999.f;
      if (det < 0) {
        dcaXYinSigma = 999.f;
      } else {
        dcaXYinSigma = std::sqrt(std::fabs((dcaX * dcaX * cYY + dcaY * dcaY * cXX - 2.f * dcaX * dcaY * cXY) / det / 2.f)); // dca xy in sigma
      }
      sigma_dcaXY = dcaXY / dcaXYinSigma;
    } else {
      return false;
    }

    if (!isSelected(pt, eta, rAtAbsorberEnd, pDCA, fwdtrack.chi2() / ndf_mchmft, fwdtrack.trackType(), dcaXY)) {
      return false;
    }

    if constexpr (fillTable) {
      float dpt = (ptMatchedMCHMID - pt) / pt;

      float detaMP = etaMatchedMCHMIDatMP - etaMatchedMFTatMP;
      float dphiMP = phiMatchedMCHMIDatMP - phiMatchedMFTatMP;
      o2::math_utils::bringToPMPi(dphiMP);

      bool isAssociatedToMPC = fwdtrack.collisionId() == collision.globalIndex();
      // LOGF(info, "isAmbiguous = %d, isAssociatedToMPC = %d, fwdtrack.globalIndex() = %d, fwdtrack.collisionId() = %d, collision.globalIndex() = %d", isAmbiguous, isAssociatedToMPC, fwdtrack.globalIndex(), fwdtrack.collisionId(), collision.globalIndex());

      emprimarymuons(collision.globalIndex(), fwdtrack.globalIndex(), fwdtrack.matchMFTTrackId(), fwdtrack.matchMCHTrackId(), fwdtrack.trackType(),
                     pt, eta, phi, fwdtrack.sign(), dcaX, dcaY, cXX, cYY, cXY, ptMatchedMCHMID, etaMatchedMCHMID, phiMatchedMCHMID,
                     etaMatchedMCHMIDatMP, phiMatchedMCHMIDatMP, etaMatchedMFTatMP, phiMatchedMFTatMP,
                     fwdtrack.nClusters(), pDCA, rAtAbsorberEnd, fwdtrack.chi2(), fwdtrack.chi2MatchMCHMID(), fwdtrack.chi2MatchMCHMFT(),
                     fwdtrack.mchBitMap(), fwdtrack.midBitMap(), fwdtrack.midBoards(), mftClusterSizesAndTrackFlags, chi2mft, isAssociatedToMPC, isAmbiguous);

      const auto& fwdcov = propmuonAtPV.getCovariances(); // covatiance matrix at PV
      emprimarymuonscov(
        fwdcov(0, 0),
        fwdcov(0, 1), fwdcov(1, 1),
        fwdcov(2, 0), fwdcov(2, 1), fwdcov(2, 2),
        fwdcov(3, 0), fwdcov(3, 1), fwdcov(3, 2), fwdcov(3, 3),
        fwdcov(4, 0), fwdcov(4, 1), fwdcov(4, 2), fwdcov(4, 3), fwdcov(4, 4));

      // See definition DataFormats/Reconstruction/include/ReconstructionDataFormats/TrackFwd.h
      // Covariance matrix of track parameters, ordered as follows:
      //  <X,X>         <Y,X>           <PHI,X>       <TANL,X>        <INVQPT,X>
      //  <X,Y>         <Y,Y>           <PHI,Y>       <TANL,Y>        <INVQPT,Y>
      // <X,PHI>       <Y,PHI>         <PHI,PHI>     <TANL,PHI>      <INVQPT,PHI>
      // <X,TANL>      <Y,TANL>       <PHI,TANL>     <TANL,TANL>     <INVQPT,TANL>
      // <X,INVQPT>   <Y,INVQPT>     <PHI,INVQPT>   <TANL,INVQPT>   <INVQPT,INVQPT>

      emprimarymuonsmatchmc(isCorrectMatchMFTMCH);

      if (fillQAHistograms) {
        fRegistry.fill(HIST("hMuonType"), fwdtrack.trackType());
        if (fwdtrack.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack) {
          fRegistry.fill(HIST("MFTMCHMID/hPt"), pt);
          fRegistry.fill(HIST("MFTMCHMID/hEtaPhi"), phi, eta);
          fRegistry.fill(HIST("MFTMCHMID/hEtaPhi_MatchedMCHMID"), phiMatchedMCHMID, etaMatchedMCHMID);
          fRegistry.fill(HIST("MFTMCHMID/hDeltaPt_Pt"), pt, dpt);
          fRegistry.fill(HIST("MFTMCHMID/hDeltaEta_Pt"), pt, deta);
          fRegistry.fill(HIST("MFTMCHMID/hDeltaPhi_Pt"), pt, dphi);
          fRegistry.fill(HIST("MFTMCHMID/hDeltaEtaAtMP_Pt"), pt, detaMP);
          fRegistry.fill(HIST("MFTMCHMID/hDeltaPhiAtMP_Pt"), pt, dphiMP);
          fRegistry.fill(HIST("MFTMCHMID/hSign"), fwdtrack.sign());
          fRegistry.fill(HIST("MFTMCHMID/hNclusters"), fwdtrack.nClusters());
          fRegistry.fill(HIST("MFTMCHMID/hNclustersMFT"), nClustersMFT);
          fRegistry.fill(HIST("MFTMCHMID/hPDCA_Rabs"), rAtAbsorberEnd, pDCA);
          fRegistry.fill(HIST("MFTMCHMID/hRatAbsorberEnd"), rAtAbsorberEnd);
          fRegistry.fill(HIST("MFTMCHMID/hChi2"), fwdtrack.chi2() / ndf_mchmft);
          fRegistry.fill(HIST("MFTMCHMID/hChi2MFT"), chi2mft / ndf_mft);
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
          fRegistry.fill(HIST("MFTMCHMID/hDCAx_PosZ"), collision.posZ(), dcaX);
          fRegistry.fill(HIST("MFTMCHMID/hDCAy_PosZ"), collision.posZ(), dcaY);
          fRegistry.fill(HIST("MFTMCHMID/hDCAx_Phi"), phi, dcaX);
          fRegistry.fill(HIST("MFTMCHMID/hDCAy_Phi"), phi, dcaY);
        } else if (fwdtrack.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack) {
          fRegistry.fill(HIST("MCHMID/hPt"), pt);
          fRegistry.fill(HIST("MCHMID/hEtaPhi"), phi, eta);
          fRegistry.fill(HIST("MCHMID/hEtaPhi_MatchedMCHMID"), phiMatchedMCHMID, etaMatchedMCHMID);
          fRegistry.fill(HIST("MCHMID/hDeltaPt_Pt"), pt, dpt);
          fRegistry.fill(HIST("MCHMID/hDeltaEta_Pt"), pt, deta);
          fRegistry.fill(HIST("MCHMID/hDeltaPhi_Pt"), pt, dphi);
          fRegistry.fill(HIST("MCHMID/hDeltaEtaAtMP_Pt"), pt, detaMP);
          fRegistry.fill(HIST("MCHMID/hDeltaPhiAtMP_Pt"), pt, dphiMP);
          fRegistry.fill(HIST("MCHMID/hSign"), fwdtrack.sign());
          fRegistry.fill(HIST("MCHMID/hNclusters"), fwdtrack.nClusters());
          fRegistry.fill(HIST("MCHMID/hNclustersMFT"), nClustersMFT);
          fRegistry.fill(HIST("MCHMID/hPDCA_Rabs"), rAtAbsorberEnd, pDCA);
          fRegistry.fill(HIST("MCHMID/hRatAbsorberEnd"), rAtAbsorberEnd);
          fRegistry.fill(HIST("MCHMID/hChi2"), fwdtrack.chi2());
          fRegistry.fill(HIST("MCHMID/hChi2MFT"), chi2mft / ndf_mft);
          fRegistry.fill(HIST("MCHMID/hChi2MatchMCHMID"), fwdtrack.chi2MatchMCHMID());
          fRegistry.fill(HIST("MCHMID/hChi2MatchMCHMFT"), fwdtrack.chi2MatchMCHMFT());
          fRegistry.fill(HIST("MCHMID/hDCAxy2D"), dcaX, dcaY);
          fRegistry.fill(HIST("MCHMID/hDCAxy2DinSigma"), dcaX / std::sqrt(cXX), dcaY / std::sqrt(cYY));
          fRegistry.fill(HIST("MCHMID/hDCAxy"), dcaXY);
          fRegistry.fill(HIST("MCHMID/hDCAxyz"), dcaXY, dcaZ);
          fRegistry.fill(HIST("MCHMID/hDCAxyinSigma"), dcaXYinSigma);
          fRegistry.fill(HIST("MCHMID/hDCAxResolutionvsPt"), pt, std::sqrt(cXX) * 1e+4); // convert cm to um
          fRegistry.fill(HIST("MCHMID/hDCAyResolutionvsPt"), pt, std::sqrt(cYY) * 1e+4); // convert cm to um
          fRegistry.fill(HIST("MCHMID/hDCAxyResolutionvsPt"), pt, sigma_dcaXY * 1e+4);   // convert cm to um
        }
      }
    }
    return true;
  }

  std::unordered_map<int, int> map_mfttrackcovs;
  std::vector<std::tuple<int, int, int>> vec_min_chi2MatchMCHMFT; // std::pair<globalIndex of global muon, globalIndex of matched MCH-MID, globalIndex of MFT> -> chi2MatchMCHMFT;
  // std::vector<std::tuple<int, int, int>> vec_min_dr;              // std::pair<globalIndex of global muon, globalIndex of matched MCH-MID, globalIndex of MFT> -> dr;
  // std::vector<std::tuple<int, int, int>> vec_min_2d;              // std::pair<globalIndex of global muon, globalIndex of matched MCH-MID, globalIndex of MFT> -> dr + chi2MatchMCHMFT;

  template <bool isMC, typename TCollision, typename TFwdTrack, typename TFwdTracks, typename TMFTTracks>
  void findBestMatchPerMCHMID(TCollision const& collision, TFwdTrack const& fwdtrack, TFwdTracks const& fwdtracks, TMFTTracks const&)
  {
    if (fwdtrack.trackType() != o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack) {
      return;
    }
    if constexpr (isMC) {
      if (!fwdtrack.has_mcParticle()) {
        return;
      }
    }

    const auto& muons_per_MCHMID = fwdtracks.sliceBy(fwdtracksPerMCHTrack, fwdtrack.globalIndex());
    // LOGF(info, "stanadalone: muon.globalIndex() = %d, muon.chi2MatchMCHMFT() = %f", muon.globalIndex(), muon.chi2MatchMCHMFT());
    // LOGF(info, "muons_per_MCHMID.size() = %d", muons_per_MCHMID.size());

    o2::dataformats::GlobalFwdTrack propmuonAtPV_Matched = propagateMuon(fwdtrack, fwdtrack, collision, propagationPoint::kToVertex, matchingZ, mBz, mZShift);
    float etaMatchedMCHMID = propmuonAtPV_Matched.getEta();
    float phiMatchedMCHMID = propmuonAtPV_Matched.getPhi();
    o2::math_utils::bringTo02Pi(phiMatchedMCHMID);

    // float min_chi2MatchMCHMFT = 1e+10, min_dr = 1e+10, min_distance_2d = 1e+10;
    float min_chi2MatchMCHMFT = 1e+10;
    std::tuple<int, int, int> tupleIds_at_min_chi2mftmch;
    // std::tuple<int, int, int> tupleIds_at_min_dr;
    // std::tuple<int, int, int> tupleIds_at_min_distance_2d;
    for (const auto& muon_tmp : muons_per_MCHMID) {
      if (muon_tmp.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack) {
        auto tupleId = std::make_tuple(muon_tmp.globalIndex(), muon_tmp.matchMCHTrackId(), muon_tmp.matchMFTTrackId());
        auto mchtrack = muon_tmp.template matchMCHTrack_as<TFwdTracks>(); // MCH-MID
        auto mfttrack = muon_tmp.template matchMFTTrack_as<TMFTTracks>(); // MFTsa

        o2::dataformats::GlobalFwdTrack propmuonAtPV = propagateMuon(muon_tmp, muon_tmp, collision, propagationPoint::kToVertex, matchingZ, mBz, mZShift);
        float eta = propmuonAtPV.getEta();
        float phi = propmuonAtPV.getPhi();
        o2::math_utils::bringTo02Pi(phi);

        float deta = etaMatchedMCHMID - eta;
        float dphi = phiMatchedMCHMID - phi;
        o2::math_utils::bringToPMPi(dphi);
        float dr = std::sqrt(deta * deta + dphi * dphi);
        int ndf = 2 * (mchtrack.nClusters() + mfttrack.nClusters()) - 5;

        if (cfgApplyPreselectionInBestMatch && cfgSlope_dr_chi2MatchMFTMCH * muon_tmp.chi2MatchMCHMFT() + cfgIntercept_dr_chi2MatchMFTMCH < dr) {
          continue;
        }

        fRegistry.fill(HIST("MFTMCHMID/hdR_Chi2MatchMCHMFT"), muon_tmp.chi2MatchMCHMFT(), dr);
        fRegistry.fill(HIST("MFTMCHMID/hdR_Chi2"), muon_tmp.chi2() / ndf, dr);
        fRegistry.fill(HIST("MFTMCHMID/hChi2_Chi2MatchMCHMFT"), muon_tmp.chi2MatchMCHMFT(), muon_tmp.chi2() / ndf);

        // LOGF(info, "muon_tmp.globalIndex() = %d, muon_tmp.matchMCHTrackId() = %d, muon_tmp.matchMFTTrackId() = %d, muon_tmp.chi2MatchMCHMFT() = %f", muon_tmp.globalIndex(), muon_tmp.matchMCHTrackId(), muon_tmp.matchMFTTrackId(), muon_tmp.chi2MatchMCHMFT());

        if (0.f < muon_tmp.chi2MatchMCHMFT() && std::sqrt(std::pow(muon_tmp.chi2MatchMCHMFT() - cfgPeakPosition_chi2MatchMFTMCH, 2)) < min_chi2MatchMCHMFT) {
          min_chi2MatchMCHMFT = std::sqrt(std::pow(muon_tmp.chi2MatchMCHMFT() - cfgPeakPosition_chi2MatchMFTMCH, 2));
          tupleIds_at_min_chi2mftmch = tupleId;
        }

        // if (std::sqrt(std::pow(dr - cfgPeakPosition_dr, 2)) < min_dr) {
        //   min_dr = std::sqrt(std::pow(dr - cfgPeakPosition_dr, 2));
        //   tupleIds_at_min_dr = tupleId;
        // }

        // float distance_2d = std::sqrt(std::pow(muon_tmp.chi2MatchMCHMFT() - cfgPeakPosition_chi2MatchMFTMCH, 2) + std::pow(dr - cfgPeakPosition_dr, 2));
        // if (distance_2d < min_distance_2d) {
        //   min_distance_2d = distance_2d;
        //   tupleIds_at_min_distance_2d = tupleId;
        // }
      }
    }
    vec_min_chi2MatchMCHMFT.emplace_back(tupleIds_at_min_chi2mftmch);
    // vec_min_dr.emplace_back(tupleIds_at_min_dr);
    // vec_min_2d.emplace_back(tupleIds_at_min_distance_2d);

    // LOGF(info, "min: muon_tmp.globalIndex() = %d, muon_tmp.matchMCHTrackId() = %d, muon_tmp.matchMFTTrackId() = %d, muon_tmp.chi2MatchMCHMFT() = %f", std::get<0>(tupleIds_at_min), std::get<1>(tupleIds_at_min), std::get<2>(tupleIds_at_min), min_chi2MatchMCHMFT);
  }

  SliceCache cache;
  Preslice<aod::FwdTracks> perCollision = o2::aod::fwdtrack::collisionId;
  Preslice<aod::FwdTrackAssoc> fwdtrackIndicesPerCollision = aod::track_association::collisionId;
  PresliceUnsorted<aod::FwdTrackAssoc> fwdtrackIndicesPerFwdTrack = aod::track_association::fwdtrackId;
  PresliceUnsorted<aod::FwdTracks> fwdtracksPerMCHTrack = aod::fwdtrack::matchMCHTrackId;
  std::unordered_multimap<int, int> multiMapSAMuonsPerCollision; // collisionId -> trackIds
  std::unordered_multimap<int, int> multiMapGLMuonsPerCollision; // collisionId -> trackIds

  void processRec_SA(MyCollisions const& collisions, MyFwdTracks const& fwdtracks, aod::MFTTracks const& mfttracks, aod::BCsWithTimestamps const&)
  {
    vec_min_chi2MatchMCHMFT.reserve(fwdtracks.size());
    // vec_min_dr.reserve(fwdtracks.size());
    // vec_min_2d.reserve(fwdtracks.size());

    for (const auto& collision : collisions) {
      auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);
      const auto& fwdtracks_per_coll = fwdtracks.sliceBy(perCollision, collision.globalIndex());
      for (const auto& fwdtrack : fwdtracks_per_coll) {
        findBestMatchPerMCHMID<false>(collision, fwdtrack, fwdtracks, mfttracks);
      } // end of fwdtrack loop
    } // end of collision loop

    for (const auto& collision : collisions) {
      const auto& bc = collision.template bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);

      if (!collision.isSelected()) {
        continue;
      }

      const auto& fwdtracks_per_coll = fwdtracks.sliceBy(perCollision, collision.globalIndex());
      for (const auto& fwdtrack : fwdtracks_per_coll) {
        if (fwdtrack.trackType() != o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack && fwdtrack.trackType() != o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack) {
          continue;
        }

        if (fwdtrack.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack && std::find(vec_min_chi2MatchMCHMFT.begin(), vec_min_chi2MatchMCHMFT.end(), std::make_tuple(fwdtrack.globalIndex(), fwdtrack.matchMCHTrackId(), fwdtrack.matchMFTTrackId())) == vec_min_chi2MatchMCHMFT.end()) {
          continue;
        }

        if (!fillFwdTrackTable<false, false, MyFwdTracks, aod::MFTTracks, false>(collision, fwdtrack, nullptr, false)) {
          continue;
        }

        if (fwdtrack.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack) {
          multiMapGLMuonsPerCollision.insert(std::make_pair(collision.globalIndex(), fwdtrack.globalIndex()));
        }
        if (fwdtrack.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack) {
          multiMapSAMuonsPerCollision.insert(std::make_pair(collision.globalIndex(), fwdtrack.globalIndex()));
        }

      } // end of fwdtrack loop
    } // end of collision loop

    for (const auto& collision : collisions) {
      int count_samuons = multiMapSAMuonsPerCollision.count(collision.globalIndex());
      int count_glmuons = multiMapGLMuonsPerCollision.count(collision.globalIndex());
      if (fillQAHistograms) {
        fRegistry.fill(HIST("MCHMID/hNmu"), count_samuons);
        fRegistry.fill(HIST("MFTMCHMID/hNmu"), count_glmuons);
      }
      if (count_samuons >= minNmuon) {
        auto range_samuons = multiMapSAMuonsPerCollision.equal_range(collision.globalIndex());
        for (auto it = range_samuons.first; it != range_samuons.second; it++) {
          auto fwdtrack = fwdtracks.rawIteratorAt(it->second);
          fillFwdTrackTable<false, false, MyFwdTracks, aod::MFTTracks, true>(collision, fwdtrack, nullptr, false);
        }
      }
      if (count_glmuons >= minNmuon) {
        auto range_glmuons = multiMapGLMuonsPerCollision.equal_range(collision.globalIndex());
        for (auto it = range_glmuons.first; it != range_glmuons.second; it++) {
          auto fwdtrack = fwdtracks.rawIteratorAt(it->second);
          fillFwdTrackTable<false, false, MyFwdTracks, aod::MFTTracks, true>(collision, fwdtrack, nullptr, false);
        }
      }
    } // end of collision loop

    multiMapSAMuonsPerCollision.clear();
    multiMapGLMuonsPerCollision.clear();
    map_mfttrackcovs.clear();
    vec_min_chi2MatchMCHMFT.clear();
    vec_min_chi2MatchMCHMFT.shrink_to_fit();
    // vec_min_dr.clear();
    // vec_min_dr.shrink_to_fit();
    // vec_min_2d.clear();
    // vec_min_2d.shrink_to_fit();
  }
  PROCESS_SWITCH(skimmerPrimaryMuon, processRec_SA, "process reconstructed info", false);

  void processRec_TTCA(MyCollisions const& collisions, MyFwdTracks const& fwdtracks, aod::MFTTracks const& mfttracks, aod::BCsWithTimestamps const&, aod::FwdTrackAssoc const& fwdtrackIndices)
  {
    vec_min_chi2MatchMCHMFT.reserve(fwdtracks.size());
    // vec_min_dr.reserve(fwdtracks.size());
    // vec_min_2d.reserve(fwdtracks.size());
    for (const auto& collision : collisions) {
      auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);
      auto fwdtrackIdsThisCollision = fwdtrackIndices.sliceBy(fwdtrackIndicesPerCollision, collision.globalIndex());
      for (const auto& fwdtrackId : fwdtrackIdsThisCollision) {
        auto fwdtrack = fwdtrackId.template fwdtrack_as<MyFwdTracks>();
        findBestMatchPerMCHMID<false>(collision, fwdtrack, fwdtracks, mfttracks);
      } // end of fwdtrack loop
    } // end of collision loop

    std::unordered_map<int64_t, bool> mapAmb; // fwdtrack.globalIndex() -> bool isAmb;
    for (const auto& fwdtrack : fwdtracks) {
      const auto& fwdtrackIdsPerFwdTrack = fwdtrackIndices.sliceBy(fwdtrackIndicesPerFwdTrack, fwdtrack.globalIndex());
      mapAmb[fwdtrack.globalIndex()] = fwdtrackIdsPerFwdTrack.size() > 1;
      // LOGF(info, "fwdtrack.globalIndex() = %d, ntimes = %d, isAmbiguous = %d", fwdtrack.globalIndex(), fwdtrackIdsPerFwdTrack.size(), mapAmb[fwdtrack.globalIndex()]);
    } // end of fwdtrack loop

    for (const auto& collision : collisions) {
      const auto& bc = collision.template bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);

      if (!collision.isSelected()) {
        continue;
      }

      const auto& fwdtrackIdsThisCollision = fwdtrackIndices.sliceBy(fwdtrackIndicesPerCollision, collision.globalIndex());
      for (const auto& fwdtrackId : fwdtrackIdsThisCollision) {
        const auto& fwdtrack = fwdtrackId.template fwdtrack_as<MyFwdTracks>();
        if (fwdtrack.trackType() != o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack && fwdtrack.trackType() != o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack) {
          continue;
        }
        if (fwdtrack.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack && std::find(vec_min_chi2MatchMCHMFT.begin(), vec_min_chi2MatchMCHMFT.end(), std::make_tuple(fwdtrack.globalIndex(), fwdtrack.matchMCHTrackId(), fwdtrack.matchMFTTrackId())) == vec_min_chi2MatchMCHMFT.end()) {
          continue;
        }

        if (!fillFwdTrackTable<false, false, MyFwdTracks, aod::MFTTracks, false>(collision, fwdtrack, nullptr, mapAmb[fwdtrack.globalIndex()])) {
          continue;
        }

        if (fwdtrack.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack) {
          multiMapGLMuonsPerCollision.insert(std::make_pair(collision.globalIndex(), fwdtrack.globalIndex()));
        }
        if (fwdtrack.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack) {
          multiMapSAMuonsPerCollision.insert(std::make_pair(collision.globalIndex(), fwdtrack.globalIndex()));
        }

      } // end of fwdtrack loop
    } // end of collision loop

    for (const auto& collision : collisions) {
      int count_samuons = multiMapSAMuonsPerCollision.count(collision.globalIndex());
      int count_glmuons = multiMapGLMuonsPerCollision.count(collision.globalIndex());
      if (fillQAHistograms) {
        fRegistry.fill(HIST("MCHMID/hNmu"), count_samuons);
        fRegistry.fill(HIST("MFTMCHMID/hNmu"), count_glmuons);
      }
      if (count_samuons >= minNmuon) {
        auto range_samuons = multiMapSAMuonsPerCollision.equal_range(collision.globalIndex());
        for (auto it = range_samuons.first; it != range_samuons.second; it++) {
          auto fwdtrack = fwdtracks.rawIteratorAt(it->second);
          fillFwdTrackTable<false, false, MyFwdTracks, aod::MFTTracks, true>(collision, fwdtrack, nullptr, mapAmb[fwdtrack.globalIndex()]);
        }
      }
      if (count_glmuons >= minNmuon) {
        auto range_glmuons = multiMapGLMuonsPerCollision.equal_range(collision.globalIndex());
        for (auto it = range_glmuons.first; it != range_glmuons.second; it++) {
          auto fwdtrack = fwdtracks.rawIteratorAt(it->second);
          fillFwdTrackTable<false, false, MyFwdTracks, aod::MFTTracks, true>(collision, fwdtrack, nullptr, mapAmb[fwdtrack.globalIndex()]);
        }
      }
    } // end of collision loop

    multiMapSAMuonsPerCollision.clear();
    multiMapGLMuonsPerCollision.clear();
    mapAmb.clear();
    map_mfttrackcovs.clear();
    vec_min_chi2MatchMCHMFT.clear();
    vec_min_chi2MatchMCHMFT.shrink_to_fit();
    // vec_min_dr.clear();
    // vec_min_dr.shrink_to_fit();
    // vec_min_2d.clear();
    // vec_min_2d.shrink_to_fit();
  }
  PROCESS_SWITCH(skimmerPrimaryMuon, processRec_TTCA, "process reconstructed info", false);

  void processRec_TTCA_withMFTCov(MyCollisions const& collisions, MyFwdTracks const& fwdtracks, aod::MFTTracks const& mfttracks, aod::BCsWithTimestamps const&, aod::FwdTrackAssoc const& fwdtrackIndices, aod::MFTTracksCov const& mftCovs)
  {
    for (const auto& mfttrackConv : mftCovs) {
      map_mfttrackcovs[mfttrackConv.matchMFTTrackId()] = mfttrackConv.globalIndex();
    }

    vec_min_chi2MatchMCHMFT.reserve(fwdtracks.size());
    // vec_min_dr.reserve(fwdtracks.size());
    // vec_min_2d.reserve(fwdtracks.size());
    for (const auto& collision : collisions) {
      auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);
      auto fwdtrackIdsThisCollision = fwdtrackIndices.sliceBy(fwdtrackIndicesPerCollision, collision.globalIndex());
      for (const auto& fwdtrackId : fwdtrackIdsThisCollision) {
        auto fwdtrack = fwdtrackId.template fwdtrack_as<MyFwdTracks>();
        findBestMatchPerMCHMID<false>(collision, fwdtrack, fwdtracks, mfttracks);
      } // end of fwdtrack loop
    } // end of collision loop

    std::unordered_map<int64_t, bool> mapAmb; // fwdtrack.globalIndex() -> bool isAmb;
    for (const auto& fwdtrack : fwdtracks) {
      const auto& fwdtrackIdsPerFwdTrack = fwdtrackIndices.sliceBy(fwdtrackIndicesPerFwdTrack, fwdtrack.globalIndex());
      mapAmb[fwdtrack.globalIndex()] = fwdtrackIdsPerFwdTrack.size() > 1;
      // LOGF(info, "fwdtrack.globalIndex() = %d, ntimes = %d, isAmbiguous = %d", fwdtrack.globalIndex(), fwdtrackIdsPerFwdTrack.size(), mapAmb[fwdtrack.globalIndex()]);
    } // end of fwdtrack loop

    for (const auto& collision : collisions) {
      const auto& bc = collision.template bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);

      if (!collision.isSelected()) {
        continue;
      }

      const auto& fwdtrackIdsThisCollision = fwdtrackIndices.sliceBy(fwdtrackIndicesPerCollision, collision.globalIndex());
      for (const auto& fwdtrackId : fwdtrackIdsThisCollision) {
        const auto& fwdtrack = fwdtrackId.template fwdtrack_as<MyFwdTracks>();
        if (fwdtrack.trackType() != o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack && fwdtrack.trackType() != o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack) {
          continue;
        }
        if (fwdtrack.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack && std::find(vec_min_chi2MatchMCHMFT.begin(), vec_min_chi2MatchMCHMFT.end(), std::make_tuple(fwdtrack.globalIndex(), fwdtrack.matchMCHTrackId(), fwdtrack.matchMFTTrackId())) == vec_min_chi2MatchMCHMFT.end()) {
          continue;
        }

        if (!fillFwdTrackTable<false, true, MyFwdTracks, aod::MFTTracks, false>(collision, fwdtrack, mftCovs, mapAmb[fwdtrack.globalIndex()])) {
          continue;
        }

        if (fwdtrack.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack) {
          multiMapGLMuonsPerCollision.insert(std::make_pair(collision.globalIndex(), fwdtrack.globalIndex()));
        }
        if (fwdtrack.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack) {
          multiMapSAMuonsPerCollision.insert(std::make_pair(collision.globalIndex(), fwdtrack.globalIndex()));
        }

      } // end of fwdtrack loop
    } // end of collision loop

    for (const auto& collision : collisions) {
      int count_samuons = multiMapSAMuonsPerCollision.count(collision.globalIndex());
      int count_glmuons = multiMapGLMuonsPerCollision.count(collision.globalIndex());
      if (fillQAHistograms) {
        fRegistry.fill(HIST("MCHMID/hNmu"), count_samuons);
        fRegistry.fill(HIST("MFTMCHMID/hNmu"), count_glmuons);
      }
      if (count_samuons >= minNmuon) {
        auto range_samuons = multiMapSAMuonsPerCollision.equal_range(collision.globalIndex());
        for (auto it = range_samuons.first; it != range_samuons.second; it++) {
          auto fwdtrack = fwdtracks.rawIteratorAt(it->second);
          fillFwdTrackTable<false, true, MyFwdTracks, aod::MFTTracks, true>(collision, fwdtrack, mftCovs, mapAmb[fwdtrack.globalIndex()]);
        }
      }
      if (count_glmuons >= minNmuon) {
        auto range_glmuons = multiMapGLMuonsPerCollision.equal_range(collision.globalIndex());
        for (auto it = range_glmuons.first; it != range_glmuons.second; it++) {
          auto fwdtrack = fwdtracks.rawIteratorAt(it->second);
          fillFwdTrackTable<false, true, MyFwdTracks, aod::MFTTracks, true>(collision, fwdtrack, mftCovs, mapAmb[fwdtrack.globalIndex()]);
        }
      }
    } // end of collision loop

    multiMapSAMuonsPerCollision.clear();
    multiMapGLMuonsPerCollision.clear();
    mapAmb.clear();
    map_mfttrackcovs.clear();
    vec_min_chi2MatchMCHMFT.clear();
    vec_min_chi2MatchMCHMFT.shrink_to_fit();
    // vec_min_dr.clear();
    // vec_min_dr.shrink_to_fit();
    // vec_min_2d.clear();
    // vec_min_2d.shrink_to_fit();
  }
  PROCESS_SWITCH(skimmerPrimaryMuon, processRec_TTCA_withMFTCov, "process reconstructed info", false);

  void processRec_SA_SWT(MyCollisionsWithSWT const& collisions, MyFwdTracks const& fwdtracks, aod::MFTTracks const& mfttracks, aod::BCsWithTimestamps const&)
  {
    vec_min_chi2MatchMCHMFT.reserve(fwdtracks.size());
    // vec_min_dr.reserve(fwdtracks.size());
    // vec_min_2d.reserve(fwdtracks.size());
    for (const auto& collision : collisions) {
      auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);
      const auto& fwdtracks_per_coll = fwdtracks.sliceBy(perCollision, collision.globalIndex());
      for (const auto& fwdtrack : fwdtracks_per_coll) {
        findBestMatchPerMCHMID<false>(collision, fwdtrack, fwdtracks, mfttracks);
      } // end of fwdtrack loop
    } // end of collision loop

    for (const auto& collision : collisions) {
      const auto& bc = collision.template bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);

      if (!collision.isSelected()) {
        continue;
      }

      if (collision.swtaliastmp_raw() == 0) {
        continue;
      }

      const auto& fwdtracks_per_coll = fwdtracks.sliceBy(perCollision, collision.globalIndex());
      for (const auto& fwdtrack : fwdtracks_per_coll) {
        if (fwdtrack.trackType() != o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack && fwdtrack.trackType() != o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack) {
          continue;
        }
        if (fwdtrack.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack && std::find(vec_min_chi2MatchMCHMFT.begin(), vec_min_chi2MatchMCHMFT.end(), std::make_tuple(fwdtrack.globalIndex(), fwdtrack.matchMCHTrackId(), fwdtrack.matchMFTTrackId())) == vec_min_chi2MatchMCHMFT.end()) {
          continue;
        }

        if (!fillFwdTrackTable<false, false, MyFwdTracks, aod::MFTTracks, false>(collision, fwdtrack, nullptr, false)) {
          continue;
        }

        if (fwdtrack.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack) {
          multiMapGLMuonsPerCollision.insert(std::make_pair(collision.globalIndex(), fwdtrack.globalIndex()));
        }
        if (fwdtrack.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack) {
          multiMapSAMuonsPerCollision.insert(std::make_pair(collision.globalIndex(), fwdtrack.globalIndex()));
        }

      } // end of fwdtrack loop
    } // end of collision loop

    for (const auto& collision : collisions) {
      int count_samuons = multiMapSAMuonsPerCollision.count(collision.globalIndex());
      int count_glmuons = multiMapGLMuonsPerCollision.count(collision.globalIndex());
      if (fillQAHistograms) {
        fRegistry.fill(HIST("MCHMID/hNmu"), count_samuons);
        fRegistry.fill(HIST("MFTMCHMID/hNmu"), count_glmuons);
      }
      if (count_samuons >= minNmuon) {
        auto range_samuons = multiMapSAMuonsPerCollision.equal_range(collision.globalIndex());
        for (auto it = range_samuons.first; it != range_samuons.second; it++) {
          auto fwdtrack = fwdtracks.rawIteratorAt(it->second);
          fillFwdTrackTable<false, false, MyFwdTracks, aod::MFTTracks, true>(collision, fwdtrack, nullptr, false);
        }
      }
      if (count_glmuons >= minNmuon) {
        auto range_glmuons = multiMapGLMuonsPerCollision.equal_range(collision.globalIndex());
        for (auto it = range_glmuons.first; it != range_glmuons.second; it++) {
          auto fwdtrack = fwdtracks.rawIteratorAt(it->second);
          fillFwdTrackTable<false, false, MyFwdTracks, aod::MFTTracks, true>(collision, fwdtrack, nullptr, false);
        }
      }
    } // end of collision loop

    multiMapSAMuonsPerCollision.clear();
    multiMapGLMuonsPerCollision.clear();
    map_mfttrackcovs.clear();
    vec_min_chi2MatchMCHMFT.clear();
    vec_min_chi2MatchMCHMFT.shrink_to_fit();
    // vec_min_dr.clear();
    // vec_min_dr.shrink_to_fit();
    // vec_min_2d.clear();
    // vec_min_2d.shrink_to_fit();
  }
  PROCESS_SWITCH(skimmerPrimaryMuon, processRec_SA_SWT, "process reconstructed info only with standalone", false);

  void processRec_TTCA_SWT(MyCollisionsWithSWT const& collisions, MyFwdTracks const& fwdtracks, aod::MFTTracks const& mfttracks, aod::BCsWithTimestamps const&, aod::FwdTrackAssoc const& fwdtrackIndices)
  {
    vec_min_chi2MatchMCHMFT.reserve(fwdtracks.size());
    // vec_min_dr.reserve(fwdtracks.size());
    // vec_min_2d.reserve(fwdtracks.size());
    for (const auto& collision : collisions) {
      auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);
      auto fwdtrackIdsThisCollision = fwdtrackIndices.sliceBy(fwdtrackIndicesPerCollision, collision.globalIndex());
      for (const auto& fwdtrackId : fwdtrackIdsThisCollision) {
        auto fwdtrack = fwdtrackId.template fwdtrack_as<MyFwdTracks>();
        findBestMatchPerMCHMID<false>(collision, fwdtrack, fwdtracks, mfttracks);
      } // end of fwdtrack loop
    } // end of collision loop

    std::unordered_map<int64_t, bool> mapAmb; // fwdtrack.globalIndex() -> bool isAmb;
    for (const auto& fwdtrack : fwdtracks) {
      const auto& fwdtrackIdsPerFwdTrack = fwdtrackIndices.sliceBy(fwdtrackIndicesPerFwdTrack, fwdtrack.globalIndex());
      mapAmb[fwdtrack.globalIndex()] = fwdtrackIdsPerFwdTrack.size() > 1;
      // LOGF(info, "fwdtrack.globalIndex() = %d, ntimes = %d, isAmbiguous = %d", fwdtrack.globalIndex(), fwdtrackIdsPerFwdTrack.size(), mapAmb[fwdtrack.globalIndex()]);
    } // end of fwdtrack loop

    for (const auto& collision : collisions) {
      const auto& bc = collision.template bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);
      if (!collision.isSelected()) {
        continue;
      }
      if (collision.swtaliastmp_raw() == 0) {
        continue;
      }

      const auto& fwdtrackIdsThisCollision = fwdtrackIndices.sliceBy(fwdtrackIndicesPerCollision, collision.globalIndex());
      for (const auto& fwdtrackId : fwdtrackIdsThisCollision) {
        const auto& fwdtrack = fwdtrackId.template fwdtrack_as<MyFwdTracks>();
        if (fwdtrack.trackType() != o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack && fwdtrack.trackType() != o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack) {
          continue;
        }
        if (fwdtrack.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack && std::find(vec_min_chi2MatchMCHMFT.begin(), vec_min_chi2MatchMCHMFT.end(), std::make_tuple(fwdtrack.globalIndex(), fwdtrack.matchMCHTrackId(), fwdtrack.matchMFTTrackId())) == vec_min_chi2MatchMCHMFT.end()) {
          continue;
        }

        if (!fillFwdTrackTable<false, false, MyFwdTracks, aod::MFTTracks, false>(collision, fwdtrack, nullptr, mapAmb[fwdtrack.globalIndex()])) {
          continue;
        }

        if (fwdtrack.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack) {
          multiMapGLMuonsPerCollision.insert(std::make_pair(collision.globalIndex(), fwdtrack.globalIndex()));
        }
        if (fwdtrack.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack) {
          multiMapSAMuonsPerCollision.insert(std::make_pair(collision.globalIndex(), fwdtrack.globalIndex()));
        }

      } // end of fwdtrack loop
    } // end of collision loop

    for (const auto& collision : collisions) {
      int count_samuons = multiMapSAMuonsPerCollision.count(collision.globalIndex());
      int count_glmuons = multiMapGLMuonsPerCollision.count(collision.globalIndex());
      if (fillQAHistograms) {
        fRegistry.fill(HIST("MCHMID/hNmu"), count_samuons);
        fRegistry.fill(HIST("MFTMCHMID/hNmu"), count_glmuons);
      }
      if (count_samuons >= minNmuon) {
        auto range_samuons = multiMapSAMuonsPerCollision.equal_range(collision.globalIndex());
        for (auto it = range_samuons.first; it != range_samuons.second; it++) {
          auto fwdtrack = fwdtracks.rawIteratorAt(it->second);
          fillFwdTrackTable<false, false, MyFwdTracks, aod::MFTTracks, true>(collision, fwdtrack, nullptr, mapAmb[fwdtrack.globalIndex()]);
        }
      }
      if (count_glmuons >= minNmuon) {
        auto range_glmuons = multiMapGLMuonsPerCollision.equal_range(collision.globalIndex());
        for (auto it = range_glmuons.first; it != range_glmuons.second; it++) {
          auto fwdtrack = fwdtracks.rawIteratorAt(it->second);
          fillFwdTrackTable<false, false, MyFwdTracks, aod::MFTTracks, true>(collision, fwdtrack, nullptr, mapAmb[fwdtrack.globalIndex()]);
        }
      }
    } // end of collision loop

    multiMapSAMuonsPerCollision.clear();
    multiMapGLMuonsPerCollision.clear();
    mapAmb.clear();
    map_mfttrackcovs.clear();
    vec_min_chi2MatchMCHMFT.clear();
    vec_min_chi2MatchMCHMFT.shrink_to_fit();
    // vec_min_dr.clear();
    // vec_min_dr.shrink_to_fit();
    // vec_min_2d.clear();
    // vec_min_2d.shrink_to_fit();
  }
  PROCESS_SWITCH(skimmerPrimaryMuon, processRec_TTCA_SWT, "process reconstructed info", false);

  void processRec_TTCA_SWT_withMFTCov(MyCollisionsWithSWT const& collisions, MyFwdTracks const& fwdtracks, aod::MFTTracks const& mfttracks, aod::BCsWithTimestamps const&, aod::FwdTrackAssoc const& fwdtrackIndices, aod::MFTTracksCov const& mftCovs)
  {
    for (const auto& mfttrackConv : mftCovs) {
      map_mfttrackcovs[mfttrackConv.matchMFTTrackId()] = mfttrackConv.globalIndex();
    }
    vec_min_chi2MatchMCHMFT.reserve(fwdtracks.size());
    // vec_min_dr.reserve(fwdtracks.size());
    // vec_min_2d.reserve(fwdtracks.size());
    for (const auto& collision : collisions) {
      auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);
      auto fwdtrackIdsThisCollision = fwdtrackIndices.sliceBy(fwdtrackIndicesPerCollision, collision.globalIndex());
      for (const auto& fwdtrackId : fwdtrackIdsThisCollision) {
        auto fwdtrack = fwdtrackId.template fwdtrack_as<MyFwdTracks>();
        findBestMatchPerMCHMID<false>(collision, fwdtrack, fwdtracks, mfttracks);
      } // end of fwdtrack loop
    } // end of collision loop

    std::unordered_map<int64_t, bool> mapAmb; // fwdtrack.globalIndex() -> bool isAmb;
    for (const auto& fwdtrack : fwdtracks) {
      auto fwdtrackIdsPerFwdTrack = fwdtrackIndices.sliceBy(fwdtrackIndicesPerFwdTrack, fwdtrack.globalIndex());
      mapAmb[fwdtrack.globalIndex()] = fwdtrackIdsPerFwdTrack.size() > 1;
      // LOGF(info, "fwdtrack.globalIndex() = %d, ntimes = %d, isAmbiguous = %d", fwdtrack.globalIndex(), fwdtrackIdsPerFwdTrack.size(), mapAmb[fwdtrack.globalIndex()]);
    } // end of fwdtrack loop

    for (const auto& collision : collisions) {
      auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);
      if (!collision.isSelected()) {
        continue;
      }
      if (collision.swtaliastmp_raw() == 0) {
        continue;
      }

      auto fwdtrackIdsThisCollision = fwdtrackIndices.sliceBy(fwdtrackIndicesPerCollision, collision.globalIndex());
      for (const auto& fwdtrackId : fwdtrackIdsThisCollision) {
        auto fwdtrack = fwdtrackId.template fwdtrack_as<MyFwdTracks>();
        if (fwdtrack.trackType() != o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack && fwdtrack.trackType() != o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack) {
          continue;
        }
        if (fwdtrack.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack && std::find(vec_min_chi2MatchMCHMFT.begin(), vec_min_chi2MatchMCHMFT.end(), std::make_tuple(fwdtrack.globalIndex(), fwdtrack.matchMCHTrackId(), fwdtrack.matchMFTTrackId())) == vec_min_chi2MatchMCHMFT.end()) {
          continue;
        }

        if (!fillFwdTrackTable<false, true, MyFwdTracks, aod::MFTTracks, false>(collision, fwdtrack, mftCovs, mapAmb[fwdtrack.globalIndex()])) {
          continue;
        }

        if (fwdtrack.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack) {
          multiMapGLMuonsPerCollision.insert(std::make_pair(collision.globalIndex(), fwdtrack.globalIndex()));
        }
        if (fwdtrack.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack) {
          multiMapSAMuonsPerCollision.insert(std::make_pair(collision.globalIndex(), fwdtrack.globalIndex()));
        }

      } // end of fwdtrack loop
    } // end of collision loop

    for (const auto& collision : collisions) {
      int count_samuons = multiMapSAMuonsPerCollision.count(collision.globalIndex());
      int count_glmuons = multiMapGLMuonsPerCollision.count(collision.globalIndex());
      if (fillQAHistograms) {
        fRegistry.fill(HIST("MCHMID/hNmu"), count_samuons);
        fRegistry.fill(HIST("MFTMCHMID/hNmu"), count_glmuons);
      }
      if (count_samuons >= minNmuon) {
        auto range_samuons = multiMapSAMuonsPerCollision.equal_range(collision.globalIndex());
        for (auto it = range_samuons.first; it != range_samuons.second; it++) {
          auto fwdtrack = fwdtracks.rawIteratorAt(it->second);
          fillFwdTrackTable<false, true, MyFwdTracks, aod::MFTTracks, true>(collision, fwdtrack, mftCovs, mapAmb[fwdtrack.globalIndex()]);
        }
      }
      if (count_glmuons >= minNmuon) {
        auto range_glmuons = multiMapGLMuonsPerCollision.equal_range(collision.globalIndex());
        for (auto it = range_glmuons.first; it != range_glmuons.second; it++) {
          auto fwdtrack = fwdtracks.rawIteratorAt(it->second);
          fillFwdTrackTable<false, true, MyFwdTracks, aod::MFTTracks, true>(collision, fwdtrack, mftCovs, mapAmb[fwdtrack.globalIndex()]);
        }
      }
    } // end of collision loop

    multiMapSAMuonsPerCollision.clear();
    multiMapGLMuonsPerCollision.clear();
    mapAmb.clear();
    map_mfttrackcovs.clear();
    vec_min_chi2MatchMCHMFT.clear();
    vec_min_chi2MatchMCHMFT.shrink_to_fit();
    // vec_min_dr.clear();
    // vec_min_dr.shrink_to_fit();
    // vec_min_2d.clear();
    // vec_min_2d.shrink_to_fit();
  }
  PROCESS_SWITCH(skimmerPrimaryMuon, processRec_TTCA_SWT_withMFTCov, "process reconstructed info", false);

  void processMC_SA(soa::Join<MyCollisions, aod::McCollisionLabels> const& collisions, MyFwdTracksMC const& fwdtracks, MFTTracksMC const& mfttracks, aod::BCsWithTimestamps const&, aod::McParticles const&)
  {
    vec_min_chi2MatchMCHMFT.reserve(fwdtracks.size());
    // vec_min_dr.reserve(fwdtracks.size());
    // vec_min_2d.reserve(fwdtracks.size());
    for (const auto& collision : collisions) {
      auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);
      auto fwdtracks_per_coll = fwdtracks.sliceBy(perCollision, collision.globalIndex());
      for (const auto& fwdtrack : fwdtracks_per_coll) {
        findBestMatchPerMCHMID<true>(collision, fwdtrack, fwdtracks, mfttracks);
      } // end of fwdtrack loop
    } // end of collision loop

    for (const auto& collision : collisions) {
      auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);
      if (!collision.isSelected()) {
        continue;
      }
      if (!collision.has_mcCollision()) {
        continue;
      }

      auto fwdtracks_per_coll = fwdtracks.sliceBy(perCollision, collision.globalIndex());
      for (const auto& fwdtrack : fwdtracks_per_coll) {
        if (!fwdtrack.has_mcParticle()) {
          continue;
        }
        if (fwdtrack.trackType() != o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack && fwdtrack.trackType() != o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack) {
          continue;
        }
        if (fwdtrack.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack && std::find(vec_min_chi2MatchMCHMFT.begin(), vec_min_chi2MatchMCHMFT.end(), std::make_tuple(fwdtrack.globalIndex(), fwdtrack.matchMCHTrackId(), fwdtrack.matchMFTTrackId())) == vec_min_chi2MatchMCHMFT.end()) {
          continue;
        }

        if (!fillFwdTrackTable<true, false, MyFwdTracksMC, MFTTracksMC, false>(collision, fwdtrack, nullptr, false)) {
          continue;
        }

        if (fwdtrack.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack) {
          multiMapGLMuonsPerCollision.insert(std::make_pair(collision.globalIndex(), fwdtrack.globalIndex()));
        }
        if (fwdtrack.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack) {
          multiMapSAMuonsPerCollision.insert(std::make_pair(collision.globalIndex(), fwdtrack.globalIndex()));
        }

      } // end of fwdtrack loop
    } // end of collision loop

    for (const auto& collision : collisions) {
      int count_samuons = multiMapSAMuonsPerCollision.count(collision.globalIndex());
      int count_glmuons = multiMapGLMuonsPerCollision.count(collision.globalIndex());
      if (fillQAHistograms) {
        fRegistry.fill(HIST("MCHMID/hNmu"), count_samuons);
        fRegistry.fill(HIST("MFTMCHMID/hNmu"), count_glmuons);
      }
      if (count_samuons >= minNmuon) {
        auto range_samuons = multiMapSAMuonsPerCollision.equal_range(collision.globalIndex());
        for (auto it = range_samuons.first; it != range_samuons.second; it++) {
          auto fwdtrack = fwdtracks.rawIteratorAt(it->second);
          fillFwdTrackTable<false, false, MyFwdTracksMC, MFTTracksMC, true>(collision, fwdtrack, nullptr, false);
        }
      }
      if (count_glmuons >= minNmuon) {
        auto range_glmuons = multiMapGLMuonsPerCollision.equal_range(collision.globalIndex());
        for (auto it = range_glmuons.first; it != range_glmuons.second; it++) {
          auto fwdtrack = fwdtracks.rawIteratorAt(it->second);
          fillFwdTrackTable<false, false, MyFwdTracksMC, MFTTracksMC, true>(collision, fwdtrack, nullptr, false);
        }
      }
    } // end of collision loop

    multiMapSAMuonsPerCollision.clear();
    multiMapGLMuonsPerCollision.clear();
    map_mfttrackcovs.clear();
    vec_min_chi2MatchMCHMFT.clear();
    vec_min_chi2MatchMCHMFT.shrink_to_fit();
    // vec_min_dr.clear();
    // vec_min_dr.shrink_to_fit();
    // vec_min_2d.clear();
    // vec_min_2d.shrink_to_fit();
  }
  PROCESS_SWITCH(skimmerPrimaryMuon, processMC_SA, "process reconstructed and MC info", false);

  void processMC_TTCA(soa::Join<MyCollisions, aod::McCollisionLabels> const& collisions, MyFwdTracksMC const& fwdtracks, MFTTracksMC const& mfttracks, aod::BCsWithTimestamps const&, aod::FwdTrackAssoc const& fwdtrackIndices, aod::McParticles const&)
  {
    vec_min_chi2MatchMCHMFT.reserve(fwdtracks.size());
    // vec_min_dr.reserve(fwdtracks.size());
    // vec_min_2d.reserve(fwdtracks.size());
    for (const auto& collision : collisions) {
      auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);
      auto fwdtrackIdsThisCollision = fwdtrackIndices.sliceBy(fwdtrackIndicesPerCollision, collision.globalIndex());
      for (const auto& fwdtrackId : fwdtrackIdsThisCollision) {
        auto fwdtrack = fwdtrackId.template fwdtrack_as<MyFwdTracksMC>();
        findBestMatchPerMCHMID<true>(collision, fwdtrack, fwdtracks, mfttracks);
      } // end of fwdtrack loop
    } // end of collision loop

    std::unordered_map<int64_t, bool> mapAmb; // fwdtrack.globalIndex() -> bool isAmb;
    for (const auto& fwdtrack : fwdtracks) {
      auto fwdtrackIdsPerFwdTrack = fwdtrackIndices.sliceBy(fwdtrackIndicesPerFwdTrack, fwdtrack.globalIndex());
      mapAmb[fwdtrack.globalIndex()] = fwdtrackIdsPerFwdTrack.size() > 1;
      // LOGF(info, "fwdtrack.globalIndex() = %d, ntimes = %d, isAmbiguous = %d", fwdtrack.globalIndex(), fwdtrackIdsPerFwdTrack.size(), mapAmb[fwdtrack.globalIndex()]);
    } // end of fwdtrack loop

    for (const auto& collision : collisions) {
      auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);
      if (!collision.isSelected()) {
        continue;
      }
      if (!collision.has_mcCollision()) {
        continue;
      }

      auto fwdtrackIdsThisCollision = fwdtrackIndices.sliceBy(fwdtrackIndicesPerCollision, collision.globalIndex());
      for (const auto& fwdtrackId : fwdtrackIdsThisCollision) {
        auto fwdtrack = fwdtrackId.template fwdtrack_as<MyFwdTracksMC>();
        if (!fwdtrack.has_mcParticle()) {
          continue;
        }
        if (fwdtrack.trackType() != o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack && fwdtrack.trackType() != o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack) {
          continue;
        }
        if (fwdtrack.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack && std::find(vec_min_chi2MatchMCHMFT.begin(), vec_min_chi2MatchMCHMFT.end(), std::make_tuple(fwdtrack.globalIndex(), fwdtrack.matchMCHTrackId(), fwdtrack.matchMFTTrackId())) == vec_min_chi2MatchMCHMFT.end()) {
          continue;
        }

        if (!fillFwdTrackTable<true, false, MyFwdTracksMC, MFTTracksMC, false>(collision, fwdtrack, nullptr, mapAmb[fwdtrack.globalIndex()])) {
          continue;
        }

        if (fwdtrack.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack) {
          multiMapGLMuonsPerCollision.insert(std::make_pair(collision.globalIndex(), fwdtrack.globalIndex()));
        }
        if (fwdtrack.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack) {
          multiMapSAMuonsPerCollision.insert(std::make_pair(collision.globalIndex(), fwdtrack.globalIndex()));
        }

      } // end of fwdtrack loop
    } // end of collision loop

    for (const auto& collision : collisions) {
      int count_samuons = multiMapSAMuonsPerCollision.count(collision.globalIndex());
      int count_glmuons = multiMapGLMuonsPerCollision.count(collision.globalIndex());
      if (fillQAHistograms) {
        fRegistry.fill(HIST("MCHMID/hNmu"), count_samuons);
        fRegistry.fill(HIST("MFTMCHMID/hNmu"), count_glmuons);
      }
      if (count_samuons >= minNmuon) {
        auto range_samuons = multiMapSAMuonsPerCollision.equal_range(collision.globalIndex());
        for (auto it = range_samuons.first; it != range_samuons.second; it++) {
          auto fwdtrack = fwdtracks.rawIteratorAt(it->second);
          fillFwdTrackTable<true, false, MyFwdTracksMC, MFTTracksMC, true>(collision, fwdtrack, nullptr, mapAmb[fwdtrack.globalIndex()]);
        }
      }
      if (count_glmuons >= minNmuon) {
        auto range_glmuons = multiMapGLMuonsPerCollision.equal_range(collision.globalIndex());
        for (auto it = range_glmuons.first; it != range_glmuons.second; it++) {
          auto fwdtrack = fwdtracks.rawIteratorAt(it->second);
          fillFwdTrackTable<true, false, MyFwdTracksMC, MFTTracksMC, true>(collision, fwdtrack, nullptr, mapAmb[fwdtrack.globalIndex()]);
        }
      }
    } // end of collision loop

    multiMapSAMuonsPerCollision.clear();
    multiMapGLMuonsPerCollision.clear();
    mapAmb.clear();
    map_mfttrackcovs.clear();
    vec_min_chi2MatchMCHMFT.clear();
    vec_min_chi2MatchMCHMFT.shrink_to_fit();
    // vec_min_dr.clear();
    // vec_min_dr.shrink_to_fit();
    // vec_min_2d.clear();
    // vec_min_2d.shrink_to_fit();
  }
  PROCESS_SWITCH(skimmerPrimaryMuon, processMC_TTCA, "process reconstructed and MC info", false);

  void processMC_TTCA_withMFTCov(soa::Join<MyCollisions, aod::McCollisionLabels> const& collisions, MyFwdTracksMC const& fwdtracks, MFTTracksMC const& mfttracks, aod::BCsWithTimestamps const&, aod::FwdTrackAssoc const& fwdtrackIndices, aod::MFTTracksCov const& mftCovs, aod::McParticles const&)
  {
    for (const auto& mfttrackConv : mftCovs) {
      map_mfttrackcovs[mfttrackConv.matchMFTTrackId()] = mfttrackConv.globalIndex();
    }
    vec_min_chi2MatchMCHMFT.reserve(fwdtracks.size());
    // vec_min_dr.reserve(fwdtracks.size());
    // vec_min_2d.reserve(fwdtracks.size());
    for (const auto& collision : collisions) {
      auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);
      auto fwdtrackIdsThisCollision = fwdtrackIndices.sliceBy(fwdtrackIndicesPerCollision, collision.globalIndex());
      for (const auto& fwdtrackId : fwdtrackIdsThisCollision) {
        auto fwdtrack = fwdtrackId.template fwdtrack_as<MyFwdTracksMC>();
        findBestMatchPerMCHMID<true>(collision, fwdtrack, fwdtracks, mfttracks);
      } // end of fwdtrack loop
    } // end of collision loop

    std::unordered_map<int64_t, bool> mapAmb; // fwdtrack.globalIndex() -> bool isAmb;
    for (const auto& fwdtrack : fwdtracks) {
      auto fwdtrackIdsPerFwdTrack = fwdtrackIndices.sliceBy(fwdtrackIndicesPerFwdTrack, fwdtrack.globalIndex());
      mapAmb[fwdtrack.globalIndex()] = fwdtrackIdsPerFwdTrack.size() > 1;
      // LOGF(info, "fwdtrack.globalIndex() = %d, ntimes = %d, isAmbiguous = %d", fwdtrack.globalIndex(), fwdtrackIdsPerFwdTrack.size(), mapAmb[fwdtrack.globalIndex()]);
    } // end of fwdtrack loop

    for (const auto& collision : collisions) {
      auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);
      if (!collision.isSelected()) {
        continue;
      }
      if (!collision.has_mcCollision()) {
        continue;
      }

      auto fwdtrackIdsThisCollision = fwdtrackIndices.sliceBy(fwdtrackIndicesPerCollision, collision.globalIndex());
      for (const auto& fwdtrackId : fwdtrackIdsThisCollision) {
        auto fwdtrack = fwdtrackId.template fwdtrack_as<MyFwdTracksMC>();
        if (!fwdtrack.has_mcParticle()) {
          continue;
        }
        if (fwdtrack.trackType() != o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack && fwdtrack.trackType() != o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack) {
          continue;
        }
        if (fwdtrack.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack && std::find(vec_min_chi2MatchMCHMFT.begin(), vec_min_chi2MatchMCHMFT.end(), std::make_tuple(fwdtrack.globalIndex(), fwdtrack.matchMCHTrackId(), fwdtrack.matchMFTTrackId())) == vec_min_chi2MatchMCHMFT.end()) {
          continue;
        }

        if (!fillFwdTrackTable<true, true, MyFwdTracksMC, MFTTracksMC, false>(collision, fwdtrack, mftCovs, mapAmb[fwdtrack.globalIndex()])) {
          continue;
        }

        if (fwdtrack.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack) {
          multiMapGLMuonsPerCollision.insert(std::make_pair(collision.globalIndex(), fwdtrack.globalIndex()));
        }
        if (fwdtrack.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack) {
          multiMapSAMuonsPerCollision.insert(std::make_pair(collision.globalIndex(), fwdtrack.globalIndex()));
        }

      } // end of fwdtrack loop
    } // end of collision loop

    for (const auto& collision : collisions) {
      int count_samuons = multiMapSAMuonsPerCollision.count(collision.globalIndex());
      int count_glmuons = multiMapGLMuonsPerCollision.count(collision.globalIndex());
      if (fillQAHistograms) {
        fRegistry.fill(HIST("MCHMID/hNmu"), count_samuons);
        fRegistry.fill(HIST("MFTMCHMID/hNmu"), count_glmuons);
      }
      if (count_samuons >= minNmuon) {
        auto range_samuons = multiMapSAMuonsPerCollision.equal_range(collision.globalIndex());
        for (auto it = range_samuons.first; it != range_samuons.second; it++) {
          auto fwdtrack = fwdtracks.rawIteratorAt(it->second);
          fillFwdTrackTable<true, true, MyFwdTracksMC, MFTTracksMC, true>(collision, fwdtrack, mftCovs, mapAmb[fwdtrack.globalIndex()]);
        }
      }
      if (count_glmuons >= minNmuon) {
        auto range_glmuons = multiMapGLMuonsPerCollision.equal_range(collision.globalIndex());
        for (auto it = range_glmuons.first; it != range_glmuons.second; it++) {
          auto fwdtrack = fwdtracks.rawIteratorAt(it->second);
          fillFwdTrackTable<true, true, MyFwdTracksMC, MFTTracksMC, true>(collision, fwdtrack, mftCovs, mapAmb[fwdtrack.globalIndex()]);
        }
      }
    } // end of collision loop

    multiMapSAMuonsPerCollision.clear();
    multiMapGLMuonsPerCollision.clear();
    mapAmb.clear();
    map_mfttrackcovs.clear();
    vec_min_chi2MatchMCHMFT.clear();
    vec_min_chi2MatchMCHMFT.shrink_to_fit();
    // vec_min_dr.clear();
    // vec_min_dr.shrink_to_fit();
    // vec_min_2d.clear();
    // vec_min_2d.shrink_to_fit();
  }
  PROCESS_SWITCH(skimmerPrimaryMuon, processMC_TTCA_withMFTCov, "process reconstructed and MC with MFTCov info", false);

  void processDummy(aod::Collisions const&) {}
  PROCESS_SWITCH(skimmerPrimaryMuon, processDummy, "process dummy", true);
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
          self_Ids.emplace_back(global_muon.globalIndex());

          // if (global_muon.collisionId() == muon.collisionId()) {
          //   self_Ids.emplace_back(global_muon.globalIndex());
          // }
        }
        em_same_mft_ids(self_Ids);
        self_Ids.clear();
        self_Ids.shrink_to_fit();
      } else {
        em_same_mft_ids(std::vector<int>{}); // empty for standalone muons
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
