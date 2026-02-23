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
// This code will create data table for inputs to machine learning for electrons.
//    Please write to: daiki.sekihata@cern.ch

#include "PWGEM/Dilepton/DataModel/lmeeMLTables.h"

#include "Common/CCDB/RCTSelectionFlags.h"
#include "Common/CCDB/ctpRateFetcher.h"
#include "Common/Core/CollisionTypeHelper.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/fwdtrackUtilities.h"
#include "Common/DataModel/EventSelection.h"

#include "CCDB/BasicCCDBManager.h"
#include "CommonConstants/PhysicsConstants.h"
#include "DCAFitter/DCAFitterN.h"
#include "DataFormatsParameters/GRPLHCIFData.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsBase/Propagator.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "GlobalTracking/MatchGlobalFwd.h"
#include "MCHTracking/TrackExtrap.h"
#include "MCHTracking/TrackParam.h"
#include "ReconstructionDataFormats/TrackFwd.h"

#include <map>
#include <random>
#include <string>
#include <unordered_map>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::constants::physics;
using namespace o2::aod::fwdtrackutils;

struct TreeCreatorMuonML {
  using MyCollisionsMC = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>;
  using MyCollisionMC = MyCollisionsMC::iterator;

  using MyFwdTracksMC = soa::Join<aod::FwdTracks, aod::FwdTracksCov, aod::McFwdTrackLabels>;
  using MyFwdTrackMC = MyFwdTracksMC::iterator;

  using MyMFTTracksMC = soa::Join<o2::aod::MFTTracks, aod::McMFTTrackLabels>;
  using MyMFTTrackMC = MyMFTTracksMC::iterator;

  Produces<aod::EMFwdTracksForML> mltable;

  // Configurables
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
  // Configurable<std::string> irSource{"irSource", "ZNC hadronic", "Estimator of the interaction rate (Recommended: pp/OO --> T0VTX, Pb-Pb --> ZNC hadronic)"};

  // for z shift for propagation
  Configurable<bool> cfgApplyZShiftFromCCDB{"cfgApplyZShiftFromCCDB", false, "flag to apply z shift"};
  Configurable<std::string> cfgZShiftPath{"cfgZShiftPath", "Users/m/mcoquet/ZShift", "CCDB path for z shift to apply to forward tracks"};
  Configurable<float> cfgManualZShift{"cfgManualZShift", 0, "manual z-shift for propagation of global muon to PV"};
  Configurable<float> cfgDownSampling{"cfgDownSampling", 1.1, "down sampling for fake matches"};

  struct : ConfigurableGroup {
    std::string prefix = "eventCutGroup";
    Configurable<float> cfgMinZvtx{"cfgMinZvtx", -10.f, "min. Zvtx of collision"};
    Configurable<float> cfgMaxZvtx{"cfgMaxZvtx", +10.f, "max. Zvtx of collision"};
    Configurable<bool> cfgRequireFT0AND{"cfgRequireFT0AND", true, "require FT0AND"};
    Configurable<bool> cfgRequireNoTFB{"cfgRequireNoTFB", true, "require No time frame border"};
    Configurable<bool> cfgRequireNoITSROFB{"cfgRequireNoITSROFB", false, "require no ITS readout frame border"};
    Configurable<bool> cfgRequireNoSameBunchPileup{"cfgRequireNoSameBunchPileup", false, "require no same bunch pileup in event cut"};
    Configurable<bool> cfgRequireGoodZvtxFT0vsPV{"cfgRequireGoodZvtxFT0vsPV", false, "require good Zvtx between FT0 vs. PV in event cut"};
  } eventCutGroup;

  struct : ConfigurableGroup {
    std::string prefix = "glMuonCutGroup";
    Configurable<float> minEta{"minEta", -3.6, "min. eta acceptance for MFT-MCH-MID"};
    Configurable<float> maxEta{"maxEta", -2.5, "max. eta acceptance for MFT-MCH-MID"};
    Configurable<float> maxMatchingChi2MCHMFT{"maxMatchingChi2MCHMFT", 100.f, "max. chi2 for MCH-MFT matching"};
    Configurable<float> maxChi2{"maxChi2", 20.f, "max. chi2 for global muon"};
    Configurable<bool> refitGlobalMuon{"refitGlobalMuon", true, "flag to refit global muon"};
    Configurable<float> matchingZ{"matchingZ", -77.5, "z position where matching is performed"};
  } glMuonCutGroup;

  o2::ccdb::CcdbApi ccdbApi;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  std::mt19937 engine;
  std::uniform_real_distribution<float> dist01;

  ctpRateFetcher mRateFetcher;
  std::string irSourceForCptFetcher{""};

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

    addHistograms();

    mRunNumber = 0;
    mBz = 0;
    mZShift = 0;

    std::random_device seed_gen;
    engine = std::mt19937(seed_gen());
    dist01 = std::uniform_real_distribution<float>(0.0f, 1.0f);
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

    o2::parameters::GRPLHCIFData* grplhcif = ccdb.service->getSpecificForRun<o2::parameters::GRPLHCIFData>("GLO/Config/GRPLHCIF", bc.runNumber());
    auto collsys = o2::common::core::CollisionSystemType::getCollisionTypeFromGrp(grplhcif);
    if (collsys == o2::common::core::CollisionSystemType::kCollSyspp || collsys == o2::common::core::CollisionSystemType::kCollSysOO) {
      irSourceForCptFetcher = std::string("T0VTX");
    } else {
      irSourceForCptFetcher = std::string("ZNC hadronic");
    }
    LOGF(info, "irSourceForCptFetcher = %s", irSourceForCptFetcher.data());
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

  template <typename TCollision>
  bool isSelectedCollision(TCollision const& collision)
  {
    if (!(eventCutGroup.cfgMinZvtx < collision.posZ() && collision.posZ() < eventCutGroup.cfgMaxZvtx)) {
      return false;
    }
    if (eventCutGroup.cfgRequireFT0AND && !collision.selection_bit(o2::aod::evsel::kIsTriggerTVX)) {
      return false;
    }
    if (eventCutGroup.cfgRequireNoTFB && !collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
      return false;
    }
    if (eventCutGroup.cfgRequireNoITSROFB && !collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
      return false;
    }
    if (eventCutGroup.cfgRequireNoSameBunchPileup && !collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      return false;
    }
    if (eventCutGroup.cfgRequireGoodZvtxFT0vsPV && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return false;
    }
    return true;
  }

  template <bool withMFTCov, typename TCollision, typename TFwdTrack, typename TFwdTracks, typename TMFTTracks, typename TMFTTracksCov>
  bool fillFwdTrackTable(TCollision const& collision, TFwdTrack const& fwdtrack, TFwdTracks const&, TMFTTracks const&, TMFTTracksCov const& mftCovs, const float hadronicRate)
  {
    if (fwdtrack.trackType() != o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack) {
      return false;
    }

    auto mchtrack = fwdtrack.template matchMCHTrack_as<TFwdTracks>(); // MCH-MID
    auto mfttrack = fwdtrack.template matchMFTTrack_as<TMFTTracks>(); // MFTsa
    if (!mfttrack.has_mcParticle() || !mchtrack.has_mcParticle() || !fwdtrack.has_mcParticle()) {
      return false;
    }
    float chi2 = fwdtrack.chi2() / (2.f * (mchtrack.nClusters() + mfttrack.nClusters()) - 5.f);
    float chi2mft = mfttrack.chi2() / (2.f * mfttrack.nClusters() - 5.f);

    // auto mcParticle_MFTMCHMID = fwdtrack.template mcParticle_as<aod::McParticles>(); // this is identical to mcParticle_MCHMID
    auto mcParticle_MCHMID = mchtrack.template mcParticle_as<aod::McParticles>(); // this is identical to mcParticle_MFTMCHMID
    auto mcParticle_MFT = mfttrack.template mcParticle_as<aod::McParticles>();

    int pdgCode = mcParticle_MCHMID.pdgCode();
    bool isPrimary = mcParticle_MCHMID.isPhysicalPrimary() || mcParticle_MCHMID.producedByGenerator();
    bool isMatched = (mcParticle_MFT.globalIndex() == mcParticle_MCHMID.globalIndex()) && (mcParticle_MFT.mcCollisionId() == mcParticle_MCHMID.mcCollisionId());

    if (!isMatched && dist01(engine) > cfgDownSampling) {
      return false;
    }

    if (fwdtrack.chi2MatchMCHMID() < 0.f) { // this should never happen. only for protection.
      return false;
    }
    if (fwdtrack.chi2MatchMCHMFT() < 0.f || glMuonCutGroup.maxMatchingChi2MCHMFT < fwdtrack.chi2MatchMCHMFT()) { // this should never happen. only for protection.
      return false;
    }
    if (fwdtrack.chi2() < 0.f || glMuonCutGroup.maxChi2 < chi2) { // this should never happen. only for protection.
      return false;
    }
    if (mfttrack.chi2() < 0.f) { // this should never happen. only for protection.
      return false;
    }

    o2::dataformats::GlobalFwdTrack propmuonAtPV = propagateMuon(fwdtrack, fwdtrack, collision, propagationPoint::kToVertex, glMuonCutGroup.matchingZ, mBz, mZShift);
    float pt = propmuonAtPV.getPt();
    float eta = propmuonAtPV.getEta();
    float phi = propmuonAtPV.getPhi();
    // o2::math_utils::bringTo02Pi(phi);
    phi = RecoDecay::constrainAngle(phi, 0, 1U);

    if (eta < glMuonCutGroup.minEta || glMuonCutGroup.maxEta < eta) {
      return false;
    }

    float dcaX = propmuonAtPV.getX() - collision.posX();
    float dcaY = propmuonAtPV.getY() - collision.posY();
    // float dcaZ = propmuonAtPV.getZ() - collision.posZ();
    float dcaXY = std::sqrt(dcaX * dcaX + dcaY * dcaY);

    o2::dataformats::GlobalFwdTrack propmuonAtPV_Matched = propagateMuon(mchtrack, mchtrack, collision, propagationPoint::kToVertex, glMuonCutGroup.matchingZ, mBz, mZShift);
    float ptMatchedMCHMID = propmuonAtPV_Matched.getPt();
    float etaMatchedMCHMID = propmuonAtPV_Matched.getEta();
    float phiMatchedMCHMID = propmuonAtPV_Matched.getPhi();
    // o2::math_utils::bringTo02Pi(phiMatchedMCHMID);
    phiMatchedMCHMID = RecoDecay::constrainAngle(phiMatchedMCHMID, 0, 1U);
    if (glMuonCutGroup.refitGlobalMuon) {
      pt = propmuonAtPV_Matched.getP() * std::sin(2.f * std::atan(std::exp(-eta)));
    }

    o2::dataformats::GlobalFwdTrack propmuonAtDCA_Matched = propagateMuon(mchtrack, mchtrack, collision, propagationPoint::kToDCA, glMuonCutGroup.matchingZ, mBz, mZShift);
    float dcaX_Matched = propmuonAtDCA_Matched.getX() - collision.posX();
    float dcaY_Matched = propmuonAtDCA_Matched.getY() - collision.posY();
    float dcaXY_Matched = std::sqrt(dcaX_Matched * dcaX_Matched + dcaY_Matched * dcaY_Matched);
    float pDCA = mchtrack.p() * dcaXY_Matched;
    float rAtAbsorberEnd = fwdtrack.rAtAbsorberEnd(); // this works only for GlobalMuonTrack

    float xMatchedMFTatMP = 999.f;
    float yMatchedMFTatMP = 999.f;
    float xMatchedMCHMIDatMP = 999.f;
    float yMatchedMCHMIDatMP = 999.f;

    if constexpr (withMFTCov) {
      auto mfttrackcov = mftCovs.rawIteratorAt(map_mfttrackcovs[mfttrack.globalIndex()]);
      o2::track::TrackParCovFwd mftsaAtMP = getTrackParCovFwdShift(mfttrack, mZShift, mfttrackcov); // values at innermost update
      mftsaAtMP.propagateToZhelix(glMuonCutGroup.matchingZ, mBz);                                   // propagated to matching plane
      xMatchedMFTatMP = mftsaAtMP.getX();
      yMatchedMFTatMP = mftsaAtMP.getY();

      auto muonAtMP = propagateMuon(mchtrack, mchtrack, collision, propagationPoint::kToMatchingPlane, glMuonCutGroup.matchingZ, mBz, mZShift); // propagated to matching plane
      xMatchedMCHMIDatMP = muonAtMP.getX();
      yMatchedMCHMIDatMP = muonAtMP.getY();
    }

    float deta = etaMatchedMCHMID - eta;
    float dphi = phiMatchedMCHMID - phi;
    o2::math_utils::bringToPMPi(dphi);

    mltable(collision.posZ(), collision.numContrib(), collision.trackOccupancyInTimeRange(), collision.ft0cOccupancyInTimeRange(), hadronicRate,
            fwdtrack.trackType(), pt, eta, phi, fwdtrack.sign(), dcaX, dcaY, ptMatchedMCHMID, etaMatchedMCHMID, phiMatchedMCHMID,
            xMatchedMCHMIDatMP, yMatchedMCHMIDatMP, xMatchedMFTatMP, yMatchedMFTatMP,
            fwdtrack.nClusters(), pDCA, rAtAbsorberEnd, chi2, fwdtrack.chi2MatchMCHMID(), fwdtrack.chi2MatchMCHMFT(),
            // fwdtrack.mchBitMap(), fwdtrack.midBitMap(), fwdtrack.midBoards(),
            mfttrack.mftClusterSizesAndTrackFlags(), chi2mft, mfttrack.nClusters(), pdgCode, isPrimary, isMatched);

    fRegistry.fill(HIST("hMuonType"), fwdtrack.trackType());
    fRegistry.fill(HIST("MFTMCHMID/hPt"), pt);
    fRegistry.fill(HIST("MFTMCHMID/hEtaPhi"), phi, eta);
    fRegistry.fill(HIST("MFTMCHMID/hEtaPhi_MatchedMCHMID"), phiMatchedMCHMID, etaMatchedMCHMID);
    fRegistry.fill(HIST("MFTMCHMID/hDeltaEta_Pt"), pt, deta);
    fRegistry.fill(HIST("MFTMCHMID/hDeltaPhi_Pt"), pt, dphi);
    fRegistry.fill(HIST("MFTMCHMID/hSign"), fwdtrack.sign());
    fRegistry.fill(HIST("MFTMCHMID/hNclusters"), fwdtrack.nClusters());
    fRegistry.fill(HIST("MFTMCHMID/hNclustersMFT"), mfttrack.nClusters());
    fRegistry.fill(HIST("MFTMCHMID/hPDCA_Rabs"), rAtAbsorberEnd, pDCA);
    fRegistry.fill(HIST("MFTMCHMID/hRatAbsorberEnd"), rAtAbsorberEnd);
    fRegistry.fill(HIST("MFTMCHMID/hChi2"), chi2);
    fRegistry.fill(HIST("MFTMCHMID/hChi2MFT"), chi2mft);
    fRegistry.fill(HIST("MFTMCHMID/hChi2MatchMCHMID"), fwdtrack.chi2MatchMCHMID());
    fRegistry.fill(HIST("MFTMCHMID/hChi2MatchMCHMFT"), fwdtrack.chi2MatchMCHMFT());
    fRegistry.fill(HIST("MFTMCHMID/hDCAxy2D"), dcaX, dcaY);
    fRegistry.fill(HIST("MFTMCHMID/hDCAxy"), dcaXY);
    fRegistry.fill(HIST("MFTMCHMID/hDCAx_PosZ"), collision.posZ(), dcaX);
    fRegistry.fill(HIST("MFTMCHMID/hDCAy_PosZ"), collision.posZ(), dcaY);
    fRegistry.fill(HIST("MFTMCHMID/hDCAx_Phi"), phi, dcaX);
    fRegistry.fill(HIST("MFTMCHMID/hDCAy_Phi"), phi, dcaY);
    return true;
  }

  SliceCache cache;
  Preslice<aod::FwdTracks> perCollision = o2::aod::fwdtrack::collisionId;

  std::unordered_map<int, int> map_mfttrackcovs;
  void processWithMFTCov(MyCollisionsMC const& collisions, aod::BCsWithTimestamps const&, MyFwdTracksMC const& fwdtracks, MyMFTTracksMC const& mfttracks, aod::MFTTracksCov const& mftCovs, aod::McParticles const&, aod::McCollisions const&)
  {
    for (const auto& mfttrackConv : mftCovs) {
      map_mfttrackcovs[mfttrackConv.matchMFTTrackId()] = mfttrackConv.globalIndex();
    }

    for (const auto& collision : collisions) {
      auto bc = collision.template foundBC_as<aod::BCsWithTimestamps>();
      initCCDB(bc);

      if (!collision.has_mcCollision()) {
        continue;
      }

      if (!isSelectedCollision(collision)) {
        continue;
      }
      float hadronicRate = mRateFetcher.fetch(ccdb.service, bc.timestamp(), bc.runNumber(), irSourceForCptFetcher) * 1.e-3; // kHz

      auto fwdtracks_coll = fwdtracks.sliceBy(perCollision, collision.globalIndex());
      for (const auto& fwdtrack : fwdtracks_coll) {
        if (!fwdtrack.has_mcParticle()) {
          continue;
        }
        if (fwdtrack.trackType() != o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack) {
          continue;
        }

        fillFwdTrackTable<true>(collision, fwdtrack, fwdtracks, mfttracks, mftCovs, hadronicRate);

      } // end of fwdtrack loop
    } // end of collision loop

    map_mfttrackcovs.clear();
  }
  PROCESS_SWITCH(TreeCreatorMuonML, processWithMFTCov, "produce ML input for single track level", true);

  void processWithoutMFTCov(MyCollisionsMC const& collisions, aod::BCsWithTimestamps const&, MyFwdTracksMC const& fwdtracks, MyMFTTracksMC const& mfttracks, aod::McParticles const&, aod::McCollisions const&)
  {
    for (const auto& collision : collisions) {
      auto bc = collision.template foundBC_as<aod::BCsWithTimestamps>();
      initCCDB(bc);

      if (!collision.has_mcCollision()) {
        continue;
      }

      if (!isSelectedCollision(collision)) {
        continue;
      }
      float hadronicRate = mRateFetcher.fetch(ccdb.service, bc.timestamp(), bc.runNumber(), irSourceForCptFetcher) * 1.e-3; // kHz

      auto fwdtracks_coll = fwdtracks.sliceBy(perCollision, collision.globalIndex());
      for (const auto& fwdtrack : fwdtracks_coll) {
        if (!fwdtrack.has_mcParticle()) {
          continue;
        }
        if (fwdtrack.trackType() != o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack) {
          continue;
        }

        fillFwdTrackTable<false>(collision, fwdtrack, fwdtracks, mfttracks, nullptr, hadronicRate);

      } // end of fwdtrack loop
    } // end of collision loop
  }
  PROCESS_SWITCH(TreeCreatorMuonML, processWithoutMFTCov, "produce ML input for single track level", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<TreeCreatorMuonML>(cfgc, TaskName{"tree-creator-muon-ml"})};
}
