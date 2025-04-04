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

#include <string>
#include <map>
#include <unordered_map>

#include "Framework/DataTypes.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "TableHelper.h"
#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "TGeoGlobalMagField.h"
#include "Field/MagneticField.h"

#include "DetectorsBase/Propagator.h"
#include "GlobalTracking/MatchGlobalFwd.h"
#include "MCHTracking/TrackExtrap.h"
#include "MCHTracking/TrackParam.h"
#include "ReconstructionDataFormats/TrackFwd.h"

#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/DataModel/PropagatedFwdTrackTables.h"
#include "Common/Core/fwdtrackUtilities.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"

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
  Configurable<float> maxChi2SA{"maxChi2SA", 1e+6, "max. chi2 for standalone muon"};
  Configurable<float> maxChi2GL{"maxChi2GL", 1e+6, "max. chi2 for global muon"};
  Configurable<bool> refitGlobalMuon{"refitGlobalMuon", true, "flag to refit global muon"};
  Configurable<bool> applyEtaCutToSAinGL{"applyEtaCutToSAinGL", false, "flag to apply eta cut to samuon in global muon"};
  Configurable<bool> requireTrueAssociation{"requireTrueAssociation", false, "flag to require true mc collision association"};

  Configurable<int> cfgCentEstimator{"cfgCentEstimator", 2, "FT0M:0, FT0A:1, FT0C:2"};
  Configurable<float> cfgCentMin{"cfgCentMin", -1.f, "min. centrality"};
  Configurable<float> cfgCentMax{"cfgCentMax", 999.f, "max. centrality"};
  Configurable<float> cfgZvtxMin{"cfgZvtxMin", -10.f, "min. Zvtx"};
  Configurable<float> cfgZvtxMax{"cfgZvtxMax", 10.f, "max. Zvtx"};

  HistogramRegistry fRegistry{"fRegistry"};
  static constexpr std::string_view muon_types[5] = {"MFTMCHMID/", "MFTMCHMIDOtherMatch/", "MFTMCH/", "MCHMID/", "MCH/"};

  void init(o2::framework::InitContext&)
  {
    if (doprocessWithoutFTTCA && doprocessWithFTTCA) {
      LOGF(fatal, "Cannot enable doprocessWithoutFTTCA and doprocessWithFTTCA at the same time. Please choose one.");
    }

    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);
    ccdbApi.init(ccdburl);

    addHistograms();
  }

  o2::ccdb::CcdbApi ccdbApi;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  int mRunNumber = -1;

  template <typename TBC>
  void initCCDB(TBC const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }
    mRunNumber = bc.runNumber();
    LOGF(info, "mRunNumber = %d", mRunNumber);
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

    fRegistry.add("MFTMCHMID/primary/correct/hPt", "pT;p_{T} (GeV/c)", kTH1F, {{100, 0.0f, 10}}, false);
    fRegistry.add("MFTMCHMID/primary/correct/hEtaPhi", "#eta vs. #varphi;#varphi (rad.);#eta", kTH2F, {{180, 0, 2 * M_PI}, {100, -6.f, -1.f}}, false);
    fRegistry.add("MFTMCHMID/primary/correct/hEtaPhi_MatchedMCHMID", "#eta vs. #varphi;#varphi (rad.);#eta", kTH2F, {{180, 0, 2 * M_PI}, {100, -6.f, -1.f}}, false);
    fRegistry.add("MFTMCHMID/primary/correct/hDeltaEtaDeltaPhi", "#Delta#eta vs. #Delta#varphi;#Delta#varphi = #varphi_{sa} - #varphi_{gl} (rad.);#Delta#eta = #eta_{sa} - #eta_{gl}", kTH2F, {{180, -M_PI, M_PI}, {400, -2.f, 2.f}}, false);
    fRegistry.add("MFTMCHMID/primary/correct/hDiffCollId", "difference in collision index;collisionId_{TTCA} - collisionId_{MP}", kTH1F, {{41, -20.5, +20.5}}, false);
    fRegistry.add("MFTMCHMID/primary/correct/hSign", "sign;sign", kTH1F, {{3, -1.5, +1.5}}, false);
    fRegistry.add("MFTMCHMID/primary/correct/hNclusters", "Nclusters;Nclusters", kTH1F, {{21, -0.5f, 20.5}}, false);
    fRegistry.add("MFTMCHMID/primary/correct/hNclustersMFT", "NclustersMFT;Nclusters MFT", kTH1F, {{11, -0.5f, 10.5}}, false);
    fRegistry.add("MFTMCHMID/primary/correct/hRatAbsorberEnd", "R at absorber end;R at absorber end (cm)", kTH1F, {{100, 0.0f, 100}}, false);
    fRegistry.add("MFTMCHMID/primary/correct/hPDCA_Rabs", "pDCA vs. Rabs;R at absorber end (cm);p #times DCA (GeV/c #upoint cm)", kTH2F, {{100, 0, 100}, {100, 0.0f, 1000}}, false);
    fRegistry.add("MFTMCHMID/primary/correct/hChi2", "chi2;chi2", kTH1F, {{100, 0.0f, 100}}, false);
    fRegistry.add("MFTMCHMID/primary/correct/hChi2MFT", "chi2 MFT;chi2 MFT", kTH1F, {{100, 0.0f, 100}}, false);
    fRegistry.add("MFTMCHMID/primary/correct/hChi2MatchMCHMID", "chi2 match MCH-MID;chi2", kTH1F, {{100, 0.0f, 100}}, false);
    fRegistry.add("MFTMCHMID/primary/correct/hChi2MatchMCHMFT", "chi2 match MCH-MFT;chi2", kTH1F, {{100, 0.0f, 100}}, false);
    fRegistry.add("MFTMCHMID/primary/correct/hDCAxy2D", "DCA x vs. y;DCA_{x} (cm);DCA_{y} (cm)", kTH2F, {{200, -0.5, 0.5}, {200, -0.5, +0.5}}, false);
    fRegistry.add("MFTMCHMID/primary/correct/hDCAxy2DinSigma", "DCA x vs. y in sigma;DCA_{x} (#sigma);DCA_{y} (#sigma)", kTH2F, {{200, -10, 10}, {200, -10, +10}}, false);
    fRegistry.add("MFTMCHMID/primary/correct/hDCAxy", "DCAxy;DCA_{xy} (cm);", kTH1F, {{100, 0, 1}}, false);
    fRegistry.add("MFTMCHMID/primary/correct/hDCAxyinSigma", "DCAxy in sigma;DCA_{xy} (#sigma);", kTH1F, {{100, 0, 10}}, false);
    fRegistry.add("MFTMCHMID/primary/correct/hDCAxResolutionvsPt", "DCA_{x} vs. p_{T};p_{T} (GeV/c);DCA_{x} resolution (#mum);", kTH2F, {{100, 0, 10.f}, {500, 0, 500}}, false);
    fRegistry.add("MFTMCHMID/primary/correct/hDCAyResolutionvsPt", "DCA_{y} vs. p_{T};p_{T} (GeV/c);DCA_{y} resolution (#mum);", kTH2F, {{100, 0, 10.f}, {500, 0, 500}}, false);
    fRegistry.add("MFTMCHMID/primary/correct/hRelDeltaPt", "pT resolution;p_{T}^{gen} (GeV/c);(p_{T}^{rec} - p_{T}^{gen})/p_{T}^{gen}", kTH2F, {{100, 0, 10}, {400, -1, +1}}, true);
    fRegistry.add("MFTMCHMID/primary/correct/hDeltaEta_Pos", "#eta resolution;p_{T}^{gen} (GeV/c);#eta^{rec} - #eta^{gen}", kTH2F, {{100, 0, 10}, {400, -0.2, +0.2}}, true);
    fRegistry.add("MFTMCHMID/primary/correct/hDeltaEta_Neg", "#eta resolution;p_{T}^{gen} (GeV/c);#eta^{rec} - #eta^{gen}", kTH2F, {{100, 0, 10}, {400, -0.2, +0.2}}, true);
    fRegistry.add("MFTMCHMID/primary/correct/hDeltaPhi_Pos", "#varphi resolution;p_{T}^{gen} (GeV/c);#varphi^{rec} - #varphi^{gen} (rad.)", kTH2F, {{100, 0, 10}, {400, -0.2, +0.2}}, true);
    fRegistry.add("MFTMCHMID/primary/correct/hDeltaPhi_Neg", "#varphi resolution;p_{T}^{gen} (GeV/c);#varphi^{rec} - #varphi^{gen} (rad.)", kTH2F, {{100, 0, 10}, {400, -0.2, +0.2}}, true);
    fRegistry.addClone("MFTMCHMID/primary/correct/", "MFTMCHMID/primary/wrong/");
    fRegistry.addClone("MFTMCHMID/primary/", "MFTMCHMID/secondary/");
  }

  bool isSelected(const float pt, const float eta, const float rAtAbsorberEnd, const float pDCA, const float chi2, const uint8_t trackType, const float etaMatchedMCHMID, const float dcaXY)
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
      if (applyEtaCutToSAinGL && (etaMatchedMCHMID < minEtaGL || maxEtaGL < etaMatchedMCHMID)) {
        return false;
      }
      if (maxDCAxy < dcaXY) {
        return false;
      }
      if (chi2 < 0.f || maxChi2GL < chi2) {
        return false;
      }
    } else if (trackType == static_cast<uint8_t>(o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack)) {
      if (eta < minEtaSA || maxEtaSA < eta) {
        return false;
      }
      if (chi2 < 0.f || maxChi2SA < chi2) {
        return false;
      }
    } else {
      return false;
    }

    return true;
  }

  template <typename TCollision, typename TFwdTrack, typename TFwdTracks, typename TMFTTracks>
  void fillHistograms(TCollision const& collision, TFwdTrack fwdtrack, TFwdTracks const&, TMFTTracks const&)
  {
    const auto& mchtrack = fwdtrack.template matchMCHTrack_as<TFwdTracks>(); // MCH-MID
    const auto& mfttrack = fwdtrack.template matchMFTTrack_as<TMFTTracks>();

    if (!fwdtrack.has_mcParticle() || !mchtrack.has_mcParticle() || !mfttrack.has_mcParticle()) {
      return;
    }

    const auto& mcParticle_MFTMCHMID = fwdtrack.template mcParticle_as<aod::McParticles>(); // this is identical to mcParticle_MCHMID
    const auto& mcParticle_MCHMID = mchtrack.template mcParticle_as<aod::McParticles>();    // this is identical to mcParticle_MFTMCHMID
    const auto& mcParticle_MFT = mfttrack.template mcParticle_as<aod::McParticles>();
    // LOGF(info, "mcParticle_MFTMCHMID.pdgCode() = %d, mcParticle_MCHMID.pdgCode() = %d, mcParticle_MFT.pdgCode() = %d", mcParticle_MFTMCHMID.pdgCode(), mcParticle_MCHMID.pdgCode(), mcParticle_MFT.pdgCode());
    // LOGF(info, "mcParticle_MFTMCHMID.globalIndex() = %d, mcParticle_MCHMID.globalIndex() = %d, mcParticle_MFT.globalIndex() = %d", mcParticle_MFTMCHMID.globalIndex(), mcParticle_MCHMID.globalIndex(), mcParticle_MFT.globalIndex());

    if (std::abs(mcParticle_MCHMID.pdgCode()) != 13) { // select true muon
      return;
    }

    if (requireTrueAssociation && (mcParticle_MCHMID.mcCollisionId() != collision.mcCollisionId())) {
      return;
    }

    bool isPrimary = mcParticle_MCHMID.isPhysicalPrimary() || mcParticle_MCHMID.producedByGenerator();
    bool isMatched = (mcParticle_MFT.globalIndex() == mcParticle_MCHMID.globalIndex()) && (mcParticle_MFT.mcCollisionId() == mcParticle_MCHMID.mcCollisionId());

    o2::dataformats::GlobalFwdTrack propmuonAtPV = propagateMuon(fwdtrack, collision, propagationPoint::kToVertex);
    o2::dataformats::GlobalFwdTrack propmuonAtDCA = propagateMuon(fwdtrack, collision, propagationPoint::kToDCA);

    float pt = propmuonAtPV.getPt();
    float eta = propmuonAtPV.getEta();
    float phi = propmuonAtPV.getPhi();
    o2::math_utils::bringTo02Pi(phi);

    float cXXatDCA = propmuonAtDCA.getSigma2X();
    float cYYatDCA = propmuonAtDCA.getSigma2Y();
    float cXYatDCA = propmuonAtDCA.getSigmaXY();

    float dcaX = propmuonAtDCA.getX() - collision.posX();
    float dcaY = propmuonAtDCA.getY() - collision.posY();
    float rAtAbsorberEnd = fwdtrack.rAtAbsorberEnd(); // this works only for GlobalMuonTrack
    float dcaXY = std::sqrt(dcaX * dcaX + dcaY * dcaY);

    float dFdx = 2.f * dcaX / dcaXY;
    float dFdy = 2.f * dcaY / dcaXY;
    float sigma_dcaXY = std::sqrt(cXXatDCA * dFdx * dFdx + cYYatDCA * dFdy * dFdy + 2.f * cXYatDCA * dFdx * dFdy);

    o2::dataformats::GlobalFwdTrack propmuonAtPV_Matched = propagateMuon(mchtrack, collision, propagationPoint::kToVertex);
    float etaMatchedMCHMID = propmuonAtPV_Matched.getEta();
    float phiMatchedMCHMID = propmuonAtPV_Matched.getPhi();
    o2::math_utils::bringTo02Pi(phiMatchedMCHMID);

    o2::dataformats::GlobalFwdTrack propmuonAtDCA_Matched = propagateMuon(mchtrack, collision, propagationPoint::kToDCA);
    float dcaX_Matched = propmuonAtDCA_Matched.getX() - collision.posX();
    float dcaY_Matched = propmuonAtDCA_Matched.getY() - collision.posY();
    float dcaXY_Matched = std::sqrt(dcaX_Matched * dcaX_Matched + dcaY_Matched * dcaY_Matched);
    float pDCA = mchtrack.p() * dcaXY_Matched;

    int nClustersMFT = mfttrack.nClusters();
    float chi2mft = mfttrack.chi2();

    if (refitGlobalMuon) {
      eta = mfttrack.eta();
      phi = mfttrack.phi();
      o2::math_utils::bringTo02Pi(phi);
      pt = propmuonAtPV_Matched.getP() * std::sin(2.f * std::atan(std::exp(-eta)));
    }

    if (!isSelected(pt, eta, rAtAbsorberEnd, pDCA, fwdtrack.chi2(), fwdtrack.trackType(), etaMatchedMCHMID, dcaXY)) {
      return;
    }

    float deta = etaMatchedMCHMID - eta;
    float dphi = phiMatchedMCHMID - phi;
    o2::math_utils::bringToPMPi(dphi);

    fRegistry.fill(HIST("hMuonType"), fwdtrack.trackType());
    if (isPrimary) {
      if (isMatched) {
        fRegistry.fill(HIST("MFTMCHMID/primary/correct/hPt"), pt);
        fRegistry.fill(HIST("MFTMCHMID/primary/correct/hEtaPhi"), phi, eta);
        fRegistry.fill(HIST("MFTMCHMID/primary/correct/hEtaPhi_MatchedMCHMID"), phiMatchedMCHMID, etaMatchedMCHMID);
        fRegistry.fill(HIST("MFTMCHMID/primary/correct/hDeltaEtaDeltaPhi"), dphi, deta);
        fRegistry.fill(HIST("MFTMCHMID/primary/correct/hDiffCollId"), collision.globalIndex() - fwdtrack.collisionId());
        fRegistry.fill(HIST("MFTMCHMID/primary/correct/hSign"), fwdtrack.sign());
        fRegistry.fill(HIST("MFTMCHMID/primary/correct/hNclusters"), fwdtrack.nClusters());
        fRegistry.fill(HIST("MFTMCHMID/primary/correct/hNclustersMFT"), nClustersMFT);
        fRegistry.fill(HIST("MFTMCHMID/primary/correct/hPDCA_Rabs"), rAtAbsorberEnd, pDCA);
        fRegistry.fill(HIST("MFTMCHMID/primary/correct/hRatAbsorberEnd"), rAtAbsorberEnd);
        fRegistry.fill(HIST("MFTMCHMID/primary/correct/hChi2"), fwdtrack.chi2());
        fRegistry.fill(HIST("MFTMCHMID/primary/correct/hChi2MFT"), chi2mft);
        fRegistry.fill(HIST("MFTMCHMID/primary/correct/hChi2MatchMCHMID"), fwdtrack.chi2MatchMCHMID());
        fRegistry.fill(HIST("MFTMCHMID/primary/correct/hChi2MatchMCHMFT"), fwdtrack.chi2MatchMCHMFT());
        fRegistry.fill(HIST("MFTMCHMID/primary/correct/hDCAxy2D"), dcaX, dcaY);
        fRegistry.fill(HIST("MFTMCHMID/primary/correct/hDCAxy2DinSigma"), dcaX / std::sqrt(cXXatDCA), dcaY / std::sqrt(cYYatDCA));
        fRegistry.fill(HIST("MFTMCHMID/primary/correct/hDCAxy"), dcaXY);
        fRegistry.fill(HIST("MFTMCHMID/primary/correct/hDCAxyinSigma"), dcaXY / sigma_dcaXY);
        fRegistry.fill(HIST("MFTMCHMID/primary/correct/hDCAxResolutionvsPt"), pt, std::sqrt(cXXatDCA) * 1e+4); // convert cm to um
        fRegistry.fill(HIST("MFTMCHMID/primary/correct/hDCAyResolutionvsPt"), pt, std::sqrt(cYYatDCA) * 1e+4); // convert cm to um
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
        fRegistry.fill(HIST("MFTMCHMID/primary/wrong/hDeltaEtaDeltaPhi"), dphi, deta);
        fRegistry.fill(HIST("MFTMCHMID/primary/wrong/hDiffCollId"), collision.globalIndex() - fwdtrack.collisionId());
        fRegistry.fill(HIST("MFTMCHMID/primary/wrong/hSign"), fwdtrack.sign());
        fRegistry.fill(HIST("MFTMCHMID/primary/wrong/hNclusters"), fwdtrack.nClusters());
        fRegistry.fill(HIST("MFTMCHMID/primary/wrong/hNclustersMFT"), nClustersMFT);
        fRegistry.fill(HIST("MFTMCHMID/primary/wrong/hPDCA_Rabs"), rAtAbsorberEnd, pDCA);
        fRegistry.fill(HIST("MFTMCHMID/primary/wrong/hRatAbsorberEnd"), rAtAbsorberEnd);
        fRegistry.fill(HIST("MFTMCHMID/primary/wrong/hChi2"), fwdtrack.chi2());
        fRegistry.fill(HIST("MFTMCHMID/primary/wrong/hChi2MFT"), chi2mft);
        fRegistry.fill(HIST("MFTMCHMID/primary/wrong/hChi2MatchMCHMID"), fwdtrack.chi2MatchMCHMID());
        fRegistry.fill(HIST("MFTMCHMID/primary/wrong/hChi2MatchMCHMFT"), fwdtrack.chi2MatchMCHMFT());
        fRegistry.fill(HIST("MFTMCHMID/primary/wrong/hDCAxy2D"), dcaX, dcaY);
        fRegistry.fill(HIST("MFTMCHMID/primary/wrong/hDCAxy2DinSigma"), dcaX / std::sqrt(cXXatDCA), dcaY / std::sqrt(cYYatDCA));
        fRegistry.fill(HIST("MFTMCHMID/primary/wrong/hDCAxy"), dcaXY);
        fRegistry.fill(HIST("MFTMCHMID/primary/wrong/hDCAxyinSigma"), dcaXY / sigma_dcaXY);
        fRegistry.fill(HIST("MFTMCHMID/primary/wrong/hDCAxResolutionvsPt"), pt, std::sqrt(cXXatDCA) * 1e+4); // convert cm to um
        fRegistry.fill(HIST("MFTMCHMID/primary/wrong/hDCAyResolutionvsPt"), pt, std::sqrt(cYYatDCA) * 1e+4); // convert cm to um
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
        fRegistry.fill(HIST("MFTMCHMID/secondary/correct/hDeltaEtaDeltaPhi"), dphi, deta);
        fRegistry.fill(HIST("MFTMCHMID/secondary/correct/hDiffCollId"), collision.globalIndex() - fwdtrack.collisionId());
        fRegistry.fill(HIST("MFTMCHMID/secondary/correct/hSign"), fwdtrack.sign());
        fRegistry.fill(HIST("MFTMCHMID/secondary/correct/hNclusters"), fwdtrack.nClusters());
        fRegistry.fill(HIST("MFTMCHMID/secondary/correct/hNclustersMFT"), nClustersMFT);
        fRegistry.fill(HIST("MFTMCHMID/secondary/correct/hPDCA_Rabs"), rAtAbsorberEnd, pDCA);
        fRegistry.fill(HIST("MFTMCHMID/secondary/correct/hRatAbsorberEnd"), rAtAbsorberEnd);
        fRegistry.fill(HIST("MFTMCHMID/secondary/correct/hChi2"), fwdtrack.chi2());
        fRegistry.fill(HIST("MFTMCHMID/secondary/correct/hChi2MFT"), chi2mft);
        fRegistry.fill(HIST("MFTMCHMID/secondary/correct/hChi2MatchMCHMID"), fwdtrack.chi2MatchMCHMID());
        fRegistry.fill(HIST("MFTMCHMID/secondary/correct/hChi2MatchMCHMFT"), fwdtrack.chi2MatchMCHMFT());
        fRegistry.fill(HIST("MFTMCHMID/secondary/correct/hDCAxy2D"), dcaX, dcaY);
        fRegistry.fill(HIST("MFTMCHMID/secondary/correct/hDCAxy2DinSigma"), dcaX / std::sqrt(cXXatDCA), dcaY / std::sqrt(cYYatDCA));
        fRegistry.fill(HIST("MFTMCHMID/secondary/correct/hDCAxy"), dcaXY);
        fRegistry.fill(HIST("MFTMCHMID/secondary/correct/hDCAxyinSigma"), dcaXY / sigma_dcaXY);
        fRegistry.fill(HIST("MFTMCHMID/secondary/correct/hDCAxResolutionvsPt"), pt, std::sqrt(cXXatDCA) * 1e+4); // convert cm to um
        fRegistry.fill(HIST("MFTMCHMID/secondary/correct/hDCAyResolutionvsPt"), pt, std::sqrt(cYYatDCA) * 1e+4); // convert cm to um
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
        fRegistry.fill(HIST("MFTMCHMID/secondary/wrong/hDeltaEtaDeltaPhi"), dphi, deta);
        fRegistry.fill(HIST("MFTMCHMID/secondary/wrong/hDiffCollId"), collision.globalIndex() - fwdtrack.collisionId());
        fRegistry.fill(HIST("MFTMCHMID/secondary/wrong/hSign"), fwdtrack.sign());
        fRegistry.fill(HIST("MFTMCHMID/secondary/wrong/hNclusters"), fwdtrack.nClusters());
        fRegistry.fill(HIST("MFTMCHMID/secondary/wrong/hNclustersMFT"), nClustersMFT);
        fRegistry.fill(HIST("MFTMCHMID/secondary/wrong/hPDCA_Rabs"), rAtAbsorberEnd, pDCA);
        fRegistry.fill(HIST("MFTMCHMID/secondary/wrong/hRatAbsorberEnd"), rAtAbsorberEnd);
        fRegistry.fill(HIST("MFTMCHMID/secondary/wrong/hChi2"), fwdtrack.chi2());
        fRegistry.fill(HIST("MFTMCHMID/secondary/wrong/hChi2MFT"), chi2mft);
        fRegistry.fill(HIST("MFTMCHMID/secondary/wrong/hChi2MatchMCHMID"), fwdtrack.chi2MatchMCHMID());
        fRegistry.fill(HIST("MFTMCHMID/secondary/wrong/hChi2MatchMCHMFT"), fwdtrack.chi2MatchMCHMFT());
        fRegistry.fill(HIST("MFTMCHMID/secondary/wrong/hDCAxy2D"), dcaX, dcaY);
        fRegistry.fill(HIST("MFTMCHMID/secondary/wrong/hDCAxy2DinSigma"), dcaX / std::sqrt(cXXatDCA), dcaY / std::sqrt(cYYatDCA));
        fRegistry.fill(HIST("MFTMCHMID/secondary/wrong/hDCAxy"), dcaXY);
        fRegistry.fill(HIST("MFTMCHMID/secondary/wrong/hDCAxyinSigma"), dcaXY / sigma_dcaXY);
        fRegistry.fill(HIST("MFTMCHMID/secondary/wrong/hDCAxResolutionvsPt"), pt, std::sqrt(cXXatDCA) * 1e+4); // convert cm to um
        fRegistry.fill(HIST("MFTMCHMID/secondary/wrong/hDCAyResolutionvsPt"), pt, std::sqrt(cYYatDCA) * 1e+4); // convert cm to um
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

  SliceCache cache;
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
      const auto& bc = collision.template bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);
      fRegistry.fill(HIST("Event/hCollisionCounter"), 0);
      if (!collision.has_mcCollision()) {
        continue;
      }
      const float centralities[3] = {collision.centFT0M(), collision.centFT0A(), collision.centFT0C()};
      if (centralities[cfgCentEstimator] < cfgCentMin || cfgCentMax < centralities[cfgCentEstimator]) {
        continue;
      }
      fRegistry.fill(HIST("Event/hCollisionCounter"), 1);
      fillEventHistograms(collision);

      const auto& fwdtracks_per_coll = fwdtracks.sliceBy(perCollision, collision.globalIndex());
      for (const auto& fwdtrack : fwdtracks_per_coll) {
        if (fwdtrack.trackType() != o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack) {
          continue;
        }
        fillHistograms(collision, fwdtrack, fwdtracks, mfttracks);
      } // end of fwdtrack loop
    } // end of collision loop
  }
  PROCESS_SWITCH(matchingMFT, processWithoutFTTCA, "process without FTTCA", false);

  void processWithFTTCA(FilteredMyCollisions const& collisions, MyFwdTracks const& fwdtracks, MyMFTTracks const& mfttracks, aod::BCsWithTimestamps const&, aod::FwdTrackAssoc const& fwdtrackIndices, aod::McParticles const&)
  {
    for (const auto& collision : collisions) {
      const auto& bc = collision.template bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);
      fRegistry.fill(HIST("Event/hCollisionCounter"), 0);
      if (!collision.has_mcCollision()) {
        continue;
      }
      const float centralities[3] = {collision.centFT0M(), collision.centFT0A(), collision.centFT0C()};
      if (centralities[cfgCentEstimator] < cfgCentMin || cfgCentMax < centralities[cfgCentEstimator]) {
        continue;
      }
      fRegistry.fill(HIST("Event/hCollisionCounter"), 1);
      fillEventHistograms(collision);

      const auto& fwdtrackIdsThisCollision = fwdtrackIndices.sliceBy(fwdtrackIndicesPerCollision, collision.globalIndex());
      for (const auto& fwdtrackId : fwdtrackIdsThisCollision) {
        const auto& fwdtrack = fwdtrackId.template fwdtrack_as<MyFwdTracks>();
        if (fwdtrack.trackType() != o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack && fwdtrack.trackType()) {
          continue;
        }
        fillHistograms(collision, fwdtrack, fwdtracks, mfttracks);
      } // end of fwdtrack loop
    } // end of collision loop
  }
  PROCESS_SWITCH(matchingMFT, processWithFTTCA, "process with FTTCA", true);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<matchingMFT>(cfgc, TaskName{"matching-mft"})};
}
