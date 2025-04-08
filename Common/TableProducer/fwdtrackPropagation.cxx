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

/// \file fwdtrackPropagator.cxx
/// \brief Common task to produce propagated forward tracks
/// \author Maurice Coquet <maurice.louis.coquet@cern.ch>
/// \author Luca Micheletti <luca.micheletti@cern.ch>
/// \author Daiki Sekihata <daiki.sekihata@cern.ch>

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

using namespace o2;
using namespace o2::aod;
using namespace o2::soa;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::fwdtrackutils;

struct FwdTrackPropagation {
  using MyFwdTracks = soa::Join<aod::FwdTracks, aod::FwdTracksCov>;

  Produces<aod::StoredPropagatedFwdTracks> propfwdtracks;
  Produces<aod::StoredPropagatedFwdTracksCov> propfwdtrackscov;

  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
  Configurable<bool> fillQAHistograms{"fillQAHistograms", true, "flag to fill QA histograms"};
  Configurable<float> minPt{"minPt", 0.2, "min pt for muon"};
  Configurable<float> maxPt{"maxPt", 1e+10, "max pt for muon"};
  Configurable<float> minEtaSA{"minEtaSA", -4.0, "min. eta acceptance for MCH-MID"};
  Configurable<float> maxEtaSA{"maxEtaSA", -2.5, "max. eta acceptance for MCH-MID"};
  Configurable<float> minEtaGL{"minEtaGL", -3.6, "min. eta acceptance for MFT-MCH-MID"};
  Configurable<float> maxEtaGL{"maxEtaGL", -2.5, "max. eta acceptance for MFT-MCH-MID"};
  Configurable<float> minRabsGL{"minRabsGL", 27.6, "min. R at absorber end for global muon (min. eta = -3.6)"}; // std::tan(2.f * std::atan(std::exp(- -3.6)) ) * -505.
  Configurable<float> minRabs{"minRabs", 17.6, "min. R at absorber end"};
  Configurable<float> midRabs{"midRabs", 26.5, "middle R at absorber end for pDCA cut"};
  Configurable<float> maxRabs{"maxRabs", 89.5, "max. R at absorber end"};
  Configurable<float> maxDCAxy{"maxDCAxy", 1e+10, "max. DCAxy for global muons"};
  Configurable<float> maxPDCAforLargeR{"maxPDCAforLargeR", 324.f, "max. pDCA for large R at absorber end"};
  Configurable<float> maxPDCAforSmallR{"maxPDCAforSmallR", 594.f, "max. pDCA for small R at absorber end"};
  Configurable<float> maxMatchingChi2MCHMFT{"maxMatchingChi2MCHMFT", 50.f, "max. chi2 for MCH-MFT matching"};
  Configurable<float> maxChi2SA{"maxChi2SA", 1e+6, "max. chi2 for standalone muon"};
  Configurable<float> maxChi2GL{"maxChi2GL", 50.f, "max. chi2 for global muon"};
  Configurable<bool> refitGlobalMuon{"refitGlobalMuon", true, "flag to refit global muon"};

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

    if (fillQAHistograms) {
      addHistograms();
    }
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
    auto hMuonType = fRegistry.add<TH1>("hMuonType", "muon type", kTH1F, {{5, -0.5f, 4.5f}}, false);
    hMuonType->GetXaxis()->SetBinLabel(1, "MFT-MCH-MID (global muon)");
    hMuonType->GetXaxis()->SetBinLabel(2, "MFT-MCH-MID (global muon other match)");
    hMuonType->GetXaxis()->SetBinLabel(3, "MFT-MCH");
    hMuonType->GetXaxis()->SetBinLabel(4, "MCH-MID");
    hMuonType->GetXaxis()->SetBinLabel(5, "MCH standalone");

    fRegistry.add("MFTMCHMID/hPt", "pT;p_{T} (GeV/c)", kTH1F, {{100, 0.0f, 10}}, false);
    fRegistry.add("MFTMCHMID/hRelDiffPt", "pT resolution;p_{T} (GeV/c);#Deltap_{T}/p_{T}", kTH2F, {{100, 0.0f, 10}, {200, 0, 0.2}}, false);
    fRegistry.add("MFTMCHMID/hEtaPhi", "#eta vs. #varphi;#varphi (rad.);#eta", kTH2F, {{180, 0, 2 * M_PI}, {60, -5.f, -2.f}}, false);
    fRegistry.add("MFTMCHMID/hEtaPhi_MatchedMCHMID", "#eta vs. #varphi;#varphi (rad.);#eta", kTH2F, {{180, 0, 2 * M_PI}, {60, -5.f, -2.f}}, false);
    fRegistry.add("MFTMCHMID/hDiffCollId", "difference in collision index;collisionId_{TTCA} - collisionId_{MP}", kTH1F, {{41, -20.5, +20.5}}, false);
    fRegistry.add("MFTMCHMID/hSign", "sign;sign", kTH1F, {{3, -1.5, +1.5}}, false);
    fRegistry.add("MFTMCHMID/hNclusters", "Nclusters;Nclusters", kTH1F, {{21, -0.5f, 20.5}}, false);
    fRegistry.add("MFTMCHMID/hNclustersMFT", "NclustersMFT;Nclusters MFT", kTH1F, {{11, -0.5f, 10.5}}, false);
    fRegistry.add("MFTMCHMID/hRatAbsorberEnd", "R at absorber end;R at absorber end (cm)", kTH1F, {{100, 0.0f, 100}}, false);
    fRegistry.add("MFTMCHMID/hPDCA_Rabs", "pDCA vs. Rabs;R at absorber end (cm);p #times DCA (GeV/c #upoint cm)", kTH2F, {{100, 0, 100}, {100, 0.0f, 1000}}, false);
    fRegistry.add("MFTMCHMID/hChi2", "chi2;chi2", kTH1F, {{100, 0.0f, 100}}, false);
    fRegistry.add("MFTMCHMID/hChi2MFT", "chi2 MFT;chi2 MFT", kTH1F, {{100, 0.0f, 100}}, false);
    fRegistry.add("MFTMCHMID/hChi2MatchMCHMID", "chi2 match MCH-MID;chi2", kTH1F, {{100, 0.0f, 100}}, false);
    fRegistry.add("MFTMCHMID/hChi2MatchMCHMFT", "chi2 match MCH-MFT;chi2", kTH1F, {{100, 0.0f, 100}}, false);
    fRegistry.add("MFTMCHMID/hMatchScoreMCHMFT", "match score MCH-MFT;score", kTH1F, {{100, 0.0f, 100}}, false);
    fRegistry.add("MFTMCHMID/hDCAxy2D", "DCA x vs. y;DCA_{x} (cm);DCA_{y} (cm)", kTH2F, {{200, -0.5, 0.5}, {200, -0.5, +0.5}}, false);
    fRegistry.add("MFTMCHMID/hDCAxy2DinSigma", "DCA x vs. y in sigma;DCA_{x} (#sigma);DCA_{y} (#sigma)", kTH2F, {{200, -10, 10}, {200, -10, +10}}, false);
    fRegistry.add("MFTMCHMID/hDCAxy", "DCAxy;DCA_{xy} (cm);", kTH1F, {{100, 0, 1}}, false);
    fRegistry.add("MFTMCHMID/hDCAxyinSigma", "DCAxy in sigma;DCA_{xy} (#sigma);", kTH1F, {{100, 0, 10}}, false);
    fRegistry.addClone("MFTMCHMID/", "MCHMID/");
    fRegistry.add("MFTMCHMID/hDCAxResolutionvsPt", "DCA_{x} vs. p_{T};p_{T} (GeV/c);DCA_{x} resolution (#mum);", kTH2F, {{100, 0, 10.f}, {500, 0, 500}}, false);
    fRegistry.add("MFTMCHMID/hDCAyResolutionvsPt", "DCA_{y} vs. p_{T};p_{T} (GeV/c);DCA_{y} resolution (#mum);", kTH2F, {{100, 0, 10.f}, {500, 0, 500}}, false);
    fRegistry.add("MCHMID/hDCAxResolutionvsPt", "DCA_{x} vs. p_{T};p_{T} (GeV/c);DCA_{x} resolution (#mum);", kTH2F, {{100, 0, 10.f}, {500, 0, 5e+5}}, false);
    fRegistry.add("MCHMID/hDCAyResolutionvsPt", "DCA_{y} vs. p_{T};p_{T} (GeV/c);DCA_{y} resolution (#mum);", kTH2F, {{100, 0, 10.f}, {500, 0, 5e+5}}, false);
  }

  bool isSelected(const float pt, const float eta, const float rAtAbsorberEnd, const float pDCA, const float chi2, const uint8_t trackType, const float dcaXY)
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
      if (chi2 < 0.f || maxChi2GL < chi2) {
        return false;
      }
      if (rAtAbsorberEnd < minRabsGL || maxRabs < rAtAbsorberEnd) {
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
  void fillFwdTrackTable(TCollision const& collision, TFwdTrack fwdtrack, TFwdTracks const&, TMFTTracks const&, const bool isAmbiguous)
  {
    if (fwdtrack.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack && (fwdtrack.chi2MatchMCHMFT() > maxMatchingChi2MCHMFT || fwdtrack.chi2() > maxChi2GL)) {
      return;
    } // Users have to decide the best match between MFT and MCH-MID at analysis level. The same global muon is repeatedly stored.

    if (fwdtrack.chi2MatchMCHMID() < 0.f) { // this should never happen. only for protection.
      return;
    }

    o2::dataformats::GlobalFwdTrack propmuonAtPV = propagateMuon(fwdtrack, collision, propagationPoint::kToVertex);
    o2::dataformats::GlobalFwdTrack propmuonAtDCA = propagateMuon(fwdtrack, collision, propagationPoint::kToDCA);

    float pt = propmuonAtPV.getPt();
    float eta = propmuonAtPV.getEta();
    float phi = propmuonAtPV.getPhi();
    float tgl = propmuonAtPV.getTgl();
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

    float pDCA = fwdtrack.p() * dcaXY;
    int nClustersMFT = 0;
    float etaMatchedMCHMID = propmuonAtPV.getEta();
    float phiMatchedMCHMID = propmuonAtPV.getPhi();
    o2::math_utils::bringTo02Pi(phiMatchedMCHMID);
    float x = fwdtrack.x();
    float y = fwdtrack.y();
    float z = fwdtrack.z();
    float chi2mft = 0.f;

    if (fwdtrack.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack) {
      const auto& mchtrack = fwdtrack.template matchMCHTrack_as<TFwdTracks>(); // MCH-MID
      o2::dataformats::GlobalFwdTrack propmuonAtPV_Matched = propagateMuon(mchtrack, collision, propagationPoint::kToVertex);
      etaMatchedMCHMID = propmuonAtPV_Matched.getEta();
      phiMatchedMCHMID = propmuonAtPV_Matched.getPhi();
      o2::math_utils::bringTo02Pi(phiMatchedMCHMID);
      o2::dataformats::GlobalFwdTrack propmuonAtDCA_Matched = propagateMuon(mchtrack, collision, propagationPoint::kToDCA);
      float dcaX_Matched = propmuonAtDCA_Matched.getX() - collision.posX();
      float dcaY_Matched = propmuonAtDCA_Matched.getY() - collision.posY();
      float dcaXY_Matched = std::sqrt(dcaX_Matched * dcaX_Matched + dcaY_Matched * dcaY_Matched);
      pDCA = mchtrack.p() * dcaXY_Matched;

      const auto& mfttrack = fwdtrack.template matchMFTTrack_as<TMFTTracks>();
      nClustersMFT = mfttrack.nClusters();
      chi2mft = mfttrack.chi2();
      if (refitGlobalMuon) {
        eta = mfttrack.eta();
        phi = mfttrack.phi();
        o2::math_utils::bringTo02Pi(phi);
        pt = propmuonAtPV_Matched.getP() * std::sin(2.f * std::atan(std::exp(-eta)));

        x = mfttrack.x();
        y = mfttrack.y();
        z = mfttrack.z();
        tgl = mfttrack.tgl();
      }
    } else if (fwdtrack.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack) {
      o2::dataformats::GlobalFwdTrack propmuonAtRabs = propagateMuon(fwdtrack, collision, propagationPoint::kToRabs); // this is necessary only for MuonStandaloneTrack
      float xAbs = propmuonAtRabs.getX();
      float yAbs = propmuonAtRabs.getY();
      rAtAbsorberEnd = std::sqrt(xAbs * xAbs + yAbs * yAbs); // Redo propagation only for muon tracks // propagation of MFT tracks alredy done in reconstruction
    } else {
      return;
    }

    if (!isSelected(pt, eta, rAtAbsorberEnd, pDCA, fwdtrack.chi2(), fwdtrack.trackType(), dcaXY)) {
      return;
    }

    const auto& fwdcov = propmuonAtPV.getCovariances(); // covatiant matrix at PV
    const float sigX = std::sqrt(fwdcov(0, 0));
    const float sigY = std::sqrt(fwdcov(1, 1));
    const float sigPhi = std::sqrt(fwdcov(2, 2));
    const float sigTgl = std::sqrt(fwdcov(3, 3));
    const float sig1Pt = std::sqrt(fwdcov(4, 4));
    const float rhoXY = 128.f * fwdcov(0, 1) / (sigX * sigY);
    const float rhoPhiX = 128.f * fwdcov(0, 2) / (sigPhi * sigX);
    const float rhoPhiY = 128.f * fwdcov(1, 2) / (sigPhi * sigY);
    const float rhoTglX = 128.f * fwdcov(0, 3) / (sigTgl * sigX);
    const float rhoTglY = 128.f * fwdcov(1, 3) / (sigTgl * sigY);
    const float rhoTglPhi = 128.f * fwdcov(2, 3) / (sigTgl * sigPhi);
    const float rho1PtX = 128.f * fwdcov(0, 4) / (sig1Pt * sigX);
    const float rho1PtY = 128.f * fwdcov(1, 4) / (sig1Pt * sigY);
    const float rho1PtPhi = 128.f * fwdcov(2, 4) / (sig1Pt * sigPhi);
    const float rho1PtTgl = 128.f * fwdcov(3, 4) / (sig1Pt * sigTgl);

    bool isAssociatedToMPC = fwdtrack.collisionId() == collision.globalIndex();
    // LOGF(info, "isAmbiguous = %d, isAssociatedToMPC = %d, fwdtrack.globalIndex() = %d, fwdtrack.collisionId() = %d, collision.globalIndex() = %d", isAmbiguous, isAssociatedToMPC, fwdtrack.globalIndex(), fwdtrack.collisionId(), collision.globalIndex());

    propfwdtracks(
      collision.globalIndex(), fwdtrack.trackType(),
      x, y, z, phi, tgl,
      fwdtrack.sign() / pt, fwdtrack.nClusters(), pDCA, rAtAbsorberEnd,
      fwdtrack.chi2(), fwdtrack.chi2MatchMCHMID(), fwdtrack.chi2MatchMCHMFT(),
      fwdtrack.matchScoreMCHMFT(), fwdtrack.globalIndex(), fwdtrack.matchMFTTrackId(), fwdtrack.matchMCHTrackId(),
      fwdtrack.mchBitMap(), fwdtrack.midBitMap(), fwdtrack.midBoards(), fwdtrack.trackTime(), fwdtrack.trackTimeRes(), dcaX, dcaY,
      cXXatDCA, cYYatDCA, cXYatDCA, etaMatchedMCHMID, phiMatchedMCHMID, isAssociatedToMPC, isAmbiguous);

    propfwdtrackscov(
      sigX, sigY, sigPhi, sigTgl, sig1Pt,
      rhoXY, rhoPhiX, rhoPhiY, rhoTglX, rhoTglY,
      rhoTglPhi, rho1PtX, rho1PtY, rho1PtPhi, rho1PtTgl);

    if (fillQAHistograms) {
      fRegistry.fill(HIST("hMuonType"), fwdtrack.trackType());
      if (fwdtrack.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack) {
        fRegistry.fill(HIST("MFTMCHMID/hPt"), pt);
        fRegistry.fill(HIST("MFTMCHMID/hRelDiffPt"), pt, sig1Pt * pt);
        fRegistry.fill(HIST("MFTMCHMID/hEtaPhi"), phi, eta);
        fRegistry.fill(HIST("MFTMCHMID/hEtaPhi_MatchedMCHMID"), phiMatchedMCHMID, etaMatchedMCHMID);
        fRegistry.fill(HIST("MFTMCHMID/hDiffCollId"), collision.globalIndex() - fwdtrack.collisionId());
        fRegistry.fill(HIST("MFTMCHMID/hSign"), fwdtrack.sign());
        fRegistry.fill(HIST("MFTMCHMID/hNclusters"), fwdtrack.nClusters());
        fRegistry.fill(HIST("MFTMCHMID/hNclustersMFT"), nClustersMFT);
        fRegistry.fill(HIST("MFTMCHMID/hPDCA_Rabs"), rAtAbsorberEnd, pDCA);
        fRegistry.fill(HIST("MFTMCHMID/hRatAbsorberEnd"), rAtAbsorberEnd);
        fRegistry.fill(HIST("MFTMCHMID/hChi2"), fwdtrack.chi2());
        fRegistry.fill(HIST("MFTMCHMID/hChi2MFT"), chi2mft);
        fRegistry.fill(HIST("MFTMCHMID/hChi2MatchMCHMID"), fwdtrack.chi2MatchMCHMID());
        fRegistry.fill(HIST("MFTMCHMID/hChi2MatchMCHMFT"), fwdtrack.chi2MatchMCHMFT());
        fRegistry.fill(HIST("MFTMCHMID/hMatchScoreMCHMFT"), fwdtrack.matchScoreMCHMFT());
        fRegistry.fill(HIST("MFTMCHMID/hDCAxy2D"), dcaX, dcaY);
        fRegistry.fill(HIST("MFTMCHMID/hDCAxy2DinSigma"), dcaX / std::sqrt(cXXatDCA), dcaY / std::sqrt(cYYatDCA));
        fRegistry.fill(HIST("MFTMCHMID/hDCAxy"), dcaXY);
        fRegistry.fill(HIST("MFTMCHMID/hDCAxyinSigma"), dcaXY / sigma_dcaXY);
        fRegistry.fill(HIST("MFTMCHMID/hDCAxResolutionvsPt"), pt, std::sqrt(cXXatDCA) * 1e+4); // convert cm to um
        fRegistry.fill(HIST("MFTMCHMID/hDCAyResolutionvsPt"), pt, std::sqrt(cYYatDCA) * 1e+4); // convert cm to um
      } else if (fwdtrack.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack) {
        fRegistry.fill(HIST("MCHMID/hPt"), pt);
        fRegistry.fill(HIST("MCHMID/hRelDiffPt"), pt, sig1Pt * pt);
        fRegistry.fill(HIST("MCHMID/hEtaPhi"), phi, eta);
        fRegistry.fill(HIST("MCHMID/hEtaPhi_MatchedMCHMID"), phiMatchedMCHMID, etaMatchedMCHMID);
        fRegistry.fill(HIST("MCHMID/hDiffCollId"), collision.globalIndex() - fwdtrack.collisionId());
        fRegistry.fill(HIST("MCHMID/hSign"), fwdtrack.sign());
        fRegistry.fill(HIST("MCHMID/hNclusters"), fwdtrack.nClusters());
        fRegistry.fill(HIST("MCHMID/hNclustersMFT"), nClustersMFT);
        fRegistry.fill(HIST("MCHMID/hPDCA_Rabs"), rAtAbsorberEnd, pDCA);
        fRegistry.fill(HIST("MCHMID/hRatAbsorberEnd"), rAtAbsorberEnd);
        fRegistry.fill(HIST("MCHMID/hChi2"), fwdtrack.chi2());
        fRegistry.fill(HIST("MCHMID/hChi2MFT"), chi2mft);
        fRegistry.fill(HIST("MCHMID/hChi2MatchMCHMID"), fwdtrack.chi2MatchMCHMID());
        fRegistry.fill(HIST("MCHMID/hChi2MatchMCHMFT"), fwdtrack.chi2MatchMCHMFT());
        fRegistry.fill(HIST("MCHMID/hMatchScoreMCHMFT"), fwdtrack.matchScoreMCHMFT());
        fRegistry.fill(HIST("MCHMID/hDCAxy2D"), dcaX, dcaY);
        fRegistry.fill(HIST("MCHMID/hDCAxy2DinSigma"), dcaX / std::sqrt(cXXatDCA), dcaY / std::sqrt(cYYatDCA));
        fRegistry.fill(HIST("MCHMID/hDCAxy"), dcaXY);
        fRegistry.fill(HIST("MCHMID/hDCAxyinSigma"), dcaXY / sigma_dcaXY);
        fRegistry.fill(HIST("MCHMID/hDCAxResolutionvsPt"), pt, std::sqrt(cXXatDCA) * 1e+4); // convert cm to um
        fRegistry.fill(HIST("MCHMID/hDCAyResolutionvsPt"), pt, std::sqrt(cYYatDCA) * 1e+4); // convert cm to um
      }
    }
  }

  SliceCache cache;
  PresliceUnsorted<aod::FwdTracks> perMFTTrack = o2::aod::fwdtrack::matchMFTTrackId;
  Preslice<aod::FwdTracks> perCollision = o2::aod::fwdtrack::collisionId;
  // Preslice<aod::MFTTracks> perCollisionMFT = o2::aod::fwdtrack::collisionId;
  Preslice<aod::FwdTrackAssoc> fwdtrackIndicesPerCollision = aod::track_association::collisionId;
  PresliceUnsorted<aod::FwdTrackAssoc> fwdtrackIndicesPerFwdTrack = aod::track_association::fwdtrackId;

  void processWithoutFTTCA(aod::Collisions const& collisions, MyFwdTracks const& fwdtracks, aod::MFTTracks const& mfttracks, aod::BCsWithTimestamps const&)
  {
    for (const auto& collision : collisions) {
      const auto& bc = collision.template bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);

      const auto& fwdtracks_per_coll = fwdtracks.sliceBy(perCollision, collision.globalIndex());
      for (const auto& fwdtrack : fwdtracks_per_coll) {
        if (fwdtrack.trackType() != o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack && fwdtrack.trackType() != o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack) {
          continue;
        }
        fillFwdTrackTable(collision, fwdtrack, fwdtracks, mfttracks, false);
      } // end of fwdtrack loop
    } // end of collision loop
  }
  PROCESS_SWITCH(FwdTrackPropagation, processWithoutFTTCA, "process without FTTCA", true);

  void processWithFTTCA(aod::Collisions const& collisions, MyFwdTracks const& fwdtracks, aod::MFTTracks const& mfttracks, aod::BCsWithTimestamps const&, aod::FwdTrackAssoc const& fwdtrackIndices)
  {
    std::unordered_map<int64_t, bool> mapAmb; // fwdtrack.globalIndex() -> bool isAmb;
    for (const auto& fwdtrack : fwdtracks) {
      const auto& fwdtrackIdsPerFwdTrack = fwdtrackIndices.sliceBy(fwdtrackIndicesPerFwdTrack, fwdtrack.globalIndex());
      mapAmb[fwdtrack.globalIndex()] = fwdtrackIdsPerFwdTrack.size() > 1;
      // LOGF(info, "fwdtrack.globalIndex() = %d, ntimes = %d, isAmbiguous = %d", fwdtrack.globalIndex(), fwdtrackIdsPerFwdTrack.size(), mapAmb[fwdtrack.globalIndex()]);
    } // end of fwdtrack loop

    for (const auto& collision : collisions) {
      const auto& bc = collision.template bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);

      const auto& fwdtrackIdsThisCollision = fwdtrackIndices.sliceBy(fwdtrackIndicesPerCollision, collision.globalIndex());
      for (const auto& fwdtrackId : fwdtrackIdsThisCollision) {
        const auto& fwdtrack = fwdtrackId.template fwdtrack_as<MyFwdTracks>();
        if (fwdtrack.trackType() != o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack && fwdtrack.trackType() != o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack) {
          continue;
        }
        fillFwdTrackTable(collision, fwdtrack, fwdtracks, mfttracks, mapAmb[fwdtrack.globalIndex()]);
      } // end of fwdtrack loop
    } // end of collision loop
    mapAmb.clear();
  }
  PROCESS_SWITCH(FwdTrackPropagation, processWithFTTCA, "process with FTTCA", false);
};

// Extends the PropagatedFwdTracks table for expression columns
struct PropagatedFwdTrackSpawner {
  Spawns<aod::PropagatedFwdTracks> propFwdTracks;
  Spawns<aod::PropagatedFwdTracksCov> propFwdTracksCov;
  void init(InitContext const&) {}
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<FwdTrackPropagation>(cfgc, TaskName{"fwdtrack-propagation"}),
    adaptAnalysisTask<PropagatedFwdTrackSpawner>(cfgc, TaskName{"propagated-fwdtrack-spawner"}),
  };
}
