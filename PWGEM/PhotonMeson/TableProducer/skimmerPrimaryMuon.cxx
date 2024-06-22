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

/// \brief write relevant information for dalitz ee analysis to an AO2D.root file. This file is then the only necessary input to perform pcm analysis.
/// \author daiki.sekihata@cern.ch

#include "Math/Vector4D.h"
#include "Math/SMatrix.h"

#include "Framework/DataTypes.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "CommonConstants/PhysicsConstants.h"

#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "TGeoGlobalMagField.h"
#include "Field/MagneticField.h"

#include "DetectorsBase/Propagator.h"
#include "GlobalTracking/MatchGlobalFwd.h"
#include "MCHTracking/TrackExtrap.h"
#include "MCHTracking/TrackParam.h"
#include "ReconstructionDataFormats/TrackFwd.h"

#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"

using namespace o2;
using namespace o2::soa;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;

using MyTracks = soa::Join<aod::FwdTracks, aod::FwdTracksCov>; // muon tracks are repeated. i.e. not exclusive.
using MyTracksMC = soa::Join<MyTracks, aod::McFwdTrackLabels>;

struct skimmerPrimaryMuon {
  enum class EM_MuMuPairType : int {
    kULS = 0,
    kLSpp = +1,
    kLSnn = -1,
  };

  // Index used to set different options for Muon propagation
  enum class MuonExtrapolation : int {
    kToVertex = 0, // propagtion to vertex by default
    kToDCA = 1,
    kToRabs = 2,
  };

  using SMatrix55 = ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepSym<double, 5>>;
  using SMatrix5 = ROOT::Math::SVector<double, 5>;

  SliceCache cache;
  Preslice<aod::FwdTracks> perCollision = o2::aod::fwdtrack::collisionId;
  Preslice<aod::MFTTracks> perCollision_mft = o2::aod::fwdtrack::collisionId;
  Produces<aod::EMPrimaryMuons> emprimarymuons;
  Produces<aod::EMPrimaryMuonsCov> emprimarymuonscov;

  // Configurables
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
  Configurable<float> minpt{"minpt", 0.1, "min pt for track"};
  Configurable<float> mineta{"mineta", -4.0, "eta acceptance"};
  Configurable<float> maxeta{"maxeta", -2.5, "eta acceptance"};
  Configurable<float> mineta_mft{"mineta_mft", -3.6, "eta acceptance"};
  Configurable<float> maxeta_mft{"maxeta_mft", -2.5, "eta acceptance"};
  Configurable<float> minRabs{"minRabs", 17.6, "min. R at absorber end"};
  Configurable<float> maxRabs{"maxRabs", 89.5, "max. R at absorber end"};
  Configurable<float> maxPDCA{"maxPDCA", 1e+3, "max. p DCA to reject beam-gas background"};

  o2::ccdb::CcdbApi ccdbApi;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  int mRunNumber;

  o2::globaltracking::MatchGlobalFwd mMatching;
  HistogramRegistry fRegistry{"output", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};
  static constexpr std::string_view muon_types[5] = {"MFTMCHMID/", "MFTMCHMIDOtherMatch/", "MFTMCH/", "MCHMID/", "MCH/"};

  void init(InitContext const&)
  {
    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);
    ccdbApi.init(ccdburl);

    o2::mch::TrackExtrap::setField();
    addHistograms();

    mRunNumber = 0;
  }

  void initCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }

    mRunNumber = bc.runNumber();
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
    fRegistry.add("Event/hCollisionCounter", "collision counter", kTH1F, {{2, -0.5f, 1.5}}, false);
    fRegistry.add("Event/hNmuon", "Number of #mu per event;N_{#mu^{#minus}};N_{#mu^{+}}", kTH2F, {{101, -0.5f, 100.5}, {101, -0.5f, 100.5}}, false);

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
    fRegistry.add("Track/MFTMCHMID/hNclustersMFT", "NclustersMFT;Nclusters MFT", kTH1F, {{21, -0.5f, 20.5}}, false);
    fRegistry.add("Track/MFTMCHMID/hRatAbsorberEnd", "R at absorber end;R at absorber end (cm)", kTH1F, {{200, 0.0f, 200}}, false);
    fRegistry.add("Track/MFTMCHMID/hPDCA", "pDCA;p_{T,#mu} at PV (GeV/c);p #times DCA (GeV/c #upoint cm)", kTH2F, {{100, 0, 10}, {100, 0.0f, 1000}}, false);
    fRegistry.add("Track/MFTMCHMID/hPDCA_recalc", "pDCA relcalculated;p_{T,#mu} at PV (GeV/c);p #times DCA (GeV/c #upoint cm)", kTH2F, {{100, 0, 10}, {100, 0.0f, 1000}}, false);
    fRegistry.add("Track/MFTMCHMID/hChi2", "chi2;chi2", kTH1F, {{100, 0.0f, 100}}, false);
    fRegistry.add("Track/MFTMCHMID/hChi2MatchMCHMID", "chi2 match MCH-MID;chi2", kTH1F, {{100, 0.0f, 100}}, false);
    fRegistry.add("Track/MFTMCHMID/hChi2MatchMCHMFT", "chi2 match MCH-MFT;chi2", kTH1F, {{100, 0.0f, 100}}, false);
    fRegistry.add("Track/MFTMCHMID/hMatchScoreMCHMFT", "match score MCH-MFT;score", kTH1F, {{100, 0.0f, 100}}, false);
    fRegistry.add("Track/MFTMCHMID/hDCAxy2D", "DCA XY;DCA_{x} (cm);DCA_{y} (cm)", kTH2F, {{200, -10, 10}, {200, -10, +10}}, false);
    fRegistry.add("Track/MFTMCHMID/hDCAxy2DinSigma", "DCA XY;DCA_{x} (#sigma);DCA_{y} (#sigma)", kTH2F, {{200, -10, 10}, {200, -10, +10}}, false);
    fRegistry.add("Track/MFTMCHMID/hDCAxySigma", "DCA_{XY};DCA_{XY} (#sigma);", kTH1F, {{100, 0, 10}}, false);
    fRegistry.add("Track/MFTMCHMID/hDCAxResolutionvsPt", "DCA_{X} vs. p_{T,#mu};p_{T,#mu} (GeV/c);DCA_{X} resolution (#mum);", kTH2F, {{100, 0, 10.f}, {100, 0, 1000}}, false);
    fRegistry.add("Track/MFTMCHMID/hDCAyResolutionvsPt", "DCA_{Y} vs. p_{T,#mu};p_{T,#mu} (GeV/c);DCA_{Y} resolution (#mum);", kTH2F, {{100, 0, 10.f}, {100, 0, 1000}}, false);
    fRegistry.addClone("Track/MFTMCHMID/", "Track/MCHMID/");

    // for pair
    fRegistry.add("Pair/MFTMCHMID/uls/hMvsPt", "m_{#mu#mu} vs. p_{T,#mu#mu};m_{#mu#mu} (GeV/c^{2});p_{T,#mu#mu} (GeV/c)", kTH2F, {{1200, 0.f, 12.f}, {100, 0, 10}}, false);
    fRegistry.addClone("Pair/MFTMCHMID/uls/", "Pair/MFTMCHMID/lspp/");
    fRegistry.addClone("Pair/MFTMCHMID/uls/", "Pair/MFTMCHMID/lsmm/");
    fRegistry.addClone("Pair/MFTMCHMID/", "Pair/MCHMID/");
  }

  template <int mu_id, typename TTrack, typename TCollision>
  void fillTrackHistogram(TTrack const& track, TCollision const& collision)
  {
    o2::dataformats::GlobalFwdTrack propmuonAtPV = PropagateMuon(track, collision, skimmerPrimaryMuon::MuonExtrapolation::kToVertex); // this is for MCH-MID tracks that cannot see the primary vertex.
    o2::dataformats::GlobalFwdTrack propmuonAtDCA = PropagateMuon(track, collision, skimmerPrimaryMuon::MuonExtrapolation::kToDCA);
    o2::dataformats::GlobalFwdTrack propmuonAtRabs = PropagateMuon(track, collision, skimmerPrimaryMuon::MuonExtrapolation::kToRabs);

    float p = propmuonAtDCA.getP();
    float pt = propmuonAtDCA.getPt();
    float eta = propmuonAtDCA.getEta();
    float phi = propmuonAtDCA.getPhi();
    if (static_cast<int>(track.trackType()) > 2) { // only for MUON standalone
      p = propmuonAtPV.getP();
      pt = propmuonAtPV.getPt();
      eta = propmuonAtPV.getEta();
      phi = propmuonAtPV.getPhi();
    }

    o2::math_utils::bringTo02Pi(phi);
    if (phi < 0.f || 2.f * M_PI < phi) {
      return;
    }

    if (eta < mineta || maxeta < eta) {
      return;
    }

    if (track.trackType() == static_cast<uint8_t>(o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack)) {
      if (eta < mineta_mft || maxeta_mft < eta) {
        return;
      }
    }

    float rAtAbsorberEnd = track.rAtAbsorberEnd();
    if (static_cast<int>(track.trackType()) > 2) { // only for MUON standalone
      float xAbs = propmuonAtRabs.getX();
      float yAbs = propmuonAtRabs.getY();
      rAtAbsorberEnd = std::sqrt(xAbs * xAbs + yAbs * yAbs);
      // Redo propagation only for muon tracks // propagation of MFT tracks alredy done in reconstruction
    }

    if (rAtAbsorberEnd < minRabs || maxRabs < rAtAbsorberEnd) {
      return;
    }

    float dcaX = (propmuonAtDCA.getX() - collision.posX());
    float dcaY = (propmuonAtDCA.getY() - collision.posY());
    float dcaXY = std::sqrt(dcaX * dcaX + dcaY * dcaY);
    float cXX = propmuonAtDCA.getSigma2X();
    float cYY = propmuonAtDCA.getSigma2Y();
    float cXY = propmuonAtDCA.getSigmaXY();

    float dcaXYinSigma = 999.f;
    float det = cXX * cYY - cXY * cXY;
    if (det < 0) {
      dcaXYinSigma = 999.f;
    } else {
      float chi2 = (dcaX * dcaX * cYY + dcaY * dcaY * cXX - 2. * dcaX * dcaY * cXY) / det;
      dcaXYinSigma = std::sqrt(std::abs(chi2) / 2.); // in sigma
    }

    fRegistry.fill(HIST("Track/") + HIST(muon_types[mu_id]) + HIST("hPt"), pt);
    fRegistry.fill(HIST("Track/") + HIST(muon_types[mu_id]) + HIST("hEtaPhi"), phi, eta);
    fRegistry.fill(HIST("Track/") + HIST(muon_types[mu_id]) + HIST("hNclusters"), track.nClusters());
    fRegistry.fill(HIST("Track/") + HIST(muon_types[mu_id]) + HIST("hPDCA"), pt, track.pDca());
    fRegistry.fill(HIST("Track/") + HIST(muon_types[mu_id]) + HIST("hPDCA_recalc"), pt, p * dcaXY);
    fRegistry.fill(HIST("Track/") + HIST(muon_types[mu_id]) + HIST("hRatAbsorberEnd"), rAtAbsorberEnd);
    fRegistry.fill(HIST("Track/") + HIST(muon_types[mu_id]) + HIST("hChi2"), track.chi2());
    fRegistry.fill(HIST("Track/") + HIST(muon_types[mu_id]) + HIST("hChi2MatchMCHMID"), track.chi2MatchMCHMID());
    fRegistry.fill(HIST("Track/") + HIST(muon_types[mu_id]) + HIST("hChi2MatchMCHMFT"), track.chi2MatchMCHMFT());
    fRegistry.fill(HIST("Track/") + HIST(muon_types[mu_id]) + HIST("hMatchScoreMCHMFT"), track.matchScoreMCHMFT());
    fRegistry.fill(HIST("Track/") + HIST(muon_types[mu_id]) + HIST("hDCAxy2D"), dcaX, dcaY);
    fRegistry.fill(HIST("Track/") + HIST(muon_types[mu_id]) + HIST("hDCAxy2DinSigma"), dcaX / std::sqrt(cXX), dcaY / std::sqrt(cYY));
    fRegistry.fill(HIST("Track/") + HIST(muon_types[mu_id]) + HIST("hDCAxySigma"), dcaXYinSigma);
    fRegistry.fill(HIST("Track/") + HIST(muon_types[mu_id]) + HIST("hDCAxResolutionvsPt"), pt, std::sqrt(cXX) * 1e+4); // convert cm to um
    fRegistry.fill(HIST("Track/") + HIST(muon_types[mu_id]) + HIST("hDCAyResolutionvsPt"), pt, std::sqrt(cYY) * 1e+4); // convert cm to um
  }

  template <typename TTrack>
  bool isGlobalMuon(TTrack const& track)
  {
    if (track.trackType() != static_cast<uint8_t>(o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack)) {
      return false;
    }
    if (track.eta() < mineta_mft || maxeta_mft < track.eta()) {
      return false;
    }
    return true;
  }

  template <typename TTrack>
  bool isStandaloneMuon(TTrack const& track)
  {
    if (track.trackType() != static_cast<uint8_t>(o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack)) {
      return false;
    }
    return true;
  }

  template <typename T, typename C>
  o2::dataformats::GlobalFwdTrack PropagateMuon(T const& muon, C const& collision, const skimmerPrimaryMuon::MuonExtrapolation endPoint)
  {
    double chi2 = muon.chi2();
    SMatrix5 tpars(muon.x(), muon.y(), muon.phi(), muon.tgl(), muon.signed1Pt());
    std::vector<double> v1{muon.cXX(), muon.cXY(), muon.cYY(), muon.cPhiX(), muon.cPhiY(),
                           muon.cPhiPhi(), muon.cTglX(), muon.cTglY(), muon.cTglPhi(), muon.cTglTgl(),
                           muon.c1PtX(), muon.c1PtY(), muon.c1PtPhi(), muon.c1PtTgl(), muon.c1Pt21Pt2()};
    SMatrix55 tcovs(v1.begin(), v1.end());
    o2::track::TrackParCovFwd fwdtrack{muon.z(), tpars, tcovs, chi2};
    o2::dataformats::GlobalFwdTrack propmuon;

    if (static_cast<int>(muon.trackType()) > 2) { // MCH-MID or MCH standalone
      o2::dataformats::GlobalFwdTrack track;
      track.setParameters(tpars);
      track.setZ(fwdtrack.getZ());
      track.setCovariances(tcovs);
      auto mchTrack = mMatching.FwdtoMCH(track);

      if (endPoint == skimmerPrimaryMuon::MuonExtrapolation::kToVertex) {
        o2::mch::TrackExtrap::extrapToVertex(mchTrack, collision.posX(), collision.posY(), collision.posZ(), collision.covXX(), collision.covYY());
      }
      if (endPoint == skimmerPrimaryMuon::MuonExtrapolation::kToDCA) {
        o2::mch::TrackExtrap::extrapToVertexWithoutBranson(mchTrack, collision.posZ());
      }
      if (endPoint == skimmerPrimaryMuon::MuonExtrapolation::kToRabs) {
        o2::mch::TrackExtrap::extrapToZ(mchTrack, -505.);
      }

      auto proptrack = mMatching.MCHtoFwd(mchTrack);
      propmuon.setParameters(proptrack.getParameters());
      propmuon.setZ(proptrack.getZ());
      propmuon.setCovariances(proptrack.getCovariances());
    } else if (static_cast<int>(muon.trackType()) < 2) { // MFT-MCH-MID
      double centerMFT[3] = {0, 0, -61.4};
      o2::field::MagneticField* field = static_cast<o2::field::MagneticField*>(TGeoGlobalMagField::Instance()->GetField());
      auto Bz = field->getBz(centerMFT); // Get field at centre of MFT
      auto geoMan = o2::base::GeometryManager::meanMaterialBudget(muon.x(), muon.y(), muon.z(), collision.posX(), collision.posY(), collision.posZ());
      auto x2x0 = static_cast<float>(geoMan.meanX2X0);
      fwdtrack.propagateToVtxhelixWithMCS(collision.posZ(), {collision.posX(), collision.posY()}, {collision.covXX(), collision.covYY()}, Bz, x2x0);
      propmuon.setParameters(fwdtrack.getParameters());
      propmuon.setZ(fwdtrack.getZ());
      propmuon.setCovariances(fwdtrack.getCovariances());
    }
    return propmuon;
  }

  template <typename TMuon, typename TCollision>
  void fillMuonTable(TMuon const& muon, TCollision const& collision)
  {
    o2::dataformats::GlobalFwdTrack propmuonAtPV = PropagateMuon(muon, collision, skimmerPrimaryMuon::MuonExtrapolation::kToVertex);
    o2::dataformats::GlobalFwdTrack propmuonAtDCA = PropagateMuon(muon, collision, skimmerPrimaryMuon::MuonExtrapolation::kToDCA);
    o2::dataformats::GlobalFwdTrack propmuonAtRabs = PropagateMuon(muon, collision, skimmerPrimaryMuon::MuonExtrapolation::kToRabs);

    float pt = propmuonAtDCA.getPt();
    float eta = propmuonAtDCA.getEta();
    float phi = propmuonAtDCA.getPhi();
    if (static_cast<int>(muon.trackType()) > 2) { // only for MUON standalone
      pt = propmuonAtPV.getPt();
      eta = propmuonAtPV.getEta();
      phi = propmuonAtPV.getPhi();
    }

    if (eta < mineta || maxeta < eta) {
      return;
    }

    o2::math_utils::bringTo02Pi(phi);
    if (phi < 0.f || 2.f * M_PI < phi) {
      return;
    }

    if (muon.pDca() > maxPDCA) {
      return;
    }

    float dcaX = propmuonAtDCA.getX() - collision.posX();
    float dcaY = propmuonAtDCA.getY() - collision.posY();
    float rAtAbsorberEnd = muon.rAtAbsorberEnd();

    int nClustersMFT = 0;
    if (muon.trackType() == static_cast<uint8_t>(o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack)) {
      auto mftsa_track = muon.template matchMFTTrack_as<aod::MFTTracks>();
      nClustersMFT = mftsa_track.nClusters();
      if (eta < mineta_mft || maxeta_mft < eta) {
        return;
      }
    } else if (muon.trackType() == static_cast<uint8_t>(o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack)) {
      nClustersMFT = 0;
      float xAbs = propmuonAtRabs.getX();
      float yAbs = propmuonAtRabs.getY();
      rAtAbsorberEnd = std::sqrt(xAbs * xAbs + yAbs * yAbs); // Redo propagation only for muon tracks // propagation of MFT tracks alredy done in reconstruction
    }

    if (rAtAbsorberEnd < minRabs || maxRabs < rAtAbsorberEnd) {
      return;
    }

    emprimarymuons(muon.collisionId(), muon.globalIndex(), muon.trackType(), pt, eta, phi, muon.sign(), dcaX, dcaY,
                   propmuonAtDCA.getX(), propmuonAtDCA.getY(), propmuonAtDCA.getZ(), propmuonAtDCA.getTgl(),
                   muon.nClusters(), muon.pDca(), rAtAbsorberEnd, muon.chi2(), muon.chi2MatchMCHMID(), muon.chi2MatchMCHMFT(),
                   muon.matchScoreMCHMFT(), muon.matchMFTTrackId(), muon.matchMCHTrackId(), muon.mchBitMap(), muon.midBitMap(), muon.midBoards(), nClustersMFT);

    auto fwdcov = propmuonAtDCA.getCovariances();
    // LOGF(info, "fwdcov(0,0) = %f, propmuonAtDCA.getSigma2X() = %f", fwdcov(0,0), propmuonAtDCA.getSigma2X());

    // See definition DataFormats/Reconstruction/include/ReconstructionDataFormats/TrackFwd.h
    /// Covariance matrix of track parameters, ordered as follows:    <pre>
    ///  <X,X>         <Y,X>           <PHI,X>       <TANL,X>        <INVQPT,X>
    ///  <X,Y>         <Y,Y>           <PHI,Y>       <TANL,Y>        <INVQPT,Y>
    /// <X,PHI>       <Y,PHI>         <PHI,PHI>     <TANL,PHI>      <INVQPT,PHI>
    /// <X,TANL>      <Y,TANL>       <PHI,TANL>     <TANL,TANL>     <INVQPT,TANL>
    /// <X,INVQPT>   <Y,INVQPT>     <PHI,INVQPT>   <TANL,INVQPT>   <INVQPT,INVQPT>  </pre
    emprimarymuonscov(
      fwdcov(0, 0),
      fwdcov(0, 1),
      fwdcov(1, 1),
      fwdcov(2, 0),
      fwdcov(2, 1),
      fwdcov(2, 2),
      fwdcov(3, 0),
      fwdcov(3, 1),
      fwdcov(3, 2),
      fwdcov(3, 3),
      fwdcov(4, 0),
      fwdcov(4, 1),
      fwdcov(4, 2),
      fwdcov(4, 3),
      fwdcov(4, 4));
  }

  // Filter muonTypeFilter = o2::aod::fwdtrack::trackType == static_cast<uint8_t>(o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack) || o2::aod::fwdtrack::trackType == static_cast<uint8_t>(o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalForwardTrack) || o2::aod::fwdtrack::trackType == static_cast<uint8_t>(o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack);
  Filter muonTypeFilter = o2::aod::fwdtrack::trackType == static_cast<uint8_t>(o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack) || o2::aod::fwdtrack::trackType == static_cast<uint8_t>(o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack);
  using MyFilteredTracks = soa::Filtered<MyTracks>;

  // Partition<MyFilteredTracks> posTracks = o2::aod::fwdtrack::signed1Pt > 0.f;
  // Partition<MyFilteredTracks> negTracks = o2::aod::fwdtrack::signed1Pt < 0.f;

  void processRec(aod::Collisions const& collisions, aod::BCsWithTimestamps const&, MyFilteredTracks const& tracks, aod::MFTTracks const& mfttracks)
  {
    for (auto& collision : collisions) {
      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);
      fRegistry.fill(HIST("Event/hCollisionCounter"), 0.f);

      auto tracks_per_coll = tracks.sliceBy(perCollision, collision.globalIndex());
      // auto mfttracks_per_coll = mfttracks.sliceBy(perCollision_mft, collision.globalIndex());
      // LOGF(info, "tracks_per_coll.size() = %d, mfttracks_per_coll.size() = %d", tracks_per_coll.size(), mfttracks_per_coll.size());

      for (auto& track : tracks_per_coll) {
        fRegistry.fill(HIST("Track/hMuonType"), track.trackType());
        // LOGF(info, "track.globalIndex() = %d, track.trackType() = %d, track.matchMFTTrackId() = %d, track.matchMCHTrackId() = %d, track.offsets() = %d", track.globalIndex(), track.trackType(), track.matchMFTTrackId(), track.matchMCHTrackId(), track.offsets());

        if (track.trackType() == static_cast<uint8_t>(o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack)) {
          fillTrackHistogram<0>(track, collision);
          auto mftsa_track = track.template matchMFTTrack_as<aod::MFTTracks>();
          fRegistry.fill(HIST("Track/MFTMCHMID/hNclustersMFT"), mftsa_track.nClusters());
          fillMuonTable(track, collision);
        } else if (track.trackType() == static_cast<uint8_t>(o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack)) {
          fillTrackHistogram<3>(track, collision);
          fRegistry.fill(HIST("Track/MCHMID/hNclustersMFT"), 0);
          fillMuonTable(track, collision);
        }
      } // end of track loop

      // auto posTracks_per_coll = posTracks->sliceByCached(o2::aod::fwdtrack::collisionId, collision.globalIndex(), cache);
      // auto negTracks_per_coll = negTracks->sliceByCached(o2::aod::fwdtrack::collisionId, collision.globalIndex(), cache);
      // fRegistry.fill(HIST("Event/hNmuon"), negTracks_per_coll.size(), posTracks_per_coll.size());
      //// LOGF(info, "collision.globalIndex() = %d, negTracks_per_coll.size() = %d, posTracks_per_coll.size() = %d", collision.globalIndex(), negTracks_per_coll.size(), posTracks_per_coll.size());

      // for (auto& [pos, neg] : combinations(CombinationsFullIndexPolicy(posTracks_per_coll, negTracks_per_coll))) { // ULS
      //   o2::dataformats::GlobalFwdTrack propmuonAtPV_pos = PropagateMuon(pos, collision, skimmerPrimaryMuon::MuonExtrapolation::kToVertex);
      //   o2::dataformats::GlobalFwdTrack propmuonAtPV_neg = PropagateMuon(neg, collision, skimmerPrimaryMuon::MuonExtrapolation::kToVertex);

      //  ROOT::Math::PtEtaPhiMVector v1(propmuonAtPV_pos.getPt(), propmuonAtPV_pos.getEta(), propmuonAtPV_pos.getPhi(), o2::constants::physics::MassMuon);
      //  ROOT::Math::PtEtaPhiMVector v2(propmuonAtPV_neg.getPt(), propmuonAtPV_neg.getEta(), propmuonAtPV_neg.getPhi(), o2::constants::physics::MassMuon);
      //  ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;

      //  if (v12.Rapidity() < mineta || maxeta < v12.Rapidity()) {
      //    continue;
      //  }
      //  if (isGlobalMuon(pos) && isGlobalMuon(neg) && (mineta_mft < v12.Rapidity() || v12.Rapidity() < maxeta_mft)) {
      //    fRegistry.fill(HIST("Pair/MFTMCHMID/uls/hMvsPt"), v12.M(), v12.Pt());
      //  }
      //  if (isStandaloneMuon(pos) && isStandaloneMuon(neg) && (mineta < v12.Rapidity() || v12.Rapidity() < maxeta)) {
      //    fRegistry.fill(HIST("Pair/MCHMID/uls/hMvsPt"), v12.M(), v12.Pt());
      //  }
      //}

      // for (auto& [pos1, pos2] : combinations(CombinationsStrictlyUpperIndexPolicy(posTracks_per_coll, posTracks_per_coll))) { // LS++
      //   o2::dataformats::GlobalFwdTrack propmuonAtPV_pos1 = PropagateMuon(pos1, collision, skimmerPrimaryMuon::MuonExtrapolation::kToVertex);
      //   o2::dataformats::GlobalFwdTrack propmuonAtPV_pos2 = PropagateMuon(pos2, collision, skimmerPrimaryMuon::MuonExtrapolation::kToVertex);
      //   ROOT::Math::PtEtaPhiMVector v1(propmuonAtPV_pos1.getPt(), propmuonAtPV_pos1.getEta(), propmuonAtPV_pos1.getPhi(), o2::constants::physics::MassMuon);
      //   ROOT::Math::PtEtaPhiMVector v2(propmuonAtPV_pos2.getPt(), propmuonAtPV_pos2.getEta(), propmuonAtPV_pos2.getPhi(), o2::constants::physics::MassMuon);
      //   ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;

      //  if (v12.Rapidity() < mineta || maxeta < v12.Rapidity()) {
      //    continue;
      //  }
      //  if (isGlobalMuon(pos1) && isGlobalMuon(pos2) && (mineta_mft < v12.Rapidity() || v12.Rapidity() < maxeta_mft)) {
      //    fRegistry.fill(HIST("Pair/MFTMCHMID/lspp/hMvsPt"), v12.M(), v12.Pt());
      //  }
      //  if (isStandaloneMuon(pos1) && isStandaloneMuon(pos2) && (mineta < v12.Rapidity() || v12.Rapidity() < maxeta)) {
      //    fRegistry.fill(HIST("Pair/MCHMID/lspp/hMvsPt"), v12.M(), v12.Pt());
      //  }
      //}

      // for (auto& [neg1, neg2] : combinations(CombinationsStrictlyUpperIndexPolicy(negTracks_per_coll, negTracks_per_coll))) { // LS--
      //   o2::dataformats::GlobalFwdTrack propmuonAtPV_neg1 = PropagateMuon(neg1, collision, skimmerPrimaryMuon::MuonExtrapolation::kToVertex);
      //   o2::dataformats::GlobalFwdTrack propmuonAtPV_neg2 = PropagateMuon(neg2, collision, skimmerPrimaryMuon::MuonExtrapolation::kToVertex);
      //   ROOT::Math::PtEtaPhiMVector v1(propmuonAtPV_neg1.getPt(), propmuonAtPV_neg1.getEta(), propmuonAtPV_neg1.getPhi(), o2::constants::physics::MassMuon);
      //   ROOT::Math::PtEtaPhiMVector v2(propmuonAtPV_neg2.getPt(), propmuonAtPV_neg2.getEta(), propmuonAtPV_neg2.getPhi(), o2::constants::physics::MassMuon);
      //   ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;

      //  if (v12.Rapidity() < mineta || maxeta < v12.Rapidity()) {
      //    continue;
      //  }
      //  if (isGlobalMuon(neg1) && isGlobalMuon(neg2) && (mineta_mft < v12.Rapidity() || v12.Rapidity() < maxeta_mft)) {
      //    fRegistry.fill(HIST("Pair/MFTMCHMID/lsmm/hMvsPt"), v12.M(), v12.Pt());
      //  }
      //  if (isStandaloneMuon(neg1) && isStandaloneMuon(neg2) && (mineta < v12.Rapidity() || v12.Rapidity() < maxeta)) {
      //    fRegistry.fill(HIST("Pair/MCHMID/lsmm/hMvsPt"), v12.M(), v12.Pt());
      //  }
      //}

    } // end of collision loop
  }
  PROCESS_SWITCH(skimmerPrimaryMuon, processRec, "process reconstructed info only", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<skimmerPrimaryMuon>(cfgc, TaskName{"skimmer-primary-muon"})};
}
