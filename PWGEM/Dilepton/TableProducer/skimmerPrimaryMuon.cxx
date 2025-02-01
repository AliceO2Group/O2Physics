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

using MyCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::EMEvSels>;
using MyCollisionsWithSWT = soa::Join<MyCollisions, aod::EMSWTriggerInfosTMP>;

using MyTracks = soa::Join<aod::FwdTracks, aod::FwdTracksCov>; // muon tracks are repeated. i.e. not exclusive.
using MyTrack = MyTracks::iterator;

using MyTracksMC = soa::Join<MyTracks, aod::McFwdTrackLabels>;
using MyTrackMC = MyTracksMC::iterator;

using MFTTracksMC = soa::Join<o2::aod::MFTTracks, aod::McMFTTrackLabels>;
using MFTTrackMC = MFTTracksMC::iterator;

struct skimmerPrimaryMuon {
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
  Configurable<bool> fillQAHistogram{"fillQAHistogram", false, "flag to fill QA histograms"};
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
  Configurable<float> minpt{"minpt", 0.2, "min pt for muon"};
  Configurable<float> mineta{"mineta", -4.0, "eta acceptance"};
  Configurable<float> maxeta{"maxeta", -2.5, "eta acceptance"};
  Configurable<float> mineta_mft{"mineta_mft", -3.6, "eta acceptance"};
  Configurable<float> maxeta_mft{"maxeta_mft", -2.5, "eta acceptance"};
  Configurable<float> minRabs{"minRabs", 17.6, "min. R at absorber end"};
  Configurable<float> maxRabs{"maxRabs", 89.5, "max. R at absorber end"};

  o2::ccdb::CcdbApi ccdbApi;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  int mRunNumber;

  o2::globaltracking::MatchGlobalFwd mMatching;
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
    if (fillQAHistogram) {
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
      fRegistry.add("Track/MFTMCHMID/hRatAbsorberEnd", "R at absorber end;R at absorber end (cm)", kTH1F, {{200, 0.0f, 200}}, false);
      fRegistry.add("Track/MFTMCHMID/hPDCA", "pDCA;p_{T,#mu} at PV (GeV/c);p #times DCA (GeV/c #upoint cm)", kTH2F, {{100, 0, 10}, {100, 0.0f, 1000}}, false);
      fRegistry.add("Track/MFTMCHMID/hPDCA_recalc", "pDCA relcalculated;p_{T,#mu} at PV (GeV/c);p #times DCA (GeV/c #upoint cm)", kTH2F, {{100, 0, 10}, {100, 0.0f, 1000}}, false);
      fRegistry.add("Track/MFTMCHMID/hChi2", "chi2;chi2", kTH1F, {{100, 0.0f, 100}}, false);
      fRegistry.add("Track/MFTMCHMID/hChi2MatchMCHMID", "chi2 match MCH-MID;chi2", kTH1F, {{100, 0.0f, 100}}, false);
      fRegistry.add("Track/MFTMCHMID/hChi2MatchMCHMFT", "chi2 match MCH-MFT;chi2", kTH1F, {{100, 0.0f, 100}}, false);
      fRegistry.add("Track/MFTMCHMID/hMatchScoreMCHMFT", "match score MCH-MFT;score", kTH1F, {{100, 0.0f, 100}}, false);
      fRegistry.add("Track/MFTMCHMID/hDCAxy2D", "DCA xy;DCA_{x} (cm);DCA_{y} (cm)", kTH2F, {{200, -10, 10}, {200, -10, +10}}, false);
      fRegistry.add("Track/MFTMCHMID/hDCAxy2DinSigma", "DCA xy;DCA_{x} (#sigma);DCA_{y} (#sigma)", kTH2F, {{200, -10, 10}, {200, -10, +10}}, false);
      fRegistry.add("Track/MFTMCHMID/hDCAxySigma", "DCA_{xy};DCA_{xy} (#sigma);", kTH1F, {{100, 0, 10}}, false);
      fRegistry.add("Track/MFTMCHMID/hDCAxResolutionvsPt", "DCA_{x} vs. p_{T,#mu};p_{T,#mu} (GeV/c);DCA_{x} resolution (#mum);", kTH2F, {{100, 0, 10.f}, {100, 0, 100}}, false);
      fRegistry.add("Track/MFTMCHMID/hDCAyResolutionvsPt", "DCA_{y} vs. p_{T,#mu};p_{T,#mu} (GeV/c);DCA_{y} resolution (#mum);", kTH2F, {{100, 0, 10.f}, {100, 0, 100}}, false);
      fRegistry.add("Track/MFTMCHMID/hChi2MatchMCHMFT_DCAxySigma", "chi2 match MCH-MFT;DCA_{xy,#mu} (#sigma);MFT-MCH matching chi2", kTH2F, {{100, 0, 10}, {100, 0.0f, 100}}, false);
      fRegistry.add("Track/MFTMCHMID/hChi2MatchMCHMFT_Pt", "chi2 match MCH-MFT;p_{T,#mu} (GeV/c);MFT-MCH matching chi2", kTH2F, {{100, 0, 10}, {100, 0.0f, 100}}, false);
      fRegistry.addClone("Track/MFTMCHMID/", "Track/MCHMID/");
    }
  }

  template <int mu_id, typename TTrack, typename TCollision>
  void fillTrackHistogram(TTrack const& track, TCollision const& collision)
  {
    o2::dataformats::GlobalFwdTrack propmuonAtPV = PropagateMuon(track, collision, skimmerPrimaryMuon::MuonExtrapolation::kToVertex); // this is for MCH-MID tracks that cannot see the primary vertex.
    o2::dataformats::GlobalFwdTrack propmuonAtDCA = PropagateMuon(track, collision, skimmerPrimaryMuon::MuonExtrapolation::kToDCA);

    float p = propmuonAtPV.getP();
    float pt = propmuonAtPV.getPt();
    float eta = propmuonAtPV.getEta();
    float phi = propmuonAtPV.getPhi();

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
      o2::dataformats::GlobalFwdTrack propmuonAtRabs = PropagateMuon(track, collision, skimmerPrimaryMuon::MuonExtrapolation::kToRabs);
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
      dcaXYinSigma = std::sqrt(std::fabs(chi2) / 2.); // in sigma
    }

    fRegistry.fill(HIST("Track/hMuonType"), track.trackType());
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
    fRegistry.fill(HIST("Track/") + HIST(muon_types[mu_id]) + HIST("hChi2MatchMCHMFT_DCAxySigma"), dcaXYinSigma, track.chi2MatchMCHMFT());
    fRegistry.fill(HIST("Track/") + HIST(muon_types[mu_id]) + HIST("hChi2MatchMCHMFT_Pt"), pt, track.chi2MatchMCHMFT());
    fRegistry.fill(HIST("Track/") + HIST(muon_types[mu_id]) + HIST("hDCAxResolutionvsPt"), pt, std::sqrt(cXX) * 1e+4); // convert cm to um
    fRegistry.fill(HIST("Track/") + HIST(muon_types[mu_id]) + HIST("hDCAyResolutionvsPt"), pt, std::sqrt(cYY) * 1e+4); // convert cm to um
  }

  template <typename T, typename C>
  o2::dataformats::GlobalFwdTrack PropagateMuon(T const& muon, C const& collision, const skimmerPrimaryMuon::MuonExtrapolation endPoint)
  {
    double chi2 = muon.chi2();
    SMatrix5 tpars(muon.x(), muon.y(), muon.phi(), muon.tgl(), muon.signed1Pt());
    std::vector<float> v1{muon.cXX(), muon.cXY(), muon.cYY(), muon.cPhiX(), muon.cPhiY(),
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

    v1.clear();
    v1.shrink_to_fit();

    return propmuon;
  }

  template <bool isMC, typename TMuon, typename TCollision>
  bool isSelected(TMuon const& muon, TCollision const& collision)
  {
    if constexpr (isMC) {
      if (!muon.has_mcParticle()) {
        return false;
      }
    }

    o2::dataformats::GlobalFwdTrack propmuonAtPV = PropagateMuon(muon, collision, skimmerPrimaryMuon::MuonExtrapolation::kToVertex);

    float pt = propmuonAtPV.getPt();
    float eta = propmuonAtPV.getEta();
    float phi = propmuonAtPV.getPhi();

    if (pt < minpt) {
      return false;
    }

    if (eta < mineta || maxeta < eta) {
      return false;
    }

    o2::math_utils::bringTo02Pi(phi);
    if (phi < 0.f || 2.f * M_PI < phi) {
      return false;
    }

    // float dcaX = propmuonAtDCA.getX() - collision.posX();
    // float dcaY = propmuonAtDCA.getY() - collision.posY();
    float rAtAbsorberEnd = muon.rAtAbsorberEnd();

    if (muon.trackType() == static_cast<uint8_t>(o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack)) {
      if (eta < mineta_mft || maxeta_mft < eta) {
        return false;
      }
    } else if (muon.trackType() == static_cast<uint8_t>(o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack)) {
      o2::dataformats::GlobalFwdTrack propmuonAtRabs = PropagateMuon(muon, collision, skimmerPrimaryMuon::MuonExtrapolation::kToRabs);
      float xAbs = propmuonAtRabs.getX();
      float yAbs = propmuonAtRabs.getY();
      rAtAbsorberEnd = std::sqrt(xAbs * xAbs + yAbs * yAbs); // Redo propagation only for muon tracks // propagation of MFT tracks alredy done in reconstruction
    } else {
      return false;
    }

    if (rAtAbsorberEnd < minRabs || maxRabs < rAtAbsorberEnd) {
      return false;
    }
    return true;
  }

  template <typename TMFTsaTracks, typename TMuon, typename TCollision>
  void fillMuonTable(TMuon const& muon, TCollision const& collision, const int new_muon_sa_id)
  {
    o2::dataformats::GlobalFwdTrack propmuonAtPV = PropagateMuon(muon, collision, skimmerPrimaryMuon::MuonExtrapolation::kToVertex);
    o2::dataformats::GlobalFwdTrack propmuonAtDCA = PropagateMuon(muon, collision, skimmerPrimaryMuon::MuonExtrapolation::kToDCA);

    float pt = propmuonAtPV.getPt();
    float eta = propmuonAtPV.getEta();
    float phi = propmuonAtPV.getPhi();

    o2::math_utils::bringTo02Pi(phi);

    float dcaX = propmuonAtDCA.getX() - collision.posX();
    float dcaY = propmuonAtDCA.getY() - collision.posY();
    float rAtAbsorberEnd = muon.rAtAbsorberEnd();

    if (muon.trackType() == static_cast<uint8_t>(o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack)) {
      o2::dataformats::GlobalFwdTrack propmuonAtRabs = PropagateMuon(muon, collision, skimmerPrimaryMuon::MuonExtrapolation::kToRabs);
      float xAbs = propmuonAtRabs.getX();
      float yAbs = propmuonAtRabs.getY();
      rAtAbsorberEnd = std::sqrt(xAbs * xAbs + yAbs * yAbs); // Redo propagation only for muon tracks // propagation of MFT tracks alredy done in reconstruction
    }

    bool isAssociatedToMPC = collision.globalIndex() == muon.collisionId();
    auto fwdcov = propmuonAtDCA.getCovariances();
    if (muon.trackType() == static_cast<uint8_t>(o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack)) {
      auto mftsa_track = muon.template matchMFTTrack_as<TMFTsaTracks>();

      emprimarymuons(collision.globalIndex(), muon.globalIndex(), muon.trackType(), pt, eta, phi, muon.sign(), dcaX, dcaY,
                     propmuonAtDCA.getX(), propmuonAtDCA.getY(), propmuonAtDCA.getZ(), propmuonAtDCA.getTgl(),
                     muon.nClusters(), muon.pDca(), rAtAbsorberEnd, muon.chi2(), muon.chi2MatchMCHMID(), muon.chi2MatchMCHMFT(),
                     new_muon_sa_id,
                     muon.mchBitMap(), muon.midBitMap(), muon.midBoards(), mftsa_track.mftClusterSizesAndTrackFlags(), mftsa_track.chi2(), isAssociatedToMPC);
      emprimarymuonscov(
        fwdcov(0, 0),
        fwdcov(0, 1), fwdcov(1, 1),
        fwdcov(2, 0), fwdcov(2, 1), fwdcov(2, 2),
        fwdcov(3, 0), fwdcov(3, 1), fwdcov(3, 2), fwdcov(3, 3),
        fwdcov(4, 0), fwdcov(4, 1), fwdcov(4, 2), fwdcov(4, 3), fwdcov(4, 4));
    } else if (muon.trackType() == static_cast<uint8_t>(o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack)) {
      emprimarymuons(collision.globalIndex(), muon.globalIndex(), muon.trackType(), pt, eta, phi, muon.sign(), dcaX, dcaY,
                     propmuonAtDCA.getX(), propmuonAtDCA.getY(), propmuonAtDCA.getZ(), propmuonAtDCA.getTgl(),
                     muon.nClusters(), muon.pDca(), rAtAbsorberEnd, muon.chi2(), muon.chi2MatchMCHMID(), muon.chi2MatchMCHMFT(),
                     new_muon_sa_id,
                     muon.mchBitMap(), muon.midBitMap(), muon.midBoards(), 0, 999999.f, isAssociatedToMPC);
      emprimarymuonscov(
        fwdcov(0, 0),
        fwdcov(0, 1), fwdcov(1, 1),
        fwdcov(2, 0), fwdcov(2, 1), fwdcov(2, 2),
        fwdcov(3, 0), fwdcov(3, 1), fwdcov(3, 2), fwdcov(3, 3),
        fwdcov(4, 0), fwdcov(4, 1), fwdcov(4, 2), fwdcov(4, 3), fwdcov(4, 4));
    }

    // See definition DataFormats/Reconstruction/include/ReconstructionDataFormats/TrackFwd.h
    /// Covariance matrix of track parameters, ordered as follows:    <pre>
    ///  <X,X>         <Y,X>           <PHI,X>       <TANL,X>        <INVQPT,X>
    ///  <X,Y>         <Y,Y>           <PHI,Y>       <TANL,Y>        <INVQPT,Y>
    /// <X,PHI>       <Y,PHI>         <PHI,PHI>     <TANL,PHI>      <INVQPT,PHI>
    /// <X,TANL>      <Y,TANL>       <PHI,TANL>     <TANL,TANL>     <INVQPT,TANL>
    /// <X,INVQPT>   <Y,INVQPT>     <PHI,INVQPT>   <TANL,INVQPT>   <INVQPT,INVQPT>  </pre
  }

  std::map<std::pair<int, int>, int> map_new_sa_muon_index; // new standalone muon index

  // Preslice<aod::FwdTracks> fwdtrackIndicesPerMFTsa = aod::fwdtrack::matchMFTTrackId;

  Preslice<aod::FwdTrackAssoc> fwdtrackIndicesPerCollision = aod::track_association::collisionId;
  Partition<MyTracks> global_muons = o2::aod::fwdtrack::trackType == uint8_t(o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack); // MFT-MCH-MID
  Partition<MyTracks> sa_muons = o2::aod::fwdtrack::trackType == uint8_t(o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack); // MCH-MID

  void processRec_SA(MyCollisions const& collisions, aod::BCsWithTimestamps const&, MyTracks const& tracks, aod::MFTTracks const&)
  {
    for (auto& collision : collisions) {
      auto bc = collision.template foundBC_as<aod::BCsWithTimestamps>();
      initCCDB(bc);
      fRegistry.fill(HIST("Event/hCollisionCounter"), 0.f);
      if (!collision.isSelected()) {
        continue;
      }

      auto sa_muons_per_coll = sa_muons->sliceByCached(o2::aod::fwdtrack::collisionId, collision.globalIndex(), cache);
      auto global_muons_per_coll = global_muons->sliceByCached(o2::aod::fwdtrack::collisionId, collision.globalIndex(), cache);

      int counter = 0;
      int offset = emprimarymuons.lastIndex() + 1;

      for (auto& track : sa_muons_per_coll) {
        if (fillQAHistogram) {
          fillTrackHistogram<3>(track, collision);
        }
        if (isSelected<false>(track, collision)) {
          map_new_sa_muon_index[std::make_pair(collision.globalIndex(), track.globalIndex())] = counter + offset;
          counter++;
        }
      } // end of standalone muon loop
      for (auto& track : global_muons_per_coll) {
        if (fillQAHistogram) {
          fillTrackHistogram<0>(track, collision);
        }

        if (map_new_sa_muon_index.find(std::make_pair(collision.globalIndex(), track.matchMCHTrackId())) == map_new_sa_muon_index.end()) { // don't apply muon selection to MCH-MID track in MFT-MCH-MID track
          map_new_sa_muon_index[std::make_pair(collision.globalIndex(), track.matchMCHTrackId())] = counter + offset;
          counter++;
        }
        if (isSelected<false>(track, collision)) {
          map_new_sa_muon_index[std::make_pair(collision.globalIndex(), track.globalIndex())] = counter + offset;
          counter++;
        }
      } // end of global muon loop

      // fill table after mapping
      for (const auto& [key, value] : map_new_sa_muon_index) {
        // int collisionId = std::get<0>(key);
        // int fwdtrackId = std::get<1>(key);
        // int new_fwdtrackId = value;
        // LOGF(info, "collisionId = %d, fwdtrackId = %d, new_fwdtrackId = %d", collisionId, fwdtrackId, new_fwdtrackId);
        auto track = tracks.iteratorAt(std::get<1>(key));
        if (track.trackType() == static_cast<uint8_t>(o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack)) {
          fillMuonTable<aod::MFTTracks>(track, collision, value);
        } else if (track.trackType() == static_cast<uint8_t>(o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack)) {
          fillMuonTable<aod::MFTTracks>(track, collision, map_new_sa_muon_index[std::make_pair(collision.globalIndex(), track.matchMCHTrackId())]);
        }
      }

      map_new_sa_muon_index.clear();
    } // end of collision loop
  }
  PROCESS_SWITCH(skimmerPrimaryMuon, processRec_SA, "process reconstructed info only with standalone", true);

  void processRec_TTCA(MyCollisions const& collisions, aod::BCsWithTimestamps const&, MyTracks const& tracks, aod::MFTTracks const&, aod::FwdTrackAssoc const& fwdtrackIndices)
  {
    for (auto& collision : collisions) {
      auto bc = collision.template foundBC_as<aod::BCsWithTimestamps>();
      initCCDB(bc);
      fRegistry.fill(HIST("Event/hCollisionCounter"), 0.f);
      if (!collision.isSelected()) {
        continue;
      }

      int counter = 0;
      int offset = emprimarymuons.lastIndex() + 1;

      auto fwdtrackIdsThisCollision = fwdtrackIndices.sliceBy(fwdtrackIndicesPerCollision, collision.globalIndex());
      for (auto& fwdtrackId : fwdtrackIdsThisCollision) {
        auto track = fwdtrackId.template fwdtrack_as<MyTracks>();
        // LOGF(info, "TTCA | collision.globalIndex() = %d, track.globalIndex() = %d, track.trackType() = %d, track.matchMFTTrackId() = %d, track.matchMCHTrackId() = %d, track.offsets() = %d", collision.globalIndex(), track.globalIndex(), track.trackType(), track.matchMFTTrackId(), track.matchMCHTrackId(), track.offsets());

        // auto collision_in_track = track.collision_as<MyCollisions>();
        // auto bc_in_track = collision_in_track.bc_as<aod::BCsWithTimestamps>();
        // LOGF(info, "track.globalIndex() = %d , bc_in_track.globalBC() = %lld", track.globalIndex(), bc_in_track.globalBC());

        if (track.trackType() == static_cast<uint8_t>(o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack)) {
          if (fillQAHistogram) {
            fillTrackHistogram<3>(track, collision);
          }
          if (isSelected<false>(track, collision)) {
            map_new_sa_muon_index[std::make_pair(collision.globalIndex(), track.globalIndex())] = counter + offset;
            counter++;
          }
        } else if (track.trackType() == static_cast<uint8_t>(o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack)) {
          if (fillQAHistogram) {
            fillTrackHistogram<0>(track, collision);
          }
          if (map_new_sa_muon_index.find(std::make_pair(collision.globalIndex(), track.matchMCHTrackId())) == map_new_sa_muon_index.end()) {
            map_new_sa_muon_index[std::make_pair(collision.globalIndex(), track.matchMCHTrackId())] = counter + offset;
            counter++;
          }
          if (isSelected<false>(track, collision)) {
            map_new_sa_muon_index[std::make_pair(collision.globalIndex(), track.globalIndex())] = counter + offset;
            counter++;
          }
        }
      } // end of track loop

      for (const auto& [key, value] : map_new_sa_muon_index) {
        // int collisionId = std::get<0>(key);
        // int fwdtrackId = std::get<1>(key);
        // int new_fwdtrackId = value;
        auto track = tracks.iteratorAt(std::get<1>(key));
        if (track.trackType() == static_cast<uint8_t>(o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack)) {
          fillMuonTable<aod::MFTTracks>(track, collision, value);
        } else if (track.trackType() == static_cast<uint8_t>(o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack)) {
          fillMuonTable<aod::MFTTracks>(track, collision, map_new_sa_muon_index[std::make_pair(collision.globalIndex(), track.matchMCHTrackId())]);
        }
      }

      map_new_sa_muon_index.clear();
    } // end of collision loop
  }
  PROCESS_SWITCH(skimmerPrimaryMuon, processRec_TTCA, "process reconstructed info only with TTCA", false);

  void processRec_SA_SWT(MyCollisionsWithSWT const& collisions, aod::BCsWithTimestamps const&, MyTracks const& tracks, aod::MFTTracks const&)
  {
    for (auto& collision : collisions) {
      auto bc = collision.template foundBC_as<aod::BCsWithTimestamps>();
      initCCDB(bc);
      fRegistry.fill(HIST("Event/hCollisionCounter"), 0.f);
      if (!collision.isSelected()) {
        continue;
      }
      if (collision.swtaliastmp_raw() == 0) {
        continue;
      }

      auto sa_muons_per_coll = sa_muons->sliceByCached(o2::aod::fwdtrack::collisionId, collision.globalIndex(), cache);
      auto global_muons_per_coll = global_muons->sliceByCached(o2::aod::fwdtrack::collisionId, collision.globalIndex(), cache);

      int counter = 0;
      int offset = emprimarymuons.lastIndex() + 1;

      for (auto& track : sa_muons_per_coll) {
        if (fillQAHistogram) {
          fillTrackHistogram<3>(track, collision);
        }
        if (isSelected<false>(track, collision)) {
          map_new_sa_muon_index[std::make_pair(collision.globalIndex(), track.globalIndex())] = counter + offset;
          counter++;
        }
      } // end of standalone muon loop
      for (auto& track : global_muons_per_coll) {
        if (fillQAHistogram) {
          fillTrackHistogram<0>(track, collision);
        }

        if (map_new_sa_muon_index.find(std::make_pair(collision.globalIndex(), track.matchMCHTrackId())) == map_new_sa_muon_index.end()) { // don't apply muon selection to MCH-MID track in MFT-MCH-MID track
          map_new_sa_muon_index[std::make_pair(collision.globalIndex(), track.matchMCHTrackId())] = counter + offset;
          counter++;
        }
        if (isSelected<false>(track, collision)) {
          map_new_sa_muon_index[std::make_pair(collision.globalIndex(), track.globalIndex())] = counter + offset;
          counter++;
        }
      } // end of global muon loop

      // fill table after mapping
      for (const auto& [key, value] : map_new_sa_muon_index) {
        // int collisionId = std::get<0>(key);
        // int fwdtrackId = std::get<1>(key);
        // int new_fwdtrackId = value;
        // LOGF(info, "collisionId = %d, fwdtrackId = %d, new_fwdtrackId = %d", collisionId, fwdtrackId, new_fwdtrackId);
        auto track = tracks.iteratorAt(std::get<1>(key));
        if (track.trackType() == static_cast<uint8_t>(o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack)) {
          fillMuonTable<aod::MFTTracks>(track, collision, value);
        } else if (track.trackType() == static_cast<uint8_t>(o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack)) {
          fillMuonTable<aod::MFTTracks>(track, collision, map_new_sa_muon_index[std::make_pair(collision.globalIndex(), track.matchMCHTrackId())]);
        }
      }

      map_new_sa_muon_index.clear();
    } // end of collision loop
  }
  PROCESS_SWITCH(skimmerPrimaryMuon, processRec_SA_SWT, "process reconstructed info only with standalone", false);

  void processRec_TTCA_SWT(MyCollisionsWithSWT const& collisions, aod::BCsWithTimestamps const&, MyTracks const& tracks, aod::MFTTracks const&, aod::FwdTrackAssoc const& fwdtrackIndices)
  {
    for (auto& collision : collisions) {
      auto bc = collision.template foundBC_as<aod::BCsWithTimestamps>();
      initCCDB(bc);
      fRegistry.fill(HIST("Event/hCollisionCounter"), 0.f);
      if (!collision.isSelected()) {
        continue;
      }
      if (collision.swtaliastmp_raw() == 0) {
        continue;
      }

      int counter = 0;
      int offset = emprimarymuons.lastIndex() + 1;

      auto fwdtrackIdsThisCollision = fwdtrackIndices.sliceBy(fwdtrackIndicesPerCollision, collision.globalIndex());
      for (auto& fwdtrackId : fwdtrackIdsThisCollision) {
        auto track = fwdtrackId.template fwdtrack_as<MyTracks>();
        // LOGF(info, "TTCA | collision.globalIndex() = %d, track.globalIndex() = %d, track.trackType() = %d, track.matchMFTTrackId() = %d, track.matchMCHTrackId() = %d, track.offsets() = %d", collision.globalIndex(), track.globalIndex(), track.trackType(), track.matchMFTTrackId(), track.matchMCHTrackId(), track.offsets());

        // auto collision_in_track = track.collision_as<MyCollisions>();
        // auto bc_in_track = collision_in_track.bc_as<aod::BCsWithTimestamps>();
        // LOGF(info, "track.globalIndex() = %d , bc_in_track.globalBC() = %lld", track.globalIndex(), bc_in_track.globalBC());

        if (track.trackType() == static_cast<uint8_t>(o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack)) {
          if (fillQAHistogram) {
            fillTrackHistogram<3>(track, collision);
          }
          if (isSelected<false>(track, collision)) {
            map_new_sa_muon_index[std::make_pair(collision.globalIndex(), track.globalIndex())] = counter + offset;
            counter++;
          }
        } else if (track.trackType() == static_cast<uint8_t>(o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack)) {
          if (fillQAHistogram) {
            fillTrackHistogram<0>(track, collision);
          }
          if (map_new_sa_muon_index.find(std::make_pair(collision.globalIndex(), track.matchMCHTrackId())) == map_new_sa_muon_index.end()) {
            map_new_sa_muon_index[std::make_pair(collision.globalIndex(), track.matchMCHTrackId())] = counter + offset;
            counter++;
          }
          if (isSelected<false>(track, collision)) {
            map_new_sa_muon_index[std::make_pair(collision.globalIndex(), track.globalIndex())] = counter + offset;
            counter++;
          }
        }
      } // end of track loop

      for (const auto& [key, value] : map_new_sa_muon_index) {
        // int collisionId = std::get<0>(key);
        // int fwdtrackId = std::get<1>(key);
        // int new_fwdtrackId = value;
        auto track = tracks.iteratorAt(std::get<1>(key));
        if (track.trackType() == static_cast<uint8_t>(o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack)) {
          fillMuonTable<aod::MFTTracks>(track, collision, value);
        } else if (track.trackType() == static_cast<uint8_t>(o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack)) {
          fillMuonTable<aod::MFTTracks>(track, collision, map_new_sa_muon_index[std::make_pair(collision.globalIndex(), track.matchMCHTrackId())]);
        }
      }

      map_new_sa_muon_index.clear();
    } // end of collision loop
  }
  PROCESS_SWITCH(skimmerPrimaryMuon, processRec_TTCA_SWT, "process reconstructed info only with TTCA", false);

  Partition<MyTracksMC> global_muons_mc = o2::aod::fwdtrack::trackType == uint8_t(o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack); // MFT-MCH-MID
  Partition<MyTracksMC> sa_muons_mc = o2::aod::fwdtrack::trackType == uint8_t(o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack); // MCH-MID

  void processMC_SA(soa::Join<MyCollisions, aod::McCollisionLabels> const& collisions, aod::BCsWithTimestamps const&, MyTracksMC const& tracks, MFTTracksMC const&)
  {
    for (auto& collision : collisions) {
      if (!collision.has_mcCollision()) {
        continue;
      }
      auto bc = collision.template foundBC_as<aod::BCsWithTimestamps>();
      initCCDB(bc);
      fRegistry.fill(HIST("Event/hCollisionCounter"), 0.f);
      if (!collision.isSelected()) {
        continue;
      }

      auto sa_muons_mc_per_coll = sa_muons_mc->sliceByCached(o2::aod::fwdtrack::collisionId, collision.globalIndex(), cache);
      auto global_muons_mc_per_coll = global_muons_mc->sliceByCached(o2::aod::fwdtrack::collisionId, collision.globalIndex(), cache);

      int counter = 0;
      int offset = emprimarymuons.lastIndex() + 1;

      for (auto& track : sa_muons_mc_per_coll) {
        if (fillQAHistogram) {
          fillTrackHistogram<3>(track, collision);
        }
        if (isSelected<true>(track, collision)) {
          map_new_sa_muon_index[std::make_pair(collision.globalIndex(), track.globalIndex())] = counter + offset;
          counter++;
        }
      } // end of standalone muon loop
      for (auto& track : global_muons_mc_per_coll) {
        auto sa_muon = tracks.iteratorAt(track.matchMCHTrackId());
        if (!sa_muon.has_mcParticle()) {
          continue;
        }
        if (fillQAHistogram) {
          fillTrackHistogram<0>(track, collision);
        }
        if (map_new_sa_muon_index.find(std::make_pair(collision.globalIndex(), track.matchMCHTrackId())) == map_new_sa_muon_index.end()) { // don't apply muon selection to MCH-MID track in MFT-MCH-MID track
          map_new_sa_muon_index[std::make_pair(collision.globalIndex(), track.matchMCHTrackId())] = counter + offset;
          counter++;
        }
        if (isSelected<true>(track, collision)) {
          map_new_sa_muon_index[std::make_pair(collision.globalIndex(), track.globalIndex())] = counter + offset;
          counter++;
        }
      } // end of global muon loop

      // fill table after mapping
      for (const auto& [key, value] : map_new_sa_muon_index) {
        // int collisionId = std::get<0>(key);
        // int fwdtrackId = std::get<1>(key);
        // int new_fwdtrackId = value;
        // LOGF(info, "collisionId = %d, fwdtrackId = %d, new_fwdtrackId = %d", collisionId, fwdtrackId, new_fwdtrackId);
        auto track = tracks.iteratorAt(std::get<1>(key));
        if (track.trackType() == static_cast<uint8_t>(o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack)) {
          fillMuonTable<MFTTracksMC>(track, collision, value);
        } else if (track.trackType() == static_cast<uint8_t>(o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack)) {
          fillMuonTable<MFTTracksMC>(track, collision, map_new_sa_muon_index[std::make_pair(collision.globalIndex(), track.matchMCHTrackId())]);
        }
      }

      map_new_sa_muon_index.clear();
    } // end of collision loop
  }
  PROCESS_SWITCH(skimmerPrimaryMuon, processMC_SA, "process reconstructed and MC info", false);

  void processMC_TTCA(soa::Join<MyCollisions, aod::McCollisionLabels> const& collisions, aod::BCsWithTimestamps const&, MyTracksMC const& tracks, MFTTracksMC const&, aod::FwdTrackAssoc const& fwdtrackIndices)
  {
    for (auto& collision : collisions) {
      if (!collision.has_mcCollision()) {
        continue;
      }
      auto bc = collision.template foundBC_as<aod::BCsWithTimestamps>();
      initCCDB(bc);
      fRegistry.fill(HIST("Event/hCollisionCounter"), 0.f);
      if (!collision.isSelected()) {
        continue;
      }

      int counter = 0;
      int offset = emprimarymuons.lastIndex() + 1;

      auto fwdtrackIdsThisCollision = fwdtrackIndices.sliceBy(fwdtrackIndicesPerCollision, collision.globalIndex());
      for (auto& fwdtrackId : fwdtrackIdsThisCollision) {
        auto track = fwdtrackId.template fwdtrack_as<MyTracksMC>();
        // LOGF(info, "TTCA | track.globalIndex() = %d, track.trackType() = %d, track.matchMFTTrackId() = %d, track.matchMCHTrackId() = %d, track.offsets() = %d", track.globalIndex(), track.trackType(), track.matchMFTTrackId(), track.matchMCHTrackId(), track.offsets());

        if (track.trackType() == static_cast<uint8_t>(o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack)) {
          if (fillQAHistogram) {
            fillTrackHistogram<3>(track, collision);
          }
          if (isSelected<true>(track, collision)) {
            map_new_sa_muon_index[std::make_pair(collision.globalIndex(), track.globalIndex())] = counter + offset;
            counter++;
          }
        } else if (track.trackType() == static_cast<uint8_t>(o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack)) {
          auto sa_muon = tracks.iteratorAt(track.matchMCHTrackId());
          if (!sa_muon.has_mcParticle()) {
            continue;
          }
          if (fillQAHistogram) {
            fillTrackHistogram<0>(track, collision);
          }
          if (map_new_sa_muon_index.find(std::make_pair(collision.globalIndex(), track.matchMCHTrackId())) == map_new_sa_muon_index.end()) {
            map_new_sa_muon_index[std::make_pair(collision.globalIndex(), track.matchMCHTrackId())] = counter + offset;
            counter++;
          }
          if (isSelected<true>(track, collision)) {
            map_new_sa_muon_index[std::make_pair(collision.globalIndex(), track.globalIndex())] = counter + offset;
            counter++;
          }
        }
      } // end of track loop

      for (const auto& [key, value] : map_new_sa_muon_index) {
        // int collisionId = std::get<0>(key);
        // int fwdtrackId = std::get<1>(key);
        // int new_fwdtrackId = value;
        auto track = tracks.iteratorAt(std::get<1>(key));
        if (track.trackType() == static_cast<uint8_t>(o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack)) {
          fillMuonTable<MFTTracksMC>(track, collision, value);
        } else if (track.trackType() == static_cast<uint8_t>(o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack)) {
          fillMuonTable<MFTTracksMC>(track, collision, map_new_sa_muon_index[std::make_pair(collision.globalIndex(), track.matchMCHTrackId())]);
        }
      }

      map_new_sa_muon_index.clear();
    } // end of collision loop
  }
  PROCESS_SWITCH(skimmerPrimaryMuon, processMC_TTCA, "process reconstructed and MC info with TTCA", false);
};
struct associateAmbiguousMuon {
  Produces<aod::EMAmbiguousMuonSelfIds> em_amb_muon_ids;

  SliceCache cache;
  PresliceUnsorted<aod::EMPrimaryMuons> perTrack = o2::aod::emprimarymuon::fwdtrackId;
  std::vector<int> ambmuon_self_Ids;

  void process(aod::EMPrimaryMuons const& muons)
  {
    for (auto& muon : muons) {
      auto muons_with_same_trackId = muons.sliceBy(perTrack, muon.fwdtrackId());
      ambmuon_self_Ids.reserve(muons_with_same_trackId.size());
      for (auto& amp_muon : muons_with_same_trackId) {
        if (amp_muon.globalIndex() == muon.globalIndex()) { // don't store myself.
          continue;
        }
        ambmuon_self_Ids.emplace_back(amp_muon.globalIndex());
      }
      em_amb_muon_ids(ambmuon_self_Ids);
      ambmuon_self_Ids.clear();
      ambmuon_self_Ids.shrink_to_fit();
    }

    // for (auto& muon : muons) {
    //   auto sa_muon = muon.template matchMCHTrack_as<aod::EMPrimaryMuons>();
    //   LOGF(info, "muon.collisionId() = %d , muon.globalIndex() = %d, muon.fwdtrackId() = %d , muon.trackType() = %d, muon.matchMCHTrackId() = %d, sa_muon.fwdtrackId() = %d", muon.collisionId(), muon.globalIndex(), muon.fwdtrackId(), muon.trackType(), muon.matchMCHTrackId(), sa_muon.fwdtrackId());
    // }
  }
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<skimmerPrimaryMuon>(cfgc, TaskName{"skimmer-primary-muon"}),
                      adaptAnalysisTask<associateAmbiguousMuon>(cfgc, TaskName{"associate-ambiguous-muon"})};
}
