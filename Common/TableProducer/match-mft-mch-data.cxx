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
#include <map>
#include <string>
#include <unordered_set>
#include <unordered_map>
#include <vector>

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/DataTypes.h"
#include "Framework/runDataProcessing.h"
#include "CCDB/BasicCCDBManager.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/CCDB/TriggerAliases.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/MftmchMatchingML.h"
#include "Common/Core/trackUtilities.h"
#include "PWGDQ/DataModel/ReducedInfoTables.h"
#include "PWGDQ/Core/VarManager.h"
#include "PWGDQ/Core/HistogramManager.h"
#include "PWGDQ/Core/AnalysisCut.h"
#include "PWGDQ/Core/AnalysisCompositeCut.h"
#include "PWGDQ/Core/HistogramsLibrary.h"
#include "PWGDQ/Core/CutsLibrary.h"
#include "DataFormatsGlobalTracking/RecoContainerCreateTracksVariadic.h"
#include "DetectorsVertexing/VertexTrackMatcher.h"
#include "ReconstructionDataFormats/PrimaryVertex.h"
#include "ReconstructionDataFormats/VtxTrackIndex.h"
#include "ReconstructionDataFormats/VtxTrackRef.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "CommonDataFormat/InteractionRecord.h"
#include "DetectorsVertexing/PVertexerParams.h"
#include "MathUtils/Primitive2D.h"
#include "DataFormatsGlobalTracking/RecoContainer.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/DataModel/MatchMFTFT0.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "Field/MagneticField.h"
#include "TGeoGlobalMagField.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "EventFiltering/Zorro.h"
#include "ReconstructionDataFormats/TrackFwd.h"
#include "Math/MatrixFunctions.h"
#include "Math/SMatrix.h"
#include "MFTTracking/Tracker.h"
#include "MCHTracking/TrackParam.h"
#include "MCHTracking/TrackExtrap.h"
#include "GlobalTracking/MatchGlobalFwd.h"
#include <TLorentzVector.h>
#include <TVector3.h>
#include "TDatabasePDG.h"

using namespace o2;
using namespace o2::soa;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;

using SMatrix55 = ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepSym<double, 5>>;
using SMatrix5 = ROOT::Math::SVector<Double_t, 5>;

// using MyEvents = soa::Join<aod::Collisions, aod::EvSels>;
// using MyEventsWithMults = soa::Join<aod::Collisions, aod::EvSels, aod::Mults>;
// using MyEventsWithFilter = soa::Join<aod::Collisions, aod::EvSels, aod::DQEventFilter>;
// using MyEventsWithMultsAndFilter = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::MultsExtra, aod::DQEventFilter>;
// using MyEventsWithCent = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs>;
// using MyEventsWithCentAndMults = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs, aod::Mults>;
using MyMuons = soa::Join<aod::FwdTracks, aod::FwdTracksCov, aod::FwdTracksDCA>;
using MyMFTs = aod::MFTTracks;
// using MyMuonsWithCov = soa::Join<aod::FwdTracks, aod::FwdTracksCov, aod::FwdTracksDCA>;
// using MyMuonsColl = soa::Join<aod::FwdTracks, aod::FwdTracksDCA, aod::FwdTrkCompColls>;
// using MyMuonsCollWithCov = soa::Join<aod::FwdTracks, aod::FwdTracksCov, aod::FwdTracksDCA, aod::FwdTrkCompColls>;
using MyBCs = soa::Join<aod::BCs, aod::Timestamps, aod::Run3MatchedToBCSparse>;
using ExtBCs = soa::Join<aod::BCs, aod::Timestamps, aod::MatchedBCCollisionsSparseMulti, aod::MatchedToFT0>;

float mMu = TDatabasePDG::Instance()->GetParticle(13)->Mass();

TRandom* rnd = new TRandom();

TLorentzVector muon1LV;
TLorentzVector muon2LV;
TLorentzVector dimuonLV;

TVector3 V1;
TVector3 V2;

namespace o2::aod
{

namespace muon_params
{
DECLARE_SOA_COLUMN(TRACKCHI2, trackChi2, float);
DECLARE_SOA_COLUMN(RABS, rabs, float);
DECLARE_SOA_COLUMN(Q, q, int16_t);
DECLARE_SOA_COLUMN(PT, pt, float);
DECLARE_SOA_COLUMN(ETA, eta, float);
DECLARE_SOA_COLUMN(PHI, phi, float);
DECLARE_SOA_COLUMN(HASMFT, has_mft, bool);
} // namespace muon_params

DECLARE_SOA_TABLE(MUONParams, "AOD", "MUON",
                  muon_params::TRACKCHI2,
                  muon_params::RABS,
                  muon_params::Q,
                  muon_params::PT,
                  muon_params::ETA,
                  muon_params::PHI,
                  muon_params::HASMFT);

namespace mft_params
{
DECLARE_SOA_COLUMN(NCLUST, nclust, int);
DECLARE_SOA_COLUMN(ISCA, isCA, bool);
DECLARE_SOA_COLUMN(TRACKCHI2, trackChi2, float);
DECLARE_SOA_COLUMN(Q, q, int16_t);
DECLARE_SOA_COLUMN(PT_AT_DCA, pt_dca, float);
DECLARE_SOA_COLUMN(ETA_AT_DCA, eta_dca, float);
DECLARE_SOA_COLUMN(PHI_AT_DCA, phi_dca, float);
DECLARE_SOA_COLUMN(DCAx, dcax, float);
DECLARE_SOA_COLUMN(DCAy, dcay, float);
} // namespace mft_params

DECLARE_SOA_TABLE(MFTParams, "AOD", "MFT",
                  mft_params::NCLUST,
                  mft_params::ISCA,
                  mft_params::TRACKCHI2,
                  mft_params::Q,
                  mft_params::PT_AT_DCA,
                  mft_params::ETA_AT_DCA,
                  mft_params::PHI_AT_DCA,
                  mft_params::DCAx,
                  mft_params::DCAy);

namespace matching_params
{
// matching parameters
DECLARE_SOA_COLUMN(NClustMFTTracks, nClustMFT, int);
DECLARE_SOA_COLUMN(Chi2MFTTracks, chi2MFT, float);

DECLARE_SOA_COLUMN(DeltaP, dp_mchplane, float);
DECLARE_SOA_COLUMN(DeltaPt, dpt_mchplane, float);
DECLARE_SOA_COLUMN(DeltaEta, deta_mchplane, float);
DECLARE_SOA_COLUMN(DeltaPhi, dphi_mchplane, float);
DECLARE_SOA_COLUMN(DeltaX, dx_mchplane, float);
DECLARE_SOA_COLUMN(DeltaY, dy_mchplane, float);

DECLARE_SOA_COLUMN(MchPt, mchpt, float);
DECLARE_SOA_COLUMN(MchEta, mcheta, float);
DECLARE_SOA_COLUMN(MchPhi, mchphi, float);
DECLARE_SOA_COLUMN(MchQ, mchq, float);

DECLARE_SOA_COLUMN(MftPt, mftpt, float);
DECLARE_SOA_COLUMN(MftEta, mfteta, float);
DECLARE_SOA_COLUMN(MftPhi, mftphi, float);
DECLARE_SOA_COLUMN(MftQ, mftq, float);

DECLARE_SOA_COLUMN(MftDCA, mftdca, float);

} // namespace matching_params

DECLARE_SOA_TABLE(MatchParams, "AOD", "MATCHING",
                  matching_params::DeltaPt,
                  matching_params::DeltaEta,
                  matching_params::DeltaPhi,
                  matching_params::DeltaX,
                  matching_params::DeltaY,
                  matching_params::MchPt,
                  matching_params::MchEta,
                  matching_params::MchQ,
                  matching_params::MftEta,
                  matching_params::MftQ);

namespace mix_matching_params
{
// matching parameters
DECLARE_SOA_COLUMN(NClustMFTTracks, nClustMFT, int);
DECLARE_SOA_COLUMN(Chi2MFTTracks, chi2MFT, float);

DECLARE_SOA_COLUMN(DeltaP, dp_mchplane, float);
DECLARE_SOA_COLUMN(DeltaPt, dpt_mchplane, float);
DECLARE_SOA_COLUMN(DeltaEta, deta_mchplane, float);
DECLARE_SOA_COLUMN(DeltaPhi, dphi_mchplane, float);
DECLARE_SOA_COLUMN(DeltaX, dx_mchplane, float);
DECLARE_SOA_COLUMN(DeltaY, dy_mchplane, float);

DECLARE_SOA_COLUMN(MchPt, mchpt, float);
DECLARE_SOA_COLUMN(MchEta, mcheta, float);
DECLARE_SOA_COLUMN(MchPhi, mchphi, float);
DECLARE_SOA_COLUMN(MchQ, mchq, float);

DECLARE_SOA_COLUMN(MftPt, mftpt, float);
DECLARE_SOA_COLUMN(MftEta, mfteta, float);
DECLARE_SOA_COLUMN(MftPhi, mftphi, float);
DECLARE_SOA_COLUMN(MftQ, mftq, float);
DECLARE_SOA_COLUMN(MftDCA, mftdca, float);
} // namespace mix_matching_params

DECLARE_SOA_TABLE(MixMatchParams, "AOD", "MIXMATCHING",
                  mix_matching_params::DeltaPt,
                  mix_matching_params::DeltaEta,
                  mix_matching_params::DeltaPhi,
                  mix_matching_params::DeltaX,
                  mix_matching_params::DeltaY,
                  mix_matching_params::MchPt,
                  mix_matching_params::MchEta,
                  mix_matching_params::MchQ,
                  mix_matching_params::MftEta,
                  mix_matching_params::MftQ);

namespace tag_matching_params
{
// matching parameters
DECLARE_SOA_COLUMN(DeltaPt, dpt_mchplane, float);
DECLARE_SOA_COLUMN(DeltaEta, deta_mchplane, float);
DECLARE_SOA_COLUMN(DeltaPhi, dphi_mchplane, float);
DECLARE_SOA_COLUMN(DeltaX, dx_mchplane, float);
DECLARE_SOA_COLUMN(DeltaY, dy_mchplane, float);

DECLARE_SOA_COLUMN(MchPt, mchpt, float);
DECLARE_SOA_COLUMN(MchEta, mcheta, float);
DECLARE_SOA_COLUMN(MchQ, mchq, float);

DECLARE_SOA_COLUMN(MftEta, mfteta, float);
DECLARE_SOA_COLUMN(MftQ, mftq, float);
DECLARE_SOA_COLUMN(IsTaged, isTaged, bool);
} // namespace tag_matching_params

DECLARE_SOA_TABLE(TagMatchParams, "AOD", "TAGMATCHING",
                  tag_matching_params::DeltaPt,
                  tag_matching_params::DeltaEta,
                  tag_matching_params::DeltaPhi,
                  tag_matching_params::DeltaX,
                  tag_matching_params::DeltaY,
                  tag_matching_params::MchPt,
                  tag_matching_params::MchEta,
                  tag_matching_params::MchQ,
                  tag_matching_params::MftEta,
                  tag_matching_params::MftQ,
                  tag_matching_params::IsTaged);

namespace probe_matching_params
{
// matching parameters
DECLARE_SOA_COLUMN(NMFTCandTagMuon, nTagMFT, int);
DECLARE_SOA_COLUMN(NMFTCandProbeMuon, nProbeMFT, int);
DECLARE_SOA_COLUMN(TagMuonP, tagmuonp, float);

DECLARE_SOA_COLUMN(NClustMFTTracks, nClustMFT, int);
DECLARE_SOA_COLUMN(Chi2MFTTracks, chi2MFT, float);

DECLARE_SOA_COLUMN(DeltaP, dp_mchplane, float);
DECLARE_SOA_COLUMN(DeltaPt, dpt_mchplane, float);
DECLARE_SOA_COLUMN(DeltaEta, deta_mchplane, float);
DECLARE_SOA_COLUMN(DeltaPhi, dphi_mchplane, float);
DECLARE_SOA_COLUMN(DeltaX, dx_mchplane, float);
DECLARE_SOA_COLUMN(DeltaY, dy_mchplane, float);

DECLARE_SOA_COLUMN(MchPt, mchpt, float);
DECLARE_SOA_COLUMN(MchEta, mcheta, float);
DECLARE_SOA_COLUMN(MchPhi, mchphi, float);
DECLARE_SOA_COLUMN(MchQ, mchq, float);

DECLARE_SOA_COLUMN(MftPt, mftpt, float);
DECLARE_SOA_COLUMN(MftEta, mfteta, float);
DECLARE_SOA_COLUMN(MftPhi, mftphi, float);
DECLARE_SOA_COLUMN(MftQ, mftq, float);
DECLARE_SOA_COLUMN(MftDCA, mftdca, float);
} // namespace probe_matching_params

DECLARE_SOA_TABLE(ProbeMatchParams, "AOD", "PROBEMATCHING",
                  probe_matching_params::NMFTCandTagMuon,
                  probe_matching_params::NMFTCandProbeMuon,
                  probe_matching_params::TagMuonP,
                  probe_matching_params::DeltaPt,
                  probe_matching_params::DeltaEta,
                  probe_matching_params::DeltaPhi,
                  probe_matching_params::DeltaX,
                  probe_matching_params::DeltaY,
                  probe_matching_params::MchPt,
                  probe_matching_params::MchEta,
                  probe_matching_params::MchQ,
                  probe_matching_params::MftEta,
                  probe_matching_params::MftQ);

namespace muon_pair
{
DECLARE_SOA_COLUMN(NMFT, nMft, int);
DECLARE_SOA_COLUMN(Q, q, int16_t);
DECLARE_SOA_COLUMN(M, m, float);
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Rap, rap, float);
} // namespace muon_pair

DECLARE_SOA_TABLE(MuonPair, "AOD", "DIMUON", muon_pair::NMFT, muon_pair::Q, muon_pair::M, muon_pair::Pt, muon_pair::Rap);

namespace tag_muon_pair
{
DECLARE_SOA_COLUMN(NMFT, nMft, int);
DECLARE_SOA_COLUMN(Q, q, int16_t);
DECLARE_SOA_COLUMN(M, m, float);
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Rap, rap, float);
} // namespace tag_muon_pair

DECLARE_SOA_TABLE(TagMuonPair, "AOD", "TAGDIMUON", tag_muon_pair::NMFT, tag_muon_pair::Q, tag_muon_pair::M, tag_muon_pair::Pt, tag_muon_pair::Rap);
} // namespace o2::aod

struct match_mft_mch_data {

  Produces<o2::aod::MatchParams> matchingParams;
  Produces<o2::aod::TagMatchParams> tagmatchingParams;
  Produces<o2::aod::ProbeMatchParams> probematchingParams;
  Produces<o2::aod::MixMatchParams> mixmatchingParams;
  Produces<o2::aod::MuonPair> muonPairs;
  Produces<o2::aod::TagMuonPair> tagmuonPairs;
  Produces<o2::aod::MUONParams> muonParams;
  Produces<o2::aod::MFTParams> mftParams;

  HistogramRegistry registry{
    "registry",
    {{"hMchP", "MCH track total momentum (at the first station); p [GeV/c]; Counts", {HistType::kTH1F, {{2000, 0, 200}}}},
     {"hMchCorrP", "MCH track total momentum (propagated to PV); p [GeV/c]; Counts", {HistType::kTH1F, {{2000, 0, 200}}}},
     {"hMassCorrMchPair", "Corrected MCH track pair mass (propagated to PV); m [GeV/c^{2}]; Counts", {HistType::kTH1F, {{1000, 0, 10}}}}}};

  ////   Variables for selecting muon tracks
  Configurable<float> fEtaMchLow{"cfgEtaMchLow", -4.0f, ""};
  Configurable<float> fEtaMchUp{"cfgEtaMchUp", -2.5f, ""};
  Configurable<float> fRabsLow1{"cfgRabsLow1", 17.6f, ""};
  Configurable<float> fRabsUp1{"cfgRabsUp1", 26.5f, ""};
  Configurable<float> fRabsLow2{"cfgRabsLow2", 26.5f, ""};
  Configurable<float> fRabsUp2{"cfgRabsUp2", 89.5f, ""};
  Configurable<float> fPdcaUp1{"cfgPdcaUp1", 594.f, ""};
  Configurable<float> fPdcaUp2{"cfgPdcaUp2", 324.f, ""};
  Configurable<float> fTrackChi2MchUp{"cfgTrackChi2MchUp", 5.f, ""};
  Configurable<float> fMatchingChi2MchMidUp{"cfgMatchingChi2MchMidUp", 999.f, ""};

  ////   Variables for selecting mft tracks
  Configurable<float> fEtaMftLow{"cfgEtaMftlow", -3.6f, ""};
  Configurable<float> fEtaMftUp{"cfgEtaMftup", -2.5f, ""};
  Configurable<int> fTrackNClustMftLow{"cfgTrackNClustMftLow", 7, ""};
  Configurable<float> fTrackChi2MftUp{"cfgTrackChi2MftUp", 999.f, ""};

  ///    Variables to add preselection for the matching table
  Configurable<float> fPreselectMatchingX{"cfgPreselectMatchingX", 999.f, ""};
  Configurable<float> fPreselectMatchingY{"cfgPreselectMatchingY", 999.f, ""};

  ///    Variables to event mixing criteria
  Configurable<float> fSaveMixedMatchingParamsRate{"cfgSaveMixedMatchingParamsRate", 0.002f, ""};
  Configurable<int> fEventMaxDeltaNMFT{"cfgEventMaxDeltaNMFT", 1, ""};
  Configurable<float> fEventMaxDeltaVtxZ{"cfgEventMaxDeltaVtxZ", 1.f, ""};

  ////   Variables for selecting tag muon
  Configurable<float> fTagMassWindowMin{"cfgTagMassWindowMin", 2.8f, ""};
  Configurable<float> fTagMassWindowMax{"cfgTagMassWindowMax", 3.3f, ""};
  Configurable<float> fTagXWindowLow{"cfgTagXWindowLow", -0.80f, ""};
  Configurable<float> fTagXWindowUp{"cfgTagXWindowUp", 0.84f, ""};
  Configurable<float> fTagYWindowLow{"cfgTagYWindowLow", -0.64f, ""};
  Configurable<float> fTagYWindowUp{"cfgTagYWindowUp", 0.95f, ""};
  Configurable<float> fTagPhiWindowLow{"cfgTagPhiWindowLow", -0.041f, ""};
  Configurable<float> fTagPhiWindowUp{"cfgTagPhiWindowUp", 0.039f, ""};
  Configurable<float> fTagEtaWindowLow{"cfgTagEtaWindowLow", -0.045f, ""};
  Configurable<float> fTagEtaWindowUp{"cfgTagEtaWindowUp", 0.033f, ""};

  ////   Variables for selecting probe muon
  Configurable<float> fProbeXWindowLow{"cfgProbeXWindowLow", -10.f, ""};
  Configurable<float> fProbeXWindowUp{"cfgProbeXWindowUp", 10.f, ""};
  Configurable<float> fProbeYWindowLow{"cfgProbeYWindowLow", -10.f, ""};
  Configurable<float> fProbeYWindowUp{"cfgProbeYWindowUp", 10.f, ""};
  Configurable<float> fProbePhiWindowLow{"cfgProbePhiWindowLow", -1.0f, ""};
  Configurable<float> fProbePhiWindowUp{"cfgProbePhiWindowUp", 1.0f, ""};
  Configurable<float> fProbeEtaWindowLow{"cfgProbeEtaWindowLow", -1.0f, ""};
  Configurable<float> fProbeEtaWindowUp{"cfgProbeEtaWindowUp", 1.0f, ""};

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};

  // o2::parameters::GRPMagField* grpmag = nullptr;
  o2::globaltracking::MatchGlobalFwd mMatching;
  o2::field::MagneticField* fieldB;

  o2::ccdb::CcdbApi ccdbApi;
  int mRunNumber;

  void init(o2::framework::InitContext&)
  {
    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);
    ccdbApi.init(ccdburl);
    mRunNumber = 0;
  }

  void initCCDB(ExtBCs::iterator const& bc)
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
    fieldB = static_cast<o2::field::MagneticField*>(TGeoGlobalMagField::Instance()->GetField());
  }

  enum ProagationPoint { ToVtx,
                         ToDCA };

  template <typename T>
  o2::dataformats::GlobalFwdTrack PropagateMUONTrack(T const& muon, int PropType)
  {

    auto collision = muon.collision();

    o2::dataformats::GlobalFwdTrack propmuon;

    double chi2 = muon.chi2();

    SMatrix5 tpars(muon.x(), muon.y(), muon.phi(), muon.tgl(), muon.signed1Pt());
    std::vector<float> v1{muon.cXX(), muon.cXY(), muon.cYY(), muon.cPhiX(), muon.cPhiY(),
                          muon.cPhiPhi(), muon.cTglX(), muon.cTglY(), muon.cTglPhi(), muon.cTglTgl(),
                          muon.c1PtX(), muon.c1PtY(), muon.c1PtPhi(), muon.c1PtTgl(), muon.c1Pt21Pt2()};
    if (isGoodMUONTrack(muon)) {
      SMatrix55 tcovs(v1.begin(), v1.end());
      o2::track::TrackParCovFwd fwdtrack{muon.z(), tpars, tcovs, chi2};

      o2::dataformats::GlobalFwdTrack track;
      track.setParameters(tpars);
      track.setZ(fwdtrack.getZ());
      track.setCovariances(tcovs);
      auto mchTrack = mMatching.FwdtoMCH(track);
      if (PropType == ProagationPoint::ToVtx)
        o2::mch::TrackExtrap::extrapToVertex(mchTrack, collision.posX(), collision.posY(), collision.posZ(), collision.covXX(), collision.covYY());
      else if (PropType == ProagationPoint::ToDCA)
        o2::mch::TrackExtrap::extrapToVertexWithoutBranson(mchTrack, collision.posZ());

      auto proptrack = mMatching.MCHtoFwd(mchTrack);
      propmuon.setParameters(proptrack.getParameters());
      propmuon.setZ(proptrack.getZ());
      propmuon.setCovariances(proptrack.getCovariances());
    }

    v1.clear();
    v1.shrink_to_fit();

    return propmuon;
  }

  template <typename T>
  o2::track::TrackParCovFwd PropagateMFT(T const& mfttrack, int PropType)
  {
    std::vector<double> mftv1;
    SMatrix55 mftcovs{mftv1.begin(), mftv1.end()};
    SMatrix5 mftpars = {mfttrack.x(), mfttrack.y(), mfttrack.phi(), mfttrack.tgl(), mfttrack.signed1Pt()};
    o2::track::TrackParCovFwd mftpartrack = {mfttrack.z(), mftpars, mftcovs, mfttrack.chi2()};
    if (PropType == ProagationPoint::ToDCA) {
      auto collision = mfttrack.collision();
      double propVec[3] = {fabs(mfttrack.x() - collision.posX()), fabs(mfttrack.y() - collision.posY()), fabs(mfttrack.z() - collision.posZ())};
      double centerZ[3] = {mfttrack.x() - propVec[0] / 2., mfttrack.y() - propVec[1] / 2., mfttrack.z() - propVec[2] / 2.};
      float Bz = fieldB->getBz(centerZ);
      mftpartrack.propagateToZ(collision.posZ(), Bz);
    }
    return mftpartrack;
  }

  template <typename MFT, typename FWD>
  o2::track::TrackParCovFwd PropagateMFTtoMatchingPlane(MFT const& mfttrack, FWD const& fwdtrack)
  {
    std::vector<double> v1; // Temporary null vector for the computation of the covariance matrix
    double propVec[3] = {fwdtrack.x() - mfttrack.x(), fwdtrack.y() - mfttrack.y(), fwdtrack.z() - mfttrack.z()};
    double centerZ[3] = {mfttrack.x() + propVec[0] / 2., mfttrack.y() + propVec[1] / 2., mfttrack.z() + propVec[2] / 2.};
    float Bz = fieldB->getBz(centerZ); // gives error if the propagator is not initFielded
    SMatrix55 tmftcovs(v1.begin(), v1.end());
    SMatrix5 tmftpars(mfttrack.x(), mfttrack.y(), mfttrack.phi(), mfttrack.tgl(), mfttrack.signed1Pt());
    o2::track::TrackParCovFwd extrap_mfttrack{mfttrack.z(), tmftpars, tmftcovs, mfttrack.chi2()};
    extrap_mfttrack.propagateToZ(fwdtrack.z(), Bz); // z in cm
    return extrap_mfttrack;
  }

  template <typename T>
  bool isGoodMUONTrack(T track)
  {
    if (!track.has_collision())
      return false;
    if (track.trackType() != o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack)
      return false;
    if (track.chi2() > fTrackChi2MchUp)
      return false;
    if (fRabsLow1 > track.rAtAbsorberEnd() || track.rAtAbsorberEnd() > fRabsUp2)
      return false;
    if (track.rAtAbsorberEnd() < fRabsUp1 && fPdcaUp1 < track.pDca())
      return false;
    if (track.rAtAbsorberEnd() > fRabsLow2 && fPdcaUp2 < track.pDca())
      return false;
    return true;
  }

  template <typename T>
  int selectTagMUON(T track1, T track2)
  {
    if (track1.pt() > track2.pt()) {
      return track1.globalIndex();
    } else {
      return track2.globalIndex();
    }
  }

  template <typename T>
  int selectProbeMUONTrack(T track1, T track2)
  {
    if (track1.pt() < track2.pt()) {
      return track1.globalIndex();
    } else {
      return track2.globalIndex();
    }
  }

  bool isGoodKenematicMUONTrack(o2::dataformats::GlobalFwdTrack track)
  {
    if (fEtaMchLow > track.getEta() || track.getEta() > fEtaMchUp)
      return false;
    return true;
  }

  template <typename T>
  bool isGoodMFTTrack(T track)
  {
    if (!track.has_collision())
      return false;
    if (track.chi2() > fTrackChi2MftUp)
      return false;
    if (track.nClusters() < fTrackNClustMftLow)
      return false;
    return true;
  }

  bool isGoodKenematicMFTTrack(o2::track::TrackParCovFwd track)
  {
    if (fEtaMftLow > track.getEta() || track.getEta() > fEtaMftUp) {
      return false;
    }
    return true;
  }

  void process(aod::Collisions const&, ExtBCs const& ebcs,
               MyMuons const& fwdtracks, MyMFTs const& mfttracks)
  {
    initCCDB(ebcs.begin());

    std::unordered_set<int> bcs_mfttrack;
    std::unordered_map<int, float> map_vtxZ;
    std::unordered_map<int, int> nmfttracks;
    std::unordered_map<int, std::vector<int64_t>> map_mfttraks;

    for (const auto& mfttrack : mfttracks) {

      if (!isGoodMFTTrack(mfttrack)) {
        continue;
      }

      bcs_mfttrack.insert(mfttrack.collisionId());
      std::vector<int64_t>& tracks = map_mfttraks[mfttrack.collisionId()];
      tracks.push_back(mfttrack.globalIndex());

      o2::track::TrackParCovFwd mftpartrack = PropagateMFT(mfttrack, ProagationPoint::ToDCA);

      if (!isGoodKenematicMFTTrack(mftpartrack)) {
        continue;
      }

      auto collision = mfttrack.collision();

      map_vtxZ[mfttrack.collisionId()] = collision.posZ();

      float dx = mftpartrack.getX() - collision.posX();
      float dy = mftpartrack.getY() - collision.posY();

      mftParams(mfttrack.nClusters(), mfttrack.isCA(), mfttrack.chi2(), mfttrack.sign(), mftpartrack.getPt(), mftpartrack.getEta(), mftpartrack.getPhi(), dx, dy);
      nmfttracks[mfttrack.collisionId()]++;
    }

    std::unordered_map<int, int> nfwdtracks;
    std::unordered_map<int, std::vector<int64_t>> map_fwdtraks;

    for (auto fwdtrack : fwdtracks) {

      if (!isGoodMUONTrack(fwdtrack)) {
        continue;
      }

      o2::dataformats::GlobalFwdTrack propmuonAtPV = PropagateMUONTrack(fwdtrack, ProagationPoint::ToVtx);
      if (!isGoodKenematicMUONTrack(propmuonAtPV)) {
        continue;
      }

      std::vector<int64_t>& tracks = map_fwdtraks[fwdtrack.collisionId()];
      tracks.push_back(fwdtrack.globalIndex());

      bool hasMFT = false;

      std::vector<int64_t>& mfttracks = map_mfttraks[fwdtrack.collisionId()];

      if (mfttracks.size() > 0) {
        hasMFT = true;
      }

      muonParams(fwdtrack.chi2(), fwdtrack.rAtAbsorberEnd(), fwdtrack.sign(), propmuonAtPV.getPt(), propmuonAtPV.getEta(), propmuonAtPV.getPhi(), hasMFT);
      nfwdtracks[fwdtrack.collisionId()]++;
    }

    for (auto fwdtrack1 : fwdtracks) {

      if (!isGoodMUONTrack(fwdtrack1)) {
        continue;
      }

      int ibc = fwdtrack1.collisionId();
      auto collision = fwdtrack1.collision();

      o2::dataformats::GlobalFwdTrack fwdtrackAtPV1 = PropagateMUONTrack(fwdtrack1, ProagationPoint::ToVtx);
      if (!isGoodKenematicMUONTrack(fwdtrackAtPV1)) {
        continue;
      }

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      ///////////////                                    MIXED EVENT                                     ///////////////
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      for (auto bc_mfttrack : bcs_mfttrack) {

        if (ibc == bc_mfttrack)
          continue;
        if (fabs(nmfttracks[ibc] - nmfttracks[bc_mfttrack]) > fEventMaxDeltaNMFT) {
          continue;
        }
        if (fabs(map_vtxZ[bc_mfttrack] - collision.posZ()) > fEventMaxDeltaVtxZ) {
          continue;
        }

        std::vector<int64_t>& mfttrackGlobalIndex = map_mfttraks[bc_mfttrack];

        for (int idmfttrack1 = 0; idmfttrack1 < static_cast<int>(mfttrackGlobalIndex.size()); ++idmfttrack1) {

          auto mfttrack1 = mfttracks.rawIteratorAt(mfttrackGlobalIndex[idmfttrack1]);
          o2::track::TrackParCovFwd mfttrack_at_matching = PropagateMFTtoMatchingPlane(mfttrack1, fwdtrack1);

          V1.SetPtEtaPhi(mfttrack_at_matching.getPt(), mfttrack_at_matching.getEta(), mfttrack_at_matching.getPhi());
          V2.SetPtEtaPhi(fwdtrack1.pt(), fwdtrack1.eta(), fwdtrack1.phi());

          double deltaPt = mfttrack_at_matching.getPt() - fwdtrack1.pt();
          double deltaX = mfttrack_at_matching.getX() - fwdtrack1.x();
          double deltaY = mfttrack_at_matching.getY() - fwdtrack1.y();
          double deltaPhi = V1.DeltaPhi(V2);
          double deltaEta = mfttrack_at_matching.getEta() - fwdtrack1.eta();

          if (fabs(deltaX) > fPreselectMatchingX) {
            continue;
          }
          if (fabs(deltaY) > fPreselectMatchingY) {
            continue;
          }

          o2::track::TrackParCovFwd mfttrack_at_dca = PropagateMFT(mfttrack1, ProagationPoint::ToDCA);

          mixmatchingParams(deltaPt, deltaEta, deltaPhi, deltaX, deltaY, fwdtrackAtPV1.getPt(), fwdtrackAtPV1.getEta(), fwdtrack1.sign(), mfttrack_at_dca.getEta(), mfttrack1.sign());
        }
      } // end of loop bc_mfttrack

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      ///////////////                                 SAME EVENT                                         ///////////////
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      std::vector<int64_t>& mfttrackGlobalIndex = map_mfttraks[ibc];

      for (int idmfttrack1 = 0; idmfttrack1 < static_cast<int>(mfttrackGlobalIndex.size()); ++idmfttrack1) {

        auto mfttrack1 = mfttracks.rawIteratorAt(mfttrackGlobalIndex[idmfttrack1]);
        o2::track::TrackParCovFwd mfttrack_at_matching = PropagateMFTtoMatchingPlane(mfttrack1, fwdtrack1);
        V1.SetPtEtaPhi(mfttrack_at_matching.getPt(), mfttrack_at_matching.getEta(), mfttrack_at_matching.getPhi());
        V2.SetPtEtaPhi(fwdtrack1.pt(), fwdtrack1.eta(), fwdtrack1.phi());

        double deltaPt = mfttrack_at_matching.getPt() - fwdtrack1.pt();
        double deltaX = mfttrack_at_matching.getX() - fwdtrack1.x();
        double deltaY = mfttrack_at_matching.getY() - fwdtrack1.y();
        double deltaPhi = V1.DeltaPhi(V2);
        double deltaEta = mfttrack_at_matching.getEta() - fwdtrack1.eta();

        if (fabs(deltaX) > fPreselectMatchingX) {
          continue;
        }
        if (fabs(deltaY) > fPreselectMatchingY) {
          continue;
        }

        o2::track::TrackParCovFwd mfttrack_at_dca = PropagateMFT(mfttrack1, ProagationPoint::ToDCA);

        matchingParams(deltaPt, deltaEta, deltaPhi, deltaX, deltaY,
                       fwdtrackAtPV1.getPt(), fwdtrackAtPV1.getEta(), fwdtrack1.sign(),
                       mfttrack_at_dca.getEta(), mfttrack1.sign());

      } // end of loop idmfttrack1

      std::vector<int64_t>& fwdtrackGlobalIndex = map_fwdtraks[ibc];

      for (int idfwdtrack2 = 0; idfwdtrack2 < static_cast<int>(fwdtrackGlobalIndex.size()); ++idfwdtrack2) {

        if (fwdtrack1.globalIndex() == fwdtrackGlobalIndex[idfwdtrack2]) {
          continue;
        }

        auto fwdtrack2 = fwdtracks.rawIteratorAt(fwdtrackGlobalIndex[idfwdtrack2]);

        if (!isGoodMUONTrack(fwdtrack2)) {
          continue;
        }

        o2::dataformats::GlobalFwdTrack fwdtrackAtPV2 = PropagateMUONTrack(fwdtrack2, ProagationPoint::ToVtx);
        if (!isGoodKenematicMUONTrack(fwdtrackAtPV2)) {
          continue;
        }

        muon1LV.SetPtEtaPhiM(fwdtrackAtPV1.getPt(), fwdtrackAtPV1.getEta(), fwdtrackAtPV1.getPhi(), mMu);
        muon2LV.SetPtEtaPhiM(fwdtrackAtPV2.getPt(), fwdtrackAtPV2.getEta(), fwdtrackAtPV2.getPhi(), mMu);
        dimuonLV = muon1LV + muon2LV;

        muonPairs(nmfttracks[ibc], fwdtrack1.sign() + fwdtrack2.sign(), dimuonLV.M(), dimuonLV.Pt(), dimuonLV.Rapidity());

        if (fabs(fwdtrack1.sign() + fwdtrack2.sign()) > 0) {
          continue;
        }
        if (fTagMassWindowMin > dimuonLV.M() || dimuonLV.M() > fTagMassWindowMax) {
          continue;
        }
        if (nmfttracks[ibc] < 1) {
          continue;
        }

        tagmuonPairs(nmfttracks[ibc], fwdtrack1.sign() + fwdtrack2.sign(), dimuonLV.M(), dimuonLV.Pt(), dimuonLV.Rapidity());

        bool isGoodTag = false;
        int nMFTCandsTagMUON = 0;

        auto tagfwdtrack = fwdtracks.rawIteratorAt(selectTagMUON(fwdtrack1, fwdtrack2));
        o2::dataformats::GlobalFwdTrack tagfwdtrackAtPV = PropagateMUONTrack(tagfwdtrack, ProagationPoint::ToVtx);

        for (int idmfttrack1 = 0; idmfttrack1 < static_cast<int>(mfttrackGlobalIndex.size()); ++idmfttrack1) {

          auto mfttrack1 = mfttracks.rawIteratorAt(mfttrackGlobalIndex[idmfttrack1]);
          o2::track::TrackParCovFwd mfttrack_at_matching = PropagateMFTtoMatchingPlane(mfttrack1, tagfwdtrack);

          V1.SetPtEtaPhi(mfttrack_at_matching.getPt(), mfttrack_at_matching.getEta(), mfttrack_at_matching.getPhi());
          V2.SetPtEtaPhi(tagfwdtrack.pt(), tagfwdtrack.eta(), tagfwdtrack.phi());

          double deltaPt = mfttrack_at_matching.getPt() - tagfwdtrack.pt();
          double deltaX = mfttrack_at_matching.getX() - tagfwdtrack.x();
          double deltaY = mfttrack_at_matching.getY() - tagfwdtrack.y();
          double deltaPhi = V1.DeltaPhi(V2);
          double deltaEta = mfttrack_at_matching.getEta() - tagfwdtrack.eta();

          if (fabs(deltaX) > fPreselectMatchingX) {
            continue;
          }
          if (fabs(deltaY) > fPreselectMatchingY) {
            continue;
          }

          o2::track::TrackParCovFwd mfttrack_at_dca = PropagateMFT(mfttrack1, ProagationPoint::ToDCA);

          bool dummyTag = false;

          if (fTagXWindowLow < deltaX && deltaX < fTagXWindowUp &&
              fTagXWindowLow < deltaY && deltaY < fTagXWindowUp &&
              fTagPhiWindowLow < deltaPhi && deltaPhi < fTagPhiWindowUp &&
              fTagEtaWindowLow < deltaEta && deltaEta < fTagEtaWindowUp) {
            isGoodTag = true;
            dummyTag = true;
            nMFTCandsTagMUON++;
          }

          tagmatchingParams(deltaPt, deltaEta, deltaPhi, deltaX, deltaY,
                            tagfwdtrackAtPV.getPt(), tagfwdtrackAtPV.getEta(), tagfwdtrack.sign(), mfttrack_at_dca.getEta(), mfttrack1.sign(), dummyTag);
        }

        if (!isGoodTag) {
          continue;
        }

        auto probefwdtrack = fwdtracks.rawIteratorAt(selectProbeMUONTrack(fwdtrack1, fwdtrack2));
        o2::dataformats::GlobalFwdTrack probefwdtrackAtPV = PropagateMUONTrack(probefwdtrack, ProagationPoint::ToVtx);

        int nMFTCandsProbeMUON = 0;
        /// Counting the number of MFT tracks in a probe muon track
        for (int idmfttrack1 = 0; idmfttrack1 < static_cast<int>(mfttrackGlobalIndex.size()); ++idmfttrack1) {

          auto mfttrack1 = mfttracks.rawIteratorAt(mfttrackGlobalIndex[idmfttrack1]);
          o2::track::TrackParCovFwd mfttrack_at_matching = PropagateMFTtoMatchingPlane(mfttrack1, probefwdtrack);
          V1.SetPtEtaPhi(mfttrack_at_matching.getPt(), mfttrack_at_matching.getEta(), mfttrack_at_matching.getPhi());
          V2.SetPtEtaPhi(probefwdtrack.pt(), probefwdtrack.eta(), probefwdtrack.phi());

          double deltaX = mfttrack_at_matching.getX() - probefwdtrack.x();
          double deltaY = mfttrack_at_matching.getY() - probefwdtrack.y();
          double deltaPhi = V1.DeltaPhi(V2);
          double deltaEta = mfttrack_at_matching.getEta() - probefwdtrack.eta();

          if (fProbeXWindowLow < deltaX && deltaX < fProbeXWindowUp &&
              fProbeXWindowLow < deltaY && deltaY < fProbeXWindowUp &&
              fProbePhiWindowLow < deltaPhi && deltaPhi < fProbePhiWindowUp &&
              fProbeEtaWindowLow < deltaEta && deltaEta < fProbeEtaWindowUp) {
            nMFTCandsProbeMUON++;
          }
        }

        for (int idmfttrack1 = 0; idmfttrack1 < static_cast<int>(mfttrackGlobalIndex.size()); ++idmfttrack1) {

          auto mfttrack1 = mfttracks.rawIteratorAt(mfttrackGlobalIndex[idmfttrack1]);
          o2::track::TrackParCovFwd mfttrack_at_matching = PropagateMFTtoMatchingPlane(mfttrack1, probefwdtrack);

          V1.SetPtEtaPhi(mfttrack_at_matching.getPt(), mfttrack_at_matching.getEta(), mfttrack_at_matching.getPhi());
          V2.SetPtEtaPhi(probefwdtrack.pt(), probefwdtrack.eta(), probefwdtrack.phi());

          double deltaPt = mfttrack_at_matching.getPt() - probefwdtrack.pt();
          double deltaX = mfttrack_at_matching.getX() - probefwdtrack.x();
          double deltaY = mfttrack_at_matching.getY() - probefwdtrack.y();
          double deltaPhi = V1.DeltaPhi(V2);
          double deltaEta = mfttrack_at_matching.getEta() - probefwdtrack.eta();

          if (fabs(deltaX) > fPreselectMatchingX) {
            continue;
          }
          if (fabs(deltaY) > fPreselectMatchingY) {
            continue;
          }

          o2::track::TrackParCovFwd mfttrack_at_dca = PropagateMFT(mfttrack1, ProagationPoint::ToDCA);

          probematchingParams(nMFTCandsTagMUON, nMFTCandsProbeMUON, tagfwdtrack.p(), deltaPt, deltaEta, deltaPhi, deltaX, deltaY, probefwdtrackAtPV.getPt(), probefwdtrackAtPV.getEta(), probefwdtrack.sign(), mfttrack_at_dca.getEta(), mfttrack1.sign());
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<match_mft_mch_data>(cfgc)};
}
