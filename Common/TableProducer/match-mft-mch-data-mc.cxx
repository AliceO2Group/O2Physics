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

#include "CCDB/BasicCCDBManager.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/CCDB/TriggerAliases.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/MftmchMatchingML.h"
#include "Common/DataModel/MatchMFTMuonData.h"
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

using namespace std;

using namespace o2;
using namespace o2::soa;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/DataTypes.h"
#include "Framework/runDataProcessing.h"

using MyCollisions = aod::Collisions;
using MyBCs = soa::Join<aod::BCs, aod::Timestamps, aod::MatchedToFT0>;
using MyMUONs = soa::Join<aod::FwdTracks, aod::FwdTracksCov, aod::McFwdTrackLabels>;
using MyMFTs = soa::Join<aod::MFTTracks, aod::McMFTTrackLabels>;

using MyCollision = MyCollisions::iterator;
using MyBC = MyBCs::iterator;
using MyMUON = MyMUONs::iterator;
using MyMFT = MyMFTs::iterator;

using SMatrix55 = ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepSym<double, 5>>;
using SMatrix5 = ROOT::Math::SVector<Double_t, 5>;

float mMu = TDatabasePDG::Instance()->GetParticle(13)->Mass();
int mRunNumber;
/*
  TLorentzVector muon1LV;
  TLorentzVector muon2LV;
  TLorentzVector dimuonLV;
*/
unordered_map<int, vector<int64_t>> map_mfttracks;
unordered_map<int, vector<int64_t>> map_muontracks;
unordered_map<int, bool> map_collisions;
unordered_map<int, bool> map_has_mfttracks_collisions;
unordered_map<int, bool> map_has_muontracks_collisions;
unordered_map<int, float> map_vtxz;
unordered_map<int, int> map_nmfttrack;

struct match_mft_mch_data_mc {

  ////  Variables for matching method
  Configurable<int> fMatchingMethod{"cfgMatchingMethod", 0, ""};

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
  Configurable<float> fPreselectMatchingX{"cfgPreselectMatchingX", 15.f, ""};
  Configurable<float> fPreselectMatchingY{"cfgPreselectMatchingY", 15.f, ""};

  ///    Variables to event mixing criteria
  Configurable<float> fSaveMixedMatchingParamsRate{"cfgSaveMixedMatchingParamsRate", 0.002f, ""};
  Configurable<int> fEventMaxDeltaNMFT{"cfgEventMaxDeltaNMFT", 1, ""};
  Configurable<float> fEventMaxDeltaVtxZ{"cfgEventMaxDeltaVtxZ", 1.f, ""};

  ////   Variables for selecting tag muon
  Configurable<float> fTagMassWindowMin{"cfgTagMassWindowMin", 2.8f, ""};
  Configurable<float> fTagMassWindowMax{"cfgTagMassWindowMax", 3.3f, ""};

  ////   Variables for ccdb
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};

  ////   Variables for Tag matching criteria
  Configurable<float> fSigmaXTagMuonCut{"cfgSigmaXTagMuonCut", 1.f, ""};
  Configurable<float> fMeanXTagMuonCut{"cfgMeanXTagMuonCut", 0.f, ""};
  Configurable<float> fSigmaYTagMuonCut{"cfgSigmaYTagMuonCut", 1.f, ""};
  Configurable<float> fMeanYTagMuonCut{"cfgMeanYTagMuonCut", 1.f, ""};

  Configurable<float> fSigmaEtaTagMuonCut{"cfgSigmaEtaTagMuonCut", 0.2f, ""};
  Configurable<float> fMeanEtaTagMuonCut{"cfgMeanEtaTagMuonCut", 0.f, ""};
  Configurable<float> fSigmaPhiTagMuonCut{"cfgSigmaPhiTagMuonCut", 0.2f, ""};
  Configurable<float> fMeanPhiTagMuonCut{"cfgMeanPhiTagMuonCut", 0.f, ""};

  template <typename MUON, typename Collision>
  class FindTagAndProbe
  {
   private:
    o2::dataformats::GlobalFwdTrack muontrack_at_pv[2];

    TLorentzVector mDimuon;
    MUON muontrack1;
    MUON muontrack2;
    Collision collision;
    int tagIdx, probeIdx;

    int16_t mQ;

    inline void fillCovarianceArray(MUON const& muontrack, float cov[15]) const
    {
      cov[0] = muontrack.cXX();
      cov[1] = muontrack.cXY();
      cov[2] = muontrack.cYY();
      cov[3] = muontrack.cPhiX();
      cov[4] = muontrack.cPhiY();
      cov[5] = muontrack.cPhiPhi();
      cov[6] = muontrack.cTglX();
      cov[7] = muontrack.cTglY();
      cov[8] = muontrack.cTglPhi();
      cov[9] = muontrack.cTglTgl();
      cov[10] = muontrack.c1PtX();
      cov[11] = muontrack.c1PtY();
      cov[12] = muontrack.c1PtPhi();
      cov[13] = muontrack.c1PtTgl();
      cov[14] = muontrack.c1Pt21Pt2();
    }

    inline o2::dataformats::GlobalFwdTrack propagateMUONtoPV(MUON const& muontrack) const
    {
      const double mz = muontrack.z();
      const double mchi2 = muontrack.chi2();
      const float mx = muontrack.x();
      const float my = muontrack.y();
      const float mphi = muontrack.phi();
      const float mtgl = muontrack.tgl();
      const float m1pt = muontrack.signed1Pt();

      float cov[15];
      fillCovarianceArray(muontrack, cov);
      SMatrix5 tpars(mx, my, mphi, mtgl, m1pt);
      SMatrix55 tcovs(cov, cov + 15);

      o2::track::TrackParCovFwd parcovmuontrack{mz, tpars, tcovs, mchi2};

      o2::dataformats::GlobalFwdTrack gtrack;
      gtrack.setParameters(tpars);
      gtrack.setZ(parcovmuontrack.getZ());
      gtrack.setCovariances(tcovs);

      o2::globaltracking::MatchGlobalFwd mMatching;
      auto mchtrack = mMatching.FwdtoMCH(gtrack);

      o2::mch::TrackExtrap::extrapToVertex(mchtrack, collision.posX(), collision.posY(), collision.posZ(), collision.covXX(), collision.covYY());

      auto fwdtrack = mMatching.MCHtoFwd(mchtrack);
      o2::dataformats::GlobalFwdTrack extrap_muontrack;
      extrap_muontrack.setParameters(fwdtrack.getParameters());
      extrap_muontrack.setZ(fwdtrack.getZ());
      extrap_muontrack.setCovariances(fwdtrack.getCovariances());

      return extrap_muontrack;
    }

    inline void setTagAndProbe()
    {
      if (muontrack1.pt() > muontrack2.pt()) {
        tagIdx = 0;
        probeIdx = 1;
      } else {
        tagIdx = 1;
        probeIdx = 0;
      }
    }

   public:
    inline FindTagAndProbe(const MUON& muon1, const MUON& muon2, const Collision& coll)
      : muontrack_at_pv(), mDimuon(), muontrack1(muon1), muontrack2(muon2), collision(coll), tagIdx(-1), probeIdx(-1), mQ(0)
    {
      mQ = muontrack1.sign() + muontrack2.sign();
      setTagAndProbe();
    }

    void calcMuonPairAtPV()
    {
      muontrack_at_pv[0] = propagateMUONtoPV(muontrack1);
      muontrack_at_pv[1] = propagateMUONtoPV(muontrack2);
      TLorentzVector vMuon1, vMuon2;
      vMuon1.SetPtEtaPhiM(muontrack_at_pv[0].getPt(), muontrack_at_pv[0].getEta(), muontrack_at_pv[0].getPhi(), mMu);
      vMuon2.SetPtEtaPhiM(muontrack_at_pv[1].getPt(), muontrack_at_pv[1].getEta(), muontrack_at_pv[1].getPhi(), mMu);
      mDimuon = vMuon1 + vMuon2;
    }
    inline int getTagMuonIndex() const { return tagIdx; }
    inline int getProbeMuonIndex() const { return probeIdx; }
    inline float getMass() const { return mDimuon.M(); }
    inline float getPt() const { return mDimuon.Pt(); }
    inline float getRap() const { return mDimuon.Rapidity(); }
    inline int16_t getCharge() const { return mQ; }
    inline const o2::dataformats::GlobalFwdTrack& getMuonAtPV(int idx) const { return muontrack_at_pv[idx]; }
  }; // end of class FindTagAndProbe

  template <typename MUON, typename MFT, typename Collision>
  class MatchingParamsML
  {
   private:
    MUON muontrack;
    MFT mfttrack;
    Collision collision;

    float mDX, mDY, mDPt, mDPhi, mDEta;
    float mGlobalMuonPtAtDCA, mGlobalMuonEtaAtDCA, mGlobalMuonPhiAtDCA, mGlobalMuonDCAx, mGlobalMuonDCAy, mGlobalMuonQ;
    int mMatchingType;

    o2::field::MagneticField* fieldB;
    o2::globaltracking::MatchGlobalFwd mMatching;

    inline o2::track::TrackParCovFwd propagateMFTtoMatchingPlane()
    {
      double covArr[15]{0.0};
      SMatrix55 tmftcovs(covArr, covArr + 15);

      SMatrix5 tmftpars(mfttrack.x(), mfttrack.y(), mfttrack.phi(), mfttrack.tgl(), mfttrack.signed1Pt());
      o2::track::TrackParCovFwd extrap_mfttrack{mfttrack.z(), tmftpars, tmftcovs, mfttrack.chi2()};

      double propVec[3] = {0.};
      float zPlane = 0.f;
      if (mMatchingType == MCH_FIRST_CLUSTER) {
        propVec[0] = muontrack.x() - mfttrack.x();
        propVec[1] = muontrack.y() - mfttrack.y();
        propVec[2] = muontrack.z() - mfttrack.z();
        zPlane = muontrack.z();
      } else if (mMatchingType == END_OF_ABSORBER || mMatchingType == BEGINING_OF_ABSORBER) {
        auto extrap_muontrack = propagateMUONtoMatchingPlane();
        propVec[0] = extrap_muontrack.getX() - mfttrack.x();
        propVec[1] = extrap_muontrack.getY() - mfttrack.y();
        propVec[2] = extrap_muontrack.getZ() - mfttrack.z();
        zPlane = (mMatchingType == END_OF_ABSORBER) ? -505.f : -90.f;
      } else {
        zPlane = mfttrack.z();
      }

      double centerZ[3] = {mfttrack.x() + propVec[0] / 2., mfttrack.y() + propVec[1] / 2., mfttrack.z() + propVec[2] / 2.};
      float Bz = fieldB->getBz(centerZ);
      extrap_mfttrack.propagateToZ(zPlane, Bz); // z in cm
      return extrap_mfttrack;
    }

    inline o2::dataformats::GlobalFwdTrack propagateMUONtoMatchingPlane()
    {
      float cov[15] = {
        muontrack.cXX(), muontrack.cXY(), muontrack.cYY(),
        muontrack.cPhiX(), muontrack.cPhiY(), muontrack.cPhiPhi(),
        muontrack.cTglX(), muontrack.cTglY(), muontrack.cTglPhi(),
        muontrack.cTglTgl(), muontrack.c1PtX(), muontrack.c1PtY(),
        muontrack.c1PtPhi(), muontrack.c1PtTgl(), muontrack.c1Pt21Pt2()};

      SMatrix5 tpars(muontrack.x(), muontrack.y(), muontrack.phi(), muontrack.tgl(), muontrack.signed1Pt());
      SMatrix55 tcovs(cov, cov + 15);
      double chi2 = muontrack.chi2();

      o2::track::TrackParCovFwd parcovmuontrack{muontrack.z(), tpars, tcovs, chi2};

      o2::dataformats::GlobalFwdTrack gtrack;
      gtrack.setParameters(tpars);
      gtrack.setZ(parcovmuontrack.getZ());
      gtrack.setCovariances(tcovs);

      auto mchtrack = mMatching.FwdtoMCH(gtrack);

      if (mMatchingType == MFT_LAST_CLUSTR) {
        o2::mch::TrackExtrap::extrapToVertexWithoutBranson(mchtrack, mfttrack.z());
      } else if (mMatchingType == END_OF_ABSORBER) {
        o2::mch::TrackExtrap::extrapToVertexWithoutBranson(mchtrack, -505.);
      } else if (mMatchingType == BEGINING_OF_ABSORBER) {
        o2::mch::TrackExtrap::extrapToVertexWithoutBranson(mchtrack, -90.);
      }

      auto fwdtrack = mMatching.MCHtoFwd(mchtrack);

      o2::dataformats::GlobalFwdTrack extrap_muontrack;
      extrap_muontrack.setParameters(fwdtrack.getParameters());
      extrap_muontrack.setZ(fwdtrack.getZ());
      extrap_muontrack.setCovariances(fwdtrack.getCovariances());
      return extrap_muontrack;
    }

    inline o2::track::TrackParCovFwd propagateMFTtoDCA()
    {
      double covArr[15]{0.0};
      SMatrix55 tmftcovs(covArr, covArr + 15);

      SMatrix5 tmftpars(mfttrack.x(), mfttrack.y(), mfttrack.phi(), mfttrack.tgl(), mfttrack.signed1Pt());
      o2::track::TrackParCovFwd extrap_mfttrack{mfttrack.z(), tmftpars, tmftcovs, mfttrack.chi2()};

      double propVec[3] = {};
      propVec[0] = collision.posX() - mfttrack.x();
      propVec[1] = collision.posY() - mfttrack.y();
      propVec[2] = collision.posZ() - mfttrack.z();

      double centerZ[3] = {mfttrack.x() + propVec[0] / 2., mfttrack.y() + propVec[1] / 2., mfttrack.z() + propVec[2] / 2.};
      float Bz = fieldB->getBz(centerZ);
      extrap_mfttrack.propagateToZ(collision.posZ(), Bz); // z in cm
      return extrap_mfttrack;
    }

    inline o2::dataformats::GlobalFwdTrack propagateMUONtoPV()
    {
      float cov[15] = {
        muontrack.cXX(), muontrack.cXY(), muontrack.cYY(),
        muontrack.cPhiX(), muontrack.cPhiY(), muontrack.cPhiPhi(),
        muontrack.cTglX(), muontrack.cTglY(), muontrack.cTglPhi(),
        muontrack.cTglTgl(), muontrack.c1PtX(), muontrack.c1PtY(),
        muontrack.c1PtPhi(), muontrack.c1PtTgl(), muontrack.c1Pt21Pt2()};

      SMatrix5 tpars(muontrack.x(), muontrack.y(), muontrack.phi(), muontrack.tgl(), muontrack.signed1Pt());
      SMatrix55 tcovs(cov, cov + 15);
      double chi2 = muontrack.chi2();

      o2::track::TrackParCovFwd parcovmuontrack{muontrack.z(), tpars, tcovs, chi2};

      o2::dataformats::GlobalFwdTrack gtrack;
      gtrack.setParameters(tpars);
      gtrack.setZ(parcovmuontrack.getZ());
      gtrack.setCovariances(tcovs);

      auto mchtrack = mMatching.FwdtoMCH(gtrack);
      o2::mch::TrackExtrap::extrapToVertex(mchtrack, collision.posX(), collision.posY(), collision.posZ(), collision.covXX(), collision.covYY());

      auto fwdtrack = mMatching.MCHtoFwd(mchtrack);
      o2::dataformats::GlobalFwdTrack extrap_muontrack;
      extrap_muontrack.setParameters(fwdtrack.getParameters());
      extrap_muontrack.setZ(fwdtrack.getZ());
      extrap_muontrack.setCovariances(fwdtrack.getCovariances());

      return extrap_muontrack;
    }

   public:
    enum MATCHING_TYPE { MCH_FIRST_CLUSTER,
                         MFT_LAST_CLUSTR,
                         END_OF_ABSORBER,
                         BEGINING_OF_ABSORBER };

    MatchingParamsML(MUON const& muon, MFT const& mft, Collision const& coll, int MType, o2::field::MagneticField* field) : muontrack(muon), mfttrack(mft), collision(coll), mDX(0.f), mDY(0.f), mDPt(0.f), mDPhi(0.f), mDEta(0.f), mGlobalMuonPtAtDCA(0.f), mGlobalMuonEtaAtDCA(0.f), mGlobalMuonPhiAtDCA(0.f), mGlobalMuonDCAx(0.f), mGlobalMuonDCAy(0.f), mGlobalMuonQ(0.f), mMatchingType(MType), fieldB(field) {}
    void calcMatchingParams()
    {
      auto mfttrack_on_matchingP = propagateMFTtoMatchingPlane();
      auto muontrack_on_matchingP = propagateMUONtoMatchingPlane();

      float dphiRaw = mfttrack_on_matchingP.getPhi() - muontrack_on_matchingP.getPhi();
      float dphi = TVector2::Phi_mpi_pi(dphiRaw);
      float deta = mfttrack_on_matchingP.getEta() - muontrack_on_matchingP.getEta();

      mDX = mfttrack_on_matchingP.getX() - muontrack_on_matchingP.getX();
      mDY = mfttrack_on_matchingP.getY() - muontrack_on_matchingP.getY();
      mDPt = mfttrack_on_matchingP.getPt() - muontrack_on_matchingP.getPt();
      mDPhi = dphi;
      mDEta = deta;
    }

    void calcGlobalMuonParams()
    {
      auto mfttrack_at_dca = propagateMFTtoDCA();
      auto muontrack_at_pv = propagateMUONtoPV();

      float momentum = muontrack_at_pv.getP();
      float theta = mfttrack_at_dca.getTheta();
      float phiTrack = mfttrack_at_dca.getPhi();
      float px = momentum * std::sin(theta) * std::cos(phiTrack);
      float py = momentum * std::sin(theta) * std::sin(phiTrack);

      mGlobalMuonQ = muontrack.sign() + mfttrack.sign();
      mGlobalMuonPtAtDCA = std::sqrt(px * px + py * py);
      mGlobalMuonEtaAtDCA = mfttrack_at_dca.getEta();
      mGlobalMuonPhiAtDCA = mfttrack_at_dca.getPhi();
      mGlobalMuonDCAx = mfttrack_at_dca.getX() - collision.posX();
      mGlobalMuonDCAy = mfttrack_at_dca.getY() - collision.posY();
    }

    inline float getDx() const { return mDX; }
    inline float getDy() const { return mDY; }
    inline float getDphi() const { return mDPhi; }
    inline float getDeta() const { return mDEta; }
    inline float getDpt() const { return mDPt; }
    inline float getGMPtAtDCA() const { return mGlobalMuonPtAtDCA; }
    inline float getGMEtaAtDCA() const { return mGlobalMuonEtaAtDCA; }
    inline float getGMPhiAtDCA() const { return mGlobalMuonPhiAtDCA; }
    inline float getGMDcaX() const { return mGlobalMuonDCAx; }
    inline float getGMDcaY() const { return mGlobalMuonDCAy; }
    inline float getGMDcaXY() const { return std::sqrt(mGlobalMuonDCAx * mGlobalMuonDCAx + mGlobalMuonDCAy * mGlobalMuonDCAy); }
    inline int16_t getGMQ() const { return static_cast<int16_t>(mGlobalMuonQ); }

  }; // end of class MatchingParamsML

  template <typename MUON, typename Collisions>
  o2::dataformats::GlobalFwdTrack propagateMUONtoPV(MUON const& muontrack, Collisions const& collisions)
  {
    auto collision = collisions.rawIteratorAt(muontrack.collisionId());
    o2::globaltracking::MatchGlobalFwd mMatching;
    o2::dataformats::GlobalFwdTrack extrap_muontrack;

    SMatrix5 tpars(muontrack.x(), muontrack.y(), muontrack.phi(), muontrack.tgl(), muontrack.signed1Pt());
    std::vector<float> v1{muontrack.cXX(), muontrack.cXY(), muontrack.cYY(),
                          muontrack.cPhiX(), muontrack.cPhiY(), muontrack.cPhiPhi(),
                          muontrack.cTglX(), muontrack.cTglY(), muontrack.cTglPhi(),
                          muontrack.cTglTgl(), muontrack.c1PtX(), muontrack.c1PtY(),
                          muontrack.c1PtPhi(), muontrack.c1PtTgl(), muontrack.c1Pt21Pt2()};
    SMatrix55 tcovs(v1.begin(), v1.end());
    double chi2 = muontrack.chi2();
    o2::track::TrackParCovFwd parcovmuontrack{muontrack.z(), tpars, tcovs, chi2};

    o2::dataformats::GlobalFwdTrack gtrack;
    gtrack.setParameters(tpars);
    gtrack.setZ(parcovmuontrack.getZ());
    gtrack.setCovariances(tcovs);
    auto mchtrack = mMatching.FwdtoMCH(gtrack);

    o2::mch::TrackExtrap::extrapToVertex(mchtrack, collision.posX(), collision.posY(), collision.posZ(), collision.covXX(), collision.covYY());

    auto fwdtrack = mMatching.MCHtoFwd(mchtrack);
    extrap_muontrack.setParameters(fwdtrack.getParameters());
    extrap_muontrack.setZ(fwdtrack.getZ());
    extrap_muontrack.setCovariances(fwdtrack.getCovariances());

    return extrap_muontrack;
  }

  inline bool isGoodTagDimuon(float M)
  {
    return !(M < fTagMassWindowMin || M > fTagMassWindowMax);
  }

  inline bool isGoodTagMatching(float mDX, float mDY, float mDEta, float mDPhi)
  {
    float dxNorm = (mDX - fMeanXTagMuonCut) / (fSigmaXTagMuonCut * 3);
    float dyNorm = (mDY - fMeanYTagMuonCut) / (fSigmaYTagMuonCut * 3);
    float detaNorm = (mDEta - fMeanEtaTagMuonCut) / (fSigmaEtaTagMuonCut * 3);
    float dphiNorm = (mDPhi - fMeanPhiTagMuonCut) / (fSigmaPhiTagMuonCut * 3);

    float rTagXY = dxNorm * dxNorm + dyNorm * dyNorm;
    float rTagEtaPhi = detaNorm * detaNorm + dphiNorm * dphiNorm;

    return (rTagXY < 1.f && rTagEtaPhi > 0.f);
  }

  template <typename MUON, typename MFT>
  bool isCorrectMatching(MUON const& muontrack, MFT const& mfttrack)
  {

    int idmuon = muontrack.mcParticleId();
    int idmft = mfttrack.mcParticleId();

    if (idmuon == -1 || idmft == -1)
      return false;
    if (idmuon != idmft)
      return false;
    else
      return true;
  };

  template <typename MUON>
  bool isGoodMuonQuality(MUON muontrack)
  {
    if (!muontrack.has_collision())
      return false;
    if (muontrack.trackType() != o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack)
      return false;
    if (muontrack.chi2() > fTrackChi2MchUp)
      return false;
    if (fRabsLow1 > muontrack.rAtAbsorberEnd() || muontrack.rAtAbsorberEnd() > fRabsUp2)
      return false;
    if (muontrack.rAtAbsorberEnd() < fRabsUp1 && fPdcaUp1 < muontrack.pDca())
      return false;
    if (muontrack.rAtAbsorberEnd() > fRabsLow2 && fPdcaUp2 < muontrack.pDca())
      return false;
    return true;
  }

  template <typename MUON>
  bool isGoodMuonKine(MUON muontrack)
  {
    if (fEtaMchLow > muontrack.getEta() || muontrack.getEta() > fEtaMchUp)
      return false;
    return true;
  }

  template <typename MFT>
  bool isGoodMFTQuality(MFT mfttrack)
  {
    if (!mfttrack.has_collision())
      return false;
    if (mfttrack.chi2() > fTrackChi2MftUp)
      return false;
    if (mfttrack.nClusters() < fTrackNClustMftLow)
      return false;
    return true;
  }

  template <typename MFT>
  bool isGoodMFTKine(MFT mfttrack)
  {
    if (fEtaMftLow > mfttrack.getEta() || mfttrack.getEta() > fEtaMftUp)
      return false;
    return true;
  }

  inline bool isPassMatchingPreselection(float Dx, float Dy)
  {
    return !(std::abs(Dx) > fPreselectMatchingX || std::abs(Dy) > fPreselectMatchingY);
  }

  template <typename MUONs, typename Collisions>
  void setMUONs(MUONs const& muontracks, Collisions const& collisions)
  {
    for (auto muontrack : muontracks) {
      if (!isGoodMuonQuality(muontrack))
        continue;
      o2::dataformats::GlobalFwdTrack muontrack_at_pv = propagateMUONtoPV(muontrack, collisions);
      if (!isGoodMuonKine(muontrack_at_pv))
        continue;

      auto collision = collisions.rawIteratorAt(muontrack.collisionId());

      bool& has = map_has_muontracks_collisions[muontrack.collisionId()];
      has = true;

      vector<int64_t>& arr_muontracks = map_muontracks[collision.globalIndex()];
      arr_muontracks.push_back(muontrack.globalIndex());
    }
  }

  template <typename MFT, typename Collisions>
  o2::track::TrackParCovFwd PropagateMFTtoDCA(MFT const& mfttrack, Collisions const& collisions, o2::field::MagneticField* field)
  {
    auto collision = collisions.rawIteratorAt(mfttrack.collisionId());
    std::vector<double> mftv1;
    SMatrix55 mftcovs{mftv1.begin(), mftv1.end()};
    SMatrix5 mftpars = {mfttrack.x(), mfttrack.y(), mfttrack.phi(), mfttrack.tgl(), mfttrack.signed1Pt()};
    o2::track::TrackParCovFwd mftpartrack = {mfttrack.z(), mftpars, mftcovs, mfttrack.chi2()};
    double propVec[3] = {fabs(mfttrack.x() - collision.posX()),
                         fabs(mfttrack.y() - collision.posY()),
                         fabs(mfttrack.z() - collision.posZ())};
    double centerZ[3] = {mfttrack.x() - propVec[0] / 2.,
                         mfttrack.y() - propVec[1] / 2.,
                         mfttrack.z() - propVec[2] / 2.};
    float Bz = field->getBz(centerZ);
    mftpartrack.propagateToZ(collision.posZ(), Bz);
    return mftpartrack;
  }

  template <typename MFTs, typename Collisions>
  void setMFTs(MFTs const& mfttracks, Collisions const& collisions, o2::field::MagneticField* field)
  {
    for (auto mfttrack : mfttracks) {
      if (!isGoodMFTQuality(mfttrack))
        continue;

      o2::track::TrackParCovFwd mfttrack_at_dca = PropagateMFTtoDCA(mfttrack, collisions, field);
      if (!isGoodMFTKine(mfttrack_at_dca))
        continue;

      auto collision = collisions.rawIteratorAt(mfttrack.collisionId());

      map_vtxz[mfttrack.collisionId()] = collision.posZ();
      map_nmfttrack[mfttrack.collisionId()] += 1;

      bool& has = map_has_mfttracks_collisions[mfttrack.collisionId()];
      has = true;

      vector<int64_t>& arr_mfttracks = map_mfttracks[collision.globalIndex()];
      arr_mfttracks.push_back(mfttrack.globalIndex());
    }
  }

  Produces<o2::aod::MatchParams> tableMatchingParams;
  Produces<o2::aod::TagMatchParams> tableTagMatchingParams;
  Produces<o2::aod::ProbeMatchParams> tableProbeMatchingParams;
  Produces<o2::aod::MixMatchParams> tableMixMatchingParams;
  Produces<o2::aod::MuonPair> tableMuonPair;

  Service<o2::ccdb::BasicCCDBManager> ccdbManager;

  o2::field::MagneticField* fieldB;
  o2::ccdb::CcdbApi ccdbApi;

  template <typename BC>
  void initCCDB(BC const& bc)
  {
    if (mRunNumber == bc.runNumber())
      return;

    mRunNumber = bc.runNumber();
    std::map<string, string> metadata;
    auto soreor = o2::ccdb::BasicCCDBManager::getRunDuration(ccdbApi, mRunNumber);
    auto ts = soreor.first;
    auto grpmag = ccdbApi.retrieveFromTFileAny<o2::parameters::GRPMagField>(grpmagPath, metadata, ts);
    o2::base::Propagator::initFieldFromGRP(grpmag);
    if (!o2::base::GeometryManager::isGeometryLoaded()) {
      ccdbManager->get<TGeoManager>(geoPath);
    }
    o2::mch::TrackExtrap::setField();
    fieldB = static_cast<o2::field::MagneticField*>(TGeoGlobalMagField::Instance()->GetField());
  }

  void init(o2::framework::InitContext&)
  {
    ccdbManager->setURL(ccdburl);
    ccdbManager->setCaching(true);
    ccdbManager->setLocalObjectValidityChecking();
    ccdbManager->setFatalWhenNull(false);
    ccdbApi.init(ccdburl);
    mRunNumber = 0;
  }

  void process(MyCollisions const& collisions,
               MyBCs const& bcs,
               MyMUONs const& muontracks,
               MyMFTs const& mfttracks)
  {
    LOG(info) << "Process()  ";
    map_muontracks.clear();
    map_mfttracks.clear();
    map_collisions.clear();
    map_has_muontracks_collisions.clear();
    map_has_mfttracks_collisions.clear();

    initCCDB(bcs.begin());
    setMUONs(muontracks, collisions);
    setMFTs(mfttracks, collisions, fieldB);

    for (auto map_has_muontracks_collision : map_has_muontracks_collisions) {
      auto idmuontrack_collisions = map_has_muontracks_collision.first;
      for (auto map_has_mfttracks_collision : map_has_mfttracks_collisions) {
        auto idmfttrack_collisions = map_has_mfttracks_collision.first;
        if (idmuontrack_collisions != idmfttrack_collisions)
          continue;
        map_collisions[idmfttrack_collisions] = true;
      }
    }

    for (auto const& map_collision : map_collisions) {
      auto const& collision = collisions.rawIteratorAt(map_collision.first);

      for (auto const& imuontrack1 : map_muontracks[map_collision.first]) {
        auto const& muontrack1 = muontracks.rawIteratorAt(imuontrack1);

        for (auto const& imfttrack1 : map_mfttracks[map_collision.first]) {
          auto const& mfttrack1 = mfttracks.rawIteratorAt(imfttrack1);

          MatchingParamsML<MyMUON, MyMFT, MyCollision> matching(muontrack1, mfttrack1, collision, fMatchingMethod, fieldB);
          matching.calcMatchingParams();

          if (!isPassMatchingPreselection(matching.getDx(), matching.getDy()))
            continue;

          matching.calcGlobalMuonParams();

          bool isTrue = isCorrectMatching(muontrack1, mfttrack1);

          tableMatchingParams(matching.getGMPtAtDCA(),
                              matching.getGMEtaAtDCA(),
                              static_cast<int16_t>(matching.getGMQ()),
                              matching.getDpt(),
                              matching.getDx(),
                              matching.getDy(),
                              matching.getDeta(),
                              matching.getDphi(),
                              isTrue);
        }

        for (auto const& map_mfttrack : map_mfttracks) {
          if (map_mfttrack.first == map_collision.first)
            continue;
          if (fabs(map_vtxz[map_mfttrack.first] - map_vtxz[map_collision.first]) > fEventMaxDeltaVtxZ)
            continue;
          if (fabs(map_nmfttrack[map_mfttrack.first] - map_nmfttrack[map_collision.first]) > fEventMaxDeltaNMFT)
            continue;

          for (auto const& imfttrack1 : map_mfttrack.second) {
            auto const& mfttrack1 = mfttracks.rawIteratorAt(imfttrack1);
            MatchingParamsML<MyMUON, MyMFT, MyCollision> matching(muontrack1, mfttrack1, collision, fMatchingMethod, fieldB);
            matching.calcMatchingParams();
            if (!isPassMatchingPreselection(matching.getDx(), matching.getDy()))
              continue;
            matching.calcGlobalMuonParams();

            bool isTrue = isCorrectMatching(muontrack1, mfttrack1);

            tableMixMatchingParams(matching.getGMPtAtDCA(),
                                   matching.getGMEtaAtDCA(),
                                   static_cast<int16_t>(matching.getGMQ()),
                                   matching.getDpt(),
                                   matching.getDx(),
                                   matching.getDy(),
                                   matching.getDeta(),
                                   matching.getDphi(),
                                   isTrue);
          }
        }

        for (auto const& imuontrack2 : map_muontracks[map_collision.first]) {

          if (imuontrack1 >= imuontrack2)
            continue;

          auto const& muontrack2 = muontracks.rawIteratorAt(imuontrack2);

          FindTagAndProbe<MyMUON, MyCollision> tagdimuon(muontrack1, muontrack2, collision);
          tagdimuon.calcMuonPairAtPV();
          tableMuonPair(tagdimuon.getCharge(), tagdimuon.getMass(), tagdimuon.getPt(), tagdimuon.getRap());

          if (!isGoodTagDimuon(tagdimuon.getMass()))
            continue;

          auto tagmuontrack = muontrack1;
          auto probemuontrack = muontrack2;

          if (tagdimuon.getTagMuonIndex() == 1) {
            tagmuontrack = muontrack2;
            probemuontrack = muontrack1;
          }

          int nTagMFTCand = 0;
          int nProbeMFTCand = 0;

          int IndexTagMFTCand = -1;
          float tagGMPtAtDCA = 0;
          // float tagGMEtaAtDCA = 0;

          float minimumR = 9999.;
          int minimumIndexProbeMFTCand = -1;

          unordered_map<int, vector<float>> map_tagMatchingParams;
          unordered_map<int, vector<float>> map_probeMatchingParams;

          for (auto const& imfttrack1 : map_mfttracks[map_collision.first]) {
            auto const& mfttrack1 = mfttracks.rawIteratorAt(imfttrack1);
            MatchingParamsML<MyMUON, MyMFT, MyCollision> matchingTag(tagmuontrack, mfttrack1, collision, fMatchingMethod, fieldB);
            matchingTag.calcMatchingParams();
            matchingTag.calcGlobalMuonParams();
            if (isGoodTagMatching(matchingTag.getDx(), matchingTag.getDy(), matchingTag.getDeta(), matchingTag.getDphi()) &&
                isPassMatchingPreselection(matchingTag.getDx(), matchingTag.getDy())) {
              bool isTrue = isCorrectMatching(tagmuontrack, mfttrack1);
              tableTagMatchingParams(matchingTag.getGMPtAtDCA(),
                                     matchingTag.getGMEtaAtDCA(),
                                     matchingTag.getGMQ(),
                                     matchingTag.getDpt(),
                                     matchingTag.getDx(),
                                     matchingTag.getDy(),
                                     matchingTag.getDeta(),
                                     matchingTag.getDphi(),
                                     isTrue);
              IndexTagMFTCand = mfttrack1.globalIndex();
              tagGMPtAtDCA = matchingTag.getGMPtAtDCA();
              // tagGMEtaAtDCA = matchingTag.getGMEtaAtDCA();
              ++nTagMFTCand;
            }
          } // end of loop imfttrack1

          if (nTagMFTCand != 1)
            continue;
          for (auto const& imfttrack1 : map_mfttracks[map_collision.first]) {
            auto const& mfttrack1 = mfttracks.rawIteratorAt(imfttrack1);
            if (mfttrack1.globalIndex() == IndexTagMFTCand)
              continue;
            MatchingParamsML<MyMUON, MyMFT, MyCollision> matchingProbe(probemuontrack, mfttrack1, collision, fMatchingMethod, fieldB);
            matchingProbe.calcMatchingParams();
            if (isPassMatchingPreselection(matchingProbe.getDx(), matchingProbe.getDy())) {
              float R = sqrt(matchingProbe.getDx() * matchingProbe.getDx() + matchingProbe.getDy() * matchingProbe.getDy());
              bool isTrue = isCorrectMatching(probemuontrack, mfttrack1);
              matchingProbe.calcGlobalMuonParams();
              vector<float>& probeMatchingParams = map_probeMatchingParams[nProbeMFTCand];
              probeMatchingParams.push_back(tagGMPtAtDCA);
              probeMatchingParams.push_back(matchingProbe.getGMPtAtDCA());
              probeMatchingParams.push_back(matchingProbe.getGMEtaAtDCA());
              probeMatchingParams.push_back(matchingProbe.getGMQ());
              probeMatchingParams.push_back(matchingProbe.getDpt());
              probeMatchingParams.push_back(matchingProbe.getDx());
              probeMatchingParams.push_back(matchingProbe.getDy());
              probeMatchingParams.push_back(matchingProbe.getDeta());
              probeMatchingParams.push_back(matchingProbe.getDphi());
              probeMatchingParams.push_back(isTrue);
              if (R < minimumR) {
                minimumIndexProbeMFTCand = nProbeMFTCand;
                minimumR = R;
              }
              ++nProbeMFTCand;
            }
          } // end of loop imfttrack1

          if (nProbeMFTCand < 1)
            continue;

          if (minimumIndexProbeMFTCand > -1) {
            vector<float>& probeMatchingParams = map_probeMatchingParams[minimumIndexProbeMFTCand];
            tableProbeMatchingParams(probeMatchingParams[0],
                                     probeMatchingParams[1],
                                     probeMatchingParams[2],
                                     static_cast<int16_t>(probeMatchingParams[3]),
                                     probeMatchingParams[4],
                                     probeMatchingParams[5],
                                     probeMatchingParams[6],
                                     probeMatchingParams[7],
                                     probeMatchingParams[8],
                                     static_cast<bool>(probeMatchingParams[9]));
          }

        } // end of loop imuontrack2
      } // end of loop imuontrack1
    } // end of loop map_collision

  } // end of processMC
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<match_mft_mch_data_mc>(cfgc)};
}
