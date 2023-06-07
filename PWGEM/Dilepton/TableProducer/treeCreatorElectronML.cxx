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

#include "Math/Vector4D.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/trackUtilities.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/PIDResponse.h"
#include "CommonConstants/PhysicsConstants.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "CCDB/BasicCCDBManager.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::constants::physics;

using FullTracksExt = soa::Join<aod::Tracks, aod::TracksCov, aod::TracksExtra, aod::TracksDCA,
                                aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr,
                                aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFbeta>;

using FullTrackExt = FullTracksExt::iterator;

using FullTracksExtMC = soa::Join<FullTracksExt, aod::McTrackLabels>;
using FullTrackExtMC = FullTracksExtMC::iterator;

namespace o2::aod
{

namespace mycollision // reconstructed collision information
{
DECLARE_SOA_COLUMN(Bz, bz, float);         //! in unit of kG.
DECLARE_SOA_COLUMN(MCPosX, mcposX, float); //!
DECLARE_SOA_COLUMN(MCPosY, mcposY, float); //!
DECLARE_SOA_COLUMN(MCPosZ, mcposZ, float); //!
} // namespace mycollision
DECLARE_SOA_TABLE(MyCollisions, "AOD", "MYCOLLISION", //! vertex information of collision
                  o2::soa::Index<>, bc::GlobalBC, bc::RunNumber, collision::PosX, collision::PosY, collision::PosZ, collision::NumContrib, evsel::Sel8, mycollision::Bz,
                  mccollision::GeneratorsID, mycollision::MCPosX, mycollision::MCPosY, mycollision::MCPosZ,
                  mult::MultTPC, mult::MultFV0A, mult::MultFV0C, mult::MultFT0A, mult::MultFT0C,
                  mult::MultFDDA, mult::MultFDDC, mult::MultZNA, mult::MultZNC, mult::MultTracklets, mult::MultNTracksPV, mult::MultNTracksPVeta1);
using MyCollision = MyCollisions::iterator;

namespace mytrack
{
DECLARE_SOA_INDEX_COLUMN(MyCollision, mycollision);              //!
DECLARE_SOA_COLUMN(Sign, sign, int);                             //!
DECLARE_SOA_COLUMN(TPCNClsFound, tpcNClsFound, int);             //!
DECLARE_SOA_COLUMN(TPCNClsCrossedRows, tpcNClsCrossedRows, int); //!
DECLARE_SOA_COLUMN(DCAresXY, dcaresXY, float);                   //!
DECLARE_SOA_COLUMN(DCAresZ, dcaresZ, float);                     //!
DECLARE_SOA_COLUMN(MCPt, mcpt, float);                           //!
DECLARE_SOA_COLUMN(MCEta, mceta, float);                         //!
DECLARE_SOA_COLUMN(MCPhi, mcphi, float);                         //!
DECLARE_SOA_COLUMN(MCVx, mcvx, float);                           //!
DECLARE_SOA_COLUMN(MCVy, mcvy, float);                           //!
DECLARE_SOA_COLUMN(MCVz, mcvz, float);                           //!
DECLARE_SOA_COLUMN(IsPhysicalPrimary, isPhysicalPrimary, bool);  //!
DECLARE_SOA_COLUMN(MotherPdgCode, motherpdgCode, int);           //!
DECLARE_SOA_COLUMN(GrandMotherPdgCode, grandmotherpdgCode, int); //!
} // namespace mytrack

// reconstructed track information
DECLARE_SOA_TABLE(MyTracks, "AOD", "MYTRACK", //!
                  o2::soa::Index<>, mytrack::MyCollisionId, mytrack::Sign,
                  track::Pt, track::Eta, track::Phi,
                  track::DcaXY, track::DcaZ, mytrack::DCAresXY, mytrack::DCAresZ,
                  track::TPCNClsFindable, track::TPCNClsFindableMinusFound, track::TPCNClsFindableMinusCrossedRows,
                  mytrack::TPCNClsFound, mytrack::TPCNClsCrossedRows,
                  track::TPCChi2NCl, track::TPCInnerParam,
                  track::TPCSignal, pidtpc::TPCNSigmaEl, pidtpc::TPCNSigmaMu, pidtpc::TPCNSigmaPi, pidtpc::TPCNSigmaKa, pidtpc::TPCNSigmaPr,
                  pidtofbeta::Beta, pidtof::TOFNSigmaEl, pidtof::TOFNSigmaMu, pidtof::TOFNSigmaPi, pidtof::TOFNSigmaKa, pidtof::TOFNSigmaPr,
                  track::ITSClusterMap, track::ITSChi2NCl, track::DetectorMap,
                  mytrack::MCPt, mytrack::MCEta, mytrack::MCPhi,
                  mytrack::MCVx, mytrack::MCVy, mytrack::MCVz,
                  mcparticle::PdgCode, mytrack::IsPhysicalPrimary, mytrack::MotherPdgCode, mytrack::GrandMotherPdgCode,
                  // dynamic column
                  track::TPCCrossedRowsOverFindableCls<track::TPCNClsFindable, track::TPCNClsFindableMinusCrossedRows>,
                  track::TPCFoundOverFindableCls<track::TPCNClsFindable, track::TPCNClsFindableMinusFound>,
                  track::ITSNCls<track::ITSClusterMap>,
                  track::HasITS<track::DetectorMap>, track::HasTPC<track::DetectorMap>,
                  track::HasTRD<track::DetectorMap>, track::HasTOF<track::DetectorMap>);

// iterators
using MyTrack = MyTracks::iterator;

namespace mypair
{
DECLARE_SOA_INDEX_COLUMN(MyCollision, mycollision);                       //!
DECLARE_SOA_INDEX_COLUMN_FULL(PosTrack, posTrack, int, MyTracks, "_Pos"); //!
DECLARE_SOA_INDEX_COLUMN_FULL(NegTrack, negTrack, int, MyTracks, "_Neg"); //!

DECLARE_SOA_COLUMN(M, m, float);                 //!
DECLARE_SOA_COLUMN(Pt, pt, float);               //!
DECLARE_SOA_COLUMN(Eta, eta, float);             //!
DECLARE_SOA_COLUMN(Phi, phi, float);             //!
DECLARE_SOA_COLUMN(PhiV, phiv, float);           //!
DECLARE_SOA_COLUMN(PairDCAxy, pairDCAxy, float); //!
DECLARE_SOA_COLUMN(PairDCAz, pairDCAz, float);   //!

DECLARE_SOA_COLUMN(IsSM, isSM, bool);        //!
DECLARE_SOA_COLUMN(IsHF, isHF, int);         //!
DECLARE_SOA_COLUMN(PairType, pairtype, int); //!
} // namespace mypair

// reconstructed track information
DECLARE_SOA_TABLE(MyPairs, "AOD", "MYPAIR", //!
                  o2::soa::Index<>, mypair::MyCollisionId, mypair::PosTrackId, mypair::NegTrackId,
                  mypair::M, mypair::Pt, mypair::Eta, mypair::Phi, mypair::PhiV, mypair::PairDCAxy, mypair::PairDCAz,
                  mypair::IsSM, mypair::IsHF, mypair::PairType,
                  mcparticle::PdgCode, mcparticle::StatusCode, mcparticle::Flags,
                  mcparticle::Vx, mcparticle::Vy, mcparticle::Vz,

                  // dynamic column
                  mcparticle::ProducedByGenerator<mcparticle::Flags>,
                  mcparticle::FromBackgroundEvent<mcparticle::Flags>,
                  mcparticle::GetGenStatusCode<mcparticle::Flags, mcparticle::StatusCode>,
                  mcparticle::GetProcess<mcparticle::Flags, mcparticle::StatusCode>,
                  mcparticle::IsPhysicalPrimary<mcparticle::Flags>);

// iterators
using MyPair = MyPairs::iterator;

} // namespace o2::aod

struct TreeCreatorElectronML {
  enum EM_HFeeType {
    kUndef = -1,
    kCe_Ce = 0,        // ULS
    kBe_Be = 1,        // ULS
    kBCe_BCe = 2,      // ULS
    kBCe_Be_SameB = 3, // ULS
    kBCe_Be_DiffB = 4, // LS
  };
  enum EM_EEPairType {
    kULS = 0,
    kLSpp = 1,
    kLSnn = 2,
  };

  SliceCache cache;
  Produces<o2::aod::MyCollisions> mycollision;
  Produces<o2::aod::MyPairs> mypair;
  Produces<o2::aod::MyTracks> mytrack;

  // Basic checks
  HistogramRegistry registry{
    "registry",
    {
      {"hEventCounter", "hEventCounter", {HistType::kTH1F, {{5, 0.5f, 5.5f}}}},
    },
  };

  // Configurables
  Configurable<int> mincrossedrows{"mincrossedrows", 70, "min. crossed rows"};
  Configurable<float> maxchi2tpc{"maxchi2tpc", 4.0, "max. chi2/NclsTPC"};
  Configurable<float> maxeta{"maxeta", 0.9, "eta acceptance"};

  int mRunNumber;
  float d_bz;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  void init(InitContext& context)
  {
    mRunNumber = 0;
    d_bz = 0;
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);
  }

  void initCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }

    auto run3grp_timestamp = bc.timestamp();
    o2::parameters::GRPObject* grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>("GLO/GRP/GRP", run3grp_timestamp);
    o2::parameters::GRPMagField* grpmag = 0x0;
    if (grpo) {
      o2::base::Propagator::initFieldFromGRP(grpo);
      // Fetch magnetic field from ccdb for current collision
      d_bz = grpo->getNominalL3Field();
      LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << d_bz << " kZG";
    } else {
      grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>("GLO/Config/GRPMagField", run3grp_timestamp);
      if (!grpmag) {
        LOG(fatal) << "Got nullptr from CCDB for path "
                   << "GLO/Config/GRPMagField"
                   << " of object GRPMagField and "
                   << "GLO/GRP/GRP"
                   << " of object GRPObject for timestamp " << run3grp_timestamp;
      }
      o2::base::Propagator::initFieldFromGRP(grpmag);
      // Fetch magnetic field from ccdb for current collision
      d_bz = std::lround(5.f * grpmag->getL3Current() / 30000.f);
      LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << d_bz << " kZG";
    }
    mRunNumber = bc.runNumber();
  }

  template <typename TTrack>
  float get_phiv(TTrack const& t1, TTrack const& t2)
  {
    // cos(phiv) = w*a /|w||a|
    // with w = u x v
    // and  a = u x z / |u x z|   , unit vector perpendicular to v12 and z-direction (magnetic field)
    // u = v12 / |v12|            , the unit vector of v12
    // v = v1 x v2 / |v1 x v2|    , unit vector perpendicular to v1 and v2

    // float bz = fgFitterTwoProngBarrel.getBz();

    ROOT::Math::PtEtaPhiMVector v1(t1.pt(), t1.eta(), t1.phi(), o2::constants::physics::MassElectron);
    ROOT::Math::PtEtaPhiMVector v2(t2.pt(), t2.eta(), t2.phi(), o2::constants::physics::MassElectron);
    ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;

    bool swapTracks = false;
    if (v1.Pt() < v2.Pt()) { // ordering of track, pt1 > pt2
      ROOT::Math::PtEtaPhiMVector v3 = v1;
      v1 = v2;
      v2 = v3;
      swapTracks = true;
    }

    // momentum of e+ and e- in (ax,ay,az) axis. Note that az=0 by definition.
    // vector product of pep X pem
    float vpx = 0, vpy = 0, vpz = 0;
    if (t1.sign() * t2.sign() > 0) { // Like Sign
      if (!swapTracks) {
        if (d_bz * t1.sign() < 0) {
          vpx = v1.Py() * v2.Pz() - v1.Pz() * v2.Py();
          vpy = v1.Pz() * v2.Px() - v1.Px() * v2.Pz();
          vpz = v1.Px() * v2.Py() - v1.Py() * v2.Px();
        } else {
          vpx = v2.Py() * v1.Pz() - v2.Pz() * v1.Py();
          vpy = v2.Pz() * v1.Px() - v2.Px() * v1.Pz();
          vpz = v2.Px() * v1.Py() - v2.Py() * v1.Px();
        }
      } else { // swaped tracks
        if (d_bz * t2.sign() < 0) {
          vpx = v1.Py() * v2.Pz() - v1.Pz() * v2.Py();
          vpy = v1.Pz() * v2.Px() - v1.Px() * v2.Pz();
          vpz = v1.Px() * v2.Py() - v1.Py() * v2.Px();
        } else {
          vpx = v2.Py() * v1.Pz() - v2.Pz() * v1.Py();
          vpy = v2.Pz() * v1.Px() - v2.Px() * v1.Pz();
          vpz = v2.Px() * v1.Py() - v2.Py() * v1.Px();
        }
      }
    } else { // Unlike Sign
      if (!swapTracks) {
        if (d_bz * t1.sign() > 0) {
          vpx = v1.Py() * v2.Pz() - v1.Pz() * v2.Py();
          vpy = v1.Pz() * v2.Px() - v1.Px() * v2.Pz();
          vpz = v1.Px() * v2.Py() - v1.Py() * v2.Px();
        } else {
          vpx = v2.Py() * v1.Pz() - v2.Pz() * v1.Py();
          vpy = v2.Pz() * v1.Px() - v2.Px() * v1.Pz();
          vpz = v2.Px() * v1.Py() - v2.Py() * v1.Px();
        }
      } else { // swaped tracks
        if (d_bz * t2.sign() > 0) {
          vpx = v1.Py() * v2.Pz() - v1.Pz() * v2.Py();
          vpy = v1.Pz() * v2.Px() - v1.Px() * v2.Pz();
          vpz = v1.Px() * v2.Py() - v1.Py() * v2.Px();
        } else {
          vpx = v2.Py() * v1.Pz() - v2.Pz() * v1.Py();
          vpy = v2.Pz() * v1.Px() - v2.Px() * v1.Pz();
          vpz = v2.Px() * v1.Py() - v2.Py() * v1.Px();
        }
      }
    }

    // unit vector of pep X pem
    float vx = vpx / TMath::Sqrt(vpx * vpx + vpy * vpy + vpz * vpz);
    float vy = vpy / TMath::Sqrt(vpx * vpx + vpy * vpy + vpz * vpz);
    float vz = vpz / TMath::Sqrt(vpx * vpx + vpy * vpy + vpz * vpz);

    float px = v12.Px();
    float py = v12.Py();
    float pz = v12.Pz();

    // unit vector of (pep+pem)
    float ux = px / TMath::Sqrt(px * px + py * py + pz * pz);
    float uy = py / TMath::Sqrt(px * px + py * py + pz * pz);
    float uz = pz / TMath::Sqrt(px * px + py * py + pz * pz);
    float ax = uy / TMath::Sqrt(ux * ux + uy * uy);
    float ay = -ux / TMath::Sqrt(ux * ux + uy * uy);

    // The third axis defined by vector product (ux,uy,uz)X(vx,vy,vz)
    float wx = uy * vz - uz * vy;
    float wy = uz * vx - ux * vz;
    // by construction, (wx,wy,wz) must be a unit vector. Measure angle between (wx,wy,wz) and (ax,ay,0).
    // The angle between them should be small if the pair is conversion. This function then returns values close to pi!
    return TMath::ACos(wx * ax + wy * ay); // phiv in [0,pi] //cosPhiV = wx * ax + wy * ay;
  }

  template <typename TMCParticle1, typename TMCParticle2, typename TMCParticles>
  int FindCommonMotherFrom2Prongs(TMCParticle1 const& p1, TMCParticle2 const& p2, TMCParticles const& mcparticles)
  {
    if (!p1.has_mothers())
      return -1;
    if (!p2.has_mothers())
      return -1;

    // LOGF(info,"original motherid1 = %d , motherid2 = %d", p1.mothersIds()[0], p2.mothersIds()[0]);

    int motherid1 = p1.mothersIds()[0];
    auto mother1 = mcparticles.iteratorAt(motherid1);
    int mother1_pdg = mother1.pdgCode();

    int motherid2 = p2.mothersIds()[0];
    auto mother2 = mcparticles.iteratorAt(motherid2);
    int mother2_pdg = mother2.pdgCode();

    // LOGF(info,"motherid1 = %d , motherid2 = %d", motherid1, motherid2);

    if (motherid1 != motherid2)
      return -1;
    if (mother1_pdg != mother2_pdg)
      return -1;

    if (abs(mother1_pdg) != 22        // photon
        && abs(mother1_pdg) != 111    // pi0
        && abs(mother1_pdg) != 221    // eta
        && abs(mother1_pdg) != 331    // eta'
        && abs(mother1_pdg) != 113    // rho
        && abs(mother1_pdg) != 223    // omega
        && abs(mother1_pdg) != 333    // phi
        && abs(mother1_pdg) != 443    // Jpsi
        && abs(mother1_pdg) != 100443 // psi2S
    ) {
      return -1;
    }
    return motherid1;
  }

  template <typename TMCParticle1, typename TMCParticle2, typename TMCParticles>
  int IsHFULS(TMCParticle1 const& p1, TMCParticle2 const& p2, TMCParticles const& mcparticles)
  {
    // in total, 4 cases for ULS pairs
    // 0. prompt c->e+ and cbar->e-
    // 1. b->e- and bbar->e+ (different b and bbar)
    // 2. b->c->e+ and bbar->cbar->e- (different b and bbar)
    // 3. b->c->e+ and b->e- (1 same b (or bbar))
    if (!p1.has_mothers())
      return EM_HFeeType::kUndef;
    if (!p2.has_mothers())
      return EM_HFeeType::kUndef;

    int motherid_p1 = p1.mothersIds()[0];
    int motherid_p2 = p2.mothersIds()[0];
    if (motherid_p1 == motherid_p2) { // different mother
      return EM_HFeeType::kUndef;     // this never happens in correlated HF->ee decays
    }

    auto mother_p1 = mcparticles.iteratorAt(motherid_p1);
    auto mother_p2 = mcparticles.iteratorAt(motherid_p2);
    int mother1_pdg = mother_p1.pdgCode();
    int mother2_pdg = mother_p2.pdgCode();

    if (((500 < abs(mother1_pdg) && abs(mother1_pdg) < 599) || (5000 < abs(mother1_pdg) && abs(mother1_pdg) < 5999)) && ((500 < abs(mother2_pdg) && abs(mother2_pdg) < 599) || (5000 < abs(mother2_pdg) && abs(mother2_pdg) < 5999))) {
      return EM_HFeeType::kBe_Be; // bb->ee, decay type = 2
    }

    if (mother_p1.has_mothers() && mother_p2.has_mothers()) { // search for decay type 1,3,4
      int grand_motherid_p1 = mother_p1.mothersIds()[0];
      int grand_motherid_p2 = mother_p2.mothersIds()[0];
      auto grand_mother_p1 = mcparticles.iteratorAt(grand_motherid_p1);
      auto grand_mother_p2 = mcparticles.iteratorAt(grand_motherid_p2);
      int grand_mother1_pdg = grand_mother_p1.pdgCode();
      int grand_mother2_pdg = grand_mother_p2.pdgCode();

      if (((400 < abs(mother1_pdg) && abs(mother1_pdg) < 499) || (4000 < abs(mother1_pdg) && abs(mother1_pdg) < 4999)) && ((400 < abs(mother2_pdg) && abs(mother2_pdg) < 499) || (4000 < abs(mother2_pdg) && abs(mother2_pdg) < 4999))) { // mother is charm

        if (((500 < abs(grand_mother1_pdg) && abs(grand_mother1_pdg) < 599) || (5000 < abs(grand_mother1_pdg) && abs(grand_mother1_pdg) < 5999)) && ((500 < abs(grand_mother2_pdg) && abs(grand_mother2_pdg) < 599) || (5000 < abs(grand_mother2_pdg) && abs(grand_mother2_pdg) < 5999))) { // grand mother is beauty
          return kBCe_BCe;                                                                                                                                                                                                                                                                  // b->c->e and b->c->e, decay type = 1
        } else {
          return kCe_Ce; // prompt cc->ee, decay type = 0
        }
      }

      if (motherid_p1 == grand_motherid_p2 || grand_motherid_p1 == motherid_p2) {
        if (
          (((500 < abs(mother1_pdg) && abs(mother1_pdg) < 599) || (5000 < abs(mother1_pdg) && abs(mother1_pdg) < 5999)) && ((500 < abs(grand_mother2_pdg) && abs(grand_mother2_pdg) < 599) || (5000 < abs(grand_mother2_pdg) && abs(grand_mother2_pdg) < 5999))) ||
          (((500 < abs(mother2_pdg) && abs(mother2_pdg) < 599) || (5000 < abs(mother2_pdg) && abs(mother2_pdg) < 5999)) && ((500 < abs(grand_mother1_pdg) && abs(grand_mother1_pdg) < 599) || (5000 < abs(grand_mother1_pdg) && abs(grand_mother1_pdg) < 5999)))) {
          return EM_HFeeType::kBCe_Be_SameB; // b->c->e and c->e, decay type = 3
        }
      }
    }
    return EM_HFeeType::kUndef;
  }

  template <typename TMCParticle1, typename TMCParticle2, typename TMCParticles>
  int IsHFLS(TMCParticle1 const& p1, TMCParticle2 const& p2, TMCParticles const& mcparticles)
  {
    // in total, 1 case for LS pairs
    // 4. b->c->e+ and bbar->e+
    if (!p1.has_mothers())
      return EM_HFeeType::kUndef;
    if (!p2.has_mothers())
      return EM_HFeeType::kUndef;

    int motherid_p1 = p1.mothersIds()[0];
    int motherid_p2 = p2.mothersIds()[0];
    if (motherid_p1 == motherid_p2) { // different mother
      return EM_HFeeType::kUndef;     // this never happens in correlated HF->ee decays
    }

    auto mother_p1 = mcparticles.iteratorAt(motherid_p1);
    auto mother_p2 = mcparticles.iteratorAt(motherid_p2);
    int mother1_pdg = mother_p1.pdgCode();
    int mother2_pdg = mother_p2.pdgCode();

    if (mother_p1.has_mothers() && mother_p2.has_mothers()) { // search for decay type 4
      int grand_motherid_p1 = mother_p1.mothersIds()[0];
      int grand_motherid_p2 = mother_p2.mothersIds()[0];
      auto grand_mother_p1 = mcparticles.iteratorAt(grand_motherid_p1);
      auto grand_mother_p2 = mcparticles.iteratorAt(grand_motherid_p2);
      int grand_mother1_pdg = grand_mother_p1.pdgCode();
      int grand_mother2_pdg = grand_mother_p2.pdgCode();

      if (motherid_p1 != grand_motherid_p2 && grand_motherid_p1 != motherid_p2) { // different mother and grand mother
        if (
          (((500 < abs(mother1_pdg) && abs(mother1_pdg) < 599) || (5000 < abs(mother1_pdg) && abs(mother1_pdg) < 5999)) && ((500 < abs(grand_mother2_pdg) && abs(grand_mother2_pdg) < 599) || (5000 < abs(grand_mother2_pdg) && abs(grand_mother2_pdg) < 5999))) ||
          (((500 < abs(mother2_pdg) && abs(mother2_pdg) < 599) || (5000 < abs(mother2_pdg) && abs(mother2_pdg) < 5999)) && ((500 < abs(grand_mother1_pdg) && abs(grand_mother1_pdg) < 599) || (5000 < abs(grand_mother1_pdg) && abs(grand_mother1_pdg) < 5999)))) {
          return EM_HFeeType::kBCe_Be_DiffB; // b->c->e and c->e, decay type = 4
        }
      }
    }
    return EM_HFeeType::kUndef;
  }

  template <typename TTrack>
  bool IsSelected(TTrack const& track)
  {
    if (!track.hasITS()) {
      return false;
    }
    if (!track.hasTPC()) {
      return false;
    }
    if (track.tpcNClsCrossedRows() < mincrossedrows) {
      return false;
    }
    if (track.itsChi2NCl() < -1) { // if tracks are not reconstructed properly, chi2/ITSncls is set to -999;
      return false;
    }
    return true;
  }

  Filter trackFilter = nabs(o2::aod::track::eta) < maxeta && o2::aod::track::tpcChi2NCl < maxchi2tpc && nabs(o2::aod::track::dcaXY) < 3.f;
  using MyFilteredTracksMC = soa::Filtered<FullTracksExtMC>;
  Preslice<MyFilteredTracksMC> perCollision = aod::track::collisionId;

  void processSingleTrack(soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::McCollisionLabels> const& collisions, aod::BCsWithTimestamps const&, MyFilteredTracksMC const& tracks, aod::McParticles const& mctracks, aod::McCollisions const&)
  {
    for (auto& collision : collisions) {
      // TODO: investigate the collisions without corresponding mcCollision
      if (!collision.has_mcCollision()) {
        continue;
      }

      registry.fill(HIST("hEventCounter"), 1.0); // all

      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      auto mccollision = collision.mcCollision();
      mycollision(bc.globalBC(), bc.runNumber(), collision.posX(), collision.posY(), collision.posZ(), collision.numContrib(), collision.sel8(), d_bz,
                  mccollision.generatorsID(), mccollision.posX(), mccollision.posY(), mccollision.posZ(),
                  collision.multTPC(), collision.multFV0A(), collision.multFV0C(), collision.multFT0A(), collision.multFT0C(),
                  collision.multFDDA(), collision.multFDDC(), collision.multZNA(), collision.multZNC(), collision.multTracklets(), collision.multNTracksPV(), collision.multNTracksPVeta1());

      auto tracks_coll = tracks.sliceBy(perCollision, collision.globalIndex());
      for (auto& track : tracks_coll) {

        if (track.tpcNClsCrossedRows() < mincrossedrows) {
          continue;
        }
        int mother_pdg = 0;
        int grand_mother_pdg = 0;
        auto mcparticle = track.mcParticle_as<aod::McParticles>();

        if (mcparticle.has_mothers()) {
          auto mother = mcparticle.mothers_as<aod::McParticles>().front(); // first mother
          mother_pdg = mother.pdgCode();
          if (mother.has_mothers()) {
            auto grand_mother = mother.mothers_as<aod::McParticles>().front(); // first mother of mother
            grand_mother_pdg = grand_mother.pdgCode();
          }
        }

        mytrack(mycollision.lastIndex(),
                track.sign(), track.pt(), track.eta(), track.phi(), track.dcaXY(), track.dcaZ(), sqrt(track.cYY()), sqrt(track.cZZ()),
                track.tpcNClsFindable(), track.tpcNClsFindableMinusFound(), track.tpcNClsFindableMinusCrossedRows(),
                track.tpcNClsFound(), track.tpcNClsCrossedRows(),
                track.tpcChi2NCl(), track.tpcInnerParam(),
                track.tpcSignal(), track.tpcNSigmaEl(), track.tpcNSigmaMu(), track.tpcNSigmaPi(), track.tpcNSigmaKa(), track.tpcNSigmaPr(),
                track.beta(), track.tofNSigmaEl(), track.tofNSigmaMu(), track.tofNSigmaPi(), track.tofNSigmaKa(), track.tofNSigmaPr(),
                track.itsClusterMap(), track.itsChi2NCl(), track.detectorMap(),
                mcparticle.pt(), mcparticle.eta(), mcparticle.phi(),
                mcparticle.vx(), mcparticle.vy(), mcparticle.vz(),
                mcparticle.pdgCode(), mcparticle.isPhysicalPrimary(), mother_pdg, grand_mother_pdg);

      } // end of track loop
    }   // end of collision loop
  }     // end of process
  PROCESS_SWITCH(TreeCreatorElectronML, processSingleTrack, "produce ML input for single track level", false);

  Partition<MyFilteredTracksMC> posTracks = o2::aod::track::signed1Pt > 0.f;
  Partition<MyFilteredTracksMC> negTracks = o2::aod::track::signed1Pt < 0.f;

  void processPair(soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::McCollisionLabels> const& collisions, aod::BCsWithTimestamps const&, MyFilteredTracksMC const& tracks, aod::McParticles const& mctracks, aod::McCollisions const&)
  {
    std::map<uint64_t, int> fNewLabels;
    int fCounter = 0;

    for (auto& collision : collisions) {
      // TODO: investigate the collisions without corresponding mcCollision
      if (!collision.has_mcCollision()) {
        continue;
      }
      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      auto mccollision = collision.mcCollision();
      initCCDB(bc);

      registry.fill(HIST("hEventCounter"), 1.0); // all

      mycollision(bc.globalBC(), bc.runNumber(), collision.posX(), collision.posY(), collision.posZ(), collision.numContrib(), collision.sel8(), d_bz,
                  mccollision.generatorsID(), mccollision.posX(), mccollision.posY(), mccollision.posZ(),
                  collision.multTPC(), collision.multFV0A(), collision.multFV0C(), collision.multFT0A(), collision.multFT0C(),
                  collision.multFDDA(), collision.multFDDC(), collision.multZNA(), collision.multZNC(), collision.multTracklets(), collision.multNTracksPV(), collision.multNTracksPVeta1());

      auto negTracks_coll = negTracks->sliceByCached(o2::aod::track::collisionId, collision.globalIndex(), cache);
      auto posTracks_coll = posTracks->sliceByCached(o2::aod::track::collisionId, collision.globalIndex(), cache);

      for (auto& track : tracks) {
        if (!IsSelected(track)) {
          continue;
        }

        if (!track.has_mcParticle()) {
          continue; // If no MC particle is found, skip the track
        }

        auto mctrack = track.mcParticle_as<aod::McParticles>();

        if (abs(mctrack.pdgCode()) != 11) {
          continue;
        }

        if (!(fNewLabels.find(track.globalIndex()) != fNewLabels.end())) {
          fNewLabels[track.globalIndex()] = fCounter;

          int mother_pdg = 0;
          int grand_mother_pdg = 0;
          if (mctrack.has_mothers()) {
            auto mother = mctrack.mothers_as<aod::McParticles>().front(); // first mother
            mother_pdg = mother.pdgCode();
            if (mother.has_mothers()) {
              auto grand_mother = mother.mothers_as<aod::McParticles>().front(); // first mother of mother
              grand_mother_pdg = grand_mother.pdgCode();
            }
          }

          mytrack(mycollision.lastIndex(),
                  track.sign(), track.pt(), track.eta(), track.phi(), track.dcaXY(), track.dcaZ(), sqrt(track.cYY()), sqrt(track.cZZ()),
                  track.tpcNClsFindable(), track.tpcNClsFindableMinusFound(), track.tpcNClsFindableMinusCrossedRows(),
                  track.tpcNClsFound(), track.tpcNClsCrossedRows(),
                  track.tpcChi2NCl(), track.tpcInnerParam(),
                  track.tpcSignal(), track.tpcNSigmaEl(), track.tpcNSigmaMu(), track.tpcNSigmaPi(), track.tpcNSigmaKa(), track.tpcNSigmaPr(),
                  track.beta(), track.tofNSigmaEl(), track.tofNSigmaMu(), track.tofNSigmaPi(), track.tofNSigmaKa(), track.tofNSigmaPr(),
                  track.itsClusterMap(), track.itsChi2NCl(), track.detectorMap(),
                  mctrack.pt(), mctrack.eta(), mctrack.phi(),
                  mctrack.vx(), mctrack.vy(), mctrack.vz(),
                  mctrack.pdgCode(), mctrack.isPhysicalPrimary(), mother_pdg, grand_mother_pdg);

          fCounter++;
        }

      } // end of track loop

      for (auto& [ele, pos] : combinations(CombinationsFullIndexPolicy(negTracks_coll, posTracks_coll))) {
        if (!IsSelected(ele) || !IsSelected(pos)) {
          continue;
        }

        if (!ele.has_mcParticle() || !pos.has_mcParticle()) {
          continue;
        }

        auto mcpos = pos.mcParticle_as<aod::McParticles>();
        auto mcele = ele.mcParticle_as<aod::McParticles>();

        if (abs(mcele.pdgCode()) != 11 || abs(mcpos.pdgCode()) != 11) {
          continue;
        } // charge swap is accepted.

        ROOT::Math::PtEtaPhiMVector v1(pos.pt(), pos.eta(), pos.phi(), o2::constants::physics::MassElectron);
        ROOT::Math::PtEtaPhiMVector v2(ele.pt(), ele.eta(), ele.phi(), o2::constants::physics::MassElectron);
        ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;

        float phiv = get_phiv(ele, pos);
        float pair_dca_xy = sqrt((pow(pos.dcaXY() / sqrt(pos.cYY()), 2) + pow(ele.dcaXY() / sqrt(ele.cYY()), 2)) / 2.);
        float pair_dca_z = sqrt((pow(pos.dcaZ() / sqrt(pos.cZZ()), 2) + pow(ele.dcaZ() / sqrt(ele.cZZ()), 2)) / 2.);

        bool isSM = false;
        int common_mother_id = FindCommonMotherFrom2Prongs(mcpos, mcele, mctracks);
        if (common_mother_id > 0) {
          isSM = true;
        }

        int isHF = IsHFULS(mcpos, mcele, mctracks); // if isHF == true, pdgCode is set to 0, because this pair is correlated HF ee pair decayed from different 2 mothers. Check pdg code of legs.

        if (isSM) {
          auto mcpair = mctracks.iteratorAt(common_mother_id);
          mypair(mycollision.lastIndex(), fNewLabels[pos.globalIndex()], fNewLabels[ele.globalIndex()],
                 v12.M(), v12.Pt(), v12.Eta(), v12.Phi(), phiv, pair_dca_xy, pair_dca_z,
                 isSM, isHF, EM_EEPairType::kULS, mcpair.pdgCode(), mcpair.statusCode(), mcpair.flags(),
                 mcpair.vx(), mcpair.vy(), mcpair.vz());
          // LOGF(info, "mcpair.pdgCode() = %d", mcpair.pdgCode());
        } else if (isHF > -1) {
          // isSM and isHF are not satisfied at the same time.
          mypair(mycollision.lastIndex(), fNewLabels[pos.globalIndex()], fNewLabels[ele.globalIndex()],
                 v12.M(), v12.Pt(), v12.Eta(), v12.Phi(), phiv, pair_dca_xy, pair_dca_z,
                 isSM, isHF, EM_EEPairType::kULS, 0, 0, 0,
                 0, 0, 0);
        } else { // this is combinatorial bkg
          mypair(mycollision.lastIndex(), fNewLabels[pos.globalIndex()], fNewLabels[ele.globalIndex()],
                 v12.M(), v12.Pt(), v12.Eta(), v12.Phi(), phiv, pair_dca_xy, pair_dca_z,
                 isSM, isHF, EM_EEPairType::kULS, 0, 0, 0,
                 0, 0, 0);
        }

      } // end of uls pair loop

      for (auto& [pos1, pos2] : combinations(CombinationsStrictlyUpperIndexPolicy(posTracks_coll, posTracks_coll))) {
        if (!IsSelected(pos1) || !IsSelected(pos2)) {
          continue;
        }

        if (!pos1.has_mcParticle() || !pos2.has_mcParticle()) {
          continue;
        }

        auto mcpos1 = pos1.mcParticle_as<aod::McParticles>();
        auto mcpos2 = pos2.mcParticle_as<aod::McParticles>();

        if (abs(mcpos1.pdgCode()) != 11 || abs(mcpos2.pdgCode()) != 11) {
          continue;
        } // charge swap is accepted.

        ROOT::Math::PtEtaPhiMVector v1(pos1.pt(), pos1.eta(), pos1.phi(), o2::constants::physics::MassElectron);
        ROOT::Math::PtEtaPhiMVector v2(pos2.pt(), pos2.eta(), pos2.phi(), o2::constants::physics::MassElectron);
        ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;

        float phiv = get_phiv(pos1, pos2);
        float pair_dca_xy = sqrt((pow(pos2.dcaXY() / sqrt(pos2.cYY()), 2) + pow(pos1.dcaXY() / sqrt(pos1.cYY()), 2)) / 2.);
        float pair_dca_z = sqrt((pow(pos2.dcaZ() / sqrt(pos2.cZZ()), 2) + pow(pos1.dcaZ() / sqrt(pos1.cZZ()), 2)) / 2.);

        bool isSM = false;
        int isHF = IsHFLS(mcpos1, mcpos2, mctracks);
        if (isHF > -1) {
          // isSM and isHF are not satisfied at the same time.
          mypair(mycollision.lastIndex(), fNewLabels[pos1.globalIndex()], fNewLabels[pos2.globalIndex()],
                 v12.M(), v12.Pt(), v12.Eta(), v12.Phi(), phiv, pair_dca_xy, pair_dca_z,
                 isSM, isHF, EM_EEPairType::kLSpp, 0, 0, 0,
                 0, 0, 0);
        } else { // this is combinatorial bkg
          mypair(mycollision.lastIndex(), fNewLabels[pos1.globalIndex()], fNewLabels[pos2.globalIndex()],
                 v12.M(), v12.Pt(), v12.Eta(), v12.Phi(), phiv, pair_dca_xy, pair_dca_z,
                 isSM, isHF, EM_EEPairType::kLSpp, 0, 0, 0,
                 0, 0, 0);
        }

      } // end of lspp pair loop

      for (auto& [ele1, ele2] : combinations(CombinationsStrictlyUpperIndexPolicy(negTracks_coll, negTracks_coll))) {
        if (!IsSelected(ele1) || !IsSelected(ele2)) {
          continue;
        }
        if (!ele1.has_mcParticle() || !ele2.has_mcParticle()) {
          continue;
        }

        auto mcele1 = ele1.mcParticle_as<aod::McParticles>();
        auto mcele2 = ele2.mcParticle_as<aod::McParticles>();

        if (abs(mcele1.pdgCode()) != 11 || abs(mcele2.pdgCode()) != 11) {
          continue;
        } // charge swap is accepted.

        ROOT::Math::PtEtaPhiMVector v1(ele1.pt(), ele1.eta(), ele1.phi(), o2::constants::physics::MassElectron);
        ROOT::Math::PtEtaPhiMVector v2(ele2.pt(), ele2.eta(), ele2.phi(), o2::constants::physics::MassElectron);
        ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;

        float phiv = get_phiv(ele1, ele2);
        float pair_dca_xy = sqrt((pow(ele2.dcaXY() / sqrt(ele2.cYY()), 2) + pow(ele1.dcaXY() / sqrt(ele1.cYY()), 2)) / 2.);
        float pair_dca_z = sqrt((pow(ele2.dcaZ() / sqrt(ele2.cZZ()), 2) + pow(ele1.dcaZ() / sqrt(ele1.cZZ()), 2)) / 2.);

        bool isSM = false;
        int isHF = IsHFLS(mcele1, mcele2, mctracks);
        if (isHF > -1) {
          // isSM and isHF are not satisfied at the same time.
          mypair(mycollision.lastIndex(), fNewLabels[ele1.globalIndex()], fNewLabels[ele2.globalIndex()],
                 v12.M(), v12.Pt(), v12.Eta(), v12.Phi(), phiv, pair_dca_xy, pair_dca_z,
                 isSM, isHF, EM_EEPairType::kLSnn, 0, 0, 0,
                 0, 0, 0);
        } else { // this is combinatorial bkg
          mypair(mycollision.lastIndex(), fNewLabels[ele1.globalIndex()], fNewLabels[ele2.globalIndex()],
                 v12.M(), v12.Pt(), v12.Eta(), v12.Phi(), phiv, pair_dca_xy, pair_dca_z,
                 isSM, isHF, EM_EEPairType::kLSnn, 0, 0, 0,
                 0, 0, 0);
        }

      } // end of lsnn pair loop

    } // end of collision loop
    fNewLabels.clear();
    fCounter = 0;
  } // end of process
  PROCESS_SWITCH(TreeCreatorElectronML, processPair, "produce ML input for pair level", false);

  void processDummy(soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::McCollisionLabels> const& collisions) {}
  PROCESS_SWITCH(TreeCreatorElectronML, processDummy, "process dummy", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<TreeCreatorElectronML>(cfgc, TaskName{"tree-creator-ele-ml"})};
}
