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
DECLARE_SOA_INDEX_COLUMN(MyCollision, mycollision);                   //!
DECLARE_SOA_COLUMN(Sign, sign, int);                                  //!
DECLARE_SOA_COLUMN(TPCNClsFound, tpcNClsFound, int);                  //!
DECLARE_SOA_COLUMN(TPCNClsCrossedRows, tpcNClsCrossedRows, int);      //!
DECLARE_SOA_COLUMN(DCAresXY, dcaresXY, float);                        //!
DECLARE_SOA_COLUMN(DCAresZ, dcaresZ, float);                          //!
DECLARE_SOA_COLUMN(MCVx, mcvx, float);                                //!
DECLARE_SOA_COLUMN(MCVy, mcvy, float);                                //!
DECLARE_SOA_COLUMN(MCVz, mcvz, float);                                //!
DECLARE_SOA_COLUMN(IsPhysicalPrimary, isPhysicalPrimary, bool);       //!
DECLARE_SOA_COLUMN(MotherIds, motherIds, std::vector<int>);           //! eta values of the matched tracks
DECLARE_SOA_COLUMN(MotherPdgCodes, motherpdgCodes, std::vector<int>); //! eta values of the matched tracks
} // namespace mytrack

// reconstructed track information
DECLARE_SOA_TABLE(MyTracks, "AOD", "MYTRACK", //!
                  o2::soa::Index<>, mytrack::MyCollisionId, mytrack::Sign,
                  track::Pt, track::Eta, track::Phi,
                  track::DcaXY, track::DcaZ, mytrack::DCAresXY, mytrack::DCAresZ,
                  track::TPCNClsFindable, mytrack::TPCNClsFound, mytrack::TPCNClsCrossedRows,
                  track::TPCChi2NCl, track::TPCInnerParam,
                  track::TPCSignal, pidtpc::TPCNSigmaEl, pidtpc::TPCNSigmaMu, pidtpc::TPCNSigmaPi, pidtpc::TPCNSigmaKa, pidtpc::TPCNSigmaPr,
                  pidtofbeta::Beta, pidtof::TOFNSigmaEl, pidtof::TOFNSigmaMu, pidtof::TOFNSigmaPi, pidtof::TOFNSigmaKa, pidtof::TOFNSigmaPr,
                  track::ITSClusterMap, track::ITSChi2NCl,
                  mytrack::MCVx, mytrack::MCVy, mytrack::MCVz,
                  mcparticle::PdgCode, mytrack::IsPhysicalPrimary, mytrack::MotherIds, mytrack::MotherPdgCodes);

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

DECLARE_SOA_COLUMN(IsSM, isSM, bool);         //!
DECLARE_SOA_COLUMN(IsHF, isHF, int);          //!
DECLARE_SOA_COLUMN(PairType, pairtype, int);  //!
DECLARE_SOA_COLUMN(IsPrompt, isPrompt, bool); //!
} // namespace mypair

// reconstructed track information
DECLARE_SOA_TABLE(MyPairs, "AOD", "MYPAIR", //!
                  o2::soa::Index<>, mypair::MyCollisionId, mypair::PosTrackId, mypair::NegTrackId,
                  mypair::M, mypair::Pt, mypair::Eta, mypair::Phi, mypair::PhiV, mypair::PairDCAxy, mypair::PairDCAz,
                  mypair::IsSM, mypair::IsHF, mypair::PairType, mypair::IsPrompt,
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
  int IsHF(TMCParticle1 const& p1, TMCParticle2 const& p2, TMCParticles const& mcparticles)
  {
    // in total, 4 cases for ULS pairs
    // 0. prompt c->e+ and cbar->e-
    // 1. b->e- and bbar->e+ (different b and bbar)
    // 2. b->c->e+ and bbar->cbar->e- (different b and bbar)
    // 3. b->c->e+ and b->e- (1 same b (or bbar))
    if (!p1.has_mothers() || !p2.has_mothers()) {
      return EM_HFeeType::kUndef;
    }

    if (p1.mothersIds()[0] == p2.mothersIds()[0]) { // same mother
      return EM_HFeeType::kUndef;                   // this never happens in correlated HF->ee decays
    }

    // store all mother1 relation
    std::vector<int> mothers_id1;
    std::vector<int> mothers_pdg1;
    int motherid1 = p1.mothersIds()[0]; // first mother index
    while (motherid1 > -1) {
      if (motherid1 < mcparticles.size()) { // protect against bad mother indices. why is this needed?
        auto mp = mcparticles.iteratorAt(motherid1);
        mothers_id1.emplace_back(motherid1);
        mothers_pdg1.emplace_back(mp.pdgCode());

        if (mp.has_mothers()) {
          motherid1 = mp.mothersIds()[0];
        } else {
          motherid1 = -999;
        }
      } else {
        LOGF(info, "Mother label(%d) exceeds the McParticles size(%d)", motherid1, mcparticles.size());
      }
    }

    // store all mother2 relation
    std::vector<int> mothers_id2;
    std::vector<int> mothers_pdg2;
    int motherid2 = p2.mothersIds()[0]; // first mother index
    while (motherid2 > -1) {
      if (motherid2 < mcparticles.size()) { // protect against bad mother indices. why is this needed?
        auto mp = mcparticles.iteratorAt(motherid2);
        mothers_id2.emplace_back(motherid2);
        mothers_pdg2.emplace_back(mp.pdgCode());

        if (mp.has_mothers()) {
          motherid2 = mp.mothersIds()[0];
        } else {
          motherid2 = -999;
        }
      } else {
        LOGF(info, "Mother label(%d) exceeds the McParticles size(%d)", motherid2, mcparticles.size());
      }
    }

    if (std::to_string(mothers_pdg1[0]).find("5") != std::string::npos && std::to_string(mothers_pdg2[0]).find("5") != std::string::npos) {
      return EM_HFeeType::kBe_Be; // bb->ee, decay type = 2
      // this is easy. first mother is b hadron for both leg.
    }

    if (std::to_string(mothers_pdg1[0]).find("4") != std::string::npos && std::to_string(mothers_pdg2[0]).find("4") != std::string::npos) {
      // mother is c hadron. next, check c is prompt or non-prompt.

      bool is_c_from_b1 = false;
      for (unsigned int i1 = 1; i1 < mothers_pdg1.size(); i1++) {
        if (std::to_string(mothers_pdg1[i1]).find("5") != std::string::npos) {
          is_c_from_b1 = true;
          break;
        }
      }
      bool is_c_from_b2 = false;
      for (unsigned int i2 = 1; i2 < mothers_pdg2.size(); i2++) {
        if (std::to_string(mothers_pdg2[i2]).find("5") != std::string::npos) {
          is_c_from_b2 = true;
          break;
        }
      }

      if (!is_c_from_b1 && !is_c_from_b2) {
        return EM_HFeeType::kCe_Ce; // prompt cc->ee, decay type = 0
      } else if (is_c_from_b1 && is_c_from_b2) {
        return EM_HFeeType::kBCe_BCe; // b->c->e and b->c->e, decay type = 1
      } else {
        for (auto& mid1 : mothers_id1) {
          for (auto& mid2 : mothers_id2) {
            if (mid1 == mid2) {
              return EM_HFeeType::kBCe_Be_SameB; // b->c->e and c->e, decay type = 3. this should happen only in ULS.
            }
          }                                // end of mother id 2
        }                                  // end of mother id 1
        return EM_HFeeType::kBCe_Be_DiffB; // b->c->e and c->e, decay type = 4. this should happen only in LS. But, this may happen, when ele/pos is reconstructed as pos/ele wrongly. and create LS pair
      }
    }
    mothers_id1.shrink_to_fit();
    mothers_pdg1.shrink_to_fit();
    mothers_id2.shrink_to_fit();
    mothers_pdg2.shrink_to_fit();
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

  Filter trackFilter = nabs(o2::aod::track::eta) < maxeta && o2::aod::track::tpcChi2NCl < maxchi2tpc && nabs(o2::aod::track::dcaXY) < 1.f && nabs(o2::aod::track::dcaZ) < 1.f;
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

        if (!IsSelected(track)) {
          continue;
        }

        if (!track.has_mcParticle()) {
          continue; // If no MC particle is found, skip the track
        }
        auto mctrack = track.mcParticle_as<aod::McParticles>();

        // store all mother relation
        std::vector<int> mothers_id;
        std::vector<int> mothers_pdg;
        int motherid = mctrack.mothersIds()[0]; // first mother index
        while (motherid > -1) {
          if (motherid < mctracks.size()) { // protect against bad mother indices. why is this needed?
            auto mp = mctracks.iteratorAt(motherid);
            mothers_id.emplace_back(motherid);
            mothers_pdg.emplace_back(mp.pdgCode());

            if (mp.has_mothers()) {
              motherid = mp.mothersIds()[0];
            } else {
              motherid = -999;
            }
          } else {
            LOGF(info, "Mother label(%d) exceeds the McParticles size(%d)", motherid, mctracks.size());
          }
        }

        mytrack(mycollision.lastIndex(),
                track.sign(), track.pt(), track.eta(), track.phi(), track.dcaXY(), track.dcaZ(), sqrt(track.cYY()), sqrt(track.cZZ()),
                track.tpcNClsFindable(), track.tpcNClsFound(), track.tpcNClsCrossedRows(),
                track.tpcChi2NCl(), track.tpcInnerParam(),
                track.tpcSignal(), track.tpcNSigmaEl(), track.tpcNSigmaMu(), track.tpcNSigmaPi(), track.tpcNSigmaKa(), track.tpcNSigmaPr(),
                track.beta(), track.tofNSigmaEl(), track.tofNSigmaMu(), track.tofNSigmaPi(), track.tofNSigmaKa(), track.tofNSigmaPr(),
                track.itsClusterMap(), track.itsChi2NCl(),
                mctrack.vx(), mctrack.vy(), mctrack.vz(),
                mctrack.pdgCode(), mctrack.isPhysicalPrimary(), mothers_id, mothers_pdg);

        mothers_id.shrink_to_fit();
        mothers_pdg.shrink_to_fit();

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

          // store all mother relation
          std::vector<int> mothers_id;
          std::vector<int> mothers_pdg;
          int motherid = mctrack.mothersIds()[0]; // first mother index
          while (motherid > -1) {
            if (motherid < mctracks.size()) { // protect against bad mother indices. why is this needed?
              auto mp = mctracks.iteratorAt(motherid);
              mothers_id.emplace_back(motherid);
              mothers_pdg.emplace_back(mp.pdgCode());

              if (mp.has_mothers()) {
                motherid = mp.mothersIds()[0];
              } else {
                motherid = -999;
              }
            } else {
              LOGF(info, "Mother label(%d) exceeds the McParticles size(%d)", motherid, mctracks.size());
            }
          }

          mytrack(mycollision.lastIndex(),
                  track.sign(), track.pt(), track.eta(), track.phi(), track.dcaXY(), track.dcaZ(), sqrt(track.cYY()), sqrt(track.cZZ()),
                  track.tpcNClsFindable(), track.tpcNClsFound(), track.tpcNClsCrossedRows(),
                  track.tpcChi2NCl(), track.tpcInnerParam(),
                  track.tpcSignal(), track.tpcNSigmaEl(), track.tpcNSigmaMu(), track.tpcNSigmaPi(), track.tpcNSigmaKa(), track.tpcNSigmaPr(),
                  track.beta(), track.tofNSigmaEl(), track.tofNSigmaMu(), track.tofNSigmaPi(), track.tofNSigmaKa(), track.tofNSigmaPr(),
                  track.itsClusterMap(), track.itsChi2NCl(),
                  mctrack.vx(), mctrack.vy(), mctrack.vz(),
                  mctrack.pdgCode(), mctrack.isPhysicalPrimary(), mothers_id, mothers_pdg);

          mothers_id.shrink_to_fit();
          mothers_pdg.shrink_to_fit();
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
        int isHF = IsHF(mcpos, mcele, mctracks); // if isHF == true, pdgCode is set to 0, because this pair is correlated HF ee pair decayed from different 2 mothers. Check pdg code of legs.

        if (isSM) {
          auto mcpair = mctracks.iteratorAt(common_mother_id);
          bool is_prompt = true; // only relevant for prompt jpsi

          int motherid_pair = mcpair.mothersIds()[0]; // first mother index
          while (motherid_pair > -1) {
            if (motherid_pair < mctracks.size()) { // protect against bad mother indices. why is this needed?
              auto mp = mctracks.iteratorAt(motherid_pair);

              if (std::to_string(mp.pdgCode()).find("5") != std::string::npos) {
                is_prompt = false;
                break;
              }

              if (mp.has_mothers()) {
                motherid_pair = mp.mothersIds()[0];
              } else {
                motherid_pair = -999;
              }
            } else {
              LOGF(info, "Mother label(%d) exceeds the McParticles size(%d)", motherid_pair, mctracks.size());
            }
          }

          mypair(mycollision.lastIndex(), fNewLabels[pos.globalIndex()], fNewLabels[ele.globalIndex()],
                 v12.M(), v12.Pt(), v12.Eta(), v12.Phi(), phiv, pair_dca_xy, pair_dca_z,
                 isSM, isHF, EM_EEPairType::kULS, is_prompt, mcpair.pdgCode(), mcpair.statusCode(), mcpair.flags(),
                 mcpair.vx(), mcpair.vy(), mcpair.vz());
          // LOGF(info, "mcpair.pdgCode() = %d", mcpair.pdgCode());
        } else if (isHF > -1) {
          // isSM and isHF are not satisfied at the same time.
          mypair(mycollision.lastIndex(), fNewLabels[pos.globalIndex()], fNewLabels[ele.globalIndex()],
                 v12.M(), v12.Pt(), v12.Eta(), v12.Phi(), phiv, pair_dca_xy, pair_dca_z,
                 isSM, isHF, EM_EEPairType::kULS, false, 0, 0, 0,
                 0, 0, 0);
        } else { // this is combinatorial bkg
          mypair(mycollision.lastIndex(), fNewLabels[pos.globalIndex()], fNewLabels[ele.globalIndex()],
                 v12.M(), v12.Pt(), v12.Eta(), v12.Phi(), phiv, pair_dca_xy, pair_dca_z,
                 isSM, isHF, EM_EEPairType::kULS, false, 0, 0, 0,
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
        int isHF = IsHF(mcpos1, mcpos2, mctracks);
        if (isHF != EM_HFeeType::kUndef) {
          // isSM and isHF are not satisfied at the same time.
          mypair(mycollision.lastIndex(), fNewLabels[pos1.globalIndex()], fNewLabels[pos2.globalIndex()],
                 v12.M(), v12.Pt(), v12.Eta(), v12.Phi(), phiv, pair_dca_xy, pair_dca_z,
                 isSM, isHF, EM_EEPairType::kLSpp, false, 0, 0, 0,
                 0, 0, 0);
        } else { // this is combinatorial bkg
          mypair(mycollision.lastIndex(), fNewLabels[pos1.globalIndex()], fNewLabels[pos2.globalIndex()],
                 v12.M(), v12.Pt(), v12.Eta(), v12.Phi(), phiv, pair_dca_xy, pair_dca_z,
                 isSM, isHF, EM_EEPairType::kLSpp, false, 0, 0, 0,
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
        int isHF = IsHF(mcele1, mcele2, mctracks);
        if (isHF != EM_HFeeType::kUndef) {
          // isSM and isHF are not satisfied at the same time.
          mypair(mycollision.lastIndex(), fNewLabels[ele1.globalIndex()], fNewLabels[ele2.globalIndex()],
                 v12.M(), v12.Pt(), v12.Eta(), v12.Phi(), phiv, pair_dca_xy, pair_dca_z,
                 isSM, isHF, EM_EEPairType::kLSnn, false, 0, 0, 0,
                 0, 0, 0);
        } else { // this is combinatorial bkg
          mypair(mycollision.lastIndex(), fNewLabels[ele1.globalIndex()], fNewLabels[ele2.globalIndex()],
                 v12.M(), v12.Pt(), v12.Eta(), v12.Phi(), phiv, pair_dca_xy, pair_dca_z,
                 isSM, isHF, EM_EEPairType::kLSnn, false, 0, 0, 0,
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
