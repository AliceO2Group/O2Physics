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

#include <vector>
#include <TMath.h>

#include "Common/Core/RecoDecay.h"
#include "Framework/AnalysisDataModel.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/CaloClusters.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/Qvectors.h"
#include "PWGJE/DataModel/EMCALClusters.h"

#ifndef PWGEM_PHOTONMESON_DATAMODEL_GAMMATABLES_H_
#define PWGEM_PHOTONMESON_DATAMODEL_GAMMATABLES_H_

namespace o2::aod
{

namespace emevent
{
DECLARE_SOA_COLUMN(CollisionId, collisionId, int);
DECLARE_SOA_COLUMN(NgammaPCM, ngpcm, int);
DECLARE_SOA_COLUMN(NgammaPHOS, ngphos, int);
DECLARE_SOA_COLUMN(NgammaEMC, ngemc, int);
DECLARE_SOA_COLUMN(NeeULS, neeuls, int);
DECLARE_SOA_COLUMN(NeeLSpp, neelspp, int);
DECLARE_SOA_COLUMN(NeeLSmm, neelsmm, int);
DECLARE_SOA_COLUMN(NmumuULS, nmumuuls, int);
DECLARE_SOA_COLUMN(NmumuLSpp, nmumulspp, int);
DECLARE_SOA_COLUMN(NmumuLSmm, nmumulsmm, int);
DECLARE_SOA_COLUMN(NcollsPerBC, ncollsPerBC, int);
DECLARE_SOA_COLUMN(Bz, bz, float);           //! kG
DECLARE_SOA_COLUMN(Q2xFT0M, q2xft0m, float); //! Qx for 2nd harmonics in FT0M (i.e. positive eta)
DECLARE_SOA_COLUMN(Q2yFT0M, q2yft0m, float); //! Qy for 2nd harmonics in FT0M (i.e. positive eta)
DECLARE_SOA_COLUMN(Q2xFT0A, q2xft0a, float); //! Qx for 2nd harmonics in FT0A (i.e. positive eta)
DECLARE_SOA_COLUMN(Q2yFT0A, q2yft0a, float); //! Qy for 2nd harmonics in FT0A (i.e. positive eta)
DECLARE_SOA_COLUMN(Q2xFT0C, q2xft0c, float); //! Qx for 2nd harmonics in FT0C (i.e. negative eta)
DECLARE_SOA_COLUMN(Q2yFT0C, q2yft0c, float); //! Qy for 2nd harmonics in FT0C (i.e. negative eta)
DECLARE_SOA_COLUMN(Q2xFV0A, q2xfv0a, float); //! Qx for 2nd harmonics in FV0A (i.e. positive eta)
DECLARE_SOA_COLUMN(Q2yFV0A, q2yfv0a, float); //! Qy for 2nd harmonics in FV0A (i.e. positive eta)
DECLARE_SOA_COLUMN(Q2xBPos, q2xbpos, float); //! Qx for 2nd harmonics in Barrel positive eta region
DECLARE_SOA_COLUMN(Q2yBPos, q2ybpos, float); //! Qy for 2nd harmonics in Barrel positive eta region
DECLARE_SOA_COLUMN(Q2xBNeg, q2xbneg, float); //! Qx for 2nd harmonics in Barrel negative eta region
DECLARE_SOA_COLUMN(Q2yBNeg, q2ybneg, float); //! Qy for 2nd harmonics in Barrel negative eta region
DECLARE_SOA_COLUMN(Q3xFT0M, q3xft0m, float); //! Qx for 3rd harmonics in FT0M (i.e. positive eta)
DECLARE_SOA_COLUMN(Q3yFT0M, q3yft0m, float); //! Qy for 3rd harmonics in FT0M (i.e. positive eta)
DECLARE_SOA_COLUMN(Q3xFT0A, q3xft0a, float); //! Qx for 3rd harmonics in FT0A (i.e. positive eta)
DECLARE_SOA_COLUMN(Q3yFT0A, q3yft0a, float); //! Qy for 3rd harmonics in FT0A (i.e. positive eta)
DECLARE_SOA_COLUMN(Q3xFT0C, q3xft0c, float); //! Qx for 3rd harmonics in FT0C (i.e. negative eta)
DECLARE_SOA_COLUMN(Q3yFT0C, q3yft0c, float); //! Qy for 3rd harmonics in FT0C (i.e. negative eta)
DECLARE_SOA_COLUMN(Q3xFV0A, q3xfv0a, float); //! Qx for 3rd harmonics in FV0A (i.e. positive eta)
DECLARE_SOA_COLUMN(Q3yFV0A, q3yfv0a, float); //! Qy for 3rd harmonics in FV0A (i.e. positive eta)
DECLARE_SOA_COLUMN(Q3xBPos, q3xbpos, float); //! Qx for 3rd harmonics in Barrel positive eta region
DECLARE_SOA_COLUMN(Q3yBPos, q3ybpos, float); //! Qy for 3rd harmonics in Barrel positive eta region
DECLARE_SOA_COLUMN(Q3xBNeg, q3xbneg, float); //! Qx for 3rd harmonics in Barrel negative eta region
DECLARE_SOA_COLUMN(Q3yBNeg, q3ybneg, float); //! Qy for 3rd harmonics in Barrel negative eta region

DECLARE_SOA_DYNAMIC_COLUMN(EP2FT0M, ep2ft0m, [](float q2x, float q2y) -> float { return std::atan2(q2y, q2x) / 2.0; });
DECLARE_SOA_DYNAMIC_COLUMN(EP2FT0A, ep2ft0a, [](float q2x, float q2y) -> float { return std::atan2(q2y, q2x) / 2.0; });
DECLARE_SOA_DYNAMIC_COLUMN(EP2FT0C, ep2ft0c, [](float q2x, float q2y) -> float { return std::atan2(q2y, q2x) / 2.0; });
DECLARE_SOA_DYNAMIC_COLUMN(EP2FV0A, ep2fv0a, [](float q2x, float q2y) -> float { return std::atan2(q2y, q2x) / 2.0; });
DECLARE_SOA_DYNAMIC_COLUMN(EP2BPos, ep2bpos, [](float q2x, float q2y) -> float { return std::atan2(q2y, q2x) / 2.0; });
DECLARE_SOA_DYNAMIC_COLUMN(EP2BNeg, ep2bneg, [](float q2x, float q2y) -> float { return std::atan2(q2y, q2x) / 2.0; });
} // namespace emevent
DECLARE_SOA_TABLE(EMEvents, "AOD", "EMEVENT", //!   Main event information table
                  o2::soa::Index<>, emevent::CollisionId, bc::GlobalBC, bc::RunNumber, evsel::Sel8, evsel::Alias, evsel::Selection, emevent::NcollsPerBC,
                  collision::PosX, collision::PosY, collision::PosZ,
                  collision::NumContrib, collision::CollisionTime, collision::CollisionTimeRes);
using EMEvent = EMEvents::iterator;

DECLARE_SOA_TABLE(EMEventsBz, "AOD", "EMEVENTBZ", emevent::Bz); // joinable to EMEvents
using EMEventBz = EMEventsBz::iterator;

DECLARE_SOA_TABLE(EMEventsMult, "AOD", "EMEVENTMULT", //!   event multiplicity table, joinable to EMEvents
                  mult::MultFV0A, mult::MultFV0C, mult::MultFT0A, mult::MultFT0C,
                  mult::MultFDDA, mult::MultFDDC, mult::MultZNA, mult::MultZNC,
                  mult::MultTPC, mult::MultTracklets, mult::MultNTracksPV, mult::MultNTracksPVeta1, mult::MultNTracksPVetaHalf,
                  mult::IsInelGt0<mult::MultNTracksPVeta1>, mult::IsInelGt1<mult::MultNTracksPVeta1>);
using EMEventMult = EMEventsMult::iterator;

DECLARE_SOA_TABLE(EMEventsCent, "AOD", "EMEVENTCENT", //!   event centrality table, joinable to EMEvents
                  cent::CentFT0M, cent::CentFT0A, cent::CentFT0C, cent::CentNTPV);
using EMEventCent = EMEventsCent::iterator;

DECLARE_SOA_TABLE(EMEventsQvec, "AOD", "EMEVENTQVEC", //!   event q vector table, joinable to EMEvents
                  emevent::Q2xFT0M, emevent::Q2yFT0M, emevent::Q2xFT0A, emevent::Q2yFT0A, emevent::Q2xFT0C, emevent::Q2yFT0C, emevent::Q2xFV0A, emevent::Q2yFV0A,
                  emevent::Q2xBPos, emevent::Q2yBPos, emevent::Q2xBNeg, emevent::Q2yBNeg,
                  emevent::Q3xFT0M, emevent::Q3yFT0M, emevent::Q3xFT0A, emevent::Q3yFT0A, emevent::Q3xFT0C, emevent::Q3yFT0C, emevent::Q3xFV0A, emevent::Q3yFV0A,
                  emevent::Q3xBPos, emevent::Q3yBPos, emevent::Q3xBNeg, emevent::Q3yBNeg,

                  // Dynamic columns
                  emevent::EP2FT0M<emevent::Q2xFT0M, emevent::Q2yFT0M>,
                  emevent::EP2FT0A<emevent::Q2xFT0A, emevent::Q2yFT0A>,
                  emevent::EP2FT0C<emevent::Q2xFT0C, emevent::Q2yFT0C>,
                  emevent::EP2FV0A<emevent::Q2xFV0A, emevent::Q2yFV0A>,
                  emevent::EP2BPos<emevent::Q2xBPos, emevent::Q2yBPos>,
                  emevent::EP2BNeg<emevent::Q2xBNeg, emevent::Q2yBNeg>);
using EMEventQvec = EMEventsQvec::iterator;

DECLARE_SOA_TABLE(EMEventsNgPCM, "AOD", "EMEVENTNGPCM", emevent::NgammaPCM); // joinable to EMEvents
using EMEventNgPCM = EMEventsNgPCM::iterator;

DECLARE_SOA_TABLE(EMEventsNgPHOS, "AOD", "EMEVENTNGPHOS", emevent::NgammaPHOS); // joinable to EMEvents
using EMEventNgPHOS = EMEventsNgPHOS::iterator;

DECLARE_SOA_TABLE(EMEventsNgEMC, "AOD", "EMEVENTNGEMC", emevent::NgammaEMC); // joinable to EMEvents
using EMEventNgEMC = EMEventsNgEMC::iterator;

DECLARE_SOA_TABLE(EMEventsNee, "AOD", "EMEVENTNEE", emevent::NeeULS, emevent::NeeLSpp, emevent::NeeLSmm); // joinable to EMEvents
using EMEventNee = EMEventsNee::iterator;

DECLARE_SOA_TABLE(EMEventsNmumu, "AOD", "EMEVENTNMUMU", emevent::NmumuULS, emevent::NmumuLSpp, emevent::NmumuLSmm); // joinable to EMEvents
using EMEventNmumu = EMEventsNmumu::iterator;

namespace emmcevent
{
DECLARE_SOA_COLUMN(McCollisionId, mcCollisionId, int);
}

DECLARE_SOA_TABLE(EMMCEvents, "AOD", "EMMCEVENT", //!   MC event information table
                  o2::soa::Index<>, emmcevent::McCollisionId, mccollision::GeneratorsID,
                  mccollision::PosX, mccollision::PosY, mccollision::PosZ,
                  mccollision::T, mccollision::ImpactParameter,

                  // dynamic column
                  mccollision::GetGeneratorId<mccollision::GeneratorsID>,
                  mccollision::GetSubGeneratorId<mccollision::GeneratorsID>,
                  mccollision::GetSourceId<mccollision::GeneratorsID>);

using EMMCEvent = EMMCEvents::iterator;

namespace emmceventlabel
{
DECLARE_SOA_INDEX_COLUMN(EMMCEvent, emmcevent); //! MC collision
DECLARE_SOA_COLUMN(McMask, mcMask, uint16_t);   //! Bit mask to indicate collision mismatches (bit ON means mismatch). Bit 15: indicates negative label
} // namespace emmceventlabel

DECLARE_SOA_TABLE(EMMCEventLabels, "AOD", "EMMCEVENTLABEL", //! Table joined to the EMEvents table containing the MC index
                  emmceventlabel::EMMCEventId, emmceventlabel::McMask);
using EMMCEventLabel = EMMCEventLabels::iterator;

namespace emmcparticle
{
DECLARE_SOA_INDEX_COLUMN(EMMCEvent, emmcevent);
DECLARE_SOA_SELF_ARRAY_INDEX_COLUMN(Mothers, mothers);     //! Mother tracks (possible empty) array. Iterate over mcParticle.mothers_as<aod::McParticles>())
DECLARE_SOA_SELF_ARRAY_INDEX_COLUMN(Daughters, daughters); //! Daughter tracks (possibly empty) array. Check for non-zero with mcParticle.has_daughters(). Iterate over mcParticle.daughters_as<aod::McParticles>())
DECLARE_SOA_DYNAMIC_COLUMN(Pt, pt, [](float px, float py) -> float { return RecoDecay::sqrtSumOfSquares(px, py); });
DECLARE_SOA_DYNAMIC_COLUMN(Eta, eta, [](float px, float py, float pz) -> float { return RecoDecay::eta(std::array{px, py, pz}); });
DECLARE_SOA_DYNAMIC_COLUMN(Phi, phi, [](float px, float py) -> float { return RecoDecay::phi(px, py); });
DECLARE_SOA_DYNAMIC_COLUMN(P, p, [](float px, float py, float pz) -> float { return RecoDecay::sqrtSumOfSquares(px, py, pz); });
DECLARE_SOA_DYNAMIC_COLUMN(Y, y, //! Particle rapidity
                           [](float pz, float e) -> float {
                             if ((e - pz) > static_cast<float>(1e-7)) {
                               return 0.5f * std::log((e + pz) / (e - pz));
                             } else {
                               return -999.0f;
                             }
                           });
} // namespace emmcparticle

// This table contains all MC truth tracks (both v0 and calos)
DECLARE_SOA_TABLE_FULL(EMMCParticles, "EMMCParticles", "AOD", "EMMCPARTICLE", //!  MC track information (on disk)
                       o2::soa::Index<>, emmcparticle::EMMCEventId,
                       mcparticle::PdgCode, mcparticle::Flags,
                       emmcparticle::MothersIds, emmcparticle::DaughtersIds,
                       mcparticle::Px, mcparticle::Py, mcparticle::Pz, mcparticle::E,
                       mcparticle::Vx, mcparticle::Vy, mcparticle::Vz, mcparticle::Vt,

                       // dynamic column
                       emmcparticle::Pt<mcparticle::Px, mcparticle::Py>,
                       emmcparticle::Eta<mcparticle::Px, mcparticle::Py, mcparticle::Pz>,
                       emmcparticle::Phi<mcparticle::Px, mcparticle::Py>,

                       emmcparticle::P<mcparticle::Px, mcparticle::Py, mcparticle::Pz>,
                       emmcparticle::Y<mcparticle::Pz, mcparticle::E>,
                       mcparticle::ProducedByGenerator<mcparticle::Flags>,
                       mcparticle::FromBackgroundEvent<mcparticle::Flags>,
                       mcparticle::IsPhysicalPrimary<mcparticle::Flags>);

using EMMCParticle = EMMCParticles::iterator;

namespace v0legmclabel
{
DECLARE_SOA_INDEX_COLUMN(EMMCParticle, emmcparticle); //!
DECLARE_SOA_COLUMN(McMask, mcMask, uint16_t);
} // namespace v0legmclabel

// NOTE: MC labels. This table has one entry for each reconstructed track (joinable with v0leg table)
DECLARE_SOA_TABLE(V0LegMCLabels, "AOD", "V0LEGMCLABEL", //!
                  v0legmclabel::EMMCParticleId, v0legmclabel::McMask);
using V0LegMCLabel = V0LegMCLabels::iterator;

namespace emprimaryelectronmclabel
{
DECLARE_SOA_INDEX_COLUMN(EMMCParticle, emmcparticle); //!
DECLARE_SOA_COLUMN(McMask, mcMask, uint16_t);
} // namespace emprimaryelectronmclabel

// NOTE: MC labels. This table has one entry for each reconstructed track (joinable with EMPrimaryElectrons table)
DECLARE_SOA_TABLE(EMPrimaryElectronMCLabels, "AOD", "EMPRMELMCLABEL", //!
                  emprimaryelectronmclabel::EMMCParticleId, emprimaryelectronmclabel::McMask);
using EMPrimaryElectronMCLabel = EMPrimaryElectronMCLabels::iterator;

namespace emprimarymuonmclabel
{
DECLARE_SOA_INDEX_COLUMN(EMMCParticle, emmcparticle); //!
DECLARE_SOA_COLUMN(McMask, mcMask, uint16_t);
} // namespace emprimarymuonmclabel

// NOTE: MC labels. This table has one entry for each reconstructed track (joinable with EMPrimaryMuons table)
DECLARE_SOA_TABLE(EMPrimaryMuonMCLabels, "AOD", "EMPRMMUMCLABEL", //!
                  emprimarymuonmclabel::EMMCParticleId, emprimarymuonmclabel::McMask);
using EMPrimaryMuonMCLabel = EMPrimaryMuonMCLabels::iterator;

// *  EMC cluster mc label tables:
// 1. EMCALMCClusters in EMCalClusters.h: Vectors of global mc particle ids and energy fractions of the cluster
// 2. EMCClusterMCLabels: Vector of global mc particle ids
// 3. EMEMCClusterMCLabels: EM MC particle ID of largest contributor to cluster
namespace emcclustermclabel
{
DECLARE_SOA_ARRAY_INDEX_COLUMN(EMMCParticle, emmcparticle); //!
} // namespace emcclustermclabel

// NOTE: MC labels. This table has one vector of global mc particle ids for each reconstructed emc cluster (joinable with emccluster table)
DECLARE_SOA_TABLE(EMCClusterMCLabels, "AOD", "EMCClsMCLABEL", //!
                  emcclustermclabel::EMMCParticleIds);
using EMCClusterMCLabel = EMCClusterMCLabels::iterator;

namespace ememcclustermclabel
{
DECLARE_SOA_INDEX_COLUMN(EMMCParticle, emmcparticle); //!
} // namespace ememcclustermclabel

// NOTE: MC labels. This table has one entry for each reconstructed emc cluster (joinable with emccluster table)
DECLARE_SOA_TABLE(EMEMCClusterMCLabels, "AOD", "EMEMCClsMCLABEL", //!
                  ememcclustermclabel::EMMCParticleId);
using EMEMCClusterMCLabel = EMEMCClusterMCLabels::iterator;

namespace v0leg
{
DECLARE_SOA_COLUMN(CollisionId, collisionId, int); //!
DECLARE_SOA_COLUMN(TrackId, trackId, int);         //!
DECLARE_SOA_COLUMN(Sign, sign, int8_t);            //!
DECLARE_SOA_COLUMN(Px, px, float);                 //! Px at SV
DECLARE_SOA_COLUMN(Py, py, float);                 //! Py at SV
DECLARE_SOA_COLUMN(Pz, pz, float);                 //! Pz at SV
DECLARE_SOA_DYNAMIC_COLUMN(P, p, [](float px, float py, float pz) -> float { return RecoDecay::sqrtSumOfSquares(px, py, pz); });
DECLARE_SOA_DYNAMIC_COLUMN(Pt, pt, [](float px, float py) -> float { return RecoDecay::sqrtSumOfSquares(px, py); });
DECLARE_SOA_DYNAMIC_COLUMN(Eta, eta, [](float px, float py, float pz) -> float { return RecoDecay::eta(std::array{px, py, pz}); });
DECLARE_SOA_DYNAMIC_COLUMN(Phi, phi, [](float px, float py) -> float { return RecoDecay::phi(px, py); });
DECLARE_SOA_DYNAMIC_COLUMN(MeanClusterSizeITS, meanClusterSizeITS, [](uint32_t itsClusterSizes) -> float {
  int total_cluster_size = 0, nl = 0;
  for (unsigned int layer = 0; layer < 7; layer++) {
    int cluster_size_per_layer = (itsClusterSizes >> (layer * 4)) & 0xf;
    if (cluster_size_per_layer > 0) {
      nl++;
    }
    total_cluster_size += cluster_size_per_layer;
  }
  if (nl > 0) {
    return static_cast<float>(total_cluster_size) / static_cast<float>(nl);
  } else {
    return 0;
  }
});
DECLARE_SOA_DYNAMIC_COLUMN(MeanClusterSizeITSib, meanClusterSizeITSib, [](uint32_t itsClusterSizes) -> float {
  int total_cluster_size = 0, nl = 0;
  for (unsigned int layer = 0; layer < 3; layer++) {
    int cluster_size_per_layer = (itsClusterSizes >> (layer * 4)) & 0xf;
    if (cluster_size_per_layer > 0) {
      nl++;
    }
    total_cluster_size += cluster_size_per_layer;
  }
  if (nl > 0) {
    return static_cast<float>(total_cluster_size) / static_cast<float>(nl);
  } else {
    return 0;
  }
});
DECLARE_SOA_DYNAMIC_COLUMN(MeanClusterSizeITSob, meanClusterSizeITSob, [](uint32_t itsClusterSizes) -> float {
  int total_cluster_size = 0, nl = 0;
  for (unsigned int layer = 3; layer < 7; layer++) {
    int cluster_size_per_layer = (itsClusterSizes >> (layer * 4)) & 0xf;
    if (cluster_size_per_layer > 0) {
      nl++;
    }
    total_cluster_size += cluster_size_per_layer;
  }
  if (nl > 0) {
    return static_cast<float>(total_cluster_size) / static_cast<float>(nl);
  } else {
    return 0;
  }
});
} // namespace v0leg
DECLARE_SOA_TABLE(V0Legs, "AOD", "V0LEG", //!
                  o2::soa::Index<>, v0leg::CollisionId, v0leg::TrackId, v0leg::Sign,
                  v0leg::Px, v0leg::Py, v0leg::Pz,
                  track::DcaXY, track::DcaZ,
                  track::TPCNClsFindable, track::TPCNClsFindableMinusFound, track::TPCNClsFindableMinusCrossedRows,
                  track::TPCChi2NCl, track::TPCInnerParam,
                  track::TPCSignal, pidtpc::TPCNSigmaEl, pidtpc::TPCNSigmaPi,
                  track::ITSClusterSizes, track::ITSChi2NCl, track::DetectorMap,
                  track::X, track::Y, track::Z, track::Tgl,

                  // dynamic column
                  v0leg::P<v0leg::Px, v0leg::Py, v0leg::Pz>,
                  v0leg::Pt<v0leg::Px, v0leg::Py>,
                  v0leg::Eta<v0leg::Px, v0leg::Py, v0leg::Pz>,
                  v0leg::Phi<v0leg::Px, v0leg::Py>,
                  track::TPCNClsFound<track::TPCNClsFindable, track::TPCNClsFindableMinusFound>,
                  track::TPCNClsCrossedRows<track::TPCNClsFindable, track::TPCNClsFindableMinusCrossedRows>,
                  track::TPCCrossedRowsOverFindableCls<track::TPCNClsFindable, track::TPCNClsFindableMinusCrossedRows>,
                  track::TPCFoundOverFindableCls<track::TPCNClsFindable, track::TPCNClsFindableMinusFound>,
                  track::v001::ITSClusterMap<track::ITSClusterSizes>, track::v001::ITSNCls<track::ITSClusterSizes>, track::v001::ITSNClsInnerBarrel<track::ITSClusterSizes>,
                  track::HasITS<track::DetectorMap>, track::HasTPC<track::DetectorMap>,
                  track::HasTRD<track::DetectorMap>, track::HasTOF<track::DetectorMap>,
                  v0leg::MeanClusterSizeITS<track::ITSClusterSizes>,
                  v0leg::MeanClusterSizeITSib<track::ITSClusterSizes>,
                  v0leg::MeanClusterSizeITSob<track::ITSClusterSizes>);

// iterators
using V0Leg = V0Legs::iterator;

namespace v0photonkf
{
DECLARE_SOA_INDEX_COLUMN(EMEvent, emevent);                             //!
DECLARE_SOA_COLUMN(CollisionId, collisionId, int);                      //!
DECLARE_SOA_INDEX_COLUMN_FULL(PosTrack, posTrack, int, V0Legs, "_Pos"); //!
DECLARE_SOA_INDEX_COLUMN_FULL(NegTrack, negTrack, int, V0Legs, "_Neg"); //!
DECLARE_SOA_COLUMN(Vx, vx, float);                                      //! secondary vertex x
DECLARE_SOA_COLUMN(Vy, vy, float);                                      //! secondary vertex y
DECLARE_SOA_COLUMN(Vz, vz, float);                                      //! secondary vertex z
DECLARE_SOA_COLUMN(Px, px, float);                                      //! px for photon kf
DECLARE_SOA_COLUMN(Py, py, float);                                      //! py for photon kf
DECLARE_SOA_COLUMN(Pz, pz, float);                                      //! pz for photon kf
DECLARE_SOA_COLUMN(MGamma, mGamma, float);                              //! invariant mass of dielectron at SV
DECLARE_SOA_COLUMN(DCAxyToPV, dcaXYtopv, float);                        //! DCAxy of V0 to PV
DECLARE_SOA_COLUMN(DCAzToPV, dcaZtopv, float);                          //! DCAz of V0 to PV
DECLARE_SOA_COLUMN(CosPA, cospa, float);                                //!
DECLARE_SOA_COLUMN(PCA, pca, float);                                    //!
DECLARE_SOA_COLUMN(Alpha, alpha, float);                                //!
DECLARE_SOA_COLUMN(QtArm, qtarm, float);                                //!
DECLARE_SOA_COLUMN(ChiSquareNDF, chiSquareNDF, float);                  // Chi2 / NDF of the reconstructed V0

DECLARE_SOA_DYNAMIC_COLUMN(E, e, [](float px, float py, float pz, float m = 0) -> float { return RecoDecay::sqrtSumOfSquares(px, py, pz, m); }); //! energy of v0 photn, mass to be given as argument when getter is called!
DECLARE_SOA_DYNAMIC_COLUMN(Pt, pt, [](float px, float py) -> float { return RecoDecay::sqrtSumOfSquares(px, py); });
DECLARE_SOA_DYNAMIC_COLUMN(Eta, eta, [](float px, float py, float pz) -> float { return RecoDecay::eta(std::array{px, py, pz}); });
DECLARE_SOA_DYNAMIC_COLUMN(Phi, phi, [](float px, float py) -> float { return RecoDecay::phi(px, py); });
DECLARE_SOA_DYNAMIC_COLUMN(P, p, [](float px, float py, float pz) -> float { return RecoDecay::sqrtSumOfSquares(px, py, pz); });
DECLARE_SOA_DYNAMIC_COLUMN(V0Radius, v0radius, [](float vx, float vy) -> float { return RecoDecay::sqrtSumOfSquares(vx, vy); });
} // namespace v0photonkf
DECLARE_SOA_TABLE(V0PhotonsKF, "AOD", "V0PHOTONKF", //!
                  o2::soa::Index<>, v0photonkf::CollisionId, v0photonkf::PosTrackId, v0photonkf::NegTrackId,
                  v0photonkf::Vx, v0photonkf::Vy, v0photonkf::Vz,
                  v0photonkf::Px, v0photonkf::Py, v0photonkf::Pz,
                  v0photonkf::MGamma,
                  v0photonkf::DCAxyToPV, v0photonkf::DCAzToPV,
                  v0photonkf::CosPA, v0photonkf::PCA,
                  v0photonkf::Alpha, v0photonkf::QtArm,
                  v0photonkf::ChiSquareNDF,

                  // dynamic column
                  v0photonkf::E<v0photonkf::Px, v0photonkf::Py, v0photonkf::Pz>,
                  v0photonkf::Pt<v0photonkf::Px, v0photonkf::Py>,
                  v0photonkf::Eta<v0photonkf::Px, v0photonkf::Py, v0photonkf::Pz>,
                  v0photonkf::Phi<v0photonkf::Px, v0photonkf::Py>,
                  v0photonkf::P<v0photonkf::Px, v0photonkf::Py, v0photonkf::Pz>,
                  v0photonkf::V0Radius<v0photonkf::Vx, v0photonkf::Vy>);
// iterators
using V0PhotonKF = V0PhotonsKF::iterator;

DECLARE_SOA_TABLE(V0KFEMEventIds, "AOD", "V0KFEMEVENTID", v0photonkf::EMEventId); // To be joined with V0PhotonsKF table at analysis level.
// iterators
using V0KFEMEventId = V0KFEMEventIds::iterator;

namespace emprimaryelectron
{
DECLARE_SOA_INDEX_COLUMN(EMEvent, emevent);        //!
DECLARE_SOA_COLUMN(CollisionId, collisionId, int); //!
DECLARE_SOA_COLUMN(TrackId, trackId, int);         //!
DECLARE_SOA_COLUMN(Sign, sign, int8_t);            //!
DECLARE_SOA_COLUMN(PrefilterBit, pfb, uint8_t);    //!
DECLARE_SOA_DYNAMIC_COLUMN(Px, px, [](float pt, float phi) -> float { return pt * std::cos(phi); });
DECLARE_SOA_DYNAMIC_COLUMN(Py, py, [](float pt, float phi) -> float { return pt * std::sin(phi); });
DECLARE_SOA_DYNAMIC_COLUMN(Pz, pz, [](float pt, float eta) -> float { return pt * std::sinh(eta); });
DECLARE_SOA_DYNAMIC_COLUMN(DcaXYinSigma, dcaXYinSigma, [](float dcaXY, float cYY) -> float { return dcaXY / std::sqrt(cYY); });
DECLARE_SOA_DYNAMIC_COLUMN(DcaZinSigma, dcaZinSigma, [](float dcaZ, float cZZ) -> float { return dcaZ / std::sqrt(cZZ); });
DECLARE_SOA_DYNAMIC_COLUMN(Dca3DinSigma, dca3DinSigma, [](float dcaXY, float dcaZ, float cYY, float cZZ, float cZY) -> float {
  float det = cYY * cZZ - cZY * cZY; // determinant
  if (det < 0) {
    return 999.f;
  } else {
    return std::sqrt(std::abs((dcaXY * dcaXY * cZZ + dcaZ * dcaZ * cYY - 2. * dcaXY * dcaZ * cZY) / det / 2.)); // dca 3d in sigma
  }
});
DECLARE_SOA_DYNAMIC_COLUMN(MeanClusterSizeITS, meanClusterSizeITS, [](uint32_t itsClusterSizes) -> float {
  int total_cluster_size = 0, nl = 0;
  for (unsigned int layer = 0; layer < 7; layer++) {
    int cluster_size_per_layer = (itsClusterSizes >> (layer * 4)) & 0xf;
    if (cluster_size_per_layer > 0) {
      nl++;
    }
    total_cluster_size += cluster_size_per_layer;
  }
  if (nl > 0) {
    return static_cast<float>(total_cluster_size) / static_cast<float>(nl);
  } else {
    return 0;
  }
});
DECLARE_SOA_DYNAMIC_COLUMN(MeanClusterSizeITSib, meanClusterSizeITSib, [](uint32_t itsClusterSizes) -> float {
  int total_cluster_size = 0, nl = 0;
  for (unsigned int layer = 0; layer < 3; layer++) {
    int cluster_size_per_layer = (itsClusterSizes >> (layer * 4)) & 0xf;
    if (cluster_size_per_layer > 0) {
      nl++;
    }
    total_cluster_size += cluster_size_per_layer;
  }
  if (nl > 0) {
    return static_cast<float>(total_cluster_size) / static_cast<float>(nl);
  } else {
    return 0;
  }
});
DECLARE_SOA_DYNAMIC_COLUMN(MeanClusterSizeITSob, meanClusterSizeITSob, [](uint32_t itsClusterSizes) -> float {
  int total_cluster_size = 0, nl = 0;
  for (unsigned int layer = 3; layer < 7; layer++) {
    int cluster_size_per_layer = (itsClusterSizes >> (layer * 4)) & 0xf;
    if (cluster_size_per_layer > 0) {
      nl++;
    }
    total_cluster_size += cluster_size_per_layer;
  }
  if (nl > 0) {
    return static_cast<float>(total_cluster_size) / static_cast<float>(nl);
  } else {
    return 0;
  }
});
} // namespace emprimaryelectron
DECLARE_SOA_TABLE(EMPrimaryElectrons, "AOD", "EMPRIMARYEL", //!
                  o2::soa::Index<>, emprimaryelectron::CollisionId,
                  emprimaryelectron::TrackId, emprimaryelectron::Sign,
                  track::Pt, track::Eta, track::Phi, track::DcaXY, track::DcaZ,
                  track::TPCNClsFindable, track::TPCNClsFindableMinusFound, track::TPCNClsFindableMinusCrossedRows,
                  track::TPCChi2NCl, track::TPCInnerParam,
                  track::TPCSignal, pidtpc::TPCNSigmaEl, pidtpc::TPCNSigmaMu, pidtpc::TPCNSigmaPi, pidtpc::TPCNSigmaKa, pidtpc::TPCNSigmaPr,
                  pidtofbeta::Beta, pidtof::TOFNSigmaEl, pidtof::TOFNSigmaMu, pidtof::TOFNSigmaPi, pidtof::TOFNSigmaKa, pidtof::TOFNSigmaPr,
                  track::ITSClusterSizes, track::ITSChi2NCl, track::DetectorMap, track::Tgl, track::CYY, track::CZZ, track::CZY,

                  // dynamic column
                  track::TPCNClsFound<track::TPCNClsFindable, track::TPCNClsFindableMinusFound>,
                  track::TPCNClsCrossedRows<track::TPCNClsFindable, track::TPCNClsFindableMinusCrossedRows>,
                  track::TPCCrossedRowsOverFindableCls<track::TPCNClsFindable, track::TPCNClsFindableMinusCrossedRows>,
                  track::TPCFoundOverFindableCls<track::TPCNClsFindable, track::TPCNClsFindableMinusFound>,
                  track::v001::ITSClusterMap<track::ITSClusterSizes>, track::v001::ITSNCls<track::ITSClusterSizes>, track::v001::ITSNClsInnerBarrel<track::ITSClusterSizes>,
                  track::HasITS<track::DetectorMap>, track::HasTPC<track::DetectorMap>,
                  track::HasTRD<track::DetectorMap>, track::HasTOF<track::DetectorMap>,
                  emprimaryelectron::Px<track::Pt, track::Phi>,
                  emprimaryelectron::Py<track::Pt, track::Phi>,
                  emprimaryelectron::Pz<track::Pt, track::Eta>,
                  emprimaryelectron::DcaXYinSigma<track::DcaXY, track::CYY>,
                  emprimaryelectron::DcaZinSigma<track::DcaZ, track::CZZ>,
                  emprimaryelectron::Dca3DinSigma<track::DcaXY, track::DcaZ, track::CYY, track::CZZ, track::CZY>,
                  emprimaryelectron::MeanClusterSizeITS<track::ITSClusterSizes>,
                  emprimaryelectron::MeanClusterSizeITSib<track::ITSClusterSizes>,
                  emprimaryelectron::MeanClusterSizeITSob<track::ITSClusterSizes>);
// iterators
using EMPrimaryElectron = EMPrimaryElectrons::iterator;

DECLARE_SOA_TABLE(EMPrimaryElectronEMEventIds, "AOD", "PRMELEMEVENTID", emprimaryelectron::EMEventId); // To be joined with EMPrimaryElectrons table at analysis level.
// iterators
using EMPrimaryElectronEMEventId = EMPrimaryElectronEMEventIds::iterator;

DECLARE_SOA_TABLE(EMPrimaryElectronsPrefilterBit, "AOD", "PRMELEPFB", emprimaryelectron::PrefilterBit); // To be joined with EMPrimaryElectrons table at analysis level.
// iterators
using EMPrimaryElectronPrefilterBit = EMPrimaryElectronsPrefilterBit::iterator;

namespace dalitzee
{
DECLARE_SOA_INDEX_COLUMN(EMEvent, emevent);                                         //!
DECLARE_SOA_INDEX_COLUMN_FULL(PosTrack, posTrack, int, EMPrimaryElectrons, "_Pos"); //!
DECLARE_SOA_INDEX_COLUMN_FULL(NegTrack, negTrack, int, EMPrimaryElectrons, "_Neg"); //!
DECLARE_SOA_COLUMN(CollisionId, collisionId, int);                                  //!
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(Mass, mass, float);
DECLARE_SOA_COLUMN(Rapidity, rapidity, float);
DECLARE_SOA_COLUMN(PhiV, phiv, float);
DECLARE_SOA_COLUMN(OpeningAngle, opangle, float);
DECLARE_SOA_COLUMN(Sign, sign, int8_t);                                                                                                  //!
DECLARE_SOA_DYNAMIC_COLUMN(Energy, e, [](float pt, float eta, float m) { return RecoDecay::sqrtSumOfSquares(pt * std::cosh(eta), m); }); // e = sqrt(p*p + m*m)
} // namespace dalitzee
DECLARE_SOA_TABLE(DalitzEEs, "AOD", "DALITZEE", //!
                  o2::soa::Index<>, dalitzee::CollisionId, dalitzee::PosTrackId, dalitzee::NegTrackId,
                  dalitzee::Pt, dalitzee::Eta, dalitzee::Phi, dalitzee::Mass, dalitzee::Rapidity,
                  dalitzee::PhiV, dalitzee::OpeningAngle, dalitzee::Sign,
                  dalitzee::Energy<o2::aod::dalitzee::Pt, o2::aod::dalitzee::Eta, o2::aod::dalitzee::Mass>);
// iterators
using DalitzEE = DalitzEEs::iterator;

DECLARE_SOA_TABLE(DalitzEEEMEventIds, "AOD", "EEEMEVENTID", dalitzee::EMEventId); // To be joined with DalitzEEs table at analysis level.
// iterators
using DalitzEEEMEventId = DalitzEEEMEventIds::iterator;

namespace emprimarymuon
{
DECLARE_SOA_INDEX_COLUMN(EMEvent, emevent);        //!
DECLARE_SOA_COLUMN(CollisionId, collisionId, int); //!
DECLARE_SOA_COLUMN(TrackId, trackId, int);         //!
DECLARE_SOA_COLUMN(Sign, sign, int8_t);            //!
DECLARE_SOA_COLUMN(PrefilterBit, pfb, uint8_t);    //!
DECLARE_SOA_DYNAMIC_COLUMN(Px, px, [](float pt, float phi) -> float { return pt * std::cos(phi); });
DECLARE_SOA_DYNAMIC_COLUMN(Py, py, [](float pt, float phi) -> float { return pt * std::sin(phi); });
DECLARE_SOA_DYNAMIC_COLUMN(Pz, pz, [](float pt, float eta) -> float { return pt * std::sinh(eta); });
DECLARE_SOA_DYNAMIC_COLUMN(DcaXYinSigma, dcaXYinSigma, [](float dcaXY, float cYY) -> float { return dcaXY / std::sqrt(cYY); });
DECLARE_SOA_DYNAMIC_COLUMN(DcaZinSigma, dcaZinSigma, [](float dcaZ, float cZZ) -> float { return dcaZ / std::sqrt(cZZ); });
DECLARE_SOA_DYNAMIC_COLUMN(Dca3DinSigma, dca3DinSigma, [](float dcaXY, float dcaZ, float cYY, float cZZ, float cZY) -> float {
  float det = cYY * cZZ - cZY * cZY; // determinant
  if (det < 0) {
    return 999.f;
  } else {
    return std::sqrt(std::abs((dcaXY * dcaXY * cZZ + dcaZ * dcaZ * cYY - 2. * dcaXY * dcaZ * cZY) / det / 2.)); // dca 3d in sigma
  }
});
DECLARE_SOA_DYNAMIC_COLUMN(MeanClusterSizeITS, meanClusterSizeITS, [](uint32_t itsClusterSizes) -> float {
  int total_cluster_size = 0, nl = 0;
  for (unsigned int layer = 0; layer < 7; layer++) {
    int cluster_size_per_layer = (itsClusterSizes >> (layer * 4)) & 0xf;
    if (cluster_size_per_layer > 0) {
      nl++;
    }
    total_cluster_size += cluster_size_per_layer;
  }
  if (nl > 0) {
    return static_cast<float>(total_cluster_size) / static_cast<float>(nl);
  } else {
    return 0;
  }
});
DECLARE_SOA_DYNAMIC_COLUMN(MeanClusterSizeITSib, meanClusterSizeITSib, [](uint32_t itsClusterSizes) -> float {
  int total_cluster_size = 0, nl = 0;
  for (unsigned int layer = 0; layer < 3; layer++) {
    int cluster_size_per_layer = (itsClusterSizes >> (layer * 4)) & 0xf;
    if (cluster_size_per_layer > 0) {
      nl++;
    }
    total_cluster_size += cluster_size_per_layer;
  }
  if (nl > 0) {
    return static_cast<float>(total_cluster_size) / static_cast<float>(nl);
  } else {
    return 0;
  }
});
DECLARE_SOA_DYNAMIC_COLUMN(MeanClusterSizeITSob, meanClusterSizeITSob, [](uint32_t itsClusterSizes) -> float {
  int total_cluster_size = 0, nl = 0;
  for (unsigned int layer = 3; layer < 7; layer++) {
    int cluster_size_per_layer = (itsClusterSizes >> (layer * 4)) & 0xf;
    if (cluster_size_per_layer > 0) {
      nl++;
    }
    total_cluster_size += cluster_size_per_layer;
  }
  if (nl > 0) {
    return static_cast<float>(total_cluster_size) / static_cast<float>(nl);
  } else {
    return 0;
  }
});
} // namespace emprimarymuon
DECLARE_SOA_TABLE(EMPrimaryMuons, "AOD", "EMPRIMARYMU", //!
                  o2::soa::Index<>, emprimarymuon::CollisionId,
                  emprimarymuon::TrackId, emprimarymuon::Sign,
                  track::Pt, track::Eta, track::Phi, track::DcaXY, track::DcaZ,
                  track::TPCNClsFindable, track::TPCNClsFindableMinusFound, track::TPCNClsFindableMinusCrossedRows,
                  track::TPCChi2NCl, track::TPCInnerParam,
                  track::TPCSignal, pidtpc::TPCNSigmaEl, pidtpc::TPCNSigmaMu, pidtpc::TPCNSigmaPi, pidtpc::TPCNSigmaKa, pidtpc::TPCNSigmaPr,
                  pidtofbeta::Beta, pidtof::TOFNSigmaEl, pidtof::TOFNSigmaMu, pidtof::TOFNSigmaPi, pidtof::TOFNSigmaKa, pidtof::TOFNSigmaPr,
                  track::ITSClusterSizes, track::ITSChi2NCl, track::DetectorMap, track::Tgl, track::CYY, track::CZZ, track::CZY,

                  // dynamic column
                  track::TPCNClsFound<track::TPCNClsFindable, track::TPCNClsFindableMinusFound>,
                  track::TPCNClsCrossedRows<track::TPCNClsFindable, track::TPCNClsFindableMinusCrossedRows>,
                  track::TPCCrossedRowsOverFindableCls<track::TPCNClsFindable, track::TPCNClsFindableMinusCrossedRows>,
                  track::TPCFoundOverFindableCls<track::TPCNClsFindable, track::TPCNClsFindableMinusFound>,
                  track::v001::ITSClusterMap<track::ITSClusterSizes>, track::v001::ITSNCls<track::ITSClusterSizes>, track::v001::ITSNClsInnerBarrel<track::ITSClusterSizes>,
                  track::HasITS<track::DetectorMap>, track::HasTPC<track::DetectorMap>,
                  track::HasTRD<track::DetectorMap>, track::HasTOF<track::DetectorMap>,
                  emprimarymuon::Px<track::Pt, track::Phi>,
                  emprimarymuon::Py<track::Pt, track::Phi>,
                  emprimarymuon::Pz<track::Pt, track::Eta>,
                  emprimarymuon::DcaXYinSigma<track::DcaXY, track::CYY>,
                  emprimarymuon::DcaZinSigma<track::DcaZ, track::CZZ>,
                  emprimarymuon::Dca3DinSigma<track::DcaXY, track::DcaZ, track::CYY, track::CZZ, track::CZY>,
                  emprimarymuon::MeanClusterSizeITS<track::ITSClusterSizes>,
                  emprimarymuon::MeanClusterSizeITSib<track::ITSClusterSizes>,
                  emprimarymuon::MeanClusterSizeITSob<track::ITSClusterSizes>);
// iterators
using EMPrimaryMuon = EMPrimaryMuons::iterator;

DECLARE_SOA_TABLE(EMPrimaryMuonEMEventIds, "AOD", "PRMMUEMEVENTID", emprimarymuon::EMEventId); // To be joined with EMPrimaryMuons table at analysis level.
// iterators
using EMPrimaryMuonEMEventId = EMPrimaryMuonEMEventIds::iterator;

DECLARE_SOA_TABLE(EMPrimaryMuonsPrefilterBit, "AOD", "PRMMUPFB", emprimarymuon::PrefilterBit); // To be joined with EMPrimaryMuons table at analysis level.
// iterators
using EMPrimaryMuonPrefilterBit = EMPrimaryMuonsPrefilterBit::iterator;

namespace dalitzmumu
{
DECLARE_SOA_INDEX_COLUMN(EMEvent, emevent);                                     //!
DECLARE_SOA_INDEX_COLUMN_FULL(PosTrack, posTrack, int, EMPrimaryMuons, "_Pos"); //!
DECLARE_SOA_INDEX_COLUMN_FULL(NegTrack, negTrack, int, EMPrimaryMuons, "_Neg"); //!
DECLARE_SOA_COLUMN(CollisionId, collisionId, int);                              //!
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(Mass, mass, float);
DECLARE_SOA_COLUMN(Rapidity, rapidity, float);
DECLARE_SOA_COLUMN(PhiV, phiv, float);
DECLARE_SOA_COLUMN(OpeningAngle, opangle, float);
DECLARE_SOA_COLUMN(Sign, sign, int8_t);                                                                                                  //!
DECLARE_SOA_DYNAMIC_COLUMN(Energy, e, [](float pt, float eta, float m) { return RecoDecay::sqrtSumOfSquares(pt * std::cosh(eta), m); }); // e = sqrt(p*p + m*m)
} // namespace dalitzmumu
DECLARE_SOA_TABLE(DalitzMuMus, "AOD", "DALITZMUMU", //!
                  o2::soa::Index<>, dalitzmumu::CollisionId, dalitzmumu::PosTrackId, dalitzmumu::NegTrackId,
                  dalitzmumu::Pt, dalitzmumu::Eta, dalitzmumu::Phi, dalitzmumu::Mass, dalitzmumu::Rapidity,
                  dalitzmumu::PhiV, dalitzmumu::OpeningAngle, dalitzmumu::Sign,
                  dalitzmumu::Energy<o2::aod::dalitzmumu::Pt, o2::aod::dalitzmumu::Eta, o2::aod::dalitzmumu::Mass>);
// iterators
using DalitzMuMu = DalitzMuMus::iterator;

DECLARE_SOA_TABLE(DalitzMuMuEMEventIds, "AOD", "MUMUEMEVENTID", dalitzmumu::EMEventId); // To be joined with DalitzMuMus table at analysis level.
// iterators
using DalitzMuMuEMEventId = DalitzMuMuEMEventIds::iterator;

namespace pwgem::photon::swtinfo
{
DECLARE_SOA_INDEX_COLUMN(EMEvent, emevent);                                                                              //!
DECLARE_SOA_COLUMN(CollisionId, collisionId, int);                                                                       //!
DECLARE_SOA_INDEX_COLUMN_FULL(TriggerV0PhotonHighPt, triggerV0PhotonHighPt, int, V0PhotonsKF, "_TriggerV0PhotonHighPt"); //! high pT PCM trigger is fired by this v0 photon
DECLARE_SOA_INDEX_COLUMN_FULL(TriggerV0PhotonPair, triggerV0PhotonPair, int, V0PhotonsKF, "_TriggerV0PhotonPair");       //! PCM+EE trigger is fired by this v0 photon and dielectron
DECLARE_SOA_INDEX_COLUMN_FULL(TriggerDielectronPair, triggerDielectronPair, int, DalitzEEs, "_TriggerDielectronPair");   //! PCM+EE trigger is fired by this v0 photon and dielectron
} // namespace pwgem::photon::swtinfo
DECLARE_SOA_TABLE(EMSwtInfosPCM, "AOD", "SWTINFOPCM", //!
                  o2::soa::Index<>, pwgem::photon::swtinfo::CollisionId, pwgem::photon::swtinfo::TriggerV0PhotonHighPtId);
using EMSwtInfoPCM = EMSwtInfosPCM::iterator;

DECLARE_SOA_TABLE(EMSwtInfoPCMEMEventIds, "AOD", "SWTPCMEVENTID", pwgem::photon::swtinfo::EMEventId, o2::soa::Marker<1>); // To be joined with EMSwtInfosPCM table at analysis level.
// iterators
using EMSwtInfoPCMEMEventId = EMSwtInfoPCMEMEventIds::iterator;

DECLARE_SOA_TABLE(EMSwtInfosPair, "AOD", "SWTINFOPAIR", //!
                  o2::soa::Index<>, pwgem::photon::swtinfo::CollisionId, pwgem::photon::swtinfo::TriggerV0PhotonPairId, pwgem::photon::swtinfo::TriggerDielectronPairId);
using EMSwtInfoPair = EMSwtInfosPair::iterator;

DECLARE_SOA_TABLE(EMSwtInfoPairEMEventIds, "AOD", "SWTPAIREVENTID", pwgem::photon::swtinfo::EMEventId, o2::soa::Marker<2>); // To be joined with EMSwtInfosPair table at analysis level.
// iterators
using EMSwtInfoPairEMEventId = EMSwtInfoPairEMEventIds::iterator;

namespace MCTracksTrue
{
DECLARE_SOA_COLUMN(SameMother, sameMother, bool); // Do the tracks have the same mother particle?
} // namespace MCTracksTrue

DECLARE_SOA_TABLE(V0DaughterMcParticles, "AOD", "MCTRACKTRUE",
                  mcparticle::PdgCode,
                  mcparticle::Px,
                  mcparticle::Py,
                  mcparticle::Pz,
                  MCTracksTrue::SameMother);

// Index table to associate V0DaughterTracks to the corresponding V0MCDaughterParticles if an entrie exists. This table contains an index column which contains the corresponding row number or -1, if there is no corresponding V0MCDaughterParticles
// This table is one column that is joinable with V0DaughterTracks as far as i understand
namespace MCParticleTrueIndex
{
DECLARE_SOA_INDEX_COLUMN(V0DaughterMcParticle, v0DaughterMcParticle);
} // namespace MCParticleTrueIndex

// DECLARE_SOA_INDEX_TABLE_USER(MCTrackIndex, V0MCDaughterParticles, "MCTRACKINDEX", MCParticleTrueIndex::V0DaughterTrackId);
DECLARE_SOA_TABLE(MCParticleIndex, "AOD", "MCPARTICLEINDEX", MCParticleTrueIndex::V0DaughterMcParticleId);

namespace gammamctrue
{
DECLARE_SOA_COLUMN(P, p, float); //! Absolute momentum in GeV/c
} // namespace gammamctrue

DECLARE_SOA_TABLE(McDaughterTrue, "AOD", "MCDAUTRUE",
                  gammamctrue::P);

namespace gammamctrue
{
DECLARE_SOA_INDEX_COLUMN_FULL(V0PhotonKF, v0photonkf, int, V0PhotonsKF, "_V0Photon"); //!
DECLARE_SOA_COLUMN(Gamma, gamma, int64_t);                                            //! Used as reference for the daughters
DECLARE_SOA_COLUMN(NDaughters, nDaughters, int);                                      //! Number of daughters
DECLARE_SOA_COLUMN(Eta, eta, float);                                                  //! Pseudorapidity
DECLARE_SOA_COLUMN(Phi, phi, float);                                                  //! Angle phi in rad
DECLARE_SOA_COLUMN(Pt, pt, float);                                                    //! Transversal momentum in GeV/c
DECLARE_SOA_COLUMN(Y, y, float);                                                      //! Rapidity

DECLARE_SOA_COLUMN(ConversionX, conversionX, float); //! x of conversion point in cm
DECLARE_SOA_COLUMN(ConversionY, conversionY, float); //! y of conversion point in cm
DECLARE_SOA_COLUMN(ConversionZ, conversionZ, float); //! z of conversion point in cm
DECLARE_SOA_COLUMN(V0Radius, v0Radius, float);       //! 2d radius of conversion point
// DECLARE_SOA_INDEX_COLUMN(McDaughterTrue, mcDaughterTrue);
DECLARE_SOA_INDEX_COLUMN_FULL(McDaughterTrueOne, mcDaughterTrueOne, int, McDaughterTrue, "_One"); // this is a reference that points to the entry in the McDaughterTrues table
DECLARE_SOA_INDEX_COLUMN_FULL(McDaughterTrueTwo, mcDaughterTrueTwo, int, McDaughterTrue, "_Two"); // this is a reference that points to the entry in the McDaughterTrues table
} // namespace gammamctrue

DECLARE_SOA_TABLE(McGammasTrue, "AOD", "MCGATRUE",
                  o2::soa::Index<>,
                  mcparticle::McCollisionId,
                  gammamctrue::Gamma,
                  // v0data::V0Id, // reference to reconstructed v0 (if its a task with reconstucted info)
                  gammamctrue::V0PhotonKFId, // reference to reconstructed v0 (if its a task with reconstucted info)
                  mcparticle::PdgCode, mcparticle::StatusCode, mcparticle::Flags,
                  mcparticle::Px, mcparticle::Py, mcparticle::Pz,
                  mcparticle::Vx, mcparticle::Vy, mcparticle::Vz, mcparticle::Vt,
                  gammamctrue::NDaughters,
                  // todo: make those expression columns in an extended table
                  gammamctrue::Eta, gammamctrue::Phi, gammamctrue::P, gammamctrue::Pt, gammamctrue::Y,
                  gammamctrue::ConversionX, gammamctrue::ConversionY, gammamctrue::ConversionZ,
                  gammamctrue::V0Radius,

                  // Index columns
                  gammamctrue::McDaughterTrueOneId,
                  gammamctrue::McDaughterTrueTwoId,

                  // Dynamic columns
                  mcparticle::ProducedByGenerator<mcparticle::Flags>,
                  mcparticle::FromBackgroundEvent<mcparticle::Flags>,
                  mcparticle::GetGenStatusCode<mcparticle::Flags, mcparticle::StatusCode>,
                  mcparticle::GetProcess<mcparticle::Flags, mcparticle::StatusCode>,
                  mcparticle::IsPhysicalPrimary<mcparticle::Flags>);

namespace skimmedcluster
{
DECLARE_SOA_INDEX_COLUMN(Collision, collision);                        //! collisionID used as index for matched clusters
DECLARE_SOA_INDEX_COLUMN(BC, bc);                                      //! bunch crossing ID used as index for ambiguous clusters
DECLARE_SOA_COLUMN(ID, id, int);                                       //! cluster ID identifying cluster in event
DECLARE_SOA_COLUMN(E, e, float);                                       //! cluster energy (GeV)
DECLARE_SOA_COLUMN(Eta, eta, float);                                   //! cluster pseudorapidity (calculated using vertex)
DECLARE_SOA_COLUMN(Phi, phi, float);                                   //! cluster azimuthal angle (calculated using vertex)
DECLARE_SOA_COLUMN(M02, m02, float);                                   //! shower shape long axis
DECLARE_SOA_COLUMN(M20, m20, float);                                   //! shower shape short axis
DECLARE_SOA_COLUMN(NCells, nCells, int);                               //! number of cells in cluster
DECLARE_SOA_COLUMN(Time, time, float);                                 //! cluster time (ns)
DECLARE_SOA_COLUMN(DistanceToBadChannel, distanceToBadChannel, float); //! distance to bad channel
DECLARE_SOA_COLUMN(NLM, nlm, int);                                     //! number of local maxima
} // namespace skimmedcluster

namespace emccluster
{
DECLARE_SOA_INDEX_COLUMN(EMEvent, emevent);                                                                                   //!
DECLARE_SOA_COLUMN(CoreEnergy, coreEnergy, float);                                                                            //! cluster core energy (GeV)
DECLARE_SOA_COLUMN(Time, time, float);                                                                                        //! cluster time (ns)
DECLARE_SOA_COLUMN(IsExotic, isExotic, bool);                                                                                 //! flag to mark cluster as exotic
DECLARE_SOA_COLUMN(Definition, definition, int);                                                                              //! cluster definition, see EMCALClusterDefinition.h
DECLARE_SOA_ARRAY_INDEX_COLUMN(Track, track);                                                                                 //! TrackIds
DECLARE_SOA_COLUMN(TrackEta, tracketa, std::vector<float>);                                                                   //! eta values of the matched tracks
DECLARE_SOA_COLUMN(TrackPhi, trackphi, std::vector<float>);                                                                   //! phi values of the matched tracks
DECLARE_SOA_COLUMN(TrackP, trackp, std::vector<float>);                                                                       //! momentum values of the matched tracks
DECLARE_SOA_COLUMN(TrackPt, trackpt, std::vector<float>);                                                                     //! pt values of the matched tracks
DECLARE_SOA_DYNAMIC_COLUMN(Pt, pt, [](float e, float eta, float m = 0) -> float { return sqrt(e * e - m * m) / cosh(eta); }); //! cluster pt, mass to be given as argument when getter is called!
} // namespace emccluster
DECLARE_SOA_TABLE(SkimEMCClusters, "AOD", "SKIMEMCCLUSTERS", //! table of skimmed EMCal clusters
                  o2::soa::Index<>, skimmedcluster::CollisionId, skimmedcluster::BCId, skimmedcluster::E, emccluster::CoreEnergy,
                  skimmedcluster::Eta, skimmedcluster::Phi, skimmedcluster::M02, skimmedcluster::M20, skimmedcluster::NCells, skimmedcluster::Time,
                  emccluster::IsExotic, skimmedcluster::DistanceToBadChannel, skimmedcluster::NLM, emccluster::Definition,
                  emccluster::TrackIds, emccluster::TrackEta, emccluster::TrackPhi, emccluster::TrackP, emccluster::TrackPt,
                  // dynamic column
                  emccluster::Pt<skimmedcluster::E, skimmedcluster::Eta>);
using SkimEMCCluster = SkimEMCClusters::iterator;

DECLARE_SOA_TABLE(EMCEMEventIds, "AOD", "EMCEMEVENTID", emccluster::EMEventId); // To be joined with SkimEMCClusters table at analysis level.
// iterators
using EMCEMEventId = EMCEMEventIds::iterator;

namespace phoscluster
{
DECLARE_SOA_INDEX_COLUMN(EMEvent, emevent);                                         //!
DECLARE_SOA_INDEX_COLUMN_FULL(MatchedTrack, matchedTrack, int, Tracks, "_Matched"); //! matched track index
DECLARE_SOA_COLUMN(X, x, float);                                                    //! cluster hit position in ALICE global coordinate
DECLARE_SOA_COLUMN(Y, y, float);                                                    //! cluster hit position in ALICE global coordinate
DECLARE_SOA_COLUMN(Z, z, float);                                                    //! cluster hit position in ALICE global coordinate
DECLARE_SOA_COLUMN(CellX, cellx, int);                                              //! cell index x of cluster hit position
DECLARE_SOA_COLUMN(CellZ, cellz, int);                                              //! cell index z of cluster hit position
// DECLARE_SOA_COLUMN(TrackEta, tracketa, float);                                      //! eta of the matched track
// DECLARE_SOA_COLUMN(TrackPhi, trackphi, float);                                      //! phi of the matched track
// DECLARE_SOA_COLUMN(TrackP, trackp, float);                                          //! momentum of the matched track
// DECLARE_SOA_COLUMN(TrackPt, trackpt, float);                                        //! pt of the matched track
DECLARE_SOA_DYNAMIC_COLUMN(Px, px, [](float e, float x, float y, float z, float m = 0) -> float { return x / RecoDecay::sqrtSumOfSquares(x, y, z) * sqrt(e * e - m * m); });
DECLARE_SOA_DYNAMIC_COLUMN(Py, py, [](float e, float x, float y, float z, float m = 0) -> float { return y / RecoDecay::sqrtSumOfSquares(x, y, z) * sqrt(e * e - m * m); });
DECLARE_SOA_DYNAMIC_COLUMN(Pz, pz, [](float e, float x, float y, float z, float m = 0) -> float { return z / RecoDecay::sqrtSumOfSquares(x, y, z) * sqrt(e * e - m * m); });
DECLARE_SOA_DYNAMIC_COLUMN(Pt, pt, [](float e, float x, float y, float z, float m = 0) -> float { return RecoDecay::sqrtSumOfSquares(x, y) / RecoDecay::sqrtSumOfSquares(x, y, z) * sqrt(e * e - m * m); });
DECLARE_SOA_DYNAMIC_COLUMN(Eta, eta, [](float x, float y, float z) -> float { return RecoDecay::eta(std::array{x, y, z}); });
DECLARE_SOA_DYNAMIC_COLUMN(Phi, phi, [](float x, float y) -> float { return RecoDecay::phi(x, y); });
} // namespace phoscluster

DECLARE_SOA_TABLE(PHOSClusters, "AOD", "PHOSCLUSTERS", //!
                  o2::soa::Index<>, skimmedcluster::CollisionId, phoscluster::MatchedTrackId,
                  skimmedcluster::E, phoscluster::X, phoscluster::Y, phoscluster::Z,
                  skimmedcluster::M02, skimmedcluster::M20, skimmedcluster::NCells,
                  skimmedcluster::Time, skimmedcluster::DistanceToBadChannel, skimmedcluster::NLM,
                  calocluster::Module, phoscluster::CellX, phoscluster::CellZ,
                  // phoscluster::TrackEta, phoscluster::TrackPhi, phoscluster::TrackP, phoscluster::TrackPt,
                  // dynamic column
                  phoscluster::Px<skimmedcluster::E, phoscluster::X, phoscluster::Y, phoscluster::Z>,
                  phoscluster::Py<skimmedcluster::E, phoscluster::X, phoscluster::Y, phoscluster::Z>,
                  phoscluster::Pz<skimmedcluster::E, phoscluster::X, phoscluster::Y, phoscluster::Z>,
                  phoscluster::Pt<skimmedcluster::E, phoscluster::X, phoscluster::Y, phoscluster::Z>,
                  phoscluster::Eta<phoscluster::X, phoscluster::Y, phoscluster::Z>,
                  phoscluster::Phi<phoscluster::X, phoscluster::Y>);
using PHOSCluster = PHOSClusters::iterator;

DECLARE_SOA_TABLE(PHOSEMEventIds, "AOD", "PHOSEMEVENTID", phoscluster::EMEventId); // To be joined with PHOSClusters table at analysis level.
// iterators
using PHOSEMEventId = PHOSEMEventIds::iterator;

namespace caloextra
{
DECLARE_SOA_INDEX_COLUMN_FULL(Cluster, cluster, int, SkimEMCClusters, ""); //! reference to the gamma in the skimmed EMCal table
DECLARE_SOA_INDEX_COLUMN_FULL(Cell, cell, int, Calos, "");                 //! reference to the gamma in the skimmed EMCal table
// DECLARE_SOA_INDEX_COLUMN(Track, track);                   //! TrackID
DECLARE_SOA_COLUMN(TrackEta, tracketa, float); //! eta of the matched track
DECLARE_SOA_COLUMN(TrackPhi, trackphi, float); //! phi of the matched track
DECLARE_SOA_COLUMN(TrackP, trackp, float);     //! momentum of the matched track
DECLARE_SOA_COLUMN(TrackPt, trackpt, float);   //! pt of the matched track
} // namespace caloextra

DECLARE_SOA_TABLE(SkimEMCCells, "AOD", "SKIMEMCCELLS",                        //! table of link between skimmed EMCal clusters and their cells
                  o2::soa::Index<>, caloextra::ClusterId, caloextra::CellId); //!
using SkimEMCCell = SkimEMCCells::iterator;

DECLARE_SOA_TABLE(SkimEMCMTs, "AOD", "SKIMEMCMTS", //! table of link between skimmed EMCal clusters and their matched tracks
                  o2::soa::Index<>, caloextra::ClusterId, caloextra::TrackEta,
                  caloextra::TrackPhi, caloextra::TrackP, caloextra::TrackPt);
using SkimEMCMT = SkimEMCMTs::iterator;

namespace gammareco
{
DECLARE_SOA_COLUMN(Method, method, int);                                         //! cut bit for PCM photon candidates
DECLARE_SOA_INDEX_COLUMN_FULL(SkimmedPCM, skimmedPCM, int, V0PhotonsKF, "");     //! reference to the gamma in the skimmed PCM table
DECLARE_SOA_INDEX_COLUMN_FULL(SkimmedPHOS, skimmedPHOS, int, PHOSClusters, "");  //! reference to the gamma in the skimmed PHOS table
DECLARE_SOA_INDEX_COLUMN_FULL(SkimmedEMC, skimmedEMC, int, SkimEMCClusters, ""); //! reference to the gamma in the skimmed EMCal table
DECLARE_SOA_COLUMN(PCMCutBit, pcmcutbit, uint64_t);                              //! cut bit for PCM photon candidates
DECLARE_SOA_COLUMN(PHOSCutBit, phoscutbit, uint64_t);                            //! cut bit for PHOS photon candidates
DECLARE_SOA_COLUMN(EMCCutBit, emccutbit, uint64_t);                              //! cut bit for EMCal photon candidates
} // namespace gammareco
DECLARE_SOA_TABLE(SkimGammas, "AOD", "SKIMGAMMAS", //! table of all gamma candidates (PCM, EMCal and PHOS) after cuts
                  o2::soa::Index<>, skimmedcluster::CollisionId, gammareco::Method,
                  skimmedcluster::E, skimmedcluster::Eta, skimmedcluster::Phi,
                  gammareco::SkimmedEMCId, gammareco::SkimmedPHOSId);
DECLARE_SOA_TABLE(SkimPCMCuts, "AOD", "SKIMPCMCUTS",                                  //! table of link between skimmed PCM photon candidates and their cuts
                  o2::soa::Index<>, gammareco::SkimmedPCMId, gammareco::PCMCutBit);   //!
DECLARE_SOA_TABLE(SkimPHOSCuts, "AOD", "SKIMPHOSCUTS",                                //! table of link between skimmed PHOS photon candidates and their cuts
                  o2::soa::Index<>, gammareco::SkimmedPHOSId, gammareco::PHOSCutBit); //!
DECLARE_SOA_TABLE(SkimEMCCuts, "AOD", "SKIMEMCCUTS",                                  //! table of link between skimmed EMCal photon candidates and their cuts
                  o2::soa::Index<>, gammareco::SkimmedEMCId, gammareco::EMCCutBit);   //!
} // namespace o2::aod

#endif // PWGEM_PHOTONMESON_DATAMODEL_GAMMATABLES_H_
