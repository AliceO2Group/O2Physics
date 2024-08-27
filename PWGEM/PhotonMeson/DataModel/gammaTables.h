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
#include "Common/DataModel/CaloClusters.h"

#include "PWGEM/Dilepton/DataModel/dileptonTables.h"

#include "PWGJE/DataModel/EMCALClusters.h"

#ifndef PWGEM_PHOTONMESON_DATAMODEL_GAMMATABLES_H_
#define PWGEM_PHOTONMESON_DATAMODEL_GAMMATABLES_H_

namespace o2::aod
{

namespace emmcbinnedgen
{
DECLARE_SOA_COLUMN(GeneratedGamma, generatedGamma, std::vector<uint16_t>);             //! gamma binned generated data
DECLARE_SOA_COLUMN(GeneratedPi0, generatedPi0, std::vector<uint16_t>);                 //! pi0 binned generated data
DECLARE_SOA_COLUMN(GeneratedEta, generatedEta, std::vector<uint16_t>);                 //! eta binned generated data
DECLARE_SOA_COLUMN(GeneratedOmega, generatedOmega, std::vector<uint16_t>);             //! omega(782) binned generated data
DECLARE_SOA_COLUMN(GeneratedPhi, generatedPhi, std::vector<uint16_t>);                 //! phi(1020) binned generated data
DECLARE_SOA_COLUMN(GeneratedChargedPion, generatedChargedPion, std::vector<uint16_t>); //! charged pion binned generated data
DECLARE_SOA_COLUMN(GeneratedChargedKaon, generatedChargedKaon, std::vector<uint16_t>); //! charged kaon binned generated data
DECLARE_SOA_COLUMN(GeneratedK0S, generatedK0S, std::vector<uint16_t>);                 //! K0S binned generated data
DECLARE_SOA_COLUMN(GeneratedLambda, generatedLambda, std::vector<uint16_t>);           //! Lambda binned generated data

// DECLARE_SOA_COLUMN(GeneratedPi0_Acc_gg, generatedPi0_acc_gg, std::vector<uint16_t>);       //! pi0 -> gg binned generated data
// DECLARE_SOA_COLUMN(GeneratedPi0_Acc_eeg, generatedPi0_acc_eeg, std::vector<uint16_t>);     //! pi0 -> eeg binned generated data
// DECLARE_SOA_COLUMN(GeneratedEta_Acc_gg, generatedEta_acc_gg, std::vector<uint16_t>);       //! eta -> gg binned generated data
// DECLARE_SOA_COLUMN(GeneratedEta_Acc_eeg, generatedEta_acc_eeg, std::vector<uint16_t>);     //! eta -> eeg binned generated data
// DECLARE_SOA_COLUMN(GeneratedEta_Acc_mumug, generatedEta_acc_mumug, std::vector<uint16_t>); //! eta -> mumug binned generated data
// DECLARE_SOA_COLUMN(GeneratedOmega_Acc_ee, generatedOmega_acc_ee, std::vector<uint16_t>);   //! omega(782) -> ee binned generated data
// DECLARE_SOA_COLUMN(GeneratedPhi_Acc_ee, generatedPhi_acc_ee, std::vector<uint16_t>);       //! phi(1020) -> ee binned generated data
// DECLARE_SOA_COLUMN(GeneratedOmega_Acc_mumu, generatedOmega_acc_mumu, std::vector<uint16_t>);   //! omega(782) -> mumu binned generated data
// DECLARE_SOA_COLUMN(GeneratedPhi_Acc_mumu, generatedPhi_acc_mumu, std::vector<uint16_t>);       //! phi(1020) -> mumu binned generated data
} // namespace emmcbinnedgen

DECLARE_SOA_TABLE(BinnedGenPts, "AOD", "BINNEDGENPT", // To be joined with EMMCEvents table at analysis level.
                  emmcbinnedgen::GeneratedGamma,
                  emmcbinnedgen::GeneratedPi0,
                  emmcbinnedgen::GeneratedEta
                  // emmcbinnedgen::GeneratedOmega,
                  // emmcbinnedgen::GeneratedPhi,
                  // emmcbinnedgen::GeneratedChargedPion,
                  // emmcbinnedgen::GeneratedChargedKaon,
                  // emmcbinnedgen::GeneratedK0S,
                  // emmcbinnedgen::GeneratedLambda
);
using BinnedGenPt = BinnedGenPts::iterator;

// DECLARE_SOA_TABLE(BinnedGenPtAccs, "AOD", "BINNEDGENPTACC", // To be joined with EMMCEvents table at analysis level.
//                   emmcbinnedgen::GeneratedPi0_Acc_gg,
//                   emmcbinnedgen::GeneratedPi0_Acc_eeg,
//                   emmcbinnedgen::GeneratedEta_Acc_gg,
//                   emmcbinnedgen::GeneratedEta_Acc_eeg,
//                   emmcbinnedgen::GeneratedEta_Acc_mumug,
//                   emmcbinnedgen::GeneratedOmega_Acc_ee,
//                   emmcbinnedgen::GeneratedPhi_Acc_ee)
////                  emmcbinnedgen::GeneratedOmega_Acc_mumu,
////                  emmcbinnedgen::GeneratedPhi_Acc_mumu);
// using BinnedGenPtAcc = BinnedGenPtAccs::iterator;

namespace v0legmclabel
{
DECLARE_SOA_INDEX_COLUMN(EMMCParticle, emmcparticle); //!
DECLARE_SOA_COLUMN(McMask, mcMask, uint16_t);
} // namespace v0legmclabel

// NOTE: MC labels. This table has one entry for each reconstructed track (joinable with v0leg table)
DECLARE_SOA_TABLE(V0LegMCLabels, "AOD", "V0LEGMCLABEL", //!
                  v0legmclabel::EMMCParticleId, v0legmclabel::McMask);
using V0LegMCLabel = V0LegMCLabels::iterator;

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

namespace emevent
{
DECLARE_SOA_COLUMN(NgPCM, ngpcm, int);
} // namespace emevent

DECLARE_SOA_TABLE(EMEventsNgPCM, "AOD", "EMEVENTNGPCM", emevent::NgPCM); // joinable to EMEvents or aod::Collisions
using EMEventNgPCM = EMEventsNgPCM::iterator;

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

DECLARE_SOA_TABLE(EMPrimaryElectronsFromDalitz, "AOD", "EMPRIMARYELDA", //!
                  o2::soa::Index<>, emprimaryelectron::CollisionId,
                  emprimaryelectron::TrackId, emprimaryelectron::Sign,
                  track::Pt, track::Eta, track::Phi, track::DcaXY, track::DcaZ,
                  track::TPCNClsFindable, track::TPCNClsFindableMinusFound, track::TPCNClsFindableMinusCrossedRows,
                  track::TPCChi2NCl, track::TPCInnerParam,
                  track::TPCSignal, pidtpc::TPCNSigmaEl, pidtpc::TPCNSigmaPi,
                  track::ITSClusterSizes, track::ITSChi2NCl, track::DetectorMap, track::Tgl,

                  // dynamic column
                  track::TPCNClsFound<track::TPCNClsFindable, track::TPCNClsFindableMinusFound>,
                  track::TPCNClsCrossedRows<track::TPCNClsFindable, track::TPCNClsFindableMinusCrossedRows>,
                  track::TPCCrossedRowsOverFindableCls<track::TPCNClsFindable, track::TPCNClsFindableMinusCrossedRows>,
                  track::TPCFoundOverFindableCls<track::TPCNClsFindable, track::TPCNClsFindableMinusFound>,
                  track::v001::ITSClusterMap<track::ITSClusterSizes>, track::v001::ITSNCls<track::ITSClusterSizes>, track::v001::ITSNClsInnerBarrel<track::ITSClusterSizes>,
                  track::HasITS<track::DetectorMap>, track::HasTPC<track::DetectorMap>,
                  track::HasTRD<track::DetectorMap>, track::HasTOF<track::DetectorMap>,
                  emprimaryelectron::Signed1Pt<track::Pt, emprimaryelectron::Sign>,
                  emprimaryelectron::P<track::Pt, track::Eta>,
                  emprimaryelectron::Px<track::Pt, track::Phi>,
                  emprimaryelectron::Py<track::Pt, track::Phi>,
                  emprimaryelectron::Pz<track::Pt, track::Eta>,
                  emprimaryelectron::MeanClusterSizeITS<track::ITSClusterSizes>,
                  emprimaryelectron::MeanClusterSizeITSib<track::ITSClusterSizes>,
                  emprimaryelectron::MeanClusterSizeITSob<track::ITSClusterSizes>);
// iterators
using EMPrimaryElectronFromDalitz = EMPrimaryElectronsFromDalitz::iterator;

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
DECLARE_SOA_INDEX_COLUMN(EMEvent, emevent);                            //!
DECLARE_SOA_COLUMN(CollisionId, collisionId, int);                     //!
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
                  o2::soa::Index<>, skimmedcluster::CollisionId, skimmedcluster::E, skimmedcluster::Eta, skimmedcluster::Phi,
                  skimmedcluster::M02, skimmedcluster::NCells, skimmedcluster::Time, emccluster::IsExotic, emccluster::TrackEta,
                  emccluster::TrackPhi, emccluster::TrackP, emccluster::TrackPt, emccluster::Pt<skimmedcluster::E, skimmedcluster::Eta>);
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
