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

/// \file gammaTables.h
/// \brief This header provides the table definitions for the soft photon and neutral meson tables
/// \author D. Sekihata, daiki.sekihata@cern.ch; Marvin Hemmer (marvin.hemmer@cern.ch) - Goethe University Frankfurt

#include "PWGEM/Dilepton/DataModel/dileptonTables.h"
#include "PWGEM/PhotonMeson/DataModel/EventTables.h" // IWYU pragma: keep (needed)

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/CaloClusters.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CommonConstants/MathConstants.h>
#include <Framework/AnalysisDataModel.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <vector>

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
// 3. EMEMCClusterMCLabels: Vector of EM MC particle ID of largest contributor to cluster
namespace emcclustermclabel
{
DECLARE_SOA_ARRAY_INDEX_COLUMN(McParticle, mcParticle)        //! Array of indice of McParticles
DECLARE_SOA_COLUMN(Amplitude, amplitude, std::vector<float>); //! vector of aplitudes stored in cluster. Ordering is identical to the ordering of EMMCParticleId
} // namespace emcclustermclabel

// NOTE: MC labels. This table has one vector of global mc particle ids for each reconstructed emc cluster (joinable with emccluster table)
DECLARE_SOA_TABLE(EMCClusterMCLabels, "AOD", "EMCClsMCLABEL", //!
                  emcclustermclabel::McParticleIds);

DECLARE_SOA_TABLE_VERSIONED(EMCClusterMCLabels_001, "AOD", "EMCClsMCLABEL", 1, //!
                            emcclustermclabel::McParticleIds, emcclustermclabel::Amplitude);
using EMCClusterMCLabel = EMCClusterMCLabels_001::iterator;

namespace ememcclustermclabel
{
DECLARE_SOA_ARRAY_INDEX_COLUMN(EMMCParticle, emmcparticle);   //! Array of indice of EmMcParticles
DECLARE_SOA_COLUMN(Amplitude, amplitude, std::vector<float>); //! vector of aplitudes stored in cluster. Ordering is identical to the ordering of EMMCParticleId
} // namespace ememcclustermclabel

// NOTE: MC labels. This table has a vector of entries for each reconstructed emc cluster (joinable with emccluster table)
DECLARE_SOA_TABLE(EMEMCClusterMCLabels, "AOD", "EMEMCClsMCLABEL", //!
                  ememcclustermclabel::EMMCParticleIds);

DECLARE_SOA_TABLE_VERSIONED(EMEMCClusterMCLabels_001, "AOD", "EMEMCClsMCLABEL", 1, //!
                            ememcclustermclabel::EMMCParticleIds, ememcclustermclabel::Amplitude);
using EMEMCClusterMCLabel = EMEMCClusterMCLabels_001::iterator;

namespace v0leg
{
DECLARE_SOA_COLUMN(CollisionId, collisionId, int); //!
DECLARE_SOA_COLUMN(TrackId, trackId, int);         //!
DECLARE_SOA_COLUMN(Sign, sign, int8_t);            //!
DECLARE_SOA_COLUMN(IsMoved, isMoved, bool);        //! moved by drift manager. relevant to TPConly tracks
DECLARE_SOA_COLUMN(Px, px, float);                 //! Px at SV
DECLARE_SOA_COLUMN(Py, py, float);                 //! Py at SV
DECLARE_SOA_COLUMN(Pz, pz, float);                 //! Pz at SV

DECLARE_SOA_DYNAMIC_COLUMN(P, p, [](float px, float py, float pz) -> float { return RecoDecay::sqrtSumOfSquares(px, py, pz); });
DECLARE_SOA_DYNAMIC_COLUMN(Pt, pt, [](float px, float py) -> float { return RecoDecay::sqrtSumOfSquares(px, py); });
DECLARE_SOA_DYNAMIC_COLUMN(Eta, eta, [](float px, float py, float pz) -> float { return RecoDecay::eta(std::array{px, py, pz}); });
DECLARE_SOA_DYNAMIC_COLUMN(Phi, phi, [](float px, float py) -> float { return RecoDecay::phi(px, py); });
DECLARE_SOA_DYNAMIC_COLUMN(Tgl, tgl, [](float px, float py, float pz) -> float { return std::tan(o2::constants::math::PIHalf - 2 * std::atan(std::exp(-RecoDecay::eta(std::array{px, py, pz})))); });
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
  }
  return 0;
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
  }
  return 0;
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
  }
  return 0;
});
} // namespace v0leg
DECLARE_SOA_TABLE(V0Legs_000, "AOD", "V0LEG", //!
                  o2::soa::Index<>, v0leg::CollisionId, v0leg::TrackId, v0leg::Sign,
                  v0leg::Px, v0leg::Py, v0leg::Pz,
                  track::DcaXY, track::DcaZ,
                  track::TPCNClsFindable, track::TPCNClsFindableMinusFound, track::TPCNClsFindableMinusCrossedRows, track::TPCNClsShared,
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
                  track::TPCFractionSharedCls<track::TPCNClsShared, track::TPCNClsFindable, track::TPCNClsFindableMinusFound>,
                  track::v001::ITSClusterMap<track::ITSClusterSizes>, track::v001::ITSNCls<track::ITSClusterSizes>, track::v001::ITSNClsInnerBarrel<track::ITSClusterSizes>,
                  track::HasITS<track::DetectorMap>, track::HasTPC<track::DetectorMap>, track::HasTRD<track::DetectorMap>, track::HasTOF<track::DetectorMap>,
                  v0leg::MeanClusterSizeITS<track::ITSClusterSizes>,
                  v0leg::MeanClusterSizeITSib<track::ITSClusterSizes>,
                  v0leg::MeanClusterSizeITSob<track::ITSClusterSizes>);

DECLARE_SOA_TABLE_VERSIONED(V0Legs_001, "AOD", "V0LEG", 1, //!
                            o2::soa::Index<>, v0leg::CollisionId, v0leg::TrackId, v0leg::Sign,
                            v0leg::Px, v0leg::Py, v0leg::Pz,
                            track::DcaXY, track::DcaZ,
                            track::TPCNClsFindable, track::TPCNClsFindableMinusFound, track::TPCNClsFindableMinusCrossedRows, track::TPCNClsShared,
                            track::TPCChi2NCl, track::TPCInnerParam,
                            track::TPCSignal, pidtpc::TPCNSigmaEl, pidtpc::TPCNSigmaPi,
                            track::ITSClusterSizes, track::ITSChi2NCl, track::DetectorMap,

                            // dynamic column
                            v0leg::P<v0leg::Px, v0leg::Py, v0leg::Pz>,
                            v0leg::Pt<v0leg::Px, v0leg::Py>,
                            v0leg::Eta<v0leg::Px, v0leg::Py, v0leg::Pz>,
                            v0leg::Phi<v0leg::Px, v0leg::Py>,
                            v0leg::Tgl<v0leg::Px, v0leg::Py, v0leg::Pz>,

                            track::TPCNClsFound<track::TPCNClsFindable, track::TPCNClsFindableMinusFound>,
                            track::TPCNClsCrossedRows<track::TPCNClsFindable, track::TPCNClsFindableMinusCrossedRows>,
                            track::TPCCrossedRowsOverFindableCls<track::TPCNClsFindable, track::TPCNClsFindableMinusCrossedRows>,
                            track::TPCFoundOverFindableCls<track::TPCNClsFindable, track::TPCNClsFindableMinusFound>,
                            track::TPCFractionSharedCls<track::TPCNClsShared, track::TPCNClsFindable, track::TPCNClsFindableMinusFound>,
                            track::v001::ITSClusterMap<track::ITSClusterSizes>, track::v001::ITSNCls<track::ITSClusterSizes>, track::v001::ITSNClsInnerBarrel<track::ITSClusterSizes>,
                            track::HasITS<track::DetectorMap>, track::HasTPC<track::DetectorMap>, track::HasTRD<track::DetectorMap>, track::HasTOF<track::DetectorMap>,
                            v0leg::MeanClusterSizeITS<track::ITSClusterSizes>,
                            v0leg::MeanClusterSizeITSib<track::ITSClusterSizes>,
                            v0leg::MeanClusterSizeITSob<track::ITSClusterSizes>);

using V0Legs = V0Legs_001;
// iterators
using V0Leg = V0Legs::iterator;

// copy of V0Legs table
DECLARE_SOA_TABLE(V0LegsTmp, "AOD", "V0LEGTMP", //!
                  o2::soa::Index<>, v0leg::CollisionId, v0leg::TrackId, v0leg::Sign,
                  v0leg::Px, v0leg::Py, v0leg::Pz,
                  track::DcaXY, track::DcaZ,
                  track::TPCNClsFindable, track::TPCNClsFindableMinusFound, track::TPCNClsFindableMinusCrossedRows, track::TPCNClsShared,
                  track::TPCChi2NCl, track::TPCInnerParam,
                  track::TPCSignal, pidtpc::TPCNSigmaEl, pidtpc::TPCNSigmaPi,
                  track::ITSClusterSizes, track::ITSChi2NCl, track::DetectorMap,

                  // dynamic column
                  v0leg::P<v0leg::Px, v0leg::Py, v0leg::Pz>,
                  v0leg::Pt<v0leg::Px, v0leg::Py>,
                  v0leg::Eta<v0leg::Px, v0leg::Py, v0leg::Pz>,
                  v0leg::Phi<v0leg::Px, v0leg::Py>,
                  v0leg::Tgl<v0leg::Px, v0leg::Py, v0leg::Pz>,

                  track::TPCNClsFound<track::TPCNClsFindable, track::TPCNClsFindableMinusFound>,
                  track::TPCNClsCrossedRows<track::TPCNClsFindable, track::TPCNClsFindableMinusCrossedRows>,
                  track::TPCCrossedRowsOverFindableCls<track::TPCNClsFindable, track::TPCNClsFindableMinusCrossedRows>,
                  track::TPCFoundOverFindableCls<track::TPCNClsFindable, track::TPCNClsFindableMinusFound>,
                  track::TPCFractionSharedCls<track::TPCNClsShared, track::TPCNClsFindable, track::TPCNClsFindableMinusFound>,
                  track::v001::ITSClusterMap<track::ITSClusterSizes>, track::v001::ITSNCls<track::ITSClusterSizes>, track::v001::ITSNClsInnerBarrel<track::ITSClusterSizes>,
                  track::HasITS<track::DetectorMap>, track::HasTPC<track::DetectorMap>, track::HasTRD<track::DetectorMap>, track::HasTOF<track::DetectorMap>,
                  v0leg::MeanClusterSizeITS<track::ITSClusterSizes>,
                  v0leg::MeanClusterSizeITSib<track::ITSClusterSizes>,
                  v0leg::MeanClusterSizeITSob<track::ITSClusterSizes>);
// iterators
using V0LegTmp = V0LegsTmp::iterator;

DECLARE_SOA_TABLE_VERSIONED(V0LegsXYZ_000, "AOD", "V0LEGXYZ", 0, track::X, track::Y, track::Z);
using V0LegsXYZ = V0LegsXYZ_000;
// iterators
using V0LegXYZ = V0LegsXYZ::iterator;

// copy of V0Legs table
DECLARE_SOA_TABLE(V0LegsXYZTmp, "AOD", "V0LEGXYZTMP", track::X, track::Y, track::Z);
// iterators
using V0LegXYZTmp = V0LegsXYZTmp::iterator;

DECLARE_SOA_TABLE_VERSIONED(V0LegsDeDxMC_000, "AOD", "V0LEGDEDXMC", 0, mcpidtpc::DeDxTunedMc, o2::soa::Marker<2>);
using V0LegsDeDxMC = V0LegsDeDxMC_000;
// iterators
using V0LegDeDxMC = V0LegsDeDxMC::iterator;

// copy of V0Legs table
DECLARE_SOA_TABLE_VERSIONED(V0LegsDeDxMCTmp, "AOD", "V0LEGDEDXMCTMP", 0, mcpidtpc::DeDxTunedMc, o2::soa::Marker<2>);
// iterators
using V0LegDeDxMCTmp = V0LegsDeDxMCTmp::iterator;

namespace v0legcov
{
DECLARE_SOA_COLUMN(X, x, float);                 //! x of the leg parameters given to the V0 fit (innermost update, cm)
DECLARE_SOA_COLUMN(Alpha, alpha, float);         //! frame angle
DECLARE_SOA_COLUMN(Y, y, float);                 //! local y (cm)
DECLARE_SOA_COLUMN(Z, z, float);                 //! z (cm), TPC-only legs after the drift correction to the V0's collision
DECLARE_SOA_COLUMN(Snp, snp, float);             //! sin(phi) in the local frame
DECLARE_SOA_COLUMN(Tgl, tgl, float);             //! tan(lambda)
DECLARE_SOA_COLUMN(Signed1Pt, signed1Pt, float); //! q/pT (c/GeV)
DECLARE_SOA_COLUMN(Cov, cov, float[15]);         // NOLINT(modernize-avoid-c-arrays) //! covariance matrix, o2::track::TrackParCov order
} // namespace v0legcov

// leg track parameters as used in the V0 fit, one row per V0Legs row (optional builder output, joinable with V0Legs)
DECLARE_SOA_TABLE_VERSIONED(V0LegsCov_000, "AOD", "V0LEGCOV", 0, //!
                            v0legcov::X, v0legcov::Alpha, v0legcov::Y, v0legcov::Z, v0legcov::Snp, v0legcov::Tgl, v0legcov::Signed1Pt, v0legcov::Cov);
using V0LegsCov = V0LegsCov_000;
// iterators
using V0LegCov = V0LegsCov::iterator;

// copy of V0LegsCov table
DECLARE_SOA_TABLE_VERSIONED(V0LegsCovTmp, "AOD", "V0LEGCOVTMP", 0, //!
                            v0legcov::X, v0legcov::Alpha, v0legcov::Y, v0legcov::Z, v0legcov::Snp, v0legcov::Tgl, v0legcov::Signed1Pt, v0legcov::Cov);
// iterators
using V0LegCovTmp = V0LegsCovTmp::iterator;

namespace emevent
{
DECLARE_SOA_COLUMN(NgPCM, ngpcm, int);
DECLARE_SOA_COLUMN(Weight, weight, float); //! Weight of the event (e.g. for JJ MCs). Set to 1 for data and non-weighted MCs.
} // namespace emevent

DECLARE_SOA_TABLE(EMEventsNgPCM, "AOD", "EMEVENTNGPCM", emevent::NgPCM); // joinable to EMEvents or o2::aod::Collisions
using EMEventNgPCM = EMEventsNgPCM::iterator;

DECLARE_SOA_TABLE(EMEventsWeight, "AOD", "EMEVENTWEIGHT", //! table contanint the weight for eache event (for JJ MCs), joinable to EMEvents
                  emevent::Weight);
using EMEventWeight = EMEventsWeight::iterator;

namespace v0photonkf
{
DECLARE_SOA_INDEX_COLUMN(EMEvent, emevent);                             //!
DECLARE_SOA_INDEX_COLUMN(PMEvent, pmevent);                             //!
DECLARE_SOA_COLUMN(CollisionId, collisionId, int);                      //!
DECLARE_SOA_COLUMN(V0Id, v0Id, int);                                    //!
DECLARE_SOA_INDEX_COLUMN_FULL(PosTrack, posTrack, int, V0Legs, "_Pos"); //!
DECLARE_SOA_INDEX_COLUMN_FULL(NegTrack, negTrack, int, V0Legs, "_Neg"); //!
DECLARE_SOA_COLUMN(Vx, vx, float);                                      //! secondary vertex x
DECLARE_SOA_COLUMN(Vy, vy, float);                                      //! secondary vertex y
DECLARE_SOA_COLUMN(Vz, vz, float);                                      //! secondary vertex z
DECLARE_SOA_COLUMN(Px, px, float);                                      //! px for photon kf
DECLARE_SOA_COLUMN(Py, py, float);                                      //! py for photon kf
DECLARE_SOA_COLUMN(Pz, pz, float);                                      //! pz for photon kf

DECLARE_SOA_COLUMN(MGamma, mGamma, float);             //! invariant mass of dielectron at SV
DECLARE_SOA_COLUMN(DCAxyToPV, dcaXYtopv, float);       //! DCAxy of V0 to PV
DECLARE_SOA_COLUMN(DCAzToPV, dcaZtopv, float);         //! DCAz of V0 to PV
DECLARE_SOA_COLUMN(CosPA, cospa, float);               //!
DECLARE_SOA_COLUMN(CosPAXY, cospaXY, float);           //!
DECLARE_SOA_COLUMN(CosPARZ, cospaRZ, float);           //!
DECLARE_SOA_COLUMN(PCA, pca, float);                   //!
DECLARE_SOA_COLUMN(Alpha, alpha, float);               //!
DECLARE_SOA_COLUMN(QtArm, qtarm, float);               //!
DECLARE_SOA_COLUMN(ChiSquareNDF, chiSquareNDF, float); //! Chi2 / NDF of the reconstructed V0

DECLARE_SOA_COLUMN(SigmaPx2, sigmaPx2, float);                 //! error^2 of px in covariant matrix
DECLARE_SOA_COLUMN(SigmaPy2, sigmaPy2, float);                 //! error^2 of py in covariant matrix
DECLARE_SOA_COLUMN(SigmaPz2, sigmaPz2, float);                 //! error^2 of pz in covariant matrix
DECLARE_SOA_COLUMN(SigmaPxPy, sigmaPxPy, float);               //! error of px x py in covariant matrix
DECLARE_SOA_COLUMN(SigmaPyPz, sigmaPyPz, float);               //! error of py x pz in covariant matrix
DECLARE_SOA_COLUMN(SigmaPzPx, sigmaPzPx, float);               //! error of pz x px in covariant matrix
DECLARE_SOA_COLUMN(PrefilterBitDerived, pfbderived, uint16_t); //!

DECLARE_SOA_DYNAMIC_COLUMN(E, e, [](float px, float py, float pz, float m = 0) -> float { return RecoDecay::sqrtSumOfSquares(px, py, pz, m); }); //! energy of v0 photn, mass to be given as argument when getter is called!
DECLARE_SOA_DYNAMIC_COLUMN(Pt, pt, [](float px, float py) -> float { return RecoDecay::sqrtSumOfSquares(px, py); });
DECLARE_SOA_DYNAMIC_COLUMN(Eta, eta, [](float px, float py, float pz) -> float { return RecoDecay::eta(std::array{px, py, pz}); });
DECLARE_SOA_DYNAMIC_COLUMN(Phi, phi, [](float px, float py) -> float { return RecoDecay::phi(px, py); });
DECLARE_SOA_DYNAMIC_COLUMN(P, p, [](float px, float py, float pz) -> float { return RecoDecay::sqrtSumOfSquares(px, py, pz); });
DECLARE_SOA_DYNAMIC_COLUMN(V0Radius, v0radius, [](float vx, float vy) -> float { return RecoDecay::sqrtSumOfSquares(vx, vy); });
} // namespace v0photonkf
DECLARE_SOA_TABLE(V0PhotonsKF_000, "AOD", "V0PHOTONKF", //!
                  o2::soa::Index<>, v0photonkf::CollisionId, v0photonkf::V0Id, v0photonkf::PosTrackId, v0photonkf::NegTrackId,
                  v0photonkf::Vx, v0photonkf::Vy, v0photonkf::Vz,
                  v0photonkf::Px, v0photonkf::Py, v0photonkf::Pz,
                  v0photonkf::MGamma,
                  v0photonkf::DCAxyToPV, v0photonkf::DCAzToPV,
                  v0photonkf::CosPA, v0photonkf::CosPAXY, v0photonkf::CosPARZ, v0photonkf::PCA,
                  v0photonkf::Alpha, v0photonkf::QtArm,
                  v0photonkf::ChiSquareNDF,

                  // dynamic column
                  v0photonkf::E<v0photonkf::Px, v0photonkf::Py, v0photonkf::Pz>,
                  v0photonkf::Pt<v0photonkf::Px, v0photonkf::Py>,
                  v0photonkf::Eta<v0photonkf::Px, v0photonkf::Py, v0photonkf::Pz>,
                  v0photonkf::Phi<v0photonkf::Px, v0photonkf::Py>,
                  v0photonkf::P<v0photonkf::Px, v0photonkf::Py, v0photonkf::Pz>,
                  v0photonkf::V0Radius<v0photonkf::Vx, v0photonkf::Vy>);

DECLARE_SOA_TABLE_VERSIONED(V0PhotonsKF_001, "AOD", "V0PHOTONKF", 1, //!
                            o2::soa::Index<>, v0photonkf::CollisionId, v0photonkf::V0Id, v0photonkf::PosTrackId, v0photonkf::NegTrackId,
                            v0photonkf::Vx, v0photonkf::Vy, v0photonkf::Vz,
                            v0photonkf::Px, v0photonkf::Py, v0photonkf::Pz,
                            v0photonkf::MGamma,
                            v0photonkf::DCAxyToPV, v0photonkf::DCAzToPV,
                            v0photonkf::CosPA, v0photonkf::CosPAXY, v0photonkf::CosPARZ, v0photonkf::PCA,
                            v0photonkf::Alpha, v0photonkf::QtArm,
                            v0photonkf::ChiSquareNDF,

                            // dynamic column
                            v0photonkf::E<v0photonkf::Px, v0photonkf::Py, v0photonkf::Pz>,
                            v0photonkf::Pt<v0photonkf::Px, v0photonkf::Py>,
                            v0photonkf::Eta<v0photonkf::Px, v0photonkf::Py, v0photonkf::Pz>,
                            v0photonkf::Phi<v0photonkf::Px, v0photonkf::Py>,
                            v0photonkf::P<v0photonkf::Px, v0photonkf::Py, v0photonkf::Pz>,
                            v0photonkf::V0Radius<v0photonkf::Vx, v0photonkf::Vy>);

using V0PhotonsKF = V0PhotonsKF_001;

// iterators
using V0PhotonKF = V0PhotonsKF::iterator;

// Clone of V0PhotonKF table for intermediate processing needed for dalitz
DECLARE_SOA_TABLE(V0PhotonsKFTmp, "AOD", "V0PHOTONKFTMP", //!
                  o2::soa::Index<>, v0photonkf::CollisionId, v0photonkf::V0Id, v0photonkf::PosTrackId, v0photonkf::NegTrackId,
                  v0photonkf::Vx, v0photonkf::Vy, v0photonkf::Vz,
                  v0photonkf::Px, v0photonkf::Py, v0photonkf::Pz,
                  v0photonkf::MGamma,
                  v0photonkf::DCAxyToPV, v0photonkf::DCAzToPV,
                  v0photonkf::CosPA, v0photonkf::CosPAXY, v0photonkf::CosPARZ, v0photonkf::PCA,
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
using V0PhotonKFTmp = V0PhotonsKFTmp::iterator;

DECLARE_SOA_TABLE(V0KFEMEventIds_000, "AOD", "V0KFEMEVENTID", v0photonkf::EMEventId);              // To be joined with V0PhotonsKF table at analysis level.
DECLARE_SOA_TABLE_VERSIONED(V0KFEMEventIds_001, "AOD", "V0KFEMEVENTID", 1, v0photonkf::PMEventId); // To be joined with V0PhotonsKF table at analysis level.
using V0KFEMEventIds = V0KFEMEventIds_001;
using V0KFEMEventId = V0KFEMEventIds::iterator;

DECLARE_SOA_TABLE(V0PhotonsKFCov, "AOD", "V0PHOTONKFCOV", //! To be joined with V0PhotonsKF table at analysis level.
                  v0photonkf::SigmaPx2, v0photonkf::SigmaPy2, v0photonkf::SigmaPz2, v0photonkf::SigmaPxPy, v0photonkf::SigmaPyPz, v0photonkf::SigmaPzPx);
// iterators
using V0PhotonKFCov = V0PhotonsKFCov::iterator;

DECLARE_SOA_TABLE(V0PhotonsKFPrefilterBitDerived, "AOD", "V0PHOTONKFPFBPI0", v0photonkf::PrefilterBitDerived); // To be joined with V0PhotonsKF table at analysis level.
// iterators
using V0PhotonKFPrefilterBitDerived = V0PhotonsKFPrefilterBitDerived::iterator;

DECLARE_SOA_TABLE(EMPrimaryElectronsFromDalitz_000, "AOD", "EMPRIMARYELDA", //!
                  o2::soa::Index<>, emprimaryelectron::CollisionId,
                  emprimaryelectron::TrackId, emprimaryelectron::Sign,
                  track::Pt, track::Eta, track::Phi, track::DcaXY, track::DcaZ, track::CYY, track::CZY, track::CZZ,
                  track::TPCNClsFindable, track::TPCNClsFindableMinusFound, track::TPCNClsFindableMinusCrossedRows, track::TPCNClsShared,
                  track::TPCChi2NCl, track::TPCInnerParam,
                  track::TPCSignal, pidtpc::TPCNSigmaEl, pidtpc::TPCNSigmaPi,
                  pidtofbeta::Beta, pidtof::TOFNSigmaEl, pidtof::TOFNSigmaPi,
                  track::ITSClusterSizes, track::ITSChi2NCl, track::TOFChi2, track::DetectorMap, track::Tgl,

                  // dynamic column
                  track::TPCNClsFound<track::TPCNClsFindable, track::TPCNClsFindableMinusFound>,
                  track::TPCNClsCrossedRows<track::TPCNClsFindable, track::TPCNClsFindableMinusCrossedRows>,
                  track::TPCCrossedRowsOverFindableCls<track::TPCNClsFindable, track::TPCNClsFindableMinusCrossedRows>,
                  track::TPCFoundOverFindableCls<track::TPCNClsFindable, track::TPCNClsFindableMinusFound>,
                  track::v001::ITSClusterMap<track::ITSClusterSizes>, track::v001::ITSNCls<track::ITSClusterSizes>, track::v001::ITSNClsInnerBarrel<track::ITSClusterSizes>,
                  track::TPCFractionSharedCls<track::TPCNClsShared, track::TPCNClsFindable, track::TPCNClsFindableMinusFound>,
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

DECLARE_SOA_TABLE_VERSIONED(EMPrimaryElectronsFromDalitz_001, "AOD", "EMPRIMARYELDA", 1, //!
                            o2::soa::Index<>, emprimaryelectron::CollisionId,
                            emprimaryelectron::TrackId, emprimaryelectron::Sign,
                            track::Pt, track::Eta, track::Phi, track::DcaXY, track::DcaZ, track::CYY, track::CZY, track::CZZ,
                            track::TPCNClsFindable, track::TPCNClsFindableMinusFound, track::TPCNClsFindableMinusCrossedRows, track::TPCNClsShared,
                            track::TPCChi2NCl, track::TPCInnerParam,
                            track::TPCSignal, pidtpc::TPCNSigmaEl, pidtpc::TPCNSigmaPi,
                            pidtofbeta::Beta, pidtof::TOFNSigmaEl,
                            track::ITSClusterSizes, track::ITSChi2NCl, track::TOFChi2, track::DetectorMap,

                            // dynamic column
                            track::TPCNClsFound<track::TPCNClsFindable, track::TPCNClsFindableMinusFound>,
                            track::TPCNClsCrossedRows<track::TPCNClsFindable, track::TPCNClsFindableMinusCrossedRows>,
                            track::TPCCrossedRowsOverFindableCls<track::TPCNClsFindable, track::TPCNClsFindableMinusCrossedRows>,
                            track::TPCFoundOverFindableCls<track::TPCNClsFindable, track::TPCNClsFindableMinusFound>,
                            track::v001::ITSClusterMap<track::ITSClusterSizes>, track::v001::ITSNCls<track::ITSClusterSizes>, track::v001::ITSNClsInnerBarrel<track::ITSClusterSizes>,
                            track::TPCFractionSharedCls<track::TPCNClsShared, track::TPCNClsFindable, track::TPCNClsFindableMinusFound>,
                            track::HasITS<track::DetectorMap>, track::HasTPC<track::DetectorMap>, track::HasTRD<track::DetectorMap>, track::HasTOF<track::DetectorMap>,

                            emprimaryelectron::Signed1Pt<track::Pt, emprimaryelectron::Sign>,
                            emprimaryelectron::P<track::Pt, track::Eta>,
                            emprimaryelectron::Px<track::Pt, track::Phi>,
                            emprimaryelectron::Py<track::Pt, track::Phi>,
                            emprimaryelectron::Pz<track::Pt, track::Eta>,
                            emprimaryelectron::Tgl<track::Eta>,
                            emprimaryelectron::MeanClusterSizeITS<track::ITSClusterSizes>,
                            emprimaryelectron::MeanClusterSizeITSib<track::ITSClusterSizes>,
                            emprimaryelectron::MeanClusterSizeITSob<track::ITSClusterSizes>);

using EMPrimaryElectronsFromDalitz = EMPrimaryElectronsFromDalitz_001;

// iterators
using EMPrimaryElectronFromDalitz = EMPrimaryElectronsFromDalitz::iterator;

// Table that is a clone of EMPrimaryELectronsFromDalitz but needed for event selection purposes
DECLARE_SOA_TABLE(EMPrimaryElectronsFromDalitzTmp, "AOD", "EMPRIMARYELDATMP", //!
                  o2::soa::Index<>, emprimaryelectron::CollisionId,
                  emprimaryelectron::TrackId, emprimaryelectron::Sign,
                  track::Pt, track::Eta, track::Phi, track::DcaXY, track::DcaZ, track::CYY, track::CZY, track::CZZ,
                  track::TPCNClsFindable, track::TPCNClsFindableMinusFound, track::TPCNClsFindableMinusCrossedRows, track::TPCNClsShared,
                  track::TPCChi2NCl, track::TPCInnerParam,
                  track::TPCSignal, pidtpc::TPCNSigmaEl, pidtpc::TPCNSigmaPi,
                  pidtofbeta::Beta, pidtof::TOFNSigmaEl,
                  track::ITSClusterSizes, track::ITSChi2NCl, track::TOFChi2, track::DetectorMap,

                  // dynamic column
                  track::TPCNClsFound<track::TPCNClsFindable, track::TPCNClsFindableMinusFound>,
                  track::TPCNClsCrossedRows<track::TPCNClsFindable, track::TPCNClsFindableMinusCrossedRows>,
                  track::TPCCrossedRowsOverFindableCls<track::TPCNClsFindable, track::TPCNClsFindableMinusCrossedRows>,
                  track::TPCFoundOverFindableCls<track::TPCNClsFindable, track::TPCNClsFindableMinusFound>,
                  track::v001::ITSClusterMap<track::ITSClusterSizes>, track::v001::ITSNCls<track::ITSClusterSizes>, track::v001::ITSNClsInnerBarrel<track::ITSClusterSizes>,
                  track::TPCFractionSharedCls<track::TPCNClsShared, track::TPCNClsFindable, track::TPCNClsFindableMinusFound>,
                  track::HasITS<track::DetectorMap>, track::HasTPC<track::DetectorMap>, track::HasTRD<track::DetectorMap>, track::HasTOF<track::DetectorMap>,

                  emprimaryelectron::Signed1Pt<track::Pt, emprimaryelectron::Sign>,
                  emprimaryelectron::P<track::Pt, track::Eta>,
                  emprimaryelectron::Px<track::Pt, track::Phi>,
                  emprimaryelectron::Py<track::Pt, track::Phi>,
                  emprimaryelectron::Pz<track::Pt, track::Eta>,
                  emprimaryelectron::Tgl<track::Eta>,
                  emprimaryelectron::MeanClusterSizeITS<track::ITSClusterSizes>,
                  emprimaryelectron::MeanClusterSizeITSib<track::ITSClusterSizes>,
                  emprimaryelectron::MeanClusterSizeITSob<track::ITSClusterSizes>);

// iterators
using EMPrimaryElectronFromDalitzTmp = EMPrimaryElectronsFromDalitzTmp::iterator;

// copy of EMPrimaryElectronDeDxMC to use as intermediate step for Dalitz analysis
DECLARE_SOA_TABLE(EMPrimaryElectronsDeDxMCTmp, "AOD", "EMELDEDXMCTMP", mcpidtpc::DeDxTunedMc, o2::soa::Marker<1>);
using EMPrimaryElectronsDeDxMCTmp = EMPrimaryElectronsDeDxMCTmp;
// iterators
using EMPrimaryElectronDeDxMCTmp = EMPrimaryElectronsDeDxMCTmp::iterator;

namespace emprimaryelectronda
{
DECLARE_SOA_INDEX_COLUMN(EMEvent, emevent); //!
DECLARE_SOA_INDEX_COLUMN(PMEvent, pmevent); //!
} // namespace emprimaryelectronda

DECLARE_SOA_TABLE_VERSIONED(EMPrimaryElectronDaEMEventIds_001, "AOD", "PRMELDAEVENTID", 1, emprimaryelectronda::PMEventId); // To be joined with EMPrimaryElectronsFromDalitz table at analysis level.
using EMPrimaryElectronDaEMEventIds = EMPrimaryElectronDaEMEventIds_001;
using EMPrimaryElectronDaEMEventId = EMPrimaryElectronDaEMEventIds::iterator;

namespace v0photonsphivpsi
{
DECLARE_SOA_INDEX_COLUMN(EMEvent, emevent); //!
DECLARE_SOA_INDEX_COLUMN(PMEvent, pmevent); //!
DECLARE_SOA_COLUMN(PhiV, phiv, float);      //!
DECLARE_SOA_COLUMN(PsiPair, psipair, float);
} // namespace v0photonsphivpsi
// KF quality of the photon: the stored v0photonkf::ChiSquareNDF is the chi2/ndf AFTER the topological constraint to the PV;
// these two separate the conversion-vertex fit from the pointing to the PV (row-by-row extension of the photon table, join by index)
namespace v0photonkfchi2
{
DECLARE_SOA_COLUMN(Chi2NDFVtx, chi2NDFVtx, float); //! chi2/ndf of the two-leg conversion-vertex fit (incl. the mass constraint if applied), before the PV constraint
DECLARE_SOA_COLUMN(Chi2PV, chi2PV, float);         //! chi2 deviation of the photon from the primary vertex (KFParticle::GetDeviationFromVertex)
} // namespace v0photonkfchi2
DECLARE_SOA_TABLE(V0PhotonsKFChi2, "AOD", "V0PHKFCHI2", o2::soa::Index<>, v0photonkfchi2::Chi2NDFVtx, v0photonkfchi2::Chi2PV);       //!
DECLARE_SOA_TABLE(V0PhotonsKFChi2Tmp, "AOD", "V0PHKFCHI2TMP", o2::soa::Index<>, v0photonkfchi2::Chi2NDFVtx, v0photonkfchi2::Chi2PV); //!
DECLARE_SOA_TABLE(V0PhotonsKFChi2Mix, "AOD", "V0PHKFCHI2MIX", o2::soa::Index<>, v0photonkfchi2::Chi2NDFVtx, v0photonkfchi2::Chi2PV); //! builder run 2

DECLARE_SOA_TABLE(V0PhotonsPhiVPsi, "AOD", "V0PHOTONPHIVPSI", //!
                  o2::soa::Index<>, v0photonsphivpsi::PhiV, v0photonsphivpsi::PsiPair);
// iterators
using V0PhotonsPhiVPsi = V0PhotonsPhiVPsi;

// clone for temporary table
DECLARE_SOA_TABLE(V0PhotonsPhiVPsiTmp, "AOD", "V0PHOTONPHIVTMP", //!
                  o2::soa::Index<>, v0photonsphivpsi::PhiV, v0photonsphivpsi::PsiPair);
// iterators
using V0PhotonsPhiVPsiTmp = V0PhotonsPhiVPsiTmp;

// co-reconstructed mixed event: artificial events made by photon-mixed-event-producer from the photons of builder run 1,
// rebuilt from their legs by builder run 2 (photon-conversion-builder-mixed)
namespace pcmmixcol
{
DECLARE_SOA_COLUMN(IsMixed, isMixed, bool);                    //! false: the photons of one event, true: N photons from N different events
DECLARE_SOA_COLUMN(SourceCollisionId, sourceCollisionId, int); //! collision (this DF) whose vertex the artificial event has
DECLARE_SOA_COLUMN(SourceEventUid, sourceEventUid, int64_t);   //! mixer event counter of that collision
DECLARE_SOA_COLUMN(ZBin, zBin, int);                           //! pool z bin
DECLARE_SOA_COLUMN(OccBin, occBin, int);                       //! pool occupancy bin
DECLARE_SOA_COLUMN(Occupancy, occupancy, float);               //! occupancy of the source collision
DECLARE_SOA_COLUMN(CentFT0M, centFT0M, float);                 //! centralities of the source collision (V0 fit input)
DECLARE_SOA_COLUMN(CentFT0A, centFT0A, float);                 //!
DECLARE_SOA_COLUMN(CentFT0C, centFT0C, float);                 //!
DECLARE_SOA_COLUMN(FirstLeg, firstLeg, int);                   //! first PCMMixLegs row of this artificial event
DECLARE_SOA_COLUMN(NLegs, nLegs, int);                         //! number of PCMMixLegs rows (2 per photon, contiguous)
DECLARE_SOA_COLUMN(RunNumber, runNumber, int);                 //! run of the source collision (CCDB: B field, TPC drift velocity)
DECLARE_SOA_COLUMN(Timestamp, timestamp, uint64_t);            //! CCDB timestamp of the source collision (builder run 2 on derived data: no BCs)
DECLARE_SOA_COLUMN(MultNTracksPVeta1, multNTracksPVeta1, int); //! PV contributors with |eta| < 1 of the source collision
DECLARE_SOA_COLUMN(MultFT0M, multFT0M, float);                 //! FT0A + FT0C amplitude of the source collision
} // namespace pcmmixcol
DECLARE_SOA_TABLE(PCMMixCollisions, "AOD", "PCMMIXCOL", //!
                  o2::soa::Index<>, collision::PosX, collision::PosY, collision::PosZ,
                  collision::CovXX, collision::CovXY, collision::CovXZ, collision::CovYY, collision::CovYZ, collision::CovZZ,
                  collision::Chi2, collision::NumContrib,
                  pcmmixcol::CentFT0M, pcmmixcol::CentFT0A, pcmmixcol::CentFT0C,
                  pcmmixcol::IsMixed, pcmmixcol::SourceCollisionId, pcmmixcol::SourceEventUid,
                  pcmmixcol::ZBin, pcmmixcol::OccBin, pcmmixcol::Occupancy, pcmmixcol::FirstLeg, pcmmixcol::NLegs,
                  pcmmixcol::RunNumber, pcmmixcol::MultNTracksPVeta1, pcmmixcol::MultFT0M, pcmmixcol::Timestamp);
using PCMMixCollision = PCMMixCollisions::iterator;

namespace pcmmixleg
{
DECLARE_SOA_COLUMN(CollisionId, collisionId, int);                //! PCMMixCollisions row
DECLARE_SOA_COLUMN(Pt, pt, float);                                //! pT of the AO2D track (track-like interface for the builder)
DECLARE_SOA_COLUMN(Phi, phi, float);                              //! phi of the AO2D track
DECLARE_SOA_COLUMN(CYY, cYY, float);                              //! covariance, o2::track::TrackParCov order (from V0LegsCov)
DECLARE_SOA_COLUMN(CZY, cZY, float);                              //!
DECLARE_SOA_COLUMN(CZZ, cZZ, float);                              //!
DECLARE_SOA_COLUMN(CSnpY, cSnpY, float);                          //!
DECLARE_SOA_COLUMN(CSnpZ, cSnpZ, float);                          //!
DECLARE_SOA_COLUMN(CSnpSnp, cSnpSnp, float);                      //!
DECLARE_SOA_COLUMN(CTglY, cTglY, float);                          //!
DECLARE_SOA_COLUMN(CTglZ, cTglZ, float);                          //!
DECLARE_SOA_COLUMN(CTglSnp, cTglSnp, float);                      //!
DECLARE_SOA_COLUMN(CTglTgl, cTglTgl, float);                      //!
DECLARE_SOA_COLUMN(C1PtY, c1PtY, float);                          //!
DECLARE_SOA_COLUMN(C1PtZ, c1PtZ, float);                          //!
DECLARE_SOA_COLUMN(C1PtSnp, c1PtSnp, float);                      //!
DECLARE_SOA_COLUMN(C1PtTgl, c1PtTgl, float);                      //!
DECLARE_SOA_COLUMN(C1Pt21Pt2, c1Pt21Pt2, float);                  //!
DECLARE_SOA_COLUMN(SourceEventUid, sourceEventUid, int64_t);      //! mixer event counter of the event the leg comes from
DECLARE_SOA_COLUMN(SourcePhotonUid, sourcePhotonUid, int64_t);    //! mixer photon counter of the run-1 photon the leg belonged to
DECLARE_SOA_COLUMN(SourceV0Id, sourceV0Id, int);                  //! v0Id of that photon in run 1 (< 0: not an SVertexer V0)
DECLARE_SOA_COLUMN(SourcePx, sourcePx, float);                    //! run-1 photon momentum at the PV
DECLARE_SOA_COLUMN(SourcePy, sourcePy, float);                    //!
DECLARE_SOA_COLUMN(SourcePz, sourcePz, float);                    //!
DECLARE_SOA_COLUMN(SourceVx, sourceVx, float);                    //! run-1 conversion point
DECLARE_SOA_COLUMN(SourceVy, sourceVy, float);                    //!
DECLARE_SOA_COLUMN(SourceVz, sourceVz, float);                    //!
DECLARE_SOA_COLUMN(SourcePvX, sourcePvX, float);                  //! primary vertex of the source event (run 2: refit of the original at its own vertex)
DECLARE_SOA_COLUMN(SourcePvY, sourcePvY, float);                  //!
DECLARE_SOA_COLUMN(SourcePvZ, sourcePvZ, float);                  //!
DECLARE_SOA_COLUMN(SourcePvCovXX, sourcePvCovXX, float);          //!
DECLARE_SOA_COLUMN(SourcePvCovXY, sourcePvCovXY, float);          //!
DECLARE_SOA_COLUMN(SourcePvCovXZ, sourcePvCovXZ, float);          //!
DECLARE_SOA_COLUMN(SourcePvCovYY, sourcePvCovYY, float);          //!
DECLARE_SOA_COLUMN(SourcePvCovYZ, sourcePvCovYZ, float);          //!
DECLARE_SOA_COLUMN(SourcePvCovZZ, sourcePvCovZZ, float);          //!
DECLARE_SOA_COLUMN(SourcePvChi2, sourcePvChi2, float);            //!
DECLARE_SOA_COLUMN(SourcePvNContrib, sourcePvNContrib, uint16_t); //!
DECLARE_SOA_COLUMN(OrigTrackId, origTrackId, int);                //! AO2D track (DF of the source event)
DECLARE_SOA_COLUMN(McParticleId, mcParticleId, int);              //! MC: particle of the leg (DF of the source event), -1 without
DECLARE_SOA_COLUMN(McMotherId, mcMotherId, int);                  //! MC: first mother (DF of the source event), -1 without
DECLARE_SOA_COLUMN(McPdg, mcPdg, int);                            //! MC: PDG of the leg
DECLARE_SOA_COLUMN(McMotherPdg, mcMotherPdg, int);                //! MC: PDG of the mother
DECLARE_SOA_COLUMN(McMotherPx, mcMotherPx, float);                //! MC: true momentum of the mother
DECLARE_SOA_COLUMN(McMotherPy, mcMotherPy, float);                //!
DECLARE_SOA_COLUMN(McMotherPz, mcMotherPz, float);                //!
DECLARE_SOA_COLUMN(McGrandMotherId, mcGrandMotherId, int);        //! MC: first mother of the mother (DF of the source event), -1 without
DECLARE_SOA_COLUMN(McGrandMotherPdg, mcGrandMotherPdg, int);      //! MC: PDG of the grandmother
DECLARE_SOA_COLUMN(McAncestorId, mcAncestorId, int);              //! MC: first physical primary along the mothers of the leg's mother (itself if primary), -1 without
DECLARE_SOA_COLUMN(McGGMotherId, mcGGMotherId, int);              //! MC: mother of the grandmother (e.g. the eta of eta -> 3 pi0), -1 without
DECLARE_SOA_COLUMN(McGGMotherPdg, mcGGMotherPdg, int);            //! MC: PDG of that particle
DECLARE_SOA_COLUMN(McRootId, mcRootId, int);                      //! MC: last particle along the first mothers of the leg's mother (no mother), -1 without
DECLARE_SOA_COLUMN(McRootPdg, mcRootPdg, int);                    //! MC: PDG of the root
DECLARE_SOA_COLUMN(McRootDepth, mcRootDepth, int);                //! MC: generations from the leg's mother to the root
DECLARE_SOA_COLUMN(McVx, mcVx, float);                            //! MC: production vertex of the leg = true conversion point (lab frame of the source event; translated with the photon), 0 without
DECLARE_SOA_COLUMN(McVy, mcVy, float);                            //!
DECLARE_SOA_COLUMN(McVz, mcVz, float);                            //!
} // namespace pcmmixleg
// one row per leg of an artificial event; the columns the builder reads from a track carry the same names as in the AO2D
DECLARE_SOA_TABLE(PCMMixLegs, "AOD", "PCMMIXLEG", //!
                  o2::soa::Index<>, pcmmixleg::CollisionId, v0leg::Sign,
                  track::X, track::Alpha, track::Y, track::Z, track::Snp, track::Tgl, track::Signed1Pt, pcmmixleg::Pt, pcmmixleg::Phi,
                  pcmmixleg::CYY, pcmmixleg::CZY, pcmmixleg::CZZ, pcmmixleg::CSnpY, pcmmixleg::CSnpZ, pcmmixleg::CSnpSnp,
                  pcmmixleg::CTglY, pcmmixleg::CTglZ, pcmmixleg::CTglSnp, pcmmixleg::CTglTgl,
                  pcmmixleg::C1PtY, pcmmixleg::C1PtZ, pcmmixleg::C1PtSnp, pcmmixleg::C1PtTgl, pcmmixleg::C1Pt21Pt2,
                  track::TPCNClsFindable, track::TPCNClsFindableMinusFound, track::TPCNClsFindableMinusCrossedRows, track::TPCNClsShared,
                  track::TPCChi2NCl, track::TPCInnerParam, track::TPCSignal, pidtpc::TPCNSigmaEl, pidtpc::TPCNSigmaPi,
                  track::ITSClusterSizes, track::ITSChi2NCl, track::DetectorMap,
                  pcmmixleg::SourceEventUid, pcmmixleg::SourcePhotonUid, pcmmixleg::SourceV0Id,
                  pcmmixleg::SourcePx, pcmmixleg::SourcePy, pcmmixleg::SourcePz, pcmmixleg::SourceVx, pcmmixleg::SourceVy, pcmmixleg::SourceVz,
                  pcmmixleg::SourcePvX, pcmmixleg::SourcePvY, pcmmixleg::SourcePvZ, pcmmixleg::SourcePvCovXX, pcmmixleg::SourcePvCovXY, pcmmixleg::SourcePvCovXZ,
                  pcmmixleg::SourcePvCovYY, pcmmixleg::SourcePvCovYZ, pcmmixleg::SourcePvCovZZ, pcmmixleg::SourcePvChi2, pcmmixleg::SourcePvNContrib,
                  pcmmixleg::OrigTrackId, pcmmixleg::McParticleId, pcmmixleg::McMotherId, pcmmixleg::McPdg, pcmmixleg::McMotherPdg,
                  pcmmixleg::McMotherPx, pcmmixleg::McMotherPy, pcmmixleg::McMotherPz, pcmmixleg::McGrandMotherId, pcmmixleg::McGrandMotherPdg, pcmmixleg::McAncestorId,
                  pcmmixleg::McGGMotherId, pcmmixleg::McGGMotherPdg, pcmmixleg::McRootId, pcmmixleg::McRootPdg, pcmmixleg::McRootDepth, pcmmixleg::McVx, pcmmixleg::McVy, pcmmixleg::McVz,

                  // dynamic column
                  track::TPCNClsFound<track::TPCNClsFindable, track::TPCNClsFindableMinusFound>,
                  track::TPCNClsCrossedRows<track::TPCNClsFindable, track::TPCNClsFindableMinusCrossedRows>,
                  track::TPCCrossedRowsOverFindableCls<track::TPCNClsFindable, track::TPCNClsFindableMinusCrossedRows>,
                  track::TPCFoundOverFindableCls<track::TPCNClsFindable, track::TPCNClsFindableMinusFound>,
                  track::TPCFractionSharedCls<track::TPCNClsShared, track::TPCNClsFindable, track::TPCNClsFindableMinusFound>,
                  track::v001::ITSClusterMap<track::ITSClusterSizes>, track::v001::ITSNCls<track::ITSClusterSizes>, track::v001::ITSNClsInnerBarrel<track::ITSClusterSizes>,
                  track::HasITS<track::DetectorMap>, track::HasTPC<track::DetectorMap>, track::HasTRD<track::DetectorMap>, track::HasTOF<track::DetectorMap>);
using PCMMixLeg = PCMMixLegs::iterator;

// free electron legs of builder run 1 (fillFreeLegs): tracks that pass the leg cuts and are the leg of no stored photon; the
// parameters as the V0 fit would use them (TPC-only legs moved to their own collision). The mixer pools them with the photons
// (cfgPoolFreeLegs), so that run 2 can pair them in the same and the mixed event alike.
namespace pcmfreeleg
{
DECLARE_SOA_COLUMN(CollisionId, collisionId, int); //! collision of the track
DECLARE_SOA_COLUMN(TrackId, trackId, int);         //! AO2D track
} // namespace pcmfreeleg
DECLARE_SOA_TABLE(PCMFreeLegs, "AOD", "PCMFREELEG", //!
                  o2::soa::Index<>, pcmfreeleg::CollisionId, pcmfreeleg::TrackId, v0leg::Sign,
                  v0legcov::X, v0legcov::Alpha, v0legcov::Y, v0legcov::Z, v0legcov::Snp, v0legcov::Tgl, v0legcov::Signed1Pt, v0legcov::Cov,
                  track::TPCNClsFindable, track::TPCNClsFindableMinusFound, track::TPCNClsFindableMinusCrossedRows, track::TPCNClsShared,
                  track::TPCChi2NCl, track::TPCInnerParam, track::TPCSignal, pidtpc::TPCNSigmaEl, pidtpc::TPCNSigmaPi,
                  track::ITSClusterSizes, track::ITSChi2NCl, track::DetectorMap,

                  // dynamic column
                  track::HasITS<track::DetectorMap>, track::HasTPC<track::DetectorMap>);
using PCMFreeLeg = PCMFreeLegs::iterator;

// free legs as AO2D track indices only (linked derived data): run 1 writes (collision, track) of every free leg of a selected
// collision; associate-photon-to-emevent adds the derived event row. The mixer on linked derived data reads the tracks from the parent AO2D.
namespace pcmfreelegid
{
DECLARE_SOA_INDEX_COLUMN(PMEvent, pmevent);          //! derived event of the free leg (-1: collision not in PMEvents)
DECLARE_SOA_COLUMN(TrackIdDelta, trackIdDelta, int); //! AO2D track index minus the one of the previous row of this DF (first row: the index): small numbers, ~5x smaller on disk
} // namespace pcmfreelegid
DECLARE_SOA_TABLE(PCMFreeLegIds, "AOD", "PCMFREELEGID", //! rows in increasing track order; track index = running sum of TrackIdDelta over the DF
                  o2::soa::Index<>, pcmfreeleg::CollisionId, pcmfreelegid::TrackIdDelta);
using PCMFreeLegId = PCMFreeLegIds::iterator;
DECLARE_SOA_TABLE(PCMFreeLegEMEventIds, "AOD", "PCMFREELEGEVID", pcmfreelegid::PMEventId); //! to be joined with PCMFreeLegIds

// builder run 2 output (same columns as the run-1 tables; collisionId = PCMMixCollisions row, trackId = PCMMixLegs row,
// v0Id >= 0: original pairing of run 1, -1: crossed pairing)
DECLARE_SOA_TABLE(V0LegsMix, "AOD", "V0LEGMIX", //!
                  o2::soa::Index<>, v0leg::CollisionId, v0leg::TrackId, v0leg::Sign,
                  v0leg::Px, v0leg::Py, v0leg::Pz,
                  track::DcaXY, track::DcaZ,
                  track::TPCNClsFindable, track::TPCNClsFindableMinusFound, track::TPCNClsFindableMinusCrossedRows, track::TPCNClsShared,
                  track::TPCChi2NCl, track::TPCInnerParam,
                  track::TPCSignal, pidtpc::TPCNSigmaEl, pidtpc::TPCNSigmaPi,
                  track::ITSClusterSizes, track::ITSChi2NCl, track::DetectorMap,

                  // dynamic column
                  v0leg::P<v0leg::Px, v0leg::Py, v0leg::Pz>,
                  v0leg::Pt<v0leg::Px, v0leg::Py>,
                  v0leg::Eta<v0leg::Px, v0leg::Py, v0leg::Pz>,
                  v0leg::Phi<v0leg::Px, v0leg::Py>,
                  v0leg::Tgl<v0leg::Px, v0leg::Py, v0leg::Pz>,

                  track::TPCNClsFound<track::TPCNClsFindable, track::TPCNClsFindableMinusFound>,
                  track::TPCNClsCrossedRows<track::TPCNClsFindable, track::TPCNClsFindableMinusCrossedRows>,
                  track::TPCCrossedRowsOverFindableCls<track::TPCNClsFindable, track::TPCNClsFindableMinusCrossedRows>,
                  track::TPCFoundOverFindableCls<track::TPCNClsFindable, track::TPCNClsFindableMinusFound>,
                  track::TPCFractionSharedCls<track::TPCNClsShared, track::TPCNClsFindable, track::TPCNClsFindableMinusFound>,
                  track::v001::ITSClusterMap<track::ITSClusterSizes>, track::v001::ITSNCls<track::ITSClusterSizes>, track::v001::ITSNClsInnerBarrel<track::ITSClusterSizes>,
                  track::HasITS<track::DetectorMap>, track::HasTPC<track::DetectorMap>, track::HasTRD<track::DetectorMap>, track::HasTOF<track::DetectorMap>,
                  v0leg::MeanClusterSizeITS<track::ITSClusterSizes>,
                  v0leg::MeanClusterSizeITSib<track::ITSClusterSizes>,
                  v0leg::MeanClusterSizeITSob<track::ITSClusterSizes>);
using V0LegMix = V0LegsMix::iterator;

namespace v0photonkfmix
{
DECLARE_SOA_INDEX_COLUMN_FULL(PosTrack, posTrack, int, V0LegsMix, "_Pos");                            //!
DECLARE_SOA_INDEX_COLUMN_FULL(NegTrack, negTrack, int, V0LegsMix, "_Neg");                            //!
DECLARE_SOA_DYNAMIC_COLUMN(PMEventId, pmeventId, [](int collisionId) -> int { return collisionId; }); //! event index as in V0KFEMEventIds (= PCMMixCollisions row), for V0PhotonCut
} // namespace v0photonkfmix
DECLARE_SOA_TABLE(V0PhotonsKFMix, "AOD", "V0PHOTONKFMIX", //!
                  o2::soa::Index<>, v0photonkf::CollisionId, v0photonkf::V0Id, v0photonkfmix::PosTrackId, v0photonkfmix::NegTrackId,
                  v0photonkf::Vx, v0photonkf::Vy, v0photonkf::Vz,
                  v0photonkf::Px, v0photonkf::Py, v0photonkf::Pz,
                  v0photonkf::MGamma,
                  v0photonkf::DCAxyToPV, v0photonkf::DCAzToPV,
                  v0photonkf::CosPA, v0photonkf::CosPAXY, v0photonkf::CosPARZ, v0photonkf::PCA,
                  v0photonkf::Alpha, v0photonkf::QtArm,
                  v0photonkf::ChiSquareNDF,

                  // dynamic column
                  v0photonkf::E<v0photonkf::Px, v0photonkf::Py, v0photonkf::Pz>,
                  v0photonkf::Pt<v0photonkf::Px, v0photonkf::Py>,
                  v0photonkf::Eta<v0photonkf::Px, v0photonkf::Py, v0photonkf::Pz>,
                  v0photonkf::Phi<v0photonkf::Px, v0photonkf::Py>,
                  v0photonkf::P<v0photonkf::Px, v0photonkf::Py, v0photonkf::Pz>,
                  v0photonkf::V0Radius<v0photonkf::Vx, v0photonkf::Vy>,
                  v0photonkfmix::PMEventId<v0photonkf::CollisionId>);
using V0PhotonKFMix = V0PhotonsKFMix::iterator;

DECLARE_SOA_TABLE(V0LegsXYZMix, "AOD", "V0LEGXYZMIX", track::X, track::Y, track::Z);                                                   //!
DECLARE_SOA_TABLE(V0LegsDeDxMCMix, "AOD", "V0LEGDEDXMCMIX", mcpidtpc::DeDxTunedMc, o2::soa::Marker<3>);                                //! stays empty (run 2 runs without MC)
DECLARE_SOA_TABLE(V0PhotonsPhiVPsiMix, "AOD", "V0PHOTONPHIVMIX", o2::soa::Index<>, v0photonsphivpsi::PhiV, v0photonsphivpsi::PsiPair); //!
DECLARE_SOA_TABLE(V0LegsCovMix, "AOD", "V0LEGCOVMIX",                                                                                  //!
                  v0legcov::X, v0legcov::Alpha, v0legcov::Y, v0legcov::Z, v0legcov::Snp, v0legcov::Tgl, v0legcov::Signed1Pt, v0legcov::Cov);

namespace dalitzee
{
DECLARE_SOA_INDEX_COLUMN(EMEvent, emevent);                                         //!
DECLARE_SOA_INDEX_COLUMN(PMEvent, pmevent);                                         //!
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

DECLARE_SOA_TABLE(DalitzEEEMEventIds_000, "AOD", "EEEMEVENTID", dalitzee::EMEventId);              // To be joined with DalitzEEs table at analysis level.
DECLARE_SOA_TABLE_VERSIONED(DalitzEEEMEventIds_001, "AOD", "EEEMEVENTID", 1, dalitzee::PMEventId); // To be joined with DalitzEEs table at analysis level.
using DalitzEEEMEventIds = DalitzEEEMEventIds_001;
using DalitzEEEMEventId = DalitzEEEMEventIds::iterator;

namespace pwgem::photon::swtinfo
{
DECLARE_SOA_INDEX_COLUMN(EMEvent, emevent);                                                                              //!
DECLARE_SOA_INDEX_COLUMN(PMEvent, pmevent);                                                                              //!
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
DECLARE_SOA_INDEX_COLUMN(PMEvent, pmevent);                            //!
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
DECLARE_SOA_INDEX_COLUMN(EMEvent, emevent);                                                                                             //!
DECLARE_SOA_INDEX_COLUMN(PMEvent, pmevent);                                                                                             //!
DECLARE_SOA_COLUMN(CoreEnergy, coreEnergy, float);                                                                                      //! cluster core energy (GeV)
DECLARE_SOA_COLUMN(Time, time, float);                                                                                                  //! cluster time (ns)
DECLARE_SOA_COLUMN(IsExotic, isExotic, bool);                                                                                           //! flag to mark cluster as exotic
DECLARE_SOA_COLUMN(Definition, definition, int);                                                                                        //! cluster definition, see EMCALClusterDefinition.h
DECLARE_SOA_ARRAY_INDEX_COLUMN(Track, track);                                                                                           //! TrackIds
DECLARE_SOA_COLUMN(DeltaPhi, deltaPhi, std::vector<float>);                                                                             //! phi values of the matched tracks
DECLARE_SOA_COLUMN(DeltaEta, deltaEta, std::vector<float>);                                                                             //! eta values of the matched tracks
DECLARE_SOA_COLUMN(TrackP, trackp, std::vector<float>);                                                                                 //! momentum values of the matched tracks
DECLARE_SOA_COLUMN(TrackPt, trackpt, std::vector<float>);                                                                               //! pt values of the matched tracks
DECLARE_SOA_COLUMN(DeltaPhiSec, deltaPhiSec, std::vector<float>);                                                                       //! phi values of the matched secondary tracks
DECLARE_SOA_COLUMN(DeltaEtaSec, deltaEtaSec, std::vector<float>);                                                                       //! eta values of the matched secondary tracks
DECLARE_SOA_COLUMN(TrackPSec, trackpSec, std::vector<float>);                                                                           //! momentum values of the matched secondary tracks
DECLARE_SOA_COLUMN(TrackPtSec, trackptSec, std::vector<float>);                                                                         //! pt values of the matched secondary tracks
DECLARE_SOA_DYNAMIC_COLUMN(Pt, pt, [](float e, float eta, float m = 0) -> float { return std::sqrt(e * e - m * m) / std::cosh(eta); }); //! cluster pt, mass to be given as argument when getter is called!
} // namespace emccluster

namespace emctm
{
DECLARE_SOA_COLUMN(DeltaPhi, deltaPhi, float); //! phi values of the matched tracks
DECLARE_SOA_COLUMN(DeltaEta, deltaEta, float); //! eta values of the matched tracks
DECLARE_SOA_COLUMN(TrackP, trackp, float);     //! momentum values of the matched tracks
DECLARE_SOA_COLUMN(TrackPt, trackpt, float);   //! pt values of the matched tracks
} // namespace emctm
DECLARE_SOA_TABLE_VERSIONED(SkimEMCClusters_000, "AOD", "SKIMEMCCLUSTER", 0, //! table of skimmed EMCal clusters
                            o2::soa::Index<>, skimmedcluster::CollisionId, emccluster::Definition, skimmedcluster::E, skimmedcluster::Eta, skimmedcluster::Phi,
                            skimmedcluster::M02, skimmedcluster::NCells, skimmedcluster::Time, emccluster::IsExotic, emccluster::DeltaPhi,
                            emccluster::DeltaEta, emccluster::TrackP, emccluster::TrackPt, emccluster::Pt<skimmedcluster::E, skimmedcluster::Eta>);

DECLARE_SOA_TABLE_VERSIONED(SkimEMCClusters_001, "AOD", "SKIMEMCCLUSTER", 1, //! table of skimmed EMCal clusters version 002 - including secondary matched tracks
                            o2::soa::Index<>, skimmedcluster::CollisionId, emccluster::Definition, skimmedcluster::E, skimmedcluster::Eta, skimmedcluster::Phi,
                            skimmedcluster::M02, skimmedcluster::NCells, skimmedcluster::Time, emccluster::IsExotic,
                            emccluster::DeltaPhi, emccluster::DeltaEta, emccluster::TrackP, emccluster::TrackPt,
                            emccluster::DeltaPhiSec, emccluster::DeltaEtaSec, emccluster::TrackPSec, emccluster::TrackPtSec,
                            emccluster::Pt<skimmedcluster::E, skimmedcluster::Eta>);

using SkimEMCClusters = SkimEMCClusters_001;
using SkimEMCCluster = SkimEMCClusters_001::iterator;

DECLARE_SOA_TABLE_VERSIONED(EmEmcClusters_000, "AOD", "EMEMCCLUSTER", 0, //! table of skimmed EMCal clusters
                            o2::soa::Index<>, skimmedcluster::CollisionId, emccluster::Definition, skimmedcluster::E, skimmedcluster::Eta, skimmedcluster::Phi,
                            skimmedcluster::M02, skimmedcluster::NCells, skimmedcluster::Time, emccluster::IsExotic, emccluster::Pt<skimmedcluster::E, skimmedcluster::Eta>);

using EmEmcClusters = EmEmcClusters_000;
using EmEmcCluster = EmEmcClusters_000::iterator;

namespace trackmatching
{
DECLARE_SOA_INDEX_COLUMN(EmEmcCluster, emEmcCluster); //!
} // namespace trackmatching

DECLARE_SOA_TABLE(EmEmcMTracks, "AOD", "EMEMCMTRACK", //!
                  o2::soa::Index<>, trackmatching::EmEmcClusterId, emctm::DeltaPhi, emctm::DeltaEta, emctm::TrackP, emctm::TrackPt);

DECLARE_SOA_TABLE(EmEmcMSTracks, "AOD", "EMEMCMSTRACK", //!
                  o2::soa::Index<>, trackmatching::EmEmcClusterId, emctm::DeltaPhi, emctm::DeltaEta, emctm::TrackP, emctm::TrackPt);

DECLARE_SOA_TABLE(EMCEMEventIds_000, "AOD", "EMCEMEVENTID", emccluster::EMEventId);              // To be joined with SkimEMCClusters table at analysis level.
DECLARE_SOA_TABLE_VERSIONED(EMCEMEventIds_001, "AOD", "EMCEMEVENTID", 1, emccluster::PMEventId); // To be joined with SkimEMCClusters table at analysis level.
using EMCEMEventIds = EMCEMEventIds_001;
using EMCEMEventId = EMCEMEventIds::iterator;

namespace phoscluster
{
DECLARE_SOA_INDEX_COLUMN(EMEvent, emevent);                                         //!
DECLARE_SOA_INDEX_COLUMN(PMEvent, pmevent);                                         //!
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
DECLARE_SOA_DYNAMIC_COLUMN(Px, px, [](float e, float x, float y, float z, float m = 0) -> float { return x / RecoDecay::sqrtSumOfSquares(x, y, z) * std::sqrt(e * e - m * m); });
DECLARE_SOA_DYNAMIC_COLUMN(Py, py, [](float e, float x, float y, float z, float m = 0) -> float { return y / RecoDecay::sqrtSumOfSquares(x, y, z) * std::sqrt(e * e - m * m); });
DECLARE_SOA_DYNAMIC_COLUMN(Pz, pz, [](float e, float x, float y, float z, float m = 0) -> float { return z / RecoDecay::sqrtSumOfSquares(x, y, z) * std::sqrt(e * e - m * m); });
DECLARE_SOA_DYNAMIC_COLUMN(Pt, pt, [](float e, float x, float y, float z, float m = 0) -> float { return RecoDecay::sqrtSumOfSquares(x, y) / RecoDecay::sqrtSumOfSquares(x, y, z) * std::sqrt(e * e - m * m); });
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

DECLARE_SOA_TABLE(PHOSEMEventIds_000, "AOD", "PHOSEMEVENTID", phoscluster::EMEventId);              // To be joined with PHOSClusters table at analysis level.
DECLARE_SOA_TABLE_VERSIONED(PHOSEMEventIds_001, "AOD", "PHOSEMEVENTID", 1, phoscluster::PMEventId); // To be joined with PHOSClusters table at analysis level.
using PHOSEMEventIds = PHOSEMEventIds_001;
using PHOSEMEventId = PHOSEMEventIds::iterator;

namespace caloextra
{
DECLARE_SOA_INDEX_COLUMN_FULL(Cluster, cluster, int, SkimEMCClusters, ""); //! reference to the gamma in the skimmed EMCal table
DECLARE_SOA_INDEX_COLUMN_FULL(Cell, cell, int, Calos, "");                 //! reference to the gamma in the skimmed EMCal table
} // namespace caloextra

DECLARE_SOA_TABLE(SkimEMCCells, "AOD", "SKIMEMCCELLS",                        //! table of link between skimmed EMCal clusters and their cells
                  o2::soa::Index<>, caloextra::ClusterId, caloextra::CellId); //!
using SkimEMCCell = SkimEMCCells::iterator;

namespace nonlin
{
DECLARE_SOA_COLUMN(CorrE, corrE, float);   //! non lin corrected energy
DECLARE_SOA_COLUMN(CorrP, corrP, float);   //! non lin corrected momentum
DECLARE_SOA_COLUMN(CorrPt, corrPt, float); //! non lin corrected transverse momentum
} // namespace nonlin

DECLARE_SOA_TABLE(NonLinV0s, "AOD", "NONLINV0", //! table of non lin corrected values for V0 Photons (so far only pt)
                  nonlin::CorrPt);              //!

DECLARE_SOA_TABLE(NonLinEmcClusters, "AOD", "NONLINEMCCLUSTER", //! table of non lin corrected values for EMCal Photons (so far only E and pT)
                  nonlin::CorrE, nonlin::CorrPt);               //!

namespace v0photonMBweights
{
DECLARE_SOA_COLUMN(OmegaMBWeight, omegaMBWeight, float);
}

DECLARE_SOA_TABLE(V0PhotonOmegaMBWeights, "AOD", "V0PHOTONMBW", v0photonMBweights::OmegaMBWeight); // store MB weights. To be joined with V0PhotonsKF table at analysis level.

using V0PhotonOmegaMBWeight = V0PhotonOmegaMBWeights::iterator;

} // namespace o2::aod

#endif // PWGEM_PHOTONMESON_DATAMODEL_GAMMATABLES_H_
