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

///
/// \file   JetTrackQa.h
/// \author Johanna Lömker
/// \since  2023-10-02
/// \brief  Header for the trackJetQa task for the analysis of the tracks for jets.
///

#ifndef PWGJE_DATAMODEL_TRACKJETQA_H_
#define PWGJE_DATAMODEL_TRACKJETQA_H_

// O2 includes
#include "ReconstructionDataFormats/Track.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/AnalysisDataModel.h"
#include "PWGJE/DataModel/Jet.h"
#include "Framework/StaticFor.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"

using namespace o2;
using namespace o2::track;
using namespace o2::framework;
using namespace o2::framework::expressions;

// Derived data model for cut variation
namespace o2::aod
{
namespace jetspectra
{

template <typename binningType>
typename binningType::binned_t packInTable(const float& valueToBin)
{
  if (valueToBin <= binningType::binned_min) {
    return (binningType::underflowBin);
  } else if (valueToBin >= binningType::binned_max) {
    return (binningType::overflowBin);
  } else if (valueToBin >= 0) {
    return (static_cast<typename binningType::binned_t>((valueToBin / binningType::bin_width) + 0.5f));
  } else {
    return (static_cast<typename binningType::binned_t>((valueToBin / binningType::bin_width) - 0.5f));
  }
}
// Function to unpack a binned value into a float
template <typename binningType>
float unPack(const typename binningType::binned_t& valueToUnpack)
{
  return binningType::bin_width * static_cast<float>(valueToUnpack);
}

struct binningDCA {
 public:
  typedef int16_t binned_t;
  static constexpr int nbins = (1 << 8 * sizeof(binned_t)) - 2;
  static constexpr binned_t overflowBin = nbins >> 1;
  static constexpr binned_t underflowBin = -(nbins >> 1);
  static constexpr float binned_max = 6.0;
  static constexpr float binned_min = -6.0;
  static constexpr float bin_width = (binned_max - binned_min) / nbins;
};

// Collision info
DECLARE_SOA_INDEX_COLUMN(BC, bc); //! Most probably BC to where this collision has occurred
DECLARE_SOA_COLUMN(RunNumber, runNumber, int);
DECLARE_SOA_COLUMN(Sel8, sel8, bool);
DECLARE_SOA_COLUMN(MultNTracksPV, multNTracksPV, int);
DECLARE_SOA_COLUMN(MultTracklets, multTracklets, int);
DECLARE_SOA_COLUMN(MultFT0M, multFT0M, float);
DECLARE_SOA_COLUMN(CentFT0M, centFT0M, float);
//  Track info
DECLARE_SOA_INDEX_COLUMN(Collision, collision); //! Index to the collision
DECLARE_SOA_COLUMN(Signed1Pt, signed1Pt, float); //! Pt (signed) of the track
DECLARE_SOA_COLUMN(Eta, eta, float);            //! Eta of the track
DECLARE_SOA_COLUMN(Phi, phi, float);            //! Phi of the track
DECLARE_SOA_COLUMN(Sigma1Pt, sigma1Pt, float);
DECLARE_SOA_COLUMN(Alpha, alpha, float);
DECLARE_SOA_COLUMN(X, x, float);
DECLARE_SOA_COLUMN(Y, y, float);
DECLARE_SOA_COLUMN(Z, z, float);
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Snp, snp, float);
DECLARE_SOA_COLUMN(Tgl, tgl, float);
// DECLARE_SOA_COLUMN(EvTimeT0AC, evTimeT0AC, float);                               //! Event time of the track computed with the T0AC
// DECLARE_SOA_COLUMN(EvTimeT0ACErr, evTimeT0ACErr, float);                         //! Resolution of the event time of the track computed with the T0AC
DECLARE_SOA_COLUMN(IsPVContributor, isPVContributor, bool); //! IsPVContributor
DECLARE_SOA_COLUMN(HasTRD, hasTRD, bool); //! Has or not the TRD match

DECLARE_SOA_COLUMN(DCAxyStore, dcaxyStore, binningDCA::binned_t); //! Stored binned dcaxy
DECLARE_SOA_COLUMN(DCAzStore, dcazStore, binningDCA::binned_t);   //! Stored binned dcaz
DECLARE_SOA_DYNAMIC_COLUMN(DCAxy, dcaXY,                          //! Unpacked dcaxy
                           [](binningDCA::binned_t binned) -> float { return unPack<binningDCA>(binned); });
DECLARE_SOA_DYNAMIC_COLUMN(DCAz, dcaZ, //! Unpacked dcaz
                           [](binningDCA::binned_t binned) -> float { return unPack<binningDCA>(binned); });
// DECLARE_SOA_DYNAMIC_COLUMN(Pt, pt, //! Absolute value of signed pT
//                            [](float signedPt) -> float { return std::abs(signedPt); });
DECLARE_SOA_DYNAMIC_COLUMN(P, p, [](float signed1pt, float eta) -> float { return std::abs(signed1pt) * cosh(eta); });
DECLARE_SOA_DYNAMIC_COLUMN(TrackType, trackType, [](float v) -> uint8_t { return o2::aod::track::TrackTypeEnum::Track; });
DECLARE_SOA_COLUMN(IsGlobalTrack, isGlobalTrack, bool);                                   // if a track passed the isGlobalTrack requirement
DECLARE_SOA_COLUMN(IsGlobalTrackWoDCA, isGlobalTrackWoDCA, bool);                         // if a track passed the isGlobalTrackWoDCA requirement
DECLARE_SOA_COLUMN(IsGlobalTrackWoPtEta, isGlobalTrackWoPtEta, bool);                     // if a track passed the isGlobalTrackWoDCA requirement
DECLARE_SOA_DYNAMIC_COLUMN(Flags, flags, [](float v) -> uint32_t { return 0; });          // Dummy
DECLARE_SOA_DYNAMIC_COLUMN(Rapidity, rapidity,                                            //! Track rapidity, computed under the mass assumption given as input
                           [](float signed1Pt, float eta, float mass) -> float {
                             const auto pt = std::abs(signed1Pt);
                             const auto p = std::abs(signed1Pt) * cosh(eta);
                             const auto pz = std::sqrt(p * p - pt * pt);
                             const auto energy = sqrt(p * p + mass * mass);
                             return 0.5f * log((energy + pz) / (energy - pz));
                           });
DECLARE_SOA_DYNAMIC_COLUMN(IsQualityTrackITS, isQualityTrackITS, [](float v) -> bool { return false; }); // Dummy
DECLARE_SOA_DYNAMIC_COLUMN(IsQualityTrackTPC, isQualityTrackTPC, [](float v) -> bool { return false; }); // Dummy

DECLARE_SOA_COLUMN(Length, length, float);
DECLARE_SOA_COLUMN(TPCChi2NCl, tpcChi2NCl, float);
DECLARE_SOA_COLUMN(ITSChi2NCl, itsChi2NCl, float);
DECLARE_SOA_COLUMN(TOFChi2, tofChi2, float);
DECLARE_SOA_COLUMN(TPCNClsShared, tpcNClsShared, int);
DECLARE_SOA_COLUMN(TPCNClsFindable, tpcNClsFindable, int);
DECLARE_SOA_COLUMN(TPCNClsFindableMinusFound, tpcNClsFindableMinusFound, int);
DECLARE_SOA_COLUMN(TPCNClsFindableMinusCrossedRows, tpcNClsFindableMinusCrossedRows, int);
DECLARE_SOA_COLUMN(ITSNCls, itsNCls, int);
DECLARE_SOA_COLUMN(ITSNClsInnerBarrel, itsNClsInnerBarrel, int);
DECLARE_SOA_COLUMN(TPCFractionSharedCls, tpcFractionSharedCls, float);
DECLARE_SOA_COLUMN(ITSClusterMap, itsClusterMap, int);
DECLARE_SOA_COLUMN(TPCNClsFound, tpcNClsFound, int);
DECLARE_SOA_COLUMN(TPCNClsCrossedRows, tpcNClsCrossedRows, int);
DECLARE_SOA_COLUMN(TPCCrossedRowsOverFindableCls, tpcCrossedRowsOverFindableCls, float);
DECLARE_SOA_COLUMN(TPCFoundOverFindableCls, tpcFoundOverFindableCls, float);

} // namespace jetspectra

DECLARE_SOA_TABLE(JeColls, "AOD", "JECOLLS",
                  o2::soa::Index<>,
                  collision::NumContrib,
                  collision::PosX,
                  collision::PosY,
                  collision::PosZ,
                  jetspectra::Sel8,
                  jetspectra::MultNTracksPV,
                  jetspectra::MultTracklets,
                  jetspectra::MultFT0M,
                  jetspectra::CentFT0M,
                  jetspectra::RunNumber);
using JeColl = JeColls::iterator;

DECLARE_SOA_TABLE(JeTracks, "AOD", "JETRACKS",
                  o2::soa::Index<>,
                  jetspectra::CollisionId,
                  jetspectra::Signed1Pt, jetspectra::Eta, jetspectra::Phi, jetspectra::Pt,
                  jetspectra::Sigma1Pt,
                  jetspectra::Alpha,
                  jetspectra::X,
                  jetspectra::Y,
                  jetspectra::Z,
                  jetspectra::Snp,
                  jetspectra::Tgl,
                  jetspectra::IsPVContributor,
                  jetspectra::HasTRD,
                  jetspectra::DCAxyStore,
                  jetspectra::DCAzStore,
                  jetspectra::DCAxy<jetspectra::DCAxyStore>,
                  jetspectra::DCAz<jetspectra::DCAzStore>,
                  jetspectra::IsGlobalTrack,
                  jetspectra::IsGlobalTrackWoDCA,
                  jetspectra::IsGlobalTrackWoPtEta,
                  // jetspectra::Pt<jetspectra::PtSigned>,
                  // jetspectra::P<jetspectra::PtSigned, jetspectra::Eta>,
                  // jetspectra::Rapidity<jetspectra::PtSigned, jetspectra::Eta>,
                  jetspectra::Flags<track::TOFChi2>,
                  jetspectra::TrackType<track::TOFChi2>,
                  jetspectra::IsQualityTrackITS<track::TOFChi2>,
                  jetspectra::IsQualityTrackTPC<track::TOFChi2>,
                  // track::Sign<jetspectra::PtSigned>,
                  jetspectra::Length,
                  jetspectra::TPCChi2NCl, jetspectra::ITSChi2NCl, track::TOFChi2,
                  jetspectra::TPCNClsShared,
                  jetspectra::TPCNClsFindable,
                  jetspectra::TPCNClsFindableMinusFound,
                  jetspectra::TPCNClsFindableMinusCrossedRows,
                  jetspectra::ITSClusterMap,
                  jetspectra::ITSNCls,
                  jetspectra::ITSNClsInnerBarrel,
                  jetspectra::TPCFractionSharedCls,
                  jetspectra::TPCNClsFound,
                  jetspectra::TPCNClsCrossedRows,
                  jetspectra::TPCCrossedRowsOverFindableCls,
                  jetspectra::TPCFoundOverFindableCls);
} // namespace o2::aod

#endif // PWGJE_DATAMODEL_TRACKJETQA_H_
