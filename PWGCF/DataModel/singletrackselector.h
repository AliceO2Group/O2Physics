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
/// \brief  create a table for single track selection.
/// \author Sofia Tomassini
/// \since 30 May 2023

#ifndef PWGCF_DATAMODEL_SINGLETRACKSELECTOR_H_
#define PWGCF_DATAMODEL_SINGLETRACKSELECTOR_H_

#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "Common/DataModel/PIDResponse.h"
#include "Framework/Logger.h"
#include "Common/DataModel/Multiplicity.h"

namespace o2::aod
{
namespace singletrackselector
{
template <typename binningType>
typename binningType::binned_t packInTable(const float& valueToBin)
{
  if (valueToBin <= binningType::binned_min) {
    return binningType::underflowBin;
  } else if (valueToBin >= binningType::binned_max) {
    return binningType::overflowBin;
  } else {
    return static_cast<typename binningType::binned_t>(valueToBin / binningType::bin_width);
  }
}

template <typename binningType>
float unPack(const typename binningType::binned_t& b)
{
  return static_cast<float>(binningType::bin_width * b);
}

template <typename binningType>
typename binningType::binned_t packInTableOffset(const float& valueToBin)
{
  if (valueToBin <= binningType::binned_min) {
    return binningType::underflowBin;
  } else if (valueToBin >= binningType::binned_max) {
    return binningType::overflowBin;
  } else {
    return static_cast<typename binningType::binned_t>(((valueToBin - (binningType::binned_max - binningType::binned_min) * 0.5) / binningType::bin_width));
  }
}

template <typename binningType>
float unPackOffset(const typename binningType::binned_t& b)
{
  return static_cast<float>((binningType::binned_max - binningType::binned_min) * 0.5 + binningType::bin_width * b);
}

namespace storedcrossedrows
{
struct binning {
 public:
  typedef int8_t binned_t;
  static constexpr int nbins = (1 << 8 * sizeof(binned_t)) - 2;
  static constexpr binned_t overflowBin = nbins >> 1;
  static constexpr binned_t underflowBin = -(nbins >> 1);
  static constexpr float binned_max = 253.5;
  static constexpr float binned_min = -0.5;
  static constexpr float bin_width = (binned_max - binned_min) / nbins;
};
} // namespace storedcrossedrows

namespace nsigma
{
struct binning {
 public:
  typedef int8_t binned_t;
  static constexpr int nbins = (1 << 8 * sizeof(binned_t)) - 2;
  static constexpr binned_t overflowBin = nbins >> 1;
  static constexpr binned_t underflowBin = -(nbins >> 1);
  static constexpr float binned_max = 10.0;
  static constexpr float binned_min = -10.0;
  static constexpr float bin_width = (binned_max - binned_min) / nbins;
};
} // namespace nsigma

DECLARE_SOA_INDEX_COLUMN(Collision, collision); // Index to the collision
DECLARE_SOA_COLUMN(HasTOF, hasTOF, bool);
DECLARE_SOA_COLUMN(HasITS, hasITS, bool);
DECLARE_SOA_COLUMN(Px, px, float); // Momentum of the track
DECLARE_SOA_COLUMN(Py, py, float); // Momentum of the track
DECLARE_SOA_COLUMN(Pz, pz, float); // Momentum of the track
DECLARE_SOA_DYNAMIC_COLUMN(P, p,
                           [](float px, float py, float pz) -> float { return std::sqrt(px * px + py * py + pz * pz); }); // Momentum of the track
DECLARE_SOA_DYNAMIC_COLUMN(Pt, pt,
                           [](float px, float py) -> float { return std::sqrt(px * px + py * py); }); // Momentum of the track
DECLARE_SOA_COLUMN(TPCInnerParam, tpcInnerParam, float);                                              // vertex position along z
DECLARE_SOA_COLUMN(TPCSignal, tpcSignal, float);                                                      // vertex position along z
DECLARE_SOA_COLUMN(Beta, beta, float);
DECLARE_SOA_COLUMN(DcaXY, dcaXY, float);               // impact parameter of the track
DECLARE_SOA_COLUMN(DcaZ, dcaZ, float);                 // impact parameter of the track
DECLARE_SOA_COLUMN(TPCNClsFound, tpcNClsFound, float); // Number of TPC clusters
DECLARE_SOA_COLUMN(TPCCrossedRowsOverFindableCls, tpcCrossedRowsOverFindableCls, float);
DECLARE_SOA_COLUMN(TPCChi2NCl, tpcChi2NCl, float); // TPC chi2
DECLARE_SOA_COLUMN(ITSNCls, itsNCls, float);       // Number of ITS clusters
DECLARE_SOA_COLUMN(ITSChi2NCl, itsChi2NCl, float); // ITS chi2
DECLARE_SOA_COLUMN(Sign, sign, int8_t);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(StoredCrossedRows, storedCrossedRows, storedcrossedrows::binning::binned_t);
DECLARE_SOA_COLUMN(StoredTOFNSigmaPr, storedTofNSigmaPr, nsigma::binning::binned_t);
DECLARE_SOA_COLUMN(StoredTPCNSigmaPr, storedTpcNSigmaPr, nsigma::binning::binned_t);
DECLARE_SOA_COLUMN(StoredTOFNSigmaDe, storedTofNSigmaDe, nsigma::binning::binned_t);
DECLARE_SOA_COLUMN(StoredTPCNSigmaDe, storedTpcNSigmaDe, nsigma::binning::binned_t);

DECLARE_SOA_DYNAMIC_COLUMN(CrossedRows, tpcNClsCrossedRows,
                           [](storedcrossedrows::binning::binned_t binned) -> float { return singletrackselector::unPackOffset<storedcrossedrows::binning>(binned); });
DECLARE_SOA_DYNAMIC_COLUMN(TOFNSigmaPr, tofNSigmaPr,
                           [](nsigma::binning::binned_t nsigma_binned) -> float { return singletrackselector::unPack<nsigma::binning>(nsigma_binned); });
DECLARE_SOA_DYNAMIC_COLUMN(TPCNSigmaPr, tpcNSigmaPr,
                           [](nsigma::binning::binned_t nsigma_binned) -> float { return singletrackselector::unPack<nsigma::binning>(nsigma_binned); });
DECLARE_SOA_DYNAMIC_COLUMN(TOFNSigmaDe, tofNSigmaDe,
                           [](nsigma::binning::binned_t nsigma_binned) -> float { return singletrackselector::unPack<nsigma::binning>(nsigma_binned); });
DECLARE_SOA_DYNAMIC_COLUMN(TPCNSigmaDe, tpcNSigmaDe,
                           [](nsigma::binning::binned_t nsigma_binned) -> float { return singletrackselector::unPack<nsigma::binning>(nsigma_binned); });

DECLARE_SOA_DYNAMIC_COLUMN(Energy, energy,
                           [](float px, float py, float pz, float mass) -> float { return sqrt(px * px + py * py + pz * pz + mass * mass); });

DECLARE_SOA_COLUMN(GlobalIndex, globalIndex, int64_t); // Index to the collision
DECLARE_SOA_COLUMN(Mult, mult, int);                   // Multiplicity of the collision
DECLARE_SOA_COLUMN(PosZ, posZ, int);                   // Vertex of the collision

} // namespace singletrackselector

DECLARE_SOA_TABLE(SingleTrackSel, "AOD", "STSEL", // Table of the variables for single track selection.
                  o2::soa::Index<>,
                  singletrackselector::CollisionId,
                  singletrackselector::HasITS,
                  singletrackselector::HasTOF,
                  singletrackselector::Px,
                  singletrackselector::Py,
                  singletrackselector::Pz,
                  singletrackselector::TPCInnerParam,
                  singletrackselector::TPCSignal,
                  singletrackselector::Beta,
                  singletrackselector::DcaXY,
                  singletrackselector::DcaZ,
                  singletrackselector::TPCNClsFound,
                  singletrackselector::TPCCrossedRowsOverFindableCls,
                  singletrackselector::TPCChi2NCl,
                  singletrackselector::ITSNCls,
                  singletrackselector::ITSChi2NCl,
                  singletrackselector::Sign,
                  singletrackselector::Eta,
                  singletrackselector::Phi,
                  singletrackselector::StoredCrossedRows,
                  singletrackselector::StoredTOFNSigmaPr,
                  singletrackselector::StoredTPCNSigmaPr,
                  singletrackselector::StoredTOFNSigmaDe,
                  singletrackselector::StoredTPCNSigmaDe,
                  singletrackselector::P<singletrackselector::Px, singletrackselector::Py, singletrackselector::Pz>,
                  singletrackselector::Pt<singletrackselector::Px, singletrackselector::Py>,
                  singletrackselector::CrossedRows<singletrackselector::StoredCrossedRows>,
                  singletrackselector::TOFNSigmaPr<singletrackselector::StoredTOFNSigmaPr>,
                  singletrackselector::TPCNSigmaPr<singletrackselector::StoredTPCNSigmaPr>,
                  singletrackselector::TOFNSigmaDe<singletrackselector::StoredTOFNSigmaDe>,
                  singletrackselector::TPCNSigmaDe<singletrackselector::StoredTPCNSigmaDe>,
                  singletrackselector::Energy<singletrackselector::Px, singletrackselector::Py, singletrackselector::Pz>);

DECLARE_SOA_TABLE(SingleCollSel, "AOD", "SCSEL", // Table of the variables for single track selection.
                  singletrackselector::GlobalIndex,
                  singletrackselector::Mult,
                  singletrackselector::PosZ);
} // namespace o2::aod
#endif // PWGCF_DATAMODEL_SINGLETRACKSELECTOR_H_
