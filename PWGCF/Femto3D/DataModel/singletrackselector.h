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
/// \author Sofia Tomassini, Gleb Romanenko, Nicolò Jacazio
/// \since 30 May 2023

#ifndef PWGCF_FEMTO3D_DATAMODEL_SINGLETRACKSELECTOR_H_
#define PWGCF_FEMTO3D_DATAMODEL_SINGLETRACKSELECTOR_H_

#include <experimental/type_traits>
#include <utility>
#include <vector>

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
inline typename binningType::binned_t packInTable(const float& valueToBin)
{
  if (valueToBin <= binningType::binned_min) {
    return binningType::underflowBin;
  } else if (valueToBin >= binningType::binned_max) {
    return binningType::overflowBin;
  } else {
    return static_cast<typename binningType::binned_t>((valueToBin - binningType::binned_center) / binningType::bin_width);
  }
}

template <typename binningType>
inline float unPack(const typename binningType::binned_t& b)
{
  return binningType::bin_width * b + binningType::binned_center;
}

template <typename binningType, typename T>
inline o2::framework::expressions::Node unPack(const T& b)
{
  return binningType::bin_width * b + binningType::binned_center;
}

template <typename binningType>
inline typename binningType::binned_t packSymmetric(const float& valueToBin)
{
  if (valueToBin <= binningType::binned_min) {
    return (binningType::underflowBin);
  } else if (valueToBin >= binningType::binned_max) {
    return (binningType::overflowBin);
  } else if (valueToBin >= 0) {
    return (static_cast<typename binningType::binned_t>((valueToBin * binningType::inv_bin_width) + 0.5f));
  } else {
    return (static_cast<typename binningType::binned_t>((valueToBin * binningType::inv_bin_width) - 0.5f));
  }
}

template <typename binningType>
inline float unPackSymmetric(const typename binningType::binned_t& b)
{
  return binningType::bin_width * static_cast<float>(b);
}

template <typename binningType, typename T>
inline o2::framework::expressions::Node unPackSymmetric(const T& b)
{
  return binningType::bin_width * static_cast<float>(b);
}

namespace binning
{

template <std::pair<float, float> lim, typename binVariable = int8_t>
struct binningParent {
 public:
  typedef binVariable binned_t;
  static constexpr int nbins = (1 << 8 * sizeof(binned_t)) - 2;
  static constexpr binned_t overflowBin = nbins >> 1;
  static constexpr binned_t underflowBin = -(nbins >> 1);

  static constexpr float binned_min = lim.first;
  static constexpr float binned_max = lim.second;
  static constexpr float binned_center = 0.5 * (binned_min + binned_max);
  static constexpr float bin_width = (binned_max - binned_min) / nbins;
  static constexpr float inv_bin_width = 1. / bin_width;
  static_assert(binned_min < binned_max, "Invalid binning range");
  static void print()
  {
    LOG(info) << "Binning: " << binned_min << " - " << binned_max << " with " << nbins << " bins, width = "
              << bin_width << ". Overflow bin " << static_cast<int>(overflowBin) << " Underflow bin " << static_cast<int>(underflowBin);
  }
};

using nsigma_v0 = binningParent<std::pair<float, float>(-10.f, 10.f)>;
using nsigma_v1 = binningParent<std::pair<float, float>(-6.35f, 6.35f)>; // Width 0.05 symmetric around 0
using nsigma = nsigma_v1;

using dca_v0 = binningParent<std::pair<float, float>(-1.f, 1.f)>;
using dca_v1 = binningParent<std::pair<float, float>(-1.f, 1.f), int16_t>;
using dca_v2 = binningParent<std::pair<float, float>(-3.2767f, 3.2767f), int16_t>; // Width 0.0001 symmetric around 0
using dca = dca_v2;

using chi2 = binningParent<std::pair<float, float>(0.f, 10.f)>;
using rowsOverFindable = binningParent<std::pair<float, float>(0.f, 3.f)>;

} // namespace binning

DECLARE_SOA_COLUMN(Mult, mult, int);                 // Multiplicity of the collision
DECLARE_SOA_COLUMN(MultPercentile, multPerc, float); // Percentiles of multiplicity of the collision
DECLARE_SOA_COLUMN(PosZ, posZ, float);               // Vertex of the collision
DECLARE_SOA_COLUMN(MagField, magField, float);       // Magnetic field corresponding to a collision (in T)

DECLARE_SOA_COLUMN(IsNoSameBunchPileup, isNoSameBunchPileup, bool);
DECLARE_SOA_COLUMN(IsGoodZvtxFT0vsPV, isGoodZvtxFT0vsPV, bool);
DECLARE_SOA_COLUMN(IsVertexITSTPC, isVertexITSTPC, bool);
DECLARE_SOA_COLUMN(HadronicRate, hadronicRate, double);
DECLARE_SOA_COLUMN(Occupancy, occupancy, int);

DECLARE_SOA_DYNAMIC_COLUMN(dIsNoSameBunchPileup, isNoSameBunchPileup, [](uint64_t selBit) -> bool { return TESTBIT(selBit, evsel::kNoSameBunchPileup); });
DECLARE_SOA_DYNAMIC_COLUMN(dIsGoodZvtxFT0vsPV, isGoodZvtxFT0vsPV, [](uint64_t selBit) -> bool { return TESTBIT(selBit, evsel::kIsGoodZvtxFT0vsPV); });
DECLARE_SOA_DYNAMIC_COLUMN(dIsVertexITSTPC, isVertexITSTPC, [](uint64_t selBit) -> bool { return TESTBIT(selBit, evsel::kIsVertexITSTPC); });
DECLARE_SOA_DYNAMIC_COLUMN(dIsVertexTOForTRDmatched, isVertexTOForTRDmatched, [](uint64_t selBit) -> int { return static_cast<int>(TESTBIT(selBit, evsel::kIsVertexTOFmatched)) + static_cast<int>(TESTBIT(selBit, evsel::kIsVertexTRDmatched)); });
DECLARE_SOA_DYNAMIC_COLUMN(dNoCollInTimeRangeStandard, noCollInTimeRangeStandard, [](uint64_t selBit) -> bool { return TESTBIT(selBit, evsel::kNoCollInTimeRangeStandard); });

} // namespace singletrackselector

DECLARE_SOA_TABLE(SingleCollSels, "AOD", "SINGLECOLLSEL", // Table of the variables for single track selection.
                  o2::soa::Index<>,
                  singletrackselector::Mult,
                  singletrackselector::MultPercentile,
                  singletrackselector::PosZ,
                  singletrackselector::MagField);

DECLARE_SOA_TABLE(SingleCollExtras_v0, "AOD", "SINGLECOLLEXTRA", // Joinable collision table with Pile-Up flags
                  singletrackselector::IsNoSameBunchPileup,
                  singletrackselector::IsGoodZvtxFT0vsPV,
                  singletrackselector::IsVertexITSTPC,
                  singletrackselector::HadronicRate);

DECLARE_SOA_TABLE(SingleCollExtras_v1, "AOD", "SINGLECOLLEXTR1", // Joinable collision table with Pile-Up flags
                  evsel::Selection,
                  singletrackselector::HadronicRate,
                  singletrackselector::Occupancy,

                  singletrackselector::dIsNoSameBunchPileup<evsel::Selection>,
                  singletrackselector::dIsGoodZvtxFT0vsPV<evsel::Selection>,
                  singletrackselector::dIsVertexITSTPC<evsel::Selection>,
                  singletrackselector::dIsVertexTOForTRDmatched<evsel::Selection>,
                  singletrackselector::dNoCollInTimeRangeStandard<evsel::Selection>);

using SingleCollExtras = SingleCollExtras_v1;

namespace singletrackselector
{
DECLARE_SOA_INDEX_COLUMN(SingleCollSel, singleCollSel); // Index to the collision
DECLARE_SOA_COLUMN(P, p, float);                        // Momentum of the track
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(Sign, sign, int8_t);
DECLARE_SOA_COLUMN(TPCNClsFound, tpcNClsFound, int16_t);   // Number of TPC clusters
DECLARE_SOA_COLUMN(TPCNClsShared, tpcNClsShared, uint8_t); // Number of shared TPC clusters
DECLARE_SOA_COLUMN(ITSNCls, itsNCls, uint8_t);             // Number of ITS clusters (only stored in v0)
DECLARE_SOA_DYNAMIC_COLUMN(ITSNClsDyn, itsNCls, [](uint32_t itsClusterSizes) -> uint8_t {
  uint8_t itsNcls = 0;
  for (int layer = 0; layer < 7; layer++) {
    if ((itsClusterSizes >> (layer * 4)) & 0xf)
      itsNcls++;
  }
  return itsNcls;
});
DECLARE_SOA_COLUMN(ITSclsMap, itsClsMap, uint8_t);
DECLARE_SOA_COLUMN(ITSclusterSizes, itsClusterSizes, uint32_t);

DECLARE_SOA_COLUMN(StoredDcaXY, storedDcaXY, binning::dca_v0::binned_t);                                                           // impact parameter of the track with 8 bits (v0)
DECLARE_SOA_COLUMN(StoredDcaZ, storedDcaZ, binning::dca_v0::binned_t);                                                             // impact parameter of the track with 8 bits (v0)
DECLARE_SOA_COLUMN(StoredDcaXY_v1, storedDcaXY_v1, binning::dca_v1::binned_t);                                                     // impact parameter of the track with 16 bits (v1)
DECLARE_SOA_COLUMN(StoredDcaZ_v1, storedDcaZ_v1, binning::dca_v1::binned_t);                                                       // impact parameter of the track with 16 bits (v1)
DECLARE_SOA_COLUMN(StoredDcaXY_v2, storedDcaXY_v2, binning::dca_v2::binned_t);                                                     // impact parameter of the track with 16 bits (v2, larger range)
DECLARE_SOA_COLUMN(StoredDcaZ_v2, storedDcaZ_v2, binning::dca_v2::binned_t);                                                       // impact parameter of the track with 16 bits (v2, larger range)
DECLARE_SOA_COLUMN(StoredTPCChi2NCl, storedTpcChi2NCl, binning::chi2::binned_t);                                                   // TPC chi2
DECLARE_SOA_COLUMN(StoredITSChi2NCl, storedItsChi2NCl, binning::chi2::binned_t);                                                   // ITS chi2
DECLARE_SOA_COLUMN(StoredTPCCrossedRowsOverFindableCls, storedTpcCrossedRowsOverFindableCls, binning::rowsOverFindable::binned_t); // Ratio of found over findable clusters

DECLARE_SOA_COLUMN(StoredTOFNSigmaPi, storedTofNSigmaPi, binning::nsigma::binned_t); // (v0)
DECLARE_SOA_COLUMN(StoredTPCNSigmaPi, storedTpcNSigmaPi, binning::nsigma::binned_t); // (v0)
DECLARE_SOA_COLUMN(StoredTOFNSigmaKa, storedTofNSigmaKa, binning::nsigma::binned_t); // (v0)
DECLARE_SOA_COLUMN(StoredTPCNSigmaKa, storedTpcNSigmaKa, binning::nsigma::binned_t); // (v0)
DECLARE_SOA_COLUMN(StoredTOFNSigmaPr, storedTofNSigmaPr, binning::nsigma::binned_t); // (v0)
DECLARE_SOA_COLUMN(StoredTPCNSigmaPr, storedTpcNSigmaPr, binning::nsigma::binned_t); // (v0)
DECLARE_SOA_COLUMN(StoredTOFNSigmaDe, storedTofNSigmaDe, binning::nsigma::binned_t); // (v0)
DECLARE_SOA_COLUMN(StoredTPCNSigmaDe, storedTpcNSigmaDe, binning::nsigma::binned_t); // (v0)
DECLARE_SOA_COLUMN(StoredTOFNSigmaHe, storedTofNSigmaHe, binning::nsigma::binned_t); // (v0)
DECLARE_SOA_COLUMN(StoredTPCNSigmaHe, storedTpcNSigmaHe, binning::nsigma::binned_t); // (v0)

DECLARE_SOA_COLUMN(StoredTOFNSigmaPi_v1, storedTofNSigmaPi_v1, binning::nsigma::binned_t); // (v1)
DECLARE_SOA_COLUMN(StoredTPCNSigmaPi_v1, storedTpcNSigmaPi_v1, binning::nsigma::binned_t); // (v1)
DECLARE_SOA_COLUMN(StoredTOFNSigmaKa_v1, storedTofNSigmaKa_v1, binning::nsigma::binned_t); // (v1)
DECLARE_SOA_COLUMN(StoredTPCNSigmaKa_v1, storedTpcNSigmaKa_v1, binning::nsigma::binned_t); // (v1)
DECLARE_SOA_COLUMN(StoredTOFNSigmaPr_v1, storedTofNSigmaPr_v1, binning::nsigma::binned_t); // (v1)
DECLARE_SOA_COLUMN(StoredTPCNSigmaPr_v1, storedTpcNSigmaPr_v1, binning::nsigma::binned_t); // (v1)
DECLARE_SOA_COLUMN(StoredTOFNSigmaDe_v1, storedTofNSigmaDe_v1, binning::nsigma::binned_t); // (v1)
DECLARE_SOA_COLUMN(StoredTPCNSigmaDe_v1, storedTpcNSigmaDe_v1, binning::nsigma::binned_t); // (v1)
DECLARE_SOA_COLUMN(StoredTOFNSigmaHe_v1, storedTofNSigmaHe_v1, binning::nsigma::binned_t); // (v1)
DECLARE_SOA_COLUMN(StoredTPCNSigmaHe_v1, storedTpcNSigmaHe_v1, binning::nsigma::binned_t); // (v1)

DECLARE_SOA_DYNAMIC_COLUMN(Energy, energy, [](float p, float mass) -> float { return sqrt(p * p + mass * mass); });
DECLARE_SOA_DYNAMIC_COLUMN(Pt, pt, [](float p, float eta) -> float { return p / std::cosh(eta); });
DECLARE_SOA_DYNAMIC_COLUMN(Px, px, [](float p, float eta, float phi) -> float { return (p / std::cosh(eta)) * std::sin(phi); });
DECLARE_SOA_DYNAMIC_COLUMN(Py, py, [](float p, float eta, float phi) -> float { return (p / std::cosh(eta)) * std::cos(phi); });
DECLARE_SOA_DYNAMIC_COLUMN(Pz, pz, [](float p, float eta) -> float { return p * std::tanh(eta); });
DECLARE_SOA_DYNAMIC_COLUMN(PhiStar, phiStar,
                           [](float p, float eta, float sign, float phi, float magfield = 0.0, float radius = 1.6) -> float {
                             if (magfield == 0.0) {
                               return -1000.0;
                             } else {
                               return phi + std::asin(-0.3 * magfield * sign * radius / (2.0 * p / std::cosh(eta)));
                             }
                           });

DECLARE_SOA_DYNAMIC_COLUMN(DcaXY_v0, dcaXY,
                           [](binning::dca_v0::binned_t dca_binned) -> float { return singletrackselector::unPack<binning::dca_v0>(dca_binned); });
DECLARE_SOA_DYNAMIC_COLUMN(DcaZ_v0, dcaZ,
                           [](binning::dca_v0::binned_t dca_binned) -> float { return singletrackselector::unPack<binning::dca_v0>(dca_binned); });
DECLARE_SOA_DYNAMIC_COLUMN(DcaXY_v1, dcaXY,
                           [](binning::dca_v1::binned_t dca_binned) -> float { return singletrackselector::unPack<binning::dca_v1>(dca_binned); });
DECLARE_SOA_DYNAMIC_COLUMN(DcaZ_v1, dcaZ,
                           [](binning::dca_v1::binned_t dca_binned) -> float { return singletrackselector::unPack<binning::dca_v1>(dca_binned); });

DECLARE_SOA_DYNAMIC_COLUMN(DcaXY_v2, dcaXY,
                           [](binning::dca_v2::binned_t dca_binned) -> float { return singletrackselector::unPackSymmetric<binning::dca_v2>(dca_binned); });
DECLARE_SOA_DYNAMIC_COLUMN(DcaZ_v2, dcaZ,
                           [](binning::dca_v2::binned_t dca_binned) -> float { return singletrackselector::unPackSymmetric<binning::dca_v2>(dca_binned); });

DECLARE_SOA_DYNAMIC_COLUMN(TPCChi2NCl, tpcChi2NCl,
                           [](binning::chi2::binned_t chi2_binned) -> float { return singletrackselector::unPack<binning::chi2>(chi2_binned); });
DECLARE_SOA_DYNAMIC_COLUMN(ITSChi2NCl, itsChi2NCl,
                           [](binning::chi2::binned_t chi2_binned) -> float { return singletrackselector::unPack<binning::chi2>(chi2_binned); });

DECLARE_SOA_DYNAMIC_COLUMN(TPCCrossedRowsOverFindableCls, tpcCrossedRowsOverFindableCls,
                           [](binning::rowsOverFindable::binned_t rowsOverFindable_binned) -> float { return singletrackselector::unPack<binning::rowsOverFindable>(rowsOverFindable_binned); });

DECLARE_SOA_DYNAMIC_COLUMN(TPCFractionSharedCls, tpcFractionSharedCls, //! Fraction of shared TPC clusters
                           [](uint8_t tpcNClsShared, int16_t tpcNClsFound) -> float { return (float)tpcNClsShared / (float)tpcNClsFound; });

DECLARE_SOA_DYNAMIC_COLUMN(TOFNSigmaPi_v0, tofNSigmaPi,
                           [](binning::nsigma_v0::binned_t nsigma_binned) -> float { return singletrackselector::unPack<binning::nsigma_v0>(nsigma_binned); });
DECLARE_SOA_DYNAMIC_COLUMN(TPCNSigmaPi_v0, tpcNSigmaPi,
                           [](binning::nsigma_v0::binned_t nsigma_binned) -> float { return singletrackselector::unPack<binning::nsigma_v0>(nsigma_binned); });
DECLARE_SOA_DYNAMIC_COLUMN(TOFNSigmaKa_v0, tofNSigmaKa,
                           [](binning::nsigma_v0::binned_t nsigma_binned) -> float { return singletrackselector::unPack<binning::nsigma_v0>(nsigma_binned); });
DECLARE_SOA_DYNAMIC_COLUMN(TPCNSigmaKa_v0, tpcNSigmaKa,
                           [](binning::nsigma_v0::binned_t nsigma_binned) -> float { return singletrackselector::unPack<binning::nsigma_v0>(nsigma_binned); });
DECLARE_SOA_DYNAMIC_COLUMN(TOFNSigmaPr_v0, tofNSigmaPr,
                           [](binning::nsigma_v0::binned_t nsigma_binned) -> float { return singletrackselector::unPack<binning::nsigma_v0>(nsigma_binned); });
DECLARE_SOA_DYNAMIC_COLUMN(TPCNSigmaPr_v0, tpcNSigmaPr,
                           [](binning::nsigma_v0::binned_t nsigma_binned) -> float { return singletrackselector::unPack<binning::nsigma_v0>(nsigma_binned); });
DECLARE_SOA_DYNAMIC_COLUMN(TOFNSigmaDe_v0, tofNSigmaDe,
                           [](binning::nsigma_v0::binned_t nsigma_binned) -> float { return singletrackselector::unPack<binning::nsigma_v0>(nsigma_binned); });
DECLARE_SOA_DYNAMIC_COLUMN(TPCNSigmaDe_v0, tpcNSigmaDe,
                           [](binning::nsigma_v0::binned_t nsigma_binned) -> float { return singletrackselector::unPack<binning::nsigma_v0>(nsigma_binned); });
DECLARE_SOA_DYNAMIC_COLUMN(TOFNSigmaHe_v0, tofNSigmaHe,
                           [](binning::nsigma_v0::binned_t nsigma_binned) -> float { return singletrackselector::unPack<binning::nsigma_v0>(nsigma_binned); });
DECLARE_SOA_DYNAMIC_COLUMN(TPCNSigmaHe_v0, tpcNSigmaHe,
                           [](binning::nsigma_v0::binned_t nsigma_binned) -> float { return singletrackselector::unPack<binning::nsigma_v0>(nsigma_binned); });

DECLARE_SOA_DYNAMIC_COLUMN(TOFNSigmaPi_v1, tofNSigmaPi,
                           [](binning::nsigma_v1::binned_t nsigma_binned) -> float { return singletrackselector::unPackSymmetric<binning::nsigma_v1>(nsigma_binned); });
DECLARE_SOA_DYNAMIC_COLUMN(TPCNSigmaPi_v1, tpcNSigmaPi,
                           [](binning::nsigma_v1::binned_t nsigma_binned) -> float { return singletrackselector::unPackSymmetric<binning::nsigma_v1>(nsigma_binned); });
DECLARE_SOA_DYNAMIC_COLUMN(TOFNSigmaKa_v1, tofNSigmaKa,
                           [](binning::nsigma_v1::binned_t nsigma_binned) -> float { return singletrackselector::unPackSymmetric<binning::nsigma_v1>(nsigma_binned); });
DECLARE_SOA_DYNAMIC_COLUMN(TPCNSigmaKa_v1, tpcNSigmaKa,
                           [](binning::nsigma_v1::binned_t nsigma_binned) -> float { return singletrackselector::unPackSymmetric<binning::nsigma_v1>(nsigma_binned); });
DECLARE_SOA_DYNAMIC_COLUMN(TOFNSigmaPr_v1, tofNSigmaPr,
                           [](binning::nsigma_v1::binned_t nsigma_binned) -> float { return singletrackselector::unPackSymmetric<binning::nsigma_v1>(nsigma_binned); });
DECLARE_SOA_DYNAMIC_COLUMN(TPCNSigmaPr_v1, tpcNSigmaPr,
                           [](binning::nsigma_v1::binned_t nsigma_binned) -> float { return singletrackselector::unPackSymmetric<binning::nsigma_v1>(nsigma_binned); });
DECLARE_SOA_DYNAMIC_COLUMN(TOFNSigmaDe_v1, tofNSigmaDe,
                           [](binning::nsigma_v1::binned_t nsigma_binned) -> float { return singletrackselector::unPackSymmetric<binning::nsigma_v1>(nsigma_binned); });
DECLARE_SOA_DYNAMIC_COLUMN(TPCNSigmaDe_v1, tpcNSigmaDe,
                           [](binning::nsigma_v1::binned_t nsigma_binned) -> float { return singletrackselector::unPackSymmetric<binning::nsigma_v1>(nsigma_binned); });
DECLARE_SOA_DYNAMIC_COLUMN(TOFNSigmaHe_v1, tofNSigmaHe,
                           [](binning::nsigma_v1::binned_t nsigma_binned) -> float { return singletrackselector::unPackSymmetric<binning::nsigma_v1>(nsigma_binned); });
DECLARE_SOA_DYNAMIC_COLUMN(TPCNSigmaHe_v1, tpcNSigmaHe,
                           [](binning::nsigma_v1::binned_t nsigma_binned) -> float { return singletrackselector::unPackSymmetric<binning::nsigma_v1>(nsigma_binned); });

DECLARE_SOA_COLUMN(TPCInnerParam, tpcInnerParam, float); // Momentum at inner wall of the TPC
DECLARE_SOA_COLUMN(TPCSignal, tpcSignal, float);         // dE/dx TPC
DECLARE_SOA_COLUMN(Beta, beta, float);                   // TOF beta

DECLARE_SOA_COLUMN(StoredTPCNSigmaEl, storedTpcNSigmaEl, binning::nsigma::binned_t);
DECLARE_SOA_DYNAMIC_COLUMN(TPCNSigmaEl, tpcNSigmaEl,
                           [](binning::nsigma_v1::binned_t nsigma_binned) -> float { return singletrackselector::unPackSymmetric<binning::nsigma_v1>(nsigma_binned); });

DECLARE_SOA_DYNAMIC_COLUMN(Rapidity, rapidity, //! Track rapidity, computed under the mass assumption given as input
                           [](float p, float eta, float mass) -> float {
                             const auto pz = p * std::tanh(eta);
                             const auto energy = std::sqrt(p * p + mass * mass);
                             return 0.5f * log((energy + pz) / (energy - pz));
                           });

} // namespace singletrackselector

DECLARE_SOA_TABLE_FULL(SingleTrackSels_v0, "SelTracks", "AOD", "SINGLETRACKSEL", // Table of the variables for single track selection.
                       o2::soa::Index<>,
                       singletrackselector::SingleCollSelId,
                       singletrackselector::P,
                       singletrackselector::Eta,
                       singletrackselector::Phi,
                       singletrackselector::Sign,
                       singletrackselector::TPCNClsFound,
                       singletrackselector::TPCNClsShared,
                       singletrackselector::ITSNCls,
                       singletrackselector::StoredDcaXY,
                       singletrackselector::StoredDcaZ,
                       singletrackselector::StoredTPCChi2NCl,
                       singletrackselector::StoredITSChi2NCl,
                       singletrackselector::StoredTPCCrossedRowsOverFindableCls,

                       singletrackselector::StoredTOFNSigmaPi,
                       singletrackselector::StoredTPCNSigmaPi,
                       singletrackselector::StoredTOFNSigmaKa,
                       singletrackselector::StoredTPCNSigmaKa,
                       singletrackselector::StoredTOFNSigmaPr,
                       singletrackselector::StoredTPCNSigmaPr,
                       singletrackselector::StoredTOFNSigmaDe,
                       singletrackselector::StoredTPCNSigmaDe,

                       singletrackselector::DcaXY_v0<singletrackselector::StoredDcaXY>,
                       singletrackselector::DcaZ_v0<singletrackselector::StoredDcaZ>,
                       singletrackselector::TPCChi2NCl<singletrackselector::StoredTPCChi2NCl>,
                       singletrackselector::ITSChi2NCl<singletrackselector::StoredITSChi2NCl>,
                       singletrackselector::TPCCrossedRowsOverFindableCls<singletrackselector::StoredTPCCrossedRowsOverFindableCls>,
                       singletrackselector::TPCFractionSharedCls<singletrackselector::TPCNClsShared, singletrackselector::TPCNClsFound>,

                       singletrackselector::TOFNSigmaPi_v0<singletrackselector::StoredTOFNSigmaPi>,
                       singletrackselector::TPCNSigmaPi_v0<singletrackselector::StoredTPCNSigmaPi>,
                       singletrackselector::TOFNSigmaKa_v0<singletrackselector::StoredTOFNSigmaKa>,
                       singletrackselector::TPCNSigmaKa_v0<singletrackselector::StoredTPCNSigmaKa>,
                       singletrackselector::TOFNSigmaPr_v0<singletrackselector::StoredTOFNSigmaPr>,
                       singletrackselector::TPCNSigmaPr_v0<singletrackselector::StoredTPCNSigmaPr>,
                       singletrackselector::TOFNSigmaDe_v0<singletrackselector::StoredTOFNSigmaDe>,
                       singletrackselector::TPCNSigmaDe_v0<singletrackselector::StoredTPCNSigmaDe>,

                       singletrackselector::Rapidity<singletrackselector::P, singletrackselector::Eta>,
                       singletrackselector::Energy<singletrackselector::P>,
                       singletrackselector::Pt<singletrackselector::P, singletrackselector::Eta>,
                       singletrackselector::Px<singletrackselector::P, singletrackselector::Eta, singletrackselector::Phi>,
                       singletrackselector::Py<singletrackselector::P, singletrackselector::Eta, singletrackselector::Phi>,
                       singletrackselector::Pz<singletrackselector::P, singletrackselector::Eta>,
                       singletrackselector::PhiStar<singletrackselector::P, singletrackselector::Eta, singletrackselector::Sign, singletrackselector::Phi>);

DECLARE_SOA_TABLE_FULL(SingleTrackSels_v1, "SelTracks", "AOD", "SINGLETRACKSEL1", // Table of the variables for single track selection.
                       o2::soa::Index<>,
                       singletrackselector::SingleCollSelId,
                       singletrackselector::P,
                       singletrackselector::Eta,
                       singletrackselector::Phi,
                       singletrackselector::Sign,
                       singletrackselector::TPCNClsFound,
                       singletrackselector::TPCNClsShared,
                       singletrackselector::ITSclsMap,
                       singletrackselector::ITSclusterSizes,
                       singletrackselector::StoredDcaXY_v1,
                       singletrackselector::StoredDcaZ_v1,
                       singletrackselector::StoredTPCChi2NCl,
                       singletrackselector::StoredITSChi2NCl,
                       singletrackselector::StoredTPCCrossedRowsOverFindableCls,

                       singletrackselector::StoredTOFNSigmaPi,
                       singletrackselector::StoredTPCNSigmaPi,
                       singletrackselector::StoredTOFNSigmaKa,
                       singletrackselector::StoredTPCNSigmaKa,
                       singletrackselector::StoredTOFNSigmaPr,
                       singletrackselector::StoredTPCNSigmaPr,
                       singletrackselector::StoredTOFNSigmaDe,
                       singletrackselector::StoredTPCNSigmaDe,
                       singletrackselector::StoredTOFNSigmaHe,
                       singletrackselector::StoredTPCNSigmaHe,

                       singletrackselector::ITSNClsDyn<singletrackselector::ITSclusterSizes>,
                       track::v001::ITSClsSizeInLayer<singletrackselector::ITSclusterSizes>,
                       singletrackselector::DcaXY_v1<singletrackselector::StoredDcaXY_v1>,
                       singletrackselector::DcaZ_v1<singletrackselector::StoredDcaZ_v1>,
                       singletrackselector::TPCChi2NCl<singletrackselector::StoredTPCChi2NCl>,
                       singletrackselector::ITSChi2NCl<singletrackselector::StoredITSChi2NCl>,
                       singletrackselector::TPCCrossedRowsOverFindableCls<singletrackselector::StoredTPCCrossedRowsOverFindableCls>,
                       singletrackselector::TPCFractionSharedCls<singletrackselector::TPCNClsShared, singletrackselector::TPCNClsFound>,

                       singletrackselector::TOFNSigmaPi_v0<singletrackselector::StoredTOFNSigmaPi>,
                       singletrackselector::TPCNSigmaPi_v0<singletrackselector::StoredTPCNSigmaPi>,
                       singletrackselector::TOFNSigmaKa_v0<singletrackselector::StoredTOFNSigmaKa>,
                       singletrackselector::TPCNSigmaKa_v0<singletrackselector::StoredTPCNSigmaKa>,
                       singletrackselector::TOFNSigmaPr_v0<singletrackselector::StoredTOFNSigmaPr>,
                       singletrackselector::TPCNSigmaPr_v0<singletrackselector::StoredTPCNSigmaPr>,
                       singletrackselector::TOFNSigmaDe_v0<singletrackselector::StoredTOFNSigmaDe>,
                       singletrackselector::TPCNSigmaDe_v0<singletrackselector::StoredTPCNSigmaDe>,
                       singletrackselector::TOFNSigmaHe_v0<singletrackselector::StoredTOFNSigmaHe>,
                       singletrackselector::TPCNSigmaHe_v0<singletrackselector::StoredTPCNSigmaHe>,

                       singletrackselector::Rapidity<singletrackselector::P, singletrackselector::Eta>,
                       singletrackselector::Energy<singletrackselector::P>,
                       singletrackselector::Pt<singletrackselector::P, singletrackselector::Eta>,
                       singletrackselector::Px<singletrackselector::P, singletrackselector::Eta, singletrackselector::Phi>,
                       singletrackselector::Py<singletrackselector::P, singletrackselector::Eta, singletrackselector::Phi>,
                       singletrackselector::Pz<singletrackselector::P, singletrackselector::Eta>,
                       singletrackselector::PhiStar<singletrackselector::P, singletrackselector::Eta, singletrackselector::Sign, singletrackselector::Phi>);

DECLARE_SOA_TABLE_FULL(SingleTrackSels_v2, "SelTracks", "AOD", "SINGLETRACKSEL2", // Table of the variables for single track selection.
                       o2::soa::Index<>,
                       singletrackselector::SingleCollSelId,
                       singletrackselector::P,
                       singletrackselector::Eta,
                       singletrackselector::Phi,
                       singletrackselector::Sign,
                       singletrackselector::TPCNClsFound,
                       singletrackselector::TPCNClsShared,
                       singletrackselector::ITSclsMap,
                       singletrackselector::ITSclusterSizes,
                       singletrackselector::StoredDcaXY_v2,
                       singletrackselector::StoredDcaZ_v2,
                       singletrackselector::StoredTPCChi2NCl,
                       singletrackselector::StoredITSChi2NCl,
                       singletrackselector::StoredTPCCrossedRowsOverFindableCls,

                       singletrackselector::StoredTOFNSigmaPi_v1,
                       singletrackselector::StoredTPCNSigmaPi_v1,
                       singletrackselector::StoredTOFNSigmaKa_v1,
                       singletrackselector::StoredTPCNSigmaKa_v1,
                       singletrackselector::StoredTOFNSigmaPr_v1,
                       singletrackselector::StoredTPCNSigmaPr_v1,
                       singletrackselector::StoredTOFNSigmaDe_v1,
                       singletrackselector::StoredTPCNSigmaDe_v1,
                       singletrackselector::StoredTOFNSigmaHe_v1,
                       singletrackselector::StoredTPCNSigmaHe_v1,

                       singletrackselector::ITSNClsDyn<singletrackselector::ITSclusterSizes>,
                       track::v001::ITSClsSizeInLayer<singletrackselector::ITSclusterSizes>,
                       singletrackselector::DcaXY_v2<singletrackselector::StoredDcaXY_v2>,
                       singletrackselector::DcaZ_v2<singletrackselector::StoredDcaZ_v2>,
                       singletrackselector::TPCChi2NCl<singletrackselector::StoredTPCChi2NCl>,
                       singletrackselector::ITSChi2NCl<singletrackselector::StoredITSChi2NCl>,
                       singletrackselector::TPCCrossedRowsOverFindableCls<singletrackselector::StoredTPCCrossedRowsOverFindableCls>,
                       singletrackselector::TPCFractionSharedCls<singletrackselector::TPCNClsShared, singletrackselector::TPCNClsFound>,

                       singletrackselector::TOFNSigmaPi_v1<singletrackselector::StoredTOFNSigmaPi_v1>,
                       singletrackselector::TPCNSigmaPi_v1<singletrackselector::StoredTPCNSigmaPi_v1>,
                       singletrackselector::TOFNSigmaKa_v1<singletrackselector::StoredTOFNSigmaKa_v1>,
                       singletrackselector::TPCNSigmaKa_v1<singletrackselector::StoredTPCNSigmaKa_v1>,
                       singletrackselector::TOFNSigmaPr_v1<singletrackselector::StoredTOFNSigmaPr_v1>,
                       singletrackselector::TPCNSigmaPr_v1<singletrackselector::StoredTPCNSigmaPr_v1>,
                       singletrackselector::TOFNSigmaDe_v1<singletrackselector::StoredTOFNSigmaDe_v1>,
                       singletrackselector::TPCNSigmaDe_v1<singletrackselector::StoredTPCNSigmaDe_v1>,
                       singletrackselector::TOFNSigmaHe_v1<singletrackselector::StoredTOFNSigmaHe_v1>,
                       singletrackselector::TPCNSigmaHe_v1<singletrackselector::StoredTPCNSigmaHe_v1>,

                       singletrackselector::Rapidity<singletrackselector::P, singletrackselector::Eta>,
                       singletrackselector::Energy<singletrackselector::P>,
                       singletrackselector::Pt<singletrackselector::P, singletrackselector::Eta>,
                       singletrackselector::Px<singletrackselector::P, singletrackselector::Eta, singletrackselector::Phi>,
                       singletrackselector::Py<singletrackselector::P, singletrackselector::Eta, singletrackselector::Phi>,
                       singletrackselector::Pz<singletrackselector::P, singletrackselector::Eta>,
                       singletrackselector::PhiStar<singletrackselector::P, singletrackselector::Eta, singletrackselector::Sign, singletrackselector::Phi>);

using SingleTrackSels = SingleTrackSels_v2;

DECLARE_SOA_TABLE(SingleTrkExtras, "AOD", "SINGLETRKEXTRA",
                  singletrackselector::TPCInnerParam,
                  singletrackselector::TPCSignal,
                  singletrackselector::Beta);

DECLARE_SOA_TABLE(SinglePIDEls, "AOD", "SINGLEPIDEL",
                  singletrackselector::StoredTPCNSigmaEl,
                  singletrackselector::TPCNSigmaEl<singletrackselector::StoredTPCNSigmaEl>);

namespace singletrackselector
{
DECLARE_SOA_COLUMN(PdgCode, pdgCode, int);
DECLARE_SOA_COLUMN(Origin, origin, int); // 0 - prymary; 1 - weak decay; 2 - material
DECLARE_SOA_COLUMN(P_MC, p_MC, float);
DECLARE_SOA_COLUMN(Eta_MC, eta_MC, float);
DECLARE_SOA_COLUMN(Phi_MC, phi_MC, float);

DECLARE_SOA_DYNAMIC_COLUMN(Pt_MC, pt_MC, [](float p, float eta) -> float { return p / std::cosh(eta); });
DECLARE_SOA_DYNAMIC_COLUMN(Px_MC, px_MC, [](float p, float eta, float phi) -> float { return (p / std::cosh(eta)) * std::sin(phi); });
DECLARE_SOA_DYNAMIC_COLUMN(Py_MC, py_MC, [](float p, float eta, float phi) -> float { return (p / std::cosh(eta)) * std::cos(phi); });
DECLARE_SOA_DYNAMIC_COLUMN(Pz_MC, pz_MC, [](float p, float eta) -> float { return p * std::tanh(eta); });

} // namespace singletrackselector

DECLARE_SOA_TABLE(SingleTrkMCs, "AOD", "SINGLETRKMC", // Table with generatad info from MC
                  singletrackselector::PdgCode,
                  singletrackselector::Origin,
                  singletrackselector::P_MC,
                  singletrackselector::Eta_MC,
                  singletrackselector::Phi_MC,
                  singletrackselector::Pt_MC<singletrackselector::P_MC, singletrackselector::Eta>,
                  singletrackselector::Px_MC<singletrackselector::P_MC, singletrackselector::Eta_MC, singletrackselector::Phi_MC>,
                  singletrackselector::Py_MC<singletrackselector::P_MC, singletrackselector::Eta_MC, singletrackselector::Phi_MC>,
                  singletrackselector::Pz_MC<singletrackselector::P_MC, singletrackselector::Eta_MC>);

} // namespace o2::aod

#endif // PWGCF_FEMTO3D_DATAMODEL_SINGLETRACKSELECTOR_H_

namespace o2::aod::singletrackselector
{

template <typename TrackType>
inline bool TPCselection(TrackType const& track, std::pair<int, std::vector<float>> const& PIDcuts)
{
  int PDG = PIDcuts.first;
  float Nsigma = -1000;
  switch (PDG) {
    case 2212:
      Nsigma = track.tpcNSigmaPr();
      break;
    case 1000010020:
      Nsigma = track.tpcNSigmaDe();
      break;
    case 1000020030:
      Nsigma = track.tpcNSigmaHe();
      break;
    case 211:
      Nsigma = track.tpcNSigmaPi();
      break;
    case 321:
      Nsigma = track.tpcNSigmaKa();
      break;
    case 0:
      return false;
    default:
      LOG(fatal) << "Cannot interpret PDG for TPC selection: " << PIDcuts.first;
  }

  if (Nsigma > PIDcuts.second[0] && Nsigma < PIDcuts.second[1]) {
    return true;
  }
  return false;
}

template <typename TrackType>
inline bool TOFselection(TrackType const& track, std::pair<int, std::vector<float>> const& PIDcuts, std::vector<float> const& TPCresidualCut = std::vector<float>{-5.0f, 5.0f})
{
  int PDG = PIDcuts.first;
  if (!TPCselection(track, std::make_pair(PDG, TPCresidualCut)))
    return false;

  float Nsigma = -1000;
  switch (PDG) {
    case 2212:
      Nsigma = track.tofNSigmaPr();
      break;
    case 1000010020:
      Nsigma = track.tofNSigmaDe();
      break;
    case 1000020030:
      Nsigma = track.tofNSigmaHe();
      break;
    case 211:
      Nsigma = track.tofNSigmaPi();
      break;
    case 321:
      Nsigma = track.tofNSigmaKa();
      break;
    case 0:
      return false;
    default:
      LOG(fatal) << "Cannot interpret PDG for TOF selection: " << PIDcuts.first;
  }

  if (Nsigma > PIDcuts.second[0] && Nsigma < PIDcuts.second[1]) {
    return true;
  }
  return false;
}

} // namespace o2::aod::singletrackselector
