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
#include "Common/DataModel/PIDResponseITS.h"
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

// using nsigma_v0 = binningParent<std::pair<float, float>(-10.f, 10.f)>;
using nsigma_v1 = binningParent<std::pair<float, float>(-6.35f, 6.35f)>; // Width 0.05 symmetric around 0
using nsigma = nsigma_v1;

// using dca_v0 = binningParent<std::pair<float, float>(-1.f, 1.f)>;
// using dca_v1 = binningParent<std::pair<float, float>(-1.f, 1.f), int16_t>;
using dca_v2 = binningParent<std::pair<float, float>(-3.2767f, 3.2767f), int16_t>; // Width 0.0001 symmetric around 0
using dca = dca_v2;

using chi2 = binningParent<std::pair<float, float>(0.f, 10.f)>;
using rowsOverFindable = binningParent<std::pair<float, float>(0.f, 3.f)>;

} // namespace binning

//==================================== base event characteristics ====================================
DECLARE_SOA_COLUMN(Mult, mult, int);                 // Multiplicity of the collision
DECLARE_SOA_COLUMN(MultPercentile, multPerc, float); // Percentiles of multiplicity of the collision
DECLARE_SOA_COLUMN(PosZ, posZ, float);               // Vertex of the collision
DECLARE_SOA_COLUMN(MagField, magField, float);       // Magnetic field corresponding to a collision (in T)

//==================================== extra event characteristics ====================================
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
//==================================== track characteristics ====================================

DECLARE_SOA_INDEX_COLUMN(SingleCollSel, singleCollSel); // Index to the collision
DECLARE_SOA_COLUMN(P, p, float);                        // Momentum of the track
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(Sign, sign, int8_t);
DECLARE_SOA_COLUMN(TPCNClsFound, tpcNClsFound, int16_t);               // Number of TPC clusters
DECLARE_SOA_COLUMN(TPCNClsShared, tpcNClsShared, uint8_t);             // Number of shared TPC clusters
DECLARE_SOA_DYNAMIC_COLUMN(TPCFractionSharedCls, tpcFractionSharedCls, //! Fraction of shared TPC clusters
                           [](uint8_t tpcNClsShared, int16_t tpcNClsFound) -> float { return (float)tpcNClsShared / (float)tpcNClsFound; });

DECLARE_SOA_COLUMN(ITSclsMap, itsClsMap, uint8_t);
DECLARE_SOA_COLUMN(ITSclusterSizes, itsClusterSizes, uint32_t);
DECLARE_SOA_DYNAMIC_COLUMN(ITSNClsDyn, itsNCls, [](uint32_t itsClusterSizes) -> uint8_t {
  uint8_t itsNcls = 0;
  for (int layer = 0; layer < 7; layer++) {
    if ((itsClusterSizes >> (layer * 4)) & 0xf)
      itsNcls++;
  }
  return itsNcls;
});

DECLARE_SOA_COLUMN(StoredTPCChi2NCl, storedTpcChi2NCl, binning::chi2::binned_t); // TPC chi2
DECLARE_SOA_DYNAMIC_COLUMN(TPCChi2NCl, tpcChi2NCl,
                           [](binning::chi2::binned_t chi2_binned) -> float { return singletrackselector::unPack<binning::chi2>(chi2_binned); });

DECLARE_SOA_COLUMN(StoredITSChi2NCl, storedItsChi2NCl, binning::chi2::binned_t); // ITS chi2
DECLARE_SOA_DYNAMIC_COLUMN(ITSChi2NCl, itsChi2NCl,
                           [](binning::chi2::binned_t chi2_binned) -> float { return singletrackselector::unPack<binning::chi2>(chi2_binned); });

DECLARE_SOA_COLUMN(StoredTPCCrossedRowsOverFindableCls, storedTpcCrossedRowsOverFindableCls, binning::rowsOverFindable::binned_t); // Ratio of found over findable clusters
DECLARE_SOA_DYNAMIC_COLUMN(TPCCrossedRowsOverFindableCls, tpcCrossedRowsOverFindableCls,
                           [](binning::rowsOverFindable::binned_t rowsOverFindable_binned) -> float { return singletrackselector::unPack<binning::rowsOverFindable>(rowsOverFindable_binned); });

//==================================== DCA ====================================

DECLARE_SOA_COLUMN(StoredDcaXY, storedDcaXY, binning::dca::binned_t); // impact parameter of the track with 16 bits (v2, larger range)
DECLARE_SOA_DYNAMIC_COLUMN(DcaXY, dcaXY,
                           [](binning::dca::binned_t dca_binned) -> float { return singletrackselector::unPackSymmetric<binning::dca>(dca_binned); });

DECLARE_SOA_COLUMN(StoredDcaZ, storedDcaZ, binning::dca::binned_t); // impact parameter of the track with 16 bits (v2, larger range)
DECLARE_SOA_DYNAMIC_COLUMN(DcaZ, dcaZ,
                           [](binning::dca::binned_t dca_binned) -> float { return singletrackselector::unPackSymmetric<binning::dca>(dca_binned); });

using StoredDcaXY_v2 = StoredDcaXY; // compatibility with the old tables of version 2 -- to be removed later
using StoredDcaZ_v2 = StoredDcaZ;   // compatibility with the old tables of version 2 -- to be removed later
//==================================== PID ====================================

//------------------------------------ Electrons ------------------------------------

DECLARE_SOA_COLUMN(StoredTOFNSigmaEl, storedTofNSigmaEl, binning::nsigma::binned_t); // (v1) TOF
DECLARE_SOA_DYNAMIC_COLUMN(TOFNSigmaEl, tofNSigmaEl,
                           [](binning::nsigma::binned_t nsigma_binned) -> float { return singletrackselector::unPackSymmetric<binning::nsigma>(nsigma_binned); });

DECLARE_SOA_COLUMN(StoredTPCNSigmaEl, storedTpcNSigmaEl, binning::nsigma::binned_t); // (v1) TPC
DECLARE_SOA_DYNAMIC_COLUMN(TPCNSigmaEl, tpcNSigmaEl,
                           [](binning::nsigma::binned_t nsigma_binned) -> float { return singletrackselector::unPackSymmetric<binning::nsigma>(nsigma_binned); });

//------------------------------------ Pions ------------------------------------
DECLARE_SOA_COLUMN(StoredTOFNSigmaPi, storedTofNSigmaPi, binning::nsigma::binned_t); // (v1) TOF
DECLARE_SOA_DYNAMIC_COLUMN(TOFNSigmaPi, tofNSigmaPi,
                           [](binning::nsigma::binned_t nsigma_binned) -> float { return singletrackselector::unPackSymmetric<binning::nsigma>(nsigma_binned); });

DECLARE_SOA_COLUMN(StoredTPCNSigmaPi, storedTpcNSigmaPi, binning::nsigma::binned_t); // (v1) TPC
DECLARE_SOA_DYNAMIC_COLUMN(TPCNSigmaPi, tpcNSigmaPi,
                           [](binning::nsigma::binned_t nsigma_binned) -> float { return singletrackselector::unPackSymmetric<binning::nsigma>(nsigma_binned); });

//------------------------------------ Kaons ------------------------------------
DECLARE_SOA_COLUMN(StoredTOFNSigmaKa, storedTofNSigmaKa, binning::nsigma::binned_t); // (v1) TOF
DECLARE_SOA_DYNAMIC_COLUMN(TOFNSigmaKa, tofNSigmaKa,
                           [](binning::nsigma::binned_t nsigma_binned) -> float { return singletrackselector::unPackSymmetric<binning::nsigma>(nsigma_binned); });

DECLARE_SOA_COLUMN(StoredTPCNSigmaKa, storedTpcNSigmaKa, binning::nsigma::binned_t); // (v1) TPC
DECLARE_SOA_DYNAMIC_COLUMN(TPCNSigmaKa, tpcNSigmaKa,
                           [](binning::nsigma::binned_t nsigma_binned) -> float { return singletrackselector::unPackSymmetric<binning::nsigma>(nsigma_binned); });

//------------------------------------ Protons ------------------------------------
DECLARE_SOA_COLUMN(StoredTOFNSigmaPr, storedTofNSigmaPr, binning::nsigma::binned_t); // (v1) TOF
DECLARE_SOA_DYNAMIC_COLUMN(TOFNSigmaPr, tofNSigmaPr,
                           [](binning::nsigma::binned_t nsigma_binned) -> float { return singletrackselector::unPackSymmetric<binning::nsigma>(nsigma_binned); });

DECLARE_SOA_COLUMN(StoredTPCNSigmaPr, storedTpcNSigmaPr, binning::nsigma::binned_t); // (v1) TPC
DECLARE_SOA_DYNAMIC_COLUMN(TPCNSigmaPr, tpcNSigmaPr,
                           [](binning::nsigma::binned_t nsigma_binned) -> float { return singletrackselector::unPackSymmetric<binning::nsigma>(nsigma_binned); });

//------------------------------------ Deutrons ------------------------------------
DECLARE_SOA_COLUMN(StoredTOFNSigmaDe, storedTofNSigmaDe, binning::nsigma::binned_t); // (v1) TOF
DECLARE_SOA_DYNAMIC_COLUMN(TOFNSigmaDe, tofNSigmaDe,
                           [](binning::nsigma::binned_t nsigma_binned) -> float { return singletrackselector::unPackSymmetric<binning::nsigma>(nsigma_binned); });

DECLARE_SOA_COLUMN(StoredTPCNSigmaDe, storedTpcNSigmaDe, binning::nsigma::binned_t); // (v1) TPC
DECLARE_SOA_DYNAMIC_COLUMN(TPCNSigmaDe, tpcNSigmaDe,
                           [](binning::nsigma::binned_t nsigma_binned) -> float { return singletrackselector::unPackSymmetric<binning::nsigma>(nsigma_binned); });

//------------------------------------ Triton ------------------------------------
DECLARE_SOA_COLUMN(StoredTOFNSigmaTr, storedTofNSigmaTr, binning::nsigma::binned_t); // (v1) TOF
DECLARE_SOA_DYNAMIC_COLUMN(TOFNSigmaTr, tofNSigmaTr,
                           [](binning::nsigma::binned_t nsigma_binned) -> float { return singletrackselector::unPackSymmetric<binning::nsigma>(nsigma_binned); });

DECLARE_SOA_COLUMN(StoredTPCNSigmaTr, storedTpcNSigmaTr, binning::nsigma::binned_t); // (v1) TPC
DECLARE_SOA_DYNAMIC_COLUMN(TPCNSigmaTr, tpcNSigmaTr,
                           [](binning::nsigma::binned_t nsigma_binned) -> float { return singletrackselector::unPackSymmetric<binning::nsigma>(nsigma_binned); });

//------------------------------------ Helium3 ------------------------------------
DECLARE_SOA_COLUMN(StoredTOFNSigmaHe, storedTofNSigmaHe, binning::nsigma::binned_t); // (v1) TOF
DECLARE_SOA_DYNAMIC_COLUMN(TOFNSigmaHe, tofNSigmaHe,
                           [](binning::nsigma::binned_t nsigma_binned) -> float { return singletrackselector::unPackSymmetric<binning::nsigma>(nsigma_binned); });

DECLARE_SOA_COLUMN(StoredTPCNSigmaHe, storedTpcNSigmaHe, binning::nsigma::binned_t); // (v1) TPC
DECLARE_SOA_DYNAMIC_COLUMN(TPCNSigmaHe, tpcNSigmaHe,
                           [](binning::nsigma::binned_t nsigma_binned) -> float { return singletrackselector::unPackSymmetric<binning::nsigma>(nsigma_binned); });

using StoredTOFNSigmaPi_v1 = StoredTOFNSigmaPi; // compatibility with the old tables of version 2 -- to be removed later
using StoredTPCNSigmaPi_v1 = StoredTPCNSigmaPi; // compatibility with the old tables of version 2 -- to be removed later

using StoredTOFNSigmaKa_v1 = StoredTOFNSigmaKa; // compatibility with the old tables of version 2 -- to be removed later
using StoredTPCNSigmaKa_v1 = StoredTPCNSigmaKa; // compatibility with the old tables of version 2 -- to be removed later

using StoredTOFNSigmaPr_v1 = StoredTOFNSigmaPr; // compatibility with the old tables of version 2 -- to be removed later
using StoredTPCNSigmaPr_v1 = StoredTPCNSigmaPr; // compatibility with the old tables of version 2 -- to be removed later

using StoredTOFNSigmaDe_v1 = StoredTOFNSigmaDe; // compatibility with the old tables of version 2 -- to be removed later
using StoredTPCNSigmaDe_v1 = StoredTPCNSigmaDe; // compatibility with the old tables of version 2 -- to be removed later

using StoredTOFNSigmaHe_v1 = StoredTOFNSigmaHe; // compatibility with the old tables of version 2 -- to be removed later
using StoredTPCNSigmaHe_v1 = StoredTPCNSigmaHe; // compatibility with the old tables of version 2 -- to be removed later

//==================================== Dynamic cols for kinematics ====================================

DECLARE_SOA_DYNAMIC_COLUMN(Energy, energy, [](float p, float mass) -> float { return sqrt(p * p + mass * mass); });
DECLARE_SOA_DYNAMIC_COLUMN(Rapidity, rapidity, //! Track rapidity, computed under the mass assumption given as input
                           [](float p, float eta, float mass) -> float {
                             const auto pz = p * std::tanh(eta);
                             const auto energy = std::sqrt(p * p + mass * mass);
                             return 0.5f * log((energy + pz) / (energy - pz));
                           });
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

//==================================== EXtra info ====================================

DECLARE_SOA_COLUMN(TPCInnerParam, tpcInnerParam, float); // Momentum at inner wall of the TPC
DECLARE_SOA_COLUMN(TPCSignal, tpcSignal, float);         // dE/dx TPC
DECLARE_SOA_COLUMN(Beta, beta, float);                   // TOF beta

} // namespace singletrackselector

DECLARE_SOA_TABLE(SingleTrackSels_v3, "AOD", "SINGLETRACKSEL", // Table of the variables for single track selection.
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
                  singletrackselector::StoredDcaXY,
                  singletrackselector::StoredDcaZ,
                  singletrackselector::StoredTPCChi2NCl,
                  singletrackselector::StoredITSChi2NCl,
                  singletrackselector::StoredTPCCrossedRowsOverFindableCls,

                  singletrackselector::ITSNClsDyn<singletrackselector::ITSclusterSizes>,
                  track::v001::ITSClsSizeInLayer<singletrackselector::ITSclusterSizes>,
                  singletrackselector::DcaXY<singletrackselector::StoredDcaXY>,
                  singletrackselector::DcaZ<singletrackselector::StoredDcaZ>,
                  singletrackselector::TPCChi2NCl<singletrackselector::StoredTPCChi2NCl>,
                  singletrackselector::ITSChi2NCl<singletrackselector::StoredITSChi2NCl>,
                  singletrackselector::TPCCrossedRowsOverFindableCls<singletrackselector::StoredTPCCrossedRowsOverFindableCls>,
                  singletrackselector::TPCFractionSharedCls<singletrackselector::TPCNClsShared, singletrackselector::TPCNClsFound>,

                  singletrackselector::Energy<singletrackselector::P>,
                  singletrackselector::Rapidity<singletrackselector::P, singletrackselector::Eta>,
                  singletrackselector::Pt<singletrackselector::P, singletrackselector::Eta>,
                  singletrackselector::Px<singletrackselector::P, singletrackselector::Eta, singletrackselector::Phi>,
                  singletrackselector::Py<singletrackselector::P, singletrackselector::Eta, singletrackselector::Phi>,
                  singletrackselector::Pz<singletrackselector::P, singletrackselector::Eta>,
                  singletrackselector::PhiStar<singletrackselector::P, singletrackselector::Eta, singletrackselector::Sign, singletrackselector::Phi>,

                  // PID with ITS (from PIDResponseITS.h)
                  o2::aod::pidits::ITSNSigmaElImp<singletrackselector::ITSclusterSizes, singletrackselector::P, singletrackselector::Eta>,
                  o2::aod::pidits::ITSNSigmaPiImp<singletrackselector::ITSclusterSizes, singletrackselector::P, singletrackselector::Eta>,
                  o2::aod::pidits::ITSNSigmaKaImp<singletrackselector::ITSclusterSizes, singletrackselector::P, singletrackselector::Eta>,
                  o2::aod::pidits::ITSNSigmaPrImp<singletrackselector::ITSclusterSizes, singletrackselector::P, singletrackselector::Eta>,
                  o2::aod::pidits::ITSNSigmaDeImp<singletrackselector::ITSclusterSizes, singletrackselector::P, singletrackselector::Eta>,
                  o2::aod::pidits::ITSNSigmaTrImp<singletrackselector::ITSclusterSizes, singletrackselector::P, singletrackselector::Eta>,
                  o2::aod::pidits::ITSNSigmaHeImp<singletrackselector::ITSclusterSizes, singletrackselector::P, singletrackselector::Eta>);

DECLARE_SOA_TABLE_VERSIONED(SingleTrackSels_v2, "AOD", "SINGLETRACKSEL2", 2, // Table of the variables for single track selection.
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
                            singletrackselector::DcaXY<singletrackselector::StoredDcaXY_v2>,
                            singletrackselector::DcaZ<singletrackselector::StoredDcaZ_v2>,
                            singletrackselector::TPCChi2NCl<singletrackselector::StoredTPCChi2NCl>,
                            singletrackselector::ITSChi2NCl<singletrackselector::StoredITSChi2NCl>,
                            singletrackselector::TPCCrossedRowsOverFindableCls<singletrackselector::StoredTPCCrossedRowsOverFindableCls>,
                            singletrackselector::TPCFractionSharedCls<singletrackselector::TPCNClsShared, singletrackselector::TPCNClsFound>,

                            singletrackselector::TOFNSigmaPi<singletrackselector::StoredTOFNSigmaPi_v1>,
                            singletrackselector::TPCNSigmaPi<singletrackselector::StoredTPCNSigmaPi_v1>,
                            singletrackselector::TOFNSigmaKa<singletrackselector::StoredTOFNSigmaKa_v1>,
                            singletrackselector::TPCNSigmaKa<singletrackselector::StoredTPCNSigmaKa_v1>,
                            singletrackselector::TOFNSigmaPr<singletrackselector::StoredTOFNSigmaPr_v1>,
                            singletrackselector::TPCNSigmaPr<singletrackselector::StoredTPCNSigmaPr_v1>,
                            singletrackselector::TOFNSigmaDe<singletrackselector::StoredTOFNSigmaDe_v1>,
                            singletrackselector::TPCNSigmaDe<singletrackselector::StoredTPCNSigmaDe_v1>,
                            singletrackselector::TOFNSigmaHe<singletrackselector::StoredTOFNSigmaHe_v1>,
                            singletrackselector::TPCNSigmaHe<singletrackselector::StoredTPCNSigmaHe_v1>,

                            singletrackselector::Rapidity<singletrackselector::P, singletrackselector::Eta>,
                            singletrackselector::Energy<singletrackselector::P>,
                            singletrackselector::Pt<singletrackselector::P, singletrackselector::Eta>,
                            singletrackselector::Px<singletrackselector::P, singletrackselector::Eta, singletrackselector::Phi>,
                            singletrackselector::Py<singletrackselector::P, singletrackselector::Eta, singletrackselector::Phi>,
                            singletrackselector::Pz<singletrackselector::P, singletrackselector::Eta>,
                            singletrackselector::PhiStar<singletrackselector::P, singletrackselector::Eta, singletrackselector::Sign, singletrackselector::Phi>);

using SingleTrackSels = SingleTrackSels_v3;

DECLARE_SOA_TABLE(SinglePIDEls_v0, "AOD", "SINGLEPIDEL0",
                  singletrackselector::StoredTPCNSigmaEl,
                  singletrackselector::TPCNSigmaEl<singletrackselector::StoredTPCNSigmaEl>);

DECLARE_SOA_TABLE(SinglePIDEls, "AOD", "SINGLEPIDEL",
                  singletrackselector::StoredTOFNSigmaEl,
                  singletrackselector::StoredTPCNSigmaEl,

                  singletrackselector::TOFNSigmaEl<singletrackselector::StoredTOFNSigmaEl>,
                  singletrackselector::TPCNSigmaEl<singletrackselector::StoredTPCNSigmaEl>);

DECLARE_SOA_TABLE(SinglePIDPis, "AOD", "SINGLEPIDPI",
                  singletrackselector::StoredTOFNSigmaPi,
                  singletrackselector::StoredTPCNSigmaPi,

                  singletrackselector::TOFNSigmaPi<singletrackselector::StoredTOFNSigmaPi>,
                  singletrackselector::TPCNSigmaPi<singletrackselector::StoredTPCNSigmaPi>);

DECLARE_SOA_TABLE(SinglePIDKas, "AOD", "SINGLEPIDKA",
                  singletrackselector::StoredTOFNSigmaKa,
                  singletrackselector::StoredTPCNSigmaKa,

                  singletrackselector::TOFNSigmaKa<singletrackselector::StoredTOFNSigmaKa>,
                  singletrackselector::TPCNSigmaKa<singletrackselector::StoredTPCNSigmaKa>);

DECLARE_SOA_TABLE(SinglePIDPrs, "AOD", "SINGLEPIDPR",
                  singletrackselector::StoredTOFNSigmaPr,
                  singletrackselector::StoredTPCNSigmaPr,

                  singletrackselector::TOFNSigmaPr<singletrackselector::StoredTOFNSigmaPr>,
                  singletrackselector::TPCNSigmaPr<singletrackselector::StoredTPCNSigmaPr>);

DECLARE_SOA_TABLE(SinglePIDDes, "AOD", "SINGLEPIDDE",
                  singletrackselector::StoredTOFNSigmaDe,
                  singletrackselector::StoredTPCNSigmaDe,

                  singletrackselector::TOFNSigmaDe<singletrackselector::StoredTOFNSigmaDe>,
                  singletrackselector::TPCNSigmaDe<singletrackselector::StoredTPCNSigmaDe>);

DECLARE_SOA_TABLE(SinglePIDTrs, "AOD", "SINGLEPIDTR",
                  singletrackselector::StoredTOFNSigmaTr,
                  singletrackselector::StoredTPCNSigmaTr,

                  singletrackselector::TOFNSigmaTr<singletrackselector::StoredTOFNSigmaTr>,
                  singletrackselector::TPCNSigmaTr<singletrackselector::StoredTPCNSigmaTr>);

DECLARE_SOA_TABLE(SinglePIDHes, "AOD", "SINGLEPIDHE",
                  singletrackselector::StoredTOFNSigmaHe,
                  singletrackselector::StoredTPCNSigmaHe,

                  singletrackselector::TOFNSigmaHe<singletrackselector::StoredTOFNSigmaHe>,
                  singletrackselector::TPCNSigmaHe<singletrackselector::StoredTPCNSigmaHe>);

DECLARE_SOA_TABLE(SingleTrkExtras, "AOD", "SINGLETRKEXTRA",
                  singletrackselector::TPCInnerParam,
                  singletrackselector::TPCSignal,
                  singletrackselector::Beta);

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

// DECLARE_SOA_COLUMN(Vx_MC, vx_MC, float);
// DECLARE_SOA_COLUMN(Vy_MC, vy_MC, float);
// DECLARE_SOA_COLUMN(Vz_MC, vz_MC, float);

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

// DECLARE_SOA_TABLE(SingleTrkMCExtras, "AOD", "SINGLETRKMCEX", // Table with generatad info from MC
//                   singletrackselector::Vx_MC,
//                   singletrackselector::Vy_MC,
//                   singletrackselector::Vz_MC);

} // namespace o2::aod

#endif // PWGCF_FEMTO3D_DATAMODEL_SINGLETRACKSELECTOR_H_

namespace o2::aod::singletrackselector
{
template <typename TrackType>
inline bool ITSselection(TrackType const& track, std::pair<int, std::vector<float>> const& PIDcuts)
{
  int PDG = PIDcuts.first;

  float Nsigma = -1000;
  switch (PDG) {
    case 2212:
      Nsigma = track.itsNSigmaPr();
      break;
    case 1000010020:
      Nsigma = track.itsNSigmaDe();
      break;
    case 1000020030:
      Nsigma = track.itsNSigmaHe();
      break;
    case 1000010030:
      Nsigma = track.itsNSigmaTr();
      break;
    case 211:
      Nsigma = track.itsNSigmaPi();
      break;
    case 321:
      Nsigma = track.itsNSigmaKa();
      break;
    case 0:
      return false;
    default:
      LOG(fatal) << "Cannot interpret PDG for ITS selection: " << PIDcuts.first;
  }

  if (Nsigma > PIDcuts.second[0] && Nsigma < PIDcuts.second[1]) {
    return true;
  }
  return false;
}

template <bool useITS, typename TrackType>
inline bool TPCselection(TrackType const& track, std::pair<int, std::vector<float>> const& PIDcuts, std::vector<float> const& ITSCut = std::vector<float>{})
{
  int PDG = PIDcuts.first;

  if constexpr (useITS) {
    if (ITSCut.size() != 0 && !ITSselection(track, std::make_pair(PDG, ITSCut)))
      return false;
  }

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
    case 1000010030:
      Nsigma = track.tpcNSigmaTr();
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
  if (!TPCselection<false>(track, std::make_pair(PDG, TPCresidualCut)))
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
    case 1000010030:
      Nsigma = track.tofNSigmaTr();
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
