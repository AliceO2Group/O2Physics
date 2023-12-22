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

namespace binning
{
struct binningParent {
 public:
  typedef int8_t binned_t;
  static constexpr int nbins = (1 << 8 * sizeof(binned_t)) - 2;
  static constexpr binned_t overflowBin = nbins >> 1;
  static constexpr binned_t underflowBin = -(nbins >> 1);
};

struct nsigma : binningParent {
 public:
  static constexpr float binned_min = -10.0;
  static constexpr float binned_max = 10.0;
  static constexpr float binned_center = 0.5 * (binned_min + binned_max);
  static constexpr float bin_width = (binned_max - binned_min) / nbins;
};

struct dca : binningParent {
 public:
  static constexpr float binned_min = -1.0;
  static constexpr float binned_max = 1.0;
  static constexpr float binned_center = 0.5 * (binned_min + binned_max);
  static constexpr float bin_width = (binned_max - binned_min) / nbins;
};

struct chi2 : binningParent {
 public:
  static constexpr float binned_min = 0.0;
  static constexpr float binned_max = 10.0;
  static constexpr float binned_center = 0.5 * (binned_min + binned_max);
  static constexpr float bin_width = (binned_max - binned_min) / nbins;
};

struct rowsOverFindable : binningParent {
 public:
  static constexpr float binned_min = 0.0;
  static constexpr float binned_max = 3.0;
  static constexpr float binned_center = 0.5 * (binned_min + binned_max);
  static constexpr float bin_width = (binned_max - binned_min) / nbins;
};

} // namespace binning

DECLARE_SOA_COLUMN(Mult, mult, int);                 // Multiplicity of the collision
DECLARE_SOA_COLUMN(MultPercentile, multPerc, float); // Percentiles of multiplicity of the collision
DECLARE_SOA_COLUMN(PosZ, posZ, float);               // Vertex of the collision
DECLARE_SOA_COLUMN(MagField, magField, float);       // Magnetic field corresponding to a collision (in T)

} // namespace singletrackselector

DECLARE_SOA_TABLE(SingleCollSels, "AOD", "SINGLECOLLSEL", // Table of the variables for single track selection.
                  o2::soa::Index<>,
                  singletrackselector::Mult,
                  singletrackselector::MultPercentile,
                  singletrackselector::PosZ,
                  singletrackselector::MagField);

namespace singletrackselector
{
DECLARE_SOA_INDEX_COLUMN(SingleCollSel, singleCollSel); // Index to the collision
DECLARE_SOA_COLUMN(P, p, float);                        // Momentum of the track
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(Sign, sign, int8_t);
DECLARE_SOA_COLUMN(TPCNClsFound, tpcNClsFound, int16_t);   // Number of TPC clusters
DECLARE_SOA_COLUMN(TPCNClsShared, tpcNClsShared, uint8_t); // Number of shared TPC clusters
DECLARE_SOA_COLUMN(ITSNCls, itsNCls, uint8_t);             // Number of ITS clusters

DECLARE_SOA_COLUMN(StoredDcaXY, storedDcaXY, binning::dca::binned_t);                                                              // impact parameter of the track
DECLARE_SOA_COLUMN(StoredDcaZ, storedDcaZ, binning::dca::binned_t);                                                                // impact parameter of the track
DECLARE_SOA_COLUMN(StoredTPCChi2NCl, storedTpcChi2NCl, binning::chi2::binned_t);                                                   // TPC chi2
DECLARE_SOA_COLUMN(StoredITSChi2NCl, storedItsChi2NCl, binning::chi2::binned_t);                                                   // ITS chi2
DECLARE_SOA_COLUMN(StoredTPCCrossedRowsOverFindableCls, storedTpcCrossedRowsOverFindableCls, binning::rowsOverFindable::binned_t); // Ratio of found over findable clusters

DECLARE_SOA_COLUMN(StoredTOFNSigmaPi, storedTofNSigmaPi, binning::nsigma::binned_t);
DECLARE_SOA_COLUMN(StoredTPCNSigmaPi, storedTpcNSigmaPi, binning::nsigma::binned_t);
DECLARE_SOA_COLUMN(StoredTOFNSigmaKa, storedTofNSigmaKa, binning::nsigma::binned_t);
DECLARE_SOA_COLUMN(StoredTPCNSigmaKa, storedTpcNSigmaKa, binning::nsigma::binned_t);
DECLARE_SOA_COLUMN(StoredTOFNSigmaPr, storedTofNSigmaPr, binning::nsigma::binned_t);
DECLARE_SOA_COLUMN(StoredTPCNSigmaPr, storedTpcNSigmaPr, binning::nsigma::binned_t);
DECLARE_SOA_COLUMN(StoredTOFNSigmaDe, storedTofNSigmaDe, binning::nsigma::binned_t);
DECLARE_SOA_COLUMN(StoredTPCNSigmaDe, storedTpcNSigmaDe, binning::nsigma::binned_t);

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

DECLARE_SOA_DYNAMIC_COLUMN(DcaXY, dcaXY,
                           [](binning::dca::binned_t dca_binned) -> float { return singletrackselector::unPack<binning::dca>(dca_binned); });
DECLARE_SOA_DYNAMIC_COLUMN(DcaZ, dcaZ,
                           [](binning::dca::binned_t dca_binned) -> float { return singletrackselector::unPack<binning::dca>(dca_binned); });
DECLARE_SOA_DYNAMIC_COLUMN(TPCChi2NCl, tpcChi2NCl,
                           [](binning::chi2::binned_t chi2_binned) -> float { return singletrackselector::unPack<binning::chi2>(chi2_binned); });
DECLARE_SOA_DYNAMIC_COLUMN(ITSChi2NCl, itsChi2NCl,
                           [](binning::chi2::binned_t chi2_binned) -> float { return singletrackselector::unPack<binning::chi2>(chi2_binned); });

DECLARE_SOA_DYNAMIC_COLUMN(TPCCrossedRowsOverFindableCls, tpcCrossedRowsOverFindableCls,
                           [](binning::rowsOverFindable::binned_t rowsOverFindable_binned) -> float { return singletrackselector::unPack<binning::rowsOverFindable>(rowsOverFindable_binned); });

DECLARE_SOA_DYNAMIC_COLUMN(TOFNSigmaPi, tofNSigmaPi,
                           [](binning::nsigma::binned_t nsigma_binned) -> float { return singletrackselector::unPack<binning::nsigma>(nsigma_binned); });
DECLARE_SOA_DYNAMIC_COLUMN(TPCNSigmaPi, tpcNSigmaPi,
                           [](binning::nsigma::binned_t nsigma_binned) -> float { return singletrackselector::unPack<binning::nsigma>(nsigma_binned); });
DECLARE_SOA_DYNAMIC_COLUMN(TOFNSigmaKa, tofNSigmaKa,
                           [](binning::nsigma::binned_t nsigma_binned) -> float { return singletrackselector::unPack<binning::nsigma>(nsigma_binned); });
DECLARE_SOA_DYNAMIC_COLUMN(TPCNSigmaKa, tpcNSigmaKa,
                           [](binning::nsigma::binned_t nsigma_binned) -> float { return singletrackselector::unPack<binning::nsigma>(nsigma_binned); });

DECLARE_SOA_DYNAMIC_COLUMN(TOFNSigmaPr, tofNSigmaPr,
                           [](binning::nsigma::binned_t nsigma_binned) -> float { return singletrackselector::unPack<binning::nsigma>(nsigma_binned); });
DECLARE_SOA_DYNAMIC_COLUMN(TPCNSigmaPr, tpcNSigmaPr,
                           [](binning::nsigma::binned_t nsigma_binned) -> float { return singletrackselector::unPack<binning::nsigma>(nsigma_binned); });
DECLARE_SOA_DYNAMIC_COLUMN(TOFNSigmaDe, tofNSigmaDe,
                           [](binning::nsigma::binned_t nsigma_binned) -> float { return singletrackselector::unPack<binning::nsigma>(nsigma_binned); });
DECLARE_SOA_DYNAMIC_COLUMN(TPCNSigmaDe, tpcNSigmaDe,
                           [](binning::nsigma::binned_t nsigma_binned) -> float { return singletrackselector::unPack<binning::nsigma>(nsigma_binned); });

DECLARE_SOA_COLUMN(TPCInnerParam, tpcInnerParam, float); // Momentum at inner wall of the TPC
DECLARE_SOA_COLUMN(TPCSignal, tpcSignal, float);         // dE/dx TPC
DECLARE_SOA_COLUMN(Beta, beta, float);                   // TOF beta

} // namespace singletrackselector

DECLARE_SOA_TABLE_FULL(SingleTrackSels, "SelTracks", "AOD", "SINGLETRACKSEL", // Table of the variables for single track selection.
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

                       singletrackselector::DcaXY<singletrackselector::StoredDcaXY>,
                       singletrackselector::DcaZ<singletrackselector::StoredDcaZ>,
                       singletrackselector::TPCChi2NCl<singletrackselector::StoredTPCChi2NCl>,
                       singletrackselector::ITSChi2NCl<singletrackselector::StoredITSChi2NCl>,
                       singletrackselector::TPCCrossedRowsOverFindableCls<singletrackselector::StoredTPCCrossedRowsOverFindableCls>,

                       singletrackselector::TOFNSigmaPi<singletrackselector::StoredTOFNSigmaPi>,
                       singletrackselector::TPCNSigmaPi<singletrackselector::StoredTPCNSigmaPi>,
                       singletrackselector::TOFNSigmaKa<singletrackselector::StoredTOFNSigmaKa>,
                       singletrackselector::TPCNSigmaKa<singletrackselector::StoredTPCNSigmaKa>,
                       singletrackselector::TOFNSigmaPr<singletrackselector::StoredTOFNSigmaPr>,
                       singletrackselector::TPCNSigmaPr<singletrackselector::StoredTPCNSigmaPr>,
                       singletrackselector::TOFNSigmaDe<singletrackselector::StoredTOFNSigmaDe>,
                       singletrackselector::TPCNSigmaDe<singletrackselector::StoredTPCNSigmaDe>,

                       singletrackselector::Energy<singletrackselector::P>,
                       singletrackselector::Pt<singletrackselector::P, singletrackselector::Eta>,
                       singletrackselector::Px<singletrackselector::P, singletrackselector::Eta, singletrackselector::Phi>,
                       singletrackselector::Py<singletrackselector::P, singletrackselector::Eta, singletrackselector::Phi>,
                       singletrackselector::Pz<singletrackselector::P, singletrackselector::Eta>,
                       singletrackselector::PhiStar<singletrackselector::P, singletrackselector::Eta, singletrackselector::Sign, singletrackselector::Phi>);

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
inline bool TOFselection(TrackType const& track, std::pair<int, std::vector<float>> const& PIDcuts, float const& TPCresidualCut = 5.0f)
{
  int PDG = PIDcuts.first;
  if (!TPCselection(track, std::make_pair(PDG, std::vector<float>{-TPCresidualCut, TPCresidualCut})))
    return false;

  float Nsigma = -1000;
  switch (PDG) {
    case 2212:
      Nsigma = track.tofNSigmaPr();
      break;
    case 1000010020:
      Nsigma = track.tofNSigmaDe();
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
