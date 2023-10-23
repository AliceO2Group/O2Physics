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

#ifndef PWGCF_DATAMODEL_SINGLETRACKSELECTOR_H_
#define PWGCF_DATAMODEL_SINGLETRACKSELECTOR_H_

#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "Common/DataModel/PIDResponse.h"
#include "Framework/Logger.h"
#include "Common/DataModel/Multiplicity.h"

#include <experimental/type_traits>

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
    // return static_cast<typename binningType::binned_t>(((valueToBin - (binningType::binned_max - binningType::binned_min) * 0.5) / binningType::bin_width));
  }
}

template <typename binningType>
float unPack(const typename binningType::binned_t& b)
{
  return static_cast<float>(binningType::bin_width * b);

  // return static_cast<float>((binningType::binned_max - binningType::binned_min) * 0.5 + binningType::bin_width * b);
}

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

DECLARE_SOA_COLUMN(Mult, mult, int);           // Multiplicity of the collision
DECLARE_SOA_COLUMN(PosZ, posZ, float);         // Vertex of the collision
DECLARE_SOA_COLUMN(MagField, magField, float); // Magnetic field corresponding to a collision (in T)

} // namespace singletrackselector

DECLARE_SOA_TABLE(SingleCollSels, "AOD", "SCSEL", // Table of the variables for single track selection.
                  o2::soa::Index<>,
                  singletrackselector::Mult,
                  singletrackselector::PosZ,
                  singletrackselector::MagField);

namespace singletrackselector
{

DECLARE_SOA_INDEX_COLUMN(SingleCollSel, singleCollSel);                                  // Index to the collision
DECLARE_SOA_COLUMN(P, p, float);                                                         // Momentum of the track
DECLARE_SOA_COLUMN(DcaXY, dcaXY, float);                                                 // impact parameter of the track
DECLARE_SOA_COLUMN(DcaZ, dcaZ, float);                                                   // impact parameter of the track
DECLARE_SOA_COLUMN(TPCInnerParam, tpcInnerParam, float);                                 // Momentum at inner wall of the TPC
DECLARE_SOA_COLUMN(TPCSignal, tpcSignal, float);                                         // dE/dx TPC
DECLARE_SOA_COLUMN(Beta, beta, float);                                                   // TOF beta
DECLARE_SOA_COLUMN(TPCNClsFound, tpcNClsFound, int16_t);                                 // Number of TPC clusters
DECLARE_SOA_COLUMN(TPCChi2NCl, tpcChi2NCl, float);                                       // TPC chi2
DECLARE_SOA_COLUMN(TPCCrossedRowsOverFindableCls, tpcCrossedRowsOverFindableCls, float); // Ratio of found over findable clusters
DECLARE_SOA_COLUMN(TPCNClsShared, tpcNClsShared, uint8_t);                               // Number of shared TPC clusters
DECLARE_SOA_COLUMN(ITSNCls, itsNCls, uint8_t);                                           // Number of ITS clusters
DECLARE_SOA_COLUMN(ITSChi2NCl, itsChi2NCl, float);                                       // ITS chi2
DECLARE_SOA_COLUMN(Sign, sign, int8_t);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(StoredTOFNSigmaPr, storedTofNSigmaPr, nsigma::binning::binned_t);
DECLARE_SOA_COLUMN(StoredTPCNSigmaPr, storedTpcNSigmaPr, nsigma::binning::binned_t);
DECLARE_SOA_COLUMN(StoredTOFNSigmaDe, storedTofNSigmaDe, nsigma::binning::binned_t);
DECLARE_SOA_COLUMN(StoredTPCNSigmaDe, storedTpcNSigmaDe, nsigma::binning::binned_t);

DECLARE_SOA_DYNAMIC_COLUMN(Energy, energy,
                           [](float p, float mass) -> float { return sqrt(p * p + mass * mass); });
DECLARE_SOA_DYNAMIC_COLUMN(Pt, pt,
                           [](float p, float eta) -> float { return p / std::cosh(eta); });
DECLARE_SOA_DYNAMIC_COLUMN(Px, px,
                           [](float p, float eta, float phi) -> float { return (p / std::cosh(eta)) * std::sin(phi); });
DECLARE_SOA_DYNAMIC_COLUMN(Py, py,
                           [](float p, float eta, float phi) -> float { return (p / std::cosh(eta)) * std::cos(phi); });
DECLARE_SOA_DYNAMIC_COLUMN(Pz, pz,
                           [](float p, float eta) -> float { return p * std::tanh(eta); });
DECLARE_SOA_DYNAMIC_COLUMN(PhiStar, phiStar,
                           [](float p, float eta, float sign, float phi, float magfield = 0.0, float radius = 1.6) -> float {
                            if(magfield==0.0) return -1000.0;
                            else return phi + std::asin( -0.3*magfield*sign*radius/(2.0*p/std::cosh(eta)) ); });

DECLARE_SOA_DYNAMIC_COLUMN(TOFNSigmaPr, tofNSigmaPr,
                           [](nsigma::binning::binned_t nsigma_binned) -> float { return singletrackselector::unPack<nsigma::binning>(nsigma_binned); });
DECLARE_SOA_DYNAMIC_COLUMN(TPCNSigmaPr, tpcNSigmaPr,
                           [](nsigma::binning::binned_t nsigma_binned) -> float { return singletrackselector::unPack<nsigma::binning>(nsigma_binned); });
DECLARE_SOA_DYNAMIC_COLUMN(TOFNSigmaDe, tofNSigmaDe,
                           [](nsigma::binning::binned_t nsigma_binned) -> float { return singletrackselector::unPack<nsigma::binning>(nsigma_binned); });
DECLARE_SOA_DYNAMIC_COLUMN(TPCNSigmaDe, tpcNSigmaDe,
                           [](nsigma::binning::binned_t nsigma_binned) -> float { return singletrackselector::unPack<nsigma::binning>(nsigma_binned); });

} // namespace singletrackselector

DECLARE_SOA_TABLE_FULL(SingleTrackSels, "SelTracks", "AOD", "STSEL", // Table of the variables for single track selection.
                       o2::soa::Index<>,
                       singletrackselector::SingleCollSelId,
                       singletrackselector::P,
                       singletrackselector::DcaXY,
                       singletrackselector::DcaZ,
                       singletrackselector::TPCInnerParam,
                       singletrackselector::TPCSignal,
                       singletrackselector::Beta,
                       singletrackselector::TPCNClsFound,
                       singletrackselector::TPCChi2NCl,
                       singletrackselector::TPCCrossedRowsOverFindableCls,
                       singletrackselector::TPCNClsShared,
                       singletrackselector::ITSNCls,
                       singletrackselector::ITSChi2NCl,
                       singletrackselector::Sign,
                       singletrackselector::Eta,
                       singletrackselector::Phi,
                       singletrackselector::StoredTOFNSigmaPr,
                       singletrackselector::StoredTPCNSigmaPr,
                       singletrackselector::StoredTOFNSigmaDe,
                       singletrackselector::StoredTPCNSigmaDe,
                       singletrackselector::Energy<singletrackselector::P>,
                       singletrackselector::Pt<singletrackselector::P, singletrackselector::Eta>,
                       singletrackselector::Px<singletrackselector::P, singletrackselector::Eta, singletrackselector::Phi>,
                       singletrackselector::Py<singletrackselector::P, singletrackselector::Eta, singletrackselector::Phi>,
                       singletrackselector::Pz<singletrackselector::P, singletrackselector::Eta>,
                       singletrackselector::PhiStar<singletrackselector::P, singletrackselector::Eta, singletrackselector::Sign, singletrackselector::Phi>,
                       singletrackselector::TOFNSigmaPr<singletrackselector::StoredTOFNSigmaPr>,
                       singletrackselector::TPCNSigmaPr<singletrackselector::StoredTPCNSigmaPr>,
                       singletrackselector::TOFNSigmaDe<singletrackselector::StoredTOFNSigmaDe>,
                       singletrackselector::TPCNSigmaDe<singletrackselector::StoredTPCNSigmaDe>);

} // namespace o2::aod
#endif // PWGCF_DATAMODEL_SINGLETRACKSELECTOR_H_

namespace o2::aod::singletrackselector
{

template <typename TrackType>
inline bool TPCselection(TrackType const& track, std::pair<int, std::vector<float>> const& PIDcuts)
{
  // add check for the size of the vector and order of the values??? must be 2 valies in order => [down, up]

  float Nsigma = -1000;

  if (PIDcuts.first == 2212)
    Nsigma = track.tpcNSigmaPr();
  if (PIDcuts.first == 1000010020)
    Nsigma = track.tpcNSigmaDe();

  if (Nsigma > PIDcuts.second[0] && Nsigma < PIDcuts.second[1])
    return true;
  else
    return false;
}

template <typename TrackType>
inline bool TOFselection(TrackType const& track, std::pair<int, std::vector<float>> const& PIDcuts)
{
  // add check for the size of the vector and order of the values??? must be 2 valies in order => [down, up]

  float Nsigma = -1000;

  if (PIDcuts.first == 2212)
    Nsigma = track.tofNSigmaPr();
  else if (PIDcuts.first == 1000010020)
    Nsigma = track.tofNSigmaDe();

  else if (PIDcuts.first == 211) {
    if constexpr (std::experimental::is_detected<o2::aod::pidutils::hasTOFPi, TrackType>::value)
      Nsigma = track.tofNSigmaPi();
  } else if (PIDcuts.first == 321) {
    if constexpr (std::experimental::is_detected<o2::aod::pidutils::hasTOFKa, TrackType>::value)
      Nsigma = track.tofNSigmaKa();
  }

  if (Nsigma > PIDcuts.second[0] && Nsigma < PIDcuts.second[1])
    return true;
  else
    return false;
}
} // namespace o2::aod::singletrackselector