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
    // return static_cast<typename binningType::binned_t>(((valueToBin - (binningType::binned_max - binningType::binned_min) * 0.5) / binningType::bin_width));
  }
}

template <typename binningType>
float unPack(const typename binningType::binned_t& b)
{
  return static_cast<float>(binningType::bin_width * b);

  // return static_cast<float>((binningType::binned_max - binningType::binned_min) * 0.5 + binningType::bin_width * b);
}

template <typename binningType>
typename binningType::binned_t packInTableInt(const float& valueToBin)
{
  if (valueToBin <= binningType::binned_min) {
    return binningType::underflowBin;
  } else if (valueToBin >= binningType::binned_max) {
    return binningType::overflowBin;
  } else {
    // return static_cast<typename binningType::binned_t>(valueToBin / binningType::bin_width);
    return static_cast<typename binningType::binned_t>(((valueToBin - (binningType::binned_max - binningType::binned_min) * 0.5) / binningType::bin_width));
  }
}

template <typename binningType>
float unPackInt(const typename binningType::binned_t& b)
{
  // return static_cast<float>(binningType::bin_width * b);
  return static_cast<float>((binningType::binned_max - binningType::binned_min) * 0.5 + binningType::bin_width * b);
}
/*
// Function to pack a float into a binned value in table
template <typename binningType>
typename binningType::binned_t packInTable(const float& valueToBin)
{
  if (valueToBin <= binningType::binned_min) {
    return binningType::underflowBin;
  } else if (valueToBin >= binningType::binned_max) {
    return binningType::overflowBin;
  } else if (valueToBin >= 0) {
    return static_cast<typename binningType::binned_t>((valueToBin / binningType::bin_width) + 0.5f);
  } else {
    return static_cast<typename binningType::binned_t>((valueToBin / binningType::bin_width) - 0.5f);
  }*/

namespace storedcrossedrows
{
struct binning {
 public:
  typedef int8_t binned_t;
  // static constexpr int nbins = 160;
  //(1 << 8 * sizeof(binned_t)) - 2;
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
DECLARE_SOA_COLUMN(Px, px, float);                       // Momentum of the track
DECLARE_SOA_COLUMN(Py, py, float);                       // Momentum of the track
DECLARE_SOA_COLUMN(Pz, pz, float);                       // Momentum of the track
DECLARE_SOA_COLUMN(P, p, float);                         // Momentum of the track
DECLARE_SOA_COLUMN(Pt, pt, float);                       // Momentum of the track
DECLARE_SOA_COLUMN(TPCInnerParam, tpcInnerParam, float); // vertex position along z
DECLARE_SOA_COLUMN(TPCSignal, tpcSignal, float);         // vertex position along z
DECLARE_SOA_COLUMN(Beta, beta, float);
DECLARE_SOA_COLUMN(DcaXY, dcaXY, float);               // impact parameter of the track
DECLARE_SOA_COLUMN(DcaZ, dcaZ, float);                 // impact parameter of the track
DECLARE_SOA_COLUMN(TPCNClsFound, tpcNClsFound, float); // Number of TPC clusters
DECLARE_SOA_COLUMN(TPCCrossedRowsOverFindableCls, tpcCrossedRowsOverFindableCls, float);
DECLARE_SOA_COLUMN(TPCChi2NCl, tpcChi2NCl, float); // TPC chi2
DECLARE_SOA_COLUMN(ITSNCls, itsNCls, float);       // Number of ITS clusters
DECLARE_SOA_COLUMN(ITSChi2NCl, itsChi2NCl, float); // ITS chi2
DECLARE_SOA_COLUMN(Sign, sign, short);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(StoredCrossedRows, storedCrossedRows, storedcrossedrows::binning::binned_t);
DECLARE_SOA_COLUMN(StoredTOFNSigmaPr, storedTofNSigmaPr, nsigma::binning::binned_t);
DECLARE_SOA_COLUMN(StoredTPCNSigmaPr, storedTpcNSigmaPr, nsigma::binning::binned_t);
DECLARE_SOA_COLUMN(StoredTOFNSigmaDe, storedTofNSigmaDe, nsigma::binning::binned_t);
DECLARE_SOA_COLUMN(StoredTPCNSigmaDe, storedTpcNSigmaDe, nsigma::binning::binned_t);

DECLARE_SOA_DYNAMIC_COLUMN(CrossedRows, tpcNClsCrossedRows,
                           [](storedcrossedrows::binning::binned_t binned) -> float { return singletrackselector::unPackInt<storedcrossedrows::binning>(binned); });
DECLARE_SOA_DYNAMIC_COLUMN(TOFNSigmaPr, tofNSigmaPr,
                           [](nsigma::binning::binned_t nsigma_binned) -> float { return singletrackselector::unPack<nsigma::binning>(nsigma_binned); });
DECLARE_SOA_DYNAMIC_COLUMN(TPCNSigmaPr, tpcNSigmaPr,
                           [](nsigma::binning::binned_t nsigma_binned) -> float { return singletrackselector::unPack<nsigma::binning>(nsigma_binned); });
DECLARE_SOA_DYNAMIC_COLUMN(TOFNSigmaDe, tofNSigmaDe,
                           [](nsigma::binning::binned_t nsigma_binned) -> float { return singletrackselector::unPack<nsigma::binning>(nsigma_binned); });
DECLARE_SOA_DYNAMIC_COLUMN(TPCNSigmaDe, tpcNSigmaDe,
                           [](nsigma::binning::binned_t nsigma_binned) -> float { return singletrackselector::unPack<nsigma::binning>(nsigma_binned); });

DECLARE_SOA_DYNAMIC_COLUMN(Energy, energy,
                           [](float p, float mass) -> float { return sqrt(p * p + mass * mass); });

DECLARE_SOA_DYNAMIC_COLUMN(TrackCuts, trackCuts,
                           [](float p, float eta, float dcaXY, float dcaZ,
                              float tpcNClsFound, float tpcChi2NCl, float itsNCls, float tpcCrossedRowsOverFindableCls, storedcrossedrows::binning::binned_t storedCrossedRows,
                              std::map<std::string, float>* track_cuts) -> bool {
                            if(p < (*track_cuts)["min_P"] || p > (*track_cuts)["max_P"]) return false;
                            if(abs(eta) > (*track_cuts)["eta"]) return false;
                            if(tpcNClsFound < (*track_cuts)["tpcNClsFound"] || tpcChi2NCl > (*track_cuts)["tpcChi2NCl"]) return false;
                            if(abs(dcaXY) > (*track_cuts)["dcaXY"] || abs(dcaZ) > (*track_cuts)["dcaZ"]) return false;
                            if(itsNCls < (*track_cuts)["itsNCls"]) return false;
                            if(singletrackselector::unPackInt<storedcrossedrows::binning>(storedCrossedRows)<(*track_cuts)["crossedrows"]) return false;
                            if(tpcCrossedRowsOverFindableCls < (*track_cuts)["crossedRows/findableCls"]) return false;
                           // if(singletrackselector::unPack<storedcrossedrows::binning>(storedCrossedRows)> (*track_cuts)["crossedrows"]) return false;

                            return true; });

DECLARE_SOA_DYNAMIC_COLUMN(PIDCuts, pidCuts,
                           [](float pt, float sign, nsigma::binning::binned_t storedTpcNSigmaPr, nsigma::binning::binned_t storedTofNSigmaPr,
                              nsigma::binning::binned_t storedTpcNSigmaDe, nsigma::binning::binned_t storedTofNSigmaDe,
                              std::map<std::string, double>* PID_cuts) -> bool {
                            if(sign != (*PID_cuts)["sign"]) return false;
                            if(pt < (*PID_cuts)["PIDtrshld"]){
                              if((*PID_cuts)["particlePDG"] == 2212 && abs(singletrackselector::unPack<nsigma::binning>(storedTpcNSigmaPr)) > (*PID_cuts)["tpcNSigma"]) return false;
                              if((*PID_cuts)["particlePDG"] == 1000010020 && abs(singletrackselector::unPack<nsigma::binning>(storedTpcNSigmaDe)) > (*PID_cuts)["tpcNSigma"]) return false;
                            }
                            else{
                              if((*PID_cuts)["particlePDG"] == 2212 && sqrt((singletrackselector::unPack<nsigma::binning>(storedTpcNSigmaPr) * singletrackselector::unPack<nsigma::binning>(storedTpcNSigmaPr)) + 
                                                                            (singletrackselector::unPack<nsigma::binning>(storedTofNSigmaPr)*singletrackselector::unPack<nsigma::binning>(storedTofNSigmaPr)))> (*PID_cuts)["tpctofNSigma"]) return false;
                              //if((*PID_cuts)["particlePDG"] == 2212 && (abs(singletrackselector::unPack<nsigma::binning>(storedTpcNSigmaPr)) > (*PID_cuts)["tpctofNSigma"]
                              //                                || abs(singletrackselector::unPack<nsigma::binning>(storedTofNSigmaPr)) > (*PID_cuts)["tpctofNSigma"])) return false;
                              //if((*PID_cuts)["particlePDG"] == 1000010020 && (abs(singletrackselector::unPack<nsigma::binning>(storedTpcNSigmaDe)) > (*PID_cuts)["tpctofNSigma"]
                              //                                || abs(singletrackselector::unPack<nsigma::binning>(storedTofNSigmaDe)) > (*PID_cuts)["tpctofNSigma"])) return false;
                            if((*PID_cuts)["particlePDG"] == 1000010020 && sqrt((singletrackselector::unPack<nsigma::binning>(storedTpcNSigmaDe) * singletrackselector::unPack<nsigma::binning>(storedTpcNSigmaDe)) + 
                                                                            (singletrackselector::unPack<nsigma::binning>(storedTofNSigmaDe)*singletrackselector::unPack<nsigma::binning>(storedTofNSigmaDe))) > (*PID_cuts)["tpctofNSigma"]) return false;
                            }

                            return true; });

// DECLARE_SOA_DYNAMIC_COLUMN(TOFNSigmaPr, tofNSigmaPr,
//                            [](nsigma::binning::binned_t nsigma_binned) -> float { return nsigma::binning::bin_width * static_cast<float>(nsigma_binned); });
// DECLARE_SOA_DYNAMIC_COLUMN(TPCNSigmaPr, tpcNSigmaPr,
//                            [](nsigma::binning::binned_t nsigma_binned) -> float { return nsigma::binning::bin_width * static_cast<float>(nsigma_binned); });

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
                  singletrackselector::P,
                  singletrackselector::Pt,
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
                  singletrackselector::CrossedRows<singletrackselector::StoredCrossedRows>,
                  singletrackselector::TOFNSigmaPr<singletrackselector::StoredTOFNSigmaPr>,
                  singletrackselector::TPCNSigmaPr<singletrackselector::StoredTPCNSigmaPr>,
                  singletrackselector::TOFNSigmaDe<singletrackselector::StoredTOFNSigmaDe>,
                  singletrackselector::TPCNSigmaDe<singletrackselector::StoredTPCNSigmaDe>,
                  singletrackselector::Energy<singletrackselector::P>,
                  singletrackselector::TrackCuts<singletrackselector::P, singletrackselector::Eta, singletrackselector::DcaXY,
                                                 singletrackselector::DcaZ, singletrackselector::TPCNClsFound, singletrackselector::TPCChi2NCl,
                                                 singletrackselector::ITSNCls,
                                                 singletrackselector::TPCCrossedRowsOverFindableCls,
                                                 singletrackselector::StoredCrossedRows>,
                  singletrackselector::PIDCuts<singletrackselector::Pt, singletrackselector::Sign,
                                               singletrackselector::StoredTPCNSigmaPr,
                                               singletrackselector::StoredTOFNSigmaPr,
                                               singletrackselector::StoredTPCNSigmaDe,
                                               singletrackselector::StoredTOFNSigmaDe>);

DECLARE_SOA_TABLE(SingleCollSel, "AOD", "SCSEL", // Table of the variables for single track selection.
                                                 // o2::soa::Index<>,
                  singletrackselector::GlobalIndex,
                  singletrackselector::Mult,
                  singletrackselector::PosZ);
} // namespace o2::aod
#endif // PWGCF_DATAMODEL_SINGLETRACKSELECTOR_H_
