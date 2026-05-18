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
/// \file LongRangeDerived.h
///
/// \brief task derived table definition for long range correlation
/// \author Abhi Modak (abhi.modak@cern.ch)
/// \since October 28, 2025

#ifndef PWGCF_TWOPARTICLECORRELATIONS_DATAMODEL_LONGRANGEDERIVED_H_
#define PWGCF_TWOPARTICLECORRELATIONS_DATAMODEL_LONGRANGEDERIVED_H_

#include "PWGMM/Mult/DataModel/bestCollisionTable.h"

#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CommonConstants/MathConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>

#include <cstdint>

namespace o2::aod
{
namespace lrcorrmccolltable
{
DECLARE_SOA_COLUMN(Multiplicity, multiplicity, float);
} // namespace lrcorrmccolltable
DECLARE_SOA_TABLE(LRMcCollisions, "AOD", "LRMCCOLLISION",
                  o2::soa::Index<>,
                  mccollision::PosZ,
                  lrcorrmccolltable::Multiplicity,
                  mult::MultMCFT0A,
                  mult::MultMCFT0C);
using LRMcCollision = LRMcCollisions::iterator;

namespace lrcorrmctrktable
{
DECLARE_SOA_INDEX_COLUMN(LRMcCollision, lrMcCollision);
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
} // namespace lrcorrmctrktable

DECLARE_SOA_TABLE(LRMidMcTracks, "AOD", "LRMIDMCTRACK",
                  o2::soa::Index<>,
                  lrcorrmctrktable::LRMcCollisionId,
                  lrcorrmctrktable::Pt,
                  lrcorrmctrktable::Eta,
                  lrcorrmctrktable::Phi,
                  mcparticle::PdgCode,
                  mcparticle::Flags,
                  mcparticle::IsPhysicalPrimary<mcparticle::Flags>);
using LRMidMcTrack = LRMidMcTracks::iterator;

DECLARE_SOA_TABLE(LRFt0aMcTracks, "AOD", "LRFT0AMCTRACK",
                  o2::soa::Index<>,
                  lrcorrmctrktable::LRMcCollisionId,
                  lrcorrmctrktable::Pt,
                  lrcorrmctrktable::Eta,
                  lrcorrmctrktable::Phi);
using LRFt0aMcTrack = LRFt0aMcTracks::iterator;

DECLARE_SOA_TABLE(LRFt0cMcTracks, "AOD", "LRFT0CMCTRACK",
                  o2::soa::Index<>,
                  lrcorrmctrktable::LRMcCollisionId,
                  lrcorrmctrktable::Pt,
                  lrcorrmctrktable::Eta,
                  lrcorrmctrktable::Phi);
using LRFt0cMcTrack = LRFt0cMcTracks::iterator;

DECLARE_SOA_TABLE(LRMftMcTracks, "AOD", "LRMFTMCTRACK",
                  o2::soa::Index<>,
                  lrcorrmctrktable::LRMcCollisionId,
                  lrcorrmctrktable::Pt,
                  lrcorrmctrktable::Eta,
                  lrcorrmctrktable::Phi);
using LRMftMcTrack = LRMftMcTracks::iterator;

namespace lrcorrcolltable
{
DECLARE_SOA_INDEX_COLUMN(LRMcCollision, lrMcCollision);
DECLARE_SOA_COLUMN(Multiplicity, multiplicity, float);
DECLARE_SOA_COLUMN(Centrality, centrality, float);
DECLARE_SOA_COLUMN(TotalFT0AmplitudeA, totalFT0AmplitudeA, float); //! sum of amplitudes on A side of FT0
DECLARE_SOA_COLUMN(TotalFT0AmplitudeC, totalFT0AmplitudeC, float); //! sum of amplitudes on C side of FT0
DECLARE_SOA_COLUMN(TotalFV0AmplitudeA, totalFV0AmplitudeA, float); //! sum of amplitudes on A side of FDD
DECLARE_SOA_COLUMN(GapSide, gapSide, uint8_t);                     // 0 for side A, 1 for side C, 2 for both sides
} // namespace lrcorrcolltable

DECLARE_SOA_TABLE(LRCollisions, "AOD", "LRCOLLISION",
                  o2::soa::Index<>,
                  bc::RunNumber,
                  collision::PosZ,
                  lrcorrcolltable::Multiplicity,
                  lrcorrcolltable::Centrality,
                  timestamp::Timestamp);
DECLARE_SOA_TABLE(LRCollLabels, "AOD", "LRCOLLLABEL",
                  lrcorrcolltable::LRMcCollisionId);
using LRCollision = LRCollisions::iterator;
using LRCollLabel = LRCollLabels::iterator;
using LRCollisionsWithLabel = soa::Join<LRCollisions, LRCollLabels>;
using LRCollisionWithLabel = LRCollisionsWithLabel::iterator;

DECLARE_SOA_TABLE(UpcLRCollisions, "AOD", "UPCLRCOLLISION",
                  o2::soa::Index<>,
                  bc::GlobalBC,
                  bc::RunNumber,
                  collision::PosZ,
                  lrcorrcolltable::Multiplicity,
                  lrcorrcolltable::TotalFT0AmplitudeA,
                  lrcorrcolltable::TotalFT0AmplitudeC,
                  lrcorrcolltable::TotalFV0AmplitudeA);
using UpcLRCollision = UpcLRCollisions::iterator;

DECLARE_SOA_TABLE(UpcSgLRCollisions, "AOD", "UPCSGLRCOLLISION",
                  lrcorrcolltable::GapSide);
using UpcSgLRCollision = UpcSgLRCollisions::iterator;

namespace lrcorrzdctable
{
DECLARE_SOA_INDEX_COLUMN(UpcLRCollision, upcLRCollision);
DECLARE_SOA_COLUMN(EnergyCommonZNA, energyCommonZNA, float);
DECLARE_SOA_COLUMN(EnergyCommonZNC, energyCommonZNC, float);
} // namespace lrcorrzdctable

DECLARE_SOA_TABLE(LRZdcs, "AOD", "LRZDC",
                  o2::soa::Index<>,
                  lrcorrzdctable::UpcLRCollisionId,
                  lrcorrzdctable::EnergyCommonZNA,
                  lrcorrzdctable::EnergyCommonZNC);
using LRZdc = LRZdcs::iterator;

namespace lrcorrtrktable
{

template <typename binningType>
inline typename binningType::binned_t packInTable(const float& valueToBin)
{
  if (valueToBin <= binningType::binned_min) {
    return binningType::underflowBin;
  } else if (valueToBin >= binningType::binned_max) {
    return binningType::overflowBin;
  } else {
    return static_cast<typename binningType::binned_t>((valueToBin - binningType::binned_min) / binningType::bin_width);
  }
}

template <typename binningType>
inline float unPack(const typename binningType::binned_t& b)
{
  return binningType::bin_width * b + binningType::binned_min;
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

namespace binning
{

template <std::pair<float, float> lim, typename binVariable = int8_t>
struct binningParent {
 public:
  typedef binVariable binned_t;

  // Reserve two bins: one for overflow and one for underflow
  static constexpr int nbins = (1 << (8 * sizeof(binned_t))) - 2;
  static constexpr binned_t overflowBin = nbins;
  static constexpr binned_t underflowBin = -1;
  static constexpr float binned_min = lim.first;
  static constexpr float binned_max = lim.second;
  static constexpr float binned_center = 0.5 * (binned_min + binned_max);
  static constexpr float bin_width = (binned_max - binned_min) / static_cast<float>(nbins);
  static constexpr float inv_bin_width = 1. / bin_width;
  static_assert(binned_min < binned_max, "Invalid binning range");
  static void print()
  {
    LOG(info) << "Binning: " << binned_min << " - " << binned_max << " with " << nbins << " bins, width = "
              << bin_width << ". Overflow bin " << static_cast<int>(overflowBin) << " Underflow bin " << static_cast<int>(underflowBin);
  }
};

using trkdca_v0 = binningParent<std::pair<float, float>(-2.0f, 2.0f), int8_t>;
using trkphi_v0 = binningParent<std::pair<float, float>(0.0f, o2::constants::math::TwoPI), uint16_t>;
using trkamp_v0 = binningParent<std::pair<float, float>(0.0f, 5000.0f), uint16_t>;
using trkpt_v0 = binningParent<std::pair<float, float>(0.0f, 10.0f), uint8_t>;
using trketa_v0 = binningParent<std::pair<float, float>(-5.0f, 5.0f), int16_t>;
using trkchi2_v0 = binningParent<std::pair<float, float>(0.0f, 10.0f), int8_t>;

using trkdca = trkdca_v0;
using trkchi2 = trkchi2_v0;
using trkpt = trkpt_v0;
using trketa = trketa_v0;
using trkphi = trkphi_v0;
using trkamp = trkamp_v0;

} // namespace binning

DECLARE_SOA_INDEX_COLUMN(LRCollision, lrCollision);
DECLARE_SOA_INDEX_COLUMN(UpcLRCollision, upcLRCollision);
DECLARE_SOA_COLUMN(ChannelID, channelID, uint8_t);
DECLARE_SOA_COLUMN(AmplitudeStore, amplitudeStore, binning::trkamp::binned_t);
DECLARE_SOA_COLUMN(InvMass, invMass, float);
DECLARE_SOA_COLUMN(IdPos, idPos, int64_t);
DECLARE_SOA_COLUMN(IdNeg, idNeg, int64_t);
DECLARE_SOA_COLUMN(TrackType, trackType, uint8_t);
DECLARE_SOA_COLUMN(V0Type, v0Type, uint8_t);
DECLARE_SOA_COLUMN(AmbDegree, ambDegree, uint8_t);
DECLARE_SOA_COLUMN(TPCNClsFound, tpcNClsFound, uint8_t);
DECLARE_SOA_COLUMN(TPCNClsCrossedRows, tpcNClsCrossedRows, uint8_t);
DECLARE_SOA_COLUMN(TPCChi2NClStore, tpcChi2NClStore, binning::trkchi2::binned_t); //! Stored binned chi2
DECLARE_SOA_COLUMN(DCAzStore, dcazStore, binning::trkdca::binned_t);              //! Stored binned dcaz
DECLARE_SOA_COLUMN(BestDCAxyStore, bestdcaxyStore, binning::trkdca::binned_t);    //! Stored binned best dcaxy
DECLARE_SOA_COLUMN(BestDCAzStore, bestdcazStore, binning::trkdca::binned_t);      //! Stored binned best dcaz
DECLARE_SOA_COLUMN(PtStore, ptStore, binning::trkpt::binned_t);                   //! Stored binned pt
DECLARE_SOA_COLUMN(EtaStore, etaStore, binning::trketa::binned_t);                //! Stored binned eta
DECLARE_SOA_COLUMN(PhiStore, phiStore, binning::trkphi::binned_t);                //! Stored binned phi
DECLARE_SOA_DYNAMIC_COLUMN(TPCChi2NCl, tpcChi2NCl,
                           [](binning::trkchi2::binned_t chi2_binned) -> float { return unPack<binning::trkchi2>(chi2_binned); });
DECLARE_SOA_DYNAMIC_COLUMN(DCAz, dcaZ,
                           [](binning::trkdca::binned_t dcaz_binned) -> float { return unPackSymmetric<binning::trkdca>(dcaz_binned); });
DECLARE_SOA_DYNAMIC_COLUMN(Pt, pt,
                           [](binning::trkpt::binned_t pt_binned) -> float { return unPack<binning::trkpt>(pt_binned); });
DECLARE_SOA_DYNAMIC_COLUMN(Eta, eta,
                           [](binning::trketa::binned_t eta_binned) -> float { return unPackSymmetric<binning::trketa>(eta_binned); });
DECLARE_SOA_DYNAMIC_COLUMN(Phi, phi,
                           [](binning::trkphi::binned_t phi_binned) -> float { return unPackSymmetric<binning::trkphi>(phi_binned); });
DECLARE_SOA_DYNAMIC_COLUMN(BestDCAXY, bestDCAXY,
                           [](binning::trkdca::binned_t bestdcaxy_binned) -> float { return unPackSymmetric<binning::trkdca>(bestdcaxy_binned); });
DECLARE_SOA_DYNAMIC_COLUMN(BestDCAZ, bestDCAZ,
                           [](binning::trkdca::binned_t bestdcaz_binned) -> float { return unPackSymmetric<binning::trkdca>(bestdcaz_binned); });
DECLARE_SOA_DYNAMIC_COLUMN(Amplitude, amplitude,
                           [](binning::trkamp::binned_t amp_binned) -> float { return unPack<binning::trkamp>(amp_binned); });
enum TrackPid {
  kSpCharge,
  kSpPion,
  kSpKaon,
  kSpProton,
  kNoPid
};
enum V0TrackPid {
  kSpK0short,
  kSpLambda,
  kSpALambda
};
} // namespace lrcorrtrktable

DECLARE_SOA_TABLE(LRMidTracks, "AOD", "LRMIDTRACK",
                  o2::soa::Index<>,
                  lrcorrtrktable::LRCollisionId,
                  lrcorrtrktable::TPCNClsFound,
                  lrcorrtrktable::TPCNClsCrossedRows,
                  lrcorrtrktable::TPCChi2NClStore,
                  lrcorrtrktable::PtStore,
                  lrcorrtrktable::EtaStore,
                  lrcorrtrktable::PhiStore,
                  lrcorrtrktable::DCAzStore,
                  lrcorrtrktable::TrackType,
                  lrcorrtrktable::TPCChi2NCl<lrcorrtrktable::TPCChi2NClStore>,
                  lrcorrtrktable::Pt<lrcorrtrktable::PtStore>,
                  lrcorrtrktable::Eta<lrcorrtrktable::EtaStore>,
                  lrcorrtrktable::Phi<lrcorrtrktable::PhiStore>,
                  lrcorrtrktable::DCAz<lrcorrtrktable::DCAzStore>);
using LRMidTrack = LRMidTracks::iterator;

DECLARE_SOA_TABLE(LRFt0aTracks, "AOD", "LRFT0ATRACK",
                  o2::soa::Index<>,
                  lrcorrtrktable::LRCollisionId,
                  lrcorrtrktable::ChannelID,
                  lrcorrtrktable::AmplitudeStore,
                  lrcorrtrktable::EtaStore,
                  lrcorrtrktable::PhiStore,
                  lrcorrtrktable::Amplitude<lrcorrtrktable::AmplitudeStore>,
                  lrcorrtrktable::Eta<lrcorrtrktable::EtaStore>,
                  lrcorrtrktable::Phi<lrcorrtrktable::PhiStore>);
using LRFt0aTrack = LRFt0aTracks::iterator;

DECLARE_SOA_TABLE(LRFt0cTracks, "AOD", "LRFT0CTRACK",
                  o2::soa::Index<>,
                  lrcorrtrktable::LRCollisionId,
                  lrcorrtrktable::ChannelID,
                  lrcorrtrktable::AmplitudeStore,
                  lrcorrtrktable::EtaStore,
                  lrcorrtrktable::PhiStore,
                  lrcorrtrktable::Amplitude<lrcorrtrktable::AmplitudeStore>,
                  lrcorrtrktable::Eta<lrcorrtrktable::EtaStore>,
                  lrcorrtrktable::Phi<lrcorrtrktable::PhiStore>);
using LRFt0cTrack = LRFt0cTracks::iterator;

DECLARE_SOA_TABLE(LRV0Tracks, "AOD", "LRV0TRACK",
                  o2::soa::Index<>,
                  lrcorrtrktable::LRCollisionId,
                  lrcorrtrktable::IdPos,
                  lrcorrtrktable::IdNeg,
                  lrcorrtrktable::PtStore,
                  lrcorrtrktable::EtaStore,
                  lrcorrtrktable::PhiStore,
                  lrcorrtrktable::InvMass,
                  lrcorrtrktable::V0Type,
                  lrcorrtrktable::Pt<lrcorrtrktable::PtStore>,
                  lrcorrtrktable::Eta<lrcorrtrktable::EtaStore>,
                  lrcorrtrktable::Phi<lrcorrtrktable::PhiStore>);
using LRV0Track = LRV0Tracks::iterator;

DECLARE_SOA_TABLE(LRMftTracks, "AOD", "LRMFTTRACK",
                  o2::soa::Index<>,
                  lrcorrtrktable::LRCollisionId,
                  lrcorrtrktable::AmbDegree,
                  lrcorrtrktable::PtStore,
                  lrcorrtrktable::EtaStore,
                  lrcorrtrktable::PhiStore,
                  fwdtrack::NClusters,
                  lrcorrtrktable::BestDCAxyStore,
                  lrcorrtrktable::BestDCAzStore,
                  lrcorrtrktable::BestDCAXY<lrcorrtrktable::BestDCAxyStore>,
                  lrcorrtrktable::BestDCAZ<lrcorrtrktable::BestDCAzStore>,
                  lrcorrtrktable::Pt<lrcorrtrktable::PtStore>,
                  lrcorrtrktable::Eta<lrcorrtrktable::EtaStore>,
                  lrcorrtrktable::Phi<lrcorrtrktable::PhiStore>);
using LRMftTrack = LRMftTracks::iterator;

DECLARE_SOA_TABLE(UpcLRMidTracks, "AOD", "UPCLRMIDTRACK",
                  o2::soa::Index<>,
                  lrcorrtrktable::UpcLRCollisionId,
                  lrcorrtrktable::TPCNClsFound,
                  lrcorrtrktable::TPCNClsCrossedRows,
                  lrcorrtrktable::TPCChi2NClStore,
                  lrcorrtrktable::PtStore,
                  lrcorrtrktable::EtaStore,
                  lrcorrtrktable::PhiStore,
                  lrcorrtrktable::DCAzStore,
                  lrcorrtrktable::TrackType,
                  lrcorrtrktable::TPCChi2NCl<lrcorrtrktable::TPCChi2NClStore>,
                  lrcorrtrktable::Pt<lrcorrtrktable::PtStore>,
                  lrcorrtrktable::Eta<lrcorrtrktable::EtaStore>,
                  lrcorrtrktable::Phi<lrcorrtrktable::PhiStore>,
                  lrcorrtrktable::DCAz<lrcorrtrktable::DCAzStore>);
using UpcLRMidTrack = UpcLRMidTracks::iterator;

DECLARE_SOA_TABLE(UpcLRFt0aTracks, "AOD", "UPCLRFT0ATRACK",
                  o2::soa::Index<>,
                  lrcorrtrktable::UpcLRCollisionId,
                  lrcorrtrktable::ChannelID,
                  lrcorrtrktable::AmplitudeStore,
                  lrcorrtrktable::EtaStore,
                  lrcorrtrktable::PhiStore,
                  lrcorrtrktable::Amplitude<lrcorrtrktable::AmplitudeStore>,
                  lrcorrtrktable::Eta<lrcorrtrktable::EtaStore>,
                  lrcorrtrktable::Phi<lrcorrtrktable::PhiStore>);
using UpcLRFt0aTrack = UpcLRFt0aTracks::iterator;

DECLARE_SOA_TABLE(UpcLRFt0cTracks, "AOD", "UPCLRFT0CTRACK",
                  o2::soa::Index<>,
                  lrcorrtrktable::UpcLRCollisionId,
                  lrcorrtrktable::ChannelID,
                  lrcorrtrktable::AmplitudeStore,
                  lrcorrtrktable::EtaStore,
                  lrcorrtrktable::PhiStore,
                  lrcorrtrktable::Amplitude<lrcorrtrktable::AmplitudeStore>,
                  lrcorrtrktable::Eta<lrcorrtrktable::EtaStore>,
                  lrcorrtrktable::Phi<lrcorrtrktable::PhiStore>);
using UpcLRFt0cTrack = UpcLRFt0cTracks::iterator;

DECLARE_SOA_TABLE(UpcLRV0Tracks, "AOD", "UPCLRV0TRACK",
                  o2::soa::Index<>,
                  lrcorrtrktable::UpcLRCollisionId,
                  lrcorrtrktable::IdPos,
                  lrcorrtrktable::IdNeg,
                  lrcorrtrktable::PtStore,
                  lrcorrtrktable::EtaStore,
                  lrcorrtrktable::PhiStore,
                  lrcorrtrktable::InvMass,
                  lrcorrtrktable::V0Type,
                  lrcorrtrktable::Pt<lrcorrtrktable::PtStore>,
                  lrcorrtrktable::Eta<lrcorrtrktable::EtaStore>,
                  lrcorrtrktable::Phi<lrcorrtrktable::PhiStore>);
using UpcLRV0Track = UpcLRV0Tracks::iterator;

DECLARE_SOA_TABLE(UpcLRMftTracks, "AOD", "UPCLRMFTTRACK",
                  o2::soa::Index<>,
                  lrcorrtrktable::UpcLRCollisionId,
                  lrcorrtrktable::AmbDegree,
                  lrcorrtrktable::PtStore,
                  lrcorrtrktable::EtaStore,
                  lrcorrtrktable::PhiStore,
                  fwdtrack::NClusters,
                  lrcorrtrktable::BestDCAxyStore,
                  lrcorrtrktable::BestDCAzStore,
                  lrcorrtrktable::BestDCAXY<lrcorrtrktable::BestDCAxyStore>,
                  lrcorrtrktable::BestDCAZ<lrcorrtrktable::BestDCAzStore>,
                  lrcorrtrktable::Pt<lrcorrtrktable::PtStore>,
                  lrcorrtrktable::Eta<lrcorrtrktable::EtaStore>,
                  lrcorrtrktable::Phi<lrcorrtrktable::PhiStore>);
using UpcLRMftTrack = UpcLRMftTracks::iterator;
} // namespace o2::aod

#endif // PWGCF_TWOPARTICLECORRELATIONS_DATAMODEL_LONGRANGEDERIVED_H_
