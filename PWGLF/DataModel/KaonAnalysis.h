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
/// \file   spectraTOF.h
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
/// \since  2022-12-03
/// \brief  Header for the spectraTOF task for the analysis of the spectra with the TOF and TPC detectors.
///

#ifndef PWGLF_DATAMODEL_KAONANALYSIS_H_
#define PWGLF_DATAMODEL_KAONANALYSIS_H_

#include <memory>

// O2 includes
#include "ReconstructionDataFormats/Track.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/StaticFor.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/FT0Corrected.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "PWGLF/DataModel/LFParticleIdentification.h"

#include "TPDGCode.h"

static constexpr o2::track::PID::ID Np = 1;
static constexpr int NCharges = 2;
static constexpr o2::track::PID::ID NpCharge = Np * NCharges;
static constexpr const char* pT[Np] = {"K"};
static constexpr const char* pN[Np] = {"ka"};
static constexpr const char* cN[NCharges] = {"pos", "neg"};
static constexpr const char* pTCharge[NpCharge] = {"K^{+}",
                                                   "K^{-}"};
static constexpr int PDGs[NpCharge] = {kKPlus,
                                       -kKPlus};

std::shared_ptr<TH2> hMultiplicityvsPercentile;
static constexpr std::string_view hnsigmatpctof[NpCharge] = {"nsigmatpctof/pos/ka",
                                                             "nsigmatpctof/neg/ka"};
static constexpr std::string_view hnsigmatof[NpCharge] = {"nsigmatof/pos/ka",
                                                          "nsigmatof/neg/ka"};
static constexpr std::string_view hnsigmatpc[NpCharge] = {"nsigmatpc/pos/ka",
                                                          "nsigmatpc/neg/ka"};
static constexpr std::string_view hdeltatof[NpCharge] = {"deltatof/pos/ka",
                                                         "deltatof/neg/ka"};
static constexpr std::string_view hdeltatpc[NpCharge] = {"deltatpc/pos/ka",
                                                         "deltatpc/neg/ka"};
static constexpr std::string_view hdcaxy[NpCharge] = {"dcaxy/pos/ka",
                                                      "dcaxy/neg/ka"};
static constexpr std::string_view hdcaxytot[NpCharge] = {"dcaxytot/pos/ka",
                                                         "dcaxytot/neg/ka"};
static constexpr std::string_view hdcaz[NpCharge] = {"dcaz/pos/ka",
                                                     "dcaz/neg/ka"};
static constexpr std::string_view hdcaztot[NpCharge] = {"dcaztot/pos/ka",
                                                        "dcaztot/neg/ka"};
static constexpr std::string_view hdcaxyphi[NpCharge] = {"dcaxyphi/pos/ka",
                                                         "dcaxyphi/neg/ka"};
// MC
static constexpr std::string_view hpt_mism_its_prm[NpCharge] = {"MC/ka/pos/prm/pt/mismITS",
                                                                "MC/ka/neg/prm/pt/mismITS"};
static constexpr std::string_view hpt_mism_tpc_prm[NpCharge] = {"MC/ka/pos/prm/pt/mismTPC",
                                                                "MC/ka/neg/prm/pt/mismTPC"};
static constexpr std::string_view hpt_mism_trd_prm[NpCharge] = {"MC/ka/pos/prm/pt/mismTRD",
                                                                "MC/ka/neg/prm/pt/mismTRD"};
static constexpr std::string_view hpt_mism_tof_prm[NpCharge] = {"MC/ka/pos/prm/pt/mismTOF",
                                                                "MC/ka/neg/prm/pt/mismTOF"};
static constexpr std::string_view hpt_num_prm[NpCharge] = {"MC/ka/pos/prm/pt/num",
                                                           "MC/ka/neg/prm/pt/num"};
static constexpr std::string_view hpt_numtof_prm[NpCharge] = {"MC/ka/pos/prm/pt/numtof",
                                                              "MC/ka/neg/prm/pt/numtof"};
static constexpr std::string_view hpt_numtofgoodmatch_prm[NpCharge] = {"MC/ka/pos/prm/pt/numtofgoodmatch",
                                                                       "MC/ka/neg/prm/pt/numtofgoodmatch"};

//********************************************RD**********************************************************************************************
static constexpr std::string_view hpt_numtof_str[NpCharge] = {"MC/ka/pos/str/pt/numtof",
                                                              "MC/ka/neg/str/pt/numtof"};
static constexpr std::string_view hpt_numtof_mat[NpCharge] = {"MC/ka/pos/mat/pt/numtof",
                                                              "MC/ka/neg/mat/pt/numtof"};
static constexpr std::string_view hpt_den_prm[NpCharge] = {"MC/ka/pos/prm/pt/den",
                                                           "MC/ka/neg/prm/pt/den"};
static constexpr std::string_view hpt_den_prm_recoev[NpCharge] = {"MC/ka/pos/prm/pt/denrecoev",
                                                                  "MC/ka/neg/prm/pt/denrecoev"};
static constexpr std::string_view hpt_den_prm_evsel[NpCharge] = {"MC/ka/pos/prm/pt/denevsel",
                                                                 "MC/ka/neg/prm/pt/denevsel"};
static constexpr std::string_view hpt_den_prm_goodev[NpCharge] = {"MC/ka/pos/prm/pt/dengoodev",
                                                                  "MC/ka/neg/prm/pt/dengoodev"};
static constexpr std::string_view hpt_den_prm_mcgoodev[NpCharge] = {"MC/ka/pos/prm/pt/denmcgoodev",
                                                                    "MC/ka/neg/prm/pt/denmcgoodev"};
static constexpr std::string_view hpt_den_prm_mcbadev[NpCharge] = {"MC/ka/pos/prm/pt/denmcbadev",
                                                                   "MC/ka/neg/prm/pt/denmcbadev"};
static constexpr std::string_view hpt_num_str[NpCharge] = {"MC/ka/pos/str/pt/num",
                                                           "MC/ka/neg/str/pt/num"};
static constexpr std::string_view hpt_den_str[NpCharge] = {"MC/ka/pos/str/pt/den",
                                                           "MC/ka/neg/str/pt/den"};
static constexpr std::string_view hpt_num_mat[NpCharge] = {"MC/ka/pos/mat/pt/num",
                                                           "MC/ka/neg/mat/pt/num"};
static constexpr std::string_view hpt_den_mat[NpCharge] = {"MC/ka/pos/mat/pt/den",
                                                           "MC/ka/neg/mat/pt/den"};
static constexpr std::string_view hdcaxyprm[NpCharge] = {"dcaxyprm/pos/ka",
                                                         "dcaxyprm/neg/ka"};
static constexpr std::string_view hdcaxyprm2[NpCharge] = {"dcaxyprm2/pos/ka",
                                                          "dcaxyprm2/neg/ka"};
static constexpr std::string_view hdcaxyD0[NpCharge] = {"dcaxyD0/pos/ka",
                                                        "dcaxyD0/neg/ka"};
static constexpr std::string_view hdcaxyprmgoodevs[NpCharge] = {"dcaxyprmgoodevs/pos/ka",
                                                                "dcaxyprmgoodevs/neg/ka"};
static constexpr std::string_view hdcazprm[NpCharge] = {"dcazprm/pos/ka",
                                                        "dcazprm/neg/ka"};
static constexpr std::string_view hdcazprm2[NpCharge] = {"dcazprm2/pos/ka",
                                                         "dcazprm2/neg/ka"};
static constexpr std::string_view hdcazD0[NpCharge] = {"dcazD0/pos/ka",
                                                       "dcazD0/neg/ka"};
static constexpr std::string_view hdcazprmgoodevs[NpCharge] = {"dcazprmgoodevs/pos/ka",
                                                               "dcazprmgoodevs/neg/ka"};
static constexpr std::string_view hdcaxystr[NpCharge] = {"dcaxystr/pos/ka",
                                                         "dcaxystr/neg/ka"};
static constexpr std::string_view hdcaxycharm[NpCharge] = {"dcaxycharm/pos/ka",
                                                           "dcaxycharm/neg/ka"};
static constexpr std::string_view hdcaxybeauty[NpCharge] = {"dcaxybeauty/pos/ka",
                                                            "dcaxybeauty/neg/ka"};
static constexpr std::string_view hdcaxymat[NpCharge] = {"dcaxymat/pos/ka",
                                                         "dcaxymat/neg/ka"};
static constexpr std::string_view hdcazstr[NpCharge] = {"dcazstr/pos/ka",
                                                        "dcazstr/neg/ka"};
static constexpr std::string_view hdcazcharm[NpCharge] = {"dcazcharm/pos/ka",
                                                          "dcazcharm/neg/ka"};
static constexpr std::string_view hdcazbeauty[NpCharge] = {"dcazbeauty/pos/ka",
                                                           "dcazbeauty/neg/ka"};
static constexpr std::string_view hdcazmat[NpCharge] = {"dcazmat/pos/ka",
                                                        "dcazmat/neg/ka"};

// Derived data model for cut variation
namespace o2::aod
{
namespace spectra
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

struct binningNSigma {
 public:
  typedef int16_t binned_t;
  static constexpr int nbins = (1 << 8 * sizeof(binned_t)) - 2;
  static constexpr binned_t overflowBin = nbins >> 1;
  static constexpr binned_t underflowBin = -(nbins >> 1);
  static constexpr float binned_max = 10.0;
  static constexpr float binned_min = -10.0;
  static constexpr float bin_width = (binned_max - binned_min) / nbins;
};

// Collision info
DECLARE_SOA_INDEX_COLUMN(BC, bc); //! Most probably BC to where this collision has occurred
DECLARE_SOA_COLUMN(IsEventReject, isEventReject, int);
DECLARE_SOA_COLUMN(RunNumber, runNumber, int);
DECLARE_SOA_COLUMN(CentFT0M, centFT0M, float);
DECLARE_SOA_COLUMN(Sel8, sel8, bool);
DECLARE_SOA_COLUMN(MultNTracksPVeta1, multNTracksPVeta1, int);
DECLARE_SOA_DYNAMIC_COLUMN(CentFV0A, centFV0A, //! Dummy
                           [](bool /*v*/) -> float { return 0.f; });
DECLARE_SOA_DYNAMIC_COLUMN(CentFT0A, centFT0A, //! Dummy
                           [](bool /*v*/) -> float { return 0.f; });
DECLARE_SOA_DYNAMIC_COLUMN(CentFT0C, centFT0C, //! Dummy
                           [](bool /*v*/) -> float { return 0.f; });
DECLARE_SOA_DYNAMIC_COLUMN(MultZeqFV0A, multZeqFV0A, //! Dummy
                           [](bool /*v*/) -> float { return 0.f; });
DECLARE_SOA_DYNAMIC_COLUMN(MultZeqFT0A, multZeqFT0A, //! Dummy
                           [](bool /*v*/) -> float { return 0.f; });
DECLARE_SOA_DYNAMIC_COLUMN(MultZeqFT0C, multZeqFT0C, //! Dummy
                           [](bool /*v*/) -> float { return 0.f; });
DECLARE_SOA_DYNAMIC_COLUMN(MultZeqFDDA, multZeqFDDA, //! Dummy
                           [](bool /*v*/) -> float { return 0.f; });
DECLARE_SOA_DYNAMIC_COLUMN(MultZeqFDDC, multZeqFDDC, //! Dummy
                           [](bool /*v*/) -> float { return 0.f; });
DECLARE_SOA_DYNAMIC_COLUMN(MultZeqNTracksPV, multZeqNTracksPV, //! Dummy
                           [](bool /*v*/) -> float { return 0.f; });
DECLARE_SOA_DYNAMIC_COLUMN(MultTracklets, multTracklets, //! Dummy
                           [](bool /*v*/) -> float { return 0.f; });
DECLARE_SOA_DYNAMIC_COLUMN(MultTPC, multTPC, //! Dummy
                           [](bool /*v*/) -> float { return 0.f; });
DECLARE_SOA_DYNAMIC_COLUMN(SelectionBit, selection_bit, //! Dummy
                           [](aod::evsel::EventSelectionFlags /*v*/) -> bool { return true; });
DECLARE_SOA_DYNAMIC_COLUMN(IsInelGt0, isInelGt0, //! is INEL > 0
                           [](int multPveta1) -> bool { return multPveta1 > 0; });
DECLARE_SOA_DYNAMIC_COLUMN(IsInelGt1, isInelGt1, //! is INEL > 1
                           [](int multPveta1) -> bool { return multPveta1 > 1; });

// Track info
DECLARE_SOA_INDEX_COLUMN(Collision, collision);                                  //! Index to the collision
DECLARE_SOA_COLUMN(PtSigned, ptSigned, float);                                   //! Pt (signed) of the track
DECLARE_SOA_COLUMN(Eta, eta, float);                                             //! Eta of the track
DECLARE_SOA_COLUMN(Phi, phi, float);                                             //! Phi of the track
DECLARE_SOA_COLUMN(EvTimeT0AC, evTimeT0AC, float);                               //! Event time of the track computed with the T0AC
DECLARE_SOA_COLUMN(EvTimeT0ACErr, evTimeT0ACErr, float);                         //! Resolution of the event time of the track computed with the T0AC
DECLARE_SOA_COLUMN(IsPVContributor, isPVContributor, bool);                      //! IsPVContributor
DECLARE_SOA_COLUMN(DetectorMap, detectorMap, uint8_t);                           //! Detector map: see enum DetectorMapEnum
DECLARE_SOA_COLUMN(LastTRDCluster, lastTRDCluster, int8_t);                      //! Index of the last cluster in the TRD, -1 if no TRD information
DECLARE_SOA_COLUMN(HasTRD, hasTRD, bool);                                        //! Has or not the TRD match
DECLARE_SOA_COLUMN(TPCNSigmaStoreKa, tpcNSigmaStoreKa, binningNSigma::binned_t); //! Stored binned nsigma with the TPC detector for kaon
DECLARE_SOA_COLUMN(TOFNSigmaStoreKa, tofNSigmaStoreKa, binningNSigma::binned_t); //! Stored binned nsigma with the TOF detector for kaon
DECLARE_SOA_DYNAMIC_COLUMN(TPCNSigmaKa, tpcNSigmaKa,                             //! Unpacked NSigma TPC Ka
                           [](binningNSigma::binned_t binned) -> float { return unPack<binningNSigma>(binned); });
DECLARE_SOA_DYNAMIC_COLUMN(TOFNSigmaKa, tofNSigmaKa, //! Unpacked NSigma TOF Ka
                           [](binningNSigma::binned_t binned) -> float { return unPack<binningNSigma>(binned); });

DECLARE_SOA_COLUMN(DCAxyStore, dcaxyStore, binningDCA::binned_t); //! Stored binned dcaxy
DECLARE_SOA_COLUMN(DCAzStore, dcazStore, binningDCA::binned_t);   //! Stored binned dcaz
DECLARE_SOA_DYNAMIC_COLUMN(DCAxy, dcaXY,                          //! Unpacked dcaxy
                           [](binningDCA::binned_t binned) -> float { return unPack<binningDCA>(binned); });
DECLARE_SOA_DYNAMIC_COLUMN(DCAz, dcaZ, //! Unpacked dcaz
                           [](binningDCA::binned_t binned) -> float { return unPack<binningDCA>(binned); });
DECLARE_SOA_DYNAMIC_COLUMN(Pt, pt, //! Absolute value of signed pT
                           [](float signedPt) -> float { return std::abs(signedPt); });
DECLARE_SOA_DYNAMIC_COLUMN(HasITS, hasITS, //! Dummy
                           [](float /*v*/) -> bool { return true; });
DECLARE_SOA_DYNAMIC_COLUMN(HasTPC, hasTPC, //! Dummy
                           [](float /*v*/) -> bool { return true; });
DECLARE_SOA_DYNAMIC_COLUMN(HasTOF, hasTOF, //! Flag to check if track has a TOF measurement
                           [](float tofSignal) -> bool { return tofSignal > 0; });
DECLARE_SOA_DYNAMIC_COLUMN(TRDSignal, trdSignal, //! Dummy
                           [](float /*v*/) -> float { return 0.f; });
DECLARE_SOA_DYNAMIC_COLUMN(P, p, [](float signedpt, float eta) -> float { return std::abs(signedpt) * cosh(eta); });
DECLARE_SOA_DYNAMIC_COLUMN(TrackType, trackType, [](float /*v*/) -> uint8_t { return o2::aod::track::TrackTypeEnum::Track; });
DECLARE_SOA_COLUMN(IsGlobalTrack, isGlobalTrack, bool);                                       // if a track passed the isGlobalTrack requirement
DECLARE_SOA_COLUMN(IsGlobalTrackWoDCA, isGlobalTrackWoDCA, bool);                             // if a track passed the isGlobalTrackWoDCA requirement
DECLARE_SOA_DYNAMIC_COLUMN(Flags, flags, [](float /*v*/) -> uint32_t { return 0; });          // Dummy
DECLARE_SOA_DYNAMIC_COLUMN(TRDPattern, trdPattern, [](float /*v*/) -> uint8_t { return 0; }); // Dummy
DECLARE_SOA_DYNAMIC_COLUMN(Rapidity, rapidity,                                                //! Track rapidity, computed under the mass assumption given as input
                           [](float signedPt, float eta, float mass) -> float {
                             const auto pt = std::abs(signedPt);
                             const auto p = std::abs(signedPt) * cosh(eta);
                             const auto pz = std::sqrt(p * p - pt * pt);
                             const auto energy = sqrt(p * p + mass * mass);
                             return 0.5f * log((energy + pz) / (energy - pz));
                           });
DECLARE_SOA_DYNAMIC_COLUMN(IsInAcceptanceTrack, isInAcceptanceTrack, [](float /*v*/) -> bool { return false; }); // Dummy
DECLARE_SOA_DYNAMIC_COLUMN(IsQualityTrackITS, isQualityTrackITS, [](float /*v*/) -> bool { return false; });     // Dummy
DECLARE_SOA_DYNAMIC_COLUMN(IsQualityTrackTPC, isQualityTrackTPC, [](float /*v*/) -> bool { return false; });     // Dummy

} // namespace spectra

DECLARE_SOA_TABLE(SpColls, "AOD", "SPCOLLS",
                  o2::soa::Index<>,
                  collision::NumContrib,
                  collision::PosX,
                  collision::PosY,
                  collision::PosZ,
                  spectra::CentFT0M,
                  spectra::Sel8,
                  spectra::MultNTracksPVeta1,
                  spectra::RunNumber,
                  spectra::IsInelGt0<spectra::MultNTracksPVeta1>,
                  spectra::IsInelGt1<spectra::MultNTracksPVeta1>,
                  spectra::CentFV0A<spectra::Sel8>,
                  spectra::CentFT0A<spectra::Sel8>,
                  spectra::CentFT0C<spectra::Sel8>,
                  spectra::MultZeqFV0A<spectra::Sel8>,
                  spectra::MultZeqFT0A<spectra::Sel8>,
                  spectra::MultZeqFT0C<spectra::Sel8>,
                  spectra::MultZeqFDDA<spectra::Sel8>,
                  spectra::MultZeqFDDC<spectra::Sel8>,
                  spectra::MultZeqNTracksPV<spectra::Sel8>,
                  spectra::MultTracklets<spectra::Sel8>,
                  spectra::MultTPC<spectra::Sel8>,
                  spectra::SelectionBit<>);
using SpColl = SpColls::iterator;

DECLARE_SOA_TABLE(SpTracks, "AOD", "SPTRACKS",
                  o2::soa::Index<>,
                  spectra::CollisionId,
                  spectra::TPCNSigmaStoreKa,
                  spectra::TOFNSigmaStoreKa,
                  spectra::PtSigned, spectra::Eta, spectra::Phi,
                  track::Length,
                  track::TPCSignal,
                  track::TPCChi2NCl, track::ITSChi2NCl, track::TOFChi2,
                  track::TPCNClsShared,
                  track::TPCNClsFindable,
                  track::TPCNClsFindableMinusFound,
                  track::TPCNClsFindableMinusCrossedRows,
                  spectra::IsPVContributor,
                  track::ITSClusterSizes,
                  spectra::HasTRD,
                  //   pidtofevtime::EvTimeTOF,
                  //   pidtofevtime::EvTimeTOFErr,
                  //   pidtofevtime::EvTimeTOFMult,
                  // spectra::EvTimeT0AC,
                  // spectra::EvTimeT0ACErr,
                  // collision::CollisionTime,
                  // collision::CollisionTimeRes,
                  pidflags::TOFFlags,
                  spectra::DCAxyStore,
                  spectra::DCAzStore,
                  spectra::IsGlobalTrack,
                  spectra::IsGlobalTrackWoDCA,
                  spectra::DCAxy<spectra::DCAxyStore>,
                  spectra::DCAz<spectra::DCAzStore>,
                  spectra::Pt<spectra::PtSigned>,
                  track::Sign<spectra::PtSigned>,
                  spectra::P<spectra::PtSigned, spectra::Eta>,
                  spectra::Rapidity<spectra::PtSigned, spectra::Eta>,
                  spectra::HasITS<track::ITSClusterSizes>,
                  spectra::HasTPC<track::TPCChi2NCl>,
                  spectra::HasTOF<track::TOFChi2>,
                  spectra::TRDSignal<track::TOFChi2>,
                  spectra::Flags<track::TOFChi2>,
                  spectra::TrackType<track::TOFChi2>,
                  spectra::TRDPattern<track::TOFChi2>,
                  spectra::IsInAcceptanceTrack<track::TOFChi2>, // Dummy
                  spectra::IsQualityTrackITS<track::TOFChi2>,   // Dummy
                  spectra::IsQualityTrackTPC<track::TOFChi2>,   // Dummy
                  track::v001::ITSClusterMap<track::ITSClusterSizes>, track::v001::ITSNCls<track::ITSClusterSizes>, track::v001::ITSNClsInnerBarrel<track::ITSClusterSizes>,
                  track::TPCFractionSharedCls<track::TPCNClsShared, track::TPCNClsFindable, track::TPCNClsFindableMinusFound>,
                  track::TPCNClsFound<track::TPCNClsFindable, track::TPCNClsFindableMinusFound>,
                  track::TPCNClsCrossedRows<track::TPCNClsFindable, track::TPCNClsFindableMinusCrossedRows>,
                  track::TPCCrossedRowsOverFindableCls<track::TPCNClsFindable, track::TPCNClsFindableMinusCrossedRows>,
                  track::TPCFoundOverFindableCls<track::TPCNClsFindable, track::TPCNClsFindableMinusFound>,
                  pidflags::IsEvTimeDefined<pidflags::TOFFlags>,
                  pidflags::IsEvTimeTOF<pidflags::TOFFlags>,
                  pidflags::IsEvTimeT0AC<pidflags::TOFFlags>,
                  pidflags::IsEvTimeTOFT0AC<pidflags::TOFFlags>,
                  spectra::TOFNSigmaKa<spectra::TOFNSigmaStoreKa>,
                  spectra::TPCNSigmaKa<spectra::TPCNSigmaStoreKa>);
} // namespace o2::aod

struct MultCodes {
  static constexpr int kNoMultiplicity = 0;
  static constexpr int kMultFV0M = 1;
  static constexpr int kMultFT0M = 2;
  static constexpr int kMultFDDM = 3;
  static constexpr int kMultTracklets = 4;
  static constexpr int kMultTPC = 5;
  static constexpr int kMultNTracksPV = 6;
  static constexpr int kMultNTracksPVeta1 = 7;
  static constexpr int kCentralityFT0C = 8;
  static constexpr int kCentralityFT0M = 9;
  static constexpr int kCentralityFV0A = 10;
  static constexpr int kNMults = 10;
};

#endif // PWGLF_DATAMODEL_KAONANALYSIS_H_
