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

/// \file CandidateReconstructionTables.h
/// \brief Definitions of tables produced by candidate reconstruction workflows
///
/// \author Gian Michele Innocenti <gian.michele.innocenti@cern.ch>, CERN
/// \author Vít Kučera <vit.kucera@cern.ch>, CERN

#ifndef PWGHF_DATAMODEL_CANDIDATERECONSTRUCTIONTABLES_H_
#define PWGHF_DATAMODEL_CANDIDATERECONSTRUCTIONTABLES_H_

#include "PWGHF/DataModel/TrackIndexSkimmingTables.h"
#include "PWGHF/Utils/utilsPid.h"
//
#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "ALICE3/DataModel/ECAL.h"
#include "Common/Core/RecoDecay.h"

#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>

#include <array>
#include <cstdint>

namespace o2::aod
{
namespace pid_tpc_tof_static_full
{
// Combined TPC and TOF NSigma
DECLARE_SOA_COLUMN(TpcTofNSigmaEl, tpcTofNSigmaEl, float); //! Combined NSigma separation with the TPC & TOF detectors for electron
DECLARE_SOA_COLUMN(TpcTofNSigmaMu, tpcTofNSigmaMu, float); //! Combined NSigma separation with the TPC & TOF detectors for muon
DECLARE_SOA_COLUMN(TpcTofNSigmaPi, tpcTofNSigmaPi, float); //! Combined NSigma separation with the TPC & TOF detectors for pion
DECLARE_SOA_COLUMN(TpcTofNSigmaKa, tpcTofNSigmaKa, float); //! Combined NSigma separation with the TPC & TOF detectors for kaon
DECLARE_SOA_COLUMN(TpcTofNSigmaPr, tpcTofNSigmaPr, float); //! Combined NSigma separation with the TPC & TOF detectors for proton
DECLARE_SOA_COLUMN(TpcTofNSigmaDe, tpcTofNSigmaDe, float); //! Combined NSigma separation with the TPC & TOF detectors for deuteron
DECLARE_SOA_COLUMN(TpcTofNSigmaTr, tpcTofNSigmaTr, float); //! Combined NSigma separation with the TPC & TOF detectors for triton
DECLARE_SOA_COLUMN(TpcTofNSigmaHe, tpcTofNSigmaHe, float); //! Combined NSigma separation with the TPC & TOF detectors for helium
} // namespace pid_tpc_tof_static_full

namespace pid_tpc_tof_static_tiny
{
// Combined TPC and TOF NSigma
DECLARE_SOA_COLUMN(TpcTofNSigmaEl, tpcTofNSigmaEl, float); //! Combined NSigma separation with the TPC & TOF detectors for electron
DECLARE_SOA_COLUMN(TpcTofNSigmaMu, tpcTofNSigmaMu, float); //! Combined NSigma separation with the TPC & TOF detectors for muon
DECLARE_SOA_COLUMN(TpcTofNSigmaPi, tpcTofNSigmaPi, float); //! Combined NSigma separation with the TPC & TOF detectors for pion
DECLARE_SOA_COLUMN(TpcTofNSigmaKa, tpcTofNSigmaKa, float); //! Combined NSigma separation with the TPC & TOF detectors for kaon
DECLARE_SOA_COLUMN(TpcTofNSigmaPr, tpcTofNSigmaPr, float); //! Combined NSigma separation with the TPC & TOF detectors for proton
DECLARE_SOA_COLUMN(TpcTofNSigmaDe, tpcTofNSigmaDe, float); //! Combined NSigma separation with the TPC & TOF detectors for deuteron
DECLARE_SOA_COLUMN(TpcTofNSigmaTr, tpcTofNSigmaTr, float); //! Combined NSigma separation with the TPC & TOF detectors for triton
DECLARE_SOA_COLUMN(TpcTofNSigmaHe, tpcTofNSigmaHe, float); //! Combined NSigma separation with the TPC & TOF detectors for helium
} // namespace pid_tpc_tof_static_tiny

// Extension of per particle tables
DECLARE_SOA_TABLE(PidTpcTofFullEl, "AOD", "PIDTPCTOFFULLEL", //! Table of the TPC & TOF Combined NSigma for electron
                  pid_tpc_tof_static_full::TpcTofNSigmaEl);
DECLARE_SOA_TABLE(PidTpcTofFullMu, "AOD", "PIDTPCTOFFULLMU", //! Table of the TPC & TOF Combined NSigma for muon
                  pid_tpc_tof_static_full::TpcTofNSigmaMu);
DECLARE_SOA_TABLE(PidTpcTofFullPi, "AOD", "PIDTPCTOFFULLPI", //! Table of the TPC & TOF Combined NSigma for pion
                  pid_tpc_tof_static_full::TpcTofNSigmaPi);
DECLARE_SOA_TABLE(PidTpcTofFullKa, "AOD", "PIDTPCTOFFULLKA", //! Table of the TPC & TOF Combined NSigma for kaon
                  pid_tpc_tof_static_full::TpcTofNSigmaKa);
DECLARE_SOA_TABLE(PidTpcTofFullPr, "AOD", "PIDTPCTOFFULLPR", //! Table of the TPC & TOF Combined NSigma for proton
                  pid_tpc_tof_static_full::TpcTofNSigmaPr);
DECLARE_SOA_TABLE(PidTpcTofFullDe, "AOD", "PIDTPCTOFFULLDe", //! Table of the TPC & TOF Combined NSigma for deuteron
                  pid_tpc_tof_static_full::TpcTofNSigmaDe);
DECLARE_SOA_TABLE(PidTpcTofFullTr, "AOD", "PIDTPCTOFFULLTr", //! Table of the TPC & TOF Combined NSigma for triton
                  pid_tpc_tof_static_full::TpcTofNSigmaTr);
DECLARE_SOA_TABLE(PidTpcTofFullHe, "AOD", "PIDTPCTOFFULLHe", //! Table of the TPC & TOF Combined NSigma for helium
                  pid_tpc_tof_static_full::TpcTofNSigmaHe);
// Extension of per particle tables
DECLARE_SOA_TABLE(PidTpcTofTinyEl, "AOD", "PIDTPCTOFTINYEL", //! Table of the TPC & TOF Combined NSigma for electron
                  pid_tpc_tof_static_tiny::TpcTofNSigmaEl);
DECLARE_SOA_TABLE(PidTpcTofTinyMu, "AOD", "PIDTPCTOFTINYMU", //! Table of the TPC & TOF Combined NSigma for muon
                  pid_tpc_tof_static_tiny::TpcTofNSigmaMu);
DECLARE_SOA_TABLE(PidTpcTofTinyPi, "AOD", "PIDTPCTOFTINYPI", //! Table of the TPC & TOF Combined NSigma for pion
                  pid_tpc_tof_static_tiny::TpcTofNSigmaPi);
DECLARE_SOA_TABLE(PidTpcTofTinyKa, "AOD", "PIDTPCTOFTINYKA", //! Table of the TPC & TOF Combined NSigma for kaon
                  pid_tpc_tof_static_tiny::TpcTofNSigmaKa);
DECLARE_SOA_TABLE(PidTpcTofTinyPr, "AOD", "PIDTPCTOFTINYPR", //! Table of the TPC & TOF Combined NSigma for proton
                  pid_tpc_tof_static_tiny::TpcTofNSigmaPr);
DECLARE_SOA_TABLE(PidTpcTofTinyDe, "AOD", "PIDTPCTOFTINYDE", //! Table of the TPC & TOF Combined NSigma for deuteron
                  pid_tpc_tof_static_tiny::TpcTofNSigmaDe);
DECLARE_SOA_TABLE(PidTpcTofTinyTr, "AOD", "PIDTPCTOFTINYTr", //! Table of the TPC & TOF Combined NSigma for triton
                  pid_tpc_tof_static_tiny::TpcTofNSigmaTr);
DECLARE_SOA_TABLE(PidTpcTofTinyHe, "AOD", "PIDTPCTOFTINYHe", //! Table of the TPC & TOF Combined NSigma for helium
                  pid_tpc_tof_static_tiny::TpcTofNSigmaHe);
// general decay properties
namespace hf_cand
{
// collision properties
DECLARE_SOA_INDEX_COLUMN(Collision, collision); //!
// secondary vertex
DECLARE_SOA_COLUMN(XSecondaryVertex, xSecondaryVertex, float); //!
DECLARE_SOA_COLUMN(YSecondaryVertex, ySecondaryVertex, float); //!
DECLARE_SOA_COLUMN(ZSecondaryVertex, zSecondaryVertex, float); //!
DECLARE_SOA_DYNAMIC_COLUMN(RSecondaryVertex, rSecondaryVertex, //!
                           [](float xVtxS, float yVtxS) -> float { return RecoDecay::sqrtSumOfSquares(xVtxS, yVtxS); });
DECLARE_SOA_COLUMN(Chi2PCA, chi2PCA, float); //! sum of (non-weighted) distances of the secondary vertex to its prongs
// prong properties
DECLARE_SOA_COLUMN(PxProng0, pxProng0, float); //!
DECLARE_SOA_COLUMN(PyProng0, pyProng0, float); //!
DECLARE_SOA_COLUMN(PzProng0, pzProng0, float); //!
DECLARE_SOA_DYNAMIC_COLUMN(PtProng0, ptProng0, //!
                           [](float px, float py) -> float { return RecoDecay::pt(px, py); });
DECLARE_SOA_DYNAMIC_COLUMN(Pt2Prong0, pt2Prong0, //!
                           [](float px, float py) -> float { return RecoDecay::pt2(px, py); });
DECLARE_SOA_DYNAMIC_COLUMN(PVectorProng0, pVectorProng0, //!
                           [](float px, float py, float pz) -> std::array<float, 3> { return std::array{px, py, pz}; });
DECLARE_SOA_COLUMN(ImpactParameter0, impactParameter0, float);                     //!
DECLARE_SOA_COLUMN(ErrorImpactParameter0, errorImpactParameter0, float);           //!
DECLARE_SOA_DYNAMIC_COLUMN(ImpactParameterNormalised0, impactParameterNormalised0, //!
                           [](float dca, float err) -> float { return dca / err; });
DECLARE_SOA_COLUMN(ImpactParameterZ0, impactParameterZ0, float);                     //!
DECLARE_SOA_COLUMN(ErrorImpactParameterZ0, errorImpactParameterZ0, float);           //!
DECLARE_SOA_DYNAMIC_COLUMN(ImpactParameterZNormalised0, impactParameterZNormalised0, //!
                           [](float dca, float err) -> float { return dca / err; });
DECLARE_SOA_COLUMN(PxProng1, pxProng1, float); //!
DECLARE_SOA_COLUMN(PyProng1, pyProng1, float); //!
DECLARE_SOA_COLUMN(PzProng1, pzProng1, float); //!
DECLARE_SOA_DYNAMIC_COLUMN(PtProng1, ptProng1, //!
                           [](float px, float py) -> float { return RecoDecay::pt(px, py); });
DECLARE_SOA_DYNAMIC_COLUMN(Pt2Prong1, pt2Prong1, //!
                           [](float px, float py) -> float { return RecoDecay::pt2(px, py); });
DECLARE_SOA_DYNAMIC_COLUMN(PVectorProng1, pVectorProng1, //!
                           [](float px, float py, float pz) -> std::array<float, 3> { return std::array{px, py, pz}; });
DECLARE_SOA_COLUMN(ImpactParameter1, impactParameter1, float);                     //!
DECLARE_SOA_COLUMN(ErrorImpactParameter1, errorImpactParameter1, float);           //!
DECLARE_SOA_DYNAMIC_COLUMN(ImpactParameterNormalised1, impactParameterNormalised1, //!
                           [](float dca, float err) -> float { return dca / err; });
DECLARE_SOA_COLUMN(ImpactParameterZ1, impactParameterZ1, float);                     //!
DECLARE_SOA_COLUMN(ErrorImpactParameterZ1, errorImpactParameterZ1, float);           //!
DECLARE_SOA_DYNAMIC_COLUMN(ImpactParameterZNormalised1, impactParameterZNormalised1, //!
                           [](float dca, float err) -> float { return dca / err; });
DECLARE_SOA_COLUMN(PxProng2, pxProng2, float); //!
DECLARE_SOA_COLUMN(PyProng2, pyProng2, float); //!
DECLARE_SOA_COLUMN(PzProng2, pzProng2, float); //!
DECLARE_SOA_DYNAMIC_COLUMN(PtProng2, ptProng2, //!
                           [](float px, float py) -> float { return RecoDecay::pt(px, py); });
DECLARE_SOA_DYNAMIC_COLUMN(Pt2Prong2, pt2Prong2, //!
                           [](float px, float py) -> float { return RecoDecay::pt2(px, py); });
DECLARE_SOA_DYNAMIC_COLUMN(PVectorProng2, pVectorProng2, //!
                           [](float px, float py, float pz) -> std::array<float, 3> { return std::array{px, py, pz}; });
DECLARE_SOA_COLUMN(ImpactParameter2, impactParameter2, float);                     //!
DECLARE_SOA_COLUMN(ErrorImpactParameter2, errorImpactParameter2, float);           //!
DECLARE_SOA_DYNAMIC_COLUMN(ImpactParameterNormalised2, impactParameterNormalised2, //!
                           [](float dca, float err) -> float { return dca / err; });
DECLARE_SOA_COLUMN(ImpactParameterZ2, impactParameterZ2, float);                     //!
DECLARE_SOA_COLUMN(ErrorImpactParameterZ2, errorImpactParameterZ2, float);           //!
DECLARE_SOA_DYNAMIC_COLUMN(ImpactParameterZNormalised2, impactParameterZNormalised2, //!
                           [](float dca, float err) -> float { return dca / err; });
DECLARE_SOA_COLUMN(PxProng3, pxProng3, float); //!
DECLARE_SOA_COLUMN(PyProng3, pyProng3, float); //!
DECLARE_SOA_COLUMN(PzProng3, pzProng3, float); //!
DECLARE_SOA_DYNAMIC_COLUMN(PtProng3, ptProng3, //!
                           [](float px, float py) -> float { return RecoDecay::pt(px, py); });
DECLARE_SOA_DYNAMIC_COLUMN(Pt2Prong3, pt2Prong3, //!
                           [](float px, float py) -> float { return RecoDecay::pt2(px, py); });
DECLARE_SOA_DYNAMIC_COLUMN(PVectorProng3, pVectorProng3, //!
                           [](float px, float py, float pz) -> std::array<float, 3> { return std::array{px, py, pz}; });
DECLARE_SOA_COLUMN(ImpactParameter3, impactParameter3, float);                     //!
DECLARE_SOA_COLUMN(ErrorImpactParameter3, errorImpactParameter3, float);           //!
DECLARE_SOA_DYNAMIC_COLUMN(ImpactParameterNormalised3, impactParameterNormalised3, //!
                           [](float dca, float err) -> float { return dca / err; });
DECLARE_SOA_COLUMN(ImpactParameterZ3, impactParameterZ3, float);                     //!
DECLARE_SOA_COLUMN(ErrorImpactParameterZ3, errorImpactParameterZ3, float);           //!
DECLARE_SOA_DYNAMIC_COLUMN(ImpactParameterZNormalised3, impactParameterZNormalised3, //!
                           [](float dca, float err) -> float { return dca / err; });
DECLARE_SOA_COLUMN(NProngsContributorsPV, nProngsContributorsPV, uint8_t);           //! number of prongs contributing to the primary-vertex reconstruction
DECLARE_SOA_COLUMN(BitmapProngsContributorsPV, bitmapProngsContributorsPV, uint8_t); //! bitmap with booleans indicating prongs contributing to the primary-vertex reconstruction
/// prong PID nsigma
DECLARE_SOA_COLUMN(NSigTpcPi0, nSigTpcPi0, float);           //! TPC nSigma for pion hypothesis - prong 0
DECLARE_SOA_COLUMN(NSigTpcPi1, nSigTpcPi1, float);           //! TPC nSigma for pion hypothesis - prong 1
DECLARE_SOA_COLUMN(NSigTpcPi2, nSigTpcPi2, float);           //! TPC nSigma for pion hypothesis - prong 2
DECLARE_SOA_COLUMN(NSigTpcKa0, nSigTpcKa0, float);           //! TPC nSigma for kaon hypothesis - prong 0
DECLARE_SOA_COLUMN(NSigTpcKa1, nSigTpcKa1, float);           //! TPC nSigma for kaon hypothesis - prong 1
DECLARE_SOA_COLUMN(NSigTpcKa2, nSigTpcKa2, float);           //! TPC nSigma for kaon hypothesis - prong 2
DECLARE_SOA_COLUMN(NSigTpcPr0, nSigTpcPr0, float);           //! TPC nSigma for proton hypothesis - prong 0
DECLARE_SOA_COLUMN(NSigTpcPr1, nSigTpcPr1, float);           //! TPC nSigma for proton hypothesis - prong 1
DECLARE_SOA_COLUMN(NSigTpcPr2, nSigTpcPr2, float);           //! TPC nSigma for proton hypothesis - prong 2
DECLARE_SOA_COLUMN(NSigTpcDe0, nSigTpcDe0, float);           //! TPC nSigma for deuteron hypothesis - prong 0
DECLARE_SOA_COLUMN(NSigTpcDe1, nSigTpcDe1, float);           //! TPC nSigma for deuteron hypothesis - prong 1
DECLARE_SOA_COLUMN(NSigTpcDe2, nSigTpcDe2, float);           //! TPC nSigma for deuteron hypothesis - prong 2
DECLARE_SOA_COLUMN(NSigTpcTr0, nSigTpcTr0, float);           //! TPC nSigma for triton hypothesis - prong 0
DECLARE_SOA_COLUMN(NSigTpcTr1, nSigTpcTr1, float);           //! TPC nSigma for triton hypothesis - prong 1
DECLARE_SOA_COLUMN(NSigTpcTr2, nSigTpcTr2, float);           //! TPC nSigma for triton hypothesis - prong 2
DECLARE_SOA_COLUMN(NSigTpcHe0, nSigTpcHe0, float);           //! TPC nSigma for helium hypothesis - prong 0
DECLARE_SOA_COLUMN(NSigTpcHe1, nSigTpcHe1, float);           //! TPC nSigma for helium hypothesis - prong 1
DECLARE_SOA_COLUMN(NSigTpcHe2, nSigTpcHe2, float);           //! TPC nSigma for helium hypothesis - prong 2
DECLARE_SOA_COLUMN(NSigTofPi0, nSigTofPi0, float);           //! TOF nSigma for pion hypothesis - prong 0
DECLARE_SOA_COLUMN(NSigTofPi1, nSigTofPi1, float);           //! TOF nSigma for pion hypothesis - prong 1
DECLARE_SOA_COLUMN(NSigTofPi2, nSigTofPi2, float);           //! TOF nSigma for pion hypothesis - prong 2
DECLARE_SOA_COLUMN(NSigTofKa0, nSigTofKa0, float);           //! TOF nSigma for kaon hypothesis - prong 0
DECLARE_SOA_COLUMN(NSigTofKa1, nSigTofKa1, float);           //! TOF nSigma for kaon hypothesis - prong 1
DECLARE_SOA_COLUMN(NSigTofKa2, nSigTofKa2, float);           //! TOF nSigma for kaon hypothesis - prong 2
DECLARE_SOA_COLUMN(NSigTofPr0, nSigTofPr0, float);           //! TOF nSigma for proton hypothesis - prong 0
DECLARE_SOA_COLUMN(NSigTofPr1, nSigTofPr1, float);           //! TOF nSigma for proton hypothesis - prong 1
DECLARE_SOA_COLUMN(NSigTofPr2, nSigTofPr2, float);           //! TOF nSigma for proton hypothesis - prong 2
DECLARE_SOA_COLUMN(NSigTofDe0, nSigTofDe0, float);           //! TOF nSigma for deuteron hypothesis - prong 0
DECLARE_SOA_COLUMN(NSigTofDe1, nSigTofDe1, float);           //! TOF nSigma for deuteron hypothesis - prong 1
DECLARE_SOA_COLUMN(NSigTofDe2, nSigTofDe2, float);           //! TOF nSigma for deuteron hypothesis - prong 2
DECLARE_SOA_COLUMN(NSigTofTr0, nSigTofTr0, float);           //! TOF nSigma for triton hypothesis - prong 0
DECLARE_SOA_COLUMN(NSigTofTr1, nSigTofTr1, float);           //! TOF nSigma for triton hypothesis - prong 1
DECLARE_SOA_COLUMN(NSigTofTr2, nSigTofTr2, float);           //! TOF nSigma for triton hypothesis - prong 2
DECLARE_SOA_COLUMN(NSigTofHe0, nSigTofHe0, float);           //! TOF nSigma for helium hypothesis - prong 0
DECLARE_SOA_COLUMN(NSigTofHe1, nSigTofHe1, float);           //! TOF nSigma for helium hypothesis - prong 1
DECLARE_SOA_COLUMN(NSigTofHe2, nSigTofHe2, float);           //! TOF nSigma for helium hypothesis - prong 2
DECLARE_SOA_DYNAMIC_COLUMN(TpcTofNSigmaPi0, tpcTofNSigmaPi0, //! Combined NSigma separation with the TPC & TOF detectors for pion - prong 0
                           [](float tpcNSigmaPi0, float tofNSigmaPi0) -> float { return pid_tpc_tof_utils::combineNSigma<false /*tiny*/>(tpcNSigmaPi0, tofNSigmaPi0); });
DECLARE_SOA_DYNAMIC_COLUMN(TpcTofNSigmaPi1, tpcTofNSigmaPi1, //! Combined NSigma separation with the TPC & TOF detectors for pion - prong 1
                           [](float tpcNSigmaPi1, float tofNSigmaPi1) -> float { return pid_tpc_tof_utils::combineNSigma<false /*tiny*/>(tpcNSigmaPi1, tofNSigmaPi1); });
DECLARE_SOA_DYNAMIC_COLUMN(TpcTofNSigmaPi2, tpcTofNSigmaPi2, //! Combined NSigma separation with the TPC & TOF detectors for pion - prong 2
                           [](float tpcNSigmaPi2, float tofNSigmaPi2) -> float { return pid_tpc_tof_utils::combineNSigma<false /*tiny*/>(tpcNSigmaPi2, tofNSigmaPi2); });
DECLARE_SOA_DYNAMIC_COLUMN(TpcTofNSigmaKa0, tpcTofNSigmaKa0, //! Combined NSigma separation with the TPC & TOF detectors for kaon - prong 0
                           [](float tpcNSigmaKa0, float tofNSigmaKa0) -> float { return pid_tpc_tof_utils::combineNSigma<false /*tiny*/>(tpcNSigmaKa0, tofNSigmaKa0); });
DECLARE_SOA_DYNAMIC_COLUMN(TpcTofNSigmaKa1, tpcTofNSigmaKa1, //! Combined NSigma separation with the TPC & TOF detectors for kaon - prong 1
                           [](float tpcNSigmaKa1, float tofNSigmaKa1) -> float { return pid_tpc_tof_utils::combineNSigma<false /*tiny*/>(tpcNSigmaKa1, tofNSigmaKa1); });
DECLARE_SOA_DYNAMIC_COLUMN(TpcTofNSigmaKa2, tpcTofNSigmaKa2, //! Combined NSigma separation with the TPC & TOF detectors for kaon - prong 2
                           [](float tpcNSigmaKa2, float tofNSigmaKa2) -> float { return pid_tpc_tof_utils::combineNSigma<false /*tiny*/>(tpcNSigmaKa2, tofNSigmaKa2); });
DECLARE_SOA_DYNAMIC_COLUMN(TpcTofNSigmaPr0, tpcTofNSigmaPr0, //! Combined NSigma separation with the TPC & TOF detectors for proton - prong 0
                           [](float tpcNSigmaPr0, float tofNSigmaPr0) -> float { return pid_tpc_tof_utils::combineNSigma<false /*tiny*/>(tpcNSigmaPr0, tofNSigmaPr0); });
DECLARE_SOA_DYNAMIC_COLUMN(TpcTofNSigmaPr1, tpcTofNSigmaPr1, //! Combined NSigma separation with the TPC & TOF detectors for proton - prong 1
                           [](float tpcNSigmaPr1, float tofNSigmaPr1) -> float { return pid_tpc_tof_utils::combineNSigma<false /*tiny*/>(tpcNSigmaPr1, tofNSigmaPr1); });
DECLARE_SOA_DYNAMIC_COLUMN(TpcTofNSigmaPr2, tpcTofNSigmaPr2, //! Combined NSigma separation with the TPC & TOF detectors for proton - prong 2
                           [](float tpcNSigmaPr2, float tofNSigmaPr2) -> float { return pid_tpc_tof_utils::combineNSigma<false /*tiny*/>(tpcNSigmaPr2, tofNSigmaPr2); });
DECLARE_SOA_DYNAMIC_COLUMN(TpcTofNSigmaDe0, tpcTofNSigmaDe0, //! Combined NSigma separation with the TPC & TOF detectors for deuteron - prong 0
                           [](float tpcNSigmaDe0, float tofNSigmaDe0) -> float { return pid_tpc_tof_utils::combineNSigma<false /*tiny*/>(tpcNSigmaDe0, tofNSigmaDe0); });
DECLARE_SOA_DYNAMIC_COLUMN(TpcTofNSigmaDe1, tpcTofNSigmaDe1, //! Combined NSigma separation with the TPC & TOF detectors for deuteron - prong 1
                           [](float tpcNSigmaDe1, float tofNSigmaDe1) -> float { return pid_tpc_tof_utils::combineNSigma<false /*tiny*/>(tpcNSigmaDe1, tofNSigmaDe1); });
DECLARE_SOA_DYNAMIC_COLUMN(TpcTofNSigmaDe2, tpcTofNSigmaDe2, //! Combined NSigma separation with the TPC & TOF detectors for deuteron - prong 2
                           [](float tpcNSigmaDe2, float tofNSigmaDe2) -> float { return pid_tpc_tof_utils::combineNSigma<false /*tiny*/>(tpcNSigmaDe2, tofNSigmaDe2); });
DECLARE_SOA_DYNAMIC_COLUMN(TpcTofNSigmaTr0, tpcTofNSigmaTr0, //! Combined NSigma separation with the TPC & TOF detectors for triton - prong 0
                           [](float tpcNSigmaTr0, float tofNSigmaTr0) -> float { return pid_tpc_tof_utils::combineNSigma<false /*tiny*/>(tpcNSigmaTr0, tofNSigmaTr0); });
DECLARE_SOA_DYNAMIC_COLUMN(TpcTofNSigmaTr1, tpcTofNSigmaTr1, //! Combined NSigma separation with the TPC & TOF detectors for triton - prong 1
                           [](float tpcNSigmaTr1, float tofNSigmaTr1) -> float { return pid_tpc_tof_utils::combineNSigma<false /*tiny*/>(tpcNSigmaTr1, tofNSigmaTr1); });
DECLARE_SOA_DYNAMIC_COLUMN(TpcTofNSigmaTr2, tpcTofNSigmaTr2, //! Combined NSigma separation with the TPC & TOF detectors for triton - prong 2
                           [](float tpcNSigmaTr2, float tofNSigmaTr2) -> float { return pid_tpc_tof_utils::combineNSigma<false /*tiny*/>(tpcNSigmaTr2, tofNSigmaTr2); });
DECLARE_SOA_DYNAMIC_COLUMN(TpcTofNSigmaHe0, tpcTofNSigmaHe0, //! Combined NSigma separation with the TPC & TOF detectors for helium - prong 0
                           [](float tpcNSigmaHe0, float tofNSigmaHe0) -> float { return pid_tpc_tof_utils::combineNSigma<false /*tiny*/>(tpcNSigmaHe0, tofNSigmaHe0); });
DECLARE_SOA_DYNAMIC_COLUMN(TpcTofNSigmaHe1, tpcTofNSigmaHe1, //! Combined NSigma separation with the TPC & TOF detectors for helium - prong 1
                           [](float tpcNSigmaHe1, float tofNSigmaHe1) -> float { return pid_tpc_tof_utils::combineNSigma<false /*tiny*/>(tpcNSigmaHe1, tofNSigmaHe1); });
DECLARE_SOA_DYNAMIC_COLUMN(TpcTofNSigmaHe2, tpcTofNSigmaHe2, //! Combined NSigma separation with the TPC & TOF detectors for helium - prong 2
                           [](float tpcNSigmaHe2, float tofNSigmaHe2) -> float { return pid_tpc_tof_utils::combineNSigma<false /*tiny*/>(tpcNSigmaHe2, tofNSigmaHe2); });
// tiny (binned) option
DECLARE_SOA_DYNAMIC_COLUMN(TpcTofNSigmaTinyPi0, tpcTofNSigmaTinyPi0, //! Combined NSigma separation with the TPC & TOF detectors for pion - prong 0
                           [](float tpcNSigmaPi0, float tofNSigmaPi0) -> float { return pid_tpc_tof_utils::combineNSigma<true /*tiny*/>(tpcNSigmaPi0, tofNSigmaPi0); });
DECLARE_SOA_DYNAMIC_COLUMN(TpcTofNSigmaTinyPi1, tpcTofNSigmaTinyPi1, //! Combined NSigma separation with the TPC & TOF detectors for pion - prong 1
                           [](float tpcNSigmaPi1, float tofNSigmaPi1) -> float { return pid_tpc_tof_utils::combineNSigma<true /*tiny*/>(tpcNSigmaPi1, tofNSigmaPi1); });
DECLARE_SOA_DYNAMIC_COLUMN(TpcTofNSigmaTinyPi2, tpcTofNSigmaTinyPi2, //! Combined NSigma separation with the TPC & TOF detectors for pion - prong 2
                           [](float tpcNSigmaPi2, float tofNSigmaPi2) -> float { return pid_tpc_tof_utils::combineNSigma<true /*tiny*/>(tpcNSigmaPi2, tofNSigmaPi2); });
DECLARE_SOA_DYNAMIC_COLUMN(TpcTofNSigmaTinyKa0, tpcTofNSigmaTinyKa0, //! Combined NSigma separation with the TPC & TOF detectors for kaon - prong 0
                           [](float tpcNSigmaKa0, float tofNSigmaKa0) -> float { return pid_tpc_tof_utils::combineNSigma<true /*tiny*/>(tpcNSigmaKa0, tofNSigmaKa0); });
DECLARE_SOA_DYNAMIC_COLUMN(TpcTofNSigmaTinyKa1, tpcTofNSigmaTinyKa1, //! Combined NSigma separation with the TPC & TOF detectors for kaon - prong 1
                           [](float tpcNSigmaKa1, float tofNSigmaKa1) -> float { return pid_tpc_tof_utils::combineNSigma<true /*tiny*/>(tpcNSigmaKa1, tofNSigmaKa1); });
DECLARE_SOA_DYNAMIC_COLUMN(TpcTofNSigmaTinyKa2, tpcTofNSigmaTinyKa2, //! Combined NSigma separation with the TPC & TOF detectors for kaon - prong 2
                           [](float tpcNSigmaKa2, float tofNSigmaKa2) -> float { return pid_tpc_tof_utils::combineNSigma<true /*tiny*/>(tpcNSigmaKa2, tofNSigmaKa2); });
DECLARE_SOA_DYNAMIC_COLUMN(TpcTofNSigmaTinyPr0, tpcTofNSigmaTinyPr0, //! Combined NSigma separation with the TPC & TOF detectors for proton - prong 0
                           [](float tpcNSigmaPr0, float tofNSigmaPr0) -> float { return pid_tpc_tof_utils::combineNSigma<true /*tiny*/>(tpcNSigmaPr0, tofNSigmaPr0); });
DECLARE_SOA_DYNAMIC_COLUMN(TpcTofNSigmaTinyPr1, tpcTofNSigmaTinyPr1, //! Combined NSigma separation with the TPC & TOF detectors for proton - prong 1
                           [](float tpcNSigmaPr1, float tofNSigmaPr1) -> float { return pid_tpc_tof_utils::combineNSigma<true /*tiny*/>(tpcNSigmaPr1, tofNSigmaPr1); });
DECLARE_SOA_DYNAMIC_COLUMN(TpcTofNSigmaTinyPr2, tpcTofNSigmaTinyPr2, //! Combined NSigma separation with the TPC & TOF detectors for proton - prong 2
                           [](float tpcNSigmaPr2, float tofNSigmaPr2) -> float { return pid_tpc_tof_utils::combineNSigma<true /*tiny*/>(tpcNSigmaPr2, tofNSigmaPr2); });

// candidate properties
DECLARE_SOA_DYNAMIC_COLUMN(Pt, pt, //!
                           [](float px, float py) -> float { return RecoDecay::pt(px, py); });
DECLARE_SOA_DYNAMIC_COLUMN(Pt2, pt2, //!
                           [](float px, float py) -> float { return RecoDecay::pt2(px, py); });
DECLARE_SOA_DYNAMIC_COLUMN(P, p, //!
                           [](float px, float py, float pz) -> float { return RecoDecay::p(px, py, pz); });
DECLARE_SOA_DYNAMIC_COLUMN(P2, p2, //!
                           [](float px, float py, float pz) -> float { return RecoDecay::p2(px, py, pz); });
DECLARE_SOA_DYNAMIC_COLUMN(PVector, pVector, //!
                           [](float px, float py, float pz) -> std::array<float, 3> { return std::array{px, py, pz}; });
DECLARE_SOA_DYNAMIC_COLUMN(Eta, eta, //!
                           [](float px, float py, float pz) -> float { return RecoDecay::eta(std::array{px, py, pz}); });
DECLARE_SOA_DYNAMIC_COLUMN(Phi, phi, //!
                           [](float px, float py) -> float { return RecoDecay::phi(px, py); });
DECLARE_SOA_DYNAMIC_COLUMN(Y, y, //!
                           [](float px, float py, float pz, double m) -> float { return RecoDecay::y(std::array{px, py, pz}, m); });
DECLARE_SOA_DYNAMIC_COLUMN(E, e, //!
                           [](float px, float py, float pz, double m) -> float { return RecoDecay::e(px, py, pz, m); });
DECLARE_SOA_DYNAMIC_COLUMN(E2, e2, //!
                           [](float px, float py, float pz, double m) -> float { return RecoDecay::e2(px, py, pz, m); });
DECLARE_SOA_DYNAMIC_COLUMN(DecayLength, decayLength, //!
                           [](float xVtxP, float yVtxP, float zVtxP, float xVtxS, float yVtxS, float zVtxS) -> float { return RecoDecay::distance(std::array{xVtxP, yVtxP, zVtxP}, std::array{xVtxS, yVtxS, zVtxS}); });
DECLARE_SOA_DYNAMIC_COLUMN(DecayLengthXY, decayLengthXY, //!
                           [](float xVtxP, float yVtxP, float xVtxS, float yVtxS) -> float { return RecoDecay::distanceXY(std::array{xVtxP, yVtxP}, std::array{xVtxS, yVtxS}); });
DECLARE_SOA_DYNAMIC_COLUMN(DecayLengthNormalised, decayLengthNormalised, //!
                           [](float xVtxP, float yVtxP, float zVtxP, float xVtxS, float yVtxS, float zVtxS, float err) -> float { return RecoDecay::distance(std::array{xVtxP, yVtxP, zVtxP}, std::array{xVtxS, yVtxS, zVtxS}) / err; });
DECLARE_SOA_DYNAMIC_COLUMN(DecayLengthXYNormalised, decayLengthXYNormalised, //!
                           [](float xVtxP, float yVtxP, float xVtxS, float yVtxS, float err) -> float { return RecoDecay::distanceXY(std::array{xVtxP, yVtxP}, std::array{xVtxS, yVtxS}) / err; });
DECLARE_SOA_COLUMN(ErrorDecayLength, errorDecayLength, float);     //!
DECLARE_SOA_COLUMN(ErrorDecayLengthXY, errorDecayLengthXY, float); //!
DECLARE_SOA_DYNAMIC_COLUMN(Cpa, cpa,                               //!
                           [](float xVtxP, float yVtxP, float zVtxP, float xVtxS, float yVtxS, float zVtxS, float px, float py, float pz) -> float { return RecoDecay::cpa(std::array{xVtxP, yVtxP, zVtxP}, std::array{xVtxS, yVtxS, zVtxS}, std::array{px, py, pz}); });
DECLARE_SOA_DYNAMIC_COLUMN(CpaXY, cpaXY, //!
                           [](float xVtxP, float yVtxP, float xVtxS, float yVtxS, float px, float py) -> float { return RecoDecay::cpaXY(std::array{xVtxP, yVtxP}, std::array{xVtxS, yVtxS}, std::array{px, py}); });
DECLARE_SOA_DYNAMIC_COLUMN(Ct, ct, //!
                           [](float xVtxP, float yVtxP, float zVtxP, float xVtxS, float yVtxS, float zVtxS, float px, float py, float pz, double m) -> float { return RecoDecay::ct(std::array{px, py, pz}, RecoDecay::distance(std::array{xVtxP, yVtxP, zVtxP}, std::array{xVtxS, yVtxS, zVtxS}), m); });
DECLARE_SOA_DYNAMIC_COLUMN(ImpactParameterXY, impactParameterXY, //!
                           [](float xVtxP, float yVtxP, float zVtxP, float xVtxS, float yVtxS, float zVtxS, float px, float py, float pz) -> float { return RecoDecay::impParXY(std::array{xVtxP, yVtxP, zVtxP}, std::array{xVtxS, yVtxS, zVtxS}, std::array{px, py, pz}); });
DECLARE_SOA_COLUMN(KfTopolChi2OverNdf, kfTopolChi2OverNdf, float); //! chi2overndf of the KFParticle topological constraint
// B-hadron mother information
DECLARE_SOA_COLUMN(PtBhadMotherPart, ptBhadMotherPart, float); //! pt of the first B-hadron mother particle (only in case of non-prompt)
DECLARE_SOA_COLUMN(PdgBhadMotherPart, pdgBhadMotherPart, int); //! pdg of the first B-hadron mother particle (only in case of non-prompt)
DECLARE_SOA_COLUMN(IdxBhadMotherPart, idxBhadMotherPart, int); //! index of the first B-hadron mother particle (only in case of non-prompt)
// Kink topology and material interaction mc flags
DECLARE_SOA_COLUMN(NTracksDecayed, nTracksDecayed, int8_t);                       //! number of tracks matched with kinked decay topology
DECLARE_SOA_COLUMN(NInteractionsWithMaterial, nInteractionsWithMaterial, int8_t); //! number of tracks matched after interaction with material

// method of secondary-vertex reconstruction
enum VertexerType { DCAFitter = 0,
                    KfParticle };
} // namespace hf_cand

// specific 2-prong decay properties
namespace hf_cand_2prong
{
DECLARE_SOA_EXPRESSION_COLUMN(Px, px, //!
                              float, 1.f * aod::hf_cand::pxProng0 + 1.f * aod::hf_cand::pxProng1);
DECLARE_SOA_EXPRESSION_COLUMN(Py, py, //!
                              float, 1.f * aod::hf_cand::pyProng0 + 1.f * aod::hf_cand::pyProng1);
DECLARE_SOA_EXPRESSION_COLUMN(Pz, pz, //!
                              float, 1.f * aod::hf_cand::pzProng0 + 1.f * aod::hf_cand::pzProng1);
DECLARE_SOA_DYNAMIC_COLUMN(ImpactParameterProduct, impactParameterProduct, //!
                           [](float dca1, float dca2) -> float { return dca1 * dca2; });
DECLARE_SOA_DYNAMIC_COLUMN(M, m, //!
                           [](float px0, float py0, float pz0, float px1, float py1, float pz1, const std::array<double, 2>& m) -> float { return RecoDecay::m(std::array{std::array{px0, py0, pz0}, std::array{px1, py1, pz1}}, m); });
DECLARE_SOA_DYNAMIC_COLUMN(M2, m2, //!
                           [](float px0, float py0, float pz0, float px1, float py1, float pz1, const std::array<double, 2>& m) -> float { return RecoDecay::m2(std::array{std::array{px0, py0, pz0}, std::array{px1, py1, pz1}}, m); });
DECLARE_SOA_DYNAMIC_COLUMN(CosThetaStar, cosThetaStar, //!
                           [](float px0, float py0, float pz0, float px1, float py1, float pz1, const std::array<double, 2>& m, double mTot, int iProng) -> float { return RecoDecay::cosThetaStar(std::array{std::array{px0, py0, pz0}, std::array{px1, py1, pz1}}, m, mTot, iProng); });
DECLARE_SOA_DYNAMIC_COLUMN(ImpactParameterProngSqSum, impactParameterProngSqSum, //!
                           [](float impParProng0, float impParProng1) -> float { return RecoDecay::sumOfSquares(impParProng0, impParProng1); });
DECLARE_SOA_DYNAMIC_COLUMN(MaxNormalisedDeltaIP, maxNormalisedDeltaIP, //!
                           [](float xVtxP, float yVtxP, float xVtxS, float yVtxS, float errDlxy, float pxM, float pyM, float ip0, float errIp0, float ip1, float errIp1, float px0, float py0, float px1, float py1) -> float { return RecoDecay::maxNormalisedDeltaIP(std::array{xVtxP, yVtxP}, std::array{xVtxS, yVtxS}, errDlxy, std::array{pxM, pyM}, std::array{ip0, ip1}, std::array{errIp0, errIp1}, std::array{std::array{px0, py0}, std::array{px1, py1}}); });
DECLARE_SOA_DYNAMIC_COLUMN(CtXY, ctXY, //!
                           [](float px0, float py0, float pz0, float px1, float py1, float pz1, float xVtxP, float yVtxP, float xVtxS, float yVtxS, const std::array<double, 2>& m) -> float { return RecoDecay::ctXY(std::array{xVtxP, yVtxP}, std::array{xVtxS, yVtxS}, std::array{std::array{px0, py0, pz0}, std::array{px1, py1, pz1}}, m); });
// MC matching result:
DECLARE_SOA_COLUMN(FlagMcMatchRec, flagMcMatchRec, int8_t);         //! reconstruction level
DECLARE_SOA_COLUMN(FlagMcMatchGen, flagMcMatchGen, int8_t);         //! generator level
DECLARE_SOA_COLUMN(OriginMcRec, originMcRec, int8_t);               //! particle origin, reconstruction level
DECLARE_SOA_COLUMN(OriginMcGen, originMcGen, int8_t);               //! particle origin, generator level
DECLARE_SOA_COLUMN(FlagMcDecayChanRec, flagMcDecayChanRec, int8_t); //! resonant decay channel flag, reconstruction level
DECLARE_SOA_COLUMN(FlagMcDecayChanGen, flagMcDecayChanGen, int8_t); //! resonant decay channel flag, reconstruction level

// KF related properties
DECLARE_SOA_COLUMN(KfGeoMassD0, kfGeoMassD0, float);       //! mass of the D0 candidate from the KFParticle geometric fit
DECLARE_SOA_COLUMN(KfGeoMassD0bar, kfGeoMassD0bar, float); //! mass of the D0bar candidate from the KFParticle geometric fit

} // namespace hf_cand_2prong

// general columns
#define HFCAND_COLUMNS                                                                                                                                                                             \
  hf_cand::CollisionId,                                                                                                                                                                            \
    collision::PosX, collision::PosY, collision::PosZ,                                                                                                                                             \
    hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex,                                                                                                               \
    hf_cand::ErrorDecayLength, hf_cand::ErrorDecayLengthXY,                                                                                                                                        \
    hf_cand::Chi2PCA,                                                                                                                                                                              \
    /* dynamic columns */ hf_cand::RSecondaryVertex<hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex>,                                                                                         \
    hf_cand::DecayLength<collision::PosX, collision::PosY, collision::PosZ, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex>,                                      \
    hf_cand::DecayLengthXY<collision::PosX, collision::PosY, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex>,                                                                                \
    hf_cand::DecayLengthNormalised<collision::PosX, collision::PosY, collision::PosZ, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex, hf_cand::ErrorDecayLength>, \
    hf_cand::DecayLengthXYNormalised<collision::PosX, collision::PosY, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ErrorDecayLengthXY>,                                         \
    /* prong 0 */ hf_cand::ImpactParameterNormalised0<hf_cand::ImpactParameter0, hf_cand::ErrorImpactParameter0>,                                                                                  \
    hf_cand::PtProng0<hf_cand::PxProng0, hf_cand::PyProng0>,                                                                                                                                       \
    hf_cand::Pt2Prong0<hf_cand::PxProng0, hf_cand::PyProng0>,                                                                                                                                      \
    hf_cand::PVectorProng0<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0>,                                                                                                               \
    /* prong 1 */ hf_cand::ImpactParameterNormalised1<hf_cand::ImpactParameter1, hf_cand::ErrorImpactParameter1>,                                                                                  \
    hf_cand::PtProng1<hf_cand::PxProng1, hf_cand::PyProng1>,                                                                                                                                       \
    hf_cand::Pt2Prong1<hf_cand::PxProng1, hf_cand::PyProng1>,                                                                                                                                      \
    hf_cand::PVectorProng1<hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>

// 2-prong decay candidate table
DECLARE_SOA_TABLE(HfCand2ProngBase, "AOD", "HFCAND2PBASE", //!
                  o2::soa::Index<>,
                  // general columns
                  HFCAND_COLUMNS,
                  // 2-prong specific columns
                  hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0,
                  hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1,
                  hf_cand::ImpactParameter0, hf_cand::ImpactParameter1,
                  hf_cand::ErrorImpactParameter0, hf_cand::ErrorImpactParameter1,
                  hf_cand::ImpactParameterZ0, hf_cand::ImpactParameterZ1,
                  hf_cand::ErrorImpactParameterZ0, hf_cand::ErrorImpactParameterZ1,
                  hf_track_index::Prong0Id, hf_track_index::Prong1Id, hf_cand::NProngsContributorsPV, hf_cand::BitmapProngsContributorsPV,
                  hf_track_index::HFflag,
                  /* dynamic columns */
                  hf_cand_2prong::M<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>,
                  hf_cand_2prong::M2<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>,
                  hf_cand_2prong::ImpactParameterProduct<hf_cand::ImpactParameter0, hf_cand::ImpactParameter1>,
                  hf_cand_2prong::CosThetaStar<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>,
                  hf_cand_2prong::ImpactParameterProngSqSum<hf_cand::ImpactParameter0, hf_cand::ImpactParameter1>,
                  /* dynamic columns that use candidate momentum components */
                  hf_cand::Pt<hf_cand_2prong::Px, hf_cand_2prong::Py>,
                  hf_cand::Pt2<hf_cand_2prong::Px, hf_cand_2prong::Py>,
                  hf_cand::P<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::P2<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::PVector<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::Cpa<collision::PosX, collision::PosY, collision::PosZ, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex, hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::CpaXY<collision::PosX, collision::PosY, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand_2prong::Px, hf_cand_2prong::Py>,
                  hf_cand::Ct<collision::PosX, collision::PosY, collision::PosZ, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex, hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::ImpactParameterXY<collision::PosX, collision::PosY, collision::PosZ, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex, hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand_2prong::MaxNormalisedDeltaIP<collision::PosX, collision::PosY, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ErrorDecayLengthXY, hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand::ImpactParameter0, hf_cand::ErrorImpactParameter0, hf_cand::ImpactParameter1, hf_cand::ErrorImpactParameter1, hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PxProng1, hf_cand::PyProng1>,
                  hf_cand::Eta<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::Phi<hf_cand_2prong::Px, hf_cand_2prong::Py>,
                  hf_cand::Y<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::E<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::E2<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>);

// extended table with expression columns that can be used as arguments of dynamic columns
DECLARE_SOA_EXTENDED_TABLE_USER(HfCand2ProngExt, HfCand2ProngBase, "HFCAND2PEXT", //!
                                hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz);

DECLARE_SOA_TABLE(HfCand2Prong0PidPi, "AOD", "HFCAND2P0PIDPI", //!
                  hf_cand::NSigTpcPi0, hf_cand::NSigTofPi0,
                  hf_cand::TpcTofNSigmaPi0<hf_cand::NSigTpcPi0, hf_cand::NSigTofPi0>);
DECLARE_SOA_TABLE(HfCand2Prong1PidPi, "AOD", "HFCAND2P1PIDPI", //!
                  hf_cand::NSigTpcPi1, hf_cand::NSigTofPi1,
                  hf_cand::TpcTofNSigmaPi1<hf_cand::NSigTpcPi1, hf_cand::NSigTofPi1>);
DECLARE_SOA_TABLE(HfCand2Prong0PidKa, "AOD", "HFCAND2P0PIDKA", //!
                  hf_cand::NSigTpcKa0, hf_cand::NSigTofKa0,
                  hf_cand::TpcTofNSigmaKa0<hf_cand::NSigTpcKa0, hf_cand::NSigTofKa0>);
DECLARE_SOA_TABLE(HfCand2Prong1PidKa, "AOD", "HFCAND2P1PIDKA", //!
                  hf_cand::NSigTpcKa1, hf_cand::NSigTofKa1,
                  hf_cand::TpcTofNSigmaKa1<hf_cand::NSigTpcKa1, hf_cand::NSigTofKa1>);

using HfCand2Prong = HfCand2ProngExt;
using HfCand2ProngWPid = soa::Join<HfCand2Prong, HfCand2Prong0PidPi, HfCand2Prong0PidKa, HfCand2Prong1PidPi, HfCand2Prong1PidKa>;

DECLARE_SOA_TABLE(HfCand2ProngKF, "AOD", "HFCAND2PKF",
                  hf_cand::KfTopolChi2OverNdf,
                  hf_cand_2prong::KfGeoMassD0, hf_cand_2prong::KfGeoMassD0bar);

// table with results of reconstruction level MC matching
DECLARE_SOA_TABLE(HfCand2ProngMcRec, "AOD", "HFCAND2PMCREC", //!
                  hf_cand_2prong::FlagMcMatchRec,
                  hf_cand_2prong::OriginMcRec,
                  hf_cand_2prong::FlagMcDecayChanRec,
                  hf_cand::PtBhadMotherPart,
                  hf_cand::PdgBhadMotherPart,
                  hf_cand::NTracksDecayed,
                  hf_cand::NInteractionsWithMaterial);

// table with results of generator level MC matching
DECLARE_SOA_TABLE(HfCand2ProngMcGen, "AOD", "HFCAND2PMCGEN", //!
                  hf_cand_2prong::FlagMcMatchGen,
                  hf_cand_2prong::OriginMcGen,
                  hf_cand_2prong::FlagMcDecayChanGen,
                  hf_cand::IdxBhadMotherPart);

// cascade decay candidate table

namespace hf_cand_casc
{
DECLARE_SOA_EXPRESSION_COLUMN(Px, px, //! px of candidate
                              float, 1.f * aod::hf_cand::pxProng0 + 1.f * aod::hf_cand::pxProng1);
DECLARE_SOA_EXPRESSION_COLUMN(Py, py, //! py of candidate
                              float, 1.f * aod::hf_cand::pyProng0 + 1.f * aod::hf_cand::pyProng1);
DECLARE_SOA_EXPRESSION_COLUMN(Pz, pz, //! pz of candidate
                              float, 1.f * aod::hf_cand::pzProng0 + 1.f * aod::hf_cand::pzProng1);
// DECLARE_SOA_DYNAMIC_COLUMN(M, m, [](float px0, float py0, float pz0, float px1, float py1, float pz1, const std::array<double, 2>& m) { return RecoDecay::m(std::array{std::array{px0, py0, pz0}, std::array{px1, py1, pz1}}, m); });
DECLARE_SOA_DYNAMIC_COLUMN(PtV0Pos, ptV0Pos, //! pt of the positive V0 daughter
                           [](float px, float py) { return RecoDecay::pt(px, py); });
DECLARE_SOA_DYNAMIC_COLUMN(PtV0Neg, ptV0Neg, //! pt of the negative V0 daughter
                           [](float px, float py) { return RecoDecay::pt(px, py); });
DECLARE_SOA_DYNAMIC_COLUMN(CtV0, ctV0, //! c*t of the V0
                           [](float xVtxP, float yVtxP, float zVtxP, float xVtxS, float yVtxS, float zVtxS, float px, float py, float pz, double m) -> float { return RecoDecay::ct(std::array{px, py, pz}, RecoDecay::distance(std::array{xVtxP, yVtxP, zVtxP}, std::array{xVtxS, yVtxS, zVtxS}), m); });
DECLARE_SOA_COLUMN(FlagMcMatchRec, flagMcMatchRec, int8_t); //! reconstruction level
DECLARE_SOA_COLUMN(FlagMcMatchGen, flagMcMatchGen, int8_t); //! generator level
DECLARE_SOA_COLUMN(OriginMcRec, originMcRec, int8_t);       //! particle origin, reconstruction level
DECLARE_SOA_COLUMN(OriginMcGen, originMcGen, int8_t);       //! particle origin, generator level
DECLARE_SOA_COLUMN(V0X, v0X, float);                        //! X position of V0 decay
DECLARE_SOA_COLUMN(V0Y, v0Y, float);                        //! Y position of V0 decay
DECLARE_SOA_COLUMN(V0Z, v0Z, float);                        //! Z position of V0 decay
} // namespace hf_cand_casc

DECLARE_SOA_TABLE(HfCandCascBase, "AOD", "HFCANDCASCBASE", //!
                                                           // general columns
                  HFCAND_COLUMNS,
                  // cascade specific columns
                  hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0,
                  hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1,
                  hf_cand::ImpactParameter0, hf_cand::ImpactParameter1,
                  hf_cand::ErrorImpactParameter0, hf_cand::ErrorImpactParameter1,
                  hf_track_index::Prong0Id,
                  hf_track_index::V0Id, // V0 index
                  // V0
                  hf_cand_casc::V0X, hf_cand_casc::V0Y, hf_cand_casc::V0Z,
                  v0data::PosTrackId, v0data::NegTrackId, // indices of V0 tracks in FullTracks table
                  v0data::PxPos, v0data::PyPos, v0data::PzPos, v0data::PxNeg, v0data::PyNeg, v0data::PzNeg,
                  v0data::DCAV0Daughters,
                  v0data::DCAPosToPV, // this is the impact param wrt prim vtx in xy!
                  v0data::DCANegToPV, // this is the impact param wrt prim vtx in xy!
                  v0data::V0CosPA,
                  /* dynamic columns */
                  hf_cand_2prong::M<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>,
                  hf_cand_2prong::M2<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>,
                  hf_cand_2prong::ImpactParameterProduct<hf_cand::ImpactParameter0, hf_cand::ImpactParameter1>,
                  hf_cand_2prong::CosThetaStar<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>,
                  hf_cand_2prong::ImpactParameterProngSqSum<hf_cand::ImpactParameter0, hf_cand::ImpactParameter1>,
                  /* dynamic columns that use candidate momentum components */
                  hf_cand::Pt<hf_cand_casc::Px, hf_cand_casc::Py>,
                  hf_cand::Pt2<hf_cand_casc::Px, hf_cand_casc::Py>,
                  hf_cand::P<hf_cand_casc::Px, hf_cand_casc::Py, hf_cand_casc::Pz>,
                  hf_cand::P2<hf_cand_casc::Px, hf_cand_casc::Py, hf_cand_casc::Pz>,
                  hf_cand::PVector<hf_cand_casc::Px, hf_cand_casc::Py, hf_cand_casc::Pz>,
                  hf_cand::Cpa<collision::PosX, collision::PosY, collision::PosZ, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex, hf_cand_casc::Px, hf_cand_casc::Py, hf_cand_casc::Pz>,
                  hf_cand::CpaXY<collision::PosX, collision::PosY, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand_casc::Px, hf_cand_casc::Py>,
                  hf_cand::Ct<collision::PosX, collision::PosY, collision::PosZ, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex, hf_cand_casc::Px, hf_cand_casc::Py, hf_cand_casc::Pz>,
                  hf_cand::ImpactParameterXY<collision::PosX, collision::PosY, collision::PosZ, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex, hf_cand_casc::Px, hf_cand_casc::Py, hf_cand_casc::Pz>,
                  hf_cand_2prong::MaxNormalisedDeltaIP<collision::PosX, collision::PosY, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ErrorDecayLengthXY, hf_cand_casc::Px, hf_cand_casc::Py, hf_cand::ImpactParameter0, hf_cand::ErrorImpactParameter0, hf_cand::ImpactParameter1, hf_cand::ErrorImpactParameter1, hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PxProng1, hf_cand::PyProng1>,
                  hf_cand::Eta<hf_cand_casc::Px, hf_cand_casc::Py, hf_cand_casc::Pz>,
                  hf_cand::Phi<hf_cand_casc::Px, hf_cand_casc::Py>,
                  hf_cand::Y<hf_cand_casc::Px, hf_cand_casc::Py, hf_cand_casc::Pz>,
                  hf_cand::E<hf_cand_casc::Px, hf_cand_casc::Py, hf_cand_casc::Pz>,
                  hf_cand::E2<hf_cand_casc::Px, hf_cand_casc::Py, hf_cand_casc::Pz>,
                  // dynamic columns from V0
                  hf_cand_casc::PtV0Pos<v0data::PxPos, v0data::PyPos>, // pT of positive V0 daughter
                  hf_cand_casc::PtV0Neg<v0data::PxNeg, v0data::PyNeg>, // pT of negative V0 daughter
                  v0data::V0Radius<hf_cand_casc::V0X, hf_cand_casc::V0Y>,
                  v0data::legacy::MLambda<v0data::PxPos, v0data::PyPos, v0data::PzPos, v0data::PxNeg, v0data::PyNeg, v0data::PzNeg>,
                  v0data::legacy::MAntiLambda<v0data::PxPos, v0data::PyPos, v0data::PzPos, v0data::PxNeg, v0data::PyNeg, v0data::PzNeg>,
                  v0data::legacy::MK0Short<v0data::PxPos, v0data::PyPos, v0data::PzPos, v0data::PxNeg, v0data::PyNeg, v0data::PzNeg>,
                  v0data::MGamma<v0data::PxPos, v0data::PyPos, v0data::PzPos, v0data::PxNeg, v0data::PyNeg, v0data::PzNeg>,
                  hf_cand_casc::CtV0<hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex, hf_cand_casc::V0X, hf_cand_casc::V0Y, hf_cand_casc::V0Z, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>);
//                  ,
//                  v0data::MLambda<v0data::PxPos, v0data::PyPos, v0data::PzPos, v0data::PxNeg, v0data::PyNeg, v0data::PzNeg>,
//                  v0data::MAntiLambda<v0data::PxPos, v0data::PyPos, v0data::PzPos, v0data::PxNeg, v0data::PyNeg, v0data::PzNeg>,
//                  v0data::MK0Short<v0data::PxPos, v0data::PyPos, v0data::PzPos, v0data::PxNeg, v0data::PyNeg, v0data::PzNeg>);

// extended table with expression columns that can be used as arguments of dynamic columns
DECLARE_SOA_EXTENDED_TABLE_USER(HfCandCascExt, HfCandCascBase, "HFCANDCASCEXT", //!
                                hf_cand_casc::Px, hf_cand_casc::Py, hf_cand_casc::Pz);

using HfCandCascade = HfCandCascExt;

// table with results of reconstruction level MC matching for Cascade
DECLARE_SOA_TABLE(HfCandCascadeMcRec, "AOD", "HFCANDCASCMCREC", //!
                  hf_cand_casc::FlagMcMatchRec,
                  hf_cand_casc::OriginMcRec,
                  hf_cand::PtBhadMotherPart,
                  hf_cand::PdgBhadMotherPart);

// table with results of generator level MC matching
DECLARE_SOA_TABLE(HfCandCascadeMcGen, "AOD", "HFCANDCASCMCGEN", //!
                  hf_cand_casc::FlagMcMatchGen,
                  hf_cand_casc::OriginMcGen,
                  hf_cand::IdxBhadMotherPart);

// specific BPlus candidate properties
namespace hf_cand_bplus
{
DECLARE_SOA_INDEX_COLUMN_FULL(Prong0, prong0, int, HfCand2Prong, "_0"); // D0 index
// MC matching result:
DECLARE_SOA_COLUMN(FlagMcMatchRec, flagMcMatchRec, int8_t);                // main decay channel, reconstruction level
DECLARE_SOA_COLUMN(FlagMcMatchGen, flagMcMatchGen, int8_t);                // main decay channel, generator level
DECLARE_SOA_COLUMN(FlagMcDecayChanRec, flagMcDecayChanRec, int8_t);        // resonant decay channel, reconstruction level
DECLARE_SOA_COLUMN(FlagMcDecayChanGen, flagMcDecayChanGen, int8_t);        // resonant decay channel, generator level
DECLARE_SOA_COLUMN(OriginMcRec, originMcRec, int8_t);                      // particle origin, reconstruction level
DECLARE_SOA_COLUMN(OriginMcGen, originMcGen, int8_t);                      // particle origin, generator level
DECLARE_SOA_COLUMN(FlagWrongCollision, flagWrongCollision, int8_t);        // reconstruction level
DECLARE_SOA_COLUMN(DebugMcRec, debugMcRec, int8_t);                        // debug flag for mis-association reconstruction level
DECLARE_SOA_DYNAMIC_COLUMN(ImpactParameterProduct, impactParameterProduct, // Impact parameter product for B+ -> J/Psi K
                           [](float px0, float py0, float pz0, float px1, float py1, float pz1, float xVtxP, float yVtxP, float zVtxP, float xVtxS, float yVtxS, float zVtxS, float impParK) -> float { return impParK * RecoDecay::impParXY(std::array{xVtxP, yVtxP, zVtxP}, std::array{xVtxS, yVtxS, zVtxS}, RecoDecay::pVec(std::array{px0, py0, pz0}, std::array{px1, py1, pz1})); });
DECLARE_SOA_DYNAMIC_COLUMN(ImpactParameterProductJpsi, impactParameterProductJpsi, // J/Psi impact parameter for B+ -> J/Psi K
                           [](float dcaDauPos, float dcaDauNeg) -> float { return dcaDauPos * dcaDauNeg; });

enum DecayTypeMc : uint8_t { BplusToD0PiToKPiPi = 0,
                             BplusToD0KToKPiK,
                             PartlyRecoDecay,
                             OtherDecay,
                             NDecayTypeMc };

enum class DecayTypeBToJpsiMc : uint8_t { BplusToJpsiKToMuMuK = 0,
                                          PartlyRecoDecay,
                                          OtherDecay,
                                          NDecayTypeMc };
} // namespace hf_cand_bplus

// declare dedicated BPlus decay candidate table
DECLARE_SOA_TABLE(HfCandBplusBase, "AOD", "HFCANDBPLUSBASE",
                  // general columns
                  HFCAND_COLUMNS,
                  // 2-prong specific columns
                  o2::soa::Index<>,
                  hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0,
                  hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1,
                  hf_cand::ImpactParameter0, hf_cand::ImpactParameter1,
                  hf_cand::ErrorImpactParameter0, hf_cand::ErrorImpactParameter1,
                  /* dynamic columns */
                  hf_cand_2prong::M<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>,
                  hf_cand_2prong::M2<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>,
                  hf_cand_2prong::ImpactParameterProduct<hf_cand::ImpactParameter0, hf_cand::ImpactParameter1>,
                  hf_cand_2prong::CosThetaStar<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>,
                  hf_cand_2prong::ImpactParameterProngSqSum<hf_cand::ImpactParameter0, hf_cand::ImpactParameter1>,
                  /* dynamic columns that use candidate momentum components */
                  hf_cand::Pt<hf_cand_2prong::Px, hf_cand_2prong::Py>,
                  hf_cand::Pt2<hf_cand_2prong::Px, hf_cand_2prong::Py>,
                  hf_cand::P<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::P2<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::PVector<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::Cpa<collision::PosX, collision::PosY, collision::PosZ, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex, hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::CpaXY<collision::PosX, collision::PosY, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand_2prong::Px, hf_cand_2prong::Py>,
                  hf_cand::Ct<collision::PosX, collision::PosY, collision::PosZ, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex, hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::ImpactParameterXY<collision::PosX, collision::PosY, collision::PosZ, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex, hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand_2prong::MaxNormalisedDeltaIP<collision::PosX, collision::PosY, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ErrorDecayLengthXY, hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand::ImpactParameter0, hf_cand::ErrorImpactParameter0, hf_cand::ImpactParameter1, hf_cand::ErrorImpactParameter1, hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PxProng1, hf_cand::PyProng1>,
                  hf_cand::Eta<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::Phi<hf_cand_2prong::Px, hf_cand_2prong::Py>,
                  hf_cand::Y<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::E<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::E2<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand_2prong::CtXY<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1, collision::PosX, collision::PosY, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex>);

// extended table with expression columns that can be used as arguments of dynamic columns
DECLARE_SOA_EXTENDED_TABLE_USER(HfCandBplusExt, HfCandBplusBase, "HFCANDBPLUSEXT",
                                hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz);

DECLARE_SOA_TABLE(HfCandBplusProngs, "AOD", "HFCANDBPPRONGS",
                  hf_cand_bplus::Prong0Id, hf_track_index::Prong1Id);

using HfCandBplus = soa::Join<HfCandBplusExt, HfCandBplusProngs>;

// table with results of reconstruction level MC matching
DECLARE_SOA_TABLE(HfCandBplusMcRec, "AOD", "HFCANDBPMCREC",
                  hf_cand_bplus::FlagMcMatchRec,
                  hf_cand_bplus::FlagMcDecayChanRec,
                  hf_cand_bplus::OriginMcRec);

// table with results of generator level MC matching
DECLARE_SOA_TABLE(HfCandBplusMcGen, "AOD", "HFCANDBPMCGEN",
                  hf_cand_bplus::FlagMcMatchGen,
                  hf_cand_bplus::FlagMcDecayChanGen,
                  hf_cand_bplus::OriginMcGen);

// specific 3-prong decay properties
namespace hf_cand_3prong
{
DECLARE_SOA_EXPRESSION_COLUMN(Px, px, //!
                              float, 1.f * aod::hf_cand::pxProng0 + 1.f * aod::hf_cand::pxProng1 + 1.f * aod::hf_cand::pxProng2);
DECLARE_SOA_EXPRESSION_COLUMN(Py, py, //!
                              float, 1.f * aod::hf_cand::pyProng0 + 1.f * aod::hf_cand::pyProng1 + 1.f * aod::hf_cand::pyProng2);
DECLARE_SOA_EXPRESSION_COLUMN(Pz, pz, //!
                              float, 1.f * aod::hf_cand::pzProng0 + 1.f * aod::hf_cand::pzProng1 + 1.f * aod::hf_cand::pzProng2);
DECLARE_SOA_DYNAMIC_COLUMN(M, m, //!
                           [](float px0, float py0, float pz0, float px1, float py1, float pz1, float px2, float py2, float pz2, const std::array<double, 3>& m) -> float { return RecoDecay::m(std::array{std::array{px0, py0, pz0}, std::array{px1, py1, pz1}, std::array{px2, py2, pz2}}, m); });
DECLARE_SOA_DYNAMIC_COLUMN(M2, m2, //!
                           [](float px0, float py0, float pz0, float px1, float py1, float pz1, float px2, float py2, float pz2, const std::array<double, 3>& m) -> float { return RecoDecay::m2(std::array{std::array{px0, py0, pz0}, std::array{px1, py1, pz1}, std::array{px2, py2, pz2}}, m); });
DECLARE_SOA_DYNAMIC_COLUMN(ImpactParameterProngSqSum, impactParameterProngSqSum, //!
                           [](float impParProng0, float impParProng1, float impParProng2) -> float { return RecoDecay::sumOfSquares(impParProng0, impParProng1, impParProng2); });
DECLARE_SOA_DYNAMIC_COLUMN(MaxNormalisedDeltaIP, maxNormalisedDeltaIP, //!
                           [](float xVtxP, float yVtxP, float xVtxS, float yVtxS, float errDlxy, float pxM, float pyM, float ip0, float errIp0, float ip1, float errIp1, float ip2, float errIp2, float px0, float py0, float px1, float py1, float px2, float py2) -> float { return RecoDecay::maxNormalisedDeltaIP(std::array{xVtxP, yVtxP}, std::array{xVtxS, yVtxS}, errDlxy, std::array{pxM, pyM}, std::array{ip0, ip1, ip2}, std::array{errIp0, errIp1, errIp2}, std::array{std::array{px0, py0}, std::array{px1, py1}, std::array{px2, py2}}); });
DECLARE_SOA_DYNAMIC_COLUMN(CtXY, ctXY, //!
                           [](float px0, float py0, float pz0, float px1, float py1, float pz1, float px2, float py2, float pz2, float xVtxP, float yVtxP, float xVtxS, float yVtxS, const std::array<double, 3>& m) -> float {
                             return RecoDecay::ctXY(std::array{xVtxP, yVtxP}, std::array{xVtxS, yVtxS}, std::array{std::array{px0, py0, pz0}, std::array{px1, py1, pz1}, std::array{px2, py2, pz2}}, m);
                           });
// MC matching result:
DECLARE_SOA_COLUMN(FlagMcMatchRec, flagMcMatchRec, int8_t);         //! reconstruction level
DECLARE_SOA_COLUMN(FlagMcMatchGen, flagMcMatchGen, int8_t);         //! generator level
DECLARE_SOA_COLUMN(OriginMcRec, originMcRec, int8_t);               //! particle origin, reconstruction level
DECLARE_SOA_COLUMN(OriginMcGen, originMcGen, int8_t);               //! particle origin, generator level
DECLARE_SOA_COLUMN(IsCandidateSwapped, isCandidateSwapped, int8_t); //! swapping of the prongs order
DECLARE_SOA_COLUMN(FlagMcDecayChanRec, flagMcDecayChanRec, int8_t); //! resonant decay channel flag, reconstruction level
DECLARE_SOA_COLUMN(FlagMcDecayChanGen, flagMcDecayChanGen, int8_t); //! resonant decay channel flag, generator level

// Ds± → K± K∓ π± or D± → K± K∓ π±

enum DecayChannelDToKKPi {
  DsToPhiPi = 1,
  DsToK0starK,
  DplusToPhiPi,  // used to describe D+ in MC production for Ds analysis
  DplusToK0starK // used to describe D+ in MC production for Ds analysis
};

// KF related properties
DECLARE_SOA_COLUMN(KfXPVError, kfXPVError, float);                       //! error of X coordinate of the event's primary vertex
DECLARE_SOA_COLUMN(KfYPVError, kfYPVError, float);                       //! error of Y coordinate of the event's primary vertex
DECLARE_SOA_COLUMN(KfZPVError, kfZPVError, float);                       //! error of Z coordinate of the event's primary vertex
DECLARE_SOA_COLUMN(KfXError, kfXError, float);                           //! error of candidate's decay point X coordinate from the KFParticle fit
DECLARE_SOA_COLUMN(KfYError, kfYError, float);                           //! error of candidate's decay point Y coordinate from the KFParticle fit
DECLARE_SOA_COLUMN(KfZError, kfZError, float);                           //! error of candidate's decay point Z coordinate from the KFParticle fit
DECLARE_SOA_COLUMN(KfMassPKPi, kfMassPKPi, float);                       //! mass of the PKPi candidate from the KFParticle fit
DECLARE_SOA_COLUMN(KfMassPiKP, kfMassPiKP, float);                       //! mass of the PiKP candidate from the KFParticle fit
DECLARE_SOA_COLUMN(KfMassPiKPi, kfMassPiKPi, float);                     //! mass of the PiKPi candidate from the KFParticle fit
DECLARE_SOA_COLUMN(KfMassKKPi, kfMassKKPi, float);                       //! mass of the KKPi candidate from the KFParticle fit
DECLARE_SOA_COLUMN(KfMassPiKK, kfMassPiKK, float);                       //! mass of the PiKK candidate from the KFParticle fit
DECLARE_SOA_COLUMN(KfMassKPi, kfMassKPi, float);                         //! mass of the KPi pair from the KFParticle fit
DECLARE_SOA_COLUMN(KfMassPiK, kfMassPiK, float);                         //! mass of the PiK pair from the KFParticle fit
DECLARE_SOA_COLUMN(KfPx, kfPx, float);                                   //! Px of the candidate from the KFParticle fit
DECLARE_SOA_COLUMN(KfPy, kfPy, float);                                   //! Py of the candidate from the KFParticle fit
DECLARE_SOA_COLUMN(KfPz, kfPz, float);                                   //! Pz of the candidate from the KFParticle fit
DECLARE_SOA_COLUMN(KfErrorPx, kfErrorPx, float);                         //! Px error of the candidate from the KFParticle fit
DECLARE_SOA_COLUMN(KfErrorPy, kfErrorPy, float);                         //! Py error of the candidate from the KFParticle fit
DECLARE_SOA_COLUMN(KfErrorPz, kfErrorPz, float);                         //! Pz error of the candidate from the KFParticle fit
DECLARE_SOA_COLUMN(KfChi2PrimProng0, kfChi2PrimProng0, float);           //! chi2 primary of the first prong
DECLARE_SOA_COLUMN(KfChi2PrimProng1, kfChi2PrimProng1, float);           //! chi2 primary of the second prong
DECLARE_SOA_COLUMN(KfChi2PrimProng2, kfChi2PrimProng2, float);           //! chi2 primary of the third prong
DECLARE_SOA_COLUMN(KfDcaProng0Prong1, kfDcaProng0Prong1, float);         //! DCA between first and second prongs
DECLARE_SOA_COLUMN(KfDcaProng0Prong2, kfDcaProng0Prong2, float);         //! DCA between first and third prongs
DECLARE_SOA_COLUMN(KfDcaProng1Prong2, kfDcaProng1Prong2, float);         //! DCA between second and third prongs
DECLARE_SOA_COLUMN(KfChi2GeoProng0Prong1, kfChi2GeoProng0Prong1, float); //! chi2 geo between first and second prongs
DECLARE_SOA_COLUMN(KfChi2GeoProng0Prong2, kfChi2GeoProng0Prong2, float); //! chi2 geo between first and third prongs
DECLARE_SOA_COLUMN(KfChi2GeoProng1Prong2, kfChi2GeoProng1Prong2, float); //! chi2 geo between second and third prongs
DECLARE_SOA_COLUMN(KfChi2Geo, kfChi2Geo, float);                         //! chi2 geo of the full candidate
DECLARE_SOA_COLUMN(KfChi2Topo, kfChi2Topo, float);                       //! chi2 topo of the full candidate (chi2prim of candidate to PV)
DECLARE_SOA_COLUMN(KfDecayLength, kfDecayLength, float);                 //! decay length
DECLARE_SOA_COLUMN(KfDecayLengthError, kfDecayLengthError, float);       //! decay length error

} // namespace hf_cand_3prong

// 3-prong decay candidate table
DECLARE_SOA_TABLE(HfCand3ProngBase, "AOD", "HFCAND3PBASE", //!
                  o2::soa::Index<>,
                  // general columns
                  HFCAND_COLUMNS,
                  // 3-prong specific columns
                  hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0,
                  hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1,
                  hf_cand::PxProng2, hf_cand::PyProng2, hf_cand::PzProng2,
                  hf_cand::ImpactParameter0, hf_cand::ImpactParameter1, hf_cand::ImpactParameter2,
                  hf_cand::ErrorImpactParameter0, hf_cand::ErrorImpactParameter1, hf_cand::ErrorImpactParameter2,
                  hf_cand::ImpactParameterZ0, hf_cand::ImpactParameterZ1, hf_cand::ImpactParameterZ2,
                  hf_cand::ErrorImpactParameterZ0, hf_cand::ErrorImpactParameterZ1, hf_cand::ErrorImpactParameterZ2,
                  hf_track_index::Prong0Id, hf_track_index::Prong1Id, hf_track_index::Prong2Id, hf_cand::NProngsContributorsPV, hf_cand::BitmapProngsContributorsPV,
                  hf_track_index::HFflag,
                  /* dynamic columns */
                  hf_cand_3prong::M<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1, hf_cand::PxProng2, hf_cand::PyProng2, hf_cand::PzProng2>,
                  hf_cand_3prong::M2<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1, hf_cand::PxProng2, hf_cand::PyProng2, hf_cand::PzProng2>,
                  hf_cand_3prong::ImpactParameterProngSqSum<hf_cand::ImpactParameter0, hf_cand::ImpactParameter1, hf_cand::ImpactParameter2>,
                  /* prong 2 */
                  hf_cand::ImpactParameterNormalised2<hf_cand::ImpactParameter2, hf_cand::ErrorImpactParameter2>,
                  hf_cand::PtProng2<hf_cand::PxProng2, hf_cand::PyProng2>,
                  hf_cand::Pt2Prong2<hf_cand::PxProng2, hf_cand::PyProng2>,
                  hf_cand::PVectorProng2<hf_cand::PxProng2, hf_cand::PyProng2, hf_cand::PzProng2>,
                  /* dynamic columns that use candidate momentum components */
                  hf_cand::Pt<hf_cand_3prong::Px, hf_cand_3prong::Py>,
                  hf_cand::Pt2<hf_cand_3prong::Px, hf_cand_3prong::Py>,
                  hf_cand::P<hf_cand_3prong::Px, hf_cand_3prong::Py, hf_cand_3prong::Pz>,
                  hf_cand::P2<hf_cand_3prong::Px, hf_cand_3prong::Py, hf_cand_3prong::Pz>,
                  hf_cand::PVector<hf_cand_3prong::Px, hf_cand_3prong::Py, hf_cand_3prong::Pz>,
                  hf_cand::Cpa<collision::PosX, collision::PosY, collision::PosZ, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex, hf_cand_3prong::Px, hf_cand_3prong::Py, hf_cand_3prong::Pz>,
                  hf_cand::CpaXY<collision::PosX, collision::PosY, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand_3prong::Px, hf_cand_3prong::Py>,
                  hf_cand::Ct<collision::PosX, collision::PosY, collision::PosZ, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex, hf_cand_3prong::Px, hf_cand_3prong::Py, hf_cand_3prong::Pz>,
                  hf_cand::ImpactParameterXY<collision::PosX, collision::PosY, collision::PosZ, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex, hf_cand_3prong::Px, hf_cand_3prong::Py, hf_cand_3prong::Pz>,
                  hf_cand_3prong::MaxNormalisedDeltaIP<collision::PosX, collision::PosY, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ErrorDecayLengthXY, hf_cand_3prong::Px, hf_cand_3prong::Py, hf_cand::ImpactParameter0, hf_cand::ErrorImpactParameter0, hf_cand::ImpactParameter1, hf_cand::ErrorImpactParameter1, hf_cand::ImpactParameter2, hf_cand::ErrorImpactParameter2, hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PxProng2, hf_cand::PyProng2>,
                  hf_cand::Eta<hf_cand_3prong::Px, hf_cand_3prong::Py, hf_cand_3prong::Pz>,
                  hf_cand::Phi<hf_cand_3prong::Px, hf_cand_3prong::Py>,
                  hf_cand::Y<hf_cand_3prong::Px, hf_cand_3prong::Py, hf_cand_3prong::Pz>,
                  hf_cand::E<hf_cand_3prong::Px, hf_cand_3prong::Py, hf_cand_3prong::Pz>,
                  hf_cand::E2<hf_cand_3prong::Px, hf_cand_3prong::Py, hf_cand_3prong::Pz>);

// extended table with expression columns that can be used as arguments of dynamic columns
DECLARE_SOA_EXTENDED_TABLE_USER(HfCand3ProngExt, HfCand3ProngBase, "HFCAND3PEXT", //!
                                hf_cand_3prong::Px, hf_cand_3prong::Py, hf_cand_3prong::Pz);

DECLARE_SOA_TABLE(HfCand3Prong0PidPi, "AOD", "HFCAND3P0PIDPI", //!
                  hf_cand::NSigTpcPi0, hf_cand::NSigTofPi0,
                  hf_cand::TpcTofNSigmaPi0<hf_cand::NSigTpcPi0, hf_cand::NSigTofPi0>);
DECLARE_SOA_TABLE(HfCand3Prong1PidPi, "AOD", "HFCAND3P1PIDPI", //!
                  hf_cand::NSigTpcPi1, hf_cand::NSigTofPi1,
                  hf_cand::TpcTofNSigmaPi1<hf_cand::NSigTpcPi1, hf_cand::NSigTofPi1>);
DECLARE_SOA_TABLE(HfCand3Prong2PidPi, "AOD", "HFCAND3P2PIDPI", //!
                  hf_cand::NSigTpcPi2, hf_cand::NSigTofPi2,
                  hf_cand::TpcTofNSigmaPi2<hf_cand::NSigTpcPi2, hf_cand::NSigTofPi2>);
DECLARE_SOA_TABLE(HfCand3Prong0PidKa, "AOD", "HFCAND3P0PIDKA", //!
                  hf_cand::NSigTpcKa0, hf_cand::NSigTofKa0,
                  hf_cand::TpcTofNSigmaKa0<hf_cand::NSigTpcKa0, hf_cand::NSigTofKa0>);
DECLARE_SOA_TABLE(HfCand3Prong1PidKa, "AOD", "HFCAND3P1PIDKA", //!
                  hf_cand::NSigTpcKa1, hf_cand::NSigTofKa1,
                  hf_cand::TpcTofNSigmaKa1<hf_cand::NSigTpcKa1, hf_cand::NSigTofKa1>);
DECLARE_SOA_TABLE(HfCand3Prong2PidKa, "AOD", "HFCAND3P2PIDKA", //!
                  hf_cand::NSigTpcKa2, hf_cand::NSigTofKa2,
                  hf_cand::TpcTofNSigmaKa2<hf_cand::NSigTpcKa2, hf_cand::NSigTofKa2>);
DECLARE_SOA_TABLE(HfCand3Prong0PidPr, "AOD", "HFCAND3P0PIDPR", //!
                  hf_cand::NSigTpcPr0, hf_cand::NSigTofPr0,
                  hf_cand::TpcTofNSigmaPr0<hf_cand::NSigTpcPr0, hf_cand::NSigTofPr0>);
DECLARE_SOA_TABLE(HfCand3Prong1PidPr, "AOD", "HFCAND3P1PIDPR", //!
                  hf_cand::NSigTpcPr1, hf_cand::NSigTofPr1,
                  hf_cand::TpcTofNSigmaPr1<hf_cand::NSigTpcPr1, hf_cand::NSigTofPr1>);
DECLARE_SOA_TABLE(HfCand3Prong2PidPr, "AOD", "HFCAND3P2PIDPR", //!
                  hf_cand::NSigTpcPr2, hf_cand::NSigTofPr2,
                  hf_cand::TpcTofNSigmaPr2<hf_cand::NSigTpcPr2, hf_cand::NSigTofPr2>);
DECLARE_SOA_TABLE(HfCand3Prong0PidDe, "AOD", "HFCAND3P0PIDDE", //!
                  hf_cand::NSigTpcDe0, hf_cand::NSigTofDe0,
                  hf_cand::TpcTofNSigmaDe0<hf_cand::NSigTpcDe0, hf_cand::NSigTofDe0>);
DECLARE_SOA_TABLE(HfCand3Prong1PidDe, "AOD", "HFCAND3P1PIDDE", //!
                  hf_cand::NSigTpcDe1, hf_cand::NSigTofDe1,
                  hf_cand::TpcTofNSigmaDe1<hf_cand::NSigTpcDe1, hf_cand::NSigTofDe1>);
DECLARE_SOA_TABLE(HfCand3Prong2PidDe, "AOD", "HFCAND3P2PIDDE", //!
                  hf_cand::NSigTpcDe2, hf_cand::NSigTofDe2,
                  hf_cand::TpcTofNSigmaDe2<hf_cand::NSigTpcDe2, hf_cand::NSigTofDe2>);
DECLARE_SOA_TABLE(HfCand3Prong0PidTr, "AOD", "HFCAND3P0PIDTR", //!
                  hf_cand::NSigTpcTr0, hf_cand::NSigTofTr0,
                  hf_cand::TpcTofNSigmaTr0<hf_cand::NSigTpcTr0, hf_cand::NSigTofTr0>);
DECLARE_SOA_TABLE(HfCand3Prong1PidTr, "AOD", "HFCAND3P1PIDTR", //!
                  hf_cand::NSigTpcTr1, hf_cand::NSigTofTr1,
                  hf_cand::TpcTofNSigmaTr1<hf_cand::NSigTpcTr1, hf_cand::NSigTofTr1>);
DECLARE_SOA_TABLE(HfCand3Prong2PidTr, "AOD", "HFCAND3P2PIDTR", //!
                  hf_cand::NSigTpcTr2, hf_cand::NSigTofTr2,
                  hf_cand::TpcTofNSigmaTr2<hf_cand::NSigTpcTr2, hf_cand::NSigTofTr2>);
DECLARE_SOA_TABLE(HfCand3Prong0PidHe, "AOD", "HFCAND3P0PIDHe", //!
                  hf_cand::NSigTpcHe0, hf_cand::NSigTofHe0,
                  hf_cand::TpcTofNSigmaTr0<hf_cand::NSigTpcHe0, hf_cand::NSigTofHe0>);
DECLARE_SOA_TABLE(HfCand3Prong1PidHe, "AOD", "HFCAND3P1PIDHe", //!
                  hf_cand::NSigTpcHe1, hf_cand::NSigTofHe1,
                  hf_cand::TpcTofNSigmaHe1<hf_cand::NSigTpcHe1, hf_cand::NSigTofHe1>);
DECLARE_SOA_TABLE(HfCand3Prong2PidHe, "AOD", "HFCAND3P2PIDHe", //!
                  hf_cand::NSigTpcHe2, hf_cand::NSigTofHe2,
                  hf_cand::TpcTofNSigmaHe2<hf_cand::NSigTpcHe2, hf_cand::NSigTofHe2>);

using HfCand3Prong = HfCand3ProngExt;
using HfCand3ProngWPidPiKaPr = soa::Join<HfCand3Prong, HfCand3Prong0PidPi, HfCand3Prong0PidPr, HfCand3Prong0PidKa, HfCand3Prong1PidPi, HfCand3Prong1PidPr, HfCand3Prong1PidKa, HfCand3Prong2PidPi, HfCand3Prong2PidPr, HfCand3Prong2PidKa>;
using HfCand3ProngWPidPiKa = soa::Join<HfCand3Prong, HfCand3Prong0PidPi, HfCand3Prong0PidKa, HfCand3Prong1PidPi, HfCand3Prong1PidKa, HfCand3Prong2PidPi, HfCand3Prong2PidKa>;
using HfCand3ProngWPidPiKaDe = soa::Join<HfCand3Prong, HfCand3Prong0PidPi, HfCand3Prong0PidDe, HfCand3Prong0PidKa, HfCand3Prong1PidPi, HfCand3Prong1PidDe, HfCand3Prong1PidKa, HfCand3Prong2PidPi, HfCand3Prong2PidDe, HfCand3Prong2PidKa>;
using HfCand3ProngWPidPiKaTr = soa::Join<HfCand3Prong, HfCand3Prong0PidPi, HfCand3Prong0PidTr, HfCand3Prong0PidKa, HfCand3Prong1PidPi, HfCand3Prong1PidTr, HfCand3Prong1PidKa, HfCand3Prong2PidPi, HfCand3Prong2PidTr, HfCand3Prong2PidKa>;
using HfCand3ProngWPidPiKaHe = soa::Join<HfCand3Prong, HfCand3Prong0PidPi, HfCand3Prong0PidHe, HfCand3Prong0PidKa, HfCand3Prong1PidPi, HfCand3Prong1PidHe, HfCand3Prong1PidKa, HfCand3Prong2PidPi, HfCand3Prong2PidHe, HfCand3Prong2PidKa>;

DECLARE_SOA_TABLE(HfCand3ProngKF, "AOD", "HFCAND3PKF",
                  hf_cand_3prong::KfXError, hf_cand_3prong::KfYError, hf_cand_3prong::KfZError,
                  hf_cand_3prong::KfXPVError, hf_cand_3prong::KfYPVError, hf_cand_3prong::KfZPVError,
                  hf_cand_3prong::KfMassPKPi, hf_cand_3prong::KfMassPiKP, hf_cand_3prong::KfMassPiKPi, hf_cand_3prong::KfMassKKPi, hf_cand_3prong::KfMassPiKK, hf_cand_3prong::KfMassKPi, hf_cand_3prong::KfMassPiK,
                  hf_cand_3prong::KfPx, hf_cand_3prong::KfPy, hf_cand_3prong::KfPz,
                  hf_cand_3prong::KfErrorPx, hf_cand_3prong::KfErrorPy, hf_cand_3prong::KfErrorPz,
                  hf_cand_3prong::KfChi2PrimProng0, hf_cand_3prong::KfChi2PrimProng1, hf_cand_3prong::KfChi2PrimProng2,
                  hf_cand_3prong::KfDcaProng1Prong2, hf_cand_3prong::KfDcaProng0Prong2, hf_cand_3prong::KfDcaProng0Prong1,
                  hf_cand_3prong::KfChi2GeoProng1Prong2, hf_cand_3prong::KfChi2GeoProng0Prong2, hf_cand_3prong::KfChi2GeoProng0Prong1,
                  hf_cand_3prong::KfChi2Geo, hf_cand_3prong::KfDecayLength, hf_cand_3prong::KfDecayLengthError, hf_cand_3prong::KfChi2Topo);

// table with results of reconstruction level MC matching
DECLARE_SOA_TABLE(HfCand3ProngMcRec, "AOD", "HFCAND3PMCREC", //!
                  hf_cand_3prong::FlagMcMatchRec,
                  hf_cand_3prong::OriginMcRec,
                  hf_cand_3prong::IsCandidateSwapped,
                  hf_cand_3prong::FlagMcDecayChanRec,
                  hf_cand::PtBhadMotherPart,
                  hf_cand::PdgBhadMotherPart,
                  hf_cand::NTracksDecayed,
                  hf_cand::NInteractionsWithMaterial);

// table with results of generator level MC matching
DECLARE_SOA_TABLE(HfCand3ProngMcGen, "AOD", "HFCAND3PMCGEN", //!
                  hf_cand_3prong::FlagMcMatchGen,
                  hf_cand_3prong::OriginMcGen,
                  hf_cand_3prong::FlagMcDecayChanGen,
                  hf_cand::IdxBhadMotherPart);

// declare dedicated BPlus -> J/Psi K decay candidate table
// convention: prongs 0 and 1 should be J/Psi decay products
DECLARE_SOA_TABLE(HfCandBpJPBase, "AOD", "HFCANDBPJPBASE",
                  // general columns
                  HFCAND_COLUMNS,
                  /* prong 2 */ hf_cand::ImpactParameterNormalised2<hf_cand::ImpactParameter2, hf_cand::ErrorImpactParameter2>,
                  hf_cand::PtProng2<hf_cand::PxProng2, hf_cand::PyProng2>,
                  hf_cand::Pt2Prong2<hf_cand::PxProng2, hf_cand::PyProng2>,
                  hf_cand::PVectorProng2<hf_cand::PxProng2, hf_cand::PyProng2, hf_cand::PzProng2>,
                  // 3-prong specific columns
                  o2::soa::Index<>,
                  hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0,
                  hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1,
                  hf_cand::PxProng2, hf_cand::PyProng2, hf_cand::PzProng2,
                  hf_cand::ImpactParameter0, hf_cand::ImpactParameter1, hf_cand::ImpactParameter2,
                  hf_cand::ErrorImpactParameter0, hf_cand::ErrorImpactParameter1, hf_cand::ErrorImpactParameter2,
                  /* dynamic columns */
                  hf_cand_3prong::M<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1, hf_cand::PxProng2, hf_cand::PyProng2, hf_cand::PzProng2>,
                  hf_cand_3prong::M2<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1, hf_cand::PxProng2, hf_cand::PyProng2, hf_cand::PzProng2>,
                  hf_cand_3prong::ImpactParameterProngSqSum<hf_cand::ImpactParameter0, hf_cand::ImpactParameter1>,
                  hf_cand_bplus::ImpactParameterProduct<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1, collision::PosX, collision::PosY, collision::PosZ, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex, hf_cand::ImpactParameter2>,
                  hf_cand_bplus::ImpactParameterProductJpsi<hf_cand::ImpactParameter0, hf_cand::ImpactParameter1>,
                  /* dynamic columns that use candidate momentum components */
                  hf_cand::Pt<hf_cand_3prong::Px, hf_cand_3prong::Py>,
                  hf_cand::Pt2<hf_cand_3prong::Px, hf_cand_3prong::Py>,
                  hf_cand::P<hf_cand_3prong::Px, hf_cand_3prong::Py, hf_cand_3prong::Pz>,
                  hf_cand::P2<hf_cand_3prong::Px, hf_cand_3prong::Py, hf_cand_3prong::Pz>,
                  hf_cand::PVector<hf_cand_3prong::Px, hf_cand_3prong::Py, hf_cand_3prong::Pz>,
                  hf_cand::Cpa<collision::PosX, collision::PosY, collision::PosZ, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex, hf_cand_3prong::Px, hf_cand_3prong::Py, hf_cand_3prong::Pz>,
                  hf_cand::CpaXY<collision::PosX, collision::PosY, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand_3prong::Px, hf_cand_3prong::Py>,
                  hf_cand::Ct<collision::PosX, collision::PosY, collision::PosZ, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex, hf_cand_3prong::Px, hf_cand_3prong::Py, hf_cand_3prong::Pz>,
                  hf_cand::ImpactParameterXY<collision::PosX, collision::PosY, collision::PosZ, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex, hf_cand_3prong::Px, hf_cand_3prong::Py, hf_cand_3prong::Pz>,
                  hf_cand_3prong::MaxNormalisedDeltaIP<collision::PosX, collision::PosY, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ErrorDecayLengthXY, hf_cand_3prong::Px, hf_cand_3prong::Py, hf_cand::ImpactParameter0, hf_cand::ErrorImpactParameter0, hf_cand::ImpactParameter1, hf_cand::ErrorImpactParameter1, hf_cand::ImpactParameter2, hf_cand::ErrorImpactParameter2, hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PxProng2, hf_cand::PyProng2>,
                  hf_cand::Eta<hf_cand_3prong::Px, hf_cand_3prong::Py, hf_cand_3prong::Pz>,
                  hf_cand::Phi<hf_cand_3prong::Px, hf_cand_3prong::Py>,
                  hf_cand::Y<hf_cand_3prong::Px, hf_cand_3prong::Py, hf_cand_3prong::Pz>,
                  hf_cand::E<hf_cand_3prong::Px, hf_cand_3prong::Py, hf_cand_3prong::Pz>,
                  hf_cand::E2<hf_cand_3prong::Px, hf_cand_3prong::Py, hf_cand_3prong::Pz>,
                  hf_cand_3prong::CtXY<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1, hf_cand::PxProng2, hf_cand::PyProng2, hf_cand::PzProng2, collision::PosX, collision::PosY, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex>);

// extended table with expression columns that can be used as arguments of dynamic columns
DECLARE_SOA_EXTENDED_TABLE_USER(HfCandBpJPExt, HfCandBpJPBase, "HFCANDBPJPEXT",
                                hf_cand_3prong::Px, hf_cand_3prong::Py, hf_cand_3prong::Pz);

DECLARE_SOA_TABLE(HfCandBpJPDaus, "AOD", "HFCANDBPJPDAUS",
                  hf_track_index::Prong0Id, hf_track_index::Prong1Id, hf_track_index::Prong2Id);

using HfCandBplusToJpsi = soa::Join<HfCandBpJPExt, HfCandBpJPDaus>;

namespace hf_cand_casc_lf
{
// mapping of construct method
enum ConstructMethod { DcaFitter = 0,
                       KfParticle };
} // namespace hf_cand_casc_lf

namespace hf_cand_x
{
DECLARE_SOA_INDEX_COLUMN_FULL(Prong0, prong0, int, HfCand2Prong, "_0"); // Jpsi index
// MC matching result:
DECLARE_SOA_COLUMN(FlagMcMatchRec, flagMcMatchRec, int8_t);         // reconstruction level
DECLARE_SOA_COLUMN(FlagMcMatchGen, flagMcMatchGen, int8_t);         // generator level
DECLARE_SOA_COLUMN(OriginMcRec, originMcRec, int8_t);               // particle origin, reconstruction level
DECLARE_SOA_COLUMN(OriginMcGen, originMcGen, int8_t);               // particle origin, generator level
DECLARE_SOA_COLUMN(FlagMcDecayChanRec, flagMcDecayChanRec, int8_t); // resonant decay channel flag, reconstruction level
DECLARE_SOA_COLUMN(FlagMcDecayChanGen, flagMcDecayChanGen, int8_t); // resonant decay channel flag, generator level

} // namespace hf_cand_x

// declare dedicated X candidate table
DECLARE_SOA_TABLE(HfCandXBase, "AOD", "HFCANDXBASE",
                  // general columns
                  HFCAND_COLUMNS,
                  // 3-prong specific columns
                  hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0,
                  hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1,
                  hf_cand::PxProng2, hf_cand::PyProng2, hf_cand::PzProng2,
                  hf_cand::ImpactParameter0, hf_cand::ImpactParameter1, hf_cand::ImpactParameter2,
                  hf_cand::ErrorImpactParameter0, hf_cand::ErrorImpactParameter1, hf_cand::ErrorImpactParameter2,
                  hf_cand_x::Prong0Id, hf_track_index::Prong1Id, hf_track_index::Prong2Id, // note the difference between Jpsi and pion indices
                  hf_track_index::HFflag,
                  /* dynamic columns */
                  hf_cand_3prong::M<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1, hf_cand::PxProng2, hf_cand::PyProng2, hf_cand::PzProng2>,
                  hf_cand_3prong::M2<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1, hf_cand::PxProng2, hf_cand::PyProng2, hf_cand::PzProng2>,
                  /* prong 2 */
                  hf_cand::PtProng2<hf_cand::PxProng2, hf_cand::PyProng2>,
                  hf_cand::Pt2Prong2<hf_cand::PxProng2, hf_cand::PyProng2>,
                  hf_cand::PVectorProng2<hf_cand::PxProng2, hf_cand::PyProng2, hf_cand::PzProng2>,
                  /* dynamic columns that use candidate momentum components */
                  hf_cand::Pt<hf_cand_3prong::Px, hf_cand_3prong::Py>,
                  hf_cand::Pt2<hf_cand_3prong::Px, hf_cand_3prong::Py>,
                  hf_cand::P<hf_cand_3prong::Px, hf_cand_3prong::Py, hf_cand_3prong::Pz>,
                  hf_cand::P2<hf_cand_3prong::Px, hf_cand_3prong::Py, hf_cand_3prong::Pz>,
                  hf_cand::PVector<hf_cand_3prong::Px, hf_cand_3prong::Py, hf_cand_3prong::Pz>,
                  hf_cand::Cpa<collision::PosX, collision::PosY, collision::PosZ, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex, hf_cand_3prong::Px, hf_cand_3prong::Py, hf_cand_3prong::Pz>,
                  hf_cand::CpaXY<collision::PosX, collision::PosY, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand_3prong::Px, hf_cand_3prong::Py>,
                  hf_cand::Ct<collision::PosX, collision::PosY, collision::PosZ, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex, hf_cand_3prong::Px, hf_cand_3prong::Py, hf_cand_3prong::Pz>,
                  hf_cand::ImpactParameterXY<collision::PosX, collision::PosY, collision::PosZ, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex, hf_cand_3prong::Px, hf_cand_3prong::Py, hf_cand_3prong::Pz>,
                  hf_cand_3prong::MaxNormalisedDeltaIP<collision::PosX, collision::PosY, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ErrorDecayLengthXY, hf_cand_3prong::Px, hf_cand_3prong::Py, hf_cand::ImpactParameter0, hf_cand::ErrorImpactParameter0, hf_cand::ImpactParameter1, hf_cand::ErrorImpactParameter1, hf_cand::ImpactParameter2, hf_cand::ErrorImpactParameter2, hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PxProng2, hf_cand::PyProng2>,
                  hf_cand::Eta<hf_cand_3prong::Px, hf_cand_3prong::Py, hf_cand_3prong::Pz>,
                  hf_cand::Phi<hf_cand_3prong::Px, hf_cand_3prong::Py>,
                  hf_cand::Y<hf_cand_3prong::Px, hf_cand_3prong::Py, hf_cand_3prong::Pz>,
                  hf_cand::E<hf_cand_3prong::Px, hf_cand_3prong::Py, hf_cand_3prong::Pz>,
                  hf_cand::E2<hf_cand_3prong::Px, hf_cand_3prong::Py, hf_cand_3prong::Pz>);

// extended table with expression columns that can be used as arguments of dynamic columns
DECLARE_SOA_EXTENDED_TABLE_USER(HfCandXExt, HfCandXBase, "HFCANDXEXT",
                                hf_cand_3prong::Px, hf_cand_3prong::Py, hf_cand_3prong::Pz);

using HfCandX = HfCandXExt;

// table with results of reconstruction level MC matching
DECLARE_SOA_TABLE(HfCandXMcRec, "AOD", "HFCANDXMCREC", //!
                  hf_cand_x::FlagMcMatchRec,
                  hf_cand_x::OriginMcRec,
                  hf_cand_x::FlagMcDecayChanRec);

// table with results of generator level MC matching
DECLARE_SOA_TABLE(HfCandXMcGen, "AOD", "HFCANDXMCGEN", //!
                  hf_cand_x::FlagMcMatchGen,
                  hf_cand_x::OriginMcGen,
                  hf_cand_x::FlagMcDecayChanGen);

// specific Xicc candidate properties
namespace hf_cand_xicc
{
DECLARE_SOA_INDEX_COLUMN_FULL(Prong0, prong0, int, HfCand3Prong, "_0"); // Xic index
// MC matching result:
DECLARE_SOA_COLUMN(FlagMcMatchRec, flagMcMatchRec, int8_t); // reconstruction level
DECLARE_SOA_COLUMN(FlagMcMatchGen, flagMcMatchGen, int8_t); // generator level
DECLARE_SOA_COLUMN(OriginMcRec, originMcRec, int8_t);       // particle origin, reconstruction level
DECLARE_SOA_COLUMN(OriginMcGen, originMcGen, int8_t);       // particle origin, generator level
DECLARE_SOA_COLUMN(DebugMcRec, debugMcRec, int8_t);         // debug flag for mis-association reconstruction level
} // namespace hf_cand_xicc

// declare dedicated Xicc candidate table
DECLARE_SOA_TABLE(HfCandXiccBase, "AOD", "HFCANDXICCBASE",
                  // general columns
                  HFCAND_COLUMNS,
                  // 2-prong specific columns
                  hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0,
                  hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1,
                  hf_cand::ImpactParameter0, hf_cand::ImpactParameter1,
                  hf_cand::ErrorImpactParameter0, hf_cand::ErrorImpactParameter1,
                  hf_cand_xicc::Prong0Id, hf_track_index::Prong1Id,
                  hf_track_index::HFflag,
                  /* dynamic columns */
                  hf_cand_2prong::M<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>,
                  hf_cand_2prong::M2<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>,
                  hf_cand_2prong::ImpactParameterProduct<hf_cand::ImpactParameter0, hf_cand::ImpactParameter1>,
                  /* dynamic columns that use candidate momentum components */
                  hf_cand::Pt<hf_cand_2prong::Px, hf_cand_2prong::Py>,
                  hf_cand::Pt2<hf_cand_2prong::Px, hf_cand_2prong::Py>,
                  hf_cand::P<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::P2<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::PVector<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::Cpa<collision::PosX, collision::PosY, collision::PosZ, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex, hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::CpaXY<collision::PosX, collision::PosY, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand_2prong::Px, hf_cand_2prong::Py>,
                  hf_cand::Ct<collision::PosX, collision::PosY, collision::PosZ, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex, hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::ImpactParameterXY<collision::PosX, collision::PosY, collision::PosZ, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex, hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand_2prong::MaxNormalisedDeltaIP<collision::PosX, collision::PosY, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ErrorDecayLengthXY, hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand::ImpactParameter0, hf_cand::ErrorImpactParameter0, hf_cand::ImpactParameter1, hf_cand::ErrorImpactParameter1, hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PxProng1, hf_cand::PyProng1>,
                  hf_cand::Eta<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::Phi<hf_cand_2prong::Px, hf_cand_2prong::Py>,
                  hf_cand::Y<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::E<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::E2<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>);

// extended table with expression columns that can be used as arguments of dynamic columns
DECLARE_SOA_EXTENDED_TABLE_USER(HfCandXiccExt, HfCandXiccBase, "HFCANDXICCEXT",
                                hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz);

using HfCandXicc = HfCandXiccExt;

// table with results of reconstruction level MC matching
DECLARE_SOA_TABLE(HfCandXiccMcRec, "AOD", "HFCANDXICCMCREC", //!
                  hf_cand_xicc::FlagMcMatchRec,
                  hf_cand_xicc::OriginMcRec,
                  hf_cand_xicc::DebugMcRec);

// table with results of generator level MC matching
DECLARE_SOA_TABLE(HfCandXiccMcGen, "AOD", "HFCANDXICCMCGEN", //!
                  hf_cand_xicc::FlagMcMatchGen,
                  hf_cand_xicc::OriginMcGen);

// specific Omegac and Xic to Xi Pi candidate properties
namespace hf_cand_xic0_omegac0
{
// Data processing results:
DECLARE_SOA_INDEX_COLUMN(Collision, collision);
DECLARE_SOA_COLUMN(XPv, xPv, float);
DECLARE_SOA_COLUMN(YPv, yPv, float);
DECLARE_SOA_COLUMN(ZPv, zPv, float);
DECLARE_SOA_COLUMN(XDecayVtxCharmBaryon, xDecayVtxCharmBaryon, float);
DECLARE_SOA_COLUMN(YDecayVtxCharmBaryon, yDecayVtxCharmBaryon, float);
DECLARE_SOA_COLUMN(ZDecayVtxCharmBaryon, zDecayVtxCharmBaryon, float);
DECLARE_SOA_COLUMN(XDecayVtxCascade, xDecayVtxCascade, float);
DECLARE_SOA_COLUMN(YDecayVtxCascade, yDecayVtxCascade, float);
DECLARE_SOA_COLUMN(ZDecayVtxCascade, zDecayVtxCascade, float);
DECLARE_SOA_COLUMN(XDecayVtxV0, xDecayVtxV0, float);
DECLARE_SOA_COLUMN(YDecayVtxV0, yDecayVtxV0, float);
DECLARE_SOA_COLUMN(ZDecayVtxV0, zDecayVtxV0, float);
DECLARE_SOA_COLUMN(SignDecay, signDecay, int8_t); // sign of pi <- xi
DECLARE_SOA_COLUMN(CovVtxCharmBaryon0, covVtxCharmBaryon0, float);
DECLARE_SOA_COLUMN(CovVtxCharmBaryon1, covVtxCharmBaryon1, float);
DECLARE_SOA_COLUMN(CovVtxCharmBaryon2, covVtxCharmBaryon2, float);
DECLARE_SOA_COLUMN(CovVtxCharmBaryon3, covVtxCharmBaryon3, float);
DECLARE_SOA_COLUMN(CovVtxCharmBaryon4, covVtxCharmBaryon4, float);
DECLARE_SOA_COLUMN(CovVtxCharmBaryon5, covVtxCharmBaryon5, float);
DECLARE_SOA_COLUMN(PxCharmBaryon, pxCharmBaryon, float);
DECLARE_SOA_COLUMN(PyCharmBaryon, pyCharmBaryon, float);
DECLARE_SOA_COLUMN(PzCharmBaryon, pzCharmBaryon, float);
DECLARE_SOA_COLUMN(PxCasc, pxCasc, float);
DECLARE_SOA_COLUMN(PyCasc, pyCasc, float);
DECLARE_SOA_COLUMN(PzCasc, pzCasc, float);
DECLARE_SOA_COLUMN(PxBachFromCharmBaryon, pxBachFromCharmBaryon, float);
DECLARE_SOA_COLUMN(PyBachFromCharmBaryon, pyBachFromCharmBaryon, float);
DECLARE_SOA_COLUMN(PzBachFromCharmBaryon, pzBachFromCharmBaryon, float);
DECLARE_SOA_COLUMN(PxLambda, pxLambda, float);
DECLARE_SOA_COLUMN(PyLambda, pyLambda, float);
DECLARE_SOA_COLUMN(PzLambda, pzLambda, float);
DECLARE_SOA_COLUMN(PxBachFromCasc, pxBachFromCasc, float);
DECLARE_SOA_COLUMN(PyBachFromCasc, pyBachFromCasc, float);
DECLARE_SOA_COLUMN(PzBachFromCasc, pzBachFromCasc, float);
DECLARE_SOA_COLUMN(PxPosV0Dau, pxPosV0Dau, float);
DECLARE_SOA_COLUMN(PyPosV0Dau, pyPosV0Dau, float);
DECLARE_SOA_COLUMN(PzPosV0Dau, pzPosV0Dau, float);
DECLARE_SOA_COLUMN(PxNegV0Dau, pxNegV0Dau, float);
DECLARE_SOA_COLUMN(PyNegV0Dau, pyNegV0Dau, float);
DECLARE_SOA_COLUMN(PzNegV0Dau, pzNegV0Dau, float);
DECLARE_SOA_COLUMN(ImpactParCascXY, impactParCascXY, float);
DECLARE_SOA_COLUMN(ImpactParBachFromCharmBaryonXY, impactParBachFromCharmBaryonXY, float);
DECLARE_SOA_COLUMN(ImpactParCascZ, impactParCascZ, float);
DECLARE_SOA_COLUMN(ImpactParBachFromCharmBaryonZ, impactParBachFromCharmBaryonZ, float);
DECLARE_SOA_COLUMN(ErrImpactParCascXY, errImpactParCascXY, float);
DECLARE_SOA_COLUMN(ErrImpactParBachFromCharmBaryonXY, errImpactParBachFromCharmBaryonXY, float);
DECLARE_SOA_INDEX_COLUMN(V0, v0);
DECLARE_SOA_INDEX_COLUMN(Cascade, cascade);
DECLARE_SOA_INDEX_COLUMN_FULL(BachelorFromCharmBaryon, bachelorFromCharmBaryon, int, Tracks, "_bachelorfromcharmbaryon");
DECLARE_SOA_COLUMN(InvMassLambda, invMassLambda, float);
DECLARE_SOA_COLUMN(InvMassCascade, invMassCascade, float);
DECLARE_SOA_COLUMN(InvMassCharmBaryon, invMassCharmBaryon, float);
DECLARE_SOA_COLUMN(CosPAV0, cosPAV0, float);
DECLARE_SOA_COLUMN(CosPACharmBaryon, cosPACharmBaryon, float);
DECLARE_SOA_COLUMN(CosPACasc, cosPACasc, float);
DECLARE_SOA_COLUMN(CosPAXYV0, cosPAXYV0, float);
DECLARE_SOA_COLUMN(CosPAXYCharmBaryon, cosPAXYCharmBaryon, float);
DECLARE_SOA_COLUMN(CosPAXYCasc, cosPAXYCasc, float);
DECLARE_SOA_COLUMN(CTauOmegac, cTauOmegac, float);
DECLARE_SOA_COLUMN(CTauCascade, cTauCascade, float);
DECLARE_SOA_COLUMN(CTauV0, cTauV0, float);
DECLARE_SOA_COLUMN(CTauXic, cTauXic, float);
DECLARE_SOA_COLUMN(EtaV0PosDau, etaV0PosDau, float);
DECLARE_SOA_COLUMN(EtaV0NegDau, etaV0NegDau, float);
DECLARE_SOA_COLUMN(EtaBachFromCasc, etaBachFromCasc, float);
DECLARE_SOA_COLUMN(EtaBachFromCharmBaryon, etaBachFromCharmBaryon, float);
DECLARE_SOA_COLUMN(EtaCharmBaryon, etaCharmBaryon, float);
DECLARE_SOA_COLUMN(EtaCascade, etaCascade, float);
DECLARE_SOA_COLUMN(EtaV0, etaV0, float);
DECLARE_SOA_COLUMN(DcaXYToPvV0Dau0, dcaXYToPvV0Dau0, float); // pos dau
DECLARE_SOA_COLUMN(DcaXYToPvV0Dau1, dcaXYToPvV0Dau1, float); // neg dau
DECLARE_SOA_COLUMN(DcaXYToPvCascDau, dcaXYToPvCascDau, float);
DECLARE_SOA_COLUMN(DcaZToPvV0Dau0, dcaZToPvV0Dau0, float);
DECLARE_SOA_COLUMN(DcaZToPvV0Dau1, dcaZToPvV0Dau1, float);
DECLARE_SOA_COLUMN(DcaZToPvCascDau, dcaZToPvCascDau, float);
DECLARE_SOA_COLUMN(DcaCascDau, dcaCascDau, float);
DECLARE_SOA_COLUMN(DcaV0Dau, dcaV0Dau, float);
DECLARE_SOA_COLUMN(DcaCharmBaryonDau, dcaCharmBaryonDau, float);
DECLARE_SOA_COLUMN(DecLenCharmBaryon, decLenCharmBaryon, float);
DECLARE_SOA_COLUMN(DecLenCascade, decLenCascade, float);
DECLARE_SOA_COLUMN(DecLenV0, decLenV0, float);
DECLARE_SOA_COLUMN(ErrorDecayLengthCharmBaryon, errorDecayLengthCharmBaryon, float);
DECLARE_SOA_COLUMN(ErrorDecayLengthXYCharmBaryon, errorDecayLengthXYCharmBaryon, float);

// KFParticle results
DECLARE_SOA_COLUMN(XPvKf, xPvKf, float);
DECLARE_SOA_COLUMN(YPvKf, yPvKf, float);
DECLARE_SOA_COLUMN(ZPvKf, zPvKf, float);
DECLARE_SOA_COLUMN(XDecayVtxV0Kf, xDecayVtxV0Kf, float);
DECLARE_SOA_COLUMN(YDecayVtxV0Kf, yDecayVtxV0Kf, float);
DECLARE_SOA_COLUMN(ZDecayVtxV0Kf, zDecayVtxV0Kf, float);
DECLARE_SOA_COLUMN(XDecayVtxCascadeKf, xDecayVtxCascadeKf, float);
DECLARE_SOA_COLUMN(YDecayVtxCascadeKf, yDecayVtxCascadeKf, float);
DECLARE_SOA_COLUMN(ZDecayVtxCascadeKf, zDecayVtxCascadeKf, float);
DECLARE_SOA_COLUMN(PxLambdaKf, pxLambdaKf, float);
DECLARE_SOA_COLUMN(PyLambdaKf, pyLambdaKf, float);
DECLARE_SOA_COLUMN(PzLambdaKf, pzLambdaKf, float);
DECLARE_SOA_COLUMN(PxCascKf, pxCascKf, float);
DECLARE_SOA_COLUMN(PyCascKf, pyCascKf, float);
DECLARE_SOA_COLUMN(PzCascKf, pzCascKf, float);
DECLARE_SOA_COLUMN(XDecayVtxOmegaKaKf, xDecayVtxOmegaKaKf, float);
DECLARE_SOA_COLUMN(YDecayVtxOmegaKaKf, yDecayVtxOmegaKaKf, float);
DECLARE_SOA_COLUMN(ZDecayVtxOmegaKaKf, zDecayVtxOmegaKaKf, float);
DECLARE_SOA_COLUMN(PxOmegaKaKf, pxOmegaKaKf, float);
DECLARE_SOA_COLUMN(PyOmegaKaKf, pyOmegaKaKf, float);
DECLARE_SOA_COLUMN(PzOmegaKaKf, pzOmegaKaKf, float);
DECLARE_SOA_COLUMN(EtaV0DauPr, etaV0DauPr, float);
DECLARE_SOA_COLUMN(EtaV0DauPi, etaV0DauPi, float);
DECLARE_SOA_COLUMN(InvMassCascadeRej, invMassCascadeRej, float);
DECLARE_SOA_COLUMN(InvMassLambdaErr, invMassLambdaErr, float);
DECLARE_SOA_COLUMN(InvMassCascadeErr, invMassCascadeErr, float);
DECLARE_SOA_COLUMN(InvMassCascadeRejErr, invMassCascadeRejErr, float);
DECLARE_SOA_COLUMN(InvMassCharmBaryonErr, invMassCharmBaryonErr, float);
DECLARE_SOA_COLUMN(CTauOmegaKa, cTauOmegaKa, float);
DECLARE_SOA_COLUMN(Chi2GeoOmegaKa, chi2GeoOmegaKa, float);
DECLARE_SOA_COLUMN(OmegaKaldl, omegaKaldl, float);
DECLARE_SOA_COLUMN(Chi2TopoKaFromOmegaKaToPv, chi2TopoKaFromOmegaKaToPv, float);
DECLARE_SOA_COLUMN(Chi2TopoOmegaKaToPv, chi2TopoOmegaKaToPv, float);
DECLARE_SOA_COLUMN(Chi2TopoKaToCasc, chi2TopoKaToCasc, float);
DECLARE_SOA_COLUMN(Chi2TopoCascToOmegaKa, chi2TopoCascToOmegaKa, float);
DECLARE_SOA_COLUMN(Chi2TopoKaToOmegaKa, chi2TopoKaToOmegaKa, float);
DECLARE_SOA_COLUMN(CosPaCascToOmegaKa, cosPaCascToOmegaKa, float);
DECLARE_SOA_COLUMN(CosPaXYCascToOmegaKa, cosPaXYCascToOmegaKa, float);
DECLARE_SOA_COLUMN(KfRapOmegaKa, kfRapOmegaKa, float);
DECLARE_SOA_COLUMN(KfPtKaFromOmegaKa, kfPtKaFromOmegaKa, float);
DECLARE_SOA_COLUMN(KfPtOmega, kfPtOmega, float);
DECLARE_SOA_COLUMN(KfPtOmegaKa, kfPtOmegaKa, float);
DECLARE_SOA_COLUMN(CosThetaStarKaFromOmegac, cosThetaStarKaFromOmegac, float);
DECLARE_SOA_COLUMN(CosThetaStarKaFromXic, cosThetaStarKaFromXic, float);
DECLARE_SOA_COLUMN(KfDcaXYPiFromOmegac, kfDcaXYPiFromOmegac, float);
DECLARE_SOA_COLUMN(KfDcaXYPiFromXic, kfDcaXYPiFromXic, float);
DECLARE_SOA_COLUMN(KfDcaXYCascToPv, kfDcaXYCascToPv, float);
DECLARE_SOA_COLUMN(Chi2GeoV0, chi2GeoV0, float);
DECLARE_SOA_COLUMN(Chi2GeoCasc, chi2GeoCasc, float);
DECLARE_SOA_COLUMN(Chi2GeoOmegac, chi2GeoOmegac, float);
DECLARE_SOA_COLUMN(Chi2GeoXic, chi2GeoXic, float);
DECLARE_SOA_COLUMN(Chi2MassV0, chi2MassV0, float);
DECLARE_SOA_COLUMN(Chi2MassCasc, chi2MassCasc, float);
DECLARE_SOA_COLUMN(V0ldl, v0ldl, float);
DECLARE_SOA_COLUMN(Cascldl, cascldl, float);
DECLARE_SOA_COLUMN(Omegacldl, omegacldl, float);
DECLARE_SOA_COLUMN(Xicldl, xicldl, float);
DECLARE_SOA_COLUMN(Chi2TopoV0ToPv, chi2TopoV0ToPv, float);
DECLARE_SOA_COLUMN(Chi2TopoCascToPv, chi2TopoCascToPv, float);
DECLARE_SOA_COLUMN(Chi2TopoPiFromOmegacToPv, chi2TopoPiFromOmegacToPv, float);
DECLARE_SOA_COLUMN(Chi2TopoPiFromXicToPv, chi2TopoPiFromXicToPv, float);
DECLARE_SOA_COLUMN(Chi2TopoOmegacToPv, chi2TopoOmegacToPv, float);
DECLARE_SOA_COLUMN(DeviationPiFromOmegacToPv, deviationPiFromOmegacToPv, float);
DECLARE_SOA_COLUMN(Chi2TopoXicToPv, chi2TopoXicToPv, float);
DECLARE_SOA_COLUMN(Chi2TopoV0ToCasc, chi2TopoV0ToCasc, float);
DECLARE_SOA_COLUMN(Chi2TopoCascToOmegac, chi2TopoCascToOmegac, float);
DECLARE_SOA_COLUMN(Chi2TopoCascToXic, chi2TopoCascToXic, float);
DECLARE_SOA_COLUMN(DecayLenXYLambda, decayLenXYLambda, float);
DECLARE_SOA_COLUMN(DecayLenXYCasc, decayLenXYCasc, float);
DECLARE_SOA_COLUMN(DecayLenXYOmegac, decayLenXYOmegac, float);
DECLARE_SOA_COLUMN(DecayLenXYXic, decayLenXYXic, float);
DECLARE_SOA_COLUMN(CosPaV0ToCasc, cosPaV0ToCasc, float);
DECLARE_SOA_COLUMN(CosPaCascToOmegac, cosPaCascToOmegac, float);
DECLARE_SOA_COLUMN(CosPaCascToXic, cosPaCascToXic, float);
DECLARE_SOA_COLUMN(CosPaXYV0ToCasc, cosPaXYV0ToCasc, float);
DECLARE_SOA_COLUMN(CosPaXYCascToOmegac, cosPaXYCascToOmegac, float);
DECLARE_SOA_COLUMN(CosPaXYCascToXic, cosPaXYCascToXic, float);
DECLARE_SOA_COLUMN(KfMassV0, kfMassV0, float);
DECLARE_SOA_COLUMN(KfMassCasc, kfMassCasc, float);
DECLARE_SOA_COLUMN(KfMassOmegac, kfMassOmegac, float);
DECLARE_SOA_COLUMN(KfRapOmegac, kfRapOmegac, float);
DECLARE_SOA_COLUMN(KfRapXic, kfRapXic, float);
DECLARE_SOA_COLUMN(KfptPiFromOmegac, kfptPiFromOmegac, float);
DECLARE_SOA_COLUMN(KfptOmegac, kfptOmegac, float);
DECLARE_SOA_COLUMN(CosThetaStarPiFromOmegac, cosThetaStarPiFromOmegac, float);
DECLARE_SOA_COLUMN(CosThetaStarPiFromXic, cosThetaStarPiFromXic, float);
DECLARE_SOA_COLUMN(V0Ndf, v0Ndf, float);
DECLARE_SOA_COLUMN(CascNdf, cascNdf, float);
DECLARE_SOA_COLUMN(OmegacNdf, omegacNdf, float);
DECLARE_SOA_COLUMN(XicNdf, xicNdf, float);
DECLARE_SOA_COLUMN(MassV0Ndf, massV0Ndf, float);
DECLARE_SOA_COLUMN(MassCascNdf, massCascNdf, float);
DECLARE_SOA_COLUMN(V0Chi2OverNdf, v0Chi2OverNdf, float);
DECLARE_SOA_COLUMN(CascChi2OverNdf, cascChi2OverNdf, float);
DECLARE_SOA_COLUMN(OmegacChi2OverNdf, omegacChi2OverNdf, float);
DECLARE_SOA_COLUMN(XicChi2OverNdf, xicChi2OverNdf, float);
DECLARE_SOA_COLUMN(MassV0Chi2OverNdf, massV0Chi2OverNdf, float);
DECLARE_SOA_COLUMN(MassCascChi2OverNdf, massCascChi2OverNdf, float);
DECLARE_SOA_COLUMN(CascRejectInvmass, cascRejectInvmass, float);

// Kf QA results:
DECLARE_SOA_COLUMN(InvMassV0Err, invMassV0Err, float);
DECLARE_SOA_COLUMN(InvMassXiErr, invMassXiErr, float);
DECLARE_SOA_COLUMN(InvMassXic0Err, invMassXic0Err, float);
DECLARE_SOA_COLUMN(V0DauPosX, v0DauPosX, float);
DECLARE_SOA_COLUMN(V0DauPosY, v0DauPosY, float);
DECLARE_SOA_COLUMN(V0DauPosZ, v0DauPosZ, float);
DECLARE_SOA_COLUMN(V0DauPosXError, v0DauPosXError, float);
DECLARE_SOA_COLUMN(V0DauPosYError, v0DauPosYError, float);
DECLARE_SOA_COLUMN(V0DauPosZError, v0DauPosZError, float);
DECLARE_SOA_COLUMN(V0DauPosPt, v0DauPosPt, float);
DECLARE_SOA_COLUMN(V0DauNegX, v0DauNegX, float);
DECLARE_SOA_COLUMN(V0DauNegY, v0DauNegY, float);
DECLARE_SOA_COLUMN(V0DauNegZ, v0DauNegZ, float);
DECLARE_SOA_COLUMN(V0DauNegXError, v0DauNegXError, float);
DECLARE_SOA_COLUMN(V0DauNegYError, v0DauNegYError, float);
DECLARE_SOA_COLUMN(V0DauNegZError, v0DauNegZError, float);
DECLARE_SOA_COLUMN(V0DauNegPt, v0DauNegPt, float);

DECLARE_SOA_COLUMN(V0VtxX, v0VtxX, float);
DECLARE_SOA_COLUMN(V0VtxY, v0VtxY, float);
DECLARE_SOA_COLUMN(V0VtxZ, v0VtxZ, float);
DECLARE_SOA_COLUMN(V0XError, v0XError, float);
DECLARE_SOA_COLUMN(V0YError, v0YError, float);
DECLARE_SOA_COLUMN(V0ZError, v0ZError, float);
DECLARE_SOA_COLUMN(V0Pt, v0Pt, float);
DECLARE_SOA_COLUMN(XiBachelorX, xiBachelorX, float);
DECLARE_SOA_COLUMN(XiBachelorY, xiBachelorY, float);
DECLARE_SOA_COLUMN(XiBachelorZ, xiBachelorZ, float);
DECLARE_SOA_COLUMN(XiBachelorPt, xiBachelorPt, float);
DECLARE_SOA_COLUMN(XiBachelorXError, xiBachelorXError, float);
DECLARE_SOA_COLUMN(XiBachelorYError, xiBachelorYError, float);
DECLARE_SOA_COLUMN(XiBachelorZError, xiBachelorZError, float);
DECLARE_SOA_COLUMN(XiX, xiX, float);
DECLARE_SOA_COLUMN(XiY, xiY, float);
DECLARE_SOA_COLUMN(XiZ, xiZ, float);
DECLARE_SOA_COLUMN(XiXError, xiXError, float);
DECLARE_SOA_COLUMN(XiYError, xiYError, float);
DECLARE_SOA_COLUMN(XiZError, xiZError, float);
DECLARE_SOA_COLUMN(XiPt, xiPt, float);
DECLARE_SOA_COLUMN(Xic0BachelorX, xic0BachelorX, float);
DECLARE_SOA_COLUMN(Xic0BachelorY, xic0BachelorY, float);
DECLARE_SOA_COLUMN(Xic0BachelorZ, xic0BachelorZ, float);
DECLARE_SOA_COLUMN(Xic0BachelorPt, xic0BachelorPt, float);
DECLARE_SOA_COLUMN(Xic0BachelorXError, xic0BachelorXError, float);
DECLARE_SOA_COLUMN(Xic0BachelorYError, xic0BachelorYError, float);
DECLARE_SOA_COLUMN(Xic0BachelorZError, xic0BachelorZError, float);
DECLARE_SOA_COLUMN(Xic0Pt, xic0Pt, float);
DECLARE_SOA_COLUMN(Xic0XError, xic0XError, float);
DECLARE_SOA_COLUMN(Xic0YError, xic0YError, float);
DECLARE_SOA_COLUMN(Xic0ZError, xic0ZError, float);

// MC matching result:
DECLARE_SOA_COLUMN(FlagMcMatchRec, flagMcMatchRec, int8_t); // reconstruction level
DECLARE_SOA_COLUMN(DebugMcRec, debugMcRec, int8_t);         // debug flag for mis-association reconstruction level
DECLARE_SOA_COLUMN(FlagMcMatchGen, flagMcMatchGen, int8_t); // generator level
DECLARE_SOA_COLUMN(CollisionMatched, collisionMatched, bool);
DECLARE_SOA_COLUMN(DebugGenCharmBar, debugGenCharmBar, int8_t);
DECLARE_SOA_COLUMN(DebugGenCasc, debugGenCasc, int8_t);
DECLARE_SOA_COLUMN(DebugGenLambda, debugGenLambda, int8_t);
DECLARE_SOA_COLUMN(OriginMcRec, originMcRec, int8_t);
DECLARE_SOA_COLUMN(OriginMcGen, originMcGen, int8_t);
DECLARE_SOA_COLUMN(PtCharmBaryonGen, ptCharmBaryonGen, float);
DECLARE_SOA_COLUMN(RapidityCharmBaryonGen, rapidityCharmBaryonGen, float);

// dynamic columns
DECLARE_SOA_DYNAMIC_COLUMN(PtCharmBaryon, ptCharmBaryon,
                           [](float px, float py) -> float { return RecoDecay::pt(px, py); });
DECLARE_SOA_DYNAMIC_COLUMN(PtCasc, ptCasc,
                           [](float px, float py) -> float { return RecoDecay::pt(px, py); });
DECLARE_SOA_DYNAMIC_COLUMN(PtLambda, ptLambda,
                           [](float px, float py) -> float { return RecoDecay::pt(px, py); });
DECLARE_SOA_DYNAMIC_COLUMN(PtPiFromCharmBaryon, ptPiFromCharmBaryon,
                           [](float px, float py) -> float { return RecoDecay::pt(px, py); });
DECLARE_SOA_DYNAMIC_COLUMN(PtKaFromCasc, ptKaFromCasc,
                           [](float px, float py) -> float { return RecoDecay::pt(px, py); });

} // end of namespace hf_cand_xic0_omegac0

// declare dedicated Omegac and Xic to Xi Pi candidate table
DECLARE_SOA_TABLE(HfCandToXiPi, "AOD", "HFCANDTOXIPI",
                  o2::soa::Index<>,
                  hf_cand_xic0_omegac0::CollisionId, hf_cand_xic0_omegac0::XPv, hf_cand_xic0_omegac0::YPv, hf_cand_xic0_omegac0::ZPv,
                  hf_cand_xic0_omegac0::XDecayVtxCharmBaryon, hf_cand_xic0_omegac0::YDecayVtxCharmBaryon, hf_cand_xic0_omegac0::ZDecayVtxCharmBaryon,
                  hf_cand_xic0_omegac0::XDecayVtxCascade, hf_cand_xic0_omegac0::YDecayVtxCascade, hf_cand_xic0_omegac0::ZDecayVtxCascade,
                  hf_cand_xic0_omegac0::XDecayVtxV0, hf_cand_xic0_omegac0::YDecayVtxV0, hf_cand_xic0_omegac0::ZDecayVtxV0,
                  hf_cand_xic0_omegac0::SignDecay, // charge pi<-cascade (neg -> omegac, pos -> antiomegac)
                  hf_cand_xic0_omegac0::CovVtxCharmBaryon0, hf_cand_xic0_omegac0::CovVtxCharmBaryon1, hf_cand_xic0_omegac0::CovVtxCharmBaryon2, hf_cand_xic0_omegac0::CovVtxCharmBaryon3, hf_cand_xic0_omegac0::CovVtxCharmBaryon4, hf_cand_xic0_omegac0::CovVtxCharmBaryon5,
                  hf_cand_xic0_omegac0::PxCharmBaryon, hf_cand_xic0_omegac0::PyCharmBaryon, hf_cand_xic0_omegac0::PzCharmBaryon,
                  hf_cand_xic0_omegac0::PxCasc, hf_cand_xic0_omegac0::PyCasc, hf_cand_xic0_omegac0::PzCasc,
                  hf_cand_xic0_omegac0::PxBachFromCharmBaryon, hf_cand_xic0_omegac0::PyBachFromCharmBaryon, hf_cand_xic0_omegac0::PzBachFromCharmBaryon,
                  hf_cand_xic0_omegac0::PxLambda, hf_cand_xic0_omegac0::PyLambda, hf_cand_xic0_omegac0::PzLambda,
                  hf_cand_xic0_omegac0::PxBachFromCasc, hf_cand_xic0_omegac0::PyBachFromCasc, hf_cand_xic0_omegac0::PzBachFromCasc,
                  hf_cand_xic0_omegac0::PxPosV0Dau, hf_cand_xic0_omegac0::PyPosV0Dau, hf_cand_xic0_omegac0::PzPosV0Dau,
                  hf_cand_xic0_omegac0::PxNegV0Dau, hf_cand_xic0_omegac0::PyNegV0Dau, hf_cand_xic0_omegac0::PzNegV0Dau,
                  hf_cand_xic0_omegac0::ImpactParCascXY, hf_cand_xic0_omegac0::ImpactParBachFromCharmBaryonXY, hf_cand_xic0_omegac0::ImpactParCascZ, hf_cand_xic0_omegac0::ImpactParBachFromCharmBaryonZ,
                  hf_cand_xic0_omegac0::ErrImpactParCascXY, hf_cand_xic0_omegac0::ErrImpactParBachFromCharmBaryonXY,
                  hf_cand_xic0_omegac0::V0Id, v0data::PosTrackId, v0data::NegTrackId, hf_cand_xic0_omegac0::CascadeId, hf_cand_xic0_omegac0::BachelorFromCharmBaryonId, cascdata::BachelorId,
                  hf_cand_xic0_omegac0::InvMassLambda, hf_cand_xic0_omegac0::InvMassCascade, hf_cand_xic0_omegac0::InvMassCharmBaryon,
                  hf_cand_xic0_omegac0::CosPAV0, hf_cand_xic0_omegac0::CosPACharmBaryon, hf_cand_xic0_omegac0::CosPACasc, hf_cand_xic0_omegac0::CosPAXYV0, hf_cand_xic0_omegac0::CosPAXYCharmBaryon, hf_cand_xic0_omegac0::CosPAXYCasc,
                  hf_cand_xic0_omegac0::CTauOmegac, hf_cand_xic0_omegac0::CTauCascade, hf_cand_xic0_omegac0::CTauV0, hf_cand_xic0_omegac0::CTauXic,
                  hf_cand_xic0_omegac0::EtaV0PosDau, hf_cand_xic0_omegac0::EtaV0NegDau, hf_cand_xic0_omegac0::EtaBachFromCasc, hf_cand_xic0_omegac0::EtaBachFromCharmBaryon,
                  hf_cand_xic0_omegac0::EtaCharmBaryon, hf_cand_xic0_omegac0::EtaCascade, hf_cand_xic0_omegac0::EtaV0,
                  hf_cand_xic0_omegac0::DcaXYToPvV0Dau0, hf_cand_xic0_omegac0::DcaXYToPvV0Dau1, hf_cand_xic0_omegac0::DcaXYToPvCascDau,
                  hf_cand_xic0_omegac0::DcaZToPvV0Dau0, hf_cand_xic0_omegac0::DcaZToPvV0Dau1, hf_cand_xic0_omegac0::DcaZToPvCascDau,
                  hf_cand_xic0_omegac0::DcaCascDau, hf_cand_xic0_omegac0::DcaV0Dau, hf_cand_xic0_omegac0::DcaCharmBaryonDau,
                  hf_cand_xic0_omegac0::DecLenCharmBaryon, hf_cand_xic0_omegac0::DecLenCascade, hf_cand_xic0_omegac0::DecLenV0, hf_cand_xic0_omegac0::ErrorDecayLengthCharmBaryon, hf_cand_xic0_omegac0::ErrorDecayLengthXYCharmBaryon,
                  // dynamic
                  hf_cand::Y<hf_cand_xic0_omegac0::PxCharmBaryon, hf_cand_xic0_omegac0::PyCharmBaryon, hf_cand_xic0_omegac0::PzCharmBaryon>);

DECLARE_SOA_TABLE(HfCandToOmegaPi, "AOD", "HFCANDTOOMEGAPI",
                  o2::soa::Index<>,
                  hf_cand_xic0_omegac0::CollisionId, hf_cand_xic0_omegac0::XPv, hf_cand_xic0_omegac0::YPv, hf_cand_xic0_omegac0::ZPv,
                  hf_cand_xic0_omegac0::XDecayVtxCharmBaryon, hf_cand_xic0_omegac0::YDecayVtxCharmBaryon, hf_cand_xic0_omegac0::ZDecayVtxCharmBaryon,
                  hf_cand_xic0_omegac0::XDecayVtxCascade, hf_cand_xic0_omegac0::YDecayVtxCascade, hf_cand_xic0_omegac0::ZDecayVtxCascade,
                  hf_cand_xic0_omegac0::XDecayVtxV0, hf_cand_xic0_omegac0::YDecayVtxV0, hf_cand_xic0_omegac0::ZDecayVtxV0,
                  hf_cand_xic0_omegac0::SignDecay, // charge pi<-cascade (neg -> omegac, pos -> antiomegac)
                  hf_cand_xic0_omegac0::CovVtxCharmBaryon0, hf_cand_xic0_omegac0::CovVtxCharmBaryon1, hf_cand_xic0_omegac0::CovVtxCharmBaryon2, hf_cand_xic0_omegac0::CovVtxCharmBaryon3, hf_cand_xic0_omegac0::CovVtxCharmBaryon4, hf_cand_xic0_omegac0::CovVtxCharmBaryon5,
                  hf_cand_xic0_omegac0::PxCharmBaryon, hf_cand_xic0_omegac0::PyCharmBaryon, hf_cand_xic0_omegac0::PzCharmBaryon,
                  hf_cand_xic0_omegac0::PxCasc, hf_cand_xic0_omegac0::PyCasc, hf_cand_xic0_omegac0::PzCasc,
                  hf_cand_xic0_omegac0::PxBachFromCharmBaryon, hf_cand_xic0_omegac0::PyBachFromCharmBaryon, hf_cand_xic0_omegac0::PzBachFromCharmBaryon,
                  hf_cand_xic0_omegac0::PxLambda, hf_cand_xic0_omegac0::PyLambda, hf_cand_xic0_omegac0::PzLambda,
                  hf_cand_xic0_omegac0::PxBachFromCasc, hf_cand_xic0_omegac0::PyBachFromCasc, hf_cand_xic0_omegac0::PzBachFromCasc,
                  hf_cand_xic0_omegac0::PxPosV0Dau, hf_cand_xic0_omegac0::PyPosV0Dau, hf_cand_xic0_omegac0::PzPosV0Dau,
                  hf_cand_xic0_omegac0::PxNegV0Dau, hf_cand_xic0_omegac0::PyNegV0Dau, hf_cand_xic0_omegac0::PzNegV0Dau,
                  // dynamic
                  hf_cand_xic0_omegac0::PtCharmBaryon<hf_cand_xic0_omegac0::PxCharmBaryon, hf_cand_xic0_omegac0::PyCharmBaryon>,
                  hf_cand_xic0_omegac0::PtCasc<hf_cand_xic0_omegac0::PxCasc, hf_cand_xic0_omegac0::PyCasc>,
                  hf_cand_xic0_omegac0::PtPiFromCharmBaryon<hf_cand_xic0_omegac0::PxBachFromCharmBaryon, hf_cand_xic0_omegac0::PyBachFromCharmBaryon>,
                  hf_cand_xic0_omegac0::PtLambda<hf_cand_xic0_omegac0::PxLambda, hf_cand_xic0_omegac0::PyLambda>,
                  hf_cand_xic0_omegac0::PtKaFromCasc<hf_cand_xic0_omegac0::PxBachFromCasc, hf_cand_xic0_omegac0::PyBachFromCasc>,

                  hf_cand_xic0_omegac0::ImpactParCascXY, hf_cand_xic0_omegac0::ImpactParBachFromCharmBaryonXY, hf_cand_xic0_omegac0::ImpactParCascZ, hf_cand_xic0_omegac0::ImpactParBachFromCharmBaryonZ,
                  hf_cand_xic0_omegac0::ErrImpactParCascXY, hf_cand_xic0_omegac0::ErrImpactParBachFromCharmBaryonXY,
                  hf_cand_xic0_omegac0::V0Id, v0data::PosTrackId, v0data::NegTrackId, hf_cand_xic0_omegac0::CascadeId, hf_cand_xic0_omegac0::BachelorFromCharmBaryonId, cascdata::BachelorId,
                  hf_cand_xic0_omegac0::InvMassLambda, hf_cand_xic0_omegac0::InvMassCascade, hf_cand_xic0_omegac0::InvMassCharmBaryon,
                  hf_cand_xic0_omegac0::CosPAV0, hf_cand_xic0_omegac0::CosPACharmBaryon, hf_cand_xic0_omegac0::CosPACasc, hf_cand_xic0_omegac0::CosPAXYV0, hf_cand_xic0_omegac0::CosPAXYCharmBaryon, hf_cand_xic0_omegac0::CosPAXYCasc,
                  hf_cand_xic0_omegac0::CTauOmegac, hf_cand_xic0_omegac0::CTauCascade, hf_cand_xic0_omegac0::CTauV0,
                  hf_cand_xic0_omegac0::EtaV0PosDau, hf_cand_xic0_omegac0::EtaV0NegDau, hf_cand_xic0_omegac0::EtaBachFromCasc, hf_cand_xic0_omegac0::EtaBachFromCharmBaryon,
                  hf_cand_xic0_omegac0::EtaCharmBaryon, hf_cand_xic0_omegac0::EtaCascade, hf_cand_xic0_omegac0::EtaV0,
                  hf_cand_xic0_omegac0::DcaXYToPvV0Dau0, hf_cand_xic0_omegac0::DcaXYToPvV0Dau1, hf_cand_xic0_omegac0::DcaXYToPvCascDau,
                  hf_cand_xic0_omegac0::DcaZToPvV0Dau0, hf_cand_xic0_omegac0::DcaZToPvV0Dau1, hf_cand_xic0_omegac0::DcaZToPvCascDau,
                  hf_cand_xic0_omegac0::DcaCascDau, hf_cand_xic0_omegac0::DcaV0Dau, hf_cand_xic0_omegac0::DcaCharmBaryonDau,
                  hf_cand_xic0_omegac0::DecLenCharmBaryon, hf_cand_xic0_omegac0::DecLenCascade, hf_cand_xic0_omegac0::DecLenV0, hf_cand_xic0_omegac0::ErrorDecayLengthCharmBaryon, hf_cand_xic0_omegac0::ErrorDecayLengthXYCharmBaryon, hf_track_index::HFflag,
                  o2::soa::Marker<1>);

DECLARE_SOA_TABLE(HfCandToOmegaK, "AOD", "HFCANDTOOMEGAK",
                  o2::soa::Index<>,
                  hf_cand_xic0_omegac0::CollisionId, hf_cand_xic0_omegac0::XPv, hf_cand_xic0_omegac0::YPv, hf_cand_xic0_omegac0::ZPv,
                  hf_cand_xic0_omegac0::XDecayVtxCharmBaryon, hf_cand_xic0_omegac0::YDecayVtxCharmBaryon, hf_cand_xic0_omegac0::ZDecayVtxCharmBaryon,
                  hf_cand_xic0_omegac0::XDecayVtxCascade, hf_cand_xic0_omegac0::YDecayVtxCascade, hf_cand_xic0_omegac0::ZDecayVtxCascade,
                  hf_cand_xic0_omegac0::XDecayVtxV0, hf_cand_xic0_omegac0::YDecayVtxV0, hf_cand_xic0_omegac0::ZDecayVtxV0,
                  hf_cand_xic0_omegac0::SignDecay, // charge pi<-cascade (neg -> omegac, pos -> antiomegac)
                  hf_cand_xic0_omegac0::CovVtxCharmBaryon0, hf_cand_xic0_omegac0::CovVtxCharmBaryon1, hf_cand_xic0_omegac0::CovVtxCharmBaryon2, hf_cand_xic0_omegac0::CovVtxCharmBaryon3, hf_cand_xic0_omegac0::CovVtxCharmBaryon4, hf_cand_xic0_omegac0::CovVtxCharmBaryon5,
                  hf_cand_xic0_omegac0::PxCharmBaryon, hf_cand_xic0_omegac0::PyCharmBaryon, hf_cand_xic0_omegac0::PzCharmBaryon,
                  hf_cand_xic0_omegac0::PxCasc, hf_cand_xic0_omegac0::PyCasc, hf_cand_xic0_omegac0::PzCasc,
                  hf_cand_xic0_omegac0::PxBachFromCharmBaryon, hf_cand_xic0_omegac0::PyBachFromCharmBaryon, hf_cand_xic0_omegac0::PzBachFromCharmBaryon,
                  hf_cand_xic0_omegac0::PxLambda, hf_cand_xic0_omegac0::PyLambda, hf_cand_xic0_omegac0::PzLambda,
                  hf_cand_xic0_omegac0::PxBachFromCasc, hf_cand_xic0_omegac0::PyBachFromCasc, hf_cand_xic0_omegac0::PzBachFromCasc,
                  hf_cand_xic0_omegac0::PxPosV0Dau, hf_cand_xic0_omegac0::PyPosV0Dau, hf_cand_xic0_omegac0::PzPosV0Dau,
                  hf_cand_xic0_omegac0::PxNegV0Dau, hf_cand_xic0_omegac0::PyNegV0Dau, hf_cand_xic0_omegac0::PzNegV0Dau,
                  hf_cand_xic0_omegac0::ImpactParCascXY, hf_cand_xic0_omegac0::ImpactParBachFromCharmBaryonXY, hf_cand_xic0_omegac0::ImpactParCascZ, hf_cand_xic0_omegac0::ImpactParBachFromCharmBaryonZ,
                  hf_cand_xic0_omegac0::ErrImpactParCascXY, hf_cand_xic0_omegac0::ErrImpactParBachFromCharmBaryonXY,
                  hf_cand_xic0_omegac0::V0Id, v0data::PosTrackId, v0data::NegTrackId, hf_cand_xic0_omegac0::CascadeId, hf_cand_xic0_omegac0::BachelorFromCharmBaryonId, cascdata::BachelorId,
                  hf_cand_xic0_omegac0::InvMassLambda, hf_cand_xic0_omegac0::InvMassCascade, hf_cand_xic0_omegac0::InvMassCharmBaryon,
                  hf_cand_xic0_omegac0::CosPAV0, hf_cand_xic0_omegac0::CosPACharmBaryon, hf_cand_xic0_omegac0::CosPACasc, hf_cand_xic0_omegac0::CosPAXYV0, hf_cand_xic0_omegac0::CosPAXYCharmBaryon, hf_cand_xic0_omegac0::CosPAXYCasc,
                  hf_cand_xic0_omegac0::CTauOmegac, hf_cand_xic0_omegac0::CTauCascade, hf_cand_xic0_omegac0::CTauV0,
                  hf_cand_xic0_omegac0::EtaV0PosDau, hf_cand_xic0_omegac0::EtaV0NegDau, hf_cand_xic0_omegac0::EtaBachFromCasc, hf_cand_xic0_omegac0::EtaBachFromCharmBaryon,
                  hf_cand_xic0_omegac0::EtaCharmBaryon, hf_cand_xic0_omegac0::EtaCascade, hf_cand_xic0_omegac0::EtaV0,
                  hf_cand_xic0_omegac0::DcaXYToPvV0Dau0, hf_cand_xic0_omegac0::DcaXYToPvV0Dau1, hf_cand_xic0_omegac0::DcaXYToPvCascDau,
                  hf_cand_xic0_omegac0::DcaZToPvV0Dau0, hf_cand_xic0_omegac0::DcaZToPvV0Dau1, hf_cand_xic0_omegac0::DcaZToPvCascDau,
                  hf_cand_xic0_omegac0::DcaCascDau, hf_cand_xic0_omegac0::DcaV0Dau, hf_cand_xic0_omegac0::DcaCharmBaryonDau,
                  hf_cand_xic0_omegac0::DecLenCharmBaryon, hf_cand_xic0_omegac0::DecLenCascade, hf_cand_xic0_omegac0::DecLenV0, hf_cand_xic0_omegac0::ErrorDecayLengthCharmBaryon, hf_cand_xic0_omegac0::ErrorDecayLengthXYCharmBaryon,
                  o2::soa::Marker<2>);

// table with results of KFParticle
DECLARE_SOA_TABLE(HfOmegacKf, "AOD", "HFOMEGACKF", //!
                  hf_cand_xic0_omegac0::KfDcaXYPiFromOmegac, hf_cand_xic0_omegac0::KfDcaXYCascToPv,
                  hf_cand_xic0_omegac0::Chi2GeoV0, hf_cand_xic0_omegac0::Chi2GeoCasc, hf_cand_xic0_omegac0::Chi2GeoOmegac,
                  hf_cand_xic0_omegac0::Chi2MassV0, hf_cand_xic0_omegac0::Chi2MassCasc,
                  hf_cand_xic0_omegac0::V0ldl, hf_cand_xic0_omegac0::Cascldl, hf_cand_xic0_omegac0::Omegacldl,
                  hf_cand_xic0_omegac0::Chi2TopoV0ToPv, hf_cand_xic0_omegac0::Chi2TopoCascToPv, hf_cand_xic0_omegac0::Chi2TopoPiFromOmegacToPv, hf_cand_xic0_omegac0::Chi2TopoOmegacToPv, hf_cand_xic0_omegac0::DeviationPiFromOmegacToPv,
                  hf_cand_xic0_omegac0::Chi2TopoV0ToCasc, hf_cand_xic0_omegac0::Chi2TopoCascToOmegac,
                  hf_cand_xic0_omegac0::DecayLenXYLambda, hf_cand_xic0_omegac0::DecayLenXYCasc, hf_cand_xic0_omegac0::DecayLenXYOmegac,
                  hf_cand_xic0_omegac0::CosPaV0ToCasc, hf_cand_xic0_omegac0::CosPaCascToOmegac, hf_cand_xic0_omegac0::CosPaXYV0ToCasc, hf_cand_xic0_omegac0::CosPaXYCascToOmegac,
                  hf_cand_xic0_omegac0::KfRapOmegac,
                  hf_cand_xic0_omegac0::KfptPiFromOmegac, hf_cand_xic0_omegac0::KfptOmegac,
                  hf_cand_xic0_omegac0::CosThetaStarPiFromOmegac,
                  hf_cand_xic0_omegac0::V0Ndf, hf_cand_xic0_omegac0::CascNdf, hf_cand_xic0_omegac0::OmegacNdf,
                  hf_cand_xic0_omegac0::MassV0Ndf, hf_cand_xic0_omegac0::MassCascNdf,
                  hf_cand_xic0_omegac0::V0Chi2OverNdf, hf_cand_xic0_omegac0::CascChi2OverNdf, hf_cand_xic0_omegac0::OmegacChi2OverNdf,
                  hf_cand_xic0_omegac0::MassV0Chi2OverNdf, hf_cand_xic0_omegac0::MassCascChi2OverNdf, hf_cand_xic0_omegac0::CascRejectInvmass);

// OmegaKa reconstruct by KFParticle
DECLARE_SOA_TABLE(HfCandToOmegaKaKf, "AOD", "HFCANDTOOMEGAKAKF",
                  o2::soa::Index<>,
                  hf_cand_xic0_omegac0::CollisionId, hf_cand_xic0_omegac0::XPv, hf_cand_xic0_omegac0::YPv, hf_cand_xic0_omegac0::ZPv,
                  hf_cand_xic0_omegac0::XPvKf, hf_cand_xic0_omegac0::YPvKf, hf_cand_xic0_omegac0::ZPvKf,
                  hf_cand_xic0_omegac0::XDecayVtxV0, hf_cand_xic0_omegac0::YDecayVtxV0, hf_cand_xic0_omegac0::ZDecayVtxV0,
                  hf_cand_xic0_omegac0::PxLambda, hf_cand_xic0_omegac0::PyLambda, hf_cand_xic0_omegac0::PzLambda,
                  hf_cand_xic0_omegac0::XDecayVtxCascade, hf_cand_xic0_omegac0::YDecayVtxCascade, hf_cand_xic0_omegac0::ZDecayVtxCascade,
                  hf_cand_xic0_omegac0::PxCasc, hf_cand_xic0_omegac0::PyCasc, hf_cand_xic0_omegac0::PzCasc,
                  hf_cand_xic0_omegac0::XDecayVtxV0Kf, hf_cand_xic0_omegac0::YDecayVtxV0Kf, hf_cand_xic0_omegac0::ZDecayVtxV0Kf,
                  hf_cand_xic0_omegac0::PxLambdaKf, hf_cand_xic0_omegac0::PyLambdaKf, hf_cand_xic0_omegac0::PzLambdaKf,
                  hf_cand_xic0_omegac0::XDecayVtxCascadeKf, hf_cand_xic0_omegac0::YDecayVtxCascadeKf, hf_cand_xic0_omegac0::ZDecayVtxCascadeKf,
                  hf_cand_xic0_omegac0::PxCascKf, hf_cand_xic0_omegac0::PyCascKf, hf_cand_xic0_omegac0::PzCascKf,
                  hf_cand_xic0_omegac0::XDecayVtxOmegaKaKf, hf_cand_xic0_omegac0::YDecayVtxOmegaKaKf, hf_cand_xic0_omegac0::ZDecayVtxOmegaKaKf,
                  hf_cand_xic0_omegac0::PxOmegaKaKf, hf_cand_xic0_omegac0::PyOmegaKaKf, hf_cand_xic0_omegac0::PzOmegaKaKf,
                  hf_cand_xic0_omegac0::SignDecay, hf_cand_xic0_omegac0::EtaV0DauPr, hf_cand_xic0_omegac0::EtaV0DauPi, hf_cand_xic0_omegac0::EtaBachFromCasc,
                  hf_cand_xic0_omegac0::EtaBachFromCharmBaryon, hf_cand_xic0_omegac0::EtaV0, hf_cand_xic0_omegac0::EtaCascade, hf_cand_xic0_omegac0::EtaCharmBaryon,
                  hf_cand_xic0_omegac0::KfRapOmegaKa, hf_cand_xic0_omegac0::ImpactParBachFromCharmBaryonXY, hf_cand_xic0_omegac0::ErrImpactParBachFromCharmBaryonXY, hf_cand_xic0_omegac0::ImpactParCascXY, hf_cand_xic0_omegac0::ErrImpactParCascXY,
                  hf_cand_xic0_omegac0::DcaV0Dau, hf_cand_xic0_omegac0::DcaCascDau, hf_cand_xic0_omegac0::DcaCharmBaryonDau,
                  hf_cand_xic0_omegac0::CosPAV0, hf_cand_xic0_omegac0::CosPACasc, hf_cand_xic0_omegac0::CosPACharmBaryon, hf_cand_xic0_omegac0::CosPAXYV0, hf_cand_xic0_omegac0::CosPAXYCasc, hf_cand_xic0_omegac0::CosPAXYCharmBaryon,
                  hf_cand_xic0_omegac0::CosPaV0ToCasc, hf_cand_xic0_omegac0::CosPaCascToOmegaKa, hf_cand_xic0_omegac0::CosPaXYV0ToCasc, hf_cand_xic0_omegac0::CosPaXYCascToOmegaKa,
                  hf_cand_xic0_omegac0::Chi2GeoV0, hf_cand_xic0_omegac0::Chi2GeoCasc, hf_cand_xic0_omegac0::Chi2GeoOmegaKa,
                  hf_cand_xic0_omegac0::MassV0Chi2OverNdf, hf_cand_xic0_omegac0::MassCascChi2OverNdf,
                  hf_cand_xic0_omegac0::Chi2TopoV0ToCasc, hf_cand_xic0_omegac0::Chi2TopoKaToCasc, hf_cand_xic0_omegac0::Chi2TopoKaToOmegaKa, hf_cand_xic0_omegac0::Chi2TopoCascToOmegaKa,
                  hf_cand_xic0_omegac0::Chi2TopoV0ToPv, hf_cand_xic0_omegac0::Chi2TopoCascToPv, hf_cand_xic0_omegac0::Chi2TopoKaFromOmegaKaToPv, hf_cand_xic0_omegac0::Chi2TopoOmegaKaToPv,
                  hf_cand_xic0_omegac0::V0ldl, hf_cand_xic0_omegac0::Cascldl, hf_cand_xic0_omegac0::OmegaKaldl,
                  hf_cand_xic0_omegac0::DecLenV0, hf_cand_xic0_omegac0::DecLenCascade, hf_cand_xic0_omegac0::DecLenCharmBaryon,
                  hf_cand_xic0_omegac0::InvMassLambda, hf_cand_xic0_omegac0::InvMassLambdaErr, hf_cand_xic0_omegac0::InvMassCascade, hf_cand_xic0_omegac0::InvMassCascadeErr, hf_cand_xic0_omegac0::InvMassCascadeRej, hf_cand_xic0_omegac0::InvMassCascadeRejErr, hf_cand_xic0_omegac0::InvMassCharmBaryon, hf_cand_xic0_omegac0::InvMassCharmBaryonErr,
                  hf_cand_xic0_omegac0::KfPtOmegaKa, hf_cand_xic0_omegac0::KfPtKaFromOmegaKa, hf_cand_xic0_omegac0::KfPtOmega,
                  hf_cand_xic0_omegac0::CosThetaStarKaFromOmegac, hf_cand_xic0_omegac0::CosThetaStarKaFromXic, hf_cand_xic0_omegac0::CTauV0, hf_cand_xic0_omegac0::CTauCascade, hf_cand_xic0_omegac0::CTauOmegaKa,
                  hf_cand_xic0_omegac0::V0Id, v0data::PosTrackId, v0data::NegTrackId, hf_cand_xic0_omegac0::CascadeId, cascdata::BachelorId, hf_cand_xic0_omegac0::BachelorFromCharmBaryonId);

DECLARE_SOA_TABLE(HfCandToXiPiKf, "AOD", "HFCANDTOXIPIKF", //!
                  o2::soa::Index<>,
                  hf_cand_xic0_omegac0::CollisionId, hf_cand_xic0_omegac0::XPv, hf_cand_xic0_omegac0::YPv, hf_cand_xic0_omegac0::ZPv,
                  hf_cand_xic0_omegac0::XDecayVtxCascade, hf_cand_xic0_omegac0::YDecayVtxCascade, hf_cand_xic0_omegac0::ZDecayVtxCascade,
                  hf_cand_xic0_omegac0::XDecayVtxV0, hf_cand_xic0_omegac0::YDecayVtxV0, hf_cand_xic0_omegac0::ZDecayVtxV0,
                  hf_cand_xic0_omegac0::SignDecay,
                  hf_cand_xic0_omegac0::CovVtxCharmBaryon0, hf_cand_xic0_omegac0::CovVtxCharmBaryon1, hf_cand_xic0_omegac0::CovVtxCharmBaryon2, hf_cand_xic0_omegac0::CovVtxCharmBaryon3, hf_cand_xic0_omegac0::CovVtxCharmBaryon4, hf_cand_xic0_omegac0::CovVtxCharmBaryon5,
                  hf_cand_xic0_omegac0::PxCharmBaryon, hf_cand_xic0_omegac0::PyCharmBaryon, hf_cand_xic0_omegac0::PzCharmBaryon,
                  hf_cand_xic0_omegac0::PxCasc, hf_cand_xic0_omegac0::PyCasc, hf_cand_xic0_omegac0::PzCasc,
                  hf_cand_xic0_omegac0::PxBachFromCharmBaryon, hf_cand_xic0_omegac0::PyBachFromCharmBaryon, hf_cand_xic0_omegac0::PzBachFromCharmBaryon,
                  hf_cand_xic0_omegac0::PxLambda, hf_cand_xic0_omegac0::PyLambda, hf_cand_xic0_omegac0::PzLambda,
                  hf_cand_xic0_omegac0::PxBachFromCasc, hf_cand_xic0_omegac0::PyBachFromCasc, hf_cand_xic0_omegac0::PzBachFromCasc,
                  hf_cand_xic0_omegac0::PxPosV0Dau, hf_cand_xic0_omegac0::PyPosV0Dau, hf_cand_xic0_omegac0::PzPosV0Dau,
                  hf_cand_xic0_omegac0::PxNegV0Dau, hf_cand_xic0_omegac0::PyNegV0Dau, hf_cand_xic0_omegac0::PzNegV0Dau,
                  hf_cand_xic0_omegac0::V0Id, v0data::PosTrackId, v0data::NegTrackId, hf_cand_xic0_omegac0::CascadeId, hf_cand_xic0_omegac0::BachelorFromCharmBaryonId, cascdata::BachelorId,
                  hf_cand_xic0_omegac0::InvMassLambda, hf_cand_xic0_omegac0::InvMassCascade, hf_cand_xic0_omegac0::InvMassCharmBaryon,
                  hf_cand_xic0_omegac0::CosPAV0, hf_cand_xic0_omegac0::CosPACasc,
                  hf_cand_xic0_omegac0::CTauCascade, hf_cand_xic0_omegac0::CTauV0, hf_cand_xic0_omegac0::CTauXic,
                  hf_cand_xic0_omegac0::EtaV0PosDau, hf_cand_xic0_omegac0::EtaV0NegDau, hf_cand_xic0_omegac0::EtaBachFromCasc, hf_cand_xic0_omegac0::EtaBachFromCharmBaryon,
                  hf_cand_xic0_omegac0::EtaCharmBaryon, hf_cand_xic0_omegac0::EtaCascade, hf_cand_xic0_omegac0::EtaV0,
                  hf_cand_xic0_omegac0::DcaXYToPvV0Dau0, hf_cand_xic0_omegac0::DcaXYToPvV0Dau1, hf_cand_xic0_omegac0::DcaXYToPvCascDau,
                  hf_cand_xic0_omegac0::DcaCascDau, hf_cand_xic0_omegac0::DcaV0Dau, hf_cand_xic0_omegac0::DcaCharmBaryonDau,
                  hf_cand_xic0_omegac0::KfDcaXYPiFromXic, hf_cand_xic0_omegac0::KfDcaXYCascToPv,
                  hf_cand_xic0_omegac0::Chi2GeoV0, hf_cand_xic0_omegac0::Chi2GeoCasc, hf_cand_xic0_omegac0::Chi2GeoXic,
                  hf_cand_xic0_omegac0::Chi2MassV0, hf_cand_xic0_omegac0::Chi2MassCasc,
                  hf_cand_xic0_omegac0::V0ldl, hf_cand_xic0_omegac0::Cascldl,
                  hf_cand_xic0_omegac0::Chi2TopoV0ToPv, hf_cand_xic0_omegac0::Chi2TopoCascToPv, hf_cand_xic0_omegac0::Chi2TopoPiFromXicToPv, hf_cand_xic0_omegac0::Chi2TopoXicToPv,
                  hf_cand_xic0_omegac0::Chi2TopoV0ToCasc, hf_cand_xic0_omegac0::Chi2TopoCascToXic,
                  hf_cand_xic0_omegac0::DecayLenXYLambda, hf_cand_xic0_omegac0::DecayLenXYCasc, hf_cand_xic0_omegac0::DecayLenXYXic,
                  hf_cand_xic0_omegac0::CosPaV0ToCasc, hf_cand_xic0_omegac0::CosPaCascToXic,
                  hf_cand_xic0_omegac0::KfRapXic,
                  hf_cand_xic0_omegac0::CosThetaStarPiFromXic,
                  hf_cand_xic0_omegac0::V0Ndf, hf_cand_xic0_omegac0::CascNdf, hf_cand_xic0_omegac0::XicNdf,
                  hf_cand_xic0_omegac0::MassV0Ndf, hf_cand_xic0_omegac0::MassCascNdf,
                  hf_cand_xic0_omegac0::V0Chi2OverNdf, hf_cand_xic0_omegac0::CascChi2OverNdf, hf_cand_xic0_omegac0::XicChi2OverNdf,
                  hf_cand_xic0_omegac0::MassV0Chi2OverNdf, hf_cand_xic0_omegac0::MassCascChi2OverNdf);

DECLARE_SOA_TABLE(HfCandToXiPiKfQa, "AOD", "HFCANDTOXIPIKFQA",
                  o2::soa::Index<>,
                  hf_cand_xic0_omegac0::InvMassLambda, hf_cand_xic0_omegac0::InvMassCascade, hf_cand_xic0_omegac0::InvMassCharmBaryon, hf_cand_xic0_omegac0::InvMassV0Err, hf_cand_xic0_omegac0::InvMassXiErr, hf_cand_xic0_omegac0::InvMassXic0Err,
                  hf_cand_xic0_omegac0::CollisionId, hf_track_index::V0Id, v0data::PosTrackId, v0data::NegTrackId, hf_cand_xic0_omegac0::CascadeId, hf_cand_xic0_omegac0::BachelorFromCharmBaryonId, cascdata::BachelorId,
                  hf_cand_xic0_omegac0::V0DauPosX, hf_cand_xic0_omegac0::V0DauPosY, hf_cand_xic0_omegac0::V0DauPosZ, hf_cand_xic0_omegac0::V0DauPosXError, hf_cand_xic0_omegac0::V0DauPosYError, hf_cand_xic0_omegac0::V0DauPosZError, hf_cand_xic0_omegac0::V0DauPosPt,
                  hf_cand_xic0_omegac0::V0DauNegX, hf_cand_xic0_omegac0::V0DauNegY, hf_cand_xic0_omegac0::V0DauNegZ, hf_cand_xic0_omegac0::V0DauNegXError, hf_cand_xic0_omegac0::V0DauNegYError, hf_cand_xic0_omegac0::V0DauNegZError, hf_cand_xic0_omegac0::V0DauNegPt,
                  hf_cand_xic0_omegac0::V0VtxX, hf_cand_xic0_omegac0::V0VtxY, hf_cand_xic0_omegac0::V0VtxZ, hf_cand_xic0_omegac0::V0XError, hf_cand_xic0_omegac0::V0YError, hf_cand_xic0_omegac0::V0ZError, hf_cand_xic0_omegac0::V0Pt,
                  hf_cand_xic0_omegac0::XiBachelorX, hf_cand_xic0_omegac0::XiBachelorY, hf_cand_xic0_omegac0::XiBachelorZ, hf_cand_xic0_omegac0::XiBachelorXError, hf_cand_xic0_omegac0::XiBachelorYError, hf_cand_xic0_omegac0::XiBachelorZError, hf_cand_xic0_omegac0::XiBachelorPt,
                  hf_cand_xic0_omegac0::XiX, hf_cand_xic0_omegac0::XiY, hf_cand_xic0_omegac0::XiZ, hf_cand_xic0_omegac0::XiXError, hf_cand_xic0_omegac0::XiYError, hf_cand_xic0_omegac0::XiZError, hf_cand_xic0_omegac0::XiPt,
                  hf_cand_xic0_omegac0::Xic0BachelorX, hf_cand_xic0_omegac0::Xic0BachelorY, hf_cand_xic0_omegac0::Xic0BachelorZ, hf_cand_xic0_omegac0::Xic0BachelorXError, hf_cand_xic0_omegac0::Xic0BachelorYError, hf_cand_xic0_omegac0::Xic0BachelorZError, hf_cand_xic0_omegac0::Xic0BachelorPt,
                  hf_cand_xic0_omegac0::XDecayVtxCharmBaryon, hf_cand_xic0_omegac0::YDecayVtxCharmBaryon, hf_cand_xic0_omegac0::ZDecayVtxCharmBaryon, hf_cand_xic0_omegac0::Xic0XError, hf_cand_xic0_omegac0::Xic0YError, hf_cand_xic0_omegac0::Xic0ZError, hf_cand_xic0_omegac0::Xic0Pt,
                  hf_cand_casc::V0X, hf_cand_casc::V0Y, hf_cand_casc::V0Z, hf_cand_xic0_omegac0::XDecayVtxCascade, hf_cand_xic0_omegac0::YDecayVtxCascade, hf_cand_xic0_omegac0::ZDecayVtxCascade);

// table with results of reconstruction level MC matching
DECLARE_SOA_TABLE(HfXicToXiPiMCRec, "AOD", "HFXICXIPIMCREC", //!
                  hf_cand_xic0_omegac0::FlagMcMatchRec,
                  hf_cand_xic0_omegac0::DebugMcRec,
                  hf_cand_xic0_omegac0::OriginMcRec,
                  hf_cand_xic0_omegac0::CollisionMatched,
                  hf_cand::PtBhadMotherPart,
                  hf_cand::PdgBhadMotherPart,
                  o2::soa::Marker<1>);
DECLARE_SOA_TABLE(HfOmegacToXiPiMCRec, "AOD", "HFOMCXIPIMCREC", //!
                  hf_cand_xic0_omegac0::FlagMcMatchRec,
                  hf_cand_xic0_omegac0::DebugMcRec,
                  hf_cand_xic0_omegac0::OriginMcRec,
                  hf_cand_xic0_omegac0::CollisionMatched,
                  hf_cand::PtBhadMotherPart,
                  hf_cand::PdgBhadMotherPart,
                  o2::soa::Marker<2>);
DECLARE_SOA_TABLE(HfToOmegaPiMCRec, "AOD", "HFTOOMEPIMCREC", //!
                  hf_cand_xic0_omegac0::FlagMcMatchRec,
                  hf_cand_xic0_omegac0::DebugMcRec,
                  hf_cand_xic0_omegac0::OriginMcRec,
                  hf_cand_xic0_omegac0::CollisionMatched,
                  hf_cand::PtBhadMotherPart,
                  hf_cand::PdgBhadMotherPart,
                  o2::soa::Marker<3>);
DECLARE_SOA_TABLE(HfToOmegaKMCRec, "AOD", "HFTOOMEKMCREC", //!
                  hf_cand_xic0_omegac0::FlagMcMatchRec,
                  hf_cand_xic0_omegac0::DebugMcRec,
                  hf_cand_xic0_omegac0::OriginMcRec,
                  hf_cand_xic0_omegac0::CollisionMatched,
                  hf_cand::PtBhadMotherPart,
                  hf_cand::PdgBhadMotherPart,
                  o2::soa::Marker<4>);

// table with results of generator level MC matching
DECLARE_SOA_TABLE(HfXicToXiPiMCGen, "AOD", "HFXICXIPIMCGEN", //!
                  hf_cand_xic0_omegac0::FlagMcMatchGen, hf_cand_xic0_omegac0::DebugGenCharmBar, hf_cand_xic0_omegac0::DebugGenCasc, hf_cand_xic0_omegac0::DebugGenLambda,
                  hf_cand_xic0_omegac0::PtCharmBaryonGen, hf_cand_xic0_omegac0::RapidityCharmBaryonGen, hf_cand_xic0_omegac0::OriginMcGen, hf_cand::IdxBhadMotherPart, o2::soa::Marker<1>);
DECLARE_SOA_TABLE(HfOmegacToXiPiMCGen, "AOD", "HFOMECXIPIMCGEN", //!
                  hf_cand_xic0_omegac0::FlagMcMatchGen, hf_cand_xic0_omegac0::DebugGenCharmBar, hf_cand_xic0_omegac0::DebugGenCasc, hf_cand_xic0_omegac0::DebugGenLambda,
                  hf_cand_xic0_omegac0::PtCharmBaryonGen, hf_cand_xic0_omegac0::RapidityCharmBaryonGen, hf_cand_xic0_omegac0::OriginMcGen, hf_cand::IdxBhadMotherPart, o2::soa::Marker<2>);
DECLARE_SOA_TABLE(HfToOmegaPiMCGen, "AOD", "HFTOOMEPIMCGEN", //!
                  hf_cand_xic0_omegac0::FlagMcMatchGen, hf_cand_xic0_omegac0::DebugGenCharmBar, hf_cand_xic0_omegac0::DebugGenCasc, hf_cand_xic0_omegac0::DebugGenLambda,
                  hf_cand_xic0_omegac0::PtCharmBaryonGen, hf_cand_xic0_omegac0::RapidityCharmBaryonGen, hf_cand_xic0_omegac0::OriginMcGen, hf_cand::IdxBhadMotherPart, o2::soa::Marker<3>);
DECLARE_SOA_TABLE(HfToOmegaKMCGen, "AOD", "HFTOOMEKMCGEN", //!
                  hf_cand_xic0_omegac0::FlagMcMatchGen, hf_cand_xic0_omegac0::DebugGenCharmBar, hf_cand_xic0_omegac0::DebugGenCasc, hf_cand_xic0_omegac0::DebugGenLambda,
                  hf_cand_xic0_omegac0::PtCharmBaryonGen, hf_cand_xic0_omegac0::RapidityCharmBaryonGen, hf_cand_xic0_omegac0::OriginMcGen, hf_cand::IdxBhadMotherPart, o2::soa::Marker<4>);

// specific Xic to Xi Pi Pi candidate properties
namespace hf_cand_xic_to_xi_pi_pi
{
DECLARE_SOA_INDEX_COLUMN_FULL(Pi0, pi0, int, Tracks, "_pi0");
DECLARE_SOA_INDEX_COLUMN_FULL(Pi1, pi1, int, Tracks, "_pi1");
DECLARE_SOA_COLUMN(Sign, sign, int8_t);
DECLARE_SOA_COLUMN(InvMassXicPlus, invMassXicPlus, float);
DECLARE_SOA_COLUMN(InvMassXi, invMassXi, float);
DECLARE_SOA_COLUMN(InvMassLambda, invMassLambda, float);
DECLARE_SOA_COLUMN(InvMassXiPi0, invMassXiPi0, float);
DECLARE_SOA_COLUMN(InvMassXiPi1, invMassXiPi1, float);
DECLARE_SOA_COLUMN(XPvErr, xPvErr, float);
DECLARE_SOA_COLUMN(YPvErr, yPvErr, float);
DECLARE_SOA_COLUMN(ZPvErr, zPvErr, float);
DECLARE_SOA_COLUMN(XSvErr, xSvErr, float);
DECLARE_SOA_COLUMN(YSvErr, ySvErr, float);
DECLARE_SOA_COLUMN(ZSvErr, zSvErr, float);
DECLARE_SOA_COLUMN(CpaXi, cpaXi, float);
DECLARE_SOA_COLUMN(CpaXYXi, cpaXYXi, float);
DECLARE_SOA_COLUMN(CpaLambda, cpaLambda, float);
DECLARE_SOA_COLUMN(CpaXYLambda, cpaXYLambda, float);
DECLARE_SOA_COLUMN(CpaLambdaToXi, cpaLambdaToXi, float);
DECLARE_SOA_COLUMN(CpaXYLambdaToXi, cpaXYLambdaToXi, float);
DECLARE_SOA_COLUMN(PBachelorPi, pBachelorPi, float);
DECLARE_SOA_COLUMN(PPiFromLambda, pPiFromLambda, float);
DECLARE_SOA_COLUMN(PPrFromLambda, pPrFromLambda, float);
DECLARE_SOA_COLUMN(DcaXiDaughters, dcaXiDaughters, float);
DECLARE_SOA_COLUMN(DcaV0Daughters, dcaV0Daughters, float);
DECLARE_SOA_COLUMN(DcaPosToPV, dcaPosToPV, float);
DECLARE_SOA_COLUMN(DcaNegToPV, dcaNegToPV, float);
DECLARE_SOA_COLUMN(DcaBachelorToPV, dcaBachelorToPV, float);
DECLARE_SOA_COLUMN(DcaXYCascToPV, dcaXYCascToPV, float);
DECLARE_SOA_COLUMN(DcaZCascToPV, dcaZCascToPV, float);
// KF specific columns
DECLARE_SOA_COLUMN(Chi2TopoXicPlusToPV, chi2TopoXicPlusToPV, float);
DECLARE_SOA_COLUMN(Chi2TopoXicPlusToPVBefConst, chi2TopoXicPlusToPVBefConst, float);
DECLARE_SOA_COLUMN(Chi2PrimXi, chi2PrimXi, float);
DECLARE_SOA_COLUMN(Chi2PrimPi0, chi2PrimPi0, float);
DECLARE_SOA_COLUMN(Chi2PrimPi1, chi2PrimPi1, float);
DECLARE_SOA_COLUMN(Chi2DevPi0Pi1, chi2DevPi0Pi1, float);
DECLARE_SOA_COLUMN(Chi2DevPi0Xi, chi2DevPi0Xi, float);
DECLARE_SOA_COLUMN(Chi2DevPi1Xi, chi2DevPi1Xi, float);
DECLARE_SOA_COLUMN(DcaPi0Pi1, dcaPi0Pi1, float);
DECLARE_SOA_COLUMN(DcaPi0Xi, dcaPi0Xi, float);
DECLARE_SOA_COLUMN(DcaPi1Xi, dcaPi1Xi, float);
DECLARE_SOA_COLUMN(DcaXYPi0Pi1, dcaXYPi0Pi1, float);
DECLARE_SOA_COLUMN(DcaXYPi0Xi, dcaXYPi0Xi, float);
DECLARE_SOA_COLUMN(DcaXYPi1Xi, dcaXYPi1Xi, float);
DECLARE_SOA_COLUMN(KfDecayLength, kfDecayLength, float);
DECLARE_SOA_COLUMN(KfDecayLengthNormalised, kfDecayLengthNormalised, float);
DECLARE_SOA_COLUMN(KfDecayLengthXY, kfDecayLengthXY, float);
DECLARE_SOA_COLUMN(KfDecayLengthXYNormalised, kfDecayLengthXYNormalised, float);
// PID
DECLARE_SOA_COLUMN(NSigTpcPiFromXicPlus0, nSigTpcPiFromXicPlus0, float);
DECLARE_SOA_COLUMN(NSigTpcPiFromXicPlus1, nSigTpcPiFromXicPlus1, float);
DECLARE_SOA_COLUMN(NSigTpcBachelorPi, nSigTpcBachelorPi, float);
DECLARE_SOA_COLUMN(NSigTpcPiFromLambda, nSigTpcPiFromLambda, float);
DECLARE_SOA_COLUMN(NSigTpcPrFromLambda, nSigTpcPrFromLambda, float);
DECLARE_SOA_COLUMN(NSigTofPiFromXicPlus0, nSigTofPiFromXicPlus0, float);
DECLARE_SOA_COLUMN(NSigTofPiFromXicPlus1, nSigTofPiFromXicPlus1, float);
DECLARE_SOA_COLUMN(NSigTofBachelorPi, nSigTofBachelorPi, float);
DECLARE_SOA_COLUMN(NSigTofPiFromLambda, nSigTofPiFromLambda, float);
DECLARE_SOA_COLUMN(NSigTofPrFromLambda, nSigTofPrFromLambda, float);
// MC matching result:
DECLARE_SOA_COLUMN(FlagMcMatchRec, flagMcMatchRec, int8_t); // reconstruction level
DECLARE_SOA_COLUMN(FlagMcMatchGen, flagMcMatchGen, int8_t); // generator level
DECLARE_SOA_COLUMN(OriginMcRec, originMcRec, int8_t);
DECLARE_SOA_COLUMN(OriginMcGen, originMcGen, int8_t);
DECLARE_SOA_COLUMN(DecayLengthMcGen, decayLengthMcGen, float);
// Residuals and pulls
DECLARE_SOA_COLUMN(PtResidual, ptResidual, float);
DECLARE_SOA_COLUMN(PResidual, pResidual, float);
DECLARE_SOA_COLUMN(XPvResidual, xPvResidual, float);
DECLARE_SOA_COLUMN(YPvResidual, yPvResidual, float);
DECLARE_SOA_COLUMN(ZPvResidual, zPvResidual, float);
DECLARE_SOA_COLUMN(XPvPull, xPvPull, float);
DECLARE_SOA_COLUMN(YPvPull, yPvPull, float);
DECLARE_SOA_COLUMN(ZPvPull, zPvPull, float);
DECLARE_SOA_COLUMN(XSvResidual, xSvResidual, float);
DECLARE_SOA_COLUMN(YSvResidual, ySvResidual, float);
DECLARE_SOA_COLUMN(ZSvResidual, zSvResidual, float);
DECLARE_SOA_COLUMN(XSvPull, xSvPull, float);
DECLARE_SOA_COLUMN(YSvPull, ySvPull, float);
DECLARE_SOA_COLUMN(ZSvPull, zSvPull, float);
// Dynamic columns
DECLARE_SOA_DYNAMIC_COLUMN(PProng0, pProng0, //!
                           [](float px, float py, float pz) -> float { return RecoDecay::p(px, py, pz); });
DECLARE_SOA_DYNAMIC_COLUMN(PProng1, pProng1, //!
                           [](float px, float py, float pz) -> float { return RecoDecay::p(px, py, pz); });
DECLARE_SOA_DYNAMIC_COLUMN(PProng2, pProng2, //!
                           [](float px, float py, float pz) -> float { return RecoDecay::p(px, py, pz); });

} // end of namespace hf_cand_xic_to_xi_pi_pi

// declare dedicated Xic to Xi Pi Pi candidate table
DECLARE_SOA_TABLE(HfCandXicBase, "AOD", "HFCANDXICBASE",
                  hf_cand::CollisionId,
                  collision::PosX, collision::PosY, collision::PosZ,
                  hf_cand_xic_to_xi_pi_pi::XPvErr, hf_cand_xic_to_xi_pi_pi::YPvErr, hf_cand_xic_to_xi_pi_pi::ZPvErr,
                  // 3-prong specific columns
                  cascdata::CascadeId, hf_cand_xic_to_xi_pi_pi::Pi0Id, hf_cand_xic_to_xi_pi_pi::Pi1Id,
                  cascdata::BachelorId, cascdata::PosTrackId, cascdata::NegTrackId,
                  hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex,
                  hf_cand_xic_to_xi_pi_pi::XSvErr, hf_cand_xic_to_xi_pi_pi::YSvErr, hf_cand_xic_to_xi_pi_pi::ZSvErr,
                  hf_cand::ErrorDecayLength, hf_cand::ErrorDecayLengthXY,
                  hf_cand::Chi2PCA, hf_cand_xic_to_xi_pi_pi::InvMassXicPlus, hf_cand_xic_to_xi_pi_pi::Sign,
                  hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0,
                  hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1,
                  hf_cand::PxProng2, hf_cand::PyProng2, hf_cand::PzProng2,
                  hf_cand::ImpactParameter0, hf_cand::ImpactParameter1, hf_cand::ImpactParameter2,
                  hf_cand::ErrorImpactParameter0, hf_cand::ErrorImpactParameter1, hf_cand::ErrorImpactParameter2,
                  // cascade specific columns
                  hf_cand_xic_to_xi_pi_pi::PBachelorPi, hf_cand_xic_to_xi_pi_pi::PPiFromLambda, hf_cand_xic_to_xi_pi_pi::PPrFromLambda,
                  hf_cand_xic_to_xi_pi_pi::CpaXi, hf_cand_xic_to_xi_pi_pi::CpaXYXi, hf_cand_xic_to_xi_pi_pi::CpaLambda, hf_cand_xic_to_xi_pi_pi::CpaXYLambda, hf_cand_xic_to_xi_pi_pi::CpaLambdaToXi, hf_cand_xic_to_xi_pi_pi::CpaXYLambdaToXi,
                  hf_cand_xic_to_xi_pi_pi::InvMassXi, hf_cand_xic_to_xi_pi_pi::InvMassLambda, hf_cand_xic_to_xi_pi_pi::InvMassXiPi0, hf_cand_xic_to_xi_pi_pi::InvMassXiPi1,
                  // DCA
                  hf_cand_xic_to_xi_pi_pi::DcaXiDaughters, hf_cand_xic_to_xi_pi_pi::DcaV0Daughters, hf_cand_xic_to_xi_pi_pi::DcaPosToPV, hf_cand_xic_to_xi_pi_pi::DcaNegToPV, hf_cand_xic_to_xi_pi_pi::DcaBachelorToPV, hf_cand_xic_to_xi_pi_pi::DcaXYCascToPV, hf_cand_xic_to_xi_pi_pi::DcaZCascToPV,
                  // PID
                  hf_cand_xic_to_xi_pi_pi::NSigTpcPiFromXicPlus0, hf_cand_xic_to_xi_pi_pi::NSigTpcPiFromXicPlus1, hf_cand_xic_to_xi_pi_pi::NSigTpcBachelorPi, hf_cand_xic_to_xi_pi_pi::NSigTpcPiFromLambda, hf_cand_xic_to_xi_pi_pi::NSigTpcPrFromLambda,
                  hf_cand_xic_to_xi_pi_pi::NSigTofPiFromXicPlus0, hf_cand_xic_to_xi_pi_pi::NSigTofPiFromXicPlus1, hf_cand_xic_to_xi_pi_pi::NSigTofBachelorPi, hf_cand_xic_to_xi_pi_pi::NSigTofPiFromLambda, hf_cand_xic_to_xi_pi_pi::NSigTofPrFromLambda,
                  /* dynamic columns */
                  hf_cand::DecayLength<collision::PosX, collision::PosY, collision::PosZ, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex>,
                  hf_cand::DecayLengthXY<collision::PosX, collision::PosY, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex>,
                  hf_cand::DecayLengthNormalised<collision::PosX, collision::PosY, collision::PosZ, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex, hf_cand::ErrorDecayLength>,
                  hf_cand::DecayLengthXYNormalised<collision::PosX, collision::PosY, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ErrorDecayLengthXY>,
                  hf_cand::ImpactParameterNormalised0<hf_cand::ImpactParameter0, hf_cand::ErrorImpactParameter0>,
                  hf_cand::ImpactParameterNormalised1<hf_cand::ImpactParameter1, hf_cand::ErrorImpactParameter1>,
                  hf_cand::ImpactParameterNormalised2<hf_cand::ImpactParameter2, hf_cand::ErrorImpactParameter2>,
                  /* dynamic columns that use daughter momentum components */
                  hf_cand_xic_to_xi_pi_pi::PProng0<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0>,
                  hf_cand::PtProng0<hf_cand::PxProng0, hf_cand::PyProng0>,
                  hf_cand_xic_to_xi_pi_pi::PProng1<hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>,
                  hf_cand::PtProng1<hf_cand::PxProng1, hf_cand::PyProng1>,
                  hf_cand_xic_to_xi_pi_pi::PProng2<hf_cand::PxProng2, hf_cand::PyProng2, hf_cand::PzProng2>,
                  hf_cand::PtProng2<hf_cand::PxProng2, hf_cand::PyProng2>,
                  /* dynamic columns that use candidate momentum components */
                  hf_cand::Pt<hf_cand_3prong::Px, hf_cand_3prong::Py>,
                  hf_cand::P<hf_cand_3prong::Px, hf_cand_3prong::Py, hf_cand_3prong::Pz>,
                  hf_cand::PVector<hf_cand_3prong::Px, hf_cand_3prong::Py, hf_cand_3prong::Pz>,
                  hf_cand::Cpa<collision::PosX, collision::PosY, collision::PosZ, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex, hf_cand_3prong::Px, hf_cand_3prong::Py, hf_cand_3prong::Pz>,
                  hf_cand::CpaXY<collision::PosX, collision::PosY, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand_3prong::Px, hf_cand_3prong::Py>,
                  hf_cand::Ct<collision::PosX, collision::PosY, collision::PosZ, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex, hf_cand_3prong::Px, hf_cand_3prong::Py, hf_cand_3prong::Pz>,
                  hf_cand::ImpactParameterXY<collision::PosX, collision::PosY, collision::PosZ, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex, hf_cand_3prong::Px, hf_cand_3prong::Py, hf_cand_3prong::Pz>,
                  hf_cand_3prong::MaxNormalisedDeltaIP<collision::PosX, collision::PosY, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ErrorDecayLengthXY, hf_cand_3prong::Px, hf_cand_3prong::Py, hf_cand::ImpactParameter0, hf_cand::ErrorImpactParameter0, hf_cand::ImpactParameter1, hf_cand::ErrorImpactParameter1, hf_cand::ImpactParameter2, hf_cand::ErrorImpactParameter2, hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PxProng2, hf_cand::PyProng2>,
                  hf_cand::Eta<hf_cand_3prong::Px, hf_cand_3prong::Py, hf_cand_3prong::Pz>,
                  hf_cand::Phi<hf_cand_3prong::Px, hf_cand_3prong::Py>,
                  hf_cand::Y<hf_cand_3prong::Px, hf_cand_3prong::Py, hf_cand_3prong::Pz>);

// extended table with expression columns that can be used as arguments of dynamic columns
DECLARE_SOA_EXTENDED_TABLE_USER(HfCandXicExt, HfCandXicBase, "HFCANDXICEXT",
                                hf_cand_3prong::Px, hf_cand_3prong::Py, hf_cand_3prong::Pz);
using HfCandXic = HfCandXicExt;

// table with KF-specific variables
DECLARE_SOA_TABLE(HfCandXicKF, "AOD", "HFCANDXICKF",
                  hf_cand_xic_to_xi_pi_pi::KfDecayLength, hf_cand_xic_to_xi_pi_pi::KfDecayLengthNormalised, hf_cand_xic_to_xi_pi_pi::KfDecayLengthXY, hf_cand_xic_to_xi_pi_pi::KfDecayLengthXYNormalised,
                  cascdata::KFCascadeChi2, cascdata::KFV0Chi2,
                  hf_cand_xic_to_xi_pi_pi::Chi2TopoXicPlusToPVBefConst, hf_cand_xic_to_xi_pi_pi::Chi2TopoXicPlusToPV,
                  hf_cand_xic_to_xi_pi_pi::Chi2PrimXi, hf_cand_xic_to_xi_pi_pi::Chi2PrimPi0, hf_cand_xic_to_xi_pi_pi::Chi2PrimPi1,
                  hf_cand_xic_to_xi_pi_pi::Chi2DevPi0Pi1, hf_cand_xic_to_xi_pi_pi::Chi2DevPi0Xi, hf_cand_xic_to_xi_pi_pi::Chi2DevPi1Xi,
                  hf_cand_xic_to_xi_pi_pi::DcaPi0Pi1, hf_cand_xic_to_xi_pi_pi::DcaPi0Xi, hf_cand_xic_to_xi_pi_pi::DcaPi1Xi,
                  hf_cand_xic_to_xi_pi_pi::DcaXYPi0Pi1, hf_cand_xic_to_xi_pi_pi::DcaXYPi0Xi, hf_cand_xic_to_xi_pi_pi::DcaXYPi1Xi);

// table with results of reconstruction level MC matching
DECLARE_SOA_TABLE(HfCandXicMcRec, "AOD", "HFCANDXICMCREC",
                  hf_cand_xic_to_xi_pi_pi::FlagMcMatchRec,
                  hf_cand_xic_to_xi_pi_pi::OriginMcRec);

// table with results of generator level MC matching
DECLARE_SOA_TABLE(HfCandXicMcGen, "AOD", "HFCANDXICMCGEN",
                  hf_cand_xic_to_xi_pi_pi::FlagMcMatchGen,
                  hf_cand_xic_to_xi_pi_pi::OriginMcGen,
                  hf_cand::PdgBhadMotherPart,
                  hf_cand_xic_to_xi_pi_pi::DecayLengthMcGen);

// table with residuals and pulls of PV
DECLARE_SOA_TABLE(HfCandXicResid, "AOD", "HFCANDXICRESID",
                  hf_cand_xic_to_xi_pi_pi::OriginMcGen,
                  hf_cand_xic_to_xi_pi_pi::PResidual,
                  hf_cand_xic_to_xi_pi_pi::PtResidual,
                  hf_cand_xic_to_xi_pi_pi::XPvResidual,
                  hf_cand_xic_to_xi_pi_pi::YPvResidual,
                  hf_cand_xic_to_xi_pi_pi::ZPvResidual,
                  hf_cand_xic_to_xi_pi_pi::XPvPull,
                  hf_cand_xic_to_xi_pi_pi::YPvPull,
                  hf_cand_xic_to_xi_pi_pi::ZPvPull,
                  hf_cand_xic_to_xi_pi_pi::XSvResidual,
                  hf_cand_xic_to_xi_pi_pi::YSvResidual,
                  hf_cand_xic_to_xi_pi_pi::ZSvResidual,
                  hf_cand_xic_to_xi_pi_pi::XSvPull,
                  hf_cand_xic_to_xi_pi_pi::YSvPull,
                  hf_cand_xic_to_xi_pi_pi::ZSvPull);

// specific chic candidate properties
namespace hf_cand_chic
{
DECLARE_SOA_INDEX_COLUMN_FULL(Prong0, prong0, int, HfCand2Prong, "_0"); // Jpsi index
DECLARE_SOA_INDEX_COLUMN_FULL(Prong1, prong1, int, ECALs, "_1");
// MC matching result:
DECLARE_SOA_COLUMN(FlagMcMatchRec, flagMcMatchRec, int8_t);         // reconstruction level
DECLARE_SOA_COLUMN(FlagMcMatchGen, flagMcMatchGen, int8_t);         // generator level
DECLARE_SOA_COLUMN(OriginMcRec, originMcRec, int8_t);               // particle origin, reconstruction level
DECLARE_SOA_COLUMN(OriginMcGen, originMcGen, int8_t);               // particle origin, generator level
DECLARE_SOA_COLUMN(FlagMcDecayChanRec, flagMcDecayChanRec, int8_t); // resonant decay channel flag, reconstruction level
DECLARE_SOA_COLUMN(FlagMcDecayChanGen, flagMcDecayChanGen, int8_t); // resonant decay channel flag, generator level
DECLARE_SOA_COLUMN(JpsiToMuMuMass, jpsiToMuMuMass, float);          // Jpsi mass
} // namespace hf_cand_chic

// declare dedicated chi_c candidate table
DECLARE_SOA_TABLE(HfCandChicBase, "AOD", "HFCANDCHICBASE",
                  // general columns
                  HFCAND_COLUMNS,
                  // 2-prong specific columns
                  hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0,
                  hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1,
                  hf_cand::ImpactParameter0, hf_cand::ImpactParameter1,
                  hf_cand::ErrorImpactParameter0, hf_cand::ErrorImpactParameter1,
                  hf_cand_chic::Prong0Id, hf_cand_chic::Prong1Id,
                  hf_track_index::HFflag, hf_cand_chic::JpsiToMuMuMass,
                  /* dynamic columns */
                  hf_cand_2prong::M<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>,
                  hf_cand_2prong::M2<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>,
                  /* prong 2 */
                  //                  hf_cand::PtProng1<hf_cand::PxProng1, hf_cand::PyProng1>,
                  //                  hf_cand::Pt2Prong1<hf_cand::PxProng1, hf_cand::PyProng1>,
                  //                  hf_cand::PVectorProng1<hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>,
                  /* dynamic columns that use candidate momentum components */
                  hf_cand::Pt<hf_cand_2prong::Px, hf_cand_2prong::Py>,
                  hf_cand::Pt2<hf_cand_2prong::Px, hf_cand_2prong::Py>,
                  hf_cand::P<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::P2<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::PVector<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::Cpa<collision::PosX, collision::PosY, collision::PosZ, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex, hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::CpaXY<collision::PosX, collision::PosY, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand_2prong::Px, hf_cand_2prong::Py>,
                  hf_cand::Ct<collision::PosX, collision::PosY, collision::PosZ, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex, hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::ImpactParameterXY<collision::PosX, collision::PosY, collision::PosZ, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex, hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand_2prong::MaxNormalisedDeltaIP<collision::PosX, collision::PosY, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ErrorDecayLengthXY, hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand::ImpactParameter0, hf_cand::ErrorImpactParameter0, hf_cand::ImpactParameter1, hf_cand::ErrorImpactParameter1, hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PxProng1, hf_cand::PyProng1>,
                  hf_cand::Eta<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::Phi<hf_cand_2prong::Px, hf_cand_2prong::Py>,
                  hf_cand::Y<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::E<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::E2<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>);

// extended table with expression columns that can be used as arguments of dynamic columns
DECLARE_SOA_EXTENDED_TABLE_USER(HfCandChicExt, HfCandChicBase, "HFCANDCHICEXT",
                                hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz);

using HfCandChic = HfCandChicExt;

// table with results of reconstruction level MC matching
DECLARE_SOA_TABLE(HfCandChicMcRec, "AOD", "HFCANDCHICMCREC", //!
                  hf_cand_chic::FlagMcMatchRec,
                  hf_cand_chic::OriginMcRec,
                  hf_cand_chic::FlagMcDecayChanRec);

// table with results of generator level MC matching
DECLARE_SOA_TABLE(HfCandChicMcGen, "AOD", "HFCANDCHICMCGEN", //!
                  hf_cand_chic::FlagMcMatchGen,
                  hf_cand_chic::OriginMcGen,
                  hf_cand_chic::FlagMcDecayChanGen);

// specific Lb candidate properties
namespace hf_cand_lb
{
DECLARE_SOA_INDEX_COLUMN_FULL(Prong0, prong0, int, HfCand3Prong, "_0"); // Lb index
// MC matching result:
DECLARE_SOA_COLUMN(FlagMcMatchRec, flagMcMatchRec, int8_t);         // main decay channel, reconstruction level
DECLARE_SOA_COLUMN(FlagMcMatchGen, flagMcMatchGen, int8_t);         // main decay channel, generator level
DECLARE_SOA_COLUMN(FlagMcDecayChanRec, flagMcDecayChanRec, int8_t); // resonant decay channel, reconstruction level
DECLARE_SOA_COLUMN(FlagMcDecayChanGen, flagMcDecayChanGen, int8_t); // resonant decay channel, generator level
DECLARE_SOA_COLUMN(OriginMcRec, originMcRec, int8_t);               // particle origin, reconstruction level
DECLARE_SOA_COLUMN(OriginMcGen, originMcGen, int8_t);               // particle origin, generator level
DECLARE_SOA_COLUMN(FlagWrongCollision, flagWrongCollision, int8_t); // reconstruction level
DECLARE_SOA_COLUMN(DebugMcRec, debugMcRec, int8_t);                 // debug flag for mis-association reconstruction level

enum DecayTypeMc : uint8_t { LbToLcPiToPKPiPi = 0,
                             LbToLcKToPKPiK,
                             B0ToDplusPiToPiKPiPi,
                             PartlyRecoDecay,
                             OtherDecay,
                             NDecayTypeMc };
} // namespace hf_cand_lb

// declare dedicated Lb candidate table
DECLARE_SOA_TABLE(HfCandLbBase, "AOD", "HFCANDLBBASE",
                  // general columns
                  HFCAND_COLUMNS,
                  // 3-prong specific columns
                  hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0,
                  hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1,
                  hf_cand::ImpactParameter0, hf_cand::ImpactParameter1,
                  hf_cand::ErrorImpactParameter0, hf_cand::ErrorImpactParameter1,
                  /* dynamic columns */
                  hf_cand_2prong::M<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>,
                  hf_cand_2prong::M2<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>,
                  hf_cand_2prong::ImpactParameterProduct<hf_cand::ImpactParameter0, hf_cand::ImpactParameter1>,
                  /* dynamic columns that use candidate momentum components */
                  hf_cand::Pt<hf_cand_2prong::Px, hf_cand_2prong::Py>,
                  hf_cand::Pt2<hf_cand_2prong::Px, hf_cand_2prong::Py>,
                  hf_cand::P<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::P2<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::PVector<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::Cpa<collision::PosX, collision::PosY, collision::PosZ, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex, hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::CpaXY<collision::PosX, collision::PosY, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand_2prong::Px, hf_cand_2prong::Py>,
                  hf_cand::Ct<collision::PosX, collision::PosY, collision::PosZ, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex, hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::ImpactParameterXY<collision::PosX, collision::PosY, collision::PosZ, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex, hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand_2prong::MaxNormalisedDeltaIP<collision::PosX, collision::PosY, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ErrorDecayLengthXY, hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand::ImpactParameter0, hf_cand::ErrorImpactParameter0, hf_cand::ImpactParameter1, hf_cand::ErrorImpactParameter1, hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PxProng1, hf_cand::PyProng1>,
                  hf_cand::Eta<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::Phi<hf_cand_2prong::Px, hf_cand_2prong::Py>,
                  hf_cand::Y<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::E<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::E2<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>);

// extended table with expression columns that can be used as arguments of dynamic columns
DECLARE_SOA_EXTENDED_TABLE_USER(HfCandLbExt, HfCandLbBase, "HFCANDLBEXT",
                                hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz);

DECLARE_SOA_TABLE(HfCandLbProngs, "AOD", "HFCANDLBPRONGS",
                  hf_cand_lb::Prong0Id, hf_track_index::Prong1Id);

using HfCandLb = soa::Join<HfCandLbExt, HfCandLbProngs>;

// table with results of reconstruction level MC matching
DECLARE_SOA_TABLE(HfCandLbMcRec, "AOD", "HFCANDLBMCREC", //!
                  hf_cand_lb::FlagMcMatchRec,
                  hf_cand_lb::FlagMcDecayChanRec,
                  hf_cand_lb::OriginMcRec,
                  hf_cand_lb::DebugMcRec);

// table with results of generator level MC matching
DECLARE_SOA_TABLE(HfCandLbMcGen, "AOD", "HFCANDLBMCGEN", //!
                  hf_cand_lb::FlagMcMatchGen,
                  hf_cand_lb::FlagMcDecayChanGen,
                  hf_cand_lb::OriginMcGen);

// specific B0 candidate properties
namespace hf_cand_b0
{
DECLARE_SOA_INDEX_COLUMN_FULL(Prong0, prong0, int, HfCand3Prong, "_0"); // D index
// MC matching result:
DECLARE_SOA_COLUMN(FlagMcMatchRec, flagMcMatchRec, int8_t);         // main decay channel, reconstruction level
DECLARE_SOA_COLUMN(FlagMcMatchGen, flagMcMatchGen, int8_t);         // main decay channel, generator level
DECLARE_SOA_COLUMN(FlagMcDecayChanRec, flagMcDecayChanRec, int8_t); // resonant decay channel, reconstruction level
DECLARE_SOA_COLUMN(FlagMcDecayChanGen, flagMcDecayChanGen, int8_t); // resonant decay channel, generator level
DECLARE_SOA_COLUMN(OriginMcRec, originMcRec, int8_t);               // particle origin, reconstruction level
DECLARE_SOA_COLUMN(OriginMcGen, originMcGen, int8_t);               // particle origin, generator level
DECLARE_SOA_COLUMN(FlagWrongCollision, flagWrongCollision, int8_t); // reconstruction level
DECLARE_SOA_COLUMN(DebugMcRec, debugMcRec, int8_t);                 // debug flag for mis-association reconstruction level

enum DecayTypeMc : uint8_t { B0ToDplusPiToPiKPiPi = 0,
                             B0ToDsPiToKKPiPi,
                             BsToDsPiToKKPiPi,
                             B0ToDplusKToPiKPiK,
                             B0ToDstarPiToD0PiPiToKPiPiPi,
                             PartlyRecoDecay,
                             OtherDecay,
                             NDecayTypeMc };
} // namespace hf_cand_b0

// declare dedicated B0 decay candidate table
DECLARE_SOA_TABLE(HfCandB0Base, "AOD", "HFCANDB0BASE",
                  // general columns
                  HFCAND_COLUMNS,
                  // 2-prong specific columns
                  hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0,
                  hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1,
                  hf_cand::ImpactParameter0, hf_cand::ImpactParameter1,
                  hf_cand::ErrorImpactParameter0, hf_cand::ErrorImpactParameter1,
                  /* dynamic columns */
                  hf_cand_2prong::M<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>,
                  hf_cand_2prong::M2<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>,
                  hf_cand_2prong::ImpactParameterProduct<hf_cand::ImpactParameter0, hf_cand::ImpactParameter1>,
                  hf_cand_2prong::CosThetaStar<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>,
                  hf_cand_2prong::ImpactParameterProngSqSum<hf_cand::ImpactParameter0, hf_cand::ImpactParameter1>,
                  /* dynamic columns that use candidate momentum components */
                  hf_cand::Pt<hf_cand_2prong::Px, hf_cand_2prong::Py>,
                  hf_cand::Pt2<hf_cand_2prong::Px, hf_cand_2prong::Py>,
                  hf_cand::P<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::P2<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::PVector<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::Cpa<collision::PosX, collision::PosY, collision::PosZ, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex, hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::CpaXY<collision::PosX, collision::PosY, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand_2prong::Px, hf_cand_2prong::Py>,
                  hf_cand::Ct<collision::PosX, collision::PosY, collision::PosZ, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex, hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::ImpactParameterXY<collision::PosX, collision::PosY, collision::PosZ, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex, hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand_2prong::MaxNormalisedDeltaIP<collision::PosX, collision::PosY, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ErrorDecayLengthXY, hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand::ImpactParameter0, hf_cand::ErrorImpactParameter0, hf_cand::ImpactParameter1, hf_cand::ErrorImpactParameter1, hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PxProng1, hf_cand::PyProng1>,
                  hf_cand::Eta<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::Phi<hf_cand_2prong::Px, hf_cand_2prong::Py>,
                  hf_cand::Y<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::E<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::E2<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>);

DECLARE_SOA_TABLE(HfCandB0DStar, "AOD", "HFCANDB0DSTAR",
                  // general columns
                  HFCAND_COLUMNS,
                  /* prong 2 */
                  hf_cand::ImpactParameterNormalised2<hf_cand::ImpactParameter2, hf_cand::ErrorImpactParameter2>,
                  hf_cand::PtProng2<hf_cand::PxProng2, hf_cand::PyProng2>,
                  hf_cand::Pt2Prong2<hf_cand::PxProng2, hf_cand::PyProng2>,
                  hf_cand::PVectorProng2<hf_cand::PxProng2, hf_cand::PyProng2, hf_cand::PzProng2>,
                  // 3-prong specific columns
                  hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0,
                  hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1,
                  hf_cand::PxProng2, hf_cand::PyProng2, hf_cand::PzProng2,
                  hf_cand::ImpactParameter0, hf_cand::ImpactParameter1, hf_cand::ImpactParameter2,
                  hf_cand::ErrorImpactParameter0, hf_cand::ErrorImpactParameter1, hf_cand::ErrorImpactParameter2,
                  /* dynamic columns */
                  hf_cand_3prong::M<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1, hf_cand::PxProng2, hf_cand::PyProng2, hf_cand::PzProng2>,
                  hf_cand_3prong::M2<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1, hf_cand::PxProng2, hf_cand::PyProng2, hf_cand::PzProng2>,
                  hf_cand_2prong::ImpactParameterProduct<hf_cand::ImpactParameter0, hf_cand::ImpactParameter2>,
                  hf_cand_3prong::ImpactParameterProngSqSum<hf_cand::ImpactParameter0, hf_cand::ImpactParameter1, hf_cand::ImpactParameter2>,
                  /* dynamic columns that use candidate momentum components */
                  hf_cand::Pt<hf_cand_3prong::Px, hf_cand_3prong::Py>,
                  hf_cand::Pt2<hf_cand_3prong::Px, hf_cand_3prong::Py>,
                  hf_cand::P<hf_cand_3prong::Px, hf_cand_3prong::Py, hf_cand_3prong::Pz>,
                  hf_cand::P2<hf_cand_3prong::Px, hf_cand_3prong::Py, hf_cand_3prong::Pz>,
                  hf_cand::PVector<hf_cand_3prong::Px, hf_cand_3prong::Py, hf_cand_3prong::Pz>,
                  hf_cand::Cpa<collision::PosX, collision::PosY, collision::PosZ, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex, hf_cand_3prong::Px, hf_cand_3prong::Py, hf_cand_3prong::Pz>,
                  hf_cand::CpaXY<collision::PosX, collision::PosY, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand_3prong::Px, hf_cand_3prong::Py>,
                  hf_cand::Ct<collision::PosX, collision::PosY, collision::PosZ, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex, hf_cand_3prong::Px, hf_cand_3prong::Py, hf_cand_3prong::Pz>,
                  hf_cand::ImpactParameterXY<collision::PosX, collision::PosY, collision::PosZ, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex, hf_cand_3prong::Px, hf_cand_3prong::Py, hf_cand_3prong::Pz>,
                  hf_cand_3prong::MaxNormalisedDeltaIP<collision::PosX, collision::PosY, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ErrorDecayLengthXY, hf_cand_3prong::Px, hf_cand_3prong::Py, hf_cand::ImpactParameter0, hf_cand::ErrorImpactParameter0, hf_cand::ImpactParameter1, hf_cand::ErrorImpactParameter1, hf_cand::ImpactParameter2, hf_cand::ErrorImpactParameter2, hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PxProng2, hf_cand::PyProng2>,
                  hf_cand::Eta<hf_cand_3prong::Px, hf_cand_3prong::Py, hf_cand_3prong::Pz>,
                  hf_cand::Phi<hf_cand_3prong::Px, hf_cand_3prong::Py>,
                  hf_cand::Y<hf_cand_3prong::Px, hf_cand_3prong::Py, hf_cand_3prong::Pz>,
                  hf_cand::E<hf_cand_3prong::Px, hf_cand_3prong::Py, hf_cand_3prong::Pz>,
                  hf_cand::E2<hf_cand_3prong::Px, hf_cand_3prong::Py, hf_cand_3prong::Pz>);

// extended table with expression columns that can be used as arguments of dynamic columns
DECLARE_SOA_EXTENDED_TABLE_USER(HfCandB0Ext, HfCandB0Base, "HFCANDB0EXT",
                                hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz);

DECLARE_SOA_TABLE(HfCandB0Prongs, "AOD", "HFCANDB0PRONGS",
                  hf_cand_b0::Prong0Id, hf_track_index::Prong1Id);

using HfCandB0 = soa::Join<HfCandB0Ext, HfCandB0Prongs>;

// extended table with expression columns that can be used as arguments of dynamic columns
DECLARE_SOA_EXTENDED_TABLE_USER(HfCandB0DStExt, HfCandB0DStar, "HFCANDB0DSTEXT",
                                hf_cand_3prong::Px, hf_cand_3prong::Py, hf_cand_3prong::Pz);

// table with results of reconstruction level MC matching
DECLARE_SOA_TABLE(HfCandB0McRec, "AOD", "HFCANDB0MCREC",
                  hf_cand_b0::FlagMcMatchRec,
                  hf_cand_b0::FlagMcDecayChanRec,
                  hf_cand_b0::OriginMcRec,
                  hf_cand_b0::DebugMcRec);

// table with results of generator level MC matching
DECLARE_SOA_TABLE(HfCandB0McGen, "AOD", "HFCANDB0MCGEN",
                  hf_cand_b0::FlagMcMatchGen,
                  hf_cand_b0::FlagMcDecayChanGen,
                  hf_cand_b0::OriginMcGen);

// specific Bs candidate properties
namespace hf_cand_bs
{
DECLARE_SOA_INDEX_COLUMN_FULL(Prong0, prong0, int, HfCand3Prong, "_0"); // Ds index
// MC matching result:
DECLARE_SOA_COLUMN(FlagMcMatchRec, flagMcMatchRec, int8_t);                // main decay channel, reconstruction level
DECLARE_SOA_COLUMN(FlagMcMatchGen, flagMcMatchGen, int8_t);                // main decay channel, generator level
DECLARE_SOA_COLUMN(FlagMcDecayChanRec, flagMcDecayChanRec, int8_t);        // resonant decay channel, reconstruction level
DECLARE_SOA_COLUMN(FlagMcDecayChanGen, flagMcDecayChanGen, int8_t);        // resonant decay channel, generator level
DECLARE_SOA_COLUMN(OriginMcRec, originMcRec, int8_t);                      // particle origin, reconstruction level
DECLARE_SOA_COLUMN(OriginMcGen, originMcGen, int8_t);                      // particle origin, generator level
DECLARE_SOA_COLUMN(FlagWrongCollision, flagWrongCollision, int8_t);        // reconstruction level
DECLARE_SOA_COLUMN(DebugMcRec, debugMcRec, int8_t);                        // debug flag for mis-association reconstruction level
DECLARE_SOA_DYNAMIC_COLUMN(ImpactParameterProduct, impactParameterProduct, // Impact parameter product for Bs -> J/Psi phi
                           [](float pxJpsiDauPos, float pyJpsiDauPos, float pzJpsiDauPos, float pxJpsiDauNeg, float pyJpsiDauNeg, float pzJpsiDauNeg, float pxLfTrack0, float pyLfTrack0, float pzLfTrack0, float pxLfTrack1, float pyLfTrack1, float pzLfTrack1, float xVtxP, float yVtxP, float zVtxP, float xVtxS, float yVtxS, float zVtxS) -> float {
                             float impParJpsi = RecoDecay::impParXY(std::array{xVtxP, yVtxP, zVtxP}, std::array{xVtxS, yVtxS, zVtxS}, RecoDecay::pVec(std::array{pxJpsiDauPos, pyJpsiDauPos, pzJpsiDauPos}, std::array{pxJpsiDauNeg, pyJpsiDauNeg, pzJpsiDauNeg}));
                             float impParPhi = RecoDecay::impParXY(std::array{xVtxP, yVtxP, zVtxP}, std::array{xVtxS, yVtxS, zVtxS}, RecoDecay::pVec(std::array{pxLfTrack0, pyLfTrack0, pzLfTrack0}, std::array{pxLfTrack1, pyLfTrack1, pzLfTrack1}));
                             return impParJpsi * impParPhi;
                           });
DECLARE_SOA_DYNAMIC_COLUMN(ImpactParameterProductJpsi, impactParameterProductJpsi, // J/Psi impact parameter for Bs -> J/Psi phi
                           [](float dcaDauPos, float dcaDauNeg) -> float { return dcaDauPos * dcaDauNeg; });
DECLARE_SOA_DYNAMIC_COLUMN(ImpactParameterProductPhi, impactParameterProductPhi, // J/Psi impact parameter for Bs -> J/Psi phi
                           [](float dcaLfTrack0, float dcaLfTrack1) -> float { return dcaLfTrack0 * dcaLfTrack1; });

enum DecayTypeMc : uint8_t { BsToDsPiToPhiPiPiToKKPiPi = 0, // Bs(bar) → Ds∓ π± → (Phi π∓) π± → (K- K+ π∓) π±
                             BsToDsPiToK0starKPiToKKPiPi,   // Bs(bar) → Ds∓ π± → (K0* K∓) π± → (K- K+ π∓) π±
                             B0ToDsPiToPhiPiPiToKKPiPi,     // B0(bar) → Ds± π∓ → (Phi π±) π∓ → (K- K+ π±) π∓
                             B0ToDsPiToK0starKPiToKKPiPi,   // B0(bar) → Ds± π∓ → (K0* K±) π∓ → (K- K+ π±) π∓
                             BsToDsKToPhiPiKToKKPiK,        // Bs(bar) → Ds± K∓ → (Phi π∓) K∓ → (K- K+ π±) K∓
                             BsToDsKToK0starKKToKKPiK,      // Bs(bar) → Ds± K∓ → (K0* K±) K∓ → (K- K+ π±) K∓
                             PartlyRecoDecay,               // 4 final state particles have another common b-hadron ancestor
                             OtherDecay,
                             NDecayTypeMc }; // counter of differentiated MC decay types

enum class DecayTypeBToJpsiMc : uint8_t { BsToJpsiPhiToMuMuKK = 0, // Bs(bar) → J/Psi Phi → (µ+ µ-) (K- K+)
                                          PartlyRecoDecay,         // 4 final state particles have another common b-hadron ancestor
                                          OtherDecay,
                                          NDecayTypeMc }; // counter of differentiated MC decay types

} // namespace hf_cand_bs

// declare dedicated Bs decay candidate table
DECLARE_SOA_TABLE(HfCandBsBase, "AOD", "HFCANDBSBASE",
                  // general columns
                  HFCAND_COLUMNS,
                  // 2-prong specific columns
                  hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0,
                  hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1,
                  hf_cand::ImpactParameter0, hf_cand::ImpactParameter1,
                  hf_cand::ErrorImpactParameter0, hf_cand::ErrorImpactParameter1,
                  /* dynamic columns */
                  hf_cand_2prong::M<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>,
                  hf_cand_2prong::M2<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>,
                  hf_cand_2prong::ImpactParameterProduct<hf_cand::ImpactParameter0, hf_cand::ImpactParameter1>,
                  hf_cand_2prong::CosThetaStar<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>,
                  hf_cand_2prong::ImpactParameterProngSqSum<hf_cand::ImpactParameter0, hf_cand::ImpactParameter1>,
                  /* dynamic columns that use candidate momentum components */
                  hf_cand::Pt<hf_cand_2prong::Px, hf_cand_2prong::Py>,
                  hf_cand::Pt2<hf_cand_2prong::Px, hf_cand_2prong::Py>,
                  hf_cand::P<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::P2<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::PVector<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::Cpa<collision::PosX, collision::PosY, collision::PosZ, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex, hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::CpaXY<collision::PosX, collision::PosY, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand_2prong::Px, hf_cand_2prong::Py>,
                  hf_cand::Ct<collision::PosX, collision::PosY, collision::PosZ, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex, hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::ImpactParameterXY<collision::PosX, collision::PosY, collision::PosZ, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex, hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand_2prong::MaxNormalisedDeltaIP<collision::PosX, collision::PosY, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ErrorDecayLengthXY, hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand::ImpactParameter0, hf_cand::ErrorImpactParameter0, hf_cand::ImpactParameter1, hf_cand::ErrorImpactParameter1, hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PxProng1, hf_cand::PyProng1>,
                  hf_cand::Eta<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::Phi<hf_cand_2prong::Px, hf_cand_2prong::Py>,
                  hf_cand::Y<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::E<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::E2<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand_2prong::CtXY<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1, collision::PosX, collision::PosY, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex>,
                  o2::soa::Marker<1>);

// extended table with expression columns that can be used as arguments of dynamic columns
DECLARE_SOA_EXTENDED_TABLE_USER(HfCandBsExt, HfCandBsBase, "HFCANDBSEXT",
                                hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz);

DECLARE_SOA_TABLE(HfCandBsProngs, "AOD", "HFCANDBSPRONGS",
                  hf_cand_bs::Prong0Id, hf_track_index::Prong1Id);

using HfCandBs = soa::Join<HfCandBsExt, HfCandBsProngs>;

// table with results of reconstruction level MC matching
DECLARE_SOA_TABLE(HfCandBsMcRec, "AOD", "HFCANDBSMCREC",
                  hf_cand_bs::FlagMcMatchRec,
                  hf_cand_bs::FlagMcDecayChanRec);

// table with results of generator level MC matching
DECLARE_SOA_TABLE(HfCandBsMcGen, "AOD", "HFCANDBSMCGEN",
                  hf_cand_bs::FlagMcMatchGen,
                  hf_cand_bs::FlagMcDecayChanGen);

namespace hf_cand_4prong
{
DECLARE_SOA_EXPRESSION_COLUMN(Px, px, //!
                              float, 1.f * aod::hf_cand::pxProng0 + 1.f * aod::hf_cand::pxProng1 + 1.f * aod::hf_cand::pxProng2 + 1.f * aod::hf_cand::pxProng3);
DECLARE_SOA_EXPRESSION_COLUMN(Py, py, //!
                              float, 1.f * aod::hf_cand::pyProng0 + 1.f * aod::hf_cand::pyProng1 + 1.f * aod::hf_cand::pyProng2 + 1.f * aod::hf_cand::pyProng3);
DECLARE_SOA_EXPRESSION_COLUMN(Pz, pz, //!
                              float, 1.f * aod::hf_cand::pzProng0 + 1.f * aod::hf_cand::pzProng1 + 1.f * aod::hf_cand::pzProng2 + 1.f * aod::hf_cand::pzProng3);
DECLARE_SOA_DYNAMIC_COLUMN(M, m, //!
                           [](float px0, float py0, float pz0, float px1, float py1, float pz1, float px2, float py2, float pz2, float px3, float py3, float pz3, const std::array<double, 4>& m) -> float { return RecoDecay::m(std::array{std::array{px0, py0, pz0}, std::array{px1, py1, pz1}, std::array{px2, py2, pz2}, std::array{px3, py3, pz3}}, m); });
DECLARE_SOA_DYNAMIC_COLUMN(M2, m2, //!
                           [](float px0, float py0, float pz0, float px1, float py1, float pz1, float px2, float py2, float pz2, float px3, float py3, float pz3, const std::array<double, 4>& m) -> float { return RecoDecay::m2(std::array{std::array{px0, py0, pz0}, std::array{px1, py1, pz1}, std::array{px2, py2, pz2}, std::array{px3, py3, pz3}}, m); });
DECLARE_SOA_DYNAMIC_COLUMN(ImpactParameterProngSqSum, impactParameterProngSqSum, //!
                           [](float impParProng0, float impParProng1, float impParProng2, float impParProng3) -> float { return RecoDecay::sumOfSquares(impParProng0, impParProng1, impParProng2, impParProng3); });
DECLARE_SOA_DYNAMIC_COLUMN(MaxNormalisedDeltaIP, maxNormalisedDeltaIP, //!
                           [](float xVtxP, float yVtxP, float xVtxS, float yVtxS, float errDlxy, float pxM, float pyM, float ip0, float errIp0, float ip1, float errIp1, float ip2, float errIp2, float ip3, float errIp3, float px0, float py0, float px1, float py1, float px2, float py2, float px3, float py3) -> float { return RecoDecay::maxNormalisedDeltaIP(std::array{xVtxP, yVtxP}, std::array{xVtxS, yVtxS}, errDlxy, std::array{pxM, pyM}, std::array{ip0, ip1, ip2, ip3}, std::array{errIp0, errIp1, errIp2, errIp3}, std::array{std::array{px0, py0}, std::array{px1, py1}, std::array{px2, py2}, std::array{px3, py3}}); });
DECLARE_SOA_DYNAMIC_COLUMN(CtXY, ctXY, //!
                           [](float px0, float py0, float pz0, float px1, float py1, float pz1, float px2, float py2, float pz2, float px3, float py3, float pz3, float xVtxP, float yVtxP, float xVtxS, float yVtxS, const std::array<double, 4>& m) -> float { return RecoDecay::ctXY(std::array{xVtxP, yVtxP}, std::array{xVtxS, yVtxS}, std::array{std::array{px0, py0, pz0}, std::array{px1, py1, pz1}, std::array{px2, py2, pz2}, std::array{px3, py3, pz3}}, m); });
} // namespace hf_cand_4prong

// declare dedicated Bs -> J/Psi phi decay candidate table
// convention: prongs 0 and 1 should be J/Psi decay products, 2 and 3 should be phi decay products
DECLARE_SOA_TABLE(HfCandBsJPBase, "AOD", "HFCANDBSJPBASE",
                  // general columns
                  HFCAND_COLUMNS,
                  /* prong 2 */ hf_cand::ImpactParameterNormalised2<hf_cand::ImpactParameter2, hf_cand::ErrorImpactParameter2>,
                  hf_cand::PtProng2<hf_cand::PxProng2, hf_cand::PyProng2>,
                  hf_cand::Pt2Prong2<hf_cand::PxProng2, hf_cand::PyProng2>,
                  hf_cand::PVectorProng2<hf_cand::PxProng2, hf_cand::PyProng2, hf_cand::PzProng2>,
                  /* prong 3 */ hf_cand::ImpactParameterNormalised3<hf_cand::ImpactParameter3, hf_cand::ErrorImpactParameter3>,
                  hf_cand::PtProng3<hf_cand::PxProng3, hf_cand::PyProng3>,
                  hf_cand::Pt2Prong3<hf_cand::PxProng3, hf_cand::PyProng3>,
                  hf_cand::PVectorProng3<hf_cand::PxProng3, hf_cand::PyProng3, hf_cand::PzProng3>,
                  // 4-prong specific columns
                  o2::soa::Index<>,
                  hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0,
                  hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1,
                  hf_cand::PxProng2, hf_cand::PyProng2, hf_cand::PzProng2,
                  hf_cand::PxProng3, hf_cand::PyProng3, hf_cand::PzProng3,
                  hf_cand::ImpactParameter0, hf_cand::ImpactParameter1, hf_cand::ImpactParameter2, hf_cand::ImpactParameter3,
                  hf_cand::ErrorImpactParameter0, hf_cand::ErrorImpactParameter1, hf_cand::ErrorImpactParameter2, hf_cand::ErrorImpactParameter3,
                  /* dynamic columns */
                  hf_cand_4prong::M<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1, hf_cand::PxProng2, hf_cand::PyProng2, hf_cand::PzProng2, hf_cand::PxProng3, hf_cand::PyProng3, hf_cand::PzProng3>,
                  hf_cand_4prong::M2<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1, hf_cand::PxProng2, hf_cand::PyProng2, hf_cand::PzProng2, hf_cand::PxProng3, hf_cand::PyProng3, hf_cand::PzProng3>,
                  hf_cand_4prong::ImpactParameterProngSqSum<hf_cand::ImpactParameter0, hf_cand::ImpactParameter1, hf_cand::ImpactParameter2, hf_cand::ImpactParameter3>,
                  hf_cand_bs::ImpactParameterProduct<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1, hf_cand::PxProng2, hf_cand::PyProng2, hf_cand::PzProng2, hf_cand::PxProng3, hf_cand::PyProng3, hf_cand::PzProng3, collision::PosX, collision::PosY, collision::PosZ, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex>,
                  hf_cand_bs::ImpactParameterProductJpsi<hf_cand::ImpactParameter0, hf_cand::ImpactParameter1>,
                  hf_cand_bs::ImpactParameterProductPhi<hf_cand::ImpactParameter2, hf_cand::ImpactParameter3>,
                  /* dynamic columns that use candidate momentum components */
                  hf_cand::Pt<hf_cand_4prong::Px, hf_cand_4prong::Py>,
                  hf_cand::Pt2<hf_cand_4prong::Px, hf_cand_4prong::Py>,
                  hf_cand::P<hf_cand_4prong::Px, hf_cand_4prong::Py, hf_cand_4prong::Pz>,
                  hf_cand::P2<hf_cand_4prong::Px, hf_cand_4prong::Py, hf_cand_4prong::Pz>,
                  hf_cand::PVector<hf_cand_4prong::Px, hf_cand_4prong::Py, hf_cand_4prong::Pz>,
                  hf_cand::Cpa<collision::PosX, collision::PosY, collision::PosZ, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex, hf_cand_4prong::Px, hf_cand_4prong::Py, hf_cand_4prong::Pz>,
                  hf_cand::CpaXY<collision::PosX, collision::PosY, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand_4prong::Px, hf_cand_4prong::Py>,
                  hf_cand::Ct<collision::PosX, collision::PosY, collision::PosZ, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex, hf_cand_4prong::Px, hf_cand_4prong::Py, hf_cand_4prong::Pz>,
                  hf_cand::ImpactParameterXY<collision::PosX, collision::PosY, collision::PosZ, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex, hf_cand_4prong::Px, hf_cand_4prong::Py, hf_cand_4prong::Pz>,
                  hf_cand_4prong::MaxNormalisedDeltaIP<collision::PosX, collision::PosY, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ErrorDecayLengthXY, hf_cand_4prong::Px, hf_cand_4prong::Py, hf_cand::ImpactParameter0, hf_cand::ErrorImpactParameter0, hf_cand::ImpactParameter1, hf_cand::ErrorImpactParameter1, hf_cand::ImpactParameter2, hf_cand::ErrorImpactParameter2, hf_cand::ImpactParameter3, hf_cand::ErrorImpactParameter3, hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PxProng2, hf_cand::PyProng2, hf_cand::PxProng3, hf_cand::PyProng3>,
                  hf_cand::Eta<hf_cand_4prong::Px, hf_cand_4prong::Py, hf_cand_4prong::Pz>,
                  hf_cand::Phi<hf_cand_4prong::Px, hf_cand_4prong::Py>,
                  hf_cand::Y<hf_cand_4prong::Px, hf_cand_4prong::Py, hf_cand_4prong::Pz>,
                  hf_cand::E<hf_cand_4prong::Px, hf_cand_4prong::Py, hf_cand_4prong::Pz>,
                  hf_cand::E2<hf_cand_4prong::Px, hf_cand_4prong::Py, hf_cand_4prong::Pz>,
                  hf_cand_4prong::CtXY<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1, hf_cand::PxProng2, hf_cand::PyProng2, hf_cand::PzProng2, hf_cand::PxProng3, hf_cand::PyProng3, hf_cand::PzProng3, collision::PosX, collision::PosY, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex>,
                  o2::soa::Marker<1>);

// extended table with expression columns that can be used as arguments of dynamic columns
DECLARE_SOA_EXTENDED_TABLE_USER(HfCandBsJPExt, HfCandBsJPBase, "HFCANDBSJPEXT",
                                hf_cand_4prong::Px, hf_cand_4prong::Py, hf_cand_4prong::Pz);

DECLARE_SOA_TABLE(HfCandBsJPDaus, "AOD", "HFCANDBSJPDAUS",
                  hf_cand_bs::Prong0Id, hf_track_index::Prong1Id, hf_track_index::Prong2Id, hf_track_index::Prong3Id);

using HfCandBsToJpsi = soa::Join<HfCandBsJPExt, HfCandBsJPDaus>;

// specific Σc0,++ candidate properties
namespace hf_cand_sigmac
{
DECLARE_SOA_INDEX_COLUMN_FULL(ProngLc, prongLc, int, HfCand3Prong, "");                //! Index to a Lc prong
DECLARE_SOA_COLUMN(Charge, charge, int8_t);                                            //! // Σc charge(either 0 or ++)
DECLARE_SOA_COLUMN(StatusSpreadLcMinvPKPiFromPDG, statusSpreadLcMinvPKPiFromPDG, int); //! // Λc Minv(pKpi) spread from PDG Λc mass
DECLARE_SOA_COLUMN(StatusSpreadLcMinvPiKPFromPDG, statusSpreadLcMinvPiKPFromPDG, int); //! // Λc Minv(piKp) spread from PDG Λc mass
DECLARE_SOA_COLUMN(SoftPiDcaXY, softPiDcaXY, float);                                   //! soft-pion impact parameter in xy
DECLARE_SOA_COLUMN(SoftPiDcaZ, softPiDcaZ, float);                                     //! soft-pion impact parameter in z
DECLARE_SOA_INDEX_COLUMN_FULL(Prong0, prong0, int, HfCand3Prong, "_0");                //! Λc index
// MC matching result:
DECLARE_SOA_COLUMN(FlagMcMatchRec, flagMcMatchRec, int8_t);             //! reconstruction level
DECLARE_SOA_COLUMN(FlagMcMatchGen, flagMcMatchGen, int8_t);             //! generator level
DECLARE_SOA_COLUMN(OriginMcRec, originMcRec, int8_t);                   //! particle origin, reconstruction level
DECLARE_SOA_COLUMN(OriginMcGen, originMcGen, int8_t);                   //! particle origin, generator level
DECLARE_SOA_COLUMN(ParticleAntiparticle, particleAntiparticle, int8_t); //! particle or antiparticle

enum Species : int { Sc2455 = 0,
                     Sc2520,
                     NSpecies };
enum Decays : int { PKPi = 0,
                    PiKP,
                    NDecays };
enum Conjugated : int { Particle = 0,
                        Antiparticle,
                        NConjugated };
constexpr int ChargeNull = 0;
constexpr int ChargePlusPlus = 2;
} // namespace hf_cand_sigmac

// declare dedicated Σc0,++ decay candidate table
// NB: no topology for Σc0,++ (strong decay)
DECLARE_SOA_TABLE(HfCandScBase, "AOD", "HFCANDSCBASE",
                  o2::soa::Index<>,
                  // general columns
                  hf_cand::CollisionId,
                  // 2-prong specific columns
                  hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0,
                  hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1,
                  // hf_cand::ImpactParameter0, hf_cand::ImpactParameter1,
                  // hf_cand::ErrorImpactParameter0, hf_cand::ErrorImpactParameter1,
                  // hf_track_index::ProngLcId, hf_track_index::Prong1Id,
                  hf_cand_sigmac::ProngLcId, hf_track_index::Prong1Id,
                  hf_track_index::HFflag,
                  /* Σc0,++ specific columns */
                  hf_cand_sigmac::Charge,
                  hf_cand_sigmac::StatusSpreadLcMinvPKPiFromPDG, hf_cand_sigmac::StatusSpreadLcMinvPiKPFromPDG,
                  hf_cand_sigmac::SoftPiDcaXY, hf_cand_sigmac::SoftPiDcaZ,
                  /* prong 0 */
                  // hf_cand::ImpactParameterNormalised0<hf_cand::ImpactParameter0, hf_cand::ErrorImpactParameter0>,
                  hf_cand::PtProng0<hf_cand::PxProng0, hf_cand::PyProng0>,
                  hf_cand::Pt2Prong0<hf_cand::PxProng0, hf_cand::PyProng0>,
                  hf_cand::PVectorProng0<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0>,
                  /* prong 1 */
                  // hf_cand::ImpactParameterNormalised1<hf_cand::ImpactParameter1, hf_cand::ErrorImpactParameter1>,
                  hf_cand::PtProng1<hf_cand::PxProng1, hf_cand::PyProng1>,
                  hf_cand::Pt2Prong1<hf_cand::PxProng1, hf_cand::PyProng1>,
                  hf_cand::PVectorProng1<hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>,
                  /* dynamic columns */
                  hf_cand_2prong::M<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>,
                  hf_cand_2prong::M2<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>,
                  /* dynamic columns that use candidate momentum components */
                  hf_cand::Pt<hf_cand_2prong::Px, hf_cand_2prong::Py>,
                  hf_cand::Pt2<hf_cand_2prong::Px, hf_cand_2prong::Py>,
                  hf_cand::P<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::P2<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::PVector<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::Eta<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::Phi<hf_cand_2prong::Px, hf_cand_2prong::Py>,
                  hf_cand::Y<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::E<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::E2<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>);

// extended table with expression columns that can be used as arguments of dynamic columns
DECLARE_SOA_EXTENDED_TABLE_USER(HfCandScExt, HfCandScBase, "HFCANDSCEXT",
                                hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz);
using HfCandSc = HfCandScExt;

// table with results of reconstruction level MC matching
DECLARE_SOA_TABLE(HfCandScMcRec, "AOD", "HFCANDSCMCREC", //!
                  hf_cand_sigmac::FlagMcMatchRec,
                  hf_cand_sigmac::OriginMcRec,
                  hf_cand::PtBhadMotherPart,
                  hf_cand::PdgBhadMotherPart,
                  hf_cand_sigmac::ParticleAntiparticle);

// table with results of generation level MC matching
DECLARE_SOA_TABLE(HfCandScMcGen, "AOD", "HFCANDSCMCGEN", //!
                  hf_cand_sigmac::FlagMcMatchGen,
                  hf_cand_sigmac::OriginMcGen,
                  hf_cand::IdxBhadMotherPart,
                  hf_cand_sigmac::ParticleAntiparticle);

// specific Σc0,++ candidate properties in cascade channel
namespace hf_cand_sigmac_to_cascade
{
DECLARE_SOA_INDEX_COLUMN_FULL(ProngLc, prongLc, int, HfCandCascade, "");               //! Index to a Lc prong
DECLARE_SOA_COLUMN(Charge, charge, int8_t);                                            //! // Σc charge(either 0 or ++)
DECLARE_SOA_COLUMN(ChargeLc, chargeLc, int8_t);                                        //! // Λc charge(+)
DECLARE_SOA_COLUMN(ChargeSoftPi, chargeSoftPi, int8_t);                                //! // pion charge(either - or +)
DECLARE_SOA_COLUMN(StatusSpreadLcMinvKs0PFromPDG, statusSpreadLcMinvKs0PFromPDG, int); //! // Λc Minv spread from PDG Λc mass
DECLARE_SOA_INDEX_COLUMN_FULL(Prong0, prong0, int, HfCandCascade, "_0");               //! Λc index
// MC matching result:
DECLARE_SOA_COLUMN(FlagMcMatchRec, flagMcMatchRec, int8_t); //! reconstruction level
DECLARE_SOA_COLUMN(FlagMcMatchGen, flagMcMatchGen, int8_t); //! generator level
DECLARE_SOA_COLUMN(OriginMcRec, originMcRec, int8_t);       //! particle origin, reconstruction level
DECLARE_SOA_COLUMN(OriginMcGen, originMcGen, int8_t);       //! particle origin, generator level
} // namespace hf_cand_sigmac_to_cascade

// declare dedicated Σc0,++ decay candidate table
// NB: no topology for Σc0,++ (strong decay)
DECLARE_SOA_TABLE(HfCandScCasBase, "AOD", "HFCANDSCCASBASE",
                  o2::soa::Index<>,
                  // general columns
                  hf_cand::CollisionId,
                  // 2-prong specific columns
                  hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0,
                  hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1,
                  hf_cand_sigmac_to_cascade::ProngLcId, hf_track_index::Prong1Id,
                  hf_cand_sigmac_to_cascade::ChargeLc,
                  hf_cand_sigmac_to_cascade::ChargeSoftPi,
                  // hf_track_index::HFflag,
                  /* Σc0,++ specific columns */
                  hf_cand_sigmac_to_cascade::Charge,
                  // hf_cand_sigmac_to_cascade::StatusSpreadLcMinvKs0PFromPDG,
                  /* prong 0 */
                  hf_cand::PtProng0<hf_cand::PxProng0, hf_cand::PyProng0>,
                  hf_cand::Pt2Prong0<hf_cand::PxProng0, hf_cand::PyProng0>,
                  hf_cand::PVectorProng0<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0>,
                  /* prong 1 */
                  hf_cand::PtProng1<hf_cand::PxProng1, hf_cand::PyProng1>,
                  hf_cand::Pt2Prong1<hf_cand::PxProng1, hf_cand::PyProng1>,
                  hf_cand::PVectorProng1<hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>,
                  /* dynamic columns */
                  hf_cand_2prong::M<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>,
                  hf_cand_2prong::M2<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>,
                  /* dynamic columns that use candidate momentum components */
                  hf_cand::Pt<hf_cand_2prong::Px, hf_cand_2prong::Py>,
                  hf_cand::Pt2<hf_cand_2prong::Px, hf_cand_2prong::Py>,
                  hf_cand::P<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::P2<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::PVector<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::Eta<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::Phi<hf_cand_2prong::Px, hf_cand_2prong::Py>,
                  hf_cand::Y<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::E<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::E2<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>);

// extended table with expression columns that can be used as arguments of dynamic columns
DECLARE_SOA_EXTENDED_TABLE_USER(HfCandScCasExt, HfCandScCasBase, "HFCANDSCCASEXT",
                                hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz);
using HfCandScCascades = HfCandScCasExt;
using HfCandScCascade = HfCandScCascades::iterator;

/// D*± → D0(bar) π±
namespace hf_cand_dstar
{
DECLARE_SOA_EXPRESSION_COLUMN(PxD0, pxD0, float, 1.f * aod::hf_cand::pxProng0 + 1.f * aod::hf_cand::pxProng1);
DECLARE_SOA_EXPRESSION_COLUMN(PyD0, pyD0, float, 1.f * aod::hf_cand::pyProng0 + 1.f * aod::hf_cand::pyProng1);
DECLARE_SOA_EXPRESSION_COLUMN(PzD0, pzD0, float, 1.f * aod::hf_cand::pzProng0 + 1.f * aod::hf_cand::pzProng1);
DECLARE_SOA_DYNAMIC_COLUMN(ImpactParameterProductD0, impactParameterProductD0,
                           [](float dca1, float dca2) -> float { return dca1 * dca2; });
// Dynamic Columns for D0 candidate using PDG masses of daughters
DECLARE_SOA_DYNAMIC_COLUMN(InvMassD0, invMassD0,
                           [](float px0, float py0, float pz0, float px1, float py1, float pz1) -> float { return RecoDecay::m(std::array{std::array{px0, py0, pz0}, std::array{px1, py1, pz1}}, std::array{constants::physics::MassPiPlus, constants::physics::MassKPlus}); });
DECLARE_SOA_DYNAMIC_COLUMN(InvMass2D0, invMass2D0,
                           [](float px0, float py0, float pz0, float px1, float py1, float pz1) -> float { return RecoDecay::m2(std::array{std::array{px0, py0, pz0}, std::array{px1, py1, pz1}}, std::array{constants::physics::MassPiPlus, constants::physics::MassKPlus}); });
DECLARE_SOA_DYNAMIC_COLUMN(CosThetaStarD0, cosThetaStarD0,
                           [](float px0, float py0, float pz0, float px1, float py1, float pz1) -> float { return RecoDecay::cosThetaStar(std::array{std::array{px0, py0, pz0}, std::array{px1, py1, pz1}}, std::array{constants::physics::MassPiPlus, constants::physics::MassKPlus}, constants::physics::MassD0, 1); });
// Dynamic Columns for D0Bar candidate using PDG masses of daughters
DECLARE_SOA_DYNAMIC_COLUMN(InvMassD0Bar, invMassD0Bar,
                           [](float px0, float py0, float pz0, float px1, float py1, float pz1) -> float { return RecoDecay::m(std::array{std::array{px0, py0, pz0}, std::array{px1, py1, pz1}}, std::array{constants::physics::MassKPlus, constants::physics::MassPiPlus}); });
DECLARE_SOA_DYNAMIC_COLUMN(InvMass2D0Bar, invMass2D0Bar,
                           [](float px0, float py0, float pz0, float px1, float py1, float pz1) -> float { return RecoDecay::m2(std::array{std::array{px0, py0, pz0}, std::array{px1, py1, pz1}}, std::array{constants::physics::MassKPlus, constants::physics::MassPiPlus}); });
DECLARE_SOA_DYNAMIC_COLUMN(CosThetaStarD0Bar, cosThetaStarD0Bar,
                           [](float px0, float py0, float pz0, float px1, float py1, float pz1) -> float { return RecoDecay::cosThetaStar(std::array{std::array{px0, py0, pz0}, std::array{px1, py1, pz1}}, std::array{constants::physics::MassKPlus, constants::physics::MassPiPlus}, constants::physics::MassD0, 0); });
DECLARE_SOA_DYNAMIC_COLUMN(ImpactParameterProngSqSumD0, impactParameterProngSqSumD0,
                           [](float impParProng0, float impParProng1) -> float { return RecoDecay::sumOfSquares(impParProng0, impParProng1); });
DECLARE_SOA_DYNAMIC_COLUMN(DeltaIPNormalisedMaxD0, deltaIPNormalisedMaxD0,
                           [](float xVtxP, float yVtxP, float xVtxS, float yVtxS, float errDlxy, float pxM, float pyM, float ip0, float errIp0, float ip1, float errIp1, float px0, float py0, float px1, float py1) -> float { return RecoDecay::maxNormalisedDeltaIP(std::array{xVtxP, yVtxP}, std::array{xVtxS, yVtxS}, errDlxy, std::array{pxM, pyM}, std::array{ip0, ip1}, std::array{errIp0, errIp1}, std::array{std::array{px0, py0}, std::array{px1, py1}}); });
DECLARE_SOA_DYNAMIC_COLUMN(PtD0, ptD0,
                           [](float px, float py) -> float { return RecoDecay::pt(px, py); });
DECLARE_SOA_DYNAMIC_COLUMN(Pt2D0, pt2D0,
                           [](float px, float py) -> float { return RecoDecay::pt2(px, py); });
DECLARE_SOA_DYNAMIC_COLUMN(PD0, pD0,
                           [](float px, float py, float pz) -> float { return RecoDecay::p(px, py, pz); });
DECLARE_SOA_DYNAMIC_COLUMN(P2D0, p2D0,
                           [](float px, float py, float pz) -> float { return RecoDecay::p2(px, py, pz); });
DECLARE_SOA_DYNAMIC_COLUMN(PVectorD0, pVectorD0,
                           [](float px, float py, float pz) -> std::array<float, 3> { return std::array{px, py, pz}; });
DECLARE_SOA_DYNAMIC_COLUMN(EtaD0, etaD0,
                           [](float px, float py, float pz) -> float { return RecoDecay::eta(std::array{px, py, pz}); });
DECLARE_SOA_DYNAMIC_COLUMN(PhiD0, phiD0,
                           [](float px, float py) -> float { return RecoDecay::phi(px, py); });
DECLARE_SOA_DYNAMIC_COLUMN(YD0, yD0,
                           [](float px, float py, float pz) -> float { return RecoDecay::y(std::array{px, py, pz}, constants::physics::MassD0); });
DECLARE_SOA_DYNAMIC_COLUMN(ED0, eD0,
                           [](float px, float py, float pz) -> float { return RecoDecay::e(px, py, pz, constants::physics::MassD0); });
DECLARE_SOA_DYNAMIC_COLUMN(E2D0, e2D0,
                           [](float px, float py, float pz) -> float { return RecoDecay::e2(px, py, pz, constants::physics::MassD0); });
// secondary vertex
DECLARE_SOA_COLUMN(Chi2PCAD0, chi2PCAD0, float);
DECLARE_SOA_COLUMN(XSecondaryVertexD0, xSecondaryVertexD0, float);
DECLARE_SOA_COLUMN(YSecondaryVertexD0, ySecondaryVertexD0, float);
DECLARE_SOA_COLUMN(ZSecondaryVertexD0, zSecondaryVertexD0, float);
DECLARE_SOA_DYNAMIC_COLUMN(RSecondaryVertexD0, rSecondaryVertexD0,
                           [](float xVtxS, float yVtxS) -> float { return RecoDecay::sqrtSumOfSquares(xVtxS, yVtxS); });
DECLARE_SOA_DYNAMIC_COLUMN(DecayLengthD0, decayLengthD0,
                           [](float xVtxP, float yVtxP, float zVtxP, float xVtxS, float yVtxS, float zVtxS) -> float { return RecoDecay::distance(std::array{xVtxP, yVtxP, zVtxP}, std::array{xVtxS, yVtxS, zVtxS}); });
DECLARE_SOA_DYNAMIC_COLUMN(DecayLengthXYD0, decayLengthXYD0,
                           [](float xVtxP, float yVtxP, float xVtxS, float yVtxS) -> float { return RecoDecay::distanceXY(std::array{xVtxP, yVtxP}, std::array{xVtxS, yVtxS}); });
DECLARE_SOA_DYNAMIC_COLUMN(DecayLengthNormalisedD0, decayLengthNormalisedD0,
                           [](float xVtxP, float yVtxP, float zVtxP, float xVtxS, float yVtxS, float zVtxS, float err) -> float { return RecoDecay::distance(std::array{xVtxP, yVtxP, zVtxP}, std::array{xVtxS, yVtxS, zVtxS}) / err; });
DECLARE_SOA_DYNAMIC_COLUMN(DecayLengthXYNormalisedD0, decayLengthXYNormalisedD0,
                           [](float xVtxP, float yVtxP, float xVtxS, float yVtxS, float err) -> float { return RecoDecay::distanceXY(std::array{xVtxP, yVtxP}, std::array{xVtxS, yVtxS}) / err; });
DECLARE_SOA_COLUMN(ErrorDecayLengthD0, errorDecayLengthD0, float);
DECLARE_SOA_COLUMN(ErrorDecayLengthXYD0, errorDecayLengthXYD0, float);
DECLARE_SOA_DYNAMIC_COLUMN(CpaD0, cpaD0,
                           [](float xVtxP, float yVtxP, float zVtxP, float xVtxS, float yVtxS, float zVtxS, float px, float py, float pz) -> float { return RecoDecay::cpa(std::array{xVtxP, yVtxP, zVtxP}, std::array{xVtxS, yVtxS, zVtxS}, std::array{px, py, pz}); });
DECLARE_SOA_DYNAMIC_COLUMN(CpaXYD0, cpaXYD0,
                           [](float xVtxP, float yVtxP, float xVtxS, float yVtxS, float px, float py) -> float { return RecoDecay::cpaXY(std::array{xVtxP, yVtxP}, std::array{xVtxS, yVtxS}, std::array{px, py}); });
DECLARE_SOA_DYNAMIC_COLUMN(CtD0, ctD0,
                           [](float xVtxP, float yVtxP, float zVtxP, float xVtxS, float yVtxS, float zVtxS, float px, float py, float pz) -> float { return RecoDecay::ct(std::array{px, py, pz}, RecoDecay::distance(std::array{xVtxP, yVtxP, zVtxP}, std::array{xVtxS, yVtxS, zVtxS}), constants::physics::MassD0); });
DECLARE_SOA_DYNAMIC_COLUMN(ImpactParameterXYD0, impactParameterXYD0,
                           [](float xVtxP, float yVtxP, float zVtxP, float xVtxS, float yVtxS, float zVtxS, float px, float py, float pz) -> float { return RecoDecay::impParXY(std::array{xVtxP, yVtxP, zVtxP}, std::array{xVtxS, yVtxS, zVtxS}, std::array{px, py, pz}); });

// Columns only for D* properties
DECLARE_SOA_INDEX_COLUMN_FULL(ProngPi, prongPi, int, Tracks, ""); //! soft-pion index

// soft pion prong
DECLARE_SOA_COLUMN(ImpParamSoftPi, impParamSoftPi, float);
DECLARE_SOA_COLUMN(ImpParamZSoftPi, impParamZSoftPi, float);
DECLARE_SOA_COLUMN(ErrorImpParamSoftPi, errorImpParamSoftPi, float);
DECLARE_SOA_COLUMN(ErrorImpParamZSoftPi, errorImpParamZSoftPi, float);
DECLARE_SOA_DYNAMIC_COLUMN(NormalisedImpParamSoftPi, normalisedImpParamSoftPi,
                           [](float dca, float err) -> float { return dca / err; });
DECLARE_SOA_DYNAMIC_COLUMN(NormalisedImpParamZSoftPi, normalisedImpParamZSoftPi,
                           [](float dca, float err) -> float { return dca / err; });
DECLARE_SOA_COLUMN(PxSoftPi, pxSoftPi, float);
DECLARE_SOA_COLUMN(PySoftPi, pySoftPi, float);
DECLARE_SOA_COLUMN(PzSoftPi, pzSoftPi, float);
DECLARE_SOA_COLUMN(DcaYSoftPi, dcaYSoftPi, float);
DECLARE_SOA_COLUMN(SigmaYSoftPi, sigmaYSoftPi, float);
DECLARE_SOA_COLUMN(SignSoftPi, signSoftPi, int8_t);
DECLARE_SOA_COLUMN(TPCNSigmaPiSoftPi, tpcNSigmaPiSoftPi, float); //! NsigmaTPCPi for soft pi, o2-linter: disable=name/o2-column (written to disk)
DECLARE_SOA_COLUMN(TOFNSigmaPiSoftPi, tofNSigmaPiSoftPi, float); //! NsigmaTOFPi for soft pi, o2-linter: disable=name/o2-column (written to disk)
// Dstar momenta
DECLARE_SOA_EXPRESSION_COLUMN(PxDstar, pxDstar, float, 1.f * aod::hf_cand::pxProng0 + 1.f * aod::hf_cand::pxProng1 + 1.f * aod::hf_cand_dstar::pxSoftPi);
DECLARE_SOA_EXPRESSION_COLUMN(PyDstar, pyDstar, float, 1.f * aod::hf_cand::pyProng0 + 1.f * aod::hf_cand::pyProng1 + 1.f * aod::hf_cand_dstar::pySoftPi);
DECLARE_SOA_EXPRESSION_COLUMN(PzDstar, pzDstar, float, 1.f * aod::hf_cand::pzProng0 + 1.f * aod::hf_cand::pzProng1 + 1.f * aod::hf_cand_dstar::pzSoftPi);
// Inv Mass (accept mass array of size 3 {π , π, k})
DECLARE_SOA_DYNAMIC_COLUMN(InvMassDstar, invMassDstar,
                           [](float pxSoftPi, float pySoftPi, float pzSoftPi, float pxProng0, float pyProng0, float pzProng0, float pxProng1, float pyProng1, float pzProng1)
                             -> float { return RecoDecay::m(std::array{std::array{pxSoftPi, pySoftPi, pzSoftPi}, std::array{pxProng0, pyProng0, pzProng0}, std::array{pxProng1, pyProng1, pzProng1}}, std::array{constants::physics::MassPiPlus, constants::physics::MassPiPlus, constants::physics::MassKPlus}); });

DECLARE_SOA_DYNAMIC_COLUMN(InvMassAntiDstar, invMassAntiDstar,
                           [](float pxSoftPi, float pySoftPi, float pzSoftPi, float pxProng0, float pyProng0, float pzProng0, float pxProng1, float pyProng1, float pzProng1)
                             -> float { return RecoDecay::m(std::array{std::array{pxSoftPi, pySoftPi, pzSoftPi}, std::array{pxProng0, pyProng0, pzProng0}, std::array{pxProng1, pyProng1, pzProng1}}, std::array{constants::physics::MassPiPlus, constants::physics::MassKPlus, constants::physics::MassPiPlus}); });

DECLARE_SOA_DYNAMIC_COLUMN(PtSoftPi, ptSoftPi, [](float pxSoftPi, float pySoftPi) -> float { return RecoDecay::pt(pxSoftPi, pySoftPi); });
DECLARE_SOA_DYNAMIC_COLUMN(PVecSoftPi, pVecSoftPi, [](float px, float py, float pz) -> std::array<float, 3> { return std::array{px, py, pz}; });
DECLARE_SOA_DYNAMIC_COLUMN(TPCTOFNSigmaPiSoftPi, tpcTofNSigmaPiSoftPi, //! Combination of NsigmaTPC and NsigmaTOF, o2-linter: disable=name/o2-column (written to disk)
                           [](float tpcNSigmaPiSoftPi, float TOFNSigmaPiSoftPi) -> float { return pid_tpc_tof_utils::combineNSigma<false /*tiny*/>(tpcNSigmaPiSoftPi, TOFNSigmaPiSoftPi); });
// MC matching result:
DECLARE_SOA_COLUMN(FlagMcMatchRec, flagMcMatchRec, int8_t);     //! reconstruction level
DECLARE_SOA_COLUMN(FlagMcMatchGen, flagMcMatchGen, int8_t);     //! generator level
DECLARE_SOA_COLUMN(FlagMcMatchRecD0, flagMcMatchRecD0, int8_t); //! reconstruction level
DECLARE_SOA_COLUMN(FlagMcMatchGenD0, flagMcMatchGenD0, int8_t); //! generator level

DECLARE_SOA_COLUMN(OriginMcRec, originMcRec, int8_t); //! particle origin, reconstruction level
DECLARE_SOA_COLUMN(OriginMcGen, originMcGen, int8_t); //! particle origin, generator level

} // namespace hf_cand_dstar

/// D0 (table) from DStar
DECLARE_SOA_TABLE(HfD0FromDstarBase, "AOD", "HFD0FRMDSTR",
                  o2::soa::Index<>,
                  // gener columns
                  hf_cand::CollisionId,
                  collision::PosX, collision::PosY, collision::PosZ,
                  hf_cand_dstar::XSecondaryVertexD0, hf_cand_dstar::YSecondaryVertexD0, hf_cand_dstar::ZSecondaryVertexD0,
                  hf_cand_dstar::ErrorDecayLengthD0, hf_cand_dstar::ErrorDecayLengthXYD0,
                  hf_cand_dstar::Chi2PCAD0,
                  /* dynamic columns */ hf_cand_dstar::RSecondaryVertexD0<hf_cand_dstar::XSecondaryVertexD0, hf_cand_dstar::YSecondaryVertexD0>,
                  hf_cand_dstar::DecayLengthD0<collision::PosX, collision::PosY, collision::PosZ, hf_cand_dstar::XSecondaryVertexD0, hf_cand_dstar::YSecondaryVertexD0, hf_cand_dstar::ZSecondaryVertexD0>,
                  hf_cand_dstar::DecayLengthXYD0<collision::PosX, collision::PosY, hf_cand_dstar::XSecondaryVertexD0, hf_cand_dstar::YSecondaryVertexD0>,
                  hf_cand_dstar::DecayLengthNormalisedD0<collision::PosX, collision::PosY, collision::PosZ, hf_cand_dstar::XSecondaryVertexD0, hf_cand_dstar::YSecondaryVertexD0, hf_cand_dstar::ZSecondaryVertexD0, hf_cand_dstar::ErrorDecayLengthD0>,
                  hf_cand_dstar::DecayLengthXYNormalisedD0<collision::PosX, collision::PosY, hf_cand_dstar::XSecondaryVertexD0, hf_cand_dstar::YSecondaryVertexD0, hf_cand_dstar::ErrorDecayLengthXYD0>,
                  /* prong 0 */ hf_cand::ImpactParameterNormalised0<hf_cand::ImpactParameter0, hf_cand::ErrorImpactParameter0>,
                  /* prong 1 */ hf_cand::ImpactParameterNormalised1<hf_cand::ImpactParameter1, hf_cand::ErrorImpactParameter1>,
                  /* prong 0 */ hf_cand::ImpactParameterZNormalised0<hf_cand::ImpactParameterZ0, hf_cand::ErrorImpactParameterZ0>,
                  /* prong 1 */ hf_cand::ImpactParameterZNormalised1<hf_cand::ImpactParameterZ1, hf_cand::ErrorImpactParameterZ1>,

                  // HFCAND_COLUMNS,
                  // 2-prong specific columns
                  hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0,
                  hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1,
                  hf_cand::ImpactParameter0, hf_cand::ImpactParameter1,
                  hf_cand::ImpactParameterZ0, hf_cand::ImpactParameterZ1,
                  hf_cand::ErrorImpactParameter0, hf_cand::ErrorImpactParameter1,
                  hf_cand::ErrorImpactParameterZ0, hf_cand::ErrorImpactParameterZ1,
                  hf_track_index::Prong0Id, hf_track_index::Prong1Id,
                  hf_track_index::HFflag,
                  /* dynamic columns */
                  hf_cand_dstar::ImpactParameterProductD0<hf_cand::ImpactParameter0, hf_cand::ImpactParameter1>,
                  hf_cand_dstar::ImpactParameterProngSqSumD0<hf_cand::ImpactParameter0, hf_cand::ImpactParameter1>,
                  /* dynamic columns that use candidate momentum components */
                  hf_cand_dstar::PtD0<hf_cand_dstar::PxD0, hf_cand_dstar::PyD0>,
                  hf_cand_dstar::Pt2D0<hf_cand_dstar::PxD0, hf_cand_dstar::PyD0>,
                  hf_cand_dstar::PD0<hf_cand_dstar::PxD0, hf_cand_dstar::PyD0, hf_cand_dstar::PzD0>,
                  hf_cand_dstar::P2D0<hf_cand_dstar::PxD0, hf_cand_dstar::PyD0, hf_cand_dstar::PzD0>,
                  hf_cand_dstar::PVectorD0<hf_cand_dstar::PxD0, hf_cand_dstar::PyD0, hf_cand_dstar::PzD0>,
                  hf_cand_dstar::CpaD0<collision::PosX, collision::PosY, collision::PosZ, hf_cand_dstar::XSecondaryVertexD0, hf_cand_dstar::YSecondaryVertexD0, hf_cand_dstar::ZSecondaryVertexD0, hf_cand_dstar::PxD0, hf_cand_dstar::PyD0, hf_cand_dstar::PzD0>,
                  hf_cand_dstar::CpaXYD0<collision::PosX, collision::PosY, hf_cand_dstar::XSecondaryVertexD0, hf_cand_dstar::YSecondaryVertexD0, hf_cand_dstar::PxD0, hf_cand_dstar::PyD0>,
                  hf_cand_dstar::CtD0<collision::PosX, collision::PosY, collision::PosZ, hf_cand_dstar::XSecondaryVertexD0, hf_cand_dstar::YSecondaryVertexD0, hf_cand_dstar::ZSecondaryVertexD0, hf_cand_dstar::PxD0, hf_cand_dstar::PyD0, hf_cand_dstar::PzD0>,
                  hf_cand_dstar::ImpactParameterXYD0<collision::PosX, collision::PosY, collision::PosZ, hf_cand_dstar::XSecondaryVertexD0, hf_cand_dstar::YSecondaryVertexD0, hf_cand_dstar::ZSecondaryVertexD0, hf_cand_dstar::PxD0, hf_cand_dstar::PyD0, hf_cand_dstar::PzD0>,
                  hf_cand_dstar::DeltaIPNormalisedMaxD0<collision::PosX, collision::PosY, hf_cand_dstar::XSecondaryVertexD0, hf_cand_dstar::YSecondaryVertexD0, hf_cand_dstar::ErrorDecayLengthXYD0, hf_cand_dstar::PxD0, hf_cand_dstar::PyD0, hf_cand::ImpactParameter0, hf_cand::ErrorImpactParameter0, hf_cand::ImpactParameter1, hf_cand::ErrorImpactParameter1, hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PxProng1, hf_cand::PyProng1>,
                  hf_cand_dstar::EtaD0<hf_cand_dstar::PxD0, hf_cand_dstar::PyD0, hf_cand_dstar::PzD0>,
                  hf_cand_dstar::PhiD0<hf_cand_dstar::PxD0, hf_cand_dstar::PyD0>,
                  hf_cand_dstar::YD0<hf_cand_dstar::PxD0, hf_cand_dstar::PyD0, hf_cand_dstar::PzD0>,
                  hf_cand_dstar::ED0<hf_cand_dstar::PxD0, hf_cand_dstar::PyD0, hf_cand_dstar::PzD0>,
                  hf_cand_dstar::E2D0<hf_cand_dstar::PxD0, hf_cand_dstar::PyD0, hf_cand_dstar::PzD0>);

// extended table with expression columns that can be used as arguments of dynamic columns
DECLARE_SOA_EXTENDED_TABLE_USER(HfD0FromDstarExt, HfD0FromDstarBase, "HFD0FRMDSTREXT",
                                hf_cand_dstar::PxD0, hf_cand_dstar::PyD0, hf_cand_dstar::PzD0);

DECLARE_SOA_TABLE(HfCandDstarProng0PidPi, "AOD", "HFDSTRP0PIDPI", //!
                  hf_cand::NSigTpcPi0, hf_cand::NSigTofPi0,
                  hf_cand::TpcTofNSigmaPi0<hf_cand::NSigTpcPi0, hf_cand::NSigTofPi0>);
DECLARE_SOA_TABLE(HfCandDstarProng1PidPi, "AOD", "HFDSTRP1PIDPI", //!
                  hf_cand::NSigTpcPi1, hf_cand::NSigTofPi1,
                  hf_cand::TpcTofNSigmaPi1<hf_cand::NSigTpcPi1, hf_cand::NSigTofPi1>);
DECLARE_SOA_TABLE(HfCandDstarProng2PidPi, "AOD", "HFDSTRP2PIDPI", //!
                  hf_cand::NSigTpcPi2, hf_cand::NSigTofPi2,
                  hf_cand::TpcTofNSigmaPi2<hf_cand::NSigTpcPi2, hf_cand::NSigTofPi2>);
DECLARE_SOA_TABLE(HfCandDstarProng0PidKa, "AOD", "HFDSTRP0PIDKA", //!
                  hf_cand::NSigTpcKa0, hf_cand::NSigTofKa0,
                  hf_cand::TpcTofNSigmaKa0<hf_cand::NSigTpcKa0, hf_cand::NSigTofKa0>);
DECLARE_SOA_TABLE(HfCandDstarProng1PidKa, "AOD", "HFDSTRP1PIDKA", //!
                  hf_cand::NSigTpcKa1, hf_cand::NSigTofKa1,
                  hf_cand::TpcTofNSigmaKa1<hf_cand::NSigTpcKa1, hf_cand::NSigTofKa1>);
DECLARE_SOA_TABLE(HfCandDstarProng2PidKa, "AOD", "HFDSTRP2PIDKA", //!
                  hf_cand::NSigTpcKa2, hf_cand::NSigTofKa2,
                  hf_cand::TpcTofNSigmaKa2<hf_cand::NSigTpcKa2, hf_cand::NSigTofKa2>);

using HfD0FromDstar = HfD0FromDstarExt;
using HfD0FromDstarWPid = soa::Join<HfD0FromDstar, HfCandDstarProng0PidPi, HfCandDstarProng0PidKa, HfCandDstarProng1PidPi, HfCandDstarProng1PidKa>;

DECLARE_SOA_TABLE(HfCandDstarBase, "AOD", "HFCANDDSTRBASE",
                  o2::soa::Index<>,
                  hf_cand::CollisionId,
                  // Primary vertex
                  collision::PosX, collision::PosY, collision::PosZ,
                  hf_cand_dstar::ProngPiId,  // Index column to softPi track table
                  hf_track_index::ProngD0Id, // Index column points to Hf2Prongs table filled by indexSkimcreator
                  // Softpi
                  hf_cand_dstar::PxSoftPi, hf_cand_dstar::PySoftPi, hf_cand_dstar::PzSoftPi,
                  hf_cand_dstar::SignSoftPi,
                  hf_cand_dstar::ImpParamSoftPi, hf_cand_dstar::ImpParamZSoftPi,
                  hf_cand_dstar::ErrorImpParamSoftPi, hf_cand_dstar::ErrorImpParamZSoftPi,
                  // Two pronges of D0
                  hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0,
                  hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1,
                  hf_track_index::Prong0Id, hf_track_index::Prong1Id,
                  // Dynamic
                  hf_cand_dstar::PtSoftPi<hf_cand_dstar::PxSoftPi, hf_cand_dstar::PySoftPi>,
                  hf_cand_dstar::PVecSoftPi<hf_cand_dstar::PxSoftPi, hf_cand_dstar::PySoftPi, hf_cand_dstar::PzSoftPi>,
                  hf_cand_dstar::NormalisedImpParamSoftPi<hf_cand_dstar::ImpParamSoftPi, hf_cand_dstar::ErrorImpParamSoftPi>,
                  hf_cand_dstar::NormalisedImpParamZSoftPi<hf_cand_dstar::ImpParamZSoftPi, hf_cand_dstar::ErrorImpParamZSoftPi>,
                  hf_cand::Pt<hf_cand_dstar::PxDstar, hf_cand_dstar::PyDstar>,
                  hf_cand::P<hf_cand_dstar::PxDstar, hf_cand_dstar::PyDstar, hf_cand_dstar::PzDstar>,
                  hf_cand::PVector<hf_cand_dstar::PxDstar, hf_cand_dstar::PyDstar, hf_cand_dstar::PzDstar>,
                  hf_cand::Eta<hf_cand_dstar::PxDstar, hf_cand_dstar::PyDstar, hf_cand_dstar::PzDstar>,
                  hf_cand::Phi<hf_cand_dstar::PxDstar, hf_cand_dstar::PyDstar>,
                  hf_cand::Y<hf_cand_dstar::PxDstar, hf_cand_dstar::PyDstar, hf_cand_dstar::PzDstar>,
                  hf_cand::E<hf_cand_dstar::PxDstar, hf_cand_dstar::PyDstar, hf_cand_dstar::PzDstar>,
                  hf_cand::PtProng0<hf_cand::PxProng0, hf_cand::PyProng0>,
                  hf_cand::Pt2Prong0<hf_cand::PxProng0, hf_cand::PyProng0>,
                  hf_cand::PVectorProng0<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0>,
                  hf_cand::PtProng1<hf_cand::PxProng1, hf_cand::PyProng1>,
                  hf_cand::Pt2Prong1<hf_cand::PxProng1, hf_cand::PyProng1>,
                  hf_cand::PVectorProng1<hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>,
                  hf_cand_dstar::InvMassDstar<hf_cand_dstar::PxSoftPi, hf_cand_dstar::PySoftPi, hf_cand_dstar::PzSoftPi, hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>,
                  hf_cand_dstar::InvMassAntiDstar<hf_cand_dstar::PxSoftPi, hf_cand_dstar::PySoftPi, hf_cand_dstar::PzSoftPi, hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>,
                  hf_cand_dstar::InvMassD0<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>,
                  hf_cand_dstar::InvMassD0Bar<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>,
                  hf_cand_dstar::InvMass2D0<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>,
                  hf_cand_dstar::InvMass2D0Bar<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>,
                  hf_cand_dstar::CosThetaStarD0<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>,
                  hf_cand_dstar::CosThetaStarD0Bar<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>);

// extended table with expression columns that can be used as arguments of dynamic columns
DECLARE_SOA_EXTENDED_TABLE_USER(HfCandDstarExt, HfCandDstarBase, "HFCANDDSTREXT",
                                hf_cand_dstar::PxDstar, hf_cand_dstar::PyDstar, hf_cand_dstar::PzDstar);

using HfCandDstars = HfCandDstarExt;
using HfCandDstar = HfCandDstars::iterator;
using HfCandDstarsWPid = soa::Join<HfCandDstars, HfCandDstarProng0PidPi, HfCandDstarProng0PidKa, HfCandDstarProng1PidPi, HfCandDstarProng1PidKa, /*soft pion*/ HfCandDstarProng2PidPi, /*soft pion*/ HfCandDstarProng2PidKa>;

// table with results of reconstruction level MC matching
DECLARE_SOA_TABLE(HfCandDstarMcRec, "AOD", "HFCANDDSTRMCREC",
                  hf_cand_dstar::FlagMcMatchRec,
                  hf_cand_dstar::FlagMcMatchRecD0,
                  hf_cand_dstar::OriginMcRec,
                  hf_cand::PtBhadMotherPart,
                  hf_cand::PdgBhadMotherPart,
                  hf_cand::NTracksDecayed,
                  hf_cand::NInteractionsWithMaterial);

// table with results of generator level MC matching
DECLARE_SOA_TABLE(HfCandDstarMcGen, "AOD", "HFCANDDSTRMCGEN",
                  hf_cand_dstar::FlagMcMatchGen,
                  hf_cand_dstar::FlagMcMatchGenD0,
                  hf_cand_dstar::OriginMcGen,
                  hf_cand::IdxBhadMotherPart);

#undef HFCAND_COLUMNS

} // namespace o2::aod

#endif // PWGHF_DATAMODEL_CANDIDATERECONSTRUCTIONTABLES_H_
