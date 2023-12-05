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
/// \file LFHypernucleiTables.h
/// \brief Slim hypernuclei tables
///

#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"

#ifndef PWGLF_DATAMODEL_LFHYPERNUCLEITABLES_H_
#define PWGLF_DATAMODEL_LFHYPERNUCLEITABLES_H_

namespace o2::aod
{
namespace hyperrec
{
DECLARE_SOA_COLUMN(IsMatter, isMatter, bool);                       // bool: true for matter
DECLARE_SOA_COLUMN(CentralityFT0A, centralityFT0A, float);          // centrality with FT0A estimator
DECLARE_SOA_COLUMN(CentralityFT0C, centralityFT0C, float);          // centrality with FT0C estimator
DECLARE_SOA_COLUMN(CentralityFT0M, centralityFT0M, float);          // centrality with FT0M estimator
DECLARE_SOA_COLUMN(PtHe3, ptHe3, float);                            // Pt of the He daughter
DECLARE_SOA_COLUMN(PhiHe3, phiHe3, float);                          // Phi of the He daughter
DECLARE_SOA_COLUMN(EtaHe3, etaHe3, float);                          // Eta of the He daughter
DECLARE_SOA_COLUMN(PtPi, ptPi, float);                              // Pt of the Pi daughter
DECLARE_SOA_COLUMN(PhiPi, phiPi, float);                            // Phi of the Pi daughter
DECLARE_SOA_COLUMN(EtaPi, etaPi, float);                            // Eta of the Pi daughter
DECLARE_SOA_COLUMN(XPrimVtx, xPrimVtx, float);                      // Decay vertex of the candidate (x direction)
DECLARE_SOA_COLUMN(YPrimVtx, yPrimVtx, float);                      // Decay vertex of the candidate (y direction)
DECLARE_SOA_COLUMN(ZPrimVtx, zPrimVtx, float);                      // Decay vertex of the candidate (z direction)
DECLARE_SOA_COLUMN(XDecVtx, xDecVtx, float);                        // Decay vertex of the candidate (x direction)
DECLARE_SOA_COLUMN(YDecVtx, yDecVtx, float);                        // Decay vertex of the candidate (y direction)
DECLARE_SOA_COLUMN(ZDecVtx, zDecVtx, float);                        // Decay vertex of the candidate (z direction)
DECLARE_SOA_COLUMN(MassH3L, massH3L, float);                        // Squared mass w/ hypertriton mass hypo
DECLARE_SOA_COLUMN(MassH4L, massH4L, float);                        // Squared mass w/ H4L mass hypo
DECLARE_SOA_COLUMN(DcaV0Daug, dcaV0Daug, float);                    // DCA between daughters
DECLARE_SOA_COLUMN(CosPA, cosPA, double);                           // Cosine of the pointing angle
DECLARE_SOA_COLUMN(NSigmaHe, nSigmaHe, float);                      // Number of sigmas of the He daughter
DECLARE_SOA_COLUMN(NTPCclusHe, nTPCclusHe, uint8_t);                // Number of TPC clusters of the He daughter
DECLARE_SOA_COLUMN(NTPCclusPi, nTPCclusPi, uint8_t);                // Number of TPC clusters of the Pi daughter
DECLARE_SOA_COLUMN(TPCsignalHe, tpcSignalHe, uint16_t);             // TPC signal of the He daughter
DECLARE_SOA_COLUMN(TPCsignalPi, tpcSignalPi, uint16_t);             // TPC signal of the Pi daughter
DECLARE_SOA_COLUMN(Flags, flags, uint8_t);                          // Flags for PID in tracking (bits [0, 3] for negative daughter, [4,7] for positive daughter)
DECLARE_SOA_COLUMN(TPCmomHe, tpcMomHe, float);                      // TPC momentum of the He daughter
DECLARE_SOA_COLUMN(TPCmomPi, tpcMomPi, float);                      // TPC momentum of the Pi daughter
DECLARE_SOA_COLUMN(ITSclusterSizesHe, itsClusterSizesHe, uint32_t); // ITS cluster size of the He daughter
DECLARE_SOA_COLUMN(ITSclusterSizesPi, itsClusterSizesPi, uint32_t); // ITS cluster size of the Pi daughter
DECLARE_SOA_COLUMN(DcaHe, dcaHe, float);                            // DCA between He daughter and V0
DECLARE_SOA_COLUMN(DcaPi, dcaPi, float);                            // DCA between pi daughter and V0
DECLARE_SOA_COLUMN(GenPt, genPt, float);                            // Pt of the hypertriton
DECLARE_SOA_COLUMN(GenPhi, genPhi, float);                          // Phi of the hypertriton
DECLARE_SOA_COLUMN(GenEta, genEta, float);                          // Eta of the hypertriton
DECLARE_SOA_COLUMN(GenPtHe3, genPtHe3, float);                      // Pt of the He daughter (to be used for the recalibration)
DECLARE_SOA_COLUMN(GenXDecVtx, genXDecVtx, float);                  // Decay vertex of the candidate (x direction)
DECLARE_SOA_COLUMN(GenYDecVtx, genYDecVtx, float);                  // Decay vertex of the candidate (y direction)
DECLARE_SOA_COLUMN(GenZDecVtx, genZDecVtx, float);                  // Decay vertex of the candidate (z direction)
DECLARE_SOA_COLUMN(IsReco, isReco, bool);                           // bool: true for reco
DECLARE_SOA_COLUMN(IsSignal, isSignal, bool);                       // bool: true for signal
} // namespace hyperrec

DECLARE_SOA_TABLE(DataHypCands, "AOD", "DATAHYPCANDS",
                  o2::soa::Index<>,
                  hyperrec::IsMatter,
                  hyperrec::CentralityFT0A,
                  hyperrec::CentralityFT0C,
                  hyperrec::CentralityFT0M,
                  hyperrec::PtHe3,
                  hyperrec::PhiHe3,
                  hyperrec::EtaHe3,
                  hyperrec::PtPi,
                  hyperrec::PhiPi,
                  hyperrec::EtaPi,
                  hyperrec::XPrimVtx,
                  hyperrec::YPrimVtx,
                  hyperrec::ZPrimVtx,
                  hyperrec::XDecVtx,
                  hyperrec::YDecVtx,
                  hyperrec::ZDecVtx,
                  hyperrec::DcaV0Daug,
                  hyperrec::DcaHe,
                  hyperrec::DcaPi,
                  hyperrec::NSigmaHe,
                  hyperrec::NTPCclusHe,
                  hyperrec::NTPCclusPi,
                  hyperrec::TPCmomHe,
                  hyperrec::TPCmomPi,
                  hyperrec::TPCsignalHe,
                  hyperrec::TPCsignalPi,
                  hyperrec::ITSclusterSizesHe,
                  hyperrec::ITSclusterSizesPi,
                  hyperrec::Flags);

DECLARE_SOA_TABLE(MCHypCands, "AOD", "MCHYPCANDS",
                  o2::soa::Index<>,
                  hyperrec::IsMatter,
                  hyperrec::CentralityFT0A,
                  hyperrec::CentralityFT0C,
                  hyperrec::CentralityFT0M,
                  hyperrec::PtHe3,
                  hyperrec::PhiHe3,
                  hyperrec::EtaHe3,
                  hyperrec::PtPi,
                  hyperrec::PhiPi,
                  hyperrec::EtaPi,
                  hyperrec::XPrimVtx,
                  hyperrec::YPrimVtx,
                  hyperrec::ZPrimVtx,
                  hyperrec::XDecVtx,
                  hyperrec::YDecVtx,
                  hyperrec::ZDecVtx,
                  hyperrec::DcaV0Daug,
                  hyperrec::DcaHe,
                  hyperrec::DcaPi,
                  hyperrec::NSigmaHe,
                  hyperrec::NTPCclusHe,
                  hyperrec::NTPCclusPi,
                  hyperrec::TPCmomHe,
                  hyperrec::TPCmomPi,
                  hyperrec::TPCsignalHe,
                  hyperrec::TPCsignalPi,
                  hyperrec::ITSclusterSizesHe,
                  hyperrec::ITSclusterSizesPi,
                  hyperrec::Flags,
                  hyperrec::GenPt,
                  hyperrec::GenPhi,
                  hyperrec::GenEta,
                  hyperrec::GenPtHe3,
                  hyperrec::GenXDecVtx,
                  hyperrec::GenYDecVtx,
                  hyperrec::GenZDecVtx,
                  hyperrec::IsReco,
                  hyperrec::IsSignal);

using DataHypCand = DataHypCands::iterator;
using MCHypCand = MCHypCands::iterator;

} // namespace o2::aod

#endif // PWGLF_DATAMODEL_LFHYPERNUCLEITABLES_H_
