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
/// \file LFLnnTables.h
/// \brief Slim lnn tables
///

#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"

#ifndef PWGLF_DATAMODEL_LFLNNTABLES_H_
#define PWGLF_DATAMODEL_LFLNNTABLES_H_

namespace o2::aod
{
namespace lnnrec
{
DECLARE_SOA_COLUMN(CentralityFT0A, centralityFT0A, float); // centrality with FT0A estimator
DECLARE_SOA_COLUMN(CentralityFT0C, centralityFT0C, float); // centrality with FT0C estimator
DECLARE_SOA_COLUMN(CentralityFT0M, centralityFT0M, float); // centrality with FT0M estimator
DECLARE_SOA_COLUMN(PsiFT0A, psiFT0A, float);               // Psi with FT0A estimator
DECLARE_SOA_COLUMN(MultFT0A, multFT0A, float);             // Multiplicity with FT0A estimator
DECLARE_SOA_COLUMN(PsiFT0C, psiFT0C, float);               // Psi with FT0C estimator
DECLARE_SOA_COLUMN(MultFT0C, multFT0C, float);             // Multiplicity with FT0C estimator
DECLARE_SOA_COLUMN(PsiTPC, psiTPC, float);                 // Psi with TPC estimator
DECLARE_SOA_COLUMN(MultTPC, multTPC, float);               // Multiplicity with TPC estimator

DECLARE_SOA_COLUMN(IsMatter, isMatter, bool);                             // bool: true for matter
DECLARE_SOA_COLUMN(Pt3H, pt3H, float);                                    // Pt of the 3H daughter
DECLARE_SOA_COLUMN(Phi3H, phi3H, float);                                  // Phi of the 3H daughter
DECLARE_SOA_COLUMN(Eta3H, eta3H, float);                                  // Eta of the 3H daughter
DECLARE_SOA_COLUMN(PtPi, ptPi, float);                                    // Pt of the Pi daughter
DECLARE_SOA_COLUMN(PhiPi, phiPi, float);                                  // Phi of the Pi daughter
DECLARE_SOA_COLUMN(EtaPi, etaPi, float);                                  // Eta of the Pi daughter
DECLARE_SOA_COLUMN(XPrimVtx, xPrimVtx, float);                            // Decay vertex of the candidate (x direction)
DECLARE_SOA_COLUMN(YPrimVtx, yPrimVtx, float);                            // Decay vertex of the candidate (y direction)
DECLARE_SOA_COLUMN(ZPrimVtx, zPrimVtx, float);                            // Decay vertex of the candidate (z direction)
DECLARE_SOA_COLUMN(XDecVtx, xDecVtx, float);                              // Decay vertex of the candidate (x direction)
DECLARE_SOA_COLUMN(YDecVtx, yDecVtx, float);                              // Decay vertex of the candidate (y direction)
DECLARE_SOA_COLUMN(ZDecVtx, zDecVtx, float);                              // Decay vertex of the candidate (z direction)
DECLARE_SOA_COLUMN(MassLNNL, massLNNL, float);                            // Squared mass w/ lnn mass hypotenuse
DECLARE_SOA_COLUMN(DcaV0Daug, dcaV0Daug, float);                          // DCA between daughters
DECLARE_SOA_COLUMN(CosPA, cosPA, double);                                 // Cosine of the pointing angle
DECLARE_SOA_COLUMN(NSigma3H, nSigma3H, float);                            // Number of sigmas of the 3H daughter
DECLARE_SOA_COLUMN(NTPCclus3H, nTPCclus3H, uint8_t);                      // Number of TPC clusters of the 3H daughter
DECLARE_SOA_COLUMN(NTPCclusPi, nTPCclusPi, uint8_t);                      // Number of TPC clusters of the Pi daughter
DECLARE_SOA_COLUMN(TPCsignal3H, tpcSignal3H, uint16_t);                   // TPC signal of the 3H daughter
DECLARE_SOA_COLUMN(TPCsignalPi, tpcSignalPi, uint16_t);                   // TPC signal of the Pi daughter
DECLARE_SOA_COLUMN(Flags, flags, uint8_t);                                // Flags for PID in tracking (bits [0, 3] for negative daughter, [4,7] for positive daughter)
DECLARE_SOA_COLUMN(TPCmom3H, tpcMom3H, float);                            // TPC momentum of the 3H daughter
DECLARE_SOA_COLUMN(TPCmomPi, tpcMomPi, float);                            // TPC momentum of the Pi daughter
DECLARE_SOA_COLUMN(MassTrTOF, mass2TrTOF, float);                         // TOF 3H mass
DECLARE_SOA_COLUMN(TPCchi3H, tpcChi3H, float);                            // tpcChi3H
DECLARE_SOA_COLUMN(ITSclusterSizes3H, itsClusterSizes3H, uint32_t);       // ITS cluster size of the 3H daughter
DECLARE_SOA_COLUMN(ITSclusterSizesPi, itsClusterSizesPi, uint32_t);       // ITS cluster size of the Pi daughter
DECLARE_SOA_COLUMN(Dca3H, dca3H, float);                                  // DCA between 3H daughter and V0
DECLARE_SOA_COLUMN(DcaPi, dcaPi, float);                                  // DCA between pi daughter and V0
DECLARE_SOA_COLUMN(GenPt, genPt, float);                                  // Pt of the lnn
DECLARE_SOA_COLUMN(GenPhi, genPhi, float);                                // Phi of the lnn
DECLARE_SOA_COLUMN(GenEta, genEta, float);                                // Eta of the lnn
DECLARE_SOA_COLUMN(GenPt3H, genPt3H, float);                              // Pt of the 3H daughter (to be used for the recalibration)
DECLARE_SOA_COLUMN(GenXDecVtx, genXDecVtx, float);                        // Decay vertex of the candidate (x direction)
DECLARE_SOA_COLUMN(GenYDecVtx, genYDecVtx, float);                        // Decay vertex of the candidate (y direction)
DECLARE_SOA_COLUMN(GenZDecVtx, genZDecVtx, float);                        // Decay vertex of the candidate (z direction)
DECLARE_SOA_COLUMN(IsReco, isReco, bool);                                 // bool: true for reco
DECLARE_SOA_COLUMN(IsSignal, isSignal, bool);                             // bool: true for signal
DECLARE_SOA_COLUMN(IsRecoMcCollision, isRecoMcCollision, bool); // bool: true for survived event selection
DECLARE_SOA_COLUMN(SurvivedEventSelection, survivedEventSelection, bool); // bool: true for survived event selection
} // namespace lnnrec

// Declaration of the table DataLnnCands which contains information of TPC, ITS and geometric variables
DECLARE_SOA_TABLE(DataLnnCands, "AOD", "LNNCANDS",
                  o2::soa::Index<>,
                  lnnrec::CentralityFT0A, lnnrec::CentralityFT0C, lnnrec::CentralityFT0M,
                  lnnrec::XPrimVtx, lnnrec::YPrimVtx, lnnrec::ZPrimVtx,

                  lnnrec::IsMatter,
                  lnnrec::Pt3H, lnnrec::Phi3H, lnnrec::Eta3H,
                  lnnrec::PtPi, lnnrec::PhiPi, lnnrec::EtaPi,
                  lnnrec::XDecVtx, lnnrec::YDecVtx, lnnrec::ZDecVtx,
                  lnnrec::DcaV0Daug, lnnrec::Dca3H, lnnrec::DcaPi,
                  lnnrec::NSigma3H, lnnrec::NTPCclus3H, lnnrec::NTPCclusPi,
                  lnnrec::TPCmom3H, lnnrec::TPCmomPi, lnnrec::TPCsignal3H, lnnrec::TPCsignalPi,
                  lnnrec::MassTrTOF, lnnrec::TPCchi3H,
                  lnnrec::ITSclusterSizes3H, lnnrec::ITSclusterSizesPi,
                  lnnrec::Flags);

DECLARE_SOA_TABLE(MCLnnCands, "AOD", "MCLNNCANDS",
                  o2::soa::Index<>,
                  lnnrec::CentralityFT0A, lnnrec::CentralityFT0C, lnnrec::CentralityFT0M,
                  lnnrec::XPrimVtx, lnnrec::YPrimVtx, lnnrec::ZPrimVtx,

                  lnnrec::IsMatter,
                  lnnrec::Pt3H, lnnrec::Phi3H, lnnrec::Eta3H,
                  lnnrec::PtPi, lnnrec::PhiPi, lnnrec::EtaPi,
                  lnnrec::XDecVtx, lnnrec::YDecVtx, lnnrec::ZDecVtx,
                  lnnrec::DcaV0Daug, lnnrec::Dca3H, lnnrec::DcaPi,
                  lnnrec::NSigma3H, lnnrec::NTPCclus3H, lnnrec::NTPCclusPi,
                  lnnrec::TPCmom3H, lnnrec::TPCmomPi, lnnrec::TPCsignal3H, lnnrec::TPCsignalPi,
                  lnnrec::MassTrTOF, lnnrec::TPCchi3H,
                  lnnrec::ITSclusterSizes3H, lnnrec::ITSclusterSizesPi,
                  lnnrec::Flags,
                  lnnrec::GenPt,
                  lnnrec::GenPhi,
                  lnnrec::GenEta,
                  lnnrec::GenPt3H,
                  lnnrec::GenXDecVtx,
                  lnnrec::GenYDecVtx,
                  lnnrec::GenZDecVtx,
                  lnnrec::IsReco,
                  lnnrec::IsSignal, lnnrec::IsRecoMcCollision,
                  lnnrec::SurvivedEventSelection);

using DataLnnCand = DataLnnCands::iterator;
using MCLnnCand = MCLnnCands::iterator;
} // namespace o2::aod

#endif // PWGLF_DATAMODEL_LFLNNTABLES_H_
