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

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"

#ifndef PWGLF_DATAMODEL_LFHYPERNUCLEITABLES_H_
#define PWGLF_DATAMODEL_LFHYPERNUCLEITABLES_H_

namespace o2::aod
{
namespace hyperrec
{
DECLARE_SOA_COLUMN(CentralityFT0A, centralityFT0A, float); // centrality with FT0A estimator
DECLARE_SOA_COLUMN(CentralityFT0C, centralityFT0C, float); // centrality with FT0C estimator
DECLARE_SOA_COLUMN(CentralityFT0M, centralityFT0M, float); // centrality with FT0M estimator
DECLARE_SOA_COLUMN(PsiFT0A, psiFT0A, float);               // Psi with FT0A estimator
DECLARE_SOA_COLUMN(MultFT0A, multFT0A, float);             // Multiplicity with FT0A estimator
DECLARE_SOA_COLUMN(PsiFT0C, psiFT0C, float);               // Psi with FT0C estimator
DECLARE_SOA_COLUMN(QFT0C, qFT0C, float);                   // Amplitude with FT0C estimator
DECLARE_SOA_COLUMN(MultFT0C, multFT0C, float);             // Multiplicity with FT0C estimator
DECLARE_SOA_COLUMN(PsiTPC, psiTPC, float);                 // Psi with TPC estimator
DECLARE_SOA_COLUMN(MultTPC, multTPC, float);               // Multiplicity with TPC estimator
DECLARE_SOA_COLUMN(CollisionId, collisionId, int64_t);     // CollisionID

DECLARE_SOA_COLUMN(RunNumber, runNumber, int32_t);                    // Run number
DECLARE_SOA_COLUMN(IsMatter, isMatter, bool);                         // bool: true for matter
DECLARE_SOA_COLUMN(PtHe3, ptHe3, float);                              // Pt of the He daughter
DECLARE_SOA_COLUMN(PhiHe3, phiHe3, float);                            // Phi of the He daughter
DECLARE_SOA_COLUMN(EtaHe3, etaHe3, float);                            // Eta of the He daughter
DECLARE_SOA_COLUMN(PtPi, ptPi, float);                                // Pt of the Pi daughter
DECLARE_SOA_COLUMN(PhiPi, phiPi, float);                              // Phi of the Pi daughter
DECLARE_SOA_COLUMN(EtaPi, etaPi, float);                              // Eta of the Pi daughter
DECLARE_SOA_COLUMN(XPrimVtx, xPrimVtx, float);                        // Decay vertex of the candidate (x direction)
DECLARE_SOA_COLUMN(YPrimVtx, yPrimVtx, float);                        // Decay vertex of the candidate (y direction)
DECLARE_SOA_COLUMN(ZPrimVtx, zPrimVtx, float);                        // Decay vertex of the candidate (z direction)
DECLARE_SOA_COLUMN(XDecVtx, xDecVtx, float);                          // Decay vertex of the candidate (x direction)
DECLARE_SOA_COLUMN(YDecVtx, yDecVtx, float);                          // Decay vertex of the candidate (y direction)
DECLARE_SOA_COLUMN(ZDecVtx, zDecVtx, float);                          // Decay vertex of the candidate (z direction)
DECLARE_SOA_COLUMN(MassH3L, massH3L, float);                          // Squared mass w/ hypertriton mass hypo
DECLARE_SOA_COLUMN(MassH4L, massH4L, float);                          // Squared mass w/ H4L mass hypo
DECLARE_SOA_COLUMN(DcaV0Daug, dcaV0Daug, float);                      // DCA between daughters
DECLARE_SOA_COLUMN(CosPA, cosPA, double);                             // Cosine of the pointing angle
DECLARE_SOA_COLUMN(NSigmaHe, nSigmaHe, float);                        // Number of sigmas of the He daughter
DECLARE_SOA_COLUMN(NTPCclusHe, nTPCclusHe, uint8_t);                  // Number of TPC clusters of the He daughter
DECLARE_SOA_COLUMN(NTPCclusPi, nTPCclusPi, uint8_t);                  // Number of TPC clusters of the Pi daughter
DECLARE_SOA_COLUMN(NTPCpidClusHe, nTPCpidClusHe, uint8_t);            // Number of TPC clusters with PID information of the He daughter
DECLARE_SOA_COLUMN(NTPCpidClusPi, nTPCpidClusPi, uint8_t);            // Number of TPC clusters with PID information of the Pi daughter
DECLARE_SOA_COLUMN(TPCsignalHe, tpcSignalHe, uint16_t);               // TPC signal of the He daughter
DECLARE_SOA_COLUMN(TPCsignalPi, tpcSignalPi, uint16_t);               // TPC signal of the Pi daughter
DECLARE_SOA_COLUMN(TPCChi2He, tpcChi2He, float);                      // TPC chi2 of the He daughter
DECLARE_SOA_COLUMN(ITSChi2He, itsChi2He, float);                      // ITS chi2 of the He daughter
DECLARE_SOA_COLUMN(ITSChi2Pi, itsChi2Pi, float);                      // ITS chi2 of the Pi daughter
DECLARE_SOA_COLUMN(TrackedClSize, trackedClSize, int);                // int: zero for non-tracked candidates
DECLARE_SOA_COLUMN(Flags, flags, uint8_t);                            // Flags for PID in tracking (bits [0, 3] for negative daughter, [4,7] for positive daughter)
DECLARE_SOA_COLUMN(TPCmomHe, tpcMomHe, float);                        // TPC momentum of the He daughter
DECLARE_SOA_COLUMN(TPCmomPi, tpcMomPi, float);                        // TPC momentum of the Pi daughter
DECLARE_SOA_COLUMN(TOFMass, tofMass, float);                          // TOF mass of the candidate
DECLARE_SOA_COLUMN(ITSclusterSizesHe, itsClusterSizesHe, uint32_t);   // ITS cluster size of the He daughter
DECLARE_SOA_COLUMN(ITSclusterSizesPi, itsClusterSizesPi, uint32_t);   // ITS cluster size of the Pi daughter
DECLARE_SOA_COLUMN(ITSclusterSizesHyp, itsClusterSizesHyp, uint32_t); // ITS cluster size of the Pi daughter
DECLARE_SOA_COLUMN(DcaHe, dcaHe, float);                              // DCA between He daughter and V0
DECLARE_SOA_COLUMN(DcaPi, dcaPi, float);                              // DCA between pi daughter and V0
DECLARE_SOA_COLUMN(GenPt, genPt, float);                              // Pt of the hypertriton
DECLARE_SOA_COLUMN(GenPhi, genPhi, float);                            // Phi of the hypertriton
DECLARE_SOA_COLUMN(GenEta, genEta, float);                            // Eta of the hypertriton
DECLARE_SOA_COLUMN(GenPtHe3, genPtHe3, float);                        // Pt of the He daughter (to be used for the recalibration)
DECLARE_SOA_COLUMN(GenXDecVtx, genXDecVtx, float);                    // Decay vertex of the candidate (x direction)
DECLARE_SOA_COLUMN(GenYDecVtx, genYDecVtx, float);                    // Decay vertex of the candidate (y direction)
DECLARE_SOA_COLUMN(GenZDecVtx, genZDecVtx, float);                    // Decay vertex of the candidate (z direction)
DECLARE_SOA_COLUMN(IsReco, isReco, bool);                             // bool: true for reco
DECLARE_SOA_COLUMN(IsFakeHeOnITSLayer, isFakeHeOnITSLayer, uint8_t);  // uint8_t: bit map for fake He on ITS layers
DECLARE_SOA_COLUMN(IsSignal, isSignal, bool);                         // bool: true for signal
DECLARE_SOA_COLUMN(IsRecoMCCollision, isRecoMCCollision, bool);       // bool: true for reco MC collision
DECLARE_SOA_COLUMN(IsSurvEvSel, isSurvEvSel, bool);                   // bool: true for survived event selection
DECLARE_SOA_COLUMN(IsTwoBodyDecay, isTwoBodyDecay, bool);             // bool: true for two body decay
} // namespace hyperrec

DECLARE_SOA_TABLE(DataHypCands, "AOD", "HYPCANDS",
                  o2::soa::Index<>,
                  hyperrec::CentralityFT0A, hyperrec::CentralityFT0C, hyperrec::CentralityFT0M,
                  hyperrec::XPrimVtx, hyperrec::YPrimVtx, hyperrec::ZPrimVtx,

                  hyperrec::RunNumber, hyperrec::IsMatter,
                  hyperrec::PtHe3, hyperrec::PhiHe3, hyperrec::EtaHe3,
                  hyperrec::PtPi, hyperrec::PhiPi, hyperrec::EtaPi,
                  hyperrec::XDecVtx, hyperrec::YDecVtx, hyperrec::ZDecVtx,
                  hyperrec::DcaV0Daug, hyperrec::DcaHe, hyperrec::DcaPi,
                  hyperrec::NSigmaHe, hyperrec::NTPCclusHe, hyperrec::NTPCclusPi, hyperrec::NTPCpidClusHe, hyperrec::NTPCpidClusPi,
                  hyperrec::TPCmomHe, hyperrec::TPCmomPi, hyperrec::TPCsignalHe, hyperrec::TPCsignalPi, hyperrec::TPCChi2He, hyperrec::ITSChi2He, hyperrec::ITSChi2Pi,
                  hyperrec::TOFMass,
                  hyperrec::ITSclusterSizesHe, hyperrec::ITSclusterSizesPi,
                  hyperrec::Flags, hyperrec::TrackedClSize);

DECLARE_SOA_TABLE(DataHypCandsFlow, "AOD", "HYPCANDSFLOW",
                  o2::soa::Index<>,
                  hyperrec::CentralityFT0A, hyperrec::CentralityFT0C, hyperrec::CentralityFT0M,
                  hyperrec::PsiFT0A, hyperrec::MultFT0A,
                  hyperrec::PsiFT0C, hyperrec::MultFT0C, hyperrec::QFT0C,
                  hyperrec::PsiTPC, hyperrec::MultTPC,
                  hyperrec::XPrimVtx, hyperrec::YPrimVtx, hyperrec::ZPrimVtx,

                  hyperrec::RunNumber, hyperrec::IsMatter,
                  hyperrec::PtHe3, hyperrec::PhiHe3, hyperrec::EtaHe3,
                  hyperrec::PtPi, hyperrec::PhiPi, hyperrec::EtaPi,
                  hyperrec::XDecVtx, hyperrec::YDecVtx, hyperrec::ZDecVtx,
                  hyperrec::DcaV0Daug, hyperrec::DcaHe, hyperrec::DcaPi,
                  hyperrec::NSigmaHe, hyperrec::NTPCclusHe, hyperrec::NTPCclusPi, hyperrec::NTPCpidClusHe, hyperrec::NTPCpidClusPi,
                  hyperrec::TPCmomHe, hyperrec::TPCmomPi, hyperrec::TPCsignalHe, hyperrec::TPCsignalPi, hyperrec::TPCChi2He, hyperrec::ITSChi2He, hyperrec::ITSChi2Pi,
                  hyperrec::TOFMass,
                  hyperrec::ITSclusterSizesHe, hyperrec::ITSclusterSizesPi,
                  hyperrec::Flags, hyperrec::TrackedClSize);

DECLARE_SOA_TABLE(MCHypCands, "AOD", "MCHYPCANDS",
                  o2::soa::Index<>,
                  hyperrec::CentralityFT0A, hyperrec::CentralityFT0C, hyperrec::CentralityFT0M,
                  hyperrec::XPrimVtx, hyperrec::YPrimVtx, hyperrec::ZPrimVtx,

                  hyperrec::RunNumber, hyperrec::IsMatter,
                  hyperrec::PtHe3, hyperrec::PhiHe3, hyperrec::EtaHe3,
                  hyperrec::PtPi, hyperrec::PhiPi, hyperrec::EtaPi,
                  hyperrec::XDecVtx, hyperrec::YDecVtx, hyperrec::ZDecVtx,
                  hyperrec::DcaV0Daug, hyperrec::DcaHe, hyperrec::DcaPi,
                  hyperrec::NSigmaHe, hyperrec::NTPCclusHe, hyperrec::NTPCclusPi, hyperrec::NTPCpidClusHe, hyperrec::NTPCpidClusPi,
                  hyperrec::TPCmomHe, hyperrec::TPCmomPi, hyperrec::TPCsignalHe, hyperrec::TPCsignalPi, hyperrec::TPCChi2He, hyperrec::ITSChi2He, hyperrec::ITSChi2Pi,
                  hyperrec::TOFMass,
                  hyperrec::ITSclusterSizesHe, hyperrec::ITSclusterSizesPi,
                  hyperrec::Flags, hyperrec::TrackedClSize,
                  hyperrec::GenPt,
                  hyperrec::GenPhi,
                  hyperrec::GenEta,
                  hyperrec::GenPtHe3,
                  hyperrec::GenXDecVtx,
                  hyperrec::GenYDecVtx,
                  hyperrec::GenZDecVtx,
                  hyperrec::IsReco,
                  hyperrec::IsFakeHeOnITSLayer,
                  hyperrec::IsSignal,
                  hyperrec::IsRecoMCCollision,
                  hyperrec::IsSurvEvSel,
                  hyperrec::IsTwoBodyDecay, aod::mcparticle::StatusCode);

DECLARE_SOA_TABLE(DataHypCandsWColl, "AOD", "HYPCANDSWCOLL",
                  o2::soa::Index<>,
                  hyperrec::CollisionId, hyperrec::CentralityFT0A, hyperrec::CentralityFT0C, hyperrec::CentralityFT0M,
                  hyperrec::XPrimVtx, hyperrec::YPrimVtx, hyperrec::ZPrimVtx,

                  hyperrec::RunNumber, hyperrec::IsMatter,
                  hyperrec::PtHe3, hyperrec::PhiHe3, hyperrec::EtaHe3,
                  hyperrec::PtPi, hyperrec::PhiPi, hyperrec::EtaPi,
                  hyperrec::XDecVtx, hyperrec::YDecVtx, hyperrec::ZDecVtx,
                  hyperrec::DcaV0Daug, hyperrec::DcaHe, hyperrec::DcaPi,
                  hyperrec::NSigmaHe, hyperrec::NTPCclusHe, hyperrec::NTPCclusPi, hyperrec::NTPCpidClusHe, hyperrec::NTPCpidClusPi,
                  hyperrec::TPCmomHe, hyperrec::TPCmomPi, hyperrec::TPCsignalHe, hyperrec::TPCsignalPi, hyperrec::TPCChi2He, hyperrec::ITSChi2He, hyperrec::ITSChi2Pi,
                  hyperrec::TOFMass,
                  hyperrec::ITSclusterSizesHe, hyperrec::ITSclusterSizesPi,
                  hyperrec::Flags, hyperrec::TrackedClSize);

using DataHypCand = DataHypCands::iterator;
using DataHypCandFlow = DataHypCandsFlow::iterator;
using MCHypCand = MCHypCands::iterator;
using DataHypCandWColl = DataHypCandsWColl::iterator;

namespace hyperkink
{
DECLARE_SOA_COLUMN(PtHyper, ptHyper, float);                              // Pt of the hypertriton
DECLARE_SOA_COLUMN(PhiHyper, phiHyper, float);                            // Phi of the hypertriton
DECLARE_SOA_COLUMN(EtaHyper, etaHyper, float);                            // Eta of the hypertriton
DECLARE_SOA_COLUMN(PtTrit, ptTrit, float);                                // Pt of the triton kink
DECLARE_SOA_COLUMN(PhiTrit, phiTrit, float);                              // Phi of the triton kink
DECLARE_SOA_COLUMN(EtaTrit, etaTrit, float);                              // Eta of the triton kink
DECLARE_SOA_COLUMN(DcaHyperPv, dcaHyperPv, float);                        // DCA of the hypertriton to the primary vertex
DECLARE_SOA_COLUMN(DcaTritPv, dcaTritPv, float);                          // DCA of the triton kink to the primary vertex
DECLARE_SOA_COLUMN(DCAKinkTopo, dcaKinkTopo, float);                      // DCA of the kink topology
DECLARE_SOA_COLUMN(ITSclusterSizesHyper, itsClusterSizesHyper, uint32_t); // ITS cluster size of the hypertriton
DECLARE_SOA_COLUMN(ITSclusterSizesTrit, itsClusterSizesTrit, uint32_t);   // ITS cluster size of the triton kink
DECLARE_SOA_COLUMN(PIDinTrackTrit, pidInTrackTrit, uint8_t);              // PID in track for the triton kink

DECLARE_SOA_COLUMN(TPCmomTrit, tpcMomTrit, float);          // TPC momentum of the triton kink
DECLARE_SOA_COLUMN(TPCsignalTrit, tpcSignalTrit, uint16_t); // TPC signal of the triton kink
DECLARE_SOA_COLUMN(NSigmaTPCTrit, nSigmaTPCTrit, float);    // Number of tpc sigmas of the triton kink
DECLARE_SOA_COLUMN(NSigmaTOFTrit, nSigmaTOFTrit, float);    // Number of tof sigmas of the triton kink

// MC additional info
DECLARE_SOA_COLUMN(GenPtTrit, genPtTrit, float);   // Pt of the triton kink
DECLARE_SOA_COLUMN(HyperPtITS, hyperPtITS, float); // Pt of the hypertriton from ITS standalone, hypertriton tagged with MC truth
DECLARE_SOA_COLUMN(MCMask, mcMask, bool);          // bool: true for fake triton

} // namespace hyperkink

DECLARE_SOA_TABLE(DataHypKinkCands, "AOD", "HYPKINKCANDS",
                  o2::soa::Index<>,
                  hyperrec::XPrimVtx, hyperrec::YPrimVtx, hyperrec::ZPrimVtx,
                  hyperrec::XDecVtx, hyperrec::YDecVtx, hyperrec::ZDecVtx,
                  hyperrec::IsMatter, hyperkink::PtHyper, hyperkink::PhiHyper, hyperkink::EtaHyper,
                  hyperkink::PtTrit, hyperkink::PhiTrit, hyperkink::EtaTrit,
                  hyperkink::DcaHyperPv, hyperkink::DcaTritPv, hyperkink::DCAKinkTopo,
                  hyperkink::ITSclusterSizesHyper, hyperkink::ITSclusterSizesTrit, hyperkink::PIDinTrackTrit,
                  hyperkink::TPCmomTrit, hyperkink::TPCsignalTrit, hyperkink::NSigmaTPCTrit, hyperkink::NSigmaTOFTrit);

DECLARE_SOA_TABLE(MCHypKinkCands, "AOD", "MCHYPKINKCANDS",
                  o2::soa::Index<>,
                  hyperrec::XPrimVtx, hyperrec::YPrimVtx, hyperrec::ZPrimVtx,
                  hyperrec::XDecVtx, hyperrec::YDecVtx, hyperrec::ZDecVtx,
                  hyperrec::IsMatter, hyperkink::PtHyper, hyperkink::PhiHyper, hyperkink::EtaHyper,
                  hyperkink::PtTrit, hyperkink::PhiTrit, hyperkink::EtaTrit,
                  hyperkink::DcaHyperPv, hyperkink::DcaTritPv, hyperkink::DCAKinkTopo,
                  hyperkink::ITSclusterSizesHyper, hyperkink::ITSclusterSizesTrit, hyperkink::PIDinTrackTrit,
                  hyperkink::TPCmomTrit, hyperkink::TPCsignalTrit, hyperkink::NSigmaTPCTrit, hyperkink::NSigmaTOFTrit,
                  hyperrec::GenXDecVtx, hyperrec::GenYDecVtx, hyperrec::GenZDecVtx,
                  hyperrec::GenPt, hyperkink::GenPtTrit,
                  hyperrec::IsReco, hyperrec::IsSignal, hyperkink::MCMask, hyperkink::HyperPtITS,
                  hyperrec::IsRecoMCCollision, hyperrec::IsSurvEvSel);

} // namespace o2::aod

#endif // PWGLF_DATAMODEL_LFHYPERNUCLEITABLES_H_
