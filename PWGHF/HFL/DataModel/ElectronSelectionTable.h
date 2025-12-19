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

/// \file ElectronSelectionTable.h
/// \brief Definitions of tables produced by Electron Selection

/// \author Rashi Gupta <rashi.gupta@cern.ch>, IIT Indore
/// \author Ravindra Singh <ravindra.singh@cern.ch>, IIT Indore

#ifndef PWGHF_HFL_DATAMODEL_ELECTRONSELECTIONTABLE_H_
#define PWGHF_HFL_DATAMODEL_ELECTRONSELECTIONTABLE_H_

#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>

#include <cstdint>
#include <vector>

namespace o2::aod
{
// definition of columns and tables for electron selection
namespace hf_sel_electron
{
DECLARE_SOA_INDEX_COLUMN(Collision, collision);                //! collisioniD of the electron track
DECLARE_SOA_INDEX_COLUMN(Track, track);                        //! trackid of of the electron track
DECLARE_SOA_COLUMN(EtaTrack, etaTrack, float);                 //! pseudorapidity of the electron track
DECLARE_SOA_COLUMN(PhiTrack, phiTrack, float);                 //! azimuth of the electron track
DECLARE_SOA_COLUMN(PtTrack, ptTrack, float);                   //! transverse momentum of the electron track
DECLARE_SOA_COLUMN(PTrack, pTrack, float);                     //! momentum of the electron track
DECLARE_SOA_COLUMN(RapidityTrack, rapidityTrack, float);       //! rapidity of the electron track
DECLARE_SOA_COLUMN(DcaXYTrack, dcaXYTrack, float);             //! dca of the electron in xy direction
DECLARE_SOA_COLUMN(DcaZTrack, dcaZTrack, float);               //! dca of the electron in z direction
DECLARE_SOA_COLUMN(TpcNSigmaElTrack, tpcNSigmaElTrack, float); //! tpcNSigma of the electron track(TPC PID)
DECLARE_SOA_COLUMN(TofNSigmaElTrack, tofNSigmaElTrack, float); //! tofNSigma of the electron track(TOF PID)

// EMCal cluster values
DECLARE_SOA_COLUMN(EnergyEmcCluster, energyEmcCluster, float);   //! energy of the EMCal cluster
DECLARE_SOA_COLUMN(EtaEmcCluster, etaEmcCluster, float);         //! pseudorapidity of the EMCal cluster
DECLARE_SOA_COLUMN(PhiEmcCluster, phiEmcCluster, float);         //! azimuth of the EMCal cluster
DECLARE_SOA_COLUMN(M02EmcCluster, m02EmcCluster, float);         //! shower shape long axis of the EMCal cluster
DECLARE_SOA_COLUMN(M20EmcCluster, m20EmcCluster, float);         //! shower shape short axis of the EMCal cluster
DECLARE_SOA_COLUMN(NCellsEmcCluster, nCellsEmcCluster, uint8_t); //! number of cells of the EMCal cluster
DECLARE_SOA_COLUMN(TimeEmcCluster, timeEmcCluster, float);       //! time of the EMCal cluster (ns)

DECLARE_SOA_COLUMN(DeltaEtaMatch, deltaEtaMatch, float); //! dEta matched track to EMCal cluster
DECLARE_SOA_COLUMN(DeltaPhiMatch, deltaPhiMatch, float); //! dPhi matched track to EMCal cluster
DECLARE_SOA_COLUMN(IsEmcal, isEmcal, bool);              //! electron information with Emcal
} // namespace hf_sel_electron
DECLARE_SOA_TABLE(HfSelEl, "AOD", "HFSELEL", //! Electron Informations
                  hf_sel_electron::CollisionId,
                  hf_sel_electron::TrackId,
                  hf_sel_electron::EtaTrack,
                  hf_sel_electron::PhiTrack,
                  hf_sel_electron::PtTrack,
                  hf_sel_electron::PTrack,
                  hf_sel_electron::RapidityTrack,
                  hf_sel_electron::DcaXYTrack,
                  hf_sel_electron::DcaZTrack,
                  hf_sel_electron::TpcNSigmaElTrack,
                  hf_sel_electron::TofNSigmaElTrack,
                  hf_sel_electron::EnergyEmcCluster,
                  hf_sel_electron::EtaEmcCluster,
                  hf_sel_electron::PhiEmcCluster,
                  hf_sel_electron::M02EmcCluster,
                  hf_sel_electron::M20EmcCluster,
                  hf_sel_electron::NCellsEmcCluster,
                  hf_sel_electron::TimeEmcCluster,
                  hf_sel_electron::DeltaEtaMatch,
                  hf_sel_electron::DeltaPhiMatch,
                  hf_sel_electron::IsEmcal);
// definition of columns and tables for HfcorrElectron Selection
namespace hf_corr_sel_electron
{
DECLARE_SOA_INDEX_COLUMN(Collision, collision);                //! collisioniD of the electron track
DECLARE_SOA_INDEX_COLUMN(Track, track);                        //! trackid of of the electron track
DECLARE_SOA_COLUMN(EtaTrack, etaTrack, float);                 //! pseudorapidity of the electron track
DECLARE_SOA_COLUMN(PhiTrack, phiTrack, float);                 //! azimuth of the electron track
DECLARE_SOA_COLUMN(PtTrack, ptTrack, float);                   //! transverse momentum of the electron track
DECLARE_SOA_COLUMN(TpcNSigmaElTrack, tpcNSigmaElTrack, float); //! tpcNSigma of the electron track(TPC PID)
DECLARE_SOA_COLUMN(TofNSigmaElTrack, tofNSigmaElTrack, float); //! tofNSigma of the electron track(TOF PID)
DECLARE_SOA_COLUMN(LSMassEE, lSMassEE, std::vector<float>);    //!  mass of the Like sign electron pair
DECLARE_SOA_COLUMN(ULSMassEE, uLSMassEE, std::vector<float>);  //! mass of UnLike sign electron pair
DECLARE_SOA_COLUMN(NElPairLS, nElPairLS, int);                 //! Number of Like sign electron pair
DECLARE_SOA_COLUMN(NElPairUS, nElPairUS, int);                 //! Number of UnLike sign electron pair
DECLARE_SOA_COLUMN(IsEmcal, isEmcal, bool);                    //! electron information
} // namespace hf_corr_sel_electron

DECLARE_SOA_TABLE(HfCorrSelEl, "AOD", "HfCORRSELEL", //! Electron Informations
                  hf_corr_sel_electron::CollisionId,
                  hf_corr_sel_electron::TrackId,
                  hf_corr_sel_electron::EtaTrack,
                  hf_corr_sel_electron::PhiTrack,
                  hf_corr_sel_electron::PtTrack,
                  hf_corr_sel_electron::TpcNSigmaElTrack,
                  hf_corr_sel_electron::TofNSigmaElTrack,
                  hf_corr_sel_electron::LSMassEE,
                  hf_corr_sel_electron::ULSMassEE,
                  hf_corr_sel_electron::NElPairLS,
                  hf_corr_sel_electron::NElPairUS,
                  hf_corr_sel_electron::IsEmcal);

// definition of columns and tables for Mc Gen HfElectron Selection
namespace hf_mcgen_sel_electron
{
DECLARE_SOA_INDEX_COLUMN(McCollision, mcCollision); //! collisioniD of the electron track
DECLARE_SOA_INDEX_COLUMN(Track, track);             //! trackid of of the electron track
DECLARE_SOA_COLUMN(EtaTrackMc, etaTrackMc, float);  //! pseudorapidity of the electron track
DECLARE_SOA_COLUMN(PhiTrackMc, phiTrackMc, float);  //! azimuth of the electron track
DECLARE_SOA_COLUMN(PtTrackMc, ptTrackMc, float);    //! transverse momentum of the electron track
DECLARE_SOA_COLUMN(IsNonHfeMc, isNonHfeMc, bool);   //! Non-Heavy flavour  electron information

} // namespace hf_mcgen_sel_electron

DECLARE_SOA_TABLE(HfMcGenSelEl, "AOD", "HFMCGENSELEL", //! Electron Informations
                  hf_mcgen_sel_electron::McCollisionId,
                  hf_mcgen_sel_electron::TrackId,
                  hf_mcgen_sel_electron::EtaTrackMc,
                  hf_mcgen_sel_electron::PhiTrackMc,
                  hf_mcgen_sel_electron::PtTrackMc,
                  hf_mcgen_sel_electron::IsNonHfeMc);
} // namespace o2::aod

#endif // PWGHF_HFL_DATAMODEL_ELECTRONSELECTIONTABLE_H_
