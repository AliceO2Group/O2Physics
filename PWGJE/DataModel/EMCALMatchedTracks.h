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

// Table definitions for EMCAL matched tracks analysis
//
/// \author Marvin Hemmer <marvin.hemmer@cern.ch>, Goethe University Franfkurt

#ifndef PWGJE_DATAMODEL_EMCALMATCHEDTRACKS_H_
#define PWGJE_DATAMODEL_EMCALMATCHEDTRACKS_H_

#include <string>
#include "Framework/AnalysisDataModel.h"
#include "EMCALClusterDefinition.h"

namespace o2::aod
{
namespace emcaltrackmatch
{
// Event values
DECLARE_SOA_COLUMN(Orbit, orbit, uint32_t);         //! orbit ID
DECLARE_SOA_COLUMN(Timestamp, timestamp, uint64_t); //! Timestamp of a BC in ms (epoch style)
DECLARE_SOA_COLUMN(RunNumber, runNumber, int);      //! Run number

// track 1 values
DECLARE_SOA_COLUMN(Track1X, track1X, float);                     //! x position of first matched track
DECLARE_SOA_COLUMN(Track1Alpha, track1Alpha, float);             //! alpha of first matched track
DECLARE_SOA_COLUMN(Track1P, track1P, float);                     //! momentum of first matched track
DECLARE_SOA_COLUMN(Track1QPt, track1QPt, float);                 //! q over pT of first matched track
DECLARE_SOA_COLUMN(Track1Y, track1Y, float);                     //! y position of first matched track
DECLARE_SOA_COLUMN(Track1Z, track1Z, float);                     //! z position of first matched track
DECLARE_SOA_COLUMN(Track1Snp, track1Snp, float);                 //! sin(phi) of first matched track
DECLARE_SOA_COLUMN(Track1Tgl, track1Tgl, float);                 //! tan(lambda) of first matched track
DECLARE_SOA_COLUMN(Track1Pt, track1Pt, float);                   //! transverse momentum of first matched track
DECLARE_SOA_COLUMN(Track1SigmaY, track1SigmaY, float);           //! convariance of y position of first matched track
DECLARE_SOA_COLUMN(Track1SigmaZ, track1SigmaZ, float);           //! convariance of z position of first matched track
DECLARE_SOA_COLUMN(Track1SigmaSnp, track1SigmaSnp, float);       //! convariance of sin(phi) of first matched track
DECLARE_SOA_COLUMN(Track1SigmaTgl, track1SigmaTgl, float);       //! convariance of tan(lambda) of first matched track
DECLARE_SOA_COLUMN(Track1SigmaPt, track1SigmaPt, float);         //! convariance of transverse momentum of first matched track
DECLARE_SOA_COLUMN(Track1Eta, track1Eta, float);                 //! eta position of first matched track
DECLARE_SOA_COLUMN(Track1Phi, track1Phi, float);                 //! phi position of first matched track
DECLARE_SOA_COLUMN(Track1EtaEMCAL, track1EtaEmcal, float);       //! eta position of first matched track propagated to calorimeter
DECLARE_SOA_COLUMN(Track1PhiEMCAL, track1PhiEmca, float);        //! phi position of first matched track propagated to calorimeter
DECLARE_SOA_COLUMN(Track1DEta, track1DEta, float);               //! dEta first matched track propagated to calorimeter
DECLARE_SOA_COLUMN(Track1DPhi, track1DPhi, float);               //! dPhi first matched track propagated to calorimeter
DECLARE_SOA_COLUMN(Track1ITSNCls, track1ItsNCls, uint8_t);       //! Number of ITS clusters of first matched track
DECLARE_SOA_COLUMN(Track1TOFExpMom, track1TofExpMom, float);     //! TOF expected momentum obtained in tracking, used to compute the expected times of first matched track
DECLARE_SOA_COLUMN(Track1TPCNSigmaE, track1TPCNSigmaE, float);   //! NSigma electron (TPC PID) of first matched track
DECLARE_SOA_COLUMN(Track1TPCNSigmaPi, track1TPCNSigmaPi, float); //! NSigma pion (TPC PID) of first matched track
DECLARE_SOA_COLUMN(Track1TOFNSigmaE, track1TOFNSigmaE, float);   //! NSigma electron (TOF PID) of first matched track
DECLARE_SOA_COLUMN(Track1TOFNSigmaPi, track1TOFNSigmaPi, float); //! NSigma pion (TOF PID) of first matched track
// track 2 values
DECLARE_SOA_COLUMN(Track2X, track2X, float);                     //! x position of second matched track
DECLARE_SOA_COLUMN(Track2Alpha, track2Alpha, float);             //! alpha of second matched track
DECLARE_SOA_COLUMN(Track2P, track2P, float);                     //! momentum of second matched track
DECLARE_SOA_COLUMN(Track2QPt, track2QPt, float);                 //! q over pT of second matched track
DECLARE_SOA_COLUMN(Track2Y, track2Y, float);                     //! y position of second matched track
DECLARE_SOA_COLUMN(Track2Z, track2Z, float);                     //! z position of second matched track
DECLARE_SOA_COLUMN(Track2Snp, track2Snp, float);                 //! sin(phi) of second matched track
DECLARE_SOA_COLUMN(Track2Tgl, track2Tgl, float);                 //! tan(lambda) of second matched track
DECLARE_SOA_COLUMN(Track2Pt, track2Pt, float);                   //! transverse momentum of second matched track
DECLARE_SOA_COLUMN(Track2SigmaY, track2SigmaY, float);           //! convariance of y position of second matched track
DECLARE_SOA_COLUMN(Track2SigmaZ, track2SigmaZ, float);           //! convariance of z position of second matched track
DECLARE_SOA_COLUMN(Track2SigmaSnp, track2SigmaSnp, float);       //! convariance of sin(phi) of second matched track
DECLARE_SOA_COLUMN(Track2SigmaTgl, track2SigmaTgl, float);       //! convariance of tan(lambda) of second matched track
DECLARE_SOA_COLUMN(Track2SigmaPt, track2SigmaPt, float);         //! convariance of transverse momentum of second matched track
DECLARE_SOA_COLUMN(Track2Eta, track2Eta, float);                 //! eta position of second matched track
DECLARE_SOA_COLUMN(Track2Phi, track2Phi, float);                 //! phi position of second matched track
DECLARE_SOA_COLUMN(Track2EtaEMCAL, track2EtaEmcal, float);       //! eta position of second matched track propagated to calorimeter
DECLARE_SOA_COLUMN(Track2PhiEMCAL, track2PhiEmca, float);        //! phi position of second matched track propagated to calorimeter
DECLARE_SOA_COLUMN(Track2DEta, track2DEta, float);               //! dEta second matched track propagated to calorimeter
DECLARE_SOA_COLUMN(Track2DPhi, track2DPhi, float);               //! dPhi second matched track propagated to calorimeter
DECLARE_SOA_COLUMN(Track2ITSNCls, track2ItsNCls, uint8_t);       //! Number of ITS clusters of second matched track
DECLARE_SOA_COLUMN(Track2TOFExpMom, track2TofExpMom, float);     //! TOF expected momentum obtained in tracking, used to compute the expected times of second matched track
DECLARE_SOA_COLUMN(Track2TPCNSigmaE, track2TPCNSigmaE, float);   //! NSigma electron (TPC PID) of second matched track
DECLARE_SOA_COLUMN(Track2TPCNSigmaPi, track2TPCNSigmaPi, float); //! NSigma pion (TPC PID) of second matched track
DECLARE_SOA_COLUMN(Track2TOFNSigmaE, track2TOFNSigmaE, float);   //! NSigma electron (TOF PID) of second matched track
DECLARE_SOA_COLUMN(Track2TOFNSigmaPi, track2TOFNSigmaPi, float); //! NSigma pion (TOF PID) of second matched track

// cluster values
DECLARE_SOA_COLUMN(ClusterEnergy, clusterE, float);        //! cluster energy (GeV)
DECLARE_SOA_COLUMN(ClusterEta, clusterEta, float);         //! cluster pseudorapidity (calculated using vertex)
DECLARE_SOA_COLUMN(ClusterPhi, clusterPhi, float);         //! cluster azimuthal angle (calculated using vertex)
DECLARE_SOA_COLUMN(ClusterM02, clusterM02, float);         //! cluster shower shape long axis
DECLARE_SOA_COLUMN(ClusterNCells, clusterNCells, uint8_t); //! number of cells in cluster
DECLARE_SOA_COLUMN(ClusterTime, clusterTime, float);       //! cluster time (ns)

} // namespace emcaltrackmatch
// table of clusters that could be matched to a collision
DECLARE_SOA_TABLE(EmcalMTs, "AOD", "EMCALMTS",                                                                                                                                                          //!
                  o2::soa::Index<>, emcaltrackmatch::Orbit, emcaltrackmatch::Timestamp, emcaltrackmatch::RunNumber,                                                                                     // event info
                  emcaltrackmatch::Track1X, emcaltrackmatch::Track1Alpha, emcaltrackmatch::Track1P, emcaltrackmatch::Track1QPt,                                                                         // track 1 X and P
                  emcaltrackmatch::Track1Y, emcaltrackmatch::Track1Z, emcaltrackmatch::Track1Snp, emcaltrackmatch::Track1Tgl, emcaltrackmatch::Track1Pt,                                                // track 1 5 vector
                  emcaltrackmatch::Track1SigmaY, emcaltrackmatch::Track1SigmaZ, emcaltrackmatch::Track1SigmaSnp, emcaltrackmatch::Track1SigmaTgl, emcaltrackmatch::Track1SigmaPt,                       // track 1 5 vector convariance
                  emcaltrackmatch::Track1Eta, emcaltrackmatch::Track1Phi, emcaltrackmatch::Track1EtaEMCAL, emcaltrackmatch::Track1PhiEMCAL, emcaltrackmatch::Track1DEta, emcaltrackmatch::Track1DPhi,   // track 1 eta phi info
                  emcaltrackmatch::Track1ITSNCls, emcaltrackmatch::Track1TOFExpMom,                                                                                                                     // track 1 ITS, TOF info
                  emcaltrackmatch::Track1TPCNSigmaE, emcaltrackmatch::Track1TPCNSigmaPi, emcaltrackmatch::Track1TOFNSigmaE, emcaltrackmatch::Track1TOFNSigmaPi,                                         // track 1 PID info
                  emcaltrackmatch::Track2X, emcaltrackmatch::Track2Alpha, emcaltrackmatch::Track2P, emcaltrackmatch::Track2QPt,                                                                         // track 2 X and P
                  emcaltrackmatch::Track2Y, emcaltrackmatch::Track2Z, emcaltrackmatch::Track2Snp, emcaltrackmatch::Track2Tgl, emcaltrackmatch::Track2Pt,                                                // track 2 5 vector
                  emcaltrackmatch::Track2SigmaY, emcaltrackmatch::Track2SigmaZ, emcaltrackmatch::Track2SigmaSnp, emcaltrackmatch::Track2SigmaTgl, emcaltrackmatch::Track2SigmaPt,                       // track 2 5 vector convariance
                  emcaltrackmatch::Track2Eta, emcaltrackmatch::Track2Phi, emcaltrackmatch::Track2EtaEMCAL, emcaltrackmatch::Track2PhiEMCAL, emcaltrackmatch::Track2DEta, emcaltrackmatch::Track2DPhi,   // track 2 eta phi info
                  emcaltrackmatch::Track2ITSNCls, emcaltrackmatch::Track2TOFExpMom,                                                                                                                     // track 2 ITS, TOF info
                  emcaltrackmatch::Track2TPCNSigmaE, emcaltrackmatch::Track2TPCNSigmaPi, emcaltrackmatch::Track2TOFNSigmaE, emcaltrackmatch::Track2TOFNSigmaPi,                                         // track 2 PID info
                  emcaltrackmatch::ClusterEnergy, emcaltrackmatch::ClusterEta, emcaltrackmatch::ClusterPhi, emcaltrackmatch::ClusterM02, emcaltrackmatch::ClusterNCells, emcaltrackmatch::ClusterTime); // cluster infos

using EmcalMT = EmcalMTs::iterator;

} // namespace o2::aod
#endif // PWGJE_DATAMODEL_EMCALMATCHEDTRACKS_H_
