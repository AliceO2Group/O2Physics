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
/// \file   OTFPIDTrk.h
/// \author Berkin Ulukutlu TUM
/// \author Henrik Fribert TUM
/// \author Nicolò Jacazio Università del Piemonte Orientale
/// \since  May 22, 2025
/// \brief  Set of tables for the ALICE3 Trk PID information
///

#ifndef ALICE3_DATAMODEL_OTFPIDTRK_H_
#define ALICE3_DATAMODEL_OTFPIDTRK_H_

// O2 includes
#include "Framework/AnalysisDataModel.h"

namespace o2::aod
{
namespace upgrade::trk
{

DECLARE_SOA_COLUMN(TimeOverThresholdBarrel, timeOverThresholdBarrel, float);   //! Time over threshold for the barrel layers
DECLARE_SOA_COLUMN(ClusterSizeBarrel, clusterSizeBarrel, float);               //! Cluster size for the barrel layers
DECLARE_SOA_COLUMN(TimeOverThresholdForward, timeOverThresholdForward, float); //! Time over threshold for the Forward layers
DECLARE_SOA_COLUMN(ClusterSizeForward, clusterSizeForward, float);             //! Cluster size for the barrel layers

DECLARE_SOA_COLUMN(NSigmaTrkEl, nSigmaEl, float); //! NSigma electron from the tracker layers
DECLARE_SOA_COLUMN(NSigmaTrkMu, nSigmaMu, float); //! NSigma muon from the tracker layers
DECLARE_SOA_COLUMN(NSigmaTrkPi, nSigmaPi, float); //! NSigma pion from the tracker layers
DECLARE_SOA_COLUMN(NSigmaTrkKa, nSigmaKa, float); //! NSigma kaon from the tracker layers
DECLARE_SOA_COLUMN(NSigmaTrkPr, nSigmaPr, float); //! NSigma proton from the tracker layers

DECLARE_SOA_DYNAMIC_COLUMN(NSigmaTrk, nSigmaTrk, //! General function to get the nSigma for the tracker layers
                           [](const float el,
                              const float mu,
                              const float pi,
                              const float ka,
                              const float pr,
                              const int id) -> float {
                             switch (std::abs(id)) {
                               case 0:
                                 return el;
                               case 1:
                                 return mu;
                               case 2:
                                 return pi;
                               case 3:
                                 return ka;
                               case 4:
                                 return pr;
                               default:
                                 LOG(fatal) << "Unrecognized PDG code for InnerTOF";
                                 return 999.f;
                             }
                           });

} // namespace upgrade::trk

DECLARE_SOA_TABLE(UpgradeTrkPidSignals, "AOD", "UPGRADETRKSIG",
                  upgrade::trk::TimeOverThresholdBarrel,
                  upgrade::trk::ClusterSizeBarrel);

DECLARE_SOA_TABLE(UpgradeTrkPids, "AOD", "UPGRADETRKPID",
                  upgrade::trk::NSigmaTrkEl,
                  upgrade::trk::NSigmaTrkMu,
                  upgrade::trk::NSigmaTrkPi,
                  upgrade::trk::NSigmaTrkKa,
                  upgrade::trk::NSigmaTrkPr,
                  upgrade::trk::NSigmaTrk<upgrade::trk::NSigmaTrkEl,
                                          upgrade::trk::NSigmaTrkMu,
                                          upgrade::trk::NSigmaTrkPi,
                                          upgrade::trk::NSigmaTrkKa,
                                          upgrade::trk::NSigmaTrkPr>);

using UpgradeTrkPidSignal = UpgradeTrkPidSignals::iterator;
using UpgradeTrkPid = UpgradeTrkPids::iterator;

} // namespace o2::aod

#endif // ALICE3_DATAMODEL_OTFPIDTRK_H_
