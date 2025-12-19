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

DECLARE_SOA_COLUMN(TimeOverThresholdBarrel, timeOverThresholdBarrel, float);   //! Time over threshold for the Barrel layers
DECLARE_SOA_COLUMN(TimeOverThresholdForward, timeOverThresholdForward, float); //! Time over threshold for the Forward layers

DECLARE_SOA_COLUMN(NSigmaTrkEl, nSigmaTrkEl, float); //! NSigma electron from the tracker layers
DECLARE_SOA_COLUMN(NSigmaTrkMu, nSigmaTrkMu, float); //! NSigma muon from the tracker layers
DECLARE_SOA_COLUMN(NSigmaTrkPi, nSigmaTrkPi, float); //! NSigma pion from the tracker layers
DECLARE_SOA_COLUMN(NSigmaTrkKa, nSigmaTrkKa, float); //! NSigma kaon from the tracker layers
DECLARE_SOA_COLUMN(NSigmaTrkPr, nSigmaTrkPr, float); //! NSigma proton from the tracker layers
DECLARE_SOA_COLUMN(NSigmaTrkDe, nSigmaTrkDe, float); //! NSigma deuteron from the tracker layers
DECLARE_SOA_COLUMN(NSigmaTrkTr, nSigmaTrkTr, float); //! NSigma triton from the tracker layers
DECLARE_SOA_COLUMN(NSigmaTrkHe, nSigmaTrkHe, float); //! NSigma helium-3 from the tracker layers
DECLARE_SOA_COLUMN(NSigmaTrkAl, nSigmaTrkAl, float); //! NSigma alpha from the tracker layers

DECLARE_SOA_DYNAMIC_COLUMN(NSigmaTrk, nSigmaTrk, //! General function to get the nSigma for the tracker layers
                           [](const float el,
                              const float mu,
                              const float pi,
                              const float ka,
                              const float pr,
                              const float de,
                              const float tr,
                              const float he,
                              const float al,
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
                               case 5:
                                 return de;
                               case 6:
                                 return tr;
                               case 7:
                                 return he;
                               case 8:
                                 return al;
                               default:
                                 LOG(fatal) << "Unrecognized PDG code";
                                 return 999.f;
                             }
                           });

} // namespace upgrade::trk

DECLARE_SOA_TABLE(UpgradeTrkPidSignals, "AOD", "UPGRADETRKSIG",
                  o2::soa::Index<>,
                  upgrade::trk::TimeOverThresholdBarrel);

DECLARE_SOA_TABLE(UpgradeTrkPids, "AOD", "UPGRADETRKPID",
                  o2::soa::Index<>,
                  upgrade::trk::NSigmaTrkEl,
                  upgrade::trk::NSigmaTrkMu,
                  upgrade::trk::NSigmaTrkPi,
                  upgrade::trk::NSigmaTrkKa,
                  upgrade::trk::NSigmaTrkPr,
                  upgrade::trk::NSigmaTrkDe,
                  upgrade::trk::NSigmaTrkTr,
                  upgrade::trk::NSigmaTrkHe,
                  upgrade::trk::NSigmaTrkAl,
                  upgrade::trk::NSigmaTrk<upgrade::trk::NSigmaTrkEl,
                                          upgrade::trk::NSigmaTrkMu,
                                          upgrade::trk::NSigmaTrkPi,
                                          upgrade::trk::NSigmaTrkKa,
                                          upgrade::trk::NSigmaTrkPr,
                                          upgrade::trk::NSigmaTrkDe,
                                          upgrade::trk::NSigmaTrkTr,
                                          upgrade::trk::NSigmaTrkHe,
                                          upgrade::trk::NSigmaTrkAl>);

using UpgradeTrkPidSignal = UpgradeTrkPidSignals::iterator;
using UpgradeTrkPid = UpgradeTrkPids::iterator;

} // namespace o2::aod

#endif // ALICE3_DATAMODEL_OTFPIDTRK_H_
