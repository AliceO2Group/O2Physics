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
/// \file   RICH.h
/// \author Nicolo' Jacazio
/// \since  25/02/2021
/// \brief  Set of tables for the ALICE3 RICH information
///

#ifndef O2_ANALYSIS_ALICE3_RICH_H_
#define O2_ANALYSIS_ALICE3_RICH_H_

// O2 includes
#include "Framework/AnalysisDataModel.h"
#include "ReconstructionDataFormats/PID.h"

namespace o2::aod
{
namespace alice3rich
{
DECLARE_SOA_INDEX_COLUMN(Track, track);                      //! Index to travel from track to RICH
DECLARE_SOA_COLUMN(RICHSignal, richSignal, float);           //! Signal in RICH
DECLARE_SOA_COLUMN(RICHSignalError, richSignalError, float); //! Error on the RICH signal
DECLARE_SOA_COLUMN(RICHDeltaEl, richDeltaEl, float);         //! signal - exp. signal for electrons
DECLARE_SOA_COLUMN(RICHDeltaMu, richDeltaMu, float);         //! signal - exp. signal for muons
DECLARE_SOA_COLUMN(RICHDeltaPi, richDeltaPi, float);         //! signal - exp. signal for pions
DECLARE_SOA_COLUMN(RICHDeltaKa, richDeltaKa, float);         //! signal - exp. signal for kaons
DECLARE_SOA_COLUMN(RICHDeltaPr, richDeltaPr, float);         //! signal - exp. signal for protons
DECLARE_SOA_COLUMN(RICHNsigmaEl, richNsigmaEl, float);       //! nsigma separation for electrons
DECLARE_SOA_COLUMN(RICHNsigmaMu, richNsigmaMu, float);       //! nsigma separation for muons
DECLARE_SOA_COLUMN(RICHNsigmaPi, richNsigmaPi, float);       //! nsigma separation for pions
DECLARE_SOA_COLUMN(RICHNsigmaKa, richNsigmaKa, float);       //! nsigma separation for kaons
DECLARE_SOA_COLUMN(RICHNsigmaPr, richNsigmaPr, float);       //! nsigma separation for protons
DECLARE_SOA_DYNAMIC_COLUMN(RICHDelta, richDelta,             //! Delta separation with the RICH detector for the combined species
                           [](const float& El, const float& Mu, const float& Pi,
                              const float& Ka, const float& Pr, const o2::track::PID::ID& index) -> float {
                             switch (index) {
                               case o2::track::PID::Electron:
                                 return El;
                               case o2::track::PID::Muon:
                                 return Mu;
                               case o2::track::PID::Pion:
                                 return Pi;
                               case o2::track::PID::Kaon:
                                 return Ka;
                               case o2::track::PID::Proton:
                                 return Pr;
                               default:
                                 return -999.f;
                             }
                           });
DECLARE_SOA_DYNAMIC_COLUMN(RICHNsigma, richNsigma, //! Nsigma separation with the RICH detector for the combined species
                           [](const float& El, const float& Mu, const float& Pi,
                              const float& Ka, const float& Pr, const o2::track::PID::ID& index) -> float {
                             switch (index) {
                               case o2::track::PID::Electron:
                                 return El;
                               case o2::track::PID::Muon:
                                 return Mu;
                               case o2::track::PID::Pion:
                                 return Pi;
                               case o2::track::PID::Kaon:
                                 return Ka;
                               case o2::track::PID::Proton:
                                 return Pr;
                               default:
                                 return -999.f;
                             }
                           });
} // namespace alice3rich

namespace alice3frich
{
DECLARE_SOA_INDEX_COLUMN(Track, track);                        //! Index to travel from track to FRICH
DECLARE_SOA_COLUMN(FRICHSignal, frichSignal, float);           //! Signal in RICH
DECLARE_SOA_COLUMN(FRICHSignalError, frichSignalError, float); //! Error on the RICH signal
DECLARE_SOA_COLUMN(FRICHDeltaEl, frichDeltaEl, float);         //! signal - exp. signal for electrons
DECLARE_SOA_COLUMN(FRICHDeltaMu, frichDeltaMu, float);         //! signal - exp. signal for muons
DECLARE_SOA_COLUMN(FRICHDeltaPi, frichDeltaPi, float);         //! signal - exp. signal for pions
DECLARE_SOA_COLUMN(FRICHDeltaKa, frichDeltaKa, float);         //! signal - exp. signal for kaons
DECLARE_SOA_COLUMN(FRICHDeltaPr, frichDeltaPr, float);         //! signal - exp. signal for protons
DECLARE_SOA_COLUMN(FRICHNsigmaEl, frichNsigmaEl, float);       //! nsigma separation for electrons
DECLARE_SOA_COLUMN(FRICHNsigmaMu, frichNsigmaMu, float);       //! nsigma separation for muons
DECLARE_SOA_COLUMN(FRICHNsigmaPi, frichNsigmaPi, float);       //! nsigma separation for pions
DECLARE_SOA_COLUMN(FRICHNsigmaKa, frichNsigmaKa, float);       //! nsigma separation for kaons
DECLARE_SOA_COLUMN(FRICHNsigmaPr, frichNsigmaPr, float);       //! nsigma separation for protons
DECLARE_SOA_DYNAMIC_COLUMN(FRICHDelta, frichDelta,             //! Delta separation with the FRICH detector for the combined species
                           [](const float& El, const float& Mu, const float& Pi,
                              const float& Ka, const float& Pr, const o2::track::PID::ID& index) -> float {
                             switch (index) {
                               case o2::track::PID::Electron:
                                 return El;
                               case o2::track::PID::Muon:
                                 return Mu;
                               case o2::track::PID::Pion:
                                 return Pi;
                               case o2::track::PID::Kaon:
                                 return Ka;
                               case o2::track::PID::Proton:
                                 return Pr;
                               default:
                                 return -999.f;
                             }
                           });
DECLARE_SOA_DYNAMIC_COLUMN(FRICHNsigma, frichNsigma, //! Nsigma separation with the FRICH detector for the combined species
                           [](const float& El, const float& Mu, const float& Pi,
                              const float& Ka, const float& Pr, const o2::track::PID::ID& index) -> float {
                             switch (index) {
                               case o2::track::PID::Electron:
                                 return El;
                               case o2::track::PID::Muon:
                                 return Mu;
                               case o2::track::PID::Pion:
                                 return Pi;
                               case o2::track::PID::Kaon:
                                 return Ka;
                               case o2::track::PID::Proton:
                                 return Pr;
                               default:
                                 return -999.f;
                             }
                           });
} // namespace alice3frich

DECLARE_SOA_TABLE(RICHs, "AOD", "RICH", //! Table for the ALICE3 RICH detector
                  o2::soa::Index<>,
                  alice3rich::TrackId,
                  alice3rich::RICHSignal,
                  alice3rich::RICHSignalError,
                  alice3rich::RICHDeltaEl,
                  alice3rich::RICHDeltaMu,
                  alice3rich::RICHDeltaPi,
                  alice3rich::RICHDeltaKa,
                  alice3rich::RICHDeltaPr,
                  alice3rich::RICHNsigmaEl,
                  alice3rich::RICHNsigmaMu,
                  alice3rich::RICHNsigmaPi,
                  alice3rich::RICHNsigmaKa,
                  alice3rich::RICHNsigmaPr,
                  alice3rich::RICHDelta<alice3rich::RICHDeltaEl, alice3rich::RICHDeltaMu, alice3rich::RICHDeltaPi,
                                        alice3rich::RICHDeltaKa, alice3rich::RICHDeltaPr>,
                  alice3rich::RICHNsigma<alice3rich::RICHNsigmaEl, alice3rich::RICHNsigmaMu, alice3rich::RICHNsigmaPi,
                                         alice3rich::RICHNsigmaKa, alice3rich::RICHNsigmaPr>);
using RICH = RICHs::iterator;

DECLARE_SOA_TABLE(FRICHs, "AOD", "FRICH", //! Table for the ALICE3 Forward RICH detector
                  o2::soa::Index<>,
                  alice3frich::TrackId,
                  alice3frich::FRICHSignal,
                  alice3frich::FRICHSignalError,
                  alice3frich::FRICHDeltaEl,
                  alice3frich::FRICHDeltaMu,
                  alice3frich::FRICHDeltaPi,
                  alice3frich::FRICHDeltaKa,
                  alice3frich::FRICHDeltaPr,
                  alice3frich::FRICHNsigmaEl,
                  alice3frich::FRICHNsigmaMu,
                  alice3frich::FRICHNsigmaPi,
                  alice3frich::FRICHNsigmaKa,
                  alice3frich::FRICHNsigmaPr,
                  alice3frich::FRICHDelta<alice3frich::FRICHDeltaEl, alice3frich::FRICHDeltaMu, alice3frich::FRICHDeltaPi,
                                          alice3frich::FRICHDeltaKa, alice3frich::FRICHDeltaPr>,
                  alice3frich::FRICHNsigma<alice3frich::FRICHNsigmaEl, alice3frich::FRICHNsigmaMu, alice3frich::FRICHNsigmaPi,
                                           alice3frich::FRICHNsigmaKa, alice3frich::FRICHNsigmaPr>);

using FRICH = FRICHs::iterator;

} // namespace o2::aod

#endif // O2_ANALYSIS_ALICE3_RICH_H_
