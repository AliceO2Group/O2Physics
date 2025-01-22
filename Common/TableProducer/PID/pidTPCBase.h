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
/// \file   pidTPCBase.h
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
/// \brief  Base to build tasks for TPC PID tasks.
///

#ifndef COMMON_TABLEPRODUCER_PID_PIDTPCBASE_H_
#define COMMON_TABLEPRODUCER_PID_PIDTPCBASE_H_

#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTPC.h"

namespace o2::aod
{

DECLARE_SOA_TABLE(PIDMults, "AOD", "PIDMults", //! TPC auxiliary table for the PID
                  o2::soa::Marker<1>,
                  mult::MultTPC);
using PIDMult = PIDMults::iterator;

} // namespace o2::aod

int getPIDIndex(const int pdgCode) // Get O2 PID index corresponding to MC PDG code
{
  switch (abs(pdgCode)) {
    case 11:
      return o2::track::PID::Electron;
    case 13:
      return o2::track::PID::Muon;
    case 211:
      return o2::track::PID::Pion;
    case 321:
      return o2::track::PID::Kaon;
    case 2212:
      return o2::track::PID::Proton;
    case 1000010020:
      return o2::track::PID::Deuteron;
    case 1000010030:
      return o2::track::PID::Triton;
    case 1000020030:
      return o2::track::PID::Helium3;
    case 1000020040:
      return o2::track::PID::Alpha;
    default: // treat as pion if not any of the above
      return o2::track::PID::Pion;
  }
}

#endif // COMMON_TABLEPRODUCER_PID_PIDTPCBASE_H_
