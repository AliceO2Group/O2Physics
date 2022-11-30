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
/// \file   tofSkimsTableCreator.h
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
/// \since  19/11/2022
/// \brief  Definition of the skimmed data format for the TOF skims
///

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::track;
using namespace o2::dataformats;

#include "Common/TableProducer/PID/pidTOFBase.h"
#include "Common/DataModel/PIDResponse.h"

namespace o2::aod
{
namespace tofskims
{
DECLARE_SOA_COLUMN(P, p, float);                             //! Momentum of the track
DECLARE_SOA_COLUMN(Pt, pt, float);                           //! Pt of the track
DECLARE_SOA_COLUMN(Eta, eta, float);                         //! Eta of the track
DECLARE_SOA_COLUMN(Phi, phi, float);                         //! Phi of the track
DECLARE_SOA_COLUMN(PIDForTracking, pidForTracking, uint8_t); //! Index for mass hypothesis used in tracking see PID.h for definition
DECLARE_SOA_DYNAMIC_COLUMN(HasTOF, hasTOF,                   //! Flag to check if track has a TOF measurement
                           [](float tofSignal) -> bool { return tofSignal > 0; });

} // namespace tofskims

DECLARE_SOA_TABLE(SkimmedTOF, "AOD", "SKIMMEDTOF", //! Table of the skimmed TOF data format. One entry per track.
                  tofskims::P,
                  tofskims::Pt,
                  tofskims::Eta,
                  tofskims::Phi,
                  tofskims::PIDForTracking,
                  track::TOFExpMom,
                  track::Length,
                  track::TOFChi2,
                  pidtofsignal::TOFSignal,
                  tofskims::HasTOF<pidtofsignal::TOFSignal>);

} // namespace o2::aod