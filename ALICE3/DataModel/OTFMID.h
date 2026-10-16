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
/// \file   OTFMID.h
/// \author David Dobrigkeit Chinellato
/// \author Antonio Ortiz, UNAM
/// \since  11/05/2023
/// \brief  Set of tables for the ALICE3 OTFMID information
/// TODO: remove the columns that are not usefulf for me

#ifndef ALICE3_DATAMODEL_MID_H_
#define ALICE3_DATAMODEL_MID_H_

#include <Framework/AnalysisDataModel.h>
#include <cstdlib>
#include <vector>

namespace o2::aod
{
namespace upgrade_mid
{
DECLARE_SOA_COLUMN(MIDTrackTime, midTrackTime, float);     //! Track time generated at the MID
DECLARE_SOA_COLUMN(MIDTrackLength, midTrackLength, float); //! track length for calculation of MID (generated)

DECLARE_SOA_COLUMN(MIDEventTime, midEventTime, float);                       //! Event time reconstructed with the MID
DECLARE_SOA_COLUMN(MIDEventTimeErr, midEventTimeErr, float);                 //! Uncertainty on the event time reconstructed with the MID
DECLARE_SOA_COLUMN(MIDTrackPosReco, midTrackPosReco, std::vector<float>);  //! Measured position in MID
DECLARE_SOA_COLUMN(MIDTrackPos, midTrackPos, std::vector<float>);  //! Expected position in MID
DECLARE_SOA_COLUMN(MIDCandidate, midCandidate, bool);              //! Muon candidate
DECLARE_SOA_COLUMN(MIDTrackTimeReco, midTrackTimeReco, float);     //! Track time measured at MID
DECLARE_SOA_COLUMN(MIDTrackLengthReco, midTrackLengthReco, float); //! track length for calculation of MID (reconstructed)

} // namespace upgrade_mid

DECLARE_SOA_TABLE(UpgradeMidMCs, "AOD", "UPGRADEMIDMC",
                  upgrade_mid::MIDTrackTime,
                  upgrade_mid::MIDTrackLength);

DECLARE_SOA_TABLE(UpgradeMids, "AOD", "UPGRADEMID",
                  upgrade_mid::MIDEventTime,
                  upgrade_mid::MIDEventTimeErr,
                  upgrade_mid::MIDTrackPosReco,
                  upgrade_mid::MIDTrackTimeReco,
                  upgrade_mid::MIDTrackLengthReco,
                  upgrade_mid::MIDCandidate);

DECLARE_SOA_TABLE(UpgradeMidExpecteds, "AOD", "UPGRADEMIDEXPT",
                  upgrade_mid::MIDTrackPos,
                  upgrade_mid::MIDTrackTime,
                  upgrade_mid::MIDTrackLength);



using UpgradeMidMC = UpgradeMidMCs::iterator;
using UpgradeMid = UpgradeMids::iterator;
using UpgradeMidExpected = UpgradeMidExpecteds::iterator;

} // namespace o2::aod

#endif // ALICE3_DATAMODEL_OTFMID_H_
