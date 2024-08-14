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

// Table for the event time. One entry per collision
namespace o2::aod
{
namespace evtime
{

DECLARE_SOA_COLUMN(EvTime, evTime, float);             //! Event time. Can be obtained via a combination of detectors e.g. TOF, FT0A, FT0C
DECLARE_SOA_COLUMN(EvTimeErr, evTimeErr, float);       //! Error of event time. Can be obtained via a combination of detectors e.g. TOF, FT0A, FT0C
DECLARE_SOA_COLUMN(EvTimeTOF, evTimeTOF, float);       //! Event time computed with the TOF detector
DECLARE_SOA_COLUMN(EvTimeTOFErr, evTimeTOFErr, float); //! Error of the event time computed with the TOF detector
DECLARE_SOA_COLUMN(EvTimeFT0, evTimeFT0, float);       //! Event time computed with the FT0 detector
DECLARE_SOA_COLUMN(EvTimeFT0Err, evTimeFT0Err, float); //! Error of the event time computed with the FT0 detector
} // namespace evtime

DECLARE_SOA_TABLE(EvTimeTOFFT0, "AOD", "EvTimeTOFFT0", //! Table of the TOF event time. One entry per collision.
                  evtime::EvTime,
                  evtime::EvTimeErr,
                  evtime::EvTimeTOF,
                  evtime::EvTimeTOFErr,
                  evtime::EvTimeFT0,
                  evtime::EvTimeFT0Err);
} // namespace o2::aod
