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

// **Hadronic Correction in the EMCAL framework: to avoid the double counting of the charged particles' contribution in jets**
/// \author Archita Rani Dash <archita.rani.dash@cern.ch>

#ifndef PWGJE_DATAMODEL_EMCALCLUSTERHADRONICCORRECTION_H_
#define PWGJE_DATAMODEL_EMCALCLUSTERHADRONICCORRECTION_H_

#include <string>
#include "Framework/AnalysisDataModel.h"
#include "EMCALClusterDefinition.h"

namespace o2::aod
{

namespace emcalhadroniccorrection
{
// four different columns for the 'Cluster Energies after Hadronic Correction' to make it more flexible
// 2 for the closest match and 2 for all matched tracks

// cluster values
// DECLARE_SOA_COLUMN(HadCorrEnergy, hadcorrEnergy, float);          //! cluster energy (GeV) after hadronic correction

// hadronic corrected energy values
DECLARE_SOA_COLUMN(HadCorrOneTrack1, hadCorrOneTrack1, float);   //! with hadronic correction fraction (100%) for one matched track
DECLARE_SOA_COLUMN(HadCorrOneTrack2, hadCorrOneTrack2, float);   //! with hadronic correction fraction (70%) for one matched track - systematic studies
DECLARE_SOA_COLUMN(HadCorrAllTracks1, hadCorrAllTracks1, float); //! with hadronic correction fraction (100%) for all matched tracks
DECLARE_SOA_COLUMN(HadCorrAllTracks2, hadCorrAllTracks2, float); //! with hadronic correction fraction (70%) for all matched tracks - for systematic studies

} // namespace emcalhadroniccorrection

// Table Definitions - define what needs to be written into the tables produced by this tableproducer task
DECLARE_SOA_TABLE(EmcalHCs, "AOD", "EMCALHCS",                //!
                  o2::soa::Index<>,                           //!
                  emcalhadroniccorrection::HadCorrOneTrack1,  // corrected cluster energy for 1 matched track (f = 100%)
                  emcalhadroniccorrection::HadCorrOneTrack2,  // corrected cluster energy for 1 matched track (f = 70%)
                  emcalhadroniccorrection::HadCorrAllTracks1, // corrected cluster energy for all matched tracks (f = 100%)
                  emcalhadroniccorrection::HadCorrAllTracks2  // corrected cluster energy for all matched tracks (f = 70%)
)

using EmcalHC = EmcalHCs::iterator;

} // namespace o2::aod
#endif // PWGJE_DATAMODEL_EMCALCLUSTERHADRONICCORRECTION_H_
