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
/// \file  TrackSelectionDefaults.h
/// \brief Class for the definition of standard track selection objects
/// \since 20-10-2020
///

#ifndef TrackSelectionDefaults_H
#define TrackSelectionDefaults_H

#include "Framework/DataTypes.h"
#include "Common/Core/TrackSelection.h"

// Default track selection requiring one hit in the SPD
TrackSelection getGlobalTrackSelection();

// Default track selection requiring a particular Run 3 ITS matching
TrackSelection getGlobalTrackSelectionRun3ITSMatch(int matching);

// Default track selection requiring no hit in the SPD and one in the innermost
// SDD -> complementary tracks to global selection
TrackSelection getGlobalTrackSelectionSDD();

#endif
