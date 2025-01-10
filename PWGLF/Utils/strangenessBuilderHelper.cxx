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
/// \file  strangenessBuilderHelper.cxx
/// \since 07/01/2025
/// \brief Utilities to build strange decays in a modular way. Used by strangeness builder task
///

#include "strangenessBuilderHelper.h"
#include "ReconstructionDataFormats/Track.h"
#include "DetectorsBase/GeometryManager.h"
#include "CommonConstants/PhysicsConstants.h"
#include "Framework/AnalysisDataModel.h"
#include "Common/Core/trackUtilities.h"

namespace o2
{
namespace pwglf
{
//_______________________________________________________________
strangenessBuilderHelper::strangenessBuilderHelper() {

}

//_______________________________________________________________
// builds V0 from two tracks. Does not check any conditionals
// except for DCA fitter convergence (should be minimal overhead)
// Resulting properties can be checked with strangenessBuilderHelper::v0

template <typename TTrack>
bool strangenessBuilderHelper::buildV0Candidate

//_______________________________________________________________
// internal helper to calculate DCAxy of a straight line to a given PV analytically
float strangenessBuilderHelper::CalculateDCAStraightToPV(float X, float Y, float Z, float Px, float Py, float Pz, float pvX, float pvY, float pvZ)


//_______________________________________________________________
} // namespace pwglf
} // namespace o2