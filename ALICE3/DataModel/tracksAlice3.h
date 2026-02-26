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
/// \file   tracksAlice3.h
/// \author David Dobrigkeit Chinellato
/// \since  11/05/2023
/// \brief  Table for ALICE 3 track-related info
///

#ifndef ALICE3_DATAMODEL_TRACKSALICE3_H_
#define ALICE3_DATAMODEL_TRACKSALICE3_H_

// O2 includes
#include "Framework/AnalysisDataModel.h"
#include <Framework/ASoA.h>

namespace o2::aod
{
namespace track_alice3
{
DECLARE_SOA_COLUMN(IsReconstructed, isReconstructed, bool); //! is reconstructed or not
DECLARE_SOA_COLUMN(NSiliconHits, nSiliconHits, int);        //! number of silicon hits
DECLARE_SOA_COLUMN(NTPCHits, nTPCHits, int);                //! number of tpc hits
DECLARE_SOA_COLUMN(PdgCode, pdgCode, int);                  //! PDG code of the linked truth MC particle
} // namespace track_alice3
DECLARE_SOA_TABLE(TracksAlice3, "AOD", "TRACKSALICE3",
                  track_alice3::IsReconstructed);
using TrackAlice3 = TracksAlice3::iterator;

DECLARE_SOA_TABLE(TracksAlice3Pdg, "AOD", "TRACKSALICE3PDG",
                  track_alice3::PdgCode);
using TrackAlice3Pdg = TracksAlice3Pdg::iterator;

DECLARE_SOA_TABLE(TracksExtraA3, "AOD", "TracksExtraA3",
                  track_alice3::NSiliconHits,
                  track_alice3::NTPCHits);
using TrackExtraA3 = TracksExtraA3::iterator;

namespace mcparticle_alice3
{
DECLARE_SOA_COLUMN(NHits, nHits, int);     //! number of silicon hits
DECLARE_SOA_COLUMN(Charge, charge, float); //! particle charge
} // namespace mcparticle_alice3
DECLARE_SOA_TABLE(MCParticlesExtraA3, "AOD", "MCParticlesExtraA3",
                  mcparticle_alice3::NHits,
                  mcparticle_alice3::Charge);
using MCParticleExtraA3 = MCParticlesExtraA3::iterator;
} // namespace o2::aod

#endif // ALICE3_DATAMODEL_TRACKSALICE3_H_
