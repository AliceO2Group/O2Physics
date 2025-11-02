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
/// \file LongRangeDerived.h
///
/// \brief task derived table definition for long range correlation
/// \author Abhi Modak (abhi.modak@cern.ch)
/// \since October 28, 2025

#ifndef PWGCF_TWOPARTICLECORRELATIONS_DATAMODEL_LONGRANGEDERIVED_H_
#define PWGCF_TWOPARTICLECORRELATIONS_DATAMODEL_LONGRANGEDERIVED_H_

#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"

namespace o2::aod
{
namespace lrcorrcolltable
{
DECLARE_SOA_COLUMN(Zvtx, zvtx, float);
DECLARE_SOA_COLUMN(Multiplicity, multiplicity, float);
DECLARE_SOA_COLUMN(Centrality, centrality, float);
} // namespace lrcorrcolltable

DECLARE_SOA_TABLE(CollLRTables, "AOD", "COLLLRTABLE",
                  o2::soa::Index<>,
                  bc::RunNumber,
                  lrcorrcolltable::Zvtx,
                  lrcorrcolltable::Multiplicity,
                  lrcorrcolltable::Centrality,
                  timestamp::Timestamp);
using CollLRTable = CollLRTables::iterator;

namespace lrcorrtrktable
{
DECLARE_SOA_INDEX_COLUMN(CollLRTable, collLRTable);
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(ChannelID, channelID, int);
DECLARE_SOA_COLUMN(Amplitude, amplitude, float);
DECLARE_SOA_COLUMN(InvMass, invMass, float);
DECLARE_SOA_COLUMN(IdPos, idPos, int64_t);
DECLARE_SOA_COLUMN(IdNeg, idNeg, int64_t);
DECLARE_SOA_COLUMN(TrackType, trackType, uint8_t);
DECLARE_SOA_COLUMN(V0Type, v0Type, uint8_t);
enum TrackPid {
  kSpCharge,
  kSpPion,
  kSpKaon,
  kSpProton
};
enum V0TrackPid {
  kSpK0short,
  kSpLambda,
  kSpALambda
};
} // namespace lrcorrtrktable

DECLARE_SOA_TABLE(TrkLRTables, "AOD", "TRKLRTABLE",
                  o2::soa::Index<>,
                  lrcorrtrktable::CollLRTableId,
                  lrcorrtrktable::Pt,
                  lrcorrtrktable::Eta,
                  lrcorrtrktable::Phi,
                  lrcorrtrktable::TrackType);
using TrkLRTable = TrkLRTables::iterator;

DECLARE_SOA_TABLE(Ft0aLRTables, "AOD", "FT0ALRTABLE",
                  o2::soa::Index<>,
                  lrcorrtrktable::CollLRTableId,
                  lrcorrtrktable::ChannelID,
                  lrcorrtrktable::Amplitude,
                  lrcorrtrktable::Eta,
                  lrcorrtrktable::Phi);
using Ft0aLRTable = Ft0aLRTables::iterator;

DECLARE_SOA_TABLE(Ft0cLRTables, "AOD", "FT0CLRTABLE",
                  o2::soa::Index<>,
                  lrcorrtrktable::CollLRTableId,
                  lrcorrtrktable::ChannelID,
                  lrcorrtrktable::Amplitude,
                  lrcorrtrktable::Eta,
                  lrcorrtrktable::Phi);
using Ft0cLRTable = Ft0cLRTables::iterator;

DECLARE_SOA_TABLE(V0TrkLRTables, "AOD", "V0TRKLRTABLE",
                  o2::soa::Index<>,
                  lrcorrtrktable::CollLRTableId,
                  lrcorrtrktable::IdPos,
                  lrcorrtrktable::IdNeg,
                  lrcorrtrktable::Pt,
                  lrcorrtrktable::Eta,
                  lrcorrtrktable::Phi,
                  lrcorrtrktable::InvMass,
                  lrcorrtrktable::V0Type);
using V0TrkLRTable = V0TrkLRTables::iterator;

DECLARE_SOA_TABLE(MftTrkLRTables, "AOD", "MFTTRKLRTABLE",
                  o2::soa::Index<>,
                  lrcorrtrktable::CollLRTableId,
                  lrcorrtrktable::Pt,
                  lrcorrtrktable::Eta,
                  lrcorrtrktable::Phi);
using MftTrkLRTable = MftTrkLRTables::iterator;

DECLARE_SOA_TABLE(MftBestTrkLRTables, "AOD", "MFTBESTTRKLRTABLE",
                  o2::soa::Index<>,
                  lrcorrtrktable::CollLRTableId,
                  lrcorrtrktable::Pt,
                  lrcorrtrktable::Eta,
                  lrcorrtrktable::Phi);
using MftBestTrkLRTable = MftBestTrkLRTables::iterator;

} // namespace o2::aod

#endif // PWGCF_TWOPARTICLECORRELATIONS_DATAMODEL_LONGRANGEDERIVED_H_
