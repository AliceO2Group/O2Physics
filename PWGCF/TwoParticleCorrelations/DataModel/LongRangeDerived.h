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
DECLARE_SOA_COLUMN(TotalFT0AmplitudeA, totalFT0AmplitudeA, float); //! sum of amplitudes on A side of FT0
DECLARE_SOA_COLUMN(TotalFT0AmplitudeC, totalFT0AmplitudeC, float); //! sum of amplitudes on C side of FT0
DECLARE_SOA_COLUMN(TotalFV0AmplitudeA, totalFV0AmplitudeA, float); //! sum of amplitudes on A side of FDD
DECLARE_SOA_COLUMN(GapSide, gapSide, uint8_t);                     // 0 for side A, 1 for side C, 2 for both sides
} // namespace lrcorrcolltable

DECLARE_SOA_TABLE(CollLRTables, "AOD", "COLLLRTABLE",
                  o2::soa::Index<>,
                  bc::RunNumber,
                  lrcorrcolltable::Zvtx,
                  lrcorrcolltable::Multiplicity,
                  lrcorrcolltable::Centrality,
                  timestamp::Timestamp);
using CollLRTable = CollLRTables::iterator;

DECLARE_SOA_TABLE(UpcCollLRTables, "AOD", "UPCCOLLLRTABLE",
                  o2::soa::Index<>,
                  bc::GlobalBC,
                  bc::RunNumber,
                  lrcorrcolltable::Zvtx,
                  lrcorrcolltable::Multiplicity,
                  lrcorrcolltable::TotalFT0AmplitudeA,
                  lrcorrcolltable::TotalFT0AmplitudeC,
                  lrcorrcolltable::TotalFV0AmplitudeA);
using UpcCollLRTable = UpcCollLRTables::iterator;

DECLARE_SOA_TABLE(UpcSgCollLRTables, "AOD", "UPCSGCOLLLRTABLE",
                  lrcorrcolltable::GapSide);
using UpcSgCollLRTable = UpcSgCollLRTables::iterator;

namespace lrcorrzdctable
{
DECLARE_SOA_INDEX_COLUMN(UpcCollLRTable, upcCollLRTable);
DECLARE_SOA_COLUMN(EnergyCommonZNA, energyCommonZNA, float);
DECLARE_SOA_COLUMN(EnergyCommonZNC, energyCommonZNC, float);
} // namespace lrcorrzdctable

DECLARE_SOA_TABLE(ZdcLRTables, "AOD", "ZDCLRTABLE",
                  o2::soa::Index<>,
                  lrcorrzdctable::UpcCollLRTableId,
                  lrcorrzdctable::EnergyCommonZNA,
                  lrcorrzdctable::EnergyCommonZNC);
using ZdcLRTable = ZdcLRTables::iterator;

namespace lrcorrtrktable
{
DECLARE_SOA_INDEX_COLUMN(CollLRTable, collLRTable);
DECLARE_SOA_INDEX_COLUMN(UpcCollLRTable, upcCollLRTable);
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

DECLARE_SOA_TABLE(TrkLRUpcTables, "AOD", "TRKLRUPCTABLE",
                  o2::soa::Index<>,
                  lrcorrtrktable::UpcCollLRTableId,
                  lrcorrtrktable::Pt,
                  lrcorrtrktable::Eta,
                  lrcorrtrktable::Phi,
                  lrcorrtrktable::TrackType);
using TrkLRUpcTable = TrkLRUpcTables::iterator;

DECLARE_SOA_TABLE(Ft0aLRUpcTables, "AOD", "FT0ALRUpcTABLE",
                  o2::soa::Index<>,
                  lrcorrtrktable::UpcCollLRTableId,
                  lrcorrtrktable::ChannelID,
                  lrcorrtrktable::Amplitude,
                  lrcorrtrktable::Eta,
                  lrcorrtrktable::Phi);
using Ft0aLRUpcTable = Ft0aLRUpcTables::iterator;

DECLARE_SOA_TABLE(Ft0cLRUpcTables, "AOD", "FT0CLRUpcTABLE",
                  o2::soa::Index<>,
                  lrcorrtrktable::UpcCollLRTableId,
                  lrcorrtrktable::ChannelID,
                  lrcorrtrktable::Amplitude,
                  lrcorrtrktable::Eta,
                  lrcorrtrktable::Phi);
using Ft0cLRUpcTable = Ft0cLRUpcTables::iterator;

DECLARE_SOA_TABLE(V0TrkLRUpcTables, "AOD", "V0TRKLRUPCTABLE",
                  o2::soa::Index<>,
                  lrcorrtrktable::UpcCollLRTableId,
                  lrcorrtrktable::IdPos,
                  lrcorrtrktable::IdNeg,
                  lrcorrtrktable::Pt,
                  lrcorrtrktable::Eta,
                  lrcorrtrktable::Phi,
                  lrcorrtrktable::InvMass,
                  lrcorrtrktable::V0Type);
using V0TrkLRUpcTable = V0TrkLRUpcTables::iterator;

DECLARE_SOA_TABLE(MftTrkLRUpcTables, "AOD", "MFTTRKLRUPCTABLE",
                  o2::soa::Index<>,
                  lrcorrtrktable::UpcCollLRTableId,
                  lrcorrtrktable::Pt,
                  lrcorrtrktable::Eta,
                  lrcorrtrktable::Phi);
using MftTrkLRUpcTable = MftTrkLRUpcTables::iterator;

DECLARE_SOA_TABLE(MftBestTrkLRUpcTables, "AOD", "MFTBESTTRKLRUPCTABLE",
                  o2::soa::Index<>,
                  lrcorrtrktable::UpcCollLRTableId,
                  lrcorrtrktable::Pt,
                  lrcorrtrktable::Eta,
                  lrcorrtrktable::Phi);
using MftBestTrkLRUpcTable = MftBestTrkLRUpcTables::iterator;

} // namespace o2::aod

#endif // PWGCF_TWOPARTICLECORRELATIONS_DATAMODEL_LONGRANGEDERIVED_H_
