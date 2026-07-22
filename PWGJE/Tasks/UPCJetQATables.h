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

#ifndef PWGJE_TASKS_UPCJETQATABLES_H_
#define PWGJE_TASKS_UPCJETQATABLES_H_

#include <Framework/AnalysisDataModel.h>

namespace o2::aod
{

// ============================================================
// DATA tables
// ============================================================
namespace upcjetevent
{
DECLARE_SOA_COLUMN(GapSide, gapSide, int);
DECLARE_SOA_COLUMN(PosZ, posZ, float);
DECLARE_SOA_COLUMN(CentFT0M, centFT0M, float);
DECLARE_SOA_COLUMN(EnergyCommonZNA, energyCommonZNA, float);
DECLARE_SOA_COLUMN(EnergyCommonZNC, energyCommonZNC, float);
DECLARE_SOA_COLUMN(TimeZNA, timeZNA, float);
DECLARE_SOA_COLUMN(TimeZNC, timeZNC, float);
DECLARE_SOA_COLUMN(NeutronClass, neutronClass, int);
DECLARE_SOA_COLUMN(IsGapTagged, isGapTagged, bool);
DECLARE_SOA_COLUMN(IsZDCTagged, isZDCTagged, bool);
} // namespace upcjetevent

DECLARE_SOA_TABLE(UpcJetEvents, "AOD", "UPCJETEVENT",
                  upcjetevent::GapSide,
                  upcjetevent::PosZ,
                  upcjetevent::CentFT0M,
                  upcjetevent::EnergyCommonZNA,
                  upcjetevent::EnergyCommonZNC,
                  upcjetevent::TimeZNA,
                  upcjetevent::TimeZNC,
                  upcjetevent::NeutronClass,
                  upcjetevent::IsGapTagged,
                  upcjetevent::IsZDCTagged);
using UpcJetEvent = UpcJetEvents::iterator;

namespace upcjet
{
DECLARE_SOA_INDEX_COLUMN(UpcJetEvent, upcJetEvent);
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(Area, area, float);
DECLARE_SOA_COLUMN(NConstituents, nConstituents, int);
} // namespace upcjet

DECLARE_SOA_TABLE(UpcJets, "AOD", "UPCJET",
                  upcjet::UpcJetEventId,
                  upcjet::Pt,
                  upcjet::Eta,
                  upcjet::Phi,
                  upcjet::Area,
                  upcjet::NConstituents);
using UpcJet = UpcJets::iterator;

namespace upctrack
{
DECLARE_SOA_INDEX_COLUMN_FULL(UpcJetEvent, upcJetEvent, int32_t, UpcJetEvents, "");
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
} // namespace upctrack

DECLARE_SOA_TABLE(UpcTracks, "AOD", "UPCTRACK",
                  upctrack::UpcJetEventId,
                  upctrack::Pt,
                  upctrack::Eta,
                  upctrack::Phi);
using UpcTrack = UpcTracks::iterator;

namespace upcjettrack
{
DECLARE_SOA_INDEX_COLUMN_FULL(UpcJet, upcJet, int32_t, UpcJets, "");
DECLARE_SOA_INDEX_COLUMN_FULL(UpcJetEvent, upcJetEvent, int32_t, UpcJetEvents, "");
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
} // namespace upcjettrack

DECLARE_SOA_TABLE(UpcJetTracks, "AOD", "UPCJETTRACK",
                  upcjettrack::UpcJetId,
                  upcjettrack::UpcJetEventId,
                  upcjettrack::Pt,
                  upcjettrack::Eta,
                  upcjettrack::Phi);
using UpcJetTrack = UpcJetTracks::iterator;

// ============================================================
// MCD tables
// ============================================================
namespace upcjeteventmcd
{
DECLARE_SOA_COLUMN(GapSide, gapSide, int);
DECLARE_SOA_COLUMN(PosZ, posZ, float);
DECLARE_SOA_COLUMN(CentFT0M, centFT0M, float);
DECLARE_SOA_COLUMN(EnergyCommonZNA, energyCommonZNA, float);
DECLARE_SOA_COLUMN(EnergyCommonZNC, energyCommonZNC, float);
DECLARE_SOA_COLUMN(TimeZNA, timeZNA, float);
DECLARE_SOA_COLUMN(TimeZNC, timeZNC, float);
DECLARE_SOA_COLUMN(NeutronClass, neutronClass, int);
DECLARE_SOA_COLUMN(IsGapTagged, isGapTagged, bool);
DECLARE_SOA_COLUMN(IsZDCTagged, isZDCTagged, bool);
} // namespace upcjeteventmcd

DECLARE_SOA_TABLE(UpcJetEventsMCD, "AOD", "UPCJETEVENTMCD",
                  upcjeteventmcd::GapSide,
                  upcjeteventmcd::PosZ,
                  upcjeteventmcd::CentFT0M,
                  upcjeteventmcd::EnergyCommonZNA,
                  upcjeteventmcd::EnergyCommonZNC,
                  upcjeteventmcd::TimeZNA,
                  upcjeteventmcd::TimeZNC,
                  upcjeteventmcd::NeutronClass,
                  upcjeteventmcd::IsGapTagged,
                  upcjeteventmcd::IsZDCTagged);
using UpcJetEventMCD = UpcJetEventsMCD::iterator;

namespace upctrackmcd
{
DECLARE_SOA_INDEX_COLUMN_FULL(UpcJetEventMCD, upcJetEventMCD, int32_t, UpcJetEventsMCD, "");
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
} // namespace upctrackmcd

DECLARE_SOA_TABLE(UpcTracksMCD, "AOD", "UPCTRACKMCD",
                  upctrackmcd::UpcJetEventMCDId,
                  upctrackmcd::Pt,
                  upctrackmcd::Eta,
                  upctrackmcd::Phi);
using UpcTrackMCD = UpcTracksMCD::iterator;

namespace upcjetmcd
{
DECLARE_SOA_INDEX_COLUMN_FULL(UpcJetEventMCD, upcJetEventMCD, int32_t, UpcJetEventsMCD, "");
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(Area, area, float);
DECLARE_SOA_COLUMN(NConstituents, nConstituents, int);
} // namespace upcjetmcd

DECLARE_SOA_TABLE(UpcJetsMCD, "AOD", "UPCJETMCD",
                  upcjetmcd::UpcJetEventMCDId,
                  upcjetmcd::Pt,
                  upcjetmcd::Eta,
                  upcjetmcd::Phi,
                  upcjetmcd::Area,
                  upcjetmcd::NConstituents);
using UpcJetMCD = UpcJetsMCD::iterator;

namespace upcjettrackmcd
{
DECLARE_SOA_INDEX_COLUMN_FULL(UpcJetMCD, upcJetMCD, int32_t, UpcJetsMCD, "");
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
} // namespace upcjettrackmcd

DECLARE_SOA_TABLE(UpcJetTracksMCD, "AOD", "UPCJETTRACKMCD",
                  upcjettrackmcd::UpcJetMCDId,
                  upcjettrackmcd::Pt,
                  upcjettrackmcd::Eta,
                  upcjettrackmcd::Phi);
using UpcJetTrackMCD = UpcJetTracksMCD::iterator;

// ============================================================
// MCP tables
// ============================================================
namespace upcjeteventmcp
{
DECLARE_SOA_COLUMN(PosZ, posZ, float);
} // namespace upcjeteventmcp

DECLARE_SOA_TABLE(UpcJetEventsMCP, "AOD", "UPCJETEVENTMCP",
                  upcjeteventmcp::PosZ);
using UpcJetEventMCP = UpcJetEventsMCP::iterator;

namespace upctrackmcp
{
DECLARE_SOA_INDEX_COLUMN_FULL(UpcJetEventMCP, upcJetEventMCP, int32_t, UpcJetEventsMCP, "");
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
} // namespace upctrackmcp

DECLARE_SOA_TABLE(UpcTracksMCP, "AOD", "UPCTRACKMCP",
                  upctrackmcp::UpcJetEventMCPId,
                  upctrackmcp::Pt,
                  upctrackmcp::Eta,
                  upctrackmcp::Phi);
using UpcTrackMCP = UpcTracksMCP::iterator;

namespace upcjetmcp
{
DECLARE_SOA_INDEX_COLUMN_FULL(UpcJetEventMCP, upcJetEventMCP, int32_t, UpcJetEventsMCP, "");
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(Area, area, float);
DECLARE_SOA_COLUMN(NConstituents, nConstituents, int);
} // namespace upcjetmcp

DECLARE_SOA_TABLE(UpcJetsMCP, "AOD", "UPCJETMCP",
                  upcjetmcp::UpcJetEventMCPId,
                  upcjetmcp::Pt,
                  upcjetmcp::Eta,
                  upcjetmcp::Phi,
                  upcjetmcp::Area,
                  upcjetmcp::NConstituents);
using UpcJetMCP = UpcJetsMCP::iterator;

namespace upcjettrackmcp
{
DECLARE_SOA_INDEX_COLUMN_FULL(UpcJetMCP, upcJetMCP, int32_t, UpcJetsMCP, "");
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
} // namespace upcjettrackmcp

DECLARE_SOA_TABLE(UpcJetTracksMCP, "AOD", "UPCJETTRACKMCP",
                  upcjettrackmcp::UpcJetMCPId,
                  upcjettrackmcp::Pt,
                  upcjettrackmcp::Eta,
                  upcjettrackmcp::Phi);
using UpcJetTrackMCP = UpcJetTracksMCP::iterator;

} // namespace o2::aod

#endif // PWGJE_TASKS_UPCJETQATABLES_H_
