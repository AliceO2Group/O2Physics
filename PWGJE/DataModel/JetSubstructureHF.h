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
/// \brief Table definitions for hf jet substrucuture observables
///
/// \author Nima Zardoshti

#ifndef PWGJE_DATAMODEL_JETSUBSTRUCTUREHF_H_
#define PWGJE_DATAMODEL_JETSUBSTRUCTUREHF_H_

#include <cmath>
#include "Framework/AnalysisDataModel.h"
#include "PWGJE/DataModel/EMCALClusters.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"

using namespace o2::analysis;

namespace o2::aod
{
namespace jetsubstructurehf
{
DECLARE_SOA_INDEX_COLUMN(HFJet, jet);                                                    //!
DECLARE_SOA_COLUMN(Zg, zg, float);                                                       //!
DECLARE_SOA_COLUMN(Rg, rg, float);                                                       //!
DECLARE_SOA_COLUMN(Nsd, nsd, float);                                                     //!
DECLARE_SOA_DYNAMIC_COLUMN(DataIdentifier, dataIdentifier, []() -> int { return 0; });   //!
DECLARE_SOA_DYNAMIC_COLUMN(MCDetIdentifier, mcDetIdentifier, []() -> int { return 0; }); //!
DECLARE_SOA_DYNAMIC_COLUMN(MCGenIdentifier, McGenIdentifier, []() -> int { return 0; }); //!
} // namespace jetsubstructurehf
DECLARE_SOA_TABLE(JetSubtructureHFData, "AOD", "JETSUBHFData", jetsubstructurehf::HFJetId, jetsubstructurehf::Zg, jetsubstructurehf::Rg, jetsubstructurehf::Nsd, jetsubstructurehf::DataIdentifier<>);  //!
DECLARE_SOA_TABLE(JetSubtructureHFMCDet, "AOD", "JETSUBHFMCD", jetsubstructurehf::HFJetId, jetsubstructurehf::Zg, jetsubstructurehf::Rg, jetsubstructurehf::Nsd, jetsubstructurehf::MCDetIdentifier<>); //!
DECLARE_SOA_TABLE(JetSubtructureHFMCGen, "AOD", "JETSUBHFMCP", jetsubstructurehf::HFJetId, jetsubstructurehf::Zg, jetsubstructurehf::Rg, jetsubstructurehf::Nsd, jetsubstructurehf::MCGenIdentifier<>); //!

namespace jetsubstructurehfoutput
{
DECLARE_SOA_COLUMN(JetPt, jetPt, float);                       //!
DECLARE_SOA_COLUMN(JetPhi, jetPhi, float);                     //!
DECLARE_SOA_COLUMN(JetEta, jetEta, float);                     //!
DECLARE_SOA_COLUMN(JetNConstituents, jetNConstituents, float); //!

DECLARE_SOA_COLUMN(CandPt, candPt, float);                  //!
DECLARE_SOA_COLUMN(CandPhi, candPhi, float);                //!
DECLARE_SOA_COLUMN(CandEta, candEta, float);                //!
DECLARE_SOA_COLUMN(CandY, candY, float);                    //!
DECLARE_SOA_COLUMN(CandInvMass, candInvMasss, float);       //!
DECLARE_SOA_COLUMN(CandBarInvMass, candBarInvMasss, float); //!
} // namespace jetsubstructurehfoutput
DECLARE_SOA_TABLE(JetSubstructureHFOutputData, "AOD", "JETSUBOutData", jetsubstructurehfoutput::JetPt, jetsubstructurehfoutput::JetPhi, jetsubstructurehfoutput::JetEta, jetsubstructurehfoutput::JetNConstituents,
                  jetsubstructurehfoutput::CandPt, jetsubstructurehfoutput::CandPhi, jetsubstructurehfoutput::CandEta, jetsubstructurehfoutput::CandY, jetsubstructurehfoutput::CandInvMass, jetsubstructurehfoutput::CandBarInvMass,
                  jetsubstructurehf::Zg, jetsubstructurehf::Rg, jetsubstructurehf::Nsd, jetsubstructurehf::DataIdentifier<>); //!

DECLARE_SOA_TABLE(JetSubstructureHFOutputMCDet, "AOD", "JETSUBOutMCD", jetsubstructurehfoutput::JetPt, jetsubstructurehfoutput::JetPhi, jetsubstructurehfoutput::JetEta, jetsubstructurehfoutput::JetNConstituents,
                  jetsubstructurehfoutput::CandPt, jetsubstructurehfoutput::CandPhi, jetsubstructurehfoutput::CandEta, jetsubstructurehfoutput::CandY, jetsubstructurehfoutput::CandInvMass, jetsubstructurehfoutput::CandBarInvMass,
                  jetsubstructurehf::Zg, jetsubstructurehf::Rg, jetsubstructurehf::Nsd, jetsubstructurehf::MCDetIdentifier<>); //!

DECLARE_SOA_TABLE(JetSubstructureHFOutputMCGen, "AOD", "JETSUBOutMCP", jetsubstructurehfoutput::JetPt, jetsubstructurehfoutput::JetPhi, jetsubstructurehfoutput::JetEta, jetsubstructurehfoutput::JetNConstituents,
                  jetsubstructurehfoutput::CandPt, jetsubstructurehfoutput::CandPhi, jetsubstructurehfoutput::CandEta, jetsubstructurehfoutput::CandY, jetsubstructurehfoutput::CandInvMass, jetsubstructurehfoutput::CandBarInvMass,
                  jetsubstructurehf::Zg, jetsubstructurehf::Rg, jetsubstructurehf::Nsd, jetsubstructurehf::MCGenIdentifier<>); //!

} // namespace o2::aod

#endif // PWGJE_DATAMODEL_JETSUBSTRUCTUREHF_H_
