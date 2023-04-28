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

#ifndef PWGJE_DATAMODEL_JETSUBSTRUCTURE_H_
#define PWGJE_DATAMODEL_JETSUBSTRUCTURE_H_

#include <cmath>
#include "Framework/AnalysisDataModel.h"
#include "PWGJE/DataModel/EMCALClusters.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGJE/DataModel/Jet.h"

using namespace o2::analysis;

// Defines the jet substrcuture table definition
#define JETSUBSTRUCTURE_TABLE_DEF(_jet_type_, _name_, _description_) \
  namespace _name_##substructure                                     \
  {                                                                  \
    DECLARE_SOA_INDEX_COLUMN(_jet_type_, jet);                       \
    DECLARE_SOA_DYNAMIC_COLUMN(Dummy##_jet_type_, dummy##_jet_type_, \
                               []() -> int { return 0; });           \
  }                                                                  \
  DECLARE_SOA_TABLE(_jet_type_##Substructure, "AOD", _description_, _name_##substructure::_jet_type_##Id, jetsubstructure::Zg, jetsubstructure::Rg, jetsubstructure::Nsd, _name_##substructure::Dummy##_jet_type_<>);

// Defines the jet substrcuture hf output table definition
#define JETSUBSTRUCTUREOUTPUT_TABLE_DEF(_jet_type_, _name_, _description_)                                                                                                                                     \
  namespace _name_##substructureoutput                                                                                                                                                                         \
  {                                                                                                                                                                                                            \
    DECLARE_SOA_DYNAMIC_COLUMN(Dummy##_jet_type_, dummy##_jet_type_,                                                                                                                                           \
                               []() -> int { return 0; });                                                                                                                                                     \
  }                                                                                                                                                                                                            \
  DECLARE_SOA_TABLE(_jet_type_##SubstructureOutput, "AOD", _description_, jetsubstructureoutput::JetPt, jetsubstructureoutput::JetPhi, jetsubstructureoutput::JetEta, jetsubstructureoutput::JetNConstituents, \
                    jetsubstructure::Zg, jetsubstructure::Rg, jetsubstructure::Nsd, _name_##substructureoutput::Dummy##_jet_type_<>);

// Defines the jet substrcuture hf output table definition
#define JETSUBSTRUCTUREOUTPUTHF_TABLE_DEF(_jet_type_, _name_, _description_)                                                                                                                                                            \
  namespace _name_##substructurehfoutput                                                                                                                                                                                                \
  {                                                                                                                                                                                                                                     \
    DECLARE_SOA_DYNAMIC_COLUMN(Dummy##_jet_type_, dummy##_jet_type_,                                                                                                                                                                    \
                               []() -> int { return 0; });                                                                                                                                                                              \
  }                                                                                                                                                                                                                                     \
  DECLARE_SOA_TABLE(_jet_type_##SubstructureOutput, "AOD", _description_, jetsubstructureoutput::JetPt, jetsubstructureoutput::JetPhi, jetsubstructureoutput::JetEta, jetsubstructureoutput::JetNConstituents,                          \
                    jetsubstructurehfoutput::CandPt, jetsubstructurehfoutput::CandPhi, jetsubstructurehfoutput::CandEta, jetsubstructurehfoutput::CandY, jetsubstructurehfoutput::CandInvMass, jetsubstructurehfoutput::CandBarInvMass, \
                    jetsubstructure::Zg, jetsubstructure::Rg, jetsubstructure::Nsd, _name_##substructurehfoutput::Dummy##_jet_type_<>);

namespace o2::aod
{
namespace jetsubstructure
{                                    //!
DECLARE_SOA_COLUMN(Zg, zg, float);   //!
DECLARE_SOA_COLUMN(Rg, rg, float);   //!
DECLARE_SOA_COLUMN(Nsd, nsd, float); //!
} // namespace jetsubstructure

namespace jetsubstructureoutput
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
} // namespace jetsubstructureoutput

namespace jetsubstructurehfoutput
{
DECLARE_SOA_COLUMN(CandPt, candPt, float);                  //!
DECLARE_SOA_COLUMN(CandPhi, candPhi, float);                //!
DECLARE_SOA_COLUMN(CandEta, candEta, float);                //!
DECLARE_SOA_COLUMN(CandY, candY, float);                    //!
DECLARE_SOA_COLUMN(CandInvMass, candInvMasss, float);       //!
DECLARE_SOA_COLUMN(CandBarInvMass, candBarInvMasss, float); //!
} // namespace jetsubstructurehfoutput

JETSUBSTRUCTURE_TABLE_DEF(ChargedJet, jet, "JETSS");
JETSUBSTRUCTUREOUTPUT_TABLE_DEF(Jet, jet, "JETSSO");

JETSUBSTRUCTURE_TABLE_DEF(ChargedMCDetectorLevelJet, mcdetectorleveljet, "JETMCDETSS");
JETSUBSTRUCTUREOUTPUT_TABLE_DEF(MCDetectorLevelJet, mcdetectorleveljet, "JETMCDETSSO");

JETSUBSTRUCTURE_TABLE_DEF(ChargedMCParticleLevelJet, mcparticleleveljet, "JETMCPARTSS");
JETSUBSTRUCTUREOUTPUT_TABLE_DEF(MCParticleLevelJet, mcparticleleveljet, "JETMCPARTSSO");

JETSUBSTRUCTURE_TABLE_DEF(D0ChargedJet, D0jet, "D0JETSS");
JETSUBSTRUCTUREOUTPUTHF_TABLE_DEF(D0ChargedJet, D0jet, "D0JETSSO");

JETSUBSTRUCTURE_TABLE_DEF(D0ChargedMCDetectorLevelJet, D0mcdjet, "D0JETMCDETSS");
JETSUBSTRUCTUREOUTPUTHF_TABLE_DEF(D0ChargedMCDetectorLevelJet, D0mcdjet, "D0JETMCDETSSO");

JETSUBSTRUCTURE_TABLE_DEF(D0ChargedMCParticleLevelJet, D0mcpjet, "D0JETMCPARTSS");
JETSUBSTRUCTUREOUTPUTHF_TABLE_DEF(D0ChargedMCParticleLevelJet, D0mcpjet, "D0JETMCPARTSSO");

JETSUBSTRUCTURE_TABLE_DEF(LcChargedJet, Lcjet, "LcJETSS");
JETSUBSTRUCTUREOUTPUTHF_TABLE_DEF(LcChargedJet, Lcjet, "LcJETSSO");

JETSUBSTRUCTURE_TABLE_DEF(LcChargedMCDetectorLevelJet, Lcmcdjet, "LcJETMCDETSS");
JETSUBSTRUCTUREOUTPUTHF_TABLE_DEF(LcChargedMCDetectorLevelJet, Lcmcdjet, "LcJETMCDETSSO");

JETSUBSTRUCTURE_TABLE_DEF(LcChargedMCParticleLevelJet, Lcmcpjet, "LcJETMCPARTSS");
JETSUBSTRUCTUREOUTPUTHF_TABLE_DEF(LcChargedMCParticleLevelJet, Lcmcpjet, "LcJETMCPARTSSO");

JETSUBSTRUCTURE_TABLE_DEF(BplusChargedJet, BPljet, "BPLJETSS");
JETSUBSTRUCTUREOUTPUTHF_TABLE_DEF(BplusChargedJet, BPljet, "BPLJETOut");

JETSUBSTRUCTURE_TABLE_DEF(BplusChargedMCDetectorLevelJet, BPlMCDjet, "BPLJETMCDETSS");
JETSUBSTRUCTUREOUTPUTHF_TABLE_DEF(BplusChargedMCDetectorLevelJet, BPlmcdjet, "BPLJETMCDETSSO");

JETSUBSTRUCTURE_TABLE_DEF(BplusChargedMCParticleLevelJet, BPlmcpjet, "BPLJETMCPARTSS");
JETSUBSTRUCTUREOUTPUTHF_TABLE_DEF(BplusChargedMCParticleLevelJet, BPlmcpjet, "BPLJETMCPARTSSO");

} // namespace o2::aod

#endif // PWGJE_DATAMODEL_JETSUBSTRUCTURE_H_
