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
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>

#ifndef PWGJE_DATAMODEL_JETSUBSTRUCTURE_H_
#define PWGJE_DATAMODEL_JETSUBSTRUCTURE_H_

#include <cmath>
#include <vector>
#include "Framework/AnalysisDataModel.h"
#include "PWGJE/DataModel/EMCALClusters.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGJE/DataModel/Jet.h"

using namespace o2::analysis;

namespace o2::aod
{

namespace jetcollision
{                                                  //!
DECLARE_SOA_COLUMN(PosZ, posZ, float);             //!
DECLARE_SOA_COLUMN(Centrality, centrality, float); //!
DECLARE_SOA_COLUMN(EventSel, eventSel, uint8_t);   //!
} // namespace jetcollision

namespace jetsubstructure
{                                                                   //!
DECLARE_SOA_COLUMN(EnergyMother, energyMother, std::vector<float>); //!
DECLARE_SOA_COLUMN(PtLeading, ptLeading, std::vector<float>);       //!
DECLARE_SOA_COLUMN(PtSubLeading, ptSubLeading, std::vector<float>); //!
DECLARE_SOA_COLUMN(Theta, theta, std::vector<float>);               //!
} // namespace jetsubstructure

namespace jetoutput
{
DECLARE_SOA_COLUMN(JetPt, jetPt, float);                     //!
DECLARE_SOA_COLUMN(JetPhi, jetPhi, float);                   //!
DECLARE_SOA_COLUMN(JetEta, jetEta, float);                   //!
DECLARE_SOA_COLUMN(JetR, jetR, float);                       //!
DECLARE_SOA_COLUMN(JetNConstituents, jetNConstituents, int); //!

} // namespace jetoutput

// Defines the jet substrcuture table definition
#define JETSUBSTRUCTURE_TABLE_DEF(_collision_type_, _table_type_, _jet_type_, _matched_jet_type_, _cand_type_, _name_, _description_)                                                                                                              \
                                                                                                                                                                                                                                                   \
  namespace _name_##collisionoutput                                                                                                                                                                                                                \
  {                                                                                                                                                                                                                                                \
    DECLARE_SOA_DYNAMIC_COLUMN(Dummy##_collision_type_, dummy##_collision_type_, []() -> int { return 0; });                                                                                                                                       \
  }                                                                                                                                                                                                                                                \
                                                                                                                                                                                                                                                   \
  DECLARE_SOA_TABLE(_table_type_##COs, "AOD", _description_ "CO", jetcollision::PosZ, jetcollision::Centrality, jetcollision::EventSel, _name_##collisionoutput::Dummy##_collision_type_<>);                                                       \
  using _table_type_##CO = _table_type_##COs::iterator;                                                                                                                                                                                            \
                                                                                                                                                                                                                                                   \
  namespace _name_##geomatched                                                                                                                                                                                                                     \
  {                                                                                                                                                                                                                                                \
    DECLARE_SOA_ARRAY_INDEX_COLUMN_FULL(_matched_jet_type_, matchedJetGeo, int32_t, _matched_jet_type_##s, "_geo");                                                                                                                                \
  }                                                                                                                                                                                                                                                \
                                                                                                                                                                                                                                                   \
  namespace _name_##ptmatched                                                                                                                                                                                                                      \
  {                                                                                                                                                                                                                                                \
    DECLARE_SOA_ARRAY_INDEX_COLUMN_FULL(_matched_jet_type_, matchedJetPt, int32_t, _matched_jet_type_##s, "_pt");                                                                                                                                  \
  }                                                                                                                                                                                                                                                \
                                                                                                                                                                                                                                                   \
  namespace _name_##candmatched                                                                                                                                                                                                                    \
  {                                                                                                                                                                                                                                                \
    DECLARE_SOA_ARRAY_INDEX_COLUMN_FULL(_matched_jet_type_, matchedJetCand, int32_t, _matched_jet_type_##s, "_hf");                                                                                                                                \
  }                                                                                                                                                                                                                                                \
  namespace _name_##jetoutput                                                                                                                                                                                                                      \
  {                                                                                                                                                                                                                                                \
    DECLARE_SOA_INDEX_COLUMN(_table_type_##CO, collision);                                                                                                                                                                                         \
    DECLARE_SOA_INDEX_COLUMN_FULL(Candidate, candidate, int, _cand_type_, "_0");                                                                                                                                                                   \
  }                                                                                                                                                                                                                                                \
  DECLARE_SOA_TABLE(_table_type_##Os, "AOD", _description_ "O", _name_##jetoutput::_table_type_##COId, _name_##jetoutput::CandidateId, jetoutput::JetPt, jetoutput::JetPhi, jetoutput::JetEta, jetoutput::JetR, jetoutput::JetNConstituents);      \
  using _table_type_##O = _table_type_##Os::iterator;                                                                                                                                                                                              \
  namespace _name_##substructure                                                                                                                                                                                                                   \
  {                                                                                                                                                                                                                                                \
    DECLARE_SOA_INDEX_COLUMN(_table_type_##O, outputTable);                                                                                                                                                                                        \
    DECLARE_SOA_DYNAMIC_COLUMN(Dummy##_jet_type_, dummy##_jet_type_, []() -> int { return 0; });                                                                                                                                                   \
  }                                                                                                                                                                                                                                                \
                                                                                                                                                                                                                                                   \
  DECLARE_SOA_TABLE(_table_type_##SSs, "AOD", _description_ "SS", jetsubstructure::EnergyMother, jetsubstructure::PtLeading, jetsubstructure::PtSubLeading, jetsubstructure::Theta, _name_##substructure::Dummy##_jet_type_<>);                    \
  DECLARE_SOA_TABLE(_table_type_##SSOs, "AOD", _description_ "SSO", _name_##substructure::_table_type_##OId, jetsubstructure::EnergyMother, jetsubstructure::PtLeading, jetsubstructure::PtSubLeading, jetsubstructure::Theta);                    \
  DECLARE_SOA_TABLE(_table_type_##MOs, "AOD", _description_ "MO", _name_##substructure::_table_type_##OId, _name_##geomatched::_matched_jet_type_##Ids, _name_##ptmatched::_matched_jet_type_##Ids, _name_##candmatched::_matched_jet_type_##Ids); \
                                                                                                                                                                                                                                                   \
  using _table_type_##O = _table_type_##Os::iterator;                                                                                                                                                                                              \
  using _table_type_##SSO = _table_type_##SSOs::iterator;                                                                                                                                                                                          \
  using _table_type_##MO = _table_type_##MOs::iterator;

#define JETSUBSTRUCTURE_TABLES_DEF(table_type_, _jet_type_, _cand_type_, _hfparticle_type_, _description_)                                                                                                \
  JETSUBSTRUCTURE_TABLE_DEF(JCollision, table_type_##Jet, _jet_type_##Jet, _jet_type_##EventWiseSubtractedJet, _cand_type_, _jet_type_##jet, _description_ "JET")                                         \
  JETSUBSTRUCTURE_TABLE_DEF(JCollision, table_type_##MCDJet, _jet_type_##MCDetectorLevelJet, _jet_type_##MCParticleLevelJet, _cand_type_, _jet_type_##mcdetectorleveljet, _description_ "MCDJET")         \
  JETSUBSTRUCTURE_TABLE_DEF(JMcCollision, table_type_##MCPJet, _jet_type_##MCParticleLevelJet, _jet_type_##MCDetectorLevelJet, _hfparticle_type_, _jet_type_##mcparticleleveljet, _description_ "MCPJET") \
  JETSUBSTRUCTURE_TABLE_DEF(JCollision, table_type_##EWSJet, _jet_type_##EventWiseSubtractedJet, _jet_type_##Jet, _cand_type_, _jet_type_##eventwisesubtractedjet, _description_ "EWSJET")

JETSUBSTRUCTURE_TABLES_DEF(C, Charged, HfD0Bases, HfD0PBases, "C");
JETSUBSTRUCTURE_TABLES_DEF(D0C, D0Charged, HfD0Bases, HfD0PBases, "D0C");
JETSUBSTRUCTURE_TABLES_DEF(LcC, LcCharged, HfD0Bases, HfD0PBases, "LCC");
JETSUBSTRUCTURE_TABLES_DEF(BplusC, BplusCharged, HfD0Bases, HfD0PBases, "BPLUSC");

} // namespace o2::aod

#endif // PWGJE_DATAMODEL_JETSUBSTRUCTURE_H_
