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
#include "Framework/AnalysisDataModel.h"
#include "PWGJE/DataModel/EMCALClusters.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGJE/DataModel/Jet.h"

using namespace o2::analysis;

namespace o2::aod
{

namespace jetsubstructure
{                                    //!
DECLARE_SOA_COLUMN(Zg, zg, float);   //!
DECLARE_SOA_COLUMN(Rg, rg, float);   //!
DECLARE_SOA_COLUMN(Nsd, nsd, float); //!
} // namespace jetsubstructure

namespace jetoutput
{
DECLARE_SOA_INDEX_COLUMN(Collision, collision);                //!
DECLARE_SOA_INDEX_COLUMN(McCollision, mcCollision);            //!
DECLARE_SOA_COLUMN(JetPt, jetPt, float);                       //!
DECLARE_SOA_COLUMN(JetPhi, jetPhi, float);                     //!
DECLARE_SOA_COLUMN(JetEta, jetEta, float);                     //!
DECLARE_SOA_COLUMN(JetR, jetR, float);                         //!
DECLARE_SOA_COLUMN(JetNConstituents, jetNConstituents, int);   //!

} // namespace jetoutput

// Defines the jet substrcuture table definition
#define JETSUBSTRUCTURE_TABLE_DEF(_collision_type_, _jet_type_, _matched_jet_type_, _cand_type_, _name_, _description_)                                                                                                                                                                                                                                                                                                     \
  namespace _name_##substructure                                                                                                                                                                                                                                                                                                                                                                                            \
  {                                                                                                                                                                                                                                                                                                                                                                                                                         \
    DECLARE_SOA_INDEX_COLUMN(_jet_type_, jet);                                                                                                                                                                                                                                                                                                                                                                              \
    DECLARE_SOA_INDEX_COLUMN_FULL(Candidate, candidate, int, _cand_type_, "_0");                                                                                                                                                                                                                                                                                                                                            \
    DECLARE_SOA_DYNAMIC_COLUMN(Dummy##_jet_type_, dummy##_jet_type_, []() -> int { return 0; });                                                                                                                                                                                                                                                                                                                            \
  }                                                                                                                                                                                                                                                                                                                                                                                                                         \
                                                                                                                                                                                                                                                                                                                                                                                                                            \
  namespace _name_##geomatched                                                                                                                                                                                                                                                                                                                                                                                              \
  {                                                                                                                                                                                                                                                                                                                                                                                                                         \
    DECLARE_SOA_ARRAY_INDEX_COLUMN_FULL(_matched_jet_type_, matchedJetGeo, int32_t, _matched_jet_type_##s, "_geo");                                                                                                                                                                                                                                                                                                         \
  }                                                                                                                                                                                                                                                                                                                                                                                                                         \
                                                                                                                                                                                                                                                                                                                                                                                                                            \
  namespace _name_##ptmatched                                                                                                                                                                                                                                                                                                                                                                                               \
  {                                                                                                                                                                                                                                                                                                                                                                                                                         \
    DECLARE_SOA_ARRAY_INDEX_COLUMN_FULL(_matched_jet_type_, matchedJetPt, int32_t, _matched_jet_type_##s, "_pt");                                                                                                                                                                                                                                                                                                           \
  }                                                                                                                                                                                                                                                                                                                                                                                                                         \
                                                                                                                                                                                                                                                                                                                                                                                                                            \
  namespace _name_##candmatched                                                                                                                                                                                                                                                                                                                                                                                             \
  {                                                                                                                                                                                                                                                                                                                                                                                                                         \
    DECLARE_SOA_ARRAY_INDEX_COLUMN_FULL(_matched_jet_type_, matchedJetCand, int32_t, _matched_jet_type_##s, "_hf");                                                                                                                                                                                                                                                                                                         \
  }                                                                                                                                                                                                                                                                                                                                                                                                                         \
                                                                                                                                                                                                                                                                                                                                                                                                                            \
  DECLARE_SOA_TABLE(_jet_type_##Substructures, "AOD", _description_ "SS", jetsubstructure::Zg, jetsubstructure::Rg, jetsubstructure::Nsd, _name_##substructure::Dummy##_jet_type_<>);                                                                                                                                                                                                                                       \
  DECLARE_SOA_TABLE(_jet_type_##Output, "AOD", _description_ "O", jetoutput::_collision_type_##Id, _name_##substructure::_jet_type_##Id, _name_##substructure::Candidate##Id, _name_##geomatched::_matched_jet_type_##Ids, _name_##ptmatched::_matched_jet_type_##Ids, _name_##candmatched::_matched_jet_type_##Ids, jetoutput::JetPt, jetoutput::JetPhi, jetoutput::JetEta, jetoutput::JetR, jetoutput::JetNConstituents); \
  DECLARE_SOA_TABLE(_jet_type_##SubstructureOutput, "AOD", _description_ "SSO", _name_##substructure::_jet_type_##Id, jetsubstructure::Zg, jetsubstructure::Rg, jetsubstructure::Nsd);

#define JETSUBSTRUCTURE_TABLES_DEF(_jet_type_, _cand_type_, _description_)                                                                                               \
  JETSUBSTRUCTURE_TABLE_DEF(Collision, _jet_type_##Jet, _jet_type_##Jet, _cand_type_, _jet_type_##jet, _description_)                                                    \
  JETSUBSTRUCTURE_TABLE_DEF(Collision, _jet_type_##MCDetectorLevelJet, _jet_type_##MCParticleLevelJet, _cand_type_, _jet_type_##mcdetectorleveljet, _description_ "MCD") \
  JETSUBSTRUCTURE_TABLE_DEF(McCollision, _jet_type_##MCParticleLevelJet, _jet_type_##MCDetectorLevelJet, McParticles, _jet_type_##mcparticleleveljet, _description_ "MCP")

JETSUBSTRUCTURE_TABLES_DEF(Charged, HfCand2Prong, "C");
JETSUBSTRUCTURE_TABLES_DEF(D0Charged, HfCand2Prong, "D0");
JETSUBSTRUCTURE_TABLES_DEF(LcCharged, HfCand3Prong, "Lc");
JETSUBSTRUCTURE_TABLES_DEF(BplusCharged, HfCandBplus, "BPL");

} // namespace o2::aod

#endif // PWGJE_DATAMODEL_JETSUBSTRUCTURE_H_
