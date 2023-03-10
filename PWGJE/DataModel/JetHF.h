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
/// \brief Table definitions for jets
///
/// Since the JE framework requires a set of nearly identical tables, most the tables are
/// generated via macros. Usually this would be avoided, but maintaining a collection of
/// (nearly) identical tables was judged to be more the larger maintenance burden.
///
/// \author Nima Zardoshti

#ifndef PWGJE_DATAMODEL_JETHF_H_
#define PWGJE_DATAMODEL_JETHF_H_

#include "PWGJE/DataModel/Jet.h"

namespace o2::aod
{

// HF jets
// D0 tagged jets
JET_TABLE_DEF(Collision, D0Jet, D0jet, "D0JET");
using D0Jet = D0Jets::iterator;
using MatchedD0Jet = MatchedD0Jets::iterator;
JET_CONSTITUENTS_ARRAY_TABLE_DEF(D0Jet, D0jet, "D0", Track, HfCand2Prong);
using D0JetConstituent = D0JetConstituents::iterator;
JET_CONSTITUENTS_SUB_TABLE_DEF(D0Jet, D0jet, "D0");
using D0JetConstituentSub = D0JetConstituentsSub::iterator;

// Lc tagged jets
JET_TABLE_DEF(Collision, LcJet, Lcjet, "LcJET");
using LcJet = LcJets::iterator;
using MatchedLcJet = MatchedLcJets::iterator;
JET_CONSTITUENTS_ARRAY_TABLE_DEF(LcJet, Lcjet, "Lc", Track, HfCand3Prong);
using LcJetConstituent = LcJetConstituents::iterator;
JET_CONSTITUENTS_SUB_TABLE_DEF(LcJet, Lcjet, "Lc");
using LcJetConstituentSub = LcJetConstituentsSub::iterator;

// B+ tagged jets
JET_TABLE_DEF(Collision, BPlusJet, BPlusjet, "BPlJET");
using BPlusJet = BPlusJets::iterator;
using MatchedBPlusJet = MatchedBPlusJets::iterator;
JET_CONSTITUENTS_ARRAY_TABLE_DEF(BPlusJet, BPlusjet, "BPl", Track, HfCandBplus);
using BPlusJetConstituent = BPlusJetConstituents::iterator;
JET_CONSTITUENTS_SUB_TABLE_DEF(BPlusJet, BPlusjet, "BPl");
using BPlusJetConstituentSub = BPlusJetConstituentsSub::iterator;

// HF jets (MC detector level)
// D0 tagged jets
JET_TABLE_DEF(Collision, MCDetectorLevelD0Jet, mcdetectorlevelD0jet, "D0JETMCD");
using MCDetectorLevelD0Jet = MCDetectorLevelD0Jets::iterator;
using MatchedMCDetectorLevelD0Jet = MatchedMCDetectorLevelD0Jets::iterator;
JET_CONSTITUENTS_ARRAY_TABLE_DEF(MCDetectorLevelD0Jet, mcdetectorlevelD0jet, "D0MCD", Track, HfCand2Prong);
using MCDetectorLevelD0JetConstituent = MCDetectorLevelD0JetConstituents::iterator;
JET_CONSTITUENTS_SUB_TABLE_DEF(MCDetectorLevelD0Jet, mcdetectorlevelD0jet, "D0MCD");
using MCDetectorLevelD0JetConstituentSub = MCDetectorLevelD0JetConstituentsSub::iterator;

// Lc tagged jets
JET_TABLE_DEF(Collision, MCDetectorLevelLcJet, mcdetectorlevelLcjet, "LcJETMCD");
using MCDetectorLevelLcJet = MCDetectorLevelLcJets::iterator;
using MatchedMCDetectorLevelLcJet = MatchedMCDetectorLevelLcJets::iterator;
JET_CONSTITUENTS_ARRAY_TABLE_DEF(MCDetectorLevelLcJet, mcdetectorlevelLcjet, "LcMCD", Track, HfCand3Prong);
using MCDetectorLevelLcJetConstituent = MCDetectorLevelLcJetConstituents::iterator;
JET_CONSTITUENTS_SUB_TABLE_DEF(MCDetectorLevelLcJet, mcdetectorlevelLcjet, "LcMCD");
using MCDetectorLevelLcJetConstituentSub = MCDetectorLevelLcJetConstituentsSub::iterator;

// B+ tagged jets
JET_TABLE_DEF(Collision, MCDetectorLevelBPlusJet, mcdetectorlevelBPlusjet, "BPlJETMCD");
using MCDetectorLevelBPlusJet = MCDetectorLevelBPlusJets::iterator;
using MatchedMCDetectorLevelBPlusJet = MatchedMCDetectorLevelBPlusJets::iterator;
JET_CONSTITUENTS_ARRAY_TABLE_DEF(MCDetectorLevelBPlusJet, mcdetectorlevelBPlusjet, "BPlMCD", Track, HfCandBplus);
using MCDetectorLevelBPlusJetConstituent = MCDetectorLevelBPlusJetConstituents::iterator;
JET_CONSTITUENTS_SUB_TABLE_DEF(MCDetectorLevelBPlusJet, mcdetectorlevelBPlusjet, "BPlMCD");
using MCDetectorLevelBPlusJetConstituentSub = MCDetectorLevelBPlusJetConstituentsSub::iterator;

// HF jets (MC particle level)
// D0 tagged jets
JET_TABLE_DEF(McCollision, MCParticleLevelD0Jet, mcparticlelevelD0jet, "D0JETMCP");
using MCParticleLevelD0Jet = MCParticleLevelD0Jets::iterator;
using MatchedMCParticleLevelD0Jet = MatchedMCParticleLevelD0Jets::iterator;
JET_CONSTITUENTS_ARRAY_TABLE_DEF(MCParticleLevelD0Jet, mcparticlelevelD0jet, "D0MCP", McParticle, McParticles);
using MCParticleLevelD0JetConstituent = MCParticleLevelD0JetConstituents::iterator;
JET_CONSTITUENTS_SUB_TABLE_DEF(MCParticleLevelD0Jet, mcparticlelevelD0jet, "D0MCP");
using MCParticleLevelD0JetConstituentSub = MCParticleLevelD0JetConstituentsSub::iterator;

// Lc tagged jets
JET_TABLE_DEF(McCollision, MCParticleLevelLcJet, mcparticlelevelLcjet, "LcJETMCP");
using MCParticleLevelLcJet = MCParticleLevelLcJets::iterator;
using MatchedMCParticleLevelLcJet = MatchedMCParticleLevelLcJets::iterator;
JET_CONSTITUENTS_ARRAY_TABLE_DEF(MCParticleLevelLcJet, mcparticlelevelLcjet, "LcMCP", McParticle, McParticles);
using MCParticleLevelLcJetConstituent = MCParticleLevelLcJetConstituents::iterator;
JET_CONSTITUENTS_SUB_TABLE_DEF(MCParticleLevelLcJet, mcparticlelevelLcjet, "LcMCP");
using MCParticleLevelLcJetConstituentSub = MCParticleLevelLcJetConstituentsSub::iterator;

// B+ tagged jets
JET_TABLE_DEF(McCollision, MCParticleLevelBPlusJet, mcparticlelevelBPlusjet, "BPlJETMCP");
using MCParticleLevelBPlusJet = MCParticleLevelBPlusJets::iterator;
using MatchedMCParticleLevelBPlusJet = MatchedMCParticleLevelBPlusJets::iterator;
JET_CONSTITUENTS_ARRAY_TABLE_DEF(MCParticleLevelBPlusJet, mcparticlelevelBPlusjet, "BPlMCP", McParticle, McParticles);
using MCParticleLevelBPlusJetConstituent = MCParticleLevelBPlusJetConstituents::iterator;
JET_CONSTITUENTS_SUB_TABLE_DEF(MCParticleLevelBPlusJet, mcparticlelevelBPlusjet, "BPlMCP");
using MCParticleLevelBPlusJetConstituentSub = MCParticleLevelBPlusJetConstituentsSub::iterator;

namespace mcdetectorlevelD0jetmatching2
{
DECLARE_SOA_INDEX_COLUMN(MCDetectorLevelD0Jet, matchedJet);
}
namespace mcparticlelevelD0jetmatching2
{
DECLARE_SOA_INDEX_COLUMN(MCParticleLevelD0Jet, matchedJet);
}
DECLARE_SOA_TABLE(MatchedMCParticleDetectorLevelD0Jets, "AOD", "D0JETMCPDMATCH", mcdetectorlevelD0jetmatching2::MCDetectorLevelD0JetId);
DECLARE_SOA_TABLE(MatchedMCDetectorParticleLevelD0Jets, "AOD", "D0JETMCDPMATCH", mcparticlelevelD0jetmatching2::MCParticleLevelD0JetId);

} // namespace o2::aod

#endif // PWGJE_DATAMODEL_JETHF_H_
