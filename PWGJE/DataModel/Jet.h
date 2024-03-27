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
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>
/// \author Jochen Klein <jochen.klein@cern.ch>
/// \author Raymond Ehlers <raymond.ehlers@cern.ch>, ORNL

#ifndef PWGJE_DATAMODEL_JET_H_
#define PWGJE_DATAMODEL_JET_H_

#include <cmath>
#include "Framework/AnalysisDataModel.h"
#include "PWGJE/DataModel/EMCALClusters.h"
#include "PWGJE/DataModel/JetReducedData.h"
#include "PWGJE/DataModel/JetReducedDataHF.h"
#include "PWGJE/DataModel/JetReducedDataLF.h"
#include "PWGJE/DataModel/JetSubtraction.h"

#include "PWGHF/DataModel/DerivedTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

#include "PWGLF/DataModel/LFStrangenessTables.h"

namespace o2::aod
{

namespace jet
{
DECLARE_SOA_INDEX_COLUMN(JCollision, collision);
DECLARE_SOA_INDEX_COLUMN(JMcCollision, mcCollision);
DECLARE_SOA_COLUMN(Pt, pt, float);         //!
DECLARE_SOA_COLUMN(Eta, eta, float);       //!
DECLARE_SOA_COLUMN(Phi, phi, float);       //!
DECLARE_SOA_COLUMN(Y, y, float);           //!
DECLARE_SOA_COLUMN(Energy, energy, float); //!
DECLARE_SOA_COLUMN(Mass, mass, float);     //!
DECLARE_SOA_COLUMN(Area, area, float);     //!
DECLARE_SOA_COLUMN(R, r, int);             //!
DECLARE_SOA_DYNAMIC_COLUMN(Px, px,         //!
                           [](float pt, float phi) -> float { return pt * std::cos(phi); });
DECLARE_SOA_DYNAMIC_COLUMN(Py, py, //!
                           [](float pt, float phi) -> float { return pt * std::sin(phi); });
DECLARE_SOA_DYNAMIC_COLUMN(Pz, pz, //!
                           [](float pt, float eta) -> float { return pt * std::sinh(eta); });
DECLARE_SOA_DYNAMIC_COLUMN(P, p, //! absolute p
                           [](float pt, float eta) -> float { return pt * std::cosh(eta); });
} // namespace jet
} // namespace o2::aod

// Defines the jet table definition
#define DECLARE_JET_TABLE(_collision_name_, _jet_type_, _name_, _description_)                      \
  namespace _name_##util                                                                            \
  {                                                                                                 \
    DECLARE_SOA_DYNAMIC_COLUMN(Dummy##_jet_type_##s, dummy##_jet_type##s,                           \
                               []() -> int { return 0; });                                          \
  }                                                                                                 \
  DECLARE_SOA_TABLE(_jet_type_##s, "AOD", _description_,                                            \
                    o2::soa::Index<>,                                                               \
                    jet::_collision_name_##Id,                                                      \
                    jet::Pt,                                                                        \
                    jet::Eta,                                                                       \
                    jet::Phi,                                                                       \
                    jet::Energy,                                                                    \
                    jet::Y,                                                                         \
                    jet::Mass,                                                                      \
                    jet::Area,                                                                      \
                    jet::R,                                                                         \
                    jet::Px<jet::Pt, jet::Phi>,                                                     \
                    jet::Py<jet::Pt, jet::Phi>,                                                     \
                    jet::Pz<jet::Pt, jet::Eta>,                                                     \
                    jet::P<jet::Pt, jet::Eta>,                                                      \
                    _name_##util::Dummy##_jet_type_##s<>);                                          \
  namespace _name_##matchingGeo                                                                     \
  {                                                                                                 \
    DECLARE_SOA_ARRAY_INDEX_COLUMN_FULL(_jet_type_, matchedJetGeo, int32_t, _jet_type_##s, "_geo"); \
  }                                                                                                 \
  namespace _name_##matchingPt                                                                      \
  {                                                                                                 \
    DECLARE_SOA_ARRAY_INDEX_COLUMN_FULL(_jet_type_, matchedJetPt, int32_t, _jet_type_##s, "_pt");   \
  }                                                                                                 \
  namespace _name_##matchingCand                                                                    \
  {                                                                                                 \
    DECLARE_SOA_ARRAY_INDEX_COLUMN_FULL(_jet_type_, matchedJetCand, int32_t, _jet_type_##s, "_hf"); \
  }

#define DECLARE_CONSTITUENTS_TABLE(_jet_type_, _name_, _Description_, _track_type_, _cand_type_)      \
  namespace _name_##constituents                                                                      \
  {                                                                                                   \
    DECLARE_SOA_INDEX_COLUMN(_jet_type_, jet);                                                        \
    DECLARE_SOA_ARRAY_INDEX_COLUMN(_track_type_, tracks);                                             \
    DECLARE_SOA_ARRAY_INDEX_COLUMN(JCluster, clusters);                                               \
    DECLARE_SOA_ARRAY_INDEX_COLUMN_FULL(HfCandidates, hfcandidates, int32_t, _cand_type_, "_hfcand"); \
    DECLARE_SOA_ARRAY_INDEX_COLUMN(V0Data, v0candidates);                                             \
  }                                                                                                   \
  DECLARE_SOA_TABLE(_jet_type_##Constituents, "AOD", _Description_ "C",                               \
                    _name_##constituents::_jet_type_##Id,                                             \
                    _name_##constituents::_track_type_##Ids,                                          \
                    _name_##constituents::JClusterIds,                                                \
                    _name_##constituents::HfCandidatesIds,                                            \
                    _name_##constituents::V0DataIds);

// combine definition of tables for jets, constituents
#define DECLARE_JET_TABLES(_collision_name_, _jet_type_, _track_type_, _hfcand_type_, _description_)        \
  DECLARE_JET_TABLE(_collision_name_, _jet_type_##Jet, _jet_type_##jet, _description_);                     \
  using _jet_type_##Jet = _jet_type_##Jet##s::iterator;                                                     \
  DECLARE_CONSTITUENTS_TABLE(_jet_type_##Jet, _jet_type_##jet, _description_, _track_type_, _hfcand_type_); \
  using _jet_type_##Jet##Constituent = _jet_type_##Jet##Constituents::iterator;

#define DECLARE_JETMATCHING_TABLE(_jet_type_base_, _jet_type_tag_, _description_)                 \
  namespace _jet_type_base_##jetsmatchedto##_jet_type_tag_                                        \
  {                                                                                               \
    DECLARE_SOA_DYNAMIC_COLUMN(Dummy##_jet_type_base_##s, dummy##_jet_type_base_##s,              \
                               []() -> int { return 0; });                                        \
  }                                                                                               \
  DECLARE_SOA_TABLE(_jet_type_base_##JetsMatchedTo##_jet_type_tag_##Jets, "AOD", _description_,   \
                    _jet_type_tag_##jetmatchingGeo::_jet_type_tag_##JetIds,                       \
                    _jet_type_tag_##jetmatchingPt::_jet_type_tag_##JetIds,                        \
                    _jet_type_tag_##jetmatchingCand::_jet_type_tag_##JetIds,                      \
                    _jet_type_base_##jetsmatchedto##_jet_type_tag_::Dummy##_jet_type_base_##s<>); \
  using _jet_type_base_##JetsMatchedTo##_jet_type_tag_##Jet = _jet_type_base_##JetsMatchedTo##_jet_type_tag_##Jets::iterator;

#define DECLARE_MCEVENTWEIGHT_TABLE(_jet_type_, _name_, _description_) \
  namespace _name_##eventweights                                       \
  {                                                                    \
    DECLARE_SOA_INDEX_COLUMN(_jet_type_##Jet, jet);                    \
    DECLARE_SOA_COLUMN(EventWeight, eventWeight, float);               \
  }                                                                    \
  DECLARE_SOA_TABLE(_jet_type_##JetEventWeights, "AOD", _description_, \
                    _name_##eventweights::_jet_type_##JetId,           \
                    _name_##eventweights::EventWeight);                \
  using _jet_type_##JetEventWeight = _jet_type_##JetEventWeights::iterator;

// generate tables for data-, detector- and particle-level jets
#define DECLARE_JET_TABLES_LEVELS(_jet_type_, _subtracted_track_type_, _hfcand_type_, _hfparticle_type_, _shortname_)                                 \
  DECLARE_JET_TABLES(JCollision, _jet_type_, JTrack, _hfcand_type_, _shortname_ "JET")                                                                \
  DECLARE_JET_TABLES(JCollision, _jet_type_##MCDetectorLevel, JTrack, _hfcand_type_, _shortname_ "DJET")                                              \
  DECLARE_JET_TABLES(JMcCollision, _jet_type_##MCParticleLevel, JMcParticle, _hfparticle_type_, _shortname_ "PJET")                                   \
  DECLARE_JETMATCHING_TABLE(_jet_type_##MCDetectorLevel, _jet_type_##MCParticleLevel, _shortname_ "JETD2P")                                           \
  DECLARE_JETMATCHING_TABLE(_jet_type_##MCParticleLevel, _jet_type_##MCDetectorLevel, _shortname_ "JETP2D")                                           \
  DECLARE_MCEVENTWEIGHT_TABLE(_jet_type_##MCDetectorLevel, _jet_type_##MCDetectorLevel, _shortname_ "DJETMW")                                         \
  DECLARE_MCEVENTWEIGHT_TABLE(_jet_type_##MCParticleLevel, _jet_type_##MCParticleLevel, _shortname_ "PETMPW")                                         \
  DECLARE_JET_TABLES(JCollision, _jet_type_##EventWiseSubtracted, _subtracted_track_type_, _hfcand_type_, _shortname_ "JETEWS")                       \
  DECLARE_JETMATCHING_TABLE(_jet_type_, _jet_type_##EventWiseSubtracted, _shortname_ "JET2EWS")                                                       \
  DECLARE_JETMATCHING_TABLE(_jet_type_##EventWiseSubtracted, _jet_type_, _shortname_ "JETEWS2")                                                       \
  DECLARE_JET_TABLES(JCollision, _jet_type_##MCDetectorLevelEventWiseSubtracted, _subtracted_track_type_, _hfcand_type_, _shortname_ "DJETEWS")       \
  DECLARE_MCEVENTWEIGHT_TABLE(_jet_type_##MCDetectorLevelEventWiseSubtracted, _jet_type_##MCDetectorLevelEventWiseSubtracted, _shortname_ "DJETEWSW") \
  DECLARE_JETMATCHING_TABLE(_jet_type_##MCDetectorLevel, _jet_type_##MCDetectorLevelEventWiseSubtracted, _shortname_ "DJET2DEWS")                     \
  DECLARE_JETMATCHING_TABLE(_jet_type_##MCDetectorLevelEventWiseSubtracted, _jet_type_##MCDetectorLevel, _shortname_ "JETDEWS2D")                     \
  DECLARE_JET_TABLES(JMcCollision, _jet_type_##MCParticleLevelEventWiseSubtracted, _subtracted_track_type_, _hfparticle_type_, _shortname_ "PJETEWS")
namespace o2::aod
{
DECLARE_JET_TABLES_LEVELS(Charged, JTrackSub, HfD0Bases, HfD0PBases, "C");
DECLARE_JET_TABLES_LEVELS(Full, JTrackSub, HfD0Bases, HfD0PBases, "F");
DECLARE_JET_TABLES_LEVELS(Neutral, JTrackSub, HfD0Bases, HfD0PBases, "N");
DECLARE_JET_TABLES_LEVELS(D0Charged, JTrackD0Sub, HfD0Bases, HfD0PBases, "D0");
DECLARE_JET_TABLES_LEVELS(LcCharged, JTrackLcSub, Hf3PBases, Hf3PPBases, "Lc");
DECLARE_JET_TABLES_LEVELS(BplusCharged, JTrackBplusSub, HfCandBplus, HfD0PBases, "BPl");
DECLARE_JET_TABLES_LEVELS(V0Charged, JTrackSub, V0Cores, V0Cores, "V0");

} // namespace o2::aod

using JetCollisions = o2::aod::JCollisions;
using JetCollision = JetCollisions::iterator;
using JetCollisionsMCD = o2::soa::Join<JetCollisions, o2::aod::JMcCollisionLbs>;
using JetTracks = o2::aod::JTracks;
using JetTracksMCD = o2::soa::Join<JetTracks, o2::aod::JMcTrackLbs>;
using JetTracksSub = o2::aod::JTrackSubs;
using JetClusters = o2::aod::JClusters;

using JetMcCollisions = o2::aod::JMcCollisions;
using JetMcCollision = JetMcCollisions::iterator;
using JetParticles = o2::aod::JMcParticles;

using CandidatesD0MCP = o2::soa::Join<o2::aod::HfD0PBases, o2::aod::JD0PIds>;
using CandidatesLcMCP = o2::soa::Join<o2::aod::Hf3PPBases, o2::aod::JLcPIds>;
using CandidatesBplusMCP = o2::soa::Join<o2::aod::JMcParticles, o2::aod::HfCandBplusMcGen>;

using CandidatesD0Data = o2::soa::Join<o2::aod::HfD0Bases, o2::aod::HfD0Pars, o2::aod::HfD0ParEs, o2::aod::HfD0Sels, o2::aod::HfD0Mls, o2::aod::JD0Ids>;
using CandidatesD0MCD = o2::soa::Join<o2::aod::HfD0Bases, o2::aod::HfD0Pars, o2::aod::HfD0ParEs, o2::aod::HfD0Sels, o2::aod::HfD0Mls, o2::aod::HfD0Mcs, o2::aod::JD0Ids>;
using JetTracksSubD0 = o2::aod::JTrackD0Subs;

using CandidatesBplusData = o2::soa::Join<o2::aod::HfCandBplus, o2::aod::HfSelBplusToD0Pi>;
using CandidatesBplusMCD = o2::soa::Join<o2::aod::HfCandBplus, o2::aod::HfSelBplusToD0Pi, o2::aod::HfCandBplusMcRec>;
using JetTracksSubBplus = o2::aod::JTrackBplusSubs;

using CandidatesLcData = o2::soa::Join<o2::aod::Hf3PBases, o2::aod::Hf3PPars, o2::aod::Hf3PParEs, o2::aod::Hf3PSels, o2::aod::Hf3PMls, o2::aod::JLcIds>;
using CandidatesLcMCD = o2::soa::Join<o2::aod::Hf3PBases, o2::aod::Hf3PPars, o2::aod::Hf3PParEs, o2::aod::Hf3PSels, o2::aod::Hf3PMls, o2::aod::Hf3PMcs, o2::aod::JLcIds>;
using JetTracksSubLc = o2::aod::JTrackLcSubs;

using CandidatesV0Data = o2::soa::Join<o2::aod::V0Cores, o2::aod::V0Extras, o2::aod::JV0Ids>;
using CandidatesV0MCD = o2::soa::Join<o2::aod::V0Cores, o2::aod::V0Extras, o2::aod::V0MCCores, o2::aod::JV0Ids, o2::aod::JV0Ids>; // add a table for McV0Labels with daughter prongs as well
using V0Daughters = o2::aod::DauTrackExtras;                                                                                      // linked by  V0Extras - check what this is
using CandidatesV0MCP = o2::soa::Join<o2::aod::JV0McParticles, o2::aod::JV0PIds>;

#endif // PWGJE_DATAMODEL_JET_H_
