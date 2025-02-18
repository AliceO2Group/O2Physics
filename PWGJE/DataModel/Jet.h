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
#include "PWGJE/DataModel/JetReducedDataV0.h"
#include "PWGJE/DataModel/JetReducedDataDQ.h"
#include "PWGJE/DataModel/JetSubtraction.h"

#include "PWGHF/DataModel/DerivedTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGDQ/DataModel/ReducedInfoTables.h"

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

#define DECLARE_CONSTITUENTS_TABLE(_jet_type_, _name_, _Description_, _track_type_, _cand_type_) \
  namespace _name_##constituents                                                                 \
  {                                                                                              \
    DECLARE_SOA_INDEX_COLUMN(_jet_type_, jet);                                                   \
    DECLARE_SOA_ARRAY_INDEX_COLUMN(_track_type_, tracks);                                        \
    DECLARE_SOA_ARRAY_INDEX_COLUMN(JCluster, clusters);                                          \
    DECLARE_SOA_ARRAY_INDEX_COLUMN_FULL(Candidates, candidates, int32_t, _cand_type_, "_cand");  \
  }                                                                                              \
  DECLARE_SOA_TABLE(_jet_type_##Constituents, "AOD", _Description_ "C",                          \
                    _name_##constituents::_jet_type_##Id,                                        \
                    _name_##constituents::_track_type_##Ids,                                     \
                    _name_##constituents::JClusterIds,                                           \
                    _name_##constituents::CandidatesIds);

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

#define STRINGIFY(x) #x

// add duplicate tables for each predefined jet table so that the same jets can be run with multiple settings
#define DECLARE_JET_DUPLICATE_TABLES_LEVELS(_jet_type_, _subtracted_track_type_, _hfcand_type_, _hfparticle_type_, _shortname_, _duplicatenumber_)                   \
  DECLARE_JET_TABLES_LEVELS(_jet_type_##_duplicatenumber_, _subtracted_track_type_, _hfcand_type_, _hfparticle_type_, _shortname_ STRINGIFY(_duplicatenumber_))      \
  DECLARE_JETMATCHING_TABLE(_jet_type_, _jet_type_##_duplicatenumber_, _shortname_ "JET2" STRINGIFY(_duplicatenumber_))                                              \
  DECLARE_JETMATCHING_TABLE(_jet_type_##_duplicatenumber_, _jet_type_, _shortname_ "JET" STRINGIFY(_duplicatenumber_) "2")                                           \
  DECLARE_JETMATCHING_TABLE(_jet_type_##MCDetectorLevel, _jet_type_##_duplicatenumber_##MCDetectorLevel, _shortname_ "JETD2" STRINGIFY(_duplicatenumber_))           \
  DECLARE_JETMATCHING_TABLE(_jet_type_##_duplicatenumber_##MCDetectorLevel, _jet_type_##MCDetectorLevel, _shortname_ "JET" STRINGIFY(_duplicatenumber_) "2D")        \
  DECLARE_JETMATCHING_TABLE(_jet_type_##MCParticleLevel, _jet_type_##_duplicatenumber_##MCParticleLevel, _shortname_ "JETP2" STRINGIFY(_duplicatenumber_))           \
  DECLARE_JETMATCHING_TABLE(_jet_type_##_duplicatenumber_##MCParticleLevel, _jet_type_##MCParticleLevel, _shortname_ "JET" STRINGIFY(_duplicatenumber_) "2P")        \
  DECLARE_JETMATCHING_TABLE(_jet_type_##EventWiseSubtracted, _jet_type_##_duplicatenumber_##EventWiseSubtracted, _shortname_ "JETEWS2" STRINGIFY(_duplicatenumber_)) \
  DECLARE_JETMATCHING_TABLE(_jet_type__duplicatenumber_##EventWiseSubtracted, _jet_type_##EventWiseSubtracted, _shortname_ "JET" STRINGIFY(_duplicatenumber_) "2EWS")

DECLARE_JET_TABLES_LEVELS(Charged, JTrackSub, HfD0Bases, HfD0PBases, "C");
DECLARE_JET_TABLES_LEVELS(Full, JTrackSub, HfD0Bases, HfD0PBases, "F");
DECLARE_JET_TABLES_LEVELS(Neutral, JTrackSub, HfD0Bases, HfD0PBases, "N");
DECLARE_JET_TABLES_LEVELS(D0Charged, JTrackD0Sub, HfD0Bases, HfD0PBases, "D0");
DECLARE_JET_TABLES_LEVELS(DplusCharged, JTrackDplusSub, HfDplusBases, HfDplusPBases, "DP");
DECLARE_JET_TABLES_LEVELS(LcCharged, JTrackLcSub, HfLcBases, HfLcPBases, "Lc");
DECLARE_JET_TABLES_LEVELS(BplusCharged, JTrackBplusSub, HfBplusBases, HfBplusPBases, "BP");
DECLARE_JET_TABLES_LEVELS(V0Charged, JTrackSub, V0Cores, JV0Mcs, "V0");
DECLARE_JET_TABLES_LEVELS(DielectronCharged, JTrackSub, Dielectrons, JDielectronMcs, "DIEL");

// duplicate jet tables (added as needed for analyses)
DECLARE_JET_DUPLICATE_TABLES_LEVELS(Charged, JTrackSub, HfD0Bases, HfD0PBases, "C", 1);

#undef DECLARE_JET_TABLE
#undef DECLARE_CONSTITUENTS_TABLE
#undef DECLARE_JET_TABLES
#undef DECLARE_JETMATCHING_TABLE
#undef DECLARE_MCEVENTWEIGHT_TABLE
#undef DECLARE_JET_TABLES_LEVELS
#undef STRINGIFY
#undef DECLARE_JET_DUPLICATE_TABLES_LEVELS

using JetCollisions = o2::soa::Join<JCollisions, JCollisionMcInfos>;
using JetCollision = JetCollisions::iterator;
using JetCollisionsMCD = o2::soa::Join<JetCollisions, JMcCollisionLbs>;
using JetCollisionMCD = o2::soa::Join<JetCollisions, JMcCollisionLbs>::iterator;
using JetTracks = JTracks;
using JetTracksMCD = o2::soa::Join<JetTracks, JMcTrackLbs>;
using JetTracksSub = JTrackSubs;
using JetClusters = o2::soa::Join<JClusters, JClustersCorrectedEnergies, JClusterTracks>;
using JetClustersMCD = o2::soa::Join<JClusters, JClustersCorrectedEnergies, JClusterTracks, JMcClusterLbs>;

using JetMcCollisions = JMcCollisions;
using JetMcCollision = JetMcCollisions::iterator;
using JetParticles = JMcParticles;
using JetParticlesSub = JMcParticleSubs;

using CollisionsD0 = o2::soa::Join<HfD0CollBases, JD0CollisionIds>;
using CandidatesD0Data = o2::soa::Join<HfD0Bases, HfD0Pars, HfD0ParEs, HfD0Sels, HfD0Mls, JD0Ids>;
using CandidatesD0MCD = o2::soa::Join<HfD0Bases, HfD0Pars, HfD0ParEs, HfD0Sels, HfD0Mls, HfD0Mcs, JD0Ids>;
using JetTracksSubD0 = JTrackD0Subs;
using JetParticlesSubD0 = JMcParticleD0Subs;
using McCollisionsD0 = o2::soa::Join<HfD0McCollBases, JD0McCollisionIds>;
using CandidatesD0MCP = o2::soa::Join<HfD0PBases, JD0PIds>;

using CollisionsDplus = o2::soa::Join<HfDplusCollBases, JDplusCollisionIds>;
using CandidatesDplusData = o2::soa::Join<HfDplusBases, HfDplusPars, HfDplusParEs, HfDplusSels, HfDplusMls, JDplusIds>;
using CandidatesDplusMCD = o2::soa::Join<HfDplusBases, HfDplusPars, HfDplusParEs, HfDplusSels, HfDplusMls, HfDplusMcs, JDplusIds>;
using JetTracksSubDplus = JTrackDplusSubs;
using JetParticlesSubDplus = JMcParticleDplusSubs;
using McCollisionsDplus = o2::soa::Join<HfDplusMcCollBases, JDplusMcCollisionIds>;
using CandidatesDplusMCP = o2::soa::Join<HfDplusPBases, JDplusPIds>;

using CollisionsLc = o2::soa::Join<HfLcCollBases, JLcCollisionIds>;
using CandidatesLcData = o2::soa::Join<HfLcBases, HfLcPars, HfLcParEs, HfLcSels, HfLcMls, JLcIds>;
using CandidatesLcMCD = o2::soa::Join<HfLcBases, HfLcPars, HfLcParEs, HfLcSels, HfLcMls, HfLcMcs, JLcIds>;
using JetTracksSubLc = JTrackLcSubs;
using JetParticlesSubLc = JMcParticleLcSubs;
using McCollisionsLc = o2::soa::Join<HfLcMcCollBases, JLcMcCollisionIds>;
using CandidatesLcMCP = o2::soa::Join<HfLcPBases, JLcPIds>;

using CollisionsBplus = o2::soa::Join<HfBplusCollBases, JBplusCollisionIds>;
using CandidatesBplusData = o2::soa::Join<HfBplusBases, HfBplusPars, HfBplusParEs, HfBplusParD0s, HfBplusSels, HfBplusMls, HfBplusMlD0s, JBplusIds>;
using CandidatesBplusMCD = o2::soa::Join<HfBplusBases, HfBplusPars, HfBplusParEs, HfBplusParD0s, HfBplusSels, HfBplusMls, HfBplusMlD0s, HfBplusMcs, JBplusIds>;
using JetTracksSubBplus = JTrackBplusSubs;
using JetParticlesSubBplus = JMcParticleBplusSubs;
using McCollisionsBplus = o2::soa::Join<HfBplusMcCollBases, JBplusMcCollisionIds>;
using CandidatesBplusMCP = o2::soa::Join<HfBplusPBases, JBplusPIds>;

using CandidatesV0Data = o2::soa::Join<V0Cores, JV0Ids>;
using CandidatesV0MCD = o2::soa::Join<V0Cores, V0MCCores, JV0Ids>;
// using V0Daughters = DauTrackExtras;
using McCollisionsV0 = o2::soa::Join<JV0McCollisions, JV0McCollisionIds>;
using CandidatesV0MCP = o2::soa::Join<JV0Mcs, JV0McIds>;

using CollisionsDielectron = o2::soa::Join<ReducedEvents, JDielectronCollisionIds>;
using CandidatesDielectronData = o2::soa::Join<Dielectrons, JDielectronIds>;
using CandidatesDielectronMCD = o2::soa::Join<Dielectrons, JDielectronIds>;
using JetTracksSubDielectron = JTrackDielectronSubs;
using JetParticlesSubDielectron = JMcParticleDielectronSubs;
using McCollisionsDielectron = o2::soa::Join<JDielectronMcCollisions, JDielectronMcCollisionIds>;
using CandidatesDielectronMCP = o2::soa::Join<JDielectronMcs, JDielectronMcIds>;

} // namespace o2::aod

#endif // PWGJE_DATAMODEL_JET_H_
