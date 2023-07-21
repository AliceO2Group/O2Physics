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
/// \author Jochen Klein
/// \author Nima Zardoshti
/// \author Raymond Ehlers

#ifndef PWGJE_DATAMODEL_JET_H_
#define PWGJE_DATAMODEL_JET_H_

#include <cmath>
#include "Framework/AnalysisDataModel.h"
#include "PWGJE/DataModel/EMCALClusters.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"

namespace o2::aod
{
namespace jet
{
DECLARE_SOA_INDEX_COLUMN(Collision, collision);
DECLARE_SOA_INDEX_COLUMN(McCollision, mcCollision);
DECLARE_SOA_COLUMN(Pt, pt, float);         //!
DECLARE_SOA_COLUMN(Eta, eta, float);       //!
DECLARE_SOA_COLUMN(Phi, phi, float);       //!
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

// Constituent sub
namespace constituentssub
{
// Jet index column will be added in the macro
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(Energy, energy, float);
DECLARE_SOA_COLUMN(Mass, mass, float);
DECLARE_SOA_COLUMN(Source, source, int);
DECLARE_SOA_DYNAMIC_COLUMN(Px, px,
                           [](float pt, float phi) -> float { return pt * std::cos(phi); });
DECLARE_SOA_DYNAMIC_COLUMN(Py, py,
                           [](float pt, float phi) -> float { return pt * std::sin(phi); });
DECLARE_SOA_DYNAMIC_COLUMN(Pz, pz,
                           [](float pt, float eta) -> float { return pt * std::sinh(eta); });
DECLARE_SOA_DYNAMIC_COLUMN(P, p,
                           [](float pt, float eta) -> float { return pt * std::cosh(eta); });
} // namespace constituentssub
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
  // DECLARE_SOA_INDEX_COLUMN_FULL(_jet_type_, matchedJetGeo, int32_t, _jet_type_##s, "_geo");
  // DECLARE_SOA_INDEX_COLUMN_FULL(_jet_type_, matchedJetPt, int32_t, _jet_type_##s, "_pt");
  // DECLARE_SOA_INDEX_COLUMN_FULL(_jet_type_, matchedJetCand, int32_t, _jet_type_##s, "_hf");

#define DECLARE_CONSTITUENTS_TABLE(_jet_type_, _name_, _Description_, _track_type_, _cand_type_)      \
  namespace _name_##constituents                                                                      \
  {                                                                                                   \
    DECLARE_SOA_INDEX_COLUMN(_jet_type_, jet);                                                        \
    DECLARE_SOA_ARRAY_INDEX_COLUMN(_track_type_, tracks);                                             \
    DECLARE_SOA_ARRAY_INDEX_COLUMN(EMCALCluster, clusters);                                           \
    DECLARE_SOA_ARRAY_INDEX_COLUMN_FULL(HfCandidates, hfcandidates, int32_t, _cand_type_, "_hfcand"); \
  }                                                                                                   \
  DECLARE_SOA_TABLE(_jet_type_##Constituents, "AOD", _Description_ "CONSTS",                          \
                    _name_##constituents::_jet_type_##Id,                                             \
                    _name_##constituents::_track_type_##Ids,                                          \
                    _name_##constituents::EMCALClusterIds,                                            \
                    _name_##constituents::HfCandidatesIds);

// Defines the jet constituent sub table
// NOTE: This relies on the jet index column being defined in the constituents namespace.
//       Since these are always paired together, there's no point in redefining them.
#define DECLARE_CONSTITUENTS_SUB_TABLE(_jet_type_, _name_, _Description_)           \
  DECLARE_SOA_TABLE(_jet_type_##ConstituentsSub, "AOD", _Description_ "CONSTSUB",   \
                    _name_##constituents::_jet_type_##Id,                           \
                    constituentssub::Pt,                                            \
                    constituentssub::Eta,                                           \
                    constituentssub::Phi,                                           \
                    constituentssub::Energy,                                        \
                    constituentssub::Mass,                                          \
                    constituentssub::Source,                                        \
                    constituentssub::Px<constituentssub::Pt, constituentssub::Phi>, \
                    constituentssub::Py<constituentssub::Pt, constituentssub::Phi>, \
                    constituentssub::Pz<constituentssub::Pt, constituentssub::Eta>, \
                    constituentssub::P<constituentssub::Pt, constituentssub::Eta>);

// combine definition of tables for jets, constituents, and substructure
#define DECLARE_JET_TABLES(_collision_name_, _jet_type_, _const_type_, _hfcand_type_, _description_)        \
  DECLARE_JET_TABLE(_collision_name_, _jet_type_##Jet, _jet_type_##jet, _description_);                     \
  using _jet_type_##Jet = _jet_type_##Jet##s::iterator;                                                     \
  DECLARE_CONSTITUENTS_TABLE(_jet_type_##Jet, _jet_type_##jet, _description_, _const_type_, _hfcand_type_); \
  using _jet_type_##Jet##Constituent = _jet_type_##Jet##Constituents::iterator;                             \
  DECLARE_CONSTITUENTS_SUB_TABLE(_jet_type_##Jet, _jet_type_##jet, _description_);                          \
  using _jet_type_##Jet##ConstituentSub = _jet_type_##Jet##ConstituentsSub::iterator;

#define DECLARE_JETMATCHING_TABLE(_jet_type_base_, _jet_type_tag_, _description_)               \
  DECLARE_SOA_TABLE(_jet_type_base_##JetsMatchedTo##_jet_type_tag_##Jets, "AOD", _description_, \
                    _jet_type_tag_##jetmatchingGeo::_jet_type_tag_##JetIds,                     \
                    _jet_type_tag_##jetmatchingPt::_jet_type_tag_##JetIds,                      \
                    _jet_type_tag_##jetmatchingCand::_jet_type_tag_##JetIds);                   \
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
#define DECLARE_JET_TABLES_LEVELS(_jet_type_, _hfcand_type_, _shortname_)                                      \
  DECLARE_JET_TABLES(Collision, _jet_type_, Track, _hfcand_type_, _shortname_ "JET")                           \
  DECLARE_JET_TABLES(Collision, _jet_type_##MCDetectorLevel, Track, _hfcand_type_, _shortname_ "DJET")         \
  DECLARE_JET_TABLES(McCollision, _jet_type_##MCParticleLevel, McParticle, McParticles, _shortname_ "PJET")    \
  DECLARE_JETMATCHING_TABLE(_jet_type_##MCDetectorLevel, _jet_type_##MCParticleLevel, _shortname_ "JETMD2P")   \
  DECLARE_JETMATCHING_TABLE(_jet_type_##MCParticleLevel, _jet_type_##MCDetectorLevel, _shortname_ "JETMP2D")   \
  DECLARE_MCEVENTWEIGHT_TABLE(_jet_type_##MCDetectorLevel, _jet_type_##MCDetectorLevel, _shortname_ "JETMDEW") \
  DECLARE_MCEVENTWEIGHT_TABLE(_jet_type_##MCParticleLevel, _jet_type_##MCParticleLevel, _shortname_ "JETMPEW")

namespace o2::aod
{
DECLARE_JET_TABLES_LEVELS(Charged, HfCand2Prong, "C");
DECLARE_JET_TABLES_LEVELS(Full, HfCand2Prong, "F");
DECLARE_JET_TABLES_LEVELS(Neutral, HfCand2Prong, "N");
DECLARE_JET_TABLES_LEVELS(D0Charged, HfCand2Prong, "D0");
DECLARE_JET_TABLES_LEVELS(LcCharged, HfCand3Prong, "Lc");
DECLARE_JET_TABLES_LEVELS(BplusCharged, HfCandBplus, "BPl");

// Hybrid intermediate
DECLARE_JET_TABLES(Collision, HybridIntermediate, Track, HfCand2Prong, "JEHYIN");
} // namespace o2::aod

#endif // PWGJE_DATAMODEL_JET_H_
