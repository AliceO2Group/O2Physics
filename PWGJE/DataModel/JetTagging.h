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

/// \file JetTagging.h
/// \brief Table definitions for hf jet tagging
///
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>
/// \author Hanseo Park <hanseo.park@cern.ch>

#ifndef PWGJE_DATAMODEL_JETTAGGING_H_
#define PWGJE_DATAMODEL_JETTAGGING_H_

#include "PWGJE/Core/JetTaggingUtilities.h"
#include "PWGJE/DataModel/Jet.h" // IWYU pragma: keep

#include "Common/Core/RecoDecay.h"

#include <Framework/ASoA.h>

#include <Rtypes.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <utility>
#include <vector>

namespace o2::aod
{

// Defines derived extension data inside jets for tagging
namespace jtracktag
{
DECLARE_SOA_COLUMN(DcaXYZ, dcaXYZ, float);
DECLARE_SOA_COLUMN(SigmaDcaXYZ2, sigmaDcaXYZ2, float);
} // namespace jtracktag

DECLARE_SOA_TABLE(JTracksTag, "AOD", "JTracksTag",
                  jtracktag::DcaXYZ,
                  jtracktag::SigmaDcaXYZ2);

using JTrackTag = JTracksTag::iterator;

namespace secondary_vertex_params
{
DECLARE_SOA_COLUMN(XPrimaryVertex, xPrimaryVertex, float);
DECLARE_SOA_COLUMN(YPrimaryVertex, yPrimaryVertex, float);
DECLARE_SOA_COLUMN(ZPrimaryVertex, zPrimaryVertex, float);
DECLARE_SOA_COLUMN(XSecondaryVertex, xSecondaryVertex, float);
DECLARE_SOA_COLUMN(YSecondaryVertex, ySecondaryVertex, float);
DECLARE_SOA_COLUMN(ZSecondaryVertex, zSecondaryVertex, float);
DECLARE_SOA_COLUMN(Px, px, float);
DECLARE_SOA_COLUMN(Py, py, float);
DECLARE_SOA_COLUMN(Pz, pz, float);
DECLARE_SOA_COLUMN(E, e, float);
DECLARE_SOA_COLUMN(M, m, float);
DECLARE_SOA_COLUMN(Chi2PCA, chi2PCA, float);
DECLARE_SOA_COLUMN(Dispersion, dispersion, float);
DECLARE_SOA_COLUMN(ErrorDecayLength, errorDecayLength, float);
DECLARE_SOA_COLUMN(ErrorDecayLengthXY, errorDecayLengthXY, float);
DECLARE_SOA_DYNAMIC_COLUMN(RSecondaryVertex, rSecondaryVertex, [](float xVtxS, float yVtxS) -> float { return RecoDecay::sqrtSumOfSquares(xVtxS, yVtxS); });
DECLARE_SOA_DYNAMIC_COLUMN(Pt, pt, [](float px, float py) -> float { return RecoDecay::pt(px, py); });
DECLARE_SOA_DYNAMIC_COLUMN(P, p, [](float px, float py, float pz) -> float { return RecoDecay::p(px, py, pz); });
DECLARE_SOA_DYNAMIC_COLUMN(PVector, pVector, [](float px, float py, float pz) -> std::array<float, 3> { return std::array{px, py, pz}; });
DECLARE_SOA_DYNAMIC_COLUMN(Eta, eta, [](float px, float py, float pz) -> float { return RecoDecay::eta(std::array{px, py, pz}); });
DECLARE_SOA_DYNAMIC_COLUMN(Phi, phi, [](float px, float py) -> float { return RecoDecay::phi(px, py); });
DECLARE_SOA_DYNAMIC_COLUMN(Y, y, [](float px, float py, float pz, double m) -> float { return RecoDecay::y(std::array{px, py, pz}, m); });
DECLARE_SOA_DYNAMIC_COLUMN(DecayLength, decayLength, [](float xVtxP, float yVtxP, float zVtxP, float xVtxS, float yVtxS, float zVtxS) -> float { return RecoDecay::distance(std::array{xVtxP, yVtxP, zVtxP}, std::array{xVtxS, yVtxS, zVtxS}); });
DECLARE_SOA_DYNAMIC_COLUMN(DecayLengthXY, decayLengthXY, [](float xVtxP, float yVtxP, float xVtxS, float yVtxS) -> float { return RecoDecay::distanceXY(std::array{xVtxP, yVtxP}, std::array{xVtxS, yVtxS}); });
DECLARE_SOA_DYNAMIC_COLUMN(DecayLengthNormalised, decayLengthNormalised, [](float xVtxP, float yVtxP, float zVtxP, float xVtxS, float yVtxS, float zVtxS, float err) -> float { return RecoDecay::distance(std::array{xVtxP, yVtxP, zVtxP}, std::array{xVtxS, yVtxS, zVtxS}) / err; });
DECLARE_SOA_DYNAMIC_COLUMN(DecayLengthXYNormalised, decayLengthXYNormalised, [](float xVtxP, float yVtxP, float xVtxS, float yVtxS, float err) -> float { return RecoDecay::distanceXY(std::array{xVtxP, yVtxP}, std::array{xVtxS, yVtxS}) / err; });
DECLARE_SOA_DYNAMIC_COLUMN(Cpa, cpa, [](float xVtxP, float yVtxP, float zVtxP, float xVtxS, float yVtxS, float zVtxS, float px, float py, float pz) -> float { return RecoDecay::cpa(std::array{xVtxP, yVtxP, zVtxP}, std::array{xVtxS, yVtxS, zVtxS}, std::array{px, py, pz}); });
DECLARE_SOA_DYNAMIC_COLUMN(ImpactParameterXY, impactParameterXY, [](float xVtxP, float yVtxP, float zVtxP, float xVtxS, float yVtxS, float zVtxS, float px, float py, float pz) -> float { return RecoDecay::impParXY(std::array{xVtxP, yVtxP, zVtxP}, std::array{xVtxS, yVtxS, zVtxS}, std::array{px, py, pz}); });
} // namespace secondary_vertex_params

#define DECLARE_SV_TABLE(_jet_type_, _name_, _description_, _datatype_)                                                                                                                                                                                                                                                                                                                                             \
  namespace _datatype_##_name_##parameters                                                                                                                                                                                                                                                                                                                                                                          \
  {                                                                                                                                                                                                                                                                                                                                                                                                                 \
    DECLARE_SOA_INDEX_COLUMN(_jet_type_, jetIndex);                                                                                                                                                                                                                                                                                                                                                                 \
  }                                                                                                                                                                                                                                                                                                                                                                                                                 \
  DECLARE_SOA_TABLE(_datatype_##_name_##s, "AOD", _description_,                                                                                                                                                                                                                                                                                                                                                    \
                    o2::soa::Index<>,                                                                                                                                                                                                                                                                                                                                                                               \
                    _datatype_##_name_##parameters::_jet_type_##Id,                                                                                                                                                                                                                                                                                                                                                 \
                    secondary_vertex_params::XPrimaryVertex,                                                                                                                                                                                                                                                                                                                                                        \
                    secondary_vertex_params::YPrimaryVertex,                                                                                                                                                                                                                                                                                                                                                        \
                    secondary_vertex_params::ZPrimaryVertex,                                                                                                                                                                                                                                                                                                                                                        \
                    secondary_vertex_params::XSecondaryVertex,                                                                                                                                                                                                                                                                                                                                                      \
                    secondary_vertex_params::YSecondaryVertex,                                                                                                                                                                                                                                                                                                                                                      \
                    secondary_vertex_params::ZSecondaryVertex,                                                                                                                                                                                                                                                                                                                                                      \
                    secondary_vertex_params::Px,                                                                                                                                                                                                                                                                                                                                                                    \
                    secondary_vertex_params::Py,                                                                                                                                                                                                                                                                                                                                                                    \
                    secondary_vertex_params::Pz,                                                                                                                                                                                                                                                                                                                                                                    \
                    secondary_vertex_params::E,                                                                                                                                                                                                                                                                                                                                                                     \
                    secondary_vertex_params::M,                                                                                                                                                                                                                                                                                                                                                                     \
                    secondary_vertex_params::Chi2PCA,                                                                                                                                                                                                                                                                                                                                                               \
                    secondary_vertex_params::Dispersion,                                                                                                                                                                                                                                                                                                                                                            \
                    secondary_vertex_params::ErrorDecayLength,                                                                                                                                                                                                                                                                                                                                                      \
                    secondary_vertex_params::ErrorDecayLengthXY,                                                                                                                                                                                                                                                                                                                                                    \
                    secondary_vertex_params::RSecondaryVertex<secondary_vertex_params::XSecondaryVertex, secondary_vertex_params::YSecondaryVertex>,                                                                                                                                                                                                                                                                \
                    secondary_vertex_params::Pt<secondary_vertex_params::Px, secondary_vertex_params::Py>,                                                                                                                                                                                                                                                                                                          \
                    secondary_vertex_params::P<secondary_vertex_params::Px, secondary_vertex_params::Py, secondary_vertex_params::Pz>,                                                                                                                                                                                                                                                                              \
                    secondary_vertex_params::PVector<secondary_vertex_params::Px, secondary_vertex_params::Py, secondary_vertex_params::Pz>,                                                                                                                                                                                                                                                                        \
                    secondary_vertex_params::Eta<secondary_vertex_params::Px, secondary_vertex_params::Py, secondary_vertex_params::Pz>,                                                                                                                                                                                                                                                                            \
                    secondary_vertex_params::Phi<secondary_vertex_params::Px, secondary_vertex_params::Py>,                                                                                                                                                                                                                                                                                                         \
                    secondary_vertex_params::Y<secondary_vertex_params::Px, secondary_vertex_params::Py, secondary_vertex_params::Pz, secondary_vertex_params::M>,                                                                                                                                                                                                                                                  \
                    secondary_vertex_params::DecayLength<secondary_vertex_params::XPrimaryVertex, secondary_vertex_params::YPrimaryVertex, secondary_vertex_params::ZPrimaryVertex, secondary_vertex_params::XSecondaryVertex, secondary_vertex_params::YSecondaryVertex, secondary_vertex_params::ZSecondaryVertex>,                                                                                               \
                    secondary_vertex_params::DecayLengthXY<secondary_vertex_params::XPrimaryVertex, secondary_vertex_params::YPrimaryVertex, secondary_vertex_params::XSecondaryVertex, secondary_vertex_params::YSecondaryVertex>,                                                                                                                                                                                 \
                    secondary_vertex_params::DecayLengthNormalised<secondary_vertex_params::XPrimaryVertex, secondary_vertex_params::YPrimaryVertex, secondary_vertex_params::ZPrimaryVertex, secondary_vertex_params::XSecondaryVertex, secondary_vertex_params::YSecondaryVertex, secondary_vertex_params::ZSecondaryVertex, secondary_vertex_params::ErrorDecayLength>,                                          \
                    secondary_vertex_params::DecayLengthXYNormalised<secondary_vertex_params::XPrimaryVertex, secondary_vertex_params::YPrimaryVertex, secondary_vertex_params::XSecondaryVertex, secondary_vertex_params::YSecondaryVertex, secondary_vertex_params::ErrorDecayLengthXY>,                                                                                                                          \
                    secondary_vertex_params::Cpa<secondary_vertex_params::XPrimaryVertex, secondary_vertex_params::YPrimaryVertex, secondary_vertex_params::ZPrimaryVertex, secondary_vertex_params::XSecondaryVertex, secondary_vertex_params::YSecondaryVertex, secondary_vertex_params::ZSecondaryVertex, secondary_vertex_params::Px, secondary_vertex_params::Py, secondary_vertex_params::Pz>,                \
                    secondary_vertex_params::ImpactParameterXY<secondary_vertex_params::XPrimaryVertex, secondary_vertex_params::YPrimaryVertex, secondary_vertex_params::ZPrimaryVertex, secondary_vertex_params::XSecondaryVertex, secondary_vertex_params::YSecondaryVertex, secondary_vertex_params::ZSecondaryVertex, secondary_vertex_params::Px, secondary_vertex_params::Py, secondary_vertex_params::Pz>); \
  namespace _name_##indices                                                                                                                                                                                                                                                                                                                                                                                         \
  {                                                                                                                                                                                                                                                                                                                                                                                                                 \
    DECLARE_SOA_ARRAY_INDEX_COLUMN(_datatype_##_name_, secondaryVertices);                                                                                                                                                                                                                                                                                                                                          \
  }                                                                                                                                                                                                                                                                                                                                                                                                                 \
  DECLARE_SOA_TABLE(_datatype_##_name_##Indices, "AOD", _description_ "SVs", _name_##indices::_datatype_##_name_##Ids);

#define JETSV_TABLES_DEF(_jet_type_, _name_, _description_)      \
  DECLARE_SV_TABLE(_jet_type_##Jet, _name_, _description_, Data) \
  DECLARE_SV_TABLE(_jet_type_##MCDetectorLevelJet, _name_, _description_ "MCD", MCD)

JETSV_TABLES_DEF(Charged, SecondaryVertex3Prong, "3PRONG");
JETSV_TABLES_DEF(Charged, SecondaryVertex2Prong, "2PRONG");
JETSV_TABLES_DEF(Charged, SecondaryVertexNProng, "NPRONG");

// Defines the jet tagging table definition
namespace jettagging
{
namespace flavourdef
{
DECLARE_SOA_COLUMN(Origin, origin, int8_t);
} // namespace flavourdef
DECLARE_SOA_COLUMN(BitTaggedjet, bitTaggedjet, uint16_t);
DECLARE_SOA_COLUMN(JetProb, jetProb, float);
DECLARE_SOA_COLUMN(ScoreML, scoreML, float);
DECLARE_SOA_DYNAMIC_COLUMN(IsTagged, isTagged, [](uint16_t bit, BJetTaggingMethod method) -> bool { return TESTBIT(bit, method); });
} // namespace jettagging

#define JETFLAVOURDEF_TABLE_DEF(_jet_type_, _name_, _description_) \
  DECLARE_SOA_TABLE(_jet_type_##FlavourDef, "AOD", _description_ "FlavourDef", jettagging::flavourdef::Origin);

#define JETTAGGING_TABLE_DEF(_jet_type_, _name_, _description_)                       \
  namespace _name_##util                                                              \
  {                                                                                   \
    DECLARE_SOA_DYNAMIC_COLUMN(Dummy##_jet_type_##tagging, dummy##_jet_type##tagging, \
                               []() -> int { return 0; });                            \
  }                                                                                   \
  DECLARE_SOA_TABLE(_jet_type_##Tags, "AOD", _description_ "Tags",                    \
                    jettagging::BitTaggedjet, jettagging::JetProb,                    \
                    jettagging::ScoreML, jettagging::IsTagged<jettagging::BitTaggedjet>, _name_##util::Dummy##_jet_type_##tagging<>);

#define JETTAGGING_TABLES_DEF(_jet_type_, _description_)                                                       \
  JETTAGGING_TABLE_DEF(_jet_type_##Jet, _jet_type_##jet, _description_)                                        \
  JETFLAVOURDEF_TABLE_DEF(_jet_type_##MCDetectorLevelJet, _jet_type_##mcdetectorleveljet, _description_ "MCD") \
  JETTAGGING_TABLE_DEF(_jet_type_##MCDetectorLevelJet, _jet_type_##mcdetectorleveljet, _description_ "MCD")    \
  JETFLAVOURDEF_TABLE_DEF(_jet_type_##MCParticleLevelJet, _jet_type_##mcparticleleveljet, _description_ "MCP") \
  JETTAGGING_TABLE_DEF(_jet_type_##MCParticleLevelJet, _jet_type_##mcparticleleveljet, _description_ "MCP")

JETTAGGING_TABLES_DEF(Charged, "C");
JETTAGGING_TABLES_DEF(Full, "F");

// Tables for full inference of GNN-based b-jet tagger (including track origin, vertex predictions)
namespace jetgnnindex
{
DECLARE_SOA_INDEX_COLUMN(JTrack, track);             //! index in JTracks of this GNN input track (not the parent AOD Tracks index)
DECLARE_SOA_COLUMN(TrackIndex, trackIndex, uint8_t); //! node index (0..255) in the GNN input after track selection/truncation; need not match the constituent-array position
// DECLARE_SOA_COLUMN(EdgeIndex, edgeIndex, uint8_t);    //! unique in-jet edge index for the (i,j) constituent-track pair, i<j, triangular order i*(2N-i-1)/2+(j-i-1)
// DECLARE_SOA_DYNAMIC_COLUMN(TrackIds, trackIds, //! inverse of the triangular edge indexing above: for n constituent tracks, recovers the (i,j) in-jet track index pair (i<j) this edge connects
//                           [](uint8_t edgeIndex, int n) -> std::pair<uint8_t, uint8_t> {
//                             int i = 0;
//                             int offset = 0;
//                             while (offset + (n - 1 - i) <= static_cast<int>(edgeIndex)) {
//                               offset += (n - 1 - i);
//                               i++;
//                             }
//                             int j = i + 1 + (static_cast<int>(edgeIndex) - offset);
//                             return std::pair<uint8_t, uint8_t>{static_cast<uint8_t>(i), static_cast<uint8_t>(j)};
//                           });
} // namespace jetgnnindex

namespace jetgnnlabel
{
DECLARE_SOA_COLUMN(Origin, origin, uint8_t);          //! trkOrigin from vertexClustering() (0: non-physical primary, 1: charm, 2: beauty, 3: primary vertex, 4: other secondary)
DECLARE_SOA_COLUMN(VertexIndex, vertexIndex, int8_t); //! trkVtxIndex from vertexClustering() (-1: unassigned, 0: primary vertex, >0: secondary vertex, ranked by distance to PV)
// DECLARE_SOA_COLUMN(SameVertex, sameVertex, bool);      //! true if the (i,j) constituent-track pair of this edge share the same in-jet vertex cluster (same trkVtxIndex from vertexClustering())
} // namespace jetgnnlabel

namespace jetgnnscore
{
DECLARE_SOA_COLUMN(ScoreML, scoreML, float); //! GNN discriminant: one track-origin D_b per GnnTracks row
} // namespace jetgnnscore

// GnnTracks has one row per evaluated (jet, track) pair. Its jet and track indices
// identify the pair independently of the ordering in JetConstituents.
// GnnTrackIds is a positional jet extension: write one row for every jet in
// Jets order, including an empty array for jets without evaluated tracks. Store
// GnnTracks row indices (e.g. the backbone producer's lastIndex()), not JTrackIds.
#define DECLARE_JET_GNNTRACKINDEX_TABLE(_jet_type_, _description_)                                            \
  namespace _jet_type_##gnntrackindextable                                                                    \
  {                                                                                                           \
    DECLARE_SOA_INDEX_COLUMN(_jet_type_, jet);                                                                \
  }                                                                                                           \
  DECLARE_SOA_TABLE(_jet_type_##GnnTracks, "AOD", _description_,                                              \
                    o2::soa::Index<>,                                                                         \
                    _jet_type_##gnntrackindextable::_jet_type_##Id,                                           \
                    jetgnnindex::JTrackId,                                                                    \
                    jetgnnindex::TrackIndex);                                                                 \
  namespace _jet_type_##gnntrackindices                                                                       \
  {                                                                                                           \
    DECLARE_SOA_ARRAY_INDEX_COLUMN_FULL(_jet_type_##GnnTrack, tracksGnn, int32_t, _jet_type_##GnnTracks, ""); \
  }                                                                                                           \
  DECLARE_SOA_TABLE(_jet_type_##GnnTrackIds, "AOD", _description_ "I",                                        \
                    _jet_type_##gnntrackindices::_jet_type_##GnnTrackIds);

/*
#define DECLARE_JET_GNNEDGEINDEX_TABLE(_jet_type_, _description_)          \
  namespace _jet_type_##gnnedgeindextable                                 \
  {                                                                       \
    DECLARE_SOA_INDEX_COLUMN(_jet_type_, jet);                           \
  }                                                                      \
  DECLARE_SOA_TABLE(_jet_type_##GnnEdges, "AOD", _description_,          \
                    o2::soa::Index<>,                                   \
                    _jet_type_##gnnedgeindextable::_jet_type_##Id,      \
                    jetgnnindex::EdgeIndex,                             \
                    jetgnnindex::TrackIds<jetgnnindex::EdgeIndex>);     \
  namespace _jet_type_##gnnedgeindices                                  \
  {                                                                     \
    DECLARE_SOA_ARRAY_INDEX_COLUMN(_jet_type_##GnnEdge, edges);         \
  }                                                                     \
  DECLARE_SOA_TABLE(_jet_type_##GnnEdgeIds, "AOD", _description_ "I", \
                    _jet_type_##gnnedgeindices::_jet_type_##GnnEdgeIds);
*/

// In-jet track/edge index backbone tables (own the in-jet ordering + the back-reference to the original track)
DECLARE_JET_GNNTRACKINDEX_TABLE(ChargedJet, "CJGNT");
DECLARE_JET_GNNTRACKINDEX_TABLE(ChargedMCDetectorLevelJet, "CDJGNT");
// DECLARE_JET_GNNEDGEINDEX_TABLE(ChargedJet, "CJGNE");
// DECLARE_JET_GNNEDGEINDEX_TABLE(ChargedMCDetectorLevelJet, "CDJGNE");

#undef DECLARE_JET_GNNTRACKINDEX_TABLE
// #undef DECLARE_JET_GNNEDGEINDEX_TABLE

// Truth labels from vertexClustering() (MC detector level only). Thin tables: join positionally with
// ChargedMCDetectorLevelJetGnnTracks / ChargedMCDetectorLevelJetGnnEdges (same producer loop, same row order).
// TrackOriginDef must have exactly one row per GnnTracks backbone row, aligned
// with GnnTrackOriginTags after GNN input selection/truncation. It is not a JTracks
// extension: a JTrack used by multiple jets has a truth row for each backbone row.
// tracksGnnIds() addresses this shared row position through a join containing
// the backbone, the origin definition and the prediction scores.
// VertexIndexDef uses the same GnnTracks row order and keeps the full-jet vertex
// labels, without renumbering after GNN track selection/truncation.
DECLARE_SOA_TABLE(ChargedMCDetectorLevelJetTrackOriginDef, "AOD", "CDJTOD", jetgnnlabel::Origin);
DECLARE_SOA_TABLE(ChargedMCDetectorLevelJetVertexIndexDef, "AOD", "CDJTVD", jetgnnlabel::VertexIndex);
// DECLARE_SOA_TABLE(ChargedMCDetectorLevelJetEdgeDef, "AOD", "CDJEDF", jetgnnlabel::SameVertex);

// GNN track-origin outputs, data and MC detector level. Thin tables: write one
// scalar score per GnnTracks row, in the same row order as the backbone table.
DECLARE_SOA_TABLE(ChargedJetGnnTrackOriginTags, "AOD", "CJTRKT", jetgnnscore::ScoreML);
DECLARE_SOA_TABLE(ChargedMCDetectorLevelJetGnnTrackOriginTags, "AOD", "CDJTRKT", jetgnnscore::ScoreML);

// GNN ScoreML outputs, edge level, data and MC detector level. Thin tables: join positionally with
// ChargedJetGnnEdges / ChargedMCDetectorLevelJetGnnEdges.
// DECLARE_SOA_TABLE(ChargedJetGnnEdgeTags, "AOD", "CJEDGT", jetgnnscore::ScoreML);
// DECLARE_SOA_TABLE(ChargedMCDetectorLevelJetGnnEdgeTags, "AOD", "CDJEDGT", jetgnnscore::ScoreML);

} // namespace o2::aod

#endif // PWGJE_DATAMODEL_JETTAGGING_H_
