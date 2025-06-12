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

#include <cmath>
#include <vector>
#include "Framework/AnalysisDataModel.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/Core/JetTaggingUtilities.h"

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
DECLARE_SOA_COLUMN(XPrimaryVertex, xPVertex, float); // o2-linter: disable=name/o2-column
DECLARE_SOA_COLUMN(YPrimaryVertex, yPVertex, float); // o2-linter: disable=name/o2-column
DECLARE_SOA_COLUMN(ZPrimaryVertex, zPVertex, float); // o2-linter: disable=name/o2-column
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
DECLARE_SOA_DYNAMIC_COLUMN(CPA, cpa, [](float xVtxP, float yVtxP, float zVtxP, float xVtxS, float yVtxS, float zVtxS, float px, float py, float pz) -> float { return RecoDecay::cpa(std::array{xVtxP, yVtxP, zVtxP}, std::array{xVtxS, yVtxS, zVtxS}, std::array{px, py, pz}); }); // o2-linter: disable=name/o2-column
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
                    secondary_vertex_params::CPA<secondary_vertex_params::XPrimaryVertex, secondary_vertex_params::YPrimaryVertex, secondary_vertex_params::ZPrimaryVertex, secondary_vertex_params::XSecondaryVertex, secondary_vertex_params::YSecondaryVertex, secondary_vertex_params::ZSecondaryVertex, secondary_vertex_params::Px, secondary_vertex_params::Py, secondary_vertex_params::Pz>,                \
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

} // namespace o2::aod

#endif // PWGJE_DATAMODEL_JETTAGGING_H_
