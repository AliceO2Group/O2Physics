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
/// \brief Table definitions for hf jet tagging
///
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>

#ifndef PWGJE_DATAMODEL_JETTAGGING_H_
#define PWGJE_DATAMODEL_JETTAGGING_H_

#include <cmath>
#include <vector>
#include "Framework/AnalysisDataModel.h"
#include "PWGJE/DataModel/Jet.h"

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

namespace SecondaryVertexParams
{
DECLARE_SOA_COLUMN(XPrimaryVertex, xPVertex, float);
DECLARE_SOA_COLUMN(YPrimaryVertex, yPVertex, float);
DECLARE_SOA_COLUMN(ZPrimaryVertex, zPVertex, float);
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
DECLARE_SOA_DYNAMIC_COLUMN(CPA, cpa, [](float xVtxP, float yVtxP, float zVtxP, float xVtxS, float yVtxS, float zVtxS, float px, float py, float pz) -> float { return RecoDecay::cpa(std::array{xVtxP, yVtxP, zVtxP}, std::array{xVtxS, yVtxS, zVtxS}, std::array{px, py, pz}); });
DECLARE_SOA_DYNAMIC_COLUMN(ImpactParameterXY, impactParameterXY, [](float xVtxP, float yVtxP, float zVtxP, float xVtxS, float yVtxS, float zVtxS, float px, float py, float pz) -> float { return RecoDecay::impParXY(std::array{xVtxP, yVtxP, zVtxP}, std::array{xVtxS, yVtxS, zVtxS}, std::array{px, py, pz}); });
} // namespace SecondaryVertexParams

#define DECLARE_SV_TABLE(_jet_type_, _name_, _description_, _datatype_)                                                                                                                                                                                                                                                                                                                         \
  namespace _datatype_##_name_##Parameters                                                                                                                                                                                                                                                                                                                                                      \
  {                                                                                                                                                                                                                                                                                                                                                                                             \
    DECLARE_SOA_INDEX_COLUMN(_jet_type_, jetIndex);                                                                                                                                                                                                                                                                                                                                             \
  }                                                                                                                                                                                                                                                                                                                                                                                             \
  DECLARE_SOA_TABLE(_datatype_##_name_##s, "AOD", _description_,                                                                                                                                                                                                                                                                                                                                \
                    o2::soa::Index<>,                                                                                                                                                                                                                                                                                                                                                           \
                    _datatype_##_name_##Parameters::_jet_type_##Id,                                                                                                                                                                                                                                                                                                                             \
                    SecondaryVertexParams::XPrimaryVertex,                                                                                                                                                                                                                                                                                                                                      \
                    SecondaryVertexParams::YPrimaryVertex,                                                                                                                                                                                                                                                                                                                                      \
                    SecondaryVertexParams::ZPrimaryVertex,                                                                                                                                                                                                                                                                                                                                      \
                    SecondaryVertexParams::XSecondaryVertex,                                                                                                                                                                                                                                                                                                                                    \
                    SecondaryVertexParams::YSecondaryVertex,                                                                                                                                                                                                                                                                                                                                    \
                    SecondaryVertexParams::ZSecondaryVertex,                                                                                                                                                                                                                                                                                                                                    \
                    SecondaryVertexParams::Px,                                                                                                                                                                                                                                                                                                                                                  \
                    SecondaryVertexParams::Py,                                                                                                                                                                                                                                                                                                                                                  \
                    SecondaryVertexParams::Pz,                                                                                                                                                                                                                                                                                                                                                  \
                    SecondaryVertexParams::E,                                                                                                                                                                                                                                                                                                                                                   \
                    SecondaryVertexParams::M,                                                                                                                                                                                                                                                                                                                                                   \
                    SecondaryVertexParams::Chi2PCA,                                                                                                                                                                                                                                                                                                                                             \
                    SecondaryVertexParams::Dispersion,                                                                                                                                                                                                                                                                                                                                          \
                    SecondaryVertexParams::ErrorDecayLength,                                                                                                                                                                                                                                                                                                                                    \
                    SecondaryVertexParams::ErrorDecayLengthXY,                                                                                                                                                                                                                                                                                                                                  \
                    SecondaryVertexParams::RSecondaryVertex<SecondaryVertexParams::XSecondaryVertex, SecondaryVertexParams::YSecondaryVertex>,                                                                                                                                                                                                                                                  \
                    SecondaryVertexParams::Pt<SecondaryVertexParams::Px, SecondaryVertexParams::Py>,                                                                                                                                                                                                                                                                                            \
                    SecondaryVertexParams::P<SecondaryVertexParams::Px, SecondaryVertexParams::Py, SecondaryVertexParams::Pz>,                                                                                                                                                                                                                                                                  \
                    SecondaryVertexParams::PVector<SecondaryVertexParams::Px, SecondaryVertexParams::Py, SecondaryVertexParams::Pz>,                                                                                                                                                                                                                                                            \
                    SecondaryVertexParams::Eta<SecondaryVertexParams::Px, SecondaryVertexParams::Py, SecondaryVertexParams::Pz>,                                                                                                                                                                                                                                                                \
                    SecondaryVertexParams::Phi<SecondaryVertexParams::Px, SecondaryVertexParams::Py>,                                                                                                                                                                                                                                                                                           \
                    SecondaryVertexParams::Y<SecondaryVertexParams::Px, SecondaryVertexParams::Py, SecondaryVertexParams::Pz, SecondaryVertexParams::M>,                                                                                                                                                                                                                                        \
                    SecondaryVertexParams::DecayLength<SecondaryVertexParams::XPrimaryVertex, SecondaryVertexParams::YPrimaryVertex, SecondaryVertexParams::ZPrimaryVertex, SecondaryVertexParams::XSecondaryVertex, SecondaryVertexParams::YSecondaryVertex, SecondaryVertexParams::ZSecondaryVertex>,                                                                                         \
                    SecondaryVertexParams::DecayLengthXY<SecondaryVertexParams::XPrimaryVertex, SecondaryVertexParams::YPrimaryVertex, SecondaryVertexParams::XSecondaryVertex, SecondaryVertexParams::YSecondaryVertex>,                                                                                                                                                                       \
                    SecondaryVertexParams::DecayLengthNormalised<SecondaryVertexParams::XPrimaryVertex, SecondaryVertexParams::YPrimaryVertex, SecondaryVertexParams::ZPrimaryVertex, SecondaryVertexParams::XSecondaryVertex, SecondaryVertexParams::YSecondaryVertex, SecondaryVertexParams::ZSecondaryVertex, SecondaryVertexParams::ErrorDecayLength>,                                      \
                    SecondaryVertexParams::DecayLengthXYNormalised<SecondaryVertexParams::XPrimaryVertex, SecondaryVertexParams::YPrimaryVertex, SecondaryVertexParams::XSecondaryVertex, SecondaryVertexParams::YSecondaryVertex, SecondaryVertexParams::ErrorDecayLengthXY>,                                                                                                                  \
                    SecondaryVertexParams::CPA<SecondaryVertexParams::XPrimaryVertex, SecondaryVertexParams::YPrimaryVertex, SecondaryVertexParams::ZPrimaryVertex, SecondaryVertexParams::XSecondaryVertex, SecondaryVertexParams::YSecondaryVertex, SecondaryVertexParams::ZSecondaryVertex, SecondaryVertexParams::Px, SecondaryVertexParams::Py, SecondaryVertexParams::Pz>,                \
                    SecondaryVertexParams::ImpactParameterXY<SecondaryVertexParams::XPrimaryVertex, SecondaryVertexParams::YPrimaryVertex, SecondaryVertexParams::ZPrimaryVertex, SecondaryVertexParams::XSecondaryVertex, SecondaryVertexParams::YSecondaryVertex, SecondaryVertexParams::ZSecondaryVertex, SecondaryVertexParams::Px, SecondaryVertexParams::Py, SecondaryVertexParams::Pz>); \
  namespace _name_##Indices                                                                                                                                                                                                                                                                                                                                                                     \
  {                                                                                                                                                                                                                                                                                                                                                                                             \
    DECLARE_SOA_ARRAY_INDEX_COLUMN(_datatype_##_name_, secondaryVertices);                                                                                                                                                                                                                                                                                                                      \
  }                                                                                                                                                                                                                                                                                                                                                                                             \
  DECLARE_SOA_TABLE(_datatype_##_name_##Indices, "AOD", _description_ "SVs", _name_##Indices::_datatype_##_name_##Ids);

#define JETSV_TABLES_DEF(_jet_type_, _name_, _description_)      \
  DECLARE_SV_TABLE(_jet_type_##Jet, _name_, _description_, Data) \
  DECLARE_SV_TABLE(_jet_type_##MCDetectorLevelJet, _name_, _description_ "MCD", MCD)

JETSV_TABLES_DEF(Charged, SecondaryVertex3Prong, "3PRONG");
JETSV_TABLES_DEF(Charged, SecondaryVertex2Prong, "2PRONG");

// Defines the jet substrcuture table definition
#define JETTAGGING_TABLE_DEF(_jet_type_, _name_, _description_)       \
  namespace _name_##tagging                                           \
  {                                                                   \
    DECLARE_SOA_COLUMN(Origin, origin, int);                          \
    DECLARE_SOA_COLUMN(JetProb, jetProb, std::vector<float>);         \
    DECLARE_SOA_COLUMN(FlagtaggedjetIP, flagtaggedjetIP, bool);       \
    DECLARE_SOA_COLUMN(FlagtaggedjetIPxyz, flagtaggedjetIPxyz, bool); \
    DECLARE_SOA_COLUMN(FlagtaggedjetSV, flagtaggedjetSV, bool);       \
    DECLARE_SOA_COLUMN(FlagtaggedjetSVxyz, flagtaggedjetSVxyz, bool); \
  }                                                                   \
  DECLARE_SOA_TABLE(_jet_type_##Tags, "AOD", _description_ "Tags", _name_##tagging::Origin, _name_##tagging::JetProb, _name_##tagging::FlagtaggedjetIP, _name_##tagging::FlagtaggedjetIPxyz, _name_##tagging::FlagtaggedjetSV, _name_##tagging::FlagtaggedjetSVxyz);

#define JETTAGGING_TABLES_DEF(_jet_type_, _description_)                                                    \
  JETTAGGING_TABLE_DEF(_jet_type_##Jet, _jet_type_##jet, _description_)                                     \
  JETTAGGING_TABLE_DEF(_jet_type_##MCDetectorLevelJet, _jet_type_##mcdetectorleveljet, _description_ "MCD") \
  JETTAGGING_TABLE_DEF(_jet_type_##MCParticleLevelJet, _jet_type_##mcparticleleveljet, _description_ "MCP")

JETTAGGING_TABLES_DEF(Charged, "C");
JETTAGGING_TABLES_DEF(Full, "F");

} // namespace o2::aod

#endif // PWGJE_DATAMODEL_JETTAGGING_H_
