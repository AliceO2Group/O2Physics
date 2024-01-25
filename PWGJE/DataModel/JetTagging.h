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
/// \author Hanseo Park <hanseo.park@cern.ch>

#ifndef PWGJE_DATAMODEL_JETTAGGING_H_
#define PWGJE_DATAMODEL_JETTAGGING_H_

#include <cmath>
#include <vector>
#include "Framework/AnalysisDataModel.h"
#include "PWGJE/DataModel/Jet.h"

using namespace o2::analysis;

namespace o2::aod
{
namespace jtracktagdca
{
DECLARE_SOA_COLUMN(DcaX, dcaX, float);
DECLARE_SOA_COLUMN(DcaY, dcaY, float);
DECLARE_SOA_COLUMN(DcaXY, dcaXY, float);
DECLARE_SOA_COLUMN(DcaZ, dcaZ, float);
DECLARE_SOA_COLUMN(DcaXYZ, dcaXYZ, float);
} // namespace jtracktagdca

DECLARE_SOA_TABLE(JTracksTagDca, "AOD", "JTrackTagDcas",
                  jtracktagdca::DcaX,
                  jtracktagdca::DcaY,
                  jtracktagdca::DcaXY,
                  jtracktagdca::DcaZ,
                  jtracktagdca::DcaXYZ);

using JTrackTagDca = JTracksTagDca::iterator;

namespace jtracktagdcacov
{
DECLARE_SOA_COLUMN(SigmaDcaXY2, sigmaDcaXY2, float);
DECLARE_SOA_COLUMN(SigmaDcaZ2, sigmaDcaZ2, float);
DECLARE_SOA_COLUMN(SigmaDcaXYZ2, sigmaDcaXYZ2, float);
} // namespace jtracktagdcacov

DECLARE_SOA_TABLE(JTracksTagDcaCov, "AOD", "JTracksTDC",
                  jtracktagdcacov::SigmaDcaXY2,
                  jtracktagdcacov::SigmaDcaZ2,
                  jtracktagdcacov::SigmaDcaXYZ2);

using JTrackTagDcaCov = JTracksTagDcaCov::iterator;

namespace jtracktagsub
{
DECLARE_SOA_COLUMN(GeoSign, geoSign, int);
} // namespace jtracktagsub

DECLARE_SOA_TABLE(JTracksTagSub, "AOD", "JTracksTagSub",
                  jtracktagsub::GeoSign);

using JTrackTagSub = JTracksTagSub::iterator;

using JetTracksTagTC = soa::Join<JTracksTagDca, JTracksTagDcaCov>;
using JetTrackTagTC = JetTracksTagTC::iterator;

// Defines the tagger jet table definition
#define JETTAGGING_TABLE_DEF(_jet_type_, _name_, _description_)  \
  namespace _name_##tagging                                      \
  {                                                              \
    DECLARE_SOA_COLUMN(Origin, origin, int);                     \
    DECLARE_SOA_COLUMN(Algorithm1, algorithm1, int);             \
    DECLARE_SOA_COLUMN(Algorithm2, algorithm2, int);             \
    DECLARE_SOA_COLUMN(Algorithm3, algorithm3, int);             \
  }                                                              \
  DECLARE_SOA_TABLE(_jet_type_##Tags, "AOD", _description_ "JT", \
                    _name_##tagging::Origin,                     \
                    _name_##tagging::Algorithm1,                 \
                    _name_##tagging::Algorithm2,                 \
                    _name_##tagging::Algorithm3);

// Defines tagger jet extention
#define STOREDJETTAG_NAME_DEF(_jet_type_, _name_, _description_)            \
  namespace _name_##storedjettag                                            \
  {                                                                         \
    DECLARE_SOA_COLUMN(JetFlavour, jetFlavour, float);                      \
    DECLARE_SOA_COLUMN(JetPt, jetPt, float);                                \
    DECLARE_SOA_COLUMN(JetEta, jetEta, float);                              \
    DECLARE_SOA_COLUMN(JetPhi, jetPhi, float);                              \
    DECLARE_SOA_COLUMN(VecGeoSign, vecGeoSign, std::vector<int>);           \
    DECLARE_SOA_COLUMN(VecTrackPt, vecTrackPt, std::vector<float>);         \
    DECLARE_SOA_COLUMN(VecTrackEta, vecTrackEta, std::vector<float>);       \
    DECLARE_SOA_COLUMN(VecTrackPhi, vecTrackPhi, std::vector<float>);       \
    DECLARE_SOA_COLUMN(VecSignedIP2D, vecSignedIP2D, std::vector<float>);   \
    DECLARE_SOA_COLUMN(VecSignedIP2Ds, vecSignedIP2Ds, std::vector<float>); \
    DECLARE_SOA_COLUMN(VecSignedIP3D, vecSignedIP3D, std::vector<float>);   \
    DECLARE_SOA_COLUMN(VecSignedIP3Ds, vecSignedIP3Ds, std::vector<float>); \
  } // namespace _name_##storedjettag

#define STOREDJETTAG_TABLE_DEF(_jet_type_, _name_, _description_)         \
  DECLARE_SOA_TABLE(Stored##_jet_type_##Tags, "AOD", _description_ "SJT", \
                    _name_##storedjettag::JetFlavour,                     \
                    _name_##storedjettag::JetPt,                          \
                    _name_##storedjettag::JetEta,                         \
                    _name_##storedjettag::JetPhi,                         \
                    _name_##storedjettag::VecGeoSign,                     \
                    _name_##storedjettag::VecTrackPt,                     \
                    _name_##storedjettag::VecTrackEta,                    \
                    _name_##storedjettag::VecTrackPhi,                    \
                    _name_##storedjettag::VecSignedIP2D,                  \
                    _name_##storedjettag::VecSignedIP2Ds,                 \
                    _name_##storedjettag::VecSignedIP3D,                  \
                    _name_##storedjettag::VecSignedIP3Ds);

// combine definition of tables for tagger jets
#define JETTAGGING_TABLES_DEF(_jet_type_, _track_type_, _description_)                                        \
  JETTAGGING_TABLE_DEF(_jet_type_##Jet, _jet_type_##jet, _description_);                                      \
  using _jet_type_##Jet##Tag = _jet_type_##Jet##Tag##s::iterator;                                             \
  STOREDJETTAG_NAME_DEF(_jet_type_##Jet, _jet_type_##jet, _description_)                                      \
  STOREDJETTAG_TABLE_DEF(_jet_type_##Jet, _jet_type_##jet, _description_)                                     \
  using Stored##_jet_type_##Jet##Tag = Stored##_jet_type_##Jet##Tags::iterator;                               \
  JETTAGGING_TABLE_DEF(_jet_type_##MCDetectorLevelJet, _jet_type_##mcdetectorleveljet, _description_ "MCD");  \
  using _jet_type_##MCDetectorLevelJet##Tag = _jet_type_##MCDetectorLevelJet##Tag##s::iterator;               \
  STOREDJETTAG_NAME_DEF(_jet_type_##MCDetectorLevelJet, _jet_type_##mcdetectorleveljet, _description_ "MCD")  \
  STOREDJETTAG_TABLE_DEF(_jet_type_##MCDetectorLevelJet, _jet_type_##mcdetectorleveljet, _description_ "MCD") \
  using Stored##_jet_type_##MCDetectorLevelJet##Tag = Stored##_jet_type_##MCDetectorLevelJet##Tags::iterator; \
  JETTAGGING_TABLE_DEF(_jet_type_##MCParticleLevelJet, _jet_type_##mcparticleleveljet, _description_ "MCP");  \
  using _jet_type_##MCParticleLevelJet##Tag = _jet_type_##MCParticleLevelJet##Tag##s::iterator;

JETTAGGING_TABLES_DEF(Charged, JTrack, "C");
JETTAGGING_TABLES_DEF(Full, JTrack, "F");
JETTAGGING_TABLES_DEF(Neutral, JTrack, "N");
JETTAGGING_TABLES_DEF(D0Charged, JTrack, "D0");
JETTAGGING_TABLES_DEF(LcCharged, JTrack, "Lc");
JETTAGGING_TABLES_DEF(BplusCharged, JTrack, "BPL");

} // namespace o2::aod

#endif // PWGJE_DATAMODEL_JETTAGGING_H_
