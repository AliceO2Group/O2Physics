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

using namespace o2::analysis;

namespace o2::aod
{

namespace jtracktag
{
DECLARE_SOA_COLUMN(Pt, pt, float);
}

DECLARE_SOA_TABLE(JTrackTags, "AOD", "JTrackTags",
                  jtracktag::Pt);

namespace jtracktagdca
{
DECLARE_SOA_COLUMN(DcaX, dcaX, float);
DECLARE_SOA_COLUMN(DcaY, dcaY, float);
DECLARE_SOA_COLUMN(DcaXY, dcaXY, float);
DECLARE_SOA_COLUMN(DcaZ, dcaZ, float);
DECLARE_SOA_COLUMN(DcaXYZ, dcaXYZ, float);
} // namespace jtracktagdca

DECLARE_SOA_TABLE(JTrackTagDcas, "AOD", "JTrackTagDcas",
                  jtracktagdca::DcaX,
                  jtracktagdca::DcaY,
                  jtracktagdca::DcaXY,
                  jtracktagdca::DcaZ,
                  jtracktagdca::DcaXYZ);

using JTrackTagDca = JTrackTagDcas::iterator;

namespace jtracktagdcacov
{
DECLARE_SOA_COLUMN(SigmaDcaXY2, sigmaDcaXY2, float);
DECLARE_SOA_COLUMN(SigmaDcaZ2, sigmaDcaZ2, float);
DECLARE_SOA_COLUMN(SigmaDcaXYZ2, sigmaDcaXYZ2, float);
} // namespace jtracktagdcacov

DECLARE_SOA_TABLE(JTrackTagDcaCovs, "AOD", "JTracksTDCs",
                  jtracktagdcacov::SigmaDcaXY2,
                  jtracktagdcacov::SigmaDcaZ2,
                  jtracktagdcacov::SigmaDcaXYZ2);

using JTrackTagDcaCov = JTrackTagDcaCovs::iterator;

namespace jtracktagext
{
DECLARE_SOA_COLUMN(GeoSign, geoSign, int);
}

DECLARE_SOA_TABLE(JTrackTagExts, "AOD", "JTracksExts",
                  jtracktagext::GeoSign);

namespace jtracktagsub
{
// Jet index column will be added in the macro
DECLARE_SOA_COLUMN(GeoSign, geoSign, int);
} // namespace jtracktagsub

namespace constituenttc
{
DECLARE_SOA_COLUMN(JetFlavour, jetFlavour, float);
DECLARE_SOA_COLUMN(JetPt, jetPt, float);
DECLARE_SOA_COLUMN(JetEta, jetEta, float);
DECLARE_SOA_COLUMN(JetPhi, jetPhi, float);
DECLARE_SOA_COLUMN(VecGeoSign, vecGeoSign, std::vector<int>);
DECLARE_SOA_COLUMN(VecTrackPt, vecTrackPt, std::vector<float>);
DECLARE_SOA_COLUMN(VecTrackEta, vecTrackEta, std::vector<float>);
DECLARE_SOA_COLUMN(VecTrackPhi, vecTrackPhi, std::vector<float>);
DECLARE_SOA_COLUMN(VecSignedIP2D, vecSignedIP2D, std::vector<float>);
DECLARE_SOA_COLUMN(VecSignedIP2Ds, vecSignedIP2Ds, std::vector<float>);
DECLARE_SOA_COLUMN(VecSignedIP3D, vecSignedIP3D, std::vector<float>);
DECLARE_SOA_COLUMN(VecSignedIP3Ds, vecSignedIP3Ds, std::vector<float>);
} // namespace constituenttc

// Defines the jet substrcuture table definition
#define JETTAGGING_TABLE_DEF(_jet_type_, _name_, _description_)  \
  namespace _name_##tagging                                      \
  {                                                              \
    DECLARE_SOA_COLUMN(Origin, origin, int);                     \
    DECLARE_SOA_COLUMN(Algorithm1, algorithm1, int);             \
    DECLARE_SOA_COLUMN(Algorithm2, algorithm2, int);             \
    DECLARE_SOA_COLUMN(Algorithm3, algorithm3, int);             \
  }                                                              \
  DECLARE_SOA_TABLE(_jet_type_##Tags, "AOD", _description_ "TJ", \
                    _name_##tagging::Origin,                     \
                    _name_##tagging::Algorithm1,                 \
                    _name_##tagging::Algorithm2,                 \
                    _name_##tagging::Algorithm3);

#define DECLARE_TAGTRACK_TABLE(_jet_type_, _name_, _track_type_, _description_)  \
  namespace _name_##tag##constituents                                            \
  {                                                                              \
    DECLARE_SOA_INDEX_COLUMN(_jet_type_##Tag, jettags);                          \
    DECLARE_SOA_ARRAY_INDEX_COLUMN(_track_type_##Tag, tracktags);                \
  }                                                                              \
  DECLARE_SOA_TABLE(_jet_type_##Tag##Constituents, "AOD", _description_ "CONST", \
                    _name_##tag##constituents::_jet_type_##Tag##Id,              \
                    _name_##tag##constituents::_track_type_##Tag##Ids);

// Defines tagger track extention
#define DECLARE_TAGTRACK_TC_TABLE(_jet_type_, _name_, _description_)              \
  DECLARE_SOA_TABLE(_jet_type_##Tag##Constituent##TCs, "AOD", _description_ "TC", \
                    _name_##tag##constituents::_jet_type_##Tag##Id,               \
                    constituent##tc::JetFlavour,                                  \
                    constituent##tc::JetPt,                                       \
                    constituent##tc::JetEta,                                      \
                    constituent##tc::JetPhi,                                      \
                    constituent##tc::VecGeoSign,                                  \
                    constituent##tc::VecTrackPt,                                  \
                    constituent##tc::VecTrackEta,                                 \
                    constituent##tc::VecTrackPhi,                                 \
                    constituent##tc::VecSignedIP2D,                               \
                    constituent##tc::VecSignedIP2Ds,                              \
                    constituent##tc::VecSignedIP3D,                               \
                    constituent##tc::VecSignedIP3Ds);

// combine definition of tables for jets, constituents, and substructure
#define JETTAGGING_TABLES_DEF(_jet_type_, _track_type_, _description_)                                                          \
  JETTAGGING_TABLE_DEF(_jet_type_##Jet, _jet_type_##jet, _description_);                                                        \
  using _jet_type_##Jet##Tag = _jet_type_##Jet##Tag##s::iterator;                                                               \
  DECLARE_TAGTRACK_TABLE(_jet_type_##Jet, _jet_type_##jet, _track_type_, _description_ "T")                                     \
  using _jet_type_##Jet##Tag##Constituent = _jet_type_##Jet##Tag##Constituents::iterator;                                       \
  DECLARE_TAGTRACK_TC_TABLE(_jet_type_##Jet, _jet_type_##jet, _description_)                                                    \
  using _jet_type_##Jet##Tag##Constituent##TC = _jet_type_##Jet##Tag##Constituent##TCs::iterator;                               \
  JETTAGGING_TABLE_DEF(_jet_type_##MCDetectorLevelJet, _jet_type_##mcdetectorleveljet, _description_ "MCD");                    \
  using _jet_type_##MCDetectorLevelJet##Tag = _jet_type_##MCDetectorLevelJet##Tag##s::iterator;                                 \
  DECLARE_TAGTRACK_TABLE(_jet_type_##MCDetectorLevelJet, _jet_type_##mcdetectorleveljet, _track_type_, _description_ "DTAG")    \
  using _jet_type_##MCDetectorLevel##Jet##Tag##Constituent = _jet_type_##MCDetectorLevel##Jet##Tag##Constituents::iterator;     \
  DECLARE_TAGTRACK_TC_TABLE(_jet_type_##MCDetectorLevelJet, _jet_type_##mcdetectorleveljet, _description_ "MCD")                \
  using _jet_type_##MCDetectorLevelJet##Tag##Constituent##TC = _jet_type_##MCDetectorLevelJet##Tag##Constituent##TCs::iterator; \
  JETTAGGING_TABLE_DEF(_jet_type_##MCParticleLevelJet, _jet_type_##mcparticleleveljet, _description_ "MCP");                    \
  using _jet_type_##MCParticleLevelJet##Tag = _jet_type_##MCParticleLevelJet##Tag##s::iterator;                                 \
  DECLARE_TAGTRACK_TABLE(_jet_type_##MCParticleLevelJet, _jet_type_##mcparticleleveljet, _track_type_, _description_ "PTAG")    \
  using _jet_type_##MCParticleLevel##Jet##Tag##Constituent = _jet_type_##MCParticleLevel##Jet##Tag##Constituents::iterator;

JETTAGGING_TABLES_DEF(Charged, JTrack, "C");
JETTAGGING_TABLES_DEF(Full, JTrack, "F");
JETTAGGING_TABLES_DEF(Neutral, JTrack, "N");
JETTAGGING_TABLES_DEF(D0Charged, JTrack, "D0");
JETTAGGING_TABLES_DEF(LcCharged, JTrack, "Lc");
JETTAGGING_TABLES_DEF(BplusCharged, JTrack, "BPL");

} // namespace o2::aod

#endif // PWGJE_DATAMODEL_JETTAGGING_H_
