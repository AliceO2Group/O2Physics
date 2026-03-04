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

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/Qvectors.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CommonConstants/MathConstants.h>
#include <Framework/AnalysisDataModel.h>

#include <Rtypes.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <string>
#include <unordered_map>
#include <vector>

#ifndef PWGEM_DILEPTON_DATAMODEL_DILEPTONTABLES_H_
#define PWGEM_DILEPTON_DATAMODEL_DILEPTONTABLES_H_

namespace o2::aod
{

namespace pwgem::dilepton::swt
{
enum class swtAliases : int { // software trigger aliases for EM
  kHighTrackMult = 0,
  kHighFt0cFv0Mult,
  kSingleE,
  kLMeeIMR,
  kLMeeHMR,
  kDiElectron,
  kSingleMuLow,
  kSingleMuHigh,
  kDiMuon,
  kNaliases
};

const std::unordered_map<std::string, int> aliasLabels = {
  {"fHighTrackMult", static_cast<int>(o2::aod::pwgem::dilepton::swt::swtAliases::kHighTrackMult)},
  {"fHighFt0cFv0Mult", static_cast<int>(o2::aod::pwgem::dilepton::swt::swtAliases::kHighFt0cFv0Mult)},
  {"fSingleE", static_cast<int>(o2::aod::pwgem::dilepton::swt::swtAliases::kSingleE)},
  {"fLMeeIMR", static_cast<int>(o2::aod::pwgem::dilepton::swt::swtAliases::kLMeeIMR)},
  {"fLMeeHMR", static_cast<int>(o2::aod::pwgem::dilepton::swt::swtAliases::kLMeeHMR)},
  {"fDiElectron", static_cast<int>(o2::aod::pwgem::dilepton::swt::swtAliases::kDiElectron)},
  {"fSingleMuLow", static_cast<int>(o2::aod::pwgem::dilepton::swt::swtAliases::kSingleMuLow)},
  {"fSingleMuHigh", static_cast<int>(o2::aod::pwgem::dilepton::swt::swtAliases::kSingleMuHigh)},
  {"fDiMuon", static_cast<int>(o2::aod::pwgem::dilepton::swt::swtAliases::kDiMuon)},
};
} // namespace pwgem::dilepton::swt

namespace emevsel
{
// Event selection criteria. See O2Physics/Common/CCDB/EventSelectionParams.h
enum EventSelectionFlags {
  kIsTriggerTVX = 0,          // FT0 vertex (acceptable FT0C-FT0A time difference) at trigger level
  kNoITSROFrameBorder,        // bunch crossing is far from ITS RO Frame border
  kNoTimeFrameBorder,         // bunch crossing is far from Time Frame borders
  kNoSameBunchPileup,         // reject collisions in case of pileup with another collision in the same foundBC
  kIsGoodZvtxFT0vsPV,         // small difference between z-vertex from PV and from FT0
  kIsVertexITSTPC,            // at least one ITS-TPC track (reject vertices built from ITS-only tracks)
  kIsVertexTOFmatched,        // at least one of vertex contributors is matched to TOF
  kIsVertexTRDmatched,        // at least one of vertex contributors is matched to TRD
  kNoCollInTimeRangeNarrow,   // no other collisions in specified time range (narrower than Strict)
  kNoCollInTimeRangeStrict,   // no other collisions in specified time range
  kNoCollInTimeRangeStandard, // no other collisions in specified time range with per-collision multiplicity above threshold
  kNoCollInRofStrict,         // no other collisions in this Readout Frame
  kNoCollInRofStandard,       // no other collisions in this Readout Frame with per-collision multiplicity above threshold
  kNoHighMultCollInPrevRof,   // veto an event if FT0C amplitude in previous ITS ROF is above threshold
  kIsGoodITSLayer3,           // number of inactive chips on ITS layer 3 is below maximum allowed value
  kIsGoodITSLayer0123,        // numbers of inactive chips on ITS layers 0-3 are below maximum allowed values
  kIsGoodITSLayersAll,        // numbers of inactive chips on all ITS layers are below maximum allowed values
  kNsel                       // counter
};

DECLARE_SOA_BITMAP_COLUMN(Selection, selection, 32); //! Bitmask of selection flags
DECLARE_SOA_DYNAMIC_COLUMN(Sel8, sel8, [](uint32_t selection_bit) -> bool { return (selection_bit & BIT(o2::aod::emevsel::kIsTriggerTVX)) && (selection_bit & BIT(o2::aod::emevsel::kNoTimeFrameBorder)) && (selection_bit & BIT(o2::aod::emevsel::kNoITSROFrameBorder)); });

template <typename TBC>
uint32_t reduceSelectionBit(TBC const& bc)
{
  // input should be aod::BcSels or aod::EvSels.
  uint32_t bitMap = 0;
  if (bc.selection_bit(o2::aod::evsel::kIsTriggerTVX)) {
    SETBIT(bitMap, o2::aod::emevsel::kIsTriggerTVX);
  }
  if (bc.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
    SETBIT(bitMap, o2::aod::emevsel::kNoTimeFrameBorder);
  }
  if (bc.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
    SETBIT(bitMap, o2::aod::emevsel::kNoITSROFrameBorder);
  }
  if (bc.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
    SETBIT(bitMap, o2::aod::emevsel::kNoSameBunchPileup);
  }
  if (bc.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
    SETBIT(bitMap, o2::aod::emevsel::kIsGoodZvtxFT0vsPV);
  }
  if (bc.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) {
    SETBIT(bitMap, o2::aod::emevsel::kIsVertexITSTPC);
  }
  if (bc.selection_bit(o2::aod::evsel::kIsVertexTRDmatched)) {
    SETBIT(bitMap, o2::aod::emevsel::kIsVertexTRDmatched);
  }
  if (bc.selection_bit(o2::aod::evsel::kIsVertexTOFmatched)) {
    SETBIT(bitMap, o2::aod::emevsel::kIsVertexTOFmatched);
  }
  if (bc.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
    SETBIT(bitMap, o2::aod::emevsel::kNoCollInTimeRangeStandard);
  }
  if (bc.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStrict)) {
    SETBIT(bitMap, o2::aod::emevsel::kNoCollInTimeRangeStrict);
  }
  if (bc.selection_bit(o2::aod::evsel::kNoCollInTimeRangeNarrow)) {
    SETBIT(bitMap, o2::aod::emevsel::kNoCollInTimeRangeNarrow);
  }
  if (bc.selection_bit(o2::aod::evsel::kNoCollInRofStandard)) {
    SETBIT(bitMap, o2::aod::emevsel::kNoCollInRofStandard);
  }
  if (bc.selection_bit(o2::aod::evsel::kNoCollInRofStrict)) {
    SETBIT(bitMap, o2::aod::emevsel::kNoCollInRofStrict);
  }
  if (bc.selection_bit(o2::aod::evsel::kNoHighMultCollInPrevRof)) {
    SETBIT(bitMap, o2::aod::emevsel::kNoHighMultCollInPrevRof);
  }
  if (bc.selection_bit(o2::aod::evsel::kIsGoodITSLayer3)) {
    SETBIT(bitMap, o2::aod::emevsel::kIsGoodITSLayer3);
  }
  if (bc.selection_bit(o2::aod::evsel::kIsGoodITSLayer0123)) {
    SETBIT(bitMap, o2::aod::emevsel::kIsGoodITSLayer0123);
  }
  if (bc.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll)) {
    SETBIT(bitMap, o2::aod::emevsel::kIsGoodITSLayersAll);
  }
  return bitMap;
}

} // namespace emevsel

namespace emevent
{
DECLARE_SOA_COLUMN(CollisionId, collisionId, int);
DECLARE_SOA_BITMAP_COLUMN(SWTAliasTmp, swtaliastmp, 16);             //! Bitmask of fired trigger aliases (see above for definitions) to be join to aod::Collisions for skimming
DECLARE_SOA_BITMAP_COLUMN(SWTAlias, swtalias, 16);                   //! Bitmask of fired trigger aliases (see above for definitions) to be join to aod::EMEvents for analysis
DECLARE_SOA_COLUMN(NInspectedTVX, nInspectedTVX, uint64_t);          //! the number of inspected TVX bcs per run
DECLARE_SOA_COLUMN(NScalars, nScalers, std::vector<uint64_t>);       //! the number of triggered bcs before down scaling per run
DECLARE_SOA_COLUMN(NSelections, nSelections, std::vector<uint64_t>); //! the number of triggered bcs after down scaling per run
DECLARE_SOA_BITMAP_COLUMN(IsAnalyzed, isAnalyzed, 16);
DECLARE_SOA_BITMAP_COLUMN(IsAnalyzedToI, isAnalyzedToI, 16);
DECLARE_SOA_COLUMN(NeeULS, neeuls, int);
DECLARE_SOA_COLUMN(NeeLSpp, neelspp, int);
DECLARE_SOA_COLUMN(NeeLSmm, neelsmm, int);
DECLARE_SOA_COLUMN(Bz, bz, float);                                          //! kG
DECLARE_SOA_COLUMN(Q2xFT0M, q2xft0m, float);                                //! Qx for 2nd harmonics in FT0M
DECLARE_SOA_COLUMN(Q2yFT0M, q2yft0m, float);                                //! Qy for 2nd harmonics in FT0M
DECLARE_SOA_COLUMN(Q2xFT0A, q2xft0a, float);                                //! Qx for 2nd harmonics in FT0A (i.e. positive eta)
DECLARE_SOA_COLUMN(Q2yFT0A, q2yft0a, float);                                //! Qy for 2nd harmonics in FT0A (i.e. positive eta)
DECLARE_SOA_COLUMN(Q2xFT0C, q2xft0c, float);                                //! Qx for 2nd harmonics in FT0C (i.e. negative eta)
DECLARE_SOA_COLUMN(Q2yFT0C, q2yft0c, float);                                //! Qy for 2nd harmonics in FT0C (i.e. negative eta)
DECLARE_SOA_COLUMN(Q2xFV0A, q2xfv0a, float);                                //! Qx for 2nd harmonics in FV0A (i.e. positive eta)
DECLARE_SOA_COLUMN(Q2yFV0A, q2yfv0a, float);                                //! Qy for 2nd harmonics in FV0A (i.e. positive eta)
DECLARE_SOA_COLUMN(Q2xBPos, q2xbpos, float);                                //! Qx for 2nd harmonics in Barrel positive eta region
DECLARE_SOA_COLUMN(Q2yBPos, q2ybpos, float);                                //! Qy for 2nd harmonics in Barrel positive eta region
DECLARE_SOA_COLUMN(Q2xBNeg, q2xbneg, float);                                //! Qx for 2nd harmonics in Barrel negative eta region
DECLARE_SOA_COLUMN(Q2yBNeg, q2ybneg, float);                                //! Qy for 2nd harmonics in Barrel negative eta region
DECLARE_SOA_COLUMN(Q2xBTot, q2xbtot, float);                                //! Qx for 2nd harmonics in Barrel full eta region
DECLARE_SOA_COLUMN(Q2yBTot, q2ybtot, float);                                //! Qy for 2nd harmonics in Barrel full eta region
DECLARE_SOA_COLUMN(Q3xFT0M, q3xft0m, float);                                //! Qx for 3rd harmonics in FT0M
DECLARE_SOA_COLUMN(Q3yFT0M, q3yft0m, float);                                //! Qy for 3rd harmonics in FT0M
DECLARE_SOA_COLUMN(Q3xFT0A, q3xft0a, float);                                //! Qx for 3rd harmonics in FT0A (i.e. positive eta)
DECLARE_SOA_COLUMN(Q3yFT0A, q3yft0a, float);                                //! Qy for 3rd harmonics in FT0A (i.e. positive eta)
DECLARE_SOA_COLUMN(Q3xFT0C, q3xft0c, float);                                //! Qx for 3rd harmonics in FT0C (i.e. negative eta)
DECLARE_SOA_COLUMN(Q3yFT0C, q3yft0c, float);                                //! Qy for 3rd harmonics in FT0C (i.e. negative eta)
DECLARE_SOA_COLUMN(Q3xFV0A, q3xfv0a, float);                                //! Qx for 3rd harmonics in FV0A (i.e. positive eta)
DECLARE_SOA_COLUMN(Q3yFV0A, q3yfv0a, float);                                //! Qy for 3rd harmonics in FV0A (i.e. positive eta)
DECLARE_SOA_COLUMN(Q3xBPos, q3xbpos, float);                                //! Qx for 3rd harmonics in Barrel positive eta region
DECLARE_SOA_COLUMN(Q3yBPos, q3ybpos, float);                                //! Qy for 3rd harmonics in Barrel positive eta region
DECLARE_SOA_COLUMN(Q3xBNeg, q3xbneg, float);                                //! Qx for 3rd harmonics in Barrel negative eta region
DECLARE_SOA_COLUMN(Q3yBNeg, q3ybneg, float);                                //! Qy for 3rd harmonics in Barrel negative eta region
DECLARE_SOA_COLUMN(Q3xBTot, q3xbtot, float);                                //! Qx for 3rd harmonics in Barrel full eta region
DECLARE_SOA_COLUMN(Q3yBTot, q3ybtot, float);                                //! Qy for 3rd harmonics in Barrel full eta region
DECLARE_SOA_COLUMN(Q4xFT0M, q4xft0m, float);                                //! Qx for 4th harmonics in FT0M
DECLARE_SOA_COLUMN(Q4yFT0M, q4yft0m, float);                                //! Qy for 4th harmonics in FT0M
DECLARE_SOA_COLUMN(Q4xFT0A, q4xft0a, float);                                //! Qx for 4th harmonics in FT0A (i.e. positive eta)
DECLARE_SOA_COLUMN(Q4yFT0A, q4yft0a, float);                                //! Qy for 4th harmonics in FT0A (i.e. positive eta)
DECLARE_SOA_COLUMN(Q4xFT0C, q4xft0c, float);                                //! Qx for 4th harmonics in FT0C (i.e. negative eta)
DECLARE_SOA_COLUMN(Q4yFT0C, q4yft0c, float);                                //! Qy for 4th harmonics in FT0C (i.e. negative eta)
DECLARE_SOA_COLUMN(Q4xFV0A, q4xfv0a, float);                                //! Qx for 4th harmonics in FV0A (i.e. positive eta)
DECLARE_SOA_COLUMN(Q4yFV0A, q4yfv0a, float);                                //! Qy for 4th harmonics in FV0A (i.e. positive eta)
DECLARE_SOA_COLUMN(Q4xBPos, q4xbpos, float);                                //! Qx for 4th harmonics in Barrel positive eta region
DECLARE_SOA_COLUMN(Q4yBPos, q4ybpos, float);                                //! Qy for 4th harmonics in Barrel positive eta region
DECLARE_SOA_COLUMN(Q4xBNeg, q4xbneg, float);                                //! Qx for 4th harmonics in Barrel negative eta region
DECLARE_SOA_COLUMN(Q4yBNeg, q4ybneg, float);                                //! Qy for 4th harmonics in Barrel negative eta region
DECLARE_SOA_COLUMN(Q4xBTot, q4xbtot, float);                                //! Qx for 4th harmonics in Barrel full eta region
DECLARE_SOA_COLUMN(Q4yBTot, q4ybtot, float);                                //! Qy for 4th harmonics in Barrel full eta region
DECLARE_SOA_COLUMN(SpherocityPtWeighted, spherocity_ptweighted, float);     //! transverse spherocity
DECLARE_SOA_COLUMN(SpherocityPtUnWeighted, spherocity_ptunweighted, float); //! transverse spherocity
DECLARE_SOA_COLUMN(NtrackSpherocity, ntspherocity, int);
DECLARE_SOA_COLUMN(IsSelected, isSelected, bool);                //! MB event selection info
DECLARE_SOA_COLUMN(IsEoI, isEoI, bool);                          //! lepton or photon exists in MB event (not for CEFP)
DECLARE_SOA_COLUMN(PosX, posX, float);                           //! only for treeCreatetorML.cxx
DECLARE_SOA_COLUMN(PosY, posY, float);                           //! only for treeCreatetorML.cxx
DECLARE_SOA_COLUMN(PosZint16, posZint16, int16_t);               //! this is only to reduce data size
DECLARE_SOA_COLUMN(CentFT0Cuint16, centFT0Cuint16, uint16_t);    //! this is only to reduce data size
DECLARE_SOA_COLUMN(PosZint8, posZint8, int8_t);                  //! this is only to reduce data size
DECLARE_SOA_COLUMN(CentFT0Cuint8, centFT0Cuint8, uint8_t);       //! this is only to reduce data size
DECLARE_SOA_COLUMN(CentNTPVuint8, centNTPVuint8, uint8_t);       //! this is only to reduce data size
DECLARE_SOA_COLUMN(CentNGlobaluint8, centNGlobaluint8, uint8_t); //! this is only to reduce data size
DECLARE_SOA_COLUMN(CentFT0Muint8, centFT0Muint8, uint8_t);       //! this is only to reduce data size
DECLARE_SOA_COLUMN(CentFT0Auint8, centFT0Auint8, uint8_t);       //! this is only to reduce data size
// DECLARE_SOA_COLUMN(CentFT0Cuint8, centFT0Cuint8, uint8_t);       //! this is only to reduce data size
// DECLARE_SOA_COLUMN(CentNTPVuint8, centNTPVuint8, uint8_t);       //! this is only to reduce data size
// DECLARE_SOA_COLUMN(CentNGlobaluint8, centNGlobaluint8, uint8_t); //! this is only to reduce data size

DECLARE_SOA_DYNAMIC_COLUMN(PosZ, posZ, [](int16_t posZint16) -> float { return (posZint16 < 0 ? std::nextafter(posZint16 * 0.01f, -std::numeric_limits<float>::infinity()) : std::nextafter(posZint16 * 0.01f, std::numeric_limits<float>::infinity())); }); //! poZ is multiplied by 100 in createEMEventDileton.cxx
DECLARE_SOA_DYNAMIC_COLUMN(CentFT0C, centFT0C, [](uint16_t centuint16) -> float { return std::nextafter(centuint16 * 0.002f, std::numeric_limits<float>::infinity()); });                                                                                    //! centrality is multiplied by 500 in createEMEventDilepton.cxx
DECLARE_SOA_DYNAMIC_COLUMN(Sel8, sel8, [](uint64_t selection_bit) -> bool { return (selection_bit & BIT(o2::aod::evsel::kIsTriggerTVX)) && (selection_bit & BIT(o2::aod::evsel::kNoTimeFrameBorder)) && (selection_bit & BIT(o2::aod::evsel::kNoITSROFrameBorder)); });
DECLARE_SOA_DYNAMIC_COLUMN(EP2FT0M, ep2ft0m, [](float q2x, float q2y) -> float { return std::atan2(q2y, q2x) / 2.0; });
DECLARE_SOA_DYNAMIC_COLUMN(EP2FT0A, ep2ft0a, [](float q2x, float q2y) -> float { return std::atan2(q2y, q2x) / 2.0; });
DECLARE_SOA_DYNAMIC_COLUMN(EP2FT0C, ep2ft0c, [](float q2x, float q2y) -> float { return std::atan2(q2y, q2x) / 2.0; });
DECLARE_SOA_DYNAMIC_COLUMN(EP2FV0A, ep2fv0a, [](float q2x, float q2y) -> float { return std::atan2(q2y, q2x) / 2.0; });
DECLARE_SOA_DYNAMIC_COLUMN(EP2BPos, ep2bpos, [](float q2x, float q2y) -> float { return std::atan2(q2y, q2x) / 2.0; });
DECLARE_SOA_DYNAMIC_COLUMN(EP2BNeg, ep2bneg, [](float q2x, float q2y) -> float { return std::atan2(q2y, q2x) / 2.0; });
DECLARE_SOA_DYNAMIC_COLUMN(EP2BTot, ep2btot, [](float q2x, float q2y) -> float { return std::atan2(q2y, q2x) / 2.0; });
DECLARE_SOA_DYNAMIC_COLUMN(EP3FT0M, ep3ft0m, [](float q3x, float q3y) -> float { return std::atan2(q3y, q3x) / 3.0; });
DECLARE_SOA_DYNAMIC_COLUMN(EP3FT0A, ep3ft0a, [](float q3x, float q3y) -> float { return std::atan2(q3y, q3x) / 3.0; });
DECLARE_SOA_DYNAMIC_COLUMN(EP3FT0C, ep3ft0c, [](float q3x, float q3y) -> float { return std::atan2(q3y, q3x) / 3.0; });
DECLARE_SOA_DYNAMIC_COLUMN(EP3FV0A, ep3fv0a, [](float q3x, float q3y) -> float { return std::atan2(q3y, q3x) / 3.0; });
DECLARE_SOA_DYNAMIC_COLUMN(EP3BPos, ep3bpos, [](float q3x, float q3y) -> float { return std::atan2(q3y, q3x) / 3.0; });
DECLARE_SOA_DYNAMIC_COLUMN(EP3BNeg, ep3bneg, [](float q3x, float q3y) -> float { return std::atan2(q3y, q3x) / 3.0; });
DECLARE_SOA_DYNAMIC_COLUMN(EP3BTot, ep3btot, [](float q3x, float q3y) -> float { return std::atan2(q3y, q3x) / 3.0; });
DECLARE_SOA_DYNAMIC_COLUMN(EP4FT0M, ep4ft0m, [](float q4x, float q4y) -> float { return std::atan2(q4y, q4x) / 4.0; });
DECLARE_SOA_DYNAMIC_COLUMN(EP4FT0A, ep4ft0a, [](float q4x, float q4y) -> float { return std::atan2(q4y, q4x) / 4.0; });
DECLARE_SOA_DYNAMIC_COLUMN(EP4FT0C, ep4ft0c, [](float q4x, float q4y) -> float { return std::atan2(q4y, q4x) / 4.0; });
DECLARE_SOA_DYNAMIC_COLUMN(EP4FV0A, ep4fv0a, [](float q4x, float q4y) -> float { return std::atan2(q4y, q4x) / 4.0; });
DECLARE_SOA_DYNAMIC_COLUMN(EP4BPos, ep4bpos, [](float q4x, float q4y) -> float { return std::atan2(q4y, q4x) / 4.0; });
DECLARE_SOA_DYNAMIC_COLUMN(EP4BNeg, ep4bneg, [](float q4x, float q4y) -> float { return std::atan2(q4y, q4x) / 4.0; });
DECLARE_SOA_DYNAMIC_COLUMN(EP4BTot, ep4btot, [](float q4x, float q4y) -> float { return std::atan2(q4y, q4x) / 4.0; });
} // namespace emevent

namespace emeventnorm
{
DECLARE_SOA_DYNAMIC_COLUMN(PosZ, posZ, [](int8_t posZint8) -> float { return (posZint8 < 0 ? std::nextafter(posZint8 * 0.5f, -std::numeric_limits<float>::infinity()) : std::nextafter(posZint8 * 0.5f, std::numeric_limits<float>::infinity())); });                     //! posZ is multiplied by 2 in createEMEventDileton.cxx
DECLARE_SOA_DYNAMIC_COLUMN(CentFT0M, centFT0M, [](uint8_t centuint8) -> float { return centuint8 < 100 ? std::nextafter(centuint8 * 0.01f, std::numeric_limits<float>::infinity()) : std::nextafter(centuint8 - 110.f, std::numeric_limits<float>::infinity()); });       //! centrality is multiplied by 100 in createEMEventDilepton.cxx
DECLARE_SOA_DYNAMIC_COLUMN(CentFT0A, centFT0A, [](uint8_t centuint8) -> float { return centuint8 < 100 ? std::nextafter(centuint8 * 0.01f, std::numeric_limits<float>::infinity()) : std::nextafter(centuint8 - 110.f, std::numeric_limits<float>::infinity()); });       //! centrality is multiplied by 100 in createEMEventDilepton.cxx
DECLARE_SOA_DYNAMIC_COLUMN(CentFT0C, centFT0C, [](uint8_t centuint8) -> float { return centuint8 < 100 ? std::nextafter(centuint8 * 0.01f, std::numeric_limits<float>::infinity()) : std::nextafter(centuint8 - 110.f, std::numeric_limits<float>::infinity()); });       //! centrality is multiplied by 100 in createEMEventDilepton.cxx
DECLARE_SOA_DYNAMIC_COLUMN(CentNTPV, centNTPV, [](uint8_t centuint8) -> float { return centuint8 < 100 ? std::nextafter(centuint8 * 0.01f, std::numeric_limits<float>::infinity()) : std::nextafter(centuint8 - 110.f, std::numeric_limits<float>::infinity()); });       //! centrality is multiplied by 100 in createEMEventDilepton.cxx
DECLARE_SOA_DYNAMIC_COLUMN(CentNGlobal, centNGlobal, [](uint8_t centuint8) -> float { return centuint8 < 100 ? std::nextafter(centuint8 * 0.01f, std::numeric_limits<float>::infinity()) : std::nextafter(centuint8 - 110.f, std::numeric_limits<float>::infinity()); }); //! centrality is multiplied by 100 in createEMEventDilepton.cxx
} // namespace emeventnorm

// namespace emcent
// {
// DECLARE_SOA_COLUMN(CentFT0Muint8, centFT0Muint8, uint8_t);       //! this is only to reduce data size
// DECLARE_SOA_COLUMN(CentFT0Auint8, centFT0Auint8, uint8_t);       //! this is only to reduce data size
// DECLARE_SOA_COLUMN(CentFT0Cuint8, centFT0Cuint8, uint8_t);       //! this is only to reduce data size
// DECLARE_SOA_COLUMN(CentNTPVuint8, centNTPVuint8, uint8_t);       //! this is only to reduce data size
// DECLARE_SOA_COLUMN(CentNGlobaluint8, centNGlobaluint8, uint8_t); //! this is only to reduce data size
//
// DECLARE_SOA_EXPRESSION_COLUMN(CentFT0A, centFT0A, float, 1.f * centFT0Auint8);          // this must be inverse of calculation in createEMEventDilepton.cxx
// DECLARE_SOA_EXPRESSION_COLUMN(CentFT0M, centFT0M, float, 1.f * centFT0Muint8);          // this must be inverse of calculation in createEMEventDilepton.cxx
// DECLARE_SOA_EXPRESSION_COLUMN(CentFT0C, centFT0C, float, 1.f * centFT0Cuint8);          // this must be inverse of calculation in createEMEventDilepton.cxx
// DECLARE_SOA_EXPRESSION_COLUMN(CentNTPV, centNTPV, float, 1.f * centNTPVuint8);          // this must be inverse of calculation in createEMEventDilepton.cxx
// DECLARE_SOA_EXPRESSION_COLUMN(CentNGlobal, centNGlobal, float, 1.f * centNGlobaluint8); // this must be inverse of calculation in createEMEventDilepton.cxx
//
// // DECLARE_SOA_EXPRESSION_COLUMN(CentFT0A, centFT0A, float, (centFT0Auint8 < 100) ? std::nextafter((1.f * centFT0Auint8) / 100.f, std::numeric_limits<float>::infinity()) : std::nextafter((1.f * centFT0Auint8) - 110.f, std::numeric_limits<float>::infinity())); // this must be inverse of calculation in createEMEventDilepton.cxx
// // DECLARE_SOA_EXPRESSION_COLUMN(CentFT0M, centFT0M, float, (centFT0Muint8 < 100) ? std::nextafter((1.f * centFT0Muint8) / 100.f, std::numeric_limits<float>::infinity()) : std::nextafter((1.f * centFT0Muint8) - 110.f, std::numeric_limits<float>::infinity())); // this must be inverse of calculation in createEMEventDilepton.cxx
// // DECLARE_SOA_EXPRESSION_COLUMN(CentFT0C, centFT0C, float, (centFT0Cuint8 < 100) ? std::nextafter((1.f * centFT0Cuint8) / 100.f, std::numeric_limits<float>::infinity()) : std::nextafter((1.f * centFT0Cuint8) - 110.f, std::numeric_limits<float>::infinity())); // this must be inverse of calculation in createEMEventDilepton.cxx
// // DECLARE_SOA_EXPRESSION_COLUMN(CentNTPV, centNTPV, float, (centNTPVuint8 < 100) ? std::nextafter((1.f * centNTPVuint8) / 100.f, std::numeric_limits<float>::infinity()) : std::nextafter((1.f * centNTPVuint8) - 110.f, std::numeric_limits<float>::infinity())); // this must be inverse of calculation in createEMEventDilepton.cxx
// // DECLARE_SOA_EXPRESSION_COLUMN(CentNGlobal, centNGlobal, float, (centNGlobaluint8 < 100) ? std::nextafter((1.f * centNGlobaluint8) / 100.f, std::numeric_limits<float>::infinity()) : std::nextafter((1.f * centNGlobaluint8) - 110.f, std::numeric_limits<float>::infinity())); // this must be inverse of calculation in createEMEventDilepton.cxx
//
// } // namespace emcent

DECLARE_SOA_TABLE(EMBCs_000, "AOD", "EMBC", //! bc information for normalization
                  o2::soa::Index<>, evsel::Selection, evsel::Rct);

DECLARE_SOA_TABLE_VERSIONED(EMBCs_001, "AOD", "EMBC", 1, //! bc information for normalization
                            o2::soa::Index<>, emevsel::Selection, evsel::Rct);

using EMBCs = EMBCs_001;
using EMBC = EMBCs::iterator;

DECLARE_SOA_TABLE_VERSIONED(EMEvents_001, "AOD", "EMEVENT", 1, //!   Main event information table
                            o2::soa::Index<>, emevent::CollisionId, bc::RunNumber, bc::GlobalBC, evsel::Alias, evsel::Selection, timestamp::Timestamp,
                            collision::PosX, collision::PosY, collision::PosZ,
                            collision::NumContrib, evsel::NumTracksInTimeRange, evsel::SumAmpFT0CInTimeRange, emevent::Sel8<evsel::Selection>);

DECLARE_SOA_TABLE_VERSIONED(EMEvents_002, "AOD", "EMEVENT", 2, //!   Main event information table
                            o2::soa::Index<>, emevent::CollisionId, bc::RunNumber, bc::GlobalBC, evsel::Alias, evsel::Selection, evsel::Rct, timestamp::Timestamp,
                            collision::PosX, collision::PosY, collision::PosZ,
                            collision::NumContrib, evsel::NumTracksInTimeRange, evsel::SumAmpFT0CInTimeRange, emevent::Sel8<evsel::Selection>);

DECLARE_SOA_TABLE_VERSIONED(EMEvents_003, "AOD", "EMEVENT", 3, //!   Main event information table
                            o2::soa::Index<>, emevent::CollisionId, bc::RunNumber, bc::GlobalBC, evsel::Alias, evsel::Selection, evsel::Rct, timestamp::Timestamp,
                            collision::PosZ,
                            collision::NumContrib, evsel::NumTracksInTimeRange, evsel::SumAmpFT0CInTimeRange, emevent::Sel8<evsel::Selection>);

DECLARE_SOA_TABLE_VERSIONED(EMEvents_004, "AOD", "EMEVENT", 4, //!   Main event information table
                            o2::soa::Index<>, emevent::CollisionId, bc::RunNumber, bc::GlobalBC, evsel::Selection, evsel::Rct, timestamp::Timestamp,
                            collision::PosZ,
                            collision::NumContrib, evsel::NumTracksInTimeRange, evsel::SumAmpFT0CInTimeRange, emevent::Sel8<evsel::Selection>);

DECLARE_SOA_TABLE_VERSIONED(EMEvents_005, "AOD", "EMEVENT", 5, //!   Main event information table
                            o2::soa::Index<>, emevent::CollisionId, bc::RunNumber, bc::GlobalBC, emevsel::Selection, evsel::Rct, timestamp::Timestamp,
                            collision::PosZ,
                            collision::NumContrib, evsel::NumTracksInTimeRange, evsel::SumAmpFT0CInTimeRange, emevsel::Sel8<emevsel::Selection>);

using EMEvents = EMEvents_005;
using EMEvent = EMEvents::iterator;

DECLARE_SOA_TABLE_VERSIONED(EMEventsAlias_000, "AOD", "EMEVENTALIAS", 0, evsel::Alias) //! joinable to EMEvents
using EMEventsAlias = EMEventsAlias_000;
using EMEventAlias = EMEventsAlias::iterator;

DECLARE_SOA_TABLE(EMEventsXY, "AOD", "EMEVENTXY", collision::PosX, collision::PosY); // joinable to EMEvents
using EMEventXY = EMEventsXY::iterator;

DECLARE_SOA_TABLE(EMEventsCov, "AOD", "EMEVENTCOV", //! joinable to EMEvents
                  collision::CovXX, collision::CovXY, collision::CovXZ, collision::CovYY, collision::CovYZ, collision::CovZZ, collision::Chi2, o2::soa::Marker<1>);
using EMEventCov = EMEventsCov::iterator;

DECLARE_SOA_TABLE(EMEventsBz, "AOD", "EMEVENTBZ", emevent::Bz); // joinable to EMEvents
using EMEventBz = EMEventsBz::iterator;

DECLARE_SOA_TABLE(EMEventsMult_000, "AOD", "EMEVENTMULT", //!   event multiplicity table, joinable to EMEvents
                  mult::MultFT0A, mult::MultFT0C, mult::MultNTracksPV, mult::MultNTracksPVeta1, mult::MultNTracksPVetaHalf,
                  mult::IsInelGt0<mult::MultNTracksPVeta1>, mult::IsInelGt1<mult::MultNTracksPVeta1>, mult::MultFT0M<mult::MultFT0A, mult::MultFT0C>);

DECLARE_SOA_TABLE_VERSIONED(EMEventsMult_001, "AOD", "EMEVENTMULT", 1,           //!   event multiplicity table, joinable to EMEvents
                            mult::MultFT0A, mult::MultFT0C, mult::MultNTracksPV, /*mult::MultNTracksGlobal,*/
                            mult::MultFT0M<mult::MultFT0A, mult::MultFT0C>);

using EMEventsMult = EMEventsMult_001;
using EMEventMult = EMEventsMult::iterator;

DECLARE_SOA_TABLE(EMEventsCent_000, "AOD", "EMEVENTCENT", //!   event centrality table, joinable to EMEvents
                  cent::CentFT0M, cent::CentFT0A, cent::CentFT0C);

DECLARE_SOA_TABLE_VERSIONED(EMEventsCent_001, "AOD", "EMEVENTCENT", 1, //! event centrality table stored in AO2D
                            cent::CentFT0M, cent::CentFT0A, cent::CentFT0C, cent::CentNTPV /*, cent::CentNGlobal*/);

using EMEventsCent = EMEventsCent_001;
using EMEventCent = EMEventsCent::iterator;

// DECLARE_SOA_TABLE_VERSIONED(EMEventsCentBase_001, "AOD", "EMEVENTCENT", 1,  //! event centrality table stored in AO2D
//     emcent::CentFT0Muint8, emcent::CentFT0Auint8, emcent::CentFT0Cuint8, emcent::CentNTPVuint8, emcent::CentNGlobaluint8);
//
// using EMEventsCentBase = EMEventsCentBase_001;
// using EMEventCentBase = EMEventsCentBase::iterator;
//
// // Extended table with expression columns that can be used for Filter.
// DECLARE_SOA_EXTENDED_TABLE_USER(EMEventsCent, EMEventsCentBase, "EMCENTEXT",
//     emcent::CentFT0M, emcent::CentFT0A, emcent::CentFT0C, emcent::CentNTPV, emcent::CentNGlobal);
//
// using EMEventCent = EMEventsCent::iterator;

DECLARE_SOA_TABLE_VERSIONED(EMEventsQvec_000, "AOD", "EMEVENTQVEC", 0, //!   event q vector table, joinable to EMEvents
                            emevent::Q2xFT0M, emevent::Q2yFT0M, emevent::Q2xFT0A, emevent::Q2yFT0A, emevent::Q2xFT0C, emevent::Q2yFT0C,
                            emevent::Q2xBPos, emevent::Q2yBPos, emevent::Q2xBNeg, emevent::Q2yBNeg, emevent::Q2xBTot, emevent::Q2yBTot,
                            emevent::Q3xFT0M, emevent::Q3yFT0M, emevent::Q3xFT0A, emevent::Q3yFT0A, emevent::Q3xFT0C, emevent::Q3yFT0C,
                            emevent::Q3xBPos, emevent::Q3yBPos, emevent::Q3xBNeg, emevent::Q3yBNeg, emevent::Q3xBTot, emevent::Q3yBTot,

                            // Dynamic columns
                            emevent::EP2FT0M<emevent::Q2xFT0M, emevent::Q2yFT0M>,
                            emevent::EP2FT0A<emevent::Q2xFT0A, emevent::Q2yFT0A>,
                            emevent::EP2FT0C<emevent::Q2xFT0C, emevent::Q2yFT0C>,
                            emevent::EP2BPos<emevent::Q2xBPos, emevent::Q2yBPos>,
                            emevent::EP2BNeg<emevent::Q2xBNeg, emevent::Q2yBNeg>,
                            emevent::EP2BTot<emevent::Q2xBTot, emevent::Q2yBTot>,
                            emevent::EP3FT0M<emevent::Q3xFT0M, emevent::Q3yFT0M>,
                            emevent::EP3FT0A<emevent::Q3xFT0A, emevent::Q3yFT0A>,
                            emevent::EP3FT0C<emevent::Q3xFT0C, emevent::Q3yFT0C>,
                            emevent::EP3BPos<emevent::Q3xBPos, emevent::Q3yBPos>,
                            emevent::EP3BNeg<emevent::Q3xBNeg, emevent::Q3yBNeg>,
                            emevent::EP3BTot<emevent::Q3xBTot, emevent::Q3yBTot>);

DECLARE_SOA_TABLE_VERSIONED(EMEventsQvec_001, "AOD", "EMEVENTQVEC", 1, //!   Main event information table
                            emevent::Q2xFT0M, emevent::Q2yFT0M, emevent::Q2xFT0A, emevent::Q2yFT0A, emevent::Q2xFT0C, emevent::Q2yFT0C,
                            emevent::Q2xFV0A, emevent::Q2yFV0A,
                            emevent::Q2xBPos, emevent::Q2yBPos, emevent::Q2xBNeg, emevent::Q2yBNeg, emevent::Q2xBTot, emevent::Q2yBTot,
                            emevent::Q3xFT0M, emevent::Q3yFT0M, emevent::Q3xFT0A, emevent::Q3yFT0A, emevent::Q3xFT0C, emevent::Q3yFT0C,
                            emevent::Q3xFV0A, emevent::Q3yFV0A,
                            emevent::Q3xBPos, emevent::Q3yBPos, emevent::Q3xBNeg, emevent::Q3yBNeg, emevent::Q3xBTot, emevent::Q3yBTot,

                            // Dynamic columns
                            emevent::EP2FT0M<emevent::Q2xFT0M, emevent::Q2yFT0M>,
                            emevent::EP2FT0A<emevent::Q2xFT0A, emevent::Q2yFT0A>,
                            emevent::EP2FT0C<emevent::Q2xFT0C, emevent::Q2yFT0C>,
                            emevent::EP2FV0A<emevent::Q2xFV0A, emevent::Q2yFV0A>,
                            emevent::EP2BPos<emevent::Q2xBPos, emevent::Q2yBPos>,
                            emevent::EP2BNeg<emevent::Q2xBNeg, emevent::Q2yBNeg>,
                            emevent::EP2BTot<emevent::Q2xBTot, emevent::Q2yBTot>,
                            emevent::EP3FT0M<emevent::Q3xFT0M, emevent::Q3yFT0M>,
                            emevent::EP3FT0A<emevent::Q3xFT0A, emevent::Q3yFT0A>,
                            emevent::EP3FT0C<emevent::Q3xFT0C, emevent::Q3yFT0C>,
                            emevent::EP3FV0A<emevent::Q3xFV0A, emevent::Q3yFV0A>,
                            emevent::EP3BPos<emevent::Q3xBPos, emevent::Q3yBPos>,
                            emevent::EP3BNeg<emevent::Q3xBNeg, emevent::Q3yBNeg>,
                            emevent::EP3BTot<emevent::Q3xBTot, emevent::Q3yBTot>);

using EMEventsQvec = EMEventsQvec_001;
using EMEventQvec = EMEventsQvec::iterator;

DECLARE_SOA_TABLE_VERSIONED(EMEventsQvec2_000, "AOD", "EMEVENTQVEC2", 0, //!   Main event information table
                            emevent::Q2xFT0M, emevent::Q2yFT0M, emevent::Q2xFT0A, emevent::Q2yFT0A, emevent::Q2xFT0C, emevent::Q2yFT0C,
                            emevent::Q2xFV0A, emevent::Q2yFV0A,
                            emevent::Q2xBPos, emevent::Q2yBPos, emevent::Q2xBNeg, emevent::Q2yBNeg, emevent::Q2xBTot, emevent::Q2yBTot,
                            // Dynamic columns
                            emevent::EP2FT0M<emevent::Q2xFT0M, emevent::Q2yFT0M>,
                            emevent::EP2FT0A<emevent::Q2xFT0A, emevent::Q2yFT0A>,
                            emevent::EP2FT0C<emevent::Q2xFT0C, emevent::Q2yFT0C>,
                            emevent::EP2FV0A<emevent::Q2xFV0A, emevent::Q2yFV0A>,
                            emevent::EP2BPos<emevent::Q2xBPos, emevent::Q2yBPos>,
                            emevent::EP2BNeg<emevent::Q2xBNeg, emevent::Q2yBNeg>,
                            emevent::EP2BTot<emevent::Q2xBTot, emevent::Q2yBTot>);

using EMEventsQvec2 = EMEventsQvec2_000;
using EMEventQvec2 = EMEventsQvec2::iterator;

DECLARE_SOA_TABLE_VERSIONED(EMEventsQvec3_000, "AOD", "EMEVENTQVEC3", 0, //!   Main event information table
                            emevent::Q3xFT0M, emevent::Q3yFT0M, emevent::Q3xFT0A, emevent::Q3yFT0A, emevent::Q3xFT0C, emevent::Q3yFT0C,
                            emevent::Q3xFV0A, emevent::Q3yFV0A,
                            emevent::Q3xBPos, emevent::Q3yBPos, emevent::Q3xBNeg, emevent::Q3yBNeg, emevent::Q3xBTot, emevent::Q3yBTot,
                            // Dynamic columns
                            emevent::EP3FT0M<emevent::Q3xFT0M, emevent::Q3yFT0M>,
                            emevent::EP3FT0A<emevent::Q3xFT0A, emevent::Q3yFT0A>,
                            emevent::EP3FT0C<emevent::Q3xFT0C, emevent::Q3yFT0C>,
                            emevent::EP3FV0A<emevent::Q3xFV0A, emevent::Q3yFV0A>,
                            emevent::EP3BPos<emevent::Q3xBPos, emevent::Q3yBPos>,
                            emevent::EP3BNeg<emevent::Q3xBNeg, emevent::Q3yBNeg>,
                            emevent::EP3BTot<emevent::Q3xBTot, emevent::Q3yBTot>);

using EMEventsQvec3 = EMEventsQvec3_000;
using EMEventQvec3 = EMEventsQvec3::iterator;

DECLARE_SOA_TABLE(EMSWTriggerBits, "AOD", "EMSWTBIT", emevent::SWTAlias, o2::soa::Marker<1>); //! joinable to EMEvents
using EMSWTriggerBit = EMSWTriggerBits::iterator;

DECLARE_SOA_TABLE(EMSWTriggerInfos, "AOD", "EMSWTINFO", bc::RunNumber, emevent::NInspectedTVX, emevent::NScalars, emevent::NSelections, o2::soa::Marker<1>); //! independent table. Don't join anything.
using EMSWTriggerInfo = EMSWTriggerInfos::iterator;

DECLARE_SOA_TABLE(EMSWTriggerATCounters, "AOD", "EMSWTAT", emevent::IsAnalyzed, o2::soa::Marker<1>); //! independent table. Don't join anything.
using EMSWTriggerATCounter = EMSWTriggerATCounters::iterator;

DECLARE_SOA_TABLE(EMSWTriggerTOICounters, "AOD", "EMSWTTOI", emevent::IsAnalyzedToI, o2::soa::Marker<1>); //! independent table. Don't join anything.
using EMSWTriggerTOICounter = EMSWTriggerTOICounters::iterator;

DECLARE_SOA_TABLE(EMSWTriggerBitsTMP, "AOD", "EMSWTBITTMP", emevent::SWTAliasTmp, o2::soa::Marker<2>); //! joinable to aod::Collisions
using EMSWTriggerBitTMP = EMSWTriggerBitsTMP::iterator;

DECLARE_SOA_TABLE(EMSWTriggerInfosTMP, "AOD", "EMSWTINFOTMP", bc::RunNumber, emevent::NInspectedTVX, emevent::NScalars, emevent::NSelections, o2::soa::Marker<2>);
using EMSWTriggerInfoTMP = EMSWTriggerInfosTMP::iterator;

DECLARE_SOA_TABLE(EMSWTriggerATCountersTMP, "AOD", "EMSWTATTMP", emevent::IsAnalyzed, o2::soa::Marker<2>); //! independent table. Don't join anything.
using EMSWTriggerATCounterTMP = EMSWTriggerATCountersTMP::iterator;

DECLARE_SOA_TABLE(EMSWTriggerTOICountersTMP, "AOD", "EMSWTTOITMP", emevent::IsAnalyzedToI, o2::soa::Marker<2>); //! independent table. Don't join anything.
using EMSWTriggerTOICounterTMP = EMSWTriggerTOICountersTMP::iterator;

DECLARE_SOA_TABLE(EMEventsProperty, "AOD", "EMEVENTPROP", //! joinable to EMEvents
                  emevent::SpherocityPtWeighted, emevent::SpherocityPtUnWeighted, emevent::NtrackSpherocity);
using EMEventProperty = EMEventsProperty::iterator;

DECLARE_SOA_TABLE(EMEventsNee, "AOD", "EMEVENTNEE", emevent::NeeULS, emevent::NeeLSpp, emevent::NeeLSmm); // joinable to EMEvents or aod::Collisions
using EMEventNee = EMEventsNee::iterator;

DECLARE_SOA_TABLE(EMEvSels, "AOD", "EMEVSEL", //! joinable to aod::Collisions
                  emevent::IsSelected);
using EMEvSel = EMEvSels::iterator;

DECLARE_SOA_TABLE(EMEoIs, "AOD", "EMEOI", //! joinable to aod::Collisions in createEMEventDilepton.cxx
                  emevent::IsEoI);
using EMEoI = EMEoIs::iterator;

DECLARE_SOA_TABLE_VERSIONED(EMEventNormInfos_000, "AOD", "EMEVENTNORMINFO", 0, //! event information for normalization
                            o2::soa::Index<>, evsel::Alias, evsel::Selection, evsel::Rct, emevent::PosZint16, cent::CentFT0C,
                            emevent::Sel8<evsel::Selection>);

DECLARE_SOA_TABLE_VERSIONED(EMEventNormInfos_001, "AOD", "EMEVENTNORMINFO", 1, //! event information for normalization
                            o2::soa::Index<>, evsel::Selection, evsel::Rct, emevent::PosZint16, emevent::CentFT0Cuint16,
                            emevent::Sel8<evsel::Selection>, emevent::PosZ<emevent::PosZint16>, emevent::CentFT0C<emevent::CentFT0Cuint16>, o2::soa::Marker<1>);

DECLARE_SOA_TABLE_VERSIONED(EMEventNormInfos_002, "AOD", "EMEVENTNORMINFO", 2,                                                                         //! event information for normalization
                            emevsel::Selection, evsel::Rct, emevent::PosZint8, emevent::CentFT0Muint8, emevent::CentFT0Cuint8, emevent::CentNTPVuint8, /*emevent::CentNGlobaluint8,*/
                            emevsel::Sel8<emevsel::Selection>, emeventnorm::PosZ<emevent::PosZint8>, emeventnorm::CentFT0M<emevent::CentFT0Muint8>, emeventnorm::CentFT0C<emevent::CentFT0Cuint8>, emeventnorm::CentNTPV<emevent::CentNTPVuint8>, /*emeventnorm::CentNTPV<emevent::CentNGlobaluint8>,*/ o2::soa::Marker<1>);

using EMEventNormInfos = EMEventNormInfos_002;
using EMEventNormInfo = EMEventNormInfos::iterator;

namespace emmcevent
{
DECLARE_SOA_INDEX_COLUMN(EMEvent, mpemevent); //! most propable emeventId
DECLARE_SOA_COLUMN(McCollisionId, mcCollisionId, int);
} // namespace emmcevent

DECLARE_SOA_TABLE(EMMCEvents, "AOD", "EMMCEVENT", //!   MC event information table
                  o2::soa::Index<>, emmcevent::McCollisionId, mccollision::GeneratorsID,
                  mccollision::PosX, mccollision::PosY, mccollision::PosZ,
                  mccollision::ImpactParameter, mccollision::EventPlaneAngle,

                  // dynamic column
                  mccollision::GetGeneratorId<mccollision::GeneratorsID>,
                  mccollision::GetSubGeneratorId<mccollision::GeneratorsID>,
                  mccollision::GetSourceId<mccollision::GeneratorsID>);
using EMMCEvent = EMMCEvents::iterator;

DECLARE_SOA_TABLE(MostProbableEMEventIdsInMC, "AOD", "MPEMEVENTIDINMC", emmcevent::EMEventId); // To be joined with EMMCEvents table at analysis level.
using MostProbableEMEventIdInMC = MostProbableEMEventIdsInMC::iterator;

namespace emmceventlabel
{
DECLARE_SOA_INDEX_COLUMN(EMMCEvent, emmcevent); //! MC collision
DECLARE_SOA_COLUMN(McMask, mcMask, uint16_t);   //! Bit mask to indicate collision mismatches (bit ON means mismatch). Bit 15: indicates negative label
} // namespace emmceventlabel

DECLARE_SOA_TABLE(EMMCEventLabels, "AOD", "EMMCEVENTLABEL", //! Table joined to the EMEvents table containing the MC index
                  emmceventlabel::EMMCEventId, emmceventlabel::McMask);
using EMMCEventLabel = EMMCEventLabels::iterator;

namespace emmcparticle
{
DECLARE_SOA_INDEX_COLUMN(EMMCEvent, emmcevent);
DECLARE_SOA_SELF_ARRAY_INDEX_COLUMN(Mothers, mothers);     //! Mother tracks (possible empty) array. Iterate over mcParticle.mothers_as<aod::McParticles>())
DECLARE_SOA_SELF_ARRAY_INDEX_COLUMN(Daughters, daughters); //! Daughter tracks (possibly empty) array. Check for non-zero with mcParticle.has_daughters(). Iterate over mcParticle.daughters_as<aod::McParticles>())
DECLARE_SOA_DYNAMIC_COLUMN(Pt, pt, [](float px, float py) -> float { return RecoDecay::sqrtSumOfSquares(px, py); });
DECLARE_SOA_DYNAMIC_COLUMN(Eta, eta, [](float px, float py, float pz) -> float { return RecoDecay::eta(std::array{px, py, pz}); });
DECLARE_SOA_DYNAMIC_COLUMN(Phi, phi, [](float px, float py) -> float { return RecoDecay::phi(px, py); });
DECLARE_SOA_DYNAMIC_COLUMN(P, p, [](float px, float py, float pz) -> float { return RecoDecay::sqrtSumOfSquares(px, py, pz); });
DECLARE_SOA_DYNAMIC_COLUMN(Y, y, //! Particle rapidity
                           [](float pz, float e) -> float {
                             if ((e - pz) > static_cast<float>(1e-7)) {
                               return 0.5f * std::log((e + pz) / (e - pz));
                             } else {
                               return -999.0f;
                             }
                           });
} // namespace emmcparticle

// This table contains all MC truth tracks
DECLARE_SOA_TABLE_FULL(EMMCParticles_000, "EMMCParticles", "AOD", "EMMCPARTICLE", //!  MC track information (on disk)
                       o2::soa::Index<>, emmcparticle::EMMCEventId,
                       mcparticle::PdgCode, mcparticle::Flags,
                       emmcparticle::MothersIds, emmcparticle::DaughtersIds,
                       mcparticle::Px, mcparticle::Py, mcparticle::Pz, mcparticle::E,
                       mcparticle::Vx, mcparticle::Vy, mcparticle::Vz,

                       // dynamic column
                       emmcparticle::Pt<mcparticle::Px, mcparticle::Py>,
                       emmcparticle::Eta<mcparticle::Px, mcparticle::Py, mcparticle::Pz>,
                       emmcparticle::Phi<mcparticle::Px, mcparticle::Py>,
                       emmcparticle::P<mcparticle::Px, mcparticle::Py, mcparticle::Pz>,
                       emmcparticle::Y<mcparticle::Pz, mcparticle::E>,
                       mcparticle::ProducedByGenerator<mcparticle::Flags>,
                       mcparticle::FromBackgroundEvent<mcparticle::Flags>,
                       mcparticle::IsPhysicalPrimary<mcparticle::Flags>);

DECLARE_SOA_TABLE_VERSIONED(EMMCParticles_001, "AOD", "EMMCPARTICLE", 1, //!  MC track information (on disk)
                            o2::soa::Index<>, emmcparticle::EMMCEventId,
                            mcparticle::PdgCode, mcparticle::Flags, mcparticle::StatusCode,
                            emmcparticle::MothersIds, emmcparticle::DaughtersIds,
                            mcparticle::Px, mcparticle::Py, mcparticle::Pz, mcparticle::E,
                            mcparticle::Vx, mcparticle::Vy, mcparticle::Vz,

                            // dynamic column
                            emmcparticle::Pt<mcparticle::Px, mcparticle::Py>,
                            emmcparticle::Eta<mcparticle::Px, mcparticle::Py, mcparticle::Pz>,
                            emmcparticle::Phi<mcparticle::Px, mcparticle::Py>,
                            emmcparticle::P<mcparticle::Px, mcparticle::Py, mcparticle::Pz>,
                            emmcparticle::Y<mcparticle::Pz, mcparticle::E>,
                            mcparticle::ProducedByGenerator<mcparticle::Flags>,
                            mcparticle::FromBackgroundEvent<mcparticle::Flags>,
                            mcparticle::IsPhysicalPrimary<mcparticle::Flags>,
                            mcparticle::GetGenStatusCode<mcparticle::Flags, mcparticle::StatusCode>,
                            mcparticle::GetHepMCStatusCode<mcparticle::Flags, mcparticle::StatusCode>,
                            mcparticle::GetProcess<mcparticle::Flags, mcparticle::StatusCode>);

using EMMCParticles = EMMCParticles_001;
using EMMCParticle = EMMCParticles::iterator;

namespace emmcgenvectormeson
{
DECLARE_SOA_INDEX_COLUMN(EMMCEvent, emmcevent);
DECLARE_SOA_COLUMN(DownScalingFactor, dsf, float); //! down scaling factor to store this mc particle in reduced AO2D.root
} // namespace emmcgenvectormeson

DECLARE_SOA_TABLE_FULL(EMMCGenVectorMesons, "EMMCGenVectorMesons", "AOD", "EMMCGENVM", //!  generated omega, phi information
                       o2::soa::Index<>, emmcgenvectormeson::EMMCEventId,
                       mcparticle::PdgCode, mcparticle::Flags,
                       mcparticle::Px, mcparticle::Py, mcparticle::Pz, mcparticle::E,
                       emmcgenvectormeson::DownScalingFactor,

                       // dynamic column
                       emmcparticle::Pt<mcparticle::Px, mcparticle::Py>,
                       emmcparticle::Eta<mcparticle::Px, mcparticle::Py, mcparticle::Pz>,
                       emmcparticle::Phi<mcparticle::Px, mcparticle::Py>,

                       emmcparticle::P<mcparticle::Px, mcparticle::Py, mcparticle::Pz>,
                       emmcparticle::Y<mcparticle::Pz, mcparticle::E>,
                       mcparticle::ProducedByGenerator<mcparticle::Flags>,
                       mcparticle::FromBackgroundEvent<mcparticle::Flags>,
                       mcparticle::IsPhysicalPrimary<mcparticle::Flags>);

using EMMCGenVectorMeson = EMMCGenVectorMesons::iterator;

namespace smearedtrack
{
DECLARE_SOA_COLUMN(PtSmeared, ptSmeared, float);
DECLARE_SOA_COLUMN(EtaSmeared, etaSmeared, float);
DECLARE_SOA_COLUMN(PhiSmeared, phiSmeared, float);
DECLARE_SOA_COLUMN(Efficiency, efficiency, float);
DECLARE_SOA_COLUMN(DCA, dca, float);

DECLARE_SOA_COLUMN(PtSmearedSAMuon, ptSmeared_sa_muon, float);
DECLARE_SOA_COLUMN(EtaSmearedSAMuon, etaSmeared_sa_muon, float);
DECLARE_SOA_COLUMN(PhiSmearedSAMuon, phiSmeared_sa_muon, float);
DECLARE_SOA_COLUMN(EfficiencySAMuon, efficiency_sa_muon, float);
DECLARE_SOA_COLUMN(DCASAMuon, dca_sa_muon, float);
DECLARE_SOA_COLUMN(PtSmearedGLMuon, ptSmeared_gl_muon, float);
DECLARE_SOA_COLUMN(EtaSmearedGLMuon, etaSmeared_gl_muon, float);
DECLARE_SOA_COLUMN(PhiSmearedGLMuon, phiSmeared_gl_muon, float);
DECLARE_SOA_COLUMN(EfficiencyGLMuon, efficiency_gl_muon, float);
DECLARE_SOA_COLUMN(DCAGLMuon, dca_gl_muon, float);
} // namespace smearedtrack

DECLARE_SOA_TABLE(SmearedElectrons, "AOD", "SMEAREDEL", // usage Join<aod::EMMCParitlces, aod::SmearedElectrons>
                  smearedtrack::PtSmeared, smearedtrack::EtaSmeared, smearedtrack::PhiSmeared, smearedtrack::Efficiency, smearedtrack::DCA);
using SmearedElectron = SmearedElectrons::iterator;

DECLARE_SOA_TABLE(SmearedMuons, "AOD", "SMEAREDMU", // usage Join<aod::EMMCParitlces, aod::SmearedSAMuons>
                  smearedtrack::PtSmearedSAMuon, smearedtrack::EtaSmearedSAMuon, smearedtrack::PhiSmearedSAMuon, smearedtrack::EfficiencySAMuon, smearedtrack::DCASAMuon,
                  smearedtrack::PtSmearedGLMuon, smearedtrack::EtaSmearedGLMuon, smearedtrack::PhiSmearedGLMuon, smearedtrack::EfficiencyGLMuon, smearedtrack::DCAGLMuon);
using SmearedMuon = SmearedMuons::iterator;

namespace emprimaryelectronmclabel
{
DECLARE_SOA_INDEX_COLUMN(EMMCParticle, emmcparticle); //!
DECLARE_SOA_COLUMN(McMask, mcMask, uint16_t);
} // namespace emprimaryelectronmclabel

// NOTE: MC labels. This table has one entry for each reconstructed track (joinable with EMPrimaryElectrons table)
DECLARE_SOA_TABLE(EMPrimaryElectronMCLabels, "AOD", "EMPRMELMCLABEL", //!
                  emprimaryelectronmclabel::EMMCParticleId, emprimaryelectronmclabel::McMask);
using EMPrimaryElectronMCLabel = EMPrimaryElectronMCLabels::iterator;

namespace emprimarymuonmclabel
{
DECLARE_SOA_INDEX_COLUMN(EMMCParticle, emmcparticle); //!
DECLARE_SOA_COLUMN(McMask, mcMask, uint16_t);
} // namespace emprimarymuonmclabel

// NOTE: MC labels. This table has one entry for each reconstructed track (joinable with EMPrimaryMuons table)
DECLARE_SOA_TABLE(EMPrimaryMuonMCLabels, "AOD", "EMPRMMUMCLABEL", //!
                  emprimarymuonmclabel::EMMCParticleId, emprimarymuonmclabel::McMask);
using EMPrimaryMuonMCLabel = EMPrimaryMuonMCLabels::iterator;

namespace emmftmclabel
{
// DECLARE_SOA_INDEX_COLUMN_FULL(EMMCParticle, emmcparticle);
DECLARE_SOA_INDEX_COLUMN_FULL(EMMFTMCParticle, emmftmcparticle, int, EMMCParticles, "_MFT");
DECLARE_SOA_COLUMN(McMask, mcMask, uint16_t);
} // namespace emmftmclabel

// NOTE: MC labels. This table has one entry for each reconstructed track (joinable with EMPrimaryMuons table)
DECLARE_SOA_TABLE(EMMFTMCLabels, "AOD", "EMMFTMCLABEL", //!
                  emmftmclabel::EMMFTMCParticleId, emmftmclabel::McMask);
using EMMFTMCLabel = EMMFTMCLabels::iterator;

namespace emprimaryelectron
{
DECLARE_SOA_INDEX_COLUMN(EMEvent, emevent);        //!
DECLARE_SOA_COLUMN(CollisionId, collisionId, int); //!
DECLARE_SOA_COLUMN(TrackId, trackId, int);         //!
DECLARE_SOA_SELF_ARRAY_INDEX_COLUMN(AmbiguousElectrons, ambiguousElectrons);
DECLARE_SOA_COLUMN(IsAssociatedToMPC, isAssociatedToMPC, bool); //! is associated to most probable collision
DECLARE_SOA_COLUMN(IsAmbiguous, isAmbiguous, bool);             //! is ambiguous
DECLARE_SOA_COLUMN(Sign, sign, int8_t);                         //!
DECLARE_SOA_COLUMN(PrefilterBit, pfb, uint8_t);                 //!
DECLARE_SOA_COLUMN(PrefilterBitDerived, pfbderived, uint16_t);  //!
DECLARE_SOA_COLUMN(ProbElBDT, probElBDT, float);                //!

DECLARE_SOA_COLUMN(ITSNSigmaEl, itsNSigmaEl, float); //!
DECLARE_SOA_COLUMN(ITSNSigmaMu, itsNSigmaMu, float); //!
DECLARE_SOA_COLUMN(ITSNSigmaPi, itsNSigmaPi, float); //!
DECLARE_SOA_COLUMN(ITSNSigmaKa, itsNSigmaKa, float); //!
DECLARE_SOA_COLUMN(ITSNSigmaPr, itsNSigmaPr, float); //!

DECLARE_SOA_DYNAMIC_COLUMN(Signed1Pt, signed1Pt, [](float pt, int8_t sign) -> float { return sign * 1. / pt; });
DECLARE_SOA_DYNAMIC_COLUMN(P, p, [](float pt, float eta) -> float { return pt * std::cosh(eta); });
DECLARE_SOA_DYNAMIC_COLUMN(Px, px, [](float pt, float phi) -> float { return pt * std::cos(phi); });
DECLARE_SOA_DYNAMIC_COLUMN(Py, py, [](float pt, float phi) -> float { return pt * std::sin(phi); });
DECLARE_SOA_DYNAMIC_COLUMN(Pz, pz, [](float pt, float eta) -> float { return pt * std::sinh(eta); });
DECLARE_SOA_DYNAMIC_COLUMN(Theta, theta, [](float tgl) -> float { return o2::constants::math::PIHalf - std::atan(tgl); });
DECLARE_SOA_DYNAMIC_COLUMN(Tgl, tgl, [](float eta) -> float { return std::tan(o2::constants::math::PIHalf - 2 * std::atan(std::exp(-eta))); });
DECLARE_SOA_DYNAMIC_COLUMN(MeanClusterSizeITS, meanClusterSizeITS, [](uint32_t itsClusterSizes) -> float {
  int total_cluster_size = 0, nl = 0;
  for (unsigned int layer = 0; layer < 7; layer++) {
    int cluster_size_per_layer = (itsClusterSizes >> (layer * 4)) & 0xf;
    if (cluster_size_per_layer > 0) {
      nl++;
    }
    total_cluster_size += cluster_size_per_layer;
  }
  if (nl > 0) {
    return static_cast<float>(total_cluster_size) / static_cast<float>(nl);
  } else {
    return 0;
  }
});
DECLARE_SOA_DYNAMIC_COLUMN(MeanClusterSizeITSib, meanClusterSizeITSib, [](uint32_t itsClusterSizes) -> float {
  int total_cluster_size = 0, nl = 0;
  for (unsigned int layer = 0; layer < 3; layer++) {
    int cluster_size_per_layer = (itsClusterSizes >> (layer * 4)) & 0xf;
    if (cluster_size_per_layer > 0) {
      nl++;
    }
    total_cluster_size += cluster_size_per_layer;
  }
  if (nl > 0) {
    return static_cast<float>(total_cluster_size) / static_cast<float>(nl);
  } else {
    return 0;
  }
});
DECLARE_SOA_DYNAMIC_COLUMN(MeanClusterSizeITSob, meanClusterSizeITSob, [](uint32_t itsClusterSizes) -> float {
  int total_cluster_size = 0, nl = 0;
  for (unsigned int layer = 3; layer < 7; layer++) {
    int cluster_size_per_layer = (itsClusterSizes >> (layer * 4)) & 0xf;
    if (cluster_size_per_layer > 0) {
      nl++;
    }
    total_cluster_size += cluster_size_per_layer;
  }
  if (nl > 0) {
    return static_cast<float>(total_cluster_size) / static_cast<float>(nl);
  } else {
    return 0;
  }
});
} // namespace emprimaryelectron

DECLARE_SOA_TABLE_VERSIONED(EMPrimaryElectrons_001, "AOD", "EMPRIMARYEL", 1, //!
                            o2::soa::Index<>, emprimaryelectron::CollisionId,
                            emprimaryelectron::TrackId, emprimaryelectron::Sign,
                            track::Pt, track::Eta, track::Phi, track::DcaXY, track::DcaZ,
                            track::TPCNClsFindable, track::TPCNClsFindableMinusFound, track::TPCNClsFindableMinusCrossedRows, track::TPCNClsShared,
                            track::TPCChi2NCl, track::TPCInnerParam,
                            track::TPCSignal, pidtpc::TPCNSigmaEl, pidtpc::TPCNSigmaMu, pidtpc::TPCNSigmaPi, pidtpc::TPCNSigmaKa, pidtpc::TPCNSigmaPr,
                            pidtofbeta::Beta, pidtof::TOFNSigmaEl, pidtof::TOFNSigmaMu, pidtof::TOFNSigmaPi, pidtof::TOFNSigmaKa, pidtof::TOFNSigmaPr,
                            track::ITSClusterSizes, track::ITSChi2NCl, track::TOFChi2, track::DetectorMap,
                            track::X, track::Alpha, track::Y, track::Z, track::Snp, track::Tgl, emprimaryelectron::IsAssociatedToMPC,

                            // dynamic column
                            track::TPCNClsFound<track::TPCNClsFindable, track::TPCNClsFindableMinusFound>,
                            track::TPCNClsCrossedRows<track::TPCNClsFindable, track::TPCNClsFindableMinusCrossedRows>,
                            track::TPCCrossedRowsOverFindableCls<track::TPCNClsFindable, track::TPCNClsFindableMinusCrossedRows>,
                            track::TPCFoundOverFindableCls<track::TPCNClsFindable, track::TPCNClsFindableMinusFound>,
                            track::TPCFractionSharedCls<track::TPCNClsShared, track::TPCNClsFindable, track::TPCNClsFindableMinusFound>,
                            track::v001::ITSClusterMap<track::ITSClusterSizes>, track::v001::ITSNCls<track::ITSClusterSizes>, track::v001::ITSNClsInnerBarrel<track::ITSClusterSizes>,
                            track::HasITS<track::DetectorMap>, track::HasTPC<track::DetectorMap>, track::HasTRD<track::DetectorMap>, track::HasTOF<track::DetectorMap>,
                            emprimaryelectron::Signed1Pt<track::Pt, emprimaryelectron::Sign>,
                            emprimaryelectron::P<track::Pt, track::Eta>,
                            emprimaryelectron::Px<track::Pt, track::Phi>,
                            emprimaryelectron::Py<track::Pt, track::Phi>,
                            emprimaryelectron::Pz<track::Pt, track::Eta>,
                            emprimaryelectron::Theta<track::Tgl>,
                            emprimaryelectron::MeanClusterSizeITS<track::ITSClusterSizes>,
                            emprimaryelectron::MeanClusterSizeITSib<track::ITSClusterSizes>,
                            emprimaryelectron::MeanClusterSizeITSob<track::ITSClusterSizes>);

DECLARE_SOA_TABLE_VERSIONED(EMPrimaryElectrons_002, "AOD", "EMPRIMARYEL", 2, //!
                            o2::soa::Index<>, emprimaryelectron::CollisionId,
                            emprimaryelectron::TrackId, emprimaryelectron::Sign,
                            track::Pt, track::Eta, track::Phi, track::DcaXY, track::DcaZ,
                            track::TPCNClsFindable, track::TPCNClsFindableMinusFound, track::TPCNClsFindableMinusCrossedRows, track::TPCNClsShared,
                            track::TPCChi2NCl, track::TPCInnerParam,
                            track::TPCSignal, pidtpc::TPCNSigmaEl, pidtpc::TPCNSigmaMu, pidtpc::TPCNSigmaPi, pidtpc::TPCNSigmaKa, pidtpc::TPCNSigmaPr,
                            pidtofbeta::Beta, pidtof::TOFNSigmaEl, pidtof::TOFNSigmaMu, pidtof::TOFNSigmaPi, pidtof::TOFNSigmaKa, pidtof::TOFNSigmaPr,
                            track::ITSClusterSizes, emprimaryelectron::ITSNSigmaEl, emprimaryelectron::ITSNSigmaMu, emprimaryelectron::ITSNSigmaPi, emprimaryelectron::ITSNSigmaKa, emprimaryelectron::ITSNSigmaPr,
                            track::ITSChi2NCl, track::TOFChi2, track::DetectorMap,
                            track::X, track::Alpha, track::Y, track::Z, track::Snp, track::Tgl, emprimaryelectron::IsAssociatedToMPC,

                            // dynamic column
                            track::TPCNClsFound<track::TPCNClsFindable, track::TPCNClsFindableMinusFound>,
                            track::TPCNClsCrossedRows<track::TPCNClsFindable, track::TPCNClsFindableMinusCrossedRows>,
                            track::TPCCrossedRowsOverFindableCls<track::TPCNClsFindable, track::TPCNClsFindableMinusCrossedRows>,
                            track::TPCFoundOverFindableCls<track::TPCNClsFindable, track::TPCNClsFindableMinusFound>,
                            track::TPCFractionSharedCls<track::TPCNClsShared, track::TPCNClsFindable, track::TPCNClsFindableMinusFound>,
                            track::v001::ITSClusterMap<track::ITSClusterSizes>, track::v001::ITSNCls<track::ITSClusterSizes>, track::v001::ITSNClsInnerBarrel<track::ITSClusterSizes>,
                            track::HasITS<track::DetectorMap>, track::HasTPC<track::DetectorMap>, track::HasTRD<track::DetectorMap>, track::HasTOF<track::DetectorMap>,
                            emprimaryelectron::Signed1Pt<track::Pt, emprimaryelectron::Sign>,
                            emprimaryelectron::P<track::Pt, track::Eta>,
                            emprimaryelectron::Px<track::Pt, track::Phi>,
                            emprimaryelectron::Py<track::Pt, track::Phi>,
                            emprimaryelectron::Pz<track::Pt, track::Eta>,
                            emprimaryelectron::Theta<track::Tgl>,
                            emprimaryelectron::MeanClusterSizeITS<track::ITSClusterSizes>,
                            emprimaryelectron::MeanClusterSizeITSib<track::ITSClusterSizes>,
                            emprimaryelectron::MeanClusterSizeITSob<track::ITSClusterSizes>);

DECLARE_SOA_TABLE_VERSIONED(EMPrimaryElectrons_003, "AOD", "EMPRIMARYEL", 3, //!
                            o2::soa::Index<>, emprimaryelectron::CollisionId,
                            emprimaryelectron::TrackId, emprimaryelectron::Sign,
                            track::Pt, track::Eta, track::Phi, track::DcaXY, track::DcaZ,
                            track::TPCNClsFindable, track::TPCNClsFindableMinusFound, track::TPCNClsFindableMinusCrossedRows, track::TPCNClsShared,
                            track::TPCChi2NCl, track::TPCInnerParam,
                            track::TPCSignal, pidtpc::TPCNSigmaEl, /*pidtpc::TPCNSigmaMu,*/ pidtpc::TPCNSigmaPi, pidtpc::TPCNSigmaKa, pidtpc::TPCNSigmaPr,
                            pidtofbeta::Beta, pidtof::TOFNSigmaEl, /*pidtof::TOFNSigmaMu,*/ pidtof::TOFNSigmaPi, pidtof::TOFNSigmaKa, pidtof::TOFNSigmaPr,
                            track::ITSClusterSizes,
                            // emprimaryelectron::ITSNSigmaEl, emprimaryelectron::ITSNSigmaMu, emprimaryelectron::ITSNSigmaPi, emprimaryelectron::ITSNSigmaKa, emprimaryelectron::ITSNSigmaPr,
                            track::ITSChi2NCl, track::TOFChi2, track::DetectorMap,
                            track::X, track::Alpha, track::Y, track::Z, track::Snp, track::Tgl, emprimaryelectron::IsAssociatedToMPC,
                            mcpidtpc::DeDxTunedMc,

                            // dynamic column
                            track::TPCNClsFound<track::TPCNClsFindable, track::TPCNClsFindableMinusFound>,
                            track::TPCNClsCrossedRows<track::TPCNClsFindable, track::TPCNClsFindableMinusCrossedRows>,
                            track::TPCCrossedRowsOverFindableCls<track::TPCNClsFindable, track::TPCNClsFindableMinusCrossedRows>,
                            track::TPCFoundOverFindableCls<track::TPCNClsFindable, track::TPCNClsFindableMinusFound>,
                            track::TPCFractionSharedCls<track::TPCNClsShared, track::TPCNClsFindable, track::TPCNClsFindableMinusFound>,
                            track::v001::ITSClusterMap<track::ITSClusterSizes>, track::v001::ITSNCls<track::ITSClusterSizes>, track::v001::ITSNClsInnerBarrel<track::ITSClusterSizes>,
                            track::HasITS<track::DetectorMap>, track::HasTPC<track::DetectorMap>, track::HasTRD<track::DetectorMap>, track::HasTOF<track::DetectorMap>,
                            emprimaryelectron::Signed1Pt<track::Pt, emprimaryelectron::Sign>,
                            emprimaryelectron::P<track::Pt, track::Eta>,
                            emprimaryelectron::Px<track::Pt, track::Phi>,
                            emprimaryelectron::Py<track::Pt, track::Phi>,
                            emprimaryelectron::Pz<track::Pt, track::Eta>,
                            emprimaryelectron::Theta<track::Tgl>,
                            emprimaryelectron::MeanClusterSizeITS<track::ITSClusterSizes>,
                            emprimaryelectron::MeanClusterSizeITSib<track::ITSClusterSizes>,
                            emprimaryelectron::MeanClusterSizeITSob<track::ITSClusterSizes>);

DECLARE_SOA_TABLE_VERSIONED(EMPrimaryElectrons_004, "AOD", "EMPRIMARYEL", 4, //!
                            o2::soa::Index<>, emprimaryelectron::CollisionId,
                            emprimaryelectron::TrackId, emprimaryelectron::Sign,
                            track::Pt, track::Eta, track::Phi,
                            track::DcaXY, track::DcaZ, aod::track::CYY, aod::track::CZY, aod::track::CZZ,
                            track::TPCNClsFindable, track::TPCNClsFindableMinusFound, track::TPCNClsFindableMinusCrossedRows, track::TPCNClsShared,
                            track::TPCChi2NCl, track::TPCInnerParam,
                            track::TPCSignal, pidtpc::TPCNSigmaEl, pidtpc::TPCNSigmaPi, pidtpc::TPCNSigmaKa, pidtpc::TPCNSigmaPr,
                            pidtofbeta::Beta, pidtof::TOFNSigmaEl,                                         /*pidtof::TOFNSigmaPi, pidtof::TOFNSigmaKa, pidtof::TOFNSigmaPr,*/
                            track::ITSClusterSizes, track::ITSChi2NCl, track::TOFChi2, track::DetectorMap, /*track::Tgl,*/
                            emprimaryelectron::IsAssociatedToMPC, emprimaryelectron::IsAmbiguous, emprimaryelectron::ProbElBDT,
                            mcpidtpc::DeDxTunedMc,

                            // dynamic column
                            track::TPCNClsFound<track::TPCNClsFindable, track::TPCNClsFindableMinusFound>,
                            track::TPCNClsCrossedRows<track::TPCNClsFindable, track::TPCNClsFindableMinusCrossedRows>,
                            track::TPCCrossedRowsOverFindableCls<track::TPCNClsFindable, track::TPCNClsFindableMinusCrossedRows>,
                            track::TPCFoundOverFindableCls<track::TPCNClsFindable, track::TPCNClsFindableMinusFound>,
                            track::TPCFractionSharedCls<track::TPCNClsShared, track::TPCNClsFindable, track::TPCNClsFindableMinusFound>,
                            track::v001::ITSClusterMap<track::ITSClusterSizes>, track::v001::ITSNCls<track::ITSClusterSizes>, track::v001::ITSNClsInnerBarrel<track::ITSClusterSizes>,
                            track::HasITS<track::DetectorMap>, track::HasTPC<track::DetectorMap>, track::HasTRD<track::DetectorMap>, track::HasTOF<track::DetectorMap>,

                            emprimaryelectron::Signed1Pt<track::Pt, emprimaryelectron::Sign>,
                            emprimaryelectron::P<track::Pt, track::Eta>,
                            emprimaryelectron::Px<track::Pt, track::Phi>,
                            emprimaryelectron::Py<track::Pt, track::Phi>,
                            emprimaryelectron::Pz<track::Pt, track::Eta>,
                            emprimaryelectron::Tgl<track::Eta>,
                            emprimaryelectron::MeanClusterSizeITS<track::ITSClusterSizes>,
                            emprimaryelectron::MeanClusterSizeITSib<track::ITSClusterSizes>,
                            emprimaryelectron::MeanClusterSizeITSob<track::ITSClusterSizes>);

DECLARE_SOA_TABLE_VERSIONED(EMPrimaryElectrons_005, "AOD", "EMPRIMARYEL", 5, //!
                            o2::soa::Index<>, emprimaryelectron::CollisionId,
                            emprimaryelectron::TrackId, emprimaryelectron::Sign,
                            track::Pt, track::Eta, track::Phi,
                            track::DcaXY, track::DcaZ, aod::track::CYY, aod::track::CZY, aod::track::CZZ,
                            track::TPCNClsFindable, track::TPCNClsFindableMinusFound, track::TPCNClsFindableMinusPID, track::TPCNClsFindableMinusCrossedRows, track::TPCNClsShared,
                            track::TPCChi2NCl, track::TPCInnerParam,
                            track::TPCSignal, pidtpc::TPCNSigmaEl, pidtpc::TPCNSigmaPi, pidtpc::TPCNSigmaKa, pidtpc::TPCNSigmaPr,
                            pidtofbeta::Beta, pidtof::TOFNSigmaEl,                                         /*pidtof::TOFNSigmaPi, pidtof::TOFNSigmaKa, pidtof::TOFNSigmaPr,*/
                            track::ITSClusterSizes, track::ITSChi2NCl, track::TOFChi2, track::DetectorMap, /*track::Tgl,*/
                            emprimaryelectron::IsAssociatedToMPC, emprimaryelectron::IsAmbiguous, emprimaryelectron::ProbElBDT,
                            mcpidtpc::DeDxTunedMc,

                            // dynamic column
                            track::TPCNClsFound<track::TPCNClsFindable, track::TPCNClsFindableMinusFound>,
                            track::TPCNClsPID<track::TPCNClsFindable, track::TPCNClsFindableMinusPID>,
                            track::TPCNClsCrossedRows<track::TPCNClsFindable, track::TPCNClsFindableMinusCrossedRows>,
                            track::TPCCrossedRowsOverFindableCls<track::TPCNClsFindable, track::TPCNClsFindableMinusCrossedRows>,
                            track::TPCFoundOverFindableCls<track::TPCNClsFindable, track::TPCNClsFindableMinusFound>,
                            track::TPCFractionSharedCls<track::TPCNClsShared, track::TPCNClsFindable, track::TPCNClsFindableMinusFound>,
                            track::v001::ITSClusterMap<track::ITSClusterSizes>, track::v001::ITSNCls<track::ITSClusterSizes>, track::v001::ITSNClsInnerBarrel<track::ITSClusterSizes>,
                            track::HasITS<track::DetectorMap>, track::HasTPC<track::DetectorMap>, track::HasTRD<track::DetectorMap>, track::HasTOF<track::DetectorMap>,

                            emprimaryelectron::Signed1Pt<track::Pt, emprimaryelectron::Sign>,
                            emprimaryelectron::P<track::Pt, track::Eta>,
                            emprimaryelectron::Px<track::Pt, track::Phi>,
                            emprimaryelectron::Py<track::Pt, track::Phi>,
                            emprimaryelectron::Pz<track::Pt, track::Eta>,
                            emprimaryelectron::Tgl<track::Eta>,
                            emprimaryelectron::MeanClusterSizeITS<track::ITSClusterSizes>,
                            emprimaryelectron::MeanClusterSizeITSib<track::ITSClusterSizes>,
                            emprimaryelectron::MeanClusterSizeITSob<track::ITSClusterSizes>);

using EMPrimaryElectrons = EMPrimaryElectrons_005;
// iterators
using EMPrimaryElectron = EMPrimaryElectrons::iterator;

DECLARE_SOA_TABLE(EMPrimaryElectronsCov_000, "AOD", "EMPRIMARYELCOV", //!
                  aod::track::CYY,
                  aod::track::CZY,
                  aod::track::CZZ,
                  aod::track::CSnpY,
                  aod::track::CSnpZ,
                  aod::track::CSnpSnp,
                  aod::track::CTglY,
                  aod::track::CTglZ,
                  aod::track::CTglSnp,
                  aod::track::CTglTgl,
                  aod::track::C1PtY,
                  aod::track::C1PtZ,
                  aod::track::C1PtSnp,
                  aod::track::C1PtTgl,
                  aod::track::C1Pt21Pt2, o2::soa::Marker<1>);

DECLARE_SOA_TABLE_VERSIONED(EMPrimaryElectronsCov_001, "AOD", "EMPRIMARYELCOV", 1, //!
                            aod::track::X,
                            aod::track::Alpha,
                            aod::track::Y,
                            aod::track::Z,
                            aod::track::Snp,
                            aod::track::CSnpY,
                            aod::track::CSnpZ,
                            aod::track::CSnpSnp,
                            aod::track::CTglY,
                            aod::track::CTglZ,
                            aod::track::CTglSnp,
                            aod::track::CTglTgl,
                            aod::track::C1PtY,
                            aod::track::C1PtZ,
                            aod::track::C1PtSnp,
                            aod::track::C1PtTgl,
                            aod::track::C1Pt21Pt2); // CYY, CZY, CZZ, Tgl are in the main electron table.

using EMPrimaryElectronsCov = EMPrimaryElectronsCov_001;
// iterators
using EMPrimaryElectronCov = EMPrimaryElectronsCov::iterator;

DECLARE_SOA_TABLE_VERSIONED(EMPrimaryElectronsDeDxMC_000, "AOD", "EMPRMELDEDXMC", 0, mcpidtpc::DeDxTunedMc, o2::soa::Marker<1>);
using EMPrimaryElectronsDeDxMC = EMPrimaryElectronsDeDxMC_000;
// iterators
using EMPrimaryElectronDeDxMC = EMPrimaryElectronsDeDxMC::iterator;

DECLARE_SOA_TABLE(EMPrimaryElectronEMEventIds, "AOD", "PRMELMEVENTID", emprimaryelectron::EMEventId); // To be joined with EMPrimaryElectrons table at analysis level.
// iterators
using EMPrimaryElectronEMEventId = EMPrimaryElectronEMEventIds::iterator;

DECLARE_SOA_TABLE(EMPrimaryElectronsPrefilterBit, "AOD", "PRMELPFB", emprimaryelectron::PrefilterBit); // To be joined with EMPrimaryElectrons table at analysis level.
// iterators
using EMPrimaryElectronPrefilterBit = EMPrimaryElectronsPrefilterBit::iterator;

DECLARE_SOA_TABLE(EMAmbiguousElectronSelfIds, "AOD", "EMAMBELSELFID", emprimaryelectron::AmbiguousElectronsIds); // To be joined with EMPrimaryElectrons table at analysis level.
// iterators
using EMAmbiguousElectronSelfId = EMAmbiguousElectronSelfIds::iterator;

DECLARE_SOA_TABLE(EMPrimaryElectronsPrefilterBitDerived, "AOD", "PRMELPFBPI0", emprimaryelectron::PrefilterBitDerived); // To be joined with EMPrimaryElectrons table at analysis level.
// iterators
using EMPrimaryElectronPrefilterBitDerived = EMPrimaryElectronsPrefilterBitDerived::iterator;

namespace emprimarymuon
{
DECLARE_SOA_INDEX_COLUMN(EMEvent, emevent);                                          //!
DECLARE_SOA_COLUMN(CollisionId, collisionId, int);                                   //!
DECLARE_SOA_COLUMN(FwdTrackId, fwdtrackId, int);                                     //!
DECLARE_SOA_COLUMN(MFTTrackId, mfttrackId, int);                                     //!
DECLARE_SOA_COLUMN(MCHTrackId, mchtrackId, int);                                     //!
DECLARE_SOA_SELF_ARRAY_INDEX_COLUMN(GlobalMuonsWithSameMCHMID, globalMuonsWithSameMCHMID); //! self indices to global muons that have the same MCHTrackId
DECLARE_SOA_SELF_ARRAY_INDEX_COLUMN(GlobalMuonsWithSameMFT, globalMuonsWithSameMFT); //! self indices to global muons that have the same MFTTrackId
DECLARE_SOA_SELF_ARRAY_INDEX_COLUMN(AmbiguousMuons, ambiguousMuons);
DECLARE_SOA_COLUMN(CXXatDCA, cXXatDCA, float);                         //! DCAx resolution squared at DCA
DECLARE_SOA_COLUMN(CYYatDCA, cYYatDCA, float);                         //! DCAy resolution squared at DCA
DECLARE_SOA_COLUMN(CXYatDCA, cXYatDCA, float);                         //! correlation term of DCAx,y resolution at DCA
DECLARE_SOA_COLUMN(PtMatchedMCHMID, ptMatchedMCHMID, float);           //! pt of MCH-MID track in MFT-MCH-MID track at PV
DECLARE_SOA_COLUMN(EtaMatchedMCHMID, etaMatchedMCHMID, float);         //! eta of MCH-MID track in MFT-MCH-MID track at PV
DECLARE_SOA_COLUMN(PhiMatchedMCHMID, phiMatchedMCHMID, float);         //! phi of MCH-MID track in MFT-MCH-MID track at PV
DECLARE_SOA_COLUMN(EtaMatchedMCHMIDatMP, etaMatchedMCHMIDatMP, float); //! eta of MCH-MID track in MFT-MCH-MID track at matching plane
DECLARE_SOA_COLUMN(PhiMatchedMCHMIDatMP, phiMatchedMCHMIDatMP, float); //! phi of MCH-MID track in MFT-MCH-MID track at matching plane
DECLARE_SOA_COLUMN(EtaMatchedMFTatMP, etaMatchedMFTatMP, float);       //! eta of MFT track in MFT-MCH-MID track at matching plane
DECLARE_SOA_COLUMN(PhiMatchedMFTatMP, phiMatchedMFTatMP, float);       //! phi of MFT track in MFT-MCH-MID track at matching plane
DECLARE_SOA_COLUMN(IsAssociatedToMPC, isAssociatedToMPC, bool);        //! is associated to most probable collision
DECLARE_SOA_COLUMN(IsAmbiguous, isAmbiguous, bool);                    //! is ambiguous
DECLARE_SOA_COLUMN(IsCorrectMatchMFTMCH, isCorrectMatchMFTMCH, bool);  //! is correct match between MFT and MCH, only for MC
DECLARE_SOA_COLUMN(Sign, sign, int8_t);                                //!
DECLARE_SOA_COLUMN(Chi2MFT, chi2MFT, float);                           //! chi2 of MFT standalone track
DECLARE_SOA_DYNAMIC_COLUMN(Tgl, tgl, [](float eta) -> float { return std::tan(o2::constants::math::PIHalf - 2 * std::atan(std::exp(-eta))); });
DECLARE_SOA_DYNAMIC_COLUMN(Signed1Pt, signed1Pt, [](float pt, int8_t sign) -> float { return sign * 1. / pt; });
DECLARE_SOA_DYNAMIC_COLUMN(P, p, [](float pt, float eta) -> float { return pt * std::cosh(eta); });
DECLARE_SOA_DYNAMIC_COLUMN(Px, px, [](float pt, float phi) -> float { return pt * std::cos(phi); });
DECLARE_SOA_DYNAMIC_COLUMN(Py, py, [](float pt, float phi) -> float { return pt * std::sin(phi); });
DECLARE_SOA_DYNAMIC_COLUMN(Pz, pz, [](float pt, float eta) -> float { return pt * std::sinh(eta); });
DECLARE_SOA_DYNAMIC_COLUMN(NClustersMFT, nClustersMFT, //! Number of MFT clusters
                           [](uint64_t mftClusterSizesAndTrackFlags) -> uint8_t {
                             uint8_t nClusters = 0;
                             for (int layer = 0; layer < 10; layer++) {
                               if ((mftClusterSizesAndTrackFlags >> (layer * 6)) & 0x3F) {
                                 nClusters++;
                               }
                             }
                             return nClusters;
                           });
DECLARE_SOA_DYNAMIC_COLUMN(MFTClusterMap, mftClusterMap, //! MFT cluster map, one bit per a layer, starting from the innermost
                           [](uint64_t mftClusterSizesAndTrackFlags) -> uint16_t {
                             uint16_t clmap = 0;
                             for (unsigned int layer = 0; layer < 10; layer++) {
                               if ((mftClusterSizesAndTrackFlags >> (layer * 6)) & 0x3f) {
                                 clmap |= (1 << layer);
                               }
                             }
                             return clmap;
                           });
} // namespace emprimarymuon
DECLARE_SOA_TABLE(EMPrimaryMuons_000, "AOD", "EMPRIMARYMU", //!
                  o2::soa::Index<>, emprimarymuon::CollisionId,
                  emprimarymuon::FwdTrackId, emprimarymuon::MFTTrackId, emprimarymuon::MCHTrackId, fwdtrack::TrackType,
                  fwdtrack::Pt, fwdtrack::Eta, fwdtrack::Phi, emprimarymuon::Sign,
                  fwdtrack::FwdDcaX, fwdtrack::FwdDcaY, emprimarymuon::CXXatDCA, emprimarymuon::CYYatDCA, emprimarymuon::CXYatDCA,
                  emprimarymuon::PtMatchedMCHMID, emprimarymuon::EtaMatchedMCHMID, emprimarymuon::PhiMatchedMCHMID,
                  // fwdtrack::X, fwdtrack::Y, fwdtrack::Z, fwdtrack::Tgl,

                  fwdtrack::NClusters, fwdtrack::PDca, fwdtrack::RAtAbsorberEnd,
                  fwdtrack::Chi2, fwdtrack::Chi2MatchMCHMID, fwdtrack::Chi2MatchMCHMFT,
                  fwdtrack::MCHBitMap, fwdtrack::MIDBitMap, fwdtrack::MIDBoards,
                  fwdtrack::MFTClusterSizesAndTrackFlags, emprimarymuon::Chi2MFT, emprimarymuon::IsAssociatedToMPC, emprimarymuon::IsAmbiguous,

                  // dynamic column
                  emprimarymuon::Signed1Pt<fwdtrack::Pt, emprimarymuon::Sign>,
                  emprimarymuon::NClustersMFT<fwdtrack::MFTClusterSizesAndTrackFlags>,
                  emprimarymuon::MFTClusterMap<fwdtrack::MFTClusterSizesAndTrackFlags>,
                  emprimarymuon::P<fwdtrack::Pt, fwdtrack::Eta>,
                  emprimarymuon::Px<fwdtrack::Pt, fwdtrack::Phi>,
                  emprimarymuon::Py<fwdtrack::Pt, fwdtrack::Phi>,
                  emprimarymuon::Pz<fwdtrack::Pt, fwdtrack::Eta>);

DECLARE_SOA_TABLE_VERSIONED(EMPrimaryMuons_001, "AOD", "EMPRIMARYMU", 1, //!
                            o2::soa::Index<>, emprimarymuon::CollisionId,
                            emprimarymuon::FwdTrackId, emprimarymuon::MFTTrackId, emprimarymuon::MCHTrackId, fwdtrack::TrackType,
                            fwdtrack::Pt, fwdtrack::Eta, fwdtrack::Phi, emprimarymuon::Sign,
                            fwdtrack::FwdDcaX, fwdtrack::FwdDcaY, emprimarymuon::CXXatDCA, emprimarymuon::CYYatDCA, emprimarymuon::CXYatDCA,
                            emprimarymuon::PtMatchedMCHMID, emprimarymuon::EtaMatchedMCHMID, emprimarymuon::PhiMatchedMCHMID,
                            emprimarymuon::EtaMatchedMCHMIDatMP, emprimarymuon::PhiMatchedMCHMIDatMP,
                            emprimarymuon::EtaMatchedMFTatMP, emprimarymuon::PhiMatchedMFTatMP,

                            fwdtrack::NClusters, fwdtrack::PDca, fwdtrack::RAtAbsorberEnd,
                            fwdtrack::Chi2, fwdtrack::Chi2MatchMCHMID, fwdtrack::Chi2MatchMCHMFT,
                            fwdtrack::MCHBitMap, fwdtrack::MIDBitMap, fwdtrack::MIDBoards,
                            fwdtrack::MFTClusterSizesAndTrackFlags, emprimarymuon::Chi2MFT, emprimarymuon::IsAssociatedToMPC, emprimarymuon::IsAmbiguous,

                            // dynamic column
                            emprimarymuon::Signed1Pt<fwdtrack::Pt, emprimarymuon::Sign>,
                            emprimarymuon::NClustersMFT<fwdtrack::MFTClusterSizesAndTrackFlags>,
                            fwdtrack::IsCA<fwdtrack::MFTClusterSizesAndTrackFlags>,
                            emprimarymuon::MFTClusterMap<fwdtrack::MFTClusterSizesAndTrackFlags>,
                            emprimarymuon::P<fwdtrack::Pt, fwdtrack::Eta>,
                            emprimarymuon::Px<fwdtrack::Pt, fwdtrack::Phi>,
                            emprimarymuon::Py<fwdtrack::Pt, fwdtrack::Phi>,
                            emprimarymuon::Pz<fwdtrack::Pt, fwdtrack::Eta>);

DECLARE_SOA_TABLE_VERSIONED(EMPrimaryMuons_002, "AOD", "EMPRIMARYMU", 2, //!
                            o2::soa::Index<>, emprimarymuon::CollisionId,
                            emprimarymuon::FwdTrackId, emprimarymuon::MFTTrackId, emprimarymuon::MCHTrackId, fwdtrack::TrackType,
                            fwdtrack::Pt, fwdtrack::Eta, fwdtrack::Phi, emprimarymuon::Sign,
                            fwdtrack::FwdDcaX, fwdtrack::FwdDcaY, aod::fwdtrack::CXX, aod::fwdtrack::CYY, aod::fwdtrack::CXY,
                            emprimarymuon::PtMatchedMCHMID, emprimarymuon::EtaMatchedMCHMID, emprimarymuon::PhiMatchedMCHMID,
                            // emprimarymuon::EtaMatchedMCHMIDatMP, emprimarymuon::PhiMatchedMCHMIDatMP,
                            // emprimarymuon::EtaMatchedMFTatMP, emprimarymuon::PhiMatchedMFTatMP,

                            fwdtrack::NClusters, fwdtrack::PDca, fwdtrack::RAtAbsorberEnd,
                            fwdtrack::Chi2, fwdtrack::Chi2MatchMCHMID, fwdtrack::Chi2MatchMCHMFT,
                            fwdtrack::MCHBitMap, fwdtrack::MIDBitMap, fwdtrack::MIDBoards,
                            fwdtrack::MFTClusterSizesAndTrackFlags, emprimarymuon::Chi2MFT, emprimarymuon::IsAssociatedToMPC, emprimarymuon::IsAmbiguous,

                            // dynamic column
                            emprimarymuon::Signed1Pt<fwdtrack::Pt, emprimarymuon::Sign>,
                            emprimarymuon::Tgl<fwdtrack::Eta>,
                            emprimarymuon::NClustersMFT<fwdtrack::MFTClusterSizesAndTrackFlags>,
                            fwdtrack::IsCA<fwdtrack::MFTClusterSizesAndTrackFlags>,
                            emprimarymuon::MFTClusterMap<fwdtrack::MFTClusterSizesAndTrackFlags>,
                            emprimarymuon::P<fwdtrack::Pt, fwdtrack::Eta>,
                            emprimarymuon::Px<fwdtrack::Pt, fwdtrack::Phi>,
                            emprimarymuon::Py<fwdtrack::Pt, fwdtrack::Phi>,
                            emprimarymuon::Pz<fwdtrack::Pt, fwdtrack::Eta>);

using EMPrimaryMuons = EMPrimaryMuons_002;
// iterators
using EMPrimaryMuon = EMPrimaryMuons::iterator;

DECLARE_SOA_TABLE_VERSIONED(EMPrimaryMuonsCov_002, "AOD", "EMPRIMARYMUCOV", 2, //!
                            fwdtrack::X, fwdtrack::Y, fwdtrack::Z,             // at PV. Signed1Pt, Tgl and Phi are in EMPrimaryMuons table.
                            // aod::fwdtrack::CXX,
                            // aod::fwdtrack::CXY,
                            // aod::fwdtrack::CYY,
                            aod::fwdtrack::CPhiX,
                            aod::fwdtrack::CPhiY,
                            aod::fwdtrack::CPhiPhi,
                            aod::fwdtrack::CTglX,
                            aod::fwdtrack::CTglY,
                            aod::fwdtrack::CTglPhi,
                            aod::fwdtrack::CTglTgl,
                            aod::fwdtrack::C1PtX,
                            aod::fwdtrack::C1PtY,
                            aod::fwdtrack::C1PtPhi,
                            aod::fwdtrack::C1PtTgl,
                            aod::fwdtrack::C1Pt21Pt2);
using EMPrimaryMuonsCov = EMPrimaryMuonsCov_002;
// iterators
using EMPrimaryMuonCov = EMPrimaryMuonsCov::iterator;

DECLARE_SOA_TABLE(EMPrimaryMuonEMEventIds, "AOD", "PRMMUEMEVENTID", emprimarymuon::EMEventId); // To be joined with EMPrimaryMuons table at analysis level.
// iterators
using EMPrimaryMuonEMEventId = EMPrimaryMuonEMEventIds::iterator;

DECLARE_SOA_TABLE(EMAmbiguousMuonSelfIds, "AOD", "EMAMBMUSELFID", emprimarymuon::AmbiguousMuonsIds); // To be joined with EMPrimaryMuons table at analysis level.
// iterators
using EMAmbiguousMuonSelfId = EMAmbiguousMuonSelfIds::iterator;

DECLARE_SOA_TABLE(EMGlobalMuonSelfIds_000, "AOD", "EMGLMUSELFID", emprimarymuon::GlobalMuonsWithSameMFTIds);                                                           // To be joined with EMPrimaryMuons table at analysis level.
DECLARE_SOA_TABLE_VERSIONED(EMGlobalMuonSelfIds_001, "AOD", "EMGLMUSELFID", 1, emprimarymuon::GlobalMuonsWithSameMCHMIDIds, emprimarymuon::GlobalMuonsWithSameMFTIds); // To be joined with EMPrimaryMuons table at analysis level.
using EMGlobalMuonSelfIds = EMGlobalMuonSelfIds_001;
// iterators
using EMGlobalMuonSelfId = EMGlobalMuonSelfIds::iterator;

DECLARE_SOA_TABLE(EMPrimaryMuonsMatchMC, "AOD", "EMMUONMATCHMC", emprimarymuon::IsCorrectMatchMFTMCH); // To be joined with EMPrimaryMuons table at analysis level. only for MC.
// iterators
using EMPrimaryMuonMatchMC = EMPrimaryMuonsMatchMC::iterator;

namespace oldemprimarytrack
{
DECLARE_SOA_COLUMN(Sign, sign, int8_t);
} // namespace oldemprimarytrack

namespace emprimarytrack
{
DECLARE_SOA_INDEX_COLUMN(EMEvent, emevent);        //!
DECLARE_SOA_COLUMN(CollisionId, collisionId, int); //!
DECLARE_SOA_COLUMN(TrackId, trackId, int);         //!
DECLARE_SOA_COLUMN(TrackBit, trackBit, uint16_t);  //!
DECLARE_SOA_COLUMN(Signed1Pt, signed1Pt, float);   //! (sign of charge)/Pt in c/GeV. Use pt() and sign() instead
DECLARE_SOA_COLUMN(Eta, eta, float);               //!
DECLARE_SOA_COLUMN(Phi, phi, float);               //!
DECLARE_SOA_DYNAMIC_COLUMN(Pt, pt, [](float signed1Pt) -> float { return 1.f / std::fabs(signed1Pt); });
DECLARE_SOA_DYNAMIC_COLUMN(Sign, sign, [](float signed1Pt) -> short { return (signed1Pt > 0) ? 1 : -1; }); //! Charge: positive: 1, negative: -1
} // namespace emprimarytrack

DECLARE_SOA_TABLE_VERSIONED(EMPrimaryTracks_000, "AOD", "EMPRIMARYTRACK", 0, //! primary charged track table for 2PC
                            o2::soa::Index<>, emprimarytrack::CollisionId, emprimarytrack::TrackId, oldemprimarytrack::Sign, track::Pt, track::Eta, track::Phi, emprimarytrack::TrackBit);

DECLARE_SOA_TABLE_VERSIONED(EMPrimaryTracks_001, "AOD", "EMPRIMARYTRACK", 1,   //! primary charged track table for 2PC
                            o2::soa::Index<>, /*emprimarytrack::CollisionId,*/ /*emprimarytrack::TrackId,*/
                            emprimarytrack::Signed1Pt, emprimarytrack::Eta, emprimarytrack::Phi, emprimarytrack::TrackBit,
                            // dynamic column
                            emprimarytrack::Sign<emprimarytrack::Signed1Pt>, emprimarytrack::Pt<emprimarytrack::Signed1Pt>);

using EMPrimaryTracks = EMPrimaryTracks_001;
// iterators
using EMPrimaryTrack = EMPrimaryTracks::iterator;

DECLARE_SOA_TABLE(EMPrimaryTrackEMEventIds, "AOD", "PRMTRKEMEVENTID", emprimarytrack::EMEventId); // To be joined with EMPrimaryTracks table at analysis level.
// iterators
using EMPrimaryTrackEMEventId = EMPrimaryTrackEMEventIds::iterator;

DECLARE_SOA_TABLE(EMPrimaryTrackEMEventIdsTMP, "AOD", "PRMTRKEVIDTMP", track::CollisionId); // To be joined with EMPrimaryTracks in associateDileptonToEMEvent
// iterators
using EMPrimaryTrackEMEventIdTMP = EMPrimaryTrackEMEventIdsTMP::iterator;

namespace emthinevent
{
DECLARE_SOA_COLUMN(EP2, ep2, float); //!
} // namespace emthinevent
DECLARE_SOA_TABLE_VERSIONED(EMThinEvents_000, "AOD", "EMTHINEVENT", 0, //! Thin event information table
                            o2::soa::Index<>, bc::RunNumber, bc::GlobalBC, timestamp::Timestamp, collision::PosZ,
                            evsel::NumTracksInTimeRange, evsel::SumAmpFT0CInTimeRange, cent::CentFT0C, emthinevent::EP2);
using EMThinEvents = EMThinEvents_000;
using EMThinEvent = EMThinEvents::iterator;

DECLARE_SOA_TABLE_VERSIONED(EMThinEventNormInfos_000, "AOD", "EMTHINEVENTNORM", 0,                                                                                       //! event information for normalization
                            o2::soa::Index<>, emevsel::Selection, evsel::Rct, emevent::PosZint8, emevent::CentFT0Muint8, emevent::CentFT0Cuint8, emevent::CentNTPVuint8, /*emevent::CentNGlobaluint8,*/
                            emevsel::Sel8<emevsel::Selection>, emeventnorm::PosZ<emevent::PosZint8>, emeventnorm::CentFT0M<emevent::CentFT0Muint8>, emeventnorm::CentFT0C<emevent::CentFT0Cuint8>, emeventnorm::CentNTPV<emevent::CentNTPVuint8>, /*emeventnorm::CentNGlobal<emevent::CentNGlobaluint8>,*/ o2::soa::Marker<2>);
using EMThinEventNormInfos = EMThinEventNormInfos_000;
using EMThinEventNormInfo = EMThinEventNormInfos::iterator;

namespace emdilepton
{
DECLARE_SOA_INDEX_COLUMN(EMThinEvent, emthinevent); //!
DECLARE_SOA_COLUMN(Pt1, pt1, float);                //!
DECLARE_SOA_COLUMN(Eta1, eta1, float);              //!
DECLARE_SOA_COLUMN(Phi1, phi1, float);              //!
DECLARE_SOA_COLUMN(Sign1, sign1, short);            //!
DECLARE_SOA_COLUMN(DCA1, dca1, float);              //! DCA in sigma. Users should decide 3D or XY or Z
DECLARE_SOA_COLUMN(Pt2, pt2, float);                //!
DECLARE_SOA_COLUMN(Eta2, eta2, float);              //!
DECLARE_SOA_COLUMN(Phi2, phi2, float);              //!
DECLARE_SOA_COLUMN(DCA2, dca2, float);              //! DCA in sigma. Users should decide 3D or XY or Z
DECLARE_SOA_COLUMN(Sign2, sign2, short);            //!
DECLARE_SOA_COLUMN(Weight, weight, float);          //! possible pair weight
} // namespace emdilepton

DECLARE_SOA_TABLE_VERSIONED(EMDileptons_000, "AOD", "EMDILEPTON", 0,
                            o2::soa::Index<>, emdilepton::EMThinEventId,
                            emdilepton::Pt1, emdilepton::Eta1, emdilepton::Phi1, emdilepton::Sign1, emdilepton::DCA1,
                            emdilepton::Pt2, emdilepton::Eta2, emdilepton::Phi2, emdilepton::Sign2, emdilepton::DCA2,
                            emdilepton::Weight);
using EMDileptons = EMDileptons_000;
using EMDilepton = EMDileptons::iterator;

// Dummy data for MC
namespace emdummydata
{
DECLARE_SOA_COLUMN(A, a, float);
DECLARE_SOA_COLUMN(B, b, float);
DECLARE_SOA_COLUMN(C, c, float);
DECLARE_SOA_COLUMN(D, d, float);
DECLARE_SOA_COLUMN(E, e, float);
DECLARE_SOA_COLUMN(F, f, float);
DECLARE_SOA_COLUMN(G, g, float);
DECLARE_SOA_COLUMN(H, h, float);
DECLARE_SOA_COLUMN(I, i, float);
DECLARE_SOA_COLUMN(J, j, float);
DECLARE_SOA_COLUMN(K, k, float);
DECLARE_SOA_COLUMN(L, l, float);
DECLARE_SOA_COLUMN(M, m, float);
DECLARE_SOA_COLUMN(N, n, float);
DECLARE_SOA_COLUMN(O, o, float);
DECLARE_SOA_COLUMN(P, p, float);
DECLARE_SOA_COLUMN(Q, q, float);
DECLARE_SOA_COLUMN(R, r, float);
DECLARE_SOA_COLUMN(S, s, float);
DECLARE_SOA_COLUMN(T, t, float);
DECLARE_SOA_COLUMN(U, u, float);
DECLARE_SOA_COLUMN(V, v, float);
DECLARE_SOA_COLUMN(W, w, float);
DECLARE_SOA_COLUMN(X, x, float);
DECLARE_SOA_COLUMN(Y, y, float);
DECLARE_SOA_COLUMN(Z, z, float);
} // namespace emdummydata
DECLARE_SOA_TABLE(EMDummyDatas, "AOD", "EMDUMMYDATA",
                  o2::soa::Index<>,
                  emdummydata::A, emdummydata::B, emdummydata::C, emdummydata::D, emdummydata::E,
                  emdummydata::F, emdummydata::G, emdummydata::H, emdummydata::I, emdummydata::J,
                  emdummydata::K, emdummydata::L, emdummydata::M, emdummydata::N, emdummydata::O,
                  emdummydata::P, emdummydata::Q, emdummydata::R, emdummydata::S, emdummydata::T,
                  emdummydata::U, emdummydata::V, emdummydata::W, emdummydata::X, emdummydata::Y,
                  emdummydata::Z);

// iterators
using EMDummyData = EMDummyDatas::iterator;
} // namespace o2::aod

#endif // PWGEM_DILEPTON_DATAMODEL_DILEPTONTABLES_H_
