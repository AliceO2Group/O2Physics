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
/// This task writes a collision and track table which are further used in
/// a diffusion wake analysis
///
/// author N. Wilson

/// \file tableDiffWake.cxx
/// \brief This task writes a collision and track table which are further used in a diffusion wake analysis
/// \author Nicola Wilson <nicola.wilson@cern.ch>

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/Core/EventPlaneHelper.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Qvectors.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

// Event selection: Only events that contain track above some threshold
// 
// -------------------------------------------------------------------------------------------
// TRACK DATA
// -------------------------------------------------------------------------------------------
// BEFORE COMPRESSION                              AFTER COMPRESSION
// Name            Data Type       Size(b)         Name            Data Type       Size(b)
// ColID           int32_t         4               [same]
// Charge          short           2               [same]
// Px, Py, Pz      float           3x4             P               unsigned long   8
// DEdx            float           4               DEdx            unsigned short  2
// DCAXY           float           4               DCAXY           short           2
// DCAZ            float           4               DCAZ            short           2
// Length          float           4               Length          unsigned short  2

// OVERALL COMPRESSION 34b->22b

// -------------------------------------------------------------------------------------------
// EVENT DATA
// -------------------------------------------------------------------------------------------
// GI              int64_t         8               [same]
// RN              int             4               [same]
// Cent            float           4               [same]
// Mult            int             4               [same]
// VertexX         float           4               [same]
// VertexY         float           4               [same]
// VertexZ         float           4               [same]
// Psi2            float           4               Psi2            short           2
// Psi3            float           4               Psi3            short           2

// OVERALL COMPRESSION 40b->36b
//--------------------------------------------------------
namespace o2::aod
{
namespace testcol
{
// Event properties
// DECLARE_SOA_COLUMN(Gi, gi, int64_t);
DECLARE_SOA_COLUMN(Rn, rn, int32_t);         // run number
DECLARE_SOA_COLUMN(Cent, cent, float);       // FT0C centrality
DECLARE_SOA_COLUMN(Mult, mult, int32_t);     // TPC multiplicity
DECLARE_SOA_COLUMN(Occu, occu, int32_t);     // Occupancy ITS
DECLARE_SOA_COLUMN(Occuft0, occuft0, float); // Occupancy FT0C amplitudes
DECLARE_SOA_COLUMN(VertexX, vertexX, float);
DECLARE_SOA_COLUMN(VertexY, vertexY, float);
DECLARE_SOA_COLUMN(VertexZ, vertexZ, float);
DECLARE_SOA_COLUMN(Psi2, psi2, int16_t);
DECLARE_SOA_COLUMN(Psi3, psi3, int16_t);
} // namespace testcol

DECLARE_SOA_TABLE(TableCol, "AOD", "TABLECOL",
                  o2::soa::Index<>,
                  testcol::Rn,
                  testcol::Cent,
                  testcol::Mult,
                  testcol::Occu,
                  testcol::Occuft0,
                  testcol::VertexX,
                  testcol::VertexY,
                  testcol::VertexZ,
                  testcol::Psi2,
                  testcol::Psi3);
using Collision = TableCol::iterator;

namespace testtrack
{

// Track properties
DECLARE_SOA_INDEX_COLUMN(Collision, collision);
DECLARE_SOA_COLUMN(Charge, charge, int16_t);
DECLARE_SOA_COLUMN(P, p, uint64_t);
DECLARE_SOA_COLUMN(Dedx, dedx, uint16_t);
DECLARE_SOA_COLUMN(Dcaxy, dcaxy, int16_t);
DECLARE_SOA_COLUMN(Dcaz, dcaz, int16_t);
} // namespace testtrack

DECLARE_SOA_TABLE(TableTrack, "AOD", "TABLETRACK",
                  o2::soa::Index<>,
                  testtrack::CollisionId,
                  testtrack::Charge,
                  testtrack::P,
                  testtrack::Dedx,
                  testtrack::Dcaxy,
                  testtrack::Dcaz);
} // namespace o2::aod
//--------------------------------------------------------
using namespace o2;
using namespace o2::framework;

struct TableDiffWake {

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  Configurable<int> nBinsPt{"nBinsPt", 100, "N bins in pT histo"};
  Configurable<double> ptThresh{"ptThresh", 20.0, "pT threshold"};
  Configurable<float> centMax{"centMax", 10, "centrality"};
  Configurable<float> zVertCut{"zVertCut", 10.0, "z_vertex cut"};

  Produces<o2::aod::TableCol> testcol;
  Produces<o2::aod::TableTrack> testtrack;

  EventPlaneHelper helperEP;

  void init(InitContext const&)
  {
    const AxisSpec axisEta{30, -1.5, +1.5, "#eta"};
    const AxisSpec axispT{nBinsPt, 0, 250, "p_{T}"};

    histos.add("etaHistogram", "etaHistogram", kTH1F, {axisEta});
    histos.add("pTHistogram", "pTHistogram", kTH1F, {axispT});
  }

  using Bcs = aod::BCs;
  void process(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs, aod::TPCMults, aod::QvectorFT0Cs>::iterator const& col,
               soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection> const& tracks,
               Bcs const&)
  {
    const float maxMomentum = 173.0;

    // Event selection corresponds to sel8FullPbPb
    if (!col.sel8())
      return;
    if (col.centFT0C() > centMax)
      return; // Centrality 0 - 10 %
    if (std::abs(col.posZ()) > zVertCut)
      return; // z position < 10 cm
    if (!col.selection_bit(o2::aod::evsel::kNoCollInRofStandard))
      return;
    if (!col.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard))
      return;

    //------ Get Run number ---------------------
    auto bc = col.bc_as<Bcs>();
    int run = bc.runNumber();
    //------------ EP ---------------------------
    double ep2 = 0.0;
    double ep3 = 0.0;
    ep2 = helperEP.GetEventPlane(col.qvecFT0CRe(), col.qvecFT0CIm(), 2);
    ep3 = helperEP.GetEventPlane(col.qvecFT0CRe(), col.qvecFT0CIm(), 3);

    //------- Only events with track above some thresh ----------

    bool eventHighpT = false;
    for (auto const &track : tracks) {

      if (!track.isGlobalTrack())
        continue;
      if (track.pt() > ptThresh) {
        eventHighpT = true;
        break;
      }
    }
    if (!eventHighpT)
      return;
    //------------------------------------------------------------
    // Translate values to less memory consuming values
    int16_t substituteEp2 = static_cast<int16_t>(ep2 * 1000);
    int16_t substituteEp3 = static_cast<int16_t>(ep3 * 1000);

    testcol(col.globalIndex(),
            run,
            col.centFT0C(),
            col.multTPC(),
            col.trackOccupancyInTimeRange(),
            col.ft0cOccupancyInTimeRange(),
            col.posX(),
            col.posY(),
            col.posZ(),
            substituteEp2,
            substituteEp3);

    for (auto const &track : tracks) {

      // Track cut
      if (!track.isGlobalTrack())
        continue; // General track cuts

      if (std::abs(track.px()) > maxMomentum || std::abs(track.py()) > maxMomentum || std::abs(track.pz()) > maxMomentum)
        continue; // to avoid overflow in Substitute_p

      histos.fill(HIST("etaHistogram"), track.eta());
      histos.fill(HIST("pTHistogram"), track.pt());

      //------------ Translate values to less memory consuming values --------------------
      // Px, Py, Pz
      uint64_t substituteP = 0;
      uint8_t uppermostBit = 20;
      uint8_t lowermostBit = 0;

      int64_t particlePx = (track.px() * 6000);
      if (particlePx < 0)
        substituteP |= static_cast<uint64_t>1 << uppermostBit;
      if (particlePx < 0)
        particlePx = (-1) * particlePx;
      for (int8_t i_bit = lowermostBit; i_bit < uppermostBit; i_bit++) {
        if ((particlePx & (static_cast<int64_t>1 << i_bit)))
          substituteP |= static_cast<uint64_t>1 << i_bit;
      }

      uppermostBit = 41;
      lowermostBit = 21;
      int64_t particlePy = (track.py() * 6000);
      if (particlePy < 0)
        substituteP |= static_cast<uint64_t>1 << uppermostBit;
      if (particlePy < 0)
        particlePy = (-1) * particlePy;
      for (int8_t i_bit = lowermostBit; i_bit < uppermostBit; i_bit++) {
        if ((particlePy & (static_cast<int64_t>1 << (i_bit - lowermostBit))))
          substituteP |= static_cast<uint64_t>1 << i_bit;
      }

      uppermostBit = 62;
      lowermostBit = 42;
      int64_t particlePz = (track.pz() * 6000);
      if (particlePz < 0)
        substituteP |= static_cast<uint64_t>1 << uppermostBit;
      if (particlePz < 0)
        particlePz = (-1) * particlePz;
      for (int8_t i_bit = lowermostBit; i_bit < uppermostBit; i_bit++) {
        if ((particlePz & (static_cast<int64_t>1 << (i_bit - lowermostBit))))
          substituteP |= static_cast<uint64_t>1 << i_bit;
      }

      // dEdx
      uint16_t substituteDEDX = static_cast<uint16_t>(track.tpcSignal() * 10);

      // DCA
      int16_t substituteDCAXY = static_cast<int16_t>(track.dcaXY() * 100);
      int16_t substituteDCAZ = static_cast<int16_t>(track.dcaZ() * 100);

      //--------------- Fill track table ------------------
      testtrack(track.collisionId(),
                track.sign(),
                substituteP,
                substituteDEDX,
                substituteDCAXY,
                substituteDCAZ);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<TableDiffWake>(cfgc)};
}
