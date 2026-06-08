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

#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>
// For centrality:
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
// For TPC Mult
#include "Common/DataModel/Multiplicity.h"
// For DCA and TrackSelection
#include "Common/DataModel/TrackSelectionTables.h"
// For EP
#include "Common/Core/EventPlaneHelper.h"
#include "Common/DataModel/Qvectors.h"
// For occupancy bit
#include "Common/CCDB/EventSelectionParams.h"

// Event selection: Only events that contain track above some threshold
/*
-------------------------------------------------------------------------------------------
TRACK DATA
-------------------------------------------------------------------------------------------
BEFORE COMPRESSION                              AFTER COMPRESSION
Name            Data Type       Size(b)         Name            Data Type       Size(b)
ColID           int32_t         4               [same]
Charge          short           2               [same]
Px, Py, Pz      float           3x4             P               unsigned long   8
DEdx            float           4               DEdx            unsigned short  2
DCAXY           float           4               DCAXY           short           2
DCAZ            float           4               DCAZ            short           2
Length          float           4               Length          unsigned short  2

OVERALL COMPRESSION 34b->22b

-------------------------------------------------------------------------------------------
EVENT DATA
-------------------------------------------------------------------------------------------
GI              int64_t         8               [same]
RN              int             4               [same]
Cent            float           4               [same]
Mult            int             4               [same]
VertexX         float           4               [same]
VertexY         float           4               [same]
VertexZ         float           4               [same]
Psi2            float           4               Psi2            short           2
Psi3            float           4               Psi3            short           2

OVERALL COMPRESSION 40b->36b

*/

//--------------------------------------------------------
namespace o2::aod
{
namespace testcol
{
// Event properties
DECLARE_SOA_COLUMN(Gi, gi, int64_t);
DECLARE_SOA_COLUMN(Rn, rn, int32_t);             // run number
DECLARE_SOA_COLUMN(Cent, cent, float);       // FT0C centrality
DECLARE_SOA_COLUMN(Mult, mult, int32_t);         // TPC multiplicity
DECLARE_SOA_COLUMN(Occu, occu, int32_t);         // Occupancy ITS
DECLARE_SOA_COLUMN(Occuft0, occuft0, float); // Occupancy FT0C amplitudes
DECLARE_SOA_COLUMN(VertexX, vertexX, float);
DECLARE_SOA_COLUMN(VertexY, vertexY, float);
DECLARE_SOA_COLUMN(VertexZ, vertexZ, float);
DECLARE_SOA_COLUMN(Psi2, psi2, int16_t);
DECLARE_SOA_COLUMN(Psi3, psi3, int16_t);
} // namespace testcol
namespace testtrack
{

// Track properties
DECLARE_SOA_COLUMN(Colid, colid, int32_t); // Collision ID
DECLARE_SOA_COLUMN(Charge, charge, int16_t);
DECLARE_SOA_COLUMN(P, p, uint64_t);
DECLARE_SOA_COLUMN(Dedx, dedx, uint16_t);
DECLARE_SOA_COLUMN(Dcaxy, dcaxy, int16_t);
DECLARE_SOA_COLUMN(Dcaz, dcaz, int16_t);
} // namespace testtrack
DECLARE_SOA_TABLE(TableCol, "AOD", "TABLECOL",
                  testcol::Gi,
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
DECLARE_SOA_TABLE(TableTrack, "AOD", "TABLETRACK",
                  testtrack::Colid,
                  testtrack::Charge,
                  testtrack::P,
                  testtrack::Dedx,
                  testtrack::Dcaxy,
                  testtrack::Dcaz);
} // namespace o2::aod
//--------------------------------------------------------
using namespace o2;
using namespace o2::framework;

struct tableDiffWake {

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  Configurable<int> nBinsPt{"nBinsPt", 100, "N bins in pT histo"};
  Configurable<double> pT_thresh{"pT_thresh", 20.0, "pT threshold"};
  Configurable<float> cent_max{"cent_max", 10, "centrality"};
  Configurable<float> z_vert_cut{"z_vert_cut", 10.0, "z_vertex cut"};

  Produces<o2::aod::TableCol> testcol;
  Produces<o2::aod::TableTrack> testtrack;

  EventPlaneHelper helperEP;

  void init(InitContext const&)
  {
    const AxisSpec axisEta{30, -1.5, +1.5, "#eta"};
    const AxisSpec axispT{nBinsPt, 0, 10, "p_{T}"};

    histos.add("etaHistogram", "etaHistogram", kTH1F, {axisEta});
    histos.add("pTHistogram", "pTHistogram", kTH1F, {axispT});
  }

  using bcs = aod::BCs;
  void process(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs, aod::TPCMults, aod::QvectorFT0Cs>::iterator const& col,
               soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection> const& tracks,
               bcs const&)
  {
    // Event selection corresponds to sel8FullPbPb
    if (!col.sel8())
      return;
    if (col.centFT0C() > cent_max)
      return; // Centrality 0 - 10 %
    if (std::abs(col.posZ()) > z_vert_cut)
      return; // z position < 10 cm
    if (!col.selection_bit(o2::aod::evsel::kNoCollInRofStandard))
      return;
    if (!col.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard))
      return;

    //------ Get Run number ---------------------
    auto bc = col.bc_as<bcs>();
    int run = bc.runNumber();
    //------------ EP ---------------------------
    double ep2 = 0.0;
    double ep3 = 0.0;
    ep2 = helperEP.GetEventPlane(col.qvecFT0CRe(), col.qvecFT0CIm(), 2);
    ep3 = helperEP.GetEventPlane(col.qvecFT0CRe(), col.qvecFT0CIm(), 3);

    //------- Only events with track above some thresh ----------

    bool eventHighpT = false;
    for (auto& track : tracks) {

      if (!track.isGlobalTrack())
        continue;
      if (track.pt() > pT_thresh) {
        eventHighpT = true;
        break;
      }
    }
    if (!eventHighpT)
      return;
    //------------------------------------------------------------
    // Translate values to less memory consuming values
    Short_t Substitute_ep2 = (Short_t)(ep2 * 1000);
    Short_t Substitute_ep3 = (Short_t)(ep3 * 1000);

    testcol(col.globalIndex(),
            run,
            col.centFT0C(),
            col.multTPC(),
            col.trackOccupancyInTimeRange(),
            col.ft0cOccupancyInTimeRange(),
            col.posX(),
            col.posY(),
            col.posZ(),
            Substitute_ep2,
            Substitute_ep3);

    for (auto& track : tracks) {

      // Track cut
      if (!track.isGlobalTrack())
        continue; // General track cuts

      histos.fill(HIST("etaHistogram"), track.eta());
      histos.fill(HIST("pTHistogram"), track.pt());

      //------------ Translate values to less memory consuming values --------------------
      // Px, Py, Pz
      ULong64_t Substitute_p = 0;

      Long64_t Particle_px = (track.px() * 6000);
      if (Particle_px < 0)
        Substitute_p |= (ULong64_t)1 << 20;
      if (Particle_px < 0)
        Particle_px = (-1) * Particle_px;
      for (Int_t i_bit = 0; i_bit < 20; i_bit++) {
        if ((Particle_px & ((Long64_t)1 << i_bit)))
          Substitute_p |= (ULong64_t)1 << i_bit;
      }

      Long64_t Particle_py = (track.py() * 6000);
      if (Particle_py < 0)
        Substitute_p |= (ULong64_t)1 << 41;
      if (Particle_py < 0)
        Particle_py = (-1) * Particle_py;
      for (Int_t i_bit = 21; i_bit < 41; i_bit++) {
        if ((Particle_py & ((Long64_t)1 << (i_bit - 21))))
          Substitute_p |= (ULong64_t)1 << i_bit;
      }

      Long64_t Particle_pz = (track.pz() * 6000);
      if (Particle_pz < 0)
        Substitute_p |= (ULong64_t)1 << 62;
      if (Particle_pz < 0)
        Particle_pz = (-1) * Particle_pz;
      for (Int_t i_bit = 42; i_bit < 62; i_bit++) {
        if ((Particle_pz & ((Long64_t)1 << (i_bit - 42))))
          Substitute_p |= (ULong64_t)1 << i_bit;
      }

      // dEdx
      UShort_t Substitute_dEdx = (UShort_t)(track.tpcSignal() * 10);

      // DCA
      Short_t Substitute_DCAXY = (Short_t)(track.dcaXY() * 100);
      Short_t Substitute_DCAZ = (Short_t)(track.dcaZ() * 100);

      //--------------- Fill track table ------------------
      testtrack(track.collisionId(),
                track.sign(),
                Substitute_p,
                Substitute_dEdx,
                Substitute_DCAXY,
                Substitute_DCAZ);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<tableDiffWake>(cfgc)};
}
