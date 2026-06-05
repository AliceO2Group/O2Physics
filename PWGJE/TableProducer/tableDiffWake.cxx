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
#include <Framework/AnalysisTask.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>
#include <Framework/AnalysisHelpers.h>
// For centrality:
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
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
namespace o2::aod {
namespace testcol
{
    // Event properties
    DECLARE_SOA_COLUMN(GI, gi, int64_t);
    DECLARE_SOA_COLUMN(RN, rn, int);                 // run number
    DECLARE_SOA_COLUMN(Cent, cent, float);           // FT0C centrality
    DECLARE_SOA_COLUMN(Mult, mult, int);             // TPC multiplicity
    DECLARE_SOA_COLUMN(Occu, occu, int);             // Occupancy ITS
    DECLARE_SOA_COLUMN(OccuFT0, occuft0, float);       // Occupancy FT0C amplitudes
    DECLARE_SOA_COLUMN(VertexX, vertexX, float);
    DECLARE_SOA_COLUMN(VertexY, vertexY, float);
    DECLARE_SOA_COLUMN(VertexZ, vertexZ, float);
    DECLARE_SOA_COLUMN(Psi2, psi2, short);
    DECLARE_SOA_COLUMN(Psi3, psi3, short);
    // Event properties
    /*
    DECLARE_SOA_COLUMN(GI, gi, int64_t);                 // global index of the collision
    DECLARE_SOA_COLUMN(RN, rn, int);                 // run number
    DECLARE_SOA_COLUMN(Cent, cent, float);           // FT0C centrality
    DECLARE_SOA_COLUMN(Mult, mult, int);             // TPC multiplicity
    DECLARE_SOA_COLUMN(VertexX, vertexX, float);
    DECLARE_SOA_COLUMN(VertexY, vertexY, float);
    DECLARE_SOA_COLUMN(VertexZ, vertexZ, float);
    DECLARE_SOA_COLUMN(Psi2, psi2, float);
    DECLARE_SOA_COLUMN(Psi3, psi3, float); */
}
namespace testtrack
{

    // Track properties
    DECLARE_SOA_COLUMN(ColID, colid, int32_t);        // Collision ID
    DECLARE_SOA_COLUMN(Charge, charge, short);
    DECLARE_SOA_COLUMN(P, p, unsigned long);
    DECLARE_SOA_COLUMN(DEdx, dedx, unsigned short);
    DECLARE_SOA_COLUMN(DCAXY, dcaxy, short);
    DECLARE_SOA_COLUMN(DCAZ, dcaz, short);
    //DECLARE_SOA_COLUMN(Length, length, unsigned short);
    // Track properties
    /*
    DECLARE_SOA_COLUMN(ColID, colid, int32_t);        // Collision ID
    DECLARE_SOA_COLUMN(Charge, charge, short);
    DECLARE_SOA_COLUMN(Px, px, float);
    DECLARE_SOA_COLUMN(Py, py, float);
    DECLARE_SOA_COLUMN(Pz, pz, float);
    DECLARE_SOA_COLUMN(DEdx, dedx, float);            // TPC dE/dx
    DECLARE_SOA_COLUMN(DCAXY, dcaxy, float);
    DECLARE_SOA_COLUMN(DCAZ, dcaz, float);
    DECLARE_SOA_COLUMN(Length, length, float);  */      // Track Length
}
DECLARE_SOA_TABLE(TableCol, "AOD", "TABLECOL",
                  testcol::GI,
                  testcol::RN,
                  testcol::Cent,
                  testcol::Mult,
                  testcol::Occu,
                  testcol::OccuFT0,
                  testcol::VertexX,
                  testcol::VertexY,
                  testcol::VertexZ,
                  testcol::Psi2,
                  testcol::Psi3);
DECLARE_SOA_TABLE(TableTrack, "AOD", "TABLETRACK",
                  testtrack::ColID,
                  testtrack::Charge,
                  testtrack::P,
                  testtrack::DEdx,
                  testtrack::DCAXY,
                  testtrack::DCAZ);
                  //testtrack::Length);
}
//--------------------------------------------------------
using namespace o2;
using namespace o2::framework;

struct tableDiffWake {
  // Histogram registry: an object to hold your histograms
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  Configurable<int> nBinsPt{"nBinsPt", 100, "N bins in pT histo"};
  Configurable<double>   pT_thresh{"pT_thresh",20.0,"pT threshold"};
  Configurable<float>    cent_max{"cent_max",10,"centrality"};
  Configurable<float>    z_vert_cut{"z_vert_cut",10.0,"z_vertex cut"};

  Produces<o2::aod::TableCol> testcol;
  Produces<o2::aod::TableTrack> testtrack;

  EventPlaneHelper helperEP;

  void init(InitContext const&)
  {
    // define axes you want to use
    const AxisSpec axisEta{30, -1.5, +1.5, "#eta"};
    const AxisSpec axispT{nBinsPt, 0, 10, "p_{T}"};

    // create histograms
    histos.add("etaHistogram", "etaHistogram", kTH1F, {axisEta});
    histos.add("pTHistogram", "pTHistogram", kTH1F, {axispT});
  }

  using bcs = aod::BCs;
  void process(soa::Join<aod::Collisions, aod::EvSels,aod::CentFT0Cs, aod::TPCMults, aod::QvectorFT0Cs>::iterator const& col,
               soa::Join<aod::TracksIU,aod::TracksExtra, aod::TracksDCA, aod::TrackSelection>  const& tracks,
               bcs const&)
  {
      if(!col.sel8()) return;
      if(col.centFT0C()>cent_max) return;               // Centrality 0 - 10 %
      if(std::abs(col.posZ())>z_vert_cut) return;          // z position < 8 cm
      if(!col.selection_bit(o2::aod::evsel::kNoCollInTimeRangeNarrow)) return;

      //------ Get Run number ---------------------
      auto bc = col.bc_as<bcs>();
      int run = bc.runNumber();
      //int run = 1;
      //------------ EP ---------------------------
      double ep2 = 0.0;
      double ep3 = 0.0;
      ep2 = helperEP.GetEventPlane(col.qvecFT0CRe(), col.qvecFT0CIm(), 2);
      ep3 = helperEP.GetEventPlane(col.qvecFT0CRe(), col.qvecFT0CIm(), 3);

      //------- Only events with track above some thresh ----------

      bool eventHighpT = false;
      for (auto& track : tracks) {

          float dcaXYval = 0.0105*0.035/std::pow(track.pt(),1.1);

          if(!track.hasITS())continue;
          if(!track.hasTPC())continue;
          if(track.tpcCrossedRowsOverFindableCls() < 0.8)continue;
          if(track.tpcFractionSharedCls() > 0.4)continue;
          if(track.tpcNClsCrossedRows() < 70.0)continue;
          if(std::abs(track.dcaXY()) > dcaXYval) continue;         // Different to Run 2
          if(std::abs(track.dcaZ()) > 2.0) continue;          // This is different compared to Run 2
          if(std::abs(track.eta()) > 0.9) continue;
          if(track.tpcChi2NCl() >  4.0) continue;
          if(track.itsChi2NCl() > 36.0) continue;
          if(track.pt() < 0.15) continue; 
          if(track.pt() > pT_thresh){
              eventHighpT = true;
              //LOG(info) << "+++ High pT event found +++";
              break;
          }
      }
      if(!eventHighpT) return;
      //------------------------------------------------------------
      // Translate values to less memory consuming values
      Short_t Substitute_ep2 = (Short_t)(ep2*1000);
      Short_t Substitute_ep3 = (Short_t)(ep3*1000);

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

      //LOG(info) << "Track Loop";
      for (auto& track : tracks) {

          float dcaXYval = 0.0105*0.035/std::pow(track.pt(),1.1);

          // Track cuts
          //if(!track.isGlobalTrack())continue;                // some general track cuts
          if(!track.hasITS())continue;
          if(!track.hasTPC())continue;
          if(track.tpcCrossedRowsOverFindableCls() < 0.8)continue;
          if(track.tpcFractionSharedCls() > 0.4)continue;
          if(track.tpcNClsCrossedRows() < 70.0)continue;
          if(std::abs(track.dcaXY()) > dcaXYval) continue;         // Different to Run 2
          if(std::abs(track.dcaZ()) > 2.0) continue;          // This is different compared to Run 2
          if(std::abs(track.eta()) > 0.9) continue;
          if(track.tpcChi2NCl() >  4.0) continue;
          if(track.itsChi2NCl() > 36.0) continue;
          if(track.pt() < 0.15) continue;

          histos.fill(HIST("etaHistogram"), track.eta());
          histos.fill(HIST("pTHistogram"), track.pt());

          //------------ Translate values to less memory consuming values --------------------
          //Px, Py, Pz
          ULong64_t Substitute_p = 0;
          Long64_t Particle_px = (track.px()*6000);
          if(Particle_px < 0)  Substitute_p |=(ULong64_t)1 << 20;
          if(Particle_px < 0)  Particle_px = (-1)*Particle_px;
          for(Int_t i_bit = 0; i_bit < 20; i_bit++)
          {
              if((Particle_px & ((Long64_t)1 <<  i_bit)))  Substitute_p |= (ULong64_t)1 << i_bit;
          };
          Long64_t Particle_py = (track.py()*6000);
          if(Particle_py < 0)  Substitute_p |=(ULong64_t)1 << 41;
          if(Particle_py < 0)  Particle_py = (-1)*Particle_py;
          for(Int_t i_bit = 21; i_bit < 41 ;i_bit++)
          {
              if((Particle_py & ((Long64_t)1 <<  (i_bit-21))))  Substitute_p |= (ULong64_t)1 << i_bit;
          };
          Long64_t Particle_pz = (track.pz()*6000);
          if(Particle_pz < 0)  Substitute_p |=(ULong64_t)1 << 62;
          if(Particle_pz < 0)  Particle_pz = (-1)*Particle_pz;
          for(Int_t i_bit = 42; i_bit < 62 ;i_bit++)
          {
              if((Particle_pz & ((Long64_t)1 <<  (i_bit-42))))  Substitute_p |= (ULong64_t)1 << i_bit;
          };

          //dEdx
          UShort_t Substitute_dEdx = (UShort_t)(track.tpcSignal()*10);

          //DCA
          Short_t Substitute_DCAXY = (Short_t)(track.dcaXY()*100);
          Short_t Substitute_DCAZ = (Short_t)(track.dcaZ()*100);

          //track length
          //UShort_t Substitute_Length = (UShort_t)(track.length()*10);


          //--------------- Fill track table ------------------
          testtrack(track.collisionId(),
                    track.sign(),     //charge
                    Substitute_p,
                    Substitute_dEdx,    // dE/dx in TPC
                    Substitute_DCAXY,
                    Substitute_DCAZ);
                    //Substitute_Length);
      }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<tableDiffWake>(cfgc)};
}
