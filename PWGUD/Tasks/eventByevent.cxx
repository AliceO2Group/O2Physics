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
//

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "iostream"
#include "PWGUD/DataModel/UDTables.h"
#include <TString.h>
#include <TTree.h>
#include "TLorentzVector.h"
#include "PWGUD/Core/SGSelector.h"
#include "PWGUD/Core/SGTrackSelector.h"

using namespace std;
using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;

/// \brief Event by event study of pions
/// \author Amrit Gautam
/// \author Anisa Khatun
/// \date 20.07.2024

namespace o2::aod
    {
    namespace tree
    {
    DECLARE_SOA_COLUMN(GAPSIDE,gapside,int);
    DECLARE_SOA_COLUMN(FT0AAMP, ft0Aamp,float);        
    DECLARE_SOA_COLUMN(FT0CAMP, ft0Camp,float);
    DECLARE_SOA_COLUMN(FDDAAMP, fddAamp,float);
    DECLARE_SOA_COLUMN(FDDCAMP, fddCamp,float);
    DECLARE_SOA_COLUMN(FV0AAMP, fv0Aamp,float);
    //ZDC tables
    DECLARE_SOA_COLUMN(ZAENERGY, zaenergy,float);        // namespace udzdc
    DECLARE_SOA_COLUMN(ZCENERGY, zcenergy,float);   
    DECLARE_SOA_COLUMN(ZATIME, zatime,float);        
    DECLARE_SOA_COLUMN(ZCTIME, zctime,float);   	 
        //track tables
    DECLARE_SOA_COLUMN(SIGMAPI, sigmapi,std::vector<float>);
    DECLARE_SOA_COLUMN(SIGMAMU, sigmamu,std::vector<float>);    
    DECLARE_SOA_COLUMN(SIGMAK, sigmak,std::vector<float>);
    DECLARE_SOA_COLUMN(SIGMAEL, sigmael,std::vector<float>);
    DECLARE_SOA_COLUMN(SIGMAPR, sigmapr,std::vector<float>);
   /* DECLARE_SOA_COLUMN(PT, Pt,float);
    DECLARE_SOA_COLUMN(RAP, rap,float);
    DECLARE_SOA_COLUMN(PHI, Phi,float);
    DECLARE_SOA_COLUMN(TOTSIGN, totsign,int);
    DECLARE_SOA_COLUMN(MASS , mass,float);*/
    DECLARE_SOA_COLUMN(TOTSIGN, totsign,int);    
    DECLARE_SOA_COLUMN(NPVTRACK, npvtrack,int);
    //Vertex xyz
    DECLARE_SOA_COLUMN(XVTX, xvtx,float);
    DECLARE_SOA_COLUMN(YVTX, yvtx,float);
    DECLARE_SOA_COLUMN(ZVTX, zvtx,float);
    // tracks kinematics	
    DECLARE_SOA_COLUMN(PXS, pxs,std::vector<float>);
    DECLARE_SOA_COLUMN(PYS, pys,std::vector<float>);
    DECLARE_SOA_COLUMN(PZS, pzs,std::vector<float>);
    DECLARE_SOA_COLUMN(SIGNS, Signs,std::vector<float>);
    DECLARE_SOA_COLUMN(UPCFLAG, upcflag,int);
    DECLARE_SOA_COLUMN(RUNNUM, runnum,int);    
        
   // DECLARE_SOA_COLUMN(RAWTRACKS, rawtracks,int);
   // DECLARE_SOA_COLUMN(PTRACKS, ptracks,int);
    
    }//end namespace

DECLARE_SOA_TABLE(TREE, "AOD", "Tree", //! ZDC information
                  tree::GAPSIDE,
                  tree::FT0AAMP,
                  tree::FT0CAMP,
                  tree::FDDAAMP,
                  tree::FDDCAMP,
                  tree::FV0AAMP,
                  tree::ZAENERGY,
                  tree::ZCENERGY,
		  tree::ZATIME,
		  tree::ZCTIME,
                 /* tree::PT,
                  tree::RAP,
                  tree::PHI,
                  tree::MASS,
                  tree::TOTSIGN,*/
                  tree::TOTSIGN,
                  tree::NPVTRACK,
                  tree::XVTX,
                  tree::YVTX,
                  tree::ZVTX,
                  tree::SIGMAPI,
                   tree::SIGMAMU,
                  tree::SIGMAK,
                  tree::SIGMAEL,
                  tree::SIGMAPR,
                  tree::PXS,
                  tree::PYS,
                  tree::PZS,
                  tree::SIGNS,
                  tree::UPCFLAG,
                  tree::RUNNUM
                  /*,
                  tree::RAWTRACKS,
                  tree::PTRACKS*/);

} // namespace o2::aod


struct EventByEvent {
  SGSelector sgSelector;
  Produces<aod::TREE> tree;

  Configurable<float> FV0_cut{"FV0", 100., "FV0A threshold"};
  Configurable<float> FT0A_cut{"FT0A", 200., "FT0A threshold"};
  Configurable<float> FT0C_cut{"FT0C", 100., "FT0C threshold"};
  Configurable<float> FDDA_cut{"FDDA", 10000., "FDDA threshold"};
  Configurable<float> FDDC_cut{"FDDC", 10000., "FDDC threshold"};
  Configurable<float> ZDC_cut{"ZDC", 100., "ZDC threshold"};
  Configurable<float> gap_Side{"gap", 2, "gap selection"};

  // Collision selection
  Configurable<float> collcontrib_cut{"collcontrib_cut", 10, "no. of PV contributor per collsion"};
  Configurable<float> Zvtx_cut{"Zvtx_cut", 15, "z-vertex selection"};

  // Track Selections
  Configurable<float> PV_cut{"PV_cut", 1.0, "Use Only PV tracks"};
  Configurable<float> dcaZ_cut{"dcaZ_cut", 2.0, "dcaZ cut"};
  Configurable<float> dcaXY_cut{"dcaXY_cut", 0.0, "dcaXY cut (0 for Pt-function)"};
  Configurable<float> tpcChi2_cut{"tpcChi2_cut", 4, "Max tpcChi2NCl"};
  Configurable<float> tpcNClsFindable_cut{"tpcNClsFindable_cut", 70, "Min tpcNClsFindable"};
  Configurable<float> itsChi2_cut{"itsChi2_cut", 36, "Max itsChi2NCl"};
  Configurable<float> eta_cut{"eta_cut", 0.9, "Track Pseudorapidity"};
  Configurable<float> pt_cut{"pt_cut", 0.1, "Track Pt"};
  Configurable<float> TPC_cluster{"TPC_cluster", 50, "No.of TPC cluster"};

  // Kinmatic cuts
  Configurable<float> PID_cut{"PID_cut", 5, "TPC PID"};
  Configurable<float> Rap_cut{"Rap_cut", 0.9, "Track rapidity"};
  Configurable<float> Mass_Max{"Mass_Max", 10, "Invariant Mass range high"};
  Configurable<float> Mass_Min{"Mass_Min", 0, "Invariant Mass range low"};
  Configurable<float> Pt_coherent{"Pt_coherent", 0.15, "Coherent selection"};

  // defining histograms using histogram registry
  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};

  //-----------------------------------------------------------------------------------------------------------------------
  void init(o2::framework::InitContext&)
  {

    registry.add("GapSide", "Gap Side; Entries", kTH1F, {{4, -1.5, 2.5}});
    registry.add("TrueGapSide", "Gap Side; Entries", kTH1F, {{4, -1.5, 2.5}});

    auto hSelectionCounter = registry.add<TH1>("hSelectionCounter", "hSelectionCounter;;NEvents", HistType::kTH1I, {{20, 0., 20.}});

    TString SelectionCuts[18] = {"NoSelection", "gapside", "goodtracks", "truegap", "ncollcontrib ", "zvtx", "2collcontrib", "ngoodtrk", "TPCPID", "Rap_cut", "unlikesign", "mass_cut", "coherent", "incoherent", "likesign", "mass_cut", "coherent", "incoherent"};
    // now we can set BinLabel in histogram Registry

    for (int i = 0; i < 18; i++) {
      hSelectionCounter->GetXaxis()->SetBinLabel(i + 1, SelectionCuts[i].Data());
    }
   
    // tracks
    registry.add("hTracks", "N_{tracks}", kTH1F, {{100, -0.5, 99.5}});
    registry.add("hQualityTracks", "N_{tracks}", kTH1F, {{100, -0.5, 99.5}});  
    //registry.add("h4TracksMuons", "N_{tracks}", kTH1F, {{100, -0.5, 99.5}});
    registry.add("hdEdx", "p vs dE/dx Signal", kTH2F, {{100, 0.0, 3.0}, {100, 0.0, 200.0}});
    //registry.add("hdEdxMuons", "p_{#pi} vs dE/dx Signal", kTH2F, {{100, 0.0, 3.0}, {100, 0.0, 200.0}});
    
    //using Angular Correlation method

  }

  using udtracks = soa::Join<aod::UDTracks, aod::UDTracksExtra, aod::UDTracksPID>;
  using udtracksfull = soa::Join<aod::UDTracks, aod::UDTracksPID, aod::UDTracksExtra, aod::UDTracksFlags, aod::UDTracksDCA>;
  using UDCollisionsFull = soa::Join<aod::UDCollisions, aod::SGCollisions, aod::UDCollisionsSels, aod::UDZdcsReduced>;
  //__________________________________________________________________________
  // Main process
  void process(UDCollisionsFull::iterator const& collision, udtracksfull const& tracks)
  {
    registry.fill(HIST("hSelectionCounter"), 0);
    // LOGF(info, " BC ID %d",collision.gapSide());
    int gapSide = collision.gapSide();
    if (gapSide < 0 || gapSide > 2)
      return;

    registry.fill(HIST("hSelectionCounter"), 1);

    float FIT_cut[5] = {FV0_cut, FT0A_cut, FT0C_cut, FDDA_cut, FDDC_cut};
    int truegapSide = sgSelector.trueGap(collision, FIT_cut[0], FIT_cut[1], FIT_cut[2], ZDC_cut);
    std::vector<float> parameters = {PV_cut, dcaZ_cut, dcaXY_cut, tpcChi2_cut, tpcNClsFindable_cut, itsChi2_cut, eta_cut, pt_cut};
    registry.fill(HIST("hSelectionCounter"), 2);

    registry.fill(HIST("GapSide"), gapSide);
    registry.fill(HIST("TrueGapSide"), truegapSide);
    gapSide = truegapSide;


      registry.fill(HIST("hSelectionCounter"), 3);
      //_____________________________________
      // Create pions and apply TPC Pion PID
      std::vector<TLorentzVector> allTracks;
      //std::vector<TLorentzVector> onlyPionTracks;
      //std::vector<float> onlyPionSigma;
      std::vector<decltype(tracks.begin())> rawGoodTracks;
      std::vector<float> trackpx;
      std::vector<float> trackpy;
      std::vector<float> trackpz;
      std::vector<float> tracksign;
      std::vector<float> pitpcpid;
      std::vector<float> mutpcpid;
      std::vector<float> ktpcpid;
      std::vector<float> eltpcpid;
      std::vector<float> prtpcpid;
      
      TLorentzVector p;
      
      if (collision.numContrib()<2 || collision.numContrib() > collcontrib_cut)
        return;

      
      
      registry.fill(HIST("hSelectionCounter"), 4); //2 - 10 collision con
      if ((collision.posZ() < -(Zvtx_cut)) || (collision.posZ() > Zvtx_cut))
        return;
      registry.fill(HIST("hSelectionCounter"), 5);// poz 
      registry.fill(HIST("hTracks"), tracks.size());  
      for (auto t : tracks) {
          
      
        if(!trackselector(t, parameters)) continue;  
        double dEdx = t.tpcSignal();
        registry.fill(HIST("hdEdx"), t.tpcInnerParam() / t.sign(), dEdx);
        
        TLorentzVector a;
        a.SetXYZM(t.px(), t.py(), t.pz(), o2::constants::physics::MassPionCharged);
        allTracks.push_back(a);
        auto nSigmaPi = t.tpcNSigmaPi();

        //if (fabs(nSigmaMu) > PID_cut) continue; 
          //onlyPionTracks.push_back(a);
          //onlyPionSigma.push_back(nSigmaPi);
          rawGoodTracks.push_back(t);
          //registry.fill(HIST("hdEdxMuon"), t.tpcInnerParam() / t.sign(), dEdx);
        
      }
      
      //_____________________________________
      
      registry.fill(HIST("hQualityTracks"), rawGoodTracks.size());
          
        if (rawGoodTracks.size() >= 2) {
           
            
            registry.fill(HIST("hSelectionCounter"), 7);
           
            int sign = 0;
            for(auto rtrk : rawGoodTracks) {
                //**********
                //to do
               //1.  make vector of vector which store piplus and piminus and its values
               //2.  save px py pz to add different channels and avoid adding mass.
               //******************
               
                sign += rtrk.sign();
                trackpx.push_back(rtrk.px());
                trackpy.push_back(rtrk.py());
                trackpz.push_back(rtrk.pz());
                tracksign.push_back(rtrk.sign());
                pitpcpid.push_back(rtrk.tpcNSigmaPi());
                mutpcpid.push_back(rtrk.tpcNSigmaMu());
                
                ktpcpid.push_back(rtrk.tpcNSigmaKa());
                eltpcpid.push_back(rtrk.tpcNSigmaEl());
                prtpcpid.push_back(rtrk.tpcNSigmaPr());
            }
           registry.fill(HIST("hTracks"), collision.numContrib());
                   
            tree(gapSide,collision.totalFT0AmplitudeA(),collision.totalFT0AmplitudeC(),collision.totalFDDAmplitudeA(),collision.totalFDDAmplitudeC(),collision.totalFV0AmplitudeA(),collision.energyCommonZNA(),collision.energyCommonZNC(),collision.timeZNA(),collision.timeZNC(),sign,collision.numContrib(),collision.posX(),collision.posY(),collision.posZ(),pitpcpid,mutpcpid,ktpcpid,eltpcpid,prtpcpid, trackpx,trackpy,trackpz,tracksign,collision.flags(),collision.runNumber());
        }
      }
      
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<EventByEvent>(cfgc)};
}
