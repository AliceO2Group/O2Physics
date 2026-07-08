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


/// \file  MultandPtFluctuations.cxx
/// \brief Calculate multiplicity and transverse momentum fluctuations using strongly intensive observables
/// \author Omama Rubza <omama.rubza@cern.ch>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/HistogramSpec.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/FT0Corrected.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "TProfile.h"
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::aod;
using namespace o2::framework::expressions;


struct myExampleTask {

  // ------ Histogram binning

  Configurable<int> nBinsPt{"nBinsPt", 100, "pT bins"};
  Configurable<int> ndcaXY{"ndcaXY", 100, "DCAxy bins"};
  Configurable<int> ndcaZ{"ndcaZ", 100, "DCAz bins"};
  Configurable<int> nCentBins{"nCentBins", 100, "Number of centrality bins"};


  // ------ Event Cuts
  
  Configurable<float> vtxZcut{"vertexZcut", 10.0f, "Vertex Z"};
  Configurable<bool> cfgNoSameBunchPileup{"cfgNoSameBunchPileup", true, "kNoSameBunchPileup"};
  Configurable<bool> cfgEvSelUseGoodZvtxFT0vsPV{"cfgEvSelUseGoodZvtxFT0vsPV",true, "Good Zvtx FT0 vs PV"};
  Configurable<bool> cfgUseGoodITSLayerAllCut{"cfgUseGoodITSLayerAllCut", true, "Good ITS Layers"};

 
  // ----- Centrality estimator

  Configurable<bool> cFT0M{"cFT0M", false, "FT0M centrality"}; // for pp
  Configurable<bool> cFT0C{"cFT0C", true,  "Use FT0C centrality"};


  // ------ Track cuts

  Configurable<float> ptMinCut{"ptMinCut", 0.2f, "Minimum pT"};
  Configurable<float> ptMaxCut{"ptMaxCut", 5.0f, "Maximum pT"};
  Configurable<float> etaCut{"etaCut",0.8f, "Maximum |eta|"};
  Configurable<int> itsNClsCut{"itsNClsCut", 4, "Minimum ITS clusters"};
  Configurable<float> tpcCrossCut{"tpcCrossCut", 70.0f, "TPC crossed rows"};
  Configurable<float> crossedRowsOverFindableCut{"crossedRowsOverFindableCut",  0.8f,"TPC crossed rows over findable"};
  Configurable<float> dcaZCut{"dcaZCut", 2.0f, "Maximum DCAz"};
  Configurable<float> dcaXYCut{"dcaXYCut", 0.2f, "Maximum DCAxy"};
  Configurable<float> tpcChiCut{"tpcChiCut", 4.0f, "TPC chi2/NCls"};
  Configurable<float> itsChiCut{"itsChiCut", 36.0f,"ITS chi2/NCls"};
  Configurable<bool> requireITS{"requireITS",true, "Require ITS hit"};
  Configurable<bool> requireTPC{"requireTPC", true, "Require TPC hit"};
  Configurable<bool> requireInnerITS{"requireInnerITS", true,"At least one hit in ITS layers 0,1,2"};
  


  //----- Histogram Registry

  HistogramRegistry histos{
    "histos", {}, OutputObjHandlingPolicy::AnalysisObject
  };

  
  void init(InitContext const&)
  {
    const AxisSpec axisEta{30, -1.5, +1.5, "#eta"};
    const AxisSpec axisPt{nBinsPt, 0.15, 10.0, "p_{T}"};
    const AxisSpec axisVtxz{80, -20, 20, "V_{Z} (cm)"};
    const AxisSpec axisPhi{40, -1, 7, "#phi"};
    const AxisSpec axisdcaXY{ndcaXY, -0.3 , 0.3, "DCAxy"};
    const AxisSpec axisdcaZ{ndcaZ, -3.0 , 3.0 , "DCAz"};
    const AxisSpec axisNch{500, 0, 500, "Nch"};
    const AxisSpec axisCent{ nCentBins, 0.0, 100.0, "Centrality (%)"}; // ------------- centrality bins

  
    histos.add("QA/BeforeCut/VtxZ", "Vertex Z", kTH1F, {axisVtxz});
    histos.add("QA/AfterCut/VtxZ",  "Vertex Z", kTH1F, {axisVtxz});

    histos.add("QA/BeforeCut/Cent", "Centrality", kTH1F, {axisCent});
    histos.add("QA/AfterCut/Cent",  "Centrality", kTH1F, {axisCent});

    histos.add("QA/BeforeCut/Eta",   "Eta",   kTH1F, {axisEta});
    histos.add("QA/AfterCut/Eta",    "Eta",   kTH1F, {axisEta});

    histos.add("QA/BeforeCut/Pt",    "Pt",    kTH1F, {axisPt});
    histos.add("QA/AfterCut/Pt",     "Pt",    kTH1F, {axisPt});

    histos.add("QA/BeforeCut/Phi",   "Phi",   kTH1F, {axisPhi});
    histos.add("QA/AfterCut/Phi",    "Phi",   kTH1F, {axisPhi});

    histos.add("QA/BeforeCut/DcaXY", "DCAxy", kTH1F, {axisdcaXY});
    histos.add("QA/AfterCut/DcaXY",  "DCAxy", kTH1F, {axisdcaXY});

    histos.add("QA/BeforeCut/DcaZ",  "DCAz",  kTH1F, {axisdcaZ});
    histos.add("QA/AfterCut/DcaZ",   "DCAz",  kTH1F, {axisdcaZ}); 
    
    histos.add("hNch" , "Nch" , kTH1F, {axisNch});
    histos.add("h2_DcaZ", "DCA_{Z}", kTH2D, {{axisPt}, {axisdcaZ}});
    histos.add("h2_DcaXY", "DCA_{XY}", kTH2D, {{axisPt}, {axisdcaXY}});

    histos.add("hEventCounter", "Number of events vs centrality", kTH1F, {axisCent});

    histos.add("p_a"  , "<A> vs centrality"    , kTProfile, {axisCent});
    histos.add("p_b"  , "<B> vs centrality"    , kTProfile, {axisCent});
    histos.add("p_a2" , "<A^{2}> vs centrality", kTProfile, {axisCent});
    histos.add("p_b2" , "<B^{2}> vs centrality", kTProfile, {axisCent});
    histos.add("p_ab", "<AB> vs centrality" , kTProfile, {axisCent});
    histos.add("p_asumb", "<A+B> vs centrality" , kTProfile, {axisCent});

  }

 //---------- Filters 

  Filter ptFilter =(aod::track::pt > ptMinCut) && (aod::track::pt < ptMaxCut);
  Filter etaFilter =nabs(aod::track::eta) < etaCut;
  Filter posZFilter = nabs(aod::collision::posZ) < vtxZcut;
 

  using myColsData = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms,aod::CentFT0Cs>;
  using myTracksData = soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra, aod::TracksDCA>;
  using myFilteredTracksData = soa::Filtered<myTracksData>;

   void process(myColsData::iterator const& col, myFilteredTracksData const& tracks)
  {
    histos.fill(HIST("QA/BeforeCut/VtxZ"), col.posZ());

    float cent = -1.0f;

    if (cFT0M) { 
        cent = col.centFT0M();
    }

    if (cFT0C) {
       cent = col.centFT0C();
   // if (!track.isGlobalTrack()) continue;
  }

    if (cent < 0) {
      return;
    }

    histos.fill(HIST("QA/BeforeCut/Cent"), cent);
    
    if (!col.sel8()) return;
    
    if (std::abs(col.posZ()) > vtxZcut) return;

    if (cfgNoSameBunchPileup && !col.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) return;
  
    if (cfgEvSelUseGoodZvtxFT0vsPV && !col.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) return;

    if (cfgUseGoodITSLayerAllCut && !col.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll))  return;

    histos.fill(HIST("hEventCounter"), cent);

    histos.fill(HIST("QA/AfterCut/VtxZ"), col.posZ());
    histos.fill(HIST("QA/AfterCut/Cent"), cent);
    
    double A = 0.0;
    double B = 0.0;  
 
   // if (!track.isGlobalTrack()) continue;
   for (auto& track : tracks) {

    histos.fill(HIST("QA/BeforeCut/Eta"), track.eta());
    histos.fill(HIST("QA/BeforeCut/Pt"), track.pt());
    histos.fill(HIST("QA/BeforeCut/Phi"), track.phi());
    histos.fill(HIST("QA/BeforeCut/DcaXY"), track.dcaXY());
    histos.fill(HIST("QA/BeforeCut/DcaZ"), track.dcaZ());

    // Track quality cuts
    
    if (requireITS && !track.hasITS()) continue;
    if (requireTPC && !track.hasTPC()) continue;

    if (track.itsNCls() < itsNClsCut) continue;
    if (track.tpcNClsCrossedRows() < tpcCrossCut) continue;

    if (track.tpcCrossedRowsOverFindableCls() < crossedRowsOverFindableCut) continue;

    if (std::abs(track.dcaZ()) > dcaZCut) continue;
    if (std::abs(track.dcaXY()) > dcaXYCut) continue;

    if (track.tpcChi2NCl() > tpcChiCut) continue;
    if (track.itsChi2NCl() > itsChiCut) continue;

    if (std::abs(track.eta()) > etaCut) continue;

    if (requireInnerITS) {

      auto itsMap = track.itsClusterMap();

      if (!(itsMap & (1 << 0)) &&
          !(itsMap & (1 << 1)) &&
          !(itsMap & (1 << 2))) {
      continue;
    }
  }
      

    /*  if (eta >= 0.6 && eta < 0.8)
        nf++;
      else if (eta > -0.8 && eta <= -0.6)Sigma_pTN_OO_NeNe.C
        nb++;*/

    // After cuts QA

    histos.fill(HIST("QA/AfterCut/Eta"), track.eta());
    histos.fill(HIST("QA/AfterCut/Pt"), track.pt());
    histos.fill(HIST("QA/AfterCut/Phi"), track.phi());
    histos.fill(HIST("QA/AfterCut/DcaXY"), track.dcaXY());
    histos.fill(HIST("QA/AfterCut/DcaZ"), track.dcaZ());

    histos.fill(HIST("h2_DcaZ"), track.pt(),track.dcaZ());
    histos.fill(HIST("h2_DcaXY"), track.pt(),track.dcaXY());
    
        A++;			//A-nch, B =pt
	B += track.pt();
     
    }
    //  LOG(info) << "Nch = " << nch;
  /*    histos.fill(HIST("hNch"),nch);
      eventNch(col, nch);*/


      // ---- Fill TProfiles (once per event)
      histos.fill(HIST("p_a"), cent, A);
      histos.fill(HIST("p_b"), cent, B);

      histos.fill(HIST("p_a2"),cent, A * A);
      histos.fill(HIST("p_b2"), cent, B * B);

      histos.fill(HIST("p_ab"), cent, A * B);
      histos.fill(HIST("p_asumb"), cent, A + B);


  }
};


    WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
    {
      return WorkflowSpec{
      adaptAnalysisTask<myExampleTask>(cfgc)
  };
}




