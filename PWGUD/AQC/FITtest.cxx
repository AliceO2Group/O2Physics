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
/// \brief A task for testing FIT selection for Ultra-perimpheral Collisions
/// \author Anisa Khatun, anisa.khatun@cern.ch
/// \since  04.08.2023

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "ReconstructionDataFormats/BCRange.h"
#include "CommonConstants/PhysicsConstants.h"
#include "Common/DataModel/FT0Corrected.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/PIDResponse.h"
#include "CommonConstants/LHCConstants.h"

#include "Framework/StaticFor.h"
#include "TLorentzVector.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace std;

struct FITtest {
  
    // inivinitialize HistogramRegistry
  HistogramRegistry registry{
    "registry",
    {}
      };

  // define abbreviations
  using CCs = soa::Join<aod::Collisions, aod::EvSels>;
  using CC = CCs::iterator;
  using BCs = soa::Join<aod::BCs, aod::BcSels, aod::Run3MatchedToBCSparse>;
  using TCs = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA, aod::TrackSelectionExtension, aod::pidTPCFullPi, aod::TOFSignal, aod::pidTOFbeta>;
  using FWs = aod::FwdTracks;
  using ATs = aod::AmbiguousTracks;
  using AFTs = aod::AmbiguousFwdTracks;
  

  void init(InitContext& context)
  {
    // add histograms for the different process functions
    if (context.mOptions.get<bool>("processMain")) {

// collisions
registry.add("collisions/Stat", "Cut statistics; Selection criterion; Collisions", {HistType::kTH1F, {{20, -0.5, 20.5}}});
registry.add("collisions/GlobalBC", "Relative BC no. in FIT; Relative BC; Collisions", {HistType::kTH1F, {{3564, -0.5, 3563.5}}});
registry.add("collisions/RelativeBC", "Relative BC number; Relative BC; Collisions", {HistType::kTH1F, {{3564, -0.5, 3563.5}}});
registry.add("collisions/trkmultiplicity", "Multiplicity of all tracks; Multiplicity; Tracks", {HistType::kTH1F, {{300, -0.5, 299.5}}});
registry.add("collisions/PVCFIT", "PV contributors with FIT; PV contributors; Tracks", {HistType::kTH1F, {{100, -0.5, 99.5}}});
registry.add("collisions/vtxTracks", "Number of vertex tracks; Number of contributors; Collisions", {HistType::kTH1F, {{300, -0.5, 299.5}}});
registry.add("collisions/PVTracks", "Number of PV tracks; Number of PV tracks; Collisions", {HistType::kTH1F, {{300, -0.5, 299.5}}});

// FV0
registry.add("FV0/hV0A", "Time FV0A; Time (ns); Collisions", {HistType::kTH1F, {{500, -50.0, 50.0}}});
registry.add("FV0/FV0Amp", "FV0A Amplitude", {HistType::kTH1F, {{2000, 0., 2000.}}});
registry.add("FV0/FV0A", "#FV0A; Channel; FV0A Amplitude", {HistType::kTH2F, {{48, -0.5, 47.5}, {2000, 0., 2000.}}});

// FT0
registry.add("FT0/hT0A", "Time FT0 A side; Time (ns); Collisions", {HistType::kTH1F, {{500, -50.0, 50.0}}});
registry.add("FT0/hT0C", "Time FT0 C side; Time (ns); Collisions", {HistType::kTH1F, {{500, -50.0, 50.0}}});
registry.add("FT0/hT0AC", "Time Correlation FT0; FT0A Time (ns) ; FT0C Time (ns)", {HistType::kTH2F, {{500, -50.0, 50.0},{500, -50.0, 50.0}}});
registry.add("FT0/FT0Aamp", "#FT0A Amplitude; FT0A Amplitude", {HistType::kTH1F, {{200, -0.5, 199.5}}});
registry.add("FT0/FT0Camp", "#FT0C Amplitude; FT0C Amplitude", {HistType::kTH1F, {{200, -0.5, 199.5}}});
registry.add("FT0/FT0A", "#FT0A; Channel; FT0A Amplitude", {HistType::kTH2F, {{96, -0.5, 95.5}, {1000, 0., 1000.}}});
registry.add("FT0/FT0C", "#FT0C; Channel; FT0C Amplitude", {HistType::kTH2F, {{112, -0.5, 111.5}, {1000, 0., 1000.}}});
registry.add("FT0/FT0ACCorr", "FT0 amp correlation; FT0A Amplitude; FT0C Amplitude", {HistType::kTH2F, {{2000, -4.5, 1995.5}, {2000, -4.5, 1995.5}}});

// FDD
registry.add("FDD/hFDDA", " Avg Time FDD A side; Time (ns); Collisions", {HistType::kTH1F, {{500, -50.0, 50.0}}});
registry.add("FDD/hFDDC", " Avg Time FDD C side; Time (ns); Collisions", {HistType::kTH1F, {{500, -50.0, 50.0}}});
registry.add("FDD/hFDDAC", "Time Correlation FDD; FDDA time (ns); FDDC time (ns)", {HistType::kTH2F, {{500, -50.0, 50.0},{500, -50.0, 50.0}}});
registry.add("FDD/FDDAamp", "#FDDA Amplitude; FDDA Amplitude", {HistType::kTH1F, {{200, -0.5, 199.5}}});
registry.add("FDD/FDDCamp", "#FDDC Amplitude; FDDC Amplitude", {HistType::kTH1F, {{200, -0.5, 199.5}}});
registry.add("FDD/FDDA", "#FDDA; Channel; FDDA Amplitude", {HistType::kTH2F, {{8, -0.5, 7.5}, {1000, 0., 1000.}}});
registry.add("FDD/FDDC", "#FDDC; Channel; FDDC Amplitude", {HistType::kTH2F, {{8, -0.5, 7.5}, {1000, 0., 1000.}}});
registry.add("FDD/FDDACCorr", "#FDD amp correlation; FDDA Amplitude; FDDC Amplitude", {HistType::kTH2F, {{2000, -4.5, 1995.5}, {2000, -4.5, 1995.5}}});

// Correlation plots
registry.add("collisions/FV0T0ACorr", "Correlation FV0 vs FT0 A side; FT0A Amplitude; FV0A Amplitude", {HistType::kTH2F, {{2000, -4.5, 1995.5}, {2000, -4.5, 1995.5}}});
registry.add("collisions/FV0T0CCorr", "Correlation FV0 vs FT0 C side; FT0C Amplitude; FV0A Amplitude", {HistType::kTH2F, {{2000, -4.5, 1995.5}, {2000, -4.5, 1995.5}}});
registry.add("collisions/FT0DDACorr", "Correlation FT0 vs FDD A side; FT0A Amplitude; FDDA Amplitude", {HistType::kTH2F, {{2000, -4.5, 1995.5}, {2000, -4.5, 1995.5}}});
registry.add("collisions/FT0DDCCorr", "Correlation FT0 vs FDD C side; FT0C Amplitude; FDDC Amplitude", {HistType::kTH2F, {{2000, -4.5, 1995.5}, {2000, -4.5, 1995.5}}});
registry.add("collisions/FT0AFDDC", "Correlation FT0 vs FDD AC side; FT0A Amplitude; FDDC Amplitude", {HistType::kTH2F, {{2000, -4.5, 1995.5}, {2000, -4.5, 1995.5}}});
registry.add("collisions/FT0CFDDA", "Correlation FT0 vs FDD CA side; FT0C Amplitude; FDDA Amplitude", {HistType::kTH2F, {{2000, -4.5, 1995.5}, {2000, -4.5, 1995.5}}});
registry.add("collisions/FV0AFDDA", "Correlation FV0 vs FDD A side; FV0A Amplitude; FDDA Amplitude", {HistType::kTH2F, {{2000, -4.5, 1995.5}, {2000, -4.5, 1995.5}}});
registry.add("collisions/FV0AFDDC", "Correlation FV0 vs FDD C side; FV0A Amplitude; FDDC Amplitude", {HistType::kTH2F, {{2000, -4.5, 1995.5}, {2000, -4.5, 1995.5}}});
        }
      
    if (context.mOptions.get<bool>("processHadronic")) {
          
     registry.add("collHadronic/StatHadronic", "Cut statistics; Selection criterion; Collisions", {HistType::kTH1F, {{20, -0.5, 20.5}}});
     registry.add("collHadronic/GlBCHadronic", "Relative BC no. in FIT; Relative BC; Collisions", {HistType::kTH1F, {{3564, -0.5, 3563.5}}});
     registry.add("collHadronic/RelBCHadronic", "Relative BC number; Relative BC; Collisions", {HistType::kTH1F, {{3564, -0.5, 3563.5}}});
     registry.add("collHadronic/trkmultHadronic", "Multiplicity of all tracks; Multiplicity; Tracks", {HistType::kTH1F, {{300, -0.5, 299.5}}});
     registry.add("collHadronic/PVCFITHadronic", "PV contributors with FIT; PV contributors; Tracks", {HistType::kTH1F, {{100, -0.5, 99.5}}});
     registry.add("collHadronic/vtxTrkHadronic", "Number of vertex tracks; Number of contributors; Collisions", {HistType::kTH1F, {{300, -0.5, 299.5}}});
     registry.add("collHadronic/PVTrkHadronic", "Number of PV tracks; Number of PV tracks; Collisions", {HistType::kTH1F, {{300, -0.5, 299.5}}});
          
          // FV0
          registry.add("collHadronic/FV0/hV0AHadronic", "Time FV0A; Time (ns); Collisions", {HistType::kTH1F, {{500, -50.0, 50.0}}});
          registry.add("collHadronic/FV0/FV0AmpHadronic", "FV0A Amplitude", {HistType::kTH1F, {{2000, 0., 2000.}}});
          
          // FT0
          registry.add("collHadronic/FT0/hT0AHadronic", "Time FT0 A side; Time (ns); Collisions", {HistType::kTH1F, {{500, -50.0, 50.0}}});
          registry.add("collHadronic/FT0/hT0CHadronic", "Time FT0 C side; Time (ns); Collisions", {HistType::kTH1F, {{500, -50.0, 50.0}}});
          registry.add("collHadronic/FT0/hT0ACHadronic", "Time Correlation FT0; FT0A Time (ns) ; FT0C Time (ns)", {HistType::kTH2F, {{500, -50.0, 50.0},{500, -50.0, 50.0}}});
          registry.add("collHadronic/FT0/FT0AampHadronic", "#FT0A Amplitude; FT0A Amplitude", {HistType::kTH1F, {{200, -0.5, 199.5}}});
          registry.add("collHadronic/FT0/FT0CampHadronic", "#FT0C Amplitude; FT0C Amplitude", {HistType::kTH1F, {{200, -0.5, 199.5}}});
          registry.add("collHadronic/FT0/FT0ACCorrHadronic", "FT0 amp correlation; FT0A Amplitude; FT0C Amplitude", {HistType::kTH2F, {{2000, -4.5, 1995.5}, {2000, -4.5, 1995.5}}});

          // FDD
          registry.add("collHadronic/FDD/hFDDAHadronic", " Avg Time FDD A side; Time (ns); Collisions", {HistType::kTH1F, {{500, -50.0, 50.0}}});
          registry.add("collHadronic/FDD/hFDDCHadronic", " Avg Time FDD C side; Time (ns); Collisions", {HistType::kTH1F, {{500, -50.0, 50.0}}});
          registry.add("collHadronic/FDD/hFDDACHadronic", "Time Correlation FDD; FDDA time (ns); FDDC time (ns)", {HistType::kTH2F, {{500, -50.0, 50.0},{500, -50.0, 50.0}}});
          registry.add("collHadronic/FDD/FDDAampHadronic", "#FDDA Amplitude; FDDA Amplitude", {HistType::kTH1F, {{200, -0.5, 199.5}}});
          registry.add("collHadronic/FDD/FDDCampHadronic", "#FDDC Amplitude; FDDC Amplitude", {HistType::kTH1F, {{200, -0.5, 199.5}}});
          registry.add("collHadronic/FDD/FDDACCorrHadronic", "#FDD amp correlation;FDDA Amplitude; FDDC Amplitude", {HistType::kTH2F, {{2000, -4.5, 1995.5}, {2000, -4.5, 1995.5}}});

          // Correlation plots
          registry.add("collHadronic/FV0T0ACorrHadronic", "Correlation FV0 vs FT0 A side; FT0A Amplitude; FV0A Amplitude", {HistType::kTH2F, {{2000, -4.5, 1995.5}, {2000, -4.5, 1995.5}}});
          registry.add("collHadronic/FV0T0CCorrHadronic", "Correlation FV0 vs FT0 C side; FT0C Amplitude; FV0A Amplitude", {HistType::kTH2F, {{2000, -4.5, 1995.5}, {2000, -4.5, 1995.5}}});
          registry.add("collHadronic/FT0DDACorrHadronic", "Correlation FT0 vs FDD A side; FT0A Amplitude; FDDA Amplitude", {HistType::kTH2F, {{2000, -4.5, 1995.5}, {2000, -4.5, 1995.5}}});
          registry.add("collHadronic/FT0DDCCorrHadronic", "Correlation FT0 vs FDD C side; FT0C Amplitude; FDDC Amplitude", {HistType::kTH2F, {{2000, -4.5, 1995.5}, {2000, -4.5, 1995.5}}});
          registry.add("collHadronic/FT0AFDDCHadronic", "Correlation FT0 vs FDD AC side; FT0A Amplitude; FDDC Amplitude", {HistType::kTH2F, {{2000, -4.5, 1995.5}, {2000, -4.5, 1995.5}}});
          registry.add("collHadronic/FT0CFDDAHadronic", "Correlation FT0 vs FDD CA side; FT0C Amplitude; FDDA Amplitude", {HistType::kTH2F, {{2000, -4.5, 1995.5}, {2000, -4.5, 1995.5}}});
          registry.add("collHadronic/FV0AFDDAHadronic", "Correlation FV0 vs FDD A side; FV0A Amplitude; FDDA Amplitude", {HistType::kTH2F, {{2000, -4.5, 1995.5}, {2000, -4.5, 1995.5}}});
          registry.add("collHadronic/FV0AFDDCHadronic", "Correlation FV0 vs FDD C side; FV0A Amplitude; FDDC Amplitude", {HistType::kTH2F, {{2000, -4.5, 1995.5}, {2000, -4.5, 1995.5}}});
      }
      
      if (context.mOptions.get<bool>("processInclusiveA")) {
          
          registry.add("colInclusive/Stat", "Cut statistics; Selection criterion; Collisions", {HistType::kTH1F, {{20, -0.5, 20.5}}});
          registry.add("colInclusive/GlobalBC", "Relative BC no. in FIT; Relative BC; Collisions", {HistType::kTH1F, {{3564, -0.5, 3563.5}}});
          registry.add("colInclusive/RelativeBC", "Relative BC number; Relative BC; Collisions", {HistType::kTH1F, {{3564, -0.5, 3563.5}}});
          registry.add("colInclusive/trkmultiplicity", "Multiplicity of all tracks; Multiplicity; Tracks", {HistType::kTH1F, {{300, -0.5, 299.5}}});
          registry.add("colInclusive/PVCFT0", "PV contributors with FIT; PV contributors; Tracks", {HistType::kTH1F, {{100, -0.5, 99.5}}});
          registry.add("colInclusive/vtxTracks", "Number of vertex tracks; Number of contributors; Collisions", {HistType::kTH1F, {{300, -0.5, 299.5}}});
          registry.add("colInclusive/PVTracks", "Number of PV tracks; Number of PV tracks; Tracks", {HistType::kTH1F, {{300, -0.5, 299.5}}});
    
          registry.add("colInclusive/FV0T0CPVtrk", "Number of vertex tracks; Number of contributors; Collisions", {HistType::kTH1F, {{300, -0.5, 299.5}}});
          registry.add("colInclusive/FV0T0Ctrk", "Number of PV tracks; Number of PV tracks; Tracks", {HistType::kTH1F, {{300, -0.5, 299.5}}});
            
          
          registry.add("colInclusive/FV0Amp", "FV0A Amplitude", {HistType::kTH1F, {{2000, 0., 2000.}}});
          registry.add("colInclusive/FT0Aamp", "#FT0A Amplitude; FT0A Amplitude", {HistType::kTH1F, {{2000, -0.5, 1999.5}}});
          registry.add("colInclusive/FT0Camp", "#FT0C Amplitude; FT0C Amplitude", {HistType::kTH1F, {{2000, -0.5, 1999.5}}});
          registry.add("colInclusive/hT0AC", "Time Correlation FT0; FT0A Time (ns) ; FT0C Time (ns)", {HistType::kTH2F, {{500, -50.0, 50.0},{500, -50.0, 50.0}}});
          registry.add("colInclusive/FT0ACCorr", "FT0 amp correlation; FT0A Amplitude; FT0C Amplitude", {HistType::kTH2F, {{2000, -4.5, 1995.5}, {2000, -4.5, 1995.5}}});
          
          registry.add("colInclusive/FDDAamp", "#FDDA Amplitude; FDDA Amplitude", {HistType::kTH1F, {{2000, -0.5, 1999.5}}});
          registry.add("colInclusive/FDDCamp", "#FDDC Amplitude; FDDC Amplitude", {HistType::kTH1F, {{2000, -0.5, 1999.5}}});
          registry.add("colInclusive/FT0DDAPVtrk", "Number of PV tracks; Number of PV tracks; Tracks", {HistType::kTH1F, {{300, -0.5, 299.5}}});
          
     registry.add("colInclusive/FT0DDAtrk", "Multiplicity of all tracks; Multiplicity; Tracks", {HistType::kTH1F, {{300, -0.5, 299.5}}});
     registry.add("colInclusive/hFDDAC", "Time Correlation FDD; FDDA time (ns); FDDC time (ns)", {HistType::kTH2F, {{500, -50.0, 50.0},{500, -50.0, 50.0}}});
     registry.add("colInclusive/FDDACCorr", "#FDD amp correlation; FDDA Amplitude; FDDC Amplitude", {HistType::kTH2F, {{2000, -4.5, 1995.5}, {2000, -4.5, 1995.5}}});
     registry.add("colInclusive/FT0DDCCorr", "Correlation FT0 vs FDD C side; FT0C Amplitude; FDDC Amplitude", {HistType::kTH2F, {{2000, -4.5, 1995.5}, {2000, -4.5, 1995.5}}});
     registry.add("colInclusive/FV0AFDDC", "Correlation FV0 vs FDD C side; FV0A Amplitude; FDDC Amplitude", {HistType::kTH2F, {{2000, -4.5, 1995.5}, {2000, -4.5, 1995.5}}});


      }
  }

//...............................................................................................................
void processMain(CC const& collision, BCs const& bct0s,TCs const& tracks,aod::FT0s const& ft0s, aod::FV0As const& fv0as,aod::FDDs const& fdds, aod::Zdcs& zdcs,aod::V0s const& v0s)
  {
    
    LOGF(debug, "<FITtest. Collision %d", collision.globalIndex());
      
      uint64_t bcnum = 0;
      registry.get<TH1>(HIST("collisions/Stat"))->Fill(0.);
  if (collision.has_foundBC()) {
        auto collbc = collision.foundBC_as<BCs>();
        bcnum = collbc.globalBC() % o2::constants::lhc::LHCMaxBunches;
      }
      
  if (!collision.has_foundBC()) return;
    
      registry.get<TH1>(HIST("collisions/Stat"))->Fill(1.);
      registry.get<TH1>(HIST("collisions/RelativeBC"))->Fill(bcnum, 1.);
      registry.get<TH1>(HIST("collisions/vtxTracks"))->Fill(collision.numContrib()); //vertex contributor
      
      // PV contributors
    int nPVcont = 0; int nCount = 0;
      for (auto const& trk : tracks) {
          if (trk.eta() > -0.5 && trk.eta() < 0.5) {
              if (trk.isPVContributor()) {nPVcont++;}
              nCount++;
            }
      } //track loop
    registry.get<TH1>(HIST("collisions/trkmultiplicity"))->Fill(nCount); //tracks
    registry.get<TH1>(HIST("collisions/PVTracks"))->Fill(nPVcont);//PVtracks
      
float totAmplitudeA = 0; float totAmplitudeC = 0; float totalAmplitudefv0 = 0; float totAmpFddA = 0; float totAmpFddC = 0;
  
if(collision.has_foundFT0()) {
    
    auto ft0 = collision.foundFT0();
    // side A
    for (size_t ind = 0; ind < ft0.channelA().size(); ind++) {
            registry.get<TH2>(HIST("FT0/FT0A"))->Fill((ft0.channelA())[ind], (ft0.amplitudeA())[ind]);
          }
          
    // side C
    for (size_t ind = 0; ind < ft0.channelC().size(); ind++) {
           registry.get<TH2>(HIST("FT0/FT0C"))->Fill((ft0.channelC())[ind], (ft0.amplitudeC())[ind]);
          }
      
    for (auto ampa : ft0.amplitudeA()) {
              totAmplitudeA += ampa;
          }
      
    for (auto ampc : ft0.amplitudeC()) {
              totAmplitudeC +=ampc;
          }
    
    registry.get<TH1>(HIST("collisions/Stat"))->Fill(2.);
    registry.get<TH1>(HIST("FT0/hT0A"))->Fill(ft0.timeA());
    registry.get<TH1>(HIST("FT0/hT0C"))->Fill(ft0.timeC());
    registry.get<TH2>(HIST("FT0/hT0AC"))->Fill(ft0.timeA(),ft0.timeC());
        
} else {
    if(!collision.has_foundFT0()) {
        totAmplitudeA=0;
        totAmplitudeC=0;
        }
    }//ends FT0
                    
      registry.get<TH1>(HIST("FT0/FT0Aamp"))->Fill(totAmplitudeA);
      registry.get<TH1>(HIST("FT0/FT0Camp"))->Fill(totAmplitudeC);
      registry.get<TH2>(HIST("FT0/FT0ACCorr"))->Fill(totAmplitudeA,totAmplitudeC);
      
      //FV0 information
  if (collision.has_foundFV0()) {
        
    auto fv0 = collision.foundFV0();
     registry.get<TH1>(HIST("collisions/Stat"))->Fill(3.);
     registry.get<TH1>(HIST("FV0/hV0A"))->Fill(fv0.time());
             
   for (size_t ind = 0; ind < fv0.channel().size(); ind++) {
     registry.get<TH2>(HIST("FV0/FV0A"))->Fill((fv0.channel())[ind], (fv0.amplitude())[ind]);
          }
          
   for (auto ampfv0a : fv0.amplitude()) {
               totalAmplitudefv0 +=ampfv0a;
           }
          
    } else {
        if(!collision.has_foundFV0()) {
            totalAmplitudefv0 = 0;
        }
    }
      registry.get<TH1>(HIST("FV0/FV0Amp"))->Fill(totalAmplitudefv0);
      
    
//FDD information
if (collision.has_foundFDD()) {
    auto fdd = collision.foundFDD();
        registry.get<TH1>(HIST("collisions/Stat"))->Fill(4.);
        registry.get<TH1>(HIST("FDD/hFDDA"))->Fill(fdd.timeA());
        registry.get<TH1>(HIST("FDD/hFDDC"))->Fill(fdd.timeC());
        registry.get<TH2>(HIST("FDD/hFDDAC"))->Fill(fdd.timeA(),fdd.timeC());
        
        // side A
        for (auto ind = 0; ind < 8; ind++) {
            registry.get<TH2>(HIST("FDD/FDDA"))->Fill(ind, (fdd.chargeA())[ind]);
        }
        
        // side C
        for (auto ind = 0; ind < 8; ind++) {
            registry.get<TH2>(HIST("FDD/FDDC"))->Fill(ind, (fdd.chargeC())[ind]);
        }
    
    for (auto ampfdd : fdd.chargeA()) {
        totAmpFddA += ampfdd;
    }
  
   for (auto ampfddc : fdd.chargeC()) {
            totAmpFddC +=ampfddc;
        }
} else {
    if (!collision.has_foundFDD()) {
        totAmpFddA=0;
        totAmpFddC=0;
    }
} //fdd
        registry.get<TH1>(HIST("FDD/FDDAamp"))->Fill(totAmpFddA);
        registry.get<TH1>(HIST("FDD/FDDCamp"))->Fill(totAmpFddC);
        registry.get<TH2>(HIST("FDD/FDDACCorr"))->Fill(totAmpFddA, totAmpFddC);
   
    //Correlation FV0 vs FT0
        registry.get<TH2>(HIST("collisions/FV0T0ACorr"))->Fill(totAmplitudeA, totalAmplitudefv0);
        registry.get<TH2>(HIST("collisions/FV0T0CCorr"))->Fill(totAmplitudeC, totalAmplitudefv0);
        
    //Correlation FDD vs FT0
        registry.get<TH2>(HIST("collisions/FT0DDACorr"))->Fill(totAmplitudeA, totAmpFddA);
        registry.get<TH2>(HIST("collisions/FT0DDCCorr"))->Fill(totAmplitudeC, totAmpFddC);
        registry.get<TH2>(HIST("collisions/FT0CFDDA"))->Fill(totAmplitudeC, totAmpFddA);
        registry.get<TH2>(HIST("collisions/FT0AFDDC"))->Fill(totAmplitudeA, totAmpFddC);
    //Correlation FDD vs FV0
        registry.get<TH2>(HIST("collisions/FV0AFDDA"))->Fill(totalAmplitudefv0, totAmpFddA);
        registry.get<TH2>(HIST("collisions/FV0AFDDC"))->Fill(totalAmplitudefv0, totAmpFddC);
      
  if(collision.has_foundFT0() && collision.has_foundFDD() && collision.has_foundFV0()){
        registry.get<TH1>(HIST("collisions/PVCFIT"))->Fill(nPVcont);
        registry.get<TH1>(HIST("collisions/GlobalBC"))->Fill(bcnum, 1.);
      }
  }
  
PROCESS_SWITCH(FITtest, processMain, "Process Main", true);
//...............................................................................................................
void processHadronic(CC const& collision, BCs const& bct0s,TCs const& tracks,aod::FT0s const& ft0s, aod::FV0As const& fv0as,aod::FDDs const& fdds, aod::Zdcs& zdcs,aod::V0s const& v0s)
{
     LOGF(debug, "<FITtest. Collision %d", collision.globalIndex());
              
        uint64_t bcnum = 0;
              
        registry.get<TH1>(HIST("collHadronic/StatHadronic"))->Fill(0.);
   if (collision.has_foundBC()) {
             auto collbc = collision.foundBC_as<BCs>();
             bcnum = collbc.globalBC() % o2::constants::lhc::LHCMaxBunches;
              }
              
   if (!collision.has_foundBC()) return;
              registry.get<TH1>(HIST("collHadronic/StatHadronic"))->Fill(1.);
   if (collision.numContrib()<10) return;
              
              registry.get<TH1>(HIST("collHadronic/StatHadronic"))->Fill(2.);
              registry.get<TH1>(HIST("collHadronic/RelBCHadronic"))->Fill(bcnum, 1.);
              registry.get<TH1>(HIST("collHadronic/vtxTrkHadronic"))->Fill(collision.numContrib());
              
              // PV contributors
        int nPVcont = 0; int nCont = 0;
              for (auto const& trk : tracks) {
                  if (trk.eta() > -0.5 && trk.eta() < 0.5) {
                        if(trk.pt() > 1){
                            if (trk.isPVContributor()) {nPVcont++;}
                            nCont++;
                        }
                    }
              }//tracks
       
        registry.get<TH1>(HIST("collHadronic/PVTrkHadronic"))->Fill(nPVcont);
        registry.get<TH1>(HIST("collHadronic/trkmultHadronic"))->Fill(nCont);
        
        float totAmplitudeA = 0; float totAmplitudeC = 0; float totalAmplitudefv0 = 0; float totAmpFddA = 0; float totAmpFddC = 0;
          
     if(collision.has_foundFT0()) {
            registry.get<TH1>(HIST("collHadronic/StatHadronic"))->Fill(3.);
            auto ft0 = collision.foundFT0();
            
            for (auto ampa : ft0.amplitudeA()) {
                      totAmplitudeA += ampa;
                  }
              
            for (auto ampc : ft0.amplitudeC()) {
                      totAmplitudeC +=ampc;
                  }
            
            registry.get<TH1>(HIST("collHadronic/FT0/hT0AHadronic"))->Fill(ft0.timeA());
            registry.get<TH1>(HIST("collHadronic/FT0/hT0CHadronic"))->Fill(ft0.timeC());
            registry.get<TH2>(HIST("collHadronic/FT0/hT0ACHadronic"))->Fill(ft0.timeA(),ft0.timeC());
                
        } else {
            if(!collision.has_foundFT0()) {
            totAmplitudeA=0;
            totAmplitudeC=0;
                }
            }//ends FT0 collsion
                            
              if(totAmplitudeA>0){registry.get<TH1>(HIST("collHadronic/FT0/FT0AampHadronic"))->Fill(nPVcont);}
              if(totAmplitudeC>0){registry.get<TH1>(HIST("collHadronic/FT0/FT0CampHadronic"))->Fill(nPVcont);}
              if((totAmplitudeA>0) || (totAmplitudeC>0)){registry.get<TH2>(HIST("collHadronic/FT0/FT0ACCorrHadronic"))->Fill(totAmplitudeA,totAmplitudeC);}
              
            //FV0 information
    if (collision.has_foundFV0()) {
                registry.get<TH1>(HIST("collHadronic/StatHadronic"))->Fill(4.);
               
                auto fv0 = collision.foundFV0();
                registry.get<TH1>(HIST("collHadronic/FV0/hV0AHadronic"))->Fill(fv0.time());
                
                for (auto ampfv0a : fv0.amplitude()) {
                       totalAmplitudefv0 +=ampfv0a;
                   }
                  
            } else {
                if(!collision.has_foundFV0()) {
                    totalAmplitudefv0 = 0;
                }
            }//FV0 collisions
    if(totalAmplitudefv0>0){registry.get<TH1>(HIST("collHadronic/FV0/FV0AmpHadronic"))->Fill(nPVcont);}
    
              
                       
        //FDD information
   if (collision.has_foundFDD()) {
            auto fdd = collision.foundFDD();
                registry.get<TH1>(HIST("collHadronic/StatHadronic"))->Fill(5.);
                registry.get<TH1>(HIST("collHadronic/FDD/hFDDAHadronic"))->Fill(fdd.timeA());
                registry.get<TH1>(HIST("collHadronic/FDD/hFDDCHadronic"))->Fill(fdd.timeC());
                registry.get<TH2>(HIST("collHadronic/FDD/hFDDACHadronic"))->Fill(fdd.timeA(),fdd.timeC());
                
            for (auto ampfdd : fdd.chargeA()) {
                totAmpFddA += ampfdd;
            }
          
           for (auto ampfddc : fdd.chargeC()) {
                    totAmpFddC +=ampfddc;
                }
        } else {
            if (!collision.has_foundFDD()) {
                totAmpFddA=0;
                totAmpFddC=0;
            }
        }//fdd
        
    if(totAmpFddA>0){registry.get<TH1>(HIST("collHadronic/FDD/FDDAampHadronic"))->Fill(nPVcont);}
    if(totAmpFddC>0){registry.get<TH1>(HIST("collHadronic/FDD/FDDCampHadronic"))->Fill(nPVcont);}
    if((totAmpFddA>0)||(totAmpFddC>0)){registry.get<TH2>(HIST("collHadronic/FDD/FDDACCorrHadronic"))->Fill(totAmpFddA, totAmpFddC);}
        
    //Correlation FV0 vs FT0
    if((totalAmplitudefv0>0)||(totAmplitudeA>0)){registry.get<TH2>(HIST("collHadronic/FV0T0ACorrHadronic"))->Fill(totAmplitudeA, totalAmplitudefv0);}
    if((totalAmplitudefv0>0)||(totAmplitudeC>0)){registry.get<TH2>(HIST("collHadronic/FV0T0CCorrHadronic"))->Fill(totAmplitudeC, totalAmplitudefv0);}
    
    //Correlation FDD vs FT0
    if((totAmpFddA>0)||(totAmplitudeA>0)){registry.get<TH2>(HIST("collHadronic/FT0DDACorrHadronic"))->Fill(totAmplitudeA, totAmpFddA);}
    if((totAmpFddC>0)||(totAmplitudeC>0)){registry.get<TH2>(HIST("collHadronic/FT0DDCCorrHadronic"))->Fill(totAmplitudeC, totAmpFddC);}
    if((totAmpFddC>0)||(totAmplitudeA>0)){registry.get<TH2>(HIST("collHadronic/FT0CFDDAHadronic"))->Fill(totAmplitudeC, totAmpFddA);}
    if((totAmpFddA>0)||(totAmplitudeC>0)){registry.get<TH2>(HIST("collHadronic/FT0AFDDCHadronic"))->Fill(totAmplitudeA, totAmpFddC);}
    
    //Correlation FDD vs FV0
    if((totAmpFddA>0)||(totalAmplitudefv0>0)){registry.get<TH2>(HIST("collHadronic/FV0AFDDAHadronic"))->Fill(totalAmplitudefv0, totAmpFddA);}
    if((totAmpFddC>0)||(totalAmplitudefv0>0)){registry.get<TH2>(HIST("collHadronic/FV0AFDDCHadronic"))->Fill(totalAmplitudefv0, totAmpFddC);}
            
   if(collision.has_foundFT0() && collision.has_foundFDD() && collision.has_foundFV0()){
        registry.get<TH1>(HIST("collHadronic/PVCFITHadronic"))->Fill(nPVcont);
        registry.get<TH1>(HIST("collHadronic/GlBCHadronic"))->Fill(bcnum, 1.);
              }
}
          
PROCESS_SWITCH(FITtest, processHadronic, "Process for hadroniclike events", true);
//...............................................................................................................
void processInclusiveA(CC const& collision, BCs const& bct0s,TCs const& tracks,aod::FT0s const& ft0s, aod::FV0As const& fv0as,aod::FDDs const& fdds, aod::Zdcs& zdcs,aod::V0s const& v0s)
    {
      
      LOGF(debug, "<FITtest. Collision %d", collision.globalIndex());
        
        uint64_t bcnum = 0;
        registry.get<TH1>(HIST("colInclusive/Stat"))->Fill(0.);
    
    if (collision.has_foundBC()) {
          auto collbc = collision.foundBC_as<BCs>();
          bcnum = collbc.globalBC() % o2::constants::lhc::LHCMaxBunches;
        }
        
    if (!collision.has_foundBC()) return;
        registry.get<TH1>(HIST("colInclusive/Stat"))->Fill(1.);
        
    if (collision.numContrib()>100) return;
        registry.get<TH1>(HIST("colInclusive/Stat"))->Fill(2.);
        
        registry.get<TH1>(HIST("colInclusive/RelativeBC"))->Fill(bcnum, 1.);
        registry.get<TH1>(HIST("colInclusive/vtxTracks"))->Fill(collision.numContrib());
        
        // PV contributors
    int nPVcont = 0; int ntrks =0;
        for (auto const& trk : tracks) {
            if (trk.eta() > -0.5 && trk.eta() < 0.5) {
                  if(trk.pt() < 10) {
                      if (trk.isPVContributor()) {nPVcont++;}
                      ntrks++;
                  }
              }
        }
    
        registry.get<TH1>(HIST("colInclusive/PVTracks"))->Fill(nPVcont);
        registry.get<TH1>(HIST("colInclusive/trkmultiplicity"))->Fill(ntrks);
    
  float totAmplitudeA = 0; float totAmplitudeC = 0; float totalAmplitudefv0 = 0; float totAmpFddA = 0; float totAmpFddC = 0;
  
    //FT0 information
  if(collision.has_foundFT0()) {
      registry.get<TH1>(HIST("colInclusive/Stat"))->Fill(3.);
        auto ft0 = collision.foundFT0();
      
      for (auto ampa : ft0.amplitudeA()) {
                totAmplitudeA += ampa;
            }
        
      for (auto ampc : ft0.amplitudeC()) {
                totAmplitudeC +=ampc;
            }
      registry.get<TH2>(HIST("colInclusive/hT0AC"))->Fill(ft0.timeA(),ft0.timeC());
          
  } else {
      if(!collision.has_foundFT0()) {
          totAmplitudeA=0;
          totAmplitudeC=0;
          }
      }//ends collsion
                      
        if(totAmplitudeC==0){registry.get<TH1>(HIST("colInclusive/FT0Aamp"))->Fill(totAmplitudeA);}
        if(totAmplitudeA==0){registry.get<TH1>(HIST("colInclusive/FT0Camp"))->Fill(totAmplitudeC);}
        if((totAmplitudeA<5) && ((totAmplitudeC>5)))registry.get<TH2>(HIST("colInclusive/FT0ACCorr"))->Fill(totAmplitudeA,totAmplitudeC);
        
    //FV0 information
    if(collision.has_foundFV0()) {
          registry.get<TH1>(HIST("colInclusive/Stat"))->Fill(4.);
         
          auto fv0 = collision.foundFV0();
            
          for (auto ampfv0a : fv0.amplitude()) {
                 totalAmplitudefv0 +=ampfv0a;
             }
            
       } else {
          if(!collision.has_foundFV0()) {
              totalAmplitudefv0 = 0;
          }
      }
     
      if((totalAmplitudefv0==0)&&(nPVcont<100)){registry.get<TH1>(HIST("colInclusive/FV0Amp"))->Fill(nPVcont);}
        
      if((totalAmplitudefv0==0) && (totAmplitudeA==0)){
            registry.get<TH1>(HIST("colInclusive/FV0T0CPVtrk"))->Fill(nPVcont);
            registry.get<TH1>(HIST("colInclusive/FV0T0Ctrk"))->Fill(ntrks);
        }
        
  //FDD information
  if (collision.has_foundFDD()) {
      auto fdd = collision.foundFDD();
          registry.get<TH1>(HIST("collisions/Stat"))->Fill(4.);
          registry.get<TH2>(HIST("colInclusive/hFDDAC"))->Fill(fdd.timeA(),fdd.timeC());
        
      for (auto ampfdd : fdd.chargeA()) {
          totAmpFddA += ampfdd;
      }
    
     for (auto ampfddc : fdd.chargeC()) {
              totAmpFddC +=ampfddc;
          }
  } else {
      if (!collision.has_foundFDD()) {
          totAmpFddA=0;
          totAmpFddC=0;
      }
  }//fdd
    
         if(totAmpFddA==0){registry.get<TH1>(HIST("colInclusive/FDDAamp"))->Fill(totAmpFddC);}
         if(totAmpFddC==0){registry.get<TH1>(HIST("colInclusive/FDDCamp"))->Fill(totAmpFddA);}
         if(totAmpFddA<5 && totAmpFddC>5){registry.get<TH2>(HIST("colInclusive/FDDACCorr"))->Fill(totAmpFddA, totAmpFddC);}
      
         if((totAmplitudeA==0) && (totAmpFddA==0)){registry.get<TH1>(HIST("colInclusive/FT0DDAPVtrk"))->Fill(nPVcont);}
         if((totAmplitudeA==0) && (totAmpFddA==0)){registry.get<TH1>(HIST("colInclusive/FT0DDAtrk"))->Fill(ntrks);}
         if((totAmplitudeA==0) && (totAmpFddA==0)){registry.get<TH2>(HIST("colInclusive/FT0DDCCorr"))->Fill(totAmplitudeC, totAmpFddC);}
           
    if((totalAmplitudefv0<5) && (totAmpFddA==0)){registry.get<TH2>(HIST("colInclusive/FV0AFDDC"))->Fill(totalAmplitudefv0, totAmpFddC);}
    
    if(collision.has_foundFT0() && collision.has_foundFDD() && collision.has_foundFV0()){
          if((totAmplitudeA==0) && (totAmpFddA==0)&&(totalAmplitudefv0==0)){registry.get<TH1>(HIST("colInclusive/PVCFT0"))->Fill(nPVcont);}
          if((totAmplitudeA==0) && (totAmpFddA==0)&&(totalAmplitudefv0==0)){registry.get<TH1>(HIST("colInclusive/GlobalBC"))->Fill(bcnum, 1.);}
        }
    }
    
    PROCESS_SWITCH(FITtest, processInclusiveA, "Process Inclusive veto A side", true);
    
/*void processInclusiveC(CC const& collision, BCs const& bct0s,TCs const& tracks,aod::FT0s const& ft0s, aod::FV0As const& fv0as,aod::FDDs const& fdds, aod::Zdcs& zdcs,aod::V0s const& v0s)
    {
      
      LOGF(debug, "<FITtest. Collision %d", collision.globalIndex());
        
        uint64_t bcnum = 0;
        //auto bc = collision.foundBC_as<BCs>();
        registry.get<TH1>(HIST("collisions/Stat"))->Fill(0.);
        if (collision.has_foundBC()) {
          auto collbc = collision.foundBC_as<BCs>();
          bcnum = collbc.globalBC() % o2::constants::lhc::LHCMaxBunches;
        }
        
        if (!collision.has_foundBC()) return;
        if (collision.numContrib()<10) return;
        registry.get<TH1>(HIST("collisions/Stat"))->Fill(1.);
        registry.get<TH1>(HIST("collisions/RelativeBC"))->Fill(bcnum, 1.);
        registry.get<TH1>(HIST("collisions/trkmultiplicity"))->Fill(tracks.size(), 1.);
        registry.get<TH1>(HIST("collisions/vtxTracks"))->Fill(collision.numContrib());
        
        // PV contributors
        int nPVcont = 0;
        for (auto const& trk : tracks) {
          if (trk.isPVContributor()) {
                  if(trk.pt() > 1){
                      if (trk.eta() > -0.5 && trk.eta() < 0.5) {nPVcont++;}
                  }
              }
        }
        registry.get<TH1>(HIST("collisions/PVTracks"))->Fill(nPVcont);
        
  float totAmplitudeA = 0; float totAmplitudeC = 0; float totalAmplitudefv0 = 0; float totAmpFddA = 0; float totAmpFddC = 0;
    
  if(collision.has_foundFT0()) {
      registry.get<TH1>(HIST("collisions/Stat"))->Fill(2.);
        auto ft0 = collision.foundFT0();
      // side A
      for (size_t ind = 0; ind < ft0.channelA().size(); ind++) {
              registry.get<TH2>(HIST("FT0A"))->Fill((ft0.channelA())[ind], (ft0.amplitudeA())[ind]);
            }
            
      // side C
      for (size_t ind = 0; ind < ft0.channelC().size(); ind++) {
             registry.get<TH2>(HIST("FT0C"))->Fill((ft0.channelC())[ind], (ft0.amplitudeC())[ind]);
            }
        
      for (auto ampa : ft0.amplitudeA()) {
                totAmplitudeA += ampa;
            }
        
      for (auto ampc : ft0.amplitudeC()) {
                totAmplitudeC +=ampc;
            }
      
      registry.get<TH1>(HIST("hT0A"))->Fill(ft0.timeA());
      registry.get<TH1>(HIST("hT0C"))->Fill(ft0.timeC());
      registry.get<TH2>(HIST("hT0AC"))->Fill(ft0.timeA(),ft0.timeC());
          
  } else {
      if(!collision.has_foundFT0()) {
          totAmplitudeA=0;
          totAmplitudeC=0;
          }
      }//ends collsion
                      
        if(totAmplitudeA>0){registry.get<TH1>(HIST("FT0Aamp"))->Fill(totAmplitudeA,nPVcont);}
        if(totAmplitudeC>0){registry.get<TH1>(HIST("FT0Camp"))->Fill(totAmplitudeC,nPVcont);}
        if((totAmplitudeA>0) && (totAmplitudeC>0))registry.get<TH2>(HIST("FT0ACCorr"))->Fill(totAmplitudeA,totAmplitudeC,nPVcont);
        
        //FV0 information
      if (collision.has_foundFV0()) {
          registry.get<TH1>(HIST("collisions/Stat"))->Fill(3.);
         
          auto fv0 = collision.foundFV0();
            
                registry.get<TH1>(HIST("hV0A"))->Fill(fv0.time());
               
          for (size_t ind = 0; ind < fv0.channel().size(); ind++) {
                    registry.get<TH2>(HIST("FV0A"))->Fill((fv0.channel())[ind], (fv0.amplitude())[ind]);
            }
            
          for (auto ampfv0a : fv0.amplitude()) {
                 totalAmplitudefv0 +=ampfv0a;
             }
            
      } else {
          if(!collision.has_foundFV0()) {
              totalAmplitudefv0 = 0;
          }
      }
        registry.get<TH1>(HIST("FV0Amp"))->Fill(totalAmplitudefv0,nPVcont);
        
     // if(collision.has_foundFT0() && collision.has_foundFV0()){
            registry.get<TH2>(HIST("FV0T0ACorr"))->Fill(totAmplitudeA, totalAmplitudefv0,nPVcont);
            registry.get<TH2>(HIST("FV0T0CCorr"))->Fill(totAmplitudeC, totalAmplitudefv0,nPVcont);
       // }
        
  //FDD information
  if (collision.has_foundFDD()) {
      auto fdd = collision.foundFDD();
      registry.get<TH1>(HIST("collisions/Stat"))->Fill(4.);
          registry.get<TH1>(HIST("hFDDA"))->Fill(fdd.timeA());
          registry.get<TH1>(HIST("hFDDC"))->Fill(fdd.timeC());
          registry.get<TH1>(HIST("hFDDAC"))->Fill(fdd.timeA()-fdd.timeC());
          
          // side A
          for (auto ind = 0; ind < 8; ind++) {
              registry.get<TH2>(HIST("FDDA"))->Fill(ind, (fdd.chargeA())[ind]);
          }
          
          // side C
          for (auto ind = 0; ind < 8; ind++) {
              registry.get<TH2>(HIST("FDDC"))->Fill(ind, (fdd.chargeC())[ind]);
          }
      
      for (auto ampfdd : fdd.chargeA()) {
          totAmpFddA += ampfdd;
      }
    
     for (auto ampfddc : fdd.chargeC()) {
              totAmpFddC +=ampfddc;
          }
      registry.get<TH2>(HIST("FDDACCorr"))->Fill(totAmpFddA, totAmpFddC,nPVcont);
      
  } else {
      if (!collision.has_foundFDD()) {
          totAmpFddA=0;
          totAmpFddC=0;
      }
  }//fdd

      //if(collision.has_foundFT0() && collision.has_foundFDD()){
            registry.get<TH2>(HIST("FT0DDACorr"))->Fill(totAmplitudeA, totAmpFddA,nPVcont);
            registry.get<TH2>(HIST("FT0DDCCorr"))->Fill(totAmplitudeC, totAmpFddC,nPVcont);
            registry.get<TH2>(HIST("FT0CFDDA"))->Fill(totAmplitudeC, totAmpFddA,nPVcont);
            registry.get<TH2>(HIST("FT0AFDDC"))->Fill(totAmplitudeA, totAmpFddC,nPVcont);
        //}
        
        //if(collision.has_foundFV0() && collision.has_foundFDD()){
            registry.get<TH2>(HIST("FV0AFDDA"))->Fill(totalAmplitudefv0, totAmpFddA,nPVcont);
            registry.get<TH2>(HIST("FV0AFDDC"))->Fill(totalAmplitudefv0, totAmpFddC,nPVcont);
        //}
        
if(collision.has_foundFT0() && collision.has_foundFDD() && collision.has_foundFV0()){
  registry.get<TH1>(HIST("collisions/PVCFT0"))->Fill(nPVcont);
  registry.get<TH1>(HIST("collisions/GlobalBC"))->Fill(bcnum, 1.);
  registry.get<TH2>(HIST("ACorr"))->Fill((totalAmplitudefv0+totAmplitudeA+totAmpFddA),(totAmpFddC+totalAmplitudefv0+totAmplitudeC));
  registry.get<TH1>(HIST("ACDiff"))->Fill((totalAmplitudefv0+totAmplitudeA+totAmpFddA)-(totAmpFddC+totalAmplitudefv0+totAmplitudeC));
        }
    }
    
    PROCESS_SWITCH(FITtest, processInclusiveC, "Process Inclusive veto C side", true);
    
void processUPC(CC const& collision, BCs const& bct0s,TCs const& tracks,aod::FT0s const& ft0s, aod::FV0As const& fv0as,aod::FDDs const& fdds, aod::Zdcs& zdcs,aod::V0s const& v0s)
    {
      
      LOGF(debug, "<FITtest. Collision %d", collision.globalIndex());
        
        uint64_t bcnum = 0;
        //auto bc = collision.foundBC_as<BCs>();
        registry.get<TH1>(HIST("collisions/Stat"))->Fill(0.);
        if (collision.has_foundBC()) {
          auto collbc = collision.foundBC_as<BCs>();
          bcnum = collbc.globalBC() % o2::constants::lhc::LHCMaxBunches;
        }
        
        if (!collision.has_foundBC()) return;
        if (collision.numContrib()<10) return;
        registry.get<TH1>(HIST("collisions/Stat"))->Fill(1.);
        registry.get<TH1>(HIST("collisions/RelativeBC"))->Fill(bcnum, 1.);
        registry.get<TH1>(HIST("collisions/trkmultiplicity"))->Fill(tracks.size(), 1.);
        registry.get<TH1>(HIST("collisions/vtxTracks"))->Fill(collision.numContrib());
        
        // PV contributors
        int nPVcont = 0;
        for (auto const& trk : tracks) {
          if (trk.isPVContributor()) {
                  if(trk.pt() > 1){
                      if (trk.eta() > -0.5 && trk.eta() < 0.5) {nPVcont++;}
                  }
              }
        }
        registry.get<TH1>(HIST("collisions/PVTracks"))->Fill(nPVcont);
        
  float totAmplitudeA = 0; float totAmplitudeC = 0; float totalAmplitudefv0 = 0; float totAmpFddA = 0; float totAmpFddC = 0;
    
  if(collision.has_foundFT0()) {
      registry.get<TH1>(HIST("collisions/Stat"))->Fill(2.);
        auto ft0 = collision.foundFT0();
      // side A
      for (size_t ind = 0; ind < ft0.channelA().size(); ind++) {
              registry.get<TH2>(HIST("FT0A"))->Fill((ft0.channelA())[ind], (ft0.amplitudeA())[ind]);
            }
            
      // side C
      for (size_t ind = 0; ind < ft0.channelC().size(); ind++) {
             registry.get<TH2>(HIST("FT0C"))->Fill((ft0.channelC())[ind], (ft0.amplitudeC())[ind]);
            }
        
      for (auto ampa : ft0.amplitudeA()) {
                totAmplitudeA += ampa;
            }
        
      for (auto ampc : ft0.amplitudeC()) {
                totAmplitudeC +=ampc;
            }
      
      registry.get<TH1>(HIST("hT0A"))->Fill(ft0.timeA());
      registry.get<TH1>(HIST("hT0C"))->Fill(ft0.timeC());
      registry.get<TH2>(HIST("hT0AC"))->Fill(ft0.timeA(),ft0.timeC());
          
  } else {
      if(!collision.has_foundFT0()) {
          totAmplitudeA=0;
          totAmplitudeC=0;
          }
      }//ends collsion
                      
        if(totAmplitudeA>0){registry.get<TH1>(HIST("FT0Aamp"))->Fill(totAmplitudeA,nPVcont);}
        if(totAmplitudeC>0){registry.get<TH1>(HIST("FT0Camp"))->Fill(totAmplitudeC,nPVcont);}
        if((totAmplitudeA>0) && (totAmplitudeC>0))registry.get<TH2>(HIST("FT0ACCorr"))->Fill(totAmplitudeA,totAmplitudeC,nPVcont);
        
        //FV0 information
      if (collision.has_foundFV0()) {
          registry.get<TH1>(HIST("collisions/Stat"))->Fill(3.);
         
          auto fv0 = collision.foundFV0();
            
                registry.get<TH1>(HIST("hV0A"))->Fill(fv0.time());
               
          for (size_t ind = 0; ind < fv0.channel().size(); ind++) {
                    registry.get<TH2>(HIST("FV0A"))->Fill((fv0.channel())[ind], (fv0.amplitude())[ind]);
            }
            
          for (auto ampfv0a : fv0.amplitude()) {
                 totalAmplitudefv0 +=ampfv0a;
             }
            
      } else {
          if(!collision.has_foundFV0()) {
              totalAmplitudefv0 = 0;
          }
      }
        registry.get<TH1>(HIST("FV0Amp"))->Fill(totalAmplitudefv0,nPVcont);
        
     // if(collision.has_foundFT0() && collision.has_foundFV0()){
            registry.get<TH2>(HIST("FV0T0ACorr"))->Fill(totAmplitudeA, totalAmplitudefv0,nPVcont);
            registry.get<TH2>(HIST("FV0T0CCorr"))->Fill(totAmplitudeC, totalAmplitudefv0,nPVcont);
       // }
        
  //FDD information
  if (collision.has_foundFDD()) {
      auto fdd = collision.foundFDD();
      registry.get<TH1>(HIST("collisions/Stat"))->Fill(4.);
          registry.get<TH1>(HIST("hFDDA"))->Fill(fdd.timeA());
          registry.get<TH1>(HIST("hFDDC"))->Fill(fdd.timeC());
          registry.get<TH1>(HIST("hFDDAC"))->Fill(fdd.timeA()-fdd.timeC());
          
          // side A
          for (auto ind = 0; ind < 8; ind++) {
              registry.get<TH2>(HIST("FDDA"))->Fill(ind, (fdd.chargeA())[ind]);
          }
          
          // side C
          for (auto ind = 0; ind < 8; ind++) {
              registry.get<TH2>(HIST("FDDC"))->Fill(ind, (fdd.chargeC())[ind]);
          }
      
      for (auto ampfdd : fdd.chargeA()) {
          totAmpFddA += ampfdd;
      }
    
     for (auto ampfddc : fdd.chargeC()) {
              totAmpFddC +=ampfddc;
          }
      registry.get<TH2>(HIST("FDDACCorr"))->Fill(totAmpFddA, totAmpFddC,nPVcont);
      
  } else {
      if (!collision.has_foundFDD()) {
          totAmpFddA=0;
          totAmpFddC=0;
      }
  }//fdd

      //if(collision.has_foundFT0() && collision.has_foundFDD()){
            registry.get<TH2>(HIST("FT0DDACorr"))->Fill(totAmplitudeA, totAmpFddA,nPVcont);
            registry.get<TH2>(HIST("FT0DDCCorr"))->Fill(totAmplitudeC, totAmpFddC,nPVcont);
            registry.get<TH2>(HIST("FT0CFDDA"))->Fill(totAmplitudeC, totAmpFddA,nPVcont);
            registry.get<TH2>(HIST("FT0AFDDC"))->Fill(totAmplitudeA, totAmpFddC,nPVcont);
        //}
        
        //if(collision.has_foundFV0() && collision.has_foundFDD()){
            registry.get<TH2>(HIST("FV0AFDDA"))->Fill(totalAmplitudefv0, totAmpFddA,nPVcont);
            registry.get<TH2>(HIST("FV0AFDDC"))->Fill(totalAmplitudefv0, totAmpFddC,nPVcont);
        //}
        
      if(collision.has_foundFT0() && collision.has_foundFDD() && collision.has_foundFV0()){
  registry.get<TH1>(HIST("collisions/PVCFT0"))->Fill(nPVcont);
  registry.get<TH1>(HIST("collisions/GlobalBC"))->Fill(bcnum, 1.);
  registry.get<TH2>(HIST("ACorr"))->Fill((totalAmplitudefv0+totAmplitudeA+totAmpFddA),(totAmpFddC+totalAmplitudefv0+totAmplitudeC));
  registry.get<TH1>(HIST("ACDiff"))->Fill((totalAmplitudefv0+totAmplitudeA+totAmpFddA)-(totAmpFddC+totalAmplitudefv0+totAmplitudeC));
        }
    }
    
PROCESS_SWITCH(FITtest, processUPC, "Process exclusiveUPC veto A and C sides", true);*/
        
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<FITtest>(cfgc, TaskName{"fittest"}),
  };
}
