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

// O2:
#include <CCDB/BasicCCDBManager.h>
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
//#include "Common/DataModel/TrackSelectionTables.h" TBI 20231008 needed for aod::TracksDCA
using namespace o2;
using namespace o2::framework;

// ROOT:
#include "TList.h"
#include "TSystem.h"
#include "TFile.h"
#include "TH1D.h"
#include "TGrid.h"
#include "Riostream.h"
#include "TRandom3.h"
#include <TComplex.h>
using namespace std;

// *) Enums:
#include "PWGCF/MultiparticleCorrelations/Core/MuPa-Enums.h"

// *) Global constants:
#include "PWGCF/MultiparticleCorrelations/Core/MuPa-GlobalConstants.h"

// *) Main task:
struct MultiparticleCorrelationsAB // this name is used in lower-case format to name the TDirectoryFile in AnalysisResults.root
{

  // *) CCDB:
  Service<ccdb::BasicCCDBManager> ccdb;

// *) Configurables (cuts):
#include "PWGCF/MultiparticleCorrelations/Core/MuPa-Configurables.h"

// *) Data members:
#include "PWGCF/MultiparticleCorrelations/Core/MuPa-DataMembers.h"

// *) Member functions:
#include "PWGCF/MultiparticleCorrelations/Core/MuPa-MemberFunctions.h"

  // -------------------------------------------

  // *) Initialize and book all objects:
  void init(InitContext const&)
  {
    // *) Trick to avoid name clashes, part 1;
    // *) Book base list;
    // *) Default configuration, booking, binning and cuts;
    // *) Configure the task with setters and getters;
    // *) Book all remaining objects;
    // ...
    // *) Trick to avoid name clashes, part 2;

    // *) Trick to avoid name clashes, part 1:
    Bool_t oldHistAddStatus = TH1::AddDirectoryStatus();
    TH1::AddDirectory(kFALSE);

    // *) Book base list:
    BookBaseList();

    // *) Default configuration, booking, binning and cuts:
    DefaultConfiguration();
    DefaultBooking();
    DefaultBinning();
    DefaultCuts(); // Remark: has to be called after DefaultBinning(), since some default cuts are defined through default binning, to ease bookeeping

    // *) Particle weights:
    if (pw_a.fUseWeights[wPHI]) {
      TH1D* phiWeights = GetHistogramWithWeights(fFileWithWeights.Data(), fRunNumber.Data(), "phi");
      SetWeightsHist(phiWeights, "phi");
    }
    if (pw_a.fUseWeights[wPT]) {
      TH1D* ptWeights = GetHistogramWithWeights(fFileWithWeights.Data(), fRunNumber.Data(), "pt");
      SetWeightsHist(ptWeights, "pt");
    }
    if (pw_a.fUseWeights[wETA]) {
      TH1D* etaWeights = GetHistogramWithWeights(fFileWithWeights.Data(), fRunNumber.Data(), "eta");
      SetWeightsHist(etaWeights, "eta");
    }

    // *) Book random generator:
    delete gRandom;
    gRandom = new TRandom3(fRandomSeed); // if uiSeed is 0, the seed is determined uniquely in space and time via TUUID

    // *) Book all remaining objects;
    BookAndNestAllLists();
    BookEventHistograms();
    BookParticleHistograms();
    BookQvectorHistograms();
    BookCorrelationsHistograms();
    BookWeightsHistograms();
    BookNestedLoopsHistograms();
    BookTest0Histograms();
    BookResultsHistograms();

    // *) Trick to avoid name clashes, part 2:
    TH1::AddDirectory(oldHistAddStatus);

  } // void init(InitContext const&)

  // -------------------------------------------

  // *) Process the data:
  //  void process(aod::Collision const& collision, soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA> const& tracks) // called once per collision found in the time frame
  void process(aod::Collision const& collision, aod::Tracks const& tracks) // called once per collision found in the time frame
  {

    // *) TBI 20231008 Temporarily here: If I reached max number of events, ignore the remaining collisions.
    //    But what I really need here is a graceful exit from subsequent processing (which will also dump the output file, etc.)
    if (ceh_a.fEventHistograms[eNumberOfEvents][eRec][eAfter]->GetBinContent(1) >= ceh_a.fEventCuts[eNumberOfEvents][eMax]) {
      return;
    }

    // *) Fill event histograms for reconstructed data before event cuts:
    FillEventHistograms(collision, tracks, eRec, eBefore);

    // *) Event cuts:
    if (!EventCuts(collision, tracks)) {
      return;
    }

    // *) Fill event histograms for reconstructed data after event cuts:
    FillEventHistograms(collision, tracks, eRec, eAfter);

    // *) Main loop over particles:
    Double_t dPhi = 0., wPhi = 1.; // azimuthal angle and corresponding phi weight
    Double_t dPt = 0., wPt = 1.;   // transverse momentum and corresponding pT weight
    Double_t dEta = 0., wEta = 1.; // pseudorapidity and corresponding eta weight
    Double_t wToPowerP = 1.;       // weight raised to power p
    fSelectedTracks = 0;           // reset number of selected tracks
    for (auto& track : tracks) {

      // *) Fill particle histograms for reconstructed data before particle cuts:
      FillParticleHistograms(track, eRec, eBefore);

      // *) Particle cuts:
      if (!ParticleCuts(track)) {
        continue;
      }

      // *) Fill particle histograms for reconstructed data after particle cuts:
      FillParticleHistograms(track, eRec, eAfter);

      // *) Fill Q-vectors:
      dPhi = track.phi();
      dPt = track.pt();
      dEta = track.eta();
      // Particle weights:
      if (pw_a.fUseWeights[wPHI]) {
        wPhi = Weight(dPhi, "phi"); // corresponding phi weight
        if (!(wPhi > 0.)) {
          LOGF(error, "\033[1;33m%s wPhi is not positive, skipping this particle for the time being...\033[0m", __PRETTY_FUNCTION__);
          LOGF(error, "dPhi = %f\nwPhi = %f", dPhi, wPhi);
          continue;
        }
      } // if(pw_a.fUseWeights[wPHI])
      if (pw_a.fUseWeights[wPT]) {
        wPt = Weight(dPt, "pt"); // corresponding pt weight
        if (!(wPt > 0.)) {
          LOGF(error, "\033[1;33m%s wPt is not positive, skipping this particle for the time being...\033[0m", __PRETTY_FUNCTION__);
          LOGF(error, "dPt = %f\nwPt = %f", dPt, wPt);
          continue;
        }
      } // if(pw_a.fUseWeights[wPT])
      if (pw_a.fUseWeights[wETA]) {
        wEta = Weight(dEta, "eta"); // corresponding eta weight
        if (!(wEta > 0.)) {
          LOGF(error, "\033[1;33m%s wEta is not positive, skipping this particle for the time being...\033[0m", __PRETTY_FUNCTION__);
          LOGF(error, "dEta = %f\nwEta = %f", dEta, wEta);
          continue;
        }
      } // if(pw_a.fUseWeights[wETA])

      for (Int_t h = 0; h < gMaxHarmonic * gMaxCorrelator + 1; h++) {
        for (Int_t wp = 0; wp < gMaxCorrelator + 1; wp++) { // weight power
          if (pw_a.fUseWeights[wPHI] || pw_a.fUseWeights[wPT] || pw_a.fUseWeights[wETA]) {
            wToPowerP = pow(wPhi * wPt * wEta, wp);
          }
          qv_a.fQvector[h][wp] += TComplex(wToPowerP * TMath::Cos(h * dPhi), wToPowerP * TMath::Sin(h * dPhi));
        } // for(Int_t wp=0;wp<gMaxCorrelator+1;wp++)
      }   // for(Int_t h=0;h<gMaxHarmonic*gMaxCorrelator+1;h++)

      // *) Nested loops containers:
      if (fCalculateNestedLoops || fCalculateCustomNestedLoop) {
        if (nl_a.ftaNestedLoops[0]) {
          nl_a.ftaNestedLoops[0]->AddAt(dPhi, fSelectedTracks);
        } // remember that the 2nd argument here must start from 0
        if (nl_a.ftaNestedLoops[1]) {
          nl_a.ftaNestedLoops[1]->AddAt(wPhi * wPt * wEta, fSelectedTracks);
        } // remember that the 2nd argument here must start from 0
      }   // if(fCalculateNestedLoops||fCalculateCustomNestedLoop)

      // *) Counter of selected tracks in the current event:
      fSelectedTracks++;
      if (fSelectedTracks >= cSelectedTracks_max) {
        break;
      }

    } // for (auto& track : tracks)

    // *) Fill remaining event histograms for reconstructed data after event (and particle) cuts:
    ceh_a.fEventHistograms[eSelectedTracks][eRec][eAfter]->Fill(fSelectedTracks);

    // *) Remaining event cuts:
    if ((fSelectedTracks < ceh_a.fEventCuts[eSelectedTracks][eMin]) || (fSelectedTracks > ceh_a.fEventCuts[eSelectedTracks][eMax])) {
      return;
    }

    // *) Calculate multiparticle correlations (standard, isotropic, same harmonic):
    if (fCalculateCorrelations) {
      CalculateCorrelations();
    }

    // *) Calculate nested loops:
    if (fCalculateNestedLoops) {
      CalculateNestedLoops();

      // TBI 20220823 this shall be called after all events are processed, only temporarily here called for each event:
      ComparisonNestedLoopsVsCorrelations();
    }

    // *) Reset event-by-event objects:
    ResetEventByEventQuantities();

  } // void process(...)

}; // struct MultiparticleCorrelationsAB

// -------------------------------------------

// *) ...
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<MultiparticleCorrelationsAB>(cfgc),
  };
} // WorkflowSpec...
