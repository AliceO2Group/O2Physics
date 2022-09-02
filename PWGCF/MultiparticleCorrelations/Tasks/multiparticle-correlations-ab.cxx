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
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
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

// *) Data members:
#include "PWGCF/MultiparticleCorrelations/Core/MuPa-DataMembers.h"

// *) Member functions:
#include "PWGCF/MultiparticleCorrelations/Core/MuPa-MemberFunctions.h"

// *) Configurables (cuts):
#include "PWGCF/MultiparticleCorrelations/Core/MuPa-Configurables.h"

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

    // *) Configure the task with setters and getters:
    //TH1D *phiWeights = task->GetHistogramWithWeights(Form("%s/%s/weights.root",directoryWeights.Data(),runNumber.Data()),"phi"); // original line
    TH1D* phiWeights = GetHistogramWithWeights("/alice/cern.ch/user/a/abilandz/weights.root", "phi"); // both relative and abs path shell be fine
    SetWeightsHist(phiWeights, "phi");

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
    BookResultsHistograms();

    // *) Trick to avoid name clashes, part 2:
    TH1::AddDirectory(oldHistAddStatus);

  } // void init(InitContext const&)

  // -------------------------------------------

  // *) Process the data:
  void process(aod::Collision const& collision, aod::Tracks const& tracks) // called once per collision found in the time frame
  {
    // *) Fill event histograms for reconstructed data before event cuts:
    FillEventHistograms(collision, tracks, eRec, eBefore);

    // *) Event cuts:
    if (!EventCuts(collision)) {
      return;
    }

    // *) Fill event histograms for reconstructed data after event cuts:
    FillEventHistograms(collision, tracks, eRec, eAfter);

    // *) Main loop over particles:
    Double_t dPhi = 0.; //, dPt = 0., dEta = 0.;
    // Double_t wPhi = 1., wPt = 1., wEta = 1.;
    Double_t wToPowerP = 1.; // final particle weight raised to power p
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
      // dPt  = track.pt();
      // dEta = track.eta();
      for (Int_t h = 0; h < gMaxHarmonic * gMaxCorrelator + 1; h++) {
        for (Int_t wp = 0; wp < gMaxCorrelator + 1; wp++) { // weight power
          // if (fUseWeights[0]||fUseWeights[1]||fUseWeights[2]) {
          //   wToPowerP = pow(wPhi*wPt*wEta,wp);
          // }
          qv_a.fQvector[h][wp] += TComplex(wToPowerP * TMath::Cos(h * dPhi), wToPowerP * TMath::Sin(h * dPhi));
        } // for(Int_t wp=0;wp<gMaxCorrelator+1;wp++)
      }   // for(Int_t h=0;h<gMaxHarmonic*gMaxCorrelator+1;h++)

      fResultsHist->Fill(pw_a.fWeightsHist[wPHI]->GetBinContent(pw_a.fWeightsHist[wPHI]->FindBin(track.phi()))); // TBI 20220713 meaningless, only temporarily here to check if this is feasible

    } // for (auto& track : tracks)

    // *) Calculate multiparticle correlations (standard, isotropic, same harmonic):
    if (fCalculateCorrelations) {
      CalculateCorrelations();
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
