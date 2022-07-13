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
using namespace std;

// *) Global constants:
#include "PWGCF/MultiparticleCorrelations/Core/MuPa-GlobalConstants.h"

// *) Enums:
#include "PWGCF/MultiparticleCorrelations/Core/MuPa-Enums.h"

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
    // *) Configure the task with setters and getters;
    // *) Book random generator;
    // *) Book all remaining objects;
    // ...
    // *) Trick to avoid name clashes, part 2;

    // *) Trick to avoid name clashes, part 1:
    Bool_t oldHistAddStatus = TH1::AddDirectoryStatus();
    TH1::AddDirectory(kFALSE);

    // *) Book base list:
    BookBaseList();

    // *) Configure the task with setters and getters:
    //TH1D *phiWeights = task->GetHistogramWithWeights(Form("%s/%s/weights.root",directoryWeights.Data(),runNumber.Data()),"phi"); // original line
    TH1D* phiWeights = GetHistogramWithWeights("/alice/cern.ch/user/a/abilandz/weights.root", "phi"); // both relative and abs path shell be fine
    SetWeightsHist(phiWeights, "phi");

    // *) Book random generator:
    delete gRandom;
    gRandom = new TRandom3(fRandomSeed); // if uiSeed is 0, the seed is determined uniquely in space and time via TUUID

    // *) Book all remaining objects;
    BookAndNestAllLists();
    BookControlEventHistograms();
    BookWeightsHistograms();
    BookResultsHistograms();

    // *) Trick to avoid name clashes, part 2:
    TH1::AddDirectory(oldHistAddStatus);

  } // void init(InitContext const&)

  // -------------------------------------------

  // *) Process the data:
  void process(aod::Collision const& collision, aod::Tracks const& tracks) // called once per collision found in the time frame
  {
    // *) Event cuts:
    ceh_a.fMultiplicityHist[BEFORE]->Fill(tracks.size()); // TBI 20220713 -> member function
    if (!EventCuts(collision)) {
      return;
    }
    ceh_a.fMultiplicityHist[AFTER]->Fill(tracks.size()); // TBI 20220713 -> member function

    // *) Main loop over particles:
    for (auto& track : tracks) {

      if (!ParticleCuts(track)) {
        continue;
      }
      fResultsHist->Fill(pw_a.fWeightsHist[wPHI]->GetBinContent(pw_a.fWeightsHist[wPHI]->FindBin(track.phi()))); // TBI 20220713 meaningless, only temporarily here to check if this is feasible

    } // for (auto& track : tracks)

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
