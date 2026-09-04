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

/// \file multiparticleCorrelationsMei.cxx
/// \brief Multiparticle correlation in O2 Framework
/// \author yuanjun.mei@cern.ch

#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h" // needed for aod::TracksDCA table

#include <CCDB/BasicCCDBManager.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/DataTypes.h>
#include <Framework/runDataProcessing.h>

#include <TGrid.h>
#include <TH1D.h>
#include <TSystem.h>

#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::constants;
using namespace std;

// Definitions of join tables for Run 3 analysis:
using EventSelection = soa::Join<aod::EvSels, aod::Mults, aod::CentFT0Cs, aod::CentFT0Ms, aod::CentFV0As, aod::CentNTPVs>;
using CollisionRec = soa::Join<aod::Collisions, EventSelection>::iterator; // use in json "isMC": "true" for "event-selection-task"
using CollisionRecSim = soa::Join<aod::Collisions, aod::McCollisionLabels, EventSelection>::iterator;
using CollisionSim = aod::McCollision;
using TracksRec = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection>;
using TrackRec = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection>::iterator;
using TracksRecSim = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::McTrackLabels>; // + use in json "isMC" : "true"
using TrackRecSim = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::McTrackLabels>::iterator;
using TracksSim = aod::McParticles;
using TrackSim = aod::McParticles::iterator;

// *) Define enums:
enum ECentralityEstimator {
  EFT0C = 0,
  EFT0M,
  EFV0A,
  ENTPV
};

enum EMultiplicityTables {
  EMultTPC = 0,
  EMultFV0M,
  EMultFT0C,
  EMultFT0M,
  EMultNTracksPV
};

enum ERecSim {
  ERec = 0,
  ESim,
  ERecAndSim
};

enum ECuts {
  ENoCuts = 0,
  EWithCuts
};

enum EProcess {
  EProcessRec = 0, // Run 3, only reconstructed
  EProcessRecSim,  // Run 3, both reconstructed and simulated
  EProcessSim,     // Run 3, only simulated
  EProcess_N
};

enum EParticlEHistograms {
  EParticlEHistogramsList = 0,
  EHistPt,
  EHistPhi,
  EHistEta,
  EParticlEHistograms_N
};

enum EEventHistograms {
  EHistCentrality = 0,
  EHistMultiplicity,
  EHistVertexX,
  EHistVertexY,
  EHistVertexZ,
  EHistImpactParameter,
  EHistReferencEMultiplicity,
  EEventHistograms_N
};

enum EExternalHistograms {
  EHistWeights = 0,
  EExternalHistograms_N
};

// *) Main task:
struct MultiparticleCorrelationsMei // this name is used in lower-case format to name the TDirectoryFile in AnalysisResults.root
{
  // *) Base TList to hold all output objects:
  TString sBaseListName = "Default list name"; // yes, I declare it separately, because I need it also later in BailOut() function
  OutputObj<TList> fBaseList{sBaseListName.Data(),
                             OutputObjHandlingPolicy::AnalysisObject,
                             OutputObjSourceType::OutputObjSource};

  // *) CCDB:
  Service<ccdb::BasicCCDBManager> ccdb; // support for offline callibration data base, not needed for the time being...

  // *) Define configurables:
  Configurable<int> centralityEstimator{"centralityEstimator", 0, "centrality estimator: 0=FT0C, 1=FT0M, 2=FV0A, 3=NTPV"};
  Configurable<int> multiplicityTables{"multiplicityTables", 0, "multiplicity tables: 0=multTPC, 1=multFV0M, 2=multFT0C, 3=multFT0M, 4=multNTracksPV"};
  Configurable<bool> cfDryRun{"cfDryRun", false, "book all histos and run without filling and calculating anything"};

  // external root files
  Configurable<std::string> cfFileWithWeights{"cfFileWithWeights", "/scratch3/go97loy/O2challenge/C5/weights.root", "path to external ROOT file which holds all particle weights"};

  // binnings
  Configurable<std::vector<float>> cfPtBins{"cfPtBins", {1000, 0., 100.}, "nPtBins, ptMin, ptMax"};            // example for an array
  Configurable<std::vector<float>> cfPhiBins{"cfPhiBins", {100, 0., math::TwoPI}, "nPhiBins, phiMin, phiMax"}; // example for an array
  Configurable<std::vector<float>> cfEtaBins{"cfEtaBins", {100, 0., 10.}, "nEtaBins, etaMin, etaMax"};         // example for an array

  Configurable<std::vector<float>> cfMultBinsRec{"cfMultBinsRec", {200, 0., 10000.}, "nMultBins, multMin, multMax"};
  Configurable<std::vector<float>> cfMultBinsRef{"cfMultBinsRef", {200, 0., 10000.}, "nMultBins, multMin, multMax"};
  Configurable<std::vector<float>> cfMultBinsSim{"cfMultBinsSim", {100, 0., 1000.}, "nMultBins, multMin, multMax"};
  Configurable<std::vector<float>> cfVxBins{"cfVxBins", {100, -10., 10.}, "nVxBins, vxMin, vxMax"};
  Configurable<std::vector<float>> cfVyBins{"cfVyBins", {100, -10., 10.}, "nVyBins, vyMin, vyMax"};
  Configurable<std::vector<float>> cfVzBins{"cfVzBins", {200, -20., 20.}, "nVzBins, vzMin, vzMax"};
  Configurable<std::vector<float>> cfCentBins{"cfCentBins", {100, 0., 100.}, "nCentBins, centMin, centMax"};
  Configurable<std::vector<float>> cfIpBins{"cfIpBins", {100, 0., 20.}, "nIPBins, ipMin, ipMax"};

  // Cuts
  Configurable<bool> cfVertexZSwitch{"cfVertexZSwitch", false, "switch to apply vertex z position cut"};
  Configurable<std::vector<float>> cfVertexZ{"cfVertexZ", {-10, 10.}, "vertex z position range: {min, max}[cm], with convention: min <= Vz < max"};
  Configurable<bool> cfPtSwitch{"cfPtSwitch", false, "switch to apply pt cut"};
  Configurable<std::vector<float>> cfPt{"cfPt", {0.2, 5.}, "pt range: {min, max}, with convention: min <= Vz < max"};

  // misc
  Configurable<double> sigmaInel{"sigmaInel", 7.71, "inelastic cross section in mb"};
  Configurable<bool> qualityAssurance{"qualityAssurance", false, "quality assurance"};

  // *) Define and initialize all data members to be called in the main process* functions:
  // **) Task configuration:
  struct TaskConfiguration {
    bool fProcess[EProcess_N] = {false}; // Set what to process. See enum EProcess for full description. Set via implicit variables within a PROCESS_SWITCH clause.
    bool fDryRun = false;                // book all histos and run without filling and calculating anything
  } tc;                                  // you have to prepend "tc." for all objects name in this group later in the code

  // **) Particle histograms:
  struct ParticlEHistograms {
    TList* fParticlEHistogramsList = NULL; //!<! list to hold all control particle histograms
    TH1F* fParticlEHistograms[EParticlEHistograms_N][2][2] = {{{NULL}}};
    TH1F* fHistPt[2] = {NULL}; // pt distribution of a particle [ 0 = rec, 1 = sim ]
    TH1F* fHistPhi[2] = {NULL};
  } pc; // you have to prepend "pc." for all objects name in this group later in the code

  // *) Event histograms:
  struct EventHistograms {
    TList* fEventHistogramsList = NULL;                            //!<! list to hold all event-level histograms
    TH1F* fEventHistograms[EEventHistograms_N][2][2] = {{{NULL}}}; //! [ type - see enum EEventHistograms ][reco,sim][before, after event cuts]
  } ec;                                                            // prepend "ec." for event counters

  // *) External histograms:
  struct ExternalHistograms {
    TList* fExternalHistogramsList = NULL;
    TH1D* fhistWeights = NULL;
  } ex;

  // *) Quality assurance histograms:
  struct QualityAssurance {
    TList* fQualityAssuranceList = NULL; //!<! list to hold all qualityAssurance histograms
    TH2F* fHistCentralityRecSim = NULL;
  } qa; // prepend "qa." for qa histograms

  // functions
  TObject* getObjectFromList(TList* list, const char* objectName)
  {
    // Get TObject pointer from TList, even if it's in some nested TList. Foreseen
    // to be used to fetch histograms or profiles from files directly.
    // Some ideas taken from TCollection::ls()
    // If you have added histograms directly to files (without TList's), then you can fetch them directly with
    // file->Get("hist-name").

    // Usage: TH1D *hist = (TH1D*) getObjectFromList("some-valid-TList-pointer","some-object-name");

    // Insanity checks:
    if (!list) {
      LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
    }
    if (!objectName) {
      LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
    }
    if (0 == list->GetEntries()) {
      return NULL;
    }

    // The object is in the current base list:
    TObject* objectFinal = list->FindObject(objectName); // the final object I am after
    if (objectFinal) {
      return objectFinal;
    }

    // Otherwise, search for the object recursively in the nested lists:
    TObject* objectIter; // iterator object in the loop below
    TIter next(list);
    while ((objectIter = next())) // double round braces are to silence the warnings
    {
      if (TString(objectIter->ClassName()).EqualTo("TList")) {
        objectFinal = getObjectFromList(reinterpret_cast<TList*>(objectIter), objectName);
        if (objectFinal)
          return objectFinal;
      }
    } // while(objectIter = next())

    return NULL;
  } // TObject* getObjectFromList(TList *list, char *objectName)

  TH1D* getHistogramWithWeights(const char* filePath, const char* runNumber)
  {
    // *) Return value:
    TH1D* hist = NULL;
    TList* baseList = NULL;     // base top-level list in the TFile, e.g. named "ccdb_object"
    TList* listWithRuns = NULL; // nested list with run-wise TList's holding run-specific weights

    // *) Determine from filePath if the file is on a local machine, or in home dir AliEn, or in CCDB:
    //    Algorithm: If filePath begins with "/alice/cern.ch/" then it's in the home dir AliEn;
    //                     If filePath begins with "/alice-ccdb.cern.ch/" then it's in CCDB. Therefore, files in AliEn and CCDB must be specified with abs path;
    //                     for local files both abs and relative paths are just fine.
    bool bFileIsInAliEn = false;
    bool bFileIsInCCDB = false;

    std::string path(cfFileWithWeights);

    if (path.starts_with("/alice/cern.ch/")) {
      bFileIsInAliEn = true;
    } else if (path.starts_with("/alice-ccdb.cern.ch/")) {
      bFileIsInCCDB = true;
    }

    if (bFileIsInAliEn) {
      // File you want to access is in your home dir in AliEn:
      const TGrid* alien = TGrid::Connect("alien", gSystem->Getenv("USER"), "", ""); // do not forget to add #include <TGrid.h> to the preamble of your analysis task
      if (!alien) {
        LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
      }
      TFile* weightsFile = TFile::Open(Form("alien://%s", filePath), "READ"); // yes, ROOT can open a file transparently, even if it's sitting in AliEn, with this specific syntax
      if (!weightsFile) {
        LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
      }
      weightsFile->GetObject("ccdb_object", baseList);
      if (!baseList) {
        // weightsFile->ls();
        LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
      }

      // Finally, from the top-level TList, get the desired nested TList => the technical problem here is that it can be nested at any level,
      // for that there is a helper utility function getObjectFromList(...) , see its implementation further below
      listWithRuns = reinterpret_cast<TList*>(getObjectFromList(baseList, runNumber));
      if (!listWithRuns) {
        TString runNumberWithLeadingZeroes = "000";
        runNumberWithLeadingZeroes += runNumber; // another try, with "000" prepended to run number
        listWithRuns = reinterpret_cast<TList*>(getObjectFromList(baseList, runNumberWithLeadingZeroes.Data()));
        if (!listWithRuns) {
          LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
        }
      }

      // OK, we got the desired TList with efficiency corrections, after that we can use the common code for all 3 cases (local, AliEn, CCDB, that common code is below)

    } else if (bFileIsInCCDB) {
      // File you want to access is in your home dir in CCDB:
      // Remember that here I do not access the file; instead, I directly access the object in that file.
      // My home dir in CCDB: https://alice-ccdb.cern.ch/browse/Users/a/abilandz/ => adapt for your case
      ccdb->setURL("https://alice-ccdb.cern.ch"); // to be able to use "ccdb" this object in your analysis task, see 4b/ below
      baseList = dynamic_cast<TList*>(ccdb->get<TList>(TString(filePath).ReplaceAll("/alice-ccdb.cern.ch/", "").Data()));
      if (!baseList) {
        LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
      }

      listWithRuns = reinterpret_cast<TList*>(getObjectFromList(baseList, runNumber));
      if (!listWithRuns) {
        TString runNumberWithLeadingZeroes = "000";
        runNumberWithLeadingZeroes += runNumber; // another try, with "000" prepended to run number
        listWithRuns = reinterpret_cast<TList*>(getObjectFromList(baseList, runNumberWithLeadingZeroes.Data()));
        if (!listWithRuns) {
          LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
        }
      }

      // OK, we got the desired TList with efficiency corrections, after that we
      // can use the common code for all 3 cases (local, AliEn, CCDB, that
      // common code is below)
    } else {
      // this is the local case:
      // Check if the external ROOT file exists at the specified path:
      if (gSystem->AccessPathName(filePath, kFileExists)) {
        LOGF(info, "\033[1;33m if(gSystem->AccessPathName(filePath,kFileExists)), filePath = %s \033[0m", filePath);
        LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
      }

      TFile* weightsFile = TFile::Open(filePath, "READ");
      if (!weightsFile) {
        LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
      }

      weightsFile->GetObject("ccdb_object", baseList);

      if (!baseList) {
        // weightsFile->ls();
        LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
      }

      listWithRuns = reinterpret_cast<TList*>(getObjectFromList(baseList, runNumber));
      if (!listWithRuns) {
        TString runNumberWithLeadingZeroes = "000";
        runNumberWithLeadingZeroes += runNumber; // another try, with "000" prepended to run number
        listWithRuns = reinterpret_cast<TList*>(getObjectFromList(baseList, runNumberWithLeadingZeroes.Data()));
        if (!listWithRuns) {
          // baseList->ls();
          // LOGF(fatal, "\033[1;31m%s at line %d : this crash can happen if in the output file there is no list with weights for the current run number = %s\033[0m", __FUNCTION__, __LINE__, tc.fRunNumber.Data());
        }
      }
    }

    // Here comes the common code for all three cases, where from "listWithRuns" you fetch the desired histogram with efficiency corrections:
    // listWithRuns->ls();

    if (!listWithRuns) {
      LOGF(fatal,
           "\033[1;31m%s: listWithRuns is null for run %s\033[0m",
           __FUNCTION__, runNumber);
    }

    hist = dynamic_cast<TH1D*>(listWithRuns->FindObject("h_invert"));
    if (!hist) {
      LOGF(fatal, "%s: histogram 'hist1' not found in run list", __FUNCTION__);
    }
    hist->SetDirectory(0);
    auto histClone = dynamic_cast<TH1D*>(hist->Clone());
    histClone->SetDirectory(0);

    delete baseList; // release back the memory

    return histClone;
  }

  // templates
  template <ERecSim rs, typename T1>
  bool eventCuts(T1 const& collision)
  {
    if constexpr (rs == ERec || rs == ERecAndSim) {
      if (cfVertexZSwitch) // Vertex Z cuts for Rec
      {
        if (collision.posZ() > cfVertexZ.value[1] || collision.posZ() < cfVertexZ.value[0]) {
          return false;
        }
        if constexpr (rs == ERecAndSim) // Vertex Z cuts for Sim
        {
          if (!collision.has_mcCollision()) {
            return false;
          }
          auto mcCollision = collision.mcCollision(); // corresponding MC truth simulated particle
          if (mcCollision.posZ() > cfVertexZ.value[1] || mcCollision.posZ() < cfVertexZ.value[0]) {
            return false;
          }
        }
      }
    }

    return true;
  }

  template <ERecSim rs, ECuts cuts, typename T1, typename T2>
  void eventHistFill(T1 const& collision, T2 const& tracks)
  {
    auto thisCent = collision.centFT0C(); // use auto to determine the type
    switch (centralityEstimator) {
      case EFT0C:
        thisCent = collision.centFT0C();
        break;
      case EFT0M:
        thisCent = collision.centFT0M();
        break;
      case EFV0A:
        thisCent = collision.centFV0A();
        break;
      case ENTPV:
        thisCent = collision.centNTPV();
        break;
    }

    auto thisRefMult = collision.multTPC(); // use auto to determine the type
    switch (multiplicityTables) {
      case EMultTPC:
        thisRefMult = collision.multTPC();
        break;
      case EMultFV0M:
        thisRefMult = collision.multFV0M();
        break;
      case EMultFT0C:
        thisRefMult = collision.multFT0C();
        break;
      case EMultFT0M:
        thisRefMult = collision.multFT0M();
        break;
      case EMultNTracksPV:
        thisRefMult = collision.multNTracksPV();
        break;
    }

    if constexpr (rs == ERec || rs == ERecAndSim) {
      // Fill reconstructed-level event histograms
      int multiplicityRec = static_cast<int>(tracks.size());
      if constexpr (cuts == ENoCuts) {
        ec.fEventHistograms[EHistMultiplicity][ERec][ENoCuts]->Fill(multiplicityRec);
        ec.fEventHistograms[EHistCentrality][ERec][ENoCuts]->Fill(thisCent);
        ec.fEventHistograms[EHistReferencEMultiplicity][ERec][ENoCuts]->Fill(thisRefMult);
        ec.fEventHistograms[EHistVertexX][ERec][ENoCuts]->Fill(collision.posX());
        ec.fEventHistograms[EHistVertexY][ERec][ENoCuts]->Fill(collision.posY());
        ec.fEventHistograms[EHistVertexZ][ERec][ENoCuts]->Fill(collision.posZ());
      }
      if constexpr (cuts == EWithCuts) {
        ec.fEventHistograms[EHistMultiplicity][ERec][EWithCuts]->Fill(multiplicityRec);
        ec.fEventHistograms[EHistCentrality][ERec][EWithCuts]->Fill(thisCent);
        ec.fEventHistograms[EHistReferencEMultiplicity][ERec][EWithCuts]->Fill(thisRefMult);
        ec.fEventHistograms[EHistVertexX][ERec][EWithCuts]->Fill(collision.posX());
        ec.fEventHistograms[EHistVertexY][ERec][EWithCuts]->Fill(collision.posY());
        ec.fEventHistograms[EHistVertexZ][ERec][EWithCuts]->Fill(collision.posZ());
      }

      // Fill MC simulated-level event histograms if both reconstructed and simulated data are processed
      if constexpr (rs == ERecAndSim) {
        if (!collision.has_mcCollision()) {
          return;
        }
        auto mcCollision = collision.mcCollision(); // corresponding MC truth simulated particle
        int multiplicitySim = static_cast<int>(tracks.size());
        float impactParameter = mcCollision.impactParameter();
        float centralitySim = math::PI * impactParameter * impactParameter / sigmaInel; // centrality for sim derived from impact parameter
        if constexpr (cuts == ENoCuts) {
          ec.fEventHistograms[EHistMultiplicity][ESim][ENoCuts]->Fill(multiplicitySim);
          ec.fEventHistograms[EHistCentrality][ESim][ENoCuts]->Fill(centralitySim);
          ec.fEventHistograms[EHistImpactParameter][ESim][ENoCuts]->Fill(mcCollision.impactParameter());
          ec.fEventHistograms[EHistVertexX][ESim][ENoCuts]->Fill(mcCollision.posX());
          ec.fEventHistograms[EHistVertexY][ESim][ENoCuts]->Fill(mcCollision.posY());
          ec.fEventHistograms[EHistVertexZ][ESim][ENoCuts]->Fill(mcCollision.posZ());
        }

        if constexpr (cuts == EWithCuts) {
          ec.fEventHistograms[EHistMultiplicity][ESim][EWithCuts]->Fill(multiplicitySim);
          ec.fEventHistograms[EHistCentrality][ESim][EWithCuts]->Fill(centralitySim);
          ec.fEventHistograms[EHistImpactParameter][ESim][EWithCuts]->Fill(mcCollision.impactParameter());
          ec.fEventHistograms[EHistVertexX][ESim][EWithCuts]->Fill(mcCollision.posX());
          ec.fEventHistograms[EHistVertexY][ESim][EWithCuts]->Fill(mcCollision.posY());
          ec.fEventHistograms[EHistVertexZ][ESim][EWithCuts]->Fill(mcCollision.posZ());
        }
      }
    }
  }

  template <ERecSim rs, typename T>
  bool particlECuts(T const& track)
  {
    if constexpr (rs == ERec || rs == ERecAndSim) {
      if (cfPtSwitch) // Vertex Z cuts for Rec
      {
        if (track.pt() < cfPt.value[0] || track.pt() > cfPt.value[1]) {
          return false;
        }
        if constexpr (rs == ERecAndSim) // Vertex Z cuts for Sim
        {
          if (!track.has_mcParticle()) {
            return false;
          }
          auto mcParticle = track.mcParticle(); // corresponding MC truth simulated particle
          if (mcParticle.pt() < cfPt.value[0] || mcParticle.pt() > cfPt.value[1]) {
            return false;
          }
        }
      }
    }

    return true;
  }

  template <ERecSim rs, ECuts cuts, typename T1>
  void particleHistFill(T1 const& track)
  {
    if constexpr (rs == ERec || rs == ERecAndSim) {
      if constexpr (cuts == ENoCuts) {
        pc.fParticlEHistograms[EHistPt][ERec][ENoCuts]->Fill(track.pt());
        pc.fParticlEHistograms[EHistPhi][ERec][ENoCuts]->Fill(track.phi());
        pc.fParticlEHistograms[EHistEta][ERec][ENoCuts]->Fill(track.eta());
      }

      if constexpr (cuts == EWithCuts) {
        pc.fParticlEHistograms[EHistPt][ERec][EWithCuts]->Fill(track.pt());
        pc.fParticlEHistograms[EHistPhi][ERec][EWithCuts]->Fill(track.phi());
        pc.fParticlEHistograms[EHistEta][ERec][EWithCuts]->Fill(track.eta());
      }

      // ... and corresponding MC truth simulated:
      // See https://github.com/AliceO2Group/O2Physics/blob/master/Tutorials/src/mcHistograms.cxx
      // See https://aliceo2group.github.io/analysis-framework/docs/datamodel/ao2dTables.html#montecarlo
      if constexpr (rs == ERecAndSim) {
        if (!track.has_mcParticle()) {
          LOGF(warning, "  No MC particle for this track, skip...");
          return;
        }
        auto mcParticle = track.mcParticle(); // corresponding MC truth simulated particle
        if constexpr (cuts == ENoCuts) {
          pc.fParticlEHistograms[EHistPt][ESim][ENoCuts]->Fill(mcParticle.pt());
          pc.fParticlEHistograms[EHistPhi][ESim][ENoCuts]->Fill(mcParticle.phi());
          pc.fParticlEHistograms[EHistEta][ESim][ENoCuts]->Fill(mcParticle.eta());
        }

        if constexpr (cuts == EWithCuts) {
          pc.fParticlEHistograms[EHistPt][ESim][EWithCuts]->Fill(mcParticle.pt());
          pc.fParticlEHistograms[EHistPhi][ESim][EWithCuts]->Fill(mcParticle.phi());
          pc.fParticlEHistograms[EHistEta][ESim][EWithCuts]->Fill(mcParticle.eta());
        }
      }
    }
  }

  template <ERecSim rs, typename T1>
  void qaFill(T1 const& collision)
  {
    auto thisCent = collision.centFT0C(); // use auto to determine the type
    switch (centralityEstimator) {
      case EFT0C:
        thisCent = collision.centFT0C();
        break;
      case EFT0M:
        thisCent = collision.centFT0M();
        break;
      case EFV0A:
        thisCent = collision.centFV0A();
        break;
      case ENTPV:
        thisCent = collision.centNTPV();
        break;
    }
    if constexpr (rs == ERecAndSim || rs == ESim) {
      if (!collision.has_mcCollision()) {
        return;
      }
      auto mcCollision = collision.mcCollision();
      float impactParameter = mcCollision.impactParameter();
      float centralitySim = math::PI * impactParameter * impactParameter / sigmaInel; // centrality for sim derived from impact parameter
      qa.fHistCentralityRecSim->Fill(thisCent, centralitySim);
    }
  }

  // *) Define all member functions to be called in the main process* functions:
  template <ERecSim rs, typename T1, typename T2>
  void steer(T1 const& collision, T2 const& tracks)
  {
    // Dry run:
    if (tc.fDryRun) {
      return;
    }

    // Fill Quality Assurance
    if (qualityAssurance) {
      qaFill<rs>(collision);
    }

    // Fill Event Hist
    eventHistFill<rs, ENoCuts>(collision, tracks);
    if (eventCuts<rs>(collision)) {
      eventHistFill<rs, EWithCuts>(collision, tracks);
    }

    // Print current run number:
    LOGF(info, "Run number: %d", collision.bc().runNumber());

    // Print centrality estimated with the selected estimator:
    // LOGF(info, "Centrality: %f", thisCent);

    // Print vertex X position:
    LOGF(info, "Vertex X position: %f", collision.posX());

    // Main loop over particles:
    auto track = tracks.iteratorAt(0); // set the type and scope from one instance
    for (int64_t i = 0; i < tracks.size(); i++) {
      // Print track azimuthal angle:
      // LOGF(info, "Track azimuthal angle: %f", track.phi());

      track = tracks.iteratorAt(i);
      // Fill reconstructed ...:
      particleHistFill<rs, ENoCuts>(track);
      if (eventCuts<rs>(collision) && particlECuts<rs>(track)) {
        particleHistFill<rs, EWithCuts>(track);
      }
    } // end of for (int64_t i = 0; i < tracks.size(); i++) {
  } // end of template <ERecSim rs, typename T1, typename T2> void steer(T1 const& collision, T2 const& tracks) {

  // *) Initialize and book all objects:
  void init(InitContext&)
  {
    // ... code to book and initialize all analysis objects ...

    // *) Set automatically what to process, from an implicit variable "doprocessSomEProcessName" within a PROCESS_SWITCH clause:
    tc.fProcess[EProcessRec] = doprocessRec;
    tc.fProcess[EProcessRecSim] = doprocessRecSim;
    tc.fProcess[EProcessSim] = doprocessSim;

    // *) Configure your task using configurables in the json file:
    tc.fDryRun = cfDryRun;

    // *) Book base list:
    TList* temp = new TList();
    temp->SetOwner(true);
    fBaseList.setObject(temp);

    // *) Book External Hist List
    ex.fExternalHistogramsList = new TList();
    ex.fExternalHistogramsList->SetName("ExternalHistograms");
    ex.fExternalHistogramsList->SetOwner(true);
    fBaseList->Add(ex.fExternalHistogramsList);

    // *) Book and Fill external hist
    ex.fhistWeights = getHistogramWithWeights(cfFileWithWeights.value.c_str(), "000123456");
    ex.fhistWeights->SetTitle(cfFileWithWeights.value.c_str());
    ex.fExternalHistogramsList->Add(ex.fhistWeights);

    // *) Book and nest all other TLists:
    pc.fParticlEHistogramsList = new TList();
    pc.fParticlEHistogramsList->SetName("ParticlEHistograms");
    pc.fParticlEHistogramsList->SetOwner(true);
    fBaseList->Add(pc.fParticlEHistogramsList); // any nested TList in the base TList appears as a subdir in the output ROOT file

    // *) Book pt and phi distribution with binning defined through configurables in the json file:
    vector<float> lPtBins = cfPtBins.value; // define local array and initialize it from an array set in the configurables
    int nBinsPt = static_cast<int>(lPtBins[0]);
    float minPt = lPtBins[1];
    float maxPt = lPtBins[2];

    vector<float> lPhiBins = cfPhiBins.value; // define local array and initialize it from an array set in the configurables
    int nBinsPhi = static_cast<int>(lPhiBins[0]);
    float minPhi = lPhiBins[1];
    float maxPhi = lPhiBins[2];

    vector<float> lEtaBins = cfEtaBins.value; // define local array and initialize it from an array set in the configurables
    int nBinsEta = static_cast<int>(lEtaBins[0]);
    float minEta = lEtaBins[1];
    float maxEta = lEtaBins[2];

    if (doprocessRec || doprocessRecSim) {
      pc.fParticlEHistograms[EHistPt][ERec][ENoCuts] = new TH1F("[EHistPt][ERec][ENoCuts]", "pt distribution for reconstructed particles before cuts", nBinsPt, minPt, maxPt);
      pc.fParticlEHistograms[EHistPt][ERec][ENoCuts]->GetXaxis()->SetTitle("p_{T}");
      pc.fParticlEHistograms[EHistPt][ERec][ENoCuts]->SetColors(kRed, -1, kRed);
      pc.fParticlEHistogramsList->Add(pc.fParticlEHistograms[EHistPt][ERec][ENoCuts]);

      pc.fParticlEHistograms[EHistPhi][ERec][ENoCuts] = new TH1F("[EHistPhi][ERec][ENoCuts]", "phi distribution for reconstructed particles before cuts", nBinsPhi, minPhi, maxPhi);
      pc.fParticlEHistograms[EHistPhi][ERec][ENoCuts]->GetXaxis()->SetTitle("p_{phi}");
      pc.fParticlEHistograms[EHistPhi][ERec][ENoCuts]->SetColors(kRed, -1, kRed);
      pc.fParticlEHistogramsList->Add(pc.fParticlEHistograms[EHistPhi][ERec][ENoCuts]);

      pc.fParticlEHistograms[EHistEta][ERec][ENoCuts] = new TH1F("[EHistEta][ERec][ENoCuts]", "eta distribution for reconstructed particles before cuts", nBinsEta, minEta, maxEta);
      pc.fParticlEHistograms[EHistEta][ERec][ENoCuts]->GetXaxis()->SetTitle("#eta");
      pc.fParticlEHistograms[EHistEta][ERec][ENoCuts]->SetColors(kRed, -1, kRed);
      pc.fParticlEHistogramsList->Add(pc.fParticlEHistograms[EHistEta][ERec][ENoCuts]);

      if (cfPtSwitch) {
        pc.fParticlEHistograms[EHistPt][ERec][EWithCuts] = new TH1F("[EHistPt][ERec][EWithCuts]", "pt distribution for reconstructed particles after cuts", nBinsPt, minPt, maxPt);
        pc.fParticlEHistograms[EHistPt][ERec][EWithCuts]->GetXaxis()->SetTitle("p_{T}");
        pc.fParticlEHistograms[EHistPt][ERec][EWithCuts]->SetColors(kGreen, -1, kGreen);
        pc.fParticlEHistogramsList->Add(pc.fParticlEHistograms[EHistPt][ERec][EWithCuts]);

        pc.fParticlEHistograms[EHistPhi][ERec][EWithCuts] = new TH1F("[EHistPhi][ERec][EWithCuts]", "phi distribution for reconstructed particles after cuts", nBinsPhi, minPhi, maxPhi);
        pc.fParticlEHistograms[EHistPhi][ERec][EWithCuts]->GetXaxis()->SetTitle("p_{phi}");
        pc.fParticlEHistograms[EHistPhi][ERec][EWithCuts]->SetColors(kGreen, -1, kGreen);
        pc.fParticlEHistogramsList->Add(pc.fParticlEHistograms[EHistPhi][ERec][EWithCuts]);

        pc.fParticlEHistograms[EHistEta][ERec][EWithCuts] = new TH1F("[EHistEta][ERec][EWithCuts]", "eta distribution for reconstructed particles after cuts", nBinsEta, minEta, maxEta);
        pc.fParticlEHistograms[EHistEta][ERec][EWithCuts]->GetXaxis()->SetTitle("#eta");
        pc.fParticlEHistograms[EHistEta][ERec][EWithCuts]->SetColors(kGreen, -1, kGreen);
        pc.fParticlEHistogramsList->Add(pc.fParticlEHistograms[EHistEta][ERec][EWithCuts]);
      }
    }

    if (doprocessSim || doprocessRecSim) {
      pc.fParticlEHistograms[EHistPt][ESim][ENoCuts] = new TH1F("[EHistPt][ESim][ENoCuts]", "pt distribution for simulated particles  before cuts", nBinsPt, minPt, maxPt);
      pc.fParticlEHistograms[EHistPt][ESim][ENoCuts]->GetXaxis()->SetTitle("p_{T}");
      pc.fParticlEHistograms[EHistPt][ESim][ENoCuts]->SetColors(kRed, -1, kRed);
      pc.fParticlEHistogramsList->Add(pc.fParticlEHistograms[EHistPt][ESim][ENoCuts]);

      pc.fParticlEHistograms[EHistPhi][ESim][ENoCuts] = new TH1F("[EHistPhi][ESim][ENoCuts]", "phi distribution for simulated particles before cuts", nBinsPhi, minPhi, maxPhi);
      pc.fParticlEHistograms[EHistPhi][ESim][ENoCuts]->GetXaxis()->SetTitle("p_{phi}");
      pc.fParticlEHistograms[EHistPhi][ESim][ENoCuts]->SetColors(kRed, -1, kRed);
      pc.fParticlEHistogramsList->Add(pc.fParticlEHistograms[EHistPhi][ESim][ENoCuts]);

      pc.fParticlEHistograms[EHistEta][ESim][ENoCuts] = new TH1F("[EHistEta][ESim][ENoCuts]", "eta distribution for simulated particles before cuts", nBinsEta, minEta, maxEta);
      pc.fParticlEHistograms[EHistEta][ESim][ENoCuts]->GetXaxis()->SetTitle("#eta");
      pc.fParticlEHistograms[EHistEta][ESim][ENoCuts]->SetColors(kRed, -1, kRed);
      pc.fParticlEHistogramsList->Add(pc.fParticlEHistograms[EHistEta][ESim][ENoCuts]);

      if (cfPtSwitch) {
        pc.fParticlEHistograms[EHistPt][ESim][EWithCuts] = new TH1F("[EHistPt][ESim][EWithCuts]", "pt distribution for simulated particles after cuts", nBinsPt, minPt, maxPt);
        pc.fParticlEHistograms[EHistPt][ESim][EWithCuts]->GetXaxis()->SetTitle("p_{T}");
        pc.fParticlEHistograms[EHistPt][ESim][EWithCuts]->SetColors(kGreen, -1, kGreen);
        pc.fParticlEHistogramsList->Add(pc.fParticlEHistograms[EHistPt][ESim][EWithCuts]);

        pc.fParticlEHistograms[EHistPhi][ESim][EWithCuts] = new TH1F("[EHistPhi][ESim][EWithCuts]", "phi distribution for simulated particles after cuts", nBinsPhi, minPhi, maxPhi);
        pc.fParticlEHistograms[EHistPhi][ESim][EWithCuts]->GetXaxis()->SetTitle("p_{phi}");
        pc.fParticlEHistograms[EHistPhi][ESim][EWithCuts]->SetColors(kGreen, -1, kGreen);
        pc.fParticlEHistogramsList->Add(pc.fParticlEHistograms[EHistPhi][ESim][EWithCuts]);

        pc.fParticlEHistograms[EHistEta][ESim][EWithCuts] = new TH1F("[EHistEta][ESim][EWithCuts]", "eta distribution for simulated particles after cuts", nBinsEta, minEta, maxEta);
        pc.fParticlEHistograms[EHistEta][ESim][EWithCuts]->GetXaxis()->SetTitle("#eta");
        pc.fParticlEHistograms[EHistEta][ESim][EWithCuts]->SetColors(kGreen, -1, kGreen);
        pc.fParticlEHistogramsList->Add(pc.fParticlEHistograms[EHistEta][ESim][EWithCuts]);
      }
    }

    // Book event-level histograms
    ec.fEventHistogramsList = new TList();
    ec.fEventHistogramsList->SetName("EventHistograms");
    ec.fEventHistogramsList->SetOwner(true);
    fBaseList->Add(ec.fEventHistogramsList);

    vector<float> lCent = cfCentBins.value;
    int nBinsCent = static_cast<int>(lCent[0]);
    float minCent = lCent[1];
    float maxCent = lCent[2];

    vector<float> lMultRec = cfMultBinsRec.value;
    int nBinsMultRec = static_cast<int>(lMultRec[0]);
    float minMultRec = lMultRec[1];
    float maxMultRec = lMultRec[2];

    vector<float> lMultRef = cfMultBinsRef.value;
    int nBinsMultRef = static_cast<int>(lMultRef[0]);
    float minMultRef = lMultRef[1];
    float maxMultRef = lMultRef[2];

    vector<float> lMultSim = cfMultBinsSim.value;
    int nBinsMultSim = static_cast<int>(lMultSim[0]);
    float minMultSim = lMultSim[1];
    float maxMultSim = lMultSim[2];

    vector<float> lVx = cfVxBins.value;
    int nBinsVx = static_cast<int>(lVx[0]);
    float minVx = lVx[1];
    float maxVx = lVx[2];

    vector<float> lVy = cfVyBins.value;
    int nBinsVy = static_cast<int>(lVy[0]);
    float minVy = lVy[1];
    float maxVy = lVy[2];

    vector<float> lVz = cfVzBins.value;
    int nBinsVz = static_cast<int>(lVz[0]);
    float minVz = lVz[1];
    float maxVz = lVz[2];

    vector<float> lIp = cfIpBins.value;
    int nBinsIp = static_cast<int>(lIp[0]);
    float minIp = lIp[1];
    float maxIp = lIp[2];

    if (doprocessRec || doprocessRecSim) {
      ec.fEventHistograms[EHistCentrality][ERec][ENoCuts] = new TH1F("[EHistCentrality][ERec][ENoCuts]", "Centrality (reconstructed) before cuts", nBinsCent, minCent, maxCent);
      ec.fEventHistograms[EHistCentrality][ERec][ENoCuts]->GetXaxis()->SetTitle("Centrality");
      ec.fEventHistograms[EHistCentrality][ERec][ENoCuts]->SetColors(kRed, -1, kRed);
      ec.fEventHistogramsList->Add(ec.fEventHistograms[EHistCentrality][ERec][ENoCuts]);

      ec.fEventHistograms[EHistMultiplicity][ERec][ENoCuts] = new TH1F("[EHistMultiplicity][ERec][ENoCuts]", "Multiplicity (reconstructed) before cuts", nBinsMultRec, minMultRec, maxMultRec);
      ec.fEventHistograms[EHistMultiplicity][ERec][ENoCuts]->GetXaxis()->SetTitle("Multiplicity");
      ec.fEventHistograms[EHistMultiplicity][ERec][ENoCuts]->SetColors(kRed, -1, kRed);
      ec.fEventHistogramsList->Add(ec.fEventHistograms[EHistMultiplicity][ERec][ENoCuts]);

      ec.fEventHistograms[EHistReferencEMultiplicity][ERec][ENoCuts] = new TH1F("[EHistReferencEMultiplicity][ERec][ENoCuts]", "Reference Multiplicity before cuts", nBinsMultRef, minMultRef, maxMultRef);
      ec.fEventHistograms[EHistReferencEMultiplicity][ERec][ENoCuts]->GetXaxis()->SetTitle("Reference Multiplicity");
      ec.fEventHistograms[EHistReferencEMultiplicity][ERec][ENoCuts]->SetColors(kRed, -1, kRed);
      ec.fEventHistogramsList->Add(ec.fEventHistograms[EHistReferencEMultiplicity][ERec][ENoCuts]);

      ec.fEventHistograms[EHistVertexX][ERec][ENoCuts] = new TH1F("[EHistVertexX][ERec][ENoCuts]", "Vertex X (reconstructed) before cuts", nBinsVx, minVx, maxVx);
      ec.fEventHistograms[EHistVertexX][ERec][ENoCuts]->GetXaxis()->SetTitle("Vertex X");
      ec.fEventHistograms[EHistVertexX][ERec][ENoCuts]->SetColors(kRed, -1, kRed);
      ec.fEventHistogramsList->Add(ec.fEventHistograms[EHistVertexX][ERec][ENoCuts]);

      ec.fEventHistograms[EHistVertexY][ERec][ENoCuts] = new TH1F("[EHistVertexY][ERec][ENoCuts]", "Vertex Y (reconstructed) before cuts", nBinsVy, minVy, maxVy);
      ec.fEventHistograms[EHistVertexY][ERec][ENoCuts]->GetXaxis()->SetTitle("Vertex Y");
      ec.fEventHistograms[EHistVertexY][ERec][ENoCuts]->SetColors(kRed, -1, kRed);
      ec.fEventHistogramsList->Add(ec.fEventHistograms[EHistVertexY][ERec][ENoCuts]);

      ec.fEventHistograms[EHistVertexZ][ERec][ENoCuts] = new TH1F("[EHistVertexZ][ERec][ENoCuts]", "Vertex Z (reconstructed) before cuts", nBinsVz, minVz, maxVz);
      ec.fEventHistograms[EHistVertexZ][ERec][ENoCuts]->GetXaxis()->SetTitle("Vertex Z");
      ec.fEventHistograms[EHistVertexZ][ERec][ENoCuts]->SetColors(kRed, -1, kRed);
      ec.fEventHistogramsList->Add(ec.fEventHistograms[EHistVertexZ][ERec][ENoCuts]);

      if (cfVertexZSwitch) {
        ec.fEventHistograms[EHistCentrality][ERec][EWithCuts] = new TH1F("[EHistCentrality][ERec][EWithCuts]", "Centrality (reconstructed) after cuts", nBinsCent, minCent, maxCent);
        ec.fEventHistograms[EHistCentrality][ERec][EWithCuts]->GetXaxis()->SetTitle("Centrality");
        ec.fEventHistograms[EHistCentrality][ERec][EWithCuts]->SetColors(kGreen, -1, kGreen);
        ec.fEventHistogramsList->Add(ec.fEventHistograms[EHistCentrality][ERec][EWithCuts]);

        ec.fEventHistograms[EHistMultiplicity][ERec][EWithCuts] = new TH1F("[EHistMultiplicity][ERec][EWithCuts]", "Multiplicity (reconstructed) after cuts", nBinsMultRec, minMultRec, maxMultRec);
        ec.fEventHistograms[EHistMultiplicity][ERec][EWithCuts]->GetXaxis()->SetTitle("Multiplicity");
        ec.fEventHistograms[EHistMultiplicity][ERec][EWithCuts]->SetColors(kGreen, -1, kGreen);
        ec.fEventHistogramsList->Add(ec.fEventHistograms[EHistMultiplicity][ERec][EWithCuts]);

        ec.fEventHistograms[EHistReferencEMultiplicity][ERec][EWithCuts] = new TH1F("[EHistReferencEMultiplicity][ERec][EWithCuts]", "Reference Multiplicity after cuts", nBinsMultRef, minMultRef, maxMultRef);
        ec.fEventHistograms[EHistReferencEMultiplicity][ERec][EWithCuts]->GetXaxis()->SetTitle("Reference Multiplicity");
        ec.fEventHistograms[EHistReferencEMultiplicity][ERec][EWithCuts]->SetColors(kGreen, -1, kGreen);
        ec.fEventHistogramsList->Add(ec.fEventHistograms[EHistReferencEMultiplicity][ERec][EWithCuts]);

        ec.fEventHistograms[EHistVertexX][ERec][EWithCuts] = new TH1F("[EHistVertexX][ERec][EWithCuts]", "Vertex X (reconstructed) after cuts", nBinsVx, minVx, maxVx);
        ec.fEventHistograms[EHistVertexX][ERec][EWithCuts]->GetXaxis()->SetTitle("Vertex X");
        ec.fEventHistograms[EHistVertexX][ERec][EWithCuts]->SetColors(kGreen, -1, kGreen);
        ec.fEventHistogramsList->Add(ec.fEventHistograms[EHistVertexX][ERec][EWithCuts]);

        ec.fEventHistograms[EHistVertexY][ERec][EWithCuts] = new TH1F("[EHistVertexY][ERec][EWithCuts]", "Vertex Y (reconstructed) after cuts", nBinsVy, minVy, maxVy);
        ec.fEventHistograms[EHistVertexY][ERec][EWithCuts]->GetXaxis()->SetTitle("Vertex Y");
        ec.fEventHistograms[EHistVertexY][ERec][EWithCuts]->SetColors(kGreen, -1, kGreen);
        ec.fEventHistogramsList->Add(ec.fEventHistograms[EHistVertexY][ERec][EWithCuts]);

        ec.fEventHistograms[EHistVertexZ][ERec][EWithCuts] = new TH1F("[EHistVertexZ][ERec][EWithCuts]", "Vertex Z (reconstructed) after cuts", nBinsVz, minVz, maxVz);
        ec.fEventHistograms[EHistVertexZ][ERec][EWithCuts]->GetXaxis()->SetTitle("Vertex Z");
        ec.fEventHistograms[EHistVertexZ][ERec][EWithCuts]->SetColors(kGreen, -1, kGreen);
        ec.fEventHistogramsList->Add(ec.fEventHistograms[EHistVertexZ][ERec][EWithCuts]);
      }
    }

    if (doprocessSim || doprocessRecSim) {
      ec.fEventHistograms[EHistCentrality][ESim][ENoCuts] = new TH1F("[EHistCentrality][ESim][ENoCuts]", "Centrality (simulated) before cuts", nBinsCent, minCent, maxCent);
      ec.fEventHistograms[EHistCentrality][ESim][ENoCuts]->GetXaxis()->SetTitle("Centrality");
      ec.fEventHistograms[EHistCentrality][ESim][ENoCuts]->SetColors(kRed, -1, kRed);
      ec.fEventHistogramsList->Add(ec.fEventHistograms[EHistCentrality][ESim][ENoCuts]);

      ec.fEventHistograms[EHistMultiplicity][ESim][ENoCuts] = new TH1F("[EHistMultiplicity][ESim][ENoCuts]", "Multiplicity (simulated) before cuts", nBinsMultSim, minMultSim, maxMultSim);
      ec.fEventHistograms[EHistMultiplicity][ESim][ENoCuts]->GetXaxis()->SetTitle("Multiplicity");
      ec.fEventHistograms[EHistMultiplicity][ESim][ENoCuts]->SetColors(kRed, -1, kRed);
      ec.fEventHistogramsList->Add(ec.fEventHistograms[EHistMultiplicity][ESim][ENoCuts]);

      ec.fEventHistograms[EHistVertexX][ESim][ENoCuts] = new TH1F("[EHistVertexX][ESim][ENoCuts]", "Vertex X (simulated) before cuts", nBinsVx, minVx, maxVx);
      ec.fEventHistograms[EHistVertexX][ESim][ENoCuts]->GetXaxis()->SetTitle("Vertex X");
      ec.fEventHistograms[EHistVertexX][ESim][ENoCuts]->SetColors(kRed, -1, kRed);
      ec.fEventHistogramsList->Add(ec.fEventHistograms[EHistVertexX][ESim][ENoCuts]);

      ec.fEventHistograms[EHistVertexY][ESim][ENoCuts] = new TH1F("[EHistVertexY][ESim][ENoCuts]", "Vertex Y (simulated) before cuts", nBinsVy, minVy, maxVy);
      ec.fEventHistograms[EHistVertexY][ESim][ENoCuts]->GetXaxis()->SetTitle("Vertex Y");
      ec.fEventHistograms[EHistVertexY][ESim][ENoCuts]->SetColors(kRed, -1, kRed);
      ec.fEventHistogramsList->Add(ec.fEventHistograms[EHistVertexY][ESim][ENoCuts]);

      ec.fEventHistograms[EHistVertexZ][ESim][ENoCuts] = new TH1F("[EHistVertexZ][ESim][ENoCuts]", "Vertex Z (simulated) before cuts", nBinsVz, minVz, maxVz);
      ec.fEventHistograms[EHistVertexZ][ESim][ENoCuts]->GetXaxis()->SetTitle("Vertex Z");
      ec.fEventHistograms[EHistVertexZ][ESim][ENoCuts]->SetColors(kRed, -1, kRed);
      ec.fEventHistogramsList->Add(ec.fEventHistograms[EHistVertexZ][ESim][ENoCuts]);

      ec.fEventHistograms[EHistImpactParameter][ESim][ENoCuts] = new TH1F("[EHistImpactParameter][ESim][ENoCuts]", "Impact Parameter (simulated) before cuts", nBinsIp, minIp, maxIp);
      ec.fEventHistograms[EHistImpactParameter][ESim][ENoCuts]->GetXaxis()->SetTitle("Impact Parameter");
      ec.fEventHistograms[EHistImpactParameter][ESim][ENoCuts]->SetColors(kRed, -1, kRed);
      ec.fEventHistogramsList->Add(ec.fEventHistograms[EHistImpactParameter][ESim][ENoCuts]);

      if (cfVertexZSwitch) {
        ec.fEventHistograms[EHistCentrality][ESim][EWithCuts] = new TH1F("[EHistCentrality][ESim][EWithCuts]", "Centrality (simulated) after cuts", nBinsCent, minCent, maxCent);
        ec.fEventHistograms[EHistCentrality][ESim][EWithCuts]->GetXaxis()->SetTitle("Centrality");
        ec.fEventHistograms[EHistCentrality][ESim][EWithCuts]->SetColors(kGreen, -1, kGreen);
        ec.fEventHistogramsList->Add(ec.fEventHistograms[EHistCentrality][ESim][EWithCuts]);

        ec.fEventHistograms[EHistMultiplicity][ESim][EWithCuts] = new TH1F("[EHistMultiplicity][ESim][EWithCuts]", "Multiplicity (simulated) after cuts", nBinsMultSim, minMultSim, maxMultSim);
        ec.fEventHistograms[EHistMultiplicity][ESim][EWithCuts]->GetXaxis()->SetTitle("Multiplicity");
        ec.fEventHistograms[EHistMultiplicity][ESim][EWithCuts]->SetColors(kGreen, -1, kGreen);
        ec.fEventHistogramsList->Add(ec.fEventHistograms[EHistMultiplicity][ESim][EWithCuts]);

        ec.fEventHistograms[EHistVertexX][ESim][EWithCuts] = new TH1F("[EHistVertexX][ESim][EWithCuts]", "Vertex X (simulated) after cuts", nBinsVx, minVx, maxVx);
        ec.fEventHistograms[EHistVertexX][ESim][EWithCuts]->GetXaxis()->SetTitle("Vertex X");
        ec.fEventHistograms[EHistVertexX][ESim][EWithCuts]->SetColors(kGreen, -1, kGreen);
        ec.fEventHistogramsList->Add(ec.fEventHistograms[EHistVertexX][ESim][EWithCuts]);

        ec.fEventHistograms[EHistVertexY][ESim][EWithCuts] = new TH1F("[EHistVertexY][ESim][EWithCuts]", "Vertex Y (simulated) after cuts", nBinsVy, minVy, maxVy);
        ec.fEventHistograms[EHistVertexY][ESim][EWithCuts]->GetXaxis()->SetTitle("Vertex Y");
        ec.fEventHistograms[EHistVertexY][ESim][EWithCuts]->SetColors(kGreen, -1, kGreen);
        ec.fEventHistogramsList->Add(ec.fEventHistograms[EHistVertexY][ESim][EWithCuts]);

        ec.fEventHistograms[EHistVertexZ][ESim][EWithCuts] = new TH1F("[EHistVertexZ][ESim][EWithCuts]", "Vertex Z (simulated) after cuts", nBinsVz, minVz, maxVz);
        ec.fEventHistograms[EHistVertexZ][ESim][EWithCuts]->GetXaxis()->SetTitle("Vertex Z");
        ec.fEventHistograms[EHistVertexZ][ESim][EWithCuts]->SetColors(kGreen, -1, kGreen);
        ec.fEventHistogramsList->Add(ec.fEventHistograms[EHistVertexZ][ESim][EWithCuts]);

        ec.fEventHistograms[EHistImpactParameter][ESim][EWithCuts] = new TH1F("[EHistImpactParameter][ESim][EWithCuts]", "Impact Parameter (simulated) after cuts", nBinsIp, minIp, maxIp);
        ec.fEventHistograms[EHistImpactParameter][ESim][EWithCuts]->GetXaxis()->SetTitle("Impact Parameter");
        ec.fEventHistograms[EHistImpactParameter][ESim][EWithCuts]->SetColors(kGreen, -1, kGreen);
        ec.fEventHistogramsList->Add(ec.fEventHistograms[EHistImpactParameter][ESim][EWithCuts]);
      }
    }

    if (qualityAssurance && doprocessRecSim) {
      qa.fQualityAssuranceList = new TList();
      qa.fQualityAssuranceList->SetName("QualityAssurance");
      qa.fQualityAssuranceList->SetOwner(true);
      fBaseList->Add(qa.fQualityAssuranceList);

      qa.fHistCentralityRecSim = new TH2F("fHistCentralityRecSim", "Centrality Rec vs Sim", nBinsCent, minCent, maxCent, nBinsCent, minCent, maxCent);
      qa.fHistCentralityRecSim->GetXaxis()->SetTitle("Centrality (reconstructed)");
      qa.fHistCentralityRecSim->GetYaxis()->SetTitle("Centrality (simulated)");
      qa.fQualityAssuranceList->Add(qa.fHistCentralityRecSim);
    }
  } // end of void init(InitContext&) {

  // A) Process only reconstructed data:
  void processRec(CollisionRec const& collision, aod::BCs const&, TracksRec const& tracks)
  {
    // *) steer all analysis steps:
    steer<ERec>(collision, tracks);
  }
  PROCESS_SWITCH(MultiparticleCorrelationsMei, processRec, "process only reconstructed data", true); // yes, keep always one process switch "true", so that there is default running version

  // -------------------------------------------

  // B) Process both reconstructed and corresponding MC truth simulated data:
  void processRecSim(CollisionRecSim const& collision, aod::BCs const&, TracksRecSim const& tracks, aod::McParticles const&, aod::McCollisions const&)
  {
    steer<ERecAndSim>(collision, tracks);
  }
  PROCESS_SWITCH(MultiparticleCorrelationsMei, processRecSim, "process both reconstructed and corresponding MC truth simulated data", false);

  // -------------------------------------------

  // C) Process only simulated data:
  void processSim(CollisionSim const& /*collision*/, aod::BCs const&, TracksSim const& /*tracks*/)
  {
    // steer<ESim>(collision, tracks); // TBI 20241105 not ready yet, but I do not really need this one urgently, since RecSim is working, and I need that one for efficiencies...
  }
  PROCESS_SWITCH(MultiparticleCorrelationsMei, processSim, "process only simulated data", false);

}; // struct MultiparticleCorrelationsMei {

// *) The final touch:
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<MultiparticleCorrelationsMei>(cfgc),
  };
} // WorkflowSpec...
