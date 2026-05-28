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

/// \file multiparticleCumulants.cxx
/// \brief Task for producing multiparticle cumulants
/// \author Pei-Ying Kuan, TU München, pei-ying.kuan@cern.ch

#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h" // needed for aod::TracksDCA table

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/MathConstants.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/InitContext.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <TCollection.h>
#include <TFile.h>
#include <TGrid.h>
#include <TH1.h>
#include <TH2.h>
#include <TList.h>
#include <TObject.h>
#include <TParticlePDG.h>
#include <TString.h>
#include <TSystem.h>

#include <RtypesCore.h>

#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;

// Definitions of join tables for Run 3 analysis:
using EventSelection = soa::Join<aod::EvSels, aod::Mults, aod::CentFT0Ms, aod::CentFV0As, aod::CentNTPVs, aod::CentFT0Cs>;
using CollisionRec = soa::Join<aod::Collisions, EventSelection>::iterator; // use in json "isMC": "true" for "event-selection-task"
using CollisionRecSim = soa::Join<aod::Collisions, aod::McCollisionLabels, EventSelection>::iterator;
using CollisionSim = aod::McCollision;
using TracksRec = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection>;
using TrackRec = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection>::iterator;
using TracksRecSim = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::McTrackLabels>; // + use in json "isMC" : "true"
using TrackRecSim = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::McTrackLabels>::iterator;
using TracksSim = aod::McParticles;
using TrackSim = aod::McParticles::iterator;

using namespace std;

// *) Define enums:
enum EnRlMc { eRl = 0,
              eMc };

enum EnRecSim { eRec = 0,
                eSim,
                eRecAndSim };

enum EnProcess {
  eProcessRec = 0, // Run 3, only reconstructed
  eProcessRecSim,  // Run 3, both reconstructed and simulated
  eProcessSim,     // Run 3, only simulated
  eProcess_N
};

enum EnEventHistograms {
  eCent,
  eMult,
  eVertexX,
  eVertexY,
  eVertexZ,
  eNumContrib,
  eEventHistograms_N
};

const char* eventHistNames[eEventHistograms_N] = {
  "Centrality",
  "Multiplicity",
  "VertexX",
  "VertexY",
  "VertexZ",
  "NumContrib"};

enum EnParticleHistograms {
  ePt,
  ePhi,
  eParticleHistograms_N
};

const char* particleHistNames[eParticleHistograms_N] = {
  "Pt",
  "Phi"};

enum EnQAHistograms {
  eQACent,
  eQAMultNumContrib,
  eQAHistograms_N
};

// *) Main task:
struct MultiparticleCumulants { // this name is used in lower-case format to name the TDirectoryFile in AnalysisResults.root

  // *) Base TList to hold all output objects:
  TString sBaseListName = "Default list name"; // yes, I declare it separately, because I need it also later in BailOut() function
  OutputObj<TList> fBaseList{sBaseListName.Data(),
                             OutputObjHandlingPolicy::AnalysisObject,
                             OutputObjSourceType::OutputObjSource};

  // *) Service:
  Service<ccdb::BasicCCDBManager> ccdb; // support for offline callibration data base
  Service<o2::framework::O2DatabasePDG> pdg;

  // *) Define configurables:
  Configurable<bool> cfDryRun{"cfDryRun", false, "book all histos and run without filling and calculating anything"}; // example for built-in type (float, string, etc.)
  Configurable<std::string> cfCentEstm{"cfCentEstm", "FT0M", "centrality estimator: FT0M, FV0A, NTPV, FT0C"};
  Configurable<std::string> cfMultEstm{"cfMultEstm", "FT0A", "multiplicity estimator: FT0A, FT0C, FT0M, FV0A, FV0C, NTracksPV"};

  Configurable<bool> cfQASwitch{"cfQASwitch", true, "quality assurance switch"};
  Configurable<bool> cfWeightSwitch{"cfWeightSwitch", true, "weight switch"};
  Configurable<bool> cfPrintSwitch{"cfPrintSwitch", true, "printing result switch"};

  // *) Event cut switches
  Configurable<bool> cfVertexZCutSwitch{"cfVertexZCutSwitch", true, "vertex z cut switch"};
  Configurable<bool> cfSel8CutSwitch{"cfSel8CutSwitch", true, "Sel8 cut switch"};
  Configurable<bool> cfCentCutSwitch{"cfCentCutSwitch", true, "centrality cut switch"};
  Configurable<bool> cfNumContribCutSwitch{"cfNumContribCutSwitch", true, "NContribution cut switch"};

  // *) Particle cut switches
  Configurable<bool> cfPtCutSwitch{"cfPtCutSwitch", true, "Pt cut switch"};
  Configurable<bool> cfEtaCutSwitch{"cfEtaCutSwitch", true, "Eta cut switch"};
  Configurable<bool> cfSignCutSwitch{"cfSignCutSwitch", true, "Charge cut switch"};
  Configurable<bool> cfTpcNClsFoundCutSwitch{"cfTpcNClsFoundCutSwitch", true, ""};
  Configurable<bool> cfDCAXYCutSwitch{"cfDCAXYCutSwitch", true, ""};
  Configurable<bool> cfDCAZCutSwitch{"cfDCAZCutSwitch", true, ""};

  // *) Event cut
  Configurable<std::vector<float>> cfVertexZCut{"cfVertexZCut", {-10., 10.}, "vertex z position range: {min, max}[cm]"};
  Configurable<std::vector<float>> cfCentCut{"cfCentCut", {10., 20.}, "centrality range: {min, max}[%]"};
  Configurable<std::vector<float>> cfNumContribCut{"cfNumContribCut", {0, 3000.}, "NContribution range: {min, max}"};

  // *) Particle cut
  Configurable<std::vector<float>> cfPtCut{"cfPtCut", {0.2, 5.0}, "Pt range: {min, max}[GeV], with convention: min <= Pt < max"};
  Configurable<std::vector<float>> cfEtaCut{"cfEtaCut", {-0.8, 0.8}, "Eta range: {min, max}, with convention: min <= Eta < max"};
  Configurable<std::vector<int>> cfSignCut{"cfSignCut", {1, 0, 1}, "sign of charge, 1 to keep and 0 to discard, {negative, neutral, positive}"};
  Configurable<std::vector<float>> cfTpcNClsFoundCut{"cfTpcNClsFoundCut", {70., 160.}, "range of found TPC clusters for this track geometry: {min, max}"};
  Configurable<std::vector<float>> cfDCAXYCut{"cfDCAXYCut", {-3.2, 3.2}, "range of distance-of-closest-approach (DCA) of the extrapolated track to the primary position in XY-direction: {min, max}[cm]"};
  Configurable<std::vector<float>> cfDCAZCut{"cfDCAZCut", {-2.4, 2.4}, "range of distance-of-closest-approach (DCA) of the extrapolated track to the primary position in Z-direction: {min, max}[cm]"};

  // *)
  Configurable<std::string> cfFileWithWeights{"cfFileWithWeights", "/scratch3/go52dab/O2tutorial/tutorial3-6/weights.root", "path to external ROOT file which holds all particle weights in O2 format"};
  Configurable<std::string> cfRunNumber{"cfRunNumber", "000123456", "run number"};

  // *) Bins
  Configurable<std::vector<float>> cfPtBins{"cfPtBins", {1000, 0., 100.}, "nPtBins, ptMin, ptMax"};
  Configurable<std::vector<float>> cfPhiBins{"cfPhiBins", {1000, 0., o2::constants::math::TwoPI}, "nPhiBins, phiMin, phiMax"};
  Configurable<std::vector<float>> cfCentBins{"cfCentBins", {100, 0., 100.}, "nCenBins, cenMin, cenMax"};
  Configurable<std::vector<float>> cfMultBins{"cfMultBins", {100, 0., 5000.}, "nMultBins, MultMin, MultMax"};
  Configurable<std::vector<float>> cfVerXBins{"cfVerXBins", {100, -0.05, 0.05}, "nVerXBins, VerXMin, VerXMax"};
  Configurable<std::vector<float>> cfVerYBins{"cfVerYBins", {100, -0.05, 0.05}, "nVerYBins, VerYMin, VerYMax"};
  Configurable<std::vector<float>> cfVerZBins{"cfVerZBins", {100, -50., 50.}, "nVerZBins, VerZMin, VerZMax"};
  Configurable<std::vector<float>> cfNumContribBins{"cfNumContribBins", {100, 0., 5000.}, "nNumContribBins, NumContribMin, NumContribMax"};

  // *) Define and initialize all data members to be called in the main process* functions:
  // **) Task configuration:
  struct TaskConfiguration {
    bool fProcess[eProcess_N] = {kFALSE}; // Set what to process. See enum EnProcess for full description. Set via implicit variables within a PROCESS_SWITCH clause.
    bool fDryRun = kFALSE;                // book all histos and run without filling and calculating anything
    std::string fCentEstm = "FT0M";
    std::string fMultEstm = "FT0A";

    bool fPrintSwitch = true;

    bool fVertexZCutSwitch = true;
    bool fSel8CutSwitch = true;
    bool fCentCutSwitch = true;
    bool fNumContribCutSwitch = true;

    bool fPtCutSwitch = true;
    bool fEtaCutSwitch = true;
    bool fSignCutSwitch = true;
    bool fTpcNClsFoundCutSwitch = true;
    bool fDCAXYCutSwitch = true;
    bool fDCAZCutSwitch = true;

    std::vector<float> fVertexZCut = {-10., 10.};
    std::vector<float> fCentCut = {10., 20.};
    std::vector<float> fNumContribCut = {0, 3000.};

    std::vector<float> fPtCut = {0.2, 5.0};
    std::vector<float> fEtaCut = {-0.8, 0.8};
    std::vector<int> fSignCut = {1, 0, 1};
    std::vector<float> fTpcNClsFoundCut = {70., 160.};
    std::vector<float> fDCAXYCut = {-3.2, 3.2};
    std::vector<float> fDCAZCut = {-2.4, 2.4};

    std::vector<float> fPtBins = {0., 100.};
    std::vector<float> fPhiBins = {0., o2::constants::math::TwoPI};

    std::vector<float> fCentBins = {0., 100.};
    std::vector<float> fMultBins = {0, 5000};
    std::vector<float> fVerXBins = {-0.05, 0.05};
    std::vector<float> fVerYBins = {-0.05, 0.05};
    std::vector<float> fVerZBins = {-50., 50.};
    std::vector<float> fNumContribBins = {0, 5000};

    std::string fFileWithWeights = "/scratch3/go52dab/O2tutorial/tutorial3-6/weights.root";
    std::string fRunNumber = "000123456";
  } tc; // you have to prepend "tc." for all objects name in this group later in the code

  // **) Particle histograms:
  struct ParticleHistograms {
    TList* fParticleHistogramsList = NULL; //!<! list to hold all control particle histograms
    TH1F* fParticleHistograms[eParticleHistograms_N][2][2] = {{{NULL}}};
  } pc; // you have to prepend "pc." for all objects name in this group later in the code

  // **) Event histograms:
  struct EventHistograms {
    TList* fEventHistogramsList = NULL;                            //!<! list to hold all control event histograms
    TH1F* fEventHistograms[eEventHistograms_N][2][2] = {{{NULL}}}; //! [ type - see enum EnEventHistograms ][reco,sim][before, after event cuts]
  } ev;

  struct QAHistograms {
    bool fQASwitch = kTRUE;
    TList* fQAHistogramsList = NULL;
    TH2F* fQAHistograms[eQAHistograms_N][2] = {{NULL}}; //[type][before/after cut]
  } qa;

  struct WeightHistograms {
    bool fWeightSwitch = kTRUE;
    TList* fWeightHistogramsList = NULL;
    std::vector<TH1F*> fWeightHistograms;
  } wt;

  struct EventByEventQuantities {
    float fReferenceMultiplicity = 0.;
    float fCentrality = 0.;
    float fCentralitySim = 0.;
    float fImpactParameter = 0.;
    float fNumContrib = 0.;
  } ebye;

  template <EnRlMc rm, typename T1>
  bool ctEventCuts(T1 const& collision)
  {
    bool pass = true;

    bool bVertexZCut = true;
    bool bSel8Cut = true;
    bool bCentCut = true;
    bool bNumContribCut = true;

    // *) For real event and MC event
    bVertexZCut = collision.posZ() < tc.fVertexZCut[1] && collision.posZ() > tc.fVertexZCut[0];

    // *) For real event only
    if constexpr (rm == eRl) {
      bSel8Cut = collision.sel8();
      bCentCut = ebye.fCentrality < tc.fCentCut[1] && ebye.fCentrality > tc.fCentCut[0];
      bNumContribCut = ebye.fNumContrib < tc.fNumContribCut[1] && ebye.fNumContrib > tc.fNumContribCut[0];
    }

    // *) For MC event only
    if constexpr (rm == eMc) {
      bCentCut = ebye.fCentralitySim < tc.fCentCut[1] && ebye.fCentralitySim > tc.fCentCut[0];
    }

    if (tc.fVertexZCutSwitch) {
      pass = pass && bVertexZCut;
    }
    if (tc.fSel8CutSwitch) {
      pass = pass && bSel8Cut;
    }
    if (tc.fCentCutSwitch) {
      pass = pass && bCentCut;
    }
    if (tc.fNumContribCutSwitch) {
      pass = pass && bNumContribCut;
    }

    return pass;
  }

  template <EnRlMc rm, typename T1>
  bool ctParticleCuts(T1 const& track)
  {
    bool pass = true;

    bool bPtCut = true;
    bool bEtaCut = true;
    bool bSignCut = true;
    bool bTpcNClsFoundCut = true;
    bool bDCAXYCut = true;
    bool bDCAZCut = true;

    // *) For real event and MC event
    bPtCut = track.pt() < tc.fPtCut[1] && track.pt() > tc.fPtCut[0];
    bEtaCut = track.eta() < tc.fEtaCut[1] && track.eta() > tc.fEtaCut[0];

    // *) For real event only
    if constexpr (rm == eRl) {
      bSignCut = (track.sign() == -1 && tc.fSignCut[0]) || (track.sign() == 0 && tc.fSignCut[1]) || (track.sign() == 1 && tc.fSignCut[2]);
      bTpcNClsFoundCut = track.tpcNClsFound() < tc.fTpcNClsFoundCut[1] && track.tpcNClsFound() > tc.fTpcNClsFoundCut[0];
      bDCAXYCut = track.dcaXY() < tc.fDCAXYCut[1] && track.dcaXY() > tc.fDCAXYCut[0];
      bDCAZCut = track.dcaZ() < tc.fDCAZCut[1] && track.dcaZ() > tc.fDCAZCut[0];
    }

    // *) For mc event only
    if constexpr (rm == eMc) {
      // TDatabasePDG *db= TDatabasePDG::Instance();
      TParticlePDG* particle = pdg->GetParticle(track.pdgCode());

      if (!particle) {
        // LOGF(warning, "PDG code %d not found", track.pdgCode());
        bSignCut = false;
      } else {
        // LOGF(info, "PDG code %d found", track.pdgCode());
        float charge = particle->Charge();
        bSignCut = (charge < 0 && tc.fSignCut[0]) || (charge == 0 && tc.fSignCut[1]) || (charge > 0 && tc.fSignCut[2]);
      }
    }

    if (tc.fPtCutSwitch) {
      pass = pass && bPtCut;
    }
    if (tc.fEtaCutSwitch) {
      pass = pass && bEtaCut;
    }
    if (tc.fSignCutSwitch) {
      pass = pass && bSignCut;
    }
    if (tc.fTpcNClsFoundCutSwitch) {
      pass = pass && bTpcNClsFoundCut;
    }
    if (tc.fDCAXYCutSwitch) {
      pass = pass && bDCAXYCut;
    }
    if (tc.fDCAZCutSwitch) {
      pass = pass && bDCAZCut;
    }

    return pass;
  }

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

  std::vector<TH1F*> getHistogramsWithWeights(const char* filePath, const char* runNumber)
  {
    // a) Return value:
    std::vector<TH1F*> histograms;
    TList* baseList = NULL;     // base top-level list in the TFile, e.g. named "ccdb_object"
    TList* listWithRuns = NULL; // nested list with run-wise TList's holding run-specific weights

    // c) Determine from filePath if the file in on a local machine, or in home dir AliEn, or in CCDB:
    //    Algorithm: If filePath begins with "/alice/data/CCDB/" then it's in home dir AliEn.
    //                     If filePath begins with "/alice-ccdb.cern.ch/" then it's in CCDB. Therefore, files in AliEn and CCDB must be specified with abs path,
    //                     for local files both abs and relative paths are just fine.
    bool bFileIsInAliEn = false;
    bool bFileIsInCCDB = false;

    TString filePathStr(filePath);

    if (filePathStr(0, 7) == "/alice/") {
      bFileIsInAliEn = true;
    } else if (filePathStr(0, 20) == "/alice-ccdb.cern.ch/") {
      bFileIsInCCDB = true;
    }

    if (bFileIsInAliEn) {
      // File you want to access is in your home dir in AliEn:
      TGrid* alien = TGrid::Connect("alien", gSystem->Getenv("USER"), "", ""); // do not forget to add #include <TGrid.h> to the preamble of your analysis task
      if (!alien) {
        LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
      }
      TFile* weightsFile = TFile::Open(Form("alien://%s", filePath), "READ"); // yes, ROOT can open a file transparently, even if it's sitting in AliEn, with this specific syntax
      if (!weightsFile) {
        LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
      }
      weightsFile->GetObject("ccdb_object", baseList);
      if (!baseList) {
        weightsFile->ls();
        LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
      }

      // Finally, from the top-level TList, get the desired nested TList => the technical problem here is that it can be nested at any level,
      // for thare there is a helper utility function getObjectFromList(...) , see its implementation further below
      listWithRuns = reinterpret_cast<TList*>(getObjectFromList(baseList, runNumber));
      if (!listWithRuns) {
        TString runNumberWithLeadingZeroes = "000";
        runNumberWithLeadingZeroes += runNumber; // another try, with "000" prepended to run number
        listWithRuns = reinterpret_cast<TList*>(getObjectFromList(baseList, runNumberWithLeadingZeroes.Data()));
        if (!listWithRuns) {
          LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
        }
      }
    } else if (bFileIsInCCDB) {
      // File you want to access is in your home dir in CCDB:
      // Remember that here I do not access the file; instead, I directly access the object in that file.
      // My home dir in CCDB: https://alice-ccdb.cern.ch/browse/Users/a/abilandz/   =>  adapt for your case
      ccdb->setURL("https://alice-ccdb.cern.ch"); // to be able to use "ccdb" this object in your analysis task, see 4b/ below
      baseList = reinterpret_cast<TList*>(ccdb->get<TList>(TString(filePath).ReplaceAll("/alice-ccdb.cern.ch/", "").Data()));
      baseList->ls();
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

      // OK, we got the desired TList with efficiency corrections, after that we can use the common code for all 3 cases (local, AliEn, CCDB, that common code is below)
    } else {
      // this is the local case, please handle this one now:
      // Check if the external ROOT file exists at specified path:

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
        weightsFile->ls();
        LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
      }

      listWithRuns = reinterpret_cast<TList*>(getObjectFromList(baseList, runNumber)); // baseList->FindObject(runNumber)
      if (!listWithRuns) {
        TString runNumberWithLeadingZeroes = "000";
        runNumberWithLeadingZeroes += runNumber; // another try, with "000" prepended to run number
        listWithRuns = reinterpret_cast<TList*>(getObjectFromList(baseList, runNumberWithLeadingZeroes.Data()));
        if (!listWithRuns) {
          baseList->ls();
          LOGF(fatal, "\033[1;31m%s at line %d : this crash can happen if in the output file there is no list with weights for the current run number = %s\033[0m", __FUNCTION__, __LINE__, runNumber);
        }
      }
    }

    TIter next(listWithRuns);
    TObject* object = nullptr;

    while (true) {

      object = next();
      if (!object) {
        break;
      }

      auto* hist = dynamic_cast<TH1F*>(object);
      if (!hist) {
        continue;
      }

      hist->SetDirectory(0);
      auto* histClone = dynamic_cast<TH1F*>(hist->Clone());
      if (!histClone) {
        LOGF(fatal, "Failed to clone histogram %s", hist->GetName());
      }
      histClone->SetTitle(Form("%s:%s", filePath, histClone->GetName()));
      histograms.push_back(histClone);
    }

    return histograms;
  }

  // *) Define all member functions to be called in the main process* functions:
  template <EnRecSim rs, typename T1, typename T2>
  void Steer(T1 const& collision, T2 const& tracks)
  {

    // Dry run:
    if (tc.fDryRun) {
      return;
    }

    if (tc.fPrintSwitch) {
      // Print current run number:
      LOGF(info, "Run number: %d", collision.bc().runNumber());
    }

    float rlCollisionCent = 0.;
    if (tc.fCentEstm == "FT0M")
      rlCollisionCent = collision.centFT0M();
    else if (tc.fCentEstm == "FV0A")
      rlCollisionCent = collision.centFV0A();
    else if (tc.fCentEstm == "NTPV")
      rlCollisionCent = collision.centNTPV();
    else if (tc.fCentEstm == "FT0C")
      rlCollisionCent = collision.centFT0C();

    float rlCollisionMult = 0.;
    if (tc.fMultEstm == "FT0A")
      rlCollisionMult = collision.multFT0A();
    else if (tc.fMultEstm == "FT0C")
      rlCollisionMult = collision.multFT0C();
    else if (tc.fMultEstm == "FT0M")
      rlCollisionMult = collision.multFT0M();
    else if (tc.fMultEstm == "FV0A")
      rlCollisionMult = collision.multFV0A();
    else if (tc.fMultEstm == "FV0C")
      rlCollisionMult = collision.multFV0C();
    else if (tc.fMultEstm == "NTracksPV")
      rlCollisionMult = collision.multNTracksPV();

    float rlCollisionNumContrib = 0.;
    rlCollisionNumContrib = static_cast<float>(collision.numContrib());

    if (tc.fPrintSwitch) {
      // Print centrality estimated with "FT0M" estimator:
      LOGF(info, "Centrality: %f", rlCollisionCent);

      // Print multiplicity:
      LOGF(info, "Multiplicity: %f", static_cast<float>(rlCollisionMult));

      // Print vertex position:
      LOGF(info, "Vertex X position: %f", collision.posX());
      LOGF(info, "Vertex Y position: %f", collision.posY());
      LOGF(info, "Vertex Z position: %f", collision.posZ());

      // Print NContributors
      LOGF(info, "NContributors: %f", static_cast<float>(rlCollisionNumContrib));
    }
    ebye.fCentrality = rlCollisionCent;
    ebye.fReferenceMultiplicity = rlCollisionMult;
    ebye.fNumContrib = rlCollisionNumContrib;

    if constexpr (rs == eRec || rs == eRecAndSim) {
      ev.fEventHistograms[eCent][eRec][0]->Fill(rlCollisionCent);
      ev.fEventHistograms[eMult][eRec][0]->Fill(rlCollisionMult);
      ev.fEventHistograms[eVertexX][eRec][0]->Fill(collision.posX());
      ev.fEventHistograms[eVertexY][eRec][0]->Fill(collision.posY());
      ev.fEventHistograms[eVertexZ][eRec][0]->Fill(collision.posZ());
      ev.fEventHistograms[eNumContrib][eRec][0]->Fill(rlCollisionNumContrib);

      if (ctEventCuts<eRl>(collision)) {
        ev.fEventHistograms[eCent][eRec][1]->Fill(rlCollisionCent);
        ev.fEventHistograms[eMult][eRec][1]->Fill(rlCollisionMult);
        ev.fEventHistograms[eVertexX][eRec][1]->Fill(collision.posX());
        ev.fEventHistograms[eVertexY][eRec][1]->Fill(collision.posY());
        ev.fEventHistograms[eVertexZ][eRec][1]->Fill(collision.posZ());
        ev.fEventHistograms[eNumContrib][eRec][1]->Fill(rlCollisionNumContrib);
      }

      if constexpr (rs == eRecAndSim) {

        if (!collision.has_mcCollision()) {
          if (tc.fPrintSwitch) {
            LOGF(warning, "  No MC collision for this collision, skip...");
          }
          return;
        }

        auto mccollision = collision.mcCollision(); // McCollisionLabels

        float mcCollisionCent = 0.;

        float b = mccollision.impactParameter() * std::pow(10, -15); // convert fm to m
        float xs = 7.71 * std::pow(10, -28);                         // convert barn to m^2
        mcCollisionCent = o2::constants::math::PI * b * b / xs * 100;

        ebye.fCentralitySim = mcCollisionCent;
        ebye.fImpactParameter = b;

        if (tc.fPrintSwitch) {
          LOGF(info, "mc impact param (fm): %f", mccollision.impactParameter());
          LOGF(info, "mc centrality: %f", mcCollisionCent);
        }

        ev.fEventHistograms[eCent][eSim][0]->Fill(mcCollisionCent);
        ev.fEventHistograms[eVertexX][eSim][0]->Fill(mccollision.posX());
        ev.fEventHistograms[eVertexY][eSim][0]->Fill(mccollision.posY());
        ev.fEventHistograms[eVertexZ][eSim][0]->Fill(mccollision.posZ());

        if (ctEventCuts<eMc>(mccollision)) {
          ev.fEventHistograms[eCent][eSim][1]->Fill(mcCollisionCent);
          ev.fEventHistograms[eVertexX][eSim][1]->Fill(mccollision.posX());
          ev.fEventHistograms[eVertexY][eSim][1]->Fill(mccollision.posY());
          ev.fEventHistograms[eVertexZ][eSim][1]->Fill(mccollision.posZ());
        }

        if (qa.fQASwitch) {
          qa.fQAHistograms[eQACent][0]->Fill(rlCollisionCent, mcCollisionCent);
          qa.fQAHistograms[eQAMultNumContrib][0]->Fill(rlCollisionMult, rlCollisionNumContrib);
          if (ctEventCuts<eRl>(collision) && ctEventCuts<eMc>(mccollision)) {
            qa.fQAHistograms[eQACent][1]->Fill(rlCollisionCent, mcCollisionCent);
            qa.fQAHistograms[eQAMultNumContrib][1]->Fill(rlCollisionMult, rlCollisionNumContrib);
          }
        }
      }
    }

    // Main loop over particles:
    for (auto const& track : tracks) {
      // LOGF(info, "Track azimuthal angle: %f", track.phi());
      // LOGF(info, "Transverse momentum: %f", track.pt());

      // Fill reconstructed ...:
      if constexpr (rs == eRec || rs == eRecAndSim) {

        // Fill track pt distribution:
        pc.fParticleHistograms[ePt][eRec][0]->Fill(track.pt());
        pc.fParticleHistograms[ePhi][eRec][0]->Fill(track.phi());

        if (ctParticleCuts<eRl>(track)) {
          pc.fParticleHistograms[ePt][eRec][1]->Fill(track.pt());
          pc.fParticleHistograms[ePhi][eRec][1]->Fill(track.phi());
        }
        // ...

        // ... and corresponding MC truth simulated:
        // See https://github.com/AliceO2Group/O2Physics/blob/master/Tutorials/src/mcHistograms.cxx
        // See https://aliceo2group.github.io/analysis-framework/docs/datamodel/ao2dTables.html#montecarlo
        if constexpr (rs == eRecAndSim) {
          if (!track.has_mcParticle()) {
            if (tc.fPrintSwitch) {
              LOGF(warning, "  No MC particle for this track, skip...");
            }
            return;
          }
          auto mcparticle = track.mcParticle(); // corresponding MC truth simulated particle
          pc.fParticleHistograms[ePt][eSim][0]->Fill(mcparticle.pt());
          pc.fParticleHistograms[ePhi][eSim][0]->Fill(mcparticle.phi());
          if (ctParticleCuts<eMc>(mcparticle)) {
            pc.fParticleHistograms[ePt][eSim][1]->Fill(mcparticle.pt());
            pc.fParticleHistograms[ePhi][eSim][1]->Fill(mcparticle.phi());
          }
        } // end of if constexpr (rs == eRecAndSim) {

      } // if constexpr (rs == eRec || rs == eRecAndSim) {

    } // end of for (int64_t i = 0; i < tracks.size(); i++) {

  } // end of template <EnRecSim rs, typename T1, typename T2> void Steer(T1 const& collision, T2 const& tracks) {

  template <EnParticleHistograms histType, typename T1>
  void BookParticleHistograms(T1 const& lPcBins, ParticleHistograms& pc)
  {
    // *) structure:
    // hists without cut
    // - name
    // - new hists
    //  └ rec
    //  └ sim
    // hists with cut
    // - name
    // - new hists
    //  └ rec
    //  └ sim

    std::vector<float> lPtBins = lPcBins[histType]; // define local array and initialize it from an array set in the configurables
    int nBinsPt = static_cast<int>(lPtBins[0]);
    float minPt = lPtBins[1];
    float maxPt = lPtBins[2];

    std::string nameRecNocut = std::string("fHist") + particleHistNames[histType] + std::string("[eRec][before cut]");
    std::string nameSimNocut = std::string("fHist") + particleHistNames[histType] + std::string("[eSim][before cut]");
    std::string nameRecNocutfull = particleHistNames[histType] + std::string(" distribution for reconstructed particles");
    std::string nameSimNocutfull = particleHistNames[histType] + std::string(" distribution for simulated particles");

    if (doprocessRec || doprocessRecSim) {
      pc.fParticleHistograms[histType][eRec][0] = new TH1F(nameRecNocut.c_str(), nameRecNocutfull.c_str(), nBinsPt, minPt, maxPt);
      pc.fParticleHistograms[histType][eRec][0]->GetXaxis()->SetTitle(particleHistNames[histType]);
      pc.fParticleHistogramsList->Add(pc.fParticleHistograms[histType][eRec][0]);
    }

    if (doprocessSim || doprocessRecSim) {
      pc.fParticleHistograms[histType][eSim][0] = new TH1F(nameSimNocut.c_str(), nameSimNocutfull.c_str(), nBinsPt, minPt, maxPt);
      pc.fParticleHistograms[histType][eSim][0]->GetXaxis()->SetTitle(particleHistNames[histType]);
      pc.fParticleHistogramsList->Add(pc.fParticleHistograms[histType][eSim][0]);
    }

    std::string nameRecCut = std::string("fHist") + particleHistNames[histType] + std::string("[eRec][after cut]");
    std::string nameSimCut = std::string("fHist") + particleHistNames[histType] + std::string("[eSim][after cut]");
    std::string nameRecCutfull = particleHistNames[histType] + std::string(" distribution for reconstructed particles");
    std::string nameSimCutfull = particleHistNames[histType] + std::string(" distribution for simulated particles");

    if (doprocessRec || doprocessRecSim) {
      pc.fParticleHistograms[histType][eRec][1] = new TH1F(nameRecCut.c_str(), nameRecCutfull.c_str(), nBinsPt, minPt, maxPt);
      pc.fParticleHistograms[histType][eRec][1]->GetXaxis()->SetTitle(particleHistNames[histType]);
      pc.fParticleHistogramsList->Add(pc.fParticleHistograms[histType][eRec][1]);
    }

    if (doprocessSim || doprocessRecSim) {
      pc.fParticleHistograms[histType][eSim][1] = new TH1F(nameSimCut.c_str(), nameSimCutfull.c_str(), nBinsPt, minPt, maxPt);
      pc.fParticleHistograms[histType][eSim][1]->GetXaxis()->SetTitle(particleHistNames[histType]);
      pc.fParticleHistogramsList->Add(pc.fParticleHistograms[histType][eSim][1]);
    }
  }

  template <EnEventHistograms histType, typename T1>
  void BookEventHistograms(T1 const& lEvBins, EventHistograms& ev)
  {
    // *) structure:
    // hists without cut
    // - name
    //  └ centrality + estimator
    //  └ multiplicity + estimator
    //  └ others
    // - new hists
    //  └ rec
    //  └ sim (without multiplicity and nContrib)
    // hists with cut
    // - name
    //  └ centrality + estimator
    //  └ multiplicity + estimator
    //  └ others
    // - new hists
    //  └ rec
    //  └ sim (without multiplicity and nContrib)

    std::vector<float> lCentBins = lEvBins[histType]; // define local array and initialize it from an array set in the configurables
    int nBinsCent = static_cast<int>(lCentBins[0]);
    float minCent = lCentBins[1];
    float maxCent = lCentBins[2];

    std::string nameRecNocut = std::string("fHist") + eventHistNames[histType] + std::string("[eRec][before cut]");
    std::string nameSimNocut = std::string("fHist") + eventHistNames[histType] + std::string("[eSim][before cut]");
    std::string nameRecNocutfull;
    std::string nameSimNocutfull;

    if constexpr (histType == eCent) {
      std::string nameRecNocutfull = tc.fCentEstm + eventHistNames[histType] + std::string(" distribution for reconstructed events");
      std::string nameSimNocutfull = tc.fCentEstm + eventHistNames[histType] + std::string(" distribution for simulated events");
    } else if constexpr (histType == eMult) {
      std::string nameRecNocutfull = tc.fMultEstm + eventHistNames[histType] + std::string(" distribution for reconstructed events");
      std::string nameSimNocutfull = tc.fMultEstm + eventHistNames[histType] + std::string(" distribution for simulated events");
    } else {
      std::string nameRecNocutfull = eventHistNames[histType] + std::string(" distribution for reconstructed events");
      std::string nameSimNocutfull = eventHistNames[histType] + std::string(" distribution for simulated events");
    }

    if (doprocessRec || doprocessRecSim) {
      ev.fEventHistograms[histType][eRec][0] = new TH1F(nameRecNocut.c_str(), nameRecNocutfull.c_str(), nBinsCent, minCent, maxCent);
      ev.fEventHistograms[histType][eRec][0]->GetXaxis()->SetTitle(eventHistNames[histType]);
      ev.fEventHistogramsList->Add(ev.fEventHistograms[histType][eRec][0]);
    }

    if (doprocessSim || doprocessRecSim) {
      if constexpr (histType != eNumContrib && histType != eMult) {
        ev.fEventHistograms[histType][eSim][0] = new TH1F(nameSimNocut.c_str(), nameSimNocutfull.c_str(), nBinsCent, minCent, maxCent);
        ev.fEventHistograms[histType][eSim][0]->GetXaxis()->SetTitle(eventHistNames[histType]);
        ev.fEventHistogramsList->Add(ev.fEventHistograms[histType][eSim][0]);
      } // No nContrib and multiplicity for processSim
    }

    std::string nameRecCut = std::string("fHist") + eventHistNames[histType] + std::string("[eRec][after cut]");
    std::string nameSimCut = std::string("fHist") + eventHistNames[histType] + std::string("[eSim][after cut]");
    std::string nameRecCutfull;
    std::string nameSimCutfull;

    if constexpr (histType == eCent) {
      std::string nameRecCutfull = tc.fCentEstm + eventHistNames[histType] + std::string(" distribution for reconstructed events");
      std::string nameSimCutfull = tc.fCentEstm + eventHistNames[histType] + std::string(" distribution for simulated events");
    } else if constexpr (histType == eMult) {
      std::string nameRecCutfull = tc.fMultEstm + eventHistNames[histType] + std::string(" distribution for reconstructed events");
      std::string nameSimCutfull = tc.fMultEstm + eventHistNames[histType] + std::string(" distribution for simulated events");
    } else {
      std::string nameRecCutfull = eventHistNames[histType] + std::string(" distribution for reconstructed events");
      std::string nameSimCutfull = eventHistNames[histType] + std::string(" distribution for simulated events");
    }

    if (doprocessRec || doprocessRecSim) {
      ev.fEventHistograms[histType][eRec][1] = new TH1F(nameRecCut.c_str(), nameRecCutfull.c_str(), nBinsCent, minCent, maxCent);
      ev.fEventHistograms[histType][eRec][1]->GetXaxis()->SetTitle(eventHistNames[histType]);
      ev.fEventHistogramsList->Add(ev.fEventHistograms[histType][eRec][1]);
    }

    if (doprocessSim || doprocessRecSim) {
      if constexpr (histType != eNumContrib && histType != eMult) {
        ev.fEventHistograms[histType][eSim][1] = new TH1F(nameSimCut.c_str(), nameSimCutfull.c_str(), nBinsCent, minCent, maxCent);
        ev.fEventHistograms[histType][eSim][1]->GetXaxis()->SetTitle(eventHistNames[histType]);
        ev.fEventHistogramsList->Add(ev.fEventHistograms[histType][eSim][1]);
      }
    }
  }

  template <EnEventHistograms histType, typename T1>
  void BookQAHistograms(T1 const& lQABins, QAHistograms& qa)
  {
    // *) structure:
    // hists without cut
    // └ simulated centrality vs. reconstructed centrality
    // └ nContrib vs. multiplicity
    // └ others
    // hists with cut
    // └ simulated centrality vs. reconstructed centrality
    // └ nContrib vs. multiplicity
    // └ others

    std::vector<float> lCentBins = lQABins[histType];
    int nBinsCent = static_cast<int>(lCentBins[0]);
    float minCent = lCentBins[1];
    float maxCent = lCentBins[2];

    std::string nameNocut = std::string("fHist") + eventHistNames[histType] + std::string("[before cut]");
    std::string nameNocutfull;
    if constexpr (histType == eCent) {
      std::string nameNocutfull = std::string("Quality assurance of ") + tc.fCentEstm + eventHistNames[histType];
    } else if constexpr (histType == eNumContrib) {
      std::string nameNocutfull = std::string("Quality assurance of ") + tc.fMultEstm + eventHistNames[histType] + std::string(" vs. NContributors");
    } else {
      std::string nameNocutfull = std::string("Quality assurance of ") + eventHistNames[histType];
    }

    qa.fQAHistograms[histType][0] = new TH2F(nameNocut.c_str(), nameNocutfull.c_str(), nBinsCent, minCent, maxCent, nBinsCent, minCent, maxCent);
    if constexpr (histType == eNumContrib) {
      qa.fQAHistograms[histType][0]->GetYaxis()->SetTitle("NContributors");
      qa.fQAHistograms[histType][0]->GetXaxis()->SetTitle("Reference multiplicity");
    } else {
      qa.fQAHistograms[histType][0]->GetYaxis()->SetTitle(Form("Simulated %s", eventHistNames[histType]));
      qa.fQAHistograms[histType][0]->GetXaxis()->SetTitle(Form("Reconstructed %s", eventHistNames[histType]));
    }
    qa.fQAHistogramsList->Add(qa.fQAHistograms[histType][0]);

    std::string nameCut = std::string("fHist") + eventHistNames[histType] + std::string("[after cut]");
    std::string nameCutfull;
    if constexpr (histType == eCent) {
      std::string nameCutfull = std::string("Quality assurance of ") + tc.fCentEstm + eventHistNames[histType];
    } else if constexpr (histType == eNumContrib) {
      std::string nameCutfull = std::string("Quality assurance of ") + tc.fMultEstm + eventHistNames[histType];
    } else {
      std::string nameCutfull = std::string("Quality assurance of ") + eventHistNames[histType];
    }

    qa.fQAHistograms[histType][1] = new TH2F(nameCut.c_str(), nameCutfull.c_str(), nBinsCent, minCent, maxCent, nBinsCent, minCent, maxCent);
    if constexpr (histType == eNumContrib) {
      qa.fQAHistograms[histType][1]->GetYaxis()->SetTitle("NContributors");
      qa.fQAHistograms[histType][1]->GetXaxis()->SetTitle("Reference multiplicity");
    } else {
      qa.fQAHistograms[histType][1]->GetYaxis()->SetTitle(Form("Simulated %s", eventHistNames[histType]));
      qa.fQAHistograms[histType][1]->GetXaxis()->SetTitle(Form("Reconstructed %s", eventHistNames[histType]));
    }
    qa.fQAHistogramsList->Add(qa.fQAHistograms[histType][1]);
  }

  // *) Initialize and book all objects:
  void init(InitContext&)
  {

    // ... code to book and initialize all analysis objects ...

    // *) Set automatically what to process, from an implicit variable "doprocessSomeProcessName" within a PROCESS_SWITCH clause:
    tc.fProcess[eProcessRec] = doprocessRec;
    tc.fProcess[eProcessRecSim] = doprocessRecSim;
    tc.fProcess[eProcessSim] = doprocessSim;

    // *) Configure your task using configurables in the json file:
    tc.fDryRun = cfDryRun;
    tc.fCentEstm = cfCentEstm;
    tc.fMultEstm = cfMultEstm;

    tc.fVertexZCutSwitch = cfVertexZCutSwitch;
    tc.fSel8CutSwitch = cfSel8CutSwitch;
    tc.fCentCutSwitch = cfCentCutSwitch;
    tc.fNumContribCutSwitch = cfNumContribCutSwitch;

    tc.fPtCutSwitch = cfPtCutSwitch;
    tc.fEtaCutSwitch = cfEtaCutSwitch;
    tc.fSignCutSwitch = cfSignCutSwitch;
    tc.fTpcNClsFoundCutSwitch = cfTpcNClsFoundCutSwitch;
    tc.fDCAXYCutSwitch = cfDCAXYCutSwitch;
    tc.fDCAZCutSwitch = cfDCAZCutSwitch;

    tc.fVertexZCut = cfVertexZCut;
    tc.fCentCut = cfCentCut;
    tc.fNumContribCut = cfNumContribCut;

    tc.fPtCut = cfPtCut;
    tc.fEtaCut = cfEtaCut;
    tc.fSignCut = cfSignCut;
    tc.fTpcNClsFoundCut = cfTpcNClsFoundCut;
    tc.fDCAXYCut = cfDCAXYCut;
    tc.fDCAZCut = cfDCAZCut;

    tc.fPtBins = cfPtBins;
    tc.fPhiBins = cfPhiBins;

    tc.fCentBins = cfCentBins;
    tc.fMultBins = cfMultBins;
    tc.fVerXBins = cfVerXBins;
    tc.fVerYBins = cfVerYBins;
    tc.fVerZBins = cfVerZBins;
    tc.fNumContribBins = cfNumContribBins;

    tc.fPrintSwitch = cfPrintSwitch;
    tc.fFileWithWeights = cfFileWithWeights;
    tc.fRunNumber = cfRunNumber;

    qa.fQASwitch = cfQASwitch;
    wt.fWeightSwitch = cfWeightSwitch;

    // *) Book base list:
    TList* temp = new TList();
    temp->SetOwner(kTRUE);
    fBaseList.setObject(temp);

    // *) Book and nest all other TLists:
    pc.fParticleHistogramsList = new TList();
    pc.fParticleHistogramsList->SetName("ParticleHistograms");
    pc.fParticleHistogramsList->SetOwner(kTRUE);
    fBaseList->Add(pc.fParticleHistogramsList); // any nested TList in the base TList appears as a subdir in the output ROOT file

    ev.fEventHistogramsList = new TList();
    ev.fEventHistogramsList->SetName("EventHistograms");
    ev.fEventHistogramsList->SetOwner(kTRUE);
    fBaseList->Add(ev.fEventHistogramsList);

    qa.fQAHistogramsList = new TList();
    qa.fQAHistogramsList->SetName("QualityAssuranceHistograms");
    qa.fQAHistogramsList->SetOwner(kTRUE);
    fBaseList->Add(qa.fQAHistogramsList);

    wt.fWeightHistogramsList = new TList();
    wt.fWeightHistogramsList->SetName("WeightHistograms");
    wt.fWeightHistogramsList->SetOwner(kTRUE);
    fBaseList->Add(wt.fWeightHistogramsList);

    std::vector<std::vector<float>> lPcBins = {tc.fPtBins, tc.fPhiBins};
    std::vector<std::vector<float>> lEvBins = {tc.fCentBins, tc.fMultBins, tc.fVerXBins, tc.fVerYBins, tc.fVerZBins, tc.fNumContribBins};
    std::vector<std::vector<float>> lQABins = {tc.fCentBins, tc.fMultBins, tc.fNumContribBins};

    BookParticleHistograms<ePt>(lPcBins, pc);
    BookParticleHistograms<ePhi>(lPcBins, pc);
    BookEventHistograms<eCent>(lEvBins, ev);
    BookEventHistograms<eMult>(lEvBins, ev);
    BookEventHistograms<eVertexX>(lEvBins, ev);
    BookEventHistograms<eVertexY>(lEvBins, ev);
    BookEventHistograms<eVertexZ>(lEvBins, ev);
    BookEventHistograms<eNumContrib>(lEvBins, ev);
    BookQAHistograms<eCent>(lQABins, qa);
    BookQAHistograms<eMult>(lQABins, qa);

    if (wt.fWeightSwitch) {
      wt.fWeightHistograms = getHistogramsWithWeights(tc.fFileWithWeights.c_str(), tc.fRunNumber.c_str());
      for (auto* hist : wt.fWeightHistograms) {
        wt.fWeightHistogramsList->Add(hist);
      }
    }

  } // end of void init(InitContext&) {

  // A) Process only reconstructed data:
  void processRec(CollisionRec const& collision, aod::BCs const&, TracksRec const& tracks)
  {
    // ...

    // *) Steer all analysis steps:
    Steer<eRec>(collision, tracks);
  }
  PROCESS_SWITCH(MultiparticleCumulants, processRec, "process only reconstructed data", true); // yes, keep always one process switch "true", so that there is default running version

  // -------------------------------------------

  // B) Process both reconstructed and corresponding MC truth simulated data:
  void processRecSim(CollisionRecSim const& collision, aod::BCs const&, TracksRecSim const& tracks, aod::McParticles const&, aod::McCollisions const&)
  {
    Steer<eRecAndSim>(collision, tracks);
  }
  PROCESS_SWITCH(MultiparticleCumulants, processRecSim, "process both reconstructed and corresponding MC truth simulated data", false);

  // -------------------------------------------

  // C) Process only simulated data:
  void processSim(CollisionSim const& /*collision*/, aod::BCs const&, TracksSim const& /*tracks*/)
  {
    // Steer<eSim>(collision, tracks); // TBI 20241105 not ready yet, but I do not really need this one urgently, since RecSim is working, and I need that one for efficiencies...
  }
  PROCESS_SWITCH(MultiparticleCumulants, processSim, "process only simulated data", false);

}; // struct MultiparticleCumulants {

// *) The final touch:
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<MultiparticleCumulants>(cfgc),
  };
} // WorkflowSpec...
