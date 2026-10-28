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

/// \file multiparticleCorrelationsDm.cxx
/// \brief Task to calculate multiparticle correlations and related observables
/// \author molnard@tcd.ie

// O2:

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h" // needed for aod::TracksDCA table

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/MathConstants.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/DataTypes.h>
#include <Framework/runDataProcessing.h>

using namespace o2;
using namespace o2::framework;

// Definitions of join tables for Run 3 analysis:
using EventSelection = soa::Join<aod::EvSels, aod::Mults, aod::CentFT0Cs, aod::CentFT0As, aod::CentFT0Ms>;
using CollisionRec = soa::Join<aod::Collisions, EventSelection>::iterator; // use in json "isMC": "true" for "event-selection-task"
using CollisionRecSim = soa::Join<aod::Collisions, aod::McCollisionLabels, EventSelection>::iterator;
using CollisionSim = aod::McCollision;
using TracksRec = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection>;
using TrackRec = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection>::iterator;
using TracksRecSim = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::McTrackLabels>; // + use in json "isMC" : "true"
using TrackRecSim = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::McTrackLabels>::iterator;
using TracksSim = aod::McParticles;
using TrackSim = aod::McParticles::iterator;

// *) ROOT and C++ include's:
#include <TFile.h>
#include <TGrid.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TList.h>
#include <TString.h>
#include <TSystem.h>

#include <string>
#include <vector>

// *) Define enums:
enum ERecSim { ERec = 0,
               ESim,
               ERecAndSim };

enum EProcess {
  EProcessRec = 0, // Run 3, only reconstructed
  EProcessRecSim,  // Run 3, both reconstructed and simulated
  EProcessSim,     // Run 3, only simulated
  EProcess_N
};

enum ECentralityEstimator {
  ECentFT0C = 0, // default
  ECentFT0M,
  ECentFT0A,
  ECentralityEstimator_N
};

enum EReferenceMultiplicityEstimator {
  ERefMultTPC = 0,
  ERefMultFT0C,
  ERefMultFT0M,
  EReferenceMultiplicityEstimator_N
};

enum EEventHistograms { // this is for cuts
  EMultiplicity = 0,
  EReferenceMultiplicity,
  ECentrality,
  EImpactParameter, // only simulated
  EVertexX,
  EVertexY,
  EVertexZ,
  EEventHistograms_N
};

enum EParticleHistograms { // this is for cuts
  EPt = 0,
  EPhi,
  EParticleHistograms_N
};

enum EBeforeAfter {
  EBefore = 0,
  EAfter,
  EBeforeAfter_N
};

// *) Main task:
struct StudentO2Physics { // this name is used in lower-case format to name the TDirectoryFile in AnalysisResults.root

  // *) Base TList to hold all output objects:
  TString sBaseListName = "Default list name"; // yes, I declare it separately, because I need it also later in BailOut() function
  OutputObj<TList> fBaseList{sBaseListName.Data(),
                             OutputObjHandlingPolicy::AnalysisObject,
                             OutputObjSourceType::OutputObjSource};

  // *) CCDB:
  // Service<ccdb::BasicCCDBManager> ccdb; // support for offline callibration data base, not needed for the time being...
  // to be able to use ccdb object in the analysis task
  Service<ccdb::BasicCCDBManager> ccdb;

  // *) Define configurables:

  Configurable<bool> cfDryRun{"cfDryRun", false, "book all histos and run without filling and calculating anything"};
  Configurable<bool> cfFillQA{"cfFillQA", false, "book and fill QA histograms (only for development and validation)"}; // switch for QA histograms
  Configurable<bool> cfUseVertexZCut{"cfUseVertexZCut", true, "apply the vertex Z event cut"};
  Configurable<bool> cfUsePtCut{"cfUsePtCut", true, "apply the pt particle cut"};

  Configurable<std::string> cfCentralityEstimator{"cfCentralityEstimator", "FT0C", "centrality estimator: FT0C, FT0M, FT0A"};                                  // centrality estimator with valid values: FT0C (default)
  Configurable<std::string> cfReferenceMultiplicityEstimator{"cfReferenceMultiplicityEstimator", "FT0C", "reference multiplicity estimator: TPC, FT0C, FT0M"}; // ref. multiplicity estimator with valid values: TPC (default) FT0C Ft0M

  // sigma_inel in fm^2 (1 b = 100 fm^2): Pb+Pb at 5.36 TeV = 7.71 b = 771 fm^2 (PH2278 Lecture 2, slide 33), change for other system or energy
  Configurable<float> cfSigmaInel{"cfSigmaInel", 771., "sigma_inel in fm^2 (1 b = 100 fm^2), Pb+Pb 5.36 TeV = 771"};

  Configurable<std::vector<float>> cfPtBins{"cfPtBins", {1000, 0., 100.}, "nPtBins, ptMin, ptMax"};
  Configurable<std::vector<float>> cfPhiBins{"cfPhiBins", {360, 0., o2::constants::math::TwoPI}, "nPhiBins, phiMin, phiMax"};
  Configurable<std::vector<float>> cfMultiplicityBins{"cfMultiplicityBins", {5000, 0., 5000.}, "nMultBins, multMin, multMax"};
  Configurable<std::vector<float>> cfRefmultiplicityBins{"cfRefmultiplicityBins", {5000, 0., 5000.}, "nRefMultBins, refMultMin, refMultMax"}; // reference multiplicity
  Configurable<std::vector<float>> cfImpactparameterBins{"cfImpactparameterBins", {200, 0., 20.}, "nBBins, bMin, bMax (fm)"};                 // impact parameter
  Configurable<std::vector<float>> cfCentralityBins{"cfCentralityBins", {100, 0., 100.}, "nCentBins, centMin, centMax"};
  Configurable<std::vector<float>> cfVertexXBins{"cfVertexXBins", {200, -1., 1.}, "nvxBins, vxMin, vxMax"};
  Configurable<std::vector<float>> cfVertexYBins{"cfVertexYBins", {200, -1, 1.}, "nvyBins, vyMin, vyMax"};
  Configurable<std::vector<float>> cfVertexZBins{"cfVertexZBins", {300, -15., 15.}, "nvzBins, vzMin, vzMax"};
  Configurable<std::vector<float>> cfVertexZ{"cfVertexZ", {-10., 10.}, "{min, max} [cm]"};
  Configurable<std::vector<float>> cfPt{"cfPt", {0.2, 5.0}, "{min, max} [GeV]"};

  // for weights (efficiency corrections)
  Configurable<std::string> cfFileWithWeights{"cfFileWithWeights", "/scratch8/go26dek/tutoarials/3/DatasetsO2/weights.root", "path to external ROOT file which holds all particle weights"};

  // *) Define and initialize all data members to be called in the main process* functions:
  // **) Task configuration:
  struct TaskConfiguration {
    bool fProcess[EProcess_N] = {false}; // Set what to process. See enum EProcess for full description. Set via implicit variables within a PROCESS_SWITCH clause.
    bool fDryRun = false;                // book all histos and run without filling and calculating anything
    bool fFillQA = false;                // book and fill QA histos
    bool fUseVertexZCut = true;
    bool fUsePtCut = true;

    int fCentralityEstimator = ECentFT0C;
    int fReferenceMultiplicityEstimator = ERefMultTPC;

    float fSigmaInel = 771.;         // sigma inel [fm^2]
    float fVertexZ[2] = {-10., 10.}; // {min, max} [cm]
    float fPt[2] = {0.2, 5.0};
    bool fWeightsLoaded = false; // weights are loaded once, at the first collision

  } tc; // you have to prepend "tc." for all objects name in this group later in the code

  // **) Particle histograms:
  struct ParticleHistograms {
    TList* fParticleHistogramsList = NULL; //!<! list to hold all control particle histograms

    TH1F* fParticleHistograms[EParticleHistograms_N][2][2] = {{{NULL}}}; //! [ type - see enum EParticleHistograms][rec, sim][before, after partic. cuts]

  } pc; // you have to prepend "pc." for all objects name in this group later in the code

  // **) Event histograms:
  struct EventHistograms {
    TList* fEventHistogramsList = NULL;

    TH1F* fEventHistograms[EEventHistograms_N][2][2] = {{{NULL}}}; //! [ type -see enum EEventHistograms][rec, sim][before, after event cuts]

  } ec; // you have to prepend "ec." for all objects name in this group later in the code

  // **) QA histograms (booked and filled only if cfFillQA is true):

  struct QAHistograms {
    TList* fQAHistogramsList = NULL;
    TH2F* fHistCentralityRecVsSim = NULL; // x = centrality sim (from impact parameter), y = centrality rec
  } qa;                                   // you have to prepend "qa." for all objects name in this group later in the code

  // **) Histogram with weights from an external file
  TH1F* histWeights = NULL;

  // *) Define all member functions to be called in the main process* functions:

  // returns centrality from the estimator selected in json
  template <typename T>
  float getCentrality(T const& collision)
  {
    switch (tc.fCentralityEstimator) {
      case ECentFT0C:
        return collision.centFT0C();
      case ECentFT0M:
        return collision.centFT0M();
      case ECentFT0A:
        return collision.centFT0A();
      default:
        LOGF(fatal, "Unknown centrality estimator");
        return -1.;
    }
  }

  // returns reference multiplicity from the estimator selected in json
  template <typename T>
  float getReferenceMultiplicity(T const& collision)
  {
    switch (tc.fReferenceMultiplicityEstimator) {
      case ERefMultTPC:
        return collision.multTPC();
      case ERefMultFT0C:
        return collision.multFT0C();
      case ERefMultFT0M:
        return collision.multFT0M();
      default:
        LOGF(fatal, "Unknown reference multiplicity estimator");
        return -1.;
    }
  }

  //  centrality from impact parameter: slide 33 PH2278 Lecture 2
  float centralityFromImpactParameter(float b) { return 100. * o2::constants::math::PI * b * b / tc.fSigmaInel; }

  // ______________________ weights _________________________

  TObject* getObjectFromList(TList* list, const char* objectName)
  {
    // Last update: 20261018

    // Get TObject pointer from TList, even if it's in some nested TList.
    // Foreseen to be used to fetch histograms or profiles from files directly.

    // Some ideas taken from TCollection::ls()

    // If you have added histograms directly to files (without TList's), then you can fetch them directly with
    // file->Get("hist-name").

    // Example usage: TH1D *hist = (TH1D*) getObjectFromList("some-valid-TList-pointer","some-object-name");

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

  } // TObject* getObjectFromList(...)

  TH1F* getHistogramWithWeights(const char* filePath, const char* runNumber)
  {
    // Last update: 20261018

    // *) Return value:
    TH1F* hist = NULL;
    TList* baseList = NULL;     // base top-level list in the TFile, e.g. named "ccdb_object"
    TList* listWithRuns = NULL; // nested list with run-wise TList's holding run-specific weights

    // *) Determine from filePath if the file is on a local machine, or in home dir AliEn, or in CCDB:

    //    Algorithm: a) If filePath begins with "/alice/cern.ch/" then it's in the home dir AliEn;
    //               b) If filePath begins with "/alice-ccdb.cern.ch/" then it's in CCDB. Therefore, files in AliEn and CCDB must be specified with abs path;
    //               c) For local files both absolute and relative paths are just fine.

    bool bFileIsInAliEn = false;
    bool bFileIsInCCDB = false;

    if (TString(filePath).BeginsWith("/alice/cern.ch/")) {
      bFileIsInAliEn = true; // file is in AliEn
    } else if (TString(filePath).BeginsWith("/alice-ccdb.cern.ch/")) {
      bFileIsInCCDB = true; // file is in CCDB
    }

    if (bFileIsInAliEn) {

      // File you want to access is in your home directory in AliEn.
      // In this branch, there is a special code to support file access in that case.

      const TGrid* alien = TGrid::Connect("alien", gSystem->Getenv("USER"), "", "");
      // add #include <TGrid.h> for it to work

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

      // Finally, from the top-level TList "ccdb_object", get the desired nested TList.
      // The technical problem here is that it can be nested at any level,
      // to overcome that there is a helper utility function TObject* getObjectFromList(...),
      // whose implementation is discussed and provided separately.
      listWithRuns = reinterpret_cast<TList*>(getObjectFromList(baseList, runNumber));
      if (!listWithRuns) {
        TString runNumberWithLeadingZeroes = "000";
        runNumberWithLeadingZeroes += runNumber; // another try, with "000" prepended to run number
        listWithRuns = reinterpret_cast<TList*>(getObjectFromList(baseList, runNumberWithLeadingZeroes.Data()));
        if (!listWithRuns) {
          LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
        }
      }

      // OK, we got the desired run-specific TList with efficiency corrections for AliEn case,
      // after that we can use the common code for all 3 cases (local, AliEn, CCDB),
      // that common code is further below.

    } else if (bFileIsInCCDB) {

      // File you want to access is in your home directory in CCDB.
      // In this branch, there is a special code to support file access in that case.

      // Remember that here I do not access the file; instead, I directly access the object in that file.
      // My home dir in CCDB: https://alice-ccdb.cern.ch/browse/Users/a/abilandz/ => adapt for your case
      ccdb->setURL("http://alice-ccdb.cern.ch"); // to be able to use the "ccdb" object in your analysis task
      baseList = reinterpret_cast<TList*>(ccdb->get<TList>(TString(filePath).ReplaceAll("/alice-ccdb.cern.ch/", "").Data()));

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

      // OK, we got the desired run-specific TList with efficiency corrections for CCDB,
      // after that we can use the common code for all 3 cases (local, AliEn, CCDB),
      // that common code is further below.

    } else {

      // File you want to access is on your local computer.
      // In this branch, there is a special code to support file access in that case.

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
          LOGF(fatal, "\033[1;31m%s at line %d : this crash can happen if in the output file there is no list with weights for the current run number = %s\033[0m", __FUNCTION__, __LINE__, runNumber); // changed runNumber instead of tc.fRunNumber.Data()
        }
      }

      // OK, we got the desired run-specific TList with efficiency corrections for local case,
      // after that we can use the common code for all 3 cases (local, AliEn, CCDB),
      // that common code is further below.

    } // end of else

    // The common code for all three cases, where from "listWithRuns" you fetch the desired histogram with efficiency corrections:
    // listWithRuns->ls();

    // get pointer to "hist" from "listWithRuns"
    hist = reinterpret_cast<TH1F*>(getObjectFromList(listWithRuns, "histWithEfficiencyCorrections_1"));
    if (!hist) {
      LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
    }

    hist->SetDirectory(0);                                    // remove ownerhip from an external file
    TH1F* histClone = reinterpret_cast<TH1F*>(hist->Clone()); // yes, I have to clone here

    if (!bFileIsInCCDB) { // not for CCDB, there the CCDB manager owns the object and deletes it itself
      delete baseList;    // release back the memory
    }

    return histClone;

  } // end of TH1F* getHistogramWithWeights(...)

  //__ EVENT CUTS __
  template <ERecSim rs, typename T1>
  bool eventCuts(T1 const& collision)
  {

    // *) Vertex Z cut
    if (tc.fUseVertexZCut) {
      if (!(tc.fVertexZ[0] <= collision.posZ() && collision.posZ() < tc.fVertexZ[1])) {
        return false; // if at least one event cut is not satisfied
      }
    }
    return true; // all event cuts passed
  }

  //__ PARTICLE CUTS __
  template <ERecSim rs, typename T>
  bool particleCuts(T const& track)
  {

    // *) pt cut
    if (tc.fUsePtCut) {
      if (!(tc.fPt[0] <= track.pt() && track.pt() < tc.fPt[1])) {
        return false; // if at least one event cut is not satisfied
      }
    }

    return true; // all particle cuts passed
  }

  //__ Fill all event histograms __ #ba = EBefore or EAfter
  template <ERecSim rs, typename T1, typename T2>
  void fillEventHistograms(T1 const& collision, T2 const& tracks, EBeforeAfter ba)
  {

    if constexpr (rs == ERec || rs == ERecAndSim) {

      float centRec = getCentrality(collision);

      ec.fEventHistograms[EMultiplicity][ERec][ba]->Fill(tracks.size()); // defined as total number of available particles
      ec.fEventHistograms[EReferenceMultiplicity][ERec][ba]->Fill(getReferenceMultiplicity(collision));
      ec.fEventHistograms[ECentrality][ERec][ba]->Fill(centRec);
      ec.fEventHistograms[EVertexX][ERec][ba]->Fill(collision.posX());
      ec.fEventHistograms[EVertexY][ERec][ba]->Fill(collision.posY());
      ec.fEventHistograms[EVertexZ][ERec][ba]->Fill(collision.posZ());

      if constexpr (rs == ERecAndSim) {
        if (!collision.has_mcCollision()) {
          LOGF(warning, "  No MC collision for this collision, skip simulated part...");
          return;
        }
        auto mcCollision = collision.mcCollision();

        ec.fEventHistograms[EMultiplicity][ESim][ba]->Fill(tracks.size());
        ec.fEventHistograms[EVertexX][ESim][ba]->Fill(mcCollision.posX());
        ec.fEventHistograms[EVertexY][ESim][ba]->Fill(mcCollision.posY());
        ec.fEventHistograms[EVertexZ][ESim][ba]->Fill(mcCollision.posZ());

        // impact parameter and centrality calculated from MC
        float b = mcCollision.impactParameter();
        if (b > 0.) { // only if the impact parameter is stored in the MC
          float centSim = centralityFromImpactParameter(b);
          ec.fEventHistograms[EImpactParameter][ESim][ba]->Fill(b);
          ec.fEventHistograms[ECentrality][ESim][ba]->Fill(centSim);

          // QA centrality sim vs centrality rec after the cuts
          if (tc.fFillQA && ba == EAfter) {
            qa.fHistCentralityRecVsSim->Fill(centSim, centRec);
          } // if (tc.fFillQA && ba == EAfter)

        } // if (b > 0.)

      } // if constexpr (rs == ERecAndSim)

    } // if constexpr (rs == ERec || rs == ERecAndSim)

  } // fillEventHistograms

  template <ERecSim rs, typename T1, typename T2>
  void steer(T1 const& collision, T2 const& tracks)
  {

    // Dry run:
    if (tc.fDryRun) {
      return;
    }

    // Load the weights once, at the first collision, because the run number is only known here
    if (!tc.fWeightsLoaded && !cfFileWithWeights.value.empty()) {
      histWeights = getHistogramWithWeights(cfFileWithWeights.value.c_str(), Form("%d", collision.bc().runNumber()));
      histWeights->SetTitle(cfFileWithWeights.value.c_str()); // book-keeping: title = file path
      // value takes the text as std::string and c_str makes it into a constant char*

      fBaseList->Add(histWeights); // so it is saved in AnalysisResults.root
      tc.fWeightsLoaded = true;
    }

    // Print current run number:
    LOGF(info, "Run number: %d", collision.bc().runNumber());

    // Print centrality estimated with "FT0M" estimator:
    LOGF(info, "Centrality: %f", getCentrality(collision));

    // Print vertex X position:
    LOGF(info, "Vertex X position: %f", collision.posX());

    // Fill event histograms before event cuts
    fillEventHistograms<rs>(collision, tracks, EBefore);

    // Event cuts: leave this collision if it doesn't pass
    if (!eventCuts<rs>(collision)) {
      return;
    }

    // Main loop over particles:
    auto track = tracks.iteratorAt(0); // set the type and scope from one instance
    for (int64_t i = 0; i < tracks.size(); i++) {

      track = tracks.iteratorAt(i);

      // Print track azimuthal angle:
      // LOGF(info, "Track azimuthal angle: %f", track.phi());

      // Fill reconstructed ...:
      if constexpr (rs == ERec || rs == ERecAndSim) {
        // Fill tracks before cuts:

        pc.fParticleHistograms[EPt][ERec][EBefore]->Fill(track.pt());
        pc.fParticleHistograms[EPhi][ERec][EBefore]->Fill(track.phi());

        // ... and corresponding MC truth simulated (before cuts):
        // See https://github.com/AliceO2Group/O2Physics/blob/master/Tutorials/src/mcHistograms.cxx
        // See https://aliceo2group.github.io/analysis-framework/docs/datamodel/ao2dTables.html#montecarlo

        if constexpr (rs == ERecAndSim) {
          if (track.has_mcParticle()) {
            auto mcparticle = track.mcParticle(); // corresponding MC truth simulated particle
            pc.fParticleHistograms[EPt][ESim][EBefore]->Fill(mcparticle.pt());
            pc.fParticleHistograms[EPhi][ESim][EBefore]->Fill(mcparticle.phi());
          } else {
            LOGF(warning, "  No MC particle for this track...");
          }
        }

        // Particle cuts: skip this track, not the whole event
        if (!particleCuts<rs>(track)) {
          continue;
        }

        // Fill tracks after cuts:

        pc.fParticleHistograms[EPt][ERec][EAfter]->Fill(track.pt());
        pc.fParticleHistograms[EPhi][ERec][EAfter]->Fill(track.phi());

        if constexpr (rs == ERecAndSim) {
          if (!track.has_mcParticle()) {
            continue;
          }

          auto mcparticle = track.mcParticle(); // corresponding MC truth simulated particles
          pc.fParticleHistograms[EPt][ESim][EAfter]->Fill(mcparticle.pt());
          pc.fParticleHistograms[EPhi][ESim][EAfter]->Fill(mcparticle.phi());

        } // end of if constexpr (rs == ERecAndSim) {

      } // if constexpr (rs == ERec || rs == ERecAndSim) {

    } // end of for (int64_t i = 0; i < tracks.size(); i++) {

    // Fill event histograms after event and particle cuts
    fillEventHistograms<rs>(collision, tracks, EAfter);

  } // end of template <ERecSim rs, typename T1, typename T2> void steer(T1 const& collision, T2 const& tracks) {

  // *) Initialize and book all objects:
  void init(InitContext&)
  {

    // ... code to book and initialize all analysis objects ...

    // *) Set automatically what to process, from an implicit variable "doprocessSomeProcessName" within a PROCESS_SWITCH clause:
    tc.fProcess[EProcessRec] = doprocessRec;
    tc.fProcess[EProcessRecSim] = doprocessRecSim;
    tc.fProcess[EProcessSim] = doprocessSim;

    // *) Configure your task using configurables in the json file:
    tc.fDryRun = cfDryRun;
    tc.fFillQA = cfFillQA;
    tc.fSigmaInel = cfSigmaInel;

    // vertex Z cut copy the swtich and range
    tc.fUseVertexZCut = cfUseVertexZCut;
    std::vector<float> lVzCut = cfVertexZ.value;

    tc.fVertexZ[0] = lVzCut[0];
    tc.fVertexZ[1] = lVzCut[1];

    // pt cut copy switch and range
    tc.fUsePtCut = cfUsePtCut;
    std::vector<float> lPtCut = cfPt.value;

    tc.fPt[0] = lPtCut[0];
    tc.fPt[1] = lPtCut[1];

    // translate the centrality estimator string into the enum
    if (cfCentralityEstimator.value == "FT0C") {
      tc.fCentralityEstimator = ECentFT0C;
    } else if (cfCentralityEstimator.value == "FT0M") {
      tc.fCentralityEstimator = ECentFT0M;
    } else if (cfCentralityEstimator.value == "FT0A") {
      tc.fCentralityEstimator = ECentFT0A;
    } else {
      LOGF(fatal, "Unknown centrality estimator in json: %s", cfCentralityEstimator.value.c_str());
    }

    // same for the reference multiplicity estimator
    if (cfReferenceMultiplicityEstimator.value == "TPC") {
      tc.fReferenceMultiplicityEstimator = ERefMultTPC;
    } else if (cfReferenceMultiplicityEstimator.value == "FT0C") {
      tc.fReferenceMultiplicityEstimator = ERefMultFT0C;
    } else if (cfReferenceMultiplicityEstimator.value == "FT0M") {
      tc.fReferenceMultiplicityEstimator = ERefMultFT0M;
    } else {
      LOGF(fatal, "Unknown reference multiplicity estimator in json: %s", cfReferenceMultiplicityEstimator.value.c_str());
    }

    // *) Book base list:
    TList* temp = new TList();
    temp->SetOwner(true);
    fBaseList.setObject(temp);

    // *) Book and nest all other TLists:
    pc.fParticleHistogramsList = new TList();
    pc.fParticleHistogramsList->SetName("ParticleHistograms");
    pc.fParticleHistogramsList->SetOwner(true);
    fBaseList->Add(pc.fParticleHistogramsList); // any nested TList in the base TList appears as a subdir in the output ROOT file

    ec.fEventHistogramsList = new TList();
    ec.fEventHistogramsList->SetName("EventHistograms");
    ec.fEventHistogramsList->SetOwner(true);
    fBaseList->Add(ec.fEventHistogramsList); // any nested TList in the base TList appears as a subdir in the output ROOT file

    if (tc.fFillQA) { // only if QA switched on
      qa.fQAHistogramsList = new TList();
      qa.fQAHistogramsList->SetName("QA");
      qa.fQAHistogramsList->SetOwner(true);
      fBaseList->Add(qa.fQAHistogramsList);
    }

    // *) Book (pt distribution and azimuthal-angle histograms) + event quantities with binning defined through configurables in the json file:
    std::vector<float> lPtBins = cfPtBins.value; // define local array and initialize it from an array set in the configurables
    int nBins = static_cast<int>(lPtBins[0]);
    float min = lPtBins[1];
    float max = lPtBins[2];

    std::vector<float> lPhiBins = cfPhiBins.value; // azimuthal angle
    int nPhiBins = static_cast<int>(lPhiBins[0]);
    float phiMin = lPhiBins[1];
    float phiMax = lPhiBins[2];

    std::vector<float> lMultBins = cfMultiplicityBins.value; // multiplicity
    int nMultBins = static_cast<int>(lMultBins[0]);
    float multMin = lMultBins[1];
    float multMax = lMultBins[2];

    std::vector<float> lRefmultBins = cfRefmultiplicityBins.value; // reference multiplicity
    int nRefMultBins = static_cast<int>(lRefmultBins[0]);
    float refMultMin = lRefmultBins[1];
    float refMultMax = lRefmultBins[2];

    std::vector<float> lCentBins = cfCentralityBins.value; // centrality
    int nCentBins = static_cast<int>(lCentBins[0]);
    float centMin = lCentBins[1];
    float centMax = lCentBins[2];

    std::vector<float> lBBins = cfImpactparameterBins.value; // impact parameter
    int nBBins = static_cast<int>(lBBins[0]);
    float bMin = lBBins[1];
    float bMax = lBBins[2];

    std::vector<float> lVxBins = cfVertexXBins.value; // vertex X
    int nvxBins = static_cast<int>(lVxBins[0]);
    float vxMin = lVxBins[1];
    float vxMax = lVxBins[2];

    std::vector<float> lVyBins = cfVertexYBins.value; // vertex y
    int nvyBins = static_cast<int>(lVyBins[0]);
    float vyMin = lVyBins[1];
    float vyMax = lVyBins[2];

    std::vector<float> lVzBins = cfVertexZBins.value; // vertex z
    int nvzBins = static_cast<int>(lVzBins[0]);
    float vzMin = lVzBins[1];
    float vzMax = lVzBins[2];

    // histogram booked before and after cuts
    const char* baNames[EBeforeAfter_N] = {"before", "after"};
    for (int ba = 0; ba < EBeforeAfter_N; ba++) {

      // par
      pc.fParticleHistograms[EPt][ERec][ba] = new TH1F(Form("fHistPt[ERec][%s]", baNames[ba]), Form("pt distribution for reconstructed particles (%s cuts)", baNames[ba]), nBins, min, max);
      pc.fParticleHistograms[EPt][ERec][ba]->GetXaxis()->SetTitle("p_{T}");
      pc.fParticleHistogramsList->Add(pc.fParticleHistograms[EPt][ERec][ba]);

      pc.fParticleHistograms[EPt][ESim][ba] = new TH1F(Form("fHistPt[ESim][%s]", baNames[ba]), Form("pt distribution for simulated particles (%s cuts)", baNames[ba]), nBins, min, max);
      pc.fParticleHistograms[EPt][ESim][ba]->GetXaxis()->SetTitle("p_{T}");
      pc.fParticleHistogramsList->Add(pc.fParticleHistograms[EPt][ESim][ba]);

      pc.fParticleHistograms[EPhi][ERec][ba] = new TH1F(Form("fHistPhi[ERec][%s]", baNames[ba]), Form("Azimuthal-angle distribution for reconstructed particles (%s cuts)", baNames[ba]), nPhiBins, phiMin, phiMax);
      pc.fParticleHistograms[EPhi][ERec][ba]->GetXaxis()->SetTitle("#varphi (rad)");
      pc.fParticleHistogramsList->Add(pc.fParticleHistograms[EPhi][ERec][ba]);

      pc.fParticleHistograms[EPhi][ESim][ba] = new TH1F(Form("fHistPhi[ESim][%s]", baNames[ba]), Form("Azimuthal-angle distribution for simulated particles (%s cuts)", baNames[ba]), nPhiBins, phiMin, phiMax);
      pc.fParticleHistograms[EPhi][ESim][ba]->GetXaxis()->SetTitle("#varphi (rad)");
      pc.fParticleHistogramsList->Add(pc.fParticleHistograms[EPhi][ESim][ba]);

      // ev
      ec.fEventHistograms[EMultiplicity][ERec][ba] = new TH1F(Form("fHistMultiplicity[ERec][%s]", baNames[ba]), Form("Multiplicity for reconstructed particles (%s cuts)", baNames[ba]), nMultBins, multMin, multMax);
      ec.fEventHistograms[EMultiplicity][ERec][ba]->GetXaxis()->SetTitle("total number of available particles");
      ec.fEventHistogramsList->Add(ec.fEventHistograms[EMultiplicity][ERec][ba]);

      ec.fEventHistograms[EMultiplicity][ESim][ba] = new TH1F(Form("fHistMultiplicity[ESim][%s]", baNames[ba]), Form("Multiplicity for simulated particles (%s cuts)", baNames[ba]), nMultBins, multMin, multMax);
      ec.fEventHistograms[EMultiplicity][ESim][ba]->GetXaxis()->SetTitle("total number of available particles");
      ec.fEventHistogramsList->Add(ec.fEventHistograms[EMultiplicity][ESim][ba]);

      // reference multiplicity, only reconstructed
      ec.fEventHistograms[EReferenceMultiplicity][ERec][ba] = new TH1F(Form("fHistReferenceMultiplicity[ERec][%s]", baNames[ba]), Form("Reference multiplicity, reconstructed (%s cuts)", baNames[ba]), nRefMultBins, refMultMin, refMultMax);
      ec.fEventHistograms[EReferenceMultiplicity][ERec][ba]->GetXaxis()->SetTitle("reference multiplicity");
      ec.fEventHistogramsList->Add(ec.fEventHistograms[EReferenceMultiplicity][ERec][ba]);

      ec.fEventHistograms[ECentrality][ERec][ba] = new TH1F(Form("fHistCentrality[ERec][%s]", baNames[ba]), Form("Centrality for reconstructed particles (%s cuts)", baNames[ba]), nCentBins, centMin, centMax);
      ec.fEventHistograms[ECentrality][ERec][ba]->GetXaxis()->SetTitle("centrality (%)");
      ec.fEventHistogramsList->Add(ec.fEventHistograms[ECentrality][ERec][ba]);

      // sim from impact parameter
      ec.fEventHistograms[ECentrality][ESim][ba] = new TH1F(Form("fHistCentrality[ESim][%s]", baNames[ba]), Form("Centrality for simulated particles, from impact parameter (%s cuts)", baNames[ba]), nCentBins, centMin, centMax);
      ec.fEventHistograms[ECentrality][ESim][ba]->GetXaxis()->SetTitle("centrality (%)");
      ec.fEventHistogramsList->Add(ec.fEventHistograms[ECentrality][ESim][ba]);

      // impact parameter, only simulated
      ec.fEventHistograms[EImpactParameter][ESim][ba] = new TH1F(Form("fHistImpactParameter[ESim][%s]", baNames[ba]), Form("Impact parameter, simulated (%s cuts)", baNames[ba]), nBBins, bMin, bMax);
      ec.fEventHistograms[EImpactParameter][ESim][ba]->GetXaxis()->SetTitle("b (fm)");
      ec.fEventHistogramsList->Add(ec.fEventHistograms[EImpactParameter][ESim][ba]);

      ec.fEventHistograms[EVertexX][ERec][ba] = new TH1F(Form("fHistVertexX[ERec][%s]", baNames[ba]), Form("Vertex X for reconstructed particles (%s cuts)", baNames[ba]), nvxBins, vxMin, vxMax);
      ec.fEventHistograms[EVertexX][ERec][ba]->GetXaxis()->SetTitle("vertex X");
      ec.fEventHistogramsList->Add(ec.fEventHistograms[EVertexX][ERec][ba]);

      ec.fEventHistograms[EVertexX][ESim][ba] = new TH1F(Form("fHistVertexX[ESim][%s]", baNames[ba]), Form("Vertex X for simulated particles (%s cuts)", baNames[ba]), nvxBins, vxMin, vxMax);
      ec.fEventHistograms[EVertexX][ESim][ba]->GetXaxis()->SetTitle("vertex X");
      ec.fEventHistogramsList->Add(ec.fEventHistograms[EVertexX][ESim][ba]);

      ec.fEventHistograms[EVertexY][ERec][ba] = new TH1F(Form("fHistVertexY[ERec][%s]", baNames[ba]), Form("Vertex Y for reconstructed particles (%s cuts)", baNames[ba]), nvyBins, vyMin, vyMax);
      ec.fEventHistograms[EVertexY][ERec][ba]->GetXaxis()->SetTitle("vertex Y");
      ec.fEventHistogramsList->Add(ec.fEventHistograms[EVertexY][ERec][ba]);

      ec.fEventHistograms[EVertexY][ESim][ba] = new TH1F(Form("fHistVertexY[ESim][%s]", baNames[ba]), Form("Vertex Y for simulated particles (%s cuts)", baNames[ba]), nvyBins, vyMin, vyMax);
      ec.fEventHistograms[EVertexY][ESim][ba]->GetXaxis()->SetTitle("vertex Y");
      ec.fEventHistogramsList->Add(ec.fEventHistograms[EVertexY][ESim][ba]);

      ec.fEventHistograms[EVertexZ][ERec][ba] = new TH1F(Form("fHistVertexZ[ERec][%s]", baNames[ba]), Form("Vertex Z for reconstructed particles (%s cuts)", baNames[ba]), nvzBins, vzMin, vzMax);
      ec.fEventHistograms[EVertexZ][ERec][ba]->GetXaxis()->SetTitle("vertex Z");
      ec.fEventHistogramsList->Add(ec.fEventHistograms[EVertexZ][ERec][ba]);

      ec.fEventHistograms[EVertexZ][ESim][ba] = new TH1F(Form("fHistVertexZ[ESim][%s]", baNames[ba]), Form("Vertex Z for simulated particles (%s cuts)", baNames[ba]), nvzBins, vzMin, vzMax);
      ec.fEventHistograms[EVertexZ][ESim][ba]->GetXaxis()->SetTitle("vertex Z");
      ec.fEventHistogramsList->Add(ec.fEventHistograms[EVertexZ][ESim][ba]);

    } // end of for (int ba = 0; ba < 2; ba++)

    // QA histograms, booked only if cfFillQA is true
    if (tc.fFillQA) {
      qa.fHistCentralityRecVsSim = new TH2F("fHistCentralityRecVsSim", "Centrality: reconstructed vs. simulated", nCentBins, centMin, centMax, nCentBins, centMin, centMax);
      qa.fHistCentralityRecVsSim->GetXaxis()->SetTitle("centrality simulated (%)");
      qa.fHistCentralityRecVsSim->GetYaxis()->SetTitle("centrality reconstructed (%)");
      qa.fQAHistogramsList->Add(qa.fHistCentralityRecVsSim);
    }

  } // end of void init(InitContext&) {

  // A) Process only reconstructed data:
  void processRec(CollisionRec const& collision, aod::BCs const&, TracksRec const& tracks)
  {
    // ...

    // *) steer all analysis steps:
    steer<ERec>(collision, tracks);
  }
  PROCESS_SWITCH(StudentO2Physics, processRec, "process only reconstructed data", true); // yes, keep always one process switch "true", so that there is default running version

  // -------------------------------------------

  // B) Process both reconstructed and corresponding MC truth simulated data:
  void processRecSim(CollisionRecSim const& collision, aod::BCs const&, TracksRecSim const& tracks, aod::McParticles const&, aod::McCollisions const&)
  {
    steer<ERecAndSim>(collision, tracks);
  }
  PROCESS_SWITCH(StudentO2Physics, processRecSim, "process both reconstructed and corresponding MC truth simulated data", false);

  // -------------------------------------------

  // C) Process only simulated data:
  void processSim(CollisionSim const& /*collision*/, aod::BCs const&, TracksSim const& /*tracks*/)
  {
    // steer<ESim>(collision, tracks); // TBI 20241105 not ready yet, but I do not really need this one urgently, since RecSim is working, and I need that one for efficiencies...
  }
  PROCESS_SWITCH(StudentO2Physics, processSim, "process only simulated data", false);

}; // struct StudentO2Physics {

// *) The final touch:
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<StudentO2Physics>(cfgc),
  };
} // WorkflowSpec...
