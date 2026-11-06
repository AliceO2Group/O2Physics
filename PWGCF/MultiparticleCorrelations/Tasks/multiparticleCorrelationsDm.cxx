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
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

// *) ROOT and C++ include's:

#include <TCollection.h>
#include <TFile.h>
#include <TGrid.h>
#include <TH1.h>
#include <TH2.h>
#include <TList.h>
#include <TObject.h>
#include <TProfile.h>
#include <TString.h>
#include <TSystem.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <string>
#include <vector>

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
struct MultiparticleCorrelationsDm { // this name is used in lower-case format to name the TDirectoryFile in AnalysisResults.root

  // *) Base TList to hold all output objects:
  TString sBaseListName = "Default list name"; // yes, I declare it separately, because I need it also later in BailOut() function
  OutputObj<TList> fBaseList{sBaseListName.Data(),
                             OutputObjHandlingPolicy::AnalysisObject,
                             OutputObjSourceType::OutputObjSource};

  // *) CCDB:
  // Service<ccdb::BasicCCDBManager> ccdb; // support for offline callibration data base, not needed for the time being...
  // to be able to use ccdb object in the analysis task
  Service<ccdb::BasicCCDBManager> ccdb{};

  // *) Define configurables:

  Configurable<bool> cfDryRun{"cfDryRun", false, "book all histos and run without filling and calculating anything"};
  Configurable<bool> cfFillQA{"cfFillQA", false, "book and fill QA histograms (only for development and validation)"}; // switch for QA histograms

  // Define cuts:
  Configurable<bool> cfUseVertexZCut{"cfUseVertexZCut", true, "apply the vertex Z event cut"};
  Configurable<bool> cfUsePtCut{"cfUsePtCut", true, "apply the pt particle cut"};
  Configurable<bool> cfUseTimeRange{"cfUseTimeRange", true, "rejects events with overlapping collisions"};
  Configurable<bool> cfUseRof{"cfUseRof", true, "removes events containing multiple collisions in the same readout frame"};
  Configurable<bool> cfUsePileUp{"cfUsePileUp", true, "apply pile up cuts"};
  Configurable<bool> cfUseITSTPC{"cfUseITSTPC", true, "use reconstructed from the ITS (Inner Tracking System) and TPC (Time Porjection Chamber) tracks"};
  Configurable<bool> cfUseITSLayers{"cfUseITSLayers", true, "use the ITS layers cut"};
  Configurable<bool> cfUseZFT0vsPV{"cfUseZFT0vsPV", true, "use this cut"};
  Configurable<bool> cfUsePrevRof{"cfUsePrevRof", true, " use this cut"};

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
  Configurable<std::vector<float>> cfNContribBins{"cfNContribBins", {500, 0., 5000.}, "nContribBins, contribMin, contribMax"};
  Configurable<std::vector<float>> cfVertexZ{"cfVertexZ", {-10., 10.}, "{min, max} [cm]"};
  Configurable<std::vector<float>> cfPt{"cfPt", {0.2, 5.0}, "{min, max} [GeV]"};
  Configurable<bool> cfUseEtaCut{"cfUseEtaCut", true, "apply the eta particle cut (pseudorapidity"};
  Configurable<std::vector<float>> cfEta{"cfEta", {-0.8, 0.8}, "eta range: {min, max}, with convention: min <= eta < max"};

  // for weights (efficiency corrections)
  Configurable<bool> cfUsePhiWeights{"cfUsePhiWeights", false, "use phi weights from the file in cfFileWithWeights"};
  Configurable<bool> cfUsePtWeights{"cfUsePtWeights", false, "use pt weights from the file in cfFileWithWeights"};
  Configurable<bool> cfUseEtaWeights{"cfUseEtaWeights", false, "use eta weights from the file in cfFileWithWeights"};
  Configurable<std::string> cfFileWithWeights{"cfFileWithWeights", "", "path to external ROOT file with weights; on Hyperloop only CCDB paths (/alice-ccdb.cern.ch/...) work"};

  Configurable<bool> cfVerbose{"cfVerbose", false, "print run number, centrality and vertex for each collision"};

  // *) Define and initialize all data members to be called in the main process* functions:
  // **) Task configuration:
  struct TaskConfiguration {
    std::array<bool, EProcess_N> fProcess{}; // see enum EProcess
    bool fDryRun = false;                    // book all histos and run without filling and calculating anything
    bool fFillQA = false;                    // book and fill QA histos
    bool fUseVertexZCut = true;
    bool fUsePtCut = true;
    bool fUseTimeRange = true;
    bool fUseRof = true;
    bool fUsePileUp = true;
    bool fUseITSTPC = true;
    bool fUseITSLayers = true;
    bool fUseZFT0vsPV = true;
    bool fUsePrevRof = true;

    int fCentralityEstimator = ECentFT0C;
    int fReferenceMultiplicityEstimator = ERefMultTPC;

    float fSigmaInel = 771.;                     // sigma inel [fm^2]
    std::array<float, 2> fVertexZ = {-10., 10.}; // {min, max} [cm]
    std::array<float, 2> fPt = {0.2, 5.0};       // {min, max} [GeV]
    bool fUseEtaCut = true;
    std::array<float, 2> fEta = {-0.8, 0.8}; // {min, max}
    int fCurrentRunNumber = -1;              // run number of the weights currently in use
    bool fUsePhiWeights = false;             // load phi weights
    bool fUsePtWeights = false;              // load pt weights
    bool fUseEtaWeights = false;             // load eta weights
    bool fVerbose = false;                   // print the values / collision

  } tc; // you have to prepend "tc." for all objects name in this group later in the code

  // **) Particle histograms:
  struct ParticleHistograms {
    TList* fParticleHistogramsList = nullptr; //!<! list to hold all control particle histograms

    std::array<std::array<std::array<TH1F*, EBeforeAfter_N>, 2>, EParticleHistograms_N> fParticleHistograms{}; // [type][rec, sim][before, after]

  } pc; // you have to prepend "pc." for all objects name in this group later in the code

  // **) Event histograms:
  struct EventHistograms {
    TList* fEventHistogramsList = nullptr;

    std::array<std::array<std::array<TH1F*, EBeforeAfter_N>, 2>, EEventHistograms_N> fEventHistograms{}; // [type][rec, sim][before, after]

  } ec; // you have to prepend "ec." for all objects name in this group later in the code

  // **) QA histograms (booked and filled only if cfFillQA is true):

  struct QAHistograms {
    TList* fQAHistogramsList = nullptr;
    TH2F* fHistCentralityRecVsSim = nullptr;                    // x = centrality sim (from impact parameter), y = centrality rec
    std::array<TH2F*, EBeforeAfter_N> fHistRefMultVsNContrib{}; // x = reference multiplicity, y = number of vertex contributors [before, after]
  } qa;                                                         // you have to prepend "qa." for all objects name in this group later in the code

  // **) Correlations:
  struct Correlations {
    TList* fCorrelationsList = nullptr;
    TProfile* fProfTwo = nullptr; // <<2>> = v2{2}^2 vs. centrality (rec)
  } cc;

  // **) Histogram with weights from an external file
  TH1F* histPhiWeights = nullptr; // phi weights for the current run
  TH1F* histPtWeights = nullptr;  // pt weights for the current run
  TH1F* histEtaWeights = nullptr; // eta weights for the current run

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
      return nullptr;
    }

    // The object is in the current base list:
    TObject* objectFinal = list->FindObject(objectName); // the final object I am after
    if (objectFinal) {
      return objectFinal;
    }

    // Otherwise, search for the object recursively in the nested lists:
    TObject* objectIter = nullptr; // iterator object in the loop below
    TIter next(list);
    while ((objectIter = next())) // double round braces are to silence the warnings
    {
      if (TString(objectIter->ClassName()).EqualTo("TList")) {
        objectFinal = getObjectFromList(dynamic_cast<TList*>(objectIter), objectName);
        if (objectFinal) {
          return objectFinal;
        }
      }
    } // while(objectIter = next())

    return nullptr;

  } // TObject* getObjectFromList(...)
  TH1F* getHistogramWithWeights(const char* filePath, const char* runNumber, const char* histName)
  {
    // Last update: 20261018

    // *) Return value:
    TH1F* hist = nullptr;
    TList* baseList = nullptr;     // base top-level list in the TFile, e.g. named "ccdb_object"
    TList* listWithRuns = nullptr; // nested list with run-wise TList's holding run-specific weights

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
      listWithRuns = dynamic_cast<TList*>(getObjectFromList(baseList, runNumber));
      if (!listWithRuns) {
        TString runNumberWithLeadingZeroes = "000";
        runNumberWithLeadingZeroes += runNumber; // another try, with "000" prepended to run number
        listWithRuns = dynamic_cast<TList*>(getObjectFromList(baseList, runNumberWithLeadingZeroes.Data()));
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
      baseList = dynamic_cast<TList*>(ccdb->get<TList>(TString(filePath).ReplaceAll("/alice-ccdb.cern.ch/", "").Data()));

      if (!baseList) {
        LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
      }

      listWithRuns = dynamic_cast<TList*>(getObjectFromList(baseList, runNumber));
      if (!listWithRuns) {
        TString runNumberWithLeadingZeroes = "000";
        runNumberWithLeadingZeroes += runNumber; // another try, with "000" prepended to run number
        listWithRuns = dynamic_cast<TList*>(getObjectFromList(baseList, runNumberWithLeadingZeroes.Data()));
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

      listWithRuns = dynamic_cast<TList*>(getObjectFromList(baseList, runNumber));
      if (!listWithRuns) {
        TString runNumberWithLeadingZeroes = "000";
        runNumberWithLeadingZeroes += runNumber; // another try, with "000" prepended to run number
        listWithRuns = dynamic_cast<TList*>(getObjectFromList(baseList, runNumberWithLeadingZeroes.Data()));
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
    hist = dynamic_cast<TH1F*>(getObjectFromList(listWithRuns, histName));
    if (!hist) {
      LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
    }

    hist->SetDirectory(nullptr);                          // remove ownerhip from an external file
    TH1F* histClone = dynamic_cast<TH1F*>(hist->Clone()); // yes, I have to clone here

    if (!bFileIsInCCDB) { // not for CCDB, there the CCDB manager owns the object and deletes it itself
      delete baseList;    // release back the memory
    }

    return histClone;

  } // end of TH1F* getHistogramWithWeights(...)

  // load one kind of weights for this run (or reuse it if already loaded), and save it in the output
  TH1F* loadWeights(const char* nameInFile, const char* nameInOutput, int runNumber)
  {
    TString name = Form("%s_%d", nameInOutput, runNumber);
    TH1F* hist = dynamic_cast<TH1F*>(fBaseList->FindObject(name.Data()));
    if (!hist) {
      hist = getHistogramWithWeights(cfFileWithWeights.value.c_str(), Form("%d", runNumber), nameInFile);
      hist->SetName(name.Data());                      // e.g. histPhiWeights_544116
      hist->SetTitle(cfFileWithWeights.value.c_str()); // book-keeping: title = file path
      fBaseList->Add(hist);                            // saved in AnalysisResults.root
    }
    return hist;
  }

  // weight of one particle: product of all enabled weights, 1 if none is enabled
  float getParticleWeight(float phi, float pt, float eta)
  {
    float weight = 1.;
    if (tc.fUsePhiWeights && histPhiWeights) {
      weight *= histPhiWeights->GetBinContent(histPhiWeights->FindBin(phi));
    }
    if (tc.fUsePtWeights && histPtWeights) {
      weight *= histPtWeights->GetBinContent(histPtWeights->FindBin(pt));
    }
    if (tc.fUseEtaWeights && histEtaWeights) {
      weight *= histEtaWeights->GetBinContent(histEtaWeights->FindBin(eta));
    }
    return weight;
  }

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
    // *) No collision in time range cut
    if (tc.fUseTimeRange) {
      if (!collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
        return false; // if at least one event cut is not satisfied
      }
    }

    // *) No collision in the same readout frame cut
    if (tc.fUseRof) {
      if (!collision.selection_bit(o2::aod::evsel::kNoCollInRofStandard)) {
        return false; // if at least one event cut is not satisfied
      }
    }

    // *) No same-bunch pile-up cut
    if (tc.fUsePileUp) {
      if (!collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
        return false; // if at least one event cut is not satisfied
      }
    }

    // *) Vertex from ITS and TPC tracks cut
    if (tc.fUseITSTPC) {
      if (!collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) {
        return false; // if at least one event cut is not satisfied
      }
    }

    // *) All ITS layers good cut
    if (tc.fUseITSLayers) {
      if (!collision.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll)) {
        return false; // if at least one event cut is not satisfied
      }
    }

    // *) z vertex from FT0 and from tracks agree cut
    if (tc.fUseZFT0vsPV) {
      if (!collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
        return false; // if at least one event cut is not satisfied
      }
    }

    // *) No high-multiplicity collision in previous readout frame cut
    if (tc.fUsePrevRof) {
      if (!collision.selection_bit(o2::aod::evsel::kNoHighMultCollInPrevRof)) {
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

    // *) eta cut
    if (tc.fUseEtaCut) {
      if (!(tc.fEta[0] <= track.eta() && track.eta() < tc.fEta[1])) {
        return false; // if at least one particle cut is not satisfied
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

      // QA: reference multiplicity vs. number of contributors
      if (tc.fFillQA) {
        qa.fHistRefMultVsNContrib[ba]->Fill(getReferenceMultiplicity(collision), collision.numContrib());
      }

      if constexpr (rs == ERecAndSim) {
        if (!collision.has_mcCollision()) {
          LOGF(debug, "  No MC collision for this collision, skip simulated part...");
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

    // weights are run-specific: load them again when the run number changes
    int runNumber = collision.bc().runNumber();
    if (runNumber != tc.fCurrentRunNumber) {
      if (tc.fUsePhiWeights) {
        histPhiWeights = loadWeights("histWithEfficiencyCorrections_1", "histPhiWeights", runNumber);
      }
      if (tc.fUsePtWeights) {
        histPtWeights = loadWeights("histWithEfficiencyCorrections_2", "histPtWeights", runNumber);
      }
      if (tc.fUseEtaWeights) {
        histEtaWeights = loadWeights("histWithEfficiencyCorrections_3", "histEtaWeights", runNumber);
      }
      tc.fCurrentRunNumber = runNumber;
    }

    // Print info for each collision, only if switched on in the config:
    if (tc.fVerbose) {
      LOGF(info, "Run number: %d", collision.bc().runNumber());
      LOGF(info, "Centrality: %f", getCentrality(collision));
      LOGF(info, "Vertex X position: %f", collision.posX());
    }

    // Fill event histograms before event cuts
    fillEventHistograms<rs>(collision, tracks, EBefore);

    // Event cuts: leave this collision if it doesn't pass
    if (!eventCuts<rs>(collision)) {
      return;
    }

    // Main loop over particles:
    // Q-vector of harmonic 2 and sums of weights, for <2> of this event
    double q2Re = 0., q2Im = 0., sumW = 0., sumW2 = 0.;

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
            LOGF(debug, "  No MC particle for this track...");
          }
        }

        // Particle cuts: skip this track, not the whole event
        if (!particleCuts<rs>(track)) {
          continue;
        }

        // Fill tracks after cuts:

        pc.fParticleHistograms[EPt][ERec][EAfter]->Fill(track.pt());
        pc.fParticleHistograms[EPhi][ERec][EAfter]->Fill(track.phi());

        // add this particle to the Q-vector
        const double w = getParticleWeight(track.phi(), track.pt(), track.eta());
        q2Re += w * std::cos(2. * track.phi());
        q2Im += w * std::sin(2. * track.phi());
        sumW += w;
        sumW2 += w * w;

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
    // <2> of this event, filled vs. centrality with event weight sumW^2 - sumW2 (= M(M-1) for unit weights)
    if constexpr (rs == ERec || rs == ERecAndSim) {
      const double denominator = sumW * sumW - sumW2;
      if (denominator > 0.) {
        const double two = (q2Re * q2Re + q2Im * q2Im - sumW2) / denominator;
        cc.fProfTwo->Fill(getCentrality(collision), two, denominator);
      }
    }

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
    tc.fVerbose = cfVerbose;
    tc.fSigmaInel = cfSigmaInel;

    tc.fUsePhiWeights = cfUsePhiWeights;
    tc.fUsePtWeights = cfUsePtWeights;
    tc.fUseEtaWeights = cfUseEtaWeights;
    if ((tc.fUsePhiWeights || tc.fUsePtWeights || tc.fUseEtaWeights) && cfFileWithWeights.value.empty()) {
      LOGF(fatal, "weights are switched on, but cfFileWithWeights is empty");
    }

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

    // eta cut copy switch and range
    tc.fUseEtaCut = cfUseEtaCut;
    std::vector<float> lEtaCut = cfEta.value;

    tc.fEta[0] = lEtaCut[0];
    tc.fEta[1] = lEtaCut[1];

    // technical event cuts: copy the switches
    tc.fUseTimeRange = cfUseTimeRange;
    tc.fUseRof = cfUseRof;
    tc.fUsePileUp = cfUsePileUp;
    tc.fUseITSTPC = cfUseITSTPC;
    tc.fUseITSLayers = cfUseITSLayers;
    tc.fUseZFT0vsPV = cfUseZFT0vsPV;
    tc.fUsePrevRof = cfUsePrevRof;

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

    // correlations: <<2>> vs. centrality, in the centrality classes of the reference paper
    cc.fCorrelationsList = new TList();
    cc.fCorrelationsList->SetName("Correlations");
    cc.fCorrelationsList->SetOwner(true);
    fBaseList->Add(cc.fCorrelationsList);
    const std::array<double, 10> centBoundaries = {0., 5., 10., 20., 30., 40., 50., 60., 70., 80.};
    cc.fProfTwo = new TProfile("fProfTwo", "<<2>> = v_{2}{2}^{2} vs. centrality", centBoundaries.size() - 1, centBoundaries.data());
    cc.fProfTwo->GetXaxis()->SetTitle("centrality (%)");
    cc.fProfTwo->GetYaxis()->SetTitle("#LT#LT2#GT#GT");
    cc.fProfTwo->Sumw2();
    cc.fCorrelationsList->Add(cc.fProfTwo);

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
    std::array<const char*, EBeforeAfter_N> baNames = {"before", "after"};
    // book rec histograms only if a process with reconstructed data runs, sim histograms only if one with simulated data runs
    const bool bookRec = doprocessRec || doprocessRecSim;
    const bool bookSim = doprocessRecSim || doprocessSim;
    for (int ba = 0; ba < EBeforeAfter_N; ba++) {

      if (bookRec) {
        // par
        pc.fParticleHistograms[EPt][ERec][ba] = new TH1F(Form("fHistPt[ERec][%s]", baNames[ba]), Form("pt distribution for reconstructed particles (%s cuts)", baNames[ba]), nBins, min, max);
        pc.fParticleHistograms[EPt][ERec][ba]->GetXaxis()->SetTitle("p_{T}");
        pc.fParticleHistogramsList->Add(pc.fParticleHistograms[EPt][ERec][ba]);

        pc.fParticleHistograms[EPhi][ERec][ba] = new TH1F(Form("fHistPhi[ERec][%s]", baNames[ba]), Form("Azimuthal-angle distribution for reconstructed particles (%s cuts)", baNames[ba]), nPhiBins, phiMin, phiMax);
        pc.fParticleHistograms[EPhi][ERec][ba]->GetXaxis()->SetTitle("#varphi (rad)");
        pc.fParticleHistogramsList->Add(pc.fParticleHistograms[EPhi][ERec][ba]);
        // ev
        ec.fEventHistograms[EMultiplicity][ERec][ba] = new TH1F(Form("fHistMultiplicity[ERec][%s]", baNames[ba]), Form("Multiplicity for reconstructed particles (%s cuts)", baNames[ba]), nMultBins, multMin, multMax);
        ec.fEventHistograms[EMultiplicity][ERec][ba]->GetXaxis()->SetTitle("total number of available particles");
        ec.fEventHistogramsList->Add(ec.fEventHistograms[EMultiplicity][ERec][ba]);

        // reference multiplicity (only reconstructed)
        ec.fEventHistograms[EReferenceMultiplicity][ERec][ba] = new TH1F(Form("fHistReferenceMultiplicity[ERec][%s]", baNames[ba]), Form("Reference multiplicity, reconstructed (%s cuts)", baNames[ba]), nRefMultBins, refMultMin, refMultMax);
        ec.fEventHistograms[EReferenceMultiplicity][ERec][ba]->GetXaxis()->SetTitle("reference multiplicity");
        ec.fEventHistogramsList->Add(ec.fEventHistograms[EReferenceMultiplicity][ERec][ba]);

        ec.fEventHistograms[ECentrality][ERec][ba] = new TH1F(Form("fHistCentrality[ERec][%s]", baNames[ba]), Form("Centrality for reconstructed particles (%s cuts)", baNames[ba]), nCentBins, centMin, centMax);
        ec.fEventHistograms[ECentrality][ERec][ba]->GetXaxis()->SetTitle("centrality (%)");
        ec.fEventHistogramsList->Add(ec.fEventHistograms[ECentrality][ERec][ba]);

        ec.fEventHistograms[EVertexX][ERec][ba] = new TH1F(Form("fHistVertexX[ERec][%s]", baNames[ba]), Form("Vertex X for reconstructed particles (%s cuts)", baNames[ba]), nvxBins, vxMin, vxMax);
        ec.fEventHistograms[EVertexX][ERec][ba]->GetXaxis()->SetTitle("vertex X");
        ec.fEventHistogramsList->Add(ec.fEventHistograms[EVertexX][ERec][ba]);

        ec.fEventHistograms[EVertexY][ERec][ba] = new TH1F(Form("fHistVertexY[ERec][%s]", baNames[ba]), Form("Vertex Y for reconstructed particles (%s cuts)", baNames[ba]), nvyBins, vyMin, vyMax);
        ec.fEventHistograms[EVertexY][ERec][ba]->GetXaxis()->SetTitle("vertex Y");
        ec.fEventHistogramsList->Add(ec.fEventHistograms[EVertexY][ERec][ba]);

        ec.fEventHistograms[EVertexZ][ERec][ba] = new TH1F(Form("fHistVertexZ[ERec][%s]", baNames[ba]), Form("Vertex Z for reconstructed particles (%s cuts)", baNames[ba]), nvzBins, vzMin, vzMax);
        ec.fEventHistograms[EVertexZ][ERec][ba]->GetXaxis()->SetTitle("vertex Z");
        ec.fEventHistogramsList->Add(ec.fEventHistograms[EVertexZ][ERec][ba]);
      }
      if (bookSim) {
        pc.fParticleHistograms[EPt][ESim][ba] = new TH1F(Form("fHistPt[ESim][%s]", baNames[ba]), Form("pt distribution for simulated particles (%s cuts)", baNames[ba]), nBins, min, max);
        pc.fParticleHistograms[EPt][ESim][ba]->GetXaxis()->SetTitle("p_{T}");
        pc.fParticleHistogramsList->Add(pc.fParticleHistograms[EPt][ESim][ba]);

        pc.fParticleHistograms[EPhi][ESim][ba] = new TH1F(Form("fHistPhi[ESim][%s]", baNames[ba]), Form("Azimuthal-angle distribution for simulated particles (%s cuts)", baNames[ba]), nPhiBins, phiMin, phiMax);
        pc.fParticleHistograms[EPhi][ESim][ba]->GetXaxis()->SetTitle("#varphi (rad)");
        pc.fParticleHistogramsList->Add(pc.fParticleHistograms[EPhi][ESim][ba]);
        // ev
        ec.fEventHistograms[EMultiplicity][ESim][ba] = new TH1F(Form("fHistMultiplicity[ESim][%s]", baNames[ba]), Form("Multiplicity for simulated particles (%s cuts)", baNames[ba]), nMultBins, multMin, multMax);
        ec.fEventHistograms[EMultiplicity][ESim][ba]->GetXaxis()->SetTitle("total number of available particles");
        ec.fEventHistogramsList->Add(ec.fEventHistograms[EMultiplicity][ESim][ba]);
        // sim from impact parameter
        ec.fEventHistograms[ECentrality][ESim][ba] = new TH1F(Form("fHistCentrality[ESim][%s]", baNames[ba]), Form("Centrality for simulated particles, from impact parameter (%s cuts)", baNames[ba]), nCentBins, centMin, centMax);
        ec.fEventHistograms[ECentrality][ESim][ba]->GetXaxis()->SetTitle("centrality (%)");
        ec.fEventHistogramsList->Add(ec.fEventHistograms[ECentrality][ESim][ba]);
        // impact parameter, only simulated
        ec.fEventHistograms[EImpactParameter][ESim][ba] = new TH1F(Form("fHistImpactParameter[ESim][%s]", baNames[ba]), Form("Impact parameter, simulated (%s cuts)", baNames[ba]), nBBins, bMin, bMax);
        ec.fEventHistograms[EImpactParameter][ESim][ba]->GetXaxis()->SetTitle("b (fm)");
        ec.fEventHistogramsList->Add(ec.fEventHistograms[EImpactParameter][ESim][ba]);

        ec.fEventHistograms[EVertexX][ESim][ba] = new TH1F(Form("fHistVertexX[ESim][%s]", baNames[ba]), Form("Vertex X for simulated particles (%s cuts)", baNames[ba]), nvxBins, vxMin, vxMax);
        ec.fEventHistograms[EVertexX][ESim][ba]->GetXaxis()->SetTitle("vertex X");
        ec.fEventHistogramsList->Add(ec.fEventHistograms[EVertexX][ESim][ba]);

        ec.fEventHistograms[EVertexY][ESim][ba] = new TH1F(Form("fHistVertexY[ESim][%s]", baNames[ba]), Form("Vertex Y for simulated particles (%s cuts)", baNames[ba]), nvyBins, vyMin, vyMax);
        ec.fEventHistograms[EVertexY][ESim][ba]->GetXaxis()->SetTitle("vertex Y");
        ec.fEventHistogramsList->Add(ec.fEventHistograms[EVertexY][ESim][ba]);

        ec.fEventHistograms[EVertexZ][ESim][ba] = new TH1F(Form("fHistVertexZ[ESim][%s]", baNames[ba]), Form("Vertex Z for simulated particles (%s cuts)", baNames[ba]), nvzBins, vzMin, vzMax);
        ec.fEventHistograms[EVertexZ][ESim][ba]->GetXaxis()->SetTitle("vertex Z");
        ec.fEventHistogramsList->Add(ec.fEventHistograms[EVertexZ][ESim][ba]);
      }
    } // end of for (int ba = 0; ba < 2; ba++)

    // QA histograms, booked only if cfFillQA is true
    if (tc.fFillQA) {
      qa.fHistCentralityRecVsSim = new TH2F("fHistCentralityRecVsSim", "Centrality: reconstructed vs. simulated", nCentBins, centMin, centMax, nCentBins, centMin, centMax);
      qa.fHistCentralityRecVsSim->GetXaxis()->SetTitle("centrality simulated (%)");
      qa.fHistCentralityRecVsSim->GetYaxis()->SetTitle("centrality reconstructed (%)");
      qa.fQAHistogramsList->Add(qa.fHistCentralityRecVsSim);

      // reference multiplicity vs. number of contributors, before and after cuts
      std::vector<float> lContribBins = cfNContribBins.value;
      for (int ba = 0; ba < EBeforeAfter_N; ba++) {
        qa.fHistRefMultVsNContrib[ba] = new TH2F(Form("fHistRefMultVsNContrib[%s]", baNames[ba]), Form("reference multiplicity vs. number of contributors (%s cuts)", baNames[ba]), nRefMultBins, refMultMin, refMultMax, static_cast<int>(lContribBins[0]), lContribBins[1], lContribBins[2]);
        qa.fHistRefMultVsNContrib[ba]->GetXaxis()->SetTitle("reference multiplicity");
        qa.fHistRefMultVsNContrib[ba]->GetYaxis()->SetTitle("number of contributors");
        qa.fQAHistogramsList->Add(qa.fHistRefMultVsNContrib[ba]);
      }
    }

  } // end of void init(InitContext&) {

  // A) Process only reconstructed data:
  void processRec(CollisionRec const& collision, aod::BCs const&, TracksRec const& tracks)
  {
    // ...

    // *) steer all analysis steps:
    steer<ERec>(collision, tracks);
  }
  PROCESS_SWITCH(MultiparticleCorrelationsDm, processRec, "process only reconstructed data", true); // yes, keep always one process switch "true", so that there is default running version

  // -------------------------------------------

  // B) Process both reconstructed and corresponding MC truth simulated data:
  void processRecSim(CollisionRecSim const& collision, aod::BCs const&, TracksRecSim const& tracks, aod::McParticles const&, aod::McCollisions const&)
  {
    steer<ERecAndSim>(collision, tracks);
  }
  PROCESS_SWITCH(MultiparticleCorrelationsDm, processRecSim, "process both reconstructed and corresponding MC truth simulated data", false);

  // -------------------------------------------

  // C) Process only simulated data:
  void processSim(CollisionSim const& /*collision*/, aod::BCs const&, TracksSim const& /*tracks*/)
  {
    // steer<ESim>(collision, tracks); // TBI 20241105 not ready yet, but I do not really need this one urgently, since RecSim is working, and I need that one for efficiencies...
  }
  PROCESS_SWITCH(MultiparticleCorrelationsDm, processSim, "process only simulated data", false);

}; // struct MultiparticleCorrelationsDm {

// *) The final touch:
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<MultiparticleCorrelationsDm>(cfgc),
  };
} // WorkflowSpec...
