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

#include <TCollection.h>
#include <TComplex.h>
#include <TF1.h>
#include <TFile.h>
#include <TGrid.h>
#include <TH1.h>
#include <TH2.h>
#include <TList.h>
#include <TMath.h>
#include <TObject.h>
#include <TString.h>
#include <TSystem.h>

#include <RtypesCore.h>

#include <cmath>
#include <cstring>
#include <string>
#include <unordered_map>
#include <vector>

using namespace o2;
using namespace o2::framework;

// Definitions of join tables for Run 3 analysis:
using EventSelection = soa::Join<aod::EvSels, aod::Mults, aod::CentFT0As, aod::CentFT0Cs, aod::CentFT0Ms, aod::CentFV0As, aod::CentNTPVs>;
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
enum eRecSim { eRec = 0,
               eSim,
               eRecAndSim };

enum eProcess {
  eProcessRec = 0, // Run 3, only reconstructed
  eProcessRecSim,  // Run 3, both reconstructed and simulated
  eProcessSim,     // Run 3, only simulated
  eProcess_N
};

enum eEventHistograms {
  eVertexZ = 0,
  ePt,
  eEventHistograms_N
};

enum eCut {
  eBefore = 0,
  eAfter,
  eCut_N
};

// *) Main task:
struct MultiharmonicCorrelations { // this name is used in lower-case format to name the TDirectoryFile in AnalysisResults.root
  Service<ccdb::BasicCCDBManager> ccdb;

  // *) Base TList to hold all output objects:
  TString sBaseListName = "Default list name";
  OutputObj<TList> fBaseList{sBaseListName.Data(),
                             OutputObjHandlingPolicy::AnalysisObject,
                             OutputObjSourceType::OutputObjSource};

  // *) Define configurables:
  Configurable<bool> cfDryRun{"cfDryRun", false, "book all histos and run without filling and calculating anything"}; // example for built-in type (float, string, etc.)
  Configurable<std::vector<float>> cfPtBins{"cfPtBins", {1000, 0., 8.}, "nPtBins, ptMin, ptMax"};                     // example for an array
  Configurable<std::vector<float>> cfPhiBins{"cfPhiBins", {360, 0., o2::constants::math::TwoPI}, "nPhiBins, phiMin, phiMax"};
  Configurable<std::vector<float>> cfCentrBins{"cfCentrBins", {80, 0., 80.}, "nCentrBins, centrMin, centrMax"};
  Configurable<std::vector<float>> cfXBins{"cfXBins", {1000, -0.04, -0.01}, "nXBins, xMin, xMax"};
  Configurable<std::vector<float>> cfYBins{"cfYBins", {1000, -0.01, 0.006}, "nYBins, yMin, yMax"};
  Configurable<std::vector<float>> cfZBins{"cfZBins", {1000, -20., 20.}, "nZBins, zMin, zMax"};
  Configurable<std::vector<float>> cfMultBins{"cfMultBins", {50, 0, 2e4}, "nMultBins, multMin, multMax"};
  Configurable<std::vector<float>> cfMselBins{"cfMselBins", {50, 0, 1e3}, "nMselBins, mselMin, mselMax"};
  Configurable<std::vector<float>> cfTPCnclsBins{"cfTPCnclsBins", {100, 0., 200.}, "ntpcnclsBins, tpnclsMin, tpcnclsMax"};
  Configurable<std::vector<float>> cfDCAxyBins{"cfDCAxyBins", {1000, -0.5, 0.5}, "ndcaxyBins, dcaxyMin, dcaxyMax"};
  Configurable<std::vector<float>> cfDCAzBins{"cfDCAzBins", {1000, -3., 3.}, "ndcazBins, dcazMin, dcazMax"};
  Configurable<std::vector<float>> cfNcontrBins{"cfNcontrBins", {100, 0., 10000.}, "nNContrBins, NContrMin, NContrMax"};

  Configurable<std::string> cfCent{"cfCent", "FT0C", "centrality estimator"};
  Configurable<std::string> cfMult{"cfMult", "TPC", "multiplicity"};
  Configurable<bool> cfQA{"cfQA", true, "quality assurance"};
  Configurable<bool> cfInitsim{"cfInitsim", false, "init histograms of sim"};
  Configurable<bool> cfUseWeights{"cfUseWeights", true, "use weights"};
  Configurable<bool> cfToyModel{"cfToyModel", true, "phi-distribution from toy model"};
  Configurable<bool> cfNest{"cfNest", true, "nested loops"};
  Configurable<bool> cfTechcuts{"cfTechcuts", true, "technical cuts"};

  Configurable<std::vector<float>> cfVertexZ{"cfVertexZ", {-10, 10.}, "vertex z position range: {min, max}[cm], with convention: min <= Vz < max"};
  Configurable<std::vector<float>> cfPt{"cfPt", {0.2, 5.0}, "transverse momentum range"};
  Configurable<std::vector<float>> cfEta{"cfEta", {-0.8, 0.8}, "eta range"};
  Configurable<std::vector<float>> cfDCAxy{"cfDCAxy", {-0.5, 0.5}, "dca xy range"};
  Configurable<std::vector<float>> cfDCAz{"cfDCAz", {-0.2, 0.2}, "dca z range"};
  Configurable<float> cftpcNClsFoundmin{"cftpcNClsFoundmin", 70., "tpcNClsFoundmin"};
  Configurable<float> cftpcChi2NClmax{"cftpcChi2NClmax", 4., "tpcChi2NClmax"};

  Configurable<std::vector<int>> cfRuns{"cfRuns", {544091, 544095, 544098, 544116, 544121, 544122, 544123, 544124}, "List of run numbers to analyze"};
  Configurable<std::string> cfFileWithWeights{"cfFileWithWeights", "/alice-ccdb.cern.ch/Users/p/pengchon/weightsfile06", "path to external ROOT file which holds all particle weights"};

  // *) Define and initialize all data members to be called in the main process* functions:
  // **) Task configuration:
  struct TaskConfiguration {
    bool fProcess[eProcess_N] = {false}; // Set what to process. See enum eProcess for full description. Set via implicit variables within a PROCESS_SWITCH clause.
    bool fDryRun = false;                // book all histos and run without filling and calculating anything
  } tc;                                  // you have to prepend "tc." for all objects name in this group later in the code

  // **) Particle histograms:
  struct ParticleHistograms {
    TList* fParticleHistogramsList = NULL; //!<! list to hold all control particle histograms
    TH1F* fHistPt[2] = {NULL};             // pt distribution of a particle [ 0 = rec, 1 = sim ]
    TH1F* fHistPhi[2] = {NULL};
    TH1F* fHistCharge[2] = {NULL};  // charge distribution
    TH1F* fHistTPCncls[2] = {NULL}; // TPCNClsFindable
    TH1F* fHistTracksdcaXY[2] = {NULL};
    TH1F* fHistTracksdcaZ[2] = {NULL};
    TH1F* histWeights = NULL;
  } pc; // you have to prepend "pc." for all objects name in this group later in the code

  struct EventHistograms {
    TList* fEventHistogramsList = NULL;
    TH1F* fHistCentr[2] = {NULL};
    TH1I* fHistMult[2] = {NULL};
    TH1F* fHistMsel[2] = {NULL};
    TH1F* fHistX[2] = {NULL};
    TH1F* fHistY[2] = {NULL};
    TH1F* fHistZ[2] = {NULL};
    TH1I* fHistNContr[2] = {NULL};
    TH1F* fEventHistograms[eEventHistograms_N][2][2] = {{{NULL}}}; //! [ type - see enum eEventHistograms ][reco,sim][before, after event cuts]
  } event;

  struct QA {
    TList* fQAList = NULL;
    TH2F* fQA = NULL;
    TH2F* fQAM_NC = NULL;
  } qa;

  static constexpr int maxHarmonic = 7;
  static constexpr int maxPower = 5;
  struct CorrelationVariables {
    TList* fCorrelationVariablesList = NULL;
    TProfile* pv22_centr = NULL;
    TProfile* pv32_centr = NULL;
    TProfile* pv42_centr = NULL;
    TProfile* pfour32_centr = NULL;
    TProfile* pfour42_centr = NULL;
    TComplex Qvector[maxHarmonic][maxPower];
    std::vector<float> vecphi;
    std::vector<float> vecwei;
    TProfile* nestedLoops[maxHarmonic] = {NULL};
    TProfile* pv2_nest = NULL;
    TProfile* pv3_nest = NULL;
    TProfile* pv4_nest = NULL;
  } cor;

  struct PhiHist {
    TList* fPhiHistList = NULL;
    std::unordered_map<int, TH1F*> histMap;
  } phih;

  struct WeightsHist {
    TList* fWeightsHistList = NULL;
    std::unordered_map<int, TH1F*> weightsmap;
  } wh;

  TObject* GetObjectFromList(TList* list, const char* objectName)
  {
    // Get TObject pointer from TList, even if it's in some nested TList. Foreseen
    // to be used to fetch histograms or profiles from files directly.
    // Some ideas taken from TCollection::ls()
    // If you have added histograms directly to files (without TList's), then you can fetch them directly with
    // file->Get("hist-name").

    // Usage: TH1D *hist = (TH1D*) GetObjectFromList("some-valid-TList-pointer","some-object-name");

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
        objectFinal = GetObjectFromList(reinterpret_cast<TList*>(objectIter), objectName);
        if (objectFinal)
          return objectFinal;
      }
    } // while(objectIter = next())

    return NULL;

  } // TObject* GetObjectFromList(TList *list, char *objectName)

  TH1F* GetHistogramWithWeights(const char* filePath, const char* runNumber)
  {
    // *) Return value:
    TH1F* hist = NULL;
    TList* baseList = NULL;     // base top-level list in the TFile, e.g. named "ccdb_object"
    TList* listWithRuns = NULL; // nested list with run-wise TList's holding run-specific weights

    // *) Determine from filePath if the file is on a local machine, or in home dir AliEn, or in CCDB:
    //    Algorithm: If filePath begins with "/alice/cern.ch/" then it's in the home dir AliEn;
    //                     If filePath begins with "/alice-ccdb.cern.ch/" then it's in CCDB. Therefore, files in AliEn and CCDB must be specified with abs path;
    //                     for local files both abs and relative paths are just fine.
    bool bFileIsInAliEn = false;
    bool bFileIsInCCDB = false;

    string pathstr = filePath;
    const string pathalien = "/alice/cern.ch/";
    const string pathccdb = "/alice-ccdb.cern.ch/";
    if (pathstr.starts_with(pathalien)) {
      bFileIsInAliEn = true;
    } else if (pathstr.starts_with(pathccdb)) {
      bFileIsInCCDB = true;
    }
    LOGF(info, "bFileIsInCCDB= %d", bFileIsInCCDB);

    if (bFileIsInAliEn) {
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
      // for thare there is a helper utility function GetObjectFromList(...) , see its implementation further below
      listWithRuns = reinterpret_cast<TList*>(GetObjectFromList(baseList, runNumber));
      if (!listWithRuns) {
        TString runNumberWithLeadingZeroes = "000";
        runNumberWithLeadingZeroes += runNumber; // another try, with "000" prepended to run number
        listWithRuns = reinterpret_cast<TList*>(GetObjectFromList(baseList, runNumberWithLeadingZeroes.Data()));
        if (!listWithRuns) {
          LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
        }
      }

    } else if (bFileIsInCCDB) {
      // File you want to access is in your home dir in CCDB:
      // Remember that here I do not access the file; instead, I directly access the object in that file.
      ccdb->setURL("http://alice-ccdb.cern.ch"); // to be able to use "ccdb" this object in your analysis task, see 4b/ below
      baseList = reinterpret_cast<TList*>(ccdb->get<TList>(TString(filePath).ReplaceAll("/alice-ccdb.cern.ch/", "").Data()));
      if (!baseList) {
        LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
      }

      listWithRuns = reinterpret_cast<TList*>(GetObjectFromList(baseList, runNumber));
      if (!listWithRuns) {
        TString runNumberWithLeadingZeroes = "000";
        runNumberWithLeadingZeroes += runNumber; // another try, with "000" prepended to run number
        listWithRuns = reinterpret_cast<TList*>(GetObjectFromList(baseList, runNumberWithLeadingZeroes.Data()));
        if (!listWithRuns) {
          LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
        }
      }

      // OK, we got the desired TList with efficiency corrections, after that we can use the common code for all 3 cases (local, AliEn, CCDB, that common code is below)

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

      listWithRuns = reinterpret_cast<TList*>(GetObjectFromList(baseList, runNumber));
      if (!listWithRuns) {
        TString runNumberWithLeadingZeroes = "000";
        runNumberWithLeadingZeroes += runNumber; // another try, with "000" prepended to run number
        listWithRuns = reinterpret_cast<TList*>(GetObjectFromList(baseList, runNumberWithLeadingZeroes.Data()));
        if (!listWithRuns) {
          // baseList->ls();
          LOGF(fatal, "\033[1;31m%s at line %d : this crash can happen if in the output file there is no list with weights for the current run number = %s\033[0m", __FUNCTION__, __LINE__ /*, tc.fRunNumber.Data()*/);
        }
      }

    } // end of  else

    // Here comes the common code for all three cases, where from "listWithRuns" you fetch the desired histogram with efficiency corrections:
    listWithRuns->ls();

    hist = reinterpret_cast<TH1F*>(listWithRuns->FindObject("histWithEfficiencyCorrections"));
    if (!hist) {
      LOGF(fatal, "no histweight");
    }

    // Once you have a valid pointer to "hist", add these technical lines:
    hist->SetDirectory(0);                                    // remove ownerhip from an external file
    TH1F* histClone = reinterpret_cast<TH1F*>(hist->Clone()); // yes, I have to clone here

    return histClone;

  } // end of TH1F* GetHistogramWithWeights(const char* filePath, const char* runNumber) {

  // ...

  // *) Define all member functions to be called in the main process* functions:
  template <eRecSim rs, typename T1>
  bool EventCuts(T1 const& collision)
  {
    vector<float> vertexZ = cfVertexZ.value;
    float vertexZmin = static_cast<float>(vertexZ[0]);
    float vertexZmax = static_cast<float>(vertexZ[1]);
    float posZ = collision.posZ();
    float NContrcut = 2.;
    if (posZ < vertexZmin || posZ > vertexZmax || (!collision.sel8()) || collision.numContrib() < NContrcut)
      return false;
    if (cfTechcuts) {
      if (!collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard) || !collision.selection_bit(o2::aod::evsel::kNoCollInRofStandard) || !collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup) || !collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC) || !collision.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll) || !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV) || !collision.selection_bit(o2::aod::evsel::kNoHighMultCollInPrevRof))
        return false;
    }
    return true;
  }

  template <eRecSim rs, typename T>
  bool ParticleCuts(T const& track)
  {
    vector<float> Pt = cfPt.value;
    float ptcutmin = static_cast<float>(Pt[0]);
    float ptcutmax = static_cast<float>(Pt[1]);
    float pt = track.pt();
    vector<float> Eta = cfEta.value;
    float etacutmin = static_cast<float>(Eta[0]);
    float etacutmax = static_cast<float>(Eta[1]);
    float eta = track.eta();
    vector<float> dcaXY = cfDCAxy.value;
    float dcaxycutmin = static_cast<float>(dcaXY[0]);
    float dcaxycutmax = static_cast<float>(dcaXY[1]);
    float dcaxy = track.dcaXY();
    vector<float> dcaZ = cfDCAz.value;
    float dcazcutmin = static_cast<float>(dcaZ[0]);
    float dcazcutmax = static_cast<float>(dcaZ[1]);
    float dcaz = track.dcaZ();
    if (pt < ptcutmin || pt > ptcutmax || eta < etacutmin || eta > etacutmax || dcaxy < dcaxycutmin || dcaxy > dcaxycutmax || dcaz < dcazcutmin || dcaz > dcazcutmax || (!track.isPrimaryTrack()) || (!track.isPVContributor()) || track.tpcNClsFound() < cftpcNClsFoundmin.value || track.tpcChi2NCl() > cftpcChi2NClmax)
      return false;
    return true;
  }

  TComplex Q(int n, int p)
  {
    // Using the fact that Q{-n,p} = Q{n,p}^*.
    if (n >= 0) {
      return cor.Qvector[n][p];
    }
    return TComplex::Conjugate(cor.Qvector[-n][p]);
  }

  TComplex Two(Int_t n1, Int_t n2)
  {
    // Generic two-particle correlation <exp[i(n1*phi1+n2*phi2)]>.
    TComplex two = Q(n1, 1) * Q(n2, 1) - Q(n1 + n2, 2);
    return two;

  } // TComplex Two(Int_t n1, Int_t n2)
  TComplex Four(Int_t n1, Int_t n2, Int_t n3, Int_t n4)
  { // Generic four-particle correlation <exp[i(n1*phi1+n2*phi2+n3*phi3+n4*phi4)]>.
    TComplex four = Q(n1, 1) * Q(n2, 1) * Q(n3, 1) * Q(n4, 1) - Q(n1 + n2, 2) * Q(n3, 1) * Q(n4, 1) - Q(n2, 1) * Q(n1 + n3, 2) * Q(n4, 1) - Q(n1, 1) * Q(n2 + n3, 2) * Q(n4, 1) + 2. * Q(n1 + n2 + n3, 3) * Q(n4, 1) - Q(n2, 1) * Q(n3, 1) * Q(n1 + n4, 2) + Q(n2 + n3, 2) * Q(n1 + n4, 2) - Q(n1, 1) * Q(n3, 1) * Q(n2 + n4, 2) + Q(n1 + n3, 2) * Q(n2 + n4, 2) + 2. * Q(n3, 1) * Q(n1 + n2 + n4, 3) - Q(n1, 1) * Q(n2, 1) * Q(n3 + n4, 2) + Q(n1 + n2, 2) * Q(n3 + n4, 2) + 2. * Q(n2, 1) * Q(n1 + n3 + n4, 3) + 2. * Q(n1, 1) * Q(n2 + n3 + n4, 3) - 6. * Q(n1 + n2 + n3 + n4, 4);
    return four;
  } // TComplex Four(Int_t n1, Int_t n2, Int_t n3, Int_t n4)

  static double pdf(const double* x, const double* par)
  {
    double y = 1;
    int harm = 6;
    for (int i = 0; i < harm; i = i + 1) {
      y = y + 2 * (0.04 + (i + 1.) * 0.01) * TMath::Cos((i + 1) * (x[0] - par[0]));
    }
    return y;
  }

  template <eRecSim rs, typename T1, typename T2>
  void Steer(T1 const& collision, T2 const& tracks)
  {
    // Dry run:
    if (tc.fDryRun) {
      return;
    }
    // Print current run number:
    int currentRun = collision.bc().runNumber();
    auto it = phih.histMap.find(currentRun);
    auto histweight = wh.weightsmap.find(currentRun);
    float centr = 0, M = 0., msel = 0.;
    // TF1* f = new TF1("f", pdf, 0, TMath::TwoPi(), 1);
    TF1* f = new TF1("f",
                     "1 +"
                     "2 * (0.05) * cos(1 * (x - [0])) +"
                     "2 * (0.06) * cos(2 * (x - [0])) +"
                     "2 * (0.07) * cos(3 * (x - [0])) +"
                     "2 * (0.08) * cos(4 * (x - [0])) +"
                     "2 * (0.09) * cos(5 * (x - [0])) +"
                     "2 * (0.10) * cos(6 * (x - [0]))",
                     0, TMath::TwoPi());
    f->SetParameters(0.);

    if constexpr (rs == eRec || rs == eRecAndSim) {
      if (cfCent.value == "FT0C")
        centr = collision.centFT0C();
      else if (cfCent.value == "FT0M")
        centr = collision.centFT0M();
      else if (cfCent.value == "FT0A")
        centr = collision.centFT0A();

      if (cfMult.value == "TPC")
        M = collision.multTPC();
      else if (cfMult.value == "FV0M")
        M = collision.multFV0M();
      else if (cfMult.value == "FT0C")
        M = collision.multFT0C();
      else if (cfMult.value == "FT0M")
        M = collision.multFT0M();
      else if (cfMult.value == "NTracksPV")
        M = collision.multNTracksPV();

      event.fHistX[eBefore]->Fill(collision.posX());
      event.fHistY[eBefore]->Fill(collision.posY());
      event.fHistZ[eBefore]->Fill(collision.posZ());
      event.fHistCentr[eBefore]->Fill(centr);
      event.fHistMult[eBefore]->Fill(M);
      event.fHistNContr[eBefore]->Fill(collision.numContrib());

      event.fEventHistograms[eVertexZ][eRec][0]->Fill(collision.posZ());

      // *) Event cuts:
      float centrcut = 80.;
      if (!EventCuts<rs>(collision) || centr > centrcut) { // Main call for event cuts
        return;
      }

      event.fHistX[eAfter]->Fill(collision.posX());
      event.fHistY[eAfter]->Fill(collision.posY());
      event.fHistZ[eAfter]->Fill(collision.posZ());
      event.fHistCentr[eAfter]->Fill(centr);
      event.fHistMult[eAfter]->Fill(M);
      event.fHistNContr[eAfter]->Fill(collision.numContrib());

      event.fEventHistograms[eVertexZ][eRec][1]->Fill(collision.posZ());
      qa.fQAM_NC->Fill(M, collision.numContrib());

      if constexpr (rs == eRecAndSim) {
        auto mccollision = collision.mcCollision();
        float b = mccollision.impactParameter();
        float centrsim = o2::constants::math::PI * b * b / 7.71;
        event.fHistX[eSim]->Fill(mccollision.posX());
        event.fHistY[eSim]->Fill(mccollision.posY());
        event.fHistZ[eSim]->Fill(mccollision.posZ());
        event.fEventHistograms[eVertexZ][eSim][1]->Fill(mccollision.posZ());
        event.fHistCentr[eSim]->Fill(centrsim);
        qa.fQA->Fill(centrsim, centr);
        centr = centrsim;
      }
    }

    // before loop over particles
    float phi = 0, weight = 1.;
    vector<float>().swap(cor.vecphi);
    vector<float>().swap(cor.vecwei);
    for (int ih = 0; ih < maxHarmonic; ih++) {
      for (int ip = 0; ip < maxPower; ip++) {
        cor.Qvector[ih][ip] = TComplex(0., 0.);
      }
    }

    /*
    for (int ih = 0; ih < maxHarmonic; ih++) {
      for (int ip = 0; ip < maxPower; ip++) {
        LOGF(info, "Qvector[%d][%d]=%f, ", ih, ip, cor.Qvector[ih][ip].Rho2());
      }
      LOGF(info, "\n");
    }
    LOGF(info, "\n\n");
    */

    // Main loop over particles:
    for (const auto& track : tracks) {

      // Fill reconstructed ...:
      float ptrec = 0., ptsim = 0.;
      if constexpr (rs == eRec || rs == eRecAndSim) {
        // Fill track pt distribution:

        pc.fHistPt[eBefore]->Fill(track.pt());
        pc.fHistPhi[eBefore]->Fill(track.phi());
        pc.fHistCharge[eBefore]->Fill(track.sign());
        pc.fHistTPCncls[eBefore]->Fill(track.tpcNClsFindable());
        pc.fHistTracksdcaXY[eBefore]->Fill(track.dcaXY());
        pc.fHistTracksdcaZ[eBefore]->Fill(track.dcaZ());

        event.fEventHistograms[ePt][eRec][0]->Fill(track.pt());
        ptrec = track.pt();

        // *) Particle cuts:
        if (!ParticleCuts<rs>(track)) { // Main call for particle cuts.
          continue;                     // not return!!
        }
        pc.fHistPt[eAfter]->Fill(track.pt());
        pc.fHistPhi[eAfter]->Fill(track.phi());
        pc.fHistCharge[eAfter]->Fill(track.sign());
        pc.fHistTPCncls[eAfter]->Fill(track.tpcNClsFindable());
        pc.fHistTracksdcaXY[eAfter]->Fill(track.dcaXY());
        pc.fHistTracksdcaZ[eAfter]->Fill(track.dcaZ());

        event.fEventHistograms[ePt][eRec][1]->Fill(ptrec);

        phi = track.phi();
        if (cfToyModel) {
          phi = f->GetRandom();
        }
        cor.vecphi.push_back(phi);
        if (it != phih.histMap.end()) {
          it->second->Fill(phi);
        }

        if (cfUseWeights) {
          if (histweight != wh.weightsmap.end() && histweight->second) {
            weight = histweight->second->GetBinContent(histweight->second->FindBin(phi));
          } else {
            LOG(warning) << "No weights found for run " << currentRun << ", using weight=1";
            weight = 1;
          }
        } else {
          weight = 1;
        }
        cor.vecwei.push_back(weight);

        // ... and corresponding MC truth simulated:
        // See https://github.com/AliceO2Group/O2Physics/blob/master/Tutorials/src/mcHistograms.cxx
        // See https://aliceo2group.github.io/analysis-framework/docs/datamodel/ao2dTables.html#montecarlo
        if constexpr (rs == eRecAndSim) {
          if (!track.has_mcParticle()) {
            LOGF(warning, "  No MC particle for this track, skip...");
            return;
          }
          auto mcparticle = track.mcParticle(); // corresponding MC truth simulated particle
          pc.fHistPt[eSim]->Fill(mcparticle.pt());
          event.fEventHistograms[ePt][eSim][0]->Fill(mcparticle.pt());
          ptsim = mcparticle.pt();
          event.fEventHistograms[ePt][eSim][0]->Fill(ptsim);
          phi = mcparticle.phi();
          pc.fHistPhi[eSim]->Fill(mcparticle.phi());
        } // end of if constexpr (rs == eRecAndSim)

      } // if constexpr (rs == eRec || rs == eRecAndSim)

      // analysis in the loop over particle
      msel = msel + 1;
      for (int ih = 0; ih < maxHarmonic; ih++) {
        for (int ip = 0; ip < maxPower; ip++) {
          cor.Qvector[ih][ip] += TComplex(std::pow(weight, ip) * TMath::Cos(ih * phi), std::pow(weight, ip) * TMath::Sin(ih * phi));
        }
      }
    } // end of for (auto track: tracks)
    event.fHistMsel[eAfter]->Fill(msel);

    if (cfNest) {
      float phi1 = 0., phi2 = 0., weight1 = 1., weight2 = 2.;
      for (Int_t c = 0; c < maxHarmonic; c++) {
        delete cor.nestedLoops[c];
        cor.nestedLoops[c] = new TProfile("", "", 1, 0., 1.);
        cor.nestedLoops[c]->Sumw2();
      }
      for (int i1 = 0; i1 < static_cast<int>(cor.vecphi.size()); i1++) { // nested loop of particles
        phi1 = cor.vecphi[i1];
        weight1 = cor.vecwei[i1];
        for (int i2 = 0; i2 < static_cast<int>(cor.vecphi.size()); i2++) {
          if (i2 == i1) {
            continue;
          }
          phi2 = cor.vecphi[i2];
          weight2 = cor.vecwei[i2];
          cor.nestedLoops[0]->Fill(0.5, TMath::Cos(2 * phi1 - 2 * phi2), weight1 * weight2);
          cor.nestedLoops[1]->Fill(0.5, TMath::Cos(3 * phi1 - 3 * phi2), weight1 * weight2);
          cor.nestedLoops[2]->Fill(0.5, TMath::Cos(4 * phi1 - 4 * phi2), weight1 * weight2);
        }
      } // end of two nested loop
    }

    // calculate correlations
    float Mmin = 4.;
    if (msel < Mmin)
      return;
    float wTwo = Two(0, 0).Re();
    float wFour = Four(0, 0, 0, 0).Re();
    float four32 = Four(3, 2, -3, -2).Re() / wFour;
    float four42 = Four(4, 2, -4, -2).Re() / wFour;
    float v22 = Two(2, -2).Re() / wTwo;
    float v32 = Two(3, -3).Re() / wTwo;
    float v42 = Two(4, -4).Re() / wTwo;
    if (std::isnan(v22) || std::isnan(v32) || std::isnan(v42) || std::isnan(four32) || std::isnan(four42)) {
      LOGF(info, "\033[1;31m%s std::isnan(v22) || std::isnan(v32) || std::isnan(v42) || std::isnan(four32) || std::isnan(four42)\033[0m", __FUNCTION__);
      LOGF(error, "v22 = %f\nv32 = %f\nv42 = %f\nfour32=%f\nv42 = %f\n", v22, v32, v42, four32, four42);
      return;
    }

    cor.pv22_centr->Fill(centr, v22, wTwo);
    cor.pv32_centr->Fill(centr, v32, wTwo);
    cor.pv42_centr->Fill(centr, v42, wTwo);
    cor.pfour32_centr->Fill(centr, four32, wFour);
    cor.pfour42_centr->Fill(centr, four42, wFour);

    if (cfNest) {
      cor.pv2_nest->Fill(centr, cor.nestedLoops[0]->GetBinContent(1), wTwo);
      cor.pv3_nest->Fill(centr, cor.nestedLoops[1]->GetBinContent(1), wTwo);
      cor.pv4_nest->Fill(centr, cor.nestedLoops[2]->GetBinContent(1), wTwo);

      LOGF(info, "v22=%f, v22_nest=%f", v22, cor.nestedLoops[0]->GetBinContent(1));
    }
  } // end of template <eRecSim rs, typename T1, typename T2> void Steer(T1 const& collision, T2 const& tracks)

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

    // *) Book base list:
    TList* temp = new TList();
    temp->SetOwner(true);
    fBaseList.setObject(temp);

    // *) Book and nest all other TLists:
    pc.fParticleHistogramsList = new TList();
    pc.fParticleHistogramsList->SetName("ParticleHistograms");
    pc.fParticleHistogramsList->SetOwner(true);
    fBaseList->Add(pc.fParticleHistogramsList); // any nested TList in the base TList appears as a subdir in the output ROOT file
    event.fEventHistogramsList = new TList();
    event.fEventHistogramsList->SetName("EventHistograms");
    event.fEventHistogramsList->SetOwner(true);
    fBaseList->Add(event.fEventHistogramsList);

    qa.fQAList = new TList();
    qa.fQAList->SetName("QAHistograms");
    qa.fQAList->SetOwner(true);
    fBaseList->Add(qa.fQAList);

    cor.fCorrelationVariablesList = new TList();
    cor.fCorrelationVariablesList->SetName("CorrelationVariables");
    cor.fCorrelationVariablesList->SetOwner(true);
    fBaseList->Add(cor.fCorrelationVariablesList);

    phih.fPhiHistList = new TList();
    phih.fPhiHistList->SetName("PhiHistograms");
    phih.fPhiHistList->SetOwner(true);
    fBaseList->Add(phih.fPhiHistList);

    wh.fWeightsHistList = new TList();
    wh.fWeightsHistList->SetName("WeightsHistograms");
    wh.fWeightsHistList->SetOwner(true);
    fBaseList->Add(wh.fWeightsHistList);

    // *) Book pt distribution with binning defined through configurables in the json file:
    vector<float> l_pt_bins = cfPtBins.value; // define local array and initialize it from an array set in the configurables
    vector<float> l_phi_bins = cfPhiBins.value;
    vector<float> l_centr_bins = cfCentrBins.value;
    vector<float> l_x_bins = cfXBins.value;
    vector<float> l_y_bins = cfYBins.value;
    vector<float> l_z_bins = cfZBins.value;
    vector<float> l_mult_bins = cfMultBins.value;
    vector<float> l_msel_bins = cfMselBins.value;
    vector<float> l_tpcncls_bins = cfTPCnclsBins.value;
    vector<float> l_dcaxy_bins = cfDCAxyBins.value;
    vector<float> l_dcaz_bins = cfDCAzBins.value;
    vector<float> l_ncontr_bins = cfNcontrBins.value;
    vector<int> targetRuns = cfRuns.value;
    int nBins = static_cast<int>(l_pt_bins[0]);
    int nBinsphi = static_cast<int>(l_phi_bins[0]);
    int nBinscentr = static_cast<int>(l_centr_bins[0]);
    int nBinsx = static_cast<int>(l_x_bins[0]);
    int nBinsy = static_cast<int>(l_y_bins[0]);
    int nBinsz = static_cast<int>(l_z_bins[0]);
    int nBinsmult = static_cast<int>(l_mult_bins[0]);
    int nBinsmsel = static_cast<int>(l_msel_bins[0]);
    int nBinscharge = 2;
    int nBinstpcncls = static_cast<int>(l_tpcncls_bins[0]);
    int nBinsdcaxy = static_cast<int>(l_dcaxy_bins[0]);
    int nBinsdcaz = static_cast<int>(l_dcaz_bins[0]);
    int nBinsncontr = static_cast<int>(l_ncontr_bins[0]);

    float min = l_pt_bins[1];
    float max = l_pt_bins[2];
    float minphi = l_phi_bins[1];
    float maxphi = l_phi_bins[2];
    float mincentr = l_centr_bins[1];
    float maxcentr = l_centr_bins[2];
    float minx = l_x_bins[1];
    float maxx = l_x_bins[2];
    float miny = l_y_bins[1];
    float maxy = l_y_bins[2];
    float minz = l_z_bins[1];
    float maxz = l_z_bins[2];
    float minmult = l_mult_bins[1];
    float maxmult = l_mult_bins[2];
    float minmsel = l_msel_bins[1];
    float maxmsel = l_msel_bins[2];
    float mincharge = -2.;
    float maxcharge = 2.;
    float mintpcncls = l_tpcncls_bins[1];
    float maxtpcncls = l_tpcncls_bins[2];
    float mindcaxy = l_dcaxy_bins[1];
    float maxdcaxy = l_dcaxy_bins[2];
    float mindcaz = l_dcaz_bins[1];
    float maxdcaz = l_dcaz_bins[2];
    float minncontr = l_ncontr_bins[1];
    float maxncontr = l_ncontr_bins[2];

    const char* cevent[] = {"vertexZ", "Pt"};
    const char* cpro[] = {"rec", "sim"};
    const char* ccut[] = {"before", "after"};
    for (int i = 0; i < eEventHistograms_N; i++) {
      for (int j = 0; j < eRecAndSim; j++) {
        for (int k = 0; k < eCut_N; k++) {
          TString histname = Form("fEventHistograms[%s][%s][%s]", cevent[i], cpro[j], ccut[k]);
          TString histtitle = Form("%s distribution for %s, %s cut", cevent[i], cpro[j], ccut[k]);
          if (i == 0)
            event.fEventHistograms[i][j][k] = new TH1F(histname, histtitle, nBinsz, minz, maxz);
          if (i == 1)
            event.fEventHistograms[i][j][k] = new TH1F(histname, histtitle, nBins, min, max);
          event.fEventHistograms[i][j][k]->GetXaxis()->SetTitle(Form("%s", cevent[i]));
          // event.fEventHistogramsList->Add(event.fEventHistograms[i][j][k]);
        }
      }
    }

    for (int icut = 0; icut < eCut_N; icut++) {
      pc.fHistPt[icut] = new TH1F(Form("fHistPt[%s]", ccut[icut]), Form("pt distribution %s cut for reconstructed particles", ccut[icut]), nBins, min, max);
      pc.fHistPhi[icut] = new TH1F(Form("fHistPhi[%s]", ccut[icut]), Form("phi distribution %s cut for reconstructed particles", ccut[icut]), nBinsphi, minphi, maxphi);
      pc.fHistCharge[icut] = new TH1F(Form("fHistCharge[%s]", ccut[icut]), Form("charge distribution %s cut for reconstructed particles", ccut[icut]), nBinscharge, mincharge, maxcharge);
      pc.fHistTPCncls[icut] = new TH1F(Form("fHistTPCncls[%s]", ccut[icut]), Form("tpcncls distribution %s cut for reconstructed particles", ccut[icut]), nBinstpcncls, mintpcncls, maxtpcncls);
      pc.fHistTracksdcaXY[icut] = new TH1F(Form("fHistTracksdcaXY[%s]", ccut[icut]), Form("dcaxy distribution %s cut for reconstructed particles", ccut[icut]), nBinsdcaxy, mindcaxy, maxdcaxy);
      pc.fHistTracksdcaZ[icut] = new TH1F(Form("fHistTracksdcaZ[%s]", ccut[icut]), Form("dcaz distribution %s cut for reconstructed particles", ccut[icut]), nBinsdcaz, mindcaz, maxdcaz);
      pc.fHistPt[icut]->GetXaxis()->SetTitle("p_{T}");
      pc.fHistPhi[icut]->GetXaxis()->SetTitle("phi");
      pc.fHistCharge[icut]->GetXaxis()->SetTitle("charge");
      pc.fHistTPCncls[icut]->GetXaxis()->SetTitle("TPCNClsFindable");
      pc.fHistTracksdcaXY[icut]->GetXaxis()->SetTitle("DCA XY");
      pc.fHistTracksdcaZ[icut]->GetXaxis()->SetTitle("DCA Z");
      pc.fParticleHistogramsList->Add(pc.fHistPt[icut]);
      pc.fParticleHistogramsList->Add(pc.fHistPhi[icut]);
      pc.fParticleHistogramsList->Add(pc.fHistCharge[icut]);
      pc.fParticleHistogramsList->Add(pc.fHistTPCncls[icut]);
      pc.fParticleHistogramsList->Add(pc.fHistTracksdcaXY[icut]);
      pc.fParticleHistogramsList->Add(pc.fHistTracksdcaZ[icut]);

      // init eventhist
      event.fHistCentr[icut] = new TH1F(Form("fHistCentr[%s]", ccut[icut]), Form("centrality distribution %s cut for reconstructed particles", ccut[icut]), nBinscentr, mincentr, maxcentr);
      event.fHistX[icut] = new TH1F(Form("fHistX[%s]", ccut[icut]), Form("posX distribution %s cut for reconstructed particles", ccut[icut]), nBinsx, minx, maxx);
      event.fHistY[icut] = new TH1F(Form("fHistY[%s]", ccut[icut]), Form("posY distribution %s cut for reconstructed particles", ccut[icut]), nBinsy, miny, maxy);
      event.fHistZ[icut] = new TH1F(Form("fHistZ[%s]", ccut[icut]), Form("posZ distribution %s cut for reconstructed particles", ccut[icut]), nBinsz, minz, maxz);
      event.fHistMult[icut] = new TH1I(Form("fHistMult[%s]", ccut[icut]), Form("mult distribution %s cut for reconstructed particles", ccut[icut]), nBinsmult, minmult, maxmult);
      event.fHistMsel[icut] = new TH1F(Form("fHistMsel[%s]", ccut[icut]), Form("selected tracks %s cut", ccut[icut]), nBinsmsel, minmsel, maxmsel);
      event.fHistNContr[icut] = new TH1I(Form("fHistNContr[%s]", ccut[icut]), Form("NContr distribution %s cut", ccut[icut]), nBinsncontr, minncontr, maxncontr);
      event.fHistCentr[icut]->GetXaxis()->SetTitle(Form("centrality, %s", cfCent.value.c_str()));
      event.fHistX[icut]->GetXaxis()->SetTitle("x");
      event.fHistY[icut]->GetXaxis()->SetTitle("y");
      event.fHistZ[icut]->GetXaxis()->SetTitle("z");
      event.fHistMult[icut]->GetXaxis()->SetTitle(Form("multiplicity, %s", cfMult.value.c_str()));
      event.fHistMsel[icut]->GetXaxis()->SetTitle("selected tracks");
      event.fHistNContr[icut]->GetXaxis()->SetTitle("numContrib");
      event.fEventHistogramsList->Add(event.fHistCentr[icut]);
      event.fEventHistogramsList->Add(event.fHistX[icut]);
      event.fEventHistogramsList->Add(event.fHistY[icut]);
      event.fEventHistogramsList->Add(event.fHistZ[icut]);
      event.fEventHistogramsList->Add(event.fHistMult[icut]);
      event.fEventHistogramsList->Add(event.fHistMsel[icut]);
      event.fEventHistogramsList->Add(event.fHistNContr[icut]);
    }

    // init of sim histograms
    if (cfInitsim) {
      pc.fHistPt[eSim] = new TH1F("fHistPt[eSim]", "pt distribution for simulated particles", nBins, min, max);
      pc.fHistPhi[eSim] = new TH1F("fHistPhi[eSim]", "phi distribution for simulated particles", nBinsphi, minphi, maxphi);
      pc.fHistCharge[eSim] = new TH1F("fHistCharge[eSim]", "charge distribution for simulated particles", nBinscharge, mincharge, maxcharge);
      pc.fHistTPCncls[eSim] = new TH1F("fHistTPCncls[eSim]", "tpcncls distribution for simulated particles", nBinstpcncls, minphi, maxtpcncls);
      pc.fHistTracksdcaXY[eSim] = new TH1F("fHistTracksdcaXY[eSim]", "dcaxy distribution for simulated particles", nBinsdcaxy, mindcaxy, maxdcaxy);
      pc.fHistTracksdcaZ[eSim] = new TH1F("fHistTracksdcaZ[eSim]", "dcaz distribution for simulated particles", nBinsdcaz, mindcaz, maxdcaz);
      pc.fHistPt[eSim]->GetXaxis()->SetTitle("p_{T}");
      pc.fHistPhi[eSim]->GetXaxis()->SetTitle("phi");
      pc.fHistCharge[eSim]->GetXaxis()->SetTitle("charge");
      pc.fHistTPCncls[eSim]->GetXaxis()->SetTitle("TPCNClsFindable");
      pc.fHistTracksdcaXY[eSim]->GetXaxis()->SetTitle("DCA XY");
      pc.fHistTracksdcaZ[eSim]->GetXaxis()->SetTitle("DCA Z");
      pc.fParticleHistogramsList->Add(pc.fHistPt[eSim]);
      pc.fParticleHistogramsList->Add(pc.fHistPhi[eSim]);
      pc.fParticleHistogramsList->Add(pc.fHistCharge[eSim]);
      pc.fParticleHistogramsList->Add(pc.fHistTPCncls[eSim]);
      pc.fParticleHistogramsList->Add(pc.fHistTracksdcaXY[eSim]);
      pc.fParticleHistogramsList->Add(pc.fHistTracksdcaZ[eSim]);

      event.fHistCentr[eSim] = new TH1F("fHistCentr[eSim]", "centrality distribution for simulated particles", nBinscentr, mincentr, maxcentr);
      event.fHistX[eSim] = new TH1F("fHistX[eSim]", "posX distribution for simulated particles", nBinsx, minx, maxx);
      event.fHistY[eSim] = new TH1F("fHistY[eSim]", "posY distribution for simulated particles", nBinsy, miny, maxy);
      event.fHistZ[eSim] = new TH1F("fHistZ[eSim]", "posZ distribution for simulated particles", nBinsz, minz, maxz);
      event.fHistMult[eSim] = new TH1I("fHistMult[eSim]", "mult distribution for simulated particles", nBinsmult, minmult, maxmult);
      event.fHistCentr[eSim]->GetXaxis()->SetTitle("centrality");
      event.fHistX[eSim]->GetXaxis()->SetTitle("x");
      event.fHistY[eSim]->GetXaxis()->SetTitle("y");
      event.fHistZ[eSim]->GetXaxis()->SetTitle("z");
      event.fHistMult[eSim]->GetXaxis()->SetTitle("multiplicity");
      event.fEventHistogramsList->Add(event.fHistCentr[eSim]);
      event.fEventHistogramsList->Add(event.fHistX[eSim]);
      event.fEventHistogramsList->Add(event.fHistY[eSim]);
      event.fEventHistogramsList->Add(event.fHistZ[eSim]);
      event.fEventHistogramsList->Add(event.fHistMult[eSim]);
    }

    // loading weights
    for (const int& run : targetRuns) {
      std::string runStr = std::to_string(run);
      TH1F* histweights = GetHistogramWithWeights(cfFileWithWeights.value.c_str(), runStr.c_str());
      if (!histweights) {
        LOG(fatal) << "Failed to load weights for run " << run;
        return;
      }
      histweights->SetName(Form("histWithEfficiencyCorrections_%d", run));
      wh.fWeightsHistList->Add(histweights);
      wh.weightsmap[run] = histweights;
    }

    qa.fQA = new TH2F("QA_centr", "quality assurance of centrality", nBinscentr, mincentr, maxcentr, nBinscentr, mincentr, maxcentr);
    qa.fQAM_NC = new TH2F("QAM_NC", "quality assurance of mult vs. NContributors", nBinsmult, minmult, maxmult, nBinsncontr, minncontr, maxncontr);
    if (cfQA) {
      if (cfInitsim)
        qa.fQAList->Add(qa.fQA);
      qa.fQAList->Add(qa.fQAM_NC);
    }

    // float quantiles[10] = {0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8};
    float quantiles[10] = {0, 5, 10, 20, 30, 40, 50, 60, 70, 80};
    cor.pv22_centr = new TProfile("pv22", "profile of v_{2}^{2}", 9, quantiles);
    cor.pv32_centr = new TProfile("pv32", "profile of v_{3}^{2}", 9, quantiles);
    cor.pv42_centr = new TProfile("pv42", "profile of v_{4}^{2}", 9, quantiles);
    cor.pfour32_centr = new TProfile("pfour32", "profile of v_{2}^{2}*v_{3}^{2}", 9, quantiles);
    cor.pfour42_centr = new TProfile("pfour42", "profile of v_{2}^{2}*v_{4}^{2}", 9, quantiles);
    cor.pv2_nest = new TProfile("pv2_nest", "profile of v_{2} from nest", 9, quantiles);
    cor.pv3_nest = new TProfile("pv3_nest", "profile of v_{3} from nest", 9, quantiles);
    cor.pv4_nest = new TProfile("pv4_nest", "profile of v_{4} from nest", 9, quantiles);
    cor.pv22_centr->GetYaxis()->SetTitle("v_{2}^{2}");
    cor.pv32_centr->GetYaxis()->SetTitle("v_{3}^{2}");
    cor.pv42_centr->GetYaxis()->SetTitle("v_{4}^{2}");
    cor.pv22_centr->GetXaxis()->SetTitle("centrality");
    cor.pv32_centr->GetXaxis()->SetTitle("centrality");
    cor.pv42_centr->GetXaxis()->SetTitle("centrality");
    cor.pfour32_centr->GetYaxis()->SetTitle("v_{2}^{2}v_{3}^{2}");
    cor.pfour42_centr->GetYaxis()->SetTitle("v_{2}^{2}v_{4}^{2}");
    cor.pfour32_centr->GetXaxis()->SetTitle("centrality");
    cor.pfour42_centr->GetXaxis()->SetTitle("centrality");
    cor.pv2_nest->GetYaxis()->SetTitle("v_{2}");
    cor.pv3_nest->GetYaxis()->SetTitle("v_{3}");
    cor.pv4_nest->GetYaxis()->SetTitle("v_{4}");
    cor.pv2_nest->GetXaxis()->SetTitle("centrality");
    cor.pv3_nest->GetXaxis()->SetTitle("centrality");
    cor.pv4_nest->GetXaxis()->SetTitle("centrality");
    cor.fCorrelationVariablesList->Add(cor.pv22_centr);
    cor.fCorrelationVariablesList->Add(cor.pv32_centr);
    cor.fCorrelationVariablesList->Add(cor.pv42_centr);
    cor.fCorrelationVariablesList->Add(cor.pfour32_centr);
    cor.fCorrelationVariablesList->Add(cor.pfour42_centr);
    if (cfNest) {
      cor.fCorrelationVariablesList->Add(cor.pv2_nest);
      cor.fCorrelationVariablesList->Add(cor.pv3_nest);
      cor.fCorrelationVariablesList->Add(cor.pv4_nest);
    }

    // init of phi hist for different runs
    for (const int& run : targetRuns) {
      std::string histName = "hphi_run_" + std::to_string(run);
      std::string histTitle = "Phi dis for Run " + std::to_string(run);

      TH1F* h = new TH1F(histName.c_str(), histTitle.c_str(), nBinsphi, minphi, maxphi);
      phih.fPhiHistList->Add(h);
      phih.histMap[run] = h;
    }

  } // end of void init(InitContext&) {

  // A) Process only reconstructed data:
  void processRec(CollisionRec const& collision, aod::BCs const&, TracksRec const& tracks)
  {
    // ...

    // *) Steer all analysis steps:
    Steer<eRec>(collision, tracks);
  }
  PROCESS_SWITCH(MultiharmonicCorrelations, processRec, "process only reconstructed data", true); // yes, keep always one process switch "true", so that there is default running version

  // -------------------------------------------

  // B) Process both reconstructed and corresponding MC truth simulated data:
  void processRecSim(CollisionRecSim const& collision, aod::BCs const&, TracksRecSim const& tracks, aod::McParticles const&, aod::McCollisions const&)
  {
    Steer<eRecAndSim>(collision, tracks);
  }
  PROCESS_SWITCH(MultiharmonicCorrelations, processRecSim, "process both reconstructed and corresponding MC truth simulated data", false);

  // -------------------------------------------

  // C) Process only simulated data:
  void processSim(CollisionSim const& /*collision*/, aod::BCs const&, TracksSim const& /*tracks*/)
  {
    // Steer<eSim>(collision, tracks); // TBI 20241105 not ready yet, but I do not really need this one urgently, since RecSim is working, and I need that one for efficiencies...
  }
  PROCESS_SWITCH(MultiharmonicCorrelations, processSim, "process only simulated data", false);

}; // struct MultiharmonicCorrelations {

// *) The final touch:
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<MultiharmonicCorrelations>(cfgc),
  };
} // WorkflowSpec...
