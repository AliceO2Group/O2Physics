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
#include <TComplex.h>
#include <TFile.h>
#include <TGrid.h>
#include <TH1.h>
#include <TH2.h>
#include <TIterator.h>
#include <TList.h>
#include <TObject.h>
#include <TProfile.h>
#include <TString.h>
#include <TSystem.h>

#include <Rtypes.h>

#include <array>
#include <cmath>
#include <cstdint>
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
  eFT0C = 0,
  eFT0M,
  eFV0A,
  eNTPV,
  eCentralityEstimator_N
};

enum EMultiplicityTables {
  eMultTPC = 0,
  eMultFV0M,
  eMultFT0C,
  eMultFT0M,
  eMultNTracksPV,
  eMultiplicityTables_N
};

enum ERecSim {
  eRec = 0,
  eSim,
  eRecAndSim,
  eRecSim_N
};

enum ECuts {
  eBefore = 0,
  eAfter,
  eCuts_N
};

enum ERealMC {
  eReal = 0,
  eMC
};

enum EProcess {
  eProcessRec = 0,
  eProcessRecSim,
  eProcessSim,
  eProcess_N
};

enum EParticleHistograms {
  eHistPt = 0,
  eHistPhi,
  eHistEta,
  eHistCharge,
  eParticleHistograms_N
};

enum EEventHistograms {
  eHistCentrality = 0,
  eHistMultiplicity,
  eHistVertexX,
  eHistVertexY,
  eHistVertexZ,
  eHistImpactParameter,
  eHistReferenceMultiplicity,
  eEventHistograms_N
};

enum EMiscHistograms {
  eRunNumber = 0,
  eMiscHistograms_N
};

enum EWeightsHistograms {
  ePhi = 0,
  ePt,
  eWeightsHistograms_N
};

enum EObservables {
  eProfTwo = 0,
  eObservablesHist_N
};

// *) Lists of names:
static constexpr std::array<const char*, eCentralityEstimator_N> CentralityEstimatorNames = {
  "FT0C",
  "FT0M",
  "FV0A",
  "NTPV"};

static constexpr std::array<const char*, eMultiplicityTables_N> MultiplicityTablesNames = {
  "multTPC",
  "multFV0M",
  "multFT0C",
  "multFT0M",
  "multNTracksPV"};

static constexpr std::array<const char*, eWeightsHistograms_N> WeightsNames = {
  "ePhi",
  "ePt"};

static constexpr std::array<const char*, eCuts_N> CutsNames = {
  "eBefore",
  "eAfter"};

// *) Main task:
struct MultiparticleCorrelationsMei // this name is used in lower-case format to name the TDirectoryFile in AnalysisResults.root
{
  // *) Base TList to hold all output objects:
  TString sBaseListName = "Default list name"; // yes, I declare it separately, because I need it also later in BailOut() function
  OutputObj<TList> fBaseList{sBaseListName.Data(),
                             OutputObjHandlingPolicy::AnalysisObject,
                             OutputObjSourceType::OutputObjSource};

  // *) CCDB:
  Service<ccdb::BasicCCDBManager> ccdb{}; // support for offline callibration data base, not needed for the time being...
  Service<o2::framework::O2DatabasePDG> pdg{};

  // *) Define configurables:
  Configurable<int> centralityEstimator{"centralityEstimator", 0, "centrality estimator: 0=FT0C, 1=FT0M, 2=FV0A, 3=NTPV"};
  Configurable<int> multiplicityTables{"multiplicityTables", 0, "multiplicity tables: 0=multTPC, 1=multFV0M, 2=multFT0C, 3=multFT0M, 4=multNTracksPV"};
  Configurable<bool> cfDryRun{"cfDryRun", false, "book all histos and run without filling and calculating anything"};

  // *) external root files
  Configurable<bool> cfExternalFileSwitch{"cfExternalFileSwitch", false, "choose to include external root files or not"};
  Configurable<std::string> cfFileWithWeights{"cfFileWithWeights", "/alice-ccdb.cern.ch/Users/m/mei/thesis-", "path to external ROOT file which holds all particle weights"};

  // *) binnings
  Configurable<bool> cfALICECentBinSwitch{"cfALICECentBinSwitch", true, "switch on or off to use ALICE default binning"};
  Configurable<std::vector<float>> cfCentBins{"cfCentBins", {100, 0., 100.}, "nCentBins, centMin, centMax"};
  Configurable<std::vector<float>> cfMultBins{"cfMultBins", {400, 0., 40000.}, "Multiplicity bins: nMultBins, multMin, multMax"};
  Configurable<std::vector<float>> cfMultBinsRef{"cfMultBinsRef", {400, 0., 40000.}, "Reference mult bins: nMultBins, multMin, multMax"};
  Configurable<std::vector<float>> cfVxBins{"cfVxBins", {300, -0.04, 0.04}, "Vertex X hist: nVxBins, vxMin, vxMax"};
  Configurable<std::vector<float>> cfVyBins{"cfVyBins", {300, -0.01, 0.01}, "Vertex Y hist: nVyBins, vyMin, vyMax"};
  Configurable<std::vector<float>> cfVzBins{"cfVzBins", {300, -20., 20.}, "Vertex Z hist: nVzBins, vzMin, vzMax"};
  Configurable<std::vector<float>> cfIpBins{"cfIpBins", {100, 0., 20.}, "Impact parameters hist (MC only): nIPBins, ipMin, ipMax"};

  Configurable<std::vector<float>> cfPtBins{"cfPtBins", {2000, 0., 5.}, "nPtBins, ptMin, ptMax"};
  Configurable<std::vector<float>> cfPhiBins{"cfPhiBins", {180, 0., math::TwoPI}, "nPhiBins, phiMin, phiMax"};
  Configurable<std::vector<float>> cfEtaBins{"cfEtaBins", {800, -3., 3.}, "nEtaBins, etaMin, etaMax"};

  // *) Cuts
  // event level cuts
  Configurable<bool> cfMasterCutSwitch{"cfMasterCutSwitch", true, "switch on or off all cuts"};
  Configurable<bool> cfEventCutSwitch{"cfEventCutSwitch", true, "switch to apply vertex z position cut"};
  Configurable<std::vector<float>> cfVertexZCutRange{"cfVertexZCutRange", {-10, 10.}, "vertex z position range: {min, max}[cm], with convention: min <= Vz < max"};

  // particle level cuts
  Configurable<bool> cfPtCutSwitch{"cfPtCutSwitch", true, "switch to apply pt cut"};
  Configurable<std::vector<float>> cfPtCutRange{"cfPtCutRange", {0.2, 5.}, "pt cut range: {min, max}, with convention: min <= pt < max"};
  Configurable<bool> cfEtaCutSwitch{"cfEtaCutSwitch", true, "switch to apply pt cut"};
  Configurable<std::vector<float>> cfEtaCutRange{"cfEtaCutRange", {-0.8, 0.8}, "eta cut range: {min, max}, with convention: min <= eta < max"};
  Configurable<bool> cfChargeCutSwitch{"cfChargeCutSwitch", true, "switch to apply charge cut (cut neutral particle out)"};

  // *) misc
  Configurable<double> sigmaInel{"sigmaInel", 7.71, "inelastic cross section in mb"};
  Configurable<bool> qualityAssuranceSwitch{"qualityAssuranceSwitch", false, "quality assurance switch"};

  // *) Define and initialize all data members to be called in the main process* functions:
  // **) Task configuration:
  struct TaskConfiguration {
    std::array<bool, eProcess_N> fProcess{false}; // Set what to process. See enum EProcess for full description. Set via implicit variables within a PROCESS_SWITCH clause.
    bool fDryRun = false;                         // book all histos and run without filling and calculating anything
  } tc;                                           // you have to prepend "tc." for all objects name in this group later in the code

  // **) Particle histograms:
  struct ParticleHistograms {
    TList* fParticleHistList = nullptr; //!<! list to hold all control particle histograms
    std::array<std::array<std::array<TH1F*, eCuts_N>, 2>, eParticleHistograms_N> fParticleHist{};
  } pc;

  // *) Event histograms:
  struct EventHistograms {
    TList* fEventHistList = nullptr;
    std::array<std::array<std::array<TH1F*, eCuts_N>, 2>, eEventHistograms_N> fEventHist{}; //! [ type - see enum EEventHistograms ][reco,sim][before, after event cuts]
  } ec;

  struct MiscHistograms {
    TList* fMiscHistList = nullptr;
    TH1D* fMiscHistRunNumber = nullptr;
  } misc;

  // *) External histograms:
  struct ExternalHistograms {
    TList* fExternalHistogramsList = nullptr;
    std::array<std::array<TH1F*, eCuts_N>, eWeightsHistograms_N> fWeights{}; //! [type][before, after cuts]
  } ex;

  struct Observables {
    TList* fObservablesList = nullptr;
    std::array<std::array<TProfile*, eCuts_N>, 2> fProfTwo{}; //! [reco,sim][before, after event cuts]
  } obs;

  // *) Quality assurance histograms:
  struct QualityAssurance {
    TList* fQualityAssuranceList = nullptr; //!<! list to hold all qualityAssurance histograms
    TH2F* fHistCentralityRecSim = nullptr;
  } qa;

  // *) functions
  std::vector<std::vector<TComplex>> initQVectorsTable(int maxCorrelator, const std::vector<int>& harmonic)
  {
    int sum = 0;
    for (const int& x : harmonic) {
      sum += std::abs(x);
    }
    const int maxHarmonic = sum + 1;        // rows
    const int maxPower = maxCorrelator + 1; // cols

    return std::vector<std::vector<TComplex>>(maxHarmonic, std::vector<TComplex>(maxPower, TComplex(0., 0.)));
  }

  template <typename T1>
  void updateQVectorsTable(std::vector<std::vector<TComplex>>& QVectorsTable, T1 const& phi, T1 weight = T1(1))
  {
    const int maxHarmonic = QVectorsTable.size();
    const int maxPower = QVectorsTable.empty() ? 0 : QVectorsTable[0].size();
    for (int h = 0; h < maxHarmonic; ++h) {
      for (int p = 0; p < maxPower; ++p) {
        const auto wp = std::pow(weight, p);
        QVectorsTable[h][p] += TComplex(wp * std::cos(h * phi), wp * std::sin(h * phi));
      }
    }
  }

  std::vector<TComplex> two(const std::vector<std::vector<TComplex>>& QVectorsTable, const std::vector<int>& harmonic)
  {
    int n1 = harmonic[0], n2 = harmonic[1];
    auto q = [&](int n, int p) -> TComplex {
      if (n >= 0) {
        return QVectorsTable[n][p];
      }
      return TComplex::Conjugate(QVectorsTable[-n][p]);
    };
    auto corr = [&](int n1, int n2) -> TComplex {
      return q(n1, 1) * q(n2, 1) - q(n1 + n2, 2);
    };

    TComplex n = corr(n1, n2);
    TComplex d = corr(0, 0);
    return {n, d};
  }

  std::vector<TComplex> recursion(int m, const std::vector<std::vector<TComplex>>& Qvector, std::vector<int>& harmonic, int mult = 1, int skip = 0)
  {
    // Calculate multi-particle correlators by using recursion (an improved faster version) originally developed by
    // Kristjan Gulbrandsen (gulbrand@nbi.dk).

    auto q = [&](int n, int p) -> TComplex {
      if (n >= 0) {
        return Qvector[n][p];
      }
      return TComplex::Conjugate(Qvector[-n][p]);
    };
    std::vector<TComplex> c = {q(harmonic[m - 1], mult), q(0, mult)};
    if ((m - 1) == 0) {
      return c;
    }
    std::vector<TComplex> temp = recursion(m - 1, Qvector, harmonic);
    c[0] *= temp[0];
    c[1] *= temp[1];
    if ((m - 1) == skip) {
      return c;
    }

    int counter1 = 0;
    int hhold = harmonic[counter1];
    harmonic[counter1] = harmonic[m - 2];
    harmonic[m - 2] = hhold + harmonic[m - 1];
    std::vector<TComplex> c2 = recursion(m - 1, Qvector, harmonic, mult + 1, m - 2);
    int counter2 = m - 3;
    while (counter2 >= skip) {
      harmonic[m - 2] = harmonic[counter1];
      harmonic[counter1] = hhold;
      ++counter1;
      hhold = harmonic[counter1];
      harmonic[counter1] = harmonic[m - 2];
      harmonic[m - 2] = hhold + harmonic[m - 1];
      temp = recursion(m - 1, Qvector, harmonic, mult + 1, counter2);
      c2[0] += temp[0];
      c2[1] += temp[1];
      --counter2;
    }
    harmonic[m - 2] = harmonic[counter1];
    harmonic[counter1] = hhold;

    if (mult == 1) {
      return {c[0] - c2[0], c[1] - c2[1]};
    }
    return {c[0] - static_cast<double>(mult) * c2[0], c[1] - static_cast<double>(mult) * c2[1]};
  }

  bool noneZeroDenom(std::vector<TComplex> resultMultCorr)
  {
    return resultMultCorr[1].Re() != 0;
  }

  TObject* getObjectFromList(TList* list, const char* objectName)
  {
    // Get TObject pointer from TList, even if it's in some nested TList. Foreseen
    // to be used to fetch histograms or profiles from files directly.
    // Some ideas taken from TCollection::ls()
    // If you have added histograms directly to files (without TList's), then you can fetch them directly with
    // file->Get("hist-name").

    // Usage: TH1F *hist = (TH1F*) getObjectFromList("some-valid-TList-pointer","some-object-name");

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
      if (auto* subList = dynamic_cast<TList*>(objectIter)) {
        objectFinal = getObjectFromList(subList, objectName);
        if (objectFinal) {
          return objectFinal;
        }
      }
    } // while(objectIter = next())

    return nullptr;
  } // TObject* getObjectFromList(TList *list, char *objectName)

  TH1F* getHistogramWithWeights(const char* filePath, const char* runNumber, const char* histName)
  {
    // *) Return value:
    TH1F* hist = nullptr;
    TList* baseList = nullptr;     // base top-level list in the TFile, e.g. named "ccdb_object"
    TList* listWithRuns = nullptr; // nested list with run-wise TList's holding run-specific weights

    // *) Determine from filePath if the file is on a local machine, or in home dir AliEn, or in CCDB:
    // Algorithm: If filePath begins with "/alice/cern.ch/" then it's in the home dir AliEn; If filePath begins with "/alice-ccdb.cern.ch/" then it's in CCDB. Therefore, files in AliEn and CCDB must be specified with abs path; for local files both abs and relative paths are just fine.
    bool bFileIsInAliEn = false;
    bool bFileIsInCCDB = false;

    std::string path(filePath);

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
        LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
      }

      // Finally, from the top-level TList, get the desired nested TList => the technical problem here is that it can be nested at any level, for that there is a helper utility function getObjectFromList(...) , see its implementation further below
      listWithRuns = dynamic_cast<TList*>(getObjectFromList(baseList, runNumber));
      if (!listWithRuns) {
        TString runNumberWithLeadingZeroes = "000";
        runNumberWithLeadingZeroes += runNumber; // another try, with "000" prepended to run number
        listWithRuns = dynamic_cast<TList*>(getObjectFromList(baseList, runNumberWithLeadingZeroes.Data()));
        if (!listWithRuns) {
          LOGF(error, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
          return nullptr;
        }
      }

      // OK, we got the desired TList with efficiency corrections, after that we can use the common code for all 3 cases (local, AliEn, CCDB, that common code is below)
    } else if (bFileIsInCCDB) {
      // File you want to access is in your home dir in CCDB:
      // Remember that here I do not access the file; instead, I directly access the object in that file.
      ccdb->setURL("https://alice-ccdb.cern.ch");
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
          LOGF(error, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
          return nullptr;
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
        LOGF(fatal, "\033[1;31m%s at line %d\033[0m can't open file", __FUNCTION__, __LINE__);
      }
      weightsFile->GetObject("ccdb_object", baseList);
      if (!baseList) {
        LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
      }

      listWithRuns = dynamic_cast<TList*>(getObjectFromList(baseList, runNumber));
      if (!listWithRuns) {
        TString runNumberWithLeadingZeroes = "000";
        runNumberWithLeadingZeroes += runNumber; // another try, with "000" prepended to run number
        listWithRuns = dynamic_cast<TList*>(getObjectFromList(baseList, runNumberWithLeadingZeroes.Data()));
        if (!listWithRuns) {
          LOGF(error, "\033[1;31m%s at line %d : this crash can happen if in the output file there is no list with weights for the current run number = %s\033[0m", __FUNCTION__, __LINE__, runNumber);
          return nullptr;
        }
      }
    }

    // Here comes the common code for all three cases, where from "listWithRuns" you fetch the desired histogram with efficiency corrections:
    if (!listWithRuns) {
      LOGF(fatal, "\033[1;31m%s: listWithRuns is null for run %s\033[0m", __FUNCTION__, runNumber);
      return nullptr;
    }

    hist = dynamic_cast<TH1F*>(listWithRuns->FindObject(histName));
    if (!hist) {
      LOGF(info, "%s: histogram 'hist' not found in run list", __FUNCTION__);
      return nullptr;
    }

    auto histClone = dynamic_cast<TH1F*>(hist->Clone());
    if (!histClone) {
      LOGF(error, "%s: histogram 'histClone' failed to be cloned", __FUNCTION__);
      return nullptr;
    }
    histClone->SetDirectory(nullptr);

    delete baseList;

    return histClone;
  } // end of TH1F* getHistogramWithWeights(const char* filePath, const char* runNumber, const char* histName)

  // templates
  template <typename T1>
  float chooseCent(T1 const& collision, int const& whichEstimator)
  {
    switch (whichEstimator) {
      case eFT0C:
        return static_cast<float>(collision.centFT0C());
      case eFT0M:
        return static_cast<float>(collision.centFT0M());
      case eFV0A:
        return static_cast<float>(collision.centFV0A());
      case eNTPV:
        return static_cast<float>(collision.centNTPV());
      default:
        LOG(warning) << "Unknown centrality estimator. Using FT0C as default.";
        return static_cast<float>(collision.centFT0C());
    }
  }

  template <typename T1>
  float chooseMult(T1 const& collision, int const& whichTable)
  {
    switch (whichTable) {
      case eMultTPC:
        return static_cast<float>(collision.multTPC());
      case eMultFV0M:
        return static_cast<float>(collision.multFV0M());
      case eMultFT0C:
        return static_cast<float>(collision.multFT0C());
      case eMultFT0M:
        return static_cast<float>(collision.multFT0M());
      case eMultNTracksPV:
        return static_cast<float>(collision.multNTracksPV());
      default:
        LOG(warning) << "Unknown multiplicity. Using multTPC as default.";
        return static_cast<float>(collision.multTPC());
    }
  }

  template <ERecSim rs, ERealMC rm, typename T1>
  bool eventCuts(T1 const& collision)
  {
    if constexpr (rs == eRec || rs == eRecAndSim) {
      if (cfEventCutSwitch) // event level cuts for Rec
      {
        if (rm == eReal) {
          if (collision.posZ() > cfVertexZCutRange.value[1] || collision.posZ() < cfVertexZCutRange.value[0]) {
            return false;
          } // vertex z cut
        }
        if constexpr (rs == eRecAndSim) // event level cuts for Sim
        {
          if (rm == eMC) {
            if (!collision.has_mcCollision()) {
              return false;
            }
            auto thisMCCollision = collision.mcCollision(); // corresponding MC truth simulated particle
            if (thisMCCollision.posZ() > cfVertexZCutRange.value[1] || thisMCCollision.posZ() < cfVertexZCutRange.value[0]) {
              return false;
            } // vertex z cut
          }
        }
      }
    }

    return true;
  }

  template <ERecSim rs, ERealMC rm, ECuts cuts, typename T1, typename T2>
  void eventHistFill(T1 const& collision, T2 const& tracks)
  {
    if constexpr (rs == eRec || rs == eRecAndSim) {
      // Fill reconstructed-level event histograms
      if (rm == eReal) {
        auto thisCent = chooseCent(collision, centralityEstimator);
        auto thisRefMult = chooseMult(collision, multiplicityTables);
        int multiplicityRec = static_cast<int>(tracks.size());
        if constexpr (cuts == eBefore) {
          ec.fEventHist[eHistMultiplicity][eRec][eBefore]->Fill(multiplicityRec);
          ec.fEventHist[eHistCentrality][eRec][eBefore]->Fill(thisCent);
          ec.fEventHist[eHistReferenceMultiplicity][eRec][eBefore]->Fill(thisRefMult);
          ec.fEventHist[eHistVertexX][eRec][eBefore]->Fill(collision.posX());
          ec.fEventHist[eHistVertexY][eRec][eBefore]->Fill(collision.posY());
          ec.fEventHist[eHistVertexZ][eRec][eBefore]->Fill(collision.posZ());
        }
        if constexpr (cuts == eAfter) {
          ec.fEventHist[eHistMultiplicity][eRec][eAfter]->Fill(multiplicityRec);
          ec.fEventHist[eHistCentrality][eRec][eAfter]->Fill(thisCent);
          ec.fEventHist[eHistReferenceMultiplicity][eRec][eAfter]->Fill(thisRefMult);
          ec.fEventHist[eHistVertexX][eRec][eAfter]->Fill(collision.posX());
          ec.fEventHist[eHistVertexY][eRec][eAfter]->Fill(collision.posY());
          ec.fEventHist[eHistVertexZ][eRec][eAfter]->Fill(collision.posZ());
        }
      }

      // Fill MC simulated-level event histograms if both reconstructed and simulated data are processed
      if constexpr (rs == eRecAndSim) {
        if (!collision.has_mcCollision()) {
          // LOGF(warning, "No MC collision for this collision, skip...");
          return;
        }
        if (rm == eMC) {
          auto thisMCCollision = collision.mcCollision(); // corresponding MC truth simulated particle
          int multiplicitySim = static_cast<int>(tracks.size());
          auto impactParameter = thisMCCollision.impactParameter();
          auto centralityMC = math::PI * impactParameter * impactParameter / sigmaInel; // centrality for sim derived from impact parameter
          if constexpr (cuts == eBefore) {
            ec.fEventHist[eHistMultiplicity][eSim][eBefore]->Fill(multiplicitySim);
            ec.fEventHist[eHistCentrality][eSim][eBefore]->Fill(centralityMC);
            ec.fEventHist[eHistImpactParameter][eSim][eBefore]->Fill(impactParameter);
            ec.fEventHist[eHistVertexX][eSim][eBefore]->Fill(thisMCCollision.posX());
            ec.fEventHist[eHistVertexY][eSim][eBefore]->Fill(thisMCCollision.posY());
            ec.fEventHist[eHistVertexZ][eSim][eBefore]->Fill(thisMCCollision.posZ());
          }

          if constexpr (cuts == eAfter) {
            ec.fEventHist[eHistMultiplicity][eSim][eAfter]->Fill(multiplicitySim);
            ec.fEventHist[eHistCentrality][eSim][eAfter]->Fill(centralityMC);
            ec.fEventHist[eHistImpactParameter][eSim][eAfter]->Fill(impactParameter);
            ec.fEventHist[eHistVertexX][eSim][eAfter]->Fill(thisMCCollision.posX());
            ec.fEventHist[eHistVertexY][eSim][eAfter]->Fill(thisMCCollision.posY());
            ec.fEventHist[eHistVertexZ][eSim][eAfter]->Fill(thisMCCollision.posZ());
          }
        } // end of if (rm == eMC) {
      }
    }
  }

  template <ERecSim rs, ERealMC rm, typename T>
  bool particleCuts(T const& track)
  {
    if constexpr (rs == eRec || rs == eRecAndSim) {
      if (cfPtCutSwitch) // pt cuts for Rec
      {
        if (rm == eReal) {
          if (track.pt() < cfPtCutRange.value[0] || track.pt() > cfPtCutRange.value[1]) {
            return false;
          }
        }
        if constexpr (rs == eRecAndSim) // pt cuts for Sim
        {
          if (rm == eMC) {
            if (!track.has_mcParticle()) {
              return false;
            }
            auto mcParticle = track.mcParticle();
            if (mcParticle.pt() < cfPtCutRange.value[0] || mcParticle.pt() > cfPtCutRange.value[1]) {
              return false;
            }
          } // end of if (rm == eMC) {
        }
      }
      if (cfEtaCutSwitch) // eta cuts for Rec
      {
        if (rm == eReal) {
          if (track.eta() < cfEtaCutRange.value[0] || track.eta() > cfEtaCutRange.value[1]) {
            return false;
          }
        }
        if constexpr (rs == eRecAndSim) // eta cuts for Sim
        {
          if (rm == eMC) {
            if (!track.has_mcParticle()) {
              return false;
            }
            auto mcParticle = track.mcParticle();
            if (mcParticle.eta() < cfEtaCutRange.value[0] || mcParticle.eta() > cfEtaCutRange.value[1]) {
              return false;
            }
          } // end of if (rm == eMC) {
        }
      }
      if (cfChargeCutSwitch) // charge cuts for Rec
      {
        if (rm == eReal) {
          if (track.sign() == 0) {
            return false;
          }
        }
        if constexpr (rs == eRecAndSim) // charge cuts for Sim
        {
          if (rm == eMC) {
            if (!track.has_mcParticle()) {
              return false;
            }
            auto mcParticle = track.mcParticle();
            auto chargeMC = pdg->GetParticle(mcParticle.pdgCode())->Charge();
            if (chargeMC == 0) {
              return false;
            }
          } // end of if (rm == eMC) {
        }
      } // end of if (cfChargeCutSwitch) // charge cuts for Rec
    } // end of if constexpr (rs == eRec || rs == eRecAndSim) {

    return true;
  }

  template <ERecSim rs, ERealMC rm, ECuts cuts, typename T1>
  void particleHistFill(T1 const& track)
  {
    if constexpr (rs == eRec || rs == eRecAndSim) {
      if (rm == eReal) {
        if constexpr (cuts == eBefore) {
          pc.fParticleHist[eHistPt][eRec][eBefore]->Fill(track.pt());
          pc.fParticleHist[eHistPhi][eRec][eBefore]->Fill(track.phi());
          pc.fParticleHist[eHistEta][eRec][eBefore]->Fill(track.eta());
          pc.fParticleHist[eHistCharge][eRec][eBefore]->Fill(track.sign());
        }

        if constexpr (cuts == eAfter) {
          pc.fParticleHist[eHistPt][eRec][eAfter]->Fill(track.pt());
          pc.fParticleHist[eHistPhi][eRec][eAfter]->Fill(track.phi());
          pc.fParticleHist[eHistEta][eRec][eAfter]->Fill(track.eta());
          pc.fParticleHist[eHistCharge][eRec][eAfter]->Fill(track.sign());
        }
      }

      if constexpr (rs == eRecAndSim) {
        if (rm == eMC) {
          if (!track.has_mcParticle()) {
            // LOGF(warning, "  No MC particle for this track, skip...");
            return;
          }
          auto mcParticle = track.mcParticle();
          auto chargeMC = pdg->GetParticle(mcParticle.pdgCode())->Charge();
          const int chargeMCUnit = 3;
          chargeMC /= chargeMCUnit;
          if constexpr (cuts == eBefore) {
            pc.fParticleHist[eHistPt][eSim][eBefore]->Fill(mcParticle.pt());
            pc.fParticleHist[eHistPhi][eSim][eBefore]->Fill(mcParticle.phi());
            pc.fParticleHist[eHistEta][eSim][eBefore]->Fill(mcParticle.eta());
            pc.fParticleHist[eHistCharge][eSim][eBefore]->Fill(chargeMC);
          }

          if constexpr (cuts == eAfter) {
            pc.fParticleHist[eHistPt][eSim][eAfter]->Fill(mcParticle.pt());
            pc.fParticleHist[eHistPhi][eSim][eAfter]->Fill(mcParticle.phi());
            pc.fParticleHist[eHistEta][eSim][eAfter]->Fill(mcParticle.eta());
            pc.fParticleHist[eHistCharge][eSim][eAfter]->Fill(chargeMC);
          }
        } // end of if (rm == eMC) {
      }
    }
  }

  void loadWeights(int runNumber)
  {
    for (int i = 0; i < eWeightsHistograms_N; ++i) {
      for (int j = 0; j < eCuts_N; ++j) {
        ex.fWeights[i][j] = getHistogramWithWeights(cfFileWithWeights.value.c_str(), Form("%d", runNumber), Form("[%s][%s]", WeightsNames[i], CutsNames[j]));
        if (!ex.fWeights[i][j]) {
          LOGF(info, "[%s][%s] not found", WeightsNames[i], CutsNames[j]);
          continue;
        }
        ex.fExternalHistogramsList->Add(ex.fWeights[i][j]);
      }
    }
  }

  template <ERecSim rs, typename T1>
  void qaFill(T1 const& collision)
  {
    auto thisCent = chooseCent(collision, centralityEstimator);
    if constexpr (rs == eRecAndSim || rs == eSim) {
      if (!collision.has_mcCollision()) {
        return;
      }
      auto thisMCCollision = collision.mcCollision();
      auto impactParameter = thisMCCollision.impactParameter();
      auto centralityMC = math::PI * impactParameter * impactParameter / sigmaInel; // centrality for sim derived from impact parameter
      qa.fHistCentralityRecSim->Fill(thisCent, centralityMC);
    }
  }

  // *) Misc:
  bool isFirstCollision = true; // this is used to ensure that weights hist are booked once and only once, otherwise, a crash

  // *) Define all member functions to be called in the main process* functions:
  template <ERecSim rs, typename T1, typename T2>
  void steer(T1 const& collision, T2 const& tracks)
  {
    // Dry run:
    if (tc.fDryRun) {
      return;
    }
    const int thisRunNumber = collision.bc().runNumber();
    LOGF(info, "Successfully running, run number is %d", thisRunNumber);

    if (isFirstCollision) {
      // Get run number
      misc.fMiscHistRunNumber->SetTitle(Form("%d", thisRunNumber));
      misc.fMiscHistRunNumber->SetBinContent(1, thisRunNumber);
      // Get weights
      if (cfExternalFileSwitch) {
        loadWeights(thisRunNumber);
      }
    }

    // Fill Quality Assurance
    if (qualityAssuranceSwitch) {
      qaFill<rs>(collision);
    }

    const bool passesEventCutsReal = eventCuts<rs, eReal>(collision);
    const bool passesEventCutsMC = eventCuts<rs, eMC>(collision);

    // Fill Event Hist
    eventHistFill<rs, eReal, eBefore>(collision, tracks);
    eventHistFill<rs, eMC, eBefore>(collision, tracks);

    if (cfMasterCutSwitch) {
      if (passesEventCutsReal) {
        eventHistFill<rs, eReal, eAfter>(collision, tracks);
      }
      if (passesEventCutsMC) {
        eventHistFill<rs, eMC, eAfter>(collision, tracks);
      }
    }

    std::vector<int> n2 = {-2, 2};
    auto qVectorsTableBeforeCutsReal = initQVectorsTable(2, n2);
    auto qVectorsTableAfterCutsReal = initQVectorsTable(2, n2);

    // Main loop over particles:
    auto track = tracks.iteratorAt(0); // set the type and scope from one instance
    for (int64_t i = 0; i < tracks.size(); i++) {
      track = tracks.iteratorAt(i);
      float thisPhi = track.phi();
      float thisPt = track.pt();
      std::array<float, 2> thisPhiAndPt = {thisPhi, thisPt};
      std::array<std::array<float, 2>, 2> thisWeights = {{{1.f, 1.f}, {1.f, 1.f}}}; // {{wPhiBefore, wPhiAfter},{wPtBefore,  wPtAfter}}

      // eWeightsHistograms_N
      if (cfExternalFileSwitch) {
        for (int i = 0; i < eWeightsHistograms_N; ++i) {
          for (int j = 0; j < eCuts_N; ++j) {
            if (ex.fWeights[i][j]) {
              thisWeights[i][j] = ex.fWeights[i][j]->GetBinContent(ex.fWeights[i][j]->FindBin(thisPhiAndPt[i]));
            }
          }
        }
      }

      particleHistFill<rs, eReal, eBefore>(track);
      particleHistFill<rs, eMC, eBefore>(track);
      updateQVectorsTable(qVectorsTableBeforeCutsReal, thisPhi, thisWeights[ePhi][eBefore] * thisWeights[ePt][eBefore]);

      if (cfMasterCutSwitch) {
        if (passesEventCutsReal && particleCuts<rs, eReal>(track)) {
          particleHistFill<rs, eReal, eAfter>(track);
          updateQVectorsTable(qVectorsTableAfterCutsReal, thisPhi, thisWeights[ePhi][eAfter] * thisWeights[ePt][eAfter]);
        }
        if (passesEventCutsMC && particleCuts<rs, eMC>(track)) {
          particleHistFill<rs, eMC, eAfter>(track);
        }
      }
    } // end of for (int64_t i = 0; i < tracks.size(); i++) {

    std::vector<TComplex> resultMultCorr(2, TComplex(0., 0.));
    resultMultCorr = two(qVectorsTableBeforeCutsReal, n2);
    if (noneZeroDenom(resultMultCorr)) {
      obs.fProfTwo[eRec][eBefore]->Fill(0.5, (resultMultCorr[0] / resultMultCorr[1].Re()).Re() / (1e-2));
    }
    if (cfMasterCutSwitch && passesEventCutsReal) {
      resultMultCorr = two(qVectorsTableAfterCutsReal, n2);
      if (noneZeroDenom(resultMultCorr)) {
        obs.fProfTwo[eRec][eAfter]->Fill(0.5, (resultMultCorr[0] / resultMultCorr[1].Re()).Re() / (1e-2));
      }
    }

    // Now the first collision ends
    isFirstCollision = false;
  } // end of template <ERecSim rs, typename T1, typename T2> void steer(T1 const& collision, T2 const& tracks) {

  // *) Initialize and book all objects:
  void init(InitContext&)
  {
    // *) Set colors
    const int afterCutFillColor = kGreen - 10;
    const int afterCutLineColor = kGreen;
    const int beforeCutFillColor = kRed - 10;
    const int beforeCutLineColor = kRed;

    // *) Set automatically what to process, from an implicit variable "doprocessSomEProcessName" within a PROCESS_SWITCH clause:
    tc.fProcess[eProcessRec] = doprocessRec;
    tc.fProcess[eProcessRecSim] = doprocessRecSim;
    tc.fProcess[eProcessSim] = doprocessSim;

    // *) Configure your task using configurables in the json file:
    tc.fDryRun = cfDryRun;

    // *) Book base list:
    auto* temp = new TList();
    temp->SetOwner(true);
    fBaseList.setObject(temp);

    // *) Book and nest all other TLists:
    if (cfExternalFileSwitch) {
      // *) Book External Hist List
      ex.fExternalHistogramsList = new TList();
      ex.fExternalHistogramsList->SetName("ExternalHistograms");
      ex.fExternalHistogramsList->SetOwner(true);
      fBaseList->Add(ex.fExternalHistogramsList);
    }

    // *) Book MiscHistograms TLists:
    misc.fMiscHistList = new TList();
    misc.fMiscHistList->SetName("MiscHistograms");
    misc.fMiscHistList->SetOwner(true);
    fBaseList->Add(misc.fMiscHistList);

    misc.fMiscHistRunNumber = new TH1D("eRunNumber", "to be set", 1, 0, 1);
    misc.fMiscHistList->Add(misc.fMiscHistRunNumber);

    // *) Book particle TLists:
    pc.fParticleHistList = new TList();
    pc.fParticleHistList->SetName("ParticleHistograms");
    pc.fParticleHistList->SetOwner(true);
    fBaseList->Add(pc.fParticleHistList); // any nested TList in the base TList appears as a subdir in the output ROOT file

    // *) Book pt and phi distribution with binning defined through configurables in the json file:
    std::vector<float> lPtBins = cfPtBins.value;
    const int nBinsPt = static_cast<int>(lPtBins[0]);
    const float minPt = lPtBins[1];
    const float maxPt = lPtBins[2];

    std::vector<float> lPhiBins = cfPhiBins.value;
    const int nBinsPhi = static_cast<int>(lPhiBins[0]);
    const float minPhi = lPhiBins[1];
    const float maxPhi = lPhiBins[2];

    std::vector<float> lEtaBins = cfEtaBins.value;
    const int nBinsEta = static_cast<int>(lEtaBins[0]);
    const float minEta = lEtaBins[1];
    const float maxEta = lEtaBins[2];

    const int nBinsCharge = 5;
    const float minCharge = -2.5;
    const float maxCharge = 2.5;

    if (doprocessRec || doprocessRecSim) {
      pc.fParticleHist[eHistPt][eRec][eBefore] = new TH1F("[eHistPt][eRec][eBefore]", "p_{T} distribution for reconstructed particles before cuts", nBinsPt, minPt, maxPt);
      pc.fParticleHist[eHistPt][eRec][eBefore]->GetXaxis()->SetTitle("p_{T}");

      pc.fParticleHist[eHistPhi][eRec][eBefore] = new TH1F("[eHistPhi][eRec][eBefore]", "#phi distribution for reconstructed particles before cuts", nBinsPhi, minPhi, maxPhi);
      pc.fParticleHist[eHistPhi][eRec][eBefore]->GetXaxis()->SetTitle("#varphi");

      pc.fParticleHist[eHistEta][eRec][eBefore] = new TH1F("[eHistEta][eRec][eBefore]", "#eta distribution for reconstructed particles before cuts", nBinsEta, minEta, maxEta);
      pc.fParticleHist[eHistEta][eRec][eBefore]->GetXaxis()->SetTitle("#eta");

      pc.fParticleHist[eHistCharge][eRec][eBefore] = new TH1F("[eHistCharge][eRec][eBefore]", "charge distribution for reconstructed particles before cuts", nBinsCharge, minCharge, maxCharge);
      pc.fParticleHist[eHistCharge][eRec][eBefore]->GetXaxis()->SetTitle("Particle charge");

      for (int i = 0; i < eParticleHistograms_N; ++i) {
        pc.fParticleHist[i][eRec][eBefore]->SetColors(beforeCutLineColor, -1, beforeCutFillColor);
        pc.fParticleHistList->Add(pc.fParticleHist[i][eRec][eBefore]);
      }

      if (cfMasterCutSwitch) {
        pc.fParticleHist[eHistPt][eRec][eAfter] = new TH1F("[eHistPt][eRec][eAfter]", "p_{T} distribution for reconstructed particles after cuts", nBinsPt, minPt, maxPt);
        pc.fParticleHist[eHistPt][eRec][eAfter]->GetXaxis()->SetTitle("p_{T}");

        pc.fParticleHist[eHistPhi][eRec][eAfter] = new TH1F("[eHistPhi][eRec][eAfter]", "#phi distribution for reconstructed particles after cuts", nBinsPhi, minPhi, maxPhi);
        pc.fParticleHist[eHistPhi][eRec][eAfter]->GetXaxis()->SetTitle("#varphi");

        pc.fParticleHist[eHistEta][eRec][eAfter] = new TH1F("[eHistEta][eRec][eAfter]", "#eta distribution for reconstructed particles after cuts", nBinsEta, minEta, maxEta);
        pc.fParticleHist[eHistEta][eRec][eAfter]->GetXaxis()->SetTitle("#eta");

        pc.fParticleHist[eHistCharge][eRec][eAfter] = new TH1F("[eHistCharge][eRec][eAfter]", "charge distribution for reconstructed particles after cuts", nBinsCharge, minCharge, maxCharge);
        pc.fParticleHist[eHistCharge][eRec][eAfter]->GetXaxis()->SetTitle("Particle charge");

        for (int i = 0; i < eParticleHistograms_N; ++i) {
          pc.fParticleHist[i][eRec][eAfter]->SetColors(afterCutLineColor, -1, afterCutFillColor);
          pc.fParticleHistList->Add(pc.fParticleHist[i][eRec][eAfter]);
        }
      }
    }

    if (doprocessSim || doprocessRecSim) {
      pc.fParticleHist[eHistPt][eSim][eBefore] = new TH1F("[eHistPt][eSim][eBefore]", "p_{T} distribution for simulated particles  before cuts", nBinsPt, minPt, maxPt);
      pc.fParticleHist[eHistPt][eSim][eBefore]->GetXaxis()->SetTitle("p_{T}");

      pc.fParticleHist[eHistPhi][eSim][eBefore] = new TH1F("[eHistPhi][eSim][eBefore]", "#phi distribution for simulated particles before cuts", nBinsPhi, minPhi, maxPhi);
      pc.fParticleHist[eHistPhi][eSim][eBefore]->GetXaxis()->SetTitle("#varphi");

      pc.fParticleHist[eHistEta][eSim][eBefore] = new TH1F("[eHistEta][eSim][eBefore]", "#eta distribution for simulated particles before cuts", nBinsEta, minEta, maxEta);
      pc.fParticleHist[eHistEta][eSim][eBefore]->GetXaxis()->SetTitle("#eta");

      pc.fParticleHist[eHistCharge][eSim][eBefore] = new TH1F("[eHistCharge][eSim][eBefore]", "charge distribution for simulated particles before cuts", nBinsCharge, minCharge, maxCharge);
      pc.fParticleHist[eHistCharge][eSim][eBefore]->GetXaxis()->SetTitle("Particle charge");

      for (int i = 0; i < eParticleHistograms_N; ++i) {
        pc.fParticleHist[i][eSim][eBefore]->SetColors(beforeCutLineColor, -1, beforeCutFillColor);
        pc.fParticleHistList->Add(pc.fParticleHist[i][eSim][eBefore]);
      }

      if (cfMasterCutSwitch) {
        pc.fParticleHist[eHistPt][eSim][eAfter] = new TH1F("[eHistPt][eSim][eAfter]", "p_{T} distribution for simulated particles after cuts", nBinsPt, minPt, maxPt);
        pc.fParticleHist[eHistPt][eSim][eAfter]->GetXaxis()->SetTitle("p_{T}");

        pc.fParticleHist[eHistPhi][eSim][eAfter] = new TH1F("[eHistPhi][eSim][eAfter]", "#phi distribution for simulated particles after cuts", nBinsPhi, minPhi, maxPhi);
        pc.fParticleHist[eHistPhi][eSim][eAfter]->GetXaxis()->SetTitle("#varphi");

        pc.fParticleHist[eHistEta][eSim][eAfter] = new TH1F("[eHistEta][eSim][eAfter]", "#eta distribution for simulated particles after cuts", nBinsEta, minEta, maxEta);
        pc.fParticleHist[eHistEta][eSim][eAfter]->GetXaxis()->SetTitle("#eta");

        pc.fParticleHist[eHistCharge][eSim][eAfter] = new TH1F("[eHistCharge][eSim][eAfter]", "charge distribution for simulated particles after cuts", nBinsCharge, minCharge, maxCharge);
        pc.fParticleHist[eHistCharge][eSim][eAfter]->GetXaxis()->SetTitle("Particle charge");

        for (int i = 0; i < eParticleHistograms_N; ++i) {
          pc.fParticleHist[i][eSim][eAfter]->SetColors(afterCutLineColor, -1, afterCutFillColor);
          pc.fParticleHistList->Add(pc.fParticleHist[i][eSim][eAfter]);
        }
      }
    }

    // Book event-level histograms
    ec.fEventHistList = new TList();
    ec.fEventHistList->SetName("EventHistograms");
    ec.fEventHistList->SetOwner(true);
    fBaseList->Add(ec.fEventHistList);

    float defaultBoundaries[] = {0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100};
    const int nDefaultBins = sizeof(defaultBoundaries) / sizeof(defaultBoundaries[0]) - 1;

    std::vector<float> lCent = cfCentBins.value;
    const int nBinsCent = static_cast<int>(lCent[0]);
    const float minCent = lCent[1];
    const float maxCent = lCent[2];

    std::vector<float> lMult = cfMultBins.value;
    const int nBinsMult = static_cast<int>(lMult[0]);
    const float minMult = lMult[1];
    const float maxMult = lMult[2];

    std::vector<float> lMultRef = cfMultBinsRef.value;
    const int nBinsMultRef = static_cast<int>(lMultRef[0]);
    const float minMultRef = lMultRef[1];
    const float maxMultRef = lMultRef[2];

    std::vector<float> lVx = cfVxBins.value;
    const int nBinsVx = static_cast<int>(lVx[0]);
    const float minVx = lVx[1];
    const float maxVx = lVx[2];

    std::vector<float> lVy = cfVyBins.value;
    const int nBinsVy = static_cast<int>(lVy[0]);
    const float minVy = lVy[1];
    const float maxVy = lVy[2];

    std::vector<float> lVz = cfVzBins.value;
    const int nBinsVz = static_cast<int>(lVz[0]);
    const float minVz = lVz[1];
    const float maxVz = lVz[2];

    std::vector<float> lIp = cfIpBins.value;
    const int nBinsIp = static_cast<int>(lIp[0]);
    const float minIp = lIp[1];
    const float maxIp = lIp[2];

    // eEventHistograms_N

    if (doprocessRec || doprocessRecSim) {
      if (cfALICECentBinSwitch) {
        ec.fEventHist[eHistCentrality][eRec][eBefore] = new TH1F("[eHistCentrality][eRec][eBefore]", "Centrality (reconstructed) before cuts", nDefaultBins, defaultBoundaries);
      } else {
        ec.fEventHist[eHistCentrality][eRec][eBefore] = new TH1F("[eHistCentrality][eRec][eBefore]", "Centrality (reconstructed) before cuts", nBinsCent, minCent, maxCent);
      }
      ec.fEventHist[eHistCentrality][eRec][eBefore]->GetXaxis()->SetTitle(Form("Centrality (%s)", CentralityEstimatorNames[centralityEstimator]));

      ec.fEventHist[eHistMultiplicity][eRec][eBefore] = new TH1F("[eHistMultiplicity][eRec][eBefore]", "Multiplicity (reconstructed) before cuts", nBinsMult, minMult, maxMult);
      ec.fEventHist[eHistMultiplicity][eRec][eBefore]->GetXaxis()->SetTitle("Multiplicity");

      ec.fEventHist[eHistReferenceMultiplicity][eRec][eBefore] = new TH1F("[eHistReferenceMultiplicity][eRec][eBefore]", "Reference Multiplicity before cuts", nBinsMultRef, minMultRef, maxMultRef);
      ec.fEventHist[eHistReferenceMultiplicity][eRec][eBefore]->GetXaxis()->SetTitle(Form("Reference Multiplicity (%s)", MultiplicityTablesNames[multiplicityTables]));

      ec.fEventHist[eHistVertexX][eRec][eBefore] = new TH1F("[eHistVertexX][eRec][eBefore]", "Vertex X (reconstructed) before cuts", nBinsVx, minVx, maxVx);
      ec.fEventHist[eHistVertexX][eRec][eBefore]->GetXaxis()->SetTitle("Vertex X");

      ec.fEventHist[eHistVertexY][eRec][eBefore] = new TH1F("[eHistVertexY][eRec][eBefore]", "Vertex Y (reconstructed) before cuts", nBinsVy, minVy, maxVy);
      ec.fEventHist[eHistVertexY][eRec][eBefore]->GetXaxis()->SetTitle("Vertex Y");

      ec.fEventHist[eHistVertexZ][eRec][eBefore] = new TH1F("[eHistVertexZ][eRec][eBefore]", "Vertex Z (reconstructed) before cuts", nBinsVz, minVz, maxVz);
      ec.fEventHist[eHistVertexZ][eRec][eBefore]->GetXaxis()->SetTitle("Vertex Z");

      for (int i = 0; i < eEventHistograms_N; ++i) {
        if (i != eHistImpactParameter) {
          ec.fEventHist[i][eRec][eBefore]->SetColors(beforeCutLineColor, -1, beforeCutFillColor);
          ec.fEventHistList->Add(ec.fEventHist[i][eRec][eBefore]);
        }
      }

      if (cfMasterCutSwitch) {
        if (cfALICECentBinSwitch) {
          ec.fEventHist[eHistCentrality][eRec][eAfter] = new TH1F("[eHistCentrality][eRec][eAfter]", "Centrality (reconstructed) after cuts", nDefaultBins, defaultBoundaries);
        } else {
          ec.fEventHist[eHistCentrality][eRec][eAfter] = new TH1F("[eHistCentrality][eRec][eAfter]", "Centrality (reconstructed) after cuts", nBinsCent, minCent, maxCent);
        }
        ec.fEventHist[eHistCentrality][eRec][eAfter]->GetXaxis()->SetTitle(Form("Centrality (%s)", CentralityEstimatorNames[centralityEstimator]));

        ec.fEventHist[eHistMultiplicity][eRec][eAfter] = new TH1F("[eHistMultiplicity][eRec][eAfter]", "Multiplicity (reconstructed) after cuts", nBinsMult, minMult, maxMult);
        ec.fEventHist[eHistMultiplicity][eRec][eAfter]->GetXaxis()->SetTitle("Multiplicity");

        ec.fEventHist[eHistReferenceMultiplicity][eRec][eAfter] = new TH1F("[eHistReferenceMultiplicity][eRec][eAfter]", "Reference Multiplicity after cuts", nBinsMultRef, minMultRef, maxMultRef);
        ec.fEventHist[eHistReferenceMultiplicity][eRec][eAfter]->GetXaxis()->SetTitle(Form("Reference Multiplicity (%s)", MultiplicityTablesNames[multiplicityTables]));

        ec.fEventHist[eHistVertexX][eRec][eAfter] = new TH1F("[eHistVertexX][eRec][eAfter]", "Vertex X (reconstructed) after cuts", nBinsVx, minVx, maxVx);
        ec.fEventHist[eHistVertexX][eRec][eAfter]->GetXaxis()->SetTitle("Vertex X");

        ec.fEventHist[eHistVertexY][eRec][eAfter] = new TH1F("[eHistVertexY][eRec][eAfter]", "Vertex Y (reconstructed) after cuts", nBinsVy, minVy, maxVy);
        ec.fEventHist[eHistVertexY][eRec][eAfter]->GetXaxis()->SetTitle("Vertex Y");

        ec.fEventHist[eHistVertexZ][eRec][eAfter] = new TH1F("[eHistVertexZ][eRec][eAfter]", "Vertex Z (reconstructed) after cuts", nBinsVz, minVz, maxVz);
        ec.fEventHist[eHistVertexZ][eRec][eAfter]->GetXaxis()->SetTitle("Vertex Z");

        for (int i = 0; i < eEventHistograms_N; ++i) {
          if (i != eHistImpactParameter) {
            ec.fEventHist[i][eRec][eAfter]->SetColors(afterCutLineColor, -1, afterCutFillColor);
            ec.fEventHistList->Add(ec.fEventHist[i][eRec][eAfter]);
          }
        }
      }
    }

    if (doprocessSim || doprocessRecSim) {
      if (cfALICECentBinSwitch) {
        ec.fEventHist[eHistCentrality][eSim][eBefore] = new TH1F("[eHistCentrality][eSim][eBefore]", "Centrality (simulated) before cuts", nDefaultBins, defaultBoundaries);
      } else {
        ec.fEventHist[eHistCentrality][eSim][eBefore] = new TH1F("[eHistCentrality][eSim][eBefore]", "Centrality (simulated) before cuts", nBinsCent, minCent, maxCent);
      }
      ec.fEventHist[eHistCentrality][eSim][eBefore]->GetXaxis()->SetTitle(Form("Centrality (%s)", CentralityEstimatorNames[centralityEstimator]));

      ec.fEventHist[eHistMultiplicity][eSim][eBefore] = new TH1F("[eHistMultiplicity][eSim][eBefore]", "Multiplicity (simulated) before cuts", nBinsMult, minMult, maxMult);
      ec.fEventHist[eHistMultiplicity][eSim][eBefore]->GetXaxis()->SetTitle("Multiplicity");

      ec.fEventHist[eHistVertexX][eSim][eBefore] = new TH1F("[eHistVertexX][eSim][eBefore]", "Vertex X (simulated) before cuts", nBinsVx, minVx, maxVx);
      ec.fEventHist[eHistVertexX][eSim][eBefore]->GetXaxis()->SetTitle("Vertex X");

      ec.fEventHist[eHistVertexY][eSim][eBefore] = new TH1F("[eHistVertexY][eSim][eBefore]", "Vertex Y (simulated) before cuts", nBinsVy, minVy, maxVy);
      ec.fEventHist[eHistVertexY][eSim][eBefore]->GetXaxis()->SetTitle("Vertex Y");

      ec.fEventHist[eHistVertexZ][eSim][eBefore] = new TH1F("[eHistVertexZ][eSim][eBefore]", "Vertex Z (simulated) before cuts", nBinsVz, minVz, maxVz);
      ec.fEventHist[eHistVertexZ][eSim][eBefore]->GetXaxis()->SetTitle("Vertex Z");

      ec.fEventHist[eHistImpactParameter][eSim][eBefore] = new TH1F("[eHistImpactParameter][eSim][eBefore]", "Impact Parameter (simulated) before cuts", nBinsIp, minIp, maxIp);
      ec.fEventHist[eHistImpactParameter][eSim][eBefore]->GetXaxis()->SetTitle("Impact Parameter");

      for (int i = 0; i < eEventHistograms_N; ++i) {
        if (i != eHistReferenceMultiplicity) {
          ec.fEventHist[i][eSim][eBefore]->SetColors(beforeCutLineColor, -1, beforeCutFillColor);
          ec.fEventHistList->Add(ec.fEventHist[i][eSim][eBefore]);
        }
      }

      if (cfMasterCutSwitch) {
        if (cfALICECentBinSwitch) {
          ec.fEventHist[eHistCentrality][eSim][eAfter] = new TH1F("[eHistCentrality][eSim][eAfter]", "Centrality (simulated) after cuts", nDefaultBins, defaultBoundaries);
        } else {
          ec.fEventHist[eHistCentrality][eSim][eAfter] = new TH1F("[eHistCentrality][eSim][eAfter]", "Centrality (simulated) after cuts", nBinsCent, minCent, maxCent);
        }
        ec.fEventHist[eHistCentrality][eSim][eAfter]->GetXaxis()->SetTitle(Form("Centrality (%s)", CentralityEstimatorNames[centralityEstimator]));

        ec.fEventHist[eHistMultiplicity][eSim][eAfter] = new TH1F("[eHistMultiplicity][eSim][eAfter]", "Multiplicity (simulated) after cuts", nBinsMult, minMult, maxMult);
        ec.fEventHist[eHistMultiplicity][eSim][eAfter]->GetXaxis()->SetTitle("Multiplicity");

        ec.fEventHist[eHistVertexX][eSim][eAfter] = new TH1F("[eHistVertexX][eSim][eAfter]", "Vertex X (simulated) after cuts", nBinsVx, minVx, maxVx);
        ec.fEventHist[eHistVertexX][eSim][eAfter]->GetXaxis()->SetTitle("Vertex X");

        ec.fEventHist[eHistVertexY][eSim][eAfter] = new TH1F("[eHistVertexY][eSim][eAfter]", "Vertex Y (simulated) after cuts", nBinsVy, minVy, maxVy);
        ec.fEventHist[eHistVertexY][eSim][eAfter]->GetXaxis()->SetTitle("Vertex Y");

        ec.fEventHist[eHistVertexZ][eSim][eAfter] = new TH1F("[eHistVertexZ][eSim][eAfter]", "Vertex Z (simulated) after cuts", nBinsVz, minVz, maxVz);
        ec.fEventHist[eHistVertexZ][eSim][eAfter]->GetXaxis()->SetTitle("Vertex Z");

        ec.fEventHist[eHistImpactParameter][eSim][eAfter] = new TH1F("[eHistImpactParameter][eSim][eAfter]", "Impact Parameter (simulated) after cuts", nBinsIp, minIp, maxIp);
        ec.fEventHist[eHistImpactParameter][eSim][eAfter]->GetXaxis()->SetTitle("Impact Parameter");

        for (int i = 0; i < eEventHistograms_N; ++i) {
          if (i != eHistReferenceMultiplicity) {
            ec.fEventHist[i][eSim][eAfter]->SetColors(afterCutLineColor, -1, afterCutFillColor);
            ec.fEventHistList->Add(ec.fEventHist[i][eSim][eAfter]);
          }
        }
      }
    }

    // *) Book observales TLists:
    obs.fObservablesList = new TList();
    obs.fObservablesList->SetName("Observables");
    obs.fObservablesList->SetOwner(true);
    fBaseList->Add(obs.fObservablesList);

    if (doprocessRec || doprocessRecSim) {
      obs.fProfTwo[eRec][eBefore] = new TProfile("obs.fProfTwo[eRec][eBefore]", "obs.fProfTwo[eRec][eBefore]", 1, 0., 1);
      obs.fProfTwo[eRec][eBefore]->GetYaxis()->SetTitle("#LT#LTk#GT#GT / 10^{-k}");
      obs.fObservablesList->Add(obs.fProfTwo[eRec][eBefore]);

      if (cfMasterCutSwitch) {
        obs.fProfTwo[eRec][eAfter] = new TProfile("obs.fProfTwo[eRec][eAfter]", "obs.fProfTwo[eRec][eAfter]", 1, 0., 1);
        obs.fProfTwo[eRec][eAfter]->GetYaxis()->SetTitle("#LT#LTk#GT#GT / 10^{-k}");
        obs.fObservablesList->Add(obs.fProfTwo[eRec][eAfter]);
      }
    }

    // *) Book and QA TLists:
    if (qualityAssuranceSwitch && doprocessRecSim) {
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
    steer<eRec>(collision, tracks);
  }
  PROCESS_SWITCH(MultiparticleCorrelationsMei, processRec, "process only reconstructed data", true); // yes, keep always one process switch "true", so that there is default running version

  // B) Process both reconstructed and corresponding MC truth simulated data:
  void processRecSim(CollisionRecSim const& collision, aod::BCs const&, TracksRecSim const& tracks, aod::McParticles const&, aod::McCollisions const&)
  {
    steer<eRecAndSim>(collision, tracks);
  }
  PROCESS_SWITCH(MultiparticleCorrelationsMei, processRecSim, "process both reconstructed and corresponding MC truth simulated data", false);

  // C) Process only simulated data:
  void processSim(CollisionSim const& /*collision*/, aod::BCs const&, TracksSim const& /*tracks*/)
  {
    // steer<eSim>(collision, tracks); // TBI 20241105 not ready yet, but I do not really need this one urgently, since RecSim is working, and I need that one for efficiencies...
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
