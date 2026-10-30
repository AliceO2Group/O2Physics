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
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <TCollection.h>
#include <TComplex.h>
#include <TF1.h>
#include <TFile.h>
#include <TGrid.h>
#include <TH1.h>
#include <TH2.h>
#include <TIterator.h>
#include <TList.h>
#include <TObject.h>
#include <TParticlePDG.h>
#include <TProfile.h>
#include <TRandom3.h>
#include <TString.h>
#include <TSystem.h>

#include <Rtypes.h>

#include <array>
#include <cmath>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::constants;

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
enum ERecSim {
  eRec = 0,
  eSim,
  eRecAndSim,
  eRecSim_N
};

enum EProcess {
  eProcessRec = 0,
  eProcessRecSim,
  eProcessSim,
  eProcess_N
};

enum ERealMC {
  eReal = 0,
  eMC
};

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

enum EEventCuts {
  eSel8Cut = 0,
  eIsGoodITSLayersAll,
  eNoCollInTimeRangeStandard,
  eNoCollInRofStandard,
  eIsGoodZvtxFT0vsPV,
  eNoHighMultCollInPrevRof,
  eNoSameBunchPileUp,
  eIsVertexITSTPC,
  eCentCut,
  eVertexZCut,
  eNContribCut,
  EEventCuts_N
};

enum EParticleCuts {
  ePtCut = 0,
  eEtaCut,
  eChargeCut,
  eTPCNumClusterCut,
  eDCAxyCut,
  eDCAzCut,
  EParticleCuts_N
};

enum EBeforeAfterCuts {
  eBefore = 0,
  eAfter,
  eBeforeAfterCuts_N
};

enum EEventHistograms {
  eHistCentrality = 0,
  eHistMult,
  eHistRefMult,
  eHistVertexX,
  eHistVertexY,
  eHistVertexZ,
  eHistImpactParameter,
  eHistNumContrib,
  eEventHistograms_N
};

enum EParticleHistograms {
  eHistPt = 0,
  eHistPhi,
  eHistEta,
  eHistCharge,
  eParticleHistograms_N
};

enum EWeightsHistograms {
  eWPhiRec = 0,
  eWPhiSim,
  eWPt,
  eWeightsHistograms_N
};

enum EMiscHistograms {
  eRunNumber = 0,
  eMiscHistograms_N
};

enum EObservables {
  eV22 = 0,
  eV22Gap,
  eV24,
  eObservables_N
};

enum EIntValHarmonics {
  eIV2 = 0,
  eIV3,
  eIV4,
  eIntValHarmonics_N
};

// *) Lists of names:
static constexpr std::array<const char*, eCentralityEstimator_N> CentralityEstimatorNames = {"FT0C", "FT0M", "FV0A", "NTPV"};
static constexpr std::array<const char*, eMultiplicityTables_N> MultiplicityTablesNames = {"multTPC", "multFV0M", "multFT0C", "multFT0M", "multNTracksPV"};

static constexpr std::array<const char*, eParticleHistograms_N> ParticleHistTitles = {"p_{T}", "#varphi", "#eta", "Particle charge"};
static constexpr std::array<const char*, eParticleHistograms_N> ParticleHistNames = {"eHistPt", "eHistPhi", "eHistEta", "eHistCharge"};

static constexpr std::array<const char*, eEventHistograms_N> EventHistTitles = {"Centrality", "Multiplicity (trak.size())", "Reference Multiplicity", "Vertex X", "Vertex Y", "Vertex Z", "Impact Parameter", "Number of Contributor"};
static constexpr std::array<const char*, eEventHistograms_N> EventHistNames = {"eHistCentrality", "eHistMult", "eHistRefMult", "eHistVertexX", "eHistVertexY", "eHistVertexZ", "eHistImpactParameter", "eHistNumContrib"};

static constexpr std::array<const char*, eWeightsHistograms_N> WeightsNames = {"eWPhiRec", "eWPhiSim", "eWPt"};

static constexpr std::array<const char*, eBeforeAfterCuts_N> CutsNames = {"eBefore", "eAfter"};
static constexpr std::array<const char*, eBeforeAfterCuts_N> CutsTitles = {"before", "after"};

static constexpr std::array<const char*, eObservables_N> ObservablesNames = {"eV22", "eV22Gap", "eV24"};
static constexpr std::array<const char*, eObservables_N> ObservablesTitles = {"v_{2}^{2}", "v_{2}^{2} with eta gap", "v_{2}^{4}"};

static constexpr std::array<const char*, eIntValHarmonics_N> IntValNames = {"eIV2", "eIV3", "eIV4"};
static constexpr std::array<const char*, eIntValHarmonics_N> IntValTitles = {"v_{2}^{2}", "v_{3}^{2}", "v_{4}^{2}"};

// *) Main task:
struct MultiparticleCorrelationsMei // this name is used in lower-case format to name the TDirectoryFile in AnalysisResults.root
{
  // *) Base TList to hold all output objects:
  TString sBaseListName = "Default list name"; // yes, I declare it separately, because I need it also later in BailOut() function
  OutputObj<TList> fBaseList{sBaseListName.Data(), OutputObjHandlingPolicy::AnalysisObject, OutputObjSourceType::OutputObjSource};

  // *) CCDB:
  Service<ccdb::BasicCCDBManager> ccdb{}; // support for offline callibration data base
  Service<o2::framework::O2DatabasePDG> pdg{};

  // *) Define configurables:
  Configurable<bool> cfDryRun{"cfDryRun", false, "book all histos and run without filling and calculating anything"};
  Configurable<int> cfCentralityEstimator{"cfCentralityEstimator", 0, "centrality estimator: 0=FT0C, 1=FT0M, 2=FV0A, 3=NTPV"};
  Configurable<int> cfMultiplicityTables{"cfMultiplicityTables", 2, "multiplicity tables: 0=multTPC, 1=multFV0M, 2=multFT0C, 3=multFT0M, 4=multNTracksPV"};

  // *) External root files
  Configurable<bool> cfExternalFileSwitch{"cfExternalFileSwitch", false, "choose to include external root files or not"};
  Configurable<std::string> cfFileWithWeights{"cfFileWithWeights", "/alice-ccdb.cern.ch/Users/m/mei/thesis-", "path to external ROOT file which holds all particle weights"};

  // *) Binnings
  Configurable<std::vector<float>> cfCentBins{"cfCentBins", {100, 0., 100.}, "Centrality bins: nCentBins, centMin, centMax"};
  Configurable<std::vector<float>> cfMultBins{"cfMultBins", {400, 0., 25000.}, "Multiplicity bins: nMultBins, multMin, multMax"};
  Configurable<std::vector<float>> cfRefMultBins{"cfRefMultBins", {400, 0., 25000.}, "Reference mult bins: nMultBins, multMin, multMax"};
  Configurable<std::vector<float>> cfNumContribBins{"cfNumContribBins", {400, 0., 7000.}, "Number of contributors bins: nBinsContrib, contribMin, contribMax"};
  Configurable<std::vector<float>> cfVxBins{"cfVxBins", {500, -0.04, 0.04}, "Vertex X hist: nVxBins, vxMin, vxMax"};
  Configurable<std::vector<float>> cfVyBins{"cfVyBins", {500, -0.015, 0.015}, "Vertex Y hist: nVyBins, vyMin, vyMax"};
  Configurable<std::vector<float>> cfVzBins{"cfVzBins", {500, -20., 20.}, "Vertex Z hist: nVzBins, vzMin, vzMax"};
  Configurable<std::vector<float>> cfIpBins{"cfIpBins", {100, 0., 20.}, "Impact parameters hist (MC only): nIPBins, ipMin, ipMax"};

  Configurable<std::vector<float>> cfPtBins{"cfPtBins", {2000, 0., 6.}, "nPtBins, ptMin, ptMax"};
  Configurable<std::vector<float>> cfPhiBins{"cfPhiBins", {180, 0., math::TwoPI}, "nPhiBins, phiMin, phiMax"};
  Configurable<std::vector<float>> cfEtaBins{"cfEtaBins", {800, -3., 3.}, "nEtaBins, etaMin, etaMax"};

  // *) Cuts
  Configurable<bool> cfMasterCutSwitch{"cfMasterCutSwitch", true, "switch on or off all cuts"};

  // event level cuts
  Configurable<std::vector<std::string>> cfEventCutSwitch{"cfEventCutSwitch", {"1Sel8", "1IsGoodITSLayersAll", "1NoCollInTimeRangeStandard", "1NoCollInRofStandard", "1IsGoodZvtxFT0vsPV", "1NoHighMultCollInPrevRof", "1NoSameBunchPileUp", "1IsVertexITSTPC", "1Cent", "1VertexZ", "1NContrib"}, "event cuts switch, on and off by the first number before name"};
  Configurable<std::vector<float>> cfVertexZCutRange{"cfVertexZCutRange", {-10., 10.}, "vertex z position range: {min, max}[cm], with convention: min <= Vz <= max"};
  Configurable<std::vector<float>> cfCentCutRange{"cfCentCutRange", {0., 80.}, "centrality range: {min, max}, with convention: min <= cent <= max"};
  Configurable<int> cfNContribCutInf{"cfNContribCutInf", 2, "Cuts on number of tracks used for the vertex, only numContrib > cfNContribCutInf survives"};

  // particle level cuts
  Configurable<std::vector<std::string>> cfParticleCutSwitch{"cfParticleCutSwitch", {"1Pt", "1Eta", "1Charge", "1TPCNumCluster", "1DCAxy", "1DCAz"}, "particle cuts switch, on and off by the first number before name"};
  Configurable<std::vector<float>> cfPtCutRange{"cfPtCutRange", {0.2, 5.}, "pt cut range: {min, max}[GeV/c], with convention: min <= pt <= max"};
  Configurable<std::vector<float>> cfEtaCutRange{"cfEtaCutRange", {-0.8, 0.8}, "eta cut range: {min, max}, with convention: min <= eta <= max"};
  Configurable<std::vector<int>> cfTPCNumClusterCutRange{"cfTPCNumClusterCutRange", {70, 160}, "range of cut on the number of clusters associated with a reconstructed track: {min, max}, with convention: min <= TPCNcls <= max"};
  Configurable<std::vector<float>> cfDCAxyCutRange{"cfDCAxyCutRange", {-3.2, 3.2}, "range of cut on signed transverse distance of closest approach between the reconstructed track and the primary interaction vertex in x-y plane: {min, max}, with convention: min <= dcaXY <= max"};
  Configurable<std::vector<float>> cfDCAzCutRange{"cfDCAzCutRange", {-2.4, 2.4}, "range of cut on signed transverse distance of closest approach between the reconstructed track and the primary interaction vertex in z direction: {min, max}, with convention: min <= dcaZ <= max"};

  // *) Misc
  Configurable<float> cfEtaGap{"cfEtaGap", 1.0, "eta gap for two particle correlation"};
  Configurable<float> cfSigmaInel{"cfSigmaInel", 7.71, "inelastic cross section in mb"};
  Configurable<bool> cfQualityAssuranceSwitch{"cfQualityAssuranceSwitch", false, "quality assurance switch"};
  Configurable<bool> cfRunMessageSwitch{"cfRunMessageSwitch", false, "run message switch"};
  Configurable<bool> cfInternalValidationSwitch{"cfInternalValidationSwitch", false, "internal validation switch"};
  Configurable<std::vector<float>> cfIVHarmonics{"cfIVHarmonics", {0.06, 0.07, 0.08}, "harmonics for internal validation: {v_2, v_3, v_4, ...}"};

  // *) Define and initialize all data members to be called in the main process* functions:
  // **) Task configuration:
  struct TaskConfiguration {
    float fSigmaInel = 7.71;
    float fEtaGap = 1.0;

    std::vector<std::string> fEventCutSwitch{};
    std::vector<float> fVertexZCutRange = {-10., 10.};
    std::vector<float> fCentCutRange = {0., 80.};
    std::vector<std::string> fParticleCutSwitch{};
    std::vector<float> fPtCutRange = {0.2, 5.};
    std::vector<float> fEtaCutRange = {-0.8, 0.8};
    std::vector<int> fTPCNumClusterCutRange = {70, 160};
    std::vector<float> fDCAxyCutRange = {-3.2, 3.2};
    std::vector<float> fDCAzCutRange = {-2.4, 2.4};

    std::vector<float> fIVHarmonics{};

    std::string fFileWithWeights = "/alice-ccdb.cern.ch/Users/m/mei/thesis-";

    int fCentralityEstimator = 0;
    int fMultiplicityTables = 0;
    int fNContribCutInf = 2;

    bool fDryRun = false; // book all histos and run without filling and calculating anything
    bool fExternalFileSwitch = false;
    bool fMasterCutSwitch = true;
    bool fQualityAssuranceSwitch = false;
    bool fRunMessageSwitch = false;
    bool fInternalValidationSwitch = false;

    std::array<bool, eProcess_N> fProcess{false}; // Set what to process. See enum EProcess for full description. Set via implicit variables within a PROCESS_SWITCH clause.
  } tc;

  struct ParticleHistograms {
    TList* fParticleHistList = nullptr;
    std::array<std::array<std::array<TH1F*, eBeforeAfterCuts_N>, 2>, eParticleHistograms_N> fParticleHist{};
  } pc;

  struct EventHistograms {
    TList* fEventHistList = nullptr;
    std::array<std::array<std::array<TH1F*, eBeforeAfterCuts_N>, 2>, eEventHistograms_N> fEventHist{}; //! [ type - see enum EEventHistograms ][reco,sim][before, after event cuts]
  } ec;

  struct ExternalHistograms {
    TList* fExternalHistogramsList = nullptr;
    std::array<TH1F*, eWeightsHistograms_N> fWeights{};
  } ex;

  struct Observables {
    TList* fObservablesList = nullptr;
    std::array<std::array<TProfile*, 2>, eObservables_N> fProfCorr{}; //! [type][reco,sim]
  } obs;

  struct QualityAssurance {
    TList* fQualityAssuranceList = nullptr;
    std::array<TH2F*, eBeforeAfterCuts_N> fHistCentralityRecSim{};
    std::array<TH2F*, eBeforeAfterCuts_N> fHistMultNContrib{};
  } qa;

  struct MiscHistograms {
    TList* fMiscHistList = nullptr;
    TProfile* fMiscHistRunNumber = nullptr;
    std::array<TH1D*, 2> fMiscEventCutsCounter{};    //! [reco,sim]
    std::array<TH1D*, 2> fMiscParticleCutsCounter{}; //! [reco,sim]
  } misc;

  struct InternalValidation {
    TRandom3* randGen = new TRandom3(1);
    TF1* fIntValFun = nullptr;
    TList* fIntValList = nullptr;
    std::array<TProfile*, eIntValHarmonics_N> fProfIntVal{};
  } intval;

  struct RunTimeInfo {
    bool isFirstCollision = true; // this is used to ensure that weights hist are booked once and only once, otherwise, a crash
    int collisionCounter = 0;
    int thisRunNumber = 0;
    static constexpr int RunMessagePeriod = 100;
  } rti;

  // *) functions and templates
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
  void updateQVectorsTable(std::vector<std::vector<TComplex>>& QVectorsTable, T1 const& phi, T1 const& weight)
  {
    for (std::size_t h = 0; h < QVectorsTable.size(); ++h) {
      for (std::size_t p = 0; p < QVectorsTable[h].size(); ++p) {
        const auto wp = std::pow(weight, static_cast<int>(p));
        QVectorsTable[h][p] += TComplex(wp * std::cos(h * phi), wp * std::sin(h * phi));
      }
    }
  }

  template <typename T1>
  void updateVecForGap(std::vector<int> const& harmonic, std::vector<TComplex>& qVecAB, std::vector<float>& weightAB, T1 const& phi, T1 const& eta, T1 const& weight)
  {
    if (eta < -tc.fEtaGap / 2.f) {
      qVecAB[0] += TComplex(weight * std::cos(harmonic[0] * phi), weight * std::sin(harmonic[0] * phi));
      weightAB[0] += weight;
    }
    if (eta > tc.fEtaGap / 2.f) {
      qVecAB[1] += TComplex(weight * std::cos(harmonic[1] * phi), weight * std::sin(harmonic[1] * phi));
      weightAB[1] += weight;
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

  std::vector<TComplex> recursion(int m, const std::vector<std::vector<TComplex>>& qVector, std::vector<int>& harmonic, int mult = 1, int skip = 0)
  {
    // Calculate multi-particle correlators by using recursion (an improved faster version) originally developed by
    // Kristjan Gulbrandsen (gulbrand@nbi.dk).

    auto q = [&](int n, int p) -> TComplex {
      if (n >= 0) {
        return qVector[n][p];
      }
      return TComplex::Conjugate(qVector[-n][p]);
    };
    std::vector<TComplex> c = {q(harmonic[m - 1], mult), q(0, mult)};
    if ((m - 1) == 0) {
      return c;
    }
    std::vector<TComplex> temp = recursion(m - 1, qVector, harmonic);
    c[0] *= temp[0];
    c[1] *= temp[1];
    if ((m - 1) == skip) {
      return c;
    }
    int counter1 = 0;
    int hhold = harmonic[counter1];
    harmonic[counter1] = harmonic[m - 2];
    harmonic[m - 2] = hhold + harmonic[m - 1];
    std::vector<TComplex> c2 = recursion(m - 1, qVector, harmonic, mult + 1, m - 2);
    int counter2 = m - 3;
    while (counter2 >= skip) {
      harmonic[m - 2] = harmonic[counter1];
      harmonic[counter1] = hhold;
      ++counter1;
      hhold = harmonic[counter1];
      harmonic[counter1] = harmonic[m - 2];
      harmonic[m - 2] = hhold + harmonic[m - 1];
      temp = recursion(m - 1, qVector, harmonic, mult + 1, counter2);
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
    // Get TObject pointer from TList, even if it's in some nested TList. Foreseen to be used to fetch histograms or profiles from files directly.
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

  template <ERecSim rs>
  void getHistogramWithWeights(const char* filePath, const char* runNumber)
  {
    // *) Return value:
    std::array<TH1F*, eWeightsHistograms_N> histW{};
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
    } else if (bFileIsInCCDB) {
      // File you want to access is in your home dir in CCDB: Remember that here I do not access the file; instead, I directly access the object in that file.
      ccdb->setURL("https://alice-ccdb.cern.ch");
      baseList = dynamic_cast<TList*>(ccdb->get<TList>(TString(filePath).ReplaceAll("/alice-ccdb.cern.ch/", "").Data()));
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
    }

    // Here comes the common code for all three cases, where from "listWithRuns" you fetch the desired histogram with efficiency corrections:
    if (!baseList) {
      LOGF(fatal, "\033[1;31m%s at line %d for case %s\033[0m", __FUNCTION__, __LINE__, filePath);
    }
    listWithRuns = dynamic_cast<TList*>(getObjectFromList(baseList, runNumber));
    if (!listWithRuns) {
      TString runNumberWithLeadingZeroes = "000";
      runNumberWithLeadingZeroes += runNumber; // another try, with "000" prepended to run number
      listWithRuns = dynamic_cast<TList*>(getObjectFromList(baseList, runNumberWithLeadingZeroes.Data()));
      if (!listWithRuns) {
        LOGF(warning, "\033[1;31m%s at line %d : this crash can happen if in the output file there is no list with weights for the current run number = %s\033[0m", __FUNCTION__, __LINE__, runNumber);
        return;
      }
    }

    for (int i = 0; i < eWeightsHistograms_N; ++i) {
      if constexpr (rs == eRec) {
        if (i == eWPhiSim) {
          continue;
        }
      }
      histW[i] = dynamic_cast<TH1F*>(listWithRuns->FindObject(Form("[%s]", WeightsNames[i])));
      if (!histW[i]) {
        LOGF(warning, "%s: histogram '%s' not found in run list with run number %s", __FUNCTION__, WeightsNames[i], runNumber);
        continue;
      }
      histW[i]->SetDirectory(nullptr);
      ex.fWeights[i] = dynamic_cast<TH1F*>(histW[i]->Clone());
      if (!ex.fWeights[i]) {
        LOGF(warning, "%s: failed to clone histogram '%s' with run number %s", __FUNCTION__, WeightsNames[i], runNumber);
        continue;
      }
      ex.fWeights[i]->SetStats(false);
      ex.fWeights[i]->SetTitle(filePath);
      ex.fWeights[i]->SetDirectory(nullptr);
      ex.fExternalHistogramsList->Add(ex.fWeights[i]);
    }
  } // end of void getHistogramWithWeights(const char* filePath, const char* runNumber)

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
    if constexpr ((rs == eSim || rs == eRecAndSim) && rm == eMC) {
      if (!collision.has_mcCollision()) {
        return false;
      }
    }

    for (int i = 0; i < EEventCuts_N; ++i) {
      bool doNotPassCut = false;
      switch (i) {
        default:
          doNotPassCut = false;
          break;
        case eSel8Cut:
          doNotPassCut = tc.fEventCutSwitch[eSel8Cut].starts_with("1") && !collision.sel8();
          break;
        case eNoCollInTimeRangeStandard:
          doNotPassCut = tc.fEventCutSwitch[eNoCollInTimeRangeStandard].starts_with("1") && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard);
          break;
        case eNoCollInRofStandard:
          doNotPassCut = tc.fEventCutSwitch[eNoCollInRofStandard].starts_with("1") && !collision.selection_bit(o2::aod::evsel::kNoCollInRofStandard);
          break;
        case eNoSameBunchPileUp:
          doNotPassCut = tc.fEventCutSwitch[eNoSameBunchPileUp].starts_with("1") && !collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup);
          break;
        case eIsVertexITSTPC:
          doNotPassCut = tc.fEventCutSwitch[eIsVertexITSTPC].starts_with("1") && !collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC);
          break;
        case eIsGoodITSLayersAll:
          doNotPassCut = tc.fEventCutSwitch[eIsGoodITSLayersAll].starts_with("1") && !collision.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll);
          break;
        case eIsGoodZvtxFT0vsPV:
          doNotPassCut = tc.fEventCutSwitch[eIsGoodZvtxFT0vsPV].starts_with("1") && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV);
          break;
        case eNoHighMultCollInPrevRof:
          doNotPassCut = tc.fEventCutSwitch[eNoHighMultCollInPrevRof].starts_with("1") && !collision.selection_bit(o2::aod::evsel::kNoHighMultCollInPrevRof);
          break;
        case eNContribCut:
          doNotPassCut = tc.fEventCutSwitch[eNContribCut].starts_with("1") && (collision.numContrib() < tc.fNContribCutInf);
          break;
        case eVertexZCut:
          if constexpr ((rs == eRec || rs == eRecAndSim) && rm == eReal) {
            doNotPassCut = tc.fEventCutSwitch[eVertexZCut].starts_with("1") && (collision.posZ() > tc.fVertexZCutRange[1] || collision.posZ() < tc.fVertexZCutRange[0]);
          }
          if constexpr ((rs == eSim || rs == eRecAndSim) && rm == eMC) {
            doNotPassCut = tc.fEventCutSwitch[eVertexZCut].starts_with("1") && (collision.mcCollision().posZ() > tc.fVertexZCutRange[1] || collision.mcCollision().posZ() < tc.fVertexZCutRange[0]);
          }
          break;
        case eCentCut:
          if constexpr ((rs == eRec || rs == eRecAndSim) && rm == eReal) {
            auto thisCent = chooseCent(collision, tc.fCentralityEstimator);
            doNotPassCut = tc.fEventCutSwitch[eCentCut].starts_with("1") && (thisCent > tc.fCentCutRange[1] || thisCent < tc.fCentCutRange[0]);
          }
          if constexpr ((rs == eSim || rs == eRecAndSim) && rm == eMC) {
            auto impactParameter = collision.mcCollision().impactParameter();
            auto centralityMC = math::PI * impactParameter * impactParameter / tc.fSigmaInel;
            doNotPassCut = tc.fEventCutSwitch[eCentCut].starts_with("1") && (centralityMC > tc.fCentCutRange[1] || centralityMC < tc.fCentCutRange[0]);
          }
          break;
      } // end of switch (i) {

      if (doNotPassCut) {
        if constexpr (rm == eReal) {
          misc.fMiscEventCutsCounter[eRec]->Fill(i + 0.5, 1);
        }
        if constexpr (rm == eMC) {
          misc.fMiscEventCutsCounter[eSim]->Fill(i + 0.5, 1);
        }
        return false;
      }
    } // end of for (int i = 0; i < EEventCuts_N; ++i) {

    return true;
  } // end of bool eventCuts(T1 const& collision)

  template <ERecSim rs, ERealMC rm, EBeforeAfterCuts cuts, typename T1, typename T2>
  void eventHistFill(T1 const& collision, T2 const& tracks)
  {
    if constexpr ((rs == eRec || rs == eRecAndSim) && rm == eReal) { // Fill reconstructed-level event histograms
      auto thisCent = chooseCent(collision, tc.fCentralityEstimator);
      auto thisRefMult = chooseMult(collision, tc.fMultiplicityTables);
      int multiplicityRec = static_cast<int>(tracks.size());
      if constexpr (cuts == eBefore) {
        ec.fEventHist[eHistMult][eRec][eBefore]->Fill(multiplicityRec);
        ec.fEventHist[eHistCentrality][eRec][eBefore]->Fill(thisCent);
        ec.fEventHist[eHistRefMult][eRec][eBefore]->Fill(thisRefMult);
        ec.fEventHist[eHistVertexX][eRec][eBefore]->Fill(collision.posX());
        ec.fEventHist[eHistVertexY][eRec][eBefore]->Fill(collision.posY());
        ec.fEventHist[eHistVertexZ][eRec][eBefore]->Fill(collision.posZ());
        ec.fEventHist[eHistNumContrib][eRec][eBefore]->Fill(collision.numContrib());
      }
      if constexpr (cuts == eAfter) {
        ec.fEventHist[eHistMult][eRec][eAfter]->Fill(multiplicityRec);
        ec.fEventHist[eHistCentrality][eRec][eAfter]->Fill(thisCent);
        ec.fEventHist[eHistRefMult][eRec][eAfter]->Fill(thisRefMult);
        ec.fEventHist[eHistVertexX][eRec][eAfter]->Fill(collision.posX());
        ec.fEventHist[eHistVertexY][eRec][eAfter]->Fill(collision.posY());
        ec.fEventHist[eHistVertexZ][eRec][eAfter]->Fill(collision.posZ());
        ec.fEventHist[eHistNumContrib][eRec][eAfter]->Fill(collision.numContrib());
      }
    }

    // Fill MC simulated-level event histograms if both reconstructed and simulated data are processed
    if constexpr ((rs == eRecAndSim || rs == eSim) && rm == eMC) {
      if (!collision.has_mcCollision()) {
        return;
      }
      auto thisMCCollision = collision.mcCollision();
      int multiplicitySim = static_cast<int>(tracks.size());
      auto impactParameter = thisMCCollision.impactParameter();
      auto centralityMC = math::PI * impactParameter * impactParameter / tc.fSigmaInel; // centrality for sim derived from impact parameter
      if constexpr (cuts == eBefore) {
        ec.fEventHist[eHistMult][eSim][eBefore]->Fill(multiplicitySim);
        ec.fEventHist[eHistCentrality][eSim][eBefore]->Fill(centralityMC);
        ec.fEventHist[eHistImpactParameter][eSim][eBefore]->Fill(impactParameter);
        ec.fEventHist[eHistVertexX][eSim][eBefore]->Fill(thisMCCollision.posX());
        ec.fEventHist[eHistVertexY][eSim][eBefore]->Fill(thisMCCollision.posY());
        ec.fEventHist[eHistVertexZ][eSim][eBefore]->Fill(thisMCCollision.posZ());
      }
      if constexpr (cuts == eAfter) {
        ec.fEventHist[eHistMult][eSim][eAfter]->Fill(multiplicitySim);
        ec.fEventHist[eHistCentrality][eSim][eAfter]->Fill(centralityMC);
        ec.fEventHist[eHistImpactParameter][eSim][eAfter]->Fill(impactParameter);
        ec.fEventHist[eHistVertexX][eSim][eAfter]->Fill(thisMCCollision.posX());
        ec.fEventHist[eHistVertexY][eSim][eAfter]->Fill(thisMCCollision.posY());
        ec.fEventHist[eHistVertexZ][eSim][eAfter]->Fill(thisMCCollision.posZ());
      }
    } // end of if constexpr (rs == eRecAndSim && rm == eMC) {
  } // end of void eventHistFill(T1 const& collision, T2 const& tracks)

  template <ERecSim rs, ERealMC rm, typename T>
  bool particleCuts(T const& track)
  {
    for (int i = 0; i < EParticleCuts_N; ++i) {
      bool doNotPassCut = false;
      if constexpr ((rs == eRec || rs == eRecAndSim) && rm == eReal) {
        switch (i) {
          default:
            doNotPassCut = false;
            break;
          case ePtCut:
            doNotPassCut = tc.fParticleCutSwitch[ePtCut].starts_with("1") && (track.pt() < tc.fPtCutRange[0] || track.pt() > tc.fPtCutRange[1]);
            break;
          case eEtaCut:
            doNotPassCut = tc.fParticleCutSwitch[eEtaCut].starts_with("1") && (track.eta() < tc.fEtaCutRange[0] || track.eta() > tc.fEtaCutRange[1]);
            break;
          case eChargeCut:
            doNotPassCut = tc.fParticleCutSwitch[eChargeCut].starts_with("1") && (track.sign() != 1 && track.sign() != -1);
            break;
          case eTPCNumClusterCut:
            doNotPassCut = tc.fParticleCutSwitch[eTPCNumClusterCut].starts_with("1") && (track.tpcNClsFound() < tc.fTPCNumClusterCutRange[0] || track.tpcNClsFound() > tc.fTPCNumClusterCutRange[1]);
            break;
          case eDCAxyCut:
            doNotPassCut = tc.fParticleCutSwitch[eDCAxyCut].starts_with("1") && (track.dcaXY() < tc.fDCAxyCutRange[0] || track.dcaXY() > tc.fDCAxyCutRange[1]);
            break;
          case eDCAzCut:
            doNotPassCut = tc.fParticleCutSwitch[eDCAzCut].starts_with("1") && (track.dcaZ() < tc.fDCAzCutRange[0] || track.dcaZ() > tc.fDCAzCutRange[1]);
            break;
        }
      } // end of if constexpr ((rs == eRec || rs == eRecAndSim) && rm == eReal) {

      if constexpr ((rs == eSim || rs == eRecAndSim) && rm == eMC) {
        if (!track.has_mcParticle()) {
          return false;
        }
        auto thisMCParticle = track.mcParticle();
        switch (i) {
          default:
            doNotPassCut = false;
            break;
          case ePtCut:
            doNotPassCut = tc.fParticleCutSwitch[ePtCut].starts_with("1") && (thisMCParticle.pt() < tc.fPtCutRange[0] || thisMCParticle.pt() > tc.fPtCutRange[1]);
            break;
          case eEtaCut:
            doNotPassCut = tc.fParticleCutSwitch[eEtaCut].starts_with("1") && (thisMCParticle.eta() < tc.fEtaCutRange[0] || thisMCParticle.eta() > tc.fEtaCutRange[1]);
            break;
          case eChargeCut:
            if (tc.fParticleCutSwitch[eChargeCut].starts_with("1")) {
              TParticlePDG* particlePDG = pdg->GetParticle(thisMCParticle.pdgCode());
              if (particlePDG) {
                const int chargeMCUnit = 3;
                auto chargeMC = particlePDG->Charge() / chargeMCUnit;
                if (chargeMC != 1 && chargeMC != -1) {
                  doNotPassCut = true;
                }
              } else {
                doNotPassCut = true;
              }
            } // end of if (tc.fParticleCutSwitch[eChargeCut]) {
            break;
          case eTPCNumClusterCut:
            doNotPassCut = tc.fParticleCutSwitch[eTPCNumClusterCut].starts_with("1") && (track.tpcNClsFound() < tc.fTPCNumClusterCutRange[0] || track.tpcNClsFound() > tc.fTPCNumClusterCutRange[1]);
            break;
          case eDCAxyCut:
            doNotPassCut = tc.fParticleCutSwitch[eDCAxyCut].starts_with("1") && (track.dcaXY() < tc.fDCAxyCutRange[0] || track.dcaXY() > tc.fDCAxyCutRange[1]);
            break;
          case eDCAzCut:
            doNotPassCut = tc.fParticleCutSwitch[eDCAzCut].starts_with("1") && (track.dcaZ() < tc.fDCAzCutRange[0] || track.dcaZ() > tc.fDCAzCutRange[1]);
            break;
        }
      } // end of if constexpr ((rs == eSim || rs == eRecAndSim) && rm == eMC) {

      if (doNotPassCut) {
        if constexpr (rm == eReal) {
          misc.fMiscParticleCutsCounter[eRec]->Fill(i + 0.5, 1);
        }
        if constexpr (rm == eMC) {
          misc.fMiscParticleCutsCounter[eSim]->Fill(i + 0.5, 1);
        }
        return false;
      }
    } // end of for (int i = 0; i < EParticleCuts_N; ++i) {

    return true;
  } // end of bool particleCuts(T const& track)

  template <ERecSim rs, ERealMC rm, EBeforeAfterCuts cuts, typename T1>
  void particleHistFill(T1 const& track)
  {
    if constexpr ((rs == eRec || rs == eRecAndSim) && rm == eReal) {
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

    if constexpr ((rs == eSim || rs == eRecAndSim) && rm == eMC) {
      if (!track.has_mcParticle()) {
        return;
      }
      auto thisMCParticle = track.mcParticle();
      TParticlePDG* particlePDG = pdg->GetParticle(thisMCParticle.pdgCode());
      float chargeMC = 0.f;
      const int chargeMCUnit = 3;
      if (particlePDG) {
        chargeMC = particlePDG->Charge();
        chargeMC /= chargeMCUnit;
      }
      if constexpr (cuts == eBefore) {
        pc.fParticleHist[eHistPt][eSim][eBefore]->Fill(thisMCParticle.pt());
        pc.fParticleHist[eHistPhi][eSim][eBefore]->Fill(thisMCParticle.phi());
        pc.fParticleHist[eHistEta][eSim][eBefore]->Fill(thisMCParticle.eta());
        if (particlePDG) {
          pc.fParticleHist[eHistCharge][eSim][eBefore]->Fill(chargeMC);
        }
      }
      if constexpr (cuts == eAfter) {
        pc.fParticleHist[eHistPt][eSim][eAfter]->Fill(thisMCParticle.pt());
        pc.fParticleHist[eHistPhi][eSim][eAfter]->Fill(thisMCParticle.phi());
        pc.fParticleHist[eHistEta][eSim][eAfter]->Fill(thisMCParticle.eta());
        if (particlePDG) {
          pc.fParticleHist[eHistCharge][eSim][eAfter]->Fill(chargeMC);
        }
      }
    } // end of if constexpr (rs == eRecAndSim && rm == eMC) {
  } // end of void particleHistFill(T1 const& track)

  template <ERecSim rs, ERealMC rm, EBeforeAfterCuts cuts, typename T1>
  void qaFill(T1 const& collision)
  {
    auto thisCent = chooseCent(collision, tc.fCentralityEstimator);
    auto thisRefMult = chooseMult(collision, tc.fMultiplicityTables);
    if constexpr (rm == eReal) {
      if constexpr (cuts == eBefore) {
        qa.fHistMultNContrib[eBefore]->Fill(thisRefMult, collision.numContrib());
      }
      if constexpr (cuts == eAfter) {
        qa.fHistMultNContrib[eAfter]->Fill(thisRefMult, collision.numContrib());
      }
    }

    if constexpr (rs == eRecAndSim && rm == eMC) {
      if (!collision.has_mcCollision()) {
        return;
      }
      auto thisMCCollision = collision.mcCollision();
      auto impactParameter = thisMCCollision.impactParameter();
      auto centralityMC = math::PI * impactParameter * impactParameter / tc.fSigmaInel;
      if constexpr (cuts == eBefore) {
        qa.fHistCentralityRecSim[eBefore]->Fill(thisCent, centralityMC);
      }
      if constexpr (cuts == eAfter) {
        qa.fHistCentralityRecSim[eAfter]->Fill(thisCent, centralityMC);
      }
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
    auto thisCollCentReal = chooseCent(collision, tc.fCentralityEstimator);

    if (rti.isFirstCollision) {
      rti.thisRunNumber = collision.bc().runNumber();
      misc.fMiscHistRunNumber->SetTitle(Form("%d", rti.thisRunNumber)); // Get run number
      misc.fMiscHistRunNumber->Fill(0.5, rti.thisRunNumber);
      if (tc.fExternalFileSwitch) {
        getHistogramWithWeights<rs>(tc.fFileWithWeights.c_str(), Form("%d", rti.thisRunNumber));
      }
    }

    if (tc.fRunMessageSwitch) {
      if (rti.collisionCounter % rti.RunMessagePeriod == 0) {
        LOGF(info, "Successfully running, run number is %d", rti.thisRunNumber);
        rti.collisionCounter = 0;
      }
      ++rti.collisionCounter;
    }

    // Fill Quality Assurance before cuts
    if (tc.fQualityAssuranceSwitch) {
      qaFill<rs, eReal, eBefore>(collision);
      qaFill<rs, eMC, eBefore>(collision);
    }

    // Fill Event Hist before cuts
    eventHistFill<rs, eReal, eBefore>(collision, tracks);
    eventHistFill<rs, eMC, eBefore>(collision, tracks);

    bool passEventCutsReal = eventCuts<rs, eReal>(collision);
    // bool passEventCutsMC = eventCuts<rs, eMC>(collision); // this is currently useless

    // Fill event hist and qa hist after cuts
    if (tc.fMasterCutSwitch) {
      if (passEventCutsReal) {
        eventHistFill<rs, eReal, eAfter>(collision, tracks);
        eventHistFill<rs, eMC, eAfter>(collision, tracks);
        if (tc.fQualityAssuranceSwitch) {
          qaFill<rs, eReal, eAfter>(collision);
          qaFill<rs, eMC, eAfter>(collision);
        }
      }
    }

    std::vector<int> harmonic22 = {-2, 2};
    std::vector<int> harmonic24 = {-2, -2, 2, 2};

    auto qVecTableReal = initQVectorsTable(4, harmonic24);
    auto qVecTableMC = initQVectorsTable(4, harmonic24);
    auto qVecTableIV = initQVectorsTable(4, harmonic24);

    std::vector<TComplex> qVecABReal(2, TComplex(0., 0.));
    std::vector<TComplex> qVecABMC(2, TComplex(0., 0.));
    std::vector<float> weightABReal(2, 0.f);
    std::vector<float> weightABMC(2, 0.f);

    if (tc.fInternalValidationSwitch) {
      double psi = intval.randGen->Uniform(0, math::TwoPI);
      intval.fIntValFun->SetParameter(0, psi);
    }

    // Main loop over particles:
    for (auto const& track : tracks) {
      float thisPhiRec = track.phi();
      float thisPt = track.pt();

      if (tc.fInternalValidationSwitch && passEventCutsReal) {
        const float thisPhiIV = intval.fIntValFun->GetRandom();
        updateQVectorsTable(qVecTableIV, thisPhiIV, 1.f);
      }

      std::array<float, 3> thisPhiAndPt = {thisPhiRec, 0.f, thisPt};         // {thisPhiRec, thisPhiSim, thisPt}
      std::array<float, eWeightsHistograms_N> thisWeights = {1.f, 1.f, 1.f}; // {wPhiRec, wPhiSim, wPt}

      particleHistFill<rs, eReal, eBefore>(track);
      if (tc.fMasterCutSwitch && passEventCutsReal && particleCuts<rs, eReal>(track)) {
        particleHistFill<rs, eReal, eAfter>(track);
        if (tc.fExternalFileSwitch) {
          for (int k = 0; k < eWeightsHistograms_N; ++k) {
            if (k != eWPhiSim && ex.fWeights[k]) {
              thisWeights[k] = ex.fWeights[k]->GetBinContent(ex.fWeights[k]->FindBin(thisPhiAndPt[k]));
            }
          }
        }

        float thisWeightReal = thisWeights[eWPhiRec] * thisWeights[eWPt];
        updateQVectorsTable(qVecTableReal, thisPhiRec, thisWeightReal);
        auto thisEtaRec = track.eta();
        updateVecForGap(harmonic22, qVecABReal, weightABReal, thisPhiRec, thisEtaRec, thisWeightReal);
      }

      if constexpr (rs == eRecAndSim) {
        particleHistFill<rs, eMC, eBefore>(track);
        if (tc.fMasterCutSwitch && passEventCutsReal && particleCuts<rs, eMC>(track)) {
          particleHistFill<rs, eMC, eAfter>(track);

          auto thisPhiSim = track.mcParticle().phi();
          auto thisEtaSim = track.mcParticle().eta();
          thisPhiAndPt[eWPhiSim] = thisPhiSim;
          if (tc.fExternalFileSwitch && ex.fWeights[eWPhiSim]) {
            thisWeights[eWPhiSim] = ex.fWeights[eWPhiSim]->GetBinContent(ex.fWeights[eWPhiSim]->FindBin(thisPhiSim));
          }

          float thisWeightSim = thisWeights[eWPhiSim] * thisWeights[eWPt];
          updateQVectorsTable(qVecTableMC, thisPhiSim, thisWeightSim);
          updateVecForGap(harmonic22, qVecABMC, weightABMC, thisPhiSim, thisEtaSim, thisWeightSim);
        }
      } // end of if constexpr (rs == eRecAndSim) {
    } // end of for (auto const& track : tracks) {

    if (tc.fMasterCutSwitch) {
      std::vector<TComplex> result22(2, TComplex(0., 0.));
      std::vector<TComplex> result24(2, TComplex(0., 0.));
      std::vector<float> result22Gap(2, 0.f);
      if (passEventCutsReal) {
        result22 = two(qVecTableReal, harmonic22);
        if (noneZeroDenom(result22)) {
          obs.fProfCorr[eV22][eRec]->Fill(thisCollCentReal, (result22[0] / result22[1].Re()).Re());
        }
        result24 = recursion(4, qVecTableReal, harmonic24);
        if (noneZeroDenom(result24)) {
          obs.fProfCorr[eV24][eRec]->Fill(thisCollCentReal, (result24[0] / result24[1].Re()).Re());
        }

        result22Gap = {static_cast<float>((qVecABReal[0] * qVecABReal[1]).Re()), weightABReal[0] * weightABReal[1]};
        if (result22Gap[1] != 0) {
          obs.fProfCorr[eV22Gap][eRec]->Fill(thisCollCentReal, result22Gap[0] / result22Gap[1]);
        }

        if constexpr (rs == eSim || rs == eRecAndSim) {
          auto impactParam = collision.mcCollision().impactParameter();
          auto thisCollCentMC = math::PI * impactParam * impactParam / tc.fSigmaInel;

          result22 = two(qVecTableMC, harmonic22);
          if (noneZeroDenom(result22)) {
            obs.fProfCorr[eV22][eSim]->Fill(thisCollCentMC, (result22[0] / result22[1].Re()).Re());
          }
          result24 = recursion(4, qVecTableMC, harmonic24);
          if (noneZeroDenom(result24)) {
            obs.fProfCorr[eV24][eSim]->Fill(thisCollCentMC, (result24[0] / result24[1].Re()).Re());
          }

          result22Gap = {static_cast<float>((qVecABMC[0] * qVecABMC[1]).Re()), weightABMC[0] * weightABMC[1]};
          if (result22Gap[1] != 0) {
            obs.fProfCorr[eV22Gap][eSim]->Fill(thisCollCentMC, result22Gap[0] / result22Gap[1]);
          }
        } // end of if constexpr (rs == eRecAndSim) {

        if (tc.fInternalValidationSwitch) {
          std::vector<std::vector<int>> harmonicIV = {{-2, 2}, {-3, 3}, {-4, 4}};
          for (int iv = 0; iv < eIntValHarmonics_N; ++iv) {
            std::vector<TComplex> result = two(qVecTableIV, harmonicIV[iv]);
            if (noneZeroDenom(result)) {
              intval.fProfIntVal[iv]->Fill(thisCollCentReal, (result[0] / result[1].Re()).Re());
            }
          }
        } // end of if (tc.fInternalValidationSwitch) {
      } // end of if (passEventCutsReal) {
    } // end of if (tc.fMasterCutSwitch) {

    rti.isFirstCollision = false; // Now the first collision ends
  } // end of template <ERecSim rs, typename T1, typename T2> void steer(T1 const& collision, T2 const& tracks) {

  // *) Initialize and book all objects:
  void init(InitContext&)
  {
    // *) Set automatically what to process, from an implicit variable "doprocessSomEProcessName" within a PROCESS_SWITCH clause:
    tc.fProcess[eProcessRec] = doprocessRec;
    tc.fProcess[eProcessRecSim] = doprocessRecSim;
    tc.fProcess[eProcessSim] = doprocessSim;

    // *) Set runtime info
    rti.isFirstCollision = true;
    rti.collisionCounter = 0;
    rti.thisRunNumber = 0;

    // *) Set colors
    std::vector<int> fillColor(eBeforeAfterCuts_N);
    fillColor[eBefore] = kRed - 10;
    fillColor[eAfter] = kGreen - 10;

    std::vector<int> lineColor(eBeforeAfterCuts_N);
    lineColor[eBefore] = kRed;
    lineColor[eAfter] = kGreen;

    // *) Configure your task using configurables in the json file:
    tc.fDryRun = cfDryRun;
    tc.fCentralityEstimator = cfCentralityEstimator;
    tc.fMultiplicityTables = cfMultiplicityTables;

    tc.fExternalFileSwitch = cfExternalFileSwitch;
    tc.fFileWithWeights = cfFileWithWeights;

    tc.fMasterCutSwitch = cfMasterCutSwitch;

    tc.fEventCutSwitch = cfEventCutSwitch.value;
    tc.fVertexZCutRange = cfVertexZCutRange.value;
    tc.fCentCutRange = cfCentCutRange.value;
    tc.fNContribCutInf = cfNContribCutInf;

    tc.fParticleCutSwitch = cfParticleCutSwitch.value;
    tc.fPtCutRange = cfPtCutRange.value;
    tc.fEtaCutRange = cfEtaCutRange.value;
    tc.fTPCNumClusterCutRange = cfTPCNumClusterCutRange.value;
    tc.fDCAxyCutRange = cfDCAxyCutRange.value;
    tc.fDCAzCutRange = cfDCAzCutRange.value;

    tc.fEtaGap = cfEtaGap;
    tc.fSigmaInel = cfSigmaInel;
    tc.fQualityAssuranceSwitch = cfQualityAssuranceSwitch;
    tc.fRunMessageSwitch = cfRunMessageSwitch;
    tc.fInternalValidationSwitch = cfInternalValidationSwitch;
    tc.fIVHarmonics = cfIVHarmonics.value;

    // *) Book base list:
    auto* temp = new TList();
    temp->SetOwner(true);
    fBaseList.setObject(temp);

    // *) Book and nest all other TLists:
    // *) Book External Hist List
    if (tc.fExternalFileSwitch) {
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

    misc.fMiscHistRunNumber = new TProfile("[eRunNumber]", "to be set", 1, 0, 1);
    misc.fMiscHistRunNumber->SetStats(false);
    misc.fMiscHistRunNumber->GetXaxis()->SetBinLabel(1, "RunNumber");
    misc.fMiscHistList->Add(misc.fMiscHistRunNumber);

    if (doprocessRec || doprocessRecSim) {
      misc.fMiscEventCutsCounter[eRec] = new TH1D("[eEventCutsCounter][eRec]", "Event Level Cuts Counter for Reconstructed Data (Sequential)", EEventCuts_N, 0, EEventCuts_N);
      misc.fMiscEventCutsCounter[eRec]->SetStats(false);
      for (int i = 1; i <= EEventCuts_N; ++i) {
        misc.fMiscEventCutsCounter[eRec]->GetXaxis()->SetBinLabel(i, tc.fEventCutSwitch[i - 1].c_str());
      }
      misc.fMiscHistList->Add(misc.fMiscEventCutsCounter[eRec]);

      misc.fMiscParticleCutsCounter[eRec] = new TH1D("[eParticleCutsCounter][eRec]", "Particle Level Cuts Counter for Reconstructed Data (Sequential)", EParticleCuts_N, 0, EParticleCuts_N);
      misc.fMiscParticleCutsCounter[eRec]->SetStats(false);
      for (int i = 1; i <= EParticleCuts_N; ++i) {
        misc.fMiscParticleCutsCounter[eRec]->GetXaxis()->SetBinLabel(i, tc.fParticleCutSwitch[i - 1].c_str());
      }
      misc.fMiscHistList->Add(misc.fMiscParticleCutsCounter[eRec]);
    }

    if (doprocessSim || doprocessRecSim) {
      misc.fMiscEventCutsCounter[eSim] = new TH1D("[eEventCutsCounter][eSim]", "Event Level Cuts Counter for Simulated Data (Sequential)", EEventCuts_N, 0, EEventCuts_N);
      misc.fMiscEventCutsCounter[eSim]->SetStats(false);
      for (int i = 1; i <= EEventCuts_N; ++i) {
        misc.fMiscEventCutsCounter[eSim]->GetXaxis()->SetBinLabel(i, tc.fEventCutSwitch[i - 1].c_str());
      }
      misc.fMiscHistList->Add(misc.fMiscEventCutsCounter[eSim]);

      misc.fMiscParticleCutsCounter[eSim] = new TH1D("[eParticleCutsCounter][eSim]", "Particle Level Cuts Counter for Simulated Data (Sequential)", EParticleCuts_N, 0, EParticleCuts_N);
      misc.fMiscParticleCutsCounter[eSim]->SetStats(false);
      for (int i = 1; i <= EParticleCuts_N; ++i) {
        misc.fMiscParticleCutsCounter[eSim]->GetXaxis()->SetBinLabel(i, tc.fParticleCutSwitch[i - 1].c_str());
      }
      misc.fMiscHistList->Add(misc.fMiscParticleCutsCounter[eSim]);
    }

    // *) Book particle TLists:
    pc.fParticleHistList = new TList();
    pc.fParticleHistList->SetName("ParticleHistograms");
    pc.fParticleHistList->SetOwner(true);
    fBaseList->Add(pc.fParticleHistList); // any nested TList in the base TList appears as a subdir in the output ROOT file

    std::vector<std::vector<float>> lParticleHistBins(eParticleHistograms_N);

    lParticleHistBins[eHistPt] = cfPtBins.value;
    lParticleHistBins[eHistPhi] = cfPhiBins.value;
    lParticleHistBins[eHistEta] = cfEtaBins.value;
    lParticleHistBins[eHistCharge] = {5, -2.5, 2.5};

    if (doprocessRec || doprocessRecSim) {
      for (int i = 0; i < eParticleHistograms_N; ++i) {
        const auto& bins = lParticleHistBins[i];
        for (int ba = 0; ba < eBeforeAfterCuts_N; ++ba) {
          if (!tc.fMasterCutSwitch && ba == eAfter) {
            continue;
          }
          pc.fParticleHist[i][eRec][ba] = new TH1F(Form("[%s][eRec][%s]", ParticleHistNames[i], CutsNames[ba]), Form("%s Distribution for Reconstructed Particles %s Cuts", ParticleHistTitles[i], CutsTitles[ba]), static_cast<int>(bins[0]), bins[1], bins[2]);
          pc.fParticleHist[i][eRec][ba]->GetXaxis()->SetTitle(Form("%s", ParticleHistTitles[i]));
          pc.fParticleHist[i][eRec][ba]->SetColors(lineColor[ba], -1, fillColor[ba]);
          pc.fParticleHistList->Add(pc.fParticleHist[i][eRec][ba]);
        } // end of for (int ba = 0; ba < eBeforeAfterCuts_N; ++ba) {
      }
    } // end of if (doprocessRec || doprocessRecSim) {

    if (doprocessSim || doprocessRecSim) {
      for (int i = 0; i < eParticleHistograms_N; ++i) {
        const auto& bins = lParticleHistBins[i];
        for (int ba = 0; ba < eBeforeAfterCuts_N; ++ba) {
          if (!tc.fMasterCutSwitch && ba == eAfter) {
            continue;
          }
          pc.fParticleHist[i][eSim][ba] = new TH1F(Form("[%s][eSim][%s]", ParticleHistNames[i], CutsNames[ba]), Form("%s Distribution for Simulated Particles %s Cuts", ParticleHistTitles[i], CutsTitles[ba]), static_cast<int>(bins[0]), bins[1], bins[2]);
          pc.fParticleHist[i][eSim][ba]->GetXaxis()->SetTitle(Form("%s", ParticleHistTitles[i]));
          pc.fParticleHist[i][eSim][ba]->SetColors(lineColor[ba], -1, fillColor[ba]);
          pc.fParticleHistList->Add(pc.fParticleHist[i][eSim][ba]);
        } // end of for (int ba = 0; ba < eBeforeAfterCuts_N; ++ba) {
      }
    } // end of if (doprocessSim || doprocessRecSim) {

    // *) Book event-level histograms
    ec.fEventHistList = new TList();
    ec.fEventHistList->SetName("EventHistograms");
    ec.fEventHistList->SetOwner(true);
    fBaseList->Add(ec.fEventHistList);

    std::vector<std::vector<float>> lEventHistBins(eEventHistograms_N);

    lEventHistBins[eHistCentrality] = cfCentBins.value;
    lEventHistBins[eHistMult] = cfMultBins.value;
    lEventHistBins[eHistRefMult] = cfRefMultBins.value;
    lEventHistBins[eHistVertexX] = cfVxBins.value;
    lEventHistBins[eHistVertexY] = cfVyBins.value;
    lEventHistBins[eHistVertexZ] = cfVzBins.value;
    lEventHistBins[eHistImpactParameter] = cfIpBins.value;
    lEventHistBins[eHistNumContrib] = cfNumContribBins.value;

    if (doprocessRec || doprocessRecSim) {
      for (int i = 0; i < eEventHistograms_N; ++i) {
        if (i == eHistImpactParameter) { // no impact parameter for Rec
          continue;
        }
        const auto& bins = lEventHistBins[i];
        for (int ba = 0; ba < eBeforeAfterCuts_N; ++ba) {
          if (!tc.fMasterCutSwitch && ba == eAfter) {
            continue;
          }
          ec.fEventHist[i][eRec][ba] = new TH1F(Form("[%s][eRec][%s]", EventHistNames[i], CutsNames[ba]), Form("%s (Reconstructed) %s Cuts", EventHistTitles[i], CutsTitles[ba]), static_cast<int>(bins[0]), bins[1], bins[2]);
          ec.fEventHist[i][eRec][ba]->GetXaxis()->SetTitle(Form("%s", EventHistTitles[i]));
          ec.fEventHist[i][eRec][ba]->SetColors(lineColor[ba], -1, fillColor[ba]);
          ec.fEventHistList->Add(ec.fEventHist[i][eRec][ba]);
        } // end of for (int ba = 0; ba < eBeforeAfterCuts_N; ++ba) {
      } // end of for (int i = 0; i < eEventHistograms_N; ++i) {
      ec.fEventHist[eHistCentrality][eRec][eBefore]->GetXaxis()->SetTitle(Form("%s (%s)", EventHistTitles[eHistCentrality], CentralityEstimatorNames[cfCentralityEstimator]));
      ec.fEventHist[eHistCentrality][eRec][eAfter]->GetXaxis()->SetTitle(Form("%s (%s)", EventHistTitles[eHistCentrality], CentralityEstimatorNames[cfCentralityEstimator]));
      ec.fEventHist[eHistRefMult][eRec][eBefore]->GetXaxis()->SetTitle(Form("%s (%s)", EventHistTitles[eHistRefMult], MultiplicityTablesNames[cfMultiplicityTables]));
      ec.fEventHist[eHistRefMult][eRec][eAfter]->GetXaxis()->SetTitle(Form("%s (%s)", EventHistTitles[eHistRefMult], MultiplicityTablesNames[cfMultiplicityTables]));
    } // end of if (doprocessRec || doprocessRecSim) {

    if (doprocessSim || doprocessRecSim) {
      for (int i = 0; i < eEventHistograms_N; ++i) {
        if (i == eHistRefMult || i == eHistNumContrib) { // no RefMult for Sim
          continue;
        }
        const auto& bins = lEventHistBins[i];
        for (int ba = 0; ba < eBeforeAfterCuts_N; ++ba) {
          if (!tc.fMasterCutSwitch && ba == eAfter) {
            continue;
          }
          ec.fEventHist[i][eSim][ba] = new TH1F(Form("[%s][eSim][%s]", EventHistNames[i], CutsNames[ba]), Form("%s (Simulated) %s Cuts", EventHistTitles[i], CutsTitles[ba]), static_cast<int>(bins[0]), bins[1], bins[2]);
          ec.fEventHist[i][eSim][ba]->GetXaxis()->SetTitle(Form("%s", EventHistTitles[i]));
          ec.fEventHist[i][eSim][ba]->SetColors(lineColor[ba], -1, fillColor[ba]);
          ec.fEventHistList->Add(ec.fEventHist[i][eSim][ba]);
        } // end of for (int ba = 0; ba < eBeforeAfterCuts_N; ++ba) {
      } // end of for (int i = 0; i < eEventHistograms_N; ++i) {
      ec.fEventHist[eHistCentrality][eSim][eBefore]->GetXaxis()->SetTitle(Form("%s (calculated from I.P.)", EventHistTitles[eHistCentrality]));
      ec.fEventHist[eHistCentrality][eSim][eAfter]->GetXaxis()->SetTitle(Form("%s (calculated from I.P.)", EventHistTitles[eHistCentrality]));
    } // end of if (doprocessSim || doprocessRecSim) {

    // *) Book observales TLists:
    obs.fObservablesList = new TList();
    obs.fObservablesList->SetName("Observables");
    obs.fObservablesList->SetOwner(true);
    fBaseList->Add(obs.fObservablesList);

    const std::vector<float> defaultCentBoundaries = {0.f, 5.f, 10.f, 20.f, 30.f, 40.f, 50.f, 60.f, 70.f, 80.f};
    const int nDefaultCentBins = defaultCentBoundaries.size() - 1;

    if (doprocessRec || doprocessRecSim) {
      for (int i = 0; i < eObservables_N; ++i) {
        obs.fProfCorr[i][eRec] = new TProfile(Form("[%s][eRec]", ObservablesNames[i]), Form("%s for Reconstructed Data after Cuts", ObservablesTitles[i]), nDefaultCentBins, defaultCentBoundaries.data());
        obs.fProfCorr[i][eRec]->SetStats(false);
        obs.fProfCorr[i][eRec]->GetXaxis()->SetTitle(Form("Centrality (%s)", CentralityEstimatorNames[tc.fCentralityEstimator]));
        obs.fProfCorr[i][eRec]->GetYaxis()->SetTitle(ObservablesTitles[i]);
        obs.fProfCorr[i][eRec]->Sumw2();
        obs.fObservablesList->Add(obs.fProfCorr[i][eRec]);
      }
    }

    if (doprocessSim || doprocessRecSim) {
      for (int i = 0; i < eObservables_N; ++i) {
        obs.fProfCorr[i][eSim] = new TProfile(Form("[%s][eSim]", ObservablesNames[i]), Form("%s for Simulated Data after Cuts", ObservablesTitles[i]), nDefaultCentBins, defaultCentBoundaries.data());
        obs.fProfCorr[i][eSim]->SetStats(false);
        obs.fProfCorr[i][eSim]->GetXaxis()->SetTitle("Centrality (calculated from I.P.)");
        obs.fProfCorr[i][eSim]->GetYaxis()->SetTitle(ObservablesTitles[i]);
        obs.fProfCorr[i][eSim]->Sumw2();
        obs.fObservablesList->Add(obs.fProfCorr[i][eSim]);
      }
    }

    // *) Book and QA TLists:
    if (tc.fQualityAssuranceSwitch) {
      std::vector<float> lContrib = cfNumContribBins.value;
      const int nBinsContrib = static_cast<int>(lContrib[0]);
      const float minContrib = lContrib[1];
      const float maxContrib = lContrib[2];
      qa.fQualityAssuranceList = new TList();
      qa.fQualityAssuranceList->SetName("QualityAssurance");
      qa.fQualityAssuranceList->SetOwner(true);
      fBaseList->Add(qa.fQualityAssuranceList);

      if (doprocessRec || doprocessRecSim) {
        const auto& refMultBins = lEventHistBins[eHistRefMult];
        for (int i = 0; i < eBeforeAfterCuts_N; ++i) {
          qa.fHistMultNContrib[i] = new TH2F(Form("fHistMultNContrib[%s]", CutsNames[i]), Form("refMult vs. nContributors %s cuts", CutsNames[i]), static_cast<int>(refMultBins[0]), refMultBins[1], refMultBins[2], nBinsContrib, minContrib, maxContrib);
          qa.fHistMultNContrib[i]->GetXaxis()->SetTitle(Form("Reference Multiplicity (%s)", MultiplicityTablesNames[tc.fMultiplicityTables]));
          qa.fHistMultNContrib[i]->GetYaxis()->SetTitle("Number of contributors");
          qa.fQualityAssuranceList->Add(qa.fHistMultNContrib[i]);
        }

        if (doprocessRecSim) {
          const auto& centBins = lEventHistBins[eHistCentrality];
          for (int i = 0; i < eBeforeAfterCuts_N; ++i) {
            qa.fHistCentralityRecSim[i] = new TH2F(Form("fHistCentralityRecSim[%s]", CutsNames[i]), Form("Centrality Rec vs Sim %s cuts", CutsNames[i]), static_cast<int>(centBins[0]), centBins[1], centBins[2], static_cast<int>(centBins[0]), centBins[1], centBins[2]);
            qa.fHistCentralityRecSim[i]->GetXaxis()->SetTitle("Centrality (reconstructed)");
            qa.fHistCentralityRecSim[i]->GetYaxis()->SetTitle("Centrality (simulated)");
            qa.fQualityAssuranceList->Add(qa.fHistCentralityRecSim[i]);
          }
        } // end of if (doprocessRecSim) {
      } // end of if (doprocessRec || doprocessRecSim) {
    } // end of if (tc.fQualityAssuranceSwitch) {

    // *) Book InternalValidation TLists:
    if (tc.fInternalValidationSwitch) {
      intval.fIntValFun = new TF1("fIntValFun", "1 + 2*0.04*cos(1*(x - [0])) + 2*[1]*cos(2*(x - [0])) + 2*[2]*cos(3*(x - [0])) + 2*[3]*cos(4*(x - [0]))", 0., math::TwoPI);
      for (std::size_t i = 1; i <= tc.fIVHarmonics.size(); ++i) {
        intval.fIntValFun->SetParameter(i, tc.fIVHarmonics[i - 1]);
      }

      intval.fIntValList = new TList();
      intval.fIntValList->SetName("InternalValidation");
      intval.fIntValList->SetOwner(true);
      fBaseList->Add(intval.fIntValList);

      for (int i = 0; i < eIntValHarmonics_N; ++i) {
        intval.fProfIntVal[i] = new TProfile(Form("[%s]", IntValNames[i]), Form("Internal Validation for %s with True Value %f", IntValTitles[i], std::pow(tc.fIVHarmonics[i], 2)), nDefaultCentBins, defaultCentBoundaries.data());
        intval.fProfIntVal[i]->SetStats(false);
        intval.fProfIntVal[i]->GetXaxis()->SetTitle("Centrality (calculated from I.P.)");
        intval.fProfIntVal[i]->GetYaxis()->SetTitle(IntValTitles[i]);
        intval.fProfIntVal[i]->Sumw2();
        intval.fIntValList->Add(intval.fProfIntVal[i]);
      }
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
  return WorkflowSpec{adaptAnalysisTask<MultiparticleCorrelationsMei>(cfgc)};
} // WorkflowSpec...
