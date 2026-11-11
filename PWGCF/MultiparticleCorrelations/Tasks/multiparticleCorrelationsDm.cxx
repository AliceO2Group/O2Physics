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
#include <TF1.h>
#include <TFile.h>
#include <TGrid.h>
#include <TH1.h>
#include <TH2.h>
#include <TList.h>
#include <TObject.h>
#include <TProfile.h>
#include <TString.h>
#include <TSystem.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
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
  ECharge,       // only reconstructed
  EDcaXY,        // only reconstructed, from aod::TracksDCA
  EDcaZ,         // only reconstructed, from aod::TracksDCA
  ETpcNClsFound, // only reconstructed, from aod::TracksExtra
  EItsNCls,      // only reconstructed, from aod::TracksExtra
  EParticleHistograms_N
};

enum EBeforeAfter {
  EBefore = 0,
  EAfter,
  EBeforeAfter_N
};

enum EOccupancyEstimator {
  ETrackOccupancy = 0, // default: number of ITS tracks in the TPC drift time window
  EFT0COccupancy,      // FT0C amplitude in the TPC drift time window
  EOccupancyEstimator_N
};

// *) Main task:
struct MultiparticleCorrelationsDm { // this name is used in lower-case format, to name the TDirectoryFile in AnalysisResults.root

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
  Configurable<bool> cfFillQA{"cfFillQA", false, "book and fill QA histograms (only for development and validation)"};                                            // main switch for all QA histograms
  Configurable<bool> cfCentralityRecVsSim{"cfCentralityRecVsSim", false, "book and fill QA histogram for 2D reconstructed vs. simulated centralities (MC only)"}; // switch for one QA histogram
  Configurable<bool> cfMultvsN{"cfMultvsN", false, "book and fill QA histogram for reference multiplicity vs. nr. of contributors (expect linear)"};              // switch for one QA histogram
  Configurable<bool> cfFillCentFT0CVsFT0M{"cfFillCentFT0CVsFT0M", true, "QA: centFT0C vs. centFT0M, before and after cuts"};
  Configurable<std::vector<float>> cfValidationCentrality{"cfValidationCentrality", {60., 70.}, "validation: centrality class {min, max} for the selected-track multiplicity"};

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

  // 2D correlation cuts:
  Configurable<bool> cfUseCentralityCorrelationCut{"cfUseCentralityCorrelationCut", false, "reject events in which centFT0C and centFT0M differ too much (relative)"};
  Configurable<std::vector<float>> cfCentralityCorrelationCut{"cfCentralityCorrelationCut", {0.2, 5.}, "{maximum relative difference, floor (%) for the scale}"};
  Configurable<bool> cfUseNContribVsRefMultCut{"cfUseNContribVsRefMultCut", false, "reject events outside the linear bands in nContrib vs. refMult"};
  Configurable<std::vector<float>> cfNContribVsRefMultLower{"cfNContribVsRefMultLower", {0., 0.}, "lower band: nContrib = p0 + p1 * refMult, {p0, p1}"};
  Configurable<std::vector<float>> cfNContribVsRefMultUpper{"cfNContribVsRefMultUpper", {5000., 0.}, "upper band: nContrib = p0 + p1 * refMult, {p0, p1}"};

  // Occupancy cut:
  Configurable<bool> cfUseOccupancyCut{"cfUseOccupancyCut", true, "apply the occupancy event cut"};
  Configurable<std::string> cfOccupancyEstimator{"cfOccupancyEstimator", "Track", "occupancy estimator: Track (default), FT0C"};
  Configurable<std::vector<float>> cfOccupancy{"cfOccupancy", {0., 6000.}, "occupancy range: {min, max}, with convention: min <= occupancy < max"};
  Configurable<bool> cfFillOccupancy{"cfFillOccupancy", true, "QA: both occupancy estimators, before and after cuts"};
  Configurable<std::vector<float>> cfTrackOccupancyBins{"cfTrackOccupancyBins", {300, 0., 30000.}, "nBins, min, max for the track occupancy"};
  Configurable<std::vector<float>> cfFT0COccupancyBins{"cfFT0COccupancyBins", {300, 0., 300000.}, "nBins, min, max for the FT0C occupancy"};

  Configurable<std::string> cfCentralityEstimator{"cfCentralityEstimator", "FT0C", "centrality estimator: FT0C, FT0M, FT0A"};                                  // centrality estimator with valid values: FT0C (default)
  Configurable<std::string> cfReferenceMultiplicityEstimator{"cfReferenceMultiplicityEstimator", "FT0C", "reference multiplicity estimator: TPC, FT0C, FT0M"}; // ref. multiplicity estimator with valid values: TPC (default) FT0C Ft0M

  // sigma_inel in fm^2 (1 b = 100 fm^2): Pb+Pb at 5.36 TeV = 7.71 b = 771 fm^2 (PH2278 Lecture 2, slide 33), change for other system or energy
  Configurable<float> cfSigmaInel{"cfSigmaInel", 771., "sigma_inel in fm^2 (1 b = 100 fm^2), Pb+Pb 5.36 TeV = 771"};

  Configurable<std::vector<float>> cfPtBins{"cfPtBins", {200, 0., 10.}, "nPtBins, ptMin, ptMax"};
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
  Configurable<bool> cfUseEtaCut{"cfUseEtaCut", true, "apply the eta particle cut (pseudorapidity)"};
  Configurable<std::vector<float>> cfEta{"cfEta", {-0.8, 0.8}, "eta range: {min, max}, with convention: min <= eta < max"};

  // track-quality and charge cuts (only for reconstructed tracks):
  Configurable<int> cfCharge{"cfCharge", 0, "charge: 0 = both charges (default), 1 = only positive, -1 = only negative"};
  Configurable<bool> cfUseTpcNClsFoundCut{"cfUseTpcNClsFoundCut", true, "apply the cut on the number of TPC clusters found"};
  Configurable<int> cfTpcNClsFoundMin{"cfTpcNClsFoundMin", 70, "keep tracks with tpcNClsFound >= this value"};
  Configurable<bool> cfUseItsNClsCut{"cfUseItsNClsCut", true, "apply the cut on the number of ITS clusters"};
  Configurable<int> cfItsNClsMin{"cfItsNClsMin", 5, "keep tracks with itsNCls >= this value"};
  Configurable<bool> cfUseDcaXYCut{"cfUseDcaXYCut", true, "apply the cut on |dcaXY|"};
  Configurable<float> cfDcaXYMax{"cfDcaXYMax", 2.4, "keep tracks with |dcaXY| < this value [cm]"};
  Configurable<bool> cfUseDcaZCut{"cfUseDcaZCut", true, "apply the cut on |dcaZ|"};
  Configurable<float> cfDcaZMax{"cfDcaZMax", 3.2, "keep tracks with |dcaZ| < this value [cm]"};

  // binning of the track-quality and charge histograms:
  Configurable<std::vector<float>> cfChargeBins{"cfChargeBins", {3, -1.5, 1.5}, "nBins, min, max"};
  Configurable<std::vector<float>> cfDcaXYBins{"cfDcaXYBins", {600, -3., 3.}, "nBins, min, max [cm]"};
  Configurable<std::vector<float>> cfDcaZBins{"cfDcaZBins", {800, -4., 4.}, "nBins, min, max [cm]"};
  Configurable<std::vector<float>> cfTpcNClsFoundBins{"cfTpcNClsFoundBins", {160, 0., 160.}, "nBins, min, max"};
  Configurable<std::vector<float>> cfItsNClsBins{"cfItsNClsBins", {8, 0., 8.}, "nBins, min, max"};

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
    bool fCentralityRecVsSim = false;        // QA: centrality rec vs. sim
    bool fMultvsN = false;                   // QA: reference multiplicity vs. nr. of contributors
    bool fUseVertexZCut = true;
    bool fUsePtCut = true;
    bool fUseTimeRange = true;
    bool fUseRof = true;
    bool fUsePileUp = true;
    bool fUseITSTPC = true;
    bool fUseITSLayers = true;
    bool fUseZFT0vsPV = true;
    bool fUsePrevRof = true;

    bool fUseCentralityCorrelationCut = false;
    std::array<float, 2> fCentralityCorrelationCut = {0.2, 5.}; // {maximum relative difference, floor (%)}
    bool fUseNContribVsRefMultCut = false;
    bool fFillCentFT0CVsFT0M = true;
    std::array<float, 2> fValidationCentrality = {60., 70.}; // {min, max} (%)

    bool fUseOccupancyCut = true;
    int fOccupancyEstimator = ETrackOccupancy;
    std::array<float, 2> fOccupancy = {0., 6000.}; // {min, max}
    bool fFillOccupancy = true;

    int fCentralityEstimator = ECentFT0C;
    int fReferenceMultiplicityEstimator = ERefMultTPC;

    float fSigmaInel = 771.;                     // sigma inel [fm^2]
    std::array<float, 2> fVertexZ = {-10., 10.}; // {min, max} [cm]
    std::array<float, 2> fPt = {0.2, 5.0};       // {min, max} [GeV]
    bool fUseEtaCut = true;
    std::array<float, 2> fEta = {-0.8, 0.8}; // {min, max}
    int fCharge = 0;                         // 0 = both, 1 = positive, -1 = negative
    bool fUseTpcNClsFoundCut = true;
    int fTpcNClsFoundMin = 70;
    bool fUseItsNClsCut = true;
    int fItsNClsMin = 5;
    bool fUseDcaXYCut = true;
    float fDcaXYMax = 2.4; // [cm]
    bool fUseDcaZCut = true;
    float fDcaZMax = 3.2;        // [cm]
    int fCurrentRunNumber = -1;  // run number of the weights currently in use
    bool fUsePhiWeights = false; // load phi weights
    bool fUsePtWeights = false;  // load pt weights
    bool fUseEtaWeights = false; // load eta weights
    bool fVerbose = false;       // print the values / collision

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

    TH1F* fHistNSelected = nullptr; // selected tracks per event (after cuts)

  } ec; // you have to prepend "ec." for all objects name in this group later in the code

  // **) QA histograms (booked and filled only if cfFillQA is true):

  struct QAHistograms {
    TList* fQAHistogramsList = nullptr;
    TH2F* fHistCentralityRecVsSim = nullptr;                                               // x = centrality sim (from impact parameter), y = centrality rec
    std::array<TH2F*, EBeforeAfter_N> fHistRefMultVsNContrib{};                            // x = reference multiplicity, y = number of vertex contributors [before, after]
    std::array<TH2F*, EBeforeAfter_N> fHistCentFT0CVsFT0M{};                               // x = centFT0M, y = centFT0C [before, after]
    TH1F* fHistSelectedMultiplicity = nullptr;                                             // validation: number of selected tracks per event in one centrality class
    std::array<std::array<TH1F*, EBeforeAfter_N>, EOccupancyEstimator_N> fHistOccupancy{}; // [estimator][before, after]
    std::array<TH2F*, EBeforeAfter_N> fHistOccupancyTrackVsFT0C{};                         // x = FT0C occupancy, y = track occupancy [before, after]
  } qa;                                                                                    // you have to prepend "qa." for all objects name in this group later in the code

  // **) Correlations:
  struct Correlations {
    TList* fCorrelationsList = nullptr;
    TProfile* fProfTwo = nullptr;  // <<2>> = v2{2}^2 vs. centrality (rec)
    TProfile* fProfFour = nullptr; // <<4>> vs. centrality (rec)
  } cc;

  // **) Q-vectors of the current event: [harmonic 0 - 4][weight power 0 - 4], and the random generator for the subsamples
  static constexpr int KMaxHarmonic = 4;
  static constexpr int KMaxPower = 4;
  static constexpr int NRangeValues = 2; // {min, max} or {p0, p1}
  std::array<std::array<std::complex<double>, KMaxPower + 1>, KMaxHarmonic + 1> fQvector{};

  // **) Histogram with weights from an external file
  TH1F* histPhiWeights = nullptr; // phi weights for the current run
  TH1F* histPtWeights = nullptr;  // pt weights for the current run
  TH1F* histEtaWeights = nullptr; // eta weights for the current run

  // **) Linear bands for the geometric cut in nContrib vs. refMult:
  TF1* fNContribVsRefMultLower = nullptr;
  TF1* fNContribVsRefMultUpper = nullptr;

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

  // returns occupancy from the estimator selected in json
  template <typename T>
  float getOccupancy(T const& collision)
  {
    switch (tc.fOccupancyEstimator) {
      case ETrackOccupancy:
        return collision.trackOccupancyInTimeRange();
      case EFT0COccupancy:
        return collision.ft0cOccupancyInTimeRange();
      default:
        LOGF(fatal, "Unknown occupancy estimator");
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

  // add one particle with azimuthal angle phi and weight w to the Q-vectors: Q(n,p) += w^p * exp(i n phi)
  void fillQvector(double phi, double w)
  {
    for (int n = 0; n <= KMaxHarmonic; n++) {
      const std::complex<double> e = std::polar(1., n * phi);
      double wp = w; // w^1, then multiplied by w for each higher power
      for (int p = 1; p <= KMaxPower; p++) {
        fQvector[n][p] += wp * e;
        wp *= w;
      }
    }
  }

  // Q(n,p), including negative harmonics: Q(-n,p) = conjugate of Q(n,p)
  std::complex<double> qv(int n, int p)
  {
    if (std::abs(n) > KMaxHarmonic || p < 1 || p > KMaxPower) {
      LOGF(fatal, "qv(%d, %d) is outside of the stored Q-vectors", n, p);
    }
    return n >= 0 ? fQvector[n][p] : std::conj(fQvector[-n][p]);
  }

  // generic 2-particle correlation, numerator of <exp[i(n1*phi1+n2*phi2)]> (same as Two() in the toy MC)
  std::complex<double> twoParticle(int n1, int n2)
  {
    return qv(n1, 1) * qv(n2, 1) - qv(n1 + n2, 2);
  }

  // generic 4-particle correlation, numerator of <exp[i(n1*phi1+n2*phi2+n3*phi3+n4*phi4)]>, same as Four() in the toy MC
  std::complex<double> fourParticle(int n1, int n2, int n3, int n4)
  {
    return qv(n1, 1) * qv(n2, 1) * qv(n3, 1) * qv(n4, 1) - qv(n1 + n2, 2) * qv(n3, 1) * qv(n4, 1) - qv(n2, 1) * qv(n1 + n3, 2) * qv(n4, 1) - qv(n1, 1) * qv(n2 + n3, 2) * qv(n4, 1) + 2. * qv(n1 + n2 + n3, 3) * qv(n4, 1) - qv(n2, 1) * qv(n3, 1) * qv(n1 + n4, 2) + qv(n2 + n3, 2) * qv(n1 + n4, 2) - qv(n1, 1) * qv(n3, 1) * qv(n2 + n4, 2) + qv(n1 + n3, 2) * qv(n2 + n4, 2) + 2. * qv(n3, 1) * qv(n1 + n2 + n4, 3) - qv(n1, 1) * qv(n2, 1) * qv(n3 + n4, 2) + qv(n1 + n2, 2) * qv(n3 + n4, 2) + 2. * qv(n2, 1) * qv(n1 + n3 + n4, 3) + 2. * qv(n1, 1) * qv(n2 + n3 + n4, 3) - 6. * qv(n1 + n2 + n3 + n4, 4);
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

    // *) Occupancy cut, convention: min <= occupancy < max (events without occupancy information, value -1, are rejected for min = 0)
    if (tc.fUseOccupancyCut) {
      const float occupancy = getOccupancy(collision);
      if (!(tc.fOccupancy[0] <= occupancy && occupancy < tc.fOccupancy[1])) {
        return false; // if at least one event cut is not satisfied
      }
    }

    // *) Direct correlation cut: centFT0C and centFT0M must agree within a relative difference
    //    scale = mean of both, but at least the floor, so that small normal differences in the most central events are not huge relative differences
    if (tc.fUseCentralityCorrelationCut) {
      const float centC = collision.centFT0C();
      const float centM = collision.centFT0M();
      const float scale = std::max(0.5f * (centC + centM), tc.fCentralityCorrelationCut[1]);
      if (std::abs(centC - centM) / scale > tc.fCentralityCorrelationCut[0]) {
        return false; // if at least one event cut is not satisfied
      }
    }

    // *) Geometric cut: nContrib must lie between the lower and upper linear bands at this refMult
    if (tc.fUseNContribVsRefMultCut) {
      const float refMult = getReferenceMultiplicity(collision);
      const float nContrib = collision.numContrib();
      if (nContrib < fNContribVsRefMultLower->Eval(refMult) || nContrib > fNContribVsRefMultUpper->Eval(refMult)) {
        return false; // if at least one event cut is not satisfied
      }
    }

    return true; // all event cuts passed
  }

  // fill charge and track-quality histograms of one reconstructed track, ba = EBefore or EAfter
  template <typename T>
  void fillTrackQualityHistograms(T const& track, EBeforeAfter ba)
  {
    pc.fParticleHistograms[ECharge][ERec][ba]->Fill(track.sign());
    pc.fParticleHistograms[EDcaXY][ERec][ba]->Fill(track.dcaXY());
    pc.fParticleHistograms[EDcaZ][ERec][ba]->Fill(track.dcaZ());
    pc.fParticleHistograms[ETpcNClsFound][ERec][ba]->Fill(track.tpcNClsFound());
    pc.fParticleHistograms[EItsNCls][ERec][ba]->Fill(track.itsNCls());
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

    // *) Charge and track-quality cuts: only for reconstructed tracks (MC particles have no sign, clusters or DCA)
    if constexpr (rs != ESim) {

      // *) charge cut: 0 = both charges, otherwise keep only the selected sign
      if (tc.fCharge != 0) {
        if (track.sign() != tc.fCharge) {
          return false; // if at least one particle cut is not satisfied
        }
      }

      // *) number of TPC clusters found
      if (tc.fUseTpcNClsFoundCut) {
        if (!(track.tpcNClsFound() >= tc.fTpcNClsFoundMin)) {
          return false; // if at least one particle cut is not satisfied
        }
      }

      // *) number of ITS clusters
      if (tc.fUseItsNClsCut) {
        if (!(track.itsNCls() >= tc.fItsNClsMin)) {
          return false; // if at least one particle cut is not satisfied
        }
      }

      // *) DCA to the primary vertex in the transverse plane
      if (tc.fUseDcaXYCut) {
        if (!(std::abs(track.dcaXY()) < tc.fDcaXYMax)) {
          return false; // if at least one particle cut is not satisfied
        }
      }

      // *) DCA to the primary vertex along the beam
      if (tc.fUseDcaZCut) {
        if (!(std::abs(track.dcaZ()) < tc.fDcaZMax)) {
          return false; // if at least one particle cut is not satisfied
        }
      }

    } // end of if constexpr (rs != ESim)

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
      if (tc.fFillQA && tc.fMultvsN) {
        qa.fHistRefMultVsNContrib[ba]->Fill(getReferenceMultiplicity(collision), collision.numContrib());
      }

      // QA: centFT0C vs. centFT0M
      if (tc.fFillQA && tc.fFillCentFT0CVsFT0M) {
        qa.fHistCentFT0CVsFT0M[ba]->Fill(collision.centFT0M(), collision.centFT0C());
      }

      // QA: both occupancy estimators, and their correlation
      if (tc.fFillQA && tc.fFillOccupancy) {
        qa.fHistOccupancy[ETrackOccupancy][ba]->Fill(collision.trackOccupancyInTimeRange());
        qa.fHistOccupancy[EFT0COccupancy][ba]->Fill(collision.ft0cOccupancyInTimeRange());
        qa.fHistOccupancyTrackVsFT0C[ba]->Fill(collision.ft0cOccupancyInTimeRange(), collision.trackOccupancyInTimeRange());
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
          if (tc.fFillQA && tc.fCentralityRecVsSim && ba == EAfter) {
            qa.fHistCentralityRecVsSim->Fill(centSim, centRec);
          } // if (tc.fFillQA && tc.fCentralityRecVsSim && ba == EAfter)

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
    fQvector = {};     // reset the Q-vectors for this event
    int nSelected = 0; // number of tracks that pass all particle cuts

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
        fillTrackQualityHistograms(track, EBefore);

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
        fillTrackQualityHistograms(track, EAfter);

        // add this particle to the Q-vector
        const double w = getParticleWeight(track.phi(), track.pt(), track.eta());
        fillQvector(track.phi(), w);
        nSelected++;

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

    // validation: selected-track multiplicity, only in the chosen centrality class
    if constexpr (rs == ERec || rs == ERecAndSim) {
      if (tc.fFillQA) {
        const float cent = getCentrality(collision);
        ec.fHistNSelected->Fill(nSelected); // after the track loop
        if (tc.fValidationCentrality[0] <= cent && cent < tc.fValidationCentrality[1]) {
          qa.fHistSelectedMultiplicity->Fill(nSelected);
        }
      }
    }

    // <2> and <4> of this event from the Q-vectors, filled vs. centrality, and the trending profiles
    if constexpr (rs == ERec || rs == ERecAndSim) {
      const float cent = getCentrality(collision);

      // <2>, with event weight = number of weighted pairs (= M(M-1) for unit weights)
      const double twoWeight = twoParticle(0, 0).real();
      double two = 0.;
      if (twoWeight > 0.) {
        two = twoParticle(-2, 2).real() / twoWeight;
        cc.fProfTwo->Fill(cent, two, twoWeight);
      }

      // <4>, with event weight = number of weighted quadruplets (= M(M-1)(M-2)(M-3) for unit weights)
      const double fourWeight = fourParticle(0, 0, 0, 0).real();
      if (fourWeight > 0.) {
        const double four = fourParticle(-2, -2, 2, 2).real() / fourWeight;
        cc.fProfFour->Fill(cent, four, fourWeight);
      }
    }

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
    tc.fCentralityRecVsSim = cfCentralityRecVsSim;
    tc.fMultvsN = cfMultvsN;
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

    // charge and track-quality cuts: copy the switches and values
    tc.fCharge = cfCharge;

    if (tc.fCharge != 0 && tc.fCharge != 1 && tc.fCharge != -1) {
      LOGF(fatal, "cfCharge must be 0 (both), 1 (positive) or -1 (negative), it is %d", tc.fCharge);
    }
    tc.fUseTpcNClsFoundCut = cfUseTpcNClsFoundCut;
    tc.fTpcNClsFoundMin = cfTpcNClsFoundMin;
    tc.fUseItsNClsCut = cfUseItsNClsCut;
    tc.fItsNClsMin = cfItsNClsMin;
    tc.fUseDcaXYCut = cfUseDcaXYCut;
    tc.fDcaXYMax = cfDcaXYMax;
    tc.fUseDcaZCut = cfUseDcaZCut;
    tc.fDcaZMax = cfDcaZMax;

    // technical event cuts: copy the switches
    tc.fUseTimeRange = cfUseTimeRange;
    tc.fUseRof = cfUseRof;
    tc.fUsePileUp = cfUsePileUp;
    tc.fUseITSTPC = cfUseITSTPC;
    tc.fUseITSLayers = cfUseITSLayers;
    tc.fUseZFT0vsPV = cfUseZFT0vsPV;
    tc.fUsePrevRof = cfUsePrevRof;

    // 2D correlation cuts and validation: copy the switches and values
    tc.fUseCentralityCorrelationCut = cfUseCentralityCorrelationCut;
    tc.fUseNContribVsRefMultCut = cfUseNContribVsRefMultCut;
    tc.fFillCentFT0CVsFT0M = cfFillCentFT0CVsFT0M;
    std::vector<float> lCentCorr = cfCentralityCorrelationCut.value;
    std::vector<float> lValCent = cfValidationCentrality.value;
    if (lCentCorr.size() != NRangeValues || lValCent.size() != NRangeValues || cfNContribVsRefMultLower.value.size() != NRangeValues || cfNContribVsRefMultUpper.value.size() != NRangeValues) {
      LOGF(fatal, "cfCentralityCorrelationCut, cfValidationCentrality, cfNContribVsRefMultLower and cfNContribVsRefMultUpper need exactly 2 values each");
    }
    tc.fCentralityCorrelationCut[0] = lCentCorr[0];
    tc.fCentralityCorrelationCut[1] = lCentCorr[1];
    tc.fValidationCentrality[0] = lValCent[0];
    tc.fValidationCentrality[1] = lValCent[1];

    // occupancy cut: copy the switch and range, translate the estimator string into the enum
    tc.fUseOccupancyCut = cfUseOccupancyCut;
    tc.fFillOccupancy = cfFillOccupancy;
    std::vector<float> lOccupancy = cfOccupancy.value;
    if (lOccupancy.size() != NRangeValues) {
      LOGF(fatal, "cfOccupancy needs exactly 2 values {min, max}");
    }
    tc.fOccupancy[0] = lOccupancy[0];
    tc.fOccupancy[1] = lOccupancy[1];
    if (cfOccupancyEstimator.value == "Track") {
      tc.fOccupancyEstimator = ETrackOccupancy;
    } else if (cfOccupancyEstimator.value == "FT0C") {
      tc.fOccupancyEstimator = EFT0COccupancy;
    } else {
      LOGF(fatal, "Unknown occupancy estimator in json: %s", cfOccupancyEstimator.value.c_str());
    }

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

    // *) Linear bands for the geometric cut, saved in the output for book-keeping
    fNContribVsRefMultLower = new TF1("fNContribVsRefMultLower", "pol1", 0., 1.e6);
    fNContribVsRefMultLower->SetParameters(cfNContribVsRefMultLower.value[0], cfNContribVsRefMultLower.value[1]);
    fNContribVsRefMultUpper = new TF1("fNContribVsRefMultUpper", "pol1", 0., 1.e6);
    fNContribVsRefMultUpper->SetParameters(cfNContribVsRefMultUpper.value[0], cfNContribVsRefMultUpper.value[1]);

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

    cc.fProfFour = new TProfile("fProfFour", "<<4>> vs. centrality", centBoundaries.size() - 1, centBoundaries.data());
    cc.fProfFour->GetXaxis()->SetTitle("centrality (%)");
    cc.fProfFour->GetYaxis()->SetTitle("#LT#LT4#GT#GT");
    cc.fProfFour->Sumw2();
    cc.fCorrelationsList->Add(cc.fProfFour);

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
        pc.fParticleHistograms[EPt][ERec][ba] = new TH1F(Form("fHistPt[ERec][%s]", baNames[ba]), Form("p_{T}, reconstructed (%s cuts)", baNames[ba]), nBins, min, max);
        pc.fParticleHistograms[EPt][ERec][ba]->GetXaxis()->SetTitle("p_{T} (GeV/#it{c})");
        pc.fParticleHistograms[EPt][ERec][ba]->GetYaxis()->SetTitle("counts");
        pc.fParticleHistogramsList->Add(pc.fParticleHistograms[EPt][ERec][ba]);

        pc.fParticleHistograms[EPhi][ERec][ba] = new TH1F(Form("fHistPhi[ERec][%s]", baNames[ba]), Form("#varphi, reconstructed (%s cuts)", baNames[ba]), nPhiBins, phiMin, phiMax);
        pc.fParticleHistograms[EPhi][ERec][ba]->GetXaxis()->SetTitle("#varphi (rad)");
        pc.fParticleHistograms[EPhi][ERec][ba]->GetYaxis()->SetTitle("counts");
        pc.fParticleHistogramsList->Add(pc.fParticleHistograms[EPhi][ERec][ba]);

        // charge and track-quality histograms (only reconstructed)
        const std::array<EParticleHistograms, 5> qualityTypes = {ECharge, EDcaXY, EDcaZ, ETpcNClsFound, EItsNCls};
        constexpr int NQualityHistograms = 5; // charge, dcaXY, dcaZ, tpcNClsFound, itsNCls / Ladung, dcaXY, dcaZ, TPC-Cluster, ITS-Cluster
        const std::array<const char*, 5> qualityNames = {"Charge", "DcaXY", "DcaZ", "TpcNClsFound", "ItsNCls"};
        const std::array<const char*, 5> qualityTitles = {"charge", "DCA_{xy} (cm)", "DCA_{z} (cm)", "number of TPC clusters found", "number of ITS clusters"};
        const std::array<std::vector<float>, 5> qualityBins = {cfChargeBins.value, cfDcaXYBins.value, cfDcaZBins.value, cfTpcNClsFoundBins.value, cfItsNClsBins.value};
        for (int q = 0; q < NQualityHistograms; q++) {
          pc.fParticleHistograms[qualityTypes[q]][ERec][ba] = new TH1F(Form("fHist%s[ERec][%s]", qualityNames[q], baNames[ba]), Form("%s, reconstructed (%s cuts)", qualityTitles[q], baNames[ba]), static_cast<int>(qualityBins[q][0]), qualityBins[q][1], qualityBins[q][2]);
          pc.fParticleHistograms[qualityTypes[q]][ERec][ba]->GetXaxis()->SetTitle(qualityTitles[q]);
          pc.fParticleHistograms[qualityTypes[q]][ERec][ba]->GetYaxis()->SetTitle("counts");
          pc.fParticleHistogramsList->Add(pc.fParticleHistograms[qualityTypes[q]][ERec][ba]);
        }
        // ev
        ec.fEventHistograms[EMultiplicity][ERec][ba] = new TH1F(Form("fHistMultiplicity[ERec][%s]", baNames[ba]), Form("multiplicity, reconstructed (%s cuts)", baNames[ba]), nMultBins, multMin, multMax);
        ec.fEventHistograms[EMultiplicity][ERec][ba]->GetXaxis()->SetTitle("number of tracks");
        ec.fEventHistograms[EMultiplicity][ERec][ba]->GetYaxis()->SetTitle("counts");
        ec.fEventHistogramsList->Add(ec.fEventHistograms[EMultiplicity][ERec][ba]);

        // reference multiplicity (only reconstructed)
        ec.fEventHistograms[EReferenceMultiplicity][ERec][ba] = new TH1F(Form("fHistReferenceMultiplicity[ERec][%s]", baNames[ba]), Form("reference multiplicity, reconstructed (%s cuts)", baNames[ba]), nRefMultBins, refMultMin, refMultMax);
        ec.fEventHistograms[EReferenceMultiplicity][ERec][ba]->GetXaxis()->SetTitle(Form("reference multiplicity (%s)", cfReferenceMultiplicityEstimator.value.c_str()));
        ec.fEventHistograms[EReferenceMultiplicity][ERec][ba]->GetYaxis()->SetTitle("counts");
        ec.fEventHistogramsList->Add(ec.fEventHistograms[EReferenceMultiplicity][ERec][ba]);

        ec.fEventHistograms[ECentrality][ERec][ba] = new TH1F(Form("fHistCentrality[ERec][%s]", baNames[ba]), Form("centrality, reconstructed (%s cuts)", baNames[ba]), nCentBins, centMin, centMax);
        ec.fEventHistograms[ECentrality][ERec][ba]->GetXaxis()->SetTitle("centrality (%)");
        ec.fEventHistograms[ECentrality][ERec][ba]->GetYaxis()->SetTitle("counts");
        ec.fEventHistogramsList->Add(ec.fEventHistograms[ECentrality][ERec][ba]);

        ec.fEventHistograms[EVertexX][ERec][ba] = new TH1F(Form("fHistVertexX[ERec][%s]", baNames[ba]), Form("V_{x}, reconstructed (%s cuts)", baNames[ba]), nvxBins, vxMin, vxMax);
        ec.fEventHistograms[EVertexX][ERec][ba]->GetXaxis()->SetTitle("V_{x} (cm)");
        ec.fEventHistograms[EVertexX][ERec][ba]->GetYaxis()->SetTitle("counts");
        ec.fEventHistogramsList->Add(ec.fEventHistograms[EVertexX][ERec][ba]);

        ec.fEventHistograms[EVertexY][ERec][ba] = new TH1F(Form("fHistVertexY[ERec][%s]", baNames[ba]), Form("V_{y}, reconstructed (%s cuts)", baNames[ba]), nvyBins, vyMin, vyMax);
        ec.fEventHistograms[EVertexY][ERec][ba]->GetXaxis()->SetTitle("V_{y} (cm)");
        ec.fEventHistograms[EVertexY][ERec][ba]->GetYaxis()->SetTitle("counts");
        ec.fEventHistogramsList->Add(ec.fEventHistograms[EVertexY][ERec][ba]);

        ec.fEventHistograms[EVertexZ][ERec][ba] = new TH1F(Form("fHistVertexZ[ERec][%s]", baNames[ba]), Form("V_{z}, reconstructed (%s cuts)", baNames[ba]), nvzBins, vzMin, vzMax);
        ec.fEventHistograms[EVertexZ][ERec][ba]->GetXaxis()->SetTitle("V_{z} (cm)");
        ec.fEventHistograms[EVertexZ][ERec][ba]->GetYaxis()->SetTitle("counts");
        ec.fEventHistogramsList->Add(ec.fEventHistograms[EVertexZ][ERec][ba]);
      }
      if (bookSim) {
        pc.fParticleHistograms[EPt][ESim][ba] = new TH1F(Form("fHistPt[ESim][%s]", baNames[ba]), Form("p_{T}, simulated (%s cuts)", baNames[ba]), nBins, min, max);
        pc.fParticleHistograms[EPt][ESim][ba]->GetXaxis()->SetTitle("p_{T} (GeV/#it{c})");
        pc.fParticleHistograms[EPt][ESim][ba]->GetYaxis()->SetTitle("counts");
        pc.fParticleHistogramsList->Add(pc.fParticleHistograms[EPt][ESim][ba]);

        pc.fParticleHistograms[EPhi][ESim][ba] = new TH1F(Form("fHistPhi[ESim][%s]", baNames[ba]), Form("#varphi, simulated (%s cuts)", baNames[ba]), nPhiBins, phiMin, phiMax);
        pc.fParticleHistograms[EPhi][ESim][ba]->GetXaxis()->SetTitle("#varphi (rad)");
        pc.fParticleHistograms[EPhi][ESim][ba]->GetYaxis()->SetTitle("counts");
        pc.fParticleHistogramsList->Add(pc.fParticleHistograms[EPhi][ESim][ba]);
        // ev
        ec.fEventHistograms[EMultiplicity][ESim][ba] = new TH1F(Form("fHistMultiplicity[ESim][%s]", baNames[ba]), Form("multiplicity, simulated (%s cuts)", baNames[ba]), nMultBins, multMin, multMax);
        ec.fEventHistograms[EMultiplicity][ESim][ba]->GetXaxis()->SetTitle("number of tracks");
        ec.fEventHistograms[EMultiplicity][ESim][ba]->GetYaxis()->SetTitle("counts");
        ec.fEventHistogramsList->Add(ec.fEventHistograms[EMultiplicity][ESim][ba]);
        // sim from impact parameter
        ec.fEventHistograms[ECentrality][ESim][ba] = new TH1F(Form("fHistCentrality[ESim][%s]", baNames[ba]), Form("centrality from impact parameter, simulated (%s cuts)", baNames[ba]), nCentBins, centMin, centMax);
        ec.fEventHistograms[ECentrality][ESim][ba]->GetXaxis()->SetTitle("centrality (%)");
        ec.fEventHistograms[ECentrality][ESim][ba]->GetYaxis()->SetTitle("counts");
        ec.fEventHistogramsList->Add(ec.fEventHistograms[ECentrality][ESim][ba]);
        // impact parameter, only simulated
        ec.fEventHistograms[EImpactParameter][ESim][ba] = new TH1F(Form("fHistImpactParameter[ESim][%s]", baNames[ba]), Form("impact parameter, simulated (%s cuts)", baNames[ba]), nBBins, bMin, bMax);
        ec.fEventHistograms[EImpactParameter][ESim][ba]->GetXaxis()->SetTitle("b (fm)");
        ec.fEventHistograms[EImpactParameter][ESim][ba]->GetYaxis()->SetTitle("counts");
        ec.fEventHistogramsList->Add(ec.fEventHistograms[EImpactParameter][ESim][ba]);

        ec.fEventHistograms[EVertexX][ESim][ba] = new TH1F(Form("fHistVertexX[ESim][%s]", baNames[ba]), Form("V_{x}, simulated (%s cuts)", baNames[ba]), nvxBins, vxMin, vxMax);
        ec.fEventHistograms[EVertexX][ESim][ba]->GetXaxis()->SetTitle("V_{x} (cm)");
        ec.fEventHistograms[EVertexX][ESim][ba]->GetYaxis()->SetTitle("counts");
        ec.fEventHistogramsList->Add(ec.fEventHistograms[EVertexX][ESim][ba]);

        ec.fEventHistograms[EVertexY][ESim][ba] = new TH1F(Form("fHistVertexY[ESim][%s]", baNames[ba]), Form("V_{y}, simulated (%s cuts)", baNames[ba]), nvyBins, vyMin, vyMax);
        ec.fEventHistograms[EVertexY][ESim][ba]->GetXaxis()->SetTitle("V_{y} (cm)");
        ec.fEventHistograms[EVertexY][ESim][ba]->GetYaxis()->SetTitle("counts");
        ec.fEventHistogramsList->Add(ec.fEventHistograms[EVertexY][ESim][ba]);

        ec.fEventHistograms[EVertexZ][ESim][ba] = new TH1F(Form("fHistVertexZ[ESim][%s]", baNames[ba]), Form("V_{z}, simulated (%s cuts)", baNames[ba]), nvzBins, vzMin, vzMax);
        ec.fEventHistograms[EVertexZ][ESim][ba]->GetXaxis()->SetTitle("V_{z} (cm)");
        ec.fEventHistograms[EVertexZ][ESim][ba]->GetYaxis()->SetTitle("counts");
        ec.fEventHistogramsList->Add(ec.fEventHistograms[EVertexZ][ESim][ba]);
      }
    } // end of for (int ba = 0; ba < 2; ba++)

    // number of selected tracks per event, for trending
    if (bookRec) {
      ec.fHistNSelected = new TH1F("fHistNSelected", "number of selected tracks per event (after cuts)", nMultBins, multMin, multMax);
      ec.fHistNSelected->GetXaxis()->SetTitle("number of selected tracks");
      ec.fHistNSelected->GetYaxis()->SetTitle("counts");
      ec.fEventHistogramsList->Add(ec.fHistNSelected);
    }

    // QA histograms, booked only if cfFillQA is true (and each one only if its own switch is on)
    if (tc.fFillQA) {

      // centrality rec vs. sim: needs reconstructed and simulated data in the same event, so only with processRecSim
      if (tc.fCentralityRecVsSim && doprocessRecSim) {
        qa.fHistCentralityRecVsSim = new TH2F("fHistCentralityRecVsSim", "Centrality: reconstructed vs. simulated", nCentBins, centMin, centMax, nCentBins, centMin, centMax);
        qa.fHistCentralityRecVsSim->GetXaxis()->SetTitle("centrality simulated (%)");
        qa.fHistCentralityRecVsSim->GetYaxis()->SetTitle("centrality reconstructed (%)");
        qa.fQAHistogramsList->Add(qa.fHistCentralityRecVsSim);
      }

      // reference multiplicity vs. number of contributors, before and after cuts: only reconstructed quantities
      if (tc.fMultvsN && bookRec) {
        std::vector<float> lContribBins = cfNContribBins.value;
        for (int ba = 0; ba < EBeforeAfter_N; ba++) {
          qa.fHistRefMultVsNContrib[ba] = new TH2F(Form("fHistRefMultVsNContrib[%s]", baNames[ba]), Form("reference multiplicity vs. number of contributors (%s cuts)", baNames[ba]), nRefMultBins, refMultMin, refMultMax, static_cast<int>(lContribBins[0]), lContribBins[1], lContribBins[2]);
          qa.fHistRefMultVsNContrib[ba]->GetXaxis()->SetTitle(Form("reference multiplicity (%s)", cfReferenceMultiplicityEstimator.value.c_str()));
          qa.fHistRefMultVsNContrib[ba]->GetYaxis()->SetTitle("number of contributors");
          qa.fQAHistogramsList->Add(qa.fHistRefMultVsNContrib[ba]);
        }
      }

      // centFT0C vs. centFT0M, before and after cuts
      if (tc.fFillCentFT0CVsFT0M && bookRec) {
        for (int ba = 0; ba < EBeforeAfter_N; ba++) {
          qa.fHistCentFT0CVsFT0M[ba] = new TH2F(Form("fHistCentFT0CVsFT0M[%s]", baNames[ba]), Form("centFT0C vs. centFT0M (%s cuts)", baNames[ba]), nCentBins, centMin, centMax, nCentBins, centMin, centMax);
          qa.fHistCentFT0CVsFT0M[ba]->GetXaxis()->SetTitle("centrality FT0M (%)");
          qa.fHistCentFT0CVsFT0M[ba]->GetYaxis()->SetTitle("centrality FT0C (%)");
          qa.fQAHistogramsList->Add(qa.fHistCentFT0CVsFT0M[ba]);
        }
      }

      // both occupancy estimators, before and after cuts, and their correlation
      if (tc.fFillOccupancy && bookRec) {
        std::vector<float> lTrackOccBins = cfTrackOccupancyBins.value;
        std::vector<float> lFT0COccBins = cfFT0COccupancyBins.value;
        for (int ba = 0; ba < EBeforeAfter_N; ba++) {
          qa.fHistOccupancy[ETrackOccupancy][ba] = new TH1F(Form("fHistOccupancy[ETrack][%s]", baNames[ba]), Form("track occupancy (%s cuts)", baNames[ba]), static_cast<int>(lTrackOccBins[0]), lTrackOccBins[1], lTrackOccBins[2]);
          qa.fHistOccupancy[ETrackOccupancy][ba]->GetXaxis()->SetTitle("track occupancy (ITS tracks in time range)");
          qa.fHistOccupancy[ETrackOccupancy][ba]->GetYaxis()->SetTitle("counts");
          qa.fQAHistogramsList->Add(qa.fHistOccupancy[ETrackOccupancy][ba]);
          qa.fHistOccupancy[EFT0COccupancy][ba] = new TH1F(Form("fHistOccupancy[EFT0C][%s]", baNames[ba]), Form("FT0C occupancy (%s cuts)", baNames[ba]), static_cast<int>(lFT0COccBins[0]), lFT0COccBins[1], lFT0COccBins[2]);
          qa.fHistOccupancy[EFT0COccupancy][ba]->GetXaxis()->SetTitle("FT0C occupancy (in time range)");
          qa.fQAHistogramsList->Add(qa.fHistOccupancy[EFT0COccupancy][ba]);
          qa.fHistOccupancyTrackVsFT0C[ba] = new TH2F(Form("fHistOccupancyTrackVsFT0C[%s]", baNames[ba]), Form("track occupancy vs. FT0C occupancy (%s cuts)", baNames[ba]), static_cast<int>(lFT0COccBins[0]), lFT0COccBins[1], lFT0COccBins[2], static_cast<int>(lTrackOccBins[0]), lTrackOccBins[1], lTrackOccBins[2]);
          qa.fHistOccupancyTrackVsFT0C[ba]->GetXaxis()->SetTitle("FT0C occupancy (in time range)");
          qa.fHistOccupancyTrackVsFT0C[ba]->GetYaxis()->SetTitle("track occupancy (ITS tracks in time range)");
          qa.fQAHistogramsList->Add(qa.fHistOccupancyTrackVsFT0C[ba]);
        }
      }

      // validation: number of selected tracks per event in one centrality class
      if (bookRec) {
        qa.fHistSelectedMultiplicity = new TH1F("fHistSelectedMultiplicity", Form("selected tracks per event, centrality %.0f-%.0f %%", tc.fValidationCentrality[0], tc.fValidationCentrality[1]), nMultBins, multMin, multMax);
        qa.fHistSelectedMultiplicity->GetXaxis()->SetTitle("number of selected tracks");
        qa.fHistSelectedMultiplicity->GetYaxis()->SetTitle("counts");
        qa.fQAHistogramsList->Add(qa.fHistSelectedMultiplicity);
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
