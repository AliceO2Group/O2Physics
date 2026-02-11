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

/// \file MuPa-DataMembers.h
/// \brief ... TBI 20250425
/// \author Ante.Bilandzic@cern.ch

#ifndef PWGCF_MULTIPARTICLECORRELATIONS_CORE_MUPA_DATAMEMBERS_H_
#define PWGCF_MULTIPARTICLECORRELATIONS_CORE_MUPA_DATAMEMBERS_H_

#include <vector>

// General remarks:
// 0. Starting with C++11, it's possible to initialize data members at declaration, so I do it here
// 1. Use //!<! for introducing a Doxygen comment interpreted as transient in both ROOT 5 and ROOT 6.
// ...

// Categories:
// a) Base list to hold all output objects ("grandmother" of all lists);
// *) Task configuration;
// *) QA;
// *) Event histograms;
// *) Particle histograms;
// *) Q-vectors;
// *) Multiparticle correlations (standard, isotropic, same harmonic);
// *) Particle weights;
// *) Nested loops;
// *) Results;
// ...

// a) Base list to hold all output objects ("grandmother" of all lists):

TString sBaseListName = "Default list name"; // yes, I declare it separately, because I need it also later in BailOut() function
OutputObj<TList> fBaseList{sBaseListName.Data(),
                           OutputObjHandlingPolicy::AnalysisObject,
                           OutputObjSourceType::OutputObjSource};
TProfile* fBasePro = NULL;           //!<! keeps flags relevant for the whole analysis
TObjArray* fBaseProBinLabels = NULL; // helper for fBasePro to hold bin labels, until SetBinLabel(...) large memory consumption is resolved

// *) Task configuration:
struct TaskConfiguration {
  TString fTaskIsConfiguredFromJson = "no";      // the trick to ensure that settings from JSON are taken into account, even if only one configurable is misconfigured, when everything dies silently
  TString fTaskName = "";                        // task name - this one is used to get the right weights programatically for this analysis.
                                                 // If not set, weights are fetched from TDirectoryFile whose name ends with "multiparticle-correlations-a-b" (default)
                                                 // If set to "someName", weights are fetched from TDirectoryFile whose name ends with "multiparticle-correlations-a-b_someName"
                                                 // TBI 20250122 Therefore, when running in HL, it's important to configure manually cfTaskName to be exactly the same as subwagon name.
                                                 //              Can I automate this?
  TString fRunNumber = "";                       // over which run number this task is executed
  bool fRunNumberIsDetermined = false;           // ensures that run number is determined in process() and propagated to already booked objects only once
  int64_t fRunTime[eRunTime_N] = {0};            // stores permanently start of run, end of run, and run duration
  bool fDryRun = false;                          // book all histos and run without storing and calculating anything
  bool fVerbose = false;                         // print additional info during debugging, but not for simply utility function or function calls per particle (see next)
  bool fVerboseUtility = false;                  // print additional info during debugging also for simply utility function, but not for function calls per particle (see next)
  bool fVerboseForEachParticle = false;          // print additional info during debugging, also for function calls per particle
  bool fVerboseEventCounter = true;              // print or not only event counter
  bool fVerboseEventCut = true;                  // print or not only which event cut didn't survive
  bool fPlainPrintout = false;                   // print in color or in plain (use the latter in HL)
  bool fDoAdditionalInsanityChecks = false;      // do additional insanity checks at run time, at the expense of losing a bit of performance
                                                 // (for instance, check if the run number in the current 'collision' is the same as run number in the first 'collision', etc.)
  bool fInsanityCheckForEachParticle = false;    // do additional insanity checks at run time for each particle, at the expense of losing a lot of performance. Use only during debugging.
  bool fProcess[eProcess_N] = {false};           // set what to process. See enum eProcess for full description. Set via implicit variables within a PROCESS_SWITCH clause.
  TString fWhichProcess = "ProcessRec";          // dump in this variable which process was used
  unsigned int fRandomSeed = 0;                  // argument for TRandom3 constructor. By default it is 0 (seed is guaranteed to be unique in time and space)
  bool fUseFisherYates = false;                  // algorithm used to randomize particle indices, set via configurable
  TArrayI* fRandomIndices = NULL;                // array to store random indices obtained from Fisher-Yates algorithm
  int fFixedNumberOfRandomlySelectedTracks = -1; // use a fixed number of randomly selected particles in each event, applies to all centralities. It is set and applied if > 0. Set to <=0 to ignore.

  bool fUseStopwatch = false;                                 // do some basing profiling with TStopwatch for where the execution time is going
  TStopwatch* fTimer[eTimer_N] = {NULL};                      // stopwatch, global (overal execution time) and local
  float fFloatingPointPrecision = 1.e-6;                      // two floats are the same if abs(f1 - f2) < fFloatingPointPrecision (there is configurable for it)
  int fSequentialBailout = 0;                                 // if fSequentialBailout > 0, then each fSequentialBailout events the function BailOut() is called. Can be used for real analysis and for IV.
  bool fUseSpecificCuts = false;                              // apply after DefaultCuts() also hardwired analysis-specific cuts, determined via tc.fWhichSpecificCuts
  TString fWhichSpecificCuts = "";                            // determine which set of analysis-specific cuts will be applied after DefaultCuts(). Use in combination with tc.fUseSpecificCuts
  TString fSkipTheseRuns = "";                                // comma-separated list of runs which will be skipped during analysis in hl (a.k.a. "bad runs")
  bool fSkipRun = false;                                      // based on the content of fWhichSpecificCuts, skip or not the current run
  bool fCalculateAsFunctionOf[eAsFunctionOf_N] = {false};     //! [0=integrated,1=vs. multiplicity,2=vs. centrality,3=pT,4=eta,5=vs. occupancy, ...]
                                                              // Example: tc.fCalculateAsFunctionOf[AFO_PT] = mupa.fCalculateCorrelationsAsFunctionOf[AFO_PT] || t0.fCalculateTest0AsFunctionOf[AFO_PT]
                                                              //                                              || es.fCalculateEtaSeparationsAsFunctionOf[AFO_PT]
  bool fCalculate2DAsFunctionOf[eAsFunctionOf2D_N] = {false}; //! [0=integrated,1=vs. multiplicity,2=vs. centrality,3=pT,4=eta,5=vs. occupancy, ...]. See example above for 1D case
  bool fCalculate3DAsFunctionOf[eAsFunctionOf3D_N] = {false}; //! [0=integrated,1=vs. multiplicity,2=vs. centrality,3=pT,4=eta,5=vs. occupancy, ...]. See example above for 1D case
  TDatabasePDG* fDatabasePDG = NULL;                          // booked only when MC info is available. There is a standard memory blow-up when booked, therefore I need to request also fUseDatabasePDG = true
                                                              // TBI 20250625 replace eventually with the service O2DatabasePDG, when memory consumption problem is resolved
  bool fUseSetBinLabel = false;                               // until SetBinLabel(...) large memory consumption is resolved, do not use hist->SetBinLabel(...), see ROOT Forum
                                                              // See also local executable PostprocessLabels.C
  bool fUseClone = false;                                     // until Clone(...) large memory consumption is resolved, do not use hist->Clone(...), see ROOT Forum
  bool fUseFormula = false;                                   // until TFormula large memory consumption is resolved, do not use, see ROOT Forum
  bool fUseDatabasePDG = false;                               // I use it at the moment only to retreive charge for MC particle from its PDG code, because there is no direct getter mcParticle.sign()
                                                              // But most likely I will use it to retrieve other particle proprties from PDG table. There is a standard memory blow-up when used.
} tc;                                                         // "tc" labels an instance of this group of variables.

// *) Event-by-event quantities:
struct EventByEventQuantities {
  int fSelectedTracks = 0;            // integer counter of tracks used to calculate Q-vectors, after all particle cuts have been applied
  float fMultiplicity = 0.;           // my internal multiplicity, can be set to fSelectedTracks (calculated internally), fReferenceMultiplicity (calculated outside of my code), etc.
                                      // Results "vs. mult" are plotted against fMultiplicity, whatever it is set to.
                                      // Use configurable cfMultiplicityEstimator[eMultiplicityEstimator] to define what is this multiplicity, by default it is "SelectedTracks"
  float fReferenceMultiplicity = 0.;  // reference multiplicity, calculated outside of my code. Can be "MultTPC", "MultFV0M", etc.
                                      // Use configurable cfReferenceMultiplicityEstimator[eReferenceMultiplicityEstimator]" to define what is this multiplicity,
                                      // by default it is "TBI 20241123 I do not know yet which estimator is best for ref. mult."
  float fCentrality = 0.;             // event-by-event centrality, in reconstructed data. Value of the default centrality estimator, set via configurable cfCentralityEstimator
  float fCentralitySim = 0.;          // event-by-event centrality, in simulated data. Calculated directly from IP at the moment, eventually I will access it from o2::aod::hepmcheavyion::Centrality
  float fOccupancy = 0.;              // event-by-event occupancy. Value of the default occupancy estimator, set via configurable cfOccupancyEstimator.
                                      // Remebmer that collision with occupanct 0. shall NOT be rejected, therefore in configurable I set -0.0001 for low edge by default.
  float fInteractionRate = 0.;        // event-by-event interaction rate
  float fCurrentRunDuration = 0.;     // how many seconds after start of run this collision was taken, i.e. seconds after start of run (SOR)
  float fVz = 0.;                     // vertex z position
  float fVzSim = 0.;                  // vertex z position, in simulated data
  float fFT0CAmplitudeOnFoundBC = 0.; // TBI20250331 finalize the comment here
  float fImpactParameter = 0.;        // calculated only for simulated/generated data
} ebye;                               // "ebye" is a common label for objects in this struct

// *) Particle-by-particle quantities:
//    Remark: Here I define all particle quantities, that I need across several member functions.
struct ParticleByParticleQuantities {
  double fPhi = 0.;      // azimuthal angle
  double fPt = 0.;       // transverse momentum
  double fEta = 0.;      // pseudorapidity
  double fCharge = -44.; // particle charge. Yes, never initialize charge to 0.
} pbyp;

// *) QA:
//    Remark 1: I keep new histograms in this group, until I need them permanently in the analysis. Then, they are moved to EventHistograms or ParticleHistograms (yes, even if they are 2D).
//    Remark 2: All 2D histograms book as TH2F, due to "stmem error" in terminate (see .cxx for further details)
struct QualityAssurance {
  TList* fQAList = NULL;                   //!<! base list to hold all QA output object
  TProfile* fQAHistogramsPro = NULL;       //!<! keeps flags relevant for the QA histograms
  bool fCheckUnderflowAndOverflow = false; // check and bail out if in event and particle histograms there are entries which went to underflow or overflow bins
  int fRebin = 1;                          // number of bins of selected heavy 2D histograms are devided with this number, there is a configurable cfRebin
  // ...

  TList* fQAEventList = NULL;                                            //!<! base list to hold all QA event output object
  TH2F* fQAEventHistograms2D[eQAEventHistograms2D_N][2][2] = {{{NULL}}}; //! [ type - see enum eQAEventHistograms2D ][reco,sim][before, after particle cuts]
  bool fFillQAEventHistograms2D = true;                                  // if false, all 2D event histograms are not filled. if true, the ones for which fBookQAEventHistograms2D[...] is true, are filled
  bool fBookQAEventHistograms2D[eQAEventHistograms2D_N] = {true};        // book or not this 2D histogram, see configurable cfBookQAEventHistograms2D
  float fEventHistogramsBins2D[eQAEventHistograms2D_N][2][3] = {{{0.}}}; // [type - see enum][x,y][nBins,min,max]
  TString fEventHistogramsName2D[eQAEventHistograms2D_N] = {""};         // name of fQAEventHistograms2D, determined programatically from other 1D names, to ease bookkeeping

  TList* fQAParticleList = NULL;                                               //!<! base list to hold all QA particle output object
  TH2F* fQAParticleHistograms2D[eQAParticleHistograms2D_N][2][2] = {{{NULL}}}; //! [ type - see enum eQAParticleHistograms2D ][reco,sim][before, after particle cuts]
  bool fFillQAParticleHistograms2D = true;                                     // if false, all 2D histograms in this category are not filled. If true, the ones for which fBookQAParticleHistograms2D[...] is true, are filled
  bool fBookQAParticleHistograms2D[eQAParticleHistograms2D_N] = {true};        // book or not this 2D histogram, see configurable cfBookQAParticleHistograms2D
  float fParticleHistogramsBins2D[eQAParticleHistograms2D_N][2][3] = {{{0.}}}; // [type - see enum][x,y][nBins,min,max]
  TString fParticleHistogramsName2D[eQAParticleHistograms2D_N] = {""};         // name of fQAParticleHistograms2D, determined programatically from other 1D names, to ease bookkeeping

  TList* fQAParticleEventList = NULL;                                                      //!<! base list to hold all QA particle event output object
  TH2F* fQAParticleEventHistograms2D[eQAParticleEventHistograms2D_N][2][2] = {{{NULL}}};   //! [ type - see enum eQAParticleEventHistograms2D ][reco,sim][before, after cuts]
  bool fFillQAParticleEventHistograms2D = true;                                            // if false, all 2D histograms in this category are not filled. If true, the ones for which fBookQAParticleEventHistograms2D[...] is true, are filled
  bool fBookQAParticleEventHistograms2D[eQAParticleEventHistograms2D_N] = {true};          // book or not this 2D histogram, see configurable cfBookQAParticleEventHistograms2D
  float fQAParticleEventHistogramsBins2D[eQAParticleEventHistograms2D_N][2][3] = {{{0.}}}; // [type - see enum][x,y][nBins,min,max]
  TString fQAParticleEventHistogramsName2D[eQAParticleEventHistograms2D_N] = {""};         // name of fQAParticleEventHistograms2D, determined programatically from other 1D names, to ease bookkeeping
  TProfile* fQAParticleEventProEbyE[2][2] = {{NULL}};                                      // helper profile to calculate <some-particle-property> event-by-event
                                                                                           // [reco, sim][before, after]. Type dimension is bin.

  TList* fQACorrelationsVsList = NULL;                                                                //!<! base list to hold all QA "CorrelationsVs" output object
  TH2F* fQACorrelationsVsHistograms2D[eQACorrelationsVsHistograms2D_N][gMaxHarmonic][2] = {{{NULL}}}; //! [ type - see enum eQACorrelationsVsHistograms2D ][reco,sim]. I do not have here support for [before, after], because I do not fill Q-vectors before cuts
  bool fFillQACorrelationsVsHistograms2D = true;                                                      // if false, all 2D histograms in this category are not filled. If true, the ones for which fBookQACorrelationsVsHistograms2D[...] is true, are filled
  bool fBookQACorrelationsVsHistograms2D[eQACorrelationsVsHistograms2D_N] = {true};                   // book or not this 2D histogram, see configurable cfBookQACorrelationsVsHistograms2D
  float fQACorrelationsVsHistogramsBins2D[eQACorrelationsVsHistograms2D_N][2][3] = {{{0.}}};          // [type - see enum][x,y][nBins,min,max]
  TString fQACorrelationsVsHistogramsName2D[eQACorrelationsVsHistograms2D_N] = {""};                  // name of fQACorrelationsVsHistograms2D, determined programatically from other 1D names, to ease bookkeeping
  int fQACorrelationsVsHistogramsMinMaxHarmonic[2];                                                   // book only for MinMaxHarmonic[0] <= harmonics < MinMaxHarmonic[1]

  TList* fQACorrelationsVsInteractionRateVsList = NULL;                                                                    //!<! base list to hold all QA "CorrelationsVsInteractionRateVs" output object
  TProfile2D* fQACorrVsIRVsProfiles2D[eQACorrelationsVsInteractionRateVsProfiles2D_N][gMaxHarmonic][2] = {{{NULL}}};       //! [ type - see enum eQACorrelationsVsInteractionRateVsProfiles2D_N ][reco,sim]. I do not have here support for [before, after], because I do not fill Q-vectors before cuts
  bool fFillQACorrelationsVsInteractionRateVsProfiles2D = true;                                                            // if false, all 2D profiles in this category are not filled. If true, the ones for which fBookQACorrelationsVsInteractionRateVsProfiles2D[...] is true, are filled
  bool fBookQACorrelationsVsInteractionRateVsProfiles2D[eQACorrelationsVsInteractionRateVsProfiles2D_N] = {true};          // book or not this 2D profile, see configurable cfBookQACorrelationsVsInteractionRateVsProfiles2D
  float fQACorrelationsVsInteractionRateVsProfilesBins2D[eQACorrelationsVsInteractionRateVsProfiles2D_N][2][3] = {{{0.}}}; // [type - see enum][x,y][nBins,min,max]
  TString fQACorrelationsVsInteractionRateVsProfilesName2D[eQACorrelationsVsInteractionRateVsProfiles2D_N] = {""};         // name of fQACorrelationsVsInteractionRateVsProfiles2D, determined programatically from other 1D names, to ease bookkeeping
  int fQACorrelationsVsInteractionRateVsProfilesMinMaxHarmonic[2];                                                         // book only for MinMaxHarmonic[0] <= harmonics < MinMaxHarmonic[1]

  float fReferenceMultiplicity[eReferenceMultiplicityEstimators_N] = {0.};                // used mostly in QA correlation plots
  TString fReferenceMultiplicityEstimatorName[eReferenceMultiplicityEstimators_N] = {""}; // TBI 20241123 add comment
  float fCentrality[eCentralityEstimators_N] = {0.};                                      // used mostly in QA correlation plots
  TString fCentralityEstimatorName[eCentralityEstimators_N] = {""};                       // TBI 20241123 add comment
  float fOccupancy[eOccupancyEstimators_N] = {0.};                                        // used mostly in QA correlation plots
  TString fOccupancyEstimatorName[eOccupancyEstimators_N] = {""};                         // TBI 20241123 add comment

} qa; // "qa" is a common label for objects in this struct

// *) Event histograms:
struct EventHistograms {
  TList* fEventHistogramsList = NULL;   //!<! list to hold all control event histograms
  TProfile* fEventHistogramsPro = NULL; //!<! keeps flags relevant for the control event histograms
  // 1D:
  TH1F* fEventHistograms[eEventHistograms_N][2][2] = {{{NULL}}}; //! [ type - see enum eEventHistograms ][reco,sim][before, after event cuts]
  bool fFillEventHistograms = true;                              // if false, all event histograms are not filled. if true, the ones for which fBookEventHistograms[...] is true, are filled
  bool fBookEventHistograms[eEventHistograms_N] = {true};        // book or not this histogram, see SetBookEventHistograms
  float fEventHistogramsBins[eEventHistograms_N][3] = {{0.}};    // [nBins,min,max]
  TString fEventHistogramsName[eEventHistograms_N] = {""};       // name of event histogram, used both for 1D and 2D histograms
  int fEventCounter[eEventCounter_N] = {0};                      // event counters, see enum eEventCounter for full explanation
  // 2D:
  // ...
  // Remark: All 2D event histograms are still in the QA group. Move here only the ones I will use regularly in the analysis
} eh; // "eh" labels an instance of group of histograms "EventHistograms"

// *) Event cuts:
struct EventCuts {
  TList* fEventCutsList = NULL;                                //!<! list to hold all event cuts objects
  TProfile* fEventCutsPro = NULL;                              //!<! keeps flags relevant for the event cuts
  bool fUseEventCuts[eEventCuts_N] = {false};                  // Use or do not use a cut enumerated in eEventHistograms + eEventCuts
  bool fUseEventCutCounterAbsolute = false;                    // profile and save how many times each event cut counter triggered (absolute). Use with care, as this is computationally heavy
  bool fUseEventCutCounterSequential = false;                  // profile and save how many times each event cut counter triggered (sequential). Use with care, as this is computationally heavy
  bool fEventCutCounterBinLabelingIsDone = false;              // this flag ensures that ordered labeling of bins, to resemble ordering of cut implementation, is done only once
  bool fPrintCutCounterContent = false;                        // if true, prints on the screen content of fEventCutCounterHist[][] (all which were booked)
  TString fEventCutName[eEventCuts_N] = {""};                  // event cut name, with default ordering defined by ordering in enum eEventCuts
  TExMap* fEventCutCounterMap[2] = {NULL};                     // map (key, value) = (enum eEventCuts, ordered bin number)
  TExMap* fEventCutCounterMapInverse[2] = {NULL};              // inverse of above fEventCutCounterMap, i.e. (ordered bin number, enum eEventCuts)
  int fEventCutCounterBinNumber[2] = {1, 1};                   // bin counter for set bin labels in fEventCutCounterHist
  float fdEventCuts[eEventCuts_N][2] = {{0.}};                 // event cuts defined via [min,max)
  TString fsEventCuts[eEventCuts_N] = {""};                    // event cuts defined via string
  TH1F* fEventCutCounterHist[2][eCutCounter_N] = {{NULL}};     //!<! [rec,sim][see enum eCutCounter] histogram to store how many times each event cut triggered
  int fBeforeAfterColor[2] = {kRed, kGreen};                   // color code before and after cuts
  float fCentralityCorrelationsCutTreshold = 5.;               // see bool CentralityCorrelationCut()
  TString fCentralityCorrelationsCutVersion = "Absolute";      // see bool CentralityCorrelationCut()
  float fCentralityValues[2] = {0.};                           // [0] value of first cent. estimator, [1] = value of second cent. estimator, when CentralityCorrelationsCut is requested
  TFormula* fEventCutsFormulas[eEventCutsFormulas_N] = {NULL}; // see enum TBI 20250415 I do not use it for the time being, due to memory blow-up in TFormula
  float fdEventCutsFormulas[eEventCutsFormulas_N][2] = {{0.}}; // I need this only temporarily until large memory consumption with TFormula is resolved.
                                                               // I support at the moment only linear cut in the format p0 + p1*x. Then, [0] = p0, [1] = p1.
} ec;                                                          // "ec" is a common label for objects in this struct

// *) Particle histograms:
struct ParticleHistograms {
  TList* fParticleHistogramsList = NULL;   //!<! list to hold all control particle histograms
  TProfile* fParticleHistogramsPro = NULL; //!<! keeps flags relevant for the control particle histograms
  // 1D:
  TH1F* fParticleHistograms[eParticleHistograms_N][2][2] = {{{NULL}}}; //! [ type - see enum eParticleHistograms ][reco,sim][before, after particle cuts]
  bool fFillParticleHistograms = true;                                 // if false, all 1D particle histograms are not filled.
                                                                       // if true, the ones for which fBookParticleHistograms[...] is true, are filled
  bool fBookParticleHistograms[eParticleHistograms_N] = {true};        // book or not the particular particle histogram, see configurable cfBookParticleHistograms
  float fParticleHistogramsBins[eParticleHistograms_N][3] = {{0.}};    // [nBins,min,max]
  TString fParticleHistogramsName[eParticleHistograms_N] = {""};       // name of particle histogram, used both for 1D and 2D histograms
  // 2D:
  TH2D* fParticleHistograms2D[eParticleHistograms2D_N][2][2] = {{{NULL}}};   //! [ type - see enum eParticleHistograms2D ][reco,sim][before, after particle cuts]
  bool fFillParticleHistograms2D = true;                                     // if false, all 2D particle histograms are not filled.
                                                                             // if true, the ones for which fBookParticleHistograms2D[...] is true, are filled
  bool fBookParticleHistograms2D[eParticleHistograms2D_N] = {true};          // book or not this 2D histogram, see configurable cfBookParticleHistograms2D
  float fParticleHistogramsBins2D[eParticleHistograms2D_N][2][3] = {{{0.}}}; // [type - see enum][x,y][nBins,min,max]
  TString fParticleHistogramsName2D[eParticleHistograms2D_N] = {""};         // name of particle histogram 2D, determined programatically from two 1D, in the format "%s_vs_%s"

  // **) n-dimensional sparse histograms:
  THnSparse* fParticleSparseHistograms[eDiffWeightCategory_N][2][2] = {{{NULL}}}; //! [ category of sparse histograms - see enum eDiffWeightCategory ][reco,sim][before, after particle cuts]
                                                                                  // Remark 0: I anticipate I will need this only for differential particle weights,
                                                                                  //           therefore I couple it with eDiffWeightCategory_N
                                                                                  // Remark 1: I fill these histograms only AFTER cuts, therefore no need for extra dimension
  bool fBookParticleSparseHistograms[eDiffWeightCategory_N] = {false};            // fill or not specific category of sparse histograms

  bool fFillParticleSparseHistogramsBeforeCuts = false; // by default, I fill sparse histograms only after the cuts. In rare cases, e.g. in internal validation
                                                        // when I am developing pT and eta weights, I calculate them from the ratio [sim][before] / [sim][after],
                                                        // therefore in that case I need to fill sparse also before cuts. As of 20251124, this is the only case when it's justified
                                                        // to fill sparse also before cuts

  // bool fFillParticleSparseHistogramsDimension[eDiffWeightCategory_N][gMaxNumberSparseDimensions] = {{true}}; // fill or not the specific dimension of a category of sparse histograms TBI 20250223 implement this eventually
  TString fParticleSparseHistogramsName[eDiffWeightCategory_N] = {""};                                      // name of particle sparse histogram, determined programatically from requested axes
  TString fParticleSparseHistogramsTitle[eDiffWeightCategory_N] = {""};                                     // title of particle sparse histogram, determined programatically from requested axes
  int fParticleSparseHistogramsNBins[eDiffWeightCategory_N][gMaxNumberSparseDimensions] = {{0}};            // number of bins. I do not have min and max, because for sparse I use BinEdges, see below
  TArrayD* fParticleSparseHistogramsBinEdges[eDiffWeightCategory_N][gMaxNumberSparseDimensions] = {{NULL}}; // arrays holding bin edges, see the usage of SetBinEdges for sparse histograms
  TString fParticleSparseHistogramsAxisTitle[eDiffWeightCategory_N][gMaxNumberSparseDimensions] = {{""}};   // axis title
  int fRebinSparse = 1;                                                                                     // used only for all fixed-length bins which are implemented directly for sparse histograms (i.e. not inherited from results histograms)
} ph;                                                                                                       // "ph" labels an instance of group of histograms "ParticleHistograms"

// *) Particle cuts:
struct ParticleCuts {
  TList* fParticleCutsList = NULL;                            //!<! list to hold all particle cuts objects
  TProfile* fParticleCutsPro = NULL;                          //!<! keeps flags relevant for the particle cuts
  bool fUseParticleCuts[eParticleCuts_N] = {false};           // true or false .
  bool fUseParticleCutCounterAbsolute = false;                // profile and save how many times each particle cut counter triggered (absolute). Use with care, as this is computationally heavy
  bool fUseParticleCutCounterSequential = false;              // profile and save how many times each particle cut counter triggered (sequential). Use with care, as this is computationally heavy
  bool fParticleCutCounterBinLabelingIsDone = false;          // this flag ensures that ordered labeling of bins, to resemble ordering of cut implementation, is done only once.
  TString fParticleCutName[eParticleCuts_N] = {""};           // particle cut name, as used in bin labels, etc.
  TExMap* fParticleCutCounterMap[2] = {NULL};                 // map (key, value) = (enum eParticleCuts, ordered bin number)
  TExMap* fParticleCutCounterMapInverse[2] = {NULL};          // inverse of above fParticleCutCounterMap, i.e. (ordered bin number, enum eParticleCuts)
  int fParticleCutCounterBinNumber[2] = {1, 1};               // bin counter for set bin labels in fParticleCutCounterHist
  float fdParticleCuts[eParticleCuts_N][2] = {{0.}};          // particles cuts defined via [min,max) . Remark: I use here eParticleHistograms_N , not to duplicate these enums for ParticleCuts.
  TString fsParticleCuts[eParticleCuts_N] = {""};             // particles cuts defined via booleans via string
  TH1F* fParticleCutCounterHist[2][eCutCounter_N] = {{NULL}}; //!<! [rec,sim][see enum eCutCounter] histogram to store how many any times each particle cut triggered
  TFormula* fPtDependentDCAxyFormula = NULL;                  // the actual formula, used to evaluate for a given pT, the corresponding DCAxy, where the parameterization is given by configurable cfPtDependentDCAxyParameterization
} pc;                                                         // "pc" is a common label for objects in this struct

// *) Q-vectors:
struct Qvector {
  TList* fQvectorList = NULL;                                                                // list to hold all Q-vector objects
  TProfile* fQvectorFlagsPro = NULL;                                                         // profile to hold all flags for Q-vector
  bool fCalculateQvectors = true;                                                            // to calculate or not to calculate Q-vectors, that's a Boolean...
                                                                                             // Does NOT apply to Qa, Qb, etc., vectors, needed for eta separ.
  TComplex fQ[gMaxHarmonic * gMaxCorrelator + 1][gMaxCorrelator + 1] = {{TComplex(0., 0.)}}; //! generic Q-vector, legacy code (TBI 20250718 remove, and switch to line below eventually)
  // std::vector<std::vector<std::complex<double>>> fQ; // generic Q-vector
  TComplex fQvector[gMaxHarmonic * gMaxCorrelator + 1][gMaxCorrelator + 1] = {{TComplex(0., 0.)}}; //! integrated Q-vector, legacy code (TBI 20250718 remove, and switch to line below eventually)
  // std::vector<std::vector<std::complex<double>>> fQvector; // dynamically allocated integrated Q-vector => it has to be done this way, to optimize memory usage

  bool fCalculateqvectorsKineAny = false;                              // by default, it's off. It's set to true automatically if any of kine correlators is requested,
                                                                       // either for Correlations, Test0, EtaSeparations, etc.
  bool fCalculateqvectorsKine[eqvectorKine_N] = {false};               // same as above, just specifically for each enum eqvectorKine + applies only to Correlations and Test0
  bool fCalculateqvectorsKineEtaSeparations[eqvectorKine_N] = {false}; // same as above, just specifically for each enum eqvectorKine + applies only to EtaSeparations

  std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>> fqvector; // dynamically allocated differential q-vector => it has to be done this way, to optimize memory usage
                                                                                     // dimensions: [eqvectorKine_N][gMaxNoBinsKine][gMaxHarmonic * gMaxCorrelator + 1][gMaxCorrelator + 1]
  std::vector<int> fNumberOfKineBins = {0};                                          // for each kine vector which was requested in this analysis, here I calculate and store the corresponding number of kine bins
  std::vector<std::vector<int>> fqvectorEntries;                                     // dynamically allocated number of entries for differential q-vector => it has to be done this way, to optimize memory usage. Dimensions: [eqvectorKine_N][gMaxNoBinsKine]

  // q-vectors for eta separations:
  TComplex fQabVector[2][gMaxHarmonic][gMaxNumberEtaSeparations] = {{{TComplex(0., 0.)}}};          //! integrated [-eta or +eta][harmonic][eta separation]
  float fMab[2][gMaxNumberEtaSeparations] = {{0.}};                                                 //! multiplicities in 2 eta separated intervals
  TH1F* fMabDist[2][2][2][gMaxNumberEtaSeparations] = {{{{NULL}}}};                                 // multiplicity distributions in A and B, for each eta separation [ A or B ] [rec or sim] [ before or after cuts ] [ eta separation value ]
  std::vector<std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>>> fqabVector; // dynamically allocated differential q-vector.
                                                                                                    // dimensions: [-eta or +eta][eqvectorKine_N][global binNo][harmonic][eta separation]
                                                                                                    // Remark: Unlike fqvector above, here I support only 2-p correlations,
                                                                                                    // therefore no need for "[gMaxHarmonic * gMaxCorrelator + 1][gMaxCorrelator + 1]", etc.
  std::vector<std::vector<std::vector<std::vector<float>>>> fmab;                                   //! multiplicities vs kine in 2 eta separated intervals
                                                                                                    // [-eta or +eta][eqvectorKine_N][global binNo][eta separation]
} qv;                                                                                               // "qv" is a common label for objects in this struct

// *) Multiparticle correlations (standard, isotropic, same harmonic):
struct MultiparticleCorrelations {
  TList* fCorrelationsList = NULL;                                           // list to hold all correlations objects
  TProfile* fCorrelationsFlagsPro = NULL;                                    // profile to hold all flags for correlations
  bool fCalculateCorrelations = false;                                       // calculate and store integrated correlations
  TProfile* fCorrelationsPro[4][gMaxHarmonic][eAsFunctionOf_N] = {{{NULL}}}; //! multiparticle correlations
                                                                             //  [2p=0,4p=1,6p=2,8p=3][n=1,n=2,...,n=gMaxHarmonic]
                                                                             //  [0=integrated,1=vs. multiplicity,2=vs. centrality,3=pT,4=eta,5=vs. occupancy, ...]
  bool fCalculateCorrelationsAsFunctionOf[eAsFunctionOf_N] = {false};        //! [0=integrated,1=vs. multiplicity,2=vs. centrality,3=pT,4=eta,5=vs. occupancy, ...]
                                                                             //  As of 20241111, 3=pT and 4=eta are not implemented, see void CalculateKineCorrelations(...)
} mupa;                                                                      // "mupa" is a common label for objects in this struct

// *) Particle weights:
struct ParticleWeights {
  TList* fWeightsList = NULL;                                             //!<! list to hold all particle weights
  TProfile* fWeightsFlagsPro = NULL;                                      //!<! profile to hold all flags for weights
  bool fUseWeights[eWeights_N] = {false};                                 // use weights [phi,pt,eta]
  TH1D* fWeightsHist[eWeights_N] = {NULL};                                //!<! particle weights
  bool fUseDiffWeights[eDiffWeights_N] = {false};                         // use differential weights [phipt,phieta] => TBI 20250215 this is obsolete and superseeded with fUseDiffPhiWeights, etc.
  TH1D* fDiffWeightsHist[eDiffWeights_N][gMaxBinsDiffWeights] = {{NULL}}; // histograms holding differential weights [phipt,phieta][bin number] => TBI 20250222 obsolete

  // ** sparse histograms:
  THnSparse* fDiffWeightsSparse[eDiffWeightCategory_N] = {NULL}; // multidimensional sparse histogram to hold all differential phi-weights (as a function of pt, eta, etc.).
                                                                 // each dimension has its own enum category, e.g. 0 = eDWPhi => eDiffPhiWeights, 1 = eDWPt => eDiffPtWeights, etc.
  bool fUseDiffPhiWeights[eDiffPhiWeights_N] = {false};          // use differential phi weights, see enum eDiffPhiWeights for supported dimensions
  bool fUseDiffPtWeights[eDiffPtWeights_N] = {false};            // use differential pt weights, see enum eDiffPtWeights for supported dimensions
  bool fUseDiffEtaWeights[eDiffEtaWeights_N] = {false};          // use differential eta weights, see enum eDiffEtaWeights for supported dimensions
  // ...
  int fDWdimension[eDiffWeightCategory_N] = {0};           // dimension of differential weight for each category in current analysis
  TArrayD* fFindBinVector[eDiffWeightCategory_N] = {NULL}; // this is the vector I use to find bin when I obtain weights with sparse histograms

  TString fFileWithWeights = "";           // path to external ROOT file which holds all particle weights
  bool fParticleWeightsAreFetched = false; // ensures that particle weights are fetched only once
} pw;                                      // "pw" labels an instance of this group of histograms

// *) Centrality weights:
struct CentralityWeights {
  TList* fCentralityWeightsList = NULL;        // list to hold all Q-vector objects
  TProfile* fCentralityWeightsFlagsPro = NULL; // profile to hold all flags for CentralityWeights
  bool fUseCentralityWeights = false;          // use centrality weights
  TH1D* fCentralityWeightsHist = NULL;         // histograms holding centrality weights
  TString fFileWithCentralityWeights = "";     // path to external ROOT file which holds all centrality weights
  bool fCentralityWeightsAreFetched = false;   // ensures that centrality weights are fetched only once
} cw;

// *) Nested loops:
struct NestedLoops {
  TList* fNestedLoopsList = NULL;                                              // list to hold all nested loops objects
  TProfile* fNestedLoopsFlagsPro = NULL;                                       // profile to hold all flags for nested loops
  bool fCalculateNestedLoops = false;                                          // calculate and store correlations with nested loops, as a cross-check
  bool fCalculateCustomNestedLoops = false;                                    // validate e-b-e all correlations with custom nested loop
  bool fCalculateKineCustomNestedLoops = false;                                // validate e-b-e all differential (vs pt, eta, etc.) correlations with custom nested loop
  int fMaxNestedLoop = -1;                                                     // if set to e.g. 4, all nested loops beyond that, e.g. 6-p and 8-p, are NOT calculated
  TProfile* fNestedLoopsPro[4][gMaxHarmonic][eAsFunctionOf_N] = {{{NULL}}};    //! multiparticle correlations from nested loops
                                                                               //! [2p=0,4p=1,6p=2,8p=3][n=1,n=2,...,n=gMaxHarmonic][0=integrated,1=vs.
                                                                               //! multiplicity,2=vs. centrality,3=pT,4=eta]
  TArrayD* ftaNestedLoops[2] = {NULL};                                         //! e-b-e container for nested loops [0=angles;1=product of all weights]
  TArrayD* ftaNestedLoopsKine[eqvectorKine_N][gMaxNoBinsKine][2] = {{{NULL}}}; //! e-b-e container for nested loops // [0=pT,1=eta,2=...][kine bin][0=angles;1=product of all weights]
} nl;                                                                          // "nl" labels an instance of this group of histograms

// *) Toy NUA (can be applied both in real data analysis and in analysis 'on-the-fly', e.g. when running internal validation):
struct NUA {
  TList* fNUAList = NULL;                                 // list to hold all NUA objects
  TProfile* fNUAFlagsPro = NULL;                          // profile to hold all flags for NUA objects
  bool fApplyNUAPDF[eNUAPDF_N] = {false};                 // apply NUA to particular kine variable (see the corresponding enum eNUAPDF)
  bool fUseDefaultNUAPDF[eNUAPDF_N] = {true, true, true}; // by default, use simple hardcoded expressions for NUA acceptance profile
  TF1* fDefaultNUAPDF[eNUAPDF_N] = {NULL};                // default distributions used as pdfs to simulate NUA on-the-fly
  TH1D* fCustomNUAPDF[eNUAPDF_N] = {NULL};                // custom, user-supplied distributions used to simulate NUA
  TString* fCustomNUAPDFHistNames[eNUAPDF_N] = {NULL};    // these are the names of histograms holding custom NUA in an external file. There is a configurable for this one.
  TString fFileWithCustomNUA = "";                        // path to external ROOT file which holds all histograms with custom NUA
  float fMaxValuePDF[eNUAPDF_N] = {0.};                   // see algorithm used in Accept(...). I implemented it as a data member, so that it is not calculated again and again at each particle call
} nua;

// *) Internal validation:
struct InternalValidation {
  TList* fInternalValidationList = NULL;              // list to hold all objects for internal validation
  TProfile* fInternalValidationFlagsPro = NULL;       // profile to hold all flags for internal validation
  bool fUseInternalValidation = false;                // use internal validation
  bool fInternalValidationForceBailout = false;       // force bailout in internal validation after either eNumberOfEvents or eSelectedEvents is reached.
                                                      // This is OK as long as I do not apply any event cuts in InternalValidation().
                                                      // Remember that for each real event, I do fnEventsInternalValidation events on-the-fly.
                                                      // Can be used in combination with setting fSequentialBailout > 0.
  unsigned int fnEventsInternalValidation = 0;        // how many on-the-fly events will be sampled for each real event, for internal validation
  TString* fHarmonicsOptionInternalValidation = NULL; // "constant", "correlated", "persistent", "ptDependent", "ptEtaDependent", see .cxx for full documentation
  bool fRescaleWithTheoreticalInput = false;          // if true, all measured correlators are rescaled with theoretical input, so that in profiles everything is at 1
  bool fRandomizeReactionPlane = true;                // if true, RP is randomized e-by-e. I need false basically only when validating against theoretical input non-isotropic correlators
  TArrayD* fInternalValidationVnPsin[2] = {NULL};     // 0 = { v1, v2, ... }, 1 = { Psi1, Psi2, ... }
  int fMultRangeInternalValidation[2] = {0, 0};       // min and max values for uniform multiplicity distribution in on-the-fly analysis (convention: min <= M < max)
} iv;

// *) Test0:
struct Test0 {
  TList* fTest0List = NULL;                                                           // list to hold all objects for Test0
  TProfile* fTest0FlagsPro = NULL;                                                    // store all flags for Test0
  bool fCalculateTest0 = false;                                                       // calculate or not Test0
  TProfile* fTest0Pro[gMaxCorrelator][gMaxIndex][eAsFunctionOf_N] = {{{NULL}}};       //! [order][index][0=integrated,1=vs. multiplicity,2=vs. centrality,3=pT,4=eta]
  bool fCalculate2DTest0 = false;                                                     // calculate or not 2D Test0
  TProfile2D* fTest0Pro2D[gMaxCorrelator][gMaxIndex][eAsFunctionOf2D_N] = {{{NULL}}}; //! [order][index][0=cent vs pt, ..., see enum eAsFunctionOf2D]
  bool fCalculate3DTest0 = false;                                                     // calculate or not 2D Test0
  TProfile3D* fTest0Pro3D[gMaxCorrelator][gMaxIndex][eAsFunctionOf3D_N] = {{{NULL}}}; //! [order][index][0=cent vs pt vs eta, ..., see enum eAsFunctionOf3D]
  TString* fTest0Labels[gMaxCorrelator][gMaxIndex] = {{NULL}};                        // all labels: k-p'th order is stored in k-1'th index. So yes, I also store 1-p
  bool fCalculateTest0AsFunctionOf[eAsFunctionOf_N] = {false};                        //! [0=integrated,1=vs. multiplicity,2=vs. centrality,3=pT,4=eta,5=vs. occupancy, ...]
  bool fCalculate2DTest0AsFunctionOf[eAsFunctionOf2D_N] = {false};                    //! [0=integrated,1=vs. multiplicity,2=vs. centrality,3=pT,4=eta,5=vs. occupancy, ...]
  bool fCalculate3DTest0AsFunctionOf[eAsFunctionOf3D_N] = {false};                    //! [0=integrated,1=vs. multiplicity,2=vs. centrality,3=pT,4=eta,5=vs. occupancy, ...]
  TString fFileWithLabels = "";                                                       // path to external ROOT file which specifies all labels of interest
  bool fUseDefaultLabels = false;                                                     // use default labels hardwired in GetDefaultObjArrayWithLabels(), the choice is made with cfWhichDefaultLabels
  TString fWhichDefaultLabels = "";                                                   // only for testing purposes, select one set of default labels, see GetDefaultObjArrayWithLabels for supported options
} t0;                                                                                 // "t0" labels an instance of this group of histograms

// *) Eta separations:
struct EtaSeparations {
  TList* fEtaSeparationsList;                                                            // list to hold all correlations with eta separations
  TProfile* fEtaSeparationsFlagsPro;                                                     // profile to hold all flags for correlations with eta separations
  bool fCalculateEtaSeparations;                                                         // calculate correlations with eta separations
  bool fCalculateEtaSeparationsAsFunctionOf[eAsFunctionOf_N] = {false};                  //! [0=integrated,1=vs. multiplicity,2=vs. centrality,3=pT,4=eta,5=vs. occupancy, ...]
  float fEtaSeparationsValues[gMaxNumberEtaSeparations] = {-1.};                         // this array holds eta separation interals for which 2p correlations with eta separation will be calculated
                                                                                         // See the corresponding cofigurable cfEtaSeparationsValues. If entry is -1, it's ignored
  bool fEtaSeparationsSkipHarmonics[gMaxHarmonic] = {false};                             // For calculation of 2p correlation with eta separation these harmonics will be skipped
  TProfile* fEtaSeparationsPro[gMaxHarmonic][gMaxNumberEtaSeparations][eAsFunctionOf_N]; // [harmonic, 0 = v1, 8 = v9][ different eta Separations - see that enum ] [ AFO ]
} es;

// *) Global cosmetics:
struct GlobalCosmetics {
  TString srs[2] = {"rec", "sim"};                             // used in the histogram name as index when saved to the file
  TString srsLong[2] = {"reconstructed", "simulated"};         // used in the histogram title
  TString sba[2] = {"before", "after"};                        // used in the histogram name as index when saved to the file
  TString sbaLong[2] = {"before cuts", "after cuts"};          // used in the histogram title
  TString scc[eCutCounter_N] = {"abs", "seq"};                 // used in the histogram name as index when saved to the file
  TString sccLong[eCutCounter_N] = {"absolute", "sequential"}; // used in the histogram title
} gc;

// *) Results:
struct Results {                                         // This is in addition also sort of "abstract" interface, which defines common binning, etc., for other groups of histograms.
  TList* fResultsList = NULL;                            //!<! list to hold all results
  TProfile* fResultsFlagsPro = NULL;                     //!<! profile to hold all flags for results
  bool fSaveResultsHistograms = false;                   // if results histos are used only as "abstract" interface for binning, then they do not need to be saved
  TProfile* fResultsPro[eAsFunctionOf_N] = {NULL};       //!<! example histogram to store some results + "abstract" interface, which defines common binning, etc., for other groups of histograms.
  TProfile2D* fResultsPro2D[eAsFunctionOf2D_N] = {NULL}; //!<! example histogram to store some results + "abstract" interface, which defines common binning, etc., for other groups of histograms.
  TProfile3D* fResultsPro3D[eAsFunctionOf3D_N] = {NULL}; //!<! example histogram to store some results + "abstract" interface, which defines common binning, etc., for other groups of histograms.

  // Remark: These settings apply to following categories fCorrelationsPro, fNestedLoopsPro, fTest0Pro, fResultsPro, and fParticleSparseHistograms
  TArrayD* fResultsProBinEdges[eAsFunctionOf_N] = {NULL};              // here I keep bin edges uniformly, both for fixed-length binning and variable-length binning
  float fResultsProFixedLengthBins[eAsFunctionOf_N][3] = {{0.}};       // [nBins,min,max]
  TArrayF* fResultsProVariableLengthBins[eAsFunctionOf_N] = {NULL};    // here for each variable in eAsFunctionOf I specify array holding bin boundaries
  bool fUseResultsProVariableLengthBins[eAsFunctionOf_N] = {false};    // use or not variable-length bins
  TString fResultsProVariableLengthBinsString[eAsFunctionOf_N] = {""}; // TBI 20241110 this one is obsolete, can be removed
  TString fResultsProXaxisTitle[eAsFunctionOf_N] = {""};               // keep ordering in sync with enum eAsFunctionOf
  TString fResultsProRawName[eAsFunctionOf_N] = {""};                  // this is how it appears simplified in the 1D hist name when saved to the file
} res;                                                                 // "res" labels an instance of this group of histograms

#endif // PWGCF_MULTIPARTICLECORRELATIONS_CORE_MUPA_DATAMEMBERS_H_
