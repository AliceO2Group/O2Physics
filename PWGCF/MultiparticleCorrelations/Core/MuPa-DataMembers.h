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

#ifndef PWGCF_MULTIPARTICLECORRELATIONS_CORE_MUPA_DATAMEMBERS_H_
#define PWGCF_MULTIPARTICLECORRELATIONS_CORE_MUPA_DATAMEMBERS_H_

// General remarks:
// 0. Starting with C++11, it's possible to initialize data members at declaration, so I do it here
// 1. Use //!<! for introducing a Doxygen comment interpreted as transient in both ROOT 5 and ROOT 6.

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

// a) Base list to hold all output objects ("grandmother" of all lists):

TString sBaseListName = "Default list name"; // yes, I declare it separately, because I need it also later in BailOut() function
OutputObj<TList> fBaseList{sBaseListName.Data(),
                           OutputObjHandlingPolicy::AnalysisObject,
                           OutputObjSourceType::OutputObjSource};
TProfile* fBasePro = NULL; //!<! keeps flags relevant for the whole analysis

// *) Task configuration:
struct TaskConfiguration {
  TString fTaskName = "";                          // task name - this one is used to get the right weights
                                                   // programatically for this analysis
  TString fRunNumber = "";                         // over which run number this task is executed
  Bool_t fRunNumberIsDetermined = kFALSE;          // ensures that run number is determined in process() and propagated to already booked objects only once
  Bool_t fDryRun = kFALSE;                         // book all histos and run without storing and calculating anything
  Bool_t fVerbose = kFALSE;                        // print additional info during debugging, but not for simply utility function or function calls per particle (see next)
  Bool_t fVerboseUtility = kFALSE;                 // print additional info during debugging also for simply utility function, but not for function calls per particle (see next)
  Bool_t fVerboseForEachParticle = kFALSE;         // print additional info during debugging, also for function calls per particle
  Bool_t fDoAdditionalInsanityChecks = kFALSE;     // do additional insanity checks at run time, at the expense of losing a bit of performance
                                                   // (for instance, check if the run number in the current 'collision' is the same as run number in the first 'collision', etc.)
  Bool_t fInsanityCheckForEachParticle = kFALSE;   // do additional insanity checks at run time for each particle, at the expense of losing a lot of performance. Use only during debugging.
  Bool_t fUseCCDB = kFALSE;                        // access personal files from CCDB (kTRUE, this is set as default in Configurables),
                                                   // or from home dir in AliEn (kFALSE, use with care, as this is discouraged)
  Bool_t fProcess[eProcess_N] = {kFALSE};          // set what to process. See enum eProcess for full description. Set via implicit variables within a PROCESS_SWITCH clause.
  TString fWhichProcess = "ProcessRec";            // dump in this variable which process was used
  UInt_t fRandomSeed = 0;                          // argument for TRandom3 constructor. By default it is 0 (seed is guaranteed to be unique in time and space)
  Bool_t fUseFisherYates = kFALSE;                 // algorithm used to randomize particle indices, set via configurable
  TArrayI* fRandomIndices = NULL;                  // array to store random indices obtained from Fisher-Yates algorithm
  Int_t fFixedNumberOfRandomlySelectedTracks = -1; // use a fixed number of randomly selected particles in each event. It is set and applied, if > 0. Set to <=0 to ignore.
  Bool_t fUseStopwatch = kFALSE;                   // do some basing profiling with TStopwatch for where the execution time is going
  TStopwatch* fTimer[eTimer_N] = {NULL};           // stopwatch, global (overal execution time) and local
  Float_t fFloatingPointPrecision = 1.e-6;         // two floats are the same if TMath::Abs(f1 - f2) < fFloatingPointPrecision (there is configurable for it)
  // Bool_t fRescaleWithTheoreticalInput; // if kTRUE, all measured correlators are
  // rescaled with theoretical input, so that in profiles everything is at 1. Used
  // both in OTF and internal val.
} tc; // "tc" labels an instance of this group of variables.

// *) Event-by-event quantities:
struct EventByEventQuantities {
  Int_t fSelectedTracks = 0; // integer counter of tracks used to calculate Q-vectors, after all particle cuts have been applied
  Float_t fCentrality = 0.;  // event-by-event centrality. Value of the default centrality estimator, set via configurable cfCentralityEstimator
  Float_t fOccupancy = 0.;   // event-by-event occupancy. Value of the default occupancy estimator, set via configurable cfOccupancyEstimator
} ebye;                      // "ebye" is a common label for objects in this struct

// *) QA:
//    Remark 1: I keep new histograms in this group, until I need them permanently in the analysis. Then, they are moved to EventHistograms or ParticleHistograms (yes, even if they are 2D).
//    Remark 2: All 2D histograms book as TH2F, due to "stmem error" in terminate (see .cxx for further details)
struct QualityAssurance {
  TList* fQAList = NULL;                                                   //!<! base list to hold all QA output object
  TProfile* fQAHistogramsPro = NULL;                                       //!<! keeps flags relevant for the QA histograms
  Bool_t fCheckUnderflowAndOverflow = kFALSE;                              // check and bail out if in event and particle histograms there are entries which went to underflow or overflow bins
  TH2F* fQAEventHistograms2D[eQAEventHistograms2D_N][2][2] = {{{NULL}}};   //! [ type - see enum eQAEventHistograms2D ][reco,sim][before, after particle cuts]
  Bool_t fFillQAEventHistograms2D = kTRUE;                                 // if kFALSE, all 2D event histograms are not filled. if kTRUE, the ones for which fBookQAEventHistograms2D[...] is kTRUE, are filled
  Bool_t fBookQAEventHistograms2D[eQAEventHistograms2D_N] = {kTRUE};       // book or not this 2D histogram, see configurable cfBookQAEventHistograms2D
  Float_t fEventHistogramsBins2D[eQAEventHistograms2D_N][2][3] = {{{0.}}}; // [type - see enum][x,y][nBins,min,max]
  TString fEventHistogramsName2D[eQAEventHistograms2D_N] = {""};           // name of fQAEventHistograms2D, determined programatically from other 1D names, to ease bookkeeping
  // Int_t fQAEventHistograms2DRebin[eQAEventHistograms2D_N][2] = {{1}};             // to reduce memory consumption, use this number to merge bins together (i.e. "rebinning") [type - see enum][x,y]
  TH2F* fQAParticleHistograms2D[eQAParticleHistograms2D_N][2][2] = {{{NULL}}};   //! [ type - see enum eQAParticleHistograms2D ][reco,sim][before, after particle cuts]
  Bool_t fFillQAParticleHistograms2D = kTRUE;                                    // if kFALSE, all 2D particle histograms are not filled. if kTRUE, the ones for which fBookQAParticleHistograms2D[...] is kTRUE, are filled
  Bool_t fBookQAParticleHistograms2D[eQAParticleHistograms2D_N] = {kTRUE};       // book or not this 2D histogram, see configurable cfBookQAParticleHistograms2D
  Float_t fParticleHistogramsBins2D[eQAParticleHistograms2D_N][2][3] = {{{0.}}}; // [type - see enum][x,y][nBins,min,max]
  TString fParticleHistogramsName2D[eQAParticleHistograms2D_N] = {""};           // name of fQAParticleHistograms2D, determined programatically from other 1D names, to ease bookkeeping
  // Int_t fQAParticleHistograms2DRebin[eQAParticleHistograms2D_N][2] = {{1}};       // to reduce memory consumption, use this number to merge bins together (i.e. "rebinning") [type - see enum][x,y]
  Float_t fCentrality[eCentralityEstimators_N] = {0.};              // used mostly in QA correlation plots
  TString fCentralityEstimatorName[eCentralityEstimators_N] = {""}; //
  Float_t fOccupancy[eOccupancyEstimators_N] = {0.};                // used mostly in QA correlation plots
  TString fOccupancyEstimatorName[eOccupancyEstimators_N] = {""};   //
} qa;                                                               // "qa" is a common label for objects in this struct

// *) Event histograms:
struct EventHistograms {
  TList* fEventHistogramsList = NULL;   //!<! list to hold all control event histograms
  TProfile* fEventHistogramsPro = NULL; //!<! keeps flags relevant for the control event histograms
  // 1D:
  TH1F* fEventHistograms[eEventHistograms_N][2][2] = {{{NULL}}}; //! [ type - see enum eEventHistograms ][reco,sim][before, after event cuts]
  Bool_t fFillEventHistograms = kTRUE;                           // if kFALSE, all event histograms are not filled. if kTRUE, the ones for which fBookEventHistograms[...] is kTRUE, are filled
  Bool_t fBookEventHistograms[eEventHistograms_N] = {kTRUE};     // book or not this histogram, see SetBookEventHistograms
  Float_t fEventHistogramsBins[eEventHistograms_N][3] = {{0.}};  // [nBins,min,max]
  TString fEventHistogramsName[eEventHistograms_N] = {""};       // name of event histogram, used both for 1D and 2D histograms
  // 2D:
  // ...
  // Remark: All 2D event histograms are still in the QA group. Move here only the ones I will use regularly in the analysis
} eh; // "eh" labels an instance of group of histograms "EventHistograms"

// *) Event cuts:
struct EventCuts {
  TList* fEventCutsList = NULL;                            //!<! list to hold all event cuts objects
  TProfile* fEventCutsPro = NULL;                          //!<! keeps flags relevant for the event cuts
  Bool_t fUseEventCuts[eEventCuts_N] = {kFALSE};           // Use or do not use a cut enumerated in eEventHistograms + eEventCuts
  Bool_t fUseEventCutCounterAbsolute = kFALSE;             // profile and save how many times each event cut counter triggered (absolute). Use with care, as this is computationally heavy
  Bool_t fUseEventCutCounterSequential = kFALSE;           // profile and save how many times each event cut counter triggered (sequential). Use with care, as this is computationally heavy
  Bool_t fEventCutCounterBinLabelingIsDone = kFALSE;       // this flag ensures that ordered labeling of bins, to resemble ordering of cut implementation, is done only once.
  Bool_t fPrintCutCounterContent = kFALSE;                 // if true, prints on the screen content of fEventCutCounterHist[][] (all which were booked)
  TString fEventCutName[eEventCuts_N] = {""};              // event cut name, with default ordering defined by ordering in enum eEventCuts
  TExMap* fEventCutCounterMap[2] = {NULL};                 // map (key, value) = (enum eEventCuts, ordered bin number)
  TExMap* fEventCutCounterMapInverse[2] = {NULL};          // inverse of above fEventCutCounterMap, i.e. (ordered bin number, enum eEventCuts)
  Int_t fEventCutCounterBinNumber[2] = {1, 1};             // bin counter for set bin labels in fEventCutCounterHist
  Float_t fdEventCuts[eEventCuts_N][2] = {{0.}};           // event cuts defined via [min,max)
  TString fsEventCuts[eEventCuts_N] = {""};                // event cuts defined via string
  TH1I* fEventCutCounterHist[2][eCutCounter_N] = {{NULL}}; //!<! [rec,sim][see enum eCutCounter] histogram to store how many any times each event cut triggered
  Int_t fBeforeAfterColor[2] = {kRed, kGreen};             // color code before and after cuts
} ec;                                                      // "ec" is a common label for objects in this struct

// *) Particle histograms:
struct ParticleHistograms {
  TList* fParticleHistogramsList = NULL;   //!<! list to hold all control particle histograms
  TProfile* fParticleHistogramsPro = NULL; //!<! keeps flags relevant for the control particle histograms
  // 1D:
  TH1F* fParticleHistograms[eParticleHistograms_N][2][2] = {{{NULL}}}; //! [ type - see enum eParticleHistograms ][reco,sim][before, after particle cuts]
  Bool_t fFillParticleHistograms = kTRUE;                              // if kFALSE, all 1D particle histograms are not filled.
                                                                       // if kTRUE, the ones for which fBookParticleHistograms[...] is kTRUE, are filled
  Bool_t fBookParticleHistograms[eParticleHistograms_N] = {kTRUE};     // book or not the particular particle histogram, see configurable cfBookParticleHistograms
  Float_t fParticleHistogramsBins[eParticleHistograms_N][3] = {{0.}};  // [nBins,min,max]
  TString fParticleHistogramsName[eParticleHistograms_N] = {""};       // name of particle histogram, used both for 1D and 2D histograms
  // 2D:
  TH2D* fParticleHistograms2D[eParticleHistograms2D_N][2][2] = {{{NULL}}};     //! [ type - see enum eParticleHistograms2D ][reco,sim][before, after particle cuts]
  Bool_t fFillParticleHistograms2D = kTRUE;                                    // if kFALSE, all 2D particle histograms are not filled.
                                                                               // if kTRUE, the ones for which fBookParticleHistograms2D[...] is kTRUE, are filled
  Bool_t fBookParticleHistograms2D[eParticleHistograms2D_N] = {kTRUE};         // book or not this 2D histogram, see configurable cfBookParticleHistograms2D
  Float_t fParticleHistogramsBins2D[eParticleHistograms2D_N][2][3] = {{{0.}}}; // [type - see enum][x,y][nBins,min,max]
  TString fParticleHistogramsName2D[eParticleHistograms2D_N] = {""};           // name of particle histogram 2D, determined programatically from two 1D, in the format "%s_vs_%s"
} ph;                                                                          // "ph" labels an instance of group of histograms "ParticleHistograms"

// *) Particle cuts:
struct ParticleCuts {
  TList* fParticleCutsList = NULL;                            //!<! list to hold all particle cuts objects
  TProfile* fParticleCutsPro = NULL;                          //!<! keeps flags relevant for the particle cuts
  Bool_t fUseParticleCuts[eParticleCuts_N] = {kFALSE};        // true or false .
  Bool_t fUseParticleCutCounterAbsolute = kFALSE;             // profile and save how many times each particle cut counter triggered (absolute). Use with care, as this is computationally heavy
  Bool_t fUseParticleCutCounterSequential = kFALSE;           // profile and save how many times each particle cut counter triggered (sequential). Use with care, as this is computationally heavy
  Bool_t fParticleCutCounterBinLabelingIsDone = kFALSE;       // this flag ensures that ordered labeling of bins, to resemble ordering of cut implementation, is done only once.
  TString fParticleCutName[eParticleCuts_N] = {""};           // particle cut name, as used in bin labels, etc.
  TExMap* fParticleCutCounterMap[2] = {NULL};                 // map (key, value) = (enum eParticleCuts, ordered bin number)
  TExMap* fParticleCutCounterMapInverse[2] = {NULL};          // inverse of above fParticleCutCounterMap, i.e. (ordered bin number, enum eParticleCuts)
  Int_t fParticleCutCounterBinNumber[2] = {1, 1};             // bin counter for set bin labels in fParticleCutCounterHist
  Float_t fdParticleCuts[eParticleCuts_N][2] = {{0.}};        // particles cuts defined via [min,max) . Remark: I use here eParticleHistograms_N , not to duplicate these enums for ParticleCuts.
  TString fsParticleCuts[eParticleCuts_N] = {""};             // particles cuts defined via booleans via string
  TH1I* fParticleCutCounterHist[2][eCutCounter_N] = {{NULL}}; //!<! [rec,sim][see enum eCutCounter] histogram to store how many any times each particle cut triggered
  TFormula* fPtDependentDCAxyFormula = NULL;                  // the actual formula, used to evaluate for a given pT, the corresponding DCAxy, where the parameterization is given by configurable cfPtDependentDCAxyParameterization
} pc;                                                         // "pc" is a common label for objects in this struct

// *) Q-vectors:
struct Qvector {
  TList* fQvectorList = NULL;                                                                                                          // list to hold all Q-vector objects
  TProfile* fQvectorFlagsPro = NULL;                                                                                                   // profile to hold all flags for Q-vector
  Bool_t fCalculateQvectors = kTRUE;                                                                                                   // to calculate or not to calculate Q-vectors, that's a Boolean...
  TComplex fQ[gMaxHarmonic * gMaxCorrelator + 1][gMaxCorrelator + 1] = {{TComplex(0., 0.)}};                                           //! generic Q-vector
  TComplex fQvector[gMaxHarmonic * gMaxCorrelator + 1][gMaxCorrelator + 1] = {{TComplex(0., 0.)}};                                     //! "integrated" Q-vector
  TComplex fqvector[eqvectorKine_N][gMaxNoBinsKine][gMaxHarmonic * gMaxCorrelator + 1][gMaxCorrelator + 1] = {{{{TComplex(0., 0.)}}}}; //! "differenttial" q-vector [kine var.][binNo][fMaxHarmonic*fMaxCorrelator+1][fMaxCorrelator+1] = [6*12+1][12+1]
  Int_t fqVectorEntries[eqvectorKine_N][gMaxNoBinsKine] = {{0}};                                                                       // count number of entries in each differential q-vector
} qv;                                                                                                                                  // "qv" is a common label for objects in this struct

// *) Multiparticle correlations (standard, isotropic, same harmonic):
struct MultiparticleCorrelations {
  TList* fCorrelationsList = NULL;                                                                     // list to hold all correlations objects
  TProfile* fCorrelationsFlagsPro = NULL;                                                              // profile to hold all flags for correlations
  Bool_t fCalculateCorrelations = kTRUE;                                                               // calculate and store integrated correlations
  TProfile* fCorrelationsPro[4][gMaxHarmonic][eAsFunctionOf_N] = {{{NULL}}};                           //! multiparticle correlations
                                                                                                       //  [2p=0,4p=1,6p=2,8p=3][n=1,n=2,...,n=gMaxHarmonic]
                                                                                                       //  [0=integrated,1=vs. multiplicity,2=vs. centrality,3=pT,4=eta,5=vs. occupancy]
  Bool_t fCalculateCorrelationsAsFunctionOf[eAsFunctionOf_N] = {true, true, true, false, false, true}; //! [0=integrated,1=vs. multiplicity,2=vs. centrality,3=pT,4=eta,5=vs. occupancy]
                                                                                                       //  As of 20241111, 3=pT and 4=eta are not implemented, see void CalculateKineCorrelations(...)
} mupa;                                                                                                // "mupa" is a common label for objects in this struct

// *) Particle weights:
struct ParticleWeights {
  TList* fWeightsList = NULL;                                             //!<! list to hold all particle weights
  TProfile* fWeightsFlagsPro = NULL;                                      //!<! profile to hold all flags for weights
  Bool_t fUseWeights[eWeights_N] = {false};                               // use weights [phi,pt,eta]
  TH1D* fWeightsHist[eWeights_N] = {NULL};                                //!<! particle weights
  Bool_t fUseDiffWeights[eDiffWeights_N] = {false};                       // use differential weights [phipt,phieta]
  TH1D* fDiffWeightsHist[eDiffWeights_N][fMaxBinsDiffWeights] = {{NULL}}; // histograms holding differential weights [phipt,phieta][bin number]
  TString fFileWithWeights = "";                                          // path to external ROOT file which holds all particle weights
  Bool_t fParticleWeightsAreFetched = kFALSE;                             // ensures that particle weights are fetched only once
} pw;                                                                     // "pw" labels an instance of this group of histograms

// *) Nested loops:
struct NestedLoops {
  TList* fNestedLoopsList = NULL;                                              // list to hold all nested loops objects
  TProfile* fNestedLoopsFlagsPro = NULL;                                       // profile to hold all flags for nested loops
  Bool_t fCalculateNestedLoops = kFALSE;                                       // calculate and store correlations with nested loops, as a cross-check
  Bool_t fCalculateCustomNestedLoops = kFALSE;                                 // validate e-b-e all correlations with custom nested loop
  Bool_t fCalculateKineCustomNestedLoops = kFALSE;                             // validate e-b-e all differential (vs pt, eta, etc.) correlations with custom nested loop
  Int_t fMaxNestedLoop = -1;                                                   // if set to e.g. 4, all nested loops beyond that, e.g. 6-p and 8-p, are NOT calculated
  TProfile* fNestedLoopsPro[4][gMaxHarmonic][eAsFunctionOf_N] = {{{NULL}}};    //! multiparticle correlations from nested loops
                                                                               //! [2p=0,4p=1,6p=2,8p=3][n=1,n=2,...,n=gMaxHarmonic][0=integrated,1=vs.
                                                                               //! multiplicity,2=vs. centrality,3=pT,4=eta]
  TArrayD* ftaNestedLoops[2] = {NULL};                                         //! e-b-e container for nested loops [0=angles;1=product of all weights]
  TArrayD* ftaNestedLoopsKine[eqvectorKine_N][gMaxNoBinsKine][2] = {{{NULL}}}; //! e-b-e container for nested loops // [0=pT,1=eta][kine bin][0=angles;1=product of all weights]
} nl;                                                                          // "nl" labels an instance of this group of histograms

// *) Toy NUA (can be applied both in real data analysis and in analysis 'on-the-fly'):
struct NUA {
  TList* fNUAList = NULL;                              // list to hold all NUA objects
  TProfile* fNUAFlagsPro = NULL;                       // profile to hold all flags for NUA objects
  Bool_t fApplyNUAPDF[eNUAPDF_N] = {kFALSE};           // apply NUA to particular kine variable (see the corresponding enum eNUAPDF)
  Bool_t fUseDefaultNUAPDF[eNUAPDF_N] = {kTRUE};       // by default, use simple hardcoded expressions for NUA acceptance profile
  TF1* fDefaultNUAPDF[eNUAPDF_N] = {NULL};             // default distributions used as pdfs to simulate events on-the-fly
  TH1D* fCustomNUAPDF[eNUAPDF_N] = {NULL};             // custom, user-supplied distributions used to simulate NUA
  TString* fCustomNUAPDFHistNames[eNUAPDF_N] = {NULL}; // these are the names of histograms holding custom NUA in an external file. There is a configurable for this one.
  TString fFileWithCustomNUA = "";                     // path to external ROOT file which holds all histograms with custom NUA
  Double_t fMaxValuePDF[eNUAPDF_N] = {0.};             // see algorithm used in Accept(...). I implemented it as a data member, so that it is not calculated again and again at each particle call
} nua;

// *) Internal validation:
struct InternalValidation {
  TList* fInternalValidationList = NULL;              // list to hold all objects for internal validation
  TProfile* fInternalValidationFlagsPro = NULL;       // profile to hold all flags for internal validation
  Bool_t fUseInternalValidation = kFALSE;             // use internal validation
  Bool_t fInternalValidationForceBailout = kFALSE;    // force bailout in internal validation after either eNumberOfEvents or eSelectedEvents is reached.
                                                      // This is OK as long as I do not apply any event cuts in InternalValidation().
                                                      // Remember that for each real event, I do fnEventsInternalValidation events on-the-fly.
  UInt_t fnEventsInternalValidation = 0;              // how many on-the-fly events will be sampled for each real event, for internal validation
  TString* fHarmonicsOptionInternalValidation = NULL; // "constant" or "correlated", see .cxx for full documentation
  Bool_t fRescaleWithTheoreticalInput = kFALSE;       // if kTRUE, all measured correlators are rescaled with theoretical input, so that in profiles everything is at 1
  TArrayD* fInternalValidationVnPsin[2] = {NULL};     // 0 = { v1, v2, ... }, 1 = { Psi1, Psi2, ... }
  Int_t fMultRangeInternalValidation[2] = {0, 0};     // min and max values for uniform multiplicity distribution in on-the-fly analysis (convention: min <= M < max)
} iv;

// *) Test0:
struct Test0 {
  TList* fTest0List = NULL;                                                                      // list to hold all objects for Test0
  TProfile* fTest0FlagsPro = NULL;                                                               // store all flags for Test0
  Bool_t fCalculateTest0 = kFALSE;                                                               // calculate or not Test0
  TProfile* fTest0Pro[gMaxCorrelator][gMaxIndex][eAsFunctionOf_N] = {{{NULL}}};                  //! [order][index][0=integrated,1=vs. multiplicity,2=vs. centrality,3=pT,4=eta]
  TString* fTest0Labels[gMaxCorrelator][gMaxIndex] = {{NULL}};                                   // all labels: k-p'th order is stored in k-1'th index. So yes, I also store 1-p
  Bool_t fCalculateTest0AsFunctionOf[eAsFunctionOf_N] = {true, true, true, false, false, false}; //! [0=integrated,1=vs. multiplicity,2=vs. centrality,3=pT,4=eta,5=vs. occupancy]
  TString fFileWithLabels = "";                                                                  // path to external ROOT file which specifies all labels of interest
  TH1I* fTest0LabelsPlaceholder = NULL;                                                          // store all Test0 labels in this histogram
} t0;                                                                                            // "t0" labels an instance of this group of histograms

// *) Results:
struct Results {                                   // This is in addition also sort of "abstract" interface, which defines common binning, etc., for other groups of histograms.
  TList* fResultsList = NULL;                      //!<! list to hold all results
  TProfile* fResultsFlagsPro = NULL;               //!<! profile to hold all flags for results
  Bool_t fSaveResultsHistograms = false;           // if results histos are used only as "abstract" interface for binning, then they do not need to be saved
  TProfile* fResultsPro[eAsFunctionOf_N] = {NULL}; //!<! example histogram to store some results + "abstract" interface, which defines common binning, etc., for other groups of histograms.

  // Remark: These settings apply to following categories fCorrelationsPro, fNestedLoopsPro, fTest0Pro, and fResultsHist
  Float_t fResultsProFixedLengthBins[eAsFunctionOf_N][3] = {{0.}};                                                             // [nBins,min,max]
  TArrayF* fResultsProVariableLengthBins[eAsFunctionOf_N] = {NULL};                                                            // here for each variable in eAsFunctionOf I specify array holding bin boundaries
  Bool_t fUseResultsProVariableLengthBins[eAsFunctionOf_N] = {kFALSE};                                                         // use or not variable-length bins
  TString fResultsProVariableLengthBinsString[eAsFunctionOf_N] = {""};                                                         // TBI 20241110 this one is obsolete, can be removed
  TString fResultsProXaxisTitle[eAsFunctionOf_N] = {"integrated", "multiplicity", "centrality", "p_{T}", "#eta", "occupancy"}; // keep ordering in sync with enum eAsFunctionOf
  TString fResultsProRawName[eAsFunctionOf_N] = {"int", "mult", "cent", "pt", "eta", "occu"};                                  // this is how it appears simplified in the hist name when saved to the file
} res;                                                                                                                         // "res" labels an instance of this group of histograms

#endif // PWGCF_MULTIPARTICLECORRELATIONS_CORE_MUPA_DATAMEMBERS_H_
