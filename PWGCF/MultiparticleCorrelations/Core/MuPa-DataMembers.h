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
  Bool_t fVerbose = kFALSE;                        // print additional info like Green(__PRETTY_FUNCTION__); etc., to
                                                   // be used during debugging, but not for function calls per particle
  Bool_t fVerboseForEachParticle = kFALSE;         // print additional info like Green(__PRETTY_FUNCTION__); etc., to
                                                   // be used during debugging, also for function calls per particle
  Bool_t fDoAdditionalInsanityChecks = kFALSE;     // do additional insanity checks at run time, at the expense of losing a bit of performance
                                                   // (for instance, check if the run number in the current 'collision' is the same as run number in the first 'collision', etc.)
  Bool_t fInsanityCheckForEachParticle = kFALSE;   // do additional insanity checks at run time for each particle, at the expense of losing a lot of performance. Use only during debugging.
  Bool_t fUseCCDB = kFALSE;                        // access personal files from CCDB (kTRUE, this is set as default in Configurables),
                                                   // or from home dir in AliEn (kFALSE, use with care, as this is discouraged)
  Bool_t fProcess[eProcess_N] = {kFALSE};          // Set what to process. See enum eProcess for full description. Set via implicit variables within a PROCESS_SWITCH clause.
  TString fWhichProcess = "ProcessRec";            // dump in this variable which process was used
  UInt_t fRandomSeed = 0;                          // argument for TRandom3 constructor. By default it is 0 (seed is guaranteed to be unique in time and space)
  Bool_t fUseFisherYates = kFALSE;                 // algorithm used to randomize particle indices, set via configurable
  TArrayI* fRandomIndices = NULL;                  // array to store random indices obtained from Fisher-Yates algorithm
  Int_t fFixedNumberOfRandomlySelectedTracks = -1; // use a fixed number of randomly selected particles in each event. It is set and applied, if > 0. Set to <=0 to ignore.
  Bool_t fUseStopwatch = kFALSE;                   // do some basing profiling with TStopwatch for where the execution time is going
  TStopwatch* fTimer[eTimer_N] = {NULL};           // stopwatch, global (overal execution time) and local
  // Bool_t fRescaleWithTheoreticalInput; // if kTRUE, all measured correlators are
  // rescaled with theoretical input, so that in profiles everything is at 1. Used
  // both in OTF and internal val.
} tc; // "tc" labels an instance of this group of variables.

// *) Event-by-event quantities:
struct EventByEventQuantities {
  Int_t fSelectedTracks = 0; // integer counter of tracks used to calculate Q-vectors, after all particle cuts have been applied
  Double_t fCentrality = 0.; // event-by-event centrality from default estimator
} ebye;                      // "ebye" is a common label for objects in this struct

// *) QA:
struct QualityAssurance {
  TList* fQAList = NULL; //!<! base list to hold all QA output object
} qa;                    // "qa" is a common label for objects in this struct

// *) Event histograms:
struct EventHistograms {
  TList* fEventHistogramsList = NULL;                            //!<! list to hold all control event histograms
  TProfile* fEventHistogramsPro = NULL;                          //!<! keeps flags relevant for the control event histograms
  TH1D* fEventHistograms[eEventHistograms_N][2][2] = {{{NULL}}}; //! [ type - see enum eEventHistograms ][reco,sim][before, after event cuts]
  Bool_t fBookEventHistograms[eEventHistograms_N] = {kTRUE};     // book or not this histogram, see SetBookEventHistograms
  Double_t fEventHistogramsBins[eEventHistograms_N][3] = {{0.}}; // [nBins,min,max]
} eh;                                                            // "eh" labels an instance of group of histograms "EventHistograms"

// *) Event cuts:
struct EventCuts {
  TList* fEventCutsList = NULL;                   //!<! list to hold all event cuts objects
  TProfile* fEventCutsPro = NULL;                 //!<! keeps flags relevant for the event cuts
  Bool_t fUseEventCuts[eEventCuts_N] = {kFALSE};  // Use or do not use a cut enumerated in eEventHistograms + eEventCuts
  Double_t fdEventCuts[eEventCuts_N][2] = {{0.}}; // [min,max)
  TString fsEventCuts[eEventCuts_N] = {""};       // specific option passed via string
} ec;                                             // "ec" is a common label for objects in this struct

// *) Particle histograms:
struct ParticleHistograms {
  TList* fParticleHistogramsList = NULL;                                        //!<! list to hold all control particle histograms
  TProfile* fParticleHistogramsPro = NULL;                                      //!<! keeps flags relevant for the control particle histograms
  TH1D* fParticleHistograms[eParticleHistograms_N][2][2] = {{{NULL}}};          //! [ type - see enum eParticleHistograms ][reco,sim][before, after particle cuts]
  Bool_t fBookParticleHistograms[eParticleHistograms_N] = {kTRUE};              // book or not this histogram, see configurable cfBookParticleHistograms
  Double_t fParticleHistogramsBins[eParticleHistograms_N][3] = {{0.}};          // [nBins,min,max]
  Double_t fParticleCuts[eParticleHistograms_N][2] = {{0.}};                    // [min,max]
  TH2D* fParticleHistograms2D[eParticleHistograms2D_N][2][2] = {{{NULL}}};      //! [ type - see enum eParticleHistograms2D ][reco,sim][before, after particle cuts]
  Bool_t fBookParticleHistograms2D[eParticleHistograms2D_N] = {kTRUE};          // book or not this 2D histogram, see configurable cfBookParticleHistograms2D
  Double_t fParticleHistogramsBins2D[eParticleHistograms2D_N][2][3] = {{{0.}}}; // [type - see enum][x,y][nBins,min,max]
} ph;                                                                           // "ph" labels an instance of group of histograms "ParticleHistograms"

// *) Particle cuts:
struct ParticleCuts {
  TList* fParticleCutsList = NULL;                      //!<! list to hold all particle cuts objects
  TProfile* fParticleCutsPro = NULL;                    //!<! keeps flags relevant for the particle cuts
  Bool_t fUseParticleCuts[eParticleCuts_N] = {kFALSE};  // true or false .
  Double_t fdParticleCuts[eParticleCuts_N][2] = {{0.}}; // [min,max) . Remark: I use here eParticleHistograms_N , not to duplicate these enums for ParticleCuts.
  TString fsParticleCuts[eParticleCuts_N] = {""};       // specific option passed via string
} pc;                                                   // "pc" is a common label for objects in this struct

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
  TList* fCorrelationsList = NULL;                                           // list to hold all correlations objects
  TProfile* fCorrelationsFlagsPro = NULL;                                    // profile to hold all flags for correlations
  Bool_t fCalculateCorrelations = kTRUE;                                     // calculate and store integrated correlations
  TProfile* fCorrelationsPro[4][gMaxHarmonic][eAsFunctionOf_N] = {{{NULL}}}; //! multiparticle correlations
                                                                             //! [2p=0,4p=1,6p=2,8p=3][n=1,n=2,...,n=gMaxHarmonic][0=integrated,1=vs.
                                                                             //! multiplicity,2=vs. centrality,3=pT,4=eta]
} mupa;                                                                      // "mupa" is a common label for objects in this struct

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

// *) Internal validation:
struct InternalValidation {
  TList* fInternalValidationList = NULL;              // list to hold all objects for internal validation
  TProfile* fInternalValidationFlagsPro = NULL;       // profile to hold all flags for internal validation
  Bool_t fUseInternalValidation = kFALSE;             // use internal validation
  Bool_t fInternalValidationForceBailout = kFALSE;    // force bailout after fnEventsInternalValidation is reached. In HL, for each real event, I do fnEventsInternalValidation events
  UInt_t fnEventsInternalValidation = 0;              // how many events will be sampled on-the-fly for internal validation
  TString* fHarmonicsOptionInternalValidation = NULL; // see .cxx for full documentation
  Bool_t fRescaleWithTheoreticalInput = kFALSE;       // if kTRUE, all measured correlators are rescaled with theoretical input, so that in profiles everything is at 1
  TArrayD* fInternalValidationVnPsin[2] = {NULL};     // 0 = { v1, v2, ... }, 1 = { Psi1, Psi2, ... }
  Int_t fMultRangeInternalValidation[2] = {0, 0};     // min and max values for uniform multiplicity distribution in on-the-fly analysis (convention: min <= M < max)
} iv;

// *) Test0:
struct Test0 {
  TList* fTest0List = NULL;                                                               // list to hold all objects for Test0
  TProfile* fTest0FlagsPro = NULL;                                                        // store all flags for Test0
  Bool_t fCalculateTest0 = kFALSE;                                                        // calculate or not Test0
  TProfile* fTest0Pro[gMaxCorrelator][gMaxIndex][eAsFunctionOf_N] = {{{NULL}}};           //! [order][index][0=integrated,1=vs. multiplicity,2=vs. centrality,3=pT,4=eta]
  TString* fTest0Labels[gMaxCorrelator][gMaxIndex] = {{NULL}};                            // all labels: k-p'th order is stored in k-1'th index. So yes, I also store 1-p
  Bool_t fCalculateTest0AsFunctionOf[eAsFunctionOf_N] = {true, true, true, false, false}; //! [0=integrated,1=vs. multiplicity,2=vs. centrality,3=pT,4=eta]
  TString fFileWithLabels = "";                                                           // path to external ROOT file which specifies all labels of interest
  TH1I* fTest0LabelsPlaceholder = NULL;                                                   // store all Test0 labels in this histogram
} t0;                                                                                     // "t0" labels an instance of this group of histograms

// *) Results:
struct Results {                                   // This is in addition also sort of "abstract" interface, which defines common binning, etc., for other groups of histograms.
  TList* fResultsList = NULL;                      //!<! list to hold all results
  TProfile* fResultsFlagsPro = NULL;               //!<! profile to hold all flags for results
  Bool_t fSaveResultsHistograms = false;           // if results histos are used only as "abstract" interface for binning, then they do not need to be saved
  TProfile* fResultsPro[eAsFunctionOf_N] = {NULL}; //!<! example histogram to store some results + "abstract" interface, which defines common binning, etc., for other groups of histograms.

  // Remark: These settings apply to following categories fCorrelationsPro, fNestedLoopsPro, fTest0Pro, and fResultsHist
  Float_t fResultsProFixedLengthBins[eAsFunctionOf_N][3] = {{0.}};                                                // [nBins,min,max]
  TArrayD* fResultsProVariableLengthBins[eAsFunctionOf_N] = {NULL};                                               // here for each variable in eAsFunctionOf I specify array holding bin boundaries
  Bool_t fUseResultsProVariableLengthBins[eAsFunctionOf_N] = {kFALSE};                                            // use or not variable-length bins
  TString fResultsProVariableLengthBinsString[eAsFunctionOf_N] = {""};                                            // TBI 20240113 temporary I do it this way
  TString fResultsProXaxisTitle[eAsFunctionOf_N] = {"integrated", "multiplicity", "centrality", "p_{T}", "#eta"}; // keep ordering in sync with enum eAsFunctionOf
  TString fResultsProRawName[eAsFunctionOf_N] = {"int", "mult", "cent", "pt", "eta"};                             // this is how it appears simplified in the hist name when saved to the file
} res;                                                                                                            // "res" labels an instance of this group of histograms

#endif // PWGCF_MULTIPARTICLECORRELATIONS_CORE_MUPA_DATAMEMBERS_H_
