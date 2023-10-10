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

// Remarks:
// 0. Starting with C++11, it's possible to initialize data members at
// declaration, so I do it here
// 1. Use //!<! for introducing a Doxygen comment interpreted as transient in
// both ROOT 5 and ROOT 6.

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

OutputObj<TList> fBaseList{"Default list name",
                           OutputObjHandlingPolicy::AnalysisObject,
                           OutputObjSourceType::OutputObjSource};
TProfile* fBasePro = NULL; //!<! keeps flags relevant for the whole analysis
UInt_t fRandomSeed = 0;    // argument to TRandom3 constructor. By default is 0,
                           // use SetRandomSeed(...) to change it

// *) Task configuration:
TString fTaskName = "";  // task name - this one is used to get the right weights
                         // programatically for this analysis
TString fRunNumber = ""; // over which run number this task is executed
Bool_t fVerbose =
  kFALSE; // print additional info like Green(__PRETTY_FUNCTION__); etc., to
          // be used during debugging, but not for function calls per particle
Bool_t fVerboseForEachParticle =
  kFALSE; // print additional info like Green(__PRETTY_FUNCTION__); etc., to
          // be used during debugging, also for function calls per particle
Bool_t fUseCCDB =
  kFALSE; // access personal files from CCDB (kTRUE, this is set as default in
          // Configurables), or from home dir in AliEn (kFALSE, use with care,
          // as this is discouraged)
// TString fModusOperandi; // "Rec" = process only reconstructed, "Sim" =
// process only simulated, "RecSim" = process both reconstructed and simulated
// UInt_t fRandomSeed; // argument to TRandom3 constructor. By default it is 0,
// use SetRandomSeed(...) to change it Bool_t fUseFisherYates; // use
// SetUseFisherYates(kTRUE); in the steering macro to randomize particle indices
// Bool_t fUseFixedNumberOfRandomlySelectedTracks; // use or not fixed number of
// randomly selected particles in each event. Use always in combination with
// SetUseFisherYates(kTRUE) Int_t fFixedNumberOfRandomlySelectedTracks; // set
// here a fixed number of randomly selected particles in each event. Use always
// in combination with SetUseFisherYates(kTRUE) Bool_t
// fRescaleWithTheoreticalInput; // if kTRUE, all measured correlators are
// rescaled with theoretical input, so that in profiles everything is at 1. Used
// both in OTF and internal val.

// *) Event-by-event quantities:
Int_t fSelectedTracks =
  0;                       // integer counter of tracks used to calculate Q-vectors, after all
                           // particle cuts have been applied
Double_t fCentrality = 0.; // event-by-event centrality from default estimator

// *) QA:
TList* fQAList = NULL; //!<! base list to hold all QA output object

// *) Event histograms:
TList* fEventHistogramsList =
  NULL; //!<! list to hold all control event histograms
TProfile* fEventHistogramsPro =
  NULL; //!<! keeps flags relevant for the control event histograms
struct EventHistograms_Arrays {
  TH1D* fEventHistograms[eEventHistograms_N][2][2] = {
    {{NULL}}}; //! [ type - see enum eEventHistograms ][reco,sim][before,after
               //! event cuts]
  Bool_t fBookEventHistograms[eEventHistograms_N] = {
    kTRUE}; // book or not this histogram, see SetBookEventHistograms
  Double_t fEventHistogramsBins[eEventHistograms_N][3] = {
    {0.}};                                             // [nBins,min,max]
  Double_t fEventCuts[eEventHistograms_N][2] = {{0.}}; // [min,max]
} ceh_a;                                               // "ceh_a" labels an instance of this group of histograms

// *) Particle histograms:
TList* fParticleHistogramsList =
  NULL; //!<! list to hold all control particle histograms
TProfile* fParticleHistogramsPro =
  NULL; //!<! keeps flags relevant for the control particle histograms
struct ParticleHistograms_Arrays {
  TH1D* fParticleHistograms[eParticleHistograms_N][2][2] = {
    {{NULL}}}; //! [ type - see enum eParticleHistograms
               //! ][reco,sim][before,after particle cuts]
  Bool_t fBookParticleHistograms[eParticleHistograms_N] = {
    kTRUE}; // book or not this histogram, see SetBookParticleHistograms
  Double_t fParticleHistogramsBins[eParticleHistograms_N][3] = {
    {0.}};                                                   // [nBins,min,max]
  Double_t fParticleCuts[eParticleHistograms_N][2] = {{0.}}; // [min,max]
} cph_a;                                                     // "cph_a" labels an instance of this group of histograms

// *) Q-vectors:
TList* fQvectorList = NULL;        // list to hold all Q-vector objects
TProfile* fQvectorFlagsPro = NULL; // profile to hold all flags for Q-vector
Bool_t fCalculateQvector =
  kTRUE; // to calculate or not to calculate Q-vectors, that's a Boolean...
struct Qvector_Arrays {
  TComplex fQ[gMaxHarmonic * gMaxCorrelator + 1][gMaxCorrelator + 1] = {
    {TComplex(0., 0.)}}; //! generic Q-vector
  TComplex fQvector[gMaxHarmonic * gMaxCorrelator + 1][gMaxCorrelator + 1] = {
    {TComplex(0., 0.)}}; //! "integrated" Q-vector
} qv_a;

// *) Multiparticle correlations (standard, isotropic, same harmonic):
TList* fCorrelationsList = NULL; // list to hold all correlations objects
TProfile* fCorrelationsFlagsPro =
  NULL; // profile to hold all flags for correlations
Bool_t fCalculateCorrelations =
  kTRUE; // calculate and store integrated correlations
struct Correlations_Arrays {
  TProfile* fCorrelationsPro[4][gMaxHarmonic][3] = {
    {{NULL}}}; //! multiparticle correlations
               //! [2p=0,4p=1,6p=2,8p=3][n=1,n=2,...,n=6][0=integrated,1=vs.
               //! multiplicity,2=vs. centrality]
} c_a;

// *) Particle weights:
TList* fWeightsList = NULL;        //!<! list to hold all particle weights
TProfile* fWeightsFlagsPro = NULL; //!<! profile to hold all flags for weights
struct ParticleWeights_Arrays {
  Bool_t fUseWeights[eWeights_N] = {false}; // use weights [phi,pt,eta]
  TH1D* fWeightsHist[eWeights_N] = {NULL};  //!<! particle weights
} pw_a;                                     // "pw_a" labels an instance of this group of histograms, e.g.
                                            // pw_a.fWeightsHist[0]
TString fFileWithWeights =
  ""; // path to external ROOT file which holds all particle weights

// *) Nested loops:
TList* fNestedLoopsList = NULL; // list to hold all nested loops objects
TProfile* fNestedLoopsFlagsPro =
  NULL;                               // profile to hold all flags for nested loops
Bool_t fCalculateNestedLoops = kTRUE; // calculate and store correlations with
                                      // nested loops, as a cross-check
Bool_t fCalculateCustomNestedLoop =
  kFALSE; // validate e-b-e all correlations with custom nested loop
struct NestedLoops_Arrays {
  TProfile* fNestedLoopsPro[4][6][5] = {
    {{NULL}}};                         //! multiparticle correlations from nested loops
                                       //! [2p=0,4p=1,6p=2,8p=3][n=1,n=2,...,n=6][0=integrated,1=vs.
                                       //! multiplicity,2=vs. centrality,3=pT,4=eta]
  TArrayD* ftaNestedLoops[2] = {NULL}; //! e-b-e container for nested loops
                                       //! [0=angles;1=product of all weights]
  // TArrayD *ftaNestedLoopsKine[gKineDependenceVariables][gMaxNoBinsKine][2];
  // //! e-b-e container for nested loops
  // [0=pT,1=eta][kine.bin][0=angles;1=product of all weights]
} nl_a;

// *) Test0:
TList* fTest0List = NULL;        // list to hold all objects for Test0
TProfile* fTest0FlagsPro = NULL; // store all flags for Test0
Bool_t fCalculateTest0 = kFALSE; // calculate or not Test0
struct Test0_Arrays {
  TProfile* fTest0Pro[gMaxCorrelator][gMaxIndex][5] = {
    {{NULL}}}; //! [gMaxCorrelator][gMaxIndex][3]
               //! [order][index][0=integrated,1=vs. multiplicity,2=vs.
               //! centrality,3=pT,4=eta]
  TString* fTest0Labels[gMaxCorrelator][gMaxIndex] = {
    {NULL}}; // all labels: k-p'th order is stored in k-1'th index. So yes, I
             // also store 1-p
} t0_a;
TString fFileWithLabels =
  ""; // path to external ROOT file which specifies all labels of interest
TH1I* fTest0LabelsPlaceholder =
  NULL; // store all Test0 labels in this histogram

// *) Results:
TList* fResultsList = NULL;        //!<! list to hold all results
TProfile* fResultsFlagsPro = NULL; //!<! profile to hold all flags for results
TH1D* fResultsHist = NULL;         //!<! example histogram to store some results

#endif // PWGCF_MULTIPARTICLECORRELATIONS_CORE_MUPA_DATAMEMBERS_H_
