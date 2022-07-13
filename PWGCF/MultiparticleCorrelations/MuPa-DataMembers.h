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

// Remarks:
// 0. Starting with C++11, it's possible to initialize data members at declaration, so I do it here
// 1. Use //!<! for introducing a Doxygen comment interpreted as transient in both ROOT 5 and ROOT 6.

// a) Base list to hold all output objects ("grandmother" of all lists);
// *) QA;
// *) Control event histograms;
// *) Particle weights;

// a) Base list to hold all output objects ("grandmother" of all lists):
OutputObj<TList> fBaseList{"Task => ...", OutputObjHandlingPolicy::AnalysisObject, OutputObjSourceType::OutputObjSource};
TProfile* fBasePro = NULL; //!<! keeps flags relevant for the whole analysis
UInt_t fRandomSeed = 0;    // argument to TRandom3 constructor. By default is 0, use SetRandomSeed(...) to change it

// *) QA:
TList* fQAList = NULL; //!<! base list to hold all QA output object

// *) Control event histograms:
TList* fControlEventHistogramsList = NULL;   //!<! list to hold all control event histograms
TProfile* fControlEventHistogramsPro = NULL; //!<! keeps flags relevant for the control event histograms
struct ControlEventHistograms_Arrays {
  TH1D* fMultiplicityHist[2] = {NULL}; //!<! distribution of multiplicity [before,after event cuts]
} ceh_a;                               // "ceh_a" labels an instance of this group of histograms, e.g. ceha_a.fMultiplicityHist[0]

// *) Particle weights:
TList* fWeightsList = NULL;        //!<! list to hold all particle weights
TProfile* fWeightsFlagsPro = NULL; //!<! profile to hold all flags for weights
struct ParticleWeights_Arrays {
  TH1D* fWeightsHist[gWeights] = {NULL}; //!<! particle weights
} pw_a;                                  // "pw_a" labels an instance of this group of histograms, e.g. pw_a.fWeightsHist[0]

// *) Results:
TList* fResultsList = NULL;        //!<! list to hold all results
TProfile* fResultsFlagsPro = NULL; //!<! profile to hold all flags for results
TH1D* fResultsHist = NULL;         //!<! example histogram to store some results
