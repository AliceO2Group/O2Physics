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
//
// This code calculates a centrality calibration based on output
// from the multiplicityQa task.
//
// Comments, suggestions, questions? Please write to:
// - victor.gonzalez@cern.ch
// - david.dobrigkeit.chinellato@cern.ch
//
#ifndef COMMON_TOOLS_MULTIPLICITY_MULTCALIBRATOR_H_
#define COMMON_TOOLS_MULTIPLICITY_MULTCALIBRATOR_H_

#include <TH1.h>
#include <TNamed.h>
#include <TString.h>

#include <Rtypes.h>
#include <RtypesCore.h>

#include <iostream>

class multCalibrator : public TNamed
{

 public:
  // Constructors/Destructor
  multCalibrator();
  explicit multCalibrator(const char* name, const char* title = "Multiplicity Calibration Class");
  ~multCalibrator();

  // void    Print(Option_t *option="") const;

  //_________________________________________________________________________
  // Interface: steering functions to be used in calibration macro

  // Set Filenames
  void SetInputFile(TString lFile) { fInputFileName = lFile.Data(); }
  void SetOutputFile(TString lFile) { fOutputFileName = lFile.Data(); }
  // Set Boundaries to find
  void SetBoundaries(Long_t lNB, Double_t* lB)
  {
    if (lNB < 2 || lNB > 1e+6) {
      std::cout << "Please make sure you are using a reasonable number of boundaries!" << std::endl;
      lNB = -1;
    }
    lDesiredBoundaries = lB;
    lNDesiredBoundaries = lNB;
  }

  void SetAnchorPointRaw(Float_t lRaw) { fAnchorPointValue = lRaw; }
  void SetAnchorPointPercentage(Float_t lPer) { fAnchorPointPercentage = lPer; }

  void SetStandardAdaptiveBoundaries();   // standard adaptive (pp-like)
  void SetStandardOnePercentBoundaries(); // standard 1% (Pb-Pb like)

  // Master Function in this Class: To be called once filenames are set
  Bool_t Calibrate();

  // Aux function. Keep public, accessible outside as rather useful utility
  TH1F* GetCalibrationHistogram(TH1* histoRaw, TString lHistoName = "hCalib");

  // Auxiliary functions
  Double_t GetRawMax(TH1* histo);
  Double_t GetBoundaryForPercentile(TH1* histo, Double_t lPercentileRequested, Double_t& lPrecisionEstimate);

  // Precision bookkeeping
  TH1D* GetPrecisionHistogram() { return fPrecisionHistogram; } // gets precision histogram from current object
  void ResetPrecisionHistogram();                               // Reset precision histogram, if it exists

  // Aliases for centrality estimators
  enum fCentEstim {
    kCentRawV0M = 0,
    kCentRawT0M,
    kCentRawFDD,
    kCentRawNTracks,
    kCentZeqV0M,
    kCentZeqT0M,
    kCentZeqFDD,
    kCentZeqNTracks,
    kNCentEstim
  };

  static const TString fCentEstimName[kNCentEstim]; //! name (internal)

 private:
  // Calibration Boundaries to locate
  Double_t* lDesiredBoundaries;
  Long_t lNDesiredBoundaries;
  Double_t fkPrecisionWarningThreshold;

  TString fInputFileName;  // Filename for TTree object for calibration purposes
  TString fBufferFileName; // Filename for TTree object (buffer file)
  TString fOutputFileName; // Filename for calibration OADB output

  // Anchor point functionality
  Float_t fAnchorPointValue;      // AP value (raw estimator)
  Float_t fAnchorPointPercentage; // AP percentage

  // TList object for storing histograms
  TList* fCalibHists;

  TH1D* fPrecisionHistogram; // for bookkeeping of precision report

  ClassDef(multCalibrator, 1);
  //(this classdef is only for bookkeeping, class will not usually
  // be streamed according to current workflow except in very specific
  // tests!)
};
#endif // COMMON_TOOLS_MULTIPLICITY_MULTCALIBRATOR_H_
