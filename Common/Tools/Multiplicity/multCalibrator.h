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
#ifndef MULTCALIBRATOR_H
#define MULTCALIBRATOR_H

#include <iostream>
#include "TNamed.h"
#include "TH1D.h"
#include <map>

using namespace std;

class multCalibrator : public TNamed
{

 public:
  //Constructors/Destructor
  multCalibrator();
  multCalibrator(const char* name, const char* title = "Multiplicity Calibration Class");
  ~multCalibrator();

  //void    Print(Option_t *option="") const;

  //_________________________________________________________________________
  //Interface: steering functions to be used in calibration macro

  //Set Filenames
  void SetInputFile(TString lFile) { fInputFileName = lFile.Data(); }
  void SetOutputFile(TString lFile) { fOutputFileName = lFile.Data(); }
  //Set Boundaries to find
  void SetBoundaries(Long_t lNB, Double_t* lB)
  {
    if (lNB < 2 || lNB > 1e+6) {
      cout << "Please make sure you are using a reasonable number of boundaries!" << endl;
      lNB = -1;
    }
    lDesiredBoundaries = lB;
    lNDesiredBoundaries = lNB;
  }

  void SetStandardAdaptiveBoundaries(); //standard adaptive (pp-like)

  //Master Function in this Class: To be called once filenames are set
  Bool_t Calibrate();

  //Auxiliary functions
  Double_t GetBoundaryForPercentile(TH1D* histo, Double_t lPercentileRequested, Double_t& lPrecisionEstimate);

 private:
  //Calibration Boundaries to locate
  Double_t* lDesiredBoundaries;
  Long_t lNDesiredBoundaries;
  Double_t fkPrecisionWarningThreshold;

  TString fInputFileName;  // Filename for TTree object for calibration purposes
  TString fBufferFileName; // Filename for TTree object (buffer file)
  TString fOutputFileName; // Filename for calibration OADB output

  // TList object for storing histograms
  TList* fCalibHists;

  ClassDef(multCalibrator, 1);
  //(this classdef is only for bookkeeping, class will not usually
  // be streamed according to current workflow except in very specific
  // tests!)
};
#endif
