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
#ifndef COMMON_TOOLS_MULTIPLICITY_MULTMCCALIBRATOR_H_
#define COMMON_TOOLS_MULTIPLICITY_MULTMCCALIBRATOR_H_

#include <TF1.h>
#include <TH1.h>
#include <TNamed.h>
#include <TProfile.h>
#include <TString.h>

#include <Rtypes.h>
#include <RtypesCore.h>

class multMCCalibrator : public TNamed
{

 public:
  // Constructors/Destructor
  multMCCalibrator();
  explicit multMCCalibrator(const char* name, const char* title = "MC Multiplicity Calibration Class");
  ~multMCCalibrator();

  //_________________________________________________________________________
  // Interface: steering functions to be used in calibration macro

  // Set Filenames
  void SetDataInputFile(TString lFile) { fDataInputFileName = lFile.Data(); }
  void SetSimInputFile(TString lFile) { fSimInputFileName = lFile.Data(); }
  void SetOutputFile(TString lFile) { fOutputFileName = lFile.Data(); }

  // Master Function in this Class: To be called once filenames are set
  Bool_t Calibrate();

  TF1* GetFit(TProfile* fProf, Bool_t lQuadratic = kTRUE);

  TList* GetCalibList() { return fCalibHists; }

 private:
  TString fDataInputFileName; // Filename for data input
  TString fSimInputFileName;  // Filename for simulation input
  TString fOutputFileName;    // Filename for calibration OADB output

  // TList object for storing histograms
  TList* fCalibHists;

  ClassDef(multMCCalibrator, 1);
  //(this classdef is only for bookkeeping, class will not usually
  // be streamed according to current workflow except in very specific
  // tests!)
};
#endif // COMMON_TOOLS_MULTIPLICITY_MULTMCCALIBRATOR_H_
