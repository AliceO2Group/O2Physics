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

#ifndef PWGCF_TWOPARTICLECORRELATIONS_CORE_DIHADRONCONTAINER_H_
#define PWGCF_TWOPARTICLECORRELATIONS_CORE_DIHADRONCONTAINER_H_

/// \file DihadronContainer.h
/// \brief minimum encapsulated container for di-hadron correlations
/// \author Zhiyong Lu (zhiyong.lu@cern.ch)
/// \since  July/2026

#include <Framework/HistogramSpec.h>

#include <TNamed.h>

#include <Rtypes.h>
#include <RtypesCore.h>

#include <cstdint>
#include <vector>

class TH1;
class TH1F;
class TH3;
class TH3F;
class TH2F;
class TH1D;
class TH2;
class TH2D;
class TCollection;
class THnSparse;
class THnBase;
class StepTHn;

class DihadronContainer : public TNamed
{
 public:
  DihadronContainer();
  DihadronContainer(const char* name, const char* objTitle, const std::vector<o2::framework::AxisSpec>& correlationAxis, const uint8_t nStep = 1);
  virtual ~DihadronContainer();

  StepTHn* getCorrHist() { return mCorrHist; }
  void setCorrHist(StepTHn* hist) { mCorrHist = hist; }
  void printCorrHist();

  DihadronContainer(const DihadronContainer& c);
  DihadronContainer& operator=(const DihadronContainer& corr);
  virtual void Copy(TObject& c) const; // NOLINT: Making this override breaks compilation for unknown reason
  void deepCopy(DihadronContainer* from);

  virtual Long64_t Merge(TCollection* list);
  void Reset();
  THnBase* changeToThn(THnBase* sparse);

 protected:
  StepTHn* mCorrHist; // container for n-dimension correlations
  uint8_t nCorrStep;

  ClassDef(DihadronContainer, 2) // underlying event histogram container
};

#endif // PWGCF_TWOPARTICLECORRELATIONS_CORE_DIHADRONCONTAINER_H_
