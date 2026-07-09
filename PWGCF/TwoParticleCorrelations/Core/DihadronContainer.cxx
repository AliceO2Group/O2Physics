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

/// \file DihadronContainer.h
/// \brief minimum encapsulated container for di-hadron correlations
/// \author Zhiyong Lu (zhiyong.lu@cern.ch)
/// \since  July/2026

#include "PWGCF/TwoParticleCorrelations/Core/DihadronContainer.h"

#include <Framework/HistogramSpec.h>
#include <Framework/Logger.h>
#include <Framework/StepTHn.h>

#include <TCollection.h>
#include <THn.h>
#include <TIterator.h>
#include <TList.h>
#include <TNamed.h>
#include <TObject.h>
#include <TString.h>

#include <Rtypes.h>
#include <RtypesCore.h>

#include <cstdint>
#include <cstring>
#include <vector>

using namespace o2;
using namespace o2::framework;

ClassImp(DihadronContainer);

DihadronContainer::DihadronContainer() : mCorrHist(nullptr), nCorrStep(1)
{
  // Default constructor
}

DihadronContainer::DihadronContainer(const char* name, const char* objTitle,
                                     const std::vector<o2::framework::AxisSpec>& correlationAxis,
                                     const uint8_t nStep) : TNamed(name, objTitle),
                                                            mCorrHist(nullptr),
                                                            nCorrStep(nStep)
{

  if (strlen(name) == 0) {
    return;
  }

  int64_t bins = 1;
  LOGF(info, "Creating DihadronContainer:");
  for (uint iAxis = 0; iAxis < correlationAxis.size(); iAxis++) {
    LOGF(info, "correlationAxis[%d].getNbins() = %d", iAxis, correlationAxis[iAxis].getNbins());
    bins *= correlationAxis[iAxis].getNbins();
  }
  LOGF(info, "DihadronContainer with %ld bins in the correlation histogram (approx. %ld-%ld MB of memory)", bins, bins * 4 / 1024 / 1024, bins * 8 / 1024 / 1024);

  mCorrHist = HistFactory::createHist<StepTHnF>({"mCorrHist", "d^{2}N_{ch}/d#varphid#eta", {HistType::kStepTHnF, correlationAxis, nStep}}).release();
}

//_____________________________________________________________________________
DihadronContainer::DihadronContainer(const DihadronContainer& c) : TNamed(c),
                                                                   mCorrHist(nullptr),
                                                                   nCorrStep(1)
{
  // DihadronContainer copy constructor
  static_cast<const DihadronContainer&>(c).Copy(*this);
}

//____________________________________________________________________
DihadronContainer::~DihadronContainer()
{
  // Destructor
  if (mCorrHist) {
    delete mCorrHist;
    mCorrHist = nullptr;
  }
}

//____________________________________________________________________
DihadronContainer& DihadronContainer::operator=(const DihadronContainer& c)
{
  // assigment operator
  if (this != &c) {
    static_cast<const DihadronContainer&>(c).Copy(*this);
  }
  return *this;
}

//____________________________________________________________________
void DihadronContainer::Copy(TObject& c) const
{
  // copy function
  auto target = dynamic_cast<DihadronContainer&>(c);

  if (mCorrHist) {
    target.mCorrHist = dynamic_cast<StepTHn*>(mCorrHist->Clone());
  }

  target.nCorrStep = nCorrStep;
}

//____________________________________________________________________
void DihadronContainer::deepCopy(DihadronContainer* from)
{
  // copies the entries of this object's members from the object <from> to this object
  // fills using the fill function and thus allows that the objects have different binning

  for (Int_t step = 0; step < mCorrHist->getNSteps(); step++) {
    LOGF(info, "Copying step %d", step);
    THnBase* target = mCorrHist->getTHn(step);
    THnBase* source = from->mCorrHist->getTHn(step);

    target->Reset();
    target->RebinnedAdd(source);
  }
}

//____________________________________________________________________
Long64_t DihadronContainer::Merge(TCollection* list)
{
  // Merge a list of DihadronContainer objects with this
  // Returns the number of merged objects (including this).

  if (!list) {
    return 0;
  }

  if (list->IsEmpty()) {
    return 1;
  }

  TIterator* iter = list->MakeIterator();
  TObject* obj = nullptr;

  // collections of objects
  const UInt_t kMaxLists = 1;
  auto lists = new TList*[kMaxLists];

  for (UInt_t i = 0; i < kMaxLists; i++) {
    lists[i] = new TList;
  }

  Int_t count = 0;
  while ((obj = iter->Next())) {

    auto entry = dynamic_cast<DihadronContainer*>(obj);
    if (entry == nullptr) {
      continue;
    }

    if (entry->mCorrHist) {
      lists[0]->Add(entry->mCorrHist);
    }

    count++;
  }

  if (mCorrHist) {
    mCorrHist->Merge(lists[0]);
  }

  for (UInt_t i = 0; i < kMaxLists; i++) {
    delete lists[i];
  }

  delete[] lists;

  return count + 1;
}

//____________________________________________________________________
void DihadronContainer::Reset()
{
  // resets all contained histograms

  for (Int_t step = 0; step < mCorrHist->getNSteps(); step++) {
    mCorrHist->getTHn(step)->Reset();
  }
}

THnBase* DihadronContainer::changeToThn(THnBase* sparse)
{
  // change the object to THn for faster processing

  return THn::CreateHn(Form("%s_thn", sparse->GetName()), sparse->GetTitle(), sparse);
}

ClassImp(DihadronContainer);
