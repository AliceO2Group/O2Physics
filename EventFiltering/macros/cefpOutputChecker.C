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
// O2 includes

#include <TFile.h>
#include <TTree.h>

void cefpOutputChecker(std::string histoFile = "AnalysisResults.root", std::string treeFile = "AO2D.root")
{
  TFile referenceFile(histoFile.data());
  TH1* refScalers = (TH1*)referenceFile.Get("central-event-filter-task/scalers/mScalers");
  TH1* refFilters = (TH1*)referenceFile.Get("central-event-filter-task/scalers/mFiltered");

  TFile outputFile("cefpOutputChecker.root", "recreate");
  TH1* newFilters = (TH1*)refFilters->Clone("newFilters");
  TH1* newScalers = (TH1*)refScalers->Clone("newScalers");
  newFilters->Reset();
  newScalers->Reset();

  TFile inputFile(treeFile.data());
  for (auto key : *(inputFile.GetListOfKeys())) {
    TTree* cefpTree = (TTree*)inputFile.Get(Form("%s/O2cefpdecision", key->GetName()));
    if (!cefpTree)
      continue;
    ULong64_t fCefpSelected, fCefpTriggered;
    cefpTree->SetBranchAddress("fCefpSelected", &fCefpSelected);
    cefpTree->SetBranchAddress("fCefpTriggered", &fCefpTriggered);
    for (int i = 0; i < cefpTree->GetEntries(); i++) {
      cefpTree->GetEntry(i);
      newFilters->Fill(0);
      newScalers->Fill(0);
      for (ULong64_t j = 0; j < 64; j++) {
        if (fCefpSelected & (1ull << j))
          newFilters->Fill(j + 1);
        if (fCefpTriggered & (1ull << j))
          newScalers->Fill(j + 1);
      }
      if (fCefpSelected)
        newFilters->Fill(newFilters->GetNbinsX() - 1);
      if (fCefpTriggered)
        newScalers->Fill(newScalers->GetNbinsX() - 1);
    }
  }

  outputFile.cd();
  newFilters->Write();
  newScalers->Write();
  newScalers->Divide(refScalers);
  newFilters->Divide(refFilters);
  newFilters->Write("ratioFilters");
  newScalers->Write("ratioScalers");
}
