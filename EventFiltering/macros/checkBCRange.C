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

#include "CommonDataFormat/InteractionRecord.h"
#include "CommonDataFormat/IRFrame.h"

using o2::InteractionRecord;
using o2::dataformats::IRFrame;

void checkBCRange(const char* filename = "AO2D.root")
{

  // Open the ROOT file and get the trees
  TFile* inputFile = new TFile(filename, "READ");

  // Loop over the TDirectories in the input file
  TList* directoryList = inputFile->GetListOfKeys();
  if (!directoryList) {
    std::cerr << "Error: input file has no TDirectory keys" << std::endl;
    return;
  }
  for (int iDir = 0; iDir < directoryList->GetEntries(); ++iDir) {
    TKey* directoryKey = static_cast<TKey*>(directoryList->At(iDir));

    TDirectoryFile* inputDir = dynamic_cast<TDirectoryFile*>(directoryKey->ReadObj());
    if (!inputDir) {
      continue;
    }
    // std::cout << "Processing directory " << directoryKey->GetName() << std::endl;
    auto dirName = directoryKey->GetName();
    TTree* treeRanges = dynamic_cast<TTree*>(inputFile->Get(Form("%s/O2bcranges", dirName)));
    TTree* treeDecision = dynamic_cast<TTree*>(inputFile->Get(Form("%s/O2cefpdecision", dirName)));

    // Get the branches we need from the trees
    ULong64_t bcstart, bcend, globalBCId, cefpSelected = 0;
    treeRanges->SetBranchAddress("fBCstart", &bcstart);
    treeRanges->SetBranchAddress("fBCend", &bcend);
    treeDecision->SetBranchAddress("fGlobalBCId", &globalBCId);
    treeDecision->SetBranchAddress("fCefpSelected", &cefpSelected);

    std::vector<InteractionRecord> bcids;
    // Loop over the entries in the decision tree and check if the BC range is valid
    int nEntriesDecision = treeDecision->GetEntries();
    for (int iEntryDecision = 0; iEntryDecision < nEntriesDecision; ++iEntryDecision) {
      treeDecision->GetEntry(iEntryDecision);
      if (cefpSelected == 0)
        continue;
      InteractionRecord ir;
      ir.setFromLong(globalBCId);
      bcids.push_back(ir);
    }
    std::vector<bool> found(bcids.size(), false);

    // Loop over the entries in the ranges tree and check if the BC range is valid
    int nEntriesRanges = treeRanges->GetEntries();
    for (int iEntryRanges = 0; iEntryRanges < nEntriesRanges; ++iEntryRanges) {
      treeRanges->GetEntry(iEntryRanges);
      InteractionRecord irstart, irend;
      irstart.setFromLong(bcstart);
      irend.setFromLong(bcend);
      IRFrame frame(irstart, irend);
      for (int i = 0; i < bcids.size(); i++) {
        auto& bcid = bcids[i];
        if (!frame.isOutside(bcid)) {
          found[i] = true;
        }
      }
    }
    int notFound = 0;
    for (int i = 0; i < bcids.size(); i++) {
      notFound += !found[i];
    }
    std::cout << "Found " << notFound << " BCs not in ranges out of " << bcids.size() << std::endl;
  }
}