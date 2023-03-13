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

#include <iostream>
#include <string>
#include "TFile.h"
#include "TDirectory.h"
#include "TTree.h"

void splitFile(const char* inputFileName = "bcSelection.root", const char* outputFileName1 = "bcRanges.root")
{
  // Open the input ROOT file
  TFile* inputFile = TFile::Open(inputFileName);
  if (!inputFile) {
    std::cerr << "Error: could not open input file " << inputFileName << std::endl;
    return;
  }

  // Open the output ROOT files
  TFile* outputFile1 = TFile::Open(outputFileName1, "RECREATE");
  if (!outputFile1) {
    std::cerr << "Error: could not create output files " << outputFileName1 << " and " << outputFileName2 << std::endl;
    return;
  }

  // Loop over the TDirectories in the input file
  TList* directoryList = inputFile->GetListOfKeys();
  if (!directoryList) {
    std::cerr << "Error: input file has no TDirectory keys" << std::endl;
    return;
  }
  bool first = true;
  TDirectory* outputDir1 = nullptr;
  TList coll;
  for (int iDir = 0; iDir < directoryList->GetEntries(); ++iDir) {
    TKey* directoryKey = static_cast<TKey*>(directoryList->At(iDir));

    TDirectoryFile* inputDir = dynamic_cast<TDirectoryFile*>(directoryKey->ReadObj());
    if (!inputDir) {
      continue;
    }

    if (first) {
      // Create the output directories in the output files
      outputDir1 = outputFile1->mkdir(inputDir->GetName());
      if (!outputDir1) {
        std::cerr << "Error: could not create output directories for " << inputDir->GetName() << std::endl;
        continue;
      }
      first = false;
    }
    // Read the trees in the input directory and copy them to the output directories
    TTree* tree1 = dynamic_cast<TTree*>(inputDir->Get("O2bcranges"));
    if (tree1)
      coll.Add(tree1);
  }

  outputDir1->cd();
  TTree* newTree1 = TTree::MergeTrees(&coll);
  newTree1->Write();
  // Close the input and output ROOT files
  inputFile->Close();
  outputFile1->Close();
}
