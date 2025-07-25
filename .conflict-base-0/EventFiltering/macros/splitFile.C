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

#include "TDirectory.h"
#include "TFile.h"
#include "TTree.h"

#include <iostream>
#include <string>

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
    std::cerr << "Error: could not create output files " << outputFileName1 << std::endl;
    return;
  }

  TDirectory* outputDir1 = nullptr;
  TList coll;
  for (auto key : *inputFile->GetListOfKeys()) {

    TDirectoryFile* inputDir = dynamic_cast<TDirectoryFile*>(inputFile->Get(key->GetName()));
    if (!inputDir) {
      continue;
    }
    std::cout << "Processing directory " << inputDir->GetName() << std::endl;
    if (!outputDir1) {
      // Create the output directories in the output files
      outputDir1 = outputFile1->mkdir(inputDir->GetName());
      if (!outputDir1) {
        std::cerr << "Error: could not create output directories for " << inputDir->GetName() << std::endl;
        continue;
      }
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
