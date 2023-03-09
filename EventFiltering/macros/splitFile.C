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

void splitFile(const char* inputFileName = "AO2D.root", const char* outputFileName1 = "bcRanges.root", const char* outputFileName2 = "cefpDecision.root")
{
  // Open the input ROOT file
  TFile* inputFile = TFile::Open(inputFileName);
  if (!inputFile) {
    std::cerr << "Error: could not open input file " << inputFileName << std::endl;
    return;
  }

  // Open the output ROOT files
  TFile* outputFile1 = TFile::Open(outputFileName1, "RECREATE");
  TFile* outputFile2 = TFile::Open(outputFileName2, "RECREATE");
  if (!outputFile1 || !outputFile2) {
    std::cerr << "Error: could not create output files " << outputFileName1 << " and " << outputFileName2 << std::endl;
    return;
  }

  // Loop over the TDirectories in the input file
  TList* directoryList = inputFile->GetListOfKeys();
  if (!directoryList) {
    std::cerr << "Error: input file has no TDirectory keys" << std::endl;
    return;
  }
  for (int iDir = 0; iDir < directoryList->GetEntries(); ++iDir) {
    TKey* directoryKey = static_cast<TKey*>(directoryList->At(iDir));
    std::cout << "Processing directory " << directoryKey->GetName() << std::endl;

    TDirectoryFile* inputDir = dynamic_cast<TDirectoryFile*>(directoryKey->ReadObj());
    if (!inputDir) {
      TMap* map = dynamic_cast<TMap*>(directoryKey->ReadObj());
      if (!map) {
        std::cerr << "Error: could not read input map " << directoryKey->GetName() << std::endl;
        continue;
      }
      outputFile1->cd();
      map->Write(directoryKey->GetName(), TObject::kSingleKey);
      outputFile2->cd();
      map->Write(directoryKey->GetName(), TObject::kSingleKey);
      continue;
    }

    // Create the output directories in the output files
    TDirectory* outputDir1 = outputFile1->mkdir(inputDir->GetName());
    TDirectory* outputDir2 = outputFile2->mkdir(inputDir->GetName());
    if (!outputDir1 || !outputDir2) {
      std::cerr << "Error: could not create output directories for " << inputDir->GetName() << std::endl;
      continue;
    }
    // Read the trees in the input directory and copy them to the output directories
    TTree* tree1 = static_cast<TTree*>(inputDir->Get("O2bcranges"));
    TTree* tree2 = static_cast<TTree*>(inputDir->Get("O2cefpdecision"));
    if (tree1) {
      outputDir1->cd();
      TTree* newTree1 = tree1->CloneTree();
      newTree1->Write();
    }
    if (tree2) {
      outputDir2->cd();
      TTree* newTree2 = tree2->CloneTree();
      newTree2->Write();
    }
  }

  // Close the input and output ROOT files
  inputFile->Close();
  outputFile1->Close();
  outputFile2->Close();
}
