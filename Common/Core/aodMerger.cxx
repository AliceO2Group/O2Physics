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

#include <map>
#include <list>
#include <fstream>
#include <getopt.h>

#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TList.h"
#include "TDirectory.h"
#include "TObjString.h"
#include <TGrid.h>
#include <TMap.h>
#include <TLeaf.h>

// AOD merger with correct index rewriting
// No need to know the datamodel because the branch names follow a canonical standard (identified by fIndex)
int main(int argc, char* argv[])
{
  std::string inputCollection("input.txt");
  std::string outputFileName("AO2D.root");
  long maxDirSize = 100000000;
  bool skipNonExistingFiles = false;
  int exitCode = 0; // 0: success, >0: failure

  int option_index = 0;
  static struct option long_options[] = {
    {"input", required_argument, nullptr, 0},
    {"output", required_argument, nullptr, 1},
    {"max-size", required_argument, nullptr, 2},
    {"skip-non-existing-files", no_argument, nullptr, 3},
    {"help", no_argument, nullptr, 4},
    {nullptr, 0, nullptr, 0}};

  while (true) {
    int c = getopt_long(argc, argv, "", long_options, &option_index);
    if (c == -1) {
      break;
    } else if (c == 0) {
      inputCollection = optarg;
    } else if (c == 1) {
      outputFileName = optarg;
    } else if (c == 2) {
      maxDirSize = atol(optarg);
    } else if (c == 3) {
      skipNonExistingFiles = true;
    } else if (c == 4) {
      printf("AOD merging tool. Options: \n");
      printf("  --input <inputfile.txt>      Contains path to files to be merged. Default: %s\n", inputCollection.c_str());
      printf("  --output <outputfile.root>   Target output ROOT file. Default: %s\n", outputFileName.c_str());
      printf("  --max-size <size in Bytes>   Target directory size. Default: %ld\n", maxDirSize);
      printf("  --skip-non-existing-files    Flag to allow skipping of non-existing files in the intput list.\n");
      return -1;
    } else {
      return -2;
    }
  }

  printf("AOD merger started with:\n");
  printf("  Input file: %s\n", inputCollection.c_str());
  printf("  Ouput file name: %s\n", outputFileName.c_str());
  printf("  Maximal folder size (uncompressed): %ld\n", maxDirSize);
  if (skipNonExistingFiles) {
    printf("  WARNING: Skipping non-existing files.\n");
  }

  std::map<std::string, TTree*> trees;
  std::map<std::string, int> offsets;
  std::map<std::string, int> unassignedIndexOffset;

  auto outputFile = TFile::Open(outputFileName.c_str(), "RECREATE", "", 501);
  TDirectory* outputDir = nullptr;
  long currentDirSize = 0;

  std::ifstream in;
  in.open(inputCollection);
  TString line;
  bool connectedToAliEn = false;
  TMap* metaData = nullptr;
  int totalMergedDFs = 0;
  int mergedDFs = 0;
  while (in.good() && exitCode == 0) {
    in >> line;

    if (line.Length() == 0) {
      continue;
    }

    if (line.BeginsWith("alien:") && !connectedToAliEn) {
      printf("Connecting to AliEn...");
      TGrid::Connect("alien:");
      connectedToAliEn = true; // Only try once
    }

    printf("Processing input file: %s\n", line.Data());

    auto inputFile = TFile::Open(line);
    if (!inputFile) {
      printf("Error: Could not open input file %s.\n", line.Data());
      if (skipNonExistingFiles) {
        continue;
      } else {
        printf("Aborting merge!\n");
        exitCode = 1;
        break;
      }
    }

    TList* keyList = inputFile->GetListOfKeys();
    keyList->Sort();

    for (auto key1 : *keyList) {
      if (((TObjString*)key1)->GetString().EqualTo("metaData")) {
        auto metaDataCurrentFile = (TMap*)inputFile->Get("metaData");
        if (metaData == nullptr) {
          metaData = metaDataCurrentFile;
          outputFile->cd();
          metaData->Write("metaData", TObject::kSingleKey);
        } else {
          for (auto metaDataPair : *metaData) {
            auto metaDataKey = ((TPair*)metaDataPair)->Key();
            if (metaDataCurrentFile->Contains(((TObjString*)metaDataKey)->GetString())) {
              auto value = (TObjString*)metaData->GetValue(((TObjString*)metaDataKey)->GetString());
              auto valueCurrentFile = (TObjString*)metaDataCurrentFile->GetValue(((TObjString*)metaDataKey)->GetString());
              if (!value->GetString().EqualTo(valueCurrentFile->GetString())) {
                printf("WARNING: Metadata differs between input files. Key %s : %s vs. %s\n", ((TObjString*)metaDataKey)->GetString().Data(),
                       value->GetString().Data(), valueCurrentFile->GetString().Data());
              }
            } else {
              printf("WARNING: Metadata differs between input files. Key %s is not present in current file\n", ((TObjString*)metaDataKey)->GetString().Data());
            }
          }
        }
      }

      if (!((TObjString*)key1)->GetString().BeginsWith("DF_")) {
        continue;
      }

      auto dfName = ((TObjString*)key1)->GetString().Data();

      printf("  Processing folder %s\n", dfName);
      ++mergedDFs;
      ++totalMergedDFs;
      auto folder = (TDirectoryFile*)inputFile->Get(dfName);
      auto treeList = folder->GetListOfKeys();
      std::list<std::string> foundTrees;

      for (auto key2 : *treeList) {
        auto treeName = ((TObjString*)key2)->GetString().Data();
        foundTrees.push_back(treeName);

        auto inputTree = (TTree*)inputFile->Get(Form("%s/%s", dfName, treeName));
        printf("    Processing tree %s with %lld entries\n", treeName, inputTree->GetEntries());

        if (trees.count(treeName) == 0) {
          if (mergedDFs > 1) {
            printf("    *** FATAL ***: The tree %s was not in the previous dataframe(s)\n", treeName);
            exitCode = 3;
          }

          // clone tree
          // NOTE Basket size etc. are copied in CloneTree()
          if (!outputDir) {
            outputDir = outputFile->mkdir(dfName);
            currentDirSize = 0;
            printf("Writing to output folder %s\n", dfName);
          }
          outputDir->cd();
          auto outputTree = inputTree->CloneTree(-1, "fast");
          outputTree->SetAutoFlush(0);
          trees[treeName] = outputTree;
          currentDirSize += inputTree->GetTotBytes();
        } else {
          // append tree
          auto outputTree = trees[treeName];

          outputTree->CopyAddresses(inputTree);

          // register index and connect VLA columns
          std::vector<std::pair<int*, int>> indexList;
          std::vector<char*> vlaPointers;
          TObjArray* branches = inputTree->GetListOfBranches();
          for (int i = 0; i < branches->GetEntriesFast(); ++i) {
            TBranch* br = (TBranch*)branches->UncheckedAt(i);
            TString branchName(br->GetName());

            // detect VLA
            if (((TLeaf*)br->GetListOfLeaves()->First())->GetLeafCount() != nullptr) {
              int maximum = ((TLeaf*)br->GetListOfLeaves()->First())->GetLeafCount()->GetMaximum();
              char* buffer = new char[maximum * 8]; // assume 64 bit as largest case
              printf("      Allocated VLA buffer of size %d for branch name %s\n", maximum, br->GetName());
              inputTree->SetBranchAddress(br->GetName(), buffer);
              outputTree->SetBranchAddress(br->GetName(), buffer);
              vlaPointers.push_back(buffer);
            }

            if (branchName.BeginsWith("fIndex")) {
              // Syntax: fIndex<Table>[_<Suffix>]
              branchName.Remove(0, 6);
              if (branchName.First("_") > 0) {
                branchName.Remove(branchName.First("_"));
              }
              branchName.Remove(branchName.Length() - 1); // remove s
              branchName.ToLower();
              branchName = "O2" + branchName;

              indexList.push_back({new int, offsets[branchName.Data()]});

              inputTree->SetBranchAddress(br->GetName(), indexList.back().first);
              outputTree->SetBranchAddress(br->GetName(), indexList.back().first);
            }
          }

          // on the first appending pass we need to find out the most negative index in the existing output
          // to correctly continue negative index assignment
          if (mergedDFs == 2) {
            int minIndex = -1;
            if (indexList.size() > 0) {
              outputTree->SetBranchStatus("*", 0);
              outputTree->SetBranchStatus("fIndex*", 1);
              auto outentries = outputTree->GetEntries();
              for (int i = 0; i < outentries; ++i) {
                outputTree->GetEntry(i);
                for (const auto& idx : indexList) {
                  minIndex = std::min(*(idx.first), minIndex);
                }
              }
              outputTree->SetBranchStatus("*", 1);
            }
            unassignedIndexOffset[treeName] = minIndex;
          }

          auto entries = inputTree->GetEntries();
          int minIndexOffset = unassignedIndexOffset[treeName];
          auto newMinIndexOffset = minIndexOffset;
          for (int i = 0; i < entries; i++) {
            inputTree->GetEntry(i);
            // shift index columns by offset
            for (const auto& idx : indexList) {
              // if negative, the index is unassigned. In this case, the different unassigned blocks have to get unique negative IDs
              if (*(idx.first) < 0) {
                *(idx.first) += minIndexOffset;
                newMinIndexOffset = std::min(newMinIndexOffset, *(idx.first));
              } else {
                *(idx.first) += idx.second;
              }
            }
            int nbytes = outputTree->Fill();
            if (nbytes > 0) {
              currentDirSize += nbytes;
            }
          }
          unassignedIndexOffset[treeName] = newMinIndexOffset;

          delete inputTree;

          for (const auto& idx : indexList) {
            delete idx.first;
          }
          for (auto& buffer : vlaPointers) {
            delete[] buffer;
          }
        }
      }
      if (exitCode > 0) {
        break;
      }

      // check if all trees were present
      if (mergedDFs > 1) {
        for (auto const& tree : trees) {
          bool found = (std::find(foundTrees.begin(), foundTrees.end(), tree.first) != foundTrees.end());
          if (found == false) {
            printf("  *** FATAL ***: The tree %s was not in the current dataframe\n", tree.first.c_str());
            exitCode = 4;
          }
        }
      }

      // update offsets
      for (auto const& tree : trees) {
        offsets[tree.first] = tree.second->GetEntries();
      }

      // check for not found tables
      for (auto const& offset : offsets) {
        if (trees.count(offset.first) == 0) {
          printf("ERROR: Index on %s but no tree found\n", offset.first.c_str());
        }
      }

      if (currentDirSize > maxDirSize) {
        printf("Maximum size reached: %ld. Closing folder %s.\n", currentDirSize, dfName);
        for (auto const& tree : trees) {
          // printf("Writing %s\n", tree.first.c_str());
          outputDir->cd();
          tree.second->Write();
          delete tree.second;
        }
        outputDir = nullptr;
        trees.clear();
        offsets.clear();
        mergedDFs = 0;
      }
    }
    inputFile->Close();
  }

  outputFile->Write();
  outputFile->Close();

  if (totalMergedDFs == 0) {
    printf("ERROR: Did not merge a single DF. This does not seem right.\n");
    exitCode = 2;
  }

  // in case of failure, remove the incomplete file
  if (exitCode != 0) {
    printf("Removing incomplete output file %s.\n", outputFile->GetName());
    gSystem->Unlink(outputFile->GetName());
  }

  printf("AOD merger finished.\n");

  return exitCode;
}
