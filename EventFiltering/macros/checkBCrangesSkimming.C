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
#include <TH1.h>
#include <TTree.h>
#include <cmath>
#include <vector>
#include <iostream>
#include "CommonDataFormat/InteractionRecord.h"
#include "CommonDataFormat/IRFrame.h"

using o2::InteractionRecord;
using o2::dataformats::IRFrame;

const ULong64_t trigger0Bit = 54;
const ULong64_t trigger1Bit = 2;
const int bcDiffTolerance = 0;

struct bcTuple {
  bcTuple(ULong64_t bcAO2D, ULong64_t bcEvSel) : bcAO2D(bcAO2D), bcEvSel(bcEvSel) {}
  ULong64_t bcAO2D{0ull};
  ULong64_t bcEvSel{0ull};
  bool operator==(const bcTuple& t) const
  {
    return (this->bcAO2D == t.bcAO2D && this->bcEvSel == t.bcEvSel);
  }
};

struct selectedFrames : public IRFrame {
  selectedFrames(ULong64_t bcAO2D, ULong64_t bcEvSel, ULong64_t triMask[2], ULong64_t selMask[2], const IRFrame& frame) : IRFrame(frame), bcAO2D(bcAO2D), bcEvSel(bcEvSel), triMask{triMask[0], triMask[1]}, selMask{selMask[0], selMask[1]} {}
  ULong64_t triMask[2]{0ull}, selMask[2]{0ull}, bcAO2D, bcEvSel;
};

std::vector<selectedFrames> getSelectedFrames(TFile& file)
{
  std::vector<selectedFrames> selectedFrames;
  ULong64_t bcAO2D{0ull}, bcEvSel{0ull}, triMask[2]{0ull}, selMask[2]{0ull};
  for (auto key : *file.GetListOfKeys()) {
    auto dir = dynamic_cast<TDirectory*>(file.Get(key->GetName()));
    if (!dir) {
      continue;
    }
    auto tree = dynamic_cast<TTree*>(dir->Get("selectedBC"));
    if (!tree) {
      continue;
    }
    tree->SetBranchAddress("bcAO2D", &bcAO2D);
    tree->SetBranchAddress("bcEvSel", &bcEvSel);
    if (tree->GetBranch("triMask")) {
      tree->SetBranchAddress("triMask", &triMask[0]);
      tree->SetBranchAddress("selMask", &selMask[0]);
    } else {
      tree->SetBranchAddress("triMask0", &triMask[0]);
      tree->SetBranchAddress("triMask1", &triMask[1]);
      tree->SetBranchAddress("selMask0", &selMask[0]);
      tree->SetBranchAddress("selMask1", &selMask[1]);
    }

    for (int i = 0; i < tree->GetEntries(); i++) {
      tree->GetEntry(i);
      if (!selMask[0] && !selMask[1]) {
        continue;
      }
      if (selMask[0] & BIT(trigger0Bit) || selMask[1] & BIT(trigger1Bit)) {
        InteractionRecord irstart, irend;
        irstart.setFromLong(std::min(bcAO2D, bcEvSel));
        irend.setFromLong(std::max(bcAO2D, bcEvSel));
        IRFrame frame(irstart, irend);
        selectedFrames.push_back({bcAO2D, bcEvSel, triMask, selMask, frame});
      }
    }
    std::cout << "TreeSize: " << tree->GetEntries() << std::endl;
  }
  std::cout << "Size: " << selectedFrames.size() << std::endl;
  return selectedFrames;
}

void checkBCrangesSkimming(std::string originalFileName = "bcRanges_fullrun.root", std::string skimmedFileName = "bcRanges_fullrun_skimmed.root")
{
  TFile originalFile(originalFileName.c_str(), "READ");
  TFile skimmedFile(skimmedFileName.c_str(), "READ");
  TH1F* hTriggerCounter = new TH1F("hTriggerCounter", "hTriggerCounter", 3, 0.5, 3.5);
  hTriggerCounter->GetXaxis()->SetBinLabel(1, "Original");
  hTriggerCounter->GetXaxis()->SetBinLabel(2, "Skimmed");
  TH1F* hNumCounter = new TH1F("hNumCounter", "hTriggerCounter", 10, -0.5, 9.5);
  TH1F* hPairedTriggerCounter = new TH1F("hPairedTriggerCounter", "hPairedTriggerCounter", 4, 0.5, 4.5);
  hPairedTriggerCounter->GetXaxis()->SetBinLabel(1, "Total");
  hPairedTriggerCounter->GetXaxis()->SetBinLabel(2, "Same AO2D BC");
  hPairedTriggerCounter->GetXaxis()->SetBinLabel(3, "Same EvSel BC");
  hPairedTriggerCounter->GetXaxis()->SetBinLabel(4, "Same Both BC");
  TH1F* hSinglePairCheck = new TH1F("hSinglePairCheck", "hSinglePairCheck", 4, 0.5, 4.5);
  hSinglePairCheck->GetXaxis()->SetBinLabel(1, "Total");
  hSinglePairCheck->GetXaxis()->SetBinLabel(2, "Same AO2D BC");
  hSinglePairCheck->GetXaxis()->SetBinLabel(3, "Same EvSel BC");
  hSinglePairCheck->GetXaxis()->SetBinLabel(4, "Same Both BC");
  TH1F* hMultiPairCheck = new TH1F("hMultiPairCheck", "hMultiPairCheck", 5, 0.5, 5.5);
  hMultiPairCheck->GetXaxis()->SetBinLabel(1, "Total");
  hMultiPairCheck->GetXaxis()->SetBinLabel(2, "Total Pair");
  hMultiPairCheck->GetXaxis()->SetBinLabel(3, "Same AO2D BC");
  hMultiPairCheck->GetXaxis()->SetBinLabel(4, "Same EvSel BC");
  hMultiPairCheck->GetXaxis()->SetBinLabel(5, "Same Both BC");
  TH1F* hBCDiffAO2D = new TH1F("hBCDiffAO2D", "hBCDiffAO2D", 1000, -5.e3, 5.e3);
  TH1F* hBCDiffEvSel = new TH1F("hBCDiffEvSel", "hBCDiffEvSel", 1000, -5.e3, 5.e3);

  TH1F* hBCOriginal = new TH1F("hBCOriginal", "hBCOriginal", 4, 0.5, 4.5);
  hBCOriginal->GetXaxis()->SetBinLabel(1, "Total");
  hBCOriginal->GetXaxis()->SetBinLabel(2, "Same AO2D BC");
  hBCOriginal->GetXaxis()->SetBinLabel(3, "Same EvSel BC");
  hBCOriginal->GetXaxis()->SetBinLabel(4, "Same Both BC");
  TH1F* hBCSkimmed = new TH1F("hBCSkimmed", "hBCSkimmed", 4, 0.5, 4.5);
  hBCSkimmed->GetXaxis()->SetBinLabel(1, "Total");
  hBCSkimmed->GetXaxis()->SetBinLabel(2, "Same AO2D BC");
  hBCSkimmed->GetXaxis()->SetBinLabel(3, "Same EvSel BC");
  hBCSkimmed->GetXaxis()->SetBinLabel(4, "Same Both BC");

  auto t1 = std::chrono::steady_clock::now();
  auto originalFrames = getSelectedFrames(originalFile);
  auto skimmedFrames = getSelectedFrames(skimmedFile);
  auto t2 = std::chrono::steady_clock::now();
  int d1 = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
  std::cout << "Readin Time: " << d1 << std::endl;

  auto t3 = std::chrono::steady_clock::now();
  std::sort(originalFrames.begin(), originalFrames.end(), [](const selectedFrames& a, const selectedFrames& b) {
    return a.getMin() < b.getMin();
  });
  std::sort(skimmedFrames.begin(), skimmedFrames.end(), [](const selectedFrames& a, const selectedFrames& b) {
    return a.getMin() < b.getMin();
  });
  auto t4 = std::chrono::steady_clock::now();
  int d2 = std::chrono::duration_cast<std::chrono::milliseconds>(t4 - t3).count();
  std::cout << "Sort Time: " << d2 << std::endl;

  auto t5 = std::chrono::steady_clock::now();
  std::vector<bcTuple> bcSet;
  std::vector<ULong64_t> bcAO2DSet;
  std::vector<ULong64_t> bcEvSelSet;
  for (auto frame : originalFrames) {
    if (frame.selMask[0] & BIT(trigger0Bit)) {
      bool found = false;
      hTriggerCounter->Fill(1);
      hBCOriginal->Fill(1);
      auto p1 = std::find_if(bcSet.begin(), bcSet.end(), [&](const auto& val) { return val.bcAO2D == frame.bcAO2D; });
      if (p1 != bcSet.end()) {
        hBCOriginal->Fill(2);
      }
      auto p2 = std::find_if(bcSet.begin(), bcSet.end(), [&](const auto& val) { return val.bcEvSel == frame.bcEvSel; });
      if (p2 != bcSet.end()) {
        hBCOriginal->Fill(3);
      }
      bcTuple currentBC(frame.bcAO2D, frame.bcEvSel);
      auto p3 = std::find(bcSet.begin(), bcSet.end(), currentBC);
      if (p3 == bcSet.end()) {
        bcSet.push_back(currentBC);
      } else {
        hBCOriginal->Fill(4);
      }
      // std::cout << "------------------------------------------------" << std::endl;
      frame.getMin() -= bcDiffTolerance;
      frame.getMax() += bcDiffTolerance;
      std::vector<bcTuple> skimmedbcs;
      int n = 0;
      for (auto& skimmedFrame : skimmedFrames) {
        if (skimmedFrame.getMin() > frame.getMax()) {
          break;
        }
        if (!frame.isOutside(skimmedFrame)) {
          found = frame.selMask[0] & skimmedFrame.selMask[0] || frame.selMask[1] & skimmedFrame.selMask[1];
          found = found && (frame.bcAO2D == skimmedFrame.bcAO2D || frame.bcEvSel == skimmedFrame.bcEvSel);
          if (found) {
            hPairedTriggerCounter->Fill(1);
            if (frame.bcAO2D == skimmedFrame.bcAO2D) {
              hPairedTriggerCounter->Fill(2);
            }
            if (frame.bcEvSel == skimmedFrame.bcEvSel) {
              hPairedTriggerCounter->Fill(3);
              if (frame.bcAO2D == skimmedFrame.bcAO2D) {
                hPairedTriggerCounter->Fill(4);
              }
            }
            skimmedbcs.push_back({skimmedFrame.bcAO2D, skimmedFrame.bcEvSel});
            n++;
          }
        }
      }
      if (n == 0) {
        // std::cout << "Trigger not found!!!   " << n << std::endl;
      } else if (n == 1) {
        hSinglePairCheck->Fill(1);
        hBCDiffAO2D->Fill(frame.bcAO2D - skimmedbcs[0].bcAO2D);
        hBCDiffEvSel->Fill(frame.bcEvSel - skimmedbcs[0].bcEvSel);
        if (frame.bcAO2D == skimmedbcs[0].bcAO2D) {
          hSinglePairCheck->Fill(2);
        }
        if (frame.bcEvSel == skimmedbcs[0].bcEvSel) {
          hSinglePairCheck->Fill(3);
          if (frame.bcAO2D == skimmedbcs[0].bcAO2D) {
            hSinglePairCheck->Fill(4);
          }
        }
      } else {
        // std::cout << "Unexpected trigger!!!   " << n << std::endl;
        hMultiPairCheck->Fill(1);
        for (auto skimmedbc : skimmedbcs) {
          hMultiPairCheck->Fill(2);
          if (frame.bcAO2D == skimmedbc.bcAO2D) {
            hMultiPairCheck->Fill(3);
          }
          if (frame.bcEvSel == skimmedbc.bcEvSel) {
            hMultiPairCheck->Fill(4);
            if (frame.bcAO2D == skimmedbc.bcAO2D) {
              hMultiPairCheck->Fill(5);
            }
          }
        }
      }
      hNumCounter->Fill(n);
    }
  }
  auto t6 = std::chrono::steady_clock::now();
  int d3 = std::chrono::duration_cast<std::chrono::milliseconds>(t6 - t5).count();
  std::cout << "Search Time: " << d3 << std::endl;

  bcSet.clear();
  for (auto& skimmedFrame : skimmedFrames) {
    if (skimmedFrame.selMask[0] & BIT(trigger0Bit) || skimmedFrame.selMask[1] & (trigger1Bit)) {
      hTriggerCounter->Fill(2);
      hBCSkimmed->Fill(1);
      auto p1 = std::find_if(bcSet.begin(), bcSet.end(), [&](const auto& val) { return val.bcAO2D == skimmedFrame.bcAO2D; });
      if (p1 != bcSet.end()) {
        hBCSkimmed->Fill(2);
      }
      auto p2 = std::find_if(bcSet.begin(), bcSet.end(), [&](const auto& val) { return val.bcEvSel == skimmedFrame.bcEvSel; });
      if (p2 != bcSet.end()) {
        hBCSkimmed->Fill(3);
      }
      bcTuple currentBC(skimmedFrame.bcAO2D, skimmedFrame.bcEvSel);
      auto p3 = std::find(bcSet.begin(), bcSet.end(), currentBC);
      if (p3 == bcSet.end()) {
        bcSet.push_back(currentBC);
      } else {
        hBCSkimmed->Fill(4);
      }
    }
  }

  TFile fout("check.root", "RECREATE");
  fout.cd();
  hTriggerCounter->Write();
  hBCOriginal->Write();
  hBCSkimmed->Write();
  hNumCounter->Write();
  hPairedTriggerCounter->Write();
  hSinglePairCheck->Write();
  hBCDiffAO2D->Write();
  hBCDiffEvSel->Write();
  hMultiPairCheck->Write();
  fout.Close();
}
