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
#include <regex>
#include "CommonDataFormat/InteractionRecord.h"
#include "CommonDataFormat/IRFrame.h"

using o2::InteractionRecord;
using o2::dataformats::IRFrame;

// Set the bit of trigger which need to be checked
const ULong64_t trigger0Bit = BIT(54);
const ULong64_t trigger1Bit = 0;
const int bcDiffTolerance = 0;
const char outputFileName[15] = "output.root";
bool skipNoneDuplicate = false;

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

std::vector<selectedFrames> getSelectedFrames(TFile& file, ULong64_t trigger0Bit, ULong64_t trigger1Bit)
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
      if (selMask[0] & trigger0Bit || selMask[1] & trigger1Bit) {
        InteractionRecord irstart, irend;
        irstart.setFromLong(std::min(bcAO2D, bcEvSel));
        irend.setFromLong(std::max(bcAO2D, bcEvSel));
        IRFrame frame(irstart, irend);
        selectedFrames.push_back({bcAO2D, bcEvSel, triMask, selMask, frame});
      }
    }
  }
  return selectedFrames;
}

// Calulate the ratio of duplicate triggers
void checkDuplicateTrigger(std::string AnaFileName = "AnalysisResults.root", std::string originalFileName = "bcRanges_fullrun.root", std::string skimmedFileName = "bcRanges_fullrun_skimmed.root")
{

  std::string runNumber = "";
  std::regex re("/5[0-9]*");
  std::smatch match;
  if (std::regex_search(originalFileName, match, re)) {
    // Remove the leading '/'
    runNumber = match.str().substr(1);
  }

  // Readin labels
  TFile AnaFile(AnaFileName.c_str(), "READ");
  TH1* hist0 = dynamic_cast<TH1*>(AnaFile.Get("central-event-filter-task/scalers/mFiltered;1"));
  std::vector<std::string> labels;
  std::vector<int> binNum;
  for (int i = 1; i <= hist0->GetNbinsX(); i++) {
    std::string label = hist0->GetXaxis()->GetBinLabel(i);
    if (label != "Total number of events" && label != "Filtered events") {
      labels.push_back(label);
      binNum.push_back(i);
    }
  }
  AnaFile.Close();

  TFile originalFile(originalFileName.c_str(), "READ");
  TFile skimmedFile(skimmedFileName.c_str(), "READ");
  std::vector<std::string> sel_labels;
  std::vector<double> numOriginal, numSkimmed, numOriginalDuplicate, numSkimmedDuplicate;
  for (int i = 0; i < labels.size(); i++) {
    ULong64_t trigger0Bit = 0, trigger1Bit = 0;
    int triggerBit = binNum[i] - 2;
    if (triggerBit < 64) {
      trigger0Bit = BIT(triggerBit);
    } else {
      trigger1Bit = BIT(triggerBit - 64);
    }
    // For Original dataset
    std::vector<bcTuple> bcSet;
    double noriginal{0}, nskimmed{0}, noriginalduplicate{0}, nskimmedduplicate{0};
    auto Frames = getSelectedFrames(originalFile, trigger0Bit, trigger1Bit);
    for (auto frame : Frames) {
      noriginal++;
      bcTuple currentBC(frame.bcAO2D, frame.bcEvSel);
      auto p = std::find(bcSet.begin(), bcSet.end(), currentBC);
      if (p == bcSet.end()) {
        bcSet.push_back(currentBC);
      } else {
        noriginalduplicate++;
      }
    }
    // For skimmed dataset
    bcSet.clear();
    auto skimmedFrames = getSelectedFrames(skimmedFile, trigger0Bit, trigger1Bit);
    for (auto& skimmedFrame : skimmedFrames) {
      nskimmed++;
      bcTuple currentBC(skimmedFrame.bcAO2D, skimmedFrame.bcEvSel);
      auto p = std::find(bcSet.begin(), bcSet.end(), currentBC);
      if (p == bcSet.end()) {
        bcSet.push_back(currentBC);
      } else {
        nskimmedduplicate++;
      }
    }
    if (!skipNoneDuplicate || noriginalduplicate != 0 || nskimmedduplicate != 0) {
      sel_labels.push_back(labels[i]);
      numOriginal.push_back(noriginal);
      numOriginalDuplicate.push_back(noriginalduplicate);
      numSkimmed.push_back(nskimmed);
      numSkimmedDuplicate.push_back(nskimmedduplicate);
    }
  }
  originalFile.Close();
  skimmedFile.Close();

  TH1D hOriginalTotal("hOriginalTotal", "AO2D Original;;Number of events", sel_labels.size(), 0, sel_labels.size());
  TH1D hOriginalDuplicate("hOriginalDuplicate", "Duplicate Trigger Original;;Number of events", sel_labels.size(), 0, sel_labels.size());
  TH1D hOriginalRatio("hOriginalRatio", (runNumber + " Original;;Duplicate / Total").data(), sel_labels.size(), 0, sel_labels.size());
  TH1D hSkimmedTotal("hSkimmedTotal", "AO2D Skimmed;;Number of events", sel_labels.size(), 0, sel_labels.size());
  TH1D hSkimmedDuplicate("hSkimmedDuplicate", "Duplicate Trigger Skimmed;;Number of events", sel_labels.size(), 0, sel_labels.size());
  TH1D hSkimmedRatio("hSkimmedRatio", (runNumber + " Skimmed;;Duplicate / Total").data(), sel_labels.size(), 0, sel_labels.size());

  for (int i = 0; i < sel_labels.size(); i++) {
    hOriginalTotal.GetXaxis()->SetBinLabel(i + 1, sel_labels[i].c_str());
    hOriginalDuplicate.GetXaxis()->SetBinLabel(i + 1, sel_labels[i].c_str());
    hOriginalRatio.GetXaxis()->SetBinLabel(i + 1, sel_labels[i].c_str());
    hOriginalTotal.SetBinContent(i + 1, numOriginal[i]);
    hOriginalDuplicate.SetBinContent(i + 1, numOriginalDuplicate[i]);
    if (hOriginalTotal.GetBinContent(i + 1) > 0) {
      hOriginalRatio.SetBinContent(i + 1, hOriginalDuplicate.GetBinContent(i + 1) / hOriginalTotal.GetBinContent(i + 1));
    } else {
      hOriginalRatio.SetBinContent(i + 1, 0);
    }

    hSkimmedTotal.GetXaxis()->SetBinLabel(i + 1, sel_labels[i].c_str());
    hSkimmedDuplicate.GetXaxis()->SetBinLabel(i + 1, sel_labels[i].c_str());
    hSkimmedRatio.GetXaxis()->SetBinLabel(i + 1, sel_labels[i].c_str());
    hSkimmedTotal.SetBinContent(i + 1, numSkimmed[i]);
    hSkimmedDuplicate.SetBinContent(i + 1, numSkimmedDuplicate[i]);
    if (hSkimmedTotal.GetBinContent(i + 1) > 0) {
      hSkimmedRatio.SetBinContent(i + 1, hSkimmedDuplicate.GetBinContent(i + 1) / hSkimmedTotal.GetBinContent(i + 1));
    } else {
      hSkimmedRatio.SetBinContent(i + 1, 0);
    }
  }

  TFile fout(outputFileName, "UPDATE");
  fout.cd();
  hOriginalTotal.Write();
  hOriginalDuplicate.Write();
  hOriginalRatio.Write();
  hSkimmedTotal.Write();
  hSkimmedDuplicate.Write();
  hSkimmedRatio.Write();
  fout.Close();
}

void checkBCrangesSkimming(std::string originalFileName = "bcRanges_fullrun.root", std::string skimmedFileName = "bcRanges_fullrun_skimmed.root")
{
  TH1F hTriggerCounter("hTriggerCounter", "hTriggerCounter", 3, 0.5, 3.5);
  hTriggerCounter.GetXaxis()->SetBinLabel(1, "Original");
  hTriggerCounter.GetXaxis()->SetBinLabel(2, "Skimmed");
  TH1F hNumCounter("hNumCounter", "hTriggerCounter", 10, -0.5, 9.5);
  TH1F hPairedTriggerCounter("hPairedTriggerCounter", "hPairedTriggerCounter", 4, 0.5, 4.5);
  hPairedTriggerCounter.GetXaxis()->SetBinLabel(1, "Total");
  hPairedTriggerCounter.GetXaxis()->SetBinLabel(2, "Same AO2D BC");
  hPairedTriggerCounter.GetXaxis()->SetBinLabel(3, "Same EvSel BC");
  hPairedTriggerCounter.GetXaxis()->SetBinLabel(4, "Same Both BC");
  TH1F hSinglePairCheck("hSinglePairCheck", "hSinglePairCheck", 4, 0.5, 4.5);
  hSinglePairCheck.GetXaxis()->SetBinLabel(1, "Total");
  hSinglePairCheck.GetXaxis()->SetBinLabel(2, "Same AO2D BC");
  hSinglePairCheck.GetXaxis()->SetBinLabel(3, "Same EvSel BC");
  hSinglePairCheck.GetXaxis()->SetBinLabel(4, "Same Both BC");
  TH1F hMultiPairCheck("hMultiPairCheck", "hMultiPairCheck", 5, 0.5, 5.5);
  hMultiPairCheck.GetXaxis()->SetBinLabel(1, "Total");
  hMultiPairCheck.GetXaxis()->SetBinLabel(2, "Total Pair");
  hMultiPairCheck.GetXaxis()->SetBinLabel(3, "Same AO2D BC");
  hMultiPairCheck.GetXaxis()->SetBinLabel(4, "Same EvSel BC");
  hMultiPairCheck.GetXaxis()->SetBinLabel(5, "Same Both BC");
  TH1F hBCDiffAO2D("hBCDiffAO2D", "hBCDiffAO2D", 1000, -5.e3, 5.e3);
  TH1F hBCDiffEvSel("hBCDiffEvSel", "hBCDiffEvSel", 1000, -5.e3, 5.e3);

  TH1F hBCOriginal("hBCOriginal", "hBCOriginal", 4, 0.5, 4.5);
  hBCOriginal.GetXaxis()->SetBinLabel(1, "Total");
  hBCOriginal.GetXaxis()->SetBinLabel(2, "Same AO2D BC");
  hBCOriginal.GetXaxis()->SetBinLabel(3, "Same EvSel BC");
  hBCOriginal.GetXaxis()->SetBinLabel(4, "Same Both BC");
  TH1F hBCSkimmed("hBCSkimmed", "hBCSkimmed", 4, 0.5, 4.5);
  hBCSkimmed.GetXaxis()->SetBinLabel(1, "Total");
  hBCSkimmed.GetXaxis()->SetBinLabel(2, "Same AO2D BC");
  hBCSkimmed.GetXaxis()->SetBinLabel(3, "Same EvSel BC");
  hBCSkimmed.GetXaxis()->SetBinLabel(4, "Same Both BC");

  auto t1 = std::chrono::steady_clock::now();
  TFile originalFile(originalFileName.c_str(), "READ");
  TFile skimmedFile(skimmedFileName.c_str(), "READ");
  auto originalFrames = getSelectedFrames(originalFile, trigger0Bit, trigger1Bit);
  auto skimmedFrames = getSelectedFrames(skimmedFile, trigger0Bit, trigger1Bit);
  originalFile.Close();
  skimmedFile.Close();
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
  for (auto frame : originalFrames) {
    if (frame.selMask[0] & trigger0Bit) {
      bool found = false;
      hTriggerCounter.Fill(1);
      hBCOriginal.Fill(1);
      auto p1 = std::find_if(bcSet.begin(), bcSet.end(), [&](const auto& val) { return val.bcAO2D == frame.bcAO2D; });
      if (p1 != bcSet.end()) {
        hBCOriginal.Fill(2);
      }
      auto p2 = std::find_if(bcSet.begin(), bcSet.end(), [&](const auto& val) { return val.bcEvSel == frame.bcEvSel; });
      if (p2 != bcSet.end()) {
        hBCOriginal.Fill(3);
      }
      bcTuple currentBC(frame.bcAO2D, frame.bcEvSel);
      auto p3 = std::find(bcSet.begin(), bcSet.end(), currentBC);
      if (p3 == bcSet.end()) {
        bcSet.push_back(currentBC);
      } else {
        hBCOriginal.Fill(4);
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
            hPairedTriggerCounter.Fill(1);
            if (frame.bcAO2D == skimmedFrame.bcAO2D) {
              hPairedTriggerCounter.Fill(2);
            }
            if (frame.bcEvSel == skimmedFrame.bcEvSel) {
              hPairedTriggerCounter.Fill(3);
              if (frame.bcAO2D == skimmedFrame.bcAO2D) {
                hPairedTriggerCounter.Fill(4);
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
        hSinglePairCheck.Fill(1);
        hBCDiffAO2D.Fill(frame.bcAO2D - skimmedbcs[0].bcAO2D);
        hBCDiffEvSel.Fill(frame.bcEvSel - skimmedbcs[0].bcEvSel);
        if (frame.bcAO2D == skimmedbcs[0].bcAO2D) {
          hSinglePairCheck.Fill(2);
        }
        if (frame.bcEvSel == skimmedbcs[0].bcEvSel) {
          hSinglePairCheck.Fill(3);
          if (frame.bcAO2D == skimmedbcs[0].bcAO2D) {
            hSinglePairCheck.Fill(4);
          }
        }
      } else {
        // std::cout << "Unexpected trigger!!!   " << n << std::endl;
        hMultiPairCheck.Fill(1);
        for (auto skimmedbc : skimmedbcs) {
          hMultiPairCheck.Fill(2);
          if (frame.bcAO2D == skimmedbc.bcAO2D) {
            hMultiPairCheck.Fill(3);
          }
          if (frame.bcEvSel == skimmedbc.bcEvSel) {
            hMultiPairCheck.Fill(4);
            if (frame.bcAO2D == skimmedbc.bcAO2D) {
              hMultiPairCheck.Fill(5);
            }
          }
        }
      }
      hNumCounter.Fill(n);
    }
  }
  auto t6 = std::chrono::steady_clock::now();
  int d3 = std::chrono::duration_cast<std::chrono::milliseconds>(t6 - t5).count();
  std::cout << "Search Time: " << d3 << std::endl;

  bcSet.clear();
  for (auto& skimmedFrame : skimmedFrames) {
    if (skimmedFrame.selMask[0] & trigger0Bit || skimmedFrame.selMask[1] & trigger1Bit) {
      hTriggerCounter.Fill(2);
      hBCSkimmed.Fill(1);
      auto p1 = std::find_if(bcSet.begin(), bcSet.end(), [&](const auto& val) { return val.bcAO2D == skimmedFrame.bcAO2D; });
      if (p1 != bcSet.end()) {
        hBCSkimmed.Fill(2);
      }
      auto p2 = std::find_if(bcSet.begin(), bcSet.end(), [&](const auto& val) { return val.bcEvSel == skimmedFrame.bcEvSel; });
      if (p2 != bcSet.end()) {
        hBCSkimmed.Fill(3);
      }
      bcTuple currentBC(skimmedFrame.bcAO2D, skimmedFrame.bcEvSel);
      auto p3 = std::find(bcSet.begin(), bcSet.end(), currentBC);
      if (p3 == bcSet.end()) {
        bcSet.push_back(currentBC);
      } else {
        hBCSkimmed.Fill(4);
      }
    }
  }

  TFile fout(outputFileName, "RECREATE");
  fout.cd();
  hTriggerCounter.Write();
  hBCOriginal.Write();
  hBCSkimmed.Write();
  hNumCounter.Write();
  hPairedTriggerCounter.Write();
  hSinglePairCheck.Write();
  hBCDiffAO2D.Write();
  hBCDiffEvSel.Write();
  hMultiPairCheck.Write();
  fout.Close();
}
