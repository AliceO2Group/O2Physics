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
const ULong64_t Trigger0BIT = BIT(61);
const ULong64_t Trigger1BIT = 0;
const ULong64_t bcDiffTolerance = 100;
const char outputFileName[15] = "output.root";

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
  int numSameTriggerInNearbyBCs = 0; // related to bcDiffTolerance
  bool isSingle() { return numSameTriggerInNearbyBCs == 0; }
  void SetNum(int n) { numSameTriggerInNearbyBCs = n; }
  int GetNum() { return numSameTriggerInNearbyBCs; }
};

int DoBCSubraction(ULong64_t bc1, ULong64_t bc2)
{
  if (bc1 > bc2) {
    return bc1 - bc2;
  } else {
    ULong64_t bcsub = bc2 - bc1;
    return -static_cast<int>(bcsub);
  }
}

bool isClose(selectedFrames a, selectedFrames b, ULong64_t bcDiffTolerance)
{
  if (a.getMin() > b.getMax() + bcDiffTolerance || a.getMax() < b.getMin() - bcDiffTolerance)
    return false;
  else
    return true;
}

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

void checkNearbyBCs(std::vector<selectedFrames>& frames, ULong64_t bcDiffTolerance)
{
  std::sort(frames.begin(), frames.end(), [](const selectedFrames& a, const selectedFrames& b) {
    return a.getMin() < b.getMin();
  });
  int firstID = 0;
  for (auto& currentFrame : frames) {
    int num = 0;
    bool isFirst = true;
    for (int i = firstID; i < frames.size(); i++) {
      auto& frame = frames[i];
      if (frame.getMin() > currentFrame.getMax() + bcDiffTolerance) {
        break;
      }
      if (isClose(currentFrame, frame, bcDiffTolerance)) {
        isFirst = false;
        bool found = currentFrame.selMask[0] & frame.selMask[0] || currentFrame.selMask[1] & frame.selMask[1];
        if (found) {
          num++;
        }
      } else {
        if (isFirst) {
          firstID = i;
        }
      }
    }
    currentFrame.SetNum(num);
  }
}

// Calulate the ratio of duplicate triggers
void checkDuplicateTriggerAndBCs(std::string AnaFileName = "AnalysisResults.root", std::string originalFileName = "bcRanges_fullrun.root", std::string skimmedFileName = "bcRanges_fullrun_skimmed.root")
{

  // Get RunNumber
  std::string runNumber = "";
  std::regex re("/5[0-9]*");
  std::smatch match;
  if (std::regex_search(originalFileName, match, re)) {
    // Remove the leading '/'
    runNumber = match.str().substr(1);
  }

  // Checks for BC difference between original and skimming data, and the ratio of triggers which have BCdiff==0
  TH1D hPairedNumCounterTotal("hPairedNumCounterTotal", "hPairedNumCounterTotal", 10, -0.5, 9.5);
  TH1D hBCDiffAO2DTotal("hBCDiffAO2DTotal", "hBCDiffAO2DTotal", 201, -100.5, 100.5);
  TH1D hBCDiffEvSelTotal("hBCDiffEvSelTotal", "hBCDiffEvSelTotal", 201, -100.5, 100.5);

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
  std::vector<int> numOriginal, numSkimmed, numOriginalSingle, numSkimmedSingle, numOriginalDouble, numSkimmedDouble, numOriginalMultiple, numSkimmedMultiple;
  std::vector<int> numpair, numpairedBCAO2D, numpairedBCEvSel, maxDeltaBCAO2D, maxDeltaBCEvSel;
  for (int i = 0; i < labels.size(); i++) {
    // std::cout << "i:" << i << std::endl;
    ULong64_t trigger0Bit = 0, trigger1Bit = 0;
    int triggerBit = binNum[i] - 2;
    if (triggerBit < 64) {
      trigger0Bit = BIT(triggerBit);
    } else {
      trigger1Bit = BIT(triggerBit - 64);
    }
    // Caculate singles, doubles, and multiples
    // For Original dataset
    std::vector<bcTuple> bcSet;
    std::vector<ULong64_t> bcFullSet;
    int noriginal{0}, nskimmed{0}, noriginalsingle{0}, nskimmedsingle{0}, noriginaldouble{0}, nskimmeddouble{0}, noriginalmultiple{0}, nskimmedmultiple{0}, maxdiffBCAO2D{0}, maxdiffBCEvSel{0};
    auto originalFrames = getSelectedFrames(originalFile, trigger0Bit, trigger1Bit);
    checkNearbyBCs(originalFrames, bcDiffTolerance);
    noriginal = originalFrames.size();
    for (auto originalFrame : originalFrames) {
      if (originalFrame.GetNum() == 0) {
        std::cerr << "Unexpected trigger!!! " << std::endl;
      } else if (originalFrame.GetNum() == 1) {
        noriginalsingle++;
      } else if (originalFrame.GetNum() == 2) {
        noriginaldouble++;
      } else {
        noriginalmultiple++;
      }
    }
    // For skimmed dataset
    auto skimmedFrames = getSelectedFrames(skimmedFile, trigger0Bit, trigger1Bit);
    checkNearbyBCs(skimmedFrames, bcDiffTolerance);
    nskimmed = skimmedFrames.size();
    for (auto& skimmedFrame : skimmedFrames) {
      if (skimmedFrame.GetNum() == 0) {
        std::cerr << "Unexpected trigger!!! " << std::endl;
      } else if (skimmedFrame.GetNum() == 1) {
        nskimmedsingle++;
      } else if (skimmedFrame.GetNum() == 2) {
        nskimmeddouble++;
      } else {
        nskimmedmultiple++;
      }
    }

    sel_labels.push_back(labels[i]);
    numOriginal.push_back(noriginal);
    numOriginalSingle.push_back(noriginalsingle);
    numOriginalDouble.push_back(noriginaldouble);
    numOriginalMultiple.push_back(noriginalmultiple);
    numSkimmed.push_back(nskimmed);
    numSkimmedSingle.push_back(nskimmedsingle);
    numSkimmedDouble.push_back(nskimmeddouble);
    numSkimmedMultiple.push_back(nskimmedmultiple);

    // Check BC differences
    int npair{0}, npairedBCAO2D{0}, npairedBCEvSel{0}, maxdeltaBCAO2D{0}, maxdeltaBCEvSel{0};
    int firstID = 0;
    for (auto frame : originalFrames) {
      if (frame.selMask[0] & trigger0Bit || frame.selMask[1] & trigger1Bit) {
        // std::cout << "------------------------------------------------" << std::endl;
        if (frame.GetNum() != 1) {
          continue; // Only check singles
        }
        std::vector<bcTuple> skimmedbcs;
        int n = 0;
        bool isFirst = true;
        for (int i = firstID; i < skimmedFrames.size(); i++) {
          auto& skimmedFrame = skimmedFrames[i];
          if (skimmedFrame.getMin() > frame.getMax()) {
            break;
          }
          if (skimmedFrame.GetNum() != 1) {
            continue; // Only check singles
          }
          if (isClose(frame, skimmedFrame, bcDiffTolerance)) {
            isFirst = false;
            bool found = frame.selMask[0] & skimmedFrame.selMask[0] || frame.selMask[1] & skimmedFrame.selMask[1];
            // found = found && (frame.bcAO2D == skimmedFrame.bcAO2D || frame.bcEvSel == skimmedFrame.bcEvSel);
            if (found) {
              skimmedbcs.push_back({skimmedFrame.bcAO2D, skimmedFrame.bcEvSel});
              n++;
            }
          } else {
            if (isFirst) {
              firstID = i;
            }
          }
        }
        if (n == 1) {
          npair++;
          int bcdiffAO2D = DoBCSubraction(frame.bcAO2D, skimmedbcs[0].bcAO2D);
          int bcdiffEvSel = DoBCSubraction(frame.bcEvSel, skimmedbcs[0].bcEvSel);
          hBCDiffAO2DTotal.Fill(bcdiffAO2D);
          hBCDiffEvSelTotal.Fill(bcdiffEvSel);
          maxdiffBCAO2D = std::max(std::abs(maxdiffBCAO2D), bcdiffAO2D);
          maxdiffBCEvSel = std::max(std::abs(maxdiffBCEvSel), bcdiffEvSel);
          if (frame.bcAO2D == skimmedbcs[0].bcAO2D) {
            npairedBCAO2D++;
          }
          if (frame.bcEvSel == skimmedbcs[0].bcEvSel) {
            npairedBCEvSel++;
          }
        }
        hPairedNumCounterTotal.Fill(n);
      }
    }
    numpair.push_back(npair);
    numpairedBCAO2D.push_back(npairedBCAO2D);
    numpairedBCEvSel.push_back(npairedBCEvSel);
    maxDeltaBCAO2D.push_back(maxdiffBCAO2D);
    maxDeltaBCEvSel.push_back(maxdiffBCEvSel);
  }
  originalFile.Close();
  skimmedFile.Close();

  TH1D hOriginalTotal("hOriginalTotal", (runNumber + " AO2D Original;;Number of events").data(), sel_labels.size(), 0, sel_labels.size());
  TH1D hOriginalSingles("hOriginalSingles", (runNumber + " Original;;Number of Singles").data(), sel_labels.size(), 0, sel_labels.size());
  TH1D hOriginalSinglesRatio("hOriginalSinglesRatio", (runNumber + " Original;;Singles / Total").data(), sel_labels.size(), 0, sel_labels.size());
  TH1D hOriginalDoubles("hOriginalDoubles", (runNumber + " Original;;Number of Doubles").data(), sel_labels.size(), 0, sel_labels.size());
  TH1D hOriginalDoublesRatio("hOriginalDoublesRatio", (runNumber + " Original;;Doubles / Total").data(), sel_labels.size(), 0, sel_labels.size());
  TH1D hOriginalMultiples("hOriginalMultiples", (runNumber + " Original;;Number of Multiples").data(), sel_labels.size(), 0, sel_labels.size());
  TH1D hOriginalMultiplesRatio("hOriginalMultiplesRatio", (runNumber + " Original;;Multiples / Total").data(), sel_labels.size(), 0, sel_labels.size());

  TH1D hSkimmedTotal("hSkimmedTotal", (runNumber + " AO2D Skimmed;;Number of events").data(), sel_labels.size(), 0, sel_labels.size());
  TH1D hSkimmedSingles("hSkimmedSingles", (runNumber + " Skimmed;;Number of Singles").data(), sel_labels.size(), 0, sel_labels.size());
  TH1D hSkimmedSinglesRatio("hSkimmedSinglesRatio", (runNumber + " Skimmed;;Singles / Total").data(), sel_labels.size(), 0, sel_labels.size());
  TH1D hSkimmedDoubles("hSkimmedDoubles", (runNumber + " Skimmed;;Number of Doubles").data(), sel_labels.size(), 0, sel_labels.size());
  TH1D hSkimmedDoublesRatio("hSkimmedDoublesRatio", (runNumber + " Skimmed;;Doubles / Total").data(), sel_labels.size(), 0, sel_labels.size());
  TH1D hSkimmedMultiples("hSkimmedMultiples", (runNumber + " Skimmed;;Number of Multiples").data(), sel_labels.size(), 0, sel_labels.size());
  TH1D hSkimmedMultiplesRatio("hSkimmedMultiplesRatio", (runNumber + " Skimmed;;Multiples / Total").data(), sel_labels.size(), 0, sel_labels.size());

  TH1D hPairedBCAO2DRatio("hPairedBCAO2DRatio", (runNumber + " One-to-One Pairs;; Pairs with same BCAO2D / Total").data(), sel_labels.size(), 0, sel_labels.size());
  TH1D hPairedBCEvSelRatio("hPairedBCEvSelRatio", (runNumber + " One-to-One Pairs;; Pairs with same BCEvSel / Total").data(), sel_labels.size(), 0, sel_labels.size());
  TH1D hMaxDiffBCAO2D("hMaxDiffBCAO2D", (runNumber + " One-to-One Pairs;;|#DeltaBCAO2D|_{max}").data(), sel_labels.size(), 0, sel_labels.size());
  TH1D hMaxDiffBCEvSel("hMaxDiffBCEvSel", (runNumber + " One-to-One Pairs;;|#DeltaBCEvSel|_{max}").data(), sel_labels.size(), 0, sel_labels.size());

  for (int i = 0; i < sel_labels.size(); i++) {
    hOriginalTotal.GetXaxis()->SetBinLabel(i + 1, sel_labels[i].c_str());
    hOriginalSingles.GetXaxis()->SetBinLabel(i + 1, sel_labels[i].c_str());
    hOriginalMultiples.GetXaxis()->SetBinLabel(i + 1, sel_labels[i].c_str());
    hOriginalSinglesRatio.GetXaxis()->SetBinLabel(i + 1, sel_labels[i].c_str());
    hOriginalDoublesRatio.GetXaxis()->SetBinLabel(i + 1, sel_labels[i].c_str());
    hOriginalMultiplesRatio.GetXaxis()->SetBinLabel(i + 1, sel_labels[i].c_str());

    hOriginalTotal.SetBinContent(i + 1, numOriginal[i]);
    hOriginalSingles.SetBinContent(i + 1, numOriginalSingle[i]);
    hOriginalDoubles.SetBinContent(i + 1, numOriginalDouble[i]);
    hOriginalMultiples.SetBinContent(i + 1, numOriginalMultiple[i]);
    if (numOriginal[i] > 0) {
      hOriginalSinglesRatio.SetBinContent(i + 1, static_cast<double>(numOriginalSingle[i]) / numOriginal[i]);
      hOriginalDoublesRatio.SetBinContent(i + 1, static_cast<double>(numOriginalSingle[i]) / numOriginal[i]);
      hOriginalMultiplesRatio.SetBinContent(i + 1, static_cast<double>(numOriginalMultiple[i]) / numOriginal[i]);
    } else {
      hOriginalSinglesRatio.SetBinContent(i + 1, 0);
      hOriginalDoublesRatio.SetBinContent(i + 1, 0);
      hOriginalMultiplesRatio.SetBinContent(i + 1, 0);
    }

    hSkimmedTotal.GetXaxis()->SetBinLabel(i + 1, sel_labels[i].c_str());
    hSkimmedSingles.GetXaxis()->SetBinLabel(i + 1, sel_labels[i].c_str());
    hSkimmedMultiples.GetXaxis()->SetBinLabel(i + 1, sel_labels[i].c_str());
    hSkimmedSinglesRatio.GetXaxis()->SetBinLabel(i + 1, sel_labels[i].c_str());
    hSkimmedDoublesRatio.GetXaxis()->SetBinLabel(i + 1, sel_labels[i].c_str());
    hSkimmedMultiplesRatio.GetXaxis()->SetBinLabel(i + 1, sel_labels[i].c_str());

    hSkimmedTotal.SetBinContent(i + 1, numSkimmed[i]);
    hSkimmedSingles.SetBinContent(i + 1, numSkimmedSingle[i]);
    hSkimmedDoubles.SetBinContent(i + 1, numSkimmedDouble[i]);
    hSkimmedMultiples.SetBinContent(i + 1, numSkimmedMultiple[i]);
    if (numSkimmed[i] > 0) {
      hSkimmedSinglesRatio.SetBinContent(i + 1, static_cast<double>(numSkimmedSingle[i]) / numSkimmed[i]);
      hSkimmedDoublesRatio.SetBinContent(i + 1, static_cast<double>(numSkimmedDouble[i]) / numSkimmed[i]);
      hSkimmedMultiplesRatio.SetBinContent(i + 1, static_cast<double>(numSkimmedMultiple[i]) / numSkimmed[i]);
    } else {
      hSkimmedSinglesRatio.SetBinContent(i + 1, 0);
      hSkimmedDoublesRatio.SetBinContent(i + 1, 0);
      hSkimmedMultiplesRatio.SetBinContent(i + 1, 0);
    }

    hPairedBCAO2DRatio.GetXaxis()->SetBinLabel(i + 1, sel_labels[i].c_str());
    hPairedBCEvSelRatio.GetXaxis()->SetBinLabel(i + 1, sel_labels[i].c_str());
    if (numpair[i] > 0) {
      hPairedBCAO2DRatio.SetBinContent(i + 1, static_cast<double>(numpairedBCAO2D[i]) / numpair[i]);
      hPairedBCEvSelRatio.SetBinContent(i + 1, static_cast<double>(numpairedBCEvSel[i]) / numpair[i]);
    }

    hMaxDiffBCAO2D.GetXaxis()->SetBinLabel(i + 1, sel_labels[i].c_str());
    hMaxDiffBCEvSel.GetXaxis()->SetBinLabel(i + 1, sel_labels[i].c_str());
    hMaxDiffBCAO2D.SetBinContent(i + 1, maxDeltaBCAO2D[i]);
    hMaxDiffBCEvSel.SetBinContent(i + 1, maxDeltaBCEvSel[i]);
  }

  TFile fout(outputFileName, "UPDATE");
  fout.cd();
  hOriginalTotal.Write();
  hOriginalSingles.Write();
  hOriginalDoubles.Write();
  hOriginalMultiples.Write();
  hOriginalSinglesRatio.Write();
  hOriginalDoublesRatio.Write();
  hOriginalMultiplesRatio.Write();
  hSkimmedTotal.Write();
  hSkimmedSingles.Write();
  hSkimmedDoubles.Write();
  hSkimmedMultiples.Write();
  hSkimmedSinglesRatio.Write();
  hSkimmedDoublesRatio.Write();
  hSkimmedMultiplesRatio.Write();
  hPairedNumCounterTotal.Write();
  hBCDiffAO2DTotal.Write();
  hBCDiffEvSelTotal.Write();
  hPairedBCAO2DRatio.Write();
  hPairedBCEvSelRatio.Write();
  hMaxDiffBCAO2D.Write();
  hMaxDiffBCEvSel.Write();
  fout.Close();
}

// Detailed checks for specific trigger
void checkBCrangesSkimming(std::string originalFileName = "bcRanges_fullrun.root", std::string skimmedFileName = "bcRanges_fullrun_skimmed.root")
{
  //---------------------------------For specific trigger----------------------------------
  TH1D hTriggerCounter("hTriggerCounter", "hTriggerCounter", 3, 0.5, 3.5);
  hTriggerCounter.GetXaxis()->SetBinLabel(1, "Original");
  hTriggerCounter.GetXaxis()->SetBinLabel(2, "Skimmed");
  TH1D hNumCounter("hNumCounter", "hNumCounter", 10, -0.5, 9.5);
  TH1D hSinglePairCheck("hSinglePairCheck", "hSinglePairCheck", 4, 0.5, 4.5);
  hSinglePairCheck.GetXaxis()->SetBinLabel(1, "Total");
  hSinglePairCheck.GetXaxis()->SetBinLabel(2, "Same AO2D BC");
  hSinglePairCheck.GetXaxis()->SetBinLabel(3, "Same EvSel BC");
  hSinglePairCheck.GetXaxis()->SetBinLabel(4, "Same Both BC");
  TH1D hMultiPairCheck("hMultiPairCheck", "hMultiPairCheck", 5, 0.5, 5.5);
  hMultiPairCheck.GetXaxis()->SetBinLabel(1, "Total");
  hMultiPairCheck.GetXaxis()->SetBinLabel(2, "Total Pair");
  hMultiPairCheck.GetXaxis()->SetBinLabel(3, "Same AO2D BC");
  hMultiPairCheck.GetXaxis()->SetBinLabel(4, "Same EvSel BC");
  hMultiPairCheck.GetXaxis()->SetBinLabel(5, "Same Both BC");
  TH1D hBCDiffAO2D("hBCDiffAO2D", "hBCDiffAO2D", 2001, -1000.5, 1000.5);
  TH1D hBCDiffEvSel("hBCDiffEvSel", "hBCDiffEvSel", 2001, -1000.5, 1000.5);

  TH1D hBCOriginal("hBCOriginal", "hBCOriginal", 4, 0.5, 4.5);
  hBCOriginal.GetXaxis()->SetBinLabel(1, "Total");
  hBCOriginal.GetXaxis()->SetBinLabel(2, "Same AO2D BC");
  hBCOriginal.GetXaxis()->SetBinLabel(3, "Same EvSel BC");
  hBCOriginal.GetXaxis()->SetBinLabel(4, "Same Both BC");
  TH1D hBCSkimmed("hBCSkimmed", "hBCSkimmed", 4, 0.5, 4.5);
  hBCSkimmed.GetXaxis()->SetBinLabel(1, "Total");
  hBCSkimmed.GetXaxis()->SetBinLabel(2, "Same AO2D BC");
  hBCSkimmed.GetXaxis()->SetBinLabel(3, "Same EvSel BC");
  hBCSkimmed.GetXaxis()->SetBinLabel(4, "Same Both BC");

  auto t1 = std::chrono::steady_clock::now();
  TFile originalFile(originalFileName.c_str(), "READ");
  TFile skimmedFile(skimmedFileName.c_str(), "READ");
  auto originalFrames = getSelectedFrames(originalFile, Trigger0BIT, Trigger1BIT);
  auto skimmedFrames = getSelectedFrames(skimmedFile, Trigger0BIT, Trigger1BIT);
  originalFile.Close();
  skimmedFile.Close();
  auto t2 = std::chrono::steady_clock::now();
  int d1 = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
  std::cout << "Readin Time: " << d1 << std::endl;

  auto t3 = std::chrono::steady_clock::now();
  checkNearbyBCs(originalFrames, bcDiffTolerance);
  checkNearbyBCs(skimmedFrames, bcDiffTolerance);
  auto t4 = std::chrono::steady_clock::now();
  int d2 = std::chrono::duration_cast<std::chrono::milliseconds>(t4 - t3).count();
  std::cout << "Sort Time: " << d2 << std::endl;

  auto t5 = std::chrono::steady_clock::now();
  std::vector<bcTuple> bcSet;
  for (auto frame : originalFrames) {
    if (frame.selMask[0] & Trigger0BIT || frame.selMask[1] & Trigger1BIT) {
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
      if (frame.GetNum() != 1) {
        continue; // Only check singles
      }
      std::vector<bcTuple> skimmedbcs;
      int n = 0;
      for (auto& skimmedFrame : skimmedFrames) {
        if (skimmedFrame.getMin() > frame.getMax()) {
          break;
        }
        if (skimmedFrame.GetNum() != 1) {
          continue; // Only check singles
        }
        if (isClose(frame, skimmedFrame, bcDiffTolerance)) {
          bool found = frame.selMask[0] & skimmedFrame.selMask[0] || frame.selMask[1] & skimmedFrame.selMask[1];
          // found = found && (frame.bcAO2D == skimmedFrame.bcAO2D || frame.bcEvSel == skimmedFrame.bcEvSel);
          if (found) {
            skimmedbcs.push_back({skimmedFrame.bcAO2D, skimmedFrame.bcEvSel});
            n++;
          }
        }
      }
      if (n == 0) {
        // std::cout << "Trigger not found!!!" << std::endl;
      } else if (n == 1) {
        hSinglePairCheck.Fill(1);
        hBCDiffAO2D.Fill(DoBCSubraction(frame.bcAO2D, skimmedbcs[0].bcAO2D));
        hBCDiffEvSel.Fill(DoBCSubraction(frame.bcEvSel, skimmedbcs[0].bcEvSel));
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
        // std::cout << "Unexpected trigger!!! n=" << n << std::endl;
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
    if (skimmedFrame.selMask[0] & Trigger0BIT || skimmedFrame.selMask[1] & Trigger1BIT) {
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
  hSinglePairCheck.Write();
  hBCDiffAO2D.Write();
  hBCDiffEvSel.Write();
  hMultiPairCheck.Write();
  fout.Close();
}
