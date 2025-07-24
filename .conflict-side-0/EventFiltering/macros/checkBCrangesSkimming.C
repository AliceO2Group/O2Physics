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

#include <cmath>
#include <vector>
#include <regex>
#include <iostream>
#include <TFile.h>
#include <TGrid.h>
#include <TH1.h>
#include <TTree.h>
#include "CommonDataFormat/InteractionRecord.h"
#include "CommonDataFormat/IRFrame.h"

using o2::InteractionRecord;
using o2::dataformats::IRFrame;

// Set the bit of trigger which need to be checked
const ULong64_t bcDiffTolerance = 0;
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
  selectedFrames(ULong64_t bcAO2D, ULong64_t bcEvSel, const IRFrame& frame) : IRFrame(frame), bcAO2D(bcAO2D), bcEvSel(bcEvSel), triMask{0, 0}, selMask{0, 0} {}
  selectedFrames(ULong64_t bcAO2D, ULong64_t bcEvSel, ULong64_t triMask[2], ULong64_t selMask[2], const IRFrame& frame) : IRFrame(frame), bcAO2D(bcAO2D), bcEvSel(bcEvSel), triMask{triMask[0], triMask[1]}, selMask{selMask[0], selMask[1]} {}
  ULong64_t triMask[2]{0ull}, selMask[2]{0ull}, bcAO2D, bcEvSel;
  int numSameTriggerInNearbyBCs = 0; // related to bcDiffTolerance
  bool isSingle() { return numSameTriggerInNearbyBCs == 0; }
  void SetNInNearbyBC(int n) { numSameTriggerInNearbyBCs = n; }
  int GetNInNearbyBC() { return numSameTriggerInNearbyBCs; }
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

int DoBCSubraction(selectedFrames bc1, selectedFrames bc2)
{
  if (bc1.getMin() > bc2.getMax()) {
    return DoBCSubraction(bc1.getMin().toLong(), bc2.getMax().toLong());
  } else if (bc1.getMax() < bc2.getMin()) {
    return DoBCSubraction(bc1.getMax().toLong(), bc2.getMin().toLong());
  } else {
    return 0;
  }
}

bool isClose(selectedFrames a, selectedFrames b, ULong64_t bcDiffTolerance)
{
  if (a.getMin() > b.getMax() + bcDiffTolerance || a.getMax() + bcDiffTolerance < b.getMin())
    return false;
  else
    return true;
}

std::vector<std::vector<selectedFrames>> getFrames(std::unique_ptr<TFile>& file, int trgIDStart, int N)
{
  ULong64_t bcAO2D{0ull}, bcEvSel{0ull}, triMask[2]{0ull}, selMask[2]{0ull};
  std::vector<std::vector<selectedFrames>> frames;
  frames.resize(N);
  for (auto key : *file->GetListOfKeys()) {
    auto dir = dynamic_cast<TDirectory*>(file->Get(key->GetName()));
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
      for (int trgID = trgIDStart; trgID < trgIDStart + N; trgID++) {
        ULong64_t trigger0Bit = 0, trigger1Bit = 0;
        if (trgID < 64) {
          trigger0Bit = BIT(trgID);
        } else {
          trigger1Bit = BIT(trgID - 64);
        }
        if (selMask[0] & trigger0Bit || selMask[1] & trigger1Bit) {
          InteractionRecord irstart, irend;
          irstart.setFromLong(std::min(bcAO2D, bcEvSel));
          irend.setFromLong(std::max(bcAO2D, bcEvSel));
          IRFrame frame(irstart, irend);
          int index = trgID - trgIDStart;
          frames[index].push_back({bcAO2D, bcEvSel, triMask, selMask, frame});
        }
      }
    }
  }

  return frames;
}

std::vector<selectedFrames> getSelectedFrames(std::unique_ptr<TFile>& file, int trgID)
{
  auto frames = getFrames(file, trgID, 1);
  return frames[0];
}

// Check how many other triggers are in a compatible BC window with the current one
// Ideally, most of triggers are singles (num = 1)
// which means for most triggered events, none of others is in a nearby time window
void checkNearbyBCs(std::vector<selectedFrames>& frames, ULong64_t bcDiffTolerance)
{
  std::sort(frames.begin(), frames.end(), [](const selectedFrames& a, const selectedFrames& b) {
    if (a.getMin() != b.getMin()) {
      return a.getMin() < b.getMin();
    } else {
      return a.getMax() < b.getMax();
    }
  });
  int firstTrg = 0;
  for (auto& currentFrame : frames) {
    int num = 0;
    bool shouldUpdate = true; // true if the maxBC of event in loop is smaller than the evaluating one ->  update firstTrg
    for (int i = firstTrg; i < frames.size(); i++) {
      auto& frame = frames[i];
      if (frame.getMin() > currentFrame.getMax() + bcDiffTolerance) {
        break;
      }
      if (isClose(currentFrame, frame, bcDiffTolerance)) {
        shouldUpdate = false;
        bool found = currentFrame.selMask[0] & frame.selMask[0] || currentFrame.selMask[1] & frame.selMask[1];
        if (found) {
          num++;
        }
      } else {
        if (shouldUpdate) {
          firstTrg = i;
        }
      }
    }
    currentFrame.SetNInNearbyBC(num);
  }
}

// Get RunNumber
std::string getRunNumber(std::string fileName)
{
  std::string runNumber = "";
  std::regex re("/5[0-9]*");
  std::smatch match;
  if (std::regex_search(fileName, match, re)) {
    // Remove the leading '/'
    runNumber = match.str().substr(1);
  }
  return runNumber;
}

// Detailed checks for specific trigger, not enabled by default
void checkBCForSelectedTrg(std::vector<selectedFrames>& originalFrames, std::vector<selectedFrames>& skimmedFrames, string runNumber, string triggerLabel)
{

  TH1D hTriggerCounter("hTriggerCounter", (runNumber + " " + triggerLabel + ";;Total number of trigger").data(), 2, -0.5, 1.5);
  hTriggerCounter.GetXaxis()->SetBinLabel(1, "Original");
  hTriggerCounter.GetXaxis()->SetBinLabel(2, "Skimmed");
  TH1D hBCDiffAO2D("hBCDiffAO2D", (runNumber + " " + triggerLabel + ";;#DeltaBC_{AO2D} between paired singles").data(), 201, -100.5, 100.5);
  TH1D hBCDiffEvSel("hBCDiffEvSel", (runNumber + " " + triggerLabel + ";;#DeltaBC_{EvSel} between paired singles").data(), 201, -100.5, 100.5);

  TH1D hBCOriginal("hBCOriginal", (runNumber + " " + triggerLabel + " Original;;Trigger counts").data(), 4, -0.5, 3.5);
  hBCOriginal.GetXaxis()->SetBinLabel(1, "Total");
  hBCOriginal.GetXaxis()->SetBinLabel(2, "Same AO2D BC");
  hBCOriginal.GetXaxis()->SetBinLabel(3, "Same EvSel BC");
  hBCOriginal.GetXaxis()->SetBinLabel(4, "Same Both BC");
  TH1D hBCSkimmed("hBCSkimmed", (runNumber + " " + triggerLabel + " Skimmed;;Trigger counts").data(), 4, -0.5, 3.5);
  hBCSkimmed.GetXaxis()->SetBinLabel(1, "Total");
  hBCSkimmed.GetXaxis()->SetBinLabel(2, "Same AO2D BC");
  hBCSkimmed.GetXaxis()->SetBinLabel(3, "Same EvSel BC");
  hBCSkimmed.GetXaxis()->SetBinLabel(4, "Same Both BC");

  TH1D hMatchedNumCounter("hMatchedNumCounter", (runNumber + " " + triggerLabel + ";;Number of matched triggers in skimmed data").data(), 10, -0.5, 9.5);

  checkNearbyBCs(originalFrames, bcDiffTolerance);
  checkNearbyBCs(skimmedFrames, bcDiffTolerance);

  std::vector<bcTuple> bcSet;
  int firstTrg = 0;
  for (int i = 0; i < originalFrames.size(); i++) {
    auto& frame = originalFrames[i];
    hTriggerCounter.Fill(0);
    hBCOriginal.Fill(0);
    //------------------------------ Check if there are triggers which have same BC, time-consuming! -------------------------------------------------------
    auto p1 = std::find_if(bcSet.begin(), bcSet.end(), [&](const auto& val) { return val.bcAO2D == frame.bcAO2D; });
    if (p1 != bcSet.end()) {
      hBCOriginal.Fill(1);
    }
    auto p2 = std::find_if(bcSet.begin(), bcSet.end(), [&](const auto& val) { return val.bcEvSel == frame.bcEvSel; });
    if (p2 != bcSet.end()) {
      hBCOriginal.Fill(2);
    }
    bcTuple currentBC(frame.bcAO2D, frame.bcEvSel);
    auto p3 = std::find(bcSet.begin(), bcSet.end(), currentBC);
    if (p3 == bcSet.end()) {
      bcSet.push_back(currentBC);
    } else {
      hBCOriginal.Fill(3);
    }
    //-------------------------------------------------------------------------------------

    if (frame.GetNInNearbyBC() != 1) {
      continue; // Only check singles
    }
    std::vector<bcTuple> skimmedbcs;
    int n = 0;
    bool shouldUpdate = true;
    for (int j = firstTrg; j < skimmedFrames.size(); j++) {
      auto& skimmedFrame = skimmedFrames[j];
      if (skimmedFrame.getMin() > frame.getMax()) {
        break;
      }
      if (skimmedFrame.GetNInNearbyBC() != 1) {
        continue; // Only check singles
      }
      if (isClose(frame, skimmedFrame, bcDiffTolerance)) {
        shouldUpdate = false;
        bool found = frame.selMask[0] & skimmedFrame.selMask[0] || frame.selMask[1] & skimmedFrame.selMask[1];
        if (found) {
          // Additional check to avoid match of skimmed singles and original multiplies
          if (i != 0 && isClose(originalFrames[i - 1], skimmedFrame, bcDiffTolerance)) {
            continue;
          }
          if (i != originalFrames.size() && isClose(originalFrames[i + 1], skimmedFrame, bcDiffTolerance)) {
            continue;
          }
          skimmedbcs.push_back({skimmedFrame.bcAO2D, skimmedFrame.bcEvSel});
          n++;
        }
      } else {
        if (shouldUpdate) {
          firstTrg = j;
        }
      }
    }
    if (n == 0) {
      // std::cout << "Trigger not found!!!" << std::endl;
    } else if (n == 1) {
      hBCDiffAO2D.Fill(DoBCSubraction(frame.bcAO2D, skimmedbcs[0].bcAO2D));
      hBCDiffEvSel.Fill(DoBCSubraction(frame.bcEvSel, skimmedbcs[0].bcEvSel));
    }
    hMatchedNumCounter.Fill(n);
  }

  //------------------------------ Check if there are triggers which have same BC, time-consuming! -------------------------------------------------------
  bcSet.clear();
  for (auto& skimmedFrame : skimmedFrames) {
    hTriggerCounter.Fill(1);
    hBCSkimmed.Fill(0);
    auto p1 = std::find_if(bcSet.begin(), bcSet.end(), [&](const auto& val) { return val.bcAO2D == skimmedFrame.bcAO2D; });
    if (p1 != bcSet.end()) {
      hBCSkimmed.Fill(1);
    }
    auto p2 = std::find_if(bcSet.begin(), bcSet.end(), [&](const auto& val) { return val.bcEvSel == skimmedFrame.bcEvSel; });
    if (p2 != bcSet.end()) {
      hBCSkimmed.Fill(2);
    }
    bcTuple currentBC(skimmedFrame.bcAO2D, skimmedFrame.bcEvSel);
    auto p3 = std::find(bcSet.begin(), bcSet.end(), currentBC);
    if (p3 == bcSet.end()) {
      bcSet.push_back(currentBC);
    } else {
      hBCSkimmed.Fill(3);
    }
  }
  //-------------------------------------------------------------------------------------

  TFile fout(outputFileName, "UPDATE");
  fout.cd();
  TDirectory* dir1 = fout.GetDirectory(runNumber.data());
  if (!dir1) {
    dir1 = fout.mkdir(runNumber.data());
  }
  dir1->cd();
  TDirectory* dir2 = dir1->GetDirectory(triggerLabel.data());
  if (!dir2) {
    dir2 = dir1->mkdir(triggerLabel.data());
  }
  dir2->cd();

  hTriggerCounter.Write();
  hBCOriginal.Write();
  hBCSkimmed.Write();
  hBCDiffAO2D.Write();
  hBCDiffEvSel.Write();
  hMatchedNumCounter.Write();
  fout.Close();
}

// Detailed checks for specific trigger
void checkBCForSelectedTrg(std::string AnaFileName = "AnalysisResults.root", std::string originalFileName = "bcRanges_fullrun.root", std::string skimmedFileName = "bcRanges_fullrun_skimmed.root", int triggerID = 1, bool useAlien = true)
{

  string runNumber = getRunNumber(originalFileName);
  if (useAlien) {
    TGrid::Connect("alien://");
    AnaFileName = "alien://" + AnaFileName;
    originalFileName = "alien://" + originalFileName;
    skimmedFileName = "alien://" + skimmedFileName;
  }

  // Readin labels
  std::unique_ptr<TFile> AnaFile{TFile::Open(AnaFileName.c_str(), "READ")};
  TH1* hist0 = dynamic_cast<TH1*>(AnaFile->Get("central-event-filter-task/scalers/mFiltered;1"));
  std::vector<std::string> labels;
  std::vector<int> binNum;
  for (int i = 1; i <= hist0->GetNbinsX(); i++) {
    std::string label = hist0->GetXaxis()->GetBinLabel(i);
    if (label != "Total number of events" && label != "Filtered events") {
      labels.push_back(label);
      binNum.push_back(i);
    }
  }
  AnaFile->Close();
  std::string triggerLabel = labels[triggerID];

  std::unique_ptr<TFile> originalFile{TFile::Open(originalFileName.c_str(), "READ")};
  std::unique_ptr<TFile> skimmedFile{TFile::Open(skimmedFileName.c_str(), "READ")};
  auto originalFrames = getSelectedFrames(originalFile, triggerID);
  auto skimmedFrames = getSelectedFrames(skimmedFile, triggerID);
  originalFile->Close();
  skimmedFile->Close();

  checkBCForSelectedTrg(originalFrames, skimmedFrames, runNumber, triggerLabel);
}

// Check the BCId compatibility of triggers on original and skimmmed data
void checkBCrangesSkimming(std::string AnaFileName = "AnalysisResults.root", std::string originalFileName = "bcRanges_fullrun.root", std::string skimmedFileName = "bcRanges_fullrun_skimmed.root", bool useAlien = true)
{

  string runNumber = getRunNumber(originalFileName);
  if (useAlien) {
    TGrid::Connect("alien://");
    AnaFileName = "alien://" + AnaFileName;
    originalFileName = "alien://" + originalFileName;
    skimmedFileName = "alien://" + skimmedFileName;
  }

  // Readin labels
  std::unique_ptr<TFile> AnaFile{TFile::Open(AnaFileName.c_str(), "READ")};
  TH1* hist0 = dynamic_cast<TH1*>(AnaFile->Get("central-event-filter-task/scalers/mFiltered;1"));
  std::vector<std::string> labels;
  std::vector<int> binNum;
  for (int i = 1; i <= hist0->GetNbinsX(); i++) {
    std::string label = hist0->GetXaxis()->GetBinLabel(i);
    if (label != "Total number of events" && label != "Filtered events") {
      labels.push_back(label);
      binNum.push_back(i);
      // std::cout << i - 2 << ": " << label << std::endl;
    }
  }
  AnaFile->Close();

  // Due to potential selection on triggers, histograms should be created later
  // for example: skip triggers which have no enrties
  std::vector<std::string> sel_labels;
  std::vector<int> numOriginal, numSkimmed, numOriginalSingle, numSkimmedSingle, numOriginalDouble, numSkimmedDouble, numOriginalMultiple, numSkimmedMultiple, numCloseSkimmed;
  std::vector<int> numpair, numpairedBCAO2D, numpairedBCEvSel;
  std::vector<double> avgDeltaBCAO2D, avgDeltaBCEvSel, avgDeltaBC, rmsDeltaBCAO2D, rmsDeltaBCEvSel, rmsDeltaBC;
  std::vector<double> avgNumPairedTrigger, rmsNumPairedTrigger;

  std::unique_ptr<TFile> originalFile{TFile::Open(originalFileName.c_str(), "READ")};
  std::unique_ptr<TFile> skimmedFile{TFile::Open(skimmedFileName.c_str(), "READ")};
  std::vector<std::vector<selectedFrames>> originalAllFrames = getFrames(originalFile, 0, labels.size());
  std::vector<std::vector<selectedFrames>> skimmedAllFrames = getFrames(skimmedFile, 0, labels.size());
  for (int trgID = 0; trgID < labels.size(); trgID++) {
    // Caculate singles, doubles, and multiples
    int noriginal{0}, nskimmed{0}, noriginalsingle{0}, nskimmedsingle{0}, noriginaldouble{0}, nskimmeddouble{0}, noriginalmultiple{0}, nskimmedmultiple{0};
    // Caculate mean and rms of diff BC
    TH1D hDiffBCAO2DCount("hDiffBCAO2DCount", "hDiffBCAO2DCount", 21, -10.5, 10.5);
    TH1D hDiffBCEvSelCount("hDiffBCEvSelCount", "hDiffBCEvSelCount", 21, -10.5, 10.5);
    TH1D hDiffBCCount("hDiffBCCount", "hDiffBCCount", 21, -10.5, 10.5);
    TH1D hNumPairedTriggerCount("hNumPairedTriggerCount", "hNumPairedTriggerCount", 10, -0.5, 9.5);
    // For Original dataset
    auto& originalFrames = originalAllFrames[trgID];
    checkNearbyBCs(originalFrames, bcDiffTolerance); // include sorting
    noriginal = originalFrames.size();
    for (auto originalFrame : originalFrames) {
      if (originalFrame.GetNInNearbyBC() == 0) {
        std::cerr << "Unexpected trigger!!! " << std::endl;
      } else if (originalFrame.GetNInNearbyBC() == 1) {
        noriginalsingle++;
      } else if (originalFrame.GetNInNearbyBC() == 2) {
        noriginaldouble++;
      } else {
        noriginalmultiple++;
      }
    }
    // For skimmed dataset
    auto& skimmedFrames = skimmedAllFrames[trgID];
    checkNearbyBCs(skimmedFrames, bcDiffTolerance); // include sorting
    nskimmed = skimmedFrames.size();
    for (auto& skimmedFrame : skimmedFrames) {
      if (skimmedFrame.GetNInNearbyBC() == 0) {
        std::cerr << "Unexpected trigger!!! " << std::endl;
      } else if (skimmedFrame.GetNInNearbyBC() == 1) {
        nskimmedsingle++;
      } else if (skimmedFrame.GetNInNearbyBC() == 2) {
        nskimmeddouble++;
      } else {
        nskimmedmultiple++;
      }
    }

    // Check BC differences
    int npair{0}, npairedBCAO2D{0}, npairedBCEvSel{0}, ncloseskimmed{0}, maxdeltaBCAO2D{0}, maxdeltaBCEvSel{0};
    int firstTrg = 0;
    for (int i = 0; i < originalFrames.size(); i++) {
      auto& frame = originalFrames[i];
      if (frame.GetNInNearbyBC() != 1) {
        continue; // Only check singles
      }
      std::vector<selectedFrames> skimmedbcs;
      int n = 0;
      bool shouldUpdate = true;
      for (int j = firstTrg; j < skimmedFrames.size(); j++) {
        auto& skimmedFrame = skimmedFrames[j];
        if (skimmedFrame.getMin() > frame.getMax()) {
          break;
        }
        if (skimmedFrame.GetNInNearbyBC() != 1) {
          continue; // Only check singles
        }
        if (isClose(frame, skimmedFrame, bcDiffTolerance)) {
          shouldUpdate = false;
          bool found = frame.selMask[0] & skimmedFrame.selMask[0] || frame.selMask[1] & skimmedFrame.selMask[1];
          if (found) {
            // Additional check to avoid match of skimmed singles and original multiplies
            if (i != 0 && isClose(originalFrames[i - 1], skimmedFrame, bcDiffTolerance)) {
              continue;
            }
            if (i != originalFrames.size() && isClose(originalFrames[i + 1], skimmedFrame, bcDiffTolerance)) {
              continue;
            }

            InteractionRecord irstart, irend;
            irstart.setFromLong(std::min(skimmedFrame.bcAO2D, skimmedFrame.bcEvSel));
            irend.setFromLong(std::max(skimmedFrame.bcAO2D, skimmedFrame.bcEvSel));
            IRFrame frame(irstart, irend);
            skimmedbcs.push_back({skimmedFrame.bcAO2D, skimmedFrame.bcEvSel, frame});
            n++;
          }
        } else {
          if (shouldUpdate) {
            firstTrg = j;
          }
        }
      }
      if (n == 1) {
        npair++;
        int bcdiffAO2D = DoBCSubraction(frame.bcAO2D, skimmedbcs[0].bcAO2D);
        int bcdiffEvSel = DoBCSubraction(frame.bcEvSel, skimmedbcs[0].bcEvSel);
        hDiffBCAO2DCount.Fill(std::abs(bcdiffAO2D));
        hDiffBCEvSelCount.Fill(std::abs(bcdiffEvSel));
        hDiffBCCount.Fill(std::abs(DoBCSubraction(frame, skimmedbcs[0])));
        if (frame.bcAO2D == skimmedbcs[0].bcAO2D) {
          npairedBCAO2D++;
        }
        if (frame.bcEvSel == skimmedbcs[0].bcEvSel) {
          npairedBCEvSel++;
        }
      }
      ncloseskimmed += n;
      hNumPairedTriggerCount.Fill(n);
    }

    // if (static_cast<double>(ncloseskimmed) / noriginal > 0.95 || noriginal == 0)
    if (noriginal == 0) {
      // continue;
    }
    sel_labels.push_back(labels[trgID]);
    numOriginal.push_back(noriginal);
    numOriginalSingle.push_back(noriginalsingle);
    numOriginalDouble.push_back(noriginaldouble);
    numOriginalMultiple.push_back(noriginalmultiple);
    numSkimmed.push_back(nskimmed);
    numSkimmedSingle.push_back(nskimmedsingle);
    numSkimmedDouble.push_back(nskimmeddouble);
    numSkimmedMultiple.push_back(nskimmedmultiple);

    numpair.push_back(npair);
    numpairedBCAO2D.push_back(npairedBCAO2D);
    numpairedBCEvSel.push_back(npairedBCEvSel);
    numCloseSkimmed.push_back(ncloseskimmed);
    avgDeltaBCAO2D.push_back(hDiffBCAO2DCount.GetMean());
    avgDeltaBCEvSel.push_back(hDiffBCEvSelCount.GetMean());
    avgDeltaBC.push_back(hDiffBCCount.GetMean());
    rmsDeltaBCAO2D.push_back(hDiffBCAO2DCount.GetRMS());
    rmsDeltaBCEvSel.push_back(hDiffBCEvSelCount.GetRMS());
    rmsDeltaBC.push_back(hDiffBCCount.GetRMS());
    avgNumPairedTrigger.push_back(hNumPairedTriggerCount.GetMean());
    rmsNumPairedTrigger.push_back(hNumPairedTriggerCount.GetRMS());
  }
  originalFile->Close();
  skimmedFile->Close();

  TH1D hOriginalTotal("hOriginalTotal", (runNumber + " Original;;Number of events").data(), sel_labels.size(), 0, sel_labels.size());
  TH1D hOriginalSingles("hOriginalSingles", (runNumber + " Original;;Number of singles").data(), sel_labels.size(), 0, sel_labels.size());
  TH1D hOriginalDoubles("hOriginalDoubles", (runNumber + " Original;;Number of doubles").data(), sel_labels.size(), 0, sel_labels.size());
  TH1D hOriginalMultiples("hOriginalMultiples", (runNumber + " Original;;Number of multiples").data(), sel_labels.size(), 0, sel_labels.size());

  TH1D hSkimmedTotal("hSkimmedTotal", (runNumber + " Skimmed;;Number of events").data(), sel_labels.size(), 0, sel_labels.size());
  TH1D hSkimmedSingles("hSkimmedSingles", (runNumber + " Skimmed;;Number of singles").data(), sel_labels.size(), 0, sel_labels.size());
  TH1D hSkimmedDoubles("hSkimmedDoubles", (runNumber + " Skimmed;;Number of doubles").data(), sel_labels.size(), 0, sel_labels.size());
  TH1D hSkimmedMultiples("hSkimmedMultiples", (runNumber + " Skimmed;;Number of multiples").data(), sel_labels.size(), 0, sel_labels.size());

  TH1D hTriggerMatchesRatio("hTriggerMatchesRatio", (runNumber + " Skimmed Efficiency;; Matched skimmed triggers / Original singles").data(), sel_labels.size(), 0, sel_labels.size());       // the ratio of triggers in skimmed dataset whose BC is compatible with original triggers to the number of original triggers, might be duplicate since we check it based on every trigger in unskimmed data
  TH1D hTriggerSingleMatchesRatio("hTriggerSingleMatchesRatio", (runNumber + " Skimmed Efficiency;; One-to-one matches / Original singles").data(), sel_labels.size(), 0, sel_labels.size()); // the ratio of 1-1 paired triggers to the number of original triggers
  TH1D hMatchesSameBCAO2DRatio("hMatchesSameBCAO2DRatio", (runNumber + " One-to-one matches;; Matchess with same BC_{AO2D} / Total").data(), sel_labels.size(), 0, sel_labels.size());        // In 1-1 matches, the ratio of matches who have same BCAO2D
  TH1D hMatchesSameBCEvSelRatio("hMatchesSameBCEvSelRatio", (runNumber + " One-to-one matches;; Matches with same BC_{EvSel} / Total").data(), sel_labels.size(), 0, sel_labels.size());      // In 1-1 matches, the ratio of matches who have same BCEvSel
  TH1D hDiffBCAO2D("hDiffBCAO2D", (runNumber + " One-to-one matches;;|#DeltaBC_{AO2D}|").data(), sel_labels.size(), 0, sel_labels.size());                                                    // difference in BCAO2D of 1-1 matches
  TH1D hDiffBCEvSel("hDiffBCEvSel", (runNumber + " One-to-one matches;;|#DeltaBC_{EvSel}|").data(), sel_labels.size(), 0, sel_labels.size());                                                 // difference in BCEvSel of 1-1 matches
  TH1D hDiffBC("hDiffBC", (runNumber + " One-to-one matches;;|#DeltaBC|").data(), sel_labels.size(), 0, sel_labels.size());                                                                   // difference between the BC tuple, expected to be 0 if bcDiffTolerance = 0
  TH1D hNumMatchesInSkimmed("hNumMatchesInSkimmed", (runNumber + " number of matched triggers in skimmed data;;Matched trigger count").data(), sel_labels.size(), 0, sel_labels.size());      // number of triggers in skimmed data which are compatible in the BC ranges of singles in original selection

  for (int i = 0; i < sel_labels.size(); i++) {
    // Original data
    hOriginalTotal.GetXaxis()->SetBinLabel(i + 1, sel_labels[i].c_str());
    hOriginalSingles.GetXaxis()->SetBinLabel(i + 1, sel_labels[i].c_str());
    hOriginalDoubles.GetXaxis()->SetBinLabel(i + 1, sel_labels[i].c_str());
    hOriginalMultiples.GetXaxis()->SetBinLabel(i + 1, sel_labels[i].c_str());

    hOriginalTotal.SetBinContent(i + 1, numOriginal[i]);
    hOriginalTotal.SetBinError(i + 1, std::sqrt(numOriginal[i]));
    hOriginalSingles.SetBinContent(i + 1, numOriginalSingle[i]);
    hOriginalSingles.SetBinError(i + 1, std::sqrt(numOriginalSingle[i]));
    hOriginalDoubles.SetBinContent(i + 1, numOriginalDouble[i]);
    hOriginalDoubles.SetBinError(i + 1, std::sqrt(numOriginalDouble[i]));
    hOriginalMultiples.SetBinContent(i + 1, numOriginalMultiple[i]);
    hOriginalMultiples.SetBinError(i + 1, std::sqrt(numOriginalMultiple[i]));

    // Skimmed data
    hSkimmedTotal.GetXaxis()->SetBinLabel(i + 1, sel_labels[i].c_str());
    hSkimmedSingles.GetXaxis()->SetBinLabel(i + 1, sel_labels[i].c_str());
    hSkimmedDoubles.GetXaxis()->SetBinLabel(i + 1, sel_labels[i].c_str());
    hSkimmedMultiples.GetXaxis()->SetBinLabel(i + 1, sel_labels[i].c_str());

    hSkimmedTotal.SetBinContent(i + 1, numSkimmed[i]);
    hSkimmedTotal.SetBinError(i + 1, std::sqrt(numSkimmed[i]));
    hSkimmedSingles.SetBinContent(i + 1, numSkimmedSingle[i]);
    hSkimmedSingles.SetBinError(i + 1, std::sqrt(numSkimmedSingle[i]));
    hSkimmedDoubles.SetBinContent(i + 1, numSkimmedDouble[i]);
    hSkimmedDoubles.SetBinError(i + 1, std::sqrt(numSkimmedDouble[i]));
    hSkimmedMultiples.SetBinContent(i + 1, numSkimmedMultiple[i]);
    hSkimmedMultiples.SetBinError(i + 1, std::sqrt(numSkimmedMultiple[i]));

    // Matches QA
    hTriggerMatchesRatio.GetXaxis()->SetBinLabel(i + 1, sel_labels[i].c_str());
    hTriggerSingleMatchesRatio.GetXaxis()->SetBinLabel(i + 1, sel_labels[i].c_str());
    hMatchesSameBCAO2DRatio.GetXaxis()->SetBinLabel(i + 1, sel_labels[i].c_str());
    hMatchesSameBCEvSelRatio.GetXaxis()->SetBinLabel(i + 1, sel_labels[i].c_str());
    hDiffBCAO2D.GetXaxis()->SetBinLabel(i + 1, sel_labels[i].c_str());
    hDiffBCEvSel.GetXaxis()->SetBinLabel(i + 1, sel_labels[i].c_str());
    hDiffBC.GetXaxis()->SetBinLabel(i + 1, sel_labels[i].c_str());
    hNumMatchesInSkimmed.GetXaxis()->SetBinLabel(i + 1, sel_labels[i].c_str());

    if (numpair[i] > 0) {
      hMatchesSameBCAO2DRatio.SetBinContent(i + 1, static_cast<double>(numpairedBCAO2D[i]) / numpair[i]);
      hMatchesSameBCEvSelRatio.SetBinContent(i + 1, static_cast<double>(numpairedBCEvSel[i]) / numpair[i]);
    }
    hTriggerMatchesRatio.SetBinContent(i + 1, numCloseSkimmed[i]);
    hTriggerMatchesRatio.SetBinError(i + 1, std::sqrt(numCloseSkimmed[i]));
    hTriggerSingleMatchesRatio.SetBinContent(i + 1, numpair[i]);
    hTriggerSingleMatchesRatio.SetBinError(i + 1, std::sqrt(numpair[i]));
    hDiffBCAO2D.SetBinContent(i + 1, avgDeltaBCAO2D[i]);
    hDiffBCAO2D.SetBinError(i + 1, rmsDeltaBCAO2D[i]);
    hDiffBCEvSel.SetBinContent(i + 1, avgDeltaBCEvSel[i]);
    hDiffBCEvSel.SetBinError(i + 1, rmsDeltaBCEvSel[i]);
    hDiffBC.SetBinContent(i + 1, avgDeltaBC[i]);
    hDiffBC.SetBinError(i + 1, rmsDeltaBC[i]);
    hNumMatchesInSkimmed.SetBinContent(i + 1, avgNumPairedTrigger[i]);
    hNumMatchesInSkimmed.SetBinError(i + 1, rmsNumPairedTrigger[i]);
  }

  TH1D* hTriggerEff; // Ratio of the total number of triggers in skimmed data to that in original data (not the real efficiency since the downscalings are removed in skimmed for this QA)
  TH1D *hOriginalSinglesRatio, *hOriginalDoublesRatio, *hOriginalMultiplesRatio;
  TH1D *hSkimmedSinglesRatio, *hSkimmedDoublesRatio, *hSkimmedMultiplesRatio;

  hTriggerEff = reinterpret_cast<TH1D*>(hSkimmedTotal.Clone("hTriggerEff"));
  hTriggerEff->SetTitle((runNumber + " skimmed efficiency;; Skimmed / Original").data());
  hTriggerEff->Divide(&hOriginalTotal);
  hTriggerMatchesRatio.Divide(&hOriginalSingles);
  hTriggerSingleMatchesRatio.Divide(&hOriginalSingles);
  hOriginalSinglesRatio = reinterpret_cast<TH1D*>(hOriginalSingles.Clone("hOriginalSinglesRatio"));
  hOriginalSinglesRatio->SetTitle((runNumber + " Original;;Singles / Total").data());
  hOriginalSinglesRatio->Divide(&hOriginalTotal);
  hOriginalDoublesRatio = reinterpret_cast<TH1D*>(hOriginalDoubles.Clone("hOriginalDoublesRatio"));
  hOriginalDoublesRatio->SetTitle((runNumber + " Original;;Doubles / Total").data());
  hOriginalDoublesRatio->Divide(&hOriginalTotal);
  hOriginalMultiplesRatio = reinterpret_cast<TH1D*>(hOriginalMultiples.Clone("hOriginalMultiplesRatio"));
  hOriginalMultiplesRatio->SetTitle((runNumber + " Original;;Multiples / Total").data());
  hOriginalMultiplesRatio->Divide(&hOriginalTotal);

  hSkimmedSinglesRatio = reinterpret_cast<TH1D*>(hSkimmedSingles.Clone("hSkimmedSinglesRatio"));
  hSkimmedSinglesRatio->SetTitle((runNumber + " Skimmed;;Singles / Total").data());
  hSkimmedSinglesRatio->Divide(&hSkimmedTotal);
  hSkimmedDoublesRatio = reinterpret_cast<TH1D*>(hSkimmedDoubles.Clone("hSkimmedDoublesRatio"));
  hSkimmedDoublesRatio->SetTitle((runNumber + " Skimmed;;Doubles / Total").data());
  hSkimmedDoublesRatio->Divide(&hSkimmedTotal);
  hSkimmedMultiplesRatio = reinterpret_cast<TH1D*>(hSkimmedMultiples.Clone("hSkimmedMultiplesRatio"));
  hSkimmedMultiplesRatio->SetTitle((runNumber + " Skimmed;;Multiples / Total").data());
  hSkimmedMultiplesRatio->Divide(&hSkimmedTotal);

  TFile fout(outputFileName, "UPDATE");
  fout.cd();
  TDirectory* dir = fout.mkdir(runNumber.data());
  dir->cd();
  hTriggerEff->Write();
  hTriggerMatchesRatio.Write();
  hTriggerSingleMatchesRatio.Write();
  hDiffBCAO2D.Write();
  hDiffBCEvSel.Write();
  hNumMatchesInSkimmed.Write();
  if (bcDiffTolerance > 0) {
    hDiffBC.Write();
  }
  TDirectory* dirextra = dir->mkdir("ExtraQA");
  dirextra->cd();
  hOriginalTotal.Write();
  hOriginalSingles.Write();
  hOriginalDoubles.Write();
  hOriginalMultiples.Write();
  hOriginalSinglesRatio->Write();
  hOriginalDoublesRatio->Write();
  hOriginalMultiplesRatio->Write();
  hSkimmedTotal.Write();
  hSkimmedSingles.Write();
  hSkimmedDoubles.Write();
  hSkimmedMultiples.Write();
  hSkimmedSinglesRatio->Write();
  hSkimmedDoublesRatio->Write();
  hSkimmedMultiplesRatio->Write();
  hMatchesSameBCAO2DRatio.Write();
  hMatchesSameBCEvSelRatio.Write();
  fout.Close();

  // Do checks for trigger
  for (int trgID = 0; trgID < labels.size(); trgID++) {
    // if (trgID == 77 || trgID == 78 || trgID == 79) {
    // checkBCForSelectedTrg(originalAllFrames[trgID], skimmedAllFrames[trgID], runNumber, labels[trgID]);
    //}
  }
}
