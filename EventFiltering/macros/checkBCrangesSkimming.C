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
#include "CommonDataFormat/InteractionRecord.h"
#include "CommonDataFormat/IRFrame.h"

using o2::InteractionRecord;
using o2::dataformats::IRFrame;

struct selectedFrames : public IRFrame {
  selectedFrames(ULong64_t triMask[2], ULong64_t selMask[2], const IRFrame& frame) : IRFrame(frame), triMask{triMask[0], triMask[1]}, selMask{selMask[0], selMask[1]} {}
  ULong64_t triMask[2]{0ull}, selMask[2]{0ull};
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
      InteractionRecord irstart, irend;
      irstart.setFromLong(std::min(bcAO2D, bcEvSel));
      irend.setFromLong(std::max(bcAO2D, bcEvSel));
      IRFrame frame(irstart, irend);
      selectedFrames.push_back({triMask, selMask, frame});
    }
  }
  return selectedFrames;
}

void checkBCrangesSkimming(std::string originalFileName = "bcRanges_fullrun.root", std::string skimmedFileName = "bcRanges_fullrun_skimmed.root")
{
  TFile originalFile(originalFileName.c_str(), "READ");
  TFile skimmedFile(skimmedFileName.c_str(), "READ");

  auto originalFrames = getSelectedFrames(originalFile);
  auto skimmedFrames = getSelectedFrames(skimmedFile);

  std::sort(originalFrames.begin(), originalFrames.end(), [](const selectedFrames& a, const selectedFrames& b) {
    return a.getMin() < b.getMin();
  });
  std::sort(skimmedFrames.begin(), skimmedFrames.end(), [](const selectedFrames& a, const selectedFrames& b) {
    return a.getMin() < b.getMin();
  });

  for (auto frame : originalFrames) {
    bool found = false;
    frame.getMin() -= 1000;
    frame.getMax() += 1000;
    for (auto& skimmedFrame : skimmedFrames) {
      if (skimmedFrame.getMin() > frame.getMax()) {
        break;
      }
      if (!frame.isOutside(skimmedFrame)) {
        found = frame.selMask[0] & skimmedFrame.selMask[0] || frame.selMask[1] & skimmedFrame.selMask[1];
      }
    }
  }
}