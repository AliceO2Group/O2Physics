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
// O2 includes// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
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

#include <TCanvas.h>
#include <TH1.h>
#include <TFile.h>

#include <iostream>
#include <vector>
#include <string>

void checkSkimming(std::string original_path = "AnalysisResults.root", std::string skimmed_path = "AnalysisResults_skimmed.root")
{
  // Define the labels to be compared (only triggers not downscaled)

  // Load the root files
  TFile file1(original_path.c_str());
  TFile file2(skimmed_path.c_str());

  // Extract the histograms
  TH1* hist1 = dynamic_cast<TH1*>(file1.Get("central-event-filter-task/scalers/mFiltered;1")); // Replace with the correct path
  TH1* hist2 = dynamic_cast<TH1*>(file2.Get("central-event-filter-task/scalers/mScalers;1"));  // Replace with the correct path

  std::vector<std::string> labels;
  for (int i = 1; i <= hist1->GetNbinsX(); i++) {
    std::string label = hist1->GetXaxis()->GetBinLabel(i);
    if (label != "Total number of events" && label != "Filtered events") {
      labels.push_back(label);
    }
  }

  if (!hist1 || !hist2) {
    std::cerr << "Error: Failed to extract histograms from the root files." << std::endl;
    return;
  }

  // Find the bins corresponding to the desired labels
  std::vector<double> selected_bins1, selected_bins2;
  for (auto lab : labels) {
    int bin1 = hist1->GetXaxis()->FindBin(lab.c_str());
    if (bin1 == -1) {
      std::cerr << "Error: Label " << lab << " not found in histogram 1." << std::endl;
      return;
    }
    selected_bins1.push_back(hist1->GetBinContent(bin1));

    int bin2 = hist2->GetXaxis()->FindBin(lab.c_str());
    if (bin2 == -1) {
      std::cerr << "Error: Label " << lab << " not found in histogram 2." << std::endl;
      return;
    }
    selected_bins2.push_back(hist2->GetBinContent(bin2));
  }

  TFile output("output.root", "RECREATE");
  TH1D* hOriginal = new TH1D("hOriginal", "Original AO2D;;Number of events", labels.size(), 0, labels.size());       // Histogram for the original values
  TH1D* hSkimmed = new TH1D("hSkimmed", "AO2D from skimmed CTF;;Number of events", labels.size(), 0, labels.size()); // Histogram for the skimmed values
  TH1D* hRatio = new TH1D("hRatio", ";;Skimmed / Original", labels.size(), 0, labels.size());                        // Histogram for the ratio of the two

  // Fill the histograms
  for (int i = 0; i < labels.size(); i++) {
    hOriginal->SetBinContent(i + 1, selected_bins1[i]);
    hSkimmed->SetBinContent(i + 1, selected_bins2[i]);
    hOriginal->GetXaxis()->SetBinLabel(i + 1, labels[i].c_str());
    hSkimmed->GetXaxis()->SetBinLabel(i + 1, labels[i].c_str());

    if (selected_bins1[i] < 1. || labels[i].find("Single") != std::string::npos || labels[i].find("Low") != std::string::npos || labels[i].find("Pt3P") != std::string::npos) {
      hRatio->GetXaxis()->SetBinLabel(i + 1, "Disabled");
      continue;
    }
    hRatio->GetXaxis()->SetBinLabel(i + 1, labels[i].c_str());
    hRatio->SetBinContent(i + 1, selected_bins2[i] / selected_bins1[i]);
  }
  hOriginal->Write();
  hSkimmed->Write();
  hRatio->Write();
}
