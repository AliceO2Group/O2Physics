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

#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TStyle.h>

#include <fstream>
#include <iostream>
#include <regex>
#include <string>
#include <vector>

void checkSkimming(std::string original_path = "AnalysisResults.root", std::string skimmed_path = "AnalysisResults_skimmed.root", TFile* outputFile = nullptr, bool skipDownscaled = true)
{
  gStyle->SetOptStat(0);
  std::string runNumber = "";
  std::regex re("/5[0-9]*");
  std::smatch match;
  if (std::regex_search(original_path, match, re)) {
    // Remove the leading '/'
    runNumber = match.str().substr(1);
  }

  // Load the root files
  TFile file1(original_path.c_str());
  TFile file2(skimmed_path.c_str());

  // Extract the histograms
  TH1* hist0 = dynamic_cast<TH1*>(file1.Get("central-event-filter-task/scalers/mScalers;1"));
  TH1* hist1 = dynamic_cast<TH1*>(file1.Get("central-event-filter-task/scalers/mFiltered;1"));
  TH1* hist2 = dynamic_cast<TH1*>(file2.Get("central-event-filter-task/scalers/mScalers;1"));

  if (!hist0 || !hist1 || !hist2) {
    std::cerr << "Error: Failed to extract histograms from the root files." << std::endl;
    return;
  }

  std::vector<std::string> labels;
  for (int i = 1; i <= hist1->GetNbinsX(); i++) {
    std::string label = hist1->GetXaxis()->GetBinLabel(i);
    if (label != "Total number of events" && label != "Filtered events" && (!skipDownscaled || hist0->GetBinContent(i) == hist1->GetBinContent(i)) && hist0->GetBinContent(i) > 0) {
      labels.push_back(label);
    }
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

  bool localFile = !outputFile;
  if (!outputFile) {
    outputFile = new TFile("output.root", "RECREATE");
  }
  if (!runNumber.empty()) {
    outputFile->mkdir(runNumber.data());
  }
  outputFile->cd(runNumber.data());

  TH1D hOriginal("hOriginal", "Original AO2D;;Number of events", labels.size(), 0, labels.size());       // Histogram for the original values
  TH1D hSkimmed("hSkimmed", "AO2D from skimmed CTF;;Number of events", labels.size(), 0, labels.size()); // Histogram for the skimmed values
  TH1D hRatio("hRatio", (runNumber + ";;Skimmed / Original").data(), labels.size(), 0, labels.size());   // Histogram for the ratio of the two

  // Fill the histograms
  for (int i = 0; i < labels.size(); i++) {
    hOriginal.SetBinContent(i + 1, selected_bins1[i]);
    hSkimmed.SetBinContent(i + 1, selected_bins2[i]);
    hOriginal.GetXaxis()->SetBinLabel(i + 1, labels[i].c_str());
    hSkimmed.GetXaxis()->SetBinLabel(i + 1, labels[i].c_str());
    hRatio.GetXaxis()->SetBinLabel(i + 1, labels[i].c_str());
    hRatio.SetBinContent(i + 1, selected_bins2[i] / selected_bins1[i]);
  }
  hOriginal.Write();
  hSkimmed.Write();
  hRatio.Write();

  if (localFile) {
    outputFile->Close();
  }
}

void checkSkimming(std::string listName = "period.txt", bool skipDownscaled = true)
{
  std::string periodName = listName.substr(0, listName.find_last_of('.'));
  std::ifstream file(listName);
  std::string line;
  TFile* outputFile = new TFile((periodName + ".root").data(), "RECREATE");
  TCanvas c1("c1", "c1", 800, 600);
  c1.SetGridy();
  int lineCounter = 0;
  while (std::getline(file, line)) {
    lineCounter++;
  }
  file.clear();
  file.seekg(0, std::ios::beg);
  int counter = 0;
  while (std::getline(file, line)) {
    size_t pos = line.find(",");
    std::string original_path = line.substr(0, pos);
    std::string skimmed_path = line.substr(pos + 1);
    checkSkimming(original_path, skimmed_path, outputFile, skipDownscaled);
    TH1* hRatio = (TH1*)gDirectory->Get("hRatio");
    hRatio->Draw();
    std::string suffix = counter == 0 ? "(" : (counter == lineCounter - 1 ? ")" : "");
    c1.Print((periodName + ".pdf" + suffix).data());
    counter++;
  }
  outputFile->Close();
}
 