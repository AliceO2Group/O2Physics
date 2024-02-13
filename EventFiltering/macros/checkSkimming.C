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

void checkSkimming(std::string original_path = "alice/data/2022/LHC22m/523308/apass4/0140/OfflineTriggerSelection/Stage_5/001/AnalysisResults_fullrun.root", std::string skimmed_path = "AnalysisResults.root")
{
  // Define the labels to be compared (only triggers not downscaled)
  std::vector<std::string> labels = {"fPHOSnbar", "fPHOSPair", "fPHOSElectron", "fPHOSPhoton", "fOmegaLargeRadius", "fSingleXiYN", "fQuadrupleXi", "fhadronXi", "fTripleXi", "fGammaHighPtDCAL", "fGammaHighPtEMCAL", "fJetFullHighPt", "fLD", "fPD", "fLLL", "fPLL", "fPPL", "fPPP", "fHighFt0cFv0Flat", "fHighFt0cFv0Mult", "fHighFt0Flat", "fHighFt0Mult", "fHighTrackMult", "fHfDoubleCharmMix", "fHfDoubleCharm3P", "fHfSoftGamma3P", "fHfFemto2P", "fHfBeauty4P", "fHfFemto3P", "fHfBeauty3P", "fHfSoftGamma2P", "fHfDoubleCharm2P", "fDiMuon", "fDiElectron", "fUDdiff", "fHe"};

  // Load the root files
  TFile file1(original_path.c_str());
  TFile file2(skimmed_path.c_str());

  // Extract the histograms
  TH1* hist1 = dynamic_cast<TH1*>(file1.Get("central-event-filter-task/scalers/mFiltered;1")); // Replace with the correct path
  TH1* hist2 = dynamic_cast<TH1*>(file2.Get("central-event-filter-task/scalers/mScalers;1"));  // Replace with the correct path

  if (!hist1 || !hist2) {
    std::cerr << "Error: Failed to extract histograms from the root files." << std::endl;
    return;
  }

  // Find the bins corresponding to the desired labels
  std::vector<int> selected_bins1, selected_bins2;
  for (auto lab : labels) {
    int bin1 = hist1->GetXaxis()->FindBin(lab.c_str());
    if (bin1 == -1) {
      std::cerr << "Error: Label " << lab << " not found in histogram 1." << std::endl;
      return;
    }
    selected_bins1.push_back(bin1);

    int bin2 = hist2->GetXaxis()->FindBin(lab.c_str());
    if (bin2 == -1) {
      std::cerr << "Error: Label " << lab << " not found in histogram 2." << std::endl;
      return;
    }
    selected_bins2.push_back(bin2);
  }

  TFile output("output.root", "RECREATE");
  TH1D* hOriginal = new TH1D("hOriginal", "Original AO2D;;Number of events", labels.size(), 0, labels.size());       // Histogram for the original values
  TH1D* hSkimmed = new TH1D("hSkimmed", "AO2D from skimmed CTF;;Number of events", labels.size(), 0, labels.size()); // Histogram for the skimmed values
  TH1D* hRatio = new TH1D("hRatio", ";;Skimmed / Original", labels.size(), 0, labels.size());                        // Histogram for the ratio of the two

  // Fill the histograms
  for (int i = 0; i < labels.size(); i++) {
    hOriginal->SetBinContent(i + 1, hist1->GetBinContent(selected_bins1[i]));
    hSkimmed->SetBinContent(i + 1, hist2->GetBinContent(selected_bins2[i]));
    hRatio->SetBinContent(i + 1, hist2->GetBinContent(selected_bins2[i]) / hist1->GetBinContent(selected_bins1[i]));
    hOriginal->GetXaxis()->SetBinLabel(i + 1, labels[i].c_str());
    hSkimmed->GetXaxis()->SetBinLabel(i + 1, labels[i].c_str());
    hRatio->GetXaxis()->SetBinLabel(i + 1, labels[i].c_str());
  }
  hOriginal->Write();
  hSkimmed->Write();
  hRatio->Write();
}
