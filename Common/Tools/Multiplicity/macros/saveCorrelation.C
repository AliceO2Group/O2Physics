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
//
/// \file saveCorrelation.C
/// \brief This file provides a simple macro to convert TGlauberMC output tuples
/// into the 2D correlation histogram necessary for the ALICE machinery
/// that performs Glauber + NBD fits.
/// \author ALICE

#include <iostream>

/// @brief function to save Npart x Ncoll correlation to file for glauber fits
/// @param filename input TGlauberMC ntuple file
/// @param outputFile output file for Npart x Ncoll correlation TH2D
void saveCorrelation(TString filename = "gmc-PbPb-snn68.21-md0.40-nd-1.0-rc1-smax99.0.root", TString outputFile = "basehistos.root")
{
  TFile* fin = new TFile(filename.Data(), "READ");
  TNtuple* ntup = (TNtuple*)fin->Get("nt_Pb_Pb");

  // try other Pb nuclear profiles in case "Pb" - "Pb" not found
  if (!ntup) {
    ntup = (TNtuple*)fin->Get("nt_Pbpn_Pbpn");
  }
  if (!ntup) {
    ntup = (TNtuple*)fin->Get("nt_PbpnVar1_PbpnVar1");
  }
  if (!ntup) {
    ntup = (TNtuple*)fin->Get("nt_PbpnVar2_PbpnVar2");
  }
  if (!ntup) {
    ntup = (TNtuple*)fin->Get("nt_PbpnVar3_PbpnVar3");
  }
  if (!ntup) {
    ntup = (TNtuple*)fin->Get("nt_PbpnVar4_PbpnVar4");
  }

  if (!ntup) {
    cout << "No tree found!" << endl;
    return;
  }

  cout << "Glauber tree entries: " << ntup->GetEntries() << endl;

  TFile* fout = new TFile(outputFile.Data(), "RECREATE");

  // 2D correlation plot necessary for Glauber + NBD fitting
  // The provided range should be enough for Pb-Pb
  TH2D* hNpNc = new TH2D("hNpNc", "", 500, -0.5, 499.5, 2500, -0.5, 2499.5);

  // let's draw this on screen for inspection
  TCanvas* c1 = new TCanvas("c1", "", 800, 600);
  c1->SetTicks(1, 1);
  c1->SetLogz();
  ntup->Draw("Ncoll:Npart>>hNpNc", "", "colz");
  fout->Write();
}
