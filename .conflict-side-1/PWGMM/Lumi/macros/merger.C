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
/// \author Junlee Kim (jikim1290@gmail.com)
/// \since November 2021
// merge all trees into one tree to construct event by event log-likelihood
// function.

void merger()
{
  TFile* fin = new TFile("../data/PbPb_NF/AnalysisResults_trees.root", "read");
  TTree* tadd;

  TFile* fout = new TFile("mergedOutput.root", "recreate");
  TTree* tin = new TTree("EventInfo_merged", "");

  TList* list = fin->GetListOfKeys();
  TIter next(fin->GetListOfKeys());

  ULong64_t fTimeStamp;
  double fVertexX;
  double fVertexY;
  double fVertexZ;
  double fVertexXY;
  double fVertexXX;
  double fVertexYY;
  double fVertexChi2;
  int fNContrib;

  tin->Branch("fTimeStamp", &fTimeStamp, "fTimeStamp/l");
  tin->Branch("fVertexX", &fVertexX, "fVertexX/D");
  tin->Branch("fVertexY", &fVertexY, "fVertexY/D");
  tin->Branch("fVertexZ", &fVertexZ, "fVertexZ/D");

  tin->Branch("fVertexXY", &fVertexXY, "fVertexXY/D");
  tin->Branch("fVertexXX", &fVertexXX, "fVertexXX/D");
  tin->Branch("fVertexYY", &fVertexYY, "fVertexYY/D");

  tin->Branch("fVertexChi2", &fVertexChi2, "fVertexChi2/D");
  tin->Branch("fNContrib", &fNContrib, "fNContrib/I");

  while (TDirectoryFile* dir = (TDirectoryFile*)next()) {
    tadd = (TTree*)fin->Get(Form("%s/O2eventinfo", dir->GetName()));

    tadd->SetBranchAddress("fTimeStamp", &fTimeStamp);

    tadd->SetBranchAddress("fVertexX", &fVertexX);
    tadd->SetBranchAddress("fVertexY", &fVertexY);
    tadd->SetBranchAddress("fVertexZ", &fVertexZ);

    tadd->SetBranchAddress("fVertexXY", &fVertexXY);
    tadd->SetBranchAddress("fVertexXX", &fVertexXX);
    tadd->SetBranchAddress("fVertexYY", &fVertexYY);

    tadd->SetBranchAddress("fVertexChi2", &fVertexChi2);
    tadd->SetBranchAddress("fNContrib", &fNContrib);

    for (int i = 0; i < tadd->GetEntries(); i++) {
      tadd->GetEntry(i);
      tin->Fill();
    }
  }
  fout->cd();
  tin->Write();
}
