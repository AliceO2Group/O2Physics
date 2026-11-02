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

#include "ALICE3/Core/FlatLutRootWriter.h"

#include <TFile.h>
#include <TTree.h>

#include <algorithm>
#include <string>
#include <utility>

namespace o2::fastsim
{

FlatLutRootWriter::FlatLutRootWriter(const lutHeader_t& header, bool covarianceStoresSigma)
  : mHeader(header),
    mCovarianceStoresSigma(covarianceStoresSigma)
{
}

void FlatLutRootWriter::reserve(std::size_t n)
{
  mRows.reserve(n);
}

void FlatLutRootWriter::addEntry(int entryIndex,
                                 int nchBin,
                                 int radBin,
                                 int etaBin,
                                 int ptBin,
                                 const lutEntry_t& entry)
{
  FlatLutRootRow row;
  row.entryIndex = entryIndex;
  row.nchBin = nchBin;
  row.radBin = radBin;
  row.etaBin = etaBin;
  row.ptBin = ptBin;

  row.nch = entry.nch;
  row.eta = entry.eta;
  row.pt = entry.pt;

  row.valid = entry.valid;
  row.eff = entry.eff;
  row.eff2 = entry.eff2;
  row.itof = entry.itof;
  row.otof = entry.otof;

  std::copy(entry.covm.begin(), entry.covm.end(), row.covm.begin());
  std::copy(entry.eigval.begin(), entry.eigval.end(), row.eigval.begin());
  row.eigvec = entry.eigvec;
  row.eiginv = entry.eiginv;

  mRows.push_back(std::move(row));
}

bool FlatLutRootWriter::writeToFile(const std::string& rootFileName) const
{
  TFile output(rootFileName.c_str(), "RECREATE");
  if (output.IsZombie()) {
    return false;
  }

  int version = mHeader.version;
  int pdg = mHeader.pdg;
  float mass = mHeader.mass;
  float field = mHeader.field;

  int nch_nbins = mHeader.nchmap.nbins;
  float nch_min = mHeader.nchmap.min;
  float nch_max = mHeader.nchmap.max;
  bool nch_log = mHeader.nchmap.log;

  int rad_nbins = mHeader.radmap.nbins;
  float rad_min = mHeader.radmap.min;
  float rad_max = mHeader.radmap.max;
  bool rad_log = mHeader.radmap.log;

  int eta_nbins = mHeader.etamap.nbins;
  float eta_min = mHeader.etamap.min;
  float eta_max = mHeader.etamap.max;
  bool eta_log = mHeader.etamap.log;

  int pt_nbins = mHeader.ptmap.nbins;
  float pt_min = mHeader.ptmap.min;
  float pt_max = mHeader.ptmap.max;
  bool pt_log = mHeader.ptmap.log;

  bool covarianceStoresSigma = mCovarianceStoresSigma;

  TTree headerTree("lutHeader", "LUT header");
  headerTree.Branch("version", &version, "version/I");
  headerTree.Branch("pdg", &pdg, "pdg/I");
  headerTree.Branch("mass", &mass, "mass/F");
  headerTree.Branch("field", &field, "field/F");

  headerTree.Branch("nch_nbins", &nch_nbins, "nch_nbins/I");
  headerTree.Branch("nch_min", &nch_min, "nch_min/F");
  headerTree.Branch("nch_max", &nch_max, "nch_max/F");
  headerTree.Branch("nch_log", &nch_log, "nch_log/O");

  headerTree.Branch("rad_nbins", &rad_nbins, "rad_nbins/I");
  headerTree.Branch("rad_min", &rad_min, "rad_min/F");
  headerTree.Branch("rad_max", &rad_max, "rad_max/F");
  headerTree.Branch("rad_log", &rad_log, "rad_log/O");

  headerTree.Branch("eta_nbins", &eta_nbins, "eta_nbins/I");
  headerTree.Branch("eta_min", &eta_min, "eta_min/F");
  headerTree.Branch("eta_max", &eta_max, "eta_max/F");
  headerTree.Branch("eta_log", &eta_log, "eta_log/O");

  headerTree.Branch("pt_nbins", &pt_nbins, "pt_nbins/I");
  headerTree.Branch("pt_min", &pt_min, "pt_min/F");
  headerTree.Branch("pt_max", &pt_max, "pt_max/F");
  headerTree.Branch("pt_log", &pt_log, "pt_log/O");

  headerTree.Branch("covarianceStoresSigma", &covarianceStoresSigma, "covarianceStoresSigma/O");
  headerTree.Fill();

  TTree entryTree("lutEntries", "LUT entries");

  int entryIndex = -1;
  int nchBin = -1;
  int radBin = -1;
  int etaBin = -1;
  int ptBin = -1;

  float nch = 0.f;
  float eta = 0.f;
  float pt = 0.f;

  bool valid = false;
  float eff = 0.f;
  float eff2 = 0.f;
  float itof = 0.f;
  float otof = 0.f;

  float covm[kNumCovarianceTerms] = {0.f};
  float eigval[kNumEigenModes] = {0.f};
  float eigvec[kNumEigenModes * kNumEigenModes] = {0.f};
  float eiginv[kNumEigenModes * kNumEigenModes] = {0.f};

  entryTree.Branch("entryIndex", &entryIndex, "entryIndex/I");
  entryTree.Branch("nchBin", &nchBin, "nchBin/I");
  entryTree.Branch("radBin", &radBin, "radBin/I");
  entryTree.Branch("etaBin", &etaBin, "etaBin/I");
  entryTree.Branch("ptBin", &ptBin, "ptBin/I");

  entryTree.Branch("nch", &nch, "nch/F");
  entryTree.Branch("eta", &eta, "eta/F");
  entryTree.Branch("pt", &pt, "pt/F");

  entryTree.Branch("valid", &valid, "valid/O");
  entryTree.Branch("eff", &eff, "eff/F");
  entryTree.Branch("eff2", &eff2, "eff2/F");
  entryTree.Branch("itof", &itof, "itof/F");
  entryTree.Branch("otof", &otof, "otof/F");

  entryTree.Branch("covm", covm, "covm[15]/F");
  entryTree.Branch("eigval", eigval, "eigval[5]/F");
  entryTree.Branch("eigvec", eigvec, "eigvec[25]/F");
  entryTree.Branch("eiginv", eiginv, "eiginv[25]/F");

  for (const auto& row : mRows) {
    entryIndex = row.entryIndex;
    nchBin = row.nchBin;
    radBin = row.radBin;
    etaBin = row.etaBin;
    ptBin = row.ptBin;

    nch = row.nch;
    eta = row.eta;
    pt = row.pt;

    valid = row.valid;
    eff = row.eff;
    eff2 = row.eff2;
    itof = row.itof;
    otof = row.otof;

    std::copy(row.covm.begin(), row.covm.end(), covm);
    std::copy(row.eigval.begin(), row.eigval.end(), eigval);
    for (int i = 0; i < kNumEigenModes; ++i) {
      for (int j = 0; j < kNumEigenModes; ++j) {
        eigvec[kNumEigenModes * i + j] = row.eigvec[i][j];
        eiginv[kNumEigenModes * i + j] = row.eiginv[i][j];
      }
    }

    entryTree.Fill();
  }

  headerTree.Write();
  entryTree.Write();
  output.Close();

  return true;
}

} // namespace o2::fastsim
