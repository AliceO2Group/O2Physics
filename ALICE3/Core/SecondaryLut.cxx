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

/// \file SecondaryLut.cxx
/// \brief ROOT-persistent LUT for secondary particles.

#include "ALICE3/Core/SecondaryLut.h"

#include <TBrowser.h>
#include <TFile.h>
#include <TH1.h>
#include <TNamed.h>
#include <TString.h>
#include <TTree.h>

#include <RtypesCore.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <memory>
#include <vector>

namespace o2::fastsim
{

namespace
{
constexpr std::array<const char*, SecondaryLut::kNDims> kDimNames = {"nch", "rad", "eta", "pt", "nHits"};
}

int SecondaryLutAxis::find(float val) const
{
  if (log) {
    val = val > 0.f ? std::log10(val) : min;
  }
  const int bin = static_cast<int>(std::floor((val - min) / (max - min) * nbins));
  return std::clamp(bin, 0, nbins - 1);
}

float SecondaryLutAxis::lowEdge(int bin) const
{
  const float v = min + bin * (max - min) / nbins;
  return log ? std::pow(10.f, v) : v;
}

float SecondaryLutAxis::upEdge(int bin) const { return lowEdge(bin + 1); }

float SecondaryLutAxis::center(int bin) const
{
  const float v = min + (bin + 0.5f) * (max - min) / nbins;
  return log ? std::pow(10.f, v) : v;
}

SecondaryLut::SecondaryLut(const char* name, int pdg, float mass, float field)
  : TNamed(name, "Secondary particle LUT"), mPdg(pdg), mMass(mass), mField(field)
{
}

void SecondaryLut::setAxis(Dim dim, int nbins, float min, float max, bool log)
{
  mAxes[dim] = {.nbins = std::max(nbins, 1), .min = min, .max = max, .log = log};
  mEntries.clear();
}

void SecondaryLut::init()
{
  std::size_t n = 1;
  for (const auto& a : mAxes) {
    n *= static_cast<std::size_t>(a.nbins);
  }
  mEntries.assign(n, SecondaryLutEntry{});
}

const SecondaryLutEntry* SecondaryLut::lookup(float nch, float rad, float eta, float pt, float nHits) const
{
  if (mEntries.empty()) {
    return nullptr;
  }
  return &mEntries[index(mAxes[kNch].find(nch), mAxes[kRad].find(rad), mAxes[kEta].find(eta),
                         mAxes[kPt].find(pt), mAxes[kNHits].find(nHits))];
}

bool SecondaryLut::writeToFile(const char* fileName, const char* option) const
{
  std::unique_ptr<TFile> f(TFile::Open(fileName, option));
  if (!f || f->IsZombie()) {
    return false;
  }
  const auto written = f->WriteObject(this, GetName());
  return written > 0 && !f->IsZombie();
}

SecondaryLut* SecondaryLut::readFromFile(const char* fileName, const char* objName)
{
  std::unique_ptr<TFile> f(TFile::Open(fileName, "READ"));
  if (!f || f->IsZombie()) {
    return nullptr;
  }
  auto* in = f->Get<SecondaryLut>(objName);
  return in ? new SecondaryLut(*in) : nullptr;
}

void SecondaryLut::unravel(std::size_t idx, std::array<int, kNDims>& bins) const
{
  // Reverse of index(): pt fastest, then hits, eta, rad, nch
  static constexpr std::array<Dim, kNDims> order = {kPt, kNHits, kEta, kRad, kNch};
  for (const Dim d : order) {
    bins[d] = static_cast<int>(idx % mAxes[d].nbins);
    idx /= mAxes[d].nbins;
  }
}

TTree* SecondaryLut::toTree(const char* treeName) const
{
  auto* tree = new TTree(treeName, "Secondary LUT flat view");
  tree->SetDirectory(nullptr);
  std::array<int, kNDims> bins{};
  std::array<float, kNDims> centers{};
  SecondaryLutEntry e;
  for (int d = 0; d < kNDims; ++d) {
    tree->Branch(Form("%sBin", kDimNames[d]), &bins[d]);
    tree->Branch(kDimNames[d], &centers[d]);
  }
  tree->Branch("valid", &e.valid);
  tree->Branch("eff", &e.eff);
  tree->Branch("eff2", &e.eff2);
  tree->Branch("itof", &e.itof);
  tree->Branch("otof", &e.otof);
  tree->Branch("covm", e.covm.data(), "covm[15]/F");
  tree->Branch("eigval", e.eigval.data(), "eigval[5]/F");
  tree->Branch("eigvec", e.eigvec.data(), "eigvec[25]/F");
  tree->Branch("eiginv", e.eiginv.data(), "eiginv[25]/F");
  for (std::size_t i = 0; i < mEntries.size(); ++i) {
    unravel(i, bins);
    for (int d = 0; d < kNDims; ++d) {
      centers[d] = mAxes[d].center(bins[d]);
    }
    e = mEntries[i];
    tree->Fill();
  }
  tree->ResetBranchAddresses();
  return tree;
}

TH1D* SecondaryLut::projectEfficiency(Dim dim, const char* histName) const
{
  const auto& ax = mAxes[dim];
  std::vector<double> edges(ax.nbins + 1);
  for (int i = 0; i <= ax.nbins; ++i) {
    edges[i] = ax.lowEdge(i);
  }
  auto* h = new TH1D(histName ? histName : Form("%s_eff_vs_%s", GetName(), kDimNames[dim]),
                     Form("Average efficiency vs %s;%s;<eff>", kDimNames[dim], kDimNames[dim]),
                     ax.nbins, edges.data());
  h->SetDirectory(nullptr);
  std::vector<double> sum(ax.nbins, 0.), cnt(ax.nbins, 0.);
  std::array<int, kNDims> bins{};
  for (std::size_t i = 0; i < mEntries.size(); ++i) {
    if (!mEntries[i].valid) {
      continue;
    }
    unravel(i, bins);
    sum[bins[dim]] += mEntries[i].eff;
    cnt[bins[dim]] += 1.;
  }
  for (int b = 0; b < ax.nbins; ++b) {
    h->SetBinContent(b + 1, cnt[b] > 0. ? sum[b] / cnt[b] : 0.);
  }
  return h;
}

void SecondaryLut::Print(Option_t*) const
{
  std::printf("SecondaryLut '%s': pdg=%d mass=%g field=%g entries=%zu (%.2f MB)\n", GetName(), mPdg, mMass, mField,
              mEntries.size(), mEntries.size() * sizeof(SecondaryLutEntry) / 1.e6);
  for (int d = 0; d < kNDims; ++d) {
    const auto& a = mAxes[d];
    std::printf("  %-6s nbins=%-4d [%g, %g]%s\n", kDimNames[d], a.nbins, a.min, a.max, a.log ? " (log10)" : "");
  }
  std::size_t nvalid = 0;
  for (const auto& e : mEntries) {
    nvalid += static_cast<std::size_t>(e.valid);
  }
  std::printf("  valid bins: %zu\n", nvalid);
}

void SecondaryLut::Browse(TBrowser* b)
{
  if (!b || mEntries.empty()) {
    return;
  }
  for (int d = 0; d < kNDims; ++d) {
    auto* h = projectEfficiency(static_cast<Dim>(d));
    b->Add(h, h->GetName());
  }
  b->Add(toTree(), "secondaryLut");
}

} // namespace o2::fastsim
