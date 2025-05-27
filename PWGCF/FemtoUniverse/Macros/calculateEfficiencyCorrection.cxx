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

/// \file calculateEfficiencyCorrection.cxx
/// \brief Macro for calculating efficiency corrections based on 3D histograms
/// \author Dawid Karpiński, WUT Warsaw, dawid.karpinski@cern.ch

#include <TFile.h>
#include <TFolder.h>
#include <TGrid.h>
#include <TGridResult.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TSystem.h>
#include <unistd.h>

#include <cassert>
#include <cmath>
#include <filesystem> // NOLINT
#include <iostream>
#include <format>
#include <string>

namespace fs = std::filesystem;

auto* getHistogram(TFile* file, const std::string& name)
{
  return dynamic_cast<TH3*>(file->Get(name.c_str()));
}

auto* projectHistogram(TH3* hist, const std::string& projection)
{
  return hist->Project3D(projection.c_str());
}

template <typename H>
H* cloneHistogram(H* hist, const std::string& name)
{
  return dynamic_cast<H*>(hist->Clone(name.c_str()));
}

auto forEachBin(TH1* hist, auto func) -> void
{
  if (hist->GetDimension() == 1) {
    for (auto x{1}; x <= hist->GetNbinsX(); ++x) {
      func(x, 0, 0);
    }
  } else if (hist->GetDimension() == 2) {
    for (auto x{1}; x <= hist->GetNbinsX(); ++x) {
      for (auto y{1}; y <= hist->GetNbinsY(); ++y) {
        func(x, y, 0);
      }
    }
  } else if (hist->GetDimension() == 3) {
    for (auto x{1}; x <= hist->GetNbinsX(); ++x) {
      for (auto y{1}; y <= hist->GetNbinsY(); ++y) {
        for (auto z{1}; z <= hist->GetNbinsZ(); ++z) {
          func(x, y, z);
        }
      }
    }
  } else {
    assert(false && "should not happen");
  }
}

auto setAxisTitles(TH1* hist, const std::string& projection) -> void
{
  auto* xAxis{hist->GetXaxis()};
  auto* yAxis{hist->GetYaxis()};
  auto* zAxis{hist->GetZaxis()};

  xAxis->SetTitle("#it{p}_{T} (GeV/#it{c})");

  if (hist->GetDimension() == 2) {
    if (projection == "yx") {
      yAxis->SetTitle("#it{#eta}");
    } else if (projection == "zx") {
      yAxis->SetTitle("mult");
    }
  } else if (hist->GetDimension() == 3) {
    yAxis->SetTitle("#it{#eta}");
    zAxis->SetTitle("mult");
  }
}

auto calculateEfficiencyCorrection(const fs::path& resultsPath, const fs::path& histPath, const std::string& projection) -> void
{
  assert(!resultsPath.empty() && !histPath.empty());
  if (projection != "" && projection != "x" && projection != "yx" && projection != "zx") {
    std::cerr << "Error: projection must be one of: x, yx, zx\n";
    std::exit(1);
    return;
  }

  auto isAlien{false};
  if (resultsPath.string().starts_with("alien://")) {
    TGrid::Connect("alien://");
    isAlien = true;
  }

  auto* resultFile{TFile::Open(resultsPath.c_str())};
  assert(resultFile != nullptr && !resultFile->IsZombie());

  using namespace std::chrono;
  auto now{duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count()};
  auto outputPath{isAlien ? std::filesystem::current_path() : resultsPath.parent_path()};
  outputPath /= std::format("{}-effcor-{}.root", resultsPath.stem().string(), now);

  auto* outputFile{TFile::Open(outputPath.c_str(), "RECREATE")};
  assert(outputFile != nullptr && !outputFile->IsZombie());

  auto* histTruthBase{getHistogram(resultFile, histPath / "hMCTruth")};
  auto* histPrimaryBase{getHistogram(resultFile, histPath / "hPrimary")};
  auto* histSecondaryBase{getHistogram(resultFile, histPath / "hSecondary")};

  assert(histTruthBase);
  assert(histPrimaryBase);
  assert(histSecondaryBase);

  TH1* histTruth{histTruthBase};
  TH1* histPrimary{histPrimaryBase};
  TH1* histSecondary{histSecondaryBase};

  if (projection != "") {
    histPrimary = projectHistogram(histPrimaryBase, projection);
    histSecondary = projectHistogram(histSecondaryBase, projection);
    histTruth = projectHistogram(histTruthBase, projection);
  }

  auto* histTotal{cloneHistogram(histPrimary, "hTotal")};
  histTotal->Add(histSecondary);

  auto* histEfficiency{cloneHistogram(histPrimary, "hEfficiency")};
  histEfficiency->Reset();
  setAxisTitles(histEfficiency, projection);

  auto* histWeights{cloneHistogram(histPrimary, "hWeights")};
  histWeights->Reset();
  setAxisTitles(histWeights, projection);

  auto* histCont{cloneHistogram(histPrimary, "hCont")};
  histCont->Reset();
  setAxisTitles(histCont, projection);

  forEachBin(histPrimary, [&](int x, int y, int z) {
    auto primVal{histPrimary->GetBinContent(x, y, z)};
    auto primErr{histPrimary->GetBinError(x, y, z)};

    auto secVal{histSecondary->GetBinContent(x, y, z)};
    auto secErr{histSecondary->GetBinError(x, y, z)};

    auto truthVal{histTruth->GetBinContent(x, y, z)};
    auto truthErr{histTruth->GetBinError(x, y, z)};

    auto effVal{0.};
    auto effErr{0.};
    if (truthVal > 0) {
      effVal = primVal / truthVal;
      effErr = std::sqrt(std::pow(primErr / truthVal, 2) + std::pow((primVal * truthErr / std::pow(truthVal, 2)), 2));
    }

    histEfficiency->SetBinContent(x, y, z, effVal);
    histEfficiency->SetBinError(x, y, z, effErr);

    auto totalVal{primVal + secVal};
    auto totalErr{std::hypot(primErr, secErr)};

    auto contVal{0.};
    auto contErr{0.};
    if (totalVal > 0) {
      contVal = secVal / totalVal;
      contErr = std::sqrt(std::pow(secErr / totalVal, 2) + std::pow((secVal * totalErr / std::pow(totalVal, 2)), 2));
    }

    histCont->SetBinContent(x, y, z, contVal);
    histCont->SetBinError(x, y, z, contErr);

    auto weightVal{0.};
    auto weightErr{0.};
    if (effVal > 0) {
      weightVal = (1 - contVal) / effVal;
      weightErr = std::sqrt(std::pow(contErr / effVal, 2) + std::pow((1 - contVal) * effErr / std::pow(effVal, 2), 2));
    }

    histWeights->SetBinContent(x, y, z, weightVal);
    histWeights->SetBinError(x, y, z, weightErr);
  });

  outputFile->WriteTObject(histEfficiency);
  outputFile->WriteTObject(histCont);
  outputFile->WriteTObject(histWeights);

  outputFile->Close();
  resultFile->Close();
}

auto printUsage(const char* name) -> void
{
  std::cerr << "Usage: " << name << "\n"
            << "  -f <path to results ROOT file>\n"
            << "  -d <path to directory within file>\n"
            << "  -p <projection> [optional, default: no projection, 3D histogram]\n"
            << "    Available projections:\n"
            << "      x - projection onto pT axis (1D histogram)\n"
            << "      yx - projection onto pT, eta (2D histogram)\n"
            << "      zx - projection onto pT, mult (2D histogram)\n";
}

int main(int argc, char** argv)
{
  std::string results{""}, hist{""}, proj{""};
  auto flag{0};

  while ((flag = getopt(argc, argv, "f:d:p:")) != -1) {
    switch (flag) {
      case 'f':
        results = optarg;
        break;
      case 'd':
        hist = optarg;
        break;
      case 'p':
        proj = optarg;
        break;
      default:
        printUsage(argv[0]);
        return 1;
    }
  }

  if (results.empty() || hist.empty()) {
    printUsage(argv[0]);
    return 1;
  }

  calculateEfficiencyCorrection(results, hist, proj);
  return 0;
}
