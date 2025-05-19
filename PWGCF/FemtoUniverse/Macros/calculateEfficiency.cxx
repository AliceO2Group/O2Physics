#include <TFile.h>
#include <TFolder.h>
#include <TGrid.h>
#include <TGridResult.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TSystem.h>
#include <cassert>
#include <cmath>
#include <filesystem>
#include <format>
#include <string>

namespace fs = std::filesystem;

auto* getHistogram(TFile* file, const std::string& name, const std::string& projection)
{
  auto* hist{dynamic_cast<TH3*>(file->Get(name.c_str()))};
  if (projection == "none") {
    return hist;
  }
  return dynamic_cast<TH3*>(hist->Project3D(projection.c_str()));
}

template <typename H>
H* cloneHistogram(H* hist, const std::string& name)
{
  return dynamic_cast<H*>(hist->Clone(name.c_str()));
}

void forEachBin(TH1* hist, auto func)
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

void calculateEfficiency(const fs::path& resultsPath, const fs::path& histPath, const std::string& projection)
{
  assert(projection == "x" || projection == "yx" || projection == "zx" || projection == "none");

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

  auto* histTruth{getHistogram(resultFile, histPath / "hMCTruth", projection)};
  assert(histTruth);

  auto* histPrimary{getHistogram(resultFile, histPath / "hPrimary", projection)};
  assert(histPrimary);

  auto* histSecondary{getHistogram(resultFile, histPath / "hSecondary", projection)};
  assert(histSecondary);

  auto* histTotal{cloneHistogram(histPrimary, "hTotal")};
  histTotal->Add(histSecondary);

  auto* histEfficiency{cloneHistogram(histPrimary, "hEfficiency")};
  histEfficiency->Reset();

  auto* histWeights{cloneHistogram(histPrimary, "hWeights")};
  histWeights->Reset();

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
  outputFile->WriteTObject(histWeights);

  outputFile->Close();
  resultFile->Close();
}

int main(int argc, char** argv)
{
  assert(argc == 4);
  calculateEfficiency(argv[1], argv[2], argv[3]);
  return 0;
}
