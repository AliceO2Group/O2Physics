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

#ifndef PWGCF_GENERICFRAMEWORK_CORE_PLOTTEMPLATESTEERER_H_
#define PWGCF_GENERICFRAMEWORK_CORE_PLOTTEMPLATESTEERER_H_

#include "MyCanvas.h"
#include "PlotYRangeUtils.h"
#include "RVersion.h"

#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TKey.h"

#include <cmath>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

class PlotTemplateSteerer
{
 public:
  using ObjType = MyCanvas::ObjType;

  // Add ratios, pair ratios, etc.
  enum class TemplateMode {
    TwoMeasurements,  // 2 inputs: same color, closed/open markers
    Paired,           // (0,1), (2,3), ... share color and closed/open marker
    Distinct,         // each series gets its own color/marker
    DistinctRatio,    // first pad: distinct, second pad: ratio to first series
    PairedRatio,      // first pad: paired, second pad: (0/1), (2/3), ...
    PairedDifference, // first pad: paired, second pad: (0-1), (2-3), ...
    PairedRelDev      // first pad: paired, second pad 1.0-|0-1|/0, 1.0-|2-3|/2
  };
  enum class ComparisonMode {
    Ratio,
    Difference,
    RelativeDeviation
  };
  struct SeriesInput {
    ObjType stat;
    ObjType sys{};
    TString legend{};
    const char* legendDrawOpt = "pe";
  };

  struct SeriesStyle {
    int color = kBlack;
    int marker = kFullCircle;
  };

  PlotTemplateSteerer()
    : mFallbackColorPalette{kBlack, kRed + 1, kBlue + 1, kGreen + 2, kOrange + 7, kMagenta + 2, kCyan + 2},
      mPairedMarkers{{kFullCircle, kOpenCircle},
                     {kFullSquare, kOpenSquare},
                     {kFullDiamond, kOpenDiamond},
                     {kFullDoubleDiamond, kOpenDoubleDiamond},
                     {kFullCross, kOpenCross},
                     {kFullStar, kOpenStar},
                     {kFullCrossX, kOpenCrossX},
                     {kFullFourTrianglesPlus, kOpenFourTrianglesPlus}},
      mDistinctMarkers{
        kFullCircle, kFullSquare, kFullDiamond, kFullDoubleDiamond, kFullCross, kFullStar, kFullCrossX, kFullFourTrianglesPlus,
        kOpenCircle, kOpenSquare, kOpenDiamond, kOpenDoubleDiamond, kOpenCross, kOpenStar, kOpenCrossX, kOpenFourTrianglesPlus} {}

  void setFallbackColorPalette(const std::vector<int>& colors)
  {
    if (colors.empty()) {
      throw std::invalid_argument("PlotTemplateSteerer: fallback color palette cannot be empty");
    }
    mFallbackColorPalette = colors;
  }

  void setColorPalette(const std::vector<int>& colors) { setFallbackColorPalette(colors); }

  void setPairedMarkers(const std::vector<std::pair<int, int>>& pairedMarkers)
  {
    if (pairedMarkers.empty()) {
      throw std::invalid_argument("PlotTemplateSteerer: paired marker list cannot be empty");
    }
    mPairedMarkers = pairedMarkers;
  }

  void setDistinctMarkers(const std::vector<int>& markers)
  {
    if (markers.empty()) {
      throw std::invalid_argument("PlotTemplateSteerer: distinct marker list cannot be empty");
    }
    mDistinctMarkers = markers;
  }

  std::vector<SeriesStyle> buildStyles(size_t nSeries, TemplateMode mode, int colorOffset = 0, int markerOffset = 0) const
  {
    if (nSeries == 0) {
      return {};
    }
    const auto colorPalette = selectColorPalette(nSeries, mode);
    std::vector<SeriesStyle> styles;
    styles.reserve(nSeries);
    for (size_t i = 0; i < nSeries; ++i) {
      styles.push_back(styleForIndex(i, mode, colorOffset, markerOffset, colorPalette));
    }
    return styles;
  }

  void addSeries(MyCanvas& canvas,
                 const std::vector<SeriesInput>& series,
                 TemplateMode mode,
                 int pad = 1,
                 int colorOffset = 0,
                 int markerOffset = 0,
                 int legInd = 0,
                 const std::vector<int>& legendIndices = {}) const
  {
    mComparisonObjects.clear();
    if (series.empty()) {
      return;
    }
    if (mode == TemplateMode::TwoMeasurements && series.size() != 2) {
      throw std::invalid_argument("PlotTemplateSteerer: TwoMeasurements mode requires exactly 2 series");
    }
    if ((mode == TemplateMode::PairedRatio || mode == TemplateMode::PairedDifference) && (series.size() < 2 || (series.size() % 2) != 0)) {
      throw std::invalid_argument("PlotTemplateSteerer: PairedRatio and PairedDifference mode requires an even number of series >= 2");
    }
    if (!legendIndices.empty() && legendIndices.size() != series.size()) {
      throw std::invalid_argument("PlotTemplateSteerer: legend index list must be empty or match the number of series");
    }

    if (!canvas.isInitialised()) {
      if (mode == TemplateMode::DistinctRatio || mode == TemplateMode::PairedRatio || mode == TemplateMode::PairedDifference || mode == TemplateMode::PairedRelDev) {
        canvas.BuildRatioCanvas(0.45);
      } else {
        canvas.BuildSingleCanvas();
      }
    }

    auto styles = buildStyles(series.size(), mode, colorOffset, markerOffset);
    for (size_t i = 0; i < series.size(); ++i) {
      canvas.AddObject(series[i].stat, pad, styles[i].color, styles[i].marker, series[i].sys);
      const int seriesLegInd = legendIndexForSeries(i, legInd, legendIndices);
      if (seriesLegInd >= 0 && !series[i].legend.IsNull()) {
        std::visit(
          [&canvas, &series, i, seriesLegInd](auto&& arg) { canvas.AddLegEntry(arg, series[i].legend, series[i].legendDrawOpt, seriesLegInd); },
          series[i].stat);
      }
      if (mode == TemplateMode::DistinctRatio) {
        if (i == 0) {
          continue;
        }
        ComparisonMode comparisonMode = ComparisonMode::Ratio;
        if (mode == TemplateMode::PairedDifference)
          comparisonMode = ComparisonMode::Difference;
        if (mode == TemplateMode::PairedRelDev)
          comparisonMode = ComparisonMode::RelativeDeviation;
        ObjType compStat = makeComparisonObject(series[i].stat, series[0].stat, static_cast<int>(i), comparisonMode);
        mComparisonObjects.push_back(compStat);
        std::visit([](auto&& arg) {
          if (!arg) {
            return;
          }
        },
                   compStat);
        canvas.AddObject(compStat, 2, styles[i].color, styles[i].marker, {});
        if (legInd >= 0 && !series[i].legend.IsNull()) {
          std::visit(
            [&canvas, &series, i, legInd](auto&& arg) { canvas.AddLegEntry(arg, series[i].legend + TString("/") + series[0].legend, series[i].legendDrawOpt, legInd + 1); },
            series[i].stat);
        }
      }
      if (mode == TemplateMode::PairedRatio || mode == TemplateMode::PairedDifference || mode == TemplateMode::PairedRelDev) {
        if ((i % 2) == 1) {
          continue;
        }
        const size_t denIdx = i + 1;
        ComparisonMode comparisonMode = ComparisonMode::Ratio;
        if (mode == TemplateMode::PairedDifference)
          comparisonMode = ComparisonMode::Difference;
        if (mode == TemplateMode::PairedRelDev)
          comparisonMode = ComparisonMode::RelativeDeviation;
        ObjType compStat = makeComparisonObject(series[i].stat, series[denIdx].stat, static_cast<int>(i), comparisonMode);
        mComparisonObjects.push_back(compStat);
        std::visit([](auto&& arg) {
          if (!arg) {
            return;
          }
        },
                   compStat);
        canvas.AddObject(compStat, 2, styles[i].color, styles[i].marker, {});
        if (legInd >= 0 && !series[i].legend.IsNull() && !series[denIdx].legend.IsNull()) {
          std::visit(
            [&canvas, &series, i, denIdx, legInd, comparisonMode](auto&& arg) {
              TString compStr = (comparisonMode == ComparisonMode::Difference) ? TString("-") : TString("/");
              canvas.AddLegEntry(arg, series[i].legend + compStr + series[denIdx].legend, series[i].legendDrawOpt, legInd + 1);
            },
            series[i].stat);
        }
      }
    }
  }

  void addTwoMeasurements(MyCanvas& canvas,
                          const SeriesInput& first,
                          const SeriesInput& second,
                          int pad = 1,
                          int colorOffset = 0,
                          int markerOffset = 0,
                          int legInd = -1) const
  {
    addSeries(canvas, {first, second}, TemplateMode::TwoMeasurements, pad, colorOffset, markerOffset, legInd);
  }

  void applyAutoYLimits(MyCanvas& canvas,
                        const std::vector<SeriesInput>& series,
                        TemplateMode mode,
                        double referenceLine = 1.0) const
  {
    if (series.empty()) {
      return;
    }

    const auto topRange = plotutils::determineTopPadYLimits(series);
    if (topRange.valid) {
      canvas.SetYLimits(topRange.ymin, topRange.ymax, 1);
    }

    plotutils::ComparisonMode cmpMode = plotutils::ComparisonMode::None;
    if (mode == TemplateMode::PairedRatio) {
      cmpMode = plotutils::ComparisonMode::PairedRatio;
    } else if (mode == TemplateMode::PairedDifference) {
      cmpMode = plotutils::ComparisonMode::PairedDifference;
    }

    if (cmpMode != plotutils::ComparisonMode::None) {
      const auto bottomRange = plotutils::determineBottomPadYLimits(series, cmpMode, referenceLine);
      if (bottomRange.valid) {
        canvas.SetYLimits(bottomRange.ymin, bottomRange.ymax, 2);
      }
    }
  }

  void applyText(MyCanvas& c, const std::vector<std::string>& topleft, const std::string& label) const
  {
    double leftMargin = 0.10;
    double rightMargin = 0.05;
    double topMargin = 0.05;
    if (gPad) {
      leftMargin = gPad->GetLeftMargin();
      rightMargin = gPad->GetRightMargin();
      topMargin = gPad->GetTopMargin();
    }

    const double xLeft = leftMargin + 0.03;
    const double yStart = 1.0 - topMargin - 0.05;
    const double lineStep = 0.06;

    for (size_t i = 0; i < topleft.size(); ++i) {
      const double y = yStart - static_cast<double>(i) * lineStep;
      c.AddText(xLeft, y, topleft[i].c_str(), true);
    }

    if (!label.empty()) {
      const double xRight = 1.0 - rightMargin - 0.02;
      // c.AddText(xRight, yStart, label.c_str(), true);
      TLatex l;
      l.SetTextFont(43);
      l.SetTextSize(c.GetTextSize());
      l.SetTextAlign(31); // right-aligned in x, top-aligned in y
      l.DrawLatexNDC(xRight, yStart, label.c_str());
    }
  }

  void createSeriesFromFile(std::string filename)
  {
    std::unique_ptr<TFile> inFile(TFile::Open(filename.c_str(), "READ"));
    if (!inFile || inFile->IsZombie()) {
      throw std::runtime_error("PlotTemplateSteerer: could not open file: " + filename);
    }

    TIter nextKey(inFile->GetListOfKeys());
    TKey* key = nullptr;
    while ((key = static_cast<TKey*>(nextKey()))) {
      const std::string keyName = key->GetName();
      mSeriesInputs.push_back(getSeriesInput(filename, keyName));
    }
  }
  void createSeriesFromFiles(const std::vector<std::string>& filename, const std::string& qualifier)
  {
    for (const auto& file : filename) {
      mSeriesInputs.push_back(getSeriesInput(file, qualifier));
    }
  }
  void createSeriesFromFilesWithSystematics(const std::vector<std::string>& filename,
                                            const std::string& statQualifier,
                                            const std::string& sysQualifier,
                                            const std::vector<std::string>& legends = {})
  {
    if (legends.size() != filename.size() && !legends.empty()) {
      throw std::runtime_error("PlotTemplateSteerer: filenames and legends have different sizes");
    }

    for (size_t i = 0; i < filename.size(); ++i) {
      if (!legends.empty()) {
        mSeriesInputs.push_back(getSeriesInput(filename[i], statQualifier, sysQualifier, legends[i]));
      } else {
        mSeriesInputs.push_back(getSeriesInput(filename[i], statQualifier, sysQualifier, ""));
      }
    }
  }
  void createSeriesFromFiles(const std::vector<std::string>& filename, const std::vector<std::string>& qualifiers, const std::vector<std::string>& legends = {})
  {
    if (qualifiers.size() != filename.size())
      throw std::runtime_error("PlotTemplateSteerer: filenames and qualifiers have different sizes");
    if (legends.size() != filename.size() && !legends.empty())
      throw std::runtime_error("PlotTemplateSteerer: filenames and qualifiers have different sizes");

    for (size_t i = 0; i < filename.size(); ++i) {
      if (!legends.empty())
        mSeriesInputs.push_back(getSeriesInput(filename[i], qualifiers[i], legends[i]));
      else
        mSeriesInputs.push_back(getSeriesInput(filename[i], qualifiers[i]));
    }
  }
  void createSeriesFromFilesWithSystematics(const std::vector<std::string>& filename,
                                            const std::vector<std::string>& statQualifiers,
                                            const std::vector<std::string>& sysQualifiers,
                                            const std::vector<std::string>& legends = {})
  {
    if (statQualifiers.size() != filename.size() || sysQualifiers.size() != filename.size()) {
      throw std::runtime_error("PlotTemplateSteerer: filenames and qualifiers have different sizes");
    }
    if (legends.size() != filename.size() && !legends.empty()) {
      throw std::runtime_error("PlotTemplateSteerer: filenames and legends have different sizes");
    }

    for (size_t i = 0; i < filename.size(); ++i) {
      if (!legends.empty()) {
        mSeriesInputs.push_back(getSeriesInput(filename[i], statQualifiers[i], sysQualifiers[i], legends[i]));
      } else {
        mSeriesInputs.push_back(getSeriesInput(filename[i], statQualifiers[i], sysQualifiers[i], ""));
      }
    }
  }
  void createSeriesFromFilesWithSystematics(const std::vector<std::string>& statFilenames,
                                            const std::vector<std::string>& statQualifiers,
                                            const std::vector<std::string>& sysFilenames,
                                            const std::vector<std::string>& sysQualifiers,
                                            const std::vector<std::string>& legends = {})
  {
    if (statQualifiers.size() != statFilenames.size() || sysFilenames.size() != statFilenames.size() || sysQualifiers.size() != statFilenames.size()) {
      throw std::runtime_error("PlotTemplateSteerer: stat/sys filenames and qualifiers have different sizes");
    }
    if (legends.size() != statFilenames.size() && !legends.empty()) {
      throw std::runtime_error("PlotTemplateSteerer: filenames and legends have different sizes");
    }

    for (size_t i = 0; i < statFilenames.size(); ++i) {
      if (!legends.empty()) {
        mSeriesInputs.push_back(getSeriesInput(statFilenames[i], statQualifiers[i], sysFilenames[i], sysQualifiers[i], legends[i]));
      } else {
        mSeriesInputs.push_back(getSeriesInput(statFilenames[i], statQualifiers[i], sysFilenames[i], sysQualifiers[i], ""));
      }
    }
  }
  SeriesInput getSeriesInput(const std::string& filename, const std::string& objname, const std::string& legend = "")
  {
    std::unique_ptr<TFile> inFile(TFile::Open(filename.c_str(), "READ"));
    if (!inFile || inFile->IsZombie()) {
      throw std::runtime_error("PlotTemplateSteerer: could not open file: " + filename);
    }

    ObjType stat = getObjectFromFile(*inFile, filename, objname);
    return {stat, {}, (legend.empty()) ? objectName(stat) : legend, "lp"};
  }
  SeriesInput getSeriesInput(const std::string& filename, const std::string& statObjname, const std::string& sysObjname, const std::string& legend)
  {
    std::unique_ptr<TFile> inFile(TFile::Open(filename.c_str(), "READ"));
    if (!inFile || inFile->IsZombie()) {
      throw std::runtime_error("PlotTemplateSteerer: could not open file: " + filename);
    }

    ObjType stat = getObjectFromFile(*inFile, filename, statObjname);
    ObjType sys = getObjectFromFile(*inFile, filename, sysObjname);
    return {stat, sys, (legend.empty()) ? objectName(stat) : legend, "lp"};
  }
  SeriesInput getSeriesInput(const std::string& statFilename,
                             const std::string& statObjname,
                             const std::string& sysFilename,
                             const std::string& sysObjname,
                             const std::string& legend)
  {
    std::unique_ptr<TFile> statFile(TFile::Open(statFilename.c_str(), "READ"));
    if (!statFile || statFile->IsZombie()) {
      throw std::runtime_error("PlotTemplateSteerer: could not open file: " + statFilename);
    }
    std::unique_ptr<TFile> sysFile(TFile::Open(sysFilename.c_str(), "READ"));
    if (!sysFile || sysFile->IsZombie()) {
      throw std::runtime_error("PlotTemplateSteerer: could not open file: " + sysFilename);
    }

    ObjType stat = getObjectFromFile(*statFile, statFilename, statObjname);
    ObjType sys = getObjectFromFile(*sysFile, sysFilename, sysObjname);
    return {stat, sys, (legend.empty()) ? objectName(stat) : legend, "lp"};
  }

  const std::vector<SeriesInput>& getLoadedSeriesInputs() const { return mSeriesInputs; }
  void clearLoadedSeriesInputs() { mSeriesInputs.clear(); }
  const std::vector<ObjType>& getComparisonObjects() const { return mComparisonObjects; }
  void clearComparisonObjects() { mComparisonObjects.clear(); }

 private:
  static size_t offsetIndex(size_t base, int offset, size_t modulus)
  {
    const int64_t m = static_cast<int64_t>(modulus);
    const int64_t b = static_cast<int64_t>(base);
    const int64_t o = static_cast<int64_t>(offset);
    const int64_t wrapped = ((b + o) % m + m) % m;
    return static_cast<size_t>(wrapped);
  }

  static int legendIndexForSeries(size_t index, int defaultLegendIndex, const std::vector<int>& legendIndices)
  {
    if (legendIndices.empty()) {
      return defaultLegendIndex;
    }
    return legendIndices[index];
  }

  static TObject* findObjectInFile(TFile& inFile, const std::string& objname)
  {
    TObject* obj = inFile.Get(objname.c_str());
    if (obj) {
      return obj;
    }

    TIter nextKey(inFile.GetListOfKeys());
    TKey* key = nullptr;
    while ((key = static_cast<TKey*>(nextKey()))) {
      const std::string keyName = key->GetName();
      if (keyName.find(objname) != std::string::npos) {
        return inFile.Get(keyName.c_str());
      }
    }
    return nullptr;
  }

  static ObjType getObjectFromFile(TFile& inFile, const std::string& filename, const std::string& objname)
  {
    TObject* obj = findObjectInFile(inFile, objname);
    if (!obj) {
      throw std::runtime_error("PlotTemplateSteerer: object not found in file " + filename + ": " + objname);
    }
    if (auto* gAsymm = dynamic_cast<TGraphAsymmErrors*>(obj)) {
      auto* clone = static_cast<TGraphAsymmErrors*>(gAsymm->Clone());
      return clone;
    }
    if (auto* gErr = dynamic_cast<TGraphErrors*>(obj)) {
      auto* clone = static_cast<TGraphErrors*>(gErr->Clone());
      return clone;
    }
    if (auto* h1d = dynamic_cast<TH1D*>(obj)) {
      auto* clone = static_cast<TH1D*>(h1d->Clone());
      clone->SetDirectory(nullptr);
      return clone;
    }
    if (auto* h1 = dynamic_cast<TH1*>(obj)) {
      auto* clone = static_cast<TH1*>(h1->Clone());
      clone->SetDirectory(nullptr);
      return clone;
    }

    throw std::runtime_error("PlotTemplateSteerer: unsupported object type for " + filename + ":" + objname);
  }

  static TString objectName(const ObjType& obj)
  {
    return std::visit(
      [](auto* ptr) -> TString {
        if (!ptr) {
          return {};
        }
        return ptr->GetName();
      },
      obj);
  }

  std::vector<int> selectColorPalette(size_t nSeries, TemplateMode mode) const
  {
    const size_t nColorGroups = (mode == TemplateMode::Distinct || mode == TemplateMode::DistinctRatio) ? nSeries : (nSeries + 1) / 2;
    if (nColorGroups <= 6) {
      std::vector<int> palette;
      palette.reserve(6);
      for (int i = 0; i < 6; ++i) {
        palette.push_back(kP6Blue + i);
      }
      return palette;
    }
    if (nColorGroups <= 8) {
      std::vector<int> palette;
      palette.reserve(8);
      for (int i = 0; i < 8; ++i) {
        palette.push_back(kP8Blue + i);
      }
      return palette;
    }
    std::vector<int> palette;
    palette.reserve(10);
    for (int i = 0; i < 10; ++i) {
      palette.push_back(kP10Blue + i);
    }
    return palette;
  }

  SeriesStyle styleForIndex(size_t index,
                            TemplateMode mode,
                            int colorOffset,
                            int markerOffset,
                            const std::vector<int>& colorPalette) const
  {
    if (mode == TemplateMode::TwoMeasurements || mode == TemplateMode::Paired || mode == TemplateMode::PairedRatio || mode == TemplateMode::PairedDifference) {
      const size_t pairIndex = index / 2;
      const bool isOpen = (index % 2) == 1;
      const auto markerIndex = offsetIndex(pairIndex, markerOffset, mPairedMarkers.size());
      const auto colorIndex = offsetIndex(pairIndex, colorOffset, colorPalette.size());
      const auto& markerPair = mPairedMarkers[markerIndex];
      const int marker = isOpen ? markerPair.second : markerPair.first;
      const int color = colorPalette[colorIndex];
      return {color, marker};
    }

    const auto markerIndex = offsetIndex(index, markerOffset, mDistinctMarkers.size());
    const auto colorIndex = offsetIndex(index, colorOffset, colorPalette.size());
    const int marker = mDistinctMarkers[markerIndex];
    const int color = colorPalette[colorIndex];
    return {color, marker};
  }

  static TH1* asHistPtr(const ObjType& obj)
  {
    return std::visit(
      [](auto* ptr) -> TH1* {
        using T = std::decay_t<decltype(ptr)>;
        if constexpr (std::is_same_v<T, TH1*> || std::is_same_v<T, TH1D*>) {
          return ptr;
        }
        return nullptr;
      },
      obj);
  }

  static bool isGraphType(const ObjType& obj)
  {
    return std::visit(
      [](auto* ptr) {
        using T = std::decay_t<decltype(ptr)>;
        return std::is_same_v<T, TGraphErrors*> || std::is_same_v<T, TGraphAsymmErrors*>;
      },
      obj);
  }
  static TGraphAsymmErrors* histToGraph(const TH1* h)
  {
    auto* g = new TGraphAsymmErrors(h);
    g->SetName(Form("%s_graph", h->GetName()));
    g->GetXaxis()->SetTitle(h->GetXaxis()->GetTitle());
    return g;
  }
  static TH1D* graphToHist(const TGraph* g,
                           const TH1* refHist,
                           const std::string& name)
  {
    if (!g || !refHist) {
      throw std::runtime_error("graphToHist: null input");
    }

    const int nBins = refHist->GetNbinsX();

    auto* h = dynamic_cast<TH1D*>(
      refHist->Clone(name.c_str()));

    if (!h) {
      throw std::runtime_error("graphToHist: failed to clone reference histogram");
    }

    h->Reset("ICES"); // keep binning, clear contents
    h->SetDirectory(nullptr);

    const int n = g->GetN();

    for (int i = 0; i < n; ++i) {
      double x = 0.0, y = 0.0;
      g->GetPoint(i, x, y);

      int bin = h->FindBin(x);

      if (bin < 1 || bin > nBins) {
        continue; // out of range → ignore
      }

      h->SetBinContent(bin, y);

      // --- uncertainty handling ---
      double ey = 0.0;

      if (auto* ge = dynamic_cast<const TGraphErrors*>(g)) {
        ey = ge->GetErrorY(i);
      } else if (auto* gae = dynamic_cast<const TGraphAsymmErrors*>(g)) {
        ey = 0.5 * (gae->GetErrorYlow(i) + gae->GetErrorYhigh(i));
      }

      h->SetBinError(bin, ey);
    }

    return h;
  }
  static ObjType makeComparisonObject(const ObjType& numerator, const ObjType& denominator, int index, ComparisonMode comparisonMode = ComparisonMode::Ratio)
  {
    ObjType numObj = numerator;
    ObjType denObj = denominator;

    // In case one object is TH1 and the other is TGraph
    bool numIsGraph = isGraphType(numObj);
    bool denIsGraph = isGraphType(denObj);

    TH1* numH = asHistPtr(numObj);
    TH1* denH = asHistPtr(denObj);

    if ((numH && denIsGraph) || (denH && numIsGraph)) {
      if (numH) {
        numObj = histToGraph(numH);
        numIsGraph = true;
      }
      if (denH) {
        denObj = histToGraph(denH);
        denIsGraph = true;
      }
    }

    if (numH && denH) {
      TH1* comparison = static_cast<TH1*>(numH->Clone(Form("%s_%s_%d", numH->GetName(), (comparisonMode == ComparisonMode::Difference) ? "difference" : "ratio", index)));
      comparison->SetDirectory(nullptr);
      if (comparisonMode == ComparisonMode::Ratio)
        comparison->Divide(denH);
      if (comparisonMode == ComparisonMode::Difference)
        comparison->Add(denH, -1);
      if (comparisonMode == ComparisonMode::RelativeDeviation) {
        for (int i = 1; i <= comparison->GetNbinsX(); i++) {
          double A = numH->GetBinContent(i);
          double B = denH->GetBinContent(i);
          double sA = numH->GetBinError(i);
          double sB = denH->GetBinError(i);

          if (A == 0) {
            comparison->SetBinContent(i, 0);
            comparison->SetBinError(i, 0);
            continue;
          }

          double D = std::abs(A - B);
          double R = D / A;

          // Error propagation (ignore sign)
          double dRdA = D / (A * A);
          double dRdB = 1.0 / A;

          double sR = std::sqrt(dRdA * dRdA * sA * sA +
                                dRdB * dRdB * sB * sB);

          comparison->SetBinContent(i, R);
          comparison->SetBinError(i, sR);
        }
      }
      return comparison;
    }

    if (numIsGraph && denIsGraph) {
      return std::visit(
        [index, comparisonMode](auto* numPtr, auto* denPtr) -> ObjType {
          using NumT = std::decay_t<decltype(numPtr)>;
          using DenT = std::decay_t<decltype(denPtr)>;
          constexpr bool numGraph = std::is_same_v<NumT, TGraphErrors*> || std::is_same_v<NumT, TGraphAsymmErrors*>;
          constexpr bool denGraph = std::is_same_v<DenT, TGraphErrors*> || std::is_same_v<DenT, TGraphAsymmErrors*>;
          if constexpr (!(numGraph && denGraph)) {
            throw std::runtime_error("PlotTemplateSteerer: incompatible ratio object types");
          } else {
            int nNum = numPtr->GetN();
            int nDen = denPtr->GetN();

            if (nNum != nDen) {
              throw std::runtime_error("PlotTemplateSteerer: graph ratio requires equal number of points");
            }

            auto* comparison = new TGraphErrors(nNum);
            std::string cmpStr = "ratio";
            if (comparisonMode == ComparisonMode::Difference)
              cmpStr = "difference";
            if (comparisonMode == ComparisonMode::RelativeDeviation)
              cmpStr = "Rel. Dev.";
            comparison->SetTitle(Form("; %s; %s", numPtr->GetXaxis()->GetTitle(), cmpStr.c_str()));
            comparison->SetName(Form("%s_%s_%d", numPtr->GetName(), cmpStr.c_str(), index));
            for (int ip = 0; ip < nNum; ++ip) {
              double xNum = 0.0, yNum = 0.0;
              double xDen = 0.0, yDen = 0.0;
              numPtr->GetPoint(ip, xNum, yNum);
              denPtr->GetPoint(ip, xDen, yDen);

              if (std::abs(xNum - xDen) > 1e-9) {
                throw std::runtime_error("PlotTemplateSteerer: graph comparison requires matching x points");
              }

              double eNumY = 0.0;
              double eDenY = 0.0;
              double eNumX = 0.0;

              if constexpr (std::is_same_v<NumT, TGraphErrors*>) {
                eNumY = numPtr->GetErrorY(ip);
                eNumX = numPtr->GetErrorX(ip);
              } else {
                eNumY = 0.5 * (numPtr->GetErrorYlow(ip) + numPtr->GetErrorYhigh(ip));
                eNumX = 0.5 * (numPtr->GetErrorXlow(ip) + numPtr->GetErrorXhigh(ip));
              }

              if constexpr (std::is_same_v<DenT, TGraphErrors*>) {
                eDenY = denPtr->GetErrorY(ip);
              } else {
                eDenY = 0.5 * (denPtr->GetErrorYlow(ip) + denPtr->GetErrorYhigh(ip));
              }

              double compY = 0.0;
              double compEY = 0.0;
              if (comparisonMode == ComparisonMode::Ratio) {
                if (std::abs(yDen) > 0.0) {
                  compY = yNum / yDen;
                  const double relNum = (std::abs(yNum) > 0.0) ? (eNumY / yNum) : 0.0;
                  const double relDen = eDenY / yDen;
                  compEY = std::abs(compY) * std::sqrt(std::abs(relNum * relNum - relDen * relDen));
                }
              }
              if (comparisonMode == ComparisonMode::Difference) {
                compY = yNum - yDen;
                compEY = std::sqrt(eNumY * eNumY + eDenY * eDenY);
              }
              if (comparisonMode == ComparisonMode::RelativeDeviation) {
                if (std::abs(yDen) > 0.0) {
                  double D = std::abs(yNum - yDen);
                  compY = D / yDen;
                  double dDen = D / (yDen * yDen);
                  double dNum = 1.0 / yDen;
                  compEY = std::sqrt(dNum * dNum * eNumY * eNumY +
                                     dDen * dDen * eDenY * eDenY);
                }
              }
              comparison->SetPoint(ip, xNum, compY);
              comparison->SetPointError(ip, eNumX, compEY);
            }
            return comparison;
          }
        },
        numObj, denObj);
    }

    throw std::runtime_error("PlotTemplateSteerer: unsupported types for ratio objects");
  }

  std::vector<int> mFallbackColorPalette;
  std::vector<std::pair<int, int>> mPairedMarkers;
  std::vector<int> mDistinctMarkers;
  std::vector<SeriesInput> mSeriesInputs;
  mutable std::vector<ObjType> mComparisonObjects;
};

#endif // PWGCF_GENERICFRAMEWORK_CORE_PLOTTEMPLATESTEERER_H_
