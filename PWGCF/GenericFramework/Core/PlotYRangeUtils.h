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

#ifndef PWGCF_GENERICFRAMEWORK_CORE_PLOTYRANGEUTILS_H_
#define PWGCF_GENERICFRAMEWORK_CORE_PLOTYRANGEUTILS_H_

#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH1D.h"

#include <algorithm>
#include <cmath>
#include <string>
#include <type_traits>
#include <vector>

namespace plotutils
{

struct YRange {
  double ymin;
  double ymax;
  bool valid;
};

struct YPoint {
  double x;
  double y;
  double ey;
};

enum class ComparisonMode {
  None,
  PairedRatio,
  PairedDifference,
  PairedRelDev
};

template <typename ObjVariant>
std::vector<YPoint> extractPointsFromObject(const ObjVariant& obj)
{
  std::vector<YPoint> points;
  std::visit(
    [&](auto* ptr) {
      using T = std::decay_t<decltype(ptr)>;
      if (!ptr) {
        return;
      }

      if constexpr (std::is_same_v<T, TH1*> || std::is_same_v<T, TH1D*>) {
        for (int i = 1; i <= ptr->GetNbinsX(); ++i) {
          points.push_back({ptr->GetBinCenter(i), ptr->GetBinContent(i), ptr->GetBinError(i)});
        }
      } else if constexpr (std::is_same_v<T, TGraphErrors*>) {
        for (int i = 0; i < ptr->GetN(); ++i) {
          double x = 0.0, y = 0.0;
          ptr->GetPoint(i, x, y);
          points.push_back({x, y, ptr->GetErrorY(i)});
        }
      } else if constexpr (std::is_same_v<T, TGraphAsymmErrors*>) {
        for (int i = 0; i < ptr->GetN(); ++i) {
          double x = 0.0, y = 0.0;
          ptr->GetPoint(i, x, y);
          const double ey = 0.5 * (ptr->GetErrorYlow(i) + ptr->GetErrorYhigh(i));
          points.push_back({x, y, ey});
        }
      }
    },
    obj);

  return points;
}

inline void appendYAndErrorEnvelope(std::vector<double>& ys, const std::vector<YPoint>& points)
{
  for (const auto& p : points) {
    if (std::isfinite(p.y)) {
      ys.push_back(p.y);
    }
    if (std::isfinite(p.y - p.ey)) {
      ys.push_back(p.y - p.ey);
    }
    if (std::isfinite(p.y + p.ey)) {
      ys.push_back(p.y + p.ey);
    }
  }
}

template <typename SeriesCollection>
std::vector<YPoint> collectTopPadPoints(const SeriesCollection& series)
{
  std::vector<YPoint> points;
  for (const auto& s : series) {
    const auto statPts = extractPointsFromObject(s.stat);
    points.insert(points.end(), statPts.begin(), statPts.end());
    const auto sysPts = extractPointsFromObject(s.sys);
    points.insert(points.end(), sysPts.begin(), sysPts.end());
  }
  return points;
}

template <typename SeriesCollection>
std::vector<YPoint> collectComparisonPoints(const SeriesCollection& series, ComparisonMode mode)
{
  std::vector<YPoint> points;

  if (mode == ComparisonMode::None) {
    return points;
  }

  for (size_t i = 0; i + 1 < series.size(); i += 2) {
    const auto numPoints = extractPointsFromObject(series[i].stat);
    const auto denPoints = extractPointsFromObject(series[i + 1].stat);

    const size_t n = std::min(numPoints.size(), denPoints.size());
    for (size_t j = 0; j < n; ++j) {
      const auto& npt = numPoints[j];
      const auto& dpt = denPoints[j];

      if (!std::isfinite(npt.y) || !std::isfinite(dpt.y)) {
        continue;
      }

      if (mode == ComparisonMode::PairedRatio || mode == ComparisonMode::PairedRelDev) {
        if (std::abs(dpt.y) < 1e-12) {
          continue;
        }
        const double r = npt.y / dpt.y;
        double relNum = 0.0;
        double relDen = dpt.ey / std::abs(dpt.y);
        if (std::abs(npt.y) > 1e-12) {
          relNum = npt.ey / std::abs(npt.y);
        }
        const double er = std::abs(r) * std::sqrt(relNum * relNum + relDen * relDen);
        if (std::isfinite(r) && std::isfinite(er)) {
          points.push_back({npt.x, r, er});
        }
      } else {
        const double d = npt.y - dpt.y;
        const double ed = std::sqrt(npt.ey * npt.ey + dpt.ey * dpt.ey);
        if (std::isfinite(d) && std::isfinite(ed)) {
          points.push_back({npt.x, d, ed});
        }
      }
    }
  }

  return points;
}

inline std::vector<double> yValuesFromPoints(const std::vector<YPoint>& points)
{
  std::vector<double> ys;
  appendYAndErrorEnvelope(ys, points);
  return ys;
}

inline YRange robustRangeWithOutlierRejection(const std::vector<double>& values, double madCut = 3.5, double minSpan = 1e-6)
{
  if (values.empty()) {
    return {0.0, 0.0, false};
  }

  std::vector<double> sorted = values;
  std::sort(sorted.begin(), sorted.end());

  auto percentile = [&](double q) {
    if (sorted.size() == 1) {
      return sorted[0];
    }
    const double pos = q * (sorted.size() - 1);
    const int lo = static_cast<int>(std::floor(pos));
    const int hi = static_cast<int>(std::ceil(pos));
    const double t = pos - lo;
    return sorted[lo] * (1.0 - t) + sorted[hi] * t;
  };

  const double median = percentile(0.5);
  std::vector<double> absDev;
  absDev.reserve(sorted.size());
  for (double v : sorted) {
    absDev.push_back(std::abs(v - median));
  }
  std::sort(absDev.begin(), absDev.end());

  double mad = 0.0;
  if (!absDev.empty()) {
    if (absDev.size() % 2 == 0) {
      mad = 0.5 * (absDev[absDev.size() / 2 - 1] + absDev[absDev.size() / 2]);
    } else {
      mad = absDev[absDev.size() / 2];
    }
  }

  std::vector<double> filtered;
  if (mad > 0.0) {
    const double robustSigma = 1.4826 * mad;
    for (double v : sorted) {
      if (std::abs(v - median) <= madCut * robustSigma) {
        filtered.push_back(v);
      }
    }
  }

  if (filtered.size() < 3) {
    filtered = sorted;
  }

  const auto mm = std::minmax_element(filtered.begin(), filtered.end());
  double ymin = *mm.first;
  double ymax = *mm.second;

  if (!(std::isfinite(ymin) && std::isfinite(ymax))) {
    return {0.0, 0.0, false};
  }

  if (std::abs(ymax - ymin) < minSpan) {
    const double expand = std::max(1e-3, 0.2 * std::max(std::abs(ymax), 1.0));
    ymin -= expand;
    ymax += expand;
  }

  return {ymin, ymax, true};
}

inline YRange determineIdealYLimits(const std::vector<double>& ys,
                                    int padIndex,
                                    double referenceLine)
{
  std::vector<double> values = ys;
  values.push_back(referenceLine);

  const YRange base = robustRangeWithOutlierRejection(values);
  if (!base.valid) {
    if (padIndex == 1) {
      return {-0.5, 0.5, true};
    }
    if (referenceLine == 0.0) {
      return {-0.05, 0.05, true};
    }
    return {0.8, 1.2, true};
  }

  const double span = base.ymax - base.ymin;
  const double margin = (padIndex == 1) ? 0.18 : 0.22;
  double ymin = base.ymin - margin * span;
  double ymax = base.ymax + margin * span;

  if (padIndex == 2) {
    if (referenceLine == 1.0) {
      ymin = std::min(ymin, 0.95);
      ymax = std::max(ymax, 1.05);
    } else if (referenceLine == 0.0) {
      ymin = std::min(ymin, -0.005);
      ymax = std::max(ymax, 0.005);
    }
  }

  return {ymin, ymax, true};
}

template <typename SeriesCollection>
inline YRange determineTopPadYLimits(const SeriesCollection& series)
{
  return determineIdealYLimits(yValuesFromPoints(collectTopPadPoints(series)), 1, 0.0);
}

template <typename SeriesCollection>
inline YRange determineBottomPadYLimits(const SeriesCollection& series,
                                        ComparisonMode mode,
                                        double referenceLine)
{
  return determineIdealYLimits(yValuesFromPoints(collectComparisonPoints(series, mode)), 2, referenceLine);
}

} // namespace plotutils

#endif // PWGCF_GENERICFRAMEWORK_CORE_PLOTYRANGEUTILS_H_
