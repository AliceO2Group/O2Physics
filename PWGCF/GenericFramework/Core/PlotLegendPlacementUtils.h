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

#ifndef PWGCF_GENERICFRAMEWORK_CORE_PLOTLEGENDPLACEMENTUTILS_H_
#define PWGCF_GENERICFRAMEWORK_CORE_PLOTLEGENDPLACEMENTUTILS_H_

#include "PlotYRangeUtils.h"

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

namespace plotutils
{

struct NDCBox {
  double x1;
  double y1;
  double x2;
  double y2;
};

inline bool insideBox(const NDCBox& b, double x, double y)
{
  return (x >= b.x1 && x <= b.x2 && y >= b.y1 && y <= b.y2);
}

inline bool overlapBox(const NDCBox& a, const NDCBox& b)
{
  return !(a.x2 < b.x1 || a.x1 > b.x2 || a.y2 < b.y1 || a.y1 > b.y2);
}

inline double mapXToNDC(double x, double xMin, double xMax, double left, double right)
{
  if (xMax <= xMin) {
    return 0.5;
  }
  return left + (x - xMin) / (xMax - xMin) * (1.0 - left - right);
}

inline double mapYToNDC(double y, double yMin, double yMax, double bottom, double top)
{
  if (yMax <= yMin) {
    return 0.5;
  }
  return bottom + (y - yMin) / (yMax - yMin) * (1.0 - bottom - top);
}

inline NDCBox chooseLegendBox(const std::vector<YPoint>& points,
                              double xMin,
                              double xMax,
                              double yMin,
                              double yMax,
                              double left,
                              double right,
                              double bottom,
                              double top,
                              const std::vector<std::string>& legendLabels,
                              bool avoidAppliedText,
                              const std::string& centralityLabel)
{
  const double innerW = std::max(0.05, 1.0 - left - right);
  const double innerH = std::max(0.05, 1.0 - bottom - top);

  size_t maxChars = 0;
  for (const auto& l : legendLabels) {
    maxChars = std::max(maxChars, l.size());
  }
  if (maxChars == 0) {
    maxChars = 6;
  }

  const int nEntries = std::max(1, static_cast<int>(legendLabels.size()));

  const double padGapX = 0.02;
  const double padGapY = 0.03;

  const double w = std::min(0.42 * innerW, std::max(0.14, 0.04 + 0.009 * static_cast<double>(maxChars)));
  const double h = std::min(0.42 * innerH, std::max(0.11, 0.055 + 0.048 * static_cast<double>(nEntries)));

  const double availXMin = left + padGapX;
  const double availXMax = 1.0 - right - padGapX;
  const double availYMin = bottom + padGapY;
  const double availYMax = 1.0 - top - padGapY;

  const double spanX = std::max(0.0, (availXMax - availXMin - w));
  const double spanY = std::max(0.0, (availYMax - availYMin - h));

  const double xLeft = availXMin;
  const double xRight = std::max(xLeft, availXMax - w);
  const double xCenter = availXMin + 0.5 * spanX;
  const double xQRightCenter = availXMin + 0.75 * spanX;
  const double xQLeftCenter = availXMin + 0.25 * spanX;

  const double yBottom = availYMin;
  const double yTop = std::max(yBottom, availYMax - h);
  const double yQTopCenter = availYMin + 0.75 * spanY;
  const double yQBottomCenter = availYMin + 0.25 * spanY;

  struct Candidate {
    NDCBox box;
    double prefPenalty;
  };

  auto mk = [&](double x1, double y1, double penalty) {
    return Candidate{{x1, y1, x1 + w, y1 + h}, penalty};
  };

  std::vector<Candidate> candidates = {
    mk(xQRightCenter, yQTopCenter, 0.0),
    mk(xRight, yTop, 6.0),
    mk(xCenter, yQTopCenter, 10.0),
    mk(xQLeftCenter, yQTopCenter, 14.0),
    mk(xLeft, yTop, 18.0),
    mk(xQRightCenter, yQBottomCenter, 22.0),
    mk(xCenter, yQBottomCenter, 26.0),
    mk(xQLeftCenter, yQBottomCenter, 30.0),
    mk(xRight, yBottom, 34.0),
    mk(xLeft, yBottom, 38.0)};

  std::vector<NDCBox> reservedTextBoxes;
  if (avoidAppliedText) {
    const double yStart = 1.0 - top - 0.1;
    const double lineStep = 0.08;
    const double xLeftText = left + 0.015;

    const double textBlockW = std::min(0.50 * innerW, std::max(0.24, 0.20 + 0.007 * 32.0));
    const double textBlockH = 3.0 * lineStep;
    reservedTextBoxes.push_back({xLeftText - 0.005, yStart - textBlockH - 0.02, xLeftText + textBlockW, yStart + 0.02});

    const double xRightLabel = 1.0 - right - 0.02;
    const double labelW = std::min(0.35 * innerW, std::max(0.08, 0.055 + 0.008 * static_cast<double>(centralityLabel.size())));
    const double labelH = 0.06;
    reservedTextBoxes.push_back({xRightLabel - labelW, yStart - labelH, xRightLabel, yStart + 0.01});
  }

  double bestScore = 1e99;
  NDCBox best = candidates.front().box;

  for (const auto& c : candidates) {
    double score = c.prefPenalty;

    for (const auto& r : reservedTextBoxes) {
      if (overlapBox(c.box, r)) {
        score += 1e6;
      }
    }

    if (c.box.x1 < availXMin || c.box.x2 > availXMax || c.box.y1 < availYMin || c.box.y2 > availYMax) {
      score += 1e6;
    }

    for (const auto& p : points) {
      if (!std::isfinite(p.x) || !std::isfinite(p.y) || !std::isfinite(p.ey)) {
        continue;
      }

      const double xNdc = mapXToNDC(p.x, xMin, xMax, left, right);
      const double yNdc = mapYToNDC(p.y, yMin, yMax, bottom, top);
      const double yLowNdc = mapYToNDC(p.y - p.ey, yMin, yMax, bottom, top);
      const double yHighNdc = mapYToNDC(p.y + p.ey, yMin, yMax, bottom, top);

      if (insideBox(c.box, xNdc, yNdc)) {
        score += 1000.0;
      }
      if (insideBox(c.box, xNdc, yLowNdc) || insideBox(c.box, xNdc, yHighNdc)) {
        score += 300.0;
      }
    }

    if (score < bestScore) {
      bestScore = score;
      best = c.box;
    }
  }

  return best;
}

} // namespace plotutils

#endif // PWGCF_GENERICFRAMEWORK_CORE_PLOTLEGENDPLACEMENTUTILS_H_
