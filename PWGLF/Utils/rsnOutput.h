// Copyright 2019-2024 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
///
/// \author Veronika Barbasova (veronika.barbasova@cern.ch)
/// \since April 3, 2024

#ifndef PWGLF_UTILS_RSNOUTPUT_H_
#define PWGLF_UTILS_RSNOUTPUT_H_

#include <utility>
#include <string>
#include <vector>

#include "Framework/HistogramRegistry.h"
#include "Framework/Logger.h"

using namespace o2::framework;

namespace o2::analysis
{
namespace rsn
{
enum class EventType {
  zvertex,
  all
};

enum class TrackType {
  px,
  py,
  pz,
  all
};

enum class PairType {
  unlikepm,
  unlikemp,
  likepp,
  likemm,
  unliketrue,
  unlikegen,
  mixingpm,
  mixingpp,
  mixingmm,
  mixingmp,
  all
};

enum class PairAxisType {
  im,
  pt,
  mu,
  ns1,
  ns2,
  y,
  vz,
  mum,
  vzm,
  unknown
};
namespace PariAxis
{
std::vector<std::string> names{"im", "pt", "mu", "ns1", "ns2", "y", "vz", "mum", "vzm"};
}
class Output
{
 public:
  virtual ~Output() = default;

  virtual void init(std::vector<std::string> const& sparseAxes, std::vector<AxisSpec> const& allAxes, bool /*produceTrue*/ = false, bool /*eventMixing*/ = false, bool /*produceLikesign*/ = false, HistogramRegistry* registry = nullptr)
  {
    mHistogramRegistry = registry;
    if (mHistogramRegistry == nullptr)
      mHistogramRegistry = new HistogramRegistry("registry");

    // check if all axes are added in correct order
    for (int i = 0; i < static_cast<int>(PairAxisType::unknown); i++) {
      auto aname = *std::move(allAxes[i].name);
      LOGF(debug, "Check axis '%s' %d", aname.c_str(), i);
      if (aname.compare(PariAxis::names[static_cast<int>(i)])) {
        LOGF(fatal, "rsn::Output::Error: Order in allAxes is not correct !!! Expected axis '%s' and has '%s'.", aname.c_str(), PariAxis::names[static_cast<int>(i)]);
      }
    }

    PairAxisType currentType;
    for (auto& c : sparseAxes) {
      currentType = type(c);
      if (currentType >= PairAxisType::unknown) {
        LOGF(warning, "Found unknown axis (rsn::PairAxisType = %d)!!! Skipping ...", static_cast<int>(currentType));
        continue;
      }
      LOGF(info, "Adding axis '%s' to all pair histograms", c.c_str());
      mCurrentAxes.push_back(allAxes[static_cast<int>(currentType)]);
      mCurrentAxisTypes.push_back(currentType);
    }

    if (mFillPoint != nullptr)
      delete mFillPoint;
    mFillPoint = new double[mCurrentAxisTypes.size()];

    LOGF(info, "Number of axis added: %d", mCurrentAxes.size());
    mPairHisto = new HistogramConfigSpec(HistType::kTHnSparseF, mCurrentAxes);
  };

  template <typename T>
  void fillSparse(const T& h, double* point)
  {
    int i = 0;
    for (auto& at : mCurrentAxisTypes) {
      mFillPoint[i++] = point[static_cast<int>(at)];
    }
    mHistogramRegistry->get<THnSparse>(h)->Fill(mFillPoint);
  }

  virtual void fill(EventType /*t*/, double* /*point*/)
  {
    LOGF(warning, "Abstract method : 'virtual void rsn::Output::fill(EventType t, double* point)' !!! Please implement it first.");
  };
  virtual void fill(TrackType /*t*/, double* /*point*/)
  {
    LOGF(warning, "Abstract method : 'virtual void rsn::Output::fill(TrackType t, double* point)' !!! Please implement it first.");
  };
  virtual void fill(PairType /*t*/, double* /*point*/)
  {
    LOGF(warning, "Abstract method : 'virtual void rsn::Output::fill(PairType t, double* point)' !!! Please implement it first.");
  };

  virtual void fillUnlikepm(double* point) = 0;
  virtual void fillUnlikemp(double* point) = 0;
  virtual void fillLikepp(double* point) = 0;
  virtual void fillLikemm(double* point) = 0;
  virtual void fillUnliketrue(double* point) = 0;
  virtual void fillUnlikegen(double* point) = 0;
  virtual void fillMixingpm(double* point) = 0;
  virtual void fillMixingpp(double* point) = 0;
  virtual void fillMixingmm(double* point) = 0;
  virtual void fillMixingmp(double* point) = 0;

  PairAxisType type(std::string name)
  {
    auto it = std::find(PariAxis::names.begin(), PariAxis::names.end(), name);
    if (it == PariAxis::names.end()) {
      return PairAxisType::unknown;
    }
    return static_cast<PairAxisType>(std::distance(PariAxis::names.begin(), it));
  }

  std::string name(PairAxisType type)
  {
    return PariAxis::names[(static_cast<int>(type))];
  }

  AxisSpec axis(std::vector<AxisSpec> const& allAxes, PairAxisType type)
  {
    const AxisSpec unknownAxis = {1, 0., 1., "unknown axis", "unknown"};
    if (type == PairAxisType::unknown)
      return unknownAxis;
    return allAxes[static_cast<int>(type)];
  }

 protected:
  HistogramRegistry* mHistogramRegistry = nullptr;
  HistogramConfigSpec* mPairHisto = nullptr;
  std::vector<AxisSpec> mCurrentAxes;
  std::vector<PairAxisType> mCurrentAxisTypes;
  double* mFillPoint = nullptr;
};

class OutputSparse : public Output
{
 public:
  virtual void init(std::vector<std::string> const& sparseAxes, std::vector<AxisSpec> const& allAxes, bool produceTrue = false, bool eventMixing = false, bool produceLikesign = false, HistogramRegistry* registry = nullptr)
  {
    Output::init(sparseAxes, allAxes, produceTrue, eventMixing, produceLikesign, registry);
    mHistogramRegistry->add("hVz", "; V_{z} (cm); Entries", kTH1F, {{40, -20., 20.}});

    mHistogramRegistry->add("unlikepm", "Unlike pm", *mPairHisto);
    mHistogramRegistry->add("unlikemp", "Unlike mp", *mPairHisto);
    if (produceLikesign) {
      mHistogramRegistry->add("likepp", "Like PP", *mPairHisto);
      mHistogramRegistry->add("likemm", "Like MM", *mPairHisto);
    }
    if (produceTrue) {
      mHistogramRegistry->add("unliketrue", "Unlike True", *mPairHisto);
      mHistogramRegistry->add("unlikegen", "Unlike Gen", *mPairHisto);
    }
    if (eventMixing) {
      mHistogramRegistry->add("mixingpm", "Event Mixing pm", *mPairHisto);
      if (produceLikesign) {
        mHistogramRegistry->add("mixingpp", "Event Mixing pp", *mPairHisto);
        mHistogramRegistry->add("mixingmm", "Event Mixing mm", *mPairHisto);
      }
      mHistogramRegistry->add("mixingmp", "Event Mixing mp", *mPairHisto);
    }
  }

  virtual void
    fill(EventType t, double* point)
  {
    switch (t) {
      case EventType::zvertex:
        mHistogramRegistry->fill(HIST("hVz"), point[0]);
        break;
      default:
        break;
    }
  }

  virtual void fill(PairType t, double* point)
  {
    switch (t) {
      case PairType::unlikepm:
        fillUnlikepm(point);
        break;
      case PairType::unlikemp:
        fillUnlikemp(point);
        break;
      case PairType::likepp:
        fillLikepp(point);
        break;
      case PairType::likemm:
        fillLikemm(point);
        break;
      case PairType::unliketrue:
        fillUnliketrue(point);
        break;
      case PairType::unlikegen:
        fillUnlikegen(point);
        break;
      case PairType::mixingpm:
        fillMixingpm(point);
        break;
      case PairType::mixingpp:
        fillMixingpp(point);
        break;
      case PairType::mixingmm:
        fillMixingmm(point);
        break;
      case PairType::mixingmp:
        fillMixingmp(point);
        break;
      default:
        break;
    }
  }

  virtual void fillUnlikepm(double* point)
  {
    fillSparse(HIST("unlikepm"), point);
  }
  virtual void fillUnlikemp(double* point)
  {
    fillSparse(HIST("unlikemp"), point);
  }
  virtual void fillLikepp(double* point)
  {
    fillSparse(HIST("likepp"), point);
  }

  virtual void fillLikemm(double* point)
  {
    fillSparse(HIST("likemm"), point);
  }

  virtual void fillUnliketrue(double* point)
  {
    fillSparse(HIST("unliketrue"), point);
  }

  virtual void fillUnlikegen(double* point)
  {
    fillSparse(HIST("unlikegen"), point);
  }
  virtual void fillMixingpm(double* point)
  {
    fillSparse(HIST("mixingpm"), point);
  }
  virtual void fillMixingpp(double* point)
  {
    fillSparse(HIST("mixingpp"), point);
  }
  virtual void fillMixingmm(double* point)
  {
    fillSparse(HIST("mixingmm"), point);
  }
  virtual void fillMixingmp(double* point)
  {
    fillSparse(HIST("mixingmp"), point);
  }
};
} // namespace rsn
} // namespace o2::analysis

#endif // PWGLF_UTILS_RSNOUTPUT_H_
