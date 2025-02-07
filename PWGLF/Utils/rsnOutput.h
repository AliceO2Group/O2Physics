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
  unlikegenold,
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
  ce,
  ns1,
  ns2,
  eta,
  y,
  vz,
  mum,
  cem,
  vzm,
  unknown
};

enum class MixingType {
  ce,
  mu,
  none
};

MixingType mixingTypeName(std::string name)
{
  if (name == "ce")
    return MixingType::ce;
  else if (name == "mu")
    return MixingType::mu;

  return MixingType::none;
}

enum class SystematicsAxisType {
  ncl,
  unknown
};
namespace PairAxis
{
std::vector<std::string> names{"im", "pt", "mu", "ce", "ns1", "ns2", "eta", "y", "vz", "mum", "cem", "vzm"};
}

namespace SystematicsAxis
{
std::vector<std::string> names{"ncl"};
}

class Output
{
 public:
  virtual ~Output() = default;

  virtual void init(std::vector<std::string> const& sparseAxes, std::vector<AxisSpec> const& allAxes, std::vector<std::string> const& sysAxes, std::vector<AxisSpec> const& allAxes_sys, bool /*produceTrue*/ = false, MixingType /*eventMixing*/ = MixingType::none, bool /*produceLikesign*/ = false, HistogramRegistry* registry = nullptr)
  {
    mHistogramRegistry = registry;
    if (mHistogramRegistry == nullptr)
      mHistogramRegistry = new HistogramRegistry("registry");

    // check if all axes are added in correct order
    for (int i = 0; i < static_cast<int>(PairAxisType::unknown); i++) {
      auto aname = *std::move(allAxes[i].name);
      LOGF(debug, "Check axis '%s' %d", aname.c_str(), i);
      if (aname.compare(PairAxis::names[static_cast<int>(i)])) {
        LOGF(fatal, "rsn::Output::Error: Order in allAxes is not correct !!! Expected axis '%s' and has '%s'.", aname.c_str(), PairAxis::names[static_cast<int>(i)]);
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

    // check if all systematic axes are added in correct order
    for (int i = 0; i < static_cast<int>(SystematicsAxisType::unknown); i++) {
      auto aname = *std::move(allAxes_sys[i].name);
      LOGF(debug, "Check axis '%s' %d", aname.c_str(), i);
      if (aname.compare(SystematicsAxis::names[static_cast<int>(i)])) {
        LOGF(fatal, "rsn::Output::Error: Order in allAxes_sys is not correct !!! Expected axis '%s' and has '%s'.", aname.c_str(), SystematicsAxis::names[static_cast<int>(i)]);
      }
    }

    SystematicsAxisType currentType_sys;
    for (auto& c : sysAxes) {
      currentType_sys = type_sys(c);
      if (currentType_sys >= SystematicsAxisType::unknown) {
        LOGF(warning, "Found unknown axis (rsn::SystematicsAxisType = %d)!!! Skipping ...", static_cast<int>(currentType_sys));
        continue;
      }
      LOGF(info, "Adding axis '%s' to systematic histogram", c.c_str());
      mCurrentAxesSys.push_back(allAxes_sys[static_cast<int>(currentType_sys)]);
      mCurrentAxisTypesSys.push_back(currentType_sys);
    }

    if (mFillPointSys != nullptr)
      delete mFillPointSys;
    mFillPointSys = new double[mCurrentAxisTypesSys.size()];

    LOGF(info, "Number of systematic axis added: %d", mCurrentAxesSys.size());
    mPairHistoSys = new HistogramConfigSpec(HistType::kTHnSparseF, mCurrentAxesSys);
  }

  template <typename T>
  void fillSparse(const T& h, double* point)
  {
    int i = 0;
    for (auto& at : mCurrentAxisTypes) {
      mFillPoint[i++] = point[static_cast<int>(at)];
    }
    mHistogramRegistry->get<THnSparse>(h)->Fill(mFillPoint);
  }

  template <typename T>
  void fillSparseSys(const T& h, double* point)
  {
    int i = 0;
    for (auto& at : mCurrentAxisTypesSys) {
      mFillPointSys[i++] = point[static_cast<int>(at)];
    }
    mHistogramRegistry->get<THnSparse>(h)->Fill(mFillPointSys);
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
  virtual void fillUnlikegenOld(double* point) = 0;
  virtual void fillMixingpm(double* point) = 0;
  virtual void fillMixingpp(double* point) = 0;
  virtual void fillMixingmm(double* point) = 0;
  virtual void fillMixingmp(double* point) = 0;
  virtual void fillSystematics(double* point) = 0;

  PairAxisType type(std::string name)
  {
    auto it = std::find(PairAxis::names.begin(), PairAxis::names.end(), name);
    if (it == PairAxis::names.end()) {
      return PairAxisType::unknown;
    }
    return static_cast<PairAxisType>(std::distance(PairAxis::names.begin(), it));
  }

  SystematicsAxisType type_sys(std::string name)
  {
    auto it = std::find(SystematicsAxis::names.begin(), SystematicsAxis::names.end(), name);
    if (it == SystematicsAxis::names.end()) {
      return SystematicsAxisType::unknown;
    }
    return static_cast<SystematicsAxisType>(std::distance(SystematicsAxis::names.begin(), it));
  }

  std::string name(PairAxisType type)
  {
    return PairAxis::names[(static_cast<int>(type))];
  }

  std::string name_sys(SystematicsAxisType type)
  {
    return SystematicsAxis::names[(static_cast<int>(type))];
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
  HistogramConfigSpec* mPairHistoSys = nullptr;
  std::vector<AxisSpec> mCurrentAxes;
  std::vector<PairAxisType> mCurrentAxisTypes;
  std::vector<AxisSpec> mCurrentAxesSys;
  std::vector<SystematicsAxisType> mCurrentAxisTypesSys;
  double* mFillPoint = nullptr;
  double* mFillPointSys = nullptr;
};

class OutputSparse : public Output
{
 public:
  virtual void init(std::vector<std::string> const& sparseAxes, std::vector<AxisSpec> const& allAxes, std::vector<std::string> const& sysAxes, std::vector<AxisSpec> const& allAxes_sys, bool produceTrue = false, MixingType eventMixing = MixingType::none, bool produceLikesign = false, HistogramRegistry* registry = nullptr)
  {
    Output::init(sparseAxes, allAxes, sysAxes, allAxes_sys, produceTrue, eventMixing, produceLikesign, registry);

    mHistogramRegistry->add("unlikepm", "Unlike pm", *mPairHisto);
    // mHistogramRegistry->add("unlikemp", "Unlike mp", *mPairHisto);
    if (produceLikesign) {
      mHistogramRegistry->add("likepp", "Like PP", *mPairHisto);
      mHistogramRegistry->add("likemm", "Like MM", *mPairHisto);
    }
    if (produceTrue) {
      mHistogramRegistry->add("unliketrue", "Unlike True", *mPairHisto);
      mHistogramRegistry->add("unlikegen", "Unlike Gen", *mPairHisto);
      mHistogramRegistry->add("unlikegenold", "Unlike Gen Old", *mPairHisto);
    }
    if (eventMixing != MixingType::none) {
      mHistogramRegistry->add("mixingpm", "Event Mixing pm", *mPairHisto);
      if (produceLikesign) {
        mHistogramRegistry->add("mixingpp", "Event Mixing pp", *mPairHisto);
        mHistogramRegistry->add("mixingmm", "Event Mixing mm", *mPairHisto);
      }
      mHistogramRegistry->add("mixingmp", "Event Mixing mp", *mPairHisto);
    }

    mHistogramRegistry->add("Mapping/systematics", "Systematics mapping", *mPairHistoSys);
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
      case PairType::unlikegenold:
        fillUnlikegenOld(point);
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
  virtual void fillUnlikegenOld(double* point)
  {
    fillSparse(HIST("unlikegenold"), point);
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
  virtual void fillSystematics(double* point)
  {
    fillSparse(HIST("Mapping/systematics"), point);
  }
};
} // namespace rsn
} // namespace o2::analysis

#endif // PWGLF_UTILS_RSNOUTPUT_H_
