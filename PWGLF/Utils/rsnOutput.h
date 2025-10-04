// Copyright 2019-2025 CERN and copyright holders of ALICE O2.
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
/// \file rsnOutput.h
/// \brief Resonance output class.
/// \author Veronika Barbasova (veronika.barbasova@cern.ch)

#ifndef PWGLF_UTILS_RSNOUTPUT_H_
#define PWGLF_UTILS_RSNOUTPUT_H_

#include "Framework/HistogramRegistry.h"
#include "Framework/Logger.h"

#include <string>
#include <utility>
#include <vector>

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
  mixingmp,
  rotationpm,
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
namespace pair_axis
{
std::vector<std::string> names{"im", "pt", "mu", "ce", "ns1", "ns2", "eta", "y", "vz", "mum", "cem", "vzm"};
}

namespace systematic_axis
{
std::vector<std::string> names{"ncl"};
}

class Output
{
 public:
  virtual ~Output() = default;

  virtual void init(std::vector<std::string> const& sparseAxes, std::vector<o2::framework::AxisSpec> const& allAxes, std::vector<std::string> const& sysAxes, std::vector<o2::framework::AxisSpec> const& allAxes_sys, bool /*produceTrue*/ = false, MixingType /*eventMixing*/ = MixingType::none, bool /*produceLikesign*/ = false, bool /*produceRotational*/ = false, o2::framework::HistogramRegistry* registry = nullptr)
  {
    mHistogramRegistry = registry;
    if (mHistogramRegistry == nullptr)
      mHistogramRegistry = new o2::framework::HistogramRegistry("registry");

    // check if all axes are added in correct order
    for (int i = 0; i < static_cast<int>(PairAxisType::unknown); i++) {
      auto aname = *std::move(allAxes[i].name);
      LOGF(debug, "Check axis '%s' %d", aname.c_str(), i);
      if (aname.compare(pair_axis::names[static_cast<int>(i)])) {
        LOGF(fatal, "rsn::Output::Error: Order in allAxes is not correct !!! Expected axis '%s' and has '%s'.", aname.c_str(), pair_axis::names[static_cast<int>(i)]);
      }
    }

    PairAxisType currentType;
    for (const auto& c : sparseAxes) {
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
    mPairHisto = new o2::framework::HistogramConfigSpec(o2::framework::HistType::kTHnSparseF, mCurrentAxes);

    // check if all systematic axes are added in correct order
    for (int i = 0; i < static_cast<int>(SystematicsAxisType::unknown); i++) {
      auto aname = *std::move(allAxes_sys[i].name);
      LOGF(debug, "Check axis '%s' %d", aname.c_str(), i);
      if (aname.compare(systematic_axis::names[static_cast<int>(i)])) {
        LOGF(fatal, "rsn::Output::Error: Order in allAxes_sys is not correct !!! Expected axis '%s' and has '%s'.", aname.c_str(), systematic_axis::names[static_cast<int>(i)]);
      }
    }

    SystematicsAxisType currentTypeSys;
    for (const auto& c : sysAxes) {
      currentTypeSys = typeSys(c);
      if (currentTypeSys >= SystematicsAxisType::unknown) {
        LOGF(warning, "Found unknown axis (rsn::SystematicsAxisType = %d)!!! Skipping ...", static_cast<int>(currentTypeSys));
        continue;
      }
      LOGF(info, "Adding axis '%s' to systematic histogram", c.c_str());
      mCurrentAxesSys.push_back(allAxes_sys[static_cast<int>(currentTypeSys)]);
      mCurrentAxisTypesSys.push_back(currentTypeSys);
    }

    if (mFillPointSys != nullptr)
      delete mFillPointSys;
    mFillPointSys = new double[mCurrentAxisTypesSys.size()];

    LOGF(info, "Number of systematic axis added: %d", mCurrentAxesSys.size());
    mPairHistoSys = new o2::framework::HistogramConfigSpec(o2::framework::HistType::kTHnSparseF, mCurrentAxesSys);
  }

  template <typename T>
  void fillSparse(const T& h, double* point)
  {
    int i = 0;
    for (const auto& at : mCurrentAxisTypes) {
      mFillPoint[i++] = point[static_cast<int>(at)];
    }
    mHistogramRegistry->get<THnSparse>(h)->Fill(mFillPoint);
  }

  template <typename T>
  void fillSparseSys(const T& h, double* point)
  {
    int i = 0;
    for (const auto& at : mCurrentAxisTypesSys) {
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
  virtual void fillMixingmp(double* point) = 0;
  virtual void fillRotationpm(double* point) = 0;
  virtual void fillSystematics(double* point) = 0;

  PairAxisType type(std::string name)
  {
    auto it = std::find(pair_axis::names.begin(), pair_axis::names.end(), name);
    if (it == pair_axis::names.end()) {
      return PairAxisType::unknown;
    }
    return static_cast<PairAxisType>(std::distance(pair_axis::names.begin(), it));
  }

  SystematicsAxisType typeSys(std::string name)
  {
    auto it = std::find(systematic_axis::names.begin(), systematic_axis::names.end(), name);
    if (it == systematic_axis::names.end()) {
      return SystematicsAxisType::unknown;
    }
    return static_cast<SystematicsAxisType>(std::distance(systematic_axis::names.begin(), it));
  }

  std::string name(PairAxisType type)
  {
    return pair_axis::names[(static_cast<int>(type))];
  }

  std::string nameSys(SystematicsAxisType type)
  {
    return systematic_axis::names[(static_cast<int>(type))];
  }

  o2::framework::AxisSpec axis(std::vector<o2::framework::AxisSpec> const& allAxes, PairAxisType type)
  {
    const o2::framework::AxisSpec unknownAxis = {1, 0., 1., "unknown axis", "unknown"};
    if (type == PairAxisType::unknown)
      return unknownAxis;
    return allAxes[static_cast<int>(type)];
  }

 protected:
  o2::framework::HistogramRegistry* mHistogramRegistry = nullptr;
  o2::framework::HistogramConfigSpec* mPairHisto = nullptr;
  o2::framework::HistogramConfigSpec* mPairHistoSys = nullptr;
  std::vector<o2::framework::AxisSpec> mCurrentAxes;
  std::vector<PairAxisType> mCurrentAxisTypes;
  std::vector<o2::framework::AxisSpec> mCurrentAxesSys;
  std::vector<SystematicsAxisType> mCurrentAxisTypesSys;
  double* mFillPoint = nullptr;
  double* mFillPointSys = nullptr;
};

class OutputSparse : public Output
{
 public:
  virtual void init(std::vector<std::string> const& sparseAxes, std::vector<o2::framework::AxisSpec> const& allAxes, std::vector<std::string> const& sysAxes, std::vector<o2::framework::AxisSpec> const& allAxes_sys, bool produceTrue = false, MixingType eventMixing = MixingType::none, bool produceLikesign = false, bool produceRotational = false, o2::framework::HistogramRegistry* registry = nullptr)
  {
    Output::init(sparseAxes, allAxes, sysAxes, allAxes_sys, produceTrue, eventMixing, produceLikesign, produceRotational, registry);

    mHistogramRegistry->add("unlikepm", "Unlike pm", *mPairHisto);
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
      mHistogramRegistry->add("mixingmp", "Event Mixing mp", *mPairHisto);
    }
    if (produceRotational) {
      mHistogramRegistry->add("rotationpm", "Rotational pm", *mPairHisto);
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
      case PairType::mixingmp:
        fillMixingmp(point);
        break;
      case PairType::rotationpm:
        fillRotationpm(point);
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
  virtual void fillMixingmp(double* point)
  {
    fillSparse(HIST("mixingmp"), point);
  }
  virtual void fillRotationpm(double* point)
  {
    fillSparse(HIST("rotationpm"), point);
  }
  virtual void fillSystematics(double* point)
  {
    fillSparse(HIST("Mapping/systematics"), point);
  }
};
} // namespace rsn
} // namespace o2::analysis

#endif // PWGLF_UTILS_RSNOUTPUT_H_
