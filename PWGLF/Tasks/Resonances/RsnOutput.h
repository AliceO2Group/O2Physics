#ifndef PWGLF_RSNOUTPUT_H_
#define PWGLF_RSNOUTPUT_H_

#include "Framework/HistogramRegistry.h"
#include "Framework/Logger.h"

using namespace o2::framework;

namespace o2::analysis
{

enum class RsnEventType {
  zvertex,
  all
};

enum class RsnTrackType {
  px,
  py,
  pz,
  all
};

enum class RsnPairType {
  unlike,
  likepp,
  likemm,
  unliketrue,
  unlikegen,
  all
};

enum class RsnPairAxisType {
  im = 0,
  pt = 1,
  mu = 2,
  ns1 = 3,
  ns2 = 4,
  y = 5,
  unknown
};
namespace RsnPariAxis
{
std::vector<std::string> names{"im", "pt", "mu", "ns1", "ns2", "y"};
}
class RsnOutput
{
 public:
  virtual ~RsnOutput() = default;

  virtual void init(std::vector<std::string> const& sparseAxes, std::vector<AxisSpec> const& allAxes, bool produceTrue = false, HistogramRegistry* registry = nullptr)
  {
    mHistogramRegistry = registry;
    if (mHistogramRegistry == nullptr)
      mHistogramRegistry = new HistogramRegistry("registry");

    // check if all axes are added in correct order
    for (int i = 0; i < static_cast<int>(RsnPairAxisType::unknown); i++) {
      auto aname = *std::move(allAxes[i].name);
      LOGF(debug, "Check axis '%s' %d", aname.c_str(), i);
      if (aname.compare(RsnPariAxis::names[static_cast<int>(i)])) {
        LOGF(fatal, "RsnOutput::Error: Order in allAxes is not correct !!! Expected axis '%s' and has '%s'.", aname.c_str(), RsnPariAxis::names[static_cast<int>(i)]);
      }
    }

    RsnPairAxisType currentType;
    for (auto& c : sparseAxes) {
      currentType = type(c);
      if (currentType >= RsnPairAxisType::unknown) {
        LOGF(warning, "Found unknown axis (RsnPairAxisType = %d)!!! Skipping ...", static_cast<int>(currentType));
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
    mPairHisto = new HistogramConfigSpec({HistType::kTHnSparseF}, mCurrentAxes);
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

  virtual void fill(RsnEventType t, double* point)
  {
    LOGF(warning, "Abstract method : 'virtual void RsnOutput::fill(RsnEventType t, double* point)' !!! Please implement it first.");
  };
  virtual void fill(RsnTrackType t, double* point)
  {
    LOGF(warning, "Abstract method : 'virtual void RsnOutput::fill(RsnTrackType t, double* point)' !!! Please implement it first.");
  };
  virtual void fill(RsnPairType t, double* point)
  {
    LOGF(warning, "Abstract method : 'virtual void RsnOutput::fill(RsnPairType t, double* point)' !!! Please implement it first.");
  };

  virtual void fillUnlike(double* point) = 0;
  virtual void fillLikepp(double* point) = 0;
  virtual void fillLikemm(double* point) = 0;
  virtual void fillUnliketrue(double* point) = 0;
  virtual void fillUnlikegen(double* point) = 0;
 
  RsnPairAxisType type(std::string name)
  {
    auto it = std::find(RsnPariAxis::names.begin(), RsnPariAxis::names.end(), name);
    if (it == RsnPariAxis::names.end()) {
      return RsnPairAxisType::unknown;
    }
    return static_cast<RsnPairAxisType>(std::distance(RsnPariAxis::names.begin(), it));
  }

  std::string name(RsnPairAxisType type)
  {
    return RsnPariAxis::names[(static_cast<int>(type))];
  }

  AxisSpec axis(std::vector<AxisSpec> const& allAxes, RsnPairAxisType type)
  {
    const AxisSpec unknownAxis = {1, 0., 1., "unknown axis", "unknown"};
    if (type == RsnPairAxisType::unknown)
      return unknownAxis;
    return allAxes[static_cast<int>(type)];
  }

 protected:
  HistogramRegistry* mHistogramRegistry = nullptr;
  HistogramConfigSpec* mPairHisto = nullptr;
  std::vector<AxisSpec> mCurrentAxes;
  std::vector<RsnPairAxisType> mCurrentAxisTypes;
  double* mFillPoint = nullptr;
};

class RsnOutputSparse : public RsnOutput
{
 public:
  virtual void init(std::vector<std::string> const& sparseAxes, std::vector<AxisSpec> const& allAxes, bool produceTrue = false, HistogramRegistry* registry = nullptr)
  {
    RsnOutput::init(sparseAxes, allAxes, produceTrue, registry);
    mHistogramRegistry->add("hVz", "; vtx_{z} (cm); Entries", kTH1F, {{40, -20., 20.}});

    mHistogramRegistry->add("unlike", "Unlike", *mPairHisto);

    mHistogramRegistry->add("likepp", "Like PP", *mPairHisto);
    mHistogramRegistry->add("likemm", "Like MM", *mPairHisto);
    if (produceTrue) {
      mHistogramRegistry->add("unlikeTrue", "Unlike True", *mPairHisto);
      mHistogramRegistry->add("unlikeGen", "Unlike Gen", *mPairHisto);
    }
  }

  virtual void
    fill(RsnEventType t, double* point)
  {
    switch (t) {
      case RsnEventType::zvertex:
        mHistogramRegistry->fill(HIST("hVz"), point[0]);
        break;
      default:
        break;
    }
  }

  virtual void fill(RsnPairType t, double* point)
  {
    switch (t) {
      case RsnPairType::unlike:
        fillUnlike(point);
        break;
      case RsnPairType::likepp:
        fillLikepp(point);
        break;
      case RsnPairType::likemm:
        fillLikemm(point);
        break;
      case RsnPairType::unliketrue:
        fillUnliketrue(point);
        break;
      case RsnPairType::unlikegen:
        fillUnlikegen(point);
        break;
      default:
        break;
    }
  }

  virtual void fillUnlike(double* point)
  {
    fillSparse(HIST("unlike"), point);
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
};
} // namespace o2::analysis

#endif