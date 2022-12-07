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
///
/// \brief This task demonstrates how to use configurable to wrap classes.
///        Use it with supplied configuration: "configurableObject.json".
/// \author
/// \since

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "configurableCut.h"

#include <sstream>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

template <typename T>
auto printArray(std::vector<T> const& vec)
{
  std::stringstream ss;
  ss << "[" << vec[0];
  for (auto i = 1u; i < vec.size(); ++i) {
    ss << ", " << vec[i];
  }
  ss << "]";
  return ss.str();
}

template <typename T>
auto printMatrix(Array2D<T> const& m)
{
  std::stringstream ss;
  ss << "[[" << m(0, 0);
  for (auto i = 1u; i < m.cols; ++i) {
    ss << "," << m(0, i);
  }
  for (auto i = 1u; i < m.rows; ++i) {
    ss << "], [" << m(i, 0);
    for (auto j = 1u; j < m.cols; ++j) {
      ss << "," << m(i, j);
    }
  }
  ss << "]]";
  return ss.str();
}

static constexpr float defaultm[3][4] = {{1.1, 1.2, 1.3, 1.4}, {2.1, 2.2, 2.3, 2.4}, {3.1, 3.2, 3.3, 3.4}};
static LabeledArray<float> la{&defaultm[0][0], 3, 4, {"r 1", "r 2", "r 3"}, {"c 1", "c 2", "c 3", "c 4"}};

struct ConfigurableObjectDemo {
  // Simple type configurables
  Configurable<float> min_pt{"min_pt", 0.5f, "Lower p_T cut"};
  Configurable<bool> require_tof{"require_tof", false, "Check for TOF hit"};

  // Configurable based on a struct
  Configurable<configurableCut> cut{"cut", {0.5, 1, true}, "generic cut"};
  MutableConfigurable<configurableCut> mutable_cut{"mutable_cut", {1., 2, false}, "generic cut"};

  // Array type configurables
  // note that size is fixed by this declaration - externally supplied vector needs to be the same size!
  Configurable<std::vector<int>> array{"array", {0, 0, 0, 0, 0, 0, 0}, "generic array"};
  Configurable<Array2D<float>> vmatrix{"matrix", {&defaultm[0][0], 3, 4}, "generic matrix"};
  Configurable<LabeledArray<float>> vla{"vla", {defaultm[0], 3, 4, {"r 1", "r 2", "r 3"}, {"c 1", "c 2", "c 3", "c 4"}}, "labeled array"};

  // Configurables can be grouped into `ConfigurableGroup`s.
  // Their names must be unique.
  struct : ConfigurableGroup {
    Configurable<float> max_eta{"max_eta", 0.8f, "Maximal eta"};
    Configurable<int16_t> min_clusters{"min_clusters", 70, "Minimal required number of clusters"};
  } trackcuts;

  void init(InitContext const&)
  {
    LOGP(info, "min_pt: {}; require_tof: {}", min_pt, require_tof);
    LOGF(info, "max_eta: %f; min_clusters: %d", (float)trackcuts.max_eta, (int)trackcuts.min_clusters);
    LOGF(info, "Cut1 bins: %s; Cut2 bins: %s", printArray(cut->getBins()), printArray(mutable_cut->getBins()));
    LOGF(info, "Cut1 labels: %s; Cut2 labels: %s", printArray(cut->getLabels()), printArray(mutable_cut->getLabels()));
    auto vec = (std::vector<int>)array;
    LOGF(info, "Array: %s", printArray(vec).c_str());
    LOGF(info, "Matrix: %s", printMatrix((Array2D<float>)vmatrix));
    LOGF(info, "Labeled:\n %s\n %s\n %s", printArray(vla->getLabelsRows()), printArray(vla->getLabelsCols()), printMatrix(vla->getData()));
  };

  void process(aod::Collision const&, aod::Tracks const& tracks)
  {
    std::stringstream tmpcut, tmpmutable_cut;
    tmpcut << cut;
    tmpmutable_cut << mutable_cut;
    LOGF(info, "Cut1: %s; Cut2: %s", tmpcut.str(), tmpmutable_cut.str());

    for (auto const& track : tracks) {
      if (track.globalIndex() % 500 == 0) {
        std::string decision1;
        std::string decision2;
        if (cut->method(std::abs(track.eta()))) {
          decision1 = "true";
        } else {
          decision1 = "false";
        }
        if (mutable_cut->method(std::abs(track.eta()))) {
          decision2 = "true";
        } else {
          decision2 = "false";
        }
        LOGF(info, "Cut1: %s; Cut2: %s", decision1, decision2);
        if (decision2 == "false") {
          mutable_cut->setState(-1);
        } else {
          mutable_cut->setState(1);
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<ConfigurableObjectDemo>(cfgc),
  };
}
