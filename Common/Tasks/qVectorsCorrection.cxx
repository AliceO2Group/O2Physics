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
/// \file   qVectorsCorrection.cxx
/// \author Cindy Mordasini <cindy.mordasini@cern.ch>
/// \author Anna Önnerstad <anna.onnerstad@cern.ch>
///
/// \brief  ...
///

// C++/ROOT includes.
#include <chrono>
#include <string>
#include <vector>
#include <TComplex.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TMath.h>

// o2Physics includes.
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/StaticFor.h"

#include "Common/DataModel/Qvectors.h"
#include "Common/Core/EventPlaneHelper.h"

// o2 includes.

using namespace o2;
using namespace o2::framework;

namespace qV
{
static constexpr std::string_view centClasses[] = {
  "Centrality_0-5/", "Centrality_5-10/", "Centrality_10-20/", "Centrality_20-30/",
  "Centrality_30-40/", "Centrality_40-50/", "Centrality_50-60/", "Centrality_60-80/"};

static constexpr std::string_view detNames[] = {"FT0A", "FT0C", "FV0A"};
} // namespace qV

struct qVectorsCorrection {
  // Configurables.
  Configurable<std::string> cfgCentEsti{"cfgCentEsti", // List from qVectorsTable.cxx
                                        "FT0C", "Centrality estimator (Run3): 0 = FT0M, 1 = FT0A, 2 = FT0C, 3 = FV0A"};
  Configurable<std::string> cfgCorrStep{"cfgCorrStep", // Used in the plotting.
                                        "Recenter", "Correction step to obtain: Recenter, Twist, Rescale, Final"};
  // No correction = recenter, recentered Qvectors = twist, twisted Qvectors = rescale.
  // NOTE: As of no, the twist gets both twist and rescale correction constants.

  // Histogram registry for the output QA figures and list of centrality classes for it.
  // Objects are NOT saved in alphabetical orders, and registry names are NOT saved
  // as TDirectoryFile.
  HistogramRegistry histosQA{"histosQA", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};

  // Helper variables.
  EventPlaneHelper helperEP;

  void init(InitContext const&)
  {
    // Fill the registry with the needed objects.
    const AxisSpec axisCent{100, 0., 100., fmt::format("Centrality percentile ({})", (std::string)cfgCentEsti)};
    const AxisSpec axisQvec{2000, -5, 5};
    const AxisSpec axisConst{12, 0., 12.}; // 4 constants x 3 detectors.

    histosQA.add("histCentFull", "Centrality distribution for valid events",
                 HistType::kTH1F, {axisCent});
    histosQA.add("Centrality_0-5/histCent", "Centrality distribution",
                 HistType::kTH1F, {axisCent});

    histosQA.add("Centrality_0-5/histQvecFT0A",
                 ("(Qx,Qy) for " + (std::string)qV::detNames[0] + " with corrections ").c_str(),
                 HistType::kTH2F, {axisQvec, axisQvec});
    histosQA.add("Centrality_0-5/histQvecFT0C",
                 ("(Qx,Qy) for " + (std::string)qV::detNames[1] + " with corrections ").c_str(),
                 HistType::kTH2F, {axisQvec, axisQvec});
    histosQA.add("Centrality_0-5/histQvecFV0A",
                 ("(Qx,Qy) for " + (std::string)qV::detNames[2] + " with corrections ").c_str(),
                 HistType::kTH2F, {axisQvec, axisQvec});

    histosQA.add("Centrality_0-5/histQvecFT0C_uncor", "", {HistType::kTH2F, {axisQvec, axisQvec}});
    histosQA.add("Centrality_0-5/histQvecFT0C_rectr", "", {HistType::kTH2F, {axisQvec, axisQvec}});
    histosQA.add("Centrality_0-5/histQvecFT0C_twist", "", {HistType::kTH2F, {axisQvec, axisQvec}});

    histosQA.add("Centrality_0-5/histQvecCorrConst",
                 ("Correction constants for " + (std::string)cfgCorrStep).c_str(),
                 HistType::kTH1F, {axisConst});

    for (int iBin = 1; iBin < 8; iBin++) {
      histosQA.addClone("Centrality_0-5/", qV::centClasses[iBin].data());
    }

  } // End void init(InitContext const&)

  // Definition of all the needed template functions.
  template <int bin, typename T>
  void fillHistosQA(const T& vec)
  {
    // Fill the centrality distribution per class for the given bin.
    histosQA.fill(HIST(qV::centClasses[bin]) + HIST("histCent"), vec.cent());

    // Fill the (Qx,Qy) distributions for each detector, after removing dummy values.
    /// NOTE: FV0 (and FT0C?) are not fully implemented yet
    /// --> Values are just dummy placeholders.
    if (TMath::Abs(vec.qvecFT0ARe()) < 100 && TMath::Abs(vec.qvecFT0AIm()) < 100) {
      histosQA.fill(HIST(qV::centClasses[bin]) + HIST("histQvecFT0A"),
                    vec.qvecFT0ARe(), vec.qvecFT0AIm());
    }
    if (TMath::Abs(vec.qvecFT0CRe()) < 100 && TMath::Abs(vec.qvecFT0CIm()) < 100) {
      histosQA.fill(HIST(qV::centClasses[bin]) + HIST("histQvecFT0C"),
                    vec.qvecFT0CRe(), vec.qvecFT0CIm());
      histosQA.fill(HIST(qV::centClasses[bin]) + HIST("histQvecFT0C_uncor"), vec.qvecFT0CUncorRe(), vec.qvecFT0CUncorIm());
      histosQA.fill(HIST(qV::centClasses[bin]) + HIST("histQvecFT0C_rectr"), vec.qvecFT0CRectrRe(), vec.qvecFT0CRectrIm());
      histosQA.fill(HIST(qV::centClasses[bin]) + HIST("histQvecFT0C_twist"), vec.qvecFT0CTwistRe(), vec.qvecFT0CTwistIm());
    }
    if (TMath::Abs(vec.qvecFV0ARe()) < 100 && TMath::Abs(vec.qvecFV0AIm()) < 100) {
      histosQA.fill(HIST(qV::centClasses[bin]) + HIST("histQvecFV0A"),
                    vec.qvecFV0ARe(), vec.qvecFV0AIm());
    }
  }

  void process(aod::Qvectors const& qVecs)
  {
    // Iterate over the table and fill the QA with the received qVecs.
    for (auto& qVec : qVecs) {
      // Get the centrality bin, and fill the centrality QA histograms.
      int centBin = helperEP.GetCentBin(qVec.cent());
      histosQA.fill(HIST("histCentFull"), qVec.cent());
      if (centBin < 0 || centBin > 8) {
        continue;
      }
      switch (centBin) { // LOKI: See if we can replace that with a const something like below.
        case 0:
          fillHistosQA<0>(qVec);
          break;
        case 1:
          fillHistosQA<1>(qVec);
          break;
        case 2:
          fillHistosQA<2>(qVec);
          break;
        case 3:
          fillHistosQA<3>(qVec);
          break;
        case 4:
          fillHistosQA<4>(qVec);
          break;
        case 5:
          fillHistosQA<5>(qVec);
          break;
        case 6:
          fillHistosQA<6>(qVec);
          break;
        case 7:
          fillHistosQA<7>(qVec);
          break;
      } // End switch(centBin)

    }   // Go to the next qVec.
  }     // End void process(...)
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<qVectorsCorrection>(cfgc)};
}
