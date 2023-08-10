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
    const AxisSpec axisQvec{2000, -1.0, 1.0};
    const AxisSpec axisConst{12, 0., 12.}; // 4 constants x 3 detectors.

    histosQA.add("histCentFull", "Centrality distribution for valid events",
                 HistType::kTH1F, {axisCent});
    histosQA.add("Centrality_0-5/histCent", "Centrality distribution",
                 HistType::kTH1F, {axisCent});

    histosQA.add("Centrality_0-5/histQvecFT0A",
                 ("(Qx,Qy) for " + (std::string)qV::detNames[0] + " before " + (std::string)cfgCorrStep).c_str(),
                 HistType::kTH2D, {axisQvec, axisQvec});
    histosQA.add("Centrality_0-5/histQvecFT0C",
                 ("(Qx,Qy) for " + (std::string)qV::detNames[1] + " before " + (std::string)cfgCorrStep).c_str(),
                 HistType::kTH2D, {axisQvec, axisQvec});
    histosQA.add("Centrality_0-5/histQvecFV0A",
                 ("(Qx,Qy) for " + (std::string)qV::detNames[2] + " before " + (std::string)cfgCorrStep).c_str(),
                 HistType::kTH2D, {axisQvec, axisQvec});

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
    }
    if (TMath::Abs(vec.qvecFV0ARe()) < 100 && TMath::Abs(vec.qvecFV0AIm()) < 100) {
      histosQA.fill(HIST(qV::centClasses[bin]) + HIST("histQvecFV0A"),
                    vec.qvecFV0ARe(), vec.qvecFV0AIm());
    }
    LOGF(info, "QA has been filled.");
  }

  void process(aod::Qvectors const& qVecs)
  {
    // Iterate over the table and fill the QA with the received qVecs.
    for (auto& qVec : qVecs) {
      // Get the centrality bin, and fill the centrality QA histograms.
      int centBin = helperEP.GetCentBin(qVec.cent());
      LOGF(info, "Centrality percentile = %.1f Centrality bin: %d", qVec.cent(), centBin);
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

    } // Go to the next qVec.

    // The QA histograms can now be used to obtain the selected correction
    // constants for each centrality class and detector, before filling the
    // corresponding TH1 histograms.
    static_for<0, 7>([&](auto iCent) {
      constexpr int indexCent = iCent.value;
      printf("Centrality class: %s\n", qV::centClasses[indexCent].data());

      const std::shared_ptr<TH1> hist1D =
        histosQA.get<TH1>(HIST(qV::centClasses[indexCent]) + HIST("histQvecCorrConst"));

      // corr1-4 represent the 4 correction constants to obtain in one pass.
      // If "Recenter": meanX, meanY, stdX, stdY (in order).
      // If "Twist": aPlus, aMinus, lambdaPlus, lambdaMinus (in order).
      float corrConst[12] = {0.};

      const std::shared_ptr<TH2> histFT0A =
        histosQA.get<TH2>(HIST(qV::centClasses[indexCent]) + HIST("histQvecFT0A"));
      const std::shared_ptr<TH2> histFT0C =
        histosQA.get<TH2>(HIST(qV::centClasses[indexCent]) + HIST("histQvecFT0C"));
      const std::shared_ptr<TH2> histFV0A =
        histosQA.get<TH2>(HIST(qV::centClasses[indexCent]) + HIST("histQvecFV0A"));

      if ((std::string)cfgCorrStep == "Recenter") { // Get the constants for the recentering.
        std::string corrLabel[12] = {               // Label for each element of corrConst.
                                     "FT0A-meanX", "FT0A-meanY", "FT0A-stdX", "FT0A-stdY",
                                     "FT0C-meanX", "FT0C-meanY", "FT0C-stdX", "FT0C-stdY",
                                     "FV0A-meanX", "FV0A-meanY", "FV0A-stdX", "FV0A-stdY"};

        // Get the constants for FT0A.
        helperEP.GetCorrRecentering(histFT0A, corrConst[0], corrConst[1]);
        helperEP.GetCorrWidth(histFT0A, corrConst[2], corrConst[3]);

        // Get the constants for FT0C.
        helperEP.GetCorrRecentering(histFT0C, corrConst[4], corrConst[5]);
        helperEP.GetCorrWidth(histFT0C, corrConst[6], corrConst[7]);

        // Get the constants for FV0A.
        helperEP.GetCorrRecentering(histFV0A, corrConst[8], corrConst[9]);
        helperEP.GetCorrWidth(histFV0A, corrConst[10], corrConst[11]);

        for (int i = 0; i < 12; i++) {
          hist1D->GetXaxis()->SetBinLabel(i + 1, corrLabel[i].data());
        }
      } else if ((std::string)cfgCorrStep == "Twist") { // End recentering.
        std::string corrLabel[12] = {                   // Label for each element of corrConst.
                                     "FT0A-aPlus", "FT0A-aMinus", "FT0A-lambdaPlus", "FT0A-lambdaMinus",
                                     "FT0C-aPlus", "FT0C-aMinus", "FT0C-lambdaPlus", "FT0C-lambdaMinus",
                                     "FV0A-aPlus", "FV0A-aMinus", "FV0A-lambdaPlus", "FV0A-lambdaMinus"};

        // Get the constants for FT0A.
        helperEP.GetCorrTwistRecale(histFT0A, corrConst[0], corrConst[1],
                                    corrConst[2], corrConst[3]);

        // Get the constants for FT0C.
        helperEP.GetCorrTwistRecale(histFT0C, corrConst[4], corrConst[5],
                                    corrConst[6], corrConst[7]);

        // Get the constants for FV0A.
        /// NOTE: Decomment once FV0A has been implemented.
        // helperEP.GetCorrTwistRecale(histFV0A, corrConst[8], corrConst[9],
        //                             corrConst[10], corrConst[11]);

        for (int i = 0; i < 12; i++) {
          hist1D->GetXaxis()->SetBinLabel(i + 1, corrLabel[i].data());
        }
      } // End twisting+rescaling.

      for (int i = 0; i < 12; i++) {
        hist1D->SetBinContent(i + 1, corrConst[i]);
        printf("Index: %d corrConst: %e\n", i, corrConst[i]);
      }
    }); // Go to the next centrality class.
  }     // End void process(...)
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<qVectorsCorrection>(cfgc)};
}
