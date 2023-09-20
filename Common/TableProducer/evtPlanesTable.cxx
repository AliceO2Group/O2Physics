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
/// \file   evtPlanesTable.cxx
/// \author Cindy Mordasini <cindy.mordasini@cern.ch>
/// \author Anna Önnerstad <anna.onnerstad@cern.ch>
///
/// \brief  Task calculating the Q-vectors for each collision in a bunch crossing
///         (with or without corrections) and save the results in a dedicated table.
///

// C++/ROOT includes.
#include <chrono>
#include <string>
#include <vector>
#include <TComplex.h>
#include <TMath.h>

// o2Physics includes.
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/StaticFor.h"

#include "Common/Core/EventPlaneHelper.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/FT0Corrected.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"

#include "Common/DataModel/Qvectors.h"
#include "Common/DataModel/EvtPlanes.h"

// o2 includes.
#include "CCDB/BasicCCDBManager.h"
#include "DetectorsCommonDataFormats/AlignParam.h"

using namespace o2;
using namespace o2::framework;

namespace ep
{
static constexpr std::string_view centClasses[] = {
  "Centrality_0-5/", "Centrality_5-10/", "Centrality_10-20/", "Centrality_20-30/",
  "Centrality_30-40/", "Centrality_40-50/", "Centrality_50-60/", "Centrality_60-80/"};

static constexpr std::string_view detNames[] = {
  "FT0A", "FT0C", "FV0A", "BPos", "BNeg",
  "FT0CUC", "FT0CRC", "FT0CTW"};
} // namespace ep

struct evtPlanesTable {
  // Configurables.
  Configurable<std::string> cfgCentEsti{"cfgCentEsti", // List from qVectorsTable.cxx
                                        "FT0C", "Centrality estimator (Run3): 0 = FT0M, 1 = FT0A, 2 = FT0C, 3 = FV0A"};
  Configurable<std::string> cfgCorrStep{"cfgCorrStep", // Used in the plotting.
                                        "Recentered", "Latest correction applied: Raw, Recentered, Twisted, Rescaled"};

  // Table.
  Produces<aod::EvtPlanes> evPlane;

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
    histosQA.add("histCentFull", "Centrality distribution for valid events",
                 HistType::kTH1F, {axisCent});

    const AxisSpec axisEP{200, -TMath::Pi() / 2., TMath::Pi() / 2.};

    for (int i = 0; i < 8; i++) {
      histosQA.add(("Centrality_0-5/histEP" + (std::string)ep::detNames[i]).c_str(),
                   ("#Psi_{2} for " + (std::string)ep::detNames[i] + (std::string)cfgCorrStep).c_str(),
                   HistType::kTH1F, {axisEP});
    }

    for (int iBin = 1; iBin < 8; iBin++) {
      histosQA.addClone("Centrality_0-5/", ep::centClasses[iBin].data());
    }

  } // void init(InitContext const&)

  template <int cBin, int det, typename T>
  void fillHistosQA(const T& val)
  {
    histosQA.fill(HIST(ep::centClasses[cBin]) + HIST("histEP") + HIST(ep::detNames[det]),
                  val[det]);
  }

  void process(aod::Qvector const& qVec)
  {
    // Get the centrality bin, and fill the centrality distribution.
    int centBin = helperEP.GetCentBin(qVec.cent());
    if (centBin < 0 || centBin > 8) {
      return;
    }
    histosQA.fill(HIST("histCentFull"), qVec.cent());

    // Calculate the event plane for each detector, then save them in the
    // corresponding distribution. The order is the same as in detNames[].
    // TODO: Update the calculation of the event plane for the central barrel.
    float evtPlaneValues[8] = {0.};
    evtPlaneValues[0] = helperEP.GetEventPlane(qVec.qvecFT0ARe(), qVec.qvecFT0AIm());
    evtPlaneValues[1] = helperEP.GetEventPlane(qVec.qvecFT0CRe(), qVec.qvecFT0CIm());
    evtPlaneValues[2] = helperEP.GetEventPlane(qVec.qvecFV0ARe(), qVec.qvecFV0AIm());
    evtPlaneValues[3] = helperEP.GetEventPlane(1., 2.);
    evtPlaneValues[4] = helperEP.GetEventPlane(2., 1.);

    evtPlaneValues[5] = helperEP.GetEventPlane(qVec.qvecFT0CUncorRe(), qVec.qvecFT0CUncorIm());
    evtPlaneValues[6] = helperEP.GetEventPlane(qVec.qvecFT0CRectrRe(), qVec.qvecFT0CRectrIm());
    evtPlaneValues[7] = helperEP.GetEventPlane(qVec.qvecFT0CTwistRe(), qVec.qvecFT0CTwistIm());

    static_for<0, 7>([&](auto iDet) {
      constexpr int indexDet = iDet.value;
      switch (centBin) {
        case 0:
          fillHistosQA<0, indexDet>(evtPlaneValues);
          break;
        case 1:
          fillHistosQA<1, indexDet>(evtPlaneValues);
          break;
        case 2:
          fillHistosQA<2, indexDet>(evtPlaneValues);
          break;
        case 3:
          fillHistosQA<3, indexDet>(evtPlaneValues);
          break;
        case 4:
          fillHistosQA<4, indexDet>(evtPlaneValues);
          break;
        case 5:
          fillHistosQA<5, indexDet>(evtPlaneValues);
          break;
        case 6:
          fillHistosQA<6, indexDet>(evtPlaneValues);
          break;
        case 7:
          fillHistosQA<7, indexDet>(evtPlaneValues);
          break;
      }
    });
    // Fill the columns of the evtPlane table.
    evPlane(qVec.cent(),
            evtPlaneValues[0], evtPlaneValues[1], evtPlaneValues[2],
            evtPlaneValues[3], evtPlaneValues[4],
            evtPlaneValues[5], evtPlaneValues[6], evtPlaneValues[7]);
  } // void process(aod::Qvector const& qVec)
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<evtPlanesTable>(cfgc)};
}
