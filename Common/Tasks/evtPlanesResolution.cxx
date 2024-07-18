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
/// \file   evtPlanesResolution.cxx
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
#include <TProfile.h>
#include <TMath.h>
#include <TH1F.h>

// o2Physics includes.
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/StaticFor.h"

#include "Common/DataModel/EvtPlanes.h"
#include "Common/Core/EventPlaneHelper.h"

// o2 includes.

using namespace o2;
using namespace o2::framework;

namespace ep
{
static constexpr std::string_view centClasses[] = {
  "Centrality_0-5/", "Centrality_5-10/", "Centrality_10-20/", "Centrality_20-30/",
  "Centrality_30-40/", "Centrality_40-50/", "Centrality_50-60/", "Centrality_60-80/"};
} // namespace ep

struct evtPlanesResolution {
  // Configurables.

  // Histogram registry for the output QA figures and list of centrality classes for it.
  // Objects are NOT saved in alphabetical orders, and registry names are NOT saved
  // as TDirectoryFile.
  HistogramRegistry histosQA{"histosQA", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};

  // Helper variables.
  EventPlaneHelper helperEP;

  Configurable<int> cfgMinTPCTracks{"cfgMinTPCTracks", 20, "minimum TPC tracks participating in Q-vector reconstruction"};
  Configurable<int> cfgnMod{"cfgnMod", 2, "Modulation of interest"};

  void init(InitContext const&)
  {
    // Fill the registry with the needed objects.
    const AxisSpec axisEvtPl{360, -constants::math::PI, constants::math::PI};

    histosQA.add("Centrality_0-5/histEvtPlUncor", "", {HistType::kTH1F, {axisEvtPl}});
    histosQA.add("Centrality_0-5/histEvtPlRectr", "", {HistType::kTH1F, {axisEvtPl}});
    histosQA.add("Centrality_0-5/histEvtPlTwist", "", {HistType::kTH1F, {axisEvtPl}});
    histosQA.add("Centrality_0-5/histEvtPlFinal", "", {HistType::kTH1F, {axisEvtPl}});

    histosQA.add("Centrality_0-5/histEvtPlTPCposUncor", "", {HistType::kTH1F, {axisEvtPl}});
    histosQA.add("Centrality_0-5/histEvtPlTPCposRectr", "", {HistType::kTH1F, {axisEvtPl}});
    histosQA.add("Centrality_0-5/histEvtPlTPCposTwist", "", {HistType::kTH1F, {axisEvtPl}});
    histosQA.add("Centrality_0-5/histEvtPlTPCposFinal", "", {HistType::kTH1F, {axisEvtPl}});

    histosQA.add("Centrality_0-5/histEvtPlTPCnegUncor", "", {HistType::kTH1F, {axisEvtPl}});
    histosQA.add("Centrality_0-5/histEvtPlTPCnegRectr", "", {HistType::kTH1F, {axisEvtPl}});
    histosQA.add("Centrality_0-5/histEvtPlTPCnegTwist", "", {HistType::kTH1F, {axisEvtPl}});
    histosQA.add("Centrality_0-5/histEvtPlTPCnegFinal", "", {HistType::kTH1F, {axisEvtPl}});

    histosQA.add("Centrality_0-5/histEvtPlResolution", "", {HistType::kTH1F, {axisEvtPl}});

    for (int iBin = 1; iBin < 8; iBin++) {
      histosQA.addClone("Centrality_0-5/", ep::centClasses[iBin].data());
    }
  } // End void init(InitContext const&)

  template <int cBin, typename T>
  void fillHistosEvtPl(const T& vec)
  {
    histosQA.fill(HIST(ep::centClasses[cBin]) + HIST("histEvtPlUncor"), vec.evtPlUncor());
    histosQA.fill(HIST(ep::centClasses[cBin]) + HIST("histEvtPlRectr"), vec.evtPlRectr());
    histosQA.fill(HIST(ep::centClasses[cBin]) + HIST("histEvtPlTwist"), vec.evtPlTwist());
    histosQA.fill(HIST(ep::centClasses[cBin]) + HIST("histEvtPlFinal"), vec.evtPlFinal());

    if (vec.nTrkTPCpos() < cfgMinTPCTracks || vec.nTrkTPCneg() < cfgMinTPCTracks)
      return;

    histosQA.fill(HIST(ep::centClasses[cBin]) + HIST("histEvtPlTPCposUncor"), vec.evtPlTPCposUncor());
    histosQA.fill(HIST(ep::centClasses[cBin]) + HIST("histEvtPlTPCposRectr"), vec.evtPlTPCposRectr());
    histosQA.fill(HIST(ep::centClasses[cBin]) + HIST("histEvtPlTPCposTwist"), vec.evtPlTPCposTwist());
    histosQA.fill(HIST(ep::centClasses[cBin]) + HIST("histEvtPlTPCposFinal"), vec.evtPlTPCposFinal());

    histosQA.fill(HIST(ep::centClasses[cBin]) + HIST("histEvtPlTPCnegUncor"), vec.evtPlTPCnegUncor());
    histosQA.fill(HIST(ep::centClasses[cBin]) + HIST("histEvtPlTPCnegRectr"), vec.evtPlTPCnegRectr());
    histosQA.fill(HIST(ep::centClasses[cBin]) + HIST("histEvtPlTPCnegTwist"), vec.evtPlTPCnegTwist());
    histosQA.fill(HIST(ep::centClasses[cBin]) + HIST("histEvtPlTPCnegFinal"), vec.evtPlTPCnegFinal());

    histosQA.fill(HIST(ep::centClasses[cBin]) + HIST("histEvtPlResolution"),
                  std::sqrt(std::cos((vec.evtPlFinal() - vec.evtPlTPCposFinal()) * cfgnMod) * std::cos((vec.evtPlFinal() - vec.evtPlTPCnegFinal()) * cfgnMod) /
                            std::cos((vec.evtPlTPCposFinal() - vec.evtPlTPCnegFinal()) * cfgnMod)));
  }

  void process(aod::EvtPlane const& evPl)
  {
    int centBin = helperEP.GetCentBin(evPl.cent());
    switch (centBin) {
      case 0:
        fillHistosEvtPl<0>(evPl);
        break;
      case 1:
        fillHistosEvtPl<1>(evPl);
        break;
      case 2:
        fillHistosEvtPl<2>(evPl);
        break;
      case 3:
        fillHistosEvtPl<3>(evPl);
        break;
      case 4:
        fillHistosEvtPl<4>(evPl);
        break;
      case 5:
        fillHistosEvtPl<5>(evPl);
        break;
      case 6:
        fillHistosEvtPl<6>(evPl);
        break;
      case 7:
        fillHistosEvtPl<7>(evPl);
        break;
    } // End switch(centBin)
  }   // End void process(...)
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<evtPlanesResolution>(cfgc)};
}
