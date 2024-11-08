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
//
// \Single Gap Event Analyzer
// \author Sasha Bylinkin, alexander.bylinkin@gmail.com
// \since  April 2023

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"

#include "TVector3.h"
#include "TTree.h"
#include "TFile.h"
#include <TH1F.h>
#include <TH2F.h>
#include "Common/DataModel/PIDResponse.h"
#include "PWGUD/DataModel/SGTables.h"
#include "PWGUD/Core/UDHelpers.h"
#include "PWGUD/Core/SGSelector.h"
#include "PWGUD/Core/SGTrackSelector.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct sgPIDAnalyzer {
  HistogramRegistry histos{"Histos", {}};

  ConfigurableAxis ptAxis{
    "ptAxis",
    {198, 0.1, 10.0},
    "Pt binning"};

  ConfigurableAxis sigmaAxis{"sigmaAxis", {100, -50, 50}, "nSigma TPC binning"};

  void init(InitContext&)
  {

    const AxisSpec ptBins{ptAxis, "p_{T} axis"};
    const AxisSpec nSigmaBins{sigmaAxis, "pseudo rapidity axis"};
    histos.add("TPC/pTPC_Pi", "Positive TPC Pi Tracks", {HistType::kTH2F, {ptBins, nSigmaBins}});
    histos.add("TPC/nTPC_Pi", "Negative TPC Pi Tracks", {HistType::kTH2F, {ptBins, nSigmaBins}});
    histos.add("TPC/pTPC_Ka", "Positive TPC Ka Tracks", {HistType::kTH2F, {ptBins, nSigmaBins}});
    histos.add("TPC/nTPC_Ka", "Negative TPC Ka Tracks", {HistType::kTH2F, {ptBins, nSigmaBins}});
    histos.add("TPC/pTPC_Pr", "Positive TPC Pr Tracks", {HistType::kTH2F, {ptBins, nSigmaBins}});
    histos.add("TPC/nTPC_Pr", "Negative TPC Pr Tracks", {HistType::kTH2F, {ptBins, nSigmaBins}});
    histos.add("TPC/pTPC_El", "Positive TPC El Tracks", {HistType::kTH2F, {ptBins, nSigmaBins}});
    histos.add("TPC/nTPC_El", "Negative TPC El Tracks", {HistType::kTH2F, {ptBins, nSigmaBins}});

    histos.add("TPC/pTPC_Pi_Ka", "Positive TPC Pi vs Ka", {HistType::kTH2F, {ptBins, nSigmaBins}});
    histos.add("TPC/pTPC_Pi_Pr", "Positive TPC Pi vs Pr", {HistType::kTH2F, {ptBins, nSigmaBins}});
    histos.add("TPC/pTPC_Pi_El", "Positive TPC Pi vs El", {HistType::kTH2F, {ptBins, nSigmaBins}});
    histos.add("TPC/pTPC_Ka_Pi", "Positive TPC Ka vs Pi", {HistType::kTH2F, {ptBins, nSigmaBins}});
    histos.add("TPC/pTPC_Ka_Pr", "Positive TPC Ka vs Pr", {HistType::kTH2F, {ptBins, nSigmaBins}});
    histos.add("TPC/pTPC_Ka_El", "Positive TPC Ka vs El", {HistType::kTH2F, {ptBins, nSigmaBins}});
    histos.add("TPC/pTPC_Pr_Pi", "Positive TPC Pr vs Pi", {HistType::kTH2F, {ptBins, nSigmaBins}});
    histos.add("TPC/pTPC_Pr_Ka", "Positive TPC Pr vs Ka", {HistType::kTH2F, {ptBins, nSigmaBins}});
    histos.add("TPC/pTPC_Pr_El", "Positive TPC Pr vs El", {HistType::kTH2F, {ptBins, nSigmaBins}});

    histos.add("TPC/nTPC_Pi_Ka", "Positive TPC Pi vs Ka", {HistType::kTH2F, {ptBins, nSigmaBins}});
    histos.add("TPC/nTPC_Pi_Pr", "Positive TPC Pi vs Pr", {HistType::kTH2F, {ptBins, nSigmaBins}});
    histos.add("TPC/nTPC_Pi_El", "Positive TPC Pi vs El", {HistType::kTH2F, {ptBins, nSigmaBins}});
    histos.add("TPC/nTPC_Ka_Pi", "Positive TPC Ka vs Pi", {HistType::kTH2F, {ptBins, nSigmaBins}});
    histos.add("TPC/nTPC_Ka_Pr", "Positive TPC Ka vs Pr", {HistType::kTH2F, {ptBins, nSigmaBins}});
    histos.add("TPC/nTPC_Ka_El", "Positive TPC Ka vs El", {HistType::kTH2F, {ptBins, nSigmaBins}});
    histos.add("TPC/nTPC_Pr_Pi", "Positive TPC Pr vs Pi", {HistType::kTH2F, {ptBins, nSigmaBins}});
    histos.add("TPC/nTPC_Pr_Ka", "Positive TPC Pr vs Ka", {HistType::kTH2F, {ptBins, nSigmaBins}});
    histos.add("TPC/nTPC_Pr_El", "Positive TPC Pr vs El", {HistType::kTH2F, {ptBins, nSigmaBins}});

    histos.add("TOF/pPi", "Positive TPC Pi vs TOF Pi vs pt", {HistType::kTH3F, {ptBins, nSigmaBins, nSigmaBins}});
    histos.add("TOF/nPi", "Negative TPC Pi vs TOF Pi vs pt", {HistType::kTH3F, {ptBins, nSigmaBins, nSigmaBins}});
    histos.add("TOF/pKa", "Positive TPC Ka vs TOF Ka vs pt", {HistType::kTH3F, {ptBins, nSigmaBins, nSigmaBins}});
    histos.add("TOF/nKa", "Negative TPC Ka vs TOF Ka vs pt", {HistType::kTH3F, {ptBins, nSigmaBins, nSigmaBins}});
    histos.add("TOF/pPr", "Positive TPC Pr vs TOF Pr vs pt", {HistType::kTH3F, {ptBins, nSigmaBins, nSigmaBins}});
    histos.add("TOF/nPr", "Negative TPC Pr vs TOF Pr vs pt", {HistType::kTH3F, {ptBins, nSigmaBins, nSigmaBins}});
    histos.add("TOF/pEl", "Positive TPC El vs TOF El vs pt", {HistType::kTH3F, {ptBins, nSigmaBins, nSigmaBins}});
    histos.add("TOF/nEl", "Negative TPC El vs TOF El vs pt", {HistType::kTH3F, {ptBins, nSigmaBins, nSigmaBins}});
  }

  void process(aod::SGEvents const& events, aod::SGTracks const& tracks)
  {
    for (const auto& track : tracks) {
      bool isPositive = (track.sign() > 0);
      if (track.tofpi() == -999) {
        // Directly fill histograms without a local variable for histName
        if (isPositive) {
          histos.fill(HIST("TPC/pTPC_Pi"), track.pt(), track.tpcpi());
          histos.fill(HIST("TPC/pTPC_Ka"), track.pt(), track.tpcka());
          histos.fill(HIST("TPC/pTPC_Pr"), track.pt(), track.tpcpr());
          histos.fill(HIST("TPC/pTPC_El"), track.pt(), track.tpcel());
          if (std::abs(track.tpcpi()) < 1) {
            histos.fill(HIST("TPC/pTPC_Pi_Ka"), track.pt(), track.tpcka());
            histos.fill(HIST("TPC/pTPC_Pi_Pr"), track.pt(), track.tpcpr());
            histos.fill(HIST("TPC/pTPC_Pi_El"), track.pt(), track.tpcel());
          }
          if (std::abs(track.tpcka()) < 1) {
            histos.fill(HIST("TPC/pTPC_Ka_Pi"), track.pt(), track.tpcpi());
            histos.fill(HIST("TPC/pTPC_Ka_Pr"), track.pt(), track.tpcpr());
            histos.fill(HIST("TPC/pTPC_Ka_El"), track.pt(), track.tpcel());
          }
          if (std::abs(track.tpcpr()) < 1) {
            histos.fill(HIST("TPC/pTPC_Pr_Pi"), track.pt(), track.tpcpi());
            histos.fill(HIST("TPC/pTPC_Pr_Ka"), track.pt(), track.tpcka());
            histos.fill(HIST("TPC/pTPC_Pr_El"), track.pt(), track.tpcel());
          }
        } else {
          histos.fill(HIST("TPC/nTPC_Pi"), track.pt(), track.tpcpi());
          histos.fill(HIST("TPC/nTPC_Ka"), track.pt(), track.tpcka());
          histos.fill(HIST("TPC/nTPC_Pr"), track.pt(), track.tpcpr());
          histos.fill(HIST("TPC/nTPC_El"), track.pt(), track.tpcel());
          if (std::abs(track.tpcpi()) < 1) {
            histos.fill(HIST("TPC/nTPC_Pi_Ka"), track.pt(), track.tpcka());
            histos.fill(HIST("TPC/nTPC_Pi_Pr"), track.pt(), track.tpcpr());
            histos.fill(HIST("TPC/nTPC_Pi_El"), track.pt(), track.tpcel());
          }
          if (std::abs(track.tpcka()) < 1) {
            histos.fill(HIST("TPC/nTPC_Ka_Pi"), track.pt(), track.tpcpi());
            histos.fill(HIST("TPC/nTPC_Ka_Pr"), track.pt(), track.tpcpr());
            histos.fill(HIST("TPC/nTPC_Ka_El"), track.pt(), track.tpcel());
          }
          if (std::abs(track.tpcpr()) < 1) {
            histos.fill(HIST("TPC/nTPC_Pr_Pi"), track.pt(), track.tpcpi());
            histos.fill(HIST("TPC/nTPC_Pr_Ka"), track.pt(), track.tpcka());
            histos.fill(HIST("TPC/nTPC_Pr_El"), track.pt(), track.tpcel());
          }
        }
      } else {
        if (isPositive) {
          histos.fill(HIST("TOF/pPi"), track.pt(), track.tpcpi(), track.tofpi());
          histos.fill(HIST("TOF/pKa"), track.pt(), track.tpcka(), track.tofka());
          histos.fill(HIST("TOF/pPr"), track.pt(), track.tpcpr(), track.tofpr());
          histos.fill(HIST("TOF/pEl"), track.pt(), track.tpcel(), track.tofel());
        } else {
          histos.fill(HIST("TOF/nPi"), track.pt(), track.tpcpi(), track.tofpi());
          histos.fill(HIST("TOF/nKa"), track.pt(), track.tpcka(), track.tofka());
          histos.fill(HIST("TOF/nPr"), track.pt(), track.tpcpr(), track.tofpr());
          histos.fill(HIST("TOF/nEl"), track.pt(), track.tpcel(), track.tofel());
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<sgPIDAnalyzer>(cfgc, TaskName{"sgpidanalyzer"})};
}
