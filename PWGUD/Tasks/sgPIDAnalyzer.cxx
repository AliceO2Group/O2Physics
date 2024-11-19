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
    {200, 0.0, 10.0},
    "Pt binning"};

  ConfigurableAxis sigmaAxis{"sigmaAxis", {1000, -20, 180}, "nSigma TPC binning"};
  ConfigurableAxis tofAxis{"tofAxis", {200, -10, 10}, "nSigma TOF binning"};
  Configurable<float> eta_min{"eta_min", -0.9, "Track Pseudorapidity"};
  Configurable<float> eta_max{"eta_max", 0.9, "Track Pseudorapidity"};

  void init(InitContext&)
  {

    const AxisSpec ptBins{ptAxis, "p_{T} axis"};
    const AxisSpec nSigmaBins{sigmaAxis, "pseudo rapidity axis"};
    const AxisSpec ntofBins{tofAxis, "pseudo rapidity axis"};
    histos.add("Events", "Selected Events", {HistType::kTH1F, {3, -.5, 2.5}});
    histos.add("TPC/pTPC_Pi", "Positive TPC Pi Tracks", {HistType::kTH2F, {ptBins, nSigmaBins}});
    histos.add("TPC/nTPC_Pi", "Negative TPC Pi Tracks", {HistType::kTH2F, {ptBins, nSigmaBins}});
    histos.add("TPC/pTPC_Ka", "Positive TPC Ka Tracks", {HistType::kTH2F, {ptBins, nSigmaBins}});
    histos.add("TPC/nTPC_Ka", "Negative TPC Ka Tracks", {HistType::kTH2F, {ptBins, nSigmaBins}});
    histos.add("TPC/pTPC_Pr", "Positive TPC Pr Tracks", {HistType::kTH2F, {ptBins, nSigmaBins}});
    histos.add("TPC/nTPC_Pr", "Negative TPC Pr Tracks", {HistType::kTH2F, {ptBins, nSigmaBins}});
    histos.add("TPC/pTPC_El", "Positive TPC El Tracks", {HistType::kTH2F, {ptBins, nSigmaBins}});
    histos.add("TPC/nTPC_El", "Negative TPC El Tracks", {HistType::kTH2F, {ptBins, nSigmaBins}});
    histos.add("TPC/pTPC_De", "Positive TPC De Tracks", {HistType::kTH2F, {ptBins, nSigmaBins}});
    histos.add("TPC/nTPC_De", "Negative TPC De Tracks", {HistType::kTH2F, {ptBins, nSigmaBins}});
    histos.add("TPC/pTPC_Tr", "Positive TPC Tr Tracks", {HistType::kTH2F, {ptBins, nSigmaBins}});
    histos.add("TPC/nTPC_Tr", "Negative TPC Tr Tracks", {HistType::kTH2F, {ptBins, nSigmaBins}});
    histos.add("TPC/pTPC_He", "Positive TPC He Tracks", {HistType::kTH2F, {ptBins, nSigmaBins}});
    histos.add("TPC/nTPC_He", "Negative TPC He Tracks", {HistType::kTH2F, {ptBins, nSigmaBins}});
    histos.add("TPC/pTPC_Al", "Positive TPC Al Tracks", {HistType::kTH2F, {ptBins, nSigmaBins}});
    histos.add("TPC/nTPC_Al", "Negative TPC Al Tracks", {HistType::kTH2F, {ptBins, nSigmaBins}});
    histos.add("TPC/pTPC_Mu", "Positive TPC Mu Tracks", {HistType::kTH2F, {ptBins, nSigmaBins}});
    histos.add("TPC/nTPC_Mu", "Negative TPC Mu Tracks", {HistType::kTH2F, {ptBins, nSigmaBins}});

    histos.add("TPC/pTPC_Pi_Ka", "Positive TPC Pi vs Ka", {HistType::kTH2F, {ptBins, nSigmaBins}});
    histos.add("TPC/pTPC_Pi_Pr", "Positive TPC Pi vs Pr", {HistType::kTH2F, {ptBins, nSigmaBins}});
    histos.add("TPC/pTPC_Pi_El", "Positive TPC Pi vs El", {HistType::kTH2F, {ptBins, nSigmaBins}});
    histos.add("TPC/pTPC_Ka_Pi", "Positive TPC Ka vs Pi", {HistType::kTH2F, {ptBins, nSigmaBins}});
    histos.add("TPC/pTPC_Ka_Pr", "Positive TPC Ka vs Pr", {HistType::kTH2F, {ptBins, nSigmaBins}});
    histos.add("TPC/pTPC_Ka_El", "Positive TPC Ka vs El", {HistType::kTH2F, {ptBins, nSigmaBins}});
    histos.add("TPC/pTPC_Pr_Pi", "Positive TPC Pr vs Pi", {HistType::kTH2F, {ptBins, nSigmaBins}});
    histos.add("TPC/pTPC_Pr_Ka", "Positive TPC Pr vs Ka", {HistType::kTH2F, {ptBins, nSigmaBins}});
    histos.add("TPC/pTPC_Pr_El", "Positive TPC Pr vs El", {HistType::kTH2F, {ptBins, nSigmaBins}});
    histos.add("TPC/pTPC_El_Pi", "Positive TPC Pr vs Pi", {HistType::kTH2F, {ptBins, nSigmaBins}});
    histos.add("TPC/pTPC_El_Ka", "Positive TPC Pr vs Ka", {HistType::kTH2F, {ptBins, nSigmaBins}});
    histos.add("TPC/pTPC_El_Pr", "Positive TPC Pr vs El", {HistType::kTH2F, {ptBins, nSigmaBins}});

    histos.add("TPC/nTPC_Pi_Ka", "Positive TPC Pi vs Ka", {HistType::kTH2F, {ptBins, nSigmaBins}});
    histos.add("TPC/nTPC_Pi_Pr", "Positive TPC Pi vs Pr", {HistType::kTH2F, {ptBins, nSigmaBins}});
    histos.add("TPC/nTPC_Pi_El", "Positive TPC Pi vs El", {HistType::kTH2F, {ptBins, nSigmaBins}});
    histos.add("TPC/nTPC_Ka_Pi", "Positive TPC Ka vs Pi", {HistType::kTH2F, {ptBins, nSigmaBins}});
    histos.add("TPC/nTPC_Ka_Pr", "Positive TPC Ka vs Pr", {HistType::kTH2F, {ptBins, nSigmaBins}});
    histos.add("TPC/nTPC_Ka_El", "Positive TPC Ka vs El", {HistType::kTH2F, {ptBins, nSigmaBins}});
    histos.add("TPC/nTPC_Pr_Pi", "Positive TPC Pr vs Pi", {HistType::kTH2F, {ptBins, nSigmaBins}});
    histos.add("TPC/nTPC_Pr_Ka", "Positive TPC Pr vs Ka", {HistType::kTH2F, {ptBins, nSigmaBins}});
    histos.add("TPC/nTPC_Pr_El", "Positive TPC Pr vs El", {HistType::kTH2F, {ptBins, nSigmaBins}});
    histos.add("TPC/nTPC_El_Pi", "Positive TPC Pr vs Pi", {HistType::kTH2F, {ptBins, nSigmaBins}});
    histos.add("TPC/nTPC_El_Ka", "Positive TPC Pr vs Ka", {HistType::kTH2F, {ptBins, nSigmaBins}});
    histos.add("TPC/nTPC_El_Pr", "Positive TPC Pr vs El", {HistType::kTH2F, {ptBins, nSigmaBins}});

    histos.add("TOF/pPi", "Positive TPC Pi vs TOF Pi vs pt", {HistType::kTH3F, {ptBins, ntofBins, ntofBins}});
    histos.add("TOF/nPi", "Negative TPC Pi vs TOF Pi vs pt", {HistType::kTH3F, {ptBins, ntofBins, ntofBins}});
    histos.add("TOF/pKa", "Positive TPC Ka vs TOF Ka vs pt", {HistType::kTH3F, {ptBins, ntofBins, ntofBins}});
    histos.add("TOF/nKa", "Negative TPC Ka vs TOF Ka vs pt", {HistType::kTH3F, {ptBins, ntofBins, ntofBins}});
    histos.add("TOF/pPr", "Positive TPC Pr vs TOF Pr vs pt", {HistType::kTH3F, {ptBins, ntofBins, ntofBins}});
    histos.add("TOF/nPr", "Negative TPC Pr vs TOF Pr vs pt", {HistType::kTH3F, {ptBins, ntofBins, ntofBins}});
    histos.add("TOF/pEl", "Positive TPC El vs TOF El vs pt", {HistType::kTH3F, {ptBins, ntofBins, ntofBins}});
    histos.add("TOF/nEl", "Negative TPC El vs TOF El vs pt", {HistType::kTH3F, {ptBins, ntofBins, ntofBins}});
    histos.add("TOF/pDe", "Positive TPC De vs TOF Pi vs pt", {HistType::kTH3F, {ptBins, ntofBins, ntofBins}});
    histos.add("TOF/nDe", "Negative TPC De vs TOF Pi vs pt", {HistType::kTH3F, {ptBins, ntofBins, ntofBins}});
    histos.add("TOF/pTr", "Positive TPC Tr vs TOF Ka vs pt", {HistType::kTH3F, {ptBins, ntofBins, ntofBins}});
    histos.add("TOF/nTr", "Negative TPC Tr vs TOF Ka vs pt", {HistType::kTH3F, {ptBins, ntofBins, ntofBins}});
    histos.add("TOF/pHe", "Positive TPC He vs TOF Pr vs pt", {HistType::kTH3F, {ptBins, ntofBins, ntofBins}});
    histos.add("TOF/nHe", "Negative TPC He vs TOF Pr vs pt", {HistType::kTH3F, {ptBins, ntofBins, ntofBins}});
    histos.add("TOF/pAl", "Positive TPC Al vs TOF El vs pt", {HistType::kTH3F, {ptBins, ntofBins, ntofBins}});
    histos.add("TOF/nAl", "Negative TPC Al vs TOF El vs pt", {HistType::kTH3F, {ptBins, ntofBins, ntofBins}});
    histos.add("TOF/pMu", "Positive TPC Mu vs TOF El vs pt", {HistType::kTH3F, {ptBins, ntofBins, ntofBins}});
    histos.add("TOF/nMu", "Negative TPC Mu vs TOF El vs pt", {HistType::kTH3F, {ptBins, ntofBins, ntofBins}});
  }
  using SGEvent = aod::SGEvents::iterator;
  void process(SGEvent const& event, aod::SGTracks const& tracks)
  {
    histos.fill(HIST("Events"), event.gs());
    for (const auto& track : tracks) {
      if (track.eta() < eta_min || track.eta() > eta_max)
        continue;
      bool isPositive = (track.sign() > 0);
      if (track.tofpi() == -999) {
        // Directly fill histograms without a local variable for histName
        if (isPositive) {
          histos.fill(HIST("TPC/pTPC_Pi"), track.pt(), track.tpcpi());
          histos.fill(HIST("TPC/pTPC_Ka"), track.pt(), track.tpcka());
          histos.fill(HIST("TPC/pTPC_Pr"), track.pt(), track.tpcpr());
          histos.fill(HIST("TPC/pTPC_El"), track.pt(), track.tpcel());
          histos.fill(HIST("TPC/pTPC_De"), track.pt(), track.tpcde());
          histos.fill(HIST("TPC/pTPC_Tr"), track.pt(), track.tpctr());
          histos.fill(HIST("TPC/pTPC_He"), track.pt(), track.tpche());
          histos.fill(HIST("TPC/pTPC_Al"), track.pt(), track.tpcal());
          histos.fill(HIST("TPC/pTPC_Mu"), track.pt(), track.tpcmu());
          if (std::abs(track.tpcpi()) < 1) {
            histos.fill(HIST("TPC/pTPC_Ka_Pi"), track.pt(), track.tpcka());
            histos.fill(HIST("TPC/pTPC_Pr_Pi"), track.pt(), track.tpcpr());
            histos.fill(HIST("TPC/pTPC_El_Pi"), track.pt(), track.tpcel());
          }
          if (std::abs(track.tpcka()) < 1) {
            histos.fill(HIST("TPC/pTPC_Pi_Ka"), track.pt(), track.tpcpi());
            histos.fill(HIST("TPC/pTPC_Pr_Ka"), track.pt(), track.tpcpr());
            histos.fill(HIST("TPC/pTPC_El_Ka"), track.pt(), track.tpcel());
          }
          if (std::abs(track.tpcpr()) < 1) {
            histos.fill(HIST("TPC/pTPC_Pi_Pr"), track.pt(), track.tpcpi());
            histos.fill(HIST("TPC/pTPC_Ka_Pr"), track.pt(), track.tpcka());
            histos.fill(HIST("TPC/pTPC_El_Pr"), track.pt(), track.tpcel());
          }
          if (std::abs(track.tpcel()) < 1) {
            histos.fill(HIST("TPC/pTPC_Pi_El"), track.pt(), track.tpcpi());
            histos.fill(HIST("TPC/pTPC_Ka_El"), track.pt(), track.tpcka());
            histos.fill(HIST("TPC/pTPC_Pr_El"), track.pt(), track.tpcpr());
          }
        } else {
          histos.fill(HIST("TPC/nTPC_Pi"), track.pt(), track.tpcpi());
          histos.fill(HIST("TPC/nTPC_Ka"), track.pt(), track.tpcka());
          histos.fill(HIST("TPC/nTPC_Pr"), track.pt(), track.tpcpr());
          histos.fill(HIST("TPC/nTPC_El"), track.pt(), track.tpcel());
          histos.fill(HIST("TPC/nTPC_De"), track.pt(), track.tpcde());
          histos.fill(HIST("TPC/nTPC_Tr"), track.pt(), track.tpctr());
          histos.fill(HIST("TPC/nTPC_He"), track.pt(), track.tpche());
          histos.fill(HIST("TPC/nTPC_Al"), track.pt(), track.tpcal());
          histos.fill(HIST("TPC/nTPC_Mu"), track.pt(), track.tpcmu());
          if (std::abs(track.tpcpi()) < 1) {
            histos.fill(HIST("TPC/nTPC_Ka_Pi"), track.pt(), track.tpcka());
            histos.fill(HIST("TPC/nTPC_Pr_Pi"), track.pt(), track.tpcpr());
            histos.fill(HIST("TPC/nTPC_El_Pi"), track.pt(), track.tpcel());
          }
          if (std::abs(track.tpcka()) < 1) {
            histos.fill(HIST("TPC/nTPC_Pi_Ka"), track.pt(), track.tpcpi());
            histos.fill(HIST("TPC/nTPC_Pr_Ka"), track.pt(), track.tpcpr());
            histos.fill(HIST("TPC/nTPC_El_Ka"), track.pt(), track.tpcel());
          }
          if (std::abs(track.tpcpr()) < 1) {
            histos.fill(HIST("TPC/nTPC_Pi_Pr"), track.pt(), track.tpcpi());
            histos.fill(HIST("TPC/nTPC_Ka_Pr"), track.pt(), track.tpcka());
            histos.fill(HIST("TPC/nTPC_El_Pr"), track.pt(), track.tpcel());
          }
          if (std::abs(track.tpcel()) < 1) {
            histos.fill(HIST("TPC/nTPC_Pi_El"), track.pt(), track.tpcpi());
            histos.fill(HIST("TPC/nTPC_Ka_El"), track.pt(), track.tpcka());
            histos.fill(HIST("TPC/nTPC_Pr_El"), track.pt(), track.tpcpr());
          }
        }
      } else {
        if (isPositive) {
          histos.fill(HIST("TOF/pPi"), track.pt(), track.tpcpi(), track.tofpi());
          histos.fill(HIST("TOF/pKa"), track.pt(), track.tpcka(), track.tofka());
          histos.fill(HIST("TOF/pPr"), track.pt(), track.tpcpr(), track.tofpr());
          histos.fill(HIST("TOF/pEl"), track.pt(), track.tpcel(), track.tofel());
          histos.fill(HIST("TOF/pDe"), track.pt(), track.tpcpi(), track.tofde());
          histos.fill(HIST("TOF/pTr"), track.pt(), track.tpcka(), track.toftr());
          histos.fill(HIST("TOF/pHe"), track.pt(), track.tpcpr(), track.tofhe());
          histos.fill(HIST("TOF/pAl"), track.pt(), track.tpcel(), track.tofal());
          histos.fill(HIST("TOF/pMu"), track.pt(), track.tpcel(), track.tofmu());
        } else {
          histos.fill(HIST("TOF/nPi"), track.pt(), track.tpcpi(), track.tofpi());
          histos.fill(HIST("TOF/nKa"), track.pt(), track.tpcka(), track.tofka());
          histos.fill(HIST("TOF/nPr"), track.pt(), track.tpcpr(), track.tofpr());
          histos.fill(HIST("TOF/nEl"), track.pt(), track.tpcel(), track.tofel());
          histos.fill(HIST("TOF/nDe"), track.pt(), track.tpcpi(), track.tofde());
          histos.fill(HIST("TOF/nTr"), track.pt(), track.tpcka(), track.toftr());
          histos.fill(HIST("TOF/nHe"), track.pt(), track.tpcpr(), track.tofhe());
          histos.fill(HIST("TOF/nAl"), track.pt(), track.tpcel(), track.tofal());
          histos.fill(HIST("TOF/nMu"), track.pt(), track.tpcel(), track.tofmu());
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
