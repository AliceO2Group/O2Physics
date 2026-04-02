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

/// \file taskCorrelationDstarHadrons.cxx
/// \author Deependra Sharma <deependra.sharma@cern.ch>, IITB
/// \author Fabrizio Grosa <fabrizio.grosa@cern.ch>, CERN
/// \author Shyam Kumar <shyam.kumar@cern.ch>

#include "PWGHF/Core/SelectorCuts.h"
#include "PWGHF/HFC/DataModel/CorrelationTables.h"
#include "PWGHF/Utils/utilsAnalysis.h"

#include <CCDB/CcdbApi.h>
#include <CommonConstants/MathConstants.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <TFile.h>
#include <TH1.h>
#include <TString.h>

#include <cstdint>
#include <map>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// string definitions, used for histogram axis labels
const TString stringPtD = "#it{p}_{T}^{D} (GeV/#it{c});";
const TString stringPtHadron = "#it{p}_{T}^{Hadron} (GeV/#it{c});";
const TString stringDeltaEta = "#it{#eta}^{Hadron}-#it{#eta}^{D};";
const TString stringDeltaPhi = "#it{#varphi}^{Hadron}-#it{#varphi}^{D} (rad);";
const TString stringDHadron = "D,Hadron candidates ";
const TString stringSignal = "signal region;";
const TString stringSideband = "sidebands;";
const TString stringPoolBin = "Pool Bin Number;";

const int nBinsPtCorrelation = 8;

const double binsPtCorrelationsDefault[nBinsPtCorrelation + 1] = {0., 2., 4., 6., 8., 12., 16., 24., 100.};
const auto vecBinsPtCorrelationsDefault = std::vector<double>{binsPtCorrelationsDefault, binsPtCorrelationsDefault + nBinsPtCorrelation + 1};

const double signalRegionLefBoundDefault[nBinsPtCorrelation] = {0.144, 0.144, 0.144, 0.144, 0.144, 0.144, 0.144, 0.144};
const auto vecSignalRegionLefBoundDefault = std::vector<double>{signalRegionLefBoundDefault, signalRegionLefBoundDefault + nBinsPtCorrelation};

const double signalRegionRightBoundDefault[nBinsPtCorrelation] = {0.146, 0.146, 0.146, 0.146, 0.146, 0.146, 0.146, 0.146};
const auto vecSignalRegionRightBoundDefault = std::vector<double>{signalRegionRightBoundDefault, signalRegionRightBoundDefault + nBinsPtCorrelation};

// const double sidebandLeftOuterDefault[nBinsPtCorrelation] = {1.7690, 1.7690, 1.7690, 1.7690, 1.7690, 1.7690, 1.7690, 1.7690};
// const auto vecSidebandLeftOuterDefault = std::vector<double>{sidebandLeftOuterDefault, sidebandLeftOuterDefault + nBinsPtCorrelation};

// const double sidebandLeftInnerDefault[nBinsPtCorrelation] = {1.8250, 1.8250, 1.8250, 1.8250, 1.8250, 1.8250, 1.8250, 1.8250};
// const auto vecSidebandLeftInnerDefault = std::vector<double>{sidebandLeftInnerDefault, sidebandLeftInnerDefault + nBinsPtCorrelation};

const double sidebandRightInnerDefault[nBinsPtCorrelation] = {0.147, 0.147, 0.147, 0.147, 0.147, 0.147, 0.147, 0.147};
const auto vecSidebandRightInnerDefault = std::vector<double>{sidebandRightInnerDefault, sidebandRightInnerDefault + nBinsPtCorrelation};

const double sidebandRightOuterDefault[nBinsPtCorrelation] = {0.154, 0.154, 0.154, 0.154, 0.154, 0.154, 0.154, 0.154};
const auto vecSidebandRightOuterDefault = std::vector<double>{sidebandRightOuterDefault, sidebandRightOuterDefault + nBinsPtCorrelation};

const int npTBinsEfficiency = o2::analysis::hf_cuts_dstar_to_d0_pi::NBinsPt;
const std::vector<double> vecEfficiencyDstarDefault(npTBinsEfficiency); // line # 76 in taskCorrelationDstarHadron.cxx; why (npTBinsEfficiency+1) ?

const int nPtBinsTrackEfficiency = o2::analysis::hf_cuts_single_track::NBinsPtTrack;
const std::vector<double> vecEfficiencyTracksDefault(nPtBinsTrackEfficiency);

// Dstar-Hadron correlation pair
struct HfTaskCorrelationDstarHadrons {

  Configurable<bool> applyEfficiency{"applyEfficiency", true, "Flag for applying efficiency weights"};
  Configurable<bool> useCcdbEfficiency{"useCcdbEfficiency", false, "Flag for using efficiency values from CCDB (if false, efficiency values must be provided via json files)"};
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> ccdbPathEfficiencyDstar{"ccdbPathEfficiencyDstar", "Users/d/desharma/HFC/Efficiency/Dstar", "path in ccdb for Dstar efficiency values"};
  Configurable<std::string> ccdbPathEfficiencyTracks{"ccdbPathEfficiencyTracks", "Users/d/desharma/HFC/Efficiency/Track", "path in ccdb for track efficiency values"};
  Configurable<int64_t> ccdbTimestamp{"ccdbTimestamp", -1, "timestamp for retrieving efficiency values from CCDB"};
  Configurable<std::string> efficiencyDstarFileName{"efficiencyDstarFileName", "efficiencyHFCDstar.root", "name of the efficiency file for Dstar"};
  Configurable<std::string> efficiencyTracksFileName{"efficiencyTracksFileName", "efficiencyHFCTrack.root", "name of the efficiency file for tracks"};
  Configurable<int> nEfficiencyHist{"nEfficiencyHist", 1, "if MB nEfficiencyHist = 1, if Centrality classes nEfficiencyHist = number of centrality classes (i.e. 10)"};

  // pT ranges for correlation plots: the default values are those embedded in hf_cuts_dplus_to_pi_k_pi (i.e. the mass pT bins), but can be redefined via json files
  Configurable<std::vector<double>> binsPtCorrelations{"binsPtCorrelations", std::vector<double>{vecBinsPtCorrelationsDefault}, "pT bin limits for correlation plots"};

  // efficiency configurables for candidate Dstar
  Configurable<std::vector<double>> binsPtEfficiency{"binsPtEfficiency", std::vector<double>{o2::analysis::hf_cuts_dstar_to_d0_pi::vecBinsPt}, "pT bin limits for efficiency"};
  Configurable<std::vector<double>> efficiencyDstar{"efficiencyDstar", std::vector<double>{vecEfficiencyDstarDefault}, "efficiency values for Dstar vs pT bin"};

  // efficiency configurables for associated tracks
  Configurable<std::vector<double>> binsPtEfficiencyTracks{"binsPtEfficiencyTracks", std::vector<double>{o2::analysis::hf_cuts_single_track::vecBinsPtTrack}, "pT bin limits for track efficiency"};
  Configurable<std::vector<double>> efficiencyTracks{"efficiencyTracks", std::vector<double>{vecEfficiencyTracksDefault}, "efficiency values for tracks vs pT bin"};

  Configurable<std::vector<double>> signalRegionLefBound{"signalRegionLefBound", std::vector<double>{vecSignalRegionLefBoundDefault}, "left boundary of signal region vs pT"};
  Configurable<std::vector<double>> signalRegionRightBound{"signalRegionRightBound", std::vector<double>{vecSignalRegionRightBoundDefault}, "right boundary of signal region vs pT"};
  // Configurable<std::vector<double>> leftSidebandOuterBoundary{"leftSidebandOuterBoundary", std::vector<double>{vecSidebandLeftOuterDefault}, "left sideband outer boundary vs pT"};
  // Configurable<std::vector<double>> leftSidebandInnerBoundary{"leftSidebandInnerBoundary", std::vector<double>{vecSidebandLeftInnerDefault}, "left sideband inner boundary vs pT"};
  Configurable<std::vector<double>> rightSidebandOuterBoundary{"rightSidebandOuterBoundary", std::vector<double>{vecSidebandRightOuterDefault}, "right sideband outer baoundary vs pT"};
  Configurable<std::vector<double>> rightSidebandInnerBoundary{"rightSidebandInnerBoundary", std::vector<double>{vecSidebandRightInnerDefault}, "right sideband inner boundary"};
  Configurable<int> nBinsDeltaPhi{"nBinsDeltaPhi", 64, "number of bins in delta phi axis"};

  ConfigurableAxis deltaEtaBinEdges{"deltaEtaBinEdges", {40, -2., 2.}, " Delta Eta Bins of equal width"};
  ConfigurableAxis ptHadronBinsEdges{"ptHadronBinsEdges", {11, 0., 11.}, "pT Bins of equal width for Hadrons"};

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  o2::ccdb::CcdbApi ccdbApi;
  std::vector<TH1F*> vecHistEfficiencyDstar;
  std::vector<TH1D*> vecHistEfficiencyTracks;

  void init(InitContext&)
  {

    auto axisPtDstar = (std::vector<double>)binsPtEfficiency;
    AxisSpec const axisSpecPtDstar = {axisPtDstar};
    AxisSpec const axisSpecDeltaPhi = {nBinsDeltaPhi, -o2::constants::math::PIHalf, 3. * o2::constants::math::PIHalf};
    AxisSpec const axisSpecDeltaEta = {deltaEtaBinEdges};
    AxisSpec const axisSpecPtHadron = {ptHadronBinsEdges};
    AxisSpec const axisSpecPoolBin = {9, 0., 9.};

    registry.add("hCorrel2DVsPtSignalRegion", stringDHadron + stringSignal + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtHadron + stringPoolBin + "entries", {HistType::kTHnSparseD, {axisSpecDeltaPhi, axisSpecDeltaEta, axisSpecPtDstar, axisSpecPtHadron, axisSpecPoolBin}}, true);
    registry.add("hCorrel2DPtIntSignalRegion", stringDHadron + stringSignal + stringDeltaPhi + stringDeltaEta + "entries", {HistType::kTH2D, {axisSpecDeltaPhi, axisSpecDeltaEta}}, true);
    registry.add("hDeltaEtaPtIntSignalRegion", stringDHadron + stringSignal + stringDeltaEta + "entries", {HistType::kTH1D, {axisSpecDeltaEta}}, true);
    registry.add("hDeltaPhiPtIntSignalRegion", stringDHadron + stringSignal + stringDeltaPhi + "entries", {HistType::kTH1D, {axisSpecDeltaPhi}}, true);
    registry.add("hCorrel2DVsPtSidebands", stringDHadron + stringSideband + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtHadron + stringPoolBin + "entries", {HistType::kTHnSparseD, {axisSpecDeltaPhi, axisSpecDeltaEta, axisSpecPtDstar, axisSpecPtHadron, axisSpecPoolBin}}, true);
    registry.add("hCorrel2DPtIntSidebands", stringDHadron + stringSideband + stringDeltaPhi + stringDeltaEta + "entries", {HistType::kTH2D, {axisSpecDeltaPhi, axisSpecDeltaEta}}, true);
    registry.add("hDeltaEtaPtIntSidebands", stringDHadron + stringSideband + stringDeltaEta + "entries", {HistType::kTH1D, {axisSpecDeltaEta}}, true);
    registry.add("hDeltaPhiPtIntSidebands", stringDHadron + stringSideband + stringDeltaPhi + "entries", {HistType::kTH1D, {axisSpecDeltaPhi}}, true);

    if (applyEfficiency && useCcdbEfficiency) {
      ccdbApi.init(ccdbUrl);
      std::map<std::string, std::string> const metadata;
      bool const isEfficiencyDstarfileAvailable = ccdbApi.retrieveBlob(ccdbPathEfficiencyDstar, ".", metadata, ccdbTimestamp, false, efficiencyDstarFileName);
      if (!isEfficiencyDstarfileAvailable) {
        LOGF(fatal, "Failed to retrieve efficiency file for Dstar from CCDB");
      }
      bool const isEfficiencyTracksfileAvailable = ccdbApi.retrieveBlob(ccdbPathEfficiencyTracks, ".", metadata, ccdbTimestamp, false, efficiencyTracksFileName);
      if (!isEfficiencyTracksfileAvailable) {
        LOGF(fatal, "Failed to retrieve efficiency file for tracks from CCDB");
      }

      TFile* efficiencyDstarRootFile = TFile::Open(efficiencyDstarFileName.value.c_str(), "READ");
      if (!efficiencyDstarRootFile || efficiencyDstarRootFile->IsZombie()) {
        LOGF(fatal, "Failed to open efficiency file for Dstar");
      }

      TFile* efficiencyTracksRootFile = TFile::Open(efficiencyTracksFileName.value.c_str(), "READ");
      if (!efficiencyTracksRootFile || efficiencyTracksRootFile->IsZombie()) {
        LOGF(fatal, "Failed to open efficiency file for tracks");
      }

      vecHistEfficiencyDstar.resize(nEfficiencyHist);
      vecHistEfficiencyTracks.resize(nEfficiencyHist);

      for (int iHist = 0; iHist < nEfficiencyHist; iHist++) {
        vecHistEfficiencyDstar[iHist] = dynamic_cast<TH1F*>(efficiencyDstarRootFile->Get(Form("hEfficiencyDstar_%d", iHist)));
        if (!vecHistEfficiencyDstar[iHist]) {
          LOGF(fatal, "Failed to retrieve Dstar efficiency histogram hEfficiencyDstar_%d from file", iHist);
        }

        vecHistEfficiencyTracks[iHist] = dynamic_cast<TH1D*>(efficiencyTracksRootFile->Get(Form("hEfficiencyTracks_%d", iHist)));
        if (!vecHistEfficiencyTracks[iHist]) {
          LOGF(fatal, "Failed to retrieve track efficiency histogram hEfficiencyTracks_%d from file", iHist);
        }
        vecHistEfficiencyDstar[iHist]->SetDirectory(nullptr);
        vecHistEfficiencyTracks[iHist]->SetDirectory(nullptr);
      }

      efficiencyDstarRootFile->Close();
      efficiencyTracksRootFile->Close();
      delete efficiencyDstarRootFile;
      delete efficiencyTracksRootFile;
    }
  }

  void processData(aod::DstarHadronPair const& dstarHPairs)
  {
    for (const auto& dstarHPair : dstarHPairs) {
      float const deltaPhi = dstarHPair.deltaPhi();
      float const deltaEta = dstarHPair.deltaEta();
      float const ptDstar = dstarHPair.ptDstar();
      float const ptTrack = dstarHPair.ptTrack();
      int const poolBin = dstarHPair.poolBin();
      float const deltaM = dstarHPair.deltaM();

      int const effBinPtDstar = o2::analysis::findBin(binsPtEfficiency, ptDstar);
      // LOG(info) << "efficiency index " << effBinPtDstar;
      int const corrBinPtDstar = o2::analysis::findBin(binsPtCorrelations, ptDstar);
      // LOG(info) << "correlation index " << corrBinPtDstar;

      int const effBinPtTrack = o2::analysis::findBin(binsPtEfficiencyTracks, ptTrack);
      // LOG(info) << "track efficiency index " << effBinPtTrack;

      // reject candidate if outside pT ranges of interst
      if (corrBinPtDstar < 0 || effBinPtDstar < 0) {
        continue;
      }
      // Why the follwing condition in Dplus task?
      // if (ptTrack > 10.0) {
      //   ptTrack = 10.5;
      // }
      float netEfficiencyWeight = 1.0;

      if (applyEfficiency && !useCcdbEfficiency) {
        float const efficiencyWeightDstar = efficiencyDstar->at(effBinPtDstar);
        // LOG(info)<<"efficiencyWeightDstar "<<efficiencyWeightDstar;
        float const efficiencyWeightTracks = efficiencyTracks->at(effBinPtTrack);
        // LOG(info)<<"efficiencyWeightTracks "<<efficiencyWeightTracks;
        netEfficiencyWeight = 1.0 / (efficiencyWeightDstar * efficiencyWeightTracks);
      } else if (applyEfficiency && useCcdbEfficiency && nEfficiencyHist == 1) {
        float const efficiencyWeightDstar = vecHistEfficiencyDstar[0]->GetBinContent(vecHistEfficiencyDstar[0]->GetXaxis()->FindBin(ptDstar));
        // LOG(info)<<"efficiencyWeightDstar "<<efficiencyWeightDstar;
        float const efficiencyWeightTracks = vecHistEfficiencyTracks[0]->GetBinContent(vecHistEfficiencyTracks[0]->GetXaxis()->FindBin(ptTrack));
        // LOG(info)<<"efficiencyWeightTracks "<<efficiencyWeightTracks;
        netEfficiencyWeight = 1.0 / (efficiencyWeightDstar * efficiencyWeightTracks);
      } else if (applyEfficiency && useCcdbEfficiency && nEfficiencyHist > 1) {
        // to do
        LOGF(fatal, "Using CCDB efficiency with more than 1 histogram is not implemented yet");
      }

      // check if correlation entry belongs to signal region, sidebands or is outside both, and fill correlation plots
      if (deltaM > signalRegionLefBound->at(corrBinPtDstar) && deltaM < signalRegionRightBound->at(corrBinPtDstar)) {
        // in signal region
        registry.fill(HIST("hCorrel2DVsPtSignalRegion"), deltaPhi, deltaEta, ptDstar, ptTrack, poolBin, netEfficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSignalRegion"), deltaPhi, deltaEta, netEfficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSignalRegion"), deltaEta, netEfficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSignalRegion"), deltaPhi, netEfficiencyWeight);
      } else if (/*(deltaM > leftSidebandOuterBoundary->at(corrBinPtDstar) && deltaM < leftSidebandInnerBoundary->at(corrBinPtDstar)) ||*/ (deltaM > rightSidebandInnerBoundary->at(corrBinPtDstar) && deltaM < rightSidebandOuterBoundary->at(corrBinPtDstar))) {
        registry.fill(HIST("hCorrel2DVsPtSidebands"), deltaPhi, deltaEta, ptDstar, ptTrack, poolBin, netEfficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSidebands"), deltaPhi, deltaEta, netEfficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSidebands"), deltaEta, netEfficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSidebands"), deltaPhi, netEfficiencyWeight);
      }
    }
  }
  PROCESS_SWITCH(HfTaskCorrelationDstarHadrons, processData, " process data only", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskCorrelationDstarHadrons>(cfgc)};
}
