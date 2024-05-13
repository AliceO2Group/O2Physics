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
/// \author Jasper Parkkila (jparkkil@cern.ch)
/// \author Dong Jo Kim (djkim@jyu.fi)
/// \since Sep 2022

#include <TFile.h>
#include <TH1.h>
#include <TTree.h>
#include <deque>
#include <memory>

#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/HistogramRegistry.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/Centrality.h"
#include "ReconstructionDataFormats/V0.h"

// #include "CCDB/BasicCCDBManager.h"

#include "PWGCF/JCorran/DataModel/JCatalyst.h"
#include "PWGCF/DataModel/CorrelationsDerived.h"
#include "JFFlucAnalysis.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

void customize(std::vector<o2::framework::ConfigParamSpec>& workflowOptions)
{
  ConfigParamSpec optionLoadWeights{"loadWeights", VariantType::Bool, false, {"Load correction weights"}};
  workflowOptions.push_back(optionLoadWeights);
}

#include "Framework/runDataProcessing.h"

#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};

namespace o2::aod
{
namespace jweight
{
DECLARE_SOA_COLUMN(WeightNUA, weightNUA, float); //! Non-uniform acceptance weight
DECLARE_SOA_COLUMN(WeightEff, weightEff, float); //! Non-uniform efficiency weight
} // namespace jweight
DECLARE_SOA_TABLE(JWeights, "AOD", "JWEIGHT", jweight::WeightNUA, jweight::WeightEff); //! JFluc table for weights
} // namespace o2::aod

// The standalone jfluc code expects the entire list of tracks for an event. At the same time, it expects weights together with other track attributes.
// This workflow creates a table of weights that can be joined with track tables.
struct jflucWeightsLoader {
  O2_DEFINE_CONFIGURABLE(pathPhiWeights, std::string, "", "Local (local://) or CCDB path for the phi acceptance correction histogram");

  UInt_t ncent;
  float* pcentEdges = 0;

  struct Map {
    Map(int _runNumber, UInt_t _centBin, TH1* _ph) : runNumber(_runNumber), centBin(_centBin), ph(_ph) {}
    ~Map() { delete ph; }
    int runNumber;
    UInt_t centBin;
    TH1* ph;
  };
  std::deque<Map> nuaCache;
  TFile* pf = 0;

  ~jflucWeightsLoader()
  {
    if (pcentEdges)
      delete[] pcentEdges;
    if (pf) {
      pf->Close();
      delete pf;
    }
  }

  Produces<aod::JWeights> output;
  void init(InitContext const&)
  {
    //
    if (!doprocessLoadWeights)
      return;
    if (pathPhiWeights.value.substr(0, 8) == "local://") {
      pf = new TFile(pathPhiWeights.value.substr(8).c_str(), "read");

      // TODO: who knows someday, replace the old collection of TH3Ds with a 4D axis that includes the centrality info.
      // At the same time, the centBin indexing needs to go, use actual bin edges instead.
      TTree* pt = static_cast<TTree*>(pf->Get("axes"));
      if (pt)
        LOGF(fatal, "The NUA correction root file does not have the axis information.");
      pt->SetBranchAddress("Ncentrality", &ncent);
      pt->GetEntry(0);

      pcentEdges = new float[ncent];
      pt->SetBranchAddress("centrality", pcentEdges);
      pt->GetEntry(0);
    }
  }

  void processLoadWeights(aod::JCollision const& collision, aod::JTracks const& tracks)
  {
    UInt_t centBin = 0;
    for (UInt_t i = 0, n = ncent - 1; i < n; ++i)
      if (collision.multiplicity() < pcentEdges[i + 1]) {
        centBin = i;
        break;
      }

    for (auto& track : tracks) {
      if (pf) {
        float phiWeight, effWeight;
        auto m = std::find_if(nuaCache.begin(), nuaCache.end(), [&](auto& t) -> bool {
          return t.runNumber == collision.runNumber() && t.centBin == centBin;
        });
        if (m == nuaCache.end()) {
          TH1* ph = static_cast<TH1*>(pf->Get(Form("PhiWeights_%u_%02u", collision.runNumber(), centBin)));
          if (ph) {
            ph->SetDirectory(0); // we delete when we please
            nuaCache.emplace_back(collision.runNumber(), centBin, ph);
            if (nuaCache.size() > ncent * 3)
              nuaCache.pop_front(); // keep at most maps for 3 runs
            Int_t bin = ph->FindBin(track.phi(), track.eta(), collision.posZ());
            phiWeight = ph->GetBinContent(bin);
          } else {
            phiWeight = 1.0f;
          }
        } else {
          Int_t bin = m->ph->FindBin(track.phi(), track.eta(), collision.posZ());
          phiWeight = m->ph->GetBinContent(bin);
        }

        effWeight = 1.0f; //<--- todo

        output(phiWeight, effWeight);
      }
    }
  }
  PROCESS_SWITCH(jflucWeightsLoader, processLoadWeights, "Load weights histograms", false);
};

struct jflucAnalysisTask {
  ~jflucAnalysisTask()
  {
    delete pcf;
  }

  O2_DEFINE_CONFIGURABLE(etamin, float, 0.4, "Minimal eta for tracks");
  O2_DEFINE_CONFIGURABLE(etamax, float, 0.8, "Maximal eta for tracks");
  O2_DEFINE_CONFIGURABLE(ptmin, float, 0.2, "Minimal pt for tracks");
  O2_DEFINE_CONFIGURABLE(ptmax, float, 0.5, "Maximal pt for tracks");

  ConfigurableAxis axisMultiplicity{"axisMultiplicity", {VARIABLE_WIDTH, 0, 5, 10, 20, 30, 40, 50, 100.1}, "multiplicity / centrality axis for histograms"};

  Filter jtrackFilter = (aod::jtrack::pt > ptmin) && (aod::jtrack::pt < ptmax);    // eta cuts done by jfluc
  Filter cftrackFilter = (aod::cftrack::pt > ptmin) && (aod::cftrack::pt < ptmax); // eta cuts done by jfluc

  OutputObj<TDirectory> output{"jflucO2"};

  void init(InitContext const&)
  {
    pcf = new JFFlucAnalysis("jflucAnalysis");
    pcf->SetNumBins(AxisSpec(axisMultiplicity).getNbins());
    pcf->AddFlags(JFFlucAnalysis::kFlucEbEWeighting);

    output->cd();
    pcf->UserCreateOutputObjects();
  }

  template <class CollisionT, class TrackT>
  void analyze(CollisionT const& collision, TrackT const& tracks)
  {
    pcf->Init();
    pcf->FillQA(tracks);
    qvecs.Calculate(tracks, etamin, etamax);
    pcf->SetJQVectors(&qvecs);
    const auto& edges = AxisSpec(axisMultiplicity).binEdges;
    for (UInt_t i = 0, n = AxisSpec(axisMultiplicity).getNbins(); i < n; ++i)
      if (collision.multiplicity() < edges[i + 1]) {
        pcf->SetEventCentralityAndBin(collision.multiplicity(), i);
        break;
      }
    const double fVertex[3] = {0.0f, 0.0f, collision.posZ()}; // TODO: check if posX/Y is really needed
    pcf->SetEventVertex(fVertex);
    pcf->SetEtaRange(etamin, etamax);
    pcf->UserExec("");
  }

  void process(aod::JCollision const& collision, soa::Filtered<aod::JTracks> const& tracks)
  {
    analyze(collision, tracks);
  }
  PROCESS_SWITCH(jflucAnalysisTask, process, "Process data", true);

  void processCorrected(aod::JCollision const& collision, soa::Filtered<soa::Join<aod::JTracks, aod::JWeights>> const& tracks)
  {
    analyze(collision, tracks);
  }
  PROCESS_SWITCH(jflucAnalysisTask, processCorrected, "Process data with corrections", false);

  void processCFDerived(aod::CFCollision const& collision, soa::Filtered<aod::CFTracks> const& tracks)
  {
    analyze(collision, tracks);
  }
  PROCESS_SWITCH(jflucAnalysisTask, processCFDerived, "Process CF derived data", false);

  void processCFDerivedCorrected(aod::CFCollision const& collision, soa::Filtered<soa::Join<aod::CFTracks, aod::JWeights>> const& tracks)
  {
    analyze(collision, tracks);
  }
  PROCESS_SWITCH(jflucAnalysisTask, processCFDerivedCorrected, "Process CF derived data with corrections", false);

  JFFlucAnalysis::JQVectorsT qvecs;
  JFFlucAnalysis* pcf;
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  if (cfgc.options().get<bool>("loadWeights")) {
    return WorkflowSpec{
      adaptAnalysisTask<jflucWeightsLoader>(cfgc),
      adaptAnalysisTask<jflucAnalysisTask>(cfgc)};
  } else {
    return WorkflowSpec{adaptAnalysisTask<jflucAnalysisTask>(cfgc)};
  }
}
