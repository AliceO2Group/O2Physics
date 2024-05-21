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
#include "JFFlucAnalysisO2Hist.h"

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

  struct Map {
    Map(THnF* _ph, int _runNumber) : ph(_ph), runNumber(_runNumber) {}
    ~Map() { delete ph; }
    THnF* ph;
    int runNumber;
  };
  std::deque<Map> nuaCache;
  TFile* pf = 0;

  ~jflucWeightsLoader()
  {
    if (pf) {
      nuaCache.clear();
      pf->Close();
      delete pf;
    }
  }

  Produces<aod::JWeights> output;
  void init(InitContext const&)
  {
    if (!doprocessLoadWeights && !doprocessLoadWeightsCF)
      return;
    if (doprocessLoadWeights && doprocessLoadWeightsCF)
      LOGF(fatal, "Only one weights loader process switch can be enabled at a time.");
    if (pathPhiWeights.value.substr(0, 8) == "local://") {
      pf = new TFile(pathPhiWeights.value.substr(8).c_str(), "read");
      if (!pf->IsOpen()) {
        delete pf;
        LOGF(fatal, "NUA correction weights file not found: %s", pathPhiWeights.value.substr(8).c_str());
      }
    }
  }

  template <class CollisionT, class TrackT>
  void loadWeights(CollisionT const& collision, TrackT const& tracks)
  {
    if (!pf)
      LOGF(fatal, "NUA correction weights file has not been opened.");
    for (auto& track : tracks) {
      float phiWeight, effWeight;
      auto m = std::find_if(nuaCache.begin(), nuaCache.end(), [&](auto& t) -> bool {
        return t.runNumber == collision.runNumber();
      });
      if (m == nuaCache.end()) {
        THnF* ph = static_cast<THnF*>(pf->Get(Form("NUAWeights_%u", collision.runNumber())));
        if (ph) {
          nuaCache.emplace_back(ph, collision.runNumber());
          if (nuaCache.size() > 3)
            nuaCache.pop_front(); // keep at most maps for 3 runs
          const Double_t coords[] = {collision.multiplicity(), track.phi(), track.eta(), collision.posZ()};
          auto bin = ph->GetBin(coords);
          phiWeight = ph->GetBinContent(bin);
        } else {
          phiWeight = 1.0f;
          // LOGF(error, "NUAWeights_%u not found", collision.runNumber());
        }
      } else {
        const Double_t coords[] = {collision.multiplicity(), track.phi(), track.eta(), collision.posZ()};
        auto bin = m->ph->GetBin(coords);
        phiWeight = m->ph->GetBinContent(bin);
      }

      effWeight = 1.0f; //<--- todo

      output(phiWeight, effWeight);
    }
  }

  void processLoadWeights(aod::JCollision const& collision, aod::JTracks const& tracks)
  {
    loadWeights(collision, tracks);
  }
  PROCESS_SWITCH(jflucWeightsLoader, processLoadWeights, "Load weights histograms for derived data table", false);

  void processLoadWeightsCF(aod::CFCollision const& collision, aod::CFTracks const& tracks)
  {
    loadWeights(collision, tracks);
  }
  PROCESS_SWITCH(jflucWeightsLoader, processLoadWeightsCF, "Load weights histograms for CF derived data table", true);
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

  HistogramRegistry registry{"registry"};
  // OutputObj<JFFlucAnalysis> output{JFFlucAnalysis("jflucO2")};

  void init(InitContext const&)
  {
    auto a = AxisSpec(axisMultiplicity);
    pcf = new JFFlucAnalysisO2Hist(registry, a);
    pcf->AddFlags(JFFlucAnalysis::kFlucEbEWeighting);
    pcf->UserCreateOutputObjects();
  }

  template <class CollisionT, class TrackT>
  void analyze(CollisionT const& collision, TrackT const& tracks)
  {
    pcf->Init();
    pcf->SetEventCentrality(collision.multiplicity());
    const double fVertex[3] = {0.0f, 0.0f, collision.posZ()}; // TODO: check if posX/Y is really needed
    pcf->SetEventVertex(fVertex);
    pcf->SetEtaRange(etamin, etamax);
    pcf->FillQA(tracks);
    qvecs.Calculate(tracks, etamin, etamax);
    pcf->SetJQVectors(&qvecs);
    pcf->UserExec("");
  }

  void processJDerived(aod::JCollision const& collision, soa::Filtered<aod::JTracks> const& tracks)
  {
    analyze(collision, tracks);
  }
  PROCESS_SWITCH(jflucAnalysisTask, processJDerived, "Process derived data", false);

  void processJDerivedCorrected(aod::JCollision const& collision, soa::Filtered<soa::Join<aod::JTracks, aod::JWeights>> const& tracks)
  {
    analyze(collision, tracks);
  }
  PROCESS_SWITCH(jflucAnalysisTask, processJDerivedCorrected, "Process derived data with corrections", false);

  void processCFDerived(aod::CFCollision const& collision, soa::Filtered<aod::CFTracks> const& tracks)
  {
    analyze(collision, tracks);
  }
  PROCESS_SWITCH(jflucAnalysisTask, processCFDerived, "Process CF derived data", true);

  void processCFDerivedCorrected(aod::CFCollision const& collision, soa::Filtered<soa::Join<aod::CFTracks, aod::JWeights>> const& tracks)
  {
    analyze(collision, tracks);
  }
  PROCESS_SWITCH(jflucAnalysisTask, processCFDerivedCorrected, "Process CF derived data with corrections", false);

  JFFlucAnalysis::JQVectorsT qvecs;
  JFFlucAnalysisO2Hist* pcf;
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
