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

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include <CCDB/BasicCCDBManager.h>
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/HistogramRegistry.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/Centrality.h"
#include "ReconstructionDataFormats/V0.h"

#include "PWGCF/JCorran/DataModel/JCatalyst.h"
#include "JFFlucAnalysis.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};

namespace o2::aod
{
namespace jweight
{
DECLARE_SOA_COLUMN(WeightNUA, weightNUA, float); //! Non-uniform acceptance weight
DECLARE_SOA_COLUMN(WeightEff, weightEff, float); //! Non-uniform efficiency weight
} // namespace cfmultiplicity
DECLARE_SOA_TABLE(JWeights, "AOD", "JWEIGHT", jweight::WeightNUA, jweight::WeightEff); //! JFluc table for weights
} // namespace o2::aod

struct jflucWeightsLoader {
  //
	Produces<aod::JWeights> output;
  void processLoadWeights(aod::JCollision const& collision, aod::JTracks const& tracks){
  	//In the future, this will do the job of weight loading
	for(auto &track : tracks){
		(void)track;
		output(1.0f,1.0f);
	}
	}
    PROCESS_SWITCH(jflucWeightsLoader, processLoadWeights, "Load weights histograms", true);
};

struct jflucAnalysisTask {
  ~jflucAnalysisTask()
  {
    delete pcf;
  }

  O2_DEFINE_CONFIGURABLE(etamin, double, 0.4, "Minimal eta for tracks");
  O2_DEFINE_CONFIGURABLE(etamax, double, 0.8, "Maximal eta for tracks");

  OutputObj<TDirectory> output{"jflucO2"};

  void init(InitContext const& ic)
  {
    //
    pcf = new JFFlucAnalysis("jflucAnalysis");
    //pcf->SetNumBins(sizeof(jflucCentBins) / sizeof(jflucCentBins[0]));
    pcf->SetNumBins(5);
    pcf->AddFlags(JFFlucAnalysis::kFlucEbEWeighting);

    output->cd();
    pcf->UserCreateOutputObjects();
  }

  //void process(soa::Join<aod::Collisions, aod::CollisionData>::iterator const& collision, aod::ParticleTrack const& tracks)
  void process(aod::JCollision const& collision, soa::Join<aod::JTracks, aod::JWeights> const& tracks)
  {
    const double fVertex[3] = {0.0f, 0.0f, collision.posZ()}; //TODO: check if posX/Y is really needed
    pcf->Init();
    pcf->FillQA(tracks);
    pcf->CalculateQvectorsQC(tracks);
    //pcf->SetEventCentralityAndBin(collision.cent(), collision.cbin());
    pcf->SetEventVertex(fVertex);
    pcf->SetEtaRange(etamin, etamax);
    pcf->UserExec("");
  }
  JFFlucAnalysis* pcf;
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<jflucWeightsLoader>(cfgc),
    adaptAnalysisTask<jflucAnalysisTask>(cfgc)};
}
