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
// Authors: Lars Jörgensen
// Date: 29.05.2024

#include <cmath>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TObjArray.h>

#include "ReconstructionDataFormats/Track.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Framework/runDataProcessing.h"
#include "Framework/HistogramRegistry.h"
#include "PWGDQ/DataModel/ReducedInfoTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct AngularCorrelationsInJets {

	// HistogramRegistry registryName{"folderTitle", {}, OutputObjHandlingPolicy::AnalysisObject, <sortHistograms:bool>, <createDir:bool>};
  HistogramRegistry jetOutputReg{"jetOutput", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry nucleiOutputReg{"nucleiOutput", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry correlationsOutputReg{"correlationsOutput", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
	// OutputObj<Type> histogramName{Type("histogramName", "histogramTitle;Axis", nbins,minbin,maxbin)};

  void init(o2::framework::InitContext&)
  {

    std::vector<double> ptBinning = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.8, 2.0, 2.2, 2.4, 2.8, 3.2, 3.6, 4., 5., 6., 8., 10., 12., 14.};

		// AxisSpec specName = {binningInfo, "axisLabel"};
    AxisSpec ptAxis = {ptBinning, "#it{p}_{T} (GeV/#it{c})"};

		// registryName.add("histogramName", "histogramTitle", HistType::Type, {{binningInfo}});
		
    jetOutputReg.add("hPtFullEvent", "p_{T} after basic cuts", HistType::kTH1F, {ptAxis});
		jetOutputReg.add("hNumPartInEvent", "Number of particles per event", HistType::kTH1I, {{500,0,500}});
  }

	// Configurable<Type> cfgName{"nameOnHyperloop", value, "Flag shown on hyperloop"};
	Configurable<float> ptMin{"ptMin", 0.5, "Minimum pT for particles"}
  Configurable<float> yMin_gen{"yMin_gen", -0.5, "Maximum rapidity (generated)"};
  Configurable<float> yMax_gen{"yMax_gen", 0.5, "Minimum rapidity (generated)"};
  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.9f, "Eta range for tracks"};
  Configurable<float> yMin_reco{"yMin_reco", -0.5, "Maximum rapidity (reconstructed)"};
  Configurable<float> yMax_reco{"yMax_reco", 0.5, "Minimum rapidity (reconstructed)"};
  Configurable<float> pTmin_reco{"pTmin_reco", 0.1f, "min pT (reconstructed)"};
  Configurable<float> pTmax_reco{"pTmax_reco", 1e+10f, "max pT (reconstructed)"};
  Configurable<float> minReqClusterITS{"minReqClusterITS", 1.0, "min number of clusters required in ITS"};
  Configurable<float> minReqClusterITSib{"minReqClusterITSib", 1.0, "min number of clusters required in ITS inner barrel"};
  Configurable<float> minTPCnClsFound{"minTPCnClsFound", 0.0f, "min number of crossed rows TPC"};
  Configurable<float> minNCrossedRowsTPC{"minNCrossedRowsTPC", 70.0f, "min number of crossed rows TPC"};
  Configurable<float> minRatioCrossedRowsTPC{"minRatioCrossedRowsTPC", 0.8f, "min ratio of crossed rows over findable clusters TPC"};
  Configurable<float> maxRatioCrossedRowsTPC{"maxRatioCrossedRowsTPC", 1.5f, "max ratio of crossed rows over findable clusters TPC"};
  Configurable<float> maxChi2ITS{"maxChi2ITS", 36.0f, "max chi2 per cluster ITS"};
  Configurable<float> maxChi2TPC{"maxChi2TPC", 4.0f, "max chi2 per cluster TPC"};
  Configurable<float> maxDCA_XY{"maxDCA_XY", 0.5f, "max DCA to vertex xy"};
  Configurable<float> maxDCA_Z{"maxDCA_Z", 2.0f, "max DCA to vertex z"};

  //****************************************************************************************************

  template <typename CollisionType, typename TracksType>
  void processJetReco(const CollisionType& collision, const TracksType& tracks)
  {
    jetOutputReg.fill(HIST("hNumPartInEvent"), collision.);

    for (const auto& track : tracks) {
      if (track.pt()<ptMin)
        continue;

      jetOutputReg.fill(HIST("hPtFullEvent"), track.pt());
    }
  }

  //****************************************************************************************************

  void process_jet_eco(aod::Collision const& collision, soa::Join<aod::Tracks,aod::TracksExtra> const& tracks)
  {
    processJetReco(collision, tracks);
  }
  PROCESS_SWITCH(AngularCorrelationsInJets, processJetReco, "reconstruct jets", true);
};

//****************************************************************************************************

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<AngularCorrelationsInJets>(cfgc, TaskName{"angular-correlations-in-jets"})};
}
