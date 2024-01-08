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
// This task is designed to do QA to the TOF PID applied to strangeness
// in the regular framework

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/LFStrangenessPIDTables.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/PIDResponse.h"

#include <TFile.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TLorentzVector.h>
#include <Math/Vector4D.h>
#include <TPDGCode.h>
#include <TDatabasePDG.h>
#include <cmath>
#include <array>
#include <cstdlib>
#include "Framework/ASoAHelpers.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;
using std::cout;
using std::endl;

struct strangepidqa {
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  ConfigurableAxis vertexZ{"vertexZ", {30, -15.0f, 15.0f}, ""};

  // base properties
  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 17.0f, 19.0f, 21.0f, 23.0f, 25.0f, 30.0f, 35.0f, 40.0f, 50.0f}, "p_{T} (GeV/c)"};
  ConfigurableAxis axisRadius{"axisRadius", {200, 0.0f, 100.0f}, "V0 radius (cm)"};

  // Invariant Mass
  ConfigurableAxis axisLambdaMass{"axisLambdaMass", {200, 1.08f, 1.16f}, "M_{#Lambda} (GeV/c^{2})"};

  // Delta time
  ConfigurableAxis axisDeltaTime{"axisDeltaTime", {200, -1000.0f, +1000.0f}, "#Delta time (ps)"};

  // Length axis
  ConfigurableAxis axisLength{"axisLength", {600, 0.0f, +600.0f}, "track Length (cm)"};

  void init(InitContext const&)
  {
    // Event counter
    histos.add("hEventVertexZ", "hEventVertexZ", kTH1F, {vertexZ});

    // V0 Radius
    histos.add("h2dLambdaRadiusVsPt", "hLambdaRadiusVsPt", {HistType::kTH2F, {axisPt, axisRadius}});

    // Invariant Mass
    histos.add("h2dLambdaMassVsPt", "hLambdaMassVsPt", {HistType::kTH2F, {axisPt, axisLambdaMass}});

    // radius vs prong length
    histos.add("h2dProtonLengthVsRadius", "h2dProtonLengthVsRadius", {HistType::kTH2F, {axisRadius, axisLength}});
    histos.add("h2dPionLengthVsRadius", "h2dPionLengthVsRadius", {HistType::kTH2F, {axisRadius, axisLength}});

    // recalculated vs topv lengths
    histos.add("h2dProtonLengthVsLengthToPV", "h2dProtonLengthVsRadiusToPV", {HistType::kTH2F, {axisRadius, axisRadius}});
    histos.add("h2dPionLengthVsLengthToPV", "h2dPionLengthVsLengthToPV", {HistType::kTH2F, {axisRadius, axisRadius}});

    // TOF PID testing for prongs
    histos.add("h2dProtonDeltaTimeVsPt", "h2dProtonDeltaTimeVsPt", {HistType::kTH2F, {axisPt, axisDeltaTime}});
    histos.add("h2dProtonDeltaTimeVsRadius", "h2dProtonDeltaTimeVsRadius", {HistType::kTH2F, {axisRadius, axisDeltaTime}});
    histos.add("h2dPionDeltaTimeVsPt", "h2dPionDeltaTimeVsPt", {HistType::kTH2F, {axisPt, axisDeltaTime}});
    histos.add("h2dPionDeltaTimeVsRadius", "h2dPionDeltaTimeVsRadius", {HistType::kTH2F, {axisRadius, axisDeltaTime}});
  }

  void process(aod::StraCollision const& coll, soa::Join<aod::V0Cores, aod::V0CollRefs, aod::V0Extras, aod::V0MCDatas, aod::V0TOFs, aod::V0TOFPIDs> const& v0s)
  {
    histos.fill(HIST("hEventVertexZ"), coll.posZ());

    for (auto& lambda : v0s) { // selecting photons from Sigma0
      if (lambda.pdgCode() != 3122)
        continue;

      histos.fill(HIST("h2dLambdaRadiusVsPt"), lambda.pt(), lambda.v0radius());
      histos.fill(HIST("h2dLambdaMassVsPt"), lambda.pt(), lambda.mLambda());

      histos.fill(HIST("h2dProtonDeltaTimeVsPt"), lambda.pt(), lambda.posTOFDeltaTLaPr());
      histos.fill(HIST("h2dProtonDeltaTimeVsRadius"), lambda.v0radius(), lambda.posTOFDeltaTLaPr());
      histos.fill(HIST("h2dPionDeltaTimeVsPt"), lambda.pt(), lambda.negTOFDeltaTLaPi());
      histos.fill(HIST("h2dPionDeltaTimeVsRadius"), lambda.v0radius(), lambda.negTOFDeltaTLaPi());

      histos.fill(HIST("h2dProtonLengthVsRadius"), lambda.v0radius(), lambda.posTOFLength());
      histos.fill(HIST("h2dPionLengthVsRadius"), lambda.v0radius(), lambda.negTOFLength());

      histos.fill(HIST("h2dProtonLengthVsLengthToPV"), lambda.posTOFLengthToPV(), lambda.posTOFLength());
      histos.fill(HIST("h2dPionLengthVsLengthToPV"), lambda.negTOFLengthToPV(), lambda.negTOFLength());
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<strangepidqa>(cfgc)};
}
