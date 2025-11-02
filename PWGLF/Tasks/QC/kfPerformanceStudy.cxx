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
/// \brief this task allows for the direct one-to-one comparison of
//         cascades computed with standard DCAFitter methods and the KFparticle
//         package. It is meant for the purposes of larger-scale QA of KF reco.

#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

#include <cmath>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;

// allows for candidate-by-candidate comparison using Cascade to []CascData link table
using CascadesCrossLinked = soa::Join<aod::Cascades, aod::CascDataLink, aod::KFCascDataLink>;

struct kfPerformanceStudy {
  // configurable binning of histograms
  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 17.0f, 19.0f, 21.0f, 23.0f, 25.0f, 30.0f, 35.0f, 40.0f, 50.0f}, "pt axis for QA histograms"};
  ConfigurableAxis axisLambdaMass{"axisLambdaMass", {200, 1.101f, 1.131f}, ""};
  ConfigurableAxis axisXiMass{"axisXiMass", {200, 1.222f, 1.422f}, ""};
  ConfigurableAxis axisOmegaMass{"axisOmegaMass", {200, 1.572f, 1.772f}, ""};
  ConfigurableAxis axisDCAxy{"axisDCAxy", {500, -1.0f, 1.0f}, ""};
  ConfigurableAxis axisPointingAngle{"axisPointingAngle", {800, 0.0f, 3.5f}, ""};
  ConfigurableAxis axisVertex{"axisVertex", {1000, -3.0f, 3.0f}, ""};
  ConfigurableAxis axisRadius{"axisRadius", {1000, 0.0f, 3.0f}, ""};

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext const&)
  {
    histos.add("hEventCounter", "hEventCounter", kTH1F, {{1, 0.0f, 1.0f}});
    histos.add("hChargeCounter", "hChargeCounter", kTH1F, {{3, -1.5f, 1.5f}});

    // baseline simple histograms of mass
    histos.add("hMassXiMinus", "hMassXiMinus", kTH2F, {axisPt, axisXiMass});
    histos.add("hMassXiPlus", "hMassXiPlus", kTH2F, {axisPt, axisXiMass});
    histos.add("hMassOmegaMinus", "hMassOmegaMinus", kTH2F, {axisPt, axisOmegaMass});
    histos.add("hMassOmegaPlus", "hMassOmegaPlus", kTH2F, {axisPt, axisOmegaMass});
    histos.add("hKFMassXiMinus", "hKFMassXiMinus", kTH2F, {axisPt, axisXiMass});
    histos.add("hKFMassXiPlus", "hKFMassXiPlus", kTH2F, {axisPt, axisXiMass});
    histos.add("hKFMassOmegaMinus", "hKFMassOmegaMinus", kTH2F, {axisPt, axisOmegaMass});
    histos.add("hKFMassOmegaPlus", "hKFMassOmegaPlus", kTH2F, {axisPt, axisOmegaMass});

    histos.add("hPointingAngle", "hPointingAngle", kTH1F, {axisPointingAngle});
    histos.add("hKFPointingAngle", "hKFPointingAngle", kTH1F, {axisPointingAngle});

    histos.add("hVertexX", "hVertexX", kTH1F, {axisVertex});
    histos.add("hKFVertexX", "hKFVertexX", kTH1F, {axisVertex});
    histos.add("hVertexY", "hVertexY", kTH1F, {axisVertex});
    histos.add("hKFVertexY", "hKFVertexY", kTH1F, {axisVertex});
    histos.add("hVertexZ", "hVertexZ", kTH1F, {axisVertex});
    histos.add("hKFVertexZ", "hKFVertexZ", kTH1F, {axisVertex});

    histos.add("hDCAxy", "hDCAxy", kTH1F, {axisDCAxy});
    histos.add("hKFDCAxy", "hKFDCAxy", kTH1F, {axisDCAxy});

    histos.add("hCascRadius", "hCascRadius", kTH1F, {axisRadius});
    histos.add("hKFCascRadius", "hKFCascRadius", kTH1F, {axisRadius});

    histos.add("h3dMassXiMinus", "h3dMassXiMinus", kTH3F, {axisPt, axisXiMass, axisXiMass});
    histos.add("h3dMassXiPlus", "h3dMassXiPlus", kTH3F, {axisPt, axisXiMass, axisXiMass});
    histos.add("h3dMassOmegaMinus", "h3dMassOmegaMinus", kTH3F, {axisPt, axisOmegaMass, axisOmegaMass});
    histos.add("h3dMassOmegaPlus", "h3dMassOmegaPlus", kTH3F, {axisPt, axisOmegaMass, axisOmegaMass});

    // all-in histograms
    histos.add("h3dPointingAngle", "h3dPointingAngle", kTH3F, {axisPt, axisPointingAngle, axisPointingAngle});
    histos.add("h3dMassLambda", "h3dMassLambda", kTH3F, {axisPt, axisLambdaMass, axisLambdaMass}); /// for x check only
    histos.add("h3dDCAxy", "h3dDCAxy", kTH3F, {axisPt, axisDCAxy, axisDCAxy});
    histos.add("hPtCorrelation", "hPtCorrelation", kTH2F, {axisPt, axisPt});
  }

  void process(aod::Collision const& collision, CascadesCrossLinked const& Cascades, aod::CascDatas const&, aod::KFCascDatas const&, aod::TracksIU const&)
  {
    histos.fill(HIST("hEventCounter"), 0.5);
    for (auto& cascade : Cascades) { // allows for cross-referencing everything
      float pt = 0.0f, ptKF = 0.0f;
      float massXi = 0.0f, massXiKF = 0.0f;
      float massOmega = 0.0f, massOmegaKF = 0.0f;
      float massLambda = 0.0f, massLambdaKF = 0.0f;
      float dcaXY = 0.0f, dcaXYKF = 0.0f;
      float pointingAngle = -1.0f, pointingAngleKF = -1.0f;
      float vertexX = 0.0f, vertexY = 0.0f, vertexZ = 0.0f;
      float vertexXKF = 0.0f, vertexYKF = 0.0f, vertexZKF = 0.0f;
      float cascRadius = 0.0f, cascRadiusKF = 0.0f;
      int charge = 0;

      // get charge from bachelor (unambiguous wrt to building)
      auto bachTrack = cascade.bachelor_as<aod::TracksIU>();
      if (bachTrack.sign() < 0) {
        charge = -1;
      } else {
        charge = +1;
      }

      histos.fill(HIST("hChargeCounter"), charge);

      if (cascade.has_cascData()) {
        // check aod::Cascades -> aod::CascData link
        // if present: this candidate was accepted by default DCAfitter building
        auto cascdata = cascade.cascData();
        pt = cascdata.pt();
        massLambda = cascdata.mLambda();
        massXi = cascdata.mXi();
        massOmega = cascdata.mOmega();
        dcaXY = cascdata.dcaXYCascToPV();
        pointingAngle = TMath::ACos(cascdata.casccosPA(collision.posX(), collision.posY(), collision.posZ()));
        vertexX = cascdata.x();
        vertexY = cascdata.y();
        vertexZ = cascdata.z();
        cascRadius = cascdata.cascradius();
      }
      if (cascade.has_kfCascData()) {
        // check aod::Cascades -> aod::KFCascData link
        // if present: this candidate was accepted by KF building
        auto cascdata = cascade.kfCascData();
        ptKF = cascdata.pt();
        massLambdaKF = cascdata.mLambda();
        massXiKF = cascdata.mXi();
        massOmegaKF = cascdata.mOmega();
        dcaXYKF = cascdata.dcaXYCascToPV();
        pointingAngleKF = TMath::ACos(cascdata.casccosPA(collision.posX(), collision.posY(), collision.posZ()));
        vertexXKF = cascdata.x();
        vertexYKF = cascdata.y();
        vertexZKF = cascdata.z();
        cascRadiusKF = cascdata.cascradius();
      }

      histos.fill(HIST("hPtCorrelation"), pt, ptKF);
      histos.fill(HIST("h3dMassLambda"), pt, massLambda, massLambdaKF);          // <- implicit pT choice, beware
      histos.fill(HIST("h3dDCAxy"), pt, dcaXY, dcaXYKF);                         // <- implicit pT choice, beware
      histos.fill(HIST("h3dPointingAngle"), pt, pointingAngle, pointingAngleKF); // <- implicit pT choice, beware

      histos.fill(HIST("hPointingAngle"), pointingAngle);
      histos.fill(HIST("hKFPointingAngle"), pointingAngleKF);

      histos.fill(HIST("hVertexX"), vertexX);
      histos.fill(HIST("hKFVertexX"), vertexXKF);
      histos.fill(HIST("hVertexY"), vertexY);
      histos.fill(HIST("hKFVertexY"), vertexYKF);
      histos.fill(HIST("hVertexZ"), vertexZ);
      histos.fill(HIST("hKFVertexZ"), vertexZKF);

      histos.fill(HIST("hCascRadius"), cascRadius);
      histos.fill(HIST("hKFCascRadius"), cascRadiusKF);

      histos.fill(HIST("hDCAxy"), dcaXY);
      histos.fill(HIST("hKFDCAxy"), dcaXYKF);

      if (charge < 0) {
        histos.fill(HIST("hMassXiMinus"), pt, massXi);
        histos.fill(HIST("hMassOmegaMinus"), pt, massOmega);
        histos.fill(HIST("hKFMassXiMinus"), ptKF, massXiKF);
        histos.fill(HIST("hKFMassOmegaMinus"), ptKF, massOmegaKF);
        histos.fill(HIST("h3dMassXiMinus"), pt, massXi, massXiKF);          // <- implicit pT choice, beware
        histos.fill(HIST("h3dMassOmegaMinus"), pt, massOmega, massOmegaKF); // <- implicit pT choice, beware
      }
      if (charge > 0) {
        histos.fill(HIST("hMassXiPlus"), pt, massXi);
        histos.fill(HIST("hMassOmegaPlus"), pt, massOmega);
        histos.fill(HIST("hKFMassXiPlus"), ptKF, massXiKF);
        histos.fill(HIST("hKFMassOmegaPlus"), ptKF, massOmegaKF);
        histos.fill(HIST("h3dMassXiPlus"), pt, massXi, massXiKF);          // <- implicit pT choice, beware
        histos.fill(HIST("h3dMassOmegaPlus"), pt, massOmega, massOmegaKF); // <- implicit pT choice, beware
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<kfPerformanceStudy>(cfgc)};
}
