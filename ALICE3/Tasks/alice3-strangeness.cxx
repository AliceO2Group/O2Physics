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
/// \file alice3-strangeness.cxx
///
/// \brief This task produces invariant mass distributions for strange hadrons
///
/// \author Lucia Anna Tarasovičová, Pavol Jozef Šafárik University (SK)
/// \since  November 20, 2025
///

#include "ALICE3/DataModel/OTFStrangeness.h"
#include "ALICE3/DataModel/tracksAlice3.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <DCAFitter/DCAFitterN.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <DetectorsBase/Propagator.h>
#include <DetectorsVertexing/PVertexer.h>
#include <DetectorsVertexing/PVertexerHelpers.h>
#include <Field/MagneticField.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/runDataProcessing.h>
#include <ReconstructionDataFormats/DCA.h>
#include <SimulationDataFormat/InteractionSampler.h>

#include <TGenPhaseSpace.h>
#include <TGeoGlobalMagField.h>
#include <TLorentzVector.h>
#include <TPDGCode.h>
#include <TRandom3.h>

#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;

using alice3tracks = soa::Join<aod::Tracks, aod::TracksCov, aod::McTrackLabels, aod::TracksDCA, aod::TracksExtraA3>;

struct alice3strangeness {
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 17.0f, 19.0f, 21.0f, 23.0f, 25.0f, 30.0f, 35.0f, 40.0f, 50.0f}, "pt axis for QA histograms"};
  ConfigurableAxis axisK0Mass{"axisK0Mass", {200, 0.4f, 0.6f}, ""};
  ConfigurableAxis axisVertexZ{"axisVertexZ", {40, -20, 20}, "vertex Z (cm)"};

  void init(InitContext&)
  {
    histos.add("K0/hMassAllCandidates", "", kTH2D, {axisK0Mass, axisPt});
    histos.add("K0/hMassSelected", "", kTH2D, {axisK0Mass, axisPt});
    histos.add("K0/hSelections", "", kTH1D, {{10, 0, 10}});
    histos.add("K0/hDCANegDaughter", "", kTH1D, {{200, -5, 5}});
    histos.add("K0/hDCAPosDaughter", "", kTH1D, {{200, -5, 5}});
    histos.add("hPVz", "hPVz", kTH1F, {axisVertexZ});
  }
  long int nEvents = 0;
  void process(aod::Collisions const& collisions, aod::McCollisions const& mcCollisions, aod::UpgradeV0s const& v0Recos, alice3tracks const&)
  {
    LOG(info) << "Event processed " << nEvents++ << " :" << collisions.size() << " " << mcCollisions.size();
    for (const auto& collision : collisions) {
      float collisionZ = collision.posZ();
      // std::cout << "______ process V0_______" <<  collision.size() << std::endl;
      histos.fill(HIST("hPVz"), collisionZ);
      for (const auto& v0Cand : v0Recos) {

        auto negV0Daughter = v0Cand.negTrack_as<alice3tracks>(); // de-reference neg track
        auto posV0Daughter = v0Cand.posTrack_as<alice3tracks>(); // de-reference pos track

        bool isK0 = v0Cand.mK0() > 0;
        if (isK0) {
          histos.fill(HIST("K0/hMassAllCandidates"), v0Cand.mK0(), v0Cand.pt());
          histos.fill(HIST("K0/hSelections"), 0); // all candidates
          histos.fill(HIST("K0/hDCANegDaughter"), negV0Daughter.dcaXY());
          histos.fill(HIST("K0/hDCAPosDaughter"), posV0Daughter.dcaXY());
          if (std::abs(negV0Daughter.dcaXY()) < 0.05)
            continue;
          histos.fill(HIST("K0/hSelections"), 1); // dcaXY cut
          if (std::abs(posV0Daughter.dcaXY()) < 0.05)
            continue;
          histos.fill(HIST("K0/hSelections"), 2); // dcaXY cut
          if (v0Cand.dcaV0Daughters() > 1.0)
            continue;
          histos.fill(HIST("K0/hSelections"), 3); // dca between daughters
          if (v0Cand.v0Radius() < 0.5)
            continue;
          histos.fill(HIST("K0/hSelections"), 4); // radius cut
          if (std::abs(negV0Daughter.eta()) > 0.8 || std::abs(posV0Daughter.eta()) > 0.8)
            continue;
          histos.fill(HIST("K0/hSelections"), 5); // eta cut
          histos.fill(HIST("K0/hMassSelected"), v0Cand.mK0(), v0Cand.pt());
        }
      }
    }
  }
  PROCESS_SWITCH(alice3strangeness, process, "", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<alice3strangeness>(cfgc)};
}
