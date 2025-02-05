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
/// \since May 2024

#include <experimental/type_traits>
#include <string>
#include <TFile.h>
#include <THn.h>

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
#include "Framework/runDataProcessing.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};

// The standalone jfluc code expects the entire list of tracks for an event. At the same time, it expects weights together with other track attributes.
// This workflow creates a table of weights that can be joined with track tables.
struct jflucWeightsLoader {
  O2_DEFINE_CONFIGURABLE(pathPhiWeights, std::string, "", "Local (local://) or CCDB path for the phi acceptance correction histogram");

  THnF* ph = 0;
  TFile* pf = 0;
  int runNumber = 0;

  ~jflucWeightsLoader()
  {
    if (ph)
      delete ph;
    if (pf) {
      pf->Close();
      delete pf;
    }
  }

  void init(InitContext const&)
  {
    if (!doprocessLoadWeights && !doprocessLoadWeightsCF)
      return;
    if (doprocessLoadWeights && doprocessLoadWeightsCF)
      LOGF(fatal, "Only one of JTracks or CFTracks processing can be enabled at a time.");
    if (pathPhiWeights.value.substr(0, 8) == "local://") {
      pf = new TFile(pathPhiWeights.value.substr(8).c_str(), "read");
      if (!pf->IsOpen()) {
        delete pf;
        pf = 0;
        LOGF(fatal, "NUA correction weights file not found: %s", pathPhiWeights.value.substr(8).c_str());
      }
    }
  }

  template <class T>
  using hasDecay = decltype(std::declval<T&>().decay());

  template <class ProducesT, class CollisionT, class TrackT>
  void loadWeights(Produces<ProducesT>& outputT, CollisionT const& collision, TrackT const& tracks)
  {
    if (!pf)
      LOGF(fatal, "NUA correction weights file has not been opened.");
    if (collision.runNumber() != runNumber) {
      if (ph)
        delete ph;
      if (!(ph = static_cast<THnF*>(pf->Get(Form("NUAWeights_%d", collision.runNumber())))))
        LOGF(warning, "NUA correction histogram not found for run %d.", collision.runNumber());
      else
        LOGF(info, "Loaded NUA correction histogram for run %d.", collision.runNumber());
      runNumber = collision.runNumber();
    }
    for (auto& track : tracks) {
      float phiWeight, effWeight;
      if (ph) {
        UInt_t partType = 0; // partType 0 = all charged hadrons
        // TODO: code below to be enabled
        /*if constexpr (std::experimental::is_detected<hasDecay, typename TrackT::iterator>::value) {
          switch (track.decay()) {
            case aod::cf2prongtrack::D0ToPiK:
            case aod::cf2prongtrack::D0barToKPi:
              partType = 1;
              break;
            default:
              break;
          }
        }*/
        const Double_t coords[] = {collision.multiplicity(), static_cast<Double_t>(partType), track.phi(), track.eta(), collision.posZ()};
        phiWeight = ph->GetBinContent(ph->GetBin(coords));
      } else {
        phiWeight = 1.0f;
      }

      effWeight = 1.0f; //<--- todo

      outputT(phiWeight, effWeight);
    }
  }

  Produces<aod::JWeights> output;
  void processLoadWeights(aod::JCollision const& collision, aod::JTracks const& tracks)
  {
    loadWeights(output, collision, tracks);
  }
  PROCESS_SWITCH(jflucWeightsLoader, processLoadWeights, "Load weights histograms for derived data table", false);

  void processLoadWeightsCF(aod::CFCollision const& collision, aod::CFTracks const& tracks)
  {
    loadWeights(output, collision, tracks);
  }
  PROCESS_SWITCH(jflucWeightsLoader, processLoadWeightsCF, "Load weights histograms for CF derived data table", true);

  Produces<aod::J2ProngWeights> output2p;
  void processLoadWeightsCF2Prong(aod::CFCollision const& collision, aod::CF2ProngTracks const& tracks2p)
  {
    loadWeights(output2p, collision, tracks2p);
  }
  PROCESS_SWITCH(jflucWeightsLoader, processLoadWeightsCF2Prong, "Load weights histograms for CF derived 2-prong tracks data table", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<jflucWeightsLoader>(cfgc)};
}
