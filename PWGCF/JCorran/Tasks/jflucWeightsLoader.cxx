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
// o2-linter: disable='doc/file'

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

#include "CCDB/BasicCCDBManager.h"

#include "PWGCF/JCorran/DataModel/JCatalyst.h"
#include "PWGCF/DataModel/CorrelationsDerived.h"
#include "Framework/runDataProcessing.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};

// The standalone jfluc code expects the entire list of tracks for an event. At the same time, it expects weights together with other track attributes.
// This workflow creates a table of weights that can be joined with track tables.
struct JflucWeightsLoader {
  O2_DEFINE_CONFIGURABLE(cfgPathPhiWeights, std::string, "http://alice-ccdb.cern.ch", "Local (local://) or CCDB path for the phi acceptance correction histogram");
  O2_DEFINE_CONFIGURABLE(cfgPathEffWeights, std::string, "", "Local (local://) or CCDB path for the efficiency correction histogram");
  O2_DEFINE_CONFIGURABLE(cfgForRunNumber, bool, false, "Get CCDB object by run");
  O2_DEFINE_CONFIGURABLE(cfgCCDBPath, std::string, "Users/m/mavirta/corrections/NUA/LHC23zzh", "Internal path in CCDB");

  THnF* ph = 0;
  TFile* pf = 0;
  THnD* pheff = 0;
  TFile* pfeff = 0;
  int runNumber = 0;
  int timestamp = 0;
  bool useCCDB = false;
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  ~JflucWeightsLoader()
  {
    if (ph)
      delete ph;
    if (pf) {
      pf->Close();
      delete pf;
    }
    if (pheff)
      delete pheff;
    if (pfeff) {
      pfeff->Close();
      delete pfeff;
    }
  }

  void initCCDB(int runNum, int ts)
  {
    if (cfgForRunNumber) {
      ph = ccdb->getForRun<THnF>(cfgCCDBPath, runNum);
    } else {
      ph = ccdb->getForTimeStamp<THnF>(cfgCCDBPath, ts);
    }
  }

  void init(InitContext const&)
  {
    if (!doprocessLoadWeights && !doprocessLoadWeightsCF) {
      return;
    }

    if (doprocessLoadWeights && doprocessLoadWeightsCF)
      LOGF(fatal, "Only one of JTracks or CFTracks processing can be enabled at a time.");
    if (cfgPathPhiWeights.value.find("ccdb") != std::string::npos) {
      LOGF(info, "Using corrections from: ccdb");
      useCCDB = true;
      ccdb->setURL(cfgPathPhiWeights);
      ccdb->setCaching(true);
      ccdb->setLocalObjectValidityChecking();
      ccdb->setFatalWhenNull(false);
    } else if (cfgPathPhiWeights.value.substr(0, 8) == "local://") {
      LOGF(info, "Using non-uniform acceptance corrections from: %s", cfgPathPhiWeights.value.substr(8).c_str());
      pf = new TFile(cfgPathPhiWeights.value.substr(8).c_str(), "read");
      if (!pf->IsOpen()) {
        delete pf;
        pf = 0;
        LOGF(fatal, "NUA correction weights file not found: %s", cfgPathPhiWeights.value.substr(8).c_str());
      }
      useCCDB = false;
    } else {
      LOGF(info, "Didn't find \"local://\" or \"ccdb\" for non-uniform acceptance corrections.");
      return;
    }

    if (cfgPathEffWeights.value.substr(0, 8) == "local://") {
      LOGF(info, "Using efficiency corrections from: %s", cfgPathEffWeights.value.substr(8).c_str());
      pfeff = new TFile(cfgPathEffWeights.value.substr(8).c_str(), "read");
      if (!pfeff->IsOpen()) {
        delete pfeff;
        pfeff = 0;
        LOGF(fatal, "Efficiency correction weights file not found: %s", cfgPathEffWeights.value.substr(8).c_str());
      }
      //
      if (!(pheff = pfeff->Get<THnD>("ccdb_object"))) {
        LOGF(warning, "Efficiency correction histogram not found.");
      } else {
        LOGF(info, "Loaded efficiency correction histogram locally.");
      }
    } else {
      LOGF(info, "Didn't find \"local://\" or \"ccdb\" for efficiency corrections.");
      return;
    }
  }

  template <class T>
  using HasDecay = decltype(std::declval<T&>().decay());

  template <class ProducesT, class CollisionT, class TrackT>
  void loadWeights(Produces<ProducesT>& outputT, CollisionT const& collision, TrackT const& tracks)
  {
    if (pf || useCCDB) {
      if (collision.runNumber() != runNumber) {
        if (ph)
          delete ph;
        if (!useCCDB) {
          // Check if NUA correction can be found from a local file and load it
          if (!(ph = pf->Get<THnF>(Form("NUAWeights_%d", collision.runNumber()))))
            LOGF(warning, "NUA correction histogram not found for run %d.", collision.runNumber());
          else
            LOGF(info, "Loaded NUA correction histogram locally for run %d.", collision.runNumber());
        } else {
          initCCDB(collision.runNumber(), timestamp);
          LOGF(info, "Loaded NUA correction histogram from CCDB for run %d.", collision.runNumber());
        }
        runNumber = collision.runNumber();
      }
    }
    for (const auto& track : tracks) {
      float phiWeight, effWeight;
      if (ph) {
        uint partType = 0; // partType 0 = all charged hadrons
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
        const double coords[] = {collision.multiplicity(), static_cast<double>(partType), track.phi(), track.eta(), collision.posZ()};
        phiWeight = ph->GetBinContent(ph->GetBin(coords));
      } else {
        phiWeight = 1.0f;
      }

      if (pheff) {
        const int effVars[] = {
          pheff->GetAxis(0)->FindBin(track.eta()),
          pheff->GetAxis(1)->FindBin(track.pt()),
          pheff->GetAxis(2)->FindBin(collision.multiplicity()),
          pheff->GetAxis(3)->FindBin(collision.posZ())};
        effWeight = pheff->GetBinContent(effVars);
      } else {
        effWeight = 1.0f;
      }

      outputT(phiWeight, effWeight);
    }
  }

  Produces<aod::JWeights> output;
  void processLoadWeights(aod::JCollision const& collision, aod::JTracks const& tracks)
  {
    loadWeights(output, collision, tracks);
  }
  PROCESS_SWITCH(JflucWeightsLoader, processLoadWeights, "Load weights histograms for derived data table", false);

  void processLoadWeightsCF(aod::CFCollision const& collision, aod::CFTracks const& tracks)
  {
    loadWeights(output, collision, tracks);
  }
  PROCESS_SWITCH(JflucWeightsLoader, processLoadWeightsCF, "Load weights histograms for CF derived data table", true);

  Produces<aod::J2ProngWeights> output2p;
  void processLoadWeightsCF2Prong(aod::CFCollision const& collision, aod::CF2ProngTracks const& tracks2p)
  {
    loadWeights(output2p, collision, tracks2p);
  }
  PROCESS_SWITCH(JflucWeightsLoader, processLoadWeightsCF2Prong, "Load weights histograms for CF derived 2-prong tracks data table", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<JflucWeightsLoader>(cfgc)};
}
