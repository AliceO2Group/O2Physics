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
/// \file jflucWeightsLoader.cxx
/// \brief Task to load the NUA and NUE weights from local files or CCDB.
/// \author Jasper Parkkila (jparkkil@cern.ch), Maxim Virta (maxim.virta@cern.ch)
/// \since May 2024
/// The weights are loaded from the local files or CCDB and stored in the JWeights table.

#include "PWGCF/DataModel/CorrelationsDerived.h"
#include "PWGCF/JCorran/DataModel/JCatalyst.h"

#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/V0.h"

#include <TFile.h>
#include <THn.h>

#include <experimental/type_traits>
#include <string>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};

// The standalone jfluc code expects the entire list of tracks for an event. At the same time, it expects weights together with other track attributes.
// This workflow creates a table of weights that can be joined with track tables.
struct JflucWeightsLoader {
  O2_DEFINE_CONFIGURABLE(cfgPathPhiWeights, std::string, "Users/m/mavirta/corrections/NUA/LHC23zzh", "Local (local://) or CCDB path for the phi acceptance correction histogram");
  O2_DEFINE_CONFIGURABLE(cfgPathEffWeights, std::string, "Users/m/mavirta/corrections/NUE/LHC23zzh", "Local (local://) or CCDB path for the efficiency correction histogram");
  O2_DEFINE_CONFIGURABLE(cfgForRunNumber, bool, false, "Get CCDB object by run");

  THnF* ph = 0;
  TFile* pf = 0;
  THnF* pheff = 0;
  TFile* pfeff = 0;
  int runNumber = 0;
  int timestamp = 0;
  bool useNUAFromCCDB = false;
  bool useEffFromCCDB = false;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  std::string ccdbURL = "http://alice-ccdb.cern.ch";
  enum { kNUA,
         kEFF };

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

  void initCCDB(int runNum, int ts, int NUAorEFF = kNUA)
  {
    if (cfgForRunNumber) {
      if (NUAorEFF == kNUA) {
        ph = ccdb->getForRun<THnF>(cfgPathPhiWeights, runNum);
      } else {
        pheff = ccdb->getForRun<THnF>(cfgPathEffWeights, runNum);
      }
    } else {
      if (NUAorEFF == kNUA) {
        ph = ccdb->getForTimeStamp<THnF>(cfgPathPhiWeights, ts);
      } else {
        pheff = ccdb->getForTimeStamp<THnF>(cfgPathEffWeights, ts);
      }
    }
  }

  void init(InitContext const&)
  {
    if (!doprocessLoadWeights && !doprocessLoadWeightsCF) {
      return;
    }

    if (doprocessLoadWeights && doprocessLoadWeightsCF)
      LOGF(fatal, "Only one of JTracks or CFTracks processing can be enabled at a time.");

    // NUA corrections from local file or CCDB
    if (cfgPathPhiWeights.value.substr(0, 8) == "local://") {
      LOGF(info, "Using NUA corrections locally from: %s", cfgPathPhiWeights.value.substr(8).c_str());
      pf = new TFile(cfgPathPhiWeights.value.substr(8).c_str(), "read");
      if (!pf->IsOpen()) {
        delete pf;
        pf = 0;
        LOGF(fatal, "NUA correction weights file not found: %s", cfgPathPhiWeights.value.substr(8).c_str());
      }
      useNUAFromCCDB = false;
    } else if (cfgPathPhiWeights.value == "") {
      LOGF(info, "No NUA corrections provided.");
      useNUAFromCCDB = false;
    } else {
      LOGF(info, "Assuming NUA corrections from CCDB.");
      useNUAFromCCDB = true;
      ccdb->setURL(ccdbURL.data()); // default CCDB URL
      ccdb->setCaching(true);
      ccdb->setLocalObjectValidityChecking();
      ccdb->setFatalWhenNull(false);
    }

    // Efficiency corrections from local file or CCDB
    if (cfgPathEffWeights.value.substr(0, 8) == "local://") {
      LOGF(info, "Using efficiency corrections locally from: %s", cfgPathEffWeights.value.substr(8).c_str());
      pfeff = new TFile(cfgPathEffWeights.value.substr(8).c_str(), "read");
      if (!pfeff->IsOpen()) {
        delete pfeff;
        pfeff = 0;
        LOGF(fatal, "Efficiency correction weights file not found: %s", cfgPathEffWeights.value.substr(8).c_str());
      } else {
        LOGF(info, "Loaded efficiency correction histogram locally.");
      }
      useEffFromCCDB = false;
    } else if (cfgPathEffWeights.value == "") {
      LOGF(info, "No efficiency corrections provided.");
      useEffFromCCDB = false;
    } else {
      LOGF(info, "Assuming efficiency corrections from CCDB.");
      useEffFromCCDB = true;
      // If NUA corrections are from CCDB, use the same CCDB URL for efficiency corrections
      if (!useNUAFromCCDB) {
        ccdb->setURL(ccdbURL.data()); // default CCDB URL
        ccdb->setCaching(true);
        ccdb->setLocalObjectValidityChecking();
        ccdb->setFatalWhenNull(false);
      }
    }
  }

  template <class T>
  using HasDecay = decltype(std::declval<T&>().decay());

  template <class ProducesT, class CollisionT, class TrackT>
  void loadWeights(Produces<ProducesT>& outputT, CollisionT const& collision, TrackT const& tracks)
  {
    if (pf || useNUAFromCCDB) {
      if (collision.runNumber() != runNumber) {
        if (ph)
          delete ph;
        if (!useNUAFromCCDB) {
          // Check if NUA correction can be found from a local file and load it
          if (!(ph = pf->Get<THnF>(Form("NUAWeights_%d", collision.runNumber()))))
            LOGF(warning, "NUA correction histogram not found for run %d.", collision.runNumber());
          else
            LOGF(info, "Loaded NUA correction histogram locally for run %d.", collision.runNumber());
        } else {
          initCCDB(collision.runNumber(), timestamp, kNUA);
          LOGF(info, "Loaded NUA correction histogram from CCDB for run %d.", collision.runNumber());
        }
      }
    }
    if (pfeff || useEffFromCCDB) {
      if (collision.runNumber() != runNumber) {
        if (pheff)
          delete pheff;
        if (!useEffFromCCDB) {
          if (!(pheff = pfeff->Get<THnF>("ccdb_object"))) {
            LOGF(warning, "Efficiency correction histogram not found.");
          } else {
            LOGF(info, "Loaded NUE correction histogram locally for run %d.", collision.runNumber());
          }
        } else {
          initCCDB(collision.runNumber(), timestamp, kEFF);
          LOGF(info, "Loaded efficiency correction histogram from CCDB for run %d.", collision.runNumber());
        }
      }
    }

    // Set run number after reading corrections
    runNumber = collision.runNumber();

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
        // NUA corrections are a function of multiplicity, partType, phi, eta, and z-vertex
        const double nuaCoords[] = {collision.multiplicity(), static_cast<double>(partType), track.phi(), track.eta(), collision.posZ()};
        phiWeight = ph->GetBinContent(ph->GetBin(nuaCoords));
      } else {
        phiWeight = 1.0f;
      }

      if (pheff) {
        // Efficiency corrections are a function of eta, pT, multiplicity, and z-vertex
        const double nueCoords[] = {track.eta(), track.pt(), collision.multiplicity(), collision.posZ()};

        effWeight = pheff->GetBinContent(pheff->GetBin(nueCoords));
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
