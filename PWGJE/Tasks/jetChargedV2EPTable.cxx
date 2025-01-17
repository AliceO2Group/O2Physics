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

// jet v2 task
/// \author Yubiao Wang <yubiao.wang@cern.ch>
// C++/ROOT includes.

#include <chrono>
#include <string>
#include <vector>
#include <TComplex.h>
#include <TH3F.h>
#include <cmath>
#include <TRandom3.h>
#include <TF1.h>
#include <algorithm>
#include <unordered_map>
// o2Physics includes.

#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPMagField.h"

#include "Framework/runDataProcessing.h"

#include "Common/DataModel/FT0Corrected.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"
#include "Common/CCDB/ctpRateFetcher.h"

//< evt pln .h >//
#include "Framework/ASoAHelpers.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/StaticFor.h"

#include "Common/DataModel/Qvectors.h"
#include "Common/Core/EventPlaneHelper.h"
//< evt pln .h | end >//

// o2 includes.
#include "DetectorsCommonDataFormats/AlignParam.h"

#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/HistogramRegistry.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "PWGJE/Core/FastJetUtilities.h"
#include "PWGJE/Core/JetFinder.h"
#include "PWGJE/Core/JetFindingUtilities.h"
#include "PWGJE/DataModel/Jet.h"

#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "EventFiltering/filterTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace o2::aod
{
DECLARE_SOA_TABLE(MyCollisions, "AOD", "MYCOLLISION", //! vertex information of collision
                  o2::soa::Index<>, collision::PosZ);
using MyCollision = MyCollisions::iterator;

  namespace myTable 
  {
    DECLARE_SOA_INDEX_COLUMN(MyCollision, mycollision);

    DECLARE_SOA_COLUMN(Psi2, psi2, Float_t);
    DECLARE_SOA_COLUMN(Psi3, psi3, Float_t);
    DECLARE_SOA_COLUMN(EvtPlRes, evtplres, Float_t);
  }
  DECLARE_SOA_TABLE(MyTable, "AOD", "MYTABLE", o2::soa::Index<>,
                    myTable::Psi2,
                    myTable::Psi3,
                    myTable::EvtPlRes);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
struct jetChargedV2EPTable {
  Produces<o2::aod::MyTable> myTable;
  Produces<o2::aod::MyCollision> outputCollisions;

  HistogramRegistry registry;


  Configurable<std::string> eventSelections{"eventSelections", "sel8", "choose event selection"};
  Configurable<std::string> trackSelections{"trackSelections", "globalTracks", "set track selections"};

  Configurable<float> vertexZCut{"vertexZCut", 10.0f, "Accepted z-vertex range"};
  Configurable<float> centralityMin{"centralityMin", -999.0, "minimum centrality"};
  Configurable<float> centralityMax{"centralityMax", 999.0, "maximum centrality"};
  Configurable<float> trackPtMin{"trackPtMin", 0.15, "minimum pT acceptance for tracks"};
  Configurable<float> trackPtMax{"trackPtMax", 1000., "maximum pT acceptance for tracks"};
  Configurable<float> trackEtaMin{"trackEtaMin", -0.9, "minimum eta acceptance for tracks"};
  Configurable<float> trackEtaMax{"trackEtaMax", 0.9, "maximum eta acceptance for tracks"};

  Configurable<float> jetAreaFractionMin{"jetAreaFractionMin", -99.0, "used to make a cut on the jet areas"};
  Configurable<float> leadingConstituentPtMin{"leadingConstituentPtMin", -99.0, "minimum pT selection on jet constituent"};
  Configurable<float> jetPtMin{"jetPtMin", 0.15, "minimum pT acceptance for jets"};
  Configurable<float> jetPtMax{"jetPtMax", 200.0, "maximum pT acceptance for jets"};

  Configurable<int> trackOccupancyInTimeRangeMax{"trackOccupancyInTimeRangeMax", 999999, "maximum occupancy of tracks in neighbouring collisions in a given time range; only applied to reconstructed collisions (data and mcd jets), not mc collisions (mcp jets)"};
  Configurable<int> trackOccupancyInTimeRangeMin{"trackOccupancyInTimeRangeMin", -999999, "minimum occupancy of tracks in neighbouring collisions in a given time range; only applied to reconstructed collisions (data and mcd jets), not mc collisions (mcp jets)"};

  //=====================< evt pln >=====================//
  Configurable<bool> cfgAddEvtSel{"cfgAddEvtSel", true, "event selection"};
  Configurable<std::vector<int>> cfgnMods{"cfgnMods", {2}, "Modulation of interest"};
  Configurable<int> cfgnTotalSystem{"cfgnTotalSystem", 7, "total qvect 56  or number"};
  Configurable<std::string> cfgDetName{"cfgDetName", "FT0M", "The name of detector to be analyzed"};
  Configurable<std::string> cfgRefAName{"cfgRefAName", "TPCpos", "The name of detector for reference A"};
  Configurable<std::string> cfgRefBName{"cfgRefBName", "TPCneg", "The name of detector for reference B"};

  EventPlaneHelper helperEP;
  int DetId;
  int RefAId;
  int RefBId;

  template <typename T>
  int GetDetId(const T& name)
  {
    if (name.value == "BPos" || name.value == "BNeg" || name.value == "BTot") {
      LOGF(warning, "Using deprecated label: %s. Please use TPCpos, TPCneg, TPCall instead.", name.value);
    }
    if (name.value == "FT0C") {
      return 0;
    } else if (name.value == "FT0A") {
      return 1;
    } else if (name.value == "FT0M") {
      return 2;
    } else if (name.value == "FV0A") {
      return 3;
    } else if (name.value == "TPCpos" || name.value == "BPos") {
      return 4;
    } else if (name.value == "TPCneg" || name.value == "BNeg") {
      return 5;
    } else if (name.value == "TPCall" || name.value == "BTot") {
      return 6;
    } else {
      return 0;
    }
  }
  //=====================< evt pln | end >=====================//

  Configurable<float> selectedJetsRadius{"selectedJetsRadius", 0.4, "resolution parameter for histograms without radius"};

  std::vector<double> jetPtBins;
  std::vector<double> jetPtBinsRhoAreaSub;

  int eventSelectiont = -1;
  int trackSelectiont = -1;

  void init(o2::framework::InitContext&)
  {
    DetId = GetDetId(cfgDetName);
    RefAId = GetDetId(cfgRefAName);
    RefBId = GetDetId(cfgRefBName);
    if (DetId == RefAId || DetId == RefBId || RefAId == RefBId) {
      LOGF(info, "Wrong detector configuration \n The FT0C will be used to get Q-Vector \n The TPCpos and TPCneg will be used as reference systems");
      DetId = 0;
      RefAId = 4;
      RefBId = 5;
    }
    eventSelectiont = jetderiveddatautilities::initialiseEventSelection(static_cast<std::string>(eventSelections));
    trackSelectiont = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));
  }

  Filter trackCuts = (aod::jtrack::pt >= trackPtMin && aod::jtrack::pt < trackPtMax && aod::jtrack::eta > trackEtaMin && aod::jtrack::eta < trackEtaMax);
  Filter eventCuts = (nabs(aod::jcollision::posZ) < vertexZCut && aod::jcollision::centrality >= centralityMin && aod::jcollision::centrality < centralityMax);

  void processJetV2Table(soa::Filtered<soa::Join<aod::JetCollisions, aod::BkgChargedRhos, aod::Qvectors>>::iterator const& collision,
                         soa::Join<aod::ChargedJets, aod::ChargedJetConstituents> const& jets,
                         aod::JetTracks const& tracks)
  {
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectiont)) {
      return;
    }
    outputCollisions(collision.posZ());
    //=====================< evt pln [n=2->\Psi_2, n=3->\Psi_3] >=====================//
    double ep2 = 0;
    double ep3 = 0;
    double res2 = 0;

    for (uint i = 0; i < cfgnMods->size(); i++) {
      int nmode = cfgnMods->at(i);
      int DetInd = DetId * 4 + cfgnTotalSystem * 4 * (nmode - 2);
      if (nmode == 2) {
        if (collision.qvecAmp()[DetId] > 1e-8) {
          ep2 = helperEP.GetEventPlane(collision.qvecRe()[DetInd + 3], collision.qvecIm()[DetInd + 3], nmode);
          res2 = helperEP.GetResolution(helperEP.GetEventPlane(collision.qvecRe()[DetInd + 3], collision.qvecIm()[DetInd + 3], nmode), helperEP.GetEventPlane(collision.qvecRe()[DetInd + 3], collision.qvecIm()[DetInd + 3], nmode), nmode);
        }
      } else if (nmode == 3) {
        if (collision.qvecAmp()[DetId] > 1e-8) {
          ep3 = helperEP.GetEventPlane(collision.qvecRe()[DetInd + 3], collision.qvecIm()[DetInd + 3], nmode);
        }
      }
    }
    myTable(ep2, ep3, res2);
    LOGF(info, "table producer Psi2 as %f, Psi3 as %f, Res2 as %f", ep2, ep3, res2);
  }
  PROCESS_SWITCH(jetChargedV2EPTable, processJetV2Table, "Jet V2 table for fit", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<jetChargedV2EPTable>(cfgc, TaskName{"jet-charged-v2ept"})};
}
