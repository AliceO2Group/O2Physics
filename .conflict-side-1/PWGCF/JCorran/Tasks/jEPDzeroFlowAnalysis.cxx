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
/// \author junlee.kim@cern.ch
/// \since Jul 2024

#include <experimental/type_traits>
#include <cmath>
#include <array>
#include <cstdlib>
#include <chrono>
#include <string>
#include <vector>

#include "TLorentzVector.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TVector2.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "Math/GenVector/Boost.h"
#include <TMath.h>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/StepTHn.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/StaticFor.h"

#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Qvectors.h"

#include "Common/Core/trackUtilities.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/EventPlaneHelper.h"

#include "CommonConstants/PhysicsConstants.h"

#include "ReconstructionDataFormats/Track.h"

#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"

#include "CCDB/CcdbApi.h"
#include "CCDB/BasicCCDBManager.h"

#include "PWGCF/DataModel/CorrelationsDerived.h"

using namespace std;
using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::constants::physics;

struct jEPDzeroFlowAnalysis {
  enum {
    kFT0C = 0,
    kFT0A = 1,
    kFT0M,
    kFV0A,
    kTPCpos,
    kTPCneg,
    kTPCall
  };

  using MyCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::CentFT0Cs, aod::Qvectors>;
  HistogramRegistry histos{
    "histos",
    {},
    OutputObjHandlingPolicy::AnalysisObject};

  Configurable<float> cfgCentSel{"cfgCentSel", 80., "Centrality selection"};
  Configurable<std::string> cfgCentEst{"cfgCentEst", "FT0C", "Centrality estimator; FT0M or FT0C available"};

  Configurable<bool> cfgPVSel{"cfgPVSel", false, "Additional PV selection flag for syst"};
  Configurable<float> cfgPV{"cfgPV", 8.0, "Additional PV selection range for syst"};
  Configurable<bool> cfgAddEvtSelPileup{"cfgAddEvtSelPileup", false, "flag for additional pileup selection"};
  Configurable<int> cfgMaxOccupancy{"cfgMaxOccupancy", 999999, "maximum occupancy of tracks in neighbouring collisions in a given time range"};
  Configurable<int> cfgMinOccupancy{"cfgMinOccupancy", 0, "maximum occupancy of tracks in neighbouring collisions in a given time range"};

  Configurable<int> cfgnMods{"cfgnMods", 1, "The number of modulations of interest starting from 2"};
  Configurable<int> cfgNQvec{"cfgNQvec", 7, "The number of total Qvectors for looping over the task"};

  Configurable<float> cfgEtaMax{"cfgEtaMax", 0.8, "eta selection"};
  Configurable<float> cfgPtMin{"cfgPtMin", 0.0, "pt selection"};

  Configurable<std::string> cfgDetName{"cfgDetName", "FT0C", "The name of detector to be analyzed"};
  Configurable<std::string> cfgRefAName{"cfgRefAName", "TPCPos", "The name of detector for reference A"};
  Configurable<std::string> cfgRefBName{"cfgRefBName", "TPCNeg", "The name of detector for reference B"};

  ConfigurableAxis massAxis{"massAxis", {175, 1.7, 2.05}, "Invariant mass axis"};
  ConfigurableAxis ptAxis{"ptAxis", {VARIABLE_WIDTH, 0.2, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.5, 8.0, 10.0, 100.0}, "Transverse momentum bins"};
  ConfigurableAxis centAxis{"centAxis", {VARIABLE_WIDTH, 0, 10, 20, 30, 40, 50, 60, 70, 80, 100}, "Centrality interval"};
  ConfigurableAxis cosAxis{"cosAxis", {110, -1.05, 1.05}, "Cosine axis"};

  //  Filter track2pFilter = (nabs(aod::cf2prongtrack::eta) < cfgEtaMax) && (aod::cf2prongtrack::pt > cfgPtMin);

  EventPlaneHelper helperEP;

  int DetId;
  int RefAId;
  int RefBId;

  float centrality;

  template <class T>
  using hasInvMass = decltype(std::declval<T&>().invMass());

  template <typename T>
  int GetDetId(const T& name)
  {
    if (name.value == "FT0C") {
      return kFT0C;
    } else if (name.value == "FT0A") {
      return kFT0A;
    } else if (name.value == "FT0M") {
      return kFT0M;
    } else if (name.value == "FV0A") {
      return kFV0A;
    } else if (name.value == "TPCpos") {
      return kTPCpos;
    } else if (name.value == "TPCneg") {
      return kTPCneg;
    } else if (name.value == "TPCall") {
      return kTPCall;
    } else {
      return 0;
    }
  }

  template <typename TCollision>
  bool eventSelected(TCollision collision)
  {
    if (!collision.sel8()) {
      return false;
    }
    if (cfgCentSel < centrality) {
      return false;
    }
    if (!collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return false;
    }
    if (!collision.selection_bit(aod::evsel::kNoSameBunchPileup)) {
      return false;
    }
    if (cfgPVSel && std::abs(collision.posZ()) > cfgPV) {
      return false;
    }
    if (cfgAddEvtSelPileup && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
      return false;
    }
    if (collision.trackOccupancyInTimeRange() > cfgMaxOccupancy || collision.trackOccupancyInTimeRange() < cfgMinOccupancy) {
      return false;
    }
    return true;
  } // event selection

  template <typename CollType, typename TrackType>
  void fillHistosFlow(const CollType& coll, TrackType& trks)
  {
    if (coll.qvecAmp()[DetId] < 1e-4 || coll.qvecAmp()[RefAId] < 1e-4 || coll.qvecAmp()[RefBId] < 1e-4) {
      return;
    }
    int DetInd = DetId * 4 + cfgNQvec * 4;
    //    int RefAInd = RefAId * 4 + cfgNQvec * 4;
    //    int RefBInd = RefBId * 4 + cfgNQvec * 4;
    for (auto& trk : trks) {
      histos.fill(HIST("hist_EP_cos_Det_v2"), trk.invMass(), trk.pt(), std::cos(2.0 * (trk.phi() - helperEP.GetEventPlane(coll.qvecRe()[DetInd + 3], coll.qvecIm()[DetInd + 3], 2))), centrality);
      histos.fill(HIST("hist_EP_sin_Det_v2"), trk.invMass(), trk.pt(), std::sin(2.0 * (trk.phi() - helperEP.GetEventPlane(coll.qvecRe()[DetInd + 3], coll.qvecIm()[DetInd + 3], 2))), centrality);
    }
  }

  void init(InitContext const&)
  {
    DetId = GetDetId(cfgDetName);
    RefAId = GetDetId(cfgRefAName);
    RefBId = GetDetId(cfgRefBName);

    if (DetId == RefAId || DetId == RefBId || RefAId == RefBId) {
      LOGF(fatal, "Wrong detector configuration \n set the systems correctly");
      DetId = 0;
      RefAId = 4;
      RefBId = 5;
    }

    histos.add(Form("hist_EP_cos_Det_v2"), "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});
    histos.add(Form("hist_EP_sin_Det_v2"), "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});
  }

  void processData(MyCollisions::iterator const& collision, aod::CF2ProngTracks const& p2tracks)
  {
    if (cfgCentEst.value == "FT0C") {
      centrality = collision.centFT0C();
    } else if (cfgCentEst.value == "FT0M") {
      centrality = collision.centFT0M();
    }
    if (!eventSelected(collision)) {
      return;
    }
    fillHistosFlow(collision, p2tracks);
  }
  PROCESS_SWITCH(jEPDzeroFlowAnalysis, processData, "Process Event for data", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<jEPDzeroFlowAnalysis>(cfgc)};
}
