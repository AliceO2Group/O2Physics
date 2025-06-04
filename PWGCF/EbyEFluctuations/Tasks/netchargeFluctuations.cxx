//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
//
// author Nida Malik (nida.malik@cern.ch)
// Department of Physics, Aligarh Muslim University, India
// to study the net charge fluctuations by observable, #nu_dyn

#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/Core/trackUtilities.h"
#include "Common/CCDB/EventSelectionParams.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/Centrality.h"
#include "CommonConstants/MathConstants.h"
#include "Common/DataModel/FT0Corrected.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/RunningWorkflowInfo.h"
#include "PWGCF/Core/CorrelationContainer.h"
#include "PWGCF/Core/PairCuts.h"
#include "Common/CCDB/TriggerAliases.h"
#include "Framework/HistogramRegistry.h"
#include "TProfile.h"
#include "TProfile2D.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace std;

namespace o2::aod
{
namespace netCharge
{
DECLARE_SOA_COLUMN(PosCharge, pos_charge, float);
DECLARE_SOA_COLUMN(NegCharge, neg_charge, float);
DECLARE_SOA_COLUMN(PosSqCharge, possq_charge, float);
DECLARE_SOA_COLUMN(NegSqCharge, negsq_charge, float);
DECLARE_SOA_COLUMN(TermpCharge, termp_charge, float);
DECLARE_SOA_COLUMN(TermnCharge, termn_charge, float);
DECLARE_SOA_COLUMN(PosNegCharge, posneg_charge, float);
DECLARE_SOA_COLUMN(Centrality, centrality, float);
} // namespace netCharge
DECLARE_SOA_TABLE(NetCharge, "AOD", "NETChargefluct",
                  netCharge::PosCharge,
                  netCharge::NegCharge,
                  netCharge::PosSqCharge,
                  netCharge::NegSqCharge,
                  netCharge::TermpCharge,
                  netCharge::TermnCharge,
                  netCharge::PosNegCharge,
                  netCharge::Centrality);
//------------------------------------------------------

using MyCollisionsRun2 = soa::Join<aod::Collisions, aod::EvSels, aod::CentRun2V0Ms, aod::Mults>;
using MyCollisionRun2 = MyCollisionsRun2::iterator;
using MyCollisionsRun3 = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::CentFT0Cs, aod::Mults>;
using MyCollisionRun3 = MyCollisionsRun3::iterator;
using MyTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::StoredTracks, aod::TrackSelection>;
using MyTrack = MyTracks::iterator;

} // namespace o2::aod

enum RunType {
  kRun3 = 0,
  kRun2
};

struct NetchargeFluctuations {
  Produces<aod::NetCharge> net_charge;

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // Configurables
  Configurable<float> cVtxZcut{"cVtx", 10.f, "Vertex Z"};
  Configurable<float> cEtacut{"cEta", 0.8, "Eta cut"};
  Configurable<float> cPtmincut{"cPtmincut", 0.2, "Pt min cut"};
  Configurable<float> cPtmaxcut{"cPtmaxcut", 5.0, "Pt max cut"};
  Configurable<float> cDcaXYcut{"cDcaXYcut", 0.12, "DCA XY cut"};
  Configurable<float> cDcaZcut{"cDcaZcut", 0.3, "DCA Z cut"};
  Configurable<float> cCentmincut{"cCentmincut", 0.0, "Min cent cut"};
  Configurable<float> cCentmaxcut{"cCentmaxcut", 90.0, "Max cent cut"};
  Configurable<int> cTPCcrosscut{"cTPCcrosscut", 70, "TPC crossrows cut"};
  Configurable<int> cITSchicut{"cITSchi2clustercut", 70, "ITS chi2 cluster cut"};
  Configurable<int> cTPCchicut{"cTPCchi2clustercut", 70, "TPC chi2 cluster cut"};

  // Event selections
  Configurable<bool> cSel8Trig{"cSel8Trig", true, "Sel8 (T0A + T0C) Selection Run3"};    // sel8
  Configurable<bool> cInt7Trig{"cInt7Trig", false, "kINT7 MB Trigger"};                  // kINT7
  Configurable<bool> cSel7Trig{"cSel7Trig", false, "Sel7 (V0A + V0C) Selection Run2"};   // sel7
  Configurable<bool> cTFBorder{"cTFBorder", false, "Timeframe Border Selection"};        // pileup
  Configurable<bool> cNoItsROBorder{"cNoItsROBorder", false, "No ITSRO Border Cut"};     // pileup
  Configurable<bool> cItsTpcVtx{"cItsTpcVtx", false, "ITS+TPC Vertex Selection"};        // pileup
  Configurable<bool> cPileupReject{"cPileupReject", false, "Pileup rejection"};          // pileup
  Configurable<bool> cZVtxTimeDiff{"cZVtxTimeDiff", false, "z-vtx time diff selection"}; // pileup
  Configurable<bool> cIsGoodITSLayers{"cIsGoodITSLayers", false, "Good ITS Layers All"}; // pileup

  // Initialization
  float cent = 0.;
  float mult = 0.;
  void init(o2::framework::InitContext&)
  {

    const AxisSpec vtxZAxis = {80, -20, 20, "V_{Z} (cm)"};
    const AxisSpec dcaAxis = {250, -0.5, 0.5, "DCA_{xy} (cm)"};
    const AxisSpec dcazAxis = {250, -0.5, 0.5, "DCA_{z} (cm)"};
    const AxisSpec ptAxis = {70, 0.0, 7.0, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec etaAxis = {20, -1., 1., "#eta"};
    const AxisSpec centAxis = {100, 0., 100., "centrality"};
    const AxisSpec multAxis = {200, 0., 10000., "FT0M Amplitude"};
    const AxisSpec TPCChi2Axis = {140, 0., 7., "Chi2"};
    const AxisSpec ITSChi2Axis = {80, 0., 40., "Chi2"};
    const AxisSpec CrossedrowTPCAxis = {160, 0., 160., "TPC Crossed rows"};
    const AxisSpec eventsAxis = {10, 0, 10, ""};

    // Before cuts
    histos.add("hVtxZ_before", "", kTH1F, {vtxZAxis});
    histos.add("hDcaXY_before", "", kTH1F, {dcaAxis});
    histos.add("hDcaZ_before", "", kTH1F, {dcazAxis});
    histos.add("hTPCchi2perCluster_before", "", kTH1D, {TPCChi2Axis});
    histos.add("hITSchi2perCluster_before", "", kTH1D, {ITSChi2Axis});
    histos.add("hTPCCrossedrows_before", "", kTH1D, {CrossedrowTPCAxis});
    histos.add("hPtDcaXY_before", "", kTH2D, {ptAxis, dcaAxis});
    histos.add("hPtDcaZ_before", "", kTH2D, {ptAxis, dcazAxis});

    // After cuts
    histos.add("hVtxZ_after", "", kTH1F, {vtxZAxis});
    histos.add("hDcaXY_after", "", kTH1F, {dcaAxis});
    histos.add("hDcaZ_after", "", kTH1F, {dcazAxis});
    histos.add("hTPCchi2perCluster_after", "", kTH1D, {TPCChi2Axis});
    histos.add("hITSchi2perCluster_after", "", kTH1D, {ITSChi2Axis});
    histos.add("hTPCCrossedrows_after", "", kTH1D, {CrossedrowTPCAxis});
    histos.add("hPtDcaXY_after", "", kTH2D, {ptAxis, dcaAxis});
    histos.add("hPtDcaZ_after", "", kTH2D, {ptAxis, dcazAxis});
    histos.add("hEta", "", kTH1F, {etaAxis});
    histos.add("hPt", "", kTH1F, {ptAxis});
    histos.add("hCentrality", "", kTH1F, {centAxis});
    histos.add("hMultiplicity", "", kTH1F, {multAxis});
  }

  template <RunType run, typename C>
  bool selCollision(C const& coll)
  {
    if (std::abs(coll.posZ()) > cVtxZcut) {
      return false;
    } // Reject the collisions with large vertex-z

    if constexpr (run == kRun3) {
      cent = coll.centFT0M(); // centrality for run3
      mult = coll.multFT0M(); // multiplicity for run3
      if (cSel8Trig && !coll.sel8()) {
        return false;
      } // require min bias trigger
    } else {
      cent = coll.centRun2V0M(); // centrality for run2
      mult = coll.multFV0M();    // multiplicity for run2
      if (cInt7Trig && !coll.alias_bit(kINT7)) {
        return false;
      }
      if (cSel7Trig && !coll.sel7()) {
        return false;
      }
    }

    if (cNoItsROBorder && !coll.selection_bit(aod::evsel::kNoITSROFrameBorder)) {
      return false;
    }

    if (cTFBorder && !coll.selection_bit(aod::evsel::kNoTimeFrameBorder)) {
      return false;
    }

    if (cPileupReject && !coll.selection_bit(aod::evsel::kNoSameBunchPileup)) {
      return false;
    }

    if (cZVtxTimeDiff && !coll.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return false;
    }

    if (cItsTpcVtx && !coll.selection_bit(aod::evsel::kIsVertexITSTPC)) {
      return false;
    }

    return true; // if all checks pass, accept the collision
  }

  template <typename T>
  bool selTrack(T const& track)
  {
    if (!track.isGlobalTrack()) {
      return false;
    } // accept only global tracks

    if (std::fabs(track.dcaXY()) > cDcaXYcut) {
      return false;
    }

    if (std::fabs(track.dcaZ()) > cDcaZcut) {
      return false;
    }

    if (std::fabs(track.eta()) >= cEtacut) {
      return false;
    }

    if (track.pt() <= cPtmincut || track.pt() >= cPtmaxcut) {
      return false;
    }

    if (track.tpcNClsCrossedRows() < cTPCcrosscut) {
      return false;
    }

    if (track.itsChi2NCl() > cITSchicut) {
      return false;
    }

    if (track.tpcChi2NCl() > cTPCchicut) {
      return false;
    }

    return true; // if all checks pass, accept the collision
  }

  template <RunType run, typename C, typename T>
  void calculation(C const& coll, T const& tracks)
  {
    histos.fill(HIST("hVtxZ_before"), coll.posZ());

    if (!selCollision<run>(coll)) {
      return;
    }

    histos.fill(HIST("hVtxZ_after"), coll.posZ());
    histos.fill(HIST("hCentrality"), cent);
    histos.fill(HIST("hMultiplicity"), mult);

    int fpos = 0, fneg = 0, posneg = 0, termn = 0, termp = 0;

    for (auto track : tracks) {
      histos.fill(HIST("hTPCchi2perCluster_before"), track.tpcChi2NCl());
      histos.fill(HIST("hITSchi2perCluster_before"), track.itsChi2NCl());
      histos.fill(HIST("hTPCCrossedrows_before"), track.tpcNClsCrossedRows());
      histos.fill(HIST("hDcaXY_before"), track.dcaXY());
      histos.fill(HIST("hDcaZ_before"), track.dcaZ());
      histos.fill(HIST("hPtDcaXY_before"), track.pt(), track.dcaXY());
      histos.fill(HIST("hPtDcaZ_before"), track.pt(), track.dcaZ());

      if (!selTrack(track)) {
        continue;
      }

      histos.fill(HIST("hDcaXY_after"), track.dcaXY());
      histos.fill(HIST("hDcaZ_after"), track.dcaZ());
      histos.fill(HIST("hPt"), track.pt());
      histos.fill(HIST("hEta"), track.eta());
      histos.fill(HIST("hPtDcaXY_after"), track.pt(), track.dcaXY());
      histos.fill(HIST("hPtDcaZ_after"), track.pt(), track.dcaZ());
      histos.fill(HIST("hTPCCrossedrows_after"), track.tpcNClsCrossedRows());
      histos.fill(HIST("hTPCchi2perCluster_after"), track.tpcChi2NCl());
      histos.fill(HIST("hITSchi2perCluster_after"), track.itsChi2NCl());

      if (track.sign() > 0) {
        fpos += 1;
        termp = fpos * (fpos - 1);
      }
      if (track.sign() < 0) {
        fneg += 1;
        termn = fneg * (fneg - 1);
      }
      posneg = fpos * fneg;
      net_charge(fpos, fneg, fpos * fpos, fneg * fneg, termp, termn, posneg, cent);
    } // tracks

    return;
  }

  void processDataRun3(aod::MyCollisionRun3 const& coll, aod::MyTracks const& tracks)
  {
    calculation<kRun3>(coll, tracks);
  }

  PROCESS_SWITCH(NetchargeFluctuations, processDataRun3, "Process for Run3 DATA", true);

  void processDataRun2(aod::MyCollisionRun2 const& coll, aod::MyTracks const& tracks)
  {
    calculation<kRun2>(coll, tracks);
  }

  PROCESS_SWITCH(NetchargeFluctuations, processDataRun2, "Process for Run2 DATA", false);

}; // struct

struct Netcharge_analysis {
  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};
  std::vector<std::vector<std::shared_ptr<TProfile>>> net;

  void init(o2::framework::InitContext&)
  {
    std::vector<double> centBinning = {0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100};
    AxisSpec centAxis = {centBinning, "centrality"};

    registry.add("pos_vs_cent", "", {HistType::kTProfile, {centAxis}});
    registry.add("neg_vs_cent", "", {HistType::kTProfile, {centAxis}});
    registry.add("termp_vs_cent", "", {HistType::kTProfile, {centAxis}});
    registry.add("termn_vs_cent", "", {HistType::kTProfile, {centAxis}});
    registry.add("pos_sq_vs_cent", "", {HistType::kTProfile, {centAxis}});
    registry.add("neg_sq_vs_cent", "", {HistType::kTProfile, {centAxis}});
    registry.add("posneg_vs_cent", "", {HistType::kTProfile, {centAxis}});
  }

  void process(aod::NetCharge::iterator const& event_netcharge)
  {
    registry.get<TProfile>(HIST("pos_vs_cent"))->Fill(event_netcharge.centrality(), event_netcharge.pos_charge());
    registry.get<TProfile>(HIST("neg_vs_cent"))->Fill(event_netcharge.centrality(), event_netcharge.neg_charge());
    registry.get<TProfile>(HIST("termp_vs_cent"))->Fill(event_netcharge.centrality(), event_netcharge.termp_charge());
    registry.get<TProfile>(HIST("termn_vs_cent"))->Fill(event_netcharge.centrality(), event_netcharge.termn_charge());
    registry.get<TProfile>(HIST("pos_sq_vs_cent"))->Fill(event_netcharge.centrality(), event_netcharge.possq_charge());
    registry.get<TProfile>(HIST("neg_sq_vs_cent"))->Fill(event_netcharge.centrality(), event_netcharge.negsq_charge());
    registry.get<TProfile>(HIST("posneg_vs_cent"))->Fill(event_netcharge.centrality(), event_netcharge.posneg_charge());
  }
}; // struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    {adaptAnalysisTask<NetchargeFluctuations>(cfgc)},
    {adaptAnalysisTask<Netcharge_analysis>(cfgc)},

  };
}
