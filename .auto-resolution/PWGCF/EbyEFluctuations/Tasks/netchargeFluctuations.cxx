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

/// \file netchargeFluctuations.cxx
/// \brief Calculate net-charge fluctuations using nu_dyn observable
///        For charged particles
///        For RUN-3
///
/// \author Nida Malik <nida.malik@cern.ch>
#include "PWGCF/Core/CorrelationContainer.h"
#include "PWGCF/Core/PairCuts.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/CCDB/TriggerAliases.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/FT0Corrected.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CommonConstants/MathConstants.h"
#include "CommonConstants/PhysicsConstants.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/runDataProcessing.h"

#include "TProfile.h"
#include "TProfile2D.h"
#include "TRandom3.h"

#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace std;
using namespace o2::constants::physics;

namespace o2
{
namespace aod
{
namespace net_charge
{
DECLARE_SOA_COLUMN(PosCharge, posCharge, float);
DECLARE_SOA_COLUMN(NegCharge, negCharge, float);
DECLARE_SOA_COLUMN(PosSqCharge, posSqCharge, float);
DECLARE_SOA_COLUMN(NegSqCharge, negSqCharge, float);
DECLARE_SOA_COLUMN(TermPCharge, termPCharge, float);
DECLARE_SOA_COLUMN(TermNCharge, termNCharge, float);
DECLARE_SOA_COLUMN(PosNegCharge, posNegCharge, float);
DECLARE_SOA_COLUMN(Centrality, centrality, float);
} // namespace net_charge

namespace net_charge_gen
{
DECLARE_SOA_COLUMN(PosCharge, posCharge, float);
DECLARE_SOA_COLUMN(NegCharge, negCharge, float);
DECLARE_SOA_COLUMN(PosSqCharge, posSqCharge, float);
DECLARE_SOA_COLUMN(NegSqCharge, negSqCharge, float);
DECLARE_SOA_COLUMN(TermPCharge, termPCharge, float);
DECLARE_SOA_COLUMN(TermNCharge, termNCharge, float);
DECLARE_SOA_COLUMN(PosNegCharge, posNegCharge, float);
DECLARE_SOA_COLUMN(Centrality, centrality, float);
} // namespace net_charge_gen

DECLARE_SOA_TABLE(NetCharge, "AOD", "NETChargefluct",
                  net_charge::PosCharge,
                  net_charge::NegCharge,
                  net_charge::PosSqCharge,
                  net_charge::NegSqCharge,
                  net_charge::TermPCharge,
                  net_charge::TermNCharge,
                  net_charge::PosNegCharge,
                  net_charge::Centrality);

DECLARE_SOA_TABLE(NetChargeGen, "AOD", "NETfluctGen",
                  net_charge_gen::PosCharge,
                  net_charge_gen::NegCharge,
                  net_charge_gen::PosSqCharge,
                  net_charge_gen::NegSqCharge,
                  net_charge_gen::TermPCharge,
                  net_charge_gen::TermNCharge,
                  net_charge_gen::PosNegCharge,
                  net_charge_gen::Centrality);

using MyCollisionsRun2 = soa::Join<aod::Collisions, aod::EvSels, aod::CentRun2V0Ms, aod::Mults>;
using MyCollisionRun2 = MyCollisionsRun2::iterator;
using MyCollisionsRun3 = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::CentFT0Cs, aod::Mults>;
using MyCollisionRun3 = MyCollisionsRun3::iterator;
using MyTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::StoredTracks, aod::TrackSelection>;
using MyTrack = MyTracks::iterator;

using MyMCCollisionsRun2 = soa::Join<aod::Collisions, aod::EvSels, aod::CentRun2V0Ms, aod::Mults, aod::McCollisionLabels>;
using MyMCCollisionRun2 = MyMCCollisionsRun2::iterator;

using MyMCCollisionsRun3 = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::CentFT0Cs, aod::Mults, aod::McCollisionLabels>;
using MyMCCollisionRun3 = MyMCCollisionsRun3::iterator;

using MyMCTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::StoredTracks, aod::TrackSelection, aod::McTrackLabels>;
using MyMCTrack = MyMCTracks::iterator;
} // namespace aod
} // namespace o2

enum RunType {
  kRun3 = 0,
  kRun2
};

struct NetchargeFluctuations {
  Produces<aod::NetCharge> netCharge;
  Produces<aod::NetChargeGen> netChargeGen;
  Service<o2::framework::O2DatabasePDG> pdgService;

  HistogramRegistry histogramRegistry{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // Configurables
  Configurable<float> vertexZcut{"vertexZcut", 10.f, "Vertex Z"};
  Configurable<float> etaCut{"etaCut", 0.8, "Eta cut"};
  Configurable<float> ptMinCut{"ptMinCut", 0.2, "Pt min cut"};
  Configurable<float> ptMaxCut{"ptMaxCut", 5.0, "Pt max cut"};
  Configurable<float> dcaXYCut{"dcaXYCut", 0.12, "DCA XY cut"};
  Configurable<float> dcaZCut{"dcaZCut", 0.3, "DCA Z cut"};
  Configurable<int> tpcCrossCut{"tpcCrossCut", 70, "TPC crossrows cut"};
  Configurable<int> itsChiCut{"itsChiCut", 70, "ITS chi2 cluster cut"};
  Configurable<int> tpcChiCut{"tpcChiCut", 70, "TPC chi2 cluster cut"};

  // Event selections
  Configurable<bool> cSel8Trig{"cSel8Trig", true, "Sel8 (T0A + T0C) Selection Run3"};    // sel8
  Configurable<bool> cInt7Trig{"cInt7Trig", true, "kINT7 MB Trigger"};                   // kINT7
  Configurable<bool> cSel7Trig{"cSel7Trig", true, "Sel7 (V0A + V0C) Selection Run2"};    // sel7
  Configurable<bool> cTFBorder{"cTFBorder", false, "Timeframe Border Selection"};        // pileup
  Configurable<bool> cNoItsROBorder{"cNoItsROBorder", false, "No ITSRO Border Cut"};     // pileup
  Configurable<bool> cItsTpcVtx{"cItsTpcVtx", false, "ITS+TPC Vertex Selection"};        // pileup
  Configurable<bool> cPileupReject{"cPileupReject", false, "Pileup rejection"};          // pileup
  Configurable<bool> cZVtxTimeDiff{"cZVtxTimeDiff", false, "z-vtx time diff selection"}; // pileup
  Configurable<bool> cIsGoodITSLayers{"cIsGoodITSLayers", false, "Good ITS Layers All"}; // pileup

  // Initialization
  void init(o2::framework::InitContext&)
  {
    const AxisSpec vtxzAxis = {800, -20, 20, "V_{Z} (cm)"};
    const AxisSpec dcaAxis = {250, -0.5, 0.5, "DCA_{xy} (cm)"};
    const AxisSpec dcazAxis = {250, -0.5, 0.5, "DCA_{z} (cm)"};
    const AxisSpec ptAxis = {70, 0.0, 7.0, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec etaAxis = {200, -1., 1., "#eta"};
    const AxisSpec cent1Axis = {10, 0., 100., "centrality"};
    const AxisSpec centAxis = {100, 0., 100., "centrality"};
    const AxisSpec multAxis = {200, 0., 10000., "FT0M Amplitude"};
    const AxisSpec tpcChiAxis = {1400, 0., 7., "Chi2"};
    const AxisSpec itsChiAxis = {800, 0., 40., "Chi2"};
    const AxisSpec crossedRowAxis = {1600, 0., 160., "TPC Crossed rows"};
    const AxisSpec eventsAxis = {10, 0, 10, ""};
    const AxisSpec signAxis = {20, -10, 10, ""};
    const AxisSpec nchAxis = {3000, 0, 3000, "Nch"};
    const AxisSpec nchpAxis = {50000, 0, 50000, "Nch"};

    histogramRegistry.add("hVtxZ_before", "", kTH1F, {vtxzAxis});
    histogramRegistry.add("hDcaXY_before", "", kTH1F, {dcaAxis});
    histogramRegistry.add("hDcaZ_before", "", kTH1F, {dcazAxis});
    histogramRegistry.add("hTPCchi2perCluster_before", "", kTH1D, {tpcChiAxis});
    histogramRegistry.add("hITSchi2perCluster_before", "", kTH1D, {itsChiAxis});
    histogramRegistry.add("hTPCCrossedrows_before", "", kTH1D, {crossedRowAxis});
    histogramRegistry.add("hPtDcaXY_before", "", kTH2D, {ptAxis, dcaAxis});
    histogramRegistry.add("hPtDcaZ_before", "", kTH2D, {ptAxis, dcazAxis});
    histogramRegistry.add("hVtxZ_after", "", kTH1F, {vtxzAxis});
    histogramRegistry.add("hDcaXY_after", "", kTH1F, {dcaAxis});
    histogramRegistry.add("hDcaZ_after", "", kTH1F, {dcazAxis});
    histogramRegistry.add("hTPCchi2perCluster_after", "", kTH1D, {tpcChiAxis});
    histogramRegistry.add("hITSchi2perCluster_after", "", kTH1D, {itsChiAxis});
    histogramRegistry.add("hTPCCrossedrows_after", "", kTH1D, {crossedRowAxis});
    histogramRegistry.add("hPtDcaXY_after", "", kTH2D, {ptAxis, dcaAxis});
    histogramRegistry.add("hPtDcaZ_after", "", kTH2D, {ptAxis, dcazAxis});
    histogramRegistry.add("hEta", "", kTH1F, {etaAxis});
    histogramRegistry.add("hPt", "", kTH1F, {ptAxis});
    histogramRegistry.add("hCentrality", "", kTH1F, {centAxis});
    histogramRegistry.add("hMultiplicity", "", kTH1F, {multAxis});
    histogramRegistry.add("gen_hVtxZ_before", "", kTH1F, {vtxzAxis});
    histogramRegistry.add("gen_hVtxZ_after", "", kTH1F, {vtxzAxis});
    histogramRegistry.add("gen_hEta", "", kTH1F, {etaAxis});
    histogramRegistry.add("gen_hSign", "", kTH1F, {signAxis});
    histogramRegistry.add("gen_hPt", "", kTH1F, {ptAxis});

    histogramRegistry.add("nch", "", kTH1D, {nchAxis});
    histogramRegistry.add("nch0", "", kTH1D, {nchAxis});
    histogramRegistry.add("nch1", "", kTH1D, {nchAxis});
    histogramRegistry.add("nch2", "", kTH1D, {nchAxis});
    histogramRegistry.add("nch3", "", kTH1D, {nchAxis});
    histogramRegistry.add("nch4", "", kTH1D, {nchAxis});
    histogramRegistry.add("nch5", "", kTH1D, {nchAxis});
    histogramRegistry.add("nch6", "", kTH1D, {nchAxis});
    histogramRegistry.add("nch7", "", kTH1D, {nchAxis});
    histogramRegistry.add("nch8", "", kTH1D, {nchAxis});
    histogramRegistry.add("nch9", "", kTH1D, {nchAxis});
    histogramRegistry.add("nch_pos", "", kTH1D, {nchAxis});
    histogramRegistry.add("nch_neg", "", kTH1D, {nchAxis});
    histogramRegistry.add("nch_negpos", "", kTH1D, {nchpAxis});

    histogramRegistry.add("nch_cent", "", kTH2D, {centAxis, nchAxis});
    histogramRegistry.add("nch_pos_cent", "", kTH2D, {centAxis, nchAxis});
    histogramRegistry.add("nch_neg_cent", "", kTH2D, {centAxis, nchAxis});
    histogramRegistry.add("nch_negpos_cent", "", kTH2D, {centAxis, nchpAxis});
  }
  template <RunType run, typename C>
  bool selCollision(C const& coll, float& cent, float& mult)
  {
    histogramRegistry.fill(HIST("hVtxZ_before"), coll.posZ());

    if (std::abs(coll.posZ()) > vertexZcut) {
      return false;
    }

    if constexpr (run == kRun3) {

      if (cSel8Trig && !coll.sel8()) {
        return false;
      } // require min bias trigger
      cent = coll.centFT0M(); // centrality for run3
      mult = coll.multFT0M(); // multiplicity for run3
    } else {
      if (cInt7Trig && !coll.alias_bit(kINT7)) {
        return false;
      }
      if (cSel7Trig && !coll.sel7()) {
        return false;
      }
      cent = coll.centRun2V0M(); // centrality for run2
      mult = coll.multFV0M();    // multiplicity for run2
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
    histogramRegistry.fill(HIST("hVtxZ_after"), coll.posZ());
    histogramRegistry.fill(HIST("hCentrality"), cent);
    histogramRegistry.fill(HIST("hMultiplicity"), mult);

    return true;
  }
  template <typename T>
  void fillBeforeQA(T const& track)
  {
    histogramRegistry.fill(HIST("hTPCchi2perCluster_before"), track.tpcChi2NCl());
    histogramRegistry.fill(HIST("hITSchi2perCluster_before"), track.itsChi2NCl());
    histogramRegistry.fill(HIST("hTPCCrossedrows_before"), track.tpcNClsCrossedRows());
    histogramRegistry.fill(HIST("hDcaXY_before"), track.dcaXY());
    histogramRegistry.fill(HIST("hDcaZ_before"), track.dcaZ());
    histogramRegistry.fill(HIST("hPtDcaXY_before"), track.pt(), track.dcaXY());
    histogramRegistry.fill(HIST("hPtDcaZ_before"), track.pt(), track.dcaZ());
  }

  template <typename T>
  void fillAfterQA(T const& track)
  {
    histogramRegistry.fill(HIST("hDcaXY_after"), track.dcaXY());
    histogramRegistry.fill(HIST("hDcaZ_after"), track.dcaZ());
    histogramRegistry.fill(HIST("hPt"), track.pt());
    histogramRegistry.fill(HIST("hEta"), track.eta());
    histogramRegistry.fill(HIST("hPtDcaXY_after"), track.pt(), track.dcaXY());
    histogramRegistry.fill(HIST("hPtDcaZ_after"), track.pt(), track.dcaZ());
    histogramRegistry.fill(HIST("hTPCCrossedrows_after"), track.tpcNClsCrossedRows());
    histogramRegistry.fill(HIST("hTPCchi2perCluster_after"), track.tpcChi2NCl());
    histogramRegistry.fill(HIST("hITSchi2perCluster_after"), track.itsChi2NCl());
  }

  template <typename N, typename C>
  void nchDistribution(N const& nch, C const& cent)
  {
    constexpr int kCent0Min = 0, kCent0Max = 5;
    constexpr int kCent1Min = 5, kCent1Max = 10;
    constexpr int kCent2Min = 10, kCent2Max = 20;
    constexpr int kCent3Min = 20, kCent3Max = 30;
    constexpr int kCent4Min = 30, kCent4Max = 40;
    constexpr int kCent5Min = 40, kCent5Max = 50;
    constexpr int kCent6Min = 50, kCent6Max = 60;
    constexpr int kCent7Min = 60, kCent7Max = 70;
    constexpr int kCent8Min = 70, kCent8Max = 80;
    constexpr int kCent9Min = 80, kCent9Max = 90;

    histogramRegistry.fill(HIST("nch"), nch);
    histogramRegistry.fill(HIST("nch_cent"), cent, nch);

    if (cent >= kCent0Min && cent < kCent0Max)
      histogramRegistry.fill(HIST("nch0"), nch);
    else if (cent >= kCent1Min && cent < kCent1Max)
      histogramRegistry.fill(HIST("nch1"), nch);
    else if (cent >= kCent2Min && cent < kCent2Max)
      histogramRegistry.fill(HIST("nch2"), nch);
    else if (cent >= kCent3Min && cent < kCent3Max)
      histogramRegistry.fill(HIST("nch3"), nch);
    else if (cent >= kCent4Min && cent < kCent4Max)
      histogramRegistry.fill(HIST("nch4"), nch);
    else if (cent >= kCent5Min && cent < kCent5Max)
      histogramRegistry.fill(HIST("nch5"), nch);
    else if (cent >= kCent6Min && cent < kCent6Max)
      histogramRegistry.fill(HIST("nch6"), nch);
    else if (cent >= kCent7Min && cent < kCent7Max)
      histogramRegistry.fill(HIST("nch7"), nch);
    else if (cent >= kCent8Min && cent < kCent8Max)
      histogramRegistry.fill(HIST("nch8"), nch);
    else if (cent >= kCent9Min && cent < kCent9Max)
      histogramRegistry.fill(HIST("nch9"), nch);
  }

  template <typename T>
  bool selTrack(T const& track)
  {
    if (!track.isGlobalTrack()) {
      return false;
    } // accept only global tracks

    if (std::fabs(track.eta()) >= etaCut) {
      return false;
    }

    if (track.pt() <= ptMinCut || track.pt() >= ptMaxCut) {
      return false;
    }

    if (track.sign() == 0) {
      return false;
    }
    /*    if (std::fabs(track.dcaXY()) > dcaXYCut) {
          return false;
        }

        if (std::fabs(track.dcaZ()) > dcaZCut) {
          return false ;
        }*/

    /*
    if (track.tpcNClsCrossedRows() < tpcCrossCut) {
      return false;
    }

    if (track.itsChi2NCl() > itsChiCut) {
      return false;
    }

    if (track.tpcChi2NCl() > tpcChiCut) {
      return false;
    }

*/

    return true;
  }

  template <RunType run, typename C, typename T>
  void calculationData(C const& coll, T const& tracks)
  {
    float cent = -1, mult = -1;
    if (!selCollision<run>(coll, cent, mult)) {
      return;
    }
    int fpos = 0, fneg = 0, posneg = 0, termn = 0, termp = 0;
    int nch = 0;
    for (const auto& track : tracks) {
      fillBeforeQA(track);

      if (!selTrack(track)) {
        continue;
      }
      nch += 1;

      fillAfterQA(track);
      if (track.sign() == 1) {
        fpos += 1;
      } else if (track.sign() == -1) {
        fneg += 1;
      }
    } // track
    termp = fpos * (fpos - 1);
    termn = fneg * (fneg - 1);
    posneg = fpos * fneg;

    nchDistribution(nch, cent);

    histogramRegistry.fill(HIST("nch_pos"), fpos);
    histogramRegistry.fill(HIST("nch_pos_cent"), cent, fpos);
    histogramRegistry.fill(HIST("nch_neg"), fneg);
    histogramRegistry.fill(HIST("nch_neg_cent"), cent, fneg);
    histogramRegistry.fill(HIST("nch_negpos"), posneg);
    histogramRegistry.fill(HIST("nch_negpos_cent"), cent, posneg);

    netCharge(fpos, fneg, fpos * fpos, fneg * fneg, termp, termn, posneg, cent);
  }

  template <RunType run, typename C, typename T, typename M, typename P>
  void calculationMc(C const& coll, T const& inputTracks, M const& mcCollisions, P const& mcParticles)
  {
    (void)mcCollisions;

    if (!coll.has_mcCollision()) {
      return;
    }
    histogramRegistry.fill(HIST("gen_hVtxZ_before"), coll.mcCollision().posZ());

    float cent = -1, mult = -1;

    if (!selCollision<run>(coll, cent, mult)) {
      return;
    }

    int fpos = 0, fneg = 0, posneg = 0, termn = 0, termp = 0;
    int nch = 0;

    for (const auto& track : inputTracks) {

      fillBeforeQA(track);

      if (!selTrack(track)) {
        continue;
      }
      nch += 1;

      fillAfterQA(track);

      if (track.sign() == 1) {
        fpos += 1;
      } else if (track.sign() == -1) {
        fneg += 1;
      }
    } // track
    termp = fpos * (fpos - 1);
    termn = fneg * (fneg - 1);
    posneg = fpos * fneg;

    nchDistribution(nch, cent);

    histogramRegistry.fill(HIST("nch_pos"), fpos);
    histogramRegistry.fill(HIST("nch_pos_cent"), cent, fpos);
    histogramRegistry.fill(HIST("nch_neg"), fneg);
    histogramRegistry.fill(HIST("nch_neg_cent"), cent, fneg);
    histogramRegistry.fill(HIST("nch_negpos"), posneg);
    histogramRegistry.fill(HIST("nch_negpos_cent"), cent, posneg);

    netCharge(fpos, fneg, fpos * fpos, fneg * fneg, termp, termn, posneg, cent);

    int posGen = 0, negGen = 0, posNegGen = 0, termNGen = 0, termPGen = 0;

    const auto& mccolgen = coll.template mcCollision_as<aod::McCollisions>();

    /* if (std::abs(mccolgen.posZ()) > vertexZcut) {
        return;
      }*/

    const auto& mcpartgen = mcParticles.sliceByCached(aod::mcparticle::mcCollisionId, mccolgen.globalIndex(), cache);
    histogramRegistry.fill(HIST("gen_hVtxZ_after"), mccolgen.posZ());

    for (const auto& mcpart : mcpartgen) {

      if (!mcpart.isPhysicalPrimary()) {
        continue;
      }
      int pid = mcpart.pdgCode();
      auto sign = 0;
      auto* pd = pdgService->GetParticle(pid);
      if (pd != nullptr) {
        sign = pd->Charge() / 3.;
      }

      if (sign == 0) {
        continue;
      }

      if (std::abs(pid) != kElectron && std::abs(pid) != kMuonMinus && std::abs(pid) != kPiPlus && std::abs(pid) != kKPlus && std::abs(pid) != kProton) {
        continue;
      }

      if (std::fabs(mcpart.eta()) > etaCut)
        continue;

      if ((mcpart.pt() <= ptMinCut) || (mcpart.pt() >= ptMaxCut))
        continue;

      histogramRegistry.fill(HIST("gen_hPt"), mcpart.pt());
      histogramRegistry.fill(HIST("gen_hEta"), mcpart.eta());
      histogramRegistry.fill(HIST("gen_hSign"), sign);

      if (sign == 1) {
        posGen += 1;
      }

      if (sign == -1) {
        negGen += 1;
      }
    } // particle
    termPGen = posGen * (posGen - 1);
    termNGen = negGen * (negGen - 1);
    posNegGen = posGen * negGen;

    netChargeGen(posGen, negGen, posGen * posGen, negGen * negGen,
                 termPGen, termNGen, posNegGen, cent);

  } // void

  SliceCache cache;
  Preslice<aod::McParticles> mcTrack = o2::aod::mcparticle::mcCollisionId;

  void processDataRun3(aod::MyCollisionRun3 const& coll, aod::MyTracks const& tracks)
  {
    calculationData<kRun3>(coll, tracks);
  }

  PROCESS_SWITCH(NetchargeFluctuations, processDataRun3, "Process for Run3 DATA", true);

  void processDataRun2(aod::MyCollisionRun2 const& coll, aod::MyTracks const& tracks)
  {
    calculationData<kRun2>(coll, tracks);
  }

  PROCESS_SWITCH(NetchargeFluctuations, processDataRun2, "Process for Run2 DATA", false);

  void processMcRun3(aod::MyMCCollisionRun3 const& coll, aod::MyMCTracks const& inputTracks,
                     aod::McCollisions const& mcCollisions, aod::McParticles const& mcParticles)
  {
    calculationMc<kRun3>(coll, inputTracks, mcCollisions, mcParticles);
  }

  PROCESS_SWITCH(NetchargeFluctuations, processMcRun3, "Process reconstructed", false);

  void processMcRun2(aod::MyMCCollisionRun2 const& coll, aod::MyMCTracks const& inputTracks,
                     aod::McCollisions const& mcCollisions, aod::McParticles const& mcParticles)
  {
    calculationMc<kRun2>(coll, inputTracks, mcCollisions, mcParticles);
  }

  PROCESS_SWITCH(NetchargeFluctuations, processMcRun2, "Process reconstructed", false);

}; // struct

struct NetchargeAnalysis {
  Configurable<int> cfSubSample{"cfSubSample", 30, "Number of subsamples"};
  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};
  std::vector<std::vector<std::shared_ptr<TProfile>>> net;
  std::vector<std::vector<std::shared_ptr<TProfile>>> subSample;
  std::vector<std::vector<std::shared_ptr<TProfile>>> genSubSample;
  TRandom3* fRndm = new TRandom3(0);

  void init(o2::framework::InitContext&)
  {
    std::vector<double> centBinning = {0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100};
    AxisSpec centAxis = {centBinning, "centrality"};

    registry.add("data/pos_vs_cent", "", {HistType::kTProfile, {centAxis}});
    registry.add("data/neg_vs_cent", "", {HistType::kTProfile, {centAxis}});
    registry.add("data/termp_vs_cent", "", {HistType::kTProfile, {centAxis}});
    registry.add("data/termn_vs_cent", "", {HistType::kTProfile, {centAxis}});
    registry.add("data/pos_sq_vs_cent", "", {HistType::kTProfile, {centAxis}});
    registry.add("data/neg_sq_vs_cent", "", {HistType::kTProfile, {centAxis}});
    registry.add("data/posneg_vs_cent", "", {HistType::kTProfile, {centAxis}});

    registry.add("gen/pos_vs_cent", "", {HistType::kTProfile, {centAxis}});
    registry.add("gen/neg_vs_cent", "", {HistType::kTProfile, {centAxis}});
    registry.add("gen/termp_vs_cent", "", {HistType::kTProfile, {centAxis}});
    registry.add("gen/termn_vs_cent", "", {HistType::kTProfile, {centAxis}});
    registry.add("gen/pos_sq_vs_cent", "", {HistType::kTProfile, {centAxis}});
    registry.add("gen/neg_sq_vs_cent", "", {HistType::kTProfile, {centAxis}});
    registry.add("gen/posneg_vs_cent", "", {HistType::kTProfile, {centAxis}});

    subSample.resize(cfSubSample);
    genSubSample.resize(cfSubSample);

    for (int i = 0; i < cfSubSample; ++i) {
      subSample[i].resize(7);
      genSubSample[i].resize(7);

      subSample[i][0] = std::get<std::shared_ptr<TProfile>>(registry.add(Form("data/subSample_%d/pos_vs_cent", i), "", {HistType::kTProfile, {centAxis}}));
      subSample[i][1] = std::get<std::shared_ptr<TProfile>>(registry.add(Form("data/subSample_%d/neg_vs_cent", i), "", {HistType::kTProfile, {centAxis}}));
      subSample[i][2] = std::get<std::shared_ptr<TProfile>>(registry.add(Form("data/subSample_%d/termp_vs_cent", i), "", {HistType::kTProfile, {centAxis}}));
      subSample[i][3] = std::get<std::shared_ptr<TProfile>>(registry.add(Form("data/subSample_%d/termn_vs_cent", i), "", {HistType::kTProfile, {centAxis}}));
      subSample[i][4] = std::get<std::shared_ptr<TProfile>>(registry.add(Form("data/subSample_%d/pos_sq_vs_cent", i), "", {HistType::kTProfile, {centAxis}}));
      subSample[i][5] = std::get<std::shared_ptr<TProfile>>(registry.add(Form("data/subSample_%d/neg_sq_vs_cent", i), "", {HistType::kTProfile, {centAxis}}));
      subSample[i][6] = std::get<std::shared_ptr<TProfile>>(registry.add(Form("data/subSample_%d/posneg_vs_cent", i), "", {HistType::kTProfile, {centAxis}}));

      genSubSample[i][0] = std::get<std::shared_ptr<TProfile>>(registry.add(Form("gen/genSubSample_%d/pos_vs_cent", i), "", {HistType::kTProfile, {centAxis}}));
      genSubSample[i][1] = std::get<std::shared_ptr<TProfile>>(registry.add(Form("gen/genSubSample_%d/neg_vs_cent", i), "", {HistType::kTProfile, {centAxis}}));
      genSubSample[i][2] = std::get<std::shared_ptr<TProfile>>(registry.add(Form("gen/genSubSample_%d/termp_vs_cent", i), "", {HistType::kTProfile, {centAxis}}));
      genSubSample[i][3] = std::get<std::shared_ptr<TProfile>>(registry.add(Form("gen/genSubSample_%d/termn_vs_cent", i), "", {HistType::kTProfile, {centAxis}}));
      genSubSample[i][4] = std::get<std::shared_ptr<TProfile>>(registry.add(Form("gen/genSubSample_%d/pos_sq_vs_cent", i), "", {HistType::kTProfile, {centAxis}}));
      genSubSample[i][5] = std::get<std::shared_ptr<TProfile>>(registry.add(Form("gen/genSubSample_%d/neg_sq_vs_cent", i), "", {HistType::kTProfile, {centAxis}}));
      genSubSample[i][6] = std::get<std::shared_ptr<TProfile>>(registry.add(Form("gen/genSubSample_%d/posneg_vs_cent", i), "", {HistType::kTProfile, {centAxis}}));
    }

  } // void

  void processData(aod::NetCharge::iterator const& event_netcharge)
  {
    registry.get<TProfile>(HIST("data/pos_vs_cent"))->Fill(event_netcharge.centrality(), event_netcharge.posCharge());
    registry.get<TProfile>(HIST("data/neg_vs_cent"))->Fill(event_netcharge.centrality(), event_netcharge.negCharge());
    registry.get<TProfile>(HIST("data/termp_vs_cent"))->Fill(event_netcharge.centrality(), event_netcharge.termPCharge());
    registry.get<TProfile>(HIST("data/termn_vs_cent"))->Fill(event_netcharge.centrality(), event_netcharge.termNCharge());
    registry.get<TProfile>(HIST("data/pos_sq_vs_cent"))->Fill(event_netcharge.centrality(), event_netcharge.posSqCharge());
    registry.get<TProfile>(HIST("data/neg_sq_vs_cent"))->Fill(event_netcharge.centrality(), event_netcharge.negSqCharge());
    registry.get<TProfile>(HIST("data/posneg_vs_cent"))->Fill(event_netcharge.centrality(), event_netcharge.posNegCharge());

    int sampleIndex = static_cast<int>(cfSubSample * fRndm->Rndm());
    subSample[sampleIndex][0]->Fill(event_netcharge.centrality(), event_netcharge.posCharge());
    subSample[sampleIndex][1]->Fill(event_netcharge.centrality(), event_netcharge.negCharge());
    subSample[sampleIndex][2]->Fill(event_netcharge.centrality(), event_netcharge.termPCharge());
    subSample[sampleIndex][3]->Fill(event_netcharge.centrality(), event_netcharge.termNCharge());
    subSample[sampleIndex][4]->Fill(event_netcharge.centrality(), event_netcharge.posSqCharge());
    subSample[sampleIndex][5]->Fill(event_netcharge.centrality(), event_netcharge.negSqCharge());
    subSample[sampleIndex][6]->Fill(event_netcharge.centrality(), event_netcharge.posNegCharge());
  } // void
  PROCESS_SWITCH(NetchargeAnalysis, processData, "Process reconstructed and Data", true);

  void processGen(aod::NetChargeGen::iterator const& event_netcharge)
  {
    registry.get<TProfile>(HIST("gen/pos_vs_cent"))->Fill(event_netcharge.centrality(), event_netcharge.posCharge());
    registry.get<TProfile>(HIST("gen/neg_vs_cent"))->Fill(event_netcharge.centrality(), event_netcharge.negCharge());
    registry.get<TProfile>(HIST("gen/termp_vs_cent"))->Fill(event_netcharge.centrality(), event_netcharge.termPCharge());
    registry.get<TProfile>(HIST("gen/termn_vs_cent"))->Fill(event_netcharge.centrality(), event_netcharge.termNCharge());
    registry.get<TProfile>(HIST("gen/pos_sq_vs_cent"))->Fill(event_netcharge.centrality(), event_netcharge.posSqCharge());
    registry.get<TProfile>(HIST("gen/neg_sq_vs_cent"))->Fill(event_netcharge.centrality(), event_netcharge.negSqCharge());
    registry.get<TProfile>(HIST("gen/posneg_vs_cent"))->Fill(event_netcharge.centrality(), event_netcharge.posNegCharge());

    int sampleIndex = static_cast<int>(cfSubSample * fRndm->Rndm());
    genSubSample[sampleIndex][0]->Fill(event_netcharge.centrality(), event_netcharge.posCharge());
    genSubSample[sampleIndex][1]->Fill(event_netcharge.centrality(), event_netcharge.negCharge());
    genSubSample[sampleIndex][2]->Fill(event_netcharge.centrality(), event_netcharge.termPCharge());
    genSubSample[sampleIndex][3]->Fill(event_netcharge.centrality(), event_netcharge.termNCharge());
    genSubSample[sampleIndex][4]->Fill(event_netcharge.centrality(), event_netcharge.posSqCharge());
    genSubSample[sampleIndex][5]->Fill(event_netcharge.centrality(), event_netcharge.negSqCharge());
    genSubSample[sampleIndex][6]->Fill(event_netcharge.centrality(), event_netcharge.posNegCharge());
  } // void
  PROCESS_SWITCH(NetchargeAnalysis, processGen, "Process generated", true);

}; // struct Netcharge_analysis

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    {adaptAnalysisTask<NetchargeFluctuations>(cfgc)},
    {adaptAnalysisTask<NetchargeAnalysis>(cfgc)}

  };
}
