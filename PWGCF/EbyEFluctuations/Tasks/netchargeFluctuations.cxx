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

#include "CCDB/BasicCCDBManager.h"
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

#include <string>
#include <vector> // Include for std::vector

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace std;
using namespace o2::constants::physics;

namespace o2
{
namespace aod
{
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
  Service<o2::framework::O2DatabasePDG> pdgService;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  TRandom3* fRndm = new TRandom3(0);
  HistogramRegistry histogramRegistry{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // Configurables
  Configurable<int64_t> ccdbNoLaterThan{"ccdbNoLaterThan", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};
  Configurable<std::string> cfgUrlCCDB{"cfgUrlCCDB", "http://alice-ccdb.cern.ch", "url of ccdb"};
  Configurable<std::string> cfgPathCCDB{"cfgPathCCDB", "Users/n/nimalik/efftest", "Path for ccdb-object"};
  Configurable<bool> cfgLoadEff{"cfgLoadEff", true, "Load efficiency"};

  Configurable<float> vertexZcut{"vertexZcut", 10.f, "Vertex Z"};
  Configurable<float> etaCut{"etaCut", 0.8, "Eta cut"};
  Configurable<float> ptMinCut{"ptMinCut", 0.2, "Pt min cut"};
  Configurable<float> ptMaxCut{"ptMaxCut", 5.0, "Pt max cut"};
  Configurable<float> dcaXYCut{"dcaXYCut", 0.12, "DCA XY cut"};
  Configurable<float> dcaZCut{"dcaZCut", 0.3, "DCA Z cut"};
  Configurable<int> tpcCrossCut{"tpcCrossCut", 70, "TPC crossrows cut"};
  Configurable<int> itsChiCut{"itsChiCut", 70, "ITS chi2 cluster cut"};
  Configurable<int> tpcChiCut{"tpcChiCut", 70, "TPC chi2 cluster cut"};
  Configurable<float> centMin{"centMin", 0.0f, "cenrality min for delta eta"};
  Configurable<float> centMax{"centMax", 10.0f, "cenrality max for delta eta"};
  Configurable<int> cfgNSubsample{"cfgNSubsample", 30, "Number of subsamples for Error"};
  Configurable<int> deltaEta{"deltaEta", 8, "Delta eta bin count"};
  Configurable<double> threshold{"threshold", 1e-6, "Delta eta bin count"};
  // Event selections
  Configurable<bool> cSel8Trig{"cSel8Trig", true, "Sel8 (T0A + T0C) Selection Run3"};                    // sel8
  Configurable<bool> cInt7Trig{"cInt7Trig", true, "kINT7 MB Trigger"};                                   // kINT7
  Configurable<bool> cSel7Trig{"cSel7Trig", true, "Sel7 (V0A + V0C) Selection Run2"};                    // sel7
  Configurable<bool> cTFBorder{"cTFBorder", false, "Timeframe Border Selection"};                        // pileup
  Configurable<bool> cNoItsROBorder{"cNoItsROBorder", false, "No ITSRO Border Cut"};                     // pileup
  Configurable<bool> cItsTpcVtx{"cItsTpcVtx", false, "ITS+TPC Vertex Selection"};                        // pileup
  Configurable<bool> cPileupReject{"cPileupReject", false, "Pileup rejection"};                          // pileup
  Configurable<bool> cZVtxTimeDiff{"cZVtxTimeDiff", false, "z-vtx time diff selection"};                 // pileup
  Configurable<bool> cfgUseGoodItsLayerAllCut{"cfgUseGoodItsLayerAllCut", false, "Good ITS Layers All"}; // pileup
  Configurable<bool> cDcaXy{"cDcaXy", false, "Dca XY cut"};
  Configurable<bool> cDcaZ{"cDcaZ", false, "Dca Z cut"};
  Configurable<bool> cTpcCr{"cTpcCr", false, "tpc crossrows"};
  Configurable<bool> cItsChi{"cItsChi", false, "ITS chi"};
  Configurable<bool> cTpcChi{"cTpcChi", false, "TPC chi"};

  // CCDB efficiency histograms
  TH2D* efficiency = nullptr;

  // Initialization
  void init(o2::framework::InitContext&)
  {
    const AxisSpec vtxzAxis = {800, -20, 20, "V_{Z} (cm)"};
    const AxisSpec dcaAxis = {250, -0.5, 0.5, "DCA_{xy} (cm)"};
    const AxisSpec dcazAxis = {250, -0.5, 0.5, "DCA_{z} (cm)"};
    const AxisSpec ptAxis = {70, 0.0, 7.0, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec etaAxis = {20, -1., 1., "#eta"};
    const AxisSpec deltaEtaAxis = {9, 0, 1.8, "#eta"};
    const AxisSpec centAxis = {100, 0., 100., "centrality"};
    const AxisSpec multAxis = {200, 0., 10000., "FT0M Amplitude"};
    const AxisSpec tpcChiAxis = {1400, 0., 7., "Chi2"};
    const AxisSpec itsChiAxis = {800, 0., 40., "Chi2"};
    const AxisSpec crossedRowAxis = {1600, 0., 160., "TPC Crossed rows"};
    const AxisSpec eventsAxis = {10, 0, 10, ""};
    const AxisSpec signAxis = {20, -10, 10, ""};
    const AxisSpec nchAxis = {5000, 0, 5000, "Nch"};
    const AxisSpec nch1Axis = {1500, 0, 1500, "Nch"};
    const AxisSpec nchpAxis = {50000, 0, 50000, "Nch"};

    std::vector<double> centBining = {0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100};
    AxisSpec cent1Axis = {centBining, "Multiplicity percentile from FT0M (%)"};

    auto noSubsample = static_cast<int>(cfgNSubsample);
    float maxSubsample = 1.0 * noSubsample;
    AxisSpec subsampleAxis = {noSubsample, 0.0, maxSubsample, "subsample no."};

    histogramRegistry.add("data/hVtxZ_before", "", kTH1F, {vtxzAxis});
    histogramRegistry.add("data/hDcaXY_before", "", kTH1F, {dcaAxis});
    histogramRegistry.add("data/hDcaZ_before", "", kTH1F, {dcazAxis});
    histogramRegistry.add("data/hTPCchi2perCluster_before", "", kTH1D, {tpcChiAxis});
    histogramRegistry.add("data/hITSchi2perCluster_before", "", kTH1D, {itsChiAxis});
    histogramRegistry.add("data/hTPCCrossedrows_before", "", kTH1D, {crossedRowAxis});
    histogramRegistry.add("data/hPtDcaXY_before", "", kTH2D, {ptAxis, dcaAxis});
    histogramRegistry.add("data/hPtDcaZ_before", "", kTH2D, {ptAxis, dcazAxis});
    histogramRegistry.add("data/hVtxZ_after", "", kTH1F, {vtxzAxis});
    histogramRegistry.add("data/hDcaXY_after", "", kTH1F, {dcaAxis});
    histogramRegistry.add("data/hDcaZ_after", "", kTH1F, {dcazAxis});
    histogramRegistry.add("data/hTPCchi2perCluster_after", "", kTH1D, {tpcChiAxis});
    histogramRegistry.add("data/hITSchi2perCluster_after", "", kTH1D, {itsChiAxis});
    histogramRegistry.add("data/hTPCCrossedrows_after", "", kTH1D, {crossedRowAxis});
    histogramRegistry.add("data/hPtDcaXY_after", "", kTH2D, {ptAxis, dcaAxis});
    histogramRegistry.add("data/hPtDcaZ_after", "", kTH2D, {ptAxis, dcazAxis});
    histogramRegistry.add("data/hEta", "", kTH1F, {etaAxis});
    histogramRegistry.add("data/hEta_cent", "", kTH2F, {cent1Axis, etaAxis});
    histogramRegistry.add("data/hPt", "", kTH1F, {ptAxis});
    histogramRegistry.add("data/hPt_cent", "", kTH2F, {cent1Axis, ptAxis});
    histogramRegistry.add("data/hPt_eta", "", kTH2F, {ptAxis, etaAxis});
    histogramRegistry.add("data/hCentrality", "", kTH1F, {centAxis});
    histogramRegistry.add("data/hMultiplicity", "", kTH1F, {multAxis});

    histogramRegistry.add("gen/hPt_eta", "", kTH2F, {ptAxis, etaAxis});
    histogramRegistry.add("gen/hVtxZ_before", "", kTH1F, {vtxzAxis});
    histogramRegistry.add("gen/hVtxZ_after", "", kTH1F, {vtxzAxis});
    histogramRegistry.add("gen/hEta", "", kTH1F, {etaAxis});
    histogramRegistry.add("gen/hEta_cent", "", kTH2F, {centAxis, etaAxis});
    histogramRegistry.add("gen/hSign", "", kTH1F, {signAxis});
    histogramRegistry.add("gen/hPt", "", kTH1F, {ptAxis});
    histogramRegistry.add("gen/hPt_cent", "", kTH2F, {centAxis, ptAxis});
    histogramRegistry.add("gen/nch", "", kTH1F, {nchAxis});

    histogramRegistry.add("mult_dist/nch", "", kTH1D, {nchAxis});
    histogramRegistry.add("mult_dist/nch_pos", "", kTH1D, {nchAxis});
    histogramRegistry.add("mult_dist/nch_neg", "", kTH1D, {nchAxis});
    histogramRegistry.add("mult_dist/nch_negpos", "", kTH1D, {nchpAxis});
    histogramRegistry.add("mult_dist/nch_cent", "", kTH2D, {centAxis, nchAxis});
    histogramRegistry.add("mult_dist/nch_pos_cent", "", kTH2D, {centAxis, nchAxis});
    histogramRegistry.add("mult_dist/nch_neg_cent", "", kTH2D, {centAxis, nchAxis});
    histogramRegistry.add("mult_dist/nch_negpos_cent", "", kTH2D, {centAxis, nchpAxis});

    histogramRegistry.add("delta_eta/cent", "Centrality", kTH1F, {cent1Axis});
    histogramRegistry.add("delta_eta/track_eta", "eta", kTH1F, {etaAxis});
    histogramRegistry.add("delta_eta/pos", "delta_eta vs fpos", kTProfile, {deltaEtaAxis});
    histogramRegistry.add("delta_eta/neg", "delta_eta vs fneg", kTProfile, {deltaEtaAxis});
    histogramRegistry.add("delta_eta/termp", "delta_eta vs termp", kTProfile, {deltaEtaAxis});
    histogramRegistry.add("delta_eta/termn", "delta_eta vs termn", kTProfile, {deltaEtaAxis});
    histogramRegistry.add("delta_eta/pos_sq", "delta_eta vs sqfpos", kTProfile, {deltaEtaAxis});
    histogramRegistry.add("delta_eta/neg_sq", "delta_eta vs sqfneg", kTProfile, {deltaEtaAxis});
    histogramRegistry.add("delta_eta/posneg", "delta_eta vs fpos*fneg", kTProfile, {deltaEtaAxis});

    histogramRegistry.add("cent/pos", "cent vs fpos", kTProfile, {cent1Axis});
    histogramRegistry.add("cent/neg", "cent vs fneg", kTProfile, {cent1Axis});
    histogramRegistry.add("cent/termp", "cent vs termp", kTProfile, {cent1Axis});
    histogramRegistry.add("cent/termn", "cent vs termn", kTProfile, {cent1Axis});
    histogramRegistry.add("cent/pos_sq", "cent vs sqfpos", kTProfile, {cent1Axis});
    histogramRegistry.add("cent/neg_sq", "cent vs sqfneg", kTProfile, {cent1Axis});
    histogramRegistry.add("cent/posneg", "cent vs fpos*fneg", kTProfile, {cent1Axis});

    histogramRegistry.add("cent/gen_pos", "cent vs fpos", kTProfile, {cent1Axis});
    histogramRegistry.add("cent/gen_neg", "cent vs fneg", kTProfile, {cent1Axis});
    histogramRegistry.add("cent/gen_termp", "cent vs termp", kTProfile, {cent1Axis});
    histogramRegistry.add("cent/gen_termn", "cent vs termn", kTProfile, {cent1Axis});
    histogramRegistry.add("cent/gen_pos_sq", "cent vs sqfpos", kTProfile, {cent1Axis});
    histogramRegistry.add("cent/gen_neg_sq", "cent vs sqfneg", kTProfile, {cent1Axis});
    histogramRegistry.add("cent/gen_posneg", "cent vs fpos*fneg", kTProfile, {cent1Axis});
    histogramRegistry.add("cent/gen_nch", "cent vs nch", kTProfile, {centAxis});

    histogramRegistry.add("cor/hPt_cor", "", kTH1F, {ptAxis});
    histogramRegistry.add("cor/hEta_cor", "", kTH1F, {etaAxis});
    histogramRegistry.add("cor/nch_vs_nchCor", "", kTProfile, {nchAxis});
    histogramRegistry.add("cor/nchCor", "", kTH1F, {nchAxis});
    histogramRegistry.add("cor/cent_nchCor", "", kTH2F, {centAxis, nchAxis});
    histogramRegistry.add("cor/fpos_cent", "", kTProfile, {centAxis});
    histogramRegistry.add("cor/fneg_cent", "", kTProfile, {centAxis});

    histogramRegistry.add("subsample/pos", "", kTProfile2D, {cent1Axis, subsampleAxis});
    histogramRegistry.add("subsample/neg", "", kTProfile2D, {cent1Axis, subsampleAxis});
    histogramRegistry.add("subsample/termp", "", kTProfile2D, {cent1Axis, subsampleAxis});
    histogramRegistry.add("subsample/termn", "", kTProfile2D, {cent1Axis, subsampleAxis});
    histogramRegistry.add("subsample/pos_sq", "", kTProfile2D, {cent1Axis, subsampleAxis});
    histogramRegistry.add("subsample/neg_sq", "", kTProfile2D, {cent1Axis, subsampleAxis});
    histogramRegistry.add("subsample/posneg", "", kTProfile2D, {cent1Axis, subsampleAxis});

    if (cfgLoadEff) {
      ccdb->setURL(cfgUrlCCDB.value);
      ccdb->setCaching(true);
      ccdb->setLocalObjectValidityChecking();

      // ccdb->setCreatedNotAfter(ccdbNoLaterThan.value);
      // LOGF(info, "Getting object %s", ccdbPath.value.data());

      TList* list = ccdb->getForTimeStamp<TList>(cfgPathCCDB.value, -1);
      efficiency = reinterpret_cast<TH2D*>(list->FindObject("efficiency_Run3"));
      if (!efficiency) {
        LOGF(info, "FATAL!! Could not find required histograms in CCDB");
      }
    }
  }

  template <RunType run, typename C>
  bool selCollision(C const& coll, float& cent, float& mult)
  {

    if (std::abs(coll.posZ()) > vertexZcut)
      return false;

    if constexpr (run == kRun3) {
      if (cSel8Trig && !coll.sel8()) {
        return false;
      } // require min bias trigger
      cent = coll.centFT0M(); // centrality for run3
      mult = coll.multFT0M(); // multiplicity for run3
    } else if constexpr (run == kRun2) {
      if (cInt7Trig && !coll.alias_bit(kINT7)) {
        return false;
      }
      if (cSel7Trig && !coll.sel7()) {
        return false;
      }
      cent = coll.centRun2V0M(); // centrality for run2
      mult = coll.multFV0M();    // multiplicity for run2
    }

    if (cNoItsROBorder && !coll.selection_bit(aod::evsel::kNoITSROFrameBorder))
      return false;
    if (cTFBorder && !coll.selection_bit(aod::evsel::kNoTimeFrameBorder))
      return false;
    if (cPileupReject && !coll.selection_bit(aod::evsel::kNoSameBunchPileup))
      return false;
    if (cZVtxTimeDiff && !coll.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV))
      return false;
    if (cItsTpcVtx && !coll.selection_bit(aod::evsel::kIsVertexITSTPC))
      return false;
    if (cfgUseGoodItsLayerAllCut && !(coll.selection_bit(aod::evsel::kIsGoodITSLayersAll)))
      return false;

    return true;
  }

  template <typename T>
  void fillBeforeQA(T const& track)
  {
    histogramRegistry.fill(HIST("data/hTPCchi2perCluster_before"), track.tpcChi2NCl());
    histogramRegistry.fill(HIST("data/hITSchi2perCluster_before"), track.itsChi2NCl());
    histogramRegistry.fill(HIST("data/hTPCCrossedrows_before"), track.tpcNClsCrossedRows());
    histogramRegistry.fill(HIST("data/hDcaXY_before"), track.dcaXY());
    histogramRegistry.fill(HIST("data/hDcaZ_before"), track.dcaZ());
    histogramRegistry.fill(HIST("data/hPtDcaXY_before"), track.pt(), track.dcaXY());
    histogramRegistry.fill(HIST("data/hPtDcaZ_before"), track.pt(), track.dcaZ());
  }

  template <typename T>
  void fillAfterQA(T const& track)
  {
    histogramRegistry.fill(HIST("data/hDcaXY_after"), track.dcaXY());
    histogramRegistry.fill(HIST("data/hDcaZ_after"), track.dcaZ());
    histogramRegistry.fill(HIST("data/hPt"), track.pt());
    histogramRegistry.fill(HIST("data/hEta"), track.eta());
    histogramRegistry.fill(HIST("data/hPt_eta"), track.pt(), track.eta());
    histogramRegistry.fill(HIST("data/hPtDcaXY_after"), track.pt(), track.dcaXY());
    histogramRegistry.fill(HIST("data/hPtDcaZ_after"), track.pt(), track.dcaZ());
    histogramRegistry.fill(HIST("data/hTPCCrossedrows_after"), track.tpcNClsCrossedRows());
    histogramRegistry.fill(HIST("data/hTPCchi2perCluster_after"), track.tpcChi2NCl());
    histogramRegistry.fill(HIST("data/hITSchi2perCluster_after"), track.itsChi2NCl());
  }

  template <typename T>
  bool selTrack(T const& track)
  {
    if (!track.isGlobalTrack())
      return false;
    if (std::fabs(track.eta()) >= etaCut)
      return false;
    if (track.pt() <= ptMinCut || track.pt() >= ptMaxCut)
      return false;
    if (track.sign() == 0)
      return false;
    if (cDcaXy && std::fabs(track.dcaXY()) > dcaXYCut)
      return false;
    if (cDcaZ && std::fabs(track.dcaZ()) > dcaZCut)
      return false;
    if (cTpcCr && track.tpcNClsCrossedRows() < tpcCrossCut)
      return false;
    if (cItsChi && track.itsChi2NCl() > itsChiCut)
      return false;
    if (cTpcChi && track.tpcChi2NCl() > tpcChiCut)
      return false;

    return true;
  }

  double getEfficiency(double pt, double eta, TH2D* hEff)
  {
    if (!hEff) {
      LOGF(error, "Efficiency histogram is null — check CCDB loading.");
      return 1e-6;
    }
    int binX = hEff->GetXaxis()->FindBin(pt);
    int binY = hEff->GetYaxis()->FindBin(eta);
    if (binX < 1 || binX > hEff->GetNbinsX() || binY < 1 || binY > hEff->GetNbinsY()) {
      LOGF(warn, "pt or eta out of histogram bounds: pt = %f, eta = %f", pt, eta);
      return 1e-6;
    }
    double eff = hEff->GetBinContent(binX, binY);
    return eff;
  }

  void fillHistograms(float nch, float cent, float fpos, float fneg, float posneg, float termp, float termn)
  {
    histogramRegistry.fill(HIST("mult_dist/nch"), nch);
    histogramRegistry.fill(HIST("mult_dist/nch_cent"), cent, nch);
    histogramRegistry.fill(HIST("mult_dist/nch_pos"), fpos);
    histogramRegistry.fill(HIST("mult_dist/nch_pos_cent"), cent, fpos);
    histogramRegistry.fill(HIST("mult_dist/nch_neg"), fneg);
    histogramRegistry.fill(HIST("mult_dist/nch_neg_cent"), cent, fneg);
    histogramRegistry.fill(HIST("mult_dist/nch_negpos"), posneg);
    histogramRegistry.fill(HIST("mult_dist/nch_negpos_cent"), cent, posneg);

    histogramRegistry.fill(HIST("cent/pos"), cent, fpos);
    histogramRegistry.fill(HIST("cent/neg"), cent, fneg);
    histogramRegistry.fill(HIST("cent/termp"), cent, termp);
    histogramRegistry.fill(HIST("cent/termn"), cent, termn);
    histogramRegistry.fill(HIST("cent/pos_sq"), cent, fpos * fpos);
    histogramRegistry.fill(HIST("cent/neg_sq"), cent, fneg * fneg);
    histogramRegistry.fill(HIST("cent/posneg"), cent, posneg);

    float lRandom = fRndm->Rndm();
    int sampleIndex = static_cast<int>(cfgNSubsample * lRandom);

    histogramRegistry.fill(HIST("subsample/pos"), cent, sampleIndex, fpos);
    histogramRegistry.fill(HIST("subsample/neg"), cent, sampleIndex, fneg);
    histogramRegistry.fill(HIST("subsample/termp"), cent, sampleIndex, termp);
    histogramRegistry.fill(HIST("subsample/termn"), cent, sampleIndex, termn);
    histogramRegistry.fill(HIST("subsample/pos_sq"), cent, sampleIndex, fpos * fpos);
    histogramRegistry.fill(HIST("subsample/neg_sq"), cent, sampleIndex, fneg * fneg);
    histogramRegistry.fill(HIST("subsample/posneg"), cent, sampleIndex, posneg);
  }

  template <RunType run, typename C, typename T>
  void calculationData(C const& coll, T const& tracks)
  {
    float cent = -1, mult = -1;
    histogramRegistry.fill(HIST("data/hVtxZ_before"), coll.posZ());
    if (!selCollision<run>(coll, cent, mult)) {
      return;
    }
    histogramRegistry.fill(HIST("data/hVtxZ_after"), coll.posZ());
    histogramRegistry.fill(HIST("data/hCentrality"), cent);
    histogramRegistry.fill(HIST("data/hMultiplicity"), mult);

    int fpos = 0, fneg = 0, posneg = 0, termn = 0, termp = 0;
    int nch = 0, nchCor = 0;
    double posWeight = 0, negWeight = 0;
    for (const auto& track : tracks) {
      fillBeforeQA(track);
      if (!selTrack(track))
        continue;
      nch += 1;
      fillAfterQA(track);
      histogramRegistry.fill(HIST("data/hEta_cent"), cent, track.eta());
      histogramRegistry.fill(HIST("data/hPt_cent"), cent, track.pt());

      double eff = getEfficiency(track.pt(), track.eta(), efficiency);
      if (eff < threshold)
        continue;
      double weight = 1.0 / eff;

      histogramRegistry.fill(HIST("cor/hPt_cor"), track.pt(), weight);
      histogramRegistry.fill(HIST("cor/hEta_cor"), track.eta(), weight);

      nchCor += weight;
      if (track.sign() == 1) {
        fpos += 1;
        posWeight += weight;
      } else if (track.sign() == -1) {
        fneg += 1;
        negWeight += weight;
      }

    } // track
    termp = fpos * (fpos - 1);
    termn = fneg * (fneg - 1);
    posneg = fpos * fneg;
    histogramRegistry.fill(HIST("cor/nch_vs_nchCor"), nch, nchCor);
    histogramRegistry.fill(HIST("cor/nchCor"), nchCor);
    histogramRegistry.fill(HIST("cor/cent_nchCor"), cent, nchCor);
    histogramRegistry.fill(HIST("cor/fpos_cent"), cent, posWeight);
    histogramRegistry.fill(HIST("cor/fneg_cent"), cent, negWeight);
    fillHistograms(nch, cent, fpos, fneg, posneg, termp, termn);
  }

  template <RunType run, typename C, typename T, typename M, typename P>
  void calculationMc(C const& coll, T const& inputTracks, M const& mcCollisions, P const& mcParticles)
  {
    (void)mcCollisions;
    if (!coll.has_mcCollision()) {
      return;
    }
    histogramRegistry.fill(HIST("gen/hVtxZ_before"), coll.mcCollision().posZ());
    float cent = -1, mult = -1;
    histogramRegistry.fill(HIST("data/hVtxZ_before"), coll.posZ());
    if (!selCollision<run>(coll, cent, mult)) {
      return;
    }
    histogramRegistry.fill(HIST("data/hVtxZ_after"), coll.posZ());
    histogramRegistry.fill(HIST("data/hCentrality"), cent);
    histogramRegistry.fill(HIST("data/hMultiplicity"), mult);

    int fpos = 0, fneg = 0, posneg = 0, termn = 0, termp = 0;
    int nch = 0, nchCor = 0;
    double posRecWeight = 0, negRecWeight = 0;

    for (const auto& track : inputTracks) {
      fillBeforeQA(track);
      if (!selTrack(track))
        continue;
      nch += 1;
      fillAfterQA(track);
      histogramRegistry.fill(HIST("data/hEta_cent"), cent, track.eta());
      histogramRegistry.fill(HIST("data/hPt_cent"), cent, track.pt());

      double eff = getEfficiency(track.pt(), track.eta(), efficiency);
      if (eff < threshold)
        continue;
      double weight = 1.0 / eff;
      histogramRegistry.fill(HIST("cor/hPt_cor"), track.pt(), weight);
      histogramRegistry.fill(HIST("cor/hEta_cor"), track.eta(), weight);
      nchCor += weight;

      if (track.sign() == 1) {
        fpos += 1;
        posRecWeight += weight;
      } else if (track.sign() == -1) {
        fneg += 1;
        negRecWeight += weight;
      }
    } // track
    termp = fpos * (fpos - 1);

    termn = fneg * (fneg - 1);

    posneg = fpos * fneg;
    histogramRegistry.fill(HIST("cor/nch_vs_nchCor"), nch, nchCor);
    histogramRegistry.fill(HIST("cor/nchCor"), nchCor);
    histogramRegistry.fill(HIST("cor/cent_nchCor"), cent, nchCor);
    histogramRegistry.fill(HIST("cor/fpos_cent"), cent, posRecWeight);
    histogramRegistry.fill(HIST("cor/fneg_cent"), cent, negRecWeight);

    fillHistograms(nch, cent, fpos, fneg, posneg, termp, termn);

    int posGen = 0, negGen = 0, posNegGen = 0, termNGen = 0, termPGen = 0, nchGen = 0;

    const auto& mccolgen = coll.template mcCollision_as<aod::McCollisions>();
    if (std::abs(mccolgen.posZ()) > vertexZcut)
      return;
    const auto& mcpartgen = mcParticles.sliceByCached(aod::mcparticle::mcCollisionId, mccolgen.globalIndex(), cache);
    histogramRegistry.fill(HIST("gen/hVtxZ_after"), mccolgen.posZ());
    for (const auto& mcpart : mcpartgen) {
      if (!mcpart.isPhysicalPrimary())
        continue;
      int pid = mcpart.pdgCode();
      auto sign = 0;
      auto* pd = pdgService->GetParticle(pid);
      if (pd != nullptr) {
        sign = pd->Charge() / 3.;
      }
      if (sign == 0)
        continue;
      if (std::abs(pid) != kElectron && std::abs(pid) != kMuonMinus && std::abs(pid) != kPiPlus && std::abs(pid) != kKPlus && std::abs(pid) != kProton)
        continue;
      if (std::fabs(mcpart.eta()) > etaCut)
        continue;
      if ((mcpart.pt() <= ptMinCut) || (mcpart.pt() >= ptMaxCut))
        continue;
      histogramRegistry.fill(HIST("gen/hPt"), mcpart.pt());
      histogramRegistry.fill(HIST("gen/hPt_cent"), cent, mcpart.pt());
      histogramRegistry.fill(HIST("gen/hEta"), mcpart.eta());
      histogramRegistry.fill(HIST("gen/hEta_cent"), cent, mcpart.eta());
      histogramRegistry.fill(HIST("gen/hSign"), sign);
      histogramRegistry.fill(HIST("gen/hPt_eta"), mcpart.pt(), mcpart.eta());
      nchGen += 1;
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
    histogramRegistry.fill(HIST("cent/gen_pos"), cent, posGen);
    histogramRegistry.fill(HIST("cent/gen_neg"), cent, negGen);
    histogramRegistry.fill(HIST("cent/gen_termp"), cent, termPGen);
    histogramRegistry.fill(HIST("cent/gen_termn"), cent, termNGen);
    histogramRegistry.fill(HIST("cent/gen_pos_sq"), cent, posGen * posGen);
    histogramRegistry.fill(HIST("cent/gen_neg_sq"), cent, negGen * negGen);
    histogramRegistry.fill(HIST("cent/gen_posneg"), cent, posNegGen);
    histogramRegistry.fill(HIST("cent/gen_nch"), cent, nchGen);
    histogramRegistry.fill(HIST("gen/nch"), nchGen);

  } // void

  template <RunType run, typename C, typename T>
  void calculationDeltaEta(C const& coll, T const& tracks, float deta1, float deta2)
  {
    float cent = -1, mult = -1;
    if (!selCollision<run>(coll, cent, mult))
      return;
    if (!(cent >= centMin && cent < centMax))
      return;
    histogramRegistry.fill(HIST("delta_eta/cent"), cent);

    int fpos = 0, fneg = 0, posneg = 0, termn = 0, termp = 0;
    for (const auto& track : tracks) {
      if (!selTrack(track))
        continue;
      double eta = track.eta();
      if (eta < deta1 || eta > deta2)
        continue;

      histogramRegistry.fill(HIST("delta_eta/track_eta"), eta);

      if (track.sign() == 1)
        fpos++;
      else if (track.sign() == -1)
        fneg++;
    }
    termp = fpos * (fpos - 1);
    termn = fneg * (fneg - 1);
    posneg = fpos * fneg;

    float deltaEtaWidth = deta2 - deta1 + 1e-5f;

    histogramRegistry.fill(HIST("delta_eta/pos"), deltaEtaWidth, fpos);
    histogramRegistry.fill(HIST("delta_eta/neg"), deltaEtaWidth, fneg);
    histogramRegistry.fill(HIST("delta_eta/termp"), deltaEtaWidth, termp);
    histogramRegistry.fill(HIST("delta_eta/termn"), deltaEtaWidth, termn);
    histogramRegistry.fill(HIST("delta_eta/pos_sq"), deltaEtaWidth, fpos * fpos);
    histogramRegistry.fill(HIST("delta_eta/neg_sq"), deltaEtaWidth, fneg * fneg);
    histogramRegistry.fill(HIST("delta_eta/posneg"), deltaEtaWidth, posneg);
  }

  SliceCache cache;
  Preslice<aod::McParticles> mcTrack = o2::aod::mcparticle::mcCollisionId;

  // process function for Data Run3
  void processDataRun3(aod::MyCollisionRun3 const& coll, aod::MyTracks const& tracks)
  {
    calculationData<kRun3>(coll, tracks);
    for (int ii = 0; ii < deltaEta; ii++) {
      float etaMin = -0.1f * (ii + 1);
      float etaMax = 0.1f * (ii + 1);

      calculationDeltaEta<kRun3>(coll, tracks, etaMin, etaMax);
    }
  }

  PROCESS_SWITCH(NetchargeFluctuations, processDataRun3, "Process for Run3 DATA", false);

  // process function for Data Run2
  void processDataRun2(aod::MyCollisionRun2 const& coll, aod::MyTracks const& tracks)
  {
    calculationData<kRun2>(coll, tracks);
    for (int ii = 0; ii < deltaEta; ii++) {
      float etaMin = -0.1f * (ii + 1); // -0.1, -0.2, ..., -0.8
      float etaMax = 0.1f * (ii + 1);  // +0.1, +0.2, ..., +0.8

      calculationDeltaEta<kRun2>(coll, tracks, etaMin, etaMax);
    }
  }

  PROCESS_SWITCH(NetchargeFluctuations, processDataRun2, "Process for Run2 DATA", false);

  // process function for MC Run3

  void processMcRun3(aod::MyMCCollisionRun3 const& coll, aod::MyMCTracks const& inputTracks,
                     aod::McCollisions const& mcCollisions, aod::McParticles const& mcParticles)
  {
    calculationMc<kRun3>(coll, inputTracks, mcCollisions, mcParticles);

    for (int ii = 0; ii < deltaEta; ii++) {
      float etaMin = -0.1f * (ii + 1);
      float etaMax = 0.1f * (ii + 1);

      calculationDeltaEta<kRun3>(coll, inputTracks, etaMin, etaMax);
    }
  }

  PROCESS_SWITCH(NetchargeFluctuations, processMcRun3, "Process reconstructed", false);

  // process function for MC Run2

  void processMcRun2(aod::MyMCCollisionRun2 const& coll, aod::MyMCTracks const& inputTracks,
                     aod::McCollisions const& mcCollisions, aod::McParticles const& mcParticles)
  {
    calculationMc<kRun2>(coll, inputTracks, mcCollisions, mcParticles);
    for (int ii = 0; ii < deltaEta; ii++) {
      float etaMin = -0.1f * (ii + 1); // -0.1, -0.2, ..., -0.8
      float etaMax = 0.1f * (ii + 1);  // +0.1, +0.2, ..., +0.8

      calculationDeltaEta<kRun2>(coll, inputTracks, etaMin, etaMax);
    }
  }

  PROCESS_SWITCH(NetchargeFluctuations, processMcRun2, "Process reconstructed", true);
};

// struct
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    {adaptAnalysisTask<NetchargeFluctuations>(cfgc)}};
}
