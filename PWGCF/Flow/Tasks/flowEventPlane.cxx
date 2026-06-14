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

/// \file flowEventPlane.cxx
/// \brief Flow calculation using event plane.
/// \author Yash Patley <yash.patley@cern.ch>

#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;
using namespace o2::constants::math;

namespace o2::aod
{
namespace colsp
{
DECLARE_SOA_COLUMN(RunNumber, runNumber, int);
DECLARE_SOA_COLUMN(Timestamp, timestamp, uint64_t);
DECLARE_SOA_COLUMN(Cent, cent, float);
DECLARE_SOA_COLUMN(Vx, vx, float);
DECLARE_SOA_COLUMN(Vy, vy, float);
DECLARE_SOA_COLUMN(Vz, vz, float);
DECLARE_SOA_COLUMN(ZnaEnergyCommon, znaEnergyCommon, float);
DECLARE_SOA_COLUMN(ZncEnergyCommon, zncEnergyCommon, float);
DECLARE_SOA_COLUMN(ZnaEnergy, znaEnergy, float[4]);
DECLARE_SOA_COLUMN(ZncEnergy, zncEnergy, float[4]);
} // namespace colsp
DECLARE_SOA_TABLE(ColSps, "AOD", "COLSPS", o2::soa::Index<>,
                  colsp::RunNumber,
                  colsp::Timestamp,
                  colsp::Cent,
                  colsp::Vx,
                  colsp::Vy,
                  colsp::Vz,
                  colsp::ZnaEnergyCommon,
                  colsp::ZncEnergyCommon,
                  colsp::ZnaEnergy,
                  colsp::ZncEnergy);

using ColSp = ColSps::iterator;

namespace tracksid
{
DECLARE_SOA_INDEX_COLUMN(ColSp, colSp);
DECLARE_SOA_COLUMN(Px, px, float);
DECLARE_SOA_COLUMN(Py, py, float);
DECLARE_SOA_COLUMN(Pz, pz, float);
DECLARE_SOA_COLUMN(Charge, charge, int8_t);
DECLARE_SOA_COLUMN(TpcInfo, tpcInfo, float[4]);
DECLARE_SOA_COLUMN(TofInfo, tofInfo, float[4]);
DECLARE_SOA_COLUMN(Dca, dca, float[2]);
} // namespace tracksid
DECLARE_SOA_TABLE(TracksId, "AOD", "TRACKSID", o2::soa::Index<>,
                  tracksid::ColSpId,
                  tracksid::Px,
                  tracksid::Py,
                  tracksid::Pz,
                  tracksid::Charge,
                  tracksid::TpcInfo,
                  tracksid::TofInfo,
                  tracksid::Dca);
using TrackId = TracksId::iterator;

namespace v0track
{
DECLARE_SOA_INDEX_COLUMN(ColSp, colSp);
DECLARE_SOA_COLUMN(Px, px, float);
DECLARE_SOA_COLUMN(Py, py, float);
DECLARE_SOA_COLUMN(Pz, pz, float);
DECLARE_SOA_COLUMN(CosPA, cosPA, float);
DECLARE_SOA_COLUMN(Ctau, ctau, float);
DECLARE_SOA_COLUMN(Mass, mass, float);
DECLARE_SOA_COLUMN(Species, species, uint8_t);
DECLARE_SOA_COLUMN(DcaDau, dcaDau, float[2]);
DECLARE_SOA_COLUMN(TpcDau, tpcDau, float[2]);
} // namespace v0track
DECLARE_SOA_TABLE(V0Tracks, "AOD", "V0TRACKS", o2::soa::Index<>,
                  v0track::ColSpId,
                  v0track::Px,
                  v0track::Py,
                  v0track::Pz,
                  v0track::CosPA,
                  v0track::Ctau,
                  v0track::Mass,
                  v0track::Species,
                  v0track::DcaDau,
                  v0track::TpcDau);
using V0Track = V0Tracks::iterator;
} // namespace o2::aod

enum TrackType {
  kCharged = 0,
  kPi,
  kKa,
  kPr
};

enum V0Type {
  kK0S = 0,
  kLambda,
  kAntiLambda
};

struct FlowEventPlane {
  // Table producer
  Produces<aod::ColSp> colSpTable;
  Produces<aod::TracksId> tracksIdTable;
  Produces<aod::V0Tracks> v0sTable;

  // Configurables
  // Collisions
  Configurable<float> cMinZVtx{"cMinZVtx", -10.0, "Min VtxZ cut"};
  Configurable<float> cMaxZVtx{"cMaxZVtx", 10.0, "Max VtxZ cut"};
  Configurable<float> cMinCent{"cMinCent", 0., "Minumum Centrality"};
  Configurable<float> cMaxCent{"cMaxCent", 100.0, "Maximum Centrality"};
  Configurable<bool> cSel8Trig{"cSel8Trig", true, "Sel8 (T0A + T0C) Selection Run3"};
  Configurable<bool> cPileupReject{"cPileupReject", true, "Pileup rejection"};
  Configurable<bool> cZVtxTimeDiff{"cZVtxTimeDiff", true, "z-vtx time diff selection"};
  Configurable<bool> cIsGoodITSLayers{"cIsGoodITSLayers", true, "Good ITS Layers All"};
  Configurable<float> cMinOccupancy{"cMinOccupancy", 0, "Minimum FT0C Occupancy"};
  Configurable<float> cMaxOccupancy{"cMaxOccupancy", 1e6, "Maximum FT0C Occupancy"};

  // Tracks
  Configurable<float> cTrackMinPt{"cTrackMinPt", 0.1, "p_{T} minimum"};
  Configurable<float> cTrackMaxPt{"cTrackMaxPt", 10.0, "p_{T} maximum"};
  Configurable<float> cTrackEtaCut{"cTrackEtaCut", 0.8, "Pseudorapidity cut"};
  Configurable<bool> cTrackGlobal{"cTrackGlobal", true, "Global Track"};
  Configurable<float> cTrackDcaXYCut{"cTrackDcaXYCut", 0.1, "DcaXY Cut"};
  Configurable<float> cTrackDcaZCut{"cTrackDcaZCut", 1., "DcaXY Cut"};

  // Track PID
  Configurable<float> cTpcElRejCutMin{"cTpcElRejCutMin", -3., "Electron Rejection Cut Minimum"};
  Configurable<float> cTpcElRejCutMax{"cTpcElRejCutMax", 5., "Electron Rejection Cut Maximum"};
  Configurable<float> cTpcNSigmaCut{"cTpcNSigmaCut", 5, "TPC NSigma Cut"};
  Configurable<float> cTpcRejCut{"cTpcRejCut", 3, "TPC Rej Cut"};
  Configurable<float> cTofNSigmaCut{"cTofNSigmaCut", 5, "TOF NSigma Cut"};
  Configurable<float> cTofRejCut{"cTofRejCut", 3, "TOF Rej Cut"};

  // V0s
  Configurable<int> cV0TypeSelection{"cV0TypeSelection", 1, "V0 Type Selection"};
  Configurable<float> cMinDcaProtonToPV{"cMinDcaProtonToPV", 0.02, "Minimum Proton DCAr to PV"};
  Configurable<float> cMinDcaPionToPV{"cMinDcaPionToPV", 0.06, "Minimum Pion DCAr to PV"};
  Configurable<float> cDcaV0Dau{"cDcaV0Dau", 1., "DCA between V0 daughters"};
  Configurable<float> cMinV0Radius{"cMinV0Radius", 0.5, "Minimum V0 radius from PV"};
  Configurable<float> cV0CTau{"cV0CTau", 50.0, "Decay length selection"};
  Configurable<float> cV0CosPA{"cV0CosPA", 0.95, "CosPA selection"};
  Configurable<float> cArmPodSel{"cArmPodSel", 0.2, "Armentros-Podolanski Selection for K0S"};
  Configurable<float> cV0EtaCut{"cV0EtaCut", 0.8, "V0 eta cut"};

  // Histogram Registry.
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // Run number
  uint64_t runNum = 0, timestamp = 0;
  int64_t colIdx = 0;
  float mult = 0., cent = 0., posX = 0., posY = 0., posZ = 0.;

  void init(InitContext const&)
  {
    // Axis collision
    const AxisSpec axisCent{100, 0., 100, "FT0C%"};

    const AxisSpec axisTrackPt(40, 0, 4, "p_{T} (GeV/#it{c})");
    const AxisSpec axisTrackDcaZ(220, -1.1, 1.1, "dca_{Z} (cm)");
    const AxisSpec axisTrackDcaXY(80, -0.2, 0.2, "dca_{XY} (cm)");
    const AxisSpec axisTrackdEdx(360, 20, 200, "#frac{dE}{dx}");
    const AxisSpec axisTrackTofSignal(240, 0, 1.2, "#beta");

    const AxisSpec axisAlpha(40, -1, 1, "#alpha");
    const AxisSpec axisQtarm(40, 0, 0.4, "q_{T}");

    // Histogram collision
    histos.add("hCent", "Centrality", kTH1F, {axisCent});

    // Tracks QA
    histos.add("QA/hDcaXY", "DCA_{XY} vs pT", kTH2F, {axisTrackPt, axisTrackDcaXY});
    histos.add("QA/hDcaZ", "DCA_{Z} vs p_{T}", kTH2F, {axisTrackPt, axisTrackDcaZ});
    histos.add("QA/hdEdX", "dE/dx vs pT", kTH2F, {axisTrackPt, axisTrackdEdx});
    histos.add("QA/hTOFSignal", "#beta_{TOF} vs p_{T}", kTH2F, {axisTrackPt, axisTrackTofSignal});

    // V0s QA
    histos.add("QA/h2f_armpod_before_sel", "Armentros-Podolanski Plot", kTH2F, {axisAlpha, axisQtarm});
    histos.add("QA/h2f_armpod_K0s", "Armentros-Podolanski Plot", kTH2F, {axisAlpha, axisQtarm});
    histos.add("QA/h2f_armpod_Lambda", "Armentros-Podolanski Plot", kTH2F, {axisAlpha, axisQtarm});
    histos.add("QA/h2f_armpod_AntiLambda", "Armentros-Podolanski Plot", kTH2F, {axisAlpha, axisQtarm});
  }

  template <typename C>
  bool selCollision(C const& col)
  {
    if (col.posZ() <= cMinZVtx || col.posZ() >= cMaxZVtx) { // VtxZ selection
      return false;
    }

    if (cSel8Trig && !col.sel8()) { // Sel8 selection
      return false;
    }

    cent = col.centFT0C();
    if (cent <= cMinCent || cent >= cMaxCent) { // Centrality selection
      return false;
    }

    if (col.ft0cOccupancyInTimeRange() < cMinOccupancy || col.ft0cOccupancyInTimeRange() > cMaxOccupancy) { // Occupancy cut
      return false;
    }

    if (cPileupReject && !col.selection_bit(aod::evsel::kNoSameBunchPileup)) { // Pile-up rejection
      return false;
    }

    if (cZVtxTimeDiff && !col.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV)) { // ZvtxFT0 vs PV
      return false;
    }

    if (cIsGoodITSLayers && !col.selection_bit(aod::evsel::kIsGoodITSLayersAll)) { // All ITS layer active
      return false;
    }

    // Set Multiplicity
    mult = col.multTPC();

    return true;
  }

  // Track Selection
  template <typename T>
  bool selectTrack(T const& track)
  {
    if (track.pt() <= cTrackMinPt || track.pt() >= cTrackMaxPt || std::abs(track.eta()) >= cTrackEtaCut) {
      return false;
    }

    if (cTrackGlobal && !track.isGlobalTrackWoDCA()) {
      return false;
    }

    if (std::abs(track.dcaXY()) >= cTrackDcaXYCut || std::abs(track.dcaZ()) >= cTrackDcaZCut) {
      return false;
    }

    return true;
  }

  // PID
  template <typename T>
  void identifyTrack(T const& track, std::array<float, 4>& nstpc, std::array<float, 4>& nstof)
  {
    // Electron rejection
    if (std::abs(track.tpcNSigmaPi()) > cTpcRejCut && std::abs(track.tpcNSigmaKa()) > cTpcRejCut && std::abs(track.tpcNSigmaPr()) > cTpcRejCut && track.tpcNSigmaEl() > cTpcElRejCutMin && track.tpcNSigmaEl() < cTpcElRejCutMax) {
      return;
    }

    // Check track for pion
    if (std::abs(track.tpcNSigmaPi()) < cTpcNSigmaCut) {
      nstpc[kPi] = track.tpcNSigmaPi();
      if (track.hasTOF() && std::abs(track.tofNSigmaPi()) < cTofNSigmaCut) {
        nstof[kPi] = track.tofNSigmaPi();
      }
    }

    // Check track for kaon
    if (std::abs(track.tpcNSigmaKa()) < cTpcNSigmaCut) {
      nstpc[kKa] = track.tpcNSigmaKa();
      if (track.hasTOF() && std::abs(track.tofNSigmaKa()) < cTofNSigmaCut) {
        nstof[kKa] = track.tofNSigmaKa();
      }
    }

    // Check track for proton
    if (std::abs(track.tpcNSigmaPr()) < cTpcNSigmaCut) {
      nstpc[kPr] = track.tpcNSigmaPr();
      if (track.hasTOF() && std::abs(track.tofNSigmaPr()) < cTofNSigmaCut) {
        nstof[kPr] = track.tofNSigmaPr();
      }
    }
  }

  // Id track table
  template <typename T>
  void analyzeIdHadrons(T const& tracks)
  {
    for (auto const& track : tracks) {
      // Global track selection
      if (!selectTrack(track)) {
        continue;
      }

      // Dca
      std::array<float, 2> dca = {track.dcaXY(), track.dcaZ()};

      // QA plots
      histos.fill(HIST("QA/hDcaXY"), track.pt(), track.dcaXY());
      histos.fill(HIST("QA/hDcaZ"), track.pt(), track.dcaZ());
      histos.fill(HIST("QA/hdEdX"), track.pt(), track.tpcSignal());
      if (track.hasTOF()) {
        histos.fill(HIST("QA/hTOFSignal"), track.pt(), track.beta());
      }

      // Identify track
      std::array<float, 4> tpc = {track.tpcSignal(), -99., -99., -99.};
      std::array<float, 4> tof = {track.beta(), -99., -99., -99.};
      identifyTrack(track, tpc, tof);

      // Fill track table
      tracksIdTable(colIdx, track.px(), track.py(), track.pz(), track.sign(), tpc.data(), tof.data(), dca.data());
    }
  }

  // V0 daughter track selection
  template <V0Type part, typename V, typename T>
  bool selV0DauTracks(V const& v0, T const& postrack, T const& negtrack, float& tpcNSigmaPosDau, float& tpcNSigmaNegDau)
  {
    // Kinematic selection
    if (postrack.pt() <= cTrackMinPt || negtrack.pt() <= cTrackMinPt || std::abs(postrack.eta()) >= cTrackEtaCut || std::abs(negtrack.eta()) >= cTrackEtaCut) {
      return false;
    }

    // Apply DCA Selection on Daughter Tracks Based on Lambda/AntiLambda/K0Short daughters
    if (part == kK0S) {
      if (std::abs(v0.dcapostopv()) <= cMinDcaPionToPV || std::abs(v0.dcanegtopv()) <= cMinDcaPionToPV) {
        return false;
      }
    } else if (part == kLambda) {
      if (std::abs(v0.dcapostopv()) <= cMinDcaProtonToPV || std::abs(v0.dcanegtopv()) <= cMinDcaPionToPV) {
        return false;
      }
    } else if (part == kAntiLambda) {
      if (std::abs(v0.dcapostopv()) <= cMinDcaPionToPV || std::abs(v0.dcanegtopv()) <= cMinDcaProtonToPV) {
        return false;
      }
    } else {
      return false;
    }

    // Daughter track PID
    switch (part) {
      // PosDau = Pion, NegDau = Pion
      case kK0S:
        tpcNSigmaPosDau = postrack.tpcNSigmaPi();
        tpcNSigmaNegDau = negtrack.tpcNSigmaPi();
        break;

      // PosDau = Proton, NegDau = Pion
      case kLambda:
        tpcNSigmaPosDau = postrack.tpcNSigmaPr();
        tpcNSigmaNegDau = negtrack.tpcNSigmaPi();
        break;

      // PosDau = Pion, NegDau = Proton
      case kAntiLambda:
        tpcNSigmaPosDau = postrack.tpcNSigmaPi();
        tpcNSigmaNegDau = negtrack.tpcNSigmaPr();
        break;
    }

    if (std::abs(tpcNSigmaPosDau) >= cTpcNSigmaCut || std::abs(tpcNSigmaNegDau) >= cTpcNSigmaCut) {
      return false;
    }

    // All checks passed
    return true;
  }

  // V0 table
  template <V0Type v0type, typename V, typename T>
  void fillV0Table(V const& v0, T const&)
  {
    // V0 parameters
    float mass = 0., ctau = 0., tpcNSigmaPosDau = 0., tpcNSigmaNegDau = 0.;
    auto postrack = v0.template posTrack_as<T>();
    auto negtrack = v0.template negTrack_as<T>();

    if (v0type == kK0S) {
      if (!selV0DauTracks<kK0S>(v0, postrack, negtrack, tpcNSigmaPosDau, tpcNSigmaNegDau)) {
        return;
      }
      mass = v0.mK0Short();
      ctau = v0.distovertotmom(posX, posY, posZ) * MassKaonNeutral;
    } else if (v0type == kLambda) {
      if (!selV0DauTracks<kLambda>(v0, postrack, negtrack, tpcNSigmaPosDau, tpcNSigmaNegDau)) {
        return;
      }
      mass = v0.mLambda();
      ctau = v0.distovertotmom(posX, posY, posZ) * MassLambda0;
    } else if (v0type == kAntiLambda) {
      if (!selV0DauTracks<kAntiLambda>(v0, postrack, negtrack, tpcNSigmaPosDau, tpcNSigmaNegDau)) {
        return;
      }
      mass = v0.mAntiLambda();
      ctau = v0.distovertotmom(posX, posY, posZ) * MassLambda0;
    } else {
      return;
    }

    // Apply V0 selection
    if (ctau > cV0CTau || v0.v0cosPA() < cV0CosPA) {
      return;
    }

    // Dca daughters
    std::array<float, 2> dcaDau = {v0.dcapostopv(), v0.dcanegtopv()};

    // Tpc daughters
    std::array<float, 2> tpcDau = {tpcNSigmaPosDau, tpcNSigmaNegDau};

    // Fill QA Plots
    if (v0type == kK0S) {
      histos.fill(HIST("QA/h2f_armpod_K0s"), v0.alpha(), v0.qtarm());
    } else if (v0type == kLambda) {
      histos.fill(HIST("QA/h2f_armpod_Lambda"), v0.alpha(), v0.qtarm());
    } else if (v0type == kAntiLambda) {
      histos.fill(HIST("QA/h2f_armpod_AntiLambda"), v0.alpha(), v0.qtarm());
    }

    // Fill V0 table
    v0sTable(colIdx, v0.px(), v0.py(), v0.pz(), v0.v0cosPA(), ctau, mass, (uint8_t)v0type, dcaDau.data(), tpcDau.data());
  }

  // V0s
  template <typename V, typename T>
  void analyzeV0s(V const& V0s, T const& tracks)
  {
    for (auto const& v0 : V0s) {
      // Topological and kinematic selections
      if (std::abs(v0.eta()) >= cV0EtaCut || v0.dcaV0daughters() >= cDcaV0Dau || v0.v0radius() <= cMinV0Radius || v0.v0Type() != cV0TypeSelection) {
        continue;
      }

      // Armenteros-Podolanski Selections
      histos.fill(HIST("QA/h2f_armpod_before_sel"), v0.alpha(), v0.qtarm());

      // K0s, Lambda, Anti-Lambda
      if (v0.qtarm() >= cArmPodSel * std::abs(v0.alpha())) { // K0s
        fillV0Table<kK0S>(v0, tracks);
      } else if ((v0.qtarm() < cArmPodSel * std::abs(v0.alpha())) && (v0.alpha() > 0)) { // Lambda
        fillV0Table<kLambda>(v0, tracks);
      } else if ((v0.qtarm() < cArmPodSel * std::abs(v0.alpha())) && (v0.alpha() < 0)) { // Anti-Lambda
        fillV0Table<kAntiLambda>(v0, tracks);
      } else {
        return;
      }
    }
  }

  template <typename B, typename C, typename T, typename V>
  void analyzeCollision(B const& bc, C const& collision, T const& tracks, V const& V0s)
  {
    // Event selection
    if (!selCollision(collision)) {
      return;
    }
    posX = collision.posX();
    posY = collision.posY();
    posZ = collision.posZ();

    // check zdc
    if (!bc.has_zdc()) {
      return;
    }

    auto zdc = bc.zdc();
    auto znaEnergy = zdc.energySectorZNA();
    auto zncEnergy = zdc.energySectorZNC();
    auto znaEnergyCommon = zdc.energyCommonZNA();
    auto zncEnergyCommon = zdc.energyCommonZNC();

    // check energy deposits
    if (znaEnergyCommon <= 0 || zncEnergyCommon <= 0 || znaEnergy[0] <= 0 || znaEnergy[1] <= 0 || znaEnergy[2] <= 0 || znaEnergy[3] <= 0 || zncEnergy[0] <= 0 || zncEnergy[1] <= 0 || zncEnergy[2] <= 0 || zncEnergy[3] <= 0) {
      return;
    }

    // Fill collision table
    histos.fill(HIST("hCent"), cent);
    colSpTable(runNum, timestamp, cent, posX, posY, posZ, znaEnergyCommon, zncEnergyCommon, znaEnergy.data(), zncEnergy.data());
    colIdx = colSpTable.lastIndex();

    // Analyze Tracks [charged, pi, K, p]
    analyzeIdHadrons(tracks);

    // Analyze V0s [Lambda, K0S]
    analyzeV0s(V0s, tracks);

    // Done
    return;
  }

  using BCsRun3 = soa::Join<aod::BCsWithTimestamps, aod::Run3MatchedToBCSparse>;
  using CollisionsRun3 = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs, aod::CentFT0Ms, aod::CentFV0As, aod::TPCMults, aod::PVMults, aod::MultsGlobal, aod::MultsExtra>;
  using Tracks = soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra, aod::TracksDCA, aod::TOFSignal, aod::pidTOFbeta, aod::pidTPCEl, aod::pidTPCPi, aod::pidTOFPi, aod::pidTPCKa, aod::pidTOFKa, aod::pidTPCPr, aod::pidTOFPr, aod::TrackCompColls>;

  void processDummy(CollisionsRun3::iterator const&) {}

  PROCESS_SWITCH(FlowEventPlane, processDummy, "Dummy process", true);

  void processFEP(CollisionsRun3::iterator const& collision, BCsRun3 const&, aod::Zdcs const&,
                  Tracks const& tracks, aod::V0Datas const& v0s)
  {
    // Get bunch crossing
    auto bc = collision.template foundBC_as<BCsRun3>();
    runNum = collision.template foundBC_as<BCsRun3>().runNumber();
    timestamp = collision.template foundBC_as<BCsRun3>().timestamp();
    analyzeCollision(bc, collision, tracks, v0s);
  }

  PROCESS_SWITCH(FlowEventPlane, processFEP, "Flow Event Plane Table Producer", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<FlowEventPlane>(cfgc)};
}
