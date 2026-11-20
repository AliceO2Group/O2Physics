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

/// \file derivedDataCreatorElectronWeakBoson.cxx
/// \brief Producer of the derived data read by taskElectronWeakBoson (processDerived)
///
/// What is decided here and can NOT be changed afterwards on the derived data:
///  - Zorro trigger selection (cfgTriggerName) and |zvtx| < vtxZ
///  - event skim: at least one track with pT > ptSkimMin matched to an EMCal cluster within rMatchSkim,
///    and (skimRequireIsoElectron) passing the loose isolated-electron cuts eopSkim*, m02Skim*, nsigTpcSkim*, isoEnergySkimMax
///  - which tracks are stored: (global track w/o pT-eta, |eta| < etaTrMax, pT > ptAssMin)  OR  pT > ptEleStoreMin
///  - isolation cone radius rIsolation (track and cluster isolation sums)
///  - EMCal cluster definition and the time window used to pick the closest cluster
///  - random-cone UE density (rcHardE, nRandomCones)
/// Everything else (electron ID, E/p, M02, isolation thresholds, matching cut, Z mass window, KF chi2,
/// track-quality cuts of the associated hadrons, UE subtraction on/off) stays free in the analysis.
///
/// The event normalisation (hEventCounter, zorroSummary) lives in the AnalysisResults.root of THIS task.
///
/// \author S. Sakai & S. Ito (Univ. of Tsukuba)

#include "PWGHF/HFL/DataModel/ElectronWeakBosonTables.h"
#include "PWGJE/DataModel/EMCALClusters.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/Zorro.h"
#include "Common/Core/ZorroSummary.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/MathConstants.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <DetectorsBase/Propagator.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/runDataProcessing.h>

#include <TRandom3.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct HfDerivedDataCreatorElectronWeakBoson {
  Produces<o2::aod::HfWbColls> rowColl;
  Produces<o2::aod::HfWbTracks> rowTrack;
  Produces<o2::aod::HfWbElectrons> rowEle;

  // event
  Configurable<float> vtxZ{"vtxZ", 10.f, "|zvtx| cut applied at production"};
  Configurable<bool> applySel8{"applySel8", false, "apply sel8 at production (it is stored anyway)"};
  Configurable<bool> cfgSkimmedProcessing{"cfgSkimmedProcessing", true, "select with Zorro (skimmed data)"};
  Configurable<std::string> cfgTriggerName{"cfgTriggerName", "fGammaHighPtEMCAL", "Trigger of interest (comma separated for multiple)"};
  Configurable<std::string> cfgCCDBPath{"cfgCCDBPath", "EventFiltering/Zorro/", "Path to CCDB for trigger data"};
  Configurable<std::string> ccdbPathGrpMag{"ccdbPathGrpMag", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object (Run 3)"};

  // event skim. Grouped: in the JSON these options sit in their own block,
  //   "derived-data-creator-electron-weak-boson": { ..., "skim": { "ptSkimMin": "10", ... } }
  // The isolated-electron cuts must be LOOSER than (or equal to) the analysis ones, incl. every systematic variation.
  struct : ConfigurableGroup {
    std::string prefix = "skim"; // JSON group name
    Configurable<float> ptSkimMin{"ptSkimMin", 10.f, "keep the event only if a track above this pT is matched to EMCal"};
    Configurable<float> rMatchSkim{"rMatchSkim", 0.05f, "|dEta|,|dPhi| matching window used by the skim"};
    Configurable<bool> skimRequireIsoElectron{"skimRequireIsoElectron", true, "require at least one isolated-electron candidate (cuts below)"};
    Configurable<float> eopSkimMin{"eopSkimMin", 0.6f, "skim: min E/p"};
    Configurable<float> eopSkimMax{"eopSkimMax", 999.f, "skim: max E/p"};
    Configurable<float> m02SkimMin{"m02SkimMin", 0.1f, "skim: min M02"};
    Configurable<float> m02SkimMax{"m02SkimMax", 0.9f, "skim: max M02"};
    Configurable<float> isoEnergySkimMax{"isoEnergySkimMax", 0.2f, "skim: max cluster isolation (E_cone - E_UE)/E"};
    Configurable<bool> useUEsubSkim{"useUEsubSkim", true, "skim: subtract the random-cone UE in the isolation (true is the looser choice)"};
    Configurable<float> nsigTpcSkimMin{"nsigTpcSkimMin", -999.f, "skim: min TPC n sigma electron (-999 = no cut)"};
    Configurable<float> nsigTpcSkimMax{"nsigTpcSkimMax", 999.f, "skim: max TPC n sigma electron (999 = no cut)"};
  } skimCfg;

  // what is stored
  Configurable<float> etaTrMax{"etaTrMax", 1.0f, "|eta| of the stored associated tracks"};
  Configurable<float> ptAssMin{"ptAssMin", 0.15f, "min pT of the stored associated tracks"};
  Configurable<float> ptEleStoreMin{"ptEleStoreMin", 3.0f, "every track above this pT gets an HfWbElectrons row (must be <= ptMin and ptZeeMin of the analysis)"};

  // EMCal / isolation (fixed at production)
  Configurable<int> clusterDefinition{"clusterDefinition", 10, "cluster definition, e.g. 10=kV3Default"};
  Configurable<float> timeEmcMin{"timeEmcMin", -25., "min EMC cluster time for the closest-cluster search"};
  Configurable<float> timeEmcMax{"timeEmcMax", +20., "max EMC cluster time for the closest-cluster search"};
  Configurable<float> rMatchQa{"rMatchQa", 0.05f, "window used for the stored emcNMatchInR (QA only)"};
  Configurable<float> rIsolation{"rIsolation", 0.3, "cone radius for track and cluster isolation"};
  Configurable<float> etaEmcMax{"etaEmcMax", 0.6f, "random cones: EMCal |eta|"};
  Configurable<float> phiEmcMin{"phiEmcMin", 1.39, "random cones: EMCal phi min"};
  Configurable<float> phiEmcMax{"phiEmcMax", 3.36, "random cones: EMCal phi max"};
  Configurable<int> nRandomCones{"nRandomCones", 100, "number of random cones"};
  Configurable<float> rcHardE{"rcHardE", 5.0, "hard cluster veto energy in the random cones"};

  Service<o2::ccdb::BasicCCDBManager> ccdb{};
  Zorro zorro;
  TRandom3* rnd = nullptr;
  int lastRunNumber = -1;
  float magneticField = 0.f;

  using CollisionsWithCent = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0As, aod::CentFT0Cs, aod::CentFT0Ms, aod::CentFV0As, aod::Mults>;
  using TrackEle = soa::Join<aod::Tracks, aod::TracksCov, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullEl, aod::pidTPCFullPi>;
  using Clusters = o2::aod::EMCALClusters;

  // collision selection is cheap and identical to processData (posZ). sel8 is applied by hand.
  Filter posZFilter = (nabs(o2::aod::collision::posZ) < vtxZ);

  PresliceUnsorted<o2::aod::EMCALMatchedTracks> perTrackMatched = o2::aod::emcalmatchedtrack::trackId;

  HistogramRegistry registry{"registry"};
  OutputObj<ZorroSummary> zorroSummary{"zorroSummary"};

  struct HfEmcMatch {
    bool has = false;
    float dEta = 999.f, dPhi = 999.f, energy = -1.f, eta = 0.f, phi = 0.f, m02 = -1.f, time = 0.f, coneEnergy = 0.f;
    int16_t nCells = 0;
    uint8_t nMatchInR = 0;
  };

  void init(InitContext const&)
  {
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    zorro.setBaseCCDBPath(cfgCCDBPath.value);
    if (cfgSkimmedProcessing) {
      zorroSummary.setObject(zorro.getZorroSummary());
    }
    rnd = new TRandom3(0);

    const AxisSpec axisCounter{1, 0, 1, "events"};
    const AxisSpec axisTrigger{3, -0.5, 2.5, "Trigger status of zorro"};
    const AxisSpec axisZvtx{40, -20, 20, "Zvtx"};
    const AxisSpec axisCent{101, 0, 101, "FT0M percentile"};
    const AxisSpec axisMultPV{200, 0, 200, "N_{PV contrib}"};
    const AxisSpec axisRho{100, 0, 10, "#rho_{UE}"};
    const AxisSpec axisNEle{20, -0.5, 19.5, "N rows"};
    const AxisSpec axisNTrk{200, -0.5, 199.5, "N rows"};
    registry.add("hEventCounterInit", "collisions after |zvtx| (before Zorro)", kTH1D, {axisCounter});
    registry.add("hEMCalTrigger", "EMCal trigger", kTH1D, {axisTrigger});
    registry.add("hEventCounter", "collisions after Zorro (normalisation)", kTH1D, {axisCounter});
    registry.add("hEventStored", "collisions written to the derived data", kTH1D, {axisCounter});
    const AxisSpec axisSkim{3, -0.5, 2.5, "skim step"};
    registry.add("hSkimCutflow", "0: after Zorro, 1: EMCal-matched track > ptSkimMin, 2: + isolated electron", kTH1D, {axisSkim});
    registry.add("hZvtx", "Z vertex (after Zorro)", kTH1D, {axisZvtx});
    registry.add("hCentFT0M", "FT0M percentile (after Zorro)", kTH1D, {axisCent});
    registry.add("hMultPV", "PV contributors (after Zorro)", kTH1D, {axisMultPV});
    registry.add("hRho", "UE density", kTH1D, {axisRho});
    registry.add("hNTracksStored", "stored HfWbTracks rows per event", kTH1D, {axisNTrk});
    registry.add("hNElectronsStored", "stored HfWbElectrons rows per event", kTH1D, {axisNEle});
  }

  void initRun(aod::BCsWithTimestamps::iterator const& bc)
  {
    const int runNumber = bc.runNumber();
    if (runNumber == lastRunNumber) {
      return;
    }
    const uint64_t timestamp = bc.timestamp();
    if (cfgSkimmedProcessing) {
      LOGF(info, "Initializing Zorro for run %d", runNumber);
      zorro.initCCDB(ccdb.service, runNumber, timestamp, cfgTriggerName);
    }
    // magnetic field: needed by the analysis for KFParticle, independent of Zorro
    auto* grpo = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(ccdbPathGrpMag, timestamp);
    o2::base::Propagator::initFieldFromGRP(grpo);
    magneticField = static_cast<float>(o2::base::Propagator::Instance()->getNominalBz());
    LOG(info) << "run " << runNumber << " : magneticField = " << magneticField;
    lastRunNumber = runNumber;
  }

  float estimateRhoRC(Clusters const& clusters)
  {
    const float area = o2::constants::math::PI * rIsolation * rIsolation;
    std::vector<float> sumErc;
    sumErc.reserve(nRandomCones);
    for (int i = 0; i < nRandomCones; i++) {
      const float etarc = rnd->Uniform(-etaEmcMax, etaEmcMax);
      const float phirc = rnd->Uniform(phiEmcMin, phiEmcMax);
      float sum = 0.f;
      for (const auto& c : clusters) {
        if (c.definition() != clusterDefinition || c.energy() > rcHardE) {
          continue;
        }
        const double dEta = etarc - c.eta();
        const double dPhi = RecoDecay::constrainAngle(phirc - c.phi(), -o2::constants::math::PI);
        if (std::sqrt(dEta * dEta + dPhi * dPhi) < rIsolation) {
          sum += c.energy();
        }
      }
      sumErc.push_back(sum);
    }
    if (sumErc.empty()) {
      return 0.f;
    }
    std::nth_element(sumErc.begin(), sumErc.begin() + sumErc.size() / 2, sumErc.end());
    return sumErc[sumErc.size() / 2] / area;
  }

  // energy of the clusters in the cone around the given cluster, the cluster itself excluded
  float coneEnergyExcl(float eta, float phi, float energy, Clusters const& clusters)
  {
    double sum = 0.;
    for (const auto& c : clusters) {
      if (c.definition() != clusterDefinition) {
        continue;
      }
      const double dEta = c.eta() - eta;
      const double dPhi = RecoDecay::constrainAngle(c.phi() - phi, -o2::constants::math::PI);
      if (std::sqrt(dEta * dEta + dPhi * dPhi) < rIsolation) {
        sum += c.energy();
      }
    }
    return static_cast<float>(sum - energy);
  }

  template <typename TMatches>
  HfEmcMatch findClosestCluster(TMatches const& matches, Clusters const& clusters)
  {
    HfEmcMatch res;
    double rMin = 1.e9;
    for (const auto& match : matches) {
      const auto& cluster = match.template emcalcluster_as<Clusters>();
      if (cluster.definition() != clusterDefinition) {
        continue;
      }
      if (cluster.time() < timeEmcMin || cluster.time() > timeEmcMax) {
        continue;
      }
      const auto& trk = match.template track_as<TrackEle>();
      const double dEta = trk.trackEtaEmcal() - cluster.eta();
      const double dPhi = RecoDecay::constrainAngle(trk.trackPhiEmcal() - cluster.phi(), -o2::constants::math::PI);
      if (std::abs(dPhi) <= rMatchQa && std::abs(dEta) <= rMatchQa && res.nMatchInR < UINT8_MAX) {
        res.nMatchInR++;
      }
      const double r = std::sqrt(dEta * dEta + dPhi * dPhi);
      if (r < rMin) {
        rMin = r;
        res.has = true;
        res.dEta = dEta;
        res.dPhi = dPhi;
        res.energy = cluster.energy();
        res.eta = cluster.eta();
        res.phi = cluster.phi();
        res.m02 = cluster.m02();
        res.time = cluster.time();
        res.nCells = static_cast<int16_t>(cluster.nCells());
      }
    }
    if (res.has) {
      res.coneEnergy = coneEnergyExcl(res.eta, res.phi, res.energy, clusters);
    }
    return res;
  }

  // isolated-electron candidate used by the skim. Same definitions as in the analysis:
  // E/p with the closest cluster, M02, TPC n sigma, cluster isolation (E_cone,excl - E_UE) / E.
  template <typename TTrack>
  bool isIsolatedElectronSkim(TTrack const& track, HfEmcMatch const& emc, float energyUE)
  {
    if (emc.energy <= 0.f) {
      return false;
    }
    const float eop = emc.energy / track.p();
    if (eop < skimCfg.eopSkimMin || eop > skimCfg.eopSkimMax) {
      return false;
    }
    if (emc.m02 < skimCfg.m02SkimMin || emc.m02 > skimCfg.m02SkimMax) {
      return false;
    }
    if (track.tpcNSigmaEl() < skimCfg.nsigTpcSkimMin || track.tpcNSigmaEl() > skimCfg.nsigTpcSkimMax) {
      return false;
    }
    const float isoEnergy = (emc.coneEnergy - energyUE) / emc.energy;
    return isoEnergy < skimCfg.isoEnergySkimMax;
  }

  void processData(soa::Filtered<CollisionsWithCent>::iterator const& collision,
                   aod::BCsWithTimestamps const&,
                   Clusters const& clusters,
                   TrackEle const& tracks,
                   o2::aod::EMCALMatchedTracks const& matchedTracks)
  {
    registry.fill(HIST("hEventCounterInit"), 0.5);
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    initRun(bc);

    if (cfgSkimmedProcessing) {
      const bool isTriggered = zorro.isSelected(bc.globalBC());
      registry.fill(HIST("hEMCalTrigger"), isTriggered ? 1 : 0);
      if (!isTriggered) {
        return;
      }
    }
    if (applySel8 && !collision.sel8()) {
      return;
    }
    // normalisation of the derived sample: every event that the analysis WOULD have seen on AO2D
    registry.fill(HIST("hEventCounter"), 0.5);
    registry.fill(HIST("hZvtx"), collision.posZ());
    registry.fill(HIST("hCentFT0M"), collision.centFT0M());
    registry.fill(HIST("hMultPV"), collision.multNTracksPV());

    // ---- pass 1 : electron candidates (needs the EMCal matching for the skim decision)
    struct HfEleCand {
      int64_t trackGlobalIndex = -1;
      float pSum = 0.f;
      int16_t nTrk = 0;
      HfEmcMatch emc{};
    };
    std::vector<HfEleCand> eleCands;
    bool hasMatchedTrack = false; // skim step 1
    bool hasIsoElectron = false;  // skim step 2

    // UE density: needed for the isolation in the skim, and stored for the analysis
    const float rho = estimateRhoRC(clusters);
    const float energyUESkim = skimCfg.useUEsubSkim ? rho * static_cast<float>(o2::constants::math::PI * rIsolation * rIsolation) : 0.f;
    registry.fill(HIST("hSkimCutflow"), 0);

    for (const auto& track : tracks) {
      if (track.pt() < ptEleStoreMin) {
        continue;
      }
      // track isolation: over ALL tracks of the collision, as getIsolatedTrack() in the task
      double pSum = 0.;
      int nTrk = 0;
      for (const auto& other : tracks) {
        const double dEta = other.eta() - track.eta();
        const double dPhi = RecoDecay::constrainAngle(other.phi() - track.phi(), -o2::constants::math::PI);
        if (std::sqrt(dEta * dEta + dPhi * dPhi) < rIsolation) {
          ++nTrk;
          pSum += other.p();
        }
      }
      HfEleCand cand{.trackGlobalIndex = track.globalIndex(), .pSum = static_cast<float>(pSum), .nTrk = static_cast<int16_t>(nTrk), .emc = HfEmcMatch{}};
      auto matches = matchedTracks.sliceBy(perTrackMatched, track.globalIndex());
      if (matches.size() != 0) {
        cand.emc = findClosestCluster(matches, clusters);
      }
      if (cand.emc.has && track.pt() > skimCfg.ptSkimMin && std::abs(cand.emc.dEta) <= skimCfg.rMatchSkim && std::abs(cand.emc.dPhi) <= skimCfg.rMatchSkim) {
        hasMatchedTrack = true;
        if (isIsolatedElectronSkim(track, cand.emc, energyUESkim)) {
          hasIsoElectron = true;
        }
      }
      eleCands.push_back(cand);
    }
    if (!hasMatchedTrack) {
      return;
    }
    registry.fill(HIST("hSkimCutflow"), 1);
    if (skimCfg.skimRequireIsoElectron && !hasIsoElectron) {
      return;
    }
    if (hasIsoElectron) {
      registry.fill(HIST("hSkimCutflow"), 2);
    }

    // ---- write the collision
    registry.fill(HIST("hRho"), rho);
    rowColl(collision.posZ(), collision.sel8(),
            collision.centFT0A(), collision.centFT0C(), collision.centFT0M(), collision.centFV0A(),
            collision.multFT0M(), collision.multNTracksPV(),
            bc.runNumber(), magneticField, rho);
    const auto collIdx = rowColl.lastIndex();
    registry.fill(HIST("hEventStored"), 0.5);

    // ---- pass 2 : tracks, and right after each electron-candidate track its HfWbElectrons row.
    // eleCands was filled in the same track order with the same pT condition, so it is read sequentially.
    int nTracksStored = 0;
    size_t iEle = 0;
    for (const auto& track : tracks) {
      const bool isAssociated = track.isGlobalTrackWoPtEta() && std::abs(track.eta()) <= etaTrMax && track.pt() > ptAssMin;
      const bool isEleCand = track.pt() >= ptEleStoreMin;
      if (!isAssociated && !isEleCand) {
        continue;
      }
      rowTrack(collIdx, track.pt(), track.eta(), track.phi(), static_cast<int8_t>(track.sign()),
               track.dcaXY(), track.dcaZ(),
               static_cast<int16_t>(track.tpcNClsFound()), static_cast<int16_t>(track.tpcNClsCrossedRows()),
               static_cast<uint8_t>(track.itsNCls()), track.tpcChi2NCl(), track.itsChi2NCl(),
               track.tpcNSigmaEl(), track.tpcNSigmaPi(), static_cast<bool>(track.isGlobalTrackWoPtEta()));
      ++nTracksStored;
      if (!isEleCand) {
        continue;
      }
      const auto& cand = eleCands.at(iEle++);
      if (cand.trackGlobalIndex != track.globalIndex()) {
        LOG(fatal) << "electron candidate out of sync: " << cand.trackGlobalIndex << " vs " << track.globalIndex();
      }
      const auto& e = cand.emc;
      rowEle(collIdx, rowTrack.lastIndex(),
             track.x(), track.alpha(), track.y(), track.z(), track.snp(), track.tgl(), track.signed1Pt(),
             track.cYY(), track.cZY(), track.cZZ(),
             track.cSnpY(), track.cSnpZ(), track.cSnpSnp(),
             track.cTglY(), track.cTglZ(), track.cTglSnp(), track.cTglTgl(),
             track.c1PtY(), track.c1PtZ(), track.c1PtSnp(), track.c1PtTgl(), track.c1Pt21Pt2(),
             cand.pSum, cand.nTrk,
             e.has, e.dEta, e.dPhi, e.energy, e.eta, e.phi, e.m02, e.time, e.nCells, e.nMatchInR, e.coneEnergy);
    }
    registry.fill(HIST("hNTracksStored"), nTracksStored);
    registry.fill(HIST("hNElectronsStored"), eleCands.size());
  }
  PROCESS_SWITCH(HfDerivedDataCreatorElectronWeakBoson, processData, "produce derived data for W/Z -> e (+h)", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfDerivedDataCreatorElectronWeakBoson>(cfgc)};
}
