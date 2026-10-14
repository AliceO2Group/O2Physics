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

/// \file emcalMcTask.cxx
/// \brief Analysis task for to obtain cluster properties from MC identified particles like photons and electrons
/// \author M. Hemmer, marvin.hemmer@cern.ch

#include "PWGEM/PhotonMeson/Core/EMBitFlags.h"
#include "PWGEM/PhotonMeson/Core/EMCPhotonCut.h"
#include "PWGEM/PhotonMeson/Core/EMPhotonEventCut.h"
#include "PWGEM/PhotonMeson/DataModel/EventTables.h"
#include "PWGEM/PhotonMeson/DataModel/GammaTablesRedux.h"
#include "PWGEM/PhotonMeson/Utils/EventHistograms.h"

#include <CommonConstants/PhysicsConstants.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Concepts.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/SliceCache.h>
#include <Framework/runDataProcessing.h>

#include <TPDGCode.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <string>
#include <string_view>
#include <vector>

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::aod::pwgem::photon;

namespace o2::em::emcal::mc
{

enum ParticleType : int {
  kPhoton = 0,
  kElectron,
  kPositron,
  kPi0,
  kEta,
  kOmega,
  kPion,
  kKaon,
  kOther,
  kNParticleTypes
};
} // namespace o2::em::emcal::mc

enum CentralityEstimator {
  None = 0,
  CFT0A = 1,
  CFT0C = 2,
  CFT0M = 3,
  NCentralityEstimators = 4
};

enum class MapLevel {
  kGood = 1,
  kNoBad = 2,
  kInEMC = 3,
  kAll = 4
};

struct EmcalMcTask {
  // configurable axis
  ConfigurableAxis thnConfigAxisE{"thnConfigAxisE", {400, 0., 20.}, "energy axis"};
  ConfigurableAxis thnConfigAxisEtaDiff{"thnConfigAxisEtaDiff", {300, -1., 2.}, "(eta rec - eta true)"};
  ConfigurableAxis thnConfigAxisPhiDiff{"thnConfigAxisPhiDiff", {300, -1., 2.}, "(phi rec - phi true"};
  ConfigurableAxis thnConfigAxisM02{"thnConfigAxisM02", {100, 0.0, 1.0}, "m02 mass axis"};
  ConfigurableAxis thnConfigAxisCent{"thnConfigAxisCent", {20, 0., 100.}, "centrality axis for the current event"};
  ConfigurableAxis thnConfigAxisMult{"thnConfigAxisMult", {60, 0., 60000.}, "multiplicity axis for the current event"};
  Configurable<bool> useCent{"useCent", false, "flag to enable usage of centrality instead of multiplicity as axis."};

  EMPhotonEventCut fEMEventCut;
  struct : ConfigurableGroup {
    std::string prefix = "eventcuts";
    Configurable<float> cfgZvtxMax{"cfgZvtxMax", 10.f, "max. Zvtx"};
    Configurable<bool> cfgRequireSel8{"cfgRequireSel8", true, "require sel8 in event cut"};
    Configurable<bool> cfgRequireFT0AND{"cfgRequireFT0AND", true, "require FT0AND in event cut"};
    Configurable<bool> cfgRequireNoTFB{"cfgRequireNoTFB", false, "require No time frame border in event cut"};
    Configurable<bool> cfgRequireNoITSROFB{"cfgRequireNoITSROFB", false, "require no ITS readout frame border in event cut"};
    Configurable<bool> cfgRequireNoSameBunchPileup{"cfgRequireNoSameBunchPileup", false, "require no same bunch pileup in event cut"};
    Configurable<bool> cfgRequireVertexITSTPC{"cfgRequireVertexITSTPC", false, "require Vertex ITSTPC in event cut"}; // ITS-TPC matched track contributes PV.
    Configurable<bool> cfgRequireGoodZvtxFT0vsPV{"cfgRequireGoodZvtxFT0vsPV", false, "require good Zvtx between FT0 vs. PV in event cut"};
    Configurable<bool> cfgRequireEMCReadoutInMB{"cfgRequireEMCReadoutInMB", true, "require the EMC to be read out in an MB collision (kTVXinEMC)"};
    Configurable<bool> cfgRequireEMCHardwareTriggered{"cfgRequireEMCHardwareTriggered", false, "require the EMC to be hardware triggered (kEMC7 or kDMC7)"};
    Configurable<float> cfgFT0COccupancyMin{"cfgFT0COccupancyMin", -1, "min. FT0C occupancy"};
    Configurable<float> cfgFT0COccupancyMax{"cfgFT0COccupancyMax", 1000000000, "max. FT0C occupancy"};
    Configurable<float> cfgMinCent{"cfgMinCent", 0, "min. centrality (%)"};
    Configurable<float> cfgMaxCent{"cfgMaxCent", 90, "max. centrality (%)"};
    Configurable<int> centEstimator{"centEstimator", 2, "Centrality estimation (FT0A: 1, FT0C: 2, FT0M: 3)"};
  } eventcuts;

  EMCPhotonCut fEMCCut;
  struct : ConfigurableGroup {
    std::string prefix = "emccuts";
    Configurable<std::string> clusterDefinition{"clusterDefinition", "kV3MostSplitSmallestTimeDiff", "Clusterizer to be selected, e.g. V3Default"};
    Configurable<float> cfgEMCminTime{"cfgEMCminTime", -25., "Minimum cluster time for EMCal time cut"};
    Configurable<float> cfgEMCmaxTime{"cfgEMCmaxTime", +30., "Maximum cluster time for EMCal time cut"};
    Configurable<float> cfgEMCminM02{"cfgEMCminM02", 0.1, "Minimum M02 for EMCal M02 cut"};
    Configurable<float> cfgEMCmaxM02{"cfgEMCmaxM02", 0.7, "Maximum M02 for EMCal M02 cut"};
    Configurable<float> cfgEMCminE{"cfgEMCminE", 0.7, "Minimum cluster energy for EMCal energy cut"};
    Configurable<int> cfgEMCminNCell{"cfgEMCminNCell", 1, "Minimum number of cells per cluster for EMCal NCell cut"};
    Configurable<std::vector<float>> cfgEMCTMEta{"cfgEMCTMEta", {0.01f, 4.07f, -2.5f}, "|eta| <= [0]+(pT+[1])^[2] for EMCal track matching"};
    Configurable<std::vector<float>> cfgEMCTMPhi{"cfgEMCTMPhi", {0.015f, 3.65f, -2.f}, "|phi| <= [0]+(pT+[1])^[2] for EMCal track matching"};
    Configurable<std::vector<float>> emcSecTMEta{"emcSecTMEta", {0.01f, 4.07f, -2.5f}, "|eta| <= [0]+(pT+[1])^[2] for EMCal track matching"};
    Configurable<std::vector<float>> emcSecTMPhi{"emcSecTMPhi", {0.015f, 3.65f, -2.f}, "|phi| <= [0]+(pT+[1])^[2] for EMCal track matching"};
    Configurable<float> cfgEMCEoverp{"cfgEMCEoverp", 1.75, "Minimum cluster energy over track momentum for EMCal track matching"};
    Configurable<bool> cfgEMCUseExoticCut{"cfgEMCUseExoticCut", true, "FLag to use the EMCal exotic cluster cut"};
    Configurable<bool> cfgEMCUseTM{"cfgEMCUseTM", false, "flag to use EMCal track matching cut or not"};
    Configurable<bool> emcUseSecondaryTM{"emcUseSecondaryTM", false, "flag to use EMCal secondary track matching cut or not"};
    Configurable<bool> cfgEnableQA{"cfgEnableQA", false, "flag to turn QA plots on/off"};
  } emccuts;

  SliceCache cache;

  using EMCalPhotons = soa::Join<aod::EMCEMEventIds, aod::MinClusters, aod::EMEMCClusterMCLabels_001>;

  using Colls = soa::Join<aod::PMEvents, aod::EMEventsAlias, aod::EMEventsMult_000, aod::EMEventsCent_000, aod::EMMCEventLabels, aod::EmMagFields>;

  using McColls = o2::soa::Join<o2::aod::EMMCEvents, o2::aod::BinnedGenPts>;
  using McParticles = EMMCParticles;

  PresliceOptional<EMCalPhotons> perCollisionEMC = o2::aod::emccluster::pmeventId;
  PresliceOptional<MinMTracks> perEMCClusterMT = o2::aod::mintm::minClusterId;
  PresliceOptional<MinMSTracks> perEMCClusterMS = o2::aod::mintm::minClusterId;

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};

  int mRunNumber{-1};
  float dBz{0.f};

  static constexpr std::array<std::string_view, static_cast<size_t>(o2::em::emcal::mc::ParticleType::kNParticleTypes)> kSubDirs = {
    "photon/", "electron/", "positron/", "pi0/",
    "eta/", "omega/", "pion/", "kaon/", "other/"};

  void defineEMEventCut()
  {
    fEMEventCut = EMPhotonEventCut("fEMEventCut", "fEMEventCut");
    fEMEventCut.SetRequireSel8(eventcuts.cfgRequireSel8);
    fEMEventCut.SetRequireFT0AND(eventcuts.cfgRequireFT0AND);
    fEMEventCut.SetZvtxRange(-eventcuts.cfgZvtxMax, +eventcuts.cfgZvtxMax);
    fEMEventCut.SetRequireNoTFB(eventcuts.cfgRequireNoTFB);
    fEMEventCut.SetRequireNoITSROFB(eventcuts.cfgRequireNoITSROFB);
    fEMEventCut.SetRequireNoSameBunchPileup(eventcuts.cfgRequireNoSameBunchPileup);
    fEMEventCut.SetRequireVertexITSTPC(eventcuts.cfgRequireVertexITSTPC);
    fEMEventCut.SetRequireGoodZvtxFT0vsPV(eventcuts.cfgRequireGoodZvtxFT0vsPV);
    fEMEventCut.SetRequireEMCReadoutInMB(eventcuts.cfgRequireEMCReadoutInMB);
    fEMEventCut.SetRequireEMCHardwareTriggered(eventcuts.cfgRequireEMCHardwareTriggered);
  }

  void defineEMCCut()
  {
    fEMCCut = EMCPhotonCut("fEMCCut", "fEMCCut");

    fEMCCut.SetTrackMatchingEtaParams(emccuts.cfgEMCTMEta->at(0), emccuts.cfgEMCTMEta->at(1), emccuts.cfgEMCTMEta->at(2));
    fEMCCut.SetTrackMatchingPhiParams(emccuts.cfgEMCTMPhi->at(0), emccuts.cfgEMCTMPhi->at(1), emccuts.cfgEMCTMPhi->at(2));

    fEMCCut.SetSecTrackMatchingEtaParams(emccuts.emcSecTMEta->at(0), emccuts.emcSecTMEta->at(1), emccuts.emcSecTMEta->at(2));
    fEMCCut.SetSecTrackMatchingPhiParams(emccuts.emcSecTMPhi->at(0), emccuts.emcSecTMPhi->at(1), emccuts.emcSecTMPhi->at(2));
    fEMCCut.SetMinEoverP(emccuts.cfgEMCEoverp);

    fEMCCut.SetMinE(emccuts.cfgEMCminE);
    fEMCCut.SetMinNCell(emccuts.cfgEMCminNCell);
    fEMCCut.SetM02Range(emccuts.cfgEMCminM02, emccuts.cfgEMCmaxM02);
    fEMCCut.SetTimeRange(emccuts.cfgEMCminTime, emccuts.cfgEMCmaxTime);
    fEMCCut.SetUseExoticCut(emccuts.cfgEMCUseExoticCut);
    fEMCCut.SetClusterizer(emccuts.clusterDefinition);
    fEMCCut.SetUseTM(emccuts.cfgEMCUseTM.value);                // disables or enables TM
    fEMCCut.SetUseSecondaryTM(emccuts.emcUseSecondaryTM.value); // disables or enables secondary TM
    fEMCCut.SetDoQA(emccuts.cfgEnableQA.value);
  }

  void init(InitContext&)
  {
    mRunNumber = 0;
    dBz = 0;

    defineEMEventCut();
    defineEMCCut();
    fEMCCut.addQAHistograms(&registry);
    o2::aod::pwgem::photonmeson::utils::eventhistogram::addEventHistograms(&registry);

    const AxisSpec thnAxisERec{thnConfigAxisE, "#it{E}_{Rec} (GeV)"};

    const AxisSpec thnAxisM02{thnConfigAxisM02, "#it{M}_{02}"};

    const AxisSpec thnAxisEtaDiff{thnConfigAxisEtaDiff, "#it{#eta}_{Rec} - #it{#eta}_{Gen}"};
    const AxisSpec thnAxisPhiDiff{thnConfigAxisPhiDiff, "#it{#varphi}_{Rec} - #it{#varphi}_{Gen}"};

    AxisSpec thnAxisCentOrMult{1, 0., 1., "Centrality/Multiplicity"}; // placeholder, overwritten in init
    if (useCent.value) {
      // PbPb: use centrality
      thnAxisCentOrMult = {thnConfigAxisCent, "Centrality (%)"};
    } else {
      // pp: use multiplicity
      thnAxisCentOrMult = {thnConfigAxisMult, "FT0C Multiplicity"};
    }

    registry.add("photon/hM02", "cluster m02 vs energy vs cent/mult", HistType::kTH3F, {thnAxisM02, thnAxisERec, thnAxisCentOrMult});
    registry.add("photon/hEtaRel", "relative #eta vs energy vs cent/mult", HistType::kTH3F, {thnAxisEtaDiff, thnAxisERec, thnAxisCentOrMult});
    registry.add("photon/hPhiRel", "relative #varphi vs energy vs cent/mult", HistType::kTH3F, {thnAxisPhiDiff, thnAxisERec, thnAxisCentOrMult});

    registry.addClone("photon/", "electron/");
    registry.addClone("photon/", "positron/");
    registry.addClone("photon/", "pi0/");
    registry.addClone("photon/", "eta/");
    registry.addClone("photon/", "omega/");
    registry.addClone("photon/", "pion/");
    registry.addClone("photon/", "kaon/");
    registry.addClone("photon/", "other/");

  }; // end init

  template <o2::soa::is_iterator TCollision>
  float getCentralityOrMultiplicity(TCollision const& collision)
  {
    if (useCent.value) {
      return getCentrality(collision);
    }
    // pp: use raw FT0C multiplicity
    return collision.multFT0C();
  }

  /// Get the centrality
  /// \param collision is the collision with the centrality information
  template <o2::soa::is_iterator TCollision>
  float getCentrality(TCollision const& collision)
  {
    float cent = -999.;
    switch (eventcuts.centEstimator) {
      case CentralityEstimator::CFT0M:
        cent = collision.centFT0M();
        break;
      case CentralityEstimator::CFT0A:
        cent = collision.centFT0A();
        break;
      case CentralityEstimator::CFT0C:
        cent = collision.centFT0C();
        break;
      default:
        LOG(warning) << "Centrality estimator not valid. Possible values are T0M, T0A, T0C. Fallback to T0C";
        cent = collision.centFT0C();
        break;
    }
    return cent;
  }

  /// \brief check if standard event cuts + FT0 occupancy + centrality + QVec good is
  /// \param collision collision that will be checked
  /// \return true if collision survives all checks, otherwise false
  template <o2::soa::is_iterator TCollision>
  bool isFullEventSelected(TCollision const& collision, bool fillHisto = false)
  {
    if (fillHisto) {
      o2::aod::pwgem::photonmeson::utils::eventhistogram::fillEventInfo<0>(&registry, collision);
    }
    if (!(fEMEventCut.IsSelected(collision))) {
      // general event selection
      return false;
    }
    if (!(eventcuts.cfgFT0COccupancyMin <= collision.ft0cOccupancyInTimeRange() && collision.ft0cOccupancyInTimeRange() < eventcuts.cfgFT0COccupancyMax)) {
      // occupancy selection
      return false;
    }
    float centOrMult = getCentralityOrMultiplicity(collision);
    if (useCent && (centOrMult < eventcuts.cfgMinCent || centOrMult > eventcuts.cfgMaxCent)) {
      // event selection
      return false;
    }
    if (fillHisto) {
      o2::aod::pwgem::photonmeson::utils::eventhistogram::fillEventInfo<1>(&registry, collision);
      registry.fill(HIST("Event/before/hCollisionCounter"), 12.0); // accepted
      registry.fill(HIST("Event/after/hCollisionCounter"), 12.0);  // accepted
    }
    return true;
  }

  // One templated fill function instead of 9 copy-pasted blocks
  template <const int type, o2::soa::is_iterator TCluster, o2::soa::is_iterator TMC>
  void fillClusterHistos(TCluster const& cluster, TMC const& mcPart, float centOrMult)
  {
    static constexpr std::string_view SubDir = kSubDirs[type];

    registry.fill(HIST(SubDir) + HIST("hM02"), cluster.m02(), cluster.e(), centOrMult);
    registry.fill(HIST(SubDir) + HIST("hEtaRel"), cluster.eta() - mcPart.eta(), cluster.e(), centOrMult);
    registry.fill(HIST(SubDir) + HIST("hPhiRel"), cluster.phi() - mcPart.phi(), cluster.e(), centOrMult);
  }

  // PCM-EMCal same event
  void processEmcal(Colls const& collisions, EMCalPhotons const& clusters, MinMTracks const& matchedPrims, MinMSTracks const& matchedSeconds, EMMCParticles const& mcParticles)
  {
    if (clusters.size() <= 0) {
      LOG(info) << "Skipping DF because there are not photons!";
      return;
    }
    EMBitFlags emcFlags(clusters.size());
    if (clusters.size() > 0) {
      fEMCCut.AreSelectedRunning(emcFlags, clusters, matchedPrims, matchedSeconds, &registry);
    }

    // create iterators for photon mc particles
    auto mcPhoton1 = mcParticles.begin();

    for (const auto& collision : collisions) {
      isFullEventSelected(collision, true);

      float centOrMult = getCentralityOrMultiplicity(collision);

      auto photonsEMCPerCollision = clusters.sliceBy(perCollisionEMC, collision.globalIndex());

      for (const auto& photonEMC : photonsEMCPerCollision) {
        if (!(emcFlags.test(photonEMC.globalIndex()))) {
          continue;
        }
        if (photonEMC.emmcparticleIds().empty()) {
          // this is a cluster with just noise, skip
          continue;
        }
        // we only want to look at the largest contribution
        mcPhoton1.setCursor(photonEMC.emmcparticleIds()[0]);

        if (std::abs(mcPhoton1.pdgCode()) == PDG_t::kGamma) {
          fillClusterHistos<o2::em::emcal::mc::ParticleType::kPhoton>(photonEMC, mcPhoton1, centOrMult);
        } else if (mcPhoton1.pdgCode() == PDG_t::kElectron) {
          fillClusterHistos<o2::em::emcal::mc::ParticleType::kElectron>(photonEMC, mcPhoton1, centOrMult);
        } else if (mcPhoton1.pdgCode() == PDG_t::kPositron) {
          fillClusterHistos<o2::em::emcal::mc::ParticleType::kPositron>(photonEMC, mcPhoton1, centOrMult);
        } else if (std::abs(mcPhoton1.pdgCode()) == PDG_t::kPi0) {
          fillClusterHistos<o2::em::emcal::mc::ParticleType::kPi0>(photonEMC, mcPhoton1, centOrMult);
        } else if (std::abs(mcPhoton1.pdgCode()) == o2::constants::physics::Pdg::kEta) {
          fillClusterHistos<o2::em::emcal::mc::ParticleType::kEta>(photonEMC, mcPhoton1, centOrMult);
        } else if (std::abs(mcPhoton1.pdgCode()) == o2::constants::physics::Pdg::kOmega) {
          fillClusterHistos<o2::em::emcal::mc::ParticleType::kOmega>(photonEMC, mcPhoton1, centOrMult);
        } else if (std::abs(mcPhoton1.pdgCode()) == PDG_t::kPiPlus) {
          fillClusterHistos<o2::em::emcal::mc::ParticleType::kPion>(photonEMC, mcPhoton1, centOrMult);
        } else if (std::abs(mcPhoton1.pdgCode()) == PDG_t::kKPlus) {
          fillClusterHistos<o2::em::emcal::mc::ParticleType::kKaon>(photonEMC, mcPhoton1, centOrMult);
        } else {
          fillClusterHistos<o2::em::emcal::mc::ParticleType::kOther>(photonEMC, mcPhoton1, centOrMult);
        }
      } // for (const auto& photonEMC : photonsEMCPerCollision) {
    }
  }
  PROCESS_SWITCH(EmcalMcTask, processEmcal, "Process for emcal", true);

}; // End struct EmcalMcTask

WorkflowSpec defineDataProcessing(ConfigContext const& context)
{
  return WorkflowSpec{adaptAnalysisTask<EmcalMcTask>(context)};
}
