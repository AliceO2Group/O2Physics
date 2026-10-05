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

/// \file enegyFlowGfw.cxx
/// \brief Charged-particle transverse-energy correlators in Xe-Xe and Pb-Pb collisions
/// \author Emil Gorm Dahlbæk Nielsen <emil.gorm.nielsen@cern.ch>

#include "PWGCF/GenericFramework/Core/FlowContainer.h"
#include "PWGCF/GenericFramework/Core/GFW.h"
#include "PWGCF/GenericFramework/Core/GFWConfig.h"
#include "PWGCF/GenericFramework/Core/GFWWeights.h"

#include "Common/CCDB/TriggerAliases.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/Expressions.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <TFile.h>
#include <TH1.h>
#include <THn.h>
#include <TNamed.h>
#include <TObjArray.h>
#include <TRandom3.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <string>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::analysis::genericframework;

struct EnergyFlowGfw {
  enum EnergyWeightMode : int {
    PtWeight = 0,
    MtWeight,
    FullEtWeight
  };

  enum ParticleId : int {
    Unidentified = 0,
    Pion,
    Kaon,
    Proton
  };

  enum ActivityAxisMode : int {
    CentralityAxis = 0,
    NGlobalTracksAxis,
    ZdcSpectatorAxis
  };

  enum AodSelectionStep : uint8_t {
    AllAodEvents = 0,
    TriggerSelected,
    CentralitySelected,
    HasValidZdc,
    SpectatorSelected
  };

  enum ZdcQaStep : uint8_t {
    MatchedZdc = 0,
    FiniteCommonEnergy,
    NonnegativeCommonEnergy,
    HasZnTiming,
    HasZpTiming,
    HasAllTiming
  };

  struct : ConfigurableGroup {
    Configurable<float> vertexZ{"vertexZ", 10.f, "Accepted |z-vertex| (cm)"};
    Configurable<std::pair<float, float>> ptCuts{"ptCuts", {0.2f, 5.f}, "Track pT cuts (GeV/c)"};
    Configurable<float> etaMax{"etaMax", 0.8f, "Maximum absolute track pseudorapidity"};
  } cuts;

  struct : ConfigurableGroup {
    Configurable<int> weightMode{"weightMode", FullEtWeight, "Energy weight: 0=pT, 1=mT, 2=E/cosh(eta)"};
    Configurable<float> nSigmaCut{"nSigmaCut", 3.f, "TPC or combined TPC+TOF PID selection in n-sigma"};
    Configurable<float> tofPtCut{"tofPtCut", 0.5f, "Require TOF PID above this pT (GeV/c)"};
    Configurable<bool> strictPid{"strictPid", true, "Mark tracks compatible with multiple species as unidentified"};
    Configurable<float> unidentifiedMass{"unidentifiedMass", o2::constants::physics::MassPionCharged, "Mass assigned to unidentified charged tracks (GeV/c^2)"};
  } energy;

  Configurable<GFWRegions> cfgRegions{"cfgRegions", {{"refA", "refB"}, {-0.8, 0.4}, {-0.4, 0.8}, {0, 0}, {1, 1}}, "GFW region definitions"};
  Configurable<GFWCorrConfigs> cfgCorrConfig{"cfgCorrConfig", {{"refA {0}", "refB {0}", "refA {0} refB {0}", "refA {1} refB {-1}", "refA {2} refB {-2}", "refA {3} refB {-3}", "refA {4} refB {-4}", "refA {5} refB {-5}", "refA {6} refB {-6}"}, {"F0A", "F0B", "F0AF0B", "EECGap12", "EECGap22", "EECGap32", "EECGap42", "EECGap52", "EECGap62"}, {0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0}}, "GFW correlations to calculate"};
  Configurable<int> cfgBootstrapSamples{"cfgBootstrapSamples", 10, "Number of FlowContainer bootstrap subsamples; zero disables bootstrapping"};
  Configurable<int> cfgActivityAxis{"cfgActivityAxis", CentralityAxis, "FlowContainer x-axis: 0=centrality, 1=number of filtered global tracks, 2=corrected ZDC spectator estimate"};
  Configurable<float> cfgMaxCentrality{"cfgMaxCentrality", 100.f, "Maximum accepted centrality percentile"};

  Configurable<float> cfgSpectatorEnergyPerNucleon{"cfgSpectatorEnergyPerNucleon", 1.f, "ZDC energy corresponding to one spectator nucleon"};
  Configurable<float> cfgMaxTotalSpectators{"cfgMaxTotalSpectators", 500.f, "Maximum accepted spectator count from ZNA+ZNC+ZPA+ZPC"};
  Configurable<float> cfgZdcResponseZNA{"cfgZdcResponseZNA", 1.f, "ZNA efficiency/acceptance response factor alpha; corrected energy is E_ZNA/alpha_ZNA"};
  Configurable<float> cfgZdcResponseZNC{"cfgZdcResponseZNC", 1.f, "ZNC efficiency/acceptance response factor alpha; corrected energy is E_ZNC/alpha_ZNC"};
  Configurable<float> cfgZdcResponseZPA{"cfgZdcResponseZPA", 1.f, "ZPA efficiency/acceptance response factor alpha; corrected energy is E_ZPA/alpha_ZPA"};
  Configurable<float> cfgZdcResponseZPC{"cfgZdcResponseZPC", 1.f, "ZPC efficiency/acceptance response factor alpha; corrected energy is E_ZPC/alpha_ZPC"};
  Configurable<float> cfgZdcMaxAbsTimeForQa{"cfgZdcMaxAbsTimeForQa", 100.f, "Maximum absolute ZDC time considered available in QA; excludes missing-time sentinels"};

  Configurable<std::string> efficiencyPath{"efficiencyPath", "", "CCDB path or local ROOT file containing a 4D efficiency-correction THn"};
  Configurable<bool> efficiencyFromLocalFile{"efficiencyFromLocalFile", false, "Load efficiencyPath from a local ROOT file"};
  Configurable<std::string> acceptancePath{"acceptancePath", "", "CCDB path to GFW acceptance weights"};
  Configurable<bool> acceptanceRunByRun{"acceptanceRunByRun", false, "Load acceptance weights from the RunByRun subdirectory"};

  ConfigurableAxis axisVertex{"axisVertex", {10, -10., 10.}, "z-vertex (cm)"};
  ConfigurableAxis axisCentrality{"axisCentrality", {VARIABLE_WIDTH, 0., 5., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100.1}, "centrality"};
  ConfigurableAxis axisNGlobalTracks{"axisNGlobalTracks", {500, 0., 5000.}, "number of filtered global tracks"};
  ConfigurableAxis axisPt{"axisPt", {100, 0., 5.}, "pT (GeV/c)"};
  ConfigurableAxis axisEnergy{"axisEnergy", {200, 0., 20.}, "transverse-energy weight (GeV)"};
  ConfigurableAxis axisTotalSpectators{"axisTotalSpectators", {500, 0., 500.}, "total spectator count"};
  ConfigurableAxis axisSpectatorActivity{"axisSpectatorActivity", {500, 0., 500.}, "spectator count used by the FlowContainer"};
  ConfigurableAxis axisZdcNeutronEnergy{"axisZdcNeutronEnergy", {410, -2000., 80000.}, "ZDC neutron common energy"};
  ConfigurableAxis axisZdcProtonEnergy{"axisZdcProtonEnergy", {310, -1000., 30000.}, "ZDC proton common energy"};
  ConfigurableAxis axisZdcTotalEnergy{"axisZdcTotalEnergy", {410, -5000., 200000.}, "summed ZDC common energy"};
  ConfigurableAxis axisZdcZemEnergy{"axisZdcZemEnergy", {210, -100., 2000.}, "summed ZEM energy"};
  ConfigurableAxis axisZdcTime{"axisZdcTime", {600, -15., 15.}, "ZDC time (ns)"};
  ConfigurableAxis axisZdcTimeCombination{"axisZdcTimeCombination", {400, -20., 20.}, "ZDC time sum or difference (ns)"};
  ConfigurableAxis axisZdcTrackMultiplicity{"axisZdcTrackMultiplicity", {400, 0., 4000.}, "number of filtered global tracks"};

  Filter collisionFilter = nabs(aod::collision::posZ) < cuts.vertexZ;
  Filter trackFilter = (nabs(aod::track::eta) < cuts.etaMax) && (aod::track::pt > cuts.ptCuts->first) && (aod::track::pt < cuts.ptCuts->second) &&
                       ((requireGlobalTrackInFilter()) || (aod::track::isGlobalTrackSDD == static_cast<uint8_t>(true)));

  using CollisionsRun2 = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentRun2V0Ms>>;
  using CollisionsRun3 = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs>>;
  using Tracks = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection,
                                         aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr,
                                         aod::pidTOFbeta, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr>>;
  using BCsRun2 = soa::Join<aod::BCs, aod::Timestamps, aod::BcSels, aod::Run2MatchedToBCSparse>;
  using BCsRun3 = soa::Join<aod::BCs, aod::Timestamps, aod::BcSels, aod::Run3MatchedToBCSparse>;

  OutputObj<FlowContainer> flowContainer{"FlowContainer"};
  HistogramRegistry registry{"registry"};
  Service<o2::ccdb::BasicCCDBManager> ccdb{};

  std::unique_ptr<GFW> gfw{std::make_unique<GFW>()};
  std::unique_ptr<TRandom3> random{std::make_unique<TRandom3>(0)};
  std::vector<GFW::CorrConfig> corrConfigs;
  std::vector<bool> zeroHarmonicConfigs;
  int gfwMask = 0;

  std::unique_ptr<THn> localEfficiency;
  THn* efficiency = nullptr;
  GFWWeights* acceptance = nullptr;
  int acceptanceRunNumber = -1;

  void init(InitContext const&)
  {
    if (doprocessRun2 && doprocessRun3) {
      LOGF(fatal, "Enable only one input format: Run 2 AOD or Run 3 AOD");
    }
    if (cuts.ptCuts->first >= cuts.ptCuts->second || cuts.etaMax <= 0.f) {
      LOGF(fatal, "Require ptMin < ptMax and etaMax > 0");
    }
    if (energy.weightMode < PtWeight || energy.weightMode > FullEtWeight) {
      LOGF(fatal, "weightMode must be 0 (pT), 1 (mT), or 2 (full ET)");
    }
    if (energy.nSigmaCut <= 0.f || energy.tofPtCut < 0.f || energy.unidentifiedMass < 0.f) {
      LOGF(fatal, "PID n-sigma and mass settings must be non-negative, with nSigmaCut > 0");
    }
    if (cfgBootstrapSamples < 0) {
      LOGF(fatal, "cfgBootstrapSamples must not be negative");
    }
    if (cfgActivityAxis < CentralityAxis || cfgActivityAxis > ZdcSpectatorAxis) {
      LOGF(fatal, "cfgActivityAxis must be 0 (centrality), 1 (filtered global tracks), or 2 (corrected ZDC spectators)");
    }
    if (cfgMaxCentrality < 0.f || cfgMaxCentrality > 100.f || cfgSpectatorEnergyPerNucleon <= 0.f || cfgMaxTotalSpectators < 0.f || cfgZdcResponseZNA <= 0.f || cfgZdcResponseZNC <= 0.f || cfgZdcResponseZPA <= 0.f || cfgZdcResponseZPC <= 0.f || cfgZdcMaxAbsTimeForQa <= 0.f) {
      LOGF(fatal, "The maximum centrality must be in [0, 100], the spectator energy scale, ZDC response factors, and ZDC QA time range must be positive, and the spectator threshold non-negative");
    }

    const int numberOfRegions = cfgRegions->GetSize();
    if (numberOfRegions <= 0 || cfgCorrConfig->GetSize() <= 0) {
      LOGF(fatal, "GFW regions and correlations must be non-empty and internally consistent");
    }

    float largestRegionEta = 0.f;
    for (int index = 0; index < numberOfRegions; ++index) {
      const float etaMin = cfgRegions->GetEtaMin()[index];
      const float etaMax = cfgRegions->GetEtaMax()[index];
      if (!std::isfinite(etaMin) || !std::isfinite(etaMax) || etaMin >= etaMax) {
        LOGF(fatal, "GFW region %s must have finite etaMin < etaMax", cfgRegions->GetNames()[index].c_str());
      }
      if (cfgRegions->GetpTDifs()[index] != 0 || cfgRegions->GetBitmasks()[index] == 0) {
        LOGF(fatal, "This task supports only pT-integrated regions with non-zero bit masks");
      }
      largestRegionEta = std::max({largestRegionEta, std::abs(etaMin), std::abs(etaMax)});
      gfwMask |= cfgRegions->GetBitmasks()[index];
      gfw->AddRegion(cfgRegions->GetNames()[index], etaMin, etaMax, 1, cfgRegions->GetBitmasks()[index]);
    }
    if (largestRegionEta > cuts.etaMax) {
      LOGF(fatal, "All GFW regions must be contained within the selected track eta range");
    }

    auto profileNames = std::make_unique<TObjArray>();
    profileNames->SetOwner(true);
    corrConfigs.reserve(cfgCorrConfig->GetSize());
    zeroHarmonicConfigs.reserve(cfgCorrConfig->GetSize());
    for (int index = 0; index < cfgCorrConfig->GetSize(); ++index) {
      if (cfgCorrConfig->GetpTDifs()[index] != 0) {
        LOGF(fatal, "This task supports only pT-integrated GFW correlations");
      }
      auto config = gfw->GetCorrelatorConfig(cfgCorrConfig->GetCorrs()[index], cfgCorrConfig->GetHeads()[index], false);
      if (config.Regs.empty() || config.Head.empty()) {
        LOGF(fatal, "Invalid GFW correlation configuration at index %d", index);
      }
      const bool zeroHarmonic = std::all_of(config.Hars.begin(), config.Hars.end(), [](const auto& harmonics) {
        return std::all_of(harmonics.begin(), harmonics.end(), [](int harmonic) { return harmonic == 0; });
      });
      zeroHarmonicConfigs.push_back(zeroHarmonic);
      if (!zeroHarmonic) {
        profileNames->Add(new TNamed(config.Head.c_str(), config.Head.c_str()));
      }
      const std::string rawName = config.Head + "_raw";
      profileNames->Add(new TNamed(rawName.c_str(), rawName.c_str()));
      corrConfigs.push_back(std::move(config));
    }
    gfw->CreateRegions();
    flowContainer.setObject(new FlowContainer("FlowContainer"));
    const AxisSpec activityAxis = [this]() {
      if (cfgActivityAxis.value == NGlobalTracksAxis) {
        return AxisSpec{axisNGlobalTracks, "N_{global tracks}"};
      }
      if (cfgActivityAxis.value == ZdcSpectatorAxis) {
        return AxisSpec{axisSpectatorActivity, "N_{spectators}"};
      }
      return AxisSpec{axisCentrality, "centrality (%)"};
    }();
    if (activityAxis.nBins.has_value()) {
      flowContainer->Initialize(profileNames.get(), activityAxis.nBins.value(), activityAxis.binEdges.front(), activityAxis.binEdges.back(), cfgBootstrapSamples);
    } else {
      flowContainer->Initialize(profileNames.get(), activityAxis, cfgBootstrapSamples);
    }

    const AxisSpec etaAxis{32, -cuts.etaMax.value, cuts.etaMax.value, "#eta"};
    registry.add("event/centrality", "Accepted events;centrality;events", HistType::kTH1F, {axisCentrality});
    registry.add("event/activity", "Accepted events;event activity;events", HistType::kTH1F, {activityAxis});
    registry.add("event/vertex", "Accepted events;z-vtx (cm);events", HistType::kTH1F, {axisVertex});
    registry.add("event/aodSelection", "AOD event selection;selection;events", HistType::kTH1F, {{5, -0.5, 4.5}});
    registry.add("event/totalSpectators", "ZDC spectator estimate;N_{spectators};events", HistType::kTH1F, {axisTotalSpectators});
    registry.add("track/ptEta", "Selected tracks;p_{T} (GeV/c);#eta", HistType::kTH2F, {axisPt, etaAxis});
    registry.add("track/energy", "Energy entering the GFW;E_{T} (GeV);tracks", HistType::kTH1F, {axisEnergy});
    registry.add("track/pid", "PID mass assignment;species;tracks", HistType::kTH1F, {{4, -0.5, 3.5}});

    registry.add("zdc/qaSelection", "ZDC QA availability;condition;events", HistType::kTH1F, {{6, -0.5, 5.5}});
    registry.add("zdc/energyZNA", "Triggered events with matched ZDC;E_{ZNA};events", HistType::kTH1F, {axisZdcNeutronEnergy});
    registry.add("zdc/energyZNC", "Triggered events with matched ZDC;E_{ZNC};events", HistType::kTH1F, {axisZdcNeutronEnergy});
    registry.add("zdc/energyZPA", "Triggered events with matched ZDC;E_{ZPA};events", HistType::kTH1F, {axisZdcProtonEnergy});
    registry.add("zdc/energyZPC", "Triggered events with matched ZDC;E_{ZPC};events", HistType::kTH1F, {axisZdcProtonEnergy});
    registry.add("zdc/energyZEM", "Triggered events with matched ZDC;E_{ZEM1}+E_{ZEM2};events", HistType::kTH1F, {axisZdcZemEnergy});
    registry.add("zdc/timeZNA", "Triggered events with matched ZDC;t_{ZNA} (ns);events", HistType::kTH1F, {axisZdcTime});
    registry.add("zdc/timeZNC", "Triggered events with matched ZDC;t_{ZNC} (ns);events", HistType::kTH1F, {axisZdcTime});
    registry.add("zdc/timeZPA", "Triggered events with matched ZDC;t_{ZPA} (ns);events", HistType::kTH1F, {axisZdcTime});
    registry.add("zdc/timeZPC", "Triggered events with matched ZDC;t_{ZPC} (ns);events", HistType::kTH1F, {axisZdcTime});
    registry.add("zdc/timeZnSumVsDiff", "Triggered events with ZNA and ZNC timing;t_{ZNA}-t_{ZNC} (ns);t_{ZNA}+t_{ZNC} (ns)", HistType::kTH2F, {axisZdcTimeCombination, axisZdcTimeCombination});
    registry.add("zdc/energyZnaVsZnc", "Triggered events with finite ZN energies;E_{ZNC};E_{ZNA}", HistType::kTH2F, {axisZdcNeutronEnergy, axisZdcNeutronEnergy});
    registry.add("zdc/energyZpaVsZpc", "Triggered events with finite ZP energies;E_{ZPC};E_{ZPA}", HistType::kTH2F, {axisZdcProtonEnergy, axisZdcProtonEnergy});
    registry.add("zdc/energyZnVsZem", "Triggered events with finite ZN and ZEM energies;E_{ZEM1}+E_{ZEM2};E_{ZNA}+E_{ZNC}", HistType::kTH2F, {axisZdcZemEnergy, axisZdcTotalEnergy});
    registry.add("zdc/energyZnaVsTime", "Triggered events with matched ZDC;t_{ZNA} (ns);E_{ZNA}", HistType::kTH2F, {axisZdcTime, axisZdcNeutronEnergy});
    registry.add("zdc/energyZncVsTime", "Triggered events with matched ZDC;t_{ZNC} (ns);E_{ZNC}", HistType::kTH2F, {axisZdcTime, axisZdcNeutronEnergy});
    registry.add("zdc/energyZpaVsTime", "Triggered events with matched ZDC;t_{ZPA} (ns);E_{ZPA}", HistType::kTH2F, {axisZdcTime, axisZdcProtonEnergy});
    registry.add("zdc/energyZpcVsTime", "Triggered events with matched ZDC;t_{ZPC} (ns);E_{ZPC}", HistType::kTH2F, {axisZdcTime, axisZdcProtonEnergy});
    registry.add("zdc/commonVsSectorZNA", "Triggered events with finite ZNA energies;E_{ZNA}^{common};#Sigma E_{ZNA}^{sectors}", HistType::kTH2F, {axisZdcNeutronEnergy, axisZdcNeutronEnergy});
    registry.add("zdc/commonVsSectorZNC", "Triggered events with finite ZNC energies;E_{ZNC}^{common};#Sigma E_{ZNC}^{sectors}", HistType::kTH2F, {axisZdcNeutronEnergy, axisZdcNeutronEnergy});
    registry.add("zdc/commonVsSectorZPA", "Triggered events with finite ZPA energies;E_{ZPA}^{common};#Sigma E_{ZPA}^{sectors}", HistType::kTH2F, {axisZdcProtonEnergy, axisZdcProtonEnergy});
    registry.add("zdc/commonVsSectorZPC", "Triggered events with finite ZPC energies;E_{ZPC}^{common};#Sigma E_{ZPC}^{sectors}", HistType::kTH2F, {axisZdcProtonEnergy, axisZdcProtonEnergy});
    registry.add("zdc/energyZnVsNGlobalTracks", "Triggered events with finite ZN energies;N_{global tracks};E_{ZNA}+E_{ZNC}", HistType::kTH2F, {axisZdcTrackMultiplicity, axisZdcTotalEnergy});
    registry.add("zdc/energyZpVsNGlobalTracks", "Triggered events with finite ZP energies;N_{global tracks};E_{ZPA}+E_{ZPC}", HistType::kTH2F, {axisZdcTrackMultiplicity, axisZdcTotalEnergy});
    registry.add("zdc/energyTotalVsNGlobalTracks", "Triggered events with finite common energies;N_{global tracks};#Sigma E_{ZDC}^{common}", HistType::kTH2F, {axisZdcTrackMultiplicity, axisZdcTotalEnergy});
    registry.add("zdc/spectatorEstimateVsNGlobalTracks", "Triggered events passing current nonnegative-energy requirement;N_{global tracks};current N_{spectators} estimate", HistType::kTH2F, {axisZdcTrackMultiplicity, axisTotalSpectators});
    registry.add("zdc/energyZnVsCentrality", "Triggered events with finite ZN energies;centrality (%);E_{ZNA}+E_{ZNC}", HistType::kTH2F, {axisCentrality, axisZdcTotalEnergy});
    registry.add("zdc/energyZpVsCentrality", "Triggered events with finite ZP energies;centrality (%);E_{ZPA}+E_{ZPC}", HistType::kTH2F, {axisCentrality, axisZdcTotalEnergy});
    registry.add("zdc/energyTotalVsCentrality", "Triggered events with finite common energies;centrality (%);#Sigma E_{ZDC}^{common}", HistType::kTH2F, {axisCentrality, axisZdcTotalEnergy});
    registry.add("zdc/spectatorEstimateVsCentrality", "Triggered events passing current nonnegative-energy requirement;centrality (%);current N_{spectators} estimate", HistType::kTH2F, {axisCentrality, axisTotalSpectators});
    registry.add("zdc/meanZNAByNGlobalTracks", "Mean ZNA response;N_{global tracks};#LT E_{ZNA} #GT", HistType::kTProfile, {axisZdcTrackMultiplicity});
    registry.add("zdc/meanZNCByNGlobalTracks", "Mean ZNC response;N_{global tracks};#LT E_{ZNC} #GT", HistType::kTProfile, {axisZdcTrackMultiplicity});
    registry.add("zdc/meanZPAByNGlobalTracks", "Mean ZPA response;N_{global tracks};#LT E_{ZPA} #GT", HistType::kTProfile, {axisZdcTrackMultiplicity});
    registry.add("zdc/meanZPCByNGlobalTracks", "Mean ZPC response;N_{global tracks};#LT E_{ZPC} #GT", HistType::kTProfile, {axisZdcTrackMultiplicity});
    registry.add("zdc/meanZnByNGlobalTracks", "Mean neutron-ZDC response;N_{global tracks};#LT E_{ZNA}+E_{ZNC} #GT", HistType::kTProfile, {axisZdcTrackMultiplicity});
    registry.add("zdc/meanZpByNGlobalTracks", "Mean proton-ZDC response;N_{global tracks};#LT E_{ZPA}+E_{ZPC} #GT", HistType::kTProfile, {axisZdcTrackMultiplicity});
    registry.add("zdc/meanTotalByNGlobalTracks", "Mean total ZDC response;N_{global tracks};#LT #Sigma E_{ZDC}^{common} #GT", HistType::kTProfile, {axisZdcTrackMultiplicity});
    registry.add("zdc/meanZNAByCentrality", "Mean ZNA response;centrality (%);#LT E_{ZNA} #GT", HistType::kTProfile, {axisCentrality});
    registry.add("zdc/meanZNCByCentrality", "Mean ZNC response;centrality (%);#LT E_{ZNC} #GT", HistType::kTProfile, {axisCentrality});
    registry.add("zdc/meanZPAByCentrality", "Mean ZPA response;centrality (%);#LT E_{ZPA} #GT", HistType::kTProfile, {axisCentrality});
    registry.add("zdc/meanZPCByCentrality", "Mean ZPC response;centrality (%);#LT E_{ZPC} #GT", HistType::kTProfile, {axisCentrality});
    registry.add("zdc/meanZnByCentrality", "Mean neutron-ZDC response;centrality (%);#LT E_{ZNA}+E_{ZNC} #GT", HistType::kTProfile, {axisCentrality});
    registry.add("zdc/meanZpByCentrality", "Mean proton-ZDC response;centrality (%);#LT E_{ZPA}+E_{ZPC} #GT", HistType::kTProfile, {axisCentrality});
    registry.add("zdc/meanTotalByCentrality", "Mean total ZDC response;centrality (%);#LT #Sigma E_{ZDC}^{common} #GT", HistType::kTProfile, {axisCentrality});

    auto selectionHistogram = registry.get<TH1>(HIST("event/aodSelection"));
    selectionHistogram->GetXaxis()->SetBinLabel(AllAodEvents + 1, "all");
    selectionHistogram->GetXaxis()->SetBinLabel(TriggerSelected + 1, "event selection");
    selectionHistogram->GetXaxis()->SetBinLabel(CentralitySelected + 1, "centrality threshold");
    selectionHistogram->GetXaxis()->SetBinLabel(HasValidZdc + 1, "valid ZDC");
    selectionHistogram->GetXaxis()->SetBinLabel(SpectatorSelected + 1, "spectator threshold");
    auto zdcQaHistogram = registry.get<TH1>(HIST("zdc/qaSelection"));
    zdcQaHistogram->GetXaxis()->SetBinLabel(MatchedZdc + 1, "matched ZDC");
    zdcQaHistogram->GetXaxis()->SetBinLabel(FiniteCommonEnergy + 1, "finite common energies");
    zdcQaHistogram->GetXaxis()->SetBinLabel(NonnegativeCommonEnergy + 1, "nonnegative common energies");
    zdcQaHistogram->GetXaxis()->SetBinLabel(HasZnTiming + 1, "ZNA and ZNC timing");
    zdcQaHistogram->GetXaxis()->SetBinLabel(HasZpTiming + 1, "ZPA and ZPC timing");
    zdcQaHistogram->GetXaxis()->SetBinLabel(HasAllTiming + 1, "all timing");
    auto pidHistogram = registry.get<TH1>(HIST("track/pid"));
    pidHistogram->GetXaxis()->SetBinLabel(Unidentified + 1, "unidentified");
    pidHistogram->GetXaxis()->SetBinLabel(Pion + 1, "pion");
    pidHistogram->GetXaxis()->SetBinLabel(Kaon + 1, "kaon");
    pidHistogram->GetXaxis()->SetBinLabel(Proton + 1, "proton");

    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    const auto now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    ccdb->setCreatedNotAfter(now);
  }

  template <typename TZdc>
  float totalSpectators(TZdc const& zdc) const
  {
    const std::array<float, 4> energies{zdc.energyCommonZNA(), zdc.energyCommonZNC(), zdc.energyCommonZPA(), zdc.energyCommonZPC()};
    const std::array<float, 4> responseFactors{cfgZdcResponseZNA.value, cfgZdcResponseZNC.value, cfgZdcResponseZPA.value, cfgZdcResponseZPC.value};
    float totalEnergy = 0.f;
    for (int index = 0; index < 4; ++index) {
      const float value = energies[index];
      if (!std::isfinite(value) || value < 0.f) {
        return -1.f;
      }
      totalEnergy += value / responseFactors[index];
    }
    return totalEnergy / cfgSpectatorEnergyPerNucleon.value;
  }

  template <typename TZdc>
  void fillZdcQa(TZdc const& zdc, float numberOfGlobalTracks, float centrality)
  {
    registry.fill(HIST("zdc/qaSelection"), MatchedZdc);

    const std::array<float, 4> commonEnergies{zdc.energyCommonZNA(), zdc.energyCommonZNC(), zdc.energyCommonZPA(), zdc.energyCommonZPC()};
    const std::array<float, 4> times{zdc.timeZNA(), zdc.timeZNC(), zdc.timeZPA(), zdc.timeZPC()};
    const bool validCentrality = std::isfinite(centrality) && centrality >= 0.f && centrality <= 100.f;
    const bool allFinite = std::all_of(commonEnergies.begin(), commonEnergies.end(), [](float value) { return std::isfinite(value); });
    const bool allNonnegative = allFinite && std::all_of(commonEnergies.begin(), commonEnergies.end(), [](float value) { return value >= 0.f; });
    if (allFinite) {
      registry.fill(HIST("zdc/qaSelection"), FiniteCommonEnergy);
    }
    if (allNonnegative) {
      registry.fill(HIST("zdc/qaSelection"), NonnegativeCommonEnergy);
    }

    const auto hasTiming = [this](float value) { return std::isfinite(value) && std::abs(value) < cfgZdcMaxAbsTimeForQa.value; };
    const bool hasZnTiming = hasTiming(times[0]) && hasTiming(times[1]);
    const bool hasZpTiming = hasTiming(times[2]) && hasTiming(times[3]);
    if (hasZnTiming) {
      registry.fill(HIST("zdc/qaSelection"), HasZnTiming);
      registry.fill(HIST("zdc/timeZnSumVsDiff"), times[0] - times[1], times[0] + times[1]);
    }
    if (hasZpTiming) {
      registry.fill(HIST("zdc/qaSelection"), HasZpTiming);
    }
    if (hasZnTiming && hasZpTiming) {
      registry.fill(HIST("zdc/qaSelection"), HasAllTiming);
    }

    if (std::isfinite(commonEnergies[0])) {
      registry.fill(HIST("zdc/energyZNA"), commonEnergies[0]);
      registry.fill(HIST("zdc/meanZNAByNGlobalTracks"), numberOfGlobalTracks, commonEnergies[0]);
      if (validCentrality) {
        registry.fill(HIST("zdc/meanZNAByCentrality"), centrality, commonEnergies[0]);
      }
    }
    if (std::isfinite(commonEnergies[1])) {
      registry.fill(HIST("zdc/energyZNC"), commonEnergies[1]);
      registry.fill(HIST("zdc/meanZNCByNGlobalTracks"), numberOfGlobalTracks, commonEnergies[1]);
      if (validCentrality) {
        registry.fill(HIST("zdc/meanZNCByCentrality"), centrality, commonEnergies[1]);
      }
    }
    if (std::isfinite(commonEnergies[2])) {
      registry.fill(HIST("zdc/energyZPA"), commonEnergies[2]);
      registry.fill(HIST("zdc/meanZPAByNGlobalTracks"), numberOfGlobalTracks, commonEnergies[2]);
      if (validCentrality) {
        registry.fill(HIST("zdc/meanZPAByCentrality"), centrality, commonEnergies[2]);
      }
    }
    if (std::isfinite(commonEnergies[3])) {
      registry.fill(HIST("zdc/energyZPC"), commonEnergies[3]);
      registry.fill(HIST("zdc/meanZPCByNGlobalTracks"), numberOfGlobalTracks, commonEnergies[3]);
      if (validCentrality) {
        registry.fill(HIST("zdc/meanZPCByCentrality"), centrality, commonEnergies[3]);
      }
    }

    if (hasTiming(times[0])) {
      registry.fill(HIST("zdc/timeZNA"), times[0]);
    }
    if (hasTiming(times[1])) {
      registry.fill(HIST("zdc/timeZNC"), times[1]);
    }
    if (hasTiming(times[2])) {
      registry.fill(HIST("zdc/timeZPA"), times[2]);
    }
    if (hasTiming(times[3])) {
      registry.fill(HIST("zdc/timeZPC"), times[3]);
    }
    if (hasTiming(times[0]) && std::isfinite(commonEnergies[0])) {
      registry.fill(HIST("zdc/energyZnaVsTime"), times[0], commonEnergies[0]);
    }
    if (hasTiming(times[1]) && std::isfinite(commonEnergies[1])) {
      registry.fill(HIST("zdc/energyZncVsTime"), times[1], commonEnergies[1]);
    }
    if (hasTiming(times[2]) && std::isfinite(commonEnergies[2])) {
      registry.fill(HIST("zdc/energyZpaVsTime"), times[2], commonEnergies[2]);
    }
    if (hasTiming(times[3]) && std::isfinite(commonEnergies[3])) {
      registry.fill(HIST("zdc/energyZpcVsTime"), times[3], commonEnergies[3]);
    }

    const auto sumSectors = [](const auto& sectors) {
      float sum = 0.f;
      for (const float value : sectors) {
        if (!std::isfinite(value)) {
          return std::numeric_limits<float>::quiet_NaN();
        }
        sum += value;
      }
      return sum;
    };
    const std::array<float, 4> sectorSums{sumSectors(zdc.energySectorZNA()), sumSectors(zdc.energySectorZNC()), sumSectors(zdc.energySectorZPA()), sumSectors(zdc.energySectorZPC())};
    if (std::isfinite(commonEnergies[0]) && std::isfinite(sectorSums[0])) {
      registry.fill(HIST("zdc/commonVsSectorZNA"), commonEnergies[0], sectorSums[0]);
    }
    if (std::isfinite(commonEnergies[1]) && std::isfinite(sectorSums[1])) {
      registry.fill(HIST("zdc/commonVsSectorZNC"), commonEnergies[1], sectorSums[1]);
    }
    if (std::isfinite(commonEnergies[2]) && std::isfinite(sectorSums[2])) {
      registry.fill(HIST("zdc/commonVsSectorZPA"), commonEnergies[2], sectorSums[2]);
    }
    if (std::isfinite(commonEnergies[3]) && std::isfinite(sectorSums[3])) {
      registry.fill(HIST("zdc/commonVsSectorZPC"), commonEnergies[3], sectorSums[3]);
    }

    const float zem1 = zdc.energyZEM1();
    const float zem2 = zdc.energyZEM2();
    const bool finiteZem = std::isfinite(zem1) && std::isfinite(zem2);
    const float zemSum = finiteZem ? zem1 + zem2 : std::numeric_limits<float>::quiet_NaN();
    if (finiteZem) {
      registry.fill(HIST("zdc/energyZEM"), zemSum);
    }
    if (std::isfinite(commonEnergies[0]) && std::isfinite(commonEnergies[1])) {
      const float znSum = commonEnergies[0] + commonEnergies[1];
      registry.fill(HIST("zdc/energyZnaVsZnc"), commonEnergies[1], commonEnergies[0]);
      registry.fill(HIST("zdc/energyZnVsNGlobalTracks"), numberOfGlobalTracks, znSum);
      registry.fill(HIST("zdc/meanZnByNGlobalTracks"), numberOfGlobalTracks, znSum);
      if (validCentrality) {
        registry.fill(HIST("zdc/energyZnVsCentrality"), centrality, znSum);
        registry.fill(HIST("zdc/meanZnByCentrality"), centrality, znSum);
      }
      if (finiteZem) {
        registry.fill(HIST("zdc/energyZnVsZem"), zemSum, znSum);
      }
    }
    if (std::isfinite(commonEnergies[2]) && std::isfinite(commonEnergies[3])) {
      const float zpSum = commonEnergies[2] + commonEnergies[3];
      registry.fill(HIST("zdc/energyZpaVsZpc"), commonEnergies[3], commonEnergies[2]);
      registry.fill(HIST("zdc/energyZpVsNGlobalTracks"), numberOfGlobalTracks, zpSum);
      registry.fill(HIST("zdc/meanZpByNGlobalTracks"), numberOfGlobalTracks, zpSum);
      if (validCentrality) {
        registry.fill(HIST("zdc/energyZpVsCentrality"), centrality, zpSum);
        registry.fill(HIST("zdc/meanZpByCentrality"), centrality, zpSum);
      }
    }
    if (allFinite) {
      const float totalEnergy = commonEnergies[0] + commonEnergies[1] + commonEnergies[2] + commonEnergies[3];
      registry.fill(HIST("zdc/energyTotalVsNGlobalTracks"), numberOfGlobalTracks, totalEnergy);
      registry.fill(HIST("zdc/meanTotalByNGlobalTracks"), numberOfGlobalTracks, totalEnergy);
      if (validCentrality) {
        registry.fill(HIST("zdc/energyTotalVsCentrality"), centrality, totalEnergy);
        registry.fill(HIST("zdc/meanTotalByCentrality"), centrality, totalEnergy);
      }
      if (allNonnegative) {
        const float spectators = totalSpectators(zdc);
        registry.fill(HIST("zdc/spectatorEstimateVsNGlobalTracks"), numberOfGlobalTracks, spectators);
        if (validCentrality) {
          registry.fill(HIST("zdc/spectatorEstimateVsCentrality"), centrality, spectators);
        }
      }
    }
  }

  template <typename TBCs, typename TCollision>
  bool eventSelected(TCollision const& collision, bool triggerSelected, float numberOfGlobalTracks, float centrality, float& spectators)
  {
    spectators = -1.f;
    registry.fill(HIST("event/aodSelection"), AllAodEvents);
    if (!triggerSelected) {
      return false;
    }
    registry.fill(HIST("event/aodSelection"), TriggerSelected);

    if (!collision.has_foundBC()) {
      return false;
    }
    const auto& bc = collision.template foundBC_as<TBCs>();
    if (!bc.has_zdc()) {
      return false;
    }
    const auto& zdc = bc.zdc();
    fillZdcQa(zdc, numberOfGlobalTracks, centrality);

    if (!std::isfinite(centrality) || centrality < 0.f || centrality > cfgMaxCentrality) {
      return false;
    }
    registry.fill(HIST("event/aodSelection"), CentralitySelected);

    spectators = totalSpectators(zdc);
    if (spectators < 0.f) {
      return false;
    }
    registry.fill(HIST("event/aodSelection"), HasValidZdc);
    registry.fill(HIST("event/totalSpectators"), spectators);
    if (spectators > cfgMaxTotalSpectators) {
      return false;
    }
    registry.fill(HIST("event/aodSelection"), SpectatorSelected);
    return true;
  }

  void loadEfficiency(uint64_t timestamp)
  {
    if (efficiencyPath.value.empty()) {
      efficiency = nullptr;
      return;
    }
    if (efficiencyFromLocalFile) {
      if (localEfficiency) {
        efficiency = localEfficiency.get();
        return;
      }
      std::unique_ptr<TFile> input(TFile::Open(efficiencyPath.value.c_str(), "READ"));
      if (!input || input->IsZombie()) {
        LOGF(fatal, "Could not open efficiency file %s", efficiencyPath.value.c_str());
        return;
      }
      auto* source = dynamic_cast<THn*>(input->Get("ccdb_object"));
      if (!source || source->GetNdimensions() != 4) {
        LOGF(fatal, "Efficiency object must be a 4D THn with axes eta, pT, centrality, z-vtx");
      }
      localEfficiency.reset(dynamic_cast<THn*>(source->Clone("eecEfficiency")));
      if (!localEfficiency) {
        LOGF(fatal, "Could not clone the efficiency object from %s", efficiencyPath.value.c_str());
      }
      efficiency = localEfficiency.get();
      return;
    }
    if (!efficiency || !ccdb->isCachedObjectValid(efficiencyPath.value, timestamp)) {
      efficiency = ccdb->getForTimeStamp<THnT<float>>(efficiencyPath.value, timestamp);
      if (!efficiency || efficiency->GetNdimensions() != 4) {
        LOGF(fatal, "Could not load a 4D efficiency correction from %s", efficiencyPath.value.c_str());
      }
    }
  }

  void loadAcceptance(uint64_t timestamp, int runNumber)
  {
    if (acceptancePath.value.empty()) {
      acceptance = nullptr;
      return;
    }
    if (acceptance && (!acceptanceRunByRun || acceptanceRunNumber == runNumber)) {
      return;
    }
    std::string path = acceptancePath.value;
    if (acceptanceRunByRun) {
      if (path.back() != '/') {
        path += '/';
      }
      path += "RunByRun/";
    }
    acceptance = ccdb->getForTimeStamp<GFWWeights>(path, timestamp);
    if (!acceptance || !acceptance->isDataFilled()) {
      LOGF(fatal, "Could not load populated GFW acceptance weights from %s", path.c_str());
    }
    acceptanceRunNumber = runNumber;
  }

  double efficiencyWeight(float eta, float pt, float centrality, float posZ) const
  {
    if (!efficiency) {
      return 1.;
    }
    std::array<int, 4> bins{
      efficiency->GetAxis(0)->FindBin(eta),
      efficiency->GetAxis(1)->FindBin(pt),
      efficiency->GetAxis(2)->FindBin(centrality),
      efficiency->GetAxis(3)->FindBin(posZ)};
    return efficiency->GetBinContent(bins.data());
  }

  template <typename TTrack>
  int particleId(TTrack const& track) const
  {
    const std::array<float, 3> tpc{track.tpcNSigmaPi(), track.tpcNSigmaKa(), track.tpcNSigmaPr()};
    std::array<float, 3> selected = tpc;
    if (track.pt() > energy.tofPtCut.value) {
      if (!track.hasTOF()) {
        return Unidentified;
      }
      selected = {std::hypot(tpc[0], track.tofNSigmaPi()),
                  std::hypot(tpc[1], track.tofNSigmaKa()),
                  std::hypot(tpc[2], track.tofNSigmaPr())};
    }

    int bestSpecies = Unidentified;
    int compatibleSpecies = 0;
    float bestNSigma = energy.nSigmaCut.value;
    for (int species = 0; species < 3; ++species) {
      const float absoluteNSigma = std::abs(selected[species]);
      if (absoluteNSigma >= energy.nSigmaCut.value) {
        continue;
      }
      ++compatibleSpecies;
      if (absoluteNSigma < bestNSigma) {
        bestNSigma = absoluteNSigma;
        bestSpecies = species + 1;
      }
    }
    if (energy.strictPid.value && compatibleSpecies > 1) {
      return Unidentified;
    }
    return bestSpecies;
  }

  double massForParticle(int pid) const
  {
    switch (pid) {
      case Pion:
        return o2::constants::physics::MassPionCharged;
      case Kaon:
        return o2::constants::physics::MassKPlus;
      case Proton:
        return o2::constants::physics::MassProton;
      default:
        return energy.unidentifiedMass.value;
    }
  }

  template <typename TTrack>
  double transverseEnergy(TTrack const& track, int pid) const
  {
    if (energy.weightMode.value == PtWeight) {
      return track.pt();
    }
    const double mass = massForParticle(pid);
    if (energy.weightMode.value == MtWeight) {
      return std::hypot(track.pt(), mass);
    }
    return std::hypot(track.pt(), mass / std::cosh(track.eta()));
  }

  template <typename TCollision, typename TTracks>
  void processCollision(TCollision const& collision, TTracks const& tracks, float centrality, float spectators)
  {
    float activity = centrality;
    if (cfgActivityAxis.value == NGlobalTracksAxis) {
      activity = static_cast<float>(tracks.size());
    } else if (cfgActivityAxis.value == ZdcSpectatorAxis) {
      activity = spectators;
    }
    if (!std::isfinite(activity)) {
      return;
    }

    gfw->Clear();
    for (const auto& track : tracks) {
      const int pid = energy.weightMode.value == PtWeight ? Unidentified : particleId(track);
      const double et = transverseEnergy(track, pid);
      const double weff = efficiencyWeight(track.eta(), track.pt(), centrality, collision.posZ());
      const double wacc = acceptance ? acceptance->getNUA(track.phi(), track.eta(), collision.posZ()) : 1.;
      if (!std::isfinite(et) || !std::isfinite(weff) || !std::isfinite(wacc) || et <= 0. || weff <= 0. || wacc <= 0.) {
        continue;
      }
      gfw->Fill(track.eta(), 0, track.phi(), et * weff * wacc, gfwMask);
      registry.fill(HIST("track/ptEta"), track.pt(), track.eta());
      registry.fill(HIST("track/energy"), et);
      if (energy.weightMode.value != PtWeight) {
        registry.fill(HIST("track/pid"), pid);
      }
    }

    registry.fill(HIST("event/centrality"), centrality);
    registry.fill(HIST("event/activity"), activity);
    registry.fill(HIST("event/vertex"), collision.posZ());
    fillFlowContainer(activity);
  }

  void fillFlowContainer(float activity)
  {
    const double randomNumber = random->Rndm();
    for (std::size_t index = 0; index < corrConfigs.size(); ++index) {
      const auto& config = corrConfigs[index];
      if (zeroHarmonicConfigs[index]) {
        const double raw = gfw->Calculate(config, 0, false).real();
        if (std::isfinite(raw) && raw > 0.) {
          const std::string rawName = config.Head + "_raw";
          flowContainer->FillProfile(rawName.c_str(), activity, raw, 1., randomNumber);
        }
        continue;
      }

      const double denominator = gfw->Calculate(config, 0, true).real();
      if (!std::isfinite(denominator) || denominator <= 0.) {
        continue;
      }
      const double raw = gfw->Calculate(config, 0, false).real();
      if (!std::isfinite(raw)) {
        continue;
      }
      const std::string rawName = config.Head + "_raw";
      flowContainer->FillProfile(rawName.c_str(), activity, raw, 1., randomNumber);
      flowContainer->FillProfile(config.Head.c_str(), activity, raw / denominator, denominator, randomNumber);
    }
  }

  void processRun2(CollisionsRun2::iterator const& collision, Tracks const& tracks, BCsRun2 const&, aod::Zdcs const&)
  {
    float spectators = -1.f;
    if (!eventSelected<BCsRun2>(collision, collision.alias_bit(kINT7) && collision.sel7(), static_cast<float>(tracks.size()), collision.centRun2V0M(), spectators)) {
      return;
    }
    const auto& bc = collision.foundBC_as<BCsRun2>();
    loadEfficiency(bc.timestamp());
    loadAcceptance(bc.timestamp(), bc.runNumber());
    processCollision(collision, tracks, collision.centRun2V0M(), spectators);
  }
  PROCESS_SWITCH(EnergyFlowGfw, processRun2, "Process Run 2 AOD data", false);

  void processRun3(CollisionsRun3::iterator const& collision, Tracks const& tracks, BCsRun3 const&, aod::Zdcs const&)
  {
    float spectators = -1.f;
    if (!eventSelected<BCsRun3>(collision, collision.sel8(), static_cast<float>(tracks.size()), collision.centFT0C(), spectators)) {
      return;
    }
    const auto& bc = collision.foundBC_as<BCsRun3>();
    loadEfficiency(bc.timestamp());
    loadAcceptance(bc.timestamp(), bc.runNumber());
    processCollision(collision, tracks, collision.centFT0C(), spectators);
  }
  PROCESS_SWITCH(EnergyFlowGfw, processRun3, "Process Run 3 AOD data", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<EnergyFlowGfw>(cfgc)};
}
