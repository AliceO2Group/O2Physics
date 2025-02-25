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

/// \file emcMCTmMonitor.cxx
/// \brief Simple monitoring task for EMCal clusters
/// \author Marvin Hemmer <marvin.hemmer@cern.ch>
/// \since 13.06.2025
///
/// This task is meant to be used for monitoring the matching between global tracks and EMCal clusters
/// properties, such as:
/// - cluster energy over track momentum
/// - difference in eta
/// - difference in phi
/// Simple event selection using the flag doEventSel is provided, which selects INT7 events if set to 1
/// For pilot beam data, instead of relying on the event selection, one can veto specific BC IDS using the flag
/// fDoVetoBCID.

#include <climits>
#include <cstdlib>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <algorithm>
#include <vector>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/HistogramRegistry.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "CommonConstants/MathConstants.h"

#include "EMCALBase/Geometry.h"
#include "EMCALCalib/BadChannelMap.h"
#include "PWGJE/DataModel/EMCALClusters.h"
#include "DataFormatsEMCAL/Cell.h"
#include "DataFormatsEMCAL/Constants.h"
#include "DataFormatsEMCAL/AnalysisCluster.h"

#include "CommonDataFormat/InteractionRecord.h"

using namespace o2::framework;
using namespace o2::framework::expressions;
using CollisionEvSel = o2::soa::Join<o2::aod::Collisions, o2::aod::EvSels, o2::aod::McCollisionLabel>;
using CollisionEvSelIt = o2::soa::Join<o2::aod::Collisions, o2::aod::EvSels, o2::aod::McCollisionLabel>::iterator;
using SelectedClusters = o2::soa::Filtered<o2::soa::Join<o2::aod::EMCALClusters, o2::aod::EMCALMCClusters>>;
using McTracks = o2::soa::Join<o2::aod::Tracks, o2::aod::TracksExtra, o2::aod::McTrackLabels>;

/// \brief Function to put labels on monitoring histogram
/// \param hRejection monitoring histogram
/// \param softwareTriggerLabel bin label for software trigger rejection
template <typename Histo>
void setParticleLabels(Histo& hParticle)
{
  auto axis = hParticle->GetAxis(0);
  axis->SetBinLabel(1, "e^{#pm}");
  axis->SetBinLabel(2, "#mu^{#pm}");
  axis->SetBinLabel(3, "#tau^{#pm}");
  axis->SetBinLabel(4, "#gamma");
  axis->SetBinLabel(5, "#pi^{0}");
  axis->SetBinLabel(6, "#pi^{#pm}");
  axis->SetBinLabel(7, "#eta");
  axis->SetBinLabel(8, "K");
  axis->SetBinLabel(9, "#eta\'");
  axis->SetBinLabel(10, "p");
  axis->SetBinLabel(11, "n");
  axis->SetBinLabel(12, "d");
  axis->SetBinLabel(13, "t");
  axis->SetBinLabel(14, "other");
}

template <typename Histo>
void setParticleLabels1D(Histo& hParticle)
{
  auto axis = hParticle->GetXaxis();
  axis->SetBinLabel(1, "e^{#pm}");
  axis->SetBinLabel(2, "#mu^{#pm}");
  axis->SetBinLabel(3, "#tau^{#pm}");
  axis->SetBinLabel(4, "#gamma");
  axis->SetBinLabel(5, "#pi^{0}");
  axis->SetBinLabel(6, "#pi^{#pm}");
  axis->SetBinLabel(7, "#eta");
  axis->SetBinLabel(8, "K");
  axis->SetBinLabel(9, "#eta\'");
  axis->SetBinLabel(10, "p");
  axis->SetBinLabel(11, "n");
  axis->SetBinLabel(12, "d");
  axis->SetBinLabel(13, "t");
  axis->SetBinLabel(14, "other");
}

bool isInEMCal(float eta, float phi)
{
  // Acceptable eta range for EMCal
  if (std::fabs(eta) >= 0.7f) {
    return false;
  }

  // Wrap phi into [0, 2π)
  float phiWrapped = RecoDecay::constrainAngle(phi);

  // Convert degrees to radians for cuts
  constexpr float phiMin = 80.f * o2::constants::math::Deg2Rad;
  constexpr float phiMax = 187.f * o2::constants::math::Deg2Rad;
  constexpr float dcalMin = 260.f * o2::constants::math::Deg2Rad;
  constexpr float dcalMax = 327.f * o2::constants::math::Deg2Rad;
  constexpr float transitionMin = 260.f * o2::constants::math::Deg2Rad;
  constexpr float transitionMax = 320.f * o2::constants::math::Deg2Rad;

  // Exclude phi outside EMCal sector (80° < φ < 187°)
  if (phiWrapped <= phiMin || phiWrapped >= dcalMax || (phiWrapped >= phiMax && phiWrapped <= dcalMin)) {
    return false;
  }

  // Exclude transition region with limited acceptance due to PHOS gap
  if ((phiWrapped > transitionMin && phiWrapped < transitionMax) && std::fabs(eta) <= 0.22f) {
    return false;
  }

  return true;
}

// calculate ratio of cluster energy to track momentum
double calcRcorr(double Eclus, double P)
{
  if (Eclus < P) {
    return Eclus / P;
  } else {
    // if we have more track momentum than energy in the cluster we "choose" to just output 1
    return 1.f;
  }
}

struct EmcMCTmMonitor {
  HistogramRegistry mHistManager{"TrackMatchingMonitorHistograms", {}, OutputObjHandlingPolicy::AnalysisObject};
  o2::emcal::Geometry* mGeometry = nullptr;

  Preslice<o2::aod::EMCALMatchedTracks> perClusterMatchedTracks = o2::aod::emcalclustercell::emcalclusterId;
  Preslice<o2::aod::McParticles> perCollision = o2::aod::mcparticle::mcCollisionId;
  // configurable parameters
  Configurable<bool> doEventSel{"doEventSel", 0, "demand kINT7"};
  Configurable<double> vertexCut{"vertexCut", -1, "apply z-vertex cut with value in cm"};
  Configurable<int> clusterDefinition{"clusterDefinition", 10, "cluster definition to be selected, e.g. 10=kV3Default"};
  ConfigurableAxis clustertimeBinning{"clustertimeBinning", {1500, -600, 900}, ""};
  Configurable<float> minTime{"minTime", -25., "Minimum cluster time for time cut"};
  Configurable<float> maxTime{"maxTime", +20., "Maximum cluster time for time cut"};
  Configurable<float> minDEta{"minDEta", 0.015, "Minimum dEta between track and cluster"};
  Configurable<float> minDPhi{"minDPhi", 0.03, "Minimum dPhi between track and cluster"};
  Configurable<std::vector<float>> eOverPRange{"eOverPRange", {0.9, 1.2}, "E/p range where one would search for electrons (first <= E/p <= second)"};
  Configurable<bool> useAlignmentFromCCDB{"useAlignmentFromCCDB", false, "States if alignment objects should be used from CCDB"};
  Configurable<bool> doFindFirst{"doFindFirst", false, "States if one wants to search and find the first appearance of a MCParticle that might have scattered multiple times. Heavy on CPU time!"};

  // vector to store important PDG Codes
  std::vector<int> mPDGCodes;

  /// \brief Create output histograms and initialize geometry
  void init(InitContext const&)
  {

    //           e,  mu, tau, y, pi0,  pi,  eta, K,  eta', p,    n,       d,          t
    mPDGCodes = {11, 13, 15, 22, 111, 211, 221, 321, 331, 2212, 2112, 1000010020, 1000010030};
    // create histograms
    using O2HistType = HistType;
    using O2Axis = AxisSpec;

    // load geometry just in case we need it
    const int runNumberForGeom = 300000;
    mGeometry = o2::emcal::Geometry::GetInstanceFromRunNumber(runNumberForGeom);
    if (useAlignmentFromCCDB.value) {
      mGeometry->SetMisalMatrixFromCcdb("Users/m/mhemmer/EMCAL/Config/GeometryAligned", 10000);
    }

    // create common axes
    LOG(info) << "Creating histograms";
    const O2Axis bcAxis{3501, -0.5, 3500.5};
    const O2Axis clusterEnergyAxis{makeEnergyBinningAliPhysics(), "#it{E}_{clus} (GeV)"};
    const O2Axis energyAxis{makeEnergyBinningAliPhysics(), "#it{E} (GeV)"};
    const O2Axis energyFractionAxis{110, 0., 1.1, "#it{E}/#it{E}_{clus}"};
    const O2Axis energyFractionParticleAxis{500, 0., 5., "#it{E}_{part}/(#it{E}_{clus}*frac)"};
    const O2Axis energyParticle{200, 0., 20., "#it{E}_{part} (GeV)"};
    const O2Axis energyCluster{200, 0., 20., "#it{E}_{clus} (GeV)"};
    const O2Axis energyParticleCluster{200, 0., 20., "#it{E}_{clus}*frac (GeV)"};
    const O2Axis productionRadius{1000, 0., 500., "#it{R} (cm)"};
    const O2Axis amplitudeAxisLarge{1000, 0., 100., "amplitudeLarge", "Amplitude (GeV)"};
    const O2Axis dEtaAxis{100, -1.f * minDEta, minDEta, "d#it{#eta}"};
    const O2Axis dPhiAxis{100, -1.f * minDPhi, minDPhi, "d#it{#varphi} (rad)"};
    const O2Axis dRAxis{150, 0.0, 0.015, "d#it{R}"};
    const O2Axis eoverpAxis{500, 0, 10, "#it{E}_{cluster}/#it{p}_{track}"};
    const O2Axis nSigmaAxis{400, -10., +30., "N#sigma"};
    const O2Axis trackptAxis{makePtBinning(), "#it{p}_{T,track}"};
    const O2Axis trackpAxis{200, 0, 100, "#it{p}_{track}"};
    const O2Axis clusterptAxis{makePtBinning(), "#it{p}_{T}"};
    const O2Axis etaAxis{160, -0.8, 0.8, "#eta"};
    const O2Axis phiAxis{72, 0, 2 * 3.14159, "#varphi (rad)"};
    const O2Axis smAxis{20, -0.5, 19.5, "SM"};
    const O2Axis ptParticle{200, 0., 20., "#it{p}_{T,track} (GeV)"};
    const O2Axis momentumParticle{200, 0., 20., "#it{p}_{part} (GeV/#it{c})"};
    const O2Axis particleSpecies{static_cast<int>(mPDGCodes.size()) + 1, -0.5, static_cast<float>(mPDGCodes.size()) + 0.5, "species"};
    const O2Axis rCorrAxis{100, 0.01f, 1.01f, "#it{R}_{corr}"};
    const O2Axis matchAxis{2, -0.5, 1.5, ""};
    const O2Axis mcProcessAxis{50, -0.5, 49.5, "MC process"};
    O2Axis timeAxis{clustertimeBinning, "t_{cl} (ns)"};

    int maxMatched = 20; // maximum number of matched tracks, hardcoded in emcalCorrectionTask.cxx!
    const O2Axis nmatchedtrack{maxMatched, -0.5, maxMatched + 0.5};

    // event properties
    mHistManager.add("eventsAll", "Number of events", O2HistType::kTH1F, {{1, 0.5, 1.5}});
    mHistManager.add("eventsSelected", "Number of events", O2HistType::kTH1F, {{1, 0.5, 1.5}});
    mHistManager.add("eventVertexZAll", "z-vertex of event (all events)", O2HistType::kTH1F, {{200, -20, 20}});
    mHistManager.add("eventVertexZSelected", "z-vertex of event (selected events)", O2HistType::kTH1F, {{200, -20, 20}});

    mHistManager.add("hContributors", "hContributors;contributor per cluster;#it{counts}", O2HistType::kTH1I, {{100, 0, 100}});

    mHistManager.add("hClusterEnergyParticleEnergy", "hClusterEnergyParticleEnergy", O2HistType::kTH2D, {{energyParticle, energyParticleCluster}});

    mHistManager.add("hRCorrNeutronTrue", "hRCorrNeutronTrue", O2HistType::kTH2D, {{rCorrAxis, momentumParticle}});
    mHistManager.add("hRCorrChargedReco", "hRCorrChargedReco", O2HistType::kTH2D, {{rCorrAxis, momentumParticle}});
    mHistManager.add("hRCorrChargedTrue", "hRCorrChargedTrue", O2HistType::kTH2D, {{rCorrAxis, momentumParticle}});

    auto hMatchType = mHistManager.add<TH2>("hMatchType", ";;counts", kTH2D, {{matchAxis, ptParticle}}, false);
    hMatchType->GetXaxis()->SetBinLabel(1, "match");
    hMatchType->GetXaxis()->SetBinLabel(2, "Miss match");

    // cluster properties (matched clusters)
    auto hFractionVsParticle = mHistManager.add<THnSparse>("hFractionVsParticle", ";;#it{E}/#it{E}_{clus}", O2HistType::kTHnSparseF, {{particleSpecies, energyFractionAxis, energyParticle, energyCluster}}, false);
    auto hParticleFractionVsParticle = mHistManager.add<THnSparse>("hParticleFractionVsParticle", ";;#it{E}/#it{E}_{clus}", O2HistType::kTHnSparseF, {{particleSpecies, energyFractionParticleAxis, energyParticle, energyCluster}}, false);
    auto hRadiusVsParticle = mHistManager.add<THnSparse>("hRadiusVsParticle", ";;#it{E} (GeV)", O2HistType::kTHnSparseF, {{particleSpecies, productionRadius, energyParticle, energyCluster}}, false);
    // cluster properties (matched clusters)
    auto hFractionVsParticleSingleCell = mHistManager.add<THnSparse>("hFractionVsParticleSingleCell", ";;#it{E}/#it{E}_{clus}", O2HistType::kTHnSparseF, {{particleSpecies, energyFractionAxis, energyParticle, energyCluster}}, false);
    auto hParticleFractionVsParticleSingleCell = mHistManager.add<THnSparse>("hParticleFractionVsParticleSingleCell", ";;#it{E}/#it{E}_{clus}", O2HistType::kTHnSparseF, {{particleSpecies, energyFractionParticleAxis, energyParticle, energyCluster}}, false);
    auto hEnergyVsAllParticle = mHistManager.add<TH2>("hEnergyVsAllParticle", ";;#it{E} (GeV)", kTH2D, {{particleSpecies, energyAxis}}, false);
    auto hRadiusVsAllParticle = mHistManager.add<TH2>("hRadiusVsAllParticle", ";;#it{E} (GeV)", kTH2D, {{particleSpecies, productionRadius}}, false);
    auto hEnergyVsOnEMCParticle = mHistManager.add<TH2>("hEnergyVsOnEMCParticle", ";;#it{E} (GeV)", kTH2D, {{particleSpecies, energyAxis}}, false);
    auto hParticleVsProcess = mHistManager.add<TH3>("hParticleVsProcess", ";;;#it{E} (GeV)", kTH3D, {{particleSpecies, mcProcessAxis, energyAxis}}, false);

    setParticleLabels(hFractionVsParticle);
    setParticleLabels(hParticleFractionVsParticle);
    setParticleLabels(hRadiusVsParticle);
    setParticleLabels(hFractionVsParticleSingleCell);
    setParticleLabels(hParticleFractionVsParticleSingleCell);
    setParticleLabels1D(hEnergyVsAllParticle);
    setParticleLabels1D(hRadiusVsAllParticle);
    setParticleLabels1D(hEnergyVsOnEMCParticle);
    setParticleLabels1D(hParticleVsProcess);
  }

  // define cluster filter. It selects only those clusters which are of the type
  // sadly passing of the string at runtime is not possible for technical region so cluster definition is
  // an integer instead
  Filter clusterDefinitionSelection = (o2::aod::emcalcluster::definition == clusterDefinition) && (o2::aod::emcalcluster::time >= minTime) && (o2::aod::emcalcluster::time <= maxTime);

  /// \brief Process EMCAL clusters that are matched to a collisions
  void processCollisions(CollisionEvSelIt const& theCollision, SelectedClusters const& clusters, o2::aod::EMCALMatchedTracks const& matchedtracks, McTracks const& /* alltracks */, o2::aod::McParticles const& mcParticles)
  {
    mHistManager.fill(HIST("eventsAll"), 1);

    // do event selection if doEventSel is specified
    // currently the event selection is hard coded to kINT7
    // but other selections are possible that are defined in TriggerAliases.h
    bool isSelected = true;
    int minimumRun3Number = 300000;
    if (doEventSel) {
      if (theCollision.bc().runNumber() < minimumRun3Number) {
        if (!theCollision.alias_bit(kINT7)) {
          isSelected = false;
        }
      } else {
        if (!theCollision.alias_bit(kTVXinEMC)) {
          isSelected = false;
        }
      }
    }
    if (!isSelected) {
      LOG(debug) << "Event not selected because it is not kINT7 or does not have EMCAL in readout, skipping";
      return;
    }
    mHistManager.fill(HIST("eventVertexZAll"), theCollision.posZ());
    if (vertexCut > 0 && std::abs(theCollision.posZ()) > vertexCut) {
      LOG(debug) << "Event not selected because of z-vertex cut z= " << theCollision.posZ() << " > " << vertexCut << " cm, skipping";
      return;
    }
    mHistManager.fill(HIST("eventsSelected"), 1);
    mHistManager.fill(HIST("eventVertexZSelected"), theCollision.posZ());

    auto mcParticlesInColl = mcParticles.sliceBy(perCollision, theCollision.mcCollisionId());
    for (const auto& mcParticle : mcParticlesInColl) {
      int pdgCode = std::abs(mcParticle.pdgCode());
      auto pdgIter = std::find(mPDGCodes.begin(), mPDGCodes.end(), pdgCode);
      int pdgIndex = pdgIter - mPDGCodes.begin();
      mHistManager.fill(HIST("hEnergyVsAllParticle"), pdgIndex, mcParticle.e());
      mHistManager.fill(HIST("hRadiusVsAllParticle"), pdgIndex, std::sqrt(mcParticle.vx() * mcParticle.vx() + mcParticle.vy() * mcParticle.vy()));

      if (isInEMCal(mcParticle.eta(), mcParticle.phi())) {
        mHistManager.fill(HIST("hEnergyVsOnEMCParticle"), pdgIndex, mcParticle.e());
      }
    }

    // loop over all clusters from accepted collision
    for (const auto& cluster : clusters) {
      // TODO:
      // check anti neutron is a thing?
      // 2D Radius vs Z (p, n) with antis
      // Check MCType
      // Sandro: Haben wir Material Interaction Maps? (X0 Verteilung und lambda Verteilung),
      // Gibt es ein Macro, was uns die Materialverteilung in radiation und interaction length geben kann für bestimmte eta/phi bei bestimmten R
      // F+ Verteilung vs cluster E (Figure 17 EMCal paper)
      // Figure 16 und 15
      // match / miss match vs track pT
      // E_cluster * frac vs E_cluster für E_part < 500 MeV (e + gamma vs rest)
      // E_cluster * frac vs E_part für leading particle (vs frac)
      double dEta, dPhi, trackEta, trackPhi;
      // int supermoduleID;
      std::vector<int> clusterPDG;
      std::vector<int> clusterMCId;
      std::vector<int> clusterUniqueMCId;
      clusterPDG.reserve(cluster.mcParticleIds().size());
      clusterMCId.reserve(cluster.mcParticleIds().size());
      clusterUniqueMCId.reserve(cluster.mcParticleIds().size());
      if (cluster.mcParticleIds().size() != cluster.amplitudeA().size()) {
        LOG(info) << "Number of mc particles " << cluster.mcParticleIds().size() << " does not match number of amplitude fractions " << cluster.amplitudeA().size();
      }
      mHistManager.fill(HIST("hContributors"), cluster.mcParticleIds().size());

      float pNeutron = 0.f;
      float pCharged = 0.f;
      for (size_t iParticle = 0; iParticle < cluster.amplitudeA().size(); ++iParticle) {
        auto mcParticle = mcParticles.iteratorAt(cluster.mcParticleIds()[iParticle]);
        float fraction = cluster.amplitudeA()[iParticle];
        int pdgCode = std::abs(mcParticle.pdgCode());
        clusterPDG.emplace_back(pdgCode);
        clusterMCId.emplace_back(cluster.mcParticleIds()[iParticle]);

        if (pdgCode == 2112) {
          pNeutron += mcParticle.p();
        } else if (pdgCode == 11 || pdgCode == 211 || pdgCode == 321 || pdgCode == 2212 || pdgCode == 13 || pdgCode == 1000010020 || pdgCode == 1000010030) {
          pCharged += mcParticle.p();
        }

        // Trying to figure out if the particle was scattered, if so, get the original particle:
        auto current = mcParticle;
        while (doFindFirst) {
          auto mothers = current.mothers_as<o2::aod::McParticles>();
          if (mothers.size() != 1) {
            break;
          }

          auto& mother = *mothers.begin(); // only one, so safe
          if (std::abs(mother.pdgCode()) != pdgCode) {
            break;
          }
          current = mother; // continue walking up
        }
        // use the current particle to obtain the production radius instead of the scattered particle

        auto uniqueIdIter = std::find(clusterUniqueMCId.begin(), clusterUniqueMCId.end(), current.globalIndex());
        if (uniqueIdIter == clusterUniqueMCId.end()) {
          clusterUniqueMCId.push_back(current.globalIndex());
        } else {
          LOG(info) << "MCParticle " << mcParticle.globalIndex() << " is just a rescatter of particle " << current.globalIndex() << " which already has another hit in the same cluster!";
        }

        auto pdgIter = std::find(mPDGCodes.begin(), mPDGCodes.end(), pdgCode);
        int pdgIndex = pdgIter - mPDGCodes.begin();

        mHistManager.fill(HIST("hFractionVsParticle"), pdgIndex, fraction, mcParticle.e(), cluster.energy());
        mHistManager.fill(HIST("hParticleFractionVsParticle"), pdgIndex, mcParticle.e() / (cluster.energy() * fraction), mcParticle.e(), cluster.energy());

        mHistManager.fill(HIST("hClusterEnergyParticleEnergy"), mcParticle.e(), (cluster.energy() * fraction));
        mHistManager.fill(HIST("hRadiusVsParticle"), pdgIndex, std::sqrt(current.vx() * current.vx() + current.vy() * current.vy()), mcParticle.e(), cluster.energy());

        if (cluster.nCells() == 1) {
          mHistManager.fill(HIST("hFractionVsParticleSingleCell"), pdgIndex, fraction, mcParticle.e(), cluster.energy());
          mHistManager.fill(HIST("hParticleFractionVsParticleSingleCell"), pdgIndex, mcParticle.e() / (cluster.energy() * fraction), mcParticle.e(), cluster.energy());
        }

        mHistManager.fill(HIST("hParticleVsProcess"), pdgIndex, current.statusCode(), mcParticle.e());
      } // end of loop over mcParticles in cluster

      mHistManager.fill(HIST("hRCorrNeutronTrue"), calcRcorr(cluster.energy(), pNeutron), pNeutron);
      mHistManager.fill(HIST("hRCorrChargedTrue"), calcRcorr(cluster.energy(), pCharged), pCharged);

      auto tracksofcluster = matchedtracks.sliceBy(perClusterMatchedTracks, cluster.globalIndex());
      float pChargedReco = 0.f;
      for (const auto& match : tracksofcluster) {
        // only consider closest match
        trackEta = match.track_as<McTracks>().trackEtaEmcal();
        trackPhi = match.track_as<McTracks>().trackPhiEmcal();
        dPhi = trackPhi - cluster.phi();
        dEta = trackEta - cluster.eta();
        if (std::fabs(dEta) >= minDEta || std::fabs(dPhi) >= minDPhi) { // dEta and dPhi cut
          continue;
        }
        pChargedReco += match.track_as<McTracks>().p();
        // Get the MC Particle coresponding to the track that was just matched
        // auto trackMCParticle = mcParticles.iteratorAt(match.track_as<McTracks>().mcParticleId());
        // int trackPdgCode = std::abs(trackMCParticle.pdgCode());

        // try to find the particle in the cluster
        auto trackPdgIter = std::find(clusterMCId.begin(), clusterMCId.end(), match.track_as<McTracks>().mcParticleId());

        if (trackPdgIter != clusterMCId.end()) {
          // track deposit energy into cluster
          mHistManager.fill(HIST("hMatchType"), 0, match.track_as<McTracks>().pt());
        } else {
          // track did not deposit energy into cluster
          mHistManager.fill(HIST("hMatchType"), 1, match.track_as<McTracks>().pt());
        }
      } // end of loop pver tracks of clusters
      mHistManager.fill(HIST("hRCorrChargedReco"), calcRcorr(cluster.energy(), pChargedReco), pChargedReco);
    } // end of loop over cluster
  }
  PROCESS_SWITCH(EmcMCTmMonitor, processCollisions, "Process clusters from collision", true);

  /// \brief Create binning for cluster energy axis (variable bin size)
  /// \return vector with bin limits
  std::vector<double>
    makeEnergyBinning() const
  {
    auto fillBinLimits = [](std::vector<double>& binlimits, double max, double binwidth) {
      auto current = *binlimits.rbegin();
      while (current < max) {
        current += binwidth;
        binlimits.emplace_back(current);
      }
    };
    std::vector<double> result = {0.};
    fillBinLimits(result, 2., 0.1);
    fillBinLimits(result, 5., 0.2);
    fillBinLimits(result, 10., 0.5);
    fillBinLimits(result, 20., 1.);
    fillBinLimits(result, 50., 2.);
    fillBinLimits(result, 100., 5.);
    fillBinLimits(result, 200., 10.);
    return result;
  }

  /// \brief Create binning for cluster energy axis (variable bin size)
  /// direct port from binning often used in AliPhysics for debugging
  /// \return vector with bin limits
  std::vector<double> makeEnergyBinningAliPhysics() const
  {

    std::vector<double> result;
    int nBinsClusterE = 235;
    for (int i = 0; i < nBinsClusterE + 1; i++) {
      if (i < 1)
        result.emplace_back(0.3 * i);
      else if (i < 55)
        result.emplace_back(0.3 + 0.05 * (i - 1));
      else if (i < 105)
        result.emplace_back(3. + 0.1 * (i - 55));
      else if (i < 140)
        result.emplace_back(8. + 0.2 * (i - 105));
      else if (i < 170)
        result.emplace_back(15. + 0.5 * (i - 140));
      else if (i < 190)
        result.emplace_back(30. + 1.0 * (i - 170));
      else if (i < 215)
        result.emplace_back(50. + 2.0 * (i - 190));
      else if (i < 235)
        result.emplace_back(100. + 5.0 * (i - 215));
      else if (i < 245)
        result.emplace_back(200. + 10.0 * (i - 235));
    }
    return result;
  }

  /// \brief Create binning for pT axis (variable bin size)
  /// direct port from binning often used in AliPhysics for debugging
  /// \return vector with bin limits
  std::vector<double> makePtBinning() const
  {
    std::vector<double> result;
    result.reserve(1000);
    double epsilon = 1e-6;
    double valGammaPt = 0;
    for (int i = 0; i < 1000; ++i) {
      result.push_back(valGammaPt);
      if (valGammaPt < 1.0 - epsilon)
        valGammaPt += 0.1;
      else if (valGammaPt < 5 - epsilon)
        valGammaPt += 0.2;
      else if (valGammaPt < 10 - epsilon)
        valGammaPt += 0.5;
      else if (valGammaPt < 50 - epsilon)
        valGammaPt += 1;
      else if (valGammaPt < 100 - epsilon)
        valGammaPt += 5;
      else
        break;
    }
    return result;
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<EmcMCTmMonitor>(cfgc)};
  return workflow;
}
