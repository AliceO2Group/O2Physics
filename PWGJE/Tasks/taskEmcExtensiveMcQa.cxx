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

/// \file taskEmcExtensiveMcQa.cxx
/// \brief Exensive monitoring task for EMCal clusters in MC
/// \author Marvin Hemmer <marvin.hemmer@cern.ch>, Goethe University Frankfurt
/// \since 31.07.2025

#include "PWGJE/DataModel/EMCALClusters.h"
// HF headers for event selection
#include "PWGHF/Core/CentralityEstimation.h"
#include "PWGHF/Utils/utilsEvSelHf.h"

#include "Common/CCDB/ctpRateFetcher.h"
#include "Common/DataModel/EventSelection.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <EMCALBase/Geometry.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <TPDGCode.h>

#include <algorithm>
#include <array>
#include <climits>
#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>

using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants;
using namespace o2::hf_evsel;
using namespace o2::hf_centrality;
using CollisionEvSels = o2::soa::Join<o2::aod::Collisions, o2::aod::EvSels>;
using BcEvSelIt = o2::soa::Join<o2::aod::BCs, o2::aod::BcSels>::iterator;
using SelectedClusters = o2::soa::Filtered<o2::soa::Join<o2::aod::EMCALClusters, o2::aod::EMCALMCClusters>>;

namespace poi
{
enum PoI {
  kPhoton = 0,
  kElectronPrim = 1,
  kElectronSec = 2,
  kMuon = 3,
  kHadronCharge = 4,
  kHadronNeutral = 5,
  kNPoI = 6
};
} // namespace poi

/// \struct TaskEmcExtensiveMcQa
struct TaskEmcExtensiveMcQa {

  static constexpr int NSM = 20; // there 20 supermodlues for the EMCal
  std::array<int, 6> arrPDGHadronNeutral = {kNeutron, kK0Short, kK0Long, kLambda0, physics::kXi0, kSigma0};

  SliceCache cache;
  Preslice<SelectedClusters> psClusterPerCollision = o2::aod::emcalcluster::collisionId;
  Preslice<o2::aod::EMCALClusterCells> perCluster = o2::aod::emcalclustercell::emcalclusterId;

  HistogramRegistry mHistManager{"EMCalExtensiveMCQAHistograms"};

  o2::emcal::Geometry* mGeometry = nullptr;
  o2::framework::Service<o2::ccdb::BasicCCDBManager> ccdb;

  o2::ctpRateFetcher rateFetcher;
  HfEventSelection hfEvSel;
  HfEventSelectionMc hfEvSelMc;

  // configurable parameters
  Configurable<bool> applyEvSels{"applyEvSels", true, "Flag to apply event selection."};
  Configurable<int> clusterDefinition{"clusterDefinition", 10, "cluster definition to be selected, e.g. 10=kV3Default"};
  Configurable<std::string> ctpFetcherSource{"ctpFetcherSource", "T0VTX", "Source for CTP rate fetching, e.g. T0VTX, T0CE, T0SC, ZNC (hadronic)"};
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};

  // configurable axis
  ConfigurableAxis nClustersBinning{"nClustersBinning", {201, -0.5, 200.5}, "binning for the number of clusters"};

  ConfigurableAxis clusterEnergy{"clusterEnergy", {100, 0., 10}, "binning for the cluster energy in GeV"};
  ConfigurableAxis clusterM02{"clusterM02", {100, 0., 2.0}, "binning for the cluster M02"};
  ConfigurableAxis clusterM20{"clusterM20", {100, 0., 2.0}, "binning for the cluster M20"};
  ConfigurableAxis clusterNCellBinning{"clusterNCellBinning", {100, 0.5, 100.5}, "binning for the number of cells per cluster"};
  ConfigurableAxis clusterOriginRadius{"clusterOriginRadius", {225, 0., 450}, "binning for the radial original point of the main contributor of a cluster"};
  ConfigurableAxis clusterNContributor{"clusterNContributor", {20, 0.5, 20.5}, "binning for the number of contributor of a cluster"};
  ConfigurableAxis clusterEnergyRatio{"clusterEnergyRatio", {100, 0., 10.}, "binning for ratio of the deposited energy of the leading particle to its generated momentum cluster"};
  ConfigurableAxis collisionCent{"collisionCent", {10, 0., 100.}, "binning for the event centrality"};
  ConfigurableAxis clusterEtaReso{"clusterEtaReso", {100, -0.1, 0.1}, "binning for cluster position resolution in eta"};
  ConfigurableAxis clusterPhiReso{"clusterPhiReso", {100, -0.1, 0.1}, "binning for cluster position resolution in phi"};

  std::vector<float> mCellTime;

  /// \brief Create output histograms and initialize geometry
  void init(InitContext const&)
  {
    // load geometry just in case we need it
    mGeometry = o2::emcal::Geometry::GetInstanceFromRunNumber(300000);

    // create common axes
    const AxisSpec numberClustersAxis{nClustersBinning, "#it{N}_{cl}/ #it{N}_{event}"};
    const AxisSpec axisParticle = {poi::kNPoI, -0.5f, +poi::kNPoI - 0.5f, ""};
    const AxisSpec axisEnergy{clusterEnergy, "#it{E}_{cl} (GeV)"};
    const AxisSpec axisM02{clusterM02, "#it{M}_{02}"};
    const AxisSpec axisM20{clusterM20, "#it{M}_{20}"};
    const AxisSpec axisNCell{clusterNCellBinning, "#it{N}_{cells}"};
    const AxisSpec axisRadius{clusterOriginRadius, "#it{R}_{origin} (cm)"};
    const AxisSpec axisNContributor{clusterNContributor, "#it{N}_{particles}"};
    const AxisSpec axisCent{collisionCent, "cent (%)"};
    const AxisSpec axisLeadingEnergy{clusterEnergy, "#it{E}_{lead} (GeV)"};
    const AxisSpec axisLeadingGenMomentum{clusterEnergy, "#it{p}_{lead, gen} (GeV/#it{c})"};
    const AxisSpec axisLeadingRatio{clusterEnergy, "#it{E}_{lead}/#it{p}_{lead, gen} (#it{c})"};
    const AxisSpec axisEtaReso{clusterEtaReso, "#Delta#eta"};
    const AxisSpec axisPhiReso{clusterPhiReso, "#Delta#varphi (rad)"};

    const AxisSpec axisSM{{20, -0.5, 19.5}, "SM"};

    // create histograms

    // event properties
    mHistManager.add("numberOfClustersEvents", "number of clusters per event (selected events)", HistType::kTH1D, {numberClustersAxis});

    // cluster properties (matched clusters)
    mHistManager.add("hSparseClusterQA", "THnSparse for Cluster QA", HistType::kTHnSparseF, {axisEnergy, axisM02, axisM20, axisNCell, axisRadius, axisParticle, axisNContributor, axisCent});
    mHistManager.add("hSparseClusterContributors", "THnSparse with cluster contributors and energies", HistType::kTHnSparseF, {axisEnergy, axisParticle, axisNContributor, axisLeadingEnergy, axisLeadingGenMomentum, axisLeadingRatio, axisCent});
    mHistManager.add("clusterEtaPhi", "Eta and phi of cluster", HistType::kTH2F, {{140, -0.7, 0.7}, {360, 0, o2::constants::math::TwoPI}});

    mHistManager.add("hSparsePosReso", "THnSparse for cluster position resolution", HistType::kTHnSparseF, {axisEnergy, axisEtaReso, axisPhiReso, axisNCell, axisSM, axisParticle, axisNContributor, axisCent});

    hfEvSel.addHistograms(mHistManager);

    ccdb->setURL(ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
  }

  template <typename Coll>
  bool isCollSelected(const Coll& coll, float& cent)
  {
    const auto rejectionMask = hfEvSel.getHfCollisionRejectionMask<true, o2::hf_centrality::CentralityEstimator::None, o2::aod::BCsWithTimestamps>(coll, cent, ccdb, mHistManager);
    /// monitor the satisfied event selections
    hfEvSel.fillHistograms(coll, rejectionMask, cent);
    return rejectionMask == 0;
  }

  /// \brief returns the PoI type of a mcparticle
  /// \param mcparticle is the mcparticle we want to find the PoI type
  /// \param mcparticles table containing the mcparticles
  /// \return PoI type of the given mcparticle
  template <typename T>
  int findPoIType(T const& mcparticle)
  {
    auto pdgValue = std::abs(mcparticle.pdgCode());
    switch (pdgValue) {
      case kGamma: {
        return poi::kPhoton;
      }
      case kElectron: {
        if (mcparticle.isPhysicalPrimary()) {
          return poi::kElectronPrim;
        } else {
          return poi::kElectronSec;
        }
      }
      case kMuonMinus: {
        return poi::kMuon;
      }
      default: {
        auto it = std::find(arrPDGHadronNeutral.begin(), arrPDGHadronNeutral.end(), pdgValue);
        if (it != arrPDGHadronNeutral.end()) {
          return poi::kHadronNeutral;
        }
        return poi::kHadronCharge;
      }
    }
  }

  Filter clusterDefinitionSelection = (o2::aod::emcalcluster::definition == clusterDefinition);

  /// \brief Process EMCAL clusters that are matched to a collisions
  void processCollisions(CollisionEvSels const& collisions, SelectedClusters const& clusters, McParticles const& /*mcparticles*/)
  {

    for (const auto& collision : collisions) {
      float cent = -1.f;
      if (applyEvSels && !isCollSelected(collision, cent)) {
        continue;
      }

      auto groupedClusters = clusters.sliceBy(psClusterPerCollision, collision.globalIndex());
      mHistManager.fill(HIST("numberOfClustersEvents"), groupedClusters.size());

      for (const auto& cluster : groupedClusters) {
        mHistManager.fill(HIST("clusterEtaPhi"), cluster.eta(), cluster.phi());
        // axisEnergy, axisTime, axisM02, axisNCell, axisRadius, axisParticle
        if (cluster.mcParticle().size() == 0) {
          LOG(info) << "Somehow cluster.mcParticle().size() == 0!";
          continue;
        }
        auto mainMcParticle = cluster.mcParticle_as<McParticles>()[0];
        float radius = std::hypot(mainMcParticle.vx(), mainMcParticle.vy());
        float momentum = mainMcParticle.p();
        float leadingEnergy = cluster.energy() * cluster.amplitudeA()[0];
        float leadingFraction = leadingEnergy / momentum;
        float dEta = cluster.eta() - mainMcParticle.eta();
        float dPhi = cluster.phi() - mainMcParticle.phi();
        int iSM = mGeometry->SuperModuleNumberFromEtaPhi(cluster.eta(), cluster.phi());
        mHistManager.fill(HIST("hSparseClusterQA"), cluster.energy(), cluster.m02(), cluster.m20(), cluster.nCells(), radius, findPoIType(mainMcParticle), cluster.mcParticle().size(), cent);
        mHistManager.fill(HIST("hSparseClusterContributors"), cluster.energy(), findPoIType(mainMcParticle), cluster.mcParticle().size(), leadingEnergy, momentum, leadingFraction, cent);
        mHistManager.fill(HIST("hSparsePosReso"), cluster.energy(), dEta, dPhi, cluster.nCells(), iSM, findPoIType(mainMcParticle), cluster.mcParticle().size(), cent);
      }
    }
  }
  PROCESS_SWITCH(TaskEmcExtensiveMcQa, processCollisions, "Process clusters from collision", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<TaskEmcExtensiveMcQa>(cfgc)};
  return workflow;
}
