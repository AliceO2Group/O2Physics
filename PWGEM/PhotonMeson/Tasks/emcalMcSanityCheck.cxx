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

/// \file emcalMcSanityCheck.cxx
/// \brief Task to test MC info for EMCal during different stages of out workflow
/// \author M. Hemmer, marvin.hemmer@cern.ch

#include "PWGEM/PhotonMeson/DataModel/GammaTablesRedux.h"
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGEM/PhotonMeson/Utils/emcalHistoDefinitions.h"
#include "PWGJE/DataModel/EMCALClusters.h"

#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

constexpr float PointOnePercent = 0.001f;
constexpr float OnePercent = 0.01f;
constexpr float TenPercent = 0.1f;
constexpr float NinetyPercent = 0.9f;
constexpr float OneHundredTenPercent = 1.1f;

// LSB from packing in MinClusters is 1 MeV
constexpr float EnergyPackingLsb = 1.f / o2::aod::emcdownscaling::downscalingFactors[o2::aod::emcdownscaling::kEnergy];

// Rounding error from packing should be half the LSB so 0.5 MeV
constexpr float EnergyComparisonEpsilon = 0.5f * EnergyPackingLsb;

// uint16_t saturates at 65535 so max value of energy for MinClusters should be 65.535 GeV
constexpr float EnergyPackingSaturationValue = static_cast<float>(std::numeric_limits<uint16_t>::max()) * EnergyPackingLsb;

enum class EnergyComparisonStatus : uint8_t {
  Match = 0,    // agrees within half a packing LSB
  Mismatch = 1, // disagrees beyond epsilon, and not explained by saturation
  Saturated = 2 // reference energy exceeded the uint16_t packing range
};

enum class McParticleComparisonStatus : uint8_t {
  Match = 0,       // same particle species, energy within epsilon, cluster fraction same
  DiffSpecies = 1, // not the same species
  DiffEnergy = 2,  // energy not within epsilon
  DiffFraction = 3 // cluster fraction not within epsilon
};

enum class EnergyRatioStatus : uint8_t {
  PointOnePercent = 0,
  OnePercent,
  TenPercent,
  NinetyPercent,
  Ok,
  OneHundredTenPercent
};

struct EmcalMcSanityCheck {

  using BeforeSkimmerCluster = soa::Join<EMCALClusters, EMCALMCClusters>;
  using AfterSkimmerCluster = soa::Join<MinClusters, EMCClusterMCLabels_001>;
  using AfterAssociateCluster = soa::Join<MinClusters, EMEMCClusterMCLabels_001>;
  using BeforeAfterAssociateCluster = soa::Join<MinClusters, EMCClusterMCLabels_001, EMEMCClusterMCLabels_001>;

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};

  void init(InitContext&)
  {
    const AxisSpec axisParticleClusterFracRatioStatus{6, -0.5, 5.5};

    auto hParticleClusterFracRatioStatus = registry.add<TH1>("BeforeSkimmer/hParticleClusterFracRatioStatus", "hParticleClusterFracRatioStatus;;counts", HistType::kTH1D, {axisParticleClusterFracRatioStatus});
    hParticleClusterFracRatioStatus->GetXaxis()->SetBinLabel(1, "#it{E}_{part}/(frac#it{E}_{clus}) < 0.001");
    hParticleClusterFracRatioStatus->GetXaxis()->SetBinLabel(2, "0.001 #leq #it{E}_{part}/(frac#it{E}_{clus}) < 0.01");
    hParticleClusterFracRatioStatus->GetXaxis()->SetBinLabel(3, "0.01 #leq #it{E}_{part}/(frac#it{E}_{clus}) < 0.1");
    hParticleClusterFracRatioStatus->GetXaxis()->SetBinLabel(4, "0.1 #leq #it{E}_{part}/(frac#it{E}_{clus}) < 0.9");
    hParticleClusterFracRatioStatus->GetXaxis()->SetBinLabel(5, "0.9 #leq #it{E}_{part}/(frac#it{E}_{clus}) < 1.1");
    hParticleClusterFracRatioStatus->GetXaxis()->SetBinLabel(6, "#it{E}_{part}/(frac#it{E}_{clus}) > 1.1");

    if (doprocessAfterSkimmer) {
      registry.addClone("BeforeSkimmer/", "AfterSkimmer/");
    }

    if (doprocessAfterAssociation) {
      registry.addClone("BeforeSkimmer/", "AfterAssociate/");
    }

    if (doprocessCompareBeforeAfterAssociate) {
      auto hMcParticleStatus = registry.add<TH1>("hMcParticleStatus", "hMcParticleStatus;;counts", HistType::kTH1D, {{4, -0.5, 3.5}});
      hMcParticleStatus->GetXaxis()->SetBinLabel(static_cast<int>(McParticleComparisonStatus::Match) + 1, "Match");
      hMcParticleStatus->GetXaxis()->SetBinLabel(static_cast<int>(McParticleComparisonStatus::DiffSpecies) + 1, "DiffSpecies");
      hMcParticleStatus->GetXaxis()->SetBinLabel(static_cast<int>(McParticleComparisonStatus::DiffEnergy) + 1, "DiffEnergy");
      hMcParticleStatus->GetXaxis()->SetBinLabel(static_cast<int>(McParticleComparisonStatus::DiffFraction) + 1, "DiffFraction");
    }
  }; // end init

  EnergyRatioStatus getEnergyRatioStatus(float clusterE, float mcParticleE, float frac)
  {
    float ratio = mcParticleE / (clusterE * frac);
    // the most likely outcome due to hadrons not depositing their full energy
    if (ratio > OneHundredTenPercent) [[likely]] {
      return EnergyRatioStatus::OneHundredTenPercent;
    }
    if (NinetyPercent <= ratio && ratio < OneHundredTenPercent) {
      return EnergyRatioStatus::Ok;
    }
    if (TenPercent <= ratio && ratio < NinetyPercent) {
      return EnergyRatioStatus::NinetyPercent;
    }
    if (OnePercent <= ratio && ratio < TenPercent) [[unlikely]] {
      return EnergyRatioStatus::TenPercent;
    }
    if (PointOnePercent <= ratio && ratio < OnePercent) [[unlikely]] {
      return EnergyRatioStatus::OnePercent;
    }
    // everything else now has to be below 0.1%
    return EnergyRatioStatus::PointOnePercent;
  }

  // referenceEnergy: full-precision energy (e.g. BeforeSkimmer's cluster.energy()).
  // compressedEnergy: energy unpacked from a uint16_t-packed table (e.g. MinClusters' cluster.e()).
  EnergyComparisonStatus compareClusterEnergy(float referenceEnergy, float compressedEnergy)
  {
    if (referenceEnergy > EnergyPackingSaturationValue) {
      return EnergyComparisonStatus::Saturated;
    }
    return std::abs(referenceEnergy - compressedEnergy) <= EnergyComparisonEpsilon
             ? EnergyComparisonStatus::Match
             : EnergyComparisonStatus::Mismatch;
  }

  void processBeforeSkimmer(BeforeSkimmerCluster const& clusters, McParticles_001 const& mcParticles)
  {
    if (clusters.size() == 0 || mcParticles.size() == 0) {
      return;
    }
    auto mcParticle = mcParticles.begin();

    for (const auto& cluster : clusters) {
      if (!cluster.has_mcParticle()) {
        continue;
      }
      const auto& ids = cluster.mcParticleIds();
      const auto& fracs = cluster.amplitudeA();
      for (std::size_t i = 0; i < ids.size(); ++i) {
        const auto& id = ids[i];
        const auto& frac = fracs[i];
        mcParticle.setCursor(id);
        registry.fill(HIST("BeforeSkimmer/hParticleClusterFracRatioStatus"), static_cast<int>(getEnergyRatioStatus(cluster.energy(), mcParticle.e(), frac)));
      }
    }
  }
  PROCESS_SWITCH(EmcalMcSanityCheck, processBeforeSkimmer, "Process EMCal cluster information before the skimmerGammaCalo", true);

  void processAfterSkimmer(AfterSkimmerCluster const& clusters, McParticles_001 const& mcParticles)
  {
    if (clusters.size() == 0 || mcParticles.size() == 0) {
      return;
    }
    auto mcParticle = mcParticles.begin();

    for (const auto& cluster : clusters) {
      if (!cluster.has_mcParticle()) {
        continue;
      }
      const auto& ids = cluster.mcParticleIds();
      const auto& fracs = cluster.amplitude();
      for (std::size_t i = 0; i < ids.size(); ++i) {
        const auto& id = ids[i];
        const auto& frac = fracs[i];
        mcParticle.setCursor(id);
        registry.fill(HIST("AfterSkimmer/hParticleClusterFracRatioStatus"), static_cast<int>(getEnergyRatioStatus(cluster.e(), mcParticle.e(), frac)));
      }
    }
  }
  PROCESS_SWITCH(EmcalMcSanityCheck, processAfterSkimmer, "Process EMCal cluster information after the skimmerGammaCalo", true);

  void processAfterAssociation(AfterAssociateCluster const& clusters, EMMCParticles_001 const& mcParticles)
  {
    if (clusters.size() == 0 || mcParticles.size() == 0) {
      return;
    }
    auto mcParticle = mcParticles.begin();

    for (const auto& cluster : clusters) {
      if (!cluster.has_emmcparticle()) {
        continue;
      }
      const auto& ids = cluster.emmcparticleIds();
      const auto& fracs = cluster.amplitude();
      for (std::size_t i = 0; i < ids.size(); ++i) {
        const auto& id = ids[i];
        const auto& frac = fracs[i];
        mcParticle.setCursor(id);
        registry.fill(HIST("AfterAssociate/hParticleClusterFracRatioStatus"), static_cast<int>(getEnergyRatioStatus(cluster.e(), mcParticle.e(), frac)));
      }
    }
  }
  PROCESS_SWITCH(EmcalMcSanityCheck, processAfterAssociation, "Process EMCal cluster information before the associateMCinfoPhoton", false);

  void processCompareBeforeAfterAssociate(MinClusters const& clusters,
                                          EMCClusterMCLabels_001 const& beforeLabels,
                                          EMEMCClusterMCLabels_001 const& afterLabels,
                                          McParticles_001 const& mcParticles,
                                          EMMCParticles_001 const& emMcParticles)
  {
    if (clusters.size() == 0 || mcParticles.size() == 0 || emMcParticles.size() == 0) {
      return;
    }

    auto mcParticle = mcParticles.begin();
    auto emMcParticle = emMcParticles.begin();

    for (int64_t i = 0; i < clusters.size(); ++i) {
      auto before = beforeLabels.iteratorAt(i);
      auto after = afterLabels.iteratorAt(i);

      if (!after.has_emmcparticle() || !before.has_mcParticle()) {
        continue;
      }

      if (after.emmcparticleIds().size() != before.mcParticleIds().size()) {
        LOG(fatal) << "Number of EmMcParticles and McParticles does not match (" << after.emmcparticleIds().size() << " vs " << before.mcParticleIds().size() << ")";
      }
      for (size_t id = 0; id < after.emmcparticleIds().size(); ++id) {
        bool isSameSpecies = true;
        bool isSameEnergy = true;
        bool hasSameFraction = true;

        mcParticle.setCursor(before.mcParticleIds()[id]);
        emMcParticle.setCursor(after.emmcparticleIds()[id]);

        if (mcParticle.pdgCode() != emMcParticle.pdgCode()) {
          isSameSpecies = false;
        }
        if (mcParticle.e() != emMcParticle.e()) {
          isSameEnergy = false;
        }
        if (before.amplitude()[id] != after.amplitude()[id]) {
          hasSameFraction = false;
        }

        if (isSameSpecies && isSameEnergy && hasSameFraction) {
          registry.fill(HIST("hMcParticleStatus"), 0);
        }
        if (!isSameSpecies) {
          registry.fill(HIST("hMcParticleStatus"), 1);
        }
        if (!isSameEnergy) {
          registry.fill(HIST("hMcParticleStatus"), 2);
        }
        if (!hasSameFraction) {
          registry.fill(HIST("hMcParticleStatus"), 3);
        }
      }
    }
  }
  PROCESS_SWITCH(EmcalMcSanityCheck, processCompareBeforeAfterAssociate, "Process and compare EMCal cluster and McParticle information before and after the associateMCinfoPhoton", false);
}; // End struct EmcalMcSanityCheck

WorkflowSpec defineDataProcessing(ConfigContext const& context)
{
  return WorkflowSpec{adaptAnalysisTask<EmcalMcSanityCheck>(context)};
}
