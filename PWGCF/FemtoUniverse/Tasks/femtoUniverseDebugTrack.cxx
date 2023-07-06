// Copyright 2019-2022 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file femtoUniverseDebugTrack.cxx
/// \brief Tasks that reads the particle tables and fills QA histograms for tracks
/// \author Luca Barioglio, TU München, luca.barioglio@cern.ch
/// \author Zuzanna Chochulska, WUT Warsaw, zuzanna.chochulska.stud@pw.edu.pl

#include "DataFormatsParameters/GRPObject.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/StepTHn.h"
#include "Framework/runDataProcessing.h"

#include "PWGCF/FemtoUniverse/Core/FemtoUniverseEventHisto.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseParticleHisto.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUtils.h"
#include "PWGCF/FemtoUniverse/DataModel/FemtoDerived.h"

using namespace o2;
using namespace o2::analysis::femtoUniverse;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

namespace
{
static constexpr int nCuts = 5;
static const std::vector<std::string> cutNames{"MaxPt", "PIDthr", "nSigmaTPC", "nSigmaTPCTOF", "MaxP"};
static const float cutsTable[1][nCuts] = {{4.05f, 0.75f, 3.f, 3.f, 100.f}};

} // namespace

struct femtoUniverseDebugTrack {
  SliceCache cache;

  Configurable<LabeledArray<float>> ConfCutTable{"ConfCutTable", {cutsTable[0], nCuts, cutNames}, "Particle selections"};
  Configurable<int> ConfNspecies{"ConfNspecies", 2, "Number of particle spieces with PID info"};
  Configurable<bool> ConfIsMC{"ConfIsMC", false, "Enable additional Histogramms in the case of a MonteCarlo Run"};
  Configurable<int> ConfPDGCodePartOne{"ConfPDGCodePartOne", 2212, "Particle 1 - PDG code"};
  Configurable<uint32_t> ConfCutPartOne{"ConfCutPartOne", 5542474, "Particle 1 - Selection bit from cutCulator"};
  Configurable<int> ConfPIDPartOne{"ConfPIDPartOne", 1, "Particle 1 - Read from cutCulator"};
  Configurable<std::vector<float>> ConfTrkPIDnSigmaMax{"ConfTrkPIDnSigmaMax", std::vector<float>{3.5f, 3.f, 2.5f}, "This configurable needs to be the same as the one used in the producer task"};
  ConfigurableAxis ConfTempFitVarBins{"ConfDTempFitVarBins", {300, -0.15, 0.15}, "binning of the TempFitVar in the pT vs. TempFitVar plot"};
  ConfigurableAxis ConfTempFitVarpTBins{"ConfTempFitVarpTBins", {20, 0.5, 4.05}, "pT binning of the pT vs. TempFitVar plot"};

  using FemtoFullParticles = soa::Join<aod::FDParticles, aod::FDExtParticles>;
  Partition<FemtoFullParticles> partsOne = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kTrack)) && ((aod::femtouniverseparticle::cut & ConfCutPartOne) == ConfCutPartOne);
  Preslice<FemtoFullParticles> perColReco = aod::femtouniverseparticle::fdCollisionId;

  using FemtoFullParticlesMC = soa::Join<aod::FDParticles, aod::FDExtParticles, aod::FDMCLabels>;
  Partition<FemtoFullParticlesMC> partsOneMC = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kTrack)) && ((aod::femtouniverseparticle::cut & ConfCutPartOne) == ConfCutPartOne);
  Preslice<FemtoFullParticlesMC> perColGen = aod::femtouniverseparticle::fdCollisionId;

  /// Histogramming for Event
  FemtoUniverseEventHisto eventHisto;

  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kTrack> trackHisto;

  /// The configurables need to be passed to an std::vector
  int vPIDPartOne;
  std::vector<float> kNsigma;

  /// Histogram output
  HistogramRegistry qaRegistry{"TrackQA", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext&)
  {
    eventHisto.init(&qaRegistry);
    trackHisto.init(&qaRegistry, ConfTempFitVarpTBins, ConfTempFitVarBins, ConfIsMC, ConfPDGCodePartOne.value, true);
    vPIDPartOne = ConfPIDPartOne.value;
    kNsigma = ConfTrkPIDnSigmaMax.value;
  }

  /// Porduce QA plots for sigle track selection in FemtoUniverse framework
  template <bool isMC, typename PartitionType>
  void FillDebugHistos(o2::aod::FDCollision& col, PartitionType& groupPartsOne)
  {
    eventHisto.fillQA(col);
    for (auto& part : groupPartsOne) {
      if (part.p() > ConfCutTable->get("MaxP") || part.pt() > ConfCutTable->get("MaxPt")) {
        continue;
      }
      if (!isFullPIDSelected(part.pidcut(), part.p(), ConfCutTable->get("PIDthr"), vPIDPartOne, ConfNspecies, kNsigma, ConfCutTable->get("nSigmaTPC"), ConfCutTable->get("nSigmaTPCTOF"))) {
        continue;
      }
      trackHisto.fillQA<isMC, true>(part);
    }
  }

  /// process function when runnning over data/ Monte Carlo reconstructed only
  /// \param col subscribe to FemtoUniverseCollision table
  /// \param parts subscribe to FemtoUniverseParticles table
  void processData(o2::aod::FDCollision& col, FemtoFullParticles& parts)
  {
    auto groupPartsOne = partsOne->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    FillDebugHistos<false>(col, groupPartsOne);
  }
  PROCESS_SWITCH(femtoUniverseDebugTrack, processData, "Enable Debug processing for Monte Carlo", true);

  /// process function when runnning over Monte Carlo with MC truth enabled

  /// \param col subscribe to FemtoUniverseCollision table
  /// \param parts subscribe to the joined table of FemtoUniverseParticles and FemtoUniverseMCLabels table
  /// \param FemtoDramMCParticles subscribe to the table containing the Monte Carlo Truth information
  void processMC(o2::aod::FDCollision& col, FemtoFullParticlesMC& parts, o2::aod::FDMCParticles&)
  {
    auto groupPartsOne = partsOneMC->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    FillDebugHistos<true>(col, groupPartsOne);
  }
  PROCESS_SWITCH(femtoUniverseDebugTrack, processMC, "Enable Debug processing for Monte Carlo", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<femtoUniverseDebugTrack>(cfgc),
  };
  return workflow;
}
