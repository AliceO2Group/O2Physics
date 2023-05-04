// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright
// holders. All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file femtoDreamDebugTrack.cxx
/// \brief Tasks that reads the particle tables and fills QA histograms for
/// tracks \author Luca Barioglio, TU München, luca.barioglio@cern.ch

#include "DataFormatsParameters/GRPObject.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/StepTHn.h"
#include "Framework/runDataProcessing.h"

#include "FemtoDreamEventHisto.h"
#include "FemtoDreamParticleHisto.h"
#include "FemtoUtils.h"
#include "PWGCF/DataModel/FemtoDerived.h"

using namespace o2;
using namespace o2::analysis::femtoDream;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

namespace
{
static constexpr int nCuts = 5;
static const std::vector<std::string> cutNames{"MaxPt", "PIDthr", "nSigmaTPC",
                                               "nSigmaTPCTOF", "MaxP"};
static const float cutsTable[1][nCuts] = {{4.05f, 0.75f, 3.5f, 3.5f, 100.f}};

} // namespace

struct femtoDreamDebugTrack {
  SliceCache cache;

  Configurable<LabeledArray<float>> cfgCutTable{"cfgCutTable", {cutsTable[0], nCuts, cutNames}, "Particle selections"};
  Configurable<int> cfgNspecies{"ccfgNspecies", 4, "Number of particle spieces with PID info"};
  Configurable<bool> ConfIsMC{"ConfIsMC", false, "Enable additional Histogramms in the case of a MonteCarlo Run"};
  Configurable<int> ConfPDGCodePartOne{"ConfPDGCodePartOne", 2212, "Particle 1 - PDG code"};
  Configurable<uint32_t> ConfCutPartOne{"ConfCutPartOne", 5542474, "Particle 1 - Selection bit from cutCulator"};
  Configurable<std::vector<int>> ConfPIDPartOne{"ConfPIDPartOne", std::vector<int>{2}, "Particle 1 - Read from cutCulator"};
  Configurable<std::vector<float>> ConfTrkPIDnSigmaMax{"ConfTrkPIDnSigmaMax", std::vector<float>{3.5f, 3.f, 2.5f}, "This configurable needs to be the same as the one used in the producer task"};
  ConfigurableAxis CfgTempFitVarBins{"CfgDTempFitVarBins", {300, -0.15, 0.15}, "binning of the TempFitVar in the pT vs. TempFitVar plot"};
  ConfigurableAxis CfgTempFitVarpTBins{"CfgTempFitVarpTBins", {20, 0.5, 4.05}, "pT binning of the pT vs. TempFitVar plot"};

  using FemtoFullParticles = soa::Join<aod::FemtoDreamParticles, aod::FemtoDreamDebugParticles>;
  Partition<FemtoFullParticles> partsOne = (aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kTrack)) && ((aod::femtodreamparticle::cut & ConfCutPartOne) == ConfCutPartOne);
  Preslice<FemtoFullParticles> perColReco = aod::femtodreamparticle::femtoDreamCollisionId;

  using FemtoFullParticlesMC = soa::Join<aod::FemtoDreamParticles, aod::FemtoDreamDebugParticles, aod::FemtoDreamMCLabels>;
  Partition<FemtoFullParticlesMC> partsOneMC = (aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kTrack)) && ((aod::femtodreamparticle::cut & ConfCutPartOne) == ConfCutPartOne);
  Preslice<FemtoFullParticlesMC> perColGen = aod::femtodreamparticle::femtoDreamCollisionId;

  /// Histogramming for Event
  FemtoDreamEventHisto eventHisto;

  FemtoDreamParticleHisto<aod::femtodreamparticle::ParticleType::kTrack> trackHisto;

  /// The configurables need to be passed to an std::vector
  std::vector<int> vPIDPartOne;
  std::vector<float> kNsigma;

  /// Histogram output
  HistogramRegistry qaRegistry{"TrackQA", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry FullQaRegistry{"FullTrackQA", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext&)
  {
    eventHisto.init(&qaRegistry);
    trackHisto.init(&qaRegistry, CfgTempFitVarpTBins, CfgTempFitVarBins, ConfIsMC, true);
    vPIDPartOne = ConfPIDPartOne.value;
    kNsigma = ConfTrkPIDnSigmaMax.value;
  }

  /// Porduce QA plots for sigle track selection in FemtoDream framework
  template <bool isMC, typename PartitionType>
  void FillDebugHistos(o2::aod::FemtoDreamCollision& col, PartitionType& groupPartsOne)
  {
    eventHisto.fillQA(col);
    for (auto& part : groupPartsOne) {
      if (part.p() > cfgCutTable->get("MaxP") ||
          part.pt() > cfgCutTable->get("MaxPt")) {
        continue;
      }
      if (!isFullPIDSelected(part.pidcut(), part.p(), cfgCutTable->get("PIDthr"), vPIDPartOne, cfgNspecies, kNsigma, cfgCutTable->get("nSigmaTPC"), cfgCutTable->get("nSigmaTPCTOF"))) {
        continue;
      }
      trackHisto.fillQA<isMC, true>(part);
      // if constexpr (isMC) {
      // if (part.has_femtoDreamMCParticle()) {
      // mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("_MC/hPt_MC"), part.pt());
      // mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("_MC/hEta_MC"), part.eta());
      // mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("_MC/hPhi_MC"), part.phi());
      // mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("_MC/hPDG"), part.pdgMCTruth());
      // mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("_MC/hOrigin_MC"), part.partOriginMCTruth());
      // } else {
      // FullQaRegistry.fill(HIST("FullTrackQA_MC/hNoMCtruthCounter"), 1);
      // }
      // }
    }
  }

  /// process function when runnning over data/ Monte Carlo reconstructed only
  /// \param col subscribe to FemtoDreamCollision table
  /// \param parts subscribe to FemtoDreamParticles table
  void processData(o2::aod::FemtoDreamCollision& col, FemtoFullParticles& parts)
  {
    auto groupPartsOne = partsOne->sliceByCached(aod::femtodreamparticle::femtoDreamCollisionId, col.globalIndex(), cache);
    FillDebugHistos<false>(col, groupPartsOne);
  }
  PROCESS_SWITCH(femtoDreamDebugTrack, processData, "Enable Debug processing for Monte Carlo", true);

  /// process function when runnning over Monte Carlo with MC truth enabled

  /// \param col subscribe to FemtoDreamCollision table
  /// \param parts subscribe to the joined table of FemtoDreamParticles and FemtoDreamMCLabels table
  /// \param FemtoDramMCParticles subscribe to the table containing the Monte Carlo Truth information
  void processMC(o2::aod::FemtoDreamCollision& col, FemtoFullParticlesMC& parts, o2::aod::FemtoDreamMCParticles&)
  {
    auto groupPartsOne = partsOneMC->sliceByCached(aod::femtodreamparticle::femtoDreamCollisionId, col.globalIndex(), cache);
    FillDebugHistos<true>(col, groupPartsOne);
  }
  PROCESS_SWITCH(femtoDreamDebugTrack, processMC, "Enable Debug processing for Monte Carlo", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<femtoDreamDebugTrack>(cfgc),
  };
  return workflow;
}
