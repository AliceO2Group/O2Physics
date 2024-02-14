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

/// \file femtoDreamDebugTrack.cxx
/// \brief Tasks that reads the particle tables and fills QA histograms for tracks
/// \author Luca Barioglio, TU München, luca.barioglio@cern.ch

#include <Framework/Expressions.h>
#include "DataFormatsParameters/GRPObject.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/StepTHn.h"
#include "Framework/runDataProcessing.h"

#include "PWGCF/FemtoDream/Core/femtoDreamEventHisto.h"
#include "PWGCF/FemtoDream/Core/femtoDreamParticleHisto.h"
#include "PWGCF/DataModel/FemtoDerived.h"

using namespace o2;
using namespace o2::analysis::femtoDream;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

struct femtoDreamDebugTrack {
  SliceCache cache;

  Configurable<bool> ConfIsMC{"ConfIsMC", false, "Enable additional Histogramms in the case of a MonteCarlo Run"};
  Configurable<int> ConfTrk1_PDGCode{"ConfTrk1_PDGCode", 2212, "Particle 1 - PDG code"};
  Configurable<aod::femtodreamparticle::cutContainerType> ConfTrk1_CutBit{"ConfTrk1_CutBit", 5542474, "Particle 1 - Selection bit from cutCulator"};
  Configurable<aod::femtodreamparticle::cutContainerType> ConfTrk1_TPCBit{"ConfTrk1_TPCBit", 1, "Particle 1 - Read from cutCulator"};
  Configurable<aod::femtodreamparticle::cutContainerType> ConfTrk1_TPCTOFBit{"ConfTrk1_TPCTOFBit", 1, "Particle 1 - Read from cutCulator"};
  Configurable<float> ConfTrk1_PIDThres{"ConfTrk1_PIDThres", 0.75, "Particle 1 - Read from cutCulator"};

  Configurable<bool> ConfOptCorrelatedPlots{"ConfOptCorrelatedPlots", false, "Enable additional three dimensional histogramms. High memory consumption. Use for debugging"};
  ConfigurableAxis ConfBinmult{"ConfBinmult", {VARIABLE_WIDTH, 0.0f, 4.0f, 8.0f, 12.0f, 16.0f, 20.0f, 24.0f, 28.0f, 32.0f, 36.0f, 40.0f, 44.0f, 48.0f, 52.0f, 56.0f, 60.0f, 64.0f, 68.0f, 72.0f, 76.0f, 80.0f, 84.0f, 88.0f, 92.0f, 96.0f, 100.0f, 200.0f}, "multiplicity Binning"};
  ConfigurableAxis ConfBinmultPercentile{"ConfBinmultPercentile", {10, 0.0f, 100.0f}, "multiplicity percentile Binning"};
  ConfigurableAxis ConfBinpT{"ConfBinpT", {{240, 0, 6}}, "pT binning"};
  ConfigurableAxis ConfBineta{"ConfBineta", {{200, -1.5, 1.5}}, "eta binning"};
  ConfigurableAxis ConfBinphi{"ConfBinphi", {{200, 0, TMath::TwoPi()}}, "phi binning"};

  ConfigurableAxis ConfTempFitVarBins{"ConfTempFitVarBins", {300, -0.15, 0.15}, "Binning of the TempFitVar"};
  ConfigurableAxis ConfNsigmaTPCBins{"ConfNsigmaTPCBins", {1600, -8, 8}, "Binning of Nsigma TPC plot"};
  ConfigurableAxis ConfNsigmaTOFBins{"ConfNsigmaTOFBins", {3000, -15, 15}, "Binning of the Nsigma TOF plot"};
  ConfigurableAxis ConfNsigmaTPCTOFBins{"ConfNsigmaTPCTOFBins", {3000, -15, 15}, "Binning of the Nsigma TPC+TOF plot"};
  ConfigurableAxis ConfTPCclustersBins{"ConfTPCClustersBins", {163, -0.5, 162.5}, "Binning of TPC found clusters plot"};
  ConfigurableAxis ConfTempFitVarMomentumBins{"ConfMomentumBins", {20, 0.5, 4.05}, "pT/p_reco/p_tpc binning of the Momentum vs. TempFitVar/Nsigma plot"};
  Configurable<int> ConfTempFitVarMomentum{"ConfTempFitVarMomentum", 0, "Momentum used for binning: 0 -> pt; 1 -> preco; 2 -> ptpc"};
  ConfigurableAxis ConfDummy{"ConfDummy", {1, 0, 1}, "Dummy axis for inv mass"};

  using FemtoFullParticles = soa::Join<aod::FDParticles, aod::FDExtParticles>;
  Partition<FemtoFullParticles> partsOne =
    (aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kTrack)) &&
    ncheckbit(aod::femtodreamparticle::cut, ConfTrk1_CutBit) &&
    ifnode(aod::femtodreamparticle::pt * (nexp(aod::femtodreamparticle::eta) + nexp(-1.f * aod::femtodreamparticle::eta)) / 2.f <= ConfTrk1_PIDThres, ncheckbit(aod::femtodreamparticle::pidcut, ConfTrk1_TPCBit), ncheckbit(aod::femtodreamparticle::pidcut, ConfTrk1_TPCTOFBit));
  Preslice<FemtoFullParticles> perColReco = aod::femtodreamparticle::fdCollisionId;

  using FemtoFullParticlesMC = soa::Join<aod::FDParticles, aod::FDExtParticles, aod::FDMCLabels>;
  Partition<FemtoFullParticlesMC> partsOneMC =
    (aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kTrack)) &&
    ncheckbit(aod::femtodreamparticle::cut, ConfTrk1_CutBit) &&
    ifnode(aod::femtodreamparticle::pt * (nexp(aod::femtodreamparticle::eta) + nexp(-1.f * aod::femtodreamparticle::eta)) / 2.f <= ConfTrk1_PIDThres, ncheckbit(aod::femtodreamparticle::pidcut, ConfTrk1_TPCBit), ncheckbit(aod::femtodreamparticle::pidcut, ConfTrk1_TPCTOFBit));
  Preslice<FemtoFullParticlesMC> perColGen = aod::femtodreamparticle::fdCollisionId;

  /// Histogramming for Event
  FemtoDreamEventHisto eventHisto;

  FemtoDreamParticleHisto<aod::femtodreamparticle::ParticleType::kTrack> trackHisto;

  /// Histogram output
  HistogramRegistry qaRegistry{"TrackQA", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext&)
  {
    eventHisto.init(&qaRegistry);
    trackHisto.init(&qaRegistry, ConfBinmult, ConfBinmultPercentile, ConfBinpT, ConfBineta, ConfBinphi, ConfTempFitVarBins, ConfNsigmaTPCBins, ConfNsigmaTOFBins, ConfNsigmaTPCTOFBins, ConfTPCclustersBins, ConfDummy, ConfIsMC, ConfTrk1_PDGCode.value, true, ConfOptCorrelatedPlots);
  }

  /// Porduce QA plots for sigle track selection in FemtoDream framework
  template <bool isMC, typename PartitionType>
  void FillDebugHistos(o2::aod::FDCollision& col, PartitionType& groupPartsOne)
  {
    eventHisto.fillQA(col);
    for (auto& part : groupPartsOne) {
      trackHisto.fillQA<isMC, true>(part, col, static_cast<aod::femtodreamparticle::MomentumType>(ConfTempFitVarMomentum.value), ConfOptCorrelatedPlots);
    }
  }

  /// process function when runnning over data/ Monte Carlo reconstructed only
  /// \param col subscribe to FemtoDreamCollision table
  /// \param parts subscribe to FemtoDreamParticles table
  void processData(o2::aod::FDCollision& col, FemtoFullParticles& parts)
  {
    auto groupPartsOne = partsOne->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
    FillDebugHistos<false>(col, groupPartsOne);
  }
  PROCESS_SWITCH(femtoDreamDebugTrack, processData, "Enable Debug processing for Monte Carlo", true);

  /// process function when runnning over Monte Carlo with MC truth enabled

  /// \param col subscribe to FemtoDreamCollision table
  /// \param parts subscribe to the joined table of FemtoDreamParticles and FemtoDreamMCLabels table
  /// \param FemtoDramMCParticles subscribe to the table containing the Monte Carlo Truth information
  void processMC(o2::aod::FDCollision& col, FemtoFullParticlesMC& parts, o2::aod::FDMCParticles&)
  {
    auto groupPartsOne = partsOneMC->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
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
