// Copyright 2019-2025 CERN and copyright holders of ALICE O2.
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
  Configurable<aod::femtodreamparticle::cutContainerType> ConfTrk1_TPCBit_Reject{"ConfTrk1_TPCBit_Reject", 0, "PID TPC bit from cutCulator to reject a particle hypothesis for particle 1 (set to 0 to ignore)"};
  Configurable<float> ConfTrk1_minPt{"ConfTrk1_minPt", 0., "Minimum pT of partricle 1 (Track)"};
  Configurable<float> ConfTrk1_maxPt{"ConfTrk1_maxPt", 999., "Maximum pT of partricle 1 (Track)"};
  Configurable<float> ConfTrk1_minEta{"ConfTrk1_minEta", -10., "Minimum eta of partricle 1 (Track)"};
  Configurable<float> ConfTrk1_maxEta{"ConfTrk1_maxEta", 10., "Maximum eta of partricle 1 (Track)"};
  Configurable<float> ConfTrk1_TempFitVarMin{"ConfTrk1_TempFitVarMin", -10., "Minimum DCAxy of partricle 1 (Track)"};
  Configurable<float> ConfTrk1_TempFitVarMax{"ConfTrk1_TempFitVarMax", 10., "Maximum DCAxy of partricle 1 (Track)"};
  Configurable<float> ConfTrk1_PIDThres{"ConfTrk1_PIDThres", 0.75, "Particle 1 - Read from cutCulator"};

  Configurable<bool> ConfOptDCACutPtDep{"ConfOptDCACutPtDep", false, "Use pt dependent dca cut"};
  Configurable<bool> ConfUseRun2Function{"ConfUseRun2Function", true, "Use Run2 pT dependent DCA selection function"};
  Configurable<bool> ConfOptCorrelatedPlots{"ConfOptCorrelatedPlots", false, "Enable additional three dimensional histogramms. High memory consumption. Use for debugging"};
  ConfigurableAxis ConfBinmult{"ConfBinmult", {1, 0, 1}, "multiplicity Binning"};
  ConfigurableAxis ConfBinmultPercentile{"ConfBinmultPercentile", {10, 0.0f, 100.0f}, "multiplicity percentile Binning"};
  ConfigurableAxis ConfBinpT{"ConfBinpT", {{240, 0, 6}}, "pT binning"};
  ConfigurableAxis ConfBineta{"ConfBineta", {{200, -1.5, 1.5}}, "eta binning"};
  ConfigurableAxis ConfBinphi{"ConfBinphi", {{200, 0, TMath::TwoPi()}}, "phi binning"};

  ConfigurableAxis ConfTempFitVarBins{"ConfTempFitVarBins", {300, -0.15, 0.15}, "Binning of the TempFitVar"};
  ConfigurableAxis ConfNsigmaTPCBins{"ConfNsigmaTPCBins", {1600, -8, 8}, "Binning of Nsigma TPC plot"};
  ConfigurableAxis ConfNsigmaTOFBins{"ConfNsigmaTOFBins", {3000, -15, 15}, "Binning of the Nsigma TOF plot"};
  ConfigurableAxis ConfNsigmaTPCTOFBins{"ConfNsigmaTPCTOFBins", {3000, -15, 15}, "Binning of the Nsigma TPC+TOF plot"};
  ConfigurableAxis ConfNsigmaITSBins{"ConfNsigmaITSBins", {3000, -15, 15}, "Binning of the Nsigma ITS plot"};
  ConfigurableAxis ConfTPCclustersBins{"ConfTPCClustersBins", {163, -0.5, 162.5}, "Binning of TPC found clusters plot"};
  Configurable<int> ConfTempFitVarMomentum{"ConfTempFitVarMomentum", 0, "Momentum used for binning: 0 -> pt; 1 -> preco; 2 -> ptpc"};
  ConfigurableAxis ConfDummy{"ConfDummy", {1, 0, 1}, "Dummy axis for inv mass"};

  Configurable<bool> ConfdoCentCut{"ConfdoCentCut", false, "Enable centrality cut"};
  Configurable<float> ConfCentMax{"ConfCentMax", 100., "Upper limit of centrality cut"};
  Configurable<float> ConfCentMin{"ConfCentMin", 0., "Lower limit of centrality cut"};

  using FemtoMCCollisions = Join<aod::FDCollisions, aod::FDMCCollLabels>;
  using FemtoMCCollision = FemtoMCCollisions::iterator;

  using FemtoFullParticles = soa::Join<aod::FDParticles, aod::FDExtParticles>;

  Partition<FemtoFullParticles> partsOne = (aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kTrack)) &&
                                           (ncheckbit(aod::femtodreamparticle::cut, ConfTrk1_CutBit)) &&
                                           ifnode(aod::femtodreamparticle::pt * (nexp(aod::femtodreamparticle::eta) + nexp(-1.f * aod::femtodreamparticle::eta)) / 2.f <= ConfTrk1_PIDThres, ncheckbit(aod::femtodreamparticle::pidcut, ConfTrk1_TPCBit) && ((aod::femtodreamparticle::pidcut & ConfTrk1_TPCBit_Reject) == 0u), ncheckbit(aod::femtodreamparticle::pidcut, ConfTrk1_TPCTOFBit)) &&
                                           (aod::femtodreamparticle::pt > ConfTrk1_minPt) &&
                                           (aod::femtodreamparticle::pt < ConfTrk1_maxPt) &&
                                           (aod::femtodreamparticle::eta > ConfTrk1_minEta) &&
                                           (aod::femtodreamparticle::eta < ConfTrk1_maxEta) &&
                                           ifnode(ConfOptDCACutPtDep, ifnode(ConfUseRun2Function, nabs(aod::femtodreamparticle::tempFitVar) < 0.0105f + (0.035f / npow(aod::femtodreamparticle::pt, 1.1f)), nabs(aod::femtodreamparticle::tempFitVar) < 0.004f + (0.013f / aod::femtodreamparticle::pt)),
                                                  (aod::femtodreamparticle::tempFitVar > ConfTrk1_TempFitVarMin) &&
                                                    (aod::femtodreamparticle::tempFitVar < ConfTrk1_TempFitVarMax));

  Preslice<FemtoFullParticles> perColReco = aod::femtodreamparticle::fdCollisionId;
  // adsf
  using FemtoFullParticlesMC = soa::Join<aod::FDParticles, aod::FDExtParticles, aod::FDMCLabels, aod::FDExtMCLabels>;
  Partition<FemtoFullParticlesMC> partsOneMC =
    (aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kTrack)) &&
    ncheckbit(aod::femtodreamparticle::cut, ConfTrk1_CutBit) &&
    ifnode(aod::femtodreamparticle::pt * (nexp(aod::femtodreamparticle::eta) + nexp(-1.f * aod::femtodreamparticle::eta)) / 2.f <= ConfTrk1_PIDThres, ncheckbit(aod::femtodreamparticle::pidcut, ConfTrk1_TPCBit) && ((aod::femtodreamparticle::pidcut & ConfTrk1_TPCBit_Reject) == 0u), ncheckbit(aod::femtodreamparticle::pidcut, ConfTrk1_TPCTOFBit));
  Preslice<FemtoFullParticlesMC> perColGen = aod::femtodreamparticle::fdCollisionId;

  /// Histogramming for Event
  FemtoDreamEventHisto eventHisto;

  FemtoDreamParticleHisto<aod::femtodreamparticle::ParticleType::kTrack> trackHisto;

  /// Histogram output
  HistogramRegistry qaRegistry{"TrackQA", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext&)
  {
    eventHisto.init(&qaRegistry, ConfIsMC);
    trackHisto.init(&qaRegistry, ConfBinmult, ConfBinmultPercentile, ConfBinpT, ConfBineta, ConfBinphi, ConfTempFitVarBins, ConfNsigmaTPCBins, ConfNsigmaTOFBins, ConfNsigmaTPCTOFBins, ConfNsigmaITSBins, ConfDummy, ConfDummy, ConfIsMC, ConfTrk1_PDGCode.value, true, ConfOptCorrelatedPlots);
  }

  /// Porduce QA plots for sigle track selection in FemtoDream framework
  template <bool isMC, typename CollisionType, typename PartitionType>
  void FillDebugHistos(CollisionType& col, PartitionType& groupPartsOne)
  {
    eventHisto.fillQA<isMC>(col);

    if (ConfdoCentCut.value && (col.multV0M() > ConfCentMax || col.multV0M() < ConfCentMin)) {
      return;
    }

    for (auto& part : groupPartsOne) {
      trackHisto.fillQA<isMC, true>(part, static_cast<aod::femtodreamparticle::MomentumType>(ConfTempFitVarMomentum.value), col.multNtr(), col.multV0M(), ConfOptCorrelatedPlots);
    }
  }

  /// process function when runnning over data/ Monte Carlo reconstructed only
  /// \param col subscribe to FemtoDreamCollision table
  /// \param parts subscribe to FemtoDreamParticles table
  void processData(o2::aod::FDCollision& col, FemtoFullParticles&)
  {
    auto groupPartsOne = partsOne->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
    FillDebugHistos<false>(col, groupPartsOne);
  }
  PROCESS_SWITCH(femtoDreamDebugTrack, processData, "Enable Debug processing for Monte Carlo", true);

  /// process function when runnning over Monte Carlo with MC truth enabled

  /// \param col subscribe to FemtoDreamCollision table
  /// \param parts subscribe to the joined table of FemtoDreamParticles and FemtoDreamMCLabels table
  /// \param FemtoDramMCParticles subscribe to the table containing the Monte Carlo Truth information
  void processMC(FemtoMCCollision& col, o2::aod::FDMCCollisions&, FemtoFullParticlesMC& /*parts*/, aod::FDMCParticles&, aod::FDExtMCParticles&)
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
