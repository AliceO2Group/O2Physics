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

/// \file femtoDreamDebugTrack.cxx
/// \brief Tasks that reads the particle tables and fills QA histograms for tracks
/// \author Luca Barioglio, TU München, luca.barioglio@cern.ch

#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/StepTHn.h"
#include "DataFormatsParameters/GRPObject.h"

#include "PWGCF/DataModel/FemtoDerived.h"
#include "FemtoDreamParticleHisto.h"
#include "FemtoDreamEventHisto.h"
#include "FemtoUtils.h"

using namespace o2;
using namespace o2::analysis::femtoDream;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

namespace
{
static constexpr int nCuts = 5;
static const std::vector<std::string> cutNames{"MaxPt", "PIDthr", "nSigmaTPC", "nSigmaTPCTOF", "MaxP"};
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
  Configurable<std::vector<float>> ConfPIDnSigmaMax{"ConfPIDnSigmaMax", std::vector<float>{3.5f, 3.f, 2.5f}, "This configurable needs to be the same as the one used in the producer task"};

  using FemtoFullParticles = soa::Join<aod::FemtoDreamParticles, aod::FemtoDreamDebugParticles>;
  Partition<FemtoFullParticles> partsOne = (aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kTrack)) && ((aod::femtodreamparticle::cut & ConfCutPartOne) == ConfCutPartOne);

  using FemtoFullParticlesMC = soa::Join<aod::FemtoDreamParticles, aod::FemtoDreamDebugParticles, aod::FemtoDreamMCLabels>;
  Partition<FemtoFullParticlesMC> partsOneMC = (aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kTrack)) && ((aod::femtodreamparticle::cut & ConfCutPartOne) == ConfCutPartOne);

  /// Histogramming for Event
  FemtoDreamEventHisto eventHisto;

  /// The configurables need to be passed to an std::vector
  std::vector<int> vPIDPartOne;
  std::vector<float> kNsigma;
  /// Histogram output
  HistogramRegistry qaRegistry{"TrackQA", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry FullQaRegistry{"FullTrackQA", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry FullQaRegistryMC{"FullTrackQA_MC", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext&)
  {
    eventHisto.init(&qaRegistry);

    FullQaRegistry.add("FullTrackQA/hPt", "; #it{p}_{T} (GeV/#it{c}); Entries", kTH1F, {{240, 0, 6}});
    FullQaRegistry.add("FullTrackQA/hEta", "; #eta; Entries", kTH1F, {{200, -1.5, 1.5}});
    FullQaRegistry.add("FullTrackQA/hPhi", "; #phi; Entries", kTH1F, {{200, 0, 2. * M_PI}});
    FullQaRegistry.add("FullTrackQA/hTPCfindable", "; TPC findable clusters; Entries", kTH1F, {{163, -0.5, 162.5}});
    FullQaRegistry.add("FullTrackQA/hTPCfound", "; TPC found clusters; Entries", kTH1F, {{163, -0.5, 162.5}});
    FullQaRegistry.add("FullTrackQA/hTPCcrossedOverFindalbe", "; TPC ratio findable; Entries", kTH1F, {{100, 0.5, 1.5}});
    FullQaRegistry.add("FullTrackQA/hTPCcrossedRows", "; TPC crossed rows; Entries", kTH1F, {{163, -0.5, 162.5}});
    FullQaRegistry.add("FullTrackQA/hTPCfindableVsCrossed", ";TPC findable clusters ; TPC crossed rows;", kTH2F, {{163, -0.5, 162.5}, {163, -0.5, 162.5}});
    FullQaRegistry.add("FullTrackQA/hTPCshared", "; TPC shared clusters; Entries", kTH1F, {{163, -0.5, 162.5}});
    FullQaRegistry.add("FullTrackQA/hITSclusters", "; ITS clusters; Entries", kTH1F, {{10, -0.5, 9.5}});
    FullQaRegistry.add("FullTrackQA/hITSclustersIB", "; ITS clusters in IB; Entries", kTH1F, {{10, -0.5, 9.5}});
    FullQaRegistry.add("FullTrackQA/hDCAxy", "; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", kTH2F, {{20, 0.5, 4.05}, {500, -5, 5}});
    FullQaRegistry.add("FullTrackQA/hDCAz", "; #it{p}_{T} (GeV/#it{c}); DCA_{z} (cm)", kTH2F, {{100, 0, 10}, {500, -5, 5}});
    FullQaRegistry.add("FullTrackQA/hDCA", "; #it{p}_{T} (GeV/#it{c}); DCA (cm)", kTH2F, {{100, 0, 10}, {301, 0., 1.5}});
    FullQaRegistry.add("FullTrackQA/hTPCdEdX", "; #it{p} (GeV/#it{c}); TPC Signal", kTH2F, {{100, 0, 10}, {1000, 0, 1000}});
    FullQaRegistry.add("FullTrackQA/nSigmaTPC_el", "; #it{p} (GeV/#it{c}); n#sigma_{TPC}^{e}", kTH2F, {{100, 0, 10}, {100, -5, 5}});
    FullQaRegistry.add("FullTrackQA/nSigmaTPC_pi", "; #it{p} (GeV/#it{c}); n#sigma_{TPC}^{#pi}", kTH2F, {{100, 0, 10}, {100, -5, 5}});
    FullQaRegistry.add("FullTrackQA/nSigmaTPC_K", "; #it{p} (GeV/#it{c}); n#sigma_{TPC}^{K}", kTH2F, {{100, 0, 10}, {100, -5, 5}});
    FullQaRegistry.add("FullTrackQA/nSigmaTPC_p", "; #it{p} (GeV/#it{c}); n#sigma_{TPC}^{p}", kTH2F, {{100, 0, 10}, {100, -5, 5}});
    FullQaRegistry.add("FullTrackQA/nSigmaTPC_d", "; #it{p} (GeV/#it{c}); n#sigma_{TPC}^{d}", kTH2F, {{100, 0, 10}, {100, -5, 5}});
    FullQaRegistry.add("FullTrackQA/nSigmaTOF_el", "; #it{p} (GeV/#it{c}); n#sigma_{TOF}^{e}", kTH2F, {{100, 0, 10}, {100, -5, 5}});
    FullQaRegistry.add("FullTrackQA/nSigmaTOF_pi", "; #it{p} (GeV/#it{c}); n#sigma_{TOF}^{#pi}", kTH2F, {{100, 0, 10}, {100, -5, 5}});
    FullQaRegistry.add("FullTrackQA/nSigmaTOF_K", "; #it{p} (GeV/#it{c}); n#sigma_{TOF}^{K}", kTH2F, {{100, 0, 10}, {100, -5, 5}});
    FullQaRegistry.add("FullTrackQA/nSigmaTOF_p", "; #it{p} (GeV/#it{c}); n#sigma_{TOF}^{p}", kTH2F, {{100, 0, 10}, {100, -5, 5}});
    FullQaRegistry.add("FullTrackQA/nSigmaTOF_d", "; #it{p} (GeV/#it{c}); n#sigma_{TOF}^{d}", kTH2F, {{100, 0, 10}, {100, -5, 5}});
    FullQaRegistry.add("FullTrackQA/nSigmaComb_el", "; #it{p} (GeV/#it{c}); n#sigma_{comb}^{e}", kTH2F, {{100, 0, 10}, {100, -5, 5}});
    FullQaRegistry.add("FullTrackQA/nSigmaComb_pi", "; #it{p} (GeV/#it{c}); n#sigma_{comb}^{#pi}", kTH2F, {{100, 0, 10}, {100, -5, 5}});
    FullQaRegistry.add("FullTrackQA/nSigmaComb_K", "; #it{p} (GeV/#it{c}); n#sigma_{comb}^{K}", kTH2F, {{100, 0, 10}, {100, -5, 5}});
    FullQaRegistry.add("FullTrackQA/nSigmaComb_p", "; #it{p} (GeV/#it{c}); n#sigma_{comb}^{p}", kTH2F, {{100, 0, 10}, {100, -5, 5}});
    FullQaRegistry.add("FullTrackQA/nSigmaComb_d", "; #it{p} (GeV/#it{c}); n#sigma_{comb}^{d}", kTH2F, {{100, 0, 10}, {100, -5, 5}});

    if (ConfIsMC) {
      FullQaRegistry.add("FullTrackQA_MC/hPt_MC", "; #it{p}_{T} (GeV/#it{c}); Entries", kTH1F, {{240, 0, 6}});
      FullQaRegistry.add("FullTrackQA_MC/hEta_MC", "; #eta; Entries", kTH1F, {{200, -1.5, 1.5}});
      FullQaRegistry.add("FullTrackQA_MC/hPhi_MC", "; #phi; Entries", kTH1F, {{200, 0, 2. * M_PI}});
      FullQaRegistry.add("FullTrackQA_MC/hPDG", "; PDG; Entries", kTH1I, {{6001, -3000, 3000}});
      FullQaRegistry.add("FullTrackQA_MC/hOrigin_MC", "; Origin; Entries", kTH1I, {{8, 0, 7}});
      FullQaRegistry.add("FullTrackQA_MC/hNoMCtruthCounter", "; Counter; Entries", kTH1I, {{1, 0, 1}});
    }
    vPIDPartOne = ConfPIDPartOne.value;
    kNsigma = ConfPIDnSigmaMax.value;
  }

  /// Porduce QA plots for sigle track selection in FemtoDream framework
  template <bool isMC, typename PartitionType>
  void FillDebugHistos(o2::aod::FemtoDreamCollision& col, PartitionType& groupPartsOne)
  {

    eventHisto.fillQA(col);

    for (auto& part : groupPartsOne) {
      if (part.p() > cfgCutTable->get("MaxP") || part.pt() > cfgCutTable->get("MaxPt")) {
        continue;
      }
      if (!isFullPIDSelected(part.pidcut(),
                             part.p(),
                             cfgCutTable->get("PIDthr"),
                             vPIDPartOne,
                             cfgNspecies,
                             kNsigma,
                             cfgCutTable->get("nSigmaTPC"),
                             cfgCutTable->get("nSigmaTPCTOF"))) {
        continue;
      }

      FullQaRegistry.fill(HIST("FullTrackQA/hPt"), part.pt());
      FullQaRegistry.fill(HIST("FullTrackQA/hEta"), part.eta());
      FullQaRegistry.fill(HIST("FullTrackQA/hPhi"), part.phi());
      FullQaRegistry.fill(HIST("FullTrackQA/hDCAxy"), part.pt(), part.tempFitVar());
      FullQaRegistry.fill(HIST("FullTrackQA/hTPCfindable"), part.tpcNClsFindable());
      FullQaRegistry.fill(HIST("FullTrackQA/hTPCfound"), part.tpcNClsFound());
      FullQaRegistry.fill(HIST("FullTrackQA/hTPCcrossedOverFindalbe"), part.tpcCrossedRowsOverFindableCls());
      FullQaRegistry.fill(HIST("FullTrackQA/hTPCcrossedRows"), part.tpcNClsCrossedRows());
      FullQaRegistry.fill(HIST("FullTrackQA/hTPCfindableVsCrossed"), part.tpcNClsFindable(), part.tpcNClsCrossedRows());
      FullQaRegistry.fill(HIST("FullTrackQA/hTPCshared"), part.tpcNClsShared());
      FullQaRegistry.fill(HIST("FullTrackQA/hITSclusters"), part.itsNCls());
      FullQaRegistry.fill(HIST("FullTrackQA/hITSclustersIB"), part.itsNClsInnerBarrel());
      FullQaRegistry.fill(HIST("FullTrackQA/hDCAz"), part.pt(), part.dcaZ());
      FullQaRegistry.fill(HIST("FullTrackQA/hDCA"), part.pt(), std::sqrt(pow(part.dcaXY(), 2.) + pow(part.dcaZ(), 2.)));
      FullQaRegistry.fill(HIST("FullTrackQA/hTPCdEdX"), part.p(), part.tpcSignal());
      FullQaRegistry.fill(HIST("FullTrackQA/nSigmaTPC_el"), part.p(), part.tpcNSigmaEl());
      FullQaRegistry.fill(HIST("FullTrackQA/nSigmaTPC_pi"), part.p(), part.tpcNSigmaPi());
      FullQaRegistry.fill(HIST("FullTrackQA/nSigmaTPC_K"), part.p(), part.tpcNSigmaKa());
      FullQaRegistry.fill(HIST("FullTrackQA/nSigmaTPC_p"), part.p(), part.tpcNSigmaPr());
      FullQaRegistry.fill(HIST("FullTrackQA/nSigmaTPC_d"), part.p(), part.tpcNSigmaDe());
      FullQaRegistry.fill(HIST("FullTrackQA/nSigmaTOF_el"), part.p(), part.tofNSigmaEl());
      FullQaRegistry.fill(HIST("FullTrackQA/nSigmaTOF_pi"), part.p(), part.tofNSigmaPi());
      FullQaRegistry.fill(HIST("FullTrackQA/nSigmaTOF_K"), part.p(), part.tofNSigmaKa());
      FullQaRegistry.fill(HIST("FullTrackQA/nSigmaTOF_p"), part.p(), part.tofNSigmaPr());
      FullQaRegistry.fill(HIST("FullTrackQA/nSigmaTOF_d"), part.p(), part.tofNSigmaDe());
      FullQaRegistry.fill(HIST("FullTrackQA/nSigmaComb_el"), part.p(), std::sqrt(part.tpcNSigmaEl() * part.tpcNSigmaEl() + part.tofNSigmaEl() * part.tofNSigmaEl()));
      FullQaRegistry.fill(HIST("FullTrackQA/nSigmaComb_pi"), part.p(), std::sqrt(part.tpcNSigmaPi() * part.tpcNSigmaPi() + part.tofNSigmaPi() * part.tofNSigmaPi()));
      FullQaRegistry.fill(HIST("FullTrackQA/nSigmaComb_K"), part.p(), std::sqrt(part.tpcNSigmaKa() * part.tpcNSigmaKa() + part.tofNSigmaKa() * part.tofNSigmaKa()));
      FullQaRegistry.fill(HIST("FullTrackQA/nSigmaComb_p"), part.p(), std::sqrt(part.tpcNSigmaPr() * part.tpcNSigmaPr() + part.tofNSigmaPr() * part.tofNSigmaPr()));
      FullQaRegistry.fill(HIST("FullTrackQA/nSigmaComb_d"), part.p(), std::sqrt(part.tpcNSigmaDe() * part.tpcNSigmaDe() + part.tofNSigmaDe() * part.tofNSigmaDe()));

      if constexpr (isMC) {
        if (part.has_femtoDreamMCParticle()) {
          auto partMC = part.template femtoDreamMCParticle_as<o2::aod::FemtoDreamMCParticles>();
          FullQaRegistry.fill(HIST("FullTrackQA_MC/hPt_MC"), partMC.pt());
          FullQaRegistry.fill(HIST("FullTrackQA_MC/hEta_MC"), partMC.eta());
          FullQaRegistry.fill(HIST("FullTrackQA_MC/hPhi_MC"), partMC.phi());
          FullQaRegistry.fill(HIST("FullTrackQA_MC/hPDG"), partMC.pdgMCTruth());
          FullQaRegistry.fill(HIST("FullTrackQA_MC/hOrigin_MC"), partMC.partOriginMCTruth());
        } else {
          FullQaRegistry.fill(HIST("FullTrackQA_MC/hNoMCtruthCounter"), 1);
        }
      }
    }
  }

  /// Porduce QA plots for sigle track selection in FemtoDream framework
  ///

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
  /// @param FemtoDramMCParticles subscribe to the table containing the Monte Carlo Truth information
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
