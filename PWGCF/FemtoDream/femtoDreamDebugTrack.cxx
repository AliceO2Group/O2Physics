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

/// \file femtoDreamPairTaskTrackTrack.cxx
/// \brief Tasks that reads the track tables used for the pairing and builds pairs of two tracks
/// \author Andi Mathis, TU München, andreas.mathis@ph.tum.de

#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/StepTHn.h"
#include <CCDB/BasicCCDBManager.h>
#include "DataFormatsParameters/GRPObject.h"

#include "PWGCF/DataModel/FemtoDerived.h"
#include "FemtoDreamParticleHisto.h"
#include "FemtoDreamEventHisto.h"

using namespace o2;
using namespace o2::analysis::femtoDream;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

namespace
{
static constexpr int nCuts = 4;
static const std::vector<std::string> cutNames{"MaxPt", "PIDthr", "nSigmaTPC", "nSigmaTPCTOF"};
static const float cutsTable[1][nCuts] = {{4.05f, 0.75f, 3.5f, 3.5f}};

static constexpr int nNsigma = 3;
static constexpr float kNsigma[nNsigma] = {3.5f, 3.f, 2.5f};

enum kPIDselection {
  k3d5sigma = 0,
  k3sigma = 1,
  k2d5sigma = 2
};

enum kDetector {
  kTPC = 0,
  kTPCTOF = 1,
  kNdetectors = 2
};
} // namespace

struct femtoDreamDebugTrack {

  Configurable<LabeledArray<float>> cfgCutTable{"cfgCutTable", {cutsTable[0], nCuts, cutNames}, "Particle selections"};
  Configurable<int> cfgNspecies{"ccfgNspecies", 4, "Number of particle spieces with PID info"};

  Configurable<int> ConfPDGCodePartOne{"ConfPDGCodePartOne", 2212, "Particle 1 - PDG code"};
  Configurable<uint32_t> ConfCutPartOne{"ConfCutPartOne", 5542474, "Particle 1 - Selection bit from cutCulator"};
  Configurable<std::vector<int>> ConfPIDPartOne{"ConfPIDPartOne", std::vector<int>{2}, "Particle 1 - Read from cutCulator"};

  using FemtoFullTracks = soa::Join<aod::FemtoDreamParticles, aod::FemtoDreamDebugParticles>;

  Partition<FemtoFullTracks> partsOne = (aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kTrack)) &&
                                        (aod::femtodreamparticle::pt < cfgCutTable->get("MaxPt")) &&
                                        ((aod::femtodreamparticle::cut & ConfCutPartOne) == ConfCutPartOne);

  /// Histogramming for Event
  FemtoDreamEventHisto eventHisto;

  /// The configurables need to be passed to an std::vector
  std::vector<int> vPIDPartOne;

  /// Histogram output
  HistogramRegistry qaRegistry{"TrackQA", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry FullQaRegistry{"FullTrackQA", {}, OutputObjHandlingPolicy::AnalysisObject};

  Service<o2::ccdb::BasicCCDBManager> ccdb; /// Accessing the CCDB

  void init(InitContext&)
  {
    eventHisto.init(&qaRegistry);

    FullQaRegistry.add("FullTrackQA/hPt", "; #it{p}_{T} (GeV/#it{c}); Entries", kTH1F, {{240, 0, 6}});
    FullQaRegistry.add("FullTrackQA/hEta", "; #eta; Entries", kTH1F, {{200, -1.5, 1.5}});
    FullQaRegistry.add("FullTrackQA/hPhi", "; #phi; Entries", kTH1F, {{200, 0, 2. * M_PI}});
    FullQaRegistry.add("FullTrackQA/hTPCfindable", "; TPC findable clusters; Entries", kTH1F, {{163, 0, 163}});
    FullQaRegistry.add("FullTrackQA/hTPCfound", "; TPC found clusters; Entries", kTH1F, {{163, 0, 163}});
    FullQaRegistry.add("FullTrackQA/hTPCcrossedOverFindalbe", "; TPC ratio findable; Entries", kTH1F, {{100, 0.5, 1.5}});
    FullQaRegistry.add("FullTrackQA/hTPCcrossedRows", "; TPC crossed rows; Entries", kTH1F, {{163, 0, 163}});
    FullQaRegistry.add("FullTrackQA/hTPCfindableVsCrossed", ";TPC findable clusters ; TPC crossed rows;", kTH2F, {{163, 0, 163}, {163, 0, 163}});
    FullQaRegistry.add("FullTrackQA/hTPCshared", "; TPC shared clusters; Entries", kTH1F, {{163, 0, 163}});
    FullQaRegistry.add("FullTrackQA/hDCAxy", "; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", kTH2F, {{20, 0.5, 4.05}, {500, -5, 5}});
    FullQaRegistry.add("FullTrackQA/hDCAz", "; #it{p}_{T} (GeV/#it{c}); DCA_{z} (cm)", kTH2F, {{100, 0, 10}, {500, -5, 5}});
    FullQaRegistry.add("FullTrackQA/hDCA", "; #it{p}_{T} (GeV/#it{c}); DCA (cm)", kTH2F, {{100, 0, 10}, {301, 0., 1.5}});
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

    vPIDPartOne = ConfPIDPartOne;

    /// Initializing CCDB
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    long now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    ccdb->setCreatedNotAfter(now);
  }

  /// Function to retrieve the nominal mgnetic field in kG (0.1T) and convert it directly to T
  float getMagneticFieldTesla(uint64_t timestamp)
  {
    // TODO done only once (and not per run). Will be replaced by CCDBConfigurable
    static o2::parameters::GRPObject* grpo = nullptr;
    if (grpo == nullptr) {
      grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>("GLO/GRP/GRP", timestamp);
      if (grpo == nullptr) {
        LOGF(fatal, "GRP object not found for timestamp %llu", timestamp);
        return 0;
      }
      LOGF(info, "Retrieved GRP for timestamp %llu with magnetic field of %d kG", timestamp, grpo->getNominalL3Field());
    }
    float output = 0.1 * (grpo->getNominalL3Field());
    return output;
  }

  /// internal function that returns the kPIDselection element corresponding to a specifica n-sigma value
  /// \param number of sigmas for PIF
  /// \return kPIDselection corresponing to n-sigma
  kPIDselection getPIDselection(float nSigma)
  {
    for (int i = 0; i < nNsigma; i++) {
      if (abs(nSigma - kNsigma[i]) < 1e-3) {
        return static_cast<kPIDselection>(i);
      }
    }
    LOG(info) << "Invalid value of nSigma: " << nSigma << ". Standard 3 sigma returned." << std::endl;
    return kPIDselection::k3sigma;
  }

  /// function that checks whether the PID selection specified in the vectors is fulfilled
  /// \param pidcut Bit-wise container for the PID
  /// \param vSpecies Vector with the different selections
  /// \return Whether the PID selection specified in the vectors is fulfilled
  bool isPIDSelected(aod::femtodreamparticle::cutContainerType const& pidcut, std::vector<int> const& vSpecies, float nSigma, kDetector iDet = kDetector::kTPC)
  {
    bool pidSelection = true;
    kPIDselection kNsigma = getPIDselection(nSigma);
    for (auto iSpecies : vSpecies) {
      //\todo we also need the possibility to specify whether the bit is true/false ->std>>vector<std::pair<int, int>>
      // if (!((pidcut >> it.first) & it.second)) {
      int bit_to_check = cfgNspecies * kDetector::kNdetectors * kNsigma + iSpecies * kDetector::kNdetectors + iDet;
      if (!(pidcut & (1UL << bit_to_check))) {
        pidSelection = false;
      }
    }
    return pidSelection;
  };

  /// function that checks whether the PID selection specified in the vectors is fulfilled, depending on the momentum TPC or TPC+TOF PID is conducted
  /// \param pidcut Bit-wise container for the PID
  /// \param mom Momentum of the track
  /// \param pidThresh Momentum threshold that separates between TPC and TPC+TOF PID
  /// \param vecTPC Vector with the different selections for the TPC PID
  /// \param vecComb Vector with the different selections for the TPC+TOF PID
  /// \return Whether the PID selection is fulfilled
  bool isFullPIDSelected(aod::femtodreamparticle::cutContainerType const& pidCut, float const& momentum, float const& pidThresh, std::vector<int> const& vSpecies, float nSigmaTPC = 3.5, float nSigmaTPCTOF = 3.5)
  {
    bool pidSelection = true;
    if (momentum < pidThresh) {
      /// TPC PID only
      pidSelection = isPIDSelected(pidCut, vSpecies, nSigmaTPC, kDetector::kTPC);
    } else {
      /// TPC + TOF PID
      pidSelection = isPIDSelected(pidCut, vSpecies, nSigmaTPCTOF, kDetector::kTPCTOF);
    }
    return pidSelection;
  };

  /// Porduce QA plots for sigle track selection in FemtoDream framework
  void process(o2::aod::FemtoDreamCollision& col, FemtoFullTracks& parts)
  {
    auto groupPartsOne = partsOne->sliceByCached(aod::femtodreamparticle::femtoDreamCollisionId, col.globalIndex());

    eventHisto.fillQA(col);

    for (auto& part : groupPartsOne) {

      if (!isFullPIDSelected(part.pidcut(), part.p(), cfgCutTable->get("PIDthr"), vPIDPartOne, cfgCutTable->get("nSigmaTPC"), cfgCutTable->get("nSigmaTPCTOF"))) {
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
      FullQaRegistry.fill(HIST("FullTrackQA/hDCAz"), part.pt(), part.dcaZ());
      FullQaRegistry.fill(HIST("FullTrackQA/hDCA"), part.pt(), std::sqrt(pow(part.dcaXY(), 2.) + pow(part.dcaZ(), 2.)));
      FullQaRegistry.fill(HIST("FullTrackQA/nSigmaTPC_el"), part.tpcInnerParam(), part.tpcNSigmaEl());
      FullQaRegistry.fill(HIST("FullTrackQA/nSigmaTPC_pi"), part.tpcInnerParam(), part.tpcNSigmaPi());
      FullQaRegistry.fill(HIST("FullTrackQA/nSigmaTPC_K"), part.tpcInnerParam(), part.tpcNSigmaKa());
      FullQaRegistry.fill(HIST("FullTrackQA/nSigmaTPC_p"), part.tpcInnerParam(), part.tpcNSigmaPr());
      FullQaRegistry.fill(HIST("FullTrackQA/nSigmaTPC_d"), part.tpcInnerParam(), part.tpcNSigmaDe());
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
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<femtoDreamDebugTrack>(cfgc),
  };
  return workflow;
}
