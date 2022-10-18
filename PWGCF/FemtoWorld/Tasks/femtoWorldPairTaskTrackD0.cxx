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

/// \file femtoWorldPairTaskTrackD0.cxx
/// \brief Tasks that reads the track tables and D0 mesons (from HF group)
/// \used for the pairing and builds pairs of track and D0/D0bar meson
/// \author Katarzyna Gwiździel, WUT Warsaw, kgwizdzi@cern.ch

#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/StepTHn.h"

#include "Framework/AnalysisDataModel.h"
#include "CCDB/BasicCCDBManager.h"
#include "Framework/GroupSlicer.h"
#include "CommonConstants/MathConstants.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/Multiplicity.h"
#include "PWGCF/Core/CorrelationContainer.h"
#include "PWGCF/Core/PairCuts.h"
#include "DataFormatsParameters/GRPObject.h"
#include "PWGHF/DataModel/HFSecondaryVertex.h"
#include "PWGHF/DataModel/HFCandidateSelectionTables.h"
#include "ReconstructionDataFormats/GlobalTrackID.h"

#include "PWGCF/FemtoWorld/DataModel/FemtoWorldDerived.h"
#include "PWGCF/FemtoWorld/Core/FemtoWorldParticleHisto.h"
#include "PWGCF/FemtoWorld/Core/FemtoWorldEventHisto.h"
#include "PWGCF/FemtoWorld/Core/FemtoWorldPairCleaner.h"
#include "PWGCF/FemtoWorld/Core/FemtoWorldContainer.h"
#include "PWGCF/FemtoWorld/Core/FemtoWorldDetaDphiStar.h"
#include "PWGCF/FemtoWorld/Core/FemtoWorldUtils.h"

using namespace o2;
using namespace o2::analysis::femtoWorld;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

using namespace o2::aod::hf_cand;
using namespace o2::aod::hf_cand_prong2;
using namespace o2::aod::hf_correlation_ddbar;
using namespace o2::analysis::hf_cuts_d0_topik;
using namespace o2::constants::math;

namespace
{
static constexpr int nPart = 2;                                                                         // number of particle types (for us it will be proton and phi)
static constexpr int nCuts = 5;                                                                         // number of cuts
static const std::vector<std::string> partNames{"PartOne", "PartTwo"};                                  // names of the rows
static const std::vector<std::string> cutNames{"MaxPt", "PIDthr", "nSigmaTPC", "nSigmaTPCTOF", "MaxP"}; // names of the columns
static const float cutsTable[nPart][nCuts]{                                                             // cuts table [rows = particles][columns = types of cuts]
                                           {1.5f, 1.f, 3.f, 3.f, 100.f},
                                           {1.5f, 1.f, 5.f, 5.f, 100.f}};

static const std::vector<float> kNsigma = {3.5f, 3.f, 2.5f};

} // namespace

struct femtoWorldPairTaskTrackD0 {

  /// Particle selection part

  /// Table for both particles
  Configurable<LabeledArray<float>> cfgCutTable{"cfgCutTable", {cutsTable[0], nPart, nCuts, partNames, cutNames}, "Particle selections"};
  Configurable<int> cfgNspecies{"ccfgNspecies", 4, "Number of particle spieces with PID info"};

  /// Particle 1 (track)
  Configurable<int> ConfPDGCodePartOne{"ConfPDGCodePartOne", 2212, "Particle 1 - PDG code"};                                 // proton (2212)
  Configurable<std::vector<int>> ConfPIDPartOne{"ConfPIDPartOne", std::vector<int>{2}, "Particle 1 - Read from cutCulator"}; // we also need the possibility to specify whether the bit is true/false ->std>>vector<std::pair<int, int>>int>>
  Configurable<float> cfgPtLowPart1{"cfgPtLowPart1", 0.5, "Lower limit for Pt for the first particle"};                      // change according to wrzesa cuts
  Configurable<float> cfgPtHighPart1{"cfgPtHighPart1", 1.5, "Higher limit for Pt for the first particle"};
  Configurable<float> cfgEtaLowPart1{"cfgEtaLowPart1", -0.8, "Lower limit for Eta for the first particle"};
  Configurable<float> cfgEtaHighPart1{"cfgEtaHighPart1", 0.8, "Higher limit for Eta for the first particle"};
  Configurable<float> cfgDcaXYPart1{"cfgDcaXYPart1", 2.4, "Value for DCA_XY for the first particle"};
  Configurable<float> cfgDcaZPart1{"cfgDcaZPart1", 3.2, "Value for DCA_Z for the first particle"};
  Configurable<int> cfgTpcClPart1{"cfgTpcClPart1", 88, "Number of tpc clasters for the first particle"};             // min number of found TPC clusters
  Configurable<int> cfgTpcCrosRoPart1{"cfgTpcCrosRoPart1", 70, "Number of tpc crossed rows for the first particle"}; // min number of crossed rows
  Configurable<float> cfgChi2TpcPart1{"cfgChi2TpcPart1", 4.0, "Chi2 / cluster for the TPC track segment for the first particle"};
  Configurable<float> cfgChi2ItsPart1{"cfgChi2ItsPart1", 36.0, "Chi2 / cluster for the ITS track segment for the first particle"};

  /// Partition for particle 1 (proton)
  Partition<aod::FemtoWorldParticles> partsOne = (aod::femtoworldparticle::partType == uint8_t(aod::femtoworldparticle::ParticleType::kTrack)) && (aod::femtoworldparticle::pt < cfgPtHighPart1) && (aod::femtoworldparticle::pt > cfgPtLowPart1) // simple pT cuts
                                                 && (aod::femtoworldparticle::eta < cfgEtaHighPart1) && (aod::femtoworldparticle::eta > cfgEtaLowPart1)                                                                                           // Eta cuts
                                                 && (o2::aod::track::dcaXY < cfgDcaXYPart1) && (o2::aod::track::dcaZ < cfgDcaZPart1)                                                                                                              // DCA cuts for XY and Z
                                                 && (aod::femtoworldparticle::tpcNClsFound > (uint8_t)cfgTpcClPart1)                                                                                                                              // Number of found TPC clusters
                                                 && (aod::femtoworldparticle::tpcNClsCrossedRows > (uint8_t)cfgTpcCrosRoPart1)                                                                                                                    // Crossed rows TPC
                                                 && (aod::femtoworldparticle::itsChi2NCl < cfgChi2ItsPart1) && (aod::femtoworldparticle::tpcChi2NCl < cfgChi2TpcPart1)                                                                            //&& // chi2 cuts
                                                 && (aod::femtoworldparticle::sign > int8_t(0));

  /// Histogramming for particle 1
  FemtoWorldParticleHisto<aod::femtoworldparticle::ParticleType::kTrack, 0> trackHistoPartOne;

  // Histogram for D0/D0bar candidates
  HistogramRegistry hfhfHistos{"hfhfHistos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // configurables for HF candidates
  Configurable<int> d_selectionFlagD0{"d_selectionFlagD0", 1, "Selection Flag for D0"};
  Configurable<int> d_selectionFlagD0bar{"d_selectionFlagD0bar", 1, "Selection Flag for D0bar"};
  Configurable<double> cutYCandMax{"cutYCandMax", -1., "max. cand. rapidity"};
  Configurable<std::vector<double>> bins{"pTBins", std::vector<double>{hf_cuts_d0_topik::pTBins_v}, "pT bin limits"};

  //  HF candidate filter
  Filter candidateFilter = aod::hf_selcandidate_d0::isSelD0 >= d_selectionFlagD0 || aod::hf_selcandidate_d0::isSelD0bar >= d_selectionFlagD0bar;
  using hfCandidates = soa::Filtered<soa::Join<aod::HfCandProng2, aod::HFSelD0Candidate>>;
  using hfCandidate = hfCandidates::iterator;
  // HF candidate partition
  // Partition<soa::Join<aod::HfCandProng2, aod::HFSelD0Candidate>> selectedD0Candidates = aod::hf_selcandidate_d0::isSelD0 >= d_selectionFlagD0 || aod::hf_selcandidate_d0::isSelD0bar >= d_selectionFlagD0bar;

  /// Histogramming for Event
  FemtoWorldEventHisto eventHisto;

  /// The configurables need to be passed to an std::vector
  std::vector<int> vPIDPartOne;

  /// Correlation part
  ConfigurableAxis CfgMultBins{"CfgMultBins", {VARIABLE_WIDTH, 0.0f, 20.0f, 40.0f, 60.0f, 80.0f, 100.0f, 200.0f, 99999.f}, "Mixing bins - multiplicity"};
  ConfigurableAxis CfgVtxBins{"CfgVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
  ConfigurableAxis CfgkstarBins{"CfgkstarBins", {1500, 0., 6.}, "binning kstar"};
  ConfigurableAxis CfgkTBins{"CfgkTBins", {150, 0., 9.}, "binning kT"};
  ConfigurableAxis CfgmTBins{"CfgmTBins", {225, 0., 7.5}, "binning mT"};
  Configurable<int> ConfPhiBins{"ConfPhiBins", 15, "Number of phi bins in deta dphi"};
  Configurable<int> ConfEtaBins{"ConfEtaBins", 15, "Number of eta bins in deta dphi"};
  Configurable<int> ConfMInvBins{"ConfMInvBins", 1000, "Number of bins in mInv distribution"};
  Configurable<int> ConfNEventsMix{"ConfNEventsMix", 5, "Number of events for mixing"};
  Configurable<bool> ConfIsCPR{"ConfIsCPR", true, "Close Pair Rejection"};
  Configurable<bool> ConfCPRPlotPerRadii{"ConfCPRPlotPerRadii", false, "Plot CPR per radii"};

  FemtoWorldContainer<femtoWorldContainer::EventType::same, femtoWorldContainer::Observable::kstar> sameEventCont;
  FemtoWorldContainer<femtoWorldContainer::EventType::mixed, femtoWorldContainer::Observable::kstar> mixedEventCont;
  FemtoWorldPairCleaner<aod::femtoworldparticle::ParticleType::kTrack, aod::femtoworldparticle::ParticleType::kPhi> pairCleaner;
  FemtoWorldDetaDphiStar<aod::femtoworldparticle::ParticleType::kTrack, aod::femtoworldparticle::ParticleType::kPhi> pairCloseRejection;
  /// Histogram output
  HistogramRegistry qaRegistry{"TrackQA", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry resultRegistry{"Correlations", {}, OutputObjHandlingPolicy::AnalysisObject};

  // PID for protons
  bool IsProtonNSigma(float mom, float nsigmaTPCP, float nsigmaTOFP) // from: https://github.com/alisw/AliPhysics/blob/master/PWGCF/FEMTOSCOPY/AliFemtoUser/AliFemtoMJTrackCut.cxx
  {
    bool fNsigmaTPCTOF = true;
    bool fNsigmaTPConly = true;
    double fNsigma = 3.0;
    double fNsigma2 = 3.0;

    if (fNsigmaTPCTOF) {
      if (mom > 0.5) {
        //        if (TMath::Hypot( nsigmaTOFP, nsigmaTPCP )/TMath::Sqrt(2) < 3.0)
        if (mom < 2.0) {
          if (TMath::Hypot(nsigmaTOFP, nsigmaTPCP) < fNsigma)
            return true;
        } else if (TMath::Hypot(nsigmaTOFP, nsigmaTPCP) < fNsigma2)
          return true;
      } else {
        if (TMath::Abs(nsigmaTPCP) < fNsigma)
          return true;
      }
    } else if (fNsigmaTPConly) {
      if (TMath::Abs(nsigmaTPCP) < fNsigma)
        return true;
    } else {
      if (mom > 0.8 && mom < 2.5) {
        if (TMath::Abs(nsigmaTPCP) < 3.0 && TMath::Abs(nsigmaTOFP) < 3.0)
          return true;
      } else if (mom > 2.5) {
        if (TMath::Abs(nsigmaTPCP) < 3.0 && TMath::Abs(nsigmaTOFP) < 2.0)
          return true;
      } else {
        if (TMath::Abs(nsigmaTPCP) < 3.0)
          return true;
      }
    }

    return false;
  }

  void init(o2::framework::InitContext&)
  {
    eventHisto.init(&qaRegistry);
    trackHistoPartOne.init(&qaRegistry);

    sameEventCont.init(&resultRegistry, CfgkstarBins, CfgMultBins, CfgkTBins, CfgmTBins, ConfPhiBins, ConfEtaBins, ConfMInvBins);
    sameEventCont.setPDGCodes(ConfPDGCodePartOne, ConfPDGCodePartTwo);
    mixedEventCont.init(&resultRegistry, CfgkstarBins, CfgMultBins, CfgkTBins, CfgmTBins, ConfPhiBins, ConfEtaBins, ConfMInvBins);
    mixedEventCont.setPDGCodes(ConfPDGCodePartOne, ConfPDGCodePartTwo);
    pairCleaner.init(&qaRegistry);
    if (ConfIsCPR) {
      pairCloseRejection.init(&resultRegistry, &qaRegistry, 0.01, 0.01, ConfCPRPlotPerRadii); /// \todo add config for Δη and ΔΦ cut values
    }

    vPIDPartOne = ConfPIDPartOne;

    // Define your axes
    // Constant bin width axis
    AxisSpec vtxZAxis = {100, -20, 20};
    // Variable bin width axis
    std::vector<double> ptBinning = {0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.8, 2.0, 2.2, 2.4, 2.8, 3.2, 3.6, 4.};
    AxisSpec ptAxis = {ptBinning, "#it{p}_{T} (GeV/#it{c})"};

    hfHistos.add("hZvtx", ";Z (cm)", kTH1F, {vtxZAxis});
    hfHistos.add("hP", ";#it{p} (GeV/#it{c})", kTH1F, {{35, 0.5, 4.}});
    hfHistos.add("hEta", ";#it{p} (GeV/#it{c})", kTH1F, {{100, -1.5, 1.5}});
    hfHistos.add("hPt", ";#it{p}_{T} (GeV/#it{c})", kTH1F, {ptAxis});
    hfHistos.add("hPhi", ";#varphi (rad)", kTH1F, {{10, -TMath::TwoPi(), TMath::TwoPi()}});
    hfHistos.add("hDeltaPhi", ";#Delta #varphi (rad)", kTH1F, {{10, -TMath::TwoPi(), TMath::TwoPi()}});
  }

  /// The function process takes a track and D0 meson
  void process(o2::aod::FemtoWorldCollision const& coll,
               o2::aod::FemtoWorldParticles const& parts,
               hfCandidate const& d0Cans)
  {

    hfHistos.fill(HIST("hZvtx"), coll.posZ());
    hfHistos.fill(HIST("hP"), d0Cans.p());
    hfHistos.fill(HIST("hPt"), d0Cans.pt());
    hfHistos.fill(HIST("hEta"), d0Cans.eta());
    hfHistos.fill(HIST("hPhi"), d0Cans.phi());

    for (auto part : parts) {
      for (auto d0meson : d0Cans) {
        hfHistos.fill(HIST("hDeltaPhi", part.phi() - d0meson.phi()));
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<femtoWorldPairTaskTrackD0>(cfgc),
  };
  return workflow;
}
