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
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "CCDB/BasicCCDBManager.h"
#include "Framework/GroupSlicer.h"
#include "Framework/StepTHn.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/RunningWorkflowInfo.h"
#include "CommonConstants/MathConstants.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/Multiplicity.h"
#include "PWGCF/DataModel/CorrelationsDerived.h"
#include "PWGCF/Core/CorrelationContainer.h"
#include "PWGCF/Core/PairCuts.h"
#include "DataFormatsParameters/GRPObject.h"
#include "PWGHF/DataModel/HFSecondaryVertex.h"
#include "PWGHF/DataModel/HFCandidateSelectionTables.h"
#include "ReconstructionDataFormats/GlobalTrackID.h"

#include <TH1F.h>
#include <cmath>
#include <TDirectory.h>
#include <THn.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace constants::math;
using namespace o2::aod::hf_cand;
using namespace o2::aod::hf_cand_prong2;
using namespace o2::analysis::hf_cuts_d0_topik;

struct HfTaskFlow {
  HistogramRegistry registry{"registry"};

  OutputObj<CorrelationContainer> sameTPCTPCCh{"sameEventTPCTPCChHadrons"};
  OutputObj<CorrelationContainer> sameTPCMFTCh{"sameEventTPCMFTChHadrons"};
  OutputObj<CorrelationContainer> sameHF{"sameEventHFHadrons"};
  OutputObj<CorrelationContainer> mixedTPCTPCCh{"mixedEventTPCTPCChHadrons"};
  OutputObj<CorrelationContainer> mixedHF{"mixedEventHFHadrons"};

  //  configurables for processing options
  Configurable<bool> processRun2{"processRun2", "false", "Flag to run on Run 2 data"};
  Configurable<bool> processRun3{"processRun3", "true", "Flag to run on Run 3 data"};

  Configurable<bool> processTPCTPChh{"processTPCTPChh", "true", "Flag to process TPC-TPC h-h correlations"};
  Configurable<bool> processTPCMFThh{"processTPCMFThh", "true", "Flag to process TPC-MFT h-h correlations"};
  Configurable<bool> processHFHadrons{"processHFHadrons", "false", "Flag to process HF-h correlations"};

  Configurable<int> cfgNoMixedEvents{"cfgNoMixedEvents", 5, "Number of mixed events per event"};

  //  configurables for collisions
  Configurable<float> cfgCutVertex{"cfgCutVertex", 7.0f, "Accepted z-vertex range"};

  //  configurables for associated particles
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8f, "Eta range for tracks"};
  Configurable<float> cfgCutPt{"cfgCutPt", 0.5f, "Minimum pT for tracks"};

  //  configurables for HF candidates
  Configurable<int> d_selectionFlagD0{"d_selectionFlagD0", 1, "Selection Flag for D0"};
  Configurable<int> d_selectionFlagD0bar{"d_selectionFlagD0bar", 1, "Selection Flag for D0bar"};
  Configurable<double> cutYCandMax{"cutYCandMax", -1., "max. cand. rapidity"};
  Configurable<std::vector<double>> bins{"pTBins", std::vector<double>{hf_cuts_d0_topik::pTBins_v}, "pT bin limits"};

  //  configurables for containers
  ConfigurableAxis axisVertex{"axisVertex", {14, -7, 7}, "vertex axis for histograms"};
  ConfigurableAxis axisDeltaPhi{"axisDeltaPhi", {72, -PIHalf, PIHalf * 3}, "delta phi axis for histograms"};
  ConfigurableAxis axisDeltaEta{"axisDeltaEta", {48, -2.4, 2.4}, "delta eta axis for histograms"};
  ConfigurableAxis axisPtTrigger{"axisPtTrigger", {VARIABLE_WIDTH, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0}, "pt trigger axis for histograms"};
  ConfigurableAxis axisPtAssoc{"axisPtAssoc", {VARIABLE_WIDTH, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0}, "pt associated axis for histograms"};
  ConfigurableAxis axisMultiplicity{"axisMultiplicity", {VARIABLE_WIDTH, 0, 5, 10, 20, 30, 40, 50, 100.1}, "multiplicity axis for histograms"};
  ConfigurableAxis axisVertexEfficiency{"axisVertexEfficiency", {10, -10, 10}, "vertex axis for efficiency histograms"};
  ConfigurableAxis axisEtaEfficiency{"axisEtaEfficiency", {20, -1.0, 1.0}, "eta axis for efficiency histograms"};
  ConfigurableAxis axisPtEfficiency{"axisPtEfficiency", {VARIABLE_WIDTH, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0}, "pt axis for efficiency histograms"};
  //  TODO: flow of HF will need to be done vs. invariant mass, in the signal and side-band regions
  //        either 1) add invariant mass axis or 2) define several containers for different inv. mass regions
  //        Note: don't forget to check inv. mass separately for D0 and D0bar candidate
  ConfigurableAxis axisMass{"axisMass", {30, 1.7, 2.0}, "axis of invariant mass of HF candidates"};

  //  Collision filters
  //  FIXME: The filter is applied also on the candidates! Beware!
  Filter collisionVtxZFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  using aodCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::Mults>>;

  //  Charged track filters
  Filter trackFilter = (nabs(aod::track::eta) < cfgCutEta) &&
                       (aod::track::pt > cfgCutPt) &&
                       requireGlobalTrackWoPtEtaInFilter();
  using aodTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TracksDCA, aod::TrackSelection>>;

  //  HF candidate filter
  //  TODO: use Partition instead of filter
  Filter candidateFilter = aod::hf_selcandidate_d0::isSelD0 >= d_selectionFlagD0 || aod::hf_selcandidate_d0::isSelD0bar >= d_selectionFlagD0bar;
  using hfCandidates = soa::Filtered<soa::Join<aod::HfCandProng2, aod::HFSelD0Candidate>>;

  //  =========================
  //      init()
  //  =========================
  void init(o2::framework::InitContext&)
  {
    //  EVENT HISTOGRAMS
    registry.add("hEventCounter", "hEventCounter", {HistType::kTH1F, {{3, 0.5, 3.5}}});
    //  set axes of the event counter histogram
    const int nBins = 3;
    std::string labels[nBins];
    labels[0] = "all";
    labels[1] = "after trigger selection (Run 2)";
    labels[2] = "after Physics selection";
    for (int iBin = 0; iBin < nBins; iBin++) {
      registry.get<TH1>(HIST("hEventCounter"))->GetXaxis()->SetBinLabel(iBin + 1, labels[iBin].data());
    }
    registry.add("hMultiplicity", "hMultiplicity", {HistType::kTH1F, {{500, 0, 500}}});
    registry.add("hVtxZ", "hVtxZ", {HistType::kTH1F, {{400, -50, 50}}});
    registry.add("hNtracks", "hNtracks", {HistType::kTH1F, {{500, 0, 500}}});

    //  histograms for event mixing
    const int maxMixBin = axisMultiplicity->size() * 14; // 14 bins for z-vertex
    registry.add("hEventCountMixing", "bin", {HistType::kTH1F, {{maxMixBin + 2, -2.5, -0.5 + maxMixBin, "bin"}}});
    registry.add("hEventCountHFMixing", "bin", {HistType::kTH1F, {{maxMixBin + 2, -2.5, -0.5 + maxMixBin, "bin"}}});
    registry.add("hEventCountSame", "bin", {HistType::kTH1F, {{maxMixBin + 2, -2.5, -0.5 + maxMixBin, "bin"}}});
    registry.add("hMultiplicityMixing", "hMultiplicityMixing", {HistType::kTH1F, {{500, 0, 500}}});
    registry.add("hVtxZMixing", "hVtxZMixing", {HistType::kTH1F, {{100, -10, 10}}});
    registry.add("hNtracksMixing", "hNtracksMixing", {HistType::kTH1F, {{500, 0, 500}}});
    registry.add("hMultiplicityHFMixing", "hMultiplicityHFMixing", {HistType::kTH1F, {{500, 0, 500}}});
    registry.add("hVtxZHFMixing", "hVtxZHFMixing", {HistType::kTH1F, {{100, -10, 10}}});
    registry.add("hNtracksHFMixing", "hNtracksHFMixing", {HistType::kTH1F, {{500, 0, 500}}});

    //  TRACK HISTOGRAMS
    //  histograms for associated particles
    registry.add("hYields", "multiplicity vs pT vs eta", {HistType::kTH3F, {{200, 0, 200, "multiplicity"}, {40, 0, 20, "p_{T}"}, {100, -2, 2, "#eta"}}});
    registry.add("hEtaPhi", "multiplicity vs eta vs phi", {HistType::kTH3F, {{200, 0, 200, "multiplicity"}, {100, -2, 2, "#eta"}, {200, 0, 2 * PI, "#varphi"}}});
    registry.add("hPt", "pT", {HistType::kTH1F, {{100, 0, 10, "p_{T}"}}});
    registry.add("hEta", "eta", {HistType::kTH1F, {{100, -4, 4, "#eta"}}});
    registry.add("hPhi", "phi", {HistType::kTH1F, {{100, 0, 2 * PI, "#varphi"}}});

    //  histograms for particles in event mixing
    registry.add("hPtMixing", "pT", {HistType::kTH1F, {{100, 0, 10, "p_{T}"}}});
    registry.add("hEtaMixing", "eta", {HistType::kTH1F, {{100, -4, 4, "#eta"}}});
    registry.add("hPhiMixing", "phi", {HistType::kTH1F, {{100, 0, 2 * PI, "#varphi"}}});

    //  histograms for MFT tracks
    registry.add("hEtaPhiMFT", "multiplicity vs eta vs phi in MFT", {HistType::kTH3F, {{200, 0, 200, "multiplicity"}, {100, -2, 2, "#eta"}, {200, 0, 2 * PI, "#varphi"}}});
    registry.add("hEtaMFT", "etaMFT", {HistType::kTH1F, {{100, -4, 4, "#eta"}}});
    registry.add("hPhiMFT", "phiMFT", {HistType::kTH1F, {{100, 0, 2 * PI, "#varphi"}}});

    //  histograms for candidates
    auto vbins = (std::vector<double>)bins;

    registry.add("hPtCand", "2-prong candidates;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0, 10.}}});
    registry.add("hPtProng0", "2-prong candidates;prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0, 10.}}});
    registry.add("hPtProng1", "2-prong candidates;prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0, 10.}}});
    registry.add("hMass", "2-prong candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{500, 0., 5.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLength", "2-prong candidates;decay length (cm);entries", {HistType::kTH2F, {{200, 0., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLengthXY", "2-prong candidates;decay length xy (cm);entries", {HistType::kTH2F, {{200, 0., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong0", "2-prong candidates;prong 0 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong1", "2-prong candidates;prong 1 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0d0", "2-prong candidates;product of DCAxy to prim. vertex (cm^{2});entries", {HistType::kTH2F, {{500, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCTS", "2-prong candidates;cos #it{#theta}* (D^{0});entries", {HistType::kTH2F, {{110, -1.1, 1.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCt", "2-prong candidates;proper lifetime (D^{0}) * #it{c} (cm);entries", {HistType::kTH2F, {{120, -20., 100.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCPA", "2-prong candidates;cosine of pointing angle;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hEtaCand", "2-prong candidates;candidate #it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hSelectionStatus", "2-prong candidates;selection status;entries", {HistType::kTH2F, {{5, -0.5, 4.5}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hImpParErr", "2-prong candidates;impact parameter error (cm);entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLenErr", "2-prong candidates;decay length error (cm);entries", {HistType::kTH2F, {{100, 0., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLenXYErr", "2-prong candidates;decay length xy error (cm);entries", {HistType::kTH2F, {{100, 0., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});

    //  histograms for candidates in event mixing
    registry.add("hPtHFMixing", "pT", {HistType::kTH1F, {{100, 0, 10, "p_{T}"}}});
    registry.add("hEtaHFMixing", "eta", {HistType::kTH1F, {{100, -4, 4, "#eta"}}});
    registry.add("hPhiHFMixing", "phi", {HistType::kTH1F, {{100, 0, 2 * PI, "#varphi"}}});

    //  set axes of the correlation container
    std::vector<AxisSpec> corrAxis = {{axisDeltaEta, "#Delta#eta"},
                                      {axisPtAssoc, "p_{T} (GeV/c)"},
                                      {axisPtTrigger, "p_{T} (GeV/c)"},
                                      {axisMultiplicity, "multiplicity"},
                                      {axisDeltaPhi, "#Delta#varphi (rad)"},
                                      {axisVertex, "z-vtx (cm)"}};
    std::vector<AxisSpec> effAxis = {{axisEtaEfficiency, "#eta"},
                                     {axisPtEfficiency, "p_{T} (GeV/c)"},
                                     {axisVertexEfficiency, "z-vtx (cm)"}};
    std::vector<AxisSpec> userAxis = {{axisMass, "m_{inv} (GeV/c^{2})"}};
    sameTPCTPCCh.setObject(new CorrelationContainer("sameEventTPCTPCChHadrons", "sameEventTPCTPCChHadrons", corrAxis, effAxis, {}));
    sameTPCMFTCh.setObject(new CorrelationContainer("sameEventTPCMFTChHadrons", "sameEventTPCMFTChHadrons", corrAxis, effAxis, {}));
    sameHF.setObject(new CorrelationContainer("sameEventHFHadrons", "sameEventHFHadrons", corrAxis, effAxis, userAxis));
    mixedTPCTPCCh.setObject(new CorrelationContainer("mixedEventTPCTPCChHadrons", "mixedEventTPCTPCChHadrons", corrAxis, effAxis, {}));
    mixedHF.setObject(new CorrelationContainer("mixedEventHFHadrons", "mixedEventHFHadrons", corrAxis, effAxis, userAxis));
  }

  //  ---------------
  //    templates
  //  FIXME: Some collisions are rejected here, what causes (part of) differences with the D0 task
  //  ---------------
  template <typename TCollision>
  bool isCollisionSelected(TCollision collision, bool fillHistograms = false)
  {
    if (processRun2 == true) {
      //  Run 2: trigger selection for data case
      if (fillHistograms)
        registry.fill(HIST("hEventCounter"), 1);
      if (!collision.alias()[kINT7]) {
        return false;
      }
      //  Run 2: further offline selection
      if (fillHistograms)
        registry.fill(HIST("hEventCounter"), 2);
      if (!collision.sel7()) {
        return false;
      }
      if (fillHistograms)
        registry.fill(HIST("hEventCounter"), 3);
    } else {
      //  Run 3: selection
      if (fillHistograms)
        registry.fill(HIST("hEventCounter"), 1);
      if (!collision.sel8()) {
        return false;
      }
      if (fillHistograms)
        registry.fill(HIST("hEventCounter"), 3);
    }
    return true;
  }

  template <typename TTracks>
  void fillQA(float multiplicity, TTracks tracks)
  {
    int Ntracks = 0;
    for (auto& track1 : tracks) {
      Ntracks++;
      registry.fill(HIST("hPt"), track1.pt());
      registry.fill(HIST("hEta"), track1.eta());
      registry.fill(HIST("hPhi"), track1.phi());
      registry.fill(HIST("hYields"), multiplicity, track1.pt(), track1.eta());
      registry.fill(HIST("hEtaPhi"), multiplicity, track1.eta(), track1.phi());
    }
    registry.fill(HIST("hNtracks"), Ntracks);
  }

  template <typename TTracks>
  void fillMixingQA(float multiplicity, float vz, TTracks tracks)
  {
    registry.fill(HIST("hMultiplicityMixing"), multiplicity);
    registry.fill(HIST("hVtxZMixing"), vz);

    int Ntracks = 0;
    for (auto& track1 : tracks) {
      Ntracks++;
      registry.fill(HIST("hPtMixing"), track1.pt());
      registry.fill(HIST("hEtaMixing"), track1.eta());
      registry.fill(HIST("hPhiMixing"), track1.phi());
    }
    registry.fill(HIST("hNtracksMixing"), Ntracks);
  }

  template <typename TTracks>
  void fillHFMixingQA(float multiplicity, float vz, TTracks tracks)
  {
    registry.fill(HIST("hMultiplicityHFMixing"), multiplicity);
    registry.fill(HIST("hVtxZHFMixing"), vz);

    int Ntracks = 0;
    for (auto& track1 : tracks) {
      Ntracks++;
      registry.fill(HIST("hPtHFMixing"), track1.pt());
      registry.fill(HIST("hEtaHFMixing"), track1.eta());
      registry.fill(HIST("hPhiHFMixing"), track1.phi());
    }
    registry.fill(HIST("hNtracksHFMixing"), Ntracks);
  }

  template <typename TTracks>
  void fillMFTQA(float multiplicity, TTracks tracks)
  {
    for (auto& track1 : tracks) {
      registry.fill(HIST("hEtaMFT"), track1.eta());
      float phi = track1.phi();
      o2::math_utils::bringTo02Pi(phi);
      registry.fill(HIST("hPhiMFT"), phi);
      registry.fill(HIST("hEtaPhiMFT"), multiplicity, track1.eta(), phi);
    }
  }

  //  TODO: Check how to put this into a Filter
  template <typename TTrack>
  bool isAcceptedCandidate(TTrack candidate)
  {
    if (!(candidate.hfflag() & 1 << DecayType::D0ToPiK)) {
      return false;
    }
    if (cutYCandMax >= 0. && std::abs(YD0(candidate)) > cutYCandMax) {
      return false;
    }
    return true;
  }

  //  TODO: Note: we do not need all these plots since they are in D0 and Lc task -> remove it after we are sure this works
  template <typename TTracks>
  void fillCandidateQA(TTracks candidates)
  {
    for (auto& candidate : candidates) {
      if (!isAcceptedCandidate(candidate)) {
        continue;
      }

      if (candidate.isSelD0() >= d_selectionFlagD0) {
        registry.fill(HIST("hMass"), InvMassD0(candidate), candidate.pt());
      }
      if (candidate.isSelD0bar() >= d_selectionFlagD0bar) {
        registry.fill(HIST("hMass"), InvMassD0bar(candidate), candidate.pt());
      }

      registry.fill(HIST("hPtCand"), candidate.pt());
      registry.fill(HIST("hPtProng0"), candidate.ptProng0());
      registry.fill(HIST("hPtProng1"), candidate.ptProng1());
      registry.fill(HIST("hDecLength"), candidate.decayLength(), candidate.pt());
      registry.fill(HIST("hDecLengthXY"), candidate.decayLengthXY(), candidate.pt());
      registry.fill(HIST("hd0Prong0"), candidate.impactParameter0(), candidate.pt());
      registry.fill(HIST("hd0Prong1"), candidate.impactParameter1(), candidate.pt());
      registry.fill(HIST("hd0d0"), candidate.impactParameterProduct(), candidate.pt());
      registry.fill(HIST("hCTS"), CosThetaStarD0(candidate), candidate.pt());
      registry.fill(HIST("hCt"), CtD0(candidate), candidate.pt());
      registry.fill(HIST("hCPA"), candidate.cpa(), candidate.pt());
      registry.fill(HIST("hEtaCand"), candidate.eta(), candidate.pt());
      registry.fill(HIST("hSelectionStatus"), candidate.isSelD0() + (candidate.isSelD0bar() * 2), candidate.pt());
      registry.fill(HIST("hImpParErr"), candidate.errorImpactParameter0(), candidate.pt());
      registry.fill(HIST("hImpParErr"), candidate.errorImpactParameter1(), candidate.pt());
      registry.fill(HIST("hDecLenErr"), candidate.errorDecayLength(), candidate.pt());
      registry.fill(HIST("hDecLenXYErr"), candidate.errorDecayLengthXY(), candidate.pt());
    }
  }

  template <typename TTarget, typename TTracksTrig, typename TTracksAssoc>
  void fillCorrelations(TTarget target, TTracksTrig tracks1, TTracksAssoc tracks2, float multiplicity, float posZ)
  {
    auto triggerWeight = 1;
    auto associatedWeight = 1;

    for (auto& track1 : tracks1) {

      float eta1 = track1.eta();
      float pt1 = track1.pt();
      float phi1 = track1.phi();
      o2::math_utils::bringTo02Pi(phi1);

      //  TODO: add getter for NUE trigger efficiency here

      //  calculating inv. mass to be filled into the container below
      //  Note: this is needed only in case of HF-hadron correlations
      bool fillingHFcontainer = false;
      double invmass = 0;
      if constexpr (std::is_same_v<hfCandidates, TTracksTrig>) {
        //  TODO: Check how to put this into a Filter
        if (!isAcceptedCandidate(track1)) {
          continue;
        }
        fillingHFcontainer = true;
        invmass = InvMassD0(track1);
      }

      //  fill single-track distributions
      if (!fillingHFcontainer) {
        target->getTriggerHist()->Fill(CorrelationContainer::kCFStepReconstructed, pt1, multiplicity, posZ, triggerWeight);
      } else {
        target->getTriggerHist()->Fill(CorrelationContainer::kCFStepReconstructed, pt1, multiplicity, posZ, invmass, triggerWeight);
      }

      for (auto& track2 : tracks2) {

        //  case of h-h correlations where the two types of tracks are the same
        //  this avoids autocorrelations and double counting of particle pairs
        if constexpr (std::is_same_v<TTracksAssoc, TTracksTrig>) {
          if (track1.index() <= track2.index()) {
            continue;
          }
        }

        //  TODO: in case of HF-h correlations, remove candidate daughters from the pool of associated hadrons
        //        with which the candidate is being correlated

        float eta2 = track2.eta();
        float pt2 = track2.pt();
        float phi2 = track2.phi();
        o2::math_utils::bringTo02Pi(phi2);

        //  TODO: add getter for NUE associated efficiency here

        //  TODO: add pair cuts on phi*

        float deltaPhi = phi1 - phi2;
        //  set range of delta phi in (-pi/2 , 3/2*pi)
        deltaPhi = RecoDecay::constrainAngle(deltaPhi, -0.5 * PI);

        if (!fillingHFcontainer) {
          //  fill pair correlations
          target->getPairHist()->Fill(CorrelationContainer::kCFStepReconstructed,
                                      eta1 - eta2, pt2, pt1, multiplicity, deltaPhi, posZ,
                                      triggerWeight * associatedWeight);
        } else {
          target->getPairHist()->Fill(CorrelationContainer::kCFStepReconstructed,
                                      eta1 - eta2, pt2, pt1, multiplicity, deltaPhi, posZ, invmass,
                                      triggerWeight * associatedWeight);
        }
      }
    }
  }

  template <typename TTracksTrig, typename TTracksAssoc, typename TLambda>
  void mixCollisions(aodCollisions& collisions, TTracksTrig& tracks1, TTracksAssoc& tracks2, TLambda getPartsSize, OutputObj<CorrelationContainer>& corrContainer)
  {
    using BinningType = FlexibleBinningPolicy<std::tuple<decltype(getPartsSize)>, aod::collision::PosZ, decltype(getPartsSize)>;
    BinningType binningWithTracksSize{{getPartsSize}, {axisVertex, axisMultiplicity}, true};

    auto tracksTuple = std::make_tuple(tracks1, tracks2);
    Pair<aodCollisions, TTracksTrig, TTracksAssoc, BinningType> pair{binningWithTracksSize, cfgNoMixedEvents, -1, collisions, tracksTuple};

    for (auto& [collision1, tracks1, collision2, tracks2] : pair) {

      if (!(isCollisionSelected(collision1, false))) {
        continue;
      }
      if (!(isCollisionSelected(collision2, false))) {
        continue;
      }

      auto binningValues = binningWithTracksSize.getBinningValues(collision1, collisions);
      int bin = binningWithTracksSize.getBin(binningValues);

      const auto multiplicity = tracks1.size();
      const auto vz = collision1.posZ();

      if (processTPCTPChh) {
        registry.fill(HIST("hEventCountMixing"), bin);
        fillMixingQA(multiplicity, vz, tracks1);
      } else if (processHFHadrons) {
        registry.fill(HIST("hEventCountHFMixing"), bin);
        fillHFMixingQA(multiplicity, vz, tracks1);
      }
      corrContainer->fillEvent(multiplicity, CorrelationContainer::kCFStepReconstructed);
      fillCorrelations(corrContainer, tracks1, tracks2, multiplicity, collision1.posZ());
    }
  }

  // =====================================
  //    process same event correlations
  // =====================================
  void processSameRun3(aodCollisions::iterator const& collision,
                       aodTracks const& tracks,
                       aod::MFTTracks const& mfttracks,
                       hfCandidates const& candidates)
  {
    if (!(isCollisionSelected(collision, true))) {
      return;
    }

    const auto multiplicity = tracks.size();
    registry.fill(HIST("hMultiplicity"), multiplicity);
    registry.fill(HIST("hVtxZ"), collision.posZ());

    BinningPolicyBase<2> baseBinning{{axisVertex, axisMultiplicity}, true};
    int bin = baseBinning.getBin(std::make_tuple(collision.posZ(), multiplicity));
    registry.fill(HIST("hEventCountSame"), bin);

    if (processTPCTPChh) {
      sameTPCTPCCh->fillEvent(multiplicity, CorrelationContainer::kCFStepReconstructed);
      fillQA(multiplicity, tracks);
      fillCorrelations(sameTPCTPCCh, tracks, tracks, multiplicity, collision.posZ());
    }

    if (processTPCMFThh) {
      sameTPCMFTCh->fillEvent(multiplicity, CorrelationContainer::kCFStepReconstructed);
      fillMFTQA(multiplicity, mfttracks);
      fillCorrelations(sameTPCMFTCh, tracks, mfttracks, multiplicity, collision.posZ());
    }

    if (processHFHadrons) {
      sameHF->fillEvent(multiplicity, CorrelationContainer::kCFStepReconstructed);
      fillCandidateQA(candidates);
      fillCorrelations(sameHF, candidates, tracks, multiplicity, collision.posZ());
    }
  }
  PROCESS_SWITCH(HfTaskFlow, processSameRun3, "Process same event for Run 3", true);

  // =====================================
  //    process same event correlations
  // =====================================
  void processSameRun2(aodCollisions::iterator const& collision,
                       aodTracks const& tracks,
                       hfCandidates const& candidates)
  {
    if (!(isCollisionSelected(collision, true))) {
      return;
    }

    const auto multiplicity = tracks.size();
    registry.fill(HIST("hMultiplicity"), multiplicity);
    registry.fill(HIST("hVtxZ"), collision.posZ());

    BinningPolicyBase<2> baseBinning{{axisVertex, axisMultiplicity}, true};
    int bin = baseBinning.getBin(std::make_tuple(collision.posZ(), multiplicity));
    registry.fill(HIST("hEventCountSame"), bin);

    if (processTPCTPChh) {
      sameTPCTPCCh->fillEvent(multiplicity, CorrelationContainer::kCFStepReconstructed);
      fillQA(multiplicity, tracks);
      fillCorrelations(sameTPCTPCCh, tracks, tracks, multiplicity, collision.posZ());
    }

    if (processHFHadrons) {
      sameHF->fillEvent(multiplicity, CorrelationContainer::kCFStepReconstructed);
      fillCandidateQA(candidates);
      fillCorrelations(sameHF, candidates, tracks, multiplicity, collision.posZ());
    }
  }
  PROCESS_SWITCH(HfTaskFlow, processSameRun2, "Process same event for Run 2", false);

  // =====================================
  //    process mixed event correlations
  // =====================================
  //  TODO: add also MFT and HFcandidate options->then it will have to be split into Run2/3 because there is no MFT in Run2

  void processMixed(aodCollisions& collisions,
                    aodTracks& tracks,
                    hfCandidates& candidates)
  {
    auto getTracksSize = [&tracks](aodCollisions::iterator const& col) {
      auto associatedTracks = tracks.sliceByCached(o2::aod::track::collisionId, col.globalIndex()); // it's cached, so slicing/grouping happens only once
      auto size = associatedTracks.size();
      return size;
    };
    auto getCandsSize = [&candidates](aodCollisions::iterator const& col) {
      auto associatedCands = candidates.sliceByCached(o2::aod::track::collisionId, col.globalIndex()); // it's cached, so slicing/grouping happens only once
      auto size = associatedCands.size();
      return size;
    };

    if (processTPCTPChh) {
      mixCollisions(collisions, tracks, tracks, getTracksSize, mixedTPCTPCCh);
    }
    if (processHFHadrons) {
      mixCollisions(collisions, candidates, tracks, getCandsSize, mixedHF);
    }
  }
  PROCESS_SWITCH(HfTaskFlow, processMixed, "Process mixed event", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfTaskFlow>(cfgc)};
}
