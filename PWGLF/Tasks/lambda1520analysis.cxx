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

/// \file lambda1520analysis.cxx
//
/// \brief 1. This task loops over all tracks and creates the corresponding analysis
/// tables that contain the typical information of Kaons and Protons required for analysis.
/// 2. It reconstructs track-track decay lambda(1520) resonance candidate
//
/// \author Hirak Kumar Koley <hirak.kumar.koley@cern.ch>

#include <CCDB/BasicCCDBManager.h>
#include <TLorentzVector.h>
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/runDataProcessing.h"
#include "PWGLF/DataModel/LFResonanceTables.h"
#include "DataFormatsParameters/GRPObject.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace std;

enum EventSelection { kNoSelection = 0,
                      kVertexCut = 1 };

namespace o2::aod
{
namespace tracktag
{
// Global bool
DECLARE_SOA_COLUMN(IsInteresting, isInteresting, bool); //! will this be built or not?
DECLARE_SOA_COLUMN(IsKaon, isKaon, bool);               //! will this be built or not?
DECLARE_SOA_COLUMN(IsProton, isProton, bool);           //! will this be built or not?
} // namespace tracktag
DECLARE_SOA_TABLE(TrackTags, "AOD", "TRACKTAGS",
                  tracktag::IsInteresting,
                  tracktag::IsKaon,
                  tracktag::IsProton);
} // namespace o2::aod

// Track preselector for lambda1520 study
struct TrackPreselector {
  // PID selection configurables
  Configurable<float> pidnSigmaPreSelectionCut{"pidnSigmaPreSelectionCut", 3.0f, "TPC and TOF minimal nSigma PID cut"};
  Configurable<float> pidnSigmaTPCCuthasTOF{"pidnSigmaTPCCuthasTOF", 5.0f, "nSigma cut for TPC if track has TOF hit"};
  Configurable<float> maxptcutTPCKaon{"maxptcutTPCKaon", 0.6f, "max pt cut for tpc only Kaons"};
  Configurable<float> maxptcutTPCProton{"maxptcutTPCProton", 1.1f, "max pt cut for tpc only Protons"};
  Configurable<float> minPtcutTOFKaon{"minPtcutTOFKaon", 0.1f, "min pt cut for tof Kaons"};
  Configurable<float> minPtcutTOFProton{"minPtcutTOFProton", 0.1f, "min pt cut for tof Protons"};

  Configurable<int> nBinsNSigma{"nBinsNSigma", 400, "Number of nSigma bins"};

  Produces<aod::TrackTags> tracktags; // Track tags

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(o2::framework::InitContext&)
  {
    // axes
    AxisSpec axisPt{100, 0, 10, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec axisEta{500, -1.0f, 1.0f, "#it{#eta}"};
    AxisSpec axisPhi{100, 0.f, 6.3f, "#it{#phi}"};
    AxisSpec axisDCA{1000, -5, 5};
    AxisSpec axisTPCcrossedrow{200, 0, 200};   // resolution = 1
    AxisSpec axisNSigma{nBinsNSigma, -10, 10}; // resolution = 0.05

    //  Track QA
    //  --- Track
    histos.add("QA/dcaZ", "DCA_{Z} distribution of Tracks; #it{p}_{T} (GeV/#it{c}); DCA_{Z} (cm); ", {HistType::kTH2F, {axisPt, axisDCA}});
    histos.add("QA/dcaXY", "DCA_{XY} momentum distribution of Tracks; #it{p}_{T} (GeV/#it{c}); DCA_{XY} (cm);", {HistType::kTH2F, {axisPt, axisDCA}});
    histos.add("QA/TPC_CR", "# TPC Crossedrows distribution of Tracks; #it{p}_{T} (GeV/#it{c}); TPC Crossed rows", {HistType::kTH2F, {axisPt, axisTPCcrossedrow}});
    histos.add("QA/pT", "pT distribution of Tracks; #it{p}_{T} (GeV/#it{c}); Counts;", {HistType::kTH1F, {axisPt}});
    histos.add("QA/eta", "#eta distribution of Tracks; #eta; Counts;", {HistType::kTH1F, {axisEta}});
    histos.add("QA/phi", "#phi distribution of Tracks; #phi; Counts;", {HistType::kTH1F, {axisPhi}});
    histos.add("QA/TOF_TPC_Map_Kaon", "Kaon TOF + TPC Combined PID of Tracks; #sigma_{TOF}^{Kaon}; #sigma_{TPC}^{Kaon}; Counts;", {HistType::kTH2F, {axisNSigma, axisNSigma}});
    histos.add("QA/TOF_Nsigma_Kaon", "Kaon TOF NSigma of Tracks; #it{p}_{T} (GeV/#it{c}); #sigma_{TOF}^{Kaon};", {HistType::kTH2F, {axisPt, axisNSigma}});
    histos.add("QA/TPC_Nsigma_Kaon", "Kaon TPC NSigma of Tracks; #it{p}_{T} (GeV/#it{c}); #sigma_{TPC}^{Kaon};", {HistType::kTH2F, {axisPt, axisNSigma}});
    histos.add("QA/TOF_TPC_Map_Proton", "Proton TOF + TPC Combined PID of Tracks; #sigma_{TOF}^{Proton}; #sigma_{TPC}^{Proton}; Counts;", {HistType::kTH2F, {axisNSigma, axisNSigma}});
    histos.add("QA/TOF_Nsigma_Proton", "Proton TOF NSigma of Tracks; #it{p}_{T} (GeV/#it{c}); #sigma_{TOF}^{Proton};", {HistType::kTH2F, {axisPt, axisNSigma}});
    histos.add("QA/TPC_Nsigma_Proton", "Proton TPC NSigma of Tracks; #it{p}_{T} (GeV/#it{c}); #sigma_{TPC}^{Proton};", {HistType::kTH2F, {axisPt, axisNSigma}});
  }

  // Kaon selection -> y/n
  template <bool fillHistos, typename kCurrentTrackType>
  bool IsKaonSelected(kCurrentTrackType const& kCurrentTrack)
  {
    // TPC only Tracks
    if ((kCurrentTrack.tofPIDselectionFlag() & aod::resodaughter::kHasTOF) != aod::resodaughter::kHasTOF) { // TPC track only check
      if (kCurrentTrack.pt() < maxptcutTPCKaon && std::abs(kCurrentTrack.tpcNSigmaKa()) > pidnSigmaPreSelectionCut)
        return false;
      if (kCurrentTrack.pt() > maxptcutTPCKaon)
        return false;
    }

    // TPC+TOF PID cut
    if ((kCurrentTrack.tofPIDselectionFlag() & aod::resodaughter::kHasTOF) == aod::resodaughter::kHasTOF) { // TOF track check
      if (kCurrentTrack.pt() >= minPtcutTOFKaon && std::abs(kCurrentTrack.tofNSigmaKa()) > static_cast<float_t>(pidnSigmaPreSelectionCut))
        return false;
      if (kCurrentTrack.pt() >= minPtcutTOFKaon && std::abs(kCurrentTrack.tpcNSigmaKa()) > static_cast<float_t>(pidnSigmaTPCCuthasTOF))
        return false;
    }
    if constexpr (fillHistos) {
      //  --- PID QA Kaons
    }
    return true;
  }

  // Proton selection -> y/n
  template <bool fillHistos, typename kCurrentTrackType>
  bool IsProtonSelected(kCurrentTrackType const& kCurrentTrack)
  {
    // TPC PID cut
    if ((kCurrentTrack.tofPIDselectionFlag() & aod::resodaughter::kHasTOF) != aod::resodaughter::kHasTOF) {
      if (kCurrentTrack.pt() < maxptcutTPCProton && std::abs(kCurrentTrack.tpcNSigmaPr()) > pidnSigmaPreSelectionCut)
        return false;
      if (kCurrentTrack.pt() > maxptcutTPCProton)
        return false;
    }

    // TOF PID cut
    if ((kCurrentTrack.tofPIDselectionFlag() & aod::resodaughter::kHasTOF) == aod::resodaughter::kHasTOF) {
      if (kCurrentTrack.pt() >= minPtcutTOFProton && std::abs(kCurrentTrack.tofNSigmaPr()) > static_cast<float_t>(pidnSigmaPreSelectionCut))
        return false;
      if (kCurrentTrack.pt() >= minPtcutTOFProton && std::abs(kCurrentTrack.tpcNSigmaPr()) > static_cast<float_t>(pidnSigmaTPCCuthasTOF))
        return false;
    }

    if constexpr (fillHistos) {
      //  --- PID QA Protons
    }
    return true;
  }

  template <bool fillHistos, typename kCurrentTrackType>
  bool IsTrackSelected(kCurrentTrackType const& kCurrentTrack)
  {
    //  --- Track Cuts (for future use)
    // For now don't apply any cuts here, otherwise there will be
    // table joining error between aod::ResoTracks & aod::TrackTags

    if constexpr (fillHistos) {
      //  --- Tracks QA
      histos.fill(HIST("QA/eta"), kCurrentTrack.eta());
      histos.fill(HIST("QA/phi"), kCurrentTrack.phi());
      histos.fill(HIST("QA/pT"), kCurrentTrack.pt());
      histos.fill(HIST("QA/dcaZ"), kCurrentTrack.pt(), kCurrentTrack.dcaZ());
      histos.fill(HIST("QA/dcaXY"), kCurrentTrack.pt(), kCurrentTrack.dcaXY());
      histos.fill(HIST("QA/TPC_CR"), kCurrentTrack.pt(), kCurrentTrack.tpcNClsCrossedRows());
      histos.fill(HIST("QA/TOF_TPC_Map_Kaon"), kCurrentTrack.tofNSigmaKa(), kCurrentTrack.tpcNSigmaKa());
      histos.fill(HIST("QA/TOF_Nsigma_Kaon"), kCurrentTrack.pt(), kCurrentTrack.tofNSigmaKa());
      histos.fill(HIST("QA/TPC_Nsigma_Kaon"), kCurrentTrack.pt(), kCurrentTrack.tpcNSigmaKa());
      histos.fill(HIST("QA/TOF_TPC_Map_Proton"), kCurrentTrack.tofNSigmaPr(), kCurrentTrack.tpcNSigmaPr());
      histos.fill(HIST("QA/TOF_Nsigma_Proton"), kCurrentTrack.pt(), kCurrentTrack.tofNSigmaPr());
      histos.fill(HIST("QA/TPC_Nsigma_Proton"), kCurrentTrack.pt(), kCurrentTrack.tpcNSigmaPr());
    }
    return true;
  }

  // Checking PID
  template <typename TrackType>
  void checkPID(TrackType const& kCurrentTrack, bool& lIsInteresting, bool& lIsKaon, bool& lIsProton)
  {
    // LF PID Check
    if (IsKaonSelected<false>(kCurrentTrack)) {
      lIsKaon = 1;
      lIsInteresting = 1;
    }
    if (IsProtonSelected<false>(kCurrentTrack)) {
      lIsProton = 1;
      lIsInteresting = 1;
    }
  }

  void processBuild(aod::ResoCollisions& collisions, o2::aod::ResoTracks const& resotracks)
  {
    for (auto& kCurrentTrack : resotracks) {
      bool lIsInteresting = false;
      bool lIsKaon = false;
      bool lIsProton = false;

      if (!IsTrackSelected<true>(kCurrentTrack))
        continue;
      checkPID(kCurrentTrack, lIsInteresting, lIsKaon, lIsProton);
      tracktags(lIsInteresting, lIsKaon, lIsProton);
    }
  }
  PROCESS_SWITCH(TrackPreselector, processBuild, "Switch to build tagged Tracks with PID preselection", true);
};

// Pre-selected Tracks
using TaggedTracks = soa::Join<aod::ResoTracks, aod::TrackTags>;
using TaggedMCTracks = soa::Join<aod::ResoTracks, aod::TrackTags, aod::ResoMCTracks>;

struct lambda1520analysis {
  // Define slice per Resocollision
  SliceCache cache;
  Preslice<aod::ResoTracks> perRCol = aod::resodaughter::resoCollisionId;

  framework::Service<o2::ccdb::BasicCCDBManager> ccdb; /// Accessing the CCDB
  ConfigurableAxis CfgMultBins{"CfgMultBins", {VARIABLE_WIDTH, 0.0f, 20.0f, 40.0f, 60.0f, 80.0f, 100.0f, 200.0f, 99999.f}, "Mixing bins - multiplicity"};
  ConfigurableAxis CfgVtxBins{"CfgVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // Configurables
  // Event cuts
  Configurable<float> ConfEvtZvtx{"ConfEvtZvtx", 10.f, "Evt sel: Max. z-Vertex (cm)"};

  // Pre-selection Track cuts
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8f, "Eta range for tracks"};
  Configurable<float> cMinPtcut{"cMinPtcut", 0.15f, "Minimal pT for tracks"};
  Configurable<int> cMinTPCncr{"cMinTPCncr", 70, "Minimum number of TPC crossed rows"};

  // DCA Selections
  // DCAr to PV
  Configurable<double> cMaxDCArToPVcut{"cMaxDCArToPVcut", 0.1f, "Track DCAr cut to PV Maximum"};
  // DCAz to PV
  Configurable<double> cMaxDCAzToPVcut{"cMaxDCAzToPVcut", 2.0f, "Track DCAz cut to PV Maximum"};
  Configurable<double> cMinDCAzToPVcut{"cMinDCAzToPVcut", 0.0f, "Track DCAz cut to PV Minimum"};

  Configurable<int> cfgNoMixedEvents{"cfgNoMixedEvents", 3, "Number of mixed events per event"};

  Configurable<int> nBinsNSigma{"nBinsNSigma", 400, "Number of nSigma bins"};
  Configurable<int> nBinsMult{"nBinsMult", 200, "Number of mass bins (Safe limit <= 800)"};

  Filter taggedFilter = aod::tracktag::isInteresting == true;

  void init(o2::framework::InitContext&)
  {
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    uint64_t now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    ccdb->setCreatedNotAfter(now);

    // axes
    AxisSpec axisvtxZ{100, -20, 20};
    AxisSpec axisPt{100, 0, 10, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec axisMassLambda1520{500, 1.3, 1.8, "Invariant Mass (GeV/#it{c}^2)"};
    AxisSpec axisNSigma{nBinsNSigma, -10, 10};
    AxisSpec axisMult{nBinsMult, 0.0f, 1000.0f, "mult_{TPC}"};
    AxisSpec axisDCA{1000, -5, 5};
    AxisSpec axisTPCcrossedrow{200, 0, 200};

    //  Event QA
    auto h = histos.add<TH1>("QA/Event/EnumEvents", "Event selection; ; Event Count", {HistType::kTH1F, {{10, 0, 10}}});
    h->GetXaxis()->SetBinLabel(1, "Events read and Ev. sel. passed");
    h->GetXaxis()->SetBinLabel(2, "posZ passed");
    histos.add("QA/Event/VertexZ", "Event selection; Vertex Z (cm); # Events", {HistType::kTH1F, {axisvtxZ}});
    histos.add("QA/Event/Mult", "Multiplicity Distribution; # Mult; # Events", {HistType::kTH1F, {axisMult}});
    histos.add("QA/Event/Sphericity", "Sphericity Distribution; S_{0}; # Events", {HistType::kTH1F, {{340, -0.2, 3.2}}});

    //  PID QA
    //  --- Kaon
    histos.add("QA/Kaon/TOF_TPC_Map", "TOF + TPC Combined PID for Kaons; #sigma_{TOF}^{Kaon}; #sigma_{TPC}^{Kaon}; Counts;", {HistType::kTH2F, {axisNSigma, axisNSigma}});
    histos.add("QA/Kaon/TOF_Nsigma", "TOF NSigma for Kaons; #it{p}_{T} (GeV/#it{c}); #sigma_{TOF}^{Kaon};", {HistType::kTH2F, {axisPt, axisNSigma}});
    histos.add("QA/Kaon/TPC_Nsigma", "TPC NSigma for Kaons; #it{p}_{T} (GeV/#it{c}); #sigma_{TPC}^{Kaon};", {HistType::kTH2F, {axisPt, axisNSigma}});
    histos.add("QA/Kaon/dcaZ", "DCA_{Z} distribution of selected Kaons; #it{p}_{T} (GeV/#it{c}); DCA_{Z} (cm); ", kTH2F, {axisPt, axisDCA});
    histos.add("QA/Kaon/dcaXY", "DCA_{XY} momentum distribution of selected Kaons; #it{p}_{T} (GeV/#it{c}); DCA_{XY} (cm);", kTH2F, {axisPt, axisDCA});
    histos.add("QA/Kaon/TPC_CR", "# TPC Crossedrows distribution of selected Kaons; #it{p}_{T} (GeV/#it{c}); TPC Crossed rows", kTH2F, {axisPt, axisTPCcrossedrow});
    histos.add("QA/Kaon/pT", "pT distribution of Kaons; #it{p}_{T} (GeV/#it{c}); Counts;", {HistType::kTH1F, {axisPt}});

    //  --- Proton
    histos.add("QA/Proton/TOF_TPC_Map", "TOF + TPC Combined PID for Protons; #sigma_{TOF}^{Proton}; #sigma_{TPC}^{Proton}; Counts;", {HistType::kTH2F, {axisNSigma, axisNSigma}});
    histos.add("QA/Proton/TOF_Nsigma", "TOF NSigma for Protons; #it{p}_{T} (GeV/#it{c}); #sigma_{TOF}^{Proton};", {HistType::kTH2F, {axisPt, axisNSigma}});
    histos.add("QA/Proton/TPC_Nsigma", "TPC NSigma for Protons; #it{p}_{T} (GeV/#it{c}); #sigma_{TPC}^{Proton};", {HistType::kTH2F, {axisPt, axisNSigma}});
    histos.add("QA/Proton/dcaZ", "DCA_{Z} distribution of selected Protons; #it{p}_{T} (GeV/#it{c}); DCA_{Z} (cm);", kTH2F, {axisPt, axisDCA});
    histos.add("QA/Proton/dcaXY", "DCA_{XY} momentum distribution of selected Protons; #it{p}_{T} (GeV/#it{c}); DCA_{XY} (cm);", kTH2F, {axisPt, axisDCA});
    histos.add("QA/Proton/TPC_CR", "# TPC Crossedrows distribution of selected Protons; #it{p}_{T} (GeV/#it{c}); TPC Crossed rows", kTH2F, {axisPt, axisTPCcrossedrow});
    histos.add("QA/Proton/pT", "pT distribution of Protons; #it{p}_{T} (GeV/#it{c}); Counts;", {HistType::kTH1F, {axisPt}});

    // Mass QA (quick check)
    histos.add("Analysis/lambda1520invmass", "Invariant mass of #Lambda(1520) K^{#pm}p^{#mp}; Invariant Mass (GeV/#it{c}^2); Counts;", {HistType::kTH1F, {axisMassLambda1520}});
    histos.add("Analysis/lambda1520invmassLS", "Invariant mass of #Lambda(1520) Like Sign Method K^{#pm}p^{#pm}; Invariant Mass (GeV/#it{c}^2); Counts;", {HistType::kTH1F, {axisMassLambda1520}});
    histos.add("Analysis/lambda1520invmassLSkp", "Invariant mass of #Lambda(1520) Like Sign Method K^{#plus}p^{#plus}; Invariant Mass (GeV/#it{c}^2); Counts;", {HistType::kTH1F, {axisMassLambda1520}});         // K+ + Pr
    histos.add("Analysis/lambda1520invmassLSkbarpbar", "Invariant mass of #Lambda(1520) Like Sign Method K^{#minus}p^{#minus}; Invariant Mass (GeV/#it{c}^2); Counts;", {HistType::kTH1F, {axisMassLambda1520}}); // K- + anti-Pr
    histos.add("Analysis/lambda1520invmassME", "Invariant mass of #Lambda(1520) mixed event K^{#pm}p^{#mp}; Invariant Mass (GeV/#it{c}^2); Counts;", {HistType::kTH1F, {axisMassLambda1520}});

    // 3d histogram
    histos.add("Analysis/h3lambda1520invmass", "Invariant mass of #Lambda(1520) K^{#pm}p^{#mp}", HistType::kTH3F, {axisMult, axisPt, axisMassLambda1520});
    histos.add("Analysis/h3lambda1520invmassLS", "Invariant mass of #Lambda(1520) Like Sign Method K^{#pm}p^{#pm}", HistType::kTH3F, {axisMult, axisPt, axisMassLambda1520});
    histos.add("Analysis/h3lambda1520invmassLSkp", "Invariant mass of #Lambda(1520) Like Sign Method K^{#plus}p^{#plus}", HistType::kTH3F, {axisMult, axisPt, axisMassLambda1520});         // K+ + Pr
    histos.add("Analysis/h3lambda1520invmassLSkbarpbar", "Invariant mass of #Lambda(1520) Like Sign Method K^{#minus}p^{#minus}", HistType::kTH3F, {axisMult, axisPt, axisMassLambda1520}); // K- + anti-Pr
    histos.add("Analysis/h3lambda1520invmassME", "Invariant mass of #Lambda(1520) mixed event K^{#pm}p^{#mp}", HistType::kTH3F, {axisMult, axisPt, axisMassLambda1520});

    if (doprocessMC) {
      histos.add("MC/h3recolambda1520invmass", "Invariant mass of Reconstructed MC #Lambda(1520)", HistType::kTH3F, {axisMult, axisPt, axisMassLambda1520});
      histos.add("MC/truelambda1520pt", "pT distribution of True MC #Lambda(1520); #it{p}_{T} (GeV/#it{c}); Counts;", {HistType::kTH1F, {axisPt}});
      histos.add("MC/recolambda1520pt", "pT distribution of Reconstructed MC #Lambda(1520); #it{p}_{T} (GeV/#it{c}); Counts;", {HistType::kTH1F, {axisPt}});
      histos.add("MC/recolambda1520invmass", "Inv mass distribution of Reconstructed MC #Lambda(1520); #it{p}_{T} (GeV/#it{c}); Counts;", {HistType::kTH1F, {axisMassLambda1520}});
    }
  }

  double massKa = TDatabasePDG::Instance()->GetParticle(kKMinus)->Mass();
  double massPr = TDatabasePDG::Instance()->GetParticle(kProton)->Mass();

  template <bool IsMC, typename CollisionType, typename TracksType>
  void fillHistograms(const CollisionType& collision, const TracksType& dTracks)
  {
    //  Collision QA
    histos.fill(HIST("QA/Event/EnumEvents"), EventSelection::kNoSelection);

    //  --- Event selection: Vertex position
    if (!(fabs(collision.posZ()) < ConfEvtZvtx))
      return;
    histos.fill(HIST("QA/Event/EnumEvents"), EventSelection::kVertexCut);
    histos.fill(HIST("QA/Event/VertexZ"), collision.posZ());
    histos.fill(HIST("QA/Event/Mult"), collision.multTPCtemp());
    histos.fill(HIST("QA/Event/Sphericity"), collision.sphericity());

    // PID QA
    for (auto kCurrentTrack : dTracks) {
      if (kCurrentTrack.isKaon()) {
        FillKaons<true>(kCurrentTrack);
      }
      if (kCurrentTrack.isProton()) {
        FillProtons<true>(kCurrentTrack);
      }
    }

    Partition<TracksType> kaons = aod::tracktag::isKaon == true;
    kaons.bindTable(dTracks);
    Partition<TracksType> protons = aod::tracktag::isProton == true;
    protons.bindTable(dTracks);

    TLorentzVector lDecayDaughter1, lDecayDaughter2, lResonance;
    for (auto& [trk1, trk2] : combinations(CombinationsFullIndexPolicy(kaons, protons))) {

      // Trk1: Kaon, Trk2: Proton
      lDecayDaughter1.SetXYZM(trk1.px(), trk1.py(), trk1.pz(), massKa);
      lDecayDaughter2.SetXYZM(trk2.px(), trk2.py(), trk2.pz(), massPr);
      lResonance = lDecayDaughter1 + lDecayDaughter2;

      if (lResonance.Rapidity() > 0.5 || lResonance.Rapidity() < -0.5)
        continue;

      // Un-like sign pair only
      if (trk1.sign() * trk2.sign() < 0) {
        histos.fill(HIST("Analysis/lambda1520invmass"), lResonance.M());
        histos.fill(HIST("Analysis/h3lambda1520invmass"), collision.multTPCtemp(), lResonance.Pt(), lResonance.M());
      } else {
        if (trk1.globalIndex() == trk2.globalIndex())
          continue;
        histos.fill(HIST("Analysis/lambda1520invmassLS"), lResonance.M());
        histos.fill(HIST("Analysis/h3lambda1520invmassLS"), collision.multTPCtemp(), lResonance.Pt(), lResonance.M());

        // Like sign pair ++
        if (trk1.sign() > 0 && trk2.sign() > 0) {
          histos.fill(HIST("Analysis/lambda1520invmassLSkp"), lResonance.M());
          histos.fill(HIST("Analysis/h3lambda1520invmassLSkp"), collision.multTPCtemp(), lResonance.Pt(), lResonance.M());
        }

        // Like sign pair --
        if (trk1.sign() < 0 && trk2.sign() < 0) {
          histos.fill(HIST("Analysis/lambda1520invmassLSkbarpbar"), lResonance.M());
          histos.fill(HIST("Analysis/h3lambda1520invmassLSkbarpbar"), collision.multTPCtemp(), lResonance.Pt(), lResonance.M());
        }
      }
      if constexpr (IsMC) {
        if (abs(trk1.pdgCode()) != kKPlus || abs(trk2.pdgCode()) != kProton) // check if the tracks are kaons and Protons
          continue;
        auto mother1 = trk1.motherId();
        auto mother2 = trk2.motherId();
        if (mother1 == mother2) {         // Same mother
          if (trk1.motherPDG() == 3124) { // lambda1520
            histos.fill(HIST("MC/recolambda1520invmass"), lResonance.M());
            histos.fill(HIST("MC/recolambda1520pt"), lResonance.Pt());
            histos.fill(HIST("MC/h3recolambda1520invmass"), collision.multTPCtemp(), lResonance.Pt(), lResonance.M());
          }
        }
      }
    }
  }

  void processData(aod::ResoCollisions& collisions, soa::Filtered<TaggedTracks> const& resotracks)
  {
    LOGF(debug, "[DATA] Processing %d collisions", collisions.size());
    for (auto& collision : collisions) {
      Partition<soa::Filtered<TaggedTracks>> selectedTracks = (o2::aod::track::pt > static_cast<float_t>(cMinPtcut)) && (nabs(o2::aod::track::eta) < static_cast<float_t>(cfgCutEta)) && (nabs(o2::aod::track::dcaZ) > static_cast<float_t>(cMinDCAzToPVcut)) && (nabs(o2::aod::track::dcaZ) < static_cast<float_t>(cMaxDCAzToPVcut)) && (nabs(o2::aod::track::dcaXY) <= static_cast<float_t>(cMaxDCArToPVcut)) /* ((0.0105f + 0.0350f / npow(o2::aod::track::pt, 1.1f)))) */ && (aod::resodaughter::tpcNClsCrossedRows > static_cast<uint8_t>(cMinTPCncr)); // Basic DCA cuts
      selectedTracks.bindTable(resotracks);
      auto colTracks = selectedTracks->sliceByCached(aod::resodaughter::resoCollisionId, collision.globalIndex(), cache);
      fillHistograms<false>(collision, colTracks);
    }
  }
  PROCESS_SWITCH(lambda1520analysis, processData, "Process Event for data", true);

  void processMC(aod::ResoCollisions& collisions, soa::Filtered<TaggedMCTracks> const& resomctracks, aod::McParticles const& mcParticles)
  {
    LOGF(debug, "[MC] MC events: %d", collisions.size());
    for (auto& collision : collisions) {
      Partition<soa::Filtered<TaggedMCTracks>> selectedTracks = (o2::aod::track::pt > static_cast<float_t>(cMinPtcut)) && (nabs(o2::aod::track::eta) < static_cast<float_t>(cfgCutEta)) && (nabs(o2::aod::track::dcaZ) > static_cast<float_t>(cMinDCAzToPVcut)) && (nabs(o2::aod::track::dcaZ) < static_cast<float_t>(cMaxDCAzToPVcut)) && (nabs(o2::aod::track::dcaXY) <= static_cast<float_t>(cMaxDCArToPVcut)) && (aod::resodaughter::tpcNClsCrossedRows > static_cast<uint8_t>(cMinTPCncr)); // Basic DCA cuts
      selectedTracks.bindTable(resomctracks);
      auto colMCTracks = selectedTracks->sliceByCached(aod::resodaughter::resoCollisionId, collision.globalIndex(), cache);
      fillHistograms<true>(collision, colMCTracks);
    }

    // Not related to the real collisions
    for (auto& part : mcParticles) {             // loop over all MC particles
      if (abs(part.pdgCode()) == 3124) {         // Lambda(1520)
        if (part.y() > 0.5 || part.y() < -0.5) { // rapidity cut
          // LOGF(info, "[Rapidity cut] Lambda(1520): %d, y: %f", part.pdgCode(), part.y());
          continue;
        }
        bool pass1 = false;
        bool pass2 = false;
        for (auto& dau : part.daughters_as<aod::McParticles>()) {
          if (dau.pt() < cMinPtcut || fabs(dau.eta()) > cfgCutEta)
            continue;
          //  LOGF(info,"daughter pt: %f, eta: %f", dau.pt(), dau.eta());
          if (abs(dau.pdgCode()) == kKPlus) { // Decay to Kaons
            pass2 = true;
          }
          if (abs(dau.pdgCode()) == kProton) { // Decay to Protons
            pass1 = true;
          }
        }
        if (!pass1 || !pass2)
          continue;
        histos.fill(HIST("MC/truelambda1520pt"), part.pt());
      }
    }
  }
  PROCESS_SWITCH(lambda1520analysis, processMC, "Process Event for MC", false);

  // Processing Event Mixing
  using BinningTypeVetZTPCtemp = ColumnBinningPolicy<aod::collision::PosZ, aod::resocollision::MultTPCtemp>;
  void processME(o2::aod::ResoCollisions& collisions, soa::Filtered<TaggedTracks> const& resomixedtracks)
  {
    LOGF(debug, "Event Mixing Started");
    auto tracksTuple = std::make_tuple(resomixedtracks);
    BinningTypeVetZTPCtemp colBinning{{CfgVtxBins, CfgMultBins}, true};
    SameKindPair<aod::ResoCollisions, soa::Filtered<TaggedTracks>, BinningTypeVetZTPCtemp> pairs{colBinning, cfgNoMixedEvents, -1, collisions, tracksTuple, &cache}; // -1 is the number of the bin to skip

    TLorentzVector lDecayDaughter1, lDecayDaughter2, lResonance;
    for (auto& [collision1, tracks1, collision2, tracks2] : pairs) {
      // Kaons
      Partition<soa::Filtered<TaggedTracks>> selectedKaons = (o2::aod::track::pt > static_cast<float_t>(cMinPtcut)) && (nabs(o2::aod::track::eta) < static_cast<float_t>(cfgCutEta)) && (nabs(o2::aod::track::dcaZ) > static_cast<float_t>(cMinDCAzToPVcut)) && (nabs(o2::aod::track::dcaZ) < static_cast<float_t>(cMaxDCAzToPVcut)) && (nabs(o2::aod::track::dcaXY) < static_cast<float_t>(cMaxDCArToPVcut)) && (aod::resodaughter::tpcNClsCrossedRows > static_cast<uint8_t>(cMinTPCncr)) && aod::tracktag::isKaon == true; // Basic DCA cuts
      selectedKaons.bindTable(tracks1);
      // Protons
      Partition<soa::Filtered<TaggedTracks>> selectedProtons = (o2::aod::track::pt > static_cast<float_t>(cMinPtcut)) && (nabs(o2::aod::track::eta) < static_cast<float_t>(cfgCutEta)) && (nabs(o2::aod::track::dcaZ) > static_cast<float_t>(cMinDCAzToPVcut)) && (nabs(o2::aod::track::dcaZ) < static_cast<float_t>(cMaxDCAzToPVcut)) && (nabs(o2::aod::track::dcaXY) < static_cast<float_t>(cMaxDCArToPVcut)) && (aod::resodaughter::tpcNClsCrossedRows > static_cast<uint8_t>(cMinTPCncr)) && aod::tracktag::isProton == true; // Basic DCA cuts
      selectedProtons.bindTable(tracks2);

      for (auto& [trk1, trk2] : combinations(CombinationsFullIndexPolicy(selectedKaons, selectedProtons))) {
        // Un-like sign pair only
        if (trk1.sign() * trk2.sign() > 0)
          continue;

        lDecayDaughter1.SetXYZM(trk1.px(), trk1.py(), trk1.pz(), massKa);
        lDecayDaughter2.SetXYZM(trk2.px(), trk2.py(), trk2.pz(), massPr);
        lResonance = lDecayDaughter1 + lDecayDaughter2;

        if (lResonance.Rapidity() > 0.5 || lResonance.Rapidity() < -0.5)
          continue;

        histos.fill(HIST("Analysis/lambda1520invmassME"), lResonance.M());
        histos.fill(HIST("Analysis/h3lambda1520invmassME"), collision1.multTPCtemp(), lResonance.Pt(), lResonance.M());
      }
    }
  };
  PROCESS_SWITCH(lambda1520analysis, processME, "Process EventMixing", true);

  template <bool fillHistos, typename kCurrentTrackType>
  bool FillKaons(kCurrentTrackType const& kCurrentTrack)
  {

    if constexpr (fillHistos) {
      //  --- PID QA Kaons
      histos.fill(HIST("QA/Kaon/TOF_Nsigma"), kCurrentTrack.pt(), kCurrentTrack.tofNSigmaKa());
      histos.fill(HIST("QA/Kaon/TPC_Nsigma"), kCurrentTrack.pt(), kCurrentTrack.tpcNSigmaKa());
      histos.fill(HIST("QA/Kaon/TOF_TPC_Map"), kCurrentTrack.tofNSigmaKa(), kCurrentTrack.tpcNSigmaKa());
      histos.fill(HIST("QA/Kaon/pT"), kCurrentTrack.pt());
      histos.fill(HIST("QA/Kaon/dcaZ"), kCurrentTrack.pt(), kCurrentTrack.dcaZ());
      histos.fill(HIST("QA/Kaon/dcaXY"), kCurrentTrack.pt(), kCurrentTrack.dcaXY());
      histos.fill(HIST("QA/Kaon/TPC_CR"), kCurrentTrack.pt(), kCurrentTrack.tpcNClsCrossedRows());
    }
    return true;
  }

  template <bool fillHistos, typename kCurrentTrackType>
  bool FillProtons(kCurrentTrackType const& kCurrentTrack)
  {

    if constexpr (fillHistos) {
      //  --- PID QA Protons
      histos.fill(HIST("QA/Proton/TOF_Nsigma"), kCurrentTrack.pt(), kCurrentTrack.tofNSigmaPr());
      histos.fill(HIST("QA/Proton/TPC_Nsigma"), kCurrentTrack.pt(), kCurrentTrack.tpcNSigmaPr());
      histos.fill(HIST("QA/Proton/TOF_TPC_Map"), kCurrentTrack.tofNSigmaPr(), kCurrentTrack.tpcNSigmaPr());
      histos.fill(HIST("QA/Proton/pT"), kCurrentTrack.pt());
      histos.fill(HIST("QA/Proton/dcaZ"), kCurrentTrack.pt(), kCurrentTrack.dcaZ());
      histos.fill(HIST("QA/Proton/dcaXY"), kCurrentTrack.pt(), kCurrentTrack.dcaXY());
      histos.fill(HIST("QA/Proton/TPC_CR"), kCurrentTrack.pt(), kCurrentTrack.tpcNClsCrossedRows());
    }
    return true;
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<TrackPreselector>(cfgc),
                      adaptAnalysisTask<lambda1520analysis>(cfgc, TaskName{"lf-lambda1520analysis"})};
}
