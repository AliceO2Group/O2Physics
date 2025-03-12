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

/// \file taskQaPidRejection.cxx
///
/// \author Peter Hristov <Peter.Hristov@cern.ch>, CERN
/// \author Gian Michele Innocenti <gian.michele.innocenti@cern.ch>, CERN
/// \author Henrique J C Zanoli <henrique.zanoli@cern.ch>, Utrecht University
/// \author Nicolo' Jacazio <nicolo.jacazio@cern.ch>, CERN

#include <utility>
#include <vector>

#include <TEfficiency.h>
#include <TList.h>
#include <TPDGCode.h>

#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "ReconstructionDataFormats/DCA.h"

#include "ALICE3/DataModel/MID.h"
#include "ALICE3/DataModel/RICH.h"
#include "Common/Core/TrackSelectorPID.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "PWGHF/DataModel/CandidateReconstructionTables.h"

using namespace o2;
using namespace o2::framework;

namespace o2::aod
{
namespace hf_track_index_alice3_pid
{
DECLARE_SOA_INDEX_COLUMN(Track, track); //!
DECLARE_SOA_INDEX_COLUMN(RICH, rich);   //!
DECLARE_SOA_INDEX_COLUMN(MID, mid);     //!
} // namespace hf_track_index_alice3_pid

DECLARE_SOA_INDEX_TABLE_USER(HfTrackIndexALICE3PID, Tracks, "HFTRKIDXA3PID", //!
                             hf_track_index_alice3_pid::TrackId,
                             hf_track_index_alice3_pid::RICHId,
                             hf_track_index_alice3_pid::MIDId);
} // namespace o2::aod

struct HfTaskQaPidRejectionAlice3PidIndexBuilder {
  Builds<aod::HfTrackIndexALICE3PID> index;
  void init(InitContext&) {}
};

void customize(std::vector<ConfigParamSpec>& workflowOptions)
{
  std::vector<ConfigParamSpec> options{
    {"rej-el", VariantType::Int, 1, {"Efficiency for the Electron PDG code"}},
    {"rej-mu", VariantType::Int, 1, {"Efficiency for the Muon PDG code"}},
    {"rej-pi", VariantType::Int, 1, {"Efficiency for the Pion PDG code"}},
    {"rej-ka", VariantType::Int, 1, {"Efficiency for the Kaon PDG code"}},
    {"rej-pr", VariantType::Int, 0, {"Efficiency for the Proton PDG code"}}};
  std::swap(workflowOptions, options);
}

#include "Framework/runDataProcessing.h"

/// Task to QA the efficiency of a particular particle defined by particlePDG
template <track::pid_constants::ID particle>
struct HfTaskQaPidRejection {
  // Particle selection
  Configurable<int> nBinsEta{"nBinsEta", 40, "Number of eta bins"};
  Configurable<float> etaMin{"etaMin", -2.f, "Lower limit in eta"};
  Configurable<float> etaMax{"etaMax", 2.f, "Upper limit in eta"};
  Configurable<int> nBinsPt{"nBinsPt", 400, "Number of pT bins"};
  Configurable<float> ptMin{"ptMin", 0.f, "Lower limit in pT"};
  Configurable<float> ptMax{"ptMax", 20.f, "Upper limit in pT"};
  // TPC
  Configurable<double> ptPidTpcMin{"ptPidTpcMin", 9999., "Lower bound of track pT for TPC PID"};
  Configurable<double> ptPidTpcMax{"ptPidTpcMax", 999999., "Upper bound of track pT for TPC PID"};
  Configurable<double> nSigmaTpcMax{"nSigmaTpcMax", 99999, "Nsigma cut on TPC only"};
  // TOF
  Configurable<double> ptPidTofMin{"ptPidTofMin", 0.15, "Lower bound of track pT for TOF PID"};
  Configurable<double> ptPidTofMax{"ptPidTofMax", 5., "Upper bound of track pT for TOF PID"};
  Configurable<double> nSigmaTofMax{"nSigmaTofMax", 3., "Nsigma cut on TOF only"};
  Configurable<double> nSigmaTofCombinedMax{"nSigmaTofCombinedMax", 0., "Nsigma cut on TOF combined with TPC"};
  // RICH
  Configurable<double> ptPidRichMin{"ptPidRichMin", 0.15, "Lower bound of track pT for RICH PID"};
  Configurable<double> ptPidRichMax{"ptPidRichMax", 10., "Upper bound of track pT for RICH PID"};
  Configurable<double> nSigmaRichMax{"nSigmaRichMax", 3., "Nsigma cut on RICH only"};
  Configurable<double> nSigmaRichCombinedTofMax{"nSigmaRichCombinedTofMax", 0., "Nsigma cut on RICH combined with TOF"};

  TrackSelectorEl selectorElectron;
  TrackSelectorMu selectorMuon;
  TrackSelectorPi selectorPion;
  TrackSelectorKa selectorKaon;
  TrackSelectorPr selectorProton;

  static constexpr PDG_t PDGs[5] = {kElectron, kMuonMinus, kPiPlus, kKPlus, kProton};
  static_assert(particle < 5 && "Maximum of particles reached");
  static constexpr int particlePDG = PDGs[particle];

  using TracksWPid = soa::Join<aod::Tracks, aod::pidTOFFullEl, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::HfTrackIndexALICE3PID>;

  HistogramRegistry histos{"HistogramsRejection"};

  void init(InitContext&)
  {
    selectorElectron.setRangePtTpc(ptPidTpcMin, ptPidTpcMax);
    selectorElectron.setRangeNSigmaTpc(-nSigmaTpcMax, nSigmaTpcMax);
    selectorElectron.setRangePtTof(ptPidTofMin, ptPidTofMax);
    selectorElectron.setRangeNSigmaTof(-nSigmaTofMax, nSigmaTofMax);
    selectorElectron.setRangeNSigmaTofCondTpc(-nSigmaTofCombinedMax, nSigmaTofCombinedMax);
    selectorElectron.setRangePtRich(ptPidRichMin, ptPidRichMax);
    selectorElectron.setRangeNSigmaRich(-nSigmaRichMax, nSigmaRichMax);
    selectorElectron.setRangeNSigmaRichCondTof(-nSigmaRichCombinedTofMax, nSigmaRichCombinedTofMax);
    selectorMuon = selectorElectron;
    selectorPion = selectorElectron;
    selectorKaon = selectorElectron;
    selectorProton = selectorElectron;

    AxisSpec ptAxis{nBinsPt, ptMin, ptMax};
    AxisSpec etaAxis{nBinsEta, etaMin, etaMax};

    TString commonTitle = "";
    if (particlePDG != 0) {
      commonTitle += Form("PDG %i", particlePDG);
    }

    const TString pt = "#it{p}_{T} [GeV/#it{c}]";
    const TString p = "#it{p} [GeV/#it{c}]";
    const TString eta = "#it{#eta}";
    const TString phi = "#it{#varphi} [rad]";

    histos.add("tracking/pteta", commonTitle + " Primary;" + pt, kTH2D, {ptAxis, etaAxis});
    histos.add("trackingTOFselElectron/pteta", commonTitle + " Primary;" + pt, kTH2D, {ptAxis, etaAxis});
    histos.add("trackingTOFselPion/pteta", commonTitle + " Primary;" + pt, kTH2D, {ptAxis, etaAxis});
    histos.add("trackingTOFselKaon/pteta", commonTitle + " Primary;" + pt, kTH2D, {ptAxis, etaAxis});
    histos.add("trackingTOFselProton/pteta", commonTitle + " Primary;" + pt, kTH2D, {ptAxis, etaAxis});
    histos.add("trackingRICHselElectron/pteta", commonTitle + " Primary;" + pt, kTH2D, {ptAxis, etaAxis});
    histos.add("trackingRICHselPion/pteta", commonTitle + " Primary;" + pt, kTH2D, {ptAxis, etaAxis});
    histos.add("trackingRICHselKaon/pteta", commonTitle + " Primary;" + pt, kTH2D, {ptAxis, etaAxis});
    histos.add("trackingRICHselProton/pteta", commonTitle + " Primary;" + pt, kTH2D, {ptAxis, etaAxis});
    histos.add("trackingMIDselMuon/pteta", commonTitle + " Primary;" + pt, kTH2D, {ptAxis, etaAxis});

    histos.add("tracking/peta", commonTitle + " Primary;" + p, kTH2D, {ptAxis, etaAxis});
    histos.add("trackingTOFselElectron/peta", commonTitle + " Primary;" + p, kTH2D, {ptAxis, etaAxis});
    histos.add("trackingTOFselPion/peta", commonTitle + " Primary;" + p, kTH2D, {ptAxis, etaAxis});
    histos.add("trackingTOFselKaon/peta", commonTitle + " Primary;" + p, kTH2D, {ptAxis, etaAxis});
    histos.add("trackingTOFselProton/peta", commonTitle + " Primary;" + p, kTH2D, {ptAxis, etaAxis});
    histos.add("trackingRICHselElectron/peta", commonTitle + " Primary;" + p, kTH2D, {ptAxis, etaAxis});
    histos.add("trackingRICHselPion/peta", commonTitle + " Primary;" + p, kTH2D, {ptAxis, etaAxis});
    histos.add("trackingRICHselKaon/peta", commonTitle + " Primary;" + p, kTH2D, {ptAxis, etaAxis});
    histos.add("trackingRICHselProton/peta", commonTitle + " Primary;" + p, kTH2D, {ptAxis, etaAxis});
    histos.add("trackingMIDselMuon/peta", commonTitle + " Primary;" + p, kTH2D, {ptAxis, etaAxis});
  }

  void process(const soa::Join<aod::Collisions, aod::McCollisionLabels>& collisions,
               const soa::Join<TracksWPid, aod::McTrackLabels>& tracks,
               const aod::McCollisions&,
               const aod::McParticles&,
               const aod::RICHs&,
               const aod::MIDs&)
  {
    std::vector<int64_t> recoEvt(collisions.size());
    std::vector<int64_t> recoTracks(tracks.size());
    LOGF(info, "%d", particlePDG);
    for (const auto& track : tracks) {
      const auto mcParticle = track.mcParticle();
      if (particlePDG != 0 && mcParticle.pdgCode() != particlePDG) { // Checking PDG code
        continue;
      }
      bool isTOFhpElectron = !(selectorElectron.statusTof(track) == TrackSelectorPID::Rejected);
      bool isRICHhpElectron = !(selectorElectron.statusRich(track) == TrackSelectorPID::Rejected);
      bool isTOFhpPion = !(selectorPion.statusTof(track) == TrackSelectorPID::Rejected);
      bool isRICHhpPion = !(selectorPion.statusRich(track) == TrackSelectorPID::Rejected);
      bool isTOFhpKaon = !(selectorKaon.statusTof(track) == TrackSelectorPID::Rejected);
      bool isRICHhpKaon = !(selectorKaon.statusRich(track) == TrackSelectorPID::Rejected);
      bool isTOFhpProton = !(selectorProton.statusTof(track) == TrackSelectorPID::Rejected);
      bool isRICHhpProton = !(selectorProton.statusRich(track) == TrackSelectorPID::Rejected);
      bool isMIDhpMuon = (selectorMuon.statusMid(track) == TrackSelectorPID::Accepted);

      if (mcParticle.isPhysicalPrimary()) {
        histos.fill(HIST("tracking/pteta"), track.pt(), track.eta());
        if (isTOFhpElectron) {
          histos.fill(HIST("trackingTOFselElectron/pteta"), track.pt(), track.eta());
        }
        if (isTOFhpPion) {
          histos.fill(HIST("trackingTOFselPion/pteta"), track.pt(), track.eta());
        }
        if (isTOFhpKaon) {
          histos.fill(HIST("trackingTOFselKaon/pteta"), track.pt(), track.eta());
        }
        if (isTOFhpProton) {
          histos.fill(HIST("trackingTOFselProton/pteta"), track.pt(), track.eta());
        }
        if (isRICHhpElectron) {
          histos.fill(HIST("trackingRICHselElectron/pteta"), track.pt(), track.eta());
        }
        if (isRICHhpPion) {
          histos.fill(HIST("trackingRICHselPion/pteta"), track.pt(), track.eta());
        }
        if (isRICHhpKaon) {
          histos.fill(HIST("trackingRICHselKaon/pteta"), track.pt(), track.eta());
        }
        if (isRICHhpProton) {
          histos.fill(HIST("trackingRICHselProton/pteta"), track.pt(), track.eta());
        }
        if (isMIDhpMuon) {
          histos.fill(HIST("trackingMIDselMuon/pteta"), track.pt(), track.eta());
        }

        histos.fill(HIST("tracking/peta"), track.p(), track.eta());
        if (isTOFhpElectron) {
          histos.fill(HIST("trackingTOFselElectron/peta"), track.p(), track.eta());
        }
        if (isTOFhpPion) {
          histos.fill(HIST("trackingTOFselPion/peta"), track.p(), track.eta());
        }
        if (isTOFhpKaon) {
          histos.fill(HIST("trackingTOFselKaon/peta"), track.p(), track.eta());
        }
        if (isTOFhpProton) {
          histos.fill(HIST("trackingTOFselProton/peta"), track.p(), track.eta());
        }
        if (isRICHhpElectron) {
          histos.fill(HIST("trackingRICHselElectron/peta"), track.p(), track.eta());
        }
        if (isRICHhpPion) {
          histos.fill(HIST("trackingRICHselPion/peta"), track.p(), track.eta());
        }
        if (isRICHhpKaon) {
          histos.fill(HIST("trackingRICHselKaon/peta"), track.p(), track.eta());
        }
        if (isRICHhpProton) {
          histos.fill(HIST("trackingRICHselProton/peta"), track.p(), track.eta());
        }
        if (isMIDhpMuon) {
          histos.fill(HIST("trackingMIDselMuon/peta"), track.p(), track.eta());
        }
      }
    }
  }
};

struct HfTaskQaPidRejectionGeneral {
  static constexpr PDG_t PDGs[5] = {kElectron, kMuonMinus, kPiPlus, kKPlus, kProton};
  // Cuts
  Configurable<float> etaMaxSel{"etaMaxSel", 1.44, "Max #eta single track"};
  Configurable<float> ptMinSel{"ptMinSel", 0.6, "p_{T} min single track"};
  // Particle selection
  Configurable<int> nBinsEta{"nBinsEta", 40, "Number of eta bins"};
  Configurable<float> etaMin{"etaMin", -2.f, "Lower limit in eta"};
  Configurable<float> etaMax{"etaMax", 2.f, "Upper limit in eta"};
  Configurable<int> nBinsPt{"nBinsPt", 400, "Number of pT bins"};
  Configurable<float> ptMin{"ptMin", 0.f, "Lower limit in pT"};
  Configurable<float> ptMax{"ptMax", 20.f, "Upper limit in pT"};
  // TPC
  Configurable<double> ptPidTpcMin{"ptPidTpcMin", 9999., "Lower bound of track pT for TPC PID"};
  Configurable<double> ptPidTpcMax{"ptPidTpcMax", 999999., "Upper bound of track pT for TPC PID"};
  Configurable<double> nSigmaTpcMax{"nSigmaTpcMax", 99999, "Nsigma cut on TPC only"};
  // TOF
  Configurable<double> ptPidTofMin{"ptPidTofMin", 0.15, "Lower bound of track pT for TOF PID"};
  Configurable<double> ptPidTofMax{"ptPidTofMax", 5., "Upper bound of track pT for TOF PID"};
  Configurable<double> nSigmaTofMax{"nSigmaTofMax", 3., "Nsigma cut on TOF only"};
  Configurable<double> nSigmaTofCombinedMax{"nSigmaTofCombinedMax", 0., "Nsigma cut on TOF combined with TPC"};
  // RICH
  Configurable<double> ptPidRichMin{"ptPidRichMin", 0.15, "Lower bound of track pT for RICH PID"};
  Configurable<double> ptPidRichMax{"ptPidRichMax", 10., "Upper bound of track pT for RICH PID"};
  Configurable<double> nSigmaRichMax{"nSigmaRichMax", 3., "Nsigma cut on RICH only"};
  Configurable<double> nSigmaRichCombinedTofMax{"nSigmaRichCombinedTofMax", 0., "Nsigma cut on RICH combined with TOF"};

  TrackSelectorEl selectorElectron;
  TrackSelectorMu selectorMuon;
  TrackSelectorPi selectorPion;
  TrackSelectorKa selectorKaon;
  TrackSelectorPr selectorProton;

  using TracksWPid = soa::Join<aod::Tracks, aod::pidTOFFullEl, aod::pidTOFFullPi, aod::HfTrackIndexALICE3PID>;

  HistogramRegistry histos{"HistogramsRejection"};

  void init(InitContext&)
  {
    selectorElectron.setRangePtTpc(ptPidTpcMin, ptPidTpcMax);
    selectorElectron.setRangeNSigmaTpc(-nSigmaTpcMax, nSigmaTpcMax);
    selectorElectron.setRangePtTof(ptPidTofMin, ptPidTofMax);
    selectorElectron.setRangeNSigmaTof(-nSigmaTofMax, nSigmaTofMax);
    selectorElectron.setRangeNSigmaTofCondTpc(-nSigmaTofCombinedMax, nSigmaTofCombinedMax);
    selectorElectron.setRangePtRich(ptPidRichMin, ptPidRichMax);
    selectorElectron.setRangeNSigmaRich(-nSigmaRichMax, nSigmaRichMax);
    selectorElectron.setRangeNSigmaRichCondTof(-nSigmaRichCombinedTofMax, nSigmaRichCombinedTofMax);
    selectorMuon = selectorElectron;
    selectorPion = selectorElectron;
    selectorKaon = selectorElectron;
    selectorProton = selectorElectron;

    AxisSpec ptAxis{nBinsPt, ptMin, ptMax};
    AxisSpec etaAxis{nBinsEta, etaMin, etaMax};

    TString commonTitle = "";

    const TString pt = "#it{p}_{T} [GeV/#it{c}]";
    const TString p = "#it{p} [GeV/#it{c}]";
    const TString eta = "#it{#eta}";
    const TString phi = "#it{#varphi} [rad]";

    histos.add("hAllNoSel/pteta", commonTitle + " Primary;" + pt, kTH2D, {ptAxis, etaAxis});
    histos.add("hPionNoSel/pteta", commonTitle + " Primary;" + pt, kTH2D, {ptAxis, etaAxis});
    histos.add("hElectronNoSel/pteta", commonTitle + " Primary;" + pt, kTH2D, {ptAxis, etaAxis});
    histos.add("hMuonNoSel/pteta", commonTitle + " Primary;" + pt, kTH2D, {ptAxis, etaAxis});
    histos.add("hKaonNoSel/pteta", commonTitle + " Primary;" + pt, kTH2D, {ptAxis, etaAxis});

    histos.add("hAllRICHSelHpElTight/pteta", commonTitle + " Primary;" + pt, kTH2D, {ptAxis, etaAxis});
    histos.add("hPionRICHSelHpElTight/pteta", commonTitle + " Primary;" + pt, kTH2D, {ptAxis, etaAxis});
    histos.add("hElectronRICHSelHpElTight/pteta", commonTitle + " Primary;" + pt, kTH2D, {ptAxis, etaAxis});
    histos.add("hMuonRICHSelHpElTight/pteta", commonTitle + " Primary;" + pt, kTH2D, {ptAxis, etaAxis});
    histos.add("hKaonRICHSelHpElTight/pteta", commonTitle + " Primary;" + pt, kTH2D, {ptAxis, etaAxis});

    histos.add("hAllRICHSelHpElTightAlt/pteta", commonTitle + " Primary;" + pt, kTH2D, {ptAxis, etaAxis});
    histos.add("hPionRICHSelHpElTightAlt/pteta", commonTitle + " Primary;" + pt, kTH2D, {ptAxis, etaAxis});
    histos.add("hElectronRICHSelHpElTightAlt/pteta", commonTitle + " Primary;" + pt, kTH2D, {ptAxis, etaAxis});
    histos.add("hMuonRICHSelHpElTightAlt/pteta", commonTitle + " Primary;" + pt, kTH2D, {ptAxis, etaAxis});
    histos.add("hKaonRICHSelHpElTightAlt/pteta", commonTitle + " Primary;" + pt, kTH2D, {ptAxis, etaAxis});

    histos.add("hAllRICHSelHpElTightAltDiff/pteta", commonTitle + " Primary;" + pt, kTH2D, {ptAxis, etaAxis});
    histos.add("hPionRICHSelHpElTightAltDiff/pteta", commonTitle + " Primary;" + pt, kTH2D, {ptAxis, etaAxis});
    histos.add("hElectronRICHSelHpElTightAltDiff/pteta", commonTitle + " Primary;" + pt, kTH2D, {ptAxis, etaAxis});
    histos.add("hMuonRICHSelHpElTightAltDiff/pteta", commonTitle + " Primary;" + pt, kTH2D, {ptAxis, etaAxis});
    histos.add("hKaonRICHSelHpElTightAltDiff/pteta", commonTitle + " Primary;" + pt, kTH2D, {ptAxis, etaAxis});

    histos.add("hAllRICHSelHpElLoose/pteta", commonTitle + " Primary;" + pt, kTH2D, {ptAxis, etaAxis});
    histos.add("hPionRICHSelHpElLoose/pteta", commonTitle + " Primary;" + pt, kTH2D, {ptAxis, etaAxis});
    histos.add("hElectronRICHSelHpElLoose/pteta", commonTitle + " Primary;" + pt, kTH2D, {ptAxis, etaAxis});
    histos.add("hMuonRICHSelHpElLoose/pteta", commonTitle + " Primary;" + pt, kTH2D, {ptAxis, etaAxis});
    histos.add("hKaonRICHSelHpElLoose/pteta", commonTitle + " Primary;" + pt, kTH2D, {ptAxis, etaAxis});

    histos.add("hAllMID/pteta", commonTitle + " Primary;" + pt, kTH2D, {ptAxis, etaAxis});
    histos.add("hElectronMID/pteta", commonTitle + " Primary;" + pt, kTH2D, {ptAxis, etaAxis});
    histos.add("hMuonMID/pteta", commonTitle + " Primary;" + pt, kTH2D, {ptAxis, etaAxis});
    histos.add("hPionMID/pteta", commonTitle + " Primary;" + pt, kTH2D, {ptAxis, etaAxis});
    histos.add("hKaonMID/pteta", commonTitle + " Primary;" + pt, kTH2D, {ptAxis, etaAxis});

    histos.add("hAllMID/peta", commonTitle + " Primary;" + p, kTH2D, {ptAxis, etaAxis});
    histos.add("hElectronMID/peta", commonTitle + " Primary;" + p, kTH2D, {ptAxis, etaAxis});
    histos.add("hMuonMID/peta", commonTitle + " Primary;" + p, kTH2D, {ptAxis, etaAxis});
    histos.add("hPionMID/peta", commonTitle + " Primary;" + p, kTH2D, {ptAxis, etaAxis});
    histos.add("hKaonMID/peta", commonTitle + " Primary;" + p, kTH2D, {ptAxis, etaAxis});
  }

  void process(const soa::Join<aod::Collisions, aod::McCollisionLabels>&,
               const soa::Join<TracksWPid, aod::McTrackLabels>& tracks,
               const aod::McCollisions&,
               const aod::McParticles&,
               const aod::RICHs&,
               const aod::MIDs&)
  {
    for (const auto& track : tracks) {

      if (std::abs(track.eta()) > etaMaxSel || track.pt() < ptMinSel) {
        continue;
      }
      const auto mcParticle = track.mcParticle();
      histos.fill(HIST("hAllNoSel/pteta"), track.pt(), track.eta());
      if (mcParticle.pdgCode() == kElectron) {
        histos.fill(HIST("hElectronNoSel/pteta"), track.pt(), track.eta());
      }
      if (mcParticle.pdgCode() == kMuonPlus) {
        histos.fill(HIST("hMuonNoSel/pteta"), track.pt(), track.eta());
      }
      if (mcParticle.pdgCode() == kPiPlus) {
        histos.fill(HIST("hPionNoSel/pteta"), track.pt(), track.eta());
      }
      if (mcParticle.pdgCode() == kKPlus) {
        histos.fill(HIST("hKaonNoSel/pteta"), track.pt(), track.eta());
      }

      bool isRICHhpElectron = !(selectorElectron.statusRich(track) == TrackSelectorPID::Rejected);
      // bool isRICHhpPion = !(selectorPion.statusRich(track) == TrackSelectorPID::Rejected);
      // bool isRICHhpKaon = !(selectorKaon.statusRich(track) == TrackSelectorPID::Rejected);
      // bool isRICHhpProton = !(selectorProton.statusRich(track) == TrackSelectorPID::Rejected);
      bool isMIDhpMuon = (selectorMuon.statusMid(track) == TrackSelectorPID::Accepted);

      bool isRICHElLoose = isRICHhpElectron;

      bool isRICHElTight = true;
      bool isrichel = true;
      bool istofel = true;
      bool istofpi = false;
      bool isrichpi = false;

      if (track.p() < 0.6) {
        istofel = std::abs(track.tofNSigmaEl()) < ptPidTofMax;
        istofpi = std::abs(track.tofNSigmaPi()) < ptPidTofMax;
        isrichel = std::abs(track.rich().richNsigmaEl()) < ptPidRichMax;
      } else if (track.p() > 0.6 && track.p() < 1.0) {
        isrichel = std::abs(track.rich().richNsigmaEl()) < ptPidRichMax;
      } else if (track.p() > 1.0 && track.p() < 2.0) {
        isrichel = std::abs(track.rich().richNsigmaEl()) < ptPidRichMax;
        isrichpi = std::abs(track.rich().richNsigmaPi()) < ptPidRichMax;
      } else {
        isrichel = std::abs(track.rich().richNsigmaEl()) < ptPidRichMax;
      }
      isRICHElTight = isrichel && !isrichpi && istofel && !istofpi;

      auto isRICHElTightAlt = selectorElectron.isElectronAndNotPion(track);

      if (isRICHElLoose) {
        histos.fill(HIST("hAllRICHSelHpElLoose/pteta"), track.pt(), track.eta());
        if (mcParticle.pdgCode() == kElectron) {
          histos.fill(HIST("hElectronRICHSelHpElLoose/pteta"), track.pt(), track.eta());
        }
        if (mcParticle.pdgCode() == kMuonPlus) {
          histos.fill(HIST("hMuonRICHSelHpElLoose/pteta"), track.pt(), track.eta());
        }
        if (mcParticle.pdgCode() == kPiPlus) {
          histos.fill(HIST("hPionRICHSelHpElLoose/pteta"), track.pt(), track.eta());
        }
        if (mcParticle.pdgCode() == kKPlus) {
          histos.fill(HIST("hKaonRICHSelHpElLoose/pteta"), track.pt(), track.eta());
        }
      }
      if (isRICHElTight) {
        histos.fill(HIST("hAllRICHSelHpElTight/pteta"), track.pt(), track.eta());
        if (mcParticle.pdgCode() == kElectron) {
          histos.fill(HIST("hElectronRICHSelHpElTight/pteta"), track.pt(), track.eta());
        }
        if (mcParticle.pdgCode() == kMuonPlus) {
          histos.fill(HIST("hMuonRICHSelHpElTight/pteta"), track.pt(), track.eta());
        }
        if (mcParticle.pdgCode() == kPiPlus) {
          histos.fill(HIST("hPionRICHSelHpElTight/pteta"), track.pt(), track.eta());
        }
        if (mcParticle.pdgCode() == kKPlus) {
          histos.fill(HIST("hKaonRICHSelHpElTight/pteta"), track.pt(), track.eta());
        }
      }
      if (isRICHElTightAlt) {
        histos.fill(HIST("hAllRICHSelHpElTightAlt/pteta"), track.pt(), track.eta());
        if (mcParticle.pdgCode() == kElectron) {
          histos.fill(HIST("hElectronRICHSelHpElTightAlt/pteta"), track.pt(), track.eta());
        }
        if (mcParticle.pdgCode() == kMuonPlus) {
          histos.fill(HIST("hMuonRICHSelHpElTightAlt/pteta"), track.pt(), track.eta());
        }
        if (mcParticle.pdgCode() == kPiPlus) {
          histos.fill(HIST("hPionRICHSelHpElTightAlt/pteta"), track.pt(), track.eta());
        }
        if (mcParticle.pdgCode() == kKPlus) {
          histos.fill(HIST("hKaonRICHSelHpElTightAlt/pteta"), track.pt(), track.eta());
        }
      }
      if (isRICHElTightAlt != isRICHElTight) {
        histos.fill(HIST("hAllRICHSelHpElTightAltDiff/pteta"), track.pt(), track.eta());
        if (mcParticle.pdgCode() == kElectron) {
          histos.fill(HIST("hElectronRICHSelHpElTightAltDiff/pteta"), track.pt(), track.eta());
        }
        if (mcParticle.pdgCode() == kMuonPlus) {
          histos.fill(HIST("hMuonRICHSelHpElTightAltDiff/pteta"), track.pt(), track.eta());
        }
        if (mcParticle.pdgCode() == kPiPlus) {
          histos.fill(HIST("hPionRICHSelHpElTightAltDiff/pteta"), track.pt(), track.eta());
        }
        if (mcParticle.pdgCode() == kKPlus) {
          histos.fill(HIST("hKaonRICHSelHpElTightAltDiff/pteta"), track.pt(), track.eta());
        }
      }
      if (isMIDhpMuon) {
        histos.fill(HIST("hAllMID/pteta"), track.pt(), track.eta());
        if (mcParticle.pdgCode() == kElectron) {
          histos.fill(HIST("hElectronMID/pteta"), track.pt(), track.eta());
        }
        if (mcParticle.pdgCode() == kMuonPlus) {
          histos.fill(HIST("hMuonMID/pteta"), track.pt(), track.eta());
        }
        if (mcParticle.pdgCode() == kPiPlus) {
          histos.fill(HIST("hPionMID/pteta"), track.pt(), track.eta());
        }
        if (mcParticle.pdgCode() == kKPlus) {
          histos.fill(HIST("hKaonMID/pteta"), track.pt(), track.eta());
        }

        histos.fill(HIST("hAllMID/peta"), track.p(), track.eta());
        if (mcParticle.pdgCode() == kElectron) {
          histos.fill(HIST("hElectronMID/peta"), track.p(), track.eta());
        }
        if (mcParticle.pdgCode() == kMuonPlus) {
          histos.fill(HIST("hMuonMID/peta"), track.p(), track.eta());
        }
        if (mcParticle.pdgCode() == kPiPlus) {
          histos.fill(HIST("hPionMID/peta"), track.p(), track.eta());
        }
        if (mcParticle.pdgCode() == kKPlus) {
          histos.fill(HIST("hKaonMID/peta"), track.p(), track.eta());
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec w;
  w.push_back(adaptAnalysisTask<HfTaskQaPidRejectionAlice3PidIndexBuilder>(cfgc));
  if (cfgc.options().get<int>("rej-el")) {
    w.push_back(adaptAnalysisTask<HfTaskQaPidRejection<track::PID::Electron>>(cfgc, TaskName{"hf-task-qa-pid-rejection-electron"}));
  }
  if (cfgc.options().get<int>("rej-ka")) {
    w.push_back(adaptAnalysisTask<HfTaskQaPidRejection<track::PID::Kaon>>(cfgc, TaskName{"hf-task-qa-pid-rejection-kaon"}));
  }
  if (cfgc.options().get<int>("rej-pr")) {
    w.push_back(adaptAnalysisTask<HfTaskQaPidRejection<track::PID::Proton>>(cfgc, TaskName{"hf-task-qa-pid-rejection-proton"}));
  }
  if (cfgc.options().get<int>("rej-mu")) {
    w.push_back(adaptAnalysisTask<HfTaskQaPidRejection<track::PID::Muon>>(cfgc, TaskName{"hf-task-qa-pid-rejection-mu"}));
  }
  if (cfgc.options().get<int>("rej-pi")) {
    w.push_back(adaptAnalysisTask<HfTaskQaPidRejection<track::PID::Pion>>(cfgc, TaskName{"hf-task-qa-pid-rejection-pion"}));
  }
  w.push_back(adaptAnalysisTask<HfTaskQaPidRejectionGeneral>(cfgc));
  return w;
}
