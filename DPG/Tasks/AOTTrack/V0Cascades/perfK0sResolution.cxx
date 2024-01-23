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
#include "ReconstructionDataFormats/Track.h"
#include "ReconstructionDataFormats/PID.h"
#include "Common/Core/trackUtilities.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponse.h"

using namespace o2;
using namespace o2::track;
using namespace o2::framework;
using namespace o2::framework::expressions;

using PIDTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::pidTPCFullPi, aod::pidTPCFullPr, aod::pidTOFFullPi>;
using PIDTracksIU = soa::Join<aod::TracksIU, aod::TracksExtra, aod::pidTPCFullPi, aod::pidTPCFullPr, aod::pidTOFFullPi>;
using SelectedCollisions = soa::Join<aod::Collisions, aod::EvSels>;

struct perfK0sResolution {
  ConfigurableAxis mBins{"mBins", {200, 0.4f, 0.6f}, "Mass binning"};
  ConfigurableAxis pTBins{"pTBins", {200, 0.f, 10.f}, "pT binning"};
  ConfigurableAxis etaBins{"etaBins", {2, -1.f, 1.f}, "eta binning"};
  ConfigurableAxis etaBinsDauthers{"etaBinsDauthers", {2, -1.f, 1.f}, "eta binning"};
  ConfigurableAxis phiBins{"phiBins", {4, 0.f, 6.28f}, "phi binning"};

  HistogramRegistry registry{"K0sResolution"};

  void init(InitContext const&)
  {
    const AxisSpec mAxis{mBins, "#it{m} (GeV/#it{c}^{2})"};
    const AxisSpec pTAxis{pTBins, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec etaAxis{etaBins, "#eta"};
    const AxisSpec etaAxisPosD{etaBinsDauthers, "#eta pos."};
    const AxisSpec etaAxisNegD{etaBinsDauthers, "#eta neg."};
    const AxisSpec phiAxis{phiBins, "#phi"};

    int nProc = 0;
    if (doprocessData) {
      LOG(info) << "processData enabled";
      nProc++;
    }
    if (doprocessMC) {
      LOG(info) << "processMC enabled";
      nProc++;
    }

    registry.add("h2_masspT", "h2_masspT", {HistType::kTH2F, {mAxis, pTAxis}});
    registry.add("h2_masseta", "h2_masseta", {HistType::kTH2F, {mAxis, etaAxis}});
    registry.add("h2_massphi", "h2_massphi", {HistType::kTH2F, {mAxis, phiAxis}});
    if (!useMultidimHisto) {
      return;
    }
    registry.add("thn_mass", "thn_mass", kTHnSparseF, {mAxis, pTAxis, etaAxis, phiAxis, etaAxisPosD, etaAxisNegD});
  }

  // Selection criteria
  Configurable<double> v0setting_cospa{"v0setting_cospa", 0.995, "V0 CosPA"};
  Configurable<float> v0setting_dcav0dau{"v0setting_dcav0dau", 1., "DCA V0 Daughters"};
  Configurable<float> v0setting_dcapostopv{"v0setting_dcapostopv", 0.1, "DCA Pos To PV"};
  Configurable<float> v0setting_dcanegtopv{"v0setting_dcanegtopv", 0.1, "DCA Neg To PV"};
  Configurable<float> v0setting_radius{"v0setting_radius", 0.9, "V0 Radius"};
  Configurable<float> v0lifetime{"v0lifetime", 3., "n ctau"};
  Configurable<float> rapidity{"rapidity", 0.5, "rapidity"};
  Configurable<float> nSigTPC{"nSigTPC", 10., "nSigTPC"};
  Configurable<int> trdSelectionPos{"trdSelectionPos", 0, "Flag for the TRD selection on positive daughters: -1 no TRD, 0 no selection, 1 TRD"};
  Configurable<int> trdSelectionNeg{"trdSelectionNeg", 0, "Flag for the TRD selection on negative daughters: -1 no TRD, 0 no selection, 1 TRD"};
  Configurable<int> tofSelectionPos{"tofSelectionPos", 0, "Flag for the TOF selection on positive daughters: -1 no TOF, 0 no selection, 1 TOF"};
  Configurable<int> tofSelectionNeg{"tofSelectionNeg", 0, "Flag for the TOF selection on negative daughters: -1 no TOF, 0 no selection, 1 TOF"};
  Configurable<bool> eventSelection{"eventSelection", true, "event selection"};
  Configurable<bool> useMultidimHisto{"useMultidimHisto", false, "use multidimentional histograms"};

  template <typename T1, typename T2, typename C>
  bool acceptV0(const T1& v0, const T2& ntrack, const T2& ptrack, const C& collision)
  {
    // Apply selections on V0
    if (TMath::Abs(v0.yK0Short()) > rapidity) {
      return false;
    }
    if (v0.v0cosPA() < v0setting_cospa) {
      return false;
    }
    if (v0.v0radius() < v0setting_radius) {
      return false;
    }
    if (v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * pid_constants::sMasses[PID::K0] > 2.684 * v0lifetime) {
      return false;
    }

    // Apply selections on V0 daughters
    if (!ntrack.hasTPC() || !ptrack.hasTPC()) {
      return false;
    }
    if (ntrack.tpcNSigmaPi() > nSigTPC || ptrack.tpcNSigmaPi() > nSigTPC) {
      return false;
    }
    // TOF selection
    switch (tofSelectionPos) {
      case -1:
        if (ptrack.hasTOF()) {
          return false;
        }
        break;
      case 0:
        break;
      case 1:
        if (!ptrack.hasTOF()) {
          return false;
        }
        break;
      default:
        LOG(fatal) << "Invalid TRD selection for positive daughter";
        break;
    }
    switch (tofSelectionNeg) {
      case -1:
        if (ntrack.hasTOF()) {
          return false;
        }
        break;
      case 0:
        break;
      case 1:
        if (!ntrack.hasTOF()) {
          return false;
        }
        break;
      default:
        LOG(fatal) << "Invalid TRD selection for negative daughter";
        break;
    }
    // TRD selection
    switch (trdSelectionPos) {
      case -1:
        if (ptrack.hasTRD()) {
          return false;
        }
        break;
      case 0:
        break;
      case 1:
        if (!ptrack.hasTRD()) {
          return false;
        }
        break;
      default:
        LOG(fatal) << "Invalid TRD selection for positive daughter";
        break;
    }
    switch (trdSelectionNeg) {
      case -1:
        if (ntrack.hasTRD()) {
          return false;
        }
        break;
      case 0:
        break;
      case 1:
        if (!ntrack.hasTRD()) {
          return false;
        }
        break;
      default:
        LOG(fatal) << "Invalid TRD selection for negative daughter";
        break;
    }
    return true;
  }

  Filter v0Filter = nabs(aod::v0data::dcapostopv) > v0setting_dcapostopv&& nabs(aod::v0data::dcanegtopv) > v0setting_dcanegtopv&& aod::v0data::dcaV0daughters < v0setting_dcav0dau;

  void processData(SelectedCollisions::iterator const& collision, soa::Filtered<aod::V0Datas> const& fullV0s, PIDTracks const& tracks)
  {
    if (eventSelection && !collision.sel8())
      return;

    for (auto& v0 : fullV0s) {

      const auto& posTrack = v0.posTrack_as<PIDTracks>();
      const auto& negTrack = v0.negTrack_as<PIDTracks>();
      if (!acceptV0(v0, negTrack, posTrack, collision))
        continue;

      registry.fill(HIST("h2_masspT"), v0.mK0Short(), v0.pt());
      registry.fill(HIST("h2_masseta"), v0.mK0Short(), v0.eta());
      registry.fill(HIST("h2_massphi"), v0.mK0Short(), v0.phi());
      if (useMultidimHisto) {
        registry.fill(HIST("thn_mass"), v0.mK0Short(), v0.pt(), v0.eta(), v0.phi(), posTrack.eta(), negTrack.eta());
      }
    }
  }
  PROCESS_SWITCH(perfK0sResolution, processData, "Process data", true);

  void processMC(SelectedCollisions::iterator const& collision, soa::Filtered<aod::V0Datas> const& fullV0s, soa::Join<PIDTracks, aod::McTrackLabels> const& tracks, aod::McParticles const&)
  {
    if (eventSelection && !collision.sel8())
      return;

    for (auto& v0 : fullV0s) {

      const auto& posTrack = v0.posTrack_as<soa::Join<PIDTracks, aod::McTrackLabels>>();
      const auto& negTrack = v0.negTrack_as<soa::Join<PIDTracks, aod::McTrackLabels>>();
      if (!acceptV0(v0, negTrack, posTrack, collision))
        continue;
      if (!posTrack.has_mcParticle()) {
        continue;
      }
      if (!negTrack.has_mcParticle()) {
        continue;
      }
      if (posTrack.mcParticle().pdgCode() != PID::Pion || negTrack.mcParticle().pdgCode() != PID::Pion) {
        continue;
      }
      registry.fill(HIST("h2_masspT"), v0.mK0Short(), v0.pt());
      registry.fill(HIST("h2_masseta"), v0.mK0Short(), v0.eta());
      registry.fill(HIST("h2_massphi"), v0.mK0Short(), v0.phi());
      if (useMultidimHisto) {
        registry.fill(HIST("thn_mass"), v0.mK0Short(), v0.pt(), v0.eta(), v0.phi(), posTrack.eta(), negTrack.eta());
      }
    }
  }
  PROCESS_SWITCH(perfK0sResolution, processMC, "Process MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<perfK0sResolution>(cfgc)}; }
