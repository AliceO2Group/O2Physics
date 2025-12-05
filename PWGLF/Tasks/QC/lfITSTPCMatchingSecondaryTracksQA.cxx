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
///
/// \file lfITSTPCMatchingSecondaryTracksQA.cxx
///
/// \brief task for QA of ITS-TPC matching efficiency of secondary tracks from V0s
/// \author Alberto Caliva (alberto.caliva@cern.ch), Francesca Ercolessi (francesca.ercolessi@cern.ch), Nicolò Jacazio (nicolo.jacazio@cern.ch)
/// \since Feb 11, 2025

#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

#include <TLorentzVector.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TPDGCode.h>
#include <TVector2.h>
#include <TVector3.h>

#include <cmath>
#include <string>
#include <vector>

using namespace std;
using namespace o2;
using namespace o2::soa;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;
using namespace o2::constants::math;
using std::array;

using SelCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms>;
using SimCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::McCollisionLabels>;
using StrHadronDaughterTracks = soa::Join<aod::Tracks, aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::TracksDCA, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr>;
using MCTracks = soa::Join<StrHadronDaughterTracks, aod::McTrackLabels>;

struct LfITSTPCMatchingSecondaryTracksQA {

  HistogramRegistry registryData{"registryData", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry registryMC{"registryMC", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  // Global Parameters
  Configurable<double> zVtx{"zVtx", 10.0, "Maximum zVertex"};

  // Track Parameters
  Configurable<float> minITSnCls{"minITSnCls", 1.0f, "min number of ITS clusters"};
  Configurable<float> minNCrossedRowsTPC{"minNCrossedRowsTPC", 80.0f, "min number of TPC crossed rows"};
  Configurable<float> maxChi2TPC{"maxChi2TPC", 4.0f, "max chi2 per cluster TPC"};
  Configurable<float> maxChi2ITS{"maxChi2ITS", 36.0f, "max chi2 per cluster ITS"};
  Configurable<float> etaMin{"etaMin", -0.8f, "eta min"};
  Configurable<float> etaMax{"etaMax", +0.8f, "eta max"};
  Configurable<float> nsigmaTPCmin{"nsigmaTPCmin", -3.0f, "Minimum nsigma TPC"};
  Configurable<float> nsigmaTPCmax{"nsigmaTPCmax", +3.0f, "Maximum nsigma TPC"};
  Configurable<float> nsigmaTOFmin{"nsigmaTOFmin", -3.0f, "Minimum nsigma TOF"};
  Configurable<float> nsigmaTOFmax{"nsigmaTOFmax", +3.0f, "Maximum nsigma TOF"};
  Configurable<float> dcaxyMax{"dcaxyMax", 0.1f, "dcaxy max"};
  Configurable<float> dcazMax{"dcazMax", 0.1f, "dcaz max"};
  Configurable<float> dcaMin{"dcaMin", 0.1f, "dca min"};
  Configurable<bool> requireTOF{"requireTOF", false, "require TOF hit"};
  Configurable<bool> requireItsHits{"requireItsHits", false, "require ITS hits"};
  Configurable<std::vector<float>> requiredHit{"requiredHit", {0, 0, 0, 0, 0, 0, 0}, "required ITS Hits (1=required, 0=not required)"};

  // V0 Parameters
  Configurable<float> minimumV0Radius{"minimumV0Radius", 0.0f, "Minimum V0 Radius"};
  Configurable<float> maximumV0Radius{"maximumV0Radius", 100.0f, "Maximum V0 Radius"};
  Configurable<float> dcanegtoPVmin{"dcanegtoPVmin", 0.1f, "Minimum DCA Neg To PV"};
  Configurable<float> dcapostoPVmin{"dcapostoPVmin", 0.1f, "Minimum DCA Pos To PV"};
  Configurable<float> v0cospaMin{"v0cospaMin", 0.99f, "Minimum V0 CosPA"};
  Configurable<float> dcaV0DaughtersMax{"dcaV0DaughtersMax", 0.5f, "Maximum DCA Daughters"};
  Configurable<float> mK0Min{"mK0Min", 0.48f, "K0 mass lower cut"};
  Configurable<float> mK0Max{"mK0Max", 0.52f, "K0 mass upper cut"};

  void init(InitContext const&)
  {
    // Event Counters
    if (doprocessData) {
      registryData.add("number_of_events_data", "number of events in data", HistType::kTH1D, {{20, 0, 20, "Event Cuts"}});
      registryData.add("dcaxyDatavspt", "dcaxyDatavspt", HistType::kTH2D, {{100, 0, 10, "#it{p}_{T} (GeV/#it{c})"}, {400, -2, 2, "DCA_{xy} (cm)"}});
      registryData.add("dcazDatavspt", "dcazDatavspt", HistType::kTH2D, {{100, 0, 10, "#it{p}_{T} (GeV/#it{c})"}, {400, -2, 2, "DCA_{z} (cm)"}});
      registryData.add("primPionTPC", "primPionTPC", HistType::kTH3D, {{100, 0, 10, "#it{p}_{T} (GeV/#it{c})"}, {16, -0.8, 0.8, "#eta"}, {100, 0, TwoPI, "#phi"}});
      registryData.add("primPionTPC_ITS", "primPionTPC_ITS", HistType::kTH3D, {{100, 0, 10, "#it{p}_{T} (GeV/#it{c})"}, {16, -0.8, 0.8, "#eta"}, {100, 0, TwoPI, "#phi"}});
      registryData.add("secPionTPC", "secPionTPC", HistType::kTH3D, {{100, 0, 10, "#it{p}_{T} (GeV/#it{c})"}, {16, -0.8, 0.8, "#eta"}, {100, 0, TwoPI, "#phi"}});
      registryData.add("secPionTPC_ITS", "secPionTPC_ITS", HistType::kTH3D, {{100, 0, 10, "#it{p}_{T} (GeV/#it{c})"}, {16, -0.8, 0.8, "#eta"}, {100, 0, TwoPI, "#phi"}});
      registryData.add("secPionV0TPC", "secPionV0TPC", HistType::kTH3D, {{100, 0, 10, "#it{p}_{T} (GeV/#it{c})"}, {16, -0.8, 0.8, "#eta"}, {100, 0, TwoPI, "#phi"}});
      registryData.add("secPionV0TPC_ITS", "secPionV0TPC_ITS", HistType::kTH3D, {{100, 0, 10, "#it{p}_{T} (GeV/#it{c})"}, {16, -0.8, 0.8, "#eta"}, {100, 0, TwoPI, "#phi"}});
    }

    if (doprocessMC) {
      registryMC.add("number_of_events_mc", "number of events in mc", HistType::kTH1D, {{20, 0, 20, "Event Cuts"}});
      registryMC.add("dcaxyMCvspt", "dcaxyMCvspt", HistType::kTH2D, {{100, 0, 10, "#it{p}_{T} (GeV/#it{c})"}, {400, -2, 2, "DCA_{xy} (cm)"}});
      registryMC.add("dcazMCvspt", "dcazMCvspt", HistType::kTH2D, {{100, 0, 10, "#it{p}_{T} (GeV/#it{c})"}, {400, -2, 2, "DCA_{z} (cm)"}});
      registryMC.add("primPionTPC_MC", "primPionTPC_MC", HistType::kTH3D, {{100, 0, 10, "#it{p}_{T} (GeV/#it{c})"}, {16, -0.8, 0.8, "#eta"}, {100, 0, TwoPI, "#phi"}});
      registryMC.add("primPionTPC_ITS_MC", "primPionTPC_ITS_MC", HistType::kTH3D, {{100, 0, 10, "#it{p}_{T} (GeV/#it{c})"}, {16, -0.8, 0.8, "#eta"}, {100, 0, TwoPI, "#phi"}});
      registryMC.add("secPionTPC_MC", "secPionTPC_MC", HistType::kTH3D, {{100, 0, 10, "#it{p}_{T} (GeV/#it{c})"}, {16, -0.8, 0.8, "#eta"}, {100, 0, TwoPI, "#phi"}});
      registryMC.add("secPionTPC_ITS_MC", "secPionTPC_ITS_MC", HistType::kTH3D, {{100, 0, 10, "#it{p}_{T} (GeV/#it{c})"}, {16, -0.8, 0.8, "#eta"}, {100, 0, TwoPI, "#phi"}});
      registryMC.add("secPionV0TPC_MC", "secPionV0TPC_MC", HistType::kTH3D, {{100, 0, 10, "#it{p}_{T} (GeV/#it{c})"}, {16, -0.8, 0.8, "#eta"}, {100, 0, TwoPI, "#phi"}});
      registryMC.add("secPionV0TPC_ITS_MC", "secPionV0TPC_ITS_MC", HistType::kTH3D, {{100, 0, 10, "#it{p}_{T} (GeV/#it{c})"}, {16, -0.8, 0.8, "#eta"}, {100, 0, TwoPI, "#phi"}});
    }
  }

  bool hasHitOnITSlayer(uint8_t itsClsmap, int layer)
  {
    unsigned char testBit = 1 << layer;
    return (itsClsmap & testBit);
  }

  template <typename TpcPrimTrack>
  bool passedTrackSelectionTpcPrimary(const TpcPrimTrack& track)
  {
    if (!track.hasTPC())
      return false;
    if (track.tpcNClsCrossedRows() < minNCrossedRowsTPC)
      return false;
    if (track.tpcChi2NCl() > maxChi2TPC)
      return false;
    if (track.eta() < etaMin || track.eta() > etaMax)
      return false;
    if (std::fabs(track.dcaXY()) > dcaxyMax)
      return false;
    if (std::fabs(track.dcaZ()) > dcazMax)
      return false;
    return true;
  }

  template <typename TpcSecTrack>
  bool passedTrackSelectionTpcSecondary(const TpcSecTrack& track)
  {
    if (!track.hasTPC())
      return false;
    if (track.tpcNClsCrossedRows() < minNCrossedRowsTPC)
      return false;
    if (track.tpcChi2NCl() > maxChi2TPC)
      return false;
    if (track.eta() < etaMin || track.eta() > etaMax)
      return false;
    if (std::sqrt(track.dcaXY() * track.dcaXY() + track.dcaZ() * track.dcaZ()) < dcaMin)
      return false;
    return true;
  }

  template <typename v0Track>
  bool passedTrackSelectionV0daughTPC(const v0Track& track)
  {
    if (!track.hasTPC())
      return false;
    if (track.tpcNClsCrossedRows() < minNCrossedRowsTPC)
      return false;
    if (track.tpcChi2NCl() > maxChi2TPC)
      return false;
    if (track.eta() < etaMin || track.eta() > etaMax)
      return false;
    return true;
  }

  // K0s Selections
  template <typename K0short>
  bool passedK0ShortSelection(const K0short& v0)
  {
    if (v0.v0cosPA() < v0cospaMin)
      return false;
    if (v0.v0radius() < minimumV0Radius || v0.v0radius() > maximumV0Radius)
      return false;
    if (v0.dcaV0daughters() > dcaV0DaughtersMax)
      return false;
    if (std::fabs(v0.dcapostopv()) < dcapostoPVmin)
      return false;
    if (std::fabs(v0.dcanegtopv()) < dcanegtoPVmin)
      return false;
    if (v0.mK0Short() < mK0Min || v0.mK0Short() > mK0Max)
      return false;

    return true;
  }

  template <typename pionTrack>
  bool passedPionSelection(const pionTrack& track)
  {
    // TPC Selection
    if (track.tpcNSigmaPi() < nsigmaTPCmin || track.tpcNSigmaPi() > nsigmaTPCmax)
      return false;

    // TOF Selection
    if (requireTOF) {
      if (track.tofNSigmaPi() < nsigmaTOFmin || track.tofNSigmaPi() > nsigmaTOFmax)
        return false;
    }
    return true;
  }

  template <typename ItsTrack>
  bool passedTrackSelectionIts(const ItsTrack& track)
  {
    /*
    if (!track.hasITS())
      return false;
    if (track.itsNCls() < minITSnCls)
      return false;
    if (track.itsChi2NCl() > maxChi2ITS)
      return false;
    */

    if (track.itsNCls() < minITSnCls)
      return false;

    auto requiredItsHit = static_cast<std::vector<float>>(requiredHit);
    if (requireItsHits) {
      for (int i = 0; i < 7; i++) {
        if (requiredItsHit[i] > 0 && !hasHitOnITSlayer(track.itsClusterMap(), i)) {
          return false;
        }
      }
    }
    return true;
  }

  void processData(SelCollisions::iterator const& collision, aod::V0Datas const& fullV0s, StrHadronDaughterTracks const& tracks)
  {
    registryData.fill(HIST("number_of_events_data"), 0.5);

    // Event Selection
    if (!collision.sel8() || std::fabs(collision.posZ()) > zVtx)
      return;
    registryData.fill(HIST("number_of_events_data"), 1.5);

    for (const auto& track : tracks) {

      // DCA distributions
      if (passedTrackSelectionV0daughTPC(track) && passedPionSelection(track)) {
        registryData.fill(HIST("dcaxyDatavspt"), track.pt(), track.dcaXY());
        registryData.fill(HIST("dcazDatavspt"), track.pt(), track.dcaZ());
      }

      // Primary Tracks
      if (passedTrackSelectionTpcPrimary(track) && passedPionSelection(track))
        registryData.fill(HIST("primPionTPC"), track.pt(), track.eta(), TVector2::Phi_0_2pi(track.phi()));
      if (passedTrackSelectionTpcPrimary(track) && passedPionSelection(track) && passedTrackSelectionIts(track))
        registryData.fill(HIST("primPionTPC_ITS"), track.pt(), track.eta(), TVector2::Phi_0_2pi(track.phi()));

      // Secondary Tracks
      if (passedTrackSelectionTpcSecondary(track) && passedPionSelection(track))
        registryData.fill(HIST("secPionTPC"), track.pt(), track.eta(), TVector2::Phi_0_2pi(track.phi()));
      if (passedTrackSelectionTpcSecondary(track) && passedPionSelection(track) && passedTrackSelectionIts(track))
        registryData.fill(HIST("secPionTPC_ITS"), track.pt(), track.eta(), TVector2::Phi_0_2pi(track.phi()));
    }

    for (const auto& v0 : fullV0s) {

      const auto& posTrack = v0.posTrack_as<StrHadronDaughterTracks>();
      const auto& negTrack = v0.negTrack_as<StrHadronDaughterTracks>();
      if (!passedK0ShortSelection(v0))
        continue;

      if (passedTrackSelectionV0daughTPC(posTrack) && passedPionSelection(posTrack))
        registryData.fill(HIST("secPionV0TPC"), posTrack.pt(), posTrack.eta(), TVector2::Phi_0_2pi(posTrack.phi()));
      if (passedTrackSelectionV0daughTPC(negTrack) && passedPionSelection(negTrack))
        registryData.fill(HIST("secPionV0TPC"), negTrack.pt(), negTrack.eta(), TVector2::Phi_0_2pi(negTrack.phi()));
      if (passedTrackSelectionV0daughTPC(posTrack) && passedPionSelection(posTrack) && passedTrackSelectionIts(posTrack))
        registryData.fill(HIST("secPionV0TPC_ITS"), posTrack.pt(), posTrack.eta(), TVector2::Phi_0_2pi(posTrack.phi()));
      if (passedTrackSelectionV0daughTPC(negTrack) && passedPionSelection(negTrack) && passedTrackSelectionIts(negTrack))
        registryData.fill(HIST("secPionV0TPC_ITS"), negTrack.pt(), negTrack.eta(), TVector2::Phi_0_2pi(negTrack.phi()));
    }
  }
  PROCESS_SWITCH(LfITSTPCMatchingSecondaryTracksQA, processData, "Process data", true);

  Preslice<aod::V0Datas> perCollisionV0 = o2::aod::v0data::collisionId;
  Preslice<MCTracks> perCollisionTrk = o2::aod::track::collisionId;

  void processMC(SimCollisions const& collisions, MCTracks const& mcTracks, aod::V0Datas const& fullV0s, const aod::McParticles&)
  {
    for (const auto& collision : collisions) {
      registryMC.fill(HIST("number_of_events_mc"), 0.5);

      if (!collision.sel8() || std::fabs(collision.posZ()) > zVtx)
        continue;
      registryMC.fill(HIST("number_of_events_mc"), 1.5);

      auto v0sPerColl = fullV0s.sliceBy(perCollisionV0, collision.globalIndex());
      auto tracksPerColl = mcTracks.sliceBy(perCollisionTrk, collision.globalIndex());

      for (const auto& track : tracksPerColl) {

        // DCA distributions
        if (passedTrackSelectionV0daughTPC(track) && passedPionSelection(track)) {
          registryMC.fill(HIST("dcaxyMCvspt"), track.pt(), track.dcaXY());
          registryMC.fill(HIST("dcazMCvspt"), track.pt(), track.dcaZ());
        }

        // Primary Tracks
        if (passedTrackSelectionTpcPrimary(track) && passedPionSelection(track))
          registryMC.fill(HIST("primPionTPC_MC"), track.pt(), track.eta(), TVector2::Phi_0_2pi(track.phi()));
        if (passedTrackSelectionTpcPrimary(track) && passedPionSelection(track) && passedTrackSelectionIts(track))
          registryMC.fill(HIST("primPionTPC_ITS_MC"), track.pt(), track.eta(), TVector2::Phi_0_2pi(track.phi()));

        // Secondary Tracks
        if (passedTrackSelectionTpcSecondary(track) && passedPionSelection(track))
          registryMC.fill(HIST("secPionTPC_MC"), track.pt(), track.eta(), TVector2::Phi_0_2pi(track.phi()));
        if (passedTrackSelectionTpcSecondary(track) && passedPionSelection(track) && passedTrackSelectionIts(track))
          registryMC.fill(HIST("secPionTPC_ITS_MC"), track.pt(), track.eta(), TVector2::Phi_0_2pi(track.phi()));
      }

      for (const auto& v0 : v0sPerColl) {

        const auto& posTrack = v0.posTrack_as<MCTracks>();
        const auto& negTrack = v0.negTrack_as<MCTracks>();
        if (!passedK0ShortSelection(v0))
          continue;

        if (passedTrackSelectionV0daughTPC(posTrack) && passedPionSelection(posTrack))
          registryMC.fill(HIST("secPionV0TPC_MC"), posTrack.pt(), posTrack.eta(), TVector2::Phi_0_2pi(posTrack.phi()));
        if (passedTrackSelectionV0daughTPC(negTrack) && passedPionSelection(negTrack))
          registryMC.fill(HIST("secPionV0TPC_MC"), negTrack.pt(), negTrack.eta(), TVector2::Phi_0_2pi(negTrack.phi()));
        if (passedTrackSelectionV0daughTPC(posTrack) && passedPionSelection(posTrack) && passedTrackSelectionIts(posTrack))
          registryMC.fill(HIST("secPionV0TPC_ITS_MC"), posTrack.pt(), posTrack.eta(), TVector2::Phi_0_2pi(posTrack.phi()));
        if (passedTrackSelectionV0daughTPC(negTrack) && passedPionSelection(negTrack) && passedTrackSelectionIts(negTrack))
          registryMC.fill(HIST("secPionV0TPC_ITS_MC"), negTrack.pt(), negTrack.eta(), TVector2::Phi_0_2pi(negTrack.phi()));
      }
    }
  }
  PROCESS_SWITCH(LfITSTPCMatchingSecondaryTracksQA, processMC, "Process MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<LfITSTPCMatchingSecondaryTracksQA>(cfgc)};
}
