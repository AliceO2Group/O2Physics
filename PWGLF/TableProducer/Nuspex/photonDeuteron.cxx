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
//
/// \file photonDeuteron.cxx
/// \brief Task to analyze photon-deuteron correlations
/// \author Arvind Khuntia <arvind.khuntia@cern.ch> and Francesco Noferini <francesco.noferini@cern.ch>

#include "PWGLF/DataModel/LFPhotonDeuteronTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include <TLorentzVector.h>
#include <TVector3.h>

#include <iostream>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct PhotonDeuteronCorrelation {
  // Histogram registry
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // Produce output table
  Produces<aod::PhotonDeuteronPairs> photonDeuteronTable;

  // Configurable parameters
  Configurable<float> cfgPhotonMassWindow{"cfgPhotonMassWindow", 0.05, "Photon mass window (GeV/c2)"};
  Configurable<float> cfgPhotonMinPt{"cfgPhotonMinPt", 0.05, "Minimum photon pT (GeV/c)"};
  Configurable<float> cfgPhotonMaxPt{"cfgPhotonMaxPt", 10.0, "Maximum photon pT (GeV/c)"};
  Configurable<float> cfgPhotonMaxEta{"cfgPhotonMaxEta", 0.8, "Maximum photon eta"};
  Configurable<float> cfgPhotonDaughterNSigmaElTPCMax{"cfgPhotonDaughterNSigmaElTPCMax", 5.0, "Maximum NSigma TPC electron for photon daughters (999=off)"};

  Configurable<float> cfgDeuteronMinPt{"cfgDeuteronMinPt", 0.2, "Minimum deuteron pT (GeV/c)"};
  Configurable<float> cfgDeuteronMaxPt{"cfgDeuteronMaxPt", 5.0, "Maximum deuteron pT (GeV/c)"};
  Configurable<float> cfgDeuteronMaxEta{"cfgDeuteronMaxEta", 0.8, "Maximum deuteron eta"};
  Configurable<float> cfgDeuteronNSigmaTPCMax{"cfgDeuteronNSigmaTPCMax", 3.0, "Maximum NSigma TPC for deuteron"};
  Configurable<float> cfgDeuteronNSigmaTOFMax{"cfgDeuteronNSigmaTOFMax", 3.0, "Maximum NSigma TOF for deuteron"};

  Configurable<float> cfgV0CosPA{"cfgV0CosPA", 0.995, "Minimum V0 cosine of pointing angle"};
  Configurable<float> cfgV0Radius{"cfgV0Radius", 5.0, "Minimum V0 radius (cm)"};
  Configurable<bool> cfgUsePhotonDaughterPIDTPCOnly{"cfgUsePhotonDaughterPIDTPCOnly", true, "Use TPC-only PID for photon daughters"};

  // Particle masses (in GeV/c²)
  static constexpr float massProton = o2::constants::physics::MassProton;
  static constexpr float massNeutron = o2::constants::physics::MassNeutron;
  static constexpr float massDeuteron = o2::constants::physics::MassDeuteron;

  // Initialize histograms
  void init(InitContext const&)
  {
    const AxisSpec axisPhotonMass{200, 0.0, 0.2, "M_{ee} (GeV/c^{2})"};
    const AxisSpec axisPt{100, 0.0, 10.0, "p_{T} (GeV/c)"};
    const AxisSpec axisEta{100, -1.0, 1.0, "#eta"};
    const AxisSpec axisPhi{64, 0.0, 2.0 * M_PI, "#phi (rad)"};
    const AxisSpec axisDeltaPhi{64, -0.5 * M_PI, 1.5 * M_PI, "#Delta#phi (rad)"};
    const AxisSpec axisDeltaEta{100, -2.0, 2.0, "#Delta#eta"};
    const AxisSpec axisNSigma{200, -10.0, 10.0, "n#sigma"};
    const AxisSpec axisV0CosPA{100, 0.9, 1.0, "cos(PA)"};
    const AxisSpec axisV0Radius{100, 0.0, 50.0, "R (cm)"};
    const AxisSpec axisRelativeMomentum{200, 0.0, 1.0, "k^{*}_{pn} (GeV/c)"};

    // Event histograms
    histos.add("hEventCounter", "Event counter", kTH1F, {{5, 0.5, 5.5}});

    // V0 (photon) histograms
    histos.add("hV0Mass", "V0 invariant mass", kTH1F, {axisPhotonMass});
    histos.add("hV0Pt", "V0 p_{T}", kTH1F, {axisPt});
    histos.add("hV0Eta", "V0 #eta", kTH1F, {axisEta});
    histos.add("hV0Phi", "V0 #phi", kTH1F, {axisPhi});
    histos.add("hV0CosPA", "V0 cos(PA)", kTH1F, {axisV0CosPA});
    histos.add("hV0Radius", "V0 radius", kTH1F, {axisV0Radius});
    histos.add("hPhotonPtEta", "Photon p_{T} vs #eta", kTH2F, {axisPt, axisEta});

    // Deuteron histograms
    histos.add("hDeuteronPt", "Deuteron p_{T}", kTH1F, {axisPt});
    histos.add("hDeuteronEta", "Deuteron #eta", kTH1F, {axisEta});
    histos.add("hDeuteronPhi", "Deuteron #phi", kTH1F, {axisPhi});
    histos.add("hDeuteronNSigmaTPC", "Deuteron n#sigma TPC", kTH1F, {axisNSigma});
    histos.add("hDeuteronNSigmaTOF", "Deuteron n#sigma TOF", kTH1F, {axisNSigma});
    histos.add("hDeuteronNSigmaTPCvsPt", "Deuteron n#sigma TPC vs p_{T}", kTH2F, {axisPt, axisNSigma});
    histos.add("hDeuteronNSigmaTOFvsPt", "Deuteron n#sigma TOF vs p_{T}", kTH2F, {axisPt, axisNSigma});

    // Correlation histograms
    histos.add("hPhotonDeuteronDeltaPhi", "Photon-Deuteron #Delta#phi", kTH1F, {axisDeltaPhi});
    histos.add("hPhotonDeuteronDeltaEta", "Photon-Deuteron #Delta#eta", kTH1F, {axisDeltaEta});
    histos.add("hPhotonDeuteronCorrelation", "Photon-Deuteron correlation", kTH2F, {axisDeltaPhi, axisDeltaEta});
    histos.add("hPhotonDeuteronInvMass", "Photon-Deuteron invariant mass", kTH1F, {{200, 0.0, 5.0, "M_{#gamma d} (GeV/c^{2})"}});
    histos.add("hPhotonDeuteronPtCorr", "Photon p_{T} vs Deuteron p_{T}", kTH2F, {axisPt, axisPt});
    histos.add("hRelativeMomentum", "Relative momentum k^{*}_{pn}", kTH1F, {axisRelativeMomentum});
    histos.add("hRelativeMomentumVsInvMass", "k^{*}_{pn} vs invariant mass", kTH2F, {{200, 0.0, 5.0, "M_{#gamma d} (GeV/c^{2})"}, {200, 0.0, 1.0, "k^{*}_{pn} (GeV/c)"}});
  }

  // V0 photon selection
  template <typename TV0, typename TTracks>
  bool selectPhoton(TV0 const& v0, TTracks const& tracks)
  {
    // Basic V0 quality cuts
    if (v0.v0cosPA() < cfgV0CosPA) {
      return false;
    }
    if (v0.v0radius() < cfgV0Radius) {
      return false;
    }

    // Photon mass window
    if (std::abs(v0.mGamma()) > cfgPhotonMassWindow) {
      return false;
    }

    // Kinematic cuts
    if (v0.pt() < cfgPhotonMinPt || v0.pt() > cfgPhotonMaxPt) {
      return false;
    }
    if (std::abs(v0.eta()) > cfgPhotonMaxEta) {
      return false;
    }

    // Optional electron PID cuts for daughter tracks (TPC-only)
    if (cfgUsePhotonDaughterPIDTPCOnly) {
      auto posTrack = v0.template posTrack_as<TTracks>();
      auto negTrack = v0.template negTrack_as<TTracks>();

      if (std::abs(posTrack.tpcNSigmaEl()) > cfgPhotonDaughterNSigmaElTPCMax) {
        return false;
      }
      if (std::abs(negTrack.tpcNSigmaEl()) > cfgPhotonDaughterNSigmaElTPCMax) {
        return false;
      }
    }

    return true;
  }

  // Deuteron PID selection
  template <typename TTrack>
  bool selectDeuteron(TTrack const& track)
  {
    // Kinematic cuts
    if (track.pt() < cfgDeuteronMinPt || track.pt() > cfgDeuteronMaxPt) {
      return false;
    }
    if (std::abs(track.eta()) > cfgDeuteronMaxEta) {
      return false;
    }

    // TPC PID
    if (std::abs(track.tpcNSigmaDe()) > cfgDeuteronNSigmaTPCMax) {
      return false;
    }

    // TOF PID (if available)
    if (track.hasTOF() && std::abs(track.tofNSigmaDe()) > cfgDeuteronNSigmaTOFMax) {
      return false;
    }

    return true;
  }

  // range [-pi/2, 3pi/2]
  float getDeltaPhi(float phi1, float phi2)
  {
    float dphi = phi1 - phi2;
    if (dphi > 1.5 * M_PI) {
      dphi -= 2.0 * M_PI;
    }
    if (dphi < -0.5 * M_PI) {
      dphi += 2.0 * M_PI;
    }
    return dphi;
  }

  // Calculate relative momentum k*_pn from the photon-deuteron invariant mass [Eq. 4.6]
  float calculateRelativeMomentum(float invMass)
  {
    if (invMass <= 0.0) {
      return -1.0;
    }

    float M = invMass;
    float M2 = M * M;
    float M4 = M2 * M2;
    float mn2 = massNeutron * massNeutron;
    float mp2 = massProton * massProton;
    float deltaMass2 = mn2 - mp2;
    float sumMass2 = mn2 + mp2;

    float term = M4 + deltaMass2 * deltaMass2 - 2.0 * sumMass2 * M2;

    if (term < 0.0) {
      return -1.0;
    }

    float kpn = 0.5 / M * std::sqrt(term);
    return kpn;
  }

  // Process function for V0s and tracks
  void process(aod::Collision const& collision,
               soa::Join<aod::V0Indices, aod::V0CoresBase, aod::V0TrackXs> const& V0s,
               soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection,
                         aod::pidTPCFullDe, aod::pidTOFFullDe,
                         aod::pidTPCFullEl, aod::pidTOFFullEl> const& tracks)
  {
    histos.fill(HIST("hEventCounter"), 1);

    // Loop over V0s to find photons
    std::vector<int> photonIndices;
    for (auto const& v0 : V0s) {
      // Fill V0 QA histograms
      histos.fill(HIST("hV0Mass"), v0.mGamma());
      histos.fill(HIST("hV0CosPA"), v0.v0cosPA());
      histos.fill(HIST("hV0Radius"), v0.v0radius());

      // Select photons
      if (!selectPhoton(v0, tracks)) {
        continue;
      }

      // Fill photon histograms
      histos.fill(HIST("hV0Pt"), v0.pt());
      histos.fill(HIST("hV0Eta"), v0.eta());
      histos.fill(HIST("hV0Phi"), v0.phi());
      histos.fill(HIST("hPhotonPtEta"), v0.pt(), v0.eta());

      if (v0.isPhotonTPConly())
        photonIndices.push_back(v0.index());
      if (v0.isPhotonTPConly())
        std::cout << " [main] global index photon: " << v0.globalIndex() << " v0 id: " << v0.index() << " pt " << v0.pt() << std::endl;
    }

    // Loop over tracks to find deuterons
    std::vector<int> deuteronIndices;
    for (auto const& track : tracks) {
      // Basic track quality
      if (!track.isGlobalTrack()) {
        continue;
      }

      // Select deuterons
      if (!selectDeuteron(track)) {
        continue;
      }

      // Fill deuteron histograms
      histos.fill(HIST("hDeuteronPt"), track.pt());
      histos.fill(HIST("hDeuteronEta"), track.eta());
      histos.fill(HIST("hDeuteronPhi"), track.phi());
      histos.fill(HIST("hDeuteronNSigmaTPC"), track.tpcNSigmaDe());
      histos.fill(HIST("hDeuteronNSigmaTPCvsPt"), track.pt(), track.tpcNSigmaDe());

      if (track.hasTOF()) {
        histos.fill(HIST("hDeuteronNSigmaTOF"), track.tofNSigmaDe());
        histos.fill(HIST("hDeuteronNSigmaTOFvsPt"), track.pt(), track.tofNSigmaDe());
      }
      deuteronIndices.push_back(track.index());
    }
    // Calculate correlations between photons and deuterons
    for (auto const& photonIdx : photonIndices) {
      const auto& photon = V0s.iteratorAt(photonIdx);

      for (auto const& deuteronIdx : deuteronIndices) {
        const auto& deuteron = tracks.iteratorAt(deuteronIdx);

        // Calculate angular correlations
        float deltaPhi = getDeltaPhi(photon.phi(), deuteron.phi());
        float deltaEta = photon.eta() - deuteron.eta();

        // Fill correlation histograms
        histos.fill(HIST("hPhotonDeuteronDeltaPhi"), deltaPhi);
        histos.fill(HIST("hPhotonDeuteronDeltaEta"), deltaEta);
        histos.fill(HIST("hPhotonDeuteronCorrelation"), deltaPhi, deltaEta);
        histos.fill(HIST("hPhotonDeuteronPtCorr"), photon.pt(), deuteron.pt());

        // Calculate invariant mass
        TLorentzVector photonVec, deuteronVec;
        photonVec.SetPtEtaPhiM(photon.pt(), photon.eta(), photon.phi(), 0.0);                  // Photon-mass = 0
        deuteronVec.SetPtEtaPhiM(deuteron.pt(), deuteron.eta(), deuteron.phi(), massDeuteron); // Deuteron-mass

        TLorentzVector combinedVec = photonVec + deuteronVec;
        float invMass = combinedVec.M();
        histos.fill(HIST("hPhotonDeuteronInvMass"), invMass);

        // Calculate relative momentum using Equation 4.6
        float kpn = calculateRelativeMomentum(invMass);
        if (kpn >= 0.0) {
          histos.fill(HIST("hRelativeMomentum"), kpn);
          histos.fill(HIST("hRelativeMomentumVsInvMass"), invMass, kpn);
        }

        // Fill the output table
        auto posTrack = photon.posTrack_as<soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::pidTPCFullDe, aod::pidTOFFullDe, aod::pidTPCFullEl, aod::pidTOFFullEl>>();
        auto negTrack = photon.negTrack_as<soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::pidTPCFullDe, aod::pidTOFFullDe, aod::pidTPCFullEl, aod::pidTOFFullEl>>();

        photonDeuteronTable(
          photon.pt(),
          photon.eta(),
          photon.phi(),
          photon.mGamma(),
          posTrack.pt(),
          posTrack.eta(),
          posTrack.phi(),
          posTrack.tpcNSigmaEl(),
          negTrack.pt(),
          negTrack.eta(),
          negTrack.phi(),
          negTrack.tpcNSigmaEl(),
          deuteron.pt(),
          deuteron.eta(),
          deuteron.phi(),
          deuteron.tpcNSigmaDe(),
          deuteron.hasTOF() ? deuteron.tofNSigmaDe() : -999.0f,
          deltaPhi,
          deltaEta,
          invMass,
          kpn);
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<PhotonDeuteronCorrelation>(cfgc)};
}
