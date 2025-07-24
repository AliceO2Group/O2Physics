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
/// \brief This task is to compute invariant mass of dimuon at forward rapidity for UPC events.
/// \author Anisa Khatun
/// \date 10.10.2024

// O2 headers
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"

// O2Physics headers
#include "PWGUD/DataModel/UDTables.h"
#include "PWGUD/Core/UDHelpers.h"
#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPLHCIFData.h"
#include "DataFormatsParameters/GRPECSObject.h"

// ROOT headers
#include "TSystem.h"
#include "TDatabasePDG.h"
#include "TLorentzVector.h"
#include "TLorentzVector.h"
#include "TMath.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace std;

struct UDTutorial06 {

  // Histogram registry: an object to hold your histograms
  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};
  Configurable<int> nBinsPt{"nBinsPt", 500, "N bins in pT histo"};
  Configurable<int> nBinsMass{"nBinsMass", 500, "N bins in InvMass histo"};

  using UDCollisionsFwd = soa::Join<o2::aod::UDCollisions, o2::aod::UDCollisionsSels, o2::aod::UDCollisionsSelsFwd>;
  using fwdtraks = soa::Join<aod::UDFwdTracks, aod::UDFwdTracksExtra>;

  void init(InitContext const&)
  {
    // define axes you want to use
    const AxisSpec axisCounter{1, 0, +1, ""};
    const AxisSpec axisEta{80, -5.0, -2.0, "#eta"};
    const AxisSpec axisPt{nBinsPt, 0.0, 10.0, "p_{T}"};
    const AxisSpec axisM{nBinsMass, 0.0, 10.0, "M_{#pi#pi}"};
    const AxisSpec axisPhi{120, -TMath::Pi(), -TMath::Pi(), "#phi"};
    const AxisSpec axisRapidity{80, -5.0, -2.0, "#it{y}"};

    // create histograms
    // Fill counter to see effect of each selection criteria
    auto hSelectionCounter = registry.add<TH1>("hSelectionCounter", "hSelectionCounter;;NEvents", HistType::kTH1I, {{15, 0., 15.}});
    TString SelectionCuts[13] = {"NoSelection", "V0A", "rAbs1", "rAbs2", "trackmatch", "eta1", "eta2", "trackpt", "pair rapidity", "mass_cut", "unlikesign", "likesign"};

    for (int i = 0; i < 12; i++) {
      hSelectionCounter->GetXaxis()->SetBinLabel(i + 1, SelectionCuts[i].Data());
    }

    registry.add("eventCounter", "eventCounter", kTH1F, {axisCounter});
    registry.add("hnumContrib", "hnumContrib;;#counts", kTH1D, {{100, -0.5, 99.5}});
    registry.add("hTracks1", "N_{tracks}", kTH1F, {{100, -0.5, 99.5}});
    registry.add("hTracks2", "N_{tracks}", kTH1F, {{100, -0.5, 99.5}});
    registry.add("hTracksMuons", "N_{tracks}", kTH1F, {{100, -0.5, 99.5}});

    registry.add("etaMuon1", "etaMuon1", kTH1F, {axisEta});
    registry.add("ptMuon1", "ptMuon1", kTH1F, {axisPt});
    registry.add("phiMuon1", "phiMuon1", kTH1F, {axisPhi});

    registry.add("etaMuon2", "etaMuon2", kTH1F, {axisEta});
    registry.add("ptMuon2", "ptMuon2", kTH1F, {axisPt});
    registry.add("phiMuon2", "phiMuon2", kTH1F, {axisPhi});

    registry.add("YJpsi", "YJpsi", kTH1F, {axisRapidity});
    registry.add("PhiJpsi", "PhiJpsi", kTH1F, {axisPhi});

    registry.add("MJpsiUnlike", "MJpsi", kTH1F, {axisM});
    registry.add("PtJpsiUnlike", "PtJpsi", kTH1F, {axisPt});
    registry.add("MJpsiLike", "MJpsi", kTH1F, {axisM});
    registry.add("PtJpsiLike", "PtJpsi", kTH1F, {axisPt});
  }

  //____________________________________________________________________________________________

  template <typename TTrack1, typename TTrack2>
  void processCandidate(UDCollisionsFwd::iterator const& collision, TTrack1& tr1, TTrack2& tr2)
  {
    registry.fill(HIST("eventCounter"), 0.5);
    registry.fill(HIST("hnumContrib"), collision.numContrib());
    registry.fill(HIST("hTracks1"), tr1.size());
    registry.fill(HIST("hTracks2"), tr2.size());
    registry.fill(HIST("hSelectionCounter"), 0);

    // V0A selection
    const auto& ampsV0A = collision.amplitudesV0A();
    const auto& ampsRelBCsV0A = collision.ampRelBCsV0A();
    for (unsigned int i = 0; i < ampsV0A.size(); ++i) {
      if (std::abs(ampsRelBCsV0A[i]) <= 1) {
        if (ampsV0A[i] > 100.)
          return;
      }
    }

    // T0A selection
    const auto& ampsT0A = collision.amplitudesT0A();
    const auto& ampsRelBCsT0A = collision.ampRelBCsT0A();
    for (unsigned int i = 0; i < ampsT0A.size(); ++i) {
      if (std::abs(ampsRelBCsT0A[i]) <= 1) {
        if (ampsT0A[i] > 100.)
          return;
      }
    }

    registry.fill(HIST("hSelectionCounter"), 1);

    // absorber end selection
    if (tr1.rAtAbsorberEnd() < 17.6 || tr1.rAtAbsorberEnd() > 89.5)
      return;
    registry.fill(HIST("hSelectionCounter"), 2);
    if (tr2.rAtAbsorberEnd() < 17.6 || tr2.rAtAbsorberEnd() > 89.5)
      return;
    registry.fill(HIST("hSelectionCounter"), 3);

    // MCH-MID match selection
    if ((tr1.chi2MatchMCHMID() < 0) || (tr2.chi2MatchMCHMID() < 0))
      return;
    registry.fill(HIST("hSelectionCounter"), 4);

    bool isUnlikeSign = (tr1.sign() + tr2.sign()) == 0;
    bool isLikeSign = ((tr1.sign() + tr2.sign()) != 0);

    // track selection
    TLorentzVector p1, p2;
    p1.SetXYZM(tr1.px(), tr1.py(), tr1.pz(), o2::constants::physics::MassMuon);
    p2.SetXYZM(tr2.px(), tr2.py(), tr2.pz(), o2::constants::physics::MassMuon);
    TLorentzVector p = p1 + p2;

    // eta cut on each track
    if (p1.Eta() < -4.0 || p1.Eta() > -2.5)
      return;
    registry.fill(HIST("hSelectionCounter"), 5);
    if (p2.Eta() < -4.0 || p2.Eta() > -2.5)
      return;
    registry.fill(HIST("hSelectionCounter"), 6);

    // pt cut on track
    if ((p1.Pt() < 0.0) || (p2.Pt() < 0.0))
      return;
    registry.fill(HIST("hSelectionCounter"), 7);

    if ((p.Rapidity() < -4.0) || (p.Rapidity() > -2.5))
      return;
    registry.fill(HIST("hSelectionCounter"), 8);

    // cuts on pair kinematics
    if (!(p.M() > 2.5 && p.M() < 5.0))
      return;
    registry.fill(HIST("hSelectionCounter"), 9);

    registry.fill(HIST("hTracksMuons"), collision.numContrib());
    registry.fill(HIST("ptMuon1"), p1.Pt());
    registry.fill(HIST("ptMuon2"), p2.Pt());
    registry.fill(HIST("etaMuon1"), p1.Eta());
    registry.fill(HIST("etaMuon2"), p2.Eta());
    registry.fill(HIST("phiMuon1"), p1.Phi());
    registry.fill(HIST("phiMuon2"), p2.Phi());
    registry.fill(HIST("YJpsi"), p.Rapidity());
    registry.fill(HIST("PhiJpsi"), p.Phi());

    if (isUnlikeSign) {
      registry.fill(HIST("hSelectionCounter"), 10);
      registry.fill(HIST("MJpsiUnlike"), p.M());
      registry.fill(HIST("PtJpsiUnlike"), p.Pt());
    }

    if (isLikeSign) {
      registry.fill(HIST("hSelectionCounter"), 11);
      registry.fill(HIST("MJpsiLike"), p.M());
      registry.fill(HIST("PtJpsiLike"), p.Pt());
    }
  }

  //____________________________________________________________________________________________

  // Template that collects all collision IDs and track per collision
  template <typename TTracks>
  void collectCandIDs(std::unordered_map<int32_t, std::vector<int32_t>>& tracksPerCand, TTracks& tracks)
  {
    for (const auto& tr : tracks) {
      int32_t candId = tr.udCollisionId();
      if (candId < 0) {
        continue;
      }
      tracksPerCand[candId].push_back(tr.globalIndex());
    }
  }

  //____________________________________________________________________________________________

  // process candidates with 2  forward tracks
  void process(UDCollisionsFwd const& eventCandidates,
               fwdtraks const& fwdTracks)
  {
    std::unordered_map<int32_t, std::vector<int32_t>> tracksPerCand;
    collectCandIDs(tracksPerCand, fwdTracks);

    // assuming that candidates have exatly 2 forward tracks
    for (const auto& item : tracksPerCand) {
      int32_t trId1 = item.second[0];
      int32_t trId2 = item.second[1];
      int32_t candID = item.first;
      const auto& collision = eventCandidates.iteratorAt(candID);
      const auto& tr1 = fwdTracks.iteratorAt(trId1);
      const auto& tr2 = fwdTracks.iteratorAt(trId2);
      processCandidate(collision, tr1, tr2);
    }
  }
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<UDTutorial06>(cfgc, TaskName{"udtutorial06"})};
}
