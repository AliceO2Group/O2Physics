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
// Preliminary QA analysis task for resonances
// (1) For Run3
// (2) Event and track selection need to be optimized
// (3) particle = 0 --> phi
// (4) particle = 1 --> kstar
// (6) 4 process function (a) Data same event (b) Data mixed event (c) MC generated (d) MC reconstructed

#include <TH1F.h>
#include <TDirectory.h>
#include <THn.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TFile.h>
#include <TH2F.h>
#include <TLorentzVector.h>
#include <TPDGCode.h>
#include <TDatabasePDG.h>
#include <CCDB/BasicCCDBManager.h>

#include <cmath>
#include <array>
#include <cstdlib>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/StepTHn.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Framework/ASoAHelpers.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;
struct resonanceqa {

  framework::Service<o2::ccdb::BasicCCDBManager> ccdb; /// Accessing the CCDB                                                                 

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  // events
  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};
  // track
  Configurable<float> cfgCutPT{"cfgCutPT", 0.2, "PT cut on daughter track"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8, "Eta cut on daughter track"};
  Configurable<float> cfgCutDCAxy{"cfgCutDCAxy", 2.0f, "DCAxy range for tracks"};
  Configurable<float> cfgCutDCAz{"cfgCutDCAz", 2.0f, "DCAz range for tracks"};
  Configurable<float> nsigmaCutTPC{"nsigmacutTPC", 3.0, "Value of the TPC Nsigma cut"};
  Configurable<float> nsigmaCutCombined{"nsigmaCutCombined", 3.0, "Value of the TOF Nsigma cut"};
  Configurable<int> cfgNoMixedEvents{"cfgNoMixedEvents", 5, "Number of mixed events per event"};
  // particle
  Configurable<int> cfgparticletype{"cfgparticletype", 0, "Resonance particle type: 0(kstar), 1(phi), 2(Lambdastar)"};
  // MC
  Configurable<bool> isMC{"isMC", false, "Run MC"};
  void init(o2::framework::InitContext&)
  {

    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    uint64_t now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    ccdb->setCreatedNotAfter(now);


    AxisSpec multAxis = {3000,0,3000, "mult"};
    AxisSpec ptAxiss = {100, 0.0f, 10.0f, "#it{p}_{T}  (GeV/#it{c})"};
    AxisSpec invmassAxis = {90, 0.6, 1.5, "{M}_{#pion #K^{0}} (GeV/#it{c}^2)"};
    AxisSpec invmassAxisphi = {300, 0.9, 1.2, "{M}_{KK} (GeV/#it{c}^2)"};


    histos.add("hCentrality", "Centrality distribution", kTH1F, {{500, 0.0, 500.0}});
    histos.add("hVtxZ", "Vertex distribution in Z;Z (cm)", kTH1F, {{400, -20.0, 20.0}});
    histos.add("hNcontributor", "Number of primary vertex contributor", kTH1F, {{1000, 0.0f, 1000.0f}});
    histos.add("hEta", "Eta distribution", kTH1F, {{200, -1.0f, 1.0f}});
    histos.add("hDcaxy", "Dcaxy distribution", kTH1F, {{200, -1.0f, 1.0f}});
    histos.add("hDcaz", "Dcaz distribution", kTH1F, {{200, -1.0f, 1.0f}});
    histos.add("hNsigmaKaonTPC", "NsigmaKaon TPC distribution", kTH1F, {{100, -10.0f, 10.0f}});
    histos.add("hNsigmaKaonTOF", "NsigmaKaon TOF distribution", kTH1F, {{100, -10.0f, 10.0f}});
    histos.add("hNsigmaPionTPC", "NsigmaPion TPC distribution", kTH1F, {{100, -10.0f, 10.0f}});
    histos.add("hNsigmaPionTOF", "NsigmaPion TOF distribution", kTH1F, {{100, -10.0f, 10.0f}});
    if (cfgparticletype == 1 && !isMC) {
      histos.add("h3PhiInvMassUnlikeSign", "Invariant mass of Phi meson Unlike Sign", kTHnSparseD, {multAxis,ptAxiss,invmassAxisphi});
      histos.add("h3PhiInvMassLikeSignPP", "Invariant mass of Phi meson Like Sign positive", kTHnSparseD, {multAxis,ptAxiss,invmassAxisphi});
      histos.add("h3PhiInvMassLikeSignMM", "Invariant mass of Phi meson Like Sign negative", kTHnSparseD, {multAxis,ptAxiss,invmassAxisphi});
      histos.add("h3PhiInvMassRotational", "Invariant mass of Phi meson Rotational", kTHnSparseD, {multAxis,ptAxiss,invmassAxisphi});
      histos.add("h3PhiInvMassMixed", "Invariant mass of Phi meson Mixed", kTHnSparseD, {multAxis,ptAxiss,invmassAxisphi});
    } else if (cfgparticletype == 0 && !isMC) {
      histos.add("h3KstarInvMassUnlikeSign", "Invariant mass of Kstar meson Unlike Sign", kTHnSparseD, {multAxis,ptAxiss,invmassAxis});
      histos.add("h3KstarInvMassLikeSignPP", "Invariant mass of Kstar meson Like Sign positive", kTHnSparseD,{multAxis,ptAxiss,invmassAxis});
      histos.add("h3KstarInvMassLikeSignMM", "Invariant mass of Kstar meson Like Sign negative", kTHnSparseD, {multAxis,ptAxiss,invmassAxis});
      histos.add("h3KstarInvMassRotational", "Invariant mass of Kstar meson Rotational", kTHnSparseD, {multAxis,ptAxiss,invmassAxis});
      histos.add("h3KstarInvMassMixed", "Invariant mass of Kstar meson Mixed", kTHnSparseD, {multAxis,ptAxiss,invmassAxis});
    } 
  }





  


  
  double massPi = RecoDecay::getMassPDG(kPiPlus);
  double massKa = RecoDecay::getMassPDG(kKPlus);
  double massPr = 0.938272088f;
  double rapidity;
  double genMass, recMass, resolution;
  double mass{0.};
  double mass_roat{0.};
  double pT{0.};
  double pT_roat{0.};
  array<float, 3> pvec0;
  array<float, 3> pvec1;
  array<float, 3> pvec1_roat;

  /////////////////////////////////////////////////////////////
  template <typename T>
  bool selectionTrack(const T& candidate)
  {
    if (!candidate.isGlobalTrack()) {
      return false;
    }
    return true;
  }

  //////////////////////////////////////////////////////////////



  /////////////////////////////////////////////////////////////
  
  template <typename T>
  bool selectionPID(const T& candidate, int PID)
  {
    if (candidate.hasTOF()) {
      if (PID == 0 && (candidate.tofNSigmaPi() * candidate.tofNSigmaPi() + candidate.tpcNSigmaPi() * candidate.tpcNSigmaPi()) < (2.0 * nsigmaCutCombined * nsigmaCutCombined)) {
        return true;
      } else if (PID == 1 && (candidate.tofNSigmaKa() * candidate.tofNSigmaKa() + candidate.tpcNSigmaKa() * candidate.tpcNSigmaKa()) < (2.0 * nsigmaCutCombined * nsigmaCutCombined)) {
        return true;
      } 
    }
    else {
      if (PID == 0 && std::abs(candidate.tpcNSigmaPi()) < nsigmaCutTPC) {
        return true;
      } else if (PID == 1 && std::abs(candidate.tpcNSigmaKa()) < nsigmaCutTPC) {
        return true;
      } 
    }
    return false;
  }

  /////////////////////////////////////////////////////////////


  ////////////////////////////////////////////////////////////
  
  template <typename T1, typename T2>
  void FillinvMass(const T1& candidate1, const T2& candidate2, float multiplicity, bool unlike, bool mix, bool rotational, bool likesign, float massd1, float massd2)
  {
    pvec0 = array{candidate1.px(), candidate1.py(), candidate1.pz()};
    pvec1 = array{candidate2.px(), candidate2.py(), candidate2.pz()};
    pvec1_roat = array{-candidate2.px(), -candidate2.py(), candidate2.pz()};
    auto arrMom = array{pvec0, pvec1};
    auto arrMom_roat = array{pvec0, pvec1_roat};
    int track1Sign = candidate1.sign();
    int track2Sign = candidate2.sign();
    mass = RecoDecay::m(arrMom, array{massd1, massd2}); //array of momentum and masses
    mass_roat = RecoDecay::m(arrMom_roat, array{massd1, massd2}); //array of momentum and masses
    pT = RecoDecay::pt(array{candidate1.px() + candidate2.px(), candidate1.py() + candidate2.py()});
    pT_roat = RecoDecay::pt(array{candidate1.px() - candidate2.px(), candidate1.py() - candidate2.py()});
    rapidity = RecoDecay::y(array{candidate1.px() + candidate2.px(), candidate1.py() + candidate2.py(), candidate1.pz() + candidate2.pz()}, mass);
    if (std::abs(rapidity) < 0.5) {
      if (track1Sign * track2Sign < 0 && unlike) /// unlike sign
      {
        if (cfgparticletype == 0) {
          histos.fill(HIST("h3KstarInvMassUnlikeSign"), multiplicity, pT, mass);
        } else if (cfgparticletype == 1) {
          histos.fill(HIST("h3PhiInvMassUnlikeSign"), multiplicity, pT, mass);
        } 
      }
      if (track1Sign * track2Sign < 0 && mix) /// mix
      {
        if (cfgparticletype == 0) {
          histos.fill(HIST("h3KstarInvMassMixed"), multiplicity, pT, mass);
        } else if (cfgparticletype == 1) {
          histos.fill(HIST("h3PhiInvMassMixed"), multiplicity, pT, mass);
        } 
      }
      if (track1Sign * track2Sign < 0 && rotational) // rotational
      {
        if (cfgparticletype == 0) {
          histos.fill(HIST("h3KstarInvMassRotational"), multiplicity, pT, mass_roat);
        } else if (cfgparticletype == 1) {
          histos.fill(HIST("h3PhiInvMassRotational"), multiplicity, pT, mass_roat);
        } 
      }
      if (track1Sign * track2Sign > 0 && likesign) /// like sign
      {
        if (track1Sign > 0 && track2Sign > 0) {
          if (cfgparticletype == 0) {
            histos.fill(HIST("h3KstarInvMassLikeSignPP"), multiplicity, pT, mass);
          } else if (cfgparticletype == 1) {
            histos.fill(HIST("h3PhiInvMassLikeSignPP"), multiplicity, pT, mass);
          } 
        } else {
          if (cfgparticletype == 0) {
            histos.fill(HIST("h3KstarInvMassLikeSignMM"), multiplicity, pT, mass);
          } else if (cfgparticletype == 1) {
            histos.fill(HIST("h3PhiInvMassLikeSignMM"), multiplicity, pT, mass);
          } 
        }
      }
    }
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter acceptanceFilter = (nabs(aod::track::eta) < cfgCutEta && nabs(aod::track::pt) > cfgCutPT);
  Filter DCAcutFilter = (nabs(aod::track::dcaXY) < cfgCutDCAxy) && (nabs(aod::track::dcaZ) < cfgCutDCAz);
  using EventCandidates = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::Mults>>;
  using TrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection,
                                                  aod::pidTPCFullPi, aod::pidTOFFullPi,
                                                  aod::pidTPCFullKa, aod::pidTOFFullKa,
                                                  aod::pidTPCFullPr, aod::pidTOFFullPr>>;


  ConfigurableAxis axisVertex{"axisVertex", {20, -10, 10}, "vertex axis for bin"};
  ConfigurableAxis axisMultiplicity{"axisMultiplicity", {VARIABLE_WIDTH, 0.0, 2.750, 5.250, 7.750, 12.750, 17.750, 22.750, 27.750, 32.750, 37.750, 42.750, 47.750, 52.750, 57.750, 62.750, 67.750, 72.750, 77.750, 82.750, 87.750, 92.750, 97.750, 250.1}, "multiplicity axis for histograms"};

  
  using BinningType = ColumnBinningPolicy<aod::collision::PosZ, aod::mult::MultTPC>;
  BinningType binningOnPositions{{axisVertex, axisMultiplicity}, true}; // true is for 'ignore overflows' (true by default)
  SameKindPair<EventCandidates, TrackCandidates, BinningType> pair{binningOnPositions, cfgNoMixedEvents, -1};


  
  
  void processSameEvent(EventCandidates::iterator const& collision, TrackCandidates const& tracks, aod::BCs const&)
  {
    if (!collision.sel8()) {
      return;
    }
    float multiplicity = collision.multTPC();
    histos.fill(HIST("hCentrality"), multiplicity);
    histos.fill(HIST("hNcontributor"), collision.numContrib());
    histos.fill(HIST("hVtxZ"), collision.posZ());
    for (auto track1 : tracks) {
      if (!selectionTrack(track1)) {
        continue;
      }
      histos.fill(HIST("hEta"), track1.eta());
      histos.fill(HIST("hDcaxy"), track1.dcaXY());
      histos.fill(HIST("hDcaz"), track1.dcaZ());
      histos.fill(HIST("hNsigmaKaonTPC"), track1.tpcNSigmaKa());
      histos.fill(HIST("hNsigmaKaonTOF"), track1.tofNSigmaKa());
      histos.fill(HIST("hNsigmaPionTPC"), track1.tpcNSigmaPi());
      histos.fill(HIST("hNsigmaPionTOF"), track1.tofNSigmaPi());
      auto track1ID = track1.globalIndex();
      for (auto track2 : tracks) {
        if (!selectionTrack(track2)) {
          continue;
        }
        auto track2ID = track2.globalIndex();
        if (track2ID <= track1ID) {
          continue;
        }
        bool unlike = true;
        bool mix = false;
        bool rotational = true;
        bool likesign = true;
        if (cfgparticletype == 0) {
          if (selectionPID(track1, 0) && selectionPID(track2, 1)) {
            FillinvMass(track1, track2, multiplicity, unlike, mix, rotational, likesign, massPi, massKa);
          } else if (selectionPID(track1, 1) && selectionPID(track2, 0)) {
            FillinvMass(track1, track2, multiplicity, unlike, mix, rotational, likesign, massKa, massPi);
          }
        } else if (cfgparticletype == 1) {
          if (selectionPID(track1, 1) && selectionPID(track2, 1)) {
            FillinvMass(track1, track2, multiplicity, unlike, mix, rotational, likesign, massKa, massKa);
          }
        }
      }
    }
  }

  
  PROCESS_SWITCH(resonanceqa, processSameEvent, "Process Same event", true);

  void processMixedEvent(EventCandidates const& collisions, TrackCandidates const& tracks)
  {
    for (auto& [c1, tracks1, c2, tracks2] : pair) {
      if (!c1.sel8()) {
        continue;
      }
      if (!c2.sel8()) {
        continue;
      }
      float multiplicity = c1.multTPC();
      for (auto& [t1, t2] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {
        bool unlike = false;
        bool mix = true;
        bool rotational = false;
        bool likesign = false;
        if (!selectionTrack(t1)) {
          continue;
        }
        if (!selectionTrack(t2)) {
          continue;
        }
        if (cfgparticletype == 0) {
          if (selectionPID(t1, 0) && selectionPID(t2, 1)) {
            FillinvMass(t1, t2, multiplicity, unlike, mix, rotational, likesign, massPi, massKa);
          } else if (selectionPID(t1, 1) && selectionPID(t2, 0)) {
            FillinvMass(t1, t2, multiplicity, unlike, mix, rotational, likesign, massKa, massPi);
          }
        } else if (cfgparticletype == 1) {
          if (selectionPID(t1, 1) && selectionPID(t2, 1)) {
            FillinvMass(t1, t2, multiplicity, unlike, mix, rotational, likesign, massKa, massKa);
          }
        } 
      }
    }
  }

  PROCESS_SWITCH(resonanceqa, processMixedEvent, "Process Mixed event", true);

};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<resonanceqa>(cfgc, TaskName{"resonanceqa"})};
}
