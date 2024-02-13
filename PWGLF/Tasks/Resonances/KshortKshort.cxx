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
/// \brief glueball resonance
/// \author Sawan (sawan.sawan@cern.ch)
#include <TDatabasePDG.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <THn.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TPDGCode.h>

#include <array>
#include <cmath>
#include <cstdlib>

#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/StepTHn.h"
#include "ReconstructionDataFormats/Track.h"

#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h" //
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponse.h" //
#include "Common/DataModel/TrackSelectionTables.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisTask.h"              //
#include "Framework/runDataProcessing.h"         //
#include "PWGLF/DataModel/LFStrangenessTables.h" //
// #include "CommonConstants/PhysicsConstants.h"
// #include "Framework/HistogramRegistry.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
// using namespace o2::constants::physics;
using std::array;

struct strangeness_tutorial {
  SliceCache cache;

  HistogramRegistry rEventSelection{"eventSelection", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rKzeroShort{"kzeroShort", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry hglue{"hglueball", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  // Configurable for event selection
  Configurable<float> cutzvertex{"cutzvertex", 10.0f, "Accepted z-vertex range (cm)"};
  Configurable<float> cfgETAcut{"cfgETAcut", 0.8f, "Track ETA cut"};

  // Configurable parameters for V0 selection
  Configurable<float> v0setting_dcav0dau{"v0setting_dcav0dau", 1, "DCA V0 Daughters"};
  Configurable<float> v0setting_dcapostopv{"v0setting_dcapostopv", 0.06, "DCA Pos To PV"};
  Configurable<float> v0setting_dcanegtopv{"v0setting_dcanegtopv", 0.06, "DCA Neg To PV"};
  Configurable<float> ConfV0PtMin{"ConfV0PtMin", 0.f, "Minimum transverse momentum of V0"};
  Configurable<float> ConfV0DCADaughMax{"ConfV0DCADaughMax", 1.0f, "Maximum DCA between the V0 daughters"};
  Configurable<float> ConfV0CPAMin{"ConfV0CPAMin", 0.97f, "Minimum CPA of V0"};
  Configurable<float> ConfV0TranRadV0Min{"ConfV0TranRadV0Min", 0.5f, "Minimum transverse radius"};
  Configurable<float> ConfV0TranRadV0Max{"ConfV0TranRadV0Max", 200.f, "Maximum transverse radius"};
  Configurable<double> cMaxV0LifeTime{"cMaxV0LifeTime", 15, "Maximum V0 life time"};
  Configurable<double> cMaxV0DCA{"cMaxV0DCA", 0.3, "DCA V0 to PV"};
  Configurable<double> cSigmaMassKs0{"cSigmaMassKs0", 4, "n Sigma cut on KS0 mass"};
  Configurable<double> cWidthKs0{"cWidthKs0", 0.005, "Width of KS0"};
  Configurable<float> ConfDaughEta{"ConfDaughEta", 0.8f, "V0 Daugh sel: max eta"};
  Configurable<float> ConfDaughTPCnclsMin{"ConfDaughTPCnclsMin", 70.f, "V0 Daugh sel: Min. nCls TPC"};
  Configurable<float> ConfDaughDCAMin{"ConfDaughDCAMin", 0.06f, "V0 Daugh sel:  Max. DCA Daugh to PV (cm)"};
  Configurable<float> ConfDaughPIDCuts{"ConfDaughPIDCuts", 4, "PID selections for KS0 daughters"};

  // Configurable parameters for PID selection
  //   Configurable<float> NSigmaTPCPion{"NSigmaTPCPion", 4, "NSigmaTPCPion"};

  // Configurable for track selection
  Configurable<float> cfgTrackDCAxy{"cfgTrackDCAxy", 2.0f, "Track DCAxy"};
  Configurable<float> cfgTrackDCAz{"cfgTrackDCAz", 2.0f, "Track DCAz"};
  Configurable<float> cfgPTcut{"cfgPTcut", 0.2f, "Track PT cut"};
  Configurable<float> cfgnSigmaTPCcut{"cfgnSigmaTPCcut", 3.0f, "Track nSigma TPC cut"};
  Configurable<float> cfgnSigmaTOFcut{"cfgnSigmaTOFcut", 3.0f, "Track nSigma TOF cut"};
  Configurable<float> nsigmaCutCombined{"nsigmaCutCombined", 3.0f, "Track nSigma Combined cut"};
  Configurable<int> cfgNmixedEvents{"cfgNmixedEvents", 5, "Number of mixed events"};
  Configurable<bool> cfgMultFT0{"cfgMultFT0", false, "Use FT0 multiplicity"};
  Configurable<bool> cfgMultFOTM{"cfgMultFOTM", false, "Use FOTM multiplicity"};
  Configurable<bool> cfgMultFT0C{"cfgMultFT0C", true, "Use FT0C multiplicity"};
  ConfigurableAxis binsCent{"binsCent", {VARIABLE_WIDTH, 0., 5., 10., 30., 50., 70., 100., 110., 150.}, "Binning of the centrality axis"};

  void init(InitContext const&)
  {
    // Axes
    AxisSpec K0ShortMassAxis = {200, 0.45f, 0.55f, "#it{M}_{inv} [GeV/#it{c}^{2}]"};
    AxisSpec glueballMassAxis = {150, 0.9f, 2.4f, "#it{M}_{inv} [GeV/#it{c}^{2}]"};
    AxisSpec vertexZAxis = {100, -15.f, 15.f, "vrtx_{Z} [cm]"}; // for histogram
    AxisSpec ptAxis = {150, 0.0f, 15.0f, "#it{p}_{T} (GeV/#it{c})"};
    // AxisSpec multiplicityAxis = {110, 0.0f, 150.0f, "Multiplicity Axis"};
    AxisSpec multiplicityAxis = {binsCent, "Multiplicity Axis"};

    // Histograms
    // Event selection
    rEventSelection.add("hVertexZRec", "hVertexZRec", {HistType::kTH1F, {vertexZAxis}});
    rEventSelection.add("hmultiplicity", "hmultiplicity", {HistType::kTH1F, {{150, 0.0f, 150.0f}}});

    // Invariant Mass
    rKzeroShort.add("hMassK0ShortSelected", "hMassK0ShortSelected", {HistType::kTH1F, {K0ShortMassAxis}});
    hglue.add("h1glueInvMassDS", "h1glueInvMassDS", kTH1F, {glueballMassAxis});
    hglue.add("h1glueInvMassME", "h1glueInvMassME", kTH1F, {glueballMassAxis});
    hglue.add("h3glueInvMassDS", "h3glueInvMassDS", kTH3F, {multiplicityAxis, ptAxis, glueballMassAxis});
    hglue.add("h3glueInvMassME", "h3glueInvMassME", kTH3F, {multiplicityAxis, ptAxis, glueballMassAxis});

    // K0s topological/PID cuts
    rKzeroShort.add("hDCAV0Daughters", "hDCAV0Daughters", {HistType::kTH1F, {{55, 0.0f, 2.2f}}});
    rKzeroShort.add("hV0CosPA", "hV0CosPA", {HistType::kTH1F, {{100, 0.95f, 1.f}}});
    rKzeroShort.add("hNSigmaPosPionFromK0s", "hNSigmaPosPionFromK0s", {HistType::kTH2F, {{ptAxis}, {100, -5.f, 5.f}}});
    rKzeroShort.add("hNSigmaNegPionFromK0s", "hNSigmaNegPionFromK0s", {HistType::kTH2F, {{ptAxis}, {100, -5.f, 5.f}}});
  }

  template <typename Collision, typename V0>
  bool SelectionV0(Collision const& collision, V0 const& candidate,
                   float multiplicity)
  {
    if (fabs(candidate.dcav0topv()) > cMaxV0DCA) {
      return false;
    }

    if (TMath::Abs(candidate.yK0Short()) > 0.5) {
      return false;
    }

    const float qtarm = candidate.qtarm();
    const float alph = candidate.alpha();
    float arm = qtarm / alph;
    const float pT = candidate.pt();
    const float tranRad = candidate.v0radius();
    const float dcaDaughv0 = candidate.dcaV0daughters();
    const float cpav0 = candidate.v0cosPA();

    float CtauK0s = candidate.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * TDatabasePDG::Instance()->GetParticle(kK0Short)->Mass(); // FIXME: Get from the common header
    float lowmasscutks0 = 0.497 - cWidthKs0 * cSigmaMassKs0;
    float highmasscutks0 = 0.497 + cWidthKs0 * cSigmaMassKs0;
    // float decayLength = candidate.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * RecoDecay::sqrtSumOfSquares(candidate.px(), candidate.py(), candidate.pz());

    if (pT < ConfV0PtMin) {
      return false;
    }
    if (dcaDaughv0 > ConfV0DCADaughMax) {
      return false;
    }
    if (cpav0 < ConfV0CPAMin) {
      return false;
    }
    if (tranRad < ConfV0TranRadV0Min) {
      return false;
    }
    if (tranRad > ConfV0TranRadV0Max) {
      return false;
    }
    if (fabs(CtauK0s) > cMaxV0LifeTime ||
        candidate.mK0Short() < lowmasscutks0 ||
        candidate.mK0Short() > highmasscutks0) {
      return false;
    }
    if (arm < 0.2) {
      return false;
    }

    // if (QAv0) {
    //   histos.fill(HIST("hLT"), CtauK0s);
    //   histos.fill(HIST("hMassvsptvsmult"), candidate.mK0Short(), candidate.pt(),
    //               multiplicity);
    //   histos.fill(HIST("hDCAV0Daughters"), candidate.dcaV0daughters());
    //   histos.fill(HIST("hV0CosPA"), candidate.v0cosPA());
    // }
    return true;
  }

  template <typename T>
  bool isSelectedV0Daughter(T const& track, float charge,
                            double nsigmaV0Daughter)
  {
    const auto eta = track.eta();
    const auto tpcNClsF = track.tpcNClsFound();
    const auto dcaXY = track.dcaXY(); // for this we need TrackDCA table
    const auto sign = track.sign();

    if (!track.hasTPC())
      return false;
    if (track.tpcNClsCrossedRows() < 70)
      return false;
    if (track.tpcCrossedRowsOverFindableCls() < 0.8)
      return false;

    if (charge < 0 && sign > 0) {
      return false;
    }
    if (charge > 0 && sign < 0) {
      return false;
    }
    if (std::abs(eta) > ConfDaughEta) {
      return false;
    }
    if (tpcNClsF < ConfDaughTPCnclsMin) {
      return false;
    }
    if (std::abs(dcaXY) < ConfDaughDCAMin) {
      return false;
    }
    // v0 PID selection
    if (std::abs(nsigmaV0Daughter) > ConfDaughPIDCuts) {
      return false;
    }

    return true;
  }

  template <typename T>
  bool PIDselection(T const& candidate)
  {
    if (candidate.hasTOF() && (candidate.tofNSigmaPi() * candidate.tofNSigmaPi() + candidate.tpcNSigmaPi() * candidate.tpcNSigmaPi()) < (nsigmaCutCombined * nsigmaCutCombined)) {
      return true;
    }
    if (!candidate.hasTOF() && std::abs(candidate.tpcNSigmaPi()) < cfgnSigmaTPCcut) {
      return true;
    }
    return false;
  }

  // Defining filters for events (event selection)
  // Filter eventFilter = (o2::aod::evsel::sel8 == true);
  Filter posZFilter = (nabs(o2::aod::collision::posZ) < cutzvertex);
  Filter AcceptenceFilter = (nabs(aod::track::eta) < cfgETAcut && nabs(aod::track::pt) > cfgPTcut);
  Filter DCAcutFilter = (nabs(aod::track::dcaXY) < cfgTrackDCAxy) && (nabs(aod::track::dcaZ) < cfgTrackDCAz);

  // Filters on V0s
  // Filter preFilterV0 = (nabs(aod::v0data::dcapostopv) > v0setting_dcapostopv &&
  //                       nabs(aod::v0data::dcanegtopv) > v0setting_dcanegtopv &&
  //                       aod::v0data::dcaV0daughters < v0setting_dcav0dau);

  // Defining the type of the daughter tracks
  using DaughterTracks = soa::Join<aod::TracksIU, aod::TracksExtra, aod::pidTPCFullPi, aod::pidTOFFullPi>;
  using EventCandidates = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::MultZeqs, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs>>;
  using TrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra, aod::TracksDCA, aod::pidTPCFullPi, aod::pidTOFFullPi>>;
  using V0TrackCandidate = aod::V0Datas;

  //   void processSE(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& collision,
  //                  soa::Filtered<aod::V0Datas> const& V0s,
  //                  DaughterTracks const&)

  void processSE(EventCandidates::iterator const& collision, TrackCandidates const& tracks, aod::V0Datas const& V0s)
  {
    const double massK0s = TDatabasePDG::Instance()->GetParticle(kK0Short)->Mass();
    if (!collision.sel8()) {
      return;
    }
    float multiplicity = 0.0f;
    if (cfgMultFT0)
      multiplicity = collision.multZeqFT0A() + collision.multZeqFT0C();
    else if (cfgMultFOTM)
      multiplicity = collision.centFT0M();
    else if (cfgMultFT0C)
      multiplicity = collision.centFT0C();

    rEventSelection.fill(HIST("hVertexZRec"), collision.posZ());
    rEventSelection.fill(HIST("hmultiplicity"), multiplicity);

    for (auto& [v1, v2] : combinations(CombinationsStrictlyUpperIndexPolicy(V0s, V0s))) {
      if (v1.size() == 0 || v2.size() == 0) {
        continue;
      }
      auto postrack1 = v1.template posTrack_as<TrackCandidates>();
      auto negtrack1 = v1.template negTrack_as<TrackCandidates>();
      auto postrack2 = v2.template posTrack_as<TrackCandidates>();
      auto negtrack2 = v2.template negTrack_as<TrackCandidates>();
      if (postrack1.globalIndex() == postrack2.globalIndex()) {
        continue;
      }
      if (negtrack1.globalIndex() == negtrack2.globalIndex()) {
        continue;
      }
      double nTPCSigmaPos1[1]{postrack1.tpcNSigmaPi()};
      double nTPCSigmaNeg1[1]{negtrack1.tpcNSigmaPi()};
      double nTPCSigmaPos2[1]{postrack2.tpcNSigmaPi()};
      double nTPCSigmaNeg2[1]{negtrack2.tpcNSigmaPi()};

      if (!isSelectedV0Daughter(postrack1, 1, nTPCSigmaPos1[0])) {
        continue;
      }
      if (!isSelectedV0Daughter(postrack2, 1, nTPCSigmaPos2[0])) {
        continue;
      }
      if (!isSelectedV0Daughter(negtrack1, -1, nTPCSigmaNeg1[0])) {
        continue;
      }
      if (!isSelectedV0Daughter(negtrack2, -1, nTPCSigmaNeg2[0])) {
        continue;
      }

      if (!SelectionV0(collision, v1, multiplicity)) {
        continue;
      }
      if (!SelectionV0(collision, v2, multiplicity)) {
        continue;
      }

      // if (!PIDselection(postrack1)) {
      //   continue;
      // }
      // if (!PIDselection(postrack2)) {
      //   continue;
      // }
      // if (!PIDselection(negtrack1)) {
      //   continue;
      // }
      // if (!PIDselection(negtrack2)) {
      //   continue;
      // }

      TLorentzVector lv1, lv2, lv3;
      lv1.SetPtEtaPhiM(v1.pt(), v1.eta(), v1.phi(), massK0s);
      lv2.SetPtEtaPhiM(v2.pt(), v2.eta(), v2.phi(), massK0s);
      lv3 = lv1 + lv2;

      if (TMath::Abs(lv3.Rapidity() < 0.5)) {
        hglue.fill(HIST("h3glueInvMassDS"), multiplicity, lv3.Pt(), lv3.M());
        hglue.fill(HIST("h1glueInvMassDS"), lv3.M());
        rKzeroShort.fill(HIST("hMassK0ShortSelected"), v1.mK0Short());
        rKzeroShort.fill(HIST("hMassK0ShortSelected"), v2.mK0Short());
        rKzeroShort.fill(HIST("hDCAV0Daughters"), v1.dcaV0daughters());
        rKzeroShort.fill(HIST("hDCAV0Daughters"), v2.dcaV0daughters());
        rKzeroShort.fill(HIST("hV0CosPA"), v1.v0cosPA());
        rKzeroShort.fill(HIST("hV0CosPA"), v2.v0cosPA());
      }

      // Filling the PID of the V0 daughters in the region of the K0 peak.
      // tpcInnerParam is the momentum at the inner wall of TPC. So momentum of tpc vs nsigma of tpc is plotted.
      if ((0.45 < v1.mK0Short() || 0.45 < v2.mK0Short()) && (v1.mK0Short() < 0.55 || v2.mK0Short() < 0.55)) {
        rKzeroShort.fill(HIST("hNSigmaPosPionFromK0s"), postrack1.tpcInnerParam(), postrack1.tpcNSigmaPi());
        rKzeroShort.fill(HIST("hNSigmaPosPionFromK0s"), postrack2.tpcInnerParam(), postrack2.tpcNSigmaPi());
        rKzeroShort.fill(HIST("hNSigmaNegPionFromK0s"), negtrack1.tpcInnerParam(), negtrack1.tpcNSigmaPi());
        rKzeroShort.fill(HIST("hNSigmaNegPionFromK0s"), negtrack2.tpcInnerParam(), negtrack2.tpcNSigmaPi());
      }
    }
  }

  PROCESS_SWITCH(strangeness_tutorial, processSE, "same event process", true);

  // use any one of 3 alias depending on the dataset. If pp then FT0M and if pbpb then FTOC
  using BinningTypeTPCMultiplicity = ColumnBinningPolicy<aod::collision::PosZ, aod::mult::MultTPC>;
  using BinningTypeCentralityM = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0M>;
  using BinningTypeVertexContributor = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0C>;
  ConfigurableAxis mevz = {"mevz", {10, -10., 10.}, "mixed event vertex z binning"};
  ConfigurableAxis memult = {"memult", {2, 0., 110.}, "mixed event multiplicity binning"};
  BinningTypeVertexContributor binningOnPositions1{{mevz, memult}, true};
  BinningTypeCentralityM binningOnPositions2{{mevz, memult}, true};

  SameKindPair<EventCandidates, V0TrackCandidate, BinningTypeVertexContributor> pair1{binningOnPositions1, cfgNmixedEvents, -1, &cache}; // for PbPb
  SameKindPair<EventCandidates, V0TrackCandidate, BinningTypeCentralityM> pair2{binningOnPositions2, cfgNmixedEvents, -1, &cache};       // for pp

  void processME(EventCandidates const& collisions, TrackCandidates const& tracks, V0TrackCandidate const& v0s)
  {
    const double massK0s = TDatabasePDG::Instance()->GetParticle(kK0Short)->Mass();

    for (auto& [c1, tracks1, c2, tracks2] : pair1) // two different centrality c1 and c2 and tracks corresponding to them
    {

      if (!c1.sel8()) {
        continue;
      }
      if (!c2.sel8()) {
        continue;
      }

      float multiplicity = 0.0f;
      if (cfgMultFT0)
        multiplicity = c1.multZeqFT0A() + c1.multZeqFT0C();
      else if (cfgMultFOTM)
        multiplicity = c1.centFT0M();
      else if (cfgMultFT0C)
        multiplicity = c1.centFT0C();

      for (auto& [t1, t2] : o2::soa::combinations(o2::soa::CombinationsStrictlyUpperIndexPolicy(tracks1, tracks2))) {
        if (t1.size() == 0 || t2.size() == 0) {
          continue;
        }

        if (!SelectionV0(c1, t1, multiplicity))
          continue;
        if (!SelectionV0(c2, t2, multiplicity))
          continue;

        auto postrack1 = t1.template posTrack_as<TrackCandidates>();
        auto negtrack1 = t1.template negTrack_as<TrackCandidates>();
        auto postrack2 = t2.template posTrack_as<TrackCandidates>();
        auto negtrack2 = t2.template negTrack_as<TrackCandidates>();
        if (postrack1.globalIndex() == postrack2.globalIndex()) {
          continue;
        }
        if (negtrack1.globalIndex() == negtrack2.globalIndex()) {
          continue;
        }
        double nTPCSigmaPos1[1]{postrack1.tpcNSigmaPi()};
        double nTPCSigmaNeg1[1]{negtrack1.tpcNSigmaPi()};
        double nTPCSigmaPos2[1]{postrack2.tpcNSigmaPi()};
        double nTPCSigmaNeg2[1]{negtrack2.tpcNSigmaPi()};

        if (!isSelectedV0Daughter(postrack1, 1, nTPCSigmaPos1[0])) {
          continue;
        }
        if (!isSelectedV0Daughter(postrack2, 1, nTPCSigmaPos2[0])) {
          continue;
        }
        if (!isSelectedV0Daughter(negtrack1, -1, nTPCSigmaNeg1[0])) {
          continue;
        }
        if (!isSelectedV0Daughter(negtrack2, -1, nTPCSigmaNeg2[0])) {
          continue;
        }

        // if (!PIDselection(postrack1)) {
        //   continue;
        // }
        // if (!PIDselection(postrack2)) {
        //   continue;
        // }
        // if (!PIDselection(negtrack1)) {
        //   continue;
        // }
        // if (!PIDselection(negtrack2)) {
        //   continue;
        // }

        TLorentzVector lv1, lv2, lv3;
        lv1.SetPtEtaPhiM(t1.pt(), t1.eta(), t1.phi(), massK0s);
        lv2.SetPtEtaPhiM(t2.pt(), t2.eta(), t2.phi(), massK0s);
        lv3 = lv1 + lv2;
        if (TMath::Abs(lv3.Rapidity() < 0.5)) {
          hglue.fill(HIST("h1glueInvMassME"), lv3.M());
          hglue.fill(HIST("h3glueInvMassME"), multiplicity, lv3.Pt(), lv3.M());
        }
      }
    }
  }
  PROCESS_SWITCH(strangeness_tutorial, processME, "mixed event process", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<strangeness_tutorial>(cfgc)};
}
