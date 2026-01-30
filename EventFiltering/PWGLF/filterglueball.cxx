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

/// \file filterglueball.cxx
/// \brief selection of events with more than 1 K0s for glueball resonance studies
/// \author Sawan <sawan.sawan@cern.ch>

#include "../filterTables.h"

// #include "PWGLF/DataModel/LFStrangenessPIDTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/EventSelection.h" //
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTPC.h" //
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h" //
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h" //
#include "ReconstructionDataFormats/Track.h"

#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "TF1.h"
#include <TDirectory.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <THn.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TPDGCode.h>
#include <TVector2.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <memory>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using std::array;

struct filterglueball {
  // Produce derived tables
  Produces<aod::GlueballFilters> tags;

  HistogramRegistry rEventSelection{"eventSelection", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry hglue{"hglueball", {}, OutputObjHandlingPolicy::AnalysisObject};
  OutputObj<TH1D> hProcessedEvents{TH1D("hProcessedEvents", ";; Number of events", 4, 0.0f, 4.0f)};

  struct : ConfigurableGroup {
    // New configurables
    Configurable<float> MinGlueMass{"MinGlueMass", 1.0, "Min Glue mass"};
    Configurable<float> MaxGlueMass{"MaxGlueMass", 3.0, "Max Glue mass"};
    Configurable<float> MinGluePt{"MinGluePt", 1.0, "Min Glue pT"};

    // Booleans
    Configurable<bool> isApplyEtaCutK0s{"isApplyEtaCutK0s", false, "Apply eta cut on K0s daughters"};
    Configurable<bool> isApplyDCAv0topv{"isApplyDCAv0topv", false, "DCA V0 to PV"};
    Configurable<bool> isapplyCompetingcut{"isapplyCompetingcut", false, "Competing cascade rejection cut"};
    Configurable<bool> isApplySel8{"isApplySel8", true, "Apply sel8 event selection"};
    Configurable<bool> isApplyTimeFrame{"isApplyTimeFrame", false, "Apply Time Frame border selection"};
    Configurable<bool> isApplyITSROF{"isApplyITSROF", false, "Apply ITS ROF border selection"};

    // Configurable parameters for V0 selection
    Configurable<float> cutzvertex{"cutzvertex", 10.0, "z-vertex cut"};
    Configurable<float> confV0PtMin{"confV0PtMin", 0.f, "Minimum transverse momentum of V0"};
    Configurable<float> confV0PtMax{"confV0PtMax", 100.f, "Maximum transverse momentum of V0"};
    Configurable<float> confPiPtMin{"confPiPtMin", 0.1f, "Minimum transverse momentum of pion daughter"};
    Configurable<float> confPiPtMax{"confPiPtMax", 100.f, "Maximum transverse momentum of pion daughter"};
    Configurable<float> confV0DCADaughMax{"confV0DCADaughMax", 1.0f, "DCA b/w V0 daughters"};
    Configurable<float> v0DCAtoPV{"v0DCAtoPV", 0.06, "DCA Pos To PV"};
    Configurable<double> cMaxV0DCA{"cMaxV0DCA", 0.5, "DCA V0 to PV"};
    Configurable<float> confV0CPAMin{"confV0CPAMin", 0.97f, "Minimum CPA of V0"};
    Configurable<float> confV0TranRadV0Min{"confV0TranRadV0Min", 0.5f, "Minimum transverse radius"};
    Configurable<double> cMaxV0LifeTime{"cMaxV0LifeTime", 15, "Maximum V0 life time"};
    Configurable<double> cSigmaMassKs0{"cSigmaMassKs0", 4, "n Sigma cut on Ks0 mass (Mass (Ks) - cSigmaMassKs0*cWidthKs0)"};
    Configurable<double> cWidthKs0{"cWidthKs0", 0.005, "Width of KS0"};
    Configurable<float> confDaughEta{"confDaughEta", 0.8f, "V0 Daugh sel: max eta"};
    Configurable<float> confDaughTPCnclsMin{"confDaughTPCnclsMin", 70.f, "V0 Daugh sel: Min. nCls TPC"};
    Configurable<float> confDaughPIDCutTPC{"confDaughPIDCutTPC", 5, "PID selections for KS0 daughters"};
    Configurable<float> confDaughPIDCutTOF{"confDaughPIDCutTOF", 5, "PID selections for KS0 daughters in TOF"};
    Configurable<float> confKsrapidity{"confKsrapidity", 0.5f, "Rapidity cut on K0s"};
    Configurable<float> cfgETAcut{"cfgETAcut", 0.8f, "Track ETA cut"};
    Configurable<float> cfgPTcut{"cfgPTcut", 0.2f, "Track PT cut"};
    Configurable<float> competingcascrejlambda{"competingcascrejlambda", 0.005, "rejecting competing cascade lambda"};

    ConfigurableAxis ksMassBins{"ksMassBins", {200, 0.45f, 0.55f}, "K0s invariant mass axis"};
    ConfigurableAxis cGlueMassBins{"cGlueMassBins", {200, 0.9f, 3.0f}, "Glueball invariant mass axis"};
    ConfigurableAxis cPtBins{"cPtBins", {200, 0.0f, 20.0f}, "Glueball pT axis"};

  } config;

  void init(o2::framework::InitContext&)
  {
    // Axes
    AxisSpec k0ShortMassAxis = {config.ksMassBins, "#it{M}_{inv} [GeV/#it{c}^{2}]"};
    AxisSpec glueballMassAxis = {config.cGlueMassBins, "#it{M}_{inv} [GeV/#it{c}^{2}]"};
    AxisSpec ptAxis = {config.cPtBins, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec NKsAxis = {15, -0.5, 14.5, "Number of K0s"};
    AxisSpec nSigmaAxis = {100, -10.f, 10.f, "n#sigma TPC"};

    // Event selection
    rEventSelection.add("hVertexZRec", "hVertexZRec", {HistType::kTH1F, {{90, -15.f, 15.f}}});
    rEventSelection.add("hmultiplicity", "multiplicity percentile distribution", {HistType::kTH1F, {{110, 0.0f, 110.0f}}});
    rEventSelection.add("htrackscheck_v0", "htrackscheck_v0", kTH1I, {{15, 0, 15}});
    rEventSelection.add("htrackscheck_v0_daughters", "htrackscheck_v0_daughters", kTH1I, {{15, 0, 15}});

    std::shared_ptr<TH1> hv0label = rEventSelection.get<TH1>(HIST("htrackscheck_v0"));
    hv0label->GetXaxis()->SetBinLabel(1, "All Tracks");
    hv0label->GetXaxis()->SetBinLabel(2, "DCA V0 to PV");
    hv0label->GetXaxis()->SetBinLabel(3, "y K0s");
    hv0label->GetXaxis()->SetBinLabel(4, "V0 pT cut");
    hv0label->GetXaxis()->SetBinLabel(5, "Daughter DCA");
    hv0label->GetXaxis()->SetBinLabel(6, "CosPA");
    hv0label->GetXaxis()->SetBinLabel(7, "Decay Radius");
    hv0label->GetXaxis()->SetBinLabel(8, "Lifetime");
    hv0label->GetXaxis()->SetBinLabel(9, "CompetingCascade");
    hv0label->GetXaxis()->SetBinLabel(10, "Mass Tolerance");

    std::shared_ptr<TH1> hv0DauLabel = rEventSelection.get<TH1>(HIST("htrackscheck_v0_daughters"));
    hv0DauLabel->GetXaxis()->SetBinLabel(1, "AllDau Tracks");
    hv0DauLabel->GetXaxis()->SetBinLabel(2, "TPC clusters");
    hv0DauLabel->GetXaxis()->SetBinLabel(3, "Charge");
    hv0DauLabel->GetXaxis()->SetBinLabel(4, "Eta");
    hv0DauLabel->GetXaxis()->SetBinLabel(5, "PID TPC");
    hv0DauLabel->GetXaxis()->SetBinLabel(6, "PID TOF");
    hv0DauLabel->GetXaxis()->SetBinLabel(7, "Pt cut");

    hProcessedEvents->GetXaxis()->SetBinLabel(1, "All events");
    hProcessedEvents->GetXaxis()->SetBinLabel(2, "Sel8");
    hProcessedEvents->GetXaxis()->SetBinLabel(3, "Events with double K0s without sel.");
    hProcessedEvents->GetXaxis()->SetBinLabel(4, "Selected events");

    hglue.add("hInvMassK0s", "hInvMassK0s", kTHnSparseF, {{k0ShortMassAxis}, {ptAxis}}, true);
    hglue.add("h3glueInvMassDS", "h3glueInvMassDS", kTHnSparseF, {ptAxis, glueballMassAxis}, true);
    hglue.add("NksProduced", "Number of K0s produced", kTH1I, {{NKsAxis}}, true);
    hglue.add("hNSigmaPosPionK0s_before", "hNSigmaPosPionK0s_before", {HistType::kTH2F, {{ptAxis}, {nSigmaAxis}}});
    hglue.add("hNSigmaNegPionK0s_before", "hNSigmaNegPionK0s_before", {HistType::kTH2F, {{ptAxis}, {nSigmaAxis}}});
  }

  template <typename Collision, typename V0>
  bool selectionV0(Collision const& collision, V0 const& candidate, float /*multiplicity*/)
  {
    const float pT = candidate.pt();
    const float tranRad = candidate.v0radius();
    const float dcaDaughv0 = candidate.dcaV0daughters();
    const float cpav0 = candidate.v0cosPA();

    float ctauK0s = candidate.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassK0Short;
    float lowmasscutks0 = o2::constants::physics::MassKPlus - config.cWidthKs0 * config.cSigmaMassKs0;
    float highmasscutks0 = o2::constants::physics::MassKPlus + config.cWidthKs0 * config.cSigmaMassKs0;

    rEventSelection.fill(HIST("htrackscheck_v0"), 0.5);

    if (config.isApplyDCAv0topv && std::fabs(candidate.dcav0topv()) > config.cMaxV0DCA) {
      return false;
    }
    rEventSelection.fill(HIST("htrackscheck_v0"), 1.5);

    if (std::abs(candidate.rapidity(0)) >= config.confKsrapidity) {
      return false;
    }
    rEventSelection.fill(HIST("htrackscheck_v0"), 2.5);

    if (pT < config.confV0PtMin || pT > config.confV0PtMax) {
      return false;
    }
    rEventSelection.fill(HIST("htrackscheck_v0"), 3.5);

    if (dcaDaughv0 > config.confV0DCADaughMax) {
      return false;
    }
    rEventSelection.fill(HIST("htrackscheck_v0"), 4.5);

    if (cpav0 < config.confV0CPAMin) {
      return false;
    }
    rEventSelection.fill(HIST("htrackscheck_v0"), 5.5);

    if (tranRad < config.confV0TranRadV0Min) {
      return false;
    }
    rEventSelection.fill(HIST("htrackscheck_v0"), 6.5);

    if (std::fabs(ctauK0s) > config.cMaxV0LifeTime) {
      return false;
    }
    rEventSelection.fill(HIST("htrackscheck_v0"), 7.5);

    if (config.isapplyCompetingcut && (std::abs(candidate.mLambda() - o2::constants::physics::MassLambda0) <= config.competingcascrejlambda || std::abs(candidate.mAntiLambda() - o2::constants::physics::MassLambda0) <= config.competingcascrejlambda)) {
      return false;
    }
    rEventSelection.fill(HIST("htrackscheck_v0"), 8.5);

    if (candidate.mK0Short() < lowmasscutks0 || candidate.mK0Short() > highmasscutks0) {
      return false;
    }
    rEventSelection.fill(HIST("htrackscheck_v0"), 9.5);

    return true;
  }

  template <typename T>
  bool isSelectedV0Daughter(T const& track, float charge, double nsigmaV0DaughterTPC)
  {
    const auto eta = track.eta();
    const auto tpcNClsF = track.tpcNClsFound();
    const auto sign = track.sign();

    rEventSelection.fill(HIST("htrackscheck_v0_daughters"), 0.5);

    if (tpcNClsF < config.confDaughTPCnclsMin) {
      return false;
    }
    rEventSelection.fill(HIST("htrackscheck_v0_daughters"), 1.5);

    if (charge < 0 && sign > 0) {
      return false;
    }

    if (charge > 0 && sign < 0) {
      return false;
    }
    rEventSelection.fill(HIST("htrackscheck_v0_daughters"), 2.5);

    if (std::abs(eta) > config.confDaughEta) {
      return false;
    }
    rEventSelection.fill(HIST("htrackscheck_v0_daughters"), 3.5);

    if (std::abs(nsigmaV0DaughterTPC) > config.confDaughPIDCutTPC) {
      return false;
    }
    rEventSelection.fill(HIST("htrackscheck_v0_daughters"), 4.5);

    // if (std::abs(v0candidate.tofNSigmaK0PiPlus()) > config.confDaughPIDCutTOF && v0candidate.positiveHasTOF()) {
    //   return false;
    // }

    // if (std::abs(v0candidate.tofNSigmaK0PiMinus()) > config.confDaughPIDCutTOF && v0candidate.negativeHasTOF()) {
    //   return false;
    // }
    rEventSelection.fill(HIST("htrackscheck_v0_daughters"), 5.5);

    if (track.pt() < config.confPiPtMin || track.pt() > config.confPiPtMax) {
      return false;
    }
    rEventSelection.fill(HIST("htrackscheck_v0_daughters"), 6.5);

    (charge == 1) ? hglue.fill(HIST("hNSigmaPosPionK0s_before"), track.pt(), track.tpcNSigmaPi()) : hglue.fill(HIST("hNSigmaNegPionK0s_before"), track.pt(), track.tpcNSigmaPi());

    return true;
  }

  Filter posZFilter = (nabs(o2::aod::collision::posZ) < config.cutzvertex);
  Filter acceptenceFilter = (nabs(aod::track::eta) < config.cfgETAcut && nabs(aod::track::pt) > config.cfgPTcut);

  // Filters on V0s
  Filter preFilterV0 = (nabs(aod::v0data::dcapostopv) > config.v0DCAtoPV && nabs(aod::v0data::dcanegtopv) > config.v0DCAtoPV);

  // Defining the type of the daughter tracks
  using EventCandidates = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms>>;
  using TrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra, aod::TracksDCA, aod::pidTPCFullPi>>;
  // using V0TrackCandidate = soa::Filtered<soa::Join<aod::V0Datas>>;
  using V0TrackCandidate = aod::V0Datas;

  ROOT::Math::PxPyPzMVector K0sVectorDummy, K0sVectorDummy2, GlueballVector;
  void processGlue(EventCandidates::iterator const& collision, TrackCandidates const&, V0TrackCandidate const& V0s)
  {
    hProcessedEvents->Fill(0.5);
    bool keepEventDoubleK0s = false;
    int nK0sInEvent = 0;

    std::vector<ROOT::Math::PtEtaPhiMVector> K0sParticle;
    std::vector<int64_t> PosdauIndex = {};
    std::vector<int64_t> NegdauIndex = {};

    if (config.isApplySel8) {
      if (!collision.sel8()) {
        return;
      }
    } else {
      if (!collision.selection_bit(aod::evsel::kIsTriggerTVX)) {
        return;
      }
    }

    // Independent conditions
    if (!config.isApplySel8 && config.isApplyTimeFrame && !collision.selection_bit(aod::evsel::kNoTimeFrameBorder)) {
      return;
    }

    if (!config.isApplySel8 && config.isApplyITSROF && !collision.selection_bit(aod::evsel::kNoITSROFrameBorder)) {
      return;
    }
    hProcessedEvents->Fill(1.5);

    rEventSelection.fill(HIST("hVertexZRec"), collision.posZ());
    rEventSelection.fill(HIST("hmultiplicity"), collision.centFT0M());

    for (const auto& V0Track : V0s) {

      if (V0Track.size() == 0) {
        continue;
      }

      if (!selectionV0(collision, V0Track, collision.centFT0M())) {
        continue;
      }

      auto postrack1 = V0Track.template posTrack_as<TrackCandidates>();
      auto negtrack1 = V0Track.template negTrack_as<TrackCandidates>();

      double nTPCSigmaPos1{postrack1.tpcNSigmaPi()};
      double nTPCSigmaNeg1{negtrack1.tpcNSigmaPi()};

      if (!isSelectedV0Daughter(negtrack1, -1, nTPCSigmaNeg1)) {
        continue;
      }

      if (!isSelectedV0Daughter(postrack1, 1, nTPCSigmaPos1)) {
        continue;
      }

      if (config.isApplyEtaCutK0s && (V0Track.eta() < config.confDaughEta || V0Track.eta() < config.confDaughEta)) {
        continue;
      }
      nK0sInEvent++;

      K0sParticle.push_back(ROOT::Math::PtEtaPhiMVector(V0Track.pt(), V0Track.eta(), V0Track.phi(), V0Track.mK0Short()));
      hglue.fill(HIST("hInvMassK0s"), V0Track.mK0Short(), V0Track.pt());
      PosdauIndex.push_back(postrack1.globalIndex());
      NegdauIndex.push_back(negtrack1.globalIndex());
    }
    if (nK0sInEvent > 1) {
      hglue.fill(HIST("NksProduced"), nK0sInEvent);
      hProcessedEvents->Fill(2.5);

      for (auto if1 = K0sParticle.begin(); if1 != K0sParticle.end(); ++if1) {
        auto i5 = std::distance(K0sParticle.begin(), if1);
        K0sVectorDummy = K0sParticle.at(i5);

        for (auto if2 = if1 + 1; if2 != K0sParticle.end(); ++if2) {
          auto i6 = std::distance(K0sParticle.begin(), if2);
          K0sVectorDummy2 = K0sParticle.at(i6);
          GlueballVector = K0sVectorDummy + K0sVectorDummy2;

          if (PosdauIndex.at(i5) == PosdauIndex.at(i6) || PosdauIndex.at(i5) == NegdauIndex.at(i6) || NegdauIndex.at(i5) == PosdauIndex.at(i6) || NegdauIndex.at(i5) == NegdauIndex.at(i6)) {
            continue;
          }

          if (GlueballVector.M() > config.MinGlueMass && GlueballVector.M() < config.MaxGlueMass && GlueballVector.Pt() > config.MinGluePt) {
            hglue.fill(HIST("h3glueInvMassDS"), GlueballVector.Pt(), GlueballVector.M());
            keepEventDoubleK0s = true;
          }
        }
      }
    }
    if (keepEventDoubleK0s) {
      hProcessedEvents->Fill(3.5);
    }
    tags(keepEventDoubleK0s);
  }
  PROCESS_SWITCH(filterglueball, processGlue, "same event process", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<filterglueball>(cfgc)};
}
