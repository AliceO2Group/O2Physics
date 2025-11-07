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
/// \brief sigma resonance
/// \author sarjeeta gami (sarjeeta.gami@cern.ch)
#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "CCDB/CcdbApi.h"
#include "CommonConstants/PhysicsConstants.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/StepTHn.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

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

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;
using std::array;

struct sigma {
  // Histograms are defined with HistogramRegistry
  HistogramRegistry rEventSelection{"eventSelection", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  // Configurable for histograms
  Configurable<int> nBins{"nBins", 100, "N bins in all histos"};

  // Confugrable for QA histograms
  Configurable<bool> QAbefore{"QAbefore", true, "QAbefore"};
  Configurable<bool> QAafter{"QAafter", true, "QAafter"};
  Configurable<bool> QAv0{"QAv0", true, "QAv0"};

  // Configurable for event selection
  Configurable<float> cutzvertex{"cutzvertex", 10.0f, "Accepted z-vertex range (cm)"};

  // Configurable parameters for V0 selection
  Configurable<float> ConfV0DCADaughMax{"ConfV0DCADaughMax", 1.6f, "Maximum DCA between the V0 daughters"};
  Configurable<float> ConfV0CPAMin{"ConfV0CPAMin", 0.99f, "Minimum CPA of V0"};
  Configurable<float> ConfV0TranRadV0Min{"ConfV0TranRadV0Min", 1.4f, "Minimum transverse radius"};
  Configurable<float> ConfV0TranRadV0Max{"ConfV0TranRadV0Max", 100.f, "Maximum transverse radius"};
  Configurable<double> cMaxV0LifeTime{"cMaxV0LifeTime", 30, "Maximum V0 life time"};
  Configurable<double> cMaxV0DCA{"cMaxV0DCA", 0.3, "DCA V0 to PV"};
  Configurable<float> ConfDaughEta{"ConfDaughEta", 0.8f, "V0 Daugh sel: max eta"};
  Configurable<float> ConfDaughTPCnclsMin{"ConfDaughTPCnclsMin", 70.f, "V0 Daugh sel: Min. nCls TPC"};
  Configurable<float> ConfDaughDCAMin{"ConfDaughDCAMin", 0.05f, "V0 Daugh sel:  Max. DCA Daugh to PV (cm)"};
  Configurable<float> ConfDaughPIDCuts{"ConfDaughPIDCuts", 4, "PID selections for lambda daughters"};
  // Configurables for track selections
  Configurable<float> cfgCutPT{"cfgCutPT", 0.2f, "PT cut on daughter track"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8f, "Eta cut on daughter track"};
  Configurable<float> cfgCutDCAxy{"cfgCutDCAxy", 2.0f, "DCAxy range for tracks"};
  Configurable<float> cfgCutDCAz{"cfgCutDCAz", 2.0f, "DCAz range for tracks"};
  Configurable<float> nsigmaCutTPC{"nsigmacutTPC", 3.0, "Value of the TPC Nsigma cut"};
  Configurable<float> nsigmaCutCombined{"nsigmaCutCombined", 3.0, "Value of the Combined Nsigma cut"};
  Configurable<int> cfgNoMixedEvents{"cfgNoMixedEvents", 5, "Number of mixed events per event"};
  Configurable<bool> cfgMultFT0{"cfgMultFT0", false, "cfgMultFT0"};
  Configurable<bool> cfgCentFT0C{"cfgCentFT0C", true, "cfgCentFT0C"};
  Configurable<bool> iscustomDCAcut{"iscustomDCAcut", false, "iscustomDCAcut"};
  Configurable<bool> ismanualDCAcut{"ismanualDCAcut", true, "ismanualDCAcut"};
  Configurable<int> cfgITScluster{"cfgITScluster", 0, "Number of ITS cluster"};

  void init(InitContext const&)
  {
    // Axes
    AxisSpec lambdaMassAxis = {200, 0.8f, 1.5f, "#it{M}_{inv} [GeV/#it{c}^{2}]"};
    AxisSpec invMassAxisReso = {1600, 0.9f, 2.5f, "Invariant Mass (GeV/#it{c}^2)"};
    AxisSpec vertexZAxis = {nBins, -10., 10., "vrtx_{Z} [cm]"};
    AxisSpec ptAxis = {200, 0.0f, 20.0f, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec multAxis = {100, 0.0f, 100.0f, "Multiplicity"};

    // Histograms
    // Event selection
    rEventSelection.add("hVertexZRec", "hVertexZRec", {HistType::kTH1F, {vertexZAxis}});
    rEventSelection.add("hmult", "Centrality distribution", kTH1F, {{200, 0.0f, 200.0f}});

    // for primary tracks
    if (QAbefore && QAafter) {
      histos.add("hNsigmaPionTPC_before", "NsigmaPion TPC distribution before", kTH1F, {{200, -10.0f, 10.0f}});
      histos.add("hNsigmaPionTOF_before", "NsigmaPion TOF distribution before", kTH1F, {{200, -10.0f, 10.0f}});
      histos.add("hEta_after", "Eta distribution", kTH1F, {{200, -1.0f, 1.0f}});
      histos.add("hDcaxy_after", "Dcaxy distribution", kTH1F, {{200, -10.0f, 10.0f}});
      histos.add("hDcaz_after", "Dcaz distribution", kTH1F, {{200, -10.0f, 10.0f}});
      histos.add("hNsigmaPionTPC_after", "NsigmaPion TPC distribution", kTH1F, {{200, -10.0f, 10.0f}});
      histos.add("hNsigmaPionTOF_after", "NsigmaPion TOF distribution", kTH1F, {{200, -10.0f, 10.0f}});
    }

    if (QAv0) {
      histos.add("lambda", "Invariant mass of lambda", kTH1F, {lambdaMassAxis});
      histos.add("hDCAV0Daughters", "hDCAV0Daughters", {HistType::kTH1F, {{50, 0.0f, 5.0f}}});
      histos.add("hLT", "hLT", {HistType::kTH1F, {{100, 0.0f, 50.0f}}});
      histos.add("hV0CosPA", "hV0CosPA", {HistType::kTH1F, {{100, 0.95f, 1.f}}});
    }
    histos.add("sigmamassDS", "Invariant mass of sigma", HistType::kTH3F, {multAxis, ptAxis, invMassAxisReso});
    histos.add("sigmamassMIX", "Invariant mass of sigma", HistType::kTH3F, {multAxis, ptAxis, invMassAxisReso});
    histos.add("sigmaDS", "Invariant mass sigma", kTH1F, {invMassAxisReso});
    histos.add("sigmaMix", "Invariant mass sigma", kTH1F, {invMassAxisReso});
  }

  double masslambda = MassLambda0;
  double massPi = MassPiPlus;
  ROOT::Math::PtEtaPhiMVector sigmaVector;

  template <typename T>
  bool selectionTrack(const T& candidate)
  {
    if (iscustomDCAcut &&
        !(candidate.isGlobalTrack() || candidate.isPVContributor() ||
          candidate.itsNCls() > cfgITScluster)) {
      return false;
    }
    if (ismanualDCAcut &&
        !(candidate.isGlobalTrackWoDCA() || candidate.isPVContributor() ||
          std::abs(candidate.dcaXY()) < cfgCutDCAxy ||
          std::abs(candidate.dcaZ()) < cfgCutDCAz ||
          candidate.itsNCls() > cfgITScluster)) {
      return false;
    }
    return true;
  }

  template <typename T>
  bool selectionPID(const T& candidate)
  {
    if (candidate.hasTOF() &&
        (candidate.tofNSigmaPi() * candidate.tofNSigmaPi() +
         candidate.tpcNSigmaPi() * candidate.tpcNSigmaPi()) <
          (nsigmaCutCombined * nsigmaCutCombined)) {
      return true;
    }
    if (!candidate.hasTOF() &&
        std::abs(candidate.tpcNSigmaPi()) < nsigmaCutTPC) {
      return true;
    }
    return false;
  }
  // v0 selection
  template <typename Collision, typename V0>
  bool SelectionV0(Collision const& collision, V0 const& candidate, float /*multiplicity*/)
  {
    if (fabs(candidate.dcav0topv()) > cMaxV0DCA) {
      return false;
    }
    const float tranRad = candidate.v0radius();
    const float dcaDaughv0 = candidate.dcaV0daughters();
    const float cpav0 = candidate.v0cosPA();
    float CtauLambda = candidate.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * masslambda;
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
    if (fabs(CtauLambda) > cMaxV0LifeTime) {
      return false;
    }
    if (fabs(candidate.mLambda() - masslambda) > 10) {
      return false;
    }
    if (QAv0) {
      histos.fill(HIST("hLT"), CtauLambda);
      histos.fill(HIST("lambda"), candidate.mLambda());
      histos.fill(HIST("hDCAV0Daughters"), candidate.dcaV0daughters());
      histos.fill(HIST("hV0CosPA"), candidate.v0cosPA());
    }
    return true;
  }

  template <typename T>
  bool isSelectedV0Daughter(T const& track, float charge, double nsigmaV0Daughter)
  {
    const auto eta = track.eta();
    const auto tpcNClsF = track.tpcNClsFound();
    const auto dcaXY = track.dcaXY();
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
    if (std::abs(nsigmaV0Daughter) > ConfDaughPIDCuts) {
      return false;
    }

    return true;
  }

  ConfigurableAxis axisVertex{"axisVertex", {20, -10, 10}, "vertex axis for bin"};
  ConfigurableAxis axisMultiplicityClass{"axisMultiplicityClass", {1, 0, 100}, "multiplicity percentile for bin"};
  ConfigurableAxis axisMultiplicity{"axisMultiplicity", {2, 0, 100}, "TPC multiplicity  for bin"};
  SliceCache cache;
  Filter posZFilter = (nabs(o2::aod::collision::posZ) < cutzvertex);
  Filter acceptanceFilter = (nabs(aod::track::eta) < cfgCutEta && nabs(aod::track::pt) > cfgCutPT);
  Filter DCAcutFilter = (nabs(aod::track::dcaXY) < cfgCutDCAxy) && (nabs(aod::track::dcaZ) < cfgCutDCAz);

  using EventCandidates = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::MultZeqs, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs>>;
  using TrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullPi, aod::pidTOFFullPi>>;
  using V0TrackCandidate = aod::V0Datas;

  using BinningTypeTPCMultiplicity = ColumnBinningPolicy<aod::collision::PosZ, aod::mult::MultTPC>;
  using BinningTypeCentralityM = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0M>;
  using BinningTypeVertexContributor = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0C>;

  BinningTypeVertexContributor binningOnPositions{{axisVertex, axisMultiplicity}, true};
  Pair<EventCandidates, TrackCandidates, V0TrackCandidate, BinningTypeVertexContributor> pair{binningOnPositions, cfgNoMixedEvents, -1, &cache};

  void processSE(EventCandidates::iterator const& collision, TrackCandidates const& tracks, aod::V0Datas const& V0s, aod::BCs const&)
  {

    if (!collision.sel8()) {
      return;
    }

    std::vector<ROOT::Math::PtEtaPhiMVector> pions, lambda;
    std::vector<int64_t> PionIndex = {};
    std::vector<int64_t> PioncollIndex = {};
    std::vector<int64_t> V0collIndex = {};
    std::vector<int64_t> lambdaPosDaughIndex = {};
    std::vector<int64_t> lambdaNegDaughIndex = {};

    float multiplicity = 0.0f;
    if (cfgMultFT0)
      multiplicity = collision.multZeqFT0A() + collision.multZeqFT0C();
    if (cfgMultFT0 == 0 && cfgCentFT0C == 1)
      multiplicity = collision.centFT0C();
    if (cfgMultFT0 == 0 && cfgCentFT0C == 0)
      multiplicity = collision.centFT0M();

    // Fill the event counter
    rEventSelection.fill(HIST("hVertexZRec"), collision.posZ());
    rEventSelection.fill(HIST("hmult"), multiplicity);

    for (auto track1 : tracks) {

      if (QAbefore) {
        histos.fill(HIST("hNsigmaPionTPC_before"), track1.tpcNSigmaPi());
        histos.fill(HIST("hNsigmaPionTOF_before"), track1.tofNSigmaPi());
      }

      if (!selectionPID(track1))
        continue; // for primary particle PID

      if (!selectionTrack(track1)) {
        continue;
      }

      if (QAafter) {
        histos.fill(HIST("hEta_after"), track1.eta());
        histos.fill(HIST("hDcaxy_after"), track1.dcaXY());
        histos.fill(HIST("hDcaz_after"), track1.dcaZ());
        histos.fill(HIST("hNsigmaPionTPC_after"), track1.tpcNSigmaPi());
        histos.fill(HIST("hNsigmaPionTOF_after"), track1.tofNSigmaPi());
      }

      ROOT::Math::PtEtaPhiMVector temp1(track1.pt(), track1.eta(), track1.phi(), massPi);
      pions.push_back(temp1);
      PionIndex.push_back(track1.globalIndex());
      PioncollIndex.push_back(track1.collisionId());

    } // track loop ends

    for (auto& v0 : V0s) {

      auto postrack = v0.template posTrack_as<TrackCandidates>();
      auto negtrack = v0.template negTrack_as<TrackCandidates>();
      double nTPCSigmaPos[1]{postrack.tpcNSigmaPi()};
      double nTPCSigmaNeg[1]{negtrack.tpcNSigmaPi()};

      if (!isSelectedV0Daughter(postrack, 1, nTPCSigmaPos[0])) {
        continue;
      }
      if (!isSelectedV0Daughter(negtrack, -1, nTPCSigmaNeg[0])) {
        continue;
      }

      if (!SelectionV0(collision, v0, multiplicity)) {
        continue;
      }

      ROOT::Math::PtEtaPhiMVector temp2(v0.pt(), v0.eta(), v0.phi(), masslambda);
      lambda.push_back(temp2);
      V0collIndex.push_back(v0.collisionId());
      lambdaPosDaughIndex.push_back(postrack.globalIndex());
      lambdaNegDaughIndex.push_back(negtrack.globalIndex());
    }

    if (pions.size() != 0 && lambda.size() != 0) {
      for (auto ipion = pions.begin(); ipion != pions.end(); ++ipion) {
        auto i1 = std::distance(pions.begin(), ipion);
        for (auto ilambda = lambda.begin(); ilambda != lambda.end();
             ++ilambda) {
          auto i3 = std::distance(lambda.begin(), ilambda);
          if (PionIndex.at(i1) == lambdaPosDaughIndex.at(i3))
            continue;
          if (PionIndex.at(i1) == lambdaNegDaughIndex.at(i3))
            continue;
          if (PioncollIndex.at(i1) != V0collIndex.at(i3))
            continue;
          sigmaVector = pions.at(i1) + lambda.at(i3);

          if (TMath::Abs(sigmaVector.Rapidity()) < 0.5) {
            histos.fill(HIST("sigmamassDS"), multiplicity, sigmaVector.Pt(), sigmaVector.M());
            histos.fill(HIST("sigmaDS"), sigmaVector.M());
          }
        }
      }
    }
  }

  PROCESS_SWITCH(sigma, processSE, "Process Same event", true);

  void processME(EventCandidates const& /*collisions*/,
                 TrackCandidates const& /*tracks*/, V0TrackCandidate const& /*V0s*/)

  {

    for (auto& [c1, tracks1, c2, tracks2] : pair) {

      if (!c1.sel8()) {
        continue;
      }
      if (!c2.sel8()) {
        continue;
      }

      float multiplicity = 0.0f;
      if (cfgMultFT0)
        multiplicity = c1.multZeqFT0A() + c1.multZeqFT0C();
      if (cfgMultFT0 == 0 && cfgCentFT0C == 1)
        multiplicity = c1.centFT0C();
      if (cfgMultFT0 == 0 && cfgCentFT0C == 0)
        multiplicity = c1.centFT0M();

      for (auto& [t1, t2] : o2::soa::combinations(
             o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {

        if (!(selectionTrack(t1)))
          continue;
        if (!(selectionPID(t1)))
          continue;
        if (!SelectionV0(c2, t2, multiplicity))
          continue;

        auto postrack = t2.template posTrack_as<TrackCandidates>();
        auto negtrack = t2.template negTrack_as<TrackCandidates>();
        double nTPCSigmaPos[1]{postrack.tpcNSigmaPi()};
        double nTPCSigmaNeg[1]{negtrack.tpcNSigmaPi()};

        if (!isSelectedV0Daughter(postrack, 1, nTPCSigmaPos[0])) {
          continue;
        }
        if (!isSelectedV0Daughter(negtrack, -1, nTPCSigmaNeg[0])) {
          continue;
        }

        TLorentzVector pi;
        pi.SetPtEtaPhiM(t1.pt(), t1.eta(), t1.phi(), massPi);
        TLorentzVector lm;
        lm.SetPtEtaPhiM(t2.pt(), t2.eta(), t2.phi(), masslambda);

        TLorentzVector lambdamix = pi + lm;

        if (TMath::Abs(lambdamix.Rapidity()) < 0.5) {
          histos.fill(HIST("sigmamassMIX"), multiplicity, lambdamix.Pt(), lambdamix.M());
          histos.fill(HIST("sigmaMix"), lambdamix.M());
        }
      }
    }
  }
  PROCESS_SWITCH(sigma, processME, "Process Mixed event", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<sigma>(cfgc)};
}
