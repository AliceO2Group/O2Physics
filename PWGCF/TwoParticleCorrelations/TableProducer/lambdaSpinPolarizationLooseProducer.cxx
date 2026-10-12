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

/// \file lambdaSpinPolarizationLooseProducer.cxx
/// \brief Lambda candidates with loose selections, written as derived data for lambdaSpinPolarization
/// \author Subhadeep Roy <subhadeep.roy@cern.ch>

#include "PWGCF/TwoParticleCorrelations/DataModel/LambdaSpinPolarizationDerived.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::constants::physics;

enum LooseEventCutFlow { kLeAll = 1,
                         kLeTrigger,
                         kLeTvx,
                         kLeTFBorder,
                         kLeItsRofBorder,
                         kLeItsTpcVtx,
                         kLeNoSameBunchPileup,
                         kLeGoodZvtxFT0vsPV,
                         kLeGoodItsLayers,
                         kLeCentrality,
                         kLeVz,
                         kLeOneV0,
                         kLeTwoV0,
                         kLeWritten,
                         kNLooseEventCutFlow };

enum LooseV0Count { kLvAll = 1,
                    kLvV0Type,
                    kLvAmbiguousTracks,
                    kLvKShortRejected,
                    kLvKinematics,
                    kLvTopology,
                    kLvNoHypothesis,
                    kLvLambdaRow,
                    kLvAntiLambdaRow,
                    kLvBothRows,
                    kNLooseV0Count };

enum LooseHypothesis { kLooseLambda = 0,
                       kLooseAntiLambda,
                       kNLooseHypotheses };

static constexpr std::size_t NV0ForPair = 2;

enum LooseCentEstimator { kLooseCentFT0M = 0,
                          kLooseCentFT0C };

struct LambdaSpinPolarizationLooseProducer {

  Produces<aod::LambdaLooseCollisions> lambdaLooseCollisionTable;
  Produces<aod::LambdaLooseTracks> lambdaLooseTrackTable;

  Configurable<int> cCentEstimator{"cCentEstimator", 0, "Centrality Estimator: 0=FT0M, 1=FT0C"};
  Configurable<float> cMinZVtx{"cMinZVtx", -10.0, "Min VtxZ (cm)"};
  Configurable<float> cMaxZVtx{"cMaxZVtx", 10.0, "Max VtxZ (cm)"};
  Configurable<float> cMinMult{"cMinMult", 0.0, "Min centrality percentile"};
  Configurable<float> cMaxMult{"cMaxMult", 100.0, "Max centrality percentile"};
  Configurable<bool> cSel8Trig{"cSel8Trig", true, "Sel8 (T0A+T0C) Run3"};
  Configurable<bool> cTriggerTvxSel{"cTriggerTvxSel", false, "TVX Trigger Selection"};
  Configurable<bool> cTFBorder{"cTFBorder", false, "Timeframe Border Selection"};
  Configurable<bool> cNoItsROBorder{"cNoItsROBorder", false, "No ITSRO Border Cut"};
  Configurable<bool> cItsTpcVtx{"cItsTpcVtx", false, "ITS+TPC Vertex Selection"};
  Configurable<bool> cPileupReject{"cPileupReject", false, "Pileup rejection"};
  Configurable<bool> cZVtxTimeDiff{"cZVtxTimeDiff", false, "z-vtx time diff selection"};
  Configurable<bool> cIsGoodITSLayers{"cIsGoodITSLayers", false, "Good ITS Layers All"};

  Configurable<float> cTrackMinPt{"cTrackMinPt", 0.15, "p_{T} minimum"};
  Configurable<float> cTrackMaxPt{"cTrackMaxPt", 999.0, "p_{T} maximum"};
  Configurable<float> cTrackEtaCut{"cTrackEtaCut", 0.8, "Pseudorapidity cut"};
  Configurable<int> cMinTpcCrossedRows{"cMinTpcCrossedRows", 50, "TPC Min Crossed Rows (loosest of the variations)"};
  Configurable<float> cMinTpcCROverCls{"cMinTpcCROverCls", -999, "TPC Min CR/Findable Cls"};
  Configurable<float> cMaxTpcSharedClusters{"cMaxTpcSharedClusters", 999, "TPC Max Shared Clusters"};
  Configurable<float> cMaxChi2Tpc{"cMaxChi2Tpc", 999, "Max TPC Chi2/ndf"};
  Configurable<double> cTpcNsigmaCut{"cTpcNsigmaCut", 5.0, "TPC nSigma PID cut (loosest of the variations)"};
  Configurable<bool> cRemoveAmbiguousTracks{"cRemoveAmbiguousTracks", false, "Remove Ambiguous Tracks"};

  Configurable<double> cMinDcaProtonToPV{"cMinDcaProtonToPV", 0.05, "Min proton DCA to PV (cm)"};
  Configurable<double> cMinDcaPionToPV{"cMinDcaPionToPV", 0.1, "Min pion DCA to PV (cm)"};
  Configurable<double> cMinV0DcaDaughters{"cMinV0DcaDaughters", 0., "Min DCA between V0 daughters"};
  Configurable<double> cMaxV0DcaDaughters{"cMaxV0DcaDaughters", 1.4, "Max DCA between V0 daughters"};
  Configurable<double> cMinDcaV0ToPV{"cMinDcaV0ToPV", 0.0, "Min DCA V0 to PV"};
  Configurable<double> cMaxDcaV0ToPV{"cMaxDcaV0ToPV", 0.2, "Max DCA V0 to PV (loosest of the variations)"};
  Configurable<double> cMinV0TransRadius{"cMinV0TransRadius", 0.5, "Min V0 decay radius (cm)"};
  Configurable<double> cMaxV0TransRadius{"cMaxV0TransRadius", 999.0, "Max V0 decay radius (cm)"};
  Configurable<double> cMinV0CTau{"cMinV0CTau", 0.0, "Min cTau (cm)"};
  Configurable<double> cMaxV0CTau{"cMaxV0CTau", 36.0, "Max cTau (cm) (loosest of the variations)"};
  Configurable<double> cMinV0CosPA{"cMinV0CosPA", 0.98, "Min V0 cos(PA) (loosest of the variations)"};
  Configurable<double> cKshortRejMassWindow{"cKshortRejMassWindow", 0.001, "K0s mass rejection window (narrowest of the variations)"};
  Configurable<bool> cKshortRejFlag{"cKshortRejFlag", true, "K0s mass rejection flag"};

  Configurable<float> cMinV0Mass{"cMinV0Mass", 1.09, "V0 Mass Min"};
  Configurable<float> cMaxV0Mass{"cMaxV0Mass", 1.14, "V0 Mass Max"};
  Configurable<float> cMinV0Pt{"cMinV0Pt", 0.6, "Minimum V0 pT"};
  Configurable<float> cMaxV0Pt{"cMaxV0Pt", 3.0, "Maximum V0 pT"};
  Configurable<float> cMaxV0Rap{"cMaxV0Rap", 0.5, "|rap| cut"};
  Configurable<bool> cDoEtaAnalysis{"cDoEtaAnalysis", false, "Do Eta Analysis"};
  Configurable<bool> cV0TypeSelFlag{"cV0TypeSelFlag", true, "V0 Type Selection Flag"};
  Configurable<int> cV0TypeSelection{"cV0TypeSelection", 1, "V0 Type Selection"};

  Configurable<int> cMinV0PerEvent{"cMinV0PerEvent", 2, "Collisions with fewer V0s passing the loose selection are not written"};

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  struct LooseCandidate {
    float px, py, pz, pt, eta, phi, rap, mass, massOther, prPx, prPy, prPz;
    int32_t posTrackId, negTrackId;
    int8_t v0Type;
    float cosPA, dcaDau, dcaV0ToPV, cTau, mK0Short;
    float posNSigmaPr, posNSigmaPi, negNSigmaPr, negNSigmaPi;
    int16_t posTpcCrossedRows, negTpcCrossedRows;
  };
  std::vector<LooseCandidate> candidates;

  float cent = 0., mult = 0.;

  void init(InitContext const&)
  {
    if (cMinV0PerEvent.value < 1) {
      LOGF(fatal, "cMinV0PerEvent must be at least 1, it is %d", cMinV0PerEvent.value);
    }
    const AxisSpec axisV0Mass(200, 1.08, 1.18, "M_{p#pi} (GeV/#it{c}^{2})");
    const AxisSpec axisV0Pt(100., 0., 10., "p_{T} (GeV/#it{c})");
    const AxisSpec axisNV0(20, -0.5, 19.5, "V0s written per collision");

    auto hCutFlow = histos.add<TH1>("Events/hEventCutFlow", "event cut flow;;collisions", kTH1D, {{kNLooseEventCutFlow - 1, 0.5, static_cast<double>(kNLooseEventCutFlow) - 0.5}});
    hCutFlow->GetXaxis()->SetBinLabel(kLeAll, "all");
    hCutFlow->GetXaxis()->SetBinLabel(kLeTrigger, "trigger (sel8)");
    hCutFlow->GetXaxis()->SetBinLabel(kLeTvx, "TVX");
    hCutFlow->GetXaxis()->SetBinLabel(kLeTFBorder, "no TF border");
    hCutFlow->GetXaxis()->SetBinLabel(kLeItsRofBorder, "no ITS ROF border");
    hCutFlow->GetXaxis()->SetBinLabel(kLeItsTpcVtx, "ITS-TPC vertex");
    hCutFlow->GetXaxis()->SetBinLabel(kLeNoSameBunchPileup, "no same-bunch pileup");
    hCutFlow->GetXaxis()->SetBinLabel(kLeGoodZvtxFT0vsPV, "good z_{vtx} FT0 vs PV");
    hCutFlow->GetXaxis()->SetBinLabel(kLeGoodItsLayers, "good ITS layers");
    hCutFlow->GetXaxis()->SetBinLabel(kLeCentrality, "centrality range");
    hCutFlow->GetXaxis()->SetBinLabel(kLeVz, "V_{z} range");
    hCutFlow->GetXaxis()->SetBinLabel(kLeOneV0, "#geq 1 V0 with a loose candidate");
    hCutFlow->GetXaxis()->SetBinLabel(kLeTwoV0, "#geq 2 V0s with a loose candidate");
    hCutFlow->GetXaxis()->SetBinLabel(kLeWritten, "written");

    auto hV0 = histos.add<TH1>("Tracks/hV0Count", "V0s;;V0s", kTH1D, {{kNLooseV0Count - 1, 0.5, static_cast<double>(kNLooseV0Count) - 0.5}});
    hV0->GetXaxis()->SetBinLabel(kLvAll, "all");
    hV0->GetXaxis()->SetBinLabel(kLvV0Type, "rejected: V0 type");
    hV0->GetXaxis()->SetBinLabel(kLvAmbiguousTracks, "rejected: ambiguous tracks");
    hV0->GetXaxis()->SetBinLabel(kLvKShortRejected, "rejected: K^{0}_{S} mass");
    hV0->GetXaxis()->SetBinLabel(kLvKinematics, "rejected: p_{T} / y");
    hV0->GetXaxis()->SetBinLabel(kLvTopology, "rejected: topology");
    hV0->GetXaxis()->SetBinLabel(kLvNoHypothesis, "rejected: no hypothesis");
    hV0->GetXaxis()->SetBinLabel(kLvLambdaRow, "#Lambda candidate");
    hV0->GetXaxis()->SetBinLabel(kLvAntiLambdaRow, "#bar{#Lambda} candidate");
    hV0->GetXaxis()->SetBinLabel(kLvBothRows, "both candidates");

    histos.add("Tracks/hLambdaMassVsPt", "#Lambda candidates", kTH2F, {axisV0Mass, axisV0Pt});
    histos.add("Tracks/hAntiLambdaMassVsPt", "#bar{#Lambda} candidates", kTH2F, {axisV0Mass, axisV0Pt});
    histos.add("Events/hNV0Written", "written collisions", kTH1F, {axisNV0});
  }

  template <typename C>
  bool selCollision(C const& col)
  {
    histos.fill(HIST("Events/hEventCutFlow"), kLeAll);
    if (cCentEstimator == kLooseCentFT0M) {
      cent = col.centFT0M();
    } else if (cCentEstimator == kLooseCentFT0C) {
      cent = col.centFT0C();
    }
    if (cSel8Trig && !col.sel8()) {
      return false;
    }
    histos.fill(HIST("Events/hEventCutFlow"), kLeTrigger);

    if (cTriggerTvxSel && !col.selection_bit(aod::evsel::kIsTriggerTVX)) {
      return false;
    }
    histos.fill(HIST("Events/hEventCutFlow"), kLeTvx);
    if (cTFBorder && !col.selection_bit(aod::evsel::kNoTimeFrameBorder)) {
      return false;
    }
    histos.fill(HIST("Events/hEventCutFlow"), kLeTFBorder);
    if (cNoItsROBorder && !col.selection_bit(aod::evsel::kNoITSROFrameBorder)) {
      return false;
    }
    histos.fill(HIST("Events/hEventCutFlow"), kLeItsRofBorder);
    if (cItsTpcVtx && !col.selection_bit(aod::evsel::kIsVertexITSTPC)) {
      return false;
    }
    histos.fill(HIST("Events/hEventCutFlow"), kLeItsTpcVtx);
    if (cPileupReject && !col.selection_bit(aod::evsel::kNoSameBunchPileup)) {
      return false;
    }
    histos.fill(HIST("Events/hEventCutFlow"), kLeNoSameBunchPileup);
    if (cZVtxTimeDiff && !col.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return false;
    }
    histos.fill(HIST("Events/hEventCutFlow"), kLeGoodZvtxFT0vsPV);
    if (cIsGoodITSLayers && !col.selection_bit(aod::evsel::kIsGoodITSLayersAll)) {
      return false;
    }
    histos.fill(HIST("Events/hEventCutFlow"), kLeGoodItsLayers);

    if (cent <= cMinMult || cent >= cMaxMult) {
      return false;
    }
    histos.fill(HIST("Events/hEventCutFlow"), kLeCentrality);
    if (col.posZ() <= cMinZVtx || col.posZ() >= cMaxZVtx) {
      return false;
    }
    histos.fill(HIST("Events/hEventCutFlow"), kLeVz);

    mult = col.multNTracksPV();
    return true;
  }

  bool kinCutSelection(float const& pt, float const& rap, float const& ptMin, float const& ptMax, float const& rapMax)
  {
    return pt > ptMin && pt < ptMax && rap < rapMax;
  }

  template <typename T>
  bool selTrack(T const& track)
  {
    if (!kinCutSelection(track.pt(), std::abs(track.eta()), cTrackMinPt, cTrackMaxPt, cTrackEtaCut)) {
      return false;
    }
    if (track.tpcNClsCrossedRows() <= cMinTpcCrossedRows) {
      return false;
    }
    if (track.tpcCrossedRowsOverFindableCls() < cMinTpcCROverCls) {
      return false;
    }
    if (track.tpcNClsShared() > cMaxTpcSharedClusters) {
      return false;
    }
    if (track.tpcChi2NCl() > cMaxChi2Tpc) {
      return false;
    }
    return true;
  }

  template <typename V, typename T>
  bool selDaughterTracks(V const& v0, T const&, LooseHypothesis const& v0Type)
  {
    auto posTrack = v0.template posTrack_as<T>();
    auto negTrack = v0.template negTrack_as<T>();
    if (!selTrack(posTrack) || !selTrack(negTrack)) {
      return false;
    }

    float dcaProton = 0., dcaPion = 0.;
    if (v0Type == kLooseLambda) {
      dcaProton = std::abs(v0.dcapostopv());
      dcaPion = std::abs(v0.dcanegtopv());
    } else if (v0Type == kLooseAntiLambda) {
      dcaPion = std::abs(v0.dcapostopv());
      dcaProton = std::abs(v0.dcanegtopv());
    }
    if (dcaProton < cMinDcaProtonToPV || dcaPion < cMinDcaPionToPV) {
      return false;
    }
    return true;
  }

  template <typename C, typename V, typename T>
  bool topoCutSelection(C const& col, V const& v0, T const&)
  {
    if (v0.dcaV0daughters() <= cMinV0DcaDaughters || v0.dcaV0daughters() >= cMaxV0DcaDaughters) {
      return false;
    }
    if (v0.dcav0topv() <= cMinDcaV0ToPV || v0.dcav0topv() >= cMaxDcaV0ToPV) {
      return false;
    }
    if (v0.v0radius() <= cMinV0TransRadius || v0.v0radius() >= cMaxV0TransRadius) {
      return false;
    }

    float ctau = v0.distovertotmom(col.posX(), col.posY(), col.posZ()) * MassLambda0;
    if (ctau <= cMinV0CTau || ctau >= cMaxV0CTau) {
      return false;
    }
    if (v0.v0cosPA() <= cMinV0CosPA) {
      return false;
    }
    return true;
  }

  template <LooseHypothesis part, typename T>
  bool selLambdaDauWithTpcPid(T const& postrack, T const& negtrack)
  {
    float tpcNSigmaPr = 0., tpcNSigmaPi = 0.;
    if constexpr (part == kLooseLambda) {
      tpcNSigmaPr = postrack.tpcNSigmaPr();
      tpcNSigmaPi = negtrack.tpcNSigmaPi();
    } else {
      tpcNSigmaPr = negtrack.tpcNSigmaPr();
      tpcNSigmaPi = postrack.tpcNSigmaPi();
    }
    return (std::abs(tpcNSigmaPr) < cTpcNsigmaCut && std::abs(tpcNSigmaPi) < cTpcNsigmaCut);
  }

  template <typename V, typename T>
  bool hasAmbiguousDaughters(V const& v0, T const&)
  {
    auto posTrack = v0.template posTrack_as<T>();
    auto negTrack = v0.template negTrack_as<T>();
    auto posCC = posTrack.compatibleCollIds();
    auto negCC = negTrack.compatibleCollIds();
    if (posCC.size() > 1 || negCC.size() > 1) {
      return true;
    }
    if ((posCC.size() != 0 && posCC[0] != posTrack.collisionId()) ||
        (negCC.size() != 0 && negCC[0] != negTrack.collisionId())) {
      return true;
    }
    return false;
  }

  using CollisionsRun3 = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::CentFT0Cs, aod::PVMults>;
  using Tracks = soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra,
                           aod::TracksDCA, aod::pidTPCPi, aod::pidTPCPr,
                           aod::TrackCompColls>;

  void processData(CollisionsRun3::iterator const& collision, aod::BCsWithTimestamps const&, aod::V0Datas const& V0s, Tracks const& tracks)
  {
    if (!selCollision(collision)) {
      return;
    }

    candidates.clear();
    std::size_t nV0 = 0;
    for (auto const& v0 : V0s) {
      histos.fill(HIST("Tracks/hV0Count"), kLvAll);
      if (cV0TypeSelFlag && v0.v0Type() != cV0TypeSelection) {
        histos.fill(HIST("Tracks/hV0Count"), kLvV0Type);
        continue;
      }
      if (cRemoveAmbiguousTracks && hasAmbiguousDaughters(v0, tracks)) {
        histos.fill(HIST("Tracks/hV0Count"), kLvAmbiguousTracks);
        continue;
      }
      if (cKshortRejFlag && (std::abs(v0.mK0Short() - MassK0Short) <= cKshortRejMassWindow)) {
        histos.fill(HIST("Tracks/hV0Count"), kLvKShortRejected);
        continue;
      }
      const float rap = cDoEtaAnalysis ? std::abs(v0.eta()) : std::abs(v0.yLambda());
      if (!kinCutSelection(v0.pt(), rap, cMinV0Pt, cMaxV0Pt, cMaxV0Rap)) {
        histos.fill(HIST("Tracks/hV0Count"), kLvKinematics);
        continue;
      }
      if (!topoCutSelection(collision, v0, tracks)) {
        histos.fill(HIST("Tracks/hV0Count"), kLvTopology);
        continue;
      }

      auto postrack = v0.posTrack_as<Tracks>();
      auto negtrack = v0.negTrack_as<Tracks>();
      const std::array<float, kNLooseHypotheses> masses = {v0.mLambda(), v0.mAntiLambda()};
      const std::array<bool, kNLooseHypotheses> pid = {selLambdaDauWithTpcPid<kLooseLambda>(postrack, negtrack), selLambdaDauWithTpcPid<kLooseAntiLambda>(postrack, negtrack)};
      const float ctau = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * MassLambda0;

      std::size_t nRows = 0;
      for (auto const& hyp : {kLooseLambda, kLooseAntiLambda}) {
        const bool inMassWindow = masses[hyp] > cMinV0Mass && masses[hyp] < cMaxV0Mass;
        if (!inMassWindow || !pid[hyp]) {
          continue;
        }
        if (!selDaughterTracks(v0, tracks, hyp)) {
          continue;
        }
        const bool isLambda = hyp == kLooseLambda;
        candidates.push_back(LooseCandidate{.px = v0.px(),
                                            .py = v0.py(),
                                            .pz = v0.pz(),
                                            .pt = v0.pt(),
                                            .eta = v0.eta(),
                                            .phi = v0.phi(),
                                            .rap = v0.yLambda(),
                                            .mass = masses[hyp],
                                            .massOther = masses[isLambda ? kLooseAntiLambda : kLooseLambda],
                                            .prPx = isLambda ? v0.pxpos() : v0.pxneg(),
                                            .prPy = isLambda ? v0.pypos() : v0.pyneg(),
                                            .prPz = isLambda ? v0.pzpos() : v0.pzneg(),
                                            .posTrackId = static_cast<int32_t>(postrack.index()),
                                            .negTrackId = static_cast<int32_t>(negtrack.index()),
                                            .v0Type = static_cast<int8_t>(hyp),
                                            .cosPA = v0.v0cosPA(),
                                            .dcaDau = v0.dcaV0daughters(),
                                            .dcaV0ToPV = v0.dcav0topv(),
                                            .cTau = ctau,
                                            .mK0Short = v0.mK0Short(),
                                            .posNSigmaPr = postrack.tpcNSigmaPr(),
                                            .posNSigmaPi = postrack.tpcNSigmaPi(),
                                            .negNSigmaPr = negtrack.tpcNSigmaPr(),
                                            .negNSigmaPi = negtrack.tpcNSigmaPi(),
                                            .posTpcCrossedRows = postrack.tpcNClsCrossedRows(),
                                            .negTpcCrossedRows = negtrack.tpcNClsCrossedRows()});
        if (isLambda) {
          histos.fill(HIST("Tracks/hV0Count"), kLvLambdaRow);
          histos.fill(HIST("Tracks/hLambdaMassVsPt"), masses[hyp], v0.pt());
        } else {
          histos.fill(HIST("Tracks/hV0Count"), kLvAntiLambdaRow);
          histos.fill(HIST("Tracks/hAntiLambdaMassVsPt"), masses[hyp], v0.pt());
        }
        ++nRows;
      }
      if (nRows == 0) {
        histos.fill(HIST("Tracks/hV0Count"), kLvNoHypothesis);
        continue;
      }
      if (nRows == kNLooseHypotheses) {
        histos.fill(HIST("Tracks/hV0Count"), kLvBothRows);
      }
      ++nV0;
    }

    if (nV0 >= 1) {
      histos.fill(HIST("Events/hEventCutFlow"), kLeOneV0);
    }
    if (nV0 >= NV0ForPair) {
      histos.fill(HIST("Events/hEventCutFlow"), kLeTwoV0);
    }
    if (std::cmp_less(nV0, cMinV0PerEvent.value)) {
      return;
    }
    histos.fill(HIST("Events/hEventCutFlow"), kLeWritten);
    histos.fill(HIST("Events/hNV0Written"), nV0);

    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    lambdaLooseCollisionTable(cent, mult, collision.posX(), collision.posY(), collision.posZ(), bc.timestamp());
    for (auto const& c : candidates) {
      lambdaLooseTrackTable(lambdaLooseCollisionTable.lastIndex(), c.px, c.py, c.pz, c.pt, c.eta, c.phi, c.rap, c.mass, c.massOther, c.prPx, c.prPy, c.prPz,
                            c.posTrackId, c.negTrackId, c.v0Type, c.cosPA, c.dcaDau, c.dcaV0ToPV, c.cTau, c.mK0Short,
                            c.posNSigmaPr, c.posNSigmaPi, c.negNSigmaPr, c.negNSigmaPi, c.posTpcCrossedRows, c.negTpcCrossedRows);
    }
  }
  PROCESS_SWITCH(LambdaSpinPolarizationLooseProducer, processData, "Run 3 data (or MC reconstructed as data)", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<LambdaSpinPolarizationLooseProducer>(cfgc)};
}
