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
#include "TF1.h"

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

  Configurable<bool> QAv0{"QAv0", false, "QAv0"};
  Configurable<bool> QAPID{"QAPID", false, "QAPID"};
  Configurable<bool> QAv0_daughters{"QAv0_daughters", false, "QA of v0 daughters"};
  Configurable<bool> inv_mass1D{"inv_mass1D", false, "1D invariant mass histograms"};
  Configurable<bool> DCAv0topv{"DCAv0topv", false, "DCA V0 to PV"};
  Configurable<bool> armcut{"armcut", true, "arm cut"};

  // Configurable for event selection
  Configurable<float> cutzvertex{"cutzvertex", 10.0f, "Accepted z-vertex range (cm)"};
  Configurable<float> cfgETAcut{"cfgETAcut", 0.8f, "Track ETA cut"};
  Configurable<bool> timFrameEvsel{"timFrameEvsel", false, "TPC Time frame boundary cut"};
  Configurable<bool> piluprejection{"piluprejection", false, "Pileup rejection"};
  Configurable<bool> goodzvertex{"goodzvertex", false, "removes collisions with large differences between z of PV by tracks and z of PV from FT0 A-C time difference."};
  Configurable<bool> itstpctracks{"itstpctracks", false, "selects collisions with at least one ITS-TPC track,"};
  Configurable<bool> additionalEvsel{"additionalEvsel", false, "Additional event selcection"};

  // Configurable parameters for V0 selection
  Configurable<float> ConfV0DCADaughMax{"ConfV0DCADaughMax", 1.0f, "DCA b/w V0 daughters"};
  Configurable<float> v0setting_dcapostopv{"v0setting_dcapostopv", 0.06, "DCA Pos To PV"};
  Configurable<float> v0setting_dcanegtopv{"v0setting_dcanegtopv", 0.06, "DCA Neg To PV"};
  Configurable<double> cMaxV0DCA{"cMaxV0DCA", 1, "DCA V0 to PV"};
  // Configurable<bool> isStandarv0{"isStandarv0", false, "Standard V0"};
  // Configurable<float> ConfDaughDCAMin{"ConfDaughDCAMin", 0.06f, "V0 Daugh sel:  Max. DCA Daugh to PV (cm)"}; // same as DCA pos to pv and neg to pv

  Configurable<float> ConfV0PtMin{"ConfV0PtMin", 0.f, "Minimum transverse momentum of V0"};
  Configurable<float> ConfV0CPAMin{"ConfV0CPAMin", 0.97f, "Minimum CPA of V0"};
  Configurable<float> ConfV0TranRadV0Min{"ConfV0TranRadV0Min", 0.5f, "Minimum transverse radius"};
  Configurable<float> ConfV0TranRadV0Max{"ConfV0TranRadV0Max", 200.f, "Maximum transverse radius"};
  Configurable<double> cMaxV0LifeTime{"cMaxV0LifeTime", 15, "Maximum V0 life time"};
  Configurable<double> cSigmaMassKs0{"cSigmaMassKs0", 4, "n Sigma cut on Ks0 mass (Mass (Ks) - cSigmaMassKs0*cWidthKs0)"};
  Configurable<double> cWidthKs0{"cWidthKs0", 0.005, "Width of KS0"};
  Configurable<float> ConfDaughEta{"ConfDaughEta", 0.8f, "V0 Daugh sel: max eta"};
  Configurable<float> ConfDaughTPCnclsMin{"ConfDaughTPCnclsMin", 70.f, "V0 Daugh sel: Min. nCls TPC"};
  Configurable<float> ConfDaughPIDCuts{"ConfDaughPIDCuts", 5, "PID selections for KS0 daughters"};
  Configurable<float> Confarmcut{"Confarmcut", 0.2f, "Armenteros cut"};
  Configurable<float> ConfKsrapidity{"ConfKsrapidity", 0.5f, "Rapidity cut on K0s"};
  // Configurable<float> lowmasscutks0{"lowmasscutks0", 0.497 - 4 * 0.005, "Low mass cut on K0s"};
  // Configurable<float> highmasscutks0{"highmasscutks0", 0.497 + 4 * 0.005, "High mass cut on K0s"};

  // Configurable for track selection and multiplicity
  Configurable<float> cfgPTcut{"cfgPTcut", 0.2f, "Track PT cut"};
  Configurable<int> cfgNmixedEvents{"cfgNmixedEvents", 5, "Number of mixed events"};
  Configurable<bool> cfgMultFOTM{"cfgMultFOTM", true, "Use FOTM multiplicity if pp else use 0 here for PbPb (FT0C)"};
  ConfigurableAxis binsCent{"binsCent", {VARIABLE_WIDTH, 0., 5., 10., 30., 50., 70., 100., 110., 150.}, "Binning of the centrality axis"};

  // Other cuts on Ks and glueball
  Configurable<bool> rapidityks{"rapidityks", true, "rapidity cut on K0s"};
  Configurable<bool> masslambda{"masslambda", false, "mass under lambda hypothesis"};
  Configurable<float> competingcascrejlambda{"competingcascrejlambda", 4.3, "rejecting competing cascade lambda"};
  Configurable<float> competingcascrejlambdaanti{"competingcascrejlambdaanti", 4.3, "rejecting competing cascade anti-lambda"};
  Configurable<int> tpcCrossedrows{"tpcCrossedrows", 70, "TPC crossed rows"};
  Configurable<float> tpcCrossedrowsOverfcls{"tpcCrossedrowsOverfcls", 0.8, "TPC crossed rows over findable clusters"};

  // Mass and pT axis as configurables
  Configurable<float> cPtMin{"cPtMin", 0.0f, "Minimum pT"};
  Configurable<float> cPtMax{"cPtMax", 15.0f, "Maximum pT"};
  Configurable<int> cPtBins{"cPtBins", 150, "Number of pT bins"};
  Configurable<float> cMassMin{"cMassMin", 0.9f, "Minimum mass of glueball"};
  Configurable<float> cMassMax{"cMassMax", 2.4f, "Maximum mass of glueball"};
  Configurable<int> cMassBins{"cMassBins", 150, "Number of mass bins for glueball"};
  Configurable<float> ksMassMin{"ksMassMin", 0.45f, "Minimum mass of K0s"};
  Configurable<float> ksMassMax{"ksMassMax", 0.55f, "Maximum mass of K0s"};
  Configurable<int> ksMassBins{"ksMassBins", 200, "Number of mass bins for K0s"};

  // Event selection cuts - Alex (Temporary, need to fix!)
  TF1* fMultPVCutLow = nullptr;
  TF1* fMultPVCutHigh = nullptr;
  TF1* fMultCutLow = nullptr;
  TF1* fMultCutHigh = nullptr;
  TF1* fMultMultPVCut = nullptr;

  void init(InitContext const&)
  {
    // Axes
    AxisSpec K0ShortMassAxis = {ksMassBins, ksMassMin, ksMassMax, "#it{M}_{inv} [GeV/#it{c}^{2}]"};
    AxisSpec glueballMassAxis = {cMassBins, cMassMin, cMassMax, "#it{M}_{inv} [GeV/#it{c}^{2}]"};
    AxisSpec vertexZAxis = {100, -15.f, 15.f, "vrtx_{Z} [cm]"}; // for histogram
    AxisSpec ptAxis = {cPtBins, cPtMin, cPtMax, "#it{p}_{T} (GeV/#it{c})"};
    // AxisSpec multiplicityAxis = {110, 0.0f, 150.0f, "Multiplicity Axis"};
    AxisSpec multiplicityAxis = {binsCent, "Multiplicity Axis"};

    // Histograms
    // Event selection
    rEventSelection.add("hVertexZRec", "hVertexZRec", {HistType::kTH1F, {vertexZAxis}});
    rEventSelection.add("hmultiplicity", "hmultiplicity", {HistType::kTH1F, {{150, 0.0f, 150.0f}}});

    if (inv_mass1D) {
      hglue.add("h1glueInvMassDS", "h1glueInvMassDS", kTH1F, {glueballMassAxis});
      hglue.add("h1glueInvMassME", "h1glueInvMassME", kTH1F, {glueballMassAxis});
    }
    hglue.add("h3glueInvMassDS", "h3glueInvMassDS", kTHnSparseF, {multiplicityAxis, ptAxis, glueballMassAxis}, true);
    hglue.add("h3glueInvMassME", "h3glueInvMassME", kTHnSparseF, {multiplicityAxis, ptAxis, glueballMassAxis}, true);
    hglue.add("heventscheck", "heventscheck", kTH1F, {{9, 0, 9}});

    // K0s topological/PID cuts
    if (QAv0) {
      // Invariant Mass
      rKzeroShort.add("hMassK0Shortbefore", "hMassK0Shortbefore", kTHnSparseF, {K0ShortMassAxis, ptAxis});
      rKzeroShort.add("hMassK0ShortSelected", "hMassK0ShortSelected", kTHnSparseF, {K0ShortMassAxis, ptAxis});
      // Topological histograms (after the selection)
      rKzeroShort.add("hDCAV0Daughters", "DCA between v0 daughters", {HistType::kTH1F, {{60, -3.0f, 3.0f}}});
      rKzeroShort.add("hV0CosPA", "hV0CosPA", {HistType::kTH1F, {{100, 0.96f, 1.1f}}});
      rKzeroShort.add("hLT", "hLT", {HistType::kTH1F, {{100, 0.0f, 50.0f}}});
      rKzeroShort.add("Mass_lambda", "Mass under lambda hypothesis", kTH1F, {glueballMassAxis});
      rKzeroShort.add("mass_AntiLambda", "Mass under anti-lambda hypothesis", kTH1F, {glueballMassAxis});
      rKzeroShort.add("mass_Gamma", "Mass under Gamma hypothesis", kTH1F, {glueballMassAxis});
      // rKzeroShort.add("mass_Hypertriton", "Mass under hypertriton hypothesis", kTH1F, {glueballMassAxis});
      // rKzeroShort.add("mass_AnitHypertriton", "Mass under anti-hypertriton hypothesis", kTH1F, {glueballMassAxis});
      rKzeroShort.add("rapidity", "Rapidity distribution", kTH1F, {{100, -1.0f, 1.0f}});
      rKzeroShort.add("hv0radius", "hv0radius", kTH1F, {{100, 0.0f, 200.0f}});
      rKzeroShort.add("hDCApostopv", "DCA positive daughter to PV", kTH1F, {{1000, -10.0f, 10.0f}});
      rKzeroShort.add("hDCAnegtopv", "DCA negative daughter to PV", kTH1F, {{1000, -10.0f, 10.0f}});
      rKzeroShort.add("hDCAv0topv", "DCA V0 to PV", kTH1F, {{60, -3.0f, 3.0f}});
      rKzeroShort.add("halpha", "Armenteros alpha", kTH1F, {{100, -5.0f, 5.0f}});
      rKzeroShort.add("hqtarm", "qtarm", kTH1F, {{100, 0.0f, 1.0f}});
      rKzeroShort.add("hpsipair", "psi pair angle", kTH1F, {{100, -5.0f, 5.0f}});

      // // Topological histograms (before the selection)
      // rKzeroShort.add("hDCAV0Daughters_before", "DCA between v0 daughters before the selection", {HistType::kTH1F, {{60, -3.0f, 3.0f}}});
      // rKzeroShort.add("hV0CosPA_before", "hV0CosPA_before", {HistType::kTH1F, {{200, 0.91f, 1.1f}}});
      // rKzeroShort.add("hLT_before", "hLT_before", {HistType::kTH1F, {{100, 0.0f, 50.0f}}});
    }
    if (QAPID) {
      rKzeroShort.add("hNSigmaPosPionK0s_before", "hNSigmaPosPionK0s_before", {HistType::kTH2F, {{ptAxis}, {100, -5.f, 5.f}}});
      rKzeroShort.add("hNSigmaPosPionK0s_after", "hNSigmaPosPionK0s_after", {HistType::kTH2F, {{ptAxis}, {100, -5.f, 5.f}}});
      rKzeroShort.add("hNSigmaNegPionK0s_before", "hNSigmaNegPionK0s_before", {HistType::kTH2F, {{ptAxis}, {100, -5.f, 5.f}}});
      rKzeroShort.add("hNSigmaNegPionK0s_after", "hNSigmaNegPionK0s_after", {HistType::kTH2F, {{ptAxis}, {100, -5.f, 5.f}}});
    }
    if (QAv0_daughters) {
      rKzeroShort.add("negative_pt", "Negative daughter pT", kTH1F, {ptAxis});
      rKzeroShort.add("positive_pt", "Positive daughter pT", kTH1F, {ptAxis});
      rKzeroShort.add("negative_eta", "Negative daughter eta", kTH1F, {{100, -1.0f, 1.0f}});
      rKzeroShort.add("positive_eta", "Positive daughter eta", kTH1F, {{100, -1.0f, 1.0f}});
      rKzeroShort.add("negative_phi", "Negative daughter phi", kTH1F, {{70, 0.0f, 7.0f}});
      rKzeroShort.add("positive_phi", "Positive daughter phi", kTH1F, {{70, 0.0f, 7.0f}});
    }
  }

  template <typename Collision>
  bool eventselection(Collision const& collision, const float& multiplicity)
  {
    if (!collision.sel8()) {
      return false;
    }
    if (timFrameEvsel && (!collision.selection_bit(aod::evsel::kNoTimeFrameBorder) || !collision.selection_bit(aod::evsel::kNoITSROFrameBorder))) {
      return false;
    }
    if (piluprejection && !collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      return false;
    }
    if (goodzvertex && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return false;
    }
    if (itstpctracks && !collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) {
      return false;
    }
    // if (collision.alias_bit(kTVXinTRD)) {
    //   // TRD triggered
    //   // return 0;
    // }
    auto multNTracksPV = collision.multNTracksPV();
    if (additionalEvsel && multNTracksPV < fMultPVCutLow->Eval(multiplicity)) {
      return false;
    }
    if (additionalEvsel && multNTracksPV > fMultPVCutHigh->Eval(multiplicity)) {
      return false;
    }
    // if (multTrk < fMultCutLow->Eval(multiplicity))
    //  return 0;
    // if (multTrk > fMultCutHigh->Eval(multiplicity))
    //  return 0;
    // if (multTrk > fMultMultPVCut->Eval(multNTracksPV))
    //  return 0;
    return true;
  }

  template <typename Collision, typename V0>
  bool SelectionV0(Collision const& collision, V0 const& candidate,
                   float /*multiplicity*/)
  {
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

    if (QAv0) {
      rKzeroShort.fill(HIST("hMassK0Shortbefore"), candidate.mK0Short(), candidate.pt());
      rKzeroShort.fill(HIST("hLT"), CtauK0s);
      rKzeroShort.fill(HIST("hDCAV0Daughters"), candidate.dcaV0daughters());
      rKzeroShort.fill(HIST("hV0CosPA"), candidate.v0cosPA());
      rKzeroShort.fill(HIST("Mass_lambda"), candidate.mLambda());
      rKzeroShort.fill(HIST("mass_AntiLambda"), candidate.mAntiLambda());
      rKzeroShort.fill(HIST("mass_Gamma"), candidate.mGamma());
      // rKzeroShort.fill(HIST("mass_Hypertriton"), candidate.mHypertriton());
      // rKzeroShort.fill(HIST("mass_AnitHypertriton"), candidate.mAntiHypertriton());
      rKzeroShort.fill(HIST("rapidity"), candidate.yK0Short());
      rKzeroShort.fill(HIST("hv0radius"), candidate.v0radius());
      rKzeroShort.fill(HIST("hDCApostopv"), candidate.dcapostopv());
      rKzeroShort.fill(HIST("hDCAnegtopv"), candidate.dcanegtopv());
      rKzeroShort.fill(HIST("hDCAv0topv"), candidate.dcav0topv());
      rKzeroShort.fill(HIST("halpha"), candidate.alpha());
      rKzeroShort.fill(HIST("hqtarm"), candidate.qtarm());
      rKzeroShort.fill(HIST("hpsipair"), candidate.psipair());
    }

    if (!DCAv0topv && fabs(candidate.dcav0topv()) > cMaxV0DCA) {
      return false;
    }

    if (rapidityks && TMath::Abs(candidate.yK0Short()) >= ConfKsrapidity) {
      return false;
    }

    if (masslambda && TMath::Abs(candidate.mLambda() - candidate.mK0Short()) >= competingcascrejlambda && TMath::Abs(candidate.mAntiLambda() - candidate.mK0Short()) >= competingcascrejlambdaanti) {
      return false;
    }

    // if (isStandarv0 && candidate.isStandardV0 == 0) {
    //   return false;
    // }

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
    if (fabs(CtauK0s) > cMaxV0LifeTime || candidate.mK0Short() < lowmasscutks0 || candidate.mK0Short() > highmasscutks0) {
      return false;
    }
    if (!armcut && arm < Confarmcut) {
      return false;
    }

    if (QAv0) {
      rKzeroShort.fill(HIST("hMassK0ShortSelected"), candidate.mK0Short(), candidate.pt());
    }
    return true;
  }

  template <typename T, typename V0s>
  bool isSelectedV0Daughter(T const& track, float charge,
                            double nsigmaV0Daughter, V0s const& /*candidate*/)
  {
    if (QAPID) {
      // Filling the PID of the V0 daughters in the region of the K0 peak.
      // tpcInnerParam is the momentum at the inner wall of TPC. So momentum of tpc vs nsigma of tpc is plotted.
      // if (0.45 < candidate.mK0Short() && candidate.mK0Short() < 0.55) {
      (charge == 1) ? rKzeroShort.fill(HIST("hNSigmaPosPionK0s_before"), track.tpcInnerParam(), track.tpcNSigmaPi()) : rKzeroShort.fill(HIST("hNSigmaNegPionK0s_before"), track.tpcInnerParam(), track.tpcNSigmaPi());
      // }
    }
    const auto eta = track.eta();
    const auto tpcNClsF = track.tpcNClsFound();
    // const auto dcaXY = track.dcaXY(); // for this we need TrackDCA table
    const auto sign = track.sign();

    if (!track.hasTPC())
      return false;
    if (track.tpcNClsCrossedRows() < tpcCrossedrows)
      return false;
    if (track.tpcCrossedRowsOverFindableCls() < tpcCrossedrowsOverfcls)
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
    // if (std::abs(dcaXY) < ConfDaughDCAMin) {
    //   return false;
    // }
    // v0 PID selection
    if (std::abs(nsigmaV0Daughter) > ConfDaughPIDCuts) {
      return false;
    }

    if (QAPID) {
      // if (0.45 < candidate.mK0Short() && candidate.mK0Short() < 0.55) {
      (charge == 1) ? rKzeroShort.fill(HIST("hNSigmaPosPionK0s_after"), track.tpcInnerParam(), track.tpcNSigmaPi()) : rKzeroShort.fill(HIST("hNSigmaNegPionK0s_after"), track.tpcInnerParam(), track.tpcNSigmaPi());
      // }
    }

    if (QAv0_daughters) {
      (charge == 1) ? rKzeroShort.fill(HIST("positive_pt"), track.pt()) : rKzeroShort.fill(HIST("negative_pt"), track.pt());
      (charge == 1) ? rKzeroShort.fill(HIST("positive_eta"), track.eta()) : rKzeroShort.fill(HIST("negative_eta"), track.eta());
      (charge == 1) ? rKzeroShort.fill(HIST("positive_phi"), track.phi()) : rKzeroShort.fill(HIST("negative_phi"), track.phi());
    }

    return true;
  }

  // Defining filters for events (event selection)
  // Filter eventFilter = (o2::aod::evsel::sel8 == true);
  Filter posZFilter = (nabs(o2::aod::collision::posZ) < cutzvertex);
  Filter AcceptenceFilter = (nabs(aod::track::eta) < cfgETAcut && nabs(aod::track::pt) > cfgPTcut);

  // Filters on V0s
  Filter preFilterV0 = (nabs(aod::v0data::dcapostopv) > v0setting_dcapostopv &&
                        nabs(aod::v0data::dcanegtopv) > v0setting_dcanegtopv);

  // Defining the type of the daughter tracks
  using DaughterTracks = soa::Join<aod::TracksIU, aod::TracksExtra, aod::pidTPCFullPi, aod::pidTOFFullPi>;
  using EventCandidates = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::FV0Mults, aod::MultZeqs, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::Mults>>;
  using TrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra, aod::TracksDCA, aod::pidTPCFullPi, aod::pidTOFFullPi>>;
  using V0TrackCandidate = aod::V0Datas;

  //   void processSE(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& collision,
  //                  soa::Filtered<aod::V0Datas> const& V0s,
  //                  DaughterTracks const&)

  void processSE(EventCandidates::iterator const& collision, TrackCandidates const& /*tracks*/, aod::V0Datas const& V0s)
  {
    hglue.fill(HIST("heventscheck"), 0.5);
    const double massK0s = TDatabasePDG::Instance()->GetParticle(kK0Short)->Mass();
    float multiplicity = 0.0f;
    if (cfgMultFOTM) {
      multiplicity = collision.centFT0M();
    } else {
      multiplicity = collision.centFT0C();
    }
    if (!eventselection(collision, multiplicity)) {
      return;
    }

    hglue.fill(HIST("heventscheck"), 1.5);

    rEventSelection.fill(HIST("hVertexZRec"), collision.posZ());
    rEventSelection.fill(HIST("hmultiplicity"), multiplicity);

    for (auto& [v1, v2] : combinations(CombinationsStrictlyUpperIndexPolicy(V0s, V0s))) {

      hglue.fill(HIST("heventscheck"), 2.5);

      if (v1.size() == 0 || v2.size() == 0) {
        continue;
      }

      if (!SelectionV0(collision, v1, multiplicity)) {
        continue;
      }
      if (!SelectionV0(collision, v2, multiplicity)) {
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

      if (!isSelectedV0Daughter(postrack1, 1, nTPCSigmaPos1[0], v1)) {
        continue;
      }
      if (!isSelectedV0Daughter(postrack2, 1, nTPCSigmaPos2[0], v2)) {
        continue;
      }
      if (!isSelectedV0Daughter(negtrack1, -1, nTPCSigmaNeg1[0], v1)) {
        continue;
      }
      if (!isSelectedV0Daughter(negtrack2, -1, nTPCSigmaNeg2[0], v2)) {
        continue;
      }

      hglue.fill(HIST("heventscheck"), 3.5);

      TLorentzVector lv1, lv2, lv3;
      lv1.SetPtEtaPhiM(v1.pt(), v1.eta(), v1.phi(), massK0s);
      lv2.SetPtEtaPhiM(v2.pt(), v2.eta(), v2.phi(), massK0s);
      lv3 = lv1 + lv2;

      if (TMath::Abs(lv3.Rapidity() < 0.5)) {
        if (inv_mass1D) {
          hglue.fill(HIST("h1glueInvMassDS"), lv3.M());
        }

        hglue.fill(HIST("h3glueInvMassDS"), multiplicity, lv3.Pt(), lv3.M());
      }
    }
  }

  PROCESS_SWITCH(strangeness_tutorial, processSE, "same event process", true);

  // use any one of 3 alias depending on the dataset. If pp then FT0M and if pbpb then FTOC
  using BinningTypeTPCMultiplicity = ColumnBinningPolicy<aod::collision::PosZ, aod::mult::MultTPC>;
  using BinningTypeCentralityM = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0M>;
  using BinningTypeVertexContributor = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0C>;
  ConfigurableAxis mevz = {"mevz", {10, -10., 10.}, "mixed event vertex z binning"};
  ConfigurableAxis memult = {"memult", {2000, 0, 10000}, "mixed event multiplicity binning"};

  void processME(EventCandidates const& collisions, TrackCandidates const& /*tracks*/, V0TrackCandidate const& v0s)
  {
    hglue.fill(HIST("heventscheck"), 4.5);

    const double massK0s = TDatabasePDG::Instance()->GetParticle(kK0Short)->Mass();
    auto tracksTuple = std::make_tuple(v0s);
    BinningTypeVertexContributor binningOnPositions1{{mevz, memult}, true};
    BinningTypeCentralityM binningOnPositions2{{mevz, memult}, true};

    SameKindPair<EventCandidates, V0TrackCandidate, BinningTypeVertexContributor> pair1{binningOnPositions1, cfgNmixedEvents, -1, collisions, tracksTuple, &cache}; // for PbPb
    SameKindPair<EventCandidates, V0TrackCandidate, BinningTypeCentralityM> pair2{binningOnPositions2, cfgNmixedEvents, -1, collisions, tracksTuple, &cache};       // for pp

    if (cfgMultFOTM) {
      for (auto& [c1, tracks1, c2, tracks2] : pair2) // two different centrality c1 and c2 and tracks corresponding to them
      {
        hglue.fill(HIST("heventscheck"), 5.5);

        float multiplicity = 0.0f;
        multiplicity = c1.centFT0M();

        if (!eventselection(c1, multiplicity) || !eventselection(c2, multiplicity)) {
          continue;
        }
        hglue.fill(HIST("heventscheck"), 6.5);

        for (auto& [t1, t2] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {
          hglue.fill(HIST("heventscheck"), 7.5);

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

          if (!isSelectedV0Daughter(postrack1, 1, nTPCSigmaPos1[0], t1)) {
            continue;
          }
          if (!isSelectedV0Daughter(postrack2, 1, nTPCSigmaPos2[0], t2)) {
            continue;
          }
          if (!isSelectedV0Daughter(negtrack1, -1, nTPCSigmaNeg1[0], t1)) {
            continue;
          }
          if (!isSelectedV0Daughter(negtrack2, -1, nTPCSigmaNeg2[0], t2)) {
            continue;
          }

          hglue.fill(HIST("heventscheck"), 8.5);

          TLorentzVector lv1, lv2, lv3;
          lv1.SetPtEtaPhiM(t1.pt(), t1.eta(), t1.phi(), massK0s);
          lv2.SetPtEtaPhiM(t2.pt(), t2.eta(), t2.phi(), massK0s);
          lv3 = lv1 + lv2;
          if (TMath::Abs(lv3.Rapidity() < 0.5)) {
            if (inv_mass1D) {
              hglue.fill(HIST("h1glueInvMassME"), lv3.M());
            }
            hglue.fill(HIST("h3glueInvMassME"), multiplicity, lv3.Pt(), lv3.M());
          }
        }
      }
    } else {
      for (auto& [c1, tracks1, c2, tracks2] : pair1) // two different centrality c1 and c2 and tracks corresponding to them
      {
        float multiplicity = 0.0f;
        multiplicity = c1.centFT0C();

        if (!eventselection(c1, multiplicity) || !eventselection(c2, multiplicity)) {
          continue;
        }

        for (auto& [t1, t2] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {
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

          if (!isSelectedV0Daughter(postrack1, 1, nTPCSigmaPos1[0], t1)) {
            continue;
          }
          if (!isSelectedV0Daughter(postrack2, 1, nTPCSigmaPos2[0], t2)) {
            continue;
          }
          if (!isSelectedV0Daughter(negtrack1, -1, nTPCSigmaNeg1[0], t1)) {
            continue;
          }
          if (!isSelectedV0Daughter(negtrack2, -1, nTPCSigmaNeg2[0], t2)) {
            continue;
          }

          TLorentzVector lv1, lv2, lv3;
          lv1.SetPtEtaPhiM(t1.pt(), t1.eta(), t1.phi(), massK0s);
          lv2.SetPtEtaPhiM(t2.pt(), t2.eta(), t2.phi(), massK0s);
          lv3 = lv1 + lv2;
          if (TMath::Abs(lv3.Rapidity() < 0.5)) {
            if (inv_mass1D) {
              hglue.fill(HIST("h1glueInvMassME"), lv3.M());
            }
            hglue.fill(HIST("h3glueInvMassME"), multiplicity, lv3.Pt(), lv3.M());
          }
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
