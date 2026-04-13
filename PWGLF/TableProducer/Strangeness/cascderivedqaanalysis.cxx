// Copyright 2019-2025 CERN and copyright holders of ALICE O2.
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
// Author Roman Nepeivoda & Romain Schotter

#include "PWGLF/DataModel/LFStrangenessPIDTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/cascqaanalysis.h"

#include "Framework/AnalysisTask.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"

#include <TPDGCode.h>
#include <string>
#include <algorithm>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using DauTracks = soa::Join<aod::DauTrackExtras, aod::DauTrackTPCPIDs>;
using CascCandidates = soa::Join<aod::CascCollRefs, aod::CascCores, aod::CascExtras, aod::CascTOFPIDs, aod::CascTOFNSigmas, aod::CascBBs>;
using CascMCCandidates = soa::Join<aod::CascCollRefs, aod::CascCores, aod::CascExtras, aod::CascTOFPIDs, aod::CascTOFNSigmas, aod::CascBBs, aod::CascCoreMCLabels>;

// constants
const float ctauXiPDG = 4.91;     // Xi PDG lifetime
const float ctauOmegaPDG = 2.461; // Omega PDG lifetime

struct CascDerivedQAanalysis {
  // Switch between Data/MC-dedicated histograms
  bool isMC;

  // for MC
  Configurable<bool> doMCAssociation{"doMCAssociation", true, "if MC, do MC association"};

  // Tables to produce
  Produces<aod::MyCascades> mycascades;
  Produces<aod::MyMCCascades> myMCcascades;

  HistogramRegistry histos{"histos"};

  // Event selection criteria
  struct : ConfigurableGroup {
    std::string prefix = "eventSelections"; // JSON group name
    Configurable<bool> requireSel8{"requireSel8", true, "require sel8 event selection"};
    Configurable<bool> requireTriggerTVX{"requireTriggerTVX", true, "require FT0 vertex (acceptable FT0C-FT0A time difference) at trigger level"};
    Configurable<bool> rejectITSROFBorder{"rejectITSROFBorder", true, "reject events at ITS ROF border (Run 3 only)"};
    Configurable<bool> rejectTFBorder{"rejectTFBorder", true, "reject events at TF border (Run 3 only)"};
    Configurable<bool> requireIsVertexITSTPC{"requireIsVertexITSTPC", false, "require events with at least one ITS-TPC track (Run 3 only)"};
    Configurable<bool> requireIsGoodZvtxFT0VsPV{"requireIsGoodZvtxFT0VsPV", false, "require events with PV position along z consistent (within 1 cm) between PV reconstructed using tracks and PV using FT0 A-C time difference (Run 3 only)"};
    Configurable<bool> requireIsVertexTOFmatched{"requireIsVertexTOFmatched", false, "require events with at least one of vertex contributors matched to TOF (Run 3 only)"};
    Configurable<bool> requireIsVertexTRDmatched{"requireIsVertexTRDmatched", false, "require events with at least one of vertex contributors matched to TRD (Run 3 only)"};
    Configurable<bool> rejectSameBunchPileup{"rejectSameBunchPileup", true, "reject collisions in case of pileup with another collision in the same foundBC (Run 3 only)"};
    Configurable<bool> requireNoCollInTimeRangeStd{"requireNoCollInTimeRangeStd", false, "reject collisions corrupted by the cannibalism, with other collisions within +/- 2 microseconds or mult above a certain threshold in -4 - -2 microseconds (Run 3 only)"};
    Configurable<bool> requireNoCollInTimeRangeStrict{"requireNoCollInTimeRangeStrict", false, "reject collisions corrupted by the cannibalism, with other collisions within +/- 10 microseconds (Run 3 only)"};
    Configurable<bool> requireNoCollInTimeRangeNarrow{"requireNoCollInTimeRangeNarrow", false, "reject collisions corrupted by the cannibalism, with other collisions within +/- 2 microseconds (Run 3 only)"};
    Configurable<bool> requireNoCollInROFStd{"requireNoCollInROFStd", false, "reject collisions corrupted by the cannibalism, with other collisions within the same ITS ROF with mult. above a certain threshold (Run 3 only)"};
    Configurable<bool> requireNoCollInROFStrict{"requireNoCollInROFStrict", false, "reject collisions corrupted by the cannibalism, with other collisions within the same ITS ROF (Run 3 only)"};
    Configurable<bool> requireINEL0{"requireINEL0", true, "require INEL>0 event selection"};
    Configurable<bool> requireINEL1{"requireINEL1", false, "require INEL>1 event selection"};

    Configurable<float> maxZVtxPosition{"maxZVtxPosition", 10., "max Z vtx position"};
  } eventSelections;

  struct : ConfigurableGroup {
    // Selection criteria: acceptance
    Configurable<float> rapidityCut{"cascSelections.rapidityCut", 0.5, "rapidity"};
    Configurable<float> daughterEtaCut{"cascSelections.daughterEtaCut", 0.8, "max eta for daughters"};

    // Standard 6 topological criteria on V0
    Configurable<float> v0cospa{"cascSelections.v0cospa", 0.97, "min V0 CosPA"};
    Configurable<float> dcav0dau{"cascSelections.dcav0dau", 1.0, "max DCA V0 Daughters (cm)"};
    Configurable<float> dcav0topv{"cascSelections.dcav0topv", .05, "min DCA V0 to PV (cm)"};
    Configurable<float> dcapiontopv{"cascSelections.dcapiontopv", .05, "min DCA Pion To PV (cm)"};
    Configurable<float> dcaprotontopv{"cascSelections.dcaprotontopv", .05, "min DCA Proton To PV (cm)"};
    Configurable<float> v0radius{"cascSelections.v0radius", 1.2, "minimum V0 radius (cm)"};
    Configurable<float> v0radiusMax{"cascSelections.v0radiusMax", 1E5, "maximum V0 radius (cm)"};

    // Standard 6 topological criteria on cascades
    Configurable<float> casccospa{"cascSelections.casccospa", 0.97, "min Cascade CosPA"};
    Configurable<float> dcacascdau{"cascSelections.dcacascdau", 1.0, "max DCA Cascade Daughters (cm)"};
    Configurable<float> dcaxybachbaryontopv{"cascSelections.dcaxybachbaryontopv", -1, "DCAxy Bachelor-Baryon to PV (cm)"};
    Configurable<float> bachbaryoncospa{"cascSelections.bachbaryoncospa", -1, "Bachelor-Baryon CosPA"};
    Configurable<float> dcabachtopv{"cascSelections.dcabachtopv", .05, "min DCA Bachelor To PV (cm)"};
    Configurable<float> cascradius{"cascSelections.cascradius", 0.5, "minimum Cascade radius (cm)"};
    Configurable<float> cascradiusMax{"cascSelections.cascradiusMax", 1E5, "maximum Cascade radius (cm)"};
    Configurable<float> cascProperLifeTime{"cascSelections.cascProperLifeTime", 3, "maximum lifetime (ctau)"};

    // invariant mass selection
    Configurable<float> v0MassWindow{"cascSelections.v0MassWindow", 0.008, "#Lambda mass (GeV/#it{c}^{2})"};
    Configurable<float> compMassRejection{"cascSelections.compMassRejection", 0.008, "Competing mass rejection (GeV/#it{c}^{2})"};

    // Track quality
    Configurable<int> minTPCrows{"cascSelections.minTPCrows", 70, "minimum TPC crossed rows"};
    Configurable<int> minITSclusters{"cascSelections.minITSclusters", -1, "minimum ITS clusters"};

    // PID (TPC/TOF)
    Configurable<float> tpcPidNsigmaCut{"cascSelections.tpcPidNsigmaCut", 5, "tpcPidNsigmaCut"};
    Configurable<float> tofPidNsigmaCutLaPr{"cascSelections.tofPidNsigmaCutLaPr", 1e+6, "tofPidNsigmaCutLaPr"};
    Configurable<float> tofPidNsigmaCutLaPi{"cascSelections.tofPidNsigmaCutLaPi", 1e+6, "tofPidNsigmaCutLaPi"};
    Configurable<float> tofPidNsigmaCutXiPi{"cascSelections.tofPidNsigmaCutXiPi", 1e+6, "tofPidNsigmaCutXiPi"};
    Configurable<float> tofPidNsigmaCutOmKa{"cascSelections.tofPidNsigmaCutOmKa", 1e+6, "tofPidNsigmaCutOmKa"};

    // PID (TOF)
    Configurable<float> maxDeltaTimeProton{"cascSelections.maxDeltaTimeProton", 1e+9, "check maximum allowed time"};
    Configurable<float> maxDeltaTimePion{"cascSelections.maxDeltaTimePion", 1e+9, "check maximum allowed time"};
    Configurable<float> maxDeltaTimeKaon{"cascSelections.maxDeltaTimeKaon", 1e+9, "check maximum allowed time"};
  } cascSelections;

  void init(InitContext const&)
  {
    // Determine if we are dealing with MC or not
    if (doprocessMCrec || doprocessMCgen) {
      isMC = true;
    } else {
      isMC = false;
    }

    TString hNEventsMCLabels[6] = {"All", "z vrtx", "INEL", "INEL>0", "INEL>1", "Associated with rec. collision"};
    TString hNEventsLabels[13] = {"All", "kIsTriggerTVX", "kNoTimeFrameBorder", "kNoITSROFrameBorder", "kIsVertexITSTPC", "kNoSameBunchPileup", "kIsGoodZvtxFT0vsPV", "isVertexTOFmatched", "kNoCollInTimeRangeNarrow", "z vrtx", "INEL", "INEL>0", "INEL>1"};

    histos.add("hEventSelection", "hEventSelection", kTH1D, {{21, -0.5f, +20.5f}});
    histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(1, "All collisions");
    histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(2, "sel8 cut");
    histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(3, "kIsTriggerTVX");
    histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(4, "kNoITSROFrameBorder");
    histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(5, "kNoTimeFrameBorder");
    histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(6, "posZ cut");
    histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(7, "kIsVertexITSTPC");
    histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(8, "kIsGoodZvtxFT0vsPV");
    histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(9, "kIsVertexTOFmatched");
    histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(10, "kIsVertexTRDmatched");
    histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(11, "kNoSameBunchPileup");
    histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(12, "kNoCollInTimeRangeStd");
    histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(13, "kNoCollInTimeRangeStrict");
    histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(14, "kNoCollInTimeRangeNarrow");
    histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(15, "kNoCollInRofStd");
    histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(16, "kNoCollInRofStrict");
    histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(17, "INEL>0");
    histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(18, "INEL>1");

    histos.add("hEventPVz", "hEventPVz", kTH1D, {{100, -20.0f, +20.0f}});

    if (isMC) {
      // Gen. lvl
      histos.add("hEventPVzMC", "hEventPVzMC", kTH1D, {{100, -20.0f, +20.0f}});
    }
  }

  // Manual grouping to avoid iterator subscription over StraMCCollisions
  std::vector<std::vector<int>> cascadesGrouped;

  template <typename TCascade, typename TCollision>
  bool isCascadeSelected(TCascade casc, TCollision collision, float rapidity, bool isXi)

  // precalculate this information so that a check is one mask operation, not many
  {
    //
    // Base topological variables
    //

    // v0 radius min/max selections
    if (casc.v0radius() < cascSelections.v0radius)
      return false;
    if (casc.v0radius() > cascSelections.v0radiusMax)
      return false;
    // DCA proton and pion to PV for Lambda and AntiLambda decay hypotheses
    if (casc.sign() < 0) { // Xi- or Omega- --> positive/negative daughter = proton/pion
      if (std::fabs(casc.dcapostopv()) < cascSelections.dcaprotontopv)
        return false;
      if (std::fabs(casc.dcanegtopv()) < cascSelections.dcapiontopv)
        return false;
    } else { // Xi+ or Omega+ --> positive/negative daughter = pion/proton
      if (std::fabs(casc.dcapostopv()) < cascSelections.dcapiontopv)
        return false;
      if (std::fabs(casc.dcanegtopv()) < cascSelections.dcaprotontopv)
        return false;
    }
    // V0 cosine of pointing angle
    if (casc.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) < cascSelections.v0cospa)
      return false;
    // DCA between v0 daughters
    if (casc.dcaV0daughters() > cascSelections.dcav0dau)
      return false;
    // DCA V0 to prim vtx
    if (casc.dcav0topv(collision.posX(), collision.posY(), collision.posZ()) < cascSelections.dcav0topv)
      return false;

    // casc radius min/max selections
    if (casc.cascradius() < cascSelections.cascradius)
      return false;
    if (casc.cascradius() > cascSelections.cascradiusMax)
      return false;
    // DCA bachelor selection
    if (std::fabs(casc.dcabachtopv()) < cascSelections.dcabachtopv)
      return false;
    // Bachelor-baryon cosPA selection
    if (casc.bachBaryonCosPA() < cascSelections.bachbaryoncospa)
      return false;
    // DCA bachelor-baryon selection
    if (std::fabs(casc.bachBaryonDCAxyToPV()) < cascSelections.dcaxybachbaryontopv)
      return false;
    // casc cosine of pointing angle
    if (casc.casccosPA(collision.posX(), collision.posY(), collision.posZ()) < cascSelections.casccospa)
      return false;
    // DCA between casc daughters
    if (casc.dcacascdaughters() > cascSelections.dcacascdau)
      return false;

    //
    // rapidity
    //
    if (std::fabs(rapidity) > cascSelections.rapidityCut)
      return false;

    //
    // invariant mass window
    //
    if (std::fabs(casc.mLambda() - o2::constants::physics::MassLambda0) > cascSelections.v0MassWindow)
      return false;

    //
    // competing mass rejection
    //
    if (isXi && std::fabs(casc.mOmega() - o2::constants::physics::MassOmegaMinus) < cascSelections.compMassRejection)
      return false;
    if (!isXi && std::fabs(casc.mXi() - o2::constants::physics::MassXiMinus) < cascSelections.compMassRejection)
      return false;

    auto bachTrackExtra = casc.template bachTrackExtra_as<DauTracks>();
    auto posTrackExtra = casc.template posTrackExtra_as<DauTracks>();
    auto negTrackExtra = casc.template negTrackExtra_as<DauTracks>();

    //
    // ITS quality flags
    //
    if (bachTrackExtra.itsNCls() < cascSelections.minITSclusters)
      return false;
    if (posTrackExtra.itsNCls() < cascSelections.minITSclusters)
      return false;
    if (negTrackExtra.itsNCls() < cascSelections.minITSclusters)
      return false;

    //
    // TPC quality flags
    //
    if (bachTrackExtra.tpcCrossedRows() < cascSelections.minTPCrows)
      return false;
    if (posTrackExtra.tpcCrossedRows() < cascSelections.minTPCrows)
      return false;
    if (negTrackExtra.tpcCrossedRows() < cascSelections.minTPCrows)
      return false;

    //
    // TPC PID
    //
    if (isXi && std::fabs(bachTrackExtra.tpcNSigmaPi()) > cascSelections.tpcPidNsigmaCut)
      return false;
    if (!isXi && std::fabs(bachTrackExtra.tpcNSigmaKa()) > cascSelections.tpcPidNsigmaCut)
      return false;
    if (casc.sign() < 0) { // Xi- or Omega- --> positive/negative daughter = proton/pion
      if (std::fabs(posTrackExtra.tpcNSigmaPr()) > cascSelections.tpcPidNsigmaCut)
        return false;
      if (std::fabs(negTrackExtra.tpcNSigmaPi()) > cascSelections.tpcPidNsigmaCut)
        return false;
    } else { // Xi+ or Omega+ --> positive/negative daughter = pion/proton
      if (std::fabs(posTrackExtra.tpcNSigmaPi()) > cascSelections.tpcPidNsigmaCut)
        return false;
      if (std::fabs(negTrackExtra.tpcNSigmaPr()) > cascSelections.tpcPidNsigmaCut)
        return false;
    }

    //
    // TOF PID in DeltaT
    // Bachelor track
    if (bachTrackExtra.hasTOF()) {
      if (isXi && std::fabs(casc.bachTOFDeltaTXiPi()) > cascSelections.maxDeltaTimePion)
        return false;
      if (!isXi && std::fabs(casc.bachTOFDeltaTOmKa()) > cascSelections.maxDeltaTimeKaon)
        return false;
    }
    // Positive track
    if (posTrackExtra.hasTOF()) {
      if (casc.sign() < 0) { // Xi- or Omega- --> positive daughter = proton
        if (isXi && std::fabs(casc.posTOFDeltaTXiPr()) > cascSelections.maxDeltaTimeProton)
          return false;
        if (!isXi && std::fabs(casc.posTOFDeltaTOmPr()) > cascSelections.maxDeltaTimeProton)
          return false;
      } else { // Xi+ or Omega+ --> positive daughter = pion
        if (isXi && std::fabs(casc.posTOFDeltaTXiPi()) > cascSelections.maxDeltaTimePion)
          return false;
        if (!isXi && std::fabs(casc.posTOFDeltaTOmPi()) > cascSelections.maxDeltaTimePion)
          return false;
      }
    }
    // Negative track
    if (negTrackExtra.hasTOF()) {
      if (casc.sign() < 0) { // Xi- or Omega- --> negative daughter = pion
        if (isXi && std::fabs(casc.negTOFDeltaTXiPi()) > cascSelections.maxDeltaTimePion)
          return false;
        if (!isXi && std::fabs(casc.negTOFDeltaTOmPi()) > cascSelections.maxDeltaTimePion)
          return false;
      } else { // Xi+ or Omega+ --> negative daughter = proton
        if (isXi && std::fabs(casc.negTOFDeltaTXiPr()) > cascSelections.maxDeltaTimeProton)
          return false;
        if (!isXi && std::fabs(casc.negTOFDeltaTOmPr()) > cascSelections.maxDeltaTimeProton)
          return false;
      }
    }

    //
    // TOF PID in NSigma
    // Bachelor track
    if (bachTrackExtra.hasTOF()) {
      if (isXi && std::fabs(casc.tofNSigmaXiPi()) > cascSelections.tofPidNsigmaCutXiPi)
        return false;
      if (!isXi && std::fabs(casc.tofNSigmaOmKa()) > cascSelections.tofPidNsigmaCutOmKa)
        return false;
    }
    // Positive track
    if (posTrackExtra.hasTOF()) {
      if (casc.sign() < 0) { // Xi- or Omega- --> positive daughter = proton
        if (isXi && std::fabs(casc.tofNSigmaXiLaPr()) > cascSelections.tofPidNsigmaCutLaPr)
          return false;
        if (!isXi && std::fabs(casc.tofNSigmaOmLaPr()) > cascSelections.tofPidNsigmaCutLaPr)
          return false;
      } else { // Xi+ or Omega+ --> positive daughter = pion
        if (isXi && std::fabs(casc.tofNSigmaXiLaPi()) > cascSelections.tofPidNsigmaCutLaPi)
          return false;
        if (!isXi && std::fabs(casc.tofNSigmaOmLaPi()) > cascSelections.tofPidNsigmaCutLaPi)
          return false;
      }
    }
    // Negative track
    if (negTrackExtra.hasTOF()) {
      if (casc.sign() < 0) { // Xi- or Omega- --> negative daughter = pion
        if (isXi && std::fabs(casc.tofNSigmaXiLaPr()) > cascSelections.tofPidNsigmaCutLaPi)
          return false;
        if (!isXi && std::fabs(casc.tofNSigmaOmLaPr()) > cascSelections.tofPidNsigmaCutLaPi)
          return false;
      } else { // Xi+ or Omega+ --> negative daughter = proton
        if (isXi && std::fabs(casc.tofNSigmaXiLaPi()) > cascSelections.tofPidNsigmaCutLaPr)
          return false;
        if (!isXi && std::fabs(casc.tofNSigmaOmLaPi()) > cascSelections.tofPidNsigmaCutLaPr)
          return false;
      }
    }

    //
    // proper lifetime
    float distOverTotMom = std::sqrt(std::pow(casc.x() - collision.posX(), 2) + std::pow(casc.y() - collision.posY(), 2) + std::pow(casc.z() - collision.posZ(), 2)) / (casc.p() + 1E-10);
    if (isXi && distOverTotMom * o2::constants::physics::MassXiMinus / ctauXiPDG > cascSelections.cascProperLifeTime)
      return false;
    if (!isXi && distOverTotMom * o2::constants::physics::MassOmegaMinus / ctauOmegaPDG > cascSelections.cascProperLifeTime)
      return false;

    //
    // MC association (if asked)
    if (doMCAssociation) {
      if constexpr (requires { casc.template cascMCCore_as<soa::Join<aod::CascMCCores, aod::CascMCCollRefs>>(); }) { // check if MC information is available
        auto cascMC = casc.template cascMCCore_as<soa::Join<aod::CascMCCores, aod::CascMCCollRefs>>();

        if (isXi) {
          if (casc.sign() < 0) {
            if (cascMC.pdgCode() != 3312 || cascMC.pdgCodePositive() != 2212 || cascMC.pdgCodeNegative() != -211 || cascMC.pdgCodeBachelor() != -211)
              return false;
          } else {
            if (cascMC.pdgCode() != -3312 || cascMC.pdgCodePositive() != 211 || cascMC.pdgCodeNegative() != -2212 || cascMC.pdgCodeBachelor() != 211)
              return false;
          }
        } else {
          if (casc.sign() < 0) {
            if (cascMC.pdgCode() != 3334 || cascMC.pdgCodePositive() != 2212 || cascMC.pdgCodeNegative() != -211 || cascMC.pdgCodeBachelor() != -321)
              return false;
          } else {
            if (cascMC.pdgCode() != -3334 || cascMC.pdgCodePositive() != 211 || cascMC.pdgCodeNegative() != -2212 || cascMC.pdgCodeBachelor() != 321)
              return false;
          }
        }
      }
    }

    return true;
  }

  template <typename TCollision>
  bool isEventAccepted(TCollision collision, bool fillHists)
  // check whether the collision passes our collision selections
  {
    if (fillHists) {
      histos.fill(HIST("hEventSelection"), 0. /* all collisions */);
    }

    if (eventSelections.requireSel8 && !collision.sel8()) {
      return false;
    }
    if (fillHists) {
      histos.fill(HIST("hEventSelection"), 1 /* sel8 collisions */);
    }

    if (eventSelections.requireTriggerTVX && !collision.selection_bit(aod::evsel::kIsTriggerTVX)) {
      return false;
    }
    if (fillHists) {
      histos.fill(HIST("hEventSelection"), 2 /* FT0 vertex (acceptable FT0C-FT0A time difference) collisions */);
    }

    if (eventSelections.rejectITSROFBorder && !collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
      return false;
    }
    if (fillHists) {
      histos.fill(HIST("hEventSelection"), 3 /* Not at ITS ROF border */);
    }

    if (eventSelections.rejectTFBorder && !collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
      return false;
    }
    if (fillHists) {
      histos.fill(HIST("hEventSelection"), 4 /* Not at TF border */);
    }

    if (std::abs(collision.posZ()) > eventSelections.maxZVtxPosition) {
      return false;
    }
    if (fillHists) {
      histos.fill(HIST("hEventSelection"), 5 /* vertex-Z selected */);
    }

    if (eventSelections.requireIsVertexITSTPC && !collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) {
      return false;
    }
    if (fillHists) {
      histos.fill(HIST("hEventSelection"), 6 /* Contains at least one ITS-TPC track */);
    }

    if (eventSelections.requireIsGoodZvtxFT0VsPV && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return false;
    }
    if (fillHists) {
      histos.fill(HIST("hEventSelection"), 7 /* PV position consistency check */);
    }

    if (eventSelections.requireIsVertexTOFmatched && !collision.selection_bit(o2::aod::evsel::kIsVertexTOFmatched)) {
      return false;
    }
    if (fillHists) {
      histos.fill(HIST("hEventSelection"), 8 /* PV with at least one contributor matched with TOF */);
    }

    if (eventSelections.requireIsVertexTRDmatched && !collision.selection_bit(o2::aod::evsel::kIsVertexTRDmatched)) {
      return false;
    }
    if (fillHists) {
      histos.fill(HIST("hEventSelection"), 9 /* PV with at least one contributor matched with TRD */);
    }

    if (eventSelections.rejectSameBunchPileup && !collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      return false;
    }
    if (fillHists) {
      histos.fill(HIST("hEventSelection"), 10 /* Not at same bunch pile-up */);
    }

    if (eventSelections.requireNoCollInTimeRangeStd && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
      return false;
    }
    if (fillHists) {
      histos.fill(HIST("hEventSelection"), 11 /* No other collision within +/- 2 microseconds or mult above a certain threshold in -4 - -2 microseconds*/);
    }

    if (eventSelections.requireNoCollInTimeRangeStrict && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStrict)) {
      return false;
    }
    if (fillHists) {
      histos.fill(HIST("hEventSelection"), 12 /* No other collision within +/- 10 microseconds */);
    }

    if (eventSelections.requireNoCollInTimeRangeNarrow && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeNarrow)) {
      return false;
    }
    if (fillHists) {
      histos.fill(HIST("hEventSelection"), 13 /* No other collision within +/- 2 microseconds */);
    }

    if (eventSelections.requireNoCollInROFStd && !collision.selection_bit(o2::aod::evsel::kNoCollInRofStandard)) {
      return false;
    }
    if (fillHists) {
      histos.fill(HIST("hEventSelection"), 14 /* No other collision within the same ITS ROF with mult. above a certain threshold */);
    }

    if (eventSelections.requireNoCollInROFStrict && !collision.selection_bit(o2::aod::evsel::kNoCollInRofStrict)) {
      return false;
    }
    if (fillHists) {
      histos.fill(HIST("hEventSelection"), 15 /* No other collision within the same ITS ROF */);
    }

    if (eventSelections.requireINEL0 && collision.multNTracksPVeta1() < 1) {
      return false;
    }
    if (fillHists) {
      histos.fill(HIST("hEventSelection"), 16 /* INEL > 0 */);
    }

    if (eventSelections.requireINEL1 && collision.multNTracksPVeta1() < 2) {
      return false;
    }
    if (fillHists) {
      histos.fill(HIST("hEventSelection"), 17 /* INEL > 1 */);
    }

    return true;
  }

  void processData(soa::Join<aod::StraCollisions, aod::StraEvSels, aod::StraCents>::iterator const& collision,
                   CascCandidates const& Cascades,
                   DauTracks const&)
  {
    if (!isEventAccepted(collision, true)) {
      return;
    }
    histos.fill(HIST("hEventPVz"), collision.posZ());

    for (const auto& cascade : Cascades) { // loop over Cascades
      if (std::abs(cascade.negativeeta()) > cascSelections.daughterEtaCut ||
          std::abs(cascade.positiveeta()) > cascSelections.daughterEtaCut ||
          std::abs(cascade.bacheloreta()) > cascSelections.daughterEtaCut)
        continue; // remove acceptance that's badly reproduced by MC / superfluous in future

      if (isCascadeSelected(cascade, collision, cascade.yXi(), true) ||
          isCascadeSelected(cascade, collision, cascade.yOmega(), false)) {

        // Fill table
        auto bachTrackExtra = cascade.template bachTrackExtra_as<DauTracks>();
        auto posTrackExtra = cascade.template posTrackExtra_as<DauTracks>();
        auto negTrackExtra = cascade.template negTrackExtra_as<DauTracks>();

        int bachITSNhits = bachTrackExtra.itsNCls();
        int posITSNhits = posTrackExtra.itsNCls();
        int negITSNhits = negTrackExtra.itsNCls();

        uint8_t evFlag = 0;
        evFlag |= o2::aod::mycascades::EvFlags::EvINEL;
        if (collision.multNTracksPVeta1() > 0) {
          evFlag |= o2::aod::mycascades::EvFlags::EvINELgt0;
        }
        if (collision.multNTracksPVeta1() > 1) {
          evFlag |= o2::aod::mycascades::EvFlags::EvINELgt1;
        }

        // c x tau
        float distOverTotMom = std::sqrt(std::pow(cascade.x() - collision.posX(), 2) + std::pow(cascade.y() - collision.posY(), 2) + std::pow(cascade.z() - collision.posZ(), 2)) / (cascade.p() + 1E-10);
        float ctauXi = distOverTotMom * o2::constants::physics::MassXiMinus;
        float ctauOmega = distOverTotMom * o2::constants::physics::MassOmegaMinus;

        mycascades(collision.posZ(),
                   collision.centFT0M(), collision.centFV0A(),
                   collision.multFT0A() + collision.multFT0C(), collision.multFV0A(),
                   cascade.sign(), cascade.pt(), cascade.yXi(), cascade.yOmega(), cascade.eta(),
                   cascade.mXi(), cascade.mOmega(), cascade.mLambda(), cascade.cascradius(), cascade.v0radius(),
                   cascade.casccosPA(collision.posX(), collision.posY(), collision.posZ()), cascade.v0cosPA(collision.posX(), collision.posY(), collision.posZ()),
                   cascade.dcapostopv(), cascade.dcanegtopv(), cascade.dcabachtopv(), cascade.dcacascdaughters(), cascade.dcaV0daughters(), cascade.dcav0topv(collision.posX(), collision.posY(), collision.posZ()),
                   cascade.positiveeta(), cascade.negativeeta(), cascade.bacheloreta(), posITSNhits, negITSNhits, bachITSNhits,
                   ctauXi, ctauOmega, negTrackExtra.tpcNSigmaPr(), posTrackExtra.tpcNSigmaPr(), negTrackExtra.tpcNSigmaPi(), posTrackExtra.tpcNSigmaPi(), bachTrackExtra.tpcNSigmaPi(), bachTrackExtra.tpcNSigmaKa(),
                   -999.f /*negdau.tofNSigmaPr()*/, -999.f /*posdau.tofNSigmaPr()*/, -999.f /*negdau.tofNSigmaPi()*/, -999.f /*posdau.tofNSigmaPi()*/, cascade.bachTOFDeltaTXiPi(), cascade.bachTOFDeltaTOmKa(),
                   posTrackExtra.tpcNClsFound(), negTrackExtra.tpcNClsFound(), bachTrackExtra.tpcNClsFound(),
                   posTrackExtra.tpcCrossedRows(), negTrackExtra.tpcCrossedRows(), bachTrackExtra.tpcCrossedRows(),
                   posTrackExtra.hasTOF(), negTrackExtra.hasTOF(), bachTrackExtra.hasTOF(),
                   cascade.positivept(), cascade.negativept(), cascade.bachelorpt(), -1, -1, cascade.bachBaryonCosPA(), cascade.bachBaryonDCAxyToPV(), evFlag, 1e3, 1e3);
      }
    }
  }

  Preslice<aod::McParticles> perMcCollision = aod::mcparticle::mcCollisionId;

  // SONRA:
  void processMCrec(soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels, aod::StraCollLabels> const& collisions,
                    CascMCCandidates const& Cascades,
                    DauTracks const&,
                    aod::MotherMCParts const&,
                    soa::Join<aod::StraMCCollisions, aod::StraMCCollMults> const& /*mccollisions*/,
                    soa::Join<aod::CascMCCores, aod::CascMCCollRefs> const& cascMCCores)

  // SONRA:
  {
    // Manual grouping to avoid iterator subscription which would require
    // CascMCCollRefs to have a sorted fIndexStraMCCollisions index
    cascadesGrouped.clear();
    cascadesGrouped.resize(collisions.size());
    for (const auto& cascade : Cascades) {
      cascadesGrouped[cascade.straCollisionId()].push_back(cascade.globalIndex());
    }

    for (const auto& collision : collisions) {
      if (!isEventAccepted(collision, true)) {
        continue;
      }
      histos.fill(HIST("hEventPVz"), collision.posZ());

      if (!collision.has_straMCCollision()) {
        continue;
      }

      const auto& mcCollision = collision.straMCCollision_as<soa::Join<aod::StraMCCollisions, aod::StraMCCollMults>>(); // aod::McCentFV0As to be added
      histos.fill(HIST("hEventPVzMC"), mcCollision.posZ());

      // SONRA:
      for (auto const& idx : cascadesGrouped[collision.globalIndex()]) {
        auto cascade = Cascades.rawIteratorAt(idx);

        if (std::abs(cascade.negativeeta()) > cascSelections.daughterEtaCut ||
            std::abs(cascade.positiveeta()) > cascSelections.daughterEtaCut ||
            std::abs(cascade.bacheloreta()) > cascSelections.daughterEtaCut)
          continue; // remove acceptance that's badly reproduced by MC / superfluous in future

        if (!cascade.has_cascMCCore())
          continue;

        // SONRA:
        auto cascadeMC = cascade.cascMCCore_as<soa::Join<aod::CascMCCores, aod::CascMCCollRefs>>();

        // Check mc association
        float lPDG = 1e3;
        float ptmc = 1e3;
        float ymc = 1e3;
        float isPrimary = -1;
        if (std::abs(cascadeMC.pdgCode()) == PDG_t::kXiMinus || std::abs(cascadeMC.pdgCode()) == PDG_t::kOmegaMinus) {
          lPDG = cascadeMC.pdgCode();
          ptmc = RecoDecay::pt(cascadeMC.pxMC(), cascadeMC.pyMC());
          if (std::fabs(cascadeMC.pdgCode()) == 3312)
            ymc = RecoDecay::y(std::array{cascadeMC.pxMC(), cascadeMC.pyMC(), cascadeMC.pzMC()}, o2::constants::physics::MassXiMinus);
          else if (std::fabs(cascadeMC.pdgCode()) == 3334)
            ymc = RecoDecay::y(std::array{cascadeMC.pxMC(), cascadeMC.pyMC(), cascadeMC.pzMC()}, o2::constants::physics::MassOmegaMinus);
          isPrimary = cascadeMC.isPhysicalPrimary() ? 1 : 0;
        }

        if (isCascadeSelected(cascade, collision, ymc, true) ||
            isCascadeSelected(cascade, collision, ymc, false)) {

          auto bachTrackExtra = cascade.template bachTrackExtra_as<DauTracks>();
          auto posTrackExtra = cascade.template posTrackExtra_as<DauTracks>();
          auto negTrackExtra = cascade.template negTrackExtra_as<DauTracks>();

          int bachITSNhits = bachTrackExtra.itsNCls();
          int posITSNhits = posTrackExtra.itsNCls();
          int negITSNhits = negTrackExtra.itsNCls();

          uint8_t evFlag = 0;
          evFlag |= o2::aod::mycascades::EvFlags::EvINEL;
          if (collision.multNTracksPVeta1() > 0) {
            evFlag |= o2::aod::mycascades::EvFlags::EvINELgt0;
          }
          if (collision.multNTracksPVeta1() > 1) {
            evFlag |= o2::aod::mycascades::EvFlags::EvINELgt1;
          }

          // c x tau
          float distOverTotMom = std::sqrt(std::pow(cascade.x() - collision.posX(), 2) + std::pow(cascade.y() - collision.posY(), 2) + std::pow(cascade.z() - collision.posZ(), 2)) / (cascade.p() + 1E-10);
          float ctauXi = distOverTotMom * o2::constants::physics::MassXiMinus;
          float ctauOmega = distOverTotMom * o2::constants::physics::MassOmegaMinus;

          mycascades(collision.posZ(),
                     collision.centFT0M(), collision.centFV0A(),
                     collision.multFT0A() + collision.multFT0C(), collision.multFV0A(),
                     cascade.sign(), cascade.pt(), cascade.yXi(), cascade.yOmega(), cascade.eta(),
                     cascade.mXi(), cascade.mOmega(), cascade.mLambda(), cascade.cascradius(), cascade.v0radius(),
                     cascade.casccosPA(collision.posX(), collision.posY(), collision.posZ()), cascade.v0cosPA(collision.posX(), collision.posY(), collision.posZ()),
                     cascade.dcapostopv(), cascade.dcanegtopv(), cascade.dcabachtopv(), cascade.dcacascdaughters(), cascade.dcaV0daughters(), cascade.dcav0topv(collision.posX(), collision.posY(), collision.posZ()),
                     cascade.positiveeta(), cascade.negativeeta(), cascade.bacheloreta(), posITSNhits, negITSNhits, bachITSNhits,
                     ctauXi, ctauOmega, negTrackExtra.tpcNSigmaPr(), posTrackExtra.tpcNSigmaPr(), negTrackExtra.tpcNSigmaPi(), posTrackExtra.tpcNSigmaPi(), bachTrackExtra.tpcNSigmaPi(), bachTrackExtra.tpcNSigmaKa(),
                     -999.f /*negdau.tofNSigmaPr()*/, -999.f /*posdau.tofNSigmaPr()*/, -999.f /*negdau.tofNSigmaPi()*/, -999.f /*posdau.tofNSigmaPi()*/, cascade.bachTOFDeltaTXiPi(), cascade.bachTOFDeltaTOmKa(),
                     posTrackExtra.tpcNClsFound(), negTrackExtra.tpcNClsFound(), bachTrackExtra.tpcNClsFound(),
                     posTrackExtra.tpcCrossedRows(), negTrackExtra.tpcCrossedRows(), bachTrackExtra.tpcCrossedRows(),
                     posTrackExtra.hasTOF(), negTrackExtra.hasTOF(), bachTrackExtra.hasTOF(),
                     cascadeMC.positiveptMC(), cascadeMC.negativeptMC(), cascadeMC.bachelorptMC(), lPDG, isPrimary, cascade.bachBaryonCosPA(), cascade.bachBaryonDCAxyToPV(), evFlag, ptmc, ymc);

        } // isCascadeSelected kapanışı  ← girinti 8 boşluk
      } // end cascade loop            ← girinti 6 boşluk
    } // end collision loop          ← girinti 4 boşluk
  } // processMCrec kapanışı       ← girinti 2 boşluk

  void processMCgen(soa::Join<aod::StraMCCollisions, aod::StraMCCollMults> const& mcCollisions,
                    soa::Join<aod::CascMCCores, aod::CascMCCollRefs> const& CascMCCores,
                    soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels, aod::StraStamps, aod::StraCollLabels> const& collisions)
  {

    // Manual grouping of CascMCCores per MC collision
    cascadesGrouped.clear();
    cascadesGrouped.resize(mcCollisions.size());
    for (auto const& cascMC : CascMCCores) {
      if (!cascMC.has_straMCCollision())
        continue;
      cascadesGrouped[cascMC.straMCCollisionId()].push_back(cascMC.globalIndex());
    }

    for (const auto& mcCollision : mcCollisions) {

      // Define the type of generated MC collision
      histos.fill(HIST("hEventPVzMC"), mcCollision.posZ());

      int evType = 0;
      uint8_t flagsGen = 0;
      flagsGen |= o2::aod::myMCcascades::EvFlags::EvINEL;
      // Generated collision is INEL>0
      if (mcCollision.multMCNParticlesEta10() > 0) {
        flagsGen |= o2::aod::myMCcascades::EvFlags::EvINELgt0;
        evType++;
      }
      // Generated collision is INEL>1
      if (mcCollision.multMCNParticlesEta10() > 1) {
        flagsGen |= o2::aod::myMCcascades::EvFlags::EvINELgt1;
        evType++;
      }

      uint16_t nchFT0 = mcCollision.multMCFT0A() + mcCollision.multMCFT0C();
      uint8_t flagsAssoc = 0;
      int nAssocColl = 0;
      int biggestNContribs = -1;
      float centrality = 100.5f;
      for (const auto& collision : collisions) {
        if (!isEventAccepted(collision, false)) {
          continue;
        }

        if (biggestNContribs < collision.multPVTotalContributors()) {
          biggestNContribs = collision.multPVTotalContributors();
          centrality = collision.centFT0M();
        }

        if (collision.straMCCollision_as<soa::Join<aod::StraMCCollisions, aod::StraMCCollMults>>().globalIndex() == mcCollision.globalIndex()) {
          flagsAssoc |= o2::aod::myMCcascades::EvFlags::EvINEL;
          if (collision.multNTracksPVeta1() > 0) {
            flagsAssoc |= o2::aod::myMCcascades::EvFlags::EvINELgt0;
          }
          if (collision.multNTracksPVeta1() > 1) {
            flagsAssoc |= o2::aod::myMCcascades::EvFlags::EvINELgt1;
          }
          nAssocColl++;
        }
      }

      for (auto const& idx : cascadesGrouped[mcCollision.globalIndex()]) {
        auto cascMC = CascMCCores.rawIteratorAt(idx);

        if (!cascMC.isPhysicalPrimary())
          continue;
        float sign = 0;
        if (cascMC.pdgCode() == PDG_t::kXiPlusBar || cascMC.pdgCode() == PDG_t::kOmegaPlusBar) {
          sign = 1;
        } else if (cascMC.pdgCode() == PDG_t::kXiMinus || cascMC.pdgCode() == PDG_t::kOmegaMinus) {
          sign = -1;
        } else {
          continue;
        }
        myMCcascades(mcCollision.posZ(), sign, cascMC.pdgCode(),
                     std::abs(cascMC.pdgCode()) == PDG_t::kXiMinus ? cascMC.rapidityMC(0) : cascMC.rapidityMC(1), RecoDecay::eta(std::array{cascMC.pxMC(), cascMC.pyMC(), cascMC.pzMC()}), RecoDecay::phi(cascMC.pxMC(), cascMC.pyMC()), cascMC.ptMC(),
                     cascMC.isPhysicalPrimary(), nAssocColl,
                     nchFT0, 0,
                     centrality, 0, // mcCollision.centFV0A() to be added
                     flagsAssoc,
                     flagsGen);

      } // end cascade loop
    } // end mcCollision loop
  }

  PROCESS_SWITCH(CascDerivedQAanalysis, processData, "Process Run 3 data", true);
  PROCESS_SWITCH(CascDerivedQAanalysis, processMCrec, "Process Run 3 mc, reconstructed", false);
  PROCESS_SWITCH(CascDerivedQAanalysis, processMCgen, "Process Run 3 mc, genereated", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<CascDerivedQAanalysis>(cfgc)};
}
