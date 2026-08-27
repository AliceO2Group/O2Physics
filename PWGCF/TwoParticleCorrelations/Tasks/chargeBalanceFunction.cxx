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

/// \file chargeBalanceFunction.cxx
/// \brief R2 P2 and BF of charged hadrons.
/// \author Yash Patley <yash.patley@cern.ch>

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/MathConstants.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/SliceCache.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>
#include <TH3.h>
#include <THnSparse.h>
#include <TList.h>
#include <TObject.h>
#include <TPDGCode.h>

#include <chrono>
#include <cmath>
#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;
using namespace o2::constants::math;

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;
using namespace o2::constants::math;

enum ParticlePairType {
  kPM = 0,
  kPP,
  kMM
};

enum RecGenType {
  kRec = 0,
  kGen
};

enum DMCType {
  kData = 0,
  kMC
};

struct ChargeBalanceFunction {
  // Collisions
  Configurable<float> cZVtxCut{"cZVtxCut", 10.0, "Min VtxZ cut"};
  Configurable<float> cMinCent{"cMinCent", 0., "Minumum Centrality"};
  Configurable<float> cMaxCent{"cMaxCent", 100.0, "Maximum Centrality"};
  Configurable<bool> cSel8Trig{"cSel8Trig", true, "Sel8 (T0A + T0C) Selection Run3"};
  Configurable<bool> cPileupReject{"cPileupReject", true, "Pileup rejection"};
  Configurable<bool> cZVtxTimeDiff{"cZVtxTimeDiff", true, "z-vtx time diff selection"};
  Configurable<bool> cIsGoodITSLayers{"cIsGoodITSLayers", true, "Good ITS Layers All"};
  Configurable<float> cMinOccupancy{"cMinOccupancy", 0, "Minimum FT0C Occupancy"};
  Configurable<float> cMaxOccupancy{"cMaxOccupancy", 1e6, "Maximum FT0C Occupancy"};

  // Tracks
  Configurable<int> cTrackNPtBins{"cTrackNPtBins", 20, "N pT bins"};
  Configurable<float> cTrackMinPt{"cTrackMinPt", 0.2, "p_{T} minimum"};
  Configurable<float> cTrackMaxPt{"cTrackMaxPt", 2.0, "p_{T} maximum"};
  Configurable<float> cTrackEtaCut{"cTrackEtaCut", 0.8, "Pseudorapidity cut"};
  Configurable<bool> cTrackGlobal{"cTrackGlobal", true, "Global Track"};
  Configurable<float> cTrackDcaXYCut{"cTrackDcaXYCut", 0.1, "DcaXY Cut"};
  Configurable<float> cTrackDcaZCut{"cTrackDcaZCut", 1., "DcaZ Cut"};
  Configurable<float> cTpcElRejCutMin{"cTpcElRejCutMin", -3., "Electron Rejection Cut Minimum"};
  Configurable<float> cTpcElRejCutMax{"cTpcElRejCutMax", 5., "Electron Rejection Cut Maximum"};
  Configurable<float> cTpcRejCut{"cTpcRejCut", 3, "TPC Rej Cut"};

  // Correlation
  Configurable<int> cNEtaBins{"cNEtaBins", 16, "N Eta Bins"};
  Configurable<int> cNPhiBins{"cNPhiBins", 36, "N Phi Bins"};

  // Femtoscopic correction
  Configurable<bool> cApplyFemtoSel{"cApplyFemtoSel", false, "Femto qinv selection"};
  Configurable<float> cFemtoCut{"cFemtoCut", 0.1, "Kaon--Lambda Femto qinv cut"};

  // Two-track cuts
  Configurable<bool> cApplyTwoTrackCut{"cApplyTwoTrackCut", false, "Flag for two track cut"};
  Configurable<float> cDEtaCut{"cDEtaCut", 0.02, "DEta cut"};
  Configurable<float> cDPhiStarCut{"cDPhiStarCut", 0.02, "DPhiStar cut"};

  // Efficiency Correction
  Configurable<bool> cGetCorrectionFlag{"cGetCorrectionFlag", false, "Apply correction flag"};
  Configurable<bool> cGetNuaCorrectionFlag{"cGetNuaCorrectionFlag", false, "Apply NUA correction flag"};

  // CCDB
  Configurable<std::string> cUrlCCDB{"cUrlCCDB", "http://alice-ccdb.cern.ch", "ALICE CCDB URL"};
  Configurable<std::string> cPathCCDBRecoEff{"cPathCCDBRecoEff", "Users/y/ypatley/CBF/Test/RecoEfficiency", "Path for ccdb-object for reco efficiency"};
  Configurable<std::string> cPathCCDBNuaCorr{"cPathCCDBNuaCorr", "Users/y/ypatley/CBF/Test/Nua", "Path for ccdb-object for NUA correction"};
  Configurable<int64_t> nolaterthan{"nolaterthan", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};

  // Configurable Axis
  ConfigurableAxis cCentBins{"cCentBins", {VARIABLE_WIDTH, 0.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.f, 60.0f, 70.0f, 80.0f, 90.0f, 100.f}, "Variable Centrality Bins"};
  ConfigurableAxis cPosZBins{"cPosZBins", {VARIABLE_WIDTH, -10.f, -8.0f, -6.0f, -4.0f, -2.0f, 0.f, 2.0f, 4.0f, 6.0f, 8.0f, 10.f}, "Variable Vz Bins"};

  // Ccdb service
  Service<o2::ccdb::BasicCCDBManager> ccdb{};
  o2::parameters::GRPMagField* grpo = nullptr;

  // Pdg service
  Service<o2::framework::O2DatabasePDG> pdg{};

  // Histogram Registry.
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // Initialize global variables
  float nrapbins = 0.;
  float kminrap = 0.;
  float kmaxrap = 0.;
  float nphibins = 0.;
  float kminphi = 0.;
  float kmaxphi = TwoPI;
  float rapbinwidth = 0.;
  float phibinwidth = 0.;
  float q = 0., e = 0., qinv = 0.;
  float posz = 0., cent = 0., mult = 0.;
  float magField = 0.;
  std::array<float, 9> vTpcRadii = {0.8, 0.10, 0.12, 0.14, 0.16, 0.18, 0.20, 0.22, 0.24};

  // Efficiency object
  TList* ccdbObjRecoEff = nullptr;
  TList* ccdbObjNuaCorr = nullptr;
  struct CorrHist {
    THnSparseF* hRecEffP = nullptr;
    THnSparseF* hRecEffM = nullptr;
    TH3F* hNuaP = nullptr;
    TH3F* hNuaM = nullptr;
  } corrHist;
  int nHistDim = 0;

  void init(InitContext const&)
  {
    // Set CCDB url
    ccdb->setURL(cUrlCCDB.value);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setCreatedNotAfter(nolaterthan.value);

    // Set Density Histogram Attributes
    nrapbins = static_cast<float>(cNEtaBins);
    kminrap = static_cast<float>(-cTrackEtaCut);
    kmaxrap = static_cast<float>(cTrackEtaCut);
    nphibins = static_cast<float>(cNPhiBins);

    rapbinwidth = (kmaxrap - kminrap) / nrapbins;
    phibinwidth = (kmaxphi - kminphi) / nphibins;

    int knrapphibins = static_cast<int>(cNEtaBins) * static_cast<int>(cNPhiBins);
    float kminrapphi = 0.;
    float kmaxrapphi = knrapphibins;

    // Initialize axis specifications
    const AxisSpec axisCheck(1, 0, 1, "");
    const AxisSpec axisChMult(200, 0, 200, "N_{ch}");
    const AxisSpec axisVz(220, -11, 11, "V_{z} (cm)");

    const AxisSpec axisCent(cCentBins, "FT0C (%)");
    const AxisSpec axisPosZ(cPosZBins, "V_{z} (cm)");

    const AxisSpec axisDEta(320, -1.6, 1.6, "#Delta#eta");
    const AxisSpec axisDPhi(640, -PIHalf, 3. * PIHalf, "#Delta#varphi");

    const AxisSpec axisTrackDcaXY{60, -0.15, 0.15, "DCA_{XY}"};
    const AxisSpec axisTrackDcaZ{230, -1.15, 1.15, "DCA_{Z}"};
    const AxisSpec axisdEdx(360, 20, 200, "#frac{dE}{dx}");

    const AxisSpec axisTrackPt(cTrackNPtBins, cTrackMinPt, cTrackMaxPt, "p_{T} (GeV/#it{c})");
    const AxisSpec axisTrackEta(cNEtaBins, -cTrackEtaCut, cTrackEtaCut, "#eta");
    const AxisSpec axisTrackPhi(cNPhiBins, 0., TwoPI, "#varphi (rad)");

    const AxisSpec axisRapPhi(knrapphibins, kminrapphi, kmaxrapphi, "#eta#varphi");

    // Create Histograms
    // Event histograms
    histos.add("Event/h1f_collision_cent", "FT0C(%)", kTH1F, {axisCent});
    histos.add("Event/h1f_collision_Vz", "V_{z}-distribution", kTH1F, {axisVz});

    // Track QA
    histos.add("TrackQA/hPtDcaXY", "DCA_{XY} vs p_{T}", kTH2F, {axisTrackPt, axisTrackDcaXY});
    histos.add("TrackQA/hPtDcaZ", "DCA_{Z} vs p_{T}", kTH2F, {axisTrackPt, axisTrackDcaZ});
    histos.add("TrackQA/hTrackTPCdEdX", "hTrackTPCdEdX", kTH2F, {axisTrackPt, axisdEdx});

    // Two track cut
    histos.add("QA/TwoTrackCut/Before/h2d_n2_detadphi", "#rho_{2}", kTH2D, {axisDEta, axisDPhi});
    histos.add("QA/TwoTrackCut/After/h2d_n2_detadphi", "#rho_{2}", kTH2D, {axisDEta, axisDPhi});
    histos.add("QA/FemtoCut/Before/h2d_n2_detadphi", "#rho_{2}", kTH2D, {axisDEta, axisDPhi});
    histos.add("QA/FemtoCut/After/h2d_n2_detadphi", "#rho_{2}", kTH2D, {axisDEta, axisDPhi});

    // Efficiency Histograms
    // Single Particle Efficiencies
    histos.add("Reco/Efficiency/h1f_n1_pt_P", "#rho_{1}^{#plus}", kTH1F, {axisTrackPt});
    histos.add("Reco/Efficiency/h1f_n1_pt_M", "#rho_{1}^{#minus}", kTH1F, {axisTrackPt});
    histos.add("Reco/Efficiency/h4f_n1_centvzptrap_P", "#rho_{1}^{#plus}", kTHnSparseF, {axisCent, axisPosZ, axisTrackPt, axisTrackEta});
    histos.add("Reco/Efficiency/h4f_n1_centvzptrap_M", "#rho_{1}^{#minus}", kTHnSparseF, {axisCent, axisPosZ, axisTrackPt, axisTrackEta});

    // NUA phi
    histos.add("Reco/NUA/h3f_n1_vzrapphi_P", "#rho_{1}^{#plus}", kTH3F, {axisPosZ, axisTrackEta, axisTrackPhi});
    histos.add("Reco/NUA/h3f_n1_vzrapphi_M", "#rho_{1}^{#minus}", kTH3F, {axisPosZ, axisTrackEta, axisTrackPhi});

    // Correction checks
    histos.add("Reco/h1f_n1_pt_P", "#rho_{1}^{#plus}", kTH1F, {axisTrackPt});
    histos.add("Reco/h1f_n1_pt_M", "#rho_{1}^{#minus}", kTH1F, {axisTrackPt});
    histos.add("Reco/h1f_n1_rap_P", "#rho_{1}^{#plus}", kTH1F, {axisTrackEta});
    histos.add("Reco/h1f_n1_rap_M", "#rho_{1}^{#minus}", kTH1F, {axisTrackEta});
    histos.add("Reco/h1f_n1_phi_P", "#rho_{1}^{#plus}", kTH1F, {axisTrackPhi});
    histos.add("Reco/h1f_n1_phi_M", "#rho_{1}^{#minus}", kTH1F, {axisTrackPhi});

    // Single and Two Particle Densities
    // Rho1 for R2 RapPhi
    histos.add("Reco/h3f_n1_rapphi_P", "#rho_{1}^{#plus}", kTH3F, {axisCent, axisTrackEta, axisTrackPhi});
    histos.add("Reco/h3f_n1_rapphi_M", "#rho_{1}^{#minus}", kTH3F, {axisCent, axisTrackEta, axisTrackPhi});

    // Rho1 for P2 RapPhi
    histos.add("Reco/h3f_pt_rapphi_P", "#rho_{1}^{#plus}", kTH3F, {axisCent, axisTrackEta, axisTrackPhi});
    histos.add("Reco/h3f_pt_rapphi_M", "#rho_{1}^{#minus}", kTH3F, {axisCent, axisTrackEta, axisTrackPhi});

    // Rho2 for R2 RapPhi
    histos.add("Reco/h3f_n2_rapphi_PM", "#rho_{2}^{#plus#minus}", kTH3F, {axisCent, axisRapPhi, axisRapPhi});
    histos.add("Reco/h3f_n2_rapphi_PP", "#rho_{2}^{#plus#plus}", kTH3F, {axisCent, axisRapPhi, axisRapPhi});
    histos.add("Reco/h3f_n2_rapphi_MM", "#rho_{2}^{#minus#minus}", kTH3F, {axisCent, axisRapPhi, axisRapPhi});

    // Rho2 for P2 RapPhi
    histos.add("Reco/h3f_ptpt_rapphi_PM", "#rho_{2}^{#plus#minus}", kTH3F, {axisCent, axisRapPhi, axisRapPhi});
    histos.add("Reco/h3f_ptpt_rapphi_PP", "#rho_{2}^{#plus#plus}", kTH3F, {axisCent, axisRapPhi, axisRapPhi});
    histos.add("Reco/h3f_ptpt_rapphi_MM", "#rho_{2}^{#minus#minus}", kTH3F, {axisCent, axisRapPhi, axisRapPhi});
    histos.add("Reco/h3f_npt_rapphi_PM", "#rho_{2}^{#plus#minus}", kTH3F, {axisCent, axisRapPhi, axisRapPhi});
    histos.add("Reco/h3f_npt_rapphi_PP", "#rho_{2}^{#plus#plus}", kTH3F, {axisCent, axisRapPhi, axisRapPhi});
    histos.add("Reco/h3f_npt_rapphi_MM", "#rho_{2}^{#minus#minus}", kTH3F, {axisCent, axisRapPhi, axisRapPhi});
    histos.add("Reco/h3f_ptn_rapphi_PM", "#rho_{2}^{#plus#minus}", kTH3F, {axisCent, axisRapPhi, axisRapPhi});
    histos.add("Reco/h3f_ptn_rapphi_PP", "#rho_{2}^{#plus#plus}", kTH3F, {axisCent, axisRapPhi, axisRapPhi});
    histos.add("Reco/h3f_ptn_rapphi_MM", "#rho_{2}^{#minus#minus}", kTH3F, {axisCent, axisRapPhi, axisRapPhi});

    histos.addClone("Reco/", "McGen/");

    // MC Generated Histograms
    if (doprocessMCRecoGen) {
      // McGen Histos
      histos.add("McGen/h1f_collision_recgen", "# of Reco Collision Associated to One Mc Generator Collision", kTH1F, {axisChMult});
      histos.add("McGen/h2f_collision_posZ", "V_{z}-distribution", kTH2F, {axisVz, axisVz});
      histos.add("McGen/h2f_collision_cent", "FT0M Centrality", kTH2F, {axisCent, axisCent});
    }

    // Load correction factor
    if (cGetCorrectionFlag) {
      // Set CCDB url
      ccdb->setURL(cUrlCCDB.value);
      ccdb->setCaching(true);
      ccdb->setLocalObjectValidityChecking();
      ccdb->setCreatedNotAfter(nolaterthan.value);

      // Get CCDB object
      ccdbObjRecoEff = ccdb->getForTimeStamp<TList>(cPathCCDBRecoEff.value, nolaterthan.value);
      ccdbObjNuaCorr = ccdb->getForTimeStamp<TList>(cPathCCDBNuaCorr.value, nolaterthan.value);

      // Load reco eff corrections
      loadRecoEfficiencyHistograms();
    }
  }

  // Load reco efficiency histograms
  void loadRecoEfficiencyHistograms()
  {
    // Efficiency correction histograms
    corrHist.hRecEffP = dynamic_cast<THnSparseF*>(ccdbObjRecoEff->FindObject("h_RecEff_P"));
    corrHist.hRecEffM = dynamic_cast<THnSparseF*>(ccdbObjRecoEff->FindObject("h_RecEff_M"));
    if (!corrHist.hRecEffP || !corrHist.hRecEffM) {
      LOGF(fatal, "CCDB efficiency object doesn't exist !");
    }
    nHistDim = corrHist.hRecEffP->GetNdimensions();
    LOGF(info, "Efficiency correction histogram dimensions: %d", nHistDim);

    // Nua correction histograms
    corrHist.hNuaP = dynamic_cast<TH3F*>(ccdbObjNuaCorr->FindObject("h_Nua_P"));
    corrHist.hNuaM = dynamic_cast<TH3F*>(ccdbObjNuaCorr->FindObject("h_Nua_M"));
    if (!corrHist.hNuaP || !corrHist.hNuaM) {
      LOGF(fatal, "CCDB NUA object doesn't exist !");
    }
  }

  // Get magnetic field
  float getMagneticField(int64_t const& timestamp)
  {
    if (grpo == nullptr) {
      grpo = ccdb->getForTimeStamp<o2::parameters::GRPMagField>("GLO/Config/GRPMagField", timestamp);
      if (grpo == nullptr) {
        LOGF(fatal, "GRP object not found for timestamp %llu", timestamp);
        return 1;
      }
    }
    auto field = std::lround(5.f * grpo->getL3Current() / 30000.f);
    return 0.1 * field;
  }

  template <typename C>
  bool selCollision(C const& col)
  {
    posz = col.posZ();
    if (std::abs(posz) >= cZVtxCut) { // VtxZ selection
      return false;
    }

    if (cSel8Trig && !col.sel8()) { // Sel8 selection
      return false;
    }

    cent = col.centFT0C();
    if (cent <= cMinCent || cent >= cMaxCent) { // Centrality selection
      return false;
    }

    if (col.ft0cOccupancyInTimeRange() < cMinOccupancy || col.ft0cOccupancyInTimeRange() > cMaxOccupancy) { // Occupancy cut
      return false;
    }

    if (cPileupReject && !col.selection_bit(aod::evsel::kNoSameBunchPileup)) { // Pile-up rejection
      return false;
    }

    if (cZVtxTimeDiff && !col.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV)) { // ZvtxFT0 vs PV
      return false;
    }

    if (cIsGoodITSLayers && !col.selection_bit(aod::evsel::kIsGoodITSLayersAll)) { // All ITS layer active
      return false;
    }

    return true;
  }

  // Track Selection
  template <typename T>
  bool selectTrack(T const& track)
  {
    // Kinematic selection
    if (track.pt() <= cTrackMinPt || track.pt() >= cTrackMaxPt || std::abs(track.eta()) >= cTrackEtaCut) {
      return false;
    }

    // Global track selection
    if (cTrackGlobal && !track.isGlobalTrackWoDCA()) {
      return false;
    }

    // DCA selection
    if (std::abs(track.dcaXY()) >= cTrackDcaXYCut || std::abs(track.dcaZ()) >= cTrackDcaZCut) {
      return false;
    }

    // Electron rejection
    if (std::abs(track.tpcNSigmaPi()) > cTpcRejCut && std::abs(track.tpcNSigmaKa()) > cTpcRejCut && std::abs(track.tpcNSigmaPr()) > cTpcRejCut && track.tpcNSigmaEl() > cTpcElRejCutMin && track.tpcNSigmaEl() < cTpcElRejCutMax) {
      return false;
    }

    // All selection passed
    return true;
  }

  template <RecGenType rec_gen, typename T, typename S>
  float getCorrectionFactor(T const& track, S const& sign)
  {
    if (!cGetCorrectionFlag) {
      return 1.;
    }

    if constexpr (rec_gen == kGen) {
      return 1.;
    }

    std::array<int, 4> binarray{};
    std::array<float, 4> v = {cent, posz, track.pt(), track.eta()};
    if (sign > 0) {
      for (int i = 0; i < nHistDim; ++i) {
        binarray[i] = corrHist.hRecEffP->GetAxis(i)->FindBin(v[i]);
      }
      return corrHist.hRecEffP->GetBinContent(corrHist.hRecEffP->GetBin(binarray.data()));
    }

    if (sign < 0) {
      for (int i = 0; i < nHistDim; ++i) {
        binarray[i] = corrHist.hRecEffM->GetAxis(i)->FindBin(v[i]);
      }
      return corrHist.hRecEffM->GetBinContent(corrHist.hRecEffM->GetBin(binarray.data()));
    }

    return 1.;
  }

  template <RecGenType rec_gen, typename T, typename S>
  float getNuaCorrectionFactor(T const& track, S const& sign)
  {
    if (!cGetNuaCorrectionFlag) {
      return 1.;
    }

    if constexpr (rec_gen == kGen) {
      return 1.;
    }

    if (sign > 0) {
      return corrHist.hNuaP->GetBinContent(corrHist.hNuaP->FindBin(posz, track.eta(), track.phi()));
    }

    if (sign < 0) {
      return corrHist.hNuaM->GetBinContent(corrHist.hNuaM->FindBin(posz, track.eta(), track.phi()));
    }

    return 1.;
  }

  template <typename T>
  bool isClosePair(T const& p1, T const& p2)
  {
    // Before
    histos.fill(HIST("QA/TwoTrackCut/Before/h2d_n2_detadphi"), p1.eta() - p2.eta(), RecoDecay::constrainAngle(p1.phi() - p2.phi(), -PIHalf));

    // DPhiStar average over TPC
    float dphistar = 0., n = 0.;
    for (auto const& radii : vTpcRadii) {
      float arg1 = 0.15 * magField * radii / p1.pt();
      float arg2 = 0.15 * magField * radii / p2.pt();
      if (std::abs(arg1) < 1.0 && std::abs(arg2) < 1.0) {
        dphistar += (p1.phi() - p2.phi() - (p1.sign() * std::abs(std::asin(arg1))) + (p2.sign() * std::abs(std::asin(arg2))));
        ++n;
      }
    }

    // Nan check
    if (n == 0) {
      return false;
    }

    // DPhistar
    dphistar = RecoDecay::constrainAngle(dphistar / n, -PIHalf);

    // DEta
    float deta = p1.eta() - p2.eta();

    // Return flag
    bool retFlag = (std::abs(deta) < cDEtaCut && std::abs(dphistar) < cDPhiStarCut);

    // Before
    if (!retFlag) {
      histos.fill(HIST("QA/TwoTrackCut/After/h2d_n2_detadphi"), deta, RecoDecay::constrainAngle(p1.phi() - p2.phi(), -PIHalf));
    }

    return retFlag;
  }

  template <typename T>
  bool isCloseQinv(T const& p1, T const& p2)
  {
    // Before
    histos.fill(HIST("QA/FemtoCut/Before/h2d_n2_detadphi"), p1.eta() - p2.eta(), RecoDecay::constrainAngle(p1.phi() - p2.phi(), -PIHalf));
    float dpt = std::abs(p1.pt() - p2.pt());
    bool retFlag = (dpt < cFemtoCut);

    if (!retFlag) {
      histos.fill(HIST("QA/FemtoCut/After/h2d_n2_detadphi"), p1.eta() - p2.eta(), RecoDecay::constrainAngle(p1.phi() - p2.phi(), -PIHalf));
    }

    return retFlag;
  }

  template <RecGenType rec_gen, typename T, typename S>
  void fillPairHist(T const& trk_1, T const& trk_2, S const& sign_1, S const& sign_2)
  {
    // Check for same index
    if (trk_1.index() == trk_2.index()) {
      return;
    }

    // Close pair rejection
    if constexpr (rec_gen == kRec) {
      if (trk_1.sign() * trk_2.sign() < 0) {
        // Close pair
        if (cApplyTwoTrackCut && isClosePair(trk_1, trk_2)) {
          return;
        }
        // Femto selection
        if (cApplyFemtoSel && isCloseQinv(trk_1, trk_2)) {
          return;
        }
      }
    }

    // Reco/Gen Dir
    static constexpr auto SubDirRecGen = std::array{"Reco/", "McGen/"};

    // RapPhi bins
    const auto rapbin1 = static_cast<int>((trk_1.eta() - kminrap) / rapbinwidth);
    const auto rapbin2 = static_cast<int>((trk_2.eta() - kminrap) / rapbinwidth);
    const auto phibin1 = static_cast<int>(trk_1.phi() / phibinwidth);
    const auto phibin2 = static_cast<int>(trk_2.phi() / phibinwidth);

    float corfac = getCorrectionFactor<rec_gen>(trk_1, sign_1) * getCorrectionFactor<rec_gen>(trk_2, sign_2);

    if (rapbin1 >= 0 && rapbin2 >= 0 && phibin1 >= 0 && phibin2 >= 0 && rapbin1 < nrapbins && rapbin2 < nrapbins && phibin1 < nphibins && phibin2 < nphibins) {

      int rapphix = rapbin1 * nphibins + phibin1;
      int rapphiy = rapbin2 * nphibins + phibin2;

      if ((sign_1 > 0 && sign_2 < 0) || (sign_1 < 0 && sign_2 > 0)) {
        histos.fill(HIST(SubDirRecGen[rec_gen]) + HIST("h3f_n2_rapphi_PM"), cent, rapphix + 0.5, rapphiy + 0.5, corfac);
        histos.fill(HIST(SubDirRecGen[rec_gen]) + HIST("h3f_ptpt_rapphi_PM"), cent, rapphix + 0.5, rapphiy + 0.5, trk_1.pt() * trk_2.pt() * corfac);
        histos.fill(HIST(SubDirRecGen[rec_gen]) + HIST("h3f_npt_rapphi_PM"), cent, rapphix + 0.5, rapphiy + 0.5, trk_1.pt() * corfac);
        histos.fill(HIST(SubDirRecGen[rec_gen]) + HIST("h3f_ptn_rapphi_PM"), cent, rapphix + 0.5, rapphiy + 0.5, trk_2.pt() * corfac);
      } else {
        if (sign_1 > 0 && sign_2 > 0) {
          histos.fill(HIST(SubDirRecGen[rec_gen]) + HIST("h3f_n2_rapphi_PP"), cent, rapphix + 0.5, rapphiy + 0.5, corfac);
          histos.fill(HIST(SubDirRecGen[rec_gen]) + HIST("h3f_ptpt_rapphi_PP"), cent, rapphix + 0.5, rapphiy + 0.5, trk_1.pt() * trk_2.pt() * corfac);
          histos.fill(HIST(SubDirRecGen[rec_gen]) + HIST("h3f_npt_rapphi_PP"), cent, rapphix + 0.5, rapphiy + 0.5, trk_1.pt() * corfac);
          histos.fill(HIST(SubDirRecGen[rec_gen]) + HIST("h3f_ptn_rapphi_PP"), cent, rapphix + 0.5, rapphiy + 0.5, trk_2.pt() * corfac);
        } else if (sign_1 < 0 && sign_2 < 0) {
          histos.fill(HIST(SubDirRecGen[rec_gen]) + HIST("h3f_n2_rapphi_MM"), cent, rapphix + 0.5, rapphiy + 0.5, corfac);
          histos.fill(HIST(SubDirRecGen[rec_gen]) + HIST("h3f_ptpt_rapphi_MM"), cent, rapphix + 0.5, rapphiy + 0.5, trk_1.pt() * trk_2.pt() * corfac);
          histos.fill(HIST(SubDirRecGen[rec_gen]) + HIST("h3f_npt_rapphi_MM"), cent, rapphix + 0.5, rapphiy + 0.5, trk_1.pt() * corfac);
          histos.fill(HIST(SubDirRecGen[rec_gen]) + HIST("h3f_ptn_rapphi_MM"), cent, rapphix + 0.5, rapphiy + 0.5, trk_2.pt() * corfac);
        }
      }
    }
  }

  template <typename T>
  void fillTrackQA(T const& track)
  {
    histos.fill(HIST("TrackQA/hPtDcaZ"), track.pt(), track.dcaZ());
    histos.fill(HIST("TrackQA/hPtDcaXY"), track.pt(), track.dcaXY());
    histos.fill(HIST("TrackQA/hTrackTPCdEdX"), track.pt(), track.tpcSignal());
  }

  template <RecGenType rec_gen, typename T, typename S>
  void fillSingleHist(T const& track, S const& sign)
  {
    // Hist array
    static constexpr auto SubDirRecGen = std::array{"Reco/", "McGen/"};

    // Correction factor
    float corrFact = getCorrectionFactor<rec_gen>(track, sign);
    float nuaCorr = getNuaCorrectionFactor<rec_gen>(track, sign);

    // Histograms
    if (sign > 0) {
      // Corrections
      histos.fill(HIST(SubDirRecGen[rec_gen]) + HIST("Efficiency/h1f_n1_pt_P"), track.pt());
      histos.fill(HIST(SubDirRecGen[rec_gen]) + HIST("Efficiency/h4f_n1_centvzptrap_P"), cent, posz, track.pt(), track.eta());
      histos.fill(HIST(SubDirRecGen[rec_gen]) + HIST("NUA/h3f_n1_vzrapphi_P"), posz, track.eta(), track.phi());

      // Checks
      histos.fill(HIST(SubDirRecGen[rec_gen]) + HIST("h1f_n1_pt_P"), track.pt(), corrFact);
      histos.fill(HIST(SubDirRecGen[rec_gen]) + HIST("h1f_n1_rap_P"), track.eta(), corrFact);
      histos.fill(HIST(SubDirRecGen[rec_gen]) + HIST("h1f_n1_phi_P"), track.phi(), nuaCorr);

      // R2 Rho1 (Eta,Phi)
      histos.fill(HIST(SubDirRecGen[rec_gen]) + HIST("h3f_n1_rapphi_P"), cent, track.eta(), track.phi(), corrFact);

      // P2 Rho1 (Eta,Phi)
      histos.fill(HIST(SubDirRecGen[rec_gen]) + HIST("h3f_pt_rapphi_P"), cent, track.eta(), track.phi(), track.pt() * corrFact);
    } else if (sign < 0) {
      // Corrections
      histos.fill(HIST(SubDirRecGen[rec_gen]) + HIST("Efficiency/h1f_n1_pt_M"), track.pt());
      histos.fill(HIST(SubDirRecGen[rec_gen]) + HIST("Efficiency/h4f_n1_centvzptrap_M"), cent, posz, track.pt(), track.eta());
      histos.fill(HIST(SubDirRecGen[rec_gen]) + HIST("NUA/h3f_n1_vzrapphi_M"), posz, track.eta(), track.phi());

      // Checks
      histos.fill(HIST(SubDirRecGen[rec_gen]) + HIST("h1f_n1_pt_M"), track.pt(), corrFact);
      histos.fill(HIST(SubDirRecGen[rec_gen]) + HIST("h1f_n1_rap_M"), track.eta(), corrFact);
      histos.fill(HIST(SubDirRecGen[rec_gen]) + HIST("h1f_n1_phi_M"), track.phi(), nuaCorr);

      // R2 Rho1 (Eta,Phi)
      histos.fill(HIST(SubDirRecGen[rec_gen]) + HIST("h3f_n1_rapphi_M"), cent, track.eta(), track.phi(), corrFact);

      // P2 Rho1 (Eta,Phi)
      histos.fill(HIST(SubDirRecGen[rec_gen]) + HIST("h3f_pt_rapphi_M"), cent, track.eta(), track.phi(), track.pt() * corrFact);
    }
  }

  template <DMCType dmc, typename B, typename C, typename T>
  void analyzeCollision(B const&, C const& collision, T const& tracks)
  {
    // Select Collision (Only for Data... McRec has been selected already !!!)
    if constexpr (dmc == kData) {
      if (!selCollision(collision)) {
        return;
      }
    }

    // Get magnetic field
    magField = getMagneticField(collision.template foundBC_as<B>().timestamp());

    // Fill Event QA
    histos.fill(HIST("Event/h1f_collision_cent"), cent);
    histos.fill(HIST("Event/h1f_collision_Vz"), collision.posZ());

    // Loop over tracks
    for (auto const& track_1 : tracks) {
      // Check for Mc matched particle
      if constexpr (dmc == kMC) {
        if (!track_1.has_mcParticle()) {
          continue;
        }
      }

      // Select track
      if (!selectTrack(track_1)) {
        continue;
      }

      // Fill QA
      fillTrackQA(track_1);

      // Fill singles
      fillSingleHist<kRec>(track_1, track_1.sign());

      // Pair
      for (auto const& track_2 : tracks) {
        // Check for Mc matched particle
        if constexpr (dmc == kMC) {
          if (!track_2.has_mcParticle()) {
            continue;
          }
        }

        // Select track
        if (!selectTrack(track_2)) {
          continue;
        }

        // Fill pair histograms
        fillPairHist<kRec>(track_1, track_2, track_1.sign(), track_2.sign());
      }
    }
  }

  template <typename M>
  bool selectMcParticle(M const& mcpart)
  {
    // Check for Primary Charged Particle
    if (!mcpart.isPhysicalPrimary()) {
      return false;
    }

    // Check pdg info
    auto pdgCode = mcpart.pdgCode();
    auto pdgInfo = pdg->GetParticle(pdgCode);
    if (pdgInfo == nullptr) { // particle with unknown pdg code
      return false;
    }

    // Remove electron
    if (pdgCode == kElectron || pdgCode == kPositron) {
      return false;
    }

    // Apply kinematic selection
    if (mcpart.pt() <= cTrackMinPt || mcpart.pt() >= cTrackMaxPt || std::abs(mcpart.eta()) >= cTrackEtaCut) {
      return false;
    }

    // All selection passed
    return true;
  }

  // MC Gen analysis
  template <typename M, typename P>
  void analyzeMcGenCollision(M const&, P const& mcParticles)
  {
    // Loop 1 over MC particles
    for (auto const& mcpart_1 : mcParticles) {
      // Select mc particle
      if (!selectMcParticle(mcpart_1)) {
        continue;
      }

      // Fill single particle densities
      fillSingleHist<kGen>(mcpart_1, pdg->GetParticle(mcpart_1.pdgCode())->Charge() / 3);

      // Loop 2
      for (auto const& mcpart_2 : mcParticles) {
        // Select mc particle
        if (!selectMcParticle(mcpart_2)) {
          continue;
        }

        // Fill pair densities
        fillPairHist<kGen>(mcpart_1, mcpart_2, pdg->GetParticle(mcpart_1.pdgCode())->Charge() / 3, pdg->GetParticle(mcpart_2.pdgCode())->Charge() / 3);
      }
    }
  }

  // MC Reco-Gen analysis
  template <DMCType dmc, typename M, typename C, typename B, typename T, typename P>
  void analyzeMcRecoGen(M const& mcCollision, C const& collisions, B const& bc, T const& tracks, P const& mcParticles)
  {
    // Number of Rec Collisions Associated to the McGen Collision
    int nRecCols = collisions.size();
    if (nRecCols != 0) {
      histos.fill(HIST("McGen/h1f_collision_recgen"), nRecCols);
    }
    // Do not analyze if more than one reco collision is accociated to one mc gen collision
    if (nRecCols != 1) {
      return;
    }
    // Check the reco collision
    if (!collisions.begin().has_mcCollision() || !selCollision(collisions.begin()) || collisions.begin().mcCollisionId() != mcCollision.globalIndex()) {
      return;
    }
    histos.fill(HIST("McGen/h2f_collision_posZ"), mcCollision.posZ(), collisions.begin().posZ());
    auto tracksThisCollision = tracks.sliceBy(tracksPerCollision, collisions.begin().globalIndex());
    analyzeCollision<dmc>(bc, collisions.begin(), tracksThisCollision);
    analyzeMcGenCollision(mcCollision, mcParticles);
  }

  // BF, Collision and Track Table
  using BCsRun3 = soa::Join<aod::BCsWithTimestamps, aod::Run3MatchedToBCSparse>;
  using Collisions = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs, aod::CentFT0Ms, aod::CentFV0As, aod::MultsExtra>;
  using Tracks = soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra, aod::TracksDCA, aod::TOFSignal, aod::pidTPCEl, aod::pidTPCPi, aod::pidTOFPi, aod::pidTPCKa, aod::pidTOFKa, aod::pidTPCPr, aod::pidTOFPr>;
  using TracksMC = soa::Join<Tracks, aod::McTrackLabels>;

  SliceCache cache;
  Preslice<TracksMC> tracksPerCollision = aod::track::collisionId;

  void processDummy(aod::Collisions const&) {}

  PROCESS_SWITCH(ChargeBalanceFunction, processDummy, "Dummy process", true);

  void processChargedData(Collisions::iterator const& collision, BCsRun3 const& bc, Tracks const& tracks)
  {
    analyzeCollision<kData>(bc, collision, tracks);
  }

  PROCESS_SWITCH(ChargeBalanceFunction, processChargedData, "Charged particle BF process", false);

  void processMCRecoGen(aod::McCollisions::iterator const& mcCollision,
                        soa::SmallGroups<soa::Join<Collisions, aod::McCollisionLabels>> const& collisions, BCsRun3 const& bc,
                        TracksMC const& tracks,
                        aod::McParticles const& mcParticles)
  {
    analyzeMcRecoGen<kMC>(mcCollision, collisions, bc, tracks, mcParticles);
  }

  PROCESS_SWITCH(ChargeBalanceFunction, processMCRecoGen, "Process for MC RecoGen", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<ChargeBalanceFunction>(cfgc)};
}
