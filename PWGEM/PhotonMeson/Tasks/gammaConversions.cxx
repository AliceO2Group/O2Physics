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

/// \brief perform photon conversion analysis on V0 candidates from aod::StoredV0Datas
/// dependencies: o2-analysis-lf-lambdakzerobuilder
/// \author stephan.friedrich.stiefelmaier@cern.ch

#include "PWGEM/PhotonMeson/Tasks/gammaConversions.h"
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGEM/PhotonMeson/Utils/gammaConvDefinitions.h"

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/runDataProcessing.h"

#include <TVector3.h>
#include <TMath.h> // for ATan2, Cos, Sin, Sqrt

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using V0DatasAdditional = soa::Join<aod::V0Datas, aod::V0Recalculated>;

// using collisionEvSelIt = soa::Join<aod::Collisions, aod::EvSels>::iterator;
struct GammaConversions {

  // event cuts
  /*
  Configurable<bool> fDoEventSel{"fDoEventSel", 0, "demand sel7 for events"};
  Configurable<float> fCentMin{"fCentMin", 0.0, "lower bound of centrality selection"};
  Configurable<float> fCentMax{"fCentMax", 100.0, "upper bound of centrality selection"};
  */

  // track cuts
  Configurable<float> fTrackEtaMax{"fTrackEtaMax", 0.8, "accepted track eta range"};
  Configurable<float> fTruePhotonEtaMax{"fTruePhotonEtaMax", fTrackEtaMax, "eta range for true mc photons that serve for validation"};
  Configurable<float> fTrackPtMin{"fTrackPtMin", 0.04, "minimum daughter track pt"};
  Configurable<float> fPIDnSigmaElectronMin{"fPIDnSigmaElectronMin", -3., "minimum sigma electron PID for V0 daughter tracks. Set 0 to disable."};
  Configurable<float> fPIDnSigmaElectronMax{"fPIDnSigmaElectronMax", 3., "maximum sigma electron PID for V0 daughter tracks"};
  Configurable<float> fPIDPionRejectionPMin{"fPIDPionRejectionPMin", 0.4, "minimum track momentum to apply any pion rejection"};                              // case 7:  // 0.4 GeV
  Configurable<float> fPIDPionRejectionPBoarder{"fPIDPionRejectionPBoarder", 8., "border between low and high momentum pion rejection"};                      // case 7:  // 8. GeV
  Configurable<float> fPIDnSigmaAbovePionLineLowPMin{"fPIDnSigmaAbovePionLineLowPMin", -10., "minimum sigma to be over the pion line for low momentum tracks"}; // case 4: 3.0sigma, 1.0 sigma at high momentum
  Configurable<float> fPIDnSigmaAbovePionLineHighPMin{"fPIDnSigmaAbovePionLineHighPMin", -10., "minimum sigma to be over the pion line for high momentum tracks"};
  Configurable<float> fMinTPCFoundOverFindableCls{"fMinTPCNClsFoundOverFindable", 0.3, "minimum ratio found tpc clusters over findable"}; // case 9:  // 0.6
  Configurable<float> fMinTPCCrossedRowsOverFindableCls{"fMinTPCCrossedRowsOverFindableCls", 0.0, "minimum ratio TPC crossed rows over findable clusters"};

  // V0 cuts
  Configurable<float> fV0CosPAngleMin{"fV0CosPAngleMin", 0.85, "Set negative to disable. Minimum cosinus of the pointing angle"}; // case 4
  Configurable<float> fV0RMin{"fV0RMin", 0., "minimum conversion radius of the V0s"};
  Configurable<float> fV0RMax{"fV0RMax", 180., "maximum conversion radius of the V0s"};
  Configurable<float> fV0PhotonAsymmetryMax{"fV0PhotonAsymmetryMax", 1.0, "maximum photon asymetry. Set negative do disable cut."};
  Configurable<float> fV0PsiPairMax{"fV0PsiPairMax", -0.1, "maximum psi angle of the track pair. Set negative do disable cut. "};
  Configurable<float> fV0QtPtMultiplicator{"fV0QtPtMultiplicator", 0.125, "Multiply pt of V0s by this value to get the 2nd denominator in the armenteros cut. The products maximum value is fV0QtMax."};
  Configurable<float> fV0QtMax{"fV0QtMax", 0.050, "the maximum value of the product, that is the maximum qt"};
  Configurable<float> LineCutZ0{"fLineCutZ0", 7.0, "The offset for the linecute used in the Z vs R plot"};
  Configurable<float> LineCutZRSlope{"LineCutZRSlope", (float)TMath::Tan(2 * TMath::ATan(TMath::Exp(-fTruePhotonEtaMax))), "The slope for the line cut"};

  std::map<ePhotonCuts, std::string> fPhotonCutLabels{
    {ePhotonCuts::kV0In, "kV0In"},
    {ePhotonCuts::kTrackEta, "kTrackEta"},
    {ePhotonCuts::kTrackPt, "kTrackPt"},
    {ePhotonCuts::kElectronPID, "kElectronPID"},
    {ePhotonCuts::kPionRejLowMom, "kPionRejLowMom"},
    {ePhotonCuts::kPionRejHighMom, "kPionRejHighMom"},
    {ePhotonCuts::kTPCFoundOverFindableCls, "kTPCFoundOverFindableCls"},
    {ePhotonCuts::kTPCCrossedRowsOverFindableCls, "kTPCCrossedRowsOverFindableCls"},
    {ePhotonCuts::kV0Radius, "kV0Radius"},
    {ePhotonCuts::kArmenteros, "kArmenteros"},
    {ePhotonCuts::kPsiPair, "kPsiPair"},
    {ePhotonCuts::kCosinePA, "kCosinePA"},
    {ePhotonCuts::kRZLine, "kRZLine"},
    {ePhotonCuts::kV0Out, "kV0Out"}};

  std::map<eV0McValidation, std::string> fV0McValidationLabels{
    {eV0McValidation::kV0in, "kV0in"},
    {eV0McValidation::kFakeV0, "kFakeV0"},
    {eV0McValidation::kMcMotherIn, "kMcMotherIn"},
    {eV0McValidation::kNoPhysicalPrimary, "kNoPhysicalPrimary"},
    {eV0McValidation::kNoPhoton, "kNoPhoton"},
    {eV0McValidation::kOutsideMCEtaAcc, "kOutsideMCEtaAcc"},
    {eV0McValidation::kMcValidatedPhotonOut, "kMcValidatedPhotonOut"},
    {eV0McValidation::kMcValAfterRecCuts, "kMcValAfterRecCuts"}};

  std::vector<std::string> fHistoSuffixes{"_MCRec", "_MCTrue", "_MCVal", "_Res"};

  tHistoRegistry fMyRegistry{""};
  HistogramRegistry fHistogramRegistry{"fHistogramRegistry"};

  void init(InitContext const&)
  {
    // make axis logarithmic
    gAxis_pT_log.makeLogaritmic();

    // Declarations for histogram
    // collision histograms
    std::vector<MyHistogramSpec> lCollisionHistoDefinitions{
      {"hCollisionZ", "hCollisionZ;z (cm);counts", {HistType::kTH1F, {gAxis_zColl}}}};

    // track histograms
    std::vector<MyHistogramSpec> lTrackHistoDefinitions{
      {"hTrackPt", "hTrackPt;p_{T} (GeV/c);counts", {HistType::kTH1F, {gAxis_pT}}},
      {"hTrackEta", "hTrackEta;#eta;counts", {HistType::kTH1F, {gAxis_eta}}},
      {"hTrackPhi", "hTrackPhi;#phi (rad);counts", {HistType::kTH1F, {gAxis_phi}}},
      {"hTPCFoundOverFindableCls", "hTPCFoundOverFindableCls;TPCFoundOverFindableCls;counts", {HistType::kTH1F, {{800, 0.9f, 1.01f}}}},
      {"hTPCCrossedRowsOverFindableCls", "hTPCCrossedRowsOverFindableCls;TPCCrossedRowsOverFindableCls;counts", {HistType::kTH1F, {{800, 0.8f, 1.5f}}}},
      {"hTPCdEdx", "hTPCdEdx;p (GeV/c);TPCdEdx", {HistType::kTH2F, {gAxis_pT_log, {800, 0.f, 200.f}}}},
      {"hTPCdEdxSigEl", "hTPCdEdxSigEl;p (GeV/c);TPCdEdxSigEl", {HistType::kTH2F, {gAxis_pT_log, gAxis_TPCdEdxSig}}},
      {"hTPCdEdxSigPi", "hTPCdEdxSigPi;p (GeV/c);TPCdEdxSigPi", {HistType::kTH2F, {gAxis_pT_log, gAxis_TPCdEdxSig}}}};

    // v0 histograms
    std::vector<MyHistogramSpec> lV0HistoDefinitions{
      {"hPt", "hPt;p_{T} (GeV/c);counts", {HistType::kTH1F, {gAxis_pT}}},
      {"hEta", "hEta;#eta;counts", {HistType::kTH1F, {gAxis_eta}}},
      {"hPhi", "hPhi;#phi (rad);counts", {HistType::kTH1F, {gAxis_phi}}},
      {"hConvPointR", "hConvPointR;conversion radius (cm);counts", {HistType::kTH1F, {gAxis_r}}},
      {"hConvPointZ", "hConvPointZ;conversion radius (cm);counts", {HistType::kTH1F, {gAxis_xyz}}},
      {"hArmenteros", "hArmenteros;#alpha;q_{T} (GeV/c)", {HistType::kTH2F, {{800, -1.f, 1.f}, gAxis_pT_armenteros}}},
      {"hPsiPt", "hPsiPt;#Psi;p_{T} (GeV/c)", {HistType::kTH2F, {gAxis_eta, gAxis_pT}}},
      {"hCosPAngle", "hCosPAngle;CosPAngle;counts", {HistType::kTH1F, {{800, 0.99f, 1.005f}}}},
      {"hRVsZ", "hRVsZ;R (cm);z (cm)", {HistType::kTH2F, {gAxis_r, gAxis_xyz}}},
      {"hpeDivpGamma", "hpeDivpGamma;p_{e}p_{#gamma};counts", {HistType::kTH1F, {{200, 0.f, 1.f}}}}};

    // recalculated conversion Point for V0, only Rec and MCVal need this
    std::vector<MyHistogramSpec> lV0HistoDefinitions_recalculated{
      {"hConvPointR_recalc", "hConvPointR_recalc;conversion radius (cm);counts", {HistType::kTH1F, {gAxis_r}}},
      {"hConvPointZ_recalc", "hConvPointZ_recalc;conversion radius (cm);counts", {HistType::kTH1F, {gAxis_xyz}}},
    };

    // only in mc
    // resolution histos
    std::vector<MyHistogramSpec> lV0ResolutionHistoDefinitions{
      {"hPtRes", "hPtRes_Rec-MC;#Delta p_T (GeV/c);counts", {HistType::kTH1F, {{800, -5.f, 5.f}}}},
      {"hEtaRes", "hEtaRes_Rec-MC;#Delta #eta;counts", {HistType::kTH1F, {gAxis_radRes}}},
      {"hPhiRes", "hPhiRes_Rec-MC;#Delta #phi (rad);counts", {HistType::kTH1F, {gAxis_radRes}}},
      {"hConvPointXYResVsXY", "hConvPointXYResVsXY_Rec-MC;R (cm);#Delta R (cm)", {HistType::kTH2F, {gAxis_r, gAxis_dr}}}, // maybe change axis later on
      {"hConvPointXYResVsXY_recalc", "hConvPointXYResVsXY_recalc_Rec-MC;R (cm);#Delta R (cm)", {HistType::kTH2F, {gAxis_r, gAxis_dr}}},
      {"hConvPointZResVsZ", "hConvPointZResVsZ_Rec-MC;Z (cm);#Delta Z (cm)", {HistType::kTH2F, {gAxis_xyz, gAxis_dr}}},
      {"hConvPointZResVsZ_recalc", "hConvPointZResVsZ_recalc_Rec-MC;Z (cm);#Delta Z (cm)", {HistType::kTH2F, {gAxis_xyz, gAxis_dr}}},
    };

    // think of better name
    std::vector<MyHistogramSpec> lSpecialHistoDefinitions{
      {"hV0Selection", "hV0Selection;V0 categories;counts", {HistType::kTH1I, {{14, -0.0f, 13.5f}}}},
      {"hV0McValidation", "hV0McValidation;V0 categories;counts", {HistType::kTH1I, {{13, -0.0f, 12.5f}}}, false /*callSumw2*/, true /*dataOnly_*/}};

    auto addLablesToHisto = [](auto const& theContainer, std::string const& theHistoName, auto const& theLables) {
      auto lHisto = theContainer.find(theHistoName);
      if (lHisto != theContainer.end()) {
        TAxis* lXaxis = std::get<std::shared_ptr<TH1>>(lHisto->second)->GetXaxis();
        for (auto& lPairIt : theLables) {
          lXaxis->SetBinLabel(static_cast<int>(lPairIt.first) + 1, lPairIt.second.data());
        }
      }
    };

    if (doprocessRec && doprocessMc) {
      LOGF(fatal, "Cannot enable doprocessRec and doprocessMc at the same time. Please choose one.");
    }

    if (doprocessRec) {
      fHistoSuffixes[0] = "Rec";
    }

    fMyRegistry.mV0.mSpecialHistos.addHistosToOfficalRegistry(
      fHistogramRegistry,
      lSpecialHistoDefinitions,
      nullptr /*theSuffix*/,
      doprocessRec /*theCheckDataOnly*/);

    // do some labeling
    addLablesToHisto(fMyRegistry.mV0.mSpecialHistos.mContainer, "hV0Selection", fPhotonCutLabels);
    if (doprocessMc) {
      addLablesToHisto(fMyRegistry.mV0.mSpecialHistos.mContainer, "hV0McValidation", fV0McValidationLabels);
    }

    for (size_t iBARecCuts = 0; iBARecCuts < 2; ++iBARecCuts) {
      for (size_t iMcKind = 0; iMcKind < (doprocessMc ? 4 : 1); ++iMcKind) {
        std::string const* lMcSuffix = &fHistoSuffixes[iMcKind];

        // for track and collision histos we only plot reconstructed quantities at the moment
        if (iMcKind == kRec) {

          // no cuts on collisions so far so only plot before cuts
          if (iBARecCuts == kBeforeRecCuts) {
            // collision histograms
            fMyRegistry.mCollision.mBeforeAfterRecCuts[iBARecCuts].mV0Kind[iMcKind].addHistosToOfficalRegistry(fHistogramRegistry,
                                                                                                               lCollisionHistoDefinitions,
                                                                                                               lMcSuffix);
          }

          // track histograms
          fMyRegistry.mTrack.mBeforeAfterRecCuts[iBARecCuts].mV0Kind[iMcKind].addHistosToOfficalRegistry(fHistogramRegistry,
                                                                                                         lTrackHistoDefinitions,
                                                                                                         lMcSuffix);
        }
        // v0 histograms
        fMyRegistry.mV0.mBeforeAfterRecCuts[iBARecCuts].mV0Kind[iMcKind].addHistosToOfficalRegistry(fHistogramRegistry,
                                                                                                    (iMcKind < 3) ? lV0HistoDefinitions : lV0ResolutionHistoDefinitions,
                                                                                                    lMcSuffix);
        // v0 recalculated conversion point histos
        if ((iMcKind == kRec) || (iMcKind == kMCVal)) {
          fMyRegistry.mV0.mBeforeAfterRecCuts[iBARecCuts].mV0Kind[iMcKind].addHistosToOfficalRegistry(fHistogramRegistry,
                                                                                                      lV0HistoDefinitions_recalculated,
                                                                                                      lMcSuffix);
        }

        if (doprocessMc) {
          // v0 mc rejection histos
          for (size_t iRejReason = 0; iRejReason < 2; ++iRejReason) {
            fMyRegistry.mV0.mRejectedByMc[iRejReason].mBeforeAfterRecCuts[iBARecCuts].mV0Kind[iMcKind].addHistosToOfficalRegistry(fHistogramRegistry,
                                                                                                                                  (iMcKind < 3) ? lV0HistoDefinitions : lV0ResolutionHistoDefinitions,
                                                                                                                                  lMcSuffix);
            if ((iMcKind == kRec) || (iMcKind == kMCVal)) {
              fMyRegistry.mV0.mRejectedByMc[iRejReason].mBeforeAfterRecCuts[iBARecCuts].mV0Kind[iMcKind].addHistosToOfficalRegistry(fHistogramRegistry,
                                                                                                                                    lV0HistoDefinitions_recalculated,
                                                                                                                                    lMcSuffix);
            }
          }
        }
      }
    }
  }

  void fillV0SelectionHisto(ePhotonCuts theCase)
  {
    fillTH1(fMyRegistry.mV0.mSpecialHistos.mContainer, "hV0Selection", static_cast<int>(theCase));
  }

  void fillV0McValidationHisto(eV0McValidation theCase)
  {
    fillTH1(fMyRegistry.mV0.mSpecialHistos.mContainer, "hV0McValidation", static_cast<int>(theCase));
  }

  template <typename TV0, typename TTRACKS>
  bool processV0(TV0 const& theV0, const float& theV0CosinePA, TTRACKS const& theTwoV0Daughters)
  {
    auto fillReconstructedInfoHistogramsI = [&](int theBeforeAfter) {
      fillReconstructedInfoHistograms(
        theBeforeAfter,
        theV0,
        theTwoV0Daughters,
        theV0CosinePA);
    };

    fillReconstructedInfoHistogramsI(kBeforeRecCuts);

    // apply track cuts and photon cuts
    if (!(v0PassesTrackCuts(theTwoV0Daughters) &&
          v0PassesPhotonCuts(theV0, theV0CosinePA))) {
      return false;
    }

    fillReconstructedInfoHistogramsI(kAfterRecCuts);

    return kTRUE;
  }

  template <typename TV0, typename TMCGAMMA, typename TTRACKS>
  void fillAllV0HistogramsForRejectedByMc(int theRejReason,
                                          int theBefAftRec,
                                          TMCGAMMA const& theMcPhoton,
                                          TV0 const& theV0,
                                          float const& theV0CosinePA,
                                          TTRACKS const& theTwoV0Daughters)
  {
    fillV0Histograms(
      fMyRegistry.mV0.mRejectedByMc[theRejReason].mBeforeAfterRecCuts[theBefAftRec].mV0Kind[kRec].mContainer,
      theV0,
      theTwoV0Daughters,
      theV0CosinePA);

    fillTruePhotonHistogramsForRejectedByMc(theRejReason,
                                            theBefAftRec,
                                            theMcPhoton,
                                            theV0,
                                            theV0CosinePA,
                                            theTwoV0Daughters);
  }

  template <typename TV0, typename TMCGAMMA, typename TTRACKS>
  bool v0IsGoodValidatedMcPhoton(TMCGAMMA const& theMcPhoton, TV0 const& theV0, float const& theV0CosinePA, bool theV0PassesRecCuts, TTRACKS const& theTwoV0Daughters)
  {
    auto fillRejectedV0HistosI = [&](eMcRejectedSaved theRejReason) {
      fillAllV0HistogramsForRejectedByMc(static_cast<int>(theRejReason),
                                         kBeforeRecCuts,
                                         theMcPhoton,
                                         theV0,
                                         theV0CosinePA,
                                         theTwoV0Daughters);
      if (theV0PassesRecCuts) {
        fillAllV0HistogramsForRejectedByMc(static_cast<int>(theRejReason),
                                           kAfterRecCuts,
                                           theMcPhoton,
                                           theV0,
                                           theV0CosinePA,
                                           theTwoV0Daughters);
      }
    };

    fillV0McValidationHisto(eV0McValidation::kMcMotherIn);

    if (!theMcPhoton.isPhysicalPrimary()) {
      fillV0McValidationHisto(eV0McValidation::kNoPhysicalPrimary);
      fillRejectedV0HistosI(eMcRejectedSaved::kNoPhysicalPrimary);
      return false;
    }

    if (theMcPhoton.pdgCode() != 22) {
      fillV0McValidationHisto(eV0McValidation::kNoPhoton);
      // add here a fillRejectedV0HistosI() call if interested in properties of non-gamma V0s
      return false;
    }

    if (std::abs(theMcPhoton.eta()) > fTruePhotonEtaMax) {
      fillV0McValidationHisto(eV0McValidation::kOutsideMCEtaAcc);
      fillRejectedV0HistosI(eMcRejectedSaved::kOutsideMCEtaAcc);
      return false;
    }

    fillV0McValidationHisto(eV0McValidation::kMcValidatedPhotonOut);
    return true;
  }

  template <typename TV0, typename TMCGAMMATABLE, typename TTRACKS>
  void processMcPhoton(TMCGAMMATABLE const& theMcPhotonForThisV0AsTable,
                       TV0 const& theV0,
                       float const& theV0CosinePA,
                       bool theV0PassesRecCuts,
                       TTRACKS const& theTwoV0Daughters)
  {
    fillV0McValidationHisto(eV0McValidation::kV0in);

    // is a table that might be empty
    if (theMcPhotonForThisV0AsTable.begin() == theMcPhotonForThisV0AsTable.end()) {
      fillV0McValidationHisto(eV0McValidation::kFakeV0);
      return;
    }
    auto const lMcPhoton = theMcPhotonForThisV0AsTable.begin();

    if (!v0IsGoodValidatedMcPhoton(lMcPhoton,
                                   theV0,
                                   theV0CosinePA,
                                   theV0PassesRecCuts,
                                   theTwoV0Daughters)) {
      return;
    }

    fillTruePhotonHistograms(kBeforeRecCuts,
                             lMcPhoton,
                             theV0,
                             theV0CosinePA,
                             theTwoV0Daughters);

    if (theV0PassesRecCuts) {
      fillV0McValidationHisto(eV0McValidation::kMcValAfterRecCuts);
      fillTruePhotonHistograms(kAfterRecCuts,
                               lMcPhoton,
                               theV0,
                               theV0CosinePA,
                               theTwoV0Daughters);
    }
  }

  template <typename TV0, typename TMCGAMMA>
  void fillV0ResolutionHistograms(mapStringHistPtr& theContainer,
                                  TMCGAMMA const& theMcPhoton,
                                  TV0 const& theV0)
  {
    TVector3 lConvPointTrue(theMcPhoton.conversionX(), theMcPhoton.conversionY(), theMcPhoton.conversionZ());
    TVector3 lConvPointRecalc(theV0.recalculatedVtxX(), theV0.recalculatedVtxY(), theV0.recalculatedVtxZ());

    fillTH1(theContainer, "hPtRes", theV0.pt() - theMcPhoton.pt());
    fillTH1(theContainer, "hEtaRes", theV0.eta() - theMcPhoton.eta());
    fillTH1(theContainer, "hPhiRes", theV0.phi() - theMcPhoton.phi());
    fillTH2(theContainer, "hConvPointXYResVsXY", lConvPointTrue.Perp(), theV0.v0radius() - lConvPointTrue.Perp());
    fillTH2(theContainer, "hConvPointXYResVsXY_recalc", lConvPointTrue.Perp(), lConvPointRecalc.Perp() - lConvPointTrue.Perp());
    fillTH2(theContainer, "hConvPointZResVsZ", theMcPhoton.conversionZ(), theV0.z() - theMcPhoton.conversionZ());
    fillTH2(theContainer, "hConvPointZResVsZ_recalc", theMcPhoton.conversionZ(), theV0.recalculatedVtxZ() - theMcPhoton.conversionZ());
  }

  template <typename TV0, typename TMCGAMMA, typename TTRACKS>
  void fillTruePhotonHistograms(int theBefAftRec, TMCGAMMA const& theMcPhoton, TV0 const& theV0, float const& theV0CosinePA, TTRACKS const& theTwoV0Daughters)
  {
    fillV0HistogramsMcGamma(
      fMyRegistry.mV0.mBeforeAfterRecCuts[theBefAftRec].mV0Kind[kMCTrue].mContainer,
      theMcPhoton);

    fillV0Histograms(
      fMyRegistry.mV0.mBeforeAfterRecCuts[theBefAftRec].mV0Kind[kMCVal].mContainer,
      theV0,
      theTwoV0Daughters,
      theV0CosinePA);

    fillV0Histograms_recalculated(
      fMyRegistry.mV0.mBeforeAfterRecCuts[theBefAftRec].mV0Kind[kMCVal].mContainer,
      theV0);

    fillV0ResolutionHistograms(
      fMyRegistry.mV0.mBeforeAfterRecCuts[theBefAftRec].mV0Kind[kRes].mContainer,
      theMcPhoton,
      theV0);
  }

  template <typename TV0, typename TMCGAMMA, typename TTRACKS>
  void fillTruePhotonHistogramsForRejectedByMc(int theRejReason,
                                               int theBefAftRec,
                                               TMCGAMMA const& theMcPhoton,
                                               TV0 const& theV0,
                                               float const& theV0CosinePA,
                                               TTRACKS const& theTwoV0Daughters)
  {
    fillV0HistogramsMcGamma(
      fMyRegistry.mV0.mRejectedByMc[theRejReason].mBeforeAfterRecCuts[theBefAftRec].mV0Kind[kMCTrue].mContainer,
      theMcPhoton);

    fillV0Histograms(
      fMyRegistry.mV0.mRejectedByMc[theRejReason].mBeforeAfterRecCuts[theBefAftRec].mV0Kind[kMCVal].mContainer,
      theV0,
      theTwoV0Daughters,
      theV0CosinePA);

    fillV0Histograms_recalculated(
      fMyRegistry.mV0.mRejectedByMc[theRejReason].mBeforeAfterRecCuts[theBefAftRec].mV0Kind[kMCVal].mContainer,
      theV0);

    fillV0ResolutionHistograms(
      fMyRegistry.mV0.mRejectedByMc[theRejReason].mBeforeAfterRecCuts[theBefAftRec].mV0Kind[kRes].mContainer,
      theMcPhoton,
      theV0);
  }

  Preslice<aod::V0DaughterTracks> perV0 = aod::v0data::v0Id;

  void processRec(aod::Collisions::iterator const& theCollision,
                  V0DatasAdditional const& theV0s,
                  aod::V0DaughterTracks const& theAllTracks)
  {
    fillTH1(fMyRegistry.mCollision.mBeforeAfterRecCuts[kBeforeRecCuts].mV0Kind[kRec].mContainer,
            "hCollisionZ",
            theCollision.posZ());

    for (auto& lV0 : theV0s) {

      auto lTwoV0Daughters = theAllTracks.sliceBy(perV0, lV0.v0Id());
      float lV0CosinePA = lV0.v0cosPA(theCollision.posX(), theCollision.posY(), theCollision.posZ());

      if (!processV0(lV0, lV0CosinePA, lTwoV0Daughters)) {
        continue;
      }
    }
  }
  PROCESS_SWITCH(GammaConversions, processRec, "process reconstructed info", true);

  Preslice<aod::McGammasTrue> gperV0 = aod::v0data::v0Id;

  void processMc(aod::Collisions::iterator const& theCollision,
                 V0DatasAdditional const& theV0s,
                 aod::V0DaughterTracks const& theAllTracks,
                 aod::McGammasTrue const& theV0sTrue)
  {
    fillTH1(fMyRegistry.mCollision.mBeforeAfterRecCuts[kBeforeRecCuts].mV0Kind[kRec].mContainer,
            "hCollisionZ",
            theCollision.posZ());

    for (auto& lV0 : theV0s) {

      auto lTwoV0Daughters = theAllTracks.sliceBy(perV0, lV0.v0Id());
      float lV0CosinePA = lV0.v0cosPA(theCollision.posX(), theCollision.posY(), theCollision.posZ());

      // check if V0 passes rec cuts and fill beforeRecCuts,afterRecCuts [kRec]
      bool lV0PassesRecCuts = processV0(lV0, lV0CosinePA, lTwoV0Daughters);

      // check if it comes from a true photon (lMcPhotonForThisV0AsTable is a table that might be empty)
      auto lMcPhotonForThisV0AsTable = theV0sTrue.sliceBy(gperV0, lV0.v0Id());
      processMcPhoton(lMcPhotonForThisV0AsTable,
                      lV0,
                      lV0CosinePA,
                      lV0PassesRecCuts,
                      lTwoV0Daughters);
    }
  }
  PROCESS_SWITCH(GammaConversions, processMc, "process reconstructed info and mc", false);

  template <typename T>
  std::shared_ptr<T> getTH(mapStringHistPtr const& theMap, std::string const& theName)
  {
    auto lPairIt = theMap.find(theName);
    if (lPairIt == theMap.end()) {
      LOGF(warning, "SFS getTH(): No element with key %s in map.", theName.data());
      return std::shared_ptr<T>(); // SFS todo verify this is correct
    }
    if (!std::holds_alternative<std::shared_ptr<T>>(lPairIt->second)) {
      LOGF(warning, "SFS getTH(): No shared_ptr<T> in map for key %s.", theName.data());
      return std::shared_ptr<T>();
    }
    std::shared_ptr<T> lHisto = std::get<std::shared_ptr<T>>(lPairIt->second);
    if (lHisto == nullptr) {
      LOGF(warning, "SFS getTH(): Found shared_ptr<T> for key %s but it is a nullptr.", theName.data());
      return std::shared_ptr<T>();
    }
    return lHisto;
  }

  void fillTH1(mapStringHistPtr const& theMap, std::string const& theName, float theValue)
  {
    std::shared_ptr<TH1> lHisto = getTH<TH1>(theMap, theName);
    if (lHisto != nullptr) {
      lHisto->Fill(theValue);
    }
  }

  void fillTH2(mapStringHistPtr const& theMap, std::string const& theName, float theValueX, float theValueY)
  {
    std::shared_ptr<TH2> lHisto = getTH<TH2>(theMap, theName);
    if (lHisto != nullptr) {
      lHisto->Fill(theValueX, theValueY);
    }
  }

  template <typename TTRACKS>
  void fillTrackHistograms(mapStringHistPtr& theContainer, TTRACKS const& theTwoV0Daughters)
  {
    auto fillTrackHistogramsI = [&](auto const& theTrack) {
      fillTH1(theContainer, "hTrackEta", theTrack.eta());
      fillTH1(theContainer, "hTrackPhi", theTrack.phi());
      fillTH1(theContainer, "hTrackPt", theTrack.pt());
      fillTH1(theContainer, "hTPCFoundOverFindableCls", theTrack.tpcFoundOverFindableCls());
      fillTH1(theContainer, "hTPCCrossedRowsOverFindableCls", theTrack.tpcCrossedRowsOverFindableCls());
      fillTH2(theContainer, "hTPCdEdxSigEl", theTrack.p(), theTrack.tpcNSigmaEl());
      fillTH2(theContainer, "hTPCdEdxSigPi", theTrack.p(), theTrack.tpcNSigmaPi());
      fillTH2(theContainer, "hTPCdEdx", theTrack.p(), theTrack.tpcSignal());
    };

    for (auto& lTrack : theTwoV0Daughters) {
      fillTrackHistogramsI(lTrack);
    }
  }

  //this is really ugly but it works. Would be easier if .end() worked..
  template <typename TTRACKS>
  void fillhpeDivpGammaHistograms(mapStringHistPtr& theContainer, TTRACKS const& theTwoV0Daughters)
  {
    float p_electron[2];
    Int_t counter = 0;

    for (auto& lTrack : theTwoV0Daughters) {
      p_electron[counter] = lTrack.p();
      counter++;
    }

    fillTH1(theContainer, "hpeDivpGamma", p_electron[0] / (p_electron[0] + p_electron[1]));
    fillTH1(theContainer, "hpeDivpGamma", p_electron[1] / (p_electron[0] + p_electron[1]));
  }

  template <typename TV0, typename TTRACKS>
  void fillV0Histograms(mapStringHistPtr& theContainer, TV0 const& theV0, TTRACKS const& theTwoV0Daughters, float const& theV0CosinePA)
  {
    fillTH1(theContainer, "hEta", theV0.eta());
    fillTH1(theContainer, "hPhi", theV0.phi());
    fillTH1(theContainer, "hPt", theV0.pt());
    fillTH1(theContainer, "hConvPointR", theV0.v0radius());
    fillTH1(theContainer, "hConvPointZ", theV0.z());
    fillTH1(theContainer, "hCosPAngle", theV0CosinePA);
    fillTH2(theContainer, "hArmenteros", theV0.alpha(), theV0.qtarm());
    fillTH2(theContainer, "hPsiPt", theV0.psipair(), theV0.pt());
    fillTH2(theContainer, "hRVsZ", theV0.recalculatedVtxR(), theV0.z()); // as long as z recalculation is not fixed use this

    float p_electron[2];
    Int_t counter = 0;

    // Fill momentum distribution histogram
    for (auto& lTrack : theTwoV0Daughters) {
      p_electron[counter] = lTrack.p();
      counter++;
    }
    fillTH1(theContainer, "hpeDivpGamma", p_electron[0] / theV0.pt());
    fillTH1(theContainer, "hpeDivpGamma", p_electron[1] / theV0.pt());
  }

  template <typename TV0>
  void fillV0Histograms_recalculated(mapStringHistPtr& theContainer, TV0 const& theV0)
  {
    fillTH1(theContainer, "hConvPointR_recalc", theV0.recalculatedVtxR());
    fillTH1(theContainer, "hConvPointZ_recalc", theV0.recalculatedVtxZ());
  }

  // SFS todo: combine fillV0Histograms and fillV0HistogramsMcGamma
  template <typename TMCGAMMA>
  void fillV0HistogramsMcGamma(mapStringHistPtr& theContainer, TMCGAMMA const& theMcGamma)
  {
    fillTH1(theContainer, "hEta", theMcGamma.eta());
    fillTH1(theContainer, "hPhi", theMcGamma.phi());
    fillTH1(theContainer, "hPt", theMcGamma.pt());
    fillTH1(theContainer, "hConvPointR", theMcGamma.v0Radius());
    fillTH1(theContainer, "hConvPointZ", theMcGamma.conversionZ());
    fillTH2(theContainer, "hRVsZ", theMcGamma.v0Radius(), theMcGamma.conversionZ());
  }

  template <typename T>
  bool trackPassesCuts(const T& theTrack)
  {
    // single track eta cut
    if (TMath::Abs(theTrack.eta()) > fTrackEtaMax) {
      fillV0SelectionHisto(ePhotonCuts::kTrackEta);
      return kFALSE;
    }

    // single track pt cut
    if (theTrack.pt() < fTrackPtMin) {
      fillV0SelectionHisto(ePhotonCuts::kTrackPt);
      return kFALSE;
    }

    if (!(selectionPIDTPC_track(theTrack))) {
      return kFALSE;
    }

    if (theTrack.tpcFoundOverFindableCls() < fMinTPCFoundOverFindableCls) {
      fillV0SelectionHisto(ePhotonCuts::kTPCFoundOverFindableCls);
      return kFALSE;
    }

    if (theTrack.tpcCrossedRowsOverFindableCls() < fMinTPCCrossedRowsOverFindableCls) {
      fillV0SelectionHisto(ePhotonCuts::kTPCCrossedRowsOverFindableCls);
      return kFALSE;
    }
    return kTRUE;
  }

  template <typename TTRACKS>
  bool v0PassesTrackCuts(TTRACKS const& theTwoV0Daughters)
  {
    for (auto& lTrack : theTwoV0Daughters) {
      if (!trackPassesCuts(lTrack)) {
        return kFALSE;
      }
    }
    return true;
  }

  template <typename TV0>
  bool v0PassesPhotonCuts(const TV0& theV0, float theV0CosinePA)
  {
    if (theV0.v0radius() < fV0RMin || theV0.v0radius() > fV0RMax) {
      fillV0SelectionHisto(ePhotonCuts::kV0Radius);
      return kFALSE;
    }

    if (fV0PhotonAsymmetryMax > 0. && !ArmenterosQtCut(theV0.alpha(), theV0.qtarm(), theV0.pt())) {
      fillV0SelectionHisto(ePhotonCuts::kArmenteros);
      return kFALSE;
    }

    if (fV0PsiPairMax > 0. && TMath::Abs(theV0.psipair()) > fV0PsiPairMax) {
      fillV0SelectionHisto(ePhotonCuts::kPsiPair);
      return kFALSE;
    }

    Configurable<float> LineCutZ0{"fLineCutZ0", 7.0, "The offset for the linecute used in the Z vs R plot"};
    Configurable<float> LineCutZRSlope{"LineCutZRSlope", (float)TMath::Tan(2 * TMath::ATan(TMath::Exp(-fTruePhotonEtaMax))), "The slope for the line cut"};

    if (fV0CosPAngleMin > 0. && theV0CosinePA < fV0CosPAngleMin) {
      fillV0SelectionHisto(ePhotonCuts::kCosinePA);
      return kFALSE;
    }

    if (TMath::Abs(theV0.z()) > LineCutZ0 + theV0.recalculatedVtxR() * LineCutZRSlope) { // as long as z recalculation is not fixed use this
      fillV0SelectionHisto(ePhotonCuts::kRZLine);
      return kFALSE;
    }
    return kTRUE;
  }

  template <typename TV0, typename TTRACKS>
  void fillReconstructedInfoHistograms(int theBefAftRec, TV0 const& theV0, TTRACKS const& theTwoV0Daughters, float const& theV0CosinePA)
  {
    fillV0SelectionHisto(!theBefAftRec ? ePhotonCuts::kV0In : ePhotonCuts::kV0Out);

    fillTrackHistograms(
      fMyRegistry.mTrack.mBeforeAfterRecCuts[theBefAftRec].mV0Kind[kRec].mContainer,
      theTwoV0Daughters);

    fillV0Histograms(
      fMyRegistry.mV0.mBeforeAfterRecCuts[theBefAftRec].mV0Kind[kRec].mContainer,
      theV0,
      theTwoV0Daughters,
      theV0CosinePA);

    fillV0Histograms_recalculated(
      fMyRegistry.mV0.mBeforeAfterRecCuts[theBefAftRec].mV0Kind[kRec].mContainer,
      theV0);
  }

  Bool_t ArmenterosQtCut(Double_t theAlpha, Double_t theQt, Double_t thePt)
  {
    // in AliPhysics this is the cut for if fDo2DQt && fDoQtGammaSelection == 2
    Float_t lQtMaxPtDep = fV0QtPtMultiplicator * thePt;
    if (lQtMaxPtDep > fV0QtMax) {
      lQtMaxPtDep = fV0QtMax;
    }
    if ((TMath::Power(theAlpha / fV0PhotonAsymmetryMax, 2) + TMath::Power(theQt / lQtMaxPtDep, 2)) >= 1) {
      return kFALSE;
    }
    return kTRUE;
  }

  template <typename T>
  bool selectionPIDTPC_track(const T& theTrack)
  {
    // TPC Electron Line
    if (fPIDnSigmaElectronMin && (theTrack.tpcNSigmaEl() < fPIDnSigmaElectronMin || theTrack.tpcNSigmaEl() > fPIDnSigmaElectronMax)) {
      fillV0SelectionHisto(ePhotonCuts::kElectronPID);
      return kFALSE;
    }

    // TPC Pion Line
    if (theTrack.p() > fPIDPionRejectionPMin) {
      // low pt Pion rej
      if (theTrack.p() < fPIDPionRejectionPBoarder) {
        if (theTrack.tpcNSigmaEl() > fPIDnSigmaElectronMin && theTrack.tpcNSigmaEl() < fPIDnSigmaElectronMax && theTrack.tpcNSigmaPi() < fPIDnSigmaAbovePionLineLowPMin) {
          fillV0SelectionHisto(ePhotonCuts::kPionRejLowMom);
          return kFALSE;
        }
      }
      // High Pt Pion rej
      else {
        if (theTrack.tpcNSigmaEl() > fPIDnSigmaElectronMin && theTrack.tpcNSigmaEl() < fPIDnSigmaElectronMax && theTrack.tpcNSigmaPi() < fPIDnSigmaAbovePionLineHighPMin) {
          fillV0SelectionHisto(ePhotonCuts::kPionRejHighMom);
          return kFALSE;
        }
      }
    }
    return kTRUE;
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<GammaConversions>(cfgc)};
}