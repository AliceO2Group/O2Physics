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
//
//  *+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
//   Decay finder task for ALICE 3
//  *+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
//
//    Uses specific ALICE 3 PID and performance for studying
//    HF decays. Work in progress: use at your own risk!
//

#include "PWGLF/DataModel/LFParticleIdentification.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "ALICE3/DataModel/A3DecayFinderTables.h"
#include "ALICE3/DataModel/OTFMulticharm.h"
#include "ALICE3/DataModel/OTFRICH.h"
#include "ALICE3/DataModel/OTFStrangeness.h"
#include "ALICE3/DataModel/OTFTOF.h"
#include "ALICE3/DataModel/tracksAlice3.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "CommonConstants/PhysicsConstants.h"
#include "DCAFitter/DCAFitterN.h"
#include "DataFormatsCalibration/MeanVertexObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsVertexing/PVertexer.h"
#include "DetectorsVertexing/PVertexerHelpers.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

#include <array>
#include <cmath>
#include <cstdlib>
#include <iterator>
#include <map>
#include <utility>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using multiCharmTracksPID = soa::Join<aod::MCharmCores, aod::MCharmPID>;
using multiCharmTracksFull = soa::Join<aod::MCharmCores, aod::MCharmPID, aod::MCharmExtra>;

struct alice3multicharm {
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  std::map<int, int> pdgToBin;

  ConfigurableAxis axisEta{"axisEta", {80, -4.0f, +4.0f}, "#eta"};
  ConfigurableAxis axisXiccMass{"axisXiccMass", {200, 3.521f, 3.721f}, "Xicc Inv Mass (GeV/c^{2})"};
  ConfigurableAxis axisDCA{"axisDCA", {400, 0, 400}, "DCA (#mum)"};
  ConfigurableAxis axisRadiusLarge{"axisRadiusLarge", {1000, 0, 20}, "Decay radius (cm)"};
  ConfigurableAxis axisRadius{"axisRadius", {10000, 0, 10000}, "Decay radius (#mum)"};
  ConfigurableAxis axisTofTrackDelta{"axisTofTrackDelta", {200, 0, 1000}, "TOF track time"};
  ConfigurableAxis axisNSigma{"axisNSigma", {21, -10, 10}, "nsigma"};
  ConfigurableAxis axisDecayLength{"axisDecayLength", {2000, 0, 2000}, "Decay lenght (#mum)"};
  ConfigurableAxis axisDcaDaughters{"axisDcaDaughters", {200, 0, 100}, "DCA (mum)"};
  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 17.0f, 19.0f, 21.0f, 23.0f, 25.0f, 30.0f, 35.0f, 40.0f, 50.0f}, "pt axis for QA histograms"};

  Configurable<float> xiMinDCAxy{"xiMinDCAxy", -1, "[0] in |DCAxy| > [0]+[1]/pT"};
  Configurable<float> xiMinDCAz{"xiMinDCAz", -1, "[0] in |DCAz| > [0]+[1]/pT"};
  Configurable<float> xiMinRadius{"xiMinRadius", -1, "Minimum R2D for Xic decay (cm)"};

  Configurable<float> picMinDCAxy{"picMinDCAxy", -1, "[0] in |DCAz| > [0]+[1]/pT"};
  Configurable<float> picMinDCAz{"picMinDCAz", -1, "[0] in |DCAxy| > [0]+[1]/pT"};
  Configurable<float> picMinPt{"picMinPt", -1, "Minimum pT for Xic pions"};

  Configurable<float> piccMinDCAxy{"piccMinDCAxy", -1, "[0] in |DCAxy| > [0]+[1]/pT"};
  Configurable<float> piccMinDCAz{"piccMinDCAz", -1, "[0] in |DCAz| > [0]+[1]/pT"};
  Configurable<float> piccMinPt{"piccMinPt", -1, "Minimum pT for Xicc pions"};

  Configurable<float> xicMaxDauDCA{"xicMaxDauDCA", 1e+4, "DCA between Xic daughters (cm)"};
  Configurable<float> xicMinDCAxy{"xicMinDCAxy", -1, "[0] in |DCAz| > [0]+[1]/pT"};
  Configurable<float> xicMinDCAz{"xicMinDCAz", -1, "[0] in |DCAxy| > [0]+[1]/pT"};
  Configurable<float> xiccMaxDCAxy{"xiccMaxDCAxy", 1e+4, "Maximum DCAxy"};
  Configurable<float> xiccMaxDCAz{"xiccMaxDCAz", 1e+4, "Maximum DCAz"};
  Configurable<float> xicMinRadius{"xicMinRadius", -1, "Minimum R2D for Xic decay (cm)"};
  Configurable<float> xicMinDecayDistanceFromPV{"xicMinDecayDistanceFromPV", -1, "Minimum distance for Xic decay from PV (cm)"};
  Configurable<float> xicMinProperLength{"xicMinProperLength", -1, "Minimum proper length for Xic decay (cm)"};
  Configurable<float> xicMaxProperLength{"xicMaxProperLength", 1e+4, "Minimum proper length for Xic decay (cm)"};

  Configurable<float> xiccMaxDauDCA{"xiccMaxDauDCA", 1e+4, "DCA between Xicc daughters (cm)"};
  Configurable<float> xiccMinRadius{"xiccMinRadius", -1, "Minimum R2D for Xicc decay (cm)"};
  Configurable<float> xiccMinProperLength{"xiccMinProperLength", -1, "Minimum proper length for Xicc decay (cm)"};
  Configurable<float> xiccMaxProperLength{"xiccMaxProperLength", 1e+4, "Minimum proper length for Xicc decay (cm)"};

  void init(InitContext&)
  {
    histos.add("SelectionQA/hDCAXicDaughters", "hDCAXicDaughters; DCA between Xic daughters (#mum)", kTH1D, {axisDcaDaughters});
    histos.add("SelectionQA/hDCAXiccDaughters", "hDCAXiccDaughters; DCA between Xicc daughters (#mum)", kTH1D, {axisDcaDaughters});
    histos.add("SelectionQA/hDCAxyXi", "hDCAxyXi; Xi DCAxy to PV (#mum)", kTH1D, {axisDCA});
    histos.add("SelectionQA/hDCAzXi", "hDCAzXi; Xi DCAz to PV (#mum)", kTH1D, {axisDCA});
    histos.add("SelectionQA/hDCAxyXic", "hDCAxyXic; Xic DCAxy to PV (#mum)", kTH1D, {axisDCA});
    histos.add("SelectionQA/hDCAzXic", "hDCAzXic; Xic DCAz to PV (#mum)", kTH1D, {axisDCA});
    histos.add("SelectionQA/hDCAxyXicc", "hDCAxyXicc; Xicc DCAxy to PV (#mum)", kTH1D, {axisDCA});
    histos.add("SelectionQA/hDCAzXicc", "hDCAzXicc; Xicc DCAz to PV (#mum)", kTH1D, {axisDCA});
    histos.add("SelectionQA/hDecayRadiusXic", "hDecayRadiusXic; Distance (#mum)", kTH1D, {axisRadius});
    histos.add("SelectionQA/hDecayRadiusXicc", "hDecayRadiusXicc; Distance (#mum)", kTH1D, {axisRadius});
    histos.add("SelectionQA/hDecayDistanceFromPVXic", "hDecayDistanceFromPVXic; Distance (#mum)", kTH1D, {axisDecayLength});
    histos.add("SelectionQA/hProperLengthXic", "hProperLengthXic; Distance (#mum)", kTH1D, {axisDecayLength});
    histos.add("SelectionQA/hProperLengthXicc", "hProperLengthXicc; Distance (#mum)", kTH1D, {axisDecayLength});
    histos.add("SelectionQA/hPi1cPt", "hPi1cPt; Pi1c pT (Gev/#it(c))", kTH1D, {axisPt});
    histos.add("SelectionQA/hPi2cPt", "hPi2cPt; Pi2c pT (Gev/#it(c))", kTH1D, {axisPt});
    histos.add("SelectionQA/hPiccPt", "hPiccPt; Picc pT (Gev/#it(c))", kTH1D, {axisPt});

    auto hMCharmBuilding = histos.add<TH1>("hMCharmBuilding", "hMCharmBuilding", kTH1D, {{22, -0.5, 21.5}});
    hMCharmBuilding->GetXaxis()->SetBinLabel(1, "nTotalCandidates");
    hMCharmBuilding->GetXaxis()->SetBinLabel(2, "xicMaxDauDCA");
    hMCharmBuilding->GetXaxis()->SetBinLabel(3, "xiccMaxDauDCA");
    hMCharmBuilding->GetXaxis()->SetBinLabel(4, "xiMinDCAxy");
    hMCharmBuilding->GetXaxis()->SetBinLabel(5, "xiMinDCAz");
    hMCharmBuilding->GetXaxis()->SetBinLabel(6, "picMinDCAxy");
    hMCharmBuilding->GetXaxis()->SetBinLabel(7, "picMinDCAz");
    hMCharmBuilding->GetXaxis()->SetBinLabel(8, "picMinDCAxy");
    hMCharmBuilding->GetXaxis()->SetBinLabel(9, "picMinDCAz");
    hMCharmBuilding->GetXaxis()->SetBinLabel(10, "piccMinDCAxy");
    hMCharmBuilding->GetXaxis()->SetBinLabel(11, "piccMinDCAz");
    hMCharmBuilding->GetXaxis()->SetBinLabel(12, "xicMinDCAxy");
    hMCharmBuilding->GetXaxis()->SetBinLabel(13, "xicMinDCAz");
    hMCharmBuilding->GetXaxis()->SetBinLabel(14, "xiccMaxDCAxy");
    hMCharmBuilding->GetXaxis()->SetBinLabel(15, "xiccMaxDCAz");
    hMCharmBuilding->GetXaxis()->SetBinLabel(16, "xicMinRadius");
    hMCharmBuilding->GetXaxis()->SetBinLabel(17, "xiccMinRadius");
    hMCharmBuilding->GetXaxis()->SetBinLabel(18, "xicMinProperLength");
    hMCharmBuilding->GetXaxis()->SetBinLabel(19, "xicMaxProperLength");
    hMCharmBuilding->GetXaxis()->SetBinLabel(20, "xiccMinProperLength");
    hMCharmBuilding->GetXaxis()->SetBinLabel(21, "xiccMaxProperLength");
    hMCharmBuilding->GetXaxis()->SetBinLabel(22, "xicMinDecayDistanceFromPV");

    if (doprocessXiccPID || doprocessXiccExtra) {
      auto hPdgCodes = histos.add<TH2>("PIDQA/hPdgCodes", "hPdgCodes", kTH2D, {{3, 0.5, 3.5}, {7, 0.5, 7.5}});
      hPdgCodes->GetXaxis()->SetBinLabel(1, "pi1c");
      hPdgCodes->GetXaxis()->SetBinLabel(2, "pi2c");
      hPdgCodes->GetXaxis()->SetBinLabel(3, "picc");
      hPdgCodes->GetYaxis()->SetBinLabel(1, "el");
      hPdgCodes->GetYaxis()->SetBinLabel(2, "mu");
      hPdgCodes->GetYaxis()->SetBinLabel(3, "pi");
      hPdgCodes->GetYaxis()->SetBinLabel(4, "ka");
      hPdgCodes->GetYaxis()->SetBinLabel(5, "pr");
      hPdgCodes->GetYaxis()->SetBinLabel(6, "xi");
      hPdgCodes->GetYaxis()->SetBinLabel(7, "other");
      pdgToBin.insert({kElectron, 1});
      pdgToBin.insert({kMuonMinus, 2});
      pdgToBin.insert({kPiPlus, 3});
      pdgToBin.insert({kKPlus, 4});
      pdgToBin.insert({kProton, 5});
      pdgToBin.insert({kXiMinus, 6});

      histos.add("PIDQA/hInnerTofTimeDeltaPi1c", "hInnerTofTimeDeltaPi1c; Reco - expected pion (ps)", kTH1D, {axisTofTrackDelta});
      histos.add("PIDQA/hInnerTofTimeDeltaPi2c", "hInnerTofTimeDeltaPi2c; Reco - expected pion (ps)", kTH1D, {axisTofTrackDelta});
      histos.add("PIDQA/hInnerTofTimeDeltaPicc", "hInnerTofTimeDeltaPicc; Reco - expected pion (ps)", kTH1D, {axisTofTrackDelta});
      histos.add("PIDQA/hOuterTofTimeDeltaPi1c", "hOuterTofTimeDeltaPi1c; Reco - expected pion (ps)", kTH1D, {axisTofTrackDelta});
      histos.add("PIDQA/hOuterTofTimeDeltaPi2c", "hOuterTofTimeDeltaPi2c; Reco - expected pion (ps)", kTH1D, {axisTofTrackDelta});
      histos.add("PIDQA/hOuterTofTimeDeltaPicc", "hOuterTofTimeDeltaPicc; Reco - expected pion (ps)", kTH1D, {axisTofTrackDelta});

      histos.add("PIDQA/hInnerTofNSigmaPi1c", "hInnerTofNSigmaPi1c; TOF NSigma pion", kTH2D, {axisPt, axisNSigma});
      histos.add("PIDQA/hOuterTofNSigmaPi1c", "hOuterTofNSigmaPi1c; TOF NSigma pion", kTH2D, {axisPt, axisNSigma});
      histos.add("PIDQA/hInnerTofNSigmaPi2c", "hInnerTofNSigmaPi2c; TOF NSigma pion", kTH2D, {axisPt, axisNSigma});
      histos.add("PIDQA/hOuterTofNSigmaPi2c", "hOuterTofNSigmaPi2c; TOF NSigma pion", kTH2D, {axisPt, axisNSigma});
      histos.add("PIDQA/hInnerTofNSigmaPicc", "hInnerTofNSigmaPicc; TOF NSigma pion", kTH2D, {axisPt, axisNSigma});
      histos.add("PIDQA/hOuterTofNSigmaPicc", "hOuterTofNSigmaPicc; TOF NSigma pion", kTH2D, {axisPt, axisNSigma});
      histos.add("PIDQA/hRichNSigmaPi1c", "hRichNSigmaPi1c; RICH NSigma pion", kTH2D, {axisPt, axisNSigma});
      histos.add("PIDQA/hRichNSigmaPi2c", "hRichNSigmaPi2c; RICH NSigma pion", kTH2D, {axisPt, axisNSigma});
      histos.add("PIDQA/hRichNSigmaPicc", "hRichNSigmaPicc; RICH NSigma pion", kTH2D, {axisPt, axisNSigma});
    }

    if (doprocessXiccExtra) {
      histos.add("XiccProngs/h3dPos", "h3dPos; Xicc pT (GeV/#it(c)); Pos pT (GeV/#it(c)); Pos #eta", kTH3D, {axisPt, axisPt, axisEta});
      histos.add("XiccProngs/h3dNeg", "h3dNeg; Xicc pT (GeV/#it(c)); Neg pT (GeV/#it(c)); Neg #eta", kTH3D, {axisPt, axisPt, axisEta});
      histos.add("XiccProngs/h3dBach", "h3dBach; Xicc pT (GeV/#it(c)); Bach pT (GeV/#it(c)); Bach #eta", kTH3D, {axisPt, axisPt, axisEta});
      histos.add("XiccProngs/h3dPi1c", "h3dPi1c; Xicc pT (GeV/#it(c)); Pi1c pT (GeV/#it(c)); Pi1c #eta", kTH3D, {axisPt, axisPt, axisEta});
      histos.add("XiccProngs/h3dPi2c", "h3dPi2c; Xicc pT (GeV/#it(c)); Pi2c pT (GeV/#it(c)); Pi2c #eta", kTH3D, {axisPt, axisPt, axisEta});
      histos.add("XiccProngs/h3dPicc", "h3dPicc; Xicc pT (GeV/#it(c)); Picc pT (GeV/#it(c)); Picc #eta", kTH3D, {axisPt, axisPt, axisEta});
    }
    histos.add("h3dXicc", "h3dXicc; Xicc pT (GeV/#it(c)); Xicc #eta; Xicc mass (GeV/#it(c)^{2})", kTH3D, {axisPt, axisEta, axisXiccMass});
  }

  int getBin(const std::map<int, int>& pdgToBin, int pdg)
  {
    auto it = pdgToBin.find(pdg);
    return (it != pdgToBin.end()) ? it->second : 7;
  }

  template <typename TMCharmCands>
  void genericProcessXicc(TMCharmCands xiccCands)
  {
    for (const auto& xiccCand : xiccCands) {

      histos.fill(HIST("hMCharmBuilding"), 0);
      if (xiccCand.xicDauDCA() > xicMaxDauDCA)
        continue;

      histos.fill(HIST("hMCharmBuilding"), 1);
      if (xiccCand.xiccDauDCA() > xiccMaxDauDCA)
        continue;

      histos.fill(HIST("hMCharmBuilding"), 2);
      if (std::fabs(xiccCand.xiDCAxy()) < xiMinDCAxy)
        continue;

      histos.fill(HIST("hMCharmBuilding"), 3);
      if (std::fabs(xiccCand.xiDCAz()) < xiMinDCAz)
        continue;

      histos.fill(HIST("hMCharmBuilding"), 4);
      if (std::fabs(xiccCand.pi1cDCAxy()) < picMinDCAxy)
        continue;

      histos.fill(HIST("hMCharmBuilding"), 5);
      if (std::fabs(xiccCand.pi1cDCAz()) < picMinDCAz)
        continue;

      histos.fill(HIST("hMCharmBuilding"), 6);
      if (std::fabs(xiccCand.pi2cDCAxy()) < picMinDCAxy)
        continue;

      histos.fill(HIST("hMCharmBuilding"), 7);
      if (std::fabs(xiccCand.pi2cDCAz()) < picMinDCAz)
        continue;

      histos.fill(HIST("hMCharmBuilding"), 8);
      if (std::fabs(xiccCand.piccDCAxy()) < piccMinDCAxy)
        continue;

      histos.fill(HIST("hMCharmBuilding"), 9);
      if (std::fabs(xiccCand.piccDCAz()) < piccMinDCAz)
        continue;

      histos.fill(HIST("hMCharmBuilding"), 10);
      if (std::fabs(xiccCand.xicDCAxy()) < xicMinDCAxy)
        continue;

      histos.fill(HIST("hMCharmBuilding"), 11);
      if (std::fabs(xiccCand.xicDCAz()) < xicMinDCAz)
        continue;

      histos.fill(HIST("hMCharmBuilding"), 12);
      if (std::fabs(xiccCand.xiccDCAxy()) > xiccMaxDCAxy)
        continue;

      histos.fill(HIST("hMCharmBuilding"), 13);
      if (std::fabs(xiccCand.xiccDCAz()) > xiccMaxDCAz)
        continue;

      histos.fill(HIST("hMCharmBuilding"), 14);
      if (xiccCand.xicDecayRadius2D() < xicMinRadius)
        continue;

      histos.fill(HIST("hMCharmBuilding"), 15);
      if (xiccCand.xiccDecayRadius2D() < xiccMinRadius)
        continue;

      histos.fill(HIST("hMCharmBuilding"), 16);
      if (xiccCand.xicProperLength() < xicMinProperLength)
        continue;

      histos.fill(HIST("hMCharmBuilding"), 17);
      if (xiccCand.xicProperLength() > xicMaxProperLength)
        continue;

      histos.fill(HIST("hMCharmBuilding"), 18);
      if (xiccCand.xiccProperLength() < xiccMinProperLength)
        continue;

      histos.fill(HIST("hMCharmBuilding"), 19);
      if (xiccCand.xiccProperLength() > xiccMaxProperLength)
        continue;

      histos.fill(HIST("hMCharmBuilding"), 20);
      if (xiccCand.xicDistanceFromPV() < xicMinDecayDistanceFromPV)
        continue;

      histos.fill(HIST("hMCharmBuilding"), 21);
      histos.fill(HIST("SelectionQA/hDCAXicDaughters"), xiccCand.xicDauDCA() * 1e+4);
      histos.fill(HIST("SelectionQA/hDCAXiccDaughters"), xiccCand.xiccDauDCA() * 1e+4);
      histos.fill(HIST("SelectionQA/hDCAxyXi"), std::fabs(xiccCand.xiDCAxy() * 1e+4));
      histos.fill(HIST("SelectionQA/hDCAzXi"), std::fabs(xiccCand.xiDCAz() * 1e+4));
      histos.fill(HIST("SelectionQA/hDCAxyXic"), std::fabs(xiccCand.xicDCAxy() * 1e+4));
      histos.fill(HIST("SelectionQA/hDCAzXic"), std::fabs(xiccCand.xicDCAz() * 1e+4));
      histos.fill(HIST("SelectionQA/hDCAxyXicc"), std::fabs(xiccCand.xiccDCAxy() * 1e+4));
      histos.fill(HIST("SelectionQA/hDCAzXicc"), std::fabs(xiccCand.xiccDCAz() * 1e+4));
      histos.fill(HIST("SelectionQA/hDecayRadiusXic"), xiccCand.xicDecayRadius2D() * 1e+4);
      histos.fill(HIST("SelectionQA/hDecayRadiusXicc"), xiccCand.xiccDecayRadius2D() * 1e+4);
      histos.fill(HIST("SelectionQA/hDecayDistanceFromPVXic"), xiccCand.xicDistanceFromPV() * 1e+4);
      histos.fill(HIST("SelectionQA/hProperLengthXic"), xiccCand.xicProperLength() * 1e+4);
      histos.fill(HIST("SelectionQA/hProperLengthXicc"), xiccCand.xiccProperLength() * 1e+4);
      histos.fill(HIST("SelectionQA/hPi1cPt"), xiccCand.pi1cPt());
      histos.fill(HIST("SelectionQA/hPi2cPt"), xiccCand.pi2cPt());
      histos.fill(HIST("SelectionQA/hPiccPt"), xiccCand.piccPt());

      if constexpr (requires { xiccCand.pi1cTofDeltaInner(); }) { // if pid table
        histos.fill(HIST("PIDQA/hInnerTofTimeDeltaPi1c"), xiccCand.pi1cTofDeltaInner());
        histos.fill(HIST("PIDQA/hInnerTofTimeDeltaPi2c"), xiccCand.pi2cTofDeltaInner());
        histos.fill(HIST("PIDQA/hInnerTofTimeDeltaPicc"), xiccCand.piccTofDeltaInner());
        histos.fill(HIST("PIDQA/hOuterTofTimeDeltaPi1c"), xiccCand.pi1cTofDeltaOuter());
        histos.fill(HIST("PIDQA/hOuterTofTimeDeltaPi2c"), xiccCand.pi2cTofDeltaOuter());
        histos.fill(HIST("PIDQA/hOuterTofTimeDeltaPicc"), xiccCand.piccTofDeltaOuter());
        histos.fill(HIST("PIDQA/hInnerTofNSigmaPi1c"), xiccCand.pi1cPt(), xiccCand.pi1cTofNSigmaInner());
        histos.fill(HIST("PIDQA/hOuterTofNSigmaPi1c"), xiccCand.pi1cPt(), xiccCand.pi1cTofNSigmaOuter());
        histos.fill(HIST("PIDQA/hInnerTofNSigmaPi2c"), xiccCand.pi2cPt(), xiccCand.pi2cTofNSigmaInner());
        histos.fill(HIST("PIDQA/hOuterTofNSigmaPi2c"), xiccCand.pi2cPt(), xiccCand.pi2cTofNSigmaOuter());
        histos.fill(HIST("PIDQA/hInnerTofNSigmaPicc"), xiccCand.piccPt(), xiccCand.piccTofNSigmaInner());
        histos.fill(HIST("PIDQA/hOuterTofNSigmaPicc"), xiccCand.piccPt(), xiccCand.piccTofNSigmaOuter());
        if (xiccCand.pi1cHasRichSignal()) {
          histos.fill(HIST("PIDQA/hRichNSigmaPi1c"), xiccCand.pi1cPt(), xiccCand.pi1cRichNSigma());
        }
        if (xiccCand.pi2cHasRichSignal()) {
          histos.fill(HIST("PIDQA/hRichNSigmaPi2c"), xiccCand.pi2cPt(), xiccCand.pi2cRichNSigma());
        }
        if (xiccCand.piccHasRichSignal()) {
          histos.fill(HIST("PIDQA/hRichNSigmaPicc"), xiccCand.piccPt(), xiccCand.piccRichNSigma());
        }

        histos.fill(HIST("PIDQA/hPdgCodes"), 1, getBin(pdgToBin, std::abs(xiccCand.pi1cPdgCode())));
        histos.fill(HIST("PIDQA/hPdgCodes"), 2, getBin(pdgToBin, std::abs(xiccCand.pi2cPdgCode())));
        histos.fill(HIST("PIDQA/hPdgCodes"), 3, getBin(pdgToBin, std::abs(xiccCand.piccPdgCode())));
      }

      if constexpr (requires { xiccCand.negPt(); }) { // if extra table
        histos.fill(HIST("XiccProngs/h3dNeg"), xiccCand.xiccPt(), xiccCand.negPt(), xiccCand.negEta());
        histos.fill(HIST("XiccProngs/h3dPos"), xiccCand.xiccPt(), xiccCand.posPt(), xiccCand.posEta());
        histos.fill(HIST("XiccProngs/h3dBach"), xiccCand.xiccPt(), xiccCand.bachPt(), xiccCand.bachEta());
        histos.fill(HIST("XiccProngs/h3dPi1c"), xiccCand.xiccPt(), xiccCand.pi1cPt(), xiccCand.pi1cEta());
        histos.fill(HIST("XiccProngs/h3dPi2c"), xiccCand.xiccPt(), xiccCand.pi2cPt(), xiccCand.pi2cEta());
        histos.fill(HIST("XiccProngs/h3dPicc"), xiccCand.xiccPt(), xiccCand.piccPt(), xiccCand.piccEta());
      }

      histos.fill(HIST("h3dXicc"), xiccCand.xiccPt(), xiccCand.xiccEta(), xiccCand.xiccMass());
    }
  }

  void processXicc(aod::MCharmCores const& multiCharmTracks)
  {
    genericProcessXicc(multiCharmTracks);
  }

  void processXiccPID(multiCharmTracksPID const& multiCharmTracks)
  {
    genericProcessXicc(multiCharmTracks);
  }

  void processXiccExtra(multiCharmTracksFull const& multiCharmTracks)
  {
    genericProcessXicc(multiCharmTracks);
  }

  PROCESS_SWITCH(alice3multicharm, processXicc, "find Xicc baryons", true);
  PROCESS_SWITCH(alice3multicharm, processXiccPID, "find Xicc baryons with more QA from PID information", false);
  PROCESS_SWITCH(alice3multicharm, processXiccExtra, "find Xicc baryons with all QA", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<alice3multicharm>(cfgc)};
}
