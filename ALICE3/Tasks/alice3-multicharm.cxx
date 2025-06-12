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

#include <cmath>
#include <array>
#include <cstdlib>
#include <map>
#include <iterator>
#include <utility>

#include "Framework/runDataProcessing.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "DCAFitter/DCAFitterN.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/LFParticleIdentification.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsCalibration/MeanVertexObject.h"
#include "ALICE3/DataModel/OTFTOF.h"
#include "ALICE3/DataModel/RICH.h"
#include "ALICE3/DataModel/A3DecayFinderTables.h"
#include "ALICE3/DataModel/OTFStrangeness.h"
#include "ALICE3/DataModel/OTFMulticharm.h"
#include "ALICE3/DataModel/tracksAlice3.h"
#include "DetectorsVertexing/PVertexer.h"
#include "DetectorsVertexing/PVertexerHelpers.h"
#include "CommonConstants/PhysicsConstants.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// using multicharmtracks = soa::Join<aod::MCharmIndices, aod::MCharmCores>;

struct alice3multicharm {
  SliceCache cache;

  ConfigurableAxis axisEta{"axisEta", {80, -4.0f, +4.0f}, "#eta"};
  ConfigurableAxis axisXiccMass{"axisXiccMass", {200, 3.521f, 3.721f}, "Xicc Inv Mass (GeV/c^{2})"};
  ConfigurableAxis axisDCA{"axisDCA", {400, 0, 400}, "DCA (#mum)"};
  ConfigurableAxis axisRadiusLarge{"axisRadiusLarge", {1000, 0, 20}, "Decay radius (cm)"};
  ConfigurableAxis axisRadius{"axisRadius", {10000, 0, 10000}, "Decay radius (#mum)"};
  ConfigurableAxis axisTofTrackDelta{"axisTofTrackDelta", {1000, 0, 5000}, "TOF track time"};
  ConfigurableAxis axisDecayLength{"axisDecayLength", {2000, 0, 2000}, "Decay lenght (#mum)"};
  ConfigurableAxis axisDcaDaughters{"axisDcaDaughters", {200, 0, 100}, "DCA (mum)"};
  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 17.0f, 19.0f, 21.0f, 23.0f, 25.0f, 30.0f, 35.0f, 40.0f, 50.0f}, "pt axis for QA histograms"};

  Configurable<float> xiMinDCAxy{"xiMinDCAxy", -1, "[0] in |DCAxy| > [0]+[1]/pT"};
  Configurable<float> xiMinDCAz{"xiMinDCAz", -1, "[0] in |DCAz| > [0]+[1]/pT"};
  Configurable<float> xiMinRadius{"xiMinRadius", -1, "Minimum R2D for Xic decay (cm)"};

  Configurable<float> picMinDCAxy{"picMinDCAxy", -1, "[0] in |DCAz| > [0]+[1]/pT"};
  Configurable<float> picMinDCAz{"picMinDCAz", -1, "[0] in |DCAxy| > [0]+[1]/pT"};
  Configurable<float> picMaxTofDiffInner{"picTofDiffInner", 1e+4, "|signal - expected| (ps)"};
  Configurable<float> picMinPt{"picMinPt", -1, "Minimum pT for Xic pions"};

  Configurable<float> piccMinDCAxy{"piccMinDCAxy", -1, "[0] in |DCAxy| > [0]+[1]/pT"};
  Configurable<float> piccMinDCAz{"piccMinDCAz", -1, "[0] in |DCAz| > [0]+[1]/pT"};
  Configurable<float> piccMaxTofDiffInner{"piccMaxTofDiffInner", 1e+4, "|signal - expected| (ps)"};
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

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

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
    histos.add("SelectionQA/hPi1cPt", "hPi1cPt; Pi1c pT (Gev/#it(c))", kTH1D, {axisPt});
    histos.add("SelectionQA/hPi2cPt", "hPi2cPt; Pi2c pT (Gev/#it(c))", kTH1D, {axisPt});
    histos.add("SelectionQA/hPiccPt", "hPiccPt; Picc pT (Gev/#it(c))", kTH1D, {axisPt});
    histos.add("SelectionQA/hXicDecayRadius", "hXicDecayRadius; Distance (#mum)", kTH1D, {axisRadius});
    histos.add("SelectionQA/hXiccDecayRadius", "hXiccDecayRadius; Distance (#mum)", kTH1D, {axisRadius});
    histos.add("SelectionQA/hXicDecayDistanceFromPV", "hXicDecayDistanceFromPV; Distance (#mum)", kTH1D, {axisDecayLength});
    histos.add("SelectionQA/hProperLengthXic", "hProperLengthXic; Distance (#mum)", kTH1D, {axisDecayLength});
    histos.add("SelectionQA/hProperLengthXicc", "hProperLengthXicc; Distance (#mum)", kTH1D, {axisDecayLength});
    histos.add("SelectionQA/hInnerTofTimeDeltaPi1c", "hInnerTofTimeDeltaPi1c; Reco - expected pion (ps)", kTH1D, {axisTofTrackDelta});
    histos.add("SelectionQA/hInnerTofTimeDeltaPi2c", "hInnerTofTimeDeltaPi2c; Reco - expected pion (ps)", kTH1D, {axisTofTrackDelta});
    histos.add("SelectionQA/hInnerTofTimeDeltaPicc", "hInnerTofTimeDeltaPicc; Reco - expected pion (ps)", kTH1D, {axisTofTrackDelta});

    histos.add("XiccProngs/h3dPos", "h3dPos; Xicc pT (GeV/#it(c)); Pos pT (GeV/#it(c)); Pos #eta", kTH3D, {axisPt, axisPt, axisEta});
    histos.add("XiccProngs/h3dNeg", "h3dNeg; Xicc pT (GeV/#it(c)); Neg pT (GeV/#it(c)); Neg #eta", kTH3D, {axisPt, axisPt, axisEta});
    histos.add("XiccProngs/h3dBach", "h3dBach; Xicc pT (GeV/#it(c)); Bach pT (GeV/#it(c)); Bach #eta", kTH3D, {axisPt, axisPt, axisEta});
    histos.add("XiccProngs/h3dPi1c", "h3dPi1c; Xicc pT (GeV/#it(c)); Pi1c pT (GeV/#it(c)); Pi1c #eta", kTH3D, {axisPt, axisPt, axisEta});
    histos.add("XiccProngs/h3dPi2c", "h3dPi2c; Xicc pT (GeV/#it(c)); Pi2c pT (GeV/#it(c)); Pi2c #eta", kTH3D, {axisPt, axisPt, axisEta});
    histos.add("XiccProngs/h3dPicc", "h3dPicc; Xicc pT (GeV/#it(c)); Picc pT (GeV/#it(c)); Picc #eta", kTH3D, {axisPt, axisPt, axisEta});
    histos.add("h3dXicc", "h3dXicc; Xicc pT (GeV/#it(c)); Xicc #eta; Xicc mass (GeV/#it(c)^{2})", kTH3D, {axisPt, axisEta, axisXiccMass});
  }

  void processXicc(aod::MCharmCores const& multiCharmTracks)
  {
    for (const auto& xiccCand : multiCharmTracks) {
      if (xiccCand.xicDauDCA() > xicMaxDauDCA || xiccCand.xiccDauDCA() > xiccMaxDauDCA)
        continue;

      if (std::fabs(xiccCand.xiDCAxy()) < xiMinDCAxy || std::fabs(xiccCand.xiDCAz()) < xiMinDCAz)
        continue;

      if (std::fabs(xiccCand.pi1cDCAxy()) < picMinDCAxy || std::fabs(xiccCand.pi1cDCAz()) < picMinDCAz)
        continue;

      if (std::fabs(xiccCand.pi2cDCAxy()) < picMinDCAxy || std::fabs(xiccCand.pi2cDCAz()) < picMinDCAz)
        continue;

      if (std::fabs(xiccCand.piccDCAxy()) < piccMinDCAxy || std::fabs(xiccCand.piccDCAz()) < piccMinDCAz)
        continue;

      if (std::fabs(xiccCand.xicDCAxy()) < xicMinDCAxy || std::fabs(xiccCand.xicDCAz()) < xicMinDCAz)
        continue;

      if (std::fabs(xiccCand.pi1cDCAxy()) < picMinDCAxy || std::fabs(xiccCand.pi1cDCAz()) < picMinDCAz)
        continue;

      if (std::fabs(xiccCand.pi2cDCAxy()) < picMinDCAxy || std::fabs(xiccCand.pi2cDCAz()) < picMinDCAz)
        continue;

      if (std::fabs(xiccCand.xiccDCAxy()) > xiccMaxDCAxy || std::fabs(xiccCand.xiccDCAz()) > xiccMaxDCAz)
        continue;

      // Cut on time delta as LoI for now
      if (xiccCand.pi1cTofDeltaInner() > picMaxTofDiffInner)
        continue;

      if (xiccCand.pi2cTofDeltaInner() > picMaxTofDiffInner)
        continue;

      if (xiccCand.piccTofDeltaInner() > piccMaxTofDiffInner)
        continue;

      if (xiccCand.pi1cPt() < picMinPt || xiccCand.pi2cPt() < picMinPt)
        continue;

      if (xiccCand.piccPt() < piccMinPt)
        continue;

      if (xiccCand.xicDecayRadius2D() < xicMinRadius)
        continue;

      if (xiccCand.xiccDecayRadius2D() < xiccMinRadius)
        continue;

      if (xiccCand.xicProperLength() < xicMinProperLength || xiccCand.xicProperLength() > xicMaxProperLength)
        continue;

      if (xiccCand.xiccProperLength() < xiccMinProperLength || xiccCand.xiccProperLength() > xiccMaxProperLength)
        continue;

      if (xiccCand.xicDistanceFromPV() < xicMinDecayDistanceFromPV)
        continue;

      histos.fill(HIST("SelectionQA/hDCAXicDaughters"), xiccCand.xicDauDCA() * 1e+4);
      histos.fill(HIST("SelectionQA/hDCAXiccDaughters"), xiccCand.xiccDauDCA() * 1e+4);
      histos.fill(HIST("SelectionQA/hDCAxyXi"), std::fabs(xiccCand.xiDCAxy() * 1e+4));
      histos.fill(HIST("SelectionQA/hDCAzXi"), std::fabs(xiccCand.xiDCAz() * 1e+4));
      histos.fill(HIST("SelectionQA/hDCAxyXic"), std::fabs(xiccCand.xicDCAxy() * 1e+4));
      histos.fill(HIST("SelectionQA/hDCAzXic"), std::fabs(xiccCand.xicDCAz() * 1e+4));
      histos.fill(HIST("SelectionQA/hDCAxyXicc"), std::fabs(xiccCand.xiccDCAxy() * 1e+4));
      histos.fill(HIST("SelectionQA/hDCAzXicc"), std::fabs(xiccCand.xiccDCAz() * 1e+4));
      histos.fill(HIST("SelectionQA/hPi1cPt"), xiccCand.pi1cPt());
      histos.fill(HIST("SelectionQA/hPi2cPt"), xiccCand.pi2cPt());
      histos.fill(HIST("SelectionQA/hPiccPt"), xiccCand.piccPt());
      histos.fill(HIST("SelectionQA/hXicDecayRadius"), xiccCand.xicDecayRadius2D() * 1e+4);
      histos.fill(HIST("SelectionQA/hXiccDecayRadius"), xiccCand.xiccDecayRadius2D() * 1e+4);
      histos.fill(HIST("SelectionQA/hXicDecayDistanceFromPV"), xiccCand.xicDistanceFromPV() * 1e+4);
      histos.fill(HIST("SelectionQA/hProperLengthXic"), xiccCand.xicProperLength() * 1e+4);
      histos.fill(HIST("SelectionQA/hProperLengthXicc"), xiccCand.xiccProperLength() * 1e+4);
      histos.fill(HIST("SelectionQA/hInnerTofTimeDeltaPi1c"), xiccCand.pi1cTofDeltaInner());
      histos.fill(HIST("SelectionQA/hInnerTofTimeDeltaPi2c"), xiccCand.pi2cTofDeltaInner());
      histos.fill(HIST("SelectionQA/hInnerTofTimeDeltaPicc"), xiccCand.piccTofDeltaInner());

      histos.fill(HIST("XiccProngs/h3dNeg"), xiccCand.pt(), xiccCand.negPt(), xiccCand.negEta());
      histos.fill(HIST("XiccProngs/h3dPos"), xiccCand.pt(), xiccCand.posPt(), xiccCand.posEta());
      histos.fill(HIST("XiccProngs/h3dBach"), xiccCand.pt(), xiccCand.bachPt(), xiccCand.bachEta());
      histos.fill(HIST("XiccProngs/h3dPi1c"), xiccCand.pt(), xiccCand.pi1cPt(), xiccCand.pi1cEta());
      histos.fill(HIST("XiccProngs/h3dPi2c"), xiccCand.pt(), xiccCand.pi2cPt(), xiccCand.pi2cEta());
      histos.fill(HIST("XiccProngs/h3dPicc"), xiccCand.pt(), xiccCand.piccPt(), xiccCand.piccEta());
      histos.fill(HIST("h3dXicc"), xiccCand.pt(), xiccCand.eta(), xiccCand.xiccMass());
    }
  }

  PROCESS_SWITCH(alice3multicharm, processXicc, "find Xicc baryons", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<alice3multicharm>(cfgc)};
}
