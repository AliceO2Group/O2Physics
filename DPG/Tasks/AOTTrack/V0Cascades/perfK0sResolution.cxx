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
/// \file perfK0sResolution.cxx
/// \brief V0s (K0s, Lambda and antiLambda) analysis task
///
/// \author Nicolò Jacazio <nicolo.jacazio@cern.ch>, Universita del Piemonte Orientale
/// \author Roman Nepeivoda <roman.nepeivoda@cern.ch>, Lund University
/// \author Romain Schotter <romain.schotter@cern.ch>, Austrian Academy of Sciences & MBI
//

#include "PWGLF/DataModel/LFStrangenessPIDTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Tools/TrackTuner.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/PhysicsConstants.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <DetectorsBase/MatLayerCylSet.h>
#include <DetectorsBase/Propagator.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>
#include <ReconstructionDataFormats/DCA.h>
#include <ReconstructionDataFormats/PID.h>
#include <ReconstructionDataFormats/TrackParametrizationWithError.h>
#include <ReconstructionDataFormats/Vertex.h>

#include <TH1.h>
#include <TString.h>

#include <RtypesCore.h>

#include <array>
#include <cstdint>
#include <string>

using namespace o2;
using namespace o2::track;
using namespace o2::framework;
using namespace o2::framework::expressions;

using PIDTracksIU = soa::Join<aod::TracksIU, aod::TracksExtra, aod::pidTPCFullPi, aod::pidTPCFullPr, aod::pidTOFFullPi>;
using PIDTracksIUMC = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::TracksDCA, aod::TracksDCACov, aod::McTrackLabels, aod::pidTPCFullPi, aod::pidTPCFullPr, aod::pidTOFFullPi>;
using SelectedCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::PVMults, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::CentFV0As, aod::CentFT0CVariant1s, aod::CentNGlobals, aod::CentMFTs>;

enum CentEstimator {
  kCentFT0C = 0,
  kCentFT0M,
  kCentFT0CVariant1,
  kCentMFT,
  kCentNGlobal
};

struct perfK0sResolution {
  // Configurable bins
  ConfigurableAxis mBins{"mBins", {200, 0.4f, 0.6f}, "Mass binning"};
  ConfigurableAxis pTBins{"pTBins", {240, 0.f, 12.f}, "pT binning"};
  ConfigurableAxis invpTBins{"invpTBins", {240, 0.f, 12.f}, "inverse pT binning"};
  ConfigurableAxis pTResBins{"pTResBins", {200, -1.2f, 1.2f}, "pT resolution binning"};
  ConfigurableAxis pTResRelBins{"pTResRelBins", {200, -0.2f, 0.2f}, "pT relative resolution binning"};
  ConfigurableAxis invpTResBins{"invpTResBins", {200, -1.2f, 1.2f}, "inv pT resolution binning"};
  ConfigurableAxis invpTResNormBins{"invpTResNormBins", {200, -4.f, 4.f}, "inv pT normalised resolution binning"};
  ConfigurableAxis etaBins{"etaBins", {2, -1.f, 1.f}, "eta binning"};
  ConfigurableAxis etaBinsDauthers{"etaBinsDauthers", {100, -1.f, 1.f}, "eta binning for daughters"};
  ConfigurableAxis phiBins{"phiBins", {100, 0.f, 6.28f}, "phi binning"};
  ConfigurableAxis relpTResBins{"relpTResBins", {200, 0.f, 0.5f}, "rel. pT resolution binning"};
  ConfigurableAxis centralityAxis{"centralityAxis", {VARIABLE_WIDTH, 0.0f, 5.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f, 100.0f, 101.f}, "Centrality"};
  ConfigurableAxis occupancyAxis{"occupancyAxis", {VARIABLE_WIDTH, -1.0f, 0.0f, 250.0f, 500.0f, 750.0f, 1000.0f, 1500.0f, 2000.0f, 3000.0f, 4500.0f, 6000.0f, 8000.0f, 10000.0f, 50000.0f}, "Occupancy"};

  struct : ConfigurableGroup {
    std::string prefix = "eventSelections"; // JSON group name
    Configurable<bool> requireSel8{"requireSel8", true, "require sel8 event selection"};
    Configurable<bool> requireTriggerTVX{"requireTriggerTVX", true, "require FT0 vertex (acceptable FT0C-FT0A time difference) at trigger level"};
    Configurable<bool> rejectITSROFBorder{"rejectITSROFBorder", true, "reject events at ITS ROF border (Run 3 only)"};
    Configurable<bool> rejectTFBorder{"rejectTFBorder", true, "reject events at TF border (Run 3 only)"};
    Configurable<bool> requireIsVertexITSTPC{"requireIsVertexITSTPC", false, "require events with at least one ITS-TPC track (Run 3 only)"};
    Configurable<bool> requireIsGoodZvtxFT0VsPV{"requireIsGoodZvtxFT0VsPV", true, "require events with PV position along z consistent (within 1 cm) between PV reconstructed using tracks and PV using FT0 A-C time difference (Run 3 only)"};
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

    Configurable<int> centralityEstimator{"centralityEstimator", kCentFT0C, "Run 3 centrality estimator (0:CentFT0C, 1:CentFT0M, 2:CentFT0CVariant1, 3:CentMFT, 4:CentNGlobal)"};
    Configurable<bool> useFT0CbasedOccupancy{"useFT0CbasedOccupancy", false, "Use sum of FT0-C amplitudes for estimating occupancy? (if not, use track-based definition)"};
    // fast check on occupancy
    Configurable<float> minOccupancy{"minOccupancy", -1, "minimum occupancy from neighbouring collisions"};
    Configurable<float> maxOccupancy{"maxOccupancy", -1, "maximum occupancy from neighbouring collisions"};
  } eventSelections;

  static constexpr float DefaultLifetimeCuts[1][2] = {{30., 20.}};

  struct : ConfigurableGroup {
    std::string prefix = "v0Selections"; // JSON group name
    Configurable<int> v0TypeSelection{"v0TypeSelection", 1, "select on a certain V0 type (leave negative if no selection desired)"};

    // Selection criteria: acceptance
    Configurable<float> rapidityCut{"rapidityCut", 0.5, "rapidity"};
    Configurable<float> daughterEtaCut{"daughterEtaCut", 0.8, "max eta for daughters"};

    // Standard 5 topological criteria
    Configurable<float> v0cospa{"v0cospa", 0.995, "min V0 CosPA"};
    Configurable<float> dcav0dau{"dcav0dau", 1.0, "max DCA V0 Daughters (cm)"};
    Configurable<float> dcanegtopv{"dcanegtopv", 0.1, "min DCA Neg To PV (cm)"};
    Configurable<float> dcapostopv{"dcapostopv", 0.1, "min DCA Pos To PV (cm)"};
    Configurable<float> v0radius{"v0radius", 1.2, "minimum V0 radius (cm)"};
    Configurable<float> v0radiusMax{"v0radiusMax", 1E5, "maximum V0 radius (cm)"};
    Configurable<LabeledArray<float>> lifetimecut{"lifetimecut", {DefaultLifetimeCuts[0], 2, {"lifetimecutLambda", "lifetimecutK0S"}}, "lifetimecut"};

    // Additional selection on the AP plot (exclusive for K0Short)
    // original equation: lArmPt*5>TMath::Abs(lArmAlpha)
    Configurable<float> armPodCut{"armPodCut", 5.0f, "pT * (cut) > |alpha|, AP cut. Negative: no cut"};

    // Track quality
    Configurable<int> minTPCrows{"minTPCrows", 70, "minimum TPC crossed rows"};
    Configurable<int> minITSclusters{"minITSclusters", -1, "minimum ITS clusters"};
    Configurable<float> minTPCrowsOverFindableClusters{"minTPCrowsOverFindableClusters", -1, "minimum nbr of TPC crossed rows over findable clusters"};
    Configurable<float> minTPCfoundOverFindableClusters{"minTPCfoundOverFindableClusters", -1, "minimum nbr of found over findable TPC clusters"};
    Configurable<float> maxFractionTPCSharedClusters{"maxFractionTPCSharedClusters", 1e+09, "maximum fraction of TPC shared clusters"};
    Configurable<float> maxITSchi2PerNcls{"maxITSchi2PerNcls", 1e+09, "maximum ITS chi2 per clusters"};
    Configurable<float> maxTPCchi2PerNcls{"maxTPCchi2PerNcls", 1e+09, "maximum TPC chi2 per clusters"};
    Configurable<int> requirePosITSib{"requirePosITSib", 0, "require ITS IB selection on positive daughters? -1: no ITS IB, 0: no selection, 1: ITS IB"};
    Configurable<int> requireNegITSib{"requireNegITSib", 0, "require ITS IB selection on negative daughters? -1: no ITS IB, 0: no selection, 1: ITS IB"};
    Configurable<int> requirePosITSafterburner{"requirePosITSafterburner", 0, "require positive track formed out of afterburner ITS tracks? -1: no AB, 0: no selection, 1: AB"};
    Configurable<int> requireNegITSafterburner{"requireNegITSafterburner", 0, "require negative track formed out of afterburner ITS tracks? -1: no AB, 0: no selection, 1: AB"};
    Configurable<int> requirePosTRD{"requirePosTRD", 0, "require TRD selection on positive daughters? -1: no TRD, 0: no selection, 1: TRD"};
    Configurable<int> requireNegTRD{"requireNegTRD", 0, "require TRD selection on negative daughters? -1: no TRD, 0: no selection, 1: TRD"};
    Configurable<int> requirePosTOF{"requirePosTOF", 0, "require TOF selection on positive daughters? -1: no TOF, 0: no selection, 1: TOF"};
    Configurable<int> requireNegTOF{"requireNegTOF", 0, "require TOF selection on negative daughters? -1: no TOF, 0: no selection, 1: TOF"};
    Configurable<int> requirePosPIDforTracking{"requirePosPIDforTracking", -1, "require specific PID hypothesis used in tracking for the positive daughters? -1: no selection, 0: Electron, 1: Muon, 2: Pion, 3: Kaon, 4: Proton"};
    Configurable<int> requireNegPIDforTracking{"requireNegPIDforTracking", -1, "require specific PID hypothesis used in tracking for the negative daughters? -1: no selection, 0: Electron, 1: Muon, 2: Pion, 3: Kaon, 4: Proton"};

    // PID (TPC/TOF)
    Configurable<float> tpcPidNsigmaCut{"tpcPidNsigmaCut", 10., "tpcPidNsigmaCut"};
    Configurable<float> tofPidNsigmaCut{"tofPidNsigmaCut", 1e+6, "tofPidNsigmaCut"};
  } v0Selections;

  // Configure plots to enable
  Configurable<bool> useMultidimHisto{"useMultidimHisto", false, "use multidimentional histograms"};
  Configurable<bool> enableTPCPlot{"enableTPCPlot", false, "Enable the TPC plot"};
  Configurable<bool> requireTrueK0s{"requireTrueK0s", true, "require rec. v0 to be true K0s"};
  Configurable<bool> doTreatPiToMuon{"doTreatPiToMuon", false, "Take pi decay into muon into account in MC"};

  HistogramRegistry rK0sResolution{"K0sResolution", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rK0sDauResolution{"K0sDauResolution", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  Service<o2::ccdb::BasicCCDBManager> ccdb;

  // for TrackTuner only (MC smearing)
  o2::parameters::GRPMagField* grpmag = nullptr;
  o2::base::MatLayerCylSet* lut = nullptr;
  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;
  TrackTuner trackTunerObj;

  Configurable<bool> useTrackTuner{"useTrackTuner", false, "Apply Improver/DCA corrections to MC"};
  Configurable<std::string> trackTunerParams{"trackTunerParams", "debugInfo=0|updateTrackCovMat=0|updateCurvature=1|updatePulls=0|isInputFileFromCCDB=1|pathInputFile=Users/m/mfaggin/test/inputsTrackTuner/PbPb2022|nameInputFile=trackTuner_DataLHC22sPass5_McLHC22l1b2_run529397.root|usePvRefitCorrections=0|oneOverPtCurrent=1|oneOverPtUpgr=1.2", "TrackTuner parameter initialization (format: <name>=<value>|<name>=<value>)"};
  OutputObj<TH1D> trackTunedTracks{TH1D("trackTunedTracks", "", 4, 0.5, 4.5), OutputObjHandlingPolicy::AnalysisObject};
  Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUTInner", "Path of the Lut parametrization"};
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};

  int runNumber = -1;

  void init(InitContext const&)
  {
    const AxisSpec statAxis{5, 0, 5, ""};
    const AxisSpec mAxis{mBins, "#it{m} (GeV/#it{c}^{2})"};
    const AxisSpec pTAxis{pTBins, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec invpTAxis{invpTBins, "1/#it{p}_{T}^{MC} (GeV/#it{c})^{-1}"};
    const AxisSpec pTResAxis{pTResBins, "#Delta#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec pTResRelAxis{pTResRelBins, "(#it{p}_{T}^{rec} - #it{p}_{T}^{MC})/#it{p}_{T}^{MC}"};
    const AxisSpec invpTResAxis{invpTResBins, "1/#it{p}_{T}-1/#it{p}_{T}^{MC} (GeV/#it{c})^{-1}"};
    const AxisSpec invpTResNormAxis{invpTResNormBins, "(1/#it{p}_{T}-1/#it{p}_{T}^{MC})/#sigma_{1/#it{p}_{T}}"};
    const AxisSpec relpTResAxis{relpTResBins, "#sigma(#it{p}_{T})/#it{p}_{T}"};
    const AxisSpec etaAxis{etaBins, "#eta"};
    const AxisSpec etaAxisPosD{etaBinsDauthers, "#eta pos."};
    const AxisSpec etaAxisNegD{etaBinsDauthers, "#eta neg."};
    const AxisSpec phiAxis{phiBins, "#phi"};
    const AxisSpec trueK0Axis{2, -0.5, 1.5, "True K0"};

    // Event Counters
    rK0sResolution.add("hEventSelection", "hEventSelection", kTH1D, {{21, -0.5f, +20.5f}});
    rK0sResolution.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(1, "All collisions");
    rK0sResolution.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(2, "sel8 cut");
    rK0sResolution.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(3, "kIsTriggerTVX");
    rK0sResolution.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(4, "kNoITSROFrameBorder");
    rK0sResolution.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(5, "kNoTimeFrameBorder");
    rK0sResolution.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(6, "posZ cut");
    rK0sResolution.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(7, "kIsVertexITSTPC");
    rK0sResolution.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(8, "kIsGoodZvtxFT0vsPV");
    rK0sResolution.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(9, "kIsVertexTOFmatched");
    rK0sResolution.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(10, "kIsVertexTRDmatched");
    rK0sResolution.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(11, "kNoSameBunchPileup");
    rK0sResolution.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(12, "kNoCollInTimeRangeStd");
    rK0sResolution.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(13, "kNoCollInTimeRangeStrict");
    rK0sResolution.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(14, "kNoCollInTimeRangeNarrow");
    rK0sResolution.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(15, "kNoCollInRofStd");
    rK0sResolution.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(16, "kNoCollInRofStrict");
    rK0sResolution.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(17, "INEL>0");
    rK0sResolution.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(18, "INEL>1");
    rK0sResolution.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(19, "Below min occup.");
    rK0sResolution.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(20, "Above max occup.");

    rK0sResolution.add("hEventCentrality", "hEventCentrality", kTH1D, {{101, 0.0f, 101.0f}});
    rK0sResolution.add("hEventOccupancy", "hEventOccupancy", kTH1D, {occupancyAxis});

    rK0sResolution.add("h1_stats", "h1_stats", {HistType::kTH1D, {statAxis}});
    TString hStatsLabels[5] = {"Selected Events", "All V0s", "Selected V0s", "Daughters have MC particles", "Daughters corr. rec."};
    for (int n = 1; n <= rK0sResolution.get<TH1>(HIST("h1_stats"))->GetNbinsX(); n++) {
      rK0sResolution.get<TH1>(HIST("h1_stats"))->GetXaxis()->SetBinLabel(n, hStatsLabels[n - 1]);
    }

    if (doprocessMC) {
      rK0sDauResolution.add("h2_massPosPtRes", "h2_massPosPtRes", {HistType::kTH2D, {mAxis, pTResAxis}});
      rK0sDauResolution.add("h2_massNegPtRes", "h2_massNegPtRes", {HistType::kTH2D, {mAxis, pTResAxis}});

      rK0sDauResolution.add("h2_genPtPosPtResNorm", "h2_genPtPosPtResNorm", {HistType::kTH2D, {pTResRelAxis, pTAxis}});
      rK0sDauResolution.add("h2_genPxPosPxResNorm", "h2_genPxPosPxResNorm", {HistType::kTH2D, {pTResRelAxis, pTAxis}});
      rK0sDauResolution.add("h2_genPyPosPyResNorm", "h2_genPyPosPyResNorm", {HistType::kTH2D, {pTResRelAxis, pTAxis}});
      rK0sDauResolution.add("h2_genPzPosPzResNorm", "h2_genPzPosPzResNorm", {HistType::kTH2D, {pTResRelAxis, pTAxis}});

      rK0sDauResolution.add("h2_genPtNegPtResNorm", "h2_genPtNegPtResNorm", {HistType::kTH2D, {pTResRelAxis, pTAxis}});
      rK0sDauResolution.add("h2_genPxNegPxResNorm", "h2_genPxNegPxResNorm", {HistType::kTH2D, {pTResRelAxis, pTAxis}});
      rK0sDauResolution.add("h2_genPyNegPyResNorm", "h2_genPyNegPyResNorm", {HistType::kTH2D, {pTResRelAxis, pTAxis}});
      rK0sDauResolution.add("h2_genPzNegPzResNorm", "h2_genPzNegPzResNorm", {HistType::kTH2D, {pTResRelAxis, pTAxis}});

      rK0sDauResolution.add("h2_genPtPosPtRes", "h2_genPtPosPtRes", {HistType::kTH2D, {pTResAxis, pTAxis}});
      rK0sDauResolution.add("h2_genPxPosPxRes", "h2_genPxPosPxRes", {HistType::kTH2D, {pTResAxis, pTAxis}});
      rK0sDauResolution.add("h2_genPyPosPyRes", "h2_genPyPosPyRes", {HistType::kTH2D, {pTResAxis, pTAxis}});
      rK0sDauResolution.add("h2_genPzPosPzRes", "h2_genPzPosPzRes", {HistType::kTH2D, {pTResAxis, pTAxis}});

      rK0sDauResolution.add("h2_genPtNegPtRes", "h2_genPtNegPtRes", {HistType::kTH2D, {pTResAxis, pTAxis}});
      rK0sDauResolution.add("h2_genPxNegPxRes", "h2_genPxNegPxRes", {HistType::kTH2D, {pTResAxis, pTAxis}});
      rK0sDauResolution.add("h2_genPyNegPyRes", "h2_genPyNegPyRes", {HistType::kTH2D, {pTResAxis, pTAxis}});
      rK0sDauResolution.add("h2_genPzNegPzRes", "h2_genPzNegPzRes", {HistType::kTH2D, {pTResAxis, pTAxis}});

      rK0sDauResolution.add("h2_genPtPosPulls", "h2_GenPtPosPulls", {HistType::kTH2D, {invpTResNormAxis, invpTAxis}});
      rK0sDauResolution.add("h2_genPtNegPulls", "h2_GenPtNegPulls", {HistType::kTH2D, {invpTResNormAxis, invpTAxis}});

      rK0sDauResolution.add("h2_PosRelPtRes", "h2_PosRelPtRes", {HistType::kTH2D, {pTAxis, relpTResAxis}});
      rK0sDauResolution.add("h2_NegRelPtRes", "h2_NegRelPtRes", {HistType::kTH2D, {pTAxis, relpTResAxis}});
    }
    rK0sResolution.add("h2_masspT", "h2_masspT", {HistType::kTH2D, {mAxis, pTAxis}});
    rK0sResolution.add("h2_masseta", "h2_masseta", {HistType::kTH2D, {mAxis, etaAxis}});
    rK0sResolution.add("h2_massphi", "h2_massphi", {HistType::kTH2D, {mAxis, phiAxis}});
    if (useMultidimHisto) {
      if (doprocessMC) {
        rK0sResolution.add("thn_mass", "thn_mass", kTHnSparseF, {mAxis, pTAxis, etaAxis, phiAxis, etaAxisPosD, etaAxisNegD, invpTResAxis, invpTResAxis, trueK0Axis});
      } else {
        rK0sResolution.add("thn_mass", "thn_mass", kTHnSparseF, {mAxis, pTAxis, etaAxis, phiAxis, etaAxisPosD, etaAxisNegD});
      }
    }
    rK0sResolution.add("h3_centralitypTmass", "h3_centralitypTmass", kTH3D, {centralityAxis, pTAxis, mAxis});
    rK0sResolution.add("h3_occupancypTmass", "h3_occupancypTmass", kTH3D, {occupancyAxis, pTAxis, mAxis});

    if (enableTPCPlot) {
      rK0sDauResolution.add("h3_tpc_vs_pid_hypothesis", "h3_tpc_vs_pid_hypothesis", {HistType::kTH3D, {{200, -10.f, 10.f, "#it{p}/Z (GeV/#it{c})"}, {1000, 0, 1000.f, "dE/dx (a.u.)"}, {10, -0.5, 9.5f, "PID hypothesis"}}});
    }

    /// TrackTuner initialization
    if (useTrackTuner) {
      ccdb->setURL(ccdbUrl);
      ccdb->setCaching(true);
      ccdb->setLocalObjectValidityChecking();
      ccdb->setFatalWhenNull(false);
      lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>(lutPath));
      std::string outputStringParams = trackTunerObj.configParams(trackTunerParams);
      trackTunerObj.getDcaGraphs();
      // QA is done in tuneTrackParams method
      trackTunedTracks->SetTitle(outputStringParams.c_str());
      trackTunedTracks->GetXaxis()->SetBinLabel(1, "all tracks");
      trackTunedTracks->GetXaxis()->SetBinLabel(2, "tracks tuned (no negative detXY)");
      trackTunedTracks->GetXaxis()->SetBinLabel(3, "untouched tracks due to negative detXY");
      trackTunedTracks->GetXaxis()->SetBinLabel(4, "original detXY<0");
    }

    // inspect histogram sizes, please
    rK0sResolution.print();
    rK0sDauResolution.print();
  }

  void initCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {
    if (runNumber == bc.runNumber()) {
      return;
    }
    grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, bc.timestamp());
    if (!grpmag) {
      LOG(fatal) << "Got nullptr from CCDB for path " << grpmagPath;
    }
    LOG(info) << "Setting magnetic field to current " << grpmag->getL3Current() << " A for run " << bc.runNumber() << " from its GRPMagField CCDB object";
    o2::base::Propagator::initFieldFromGRP(grpmag);
    o2::base::Propagator::Instance()->setMatLUT(lut);
    runNumber = bc.runNumber();
  }

  // ______________________________________________________
  // Return slicing output
  template <typename TCollision>
  float getCentralityRun3(TCollision const& collision)
  {
    if (eventSelections.centralityEstimator == kCentFT0C)
      return collision.centFT0C();
    else if (eventSelections.centralityEstimator == kCentFT0M)
      return collision.centFT0M();
    else if (eventSelections.centralityEstimator == kCentFT0CVariant1)
      return collision.centFT0CVariant1();
    else if (eventSelections.centralityEstimator == kCentMFT)
      return collision.centMFT();
    else if (eventSelections.centralityEstimator == kCentNGlobal)
      return collision.centNGlobal();

    return -1.f;
  }

  template <typename TCollision>
  bool isEventAccepted(TCollision collision, bool fillHists)
  // check whether the collision passes our collision selections
  {
    if (fillHists)
      rK0sResolution.fill(HIST("hEventSelection"), 0. /* all collisions */);

    if (eventSelections.requireSel8 && !collision.sel8()) {
      return false;
    }
    if (fillHists)
      rK0sResolution.fill(HIST("hEventSelection"), 1 /* sel8 collisions */);

    if (eventSelections.requireTriggerTVX && !collision.selection_bit(aod::evsel::kIsTriggerTVX)) {
      return false;
    }
    if (fillHists)
      rK0sResolution.fill(HIST("hEventSelection"), 2 /* FT0 vertex (acceptable FT0C-FT0A time difference) collisions */);

    if (eventSelections.rejectITSROFBorder && !collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
      return false;
    }
    if (fillHists)
      rK0sResolution.fill(HIST("hEventSelection"), 3 /* Not at ITS ROF border */);

    if (eventSelections.rejectTFBorder && !collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
      return false;
    }
    if (fillHists)
      rK0sResolution.fill(HIST("hEventSelection"), 4 /* Not at TF border */);

    if (std::abs(collision.posZ()) > eventSelections.maxZVtxPosition) {
      return false;
    }
    if (fillHists)
      rK0sResolution.fill(HIST("hEventSelection"), 5 /* vertex-Z selected */);

    if (eventSelections.requireIsVertexITSTPC && !collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) {
      return false;
    }
    if (fillHists)
      rK0sResolution.fill(HIST("hEventSelection"), 6 /* Contains at least one ITS-TPC track */);

    if (eventSelections.requireIsGoodZvtxFT0VsPV && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return false;
    }
    if (fillHists)
      rK0sResolution.fill(HIST("hEventSelection"), 7 /* PV position consistency check */);

    if (eventSelections.requireIsVertexTOFmatched && !collision.selection_bit(o2::aod::evsel::kIsVertexTOFmatched)) {
      return false;
    }
    if (fillHists)
      rK0sResolution.fill(HIST("hEventSelection"), 8 /* PV with at least one contributor matched with TOF */);

    if (eventSelections.requireIsVertexTRDmatched && !collision.selection_bit(o2::aod::evsel::kIsVertexTRDmatched)) {
      return false;
    }
    if (fillHists)
      rK0sResolution.fill(HIST("hEventSelection"), 9 /* PV with at least one contributor matched with TRD */);

    if (eventSelections.rejectSameBunchPileup && !collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      return false;
    }
    if (fillHists)
      rK0sResolution.fill(HIST("hEventSelection"), 10 /* Not at same bunch pile-up */);

    if (eventSelections.requireNoCollInTimeRangeStd && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
      return false;
    }
    if (fillHists)
      rK0sResolution.fill(HIST("hEventSelection"), 11 /* No other collision within +/- 2 microseconds or mult above a certain threshold in -4 - -2 microseconds*/);

    if (eventSelections.requireNoCollInTimeRangeStrict && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStrict)) {
      return false;
    }
    if (fillHists)
      rK0sResolution.fill(HIST("hEventSelection"), 12 /* No other collision within +/- 10 microseconds */);

    if (eventSelections.requireNoCollInTimeRangeNarrow && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeNarrow)) {
      return false;
    }
    if (fillHists)
      rK0sResolution.fill(HIST("hEventSelection"), 13 /* No other collision within +/- 2 microseconds */);

    if (eventSelections.requireNoCollInROFStd && !collision.selection_bit(o2::aod::evsel::kNoCollInRofStandard)) {
      return false;
    }
    if (fillHists)
      rK0sResolution.fill(HIST("hEventSelection"), 14 /* No other collision within the same ITS ROF with mult. above a certain threshold */);

    if (eventSelections.requireNoCollInROFStrict && !collision.selection_bit(o2::aod::evsel::kNoCollInRofStrict)) {
      return false;
    }
    if (fillHists)
      rK0sResolution.fill(HIST("hEventSelection"), 15 /* No other collision within the same ITS ROF */);

    if (eventSelections.requireINEL0 && collision.multNTracksPVeta1() < 1) {
      return false;
    }
    if (fillHists)
      rK0sResolution.fill(HIST("hEventSelection"), 16 /* INEL > 0 */);

    if (eventSelections.requireINEL1 && collision.multNTracksPVeta1() < 2) {
      return false;
    }
    if (fillHists)
      rK0sResolution.fill(HIST("hEventSelection"), 17 /* INEL > 1 */);

    float collisionOccupancy = eventSelections.useFT0CbasedOccupancy ? collision.ft0cOccupancyInTimeRange() : collision.trackOccupancyInTimeRange();
    if (eventSelections.minOccupancy >= 0 && collisionOccupancy < eventSelections.minOccupancy) {
      return false;
    }
    if (fillHists)
      rK0sResolution.fill(HIST("hEventSelection"), 18 /* Below min occupancy */);

    if (eventSelections.maxOccupancy >= 0 && collisionOccupancy > eventSelections.maxOccupancy) {
      return false;
    }
    if (fillHists)
      rK0sResolution.fill(HIST("hEventSelection"), 19 /* Above max occupancy */);

    return true;
  }

  template <typename T1, typename T2, typename C>
  bool acceptV0(const T1& v0, const T2& ntrack, const T2& ptrack, const C& collision)
  {
    // Apply selections on V0
    if (std::abs(v0.yK0Short()) > v0Selections.rapidityCut) {
      return false;
    }

    if (std::abs(v0.negativeeta()) > v0Selections.daughterEtaCut || std::abs(v0.positiveeta()) > v0Selections.daughterEtaCut)
      return false;

    if (v0Selections.v0TypeSelection > -1 && v0.v0Type() != v0Selections.v0TypeSelection)
      return false; // skip V0s that are not standard

    // Base topological variables
    if (v0.v0radius() < v0Selections.v0radius)
      return false;
    if (v0.v0radius() > v0Selections.v0radiusMax)
      return false;
    if (std::abs(v0.dcapostopv()) < v0Selections.dcapostopv)
      return false;
    if (std::abs(v0.dcanegtopv()) < v0Selections.dcanegtopv)
      return false;
    if (v0.v0cosPA() < v0Selections.v0cospa)
      return false;
    if (v0.dcaV0daughters() > v0Selections.dcav0dau)
      return false;

    if (v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassK0Short > v0Selections.lifetimecut->get("lifetimecutK0S")) {
      return false;
    }

    if (v0Selections.armPodCut > 1e-4 && v0.qtarm() * v0Selections.armPodCut < std::abs(v0.alpha())) {
      return false;
    }

    // Apply selections on V0 daughters

    // ITS quality flags
    // check minium ITS clusters
    if (ptrack.itsNCls() < v0Selections.minITSclusters)
      return false;
    if (ntrack.itsNCls() < v0Selections.minITSclusters)
      return false;
    // check maximum ITS chi2 per clusters
    if (ptrack.itsChi2NCl() > v0Selections.maxITSchi2PerNcls)
      return false;
    if (ntrack.itsChi2NCl() > v0Selections.maxITSchi2PerNcls)
      return false;

    // TPC quality flags
    // check minimum TPC crossed rows
    if (ptrack.tpcNClsCrossedRows() < v0Selections.minTPCrows)
      return false;
    if (ntrack.tpcNClsCrossedRows() < v0Selections.minTPCrows)
      return false;
    // check maximum TPC chi2 per clusters
    if (ptrack.tpcChi2NCl() > v0Selections.maxTPCchi2PerNcls)
      return false;
    if (ntrack.tpcChi2NCl() > v0Selections.maxTPCchi2PerNcls)
      return false;
    // check minimum fraction of TPC rows over findable
    if (ptrack.tpcCrossedRowsOverFindableCls() < v0Selections.minTPCrowsOverFindableClusters)
      return false;
    if (ntrack.tpcCrossedRowsOverFindableCls() < v0Selections.minTPCrowsOverFindableClusters)
      return false;
    // check minimum fraction of found over findable TPC clusters
    if (ptrack.tpcFoundOverFindableCls() < v0Selections.minTPCfoundOverFindableClusters)
      return false;
    if (ntrack.tpcFoundOverFindableCls() < v0Selections.minTPCfoundOverFindableClusters)
      return false;
    // check the maximum fraction of allowed shared TPC clusters
    if (ptrack.tpcChi2NCl() > v0Selections.maxFractionTPCSharedClusters)
      return false;
    if (ntrack.tpcChi2NCl() > v0Selections.maxFractionTPCSharedClusters)
      return false;

    // ITS Inner Barrel selection
    if (std::abs(v0Selections.requirePosITSib) > 0) {
      if (v0Selections.requirePosITSib < 0 && ptrack.itsNClsInnerBarrel() > 0) // require no ITS IB
        return false;
      if (v0Selections.requirePosITSib > 0 && ptrack.itsNClsInnerBarrel() < 1) // require ITS IB
        return false;
    }
    if (std::abs(v0Selections.requireNegITSib) > 0) {
      if (v0Selections.requireNegITSib < 0 && ntrack.itsNClsInnerBarrel() > 0) // require no ITS IB
        return false;
      if (v0Selections.requireNegITSib > 0 && ntrack.itsNClsInnerBarrel() < 1) // require ITS IB
        return false;
    }

    // ITS AfterBurner selection
    if (std::abs(v0Selections.requirePosITSafterburner) > 0) {
      if (v0Selections.requirePosITSafterburner < 0 && ptrack.isITSAfterburner()) // require no ITS AB
        return false;
      if (v0Selections.requirePosITSafterburner > 0 && !ptrack.isITSAfterburner()) // require ITS AB
        return false;
    }
    if (std::abs(v0Selections.requireNegITSafterburner) > 0) {
      if (v0Selections.requireNegITSafterburner < 0 && ntrack.isITSAfterburner()) // require no ITS AB
        return false;
      if (v0Selections.requireNegITSafterburner > 0 && !ntrack.isITSAfterburner()) // require ITS AB
        return false;
    }

    // TPC PID selection
    if (std::abs(ptrack.tpcNSigmaPi()) > v0Selections.tpcPidNsigmaCut) {
      return false;
    }
    if (std::abs(ntrack.tpcNSigmaPi()) > v0Selections.tpcPidNsigmaCut) {
      return false;
    }

    // TOF PID selection
    if (ptrack.hasTOF() && std::fabs(v0.tofNSigmaK0PiPlus()) > v0Selections.tofPidNsigmaCut) {
      return false;
    }
    if (ntrack.hasTOF() && std::fabs(v0.tofNSigmaK0PiMinus()) > v0Selections.tofPidNsigmaCut) {
      return false;
    }

    // TOF selection
    if (std::abs(v0Selections.requirePosTOF) > 0) {
      if (v0Selections.requirePosTOF < 0 && ptrack.hasTOF()) // require no TOF
        return false;
      if (v0Selections.requirePosTOF > 0 && !ptrack.hasTOF()) // require TOF
        return false;
    }
    if (std::abs(v0Selections.requireNegTOF) > 0) {
      if (v0Selections.requireNegTOF < 0 && ntrack.hasTOF()) // require no TOF
        return false;
      if (v0Selections.requireNegTOF > 0 && !ntrack.hasTOF()) // require TOF
        return false;
    }

    // TRD selection
    if (std::abs(v0Selections.requirePosTRD) > 0) {
      if (v0Selections.requirePosTRD < 0 && ptrack.hasTRD()) // require no TRD
        return false;
      if (v0Selections.requirePosTRD > 0 && !ptrack.hasTRD()) // require TRD
        return false;
    }
    if (std::abs(v0Selections.requireNegTRD) > 0) {
      if (v0Selections.requireNegTRD < 0 && ntrack.hasTRD()) // require no TRD
        return false;
      if (v0Selections.requireNegTRD > 0 && !ntrack.hasTRD()) // require TRD
        return false;
    }

    // Specific PID for tracking selection
    if (v0Selections.requirePosPIDforTracking > -1 && ptrack.pidForTracking() != static_cast<uint32_t>(v0Selections.requirePosPIDforTracking)) {
      return false;
    }
    if (v0Selections.requireNegPIDforTracking > -1 && ntrack.pidForTracking() != static_cast<uint32_t>(v0Selections.requireNegPIDforTracking)) {
      return false;
    }

    return true;
  }

  void processData(SelectedCollisions::iterator const& collision,
                   soa::Join<aod::V0Datas, aod::V0TOFNSigmas> const& fullV0s,
                   PIDTracksIU const&)
  {
    if (!isEventAccepted(collision, true))
      return;

    float centrality = getCentralityRun3(collision);
    float occupancy = eventSelections.useFT0CbasedOccupancy ? collision.ft0cOccupancyInTimeRange() : collision.trackOccupancyInTimeRange();
    rK0sResolution.fill(HIST("hEventCentrality"), centrality);
    rK0sResolution.fill(HIST("hEventOccupancy"), occupancy);

    rK0sResolution.fill(HIST("h1_stats"), 0.5);
    for (const auto& v0 : fullV0s) {
      rK0sResolution.fill(HIST("h1_stats"), 1.5);
      const auto& posTrack = v0.posTrack_as<PIDTracksIU>();
      const auto& negTrack = v0.negTrack_as<PIDTracksIU>();
      if (!acceptV0(v0, negTrack, posTrack, collision))
        continue;
      rK0sResolution.fill(HIST("h1_stats"), 2.5);

      float mass = v0.mK0Short();

      rK0sResolution.fill(HIST("h2_masspT"), mass, v0.pt());
      rK0sResolution.fill(HIST("h2_masseta"), mass, v0.eta());
      rK0sResolution.fill(HIST("h2_massphi"), mass, v0.phi());
      if (useMultidimHisto) {
        rK0sResolution.fill(HIST("thn_mass"), mass, v0.pt(), v0.eta(), v0.phi(), posTrack.eta(), negTrack.eta());
      }
      rK0sResolution.fill(HIST("h3_centralitypTmass"), centrality, v0.pt(), mass);
      rK0sResolution.fill(HIST("h3_occupancypTmass"), occupancy, v0.pt(), mass);
      if (enableTPCPlot) {
        rK0sDauResolution.fill(HIST("h3_tpc_vs_pid_hypothesis"), posTrack.tpcInnerParam(), posTrack.tpcSignal(), posTrack.pidForTracking());
        rK0sDauResolution.fill(HIST("h3_tpc_vs_pid_hypothesis"), -negTrack.tpcInnerParam(), negTrack.tpcSignal(), negTrack.pidForTracking());
      }
    }
  }
  PROCESS_SWITCH(perfK0sResolution, processData, "Process data", true);

  // Running variables
  o2::dataformats::VertexBase mVtx;
  o2::dataformats::DCA mDcaInfoCovPos;
  o2::dataformats::DCA mDcaInfoCovNeg;
  o2::track::TrackParametrizationWithError<float> mTrackParCovPos;
  o2::track::TrackParametrizationWithError<float> mTrackParCovNeg;

  template <typename TV0, typename TV0Track>
  void tuneV0(TV0 const& v0,
              TV0Track const& posTrack,
              TV0Track const& negTrack,
              aod::McParticles const&,
              aod::BCsWithTimestamps const& bcs)
  {
    initCCDB(bcs.begin());
    trackTunedTracks->Fill(1, 2); // tune 2 tracks
    setTrackParCov(posTrack, mTrackParCovPos);
    setTrackParCov(negTrack, mTrackParCovNeg);
    mTrackParCovPos.setPID(posTrack.pidForTracking());
    mTrackParCovNeg.setPID(negTrack.pidForTracking());
    mDcaInfoCovPos.set(999, 999, 999, 999, 999);
    mDcaInfoCovNeg.set(999, 999, 999, 999, 999);
    auto mcParticlePos = posTrack.mcParticle();
    auto mcParticleNeg = negTrack.mcParticle();

    // LOG(info) << "Inside tuneTrack: before calling tuneTrackParams trackParCov.getY(): " << mTrackParCovPos.getY();
    trackTunerObj.tuneTrackParams(mcParticlePos, mTrackParCovPos, matCorr, &mDcaInfoCovPos, trackTunedTracks);
    trackTunerObj.tuneTrackParams(mcParticleNeg, mTrackParCovNeg, matCorr, &mDcaInfoCovNeg, trackTunedTracks);
    // LOG(info) << "Inside tuneTrack: after calling tuneTrackParams trackParCov.getY(): " << mTrackParCovPos.getY();
    //  trackTunedTracks->Fill(1, 2);
    mVtx.setPos({v0.x(), v0.y(), v0.z()});
    mVtx.setCov(v0.positionCovMat()[0], v0.positionCovMat()[1], v0.positionCovMat()[2], v0.positionCovMat()[3], v0.positionCovMat()[4], v0.positionCovMat()[5]);
    o2::base::Propagator::Instance()->propagateToDCABxByBz(mVtx, mTrackParCovPos, 2.f, matCorr, &mDcaInfoCovPos);
    o2::base::Propagator::Instance()->propagateToDCABxByBz(mVtx, mTrackParCovNeg, 2.f, matCorr, &mDcaInfoCovNeg);
  }

  void processMC(SelectedCollisions::iterator const& collision,
                 soa::Join<aod::V0Datas, aod::V0Covs, aod::V0DauCovs, aod::V0TOFNSigmas, aod::McV0Labels> const& fullV0s,
                 PIDTracksIUMC const&,
                 aod::McParticles const& mcParticles,
                 aod::BCsWithTimestamps const& bcs)
  {
    if (!isEventAccepted(collision, true))
      return;

    float centrality = getCentralityRun3(collision);
    float occupancy = eventSelections.useFT0CbasedOccupancy ? collision.ft0cOccupancyInTimeRange() : collision.trackOccupancyInTimeRange();
    rK0sResolution.fill(HIST("hEventCentrality"), centrality);
    rK0sResolution.fill(HIST("hEventOccupancy"), occupancy);

    rK0sResolution.fill(HIST("h1_stats"), 0.5);
    for (const auto& v0 : fullV0s) {
      bool daughtersHaveMCParticles = false;
      bool daughtersCorrRec = false;
      rK0sResolution.fill(HIST("h1_stats"), 1.5);
      const auto& posTrack = v0.posTrack_as<PIDTracksIUMC>();
      const auto& negTrack = v0.negTrack_as<PIDTracksIUMC>();
      if (!acceptV0(v0, negTrack, posTrack, collision))
        continue;
      rK0sResolution.fill(HIST("h1_stats"), 2.5);

      if (posTrack.has_mcParticle() && negTrack.has_mcParticle()) {
        daughtersHaveMCParticles = true;
        rK0sResolution.fill(HIST("h1_stats"), 3.5);
        bool isPositivePion = posTrack.mcParticle().pdgCode() == PDG_t::kPiPlus || (doTreatPiToMuon && posTrack.mcParticle().pdgCode() == PDG_t::kMuonPlus);
        bool isNegativePion = negTrack.mcParticle().pdgCode() == PDG_t::kPiMinus || (doTreatPiToMuon && negTrack.mcParticle().pdgCode() == PDG_t::kMuonMinus);
        if (isPositivePion && isNegativePion) {
          daughtersCorrRec = true;
          rK0sResolution.fill(HIST("h1_stats"), 4.5);
        }
      }

      if (useTrackTuner && daughtersHaveMCParticles) {
        tuneV0(v0, posTrack, negTrack, mcParticles, bcs);
      }

      float mass = v0.mK0Short();

      if (useTrackTuner && daughtersHaveMCParticles) {
        std::array<float, 3> pPos{0., 0., 0.};
        std::array<float, 3> pNeg{0., 0., 0.};
        mTrackParCovPos.getPxPyPzGlo(pPos);
        mTrackParCovNeg.getPxPyPzGlo(pNeg);
        mass = RecoDecay::m(std::array{std::array{pPos[0], pPos[1], pPos[2]},
                                       std::array{pNeg[0], pNeg[1], pNeg[2]}},
                            std::array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassPionCharged});
      }

      bool isTrueK0s = (v0.has_mcParticle() && std::abs(v0.mcParticle().pdgCode()) == PDG_t::kK0Short && v0.mcParticle().isPhysicalPrimary() && daughtersCorrRec);
      if (requireTrueK0s && !isTrueK0s) {
        continue;
      }

      // QA of correctly reconstructed V0 daughters
      if (daughtersCorrRec) {
        rK0sDauResolution.fill(HIST("h2_genPtPosPtResNorm"), (v0.positivept() - posTrack.mcParticle().pt()) / posTrack.mcParticle().pt(), posTrack.mcParticle().pt());
        rK0sDauResolution.fill(HIST("h2_genPxPosPxResNorm"), (v0.pxpos() - posTrack.mcParticle().px()) / posTrack.mcParticle().px(), posTrack.mcParticle().px());
        rK0sDauResolution.fill(HIST("h2_genPyPosPyResNorm"), (v0.pypos() - posTrack.mcParticle().py()) / posTrack.mcParticle().py(), posTrack.mcParticle().py());
        rK0sDauResolution.fill(HIST("h2_genPzPosPzResNorm"), (v0.pzpos() - posTrack.mcParticle().pz()) / posTrack.mcParticle().pz(), posTrack.mcParticle().pz());

        rK0sDauResolution.fill(HIST("h2_genPtNegPtResNorm"), (v0.negativept() - negTrack.mcParticle().pt()) / negTrack.mcParticle().pt(), negTrack.mcParticle().pt());
        rK0sDauResolution.fill(HIST("h2_genPxNegPxResNorm"), (v0.pxneg() - negTrack.mcParticle().px()) / negTrack.mcParticle().px(), negTrack.mcParticle().px());
        rK0sDauResolution.fill(HIST("h2_genPyNegPyResNorm"), (v0.pyneg() - negTrack.mcParticle().py()) / negTrack.mcParticle().py(), negTrack.mcParticle().py());
        rK0sDauResolution.fill(HIST("h2_genPzNegPzResNorm"), (v0.pzneg() - negTrack.mcParticle().pz()) / negTrack.mcParticle().pz(), negTrack.mcParticle().pz());

        rK0sDauResolution.fill(HIST("h2_genPtPosPtRes"), (v0.positivept() - posTrack.mcParticle().pt()), posTrack.mcParticle().pt());
        rK0sDauResolution.fill(HIST("h2_genPxPosPxRes"), (v0.pxpos() - posTrack.mcParticle().px()), posTrack.mcParticle().px());
        rK0sDauResolution.fill(HIST("h2_genPyPosPyRes"), (v0.pypos() - posTrack.mcParticle().py()), posTrack.mcParticle().py());
        rK0sDauResolution.fill(HIST("h2_genPzPosPzRes"), (v0.pzpos() - posTrack.mcParticle().pz()), posTrack.mcParticle().pz());

        rK0sDauResolution.fill(HIST("h2_genPtNegPtRes"), (v0.negativept() - negTrack.mcParticle().pt()), negTrack.mcParticle().pt());
        rK0sDauResolution.fill(HIST("h2_genPxNegPxRes"), (v0.pxneg() - negTrack.mcParticle().px()), negTrack.mcParticle().px());
        rK0sDauResolution.fill(HIST("h2_genPyNegPyRes"), (v0.pyneg() - negTrack.mcParticle().py()), negTrack.mcParticle().py());
        rK0sDauResolution.fill(HIST("h2_genPzNegPzRes"), (v0.pzneg() - negTrack.mcParticle().pz()), negTrack.mcParticle().pz());

        rK0sDauResolution.fill(HIST("h2_massPosPtRes"), mass, v0.positivept() - posTrack.mcParticle().pt());
        rK0sDauResolution.fill(HIST("h2_massNegPtRes"), mass, v0.negativept() - negTrack.mcParticle().pt());

        rK0sDauResolution.fill(HIST("h2_genPtPosPulls"), (1. / v0.positivept() - 1. / posTrack.mcParticle().pt()) / (RecoDecay::sqrtSumOfSquares(v0.covMatPosDau()[9], v0.covMatPosDau()[14]) / RecoDecay::sq(v0.positivept())), 1. / posTrack.mcParticle().pt());
        rK0sDauResolution.fill(HIST("h2_genPtNegPulls"), (1. / v0.negativept() - 1. / negTrack.mcParticle().pt()) / (RecoDecay::sqrtSumOfSquares(v0.covMatNegDau()[9], v0.covMatNegDau()[14]) / RecoDecay::sq(v0.negativept())), 1. / negTrack.mcParticle().pt());

        if (useMultidimHisto) {
          rK0sResolution.fill(HIST("thn_mass"), mass, v0.pt(), v0.eta(), v0.phi(), posTrack.eta(), negTrack.eta(),
                              1. / v0.positivept() - 1. / posTrack.mcParticle().pt(),
                              1. / v0.negativept() - 1. / negTrack.mcParticle().pt(),
                              isTrueK0s);
        }
      }

      // QA of seleted V0s
      rK0sDauResolution.fill(HIST("h2_PosRelPtRes"), v0.positivept(), RecoDecay::sqrtSumOfSquares(v0.covMatPosDau()[9], v0.covMatPosDau()[14]) / v0.positivept());
      rK0sDauResolution.fill(HIST("h2_NegRelPtRes"), v0.negativept(), RecoDecay::sqrtSumOfSquares(v0.covMatNegDau()[9], v0.covMatNegDau()[14]) / v0.negativept());
      rK0sResolution.fill(HIST("h2_masspT"), mass, v0.pt());
      rK0sResolution.fill(HIST("h2_masseta"), mass, v0.eta());
      rK0sResolution.fill(HIST("h2_massphi"), mass, v0.phi());
      rK0sResolution.fill(HIST("h3_centralitypTmass"), centrality, v0.pt(), mass);
      rK0sResolution.fill(HIST("h3_occupancypTmass"), occupancy, v0.pt(), mass);
    }
  }
  PROCESS_SWITCH(perfK0sResolution, processMC, "Process MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<perfK0sResolution>(cfgc)};
}
