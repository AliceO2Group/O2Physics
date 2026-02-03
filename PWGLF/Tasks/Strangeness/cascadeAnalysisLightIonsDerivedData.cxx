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
/// \file cascadeAnalysisLightIonsDerivedData.cxx
/// \brief analysis of cascades (Xi, antiXi, Omega, antiOmega) in light-ion collisions using derived data
///
/// \author Sara Pucillo (sara.pucillo@cern.ch), Alberto Caliva (alberto.caliva@cern.ch)
//

#include "PWGLF/DataModel/LFStrangenessPIDTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/CCDB/ctpRateFetcher.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CommonConstants/MathConstants.h"
#include "CommonConstants/PhysicsConstants.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

#include <TFile.h>
#include <TH2D.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TPDGCode.h>
#include <TProfile.h>

#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::math;
using namespace o2::constants::physics;
using std::array;

using SelCollisions = soa::Join<aod::StraCollisions, aod::StraEvSels, aod::StraCents, aod::StraStamps>;
using SimCollisions = soa::Join<aod::StraCollisions, aod::StraEvSels, aod::StraCents, aod::StraStamps, aod::StraCollLabels>;
using CascadeCandidates = soa::Join<aod::CascCollRefs, aod::CascCores, aod::CascExtras, aod::CascTOFPIDs, aod::CascTOFNSigmas>;
using CascadeMCCandidates = soa::Join<aod::CascCollRefs, aod::CascCores, aod::CascExtras, aod::CascTOFPIDs, aod::CascTOFNSigmas, aod::CascCoreMCLabels>;
using DaughterTracks = soa::Join<aod::DauTrackExtras, aod::DauTrackTPCPIDs>;
using CollisionMCTrueTable = soa::Join<aod::StraMCCollisions, aod::StraMCCollMults>;
using CascadeMCCores = soa::Join<aod::CascMCCores, aod::CascMCCollRefs>;

struct CascadeAnalysisLightIonsDerivedData {

  // Instantiate the CCDB service and API interface
  // CCDB options
  struct : ConfigurableGroup {
    std::string prefix = "ccdbConfigurations"; // JSON group name
    Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
    Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
    Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
    Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
    Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
    Configurable<std::string> mVtxPath{"mVtxPath", "GLO/Calib/MeanVertex", "Path of the mean vertex file"};
  } ccdbConfigurations;

  o2::ccdb::CcdbApi ccdbApi;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  int mRunNumber;

  // Define histogram registries
  HistogramRegistry registryData{"registryData", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry registryMC{"registryMC", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry registryQC{"registryQC", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  // Global analysis parameters
  Configurable<double> zVtx{"zVtx", 10.0, "Maximum zVertex"};
  Configurable<double> deltaEtaEdge{"deltaEtaEdge", 0.05, "eta gap from the edge"};
  Configurable<int> nBins{"nBins", 100, "N bins in all QC histos"};

  // Event selections parameters
  Configurable<bool> applySel8{"applySel8", true, "0 - no, 1 - yes"};
  Configurable<bool> applyVtxZ{"applyVtxZ", true, "0 - no, 1 - yes"};
  Configurable<bool> rejectITSROFBorder{"rejectITSROFBorder", true, "reject events at ITS ROF border"};
  Configurable<bool> rejectTFBorder{"rejectTFBorder", true, "reject events at TF border"};
  Configurable<bool> requireVertexITSTPC{"requireVertexITSTPC", false, "require events with at least one ITS-TPC track"};
  Configurable<bool> requireIsGoodZvtxFT0VsPV{"requireIsGoodZvtxFT0VsPV", true, "require is good Zvtx FT0 vs PV"};
  Configurable<bool> requireIsVertexTOFmatched{"requireIsVertexTOFmatched", false, "require events with at least one of vertex contributors matched to TOF"};
  Configurable<bool> requireIsVertexTRDmatched{"requireIsVertexTRDmatched", false, "require events with at least one of vertex contributors matched to TRD"};
  Configurable<bool> rejectSameBunchPileup{"rejectSameBunchPileup", true, "reject collisions in case of pileup with another collision in the same foundBC"};
  Configurable<bool> applyZVtxSelOnMCPV{"applyZVtxSelOnMCPV", false, "Apply Z-vtx cut on the PV of the generated collision?"};
  Configurable<bool> requireInel0{"requireInel0", false, "Enable INEL > 0 selection"};
  Configurable<bool> requireInel0OnMC{"requireInel0OnMC", false, "Enable INEL > 0 selection for MC gen events"};

  // Track analysis Parameters
  Configurable<float> minITSnCls{"minITSnCls", 4.0f, "min number of ITS clusters"};
  Configurable<float> minTPCnClsFound{"minTPCnClsFound", 80.0f, "min number of found TPC clusters"};
  Configurable<float> minNCrossedRowsTPC{"minNCrossedRowsTPC", 80.0f, "min number of TPC crossed rows"};
  Configurable<float> etaMin{"etaMin", -0.8f, "eta min"};
  Configurable<float> etaMax{"etaMax", +0.8f, "eta max"};
  Configurable<float> rapcut{"rapcut", +0.5f, "rapidity cut"};
  Configurable<float> maxChi2TPC{"maxChi2TPC", 4.0f, "max chi2 per cluster TPC"};
  Configurable<bool> requireITS{"requireITS", false, "require ITS hit"};
  Configurable<bool> requireTOF{"requireTOF", false, "require TOF hit"};

  // Configurable parameters for PID selection
  Configurable<float> nsigmaTPCmin{"nsigmaTPCmin", -5.0f, "Minimum nsigma TPC"};
  Configurable<float> nsigmaTPCmax{"nsigmaTPCmax", +5.0f, "Maximum nsigma TPC"};

  // Configurable parameters for TOF PID selection
  Configurable<float> nsigmaTOFmin{"nsigmaTOFmin", -3.0f, "Minimum nsigma TOF"};
  Configurable<float> nsigmaTOFmax{"nsigmaTOFmax", +3.0f, "Maximum nsigma TOF"};

  // Topological Parameters
  Configurable<float> dcanegtoPVmin{"dcanegtoPVmin", 0.1f, "Minimum DCA Neg To PV"};
  Configurable<float> dcapostoPVmin{"dcapostoPVmin", 0.1f, "Minimum DCA Pos To PV"};
  Configurable<float> dcabachtopvMin{"dcabachtopvMin", 0.1f, "Minimum DCA bachelor to PV"};
  Configurable<float> dcaV0DaughtersMax{"dcaV0DaughtersMax", 0.7f, "Maximum DCA Daughters"};
  Configurable<float> dcaV0topvMin{"dcaV0topvMin", 0.02f, "Minimum DCA V0 to PV"};
  Configurable<float> dcaCascDaughtersMax{"dcaCascDaughtersMax", 0.8f, "Maximum DCA Daughters"};
  Configurable<float> v0cospaMin{"v0cospaMin", 0.99f, "Minimum V0 CosPA"};
  Configurable<float> casccospaMin{"casccospaMin", 0.99f, "Minimum Cascade CosPA"};
  Configurable<float> minimumV0Radius{"minimumV0Radius", 2.5f, "Minimum V0 Radius"};
  Configurable<float> minimumCascRadius{"minimumCascRadius", 1.1f, "Minimum Cascade Radius"};
  Configurable<float> v0masswindow{"v0masswindow", 0.005, "v0 mass window"};
  Configurable<float> competingmassrej{"competingmassrej", 0.008, "Competing mass rejection"};

  // Axes parameters
  ConfigurableAxis centEstimatorHistBin{"centEstimatorHistBin", {101, 0.0f, 101.0f}, ""};
  ConfigurableAxis centralityBinning{"centralityBinning", {VARIABLE_WIDTH, 0.0f, 5.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f}, ""};
  ConfigurableAxis axisNch{"axisNch", {500, 0.0f, +5000.0f}, "Number of charged particles"};
  ConfigurableAxis axisMult{"axisMult", {500, 0.0f, +100000.0f}, "Multiplicity"};

  // Centrality estimator
  Configurable<int> centralityEstimator{"centralityEstimator", 0, "0 = FT0C, 1 = FTOM, 2 = FV0A, 3 = NGlobal"};

  // List of estimators
  enum Option { kFT0C,
                kFT0M,
                kFV0A,
                kNGlobal };

  // For manual sliceBy
  PresliceUnsorted<soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels, aod::StraCollLabels>> perMcCollision = aod::cascdata::straMCCollisionId;

  void init(InitContext const&)
  {
    // setting CCDB service
    ccdb->setURL(ccdbConfigurations.ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setFatalWhenNull(false);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setCreatedNotAfter(
      std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::system_clock::now().time_since_epoch())
        .count());

    // Axes and Binning
    AxisSpec axisCentEstimator = {centEstimatorHistBin, "CentEstimator", "CentEstimatorAxis"};
    AxisSpec centAxis = {centralityBinning, "Centrality", "CentralityAxis"};
    const AxisSpec vertexZAxis = {nBins, -15., 15., "vrtx_{Z} [cm]"};
    const AxisSpec etaAxis{18, -0.9, 0.9, "#eta"};
    const AxisSpec ptAxis{100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec invMassXiAxis{200, 1.28, 1.36, "m_{p#pi#pi} (GeV/#it{c}^{2})"};
    const AxisSpec invMassOmegaAxis{200, 1.63, 1.71, "m_{p#piK} (GeV/#it{c}^{2})"};
    const AxisSpec nsigmaTOFAxis{200, -10, 10, "n#sigma_{TOF}"};
    const AxisSpec nsigmaTPCAxis{200, -10, 10, "n#sigma_{TPC}"};

    // Histograms for data
    if (doprocessData) {
      registryData.add("number_of_events_data", "number of events in data", HistType::kTH1D, {{20, -0.5f, +19.5f}});
      registryData.get<TH1>(HIST("number_of_events_data"))->GetXaxis()->SetBinLabel(1, "All collisions");
      registryData.get<TH1>(HIST("number_of_events_data"))->GetXaxis()->SetBinLabel(2, "sel8 cut");
      registryData.get<TH1>(HIST("number_of_events_data"))->GetXaxis()->SetBinLabel(3, "posZ cut");
      registryData.get<TH1>(HIST("number_of_events_data"))->GetXaxis()->SetBinLabel(4, "kNoITSROFrameBorder");
      registryData.get<TH1>(HIST("number_of_events_data"))->GetXaxis()->SetBinLabel(5, "kNoTimeFrameBorder");
      registryData.get<TH1>(HIST("number_of_events_data"))->GetXaxis()->SetBinLabel(6, "kIsVertexITSTPC");
      registryData.get<TH1>(HIST("number_of_events_data"))->GetXaxis()->SetBinLabel(7, "kIsGoodZvtxFT0vsPV");
      registryData.get<TH1>(HIST("number_of_events_data"))->GetXaxis()->SetBinLabel(8, "kIsVertexTOFmatched");
      registryData.get<TH1>(HIST("number_of_events_data"))->GetXaxis()->SetBinLabel(9, "kIsVertexTRDmatched");
      registryData.get<TH1>(HIST("number_of_events_data"))->GetXaxis()->SetBinLabel(10, "kNoSameBunchPileup");
      registryData.get<TH1>(HIST("number_of_events_data"))->GetXaxis()->SetBinLabel(11, "kINELgr0");

      registryData.add("number_of_events_data_vs_centrality", "number of events in data vs centrality", HistType::kTH2D, {{20, -0.5f, +19.5f}, {101, 0.0f, 101.0f}});
      registryData.get<TH2>(HIST("number_of_events_data_vs_centrality"))->GetXaxis()->SetBinLabel(1, "All collisions");
      registryData.get<TH2>(HIST("number_of_events_data_vs_centrality"))->GetXaxis()->SetBinLabel(2, "sel8 cut");
      registryData.get<TH2>(HIST("number_of_events_data_vs_centrality"))->GetXaxis()->SetBinLabel(3, "posZ cut");
      registryData.get<TH2>(HIST("number_of_events_data_vs_centrality"))->GetXaxis()->SetBinLabel(4, "kNoITSROFrameBorder");
      registryData.get<TH2>(HIST("number_of_events_data_vs_centrality"))->GetXaxis()->SetBinLabel(5, "kNoTimeFrameBorder");
      registryData.get<TH2>(HIST("number_of_events_data_vs_centrality"))->GetXaxis()->SetBinLabel(6, "kIsVertexITSTPC");
      registryData.get<TH2>(HIST("number_of_events_data_vs_centrality"))->GetXaxis()->SetBinLabel(7, "kIsGoodZvtxFT0vsPV");
      registryData.get<TH2>(HIST("number_of_events_data_vs_centrality"))->GetXaxis()->SetBinLabel(8, "kIsVertexTOFmatched");
      registryData.get<TH2>(HIST("number_of_events_data_vs_centrality"))->GetXaxis()->SetBinLabel(9, "kIsVertexTRDmatched");
      registryData.get<TH2>(HIST("number_of_events_data_vs_centrality"))->GetXaxis()->SetBinLabel(10, "kNoSameBunchPileup");
      registryData.get<TH2>(HIST("number_of_events_data_vs_centrality"))->GetXaxis()->SetBinLabel(11, "kINELgr0");
      registryData.get<TH2>(HIST("number_of_events_data_vs_centrality"))->GetYaxis()->SetTitle("Centrality (%)");

      // QC Histograms
      registryQC.add("hVertexZdata", "hVertexZdata", HistType::kTH1D, {vertexZAxis});
      registryQC.add("hv0cosPAdata", "hv0cosPAdata", HistType::kTH1D, {{nBins, 0.95f, 1.f}});
      registryQC.add("hcasccosPAdata", "hcasccosPAdata", HistType::kTH1D, {{nBins, 0.95f, 1.f}});
      registryQC.add("hv0radiusdata", "hv0radiusdata", HistType::kTH1D, {{nBins, 0.0f, 5.0f}});
      registryQC.add("hcascradiusdata", "hcascradiusdata", HistType::kTH1D, {{nBins, 0.0f, 5.0f}});
      registryQC.add("hdcaV0daughtersdata", "hdcaV0daughtersdata", HistType::kTH1D, {{nBins, 0.0f, 1.5f}});
      registryQC.add("hdcacascdaughtersdata", "hdcacascdaughtersdata", HistType::kTH1D, {{nBins, 0.0f, 1.5f}});
      registryQC.add("hdcapostopvdata", "hdcapostopvdata", HistType::kTH1D, {{nBins, 0.0f, 2.0f}});
      registryQC.add("hdcanegtopvdata", "hdcanegtopvdata", HistType::kTH1D, {{nBins, 0.0f, 2.0f}});
      registryQC.add("hdcabachtopvdata", "hdcabachtopvdata", HistType::kTH1D, {{nBins, 0.0f, 2.0f}});
      registryQC.add("hdcav0topvdata", "hdcav0topvdata", HistType::kTH1D, {{nBins, 0.0f, 2.0f}});

      // Multiplicity Histograms
      registryData.add("hCentEstimator", "hCentEstimator", HistType::kTH1D, {{101, 0.0f, 101.0f}});
      registryData.add("hCentralityVsNch", "hCentralityVsNch", HistType::kTH2D, {{101, 0.0f, 101.0f}, axisNch});
      if (centralityEstimator == Option::kNGlobal)
        registryData.add("hCentralityVsMultiplicity", "hCentralityVsMultiplicity", HistType::kTH2D, {{101, 0.0f, 101.0f}, axisNch});
      else
        registryData.add("hCentralityVsMultiplicity", "hCentralityVsMultiplicity", HistType::kTH2D, {{101, 0.0f, 101.0f}, axisMult});

      // Histograms for xi (data)
      registryData.add("hMassXipos", "hMassXipos", HistType::kTH3D, {centAxis, ptAxis, invMassXiAxis});
      registryData.add("hMassXineg", "hMassXineg", HistType::kTH3D, {centAxis, ptAxis, invMassXiAxis});
      registryData.add("hMassXiposSelected", "hMassXiposSelected", HistType::kTH3D, {centAxis, ptAxis, invMassXiAxis});
      registryData.add("hMassXinegSelected", "hMassXinegSelected", HistType::kTH3D, {centAxis, ptAxis, invMassXiAxis});

      // Histograms for omega (data)
      registryData.add("hMassOmegapos", "hMassOmegapos", HistType::kTH3D, {centAxis, ptAxis, invMassOmegaAxis});
      registryData.add("hMassOmeganeg", "hMassOmeganeg", HistType::kTH3D, {centAxis, ptAxis, invMassOmegaAxis});
      registryData.add("hMassOmegaposSelected", "hMassOmegaposSelected", HistType::kTH3D, {centAxis, ptAxis, invMassOmegaAxis});
      registryData.add("hMassOmeganegSelected", "hMassOmeganegSelected", HistType::kTH3D, {centAxis, ptAxis, invMassOmegaAxis});
    }

    if (doprocessMonteCarloRec) {
      // Histograms for mc reconstructed
      registryMC.add("number_of_events_mc_rec", "number of events in mc_rec", HistType::kTH1D, {{20, -0.5f, +19.5f}});
      registryMC.get<TH1>(HIST("number_of_events_mc_rec"))->GetXaxis()->SetBinLabel(1, "All collisions");
      registryMC.get<TH1>(HIST("number_of_events_mc_rec"))->GetXaxis()->SetBinLabel(2, "sel8 cut");
      registryMC.get<TH1>(HIST("number_of_events_mc_rec"))->GetXaxis()->SetBinLabel(3, "posZ cut");
      registryMC.get<TH1>(HIST("number_of_events_mc_rec"))->GetXaxis()->SetBinLabel(4, "kNoITSROFrameBorder");
      registryMC.get<TH1>(HIST("number_of_events_mc_rec"))->GetXaxis()->SetBinLabel(5, "kNoTimeFrameBorder");
      registryMC.get<TH1>(HIST("number_of_events_mc_rec"))->GetXaxis()->SetBinLabel(6, "kIsVertexITSTPC");
      registryMC.get<TH1>(HIST("number_of_events_mc_rec"))->GetXaxis()->SetBinLabel(7, "kIsGoodZvtxFT0vsPV");
      registryMC.get<TH1>(HIST("number_of_events_mc_rec"))->GetXaxis()->SetBinLabel(8, "kIsVertexTOFmatched");
      registryMC.get<TH1>(HIST("number_of_events_mc_rec"))->GetXaxis()->SetBinLabel(9, "kIsVertexTRDmatched");
      registryMC.get<TH1>(HIST("number_of_events_mc_rec"))->GetXaxis()->SetBinLabel(10, "kNoSameBunchPileup");
      registryMC.get<TH1>(HIST("number_of_events_mc_rec"))->GetXaxis()->SetBinLabel(11, "kINELgr0");

      registryMC.add("number_of_events_mc_rec_vs_centrality", "number of events in mc_rec vs centrality", HistType::kTH2D, {{20, -0.5f, +19.5f}, {101, 0.0f, 101.0f}});
      registryMC.get<TH2>(HIST("number_of_events_mc_rec_vs_centrality"))->GetXaxis()->SetBinLabel(1, "All collisions");
      registryMC.get<TH2>(HIST("number_of_events_mc_rec_vs_centrality"))->GetXaxis()->SetBinLabel(2, "sel8 cut");
      registryMC.get<TH2>(HIST("number_of_events_mc_rec_vs_centrality"))->GetXaxis()->SetBinLabel(3, "posZ cut");
      registryMC.get<TH2>(HIST("number_of_events_mc_rec_vs_centrality"))->GetXaxis()->SetBinLabel(4, "kNoITSROFrameBorder");
      registryMC.get<TH2>(HIST("number_of_events_mc_rec_vs_centrality"))->GetXaxis()->SetBinLabel(5, "kNoTimeFrameBorder");
      registryMC.get<TH2>(HIST("number_of_events_mc_rec_vs_centrality"))->GetXaxis()->SetBinLabel(6, "kIsVertexITSTPC");
      registryMC.get<TH2>(HIST("number_of_events_mc_rec_vs_centrality"))->GetXaxis()->SetBinLabel(7, "kIsGoodZvtxFT0vsPV");
      registryMC.get<TH2>(HIST("number_of_events_mc_rec_vs_centrality"))->GetXaxis()->SetBinLabel(8, "kIsVertexTOFmatched");
      registryMC.get<TH2>(HIST("number_of_events_mc_rec_vs_centrality"))->GetXaxis()->SetBinLabel(9, "kIsVertexTRDmatched");
      registryMC.get<TH2>(HIST("number_of_events_mc_rec_vs_centrality"))->GetXaxis()->SetBinLabel(10, "kNoSameBunchPileup");
      registryMC.get<TH2>(HIST("number_of_events_mc_rec_vs_centrality"))->GetXaxis()->SetBinLabel(11, "kINELgr0");
      registryMC.get<TH2>(HIST("number_of_events_mc_rec_vs_centrality"))->GetYaxis()->SetTitle("Centrality (%)");

      // QC Histograms
      registryQC.add("hVertexZRec", "hVertexZRec", HistType::kTH1D, {{vertexZAxis}});
      registryQC.add("hv0cosPARec", "hv0cosPARec", HistType::kTH1D, {{nBins, 0.95f, 1.f}});
      registryQC.add("hcasccosPARec", "hcasccosPARec", HistType::kTH1D, {{nBins, 0.95f, 1.f}});
      registryQC.add("hv0radiusRec", "hv0radiusRec", HistType::kTH1D, {{nBins, 0.0f, 5.0f}});
      registryQC.add("hcascradiusRec", "hcascradiusRec", HistType::kTH1D, {{nBins, 0.0f, 5.0f}});
      registryQC.add("hdcaV0daughtersRec", "hdcaV0daughtersRec", HistType::kTH1D, {{nBins, 0.0f, 1.5f}});
      registryQC.add("hdcacascdaughtersRec", "hdcacascdaughtersRec", HistType::kTH1D, {{nBins, 0.0f, 1.5f}});
      registryQC.add("hdcapostopvRec", "hdcapostopvRec", HistType::kTH1D, {{nBins, 0.0f, 2.0f}});
      registryQC.add("hdcanegtopvRec", "hdcanegtopvRec", HistType::kTH1D, {{nBins, 0.0f, 2.0f}});
      registryQC.add("hdcabachtopvRec", "hdcabachtopvRec", HistType::kTH1D, {{nBins, 0.0f, 2.0f}});
      registryQC.add("hdcav0topvRec", "hdcav0topvRec", HistType::kTH1D, {{nBins, 0.0f, 2.0f}});

      // Multiplicity Histograms
      registryMC.add("hCentEstimator_truerec", "hCentEstimator_truerec", HistType::kTH1D, {{101, 0.0f, 101.0f}});
      registryMC.add("hCentralityVsNch_truerec", "hCentralityVsNch_truerec", HistType::kTH2D, {{101, 0.0f, 101.0f}, axisNch});
      if (centralityEstimator == Option::kNGlobal)
        registryMC.add("hCentralityVsMultiplicity_truerec", "hCentralityVsMultiplicity_truerec", HistType::kTH2D, {{101, 0.0f, 101.0f}, axisNch});
      else
        registryMC.add("hCentralityVsMultiplicity_truerec", "hCentralityVsMultiplicity_truerec", HistType::kTH2D, {{101, 0.0f, 101.0f}, axisMult});

      // Histograms for xi (mc)
      registryMC.add("hMassXipos_truerec", "hMassXipos_truerec", HistType::kTH3D, {centAxis, ptAxis, invMassXiAxis});
      registryMC.add("hMassXineg_truerec", "hMassXineg_truerec", HistType::kTH3D, {centAxis, ptAxis, invMassXiAxis});
      registryMC.add("hMassXiposSelected_truerec", "hMassXiposSelected_truerec", HistType::kTH3D, {centAxis, ptAxis, invMassXiAxis});
      registryMC.add("hMassXinegSelected_truerec", "hMassXinegSelected_truerec", HistType::kTH3D, {centAxis, ptAxis, invMassXiAxis});

      // Histograms for omega (mc)
      registryMC.add("hMassOmegapos_truerec", "hMassOmegapos_truerec", HistType::kTH3D, {centAxis, ptAxis, invMassOmegaAxis});
      registryMC.add("hMassOmeganeg_truerec", "hMassOmeganeg_truerec", HistType::kTH3D, {centAxis, ptAxis, invMassOmegaAxis});
      registryMC.add("hMassOmegaposSelected_truerec", "hMassOmegaposSelected_truerec", HistType::kTH3D, {centAxis, ptAxis, invMassOmegaAxis});
      registryMC.add("hMassOmeganegSelected_truerec", "hMassOmeganegSelected_truerec", HistType::kTH3D, {centAxis, ptAxis, invMassOmegaAxis});
    }

    if (doprocessMonteCarloGen) {
      // Histograms for mc generated
      // QC Histograms
      registryQC.add("hVertexZGen", "hVertexZGen", HistType::kTH1D, {{vertexZAxis}});
      // Histograms for xi (mc)
      registryMC.add("h2dGenXiMinusVsMultMC_RecoedEvt", "h2dGenXiMinusVsMultMC_RecoedEvt", HistType::kTH2D, {axisNch, ptAxis});
      registryMC.add("h2dGenXiPlusVsMultMC_RecoedEvt", "h2dGenXiPlusVsMultMC_RecoedEvt", HistType::kTH2D, {axisNch, ptAxis});
      registryMC.add("h2dGenXiMinusVsMultMC", "h2dGenXiMinusVsMultMC", HistType::kTH2D, {axisNch, ptAxis});
      registryMC.add("h2dGenXiPlusVsMultMC", "h2dGenXiPlusVsMultMC", HistType::kTH2D, {axisNch, ptAxis});
      registryMC.add("h2dGenXiMinus", "h2dGenXiMinus", HistType::kTH2D, {centAxis, ptAxis});
      registryMC.add("h2dGenXiPlus", "h2dGenXiPlus", HistType::kTH2D, {centAxis, ptAxis});

      // Histograms for omega (mc)
      registryMC.add("h2dGenOmegaMinusVsMultMC_RecoedEvt", "h2dGenOmegaMinusVsMultMC_RecoedEvt", HistType::kTH2D, {axisNch, ptAxis});
      registryMC.add("h2dGenOmegaPlusVsMultMC_RecoedEvt", "h2dGenOmegaPlusVsMultMC_RecoedEvt", HistType::kTH2D, {axisNch, ptAxis});
      registryMC.add("h2dGenOmegaMinusVsMultMC", "h2dGenOmegaMinusVsMultMC", HistType::kTH2D, {axisNch, ptAxis});
      registryMC.add("h2dGenOmegaPlusVsMultMC", "h2dGenOmegaPlusVsMultMC", HistType::kTH2D, {axisNch, ptAxis});
      registryMC.add("h2dGenOmegaMinus", "h2dGenOmegaMinus", HistType::kTH2D, {centAxis, ptAxis});
      registryMC.add("h2dGenOmegaPlus", "h2dGenOmegaPlus", HistType::kTH2D, {centAxis, ptAxis});

      // Histograms for event loss/splitting
      registryMC.add("hGenEvents", "hGenEvents", HistType::kTH2D, {{axisNch}, {4, -0.5f, +3.5f}});
      registryMC.get<TH2>(HIST("hGenEvents"))->GetYaxis()->SetBinLabel(1, "All gen. events");
      registryMC.get<TH2>(HIST("hGenEvents"))->GetYaxis()->SetBinLabel(2, "All gen. events in INEL > 0");
      registryMC.get<TH2>(HIST("hGenEvents"))->GetYaxis()->SetBinLabel(3, "Gen. with at least 1 rec. events");
      registryMC.get<TH2>(HIST("hGenEvents"))->GetYaxis()->SetBinLabel(4, "Gen. with at least 1 rec. events in INEL > 0");
      registryMC.add("hGenEventCentrality", "hGenEventCentrality", kTH1D, {{101, 0.0f, 101.0f}});

      registryMC.add("hCentralityVsNcoll_beforeEvSel", "hCentralityVsNcoll_beforeEvSel", HistType::kTH2D, {centAxis, {50, -0.5f, 49.5f}});
      registryMC.add("hCentralityVsNcoll_afterEvSel", "hCentralityVsNcoll_afterEvSel", HistType::kTH2D, {centAxis, {50, -0.5f, 49.5f}});

      registryMC.add("hCentralityVsMultMC", "hCentralityVsMultMC", HistType::kTH2D, {{101, 0.0f, 101.0f}, axisNch});
    }
  }

  // ______________________________________________________
  // Return slicing output
  template <typename TCollisions>
  auto getGroupedCollisions(TCollisions const& collisions, int globalIndex)
  {
    return collisions.sliceBy(perMcCollision, globalIndex);
  }

  template <typename TCollision>
  void initCCDB(TCollision collision)
  {
    if (mRunNumber == collision.runNumber()) {
      return;
    }

    mRunNumber = collision.runNumber();
  }

  // Find ITS hit
  // template <typename TrackIts>
  // bool hasITSHitOnLayer(const TrackIts& track, int layer)
  // {
  //   int ibit = layer - 1;
  //   return (track.itsClusterMap() & (1 << ibit));
  // }

  // Single-Track Selection
  template <typename Track>
  bool passedSingleTrackSelection(const Track& track)
  {
    if (requireITS && (!track.hasITS()))
      return false;
    if (requireITS && track.itsNCls() < minITSnCls)
      return false;
    if (!track.hasTPC())
      return false;
    if (track.tpcNClsFound() < minTPCnClsFound)
      return false;
    if (track.tpcNClsCrossedRows() < minNCrossedRowsTPC)
      return false;
    if (track.tpcChi2NCl() > maxChi2TPC)
      return false;
    if (requireTOF && (!track.hasTOF()))
      return false;
    return true;
  }

  // Xi Selection
  template <typename Xi, typename TrackPos, typename TrackNeg,
            typename TrackBac, typename Coll>
  bool passedXiSelection(const Xi& casc, const TrackPos& ptrack,
                         const TrackNeg& ntrack, const TrackBac& btrack,
                         const Coll& coll)
  {
    if (!passedSingleTrackSelection(ptrack))
      return false;
    if (!passedSingleTrackSelection(ntrack))
      return false;
    if (!passedSingleTrackSelection(btrack))
      return false;

    if (std::abs(casc.yXi()) > rapcut)
      return false;

    // Xi+ Selection (Xi+ -> antiL + pi+)
    if (casc.sign() > 0) {

      // PID Selections (TPC)
      if (ntrack.tpcNSigmaPr() < nsigmaTPCmin ||
          ntrack.tpcNSigmaPr() > nsigmaTPCmax)
        return false;
      if (ptrack.tpcNSigmaPi() < nsigmaTPCmin ||
          ptrack.tpcNSigmaPi() > nsigmaTPCmax)
        return false;

      // PID Selections (TOF)
      if (requireTOF) {
        if (casc.tofNSigmaXiLaPr() < nsigmaTOFmin ||
            casc.tofNSigmaXiLaPr() > nsigmaTOFmax)
          return false;
        if (casc.tofNSigmaXiLaPi() < nsigmaTOFmin ||
            casc.tofNSigmaXiLaPi() > nsigmaTOFmax)
          return false;
      }
    }

    // Xi- Selection (Xi- -> L + pi-)
    if (casc.sign() < 0) {

      // PID Selections (TPC)
      if (ptrack.tpcNSigmaPr() < nsigmaTPCmin ||
          ptrack.tpcNSigmaPr() > nsigmaTPCmax)
        return false;
      if (ntrack.tpcNSigmaPi() < nsigmaTPCmin ||
          ntrack.tpcNSigmaPi() > nsigmaTPCmax)
        return false;

      // PID Selections (TOF)
      if (requireTOF) {
        if (casc.tofNSigmaXiLaPr() < nsigmaTOFmin ||
            casc.tofNSigmaXiLaPr() > nsigmaTOFmax)
          return false;
        if (casc.tofNSigmaXiLaPi() < nsigmaTOFmin ||
            casc.tofNSigmaXiLaPi() > nsigmaTOFmax)
          return false;
      }
    }

    if (casc.v0cosPA(coll.posX(), coll.posY(), coll.posZ()) < v0cospaMin)
      return false;
    if (casc.v0radius() < minimumV0Radius)
      return false;
    if (std::fabs(casc.dcaV0daughters()) > dcaV0DaughtersMax)
      return false;
    if (std::fabs(casc.dcapostopv()) < dcapostoPVmin)
      return false;
    if (std::fabs(casc.dcanegtopv()) < dcanegtoPVmin)
      return false;

    if (casc.casccosPA(coll.posX(), coll.posY(), coll.posZ()) < casccospaMin)
      return false;
    if (casc.cascradius() < minimumCascRadius)
      return false;
    if (std::fabs(casc.dcabachtopv()) < dcabachtopvMin)
      return false;
    if (std::fabs(casc.dcav0topv(coll.posX(), coll.posY(), coll.posZ())) <
        dcaV0topvMin)
      return false;
    if (std::fabs(casc.dcacascdaughters()) > dcaCascDaughtersMax)
      return false;

    // V0 mass window
    if (std::abs(casc.mLambda() - o2::constants::physics::MassLambda0) >
        v0masswindow)
      return false;

    // reject candidates compatible with omega
    if (std::abs(casc.mOmega() - o2::constants::physics::MassOmegaMinus) <
        competingmassrej)
      return false;

    // PID Selection on bachelor
    if (btrack.tpcNSigmaPi() < nsigmaTPCmin ||
        btrack.tpcNSigmaPi() > nsigmaTPCmax)
      return false;

    // PID Selections (TOF)
    if (requireTOF) {
      if (casc.tofNSigmaXiPi() < nsigmaTOFmin ||
          casc.tofNSigmaXiPi() > nsigmaTOFmax)
        return false;
    }
    return true;
  }

  // Omega Selection
  template <typename Omega, typename TrackPos, typename TrackNeg,
            typename TrackBac, typename Coll>
  bool passedOmegaSelection(const Omega& casc, const TrackPos& ptrack,
                            const TrackNeg& ntrack, const TrackBac& btrack,
                            const Coll& coll)
  {
    if (!passedSingleTrackSelection(ptrack))
      return false;
    if (!passedSingleTrackSelection(ntrack))
      return false;
    if (!passedSingleTrackSelection(btrack))
      return false;

    if (std::abs(casc.yOmega()) > rapcut)
      return false;

    // Omega+ Selection (Omega+ -> antiL + K+)
    if (casc.sign() > 0) {
      // PID Selections (TPC)
      if (ntrack.tpcNSigmaPr() < nsigmaTPCmin ||
          ntrack.tpcNSigmaPr() > nsigmaTPCmax)
        return false;
      if (ptrack.tpcNSigmaPi() < nsigmaTPCmin ||
          ptrack.tpcNSigmaPi() > nsigmaTPCmax)
        return false;

      // PID Selections (TOF)
      if (requireTOF) {
        if (casc.tofNSigmaOmLaPr() < nsigmaTOFmin ||
            casc.tofNSigmaOmLaPr() > nsigmaTOFmax)
          return false;
        if (casc.tofNSigmaOmLaPi() < nsigmaTOFmin ||
            casc.tofNSigmaOmLaPi() > nsigmaTOFmax)
          return false;
      }
    }

    // Omega- Selection (Omega- -> L + K-)
    if (casc.sign() < 0) {
      // PID Selections (TPC)
      if (ptrack.tpcNSigmaPr() < nsigmaTPCmin ||
          ptrack.tpcNSigmaPr() > nsigmaTPCmax)
        return false;
      if (ntrack.tpcNSigmaPi() < nsigmaTPCmin ||
          ntrack.tpcNSigmaPi() > nsigmaTPCmax)
        return false;

      // PID Selections (TOF)
      if (requireTOF) {
        if (casc.tofNSigmaOmLaPr() < nsigmaTOFmin ||
            casc.tofNSigmaOmLaPr() > nsigmaTOFmax)
          return false;
        if (casc.tofNSigmaOmLaPi() < nsigmaTOFmin ||
            casc.tofNSigmaOmLaPi() > nsigmaTOFmax)
          return false;
      }
    }

    if (casc.v0cosPA(coll.posX(), coll.posY(), coll.posZ()) < v0cospaMin)
      return false;
    if (casc.v0radius() < minimumV0Radius)
      return false;
    if (std::fabs(casc.dcaV0daughters()) > dcaV0DaughtersMax)
      return false;
    if (std::fabs(casc.dcapostopv()) < dcapostoPVmin)
      return false;
    if (std::fabs(casc.dcanegtopv()) < dcanegtoPVmin)
      return false;

    if (casc.casccosPA(coll.posX(), coll.posY(), coll.posZ()) < casccospaMin)
      return false;
    if (casc.cascradius() < minimumCascRadius)
      return false;
    if (std::fabs(casc.dcabachtopv()) < dcabachtopvMin)
      return false;
    if (std::fabs(casc.dcav0topv(coll.posX(), coll.posY(), coll.posZ())) <
        dcaV0topvMin)
      return false;
    if (std::fabs(casc.dcacascdaughters()) > dcaCascDaughtersMax)
      return false;

    // V0 mass window
    if (std::abs(casc.mLambda() - o2::constants::physics::MassLambda0) >
        v0masswindow)
      return false;

    // Reject candidates compatible with Xi
    if (std::abs(casc.mXi() - o2::constants::physics::MassXiMinus) <
        competingmassrej)
      return false;

    // PID Selection on bachelor
    if (btrack.tpcNSigmaKa() < nsigmaTPCmin ||
        btrack.tpcNSigmaKa() > nsigmaTPCmax)
      return false;

    // PID Selections (TOF)
    if (requireTOF) {
      if (casc.tofNSigmaOmKa() < nsigmaTOFmin ||
          casc.tofNSigmaOmKa() > nsigmaTOFmax)
        return false;
    }
    return true;
  }

  // Return the list of indices to the recoed collision associated to a given MC collision.
  template <typename TMCollisions, typename TCollisions>
  std::vector<int> getListOfRecoCollIndices(TMCollisions const& mcCollisions, TCollisions const& collisions)
  {
    std::vector<int> listBestCollisionIdx(mcCollisions.size());
    for (auto const& mcCollision : mcCollisions) {
      auto groupedCollisions = getGroupedCollisions(collisions, mcCollision.globalIndex());
      int biggestNContribs = -1;
      int bestCollisionIndex = -1;
      for (auto const& collision : groupedCollisions) {
        // event selections
        if (applySel8 && !collision.sel8())
          continue;

        if (applyVtxZ && std::fabs(collision.posZ()) > zVtx)
          continue;

        if (rejectITSROFBorder && !collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
          continue;
        }

        if (rejectTFBorder && !collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
          continue;
        }

        if (requireVertexITSTPC && !collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) {
          continue;
        }

        if (requireIsGoodZvtxFT0VsPV && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
          continue;
        }

        if (requireIsVertexTOFmatched && !collision.selection_bit(o2::aod::evsel::kIsVertexTOFmatched)) {
          continue;
        }

        if (requireIsVertexTRDmatched && !collision.selection_bit(o2::aod::evsel::kIsVertexTRDmatched)) {
          continue;
        }

        if (rejectSameBunchPileup && !collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
          continue;
        }

        if (requireInel0 && collision.multNTracksPVeta1() < 1) {
          continue;
        }

        // Find the collision with the biggest nbr of PV contributors
        // Follows what was done here: https://github.com/AliceO2Group/O2Physics/blob/master/Common/TableProducer/mcCollsExtra.cxx#L93
        if (biggestNContribs < collision.multPVTotalContributors()) {
          biggestNContribs = collision.multPVTotalContributors();
          bestCollisionIndex = collision.globalIndex();
        }
      }
      listBestCollisionIdx[mcCollision.globalIndex()] = bestCollisionIndex;
    }
    return listBestCollisionIdx;
  }

  // Fill generated event information (for event loss/splitting estimation)
  template <typename TMCCollisions, typename TCollisions>
  void fillGeneratedEventProperties(TMCCollisions const& mcCollisions, TCollisions const& collisions)
  {
    std::vector<int> listBestCollisionIdx(mcCollisions.size());
    for (auto const& mcCollision : mcCollisions) {
      // event selections
      if (applyZVtxSelOnMCPV && std::fabs(mcCollision.posZ()) > zVtx)
        continue;

      registryMC.fill(HIST("hGenEvents"), mcCollision.multMCNParticlesEta05(), 0 /* all gen. events*/);

      if (requireInel0OnMC && mcCollision.multMCNParticlesEta10() < 1) {
        continue;
      }
      registryMC.fill(HIST("hGenEvents"), mcCollision.multMCNParticlesEta05(), 1 /* all gen. events in INEL > 0*/);

      auto groupedCollisions = getGroupedCollisions(collisions, mcCollision.globalIndex());
      // Check if there is at least one of the reconstructed collisions associated to this MC collision
      // If so, we consider it
      bool atLeastOne = false;
      int biggestNContribs = -1;
      int nCollisions = 0;
      float centralitydata = 100.5f;
      for (auto const& collision : groupedCollisions) {
        // event selections
        if (applySel8 && !collision.sel8())
          continue;

        if (applyVtxZ && std::fabs(collision.posZ()) > zVtx)
          continue;

        if (rejectITSROFBorder && !collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
          continue;
        }

        if (rejectTFBorder && !collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
          continue;
        }

        if (requireVertexITSTPC && !collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) {
          continue;
        }

        if (requireIsGoodZvtxFT0VsPV && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
          continue;
        }

        if (requireIsVertexTOFmatched && !collision.selection_bit(o2::aod::evsel::kIsVertexTOFmatched)) {
          continue;
        }

        if (requireIsVertexTRDmatched && !collision.selection_bit(o2::aod::evsel::kIsVertexTRDmatched)) {
          continue;
        }

        if (rejectSameBunchPileup && !collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
          continue;
        }

        if (requireInel0 && collision.multNTracksPVeta1() < 1) {
          continue;
        }

        if (biggestNContribs < collision.multPVTotalContributors()) {
          biggestNContribs = collision.multPVTotalContributors();
          if (centralityEstimator == Option::kFT0C)
            centralitydata = collision.centFT0C();
          if (centralityEstimator == Option::kFT0M)
            centralitydata = collision.centFT0M();
          if (centralityEstimator == Option::kFV0A)
            centralitydata = collision.centFV0A();
          if (centralityEstimator == Option::kNGlobal)
            centralitydata = collision.centNGlobal();
        }
        nCollisions++;

        atLeastOne = true;
      }

      registryMC.fill(HIST("hCentralityVsNcoll_beforeEvSel"), centralitydata, groupedCollisions.size());
      registryMC.fill(HIST("hCentralityVsNcoll_afterEvSel"), centralitydata, nCollisions);
      registryMC.fill(HIST("hCentralityVsMultMC"), centralitydata, mcCollision.multMCNParticlesEta05());

      registryQC.fill(HIST("hVertexZGen"), mcCollision.posZ());

      if (atLeastOne) {
        registryMC.fill(HIST("hGenEvents"), mcCollision.multMCNParticlesEta05(), 2 /* at least 1 rec. event*/);
        registryMC.fill(HIST("hGenEvents"), mcCollision.multMCNParticlesEta05(), 3 /* at least 1 rec. event in INEL > 0*/);
        registryMC.fill(HIST("hGenEventCentrality"), centralitydata);
      }
    }
    return;
  }

  void processData(SelCollisions::iterator const& collision,
                   CascadeCandidates const& fullCascades,
                   DaughterTracks const&)
  {
    // Fill event counter before event selection
    registryData.fill(HIST("number_of_events_data"), 0);

    // Initialize CCDB objects using the BC info
    initCCDB(collision);

    // Define the event centrality using different estimators
    float centrality = -1;
    float multiplicity = -1;

    if (centralityEstimator == Option::kFT0C) {
      centrality = collision.centFT0C();
      multiplicity = collision.multFT0C();
    }
    if (centralityEstimator == Option::kFT0M) {
      centrality = collision.centFT0M();
      multiplicity = collision.multFT0C() + collision.multFT0A();
    }
    if (centralityEstimator == Option::kFV0A) {
      centrality = collision.centFV0A();
      multiplicity = collision.multFV0A();
    }
    if (centralityEstimator == Option::kNGlobal) {
      centrality = collision.centNGlobal();
      multiplicity = collision.multNTracksGlobal();
    }

    registryData.fill(HIST("number_of_events_data_vs_centrality"), 0, centrality);

    // event selections
    if (applySel8 && !collision.sel8())
      return;
    registryData.fill(HIST("number_of_events_data"), 1);
    registryData.fill(HIST("number_of_events_data_vs_centrality"), 1, centrality);

    if (applyVtxZ && std::fabs(collision.posZ()) > zVtx)
      return;
    registryData.fill(HIST("number_of_events_data"), 2);
    registryData.fill(HIST("number_of_events_data_vs_centrality"), 2, centrality);

    if (rejectITSROFBorder && !collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
      return;
    }
    registryData.fill(HIST("number_of_events_data"), 3 /* Not at ITS ROF border */);
    registryData.fill(HIST("number_of_events_data_vs_centrality"), 3, centrality);

    if (rejectTFBorder && !collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
      return;
    }
    registryData.fill(HIST("number_of_events_data"), 4 /* Not at TF border */);
    registryData.fill(HIST("number_of_events_data_vs_centrality"), 4, centrality);

    if (requireVertexITSTPC && !collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) {
      return;
    }
    registryData.fill(HIST("number_of_events_data"), 5 /* Contains at least one ITS-TPC track */);
    registryData.fill(HIST("number_of_events_data_vs_centrality"), 5, centrality);

    if (requireIsGoodZvtxFT0VsPV && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return;
    }
    registryData.fill(HIST("number_of_events_data"), 6 /* PV position consistency check */);
    registryData.fill(HIST("number_of_events_data_vs_centrality"), 6, centrality);

    if (requireIsVertexTOFmatched && !collision.selection_bit(o2::aod::evsel::kIsVertexTOFmatched)) {
      return;
    }
    registryData.fill(HIST("number_of_events_data"), 7 /* PV with at least one contributor matched with TOF */);
    registryData.fill(HIST("number_of_events_data_vs_centrality"), 7, centrality);

    if (requireIsVertexTRDmatched && !collision.selection_bit(o2::aod::evsel::kIsVertexTRDmatched)) {
      return;
    }
    registryData.fill(HIST("number_of_events_data"), 8 /* PV with at least one contributor matched with TRD */);
    registryData.fill(HIST("number_of_events_data_vs_centrality"), 8, centrality);

    if (rejectSameBunchPileup && !collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      return;
    }
    registryData.fill(HIST("number_of_events_data"), 9 /* Not at same bunch pile-up */);
    registryData.fill(HIST("number_of_events_data_vs_centrality"), 9, centrality);

    if (requireInel0 && collision.multNTracksPVeta1() < 1) {
      return;
    }
    registryData.fill(HIST("number_of_events_data"), 10 /* INEL > 0 */);
    registryData.fill(HIST("number_of_events_data_vs_centrality"), 10, centrality);

    // Store the Zvtx
    registryQC.fill(HIST("hVertexZdata"), collision.posZ());

    // Store the event centrality
    registryData.fill(HIST("hCentEstimator"), centrality);
    registryData.fill(HIST("hCentralityVsNch"), centrality, collision.multNTracksPVeta1());
    registryData.fill(HIST("hCentralityVsMultiplicity"), centrality, multiplicity);

    // Loop over cascades
    for (const auto& casc : fullCascades) {
      if (etaMin > casc.bacheloreta() || casc.bacheloreta() > etaMax ||
          etaMin > casc.negativeeta() || casc.negativeeta() > etaMax ||
          etaMin > casc.positiveeta() || casc.positiveeta() > etaMax)
        continue; // remove acceptance that's badly reproduced by MC / superfluous in future

      // Get cascade daughters
      auto bach = casc.bachTrackExtra_as<DaughterTracks>();
      auto pos = casc.posTrackExtra_as<DaughterTracks>();
      auto neg = casc.negTrackExtra_as<DaughterTracks>();

      // ------------------------------------- Store selctions distribution for QC
      registryQC.fill(HIST("hv0cosPAdata"), casc.v0cosPA(collision.posX(), collision.posY(), collision.posZ()));
      registryQC.fill(HIST("hcasccosPAdata"), casc.casccosPA(collision.posX(), collision.posY(), collision.posZ()));
      registryQC.fill(HIST("hv0radiusdata"), casc.v0radius());
      registryQC.fill(HIST("hcascradiusdata"), casc.cascradius());
      registryQC.fill(HIST("hdcaV0daughtersdata"), casc.dcaV0daughters());
      registryQC.fill(HIST("hdcacascdaughtersdata"), casc.dcacascdaughters());
      registryQC.fill(HIST("hdcapostopvdata"), casc.dcapostopv());
      registryQC.fill(HIST("hdcanegtopvdata"), casc.dcanegtopv());
      registryQC.fill(HIST("hdcabachtopvdata"), casc.dcabachtopv());
      registryQC.fill(HIST("hdcav0topvdata"), casc.dcav0topv(collision.posX(), collision.posY(), collision.posZ()));

      // ------------------------------------- Store selctions distribution for analysis
      if (casc.sign() < 0) {
        registryData.fill(HIST("hMassXineg"), centrality, casc.pt(), casc.mXi());
        registryData.fill(HIST("hMassOmeganeg"), centrality, casc.pt(), casc.mOmega());
      }
      if (casc.sign() > 0) {
        registryData.fill(HIST("hMassXipos"), centrality, casc.pt(), casc.mXi());
        registryData.fill(HIST("hMassOmegapos"), centrality, casc.pt(), casc.mOmega());
      }

      if (casc.sign() < 0 && passedXiSelection(casc, pos, neg, bach, collision)) {
        registryData.fill(HIST("hMassXinegSelected"), centrality, casc.pt(), casc.mXi());
      }
      if (casc.sign() < 0 && passedOmegaSelection(casc, pos, neg, bach, collision)) {
        registryData.fill(HIST("hMassOmeganegSelected"), centrality, casc.pt(), casc.mOmega());
      }
      if (casc.sign() > 0 && passedXiSelection(casc, pos, neg, bach, collision)) {
        registryData.fill(HIST("hMassXiposSelected"), centrality, casc.pt(), casc.mXi());
      }
      if (casc.sign() > 0 && passedOmegaSelection(casc, pos, neg, bach, collision)) {
        registryData.fill(HIST("hMassOmegaposSelected"), centrality, casc.pt(), casc.mOmega());
      }
    }
  }

  PROCESS_SWITCH(CascadeAnalysisLightIonsDerivedData, processData, "Process data", true);

  void processMonteCarloRec(SimCollisions const& RecCols, CascadeMCCandidates const& fullCascades, DaughterTracks const&, CascadeMCCores const&)
  {
    for (const auto& RecCol : RecCols) {
      // Fill event counter before event selection
      registryMC.fill(HIST("number_of_events_mc_rec"), 0);

      // Initialize CCDB objects using the BC info
      initCCDB(RecCol);

      // Define the event centrality using different estimators
      float centralityMcRec = -1;
      float multiplicityMcRec = -1;

      if (centralityEstimator == Option::kFT0C) {
        centralityMcRec = RecCol.centFT0C();
        multiplicityMcRec = RecCol.multFT0C();
      }
      if (centralityEstimator == Option::kFT0M) {
        centralityMcRec = RecCol.centFT0M();
        multiplicityMcRec = RecCol.multFT0C() + RecCol.multFT0A();
      }
      if (centralityEstimator == Option::kFV0A) {
        centralityMcRec = RecCol.centFV0A();
        multiplicityMcRec = RecCol.multFV0A();
      }
      if (centralityEstimator == Option::kNGlobal) {
        centralityMcRec = RecCol.centNGlobal();
        multiplicityMcRec = RecCol.multNTracksGlobal();
      }

      registryMC.fill(HIST("number_of_events_mc_rec_vs_centrality"), 0, centralityMcRec);

      // event selections
      if (applySel8 && !RecCol.sel8())
        continue;
      registryMC.fill(HIST("number_of_events_mc_rec"), 1);
      registryMC.fill(HIST("number_of_events_mc_rec_vs_centrality"), 1, centralityMcRec);

      if (applyVtxZ && std::fabs(RecCol.posZ()) > zVtx)
        continue;
      registryMC.fill(HIST("number_of_events_mc_rec"), 2);
      registryMC.fill(HIST("number_of_events_mc_rec_vs_centrality"), 2, centralityMcRec);

      if (rejectITSROFBorder && !RecCol.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
        continue;
      }
      registryMC.fill(HIST("number_of_events_mc_rec"), 3 /* Not at ITS ROF border */);
      registryMC.fill(HIST("number_of_events_mc_rec_vs_centrality"), 3, centralityMcRec);

      if (rejectTFBorder && !RecCol.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
        continue;
      }
      registryMC.fill(HIST("number_of_events_mc_rec"), 4 /* Not at TF border */);
      registryMC.fill(HIST("number_of_events_mc_rec_vs_centrality"), 4, centralityMcRec);

      if (requireVertexITSTPC && !RecCol.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) {
        continue;
      }
      registryMC.fill(HIST("number_of_events_mc_rec"), 5 /* Contains at least one ITS-TPC track */);
      registryMC.fill(HIST("number_of_events_mc_rec_vs_centrality"), 5, centralityMcRec);

      if (requireIsGoodZvtxFT0VsPV && !RecCol.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
        continue;
      }
      registryMC.fill(HIST("number_of_events_mc_rec"), 6 /* PV position consistency check */);
      registryMC.fill(HIST("number_of_events_mc_rec_vs_centrality"), 6, centralityMcRec);

      if (requireIsVertexTOFmatched && !RecCol.selection_bit(o2::aod::evsel::kIsVertexTOFmatched)) {
        continue;
      }
      registryMC.fill(HIST("number_of_events_mc_rec"), 7 /* PV with at least one contributor matched with TOF */);
      registryMC.fill(HIST("number_of_events_mc_rec_vs_centrality"), 7, centralityMcRec);

      if (requireIsVertexTRDmatched && !RecCol.selection_bit(o2::aod::evsel::kIsVertexTRDmatched)) {
        continue;
      }
      registryMC.fill(HIST("number_of_events_mc_rec"), 8 /* PV with at least one contributor matched with TRD */);
      registryMC.fill(HIST("number_of_events_mc_rec_vs_centrality"), 8, centralityMcRec);

      if (rejectSameBunchPileup && !RecCol.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
        continue;
      }
      registryMC.fill(HIST("number_of_events_mc_rec"), 9 /* Not at same bunch pile-up */);
      registryMC.fill(HIST("number_of_events_mc_rec_vs_centrality"), 9, centralityMcRec);

      if (requireInel0 && collision.multNTracksPVeta1() < 1) {
        continue;
      }
      registryMC.fill(HIST("number_of_events_mc_rec"), 10 /* INEL > 0 */);
      registryMC.fill(HIST("number_of_events_mc_rec_vs_centrality"), 10, centralityMcRec);

      // Store the Zvtx
      registryQC.fill(HIST("hVertexZRec"), RecCol.posZ());

      // Store the event centrality
      registryMC.fill(HIST("hCentEstimator_truerec"), centralityMcRec);
      registryMC.fill(HIST("hCentralityVsNch_truerec"), centralityMcRec, RecCol.multNTracksPVeta1());
      registryMC.fill(HIST("hCentralityVsMultiplicity_truerec"), centralityMcRec, multiplicityMcRec);

      for (const auto& casc : fullCascades) {
        if (etaMin > casc.bacheloreta() || casc.bacheloreta() > etaMax ||
            etaMin > casc.negativeeta() || casc.negativeeta() > etaMax ||
            etaMin > casc.positiveeta() || casc.positiveeta() > etaMax)
          continue; // remove acceptance that's badly reproduced by MC / superfluous in future

        if (!casc.has_cascMCCore())
          continue;

        auto cascMC = casc.template cascMCCore_as<CascadeMCCores>();

        auto bach = casc.bachTrackExtra_as<DaughterTracks>();
        auto pos = casc.posTrackExtra_as<DaughterTracks>();
        auto neg = casc.negTrackExtra_as<DaughterTracks>();

        int pdgParent = cascMC.pdgCode();
        bool isPhysPrim = cascMC.isPhysicalPrimary();
        if (pdgParent == 0)
          continue;
        if (!isPhysPrim)
          continue;

        float ptmc = RecoDecay::sqrtSumOfSquares(cascMC.pxMC(), cascMC.pyMC());

        // ------------------------------------- Store selctions distribution for QC
        registryQC.fill(HIST("hv0cosPARec"), casc.v0cosPA(RecCol.posX(), RecCol.posY(), RecCol.posZ()));
        registryQC.fill(HIST("hcasccosPARec"), casc.casccosPA(RecCol.posX(), RecCol.posY(), RecCol.posZ()));
        registryQC.fill(HIST("hv0radiusRec"), casc.v0radius());
        registryQC.fill(HIST("hcascradiusRec"), casc.cascradius());
        registryQC.fill(HIST("hdcaV0daughtersRec"), casc.dcaV0daughters());
        registryQC.fill(HIST("hdcacascdaughtersRec"), casc.dcacascdaughters());
        registryQC.fill(HIST("hdcapostopvRec"), casc.dcapostopv());
        registryQC.fill(HIST("hdcanegtopvRec"), casc.dcanegtopv());
        registryQC.fill(HIST("hdcabachtopvRec"), casc.dcabachtopv());
        registryQC.fill(HIST("hdcav0topvRec"), casc.dcav0topv(RecCol.posX(), RecCol.posY(), RecCol.posZ()));

        // ------------------------------------- Store selctions distribution for analysis
        if (casc.sign() < 0) {
          if (pdgParent == kXiMinus) {
            registryMC.fill(HIST("hMassXineg_truerec"), centralityMcRec, ptmc, casc.mXi());
          }
          if (pdgParent == kOmegaMinus) {
            registryMC.fill(HIST("hMassOmeganeg_truerec"), centralityMcRec, ptmc, casc.mOmega());
          }
        }

        if (casc.sign() > 0) {
          if (pdgParent == kXiPlusBar) {
            registryMC.fill(HIST("hMassXipos_truerec"), centralityMcRec, ptmc, casc.mXi());
          }
          if (pdgParent == kOmegaPlusBar) {
            registryMC.fill(HIST("hMassOmegapos_truerec"), centralityMcRec, ptmc, casc.mOmega());
          }
        }

        if (casc.sign() < 0 && pdgParent == kXiMinus && passedXiSelection(casc, pos, neg, bach, RecCol)) {
          registryMC.fill(HIST("hMassXinegSelected_truerec"), centralityMcRec, ptmc, casc.mXi());
        }
        if (casc.sign() < 0 && pdgParent == kOmegaMinus && passedOmegaSelection(casc, pos, neg, bach, RecCol)) {
          registryMC.fill(HIST("hMassOmeganegSelected_truerec"), centralityMcRec, ptmc, casc.mOmega());
        }
        if (casc.sign() > 0 && pdgParent == kXiPlusBar && passedXiSelection(casc, pos, neg, bach, RecCol)) {
          registryMC.fill(HIST("hMassXiposSelected_truerec"), centralityMcRec, ptmc, casc.mXi());
        }
        if (casc.sign() > 0 && pdgParent == kOmegaPlusBar && passedOmegaSelection(casc, pos, neg, bach, RecCol)) {
          registryMC.fill(HIST("hMassOmegaposSelected_truerec"), centralityMcRec, ptmc, casc.mOmega());
        }
      } // casc loop
    } // rec.collision loop
  }

  PROCESS_SWITCH(CascadeAnalysisLightIonsDerivedData, processMonteCarloRec, "Process MC Rec", false);

  void processMonteCarloGen(CollisionMCTrueTable const& mcCollisions,
                            CascadeMCCores const& CascMCCores,
                            SimCollisions const& RecCols)
  {
    // Fill generated event information (for event loss/splitting estimation)
    fillGeneratedEventProperties(mcCollisions, RecCols);
    std::vector<int> listBestCollisionIdx = getListOfRecoCollIndices(mcCollisions, RecCols);
    for (auto const& cascMC : CascMCCores) {
      int pdgParent = cascMC.pdgCode();
      bool isPhysPrim = cascMC.isPhysicalPrimary();
      if (pdgParent == 0)
        continue;
      if (!isPhysPrim)
        continue;

      float ptmc = cascMC.ptMC();

      auto mcCollision = cascMC.template straMCCollision_as<CollisionMCTrueTable>();

      // event selections
      if (applyZVtxSelOnMCPV && std::abs(mcCollision.posZ()) > zVtx)
        continue;
      // Store the Zvtx
      registryQC.fill(HIST("hVertexZGen"), mcCollision.posZ());

      if (requireInel0OnMC && mcCollision.multMCNParticlesEta10() < 1) {
        continue;
      }

      float centralityMC = 100.5f;

      if (listBestCollisionIdx[mcCollision.globalIndex()] > -1) {
        auto collision = RecCols.iteratorAt(listBestCollisionIdx[mcCollision.globalIndex()]);
        if (centralityEstimator == Option::kFT0C)
          centralityMC = collision.centFT0C();
        if (centralityEstimator == Option::kFT0M)
          centralityMC = collision.centFT0M();
        if (centralityEstimator == Option::kFV0A)
          centralityMC = collision.centFV0A();
        if (centralityEstimator == Option::kNGlobal)
          centralityMC = collision.centNGlobal();

        if (cascMC.pdgCode() == kXiMinus && std::abs(cascMC.rapidityMC(0)) < rapcut) {
          registryMC.fill(HIST("h2dGenXiMinusVsMultMC_RecoedEvt"), mcCollision.multMCNParticlesEta05(), ptmc);
        }
        if (cascMC.pdgCode() == kXiPlusBar && std::abs(cascMC.rapidityMC(0)) < rapcut) {
          registryMC.fill(HIST("h2dGenXiPlusVsMultMC_RecoedEvt"), mcCollision.multMCNParticlesEta05(), ptmc);
        }
        if (cascMC.pdgCode() == kOmegaMinus && std::abs(cascMC.rapidityMC(2)) < rapcut) {
          registryMC.fill(HIST("h2dGenOmegaMinusVsMultMC_RecoedEvt"), mcCollision.multMCNParticlesEta05(), ptmc);
        }
        if (cascMC.pdgCode() == kOmegaPlusBar && std::abs(cascMC.rapidityMC(2)) < rapcut) {
          registryMC.fill(HIST("h2dGenOmegaPlusVsMultMC_RecoedEvt"), mcCollision.multMCNParticlesEta05(), ptmc);
        }
      }

      if (cascMC.pdgCode() == kXiMinus && std::abs(cascMC.rapidityMC(0)) < rapcut) {
        registryMC.fill(HIST("h2dGenXiMinus"), centralityMC, ptmc);
        registryMC.fill(HIST("h2dGenXiMinusVsMultMC"), mcCollision.multMCNParticlesEta05(), ptmc);
      }
      if (cascMC.pdgCode() == kXiPlusBar && std::abs(cascMC.rapidityMC(0)) < rapcut) {
        registryMC.fill(HIST("h2dGenXiPlus"), centralityMC, ptmc);
        registryMC.fill(HIST("h2dGenXiPlusVsMultMC"), mcCollision.multMCNParticlesEta05(), ptmc);
      }
      if (cascMC.pdgCode() == kOmegaMinus && std::abs(cascMC.rapidityMC(2)) < rapcut) {
        registryMC.fill(HIST("h2dGenOmegaMinus"), centralityMC, ptmc);
        registryMC.fill(HIST("h2dGenOmegaMinusVsMultMC"), mcCollision.multMCNParticlesEta05(), ptmc);
      }
      if (cascMC.pdgCode() == kOmegaPlusBar && std::abs(cascMC.rapidityMC(2)) < rapcut) {
        registryMC.fill(HIST("h2dGenOmegaPlus"), centralityMC, ptmc);
        registryMC.fill(HIST("h2dGenOmegaPlusVsMultMC"), mcCollision.multMCNParticlesEta05(), ptmc);
      }
    } // cascMC loop
  }

  PROCESS_SWITCH(CascadeAnalysisLightIonsDerivedData, processMonteCarloGen, "Process MC Gen", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<CascadeAnalysisLightIonsDerivedData>(cfgc)};
}
