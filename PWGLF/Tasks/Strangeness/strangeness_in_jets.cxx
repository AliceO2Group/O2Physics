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
/// \author Alberto Caliva (alberto.caliva@cern.ch), Francesca Ercolessi (francesca.ercolessi@cern.ch)
/// \since May 22, 2024

#include <TDatabasePDG.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TPDGCode.h>
#include <TVector2.h>
#include <TVector3.h>
#include <cmath>
#include <vector>

#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "ReconstructionDataFormats/Track.h"

using namespace std;
using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;
using std::array;

using SelectedCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms>;
using SimCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::McCollisionLabels>;

using FullTracks = soa::Join<aod::Tracks, aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::TracksDCA, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr>;

using MCTracks = soa::Join<aod::Tracks, aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::TracksDCA, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::McTrackLabels>;

struct strangeness_in_jets {

  // QC Histograms
  HistogramRegistry registryQC{
    "registryQC",
    {},
    OutputObjHandlingPolicy::AnalysisObject,
    true,
    true};

  // Analysis Histograms: Data
  HistogramRegistry registryData{
    "registryData",
    {},
    OutputObjHandlingPolicy::AnalysisObject,
    true,
    true};

  // MC Histograms
  HistogramRegistry registryMC{
    "registryMC",
    {},
    OutputObjHandlingPolicy::AnalysisObject,
    true,
    true};

  // Global Parameters
  Configurable<float> ptLeadingMin{"ptLeadingMin", 5.0f, "pt leading min"};
  Configurable<float> Rjet{"Rjet", 0.3f, "jet resolution parameter R"};
  Configurable<float> Rmax{"Rmax", 0.3f, "radius of the jet and UE cones"};
  Configurable<float> zVtx{"zVtx", 10.0f, "z vertex cut"};

  // Track Parameters
  Configurable<float> minITSnCls{"minITSnCls", 4.0f, "min number of ITS clusters"};
  Configurable<float> minTPCnClsFound{"minTPCnClsFound", 80.0f, "min number of found TPC clusters"};
  Configurable<float> minNCrossedRowsTPC{"minNCrossedRowsTPC", 80.0f, "min number of TPC crossed rows"};
  Configurable<float> maxChi2TPC{"maxChi2TPC", 4.0f, "max chi2 per cluster TPC"};
  Configurable<float> maxChi2ITS{"maxChi2ITS", 36.0f, "max chi2 per cluster ITS"};
  Configurable<float> etaMin{"etaMin", -0.8f, "eta min"};
  Configurable<float> etaMax{"etaMax", +0.8f, "eta max"};
  Configurable<float> ptMin_V0_proton{"ptMin_V0_proton", 0.3f, "pt min of proton from V0"};
  Configurable<float> ptMax_V0_proton{"ptMax_V0_proton", 10.0f, "pt max of proton from V0"};
  Configurable<float> ptMin_V0_pion{"ptMin_V0_pion", 0.1f, "pt min of pion from V0"};
  Configurable<float> ptMax_V0_pion{"ptMax_V0_pion", 1.5f, "pt max of pion from V0"};
  Configurable<float> ptMin_K0_pion{"ptMin_K0_pion", 0.3f, "pt min of pion from K0"};
  Configurable<float> ptMax_K0_pion{"ptMax_K0_pion", 10.0f, "pt max of pion from K0"};
  Configurable<float> nsigmaTPCmin{"nsigmaTPCmin", -3.0f, "Minimum nsigma TPC"};
  Configurable<float> nsigmaTPCmax{"nsigmaTPCmax", +3.0f, "Maximum nsigma TPC"};
  Configurable<float> nsigmaTOFmin{"nsigmaTOFmin", -3.0f, "Minimum nsigma TOF"};
  Configurable<float> nsigmaTOFmax{"nsigmaTOFmax", +3.0f, "Maximum nsigma TOF"};
  Configurable<float> dcanegtoPVmin{"dcanegtoPVmin", 0.1f, "Minimum DCA Neg To PV"};
  Configurable<float> dcapostoPVmin{"dcapostoPVmin", 0.1f, "Minimum DCA Pos To PV"};
  Configurable<bool> requireTOF{"requireTOF", false, "require TOF hit"};
  Configurable<bool> requireITS{"requireITS", false, "require ITS hit"};
  Configurable<float> dcaxyMax{"dcaxyMax", 0.1f, "Maximum DCAxy to primary vertex"};
  Configurable<float> dcazMax{"dcazMax", 0.1f, "Maximum DCAz to primary vertex"};

  // V0 Parameters
  Configurable<float> yMin{"yMin", -0.5f, "minimum y"};
  Configurable<float> yMax{"yMax", +0.5f, "maximum y"};
  Configurable<float> minimumV0Radius{"minimumV0Radius", 0.5f, "Minimum V0 Radius"};
  Configurable<float> maximumV0Radius{"maximumV0Radius", 40.0f, "Maximum V0 Radius"};
  Configurable<float> v0cospaMin{"v0cospaMin", 0.99f, "Minimum V0 CosPA"};
  Configurable<float> dcaV0DaughtersMax{"dcaV0DaughtersMax", 0.5f, "Maximum DCA Daughters"};

  // Cascade Parameters
  Configurable<float> minimumCascRadius{"minimumCascRadius", 0.1f, "Minimum Cascade Radius"};
  Configurable<float> maximumCascRadius{"maximumCascRadius", 40.0f, "Maximum Cascade Radius"};
  Configurable<float> casccospaMin{"casccospaMin", 0.99f, "Minimum Cascade CosPA"};
  Configurable<float> dcabachtopvMin{"dcabachtopvMin", 0.1f, "Minimum DCA bachelor to PV"};
  Configurable<float> dcaV0topvMin{"dcaV0topvMin", 0.1f, "Minimum DCA V0 to PV"};
  Configurable<float> dcaCascDaughtersMax{"dcaCascDaughtersMax", 0.5f, "Maximum DCA Daughters"};

  // Debug
  Configurable<bool> debug_level0{"debug_level0", false, "debug0"};
  Configurable<bool> debug_level1{"debug_level1", false, "debug1"};
  Configurable<bool> debug_level2{"debug_level2", false, "debug2"};
  Configurable<bool> debug_level3{"debug_level3", false, "debug3"};
  Configurable<bool> debug_level4{"debug_level4", false, "debug4"};
  Configurable<bool> debug_level5{"debug_level5", false, "debug5"};
  Configurable<bool> debug_level6{"debug_level6", false, "debug6"};
  Configurable<bool> debug_level7{"debug_level7", false, "debug7"};
  Configurable<bool> debug_level8{"debug_level8", false, "debug8"};
  Configurable<bool> debug_level9{"debug_level9", false, "debug9"};
  Configurable<bool> debug_level10{"debug_level10", false, "debug10"};
  Configurable<bool> debug_level11{"debug_level11", false, "debug11"};
  Configurable<bool> debug_level12{"debug_level12", false, "debug12"};
  Configurable<bool> debug_level13{"debug_level13", false, "debug13"};
  Configurable<bool> debug_level14{"debug_level14", false, "debug14"};
  Configurable<bool> debug_level15{"debug_level15", false, "debug15"};

  void init(InitContext const&)
  {
    // Global Properties and QC
    registryQC.add("number_of_events_data", "number of events in data", HistType::kTH1F, {{15, 0, 15, "Event Cuts"}});
    registryQC.add("number_of_events_mc", "number of events in mc", HistType::kTH1F, {{15, 0, 15, "Event Cuts"}});

    // Multiplicity Binning
    std::vector<double> multBinning = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100};
    AxisSpec multAxis = {multBinning, "FT0C percentile"};

    // Histograms (Pions)
    registryData.add("piplus_tpc_in_jet", "piplus_tpc_in_jet", HistType::kTH3F, {multBinning, {100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -10, 10, "n#sigma_{TPC}"}});
    registryData.add("piplus_tof_in_jet", "piplus_tof_in_jet", HistType::kTH3F, {multBinning, {100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -10, 10, "n#sigma_{TOF}"}});
    registryData.add("piplus_tpc_in_ue", "piplus_tpc_in_ue", HistType::kTH3F, {multBinning, {100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -10, 10, "n#sigma_{TPC}"}});
    registryData.add("piplus_tof_in_ue", "piplus_tof_in_ue", HistType::kTH3F, {multBinning, {100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -10, 10, "n#sigma_{TOF}"}});
    registryData.add("piminus_tpc_in_jet", "piminus_tpc_in_jet", HistType::kTH3F, {multBinning, {100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -10, 10, "n#sigma_{TPC}"}});
    registryData.add("piminus_tof_in_jet", "piminus_tof_in_jet", HistType::kTH3F, {multBinning, {100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -10, 10, "n#sigma_{TOF}"}});
    registryData.add("piminus_tpc_in_ue", "piminus_tpc_in_ue", HistType::kTH3F, {multBinning, {100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -10, 10, "n#sigma_{TPC}"}});
    registryData.add("piminus_tof_in_ue", "piminus_tof_in_ue", HistType::kTH3F, {multBinning, {100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}, {200, -10, 10, "n#sigma_{TOF}"}});

    // Histograms (Lambda)
    registryData.add("Lambda_in_jet", "Lambda_in_jet", HistType::kTH3F, {multBinning, {100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}, {200, 1.09, 1.14, "m_{p#pi} (GeV/#it{c}^{2})"}});
    registryData.add("AntiLambda_in_jet", "AntiLambda_in_jet", HistType::kTH3F, {multBinning, {100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}, {200, 1.09, 1.14, "m_{p#pi} (GeV/#it{c}^{2})"}});
    registryData.add("Lambda_in_ue", "Lambda_in_ue", HistType::kTH3F, {multBinning, {100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}, {200, 1.09, 1.14, "m_{p#pi} (GeV/#it{c}^{2})"}});
    registryData.add("AntiLambda_in_ue", "AntiLambda_in_ue", HistType::kTH3F, {multBinning, {100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}, {200, 1.09, 1.14, "m_{p#pi} (GeV/#it{c}^{2})"}});

    // Histograms (K0s)
    registryData.add("K0s_in_jet", "K0s_in_jet", HistType::kTH3F, {multBinning, {100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}, {200, 0.44, 0.56, "m_{#pi#pi} (GeV/#it{c}^{2})"}});
    registryData.add("K0s_in_ue", "K0s_in_ue", HistType::kTH3F, {multBinning, {100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}, {200, 0.44, 0.56, "m_{#pi#pi} (GeV/#it{c}^{2})"}});

    // Histograms (Xi)
    registryData.add("XiPos_in_jet", "XiPos_in_jet", HistType::kTH3F, {multBinning, {100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}, {200, 1.28, 1.36, "m_{p#pi#pi} (GeV/#it{c}^{2})"}});
    registryData.add("XiPos_in_ue", "XiPos_in_ue", HistType::kTH3F, {multBinning, {100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}, {200, 1.28, 1.36, "m_{p#pi#pi} (GeV/#it{c}^{2})"}});
    registryData.add("XiNeg_in_jet", "XiNeg_in_jet", HistType::kTH3F, {multBinning, {100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}, {200, 1.28, 1.36, "m_{p#pi#pi} (GeV/#it{c}^{2})"}});
    registryData.add("XiNeg_in_ue", "XiNeg_in_ue", HistType::kTH3F, {multBinning, {100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}, {200, 1.28, 1.36, "m_{p#pi#pi} (GeV/#it{c}^{2})"}});

    // Histograms (Omega)
    registryData.add("OmegaPos_in_jet", "OmegaPos_in_jet", HistType::kTH3F, {multBinning, {100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}, {200, 1.63, 1.71, "m_{p#piK} (GeV/#it{c}^{2})"}});
    registryData.add("OmegaPos_in_ue", "OmegaPos_in_ue", HistType::kTH3F, {multBinning, {100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}, {200, 1.63, 1.71, "m_{p#piK} (GeV/#it{c}^{2})"}});
    registryData.add("OmegaNeg_in_jet", "OmegaNeg_in_jet", HistType::kTH3F, {multBinning, {100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}, {200, 1.63, 1.71, "m_{p#piK} (GeV/#it{c}^{2})"}});
    registryData.add("OmegaNeg_in_ue", "OmegaNeg_in_ue", HistType::kTH3F, {multBinning, {100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}, {200, 1.63, 1.71, "m_{p#piK} (GeV/#it{c}^{2})"}});

    // Histograms (MC)
    registryMC.add("K0s_generated", "K0s_generated", HistType::kTH2F, {multBinning, {100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryMC.add("Lambda_generated", "Lambda_generated", HistType::kTH2F, {multBinning, {100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryMC.add("AntiLambda_generated", "AntiLambda_generated", HistType::kTH2F, {multBinning, {100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryMC.add("XiPos_generated", "XiPos_generated", HistType::kTH2F, {multBinning, {100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryMC.add("XiNeg_generated", "XiNeg_generated", HistType::kTH2F, {multBinning, {100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryMC.add("OmegaPos_generated", "OmegaPos_generated", HistType::kTH2F, {multBinning, {100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryMC.add("OmegaNeg_generated", "OmegaNeg_generated", HistType::kTH2F, {multBinning, {100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryMC.add("K0s_reconstructed", "K0s_reconstructed", HistType::kTH2F, {multBinning, {100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryMC.add("Lambda_reconstructed", "Lambda_reconstructed", HistType::kTH2F, {multBinning, {100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryMC.add("AntiLambda_reconstructed", "AntiLambda_reconstructed", HistType::kTH2F, {multBinning, {100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryMC.add("K0s_reconstructed_incl", "K0s_reconstructed_incl", HistType::kTH2F, {multBinning, {100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryMC.add("Lambda_reconstructed_incl", "Lambda_reconstructed_incl", HistType::kTH2F, {multBinning, {100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryMC.add("AntiLambda_reconstructed_incl", "AntiLambda_reconstructed_incl", HistType::kTH2F, {multBinning, {100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryMC.add("XiPos_reconstructed", "XiPos_reconstructed", HistType::kTH2F, {multBinning, {100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryMC.add("XiNeg_reconstructed", "XiNeg_reconstructed", HistType::kTH2F, {multBinning, {100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryMC.add("OmegaPos_reconstructed", "OmegaPos_reconstructed", HistType::kTH2F, {multBinning, {100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryMC.add("OmegaNeg_reconstructed", "OmegaNeg_reconstructed", HistType::kTH2F, {multBinning, {100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}});

    // Histograms for reweighting
    registryMC.add("K0s_eta_pt_jet", "K0s_eta_pt_jet", HistType::kTH2F, {{100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}, {18, -0.9, 0.9, "#eta"}});
    registryMC.add("K0s_eta_pt_ue", "K0s_eta_pt_ue", HistType::kTH2F, {{100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}, {18, -0.9, 0.9, "#eta"}});
    registryMC.add("Lambda_eta_pt_jet", "Lambda_eta_pt_jet", HistType::kTH2F, {{100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}, {18, -0.9, 0.9, "#eta"}});
    registryMC.add("Lambda_eta_pt_ue", "Lambda_eta_pt_ue", HistType::kTH2F, {{100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}, {18, -0.9, 0.9, "#eta"}});
    registryMC.add("AntiLambda_eta_pt_jet", "AntiLambda_eta_pt_jet", HistType::kTH2F, {{100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}, {18, -0.9, 0.9, "#eta"}});
    registryMC.add("AntiLambda_eta_pt_ue", "AntiLambda_eta_pt_ue", HistType::kTH2F, {{100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}, {18, -0.9, 0.9, "#eta"}});
    registryMC.add("K0s_eta_pt_pythia", "K0s_eta_pt_pythia", HistType::kTH2F, {{100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}, {18, -0.9, 0.9, "#eta"}});
    registryMC.add("Lambda_eta_pt_pythia", "Lambda_eta_pt_pythia", HistType::kTH2F, {{100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}, {18, -0.9, 0.9, "#eta"}});
    registryMC.add("AntiLambda_eta_pt_pythia", "AntiLambda_eta_pt_pythia", HistType::kTH2F, {{100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}, {18, -0.9, 0.9, "#eta"}});
  }

  template <typename T1>
  bool passedTrackSelectionForJets(const T1& track)
  {
    if (!track.hasITS())
      return false;
    if (track.itsNCls() < 2)
      return false;
    if (!track.hasTPC())
      return false;
    if (track.tpcNClsFound() < 70)
      return false;
    if (track.tpcNClsCrossedRows() < 70)
      return false;
    if (track.tpcChi2NCl() > 4)
      return false;
    if (track.itsChi2NCl() > 36)
      return false;
    if (track.eta() < -0.8 || track.eta() > 0.8)
      return false;
    if (track.pt() < 0.15)
      return false;
    if (TMath::Abs(track.dcaXY()) > 0.5)
      return false;
    if (TMath::Abs(track.dcaZ()) > 0.5)
      return false;
    return true;
  }

  template <typename pionTrack>
  bool passedTrackSelectionForPions(const pionTrack& track)
  {
    if (!track.hasITS())
      return false;
    if (track.itsNCls() < minITSnCls)
      return false;
    if (!track.hasTPC())
      return false;
    if (track.tpcNClsFound() < minTPCnClsFound)
      return false;
    if (track.tpcNClsCrossedRows() < minNCrossedRowsTPC)
      return false;
    if (track.tpcChi2NCl() > maxChi2TPC)
      return false;
    if (track.itsChi2NCl() > maxChi2ITS)
      return false;
    if (track.eta() < etaMin || track.eta() > etaMax)
      return false;
    if (TMath::Abs(track.dcaXY()) > dcaxyMax)
      return false;
    if (TMath::Abs(track.dcaZ()) > dcazMax)
      return false;

    // Rapidity Selection
    TLorentzVector lorentzVect;
    lorentzVect.SetXYZM(track.px(), track.py(), track.pz(), 0.13957021);
    if (lorentzVect.Rapidity() < yMin || lorentzVect.Rapidity() > yMax)
      return false;

    return true;
  }

  // Lambda Selections
  template <typename V, typename T1, typename T2>
  bool passedLambdaSelection(const V& v0, const T1& ptrack, const T2& ntrack)
  {
    // Single-Track Selections
    if (!passedSingleTrackSelection(ptrack))
      return false;
    if (!passedSingleTrackSelection(ntrack))
      return false;

    // Momentum Lambda Daughters
    TVector3 proton(v0.pxpos(), v0.pypos(), v0.pzpos());
    TVector3 pion(v0.pxneg(), v0.pyneg(), v0.pzneg());

    if (proton.Pt() < ptMin_V0_proton)
      return false;
    if (proton.Pt() > ptMax_V0_proton)
      return false;
    if (pion.Pt() < ptMin_V0_pion)
      return false;
    if (pion.Pt() > ptMax_V0_pion)
      return false;

    // V0 Selections
    if (v0.v0cosPA() < v0cospaMin)
      return false;
    if (v0.v0radius() < minimumV0Radius || v0.v0radius() > maximumV0Radius)
      return false;
    if (v0.dcaV0daughters() > dcaV0DaughtersMax)
      return false;
    if (v0.dcapostopv() < dcapostoPVmin)
      return false;
    if (v0.dcanegtopv() < dcanegtoPVmin)
      return false;

    // PID Selections (TPC)
    if (ptrack.tpcNSigmaPr() < nsigmaTPCmin || ptrack.tpcNSigmaPr() > nsigmaTPCmax)
      return false;
    if (ntrack.tpcNSigmaPi() < nsigmaTPCmin || ntrack.tpcNSigmaPi() > nsigmaTPCmax)
      return false;

    // PID Selections (TOF)
    if (requireTOF) {
      if (ptrack.tofNSigmaPr() < nsigmaTOFmin || ptrack.tofNSigmaPr() > nsigmaTOFmax)
        return false;
      if (ntrack.tofNSigmaPi() < nsigmaTOFmin || ntrack.tofNSigmaPi() > nsigmaTOFmax)
        return false;
    }

    // Rapidity Selection
    TLorentzVector lorentzVect;
    lorentzVect.SetXYZM(v0.px(), v0.py(), v0.pz(), 1.115683);
    if (lorentzVect.Rapidity() < yMin || lorentzVect.Rapidity() > yMax)
      return false;
    return true;
  }

  // AntiLambda Selections
  template <typename V, typename T1, typename T2>
  bool passedAntiLambdaSelection(const V& v0, const T1& ptrack, const T2& ntrack)
  {
    // Single-Track Selections
    if (!passedSingleTrackSelection(ptrack))
      return false;
    if (!passedSingleTrackSelection(ntrack))
      return false;

    // Momentum AntiLambda Daughters
    TVector3 pion(v0.pxpos(), v0.pypos(), v0.pzpos());
    TVector3 proton(v0.pxneg(), v0.pyneg(), v0.pzneg());

    if (proton.Pt() < ptMin_V0_proton)
      return false;
    if (proton.Pt() > ptMax_V0_proton)
      return false;
    if (pion.Pt() < ptMin_V0_pion)
      return false;
    if (pion.Pt() > ptMax_V0_pion)
      return false;

    // V0 Selections
    if (v0.v0cosPA() < v0cospaMin)
      return false;
    if (v0.v0radius() < minimumV0Radius || v0.v0radius() > maximumV0Radius)
      return false;
    if (v0.dcaV0daughters() > dcaV0DaughtersMax)
      return false;
    if (v0.dcapostopv() < dcapostoPVmin)
      return false;
    if (v0.dcanegtopv() < dcanegtoPVmin)
      return false;

    // PID Selections (TPC)
    if (ptrack.tpcNSigmaPi() < nsigmaTPCmin || ptrack.tpcNSigmaPi() > nsigmaTPCmax)
      return false;
    if (ntrack.tpcNSigmaPr() < nsigmaTPCmin || ntrack.tpcNSigmaPr() > nsigmaTPCmax)
      return false;

    // PID Selections (TOF)
    if (requireTOF) {
      if (ptrack.tofNSigmaPi() < nsigmaTOFmin || ptrack.tofNSigmaPi() > nsigmaTOFmax)
        return false;
      if (ntrack.tofNSigmaPr() < nsigmaTOFmin || ntrack.tofNSigmaPr() > nsigmaTOFmax)
        return false;
    }

    // Rapidity Selection
    TLorentzVector lorentzVect;
    lorentzVect.SetXYZM(v0.px(), v0.py(), v0.pz(), 1.115683);
    if (lorentzVect.Rapidity() < yMin || lorentzVect.Rapidity() > yMax)
      return false;
    return true;
  }

  // K0s Selections
  template <typename V, typename T1, typename T2>
  bool passedK0ShortSelection(const V& v0, const T1& ptrack, const T2& ntrack)
  {
    // Single-Track Selections
    if (!passedSingleTrackSelection(ptrack))
      return false;
    if (!passedSingleTrackSelection(ntrack))
      return false;

    // Momentum K0s Daughters
    TVector3 pion_pos(v0.pxpos(), v0.pypos(), v0.pzpos());
    TVector3 pion_neg(v0.pxneg(), v0.pyneg(), v0.pzneg());

    if (pion_pos.Pt() < ptMin_K0_pion)
      return false;
    if (pion_pos.Pt() > ptMax_K0_pion)
      return false;
    if (pion_neg.Pt() < ptMin_K0_pion)
      return false;
    if (pion_neg.Pt() > ptMax_K0_pion)
      return false;

    // V0 Selections
    if (v0.v0cosPA() < v0cospaMin)
      return false;
    if (v0.v0radius() < minimumV0Radius || v0.v0radius() > maximumV0Radius)
      return false;
    if (v0.dcaV0daughters() > dcaV0DaughtersMax)
      return false;
    if (v0.dcapostopv() < dcapostoPVmin)
      return false;
    if (v0.dcanegtopv() < dcanegtoPVmin)
      return false;

    // PID Selections (TPC)
    if (ptrack.tpcNSigmaPi() < nsigmaTPCmin || ptrack.tpcNSigmaPi() > nsigmaTPCmax)
      return false;
    if (ntrack.tpcNSigmaPi() < nsigmaTPCmin || ntrack.tpcNSigmaPi() > nsigmaTPCmax)
      return false;

    // PID Selections (TOF)
    if (requireTOF) {
      if (ptrack.tofNSigmaPi() < nsigmaTOFmin || ptrack.tofNSigmaPi() > nsigmaTOFmax)
        return false;
      if (ntrack.tofNSigmaPi() < nsigmaTOFmin || ntrack.tofNSigmaPi() > nsigmaTOFmax)
        return false;
    }

    // Rapidity Selection
    TLorentzVector lorentzVect;
    lorentzVect.SetXYZM(v0.px(), v0.py(), v0.pz(), 0.497614);
    if (lorentzVect.Rapidity() < yMin || lorentzVect.Rapidity() > yMax)
      return false;
    return true;
  }

  // Xi Selections
  template <typename CascType, typename T1, typename T2, typename T3, typename C>
  bool passedXiSelection(const CascType& casc, const T1& ptrack, const T2& ntrack, const T3& btrack, const C& coll)
  {
    if (!passedSingleTrackSelection(ptrack))
      return false;
    if (!passedSingleTrackSelection(ntrack))
      return false;
    if (!passedSingleTrackSelection(btrack))
      return false;

    // Xi+ Selection (Xi+ -> antiL + pi+)
    if (btrack.sign() > 0) {
      if (ntrack.pt() < ptMin_V0_proton)
        return false;
      if (ntrack.pt() > ptMax_V0_proton)
        return false;
      if (ptrack.pt() < ptMin_V0_pion)
        return false;
      if (ptrack.pt() > ptMax_V0_pion)
        return false;

      // PID Selections (TPC)
      if (ntrack.tpcNSigmaPr() < nsigmaTPCmin || ntrack.tpcNSigmaPr() > nsigmaTPCmax)
        return false;
      if (ptrack.tpcNSigmaPi() < nsigmaTPCmin || ptrack.tpcNSigmaPi() > nsigmaTPCmax)
        return false;

      // PID Selections (TOF)
      if (requireTOF) {
        if (ntrack.tofNSigmaPr() < nsigmaTOFmin || ntrack.tofNSigmaPr() > nsigmaTOFmax)
          return false;
        if (ptrack.tofNSigmaPi() < nsigmaTOFmin || ptrack.tofNSigmaPi() > nsigmaTOFmax)
          return false;
      }
    }

    // Xi- Selection (Xi- -> L + pi-)
    if (btrack.sign() < 0) {
      if (ptrack.pt() < ptMin_V0_proton)
        return false;
      if (ptrack.pt() > ptMax_V0_proton)
        return false;
      if (ntrack.pt() < ptMin_V0_pion)
        return false;
      if (ntrack.pt() > ptMax_V0_pion)
        return false;

      // PID Selections (TPC)
      if (ptrack.tpcNSigmaPr() < nsigmaTPCmin || ptrack.tpcNSigmaPr() > nsigmaTPCmax)
        return false;
      if (ntrack.tpcNSigmaPi() < nsigmaTPCmin || ntrack.tpcNSigmaPi() > nsigmaTPCmax)
        return false;

      // PID Selections (TOF)
      if (requireTOF) {
        if (ptrack.tofNSigmaPr() < nsigmaTOFmin || ptrack.tofNSigmaPr() > nsigmaTOFmax)
          return false;
        if (ntrack.tofNSigmaPi() < nsigmaTOFmin || ntrack.tofNSigmaPi() > nsigmaTOFmax)
          return false;
      }
    }

    // V0 Selections
    if (casc.v0cosPA(coll.posX(), coll.posY(), coll.posZ()) < v0cospaMin)
      return false;
    if (casc.v0radius() < minimumV0Radius || casc.v0radius() > maximumV0Radius)
      return false;
    if (casc.dcaV0daughters() > dcaV0DaughtersMax)
      return false;
    if (casc.dcapostopv() < dcapostoPVmin)
      return false;
    if (casc.dcanegtopv() < dcanegtoPVmin)
      return false;

    // Cascade Selections
    if (casc.cascradius() < minimumCascRadius || casc.cascradius() > maximumCascRadius)
      return false;
    if (casc.casccosPA(coll.posX(), coll.posY(), coll.posZ()) < casccospaMin)
      return false;
    if (casc.dcabachtopv() < dcabachtopvMin)
      return false;
    if (casc.dcav0topv(coll.posX(), coll.posY(), coll.posZ()) < dcaV0topvMin)
      return false;
    if (casc.dcacascdaughters() > dcaCascDaughtersMax)
      return false;

    // PID Selection on bachelor
    if (btrack.tpcNSigmaPi() < nsigmaTPCmin || btrack.tpcNSigmaPi() > nsigmaTPCmax)
      return false;

    // PID Selections (TOF)
    if (requireTOF) {
      if (btrack.tofNSigmaPi() < nsigmaTOFmin || btrack.tofNSigmaPi() > nsigmaTOFmax)
        return false;
    }

    // Rapidity Selection
    TLorentzVector lorentzVect;
    lorentzVect.SetXYZM(casc.px(), casc.py(), casc.pz(), 1.32171);
    if (lorentzVect.Rapidity() < yMin || lorentzVect.Rapidity() > yMax)
      return false;
    return true;
  }

  // Omega Selections
  template <typename CascType, typename T1, typename T2, typename T3, typename C>
  bool passedOmegaSelection(const CascType& casc, const T1& ptrack, const T2& ntrack, const T3& btrack, const C& coll)
  {
    if (!passedSingleTrackSelection(ptrack))
      return false;
    if (!passedSingleTrackSelection(ntrack))
      return false;
    if (!passedSingleTrackSelection(btrack))
      return false;

    // Omega+ Selection (Omega+ -> antiL + K+)
    if (btrack.sign() > 0) {
      if (ntrack.pt() < ptMin_V0_proton)
        return false;
      if (ntrack.pt() > ptMax_V0_proton)
        return false;
      if (ptrack.pt() < ptMin_V0_pion)
        return false;
      if (ptrack.pt() > ptMax_V0_pion)
        return false;

      // PID Selections (TPC)
      if (ntrack.tpcNSigmaPr() < nsigmaTPCmin || ntrack.tpcNSigmaPr() > nsigmaTPCmax)
        return false;
      if (ptrack.tpcNSigmaPi() < nsigmaTPCmin || ptrack.tpcNSigmaPi() > nsigmaTPCmax)
        return false;

      // PID Selections (TOF)
      if (requireTOF) {
        if (ntrack.tofNSigmaPr() < nsigmaTOFmin || ntrack.tofNSigmaPr() > nsigmaTOFmax)
          return false;
        if (ptrack.tofNSigmaPi() < nsigmaTOFmin || ptrack.tofNSigmaPi() > nsigmaTOFmax)
          return false;
      }
    }

    // Omega- Selection (Omega- -> L + K-)
    if (btrack.sign() < 0) {
      if (ptrack.pt() < ptMin_V0_proton)
        return false;
      if (ptrack.pt() > ptMax_V0_proton)
        return false;
      if (ntrack.pt() < ptMin_V0_pion)
        return false;
      if (ntrack.pt() > ptMax_V0_pion)
        return false;

      // PID Selections (TPC)
      if (ptrack.tpcNSigmaPr() < nsigmaTPCmin || ptrack.tpcNSigmaPr() > nsigmaTPCmax)
        return false;
      if (ntrack.tpcNSigmaPi() < nsigmaTPCmin || ntrack.tpcNSigmaPi() > nsigmaTPCmax)
        return false;

      // PID Selections (TOF)
      if (requireTOF) {
        if (ptrack.tofNSigmaPr() < nsigmaTOFmin || ptrack.tofNSigmaPr() > nsigmaTOFmax)
          return false;
        if (ntrack.tofNSigmaPi() < nsigmaTOFmin || ntrack.tofNSigmaPi() > nsigmaTOFmax)
          return false;
      }
    }

    // V0 Selections
    if (casc.v0cosPA(coll.posX(), coll.posY(), coll.posZ()) < v0cospaMin)
      return false;
    if (casc.v0radius() < minimumV0Radius || casc.v0radius() > maximumV0Radius)
      return false;
    if (casc.dcaV0daughters() > dcaV0DaughtersMax)
      return false;
    if (casc.dcapostopv() < dcapostoPVmin)
      return false;
    if (casc.dcanegtopv() < dcanegtoPVmin)
      return false;

    // Cascade Selections
    if (casc.cascradius() < minimumCascRadius || casc.cascradius() > maximumCascRadius)
      return false;
    if (casc.casccosPA(coll.posX(), coll.posY(), coll.posZ()) < casccospaMin)
      return false;
    if (casc.dcabachtopv() < dcabachtopvMin)
      return false;
    if (casc.dcav0topv(coll.posX(), coll.posY(), coll.posZ()) < dcaV0topvMin)
      return false;
    if (casc.dcacascdaughters() > dcaCascDaughtersMax)
      return false;

    // PID Selection on bachelor
    if (btrack.tpcNSigmaKa() < nsigmaTPCmin || btrack.tpcNSigmaKa() > nsigmaTPCmax)
      return false;

    // PID Selections (TOF)
    if (requireTOF) {
      if (btrack.tofNSigmaKa() < nsigmaTOFmin || btrack.tofNSigmaKa() > nsigmaTOFmax)
        return false;
    }

    // Rapidity Selection
    TLorentzVector lorentzVect;
    lorentzVect.SetXYZM(casc.px(), casc.py(), casc.pz(), 1.6724);
    if (lorentzVect.Rapidity() < yMin || lorentzVect.Rapidity() > yMax)
      return false;
    return true;
  }

  // Single-Track Selection
  template <typename T1>
  bool passedSingleTrackSelection(const T1& track)
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
    if (track.eta() < etaMin || track.eta() > etaMax)
      return false;
    if (requireTOF && (!track.hasTOF()))
      return false;
    return true;
  }

  float Minimum(float x1, float x2)
  {
    float x_min(x1);
    if (x1 < x2)
      x_min = x1;
    if (x1 >= x2)
      x_min = x2;

    return x_min;
  }

  double GetDeltaPhi(double a1, double a2)
  {
    double delta_phi(0);

    double phi1 = TVector2::Phi_0_2pi(a1);
    double phi2 = TVector2::Phi_0_2pi(a2);
    double diff = TMath::Abs(phi1 - phi2);

    if (diff <= TMath::Pi())
      delta_phi = diff;
    if (diff > TMath::Pi())
      delta_phi = TMath::TwoPi() - diff;

    return delta_phi;
  }

  void get_perpendicular_cone(TVector3 p, TVector3& u, float sign)
  {
    // Initialization
    float ux(0), uy(0), uz(0);

    // Components of Vector p
    float px = p.X();
    float py = p.Y();
    float pz = p.Z();

    // Protection 1
    if (px == 0 && py != 0) {

      uy = -(pz * pz) / py;
      ux = sign * sqrt(py * py - (pz * pz * pz * pz) / (py * py));
      uz = pz;
      u.SetXYZ(ux, uy, uz);
      return;
    }

    // Protection 2
    if (py == 0 && px != 0) {

      ux = -(pz * pz) / px;
      uy = sign * sqrt(px * px - (pz * pz * pz * pz) / (px * px));
      uz = pz;
      u.SetXYZ(ux, uy, uz);
      return;
    }

    // Equation Parameters
    float a = px * px + py * py;
    float b = 2.0 * px * pz * pz;
    float c = pz * pz * pz * pz - py * py * py * py - px * px * py * py;
    float delta = b * b - 4.0 * a * c;

    // Protection agains delta<0
    if (delta < 0) {
      return;
    }

    // Solutions
    ux = (-b + sign * sqrt(delta)) / (2.0 * a);
    uy = (-pz * pz - px * ux) / py;
    uz = pz;
    u.SetXYZ(ux, uy, uz);
    return;
  }

  void processData(SelectedCollisions::iterator const& collision, aod::V0Datas const& fullV0s, aod::CascDataExt const& Cascades, aod::V0sLinked const& /*V0linked*/, FullTracks const& tracks)
  {
    registryQC.fill(HIST("number_of_events_data"), 0.5);
    if (!collision.sel8())
      return;

    registryQC.fill(HIST("number_of_events_data"), 1.5);
    if (abs(collision.posZ()) > zVtx)
      return;

    registryQC.fill(HIST("number_of_events_data"), 2.5);

    // Find Leading Particle
    std::vector<int> particle_ID;
    int leading_ID(0);
    float ptMax(0);

    // Track Index
    int i = -1;
    for (auto track : tracks) {

      i++;
      if (!passedTrackSelectionForJets(track))
        continue;

      if (track.pt() > ptMax) {
        leading_ID = i;
        ptMax = track.pt();
      }
      particle_ID.push_back(i);
    }

    if (ptMax < ptLeadingMin)
      return;
    registryQC.fill(HIST("number_of_events_data"), 3.5);

    if (debug_level0 == true)
      return;

    auto const& leading_track = tracks.iteratorAt(leading_ID);
    TVector3 p_leading(leading_track.px(), leading_track.py(), leading_track.pz());
    int nParticles = static_cast<int>(particle_ID.size());

    // Jet Finder
    int exit(0);
    int nPartAssociated(0);
    do {
      // Initialization
      float distance_jet_min(1e+08);
      float distance_bkg_min(1e+08);
      int label_jet_particle(0);
      int i_jet_particle(0);

      for (int i = 0; i < nParticles; i++) {

        // Skip Leading Particle & Elements already associated to the Jet
        if (particle_ID[i] == leading_ID || particle_ID[i] == -1)
          continue;

        // Get Particle Momentum
        auto stored_track = tracks.iteratorAt(particle_ID[i]);
        TVector3 p_particle(stored_track.px(), stored_track.py(), stored_track.pz());

        // Variables
        float one_over_pt2_part = 1.0 / (p_particle.Pt() * p_particle.Pt());
        float one_over_pt2_lead = 1.0 / (p_leading.Pt() * p_leading.Pt());
        float deltaEta = p_particle.Eta() - p_leading.Eta();
        float deltaPhi = GetDeltaPhi(p_particle.Phi(), p_leading.Phi());
        float min = Minimum(one_over_pt2_part, one_over_pt2_lead);
        float Delta2 = deltaEta * deltaEta + deltaPhi * deltaPhi;

        // Distances
        float distance_jet = min * Delta2 / (Rjet * Rjet);
        float distance_bkg = one_over_pt2_part;

        // Find Minimum Distance Jet
        if (distance_jet < distance_jet_min) {
          distance_jet_min = distance_jet;
          label_jet_particle = particle_ID[i];
          i_jet_particle = i;
        }

        // Find Minimum Distance Bkg
        if (distance_bkg < distance_bkg_min) {
          distance_bkg_min = distance_bkg;
        }
      }

      if (distance_jet_min <= distance_bkg_min) {

        // Add Particle to Jet
        // jet_particle_ID.push_back(label_jet_particle);

        // Update Momentum of Leading Particle
        auto jet_track = tracks.iteratorAt(label_jet_particle);
        TVector3 p_i(jet_track.px(), jet_track.py(), jet_track.pz());
        p_leading = p_leading + p_i;

        // Remove Element
        particle_ID[i_jet_particle] = -1;
        nPartAssociated++;
      }

      if (nPartAssociated >= (nParticles - 1))
        exit = 1;
      if (distance_jet_min > distance_bkg_min)
        exit = 2;

    } while (exit == 0);

    if (debug_level1 == true)
      return;

    // Jet Axis
    TVector3 jet_axis(p_leading.X(), p_leading.Y(), p_leading.Z());

    // Cut events with jet not fully inside acceptance
    if ((abs(jet_axis.Eta()) + Rmax) > etaMax)
      return;
    registryQC.fill(HIST("number_of_events_data"), 4.5);

    // Perpendicular Cones for UE
    TVector3 ue_axis1(0.0, 0.0, 0.0);
    TVector3 ue_axis2(0.0, 0.0, 0.0);
    get_perpendicular_cone(jet_axis, ue_axis1, +1.0);
    get_perpendicular_cone(jet_axis, ue_axis2, -1.0);

    // Protection against delta<0
    if (ue_axis1.X() == 0 && ue_axis1.Y() == 0 && ue_axis1.Z() == 0)
      return;
    if (ue_axis2.X() == 0 && ue_axis2.Y() == 0 && ue_axis2.Z() == 0)
      return;
    registryQC.fill(HIST("number_of_events_data"), 5.5);

    if (debug_level2 == true)
      return;

    // Event multiplicity
    float multiplicity = collision.centFT0M();

    if (debug_level3 == true)
      return;

    for (auto& v0 : fullV0s) {

      const auto& pos = v0.posTrack_as<FullTracks>();
      const auto& neg = v0.negTrack_as<FullTracks>();

      if (debug_level4 == true)
        return;

      TVector3 v0dir(pos.px() + neg.px(), pos.py() + neg.py(), pos.pz() + neg.pz());

      float deltaEta_jet = v0dir.Eta() - jet_axis.Eta();
      float deltaPhi_jet = GetDeltaPhi(v0dir.Phi(), jet_axis.Phi());
      float deltaR_jet = sqrt(deltaEta_jet * deltaEta_jet + deltaPhi_jet * deltaPhi_jet);

      float deltaEta_ue1 = v0dir.Eta() - ue_axis1.Eta();
      float deltaPhi_ue1 = GetDeltaPhi(v0dir.Phi(), ue_axis1.Phi());
      float deltaR_ue1 = sqrt(deltaEta_ue1 * deltaEta_ue1 + deltaPhi_ue1 * deltaPhi_ue1);

      float deltaEta_ue2 = v0dir.Eta() - ue_axis2.Eta();
      float deltaPhi_ue2 = GetDeltaPhi(v0dir.Phi(), ue_axis2.Phi());
      float deltaR_ue2 = sqrt(deltaEta_ue2 * deltaEta_ue2 + deltaPhi_ue2 * deltaPhi_ue2);

      // K0s
      if (passedK0ShortSelection(v0, pos, neg)) {
        if (debug_level5 == true)
          return;

        if (deltaR_jet < Rmax) {
          registryData.fill(HIST("K0s_in_jet"), multiplicity, v0.pt(), v0.mK0Short());
        }
        if (deltaR_ue1 < Rmax || deltaR_ue2 < Rmax) {
          registryData.fill(HIST("K0s_in_ue"), multiplicity, v0.pt(), v0.mK0Short());
        }
      }

      // Lambda
      if (passedLambdaSelection(v0, pos, neg)) {
        if (debug_level6 == true)
          return;

        if (deltaR_jet < Rmax) {
          registryData.fill(HIST("Lambda_in_jet"), multiplicity, v0.pt(), v0.mLambda());
        }

        if (deltaR_ue1 < Rmax || deltaR_ue2 < Rmax) {
          registryData.fill(HIST("Lambda_in_ue"), multiplicity, v0.pt(), v0.mLambda());
        }
      }

      // AntiLambda
      if (passedAntiLambdaSelection(v0, pos, neg)) {
        if (debug_level7 == true)
          return;

        if (deltaR_jet < Rmax) {
          registryData.fill(HIST("AntiLambda_in_jet"), multiplicity, v0.pt(), v0.mAntiLambda());
        }

        if (deltaR_ue1 < Rmax || deltaR_ue2 < Rmax) {
          registryData.fill(HIST("AntiLambda_in_ue"), multiplicity, v0.pt(), v0.mAntiLambda());
        }
      }
    }

    if (debug_level8 == true)
      return;

    // Cascades
    for (auto& casc : Cascades) {

      auto bach = casc.bachelor_as<FullTracks>();
      auto pos = casc.posTrack_as<FullTracks>();
      auto neg = casc.negTrack_as<FullTracks>();

      if (debug_level9 == true)
        return;

      TVector3 cascade_dir(casc.px(), casc.py(), casc.pz());

      float deltaEta_jet = cascade_dir.Eta() - jet_axis.Eta();
      float deltaPhi_jet = GetDeltaPhi(cascade_dir.Phi(), jet_axis.Phi());
      float deltaR_jet = sqrt(deltaEta_jet * deltaEta_jet + deltaPhi_jet * deltaPhi_jet);
      float deltaEta_ue1 = cascade_dir.Eta() - ue_axis1.Eta();
      float deltaPhi_ue1 = GetDeltaPhi(cascade_dir.Phi(), ue_axis1.Phi());
      float deltaR_ue1 = sqrt(deltaEta_ue1 * deltaEta_ue1 + deltaPhi_ue1 * deltaPhi_ue1);
      float deltaEta_ue2 = cascade_dir.Eta() - ue_axis2.Eta();
      float deltaPhi_ue2 = GetDeltaPhi(cascade_dir.Phi(), ue_axis2.Phi());
      float deltaR_ue2 = sqrt(deltaEta_ue2 * deltaEta_ue2 + deltaPhi_ue2 * deltaPhi_ue2);

      // Xi+
      if (passedXiSelection(casc, pos, neg, bach, collision) && bach.sign() > 0) {
        if (debug_level10 == true)
          return;

        if (deltaR_jet < Rmax) {
          registryData.fill(HIST("XiPos_in_jet"), multiplicity, casc.pt(), casc.mXi());
        }
        if (deltaR_ue1 < Rmax || deltaR_ue2 < Rmax) {
          registryData.fill(HIST("XiPos_in_ue"), multiplicity, casc.pt(), casc.mXi());
        }
      }
      // Xi-
      if (passedXiSelection(casc, pos, neg, bach, collision) && bach.sign() < 0) {
        if (debug_level11 == true)
          return;

        if (deltaR_jet < Rmax) {
          registryData.fill(HIST("XiNeg_in_jet"), multiplicity, casc.pt(), casc.mXi());
        }
        if (deltaR_ue1 < Rmax || deltaR_ue2 < Rmax) {
          registryData.fill(HIST("XiNeg_in_ue"), multiplicity, casc.pt(), casc.mXi());
        }
      }

      // Omega+
      if (passedOmegaSelection(casc, pos, neg, bach, collision) && bach.sign() > 0) {
        if (debug_level12 == true)
          return;

        if (deltaR_jet < Rmax) {
          registryData.fill(HIST("OmegaPos_in_jet"), multiplicity, casc.pt(), casc.mOmega());
        }
        if (deltaR_ue1 < Rmax || deltaR_ue2 < Rmax) {
          registryData.fill(HIST("OmegaPos_in_ue"), multiplicity, casc.pt(), casc.mOmega());
        }
      }
      // Omega-
      if (passedOmegaSelection(casc, pos, neg, bach, collision) && bach.sign() < 0) {
        if (debug_level13 == true)
          return;

        if (deltaR_jet < Rmax) {
          registryData.fill(HIST("OmegaNeg_in_jet"), multiplicity, casc.pt(), casc.mOmega());
        }
        if (deltaR_ue1 < Rmax || deltaR_ue2 < Rmax) {
          registryData.fill(HIST("OmegaNeg_in_ue"), multiplicity, casc.pt(), casc.mOmega());
        }
      }
    }

    for (auto track : tracks) {

      if (!passedTrackSelectionForPions(track))
        continue;

      TVector3 track_dir(track.px(), track.py(), track.pz());
      float deltaEta_jet = track_dir.Eta() - jet_axis.Eta();
      float deltaPhi_jet = GetDeltaPhi(track_dir.Phi(), jet_axis.Phi());
      float deltaR_jet = sqrt(deltaEta_jet * deltaEta_jet + deltaPhi_jet * deltaPhi_jet);
      float deltaEta_ue1 = track_dir.Eta() - ue_axis1.Eta();
      float deltaPhi_ue1 = GetDeltaPhi(track_dir.Phi(), ue_axis1.Phi());
      float deltaR_ue1 = sqrt(deltaEta_ue1 * deltaEta_ue1 + deltaPhi_ue1 * deltaPhi_ue1);
      float deltaEta_ue2 = track_dir.Eta() - ue_axis2.Eta();
      float deltaPhi_ue2 = GetDeltaPhi(track_dir.Phi(), ue_axis2.Phi());
      float deltaR_ue2 = sqrt(deltaEta_ue2 * deltaEta_ue2 + deltaPhi_ue2 * deltaPhi_ue2);

      // TPC
      if (deltaR_jet < Rmax && track.sign() > 0) {
        registryData.fill(HIST("piplus_tpc_in_jet"), multiplicity, track.pt(), track.tpcNSigmaPi());
      }
      if ((deltaR_ue1 < Rmax || deltaR_ue2) < Rmax && track.sign() > 0) {
        registryData.fill(HIST("piplus_tpc_in_ue"), multiplicity, track.pt(), track.tpcNSigmaPi());
      }
      if (deltaR_jet < Rmax && track.sign() < 0) {
        registryData.fill(HIST("piminus_tpc_in_jet"), multiplicity, track.pt(), track.tpcNSigmaPi());
      }
      if ((deltaR_ue1 < Rmax || deltaR_ue2) < Rmax && track.sign() < 0) {
        registryData.fill(HIST("piminus_tpc_in_ue"), multiplicity, track.pt(), track.tpcNSigmaPi());
      }

      if (track.tpcNSigmaPi() < nsigmaTPCmin || track.tpcNSigmaPi() > nsigmaTPCmax)
        continue;

      if (!track.hasTOF())
        continue;

      // TOF
      if (deltaR_jet < Rmax && track.sign() > 0) {
        registryData.fill(HIST("piplus_tof_in_jet"), multiplicity, track.pt(), track.tofNSigmaPi());
      }
      if ((deltaR_ue1 < Rmax || deltaR_ue2) < Rmax && track.sign() > 0) {
        registryData.fill(HIST("piplus_tof_in_ue"), multiplicity, track.pt(), track.tofNSigmaPi());
      }
      if (deltaR_jet < Rmax && track.sign() < 0) {
        registryData.fill(HIST("piminus_tof_in_jet"), multiplicity, track.pt(), track.tofNSigmaPi());
      }
      if ((deltaR_ue1 < Rmax || deltaR_ue2) < Rmax && track.sign() < 0) {
        registryData.fill(HIST("piminus_tof_in_ue"), multiplicity, track.pt(), track.tofNSigmaPi());
      }
    }
  }
  PROCESS_SWITCH(strangeness_in_jets, processData, "Process data", true);

  Preslice<aod::V0Datas> perCollisionV0 = o2::aod::v0data::collisionId;
  Preslice<aod::CascDataExt> perCollisionCasc = o2::aod::cascade::collisionId;
  Preslice<aod::McParticles> perMCCollision = o2::aod::mcparticle::mcCollisionId;

  void processMCefficiency(SimCollisions const& collisions, MCTracks const& /*mcTracks*/, aod::V0Datas const& fullV0s, aod::CascDataExt const& Cascades, aod::McCollisions const& /*mcCollisions*/, const aod::McParticles& mcParticles)
  {
    for (const auto& collision : collisions) {
      registryQC.fill(HIST("number_of_events_mc"), 0.5);
      if (!collision.sel8())
        continue;

      registryQC.fill(HIST("number_of_events_mc"), 1.5);
      if (abs(collision.posZ()) > 10.0)
        continue;

      registryQC.fill(HIST("number_of_events_mc"), 2.5);
      float multiplicity = collision.centFT0M();

      auto v0s_per_coll = fullV0s.sliceBy(perCollisionV0, collision.globalIndex());
      auto casc_per_coll = Cascades.sliceBy(perCollisionCasc, collision.globalIndex());
      auto mcParticles_per_coll = mcParticles.sliceBy(perMCCollision, collision.globalIndex());

      for (auto& v0 : v0s_per_coll) {

        const auto& pos = v0.posTrack_as<MCTracks>();
        const auto& neg = v0.negTrack_as<MCTracks>();
        if (!pos.has_mcParticle())
          continue;
        if (!neg.has_mcParticle())
          continue;

        auto posParticle = pos.mcParticle_as<aod::McParticles>();
        auto negParticle = neg.mcParticle_as<aod::McParticles>();
        if (!posParticle.has_mothers())
          continue;
        if (!negParticle.has_mothers())
          continue;

        int pdg_parent(0);
        bool isPhysPrim = false;
        for (auto& particleMotherOfNeg : negParticle.mothers_as<aod::McParticles>()) {
          for (auto& particleMotherOfPos : posParticle.mothers_as<aod::McParticles>()) {
            if (particleMotherOfNeg.globalIndex() == particleMotherOfPos.globalIndex()) {
              pdg_parent = particleMotherOfNeg.pdgCode();
              isPhysPrim = particleMotherOfNeg.isPhysicalPrimary();
            }
          }
        }
        if (pdg_parent == 0)
          continue;

        // K0s
        if (passedK0ShortSelection(v0, pos, neg) && pdg_parent == 310) {
          registryMC.fill(HIST("K0s_reconstructed_incl"), multiplicity, v0.pt());
        }
        if (passedLambdaSelection(v0, pos, neg) && pdg_parent == 3122) {
          registryMC.fill(HIST("Lambda_reconstructed_incl"), multiplicity, v0.pt());
        }
        if (passedAntiLambdaSelection(v0, pos, neg) && pdg_parent == -3122) {
          registryMC.fill(HIST("AntiLambda_reconstructed_incl"), multiplicity, v0.pt());
        }
        if (!isPhysPrim)
          continue;

        // K0s
        if (passedK0ShortSelection(v0, pos, neg) && pdg_parent == 310) {
          registryMC.fill(HIST("K0s_reconstructed"), multiplicity, v0.pt());
        }
        if (passedLambdaSelection(v0, pos, neg) && pdg_parent == 3122) {
          registryMC.fill(HIST("Lambda_reconstructed"), multiplicity, v0.pt());
        }
        if (passedAntiLambdaSelection(v0, pos, neg) && pdg_parent == -3122) {
          registryMC.fill(HIST("AntiLambda_reconstructed"), multiplicity, v0.pt());
        }
      }

      // Cascades
      for (auto& casc : casc_per_coll) {
        auto bach = casc.template bachelor_as<MCTracks>();
        auto neg = casc.template negTrack_as<MCTracks>();
        auto pos = casc.template posTrack_as<MCTracks>();

        if (!bach.has_mcParticle())
          continue;
        if (!pos.has_mcParticle())
          continue;
        if (!neg.has_mcParticle())
          continue;

        auto posParticle = pos.mcParticle_as<aod::McParticles>();
        auto negParticle = neg.mcParticle_as<aod::McParticles>();
        auto bachParticle = bach.mcParticle_as<aod::McParticles>();
        if (!posParticle.has_mothers())
          continue;
        if (!negParticle.has_mothers())
          continue;
        if (!bachParticle.has_mothers())
          continue;

        int pdg_parent(0);
        for (auto& particleMotherOfNeg : negParticle.mothers_as<aod::McParticles>()) {
          for (auto& particleMotherOfPos : posParticle.mothers_as<aod::McParticles>()) {
            for (auto& particleMotherOfBach : bachParticle.mothers_as<aod::McParticles>()) {
              if (particleMotherOfNeg.globalIndex() != particleMotherOfPos.globalIndex())
                continue;
              if (abs(particleMotherOfNeg.pdgCode()) != 3122)
                continue;
              if (!particleMotherOfBach.isPhysicalPrimary())
                continue;
              pdg_parent = particleMotherOfBach.pdgCode();
            }
          }
        }
        if (pdg_parent == 0)
          continue;

        // Xi+
        if (passedXiSelection(casc, pos, neg, bach, collision) && pdg_parent == -3312) {
          registryMC.fill(HIST("XiPos_reconstructed"), multiplicity, casc.pt());
        }
        // Xi-
        if (passedXiSelection(casc, pos, neg, bach, collision) && pdg_parent == 3312) {
          registryMC.fill(HIST("XiNeg_reconstructed"), multiplicity, casc.pt());
        }
        // Omega+
        if (passedOmegaSelection(casc, pos, neg, bach, collision) && pdg_parent == -3334) {
          registryMC.fill(HIST("OmegaPos_reconstructed"), multiplicity, casc.pt());
        }
        // Omega-
        if (passedOmegaSelection(casc, pos, neg, bach, collision) && pdg_parent == 3334) {
          registryMC.fill(HIST("OmegaNeg_reconstructed"), multiplicity, casc.pt());
        }
      }

      for (auto& mcParticle : mcParticles_per_coll) {

        if (mcParticle.y() < yMin || mcParticle.y() > yMax)
          continue;
        if (!mcParticle.isPhysicalPrimary())
          continue;

        // K0s
        if (mcParticle.pdgCode() == 310) {
          registryMC.fill(HIST("K0s_generated"), multiplicity, mcParticle.pt());
          registryMC.fill(HIST("K0s_eta_pt_pythia"), mcParticle.pt(), mcParticle.eta());
        }
        // Lambda
        if (mcParticle.pdgCode() == 3122) {
          registryMC.fill(HIST("Lambda_generated"), multiplicity, mcParticle.pt());
          registryMC.fill(HIST("Lambda_eta_pt_pythia"), mcParticle.pt(), mcParticle.eta());
        }
        // AntiLambda
        if (mcParticle.pdgCode() == -3122) {
          registryMC.fill(HIST("AntiLambda_generated"), multiplicity, mcParticle.pt());
          registryMC.fill(HIST("AntiLambda_eta_pt_pythia"), mcParticle.pt(), mcParticle.eta());
        }
        // Xi Pos
        if (mcParticle.pdgCode() == -3312) {
          registryMC.fill(HIST("XiPos_generated"), multiplicity, mcParticle.pt());
        }
        // Xi Neg
        if (mcParticle.pdgCode() == 3312) {
          registryMC.fill(HIST("XiNeg_generated"), multiplicity, mcParticle.pt());
        }
        // Omega Pos
        if (mcParticle.pdgCode() == -3334) {
          registryMC.fill(HIST("OmegaPos_generated"), multiplicity, mcParticle.pt());
        }
        // Omega Neg
        if (mcParticle.pdgCode() == 3334) {
          registryMC.fill(HIST("OmegaNeg_generated"), multiplicity, mcParticle.pt());
        }
      }
    }
  }
  PROCESS_SWITCH(strangeness_in_jets, processMCefficiency, "Process MC Efficiency", false);

  void processWeights(SimCollisions const& collisions, aod::McParticles const& mcParticles)
  {
    // Loop over MC Collisions
    for (const auto& collision : collisions) {

      // Selection on z_{vertex}
      if (abs(collision.posZ()) > 10)
        continue;

      // MC Particles per Collision
      auto mcParticles_per_coll = mcParticles.sliceBy(perMCCollision, collision.globalIndex());

      std::vector<int> particle_ID;
      int leading_ID = 0;
      float pt_max(0);

      for (auto& particle : mcParticles_per_coll) {

        // Global Index
        int i = particle.globalIndex();

        // Select Primary Particles
        float deltaX = particle.vx() - collision.posX();
        float deltaY = particle.vy() - collision.posY();
        float deltaZ = particle.vz() - collision.posZ();
        float deltaR = sqrt(deltaX * deltaX + deltaY * deltaY + deltaZ * deltaZ);
        if (deltaR > 0.1)
          continue;

        // Pseudorapidity Selection
        if (abs(particle.eta()) > 0.8)
          continue;

        // PDG Selection
        int pdg = abs(particle.pdgCode());
        if ((pdg != 11) && (pdg != 211) && (pdg != 321) && (pdg != 2212))
          continue;

        // Find pt Leading
        if (particle.pt() > pt_max) {
          leading_ID = i;
          pt_max = particle.pt();
        }

        // Store Array Element
        particle_ID.push_back(i);
      }

      // Skip Events with pt<pt_leading_min
      if (pt_max < ptLeadingMin)
        continue;

      // Number of Stored Particles
      int nParticles = static_cast<int>(particle_ID.size());

      // Momentum of the Leading Particle
      auto const& leading_track = mcParticles_per_coll.iteratorAt(leading_ID);
      TVector3 p_leading(leading_track.px(), leading_track.py(), leading_track.pz());

      // Labels
      int exit(0);
      int nPartAssociated(0);

      // Jet Finder
      do {
        // Initialization
        float distance_jet_min(1e+08);
        float distance_bkg_min(1e+08);
        int label_jet_particle(0);
        int i_jet_particle(0);

        for (int i = 0; i < nParticles; i++) {

          // Skip Leading Particle & Elements already associated to the Jet
          if (particle_ID[i] == leading_ID || particle_ID[i] == -1)
            continue;

          // Get Particle Momentum
          auto stored_track = mcParticles_per_coll.iteratorAt(particle_ID[i]);
          TVector3 p_particle(stored_track.px(), stored_track.py(), stored_track.pz());

          // Variables
          float one_over_pt2_part = 1.0 / (p_particle.Pt() * p_particle.Pt());
          float one_over_pt2_lead = 1.0 / (p_leading.Pt() * p_leading.Pt());
          float deltaEta = p_particle.Eta() - p_leading.Eta();
          float deltaPhi = GetDeltaPhi(p_particle.Phi(), p_leading.Phi());
          float min = Minimum(one_over_pt2_part, one_over_pt2_lead);
          float Delta2 = deltaEta * deltaEta + deltaPhi * deltaPhi;

          // Distances
          float distance_jet = min * Delta2 / (Rjet * Rjet);
          float distance_bkg = one_over_pt2_part;

          // Find Minimum Distance Jet
          if (distance_jet < distance_jet_min) {
            distance_jet_min = distance_jet;
            label_jet_particle = particle_ID[i];
            i_jet_particle = i;
          }

          // Find Minimum Distance Bkg
          if (distance_bkg < distance_bkg_min) {
            distance_bkg_min = distance_bkg;
          }
        }

        if (distance_jet_min <= distance_bkg_min) {

          // Add Particle to Jet
          // jet_particle_ID.push_back(label_jet_particle);

          // Update Momentum of Leading Particle
          auto jet_track = mcParticles_per_coll.iteratorAt(label_jet_particle);
          TVector3 p_i(jet_track.px(), jet_track.py(), jet_track.pz());
          p_leading = p_leading + p_i;

          // Remove Element
          particle_ID[i_jet_particle] = -1;
          nPartAssociated++;
        }

        if (nPartAssociated >= (nParticles - 1))
          exit = 1;
        if (distance_jet_min > distance_bkg_min)
          exit = 2;

      } while (exit == 0);

      // Jet Axis
      TVector3 jet_axis(p_leading.X(), p_leading.Y(), p_leading.Z());

      if ((abs(jet_axis.Eta()) + Rmax) > etaMax)
        return;

      // Perpendicular Cones for UE
      TVector3 ue_axis1(0.0, 0.0, 0.0);
      TVector3 ue_axis2(0.0, 0.0, 0.0);
      get_perpendicular_cone(jet_axis, ue_axis1, +1.0);
      get_perpendicular_cone(jet_axis, ue_axis2, -1.0);

      // Protection against delta<0
      if (ue_axis1.X() == 0 && ue_axis1.Y() == 0 && ue_axis1.Z() == 0)
        return;
      if (ue_axis2.X() == 0 && ue_axis2.Y() == 0 && ue_axis2.Z() == 0)
        return;

      // Generated Particles
      for (auto& particle : mcParticles_per_coll) {

        // PDG Selection
        int pdg = particle.pdgCode();

        if (!particle.isPhysicalPrimary())
          continue;
        if (particle.y() < yMin || particle.y() > yMax)
          continue;

        TVector3 p_particle(particle.px(), particle.py(), particle.pz());
        float deltaEta_jet = p_particle.Eta() - jet_axis.Eta();
        float deltaPhi_jet = GetDeltaPhi(p_particle.Phi(), jet_axis.Phi());
        float deltaR_jet = sqrt(deltaEta_jet * deltaEta_jet + deltaPhi_jet * deltaPhi_jet);
        float deltaEta_ue1 = p_particle.Eta() - ue_axis1.Eta();
        float deltaPhi_ue1 = GetDeltaPhi(p_particle.Phi(), ue_axis1.Phi());
        float deltaR_ue1 = sqrt(deltaEta_ue1 * deltaEta_ue1 + deltaPhi_ue1 * deltaPhi_ue1);
        float deltaEta_ue2 = p_particle.Eta() - ue_axis2.Eta();
        float deltaPhi_ue2 = GetDeltaPhi(p_particle.Phi(), ue_axis2.Phi());
        float deltaR_ue2 = sqrt(deltaEta_ue2 * deltaEta_ue2 + deltaPhi_ue2 * deltaPhi_ue2);

        // Fill K0s
        if (pdg == 310) {
          if (deltaR_jet < Rmax)
            registryMC.fill(HIST("K0s_eta_pt_jet"), particle.pt(), particle.eta());
          if (deltaR_ue1 < Rmax || deltaR_ue2 < Rmax)
            registryMC.fill(HIST("K0s_eta_pt_ue"), particle.pt(), particle.eta());
        }

        // Fill Lambda
        if (pdg == 3122) {
          if (deltaR_jet < Rmax)
            registryMC.fill(HIST("Lambda_eta_pt_jet"), particle.pt(), particle.eta());
          if (deltaR_ue1 < Rmax || deltaR_ue2 < Rmax)
            registryMC.fill(HIST("Lambda_eta_pt_ue"), particle.pt(), particle.eta());
        }

        // Fill AntiLambda
        if (pdg == -3122) {
          if (deltaR_jet < Rmax)
            registryMC.fill(HIST("AntiLambda_eta_pt_jet"), particle.pt(), particle.eta());
          if (deltaR_ue1 < Rmax || deltaR_ue2 < Rmax)
            registryMC.fill(HIST("AntiLambda_eta_pt_ue"), particle.pt(), particle.eta());
        }
      }
    }
  }
  PROCESS_SWITCH(strangeness_in_jets, processWeights, "Process Weights", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<strangeness_in_jets>(cfgc)};
}
