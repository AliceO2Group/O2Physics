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
/// \file phik0sanalysis.cxx
/// \brief Analysis task for the Phi and K0S rapidity correlations analysis
/// \author Stefano Cannito (stefano.cannito@cern.ch)

#include <TH1F.h>
#include <TRandom.h>
#include <TDirectory.h>
#include <THn.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TFile.h>
#include <TH2F.h>
#include <TLorentzVector.h>
#include <TPDGCode.h>
#include <Math/Vector4D.h>
#include <array>
#include <vector>
#include <cmath>
#include <cstdlib>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Common/DataModel/EventSelection.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "Common/DataModel/PIDResponse.h"
#include "Framework/ASoAHelpers.h"
#include "CommonConstants/PhysicsConstants.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/HistogramRegistry.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Core/trackUtilities.h"
#include "Common/Core/TrackSelection.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "PWGLF/Utils/inelGt.h"
#include "PWGLF/DataModel/mcCentrality.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct phik0shortanalysis {
  // Histograms are defined with HistogramRegistry
  HistogramRegistry eventHist{"eventHist", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry MCeventHist{"MCeventHist", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry PhicandHist{"PhicandHist", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry K0SHist{"K0SHist", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry PhipurHist{"PhipurHist", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry PhiK0SHist{"PhiK0SHist", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry MCPhiK0SHist{"MCPhiK0SHist", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry PhiPionHist{"PhiPionHist", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry MCPhiPionHist{"MCPhiPionHist", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry PhieffHist{"PhieffHist", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry K0SeffHist{"K0SeffHist", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry PioneffHist{"PioneffHist", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry closureMCPhipurHist{"closureMCPhipurHist", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry closureMCPhiK0SHist{"closureMCPhiK0SHist", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry closureMCPhiPionHist{"closureMCPhiPionHist", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  // Configurable for event selection
  Configurable<float> cutzvertex{"cutzvertex", 10.0f, "Accepted z-vertex range (cm)"};

  // Configurable on multiplicity bins
  Configurable<std::vector<double>> binsMult{"binsMult", {0.0, 1.0, 5.0, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 70.0, 100.0}, "Multiplicity bin limits"};

  // Configurables for V0 selection
  Configurable<int> minTPCnClsFound{"minTPCnClsFound", 70, "min number of found TPC clusters"};
  Configurable<int> minNCrossedRowsTPC{"minNCrossedRowsTPC", 80, "min number of TPC crossed rows"};
  Configurable<float> maxChi2TPC{"maxChi2TPC", 4.0f, "max chi2 per cluster TPC"};
  Configurable<float> etaMax{"etaMax", 0.8f, "eta max"};

  Configurable<float> v0setting_cospa{"v0setting_cospa", 0.98, "V0 CosPA"};
  Configurable<float> v0setting_radius{"v0setting_radius", 0.5, "v0radius"};
  Configurable<float> v0setting_dcav0dau{"v0setting_dcav0dau", 1, "DCA V0 Daughters"};
  Configurable<float> v0setting_dcapostopv{"v0setting_dcapostopv", 0.06, "DCA Pos To PV"};
  Configurable<float> v0setting_dcanegtopv{"v0setting_dcanegtopv", 0.06, "DCA Neg To PV"};
  Configurable<float> NSigmaTPCPion{"NSigmaTPCPion", 4.0, "NSigmaTPCPion"};

  // Configurables on K0S mass
  Configurable<float> lowmK0S{"lowmK0S", 0.48, "Lower limit on K0Short mass"};
  Configurable<float> upmK0S{"upmK0S", 0.52, "Upper limit on K0Short mass"};

  // Configurable on K0S pT bins
  Configurable<std::vector<double>> binspTK0S{"binspTK0S", {0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0}, "pT bin limits for K0S"};

  // Configurables on Phi mass
  Configurable<int> nBinsmPhi{"nBinsmPhi", 13, "N bins in cfgPhimassaxis"};
  Configurable<float> lowmPhi{"lowmPhiMB", 1.0095, "Upper limits on Phi mass for signal extraction"};
  Configurable<float> upmPhi{"upmPhiMB", 1.029, "Upper limits on Phi mass for signal extraction"};

  // Configurables for phi selection
  Configurable<float> cfgCutCharge{"cfgCutCharge", 0.0, "Cut on charge"};
  Configurable<bool> cfgPrimaryTrack{"cfgPrimaryTrack", false, "Primary track selection"};
  Configurable<bool> cfgGlobalWoDCATrack{"cfgGlobalWoDCATrack", true, "Global track selection without DCA"};
  Configurable<bool> cfgPVContributor{"cfgPVContributor", true, "PV contributor track selection"};
  Configurable<float> cMinKaonPtcut{"cMinKaonPtcut", 0.15f, "Track minimum pt cut"};
  Configurable<float> cMaxDCAzToPVcut{"cMaxDCAzToPVcut", 2.0f, "Track DCAz cut to PV Maximum"};
  Configurable<float> cMaxDCArToPV1{"cMaxDCArToPV1", 0.004f, "Track DCAr cut to PV config 1"};
  Configurable<float> cMaxDCArToPV2{"cMaxDCArToPV2", 0.013f, "Track DCAr cut to PV config 2"};
  Configurable<float> cMaxDCArToPV3{"cMaxDCArToPV3", 1.0f, "Track DCAr cut to PV config 3"};

  Configurable<bool> isNoTOF{"isNoTOF", false, "isNoTOF"};
  Configurable<float> nsigmaCutTPCKa{"nsigmacutTPC", 3.0, "Value of the TPC Nsigma cut"};
  Configurable<float> nsigmaCutCombinedKa{"nsigmaCutCombined", 3.0, "Value of the TOF Nsigma cut"};

  // Configurables for pions selection(extra with respect to a few of those defined for V0)
  Configurable<int> minITSnCls{"minITSnCls", 4, "min number of ITS clusters"};
  Configurable<float> maxChi2ITS{"maxChi2ITS", 36.0f, "max chi2 per cluster ITS"};
  Configurable<float> dcaxyMax{"dcaxyMax", 0.1f, "Maximum DCAxy to primary vertex"};
  Configurable<float> dcazMax{"dcazMax", 0.1f, "Maximum DCAz to primary vertex"};
  Configurable<float> NSigmaTOFPion{"NSigmaTOFPion", 5.0, "NSigmaTOFPion"};

  // Configurable on pion pT bins
  Configurable<std::vector<double>> binspTPi{"binspTPi", {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.5, 2.0, 3.0}, "pT bin limits for pions"};

  // Configurables for delta y selection
  Configurable<int> nBinsy{"nBinsy", 80, "Number of bins in y and deltay axis"};
  Configurable<float> cfgyAcceptance{"cfgyAcceptance", 0.5f, "Rapidity acceptance"};
  Configurable<float> cfgFirstCutonDeltay{"cgfFirstCutonDeltay", 0.5f, "First upper bound on Deltay selection"};
  Configurable<float> cfgSecondCutonDeltay{"cgfSecondCutonDeltay", 0.1f, "Second upper bound on Deltay selection"};
  Configurable<float> cfgyAcceptanceSmear{"cfgyAcceptanceSmear", 0.8f, "Rapidity acceptance for smearing matrix study"};

  // Configurable for RecMC
  Configurable<bool> cfgiskNoITSROFrameBorder{"cfgiskNoITSROFrameBorder", false, "kNoITSROFrameBorder request on RecMC collisions"};

  // Constants
  double massKa = o2::constants::physics::MassKPlus;
  double massPi = o2::constants::physics::MassPiPlus;

  // Defining filters for events (event selection)
  // Processed events will be already fulfilling the event selection requirements
  Filter eventFilter = (o2::aod::evsel::sel8 == true);
  Filter posZFilter = (nabs(o2::aod::collision::posZ) < cutzvertex);

  // Defining filters on V0s (cannot filter on dynamic columns)
  Filter preFilterV0 = (nabs(aod::v0data::dcapostopv) > v0setting_dcapostopv && nabs(aod::v0data::dcanegtopv) > v0setting_dcanegtopv && aod::v0data::dcaV0daughters < v0setting_dcav0dau);

  // Defining the type of the collisions for data and MC
  using SelCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::PVMults>;
  using SimCollisions = soa::Join<SelCollisions, aod::McCollisionLabels>;
  using MCCollisions = soa::Join<aod::McCollisions, aod::McCentFT0Ms>;

  // Defining the type of the V0s
  using FullV0s = soa::Filtered<aod::V0Datas>;
  using FullMCV0s = soa::Join<FullV0s, aod::McV0Labels>;

  // Defining the type of the tracks for data and MC
  using FullTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTOFFullPi, aod::pidTOFFullKa>;
  using FullMCTracks = soa::Join<FullTracks, aod::McTrackLabels>;

  using V0DauTracks = soa::Join<aod::TracksIU, aod::TracksExtra, aod::pidTPCFullPi>;
  using V0DauMCTracks = soa::Join<V0DauTracks, aod::McTrackLabels>;

  // Defining the binning policy for mixed event
  using BinningTypeVertexContributor = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0M>;

  SliceCache cache;

  Partition<FullTracks> posTracks = aod::track::signed1Pt > cfgCutCharge;
  Partition<FullTracks> negTracks = aod::track::signed1Pt < cfgCutCharge;

  Partition<FullMCTracks> posMCTracks = aod::track::signed1Pt > cfgCutCharge;
  Partition<FullMCTracks> negMCTracks = aod::track::signed1Pt < cfgCutCharge;

  // Necessary to flag INEL>0 events in GenMC
  Service<o2::framework::O2DatabasePDG> pdgDB;

  void init(InitContext&)
  {
    // Axes
    AxisSpec K0SmassAxis = {200, 0.45f, 0.55f, "#it{M}_{inv} [GeV/#it{c}^{2}]"};
    AxisSpec PhimassAxis = {200, 0.9f, 1.2f, "#it{M}_{inv} [GeV/#it{c}^{2}]"};
    AxisSpec sigPhimassAxis = {nBinsmPhi, lowmPhi, upmPhi, "#it{M}_{inv} [GeV/#it{c}^{2}]"};
    AxisSpec vertexZAxis = {100, -15.f, 15.f, "vrtx_{Z} [cm]"};
    AxisSpec yAxis = {nBinsy, -cfgyAcceptanceSmear, cfgyAcceptanceSmear, "#it{y}"};
    AxisSpec deltayAxis = {nBinsy, 0.0f, 1.6f, "|#it{#Deltay}|"};
    AxisSpec multAxis = {120, 0.0f, 120.0f, "centFT0M"};
    AxisSpec binnedmultAxis{(std::vector<double>)binsMult, "centFT0M"};
    AxisSpec ptK0SAxis = {60, 0.0f, 6.0f, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec binnedptK0SAxis{(std::vector<double>)binspTK0S, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec ptPiAxis = {30, 0.0f, 3.0f, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec binnedptPiAxis{(std::vector<double>)binspTPi, "#it{p}_{T} (GeV/#it{c})"};

    // Histograms
    // Number of events per selection
    eventHist.add("hEventSelection", "hEventSelection", kTH1F, {{5, -0.5f, 4.5f}});
    eventHist.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(1, "All collisions");
    eventHist.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(2, "sel8 cut");
    eventHist.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(3, "posZ cut");
    eventHist.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(4, "INEL>0 cut");
    eventHist.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(5, "With at least a #phi cand");

    // Event information
    eventHist.add("hVertexZ", "hVertexZ", kTH1F, {vertexZAxis});
    eventHist.add("hMultiplicityPercent", "Multiplicity Percentile", kTH1F, {multAxis});

    // Number of MC events per selection for Rec and Gen
    MCeventHist.add("hRecMCEventSelection", "hRecMCEventSelection", kTH1F, {{8, -0.5f, 7.5f}});
    MCeventHist.get<TH1>(HIST("hRecMCEventSelection"))->GetXaxis()->SetBinLabel(1, "All collisions");
    MCeventHist.get<TH1>(HIST("hRecMCEventSelection"))->GetXaxis()->SetBinLabel(2, "kIsTriggerTVX");
    MCeventHist.get<TH1>(HIST("hRecMCEventSelection"))->GetXaxis()->SetBinLabel(3, "kNoTimeFrameBorder");
    MCeventHist.get<TH1>(HIST("hRecMCEventSelection"))->GetXaxis()->SetBinLabel(4, "kNoITSROFrameBorder");
    MCeventHist.get<TH1>(HIST("hRecMCEventSelection"))->GetXaxis()->SetBinLabel(5, "posZ cut");
    MCeventHist.get<TH1>(HIST("hRecMCEventSelection"))->GetXaxis()->SetBinLabel(6, "INEL>0 cut");
    MCeventHist.get<TH1>(HIST("hRecMCEventSelection"))->GetXaxis()->SetBinLabel(7, "With at least a gen coll");
    MCeventHist.get<TH1>(HIST("hRecMCEventSelection"))->GetXaxis()->SetBinLabel(8, "With at least a #phi");

    MCeventHist.add("hGenMCEventSelection", "hGenMCEventSelection", kTH1F, {{5, -0.5f, 4.5f}});
    MCeventHist.get<TH1>(HIST("hGenMCEventSelection"))->GetXaxis()->SetBinLabel(1, "All collisions");
    MCeventHist.get<TH1>(HIST("hGenMCEventSelection"))->GetXaxis()->SetBinLabel(2, "posZ cut");
    MCeventHist.get<TH1>(HIST("hGenMCEventSelection"))->GetXaxis()->SetBinLabel(3, "INEL>0 cut");
    MCeventHist.get<TH1>(HIST("hGenMCEventSelection"))->GetXaxis()->SetBinLabel(4, "With at least a #phi");
    MCeventHist.get<TH1>(HIST("hGenMCEventSelection"))->GetXaxis()->SetBinLabel(5, "With at least a reco coll");

    // MC Event information for Rec and Gen
    MCeventHist.add("hRecMCVertexZ", "hRecMCVertexZ", kTH1F, {vertexZAxis});
    MCeventHist.add("hRecMCMultiplicityPercent", "RecMC Multiplicity Percentile", kTH1F, {multAxis});
    MCeventHist.add("hRecMCGenMultiplicityPercent", "RecMC Gen Multiplicity Percentile", kTH1F, {binnedmultAxis});

    MCeventHist.add("hGenMCVertexZ", "hGenMCVertexZ", kTH1F, {vertexZAxis});
    MCeventHist.add("hGenMCMultiplicityPercent", "GenMC Multiplicity Percentile", kTH1F, {binnedmultAxis});

    // Phi tpological/PID cuts
    PhicandHist.add("hEta", "Eta distribution", kTH1F, {{200, -1.0f, 1.0f}});
    PhicandHist.add("hDcaxy", "Dcaxy distribution", kTH1F, {{200, -1.0f, 1.0f}});
    PhicandHist.add("hDcaz", "Dcaz distribution", kTH1F, {{200, -1.0f, 1.0f}});
    PhicandHist.add("hNsigmaKaonTPC", "NsigmaKaon TPC distribution", kTH2F, {ptK0SAxis, {100, -10.0f, 10.0f}});
    PhicandHist.add("hNsigmaKaonTOF", "NsigmaKaon TOF distribution", kTH2F, {ptK0SAxis, {100, -10.0f, 10.0f}});

    // K0S topological/PID cuts
    K0SHist.add("hDCAV0Daughters", "hDCAV0Daughters", kTH1F, {{55, 0.0f, 2.2f}});
    K0SHist.add("hV0CosPA", "hV0CosPA", kTH1F, {{100, 0.95f, 1.f}});
    K0SHist.add("hNSigmaPosPionFromK0S", "hNSigmaPosPionFromK0Short", kTH2F, {ptK0SAxis, {100, -5.f, 5.f}});
    K0SHist.add("hNSigmaNegPionFromK0S", "hNSigmaNegPionFromK0Short", kTH2F, {ptK0SAxis, {100, -5.f, 5.f}});

    // Phi invariant mass for computing purities and normalisation
    PhipurHist.add("h2PhipurInvMass", "Invariant mass of Phi for Purity (no K0S/Pi)", kTH2F, {binnedmultAxis, PhimassAxis});

    PhipurHist.add("h3PhipurK0SInvMassInc", "Invariant mass of Phi for Purity (K0S) Inclusive", kTH3F, {binnedmultAxis, binnedptK0SAxis, PhimassAxis});
    PhipurHist.add("h3PhipurK0SInvMassFCut", "Invariant mass of Phi for Purity (K0S) Deltay < FirstCut", kTH3F, {binnedmultAxis, binnedptK0SAxis, PhimassAxis});
    PhipurHist.add("h3PhipurK0SInvMassSCut", "Invariant mass of Phi for Purity (K0S) Deltay < SecondCut", kTH3F, {binnedmultAxis, binnedptK0SAxis, PhimassAxis});

    PhipurHist.add("h3PhipurPiInvMassInc", "Invariant mass of Phi for Purity (Pi) Inclusive", kTH3F, {binnedmultAxis, binnedptPiAxis, PhimassAxis});
    PhipurHist.add("h3PhipurPiInvMassFCut", "Invariant mass of Phi for Purity (Pi) Deltay < FirstCut", kTH3F, {binnedmultAxis, binnedptPiAxis, PhimassAxis});
    PhipurHist.add("h3PhipurPiInvMassSCut", "Invariant mass of Phi for Purity (Pi) Deltay < SecondCut", kTH3F, {binnedmultAxis, binnedptPiAxis, PhimassAxis});

    // MCPhi invariant mass for computing purities
    closureMCPhipurHist.add("h2MCPhipurInvMass", "Invariant mass of Phi for Purity (no K0S/Pi)", kTH2F, {binnedmultAxis, PhimassAxis});

    closureMCPhipurHist.add("h3MCPhipurK0SInvMassInc", "Invariant mass of Phi for Purity (K0S) Inclusive", kTH3F, {binnedmultAxis, binnedptK0SAxis, PhimassAxis});
    closureMCPhipurHist.add("h3MCPhipurK0SInvMassFCut", "Invariant mass of Phi for Purity (K0S) Deltay < FirstCut", kTH3F, {binnedmultAxis, binnedptK0SAxis, PhimassAxis});
    closureMCPhipurHist.add("h3MCPhipurK0SInvMassSCut", "Invariant mass of Phi for Purity (K0S) Deltay < SecondCut", kTH3F, {binnedmultAxis, binnedptK0SAxis, PhimassAxis});

    closureMCPhipurHist.add("h3MCPhipurPiInvMassInc", "Invariant mass of Phi for Purity (Pi) Inclusive", kTH3F, {binnedmultAxis, binnedptPiAxis, PhimassAxis});
    closureMCPhipurHist.add("h3MCPhipurPiInvMassFCut", "Invariant mass of Phi for Purity (Pi) Deltay < FirstCut", kTH3F, {binnedmultAxis, binnedptPiAxis, PhimassAxis});
    closureMCPhipurHist.add("h3MCPhipurPiInvMassSCut", "Invariant mass of Phi for Purity (Pi) Deltay < SecondCut", kTH3F, {binnedmultAxis, binnedptPiAxis, PhimassAxis});

    // 2D mass for Phi and K0S for Data
    PhiK0SHist.add("h4PhiK0SSEInc", "2D Invariant mass of Phi and K0Short for Same Event Inclusive", kTHnSparseF, {binnedmultAxis, binnedptK0SAxis, K0SmassAxis, sigPhimassAxis});
    PhiK0SHist.add("h4PhiK0SSEFCut", "2D Invariant mass of Phi and K0Short for Same Event Deltay < FirstCut", kTHnSparseF, {binnedmultAxis, binnedptK0SAxis, K0SmassAxis, sigPhimassAxis});
    PhiK0SHist.add("h4PhiK0SSESCut", "2D Invariant mass of Phi and K0Short for Same Event Deltay < SecondCut", kTHnSparseF, {binnedmultAxis, binnedptK0SAxis, K0SmassAxis, sigPhimassAxis});

    // RecMC K0S coupled to Phi
    MCPhiK0SHist.add("h3RecMCPhiK0SSEInc", "2D Invariant mass of Phi and K0Short for RecMC Inclusive", kTH3F, {binnedmultAxis, binnedptK0SAxis, K0SmassAxis});
    MCPhiK0SHist.add("h3RecMCPhiK0SSEFCut", "2D Invariant mass of Phi and K0Short for RecMC Deltay < FirstCut", kTH3F, {binnedmultAxis, binnedptK0SAxis, K0SmassAxis});
    MCPhiK0SHist.add("h3RecMCPhiK0SSESCut", "2D Invariant mass of Phi and K0Short for RecMC Deltay < SecondCut", kTH3F, {binnedmultAxis, binnedptK0SAxis, K0SmassAxis});

    // GenMC K0S coupled to Phi
    MCPhiK0SHist.add("h2PhiK0SGenMCInc", "K0Short coupled to Phi for GenMC Inclusive", kTH2F, {binnedmultAxis, binnedptK0SAxis});
    MCPhiK0SHist.add("h2PhiK0SGenMCFCut", "K0Short coupled to Phi for GenMC Deltay < FirstCut", kTH2F, {binnedmultAxis, binnedptK0SAxis});
    MCPhiK0SHist.add("h2PhiK0SGenMCSCut", "K0Short coupled to Phi for GenMC Deltay < SecondCut", kTH2F, {binnedmultAxis, binnedptK0SAxis});

    MCPhiK0SHist.add("h2PhiK0SGenMCIncAssocReco", "K0Short coupled to Phi for GenMC Inclusive Associated Reco Collision", kTH2F, {binnedmultAxis, binnedptK0SAxis});
    MCPhiK0SHist.add("h2PhiK0SGenMCFCutAssocReco", "K0Short coupled to Phi for GenMC Deltay < FirstCut Associated Reco Collision", kTH2F, {binnedmultAxis, binnedptK0SAxis});
    MCPhiK0SHist.add("h2PhiK0SGenMCSCutAssocReco", "K0Short coupled to Phi for GenMC Deltay < SecondCut Associated Reco Collision", kTH2F, {binnedmultAxis, binnedptK0SAxis});

    // 2D mass for Phi and K0S for Closure Test
    closureMCPhiK0SHist.add("h4ClosureMCPhiK0SSEInc", "2D Invariant mass of Phi and K0Short for Same Event Inclusive for Closure Test", kTHnSparseF, {binnedmultAxis, binnedptK0SAxis, K0SmassAxis, sigPhimassAxis});
    closureMCPhiK0SHist.add("h4ClosureMCPhiK0SSEFCut", "2D Invariant mass of Phi and K0Short for Same Event Deltay < FirstCut for Closure Test", kTHnSparseF, {binnedmultAxis, binnedptK0SAxis, K0SmassAxis, sigPhimassAxis});
    closureMCPhiK0SHist.add("h4ClosureMCPhiK0SSESCut", "2D Invariant mass of Phi and K0Short for Same Event Deltay < SecondCut for Closure Test", kTHnSparseF, {binnedmultAxis, binnedptK0SAxis, K0SmassAxis, sigPhimassAxis});

    // Phi mass vs Pion NSigma dE/dx for Data
    PhiPionHist.add("h5PhiPiSEInc", "Phi Invariant mass vs Pion nSigma TPC/TOF for Same Event Inclusive", kTHnSparseF, {binnedmultAxis, binnedptPiAxis, {100, -10.0f, 10.0f}, {100, -10.0f, 10.0f}, sigPhimassAxis});
    PhiPionHist.add("h5PhiPiSEFCut", "Phi Invariant mass vs Pion nSigma TPC/TOF for Same Event Deltay < FirstCut", kTHnSparseF, {binnedmultAxis, binnedptPiAxis, {100, -10.0f, 10.0f}, {100, -10.0f, 10.0f}, sigPhimassAxis});
    PhiPionHist.add("h5PhiPiSESCut", "Phi Invariant mass vs Pion nSigma TPC/TOF for Same Event Deltay < SecondCut", kTHnSparseF, {binnedmultAxis, binnedptPiAxis, {100, -10.0f, 10.0f}, {100, -10.0f, 10.0f}, sigPhimassAxis});

    // RecMC Pion coupled to Phi
    MCPhiPionHist.add("h4RecMCPhiPiSEInc", "Phi Invariant mass vs Pion nSigma TPC/TOF for RecMC Inclusive", kTHnSparseF, {binnedmultAxis, binnedptPiAxis, {100, -10.0f, 10.0f}, {100, -10.0f, 10.0f}});
    MCPhiPionHist.add("h4RecMCPhiPiSEFCut", "Phi Invariant mass vs Pion nSigma TPC/TOF for RecMC Deltay < FirstCut", kTHnSparseF, {binnedmultAxis, binnedptPiAxis, {100, -10.0f, 10.0f}, {100, -10.0f, 10.0f}});
    MCPhiPionHist.add("h4RecMCPhiPiSESCut", "Phi Invariant mass vs Pion nSigma TPC/TOF for RecMC Deltay < SecondCut", kTHnSparseF, {binnedmultAxis, binnedptPiAxis, {100, -10.0f, 10.0f}, {100, -10.0f, 10.0f}});

    // GenMC Pion coupled to Phi
    MCPhiPionHist.add("h2PhiPiGenMCInc", "Pion coupled to Phi for GenMC Inclusive", kTH2F, {binnedmultAxis, binnedptPiAxis});
    MCPhiPionHist.add("h2PhiPiGenMCFCut", "Pion coupled to Phi for GenMC Deltay < FirstCut", kTH2F, {binnedmultAxis, binnedptPiAxis});
    MCPhiPionHist.add("h2PhiPiGenMCSCut", "Pion coupled to Phi for GenMC Deltay < SecondCut", kTH2F, {binnedmultAxis, binnedptPiAxis});

    MCPhiPionHist.add("h2PhiPiGenMCIncAssocReco", "Pion coupled to Phi for GenMC Inclusive Associated Reco Collision", kTH2F, {binnedmultAxis, binnedptPiAxis});
    MCPhiPionHist.add("h2PhiPiGenMCFCutAssocReco", "Pion coupled to Phi for GenMC Deltay < FirstCut Associated Reco Collision", kTH2F, {binnedmultAxis, binnedptPiAxis});
    MCPhiPionHist.add("h2PhiPiGenMCSCutAssocReco", "Pion coupled to Phi for GenMC Deltay < SecondCut Associated Reco Collision", kTH2F, {binnedmultAxis, binnedptPiAxis});

    // Phi mass vs Pion NSigma dE/dx for Closure Test
    closureMCPhiPionHist.add("h5ClosureMCPhiPiSEInc", "Phi Invariant mass vs Pion nSigma TPC/TOF for Same Event Inclusive for Closure Test", kTHnSparseF, {binnedmultAxis, binnedptPiAxis, {100, -10.0f, 10.0f}, {100, -10.0f, 10.0f}, sigPhimassAxis});
    closureMCPhiPionHist.add("h5ClosureMCPhiPiSEFCut", "Phi Invariant mass vs Pion nSigma TPC/TOF for Same Event Deltay < FirstCut for Closure Test", kTHnSparseF, {binnedmultAxis, binnedptPiAxis, {100, -10.0f, 10.0f}, {100, -10.0f, 10.0f}, sigPhimassAxis});
    closureMCPhiPionHist.add("h5ClosureMCPhiPiSESCut", "Phi Invariant mass vs Pion nSigma TPC/TOF for Same Event Deltay < SecondCut for Closure Test", kTHnSparseF, {binnedmultAxis, binnedptPiAxis, {100, -10.0f, 10.0f}, {100, -10.0f, 10.0f}, sigPhimassAxis});

    // MCPhi invariant mass for computing efficiencies and MCnormalisation
    PhieffHist.add("h2PhieffInvMass", "Invariant mass of Phi for Efficiency (no K0S/Pi)", kTH2F, {binnedmultAxis, PhimassAxis});

    PhieffHist.add("h3PhieffK0SInvMassInc", "Invariant mass of Phi for Efficiency (K0S) Inclusive", kTH3F, {binnedmultAxis, binnedptK0SAxis, PhimassAxis});
    PhieffHist.add("h3PhieffK0SInvMassFCut", "Invariant mass of Phi for Efficiency (K0S) Deltay < FirstCut", kTH3F, {binnedmultAxis, binnedptK0SAxis, PhimassAxis});
    PhieffHist.add("h3PhieffK0SInvMassSCut", "Invariant mass of Phi for Efficiency (K0S) Deltay < SecondCut", kTH3F, {binnedmultAxis, binnedptK0SAxis, PhimassAxis});

    PhieffHist.add("h3PhieffPiInvMassInc", "Invariant mass of Phi for Efficiency (Pi) Inclusive", kTH3F, {binnedmultAxis, binnedptPiAxis, PhimassAxis});
    PhieffHist.add("h3PhieffPiInvMassFCut", "Invariant mass of Phi for Efficiency (Pi) Deltay < FirstCut", kTH3F, {binnedmultAxis, binnedptPiAxis, PhimassAxis});
    PhieffHist.add("h3PhieffPiInvMassSCut", "Invariant mass of Phi for Efficiency (Pi) Deltay < SecondCut", kTH3F, {binnedmultAxis, binnedptPiAxis, PhimassAxis});

    // GenMC Phi and Phi coupled to K0S and Pion
    PhieffHist.add("h1PhiGenMC", "Phi for GenMC", kTH1F, {binnedmultAxis});
    PhieffHist.add("h1PhiGenMCAssocReco", "Phi for GenMC Associated Reco Collision", kTH1F, {binnedmultAxis});

    PhieffHist.add("h2PhieffK0SGenMCInc", "Phi coupled to K0Short for GenMC Inclusive", kTH2F, {binnedmultAxis, binnedptK0SAxis});
    PhieffHist.add("h2PhieffK0SGenMCFCut", "Phi coupled to K0Short for GenMC Deltay < FirstCut", kTH2F, {binnedmultAxis, binnedptK0SAxis});
    PhieffHist.add("h2PhieffK0SGenMCSCut", "Phi coupled to K0Short for GenMC Deltay < SecondCut", kTH2F, {binnedmultAxis, binnedptK0SAxis});

    PhieffHist.add("h2PhieffK0SGenMCIncAssocReco", "Phi coupled to K0Short for GenMC Inclusive Associated Reco Collision", kTH2F, {binnedmultAxis, binnedptK0SAxis});
    PhieffHist.add("h2PhieffK0SGenMCFCutAssocReco", "Phi coupled to K0Short for GenMC Deltay < FirstCut Associated Reco Collision", kTH2F, {binnedmultAxis, binnedptK0SAxis});
    PhieffHist.add("h2PhieffK0SGenMCSCutAssocReco", "Phi coupled to K0Short for GenMC Deltay < SecondCut Associated Reco Collision", kTH2F, {binnedmultAxis, binnedptK0SAxis});

    PhieffHist.add("h2PhieffPiGenMCInc", "Phi coupled to Pion for GenMC Inclusive", kTH2F, {binnedmultAxis, binnedptPiAxis});
    PhieffHist.add("h2PhieffPiGenMCFCut", "Phi coupled to Pion for GenMC Deltay < FirstCut", kTH2F, {binnedmultAxis, binnedptPiAxis});
    PhieffHist.add("h2PhieffPiGenMCSCut", "Phi coupled to Pion for GenMC Deltay < SecondCut", kTH2F, {binnedmultAxis, binnedptPiAxis});

    PhieffHist.add("h2PhieffPiGenMCIncAssocReco", "Phi coupled to Pion for GenMC Inclusive Associated Reco Collision", kTH2F, {binnedmultAxis, binnedptPiAxis});
    PhieffHist.add("h2PhieffPiGenMCFCutAssocReco", "Phi coupled to Pion for GenMC Deltay < FirstCut Associated Reco Collision", kTH2F, {binnedmultAxis, binnedptPiAxis});
    PhieffHist.add("h2PhieffPiGenMCSCutAssocReco", "Phi coupled to Pion for GenMC Deltay < SecondCut Associated Reco Collision", kTH2F, {binnedmultAxis, binnedptPiAxis});

    // Rapidity smearing matrix for Phi
    PhieffHist.add("h3PhiRapiditySmearing", "Rapidity Smearing Matrix for Phi", kTH3F, {binnedmultAxis, yAxis, yAxis});

    // MCK0S invariant mass and GenMC K0S for computing efficiencies
    K0SeffHist.add("h3K0SeffInvMass", "Invariant mass of K0Short for Efficiency", kTH3F, {binnedmultAxis, binnedptK0SAxis, K0SmassAxis});
    K0SeffHist.add("h2K0SGenMC", "K0Short for GenMC", kTH2F, {binnedmultAxis, binnedptK0SAxis});
    K0SeffHist.add("h2K0SGenMCAssocReco", "K0Short for GenMC Associated Reco Collision", kTH2F, {binnedmultAxis, binnedptK0SAxis});

    // Rapidity smearing matrix for K0S
    K0SeffHist.add("h4K0SRapiditySmearing", "Rapidity Smearing Matrix for K0Short", kTHnSparseF, {binnedmultAxis, binnedptK0SAxis, yAxis, yAxis});

    // MCPion invariant mass and GenMC Pion for computing efficiencies
    PioneffHist.add("h4PieffInvMass", "Invariant mass of Pion for Efficiency", kTHnSparseF, {binnedmultAxis, binnedptPiAxis, {100, -10.0f, 10.0f}, {100, -10.0f, 10.0f}});
    PioneffHist.add("h2PiGenMC", "Pion for GenMC", kTH2F, {binnedmultAxis, binnedptPiAxis});
    PioneffHist.add("h2PiGenMCAssocReco", "Pion for GenMC Associated Reco Collision", kTH2F, {binnedmultAxis, binnedptPiAxis});

    // Rapidity smearing matrix for Pion
    PioneffHist.add("h4PiRapiditySmearing", "Rapidity Smearing Matrix for Pion", kTHnSparseF, {binnedmultAxis, binnedptPiAxis, yAxis, yAxis});
  }

  // Event selection and QA filling
  template <bool isMC, typename T>
  bool acceptEventQA(const T& collision, bool QA)
  {
    if constexpr (!isMC) { // data event
      if (QA)
        eventHist.fill(HIST("hEventSelection"), 0); // all collisions
      if (!collision.sel8())
        return false;
      if (QA)
        eventHist.fill(HIST("hEventSelection"), 1); // sel8 collisions
      if (std::abs(collision.posZ()) >= cutzvertex)
        return false;
      if (QA) {
        eventHist.fill(HIST("hEventSelection"), 2); // vertex-Z selected
        eventHist.fill(HIST("hVertexZ"), collision.posZ());
      }
      if (!collision.isInelGt0())
        return false;
      if (QA)
        eventHist.fill(HIST("hEventSelection"), 3); // INEL>0 collisions
      return true;
    } else { // RecMC event
      if (QA)
        MCeventHist.fill(HIST("hRecMCEventSelection"), 0); // all collisions
      if (!collision.selection_bit(aod::evsel::kIsTriggerTVX))
        return false;
      if (QA)
        MCeventHist.fill(HIST("hRecMCEventSelection"), 1); // kIsTriggerTVX collisions
      if (!collision.selection_bit(aod::evsel::kNoTimeFrameBorder))
        return false;
      if (QA)
        MCeventHist.fill(HIST("hRecMCEventSelection"), 2); // kNoTimeFrameBorder collisions
      if (cfgiskNoITSROFrameBorder && !collision.selection_bit(aod::evsel::kNoITSROFrameBorder))
        return false;
      if (QA)
        MCeventHist.fill(HIST("hRecMCEventSelection"), 3); // kNoITSROFrameBorder collisions (by default not requested by the selection)
      if (std::abs(collision.posZ()) > cutzvertex)
        return false;
      if (QA) {
        MCeventHist.fill(HIST("hRecMCEventSelection"), 4); // vertex-Z selected
        MCeventHist.fill(HIST("hRecMCVertexZ"), collision.posZ());
      }
      if (!collision.isInelGt0())
        return false;
      if (QA)
        MCeventHist.fill(HIST("hRecMCEventSelection"), 5); // INEL>0 collisions
      return true;
    }
  }

  // Single track selection for strangeness sector
  template <typename T>
  bool selectionTrackStrangeness(const T& track)
  {
    if (!track.hasTPC())
      return false;
    if (track.tpcNClsFound() < minTPCnClsFound)
      return false;
    if (track.tpcNClsCrossedRows() < minNCrossedRowsTPC)
      return false;
    if (track.tpcChi2NCl() > maxChi2TPC)
      return false;
    return true;
  }

  // V0 selection
  template <typename T1, typename T2>
  bool selectionV0(const T1& v0, const T2& daughter1, const T2& daughter2)
  {
    if (!selectionTrackStrangeness(daughter1) || !selectionTrackStrangeness(daughter2))
      return false;

    if (v0.v0cosPA() < v0setting_cospa)
      return false;
    if (v0.v0radius() < v0setting_radius)
      return false;

    if (std::abs(daughter1.tpcNSigmaPi()) > NSigmaTPCPion)
      return false;
    if (std::abs(daughter2.tpcNSigmaPi()) > NSigmaTPCPion)
      return false;
    return true;
  }

  // Topological track selection
  template <typename T>
  bool selectionTrackResonance(const T& track)
  {
    if (cfgPrimaryTrack && !track.isPrimaryTrack())
      return false;
    if (cfgGlobalWoDCATrack && !track.isGlobalTrackWoDCA())
      return false;
    if (cfgPVContributor && !track.isPVContributor())
      return false;

    if (track.pt() < cMinKaonPtcut)
      return false;
    if (std::abs(track.dcaZ()) > cMaxDCAzToPVcut)
      return false;
    if (std::abs(track.dcaXY()) > cMaxDCArToPV1 + (cMaxDCArToPV2 / std::pow(track.pt(), cMaxDCArToPV3)))
      return false;
    if (track.tpcNClsFound() < minTPCnClsFound)
      return false;
    return true;
  }

  // PIDKaon track selection
  template <typename T>
  bool selectionPIDKaon(const T& candidate)
  {
    if (!isNoTOF && candidate.hasTOF() && (candidate.tofNSigmaKa() * candidate.tofNSigmaKa() + candidate.tpcNSigmaKa() * candidate.tpcNSigmaKa()) < (nsigmaCutCombinedKa * nsigmaCutCombinedKa))
      return true;
    if (!isNoTOF && !candidate.hasTOF() && std::abs(candidate.tpcNSigmaKa()) < nsigmaCutTPCKa)
      return true;
    if (isNoTOF && std::abs(candidate.tpcNSigmaKa()) < nsigmaCutTPCKa)
      return true;
    return false;
  }

  template <typename T>
  bool selectionPIDKaonpTdependent(const T& candidate)
  {
    if (candidate.pt() < 0.5 && std::abs(candidate.tpcNSigmaKa()) < nsigmaCutTPCKa)
      return true;
    if (candidate.pt() >= 0.5 && candidate.hasTOF() && ((candidate.tofNSigmaKa() * candidate.tofNSigmaKa()) + (candidate.tpcNSigmaKa() * candidate.tpcNSigmaKa())) < (nsigmaCutCombinedKa * nsigmaCutCombinedKa))
      return true;
    return false;
  }

  // Reconstruct the Phi
  template <typename T1, typename T2>
  TLorentzVector recMother(const T1& candidate1, const T2& candidate2, float masscand1, float masscand2)
  {
    TLorentzVector daughter1, daughter2, mother;

    daughter1.SetXYZM(candidate1.px(), candidate1.py(), candidate1.pz(), masscand1); // set the daughter1 4-momentum
    daughter2.SetXYZM(candidate2.px(), candidate2.py(), candidate2.pz(), masscand2); // set the daughter2 4-momentum
    mother = daughter1 + daughter2;                                                  // calculate the mother 4-momentum

    return mother;
  }

  // Topological selection for pions
  template <typename T>
  bool selectionPion(const T& track)
  {
    if (!track.hasITS())
      return false;
    if (track.itsNCls() < minITSnCls)
      return false;
    if (track.itsChi2NCl() > maxChi2ITS)
      return false;

    if (track.pt() < 0.2)
      return false;

    if (track.pt() < 0.8) {
      if (!track.hasTPC())
        return false;
      if (track.tpcNClsFound() < minTPCnClsFound)
        return false;
      if (track.tpcNClsCrossedRows() < minNCrossedRowsTPC)
        return false;
      if (track.tpcChi2NCl() > maxChi2TPC)
        return false;
    }

    if (track.pt() >= 0.5) {
      if (!track.hasTOF())
        return false;
    }

    if (std::abs(track.dcaXY()) > dcaxyMax)
      return false;
    if (std::abs(track.dcaZ()) > dcazMax)
      return false;
    return true;
  }

  // Fill 2D invariant mass histogram for V0 and Phi
  template <bool isMC, typename T>
  void fillInvMass2D(const T& V0, const std::vector<TLorentzVector>& listPhi, float multiplicity, const std::array<float, 3> weights)
  {
    for (const auto& Phi : listPhi) {
      if constexpr (!isMC) { // same event
        PhiK0SHist.fill(HIST("h4PhiK0SSEInc"), multiplicity, V0.pt(), V0.mK0Short(), Phi.M(), weights.at(0));
        if (std::abs(V0.yK0Short() - Phi.Rapidity()) > cfgFirstCutonDeltay)
          continue;
        PhiK0SHist.fill(HIST("h4PhiK0SSEFCut"), multiplicity, V0.pt(), V0.mK0Short(), Phi.M(), weights.at(1));
        if (std::abs(V0.yK0Short() - Phi.Rapidity()) > cfgSecondCutonDeltay)
          continue;
        PhiK0SHist.fill(HIST("h4PhiK0SSESCut"), multiplicity, V0.pt(), V0.mK0Short(), Phi.M(), weights.at(2));
      } else { // MC event
        closureMCPhiK0SHist.fill(HIST("h4ClosureMCPhiK0SSEInc"), multiplicity, V0.pt(), V0.mK0Short(), Phi.M(), weights.at(0));
        if (std::abs(V0.yK0Short() - Phi.Rapidity()) > cfgFirstCutonDeltay)
          continue;
        closureMCPhiK0SHist.fill(HIST("h4ClosureMCPhiK0SSEFCut"), multiplicity, V0.pt(), V0.mK0Short(), Phi.M(), weights.at(1));
        if (std::abs(V0.yK0Short() - Phi.Rapidity()) > cfgSecondCutonDeltay)
          continue;
        closureMCPhiK0SHist.fill(HIST("h4ClosureMCPhiK0SSESCut"), multiplicity, V0.pt(), V0.mK0Short(), Phi.M(), weights.at(2));
      }
    }
  }

  // Fill Phi invariant mass vs Pion nSigmadE/dx histogram
  template <bool isMC, typename T>
  void fillInvMassNSigma(const T& Pi, const std::vector<TLorentzVector>& listPhi, float multiplicity, const std::array<float, 3> weights)
  {
    float nSigmaTPCPi = (Pi.hasTPC() ? Pi.tpcNSigmaPi() : -999);
    float nSigmaTOFPi = (Pi.hasTOF() ? Pi.tofNSigmaPi() : -999);

    for (const auto& Phi : listPhi) {
      if constexpr (!isMC) { // same event
        PhiPionHist.fill(HIST("h5PhiPiSEInc"), multiplicity, Pi.pt(), nSigmaTPCPi, nSigmaTOFPi, Phi.M(), weights.at(0));
        if (std::abs(Pi.rapidity(massPi) - Phi.Rapidity()) > cfgFirstCutonDeltay)
          continue;
        PhiPionHist.fill(HIST("h5PhiPiSEFCut"), multiplicity, Pi.pt(), nSigmaTPCPi, nSigmaTOFPi, Phi.M(), weights.at(1));
        if (std::abs(Pi.rapidity(massPi) - Phi.Rapidity()) > cfgSecondCutonDeltay)
          continue;
        PhiPionHist.fill(HIST("h5PhiPiSESCut"), multiplicity, Pi.pt(), nSigmaTPCPi, nSigmaTOFPi, Phi.M(), weights.at(2));
      } else { // MC event
        closureMCPhiPionHist.fill(HIST("h5ClosureMCPhiPiSEInc"), multiplicity, Pi.pt(), nSigmaTPCPi, nSigmaTOFPi, Phi.M(), weights.at(0));
        if (std::abs(Pi.rapidity(massPi) - Phi.Rapidity()) > cfgFirstCutonDeltay)
          continue;
        closureMCPhiPionHist.fill(HIST("h5ClosureMCPhiPiSEFCut"), multiplicity, Pi.pt(), nSigmaTPCPi, nSigmaTOFPi, Phi.M(), weights.at(1));
        if (std::abs(Pi.rapidity(massPi) - Phi.Rapidity()) > cfgSecondCutonDeltay)
          continue;
        closureMCPhiPionHist.fill(HIST("h5ClosureMCPhiPiSESCut"), multiplicity, Pi.pt(), nSigmaTPCPi, nSigmaTOFPi, Phi.M(), weights.at(2));
      }
    }
  }

  void processQAPurity(SelCollisions::iterator const& collision, FullTracks const& fullTracks, FullV0s const& V0s, V0DauTracks const&)
  {
    // Check if the event selection is passed
    if (!acceptEventQA<false>(collision, true))
      return;

    float multiplicity = collision.centFT0M();
    eventHist.fill(HIST("hMultiplicityPercent"), multiplicity);

    // Defining positive and negative tracks for phi reconstruction
    auto posThisColl = posTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto negThisColl = negTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);

    bool isCountedPhi = false;
    bool isFilledhV0 = false;

    for (const auto& track1 : posThisColl) { // loop over all selected tracks
      if (!selectionTrackResonance(track1) || !selectionPIDKaonpTdependent(track1))
        continue; // topological and PID selection

      PhicandHist.fill(HIST("hEta"), track1.eta());
      PhicandHist.fill(HIST("hDcaxy"), track1.dcaXY());
      PhicandHist.fill(HIST("hDcaz"), track1.dcaZ());
      PhicandHist.fill(HIST("hNsigmaKaonTPC"), track1.tpcInnerParam(), track1.tpcNSigmaKa());
      PhicandHist.fill(HIST("hNsigmaKaonTOF"), track1.tpcInnerParam(), track1.tofNSigmaKa());

      auto track1ID = track1.globalIndex();

      // Loop over all negative candidates
      for (const auto& track2 : negThisColl) {
        if (!selectionTrackResonance(track2) || !selectionPIDKaonpTdependent(track2))
          continue; // topological and PID selection

        auto track2ID = track2.globalIndex();
        if (track2ID == track1ID)
          continue; // condition to avoid double counting of pair

        TLorentzVector recPhi = recMother(track1, track2, massKa, massKa);
        if (std::abs(recPhi.Rapidity()) > cfgyAcceptance)
          continue;

        if (!isCountedPhi) {
          eventHist.fill(HIST("hEventSelection"), 4); // at least a Phi candidate in the event
          isCountedPhi = true;
        }

        PhipurHist.fill(HIST("h2PhipurInvMass"), multiplicity, recPhi.M());

        std::array<bool, 3> isCountedK0S{false, false, false};

        // V0 already reconstructed by the builder
        for (const auto& v0 : V0s) {
          const auto& posDaughterTrack = v0.posTrack_as<V0DauTracks>();
          const auto& negDaughterTrack = v0.negTrack_as<V0DauTracks>();

          // Cut on V0 dynamic columns
          if (!selectionV0(v0, posDaughterTrack, negDaughterTrack))
            continue;

          if (!isFilledhV0) {
            K0SHist.fill(HIST("hDCAV0Daughters"), v0.dcaV0daughters());
            K0SHist.fill(HIST("hV0CosPA"), v0.v0cosPA());

            // Filling the PID of the V0 daughters in the region of the K0 peak
            if (lowmK0S < v0.mK0Short() && v0.mK0Short() < upmK0S) {
              K0SHist.fill(HIST("hNSigmaPosPionFromK0S"), posDaughterTrack.tpcInnerParam(), posDaughterTrack.tpcNSigmaPi());
              K0SHist.fill(HIST("hNSigmaNegPionFromK0S"), negDaughterTrack.tpcInnerParam(), negDaughterTrack.tpcNSigmaPi());
            }
          }

          if (std::abs(v0.yK0Short()) > cfgyAcceptance)
            continue;
          if (!isCountedK0S.at(0)) {
            PhipurHist.fill(HIST("h3PhipurK0SInvMassInc"), multiplicity, v0.pt(), recPhi.M());
            isCountedK0S.at(0) = true;
          }
          if (std::abs(v0.yK0Short() - recPhi.Rapidity()) > cfgFirstCutonDeltay)
            continue;
          if (!isCountedK0S.at(1)) {
            PhipurHist.fill(HIST("h3PhipurK0SInvMassFCut"), multiplicity, v0.pt(), recPhi.M());
            isCountedK0S.at(1) = true;
          }
          if (std::abs(v0.yK0Short() - recPhi.Rapidity()) > cfgSecondCutonDeltay)
            continue;
          if (!isCountedK0S.at(2)) {
            PhipurHist.fill(HIST("h3PhipurK0SInvMassSCut"), multiplicity, v0.pt(), recPhi.M());
            isCountedK0S.at(2) = true;
          }
        }

        isFilledhV0 = true;

        std::array<bool, 3> isCountedPi{false, false, false};

        // Loop over all primary pion candidates
        for (const auto& track : fullTracks) {
          if (!selectionPion(track))
            continue;

          if (std::abs(track.rapidity(massPi)) > cfgyAcceptance)
            continue;
          if (!isCountedPi.at(0)) {
            PhipurHist.fill(HIST("h3PhipurPiInvMassInc"), multiplicity, track.pt(), recPhi.M());
            isCountedPi.at(0) = true;
          }
          if (std::abs(track.rapidity(massPi) - recPhi.Rapidity()) > cfgFirstCutonDeltay)
            continue;
          if (!isCountedPi.at(1)) {
            PhipurHist.fill(HIST("h3PhipurPiInvMassFCut"), multiplicity, track.pt(), recPhi.M());
            isCountedPi.at(1) = true;
          }
          if (std::abs(track.rapidity(massPi) - recPhi.Rapidity()) > cfgSecondCutonDeltay)
            continue;
          if (!isCountedPi.at(2)) {
            PhipurHist.fill(HIST("h3PhipurPiInvMassSCut"), multiplicity, track.pt(), recPhi.M());
            isCountedPi.at(2) = true;
          }
        }
      }
    }
  }

  PROCESS_SWITCH(phik0shortanalysis, processQAPurity, "Process for QA and Phi Purities", true);

  void processSEPhiK0S(soa::Filtered<SelCollisions>::iterator const& collision, FullTracks const&, FullV0s const& V0s, V0DauTracks const&)
  {
    if (!collision.isInelGt0())
      return;

    float multiplicity = collision.centFT0M();
    eventHist.fill(HIST("hMultiplicityPercent"), multiplicity);

    // Defining positive and negative tracks for phi reconstruction
    auto posThisColl = posTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto negThisColl = negTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);

    // V0 already reconstructed by the builder
    for (const auto& v0 : V0s) {
      const auto& posDaughterTrack = v0.posTrack_as<V0DauTracks>();
      const auto& negDaughterTrack = v0.negTrack_as<V0DauTracks>();

      // Cut on V0 dynamic columns
      if (!selectionV0(v0, posDaughterTrack, negDaughterTrack))
        continue;

      if (std::abs(v0.yK0Short()) > cfgyAcceptance)
        continue;

      std::vector<TLorentzVector> listrecPhi;
      std::array<int, 3> counts{};

      // Phi reconstruction
      // Loop over positive candidates
      for (const auto& track1 : posThisColl) { // loop over all selected tracks
        if (!selectionTrackResonance(track1) || !selectionPIDKaonpTdependent(track1))
          continue; // topological and PID selection

        auto track1ID = track1.globalIndex();

        // Loop over all negative candidates
        for (const auto& track2 : negThisColl) {
          if (!selectionTrackResonance(track2) || !selectionPIDKaonpTdependent(track2))
            continue; // topological and PID selection

          auto track2ID = track2.globalIndex();
          if (track2ID == track1ID)
            continue; // condition to avoid double counting of pair

          TLorentzVector recPhi = recMother(track1, track2, massKa, massKa);

          if (recPhi.M() < lowmPhi || recPhi.M() > upmPhi)
            continue;

          if (std::abs(recPhi.Rapidity()) > cfgyAcceptance)
            continue;
          listrecPhi.push_back(recPhi);
          counts.at(0)++;
          if (std::abs(v0.yK0Short() - recPhi.Rapidity()) > cfgFirstCutonDeltay)
            continue;
          counts.at(1)++;
          if (std::abs(v0.yK0Short() - recPhi.Rapidity()) > cfgSecondCutonDeltay)
            continue;
          counts.at(2)++;
        }
      }

      std::array<float, 3> weights{};
      for (unsigned int i = 0; i < counts.size(); i++) {
        weights.at(i) = (counts.at(i) > 0 ? 1. / static_cast<float>(counts.at(i)) : 0);
      }

      fillInvMass2D<false>(v0, listrecPhi, multiplicity, weights);
    }
  }

  PROCESS_SWITCH(phik0shortanalysis, processSEPhiK0S, "Process Same Event for Phi-K0S Analysis", false);

  void processSEPhiPion(soa::Filtered<SelCollisions>::iterator const& collision, FullTracks const& fullTracks)
  {
    if (!collision.isInelGt0())
      return;

    float multiplicity = collision.centFT0M();
    eventHist.fill(HIST("hMultiplicityPercent"), multiplicity);

    // Defining positive and negative tracks for phi reconstruction
    auto posThisColl = posTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto negThisColl = negTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);

    // Loop over all primary pion candidates
    for (const auto& track : fullTracks) {

      // Pion selection
      if (!selectionPion(track))
        continue;

      if (std::abs(track.rapidity(massPi)) > cfgyAcceptance)
        continue;

      std::vector<TLorentzVector> listrecPhi;
      std::array<int, 3> counts{};

      // Phi reconstruction
      // Loop over positive candidates
      for (const auto& track1 : posThisColl) { // loop over all selected tracks
        if (!selectionTrackResonance(track1) || !selectionPIDKaonpTdependent(track1))
          continue; // topological and PID selection

        auto track1ID = track1.globalIndex();

        // Loop over all negative candidates
        for (const auto& track2 : negThisColl) {
          if (!selectionTrackResonance(track2) || !selectionPIDKaonpTdependent(track2))
            continue; // topological and PID selection

          auto track2ID = track2.globalIndex();
          if (track2ID == track1ID)
            continue; // condition to avoid double counting of pair

          TLorentzVector recPhi = recMother(track1, track2, massKa, massKa);

          if (recPhi.M() < lowmPhi || recPhi.M() > upmPhi)
            continue;

          if (std::abs(recPhi.Rapidity()) > cfgyAcceptance)
            continue;
          listrecPhi.push_back(recPhi);
          counts.at(0)++;
          if (std::abs(track.rapidity(massPi) - recPhi.Rapidity()) > cfgFirstCutonDeltay)
            continue;
          counts.at(1)++;
          if (std::abs(track.rapidity(massPi) - recPhi.Rapidity()) > cfgSecondCutonDeltay)
            continue;
          counts.at(2)++;
        }
      }

      std::array<float, 3> weights{};
      for (unsigned int i = 0; i < counts.size(); i++) {
        weights.at(i) = (counts.at(i) > 0 ? 1. / static_cast<float>(counts.at(i)) : 0);
      }

      fillInvMassNSigma<false>(track, listrecPhi, multiplicity, weights);
    }
  }

  PROCESS_SWITCH(phik0shortanalysis, processSEPhiPion, "Process Same Event for Phi-Pion Analysis", false);

  void processRecMCPhiQA(SimCollisions::iterator const& collision, FullMCTracks const& fullMCTracks, FullMCV0s const& V0s, V0DauMCTracks const&, MCCollisions const&, aod::McParticles const&)
  {
    if (!acceptEventQA<true>(collision, true))
      return;

    float multiplicity = collision.centFT0M();
    MCeventHist.fill(HIST("hRecMCMultiplicityPercent"), multiplicity);

    if (!collision.has_mcCollision())
      return;
    MCeventHist.fill(HIST("hRecMCEventSelection"), 6); // with at least a gen collision

    const auto& mcCollision = collision.mcCollision_as<MCCollisions>();
    float genmultiplicity = mcCollision.centFT0M();
    MCeventHist.fill(HIST("hRecMCGenMultiplicityPercent"), genmultiplicity);

    // Defining positive and negative tracks for phi reconstruction
    auto posThisColl = posMCTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto negThisColl = negMCTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);

    bool isCountedPhi = false;

    for (const auto& track1 : posThisColl) { // loop over all selected tracks
      if (!selectionTrackResonance(track1) || !selectionPIDKaonpTdependent(track1))
        continue; // topological and PID selection

      auto track1ID = track1.globalIndex();

      if (!track1.has_mcParticle())
        continue;
      auto MCtrack1 = track1.mcParticle_as<aod::McParticles>();
      if (MCtrack1.pdgCode() != 321 || !MCtrack1.isPhysicalPrimary())
        continue;

      // Loop over all negative candidates
      for (const auto& track2 : negThisColl) {
        if (!selectionTrackResonance(track2) || !selectionPIDKaonpTdependent(track2))
          continue; // topological and PID selection

        auto track2ID = track2.globalIndex();
        if (track2ID == track1ID)
          continue; // condition to avoid double counting of pair

        if (!track2.has_mcParticle())
          continue;
        auto MCtrack2 = track2.mcParticle_as<aod::McParticles>();
        if (MCtrack2.pdgCode() != -321 || !MCtrack2.isPhysicalPrimary())
          continue;

        bool isMCMotherPhi = false;
        auto MCMotherPhi = MCtrack1.mothers_as<aod::McParticles>()[0];
        for (const auto& MotherOfMCtrack1 : MCtrack1.mothers_as<aod::McParticles>()) {
          for (const auto& MotherOfMCtrack2 : MCtrack2.mothers_as<aod::McParticles>()) {
            if (MotherOfMCtrack1 == MotherOfMCtrack2 && MotherOfMCtrack1.pdgCode() == 333) {
              MCMotherPhi = MotherOfMCtrack1;
              isMCMotherPhi = true;
            }
          }
        }

        if (!isMCMotherPhi)
          continue;

        TLorentzVector recPhi = recMother(track1, track2, massKa, massKa);

        PhieffHist.fill(HIST("h3PhiRapiditySmearing"), genmultiplicity, recPhi.Rapidity(), MCMotherPhi.y());

        if (std::abs(recPhi.Rapidity()) > cfgyAcceptance)
          continue;

        if (!isCountedPhi) {
          MCeventHist.fill(HIST("hRecMCEventSelection"), 7); // at least a Phi in the event
          isCountedPhi = true;
        }

        PhieffHist.fill(HIST("h2PhieffInvMass"), genmultiplicity, recPhi.M());

        std::array<bool, 3> isCountedK0S{false, false, false};

        // V0 already reconstructed by the builder
        for (const auto& v0 : V0s) {
          if (!v0.has_mcParticle()) {
            continue;
          }

          auto v0mcparticle = v0.mcParticle();
          if (v0mcparticle.pdgCode() != 310 || !v0mcparticle.isPhysicalPrimary())
            continue;

          const auto& posDaughterTrack = v0.posTrack_as<V0DauMCTracks>();
          const auto& negDaughterTrack = v0.negTrack_as<V0DauMCTracks>();

          // Cut on V0 dynamic columns
          if (!selectionV0(v0, posDaughterTrack, negDaughterTrack))
            continue;

          if (std::abs(v0.yK0Short()) > cfgyAcceptance)
            continue;
          if (!isCountedK0S.at(0)) {
            PhieffHist.fill(HIST("h3PhieffK0SInvMassInc"), genmultiplicity, v0.pt(), recPhi.M());
            isCountedK0S.at(0) = true;
          }
          if (std::abs(v0.yK0Short() - recPhi.Rapidity()) > cfgFirstCutonDeltay)
            continue;
          if (!isCountedK0S.at(1)) {
            PhieffHist.fill(HIST("h3PhieffK0SInvMassFCut"), genmultiplicity, v0.pt(), recPhi.M());
            isCountedK0S.at(1) = true;
          }
          if (std::abs(v0.yK0Short() - recPhi.Rapidity()) > cfgSecondCutonDeltay)
            continue;
          if (!isCountedK0S.at(2)) {
            PhieffHist.fill(HIST("h3PhieffK0SInvMassSCut"), genmultiplicity, v0.pt(), recPhi.M());
            isCountedK0S.at(2) = true;
          }
        }

        std::array<bool, 3> isCountedPi{false, false, false};

        // Loop over all primary pion candidates
        for (const auto& track : fullMCTracks) {
          if (!track.has_mcParticle())
            continue;

          auto MCtrack = track.mcParticle_as<aod::McParticles>();
          if (std::abs(MCtrack.pdgCode()) != 211 || !MCtrack.isPhysicalPrimary())
            continue;

          if (!selectionPion(track))
            continue;

          if (std::abs(track.rapidity(massPi)) > cfgyAcceptance)
            continue;
          if (!isCountedPi.at(0)) {
            PhieffHist.fill(HIST("h3PhieffPiInvMassInc"), genmultiplicity, track.pt(), recPhi.M());
            isCountedPi.at(0) = true;
          }
          if (std::abs(track.rapidity(massPi) - recPhi.Rapidity()) > cfgFirstCutonDeltay)
            continue;
          if (!isCountedPi.at(1)) {
            PhieffHist.fill(HIST("h3PhieffPiInvMassFCut"), genmultiplicity, track.pt(), recPhi.M());
            isCountedPi.at(1) = true;
          }
          if (std::abs(track.rapidity(massPi) - recPhi.Rapidity()) > cfgSecondCutonDeltay)
            continue;
          if (!isCountedPi.at(2)) {
            PhieffHist.fill(HIST("h3PhieffPiInvMassSCut"), genmultiplicity, track.pt(), recPhi.M());
            isCountedPi.at(2) = true;
          }
        }
      }
    }
  }

  PROCESS_SWITCH(phik0shortanalysis, processRecMCPhiQA, "Process for ReCMCQA and Phi in RecMC", false);

  void processRecMCPhiK0S(SimCollisions::iterator const& collision, FullMCTracks const&, FullMCV0s const& V0s, V0DauMCTracks const&, MCCollisions const&, aod::McParticles const& mcParticles)
  {
    if (!acceptEventQA<true>(collision, false))
      return;

    if (!collision.has_mcCollision())
      return;

    const auto& mcCollision = collision.mcCollision_as<MCCollisions>();
    float genmultiplicity = mcCollision.centFT0M();

    // Defining positive and negative tracks for phi reconstruction
    auto posThisColl = posMCTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto negThisColl = negMCTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);

    // Defining McParticles in the collision
    auto mcParticlesThisColl = mcParticles.sliceByCached(aod::mcparticle::mcCollisionId, mcCollision.globalIndex(), cache);

    // V0 already reconstructed by the builder
    for (const auto& v0 : V0s) {
      if (!v0.has_mcParticle())
        continue;

      auto v0mcparticle = v0.mcParticle();
      if (v0mcparticle.pdgCode() != 310 || !v0mcparticle.isPhysicalPrimary())
        continue;

      const auto& posDaughterTrack = v0.posTrack_as<V0DauMCTracks>();
      const auto& negDaughterTrack = v0.negTrack_as<V0DauMCTracks>();

      if (!selectionV0(v0, posDaughterTrack, negDaughterTrack))
        continue;

      K0SeffHist.fill(HIST("h4K0SRapiditySmearing"), genmultiplicity, v0.pt(), v0.yK0Short(), v0mcparticle.y());

      if (std::abs(v0mcparticle.y()) > cfgyAcceptance)
        continue;

      K0SeffHist.fill(HIST("h3K0SeffInvMass"), genmultiplicity, v0.pt(), v0.mK0Short());

      std::array<bool, 3> isCountedMCPhi{false, false, false};

      for (const auto& mcParticle : mcParticlesThisColl) {
        if (mcParticle.pdgCode() != 333)
          continue;
        auto kDaughters = mcParticle.daughters_as<aod::McParticles>();
        if (kDaughters.size() != 2)
          continue;
        bool isPosKaon = false, isNegKaon = false;
        for (const auto& kDaughter : kDaughters) {
          if (kDaughter.pdgCode() == 321)
            isPosKaon = true;
          if (kDaughter.pdgCode() == -321)
            isNegKaon = true;
        }
        if (!isPosKaon || !isNegKaon)
          continue;
        if (std::abs(mcParticle.y()) > cfgyAcceptance)
          continue;

        if (!isCountedMCPhi.at(0)) {
          MCPhiK0SHist.fill(HIST("h3RecMCPhiK0SSEInc"), genmultiplicity, v0.pt(), v0.mK0Short());
          isCountedMCPhi.at(0) = true;
        }
        if (std::abs(v0mcparticle.y() - mcParticle.y()) > cfgFirstCutonDeltay)
          continue;
        if (!isCountedMCPhi.at(1)) {
          MCPhiK0SHist.fill(HIST("h3RecMCPhiK0SSEFCut"), genmultiplicity, v0.pt(), v0.mK0Short());
          isCountedMCPhi.at(1) = true;
        }
        if (std::abs(v0mcparticle.y() - mcParticle.y()) > cfgSecondCutonDeltay)
          continue;
        if (!isCountedMCPhi.at(2)) {
          MCPhiK0SHist.fill(HIST("h3RecMCPhiK0SSESCut"), genmultiplicity, v0.pt(), v0.mK0Short());
          isCountedMCPhi.at(2) = true;
        }
      }
    }
  }

  PROCESS_SWITCH(phik0shortanalysis, processRecMCPhiK0S, "Process RecMC for Phi-K0S Analysis", false);

  void processRecMCPhiPion(SimCollisions::iterator const& collision, FullMCTracks const& fullMCTracks, MCCollisions const&, aod::McParticles const& mcParticles)
  {
    if (!acceptEventQA<true>(collision, false))
      return;

    if (!collision.has_mcCollision())
      return;

    const auto& mcCollision = collision.mcCollision_as<MCCollisions>();
    float genmultiplicity = mcCollision.centFT0M();

    // Defining positive and negative tracks for phi reconstruction
    auto posThisColl = posMCTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto negThisColl = negMCTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);

    // Defining McParticles in the collision
    auto mcParticlesThisColl = mcParticles.sliceByCached(aod::mcparticle::mcCollisionId, mcCollision.globalIndex(), cache);

    // Loop over all primary pion candidates
    for (const auto& track : fullMCTracks) {

      if (!track.has_mcParticle())
        continue;

      auto MCtrack = track.mcParticle_as<aod::McParticles>();
      if (std::abs(MCtrack.pdgCode()) != 211 || !MCtrack.isPhysicalPrimary())
        continue;

      // Pion selection
      if (!selectionPion(track))
        continue;

      PioneffHist.fill(HIST("h4PiRapiditySmearing"), genmultiplicity, track.pt(), track.rapidity(massPi), MCtrack.y());

      if (std::abs(MCtrack.y()) > cfgyAcceptance)
        continue;

      float nsigmaTPC = (track.hasTPC() ? track.tpcNSigmaPi() : -999);
      float nsigmaTOF = (track.hasTOF() ? track.tofNSigmaPi() : -999);

      PioneffHist.fill(HIST("h4PieffInvMass"), genmultiplicity, track.pt(), nsigmaTPC, nsigmaTOF);

      std::array<bool, 3> isCountedMCPhi{false, false, false};

      for (const auto& mcParticle : mcParticlesThisColl) {
        if (mcParticle.pdgCode() != 333)
          continue;
        auto kDaughters = mcParticle.daughters_as<aod::McParticles>();
        if (kDaughters.size() != 2)
          continue;
        bool isPosKaon = false, isNegKaon = false;
        for (const auto& kDaughter : kDaughters) {
          if (kDaughter.pdgCode() == 321)
            isPosKaon = true;
          if (kDaughter.pdgCode() == -321)
            isNegKaon = true;
        }
        if (!isPosKaon || !isNegKaon)
          continue;
        if (std::abs(mcParticle.y()) > cfgyAcceptance)
          continue;

        if (!isCountedMCPhi.at(0)) {
          MCPhiPionHist.fill(HIST("h4RecMCPhiPiSEInc"), genmultiplicity, track.pt(), nsigmaTPC, nsigmaTOF);
          isCountedMCPhi.at(0) = true;
        }
        if (std::abs(MCtrack.y() - mcParticle.y()) > cfgFirstCutonDeltay)
          continue;
        if (!isCountedMCPhi.at(1)) {
          MCPhiPionHist.fill(HIST("h4RecMCPhiPiSEFCut"), genmultiplicity, track.pt(), nsigmaTPC, nsigmaTOF);
          isCountedMCPhi.at(1) = true;
        }
        if (std::abs(MCtrack.y() - mcParticle.y()) > cfgSecondCutonDeltay)
          continue;
        if (!isCountedMCPhi.at(2)) {
          MCPhiPionHist.fill(HIST("h4RecMCPhiPiSESCut"), genmultiplicity, track.pt(), nsigmaTPC, nsigmaTOF);
          isCountedMCPhi.at(2) = true;
        }
      }
    }
  }

  PROCESS_SWITCH(phik0shortanalysis, processRecMCPhiPion, "Process RecMC for Phi-Pion Analysis", false);

  void processRecMCClosurePhiQA(SimCollisions::iterator const& collision, FullMCTracks const& fullMCTracks, FullV0s const& V0s, V0DauMCTracks const&, MCCollisions const&)
  {
    if (!acceptEventQA<true>(collision, true))
      return;

    if (!collision.has_mcCollision())
      return;
    MCeventHist.fill(HIST("hRecMCEventSelection"), 6); // with at least a gen collision

    const auto& mcCollision = collision.mcCollision_as<MCCollisions>();
    float genmultiplicity = mcCollision.centFT0M();
    MCeventHist.fill(HIST("hRecMCGenMultiplicityPercent"), genmultiplicity);

    // Defining positive and negative tracks for phi reconstruction
    auto posThisColl = posMCTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto negThisColl = negMCTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);

    bool isCountedPhi = false;

    for (const auto& track1 : posThisColl) { // loop over all selected tracks
      if (!selectionTrackResonance(track1) || !selectionPIDKaonpTdependent(track1))
        continue; // topological and PID selection

      auto track1ID = track1.globalIndex();

      // Loop over all negative candidates
      for (const auto& track2 : negThisColl) {
        if (!selectionTrackResonance(track2) || !selectionPIDKaonpTdependent(track2))
          continue; // topological and PID selection

        auto track2ID = track2.globalIndex();
        if (track2ID == track1ID)
          continue; // condition to avoid double counting of pair

        TLorentzVector recPhi = recMother(track1, track2, massKa, massKa);
        if (std::abs(recPhi.Rapidity()) > cfgyAcceptance)
          continue;

        if (!isCountedPhi) {
          MCeventHist.fill(HIST("hRecMCEventSelection"), 7); // at least a Phi in the event
          isCountedPhi = true;
        }

        closureMCPhipurHist.fill(HIST("h2MCPhipurInvMass"), genmultiplicity, recPhi.M());

        std::array<bool, 3> isCountedK0S{false, false, false};

        // V0 already reconstructed by the builder
        for (const auto& v0 : V0s) {
          const auto& posDaughterTrack = v0.posTrack_as<V0DauMCTracks>();
          const auto& negDaughterTrack = v0.negTrack_as<V0DauMCTracks>();

          if (!selectionV0(v0, posDaughterTrack, negDaughterTrack))
            continue;

          if (std::abs(v0.yK0Short()) > cfgyAcceptance)
            continue;

          if (!isCountedK0S.at(0)) {
            closureMCPhipurHist.fill(HIST("h3MCPhipurK0SInvMassInc"), genmultiplicity, v0.pt(), recPhi.M());
            isCountedK0S.at(0) = true;
          }
          if (std::abs(v0.yK0Short() - recPhi.Rapidity()) > cfgFirstCutonDeltay)
            continue;
          if (!isCountedK0S.at(1)) {
            closureMCPhipurHist.fill(HIST("h3MCPhipurK0SInvMassFCut"), genmultiplicity, v0.pt(), recPhi.M());
            isCountedK0S.at(1) = true;
          }
          if (std::abs(v0.yK0Short() - recPhi.Rapidity()) > cfgSecondCutonDeltay)
            continue;
          if (!isCountedK0S.at(2)) {
            closureMCPhipurHist.fill(HIST("h3MCPhipurK0SInvMassSCut"), genmultiplicity, v0.pt(), recPhi.M());
            isCountedK0S.at(2) = true;
          }
        }

        std::array<bool, 3> isCountedPi{false, false, false};

        // Loop over all primary pion candidates
        for (const auto& track : fullMCTracks) {

          if (!selectionPion(track))
            continue;

          if (std::abs(track.rapidity(massPi)) > cfgyAcceptance)
            continue;

          if (!isCountedPi.at(0)) {
            closureMCPhipurHist.fill(HIST("h3MCPhipurPiInvMassInc"), genmultiplicity, track.pt(), recPhi.M());
            isCountedPi.at(0) = true;
          }
          if (std::abs(track.rapidity(massPi) - recPhi.Rapidity()) > cfgFirstCutonDeltay)
            continue;
          if (!isCountedPi.at(1)) {
            closureMCPhipurHist.fill(HIST("h3MCPhipurPiInvMassFCut"), genmultiplicity, track.pt(), recPhi.M());
            isCountedPi.at(1) = true;
          }
          if (std::abs(track.rapidity(massPi) - recPhi.Rapidity()) > cfgSecondCutonDeltay)
            continue;
          if (!isCountedPi.at(2)) {
            closureMCPhipurHist.fill(HIST("h3MCPhipurPiInvMassSCut"), genmultiplicity, track.pt(), recPhi.M());
            isCountedPi.at(2) = true;
          }
        }
      }
    }
  }

  PROCESS_SWITCH(phik0shortanalysis, processRecMCClosurePhiQA, "Process for ReCMCQA and Phi in RecMCClosure", false);

  void processRecMCClosurePhiK0S(SimCollisions::iterator const& collision, FullMCTracks const&, FullV0s const& V0s, V0DauMCTracks const&, MCCollisions const&)
  {
    if (!acceptEventQA<true>(collision, false))
      return;

    if (!collision.has_mcCollision())
      return;

    const auto& mcCollision = collision.mcCollision_as<MCCollisions>();
    float genmultiplicity = mcCollision.centFT0M();

    // Defining positive and negative tracks for phi reconstruction
    auto posThisColl = posMCTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto negThisColl = negMCTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);

    // V0 already reconstructed by the builder
    for (const auto& v0 : V0s) {
      const auto& posDaughterTrack = v0.posTrack_as<V0DauMCTracks>();
      const auto& negDaughterTrack = v0.negTrack_as<V0DauMCTracks>();

      if (!selectionV0(v0, posDaughterTrack, negDaughterTrack))
        continue;

      if (std::abs(v0.yK0Short()) > cfgyAcceptance)
        continue;

      std::vector<TLorentzVector> listrecPhi;
      std::array<int, 3> counts{};

      // Phi reconstruction
      for (const auto& track1 : posThisColl) { // loop over all selected tracks
        if (!selectionTrackResonance(track1) || !selectionPIDKaonpTdependent(track1))
          continue; // topological and PID selection

        auto track1ID = track1.globalIndex();

        for (const auto& track2 : negThisColl) {
          if (!selectionTrackResonance(track2) || !selectionPIDKaonpTdependent(track2))
            continue; // topological and PID selection

          auto track2ID = track2.globalIndex();
          if (track2ID == track1ID)
            continue; // condition to avoid double counting of pair

          TLorentzVector recPhi = recMother(track1, track2, massKa, massKa);

          if (recPhi.M() < lowmPhi || recPhi.M() > upmPhi)
            continue;

          if (std::abs(recPhi.Rapidity()) > cfgyAcceptance)
            continue;
          listrecPhi.push_back(recPhi);
          counts.at(0)++;
          if (std::abs(v0.yK0Short() - recPhi.Rapidity()) > cfgFirstCutonDeltay)
            continue;
          counts.at(1)++;
          if (std::abs(v0.yK0Short() - recPhi.Rapidity()) > cfgSecondCutonDeltay)
            continue;
          counts.at(2)++;
        }
      }

      std::array<float, 3> weights{};
      for (unsigned int i = 0; i < counts.size(); i++) {
        weights.at(i) = (counts.at(i) > 0 ? 1. / static_cast<float>(counts.at(i)) : 0);
      }

      fillInvMass2D<true>(v0, listrecPhi, genmultiplicity, weights);
    }
  }

  PROCESS_SWITCH(phik0shortanalysis, processRecMCClosurePhiK0S, "Process RecMC for MCClosure Phi-K0S Analysis", false);

  void processRecMCClosurePhiPion(SimCollisions::iterator const& collision, FullMCTracks const& fullMCTracks, MCCollisions const&)
  {
    if (!acceptEventQA<true>(collision, false))
      return;

    if (!collision.has_mcCollision())
      return;

    const auto& mcCollision = collision.mcCollision_as<MCCollisions>();
    float genmultiplicity = mcCollision.centFT0M();

    // Defining positive and negative tracks for phi reconstruction
    auto posThisColl = posMCTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto negThisColl = negMCTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);

    // Loop over all primary pion candidates
    for (const auto& track : fullMCTracks) {

      // Pion selection
      if (!selectionPion(track))
        continue;

      if (std::abs(track.rapidity(massPi)) > cfgyAcceptance)
        continue;

      std::vector<TLorentzVector> listrecPhi;
      std::array<int, 3> counts{};

      // Phi reconstruction
      for (const auto& track1 : posThisColl) { // loop over all selected tracks
        if (!selectionTrackResonance(track1) || !selectionPIDKaonpTdependent(track1))
          continue; // topological and PID selection

        auto track1ID = track1.globalIndex();

        for (const auto& track2 : negThisColl) {
          if (!selectionTrackResonance(track2) || !selectionPIDKaonpTdependent(track2))
            continue; // topological and PID selection

          auto track2ID = track2.globalIndex();
          if (track2ID == track1ID)
            continue; // condition to avoid double counting of pair

          TLorentzVector recPhi = recMother(track1, track2, massKa, massKa);

          if (recPhi.M() < lowmPhi || recPhi.M() > upmPhi)
            continue;

          if (std::abs(recPhi.Rapidity()) > cfgyAcceptance)
            continue;
          listrecPhi.push_back(recPhi);
          counts.at(0)++;
          if (std::abs(track.rapidity(massPi) - recPhi.Rapidity()) > cfgFirstCutonDeltay)
            continue;
          counts.at(1)++;
          if (std::abs(track.rapidity(massPi) - recPhi.Rapidity()) > cfgSecondCutonDeltay)
            continue;
          counts.at(2)++;
        }
      }

      std::array<float, 3> weights{};
      for (unsigned int i = 0; i < counts.size(); i++) {
        weights.at(i) = (counts.at(i) > 0 ? 1. / static_cast<float>(counts.at(i)) : 0);
      }

      fillInvMassNSigma<true>(track, listrecPhi, genmultiplicity, weights);
    }
  }

  PROCESS_SWITCH(phik0shortanalysis, processRecMCClosurePhiPion, "Process RecMC for MCClosure Phi-Pion Analysis", false);

  void processGenMCPhiQA(MCCollisions::iterator const& mcCollision, soa::SmallGroups<SimCollisions> const& collisions, aod::McParticles const& mcParticles)
  {
    MCeventHist.fill(HIST("hGenMCEventSelection"), 0); // all collisions
    if (std::abs(mcCollision.posZ()) > cutzvertex)
      return;
    MCeventHist.fill(HIST("hGenMCEventSelection"), 1); // vertex-Z selected
    MCeventHist.fill(HIST("hGenMCVertexZ"), mcCollision.posZ());
    if (!pwglf::isINELgtNmc(mcParticles, 0, pdgDB))
      return;
    MCeventHist.fill(HIST("hGenMCEventSelection"), 2); // INEL>0 collisions

    bool isAssocColl = false;
    for (const auto& collision : collisions) {
      if (acceptEventQA<true>(collision, false)) {
        isAssocColl = true;
        break;
      }
    }

    float genmultiplicity = mcCollision.centFT0M();
    MCeventHist.fill(HIST("hGenMCMultiplicityPercent"), genmultiplicity);

    bool isCountedPhi = false;

    for (const auto& mcParticle1 : mcParticles) {
      if (mcParticle1.pdgCode() != 333)
        continue;
      auto kDaughters = mcParticle1.daughters_as<aod::McParticles>();
      if (kDaughters.size() != 2)
        continue;
      bool isPosKaon = false, isNegKaon = false;
      for (const auto& kDaughter : kDaughters) {
        if (kDaughter.pdgCode() == 321)
          isPosKaon = true;
        if (kDaughter.pdgCode() == -321)
          isNegKaon = true;
      }
      if (!isPosKaon || !isNegKaon)
        continue;
      if (std::abs(mcParticle1.y()) > cfgyAcceptance)
        continue;

      if (!isCountedPhi) {
        MCeventHist.fill(HIST("hGenMCEventSelection"), 3); // at least a Phi in the event
        if (isAssocColl)
          MCeventHist.fill(HIST("hGenMCEventSelection"), 4); // with at least a reco collision
        isCountedPhi = true;
      }

      PhieffHist.fill(HIST("h1PhiGenMC"), genmultiplicity);
      if (isAssocColl)
        PhieffHist.fill(HIST("h1PhiGenMCAssocReco"), genmultiplicity);

      std::array<bool, 3> isCountedK0S = {false, false, false};

      for (const auto& mcParticle2 : mcParticles) {
        if (mcParticle2.pdgCode() != 310)
          continue;
        if (!mcParticle2.isPhysicalPrimary())
          continue;
        auto kDaughters2 = mcParticle2.daughters_as<aod::McParticles>();
        if (kDaughters2.size() != 2)
          continue;
        bool isPosPion = false, isNegPion = false;
        for (const auto& kDaughter2 : kDaughters2) {
          if (kDaughter2.pdgCode() == 211)
            isPosPion = true;
          if (kDaughter2.pdgCode() == -211)
            isNegPion = true;
        }
        if (!isPosPion || !isNegPion)
          continue;

        if (std::abs(mcParticle2.y()) > cfgyAcceptance)
          continue;
        if (!isCountedK0S.at(0)) {
          PhieffHist.fill(HIST("h2PhieffK0SGenMCInc"), genmultiplicity, mcParticle2.pt());
          if (isAssocColl)
            PhieffHist.fill(HIST("h2PhieffK0SGenMCFCutAssocReco"), genmultiplicity, mcParticle2.pt());
          isCountedK0S.at(0) = true;
        }
        if (std::abs(mcParticle1.y() - mcParticle2.y()) > cfgFirstCutonDeltay)
          continue;
        if (!isCountedK0S.at(1)) {
          PhieffHist.fill(HIST("h2PhieffK0SGenMCFCut"), genmultiplicity, mcParticle2.pt());
          if (isAssocColl)
            PhieffHist.fill(HIST("h2PhieffK0SGenMCFCutAssocReco"), genmultiplicity, mcParticle2.pt());
          isCountedK0S.at(1) = true;
        }
        if (std::abs(mcParticle1.y() - mcParticle2.y()) > cfgSecondCutonDeltay)
          continue;
        if (!isCountedK0S.at(2)) {
          PhieffHist.fill(HIST("h2PhieffK0SGenMCSCut"), genmultiplicity, mcParticle2.pt());
          if (isAssocColl)
            PhieffHist.fill(HIST("h2PhieffK0SGenMCSCutAssocReco"), genmultiplicity, mcParticle2.pt());
          isCountedK0S.at(2) = true;
        }
      }

      std::array<bool, 3> isCountedPi = {false, false, false};

      for (const auto& mcParticle2 : mcParticles) {
        if (std::abs(mcParticle2.pdgCode()) != 211)
          continue;
        if (!mcParticle2.isPhysicalPrimary())
          continue;

        if (std::abs(mcParticle2.y()) > cfgyAcceptance)
          continue;
        if (!isCountedPi.at(0)) {
          PhieffHist.fill(HIST("h2PhieffPiGenMCInc"), genmultiplicity, mcParticle2.pt());
          if (isAssocColl)
            PhieffHist.fill(HIST("h2PhieffPiGenMCIncAssocReco"), genmultiplicity, mcParticle2.pt());
          isCountedPi.at(0) = true;
        }
        if (std::abs(mcParticle1.y() - mcParticle2.y()) > cfgFirstCutonDeltay)
          continue;
        if (!isCountedPi.at(1)) {
          PhieffHist.fill(HIST("h2PhieffPiGenMCFCut"), genmultiplicity, mcParticle2.pt());
          if (isAssocColl)
            PhieffHist.fill(HIST("h2PhieffPiGenMCFCutAssocReco"), genmultiplicity, mcParticle2.pt());
          isCountedPi.at(1) = true;
        }
        if (std::abs(mcParticle1.y() - mcParticle2.y()) > cfgSecondCutonDeltay)
          continue;
        if (!isCountedPi.at(2)) {
          PhieffHist.fill(HIST("h2PhieffPiGenMCSCut"), genmultiplicity, mcParticle2.pt());
          if (isAssocColl)
            PhieffHist.fill(HIST("h2PhieffPiGenMCSCutAssocReco"), genmultiplicity, mcParticle2.pt());
          isCountedPi.at(2) = true;
        }
      }
    }
  }

  PROCESS_SWITCH(phik0shortanalysis, processGenMCPhiQA, "Process for ReCMCQA and Phi in RecMC", false);

  void processGenMCPhiK0S(MCCollisions::iterator const& mcCollision, soa::SmallGroups<SimCollisions> const& collisions, aod::McParticles const& mcParticles)
  {
    if (std::abs(mcCollision.posZ()) > cutzvertex)
      return;
    if (!pwglf::isINELgtNmc(mcParticles, 0, pdgDB))
      return;

    bool isAssocColl = false;
    for (const auto& collision : collisions) {
      if (acceptEventQA<true>(collision, false)) {
        isAssocColl = true;
        break;
      }
    }

    float genmultiplicity = mcCollision.centFT0M();
    MCeventHist.fill(HIST("hGenMCMultiplicityPercent"), genmultiplicity);

    for (const auto& mcParticle1 : mcParticles) {
      if (mcParticle1.pdgCode() != 310)
        continue;
      if (!mcParticle1.isPhysicalPrimary())
        continue;
      auto kDaughters1 = mcParticle1.daughters_as<aod::McParticles>();
      if (kDaughters1.size() != 2)
        continue;
      bool isPosPion = false, isNegPion = false;
      for (const auto& kDaughter1 : kDaughters1) {
        if (kDaughter1.pdgCode() == 211)
          isPosPion = true;
        if (kDaughter1.pdgCode() == -211)
          isNegPion = true;
      }
      if (!isPosPion || !isNegPion)
        continue;
      if (std::abs(mcParticle1.y()) > cfgyAcceptance)
        continue;

      K0SeffHist.fill(HIST("h2K0SGenMC"), genmultiplicity, mcParticle1.pt());
      if (isAssocColl)
        K0SeffHist.fill(HIST("h2K0SGenMCAssocReco"), genmultiplicity, mcParticle1.pt());

      std::array<bool, 3> isCountedPhi = {false, false, false};

      for (const auto& mcParticle2 : mcParticles) {
        if (mcParticle2.pdgCode() != 333)
          continue;
        auto kDaughters2 = mcParticle2.daughters_as<aod::McParticles>();
        if (kDaughters2.size() != 2)
          continue;
        bool isPosKaon = false, isNegKaon = false;
        for (const auto& kDaughter2 : kDaughters2) {
          if (kDaughter2.pdgCode() == 321)
            isPosKaon = true;
          if (kDaughter2.pdgCode() == -321)
            isNegKaon = true;
        }
        if (!isPosKaon || !isNegKaon)
          continue;

        if (std::abs(mcParticle2.y()) > cfgyAcceptance)
          continue;
        if (!isCountedPhi.at(0)) {
          MCPhiK0SHist.fill(HIST("h2PhiK0SGenMCInc"), genmultiplicity, mcParticle1.pt());
          if (isAssocColl)
            MCPhiK0SHist.fill(HIST("h2PhiK0SGenMCIncAssocReco"), genmultiplicity, mcParticle1.pt());
          isCountedPhi.at(0) = true;
        }
        if (std::abs(mcParticle1.y() - mcParticle2.y()) > cfgFirstCutonDeltay)
          continue;
        if (!isCountedPhi.at(1)) {
          MCPhiK0SHist.fill(HIST("h2PhiK0SGenMCFCut"), genmultiplicity, mcParticle1.pt());
          if (isAssocColl)
            MCPhiK0SHist.fill(HIST("h2PhiK0SGenMCFCutAssocReco"), genmultiplicity, mcParticle1.pt());
          isCountedPhi.at(1) = true;
        }
        if (std::abs(mcParticle1.y() - mcParticle2.y()) > cfgSecondCutonDeltay)
          continue;
        if (!isCountedPhi.at(2)) {
          MCPhiK0SHist.fill(HIST("h2PhiK0SGenMCSCut"), genmultiplicity, mcParticle1.pt());
          if (isAssocColl)
            MCPhiK0SHist.fill(HIST("h2PhiK0SGenMCSCutAssocReco"), genmultiplicity, mcParticle1.pt());
          isCountedPhi.at(2) = true;
        }
      }
    }
  }

  PROCESS_SWITCH(phik0shortanalysis, processGenMCPhiK0S, "Process GenMC for Phi-K0S Analysis", false);

  void processGenMCPhiPion(MCCollisions::iterator const& mcCollision, soa::SmallGroups<SimCollisions> const& collisions, aod::McParticles const& mcParticles)
  {
    if (std::abs(mcCollision.posZ()) > cutzvertex)
      return;
    if (!pwglf::isINELgtNmc(mcParticles, 0, pdgDB))
      return;

    bool isAssocColl = false;
    for (const auto& collision : collisions) {
      if (acceptEventQA<true>(collision, false)) {
        isAssocColl = true;
        break;
      }
    }

    float genmultiplicity = mcCollision.centFT0M();
    MCeventHist.fill(HIST("hGenMCMultiplicityPercent"), genmultiplicity);

    for (const auto& mcParticle1 : mcParticles) {
      if (std::abs(mcParticle1.pdgCode()) != 211)
        continue;
      if (!mcParticle1.isPhysicalPrimary())
        continue;
      if (std::abs(mcParticle1.y()) > cfgyAcceptance)
        continue;

      PioneffHist.fill(HIST("h2PiGenMC"), genmultiplicity, mcParticle1.pt());
      if (isAssocColl)
        PioneffHist.fill(HIST("h2PiGenMCAssocReco"), genmultiplicity, mcParticle1.pt());

      std::array<bool, 3> isCountedPhi = {false, false, false};

      for (const auto& mcParticle2 : mcParticles) {
        if (mcParticle2.pdgCode() != 333)
          continue;
        auto kDaughters = mcParticle2.daughters_as<aod::McParticles>();
        if (kDaughters.size() != 2)
          continue;
        bool isPosKaon = false, isNegKaon = false;
        for (const auto& kDaughter : kDaughters) {
          if (kDaughter.pdgCode() == 321)
            isPosKaon = true;
          if (kDaughter.pdgCode() == -321)
            isNegKaon = true;
        }
        if (!isPosKaon || !isNegKaon)
          continue;

        if (std::abs(mcParticle2.y()) > cfgyAcceptance)
          continue;
        if (!isCountedPhi.at(0)) {
          MCPhiPionHist.fill(HIST("h2PhiPiGenMCInc"), genmultiplicity, mcParticle1.pt());
          if (isAssocColl)
            MCPhiPionHist.fill(HIST("h2PhiPiGenMCIncAssocReco"), genmultiplicity, mcParticle1.pt());
          isCountedPhi.at(0) = true;
        }
        if (std::abs(mcParticle1.y() - mcParticle2.y()) > cfgFirstCutonDeltay)
          continue;
        if (!isCountedPhi.at(1)) {
          MCPhiPionHist.fill(HIST("h2PhiPiGenMCFCut"), genmultiplicity, mcParticle1.pt());
          if (isAssocColl)
            MCPhiPionHist.fill(HIST("h2PhiPiGenMCFCutAssocReco"), genmultiplicity, mcParticle1.pt());
          isCountedPhi.at(1) = true;
        }
        if (std::abs(mcParticle1.y() - mcParticle2.y()) > cfgSecondCutonDeltay)
          continue;
        if (!isCountedPhi.at(2)) {
          MCPhiPionHist.fill(HIST("h2PhiPiGenMCSCut"), genmultiplicity, mcParticle1.pt());
          if (isAssocColl)
            MCPhiPionHist.fill(HIST("h2PhiPiGenMCSCutAssocReco"), genmultiplicity, mcParticle1.pt());
          isCountedPhi.at(2) = true;
        }
      }
    }
  }

  PROCESS_SWITCH(phik0shortanalysis, processGenMCPhiPion, "Process GenMC for Phi-Pion Analysis", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<phik0shortanalysis>(cfgc, TaskName{"lf-phik0shortanalysis"})};
}
