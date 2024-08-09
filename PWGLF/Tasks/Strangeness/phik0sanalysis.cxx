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
#include <TDatabasePDG.h>
#include <Math/Vector4D.h>
#include <vector>
#include <cmath>
#include <array>
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

namespace
{
static constexpr int nMultBin = 10;
static constexpr int nPtBin = 3;

constexpr float flowmPhiInc[nMultBin] = {1.01074f, 1.01073f, 1.01072f, 1.01074f, 1.01075f, 1.01074f, 1.01075f, 1.01074f, 1.01073f, 1.01074f};
constexpr float fupmPhiInc[nMultBin] = {1.02778f, 1.02777f, 1.02776f, 1.02778f, 1.02779f, 1.02778f, 1.02779f, 1.02778f, 1.02777f, 1.02778f};

constexpr float flowmPhiFCut[nMultBin] = {1.01072f, 1.01073f, 1.01072f, 1.01074f, 1.01075f, 1.01076f, 1.01076f, 1.01076f, 1.01075f, 1.01073f};
constexpr float fupmPhiFCut[nMultBin] = {1.02776f, 1.02777f, 1.02776f, 1.02778f, 1.02779f, 1.02778f, 1.02778f, 1.02778f, 1.02779f, 1.02777f};

constexpr float flowmPhiSCut[nMultBin] = {1.01072f, 1.01074f, 1.01070f, 1.01076f, 1.01075f, 1.01077f, 1.01075f, 1.01075f, 1.01076f, 1.01077f};
constexpr float fupmPhiSCut[nMultBin] = {1.02776f, 1.02778f, 1.02774f, 1.02780f, 1.02779f, 1.02781f, 1.02779f, 1.02779f, 1.02780f, 1.02774f};

static constexpr float multBin[nMultBin + 1] = {0.0, 1.0, 5.0, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 70.0, 100.0};
static constexpr float pTBin[nPtBin + 1] = {0.0, 0.5, 1.2, 10.0};

static constexpr std::string_view PhiK0SSEInc[nMultBin] = {"h2PhiK0SSEInc_0_1", "h2PhiK0SSEInc_1_5", "h2PhiK0SSEInc_5_10", "h2PhiK0SSEInc_10_15", "h2PhiK0SSEInc_15_20",
                                                           "h2PhiK0SSEInc_20_30", "h2PhiK0SSEInc_30_40", "h2PhiK0SSEInc_40_50", "h2PhiK0SSEInc_50_70", "h2PhiK0SSEInc_70_100"};
static constexpr std::string_view PhiK0SSEFCut[nMultBin] = {"h2PhiK0SSEFCut_0_1", "h2PhiK0SSEFCut_1_5", "h2PhiK0SSEFCut_5_10", "h2PhiK0SSEFCut_10_15", "h2PhiK0SSEFCut_15_20",
                                                            "h2PhiK0SSEFCut_20_30", "h2PhiK0SSEFCut_30_40", "h2PhiK0SSEFCut_40_50", "h2PhiK0SSEFCut_50_70", "h2PhiK0SSEFCut_70_100"};
static constexpr std::string_view PhiK0SSESCut[nMultBin] = {"h2PhiK0SSESCut_0_1", "h2PhiK0SSESCut_1_5", "h2PhiK0SSESCut_5_10", "h2PhiK0SSESCut_10_15", "h2PhiK0SSESCut_15_20",
                                                            "h2PhiK0SSESCut_20_30", "h2PhiK0SSESCut_30_40", "h2PhiK0SSESCut_40_50", "h2PhiK0SSESCut_50_70", "h2PhiK0SSESCut_70_100"};

static constexpr std::string_view PhiPiSEInc[nMultBin] = {"h2PhiPiSEInc_0_1", "h2PhiPiSEInc_1_5", "h2PhiPiSEInc_5_10", "h2PhiPiSEInc_10_15", "h2PhiPiSEInc_15_20",
                                                          "h2PhiPiSEInc_20_30", "h2PhiPiSEInc_30_40", "h2PhiPiSEInc_40_50", "h2PhiPiSEInc_50_70", "h2PhiPiSEInc_70_100"};
static constexpr std::string_view PhiPiSEFCut[nMultBin] = {"h2PhiPiSEFCut_0_1", "h2PhiPiSEFCut_1_5", "h2PhiPiSEFCut_5_10", "h2PhiPiSEFCut_10_15", "h2PhiPiSEFCut_15_20",
                                                           "h2PhiPiSEFCut_20_30", "h2PhiPiSEFCut_30_40", "h2PhiPiSEFCut_40_50", "h2PhiPiSEFCut_50_70", "h2PhiPiSEFCut_70_100"};
static constexpr std::string_view PhiPiSESCut[nMultBin] = {"h2PhiPiSESCut_0_1", "h2PhiPiSESCut_1_5", "h2PhiPiSESCut_5_10", "h2PhiPiSESCut_10_15", "h2PhiPiSESCut_15_20",
                                                           "h2PhiPiSESCut_20_30", "h2PhiPiSESCut_30_40", "h2PhiPiSESCut_40_50", "h2PhiPiSESCut_50_70", "h2PhiPiSESCut_70_100"};

static constexpr std::string_view MCPhiK0SSEInc[nMultBin] = {"h2RecMCPhiK0SSEInc_0_1", "h2RecMCPhiK0SSEInc_1_5", "h2RecMCPhiK0SSEInc_5_10", "h2RecMCPhiK0SSEInc_10_15", "h2RecMCPhiK0SSEInc_15_20",
                                                             "h2RecMCPhiK0SSEInc_20_30", "h2RecMCPhiK0SSEInc_30_40", "h2RecMCPhiK0SSEInc_40_50", "h2RecMCPhiK0SSEInc_50_70", "h2RecMCPhiK0SSEInc_70_100"};
static constexpr std::string_view MCPhiK0SSEFCut[nMultBin] = {"h2RecMCPhiK0SSEFCut_0_1", "h2RecMCPhiK0SSEFCut_1_5", "h2RecMCPhiK0SSEFCut_5_10", "h2RecMCPhiK0SSEFCut_10_15", "h2RecMCPhiK0SSEFCut_15_20",
                                                              "h2RecMCPhiK0SSEFCut_20_30", "h2RecMCPhiK0SSEFCut_30_40", "h2RecMCPhiK0SSEFCut_40_50", "h2RecMCPhiK0SSEFCut_50_70", "h2RecMCPhiK0SSEFCut_70_100"};
static constexpr std::string_view MCPhiK0SSESCut[nMultBin] = {"h2RecMCPhiK0SSESCut_0_1", "h2RecMCPhiK0SSESCut_1_5", "h2RecMCPhiK0SSESCut_5_10", "h2RecMCPhiK0SSESCut_10_15", "h2RecMCPhiK0SSESCut_15_20",
                                                              "h2RecMCPhiK0SSESCut_20_30", "h2RecMCPhiK0SSESCut_30_40", "h2RecMCPhiK0SSESCut_40_50", "h2RecMCPhiK0SSESCut_50_70", "h2RecMCPhiK0SSESCut_70_100"};

static constexpr std::string_view MCPhiPiSEInc[nMultBin] = {"h2RecMCPhiPiSEInc_0_1", "h2RecMCPhiPiSEInc_1_5", "h2RecMCPhiPiSEInc_5_10", "h2RecMCPhiPiSEInc_10_15", "h2RecMCPhiPiSEInc_15_20",
                                                            "h2RecMCPhiPiSEInc_20_30", "h2RecMCPhiPiSEInc_30_40", "h2RecMCPhiPiSEInc_40_50", "h2RecMCPhiPiSEInc_50_70", "h2RecMCPhiPiSEInc_70_100"};
static constexpr std::string_view MCPhiPiSEFCut[nMultBin] = {"h2RecMCPhiPiSEFCut_0_1", "h2RecMCPhiPiSEFCut_1_5", "h2RecMCPhiPiSEFCut_5_10", "h2RecMCPhiPiSEFCut_10_15", "h2RecMCPhiPiSEFCut_15_20",
                                                             "h2RecMCPhiPiSEFCut_20_30", "h2RecMCPhiPiSEFCut_30_40", "h2RecMCPhiPiSEFCut_40_50", "h2RecMCPhiPiSEFCut_50_70", "h2RecMCPhiPiSEFCut_70_100"};
static constexpr std::string_view MCPhiPiSESCut[nMultBin] = {"h2RecMCPhiPiSESCut_0_1", "h2RecMCPhiPiSESCut_1_5", "h2RecMCPhiPiSESCut_5_10", "h2RecMCPhiPiSESCut_10_15", "h2RecMCPhiPiSESCut_15_20",
                                                             "h2RecMCPhiPiSESCut_20_30", "h2RecMCPhiPiSESCut_30_40", "h2RecMCPhiPiSESCut_40_50", "h2RecMCPhiPiSESCut_50_70", "h2RecMCPhiPiSESCut_70_100"};
} // namespace

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

  // Configurable for event selection
  Configurable<float> cutzvertex{"cutzvertex", 10.0f, "Accepted z-vertex range (cm)"};

  // Configurables for V0 selection
  Configurable<float> minTPCnClsFound{"minTPCnClsFound", 80.0f, "min number of found TPC clusters"};
  Configurable<float> minNCrossedRowsTPC{"minNCrossedRowsTPC", 80.0f, "min number of TPC crossed rows"};
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

  // Configurables on Phi mass
  Configurable<int> nBins{"nBins", 15, "N bins in cfgPhimassaxis"};
  Configurable<std::vector<float>> lowmPhiInc{"lowmPhiInc", std::vector<float>{flowmPhiInc, flowmPhiInc + nMultBin}, "Lower limits on Phi mass Inclusive"};
  Configurable<std::vector<float>> upmPhiInc{"upmPhiInc", std::vector<float>{fupmPhiInc, fupmPhiInc + nMultBin}, "Upper limits on Phi mass Inclusive"};
  Configurable<std::vector<float>> lowmPhiFCut{"lowmPhiFCut", std::vector<float>{flowmPhiFCut, flowmPhiFCut + nMultBin}, "Lower limits on Phi mass First Cut"};
  Configurable<std::vector<float>> upmPhiFCut{"upmPhiFCut", std::vector<float>{fupmPhiFCut, fupmPhiFCut + nMultBin}, "Upper limits on Phi mass First Cut"};
  Configurable<std::vector<float>> lowmPhiSCut{"lowmPhiSCut", std::vector<float>{flowmPhiSCut, flowmPhiSCut + nMultBin}, "Lower limits on Phi mass Second Cut"};
  Configurable<std::vector<float>> upmPhiSCut{"upmPhiSCut", std::vector<float>{fupmPhiSCut, fupmPhiSCut + nMultBin}, "Upper limits on Phi mass Second Cut"};

  // Configurables for phi selection
  Configurable<double> cMinPtcut{"cMinPtcut", 0.15, "Track minimum pt cut"};
  Configurable<float> cfgCutCharge{"cfgCutCharge", 0.0, "Cut on charge"};
  Configurable<bool> cfgPrimaryTrack{"cfgPrimaryTrack", false, "Primary track selection"};
  Configurable<bool> cfgGlobalWoDCATrack{"cfgGlobalWoDCATrack", true, "Global track selection without DCA"};
  Configurable<bool> cfgPVContributor{"cfgPVContributor", true, "PV contributor track selection"};
  Configurable<double> cMaxDCArToPVcut{"cMaxDCArToPVcut", 0.5, "Track DCAr cut to PV Maximum"};
  Configurable<double> cMaxDCAzToPVcut{"cMaxDCAzToPVcut", 2.0, "Track DCAz cut to PV Maximum"};

  Configurable<bool> isNoTOF{"isNoTOF", false, "isNoTOF"};
  Configurable<float> nsigmaCutTPCKa{"nsigmacutTPC", 3.0, "Value of the TPC Nsigma cut"};
  Configurable<float> nsigmaCutCombinedKa{"nsigmaCutCombined", 3.0, "Value of the TOF Nsigma cut"};

  // Configurables for pions(extra with respect to a few of those defined for V0)
  Configurable<float> minITSnCls{"minITSnCls", 4.0f, "min number of ITS clusters"};
  Configurable<float> maxChi2ITS{"maxChi2ITS", 36.0f, "max chi2 per cluster ITS"};
  Configurable<float> dcaxyMax{"dcaxyMax", 0.1f, "Maximum DCAxy to primary vertex"};
  Configurable<float> dcazMax{"dcazMax", 0.1f, "Maximum DCAz to primary vertex"};
  Configurable<float> NSigmaTOFPion{"NSigmaTOFPion", 5.0, "NSigmaTOFPion"};
  Configurable<std::vector<float>> binspT{"binspT", std::vector<float>{0.0, 0.5, 1.2, 10.0}, "pT bin limits for pions"};

  // Configurables for invariant mass histograms filling
  Configurable<double> cfgFirstCutonDeltay{"cgfFirstCutonDeltay", 0.5, "First upper bound on Deltay selection"};
  Configurable<double> cfgSecondCutonDeltay{"cgfSecondCutonDeltay", 0.2, "Second upper bound on Deltay selection"};

  // Configurable for event mixing
  Configurable<int> cfgNoMixedEvents{"cfgNoMixedEvents", 5, "Number of mixed events per event"};

  // Configurable for RecMC
  Configurable<bool> cfgiskNoITSROFrameBorder{"cfgiskNoITSROFrameBorder", false, "kNoITSROFrameBorder request on RecMC collisions"};

  // Configurable axis
  ConfigurableAxis axisVertex{"axisVertex", {20, -10, 10}, "vertex axis for bin"};
  ConfigurableAxis axisMultiplicityClass{"axisMultiplicityClass", {20, 0, 100}, "multiplicity percentile for bin"};
  ConfigurableAxis axisMultiplicity{"axisMultiplicity", {2000, 0, 10000}, "TPC multiplicity  for bin"};

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

  void init(InitContext const&)
  {
    // Axes
    AxisSpec K0SmassAxis = {200, 0.45f, 0.55f, "#it{M}_{inv} [GeV/#it{c}^{2}]"};
    AxisSpec PhimassAxis = {200, 0.9f, 1.2f, "#it{M}_{inv} [GeV/#it{c}^{2}]"};
    AxisSpec vertexZAxis = {100, -15.f, 15.f, "vrtx_{Z} [cm]"};
    AxisSpec deltayAxis = {16, 0.0f, 0.8f, "|#it{#Deltay}|"};
    AxisSpec multAxis = {120, 0.0f, 120.0f, "centFT0M"};
    AxisSpec binnedmultAxis{{0.0, 1.0, 5.0, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 70.0, 100.0}, "centFT0M"};
    AxisSpec ptAxis = {100, 0.0f, 10.0f, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec binnedptAxis{{0.0, 0.5, 1.2, 10.0}, "#it{p}_{T} (GeV/#it{c})"};
    std::vector<AxisSpec> cfgPhimassAxisInc;
    std::vector<AxisSpec> cfgPhimassAxisFCut;
    std::vector<AxisSpec> cfgPhimassAxisSCut;
    for (int i = 0; i < nMultBin; i++) {
      cfgPhimassAxisInc.push_back({nBins, lowmPhiInc->at(i), upmPhiInc->at(i), "#it{M}_{inv} [GeV/#it{c}^{2}]"});
      cfgPhimassAxisFCut.push_back({nBins, lowmPhiFCut->at(i), upmPhiFCut->at(i), "#it{M}_{inv} [GeV/#it{c}^{2}]"});
      cfgPhimassAxisSCut.push_back({nBins, lowmPhiSCut->at(i), upmPhiSCut->at(i), "#it{M}_{inv} [GeV/#it{c}^{2}]"});
    }

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

    // Histo to check phi numbers in data
    eventHist.add("thereisnoPhiwK0S", "thereisnoPhiwK0S", kTH1F, {{3, -0.5, 2.5}});
    eventHist.get<TH1>(HIST("thereisnoPhiwK0S"))->GetXaxis()->SetBinLabel(1, "Inclusive");
    eventHist.get<TH1>(HIST("thereisnoPhiwK0S"))->GetXaxis()->SetBinLabel(2, "|#Delta y|<0.5");
    eventHist.get<TH1>(HIST("thereisnoPhiwK0S"))->GetXaxis()->SetBinLabel(3, "|#Delta y|<0.2");

    eventHist.add("thereisnoPhiwPi", "thereisnoPhiwPi", kTH1F, {{3, -0.5, 2.5}});
    eventHist.get<TH1>(HIST("thereisnoPhiwPi"))->GetXaxis()->SetBinLabel(1, "Inclusive");
    eventHist.get<TH1>(HIST("thereisnoPhiwPi"))->GetXaxis()->SetBinLabel(2, "|#Delta y|<0.5");
    eventHist.get<TH1>(HIST("thereisnoPhiwPi"))->GetXaxis()->SetBinLabel(3, "|#Delta y|<0.2");

    // Number of MC events per selection for Rec and Gen
    MCeventHist.add("hRecMCEventSelection", "hRecMCEventSelection", kTH1F, {{7, -0.5f, 6.5f}});
    MCeventHist.get<TH1>(HIST("hRecMCEventSelection"))->GetXaxis()->SetBinLabel(1, "All collisions");
    MCeventHist.get<TH1>(HIST("hRecMCEventSelection"))->GetXaxis()->SetBinLabel(2, "kIsTriggerTVX");
    MCeventHist.get<TH1>(HIST("hRecMCEventSelection"))->GetXaxis()->SetBinLabel(3, "kNoTimeFrameBorder");
    MCeventHist.get<TH1>(HIST("hRecMCEventSelection"))->GetXaxis()->SetBinLabel(4, "kNoITSROFrameBorder");
    MCeventHist.get<TH1>(HIST("hRecMCEventSelection"))->GetXaxis()->SetBinLabel(5, "posZ cut");
    MCeventHist.get<TH1>(HIST("hRecMCEventSelection"))->GetXaxis()->SetBinLabel(6, "INEL>0 cut");
    MCeventHist.get<TH1>(HIST("hRecMCEventSelection"))->GetXaxis()->SetBinLabel(7, "With at least a #phi");

    // MC Event information for Rec and Gen
    MCeventHist.add("hRecMCVertexZ", "hRecMCVertexZ", kTH1F, {vertexZAxis});
    MCeventHist.add("hRecMCMultiplicityPercent", "RecMC Multiplicity Percentile", kTH1F, {multAxis});

    // Histo to check phi numbers in RecMC
    MCeventHist.add("thereisnoPhiwK0SMC", "thereisnoPhiwK0SMC", kTH1F, {{3, -0.5, 2.5}});
    MCeventHist.get<TH1>(HIST("thereisnoPhiwK0SMC"))->GetXaxis()->SetBinLabel(1, "Inclusive");
    MCeventHist.get<TH1>(HIST("thereisnoPhiwK0SMC"))->GetXaxis()->SetBinLabel(2, "|#Delta y|<0.5");
    MCeventHist.get<TH1>(HIST("thereisnoPhiwK0SMC"))->GetXaxis()->SetBinLabel(3, "|#Delta y|<0.2");

    MCeventHist.add("thereisnoPhiwPiMC", "thereisnoPhiwPiMC", kTH1F, {{3, -0.5, 2.5}});
    MCeventHist.get<TH1>(HIST("thereisnoPhiwPiMC"))->GetXaxis()->SetBinLabel(1, "Inclusive");
    MCeventHist.get<TH1>(HIST("thereisnoPhiwPiMC"))->GetXaxis()->SetBinLabel(2, "|#Delta y|<0.5");
    MCeventHist.get<TH1>(HIST("thereisnoPhiwPiMC"))->GetXaxis()->SetBinLabel(3, "|#Delta y|<0.2");

    MCeventHist.add("hGenMCEventSelection", "hGenMCEventSelection", kTH1F, {{5, -0.5f, 4.5f}});
    MCeventHist.get<TH1>(HIST("hGenMCEventSelection"))->GetXaxis()->SetBinLabel(1, "All collisions");
    MCeventHist.get<TH1>(HIST("hGenMCEventSelection"))->GetXaxis()->SetBinLabel(2, "posZ cut");
    MCeventHist.get<TH1>(HIST("hGenMCEventSelection"))->GetXaxis()->SetBinLabel(3, "INEL>0 cut");
    MCeventHist.get<TH1>(HIST("hGenMCEventSelection"))->GetXaxis()->SetBinLabel(4, "With at least a reco coll");
    MCeventHist.get<TH1>(HIST("hGenMCEventSelection"))->GetXaxis()->SetBinLabel(5, "With at least a #phi");

    MCeventHist.add("hGenMCVertexZ", "hGenMCVertexZ", kTH1F, {vertexZAxis});
    MCeventHist.add("hGenMCMultiplicityPercent", "GenMC Multiplicity Percentile", kTH1F, {binnedmultAxis});

    // Phi tpological/PID cuts
    PhicandHist.add("hEta", "Eta distribution", kTH1F, {{200, -1.0f, 1.0f}});
    PhicandHist.add("hDcaxy", "Dcaxy distribution", kTH1F, {{200, -1.0f, 1.0f}});
    PhicandHist.add("hDcaz", "Dcaz distribution", kTH1F, {{200, -1.0f, 1.0f}});
    PhicandHist.add("hNsigmaKaonTPC", "NsigmaKaon TPC distribution", kTH2F, {ptAxis, {100, -10.0f, 10.0f}});
    PhicandHist.add("hNsigmaKaonTOF", "NsigmaKaon TOF distribution", kTH2F, {ptAxis, {100, -10.0f, 10.0f}});

    // K0S topological/PID cuts
    K0SHist.add("hDCAV0Daughters", "hDCAV0Daughters", kTH1F, {{55, 0.0f, 2.2f}});
    K0SHist.add("hV0CosPA", "hV0CosPA", kTH1F, {{100, 0.95f, 1.f}});
    K0SHist.add("hNSigmaPosPionFromK0S", "hNSigmaPosPionFromK0Short", kTH2F, {ptAxis, {100, -5.f, 5.f}});
    K0SHist.add("hNSigmaNegPionFromK0S", "hNSigmaNegPionFromK0Short", kTH2F, {ptAxis, {100, -5.f, 5.f}});

    // Phi invariant mass for computing purities and normalisation
    PhipurHist.add("h2PhipurInvMass", "Invariant mass of Phi for Purity (no K0S/Pi)", kTH2F, {multAxis, PhimassAxis});
    PhipurHist.add("h2PhipurK0SInvMassInclusive", "Invariant mass of Phi for Purity (K0S) Inclusive", kTH2F, {multAxis, PhimassAxis});
    PhipurHist.add("h2PhipurK0SInvMassFirstCut", "Invariant mass of Phi for Purity (K0S) Deltay < FirstCut", kTH2F, {multAxis, PhimassAxis});
    PhipurHist.add("h2PhipurK0SInvMassSecondCut", "Invariant mass of Phi for Purity (K0S) Deltay < SecondCut", kTH2F, {multAxis, PhimassAxis});
    PhipurHist.add("h3PhipurPiInvMassInclusive", "Invariant mass of Phi for Purity (Pi) Inclusive", kTH3F, {multAxis, ptAxis, PhimassAxis});
    PhipurHist.add("h3PhipurPiInvMassFirstCut", "Invariant mass of Phi for Purity (Pi) Deltay < FirstCut", kTH3F, {multAxis, ptAxis, PhimassAxis});
    PhipurHist.add("h3PhipurPiInvMassSecondCut", "Invariant mass of Phi for Purity (Pi) Deltay < SecondCut", kTH3F, {multAxis, ptAxis, PhimassAxis});

    // 2D mass for Phi and K0S for Same Event and Mixed Event
    for (int i = 0; i < nMultBin; i++) {
      PhiK0SHist.add(PhiK0SSEInc[i].data(), "2D Invariant mass of Phi and K0Short for Same Event Inclusive", kTH2F, {K0SmassAxis, cfgPhimassAxisInc.at(i)});
      PhiK0SHist.add(PhiK0SSEFCut[i].data(), "2D Invariant mass of Phi and K0Short for Same Event Deltay < FirstCut", kTH2F, {K0SmassAxis, cfgPhimassAxisFCut.at(i)});
      PhiK0SHist.add(PhiK0SSESCut[i].data(), "2D Invariant mass of Phi and K0Short for Same Event Deltay < SecondCut", kTH2F, {K0SmassAxis, cfgPhimassAxisSCut.at(i)});
    }

    PhiK0SHist.add("h3PhiK0SInvMassMixedEventInclusive", "2D Invariant mass of Phi and K0Short for Mixed Event Inclusive", kTH3F, {multAxis, K0SmassAxis, PhimassAxis});
    PhiK0SHist.add("h3PhiK0SInvMassMixedEventFirstCut", "2D Invariant mass of Phi and K0Short for Mixed Event Deltay < FirstCut", kTH3F, {multAxis, K0SmassAxis, PhimassAxis});
    PhiK0SHist.add("h3PhiK0SInvMassMixedEventSecondCut", "2D Invariant mass of Phi and K0Short for Mixed Event Deltay < SecondCut", kTH3F, {multAxis, K0SmassAxis, PhimassAxis});

    // MC 2D mass for Phi and K0S
    for (int i = 0; i < nMultBin; i++) {
      MCPhiK0SHist.add(MCPhiK0SSEInc[i].data(), "2D Invariant mass of Phi and K0Short for RecMC Inclusive", kTH2F, {K0SmassAxis, cfgPhimassAxisInc.at(i)});
      MCPhiK0SHist.add(MCPhiK0SSEFCut[i].data(), "2D Invariant mass of Phi and K0Short for RecMC Deltay < FirstCut", kTH2F, {K0SmassAxis, cfgPhimassAxisFCut.at(i)});
      MCPhiK0SHist.add(MCPhiK0SSESCut[i].data(), "2D Invariant mass of Phi and K0Short for RecMC Deltay < SecondCut", kTH2F, {K0SmassAxis, cfgPhimassAxisSCut.at(i)});
    }

    // GenMC K0S coupled to Phi
    MCPhiK0SHist.add("h1PhiK0SGenMCInclusive", "K0Short coupled to Phi for GenMC Inclusive", kTH1F, {{10, -0.5f, 9.5f}});
    MCPhiK0SHist.add("h1PhiK0SGenMCFirstCut", "K0Short coupled to Phi for GenMC Deltay < FirstCut", kTH1F, {{10, -0.5f, 9.5f}});
    MCPhiK0SHist.add("h1PhiK0SGenMCSecondCut", "K0Short coupled to Phi for GenMC Deltay < SecondCut", kTH1F, {{10, -0.5f, 9.5f}});

    MCPhiK0SHist.add("h1PhiK0SGenMCInclusiveAssocReco", "K0Short coupled to Phi for GenMC Inclusive Associated Reco Collision", kTH1F, {{10, -0.5f, 9.5f}});
    MCPhiK0SHist.add("h1PhiK0SGenMCFirstCutAssocReco", "K0Short coupled to Phi for GenMC Deltay < FirstCut Associated Reco Collision", kTH1F, {{10, -0.5f, 9.5f}});
    MCPhiK0SHist.add("h1PhiK0SGenMCSecondCutAssocReco", "K0Short coupled to Phi for GenMC Deltay < SecondCut Associated Reco Collision", kTH1F, {{10, -0.5f, 9.5f}});

    // Phi mass vs Pion NSigma dE/dx for Same Event and Mixed Event
    for (int i = 0; i < nMultBin; i++) {
      PhiPionHist.add(PhiPiSEInc[i].data(), "Phi Invariant mass vs Pion nSigma TPC/TOF for Same Event Inclusive", kTHnSparseF, {binnedptAxis, {100, -10.0f, 10.0f}, {100, -10.0f, 10.0f}, cfgPhimassAxisInc.at(i)});
      PhiPionHist.add(PhiPiSEFCut[i].data(), "Phi Invariant mass vs Pion nSigma TPC/TOF for Same Event Deltay < FirstCut", kTHnSparseF, {binnedptAxis, {100, -10.0f, 10.0f}, {100, -10.0f, 10.0f}, cfgPhimassAxisFCut.at(i)});
      PhiPionHist.add(PhiPiSESCut[i].data(), "Phi Invariant mass vs Pion nSigma TPC/TOF for Same Event Deltay < SecondCut", kTHnSparseF, {binnedptAxis, {100, -10.0f, 10.0f}, {100, -10.0f, 10.0f}, cfgPhimassAxisSCut.at(i)});
    }

    PhiPionHist.add("h4PhiInvMassPiNSigmadEdxMixedEventInclusive", "Phi Invariant mass vs Pion nSigma TPC/TOF for Mixed Event Inclusive", kTHnSparseF, {multAxis, binnedptAxis, {100, -10.0f, 10.0f}, {100, -10.0f, 10.0f}, PhimassAxis});
    PhiPionHist.add("h4PhiInvMassPiNSigmadEdxMixedEventFirstCut", "Phi Invariant mass vs Pion nSigma TPC/TOF for Mixed Event Deltay < FirstCut", kTHnSparseF, {multAxis, binnedptAxis, {100, -10.0f, 10.0f}, {100, -10.0f, 10.0f}, PhimassAxis});
    PhiPionHist.add("h4PhiInvMassPiNSigmadEdxMixedEventSecondCut", "Phi Invariant mass vs Pion nSigma TPC/TOF for Mixed Event Deltay < SecondCut", kTHnSparseF, {multAxis, binnedptAxis, {100, -10.0f, 10.0f}, {100, -10.0f, 10.0f}, PhimassAxis});

    // MC Phi mass vs Pion NSigma dE/dx
    for (int i = 0; i < nMultBin; i++) {
      MCPhiPionHist.add(MCPhiPiSEInc[i].data(), "Phi Invariant mass vs Pion nSigma TPC/TOF for RecMC Inclusive", kTHnSparseF, {binnedptAxis, {100, -10.0f, 10.0f}, {100, -10.0f, 10.0f}, cfgPhimassAxisInc.at(i)});
      MCPhiPionHist.add(MCPhiPiSEFCut[i].data(), "Phi Invariant mass vs Pion nSigma TPC/TOF for RecMC Deltay < FirstCut", kTHnSparseF, {binnedptAxis, {100, -10.0f, 10.0f}, {100, -10.0f, 10.0f}, cfgPhimassAxisFCut.at(i)});
      MCPhiPionHist.add(MCPhiPiSESCut[i].data(), "Phi Invariant mass vs Pion nSigma TPC/TOF for RecMC Deltay < SecondCut", kTHnSparseF, {binnedptAxis, {100, -10.0f, 10.0f}, {100, -10.0f, 10.0f}, cfgPhimassAxisSCut.at(i)});
    }

    // GenMC Pion coupled to Phi
    MCPhiPionHist.add("h2PhiPiGenMCInclusive", "Pion coupled to Phi for GenMC Inclusive", kTH2F, {{10, -0.5f, 9.5f}, {3, -0.5f, 2.5f}});
    MCPhiPionHist.add("h2PhiPiGenMCFirstCut", "Pion coupled to Phi for GenMC Deltay < FirstCut", kTH2F, {{10, -0.5f, 9.5f}, {3, -0.5f, 2.5f}});
    MCPhiPionHist.add("h2PhiPiGenMCSecondCut", "Pion coupled to Phi for GenMC Deltay < SecondCut", kTH2F, {{10, -0.5f, 9.5f}, {3, -0.5f, 2.5f}});

    MCPhiPionHist.add("h2PhiPiGenMCInclusiveAssocReco", "Pion coupled to Phi for GenMC Inclusive Associated Reco Collision", kTH2F, {{10, -0.5f, 9.5f}, {3, -0.5f, 2.5f}});
    MCPhiPionHist.add("h2PhiPiGenMCFirstCutAssocReco", "Pion coupled to Phi for GenMC Deltay < FirstCut Associated Reco Collision", kTH2F, {{10, -0.5f, 9.5f}, {3, -0.5f, 2.5f}});
    MCPhiPionHist.add("h2PhiPiGenMCSecondCutAssocReco", "Pion coupled to Phi for GenMC Deltay < SecondCut Associated Reco Collision", kTH2F, {{10, -0.5f, 9.5f}, {3, -0.5f, 2.5f}});

    // MCPhi invariant mass for computing efficiencies and MCnormalisation
    PhieffHist.add("h2PhieffInvMass", "Invariant mass of Phi for Efficiency (no K0S/Pi)", kTH2F, {multAxis, PhimassAxis});
    PhieffHist.add("h2PhieffK0SInvMassInclusive", "Invariant mass of Phi for Efficiency (K0S) Inclusive", kTH2F, {multAxis, PhimassAxis});
    PhieffHist.add("h2PhieffK0SInvMassFirstCut", "Invariant mass of Phi for Efficiency (K0S) Deltay < FirstCut", kTH2F, {multAxis, PhimassAxis});
    PhieffHist.add("h2PhieffK0SInvMassSecondCut", "Invariant mass of Phi for Efficiency (K0S) Deltay < SecondCut", kTH2F, {multAxis, PhimassAxis});
    PhieffHist.add("h3PhieffPiInvMassInclusive", "Invariant mass of Phi for Efficiency (Pi) Inclusive", kTH3F, {multAxis, ptAxis, PhimassAxis});
    PhieffHist.add("h3PhieffPiInvMassFirstCut", "Invariant mass of Phi for Efficiency (Pi) Deltay < FirstCut", kTH3F, {multAxis, ptAxis, PhimassAxis});
    PhieffHist.add("h3PhieffPiInvMassSecondCut", "Invariant mass of Phi for Efficiency (Pi) Deltay < SecondCut", kTH3F, {multAxis, ptAxis, PhimassAxis});

    // GenMC Phi and Phi coupled to K0S and Pion
    PhieffHist.add("h1PhiGenMC", "Phi for GenMC", kTH1F, {{10, -0.5f, 9.5f}});
    PhieffHist.add("h1PhieffK0SGenMCInclusive", "Phi coupled to K0Short for GenMC Inclusive", kTH1F, {{10, -0.5f, 9.5f}});
    PhieffHist.add("h1PhieffK0SGenMCFirstCut", "Phi coupled to K0Short for GenMC Deltay < FirstCut", kTH1F, {{10, -0.5f, 9.5f}});
    PhieffHist.add("h1PhieffK0SGenMCSecondCut", "Phi coupled to K0Short for GenMC Deltay < SecondCut", kTH1F, {{10, -0.5f, 9.5f}});
    PhieffHist.add("h2PhieffPiGenMCInclusive", "Phi coupled to Pion for GenMC Inclusive", kTH2F, {{10, -0.5f, 9.5f}, {3, -0.5f, 2.5f}});
    PhieffHist.add("h2PhieffPiGenMCFirstCut", "Phi coupled to Pion for GenMC Deltay < FirstCut", kTH2F, {{10, -0.5f, 9.5f}, {3, -0.5f, 2.5f}});
    PhieffHist.add("h2PhieffPiGenMCSecondCut", "Phi coupled to Pion for GenMC Deltay < SecondCut", kTH2F, {{10, -0.5f, 9.5f}, {3, -0.5f, 2.5f}});
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
      if (std::abs(collision.posZ()) > cutzvertex)
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
      if (cfgiskNoITSROFrameBorder) {
        if (!collision.selection_bit(aod::evsel::kNoITSROFrameBorder))
          return false;
        if (QA)
          MCeventHist.fill(HIST("hRecMCEventSelection"), 3); // kNoITSROFrameBorder collisions (by default not requested by the selection)
      }
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
    if (std::abs(track.pt()) < cMinPtcut)
      return false;
    if (std::abs(track.dcaXY()) > cMaxDCArToPVcut)
      return false;
    if (std::abs(track.dcaZ()) > cMaxDCAzToPVcut)
      return false;

    if (cfgPrimaryTrack && !track.isPrimaryTrack())
      return false;
    if (cfgGlobalWoDCATrack && !track.isGlobalTrackWoDCA())
      return false;
    if (cfgPVContributor && !track.isPVContributor())
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

    if (track.pt() < 1.2) {
      if (!track.hasTPC())
        return false;
      if (track.tpcNClsFound() < minTPCnClsFound)
        return false;
      if (track.tpcNClsCrossedRows() < minNCrossedRowsTPC)
        return false;
      if (track.tpcChi2NCl() > maxChi2TPC)
        return false;
    }

    if (track.pt() > 0.5) {
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
  template <bool isMix, bool isMC, std::size_t iBin>
  void fillInvMass2D(TLorentzVector V0, const std::vector<TLorentzVector> listPhi, float multiplicity, double weightInclusive, double weightLtFirstCut, double weightLtSecondCut)
  {
    double massV0 = V0.M();
    double rapidityV0 = V0.Rapidity();

    for (unsigned int phitag = 0; phitag < listPhi.size(); phitag++) {
      double massPhi = listPhi[phitag].M();
      double rapidityPhi = listPhi[phitag].Rapidity();
      double deltay = std::abs(rapidityV0 - rapidityPhi);

      if constexpr (!isMix) { // same event
        PhiK0SHist.fill(HIST(PhiK0SSEInc[iBin]), massV0, massPhi, weightInclusive);
        if (deltay > cfgFirstCutonDeltay)
          continue;
        PhiK0SHist.fill(HIST(PhiK0SSEFCut[iBin]), massV0, massPhi, weightLtFirstCut);
        if (deltay > cfgSecondCutonDeltay)
          continue;
        PhiK0SHist.fill(HIST(PhiK0SSESCut[iBin]), massV0, massPhi, weightLtSecondCut);
      } else { // mixed event
        PhiK0SHist.fill(HIST("h3PhiK0SInvMassMixedEventInclusive"), multiplicity, massV0, massPhi, weightInclusive);
        if (deltay > cfgFirstCutonDeltay)
          continue;
        PhiK0SHist.fill(HIST("h3PhiK0SInvMassMixedEventFirstCut"), multiplicity, massV0, massPhi, weightLtFirstCut);
        if (deltay > cfgSecondCutonDeltay)
          continue;
        PhiK0SHist.fill(HIST("h3PhiK0SInvMassMixedEventSecondCut"), multiplicity, massV0, massPhi, weightLtSecondCut);
      }

      if constexpr (isMC) { // MC event
        MCPhiK0SHist.fill(HIST(MCPhiK0SSEInc[iBin]), massV0, massPhi, weightInclusive);
        if (deltay > cfgFirstCutonDeltay)
          continue;
        MCPhiK0SHist.fill(HIST(MCPhiK0SSEFCut[iBin]), massV0, massPhi, weightLtFirstCut);
        if (deltay > cfgSecondCutonDeltay)
          continue;
        MCPhiK0SHist.fill(HIST(MCPhiK0SSESCut[iBin]), massV0, massPhi, weightLtSecondCut);
      }
    }
  }

  // Fill Phi invariant mass vs Pion nSigmadE/dx histogram
  template <bool isMix, bool isMC, std::size_t iBin>
  void fillInvMassNSigma(TLorentzVector Pi, float nSigmaTPCPi, float nSigmaTOFPi, const std::vector<TLorentzVector> listPhi, float multiplicity, double weightInclusive, double weightLtFirstCut, double weightLtSecondCut)
  {
    double rapidityPi = Pi.Rapidity();
    double ptPi = Pi.Pt();

    for (unsigned int phitag = 0; phitag < listPhi.size(); phitag++) {
      double massPhi = listPhi[phitag].M();
      double rapidityPhi = listPhi[phitag].Rapidity();
      double deltay = std::abs(rapidityPi - rapidityPhi);

      if constexpr (!isMix) { // same event
        PhiPionHist.fill(HIST(PhiPiSEInc[iBin]), ptPi, nSigmaTPCPi, nSigmaTOFPi, massPhi, weightInclusive);
        if (deltay > cfgFirstCutonDeltay)
          continue;
        PhiPionHist.fill(HIST(PhiPiSEFCut[iBin]), ptPi, nSigmaTPCPi, nSigmaTOFPi, massPhi, weightLtFirstCut);
        if (deltay > cfgSecondCutonDeltay)
          continue;
        PhiPionHist.fill(HIST(PhiPiSESCut[iBin]), ptPi, nSigmaTPCPi, nSigmaTOFPi, massPhi, weightLtSecondCut);
      } else { // mixed event
        PhiPionHist.fill(HIST("h4PhiInvMassPiNSigmadEdxMixedEventInclusive"), multiplicity, ptPi, nSigmaTPCPi, nSigmaTOFPi, massPhi, weightInclusive);
        if (deltay > cfgFirstCutonDeltay)
          continue;
        PhiPionHist.fill(HIST("h4PhiInvMassPiNSigmadEdxMixedEventFirstCut"), multiplicity, ptPi, nSigmaTPCPi, nSigmaTOFPi, massPhi, weightLtFirstCut);
        if (deltay > cfgSecondCutonDeltay)
          continue;
        PhiPionHist.fill(HIST("h4PhiInvMassPiNSigmadEdxMixedEventSecondCut"), multiplicity, ptPi, nSigmaTPCPi, nSigmaTOFPi, massPhi, weightLtSecondCut);
      }

      if constexpr (isMC) { // MC event
        MCPhiPionHist.fill(HIST(MCPhiPiSEInc[iBin]), ptPi, nSigmaTPCPi, nSigmaTOFPi, massPhi, weightInclusive);
        if (deltay > cfgFirstCutonDeltay)
          continue;
        MCPhiPionHist.fill(HIST(MCPhiPiSEFCut[iBin]), ptPi, nSigmaTPCPi, nSigmaTOFPi, massPhi, weightLtFirstCut);
        if (deltay > cfgSecondCutonDeltay)
          continue;
        MCPhiPionHist.fill(HIST(MCPhiPiSESCut[iBin]), ptPi, nSigmaTPCPi, nSigmaTOFPi, massPhi, weightLtSecondCut);
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

    for (auto track1 : posThisColl) { // loop over all selected tracks
      if (!selectionTrackResonance(track1) || !selectionPIDKaon(track1))
        continue; // topological and PID selection

      PhicandHist.fill(HIST("hEta"), track1.eta());
      PhicandHist.fill(HIST("hDcaxy"), track1.dcaXY());
      PhicandHist.fill(HIST("hDcaz"), track1.dcaZ());
      PhicandHist.fill(HIST("hNsigmaKaonTPC"), track1.tpcInnerParam(), track1.tpcNSigmaKa());
      PhicandHist.fill(HIST("hNsigmaKaonTOF"), track1.tpcInnerParam(), track1.tofNSigmaKa());

      auto track1ID = track1.globalIndex();

      // Loop over all negative candidates
      for (auto track2 : negThisColl) {
        if (!selectionTrackResonance(track2) || !selectionPIDKaon(track2))
          continue; // topological and PID selection

        auto track2ID = track2.globalIndex();
        if (track2ID == track1ID)
          continue; // condition to avoid double counting of pair

        TLorentzVector recPhi;
        recPhi = recMother(track1, track2, massKa, massKa);
        if (recPhi.Rapidity() > 0.8)
          continue;

        if (!isCountedPhi) {
          eventHist.fill(HIST("hEventSelection"), 4); // at least a Phi candidate in the event
          isCountedPhi = true;
        }

        PhipurHist.fill(HIST("h2PhipurInvMass"), multiplicity, recPhi.M());

        bool isCountedK0SInclusive = false, isCountedK0SFirstCut = false, isCountedK0SSecondCut = false;

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

          TLorentzVector recK0S;
          recK0S.SetXYZM(v0.px(), v0.py(), v0.pz(), v0.mK0Short());

          if (recK0S.Rapidity() > 0.8)
            continue;
          if (!isCountedK0SInclusive) {
            PhipurHist.fill(HIST("h2PhipurK0SInvMassInclusive"), multiplicity, recPhi.M());
            isCountedK0SInclusive = true;
          }
          if (std::abs(recK0S.Rapidity() - recPhi.Rapidity()) > cfgFirstCutonDeltay)
            continue;
          if (!isCountedK0SFirstCut) {
            PhipurHist.fill(HIST("h2PhipurK0SInvMassFirstCut"), multiplicity, recPhi.M());
            isCountedK0SFirstCut = true;
          }
          if (std::abs(recK0S.Rapidity() - recPhi.Rapidity()) > cfgSecondCutonDeltay)
            continue;
          if (!isCountedK0SSecondCut) {
            PhipurHist.fill(HIST("h2PhipurK0SInvMassSecondCut"), multiplicity, recPhi.M());
            isCountedK0SSecondCut = true;
          }
        }
        isFilledhV0 = true;

        bool isCountedPiInclusive = false, isCountedPiFirstCut = false, isCountedPiSecondCut = false;

        // Loop over all primary pion candidates
        for (const auto& track : fullTracks) {
          if (!selectionPion(track))
            continue;

          TLorentzVector recPi;
          recPi.SetXYZM(track.px(), track.py(), track.pz(), massPi);

          if (recPi.Rapidity() > 0.8)
            continue;
          if (!isCountedPiInclusive) {
            PhipurHist.fill(HIST("h3PhipurPiInvMassInclusive"), multiplicity, recPi.Pt(), recPhi.M());
            isCountedPiInclusive = true;
          }
          if (std::abs(recPi.Rapidity() - recPhi.Rapidity()) > cfgFirstCutonDeltay)
            continue;
          if (!isCountedPiFirstCut) {
            PhipurHist.fill(HIST("h3PhipurPiInvMassFirstCut"), multiplicity, recPi.Pt(), recPhi.M());
            isCountedPiFirstCut = true;
          }
          if (std::abs(recPi.Rapidity() - recPhi.Rapidity()) > cfgSecondCutonDeltay)
            continue;
          if (!isCountedPiSecondCut) {
            PhipurHist.fill(HIST("h3PhipurPiInvMassSecondCut"), multiplicity, recPi.Pt(), recPhi.M());
            isCountedPiSecondCut = true;
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

    int iBin = 0;
    for (int i = 0; i < nMultBin; i++) {
      if (multBin[i] < multiplicity && multiplicity <= multBin[i + 1]) {
        iBin = i;
        break;
      }
    }

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

      TLorentzVector recK0S;
      recK0S.SetXYZM(v0.px(), v0.py(), v0.pz(), v0.mK0Short());
      if (recK0S.Rapidity() > 0.8)
        continue;

      std::vector<TLorentzVector> listrecPhi;
      int countInclusive = 0, countLtFirstCut = 0, countLtSecondCut = 0;

      // Phi reconstruction
      // Loop over positive candidates
      for (auto track1 : posThisColl) { // loop over all selected tracks
        if (!selectionTrackResonance(track1) || !selectionPIDKaon(track1))
          continue; // topological and PID selection

        auto track1ID = track1.globalIndex();

        // Loop over all negative candidates
        for (auto track2 : negThisColl) {
          if (!selectionTrackResonance(track2) || !selectionPIDKaon(track2))
            continue; // topological and PID selection

          auto track2ID = track2.globalIndex();
          if (track2ID == track1ID)
            continue; // condition to avoid double counting of pair

          TLorentzVector recPhi;
          recPhi = recMother(track1, track2, massKa, massKa);
          if (recPhi.Rapidity() > 0.8)
            continue;

          listrecPhi.push_back(recPhi);

          if (lowmPhiInc->at(iBin) <= recPhi.M() && recPhi.M() <= upmPhiInc->at(iBin))
            countInclusive++;
          if (std::abs(recK0S.Rapidity() - recPhi.Rapidity()) > cfgFirstCutonDeltay)
            continue;
          if (lowmPhiFCut->at(iBin) <= recPhi.M() && recPhi.M() <= upmPhiFCut->at(iBin))
            countLtFirstCut++;
          if (std::abs(recK0S.Rapidity() - recPhi.Rapidity()) > cfgSecondCutonDeltay)
            continue;
          if (lowmPhiSCut->at(iBin) <= recPhi.M() && recPhi.M() <= upmPhiSCut->at(iBin))
            countLtSecondCut++;
        }
      }

      float weightInclusive, weightLtFirstCut, weightLtSecondCut;
      if (countInclusive > 0) {
        weightInclusive = 1. / static_cast<float>(countInclusive);
      } else {
        weightInclusive = 0;
        eventHist.fill(HIST("thereisnoPhiwK0S"), 0);
      }
      if (countLtFirstCut > 0) {
        weightLtFirstCut = 1. / static_cast<float>(countLtFirstCut);
      } else {
        weightLtFirstCut = 0;
        eventHist.fill(HIST("thereisnoPhiwK0S"), 1);
      }
      if (countLtSecondCut > 0) {
        weightLtSecondCut = 1. / static_cast<float>(countLtSecondCut);
      } else {
        weightLtSecondCut = 0;
        eventHist.fill(HIST("thereisnoPhiwK0S"), 2);
      }

      switch (iBin) {
        case 0: {
          fillInvMass2D<false, false, 0>(recK0S, listrecPhi, multiplicity, weightInclusive, weightLtFirstCut, weightLtSecondCut);
          break;
        }
        case 1: {
          fillInvMass2D<false, false, 1>(recK0S, listrecPhi, multiplicity, weightInclusive, weightLtFirstCut, weightLtSecondCut);
          break;
        }
        case 2: {
          fillInvMass2D<false, false, 2>(recK0S, listrecPhi, multiplicity, weightInclusive, weightLtFirstCut, weightLtSecondCut);
          break;
        }
        case 3: {
          fillInvMass2D<false, false, 3>(recK0S, listrecPhi, multiplicity, weightInclusive, weightLtFirstCut, weightLtSecondCut);
          break;
        }
        case 4: {
          fillInvMass2D<false, false, 4>(recK0S, listrecPhi, multiplicity, weightInclusive, weightLtFirstCut, weightLtSecondCut);
          break;
        }
        case 5: {
          fillInvMass2D<false, false, 5>(recK0S, listrecPhi, multiplicity, weightInclusive, weightLtFirstCut, weightLtSecondCut);
          break;
        }
        case 6: {
          fillInvMass2D<false, false, 6>(recK0S, listrecPhi, multiplicity, weightInclusive, weightLtFirstCut, weightLtSecondCut);
          break;
        }
        case 7: {
          fillInvMass2D<false, false, 7>(recK0S, listrecPhi, multiplicity, weightInclusive, weightLtFirstCut, weightLtSecondCut);
          break;
        }
        case 8: {
          fillInvMass2D<false, false, 8>(recK0S, listrecPhi, multiplicity, weightInclusive, weightLtFirstCut, weightLtSecondCut);
          break;
        }
        case 9: {
          fillInvMass2D<false, false, 9>(recK0S, listrecPhi, multiplicity, weightInclusive, weightLtFirstCut, weightLtSecondCut);
          break;
        }
        default:
          break;
      }
    }
  }

  PROCESS_SWITCH(phik0shortanalysis, processSEPhiK0S, "Process Same Event for Phi-K0S Analysis", false);

  void processSEPhiPion(soa::Filtered<SelCollisions>::iterator const& collision, FullTracks const& fullTracks)
  {
    if (!collision.isInelGt0())
      return;

    float multiplicity = collision.centFT0M();
    eventHist.fill(HIST("hMultiplicityPercent"), multiplicity);

    int iBin = 0;
    for (int i = 0; i < nMultBin; i++) {
      if (multBin[i] < multiplicity && multiplicity <= multBin[i + 1]) {
        iBin = i;
        break;
      }
    }

    // Defining positive and negative tracks for phi reconstruction
    auto posThisColl = posTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto negThisColl = negTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);

    // Loop over all primary pion candidates
    for (const auto& track : fullTracks) {

      // Pion selection
      if (!selectionPion(track))
        continue;

      TLorentzVector recPi;
      recPi.SetXYZM(track.px(), track.py(), track.pz(), massPi);
      if (recPi.Rapidity() > 0.8)
        continue;

      std::vector<TLorentzVector> listrecPhi;
      int countInclusive = 0, countLtFirstCut = 0, countLtSecondCut = 0;

      // Phi reconstruction
      // Loop over positive candidates
      for (auto track1 : posThisColl) { // loop over all selected tracks
        if (!selectionTrackResonance(track1) || !selectionPIDKaon(track1))
          continue; // topological and PID selection

        auto track1ID = track1.globalIndex();

        // Loop over all negative candidates
        for (auto track2 : negThisColl) {
          if (!selectionTrackResonance(track2) || !selectionPIDKaon(track2))
            continue; // topological and PID selection

          auto track2ID = track2.globalIndex();
          if (track2ID == track1ID)
            continue; // condition to avoid double counting of pair

          TLorentzVector recPhi;
          recPhi = recMother(track1, track2, massKa, massKa);
          if (recPhi.Rapidity() > 0.8)
            continue;

          listrecPhi.push_back(recPhi);

          if (lowmPhiInc->at(iBin) <= recPhi.M() && recPhi.M() <= upmPhiInc->at(iBin))
            countInclusive++;
          if (std::abs(recPi.Rapidity() - recPhi.Rapidity()) > cfgFirstCutonDeltay)
            continue;
          if (lowmPhiFCut->at(iBin) <= recPhi.M() && recPhi.M() <= upmPhiFCut->at(iBin))
            countLtFirstCut++;
          if (std::abs(recPi.Rapidity() - recPhi.Rapidity()) > cfgSecondCutonDeltay)
            continue;
          if (lowmPhiSCut->at(iBin) <= recPhi.M() && recPhi.M() <= upmPhiSCut->at(iBin))
            countLtSecondCut++;
        }
      }

      float nsigmaTPC, nsigmaTOF;
      if (track.hasTPC())
        nsigmaTPC = track.tpcNSigmaPi();
      else
        nsigmaTPC = -9.99;
      if (track.hasTOF())
        nsigmaTOF = track.tofNSigmaPi();
      else
        nsigmaTOF = -9.99;

      float weightInclusive, weightLtFirstCut, weightLtSecondCut;
      if (countInclusive > 0) {
        weightInclusive = 1. / static_cast<float>(countInclusive);
      } else {
        weightInclusive = 0;
        eventHist.fill(HIST("thereisnoPhiwPi"), 0);
      }
      if (countLtFirstCut > 0) {
        weightLtFirstCut = 1. / static_cast<float>(countLtFirstCut);
      } else {
        weightLtFirstCut = 0;
        eventHist.fill(HIST("thereisnoPhiwPi"), 1);
      }
      if (countLtSecondCut > 0) {
        weightLtSecondCut = 1. / static_cast<float>(countLtSecondCut);
      } else {
        weightLtSecondCut = 0;
        eventHist.fill(HIST("thereisnoPhiwPi"), 2);
      }

      switch (iBin) {
        case 0: {
          fillInvMassNSigma<false, false, 0>(recPi, nsigmaTPC, nsigmaTOF, listrecPhi, multiplicity, weightInclusive, weightLtFirstCut, weightLtSecondCut);
          break;
        }
        case 1: {
          fillInvMassNSigma<false, false, 1>(recPi, nsigmaTPC, nsigmaTOF, listrecPhi, multiplicity, weightInclusive, weightLtFirstCut, weightLtSecondCut);
          break;
        }
        case 2: {
          fillInvMassNSigma<false, false, 2>(recPi, nsigmaTPC, nsigmaTOF, listrecPhi, multiplicity, weightInclusive, weightLtFirstCut, weightLtSecondCut);
          break;
        }
        case 3: {
          fillInvMassNSigma<false, false, 3>(recPi, nsigmaTPC, nsigmaTOF, listrecPhi, multiplicity, weightInclusive, weightLtFirstCut, weightLtSecondCut);
          break;
        }
        case 4: {
          fillInvMassNSigma<false, false, 4>(recPi, nsigmaTPC, nsigmaTOF, listrecPhi, multiplicity, weightInclusive, weightLtFirstCut, weightLtSecondCut);
          break;
        }
        case 5: {
          fillInvMassNSigma<false, false, 5>(recPi, nsigmaTPC, nsigmaTOF, listrecPhi, multiplicity, weightInclusive, weightLtFirstCut, weightLtSecondCut);
          break;
        }
        case 6: {
          fillInvMassNSigma<false, false, 6>(recPi, nsigmaTPC, nsigmaTOF, listrecPhi, multiplicity, weightInclusive, weightLtFirstCut, weightLtSecondCut);
          break;
        }
        case 7: {
          fillInvMassNSigma<false, false, 7>(recPi, nsigmaTPC, nsigmaTOF, listrecPhi, multiplicity, weightInclusive, weightLtFirstCut, weightLtSecondCut);
          break;
        }
        case 8: {
          fillInvMassNSigma<false, false, 8>(recPi, nsigmaTPC, nsigmaTOF, listrecPhi, multiplicity, weightInclusive, weightLtFirstCut, weightLtSecondCut);
          break;
        }
        case 9: {
          fillInvMassNSigma<false, false, 9>(recPi, nsigmaTPC, nsigmaTOF, listrecPhi, multiplicity, weightInclusive, weightLtFirstCut, weightLtSecondCut);
          break;
        }
        default:
          break;
      }
    }
  }

  PROCESS_SWITCH(phik0shortanalysis, processSEPhiPion, "Process Same Event for Phi-Pion Analysis", false);

  void processMEPhiK0S(soa::Filtered<SelCollisions> const& collisions, FullTracks const&, FullV0s const& V0s, V0DauTracks const&)
  {
    // Mixing the events with similar vertex z and multiplicity
    BinningTypeVertexContributor binningOnPositions{{axisVertex, axisMultiplicity}, true};
    for (auto const& [collision1, collision2] : o2::soa::selfCombinations(binningOnPositions, cfgNoMixedEvents, -1, collisions, collisions)) {
      if (!collision1.isInelGt0() || !collision2.isInelGt0())
        continue;

      float multiplicity = collision1.centFT0M();

      int iBin = 0;
      for (int i = 0; i < nMultBin; i++) {
        if (multBin[i] < multiplicity && multiplicity <= multBin[i + 1]) {
          iBin = i;
          break;
        }
      }

      // Defining V0s from collision1
      auto V0ThisColl = V0s.sliceByCached(aod::v0::collisionId, collision1.globalIndex(), cache);

      // Defining positive and negative tracks for phi reconstruction from collision1 and collision2, respectively
      auto posThisColl = posTracks->sliceByCached(aod::track::collisionId, collision1.globalIndex(), cache);
      auto negThisColl = negTracks->sliceByCached(aod::track::collisionId, collision2.globalIndex(), cache);

      for (const auto& v0 : V0ThisColl) {
        const auto& posDaughterTrack = v0.posTrack_as<V0DauTracks>();
        const auto& negDaughterTrack = v0.negTrack_as<V0DauTracks>();

        // Cut on V0 dynamic columns
        if (!selectionV0(v0, posDaughterTrack, negDaughterTrack))
          continue;

        TLorentzVector recK0S;
        recK0S.SetXYZM(v0.px(), v0.py(), v0.pz(), v0.mK0Short());
        if (recK0S.Rapidity() > 0.8)
          continue;

        std::vector<TLorentzVector> listrecPhi;
        int countInclusive = 0, countLtFirstCut = 0, countLtSecondCut = 0;

        // Combinatorial background simulation
        for (auto& [track1, track2] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(posThisColl, negThisColl))) {
          if (!selectionTrackResonance(track1) || !selectionPIDKaon(track1) || !selectionTrackResonance(track2) || !selectionPIDKaon(track2))
            continue; // topological and PID selection

          TLorentzVector recPhi;
          recPhi = recMother(track1, track2, massKa, massKa);
          if (recPhi.Rapidity() > 0.8)
            continue;

          listrecPhi.push_back(recPhi);

          countInclusive++;
          if (std::abs(recK0S.Rapidity() - recPhi.Rapidity()) > cfgFirstCutonDeltay)
            continue;
          countLtFirstCut++;
          if (std::abs(recK0S.Rapidity() - recPhi.Rapidity()) > cfgSecondCutonDeltay)
            continue;
          countLtSecondCut++;
        }

        float weightInclusive = 1. / static_cast<float>(countInclusive);
        float weightLtFirstCut = 1. / static_cast<float>(countLtFirstCut);
        float weightLtSecondCut = 1. / static_cast<float>(countLtSecondCut);

        switch (iBin) {
          case 0: {
            fillInvMass2D<true, false, 0>(recK0S, listrecPhi, multiplicity, weightInclusive, weightLtFirstCut, weightLtSecondCut);
            break;
          }
          case 1: {
            fillInvMass2D<true, false, 1>(recK0S, listrecPhi, multiplicity, weightInclusive, weightLtFirstCut, weightLtSecondCut);
            break;
          }
          case 2: {
            fillInvMass2D<true, false, 2>(recK0S, listrecPhi, multiplicity, weightInclusive, weightLtFirstCut, weightLtSecondCut);
            break;
          }
          case 3: {
            fillInvMass2D<true, false, 3>(recK0S, listrecPhi, multiplicity, weightInclusive, weightLtFirstCut, weightLtSecondCut);
            break;
          }
          case 4: {
            fillInvMass2D<true, false, 4>(recK0S, listrecPhi, multiplicity, weightInclusive, weightLtFirstCut, weightLtSecondCut);
            break;
          }
          case 5: {
            fillInvMass2D<true, false, 5>(recK0S, listrecPhi, multiplicity, weightInclusive, weightLtFirstCut, weightLtSecondCut);
            break;
          }
          case 6: {
            fillInvMass2D<true, false, 6>(recK0S, listrecPhi, multiplicity, weightInclusive, weightLtFirstCut, weightLtSecondCut);
            break;
          }
          case 7: {
            fillInvMass2D<true, false, 7>(recK0S, listrecPhi, multiplicity, weightInclusive, weightLtFirstCut, weightLtSecondCut);
            break;
          }
          case 8: {
            fillInvMass2D<true, false, 8>(recK0S, listrecPhi, multiplicity, weightInclusive, weightLtFirstCut, weightLtSecondCut);
            break;
          }
          case 9: {
            fillInvMass2D<true, false, 9>(recK0S, listrecPhi, multiplicity, weightInclusive, weightLtFirstCut, weightLtSecondCut);
            break;
          }
          default:
            break;
        }
      }
    }
  }

  PROCESS_SWITCH(phik0shortanalysis, processMEPhiK0S, "Process Mixed Event for Phi-K0S Analysis", false);

  void processMEPhiPion(soa::Filtered<SelCollisions> const& collisions, FullTracks const& fullTracks)
  {
    // Mixing the events with similar vertex z and multiplicity
    BinningTypeVertexContributor binningOnPositions{{axisVertex, axisMultiplicity}, true};
    for (auto const& [collision1, collision2] : o2::soa::selfCombinations(binningOnPositions, cfgNoMixedEvents, -1, collisions, collisions)) {
      if (!collision1.isInelGt0() || !collision2.isInelGt0())
        continue;

      float multiplicity = collision1.centFT0M();

      int iBin = 0;
      for (int i = 0; i < nMultBin; i++) {
        if (multBin[i] < multiplicity && multiplicity <= multBin[i + 1]) {
          iBin = i;
          break;
        }
      }

      // Defining V0s from collision1
      auto trackThisColl = fullTracks.sliceByCached(aod::track::collisionId, collision1.globalIndex(), cache);

      // Defining positive and negative tracks for phi reconstruction from collision1 and collision2, respectively
      auto posThisColl = posTracks->sliceByCached(aod::track::collisionId, collision1.globalIndex(), cache);
      auto negThisColl = negTracks->sliceByCached(aod::track::collisionId, collision2.globalIndex(), cache);

      for (const auto& track : trackThisColl) {

        if (!selectionPion(track))
          continue;

        TLorentzVector recPi;
        recPi.SetXYZM(track.px(), track.py(), track.pz(), massPi);
        if (recPi.Rapidity() > 0.8)
          continue;

        std::vector<TLorentzVector> listrecPhi;
        int countInclusive = 0, countLtFirstCut = 0, countLtSecondCut = 0;

        // Combinatorial background simulation
        for (auto& [track1, track2] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(posThisColl, negThisColl))) {
          if (!selectionTrackResonance(track1) || !selectionPIDKaon(track1) || !selectionTrackResonance(track2) || !selectionPIDKaon(track2))
            continue; // topological and PID selection

          TLorentzVector recPhi;
          recPhi = recMother(track1, track2, massKa, massKa);
          if (recPhi.Rapidity() > 0.8)
            continue;

          listrecPhi.push_back(recPhi);

          countInclusive++;
          if (std::abs(recPi.Rapidity() - recPhi.Rapidity()) > cfgFirstCutonDeltay)
            continue;
          countLtFirstCut++;
          if (std::abs(recPi.Rapidity() - recPhi.Rapidity()) > cfgSecondCutonDeltay)
            continue;
          countLtSecondCut++;
        }

        float nsigmaTPC, nsigmaTOF;
        if (track.hasTPC())
          nsigmaTPC = track.tpcNSigmaPi();
        else
          nsigmaTPC = -9.99;
        if (track.hasTOF())
          nsigmaTOF = track.tofNSigmaPi();
        else
          nsigmaTOF = -9.99;

        float weightInclusive = 1. / static_cast<float>(countInclusive);
        float weightLtFirstCut = 1. / static_cast<float>(countLtFirstCut);
        float weightLtSecondCut = 1. / static_cast<float>(countLtSecondCut);

        switch (iBin) {
          case 0: {
            fillInvMassNSigma<true, false, 0>(recPi, nsigmaTPC, nsigmaTOF, listrecPhi, multiplicity, weightInclusive, weightLtFirstCut, weightLtSecondCut);
            break;
          }
          case 1: {
            fillInvMassNSigma<true, false, 1>(recPi, nsigmaTPC, nsigmaTOF, listrecPhi, multiplicity, weightInclusive, weightLtFirstCut, weightLtSecondCut);
            break;
          }
          case 2: {
            fillInvMassNSigma<true, false, 2>(recPi, nsigmaTPC, nsigmaTOF, listrecPhi, multiplicity, weightInclusive, weightLtFirstCut, weightLtSecondCut);
            break;
          }
          case 3: {
            fillInvMassNSigma<true, false, 3>(recPi, nsigmaTPC, nsigmaTOF, listrecPhi, multiplicity, weightInclusive, weightLtFirstCut, weightLtSecondCut);
            break;
          }
          case 4: {
            fillInvMassNSigma<true, false, 4>(recPi, nsigmaTPC, nsigmaTOF, listrecPhi, multiplicity, weightInclusive, weightLtFirstCut, weightLtSecondCut);
            break;
          }
          case 5: {
            fillInvMassNSigma<true, false, 5>(recPi, nsigmaTPC, nsigmaTOF, listrecPhi, multiplicity, weightInclusive, weightLtFirstCut, weightLtSecondCut);
            break;
          }
          case 6: {
            fillInvMassNSigma<true, false, 6>(recPi, nsigmaTPC, nsigmaTOF, listrecPhi, multiplicity, weightInclusive, weightLtFirstCut, weightLtSecondCut);
            break;
          }
          case 7: {
            fillInvMassNSigma<true, false, 7>(recPi, nsigmaTPC, nsigmaTOF, listrecPhi, multiplicity, weightInclusive, weightLtFirstCut, weightLtSecondCut);
            break;
          }
          case 8: {
            fillInvMassNSigma<true, false, 8>(recPi, nsigmaTPC, nsigmaTOF, listrecPhi, multiplicity, weightInclusive, weightLtFirstCut, weightLtSecondCut);
            break;
          }
          case 9: {
            fillInvMassNSigma<true, false, 9>(recPi, nsigmaTPC, nsigmaTOF, listrecPhi, multiplicity, weightInclusive, weightLtFirstCut, weightLtSecondCut);
            break;
          }
          default:
            break;
        }
      }
    }
  }

  PROCESS_SWITCH(phik0shortanalysis, processMEPhiPion, "Process Mixed Event for Phi-Pion Analysis", false);

  void processRecMCPhiQA(SimCollisions::iterator const& collision, FullMCTracks const& fullMCTracks, FullV0s const& V0s, V0DauMCTracks const&, MCCollisions const&, aod::McParticles const&)
  {
    if (!acceptEventQA<true>(collision, true))
      return;

    float multiplicity = collision.centFT0M();
    MCeventHist.fill(HIST("hRecMCMultiplicityPercent"), multiplicity);

    // Defining positive and negative tracks for phi reconstruction
    auto posThisColl = posMCTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto negThisColl = negMCTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);

    bool isCountedPhi = false;

    for (auto track1 : posThisColl) { // loop over all selected tracks
      if (!selectionTrackResonance(track1) || !selectionPIDKaon(track1))
        continue; // topological and PID selection

      auto track1ID = track1.globalIndex();

      if (!track1.has_mcParticle())
        continue;

      // Loop over all negative candidates
      for (auto track2 : negThisColl) {
        if (!selectionTrackResonance(track2) || !selectionPIDKaon(track2))
          continue; // topological and PID selection

        auto track2ID = track2.globalIndex();
        if (track2ID == track1ID)
          continue; // condition to avoid double counting of pair

        if (!track2.has_mcParticle())
          continue;

        auto MCtrack1 = track1.mcParticle_as<aod::McParticles>();
        auto MCtrack2 = track2.mcParticle_as<aod::McParticles>();
        if (MCtrack1.pdgCode() != 321 || MCtrack2.pdgCode() != -321)
          continue;
        if (!MCtrack1.has_mothers() || !MCtrack2.has_mothers())
          continue;
        if (!MCtrack1.isPhysicalPrimary() || !MCtrack2.isPhysicalPrimary())
          continue;

        int pdgParentPhi = 0;
        for (const auto& MotherOfMCtrack1 : MCtrack1.mothers_as<aod::McParticles>()) {
          for (const auto& MotherOfMCtrack2 : MCtrack2.mothers_as<aod::McParticles>()) {
            if (MotherOfMCtrack1 == MotherOfMCtrack2) {
              pdgParentPhi = MotherOfMCtrack1.pdgCode();
            }
          }
        }

        if (pdgParentPhi != 333)
          continue;

        TLorentzVector recPhi;
        recPhi = recMother(track1, track2, massKa, massKa);
        if (recPhi.Rapidity() > 0.8)
          continue;

        if (!isCountedPhi) {
          MCeventHist.fill(HIST("hRecMCEventSelection"), 6); // at least a Phi in the event
          isCountedPhi = true;
        }

        PhieffHist.fill(HIST("h2PhieffInvMass"), multiplicity, recPhi.M());

        bool isCountedK0SInclusive = false, isCountedK0SFirstCut = false, isCountedK0SSecondCut = false;

        // V0 already reconstructed by the builder
        for (const auto& v0 : V0s) {
          const auto& posDaughterTrack = v0.posTrack_as<V0DauMCTracks>();
          const auto& negDaughterTrack = v0.negTrack_as<V0DauMCTracks>();
          if (!posDaughterTrack.has_mcParticle() || !negDaughterTrack.has_mcParticle())
            continue;

          auto posMCDaughterTrack = posDaughterTrack.mcParticle_as<aod::McParticles>();
          auto negMCDaughterTrack = negDaughterTrack.mcParticle_as<aod::McParticles>();
          if (posMCDaughterTrack.pdgCode() != 211 || negMCDaughterTrack.pdgCode() != -211)
            continue;
          if (!posMCDaughterTrack.has_mothers() || !negMCDaughterTrack.has_mothers())
            continue;

          int pdgParentv0 = 0;
          bool isPhysPrim = false;
          for (const auto& particleMotherOfNeg : negMCDaughterTrack.mothers_as<aod::McParticles>()) {
            for (const auto& particleMotherOfPos : posMCDaughterTrack.mothers_as<aod::McParticles>()) {
              if (particleMotherOfNeg == particleMotherOfPos) {
                pdgParentv0 = particleMotherOfNeg.pdgCode();
                isPhysPrim = particleMotherOfNeg.isPhysicalPrimary();
              }
            }
          }
          if (pdgParentv0 != 310 || !isPhysPrim)
            continue;

          // Cut on V0 dynamic columns
          if (!selectionV0(v0, posDaughterTrack, negDaughterTrack))
            continue;

          TLorentzVector recK0S;
          recK0S.SetXYZM(v0.px(), v0.py(), v0.pz(), v0.mK0Short());

          if (recK0S.Rapidity() > 0.8)
            continue;
          if (!isCountedK0SInclusive) {
            PhieffHist.fill(HIST("h2PhieffK0SInvMassInclusive"), multiplicity, recPhi.M());
            isCountedK0SInclusive = true;
          }
          if (std::abs(recK0S.Rapidity() - recPhi.Rapidity()) > cfgFirstCutonDeltay)
            continue;
          if (!isCountedK0SFirstCut) {
            PhieffHist.fill(HIST("h2PhieffK0SInvMassFirstCut"), multiplicity, recPhi.M());
            isCountedK0SFirstCut = true;
          }
          if (std::abs(recK0S.Rapidity() - recPhi.Rapidity()) > cfgSecondCutonDeltay)
            continue;
          if (!isCountedK0SSecondCut) {
            PhieffHist.fill(HIST("h2PhieffK0SInvMassSecondCut"), multiplicity, recPhi.M());
            isCountedK0SSecondCut = true;
          }
        }

        bool isCountedPiInclusive = false, isCountedPiFirstCut = false, isCountedPiSecondCut = false;

        // Loop over all primary pion candidates
        for (const auto& track : fullMCTracks) {
          if (!track.has_mcParticle())
            continue;

          auto MCtrack = track.mcParticle_as<aod::McParticles>();
          if (std::abs(MCtrack.pdgCode()) != 211 || !MCtrack.isPhysicalPrimary())
            continue;

          if (!selectionPion(track))
            continue;

          TLorentzVector recPi;
          recPi.SetXYZM(track.px(), track.py(), track.pz(), massPi);

          if (recPi.Rapidity() > 0.8)
            continue;
          if (!isCountedPiInclusive) {
            PhieffHist.fill(HIST("h3PhieffPiInvMassInclusive"), multiplicity, recPi.Pt(), recPhi.M());
            isCountedPiInclusive = true;
          }
          if (std::abs(recPi.Rapidity() - recPhi.Rapidity()) > cfgFirstCutonDeltay)
            continue;
          if (!isCountedPiFirstCut) {
            PhieffHist.fill(HIST("h3PhieffPiInvMassFirstCut"), multiplicity, recPi.Pt(), recPhi.M());
            isCountedPiFirstCut = true;
          }
          if (std::abs(recPi.Rapidity() - recPhi.Rapidity()) > cfgSecondCutonDeltay)
            continue;
          if (!isCountedPiSecondCut) {
            PhieffHist.fill(HIST("h3PhieffPiInvMassSecondCut"), multiplicity, recPi.Pt(), recPhi.M());
            isCountedPiSecondCut = true;
          }
        }
      }
    }
  }

  PROCESS_SWITCH(phik0shortanalysis, processRecMCPhiQA, "Process for ReCMCQA and Phi in RecMC", false);

  void processRecMCPhiK0S(SimCollisions::iterator const& collision, FullMCTracks const&, FullV0s const& V0s, V0DauMCTracks const&, MCCollisions const&, aod::McParticles const&)
  {
    if (!acceptEventQA<true>(collision, false))
      return;

    float multiplicity = collision.centFT0M();
    MCeventHist.fill(HIST("hRecMCMultiplicityPercent"), multiplicity);

    int iBin = 0;
    for (int i = 0; i < nMultBin; i++) {
      if (multBin[i] < multiplicity && multiplicity <= multBin[i + 1]) {
        iBin = i;
        break;
      }
    }

    // Defining positive and negative tracks for phi reconstruction
    auto posThisColl = posMCTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto negThisColl = negMCTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);

    // V0 already reconstructed by the builder
    for (const auto& v0 : V0s) {
      const auto& posDaughterTrack = v0.posTrack_as<V0DauMCTracks>();
      const auto& negDaughterTrack = v0.negTrack_as<V0DauMCTracks>();
      if (!posDaughterTrack.has_mcParticle() || !negDaughterTrack.has_mcParticle())
        continue;

      auto posMCDaughterTrack = posDaughterTrack.mcParticle_as<aod::McParticles>();
      auto negMCDaughterTrack = negDaughterTrack.mcParticle_as<aod::McParticles>();
      if (posMCDaughterTrack.pdgCode() != 211 || negMCDaughterTrack.pdgCode() != -211)
        continue;
      if (!posMCDaughterTrack.has_mothers() || !negMCDaughterTrack.has_mothers())
        continue;

      int pdgParentv0 = 0;
      bool isPhysPrim = false;
      for (const auto& particleMotherOfNeg : negMCDaughterTrack.mothers_as<aod::McParticles>()) {
        for (const auto& particleMotherOfPos : posMCDaughterTrack.mothers_as<aod::McParticles>()) {
          if (particleMotherOfNeg == particleMotherOfPos) {
            pdgParentv0 = particleMotherOfNeg.pdgCode();
            isPhysPrim = particleMotherOfNeg.isPhysicalPrimary();
          }
        }
      }
      if (pdgParentv0 != 310 || !isPhysPrim)
        continue;

      if (!selectionV0(v0, posDaughterTrack, negDaughterTrack))
        continue;

      TLorentzVector recK0S;
      recK0S.SetXYZM(v0.px(), v0.py(), v0.pz(), v0.mK0Short());
      if (recK0S.Rapidity() > 0.8)
        continue;

      std::vector<TLorentzVector> listrecPhi;
      int countInclusive = 0, countLtFirstCut = 0, countLtSecondCut = 0;

      // Phi reconstruction
      for (auto track1 : posThisColl) { // loop over all selected tracks
        if (!selectionTrackResonance(track1) || !selectionPIDKaon(track1))
          continue; // topological and PID selection

        auto track1ID = track1.globalIndex();

        if (!track1.has_mcParticle())
          continue;

        for (auto track2 : negThisColl) {
          if (!selectionTrackResonance(track2) || !selectionPIDKaon(track2))
            continue; // topological and PID selection

          auto track2ID = track2.globalIndex();
          if (track2ID == track1ID)
            continue; // condition to avoid double counting of pair

          if (!track2.has_mcParticle())
            continue;

          auto MCtrack1 = track1.mcParticle_as<aod::McParticles>();
          auto MCtrack2 = track2.mcParticle_as<aod::McParticles>();
          if (MCtrack1.pdgCode() != 321 || MCtrack2.pdgCode() != -321)
            continue;
          if (!MCtrack1.has_mothers() || !MCtrack2.has_mothers())
            continue;
          if (!MCtrack1.isPhysicalPrimary() || !MCtrack2.isPhysicalPrimary())
            continue;

          int pdgParentPhi = 0;
          for (const auto& MotherOfMCtrack1 : MCtrack1.mothers_as<aod::McParticles>()) {
            for (const auto& MotherOfMCtrack2 : MCtrack2.mothers_as<aod::McParticles>()) {
              if (MotherOfMCtrack1 == MotherOfMCtrack2) {
                pdgParentPhi = MotherOfMCtrack1.pdgCode();
              }
            }
          }

          if (pdgParentPhi != 333)
            continue;

          TLorentzVector recPhi;
          recPhi = recMother(track1, track2, massKa, massKa);
          if (recPhi.Rapidity() > 0.8)
            continue;

          listrecPhi.push_back(recPhi);

          countInclusive++;
          if (std::abs(recK0S.Rapidity() - recPhi.Rapidity()) > cfgFirstCutonDeltay)
            continue;
          countLtFirstCut++;
          if (std::abs(recK0S.Rapidity() - recPhi.Rapidity()) > cfgSecondCutonDeltay)
            continue;
          countLtSecondCut++;
        }
      }

      float weightInclusive, weightLtFirstCut, weightLtSecondCut;
      if (countInclusive > 0) {
        weightInclusive = 1. / static_cast<float>(countInclusive);
      } else {
        weightInclusive = 0;
        MCeventHist.fill(HIST("thereisnoPhiwK0SMC"), 0);
      }
      if (countLtFirstCut > 0) {
        weightLtFirstCut = 1. / static_cast<float>(countLtFirstCut);
      } else {
        weightLtFirstCut = 0;
        MCeventHist.fill(HIST("thereisnoPhiwK0SMC"), 1);
      }
      if (countLtSecondCut > 0) {
        weightLtSecondCut = 1. / static_cast<float>(countLtSecondCut);
      } else {
        weightLtSecondCut = 0;
        MCeventHist.fill(HIST("thereisnoPhiwK0SMC"), 2);
      }

      switch (iBin) {
        case 0: {
          fillInvMass2D<false, true, 0>(recK0S, listrecPhi, multiplicity, weightInclusive, weightLtFirstCut, weightLtSecondCut);
          break;
        }
        case 1: {
          fillInvMass2D<false, true, 1>(recK0S, listrecPhi, multiplicity, weightInclusive, weightLtFirstCut, weightLtSecondCut);
          break;
        }
        case 2: {
          fillInvMass2D<false, true, 2>(recK0S, listrecPhi, multiplicity, weightInclusive, weightLtFirstCut, weightLtSecondCut);
          break;
        }
        case 3: {
          fillInvMass2D<false, true, 3>(recK0S, listrecPhi, multiplicity, weightInclusive, weightLtFirstCut, weightLtSecondCut);
          break;
        }
        case 4: {
          fillInvMass2D<false, true, 4>(recK0S, listrecPhi, multiplicity, weightInclusive, weightLtFirstCut, weightLtSecondCut);
          break;
        }
        case 5: {
          fillInvMass2D<false, true, 5>(recK0S, listrecPhi, multiplicity, weightInclusive, weightLtFirstCut, weightLtSecondCut);
          break;
        }
        case 6: {
          fillInvMass2D<false, true, 6>(recK0S, listrecPhi, multiplicity, weightInclusive, weightLtFirstCut, weightLtSecondCut);
          break;
        }
        case 7: {
          fillInvMass2D<false, true, 7>(recK0S, listrecPhi, multiplicity, weightInclusive, weightLtFirstCut, weightLtSecondCut);
          break;
        }
        case 8: {
          fillInvMass2D<false, true, 8>(recK0S, listrecPhi, multiplicity, weightInclusive, weightLtFirstCut, weightLtSecondCut);
          break;
        }
        case 9: {
          fillInvMass2D<false, true, 9>(recK0S, listrecPhi, multiplicity, weightInclusive, weightLtFirstCut, weightLtSecondCut);
          break;
        }
        default:
          break;
      }
    }
  }

  PROCESS_SWITCH(phik0shortanalysis, processRecMCPhiK0S, "Process RecMC for Phi-K0S Analysis", false);

  void processRecMCPhiPion(SimCollisions::iterator const& collision, FullMCTracks const& fullMCTracks, MCCollisions const&, aod::McParticles const&)
  {
    if (!acceptEventQA<true>(collision, false))
      return;

    float multiplicity = collision.centFT0M();
    MCeventHist.fill(HIST("hRecMCMultiplicityPercent"), multiplicity);

    int iBin = 0;
    for (int i = 0; i < nMultBin; i++) {
      if (multBin[i] < multiplicity && multiplicity <= multBin[i + 1]) {
        iBin = i;
        break;
      }
    }

    // Defining positive and negative tracks for phi reconstruction
    auto posThisColl = posMCTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto negThisColl = negMCTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);

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

      TLorentzVector recPi;
      recPi.SetXYZM(track.px(), track.py(), track.pz(), massPi);
      if (recPi.Rapidity() > 0.8)
        continue;

      std::vector<TLorentzVector> listrecPhi;
      int countInclusive = 0, countLtFirstCut = 0, countLtSecondCut = 0;

      // Phi reconstruction
      for (auto track1 : posThisColl) { // loop over all selected tracks
        if (!selectionTrackResonance(track1) || !selectionPIDKaon(track1))
          continue; // topological and PID selection

        auto track1ID = track1.globalIndex();

        if (!track1.has_mcParticle())
          continue;

        for (auto track2 : negThisColl) {
          if (!selectionTrackResonance(track2) || !selectionPIDKaon(track2))
            continue; // topological and PID selection

          auto track2ID = track2.globalIndex();
          if (track2ID == track1ID)
            continue; // condition to avoid double counting of pair

          if (!track2.has_mcParticle())
            continue;

          auto MCtrack1 = track1.mcParticle_as<aod::McParticles>();
          auto MCtrack2 = track2.mcParticle_as<aod::McParticles>();
          if (MCtrack1.pdgCode() != 321 || MCtrack2.pdgCode() != -321)
            continue;
          if (!MCtrack1.has_mothers() || !MCtrack2.has_mothers())
            continue;
          if (!MCtrack1.isPhysicalPrimary() || !MCtrack2.isPhysicalPrimary())
            continue;

          int pdgParentPhi = 0;
          for (const auto& MotherOfMCtrack1 : MCtrack1.mothers_as<aod::McParticles>()) {
            for (const auto& MotherOfMCtrack2 : MCtrack2.mothers_as<aod::McParticles>()) {
              if (MotherOfMCtrack1 == MotherOfMCtrack2) {
                pdgParentPhi = MotherOfMCtrack1.pdgCode();
              }
            }
          }

          if (pdgParentPhi != 333)
            continue;

          TLorentzVector recPhi;
          recPhi = recMother(track1, track2, massKa, massKa);
          if (recPhi.Rapidity() > 0.8)
            continue;

          listrecPhi.push_back(recPhi);

          countInclusive++;
          if (std::abs(recPi.Rapidity() - recPhi.Rapidity()) > cfgFirstCutonDeltay)
            continue;
          countLtFirstCut++;
          if (std::abs(recPi.Rapidity() - recPhi.Rapidity()) > cfgSecondCutonDeltay)
            continue;
          countLtSecondCut++;
        }
      }

      float nsigmaTPC, nsigmaTOF;
      if (track.hasTPC())
        nsigmaTPC = track.tpcNSigmaPi();
      else
        nsigmaTPC = -9.99;
      if (track.hasTOF())
        nsigmaTOF = track.tofNSigmaPi();
      else
        nsigmaTOF = -9.99;

      float weightInclusive, weightLtFirstCut, weightLtSecondCut;
      if (countInclusive > 0) {
        weightInclusive = 1. / static_cast<float>(countInclusive);
      } else {
        weightInclusive = 0;
        MCeventHist.fill(HIST("thereisnoPhiwPiMC"), 0);
      }
      if (countLtFirstCut > 0) {
        weightLtFirstCut = 1. / static_cast<float>(countLtFirstCut);
      } else {
        weightLtFirstCut = 0;
        MCeventHist.fill(HIST("thereisnoPhiwPiMC"), 1);
      }
      if (countLtSecondCut > 0) {
        weightLtSecondCut = 1. / static_cast<float>(countLtSecondCut);
      } else {
        weightLtSecondCut = 0;
        MCeventHist.fill(HIST("thereisnoPhiwPiMC"), 2);
      }

      switch (iBin) {
        case 0: {
          fillInvMassNSigma<false, true, 0>(recPi, nsigmaTPC, nsigmaTOF, listrecPhi, multiplicity, weightInclusive, weightLtFirstCut, weightLtSecondCut);
          break;
        }
        case 1: {
          fillInvMassNSigma<false, true, 1>(recPi, nsigmaTPC, nsigmaTOF, listrecPhi, multiplicity, weightInclusive, weightLtFirstCut, weightLtSecondCut);
          break;
        }
        case 2: {
          fillInvMassNSigma<false, true, 2>(recPi, nsigmaTPC, nsigmaTOF, listrecPhi, multiplicity, weightInclusive, weightLtFirstCut, weightLtSecondCut);
          break;
        }
        case 3: {
          fillInvMassNSigma<false, true, 3>(recPi, nsigmaTPC, nsigmaTOF, listrecPhi, multiplicity, weightInclusive, weightLtFirstCut, weightLtSecondCut);
          break;
        }
        case 4: {
          fillInvMassNSigma<false, true, 4>(recPi, nsigmaTPC, nsigmaTOF, listrecPhi, multiplicity, weightInclusive, weightLtFirstCut, weightLtSecondCut);
          break;
        }
        case 5: {
          fillInvMassNSigma<false, true, 5>(recPi, nsigmaTPC, nsigmaTOF, listrecPhi, multiplicity, weightInclusive, weightLtFirstCut, weightLtSecondCut);
          break;
        }
        case 6: {
          fillInvMassNSigma<false, true, 6>(recPi, nsigmaTPC, nsigmaTOF, listrecPhi, multiplicity, weightInclusive, weightLtFirstCut, weightLtSecondCut);
          break;
        }
        case 7: {
          fillInvMassNSigma<false, true, 7>(recPi, nsigmaTPC, nsigmaTOF, listrecPhi, multiplicity, weightInclusive, weightLtFirstCut, weightLtSecondCut);
          break;
        }
        case 8: {
          fillInvMassNSigma<false, true, 8>(recPi, nsigmaTPC, nsigmaTOF, listrecPhi, multiplicity, weightInclusive, weightLtFirstCut, weightLtSecondCut);
          break;
        }
        case 9: {
          fillInvMassNSigma<false, true, 9>(recPi, nsigmaTPC, nsigmaTOF, listrecPhi, multiplicity, weightInclusive, weightLtFirstCut, weightLtSecondCut);
          break;
        }
        default:
          break;
      }
    }
  }

  PROCESS_SWITCH(phik0shortanalysis, processRecMCPhiPion, "Process RecMC for Phi-Pion Analysis", false);

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

    if (collisions.size() < 1)
      return;

    bool isAssocColl = false;
    for (auto collision : collisions) {
      if (acceptEventQA<true>(collision, false)) {
        isAssocColl = true;
        break;
      }
    }
    if (!isAssocColl)
      return;
    MCeventHist.fill(HIST("hGenMCEventSelection"), 3); // with at least a reco collision

    float multiplicity = mcCollision.centFT0M();
    MCeventHist.fill(HIST("hGenMCMultiplicityPercent"), multiplicity);

    int imultBin = 0;
    for (int i = 0; i < nMultBin; i++) {
      if (multBin[i] < multiplicity && multiplicity <= multBin[i + 1]) {
        imultBin = i;
        break;
      }
    }

    bool isCountedPhi = false;

    for (auto mcParticle1 : mcParticles) {
      if (mcParticle1.pdgCode() != 333)
        continue;
      auto kDaughters = mcParticle1.daughters_as<aod::McParticles>();
      if (kDaughters.size() != 2)
        continue;
      bool isPosKaon = false, isNegKaon = false;
      for (auto kDaughter : kDaughters) {
        if (kDaughter.pdgCode() == 321)
          isPosKaon = true;
        if (kDaughter.pdgCode() == -321)
          isNegKaon = true;
      }
      if (!isPosKaon || !isNegKaon)
        continue;
      if (mcParticle1.y() > 0.8)
        continue;

      if (!isCountedPhi) {
        MCeventHist.fill(HIST("hGenMCEventSelection"), 4); // at least a Phi in the event
        isCountedPhi = true;
      }

      bool isCountedK0SInclusive = false, isCountedK0SFirstCut = false, isCountedK0SSecondCut = false;

      for (auto mcParticle2 : mcParticles) {
        if (mcParticle2.pdgCode() != 310)
          continue;
        if (!mcParticle2.isPhysicalPrimary())
          continue;
        if (mcParticle2.y() > 0.8)
          continue;

        if (!isCountedK0SInclusive) {
          PhieffHist.fill(HIST("h1PhieffK0SGenMCInclusive"), imultBin);
          isCountedK0SInclusive = true;
        }
        if (std::abs(mcParticle1.y() - mcParticle2.y()) > cfgFirstCutonDeltay)
          continue;
        if (!isCountedK0SFirstCut) {
          PhieffHist.fill(HIST("h1PhieffK0SGenMCFirstCut"), imultBin);
          isCountedK0SFirstCut = true;
        }
        if (std::abs(mcParticle1.y() - mcParticle2.y()) > cfgSecondCutonDeltay)
          continue;
        if (!isCountedK0SSecondCut) {
          PhieffHist.fill(HIST("h1PhieffK0SGenMCSecondCut"), imultBin);
          isCountedK0SSecondCut = true;
        }
      }

      bool isCountedPiInclusive = false, isCountedPiFirstCut = false, isCountedPiSecondCut = false;

      for (auto mcParticle2 : mcParticles) {
        if (std::abs(mcParticle2.pdgCode()) != 211)
          continue;
        if (!mcParticle2.isPhysicalPrimary())
          continue;
        if (mcParticle2.y() > 0.8)
          continue;

        int ipTBin = 0;
        for (int i = 0; i < nPtBin; i++) {
          if (pTBin[i] < mcParticle1.pt() && mcParticle1.pt() <= pTBin[i + 1]) {
            ipTBin = i;
            break;
          }
        }

        if (!isCountedPiInclusive) {
          PhieffHist.fill(HIST("h2PhieffPiGenMCInclusive"), imultBin, ipTBin);
          isCountedPiInclusive = true;
        }
        if (std::abs(mcParticle1.y() - mcParticle2.y()) > cfgFirstCutonDeltay)
          continue;
        if (!isCountedPiFirstCut) {
          PhieffHist.fill(HIST("h2PhieffPiGenMCFirstCut"), imultBin, ipTBin);
          isCountedPiFirstCut = true;
        }
        if (std::abs(mcParticle1.y() - mcParticle2.y()) > cfgSecondCutonDeltay)
          continue;
        if (!isCountedPiSecondCut) {
          PhieffHist.fill(HIST("h2PhieffPiGenMCSecondCut"), imultBin, ipTBin);
          isCountedPiSecondCut = true;
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
    for (auto collision : collisions) {
      if (acceptEventQA<true>(collision, false)) {
        isAssocColl = true;
        break;
      }
    }

    float multiplicity = mcCollision.centFT0M();
    MCeventHist.fill(HIST("hGenMCMultiplicityPercent"), multiplicity);

    int iBin = 0;
    for (int i = 0; i < nMultBin; i++) {
      if (multBin[i] < multiplicity && multiplicity <= multBin[i + 1]) {
        iBin = i;
        break;
      }
    }

    for (auto mcParticle1 : mcParticles) {
      if (mcParticle1.pdgCode() != 310)
        continue;
      if (!mcParticle1.isPhysicalPrimary())
        continue;
      if (mcParticle1.y() > 0.8)
        continue;

      bool isCountedPhiInclusive = false, isCountedPhiFirstCut = false, isCountedPhiSecondCut = false;
      bool isCountedPhiInclusiveAssocReco = false, isCountedPhiFirstCutAssocReco = false, isCountedPhiSecondCutAssocReco = false;

      for (auto mcParticle2 : mcParticles) {
        if (mcParticle2.pdgCode() != 333)
          continue;
        auto kDaughters = mcParticle2.daughters_as<aod::McParticles>();
        if (kDaughters.size() != 2)
          continue;
        bool isPosKaon = false, isNegKaon = false;
        for (auto kDaughter : kDaughters) {
          if (kDaughter.pdgCode() == 321)
            isPosKaon = true;
          if (kDaughter.pdgCode() == -321)
            isNegKaon = true;
        }
        if (!isPosKaon || !isNegKaon)
          continue;
        if (mcParticle2.y() > 0.8)
          continue;

        if (isAssocColl) {
          if (!isCountedPhiInclusiveAssocReco) {
            MCPhiK0SHist.fill(HIST("h1PhiK0SGenMCInclusiveAssocReco"), iBin);
            isCountedPhiInclusiveAssocReco = true;
          }
          if (std::abs(mcParticle1.y() - mcParticle2.y()) > cfgFirstCutonDeltay)
            continue;
          if (!isCountedPhiFirstCutAssocReco) {
            MCPhiK0SHist.fill(HIST("h1PhiK0SGenMCFirstCutAssocReco"), iBin);
            isCountedPhiFirstCutAssocReco = true;
          }
          if (std::abs(mcParticle1.y() - mcParticle2.y()) > cfgSecondCutonDeltay)
            continue;
          if (!isCountedPhiSecondCutAssocReco) {
            MCPhiK0SHist.fill(HIST("h1PhiK0SGenMCSecondCutAssocReco"), iBin);
            isCountedPhiSecondCutAssocReco = true;
          }
        } else {
          if (!isCountedPhiInclusive) {
            MCPhiK0SHist.fill(HIST("h1PhiK0SGenMCInclusive"), iBin);
            isCountedPhiInclusive = true;
          }
          if (std::abs(mcParticle1.y() - mcParticle2.y()) > cfgFirstCutonDeltay)
            continue;
          if (!isCountedPhiFirstCut) {
            MCPhiK0SHist.fill(HIST("h1PhiK0SGenMCFirstCut"), iBin);
            isCountedPhiFirstCut = true;
          }
          if (std::abs(mcParticle1.y() - mcParticle2.y()) > cfgSecondCutonDeltay)
            continue;
          if (!isCountedPhiSecondCut) {
            MCPhiK0SHist.fill(HIST("h1PhiK0SGenMCSecondCut"), iBin);
            isCountedPhiSecondCut = true;
          }
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
    for (auto collision : collisions) {
      if (acceptEventQA<true>(collision, false)) {
        isAssocColl = true;
        break;
      }
    }

    float multiplicity = mcCollision.centFT0M();
    MCeventHist.fill(HIST("hGenMCMultiplicityPercent"), multiplicity);

    int imultBin = 0;
    for (int i = 0; i < nMultBin; i++) {
      if (multBin[i] < multiplicity && multiplicity <= multBin[i + 1]) {
        imultBin = i;
        break;
      }
    }

    for (auto mcParticle1 : mcParticles) {
      if (std::abs(mcParticle1.pdgCode()) != 211)
        continue;
      if (!mcParticle1.isPhysicalPrimary())
        continue;
      if (mcParticle1.y() > 0.8)
        continue;

      int ipTBin = 0;
      for (int i = 0; i < nPtBin; i++) {
        if (pTBin[i] < mcParticle1.pt() && mcParticle1.pt() <= pTBin[i + 1]) {
          ipTBin = i;
          break;
        }
      }

      bool isCountedPhiInclusive = false, isCountedPhiFirstCut = false, isCountedPhiSecondCut = false;
      bool isCountedPhiInclusiveAssocReco = false, isCountedPhiFirstCutAssocReco = false, isCountedPhiSecondCutAssocReco = false;

      for (auto mcParticle2 : mcParticles) {
        if (mcParticle2.pdgCode() != 333)
          continue;
        auto kDaughters = mcParticle2.daughters_as<aod::McParticles>();
        if (kDaughters.size() != 2)
          continue;
        bool isPosKaon = false, isNegKaon = false;
        for (auto kDaughter : kDaughters) {
          if (kDaughter.pdgCode() == 321)
            isPosKaon = true;
          if (kDaughter.pdgCode() == -321)
            isNegKaon = true;
        }
        if (!isPosKaon || !isNegKaon)
          continue;
        if (mcParticle2.y() > 0.8)
          continue;

        if (isAssocColl) {
          if (!isCountedPhiInclusiveAssocReco) {
            MCPhiPionHist.fill(HIST("h1PhiPiGenMCInclusiveAssocReco"), imultBin, ipTBin);
            isCountedPhiInclusiveAssocReco = true;
          }
          if (std::abs(mcParticle1.y() - mcParticle2.y()) > cfgFirstCutonDeltay)
            continue;
          if (!isCountedPhiFirstCutAssocReco) {
            MCPhiPionHist.fill(HIST("h1PhiPiGenMCFirstCutAssocReco"), imultBin, ipTBin);
            isCountedPhiFirstCutAssocReco = true;
          }
          if (std::abs(mcParticle1.y() - mcParticle2.y()) > cfgSecondCutonDeltay)
            continue;
          if (!isCountedPhiSecondCutAssocReco) {
            MCPhiPionHist.fill(HIST("h1PhiPiGenMCSecondCutAssocReco"), imultBin, ipTBin);
            isCountedPhiSecondCutAssocReco = true;
          }
        } else {
          if (!isCountedPhiInclusive) {
            MCPhiPionHist.fill(HIST("h1PhiPiGenMCInclusive"), imultBin, ipTBin);
            isCountedPhiInclusive = true;
          }
          if (std::abs(mcParticle1.y() - mcParticle2.y()) > cfgFirstCutonDeltay)
            continue;
          if (!isCountedPhiFirstCut) {
            MCPhiPionHist.fill(HIST("h1PhiPiGenMCFirstCut"), imultBin, ipTBin);
            isCountedPhiFirstCut = true;
          }
          if (std::abs(mcParticle1.y() - mcParticle2.y()) > cfgSecondCutonDeltay)
            continue;
          if (!isCountedPhiSecondCut) {
            MCPhiPionHist.fill(HIST("h1PhiPiGenMCSecondCut"), imultBin, ipTBin);
            isCountedPhiSecondCut = true;
          }
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
