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

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace
{
static constexpr int nMultBin = 10;

constexpr float mPhi[nMultBin] = {};

static constexpr float multBin[nMultBin + 1] = {0.0, 1.0, 5.0, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 70.0, 100.0};

static constexpr std::string_view PhiK0SSEInc[nMultBin] = {"h2PhiK0SSEInc_0_1", "h2PhiK0SSEInc_1_5", "h2PhiK0SSEInc_5_10", "h2PhiK0SSEInc_10_15", "h2PhiK0SSEInc_15_20",
                                                           "h2PhiK0SSEInc_20_30", "h2PhiK0SSEInc_30_40", "h2PhiK0SSEInc_40_50", "h2PhiK0SSEInc_50_70", "h2PhiK0SSEInc_70_100"};
static constexpr std::string_view PhiK0SSEFCut[nMultBin] = {"h2PhiK0SSEFCut_0_1", "h2PhiK0SSEFCut_1_5", "h2PhiK0SSEFCut_5_10", "h2PhiK0SSEFCut_10_15", "h2PhiK0SSEFCut_15_20",
                                                            "h2PhiK0SSEFCut_20_30", "h2PhiK0SSEFCut_30_40", "h2PhiK0SSEFCut_40_50", "h2PhiK0SSEFCut_50_70", "h2PhiK0SSEFCut_70_100"};
static constexpr std::string_view PhiK0SSESCut[nMultBin] = {"h2PhiK0SSESCut_0_1", "h2PhiK0SSESCut_1_5", "h2PhiK0SSESCut_5_10", "h2PhiK0SSESCut_10_15", "h2PhiK0SSESCut_15_20",
                                                            "h2PhiK0SSESCut_20_30", "h2PhiK0SSESCut_30_40", "h2PhiK0SSESCut_40_50", "h2PhiK0SSESCut_50_70", "h2PhiK0SSESCut_70_100"};
} // namespace

struct phik0shortanalysis {
  // Histograms are defined with HistogramRegistry
  HistogramRegistry eventHist{"eventHist", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry K0SHist{"K0SHist", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry PhicandHist{"PhicandHist", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry PhipurHist{"PhipurHist", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry PhiK0SHist{"PhiK0SHist", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry PhiPionHist{"PhiPionHist", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

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
  Configurable<int> nBins{"nBins", 50, "N bins in cfgPhimassaxis"};
  Configurable<std::vector<float>> lowmPhiInc{"lowmPhiInc", std::vector<float>{mPhi, mPhi + nMultBin}, "Lower limits on Phi mass Inclusive"};
  Configurable<std::vector<float>> upmPhiInc{"upmPhiInc", std::vector<float>{mPhi, mPhi + nMultBin}, "Upper limits on Phi mass Inclusive"};
  Configurable<std::vector<float>> lowmPhiFCut{"lowmPhiFCut", std::vector<float>{mPhi, mPhi + nMultBin}, "Lower limits on Phi mass First Cut"};
  Configurable<std::vector<float>> upmPhiFCut{"upmPhiFCut", std::vector<float>{mPhi, mPhi + nMultBin}, "Upper limits on Phi mass First Cut"};
  Configurable<std::vector<float>> lowmPhiSCut{"lowmPhiSCut", std::vector<float>{mPhi, mPhi + nMultBin}, "Lower limits on Phi mass Second Cut"};
  Configurable<std::vector<float>> upmPhiSCut{"upmPhiSCut", std::vector<float>{mPhi, mPhi + nMultBin}, "Upper limits on Phi mass Second Cut"};

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

  // Configurables for invariant mass histograms filling
  Configurable<double> cfgFirstCutonDeltay{"cgfFirstCutonDeltay", 0.5, "First upper bound on Deltay selection"};
  Configurable<double> cfgSecondCutonDeltay{"cgfSecondCutonDeltay", 0.2, "Second upper bound on Deltay selection"};

  // Configurable for event mixing
  Configurable<int> cfgNoMixedEvents{"cfgNoMixedEvents", 5, "Number of mixed events per event"};

  // Configurable axis
  ConfigurableAxis axisVertex{"axisVertex", {20, -10, 10}, "vertex axis for bin"};
  ConfigurableAxis axisMultiplicityClass{"axisMultiplicityClass", {20, 0, 100}, "multiplicity percentile for bin"};
  ConfigurableAxis axisMultiplicity{"axisMultiplicity", {2000, 0, 10000}, "TPC multiplicity  for bin"};

  void init(InitContext const&)
  {
    // Axes
    AxisSpec K0SmassAxis = {200, 0.45f, 0.55f, "#it{M}_{inv} [GeV/#it{c}^{2}]"};
    AxisSpec PhimassAxis = {200, 0.9f, 1.2f, "#it{M}_{inv} [GeV/#it{c}^{2}]"};
    AxisSpec vertexZAxis = {100, -15.f, 15.f, "vrtx_{Z} [cm]"};
    AxisSpec deltayAxis = {16, 0.0f, 0.8f, "|#it{#Deltay}|"};
    AxisSpec multAxis = {120, 0.0f, 120.0f, "centFT0M"};
    AxisSpec ptAxis = {100, 0.0f, 10.0f, "#it{p}_{T} (GeV/#it{c})"};
    std::vector<AxisSpec> cfgPhimassAxisInc;
    std::vector<AxisSpec> cfgPhimassAxisFCut;
    std::vector<AxisSpec> cfgPhimassAxisSCut;
    for (int i = 0; i < nMultBin; i++) {
      cfgPhimassAxisInc.push_back({nBins, lowmPhiInc->at(i), upmPhiInc->at(i), "#it{M}_{inv} [GeV/#it{c}^{2}]"});
      cfgPhimassAxisFCut.push_back({nBins, lowmPhiFCut->at(i), upmPhiFCut->at(i), "#it{M}_{inv} [GeV/#it{c}^{2}]"});
      cfgPhimassAxisSCut.push_back({nBins, lowmPhiSCut->at(i), upmPhiFCut->at(i), "#it{M}_{inv} [GeV/#it{c}^{2}]"});
    }

    // Histograms
    // Number of events per selection
    eventHist.add("hEventSelection", "hEVentSelection", kTH1F, {{5, -0.5f, 4.5f}});
    eventHist.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(1, "All collisions");
    eventHist.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(2, "sel8 cut");
    eventHist.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(3, "posZ cut");
    eventHist.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(4, "INEL>0 cut");
    eventHist.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(5, "With at least a #phi cand");

    // Event information
    eventHist.add("hVertexZRec", "hVertexZRec", kTH1F, {vertexZAxis});
    eventHist.add("hMultiplicityPercent", "Multiplicity Percentile", kTH1F, {multAxis});

    // K0S topological/PID cuts
    K0SHist.add("hDCAV0Daughters", "hDCAV0Daughters", kTH1F, {{55, 0.0f, 2.2f}});
    K0SHist.add("hV0CosPA", "hV0CosPA", kTH1F, {{100, 0.95f, 1.f}});
    K0SHist.add("hNSigmaPosPionFromK0S", "hNSigmaPosPionFromK0Short", kTH2F, {ptAxis, {100, -5.f, 5.f}});
    K0SHist.add("hNSigmaNegPionFromK0S", "hNSigmaNegPionFromK0Short", kTH2F, {ptAxis, {100, -5.f, 5.f}});

    // Phi tpological/PID cuts
    PhicandHist.add("hEta", "Eta distribution", kTH1F, {{200, -1.0f, 1.0f}});
    PhicandHist.add("hDcaxy", "Dcaxy distribution", kTH1F, {{200, -1.0f, 1.0f}});
    PhicandHist.add("hDcaz", "Dcaz distribution", kTH1F, {{200, -1.0f, 1.0f}});
    PhicandHist.add("hNsigmaKaonTPC", "NsigmaKaon TPC distribution", kTH2F, {ptAxis, {100, -10.0f, 10.0f}});
    PhicandHist.add("hNsigmaKaonTOF", "NsigmaKaon TOF distribution", kTH2F, {ptAxis, {100, -10.0f, 10.0f}});

    // Phi invariant mass for computing purities and normalisation
    PhipurHist.add("h2PhipurInvMass", "Invariant mass of Phi for Purity (no K0S/Pi)", kTH2F, {multAxis, PhimassAxis});
    PhipurHist.add("h2PhipurK0SInvMassInclusive", "Invariant mass of Phi for Purity (K0S) Inclusive", kTH2F, {multAxis, PhimassAxis});
    PhipurHist.add("h2PhipurK0SInvMassFirstCut", "Invariant mass of Phi for Purity (K0S) Deltay < FirstCut", kTH2F, {multAxis, PhimassAxis});
    PhipurHist.add("h2PhipurK0SInvMassSecondCut", "Invariant mass of Phi for Purity (K0S) Deltay < SecondCut", kTH2F, {multAxis, PhimassAxis});
    PhipurHist.add("h3PhipurPiInvMassInclusive", "Invariant mass of Phi for Purity (Pi) Inclusive", kTH3F, {multAxis, ptAxis, PhimassAxis});
    PhipurHist.add("h3PhipurPiInvMassFirstCut", "Invariant mass of Phi for Purity (Pi) Deltay < FirstCut", kTH3F, {multAxis, ptAxis, PhimassAxis});
    PhipurHist.add("h3PhipurPiInvMassSecondCut", "Invariant mass of Phi for Purity (Pi) Deltay < SecondCut", kTH3F, {multAxis, ptAxis, PhimassAxis});

    // 2D mass for Phi and K0S
    for (int i = 0; i < nMultBin; i++) {
      PhiK0SHist.add(PhiK0SSEInc[i].data(), "2D Invariant mass of Phi and K0Short for Same Event Inclusive", kTH2F, {K0SmassAxis, cfgPhimassAxisInc.at(i)});
      PhiK0SHist.add(PhiK0SSEFCut[i].data(), "2D Invariant mass of Phi and K0Short for Same Event Deltay < FirstCut", kTH2F, {K0SmassAxis, cfgPhimassAxisFCut.at(i)});
      PhiK0SHist.add(PhiK0SSESCut[i].data(), "2D Invariant mass of Phi and K0Short for Same Event Deltay < SecondCut", kTH2F, {K0SmassAxis, cfgPhimassAxisSCut.at(i)});
    }

    PhiK0SHist.add("h3PhiK0SInvMassMixedEventInclusive", "2D Invariant mass of Phi and K0Short for Mixed Event Inclusive", kTH3F, {multAxis, K0SmassAxis, PhimassAxis});
    PhiK0SHist.add("h3PhiK0SInvMassMixedEventFirstCut", "2D Invariant mass of Phi and K0Short for Mixed Event Deltay < FirstCut", kTH3F, {multAxis, K0SmassAxis, PhimassAxis});
    PhiK0SHist.add("h3PhiK0SInvMassMixedEventSecondCut", "2D Invariant mass of Phi and K0Short for Mixed Event Deltay < SecondCut", kTH3F, {multAxis, K0SmassAxis, PhimassAxis});

    // Phi mass vs Pion dE/dx
    PhiPionHist.add("h4PhiInvMassPidEdxSameEventInclusive", "Phi Invariant mass vs Pion dE/dx for Same Event Inclusive", kTHnSparseF, {multAxis, ptAxis, {100, -10.0f, 10.0f}, PhimassAxis});
    PhiPionHist.add("h4PhiInvMassPidEdxSameEventFirstCut", "Phi Invariant mass vs Pion dE/dx for Same Event Deltay < FirstCut", kTHnSparseF, {multAxis, ptAxis, {100, -10.0f, 10.0f}, PhimassAxis});
    PhiPionHist.add("h4PhiInvMassPidEdxSameEventSecondCut", "Phi Invariant mass vs Pion dE/dx for Same Event Deltay < SecondCut", kTHnSparseF, {multAxis, ptAxis, {100, -10.0f, 10.0f}, PhimassAxis});

    PhiPionHist.add("h4PhiInvMassPidEdxMixedEventInclusive", "Phi Invariant mass vs Pion dE/dx for Mixed Event Inclusive", kTHnSparseF, {multAxis, ptAxis, {100, -10.0f, 10.0f}, PhimassAxis});
    PhiPionHist.add("h4PhiInvMassPidEdxMixedEventFirstCut", "Phi Invariant mass vs Pion dE/dx for Mixed Event Deltay < FirstCut", kTHnSparseF, {multAxis, ptAxis, {100, -10.0f, 10.0f}, PhimassAxis});
    PhiPionHist.add("h4PhiInvMassPidEdxMixedEventSecondCut", "Phi Invariant mass vs Pion dE/dx for Mixed Event Deltay < SecondCut", kTHnSparseF, {multAxis, ptAxis, {100, -10.0f, 10.0f}, PhimassAxis});
  }

  // Constants
  double massKa = o2::constants::physics::MassKPlus;
  double massPi = o2::constants::physics::MassPiPlus;

  // Defining filters for events (event selection)
  // Processed events will be already fulfilling the event selection requirements
  Filter eventFilter = (o2::aod::evsel::sel8 == true);
  Filter posZFilter = (nabs(o2::aod::collision::posZ) < cutzvertex);

  // Defining filters on V0s (cannot filter on dynamic columns)
  Filter preFilterV0 = (nabs(aod::v0data::dcapostopv) > v0setting_dcapostopv && nabs(aod::v0data::dcanegtopv) > v0setting_dcanegtopv && aod::v0data::dcaV0daughters < v0setting_dcav0dau);

  // Defining the type of the event
  using EventCandidates = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::PVMults>;

  // Defining the type of the Phi daughter tracks
  using PhiDaughterCandidates = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullKa, aod::pidTOFFullKa>;

  // Defining the type of the V0s
  using V0Candidates = soa::Filtered<aod::V0Datas>;

  // Defining the type of the V0 daughter tracks
  using V0DaughterCandidates = soa::Join<aod::TracksIU, aod::TracksExtra, aod::pidTPCFullPi>;

  // Defining the type of the Pions
  using PionCandidates = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullPi, aod::pidTOFFullPi>;

  // Defining the binning policy for mixed event
  using BinningTypeVertexContributor = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0M>;

  SliceCache cache;
  Partition<PhiDaughterCandidates> posTracks = aod::track::signed1Pt > cfgCutCharge;
  Partition<PhiDaughterCandidates> negTracks = aod::track::signed1Pt < cfgCutCharge;

  // Event selection and QA filling
  template <typename T>
  bool acceptEventQA(const T& collision)
  {
    eventHist.fill(HIST("hEventSelection"), 0); // all collisions
    if (!collision.sel8())
      return false;
    eventHist.fill(HIST("hEventSelection"), 1); // sel8 collisions
    if (std::abs(collision.posZ()) > cutzvertex)
      return false;
    eventHist.fill(HIST("hEventSelection"), 2); // vertex-Z selected
    eventHist.fill(HIST("hVertexZRec"), collision.posZ());
    if (!collision.isInelGt0())
      return false;
    eventHist.fill(HIST("hEventSelection"), 3); // INEL>0 collisions
    return true;
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
    if (std::abs(track.dcaXY()) > dcaxyMax)
      return false;
    if (std::abs(track.dcaZ()) > dcazMax)
      return false;
    return true;
  }

  // Fill 2D invariant mass histogram for V0 and Phi
  template <bool isMix, std::size_t iBin>
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
    }
  }

  // Fill Phi invariant mass vs Pion dE/dx histogram
  template <bool isMix>
  void fillInvMassdEdx(TLorentzVector Pi, float dEdx, const std::vector<TLorentzVector> listPhi, float multiplicity, double weightInclusive, double weightLtFirstCut, double weightLtSecondCut)
  {
    double rapidityPi = Pi.Rapidity();
    double ptPi = Pi.Pt();

    for (unsigned int phitag = 0; phitag < listPhi.size(); phitag++) {
      double massPhi = listPhi[phitag].M();
      double rapidityPhi = listPhi[phitag].Rapidity();
      double deltay = std::abs(rapidityPi - rapidityPhi);

      if constexpr (!isMix) { // same event
        PhiK0SHist.fill(HIST("h4PhiInvMassPidEdxSameEventInclusive"), multiplicity, ptPi, dEdx, massPhi, weightInclusive);
        if (deltay > cfgFirstCutonDeltay)
          continue;
        PhiK0SHist.fill(HIST("h4PhiInvMassPidEdxSameEventFirstCut"), multiplicity, ptPi, dEdx, massPhi, weightLtFirstCut);
        if (deltay > cfgSecondCutonDeltay)
          continue;
        PhiK0SHist.fill(HIST("h4PhiInvMassPidEdxsSameEventSecondCut"), multiplicity, ptPi, dEdx, massPhi, weightLtSecondCut);
      } else { // mixed event
        PhiK0SHist.fill(HIST("h4PhiInvMassPidEdxMixedEventInclusive"), multiplicity, ptPi, dEdx, massPhi, weightInclusive);
        if (deltay > cfgFirstCutonDeltay)
          continue;
        PhiK0SHist.fill(HIST("h4PhiInvMassPidEdxMixedEventFirstCut"), multiplicity, ptPi, dEdx, massPhi, weightLtFirstCut);
        if (deltay > cfgSecondCutonDeltay)
          continue;
        PhiK0SHist.fill(HIST("h4PhiInvMassPidEdxMixedEventSecondCut"), multiplicity, ptPi, dEdx, massPhi, weightLtSecondCut);
      }
    }
  }

  void processQAPurity(EventCandidates::iterator const& collision, PhiDaughterCandidates const&, PionCandidates const& Pis, V0Candidates const& V0s, V0DaughterCandidates const&)
  {
    // Check if the event selection is passed
    if (!acceptEventQA(collision))
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

        bool isCountedK0SInclusive = false;
        bool isCountedK0SFirstCut = false;
        bool isCountedK0SSecondCut = false;

        // V0 already reconstructed by the builder
        for (const auto& v0 : V0s) {
          const auto& posDaughterTrack = v0.posTrack_as<V0DaughterCandidates>();
          const auto& negDaughterTrack = v0.negTrack_as<V0DaughterCandidates>();

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

        bool isCountedPiInclusive = false;
        bool isCountedPiFirstCut = false;
        bool isCountedPiSecondCut = false;

        // Loop over all primary pion candidates
        for (const auto& pi : Pis) {
          if (!selectionPion(pi))
            continue;

          TLorentzVector recPi;
          recPi.SetXYZM(pi.px(), pi.py(), pi.pz(), massPi);

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

  PROCESS_SWITCH(phik0shortanalysis, processQAPurity, "Process Same Event for QA and Phi Purities", true);

  void processSEPhiK0S(soa::Filtered<EventCandidates>::iterator const& collision, PhiDaughterCandidates const&, V0Candidates const& V0s, V0DaughterCandidates const&)
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
      const auto& posDaughterTrack = v0.posTrack_as<V0DaughterCandidates>();
      const auto& negDaughterTrack = v0.negTrack_as<V0DaughterCandidates>();

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

      float weightInclusive = 1. / static_cast<float>(countInclusive);
      float weightLtFirstCut = 1. / static_cast<float>(countLtFirstCut);
      float weightLtSecondCut = 1. / static_cast<float>(countLtSecondCut);

      switch (iBin) {
        case 0: {
          fillInvMass2D<false, 0>(recK0S, listrecPhi, multiplicity, weightInclusive, weightLtFirstCut, weightLtSecondCut);
          break;
        }
        case 1: {
          fillInvMass2D<false, 1>(recK0S, listrecPhi, multiplicity, weightInclusive, weightLtFirstCut, weightLtSecondCut);
          break;
        }
        case 2: {
          fillInvMass2D<false, 2>(recK0S, listrecPhi, multiplicity, weightInclusive, weightLtFirstCut, weightLtSecondCut);
          break;
        }
        case 3: {
          fillInvMass2D<false, 3>(recK0S, listrecPhi, multiplicity, weightInclusive, weightLtFirstCut, weightLtSecondCut);
          break;
        }
        case 4: {
          fillInvMass2D<false, 4>(recK0S, listrecPhi, multiplicity, weightInclusive, weightLtFirstCut, weightLtSecondCut);
          break;
        }
        case 5: {
          fillInvMass2D<false, 5>(recK0S, listrecPhi, multiplicity, weightInclusive, weightLtFirstCut, weightLtSecondCut);
          break;
        }
        case 6: {
          fillInvMass2D<false, 6>(recK0S, listrecPhi, multiplicity, weightInclusive, weightLtFirstCut, weightLtSecondCut);
          break;
        }
        case 7: {
          fillInvMass2D<false, 7>(recK0S, listrecPhi, multiplicity, weightInclusive, weightLtFirstCut, weightLtSecondCut);
          break;
        }
        case 8: {
          fillInvMass2D<false, 8>(recK0S, listrecPhi, multiplicity, weightInclusive, weightLtFirstCut, weightLtSecondCut);
          break;
        }
        case 9: {
          fillInvMass2D<false, 9>(recK0S, listrecPhi, multiplicity, weightInclusive, weightLtFirstCut, weightLtSecondCut);
          break;
        }
        default:
          break;
      }
    }
  }

  PROCESS_SWITCH(phik0shortanalysis, processSEPhiK0S, "Process Same Event for Analysis Phi-K0S", false);

  void processSEPhiPion(soa::Filtered<EventCandidates>::iterator const& collision, PhiDaughterCandidates const&, PionCandidates const& Pis)
  {
    if (!collision.isInelGt0())
      return;

    float multiplicity = collision.centFT0M();
    eventHist.fill(HIST("hMultiplicityPercent"), multiplicity);

    // Defining positive and negative tracks for phi reconstruction
    auto posThisColl = posTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto negThisColl = negTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);

    // Loop over all primary pion candidates
    for (const auto& pi : Pis) {

      // Pion selection
      if (!selectionPion(pi))
        continue;

      TLorentzVector recPi;
      recPi.SetXYZM(pi.px(), pi.py(), pi.pz(), massPi);
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

          countInclusive++;
          if (std::abs(recPi.Rapidity() - recPhi.Rapidity()) > cfgFirstCutonDeltay)
            continue;
          countLtFirstCut++;
          if (std::abs(recPi.Rapidity() - recPhi.Rapidity()) > cfgSecondCutonDeltay)
            continue;
          countLtSecondCut++;
        }
      }

      float weightInclusive = 1. / static_cast<float>(countInclusive);
      float weightLtFirstCut = 1. / static_cast<float>(countLtFirstCut);
      float weightLtSecondCut = 1. / static_cast<float>(countLtSecondCut);
      fillInvMassdEdx<false>(recPi, pi.tpcNSigmaPi(), listrecPhi, multiplicity, weightInclusive, weightLtFirstCut, weightLtSecondCut);
    }
  }

  PROCESS_SWITCH(phik0shortanalysis, processSEPhiPion, "Process Same Event for Analysis Phi-Pion", false);

  void processMEPhiK0S(soa::Filtered<EventCandidates> const& collisions, PhiDaughterCandidates const&, V0Candidates const& V0s, V0DaughterCandidates const&)
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
        const auto& posDaughterTrack = v0.posTrack_as<V0DaughterCandidates>();
        const auto& negDaughterTrack = v0.negTrack_as<V0DaughterCandidates>();

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
            fillInvMass2D<true, 0>(recK0S, listrecPhi, multiplicity, weightInclusive, weightLtFirstCut, weightLtSecondCut);
            break;
          }
          case 1: {
            fillInvMass2D<true, 1>(recK0S, listrecPhi, multiplicity, weightInclusive, weightLtFirstCut, weightLtSecondCut);
            break;
          }
          case 2: {
            fillInvMass2D<true, 2>(recK0S, listrecPhi, multiplicity, weightInclusive, weightLtFirstCut, weightLtSecondCut);
            break;
          }
          case 3: {
            fillInvMass2D<true, 3>(recK0S, listrecPhi, multiplicity, weightInclusive, weightLtFirstCut, weightLtSecondCut);
            break;
          }
          case 4: {
            fillInvMass2D<true, 4>(recK0S, listrecPhi, multiplicity, weightInclusive, weightLtFirstCut, weightLtSecondCut);
            break;
          }
          case 5: {
            fillInvMass2D<true, 5>(recK0S, listrecPhi, multiplicity, weightInclusive, weightLtFirstCut, weightLtSecondCut);
            break;
          }
          case 6: {
            fillInvMass2D<true, 6>(recK0S, listrecPhi, multiplicity, weightInclusive, weightLtFirstCut, weightLtSecondCut);
            break;
          }
          case 7: {
            fillInvMass2D<true, 7>(recK0S, listrecPhi, multiplicity, weightInclusive, weightLtFirstCut, weightLtSecondCut);
            break;
          }
          case 8: {
            fillInvMass2D<true, 8>(recK0S, listrecPhi, multiplicity, weightInclusive, weightLtFirstCut, weightLtSecondCut);
            break;
          }
          case 9: {
            fillInvMass2D<true, 9>(recK0S, listrecPhi, multiplicity, weightInclusive, weightLtFirstCut, weightLtSecondCut);
            break;
          }
          default:
            break;
        }
      }
    }
  }

  PROCESS_SWITCH(phik0shortanalysis, processMEPhiK0S, "Process Mixed Event for Analysis Phi-K0S", false);

  void processMEPhiPion(soa::Filtered<EventCandidates> const& collisions, PhiDaughterCandidates const&, PionCandidates const& Pis)
  {
    // Mixing the events with similar vertex z and multiplicity
    BinningTypeVertexContributor binningOnPositions{{axisVertex, axisMultiplicity}, true};
    for (auto const& [collision1, collision2] : o2::soa::selfCombinations(binningOnPositions, cfgNoMixedEvents, -1, collisions, collisions)) {
      if (!collision1.isInelGt0() || !collision2.isInelGt0())
        continue;

      float multiplicity = collision1.centFT0M();

      // Defining V0s from collision1
      auto PiThisColl = Pis.sliceByCached(aod::track::collisionId, collision1.globalIndex(), cache);

      // Defining positive and negative tracks for phi reconstruction from collision1 and collision2, respectively
      auto posThisColl = posTracks->sliceByCached(aod::track::collisionId, collision1.globalIndex(), cache);
      auto negThisColl = negTracks->sliceByCached(aod::track::collisionId, collision2.globalIndex(), cache);

      for (const auto& pi : PiThisColl) {

        if (!selectionPion(pi))
          continue;

        TLorentzVector recPi;
        recPi.SetXYZM(pi.px(), pi.py(), pi.pz(), massPi);
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

        float weightInclusive = 1. / static_cast<float>(countInclusive);
        float weightLtFirstCut = 1. / static_cast<float>(countLtFirstCut);
        float weightLtSecondCut = 1. / static_cast<float>(countLtSecondCut);
        fillInvMassdEdx<true>(recPi, pi.tpcNSigmaPi(), listrecPhi, multiplicity, weightInclusive, weightLtFirstCut, weightLtSecondCut);
      }
    }
  }

  PROCESS_SWITCH(phik0shortanalysis, processMEPhiPion, "Process Mixed Event for Analysis Phi-Pion", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<phik0shortanalysis>(cfgc, TaskName{"lf-phik0shortanalysis"})};
}
