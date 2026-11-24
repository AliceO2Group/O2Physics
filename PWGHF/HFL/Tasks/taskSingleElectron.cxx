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
/// \file taskSingleElectron.cxx
/// \brief task for electrons from heavy-flavour hadron decays
/// \author Jonghan Park (Jeonbuk National University), Seul I Jeong (Pusan National University)

#include "PWGHF/DataModel/AliasTables.h"
#include "PWGHF/Utils/utilsBfieldCCDB.h"

#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/PhysicsConstants.h>
#include <CommonUtils/ConfigurableParam.h>
#include <DetectorsBase/MatLayerCylSet.h>
#include <DetectorsBase/Propagator.h>
#include <DetectorsVertexing/PVertexer.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/DataTypes.h>
#include <Framework/Expressions.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/SliceCache.h>
#include <Framework/runDataProcessing.h>
#include <ReconstructionDataFormats/Track.h>
#include <ReconstructionDataFormats/Vertex.h>

#include <TPDGCode.h>

#include <algorithm>
#include <array>
#include <cstdint>
#include <map>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::constants::math;
using namespace o2::constants::physics;
using namespace o2::framework;
using namespace o2::framework::expressions;

enum SourceType {
  NotElec = 0,      // not electron
  DirectCharm = 1,  // electrons from prompt charm hadrons
  DirectBeauty = 2, // electrons from primary beauty hadrons
  BeautyCharm = 3,  // electrons from non-prompt charm hadrons
  DirectGamma = 4,  // electrons from direct photon
  GammaPi0 = 5,
  GammaEta = 6,
  GammaOmega = 7,
  GammaPhi = 8,
  GammaEtaPrime = 9,
  GammaRho0 = 10,
  GammaK0s = 11,
  GammaK0l = 12,
  GammaKe3 = 13,
  GammaLambda0 = 14,
  GammaSigma = 15,
  Pi0 = 16,
  Eta = 17,
  Omega = 18,
  Phi = 19,
  EtaPrime = 20,
  Rho0 = 21,
  K0s = 22,
  K0l = 23,
  Ke3 = 24,
  Lambda0 = 25,
  Sigma = 26,
  Else = 27
};

struct HfTaskSingleElectron {

  // Produces

  // Configurable
  Configurable<int> nContribMin{"nContribMin", 2, "min number of contributors"};
  Configurable<float> posZMax{"posZMax", 10., "max posZ cut"};
  Configurable<float> ptTrackMax{"ptTrackMax", 10., "max pt cut"};
  Configurable<float> ptTrackMin{"ptTrackMin", 0.5, "min pt cut"};
  Configurable<float> etaTrackMax{"etaTrackMax", 0.8, "eta cut"};
  Configurable<int> nCrossedRowTpcMin{"nCrossedRowTpcMin", 70, "min # of TPC n cluster crossed rows"};
  Configurable<float> nCrossedRowsOverFindableClsTpcMin{"nCrossedRowsOverFindableClsTpcMin", 0.8, "min ratio of TPC crossed rows over findable clusters"};
  Configurable<float> chi2PerNClTpcMax{"chi2PerNClTpcMax", 4., "max # of tpc chi2 per clusters"};
  Configurable<int> clsIbItsMin{"clsIbItsMin", 3, "min # of its clusters in IB"};
  Configurable<float> chi2PerNClItsMax{"chi2PerNClItsMax", 6., "min # of its chi2 per clusters"};
  Configurable<float> dcaxyMax{"dcaxyMax", 1., "max of track dca in xy"};
  Configurable<float> dcazMax{"dcazMax", 2., "max of track dca in z"};
  Configurable<int> nClsTpcMin{"nClsTpcMin", 0, "min # of found TPC clusters"};
  Configurable<int> nClsItsMin{"nClsItsMin", 0, "min # of total ITS clusters"};
  Configurable<bool> requireItsTpcRefit{"requireItsTpcRefit", true, "require ITS and TPC refit"};
  Configurable<float> nSigmaTofMax{"nSigmaTofMax", 3., "max of tof nsigma"};
  Configurable<float> nSigmaTpcMin{"nSigmaTpcMin", -1., "min of tpc nsigma"};
  Configurable<float> nSigmaTpcMax{"nSigmaTpcMax", 3., "max of tpc nsigma"};

  Configurable<int> nBinsP{"nBinsP", 1500, "number of bins of particle momentum"};
  Configurable<int> nBinsPt{"nBinsPt", 100, "N bins in pT histo"};

  Configurable<int> nSigmaTpcHadronMax{"nSigmaTpcHadronMax", -3, "max of tpc hadron nsigma"};
  Configurable<int> nSigmaTpcHadronMin{"nSigmaTpcHadronMin", -5, "min of tpc hadron nsigma"};

  Configurable<float> bFieldSign{"bFieldSign", 1.f, "sign of the L3 magnetic field (+1 or -1), used in DCA_xy x charge x sign(B)"};

  // PV refit (processDataPvRefit)
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "URL of the CCDB repository"};
  Configurable<std::string> ccdbPathLut{"ccdbPathLut", "GLO/Param/MatLUT", "CCDB path of the material LUT"};
  Configurable<std::string> ccdbPathGrpMag{"ccdbPathGrpMag", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};

  // SliceCache
  SliceCache cache;

  Service<o2::ccdb::BasicCCDBManager> ccdb{};
  o2::base::MatLayerCylSet* lut{};
  o2::vertexing::PVertexer vertexer;
  int runNumber{0};

  // using declarations
  using MyCollisions = soa::Join<aod::Collisions, aod::EvSels>;
  using TracksEl = soa::Join<aod::Tracks, aod::TrackSelection, aod::TrackSelectionExtension, aod::TracksExtra, aod::TracksDCA, aod::pidTOFFullEl, aod::pidTPCFullEl>;
  using McTracksEl = soa::Join<aod::Tracks, aod::TrackSelection, aod::TrackSelectionExtension, aod::TrackExtra, aod::TracksDCA, aod::pidTOFFullEl, aod::pidTPCFullEl, aod::McTrackLabels>;
  using AllTracks = aod::TracksWCovDcaExtra;

  // Filter
  Filter collZFilter = nabs(aod::collision::posZ) < posZMax;

  // Partition
  Partition<AllTracks> pvContributors = ((aod::track::flags & static_cast<uint32_t>(aod::track::PVContributor)) == static_cast<uint32_t>(aod::track::PVContributor));

  // ConfigurableAxis
  ConfigurableAxis axisPtEl{"axisPtEl", {VARIABLE_WIDTH, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.75f, 2.0f, 2.25f, 2.5f, 2.75f, 3.f, 3.5f, 4.0f, 5.0f, 6.0f, 8.0f, 10.0f}, "electron pt bins"};
  ConfigurableAxis axisDcaSigned{"axisDcaSigned", {4000, -0.2, 0.2}, "DCA_xy x charge x sign(B) bins (cm)"};
  ConfigurableAxis axisNContrib{"axisNContrib", {VARIABLE_WIDTH, 0., 5., 10., 15., 20., 25., 30., 40., 50., 60., 80., 100., 150., 250.}, "number of PV contributors bins"};

  // Histogram registry
  HistogramRegistry histos{"histos"};

  void init(InitContext const&)
  {
    // AxisSpec
    const AxisSpec axisEvt{4, 0., 4., "nEvents"};
    const AxisSpec axisNCont{100, 0., 100., "nCont"};
    const AxisSpec axisPosZ{600, -30., 30., "Z_{pos}"};
    const AxisSpec axisEta{30, -1.5, +1.5, "#eta"};
    const AxisSpec axisP{nBinsP, 0., 15., "p_{T}"};
    const AxisSpec axisPt{nBinsPt, 0., 15., "p_{T}"};
    const AxisSpec axisNsig{800, -20., 20.};
    const AxisSpec axisTrackIp{4000, -0.2, 0.2, "dca"};

    // create histograms
    histos.add("hNEvents", "Number of events", kTH1D, {{1, 0., 1.}});
    histos.add("hVtxZ", "hVtxZ; cm; entries", kTH1D, {axisPosZ});
    histos.add("hEtaTrack", "hEtaTrack; #eta; entries", kTH1D, {axisEta});
    histos.add("hPtTrack", "#it{p}_{T} distribution of selected tracks; #it{p}_{T} (GeV/#it{c}); entries", kTH1D, {axisPt});

    // QA plots for trigger track selection
    histos.add("hNClsTpcTrack", "hNClsTpcTrack", kTH1D, {{200, 0, 200}});
    histos.add("hCrossedRowsOverFindableTpcTrack", "", kTH1D, {{200, 0, 2}});
    histos.add("hChi2TpcTrack", "", kTH1D, {{100, 0, 10}});
    histos.add("hIbClsItsTrack", "", kTH1D, {{10, 0, 10}});
    histos.add("hChi2ItsTrack", "", kTH1D, {{50, 0, 50}});
    histos.add("hDcaXYTrack", "", kTH1D, {{600, -3, 3}});
    histos.add("hDcaZTrack", "", kTH1D, {{600, -3, 3}});

    // QA of track-selection variables before any track cut, differential in pT (binning follows DPG qaEventTrack)
    histos.add("hPtEtaPreCut", "before track selection;#it{p}_{T} (GeV/#it{c});#eta", kTH2D, {{axisPtEl}, {40, -1., 1.}});
    histos.add("hPtCrossedRowsTpcPreCut", "before track selection;#it{p}_{T} (GeV/#it{c});TPC crossed rows", kTH2D, {{axisPtEl}, {165, -0.5, 164.5}});
    histos.add("hPtCrossedRowsOverFindableTpcPreCut", "before track selection;#it{p}_{T} (GeV/#it{c});TPC crossed rows / findable clusters", kTH2D, {{axisPtEl}, {200, 0., 2.}});
    histos.add("hPtNClsFoundTpcPreCut", "before track selection;#it{p}_{T} (GeV/#it{c});TPC found clusters", kTH2D, {{axisPtEl}, {165, -0.5, 164.5}});
    histos.add("hPtNClsFindableTpcPreCut", "before track selection;#it{p}_{T} (GeV/#it{c});TPC findable clusters", kTH2D, {{axisPtEl}, {165, -0.5, 164.5}});
    histos.add("hPtFoundOverFindableTpcPreCut", "before track selection;#it{p}_{T} (GeV/#it{c});TPC found / findable clusters", kTH2D, {{axisPtEl}, {200, 0., 2.}});
    histos.add("hPtFracSharedClsTpcPreCut", "before track selection;#it{p}_{T} (GeV/#it{c});fraction of shared TPC clusters", kTH2D, {{axisPtEl}, {100, 0., 1.}});
    histos.add("hPtChi2TpcPreCut", "before track selection;#it{p}_{T} (GeV/#it{c});#chi^{2} / cluster TPC", kTH2D, {{axisPtEl}, {100, 0., 10.}});
    histos.add("hPtNClsItsPreCut", "before track selection;#it{p}_{T} (GeV/#it{c});# clusters ITS", kTH2D, {{axisPtEl}, {8, -0.5, 7.5}});
    histos.add("hPtIbClsItsPreCut", "before track selection;#it{p}_{T} (GeV/#it{c});# clusters ITS inner barrel", kTH2D, {{axisPtEl}, {4, -0.5, 3.5}});
    histos.add("hPtChi2ItsPreCut", "before track selection;#it{p}_{T} (GeV/#it{c});#chi^{2} / cluster ITS", kTH2D, {{axisPtEl}, {100, 0., 40.}});
    histos.add("hPtDcaXYPreCut", "before track selection;#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm)", kTH2D, {{axisPtEl}, {600, -3., 3.}});
    histos.add("hPtDcaZPreCut", "before track selection;#it{p}_{T} (GeV/#it{c});DCA_{z} (cm)", kTH2D, {{axisPtEl}, {600, -3., 3.}});

    // pid
    histos.add("hTofNSigPt", "", kTH2D, {{axisPtEl}, {axisNsig}});
    histos.add("hTofNSigPtQA", "", kTH2D, {{axisPtEl}, {axisNsig}});
    histos.add("hTpcNSigP", "", kTH2D, {{axisP}, {axisNsig}});
    histos.add("hTpcNSigPt", "", kTH2D, {{axisPtEl}, {axisNsig}});
    histos.add("hTpcNSigPAfterTofCut", "", kTH2D, {{axisP}, {axisNsig}});
    histos.add("hTpcNSigPtAfterTofCut", "", kTH2D, {{axisPtEl}, {axisNsig}});
    histos.add("hTpcNSigPtQA", "", kTH2D, {{axisPtEl}, {axisNsig}});

    // TPC nSigma(e) vs pT per MC truth species (electron, pion, kaon, proton, other), after the track selection and after the TOF electron selection
    histos.add("hTpcNSigPtRecElectron", "MC truth e^{#pm} after track selection;#it{p}_{T} (GeV/#it{c});n#sigma_{e}^{TPC}", kTH2D, {{axisPtEl}, {axisNsig}});
    histos.add("hTpcNSigPtRecPion", "MC truth #pi^{#pm} after track selection;#it{p}_{T} (GeV/#it{c});n#sigma_{e}^{TPC}", kTH2D, {{axisPtEl}, {axisNsig}});
    histos.add("hTpcNSigPtRecKaon", "MC truth K^{#pm} after track selection;#it{p}_{T} (GeV/#it{c});n#sigma_{e}^{TPC}", kTH2D, {{axisPtEl}, {axisNsig}});
    histos.add("hTpcNSigPtRecProton", "MC truth p after track selection;#it{p}_{T} (GeV/#it{c});n#sigma_{e}^{TPC}", kTH2D, {{axisPtEl}, {axisNsig}});
    histos.add("hTpcNSigPtRecOther", "MC truth other after track selection;#it{p}_{T} (GeV/#it{c});n#sigma_{e}^{TPC}", kTH2D, {{axisPtEl}, {axisNsig}});
    histos.add("hTpcNSigPtAfterTofCutRecElectron", "MC truth e^{#pm} after TOF selection;#it{p}_{T} (GeV/#it{c});n#sigma_{e}^{TPC}", kTH2D, {{axisPtEl}, {axisNsig}});
    histos.add("hTpcNSigPtAfterTofCutRecPion", "MC truth #pi^{#pm} after TOF selection;#it{p}_{T} (GeV/#it{c});n#sigma_{e}^{TPC}", kTH2D, {{axisPtEl}, {axisNsig}});
    histos.add("hTpcNSigPtAfterTofCutRecKaon", "MC truth K^{#pm} after TOF selection;#it{p}_{T} (GeV/#it{c});n#sigma_{e}^{TPC}", kTH2D, {{axisPtEl}, {axisNsig}});
    histos.add("hTpcNSigPtAfterTofCutRecProton", "MC truth p after TOF selection;#it{p}_{T} (GeV/#it{c});n#sigma_{e}^{TPC}", kTH2D, {{axisPtEl}, {axisNsig}});
    histos.add("hTpcNSigPtAfterTofCutRecOther", "MC truth other after TOF selection;#it{p}_{T} (GeV/#it{c});n#sigma_{e}^{TPC}", kTH2D, {{axisPtEl}, {axisNsig}});

    // track impact parameter
    histos.add("hDcaTrack", "", kTH2D, {{axisPtEl}, {axisTrackIp}});
    histos.add("hDcaBeauty", "", kTH2D, {{axisPtEl}, {axisTrackIp}});
    histos.add("hDcaCharm", "", kTH2D, {{axisPtEl}, {axisTrackIp}});
    histos.add("hDcaDalitz", "", kTH2D, {{axisPtEl}, {axisTrackIp}});
    histos.add("hDcaConv", "", kTH2D, {{axisPtEl}, {axisTrackIp}});
    histos.add("hDcaHadron", "", kTH2D, {{axisPtEl}, {axisTrackIp}});

    // DCA_xy x q x sign(B)
    const AxisSpec axisDcaSig{axisDcaSigned, "DCA_{xy} #times #it{q} #times sign(#it{B}) (cm)"};
    const AxisSpec axisNContribEvt{axisNContrib, "#it{N}_{contrib}"};
    histos.add("hDcaSignedTrack", "after electron identification;#it{p}_{T} (GeV/#it{c});DCA_{xy} #times #it{q} #times sign(#it{B}) (cm);#it{N}_{contrib}", kTH3D, {{axisPtEl}, {axisDcaSig}, {axisNContribEvt}});
    histos.add("hDcaSignedHadron", "TPC hadron band, before TOF selection;#it{p}_{T} (GeV/#it{c});DCA_{xy} #times #it{q} #times sign(#it{B}) (cm)", kTH2D, {{axisPtEl}, {axisDcaSig}});
    histos.add("hDcaSignedHadronAfterTofCut", "TPC hadron band, after TOF selection;#it{p}_{T} (GeV/#it{c});DCA_{xy} #times #it{q} #times sign(#it{B}) (cm)", kTH2D, {{axisPtEl}, {axisDcaSig}});

    if (doprocessData && doprocessDataPvRefit) {
      LOGF(fatal, "processData and processDataPvRefit cannot be enabled together");
    }
    if (doprocessDataPvRefit) {
      histos.add("hDcaSignedTrackPvRefit", "after electron identification, PV refit;#it{p}_{T} (GeV/#it{c});DCA_{xy} #times #it{q} #times sign(#it{B}) (cm);#it{N}_{contrib}", kTH3D, {{axisPtEl}, {axisDcaSig}, {axisNContribEvt}});
      histos.add("hDcaSignedHadronPvRefit", "TPC hadron band, before TOF selection, PV refit;#it{p}_{T} (GeV/#it{c});DCA_{xy} #times #it{q} #times sign(#it{B}) (cm)", kTH2D, {{axisPtEl}, {axisDcaSig}});
      histos.add("hDcaSignedHadronAfterTofCutPvRefit", "TPC hadron band, after TOF selection, PV refit;#it{p}_{T} (GeV/#it{c});DCA_{xy} #times #it{q} #times sign(#it{B}) (cm)", kTH2D, {{axisPtEl}, {axisDcaSig}});
      // PV refit QA
      const AxisSpec axisCollisionNContrib{1000, 0, 1000, "Number of contributors"};
      const AxisSpec axisPvRefitDeltaX{1000, -0.05f, 0.05f, "#Delta x_{PV} (cm)"};
      const AxisSpec axisPvRefitDeltaY{1000, -0.05f, 0.05f, "#Delta y_{PV} (cm)"};
      const AxisSpec axisPvRefitDeltaZ{1000, -0.05f, 0.05f, "#Delta z_{PV} (cm)"};
      histos.add("PvRefit/hVerticesPerTrack", "", kTH1D, {{3, 0.5f, 3.5f, ""}});
      histos.add("PvRefit/hPvDeltaXvsNContrib", "", kTH2D, {axisCollisionNContrib, axisPvRefitDeltaX});
      histos.add("PvRefit/hPvDeltaYvsNContrib", "", kTH2D, {axisCollisionNContrib, axisPvRefitDeltaY});
      histos.add("PvRefit/hPvDeltaZvsNContrib", "", kTH2D, {axisCollisionNContrib, axisPvRefitDeltaZ});
      histos.add("PvRefit/hChi2vsNContrib", "", kTH2D, {axisCollisionNContrib, {102, -1.5, 100.5, "#chi^{2} PV refit"}});
      histos.add("PvRefit/hPvRefitXChi2Minus1", "PV refit with #it{#chi}^{2}==#minus1", kTH2D, {{100, -20.f, 20.f}, {100, -2.f, 2.f}});
      histos.add("PvRefit/hPvRefitYChi2Minus1", "PV refit with #it{#chi}^{2}==#minus1", kTH2D, {{100, -20.f, 20.f}, {100, -2.f, 2.f}});
      histos.add("PvRefit/hPvRefitZChi2Minus1", "PV refit with #it{#chi}^{2}==#minus1", kTH2D, {{100, -20.f, 20.f}, {100, -2.f, 2.f}});
      histos.add("PvRefit/hNContribPvRefitNotDoable", "N. contrib for PV refit not doable", kTH1D, {axisCollisionNContrib});
      histos.add("PvRefit/hNContribPvRefitChi2Minus1", "N. contrib for PV refit #it{#chi}^{2}==#minus1", kTH1D, {axisCollisionNContrib});
      ccdb->setURL(ccdbUrl);
      ccdb->setCaching(true);
      ccdb->setLocalObjectValidityChecking();
      lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>(ccdbPathLut));
      o2::conf::ConfigurableParam::updateFromString("pvertexer.useMeanVertexConstraint=false");
    }

    // QA plots for MC
    histos.add("hPdgC", "", kTH1D, {{10001, -0.5, 10000.5}});
    histos.add("hPdgB", "", kTH1D, {{10001, -0.5, 10000.5}});
    histos.add("hPdgDa", "", kTH1D, {{10001, -0.5, 10000.5}});
    histos.add("hPdgCo", "", kTH1D, {{10001, -0.5, 10000.5}});
  }

  template <typename TrackType>
  void fillTrackQaPreCut(const TrackType& track)
  {
    double const pt = track.pt();
    histos.fill(HIST("hPtEtaPreCut"), pt, track.eta());
    histos.fill(HIST("hPtCrossedRowsTpcPreCut"), pt, track.tpcNClsCrossedRows());
    histos.fill(HIST("hPtCrossedRowsOverFindableTpcPreCut"), pt, track.tpcCrossedRowsOverFindableCls());
    histos.fill(HIST("hPtNClsFoundTpcPreCut"), pt, track.tpcNClsFound());
    histos.fill(HIST("hPtNClsFindableTpcPreCut"), pt, track.tpcNClsFindable());
    histos.fill(HIST("hPtFoundOverFindableTpcPreCut"), pt, track.tpcFoundOverFindableCls());
    histos.fill(HIST("hPtFracSharedClsTpcPreCut"), pt, track.tpcFractionSharedCls());
    histos.fill(HIST("hPtChi2TpcPreCut"), pt, track.tpcChi2NCl());
    histos.fill(HIST("hPtNClsItsPreCut"), pt, track.itsNCls());
    histos.fill(HIST("hPtIbClsItsPreCut"), pt, track.itsNClsInnerBarrel());
    histos.fill(HIST("hPtChi2ItsPreCut"), pt, track.itsChi2NCl());
    histos.fill(HIST("hPtDcaXYPreCut"), pt, track.dcaXY());
    histos.fill(HIST("hPtDcaZPreCut"), pt, track.dcaZ());
  }

  template <typename TrackType>
  void fillTpcNSigPtRec(const TrackType& track, bool afterTofCut)
  {
    if (!track.has_mcParticle()) {
      return;
    }
    double const pt = track.pt();
    double const nSigmaEl = track.tpcNSigmaEl();
    int const absPdg = std::abs(track.mcParticle().pdgCode());
    if (!afterTofCut) {
      if (absPdg == kElectron) {
        histos.fill(HIST("hTpcNSigPtRecElectron"), pt, nSigmaEl);
      } else if (absPdg == kPiPlus) {
        histos.fill(HIST("hTpcNSigPtRecPion"), pt, nSigmaEl);
      } else if (absPdg == kKPlus) {
        histos.fill(HIST("hTpcNSigPtRecKaon"), pt, nSigmaEl);
      } else if (absPdg == kProton) {
        histos.fill(HIST("hTpcNSigPtRecProton"), pt, nSigmaEl);
      } else {
        histos.fill(HIST("hTpcNSigPtRecOther"), pt, nSigmaEl);
      }
      return;
    }
    if (absPdg == kElectron) {
      histos.fill(HIST("hTpcNSigPtAfterTofCutRecElectron"), pt, nSigmaEl);
    } else if (absPdg == kPiPlus) {
      histos.fill(HIST("hTpcNSigPtAfterTofCutRecPion"), pt, nSigmaEl);
    } else if (absPdg == kKPlus) {
      histos.fill(HIST("hTpcNSigPtAfterTofCutRecKaon"), pt, nSigmaEl);
    } else if (absPdg == kProton) {
      histos.fill(HIST("hTpcNSigPtAfterTofCutRecProton"), pt, nSigmaEl);
    } else {
      histos.fill(HIST("hTpcNSigPtAfterTofCutRecOther"), pt, nSigmaEl);
    }
  }

  template <typename TrackType>
  bool trackSel(const TrackType& track)
  {
    if ((track.pt() > ptTrackMax) || (track.pt() < ptTrackMin)) {
      return false;
    }
    if (std::abs(track.eta()) > etaTrackMax) {
      return false;
    }

    if (track.tpcNClsCrossedRows() < nCrossedRowTpcMin) {
      return false;
    }

    if (track.tpcCrossedRowsOverFindableCls() < nCrossedRowsOverFindableClsTpcMin) {
      return false;
    }

    if (track.tpcNClsFound() < nClsTpcMin) {
      return false;
    }

    if (track.tpcChi2NCl() > chi2PerNClTpcMax) {
      return false;
    }

    if (track.itsNClsInnerBarrel() < clsIbItsMin) { // inner barrel has 3 layers: default 3 is equivalent to the previous == requirement
      return false;
    }

    if (track.itsNCls() < nClsItsMin) {
      return false;
    }

    if (track.itsChi2NCl() > chi2PerNClItsMax) {
      return false;
    }

    if (requireItsTpcRefit && !(track.passedITSRefit() && track.passedTPCRefit())) {
      return false;
    }

    if (std::abs(track.dcaXY()) > dcaxyMax) {
      return false;
    }

    if (std::abs(track.dcaZ()) > dcazMax) {
      return false;
    }

    return true;
  }

  template <typename TrackType>
  int getElecSource(const TrackType& track, double& mpt, int& mpdg)
  {
    auto mcpart = track.mcParticle();
    if (std::abs(mcpart.pdgCode()) != kElectron) {
      return NotElec;
    }

    int motherPdg = -999;
    int grmotherPdg = -999;
    int ggrmotherPdg = -999; // mother, grand mother, grand grand mother pdg
    int motherPt = -999.;
    int grmotherPt = -999;
    int ggrmotherPt = -999.; // mother, grand mother, grand grand mother pt

    auto partMother = mcpart.template mothers_as<aod::McParticles>(); // first mother particle of electron
    auto partMotherCopy = partMother;                                 // copy of the first mother
    auto mctrack = partMother;                                        // will change all the time

    motherPt = partMother.front().pt();                 // first mother pt
    motherPdg = std::abs(partMother.front().pdgCode()); // first mother pdg
    mpt = motherPt;                                     // copy of first mother pt
    mpdg = motherPdg;                                   // copy of first mother pdg

    // check if electron from charm hadrons
    if ((static_cast<int>(motherPdg / 100.) % 10) == kCharm || (static_cast<int>(motherPdg / 1000.) % 10) == kCharm) {

      // iterate until B hadron is found as an ancestor
      while (partMother.size()) {
        mctrack = partMother.front().template mothers_as<aod::McParticles>();
        if (mctrack.size()) {
          auto const& grmothersIdsVec = mctrack.front().mothersIds();

          if (grmothersIdsVec.empty()) {
            return DirectCharm;
          }
          grmotherPt = mctrack.front().pt();
          grmotherPdg = std::abs(mctrack.front().pdgCode());
          if ((static_cast<int>(grmotherPdg / 100.) % 10) == kBottom || (static_cast<int>(grmotherPdg / 1000.) % 10) == kBottom) {
            mpt = grmotherPt;
            mpdg = grmotherPdg;
            return BeautyCharm;
          }
        }
        partMother = mctrack;
      }
    } else if ((static_cast<int>(motherPdg / 100.) % 10) == kBottom || (static_cast<int>(motherPdg / 1000.) % 10) == kBottom) { // check if electron from beauty hadrons
      return DirectBeauty;
    } else if (motherPdg == kGamma) { // check if electron from photon conversion
      mctrack = partMother.front().template mothers_as<aod::McParticles>();
      if (mctrack.size()) {
        auto const& grmothersIdsVec = mctrack.front().mothersIds();
        if (grmothersIdsVec.empty()) {
          return DirectGamma;
        }
        grmotherPdg = std::abs(mctrack.front().pdgCode());
        mpdg = grmotherPdg;
        mpt = mctrack.front().pt();

        partMother = mctrack;
        mctrack = partMother.front().template mothers_as<aod::McParticles>();
        if (mctrack.size()) {
          auto const& ggrmothersIdsVec = mctrack.front().mothersIds();
          if (ggrmothersIdsVec.empty()) {
            if (grmotherPdg == kPi0) {
              return GammaPi0;
            }
            if (grmotherPdg == Pdg::kEta) {
              return GammaEta;
            }
            if (grmotherPdg == Pdg::kOmega) {
              return GammaOmega;
            }
            if (grmotherPdg == Pdg::kPhi) {
              return GammaPhi;
            }
            if (grmotherPdg == Pdg::kEtaPrime) {
              return GammaEtaPrime;
            }
            if (grmotherPdg == kRho770_0) {
              return GammaRho0;
            }
            return Else;
          }
          ggrmotherPdg = mctrack.front().pdgCode();
          ggrmotherPt = mctrack.front().pt();
          mpdg = ggrmotherPdg;
          mpt = ggrmotherPt;
          if (grmotherPdg == kPi0) {
            if (ggrmotherPdg == kK0Short) {
              return GammaK0s;
            }
            if (ggrmotherPdg == kK0Long) {
              return GammaK0l;
            }
            if (ggrmotherPdg == kKPlus) {
              return GammaKe3;
            }
            if (ggrmotherPdg == kLambda0) {
              return GammaLambda0;
            }
            if (ggrmotherPdg == kSigmaPlus) {
              return GammaSigma;
            }
            mpdg = grmotherPdg;
            mpt = grmotherPt;
            return GammaPi0;
          }
          if (grmotherPdg == Pdg::kEta) {
            mpdg = grmotherPdg;
            mpt = grmotherPt;
            return GammaEta;
          }
          if (grmotherPdg == Pdg::kOmega) {
            mpdg = grmotherPdg;
            mpt = grmotherPt;
            return GammaOmega;
          }
          if (grmotherPdg == Pdg::kPhi) {
            mpdg = grmotherPdg;
            mpt = grmotherPt;
            return GammaPhi;
          }
          if (grmotherPdg == Pdg::kEtaPrime) {
            mpdg = grmotherPdg;
            mpt = grmotherPt;
            return GammaEtaPrime;
          }
          if (grmotherPdg == kRho770_0) {
            mpdg = grmotherPdg;
            mpt = grmotherPt;
            return GammaRho0;
          }
          return Else;
        }
      }
    } else { // check if electron from Dalitz decays
      mctrack = partMother.front().template mothers_as<aod::McParticles>();
      if (mctrack.size()) {
        auto const& grmothersIdsVec = mctrack.front().mothersIds();
        if (grmothersIdsVec.empty()) {
          static const std::map<int, SourceType> pdgToSource = {
            {kPi0, Pi0},
            {Pdg::kEta, Eta},
            {Pdg::kOmega, Omega},
            {Pdg::kPhi, Phi},
            {Pdg::kEtaPrime, EtaPrime},
            {kRho770_0, Rho0},
            {kKPlus, Ke3},
            {kK0Long, K0l}};

          auto it = pdgToSource.find(motherPdg);
          if (it != pdgToSource.end()) {
            return it->second;
          }
          return Else;
        }
        if (motherPdg == kPi0) {
          grmotherPt = mctrack.front().pt();
          grmotherPdg = mctrack.front().pdgCode();
          mpt = grmotherPt;
          mpdg = grmotherPdg;
          if (grmotherPdg == kK0Short) {
            return K0s;
          }
          if (grmotherPdg == kK0Long) {
            return K0l;
          }
          if (grmotherPdg == kKPlus) {
            return Ke3;
          }
          if (grmotherPdg == kLambda0) {
            return Lambda0;
          }
          if (grmotherPdg == kSigmaPlus) {
            return Sigma;
          }
          mpt = motherPt;
          mpdg = motherPdg;
          return Pi0;
        }
        if (motherPdg == Pdg::kEta) {
          return Eta;
        }
        if (motherPdg == Pdg::kOmega) {
          return Omega;
        }
        if (motherPdg == Pdg::kPhi) {
          return Phi;
        }
        if (motherPdg == Pdg::kEtaPrime) {
          return EtaPrime;
        }
        if (motherPdg == kRho770_0) {
          return Rho0;
        }
        if (motherPdg == kKPlus) {
          return Ke3;
        }
        if (motherPdg == kK0Long) {
          return K0l;
        }
        return Else;
      }
    }

    return Else;
  }

  template <bool DoPvRefit, typename TCollision, typename TTracks>
  void runData(TCollision const& collision, TTracks const& tracks)
  {
    float const flagAnalysedEvt = 0.5;

    if (!collision.sel8()) {
      return;
    }

    if (collision.numContrib() < nContribMin) {
      return;
    }

    histos.fill(HIST("hVtxZ"), collision.posZ());
    histos.fill(HIST("hNEvents"), flagAnalysedEvt);

    // PV refit: vertexer per collision
    std::vector<int64_t> pvContribGlobId{};
    std::vector<o2::track::TrackParCov> pvContribTrackParCov{};
    o2::dataformats::VertexBase primVtx;
    bool pvRefitDoable = false;
    if constexpr (DoPvRefit) {
      initCCDB(collision.template bc_as<aod::BCsWithTimestamps>(), runNumber, ccdb, ccdbPathGrpMag, lut, false);
      auto pvContribColl = pvContributors->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
      for (const auto& contributor : pvContribColl) {
        pvContribGlobId.push_back(contributor.globalIndex());
        pvContribTrackParCov.push_back(getTrackParCov(contributor));
      }
      primVtx.setX(collision.posX());
      primVtx.setY(collision.posY());
      primVtx.setZ(collision.posZ());
      primVtx.setCov(collision.covXX(), collision.covXY(), collision.covYY(), collision.covXZ(), collision.covYZ(), collision.covZZ());
      vertexer.init();
      pvRefitDoable = vertexer.prepareVertexRefit(pvContribTrackParCov, primVtx);
    }
    std::vector<bool> pvContribUsed(pvContribGlobId.size(), true);

    for (const auto& track : tracks) {

      fillTrackQaPreCut(track);

      if (!trackSel(track)) {
        continue;
      }

      histos.fill(HIST("hEtaTrack"), track.eta());
      histos.fill(HIST("hPtTrack"), track.pt());

      histos.fill(HIST("hNClsTpcTrack"), track.tpcNClsCrossedRows());
      histos.fill(HIST("hCrossedRowsOverFindableTpcTrack"), track.tpcCrossedRowsOverFindableCls());
      histos.fill(HIST("hChi2TpcTrack"), track.tpcChi2NCl());
      histos.fill(HIST("hIbClsItsTrack"), track.itsNClsInnerBarrel());
      histos.fill(HIST("hChi2ItsTrack"), track.itsChi2NCl());
      histos.fill(HIST("hDcaXYTrack"), track.dcaXY());
      histos.fill(HIST("hDcaZTrack"), track.dcaZ());

      histos.fill(HIST("hTofNSigPt"), track.pt(), track.tofNSigmaEl());
      histos.fill(HIST("hTpcNSigP"), track.p(), track.tpcNSigmaEl());
      histos.fill(HIST("hTpcNSigPt"), track.pt(), track.tpcNSigmaEl());

      double const dcaSigned = track.dcaXY() * track.sign() * bFieldSign;
      float dcaXYPvRefit = track.dcaXY();
      bool pvRefitDone = false;
      [[maybe_unused]] auto getDcaSignedPvRefit = [&]() -> double {
        if (!pvRefitDone) {
          pvRefitDone = true;
          const auto it = std::find(pvContribGlobId.begin(), pvContribGlobId.end(), track.globalIndex());
          if (track.isPVContributor() && it != pvContribGlobId.end()) {
            histos.fill(HIST("PvRefit/hVerticesPerTrack"), 1);
            if (!pvRefitDoable) {
              histos.fill(HIST("PvRefit/hNContribPvRefitNotDoable"), collision.numContrib());
            } else {
              histos.fill(HIST("PvRefit/hVerticesPerTrack"), 2);
              const auto entry = std::distance(pvContribGlobId.begin(), it);
              pvContribUsed[entry] = false;
              const auto primVtxRefitted = vertexer.refitVertex(pvContribUsed, primVtx);
              pvContribUsed[entry] = true;
              histos.fill(HIST("PvRefit/hChi2vsNContrib"), primVtxRefitted.getNContributors(), primVtxRefitted.getChi2());
              if (primVtxRefitted.getChi2() < 0) {
                histos.fill(HIST("PvRefit/hPvRefitXChi2Minus1"), primVtxRefitted.getX(), collision.posX());
                histos.fill(HIST("PvRefit/hPvRefitYChi2Minus1"), primVtxRefitted.getY(), collision.posY());
                histos.fill(HIST("PvRefit/hPvRefitZChi2Minus1"), primVtxRefitted.getZ(), collision.posZ());
                histos.fill(HIST("PvRefit/hNContribPvRefitChi2Minus1"), collision.numContrib());
              } else {
                histos.fill(HIST("PvRefit/hVerticesPerTrack"), 3);
                histos.fill(HIST("PvRefit/hPvDeltaXvsNContrib"), primVtxRefitted.getNContributors(), primVtx.getX() - primVtxRefitted.getX());
                histos.fill(HIST("PvRefit/hPvDeltaYvsNContrib"), primVtxRefitted.getNContributors(), primVtx.getY() - primVtxRefitted.getY());
                histos.fill(HIST("PvRefit/hPvDeltaZvsNContrib"), primVtxRefitted.getNContributors(), primVtx.getZ() - primVtxRefitted.getZ());
                auto trackPar = getTrackPar(track);
                std::array<float, 2> dcaInfo{-999.f, -999.f};
                if (o2::base::Propagator::Instance()->propagateToDCABxByBz({primVtxRefitted.getX(), primVtxRefitted.getY(), primVtxRefitted.getZ()}, trackPar, 2.f, o2::base::Propagator::MatCorrType::USEMatCorrNONE, &dcaInfo)) {
                  dcaXYPvRefit = dcaInfo[0];
                }
              }
            }
          }
        }
        return dcaXYPvRefit * track.sign() * bFieldSign;
      };
      if (track.tpcNSigmaEl() < nSigmaTpcHadronMax && track.tpcNSigmaEl() > nSigmaTpcHadronMin) {
        histos.fill(HIST("hDcaSignedHadron"), track.pt(), dcaSigned);
        if constexpr (DoPvRefit) {
          histos.fill(HIST("hDcaSignedHadronPvRefit"), track.pt(), getDcaSignedPvRefit());
        }
      }

      if (std::abs(track.tofNSigmaEl()) > nSigmaTofMax) {
        continue;
      }
      histos.fill(HIST("hTofNSigPtQA"), track.pt(), track.tofNSigmaEl());
      histos.fill(HIST("hTpcNSigPAfterTofCut"), track.p(), track.tpcNSigmaEl());
      histos.fill(HIST("hTpcNSigPtAfterTofCut"), track.pt(), track.tpcNSigmaEl());
      if (track.tpcNSigmaEl() < nSigmaTpcHadronMax && track.tpcNSigmaEl() > nSigmaTpcHadronMin) {
        histos.fill(HIST("hDcaSignedHadronAfterTofCut"), track.pt(), dcaSigned);
        if constexpr (DoPvRefit) {
          histos.fill(HIST("hDcaSignedHadronAfterTofCutPvRefit"), track.pt(), getDcaSignedPvRefit());
        }
      }

      if (track.tpcNSigmaEl() < nSigmaTpcMin || track.tpcNSigmaEl() > nSigmaTpcMax) {
        continue;
      }

      if (track.tpcNSigmaEl() < nSigmaTpcHadronMax && track.tpcNSigmaEl() > nSigmaTpcHadronMin) {

        histos.fill(HIST("hDcaHadron"), track.pt(), track.dcaXY());
      }

      histos.fill(HIST("hTpcNSigPtQA"), track.pt(), track.tpcNSigmaEl());

      histos.fill(HIST("hDcaTrack"), track.pt(), track.dcaXY());
      histos.fill(HIST("hDcaSignedTrack"), track.pt(), dcaSigned, collision.numContrib());
      if constexpr (DoPvRefit) {
        histos.fill(HIST("hDcaSignedTrackPvRefit"), track.pt(), getDcaSignedPvRefit(), collision.numContrib());
      }
    }
  }

  void processData(soa::Filtered<MyCollisions>::iterator const& collision,
                   TracksEl const& tracks)
  {
    runData<false>(collision, tracks);
  }
  PROCESS_SWITCH(HfTaskSingleElectron, processData, "For real data", true);

  void processDataPvRefit(soa::Filtered<MyCollisions>::iterator const& collision,
                          TracksEl const& tracks,
                          AllTracks const&,
                          aod::BCsWithTimestamps const&)
  {
    runData<true>(collision, tracks);
  }
  PROCESS_SWITCH(HfTaskSingleElectron, processDataPvRefit, "For real data, also with the DCA to the PV refitted without the track", false);

  void processMc(soa::Filtered<MyCollisions>::iterator const& collision,
                 McTracksEl const& tracks,
                 aod::McParticles const&)
  {
    float const flagAnalysedEvt = 0.5;

    if (!collision.sel8()) {
      return;
    }

    if (collision.numContrib() < nContribMin) {
      return;
    }

    histos.fill(HIST("hVtxZ"), collision.posZ());
    histos.fill(HIST("hNEvents"), flagAnalysedEvt);

    for (const auto& track : tracks) {

      fillTrackQaPreCut(track);

      if (!trackSel(track)) {
        continue;
      }

      histos.fill(HIST("hEtaTrack"), track.eta());
      histos.fill(HIST("hPtTrack"), track.pt());

      histos.fill(HIST("hNClsTpcTrack"), track.tpcNClsCrossedRows());
      histos.fill(HIST("hCrossedRowsOverFindableTpcTrack"), track.tpcCrossedRowsOverFindableCls());
      histos.fill(HIST("hChi2TpcTrack"), track.tpcChi2NCl());
      histos.fill(HIST("hIbClsItsTrack"), track.itsNClsInnerBarrel());
      histos.fill(HIST("hDcaXYTrack"), track.dcaXY());
      histos.fill(HIST("hDcaZTrack"), track.dcaZ());

      histos.fill(HIST("hTofNSigPt"), track.pt(), track.tofNSigmaEl());
      histos.fill(HIST("hTpcNSigPt"), track.pt(), track.tpcNSigmaEl());
      fillTpcNSigPtRec(track, false);

      int mpdg{};   // electron source pdg code
      double mpt{}; // electron source pt
      int const source = getElecSource(track, mpt, mpdg);

      if (source == DirectBeauty || source == BeautyCharm) {
        histos.fill(HIST("hPdgB"), mpdg);
        histos.fill(HIST("hDcaBeauty"), track.pt(), track.dcaXY());
      }

      if (source == DirectCharm) {
        histos.fill(HIST("hPdgC"), mpdg);
        histos.fill(HIST("hDcaCharm"), track.pt(), track.dcaXY());
      }

      if (source >= GammaPi0 && source <= GammaSigma) {
        histos.fill(HIST("hPdgCo"), mpdg);
        histos.fill(HIST("hDcaConv"), track.pt(), track.dcaXY());
      }

      if (source >= Pi0 && source <= Sigma) {
        histos.fill(HIST("hPdgDa"), mpdg);
        histos.fill(HIST("hDcaDalitz"), track.pt(), track.dcaXY());
      }

      if (track.tpcNSigmaEl() < nSigmaTpcHadronMax && track.tpcNSigmaEl() > nSigmaTpcHadronMin)
        histos.fill(HIST("hDcaHadron"), track.pt(), track.dcaXY());

      if (std::abs(track.tofNSigmaEl()) > nSigmaTofMax) {
        continue;
      }
      histos.fill(HIST("hTofNSigPtQA"), track.pt(), track.tofNSigmaEl());
      histos.fill(HIST("hTpcNSigPtAfterTofCut"), track.pt(), track.tpcNSigmaEl());
      fillTpcNSigPtRec(track, true);

      if (track.tpcNSigmaEl() < nSigmaTpcMin || track.tpcNSigmaEl() > nSigmaTpcMax) {
        continue;
      }
      histos.fill(HIST("hTpcNSigPtQA"), track.pt(), track.tpcNSigmaEl());

      histos.fill(HIST("hDcaTrack"), track.pt(), track.dcaXY());
    }
  }
  PROCESS_SWITCH(HfTaskSingleElectron, processMc, "For MC simulations", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfTaskSingleElectron>(cfgc)};
}
