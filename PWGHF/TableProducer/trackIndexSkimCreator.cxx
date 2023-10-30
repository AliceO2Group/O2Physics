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

/// \file trackIndexSkimCreator.cxx
/// \brief Pre-selection of 2-prong and 3-prong secondary vertices of heavy-flavour decay candidates
///
/// \author Gian Michele Innocenti <gian.michele.innocenti@cern.ch>, CERN
/// \author Vít Kučera <vit.kucera@cern.ch>, CERN
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>, CERN
/// \author Mattia Faggin <mfaggin@cern.ch>, University and INFN Padova
/// \author Jinjoo Seo <jseo@cern.ch>, Inha University
/// \author Fabrizio Grosa <fgrosa@cern.ch>, CERN

#include <algorithm> // std::find
#include <iterator>  // std::distance
#include <string>    // std::string
#include <vector>    // std::vector

#include "CCDB/BasicCCDBManager.h"             // for PV refit
#include "DataFormatsParameters/GRPMagField.h" // for PV refit
#include "DataFormatsParameters/GRPObject.h"   // for PV refit
#include "DCAFitter/DCAFitterN.h"
#include "DetectorsBase/Propagator.h"     // for PV refit
#include "DetectorsVertexing/PVertexer.h" // for PV refit
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/V0.h"
#include "ReconstructionDataFormats/Vertex.h" // for PV refit

#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/Utils/utilsAnalysis.h"
#include "PWGHF/Utils/utilsBfieldCCDB.h"
#include "PWGHF/Utils/utilsDebugLcToK0sP.h"

using namespace o2;
using namespace o2::analysis;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;

// enum for candidate type
enum CandidateType {
  Cand2Prong = 0,
  Cand3Prong,
  CandV0bachelor,
  CandDstar,
  NCandidateTypes
};

// event rejection types
enum EventRejection {
  Trigger = 0,
  PositionX,
  PositionY,
  PositionZ,
  NContrib,
  Chi2,
  NEventRejection
};

// #define MY_DEBUG

#ifdef MY_DEBUG
using TracksWithSelAndDCA = soa::Join<aod::TracksWCovDcaExtra, aod::TrackSelection, aod::McTrackLabels>;
#define MY_DEBUG_MSG(condition, cmd) \
  if (condition) {                   \
    cmd;                             \
  }
#else
using TracksWithSelAndDCA = soa::Join<aod::TracksWCovDcaExtra, aod::TrackSelection>;
#define MY_DEBUG_MSG(condition, cmd)
#endif

/// Event selection
struct HfTrackIndexSkimCreatorTagSelCollisions {
  Produces<aod::HfSelCollision> rowSelectedCollision;

  Configurable<bool> fillHistograms{"fillHistograms", true, "fill histograms"};
  Configurable<double> xVertexMin{"xVertexMin", -100., "min. x of primary vertex [cm]"};
  Configurable<double> xVertexMax{"xVertexMax", 100., "max. x of primary vertex [cm]"};
  Configurable<double> yVertexMin{"yVertexMin", -100., "min. y of primary vertex [cm]"};
  Configurable<double> yVertexMax{"yVertexMax", 100., "max. y of primary vertex [cm]"};
  Configurable<double> zVertexMin{"zVertexMin", -100., "min. z of primary vertex [cm]"};
  Configurable<double> zVertexMax{"zVertexMax", 100., "max. z of primary vertex [cm]"};
  Configurable<int> nContribMin{"nContribMin", 0, "min. number of contributors to primary-vertex reconstruction"};
  Configurable<double> chi2Max{"chi2Max", 0., "max. chi^2 of primary-vertex reconstruction"};
  Configurable<std::string> triggerClassName{"triggerClassName", "kINT7", "trigger class"};
  Configurable<bool> useSel8Trigger{"useSel8Trigger", false, "use sel8 trigger condition, for Run3 studies"};

  ConfigurableAxis axisNumContributors{"axisNumContributors", {200, -0.5f, 199.5f}, "Number of PV contributors"};

  int triggerClass;

  HistogramRegistry registry{"registry"};

  void init(InitContext const&)
  {
    triggerClass = std::distance(aliasLabels, std::find(aliasLabels, aliasLabels + kNaliases, triggerClassName.value.data()));

    if (fillHistograms) {
      constexpr int kNBinsEvents = 2 + EventRejection::NEventRejection;
      std::string labels[kNBinsEvents];
      labels[0] = "processed";
      labels[1] = "selected";
      labels[2 + EventRejection::Trigger] = "rej. trigger";
      labels[2 + EventRejection::PositionX] = "rej. #it{x}";
      labels[2 + EventRejection::PositionY] = "rej. #it{y}";
      labels[2 + EventRejection::PositionZ] = "rej. #it{z}";
      labels[2 + EventRejection::NContrib] = "rej. # of contributors";
      labels[2 + EventRejection::Chi2] = "rej. #it{#chi}^{2}";
      AxisSpec axisEvents = {kNBinsEvents, 0.5, kNBinsEvents + 0.5, ""};
      registry.add("hEvents", "Events;;entries", HistType::kTH1F, {axisEvents});
      for (int iBin = 0; iBin < kNBinsEvents; iBin++) {
        registry.get<TH1>(HIST("hEvents"))->GetXaxis()->SetBinLabel(iBin + 1, labels[iBin].data());
      }
      // primary vertex histograms
      registry.add("hNContributors", "Number of PV contributors;entries", {HistType::kTH1F, {axisNumContributors}});
      registry.add("hPrimVtxX", "selected events;#it{x}_{prim. vtx.} (cm);entries", {HistType::kTH1F, {{200, -0.5, 0.5}}});
      registry.add("hPrimVtxY", "selected events;#it{y}_{prim. vtx.} (cm);entries", {HistType::kTH1F, {{200, -0.5, 0.5}}});
      registry.add("hPrimVtxZ", "selected events;#it{z}_{prim. vtx.} (cm);entries", {HistType::kTH1F, {{200, -20., 20.}}});
    }
  }

  /// Primary-vertex selection
  /// \param collision  collision table row
  /// \param statusCollision  bit map with rejection results
  template <typename Col>
  void selectVertex(const Col& collision, int& statusCollision)
  {
    if (fillHistograms) {
      registry.fill(HIST("hNContributors"), collision.numContrib());
    }

    // x position
    if (collision.posX() < xVertexMin || collision.posX() > xVertexMax) {
      SETBIT(statusCollision, EventRejection::PositionX);
      if (fillHistograms) {
        registry.fill(HIST("hEvents"), 3 + EventRejection::PositionX);
      }
    }
    // y position
    if (collision.posY() < yVertexMin || collision.posY() > yVertexMax) {
      SETBIT(statusCollision, EventRejection::PositionY);
      if (fillHistograms) {
        registry.fill(HIST("hEvents"), 3 + EventRejection::PositionY);
      }
    }
    // z position
    if (collision.posZ() < zVertexMin || collision.posZ() > zVertexMax) {
      SETBIT(statusCollision, EventRejection::PositionZ);
      if (fillHistograms) {
        registry.fill(HIST("hEvents"), 3 + EventRejection::PositionZ);
      }
    }
    // number of contributors
    if (collision.numContrib() < nContribMin) {
      SETBIT(statusCollision, EventRejection::NContrib);
      if (fillHistograms) {
        registry.fill(HIST("hEvents"), 3 + EventRejection::NContrib);
      }
    }
    // chi^2
    if (chi2Max > 0. && collision.chi2() > chi2Max) {
      SETBIT(statusCollision, EventRejection::Chi2);
      if (fillHistograms) {
        registry.fill(HIST("hEvents"), 3 + EventRejection::Chi2);
      }
    }
  }

  /// Event selection with trigger selection
  void processTrigSel(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision)
  {
    int statusCollision = 0;

    // processed events
    if (fillHistograms) {
      registry.fill(HIST("hEvents"), 1);
    }

    // trigger selection
    if ((!useSel8Trigger && !collision.alias_bit(triggerClass)) || (useSel8Trigger && !collision.sel8())) {
      SETBIT(statusCollision, EventRejection::Trigger);
      if (fillHistograms) {
        registry.fill(HIST("hEvents"), 3 + EventRejection::Trigger);
      }
    }

    // vertex selection
    selectVertex(collision, statusCollision);

    // TODO: add more event selection criteria

    // selected events
    if (fillHistograms && statusCollision == 0) {
      registry.fill(HIST("hEvents"), 2);
      registry.fill(HIST("hPrimVtxX"), collision.posX());
      registry.fill(HIST("hPrimVtxY"), collision.posY());
      registry.fill(HIST("hPrimVtxZ"), collision.posZ());
    }

    // fill table row
    rowSelectedCollision(statusCollision);
  }
  PROCESS_SWITCH(HfTrackIndexSkimCreatorTagSelCollisions, processTrigSel, "Use trigger selection", false);

  /// Event selection without trigger selection
  void processNoTrigSel(aod::Collision const& collision)
  {
    int statusCollision = 0;

    // processed events
    if (fillHistograms) {
      registry.fill(HIST("hEvents"), 1);
    }

    // vertex selection
    selectVertex(collision, statusCollision);

    // selected events
    if (fillHistograms && statusCollision == 0) {
      registry.fill(HIST("hEvents"), 2);
      registry.fill(HIST("hPrimVtxX"), collision.posX());
      registry.fill(HIST("hPrimVtxY"), collision.posY());
      registry.fill(HIST("hPrimVtxZ"), collision.posZ());
    }

    // fill table row
    rowSelectedCollision(statusCollision);
  }
  PROCESS_SWITCH(HfTrackIndexSkimCreatorTagSelCollisions, processNoTrigSel, "Do not use trigger selection", true);
};

/// Track selection
struct HfTrackIndexSkimCreatorTagSelTracks {
  SliceCache cache;
  Preslice<aod::Tracks> perCol = aod::track::collisionId;

  Produces<aod::HfSelTrack> rowSelectedTrack;
  Produces<aod::HfPvRefitTrack> tabPvRefitTrack;

  Configurable<bool> isRun2{"isRun2", false, "enable Run 2 or Run 3 GRP objects for magnetic field"};
  Configurable<bool> doPvRefit{"doPvRefit", false, "do PV refit excluding the considered track"};
  Configurable<bool> fillHistograms{"fillHistograms", true, "fill histograms"};
  Configurable<bool> debugPvRefit{"debugPvRefit", false, "debug lines for primary vertex refit"};
  // Configurable<double> bz{"bz", 5., "bz field"};
  // quality cut
  Configurable<bool> doCutQuality{"doCutQuality", true, "apply quality cuts"};
  Configurable<bool> useIsGlobalTrack{"useIsGlobalTrack", false, "check isGlobalTrack status for tracks, for Run3 studies"};
  Configurable<bool> useIsGlobalTrackWoDCA{"useIsGlobalTrackWoDCA", false, "check isGlobalTrackWoDCA status for tracks, for Run3 studies"};
  Configurable<int> tpcNClsFoundMin{"tpcNClsFoundMin", 70, "min. number of found TPC clusters"};
  // pT bins for single-track cuts
  Configurable<std::vector<double>> binsPtTrack{"binsPtTrack", std::vector<double>{hf_cuts_single_track::vecBinsPtTrack}, "track pT bin limits for 2-prong DCA XY pT-dependent cut"};
  // 2-prong cuts
  Configurable<double> ptMinTrack2Prong{"ptMinTrack2Prong", -1., "min. track pT for 2 prong candidate"};
  Configurable<LabeledArray<double>> cutsTrack2Prong{"cutsTrack2Prong", {hf_cuts_single_track::cutsTrack[0], hf_cuts_single_track::nBinsPtTrack, hf_cuts_single_track::nCutVarsTrack, hf_cuts_single_track::labelsPtTrack, hf_cuts_single_track::labelsCutVarTrack}, "Single-track selections per pT bin for 2-prong candidates"};
  Configurable<double> etaMinTrack2Prong{"etaMinTrack2Prong", -99999., "min. pseudorapidity for 2 prong candidate"};
  Configurable<double> etaMaxTrack2Prong{"etaMaxTrack2Prong", 4., "max. pseudorapidity for 2 prong candidate"};
  // 3-prong cuts
  Configurable<double> ptMinTrack3Prong{"ptMinTrack3Prong", -1., "min. track pT for 3 prong candidate"};
  Configurable<LabeledArray<double>> cutsTrack3Prong{"cutsTrack3Prong", {hf_cuts_single_track::cutsTrack[0], hf_cuts_single_track::nBinsPtTrack, hf_cuts_single_track::nCutVarsTrack, hf_cuts_single_track::labelsPtTrack, hf_cuts_single_track::labelsCutVarTrack}, "Single-track selections per pT bin for 3-prong candidates"};
  Configurable<double> etaMinTrack3Prong{"etaMinTrack3Prong", -99999., "min. pseudorapidity for 3 prong candidate"};
  Configurable<double> etaMaxTrack3Prong{"etaMaxTrack3Prong", 4., "max. pseudorapidity for 3 prong candidate"};
  // bachelor cuts (when using cascades)
  Configurable<double> ptMinTrackBach{"ptMinTrackBach", 0.3, "min. track pT for bachelor in cascade candidate"}; // 0.5 for PbPb 2015?
  Configurable<LabeledArray<double>> cutsTrackBach{"cutsTrackBach", {hf_cuts_single_track::cutsTrack[0], hf_cuts_single_track::nBinsPtTrack, hf_cuts_single_track::nCutVarsTrack, hf_cuts_single_track::labelsPtTrack, hf_cuts_single_track::labelsCutVarTrack}, "Single-track selections per pT bin for the bachelor of V0-bachelor candidates"};
  Configurable<double> etaMinTrackBach{"etaMinTrackBach", -99999., "min. pseudorapidity for bachelor in cascade candidate"};
  Configurable<double> etaMaxTrackBach{"etaMaxTrackBach", 0.8, "max. pseudorapidity for bachelor in cascade candidate"};
  // soft pion cuts for D*
  Configurable<double> ptMinSoftPionForDstar{"ptMinSoftPionForDstar", 0.05, "min. track pT for soft pion in D* candidate"};
  Configurable<double> ptMaxSoftPionForDstar{"ptMaxSoftPionForDstar", 2., "max. track pT for soft pion in D* candidate"};
  Configurable<double> etaMinSoftPionForDstar{"etaMinSoftPionForDstar", -99999., "min. pseudorapidity for soft pion in D* candidate"};
  Configurable<double> etaMaxSoftPionForDstar{"etaMaxSoftPionForDstar", 0.8, "max. pseudorapidity for soft pion in D* candidate"};
  Configurable<LabeledArray<double>> cutsTrackDstar{"cutsTrackDstar", {hf_cuts_single_track::cutsTrackPrimary[0], hf_cuts_single_track::nBinsPtTrack, hf_cuts_single_track::nCutVarsTrack, hf_cuts_single_track::labelsPtTrack, hf_cuts_single_track::labelsCutVarTrack}, "Single-track selections per pT bin for the soft pion of D* candidates"};
  Configurable<bool> useIsGlobalTrackForSoftPion{"useIsGlobalTrackForSoftPion", false, "check isGlobalTrack status for soft pion tracks"};
  Configurable<bool> useIsGlobalTrackWoDCAForSoftPion{"useIsGlobalTrackWoDCAForSoftPion", false, "check isGlobalTrackWoDCA status for soft pion tracks"};
  Configurable<bool> useIsQualityTrackITSForSoftPion{"useIsQualityTrackITSForSoftPion", true, "check qualityTracksITS status for soft pion tracks"};
  // CCDB
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> ccdbPathLut{"ccdbPathLut", "GLO/Param/MatLUT", "Path for LUT parametrization"};
  Configurable<std::string> ccdbPathGrp{"ccdbPathGrp", "GLO/GRP/GRP", "Path of the grp file (Run 2)"};
  Configurable<std::string> ccdbPathGrpMag{"ccdbPathGrpMag", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object (Run 3)"};
  // for debugging
#ifdef MY_DEBUG
  Configurable<std::vector<int>> indexK0Spos{"indexK0Spos", {729, 2866, 4754, 5457, 6891, 7824, 9243, 9810}, "indices of K0S positive daughters, for debug"};
  Configurable<std::vector<int>> indexK0Sneg{"indexK0Sneg", {730, 2867, 4755, 5458, 6892, 7825, 9244, 9811}, "indices of K0S negative daughters, for debug"};
  Configurable<std::vector<int>> indexProton{"indexProton", {717, 2810, 4393, 5442, 6769, 7793, 9002, 9789}, "indices of protons, for debug"};
#endif

  // Needed for PV refitting
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::base::MatLayerCylSet* lut;
  o2::base::Propagator::MatCorrType noMatCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;
  int runNumber;

  // single-track cuts
  static const int nCuts = 4;
  // array of 2-prong and 3-prong cuts
  std::array<LabeledArray<double>, CandidateType::NCandidateTypes> cutsSingleTrack;

  // QA of PV refit
  ConfigurableAxis axisPvRefitDeltaX{"axisPvRefitDeltaX", {1000, -0.5f, 0.5f}, "DeltaX binning PV refit"};
  ConfigurableAxis axisPvRefitDeltaY{"axisPvRefitDeltaY", {1000, -0.5f, 0.5f}, "DeltaY binning PV refit"};
  ConfigurableAxis axisPvRefitDeltaZ{"axisPvRefitDeltaZ", {1000, -0.5f, 0.5f}, "DeltaZ binning PV refit"};

  HistogramRegistry registry{"registry"};

  void init(InitContext const&)
  {
    cutsSingleTrack = {cutsTrack2Prong, cutsTrack3Prong, cutsTrackBach, cutsTrackDstar};

    if (etaMinTrack2Prong == -99999.) {
      etaMinTrack2Prong.value = -etaMaxTrack2Prong;
    }
    if (etaMinTrack3Prong == -99999.) {
      etaMinTrack3Prong.value = -etaMaxTrack3Prong;
    }
    if (etaMinTrackBach == -99999.) {
      etaMinTrackBach.value = -etaMaxTrackBach;
    }
    if (etaMinSoftPionForDstar == -99999.) {
      etaMinSoftPionForDstar.value = -etaMaxSoftPionForDstar;
    }

    if (fillHistograms) {
      // general tracks
      registry.add("hRejTracks", "Tracks;;entries", {HistType::kTH1F, {{20, 0.5, 20.5}}});
      registry.add("hPtNoCuts", "all tracks;#it{p}_{T}^{track} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}});

      // 2-prong histograms
      registry.add("hPtCuts2Prong", "tracks selected for 2-prong vertexing;#it{p}_{T}^{track} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}});
      registry.add("hDCAToPrimXYVsPtCuts2Prong", "tracks selected for 2-prong vertexing;#it{p}_{T}^{track} (GeV/#it{c});DCAxy to prim. vtx. (cm);entries", {HistType::kTH2F, {{360, 0., 36.}, {400, -2., 2.}}});
      registry.add("hEtaCuts2Prong", "tracks selected for 2-prong vertexing;#it{#eta};entries", {HistType::kTH1F, {{static_cast<int>(0.6 * (etaMaxTrack2Prong - etaMinTrack2Prong) * 100), -1.2 * etaMinTrack2Prong, 1.2 * etaMaxTrack2Prong}}});
      // 3-prong histograms
      registry.add("hPtCuts3Prong", "tracks selected for 3-prong vertexing;#it{p}_{T}^{track} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}});
      registry.add("hDCAToPrimXYVsPtCuts3Prong", "tracks selected for 3-prong vertexing;#it{p}_{T}^{track} (GeV/#it{c});DCAxy to prim. vtx. (cm);entries", {HistType::kTH2F, {{360, 0., 36.}, {400, -2., 2.}}});
      registry.add("hEtaCuts3Prong", "tracks selected for 3-prong vertexing;#it{#eta};entries", {HistType::kTH1F, {{static_cast<int>(0.6 * (etaMaxTrack3Prong - etaMinTrack3Prong) * 100), -1.2 * etaMinTrack3Prong, 1.2 * etaMaxTrack3Prong}}});
      // bachelor (for cascades) histograms
      registry.add("hPtCutsV0bachelor", "tracks selected for V0-bachelor vertexing;#it{p}_{T}^{track} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}});
      registry.add("hDCAToPrimXYVsPtCutsV0bachelor", "tracks selected for V0-bachelor vertexing;#it{p}_{T}^{track} (GeV/#it{c});DCAxy to prim. vtx. (cm);entries", {HistType::kTH2F, {{360, 0., 36.}, {400, -2., 2.}}});
      registry.add("hEtaCutsV0bachelor", "tracks selected for V0-bachelor vertexing;#it{#eta};entries", {HistType::kTH1F, {{static_cast<int>(0.6 * (etaMaxTrackBach - etaMinTrackBach) * 100), -1.2 * etaMinTrackBach, 1.2 * etaMaxTrackBach}}});
      // soft pion (for D*) histograms
      registry.add("hPtCutsSoftPionForDstar", "tracks selected for D* soft pion;#it{p}_{T}^{track} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}});
      registry.add("hDCAToPrimXYVsPtCutsSoftPionForDstar", "tracks selected for D* soft pion;#it{p}_{T}^{track} (GeV/#it{c});DCAxy to prim. vtx. (cm);entries", {HistType::kTH2F, {{360, 0., 36.}, {400, -2., 2.}}});
      registry.add("hEtaCutsSoftPionForDstar", "tracks selected for D* soft pion;#it{#eta};entries", {HistType::kTH1F, {{static_cast<int>(0.6 * (etaMaxSoftPionForDstar - etaMinSoftPionForDstar) * 100), -1.2 * etaMinSoftPionForDstar, 1.2 * etaMaxSoftPionForDstar}}});

      std::string cutNames[nCuts + 1] = {"selected", "rej pT", "rej eta", "rej track quality", "rej dca"};
      std::string candNames[CandidateType::NCandidateTypes] = {"2-prong", "3-prong", "bachelor", "dstar"};
      for (int iCandType = 0; iCandType < CandidateType::NCandidateTypes; iCandType++) {
        for (int iCut = 0; iCut < nCuts + 1; iCut++) {
          registry.get<TH1>(HIST("hRejTracks"))->GetXaxis()->SetBinLabel((nCuts + 1) * iCandType + iCut + 1, Form("%s %s", candNames[iCandType].data(), cutNames[iCut].data()));
        }
      }
    }

    // Needed for PV refitting
    if (doPvRefit) {
      if (fillHistograms) {
        AxisSpec axisCollisionX{100, -20.f, 20.f, "X (cm)"};
        AxisSpec axisCollisionY{100, -20.f, 20.f, "Y (cm)"};
        AxisSpec axisCollisionZ{100, -20.f, 20.f, "Z (cm)"};
        AxisSpec axisCollisionXOriginal{100, -2.f, 2.f, "X original PV (cm)"};
        AxisSpec axisCollisionYOriginal{100, -2.f, 2.f, "Y original PV (cm)"};
        AxisSpec axisCollisionZOriginal{100, -2.f, 2.f, "Z original PV (cm)"};
        AxisSpec axisCollisionNContrib{1000, 0, 1000, "Number of contributors"};
        AxisSpec axisCollisionDeltaX{axisPvRefitDeltaX, "#Delta x_{PV} (cm)"};
        AxisSpec axisCollisionDeltaY{axisPvRefitDeltaY, "#Delta y_{PV} (cm)"};
        AxisSpec axisCollisionDeltaZ{axisPvRefitDeltaZ, "#Delta z_{PV} (cm)"};

        registry.add("PvRefit/hVerticesPerTrack", "", kTH1F, {{3, 0.5f, 3.5f, ""}});
        registry.get<TH1>(HIST("PvRefit/hVerticesPerTrack"))->GetXaxis()->SetBinLabel(1, "All PV");
        registry.get<TH1>(HIST("PvRefit/hVerticesPerTrack"))->GetXaxis()->SetBinLabel(2, "PV refit doable");
        registry.get<TH1>(HIST("PvRefit/hVerticesPerTrack"))->GetXaxis()->SetBinLabel(3, "PV refit #chi^{2}!=-1");
        registry.add("PvRefit/hPvDeltaXvsNContrib", "", kTH2F, {axisCollisionNContrib, axisCollisionDeltaX});
        registry.add("PvRefit/hPvDeltaYvsNContrib", "", kTH2F, {axisCollisionNContrib, axisCollisionDeltaY});
        registry.add("PvRefit/hPvDeltaZvsNContrib", "", kTH2F, {axisCollisionNContrib, axisCollisionDeltaZ});
        registry.add("PvRefit/hChi2vsNContrib", "", kTH2F, {axisCollisionNContrib, {102, -1.5, 100.5, "#chi^{2} PV refit"}});
        registry.add("PvRefit/hPvRefitXChi2Minus1", "PV refit with #it{#chi}^{2}==#minus1", kTH2F, {axisCollisionX, axisCollisionXOriginal});
        registry.add("PvRefit/hPvRefitYChi2Minus1", "PV refit with #it{#chi}^{2}==#minus1", kTH2F, {axisCollisionY, axisCollisionYOriginal});
        registry.add("PvRefit/hPvRefitZChi2Minus1", "PV refit with #it{#chi}^{2}==#minus1", kTH2F, {axisCollisionZ, axisCollisionZOriginal});
        registry.add("PvRefit/hNContribPvRefitNotDoable", "N. contributors for PV refit not doable", kTH1F, {axisCollisionNContrib});
        registry.add("PvRefit/hNContribPvRefitChi2Minus1", "N. contributors original PV for PV refit #it{#chi}^{2}==#minus1", kTH1F, {axisCollisionNContrib});
      }

      ccdb->setURL(ccdbUrl);
      ccdb->setCaching(true);
      ccdb->setLocalObjectValidityChecking();

      lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>(ccdbPathLut));
      runNumber = 0;
    }
  }

  /// Single-track cuts for 2-prongs, 3-prongs, or cascades
  /// \param trackPt is the track pt
  /// \param dca is a 2-element array with dca in transverse and longitudinal directions
  /// \param candType is the flag to decide which cuts to be applied (either for 2-prong, 3-prong, or cascade decays)
  /// \return true if track passes all cuts
  bool isSelectedTrackDCA(const float& trackPt, const std::array<float, 2>& dca, const int candType)
  {
    auto pTBinTrack = findBin(binsPtTrack, trackPt);
    if (pTBinTrack == -1) {
      return false;
    }

    if (std::abs(dca[0]) < cutsSingleTrack[candType].get(pTBinTrack, "min_dcaxytoprimary")) {
      return false; // minimum DCAxy
    }
    if (std::abs(dca[0]) > cutsSingleTrack[candType].get(pTBinTrack, "max_dcaxytoprimary")) {
      return false; // maximum DCAxy
    }
    return true;
  }

  /// Single-track cuts for 2-prongs or 3-prongs
  /// \param hfTrack is a track
  /// \param trackPt is the track pt
  /// \param trackEta is the track eta
  /// \param dca is a 2-element array with dca in transverse and longitudinal directions
  /// \param statusProng is the selection flag
  template <typename T>
  void isSelectedTrack(const T& hfTrack, const float& trackPt, const float& trackEta, const std::array<float, 2>& dca, int& statusProng)
  {
    if (fillHistograms) {
      registry.fill(HIST("hPtNoCuts"), trackPt);
    }

    int iCut{2};
    // pT cut
    if (trackPt < ptMinTrack2Prong) {
      CLRBIT(statusProng, CandidateType::Cand2Prong); // set the nth bit to 0
      if (fillHistograms) {
        registry.fill(HIST("hRejTracks"), (nCuts + 1) * CandidateType::Cand2Prong + iCut);
      }
    }
    if (trackPt < ptMinTrack3Prong) {
      CLRBIT(statusProng, CandidateType::Cand3Prong);
      if (fillHistograms) {
        registry.fill(HIST("hRejTracks"), (nCuts + 1) * CandidateType::Cand3Prong + iCut);
      }
    }
    MY_DEBUG_MSG(isProtonFromLc, LOG(info) << "proton " << indexBach << " pt = " << trackPt << " (cut " << ptMinTrackBach << ")");

    if (trackPt < ptMinTrackBach) {
      CLRBIT(statusProng, CandidateType::CandV0bachelor);
      if (fillHistograms) {
        registry.fill(HIST("hRejTracks"), (nCuts + 1) * CandidateType::CandV0bachelor + iCut);
      }
    }
    if (trackPt < ptMinSoftPionForDstar || trackPt > ptMaxSoftPionForDstar) {
      CLRBIT(statusProng, CandidateType::CandDstar);
      if (fillHistograms) {
        registry.fill(HIST("hRejTracks"), (nCuts + 1) * CandidateType::CandDstar + iCut);
      }
    }

    iCut = 3;
    // eta cut
    if (TESTBIT(statusProng, CandidateType::Cand2Prong) && (trackEta > etaMaxTrack2Prong || trackEta < etaMinTrack2Prong)) {
      CLRBIT(statusProng, CandidateType::Cand2Prong);
      if (fillHistograms) {
        registry.fill(HIST("hRejTracks"), (nCuts + 1) * CandidateType::Cand2Prong + iCut);
      }
    }
    if (TESTBIT(statusProng, CandidateType::Cand3Prong) && (trackEta > etaMaxTrack3Prong || trackEta < etaMinTrack3Prong)) {
      CLRBIT(statusProng, CandidateType::Cand3Prong);
      if (fillHistograms) {
        registry.fill(HIST("hRejTracks"), (nCuts + 1) * CandidateType::Cand3Prong + iCut);
      }
    }
    MY_DEBUG_MSG(isProtonFromLc, LOG(info) << "proton " << indexBach << " eta = " << trackEta << " (cut " << etaMinTrackBach << " to " << etaMaxTrackBach << ")");

    if (TESTBIT(statusProng, CandidateType::CandV0bachelor) && (trackEta > etaMaxTrackBach || trackEta < etaMinTrackBach)) {
      CLRBIT(statusProng, CandidateType::CandV0bachelor);
      if (fillHistograms) {
        registry.fill(HIST("hRejTracks"), (nCuts + 1) * CandidateType::CandV0bachelor + iCut);
      }
    }

    if (TESTBIT(statusProng, CandidateType::CandDstar) && (trackEta > etaMaxSoftPionForDstar || trackEta < etaMinSoftPionForDstar)) {
      CLRBIT(statusProng, CandidateType::CandDstar);
      if (fillHistograms) {
        registry.fill(HIST("hRejTracks"), (nCuts + 1) * CandidateType::CandDstar + iCut);
      }
    }

    // quality cut
    MY_DEBUG_MSG(isProtonFromLc, LOG(info) << "proton " << indexBach << " tpcNClsFound = " << hfTrack.tpcNClsFound() << " (cut " << tpcNClsFoundMin.value << ")");

    iCut = 4;
    bool hasGoodQuality = true;
    if (doCutQuality.value && statusProng > 0) { // FIXME to make a more complete selection e.g track.flags() & o2::aod::track::TPCrefit && track.flags() & o2::aod::track::GoldenChi2 &&
      if (useIsGlobalTrack) {
        if (!hfTrack.isGlobalTrack()) {
          hasGoodQuality = false;
        }
      } else if (useIsGlobalTrackWoDCA) {
        if (!hfTrack.isGlobalTrackWoDCA()) {
          hasGoodQuality = false;
        }
      } else {
        UChar_t clustermap = hfTrack.itsClusterMap();
        if (!(hfTrack.tpcNClsFound() >= tpcNClsFoundMin.value && // is this the number of TPC clusters? It should not be used
              TESTBIT(hfTrack.flags(), o2::aod::track::ITSrefit) &&
              (TESTBIT(clustermap, 0) || TESTBIT(clustermap, 1)))) {
          hasGoodQuality = false;
        }
      }
      if (!hasGoodQuality) {
        MY_DEBUG_MSG(isProtonFromLc, LOG(info) << "proton " << indexBach << " did not pass clusters cut");
        for (int iCandType = 0; iCandType < CandidateType::NCandidateTypes; iCandType++) {
          if (iCandType == CandidateType::CandDstar) { // different quality criteria for D* soft pions
            continue;
          }
          CLRBIT(statusProng, iCandType);
          if (fillHistograms) {
            registry.fill(HIST("hRejTracks"), (nCuts + 1) * iCandType + iCut);
          }
        }
      }
    }

    // quality cut for soft pion
    hasGoodQuality = true;
    if (doCutQuality.value && TESTBIT(statusProng, CandidateType::CandDstar)) {
      if (useIsGlobalTrackForSoftPion) {
        if (!hfTrack.isGlobalTrack()) {
          hasGoodQuality = false;
        }
      } else if (useIsGlobalTrackWoDCAForSoftPion) {
        if (!hfTrack.isGlobalTrackWoDCA()) {
          hasGoodQuality = false;
        }
      } else if (useIsQualityTrackITSForSoftPion) {
        if (!hfTrack.isQualityTrackITS()) {
          hasGoodQuality = false;
        }
      } else { // selections for Run2 converted data
        UChar_t clustermap = hfTrack.itsClusterMap();
        if (!(TESTBIT(hfTrack.flags(), o2::aod::track::ITSrefit) && (TESTBIT(clustermap, 0) || TESTBIT(clustermap, 1)))) {
          hasGoodQuality = false;
        }
      }
      if (!hasGoodQuality) {
        CLRBIT(statusProng, CandidateType::CandDstar);
        if (fillHistograms) {
          registry.fill(HIST("hRejTracks"), (nCuts + 1) * CandidateType::CandDstar + iCut);
        }
      }
    }

    // DCA cut
    iCut = 5;
    if (statusProng > 0) {
      for (int iCandType = 0; iCandType < CandidateType::NCandidateTypes; ++iCandType) {
        if (TESTBIT(statusProng, iCandType) && !isSelectedTrackDCA(trackPt, dca, iCandType)) {
          CLRBIT(statusProng, iCandType);
          if (fillHistograms) {
            registry.fill(HIST("hRejTracks"), (nCuts + 1) * iCandType + iCut);
          }
        }
      }
    }
    MY_DEBUG_MSG(isProtonFromLc, LOG(info) << "statusProng = " << statusProng; printf("\n"));

    // fill histograms
    if (fillHistograms) {
      iCut = 1;
      if (TESTBIT(statusProng, CandidateType::Cand2Prong)) {
        registry.fill(HIST("hPtCuts2Prong"), trackPt);
        registry.fill(HIST("hEtaCuts2Prong"), trackEta);
        registry.fill(HIST("hDCAToPrimXYVsPtCuts2Prong"), trackPt, dca[0]);
        registry.fill(HIST("hRejTracks"), (nCuts + 1) * CandidateType::Cand2Prong + iCut);
      }
      if (TESTBIT(statusProng, CandidateType::Cand3Prong)) {
        registry.fill(HIST("hPtCuts3Prong"), trackPt);
        registry.fill(HIST("hEtaCuts3Prong"), trackEta);
        registry.fill(HIST("hDCAToPrimXYVsPtCuts3Prong"), trackPt, dca[0]);
        registry.fill(HIST("hRejTracks"), (nCuts + 1) * CandidateType::Cand3Prong + iCut);
      }
      if (TESTBIT(statusProng, CandidateType::CandV0bachelor)) {
        MY_DEBUG_MSG(isProtonFromLc, LOG(info) << "Will be kept: Proton from Lc " << indexBach);
        registry.fill(HIST("hPtCutsV0bachelor"), trackPt);
        registry.fill(HIST("hEtaCutsV0bachelor"), trackEta);
        registry.fill(HIST("hDCAToPrimXYVsPtCutsV0bachelor"), trackPt, dca[0]);
        registry.fill(HIST("hRejTracks"), (nCuts + 1) * CandidateType::CandV0bachelor + iCut);
      }
      if (TESTBIT(statusProng, CandidateType::CandDstar)) {
        registry.fill(HIST("hPtCutsSoftPionForDstar"), trackPt);
        registry.fill(HIST("hEtaCutsSoftPionForDstar"), trackEta);
        registry.fill(HIST("hDCAToPrimXYVsPtCutsSoftPionForDstar"), trackPt, dca[0]);
        registry.fill(HIST("hRejTracks"), (nCuts + 1) * CandidateType::CandDstar + iCut);
      }
    }
  }

  /// Method for the PV refit and DCA recalculation for tracks with a collision assigned
  /// \param collision is a collision
  /// \param bcWithTimeStamps is a table of bunch crossing joined with timestamps used to query the CCDB for B and material budget
  /// \param vecPvContributorGlobId is a vector containing the global ID of PV contributors for the current collision
  /// \param vecPvContributorTrackParCov is a vector containing the TrackParCov of PV contributors for the current collision
  /// \param myTrack is the track to be removed, if contributor, from the PV refit
  /// \param pvCoord is an array containing the coordinates of the refitted PV
  /// \param pvCovMatrix is an array containing the covariance matrix values of the refitted PV
  /// \param dcaXYdcaZ is an array containing the dcaXY and dcaZ of myTrack with respect to the refitted PV
  void performPvRefitTrack(aod::Collision const& collision,
                           aod::BCsWithTimestamps const& bcWithTimeStamps,
                           std::vector<int64_t> vecPvContributorGlobId,
                           std::vector<o2::track::TrackParCov> vecPvContributorTrackParCov,
                           TracksWithSelAndDCA::iterator const& myTrack,
                           std::array<float, 3>& pvCoord,
                           std::array<float, 6>& pvCovMatrix,
                           std::array<float, 2>& dcaXYdcaZ)
  {
    std::vector<bool> vecPvRefitContributorUsed(vecPvContributorGlobId.size(), true);

    /// Prepare the vertex refitting
    // set the magnetic field from CCDB
    auto bc = collision.bc_as<o2::aod::BCsWithTimestamps>();
    initCCDB(bc, runNumber, ccdb, isRun2 ? ccdbPathGrp : ccdbPathGrpMag, lut, isRun2);
    /*if (runNumber != bc.runNumber()) {

      if (isRun2) { // Run 2 GRP object
        o2::parameters::GRPObject* grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>(ccdbPathGrp, bc.timestamp());
        if (grpo != nullptr) {
          o2::base::Propagator::initFieldFromGRP(grpo);
          o2::base::Propagator::Instance()->setMatLUT(lut);
          LOGF(info, "Setting magnetic field to %d kG for run %d from its GRP CCDB object (type o2::parameters::GRPObject)", grpo->getNominalL3Field(), bc.runNumber());
        } else {
          LOGF(fatal, "Run 2 GRP object (type o2::parameters::GRPObject) is not available in CCDB for run=%d at timestamp=%llu", bc.runNumber(), bc.timestamp());
        }
      } else { // Run 3 GRP object
        o2::parameters::GRPMagField* grpo = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(ccdbPathGrpMag, bc.timestamp());
        if (grpo != nullptr) {
          o2::base::Propagator::initFieldFromGRP(grpo);
          o2::base::Propagator::Instance()->setMatLUT(lut);
          LOG(info) << "Setting magnetic field to current" << grpo->getL3Current() << " A for run" << bc.runNumber() << " from its GRP CCDB object (type o2::parameters::GRPMagField)";
        } else {
          LOGF(fatal, "Run 3 GRP object (type o2::parameters::GRPMagField) is not available in CCDB for run=%d at timestamp=%llu", bc.runNumber(), bc.timestamp());
        }
      }

      runNumber = bc.runNumber();
    }*/

    // build the VertexBase to initialize the vertexer
    o2::dataformats::VertexBase primVtx;
    primVtx.setX(collision.posX());
    primVtx.setY(collision.posY());
    primVtx.setZ(collision.posZ());
    primVtx.setCov(collision.covXX(), collision.covXY(), collision.covYY(), collision.covXZ(), collision.covYZ(), collision.covZZ());
    // configure PVertexer
    o2::vertexing::PVertexer vertexer;
    o2::conf::ConfigurableParam::updateFromString("pvertexer.useMeanVertexConstraint=false"); /// remove diamond constraint (let's keep it at the moment...)
    vertexer.init();
    bool pvRefitDoable = vertexer.prepareVertexRefit(vecPvContributorTrackParCov, primVtx);
    if (!pvRefitDoable) {
      LOG(info) << "Not enough tracks accepted for the refit";
      if (doPvRefit && fillHistograms) {
        registry.fill(HIST("PvRefit/hNContribPvRefitNotDoable"), collision.numContrib());
      }
    }
    if (debugPvRefit) {
      LOG(info) << "prepareVertexRefit = " << pvRefitDoable << " Ncontrib= " << vecPvContributorTrackParCov.size() << " Ntracks= " << collision.numContrib() << " Vtx= " << primVtx.asString();
    }

    if (fillHistograms) {
      registry.fill(HIST("PvRefit/hVerticesPerTrack"), 1);
      if (pvRefitDoable) {
        registry.fill(HIST("PvRefit/hVerticesPerTrack"), 2);
      }
    }
    /// PV refitting, if the tracks contributed to this at the beginning
    o2::dataformats::VertexBase primVtxBaseRecalc;
    bool recalcImpPar = false;
    if (doPvRefit && pvRefitDoable) {
      recalcImpPar = true;
      auto trackIterator = std::find(vecPvContributorGlobId.begin(), vecPvContributorGlobId.end(), myTrack.globalIndex()); /// track global index
      if (trackIterator != vecPvContributorGlobId.end()) {

        /// this track contributed to the PV fit: let's do the refit without it
        const int entry = std::distance(vecPvContributorGlobId.begin(), trackIterator);

        vecPvRefitContributorUsed[entry] = false; /// remove the track from the PV refitting

        auto primVtxRefitted = vertexer.refitVertex(vecPvRefitContributorUsed, primVtx); // vertex refit
        // LOG(info) << "refit " << cnt << "/" << ntr << " result = " << primVtxRefitted.asString();
        if (debugPvRefit) {
          LOG(info) << "refit for track with global index " << static_cast<int>(myTrack.globalIndex()) << " " << primVtxRefitted.asString();
        }
        if (primVtxRefitted.getChi2() < 0) {
          if (debugPvRefit) {
            LOG(info) << "---> Refitted vertex has bad chi2 = " << primVtxRefitted.getChi2();
          }
          if (fillHistograms) {
            registry.fill(HIST("PvRefit/hPvRefitXChi2Minus1"), primVtxRefitted.getX(), collision.posX());
            registry.fill(HIST("PvRefit/hPvRefitYChi2Minus1"), primVtxRefitted.getY(), collision.posY());
            registry.fill(HIST("PvRefit/hPvRefitZChi2Minus1"), primVtxRefitted.getZ(), collision.posZ());
            registry.fill(HIST("PvRefit/hNContribPvRefitChi2Minus1"), collision.numContrib());
          }
          recalcImpPar = false;
        } else if (fillHistograms) {
          registry.fill(HIST("PvRefit/hVerticesPerTrack"), 3);
        }
        if (fillHistograms) {
          registry.fill(HIST("PvRefit/hChi2vsNContrib"), primVtxRefitted.getNContributors(), primVtxRefitted.getChi2());
        }

        vecPvRefitContributorUsed[entry] = true; /// restore the track for the next PV refitting (probably not necessary here)

        if (recalcImpPar) {
          // fill the histograms for refitted PV with good Chi2
          const double deltaX = primVtx.getX() - primVtxRefitted.getX();
          const double deltaY = primVtx.getY() - primVtxRefitted.getY();
          const double deltaZ = primVtx.getZ() - primVtxRefitted.getZ();
          if (fillHistograms) {
            registry.fill(HIST("PvRefit/hPvDeltaXvsNContrib"), primVtxRefitted.getNContributors(), deltaX);
            registry.fill(HIST("PvRefit/hPvDeltaYvsNContrib"), primVtxRefitted.getNContributors(), deltaY);
            registry.fill(HIST("PvRefit/hPvDeltaZvsNContrib"), primVtxRefitted.getNContributors(), deltaZ);
          }

          // fill the newly calculated PV
          primVtxBaseRecalc.setX(primVtxRefitted.getX());
          primVtxBaseRecalc.setY(primVtxRefitted.getY());
          primVtxBaseRecalc.setZ(primVtxRefitted.getZ());
          primVtxBaseRecalc.setCov(primVtxRefitted.getSigmaX2(), primVtxRefitted.getSigmaXY(), primVtxRefitted.getSigmaY2(), primVtxRefitted.getSigmaXZ(), primVtxRefitted.getSigmaYZ(), primVtxRefitted.getSigmaZ2());
        }

        // cnt++;
      }
    } /// end 'if (doPvRefit && pvRefitDoable)'

    // updated value after PV recalculation
    if (recalcImpPar) {

      /// Track propagation to the PV refit considering also the material budget
      /// Mandatory for tracks updated at most only to the innermost ITS layer
      auto trackPar = getTrackPar(myTrack);
      o2::gpu::gpustd::array<float, 2> dcaInfo{-999., -999.};
      if (o2::base::Propagator::Instance()->propagateToDCABxByBz({primVtxBaseRecalc.getX(), primVtxBaseRecalc.getY(), primVtxBaseRecalc.getZ()}, trackPar, 2.f, noMatCorr, &dcaInfo)) {
        pvCoord[0] = primVtxBaseRecalc.getX();
        pvCoord[1] = primVtxBaseRecalc.getY();
        pvCoord[2] = primVtxBaseRecalc.getZ();
        pvCovMatrix[0] = primVtxBaseRecalc.getSigmaX2();
        pvCovMatrix[1] = primVtxBaseRecalc.getSigmaXY();
        pvCovMatrix[2] = primVtxBaseRecalc.getSigmaY2();
        pvCovMatrix[3] = primVtxBaseRecalc.getSigmaXZ();
        pvCovMatrix[4] = primVtxBaseRecalc.getSigmaYZ();
        pvCovMatrix[5] = primVtxBaseRecalc.getSigmaZ2();
        dcaXYdcaZ[0] = dcaInfo[0]; // [cm]
        dcaXYdcaZ[1] = dcaInfo[1]; // [cm]
        // TODO: add DCAxy and DCAz uncertainties?
      }

      /// Track propagation to the PV refit done only ia geometrical way
      /// Correct only if no further material budget is crossed, namely for tracks already propagated to the original PV
      // o2::dataformats::DCA impactParameter;
      // if (getTrackParCov(myTrack).propagateToDCA(primVtxBaseRecalc, o2::base::Propagator::Instance()->getNominalBz(), &impactParameter)) {
      //   if (debug) {
      //     LOG(info) << "===> nominal Bz: " << o2::base::Propagator::Instance()->getNominalBz();
      //   }
      //   pvCoord[0] = primVtxBaseRecalc.getX();
      //   pvCoord[1] = primVtxBaseRecalc.getY();
      //   pvCoord[2] = primVtxBaseRecalc.getZ();
      //   pvCovMatrix[0] = primVtxBaseRecalc.getSigmaX2();
      //   pvCovMatrix[1] = primVtxBaseRecalc.getSigmaXY();
      //   pvCovMatrix[2] = primVtxBaseRecalc.getSigmaY2();
      //   pvCovMatrix[3] = primVtxBaseRecalc.getSigmaXZ();
      //   pvCovMatrix[4] = primVtxBaseRecalc.getSigmaYZ();
      //   pvCovMatrix[5] = primVtxBaseRecalc.getSigmaZ2();
      //   dcaXYdcaZ[0] = impactParameter.getY(); // [cm]
      //   dcaXYdcaZ[1] = impactParameter.getZ();    // [cm]
      //   // TODO: add DCAxy and DCAz uncertainties?
      // }
    }

    return;
  } /// end of performPvRefitTrack function

  /// Partition for PV contributors

  Preslice<TrackAssoc> trackIndicesPerCollision = aod::track_association::collisionId;
  Partition<TracksWithSelAndDCA> pvContributors = ((aod::track::flags & (uint32_t)aod::track::PVContributor) == (uint32_t)aod::track::PVContributor);

  void process(aod::Collisions const& collisions,
               TrackAssoc const& trackIndices,
               TracksWithSelAndDCA const& tracks,
               aod::BCsWithTimestamps const& bcWithTimeStamps // for PV refit
#ifdef MY_DEBUG
               ,
               aod::McParticles& mcParticles
#endif
  )
  {
    rowSelectedTrack.reserve(tracks.size());

    // prepare vectors to cache quantities needed for PV refit
    std::vector<std::array<float, 2>> pvRefitDcaPerTrack{};
    std::vector<std::array<float, 3>> pvRefitPvCoordPerTrack{};
    std::vector<std::array<float, 6>> pvRefitPvCovMatrixPerTrack{};
    if (doPvRefit) {
      pvRefitDcaPerTrack.resize(tracks.size());
      pvRefitPvCoordPerTrack.resize(tracks.size());
      pvRefitPvCovMatrixPerTrack.resize(tracks.size());
      if (debugPvRefit) {
        LOG(info) << ">>> number of tracks: " << tracks.size();
        LOG(info) << ">>> number of collisions: " << collisions.size();
      }
      tabPvRefitTrack.reserve(tracks.size());
    }

    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      auto groupedTrackIndices = trackIndices.sliceBy(trackIndicesPerCollision, thisCollId);

      for (const auto& trackId : groupedTrackIndices) {
        int statusProng = BIT(CandidateType::NCandidateTypes) - 1; // all bits on
        auto track = trackId.track_as<TracksWithSelAndDCA>();
        auto trackIdx = track.globalIndex();
        float trackPt = track.pt();
        float trackEta = track.eta();

        std::array<float, 2> pvRefitDcaXYDcaZ{track.dcaXY(), track.dcaZ()};
        std::array<float, 3> pvRefitPvCoord{0.f, 0.f, 0.f};
        std::array<float, 6> pvRefitPvCovMatrix{1e10f, 1e10f, 1e10f, 1e10f, 1e10f, 1e10f};

#ifdef MY_DEBUG
        auto indexBach = track.mcParticleId();
        //      LOG(info) << "Checking label " << indexBach;
        bool isProtonFromLc = isProtonFromLcFunc(indexBach, indexProton);
#endif
        // PV refit and DCA recalculation only for tracks with an assigned collision
        if (doPvRefit && track.has_collision() && track.collisionId() == thisCollId && track.isPVContributor()) {
          pvRefitPvCoord = {collision.posX(), collision.posY(), collision.posZ()};
          pvRefitPvCovMatrix = {collision.covXX(), collision.covXY(), collision.covYY(), collision.covXZ(), collision.covYZ(), collision.covZZ()};

          /// retrieve PV contributors for the current collision
          std::vector<int64_t> vecPvContributorGlobId = {};
          std::vector<o2::track::TrackParCov> vecPvContributorTrackParCov = {};

          /// contributors for the current collision
          auto pvContrCollision = pvContributors->sliceByCached(aod::track::collisionId, thisCollId, cache);
          for (const auto& contributor : pvContrCollision) {
            vecPvContributorGlobId.push_back(contributor.globalIndex());
            vecPvContributorTrackParCov.push_back(getTrackParCov(contributor));
          }
          if (debugPvRefit) {
            LOG(info) << "### vecPvContributorGlobId.size()=" << vecPvContributorGlobId.size() << ", vecPvContributorTrackParCov.size()=" << vecPvContributorTrackParCov.size() << ", N. original contributors=" << collision.numContrib();
          }

          /// Perform the PV refit only for tracks with an assigned collision
          if (debugPvRefit) {
            LOG(info) << "[BEFORE performPvRefitTrack] track.collision().globalIndex(): " << collision.globalIndex();
          }
          performPvRefitTrack(collision, bcWithTimeStamps, vecPvContributorGlobId, vecPvContributorTrackParCov, track, pvRefitPvCoord, pvRefitPvCovMatrix, pvRefitDcaXYDcaZ);
          pvRefitDcaPerTrack[trackIdx] = pvRefitDcaXYDcaZ;
          pvRefitPvCoordPerTrack[trackIdx] = pvRefitPvCoord;
          pvRefitPvCovMatrixPerTrack[trackIdx] = pvRefitPvCovMatrix;
        } else if (track.collisionId() != thisCollId) {
          auto bc = collision.bc_as<o2::aod::BCsWithTimestamps>();
          initCCDB(bc, runNumber, ccdb, isRun2 ? ccdbPathGrp : ccdbPathGrpMag, lut, isRun2);
          auto trackPar = getTrackPar(track);
          o2::gpu::gpustd::array<float, 2> dcaInfo{-999., -999.};
          o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackPar, 2.f, noMatCorr, &dcaInfo);
          trackPt = trackPar.getPt();
          trackEta = trackPar.getEta();
          pvRefitDcaXYDcaZ[0] = dcaInfo[0];
          pvRefitDcaXYDcaZ[1] = dcaInfo[1];
        }

        MY_DEBUG_MSG(isProtonFromLc, LOG(info) << "\nWe found the proton " << indexBach);

        // bool cutStatus[CandidateType::NCandidateTypes][nCuts];
        // if (debug) {
        //   for (int iCandType = 0; iCandType < CandidateType::NCandidateTypes; iCandType++) {
        //     for (int iCut = 0; iCut < nCuts; iCut++) {
        //       cutStatus[iCandType][iCut] = true;
        //     }
        //   }
        // }

        isSelectedTrack(track, trackPt, trackEta, pvRefitDcaXYDcaZ, statusProng);
        rowSelectedTrack(statusProng);
      }
    }

    if (doPvRefit) { /// fill table with PV refit info (it has to be filled per track)
      for (auto iTrack{0u}; iTrack < tracks.size(); ++iTrack) {
        tabPvRefitTrack(pvRefitPvCoordPerTrack[iTrack][0], pvRefitPvCoordPerTrack[iTrack][1], pvRefitPvCoordPerTrack[iTrack][2],
                        pvRefitPvCovMatrixPerTrack[iTrack][0], pvRefitPvCovMatrixPerTrack[iTrack][1], pvRefitPvCovMatrixPerTrack[iTrack][2], pvRefitPvCovMatrixPerTrack[iTrack][3], pvRefitPvCovMatrixPerTrack[iTrack][4], pvRefitPvCovMatrixPerTrack[iTrack][5],
                        pvRefitDcaPerTrack[iTrack][0], pvRefitDcaPerTrack[iTrack][1]);
      }
    }
  }
};

//____________________________________________________________________________________________________________________________________________

/// Pre-selection of 2-prong and 3-prong secondary vertices
struct HfTrackIndexSkimCreator {
  Produces<aod::Hf2Prongs> rowTrackIndexProng2;
  Produces<aod::HfCutStatus2Prong> rowProng2CutStatus;
  Produces<aod::HfPvRefit2Prong> rowProng2PVrefit;
  Produces<aod::Hf3Prongs> rowTrackIndexProng3;
  Produces<aod::HfCutStatus3Prong> rowProng3CutStatus;
  Produces<aod::HfPvRefit3Prong> rowProng3PVrefit;
  Produces<aod::HfDstars> rowTrackIndexDstar;
  Produces<aod::HfCutStatusDstar> rowDstarCutStatus;
  Produces<aod::HfPvRefitDstar> rowDstarPVrefit;

  Configurable<bool> isRun2{"isRun2", false, "enable Run 2 or Run 3 GRP objects for magnetic field"};
  Configurable<int> do3Prong{"do3Prong", 0, "do 3 prong"};
  Configurable<bool> doDstar{"doDstar", false, "do D* candidates"};
  Configurable<bool> debug{"debug", false, "debug mode"};
  Configurable<bool> debugPvRefit{"debugPvRefit", false, "debug lines for primary vertex refit"};
  Configurable<bool> fillHistograms{"fillHistograms", true, "fill histograms"};
  ConfigurableAxis axisNumTracks{"axisNumTracks", {250, -0.5f, 249.5f}, "Number of tracks"};
  ConfigurableAxis axisNumCands{"axisNumCands", {200, -0.5f, 199.f}, "Number of candidates"};
  // Configurable<int> nCollsMax{"nCollsMax", -1, "Max collisions per file"}; //can be added to run over limited collisions per file - for tesing purposes
  // preselection
  Configurable<double> ptTolerance{"ptTolerance", 0.1, "pT tolerance in GeV/c for applying preselections before vertex reconstruction"};
  // vertexing
  // Configurable<double> bz{"bz", 5., "magnetic field kG"};
  Configurable<bool> propagateToPCA{"propagateToPCA", true, "create tracks version propagated to PCA"};
  Configurable<bool> useAbsDCA{"useAbsDCA", false, "Minimise abs. distance rather than chi2"};
  Configurable<bool> useWeightedFinalPCA{"useWeightedFinalPCA", false, "Recalculate vertex position using track covariances, effective only if useAbsDCA is true"};
  Configurable<double> maxR{"maxR", 200., "reject PCA's above this radius"};
  Configurable<double> maxDZIni{"maxDZIni", 4., "reject (if>0) PCA candidate if tracks DZ exceeds threshold"};
  Configurable<double> minParamChange{"minParamChange", 1.e-3, "stop iterations if largest change of any X is smaller than this"};
  Configurable<double> minRelChi2Change{"minRelChi2Change", 0.9, "stop iterations if chi2/chi2old > this"};
  // CCDB
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> ccdbPathLut{"ccdbPathLut", "GLO/Param/MatLUT", "Path for LUT parametrization"};
  Configurable<std::string> ccdbPathGrp{"ccdbPathGrp", "GLO/GRP/GRP", "Path of the grp file (Run 2)"};
  Configurable<std::string> ccdbPathGrpMag{"ccdbPathGrpMag", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object (Run 3)"};

  // D0 cuts
  Configurable<std::vector<double>> binsPtD0ToPiK{"binsPtD0ToPiK", std::vector<double>{hf_cuts_presel_2prong::vecBinsPt}, "pT bin limits for D0->piK pT-dependent cuts"};
  Configurable<LabeledArray<double>> cutsD0ToPiK{"cutsD0ToPiK", {hf_cuts_presel_2prong::cuts[0], hf_cuts_presel_2prong::nBinsPt, hf_cuts_presel_2prong::nCutVars, hf_cuts_presel_2prong::labelsPt, hf_cuts_presel_2prong::labelsCutVar}, "D0->piK selections per pT bin"};
  // Jpsi -> ee cuts
  Configurable<std::vector<double>> binsPtJpsiToEE{"binsPtJpsiToEE", std::vector<double>{hf_cuts_presel_2prong::vecBinsPt}, "pT bin limits for Jpsi->ee pT-dependent cuts"};
  Configurable<LabeledArray<double>> cutsJpsiToEE{"cutsJpsiToEE", {hf_cuts_presel_2prong::cuts[0], hf_cuts_presel_2prong::nBinsPt, hf_cuts_presel_2prong::nCutVars, hf_cuts_presel_2prong::labelsPt, hf_cuts_presel_2prong::labelsCutVar}, "Jpsi->ee selections per pT bin"};
  // Jpsi -> mumu cuts
  Configurable<std::vector<double>> binsPtJpsiToMuMu{"binsPtJpsiToMuMu", std::vector<double>{hf_cuts_presel_2prong::vecBinsPt}, "pT bin limits for Jpsi->mumu pT-dependent cuts"};
  Configurable<LabeledArray<double>> cutsJpsiToMuMu{"cutsJpsiToMuMu", {hf_cuts_presel_2prong::cuts[0], hf_cuts_presel_2prong::nBinsPt, hf_cuts_presel_2prong::nCutVars, hf_cuts_presel_2prong::labelsPt, hf_cuts_presel_2prong::labelsCutVar}, "Jpsi->mumu selections per pT bin"};
  // D+ cuts
  Configurable<std::vector<double>> binsPtDplusToPiKPi{"binsPtDplusToPiKPi", std::vector<double>{hf_cuts_presel_3prong::vecBinsPt}, "pT bin limits for D+->piKpi pT-dependent cuts"};
  Configurable<LabeledArray<double>> cutsDplusToPiKPi{"cutsDplusToPiKPi", {hf_cuts_presel_3prong::cuts[0], hf_cuts_presel_3prong::nBinsPt, hf_cuts_presel_3prong::nCutVars, hf_cuts_presel_3prong::labelsPt, hf_cuts_presel_3prong::labelsCutVar}, "D+->piKpi selections per pT bin"};
  // Ds+ cuts
  Configurable<std::vector<double>> binsPtDsToKKPi{"binsPtDsToKKPi", std::vector<double>{hf_cuts_presel_3prong::vecBinsPt}, "pT bin limits for Ds+->KKPi pT-dependent cuts"};
  Configurable<LabeledArray<double>> cutsDsToKKPi{"cutsDsToKKPi", {hf_cuts_presel_3prong::cuts[0], hf_cuts_presel_3prong::nBinsPt, hf_cuts_presel_3prong::nCutVars, hf_cuts_presel_3prong::labelsPt, hf_cuts_presel_3prong::labelsCutVar}, "Ds+->KKPi selections per pT bin"};
  // Lc+ cuts
  Configurable<std::vector<double>> binsPtLcToPKPi{"binsPtLcToPKPi", std::vector<double>{hf_cuts_presel_3prong::vecBinsPt}, "pT bin limits for Lc->pKpi pT-dependent cuts"};
  Configurable<LabeledArray<double>> cutsLcToPKPi{"cutsLcToPKPi", {hf_cuts_presel_3prong::cuts[0], hf_cuts_presel_3prong::nBinsPt, hf_cuts_presel_3prong::nCutVars, hf_cuts_presel_3prong::labelsPt, hf_cuts_presel_3prong::labelsCutVar}, "Lc->pKpi selections per pT bin"};
  // Xic+ cuts
  Configurable<std::vector<double>> binsPtXicToPKPi{"binsPtXicToPKPi", std::vector<double>{hf_cuts_presel_3prong::vecBinsPt}, "pT bin limits for Xic->pKpi pT-dependent cuts"};
  Configurable<LabeledArray<double>> cutsXicToPKPi{"cutsXicToPKPi", {hf_cuts_presel_3prong::cuts[0], hf_cuts_presel_3prong::nBinsPt, hf_cuts_presel_3prong::nCutVars, hf_cuts_presel_3prong::labelsPt, hf_cuts_presel_3prong::labelsCutVar}, "Xic->pKpi selections per pT bin"};
  // D*+ cuts
  Configurable<std::vector<double>> binsPtDstarToD0Pi{"binsPtDstarToD0Pi", std::vector<double>{hf_cuts_presel_dstar::vecBinsPt}, "pT bin limits for D*+->D0pi pT-dependent cuts"};
  Configurable<LabeledArray<double>> cutsDstarToD0Pi{"cutsDstarToD0Pi", {hf_cuts_presel_dstar::cuts[0], hf_cuts_presel_dstar::nBinsPt, hf_cuts_presel_dstar::nCutVars, hf_cuts_presel_dstar::labelsPt, hf_cuts_presel_dstar::labelsCutVar}, "D*+->D0pi selections per pT bin"};

  // Needed for PV refitting
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::base::MatLayerCylSet* lut;
  o2::base::Propagator::MatCorrType noMatCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;
  int runNumber;

  double massPi{0.};
  double massK{0.};
  double massProton{0.};
  double massElectron{0.};
  double massMuon{0.};
  double massDzero{0.};

  // int nColls{0}; //can be added to run over limited collisions per file - for tesing purposes

  static constexpr int kN2ProngDecays = hf_cand_2prong::DecayType::N2ProngDecays; // number of 2-prong hadron types
  static constexpr int kN3ProngDecays = hf_cand_3prong::DecayType::N3ProngDecays; // number of 3-prong hadron types
  static constexpr int kNCuts2Prong = 4;                                          // how many different selections are made on 2-prongs
  static constexpr int kNCuts3Prong = 4;                                          // how many different selections are made on 3-prongs
  static constexpr int kNCutsDstar = 3;                                           // how many different selections are made on Dstars
  std::array<std::array<std::array<double, 2>, 2>, kN2ProngDecays> arrMass2Prong;
  std::array<std::array<std::array<double, 3>, 2>, kN3ProngDecays> arrMass3Prong;
  // arrays of 2-prong and 3-prong cuts
  std::array<LabeledArray<double>, kN2ProngDecays> cut2Prong;
  std::array<std::vector<double>, kN2ProngDecays> pTBins2Prong;
  std::array<LabeledArray<double>, kN3ProngDecays> cut3Prong;
  std::array<std::vector<double>, kN3ProngDecays> pTBins3Prong;

  using SelectedCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::HfSelCollision>>;
  using TracksWithPVRefitAndDCA = soa::Join<aod::TracksWCovDcaExtra, aod::HfPvRefitTrack>;
  using FilteredTrackAssocSel = soa::Filtered<soa::Join<aod::TrackAssoc, aod::HfSelTrack>>;

  // filter collisions
  Filter filterSelectCollisions = (aod::hf_sel_collision::whyRejectColl == 0);

  // define slice of track indices per collisions
  Preslice<TracksWithPVRefitAndDCA> tracksPerCollision = aod::track::collisionId; // needed for PV refit

  // filter track indices
  Filter filterSelectTrackIds = (aod::hf_sel_track::isSelProng > 0);
  Preslice<FilteredTrackAssocSel> trackIndicesPerCollision = aod::track_association::collisionId;

  // FIXME
  // Partition<TracksWithPVRefitAndDCA> tracksPos = aod::track::signed1Pt > 0.f;
  // Partition<TracksWithPVRefitAndDCA> tracksNeg = aod::track::signed1Pt < 0.f;

  // QA of PV refit
  ConfigurableAxis axisPvRefitDeltaX{"axisPvRefitDeltaX", {1000, -0.5f, 0.5f}, "DeltaX binning PV refit"};
  ConfigurableAxis axisPvRefitDeltaY{"axisPvRefitDeltaY", {1000, -0.5f, 0.5f}, "DeltaY binning PV refit"};
  ConfigurableAxis axisPvRefitDeltaZ{"axisPvRefitDeltaZ", {1000, -0.5f, 0.5f}, "DeltaZ binning PV refit"};

  HistogramRegistry registry{"registry"};

  void init(InitContext const& context)
  {
    if (!doprocess2And3ProngsWithPvRefit && !doprocess2And3ProngsNoPvRefit) {
      return;
    }

    massPi = o2::analysis::pdg::MassPiPlus;
    massK = o2::analysis::pdg::MassKPlus;
    massProton = o2::analysis::pdg::MassProton;
    massElectron = o2::analysis::pdg::MassElectron;
    massMuon = o2::analysis::pdg::MassMuonPlus;
    massDzero = o2::analysis::pdg::MassD0;

    arrMass2Prong[hf_cand_2prong::DecayType::D0ToPiK] = std::array{std::array{massPi, massK},
                                                                   std::array{massK, massPi}};

    arrMass2Prong[hf_cand_2prong::DecayType::JpsiToEE] = std::array{std::array{massElectron, massElectron},
                                                                    std::array{massElectron, massElectron}};

    arrMass2Prong[hf_cand_2prong::DecayType::JpsiToMuMu] = std::array{std::array{massMuon, massMuon},
                                                                      std::array{massMuon, massMuon}};

    arrMass3Prong[hf_cand_3prong::DecayType::DplusToPiKPi] = std::array{std::array{massPi, massK, massPi},
                                                                        std::array{massPi, massK, massPi}};

    arrMass3Prong[hf_cand_3prong::DecayType::LcToPKPi] = std::array{std::array{massProton, massK, massPi},
                                                                    std::array{massPi, massK, massProton}};

    arrMass3Prong[hf_cand_3prong::DecayType::DsToKKPi] = std::array{std::array{massK, massK, massPi},
                                                                    std::array{massPi, massK, massK}};

    arrMass3Prong[hf_cand_3prong::DecayType::XicToPKPi] = std::array{std::array{massProton, massK, massPi},
                                                                     std::array{massPi, massK, massProton}};

    // cuts for 2-prong decays retrieved by json. the order must be then one in hf_cand_2prong::DecayType
    cut2Prong = {cutsD0ToPiK, cutsJpsiToEE, cutsJpsiToMuMu};
    pTBins2Prong = {binsPtD0ToPiK, binsPtJpsiToEE, binsPtJpsiToMuMu};
    // cuts for 3-prong decays retrieved by json. the order must be then one in hf_cand_3prong::DecayType
    cut3Prong = {cutsDplusToPiKPi, cutsLcToPKPi, cutsDsToKKPi, cutsXicToPKPi};
    pTBins3Prong = {binsPtDplusToPiKPi, binsPtLcToPKPi, binsPtDsToKKPi, binsPtXicToPKPi};

    if (fillHistograms) {
      registry.add("hNTracks", "Number of selected tracks;# of selected tracks;entries", {HistType::kTH1F, {axisNumTracks}});
      // 2-prong histograms
      registry.add("hVtx2ProngX", "2-prong candidates;#it{x}_{sec. vtx.} (cm);entries", {HistType::kTH1F, {{1000, -2., 2.}}});
      registry.add("hVtx2ProngY", "2-prong candidates;#it{y}_{sec. vtx.} (cm);entries", {HistType::kTH1F, {{1000, -2., 2.}}});
      registry.add("hVtx2ProngZ", "2-prong candidates;#it{z}_{sec. vtx.} (cm);entries", {HistType::kTH1F, {{1000, -20., 20.}}});
      registry.add("hNCand2Prong", "2-prong candidates preselected;# of candidates;entries", {HistType::kTH1F, {axisNumCands}});
      registry.add("hNCand2ProngVsNTracks", "2-prong candidates preselected;# of selected tracks;# of candidates;entries", {HistType::kTH2F, {axisNumTracks, axisNumCands}});
      registry.add("hMassD0ToPiK", "D^{0} candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{500, 0., 5.}}});
      registry.add("hMassJpsiToEE", "J/#psi candidates;inv. mass (e^{#plus} e^{#minus}) (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{500, 0., 5.}}});
      registry.add("hMassJpsiToMuMu", "J/#psi candidates;inv. mass (#mu^{#plus} #mu^{#minus}) (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{500, 0., 5.}}});
      // 3-prong histograms
      registry.add("hVtx3ProngX", "3-prong candidates;#it{x}_{sec. vtx.} (cm);entries", {HistType::kTH1F, {{1000, -2., 2.}}});
      registry.add("hVtx3ProngY", "3-prong candidates;#it{y}_{sec. vtx.} (cm);entries", {HistType::kTH1F, {{1000, -2., 2.}}});
      registry.add("hVtx3ProngZ", "3-prong candidates;#it{z}_{sec. vtx.} (cm);entries", {HistType::kTH1F, {{1000, -20., 20.}}});
      registry.add("hNCand3Prong", "3-prong candidates preselected;# of candidates;entries", {HistType::kTH1F, {axisNumCands}});
      registry.add("hNCand3ProngVsNTracks", "3-prong candidates preselected;# of selected tracks;# of candidates;entries", {HistType::kTH2F, {axisNumTracks, axisNumCands}});
      registry.add("hMassDPlusToPiKPi", "D^{#plus} candidates;inv. mass (#pi K #pi) (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{500, 0., 5.}}});
      registry.add("hMassLcToPKPi", "#Lambda_{c}^{#plus} candidates;inv. mass (p K #pi) (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{500, 0., 5.}}});
      registry.add("hMassDsToKKPi", "D_{s}^{#plus} candidates;inv. mass (K K #pi) (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{500, 0., 5.}}});
      registry.add("hMassXicToPKPi", "#Xi_{c}^{#plus} candidates;inv. mass (p K #pi) (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{500, 0., 5.}}});
      registry.add("hMassDstarToD0Pi", "D^{*#plus} candidates;inv. mass (K #pi #pi) - mass (K #pi) (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{500, 0.135, 0.185}}});

      // needed for PV refitting
      if (doprocess2And3ProngsWithPvRefit) {
        AxisSpec axisCollisionX{100, -20.f, 20.f, "X (cm)"};
        AxisSpec axisCollisionY{100, -20.f, 20.f, "Y (cm)"};
        AxisSpec axisCollisionZ{100, -20.f, 20.f, "Z (cm)"};
        AxisSpec axisCollisionXOriginal{1000, -20.f, 20.f, "X original PV (cm)"};
        AxisSpec axisCollisionYOriginal{1000, -20.f, 20.f, "Y original PV (cm)"};
        AxisSpec axisCollisionZOriginal{1000, -20.f, 20.f, "Z original PV (cm)"};
        AxisSpec axisCollisionNContrib{1000, 0, 1000, "Number of contributors"};
        AxisSpec axisCollisionDeltaX{axisPvRefitDeltaX, "#Delta x_{PV} (cm)"};
        AxisSpec axisCollisionDeltaY{axisPvRefitDeltaY, "#Delta y_{PV} (cm)"};
        AxisSpec axisCollisionDeltaZ{axisPvRefitDeltaZ, "#Delta z_{PV} (cm)"};
        registry.add("PvRefit/verticesPerCandidate", "", kTH1F, {{6, 0.5f, 6.5f, ""}});
        registry.get<TH1>(HIST("PvRefit/verticesPerCandidate"))->GetXaxis()->SetBinLabel(1, "All PV");
        registry.get<TH1>(HIST("PvRefit/verticesPerCandidate"))->GetXaxis()->SetBinLabel(2, "PV refit doable");
        registry.get<TH1>(HIST("PvRefit/verticesPerCandidate"))->GetXaxis()->SetBinLabel(3, "PV refit #chi^{2}!=-1");
        registry.get<TH1>(HIST("PvRefit/verticesPerCandidate"))->GetXaxis()->SetBinLabel(4, "PV refit #it{#chi}^{2}==#minus1");
        registry.get<TH1>(HIST("PvRefit/verticesPerCandidate"))->GetXaxis()->SetBinLabel(5, "1 daughter contr.");
        registry.get<TH1>(HIST("PvRefit/verticesPerCandidate"))->GetXaxis()->SetBinLabel(6, "no PV refit");
        registry.add("PvRefit/hPvDeltaXvsNContrib", "", kTH2F, {axisCollisionNContrib, axisCollisionDeltaX});
        registry.add("PvRefit/hPvDeltaYvsNContrib", "", kTH2F, {axisCollisionNContrib, axisCollisionDeltaY});
        registry.add("PvRefit/hPvDeltaZvsNContrib", "", kTH2F, {axisCollisionNContrib, axisCollisionDeltaZ});
        registry.add("PvRefit/hChi2vsNContrib", "", kTH2F, {axisCollisionNContrib, {102, -1.5, 100.5, "#chi^{2} PV refit"}});
        registry.add("PvRefit/hPvRefitXChi2Minus1", "PV refit with #it{#chi}^{2}==#minus1", kTH2F, {axisCollisionX, axisCollisionXOriginal});
        registry.add("PvRefit/hPvRefitYChi2Minus1", "PV refit with #it{#chi}^{2}==#minus1", kTH2F, {axisCollisionY, axisCollisionYOriginal});
        registry.add("PvRefit/hPvRefitZChi2Minus1", "PV refit with #it{#chi}^{2}==#minus1", kTH2F, {axisCollisionZ, axisCollisionZOriginal});
        registry.add("PvRefit/hNContribPvRefitNotDoable", "N. contributors for PV refit not doable", kTH1F, {axisCollisionNContrib});
        registry.add("PvRefit/hNContribPvRefitChi2Minus1", "N. contributors original PV for PV refit #it{#chi}^{2}==#minus1", kTH1F, {axisCollisionNContrib});
      }
    }

    ccdb->setURL(ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>(ccdbPathLut));
    runNumber = 0;
  }

  /// Method to perform selections for 2-prong candidates before vertex reconstruction
  /// \param pVecTrack0 is the momentum array of the first daughter track
  /// \param pVecTrack1 is the momentum array of the second daughter track
  /// \param dcaTrack0 is the dcaXY of the first daughter track
  /// \param dcaTrack1 is the dcaXY of the second daughter track
  /// \param cutStatus is a 2D array with outcome of each selection (filled only in debug mode)
  /// \param whichHypo information of the mass hypoteses that were selected
  /// \param isSelected ia s bitmap with selection outcome
  template <typename T1, typename T2, typename T3, typename T4>
  void is2ProngPreselected(T1 const& pVecTrack0, T1 const& pVecTrack1, T2 const& dcaTrack0, T2 const& dcaTrack1, T3& cutStatus, T4& whichHypo, int& isSelected)
  {
    /// FIXME: this would be better fixed by having a convention on the position of min and max in the 2D Array
    static std::vector<int> massMinIndex;
    static std::vector<int> massMaxIndex;
    static std::vector<int> d0d0Index;
    static auto cacheIndices = [](std::array<LabeledArray<double>, kN2ProngDecays>& cut2Prong, std::vector<int>& mins, std::vector<int>& maxs, std::vector<int>& d0d0) {
      mins.resize(cut2Prong.size());
      maxs.resize(cut2Prong.size());
      d0d0.resize(cut2Prong.size());
      for (size_t i = 0; i < cut2Prong.size(); ++i) {
        mins[i] = cut2Prong[i].colmap.find("massMin")->second;
        maxs[i] = cut2Prong[i].colmap.find("massMax")->second;
        d0d0[i] = cut2Prong[i].colmap.find("d0d0")->second;
      }
      return true;
    };
    cacheIndices(cut2Prong, massMinIndex, massMaxIndex, d0d0Index);

    auto arrMom = std::array{pVecTrack0, pVecTrack1};
    auto pT = RecoDecay::pt(pVecTrack0, pVecTrack1) + ptTolerance; // add tolerance because of no reco decay vertex

    for (int iDecay2P = 0; iDecay2P < kN2ProngDecays; iDecay2P++) {

      // pT
      auto pTBin = findBin(&pTBins2Prong[iDecay2P], pT);
      // return immediately if it is outside the defined pT bins
      if (pTBin == -1) {
        CLRBIT(isSelected, iDecay2P);
        if (debug) {
          cutStatus[iDecay2P][0] = false;
        }
        continue;
      }

      // invariant mass
      double massHypos[2];
      whichHypo[iDecay2P] = 3;
      double min2 = pow(cut2Prong[iDecay2P].get(pTBin, massMinIndex[iDecay2P]), 2);
      double max2 = pow(cut2Prong[iDecay2P].get(pTBin, massMaxIndex[iDecay2P]), 2);

      if ((debug || TESTBIT(isSelected, iDecay2P)) && cut2Prong[iDecay2P].get(pTBin, massMinIndex[iDecay2P]) >= 0. && cut2Prong[iDecay2P].get(pTBin, massMaxIndex[iDecay2P]) > 0.) {
        massHypos[0] = RecoDecay::m2(arrMom, arrMass2Prong[iDecay2P][0]);
        massHypos[1] = RecoDecay::m2(arrMom, arrMass2Prong[iDecay2P][1]);
        if (massHypos[0] < min2 || massHypos[0] >= max2) {
          whichHypo[iDecay2P] -= 1;
        }
        if (massHypos[1] < min2 || massHypos[1] >= max2) {
          whichHypo[iDecay2P] -= 2;
        }
        if (whichHypo[iDecay2P] == 0) {
          CLRBIT(isSelected, iDecay2P);
          if (debug) {
            cutStatus[iDecay2P][1] = false;
          }
        }
      }

      // imp. par. product cut
      if (debug || TESTBIT(isSelected, iDecay2P)) {
        auto impParProduct = dcaTrack0 * dcaTrack1;
        if (impParProduct > cut2Prong[iDecay2P].get(pTBin, d0d0Index[iDecay2P])) {
          CLRBIT(isSelected, iDecay2P);
          if (debug) {
            cutStatus[iDecay2P][2] = false;
          }
        }
      }
    }
  }

  /// Method to perform selections for 3-prong candidates before vertex reconstruction
  /// \param pVecTrack0 is the momentum array of the first daughter track
  /// \param pVecTrack1 is the momentum array of the second daughter track
  /// \param pVecTrack2 is the momentum array of the third daughter track
  /// \param cutStatus is a 2D array with outcome of each selection (filled only in debug mode)
  /// \param whichHypo information of the mass hypoteses that were selected
  /// \param isSelected ia s bitmap with selection outcome
  template <typename T1, typename T2, typename T3>
  void is3ProngPreselected(T1 const& pVecTrack0, T1 const& pVecTrack1, T1 const& pVecTrack2, T2& cutStatus, T3& whichHypo, int& isSelected)
  {
    /// FIXME: this would be better fixed by having a convention on the position of min and max in the 2D Array
    static std::vector<int> massMinIndex;
    static std::vector<int> massMaxIndex;
    static auto cacheIndices = [](std::array<LabeledArray<double>, kN3ProngDecays>& cut3Prong, std::vector<int>& mins, std::vector<int>& maxs) {
      mins.resize(cut3Prong.size());
      maxs.resize(cut3Prong.size());
      for (size_t iDecay3P = 0; iDecay3P < cut3Prong.size(); ++iDecay3P) {
        mins[iDecay3P] = cut3Prong[iDecay3P].colmap.find("massMin")->second;
        maxs[iDecay3P] = cut3Prong[iDecay3P].colmap.find("massMax")->second;
      }
      return true;
    };
    cacheIndices(cut3Prong, massMinIndex, massMaxIndex);

    auto arrMom = std::array{pVecTrack0, pVecTrack1, pVecTrack2};
    auto pT = RecoDecay::pt(pVecTrack0, pVecTrack1, pVecTrack2) + ptTolerance; // add tolerance because of no reco decay vertex

    for (int iDecay3P = 0; iDecay3P < kN3ProngDecays; iDecay3P++) {

      // pT
      auto pTBin = findBin(&pTBins3Prong[iDecay3P], pT);
      // return immediately if it is outside the defined pT bins
      if (pTBin == -1) {
        CLRBIT(isSelected, iDecay3P);
        if (debug) {
          cutStatus[iDecay3P][0] = false;
        }
        continue;
      }

      // invariant mass
      double massHypos[2];
      whichHypo[iDecay3P] = 3;
      double min2 = pow(cut3Prong[iDecay3P].get(pTBin, massMinIndex[iDecay3P]), 2);
      double max2 = pow(cut3Prong[iDecay3P].get(pTBin, massMaxIndex[iDecay3P]), 2);

      if ((debug || TESTBIT(isSelected, iDecay3P)) && cut3Prong[iDecay3P].get(pTBin, massMinIndex[iDecay3P]) >= 0. && cut3Prong[iDecay3P].get(pTBin, massMaxIndex[iDecay3P]) > 0.) { // no need to check isSelected but to avoid mistakes
        massHypos[0] = RecoDecay::m2(arrMom, arrMass3Prong[iDecay3P][0]);
        massHypos[1] = RecoDecay::m2(arrMom, arrMass3Prong[iDecay3P][1]);
        if (massHypos[0] < min2 || massHypos[0] >= max2) {
          whichHypo[iDecay3P] -= 1;
        }
        if (massHypos[1] < min2 || massHypos[1] >= max2) {
          whichHypo[iDecay3P] -= 2;
        }
        if (whichHypo[iDecay3P] == 0) {
          CLRBIT(isSelected, iDecay3P);
          if (debug) {
            cutStatus[iDecay3P][1] = false;
          }
        }
      }
    }
  }

  /// Method to perform selections for 2-prong candidates after vertex reconstruction
  /// \param pVecCand is the array for the candidate momentum after reconstruction of secondary vertex
  /// \param secVtx is the secondary vertex
  /// \param primVtx is the primary vertex
  /// \param cutStatus is a 2D array with outcome of each selection (filled only in debug mode)
  /// \param isSelected ia s bitmap with selection outcome
  template <typename T1, typename T2, typename T3, typename T4>
  void is2ProngSelected(const T1& pVecCand, const T2& secVtx, const T3& primVtx, T4& cutStatus, int& isSelected)
  {
    if (debug || isSelected > 0) {

      /// FIXME: this would be better fixed by having a convention on the position of cuts in the 2D Array
      static std::vector<int> cospIndex;
      static auto cacheIndices = [](std::array<LabeledArray<double>, kN2ProngDecays>& cut2Prong, std::vector<int>& cosp) {
        cosp.resize(cut2Prong.size());
        for (size_t iDecay2P = 0; iDecay2P < cut2Prong.size(); ++iDecay2P) {
          cosp[iDecay2P] = cut2Prong[iDecay2P].colmap.find("cosp")->second;
        }
        return true;
      };
      cacheIndices(cut2Prong, cospIndex);

      for (int iDecay2P = 0; iDecay2P < kN2ProngDecays; iDecay2P++) {

        // pT
        auto pTBin = findBin(&pTBins2Prong[iDecay2P], RecoDecay::pt(pVecCand));
        if (pTBin == -1) { // cut if it is outside the defined pT bins
          CLRBIT(isSelected, iDecay2P);
          if (debug) {
            cutStatus[iDecay2P][0] = false;
          }
          continue;
        }

        // cosp
        if (debug || TESTBIT(isSelected, iDecay2P)) {
          auto cpa = RecoDecay::cpa(primVtx, secVtx, pVecCand);
          if (cpa < cut2Prong[iDecay2P].get(pTBin, cospIndex[iDecay2P])) {
            CLRBIT(isSelected, iDecay2P);
            if (debug) {
              cutStatus[iDecay2P][3] = false;
            }
          }
        }
      }
    }
  }

  /// Method to perform selections for 3-prong candidates after vertex reconstruction
  /// \param pVecCand is the array for the candidate momentum after reconstruction of secondary vertex
  /// \param secVtx is the secondary vertex
  /// \param primVtx is the primary vertex
  /// \param cutStatus is a 2D array with outcome of each selection (filled only in debug mode)
  /// \param isSelected ia s bitmap with selection outcome
  template <typename T1, typename T2, typename T3, typename T4>
  void is3ProngSelected(const T1& pVecCand, const T2& secVtx, const T3& primVtx, T4& cutStatus, int& isSelected)
  {
    if (debug || isSelected > 0) {

      /// FIXME: this would be better fixed by having a convention on the position of cuts in the 2D Array
      static std::vector<int> cospIndex;
      static std::vector<int> decLenIndex;
      static auto cacheIndices = [](std::array<LabeledArray<double>, kN3ProngDecays>& cut3Prong, std::vector<int>& cosp, std::vector<int>& decL) {
        cosp.resize(cut3Prong.size());
        decL.resize(cut3Prong.size());
        for (size_t iDecay3P = 0; iDecay3P < cut3Prong.size(); ++iDecay3P) {
          cosp[iDecay3P] = cut3Prong[iDecay3P].colmap.find("cosp")->second;
          decL[iDecay3P] = cut3Prong[iDecay3P].colmap.find("decL")->second;
        }
        return true;
      };
      cacheIndices(cut3Prong, cospIndex, decLenIndex);

      for (int iDecay3P = 0; iDecay3P < kN3ProngDecays; iDecay3P++) {

        // pT
        auto pTBin = findBin(&pTBins3Prong[iDecay3P], RecoDecay::pt(pVecCand));
        if (pTBin == -1) { // cut if it is outside the defined pT bins
          CLRBIT(isSelected, iDecay3P);
          if (debug) {
            cutStatus[iDecay3P][0] = false;
          }
          continue;
        }

        // cosp
        if ((debug || TESTBIT(isSelected, iDecay3P))) {
          auto cpa = RecoDecay::cpa(primVtx, secVtx, pVecCand);
          if (cpa < cut3Prong[iDecay3P].get(pTBin, cospIndex[iDecay3P])) {
            CLRBIT(isSelected, iDecay3P);
            if (debug) {
              cutStatus[iDecay3P][2] = false;
            }
          }
        }

        // decay length
        if ((debug || TESTBIT(isSelected, iDecay3P))) {
          auto decayLength = RecoDecay::distance(primVtx, secVtx);
          if (decayLength < cut3Prong[iDecay3P].get(pTBin, decLenIndex[iDecay3P])) {
            CLRBIT(isSelected, iDecay3P);
            if (debug) {
              cutStatus[iDecay3P][3] = false;
            }
          }
        }
      }
    }
  }

  /// Method to perform selections for D* candidates before vertex reconstruction
  /// \param pVecTrack0 is the momentum array of the first daughter track (same charge)
  /// \param pVecTrack1 is the momentum array of the second daughter track (opposite charge)
  /// \param pVecTrack2 is the momentum array of the third daughter track (same charge)
  /// \param cutStatus is the cut status (filled only in case of debug)
  /// \param deltaMass is the M(Kpipi) - M(Kpi) value to be filled in the control histogram
  /// \return a bitmap with selection outcome
  template <typename T1, typename T2>
  uint8_t isDstarSelected(T1 const& pVecTrack0, T1 const& pVecTrack1, T1 const& pVecTrack2, uint8_t& cutStatus, T2& deltaMass)
  {
    uint8_t isSelected{1};
    /// FIXME: this would be better fixed by having a convention on the position of cuts in the 2D Array
    int deltaMassIndex, deltaMassD0Index;
    static auto cacheIndices = [](LabeledArray<double>& cutDstar, int& deltaMassIdx, int& deltaMassD0Idx) {
      deltaMassIdx = cutDstar.colmap.find("deltaMassMax")->second;
      deltaMassD0Idx = cutDstar.colmap.find("deltaMassD0")->second;
      return true;
    };
    cacheIndices(cutsDstarToD0Pi.value, deltaMassIndex, deltaMassD0Index);

    auto arrMom = std::array{pVecTrack0, pVecTrack1, pVecTrack2};
    auto arrMomD0 = std::array{pVecTrack0, pVecTrack1};
    auto pT = RecoDecay::pt(pVecTrack0, pVecTrack1, pVecTrack2) + ptTolerance; // add tolerance because of no reco decay vertex

    // pT
    auto pTBin = findBin(binsPtDstarToD0Pi, pT);
    // return immediately if it is outside the defined pT bins
    if (pTBin == -1) {
      isSelected = 0;
      if (debug) {
        CLRBIT(cutStatus, 0);
      } else {
        return isSelected;
      }
    }

    // D0 mass
    double deltaMassD0 = cutsDstarToD0Pi->get(pTBin, deltaMassD0Index);
    double invMassD0 = RecoDecay::m(arrMomD0, std::array{massPi, massK});
    if (std::abs(invMassD0 - massDzero) > deltaMassD0) {
      isSelected = 0;
      if (debug) {
        CLRBIT(cutStatus, 1);
      } else {
        return isSelected;
      }
    }

    // D*+ mass
    double maxDeltaMass = cutsDstarToD0Pi->get(pTBin, deltaMassIndex);
    double invMassDstar = RecoDecay::m(arrMom, std::array{massPi, massK, massPi});
    deltaMass = invMassDstar - invMassD0;
    if (deltaMass > maxDeltaMass) {
      isSelected = 0;
      if (debug) {
        CLRBIT(cutStatus, 1);
      } else {
        return isSelected;
      }
    }

    return isSelected;
  }

  /// Method for the PV refit excluding the candidate daughters
  /// \param collision is a collision
  /// \param bcWithTimeStamps is a table of bunch crossing joined with timestamps used to query the CCDB for B and material budget
  /// \param vecPvContributorGlobId is a vector containing the global ID of PV contributors for the current collision
  /// \param vecPvContributorTrackParCov is a vector containing the TrackParCov of PV contributors for the current collision
  /// \param vecCandPvContributorGlobId is a vector containing the global indices of daughter tracks that contributed to the original PV refit
  /// \param pvCoord is a vector where to store X, Y and Z values of refitted PV
  /// \param pvCovMatrix is a vector where to store the covariance matrix values of refitted PV
  void performPvRefitCandProngs(SelectedCollisions::iterator const& collision,
                                aod::BCsWithTimestamps const& bcWithTimeStamps,
                                std::vector<int64_t> vecPvContributorGlobId,
                                std::vector<o2::track::TrackParCov> vecPvContributorTrackParCov,
                                std::vector<int64_t> vecCandPvContributorGlobId,
                                std::array<float, 3>& pvCoord,
                                std::array<float, 6>& pvCovMatrix)
  {
    std::vector<bool> vecPvRefitContributorUsed(vecPvContributorGlobId.size(), true);

    /// Prepare the vertex refitting
    // set the magnetic field from CCDB
    auto bc = collision.bc_as<o2::aod::BCsWithTimestamps>();
    initCCDB(bc, runNumber, ccdb, isRun2 ? ccdbPathGrp : ccdbPathGrpMag, lut, isRun2);

    // build the VertexBase to initialize the vertexer
    o2::dataformats::VertexBase primVtx;
    primVtx.setX(collision.posX());
    primVtx.setY(collision.posY());
    primVtx.setZ(collision.posZ());
    primVtx.setCov(collision.covXX(), collision.covXY(), collision.covYY(), collision.covXZ(), collision.covYZ(), collision.covZZ());
    // configure PVertexer
    o2::vertexing::PVertexer vertexer;
    o2::conf::ConfigurableParam::updateFromString("pvertexer.useMeanVertexConstraint=false"); /// remove diamond constraint (let's keep it at the moment...)
    vertexer.init();
    bool pvRefitDoable = vertexer.prepareVertexRefit(vecPvContributorTrackParCov, primVtx);
    if (!pvRefitDoable) {
      LOG(info) << "Not enough tracks accepted for the refit";
      if (doprocess2And3ProngsWithPvRefit && fillHistograms) {
        registry.fill(HIST("PvRefit/hNContribPvRefitNotDoable"), collision.numContrib());
      }
    }
    if (debugPvRefit) {
      LOG(info) << "prepareVertexRefit = " << pvRefitDoable << " Ncontrib= " << vecPvContributorTrackParCov.size() << " Ntracks= " << collision.numContrib() << " Vtx= " << primVtx.asString();
    }

    /// PV refitting, if the tracks contributed to this at the beginning
    o2::dataformats::VertexBase primVtxBaseRecalc;
    bool recalcPvRefit = false;
    if (doprocess2And3ProngsWithPvRefit && pvRefitDoable) {
      if (fillHistograms) {
        registry.fill(HIST("PvRefit/verticesPerCandidate"), 2);
      }
      recalcPvRefit = true;
      int nCandContr = 0;
      for (uint64_t myGlobalID : vecCandPvContributorGlobId) {
        auto trackIterator = std::find(vecPvContributorGlobId.begin(), vecPvContributorGlobId.end(), myGlobalID); /// track global index
        if (trackIterator != vecPvContributorGlobId.end()) {
          /// this is a contributor, let's remove it for the PV refit
          const int entry = std::distance(vecPvContributorGlobId.begin(), trackIterator);
          vecPvRefitContributorUsed[entry] = false; /// remove the track from the PV refitting
          nCandContr++;
        }
      }

      /// do the PV refit excluding the candidate daughters that originally contributed to fit it
      if (debugPvRefit) {
        LOG(info) << "### PV refit after removing " << nCandContr << " tracks";
      }
      auto primVtxRefitted = vertexer.refitVertex(vecPvRefitContributorUsed, primVtx); // vertex refit
      // LOG(info) << "refit " << cnt << "/" << ntr << " result = " << primVtxRefitted.asString();
      // LOG(info) << "refit for track with global index " << static_cast<int>(myTrack.globalIndex()) << " " << primVtxRefitted.asString();
      if (primVtxRefitted.getChi2() < 0) {
        if (debugPvRefit) {
          LOG(info) << "---> Refitted vertex has bad chi2 = " << primVtxRefitted.getChi2();
        }
        if (fillHistograms) {
          registry.fill(HIST("PvRefit/verticesPerCandidate"), 4);
          registry.fill(HIST("PvRefit/hPvRefitXChi2Minus1"), primVtxRefitted.getX(), collision.posX());
          registry.fill(HIST("PvRefit/hPvRefitYChi2Minus1"), primVtxRefitted.getY(), collision.posY());
          registry.fill(HIST("PvRefit/hPvRefitZChi2Minus1"), primVtxRefitted.getZ(), collision.posZ());
          registry.fill(HIST("PvRefit/hNContribPvRefitChi2Minus1"), collision.numContrib());
        }
        recalcPvRefit = false;
      } else if (fillHistograms) {
        registry.fill(HIST("PvRefit/verticesPerCandidate"), 3);
      }
      if (fillHistograms) {
        registry.fill(HIST("PvRefit/hChi2vsNContrib"), primVtxRefitted.getNContributors(), primVtxRefitted.getChi2());
      }

      for (size_t i = 0; i < vecPvContributorGlobId.size(); i++) {
        vecPvRefitContributorUsed[i] = true; /// restore the tracks for the next PV refitting (probably not necessary here)
      }

      if (recalcPvRefit) {
        // fill the histograms for refitted PV with good Chi2
        const double deltaX = primVtx.getX() - primVtxRefitted.getX();
        const double deltaY = primVtx.getY() - primVtxRefitted.getY();
        const double deltaZ = primVtx.getZ() - primVtxRefitted.getZ();
        if (fillHistograms) {
          registry.fill(HIST("PvRefit/hPvDeltaXvsNContrib"), primVtxRefitted.getNContributors(), deltaX);
          registry.fill(HIST("PvRefit/hPvDeltaYvsNContrib"), primVtxRefitted.getNContributors(), deltaY);
          registry.fill(HIST("PvRefit/hPvDeltaZvsNContrib"), primVtxRefitted.getNContributors(), deltaZ);
        }

        // fill the newly calculated PV
        primVtxBaseRecalc.setX(primVtxRefitted.getX());
        primVtxBaseRecalc.setY(primVtxRefitted.getY());
        primVtxBaseRecalc.setZ(primVtxRefitted.getZ());
        primVtxBaseRecalc.setCov(primVtxRefitted.getSigmaX2(), primVtxRefitted.getSigmaXY(), primVtxRefitted.getSigmaY2(), primVtxRefitted.getSigmaXZ(), primVtxRefitted.getSigmaYZ(), primVtxRefitted.getSigmaZ2());

      } else {
        /// copy the original collision PV
        primVtxBaseRecalc.setX(primVtx.getX());
        primVtxBaseRecalc.setY(primVtx.getY());
        primVtxBaseRecalc.setZ(primVtx.getZ());
        primVtxBaseRecalc.setCov(primVtx.getSigmaX2(), primVtx.getSigmaXY(), primVtx.getSigmaY2(), primVtx.getSigmaXZ(), primVtx.getSigmaYZ(), primVtx.getSigmaZ2());
      }

      // fill the output
      pvCoord[0] = primVtxBaseRecalc.getX();
      pvCoord[1] = primVtxBaseRecalc.getY();
      pvCoord[2] = primVtxBaseRecalc.getZ();
      pvCovMatrix[0] = primVtxBaseRecalc.getSigmaX2();
      pvCovMatrix[1] = primVtxBaseRecalc.getSigmaXY();
      pvCovMatrix[2] = primVtxBaseRecalc.getSigmaY2();
      pvCovMatrix[3] = primVtxBaseRecalc.getSigmaXZ();
      pvCovMatrix[4] = primVtxBaseRecalc.getSigmaYZ();
      pvCovMatrix[5] = primVtxBaseRecalc.getSigmaZ2();

      // cnt++;

    } /// end 'if (doprocess2And3ProngsWithPvRefit && pvRefitDoable)'

    return;
  } /// end of performPvRefitCandProngs function

  template <bool doPvRefit = false, typename TTracks>
  void run2And3Prongs(SelectedCollisions const& collisions,
                      aod::BCsWithTimestamps const& bcWithTimeStamps,
                      FilteredTrackAssocSel const& trackIndices,
                      TTracks const& tracks)
  {

    // can be added to run over limited collisions per file - for tesing purposes
    /*
    if (nCollsMax > -1){
      if (nColls == nCollMax){
        return;
        //can be added to run over limited collisions per file - for tesing purposes
      }
      nColls++;
    }
    */

    for (const auto& collision : collisions) {

      /// retrieve PV contributors for the current collision
      std::vector<int64_t> vecPvContributorGlobId{};
      std::vector<o2::track::TrackParCov> vecPvContributorTrackParCov{};
      std::vector<bool> vecPvRefitContributorUsed{};
      if constexpr (doPvRefit) {
        auto groupedTracksUnfiltered = tracks.sliceBy(tracksPerCollision, collision.globalIndex());
        const int nTrk = groupedTracksUnfiltered.size();
        int nContrib = 0;
        int nNonContrib = 0;
        for (const auto& trackUnfiltered : groupedTracksUnfiltered) {
          if (!trackUnfiltered.isPVContributor()) {
            /// the track did not contribute to fit the primary vertex
            nNonContrib++;
            continue;
          } else {
            vecPvContributorGlobId.push_back(trackUnfiltered.globalIndex());
            vecPvContributorTrackParCov.push_back(getTrackParCov(trackUnfiltered));
            nContrib++;
            if (debugPvRefit) {
              LOG(info) << "---> a contributor! stuff saved";
              LOG(info) << "vec_contrib size: " << vecPvContributorTrackParCov.size() << ", nContrib: " << nContrib;
            }
          }
        }
        if (debugPvRefit) {
          LOG(info) << "===> nTrk: " << nTrk << ",   nContrib: " << nContrib << ",   nNonContrib: " << nNonContrib;
          if ((uint16_t)vecPvContributorTrackParCov.size() != collision.numContrib() || (uint16_t)nContrib != collision.numContrib()) {
            LOG(info) << "!!! Some problem here !!! vecPvContributorTrackParCov.size()= " << vecPvContributorTrackParCov.size() << ", nContrib=" << nContrib << ", collision.numContrib()" << collision.numContrib();
          }
        }
        vecPvRefitContributorUsed = std::vector<bool>(vecPvContributorGlobId.size(), true);
      }

      // auto centrality = collision.centV0M(); //FIXME add centrality when option for variations to the process function appears

      int n2ProngBit = BIT(kN2ProngDecays) - 1; // bit value for 2-prong candidates where each candidate is one bit and they are all set to 1
      int n3ProngBit = BIT(kN3ProngDecays) - 1; // bit value for 3-prong candidates where each candidate is one bit and they are all set to 1

      bool cutStatus2Prong[kN2ProngDecays][kNCuts2Prong];
      bool cutStatus3Prong[kN3ProngDecays][kNCuts3Prong];
      int nCutStatus2ProngBit = BIT(kNCuts2Prong) - 1; // bit value for selection status for each 2-prong candidate where each selection is one bit and they are all set to 1
      int nCutStatus3ProngBit = BIT(kNCuts3Prong) - 1; // bit value for selection status for each 3-prong candidate where each selection is one bit and they are all set to 1

      int whichHypo2Prong[kN2ProngDecays];
      int whichHypo3Prong[kN3ProngDecays];

      // set the magnetic field from CCDB
      auto bc = collision.bc_as<o2::aod::BCsWithTimestamps>();
      initCCDB(bc, runNumber, ccdb, isRun2 ? ccdbPathGrp : ccdbPathGrpMag, lut, isRun2);

      // 2-prong vertex fitter
      o2::vertexing::DCAFitterN<2> df2;
      df2.setBz(o2::base::Propagator::Instance()->getNominalBz());
      df2.setPropagateToPCA(propagateToPCA);
      df2.setMaxR(maxR);
      df2.setMaxDZIni(maxDZIni);
      df2.setMinParamChange(minParamChange);
      df2.setMinRelChi2Change(minRelChi2Change);
      df2.setUseAbsDCA(useAbsDCA);
      df2.setWeightedFinalPCA(useWeightedFinalPCA);

      // 3-prong vertex fitter
      o2::vertexing::DCAFitterN<3> df3;
      df3.setBz(o2::base::Propagator::Instance()->getNominalBz());
      df3.setPropagateToPCA(propagateToPCA);
      df3.setMaxR(maxR);
      df3.setMaxDZIni(maxDZIni);
      df3.setMinParamChange(minParamChange);
      df3.setMinRelChi2Change(minRelChi2Change);
      df3.setUseAbsDCA(useAbsDCA);
      df3.setWeightedFinalPCA(useWeightedFinalPCA);

      // used to calculate number of candidiates per event
      auto nCand2 = rowTrackIndexProng2.lastIndex();
      auto nCand3 = rowTrackIndexProng3.lastIndex();

      // if there isn't at least a positive and a negative track, continue immediately
      // if (tracksPos.size() < 1 || tracksNeg.size() < 1) {
      //  return;
      //}

      // first loop over positive tracks
      // for (auto trackPos1 = tracksPos.begin(); trackPos1 != tracksPos.end(); ++trackPos1) {

      auto thisCollId = collision.globalIndex();
      auto groupedTrackIndices = trackIndices.sliceBy(trackIndicesPerCollision, thisCollId);

      for (auto trackIndexPos1 = groupedTrackIndices.begin(); trackIndexPos1 != groupedTrackIndices.end(); ++trackIndexPos1) {
        auto trackPos1 = trackIndexPos1.template track_as<TTracks>();

        if (trackPos1.signed1Pt() < 0) {
          continue;
        }

        // retrieve the selection flag that corresponds to this collision
        auto isSelProngPos1 = trackIndexPos1.isSelProng();
        bool sel2ProngStatusPos = TESTBIT(isSelProngPos1, CandidateType::Cand2Prong);
        bool sel3ProngStatusPos1 = TESTBIT(isSelProngPos1, CandidateType::Cand3Prong);
        if (!sel2ProngStatusPos && !sel3ProngStatusPos1) {
          continue;
        }

        auto trackParVarPos1 = getTrackParCov(trackPos1);
        std::array<float, 3> pVecTrackPos1{trackPos1.px(), trackPos1.py(), trackPos1.pz()};
        o2::gpu::gpustd::array<float, 2> dcaInfoPos1{trackPos1.dcaXY(), trackPos1.dcaZ()};
        if (thisCollId != trackPos1.collisionId()) { // this is not the "default" collision for this track, we have to re-propagate it
          o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParVarPos1, 2.f, noMatCorr, &dcaInfoPos1);
          getPxPyPz(trackParVarPos1, pVecTrackPos1);
        }

        // first loop over negative tracks
        // for (auto trackNeg1 = tracksNeg.begin(); trackNeg1 != tracksNeg.end(); ++trackNeg1) {
        for (auto trackIndexNeg1 = groupedTrackIndices.begin(); trackIndexNeg1 != groupedTrackIndices.end(); ++trackIndexNeg1) {
          auto trackNeg1 = trackIndexNeg1.template track_as<TTracks>();
          if (trackNeg1.signed1Pt() > 0) {
            continue;
          }

          // retrieve the selection flag that corresponds to this collision
          auto isSelProngNeg1 = trackIndexNeg1.isSelProng();
          bool sel2ProngStatusNeg = TESTBIT(isSelProngNeg1, CandidateType::Cand2Prong);
          bool sel3ProngStatusNeg1 = TESTBIT(isSelProngNeg1, CandidateType::Cand3Prong);
          if (!sel2ProngStatusNeg && !sel3ProngStatusNeg1) {
            continue;
          }

          auto trackParVarNeg1 = getTrackParCov(trackNeg1);
          std::array<float, 3> pVecTrackNeg1{trackNeg1.px(), trackNeg1.py(), trackNeg1.pz()};
          o2::gpu::gpustd::array<float, 2> dcaInfoNeg1{trackNeg1.dcaXY(), trackNeg1.dcaZ()};
          if (thisCollId != trackNeg1.collisionId()) { // this is not the "default" collision for this track, we have to re-propagate it
            o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParVarNeg1, 2.f, noMatCorr, &dcaInfoNeg1);
            getPxPyPz(trackParVarNeg1, pVecTrackNeg1);
          }

          int isSelected2ProngCand = n2ProngBit; // bitmap for checking status of two-prong candidates (1 is true, 0 is rejected)

          if (debug) {
            for (int iDecay2P = 0; iDecay2P < kN2ProngDecays; iDecay2P++) {
              for (int iCut = 0; iCut < kNCuts2Prong; iCut++) {
                cutStatus2Prong[iDecay2P][iCut] = true;
              }
            }
          }

          // 2-prong vertex reconstruction
          if (sel2ProngStatusPos && sel2ProngStatusNeg) {

            // 2-prong preselections
            // TODO: in case of PV refit, the single-track DCA is calculated wrt two different PV vertices (only 1 track excluded)
            is2ProngPreselected(pVecTrackPos1, pVecTrackNeg1, dcaInfoPos1[0], dcaInfoNeg1[0], cutStatus2Prong, whichHypo2Prong, isSelected2ProngCand);

            // secondary vertex reconstruction and further 2-prong selections
            if (isSelected2ProngCand > 0 && df2.process(trackParVarPos1, trackParVarNeg1) > 0) { // should it be this or > 0 or are they equivalent
              // get secondary vertex
              const auto& secondaryVertex2 = df2.getPCACandidate();
              // get track momenta
              std::array<float, 3> pvec0;
              std::array<float, 3> pvec1;
              df2.getTrack(0).getPxPyPzGlo(pvec0);
              df2.getTrack(1).getPxPyPzGlo(pvec1);

              /// PV refit excluding the candidate daughters, if contributors
              std::array<float, 3> pvRefitCoord2Prong = {collision.posX(), collision.posY(), collision.posZ()}; /// initialize to the original PV
              std::array<float, 6> pvRefitCovMatrix2Prong = getPrimaryVertex(collision).getCov();               /// initialize to the original PV
              if constexpr (doPvRefit) {
                if (fillHistograms) {
                  registry.fill(HIST("PvRefit/verticesPerCandidate"), 1);
                }
                int nCandContr = 2;
                auto trackFirstIt = std::find(vecPvContributorGlobId.begin(), vecPvContributorGlobId.end(), trackPos1.globalIndex());
                auto trackSecondIt = std::find(vecPvContributorGlobId.begin(), vecPvContributorGlobId.end(), trackNeg1.globalIndex());
                bool isTrackFirstContr = true;
                bool isTrackSecondContr = true;
                if (trackFirstIt == vecPvContributorGlobId.end()) {
                  /// This track did not contribute to the original PV refit
                  if (debugPvRefit) {
                    LOG(info) << "--- [2 Prong] trackPos1 with globalIndex " << trackPos1.globalIndex() << " was not a PV contributor";
                  }
                  nCandContr--;
                  isTrackFirstContr = false;
                }
                if (trackSecondIt == vecPvContributorGlobId.end()) {
                  /// This track did not contribute to the original PV refit
                  if (debugPvRefit) {
                    LOG(info) << "--- [2 Prong] trackNeg1 with globalIndex " << trackNeg1.globalIndex() << " was not a PV contributor";
                  }
                  nCandContr--;
                  isTrackSecondContr = false;
                }
                if (nCandContr == 2) {
                  /// Both the daughter tracks were used for the original PV refit, let's refit it after excluding them
                  if (debugPvRefit) {
                    LOG(info) << "### [2 Prong] Calling performPvRefitCandProngs for HF 2 prong candidate";
                  }
                  performPvRefitCandProngs(collision, bcWithTimeStamps, vecPvContributorGlobId, vecPvContributorTrackParCov, {trackPos1.globalIndex(), trackNeg1.globalIndex()}, pvRefitCoord2Prong, pvRefitCovMatrix2Prong);
                } else if (nCandContr == 1) {
                  /// Only one daughter was a contributor, let's use then the PV recalculated by excluding only it
                  if (debugPvRefit) {
                    LOG(info) << "####### [2 Prong] nCandContr==" << nCandContr << " ---> just 1 contributor!";
                  }
                  if (fillHistograms) {
                    registry.fill(HIST("PvRefit/verticesPerCandidate"), 5);
                  }
                  if (isTrackFirstContr && !isTrackSecondContr) {
                    /// the first daughter is contributor, the second is not
                    pvRefitCoord2Prong = {trackPos1.pvRefitX(), trackPos1.pvRefitY(), trackPos1.pvRefitZ()};
                    pvRefitCovMatrix2Prong = {trackPos1.pvRefitSigmaX2(), trackPos1.pvRefitSigmaXY(), trackPos1.pvRefitSigmaY2(), trackPos1.pvRefitSigmaXZ(), trackPos1.pvRefitSigmaYZ(), trackPos1.pvRefitSigmaZ2()};
                  } else if (!isTrackFirstContr && isTrackSecondContr) {
                    ///  the second daughter is contributor, the first is not
                    pvRefitCoord2Prong = {trackNeg1.pvRefitX(), trackNeg1.pvRefitY(), trackNeg1.pvRefitZ()};
                    pvRefitCovMatrix2Prong = {trackNeg1.pvRefitSigmaX2(), trackNeg1.pvRefitSigmaXY(), trackNeg1.pvRefitSigmaY2(), trackNeg1.pvRefitSigmaXZ(), trackNeg1.pvRefitSigmaYZ(), trackNeg1.pvRefitSigmaZ2()};
                  }
                } else {
                  /// 0 contributors among the HF candidate daughters
                  if (fillHistograms) {
                    registry.fill(HIST("PvRefit/verticesPerCandidate"), 6);
                  }
                  if (debugPvRefit) {
                    LOG(info) << "####### [2 Prong] nCandContr==" << nCandContr << " ---> some of the candidate daughters did not contribute to the original PV fit, PV refit not redone";
                  }
                }
              }

              auto pVecCandProng2 = RecoDecay::pVec(pvec0, pvec1);
              // 2-prong selections after secondary vertex
              std::array<float, 3> pvCoord2Prong = {collision.posX(), collision.posY(), collision.posZ()};
              if constexpr (doPvRefit) {
                pvCoord2Prong[0] = pvRefitCoord2Prong[0];
                pvCoord2Prong[1] = pvRefitCoord2Prong[1];
                pvCoord2Prong[2] = pvRefitCoord2Prong[2];
              }
              is2ProngSelected(pVecCandProng2, secondaryVertex2, pvCoord2Prong, cutStatus2Prong, isSelected2ProngCand);

              if (isSelected2ProngCand > 0) {
                // fill table row
                rowTrackIndexProng2(thisCollId, trackPos1.globalIndex(), trackNeg1.globalIndex(), isSelected2ProngCand);

                if constexpr (doPvRefit) {
                  // fill table row with coordinates of PV refit
                  rowProng2PVrefit(pvRefitCoord2Prong[0], pvRefitCoord2Prong[1], pvRefitCoord2Prong[2],
                                   pvRefitCovMatrix2Prong[0], pvRefitCovMatrix2Prong[1], pvRefitCovMatrix2Prong[2], pvRefitCovMatrix2Prong[3], pvRefitCovMatrix2Prong[4], pvRefitCovMatrix2Prong[5]);
                }

                if (debug) {
                  int Prong2CutStatus[kN2ProngDecays];
                  for (int iDecay2P = 0; iDecay2P < kN2ProngDecays; iDecay2P++) {
                    Prong2CutStatus[iDecay2P] = nCutStatus2ProngBit;
                    for (int iCut = 0; iCut < kNCuts2Prong; iCut++) {
                      if (!cutStatus2Prong[iDecay2P][iCut]) {
                        CLRBIT(Prong2CutStatus[iDecay2P], iCut);
                      }
                    }
                  }
                  rowProng2CutStatus(Prong2CutStatus[0], Prong2CutStatus[1], Prong2CutStatus[2]); // FIXME when we can do this by looping over kN2ProngDecays
                }

                // fill histograms
                if (fillHistograms) {
                  registry.fill(HIST("hVtx2ProngX"), secondaryVertex2[0]);
                  registry.fill(HIST("hVtx2ProngY"), secondaryVertex2[1]);
                  registry.fill(HIST("hVtx2ProngZ"), secondaryVertex2[2]);
                  std::array<std::array<float, 3>, 2> arrMom = {pvec0, pvec1};
                  for (int iDecay2P = 0; iDecay2P < kN2ProngDecays; iDecay2P++) {
                    if (TESTBIT(isSelected2ProngCand, iDecay2P)) {
                      if (whichHypo2Prong[iDecay2P] == 1 || whichHypo2Prong[iDecay2P] == 3) {
                        auto mass2Prong = RecoDecay::m(arrMom, arrMass2Prong[iDecay2P][0]);
                        switch (iDecay2P) {
                          case hf_cand_2prong::DecayType::D0ToPiK:
                            registry.fill(HIST("hMassD0ToPiK"), mass2Prong);
                            break;
                          case hf_cand_2prong::DecayType::JpsiToEE:
                            registry.fill(HIST("hMassJpsiToEE"), mass2Prong);
                            break;
                          case hf_cand_2prong::DecayType::JpsiToMuMu:
                            registry.fill(HIST("hMassJpsiToMuMu"), mass2Prong);
                            break;
                        }
                      }
                      if (whichHypo2Prong[iDecay2P] >= 2) {
                        auto mass2Prong = RecoDecay::m(arrMom, arrMass2Prong[iDecay2P][1]);
                        if (iDecay2P == hf_cand_2prong::DecayType::D0ToPiK) {
                          registry.fill(HIST("hMassD0ToPiK"), mass2Prong);
                        }
                      }
                    }
                  }
                }
              }
            }
          }

          // 3-prong vertex and D* reconstruction
          if (do3Prong == 1 || doDstar) {
            if (!doDstar && (!sel3ProngStatusPos1 || !sel3ProngStatusNeg1)) {
              continue;
            }

            // second loop over positive tracks
            // for (auto trackPos2 = trackPos1 + 1; trackPos2 != tracksPos.end(); ++trackPos2) {
            for (auto trackIndexPos2 = trackIndexPos1 + 1; trackIndexPos2 != groupedTrackIndices.end(); ++trackIndexPos2) {
              auto trackPos2 = trackIndexPos2.template track_as<TTracks>();
              if (trackPos2.signed1Pt() < 0) {
                continue;
              }

              auto trackParVarPos2 = getTrackParCov(trackPos2);
              std::array<float, 3> pVecTrackPos2{trackPos2.px(), trackPos2.py(), trackPos2.pz()};
              o2::gpu::gpustd::array<float, 2> dcaInfoPos2{trackPos2.dcaXY(), trackPos2.dcaZ()};
              bool propagatedPos2{false};

              // first we build D*+ candidates if enabled
              auto isSelProngPos2 = trackIndexPos2.isSelProng();
              uint8_t isSelectedDstar{0};
              if (doDstar && TESTBIT(isSelected2ProngCand, hf_cand_2prong::DecayType::D0ToPiK) && (whichHypo2Prong[0] == 1 || whichHypo2Prong[0] == 3)) { // the 2-prong decay is compatible with a D0
                if (TESTBIT(isSelProngPos2, CandidateType::CandDstar)) {                                                                                  // compatible with a soft pion
                  if (thisCollId != trackPos2.collisionId()) {                                                                                            // this is not the "default" collision for this track, we have to re-propagate it
                    o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParVarPos2, 2.f, noMatCorr, &dcaInfoPos2);
                    getPxPyPz(trackParVarPos2, pVecTrackPos2);
                    propagatedPos2 = true;
                  }
                  uint8_t cutStatus{BIT(kNCutsDstar) - 1};
                  float deltaMass{-1.};
                  isSelectedDstar = isDstarSelected(pVecTrackPos1, pVecTrackNeg1, pVecTrackPos2, cutStatus, deltaMass); // we do not compute the D* decay vertex at this stage because we are not interested in applying topological selections
                  if (isSelectedDstar) {
                    rowTrackIndexDstar(thisCollId, trackPos2.globalIndex(), rowTrackIndexProng2.lastIndex());
                    if (fillHistograms) {
                      registry.fill(HIST("hMassDstarToD0Pi"), deltaMass);
                    }
                  }
                  if (debug) {
                    rowDstarCutStatus(cutStatus);
                  }
                }
              }

              // preselection of 3-prong candidates
              int isSelected3ProngCand = n3ProngBit;
              if (do3Prong && TESTBIT(isSelProngPos2, CandidateType::Cand3Prong) && (sel3ProngStatusPos1 && sel3ProngStatusNeg1)) {
                if (thisCollId != trackPos2.collisionId() && !propagatedPos2) { // this is not the "default" collision for this track and we still did not re-propagate it, we have to re-propagate it
                  o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParVarPos2, 2.f, noMatCorr, &dcaInfoPos2);
                  getPxPyPz(trackParVarPos2, pVecTrackPos2);
                }

                if (debug) {
                  for (int iDecay3P = 0; iDecay3P < kN3ProngDecays; iDecay3P++) {
                    for (int iCut = 0; iCut < kNCuts3Prong; iCut++) {
                      cutStatus3Prong[iDecay3P][iCut] = true;
                    }
                  }
                }

                // 3-prong preselections
                is3ProngPreselected(pVecTrackPos1, pVecTrackNeg1, pVecTrackPos2, cutStatus3Prong, whichHypo3Prong, isSelected3ProngCand);
                if (!debug && isSelected3ProngCand == 0) {
                  continue;
                }
              } else {
                isSelected3ProngCand = 0;
              }

              // if we did not preselected any D* or 3-prong candidate, continue
              if (!debug && (isSelectedDstar == 0 && isSelected3ProngCand == 0)) {
                continue;
              }

              /// PV refit excluding the candidate daughters, if contributors
              std::array<float, 3> pvRefitCoord3Prong2Pos1Neg = {collision.posX(), collision.posY(), collision.posZ()}; /// initialize to the original PV
              std::array<float, 6> pvRefitCovMatrix3Prong2Pos1Neg = getPrimaryVertex(collision).getCov();               /// initialize to the original PV
              if constexpr (doPvRefit) {
                if (fillHistograms) {
                  registry.fill(HIST("PvRefit/verticesPerCandidate"), 1);
                }
                int nCandContr = 3;
                auto trackFirstIt = std::find(vecPvContributorGlobId.begin(), vecPvContributorGlobId.end(), trackPos1.globalIndex());
                auto trackSecondIt = std::find(vecPvContributorGlobId.begin(), vecPvContributorGlobId.end(), trackNeg1.globalIndex());
                auto trackThirdIt = std::find(vecPvContributorGlobId.begin(), vecPvContributorGlobId.end(), trackPos2.globalIndex());
                bool isTrackFirstContr = true;
                bool isTrackSecondContr = true;
                bool isTrackThirdContr = true;
                if (trackFirstIt == vecPvContributorGlobId.end()) {
                  /// This track did not contribute to the original PV refit
                  if (debugPvRefit) {
                    LOG(info) << "--- [3 prong] trackPos1 with globalIndex " << trackPos1.globalIndex() << " was not a PV contributor";
                  }
                  nCandContr--;
                  isTrackFirstContr = false;
                }
                if (trackSecondIt == vecPvContributorGlobId.end()) {
                  /// This track did not contribute to the original PV refit
                  if (debugPvRefit) {
                    LOG(info) << "--- [3 prong] trackNeg1 with globalIndex " << trackNeg1.globalIndex() << " was not a PV contributor";
                  }
                  nCandContr--;
                  isTrackSecondContr = false;
                }
                if (trackThirdIt == vecPvContributorGlobId.end()) {
                  /// This track did not contribute to the original PV refit
                  if (debugPvRefit) {
                    LOG(info) << "--- [3 prong] trackPos2 with globalIndex " << trackPos2.globalIndex() << " was not a PV contributor";
                  }
                  nCandContr--;
                  isTrackThirdContr = false;
                }

                // Fill a vector with global ID of candidate daughters that are contributors
                std::vector<int64_t> vecCandPvContributorGlobId = {};
                if (isTrackFirstContr) {
                  vecCandPvContributorGlobId.push_back(trackPos1.globalIndex());
                }
                if (isTrackSecondContr) {
                  vecCandPvContributorGlobId.push_back(trackNeg1.globalIndex());
                }
                if (isTrackThirdContr) {
                  vecCandPvContributorGlobId.push_back(trackPos2.globalIndex());
                }

                if (nCandContr == 3 || nCandContr == 2) {
                  /// At least two of the daughter tracks were used for the original PV refit, let's refit it after excluding them
                  if (debugPvRefit) {
                    LOG(info) << "### [3 prong] Calling performPvRefitCandProngs for HF 3 prong candidate, removing " << nCandContr << " daughters";
                  }
                  performPvRefitCandProngs(collision, bcWithTimeStamps, vecPvContributorGlobId, vecPvContributorTrackParCov, vecCandPvContributorGlobId, pvRefitCoord3Prong2Pos1Neg, pvRefitCovMatrix3Prong2Pos1Neg);
                } else if (nCandContr == 1) {
                  /// Only one daughter was a contributor, let's use then the PV recalculated by excluding only it
                  if (debugPvRefit) {
                    LOG(info) << "####### [3 Prong] nCandContr==" << nCandContr << " ---> just 1 contributor!";
                  }
                  if (fillHistograms) {
                    registry.fill(HIST("PvRefit/verticesPerCandidate"), 5);
                  }
                  if (isTrackFirstContr && !isTrackSecondContr && !isTrackThirdContr) {
                    /// the first daughter is contributor, the second and the third are not
                    pvRefitCoord3Prong2Pos1Neg = {trackPos1.pvRefitX(), trackPos1.pvRefitY(), trackPos1.pvRefitZ()};
                    pvRefitCovMatrix3Prong2Pos1Neg = {trackPos1.pvRefitSigmaX2(), trackPos1.pvRefitSigmaXY(), trackPos1.pvRefitSigmaY2(), trackPos1.pvRefitSigmaXZ(), trackPos1.pvRefitSigmaYZ(), trackPos1.pvRefitSigmaZ2()};
                  } else if (!isTrackFirstContr && isTrackSecondContr && !isTrackThirdContr) {
                    /// the second daughter is contributor, the first and the third are not
                    pvRefitCoord3Prong2Pos1Neg = {trackNeg1.pvRefitX(), trackNeg1.pvRefitY(), trackNeg1.pvRefitZ()};
                    pvRefitCovMatrix3Prong2Pos1Neg = {trackNeg1.pvRefitSigmaX2(), trackNeg1.pvRefitSigmaXY(), trackNeg1.pvRefitSigmaY2(), trackNeg1.pvRefitSigmaXZ(), trackNeg1.pvRefitSigmaYZ(), trackNeg1.pvRefitSigmaZ2()};
                  } else if (!isTrackFirstContr && !isTrackSecondContr && isTrackThirdContr) {
                    /// the third daughter is contributor, the first and the second are not
                    pvRefitCoord3Prong2Pos1Neg = {trackPos2.pvRefitX(), trackPos2.pvRefitY(), trackPos2.pvRefitZ()};
                    pvRefitCovMatrix3Prong2Pos1Neg = {trackPos2.pvRefitSigmaX2(), trackPos2.pvRefitSigmaXY(), trackPos2.pvRefitSigmaY2(), trackPos2.pvRefitSigmaXZ(), trackPos2.pvRefitSigmaYZ(), trackPos2.pvRefitSigmaZ2()};
                  }
                } else {
                  /// 0 contributors among the HF candidate daughters
                  if (fillHistograms) {
                    registry.fill(HIST("PvRefit/verticesPerCandidate"), 6);
                  }
                  if (debugPvRefit) {
                    LOG(info) << "####### [3 prong] nCandContr==" << nCandContr << " ---> some of the candidate daughters did not contribute to the original PV fit, PV refit not redone";
                  }
                }
                if (isSelectedDstar) {
                  rowDstarPVrefit(pvRefitCoord3Prong2Pos1Neg[0], pvRefitCoord3Prong2Pos1Neg[1], pvRefitCoord3Prong2Pos1Neg[2],
                                  pvRefitCovMatrix3Prong2Pos1Neg[0], pvRefitCovMatrix3Prong2Pos1Neg[1], pvRefitCovMatrix3Prong2Pos1Neg[2], pvRefitCovMatrix3Prong2Pos1Neg[3], pvRefitCovMatrix3Prong2Pos1Neg[4], pvRefitCovMatrix3Prong2Pos1Neg[5]);
                }
              } // end of D*

              if (!debug && isSelected3ProngCand == 0) { // below only for 3-prong candidates
                continue;
              }

              // reconstruct the 3-prong secondary vertex
              if (df3.process(trackParVarPos1, trackParVarNeg1, trackParVarPos2) == 0) {
                continue;
              }
              // get secondary vertex
              const auto& secondaryVertex3 = df3.getPCACandidate();
              // get track momenta
              std::array<float, 3> pvec0;
              std::array<float, 3> pvec1;
              std::array<float, 3> pvec2;
              df3.getTrack(0).getPxPyPzGlo(pvec0);
              df3.getTrack(1).getPxPyPzGlo(pvec1);
              df3.getTrack(2).getPxPyPzGlo(pvec2);
              auto pVecCandProng3Pos = RecoDecay::pVec(pvec0, pvec1, pvec2);

              // 3-prong selections after secondary vertex
              is3ProngSelected(pVecCandProng3Pos, secondaryVertex3, pvRefitCoord3Prong2Pos1Neg, cutStatus3Prong, isSelected3ProngCand);
              if (!debug && isSelected3ProngCand == 0) {
                continue;
              }

              // fill table row
              rowTrackIndexProng3(thisCollId, trackPos1.globalIndex(), trackNeg1.globalIndex(), trackPos2.globalIndex(), isSelected3ProngCand);
              if constexpr (doPvRefit) {
                // fill table row of coordinates of PV refit
                rowProng3PVrefit(pvRefitCoord3Prong2Pos1Neg[0], pvRefitCoord3Prong2Pos1Neg[1], pvRefitCoord3Prong2Pos1Neg[2],
                                 pvRefitCovMatrix3Prong2Pos1Neg[0], pvRefitCovMatrix3Prong2Pos1Neg[1], pvRefitCovMatrix3Prong2Pos1Neg[2], pvRefitCovMatrix3Prong2Pos1Neg[3], pvRefitCovMatrix3Prong2Pos1Neg[4], pvRefitCovMatrix3Prong2Pos1Neg[5]);
              }

              if (debug) {
                int Prong3CutStatus[kN3ProngDecays];
                for (int iDecay3P = 0; iDecay3P < kN3ProngDecays; iDecay3P++) {
                  Prong3CutStatus[iDecay3P] = nCutStatus3ProngBit;
                  for (int iCut = 0; iCut < kNCuts3Prong; iCut++) {
                    if (!cutStatus3Prong[iDecay3P][iCut]) {
                      CLRBIT(Prong3CutStatus[iDecay3P], iCut);
                    }
                  }
                }
                rowProng3CutStatus(Prong3CutStatus[0], Prong3CutStatus[1], Prong3CutStatus[2], Prong3CutStatus[3]); // FIXME when we can do this by looping over kN3ProngDecays
              }

              // fill histograms
              if (fillHistograms) {
                registry.fill(HIST("hVtx3ProngX"), secondaryVertex3[0]);
                registry.fill(HIST("hVtx3ProngY"), secondaryVertex3[1]);
                registry.fill(HIST("hVtx3ProngZ"), secondaryVertex3[2]);
                std::array<std::array<float, 3>, 3> arr3Mom = {pvec0, pvec1, pvec2};
                for (int iDecay3P = 0; iDecay3P < kN3ProngDecays; iDecay3P++) {
                  if (TESTBIT(isSelected3ProngCand, iDecay3P)) {
                    if (whichHypo3Prong[iDecay3P] == 1 || whichHypo3Prong[iDecay3P] == 3) {
                      auto mass3Prong = RecoDecay::m(arr3Mom, arrMass3Prong[iDecay3P][0]);
                      switch (iDecay3P) {
                        case hf_cand_3prong::DecayType::DplusToPiKPi:
                          registry.fill(HIST("hMassDPlusToPiKPi"), mass3Prong);
                          break;
                        case hf_cand_3prong::DecayType::DsToKKPi:
                          registry.fill(HIST("hMassDsToKKPi"), mass3Prong);
                          break;
                        case hf_cand_3prong::DecayType::LcToPKPi:
                          registry.fill(HIST("hMassLcToPKPi"), mass3Prong);
                          break;
                        case hf_cand_3prong::DecayType::XicToPKPi:
                          registry.fill(HIST("hMassXicToPKPi"), mass3Prong);
                          break;
                      }
                    }
                    if (whichHypo3Prong[iDecay3P] >= 2) {
                      auto mass3Prong = RecoDecay::m(arr3Mom, arrMass3Prong[iDecay3P][1]);
                      switch (iDecay3P) {
                        case hf_cand_3prong::DecayType::DsToKKPi:
                          registry.fill(HIST("hMassDsToKKPi"), mass3Prong);
                          break;
                        case hf_cand_3prong::DecayType::LcToPKPi:
                          registry.fill(HIST("hMassLcToPKPi"), mass3Prong);
                          break;
                        case hf_cand_3prong::DecayType::XicToPKPi:
                          registry.fill(HIST("hMassXicToPKPi"), mass3Prong);
                          break;
                      }
                    }
                  }
                }
              }
            }

            // second loop over negative tracks
            // for (auto trackNeg2 = trackNeg1 + 1; trackNeg2 != tracksNeg.end(); ++trackNeg2) {
            for (auto trackIndexNeg2 = trackIndexNeg1 + 1; trackIndexNeg2 != groupedTrackIndices.end(); ++trackIndexNeg2) {
              auto trackNeg2 = trackIndexNeg2.template track_as<TTracks>();
              if (trackNeg2.signed1Pt() > 0) {
                continue;
              }

              auto trackParVarNeg2 = getTrackParCov(trackNeg2);
              std::array<float, 3> pVecTrackNeg2{trackNeg2.px(), trackNeg2.py(), trackNeg2.pz()};
              o2::gpu::gpustd::array<float, 2> dcaInfoNeg2{trackNeg2.dcaXY(), trackNeg2.dcaZ()};
              bool propagatedNeg2{false};

              // first we build D*+ candidates if enabled
              auto isSelProngNeg2 = trackIndexNeg2.isSelProng();
              uint8_t isSelectedDstar{0};
              if (doDstar && TESTBIT(isSelected2ProngCand, hf_cand_2prong::DecayType::D0ToPiK) && (whichHypo2Prong[0] >= 2)) { // the 2-prong decay is compatible with a D0bar
                if (TESTBIT(isSelProngNeg2, CandidateType::CandDstar)) {                                                       // compatible with a soft pion
                  if (thisCollId != trackNeg2.collisionId()) {                                                                 // this is not the "default" collision for this track, we have to re-propagate it
                    o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParVarNeg2, 2.f, noMatCorr, &dcaInfoNeg2);
                    getPxPyPz(trackParVarNeg2, pVecTrackNeg2);
                    propagatedNeg2 = true;
                  }
                  uint8_t cutStatus{BIT(kNCutsDstar) - 1};
                  float deltaMass{-1.};
                  isSelectedDstar = isDstarSelected(pVecTrackNeg1, pVecTrackPos1, pVecTrackNeg2, cutStatus, deltaMass); // we do not compute the D* decay vertex at this stage because we are not interested in applying topological selections
                  if (isSelectedDstar) {
                    rowTrackIndexDstar(thisCollId, trackNeg2.globalIndex(), rowTrackIndexProng2.lastIndex());
                    if (fillHistograms) {
                      registry.fill(HIST("hMassDstarToD0Pi"), deltaMass);
                    }
                  }
                  if (debug) {
                    rowDstarCutStatus(cutStatus);
                  }
                }
              }

              // preselection of 3-prong candidates
              int isSelected3ProngCand = n3ProngBit;
              if (do3Prong && TESTBIT(isSelProngNeg2, CandidateType::Cand3Prong) && (sel3ProngStatusPos1 && sel3ProngStatusNeg1)) {
                if (thisCollId != trackNeg2.collisionId() && !propagatedNeg2) { // this is not the "default" collision for this track and we still did not re-propagate it, we have to re-propagate it
                  o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParVarNeg2, 2.f, noMatCorr, &dcaInfoNeg2);
                  getPxPyPz(trackParVarNeg2, pVecTrackNeg2);
                }

                if (debug) {
                  for (int iDecay3P = 0; iDecay3P < kN3ProngDecays; iDecay3P++) {
                    for (int iCut = 0; iCut < kNCuts3Prong; iCut++) {
                      cutStatus3Prong[iDecay3P][iCut] = true;
                    }
                  }
                }

                // 3-prong preselections
                is3ProngPreselected(pVecTrackNeg1, pVecTrackPos1, pVecTrackNeg2, cutStatus3Prong, whichHypo3Prong, isSelected3ProngCand);
                if (!debug && isSelected3ProngCand == 0) {
                  continue;
                }
              } else {
                isSelected3ProngCand = 0;
              }

              // if we did not preselected any D* or 3-prong candidate, continue
              if (!debug && (isSelectedDstar == 0 && isSelected3ProngCand == 0)) {
                continue;
              }

              /// PV refit excluding the candidate daughters, if contributors
              std::array<float, 3> pvRefitCoord3Prong1Pos2Neg = {collision.posX(), collision.posY(), collision.posZ()}; /// initialize to the original PV
              std::array<float, 6> pvRefitCovMatrix3Prong1Pos2Neg = getPrimaryVertex(collision).getCov();               /// initialize to the original PV
              if constexpr (doPvRefit) {
                if (fillHistograms) {
                  registry.fill(HIST("PvRefit/verticesPerCandidate"), 1);
                }
                int nCandContr = 3;
                auto trackFirstIt = std::find(vecPvContributorGlobId.begin(), vecPvContributorGlobId.end(), trackPos1.globalIndex());
                auto trackSecondIt = std::find(vecPvContributorGlobId.begin(), vecPvContributorGlobId.end(), trackNeg1.globalIndex());
                auto trackThirdIt = std::find(vecPvContributorGlobId.begin(), vecPvContributorGlobId.end(), trackNeg2.globalIndex());
                bool isTrackFirstContr = true;
                bool isTrackSecondContr = true;
                bool isTrackThirdContr = true;
                if (trackFirstIt == vecPvContributorGlobId.end()) {
                  /// This track did not contribute to the original PV refit
                  if (debugPvRefit) {
                    LOG(info) << "--- [3 prong] trackPos1 with globalIndex " << trackPos1.globalIndex() << " was not a PV contributor";
                  }
                  nCandContr--;
                  isTrackFirstContr = false;
                }
                if (trackSecondIt == vecPvContributorGlobId.end()) {
                  /// This track did not contribute to the original PV refit
                  if (debugPvRefit) {
                    LOG(info) << "--- [3 prong] trackNeg1 with globalIndex " << trackNeg1.globalIndex() << " was not a PV contributor";
                  }
                  nCandContr--;
                  isTrackSecondContr = false;
                }
                if (trackThirdIt == vecPvContributorGlobId.end()) {
                  /// This track did not contribute to the original PV refit
                  if (debugPvRefit) {
                    LOG(info) << "--- [3 prong] trackNeg2 with globalIndex " << trackNeg2.globalIndex() << " was not a PV contributor";
                  }
                  nCandContr--;
                  isTrackThirdContr = false;
                }

                // Fill a vector with global ID of candidate daughters that are contributors
                std::vector<int64_t> vecCandPvContributorGlobId = {};
                if (isTrackFirstContr) {
                  vecCandPvContributorGlobId.push_back(trackPos1.globalIndex());
                }
                if (isTrackSecondContr) {
                  vecCandPvContributorGlobId.push_back(trackNeg1.globalIndex());
                }
                if (isTrackThirdContr) {
                  vecCandPvContributorGlobId.push_back(trackNeg2.globalIndex());
                }

                if (nCandContr == 3 || nCandContr == 2) {
                  /// At least two of the daughter tracks were used for the original PV refit, let's refit it after excluding them
                  if (debugPvRefit) {
                    LOG(info) << "### [3 prong] Calling performPvRefitCandProngs for HF 3 prong candidate, removing " << nCandContr << " daughters";
                  }
                  performPvRefitCandProngs(collision, bcWithTimeStamps, vecPvContributorGlobId, vecPvContributorTrackParCov, vecCandPvContributorGlobId, pvRefitCoord3Prong1Pos2Neg, pvRefitCovMatrix3Prong1Pos2Neg);
                } else if (nCandContr == 1) {
                  /// Only one daughter was a contributor, let's use then the PV recalculated by excluding only it
                  if (debugPvRefit) {
                    LOG(info) << "####### [3 Prong] nCandContr==" << nCandContr << " ---> just 1 contributor!";
                  }
                  if (fillHistograms) {
                    registry.fill(HIST("PvRefit/verticesPerCandidate"), 5);
                  }
                  if (isTrackFirstContr && !isTrackSecondContr && !isTrackThirdContr) {
                    /// the first daughter is contributor, the second and the third are not
                    pvRefitCoord3Prong1Pos2Neg = {trackPos1.pvRefitX(), trackPos1.pvRefitY(), trackPos1.pvRefitZ()};
                    pvRefitCovMatrix3Prong1Pos2Neg = {trackPos1.pvRefitSigmaX2(), trackPos1.pvRefitSigmaXY(), trackPos1.pvRefitSigmaY2(), trackPos1.pvRefitSigmaXZ(), trackPos1.pvRefitSigmaYZ(), trackPos1.pvRefitSigmaZ2()};
                  } else if (!isTrackFirstContr && isTrackSecondContr && !isTrackThirdContr) {
                    /// the second daughter is contributor, the first and the third are not
                    pvRefitCoord3Prong1Pos2Neg = {trackNeg1.pvRefitX(), trackNeg1.pvRefitY(), trackNeg1.pvRefitZ()};
                    pvRefitCovMatrix3Prong1Pos2Neg = {trackNeg1.pvRefitSigmaX2(), trackNeg1.pvRefitSigmaXY(), trackNeg1.pvRefitSigmaY2(), trackNeg1.pvRefitSigmaXZ(), trackNeg1.pvRefitSigmaYZ(), trackNeg1.pvRefitSigmaZ2()};
                  } else if (!isTrackFirstContr && !isTrackSecondContr && isTrackThirdContr) {
                    /// the third daughter is contributor, the first and the second are not
                    pvRefitCoord3Prong1Pos2Neg = {trackNeg2.pvRefitX(), trackNeg2.pvRefitY(), trackNeg2.pvRefitZ()};
                    pvRefitCovMatrix3Prong1Pos2Neg = {trackNeg2.pvRefitSigmaX2(), trackNeg2.pvRefitSigmaXY(), trackNeg2.pvRefitSigmaY2(), trackNeg2.pvRefitSigmaXZ(), trackNeg2.pvRefitSigmaYZ(), trackNeg2.pvRefitSigmaZ2()};
                  }
                } else {
                  /// 0 contributors among the HF candidate daughters
                  if (fillHistograms) {
                    registry.fill(HIST("PvRefit/verticesPerCandidate"), 6);
                  }
                  if (debugPvRefit) {
                    LOG(info) << "####### [3 prong] nCandContr==" << nCandContr << " ---> some of the candidate daughters did not contribute to the original PV fit, PV refit not redone";
                  }
                }
                if (isSelectedDstar) {
                  rowDstarPVrefit(pvRefitCoord3Prong1Pos2Neg[0], pvRefitCoord3Prong1Pos2Neg[1], pvRefitCoord3Prong1Pos2Neg[2],
                                  pvRefitCovMatrix3Prong1Pos2Neg[0], pvRefitCovMatrix3Prong1Pos2Neg[1], pvRefitCovMatrix3Prong1Pos2Neg[2], pvRefitCovMatrix3Prong1Pos2Neg[3], pvRefitCovMatrix3Prong1Pos2Neg[4], pvRefitCovMatrix3Prong1Pos2Neg[5]);
                }
              } // end of D*

              if (!debug && isSelected3ProngCand == 0) { // below only for 3-prong candidates
                continue;
              }

              // reconstruct the 3-prong secondary vertex
              if (df3.process(trackParVarNeg1, trackParVarPos1, trackParVarNeg2) == 0) {
                continue;
              }
              // get secondary vertex
              const auto& secondaryVertex3 = df3.getPCACandidate();
              // get track momenta
              std::array<float, 3> pvec0;
              std::array<float, 3> pvec1;
              std::array<float, 3> pvec2;
              df3.getTrack(0).getPxPyPzGlo(pvec0);
              df3.getTrack(1).getPxPyPzGlo(pvec1);
              df3.getTrack(2).getPxPyPzGlo(pvec2);
              auto pVecCandProng3Neg = RecoDecay::pVec(pvec0, pvec1, pvec2);

              // 3-prong selections after secondary vertex
              is3ProngSelected(pVecCandProng3Neg, secondaryVertex3, pvRefitCoord3Prong1Pos2Neg, cutStatus3Prong, isSelected3ProngCand);
              if (!debug && isSelected3ProngCand == 0) {
                continue;
              }

              // fill table row
              rowTrackIndexProng3(thisCollId, trackNeg1.globalIndex(), trackPos1.globalIndex(), trackNeg2.globalIndex(), isSelected3ProngCand);
              // fill table row of coordinates of PV refit
              if constexpr (doPvRefit) {
                rowProng3PVrefit(pvRefitCoord3Prong1Pos2Neg[0], pvRefitCoord3Prong1Pos2Neg[1], pvRefitCoord3Prong1Pos2Neg[2],
                                 pvRefitCovMatrix3Prong1Pos2Neg[0], pvRefitCovMatrix3Prong1Pos2Neg[1], pvRefitCovMatrix3Prong1Pos2Neg[2], pvRefitCovMatrix3Prong1Pos2Neg[3], pvRefitCovMatrix3Prong1Pos2Neg[4], pvRefitCovMatrix3Prong1Pos2Neg[5]);
              }

              if (debug) {
                int Prong3CutStatus[kN3ProngDecays];
                for (int iDecay3P = 0; iDecay3P < kN3ProngDecays; iDecay3P++) {
                  Prong3CutStatus[iDecay3P] = nCutStatus3ProngBit;
                  for (int iCut = 0; iCut < kNCuts3Prong; iCut++) {
                    if (!cutStatus3Prong[iDecay3P][iCut]) {
                      CLRBIT(Prong3CutStatus[iDecay3P], iCut);
                    }
                  }
                }
                rowProng3CutStatus(Prong3CutStatus[0], Prong3CutStatus[1], Prong3CutStatus[2], Prong3CutStatus[3]); // FIXME when we can do this by looping over kN3ProngDecays
              }

              // fill histograms
              if (fillHistograms) {
                registry.fill(HIST("hVtx3ProngX"), secondaryVertex3[0]);
                registry.fill(HIST("hVtx3ProngY"), secondaryVertex3[1]);
                registry.fill(HIST("hVtx3ProngZ"), secondaryVertex3[2]);
                std::array<std::array<float, 3>, 3> arr3Mom = {pvec0, pvec1, pvec2};
                for (int iDecay3P = 0; iDecay3P < kN3ProngDecays; iDecay3P++) {
                  if (TESTBIT(isSelected3ProngCand, iDecay3P)) {
                    if (whichHypo3Prong[iDecay3P] == 1 || whichHypo3Prong[iDecay3P] == 3) {
                      auto mass3Prong = RecoDecay::m(arr3Mom, arrMass3Prong[iDecay3P][0]);
                      switch (iDecay3P) {
                        case hf_cand_3prong::DecayType::DplusToPiKPi:
                          registry.fill(HIST("hMassDPlusToPiKPi"), mass3Prong);
                          break;
                        case hf_cand_3prong::DecayType::DsToKKPi:
                          registry.fill(HIST("hMassDsToKKPi"), mass3Prong);
                          break;
                        case hf_cand_3prong::DecayType::LcToPKPi:
                          registry.fill(HIST("hMassLcToPKPi"), mass3Prong);
                          break;
                        case hf_cand_3prong::DecayType::XicToPKPi:
                          registry.fill(HIST("hMassXicToPKPi"), mass3Prong);
                          break;
                      }
                    }
                    if (whichHypo3Prong[iDecay3P] >= 2) {
                      auto mass3Prong = RecoDecay::m(arr3Mom, arrMass3Prong[iDecay3P][1]);
                      switch (iDecay3P) {
                        case hf_cand_3prong::DecayType::DsToKKPi:
                          registry.fill(HIST("hMassDsToKKPi"), mass3Prong);
                          break;
                        case hf_cand_3prong::DecayType::LcToPKPi:
                          registry.fill(HIST("hMassLcToPKPi"), mass3Prong);
                          break;
                        case hf_cand_3prong::DecayType::XicToPKPi:
                          registry.fill(HIST("hMassXicToPKPi"), mass3Prong);
                          break;
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }

      int nTracks = 0;
      // auto nTracks = trackIndicesPerCollision.lastIndex() - trackIndicesPerCollision.firstIndex(); // number of tracks passing 2 and 3 prong selection in this collision
      nCand2 = rowTrackIndexProng2.lastIndex() - nCand2; // number of 2-prong candidates in this collision
      nCand3 = rowTrackIndexProng3.lastIndex() - nCand3; // number of 3-prong candidates in this collision

      if (fillHistograms) {
        registry.fill(HIST("hNTracks"), nTracks);
        registry.fill(HIST("hNCand2Prong"), nCand2);
        registry.fill(HIST("hNCand3Prong"), nCand3);
        registry.fill(HIST("hNCand2ProngVsNTracks"), nTracks, nCand2);
        registry.fill(HIST("hNCand3ProngVsNTracks"), nTracks, nCand3);
      }
    }
  } /// end of run2And3Prongs function

  void processNo2And3Prongs(SelectedCollisions const&)
  {
    // dummy
  }
  PROCESS_SWITCH(HfTrackIndexSkimCreator, processNo2And3Prongs, "Do not process 2-prongs and 3-prongs", false);

  void process2And3ProngsWithPvRefit( // soa::Join<aod::Collisions, aod::CentV0Ms>::iterator const& collision, //FIXME add centrality when option for variations to the process function appears
    SelectedCollisions const& collisions,
    aod::BCsWithTimestamps const& bcWithTimeStamps,
    FilteredTrackAssocSel const& trackIndices,
    TracksWithPVRefitAndDCA const& tracks)
  {
    run2And3Prongs<true>(collisions, bcWithTimeStamps, trackIndices, tracks);
  }

  PROCESS_SWITCH(HfTrackIndexSkimCreator, process2And3ProngsWithPvRefit, "Process 2-prong and 3-prong skim with PV refit", false);

  void process2And3ProngsNoPvRefit( // soa::Join<aod::Collisions, aod::CentV0Ms>::iterator const& collision, //FIXME add centrality when option for variations to the process function appears
    SelectedCollisions const& collisions,
    aod::BCsWithTimestamps const& bcWithTimeStamps,
    FilteredTrackAssocSel const& trackIndices,
    aod::TracksWCovDcaExtra const& tracks)
  {
    run2And3Prongs(collisions, bcWithTimeStamps, trackIndices, tracks);
  }

  PROCESS_SWITCH(HfTrackIndexSkimCreator, process2And3ProngsNoPvRefit, "Process 2-prong and 3-prong skim without PV refit", true);
};

//________________________________________________________________________________________________________________________

/// Pre-selection of cascade secondary vertices
/// It will produce in any case a Hf2Prongs object, but mixing a V0
/// with a track, instead of 2 tracks

/// to run: o2-analysis-weak-decay-indices --aod-file AO2D.root -b | o2-analysis-lambdakzerobuilder -b |
///         o2-analysis-trackextension -b | o2-analysis-hf-track-index-skim-creator -b

struct HfTrackIndexSkimCreatorCascades {
  Produces<aod::HfCascades> rowTrackIndexCasc;

  Configurable<bool> isRun2{"isRun2", false, "enable Run 2 or Run 3 GRP objects for magnetic field"};
  Configurable<bool> fillHistograms{"fillHistograms", true, "fill histograms"};
  // event selection
  // Configurable<int> triggerindex{"triggerindex", -1, "trigger index"};
  // vertexing
  // Configurable<double> bz{"bz", 5., "magnetic field"};
  Configurable<bool> propagateToPCA{"propagateToPCA", true, "create tracks version propagated to PCA"};
  Configurable<double> maxR{"maxR", 200., "reject PCA's above this radius"};
  Configurable<double> maxDZIni{"maxDZIni", 4., "reject (if>0) PCA candidate if tracks DZ exceeds threshold"};
  Configurable<double> minParamChange{"minParamChange", 1.e-3, "stop iterations if largest change of any X is smaller than this"};
  Configurable<double> minRelChi2Change{"minRelChi2Change", 0.9, "stop iterations if chi2/chi2old > this"};
  Configurable<bool> useAbsDCA{"useAbsDCA", true, "Minimise abs. distance rather than chi2"};
  Configurable<bool> useWeightedFinalPCA{"useWeightedFinalPCA", true, "Recalculate vertex position using track covariances, effective only if useAbsDCA is true"};
  // quality cut
  Configurable<bool> doCutQuality{"doCutQuality", true, "apply quality cuts"};
  // track cuts for bachelor
  Configurable<bool> tpcRefitBach{"tpcRefitBach", true, "request TPC refit bachelor"};
  Configurable<int> nCrossedRowsMinBach{"nCrossedRowsMinBach", 50, "min crossed rows bachelor"};
  // track cuts for V0 daughters
  Configurable<bool> tpcRefitV0Daugh{"tpcRefitV0Daugh", true, "request TPC refit V0 daughters"};
  Configurable<int> nCrossedRowsMinV0Daugh{"nCrossedRowsMinV0Daugh", 50, "min crossed rows V0 daughters"};
  Configurable<double> etaMinV0Daugh{"etaMinV0Daugh", -99999., "min. pseudorapidity V0 daughters"};
  Configurable<double> etaMaxV0Daugh{"etaMaxV0Daugh", 1.1, "max. pseudorapidity V0 daughters"};
  Configurable<double> ptMinV0Daugh{"ptMinV0Daugh", 0.05, "min. pT V0 daughters"};
  // bachelor cuts
  //  Configurable<float> dcabachtopv{"dcabachtopv", .1, "DCA Bach To PV"};
  //  Configurable<double> ptminbach{"ptminbach", -1., "min. track pT bachelor"};
  // v0 cuts
  Configurable<double> cpaV0Min{"cpaV0Min", .995, "min. cos PA V0"};                    // as in the task that create the V0s
  Configurable<double> dcaXYNegToPvMin{"dcaXYNegToPvMin", .1, "min. DCA_XY Neg To PV"}; // check: in HF Run 2, it was 0 at filtering
  Configurable<double> dcaXYPosToPvMin{"dcaXYPosToPvMin", .1, "min. DCA_XY Pos To PV"}; // check: in HF Run 2, it was 0 at filtering
  Configurable<double> cutInvMassV0{"cutInvMassV0", 0.05, "V0 candidate invariant mass difference wrt PDG"};
  // cascade cuts
  Configurable<double> ptCascCandMin{"ptCascCandMin", -1., "min. pT of the cascade candidate"};                    // PbPb 2018: use 1
  Configurable<double> cutInvMassCascLc{"cutInvMassCascLc", 1., "Lc candidate invariant mass difference wrt PDG"}; // for PbPb 2018: use 0.2
  // Configurable<double> cutCascDCADaughters{"cutCascDCADaughters", .1, "DCA between V0 and bachelor in cascade"};
  // CCDB
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> ccdbPathLut{"ccdbPathLut", "GLO/Param/MatLUT", "Path for LUT parametrization"};
  Configurable<std::string> ccdbPathGrp{"ccdbPathGrp", "GLO/GRP/GRP", "Path of the grp file (Run 2)"};
  Configurable<std::string> ccdbPathGrpMag{"ccdbPathGrpMag", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object (Run 3)"};
  // for debugging
#ifdef MY_DEBUG
  Configurable<std::vector<int>> indexK0Spos{"indexK0Spos", {729, 2866, 4754, 5457, 6891, 7824, 9243, 9810}, "indices of K0S positive daughters, for debug"};
  Configurable<std::vector<int>> indexK0Sneg{"indexK0Sneg", {730, 2867, 4755, 5458, 6892, 7825, 9244, 9811}, "indices of K0S negative daughters, for debug"};
  Configurable<std::vector<int>> indexProton{"indexProton", {717, 2810, 4393, 5442, 6769, 7793, 9002, 9789}, "indices of protons, for debug"};
#endif

  // Needed for PV refitting
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::base::MatLayerCylSet* lut;
  o2::base::Propagator::MatCorrType noMatCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;
  int runNumber;

  double massP{0.};
  double massK0s{0.};
  double massPi{0.};
  double massLc{0.};
  double mass2K0sP{0.}; // WHY HERE?

  using SelectedCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::HfSelCollision>>;
  using FilteredTrackAssocSel = soa::Filtered<soa::Join<aod::TrackAssoc, aod::HfSelTrack>>;

  Filter filterSelectCollisions = (aod::hf_sel_collision::whyRejectColl == 0);
  Filter filterSelectTrackIds = (aod::hf_sel_track::isSelProng >= 4);
  // Partition<MyTracks> TracksWithPVRefitAndDCA = aod::hf_sel_track::isSelProng >= 4;

  Preslice<FilteredTrackAssocSel> trackIndicesPerCollision = aod::track_association::collisionId;
  Preslice<aod::V0Datas> v0sPerCollision = aod::v0data::collisionId;

  // histograms
  HistogramRegistry registry{"registry"};

  void init(InitContext const& context)
  {
    if (!doprocessCascades) {
      return;
    }

    if (etaMinV0Daugh == -99999.) {
      etaMinV0Daugh.value = -etaMaxV0Daugh;
    }

    massP = o2::analysis::pdg::MassProton;
    massK0s = o2::analysis::pdg::MassK0Short;
    massPi = o2::analysis::pdg::MassPiPlus;
    massLc = o2::analysis::pdg::MassLambdaCPlus;

    ccdb->setURL(ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>(ccdbPathLut));
    runNumber = 0;

    if (fillHistograms) {
      registry.add("hVtx2ProngX", "2-prong candidates;#it{x}_{sec. vtx.} (cm);entries", {HistType::kTH1F, {{1000, -2., 2.}}});
      registry.add("hVtx2ProngY", "2-prong candidates;#it{y}_{sec. vtx.} (cm);entries", {HistType::kTH1F, {{1000, -2., 2.}}});
      registry.add("hVtx2ProngZ", "2-prong candidates;#it{z}_{sec. vtx.} (cm);entries", {HistType::kTH1F, {{1000, -2., 2.}}});
      registry.add("hMassLcToPK0S", "#Lambda_{c}^ candidates;inv. mass (p K_{S}^{0}) (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{500, 0., 5.}}});
    }
  }

  void processNoCascades(SelectedCollisions const&)
  {
    // dummy
  }
  PROCESS_SWITCH(HfTrackIndexSkimCreatorCascades, processNoCascades, "Do not skim HF -> V0 cascades", true);

  void processCascades(SelectedCollisions const& collisions,
                       aod::V0Datas const& v0s, // TODO: I am now assuming that the V0s are already filtered with my cuts (David's work to come)
                       FilteredTrackAssocSel const& trackIndices,
                       aod::TracksWCovDcaExtra const& tracks,
#ifdef MY_DEBUG
                       aod::McParticles& mcParticles,
#endif
                       aod::BCsWithTimestamps const&)
  {
    // set the magnetic field from CCDB
    for (const auto& collision : collisions) {
      auto bc = collision.bc_as<o2::aod::BCsWithTimestamps>();
      initCCDB(bc, runNumber, ccdb, isRun2 ? ccdbPathGrp : ccdbPathGrpMag, lut, isRun2);

      // Define o2 fitter, 2-prong
      o2::vertexing::DCAFitterN<2> fitter;
      fitter.setBz(o2::base::Propagator::Instance()->getNominalBz());
      fitter.setPropagateToPCA(propagateToPCA);
      fitter.setMaxR(maxR);
      fitter.setMinParamChange(minParamChange);
      fitter.setMinRelChi2Change(minRelChi2Change);
      // fitter.setMaxDZIni(1e9); // used in cascadeproducer.cxx, but not for the 2 prongs
      // fitter.setMaxChi2(1e9);  // used in cascadeproducer.cxx, but not for the 2 prongs
      fitter.setUseAbsDCA(useAbsDCA);
      fitter.setWeightedFinalPCA(useWeightedFinalPCA);

      // fist we loop over the bachelor candidate

      const auto thisCollId = collision.globalIndex();
      auto groupedBachTrackIndices = trackIndices.sliceBy(trackIndicesPerCollision, thisCollId);

      // for (const auto& bach : selectedTracks) {
      for (const auto& bachIdx : groupedBachTrackIndices) {
        auto bach = bachIdx.track_as<aod::TracksWCovDcaExtra>();

        MY_DEBUG_MSG(1, printf("\n"); LOG(info) << "Bachelor loop");
#ifdef MY_DEBUG
        auto indexBach = bach.mcParticleId();
        bool isProtonFromLc = isProtonFromLcFunc(indexBach, indexProton);
#endif
        // selections on the bachelor

        // // retrieve the selection flag that corresponds to this collision
        auto isSelProngBach = bachIdx.isSelProng();

        // pT cut
        if (!TESTBIT(isSelProngBach, CandidateType::CandV0bachelor)) {
          MY_DEBUG_MSG(isProtonFromLc, LOG(info) << "proton " << indexBach << ": rejected due to HFsel");
          continue;
        }
        if (tpcRefitBach) {
          if (!(bach.trackType() & o2::aod::track::TPCrefit)) {
            MY_DEBUG_MSG(isProtonFromLc, LOG(info) << "proton " << indexBach << ": rejected due to TPCrefit");
            continue;
          }
        }
        if (bach.tpcNClsCrossedRows() < nCrossedRowsMinBach) {
          MY_DEBUG_MSG(isProtonFromLc, LOG(info) << "proton " << indexBach << ": rejected due to minNUmberOfCrossedRows " << bach.tpcNClsCrossedRows() << " (cut " << nCrossedRowsMinBach << ")");
          continue;
        }
        MY_DEBUG_MSG(isProtonFromLc, LOG(info) << "KEPT! proton from Lc with daughters " << indexBach);

        auto trackBach = getTrackParCov(bach);

        auto groupedV0s = v0s.sliceBy(v0sPerCollision, thisCollId);
        // now we loop over the V0s
        for (const auto& v0 : groupedV0s) {
          MY_DEBUG_MSG(1, LOG(info) << "*** Checking next K0S");
          // selections on the V0 daughters
          const auto& trackV0DaughPos = v0.posTrack_as<aod::TracksWCovDcaExtra>();
          const auto& trackV0DaughNeg = v0.negTrack_as<aod::TracksWCovDcaExtra>();

          // check not to take the same track twice (as bachelor and V0 daughter)
          if (trackV0DaughPos.globalIndex() == bach.globalIndex() || trackV0DaughNeg.globalIndex() == bach.globalIndex()) {
            continue;
          }

#ifdef MY_DEBUG
          auto indexV0DaughPos = trackV0DaughPos.mcParticleId();
          auto indexV0DaughNeg = trackV0DaughNeg.mcParticleId();
          bool isK0SfromLc = isK0SfromLcFunc(indexV0DaughPos, indexV0DaughNeg, indexK0Spos, indexK0Sneg);

          bool isLc = isLcK0SpFunc(indexBach, indexV0DaughPos, indexV0DaughNeg, indexProton, indexK0Spos, indexK0Sneg);
#endif
          MY_DEBUG_MSG(isK0SfromLc, LOG(info) << "K0S from Lc found, trackV0DaughPos --> " << indexV0DaughPos << ", trackV0DaughNeg --> " << indexV0DaughNeg);

          MY_DEBUG_MSG(isK0SfromLc && isProtonFromLc,
                       LOG(info) << "ACCEPTED!!!";
                       LOG(info) << "proton belonging to a Lc found: label --> " << indexBach;
                       LOG(info) << "K0S belonging to a Lc found: trackV0DaughPos --> " << indexV0DaughPos << ", trackV0DaughNeg --> " << indexV0DaughNeg);

          MY_DEBUG_MSG(isLc, LOG(info) << "Combination of K0S and p which correspond to a Lc found!");

          if (tpcRefitV0Daugh) {
            if (!(trackV0DaughPos.trackType() & o2::aod::track::TPCrefit) ||
                !(trackV0DaughNeg.trackType() & o2::aod::track::TPCrefit)) {
              MY_DEBUG_MSG(isK0SfromLc, LOG(info) << "K0S with daughters " << indexV0DaughPos << " and " << indexV0DaughNeg << ": rejected due to TPCrefit");
              continue;
            }
          }
          if (trackV0DaughPos.tpcNClsCrossedRows() < nCrossedRowsMinV0Daugh ||
              trackV0DaughNeg.tpcNClsCrossedRows() < nCrossedRowsMinV0Daugh) {
            MY_DEBUG_MSG(isK0SfromLc, LOG(info) << "K0S with daughters " << indexV0DaughPos << " and " << indexV0DaughNeg << ": rejected due to minCrossedRows");
            continue;
          }
          //
          // if (trackV0DaughPos.dcaXY() < dcaXYPosToPvMin ||   // to the filters?
          //     trackV0DaughNeg.dcaXY() < dcaXYNegToPvMin) {
          //   continue;
          // }
          //
          if (trackV0DaughPos.pt() < ptMinV0Daugh || // to the filters? I can't for now, it is not in the tables
              trackV0DaughNeg.pt() < ptMinV0Daugh) {
            MY_DEBUG_MSG(isK0SfromLc, LOG(info) << "K0S with daughters " << indexV0DaughPos << " and " << indexV0DaughNeg << ": rejected due to minPt --> pos " << trackV0DaughPos.pt() << ", neg " << trackV0DaughNeg.pt() << " (cut " << ptMinV0Daugh << ")");
            continue;
          }
          if ((trackV0DaughPos.eta() > etaMaxV0Daugh || trackV0DaughPos.eta() < etaMinV0Daugh) || // to the filters? I can't for now, it is not in the tables
              (trackV0DaughNeg.eta() > etaMaxV0Daugh || trackV0DaughNeg.eta() < etaMinV0Daugh)) {
            MY_DEBUG_MSG(isK0SfromLc, LOG(info) << "K0S with daughters " << indexV0DaughPos << " and " << indexV0DaughNeg << ": rejected due to eta --> pos " << trackV0DaughPos.eta() << ", neg " << trackV0DaughNeg.eta() << " (cut " << etaMinV0Daugh << " to " << etaMaxV0Daugh << ")");
            continue;
          }

          // V0 invariant mass selection
          if (std::abs(v0.mK0Short() - massK0s) > cutInvMassV0) {
            MY_DEBUG_MSG(isK0SfromLc, LOG(info) << "K0S with daughters " << indexV0DaughPos << " and " << indexV0DaughNeg << ": rejected due to invMass --> " << v0.mK0Short() - massK0s << " (cut " << cutInvMassV0 << ")");
            continue; // should go to the filter, but since it is a dynamic column, I cannot use it there
          }

          // V0 cosPointingAngle selection
          if (v0.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) < cpaV0Min) {
            MY_DEBUG_MSG(isK0SfromLc, LOG(info) << "K0S with daughters " << indexV0DaughPos << " and " << indexV0DaughNeg << ": rejected due to cosPA --> " << v0.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) << " (cut " << cpaV0Min << ")");
            continue;
          }

          const std::array<float, 3> momentumV0 = {v0.px(), v0.py(), v0.pz()};

          // invariant-mass cut: we do it here, before updating the momenta of bach and V0 during the fitting to save CPU
          // TODO: but one should better check that the value here and after the fitter do not change significantly!!!
          mass2K0sP = RecoDecay::m(std::array{std::array{bach.px(), bach.py(), bach.pz()}, momentumV0}, std::array{massP, massK0s});
          if ((cutInvMassCascLc >= 0.) && (std::abs(mass2K0sP - massLc) > cutInvMassCascLc)) {
            MY_DEBUG_MSG(isK0SfromLc && isProtonFromLc, LOG(info) << "True Lc from proton " << indexBach << " and K0S pos " << indexV0DaughPos << " and neg " << indexV0DaughNeg << " rejected due to invMass cut: " << mass2K0sP << ", mass Lc " << massLc << " (cut " << cutInvMassCascLc << ")");
            continue;
          }

          MY_DEBUG_MSG(isK0SfromLc, LOG(info) << "KEPT! K0S from Lc with daughters " << indexV0DaughPos << " and " << indexV0DaughNeg);

          auto trackParCovV0DaughPos = getTrackParCov(trackV0DaughPos);
          trackParCovV0DaughPos.propagateTo(v0.posX(), o2::base::Propagator::Instance()->getNominalBz()); // propagate the track to the X closest to the V0 vertex
          auto trackParCovV0DaughNeg = getTrackParCov(trackV0DaughNeg);
          trackParCovV0DaughNeg.propagateTo(v0.negX(), o2::base::Propagator::Instance()->getNominalBz()); // propagate the track to the X closest to the V0 vertex
          std::array<float, 3> pVecV0 = {0., 0., 0.};
          std::array<float, 3> pVecBach = {0., 0., 0.};

          const std::array<float, 3> vertexV0 = {v0.x(), v0.y(), v0.z()};
          // we build the neutral track to then build the cascade
          auto trackV0 = o2::dataformats::V0(vertexV0, momentumV0, {0, 0, 0, 0, 0, 0}, trackParCovV0DaughPos, trackParCovV0DaughNeg); // build the V0 track

          // now we find the DCA between the V0 and the bachelor, for the cascade
          int nCand2 = fitter.process(trackV0, trackBach);
          MY_DEBUG_MSG(isK0SfromLc && isProtonFromLc, LOG(info) << "Fitter result = " << nCand2 << " proton = " << indexBach << " and K0S pos " << indexV0DaughPos << " and neg " << indexV0DaughNeg);
          MY_DEBUG_MSG(isLc, LOG(info) << "Fitter result for true Lc = " << nCand2);
          if (nCand2 == 0) {
            continue;
          }
          fitter.propagateTracksToVertex();        // propagate the bach and V0 to the Lc vertex
          fitter.getTrack(0).getPxPyPzGlo(pVecV0); // take the momentum at the Lc vertex
          fitter.getTrack(1).getPxPyPzGlo(pVecBach);

          // cascade candidate pT cut
          auto ptCascCand = RecoDecay::pt(pVecBach, pVecV0);
          if (ptCascCand < ptCascCandMin) {
            MY_DEBUG_MSG(isK0SfromLc && isProtonFromLc, LOG(info) << "True Lc from proton " << indexBach << " and K0S pos " << indexV0DaughPos << " and neg " << indexV0DaughNeg << " rejected due to pt cut: " << ptCascCand << " (cut " << ptCascCandMin << ")");
            continue;
          }

          // invariant mass
          // re-calculate invariant masses with updated momenta, to fill the histogram
          mass2K0sP = RecoDecay::m(std::array{pVecBach, pVecV0}, std::array{massP, massK0s});

          std::array<float, 3> posCasc = {0., 0., 0.};
          const auto& cascVtx = fitter.getPCACandidate();
          for (int i = 0; i < 3; i++) {
            posCasc[i] = cascVtx[i];
          }

          // fill table row
          rowTrackIndexCasc(thisCollId, bach.globalIndex(), v0.v0Id());
          // fill histograms
          if (fillHistograms) {
            MY_DEBUG_MSG(isK0SfromLc && isProtonFromLc && isLc, LOG(info) << "KEPT! True Lc from proton " << indexBach << " and K0S pos " << indexV0DaughPos << " and neg " << indexV0DaughNeg);
            registry.fill(HIST("hVtx2ProngX"), posCasc[0]);
            registry.fill(HIST("hVtx2ProngY"), posCasc[1]);
            registry.fill(HIST("hVtx2ProngZ"), posCasc[2]);
            registry.fill(HIST("hMassLcToPK0S"), mass2K0sP);
          }
        } // loop over V0s
      }   // loop over tracks
    }     // loop over collisions
  }       // processCascades

  PROCESS_SWITCH(HfTrackIndexSkimCreatorCascades, processCascades, "Skim HF -> V0 cascades", false);
};

struct HfTrackIndexSkimCreatorLfCascades {
  Produces<aod::HfCascLf2Prongs> rowTrackIndexCasc2Prong;
  Produces<aod::HfCascLf3Prongs> rowTrackIndexCasc3Prong;

  // whether to do or not validation plots
  Configurable<bool> fillHistograms{"fillHistograms", true, "fill histograms"};
  Configurable<bool> do3Prong{"do3Prong", false, "do 3-prong cascade"};

  // vertexing parameters
  Configurable<bool> propagateToPCA{"propagateToPCA", true, "create tracks version propagated to PCA"};
  Configurable<double> maxR{"maxR", 200., "reject PCA's above this radius"};
  Configurable<double> maxDZIni{"maxDZIni", 4., "reject (if>0) PCA candidate if tracks DZ exceeds threshold"};
  Configurable<double> minParamChange{"minParamChange", 1.e-3, "stop iterations if largest change of any X is smaller than this"};
  Configurable<double> minRelChi2Change{"minRelChi2Change", 0.9, "stop iterations if chi2/chi2old > this"};
  Configurable<bool> useAbsDCA{"useAbsDCA", true, "Minimise abs. distance rather than chi2"};
  Configurable<bool> useWeightedFinalPCA{"useWeightedFinalPCA", true, "Recalculate vertex position using track covariances, effective only if useAbsDCA is true"};
  Configurable<bool> rejDiffCollTrack{"rejDiffCollTrack", true, "Reject tracks coming from different collisions"};

  // quality cut
  Configurable<bool> doCutQuality{"doCutQuality", true, "apply quality cuts"};

  // Selection criteria
  Configurable<double> v0CosPA{"v0CosPA", 0.95, "V0 CosPA"};       // double -> N.B. dcos(x)/dx = 0 at x=0)
  Configurable<double> cascCosPA{"cascCosPA", 0.95, "Casc CosPA"}; // double -> N.B. dcos(x)/dx = 0 at x=0)
  Configurable<float> dcaV0Dau{"dcaV0Dau", 2.0, "DCA V0 Daughters"};
  Configurable<float> dcaCascDau{"dcaCascDau", 1.0, "DCA Casc Daughters"};
  Configurable<float> dcaNegToPv{"dcaNegToPv", .05, "DCA Neg To PV"};
  Configurable<float> dcaPosToPv{"dcaPosToPv", .05, "DCA Pos To PV"};
  Configurable<float> dcaBachToPv{"dcaBachToPv", .05, "DCA Bach To PV"};
  Configurable<float> dcaV0ToPv{"dcaV0ToPv", .05, "DCA V0 To PV"};
  Configurable<float> v0Radius{"v0Radius", 0.9, "V0 radius"};
  Configurable<float> cascRadius{"cascRadius", 0.5, "Casc radius"};
  Configurable<float> v0MassWindow{"v0MassWindow", 0.008, "V0 mass window"};

  // Track identification configurables
  Configurable<float> tpcNsigmaBachelor{"tpcNsigmaBachelor", 4, "TPC NSigma bachelor (>10 is no cut)"};
  Configurable<float> tpcNsigmaProton{"tpcNsigmaProton", 4, "TPC NSigma proton <- lambda (>10 is no cut)"};
  Configurable<float> tpcNsigmaPion{"tpcNsigmaPion", 4, "TPC NSigma pion <- lambda (>10 is no cut)"};

  // magnetic field setting from CCDB
  Configurable<bool> isRun2{"isRun2", false, "enable Run 2 or Run 3 GRP objects for magnetic field"};
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "URL of the CCDB repository"};
  Configurable<std::string> ccdbPathLut{"ccdbPathLut", "GLO/Param/MatLUT", "Path for LUT parametrization"};
  Configurable<std::string> ccdbPathGrp{"ccdbPathGrp", "GLO/GRP/GRP", "Path of the grp file (Run 2)"};
  Configurable<std::string> ccdbPathGrpMag{"ccdbPathGrpMag", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object (Run 3)"};

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::base::MatLayerCylSet* lut;
  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;
  int runNumber;

  static constexpr int kN2ProngDecays = hf_cand_casc_lf_2prong::DecayType::N2ProngDecays; // number of 2-prong hadron types
  static constexpr int kN3ProngDecays = hf_cand_casc_lf_3prong::DecayType::N3ProngDecays; // number of 3-prong hadron types
  std::array<std::array<std::array<double, 2>, 2>, kN2ProngDecays> arrMass2Prong;
  std::array<std::array<std::array<double, 3>, 2>, kN3ProngDecays> arrMass3Prong;

  double massP{0.};
  double massPi{0.};
  double massXi{0.};
  double massOmega{0.};
  double massXiczero{0.};
  double massXicplus{0.};

  // histograms
  HistogramRegistry registry{"registry"};

  void init(InitContext const&)
  {
    if (!doprocessLfCascades) {
      return;
    }

    massP = o2::analysis::pdg::MassProton;
    massPi = o2::analysis::pdg::MassPiPlus;
    massXi = o2::analysis::pdg::MassXiMinus;
    massOmega = o2::analysis::pdg::MassOmegaMinus;
    massXiczero = o2::analysis::pdg::MassXiCZero;
    massXicplus = o2::analysis::pdg::MassXiCPlus;

    arrMass2Prong[hf_cand_casc_lf_2prong::DecayType::XiczeroToXiPi] = std::array{std::array{massXi, massPi},
                                                                                 std::array{massPi, massXi}};

    arrMass2Prong[hf_cand_casc_lf_2prong::DecayType::OmegaczeroToOmegaPi] = std::array{std::array{massOmega, massPi},
                                                                                       std::array{massPi, massOmega}};

    arrMass3Prong[hf_cand_casc_lf_3prong::DecayType::XicplusToXiPiPi] = std::array{std::array{massXi, massPi, massPi},
                                                                                   std::array{massPi, massPi, massXi}};

    ccdb->setURL(ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>(ccdbPathLut));
    runNumber = 0;

    if (fillHistograms) {
      AxisSpec ptAxis = {200, 0.0f, 10.0f, "it{p}_{T} (GeV/c)"};
      AxisSpec massAxisXi = {200, 1.222f, 1.422f, "Inv. Mass (GeV/c^{2})"};
      AxisSpec massAxisOmega = {200, 1.572f, 1.772f, "Inv. Mass (GeV/c^{2})"};

      registry.add("hCandidateCounter", "hCandidateCounter", {HistType::kTH1F, {{10, 0.0f, 10.0f}}});

      // Cascade mass spectra
      registry.add("hMassXiMinus", "hMassXiMinus", {HistType::kTH1F, {{3000, 0.0f, 3.0f, "Inv. Mass (GeV/c^{2})"}}});
      registry.add("hMassXiPlus", "hMassXiPlus", {HistType::kTH1F, {{3000, 0.0f, 3.0f, "Inv. Mass (GeV/c^{2}²)"}}});
      registry.add("hMassOmegaMinus", "hMassOmegaMinus", {HistType::kTH1F, {{3000, 0.0f, 3.0f, "Inv. Mass (GeV/c^{2})"}}});
      registry.add("hMassOmegaPlus", "hMassOmegaPlus", {HistType::kTH1F, {{3000, 0.0f, 3.0f, "Inv. Mass (GeV/c^{2})"}}});
      registry.add("h2dMassXiMinus", "h2dMassXiMinus", {HistType::kTH2F, {ptAxis, massAxisXi}});
      registry.add("h2dMassXiPlus", "h2dMassXiPlus", {HistType::kTH2F, {ptAxis, massAxisXi}});
      registry.add("h2dMassOmegaMinus", "h2dMassOmegaMinus", {HistType::kTH2F, {ptAxis, massAxisOmega}});
      registry.add("h2dMassOmegaPlus", "h2dMassOmegaPlus", {HistType::kTH2F, {ptAxis, massAxisOmega}});

      // Cascade topology
      registry.add("hV0Radius", "hV0Radius", {HistType::kTH1F, {{1000, 0.0f, 100.0f, "cm"}}});
      registry.add("hCascRadius", "hCascRadius", {HistType::kTH1F, {{1000, 0.0f, 100.0f, "cm"}}});
      registry.add("hV0CosPA", "hV0CosPA", {HistType::kTH1F, {{1000, 0.95f, 1.0f}}});
      registry.add("hCascCosPA", "hCascCosPA", {HistType::kTH1F, {{1000, 0.95f, 1.0f}}});
      registry.add("hDCAPosToPV", "hDCAPosToPV", {HistType::kTH1F, {{1000, -10.0f, 10.0f, "cm"}}});
      registry.add("hDCANegToPV", "hDCANegToPV", {HistType::kTH1F, {{1000, -10.0f, 10.0f, "cm"}}});
      registry.add("hDCABachToPV", "hDCABachToPV", {HistType::kTH1F, {{1000, -10.0f, 10.0f, "cm"}}});
      registry.add("hDCAV0ToPV", "hDCAV0ToPV", {HistType::kTH1F, {{1000, -10.0f, 10.0f, "cm"}}});
      registry.add("hDCAV0Dau", "hDCAV0Dau", {HistType::kTH1F, {{1000, 0.0f, 10.0f, "cm^{2}"}}});
      registry.add("hDCACascDau", "hDCACascDau", {HistType::kTH1F, {{1000, 0.0f, 10.0f, "cm^{2}"}}});
      registry.add("hLambdaMass", "hLambdaMass", {HistType::kTH1F, {{1000, 0.0f, 10.0f, "Inv. Mass (GeV/c^{2})"}}});

      registry.add("hVtx2ProngX", "2-prong candidates;#it{x}_{sec. vtx.} (cm);entries", {HistType::kTH1F, {{1000, -2., 2.}}});
      registry.add("hVtx2ProngY", "2-prong candidates;#it{y}_{sec. vtx.} (cm);entries", {HistType::kTH1F, {{1000, -2., 2.}}});
      registry.add("hVtx2ProngZ", "2-prong candidates;#it{z}_{sec. vtx.} (cm);entries", {HistType::kTH1F, {{1000, -20., 20.}}});
      registry.add("hMassXicToXiPi", "2-prong candidates;inv. mass (#Xi pi) (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{100, 2., 3.}}});
      registry.add("hVtx3ProngX", "3-prong candidates;#it{x}_{sec. vtx.} (cm);entries", {HistType::kTH1F, {{1000, -2., 2.}}});
      registry.add("hVtx3ProngY", "3-prong candidates;#it{y}_{sec. vtx.} (cm);entries", {HistType::kTH1F, {{1000, -2., 2.}}});
      registry.add("hVtx3ProngZ", "3-prong candidates;#it{z}_{sec. vtx.} (cm);entries", {HistType::kTH1F, {{1000, -20., 20.}}});

      // mass spectra
      registry.add("hMassXicZeroToXiPi", "2-prong candidates;inv. mass (#Xi #pi) (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{100, 2., 3.}}});
      registry.add("hMassOmegacZeroToOmegaPi", "2-prong candidates;inv. mass (#Omega #pi) (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{100, 2., 3.}}});
      registry.add("hMassXicPlusToXiPiPi", "3-prong candidates;inv. mass (#Xi #pi #pi) (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{100, 2., 3.}}});
    }
  }

  Filter filterSelectCollisions = (aod::hf_sel_collision::whyRejectColl == 0);
  Filter filterSelectTrackIds = (aod::hf_sel_track::isSelProng > 0);

  using SelectedCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::HfSelCollision>>;
  using TracksWithPVRefitAndDCA = soa::Join<aod::TracksWDca, aod::HfPvRefitTrack>;
  using V0Full = soa::Join<aod::V0Datas, aod::V0Covs>;

  Preslice<TracksWithPVRefitAndDCA> tracksPerCollision = aod::track::collisionId; // needed for PV refit
  using FilteredTrackAssocSel = soa::Filtered<soa::Join<aod::TrackAssoc, aod::HfSelTrack>>;
  Preslice<FilteredTrackAssocSel> trackIndicesPerCollision = aod::track_association::collisionId;
  Preslice<aod::CascDataFull> cascadesPerCollision = aod::cascdata::collisionId;

  /// Single-cascade cuts for 2-prongs or 3-prongs
  /// From cascadeanalysis.cxx w/o PID and Centrality study
  ///  Function to process cascades and generate corresponding invariant-mass distributions
  template <typename TCascade, typename T1, typename T2, typename T3>
  bool isPreselectedCascade(const TCascade& casc, const T1& bachTrackCast, const T2& posTrackCast, const T3& negTrackCast, const float& pvx, const float& pvy, const float& pvz)
  {
    // Logic: either you have enough TPC clusters, OR you enabled ITSSA and have enough ITS clusters as requested
    // N.B.: This will require dedicated studies!
    registry.fill(HIST("hCandidateCounter"), 2.5); // okay track quality

    if (casc.v0cosPA(pvx, pvy, pvz) > v0CosPA &&
        casc.casccosPA(pvx, pvy, pvz) > cascCosPA &&
        casc.dcacascdaughters() > dcaCascDau &&
        casc.dcaV0daughters() > dcaV0Dau &&
        casc.dcanegtopv() > dcaNegToPv &&
        casc.dcapostopv() > dcaPosToPv &&
        casc.dcabachtopv() > dcaBachToPv &&
        casc.dcav0topv(pvx, pvy, pvz) > dcaV0ToPv &&
        casc.v0radius() > v0Radius &&
        casc.cascradius() > cascRadius &&
        std::abs(casc.mLambda() - 1.115683) < v0MassWindow) {
      registry.fill(HIST("hCandidateCounter"), 3.5); // pass cascade selections
      if (casc.sign() < 0) {                         // FIXME: could be done better...
        registry.fill(HIST("hMassXiMinus"), casc.mXi());
        registry.fill(HIST("hMassOmegaMinus"), casc.mOmega());
      } else {
        registry.fill(HIST("hMassXiPlus"), casc.mXi());
        registry.fill(HIST("hMassOmegaPlus"), casc.mOmega());
      }

      if (fillHistograms) {
        // The basic eleven!
        registry.fill(HIST("hV0Radius"), casc.v0radius());
        registry.fill(HIST("hCascRadius"), casc.cascradius());
        registry.fill(HIST("hV0CosPA"), casc.v0cosPA(pvx, pvy, pvz));
        registry.fill(HIST("hCascCosPA"), casc.casccosPA(pvx, pvy, pvz));
        registry.fill(HIST("hDCAPosToPV"), casc.dcapostopv());
        registry.fill(HIST("hDCANegToPV"), casc.dcanegtopv());
        registry.fill(HIST("hDCABachToPV"), casc.dcabachtopv());
        registry.fill(HIST("hDCAV0ToPV"), casc.dcav0topv(pvx, pvy, pvz));
        registry.fill(HIST("hDCAV0Dau"), casc.dcaV0daughters());
        registry.fill(HIST("hDCACascDau"), casc.dcacascdaughters());
        registry.fill(HIST("hLambdaMass"), casc.mLambda());
      }
      return true;
    }
    return false;
  }

  void processNoLfCascades(SelectedCollisions const&)
  {
    // dummy
  }
  PROCESS_SWITCH(HfTrackIndexSkimCreatorLfCascades, processNoLfCascades, "Do not skim LF cascades", true);

  void processLfCascades(SelectedCollisions const& collisions,
                         aod::CascDataFull const& cascades,
                         FilteredTrackAssocSel const& trackIndices,
                         aod::TracksWCovDca const& tracks,
                         aod::BCsWithTimestamps const&,
                         aod::V0sLinked const&,
                         V0Full const&)
  {
    for (const auto& collision : collisions) {

      // set the magnetic field from CCDB
      auto bc = collision.bc_as<o2::aod::BCsWithTimestamps>();
      initCCDB(bc, runNumber, ccdb, isRun2 ? ccdbPathGrp : ccdbPathGrpMag, lut, isRun2);
      auto magneticField = o2::base::Propagator::Instance()->getNominalBz(); // z component

      // Define o2 fitter, 2-prong
      o2::vertexing::DCAFitterN<2> df2;
      df2.setBz(magneticField);
      df2.setPropagateToPCA(propagateToPCA);
      df2.setMaxR(maxR);
      df2.setMaxDZIni(maxDZIni);
      df2.setMinParamChange(minParamChange);
      df2.setMinRelChi2Change(minRelChi2Change);
      df2.setUseAbsDCA(useAbsDCA);
      df2.setWeightedFinalPCA(useWeightedFinalPCA);

      // 3-prong vertex fitter
      o2::vertexing::DCAFitterN<3> df3;
      df3.setBz(magneticField);
      df3.setPropagateToPCA(propagateToPCA);
      df3.setMaxR(maxR);
      df3.setMaxDZIni(maxDZIni);
      df3.setMinParamChange(minParamChange);
      df3.setMinRelChi2Change(minRelChi2Change);
      df3.setUseAbsDCA(useAbsDCA);
      df3.setWeightedFinalPCA(useWeightedFinalPCA);

      // cascade loop
      auto thisCollId = collision.globalIndex();
      auto groupedCascades = cascades.sliceBy(cascadesPerCollision, thisCollId);
      for (const auto& casc : groupedCascades) {

        registry.fill(HIST("hCandidateCounter"), 0.5); // all candidates

        auto trackXiDauCharged = casc.bachelor_as<aod::TracksWCovDca>(); // pion <- xi track from TracksWithPVRefitAndDCA table
        // cascade daughter - V0
        if (!casc.v0_as<aod::V0sLinked>().has_v0Data()) { // check that V0 data are stored
          continue;
        }
        registry.fill(HIST("hCandidateCounter"), 1.5); // v0data exists
        auto v0index = casc.v0_as<aod::V0sLinked>();
        auto v0 = v0index.v0Data_as<V0Full>(); // V0 element from LF table containing V0 info
        // V0 positive daughter
        auto trackV0DauPos = v0.posTrack_as<aod::TracksWCovDca>(); // p <- V0 track (positive track) from TracksWithPVRefitAndDCA table
        // V0 negative daughter
        auto trackV0DauNeg = v0.negTrack_as<aod::TracksWCovDca>(); // pion <- V0 track (negative track) from TracksWithPVRefitAndDCA table

        // check that particles come from the same collision
        if (rejDiffCollTrack) {
          if (trackV0DauPos.collisionId() != trackV0DauNeg.collisionId()) {
            continue;
          }
          if (trackXiDauCharged.collisionId() != trackV0DauPos.collisionId()) {
            continue;
          }
        }

        if (!(isPreselectedCascade(casc, trackXiDauCharged, trackV0DauPos, trackV0DauNeg, collision.posX(), collision.posY(), collision.posZ()))) {
          continue;
        }

        o2::vertexing::DCAFitterN<2> dfc;
        dfc.setBz(magneticField);
        dfc.setPropagateToPCA(propagateToPCA);
        dfc.setMaxR(maxR);
        dfc.setMaxDZIni(maxDZIni);
        dfc.setMinParamChange(minParamChange);
        dfc.setMinRelChi2Change(minRelChi2Change);
        dfc.setUseAbsDCA(useAbsDCA);
        dfc.setWeightedFinalPCA(useWeightedFinalPCA);

        auto trackParVarXiDauCharged = getTrackParCov(trackXiDauCharged);

        // Set up covariance matrices (should in fact be optional)
        std::array<float, 21> covV = {0.};
        constexpr int MomInd[6] = {9, 13, 14, 18, 19, 20}; // cov matrix elements for momentum component
        for (int i = 0; i < 6; i++) {
          covV[MomInd[i]] = v0.momentumCovMat()[i];
          covV[i] = v0.positionCovMat()[i];
        }
        auto trackV0 = o2::track::TrackParCov(
          {v0.x(), v0.y(), v0.z()},
          {v0.pxpos() + v0.pxneg(), v0.pypos() + v0.pyneg(), v0.pzpos() + v0.pzneg()},
          covV, 0, true);
        trackV0.setAbsCharge(0);
        trackV0.setPID(o2::track::PID::Lambda);

        // reconstruct the cascade
        if (dfc.process(trackV0, trackParVarXiDauCharged) == 0) {
          continue;
        }

        std::array<float, 3> pvecV0;        // V0
        std::array<float, 3> pvecXiDauPion; // bach pion

        dfc.getTrack(0).getPxPyPzGlo(pvecV0);
        dfc.getTrack(1).getPxPyPzGlo(pvecXiDauPion);

        std::array<float, 3> coordVtxCasc = dfc.getPCACandidatePos();
        std::array<float, 6> covVtxCasc = dfc.calcPCACovMatrixFlat();
        std::array<float, 3> pvecCascAsM = {pvecV0[0] + pvecXiDauPion[0], pvecV0[1] + pvecXiDauPion[1], pvecV0[2] + pvecXiDauPion[2]};

        auto trackCasc = o2::dataformats::V0(coordVtxCasc, pvecCascAsM, covVtxCasc, trackV0, trackParVarXiDauCharged);

        //-------------------combining cascade and pion tracks--------------------------
        // first loop over positive tracks
        auto groupedBachTrackIndices = trackIndices.sliceBy(trackIndicesPerCollision, thisCollId);
        for (auto trackIdPion1 = groupedBachTrackIndices.begin(); trackIdPion1 != groupedBachTrackIndices.end(); ++trackIdPion1) {
          auto trackPion1 = trackIdPion1.track_as<aod::TracksWCovDca>();

          if ((rejDiffCollTrack) && (trackXiDauCharged.collisionId() != trackPion1.collisionId())) {
            continue;
          }

          // ask for opposite sign daughters
          if (trackPion1.sign() * trackXiDauCharged.sign() >= 0) {
            continue;
          }

          // check not to take the same particle twice in the decay chain
          if (trackPion1.globalIndex() == trackXiDauCharged.globalIndex() || trackPion1.globalIndex() == trackV0DauPos.globalIndex() || trackPion1.globalIndex() == trackV0DauNeg.globalIndex() || trackPion1.globalIndex() == casc.globalIndex()) {
            continue;
          }

          // primary pion track to be processed with DCAFitter
          auto trackParVarPion1 = getTrackParCov(trackPion1);

          // first loop over tracks
          if (do3Prong) {
            // second loop over positive tracks
            for (auto trackIdPion2 = trackIdPion1 + 1; trackIdPion2 != groupedBachTrackIndices.end(); ++trackIdPion2) {
              auto trackPion2 = trackIdPion2.track_as<aod::TracksWCovDca>();

              if ((rejDiffCollTrack) && (trackXiDauCharged.collisionId() != trackPion2.collisionId())) {
                continue;
              }

              // ask for same sign daughters
              if (trackPion2.sign() * trackPion1.sign() <= 0) {
                continue;
              }

              // check not to take the same particle twice in the decay chain
              if (trackPion2.globalIndex() == trackXiDauCharged.globalIndex() || trackPion2.globalIndex() == trackV0DauPos.globalIndex() || trackPion2.globalIndex() == trackV0DauNeg.globalIndex() || trackPion2.globalIndex() == casc.globalIndex()) {
                continue;
              }

              // primary pion track to be processed with DCAFitter
              auto trackParVarPion2 = getTrackParCov(trackPion2);

              // reconstruct Xic with DCAFitter
              if (df3.process(trackCasc, trackParVarPion1, trackParVarPion2) == 0) {
                continue;
              }

              std::array<float, 3> pVec1 = {0.};
              std::array<float, 3> pVec2 = {0.};
              std::array<float, 3> pVec3 = {0.};

              df3.getTrack(0).getPxPyPzGlo(pVec1); // take the momentum at the Xic vertex
              df3.getTrack(1).getPxPyPzGlo(pVec2);
              df3.getTrack(2).getPxPyPzGlo(pVec3);

              // std::array<float, 3> secondaryVertex3 = {0., 0., 0.};
              const auto& secondaryVertex3 = df3.getPCACandidate();

              // fill table row
              rowTrackIndexCasc3Prong(casc.globalIndex(),
                                      trackPion1.globalIndex(),
                                      trackPion2.globalIndex());

              // fill histograms
              if (fillHistograms) {
                registry.fill(HIST("hVtx3ProngX"), secondaryVertex3[0]);
                registry.fill(HIST("hVtx3ProngY"), secondaryVertex3[1]);
                registry.fill(HIST("hVtx3ProngZ"), secondaryVertex3[2]);

                std::array<std::array<float, 3>, 3> arr3Mom = {pVec1, pVec2, pVec3};
                for (int iDecay3P = 0; iDecay3P < kN3ProngDecays; iDecay3P++) {
                  auto mass3Prong = RecoDecay::m(arr3Mom, arrMass3Prong[iDecay3P][0]);
                  switch (iDecay3P) {
                    case hf_cand_casc_lf_3prong::DecayType::XicplusToXiPiPi:
                      registry.fill(HIST("hMassXicPlusToXiPiPi"), mass3Prong);
                      break;
                  }
                }
              }
            }
          }

          if (df2.process(trackCasc, trackParVarPion1) == 0) {
            continue;
          }

          std::array<float, 3> pVec1 = {0.};
          std::array<float, 3> pVec2 = {0.};

          df2.getTrack(0).getPxPyPzGlo(pVec1);
          df2.getTrack(1).getPxPyPzGlo(pVec2);

          const auto& secondaryVertex2 = df2.getPCACandidate();

          // fill table row
          rowTrackIndexCasc2Prong(casc.globalIndex(),
                                  trackPion1.globalIndex());

          // fill histograms
          if (fillHistograms) {
            registry.fill(HIST("hVtx2ProngX"), secondaryVertex2[0]);
            registry.fill(HIST("hVtx2ProngY"), secondaryVertex2[1]);
            registry.fill(HIST("hVtx2ProngZ"), secondaryVertex2[2]);

            std::array<std::array<float, 3>, 2> arrMom = {pVec1, pVec2};
            for (int iDecay2P = 0; iDecay2P < kN2ProngDecays; iDecay2P++) {
              auto mass2Prong = RecoDecay::m(arrMom, arrMass2Prong[iDecay2P][0]);
              switch (iDecay2P) {
                case hf_cand_casc_lf_2prong::DecayType::XiczeroToXiPi:
                  registry.fill(HIST("hMassXicZeroToXiPi"), mass2Prong);
                  break;
                case hf_cand_casc_lf_2prong::DecayType::OmegaczeroToOmegaPi:
                  registry.fill(HIST("hMassOmegacZeroToXiPi"), mass2Prong);
                  break;
              }
            }
          }
        } // loop over pion
      }   // loop over cascade
    }     // loop over collisions
  }       // processLfCascades

  PROCESS_SWITCH(HfTrackIndexSkimCreatorLfCascades, processLfCascades, "Skim LF cascades", false);
};

//________________________________________________________________________________________________________________________
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{};
  workflow.push_back(adaptAnalysisTask<HfTrackIndexSkimCreatorTagSelCollisions>(cfgc));
  workflow.push_back(adaptAnalysisTask<HfTrackIndexSkimCreatorTagSelTracks>(cfgc));
  workflow.push_back(adaptAnalysisTask<HfTrackIndexSkimCreator>(cfgc));
  workflow.push_back(adaptAnalysisTask<HfTrackIndexSkimCreatorCascades>(cfgc));
  workflow.push_back(adaptAnalysisTask<HfTrackIndexSkimCreatorLfCascades>(cfgc));
  return workflow;
}
