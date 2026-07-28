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

/// \file kstarpbpb.cxx
/// \brief Code for K*(892)^0 resonance flow and spin alignment analysis
/// \author sourav.kundu@cern.ch , sarjeeta.gami@cern.ch
///

#include "PWGLF/DataModel/EPCalibrationTables.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/CCDB/RCTSelectionFlags.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CCDB/CcdbApi.h>
#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoAHelpers.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/BinningPolicy.h>
#include <Framework/Configurable.h>
#include <Framework/GroupedCombinations.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <Math/GenVector/Boost.h>
#include <Math/Vector3Dfwd.h>
#include <Math/Vector4D.h> // IWYU pragma: keep (do not replace with Math/Vector4Dfwd.h)
#include <Math/Vector4Dfwd.h>
#include <TF1.h>
#include <THn.h>
#include <TPDGCode.h>

#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::rctsel;

struct Kstarpbpb {

  struct : ConfigurableGroup {
    Configurable<std::string> cfgURL{"cfgURL", "http://alice-ccdb.cern.ch", "Address of the CCDB to browse"};
    Configurable<int64_t> nolaterthan{"nolaterthan", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "Latest acceptable timestamp of creation for the object"};
  } cfgCcdbParam;

  // Enable access to the CCDB for the offset and correction constants and save them in dedicated variables.
  Service<o2::ccdb::BasicCCDBManager> ccdb{};
  o2::ccdb::CcdbApi ccdbApi;
  // Service<o2::framework::O2DatabasePDG> pdg;
  struct RCTCut : ConfigurableGroup {
    Configurable<bool> requireRCTFlagChecker{"requireRCTFlagChecker", true, "Check event quality in run condition table"};
    Configurable<std::string> cfgEvtRCTFlagCheckerLabel{"cfgEvtRCTFlagCheckerLabel", "CBT_hadronPID", "Evt sel: RCT flag checker label"};
    Configurable<bool> cfgEvtRCTFlagCheckerZDCCheck{"cfgEvtRCTFlagCheckerZDCCheck", false, "Evt sel: RCT flag checker ZDC check"};
    Configurable<bool> cfgEvtRCTFlagCheckerLimitAcceptAsBad{"cfgEvtRCTFlagCheckerLimitAcceptAsBad", true, "Evt sel: RCT flag checker treat Limited Acceptance As Bad"};

    RCTFlagsChecker rctChecker;
  };

  RCTCut rctCut;

  // CCDB options
  // Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  // Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  // Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  // Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
  // Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};

  // events
  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};
  Configurable<float> cfgCutCentrality{"cfgCutCentrality", 80.0f, "Accepted maximum Centrality"};
  // track
  Configurable<float> cfgCutCharge{"cfgCutCharge", 0.0, "cut on Charge"};
  Configurable<float> cfgCutPT{"cfgCutPT", 0.2, "PT cut on daughter track"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8, "Eta cut on daughter track"};
  Configurable<float> cfgCutDCAxy{"cfgCutDCAxy", 2.0f, "DCAxy range for tracks"};
  Configurable<float> cfgCutDCAz{"cfgCutDCAz", 2.0f, "DCAz range for tracks"};
  Configurable<bool> useGlobalTrack{"useGlobalTrack", true, "use Global track"};
  Configurable<bool> usepolar{"usepolar", true, "flag to fill type of SA"};
  Configurable<float> nsigmaCutTOF{"nsigmaCutTOF", 3.0, "Value of the TOF Nsigma cut"};
  Configurable<float> nsigmaCutTPC{"nsigmaCutTPC", 3.0, "Value of the TPC Nsigma cut"};
  Configurable<bool> isTOFOnly{"isTOFOnly", false, "use TOF only PID"};
  Configurable<float> nsigmaCutCombined{"nsigmaCutCombined", 3.0, "Value of the TOF Nsigma cut"};
  Configurable<int> cfgNoMixedEvents{"cfgNoMixedEvents", 1, "Number of mixed events per event"};
  Configurable<int> cfgITScluster{"cfgITScluster", 0, "Number of ITS cluster"};
  Configurable<int> cfgTPCcluster{"cfgTPCcluster", 70, "Number of TPC cluster"};
  Configurable<double> confRapidity{"confRapidity", 0.5, "Rapidity cut"};
  ConfigurableAxis configThnAxisInvMass{"configThnAxisInvMass", {180, 0.6, 1.5}, "#it{M} (GeV/#it{c}^{2})"};
  ConfigurableAxis configThnAxisPt{"configThnAxisPt", {100, 0.0, 10.}, "#it{p}_{T} (GeV/#it{c})"};
  ConfigurableAxis configThnAxisCentrality{"configThnAxisCentrality", {8, 0., 80}, "Centrality"};
  ConfigurableAxis configrapAxis{"configrapAxis", {VARIABLE_WIDTH, -0.8, -0.4, 0.4, 0.8}, "Rapidity"};
  Configurable<float> confFakeKaonCut{"confFakeKaonCut", 0.1, "Cut based on track from momentum difference"};
  ConfigurableAxis configThnAxisV2{"configThnAxisV2", {400, -16, 16}, "V2"};
  Configurable<bool> isNoTOF{"isNoTOF", true, "isNoTOF"};
  Configurable<bool> pdgcheck{"pdgcheck", true, "pdgcheck"};
  Configurable<int> strategyPID{"strategyPID", 2, "PID strategy"};
  Configurable<float> cfgCutTOFBeta{"cfgCutTOFBeta", 0.0, "cut TOF beta"};
  Configurable<bool> additionalQAplots{"additionalQAplots", true, "Additional QA plots"};
  Configurable<bool> additionalQAplots1{"additionalQAplots1", true, "Additional QA plots"};
  Configurable<float> confMinRot{"confMinRot", 5.0 * o2::constants::math::PI / 6.0, "Minimum of rotation"};
  Configurable<float> confMaxRot{"confMaxRot", 7.0 * o2::constants::math::PI / 6.0, "Maximum of rotation"};
  Configurable<int> nBkgRotations{"nBkgRotations", 9, "Number of rotated copies (background) per each original candidate"};
  Configurable<bool> fillRotation{"fillRotation", true, "fill rotation"};
  Configurable<bool> fillSA{"fillSA", true, "same event SA"};
  Configurable<bool> fillOccupancy{"fillOccupancy", false, "fill Occupancy"};
  Configurable<int> cfgOccupancyCut{"cfgOccupancyCut", 500, "Occupancy cut"};
  Configurable<bool> useWeight{"useWeight", false, "use EP dep effi weight"};
  Configurable<bool> useSP{"useSP", false, "use SP"};
  Configurable<float> cfgMinTrackPt{"cfgMinTrackPt", 0.15f,
                                    "Minimum track pT"};

  Configurable<float> cfgMaxTrackPt{"cfgMaxTrackPt", 10.0f,
                                    "Maximum track pT"};
  Configurable<bool> genacceptancecut{"genacceptancecut", true, "use acceptance cut for generated"};
  Configurable<bool> avoidsplitrackMC{"avoidsplitrackMC", false, "avoid split track in MC"};
  Configurable<bool> additionalEvSel1{"additionalEvSel1", true, "Additional evsel1"};
  Configurable<bool> additionalEvSel2{"additionalEvSel2", true, "Additional evsel2"};
  Configurable<bool> additionalEvSel3{"additionalEvSel3", true, "Additional evsel3"};
  Configurable<bool> additionalEvSel4{"additionalEvSel4", true, "Additional evsel4"};
  Configurable<std::string> confWeightPath{"confWeightPath", "Users/s/skundu/My/Object/fitweight", "Path to gain calibration"};
  ConfigurableAxis axisPtKaonWeight{"axisPtKaonWeight", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f}, "pt axis"};

  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter centralityFilter = nabs(aod::cent::centFT0C) < cfgCutCentrality;
  Filter acceptanceFilter = (nabs(aod::track::eta) < cfgCutEta && nabs(aod::track::pt) > cfgCutPT);
  Filter dcacutFilter = (nabs(aod::track::dcaXY) < cfgCutDCAxy) && (nabs(aod::track::dcaZ) < cfgCutDCAz);

  using EventCandidates = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::FV0Mults, aod::TPCMults, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0Cs, aod::CentFT0As, aod::EPCalibrationTables, aod::Mults>>;
  using TrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTOFbeta, aod::pidTPCFullKa, aod::pidTOFFullKa, aod::pidTPCFullPi, aod::pidTOFFullPi>>;

  using CollisionMCTrueTable = aod::McCollisions;
  using TrackMCTrueTable = aod::McParticles;
  using CollisionMCRecTableCentFT0C = soa::SmallGroups<soa::Join<aod::McCollisionLabels, aod::Collisions, aod::CentFT0Cs, aod::EvSels>>;
  using TrackMCRecTable = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::McTrackLabels, aod::TrackSelection, aod::pidTOFbeta, aod::pidTPCFullKa, aod::pidTOFFullKa, aod::pidTPCFullPi, aod::pidTOFFullPi>;
  using FilTrackMCRecTable = soa::Filtered<TrackMCRecTable>;

  Preslice<TrackMCRecTable> perCollision = aod::track::collisionId;

  SliceCache cache;
  // Partition<TrackCandidates> posTracks = aod::track::signed1Pt > cfgCutCharge;
  // Partition<TrackCandidates> negTracks = aod::track::signed1Pt < cfgCutCharge;

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // Event selection cuts - Alex
  TF1* fMultPVCutLow = nullptr;
  TF1* fMultPVCutHigh = nullptr;
  TF1* fMultCutLow = nullptr;
  TF1* fMultCutHigh = nullptr;
  TF1* fMultMultPVCut = nullptr;

  void init(o2::framework::InitContext&)
  {
    rctCut.rctChecker.init(
      rctCut.cfgEvtRCTFlagCheckerLabel,
      rctCut.cfgEvtRCTFlagCheckerZDCCheck,
      rctCut.cfgEvtRCTFlagCheckerLimitAcceptAsBad);

    std::vector<double> occupancyBinning = {0.0, 500.0, 1000.0, 1500.0, 2000.0, 3000.0, 4000.0, 5000.0, 50000.0};
    AxisSpec phiAxis = {500, -6.28, 6.28, "phi"};
    AxisSpec resAxis = {6000, -30, 30, "Res"};
    AxisSpec centAxis = {8, 0, 80, "V0M (%)"};
    AxisSpec occupancyAxis = {occupancyBinning, "Occupancy"};
    histos.add("hEvtSelInfo", "hEvtSelInfo", kTH1F, {{10, 0, 10.0}});
    if (!fillSA) {
      histos.add("hSparseV2SASameEvent_V2", "hSparseV2SASameEvent_V2", HistType::kTHnSparseF, {configThnAxisInvMass, configThnAxisPt, configThnAxisV2, configThnAxisCentrality});
    }
    if (fillRotation) {
      if (!fillSA) {
        histos.add("hRotation", "hRotation", kTH1F, {{360, 0.0, o2::constants::math::TwoPI}});
        histos.add("hSparseV2SASameEventRotational_V2", "hSparseV2SASameEventRotational_V2", HistType::kTHnSparseF, {configThnAxisInvMass, configThnAxisPt, configThnAxisV2, configThnAxisCentrality});
      }
    }

    if (fillSA) {
      histos.add("hSparseSAvsrapsameunlike", "hSparseSAvsrapsameunlike", HistType::kTHnSparseF, {configThnAxisInvMass, configThnAxisPt, configThnAxisV2, configrapAxis, configThnAxisCentrality}, true);
      histos.add("hSparseSAvsrapsamelike", "hSparseSAvsrapsamelike", HistType::kTHnSparseF, {configThnAxisInvMass, configThnAxisPt, configThnAxisV2, configrapAxis, configThnAxisCentrality}, true);
      histos.add("hSparseSAvsraprot", "hSparseSAvsraprot", HistType::kTHnSparseF, {configThnAxisInvMass, configThnAxisPt, configThnAxisV2, configrapAxis, configThnAxisCentrality}, true);
      histos.add("hSparseSAvsrapmix", "hSparseSAvsrapmix", HistType::kTHnSparseF, {configThnAxisInvMass, configThnAxisPt, configThnAxisV2, configrapAxis, configThnAxisCentrality}, true);
    }
    histos.add("hSparseV2SAGen_V2", "hSparseV2SAGen_V2", HistType::kTHnSparseF, {configThnAxisInvMass, configThnAxisPt, configThnAxisV2, configThnAxisCentrality});
    histos.add("hSparseV2SARec_V2", "hSparseV2SARec_V2", HistType::kTHnSparseF, {configThnAxisInvMass, configThnAxisPt, configThnAxisV2, configThnAxisCentrality});
    histos.add("hpt", "hpt", kTH1F, {configThnAxisPt});
    histos.add("hMC", "MC Event statistics", kTH1F, {{10, 0.0f, 10.0f}});
    histos.add("h1PhiRecsplit", "Phi meson Rec split", kTH1F, {{100, 0.0f, 10.0f}});
    histos.add("CentPercentileMCRecHist", "MC Centrality", kTH1F, {{100, 0.0f, 100.0f}});
    histos.add("hSparseV2SAMixedEvent_V2", "hSparseV2SAMixedEvent_V2", HistType::kTHnSparseF, {configThnAxisInvMass, configThnAxisPt, configThnAxisV2, configThnAxisCentrality});
    histos.add("h2PhiGen2", "Phi meson gen", kTH2F, {configThnAxisPt, configThnAxisCentrality});
    histos.add("h2PhiRec2", "Phi meson Rec", kTH2F, {configThnAxisPt, configThnAxisCentrality});
    histos.add("hImpactParameter", "Impact parameter", kTH1F, {{200, 0.0f, 20.0f}});
    histos.add("hEventPlaneAngle", "hEventPlaneAngle", kTH1F, {{200, -o2::constants::math::TwoPI, o2::constants::math::TwoPI}});
    histos.add("hSparseKstarMCGenWeight", "hSparseKstarMCGenWeight", HistType::kTHnSparseD, {configThnAxisCentrality, {36, 0.0f, o2::constants::math::PI}, {400, 0.0f, 1}, configThnAxisPt, {8, -0.8, 0.8}});
    histos.add("hSparseKstarMCRecWeight", "hSparseKstarMCRecWeight", HistType::kTHnSparseD, {configThnAxisCentrality, {36, 0.0f, o2::constants::math::PI}, {400, 0.0f, 1}, configThnAxisPt, {8, -0.8, 0.8}});
    histos.add("hSparseKstarMCGenKaonWeight", "hSparseKstarMCGenKaonWeight", HistType::kTHnSparseD, {configThnAxisCentrality, {36, 0.0f, o2::constants::math::PI}, {400, 0.0f, 1}, axisPtKaonWeight, {8, -0.8, 0.8}});
    histos.add("hSparseKstarMCRecKaonWeight", "hSparseKstarMCRecKaonWeight", HistType::kTHnSparseD, {configThnAxisCentrality, {36, 0.0f, o2::constants::math::PI}, {400, 0.0f, 1}, axisPtKaonWeight, {8, -0.8, 0.8}});
    histos.add("hSparseKstarMCRecKaonMissMatchWeight", "hSparseKstarMCRecKaonMissMatchWeight", HistType::kTHnSparseD, {configThnAxisCentrality, {36, 0.0f, o2::constants::math::PI}, {400, 0.0f, 1}, axisPtKaonWeight, {8, -0.8, 0.8}});
    histos.add("hSparseKstarMCGenPionWeight", "hSparseKstarMCGenPionWeight", HistType::kTHnSparseD, {configThnAxisCentrality, {36, 0.0f, o2::constants::math::PI}, {400, 0.0f, 1}, axisPtKaonWeight, {8, -0.8, 0.8}});
    histos.add("hSparseKstarMCRecPionWeight", "hSparseKstarMCRecPionWeight", HistType::kTHnSparseD, {configThnAxisCentrality, {36, 0.0f, o2::constants::math::PI}, {400, 0.0f, 1}, axisPtKaonWeight, {8, -0.8, 0.8}});
    histos.add("hSparseKstarMCRecPionMissMatchWeight", "hSparseKstarMCRecPionMissMatchWeight", HistType::kTHnSparseD, {configThnAxisCentrality, {36, 0.0f, o2::constants::math::PI}, {400, 0.0f, 1}, axisPtKaonWeight, {8, -0.8, 0.8}});
    histos.add("hSparseKstarMCGenSA", "hSparseKstarMCGenSA", HistType::kTHnSparseD, {configThnAxisInvMass, configThnAxisPt, configThnAxisV2, configrapAxis, configThnAxisCentrality});
    histos.add("hSparseKstarMCGenCosThetaStar_effy", "hSparseKstarMCGenCosThetaStar_effy", HistType::kTHnSparseD, {configThnAxisInvMass, configThnAxisPt, configThnAxisV2, configrapAxis, configThnAxisCentrality});
    histos.add("hSparseKstarMCRecSA", "hSparseKstarMCRecSA", HistType::kTHnSparseD, {configThnAxisInvMass, configThnAxisPt, configThnAxisV2, configrapAxis, configThnAxisCentrality});
    histos.add("hSparseKstarMCRecCosThetaStar_effy", "hSparseKstarMCRecCosThetaStar_effy", HistType::kTHnSparseD, {configThnAxisInvMass, configThnAxisPt, configThnAxisV2, configrapAxis, configThnAxisCentrality});
    if (additionalQAplots1) {
      histos.add("hFTOCvsTPCSelected", "Mult correlation FT0C vs. TPC after selection", kTH2F, {{80, 0.0f, 80.0f}, {100, -0.5f, 5999.5f}});
      histos.add("hCentrality", "Centrality distribution", kTH1F, {{200, 0.0, 200.0}});
      histos.add("hOccupancy", "Occupancy distribution", kTH1F, {occupancyAxis});
      histos.add("hVtxZ", "Vertex distribution in Z;Z (cm)", kTH1F, {{400, -20.0, 20.0}});
      histos.add("hPsiFT0C", "PsiFT0C", kTH2F, {centAxis, phiAxis});
      histos.add("hPsiFT0A", "PsiFT0A", kTH2F, {centAxis, phiAxis});
      histos.add("hPsiTPC", "PsiTPC", kTH2F, {centAxis, phiAxis});
      histos.add("ResFT0CTPC", "ResFT0CTPC", kTH2F, {centAxis, resAxis});
      histos.add("ResFT0CFT0A", "ResFT0CFT0A", kTH2F, {centAxis, resAxis});
      histos.add("ResFT0ATPC", "ResFT0ATPC", kTH2F, {centAxis, resAxis});
      histos.add("ResFT0CTPCSP", "ResFT0CTPCSP", kTH2F, {centAxis, resAxis});
      histos.add("ResFT0CFT0ASP", "ResFT0CFT0ASP", kTH2F, {centAxis, resAxis});
      histos.add("ResFT0ATPCSP", "ResFT0ATPCSP", kTH2F, {centAxis, resAxis});
      histos.add("ResTrackSPFT0CTPC", "ResTrackSPFT0CTPC", kTH3F, {centAxis, occupancyAxis, resAxis});
      histos.add("ResTrackSPFT0CFT0A", "ResTrackSPFT0CFT0A", kTH3F, {centAxis, occupancyAxis, resAxis});
      histos.add("ResTrackSPFT0ATPC", "ResTrackSPFT0ATPC", kTH3F, {centAxis, occupancyAxis, resAxis});
    }
    if (additionalQAplots) {
      // DCA QA
      histos.add("QAbefore/trkDCAxyka", "DCAxy distribution of kaon track candidates", HistType::kTH1F, {{150, -1.0f, 1.0f}});
      histos.add("QAbefore/trkDCAzka", "DCAz distribution of kaon track candidates", HistType::kTH1F, {{150, -1.0f, 1.0f}});
      histos.add("QAafter/trkDCAxyka", "DCAxy distribution of kaon track candidates", HistType::kTH1F, {{150, -1.0f, 1.0f}});
      histos.add("QAafter/trkDCAzka", "DCAz distribution of kaon track candidates", HistType::kTH1F, {{150, -1.0f, 1.0f}});

      // PID QA before cuts
      histos.add("QAbefore/TOF_TPC_Mapka_allka", "TOF + TPC Combined PID for Kaon;#sigma_{TOF}^{Kaon};#sigma_{TPC}^{Kaon}", {HistType::kTH2D, {{100, -6, 6}, {100, -6, 6}}});
      histos.add("QAbefore/TOF_Nsigma_allka", "TOF NSigma for Kaon;#it{p}_{T} (GeV/#it{c});#sigma_{TOF}^{Kaon};", {HistType::kTH3D, {{200, 0.0, 20.0}, {100, -6, 6}, {100, 0.0, 100.0}}});
      histos.add("QAbefore/TPC_Nsigma_allka", "TPC NSigma for Kaon;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Kaon};", {HistType::kTH3D, {{200, 0.0, 20.0}, {100, -6, 6}, {100, 0.0, 100.0}}});
      // PID QA after cuts
      histos.add("QAafter/TOF_TPC_Mapka_allka", "TOF + TPC Combined PID for Kaon;#sigma_{TOF}^{Kaon};#sigma_{TPC}^{Kaon}", {HistType::kTH2D, {{100, -6, 6}, {100, -6, 6}}});
      histos.add("QAafter/TOF_Nsigma_allka", "TOF NSigma for Kaon;#it{p}_{T} (GeV/#it{c});#sigma_{TOF}^{Kaon};", {HistType::kTH3D, {{200, 0.0, 20.0}, {100, -6, 6}, {100, 0.0, 100.0}}});
      histos.add("QAafter/TPC_Nsigma_allka", "TPC NSigma for Kaon;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Kaon};", {HistType::kTH3D, {{200, 0.0, 20.0}, {100, -6, 6}, {100, 0.0, 100.0}}});

      // DCA QA
      histos.add("QAbefore/trkDCAxypi", "DCAxy distribution of pion track candidates", HistType::kTH1F, {{150, -1.0f, 1.0f}});
      histos.add("QAbefore/trkDCAzpi", "DCAz distribution of pion track candidates", HistType::kTH1F, {{150, -1.0f, 1.0f}});
      histos.add("QAafter/trkDCAxypi", "DCAxy distribution of pion track candidates", HistType::kTH1F, {{150, -1.0f, 1.0f}});
      histos.add("QAafter/trkDCAzpi", "DCAz distribution of pion track candidates", HistType::kTH1F, {{150, -1.0f, 1.0f}});
      // PID QA before cuts
      histos.add("QAbefore/TOF_TPC_Mapka_allpi", "TOF + TPC Combined PID for pion;#sigma_{TOF}^{pion};#sigma_{TPC}^{pion}", {HistType::kTH2D, {{100, -6, 6}, {100, -6, 6}}});
      histos.add("QAbefore/TOF_Nsigma_allpi", "TOF NSigma for pion;#it{p}_{T} (GeV/#it{c});#sigma_{TOF}^{pion};", {HistType::kTH3D, {{200, 0.0, 20.0}, {100, -6, 6}, {100, 0.0, 100.0}}});
      histos.add("QAbefore/TPC_Nsigma_allpi", "TPC NSigma for pion;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{pion};", {HistType::kTH3D, {{200, 0.0, 20.0}, {100, -6, 6}, {100, 0.0, 100.0}}});
      // PID QA after cuts
      histos.add("QAafter/TOF_TPC_Mapka_allpi", "TOF + TPC Combined PID for pion;#sigma_{TOF}^{pion};#sigma_{TPC}^{pion}", {HistType::kTH2D, {{100, -6, 6}, {100, -6, 6}}});
      histos.add("QAafter/TOF_Nsigma_allpi", "TOF NSigma for pion;#it{p}_{T} (GeV/#it{c});#sigma_{TOF}^{pion};", {HistType::kTH3D, {{200, 0.0, 20.0}, {100, -6, 6}, {100, 0.0, 100.0}}});
      histos.add("QAafter/TPC_Nsigma_allpi", "TPC NSigma for pion;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{pion};", {HistType::kTH3D, {{200, 0.0, 20.0}, {100, -6, 6}, {100, 0.0, 100.0}}});
    }

    histos.add("hMassSameEventLikeNN", "Same-event like-sign (--) mass", kTH2F, {configThnAxisInvMass, configThnAxisCentrality});
    histos.add("hMassSameEventLikePP", "Same-event like-sign (++) mass", kTH2F, {configThnAxisInvMass, configThnAxisCentrality});
    histos.add("hMassMixedEventLikeNN", "Mixed-event like-sign (--) mass", kTH2F, {configThnAxisInvMass, configThnAxisCentrality});
    histos.add("hMassMixedEventLikePP", "Mixed-event like-sign (++) mass", kTH2F, {configThnAxisInvMass, configThnAxisCentrality});
    histos.add("hMassMixedEventUnlike", "Mixed-event unlike-sign mass", kTH2F, {configThnAxisInvMass, configThnAxisCentrality});

    ccdb->setURL(cfgCcdbParam.cfgURL);
    ccdbApi.init("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setCreatedNotAfter(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
  }

  double massKa = o2::constants::physics::MassKPlus;
  double massPi = o2::constants::physics::MassPiMinus;

  template <typename T>
  bool selectionTrack(const T& candidate)
  {
    if (useGlobalTrack && !(candidate.isGlobalTrack() && candidate.isPVContributor() && candidate.itsNCls() > cfgITScluster && candidate.tpcNClsFound() > cfgTPCcluster)) {
      return false;
    }
    if (!useGlobalTrack && !(candidate.isPVContributor() && candidate.itsNCls() > cfgITScluster)) {
      return false;
    }
    return true;
  }
  static constexpr float TPCOnlyPt = 0.5f;
  static constexpr int Strategy = 2;
  template <typename T>
  bool selectionPID2(const T& candidate, int PID)
  {
    if (PID == 0) {
      if (candidate.hasTOF() && candidate.beta() > cfgCutTOFBeta && std::abs(candidate.tofNSigmaKa()) < nsigmaCutTOF) {
        return true;
      }
    }
    if (PID == 1) {
      if (candidate.hasTOF() && candidate.beta() > cfgCutTOFBeta && std::abs(candidate.tofNSigmaPi()) < nsigmaCutTOF) {
        return true;
      }
    }
    return false;
  }
  template <typename T>
  bool strategySelectionPID(const T& candidate, int PID, int strategy)
  {
    if (PID == 0) {
      if (strategy == 0) {
        if (!isNoTOF && !candidate.hasTOF() && std::abs(candidate.tpcNSigmaKa()) < nsigmaCutTPC) {
          return true;
        }
        if (!isNoTOF && candidate.hasTOF() && std::abs(candidate.tpcNSigmaKa()) < nsigmaCutTPC && std::abs(candidate.tofNSigmaKa()) < nsigmaCutTOF) {
          return true;
        }
        if (isNoTOF && std::abs(candidate.tpcNSigmaKa()) < nsigmaCutTPC) {
          return true;
        }
      } else if (strategy == 1) {
        if (!isNoTOF && !candidate.hasTOF() && std::abs(candidate.tpcNSigmaKa()) < nsigmaCutTPC) {
          return true;
        }
        if (!isNoTOF && candidate.hasTOF() && ((candidate.tofNSigmaKa() * candidate.tofNSigmaKa()) + (candidate.tpcNSigmaKa() * candidate.tpcNSigmaKa())) < (nsigmaCutCombined * nsigmaCutCombined)) {
          return true;
        }
        if (isNoTOF && std::abs(candidate.tpcNSigmaKa()) < nsigmaCutTPC) {
          return true;
        }
      } else if (strategy == Strategy) {
        if (candidate.pt() < TPCOnlyPt && std::abs(candidate.tpcNSigmaKa()) < nsigmaCutTPC) {
          return true;
        }
        if (candidate.pt() >= TPCOnlyPt && std::abs(candidate.tpcNSigmaKa()) < nsigmaCutTPC && candidate.hasTOF() && std::abs(candidate.tofNSigmaKa()) < nsigmaCutTOF && candidate.beta() > cfgCutTOFBeta) {
          return true;
        }
        if (candidate.pt() >= TPCOnlyPt && std::abs(candidate.tpcNSigmaKa()) < nsigmaCutTPC && !candidate.hasTOF()) {
          return true;
        }
      }
    }
    if (PID == 1) {
      if (strategy == 0) {
        if (!isNoTOF && !candidate.hasTOF() && std::abs(candidate.tpcNSigmaPi()) < nsigmaCutTPC) {
          return true;
        }
        if (!isNoTOF && candidate.hasTOF() && std::abs(candidate.tpcNSigmaPi()) < nsigmaCutTPC && std::abs(candidate.tofNSigmaPi()) < nsigmaCutTOF) {
          return true;
        }
        if (isNoTOF && std::abs(candidate.tpcNSigmaPi()) < nsigmaCutTPC) {
          return true;
        }
      } else if (strategy == 1) {
        if (!isNoTOF && !candidate.hasTOF() && std::abs(candidate.tpcNSigmaPi()) < nsigmaCutTPC) {
          return true;
        }
        if (!isNoTOF && candidate.hasTOF() && ((candidate.tofNSigmaPi() * candidate.tofNSigmaPi()) + (candidate.tpcNSigmaPi() * candidate.tpcNSigmaPi())) < (nsigmaCutCombined * nsigmaCutCombined)) {
          return true;
        }
        if (isNoTOF && std::abs(candidate.tpcNSigmaPi()) < nsigmaCutTPC) {
          return true;
        }
      } else if (strategy == Strategy) {
        if (candidate.pt() < TPCOnlyPt && std::abs(candidate.tpcNSigmaPi()) < nsigmaCutTPC) {
          return true;
        }
        if (candidate.pt() >= TPCOnlyPt && std::abs(candidate.tpcNSigmaPi()) < nsigmaCutTPC && candidate.hasTOF() && std::abs(candidate.tofNSigmaPi()) < nsigmaCutTOF && candidate.beta() > cfgCutTOFBeta) {
          return true;
        }
        if (candidate.pt() >= TPCOnlyPt && std::abs(candidate.tpcNSigmaPi()) < nsigmaCutTPC && !candidate.hasTOF()) {
          return true;
        }
      }
    }
    return false;
  }

  double getPhiInRange(double phi)
  {
    double result = phi;
    while (result < 0) {
      result += o2::constants::math::PI;
    }
    while (result >= o2::constants::math::PI) { // >= not >
      result -= o2::constants::math::PI;
    }
    return result;
  }
  template <typename T>
  bool isFakeKaon(T const& track, int /*PID*/)
  {
    const auto pglobal = track.p();
    const auto ptpc = track.tpcInnerParam();
    return std::abs(pglobal - ptpc) > confFakeKaonCut;
  }

  static constexpr float HalfPI = o2::constants::math::PI * 0.5f;
  ConfigurableAxis axisVertex{"axisVertex", {20, -10, 10}, "vertex axis for bin"};
  ConfigurableAxis axisMultiplicityClass{"axisMultiplicityClass", {20, 0, 100}, "multiplicity percentile for bin"};
  ConfigurableAxis axisEPAngle{"axisEPAngle",
                               {6, -HalfPI, HalfPI},
                               "event plane angle"};
  ConfigurableAxis axisOccup{"axisOccup", {20, -0.5, 40000.0}, "occupancy axis"};
  double v2 = 0.;
  double v2Rot = 0.;

  using BinningTypeVertexContributor = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0C, aod::epcalibrationtable::PsiFT0C>;
  ROOT::Math::PxPyPzMVector kstarMother, fourVecDauCM, daughter1, daughter2, kaonrot, kstarrot, kaonPlus, pionMinus;
  ROOT::Math::XYZVector threeVecDauCM, threeVecDauCMXY, eventplaneVec, eventplaneVecNorm;
  ROOT::Math::PxPyPzMVector daughter2rot, fourVecDauCMrot;
  ROOT::Math::XYZVector threeVecDauCMrot, threeVecDauCMXYrot;

  int currentRunNumber = -999;
  int lastRunNumber = -999;
  TH2D* hweight = nullptr;
  void processSE(EventCandidates::iterator const& collision, TrackCandidates const& tracks, aod::BCsWithTimestamps const&)
  {
    histos.fill(HIST("hEvtSelInfo"), 0.5);
    if (rctCut.requireRCTFlagChecker && !rctCut.rctChecker(collision)) {
      return;
    }
    histos.fill(HIST("hEvtSelInfo"), 1.5);
    if (!collision.sel8()) {
      return;
    }
    if (!collision.triggereventep()) {
      return;
    }
    if (additionalEvSel1 &&
        !collision.selection_bit(aod::evsel::kNoTimeFrameBorder)) {
      return;
    }
    if (additionalEvSel2 &&
        !collision.selection_bit(aod::evsel::kNoITSROFrameBorder)) {
      return;
    }
    if (additionalEvSel3 &&
        !collision.selection_bit(aod::evsel::kNoSameBunchPileup)) {
      return;
    }
    if (additionalEvSel4 &&
        !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return;
    }
    histos.fill(HIST("hEvtSelInfo"), 2.5);
    auto centrality = collision.centFT0C();
    auto multTPC = collision.multNTracksPV();
    int occupancy = collision.trackOccupancyInTimeRange();
    auto psiFT0C = collision.psiFT0C();
    auto psiFT0A = collision.psiFT0A();
    auto psiTPC = collision.psiTPC();
    auto qFT0C = collision.qFT0C();
    auto qFT0A = collision.qFT0A();
    auto qTPC = collision.qTPC();
    if (fillOccupancy && occupancy > cfgOccupancyCut) {
      return;
    }
    histos.fill(HIST("hEvtSelInfo"), 3.5);
    if (additionalQAplots1) {
      histos.fill(HIST("hFTOCvsTPCSelected"), centrality, multTPC);
      histos.fill(HIST("hPsiFT0C"), centrality, psiFT0C);
      histos.fill(HIST("hPsiFT0A"), centrality, psiFT0A);
      histos.fill(HIST("hPsiTPC"), centrality, psiTPC);
      histos.fill(HIST("ResFT0CTPC"), centrality, std::cos(2.0 * (psiFT0C - psiTPC)));
      histos.fill(HIST("ResFT0CFT0A"), centrality, std::cos(2.0 * (psiFT0C - psiFT0A)));
      histos.fill(HIST("ResFT0ATPC"), centrality, std::cos(2.0 * (psiTPC - psiFT0A)));
      histos.fill(HIST("ResFT0CTPCSP"), centrality, qFT0C * qTPC * std::cos(2.0 * (psiFT0C - psiTPC)));
      histos.fill(HIST("ResFT0CFT0ASP"), centrality, qFT0C * qFT0A * std::cos(2.0 * (psiFT0C - psiFT0A)));
      histos.fill(HIST("ResFT0ATPCSP"), centrality, qTPC * qFT0A * std::cos(2.0 * (psiTPC - psiFT0A)));
      histos.fill(HIST("hCentrality"), centrality);
      histos.fill(HIST("hOccupancy"), occupancy);
      histos.fill(HIST("hVtxZ"), collision.posZ());
    }
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    currentRunNumber = collision.bc_as<aod::BCsWithTimestamps>().runNumber();
    if (useWeight && (currentRunNumber != lastRunNumber)) {
      hweight = ccdb->getForTimeStamp<TH2D>(confWeightPath.value, bc.timestamp());
    }
    lastRunNumber = currentRunNumber;
    float weight1 = 1.0;
    float weight2 = 1.0;
    for (const auto& track1 : tracks) {
      if (!selectionTrack(track1)) {
        continue;
      }

      auto track1ID = track1.globalIndex();
      if (!isTOFOnly && !strategySelectionPID(track1, 0, strategyPID)) {
        continue;
      }
      if (isTOFOnly && !selectionPID2(track1, 0)) {
        continue;
      }

      if (useWeight) {
        if (track1.pt() < cfgMaxTrackPt &&
            track1.pt() > cfgMinTrackPt) {
          weight1 = 1 + hweight->GetBinContent(hweight->FindBin(centrality, track1.pt() + 0.000005)) * std::cos(2.0 * getPhiInRange(track1.phi() - psiFT0C));
        } else {
          weight1 = 1;
        }
      }
      for (const auto& track2 : tracks) {
        if (!selectionTrack(track2)) {
          continue;
        }

        auto track2ID = track2.globalIndex();
        if (!isTOFOnly && !strategySelectionPID(track2, 1, strategyPID)) {
          continue;
        }
        if (isTOFOnly && !selectionPID2(track2, 1)) {
          continue;
        }

        if (track2ID == track1ID) {
          continue;
        }

        if (additionalQAplots) {
          histos.fill(HIST("QAafter/TPC_Nsigma_allka"), track1.pt(), track1.tpcNSigmaKa(), centrality);
          histos.fill(HIST("QAafter/TOF_Nsigma_allka"), track1.pt(), track1.tofNSigmaKa(), centrality);
          histos.fill(HIST("QAafter/trkDCAxyka"), track1.dcaXY());
          histos.fill(HIST("QAafter/trkDCAzka"), track1.dcaZ());
          histos.fill(HIST("QAafter/TOF_TPC_Mapka_allka"), track1.tofNSigmaKa(), track1.tpcNSigmaKa());
          histos.fill(HIST("QAafter/TOF_TPC_Mapka_allpi"), track2.tofNSigmaPi(), track2.tpcNSigmaPi());
          histos.fill(HIST("QAafter/TPC_Nsigma_allpi"), track2.pt(), track2.tpcNSigmaPi(), centrality);
          histos.fill(HIST("QAafter/TOF_Nsigma_allpi"), track2.pt(), track2.tofNSigmaPi(), centrality);
          histos.fill(HIST("QAafter/trkDCAxypi"), track2.dcaXY());
          histos.fill(HIST("QAafter/trkDCAzpi"), track2.dcaZ());
        }
        if (useWeight) {
          if (track2.pt() < cfgMaxTrackPt &&
              track2.pt() > cfgMinTrackPt) {
            weight2 = 1 + hweight->GetBinContent(hweight->FindBin(centrality, track2.pt() + 0.000005)) * std::cos(2.0 * getPhiInRange(track2.phi() - psiFT0C));
          } else {
            weight2 = 1;
          }
        }
        daughter1 = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), massKa);
        daughter2 = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), massPi);
        kstarMother = daughter1 + daughter2;
        if (std::abs(kstarMother.Rapidity()) > confRapidity) {
          continue;
        }
        auto phiMinusPsi = getPhiInRange(kstarMother.Phi() - psiFT0C);

        if (useSP) {
          v2 = std::cos(2.0 * phiMinusPsi) * qFT0C;
        }
        if (!useSP) {
          v2 = std::cos(2.0 * phiMinusPsi);
        }
        auto totalweight = weight1 * weight2;
        static constexpr float MinTotalWeight = 5e-7f;
        if (totalweight <= MinTotalWeight) {
          totalweight = 1.0;
        }
        if (additionalQAplots1) {
          histos.fill(HIST("ResTrackSPFT0CTPC"), centrality, occupancy, qFT0C * qTPC * std::cos(2.0 * (psiFT0C - psiTPC)));
          histos.fill(HIST("ResTrackSPFT0CFT0A"), centrality, occupancy, qFT0C * qFT0A * std::cos(2.0 * (psiFT0C - psiFT0A)));
          histos.fill(HIST("ResTrackSPFT0ATPC"), centrality, occupancy, qTPC * qFT0A * std::cos(2.0 * (psiTPC - psiFT0A)));
        }
        if (!fillSA) {

          if (useWeight) {
            histos.fill(HIST("hSparseV2SASameEvent_V2"), kstarMother.M(), kstarMother.Pt(), v2, centrality, 1 / totalweight);
          } else {
            histos.fill(HIST("hSparseV2SASameEvent_V2"), kstarMother.M(), kstarMother.Pt(), v2, centrality);
          }
        }
        int track1Sign = track1.sign();
        int track2Sign = track2.sign();

        if (track1Sign * track2Sign > 0) {
          if (track1Sign > 0) {
            histos.fill(HIST("hMassSameEventLikePP"), kstarMother.M(), centrality);
          } else {
            histos.fill(HIST("hMassSameEventLikeNN"), kstarMother.M(), centrality);
          }
        } else {
        }

        if (fillSA) {
          ROOT::Math::Boost boost{kstarMother.BoostToCM()};
          fourVecDauCM = boost(daughter1);
          threeVecDauCM = fourVecDauCM.Vect();
          threeVecDauCMXY = ROOT::Math::XYZVector(threeVecDauCM.X(), threeVecDauCM.Y(), 0.);
          eventplaneVec = ROOT::Math::XYZVector(std::cos(2.0 * psiFT0C), std::sin(2.0 * psiFT0C), 0);
          eventplaneVecNorm = ROOT::Math::XYZVector(std::sin(2.0 * psiFT0C), -std::cos(2.0 * psiFT0C), 0);
          auto cosPhistarminuspsi = getPhiInRange(fourVecDauCM.Phi() - psiFT0C);
          auto sa = std::cos(2.0 * cosPhistarminuspsi);
          auto cosThetaStar = eventplaneVecNorm.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2()) / std::sqrt(eventplaneVecNorm.Mag2());

          if (track1Sign * track2Sign < 0) {
            if (usepolar) {
              histos.fill(HIST("hSparseSAvsrapsameunlike"), kstarMother.M(), kstarMother.Pt(), cosThetaStar, kstarMother.Rapidity(), centrality);
            } else {
              histos.fill(HIST("hSparseSAvsrapsameunlike"), kstarMother.M(), kstarMother.Pt(), sa, kstarMother.Rapidity(), centrality);
            }
          } else if (track1Sign * track2Sign > 0) {
            if (usepolar) {
              histos.fill(HIST("hSparseSAvsrapsamelike"), kstarMother.M(), kstarMother.Pt(), cosThetaStar, kstarMother.Rapidity(), centrality);
            } else {
              histos.fill(HIST("hSparseSAvsrapsamelike"), kstarMother.M(), kstarMother.Pt(), sa, kstarMother.Rapidity(), centrality);
            }
          }
        }
        if (fillRotation) {
          for (int nrotbkg = 0; nrotbkg < nBkgRotations; nrotbkg++) {
            auto anglestart = confMinRot;
            auto angleend = confMaxRot;
            auto anglestep = (angleend - anglestart) / (1.0 * (nBkgRotations - 1));
            auto rotangle = anglestart + nrotbkg * anglestep;
            if (!fillSA) {
              histos.fill(HIST("hRotation"), rotangle);
            }
            auto rotkaonPx = track1.px() * std::cos(rotangle) - track1.py() * std::sin(rotangle);
            auto rotkaonPy = track1.px() * std::sin(rotangle) + track1.py() * std::cos(rotangle);
            kaonrot = ROOT::Math::PxPyPzMVector(rotkaonPx, rotkaonPy, track1.pz(), massKa);
            kstarrot = kaonrot + daughter2;
            if (std::abs(kstarrot.Rapidity()) > confRapidity) {
              continue;
            }
            auto phiMinusPsiRot = getPhiInRange(kstarrot.Phi() - psiFT0C);

            if (useSP) {
              v2Rot = std::cos(2.0 * phiMinusPsiRot) * qFT0C;
            }
            if (!useSP) {
              v2Rot = std::cos(2.0 * phiMinusPsiRot);
            }
            if (!fillSA) {
              histos.fill(HIST("hSparseV2SASameEventRotational_V2"), kstarrot.M(), kstarrot.Pt(), v2Rot, centrality);
            }
            if (fillSA) {
              if (track1Sign * track2Sign < 0) {
                ROOT::Math::Boost boost{kstarrot.BoostToCM()};
                fourVecDauCMrot = boost(kaonrot);
                threeVecDauCMrot = fourVecDauCMrot.Vect();
                threeVecDauCMXYrot = ROOT::Math::XYZVector(threeVecDauCMrot.X(), threeVecDauCMrot.Y(), 0.);
                auto cosPhistarminuspsirot = getPhiInRange(fourVecDauCMrot.Phi() - psiFT0C);
                auto sarot = std::cos(2.0 * cosPhistarminuspsirot);
                auto cosThetaStarrot = eventplaneVecNorm.Dot(threeVecDauCMrot) / std::sqrt(threeVecDauCMrot.Mag2()) / std::sqrt(eventplaneVecNorm.Mag2());
                if (usepolar) {
                  histos.fill(HIST("hSparseSAvsraprot"), kstarrot.M(), kstarrot.Pt(), cosThetaStarrot, kstarrot.Rapidity(), centrality);
                } else {
                  histos.fill(HIST("hSparseSAvsraprot"), kstarrot.M(), kstarrot.Pt(), sarot, kstarrot.Rapidity(), centrality);
                }
              }
            }
          }
        }
      }
    }
  }
  PROCESS_SWITCH(Kstarpbpb, processSE, "Process Same event latest", true);

  void processMixedEvent(EventCandidates const& collisions, TrackCandidates const& tracks)
  {

    auto tracksTuple = std::make_tuple(tracks);
    BinningTypeVertexContributor binningOnPositions{{axisVertex, axisMultiplicityClass, axisEPAngle}, true};
    SameKindPair<EventCandidates, TrackCandidates, BinningTypeVertexContributor> pair{binningOnPositions, cfgNoMixedEvents, -1, collisions, tracksTuple, &cache};
    for (const auto& [collision1, tracks1, collision2, tracks2] : pair) {
      if (rctCut.requireRCTFlagChecker && !rctCut.rctChecker(collision1)) {
        continue;
      }
      if (rctCut.requireRCTFlagChecker && !rctCut.rctChecker(collision2)) {
        continue;
      }

      if (!collision1.sel8() ||
          !collision2.sel8() ||

          !collision1.triggereventep() ||
          !collision2.triggereventep() ||

          (additionalEvSel1 &&
           (!collision1.selection_bit(aod::evsel::kNoTimeFrameBorder) ||
            !collision2.selection_bit(aod::evsel::kNoTimeFrameBorder))) ||

          (additionalEvSel2 &&
           (!collision1.selection_bit(aod::evsel::kNoITSROFrameBorder) ||
            !collision2.selection_bit(aod::evsel::kNoITSROFrameBorder))) ||

          (additionalEvSel3 &&
           (!collision1.selection_bit(aod::evsel::kNoSameBunchPileup) ||
            !collision2.selection_bit(aod::evsel::kNoSameBunchPileup))) ||

          (additionalEvSel4 &&
           (!collision1.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV) ||
            !collision2.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)))) {
        continue;
      }
      if (collision1.bcId() == collision2.bcId()) {
        continue;
      }
      int occupancy1 = collision1.trackOccupancyInTimeRange();
      int occupancy2 = collision2.trackOccupancyInTimeRange();
      if (fillOccupancy && occupancy1 >= cfgOccupancyCut && occupancy2 >= cfgOccupancyCut) // occupancy info is available for this collision (*)
      {
        continue;
      }
      auto centrality = collision1.centFT0C();
      auto psiFT0C1 = collision1.psiFT0C();
      auto qFT0C1 = collision1.qFT0C();
      auto psiFT0C2 = collision2.psiFT0C();
      auto qFT0C2 = collision2.qFT0C();

      for (const auto& [track1, track2] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {

        if (!selectionTrack(track1) || !selectionTrack(track2)) {
          continue;
        }
        if (!isTOFOnly && !strategySelectionPID(track1, 0, strategyPID)) {
          continue;
        }
        if (!isTOFOnly && !strategySelectionPID(track2, 1, strategyPID)) {
          continue;
        }
        if (isTOFOnly && !selectionPID2(track1, 0)) {
          continue;
        }
        if (isTOFOnly && !selectionPID2(track2, 1)) {
          continue;
        }

        daughter1 = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), massKa);
        daughter2 = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), massPi);

        kstarMother = daughter1 + daughter2;
        if (std::abs(kstarMother.Rapidity()) > confRapidity) {
          continue;
        }

        int s1 = track1.sign();
        int s2 = track2.sign();

        if (s1 * s2 < 0) {

          auto phi1 = track1.phi();
          auto phi2 = track2.phi();
          auto phiKstar = kstarMother.Phi();

          double term1 = qFT0C1 * std::cos(2.0 * getPhiInRange(phi1 - psiFT0C1)) * std::cos(2.0 * getPhiInRange(phi1 - phiKstar));
          double term2 = qFT0C2 * std::cos(2.0 * getPhiInRange(phi2 - psiFT0C2)) * std::cos(2.0 * getPhiInRange(phi2 - phiKstar));

          v2 = term1 + term2;

          if (!fillSA) {
            histos.fill(HIST("hSparseV2SAMixedEvent_V2"), kstarMother.M(), kstarMother.Pt(), v2, centrality);
          }
          histos.fill(HIST("hMassMixedEventUnlike"), kstarMother.M(), centrality);

          if (fillSA) {
            ROOT::Math::Boost boost{kstarMother.BoostToCM()};
            fourVecDauCM = boost(daughter1);
            threeVecDauCM = fourVecDauCM.Vect();
            threeVecDauCMXY = ROOT::Math::XYZVector(threeVecDauCM.X(), threeVecDauCM.Y(), 0.);
            eventplaneVec = ROOT::Math::XYZVector(std::cos(2.0 * psiFT0C1), std::sin(2.0 * psiFT0C1), 0);
            eventplaneVecNorm = ROOT::Math::XYZVector(std::sin(2.0 * psiFT0C1), -std::cos(2.0 * psiFT0C1), 0);
            auto cosPhistarminuspsi = getPhiInRange(fourVecDauCM.Phi() - psiFT0C1);
            auto sa = std::cos(2.0 * cosPhistarminuspsi);
            auto cosThetaStar = eventplaneVecNorm.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2()) / std::sqrt(eventplaneVecNorm.Mag2());
            if (usepolar) {
              histos.fill(HIST("hSparseSAvsrapmix"), kstarMother.M(), kstarMother.Pt(), cosThetaStar, kstarMother.Rapidity(), centrality);
            } else {
              histos.fill(HIST("hSparseSAvsrapmix"), kstarMother.M(), kstarMother.Pt(), sa, kstarMother.Rapidity(), centrality);
            }
          }
        } else {

          if (s1 > 0) {
            histos.fill(HIST("hMassMixedEventLikePP"), kstarMother.M(), centrality);
          } else {
            histos.fill(HIST("hMassMixedEventLikeNN"), kstarMother.M(), centrality);
          }
        }
      }
    }
  }
  PROCESS_SWITCH(Kstarpbpb, processMixedEvent, "Process Mixed event", true);

  void processMC(CollisionMCTrueTable::iterator const& /*TrueCollision*/, CollisionMCRecTableCentFT0C const& RecCollisions, TrackMCTrueTable const& GenParticles, FilTrackMCRecTable const& RecTracks)
  {
    histos.fill(HIST("hMC"), 0);
    if (RecCollisions.size() == 0) {
      histos.fill(HIST("hMC"), 1);
      return;
    }
    if (RecCollisions.size() > 1) {
      histos.fill(HIST("hMC"), 2);
      return;
    }
    for (const auto& RecCollision : RecCollisions) {
      auto psiFT0C = 0.0;
      histos.fill(HIST("hMC"), 3);
      if (!RecCollision.sel8()) {
        histos.fill(HIST("hMC"), 4);
        continue;
      }

      if (additionalEvSel1 &&
          !RecCollision.selection_bit(aod::evsel::kNoTimeFrameBorder)) {
        continue;
      }

      if (additionalEvSel2 &&
          !RecCollision.selection_bit(aod::evsel::kNoITSROFrameBorder)) {
        continue;
      }

      if (additionalEvSel3 &&
          !RecCollision.selection_bit(aod::evsel::kNoSameBunchPileup)) {
        continue;
      }

      if (additionalEvSel4 &&
          !RecCollision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
        continue;
      }
      histos.fill(HIST("hMC"), 5);
      if (std::abs(RecCollision.posZ()) > cfgCutVertex) {
        histos.fill(HIST("hMC"), 6);
        continue;
      }
      histos.fill(HIST("hMC"), 7);
      auto centrality = RecCollision.centFT0C();
      histos.fill(HIST("CentPercentileMCRecHist"), centrality);
      auto oldindex = -999;
      auto rectrackspart = RecTracks.sliceBy(perCollision, RecCollision.globalIndex());
      // loop over reconstructed particle
      for (const auto& track1 : rectrackspart) {
        if (!selectionTrack(track1)) {
          continue;
        }
        if (!isTOFOnly && !strategySelectionPID(track1, 0, strategyPID)) {
          continue;
        }
        if (isTOFOnly && !selectionPID2(track1, 0)) {
          continue;
        }
        if (!track1.has_mcParticle()) {
          continue;
        }
        auto track1ID = track1.index();
        for (const auto& track2 : rectrackspart) {
          auto track2ID = track2.index();
          if (track2ID <= track1ID) {
            continue;
          }
          if (!selectionTrack(track2)) {
            continue;
          }
          if (!isTOFOnly && !strategySelectionPID(track2, 1, strategyPID)) {
            continue;
          }
          if (isTOFOnly && !selectionPID2(track2, 1)) {
            continue;
          }
          if (!track2.has_mcParticle()) {
            continue;
          }
          if (track1.sign() * track2.sign() > 0) {
            continue;
          }
          const auto mctrack1 = track1.mcParticle();
          const auto mctrack2 = track2.mcParticle();
          int track1PDG = std::abs(mctrack1.pdgCode());
          int track2PDG = std::abs(mctrack2.pdgCode());
          if (!mctrack1.isPhysicalPrimary()) {
            continue;
          }
          if (!mctrack2.isPhysicalPrimary()) {
            continue;
          }
          if (track1PDG != PDG_t::kKPlus || track2PDG != PDG_t::kPiPlus) {
            continue;
          }
          for (const auto& mothertrack1 : mctrack1.mothers_as<aod::McParticles>()) {
            for (const auto& mothertrack2 : mctrack2.mothers_as<aod::McParticles>()) {
              if (mothertrack1.pdgCode() != mothertrack2.pdgCode()) {
                continue;
              }
              if (mothertrack1 != mothertrack2) {
                continue;
              }
              if (std::abs(mothertrack1.y()) > confRapidity) {
                continue;
              }
              if (pdgcheck && std::abs(mothertrack1.pdgCode()) != o2::constants::physics::kK0Star892) {
                continue;
              }
              if (!isTOFOnly && !(strategySelectionPID(track1, 0, strategyPID) || strategySelectionPID(track2, 1, strategyPID))) {
                continue;
              }
              if (isTOFOnly && !(selectionPID2(track1, 0) || selectionPID2(track2, 1))) {
                continue;
              }
              if (avoidsplitrackMC && oldindex == mothertrack1.globalIndex()) {
                continue;
              }
              oldindex = mothertrack1.globalIndex();
              if (track1.sign() > 0 && track2.sign() < 0) {
                kaonPlus = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), massKa);
                pionMinus = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), massPi);
              }
              if (track1.sign() < 0 && track2.sign() > 0) {
                pionMinus = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), massKa);
                kaonPlus = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), massPi);
              }
              kstarMother = kaonPlus + pionMinus;
              if (std::abs(kstarMother.Rapidity()) > confRapidity) {
                continue;
              }
              auto phiMinusPsi = getPhiInRange(kstarMother.Phi() - psiFT0C);

              v2 = std::cos(2.0 * phiMinusPsi);

              histos.fill(HIST("hSparseV2SARec_V2"), kstarMother.M(), kstarMother.Pt(), v2, centrality);
              histos.fill(HIST("h2PhiRec2"), kstarMother.pt(), centrality);
              histos.fill(HIST("hpt"), kstarMother.Pt());

              {
                ROOT::Math::Boost boost{kstarMother.BoostToCM()};
                auto fourVecDauCMRec = boost(kaonPlus);
                auto threeVecDauCMRec = fourVecDauCMRec.Vect();
                auto eventplaneVecNormRec = ROOT::Math::XYZVector(std::sin(2.0 * psiFT0C), -std::cos(2.0 * psiFT0C), 0);
                auto cosPhistarminuspsiRec = getPhiInRange(fourVecDauCMRec.Phi() - psiFT0C);
                auto saRec = std::cos(2.0 * cosPhistarminuspsiRec);
                auto cosThetaStarRec = eventplaneVecNormRec.Dot(threeVecDauCMRec) / std::sqrt(threeVecDauCMRec.Mag2()) / std::sqrt(eventplaneVecNormRec.Mag2());

                histos.fill(HIST("hSparseKstarMCRecSA"), kstarMother.M(), kstarMother.Pt(), saRec, std::abs(kstarMother.Rapidity()), centrality);
                histos.fill(HIST("hSparseKstarMCRecCosThetaStar_effy"), kstarMother.M(), kstarMother.Pt(), cosThetaStarRec, std::abs(kstarMother.Rapidity()), centrality);
              }
            }
          }
        }
      }
      // loop over generated particle
      for (const auto& mcParticle : GenParticles) {
        if (std::abs(mcParticle.y()) > confRapidity) {
          continue;
        }
        if (pdgcheck && mcParticle.pdgCode() != o2::constants::physics::kK0Star892) {
          continue;
        }
        auto kDaughters = mcParticle.daughters_as<aod::McParticles>();
        static constexpr std::size_t NumberOfDaughters = 2;

        if (kDaughters.size() != NumberOfDaughters) {
          continue;
        }
        auto daughtp = false;
        auto daughtm = false;
        for (const auto& kCurrentDaughter : kDaughters) {
          if (!kCurrentDaughter.isPhysicalPrimary()) {
            continue;
          }
          if (kCurrentDaughter.pdgCode() == +PDG_t::kKPlus) {
            if (genacceptancecut && kCurrentDaughter.pt() > cfgCutPT && std::abs(kCurrentDaughter.eta()) < cfgCutEta) {
              daughtp = true;
            }
            if (!genacceptancecut) {
              daughtp = true;
            }
            kaonPlus = ROOT::Math::PxPyPzMVector(kCurrentDaughter.px(), kCurrentDaughter.py(), kCurrentDaughter.pz(), massKa);
          } else if (kCurrentDaughter.pdgCode() == -PDG_t::kPiPlus) {
            if (genacceptancecut && kCurrentDaughter.pt() > cfgCutPT && std::abs(kCurrentDaughter.eta()) < cfgCutEta) {
              daughtm = true;
            }
            if (!genacceptancecut) {
              daughtm = true;
            }
            pionMinus = ROOT::Math::PxPyPzMVector(kCurrentDaughter.px(), kCurrentDaughter.py(), kCurrentDaughter.pz(), massPi);
          }
        }
        if (daughtp && daughtm) {
          kstarMother = kaonPlus + pionMinus;
          if (std::abs(kstarMother.Rapidity()) > confRapidity) {
            continue;
          }
          auto phiMinusPsi = getPhiInRange(kstarMother.Phi() - psiFT0C);

          v2 = std::cos(2.0 * phiMinusPsi);

          histos.fill(HIST("hSparseV2SAGen_V2"), kstarMother.M(), kstarMother.Pt(), v2, centrality);
          histos.fill(HIST("h2PhiGen2"), kstarMother.pt(), centrality);

          {
            ROOT::Math::Boost boost{kstarMother.BoostToCM()};
            auto fourVecDauCMGen = boost(kaonPlus);
            auto threeVecDauCMGen = fourVecDauCMGen.Vect();
            auto eventplaneVecNormGen = ROOT::Math::XYZVector(std::sin(2.0 * psiFT0C), -std::cos(2.0 * psiFT0C), 0);
            auto cosPhistarminuspsiGen = getPhiInRange(fourVecDauCMGen.Phi() - psiFT0C);
            auto saGen = std::cos(2.0 * cosPhistarminuspsiGen);
            auto cosThetaStarGen = eventplaneVecNormGen.Dot(threeVecDauCMGen) / std::sqrt(threeVecDauCMGen.Mag2()) / std::sqrt(eventplaneVecNormGen.Mag2());

            histos.fill(HIST("hSparseKstarMCGenSA"), kstarMother.M(), kstarMother.Pt(), saGen, std::abs(kstarMother.Rapidity()), centrality);
            histos.fill(HIST("hSparseKstarMCGenCosThetaStar_effy"), kstarMother.M(), kstarMother.Pt(), cosThetaStarGen, std::abs(kstarMother.Rapidity()), centrality);
          }
        }
      }
    } // rec collision loop

  } // process MC
  PROCESS_SWITCH(Kstarpbpb, processMC, "Process MC", false);

  void processMCkstarWeight(CollisionMCTrueTable::iterator const& TrueCollision, CollisionMCRecTableCentFT0C const& RecCollisions, TrackMCTrueTable const& GenParticles, FilTrackMCRecTable const& RecTracks)
  {
    float imp = TrueCollision.impactParameter();
    float evPhi = TrueCollision.eventPlaneAngle() / 2.0;
    static constexpr std::array<float, 10> CentEdges = {
      0.0f, 3.49f, 4.93f, 6.98f, 8.55f,
      9.87f, 11.0f, 12.1f, 13.1f, 14.0f};

    static constexpr std::array<float, 9> CentValues = {
      2.5f, 7.5f, 15.0f, 25.0f, 35.0f,
      45.0f, 55.0f, 65.0f, 75.0f};
    float centclass = -999.f;

    for (size_t i = 0; i < CentValues.size(); ++i) {
      if (imp >= CentEdges[i] && imp < CentEdges[i + 1]) {
        centclass = CentValues[i];
        break;
      }
    }
    histos.fill(HIST("hImpactParameter"), imp);
    histos.fill(HIST("hEventPlaneAngle"), evPhi);
    static constexpr float MinCentrality = 0.0f;
    static constexpr float MaxCentrality = 80.0f;

    if (centclass < MinCentrality || centclass > MaxCentrality) {
      return;
    }
    for (const auto& RecCollision : RecCollisions) {
      auto psiFT0C = TrueCollision.eventPlaneAngle();
      /*
  if (!RecCollision.sel8()) {
        continue;
      }
      if (!RecCollision.selection_bit(aod::evsel::kNoITSROFrameBorder)) {
      continue;
      }
      */
      if (std::abs(RecCollision.posZ()) > cfgCutVertex) {
        continue;
      }
      auto oldindex = -999;
      auto rectrackspart = RecTracks.sliceBy(perCollision, RecCollision.globalIndex());
      // loop over reconstructed particle
      for (const auto& track1 : rectrackspart) {
        if (!track1.has_mcParticle()) {
          continue;
        }

        const auto mctrack1 = track1.mcParticle();

        if (selectionTrack(track1) && strategySelectionPID(track1, 0, strategyPID) && std::abs(mctrack1.pdgCode()) == PDG_t::kKPlus && mctrack1.isPhysicalPrimary()) {
          histos.fill(HIST("hSparseKstarMCRecKaonWeight"), centclass, getPhiInRange(mctrack1.phi() - psiFT0C), std::pow(std::cos(2.0 * getPhiInRange(mctrack1.phi() - psiFT0C)), 2.0), mctrack1.pt(), mctrack1.eta());
        }
        if (selectionTrack(track1) && track1.pt() > TPCOnlyPt && track1.hasTOF() && std::abs(track1.tofNSigmaKa()) > nsigmaCutTOF && std::abs(track1.tpcNSigmaKa()) < nsigmaCutTPC && std::abs(mctrack1.pdgCode()) == PDG_t::kKPlus && mctrack1.isPhysicalPrimary()) {
          histos.fill(HIST("hSparseKstarMCRecKaonMissMatchWeight"), centclass, getPhiInRange(mctrack1.phi() - psiFT0C), std::pow(std::cos(2.0 * getPhiInRange(mctrack1.phi() - psiFT0C)), 2.0), mctrack1.pt(), mctrack1.eta());
        }
        if (selectionTrack(track1) && strategySelectionPID(track1, 1, strategyPID) && std::abs(mctrack1.pdgCode()) == PDG_t::kPiPlus && mctrack1.isPhysicalPrimary()) {
          histos.fill(HIST("hSparseKstarMCRecPionWeight"), centclass, getPhiInRange(mctrack1.phi() - psiFT0C), std::pow(std::cos(2.0 * getPhiInRange(mctrack1.phi() - psiFT0C)), 2.0), mctrack1.pt(), mctrack1.eta());
        }
        if (selectionTrack(track1) && track1.pt() > TPCOnlyPt && track1.hasTOF() && std::abs(track1.tofNSigmaPi()) > nsigmaCutTOF && std::abs(track1.tpcNSigmaPi()) < nsigmaCutTPC && std::abs(mctrack1.pdgCode()) == PDG_t::kPiPlus && mctrack1.isPhysicalPrimary()) {
          histos.fill(HIST("hSparseKstarMCRecPionMissMatchWeight"), centclass, getPhiInRange(mctrack1.phi() - psiFT0C), std::pow(std::cos(2.0 * getPhiInRange(mctrack1.phi() - psiFT0C)), 2.0), mctrack1.pt(), mctrack1.eta());
        }
        auto track1ID = track1.index();
        for (const auto& track2 : rectrackspart) {
          if (!track2.has_mcParticle()) {
            continue;
          }
          auto track2ID = track2.index();
          if (track2ID <= track1ID) {
            continue;
          }
          const auto mctrack2 = track2.mcParticle();
          int track1PDG = std::abs(mctrack1.pdgCode());
          int track2PDG = std::abs(mctrack2.pdgCode());
          if (!mctrack1.isPhysicalPrimary()) {
            continue;
          }
          if (!mctrack2.isPhysicalPrimary()) {
            continue;
          }
          if (track1PDG != PDG_t::kKPlus || track2PDG != PDG_t::kPiPlus) {
            continue;
          }
          if (!selectionTrack(track1) || !selectionTrack(track2) || track1.sign() * track2.sign() > 0) {
            continue;
          }
          // PID check
          if (!isTOFOnly && (!strategySelectionPID(track1, 0, strategyPID) || !strategySelectionPID(track2, 1, strategyPID))) {
            continue;
          }
          if (isTOFOnly && (!selectionPID2(track1, 0) || !selectionPID2(track2, 1))) {
            continue;
          }
          for (const auto& mothertrack1 : mctrack1.mothers_as<aod::McParticles>()) {
            for (const auto& mothertrack2 : mctrack2.mothers_as<aod::McParticles>()) {
              if (mothertrack1.pdgCode() != mothertrack2.pdgCode()) {
                continue;
              }
              if (mothertrack1 != mothertrack2) {
                continue;
              }
              if (std::abs(mothertrack1.y()) > confRapidity) {
                continue;
              }
              if (std::abs(mothertrack1.pdgCode()) != o2::constants::physics::kK0Star892) {
                continue;
              }
              // if (avoidsplitrackMC && oldindex == mothertrack1.globalIndex()) {
              if (avoidsplitrackMC && oldindex == mothertrack1.index()) {
                histos.fill(HIST("h1PhiRecsplit"), mothertrack1.pt());
                continue;
              }
              // oldindex = mothertrack1.globalIndex();
              oldindex = mothertrack1.index();
              auto phiMinusPsi = getPhiInRange(mothertrack1.phi() - psiFT0C);
              histos.fill(HIST("hSparseKstarMCRecWeight"), centclass, phiMinusPsi, std::pow(std::cos(2.0 * phiMinusPsi), 2.0), mothertrack1.pt(), mothertrack1.eta());
            }
          }
        }
      }
      // loop over generated particle
      for (const auto& mcParticle : GenParticles) {
        static constexpr float MaxEtaAcceptance = 0.8f;
        if (std::abs(mcParticle.eta()) > MaxEtaAcceptance) {
          continue;
        }
        if (std::abs(mcParticle.pdgCode()) == PDG_t::kKPlus && mcParticle.isPhysicalPrimary()) {
          histos.fill(HIST("hSparseKstarMCGenKaonWeight"), centclass, getPhiInRange(mcParticle.phi() - psiFT0C), std::pow(std::cos(2.0 * getPhiInRange(mcParticle.phi() - psiFT0C)), 2.0), mcParticle.pt(), mcParticle.eta());
        }
        if (std::abs(mcParticle.pdgCode()) == PDG_t::kPiPlus && mcParticle.isPhysicalPrimary()) {
          histos.fill(HIST("hSparseKstarMCGenPionWeight"), centclass, getPhiInRange(mcParticle.phi() - psiFT0C), std::pow(std::cos(2.0 * getPhiInRange(mcParticle.phi() - psiFT0C)), 2.0), mcParticle.pt(), mcParticle.eta());
        }
        if (std::abs(mcParticle.y()) > confRapidity) {
          continue;
        }
        if (mcParticle.pdgCode() != o2::constants::physics::kK0Star892) {
          continue;
        }
        auto kDaughters = mcParticle.daughters_as<aod::McParticles>();
        static constexpr std::size_t NumberOfDaughters = 2;

        if (kDaughters.size() != NumberOfDaughters) {
          continue;
        }
        auto daughtp = false;
        auto daughtm = false;
        for (const auto& kCurrentDaughter : kDaughters) {
          if (!kCurrentDaughter.isPhysicalPrimary()) {
            continue;
          }
          if (kCurrentDaughter.pdgCode() == +PDG_t::kKPlus) {
            if (kCurrentDaughter.pt() > cfgCutPT && std::abs(kCurrentDaughter.eta()) < cfgCutEta) {
              daughtp = true;
            }
            kaonPlus = ROOT::Math::PxPyPzMVector(kCurrentDaughter.px(), kCurrentDaughter.py(), kCurrentDaughter.pz(), massKa);
          } else if (kCurrentDaughter.pdgCode() == -PDG_t::kPiPlus) {
            if (kCurrentDaughter.pt() > cfgCutPT && std::abs(kCurrentDaughter.eta()) < cfgCutEta) {
              daughtm = true;
            }
            pionMinus = ROOT::Math::PxPyPzMVector(kCurrentDaughter.px(), kCurrentDaughter.py(), kCurrentDaughter.pz(), massPi);
          }
        }
        if (daughtp && daughtm) {
          auto phiMinusPsiGen = getPhiInRange(mcParticle.phi() - psiFT0C);
          histos.fill(HIST("hSparseKstarMCGenWeight"), centclass, phiMinusPsiGen, std::pow(std::cos(2.0 * phiMinusPsiGen), 2.0), mcParticle.pt(), mcParticle.eta());
        }
      }
    } // rec collision loop

  } // process MC
  PROCESS_SWITCH(Kstarpbpb, processMCkstarWeight, "Process MC kstar Weight", false);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<Kstarpbpb>(cfgc)};
}
