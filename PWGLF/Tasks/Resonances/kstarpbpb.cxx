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
// sourav.kundu@cern.ch , sarjeeta.gami@cern.ch

#include "PWGLF/DataModel/EPCalibrationTables.h"
#include "PWGMM/Mult/DataModel/Index.h" // for Particles2Tracks table

#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseITS.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "CCDB/CcdbApi.h"
#include "CommonConstants/PhysicsConstants.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/StepTHn.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

#include "Math/GenVector/Boost.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "TF1.h"
#include "TRandom3.h"
#include <TDirectory.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <THn.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TPDGCode.h>

#include <array>
#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;
using namespace o2::aod::rctsel;
struct kstarpbpb {

  struct : ConfigurableGroup {
    Configurable<std::string> cfgURL{"cfgURL", "http://alice-ccdb.cern.ch", "Address of the CCDB to browse"};
    Configurable<int64_t> nolaterthan{"ccdb-no-later-than", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "Latest acceptable timestamp of creation for the object"};
  } cfgCcdbParam;

  // Enable access to the CCDB for the offset and correction constants and save them in dedicated variables.
  Service<o2::ccdb::BasicCCDBManager> ccdb;
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
  Configurable<bool> additionalEvSel2{"additionalEvSel2", true, "Additional evsel2"};
  Configurable<bool> additionalEvSel3{"additionalEvSel3", true, "Additional evsel3"};
  Configurable<float> cfgCutPT{"cfgCutPT", 0.2, "PT cut on daughter track"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8, "Eta cut on daughter track"};
  Configurable<float> cfgCutDCAxy{"cfgCutDCAxy", 2.0f, "DCAxy range for tracks"};
  Configurable<float> cfgCutDCAz{"cfgCutDCAz", 2.0f, "DCAz range for tracks"};
  Configurable<bool> useGlobalTrack{"useGlobalTrack", true, "use Global track"};
  Configurable<bool> usepolar{"usepolar", true, "flag to fill type of SA"};
  Configurable<float> nsigmaCutTOF{"nsigmacutTOF", 3.0, "Value of the TOF Nsigma cut"};
  Configurable<float> nsigmaCutTPC{"nsigmacutTPC", 3.0, "Value of the TPC Nsigma cut"};
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
  Configurable<bool> removefaketrak{"removefaketrack", true, "Remove fake track from momentum difference"};
  Configurable<float> ConfFakeKaonCut{"ConfFakeKaonCut", 0.1, "Cut based on track from momentum difference"};
  ConfigurableAxis configThnAxisV2{"configThnAxisV2", {400, -16, 16}, "V2"};
  Configurable<bool> additionalEvsel{"additionalEvsel", false, "Additional event selcection"};
  Configurable<bool> timFrameEvsel{"timFrameEvsel", false, "TPC Time frame boundary cut"};
  Configurable<bool> additionalEvselITS{"additionalEvselITS", true, "Additional event selcection for ITS"};
  Configurable<bool> ispTdepPID{"ispTdepPID", true, "pT dependent PID"};
  Configurable<bool> isNoTOF{"isNoTOF", true, "isNoTOF"};
  Configurable<bool> PDGcheck{"PDGcheck", true, "PDGcheck"};
  Configurable<int> strategyPID{"strategyPID", 2, "PID strategy"};
  Configurable<float> cfgCutTOFBeta{"cfgCutTOFBeta", 0.0, "cut TOF beta"};
  Configurable<bool> additionalQAplots{"additionalQAplots", true, "Additional QA plots"};
  Configurable<bool> additionalQAplots1{"additionalQAplots1", true, "Additional QA plots"};
  Configurable<float> confMinRot{"confMinRot", 5.0 * TMath::Pi() / 6.0, "Minimum of rotation"};
  Configurable<float> confMaxRot{"confMaxRot", 7.0 * TMath::Pi() / 6.0, "Maximum of rotation"};
  Configurable<int> nBkgRotations{"nBkgRotations", 9, "Number of rotated copies (background) per each original candidate"};
  Configurable<bool> fillRotation{"fillRotation", true, "fill rotation"};
  Configurable<bool> same{"same", true, "same event"};
  Configurable<bool> like{"like", false, "like-sign"};
  Configurable<bool> fillSA{"fillSA", true, "same event SA"};
  Configurable<bool> fillOccupancy{"fillOccupancy", false, "fill Occupancy"};
  Configurable<int> cfgOccupancyCut{"cfgOccupancyCut", 500, "Occupancy cut"};
  Configurable<bool> useWeight{"useWeight", false, "use EP dep effi weight"};
  Configurable<bool> useSP{"useSP", false, "use SP"};
  Configurable<bool> genacceptancecut{"genacceptancecut", true, "use acceptance cut for generated"};
  Configurable<bool> avoidsplitrackMC{"avoidsplitrackMC", false, "avoid split track in MC"};
  Configurable<std::string> ConfWeightPath{"ConfWeightPath", "Users/s/skundu/My/Object/fitweight", "Path to gain calibration"};
  ConfigurableAxis axisPtKaonWeight{"axisPtKaonWeight", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f}, "pt axis"};

  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter centralityFilter = nabs(aod::cent::centFT0C) < cfgCutCentrality;
  Filter acceptanceFilter = (nabs(aod::track::eta) < cfgCutEta && nabs(aod::track::pt) > cfgCutPT);
  Filter DCAcutFilter = (nabs(aod::track::dcaXY) < cfgCutDCAxy) && (nabs(aod::track::dcaZ) < cfgCutDCAz);

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
      if (same) {
        histos.add("hSparseV2SASameEvent_V2", "hSparseV2SASameEvent_V2", HistType::kTHnSparseF, {configThnAxisInvMass, configThnAxisPt, configThnAxisV2, configThnAxisCentrality});
      }
      if (like) {
        histos.add("hSparseV2SAlikeEventNN_V2", "hSparseV2SAlikeEventNN_V2", HistType::kTHnSparseF, {configThnAxisInvMass, configThnAxisPt, configThnAxisV2, configThnAxisCentrality});
        histos.add("hSparseV2SAlikeEventPP_V2", "hSparseV2SAlikeEventPP_V2", HistType::kTHnSparseF, {configThnAxisInvMass, configThnAxisPt, configThnAxisV2, configThnAxisCentrality});
      }
    }
    if (fillRotation) {
      if (!fillSA) {
        histos.add("hRotation", "hRotation", kTH1F, {{360, 0.0, 2.0 * TMath::Pi()}});
        histos.add("hSparseV2SASameEventRotational_V2", "hSparseV2SASameEventRotational_V2", HistType::kTHnSparseF, {configThnAxisInvMass, configThnAxisPt, configThnAxisV2, configThnAxisCentrality});
      }
    }

    if (fillSA) {
      histos.add("hSparseSAvsrapsameunlike", "hSparseSAvsrapsameunlike", HistType::kTHnSparseF, {configThnAxisInvMass, configThnAxisPt, configThnAxisV2, configrapAxis, configThnAxisCentrality}, true);
      histos.add("hSparseSAvsrapsamelike", "hSparseSAvsrapsamelike", HistType::kTHnSparseF, {configThnAxisInvMass, configThnAxisPt, configThnAxisV2, configrapAxis, configThnAxisCentrality}, true);
      histos.add("hSparseSAvsraprot", "hSparseSAvsraprot", HistType::kTHnSparseF, {configThnAxisInvMass, configThnAxisPt, configThnAxisV2, configrapAxis, configThnAxisCentrality}, true);
      histos.add("hSparseSAvsrapmix", "hSparseSAvsrapmix", HistType::kTHnSparseF, {configThnAxisInvMass, configThnAxisPt, configThnAxisV2, configrapAxis, configThnAxisCentrality}, true);
    }
    if (!fillSA) {
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
      histos.add("hEventPlaneAngle", "hEventPlaneAngle", kTH1F, {{200, -2.0f * TMath::Pi(), 2.0f * TMath::Pi()}});
      histos.add("hSparseKstarMCGenWeight", "hSparseKstarMCGenWeight", HistType::kTHnSparseD, {configThnAxisCentrality, {36, 0.0f, TMath::Pi()}, {400, 0.0f, 1}, configThnAxisPt, {8, -0.8, 0.8}});
      histos.add("hSparseKstarMCRecWeight", "hSparseKstarMCRecWeight", HistType::kTHnSparseD, {configThnAxisCentrality, {36, 0.0f, TMath::Pi()}, {400, 0.0f, 1}, configThnAxisPt, {8, -0.8, 0.8}});
      histos.add("hSparseKstarMCGenKaonWeight", "hSparseKstarMCGenKaonWeight", HistType::kTHnSparseD, {configThnAxisCentrality, {36, 0.0f, TMath::Pi()}, {400, 0.0f, 1}, axisPtKaonWeight, {8, -0.8, 0.8}});
      histos.add("hSparseKstarMCRecKaonWeight", "hSparseKstarMCRecKaonWeight", HistType::kTHnSparseD, {configThnAxisCentrality, {36, 0.0f, TMath::Pi()}, {400, 0.0f, 1}, axisPtKaonWeight, {8, -0.8, 0.8}});
      histos.add("hSparseKstarMCRecKaonMissMatchWeight", "hSparseKstarMCRecKaonMissMatchWeight", HistType::kTHnSparseD, {configThnAxisCentrality, {36, 0.0f, TMath::Pi()}, {400, 0.0f, 1}, axisPtKaonWeight, {8, -0.8, 0.8}});
      histos.add("hSparseKstarMCGenPionWeight", "hSparseKstarMCGenPionWeight", HistType::kTHnSparseD, {configThnAxisCentrality, {36, 0.0f, TMath::Pi()}, {400, 0.0f, 1}, axisPtKaonWeight, {8, -0.8, 0.8}});
      histos.add("hSparseKstarMCRecPionWeight", "hSparseKstarMCRecPionWeight", HistType::kTHnSparseD, {configThnAxisCentrality, {36, 0.0f, TMath::Pi()}, {400, 0.0f, 1}, axisPtKaonWeight, {8, -0.8, 0.8}});
      histos.add("hSparseKstarMCRecPionMissMatchWeight", "hSparseKstarMCRecPionMissMatchWeight", HistType::kTHnSparseD, {configThnAxisCentrality, {36, 0.0f, TMath::Pi()}, {400, 0.0f, 1}, axisPtKaonWeight, {8, -0.8, 0.8}});
    }
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

    // Event selection cut additional - Alex
    if (additionalEvsel) {
      fMultPVCutLow = new TF1("fMultPVCutLow", "[0]+[1]*x+[2]*x*x+[3]*x*x*x - 2.5*([4]+[5]*x+[6]*x*x+[7]*x*x*x+[8]*x*x*x*x)", 0, 100);
      fMultPVCutLow->SetParameters(2834.66, -87.0127, 0.915126, -0.00330136, 332.513, -12.3476, 0.251663, -0.00272819, 1.12242e-05);
      fMultPVCutHigh = new TF1("fMultPVCutHigh", "[0]+[1]*x+[2]*x*x+[3]*x*x*x + 2.5*([4]+[5]*x+[6]*x*x+[7]*x*x*x+[8]*x*x*x*x)", 0, 100);
      fMultPVCutHigh->SetParameters(2834.66, -87.0127, 0.915126, -0.00330136, 332.513, -12.3476, 0.251663, -0.00272819, 1.12242e-05);
      fMultCutLow = new TF1("fMultCutLow", "[0]+[1]*x+[2]*x*x+[3]*x*x*x - 2.5*([4]+[5]*x)", 0, 100);
      fMultCutLow->SetParameters(1893.94, -53.86, 0.502913, -0.0015122, 109.625, -1.19253);
      fMultCutHigh = new TF1("fMultCutHigh", "[0]+[1]*x+[2]*x*x+[3]*x*x*x + 3.*([4]+[5]*x)", 0, 100);
      fMultCutHigh->SetParameters(1893.94, -53.86, 0.502913, -0.0015122, 109.625, -1.19253);
      fMultMultPVCut = new TF1("fMultMultPVCut", "[0]+[1]*x+[2]*x*x", 0, 5000);
      fMultMultPVCut->SetParameters(-0.1, 0.785, -4.7e-05);
    }
    ccdb->setURL(cfgCcdbParam.cfgURL);
    ccdbApi.init("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setCreatedNotAfter(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
  }

  double massKa = o2::constants::physics::MassKPlus;
  double massPi = o2::constants::physics::MassPiMinus;

  template <typename TCollision>
  bool eventSelected(TCollision collision, const float& centrality)
  {
    if (collision.alias_bit(kTVXinTRD)) {
      // TRD triggered
      // return 0;
    }
    auto multNTracksPV = collision.multNTracksPV();
    if (multNTracksPV < fMultPVCutLow->Eval(centrality))
      return 0;
    if (multNTracksPV > fMultPVCutHigh->Eval(centrality))
      return 0;
    // if (multTrk < fMultCutLow->Eval(centrality))
    //  return 0;
    // if (multTrk > fMultCutHigh->Eval(centrality))
    //  return 0;
    // if (multTrk > fMultMultPVCut->Eval(multNTracksPV))
    //  return 0;

    return 1;
  }
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

  template <typename T>
  bool selectionPIDNew(const T& candidate, int PID)
  {
    if (PID == 0) {
      if (candidate.pt() < 0.5 && TMath::Abs(candidate.tpcNSigmaKa()) < nsigmaCutTPC) {
        return true;
      }
      if (candidate.pt() >= 0.5 && TMath::Abs(candidate.tpcNSigmaKa()) < nsigmaCutTPC && candidate.hasTOF() && TMath::Abs(candidate.tofNSigmaKa()) < nsigmaCutTOF && candidate.beta() > cfgCutTOFBeta) {
        return true;
      }
      if (candidate.pt() >= 0.5 && TMath::Abs(candidate.tpcNSigmaKa()) < nsigmaCutTPC && !candidate.hasTOF()) {
        return true;
      }
    } else if (PID == 1) {
      if (candidate.pt() < 0.5 && TMath::Abs(candidate.tpcNSigmaPi()) < nsigmaCutTPC) {
        return true;
      }
      if (candidate.pt() >= 0.5 && TMath::Abs(candidate.tpcNSigmaPi()) < nsigmaCutTPC && candidate.hasTOF() && TMath::Abs(candidate.tofNSigmaPi()) < nsigmaCutTOF && candidate.beta() > cfgCutTOFBeta) {
        return true;
      }
      if (candidate.pt() >= 0.5 && TMath::Abs(candidate.tpcNSigmaPi()) < nsigmaCutTPC && !candidate.hasTOF()) {
        return true;
      }
    }
    return false;
  }
  template <typename T>
  bool selectionPID2(const T& candidate, int PID)
  {
    if (PID == 0) {
      if (candidate.hasTOF() && candidate.beta() > cfgCutTOFBeta && TMath::Abs(candidate.tofNSigmaKa()) < nsigmaCutTOF) {
        return true;
      }
    }
    if (PID == 1) {
      if (candidate.hasTOF() && candidate.beta() > cfgCutTOFBeta && TMath::Abs(candidate.tofNSigmaPi()) < nsigmaCutTOF) {
        return true;
      }
    }
    return false;
  }
  template <typename T>
  bool selectionPID(const T& candidate, int PID)
  {
    if (PID == 0) {
      if (!isNoTOF && !candidate.hasTOF() && TMath::Abs(candidate.tpcNSigmaKa()) < nsigmaCutTPC) {
        return true;
      }
      if (!isNoTOF && candidate.hasTOF() && ((candidate.tofNSigmaKa() * candidate.tofNSigmaKa()) + (candidate.tpcNSigmaKa() * candidate.tpcNSigmaKa())) < (nsigmaCutCombined * nsigmaCutCombined)) {
        return true;
      }
      if (isNoTOF && std::abs(candidate.tpcNSigmaKa()) < nsigmaCutTPC) {
        return true;
      }
    } else if (PID == 1) {
      if (!isNoTOF && !candidate.hasTOF() && TMath::Abs(candidate.tpcNSigmaPi()) < nsigmaCutTPC) {
        return true;
      }
      if (!isNoTOF && candidate.hasTOF() && ((candidate.tofNSigmaPi() * candidate.tofNSigmaPi()) + (candidate.tpcNSigmaPi() * candidate.tpcNSigmaPi())) < (nsigmaCutCombined * nsigmaCutCombined)) {
        return true;
      }
      if (isNoTOF && std::abs(candidate.tpcNSigmaPi()) < nsigmaCutTPC) {
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
        if (!isNoTOF && !candidate.hasTOF() && TMath::Abs(candidate.tpcNSigmaKa()) < nsigmaCutTPC) {
          return true;
        }
        if (!isNoTOF && candidate.hasTOF() && TMath::Abs(candidate.tpcNSigmaKa()) < nsigmaCutTPC && TMath::Abs(candidate.tofNSigmaKa()) < nsigmaCutTOF) {
          return true;
        }
        if (isNoTOF && std::abs(candidate.tpcNSigmaKa()) < nsigmaCutTPC) {
          return true;
        }
      } else if (strategy == 1) {
        if (candidate.pt() < 0.5 && TMath::Abs(candidate.tpcNSigmaKa()) < nsigmaCutTPC) {
          return true;
        }
        if (candidate.pt() >= 0.5 && TMath::Abs(candidate.tpcNSigmaKa()) < nsigmaCutTPC && candidate.hasTOF() && TMath::Abs(candidate.tofNSigmaKa()) < nsigmaCutTOF && candidate.beta() > cfgCutTOFBeta) {
          return true;
        }
        if (!useGlobalTrack && !candidate.hasTPC()) {
          return true;
        }
      } else if (strategy == 2) {
        if (candidate.pt() < 0.5 && TMath::Abs(candidate.tpcNSigmaKa()) < nsigmaCutTPC) {
          return true;
        }
        if (candidate.pt() >= 0.5 && TMath::Abs(candidate.tpcNSigmaKa()) < nsigmaCutTPC && candidate.hasTOF() && TMath::Abs(candidate.tofNSigmaKa()) < nsigmaCutTOF && candidate.beta() > cfgCutTOFBeta) {
          return true;
        }
        if (candidate.pt() >= 0.5 && TMath::Abs(candidate.tpcNSigmaKa()) < nsigmaCutTPC && !candidate.hasTOF()) {
          return true;
        }
      }
    }
    if (PID == 1) {
      if (strategy == 0) {
        if (!isNoTOF && !candidate.hasTOF() && TMath::Abs(candidate.tpcNSigmaPi()) < nsigmaCutTPC) {
          return true;
        }
        if (!isNoTOF && candidate.hasTOF() && TMath::Abs(candidate.tpcNSigmaPi()) < nsigmaCutTPC && TMath::Abs(candidate.tofNSigmaPi()) < nsigmaCutTOF) {
          return true;
        }
        if (isNoTOF && std::abs(candidate.tpcNSigmaPi()) < nsigmaCutTPC) {
          return true;
        }
      } else if (strategy == 1) {
        if (candidate.pt() < 0.5 && TMath::Abs(candidate.tpcNSigmaPi()) < nsigmaCutTPC) {
          return true;
        }
        if (candidate.pt() >= 0.5 && TMath::Abs(candidate.tpcNSigmaPi()) < nsigmaCutTPC && candidate.hasTOF() && TMath::Abs(candidate.tofNSigmaPi()) < nsigmaCutTOF && candidate.beta() > cfgCutTOFBeta) {
          return true;
        }
        if (!useGlobalTrack && !candidate.hasTPC()) {
          return true;
        }
      } else if (strategy == 2) {
        if (candidate.pt() < 0.5 && TMath::Abs(candidate.tpcNSigmaPi()) < nsigmaCutTPC) {
          return true;
        }
        if (candidate.pt() >= 0.5 && TMath::Abs(candidate.tpcNSigmaPi()) < nsigmaCutTPC && candidate.hasTOF() && TMath::Abs(candidate.tofNSigmaPi()) < nsigmaCutTOF && candidate.beta() > cfgCutTOFBeta) {
          return true;
        }
        if (candidate.pt() >= 0.5 && TMath::Abs(candidate.tpcNSigmaPi()) < nsigmaCutTPC && !candidate.hasTOF()) {
          return true;
        }
      }
    }
    return false;
  }

  double GetPhiInRange(double phi)
  {
    double result = phi;
    while (result < 0) {
      result = result + 2. * TMath::Pi() / 2;
    }
    while (result > 2. * TMath::Pi() / 2) {
      result = result - 2. * TMath::Pi() / 2;
    }
    return result;
  }
  template <typename T>
  bool isFakeKaon(T const& track, int /*PID*/)
  {
    const auto pglobal = track.p();
    const auto ptpc = track.tpcInnerParam();
    if (TMath::Abs(pglobal - ptpc) > ConfFakeKaonCut) {
      return true;
    }
    return false;
  }
  ConfigurableAxis axisVertex{"axisVertex", {20, -10, 10}, "vertex axis for bin"};
  ConfigurableAxis axisMultiplicityClass{"axisMultiplicityClass", {20, 0, 100}, "multiplicity percentile for bin"};
  ConfigurableAxis axisEPAngle{"axisEPAngle", {6, -TMath::Pi() / 2, TMath::Pi() / 2}, "event plane angle"};
  ConfigurableAxis axisOccup{"axisOccup", {20, -0.5, 40000.0}, "occupancy axis"};
  double v2, v2Rot;

  using BinningTypeVertexContributor = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0C, aod::epcalibrationtable::PsiFT0C>;
  ROOT::Math::PxPyPzMVector KstarMother, fourVecDauCM, daughter1, daughter2, kaonrot, kstarrot, KaonPlus, PionMinus;
  ROOT::Math::XYZVector threeVecDauCM, threeVecDauCMXY, eventplaneVec, eventplaneVecNorm;
  ROOT::Math::PxPyPzMVector daughter2rot, fourVecDauCMrot;
  ROOT::Math::XYZVector threeVecDauCMrot, threeVecDauCMXYrot;

  int currentRunNumber = -999;
  int lastRunNumber = -999;
  TH2D* hweight;
  void processSE(EventCandidates::iterator const& collision, TrackCandidates const& tracks, aod::BCsWithTimestamps const&)
  {
    histos.fill(HIST("hEvtSelInfo"), 0.5);
    if (rctCut.requireRCTFlagChecker && !rctCut.rctChecker(collision)) {
      return;
    }
    histos.fill(HIST("hEvtSelInfo"), 1.5);
    if (!collision.sel8() || !collision.triggereventep() || !collision.selection_bit(aod::evsel::kNoTimeFrameBorder) || !collision.selection_bit(aod::evsel::kNoITSROFrameBorder) || !collision.selection_bit(aod::evsel::kNoSameBunchPileup) || !collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV) || !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
      return;
    }
    histos.fill(HIST("hEvtSelInfo"), 2.5);
    auto centrality = collision.centFT0C();
    auto multTPC = collision.multNTracksPV();
    int occupancy = collision.trackOccupancyInTimeRange();
    auto psiFT0C = collision.psiFT0C();
    auto psiFT0A = collision.psiFT0A();
    auto psiTPC = collision.psiTPC();
    auto QFT0C = collision.qFT0C();
    auto QFT0A = collision.qFT0A();
    auto QTPC = collision.qTPC();
    if (fillOccupancy && occupancy > cfgOccupancyCut) {
      return;
    }
    histos.fill(HIST("hEvtSelInfo"), 3.5);
    if (additionalEvsel && !eventSelected(collision, centrality)) {
      return;
    }
    histos.fill(HIST("hEvtSelInfo"), 4.5);
    if (additionalEvselITS && !collision.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll)) {
      return;
    }
    histos.fill(HIST("hEvtSelInfo"), 5.5);
    if (additionalQAplots1) {
      histos.fill(HIST("hFTOCvsTPCSelected"), centrality, multTPC);
      histos.fill(HIST("hPsiFT0C"), centrality, psiFT0C);
      histos.fill(HIST("hPsiFT0A"), centrality, psiFT0A);
      histos.fill(HIST("hPsiTPC"), centrality, psiTPC);
      histos.fill(HIST("ResFT0CTPC"), centrality, TMath::Cos(2.0 * (psiFT0C - psiTPC)));
      histos.fill(HIST("ResFT0CFT0A"), centrality, TMath::Cos(2.0 * (psiFT0C - psiFT0A)));
      histos.fill(HIST("ResFT0ATPC"), centrality, TMath::Cos(2.0 * (psiTPC - psiFT0A)));
      histos.fill(HIST("ResFT0CTPCSP"), centrality, QFT0C * QTPC * TMath::Cos(2.0 * (psiFT0C - psiTPC)));
      histos.fill(HIST("ResFT0CFT0ASP"), centrality, QFT0C * QFT0A * TMath::Cos(2.0 * (psiFT0C - psiFT0A)));
      histos.fill(HIST("ResFT0ATPCSP"), centrality, QTPC * QFT0A * TMath::Cos(2.0 * (psiTPC - psiFT0A)));
      histos.fill(HIST("hCentrality"), centrality);
      histos.fill(HIST("hOccupancy"), occupancy);
      histos.fill(HIST("hVtxZ"), collision.posZ());
    }
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    currentRunNumber = collision.bc_as<aod::BCsWithTimestamps>().runNumber();
    if (useWeight && (currentRunNumber != lastRunNumber)) {
      hweight = ccdb->getForTimeStamp<TH2D>(ConfWeightPath.value, bc.timestamp());
    }
    lastRunNumber = currentRunNumber;
    float weight1 = 1.0;
    float weight2 = 1.0;
    for (auto track1 : tracks) {
      if (!selectionTrack(track1)) {
        continue;
      }
      bool track1kaon = false;
      auto track1ID = track1.globalIndex();
      if (!isTOFOnly && !strategySelectionPID(track1, 0, strategyPID)) {
        continue;
      }
      if (isTOFOnly && !selectionPID2(track1, 0)) {
        continue;
      }
      track1kaon = true;

      if (useWeight) {
        if (track1.pt() < 10.0 && track1.pt() > 0.15) {
          weight1 = 1 + hweight->GetBinContent(hweight->FindBin(centrality, track1.pt() + 0.000005)) * TMath::Cos(2.0 * GetPhiInRange(track1.phi() - psiFT0C));
        } else {
          weight1 = 1;
        }
      }
      for (auto track2 : tracks) {
        if (!selectionTrack(track2)) {
          continue;
        }
        bool track2pion = false;
        auto track2ID = track2.globalIndex();
        if (!isTOFOnly && !strategySelectionPID(track2, 1, strategyPID)) {
          continue;
        }
        if (isTOFOnly && !selectionPID2(track2, 1)) {
          continue;
        }
        track2pion = true;
        if (track2ID == track1ID) {
          continue;
        }
        if (!track1kaon || !track2pion) {
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
          if (track2.pt() < 10.0 && track2.pt() > 0.15) {
            weight2 = 1 + hweight->GetBinContent(hweight->FindBin(centrality, track2.pt() + 0.000005)) * TMath::Cos(2.0 * GetPhiInRange(track2.phi() - psiFT0C));
          } else {
            weight2 = 1;
          }
        }
        daughter1 = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), massKa);
        daughter2 = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), massPi);
        KstarMother = daughter1 + daughter2;
        if (TMath::Abs(KstarMother.Rapidity()) > confRapidity) {
          continue;
        }
        auto phiminuspsi = GetPhiInRange(KstarMother.Phi() - psiFT0C);

        if (useSP) {
          v2 = TMath::Cos(2.0 * phiminuspsi) * QFT0C;
        }
        if (!useSP) {
          v2 = TMath::Cos(2.0 * phiminuspsi);
        }
        auto totalweight = weight1 * weight2;
        if (totalweight <= 0.0000005) {
          totalweight = 1.0;
        }
        if (additionalQAplots1) {
          histos.fill(HIST("ResTrackSPFT0CTPC"), centrality, occupancy, QFT0C * QTPC * TMath::Cos(2.0 * (psiFT0C - psiTPC)));
          histos.fill(HIST("ResTrackSPFT0CFT0A"), centrality, occupancy, QFT0C * QFT0A * TMath::Cos(2.0 * (psiFT0C - psiFT0A)));
          histos.fill(HIST("ResTrackSPFT0ATPC"), centrality, occupancy, QTPC * QFT0A * TMath::Cos(2.0 * (psiTPC - psiFT0A)));
        }
        if (!fillSA) {
          if (same) {
            if (useWeight) {
              histos.fill(HIST("hSparseV2SASameEvent_V2"), KstarMother.M(), KstarMother.Pt(), v2, centrality, 1 / totalweight);
            } else {
              histos.fill(HIST("hSparseV2SASameEvent_V2"), KstarMother.M(), KstarMother.Pt(), v2, centrality);
            }
          }
        }
        int track1Sign = track1.sign();
        int track2Sign = track2.sign();

        if (fillSA) {
          ROOT::Math::Boost boost{KstarMother.BoostToCM()};
          fourVecDauCM = boost(daughter1);
          threeVecDauCM = fourVecDauCM.Vect();
          threeVecDauCMXY = ROOT::Math::XYZVector(threeVecDauCM.X(), threeVecDauCM.Y(), 0.);
          eventplaneVec = ROOT::Math::XYZVector(std::cos(2.0 * psiFT0C), std::sin(2.0 * psiFT0C), 0);
          eventplaneVecNorm = ROOT::Math::XYZVector(std::sin(2.0 * psiFT0C), -std::cos(2.0 * psiFT0C), 0);
          auto cosPhistarminuspsi = GetPhiInRange(fourVecDauCM.Phi() - psiFT0C);
          auto SA = TMath::Cos(2.0 * cosPhistarminuspsi);
          auto cosThetaStar = eventplaneVecNorm.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2()) / std::sqrt(eventplaneVecNorm.Mag2());

          if (track1Sign * track2Sign < 0) {
            if (usepolar) {
              histos.fill(HIST("hSparseSAvsrapsameunlike"), KstarMother.M(), KstarMother.Pt(), cosThetaStar, KstarMother.Rapidity(), centrality);
            } else {
              histos.fill(HIST("hSparseSAvsrapsameunlike"), KstarMother.M(), KstarMother.Pt(), SA, KstarMother.Rapidity(), centrality);
            }
          } else if (track1Sign * track2Sign > 0) {
            if (usepolar) {
              histos.fill(HIST("hSparseSAvsrapsamelike"), KstarMother.M(), KstarMother.Pt(), cosThetaStar, KstarMother.Rapidity(), centrality);
            } else {
              histos.fill(HIST("hSparseSAvsrapsamelike"), KstarMother.M(), KstarMother.Pt(), SA, KstarMother.Rapidity(), centrality);
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
            if (TMath::Abs(kstarrot.Rapidity()) > confRapidity) {
              continue;
            }
            auto phiminuspsiRot = GetPhiInRange(kstarrot.Phi() - psiFT0C);

            if (useSP) {
              v2Rot = TMath::Cos(2.0 * phiminuspsiRot) * QFT0C;
            }
            if (!useSP) {
              v2Rot = TMath::Cos(2.0 * phiminuspsiRot);
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
                auto cosPhistarminuspsirot = GetPhiInRange(fourVecDauCMrot.Phi() - psiFT0C);
                auto SArot = TMath::Cos(2.0 * cosPhistarminuspsirot);
                auto cosThetaStarrot = eventplaneVecNorm.Dot(threeVecDauCMrot) / std::sqrt(threeVecDauCMrot.Mag2()) / std::sqrt(eventplaneVecNorm.Mag2());
                if (usepolar) {
                  histos.fill(HIST("hSparseSAvsraprot"), kstarrot.M(), kstarrot.Pt(), cosThetaStarrot, kstarrot.Rapidity(), centrality);
                } else {
                  histos.fill(HIST("hSparseSAvsraprot"), kstarrot.M(), kstarrot.Pt(), SArot, kstarrot.Rapidity(), centrality);
                }
              }
            }
          }
        }
      }
    }
  }
  PROCESS_SWITCH(kstarpbpb, processSE, "Process Same event latest", true);
  /*
  void processSameEvent(EventCandidates::iterator const& collision, TrackCandidates const& , aod::BCsWithTimestamps const&)
  {
    if (!collision.sel8()) {
      return;
    }
    auto centrality = collision.centFT0C();
    auto multTPC = collision.multNTracksPV();
    auto QFT0C = collision.qFT0C();
    if (!collision.triggereventep()) {
      return;
    }
    if (timFrameEvsel && (!collision.selection_bit(aod::evsel::kNoTimeFrameBorder) || !collision.selection_bit(aod::evsel::kNoITSROFrameBorder))) {
      return;
    }
    if (additionalEvSel2 && (!collision.selection_bit(aod::evsel::kNoSameBunchPileup) || !collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV))) {
      return;
    }
    if (additionalEvSel3 && (!collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard))) {
      return;
    }
    if (additionalEvselITS && !collision.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll)) {
      return;
    }
    int occupancy = collision.trackOccupancyInTimeRange();
    auto posThisColl = posTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto negThisColl = negTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto psiFT0C = collision.psiFT0C();
    auto psiFT0A = collision.psiFT0A();
    auto psiTPC = collision.psiTPC();
    if (fillOccupancy && occupancy >= cfgOccupancyCut) // occupancy info is available for this collision (*)
    {
      return;
    }
    if (additionalEvsel && !eventSelected(collision, centrality)) {
      return;
    }
    if (additionalQAplots1) {
      histos.fill(HIST("hFTOCvsTPCSelected"), centrality, multTPC);
      histos.fill(HIST("hPsiFT0C"), centrality, psiFT0C);
      histos.fill(HIST("hPsiFT0A"), centrality, psiFT0A);
      histos.fill(HIST("hPsiTPC"), centrality, psiTPC);
      histos.fill(HIST("ResFT0CTPC"), centrality, TMath::Cos(2.0 * (psiFT0C - psiTPC)));
      histos.fill(HIST("ResFT0CFT0A"), centrality, TMath::Cos(2.0 * (psiFT0C - psiFT0A)));
      histos.fill(HIST("ResFT0ATPC"), centrality, TMath::Cos(2.0 * (psiTPC - psiFT0A)));
      histos.fill(HIST("hCentrality"), centrality);
      histos.fill(HIST("hOccupancy"), occupancy);
      histos.fill(HIST("hVtxZ"), collision.posZ());
    }
    for (auto track1 : posThisColl) {
      if (!selectionTrack(track1)) {
        continue;
      }
      if (additionalQAplots) {
        histos.fill(HIST("QAbefore/TPC_Nsigma_allka"), track1.pt(), track1.tpcNSigmaKa(), centrality);
        histos.fill(HIST("QAbefore/TOF_Nsigma_allka"), track1.pt(), track1.tofNSigmaKa(), centrality);
        histos.fill(HIST("QAbefore/trkDCAxyka"), track1.dcaXY());
        histos.fill(HIST("QAbefore/trkDCAzka"), track1.dcaZ());
        histos.fill(HIST("QAbefore/TOF_TPC_Mapka_allka"), track1.tofNSigmaKa(), track1.tpcNSigmaKa());
      }

      bool track1pion = false;
      bool track1kaon = false;
      if (ispTdepPID && !isTOFOnly && !(selectionPIDNew(track1, 0) || selectionPIDNew(track1, 1))) {
        continue;
      }
      if (!ispTdepPID && !isTOFOnly && !(selectionPID(track1, 0) || selectionPID(track1, 1))) {
        continue;
      }
      if (isTOFOnly && !(selectionPID2(track1, 0) || selectionPID2(track1, 1))) {
        continue;
      }
      auto track1ID = track1.globalIndex();
      for (auto track2 : negThisColl) {
        bool track2pion = false;
        bool track2kaon = false;
        if (!selectionTrack(track2)) {
          continue;
        }
        if (additionalQAplots) {
          histos.fill(HIST("QAbefore/TOF_TPC_Mapka_allpi"), track2.tofNSigmaPi(), track2.tpcNSigmaPi());
          histos.fill(HIST("QAbefore/TPC_Nsigma_allpi"), track2.pt(), track2.tpcNSigmaPi(), centrality);
          histos.fill(HIST("QAbefore/TOF_Nsigma_allpi"), track2.pt(), track2.tofNSigmaPi(), centrality);
          histos.fill(HIST("QAbefore/trkDCAxypi"), track2.dcaXY());
          histos.fill(HIST("QAbefore/trkDCAzpi"), track2.dcaZ());
        }
        if (ispTdepPID && !isTOFOnly && !(selectionPIDNew(track2, 0) || selectionPIDNew(track2, 1))) {
          continue;
        }
        if (!ispTdepPID && !isTOFOnly && !(selectionPID(track2, 0) || selectionPID(track2, 1))) {
          continue;
        }
        if (isTOFOnly && !(selectionPID2(track2, 0) || selectionPID2(track2, 1))) {
          continue;
        }
        auto track2ID = track2.globalIndex();
        if (track2ID == track1ID) {
          continue;
        }
        if (track1.sign() * track2.sign() > 0) {
          continue;
        }

        if (ispTdepPID && !isTOFOnly) {
          if (selectionPIDNew(track1, 1) && selectionPIDNew(track2, 0)) {
            track1pion = true;
            track2kaon = true;
            if (removefaketrak && isFakeKaon(track2, 0)) {
              continue;
            }
          }
          if (selectionPIDNew(track2, 1) && selectionPIDNew(track1, 0)) {
            track2pion = true;
            track1kaon = true;
            if (removefaketrak && isFakeKaon(track1, 0)) {
              continue;
            }
          }
        }
        if (!ispTdepPID && !isTOFOnly) {
          if (selectionPID(track1, 1) && selectionPID(track2, 0)) {
            track1pion = true;
            track2kaon = true;
            if (removefaketrak && isFakeKaon(track2, 0)) {
              continue;
            }
          }
          if (selectionPID(track2, 1) && selectionPID(track1, 0)) {
            track2pion = true;
            track1kaon = true;
            if (removefaketrak && isFakeKaon(track1, 0)) {
              continue;
            }
          }
        }
        if (isTOFOnly) {
          if (selectionPID2(track1, 1) && selectionPID2(track2, 0)) {
            track1pion = true;
            track2kaon = true;
            if (removefaketrak && isFakeKaon(track2, 0)) {
              continue;
            }
          }
          if (selectionPID2(track2, 1) && selectionPID2(track1, 0)) {
            track2pion = true;
            track1kaon = true;
            if (removefaketrak && isFakeKaon(track1, 0)) {
              continue;
            }
          }
        }
        if (same) {
          if (track1kaon && track2pion) {
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
            daughter1 = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), massKa);
            daughter2 = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), massPi);
          } else if (track1pion && track2kaon) {
            daughter1 = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), massPi);
            daughter2 = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), massKa);
          } else {
            continue;
          }

          KstarMother = daughter1 + daughter2;
          if (TMath::Abs(KstarMother.Rapidity()) > confRapidity) {
            continue;
          }
          auto phiminuspsi = GetPhiInRange(KstarMother.Phi() - psiFT0C);

          if (useSP) {
            v2 = TMath::Cos(2.0 * phiminuspsi) * QFT0C;
          }
          if (!useSP) {
            v2 = TMath::Cos(2.0 * phiminuspsi);
          }

          histos.fill(HIST("hSparseV2SASameEvent_V2"), KstarMother.M(), KstarMother.Pt(), v2, centrality);
        }
        if (fillRotation) {
          for (int nrotbkg = 0; nrotbkg < nBkgRotations; nrotbkg++) {
            auto anglestart = confMinRot;
            auto angleend = confMaxRot;
            auto anglestep = (angleend - anglestart) / (1.0 * (nBkgRotations - 1));
            auto rotangle = anglestart + nrotbkg * anglestep;
            histos.fill(HIST("hRotation"), rotangle);
            if (track1kaon && track2pion) {
              auto rotkaonPx = track1.px() * std::cos(rotangle) - track1.py() * std::sin(rotangle);
              auto rotkaonPy = track1.px() * std::sin(rotangle) + track1.py() * std::cos(rotangle);
              kaonrot = ROOT::Math::PxPyPzMVector(rotkaonPx, rotkaonPy, track1.pz(), massKa);
              daughter2 = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), massPi);
            } else if (track1pion && track2kaon) {
              auto rotkaonPx = track2.px() * std::cos(rotangle) - track2.py() * std::sin(rotangle);
              auto rotkaonPy = track2.px() * std::sin(rotangle) + track2.py() * std::cos(rotangle);
              kaonrot = ROOT::Math::PxPyPzMVector(rotkaonPx, rotkaonPy, track2.pz(), massKa);
              daughter2 = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), massPi);
            } else {
              continue;
            }
            kstarrot = kaonrot + daughter2;
            if (TMath::Abs(kstarrot.Rapidity()) > confRapidity) {
              continue;
            }
            auto phiminuspsiRot = GetPhiInRange(kstarrot.Phi() - psiFT0C);

            if (useSP) {
              v2Rot = TMath::Cos(2.0 * phiminuspsiRot) * QFT0C;
            }
            if (!useSP) {
              v2Rot = TMath::Cos(2.0 * phiminuspsiRot);
            }

            histos.fill(HIST("hSparseV2SASameEventRotational_V2"), kstarrot.M(), kstarrot.Pt(), v2Rot, centrality);
          }
        }
      }
    }
  }
  PROCESS_SWITCH(kstarpbpb, processSameEvent, "Process Same event", false);

  void processlikeEvent(EventCandidates::iterator const& collision, TrackCandidates const& tracks, aod::BCs const&)
  {
    if (!collision.sel8()) {
      return;
    }
    auto centrality = collision.centFT0C();
    auto QFT0C = collision.qFT0C();
    if (!collision.triggereventep()) {
      return;
    }
    if (timFrameEvsel && (!collision.selection_bit(aod::evsel::kNoTimeFrameBorder) || !collision.selection_bit(aod::evsel::kNoITSROFrameBorder))) {
      return;
    }
    if (additionalEvSel2 && (!collision.selection_bit(aod::evsel::kNoSameBunchPileup) || !collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV))) {
      return;
    }
    if (additionalEvSel3 && (!collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard))) {
      return;
    }
    if (additionalEvselITS && !collision.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll)) {
      return;
    }
    int occupancy = collision.trackOccupancyInTimeRange();
    auto psiFT0C = collision.psiFT0C();
    if (fillOccupancy && occupancy >= cfgOccupancyCut) // occupancy info is available for this collision (*)
    {
      return;
    }

    if (additionalEvsel && !eventSelected(collision, centrality)) {
      return;
    }
    for (auto track1 : tracks) {
      if (!selectionTrack(track1)) {
        continue;
      }
      bool track1pion = false;
      bool track1kaon = false;
      if (ispTdepPID && !isTOFOnly && !(selectionPIDNew(track1, 0) || selectionPIDNew(track1, 1))) {
        continue;
      }
      if (!ispTdepPID && !isTOFOnly && !(selectionPID(track1, 0) || selectionPID(track1, 1))) {
        continue;
      }
      if (isTOFOnly && !(selectionPID2(track1, 0) || selectionPID2(track1, 1))) {
        continue;
      }
      for (auto track2 : tracks) {
        bool track2pion = false;
        bool track2kaon = false;
        if (!selectionTrack(track2)) {
          continue;
        }
        if (ispTdepPID && !isTOFOnly && !(selectionPIDNew(track2, 0) || selectionPIDNew(track2, 1))) {
          continue;
        }
        if (!ispTdepPID && !isTOFOnly && !(selectionPID(track2, 0) || selectionPID(track2, 1))) {
          continue;
        }
        if (isTOFOnly && !(selectionPID2(track2, 0) || selectionPID2(track2, 1))) {
          continue;
        }
        if (track1.sign() * track2.sign() < 0) {
          continue;
        }

        if (ispTdepPID && !isTOFOnly) {
          if (selectionPIDNew(track1, 1) && selectionPIDNew(track2, 0)) {
            track1pion = true;
            track2kaon = true;
            if (removefaketrak && isFakeKaon(track2, 0)) {
              continue;
            }
          }
          if (selectionPIDNew(track2, 1) && selectionPIDNew(track1, 0)) {
            track2pion = true;
            track1kaon = true;
            if (removefaketrak && isFakeKaon(track1, 0)) {
              continue;
            }
          }
        }
        if (!ispTdepPID && !isTOFOnly) {
          if (selectionPID(track1, 1) && selectionPID(track2, 0)) {
            track1pion = true;
            track2kaon = true;
            if (removefaketrak && isFakeKaon(track2, 0)) {
              continue;
            }
          }
          if (selectionPID(track2, 1) && selectionPID(track1, 0)) {
            track2pion = true;
            track1kaon = true;
            if (removefaketrak && isFakeKaon(track1, 0)) {
              continue;
            }
          }
        }
        if (isTOFOnly) {
          if (selectionPID2(track1, 1) && selectionPID2(track2, 0)) {
            track1pion = true;
            track2kaon = true;
            if (removefaketrak && isFakeKaon(track2, 0)) {
              continue;
            }
          }
          if (selectionPID2(track2, 1) && selectionPID2(track1, 0)) {
            track2pion = true;
            track1kaon = true;
            if (removefaketrak && isFakeKaon(track1, 0)) {
              continue;
            }
          }
        }
        if (track1kaon && track2pion) {
          if (track1.sign() < 0 && track2.sign() < 0) {

            daughter1 = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), massKa);
            daughter2 = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), massPi);

            KstarMother = daughter1 + daughter2;
            if (TMath::Abs(KstarMother.Rapidity()) > confRapidity) {
              continue;
            }
            auto phiminuspsi = GetPhiInRange(KstarMother.Phi() - psiFT0C);

            if (useSP) {
              v2 = TMath::Cos(2.0 * phiminuspsi) * QFT0C;
            }
            if (!useSP) {
              v2 = TMath::Cos(2.0 * phiminuspsi);
            }
            histos.fill(HIST("hSparseV2SAlikeEventNN_V2"), KstarMother.M(), KstarMother.Pt(), v2, centrality);
          }
        } else if (track1pion && track2kaon) {
          if (track1.sign() > 0 && track2.sign() > 0) {
            daughter1 = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), massPi);
            daughter2 = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), massKa);
            KstarMother = daughter1 + daughter2;
            if (TMath::Abs(KstarMother.Rapidity()) > confRapidity) {
              continue;
            }
            auto phiminuspsi = GetPhiInRange(KstarMother.Phi() - psiFT0C);

            if (useSP) {
              v2 = TMath::Cos(2.0 * phiminuspsi) * QFT0C;
            }
            if (!useSP) {
              v2 = TMath::Cos(2.0 * phiminuspsi);
            }

            histos.fill(HIST("hSparseV2SAlikeEventPP_V2"), KstarMother.M(), KstarMother.Pt(), v2, centrality);
          }
        }
      }
    }
  }

  PROCESS_SWITCH(kstarpbpb, processlikeEvent, "Process like event", false);
  */

  void processMixedEvent(EventCandidates const& collisions, TrackCandidates const& tracks)
  {

    auto tracksTuple = std::make_tuple(tracks);
    BinningTypeVertexContributor binningOnPositions{{axisVertex, axisMultiplicityClass, axisOccup}, true};
    SameKindPair<EventCandidates, TrackCandidates, BinningTypeVertexContributor> pair{binningOnPositions, cfgNoMixedEvents, -1, collisions, tracksTuple, &cache};
    for (auto& [collision1, tracks1, collision2, tracks2] : pair) {
      if (rctCut.requireRCTFlagChecker && !rctCut.rctChecker(collision1)) {
        continue;
      }
      if (rctCut.requireRCTFlagChecker && !rctCut.rctChecker(collision2)) {
        continue;
      }
      if (!collision1.sel8() || !collision1.triggereventep() || !collision1.selection_bit(aod::evsel::kNoTimeFrameBorder) || !collision1.selection_bit(aod::evsel::kNoITSROFrameBorder) || !collision1.selection_bit(aod::evsel::kNoSameBunchPileup) || !collision1.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV) || !collision1.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
        continue;
      }
      if (!collision2.sel8() || !collision2.triggereventep() || !collision2.selection_bit(aod::evsel::kNoTimeFrameBorder) || !collision2.selection_bit(aod::evsel::kNoITSROFrameBorder) || !collision2.selection_bit(aod::evsel::kNoSameBunchPileup) || !collision2.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV) || !collision2.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
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
      auto centrality2 = collision2.centFT0C();
      auto psiFT0C = collision1.psiFT0C();
      auto QFT0C = collision1.qFT0C();
      if (additionalEvsel && !eventSelected(collision1, centrality)) {
        continue;
      }
      if (additionalEvsel && !eventSelected(collision2, centrality2)) {
        continue;
      }
      if (additionalEvselITS && !collision1.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll)) {
        continue;
      }
      if (additionalEvselITS && !collision2.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll)) {
        continue;
      }

      for (auto& [track1, track2] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {
        if (track1.sign() * track2.sign() > 0) {
          continue;
        }
        if (!selectionTrack(track1) || !selectionTrack(track2)) {

          continue;
        }
        if (ispTdepPID && !isTOFOnly && !(selectionPIDNew(track1, 0))) {
          continue;
        }
        if (ispTdepPID && !isTOFOnly && !(selectionPIDNew(track2, 1))) {
          continue;
        }
        if (!ispTdepPID && !isTOFOnly && !(selectionPID(track1, 0))) {
          continue;
        }
        if (!ispTdepPID && !isTOFOnly && !(selectionPID(track2, 1))) {
          continue;
        }
        if (isTOFOnly && !selectionPID2(track1, 0)) {
          continue;
        }
        if (isTOFOnly && !selectionPID2(track2, 1)) {
          continue;
        }
        // if (track1.sign() > 0 && track2.sign() < 0) {
        daughter1 = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), massKa);
        daughter2 = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), massPi);
        /*} else if (track1.sign() < 0 && track2.sign() > 0) {
              daughter2 = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), massKa);
              daughter1 = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), massPi);
        }*/
        KstarMother = daughter1 + daughter2;
        if (TMath::Abs(KstarMother.Rapidity()) > confRapidity) {
          continue;
        }
        auto phiminuspsi = GetPhiInRange(KstarMother.Phi() - psiFT0C);

        v2 = TMath::Cos(2.0 * phiminuspsi) * QFT0C;
        if (!fillSA) {
          if (track1.sign() * track2.sign() < 0)
            histos.fill(HIST("hSparseV2SAMixedEvent_V2"), KstarMother.M(), KstarMother.Pt(), v2, centrality);
        }
        if (fillSA) {
          ROOT::Math::Boost boost{KstarMother.BoostToCM()};
          fourVecDauCM = boost(daughter1);
          threeVecDauCM = fourVecDauCM.Vect();
          threeVecDauCMXY = ROOT::Math::XYZVector(threeVecDauCM.X(), threeVecDauCM.Y(), 0.);
          eventplaneVec = ROOT::Math::XYZVector(std::cos(2.0 * psiFT0C), std::sin(2.0 * psiFT0C), 0);
          eventplaneVecNorm = ROOT::Math::XYZVector(std::sin(2.0 * psiFT0C), -std::cos(2.0 * psiFT0C), 0);
          auto cosPhistarminuspsi = GetPhiInRange(fourVecDauCM.Phi() - psiFT0C);
          auto SA = TMath::Cos(2.0 * cosPhistarminuspsi);
          auto cosThetaStar = eventplaneVecNorm.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2()) / std::sqrt(eventplaneVecNorm.Mag2());
          if (usepolar) {
            histos.fill(HIST("hSparseSAvsrapmix"), KstarMother.M(), KstarMother.Pt(), cosThetaStar, KstarMother.Rapidity(), centrality);
          } else {
            histos.fill(HIST("hSparseSAvsrapmix"), KstarMother.M(), KstarMother.Pt(), SA, KstarMother.Rapidity(), centrality);
          }
        }
      }
    }
  }
  PROCESS_SWITCH(kstarpbpb, processMixedEvent, "Process Mixed event", true);
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
    for (auto& RecCollision : RecCollisions) {
      auto psiFT0C = 0.0;
      histos.fill(HIST("hMC"), 3);
      if (!RecCollision.sel8()) {
        histos.fill(HIST("hMC"), 4);
        continue;
      }

      if (!RecCollision.selection_bit(aod::evsel::kNoTimeFrameBorder) || !RecCollision.selection_bit(aod::evsel::kNoITSROFrameBorder)) {
        histos.fill(HIST("hMC"), 5);
        continue;
      }
      if (TMath::Abs(RecCollision.posZ()) > cfgCutVertex) {
        histos.fill(HIST("hMC"), 6);
        continue;
      }
      histos.fill(HIST("hMC"), 7);
      auto centrality = RecCollision.centFT0C();
      histos.fill(HIST("CentPercentileMCRecHist"), centrality);
      auto oldindex = -999;
      auto Rectrackspart = RecTracks.sliceBy(perCollision, RecCollision.globalIndex());
      // loop over reconstructed particle
      for (auto track1 : Rectrackspart) {
        if (!selectionTrack(track1)) {
          continue;
        }
        if (ispTdepPID && !(selectionPIDNew(track1, 0))) {
          continue;
        }
        if (!ispTdepPID && !(selectionPID(track1, 0))) {
          continue;
        }
        if (!track1.has_mcParticle()) {
          continue;
        }
        auto track1ID = track1.index();
        for (auto track2 : Rectrackspart) {
          auto track2ID = track2.index();
          if (track2ID <= track1ID) {
            continue;
          }
          if (!selectionTrack(track2)) {
            continue;
          }
          if (ispTdepPID && !(selectionPIDNew(track2, 1))) {
            continue;
          }
          if (!ispTdepPID && !(selectionPID(track2, 1))) {
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
          int track1PDG = TMath::Abs(mctrack1.pdgCode());
          int track2PDG = TMath::Abs(mctrack2.pdgCode());
          if (!mctrack1.isPhysicalPrimary()) {
            continue;
          }
          if (!mctrack2.isPhysicalPrimary()) {
            continue;
          }
          if (!(track1PDG == 321 && track2PDG == 211)) {
            continue;
          }
          for (auto& mothertrack1 : mctrack1.mothers_as<aod::McParticles>()) {
            for (auto& mothertrack2 : mctrack2.mothers_as<aod::McParticles>()) {
              if (mothertrack1.pdgCode() != mothertrack2.pdgCode()) {
                continue;
              }
              if (mothertrack1 != mothertrack2) {
                continue;
              }
              if (TMath::Abs(mothertrack1.y()) > confRapidity) {
                continue;
              }
              if (PDGcheck && TMath::Abs(mothertrack1.pdgCode()) != 313) {
                continue;
              }
              if (ispTdepPID && !(selectionPIDNew(track1, 0) || selectionPIDNew(track2, 1))) {
                continue;
              }
              if (!ispTdepPID && !(selectionPID(track1, 0) || selectionPID(track2, 1))) {
                continue;
              }
              if (avoidsplitrackMC && oldindex == mothertrack1.globalIndex()) {
                continue;
              }
              oldindex = mothertrack1.globalIndex();
              if (track1.sign() > 0 && track2.sign() < 0) {
                KaonPlus = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), massKa);
                PionMinus = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), massPi);
              }
              if (track1.sign() < 0 && track2.sign() > 0) {
                PionMinus = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), massKa);
                KaonPlus = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), massPi);
              }
              KstarMother = KaonPlus + PionMinus;
              if (TMath::Abs(KstarMother.Rapidity()) > confRapidity) {
                continue;
              }
              auto phiminuspsi = GetPhiInRange(KstarMother.Phi() - psiFT0C);

              v2 = TMath::Cos(2.0 * phiminuspsi);

              histos.fill(HIST("hSparseV2SARec_V2"), KstarMother.M(), KstarMother.Pt(), v2, centrality);
              histos.fill(HIST("h2PhiRec2"), KstarMother.pt(), centrality);
              histos.fill(HIST("hpt"), KstarMother.Pt());
            }
          }
        }
      }
      // loop over generated particle
      for (auto& mcParticle : GenParticles) {
        if (TMath::Abs(mcParticle.y()) > confRapidity) {
          continue;
        }
        if (PDGcheck && mcParticle.pdgCode() != 313) {
          continue;
        }
        auto kDaughters = mcParticle.daughters_as<aod::McParticles>();
        if (kDaughters.size() != 2) {
          continue;
        }
        auto daughtp = false;
        auto daughtm = false;
        for (auto kCurrentDaughter : kDaughters) {
          if (!kCurrentDaughter.isPhysicalPrimary()) {
            continue;
          }
          if (kCurrentDaughter.pdgCode() == +321) {
            if (genacceptancecut && kCurrentDaughter.pt() > cfgCutPT && TMath::Abs(kCurrentDaughter.eta()) < cfgCutEta) {
              daughtp = true;
            }
            if (!genacceptancecut) {
              daughtp = true;
            }
            KaonPlus = ROOT::Math::PxPyPzMVector(kCurrentDaughter.px(), kCurrentDaughter.py(), kCurrentDaughter.pz(), massKa);
          } else if (kCurrentDaughter.pdgCode() == -211) {
            if (genacceptancecut && kCurrentDaughter.pt() > cfgCutPT && TMath::Abs(kCurrentDaughter.eta()) < cfgCutEta) {
              daughtm = true;
            }
            if (!genacceptancecut) {
              daughtm = true;
            }
            PionMinus = ROOT::Math::PxPyPzMVector(kCurrentDaughter.px(), kCurrentDaughter.py(), kCurrentDaughter.pz(), massPi);
          }
        }
        if (daughtp && daughtm) {
          KstarMother = KaonPlus + PionMinus;
          if (TMath::Abs(KstarMother.Rapidity()) > confRapidity) {
            continue;
          }
          auto phiminuspsi = GetPhiInRange(KstarMother.Phi() - psiFT0C);

          v2 = TMath::Cos(2.0 * phiminuspsi);

          histos.fill(HIST("hSparseV2SAGen_V2"), KstarMother.M(), KstarMother.Pt(), v2, centrality);
          histos.fill(HIST("h2PhiGen2"), KstarMother.pt(), centrality);
        }
      }
    } // rec collision loop

  } // process MC
  PROCESS_SWITCH(kstarpbpb, processMC, "Process MC", false);

  void processMCkstarWeight(CollisionMCTrueTable::iterator const& TrueCollision, CollisionMCRecTableCentFT0C const& RecCollisions, TrackMCTrueTable const& GenParticles, FilTrackMCRecTable const& RecTracks)
  {
    float imp = TrueCollision.impactParameter();
    float evPhi = TrueCollision.eventPlaneAngle() / 2.0;
    float centclass = -999;
    if (imp >= 0 && imp < 3.49) {
      centclass = 2.5;
    }
    if (imp >= 3.49 && imp < 4.93) {
      centclass = 7.5;
    }
    if (imp >= 4.93 && imp < 6.98) {
      centclass = 15.0;
    }
    if (imp >= 6.98 && imp < 8.55) {
      centclass = 25.0;
    }
    if (imp >= 8.55 && imp < 9.87) {
      centclass = 35.0;
    }
    if (imp >= 9.87 && imp < 11) {
      centclass = 45.0;
    }
    if (imp >= 11 && imp < 12.1) {
      centclass = 55.0;
    }
    if (imp >= 12.1 && imp < 13.1) {
      centclass = 65.0;
    }
    if (imp >= 13.1 && imp < 14) {
      centclass = 75.0;
    }
    histos.fill(HIST("hImpactParameter"), imp);
    histos.fill(HIST("hEventPlaneAngle"), evPhi);
    if (centclass < 0.0 || centclass > 80.0) {
      return;
    }
    for (auto& RecCollision : RecCollisions) {
      auto psiFT0C = TrueCollision.eventPlaneAngle();
      /*
  if (!RecCollision.sel8()) {
        continue;
      }
      if (!RecCollision.selection_bit(aod::evsel::kNoITSROFrameBorder)) {
      continue;
      }
      */
      if (TMath::Abs(RecCollision.posZ()) > cfgCutVertex) {
        continue;
      }
      auto oldindex = -999;
      auto Rectrackspart = RecTracks.sliceBy(perCollision, RecCollision.globalIndex());
      // loop over reconstructed particle
      for (auto track1 : Rectrackspart) {
        if (!track1.has_mcParticle()) {
          continue;
        }

        const auto mctrack1 = track1.mcParticle();

        if (selectionTrack(track1) && strategySelectionPID(track1, 0, strategyPID) && TMath::Abs(mctrack1.pdgCode()) == 321 && mctrack1.isPhysicalPrimary()) {
          histos.fill(HIST("hSparseKstarMCRecKaonWeight"), centclass, GetPhiInRange(mctrack1.phi() - psiFT0C), TMath::Power(TMath::Cos(2.0 * GetPhiInRange(mctrack1.phi() - psiFT0C)), 2.0), mctrack1.pt(), mctrack1.eta());
        }
        if (selectionTrack(track1) && track1.pt() > 0.5 && track1.hasTOF() && TMath::Abs(track1.tofNSigmaKa()) > nsigmaCutTOF && TMath::Abs(track1.tpcNSigmaKa()) < nsigmaCutTPC && TMath::Abs(mctrack1.pdgCode()) == 321 && mctrack1.isPhysicalPrimary()) {
          histos.fill(HIST("hSparseKstarMCRecKaonMissMatchWeight"), centclass, GetPhiInRange(mctrack1.phi() - psiFT0C), TMath::Power(TMath::Cos(2.0 * GetPhiInRange(mctrack1.phi() - psiFT0C)), 2.0), mctrack1.pt(), mctrack1.eta());
        }
        if (selectionTrack(track1) && strategySelectionPID(track1, 1, strategyPID) && TMath::Abs(mctrack1.pdgCode()) == 211 && mctrack1.isPhysicalPrimary()) {
          histos.fill(HIST("hSparseKstarMCRecPionWeight"), centclass, GetPhiInRange(mctrack1.phi() - psiFT0C), TMath::Power(TMath::Cos(2.0 * GetPhiInRange(mctrack1.phi() - psiFT0C)), 2.0), mctrack1.pt(), mctrack1.eta());
        }
        if (selectionTrack(track1) && track1.pt() > 0.5 && track1.hasTOF() && TMath::Abs(track1.tofNSigmaPi()) > nsigmaCutTOF && TMath::Abs(track1.tpcNSigmaPi()) < nsigmaCutTPC && TMath::Abs(mctrack1.pdgCode()) == 211 && mctrack1.isPhysicalPrimary()) {
          histos.fill(HIST("hSparseKstarMCRecPionMissMatchWeight"), centclass, GetPhiInRange(mctrack1.phi() - psiFT0C), TMath::Power(TMath::Cos(2.0 * GetPhiInRange(mctrack1.phi() - psiFT0C)), 2.0), mctrack1.pt(), mctrack1.eta());
        }
        auto track1ID = track1.index();
        for (auto track2 : Rectrackspart) {
          if (!track2.has_mcParticle()) {
            continue;
          }
          auto track2ID = track2.index();
          if (track2ID <= track1ID) {
            continue;
          }
          const auto mctrack2 = track2.mcParticle();
          int track1PDG = TMath::Abs(mctrack1.pdgCode());
          int track2PDG = TMath::Abs(mctrack2.pdgCode());
          if (!mctrack1.isPhysicalPrimary()) {
            continue;
          }
          if (!mctrack2.isPhysicalPrimary()) {
            continue;
          }
          if (!(track1PDG == 321 && track2PDG == 211)) {
            continue;
          }
          if (!selectionTrack(track1) || !selectionTrack(track2) || track1.sign() * track2.sign() > 0) {
            continue;
          }
          // PID check
          if (ispTdepPID && !isTOFOnly && (!strategySelectionPID(track1, 0, strategyPID) || !strategySelectionPID(track2, 1, strategyPID))) {
            continue;
          }
          if (!ispTdepPID && !isTOFOnly && (!selectionPID(track1, 0) || !selectionPID(track2, 1))) {
            continue;
          }
          if (isTOFOnly && (!selectionPID2(track1, 0) || !selectionPID2(track2, 1))) {
            continue;
          }
          for (auto& mothertrack1 : mctrack1.mothers_as<aod::McParticles>()) {
            for (auto& mothertrack2 : mctrack2.mothers_as<aod::McParticles>()) {
              if (mothertrack1.pdgCode() != mothertrack2.pdgCode()) {
                continue;
              }
              if (mothertrack1 != mothertrack2) {
                continue;
              }
              if (TMath::Abs(mothertrack1.y()) > confRapidity) {
                continue;
              }
              if (TMath::Abs(mothertrack1.pdgCode()) != 313) {
                continue;
              }
              // if (avoidsplitrackMC && oldindex == mothertrack1.globalIndex()) {
              if (avoidsplitrackMC && oldindex == mothertrack1.index()) {
                histos.fill(HIST("h1PhiRecsplit"), mothertrack1.pt());
                continue;
              }
              // oldindex = mothertrack1.globalIndex();
              oldindex = mothertrack1.index();
              auto PhiMinusPsi = GetPhiInRange(mothertrack1.phi() - psiFT0C);
              histos.fill(HIST("hSparseKstarMCRecWeight"), centclass, PhiMinusPsi, TMath::Power(TMath::Cos(2.0 * PhiMinusPsi), 2.0), mothertrack1.pt(), mothertrack1.eta());
            }
          }
        }
      }
      // loop over generated particle
      for (auto& mcParticle : GenParticles) {
        if (TMath::Abs(mcParticle.eta()) > 0.8) // main acceptance
          continue;
        if (TMath::Abs(mcParticle.pdgCode()) == 321 && mcParticle.isPhysicalPrimary()) {
          histos.fill(HIST("hSparseKstarMCGenKaonWeight"), centclass, GetPhiInRange(mcParticle.phi() - psiFT0C), TMath::Power(TMath::Cos(2.0 * GetPhiInRange(mcParticle.phi() - psiFT0C)), 2.0), mcParticle.pt(), mcParticle.eta());
        }
        if (TMath::Abs(mcParticle.pdgCode()) == 211 && mcParticle.isPhysicalPrimary()) {
          histos.fill(HIST("hSparseKstarMCGenPionWeight"), centclass, GetPhiInRange(mcParticle.phi() - psiFT0C), TMath::Power(TMath::Cos(2.0 * GetPhiInRange(mcParticle.phi() - psiFT0C)), 2.0), mcParticle.pt(), mcParticle.eta());
        }
        if (TMath::Abs(mcParticle.y()) > confRapidity) {
          continue;
        }
        if (mcParticle.pdgCode() != 313) {
          continue;
        }
        auto kDaughters = mcParticle.daughters_as<aod::McParticles>();
        if (kDaughters.size() != 2) {
          continue;
        }
        auto daughtp = false;
        auto daughtm = false;
        for (auto kCurrentDaughter : kDaughters) {
          if (!kCurrentDaughter.isPhysicalPrimary()) {
            continue;
          }
          if (kCurrentDaughter.pdgCode() == +321) {
            if (kCurrentDaughter.pt() > cfgCutPT && TMath::Abs(kCurrentDaughter.eta()) < cfgCutEta) {
              daughtp = true;
            }
            KaonPlus = ROOT::Math::PxPyPzMVector(kCurrentDaughter.px(), kCurrentDaughter.py(), kCurrentDaughter.pz(), massKa);
          } else if (kCurrentDaughter.pdgCode() == -211) {
            if (kCurrentDaughter.pt() > cfgCutPT && TMath::Abs(kCurrentDaughter.eta()) < cfgCutEta) {
              daughtm = true;
            }
            PionMinus = ROOT::Math::PxPyPzMVector(kCurrentDaughter.px(), kCurrentDaughter.py(), kCurrentDaughter.pz(), massPi);
          }
        }
        if (daughtp && daughtm) {
          auto PhiMinusPsiGen = GetPhiInRange(mcParticle.phi() - psiFT0C);
          histos.fill(HIST("hSparseKstarMCGenWeight"), centclass, PhiMinusPsiGen, TMath::Power(TMath::Cos(2.0 * PhiMinusPsiGen), 2.0), mcParticle.pt(), mcParticle.eta());
        }
      }
    } // rec collision loop

  } // process MC
  PROCESS_SWITCH(kstarpbpb, processMCkstarWeight, "Process MC kstar Weight", false);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<kstarpbpb>(cfgc, TaskName{"kstarpbpb"})};
}
