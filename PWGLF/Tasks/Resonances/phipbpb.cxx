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
// Phi meson spin alignment task
// sourav.kundu@cern.ch

#include <TH1F.h>
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
#include <cmath>
#include <array>
#include <cstdlib>

#include "TRandom3.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "Math/GenVector/Boost.h"
#include "TF1.h"

#include "PWGLF/DataModel/EPCalibrationTables.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/StepTHn.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/Core/trackUtilities.h"
#include "CommonConstants/PhysicsConstants.h"
#include "Common/Core/TrackSelection.h"
#include "Framework/ASoAHelpers.h"
#include "ReconstructionDataFormats/Track.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "CCDB/BasicCCDBManager.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;
struct phipbpb {

  int mRunNumber;
  int multEstimator;
  float d_bz;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Service<o2::framework::O2DatabasePDG> pdg;

  // CCDB options
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};

  // events
  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};
  Configurable<float> cfgCutCentrality{"cfgCutCentrality", 80.0f, "Accepted maximum Centrality"};
  // track
  Configurable<bool> fillRapidity{"fillRapidity", true, "fill rapidity bin"};
  Configurable<bool> useGlobalTrack{"useGlobalTrack", true, "use Global track"};
  Configurable<float> cfgCutCharge{"cfgCutCharge", 0.0, "cut on Charge"};
  Configurable<float> cfgCutPT{"cfgCutPT", 0.2, "PT cut on daughter track"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8, "Eta cut on daughter track"};
  Configurable<float> cfgCutDCAxy{"cfgCutDCAxy", 2.0f, "DCAxy range for tracks"};
  Configurable<float> cfgCutDCAz{"cfgCutDCAz", 2.0f, "DCAz range for tracks"};
  Configurable<float> nsigmaCutTPC{"nsigmacutTPC", 3.0, "Value of the TPC Nsigma cut"};
  Configurable<float> nsigmaCutCombined{"nsigmaCutCombined", 3.0, "Value of the TOF Nsigma cut"};
  Configurable<int> cfgNoMixedEvents{"cfgNoMixedEvents", 1, "Number of mixed events per event"};
  Configurable<int> cfgITScluster{"cfgITScluster", 0, "Number of ITS cluster"};
  Configurable<int> cfgTPCcluster{"cfgTPCcluster", 70, "Number of TPC cluster"};
  Configurable<bool> isNoTOF{"isNoTOF", false, "isNoTOF"};
  Configurable<bool> isDeepAngle{"isDeepAngle", false, "Deep Angle cut"};
  Configurable<bool> ispTdepPID{"ispTdepPID", true, "pT dependent PID"};
  Configurable<double> cfgDeepAngle{"cfgDeepAngle", 0.04, "Deep Angle cut value"};
  Configurable<double> confRapidity{"confRapidity", 0.5, "Rapidity cut"};
  ConfigurableAxis configThnAxisInvMass{"configThnAxisInvMass", {120, 0.98, 1.1}, "#it{M} (GeV/#it{c}^{2})"};
  ConfigurableAxis configThnAxisPt{"configThnAxisPt", {100, 0.0, 10.}, "#it{p}_{T} (GeV/#it{c})"};
  ConfigurableAxis configThnAxisCosThetaStar{"configThnAxisCosThetaStar", {10, -1.0, 1.}, "cos(#vartheta)"};
  ConfigurableAxis configThnAxisCentrality{"configThnAxisCentrality", {8, 0., 80}, "Centrality"};
  ConfigurableAxis configThnAxisPhiminusPsi{"configThnAxisPhiminusPsi", {6, 0.0, TMath::Pi()}, "#phi - #psi"};
  ConfigurableAxis configThnAxisV2{"configThnAxisV2", {200, -1, 1}, "V2"};
  ConfigurableAxis configThnAxisRapidity{"configThnAxisRapidity", {8, 0, 0.8}, "Rapidity"};
  ConfigurableAxis configThnAxisSA{"configThnAxisSA", {200, -1, 1}, "SA"};
  ConfigurableAxis configThnAxiscosthetaSA{"configThnAxiscosthetaSA", {200, 0, 1}, "costhetaSA"};
  Configurable<bool> additionalEvsel{"additionalEvsel", false, "Additional event selcection"};
  Configurable<bool> timFrameEvsel{"timFrameEvsel", false, "TPC Time frame boundary cut"};
  Configurable<bool> isMC{"isMC", false, "use MC"};
  Configurable<bool> genacceptancecut{"genacceptancecut", true, "use acceptance cut for generated"};
  Configurable<bool> avoidsplitrackMC{"avoidsplitrackMC", false, "avoid split track in MC"};
  Configurable<bool> islike{"islike", false, "use like"};
  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter centralityFilter = nabs(aod::cent::centFT0C) < cfgCutCentrality;
  Filter acceptanceFilter = (nabs(aod::track::eta) < cfgCutEta && nabs(aod::track::pt) > cfgCutPT);
  Filter DCAcutFilter = (nabs(aod::track::dcaXY) < cfgCutDCAxy) && (nabs(aod::track::dcaZ) < cfgCutDCAz);
  Filter PIDcutFilter = nabs(aod::pidtpc::tpcNSigmaKa) < nsigmaCutTPC;

  using EventCandidates = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::FV0Mults, aod::TPCMults, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0Cs, aod::CentFT0As, aod::EPCalibrationTables, aod::Mults>>;
  using TrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullKa, aod::pidTOFFullKa>>;

  using CollisionMCTrueTable = aod::McCollisions;
  using TrackMCTrueTable = aod::McParticles;
  using CollisionMCRecTableCentFT0C = soa::SmallGroups<soa::Join<aod::McCollisionLabels, aod::Collisions, aod::CentFT0Cs, aod::EvSels>>;
  using TrackMCRecTable = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::McTrackLabels, aod::TrackSelection, aod::pidTPCFullKa, aod::pidTOFFullKa>;
  using FilTrackMCRecTable = soa::Filtered<TrackMCRecTable>;

  Preslice<TrackMCRecTable> perCollision = aod::track::collisionId;

  SliceCache cache;
  Partition<TrackCandidates> posTracks = aod::track::signed1Pt > cfgCutCharge;
  Partition<TrackCandidates> negTracks = aod::track::signed1Pt < cfgCutCharge;

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // Event selection cuts - Alex
  TF1* fMultPVCutLow = nullptr;
  TF1* fMultPVCutHigh = nullptr;
  TF1* fMultCutLow = nullptr;
  TF1* fMultCutHigh = nullptr;
  TF1* fMultMultPVCut = nullptr;

  void init(o2::framework::InitContext&)
  {
    const AxisSpec thnAxisInvMass{configThnAxisInvMass, "#it{M} (GeV/#it{c}^{2})"};
    const AxisSpec thnAxisPt{configThnAxisPt, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec thnAxisCosThetaStarOP{configThnAxisCosThetaStar, "cos(#vartheta_{OP})"};
    const AxisSpec thnAxisCosThetaStarIP{configThnAxisCosThetaStar, "cos(#vartheta_{IP})"};
    const AxisSpec thnAxisPhiminusPsi{configThnAxisPhiminusPsi, "#phi - #psi"};
    const AxisSpec thnAxisCentrality{configThnAxisCentrality, "Centrality (%)"};
    const AxisSpec thnAxisV2{configThnAxisV2, "V2"};
    const AxisSpec thnAxisRapidity{configThnAxisRapidity, "Rapidity"};
    const AxisSpec thnAxisSA{configThnAxisSA, "SA"};
    const AxisSpec thnAxiscosthetaSA{configThnAxiscosthetaSA, "costhetaSA"};
    AxisSpec phiAxis = {500, -6.28, 6.28, "phi"};
    AxisSpec resAxis = {400, -2, 2, "Res"};
    AxisSpec centAxis = {8, 0, 80, "V0M (%)"};

    histos.add("hpTvsRapidity", "pT vs Rapidity", kTH2F, {{100, 0.0f, 10.0f}, {300, -1.5f, 1.5f}});
    histos.add("hFTOCvsTPCNoCut", "Mult correlation FT0C vs. TPC without any cut", kTH2F, {{80, 0.0f, 80.0f}, {100, -0.5f, 5999.5f}});
    histos.add("hFTOCvsTPC", "Mult correlation FT0C vs. TPC", kTH2F, {{80, 0.0f, 80.0f}, {100, -0.5f, 5999.5f}});
    histos.add("hFTOCvsTPCSelected", "Mult correlation FT0C vs. TPC after selection", kTH2F, {{80, 0.0f, 80.0f}, {100, -0.5f, 5999.5f}});
    histos.add("hCentrality", "Centrality distribution", kTH1F, {{200, 0.0, 200.0}});
    histos.add("hVtxZ", "Vertex distribution in Z;Z (cm)", kTH1F, {{400, -20.0, 20.0}});
    histos.add("hEta", "Eta distribution", kTH1F, {{200, -1.0f, 1.0f}});
    histos.add("hDcaxy", "Dcaxy distribution", kTH1F, {{200, -1.0f, 1.0f}});
    histos.add("hDcaz", "Dcaz distribution", kTH1F, {{200, -1.0f, 1.0f}});
    histos.add("hNsigmaKaonTPC", "NsigmaKaon TPC distribution", kTH1F, {{200, -10.0f, 10.0f}});
    histos.add("hNsigmaKaonTOF", "NsigmaKaon TOF distribution", kTH1F, {{200, -10.0f, 10.0f}});
    histos.add("hPsiFT0C", "PsiFT0C", kTH2F, {centAxis, phiAxis});
    histos.add("hPsiFT0A", "PsiFT0A", kTH2F, {centAxis, phiAxis});
    histos.add("hPsiTPC", "PsiTPC", kTH2F, {centAxis, phiAxis});
    histos.add("hPsiTPCR", "PsiTPCR", kTH2F, {centAxis, phiAxis});
    histos.add("hPsiTPCL", "PsiTPCL", kTH2F, {centAxis, phiAxis});
    if (!fillRapidity) {
      histos.add("hSparseV2SASameEvent_costhetastarOP", "hSparseV2SASameEvent_costhetastarOP", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisCosThetaStarOP, thnAxisPhiminusPsi, thnAxisCentrality});
      histos.add("hSparseV2SASameEvent_costhetastarIP", "hSparseV2SASameEvent_costhetastarIP", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisCosThetaStarOP, thnAxisPhiminusPsi, thnAxisCentrality});
      histos.add("hSparseV2SASameEvent_SA", "hSparseV2SASameEvent_SA", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisSA, thnAxisPhiminusPsi, thnAxisCentrality});
      histos.add("hSparseV2SASameEvent_costheta_SA", "hSparseV2SASameEvent_costheta_SA", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxiscosthetaSA, thnAxisPhiminusPsi, thnAxisCentrality});
      histos.add("hSparseV2SASameEvent_SA_A0", "hSparseV2SASameEvent_SA_A0", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisSA, thnAxisPhiminusPsi, thnAxisCentrality});
      histos.add("hSparseV2SASameEvent_V2", "hSparseV2SASameEvent_V2", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisV2, thnAxisCentrality});
      histos.add("hSparseV2SAMixedEvent_costhetastarOP", "hSparseV2SAMixedEvent_costhetastarOP", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisCosThetaStarOP, thnAxisPhiminusPsi, thnAxisCentrality});
      histos.add("hSparseV2SAMixedEvent_costhetastarIP", "hSparseV2SAMixedEvent_costhetastarIP", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisCosThetaStarOP, thnAxisPhiminusPsi, thnAxisCentrality});
      histos.add("hSparseV2SAMixedEvent_SA", "hSparseV2SAMixedEvent_SA", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisSA, thnAxisPhiminusPsi, thnAxisCentrality});
      histos.add("hSparseV2SAMixedEvent_costheta_SA", "hSparseV2SAMixedEvent_costheta_SA", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxiscosthetaSA, thnAxisPhiminusPsi, thnAxisCentrality});
      histos.add("hSparseV2SAMixedEvent_SA_A0", "hSparseV2SAMixedEvent_SA_A0", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisSA, thnAxisPhiminusPsi, thnAxisCentrality});
      histos.add("hSparseV2SAMixedEvent_V2", "hSparseV2SAMixedEvent_V2", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisV2, thnAxisCentrality});
    }
    if (fillRapidity) {
      histos.add("hSparseV2SASameEvent_costhetastarOP", "hSparseV2SASameEvent_costhetastarOP", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisCosThetaStarOP, thnAxisRapidity, thnAxisCentrality});
      histos.add("hSparseV2SASameEvent_costhetastarIP", "hSparseV2SASameEvent_costhetastarIP", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisCosThetaStarOP, thnAxisRapidity, thnAxisCentrality});
      histos.add("hSparseV2SASameEvent_SA", "hSparseV2SASameEvent_SA", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisSA, thnAxisRapidity, thnAxisCentrality});
      histos.add("hSparseV2SASameEvent_costheta_SA", "hSparseV2SASameEvent_costheta_SA", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxiscosthetaSA, thnAxisRapidity, thnAxisCentrality});
      histos.add("hSparseV2SASameEvent_SA_A0", "hSparseV2SASameEvent_SA_A0", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisSA, thnAxisRapidity, thnAxisCentrality});
      histos.add("hSparseV2SASameEvent_V2", "hSparseV2SASameEvent_V2", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisV2, thnAxisCentrality});
      histos.add("hSparseV2SAMixedEvent_costhetastarOP", "hSparseV2SAMixedEvent_costhetastarOP", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisCosThetaStarOP, thnAxisRapidity, thnAxisCentrality});
      histos.add("hSparseV2SAMixedEvent_costhetastarIP", "hSparseV2SAMixedEvent_costhetastarIP", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisCosThetaStarOP, thnAxisRapidity, thnAxisCentrality});
      histos.add("hSparseV2SAMixedEvent_SA", "hSparseV2SAMixedEvent_SA", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisSA, thnAxisRapidity, thnAxisCentrality});
      histos.add("hSparseV2SAMixedEvent_costheta_SA", "hSparseV2SAMixedEvent_costheta_SA", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxiscosthetaSA, thnAxisRapidity, thnAxisCentrality});
      histos.add("hSparseV2SAMixedEvent_SA_A0", "hSparseV2SAMixedEvent_SA_A0", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisSA, thnAxisRapidity, thnAxisCentrality});
      histos.add("hSparseV2SAMixedEvent_V2", "hSparseV2SAMixedEvent_V2", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisV2, thnAxisCentrality});
    }
    // histogram for resolution
    histos.add("ResFT0CTPC", "ResFT0CTPC", kTH2F, {centAxis, resAxis});
    histos.add("ResFT0CTPCR", "ResFT0CTPCR", kTH2F, {centAxis, resAxis});
    histos.add("ResFT0CTPCL", "ResFT0CTPCL", kTH2F, {centAxis, resAxis});
    histos.add("ResTPCRTPCL", "ResTPCRTPCL", kTH2F, {centAxis, resAxis});
    histos.add("ResFT0CFT0A", "ResFT0CFT0A", kTH2F, {centAxis, resAxis});
    histos.add("ResFT0ATPC", "ResFT0ATPC", kTH2F, {centAxis, resAxis});

    // MC histogram
    if (isMC) {
      histos.add("hMC", "MC Event statistics", kTH1F, {{10, 0.0f, 10.0f}});
      histos.add("h1PhiRecsplit", "Phi meson Rec split", kTH1F, {{100, 0.0f, 10.0f}});
      histos.add("CentPercentileMCRecHist", "MC Centrality", kTH1F, {{100, 0.0f, 100.0f}});
      if (!fillRapidity) {
        histos.add("hSparseV2SASameEvent_costhetastarOP_beam_MCGen", "hSparseV2SASameEvent_costhetastar_beamOP_MCGen", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisCosThetaStarOP, thnAxisPhiminusPsi, thnAxisCentrality});
        histos.add("hSparseV2SASameEvent_costhetastarOP_MCGen", "hSparseV2SASameEvent_costhetastarOP_MCGen", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisCosThetaStarOP, thnAxisPhiminusPsi, thnAxisCentrality});
        histos.add("hSparseV2SASameEvent_costhetastarIP_MCGen", "hSparseV2SASameEvent_costhetastarIP_MCGen", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisCosThetaStarOP, thnAxisPhiminusPsi, thnAxisCentrality});
        histos.add("hSparseV2SASameEvent_SA_MCGen", "hSparseV2SASameEvent_SA_MCGen", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisSA, thnAxisPhiminusPsi, thnAxisCentrality});
        histos.add("hSparseV2SASameEvent_costheta_SA_MCGen", "hSparseV2SASameEvent_costheta_SA_MCGen", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxiscosthetaSA, thnAxisPhiminusPsi, thnAxisCentrality});
        histos.add("hSparseV2SASameEvent_SA_A0_MCGen", "hSparseV2SASameEvent_SA_A0_MCGen", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisSA, thnAxisPhiminusPsi, thnAxisCentrality});
        histos.add("hSparseV2SASameEvent_V2_MCGen", "hSparseV2SASameEvent_V2_MCGen", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisV2, thnAxisCentrality});

        histos.add("hSparseV2SASameEvent_costhetastarOP_beam_MCRec", "hSparseV2SASameEvent_costhetastar_beamOP_MCRec", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisCosThetaStarOP, thnAxisPhiminusPsi, thnAxisCentrality});
        histos.add("hSparseV2SASameEvent_costhetastarOP_MCRec", "hSparseV2SASameEvent_costhetastarOP_MCRec", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisCosThetaStarOP, thnAxisPhiminusPsi, thnAxisCentrality});
        histos.add("hSparseV2SASameEvent_costhetastarIP_MCRec", "hSparseV2SASameEvent_costhetastarIP_MCRec", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisCosThetaStarOP, thnAxisPhiminusPsi, thnAxisCentrality});
        histos.add("hSparseV2SASameEvent_SA_MCRec", "hSparseV2SASameEvent_SA_MCRec", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisSA, thnAxisPhiminusPsi, thnAxisCentrality});
        histos.add("hSparseV2SASameEvent_costheta_SA_MCRec", "hSparseV2SASameEvent_costheta_SA_MCRec", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxiscosthetaSA, thnAxisPhiminusPsi, thnAxisCentrality});
        histos.add("hSparseV2SASameEvent_SA_A0_MCRec", "hSparseV2SASameEvent_SA_A0_MCRec", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisSA, thnAxisPhiminusPsi, thnAxisCentrality});
        histos.add("hSparseV2SASameEvent_V2_MCRec", "hSparseV2SASameEvent_V2_MCRec", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisV2, thnAxisCentrality});
      }
      if (fillRapidity) {
        histos.add("hSparseV2SASameEvent_costhetastarOP_beam_MCGen", "hSparseV2SASameEvent_costhetastar_beamOP_MCGen", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisCosThetaStarOP, thnAxisRapidity, thnAxisCentrality});
        histos.add("hSparseV2SASameEvent_costhetastarOP_MCGen", "hSparseV2SASameEvent_costhetastarOP_MCGen", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisCosThetaStarOP, thnAxisRapidity, thnAxisCentrality});
        histos.add("hSparseV2SASameEvent_costhetastarIP_MCGen", "hSparseV2SASameEvent_costhetastarIP_MCGen", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisCosThetaStarOP, thnAxisRapidity, thnAxisCentrality});
        histos.add("hSparseV2SASameEvent_SA_MCGen", "hSparseV2SASameEvent_SA_MCGen", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisSA, thnAxisRapidity, thnAxisCentrality});
        histos.add("hSparseV2SASameEvent_costheta_SA_MCGen", "hSparseV2SASameEvent_costheta_SA_MCGen", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxiscosthetaSA, thnAxisRapidity, thnAxisCentrality});
        histos.add("hSparseV2SASameEvent_SA_A0_MCGen", "hSparseV2SASameEvent_SA_A0_MCGen", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisSA, thnAxisRapidity, thnAxisCentrality});
        histos.add("hSparseV2SASameEvent_V2_MCGen", "hSparseV2SASameEvent_V2_MCGen", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisV2, thnAxisCentrality});

        histos.add("hSparseV2SASameEvent_costhetastarOP_beam_MCRec", "hSparseV2SASameEvent_costhetastar_beamOP_MCRec", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisCosThetaStarOP, thnAxisRapidity, thnAxisCentrality});
        histos.add("hSparseV2SASameEvent_costhetastarOP_MCRec", "hSparseV2SASameEvent_costhetastarOP_MCRec", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisCosThetaStarOP, thnAxisRapidity, thnAxisCentrality});
        histos.add("hSparseV2SASameEvent_costhetastarIP_MCRec", "hSparseV2SASameEvent_costhetastarIP_MCRec", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisCosThetaStarOP, thnAxisRapidity, thnAxisCentrality});
        histos.add("hSparseV2SASameEvent_SA_MCRec", "hSparseV2SASameEvent_SA_MCRec", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisSA, thnAxisRapidity, thnAxisCentrality});
        histos.add("hSparseV2SASameEvent_costheta_SA_MCRec", "hSparseV2SASameEvent_costheta_SA_MCRec", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxiscosthetaSA, thnAxisRapidity, thnAxisCentrality});
        histos.add("hSparseV2SASameEvent_SA_A0_MCRec", "hSparseV2SASameEvent_SA_A0_MCRec", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisSA, thnAxisRapidity, thnAxisCentrality});
        histos.add("hSparseV2SASameEvent_V2_MCRec", "hSparseV2SASameEvent_V2_MCRec", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisV2, thnAxisCentrality});
      }
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
  }

  double massKa = o2::constants::physics::MassKPlus;

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
  bool selectionPIDpTdependent(const T& candidate)
  {
    if (candidate.p() < 0.5 && TMath::Abs(candidate.tpcNSigmaKa()) < nsigmaCutTPC) {
      return true;
    }
    if (candidate.p() >= 0.5 && candidate.hasTOF() && ((candidate.tofNSigmaKa() * candidate.tofNSigmaKa()) + (candidate.tpcNSigmaKa() * candidate.tpcNSigmaKa())) < (nsigmaCutCombined * nsigmaCutCombined)) {
      return true;
    }
    if (!useGlobalTrack && !candidate.hasTPC()) {
      return true;
    }
    return false;
  }
  template <typename T>
  bool selectionPID(const T& candidate)
  {
    if (!isNoTOF && !candidate.hasTOF() && TMath::Abs(candidate.tpcNSigmaKa()) < nsigmaCutTPC) {
      return true;
    }
    if (!isNoTOF && candidate.hasTOF() && ((candidate.tofNSigmaKa() * candidate.tofNSigmaKa()) + (candidate.tpcNSigmaKa() * candidate.tpcNSigmaKa())) < (nsigmaCutCombined * nsigmaCutCombined)) {
      return true;
    }
    if (isNoTOF && TMath::Abs(candidate.tpcNSigmaKa()) < nsigmaCutTPC) {
      return true;
    }
    return false;
  }
  // deep angle cut on pair to remove photon conversion
  template <typename T1, typename T2>
  bool selectionPair(const T1& candidate1, const T2& candidate2)
  {
    double pt1, pt2, pz1, pz2, p1, p2, angle;
    pt1 = candidate1.pt();
    pt2 = candidate2.pt();
    pz1 = candidate1.pz();
    pz2 = candidate2.pz();
    p1 = candidate1.p();
    p2 = candidate2.p();
    angle = TMath::ACos((pt1 * pt2 + pz1 * pz2) / (p1 * p2));
    if (isDeepAngle && angle < cfgDeepAngle) {
      return false;
    }
    return true;
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

  double GetDeltaPsiSubInRange(double psi1, double psi2)
  {
    double delta = psi1 - psi2;
    if (TMath::Abs(delta) > TMath::Pi() / 2) {
      if (delta > 0.)
        delta -= 2. * TMath::Pi() / 2;
      else
        delta += 2. * TMath::Pi() / 2;
    }
    return delta;
  }

  ConfigurableAxis axisVertex{"axisVertex", {20, -10, 10}, "vertex axis for bin"};
  ConfigurableAxis axisMultiplicityClass{"axisMultiplicityClass", {20, 0, 100}, "multiplicity percentile for bin"};
  ConfigurableAxis axisEPAngle{"axisEPAngle", {6, -TMath::Pi() / 2, TMath::Pi() / 2}, "event plane angle"};

  using BinningTypeVertexContributor = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0C, aod::epcalibrationtable::PsiFT0C>;
  ROOT::Math::PxPyPzMVector PhiMesonMother, KaonPlus, KaonMinus, fourVecDauCM;
  ROOT::Math::XYZVector threeVecDauCM, threeVecDauCMXY, eventplaneVec, eventplaneVecNorm, beamvector;

  void processSameEvent(EventCandidates::iterator const& collision, TrackCandidates const& /*tracks*/, aod::BCs const&)
  {
    if (!collision.sel8()) {
      return;
    }
    auto centrality = collision.centFT0C();
    auto multTPC = collision.multNTracksPV();
    histos.fill(HIST("hFTOCvsTPCNoCut"), centrality, multTPC);
    if (!collision.triggereventep()) {
      return;
    }
    if (timFrameEvsel && (!collision.selection_bit(aod::evsel::kNoTimeFrameBorder) || !collision.selection_bit(aod::evsel::kNoITSROFrameBorder))) {
      return;
    }
    auto posThisColl = posTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto negThisColl = negTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto psiFT0C = collision.psiFT0C();
    auto psiFT0A = collision.psiFT0A();
    auto psiTPC = collision.psiTPC();
    auto psiTPCR = collision.psiTPCR();
    auto psiTPCL = collision.psiTPCL();
    histos.fill(HIST("hFTOCvsTPC"), centrality, multTPC);
    if (additionalEvsel && !eventSelected(collision, centrality)) {
      return;
    }
    histos.fill(HIST("hFTOCvsTPCSelected"), centrality, multTPC);
    histos.fill(HIST("hPsiFT0C"), centrality, psiFT0C);
    histos.fill(HIST("hPsiFT0A"), centrality, psiFT0A);
    histos.fill(HIST("hPsiTPC"), centrality, psiTPC);
    histos.fill(HIST("hPsiTPCR"), centrality, psiTPCR);
    histos.fill(HIST("hPsiTPCL"), centrality, psiTPCL);
    histos.fill(HIST("ResFT0CTPC"), centrality, TMath::Cos(2.0 * (psiFT0C - psiTPC)));
    histos.fill(HIST("ResFT0CTPCR"), centrality, TMath::Cos(2.0 * (psiFT0C - psiTPCR)));
    histos.fill(HIST("ResFT0CTPCL"), centrality, TMath::Cos(2.0 * (psiFT0C - psiTPCL)));
    histos.fill(HIST("ResTPCRTPCL"), centrality, TMath::Cos(2.0 * (psiTPCR - psiTPCL)));
    histos.fill(HIST("ResFT0CFT0A"), centrality, TMath::Cos(2.0 * (psiFT0C - psiFT0A)));
    histos.fill(HIST("ResFT0ATPC"), centrality, TMath::Cos(2.0 * (psiTPC - psiFT0A)));
    histos.fill(HIST("hCentrality"), centrality);
    histos.fill(HIST("hVtxZ"), collision.posZ());
    for (auto track1 : posThisColl) {
      // track selection
      if (!selectionTrack(track1)) {
        continue;
      }
      // PID check
      if (ispTdepPID && !selectionPIDpTdependent(track1)) {
        continue;
      }
      if (!ispTdepPID && !selectionPID(track1)) {
        continue;
      }
      histos.fill(HIST("hEta"), track1.eta());
      histos.fill(HIST("hDcaxy"), track1.dcaXY());
      histos.fill(HIST("hDcaz"), track1.dcaZ());
      histos.fill(HIST("hNsigmaKaonTPC"), track1.tpcNSigmaKa());
      histos.fill(HIST("hNsigmaKaonTOF"), track1.tofNSigmaKa());
      auto track1ID = track1.globalIndex();
      for (auto track2 : negThisColl) {
        // track selection
        if (!selectionTrack(track2)) {
          continue;
        }
        // PID check
        if (ispTdepPID && !selectionPIDpTdependent(track2)) {
          continue;
        }
        if (!ispTdepPID && !selectionPID(track2)) {
          continue;
        }
        auto track2ID = track2.globalIndex();
        if (track2ID == track1ID) {
          continue;
        }
        if (!selectionPair(track1, track2)) {
          continue;
        }
        KaonPlus = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), massKa);
        KaonMinus = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), massKa);
        PhiMesonMother = KaonPlus + KaonMinus;
        histos.fill(HIST("hpTvsRapidity"), PhiMesonMother.Pt(), PhiMesonMother.Rapidity());
        if (TMath::Abs(PhiMesonMother.Rapidity()) > confRapidity) {
          continue;
        }
        ROOT::Math::Boost boost{PhiMesonMother.BoostToCM()};
        fourVecDauCM = boost(KaonMinus);
        threeVecDauCM = fourVecDauCM.Vect();
        threeVecDauCMXY = ROOT::Math::XYZVector(threeVecDauCM.X(), threeVecDauCM.Y(), 0.);
        eventplaneVec = ROOT::Math::XYZVector(std::cos(2.0 * psiFT0C), std::sin(2.0 * psiFT0C), 0);
        eventplaneVecNorm = ROOT::Math::XYZVector(std::sin(2.0 * psiFT0C), -std::cos(2.0 * psiFT0C), 0);

        // auto cosinephidaughterstarminuspsi = eventplaneVec.Dot(threeVecDauCMXY) / std::sqrt(threeVecDauCMXY.Mag2()) / std::sqrt(eventplaneVec.Mag2());
        // auto SA = (2.0 * cosinephidaughterstarminuspsi * cosinephidaughterstarminuspsi) - 1.0;
        auto cosPhistarminuspsi = GetPhiInRange(fourVecDauCM.Phi() - psiFT0C);
        auto SA = TMath::Cos(2.0 * cosPhistarminuspsi);
        // auto cosThetaStarOP = TMath::Abs(eventplaneVecNorm.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2()) / std::sqrt(eventplaneVecNorm.Mag2()));
        auto cosThetaStarOP = eventplaneVecNorm.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2()) / std::sqrt(eventplaneVecNorm.Mag2());
        auto SA_A0 = 1 - (cosThetaStarOP * cosThetaStarOP);
        // auto cosThetaStarIP = TMath::Abs(eventplaneVec.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2()) / std::sqrt(eventplaneVec.Mag2()));
        auto cosThetaStarIP = eventplaneVec.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2()) / std::sqrt(eventplaneVec.Mag2());
        auto phiminuspsi = GetPhiInRange(PhiMesonMother.Phi() - psiFT0C);
        auto v2 = TMath::Cos(2.0 * phiminuspsi);
        if (!fillRapidity) {
          histos.fill(HIST("hSparseV2SASameEvent_costhetastarOP"), PhiMesonMother.M(), PhiMesonMother.Pt(), cosThetaStarOP, phiminuspsi, centrality);
          histos.fill(HIST("hSparseV2SASameEvent_costheta_SA"), PhiMesonMother.M(), PhiMesonMother.Pt(), cosThetaStarOP * cosThetaStarOP, phiminuspsi, centrality);
          histos.fill(HIST("hSparseV2SASameEvent_costhetastarIP"), PhiMesonMother.M(), PhiMesonMother.Pt(), cosThetaStarIP, phiminuspsi, centrality);
          histos.fill(HIST("hSparseV2SASameEvent_SA"), PhiMesonMother.M(), PhiMesonMother.Pt(), SA, phiminuspsi, centrality);
          histos.fill(HIST("hSparseV2SASameEvent_SA_A0"), PhiMesonMother.M(), PhiMesonMother.Pt(), SA_A0, phiminuspsi, centrality);
          histos.fill(HIST("hSparseV2SASameEvent_V2"), PhiMesonMother.M(), PhiMesonMother.Pt(), v2, centrality);
        }

        if (fillRapidity) {
          histos.fill(HIST("hSparseV2SASameEvent_costhetastarOP"), PhiMesonMother.M(), PhiMesonMother.Pt(), cosThetaStarOP, TMath::Abs(PhiMesonMother.Rapidity()), centrality);
          histos.fill(HIST("hSparseV2SASameEvent_costheta_SA"), PhiMesonMother.M(), PhiMesonMother.Pt(), cosThetaStarOP * cosThetaStarOP, TMath::Abs(PhiMesonMother.Rapidity()), centrality);
          histos.fill(HIST("hSparseV2SASameEvent_costhetastarIP"), PhiMesonMother.M(), PhiMesonMother.Pt(), cosThetaStarIP, TMath::Abs(PhiMesonMother.Rapidity()), centrality);
          histos.fill(HIST("hSparseV2SASameEvent_SA"), PhiMesonMother.M(), PhiMesonMother.Pt(), SA, TMath::Abs(PhiMesonMother.Rapidity()), centrality);
          histos.fill(HIST("hSparseV2SASameEvent_SA_A0"), PhiMesonMother.M(), PhiMesonMother.Pt(), SA_A0, TMath::Abs(PhiMesonMother.Rapidity()), centrality);
          histos.fill(HIST("hSparseV2SASameEvent_V2"), PhiMesonMother.M(), PhiMesonMother.Pt(), v2, centrality);
        }
      }
    }
  }
  PROCESS_SWITCH(phipbpb, processSameEvent, "Process Same event", true);
  void processMixedEvent(EventCandidates const& collisions, TrackCandidates const& /*tracks*/)
  {
    BinningTypeVertexContributor binningOnPositions{{axisVertex, axisMultiplicityClass, axisEPAngle}, true};
    for (auto const& [collision1, collision2] : o2::soa::selfCombinations(binningOnPositions, cfgNoMixedEvents, -1, collisions, collisions)) {
      if (!collision1.sel8() || !collision2.sel8()) {
        // printf("Mix = %d\n", 1);
        continue;
      }
      if (!collision1.triggereventep() || !collision2.triggereventep()) {
        // printf("Mix = %d\n", 2);
        continue;
      }
      if (timFrameEvsel && (!collision1.selection_bit(aod::evsel::kNoTimeFrameBorder) || !collision2.selection_bit(aod::evsel::kNoTimeFrameBorder) || !collision1.selection_bit(aod::evsel::kNoITSROFrameBorder) || !collision2.selection_bit(aod::evsel::kNoITSROFrameBorder))) {
        // printf("Mix = %d\n", 3);
        continue;
      }
      auto posThisColl = posTracks->sliceByCached(aod::track::collisionId, collision1.globalIndex(), cache);
      auto negThisColl = negTracks->sliceByCached(aod::track::collisionId, collision2.globalIndex(), cache);
      auto centrality = collision1.centFT0C();
      auto centrality2 = collision2.centFT0C();
      auto psiFT0C = collision1.psiFT0C();

      if (additionalEvsel && !eventSelected(collision1, centrality)) {
        // printf("Mix = %d\n", 4);
        continue;
      }
      if (additionalEvsel && !eventSelected(collision2, centrality2)) {
        // printf("Mix = %d\n", 5);
        continue;
      }
      for (auto& [track1, track2] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(posThisColl, negThisColl))) {
        // track selection
        if (!selectionTrack(track1) || !selectionTrack(track2)) {
          // printf("Mix = %d\n", 6);
          continue;
        }
        // PID check
        if (ispTdepPID && (!selectionPIDpTdependent(track1) || !selectionPIDpTdependent(track2))) {
          // printf("Mix = %d\n", 7);
          continue;
        }
        if (!ispTdepPID && (!selectionPID(track1) || !selectionPID(track2))) {
          continue;
        }
        if (!selectionPair(track1, track2)) {
          // printf("Mix = %d\n", 8);
          continue;
        }
        KaonPlus = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), massKa);
        KaonMinus = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), massKa);
        PhiMesonMother = KaonPlus + KaonMinus;
        if (TMath::Abs(PhiMesonMother.Rapidity()) > confRapidity) {
          continue;
        }
        ROOT::Math::Boost boost{PhiMesonMother.BoostToCM()};
        fourVecDauCM = boost(KaonMinus);
        threeVecDauCM = fourVecDauCM.Vect();
        threeVecDauCMXY = ROOT::Math::XYZVector(threeVecDauCM.X(), threeVecDauCM.Y(), 0.);
        eventplaneVec = ROOT::Math::XYZVector(std::cos(2.0 * psiFT0C), std::sin(2.0 * psiFT0C), 0);
        eventplaneVecNorm = ROOT::Math::XYZVector(std::sin(2.0 * psiFT0C), -std::cos(2.0 * psiFT0C), 0);

        // auto cosinephidaughterstarminuspsi = eventplaneVec.Dot(threeVecDauCMXY) / std::sqrt(threeVecDauCMXY.Mag2()) / std::sqrt(eventplaneVec.Mag2());
        // auto SA = (2.0 * cosinephidaughterstarminuspsi * cosinephidaughterstarminuspsi) - 1.0;
        auto cosPhistarminuspsi = GetPhiInRange(fourVecDauCM.Phi() - psiFT0C);
        auto SA = TMath::Cos(2.0 * cosPhistarminuspsi);
        // auto cosThetaStarOP = TMath::Abs(eventplaneVecNorm.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2()) / std::sqrt(eventplaneVecNorm.Mag2()));
        auto cosThetaStarOP = eventplaneVecNorm.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2()) / std::sqrt(eventplaneVecNorm.Mag2());
        auto SA_A0 = 1 - (cosThetaStarOP * cosThetaStarOP);
        // auto cosThetaStarIP = TMath::Abs(eventplaneVec.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2()) / std::sqrt(eventplaneVec.Mag2()));
        auto cosThetaStarIP = eventplaneVec.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2()) / std::sqrt(eventplaneVec.Mag2());
        auto phiminuspsi = GetPhiInRange(PhiMesonMother.Phi() - psiFT0C);
        auto v2 = TMath::Cos(2.0 * phiminuspsi);
        if (!fillRapidity) {
          histos.fill(HIST("hSparseV2SAMixedEvent_costhetastarOP"), PhiMesonMother.M(), PhiMesonMother.Pt(), cosThetaStarOP, phiminuspsi, centrality);
          histos.fill(HIST("hSparseV2SAMixedEvent_costheta_SA"), PhiMesonMother.M(), PhiMesonMother.Pt(), cosThetaStarOP * cosThetaStarOP, phiminuspsi, centrality);
          histos.fill(HIST("hSparseV2SAMixedEvent_costhetastarIP"), PhiMesonMother.M(), PhiMesonMother.Pt(), cosThetaStarIP, phiminuspsi, centrality);
          histos.fill(HIST("hSparseV2SAMixedEvent_SA"), PhiMesonMother.M(), PhiMesonMother.Pt(), SA, phiminuspsi, centrality);
          histos.fill(HIST("hSparseV2SAMixedEvent_SA_A0"), PhiMesonMother.M(), PhiMesonMother.Pt(), SA_A0, phiminuspsi, centrality);
          histos.fill(HIST("hSparseV2SAMixedEvent_V2"), PhiMesonMother.M(), PhiMesonMother.Pt(), v2, centrality);
        }
        if (fillRapidity) {
          histos.fill(HIST("hSparseV2SAMixedEvent_costhetastarOP"), PhiMesonMother.M(), PhiMesonMother.Pt(), cosThetaStarOP, TMath::Abs(PhiMesonMother.Rapidity()), centrality);
          histos.fill(HIST("hSparseV2SAMixedEvent_costheta_SA"), PhiMesonMother.M(), PhiMesonMother.Pt(), cosThetaStarOP * cosThetaStarOP, TMath::Abs(PhiMesonMother.Rapidity()), centrality);
          histos.fill(HIST("hSparseV2SAMixedEvent_costhetastarIP"), PhiMesonMother.M(), PhiMesonMother.Pt(), cosThetaStarIP, TMath::Abs(PhiMesonMother.Rapidity()), centrality);
          histos.fill(HIST("hSparseV2SAMixedEvent_SA"), PhiMesonMother.M(), PhiMesonMother.Pt(), SA, TMath::Abs(PhiMesonMother.Rapidity()), centrality);
          histos.fill(HIST("hSparseV2SAMixedEvent_SA_A0"), PhiMesonMother.M(), PhiMesonMother.Pt(), SA_A0, TMath::Abs(PhiMesonMother.Rapidity()), centrality);
          histos.fill(HIST("hSparseV2SAMixedEvent_V2"), PhiMesonMother.M(), PhiMesonMother.Pt(), v2, centrality);
        }
      }
    }
  }
  PROCESS_SWITCH(phipbpb, processMixedEvent, "Process Mixed event", true);
  void processMixedEventOpti(EventCandidates const& collisions, TrackCandidates const& tracks)
  {
    auto tracksTuple = std::make_tuple(tracks);
    BinningTypeVertexContributor binningOnPositions{{axisVertex, axisMultiplicityClass, axisEPAngle}, true};
    SameKindPair<EventCandidates, TrackCandidates, BinningTypeVertexContributor> pair{binningOnPositions, cfgNoMixedEvents, -1, collisions, tracksTuple, &cache};
    for (auto& [collision1, tracks1, collision2, tracks2] : pair) {
      if (!collision1.sel8() || !collision2.sel8()) {
        continue;
      }
      if (!collision1.triggereventep() || !collision2.triggereventep()) {
        continue;
      }
      if (timFrameEvsel && (!collision1.selection_bit(aod::evsel::kNoTimeFrameBorder) || !collision2.selection_bit(aod::evsel::kNoTimeFrameBorder) || !collision1.selection_bit(aod::evsel::kNoITSROFrameBorder) || !collision2.selection_bit(aod::evsel::kNoITSROFrameBorder))) {
        continue;
      }
      auto centrality = collision1.centFT0C();
      auto centrality2 = collision2.centFT0C();
      auto psiFT0C = collision1.psiFT0C();
      if (additionalEvsel && !eventSelected(collision1, centrality)) {
        continue;
      }
      if (additionalEvsel && !eventSelected(collision2, centrality2)) {
        continue;
      }

      for (auto& [track1, track2] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {
        if (track1.sign() * track2.sign() > 0) {
          continue;
        }
        if (!selectionTrack(track1) || !selectionTrack(track2)) {
          continue;
        }
        // PID check
        if (ispTdepPID && (!selectionPIDpTdependent(track1) || !selectionPIDpTdependent(track2))) {
          continue;
        }
        if (!ispTdepPID && (!selectionPID(track1) || !selectionPID(track2))) {
          continue;
        }
        if (!selectionPair(track1, track2)) {
          continue;
        }
        if (track1.sign() > 0 && track2.sign() < 0) {
          KaonPlus = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), massKa);
          KaonMinus = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), massKa);
        } else if (track1.sign() < 0 && track2.sign() > 0) {
          KaonMinus = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), massKa);
          KaonPlus = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), massKa);
        }

        PhiMesonMother = KaonPlus + KaonMinus;
        if (TMath::Abs(PhiMesonMother.Rapidity()) > confRapidity) {
          continue;
        }
        ROOT::Math::Boost boost{PhiMesonMother.BoostToCM()};
        fourVecDauCM = boost(KaonMinus);
        threeVecDauCM = fourVecDauCM.Vect();
        threeVecDauCMXY = ROOT::Math::XYZVector(threeVecDauCM.X(), threeVecDauCM.Y(), 0.);
        eventplaneVec = ROOT::Math::XYZVector(std::cos(2.0 * psiFT0C), std::sin(2.0 * psiFT0C), 0);
        eventplaneVecNorm = ROOT::Math::XYZVector(std::sin(2.0 * psiFT0C), -std::cos(2.0 * psiFT0C), 0);

        auto cosPhistarminuspsi = GetPhiInRange(fourVecDauCM.Phi() - psiFT0C);
        auto SA = TMath::Cos(2.0 * cosPhistarminuspsi);
        // auto cosThetaStarOP = TMath::Abs(eventplaneVecNorm.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2()) / std::sqrt(eventplaneVecNorm.Mag2()));
        auto cosThetaStarOP = eventplaneVecNorm.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2()) / std::sqrt(eventplaneVecNorm.Mag2());
        auto SA_A0 = 1 - (cosThetaStarOP * cosThetaStarOP);
        // auto cosThetaStarIP = TMath::Abs(eventplaneVec.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2()) / std::sqrt(eventplaneVec.Mag2()));
        auto cosThetaStarIP = eventplaneVec.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2()) / std::sqrt(eventplaneVec.Mag2());
        auto phiminuspsi = GetPhiInRange(PhiMesonMother.Phi() - psiFT0C);
        auto v2 = TMath::Cos(2.0 * phiminuspsi);
        if (!fillRapidity) {
          histos.fill(HIST("hSparseV2SAMixedEvent_costhetastarOP"), PhiMesonMother.M(), PhiMesonMother.Pt(), cosThetaStarOP, phiminuspsi, centrality);
          histos.fill(HIST("hSparseV2SAMixedEvent_costheta_SA"), PhiMesonMother.M(), PhiMesonMother.Pt(), cosThetaStarOP * cosThetaStarOP, phiminuspsi, centrality);
          histos.fill(HIST("hSparseV2SAMixedEvent_costhetastarIP"), PhiMesonMother.M(), PhiMesonMother.Pt(), cosThetaStarIP, phiminuspsi, centrality);
          histos.fill(HIST("hSparseV2SAMixedEvent_SA"), PhiMesonMother.M(), PhiMesonMother.Pt(), SA, phiminuspsi, centrality);
          histos.fill(HIST("hSparseV2SAMixedEvent_SA_A0"), PhiMesonMother.M(), PhiMesonMother.Pt(), SA_A0, phiminuspsi, centrality);
          histos.fill(HIST("hSparseV2SAMixedEvent_V2"), PhiMesonMother.M(), PhiMesonMother.Pt(), v2, centrality);
        }
        if (fillRapidity) {
          histos.fill(HIST("hSparseV2SAMixedEvent_costhetastarOP"), PhiMesonMother.M(), PhiMesonMother.Pt(), cosThetaStarOP, TMath::Abs(PhiMesonMother.Rapidity()), centrality);
          histos.fill(HIST("hSparseV2SAMixedEvent_costheta_SA"), PhiMesonMother.M(), PhiMesonMother.Pt(), cosThetaStarOP * cosThetaStarOP, TMath::Abs(PhiMesonMother.Rapidity()), centrality);
          histos.fill(HIST("hSparseV2SAMixedEvent_costhetastarIP"), PhiMesonMother.M(), PhiMesonMother.Pt(), cosThetaStarIP, TMath::Abs(PhiMesonMother.Rapidity()), centrality);
          histos.fill(HIST("hSparseV2SAMixedEvent_SA"), PhiMesonMother.M(), PhiMesonMother.Pt(), SA, TMath::Abs(PhiMesonMother.Rapidity()), centrality);
          histos.fill(HIST("hSparseV2SAMixedEvent_SA_A0"), PhiMesonMother.M(), PhiMesonMother.Pt(), SA_A0, TMath::Abs(PhiMesonMother.Rapidity()), centrality);
          histos.fill(HIST("hSparseV2SAMixedEvent_V2"), PhiMesonMother.M(), PhiMesonMother.Pt(), v2, centrality);
        }
      }
    }
  }
  PROCESS_SWITCH(phipbpb, processMixedEventOpti, "Process Mixed event new", true);
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
      if (timFrameEvsel && (!RecCollision.selection_bit(aod::evsel::kNoTimeFrameBorder) || !RecCollision.selection_bit(aod::evsel::kNoITSROFrameBorder))) {
        histos.fill(HIST("hMC"), 5);
        continue;
      }
      if (std::abs(RecCollision.posZ()) > cfgCutVertex) {
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
        if (ispTdepPID && !selectionPIDpTdependent(track1)) {
          continue;
        }
        if (!ispTdepPID && !selectionPID(track1)) {
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
          if (ispTdepPID && !selectionPIDpTdependent(track2)) {
            continue;
          }
          if (!ispTdepPID && !selectionPID(track2)) {
            continue;
          }
          if (!track2.has_mcParticle()) {
            continue;
          }
          if (!selectionPair(track1, track2)) {
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
          if (!(track1PDG == 321 && track2PDG == 321)) {
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
              if (std::abs(mothertrack1.y()) > confRapidity) {
                continue;
              }
              if (std::abs(mothertrack1.pdgCode()) != 333) {
                continue;
              }
              if (!selectionPID(track1) || !selectionPID(track2)) {
                continue;
              }
              if (avoidsplitrackMC && oldindex == mothertrack1.globalIndex()) {
                histos.fill(HIST("h1PhiRecsplit"), mothertrack1.pt());
                continue;
              }
              oldindex = mothertrack1.globalIndex();
              if (track1.sign() > 0 && track2.sign() < 0) {
                KaonPlus = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), massKa);
                KaonMinus = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), massKa);
              }
              if (track1.sign() < 0 && track2.sign() > 0) {
                KaonMinus = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), massKa);
                KaonPlus = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), massKa);
              }
              PhiMesonMother = KaonPlus + KaonMinus;

              if (TMath::Abs(PhiMesonMother.Rapidity()) > confRapidity) {
                continue;
              }
              ROOT::Math::Boost boost{PhiMesonMother.BoostToCM()};
              fourVecDauCM = boost(KaonMinus);
              threeVecDauCM = fourVecDauCM.Vect();
              beamvector = ROOT::Math::XYZVector(0.0, 0.0, 1.0);
              threeVecDauCMXY = ROOT::Math::XYZVector(threeVecDauCM.X(), threeVecDauCM.Y(), 0.);
              eventplaneVec = ROOT::Math::XYZVector(std::cos(2.0 * psiFT0C), std::sin(2.0 * psiFT0C), 0);
              eventplaneVecNorm = ROOT::Math::XYZVector(std::sin(2.0 * psiFT0C), -std::cos(2.0 * psiFT0C), 0);
              auto cosPhistarminuspsi = GetPhiInRange(fourVecDauCM.Phi() - psiFT0C);
              auto SA = TMath::Cos(2.0 * cosPhistarminuspsi);
              // auto cosThetaStarOP = TMath::Abs(eventplaneVecNorm.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2()) / std::sqrt(eventplaneVecNorm.Mag2()));
              auto cosThetaStarIP = eventplaneVec.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2()) / std::sqrt(eventplaneVec.Mag2());
              auto cosThetaStarOP = eventplaneVecNorm.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2()) / std::sqrt(eventplaneVecNorm.Mag2());
              auto cosThetaStarOPbeam = beamvector.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2()) / std::sqrt(beamvector.Mag2());
              auto SA_A0 = 1 - (cosThetaStarOP * cosThetaStarOP);
              // auto cosThetaStarIP = TMath::Abs(eventplaneVec.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2()) / std::sqrt(eventplaneVec.Mag2()));
              auto phiminuspsi = GetPhiInRange(PhiMesonMother.Phi() - psiFT0C);
              auto v2 = TMath::Cos(2.0 * phiminuspsi);
              if (!fillRapidity) {
                histos.fill(HIST("hSparseV2SASameEvent_costhetastarOP_MCRec"), PhiMesonMother.M(), PhiMesonMother.Pt(), cosThetaStarOP, phiminuspsi, centrality);
                histos.fill(HIST("hSparseV2SASameEvent_costhetastarOP_beam_MCRec"), PhiMesonMother.M(), PhiMesonMother.Pt(), cosThetaStarOPbeam, phiminuspsi, centrality);
                histos.fill(HIST("hSparseV2SASameEvent_costheta_SA_MCRec"), PhiMesonMother.M(), PhiMesonMother.Pt(), cosThetaStarOP * cosThetaStarOP, phiminuspsi, centrality);
                histos.fill(HIST("hSparseV2SASameEvent_costhetastarIP_MCRec"), PhiMesonMother.M(), PhiMesonMother.Pt(), cosThetaStarIP, phiminuspsi, centrality);
                histos.fill(HIST("hSparseV2SASameEvent_SA_MCRec"), PhiMesonMother.M(), PhiMesonMother.Pt(), SA, phiminuspsi, centrality);
                histos.fill(HIST("hSparseV2SASameEvent_SA_A0_MCRec"), PhiMesonMother.M(), PhiMesonMother.Pt(), SA_A0, phiminuspsi, centrality);
                histos.fill(HIST("hSparseV2SASameEvent_V2_MCRec"), PhiMesonMother.M(), PhiMesonMother.Pt(), v2, centrality);
              }
              if (fillRapidity) {
                histos.fill(HIST("hSparseV2SASameEvent_costhetastarOP_MCRec"), PhiMesonMother.M(), PhiMesonMother.Pt(), cosThetaStarOP, TMath::Abs(PhiMesonMother.Rapidity()), centrality);
                histos.fill(HIST("hSparseV2SASameEvent_costhetastarOP_beam_MCRec"), PhiMesonMother.M(), PhiMesonMother.Pt(), cosThetaStarOPbeam, TMath::Abs(PhiMesonMother.Rapidity()), centrality);
                histos.fill(HIST("hSparseV2SASameEvent_costheta_SA_MCRec"), PhiMesonMother.M(), PhiMesonMother.Pt(), cosThetaStarOP * cosThetaStarOP, TMath::Abs(PhiMesonMother.Rapidity()), centrality);
                histos.fill(HIST("hSparseV2SASameEvent_costhetastarIP_MCRec"), PhiMesonMother.M(), PhiMesonMother.Pt(), cosThetaStarIP, TMath::Abs(PhiMesonMother.Rapidity()), centrality);
                histos.fill(HIST("hSparseV2SASameEvent_SA_MCRec"), PhiMesonMother.M(), PhiMesonMother.Pt(), SA, TMath::Abs(PhiMesonMother.Rapidity()), centrality);
                histos.fill(HIST("hSparseV2SASameEvent_SA_A0_MCRec"), PhiMesonMother.M(), PhiMesonMother.Pt(), SA_A0, TMath::Abs(PhiMesonMother.Rapidity()), centrality);
                histos.fill(HIST("hSparseV2SASameEvent_V2_MCRec"), PhiMesonMother.M(), PhiMesonMother.Pt(), v2, centrality);
              }
            }
          }
        }
      }
      // loop over generated particle
      for (auto& mcParticle : GenParticles) {
        if (std::abs(mcParticle.y()) > confRapidity) {
          continue;
        }
        if (mcParticle.pdgCode() != 333) {
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
          } else if (kCurrentDaughter.pdgCode() == -321) {
            if (genacceptancecut && kCurrentDaughter.pt() > cfgCutPT && TMath::Abs(kCurrentDaughter.eta()) < cfgCutEta) {
              daughtm = true;
            }
            if (!genacceptancecut) {
              daughtm = true;
            }
            KaonMinus = ROOT::Math::PxPyPzMVector(kCurrentDaughter.px(), kCurrentDaughter.py(), kCurrentDaughter.pz(), massKa);
          }
        }
        if (daughtp && daughtm) {
          PhiMesonMother = KaonPlus + KaonMinus;
          ROOT::Math::Boost boost{PhiMesonMother.BoostToCM()};
          fourVecDauCM = boost(KaonMinus);
          threeVecDauCM = fourVecDauCM.Vect();
          threeVecDauCMXY = ROOT::Math::XYZVector(threeVecDauCM.X(), threeVecDauCM.Y(), 0.);
          beamvector = ROOT::Math::XYZVector(0.0, 0.0, 1.0);
          eventplaneVec = ROOT::Math::XYZVector(std::cos(2.0 * psiFT0C), std::sin(2.0 * psiFT0C), 0);
          eventplaneVecNorm = ROOT::Math::XYZVector(std::sin(2.0 * psiFT0C), -std::cos(2.0 * psiFT0C), 0);
          auto cosPhistarminuspsi = GetPhiInRange(fourVecDauCM.Phi() - psiFT0C);
          auto SA = TMath::Cos(2.0 * cosPhistarminuspsi);
          auto cosThetaStarIP = eventplaneVec.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2()) / std::sqrt(eventplaneVec.Mag2());
          auto cosThetaStarOP = eventplaneVecNorm.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2()) / std::sqrt(eventplaneVecNorm.Mag2());
          auto cosThetaStarOPbeam = beamvector.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2()) / std::sqrt(beamvector.Mag2());
          auto SA_A0 = 1 - (cosThetaStarOP * cosThetaStarOP);
          auto phiminuspsi = GetPhiInRange(PhiMesonMother.Phi() - psiFT0C);
          auto v2 = TMath::Cos(2.0 * phiminuspsi);
          if (!fillRapidity) {
            histos.fill(HIST("hSparseV2SASameEvent_costhetastarOP_MCGen"), PhiMesonMother.M(), PhiMesonMother.Pt(), cosThetaStarOP, phiminuspsi, centrality);
            histos.fill(HIST("hSparseV2SASameEvent_costhetastarOP_beam_MCGen"), PhiMesonMother.M(), PhiMesonMother.Pt(), cosThetaStarOPbeam, phiminuspsi, centrality);
            histos.fill(HIST("hSparseV2SASameEvent_costheta_SA_MCGen"), PhiMesonMother.M(), PhiMesonMother.Pt(), cosThetaStarOP * cosThetaStarOP, phiminuspsi, centrality);
            histos.fill(HIST("hSparseV2SASameEvent_costhetastarIP_MCGen"), PhiMesonMother.M(), PhiMesonMother.Pt(), cosThetaStarIP, phiminuspsi, centrality);
            histos.fill(HIST("hSparseV2SASameEvent_SA_MCGen"), PhiMesonMother.M(), PhiMesonMother.Pt(), SA, phiminuspsi, centrality);
            histos.fill(HIST("hSparseV2SASameEvent_SA_A0_MCGen"), PhiMesonMother.M(), PhiMesonMother.Pt(), SA_A0, phiminuspsi, centrality);
            histos.fill(HIST("hSparseV2SASameEvent_V2_MCGen"), PhiMesonMother.M(), PhiMesonMother.Pt(), v2, centrality);
          }
          if (fillRapidity) {
            histos.fill(HIST("hSparseV2SASameEvent_costhetastarOP_MCGen"), PhiMesonMother.M(), PhiMesonMother.Pt(), cosThetaStarOP, TMath::Abs(PhiMesonMother.Rapidity()), centrality);
            histos.fill(HIST("hSparseV2SASameEvent_costhetastarOP_beam_MCGen"), PhiMesonMother.M(), PhiMesonMother.Pt(), cosThetaStarOPbeam, TMath::Abs(PhiMesonMother.Rapidity()), centrality);
            histos.fill(HIST("hSparseV2SASameEvent_costheta_SA_MCGen"), PhiMesonMother.M(), PhiMesonMother.Pt(), cosThetaStarOP * cosThetaStarOP, TMath::Abs(PhiMesonMother.Rapidity()), centrality);
            histos.fill(HIST("hSparseV2SASameEvent_costhetastarIP_MCGen"), PhiMesonMother.M(), PhiMesonMother.Pt(), cosThetaStarIP, TMath::Abs(PhiMesonMother.Rapidity()), centrality);
            histos.fill(HIST("hSparseV2SASameEvent_SA_MCGen"), PhiMesonMother.M(), PhiMesonMother.Pt(), SA, TMath::Abs(PhiMesonMother.Rapidity()), centrality);
            histos.fill(HIST("hSparseV2SASameEvent_SA_A0_MCGen"), PhiMesonMother.M(), PhiMesonMother.Pt(), SA_A0, TMath::Abs(PhiMesonMother.Rapidity()), centrality);
            histos.fill(HIST("hSparseV2SASameEvent_V2_MCGen"), PhiMesonMother.M(), PhiMesonMother.Pt(), v2, centrality);
          }
        }
      }
    } // rec collision loop

  } // process MC
  PROCESS_SWITCH(phipbpb, processMC, "Process MC", false);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<phipbpb>(cfgc, TaskName{"phipbpb"})};
}
