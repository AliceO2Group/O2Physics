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
#include <vector>
#include <string>

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
#include "CCDB/CcdbApi.h"
#include "Common/DataModel/PIDResponseITS.h"
#include "PWGMM/Mult/DataModel/Index.h" // for Particles2Tracks table

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;
struct phipbpb {

  int mRunNumber;
  int multEstimator;
  float d_bz;

  struct : ConfigurableGroup {
    Configurable<std::string> cfgURL{"cfgURL", "http://alice-ccdb.cern.ch", "Address of the CCDB to browse"};
    Configurable<int64_t> nolaterthan{"ccdb-no-later-than", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "Latest acceptable timestamp of creation for the object"};
  } cfgCcdbParam;

  // Enable access to the CCDB for the offset and correction constants and save them in dedicated variables.
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::ccdb::CcdbApi ccdbApi;
  Service<o2::framework::O2DatabasePDG> pdg;

  // CCDB options
  // Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  // Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  // Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  // Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
  // Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};

  // filling
  Configurable<bool> fillSA{"fillSA", false, "fill spin alignment"};
  // events
  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};
  Configurable<float> cfgCutCentrality{"cfgCutCentrality", 80.0f, "Accepted maximum Centrality"};
  Configurable<int> cfgCutOccupancy{"cfgCutOccupancy", 3000, "Occupancy cut"};
  // track
  Configurable<bool> useSP{"useSP", false, "use SP"};
  Configurable<bool> additionalEvsel{"additionalEvsel", false, "Additional event selcection"};
  Configurable<bool> additionalEvselITS{"additionalEvselITS", true, "Additional event selcection for ITS"};
  Configurable<bool> removefaketrak{"removefaketrack", true, "Remove fake track from momentum difference"};
  Configurable<float> ConfFakeKaonCut{"ConfFakeKaonCut", 0.1, "Cut based on track from momentum difference"};
  Configurable<bool> useGlobalTrack{"useGlobalTrack", true, "use Global track"};
  Configurable<float> cfgCutTOFBeta{"cfgCutTOFBeta", 0.0, "cut TOF beta"};
  Configurable<float> cfgCutCharge{"cfgCutCharge", 0.0, "cut on Charge"};
  Configurable<float> cfgCutPT{"cfgCutPT", 0.2, "PT cut on daughter track"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8, "Eta cut on daughter track"};
  Configurable<float> cfgCutDCAxy{"cfgCutDCAxy", 2.0f, "DCAxy range for tracks"};
  Configurable<float> cfgCutDCAz{"cfgCutDCAz", 2.0f, "DCAz range for tracks"};
  Configurable<float> nsigmaCutTPC{"nsigmacutTPC", 3.0, "Value of the TPC Nsigma cut"};
  Configurable<float> nsigmaCutTOF{"nsigmaCutTOF", 3.0, "Value of the TOF Nsigma cut"};
  Configurable<int> cfgNoMixedEvents{"cfgNoMixedEvents", 1, "Number of mixed events per event"};
  Configurable<int> cfgITScluster{"cfgITScluster", 0, "Number of ITS cluster"};
  Configurable<int> cfgTPCcluster{"cfgTPCcluster", 70, "Number of TPC cluster"};
  Configurable<float> cfgTPCSharedcluster{"cfgTPCSharedcluster", 0.4, "Maximum Number of TPC shared cluster"};
  Configurable<bool> isDeepAngle{"isDeepAngle", false, "Deep Angle cut"};
  Configurable<bool> ispTdepPID{"ispTdepPID", true, "pT dependent PID"};
  Configurable<bool> isTOFOnly{"isTOFOnly", false, "use TOF only PID"};
  Configurable<bool> checkAllCharge{"checkAllCharge", true, "check all charge for MC weight"};
  Configurable<double> cfgDeepAngle{"cfgDeepAngle", 0.04, "Deep Angle cut value"};
  Configurable<double> confRapidity{"confRapidity", 0.5, "Rapidity cut"};
  ConfigurableAxis configThnAxisInvMass{"configThnAxisInvMass", {120, 0.98, 1.1}, "#it{M} (GeV/#it{c}^{2})"};
  ConfigurableAxis configThnAxisPt{"configThnAxisPt", {100, 0.0, 10.}, "#it{p}_{T} (GeV/#it{c})"};
  ConfigurableAxis configThnAxisCosThetaStar{"configThnAxisCosThetaStar", {10, -1.0, 1.}, "cos(#vartheta)"};
  ConfigurableAxis configThnAxisCentrality{"configThnAxisCentrality", {8, 0., 80}, "Centrality"};
  ConfigurableAxis configThnAxisPhiminusPsi{"configThnAxisPhiminusPsi", {6, 0.0, TMath::Pi()}, "#phi - #psi"};
  ConfigurableAxis configThnAxisV2{"configThnAxisV2", {200, -6, 6}, "V2"};
  ConfigurableAxis configThnAxisRapidity{"configThnAxisRapidity", {8, 0, 0.8}, "Rapidity"};
  ConfigurableAxis configThnAxisSA{"configThnAxisSA", {200, -1, 1}, "SA"};
  ConfigurableAxis configThnAxiscosthetaSA{"configThnAxiscosthetaSA", {200, 0, 1}, "costhetaSA"};
  ConfigurableAxis axisPtKaonWeight{"axisPtKaonWeight", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f}, "pt axis"};
  Configurable<bool> isMC{"isMC", false, "use MC"};
  Configurable<bool> genacceptancecut{"genacceptancecut", true, "use acceptance cut for generated"};
  Configurable<bool> avoidsplitrackMC{"avoidsplitrackMC", false, "avoid split track in MC"};
  Configurable<bool> islike{"islike", false, "use like"};
  Configurable<bool> useWeight{"useWeight", true, "use EP dep effi weight"};
  Configurable<std::string> ConfWeightPath{"ConfWeightPath", "Users/s/skundu/My/Object/mcweight", "Path to gain calibration"};

  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter centralityFilter = nabs(aod::cent::centFT0C) < cfgCutCentrality;
  Filter acceptanceFilter = (nabs(aod::track::eta) < cfgCutEta && nabs(aod::track::pt) > cfgCutPT);
  Filter DCAcutFilter = (nabs(aod::track::dcaXY) < cfgCutDCAxy) && (nabs(aod::track::dcaZ) < cfgCutDCAz);
  Filter PIDcutFilter = nabs(aod::pidtpc::tpcNSigmaKa) < nsigmaCutTPC;

  using EventCandidates = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::FV0Mults, aod::TPCMults, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0Cs, aod::CentFT0As, aod::EPCalibrationTables, aod::Mults>>;
  using TrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTOFbeta, aod::pidTPCFullKa, aod::pidTOFFullKa>>;

  using CollisionMCTrueTable = aod::McCollisions;
  using TrackMCTrueTable = aod::McParticles;

  using CollisionMCRecTableCentFT0C = soa::SmallGroups<soa::Join<aod::McCollisionLabels, aod::Collisions, aod::CentFT0Cs, aod::EvSels>>;
  using TrackMCRecTable = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::McTrackLabels, aod::TrackSelection, aod::pidTOFbeta, aod::pidTPCFullKa, aod::pidTOFFullKa>;
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
    // std::vector<double> occupancyBinning = {0.0, 500.0, 1000.0, 3000.0, 6000.0, 50000.0};
    std::vector<double> occupancyBinning = {-0.5, 500.0, 1000.0, 1500.0, 2000.0, 3000.0, 4000.0, 5000.0, 50000.0};
    const AxisSpec thnAxisInvMass{configThnAxisInvMass, "#it{M} (GeV/#it{c}^{2})"};
    const AxisSpec thnAxisPt{configThnAxisPt, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec thnAxisCosThetaStar{configThnAxisCosThetaStar, "cos(#vartheta_{OP})"};
    const AxisSpec thnAxisCentrality{configThnAxisCentrality, "Centrality (%)"};
    const AxisSpec thnAxisV2{configThnAxisV2, "V2"};
    const AxisSpec thnAxisRapidity{configThnAxisRapidity, "Rapidity"};
    const AxisSpec thnAxisSA{configThnAxisSA, "SA"};
    AxisSpec cumulantAxis = {200, -1, 1, "phi"};
    AxisSpec itsAxis = {8, -0.5, 7.5, "its"};
    AxisSpec tpcAxis = {130, 69.5, 199.5, "its"};
    AxisSpec squareAxis = {200, 0, 1, "aossquare"};
    AxisSpec phiAxis = {500, -6.28, 6.28, "phi"};
    AxisSpec resAxis = {6000, -30, 30, "Res"};
    AxisSpec resAxisSquare = {800, -1, 1, "Res"};
    AxisSpec centAxis = {8, 0, 80, "V0M (%)"};
    AxisSpec occupancyAxis = {occupancyBinning, "Occupancy"};

    histos.add("hTPCglobalmomcorr", "Momentum correlation", kTH3F, {{200, -10.0f, 10.0f}, {200, -10.0f, 10.0f}, {8, 0.0f, 80.0f}});
    histos.add("hpTvsRapidity", "pT vs Rapidity", kTH2F, {{100, 0.0f, 10.0f}, {300, -1.5f, 1.5f}});
    histos.add("hFTOCvsTPCNoCut", "Mult correlation FT0C vs. TPC without any cut", kTH2F, {{80, 0.0f, 80.0f}, {100, -0.5f, 5999.5f}});
    histos.add("hFTOCvsTPC", "Mult correlation FT0C vs. TPC", kTH2F, {{80, 0.0f, 80.0f}, {100, -0.5f, 5999.5f}});
    histos.add("hFTOCvsTPCSelected", "Mult correlation FT0C vs. TPC after selection", kTH2F, {{80, 0.0f, 80.0f}, {100, -0.5f, 5999.5f}});
    histos.add("hCentrality", "Centrality distribution", kTH1F, {{200, 0.0, 200.0}});
    histos.add("hVtxZ", "Vertex distribution in Z;Z (cm)", kTH1F, {{400, -20.0, 20.0}});
    histos.add("hEta", "Eta distribution", kTH1F, {{200, -1.0f, 1.0f}});
    histos.add("hDcaxy", "Dcaxy distribution", kTH1F, {{200, -1.0f, 1.0f}});
    histos.add("hDcaz", "Dcaz distribution", kTH1F, {{200, -1.0f, 1.0f}});

    histos.add("hITS", "ITS cluster", kTH2F, {{10, -0.5f, 9.5f}, {200, -1.0, 1.0}});
    histos.add("hTPC", "TPC crossed rows", kTH2F, {{90, 69.5f, 159.5f}, {200, -1.0, 1.0}});
    histos.add("hTPCScls", "TPC Shared cluster", kTH2F, {{16, -0.5f, 159.5f}, {200, -1.0, 1.0}});
    histos.add("hTPCSclsFrac", "Fraction of TPC Shared cluster", kTH2F, {{100, -0.0f, 2.0f}, {200, -1.0, 1.0}});

    histos.add("hNsigmaKaonTPC", "NsigmaKaon TPC distribution", kTH1F, {{200, -10.0f, 10.0f}});
    histos.add("hNsigmaKaonTOF", "NsigmaKaon TOF distribution", kTH1F, {{200, -10.0f, 10.0f}});
    histos.add("hPsiFT0C", "PsiFT0C", kTH3F, {centAxis, occupancyAxis, phiAxis});
    histos.add("hPsiFT0A", "PsiFT0A", kTH3F, {centAxis, occupancyAxis, phiAxis});
    histos.add("hPsiTPC", "PsiTPC", kTH3F, {centAxis, occupancyAxis, phiAxis});
    histos.add("hPsiTPCR", "PsiTPCR", kTH3F, {centAxis, occupancyAxis, phiAxis});
    histos.add("hPsiTPCL", "PsiTPCL", kTH3F, {centAxis, occupancyAxis, phiAxis});

    histos.add("hSparseV2SameEventCosPhi", "hSparseV2SameEventCosPhi", HistType::kTHnSparseD, {thnAxisInvMass, thnAxisPt, cumulantAxis, thnAxisCentrality});
    histos.add("hSparseV2SameEventSinPhi", "hSparseV2SameEventSinPhi", HistType::kTHnSparseD, {thnAxisInvMass, thnAxisPt, cumulantAxis, thnAxisCentrality});
    histos.add("hSparseV2SameEventCosPsi", "hSparseV2SameEventCosPsi", HistType::kTHnSparseD, {thnAxisInvMass, thnAxisPt, cumulantAxis, thnAxisCentrality});
    histos.add("hSparseV2SameEventSinPsi", "hSparseV2SameEventSinPsi", HistType::kTHnSparseD, {thnAxisInvMass, thnAxisPt, cumulantAxis, thnAxisCentrality});

    histos.add("hSparseV2SameEventCosDeltaPhi", "hSparseV2SameEventCosDeltaPhi", HistType::kTHnSparseD, {thnAxisInvMass, thnAxisPt, thnAxisV2, thnAxisCentrality});
    histos.add("hSparseV2MixedEventCosDeltaPhi", "hSparseV2MixedEventCosDeltaPhi", HistType::kTHnSparseD, {thnAxisInvMass, thnAxisPt, thnAxisV2, thnAxisCentrality});

    histos.add("hSparseV2SameEventCos2DeltaPhi", "hSparseV2SameEventCos2DeltaPhi", HistType::kTHnSparseD, {thnAxisInvMass, thnAxisPt, thnAxisV2, thnAxisCentrality});
    histos.add("hSparseV2MixedEventCos2DeltaPhi", "hSparseV2MixedEventCos2DeltaPhi", HistType::kTHnSparseD, {thnAxisInvMass, thnAxisPt, thnAxisV2, thnAxisCentrality});

    histos.add("hSparseV2SameEventCosDeltaPhiSquare", "hSparseV2SameEventCosDeltaPhiSquare", HistType::kTHnSparseD, {thnAxisInvMass, thnAxisPt, squareAxis, thnAxisCentrality});
    histos.add("hSparseV2SameEventCosDeltaPhiCube", "hSparseV2SameEventCosDeltaPhiCube", HistType::kTHnSparseD, {thnAxisInvMass, thnAxisPt, thnAxisV2, thnAxisCentrality});
    histos.add("hSparseV2MixedEventCosDeltaPhiSquare", "hSparseV2MixedEventCosDeltaPhiSquare", HistType::kTHnSparseD, {thnAxisInvMass, thnAxisPt, squareAxis, thnAxisCentrality});

    histos.add("hSparseV2SameEventSinDeltaPhi", "hSparseV2SameEventSinDeltaPhi", HistType::kTHnSparseD, {thnAxisInvMass, thnAxisPt, thnAxisV2, thnAxisCentrality});
    histos.add("hSparseV2MixedEventSinDeltaPhi", "hSparseV2MixedEventSinDeltaPhi", HistType::kTHnSparseD, {thnAxisInvMass, thnAxisPt, thnAxisV2, thnAxisCentrality});

    if (fillSA) {
      histos.add("hSparseV2SameEventSA", "hSparseV2SameEventSA", HistType::kTHnSparseD, {thnAxisInvMass, thnAxisPt, thnAxisSA, thnAxisRapidity, thnAxisCentrality});
      histos.add("hSparseV2MixedEventSA", "hSparseV2MixedEventSA", HistType::kTHnSparseD, {thnAxisInvMass, thnAxisPt, thnAxisSA, thnAxisRapidity, thnAxisCentrality});
      histos.add("hSparseV2SameEventCosThetaStar", "hSparseV2SameEventCosThetaStar", HistType::kTHnSparseD, {thnAxisInvMass, thnAxisPt, thnAxisCosThetaStar, thnAxisRapidity, thnAxisCentrality});
      histos.add("hSparseV2MixedEventCosThetaStar", "hSparseV2MixedEventCosThetaStar", HistType::kTHnSparseD, {thnAxisInvMass, thnAxisPt, thnAxisCosThetaStar, thnAxisRapidity, thnAxisCentrality});
    }
    // histogram for resolution
    histos.add("ResFT0CTPC", "ResFT0CTPC", kTH3F, {centAxis, occupancyAxis, resAxis});
    histos.add("ResFT0CTPCR", "ResFT0CTPCR", kTH3F, {centAxis, occupancyAxis, resAxis});
    histos.add("ResFT0CTPCL", "ResFT0CTPCL", kTH3F, {centAxis, occupancyAxis, resAxis});
    histos.add("ResTPCRTPCL", "ResTPCRTPCL", kTH3F, {centAxis, occupancyAxis, resAxis});
    histos.add("ResFT0CFT0A", "ResFT0CFT0A", kTH3F, {centAxis, occupancyAxis, resAxis});
    histos.add("ResFT0ATPC", "ResFT0ATPC", kTH3F, {centAxis, occupancyAxis, resAxis});

    histos.add("Res4FT0CTPC", "Res4FT0CTPC", kTH3F, {centAxis, occupancyAxis, resAxis});
    histos.add("Res4FT0CTPCR", "Res4FT0CTPCR", kTH3F, {centAxis, occupancyAxis, resAxis});
    histos.add("Res4FT0CTPCL", "Res4FT0CTPCL", kTH3F, {centAxis, occupancyAxis, resAxis});
    histos.add("Res4TPCRTPCL", "Res4TPCRTPCL", kTH3F, {centAxis, occupancyAxis, resAxis});
    histos.add("Res4FT0CFT0A", "Res4FT0CFT0A", kTH3F, {centAxis, occupancyAxis, resAxis});
    histos.add("Res4FT0ATPC", "Res4FT0ATPC", kTH3F, {centAxis, occupancyAxis, resAxis});

    histos.add("ResFT0CTPCSquare", "ResFT0CTPCSquare", kTH3F, {centAxis, occupancyAxis, resAxisSquare});
    histos.add("ResFT0CTPCRSquare", "ResFT0CTPCRSquare", kTH3F, {centAxis, occupancyAxis, resAxisSquare});
    histos.add("ResFT0CTPCLSquare", "ResFT0CTPCLSquare", kTH3F, {centAxis, occupancyAxis, resAxisSquare});
    histos.add("ResTPCRTPCLSquare", "ResTPCRTPCLSquare", kTH3F, {centAxis, occupancyAxis, resAxisSquare});
    histos.add("ResFT0CFT0ASquare", "ResFT0CFT0ASquare", kTH3F, {centAxis, occupancyAxis, resAxisSquare});
    histos.add("ResFT0ATPCSquare", "ResFT0ATPCSquare", kTH3F, {centAxis, occupancyAxis, resAxisSquare});

    histos.add("ResSPFT0CTPC", "ResSPFT0CTPC", kTH3F, {centAxis, occupancyAxis, resAxis});
    histos.add("ResSPFT0CTPCR", "ResSPFT0CTPCR", kTH3F, {centAxis, occupancyAxis, resAxis});
    histos.add("ResSPFT0CTPCL", "ResSPFT0CTPCL", kTH3F, {centAxis, occupancyAxis, resAxis});
    histos.add("ResSPTPCRTPCL", "ResSPTPCRTPCL", kTH3F, {centAxis, occupancyAxis, resAxis});
    histos.add("ResSPFT0CFT0A", "ResSPFT0CFT0A", kTH3F, {centAxis, occupancyAxis, resAxis});
    histos.add("ResSPFT0ATPC", "ResSPFT0ATPC", kTH3F, {centAxis, occupancyAxis, resAxis});

    histos.add("Res4SPFT0CTPC", "Res4SPFT0CTPC", kTH3F, {centAxis, occupancyAxis, resAxis});
    histos.add("Res4SPFT0CTPCR", "Res4SPFT0CTPCR", kTH3F, {centAxis, occupancyAxis, resAxis});
    histos.add("Res4SPFT0CTPCL", "Res4SPFT0CTPCL", kTH3F, {centAxis, occupancyAxis, resAxis});
    histos.add("Res4SPTPCRTPCL", "Res4SPTPCRTPCL", kTH3F, {centAxis, occupancyAxis, resAxis});
    histos.add("Res4SPFT0CFT0A", "Res4SPFT0CFT0A", kTH3F, {centAxis, occupancyAxis, resAxis});
    histos.add("Res4SPFT0ATPC", "Res4SPFT0ATPC", kTH3F, {centAxis, occupancyAxis, resAxis});

    histos.add("ResTrackSPFT0CTPC", "ResTrackSPFT0CTPC", kTH3F, {centAxis, occupancyAxis, resAxis});
    histos.add("ResTrackSPFT0CTPCR", "ResTrackSPFT0CTPCR", kTH3F, {centAxis, occupancyAxis, resAxis});
    histos.add("ResTrackSPFT0CTPCL", "ResTrackSPFT0CTPCL", kTH3F, {centAxis, occupancyAxis, resAxis});
    histos.add("ResTrackSPTPCRTPCL", "ResTrackSPTPCRTPCL", kTH3F, {centAxis, occupancyAxis, resAxis});
    histos.add("ResTrackSPFT0CFT0A", "ResTrackSPFT0CFT0A", kTH3F, {centAxis, occupancyAxis, resAxis});
    histos.add("ResTrackSPFT0ATPC", "ResTrackSPFT0ATPC", kTH3F, {centAxis, occupancyAxis, resAxis});

    // MC histogram
    if (isMC) {
      histos.add("hMC", "MC Event statistics", kTH1F, {{10, 0.0f, 10.0f}});
      histos.add("h1PhiRecsplit", "Phi meson Rec split", kTH1F, {{100, 0.0f, 10.0f}});
      histos.add("CentPercentileMCRecHist", "MC Centrality", kTH1F, {{100, 0.0f, 100.0f}});

      histos.add("hSparseV2MCGenSA", "hSparseV2SameEventSA", HistType::kTHnSparseD, {thnAxisInvMass, thnAxisPt, thnAxisSA, thnAxisRapidity, thnAxisCentrality});
      histos.add("hSparseV2MCGenCosThetaStar_effy", "hSparseV2SameEventCosThetaStar_effy", HistType::kTHnSparseD, {thnAxisInvMass, thnAxisPt, thnAxisCosThetaStar, thnAxisRapidity, thnAxisCentrality});

      histos.add("hSparseV2MCRecSA", "hSparseV2SameEventSA", HistType::kTHnSparseD, {thnAxisInvMass, thnAxisPt, thnAxisSA, thnAxisRapidity, thnAxisCentrality});
      histos.add("hSparseV2MCRecCosThetaStar_effy", "hSparseV2SameEventCosThetaStar_effy", HistType::kTHnSparseD, {thnAxisInvMass, thnAxisPt, thnAxisCosThetaStar, thnAxisRapidity, thnAxisCentrality});

      // weight

      histos.add("hSparsePhiMCGenWeight", "hSparsePhiMCGenWeight", HistType::kTHnSparseD, {thnAxisCentrality, {36, 0.0f, TMath::Pi()}, {400, -1.0, 1}, thnAxisPt, {8, -0.8, 0.8}});
      histos.add("hSparsePhiMCRecWeight", "hSparsePhiMCRecWeight", HistType::kTHnSparseD, {thnAxisCentrality, {36, 0.0f, TMath::Pi()}, {400, -1.0, 1}, thnAxisPt, {8, -0.8, 0.8}});
      histos.add("hSparsePhiMCGenKaonWeight", "hSparsePhiMCGenKaonWeight", HistType::kTHnSparseD, {thnAxisCentrality, {36, 0.0f, TMath::Pi()}, {400, -1.0, 1}, axisPtKaonWeight, {8, -0.8, 0.8}});
      histos.add("hSparsePhiMCRecKaonWeight", "hSparsePhiMCRecKaonWeight", HistType::kTHnSparseD, {thnAxisCentrality, {36, 0.0f, TMath::Pi()}, {400, -1.0, 1}, axisPtKaonWeight, {8, -0.8, 0.8}});
      histos.add("hSparsePhiMCRecKaonMissMatchWeight", "hSparsePhiMCRecKaonMissMatchWeight", HistType::kTHnSparseD, {thnAxisCentrality, {36, 0.0f, TMath::Pi()}, {400, -1.0, 1}, axisPtKaonWeight, {8, -0.8, 0.8}});

      histos.add("hSparseMCGenWeight", "hSparseMCGenWeight", HistType::kTHnSparseD, {thnAxisCentrality, {36, 0.0f, TMath::Pi()}, {400, -1.0, 1}, axisPtKaonWeight, {8, -0.8, 0.8}});
      histos.add("hSparseMCRecWeight", "hSparseMCRecWeight", HistType::kTHnSparseD, {thnAxisCentrality, {36, 0.0f, TMath::Pi()}, {400, -1.0, 1}, axisPtKaonWeight, {8, -0.8, 0.8}});
      histos.add("hSparseMCRecAllTrackWeight", "hSparseMCRecAllTrackWeight", HistType::kTHnSparseD, {thnAxisCentrality, {36, 0.0, TMath::Pi()}, {400, -1.0, 1}, axisPtKaonWeight, {8, -0.8, 0.8}});

      histos.add("hImpactParameter", "Impact parameter", kTH1F, {{200, 0.0f, 20.0f}});
      histos.add("hEventPlaneAngle", "hEventPlaneAngle", kTH1F, {{200, -2.0f * TMath::Pi(), 2.0f * TMath::Pi()}});
      histos.add("hEventPlaneAngleRec", "hEventPlaneAngleRec", kTH1F, {{200, -2.0f * TMath::Pi(), 2.0f * TMath::Pi()}});
      histos.add("hNchVsImpactParameter", "hNchVsImpactParameter", kTH2F, {{200, 0.0f, 20.0f}, {500, -0.5f, 5000.5f}});
      histos.add("hSparseMCGenV2", "hSparseMCGenV2", HistType::kTHnSparseD, {thnAxisCentrality, {200, -1.0, 1.0}, axisPtKaonWeight});
      histos.add("hSparseMCRecV2", "hSparseMCRecV2", HistType::kTHnSparseD, {thnAxisCentrality, {200, -1.0, 1.0}, axisPtKaonWeight});
      // histos.add("hSparseMCGenV2Square", "hSparseMCGenV2Square", HistType::kTHnSparseD, {thnAxisCentrality, {1000, 0.0, 1.0}, axisPtKaonWeight});
      // histos.add("hSparseMCRecV2Square", "hSparseMCRecV2Square", HistType::kTHnSparseD, {thnAxisCentrality, {1000, 0.0, 1.0}, axisPtKaonWeight});
      histos.add("hSparseMCGenV2Square", "hSparseMCGenV2Square", HistType::kTH3D, {thnAxisCentrality, {1000, 0.0, 1.0}, axisPtKaonWeight});
      histos.add("hSparseMCRecV2Square", "hSparseMCRecV2Square", HistType::kTH3D, {thnAxisCentrality, {1000, 0.0, 1.0}, axisPtKaonWeight});
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
    if (useGlobalTrack && !(candidate.isGlobalTrack() && candidate.isPVContributor() && candidate.itsNCls() > cfgITScluster && candidate.tpcNClsCrossedRows() > cfgTPCcluster && candidate.tpcFractionSharedCls() < cfgTPCSharedcluster)) {
      return false;
    }
    if (!useGlobalTrack && !(candidate.tpcNClsFound() > cfgTPCcluster)) {
      return false;
    }
    return true;
  }
  template <typename T>
  bool selectionPIDpTdependent(const T& candidate)
  {
    if (candidate.pt() < 0.5 && TMath::Abs(candidate.tpcNSigmaKa()) < nsigmaCutTPC) {
      return true;
    }
    if (candidate.pt() >= 0.5 && candidate.hasTOF() && candidate.beta() > cfgCutTOFBeta && TMath::Abs(candidate.tpcNSigmaKa()) < nsigmaCutTPC && TMath::Abs(candidate.tofNSigmaKa()) < nsigmaCutTOF) {
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
    if (!candidate.hasTOF() && TMath::Abs(candidate.tpcNSigmaKa()) < nsigmaCutTPC) {
      return true;
    }
    if (candidate.hasTOF() && candidate.beta() > cfgCutTOFBeta && TMath::Abs(candidate.tpcNSigmaKa()) < nsigmaCutTPC && TMath::Abs(candidate.tofNSigmaKa()) < nsigmaCutTOF) {
      return true;
    }
    return false;
  }

  template <typename T>
  bool selectionPID2(const T& candidate)
  {
    if (candidate.hasTOF() && candidate.beta() > cfgCutTOFBeta && TMath::Abs(candidate.tofNSigmaKa()) < nsigmaCutTOF) {
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

  template <typename T>
  bool isFakeKaon(T const& track)
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

  using BinningTypeVertexContributor = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0C, o2::aod::evsel::NumTracksInTimeRange>;
  ROOT::Math::PxPyPzMVector PhiMesonMother, KaonPlus, KaonMinus, fourVecDauCM;
  ROOT::Math::XYZVector threeVecDauCM, threeVecDauCMXY, eventplaneVec, eventplaneVecNorm, beamvector;
  int currentRunNumber = -999;
  int lastRunNumber = -999;
  // TH3D* hweight;
  TH2D* hweight;
  void processSameEvent(EventCandidates::iterator const& collision, TrackCandidates const& /*tracks, aod::BCs const&*/, aod::BCsWithTimestamps const&)
  {
    // if (!collision.sel8() || !collision.triggereventep() || !collision.selection_bit(aod::evsel::kNoTimeFrameBorder) || !collision.selection_bit(aod::evsel::kNoITSROFrameBorder) || !collision.selection_bit(aod::evsel::kNoSameBunchPileup) || !collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV) || !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
    // return;
    // }
    if (!collision.sel8() || !collision.triggereventep() || !collision.selection_bit(aod::evsel::kNoSameBunchPileup)) {
      return;
    }
    auto centrality = collision.centFT0C();
    auto multTPC = collision.multNTracksPV();
    histos.fill(HIST("hFTOCvsTPCNoCut"), centrality, multTPC);
    auto posThisColl = posTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto negThisColl = negTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto psiFT0C = collision.psiFT0C();
    auto psiFT0A = collision.psiFT0A();
    auto psiTPC = collision.psiTPC();
    auto psiTPCR = collision.psiTPCR();
    auto psiTPCL = collision.psiTPCL();
    auto QFT0C = collision.qFT0C();
    auto QFT0A = collision.qFT0A();
    auto QTPC = collision.qTPC();
    auto QTPCR = collision.qTPCR();
    auto QTPCL = collision.qTPCL();
    int occupancy = collision.trackOccupancyInTimeRange();
    o2::aod::ITSResponse itsResponse;
    if (occupancy > cfgCutOccupancy) {
      return;
    }
    histos.fill(HIST("hFTOCvsTPC"), centrality, multTPC);
    if (additionalEvsel && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
      return;
    }
    if (additionalEvselITS && !collision.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll)) {
      return;
    }
    histos.fill(HIST("hFTOCvsTPCSelected"), centrality, multTPC);
    histos.fill(HIST("hPsiFT0C"), centrality, occupancy, psiFT0C);
    histos.fill(HIST("hPsiFT0A"), centrality, occupancy, psiFT0A);
    histos.fill(HIST("hPsiTPC"), centrality, occupancy, psiTPC);
    histos.fill(HIST("hPsiTPCR"), centrality, occupancy, psiTPCR);
    histos.fill(HIST("hPsiTPCL"), centrality, occupancy, psiTPCL);

    histos.fill(HIST("ResFT0CTPC"), centrality, occupancy, TMath::Cos(2.0 * (psiFT0C - psiTPC)));
    histos.fill(HIST("ResFT0CTPCR"), centrality, occupancy, TMath::Cos(2.0 * (psiFT0C - psiTPCR)));
    histos.fill(HIST("ResFT0CTPCL"), centrality, occupancy, TMath::Cos(2.0 * (psiFT0C - psiTPCL)));
    histos.fill(HIST("ResTPCRTPCL"), centrality, occupancy, TMath::Cos(2.0 * (psiTPCR - psiTPCL)));
    histos.fill(HIST("ResFT0CFT0A"), centrality, occupancy, TMath::Cos(2.0 * (psiFT0C - psiFT0A)));
    histos.fill(HIST("ResFT0ATPC"), centrality, occupancy, TMath::Cos(2.0 * (psiTPC - psiFT0A)));

    histos.fill(HIST("Res4FT0CTPC"), centrality, occupancy, TMath::Cos(4.0 * (psiFT0C - psiTPC)));
    histos.fill(HIST("Res4FT0CTPCR"), centrality, occupancy, TMath::Cos(4.0 * (psiFT0C - psiTPCR)));
    histos.fill(HIST("Res4FT0CTPCL"), centrality, occupancy, TMath::Cos(4.0 * (psiFT0C - psiTPCL)));
    histos.fill(HIST("Res4TPCRTPCL"), centrality, occupancy, TMath::Cos(4.0 * (psiTPCR - psiTPCL)));
    histos.fill(HIST("Res4FT0CFT0A"), centrality, occupancy, TMath::Cos(4.0 * (psiFT0C - psiFT0A)));
    histos.fill(HIST("Res4FT0ATPC"), centrality, occupancy, TMath::Cos(4.0 * (psiTPC - psiFT0A)));

    histos.fill(HIST("ResFT0CTPCSquare"), centrality, occupancy, TMath::Power(TMath::Cos(2.0 * (psiFT0C - psiTPC)), 2.0));
    histos.fill(HIST("ResFT0CTPCRSquare"), centrality, occupancy, TMath::Power(TMath::Cos(2.0 * (psiFT0C - psiTPCR)), 2.0));
    histos.fill(HIST("ResFT0CTPCLSquare"), centrality, occupancy, TMath::Power(TMath::Cos(2.0 * (psiFT0C - psiTPCL)), 2.0));
    histos.fill(HIST("ResTPCRTPCLSquare"), centrality, occupancy, TMath::Power(TMath::Cos(2.0 * (psiTPCR - psiTPCL)), 2.0));
    histos.fill(HIST("ResFT0CFT0ASquare"), centrality, occupancy, TMath::Power(TMath::Cos(2.0 * (psiFT0C - psiFT0A)), 2.0));
    histos.fill(HIST("ResFT0ATPCSquare"), centrality, occupancy, TMath::Power(TMath::Cos(2.0 * (psiTPC - psiFT0A)), 2.0));

    histos.fill(HIST("ResSPFT0CTPC"), centrality, occupancy, QFT0C * QTPC * TMath::Cos(2.0 * (psiFT0C - psiTPC)));
    histos.fill(HIST("ResSPFT0CTPCR"), centrality, occupancy, QFT0C * QTPCR * TMath::Cos(2.0 * (psiFT0C - psiTPCR)));
    histos.fill(HIST("ResSPFT0CTPCL"), centrality, occupancy, QFT0C * QTPCL * TMath::Cos(2.0 * (psiFT0C - psiTPCL)));
    histos.fill(HIST("ResSPTPCRTPCL"), centrality, occupancy, QTPCR * QTPCL * TMath::Cos(2.0 * (psiTPCR - psiTPCL)));
    histos.fill(HIST("ResSPFT0CFT0A"), centrality, occupancy, QFT0C * QFT0A * TMath::Cos(2.0 * (psiFT0C - psiFT0A)));
    histos.fill(HIST("ResSPFT0ATPC"), centrality, occupancy, QTPC * QFT0A * TMath::Cos(2.0 * (psiTPC - psiFT0A)));

    histos.fill(HIST("Res4SPFT0CTPC"), centrality, occupancy, QFT0C * QTPC * TMath::Cos(4.0 * (psiFT0C - psiTPC)));
    histos.fill(HIST("Res4SPFT0CTPCR"), centrality, occupancy, QFT0C * QTPCR * TMath::Cos(4.0 * (psiFT0C - psiTPCR)));
    histos.fill(HIST("Res4SPFT0CTPCL"), centrality, occupancy, QFT0C * QTPCL * TMath::Cos(4.0 * (psiFT0C - psiTPCL)));
    histos.fill(HIST("Res4SPTPCRTPCL"), centrality, occupancy, QTPCR * QTPCL * TMath::Cos(4.0 * (psiTPCR - psiTPCL)));
    histos.fill(HIST("Res4SPFT0CFT0A"), centrality, occupancy, QFT0C * QFT0A * TMath::Cos(4.0 * (psiFT0C - psiFT0A)));
    histos.fill(HIST("Res4SPFT0ATPC"), centrality, occupancy, QTPC * QFT0A * TMath::Cos(4.0 * (psiTPC - psiFT0A)));

    histos.fill(HIST("hCentrality"), centrality);
    histos.fill(HIST("hVtxZ"), collision.posZ());

    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    currentRunNumber = collision.bc_as<aod::BCsWithTimestamps>().runNumber();
    if (useWeight && (currentRunNumber != lastRunNumber)) {
      // hweight = ccdb->getForTimeStamp<TH3D>(ConfWeightPath.value, bc.timestamp());
      hweight = ccdb->getForTimeStamp<TH2D>(ConfWeightPath.value, bc.timestamp());
    }
    lastRunNumber = currentRunNumber;
    int Npostrack = 0;
    float weight1 = 1.0;
    float weight2 = 1.0;
    for (auto track1 : posThisColl) {
      // track selection
      if (!selectionTrack(track1)) {
        continue;
      }
      // PID check
      if (ispTdepPID && !isTOFOnly && !selectionPIDpTdependent(track1)) {
        continue;
      }
      if (!ispTdepPID && !isTOFOnly && !selectionPID(track1)) {
        continue;
      }
      if (isTOFOnly && !selectionPID2(track1)) {
        continue;
      }
      if (useGlobalTrack && track1.p() < 1.0 && !(itsResponse.nSigmaITS<o2::track::PID::Kaon>(track1) > -2.5 && itsResponse.nSigmaITS<o2::track::PID::Kaon>(track1) < 2.5)) {
        continue;
      }
      histos.fill(HIST("hTPCglobalmomcorr"), track1.p() / track1.sign(), track1.p() - track1.tpcInnerParam(), centrality);
      histos.fill(HIST("hEta"), track1.eta());
      histos.fill(HIST("hDcaxy"), track1.dcaXY());
      histos.fill(HIST("hDcaz"), track1.dcaZ());
      histos.fill(HIST("hNsigmaKaonTPC"), track1.tpcNSigmaKa());
      histos.fill(HIST("hNsigmaKaonTOF"), track1.tofNSigmaKa());
      auto V2Track = TMath::Cos(2.0 * GetPhiInRange(track1.phi() - psiFT0C));
      histos.fill(HIST("hITS"), track1.itsNCls(), V2Track);
      histos.fill(HIST("hTPC"), track1.tpcNClsCrossedRows(), V2Track);
      histos.fill(HIST("hTPCScls"), track1.tpcNClsShared(), V2Track);
      histos.fill(HIST("hTPCSclsFrac"), track1.tpcFractionSharedCls(), V2Track);
      auto track1ID = track1.globalIndex();
      if (useWeight) {
        if (track1.pt() < 10.0 && track1.pt() > 0.15) {
          // weight1 = hweight->GetBinContent(hweight->FindBin(centrality, GetPhiInRange(track1.phi() - psiFT0C), track1.pt() + 0.000005));
          weight1 = 1 + hweight->GetBinContent(hweight->FindBin(centrality, track1.pt() + 0.000005)) * TMath::Cos(2.0 * GetPhiInRange(track1.phi() - psiFT0C));
        } else {
          weight1 = 1;
        }
      }
      for (auto track2 : negThisColl) {
        // track selection
        if (!selectionTrack(track2)) {
          continue;
        }
        // PID check
        if (ispTdepPID && !isTOFOnly && !selectionPIDpTdependent(track2)) {
          continue;
        }
        if (!ispTdepPID && !isTOFOnly && !selectionPID(track2)) {
          continue;
        }
        if (isTOFOnly && !selectionPID2(track2)) {
          continue;
        }
        auto track2ID = track2.globalIndex();
        if (track2ID == track1ID) {
          continue;
        }
        if (!selectionPair(track1, track2)) {
          continue;
        }
        if (Npostrack == 0) {
          histos.fill(HIST("hTPCglobalmomcorr"), track2.p() / track2.sign(), track2.p() - track2.tpcInnerParam(), centrality);
        }
        if (removefaketrak && isFakeKaon(track1)) {
          continue;
        }
        if (removefaketrak && isFakeKaon(track2)) {
          continue;
        }
        if (useGlobalTrack && track2.p() < 1.0 && !(itsResponse.nSigmaITS<o2::track::PID::Kaon>(track2) > -2.5 && itsResponse.nSigmaITS<o2::track::PID::Kaon>(track2) < 2.5)) {
          continue;
        }
        if (useWeight) {
          if (track2.pt() < 10.0 && track2.pt() > 0.15) {
            // weight2 = hweight->GetBinContent(hweight->FindBin(centrality, GetPhiInRange(track2.phi() - psiFT0C), track2.pt() + 0.000005));
            weight2 = 1 + hweight->GetBinContent(hweight->FindBin(centrality, track2.pt() + 0.000005)) * TMath::Cos(2.0 * GetPhiInRange(track2.phi() - psiFT0C));
          } else {
            weight2 = 1;
          }
        }
        KaonPlus = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), massKa);
        KaonMinus = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), massKa);
        PhiMesonMother = KaonPlus + KaonMinus;
        auto phiminuspsi = GetPhiInRange(PhiMesonMother.Phi() - psiFT0C);
        auto v2 = TMath::Cos(2.0 * phiminuspsi);
        auto v2acc = TMath::Cos(4.0 * phiminuspsi);
        auto v2sin = TMath::Sin(2.0 * phiminuspsi);
        auto phimother = PhiMesonMother.Phi();
        histos.fill(HIST("hpTvsRapidity"), PhiMesonMother.Pt(), PhiMesonMother.Rapidity());
        auto totalweight = weight1 * weight2;
        if (totalweight <= 0.0000005) {
          totalweight = 1.0;
        }
        // LOGF(info, Form("weight %f    %f",weight1, weight2));
        if (TMath::Abs(PhiMesonMother.Rapidity()) < confRapidity) {
          histos.fill(HIST("ResTrackSPFT0CTPC"), centrality, occupancy, QFT0C * QTPC * TMath::Cos(2.0 * (psiFT0C - psiTPC)));
          histos.fill(HIST("ResTrackSPFT0CTPCR"), centrality, occupancy, QFT0C * QTPCR * TMath::Cos(2.0 * (psiFT0C - psiTPCR)));
          histos.fill(HIST("ResTrackSPFT0CTPCL"), centrality, occupancy, QFT0C * QTPCL * TMath::Cos(2.0 * (psiFT0C - psiTPCL)));
          histos.fill(HIST("ResTrackSPTPCRTPCL"), centrality, occupancy, QTPCR * QTPCL * TMath::Cos(2.0 * (psiTPCR - psiTPCL)));
          histos.fill(HIST("ResTrackSPFT0CFT0A"), centrality, occupancy, QFT0C * QFT0A * TMath::Cos(2.0 * (psiFT0C - psiFT0A)));
          histos.fill(HIST("ResTrackSPFT0ATPC"), centrality, occupancy, QTPC * QFT0A * TMath::Cos(2.0 * (psiTPC - psiFT0A)));
          if (useWeight) {
            histos.fill(HIST("hSparseV2SameEventCosDeltaPhi"), PhiMesonMother.M(), PhiMesonMother.Pt(), v2 * QFT0C, centrality, 1 / totalweight);
            histos.fill(HIST("hSparseV2SameEventCos2DeltaPhi"), PhiMesonMother.M(), PhiMesonMother.Pt(), v2acc * QFT0C, centrality, 1 / totalweight);

            histos.fill(HIST("hSparseV2SameEventCosDeltaPhiSquare"), PhiMesonMother.M(), PhiMesonMother.Pt(), v2 * v2, centrality, 1 / totalweight);
            histos.fill(HIST("hSparseV2SameEventSinDeltaPhi"), PhiMesonMother.M(), PhiMesonMother.Pt(), v2sin * QFT0C, centrality, 1 / totalweight);

            histos.fill(HIST("hSparseV2SameEventCosPhi"), PhiMesonMother.M(), PhiMesonMother.Pt(), TMath::Cos(2.0 * phimother), centrality, 1 / totalweight);
            histos.fill(HIST("hSparseV2SameEventSinPhi"), PhiMesonMother.M(), PhiMesonMother.Pt(), TMath::Sin(2.0 * phimother), centrality, 1 / totalweight);
            histos.fill(HIST("hSparseV2SameEventCosPsi"), PhiMesonMother.M(), PhiMesonMother.Pt(), TMath::Cos(2.0 * psiFT0C), centrality, 1 / totalweight);
            histos.fill(HIST("hSparseV2SameEventSinPsi"), PhiMesonMother.M(), PhiMesonMother.Pt(), TMath::Sin(2.0 * psiFT0C), centrality, 1 / totalweight);
          } else {
            if (useSP) {
              histos.fill(HIST("hSparseV2SameEventCosDeltaPhi"), PhiMesonMother.M(), PhiMesonMother.Pt(), v2 * QFT0C, centrality);
              histos.fill(HIST("hSparseV2SameEventCos2DeltaPhi"), PhiMesonMother.M(), PhiMesonMother.Pt(), v2acc * QFT0C, centrality);
            } else {
              histos.fill(HIST("hSparseV2SameEventCosDeltaPhi"), PhiMesonMother.M(), PhiMesonMother.Pt(), v2, centrality);
              histos.fill(HIST("hSparseV2SameEventCos2DeltaPhi"), PhiMesonMother.M(), PhiMesonMother.Pt(), v2acc, centrality);
            }
            histos.fill(HIST("hSparseV2SameEventCosDeltaPhiSquare"), PhiMesonMother.M(), PhiMesonMother.Pt(), v2 * v2, centrality);
            histos.fill(HIST("hSparseV2SameEventCosDeltaPhiCube"), PhiMesonMother.M(), PhiMesonMother.Pt(), v2 * v2 * v2, centrality);
            histos.fill(HIST("hSparseV2SameEventSinDeltaPhi"), PhiMesonMother.M(), PhiMesonMother.Pt(), v2sin * QFT0C, centrality);

            histos.fill(HIST("hSparseV2SameEventCosPhi"), PhiMesonMother.M(), PhiMesonMother.Pt(), TMath::Cos(2.0 * phimother), centrality);
            histos.fill(HIST("hSparseV2SameEventSinPhi"), PhiMesonMother.M(), PhiMesonMother.Pt(), TMath::Sin(2.0 * phimother), centrality);
            histos.fill(HIST("hSparseV2SameEventCosPsi"), PhiMesonMother.M(), PhiMesonMother.Pt(), TMath::Cos(2.0 * psiFT0C), centrality);
            histos.fill(HIST("hSparseV2SameEventSinPsi"), PhiMesonMother.M(), PhiMesonMother.Pt(), TMath::Sin(2.0 * psiFT0C), centrality);
          }
        }
        if (fillSA) {
          ROOT::Math::Boost boost{PhiMesonMother.BoostToCM()};
          fourVecDauCM = boost(KaonMinus);
          threeVecDauCM = fourVecDauCM.Vect();
          threeVecDauCMXY = ROOT::Math::XYZVector(threeVecDauCM.X(), threeVecDauCM.Y(), 0.);
          eventplaneVec = ROOT::Math::XYZVector(std::cos(2.0 * psiFT0C), std::sin(2.0 * psiFT0C), 0);
          eventplaneVecNorm = ROOT::Math::XYZVector(std::sin(2.0 * psiFT0C), -std::cos(2.0 * psiFT0C), 0);
          auto cosPhistarminuspsi = GetPhiInRange(fourVecDauCM.Phi() - psiFT0C);
          auto SA = TMath::Cos(2.0 * cosPhistarminuspsi);
          auto cosThetaStar = eventplaneVecNorm.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2()) / std::sqrt(eventplaneVecNorm.Mag2());
          if (useWeight) {
            histos.fill(HIST("hSparseV2SameEventSA"), PhiMesonMother.M(), PhiMesonMother.Pt(), SA, TMath::Abs(PhiMesonMother.Rapidity()), centrality, 1 / totalweight);
            histos.fill(HIST("hSparseV2SameEventCosThetaStar"), PhiMesonMother.M(), PhiMesonMother.Pt(), cosThetaStar, TMath::Abs(PhiMesonMother.Rapidity()), centrality, 1 / totalweight);
          } else {
            histos.fill(HIST("hSparseV2SameEventSA"), PhiMesonMother.M(), PhiMesonMother.Pt(), SA, TMath::Abs(PhiMesonMother.Rapidity()), centrality);
            histos.fill(HIST("hSparseV2SameEventCosThetaStar"), PhiMesonMother.M(), PhiMesonMother.Pt(), cosThetaStar, TMath::Abs(PhiMesonMother.Rapidity()));
          }
        }
      }
      Npostrack = Npostrack + 1;
    }
  }
  PROCESS_SWITCH(phipbpb, processSameEvent, "Process Same event", true);
  void processMixedEventOpti(EventCandidates const& collisions, TrackCandidates const& tracks)
  {
    auto tracksTuple = std::make_tuple(tracks);
    BinningTypeVertexContributor binningOnPositions{{axisVertex, axisMultiplicityClass, axisOccup}, true};
    SameKindPair<EventCandidates, TrackCandidates, BinningTypeVertexContributor> pair{binningOnPositions, cfgNoMixedEvents, -1, collisions, tracksTuple, &cache};
    for (auto& [collision1, tracks1, collision2, tracks2] : pair) {
      // if (!collision1.sel8() || !collision1.triggereventep() || !collision1.selection_bit(aod::evsel::kNoTimeFrameBorder) || !collision1.selection_bit(aod::evsel::kNoITSROFrameBorder) || !collision1.selection_bit(aod::evsel::kNoSameBunchPileup) || !collision1.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV) || !collision1.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
      // continue;
      // }
      // if (!collision2.sel8() || !collision2.triggereventep() || !collision2.selection_bit(aod::evsel::kNoTimeFrameBorder) || !collision2.selection_bit(aod::evsel::kNoITSROFrameBorder) || !collision2.selection_bit(aod::evsel::kNoSameBunchPileup) || !collision2.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV) || !collision2.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
      // continue;
      // }
      if (!collision1.sel8() || !collision1.triggereventep() || !collision1.selection_bit(aod::evsel::kNoSameBunchPileup)) {
        continue;
      }
      if (!collision2.sel8() || !collision2.triggereventep() || !collision2.selection_bit(aod::evsel::kNoSameBunchPileup)) {
        continue;
      }
      if (collision1.bcId() == collision2.bcId()) {
        continue;
      }
      o2::aod::ITSResponse itsResponse;
      int occupancy1 = collision1.trackOccupancyInTimeRange();
      int occupancy2 = collision2.trackOccupancyInTimeRange();
      if (occupancy1 > cfgCutOccupancy) {
        continue;
      }
      if (occupancy2 > cfgCutOccupancy) {
        continue;
      }
      auto centrality = collision1.centFT0C();
      // auto centrality2 = collision2.centFT0C();
      auto psiFT0C = collision1.psiFT0C();
      auto QFT0C = collision1.qFT0C();
      if (additionalEvsel && !collision1.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
        continue;
      }
      if (additionalEvsel && !collision2.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
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
        if (useGlobalTrack && track1.p() < 1.0 && !(itsResponse.nSigmaITS<o2::track::PID::Kaon>(track1) > -2.5 && itsResponse.nSigmaITS<o2::track::PID::Kaon>(track1) < 2.5)) {
          continue;
        }
        if (useGlobalTrack && track2.p() < 1.0 && !(itsResponse.nSigmaITS<o2::track::PID::Kaon>(track2) > -2.5 && itsResponse.nSigmaITS<o2::track::PID::Kaon>(track2) < 2.5)) {
          continue;
        }
        if (!selectionTrack(track1) || !selectionTrack(track2)) {
          continue;
        }
        // PID check
        if (ispTdepPID && !isTOFOnly && (!selectionPIDpTdependent(track1) || !selectionPIDpTdependent(track2))) {
          continue;
        }
        if (!ispTdepPID && !isTOFOnly && (!selectionPID(track1) || !selectionPID(track2))) {
          continue;
        }
        if (isTOFOnly && (!selectionPID2(track1) || !selectionPID2(track2))) {
          continue;
        }
        if (!selectionPair(track1, track2)) {
          continue;
        }
        if (removefaketrak && isFakeKaon(track1)) {
          continue;
        }
        if (removefaketrak && isFakeKaon(track2)) {
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
        auto phiminuspsi = GetPhiInRange(PhiMesonMother.Phi() - psiFT0C);
        auto v2 = TMath::Cos(2.0 * phiminuspsi);
        auto v2acc = TMath::Cos(4.0 * phiminuspsi);
        auto v2sin = TMath::Sin(2.0 * phiminuspsi);
        histos.fill(HIST("hpTvsRapidity"), PhiMesonMother.Pt(), PhiMesonMother.Rapidity());
        if (TMath::Abs(PhiMesonMother.Rapidity()) < confRapidity) {
          if (useSP) {
            histos.fill(HIST("hSparseV2MixedEventCosDeltaPhi"), PhiMesonMother.M(), PhiMesonMother.Pt(), v2 * QFT0C, centrality);
            histos.fill(HIST("hSparseV2MixedEventCos2DeltaPhi"), PhiMesonMother.M(), PhiMesonMother.Pt(), v2acc * QFT0C, centrality);
          } else {
            histos.fill(HIST("hSparseV2MixedEventCosDeltaPhi"), PhiMesonMother.M(), PhiMesonMother.Pt(), v2, centrality);
            histos.fill(HIST("hSparseV2MixedEventCos2DeltaPhi"), PhiMesonMother.M(), PhiMesonMother.Pt(), v2acc, centrality);
          }
          histos.fill(HIST("hSparseV2MixedEventCosDeltaPhiSquare"), PhiMesonMother.M(), PhiMesonMother.Pt(), v2 * v2, centrality);
          histos.fill(HIST("hSparseV2MixedEventSinDeltaPhi"), PhiMesonMother.M(), PhiMesonMother.Pt(), v2sin * QFT0C, centrality);
        }
        if (fillSA) {
          ROOT::Math::Boost boost{PhiMesonMother.BoostToCM()};
          fourVecDauCM = boost(KaonMinus);
          threeVecDauCM = fourVecDauCM.Vect();
          threeVecDauCMXY = ROOT::Math::XYZVector(threeVecDauCM.X(), threeVecDauCM.Y(), 0.);
          eventplaneVec = ROOT::Math::XYZVector(std::cos(2.0 * psiFT0C), std::sin(2.0 * psiFT0C), 0);
          eventplaneVecNorm = ROOT::Math::XYZVector(std::sin(2.0 * psiFT0C), -std::cos(2.0 * psiFT0C), 0);
          auto cosPhistarminuspsi = GetPhiInRange(fourVecDauCM.Phi() - psiFT0C);
          auto SA = TMath::Cos(2.0 * cosPhistarminuspsi);
          auto cosThetaStar = eventplaneVecNorm.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2()) / std::sqrt(eventplaneVecNorm.Mag2());
          histos.fill(HIST("hSparseV2MixedEventSA"), PhiMesonMother.M(), PhiMesonMother.Pt(), SA, TMath::Abs(PhiMesonMother.Rapidity()), centrality);
          histos.fill(HIST("hSparseV2MixedEventCosThetaStar"), PhiMesonMother.M(), PhiMesonMother.Pt(), cosThetaStar, TMath::Abs(PhiMesonMother.Rapidity()), centrality);
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
          int track1PDG = TMath::Abs(mctrack1.pdgCode());
          int track2PDG = TMath::Abs(mctrack2.pdgCode());
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
              if (TMath::Abs(mothertrack1.y()) > confRapidity) {
                continue;
              }
              if (TMath::Abs(mothertrack1.pdgCode()) != 333) {
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
              ROOT::Math::Boost boost{PhiMesonMother.BoostToCM()};
              fourVecDauCM = boost(KaonMinus);
              threeVecDauCM = fourVecDauCM.Vect();
              beamvector = ROOT::Math::XYZVector(0.0, 0.0, 1.0);
              threeVecDauCMXY = ROOT::Math::XYZVector(threeVecDauCM.X(), threeVecDauCM.Y(), 0.);
              eventplaneVec = ROOT::Math::XYZVector(std::cos(2.0 * psiFT0C), std::sin(2.0 * psiFT0C), 0);
              eventplaneVecNorm = ROOT::Math::XYZVector(std::sin(2.0 * psiFT0C), -std::cos(2.0 * psiFT0C), 0);
              auto cosPhistarminuspsi = GetPhiInRange(fourVecDauCM.Phi() - psiFT0C);
              auto SA = TMath::Cos(2.0 * cosPhistarminuspsi);
              auto cosThetaStar = eventplaneVec.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2()) / std::sqrt(eventplaneVec.Mag2());
              histos.fill(HIST("hSparseV2MCRecCosThetaStar_effy"), PhiMesonMother.M(), PhiMesonMother.Pt(), cosThetaStar, TMath::Abs(PhiMesonMother.Rapidity()), centrality);
              histos.fill(HIST("hSparseV2MCRecSA"), PhiMesonMother.M(), PhiMesonMother.Pt(), SA, TMath::Abs(PhiMesonMother.Rapidity()), centrality);
            }
          }
        }
      }
      // loop over generated particle
      for (auto& mcParticle : GenParticles) {
        if (TMath::Abs(mcParticle.y()) > confRapidity) {
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
          auto cosThetaStar = eventplaneVec.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2()) / std::sqrt(eventplaneVec.Mag2());
          histos.fill(HIST("hSparseV2MCGenCosThetaStar_effy"), PhiMesonMother.M(), PhiMesonMother.Pt(), cosThetaStar, TMath::Abs(PhiMesonMother.Rapidity()), centrality);
          histos.fill(HIST("hSparseV2MCGenSA"), PhiMesonMother.M(), PhiMesonMother.Pt(), SA, TMath::Abs(PhiMesonMother.Rapidity()), centrality);
        }
      }
    } // rec collision loop

  } // process MC
  PROCESS_SWITCH(phipbpb, processMC, "Process MC", false);
  using recoTracks = soa::Join<aod::TracksIU, aod::TracksExtra>;
  void processMCweight(aod::McCollision const& mcCollision, soa::Join<aod::McParticles, aod::ParticlesToTracks> const& mcParticles, recoTracks const&)
  {
    float imp = mcCollision.impactParameter();
    float evPhi = mcCollision.eventPlaneAngle() / 2.0;
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
    // if (evPhi < 0)
    //   evPhi += 2. * TMath::Pi();

    int nCh = 0;

    if (centclass > 0 && centclass < 80) {
      // event within range
      histos.fill(HIST("hImpactParameter"), imp);
      histos.fill(HIST("hEventPlaneAngle"), evPhi);
      for (auto const& mcParticle : mcParticles) {

        float deltaPhi = mcParticle.phi() - mcCollision.eventPlaneAngle();
        // focus on bulk: e, mu, pi, k, p
        int pdgCode = TMath::Abs(mcParticle.pdgCode());
        if (checkAllCharge && pdgCode != 11 && pdgCode != 13 && pdgCode != 211 && pdgCode != 321 && pdgCode != 2212)
          continue;
        if (!checkAllCharge && pdgCode != 321)
          continue;
        if (!mcParticle.isPhysicalPrimary())
          continue;
        if (TMath::Abs(mcParticle.eta()) > 0.8) // main acceptance
          continue;
        histos.fill(HIST("hSparseMCGenWeight"), centclass, GetPhiInRange(deltaPhi), TMath::Power(TMath::Cos(4.0 * GetPhiInRange(deltaPhi)), 1.0), mcParticle.pt(), mcParticle.eta());
        histos.fill(HIST("hSparseMCGenV2"), centclass, TMath::Cos(2.0 * GetPhiInRange(deltaPhi)), mcParticle.pt());
        histos.fill(HIST("hSparseMCGenV2Square"), centclass, TMath::Power(TMath::Cos(2.0 * GetPhiInRange(deltaPhi)), 2.0), mcParticle.pt());
        nCh++;
        bool validGlobal = false;
        bool validAny = false;
        if (mcParticle.has_tracks()) {
          auto const& tracks = mcParticle.tracks_as<recoTracks>();
          for (auto const& track : tracks) {
            if (track.hasTPC() && track.hasITS()) {
              validGlobal = true;
            }
            if (track.hasTPC() || track.hasITS()) {
              validAny = true;
            }
          }
        }
        // if valid global, fill
        if (validGlobal) {
          histos.fill(HIST("hSparseMCRecWeight"), centclass, GetPhiInRange(deltaPhi), TMath::Power(TMath::Cos(4.0 * GetPhiInRange(deltaPhi)), 1.0), mcParticle.pt(), mcParticle.eta());
          histos.fill(HIST("hSparseMCRecV2"), centclass, TMath::Cos(2.0 * GetPhiInRange(deltaPhi)), mcParticle.pt());
          histos.fill(HIST("hSparseMCRecV2Square"), centclass, TMath::Power(TMath::Cos(2.0 * GetPhiInRange(deltaPhi)), 2.0), mcParticle.pt());
        }
        if (validAny) {
          histos.fill(HIST("hSparseMCRecAllTrackWeight"), centclass, GetPhiInRange(deltaPhi), TMath::Power(TMath::Cos(4.0 * GetPhiInRange(deltaPhi)), 1.0), mcParticle.pt(), mcParticle.eta());
          histos.fill(HIST("hEventPlaneAngleRec"), GetPhiInRange(deltaPhi));
        }
        // if any track present, fill
      }
    }
    histos.fill(HIST("hNchVsImpactParameter"), imp, nCh);
  }
  PROCESS_SWITCH(phipbpb, processMCweight, "Process MC Weight", false);

  void processMCPhiWeight(CollisionMCTrueTable::iterator const& TrueCollision, CollisionMCRecTableCentFT0C const& RecCollisions, TrackMCTrueTable const& GenParticles, FilTrackMCRecTable const& RecTracks)
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

        if (selectionTrack(track1) && selectionPIDpTdependent(track1) && TMath::Abs(mctrack1.pdgCode()) == 321 && mctrack1.isPhysicalPrimary()) {
          histos.fill(HIST("hSparsePhiMCRecKaonWeight"), centclass, GetPhiInRange(mctrack1.phi() - psiFT0C), TMath::Power(TMath::Cos(4.0 * GetPhiInRange(mctrack1.phi() - psiFT0C)), 1.0), mctrack1.pt(), mctrack1.eta());
        }

        if (selectionTrack(track1) && track1.pt() > 0.5 && track1.hasTOF() && TMath::Abs(track1.tofNSigmaKa()) > nsigmaCutTOF && TMath::Abs(track1.tpcNSigmaKa()) < nsigmaCutTPC && TMath::Abs(mctrack1.pdgCode()) == 321 && mctrack1.isPhysicalPrimary()) {
          histos.fill(HIST("hSparsePhiMCRecKaonMissMatchWeight"), centclass, GetPhiInRange(mctrack1.phi() - psiFT0C), TMath::Power(TMath::Cos(4.0 * GetPhiInRange(mctrack1.phi() - psiFT0C)), 1.0), mctrack1.pt(), mctrack1.eta());
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
          if (!(track1PDG == 321 && track2PDG == 321)) {
            continue;
          }
          if (!selectionTrack(track1) || !selectionTrack(track2) || track1.sign() * track2.sign() > 0) {
            continue;
          }
          // PID check
          if (ispTdepPID && !isTOFOnly && (!selectionPIDpTdependent(track1) || !selectionPIDpTdependent(track2))) {
            continue;
          }
          if (!ispTdepPID && !isTOFOnly && (!selectionPID(track1) || !selectionPID(track2))) {
            continue;
          }
          if (isTOFOnly && (!selectionPID2(track1) || !selectionPID2(track2))) {
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
              if (TMath::Abs(mothertrack1.pdgCode()) != 333) {
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
              histos.fill(HIST("hSparsePhiMCRecWeight"), centclass, PhiMinusPsi, TMath::Power(TMath::Cos(4.0 * PhiMinusPsi), 1.0), mothertrack1.pt(), mothertrack1.eta());
            }
          }
        }
      }
      // loop over generated particle
      for (auto& mcParticle : GenParticles) {
        if (TMath::Abs(mcParticle.eta()) > 0.8) // main acceptance
          continue;
        if (TMath::Abs(mcParticle.pdgCode()) == 321 && mcParticle.isPhysicalPrimary()) {
          histos.fill(HIST("hSparsePhiMCGenKaonWeight"), centclass, GetPhiInRange(mcParticle.phi() - psiFT0C), TMath::Power(TMath::Cos(4.0 * GetPhiInRange(mcParticle.phi() - psiFT0C)), 1.0), mcParticle.pt(), mcParticle.eta());
        }
        if (TMath::Abs(mcParticle.y()) > confRapidity) {
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
            if (kCurrentDaughter.pt() > cfgCutPT && TMath::Abs(kCurrentDaughter.eta()) < cfgCutEta) {
              daughtp = true;
            }
            KaonPlus = ROOT::Math::PxPyPzMVector(kCurrentDaughter.px(), kCurrentDaughter.py(), kCurrentDaughter.pz(), massKa);
          } else if (kCurrentDaughter.pdgCode() == -321) {
            if (kCurrentDaughter.pt() > cfgCutPT && TMath::Abs(kCurrentDaughter.eta()) < cfgCutEta) {
              daughtm = true;
            }
            KaonMinus = ROOT::Math::PxPyPzMVector(kCurrentDaughter.px(), kCurrentDaughter.py(), kCurrentDaughter.pz(), massKa);
          }
        }
        if (daughtp && daughtm) {
          auto PhiMinusPsiGen = GetPhiInRange(mcParticle.phi() - psiFT0C);
          histos.fill(HIST("hSparsePhiMCGenWeight"), centclass, PhiMinusPsiGen, TMath::Power(TMath::Cos(4.0 * PhiMinusPsiGen), 1.0), mcParticle.pt(), mcParticle.eta());
        }
      }
    } // rec collision loop

  } // process MC
  PROCESS_SWITCH(phipbpb, processMCPhiWeight, "Process MC Phi Weight", false);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<phipbpb>(cfgc, TaskName{"phipbpb"})};
}
