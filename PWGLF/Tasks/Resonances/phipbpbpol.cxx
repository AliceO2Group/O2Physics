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
// spin alignment task for resonances : dukhishyam mallick

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
struct phipbpbpol {

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
  Configurable<bool> fillRapidity{"fillRapidity", false, "fill rapidity bin"};
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
  ConfigurableAxis configThnAxisCentrality{"configThnAxisCentrality", {8, 0., 80}, "Centrality"};
  ConfigurableAxis configThnAxisPOL{"configThnAxisPOL", {20, -1.0, 1.0}, "Costheta"};
  ConfigurableAxis configThnAxisRapidity{"configThnAxisRapidity", {8, 0, 0.8}, "Rapidity"};
  Configurable<bool> additionalEvsel{"additionalEvsel", false, "Additional event selcection"};
  Configurable<bool> timFrameEvsel{"timFrameEvsel", false, "TPC Time frame boundary cut"};
  Configurable<bool> isMC{"isMC", false, "use MC"};
  Configurable<bool> genacceptancecut{"genacceptancecut", true, "use acceptance cut for generated"};
  Configurable<bool> avoidsplitrackMC{"avoidsplitrackMC", false, "avoid split track in MC"};
  Configurable<bool> islike{"islike", false, "use like"};

  // AxisSpec axisMult{cCentBins, cCentLow, cCentHigh};
  // AxisSpec axisPt{cpTbins, cpTlow, cpThigh};
  // AxisSpec axisMass{cMassbins, cMasslow, cMasshigh};
  // AxisSpec axisMass{cCosbins, cCoslow, cCoshigh};

  // ConfigurableAxis configThnAxisInvMass{"configThnAxisInvMass", {120, 0.66, 1.2}, "#it{M} (GeV/#it{c}^{2})"};
  // ConfigurableAxis configThnAxisPt{"configThnAxisPt", {100, 0.0, 10.}, "#it{p}_{T} (GeV/#it{c})"};
  // ConfigurableAxis configThnAxisCentrality{"configThnAxisCentrality", {8, 0., 80}, "Centrality"};

  // output THnSparses
  Configurable<bool> activateTHnSparseCosThStarHelicity{"activateTHnSparseCosThStarHelicity", true, "Activate the THnSparse with cosThStar w.r.t. helicity axis"};
  Configurable<bool> activateTHnSparseCosThStarProduction{"activateTHnSparseCosThStarProduction", false, "Activate the THnSparse with cosThStar w.r.t. production axis"};
  Configurable<bool> activateTHnSparseCosThStarBeam{"activateTHnSparseCosThStarBeam", false, "Activate the THnSparse with cosThStar w.r.t. beam axis"};
  Configurable<bool> activateTHnSparseCosThStarRandom{"activateTHnSparseCosThStarRandom", false, "Activate the THnSparse with cosThStar w.r.t. random axis"};

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
    const AxisSpec thnAxisCentrality{configThnAxisCentrality, "Centrality"};
    const AxisSpec thnAxisRapidity{configThnAxisRapidity, "Rapidity"};
    const AxisSpec thnAxisPOL{configThnAxisPOL, "POL"};

    // AxisSpec phiAxis = {500, -6.28, 6.28, "phi"};
    // AxisSpec resAxis = {1000, -5, 5, "Res"};
    //  AxisSpec centAxis = {8, 0, 80, "V0M (%)"};

    /*AxisSpec axisMult{cCentBins, cCentLow, cCentHigh};
    AxisSpec axisPt{cpTbins, cpTlow, cpThigh};
    AxisSpec axisMass{cMassbins, cMasslow, cMasshigh};
    AxisSpec axisCostheta{cCosbins, cCoslow, cCoshigh};

    ConfigurableAxis configThnAxisInvMass{"configThnAxisInvMass", {cMassbins,cMasslow,cMasshigh}, "#it{M} (GeV/#it{c}^{2})"};
    ConfigurableAxis configThnAxisPt{"configThnAxisPt", {cpTbins,cpTlow, cpThigh}, "#it{p}_{T} (GeV/#it{c})"};
    ConfigurableAxis configThnAxisCentrality{"configThnAxisCentrality", {cCentBins,cCentLow,cCentHigh}, "Centrality"};
    ConfigurableAxis configThnAxisPOL{"configThnAxisPOL", {cCosbins,cCoslow,cCoshigh}, "Costheta"};
      */
    //  THnSparses
    std::array<int, 4> sparses = {activateTHnSparseCosThStarHelicity, activateTHnSparseCosThStarProduction, activateTHnSparseCosThStarBeam, activateTHnSparseCosThStarRandom};
    if (std::accumulate(sparses.begin(), sparses.end(), 0) == 0) {
      LOGP(fatal, "No output THnSparses enabled");
    } else {
      if (activateTHnSparseCosThStarHelicity) {
        LOGP(info, "THnSparse with cosThStar w.r.t. helicity axis active.");
      }
      if (activateTHnSparseCosThStarProduction) {
        LOGP(info, "THnSparse with cosThStar w.r.t. production axis active.");
      }
      if (activateTHnSparseCosThStarBeam) {
        LOGP(info, "THnSparse with cosThStar w.r.t. beam axis active.");
      }
      if (activateTHnSparseCosThStarRandom) {
        LOGP(info, "THnSparse with cosThStar w.r.t. random axis active.");
      }
    }
    histos.add("hpTvsRapidity", "pT vs Rapidity", kTH2F, {{100, 0.0f, 10.0f}, {300, -1.5f, 1.5f}});
    histos.add("hFTOCvsTPCNoCut", "Mult correlation FT0C vs. TPC without any cut", kTH2F, {{80, 0.0f, 80.0f}, {100, -0.5f, 5999.5f}});
    histos.add("hFTOCvsTPC", "Mult correlation FT0C vs. TPC", kTH2F, {{80, 0.0f, 80.0f}, {100, -0.5f, 5999.5f}});
    histos.add("hFTOCvsTPCSelected", "Mult correlation FT0C vs. TPC after selection", kTH2F, {{80, 0.0f, 80.0f}, {100, -0.5f, 5999.5f}});
    histos.add("hNsigmaKaonTPC", "NsigmaKaon TPC vs pT distribution", kTH2F, {{200, -10.0f, 10.0f}, {100, 0.0f, 10.0f}});
    histos.add("hNsigmaKaonTOF", "NsigmaKaon TOF vs pT distribution", kTH2F, {{200, -10.0f, 10.0f}, {100, 0.0f, 10.0f}});

    histos.add("hCentrality", "Centrality distribution", kTH1F, {{200, 0.0, 200.0}});
    histos.add("hVtxZ", "Vertex distribution in Z;Z (cm)", kTH1F, {{400, -20.0, 20.0}});
    histos.add("hEta", "Eta distribution", kTH1F, {{200, -1.0f, 1.0f}});
    histos.add("hDcaxy", "Dcaxy distribution", kTH1F, {{200, -1.0f, 1.0f}});
    histos.add("hDcaz", "Dcaz distribution", kTH1F, {{200, -1.0f, 1.0f}});

    // define histogram same event pair and mixed event pair
    histos.add("hSparseHESASameEvent", "hSparseHESASameEvent", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisCentrality, thnAxisPOL});
    histos.add("hSparsePPSASameEvent", "hSparsePPSASameEvent", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisCentrality, thnAxisPOL});
    histos.add("hSparseBASASameEvent", "hSparseBASASameEvent", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisCentrality, thnAxisPOL});
    histos.add("hSparseRndASASameEvent", "hSparseRndASASameEvent", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisCentrality, thnAxisPOL});

    // define histogram  mixed event pair
    histos.add("hSparseHESAMixedEvent", "hSparseHESAMixedEvent", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisCentrality, thnAxisPOL});
    histos.add("hSparsePPSAMixedEvent", "hSparsePPSAMixedEvent", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisCentrality, thnAxisPOL});
    histos.add("hSparseBASAMixedEvent", "hSparseBASAMixedEvent", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisCentrality, thnAxisPOL});
    histos.add("hSparseRndASAMixedEvent", "hSparseRndASAMixedEvent", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisCentrality, thnAxisPOL});

    // define histogram  mixed event pair
    histos.add("hSparseHESALikeSignPair", "hSparseHESALikePair", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisCentrality, thnAxisPOL});
    histos.add("hSparsePPSALikeSignPair", "hSparsePPSALikePair", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisCentrality, thnAxisPOL});
    histos.add("hSparseBASALikeSignPair", "hSparseBASALikePair", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisCentrality, thnAxisPOL});
    histos.add("hSparseRndASALikeSignPair", "hSparseRndASALikePair", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisCentrality, thnAxisPOL});

    // MC histogram
    if (isMC) {
      histos.add("hMC", "MC Event statistics", kTH1F, {{10, 0.0f, 10.0f}});
      histos.add("h1PhiRecsplit", "Phi meson Rec split", kTH1F, {{100, 0.0f, 10.0f}});
      histos.add("CentPercentileMCRecHist", "MC Centrality", kTH1F, {{100, 0.0f, 100.0f}});

      histos.add("hSparseHESASameEvent_costheta_SA_MCGen", "hSparseHESASameEvent_costheta_SA_MCGen", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisCentrality, thnAxisPOL});
      histos.add("hSparsePPSASameEvent_costheta_SA_MCGen", "hSparsePPSASameEvent_costheta_SA_MCGen", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisCentrality, thnAxisPOL});
      histos.add("hSparseBASASameEvent_costheta_SA_MCGen", "hSparseBASASameEvent_costheta_SA_MCGen", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisCentrality, thnAxisPOL});
      histos.add("hSparseRndASASameEvent_costheta_SA_MCGen", "hSparseRndASASameEvent_costheta_SA_MCGen", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisCentrality, thnAxisPOL});

      histos.add("hSparseHESASameEvent_costheta_SA_MCRec", "hSparseHESASameEvent_costheta_SA_MCRec", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisCentrality, thnAxisPOL});
      histos.add("hSparsePPSASameEvent_costheta_SA_MCRec", "hSparsePPSASameEvent_costheta_SA_MCRec", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisCentrality, thnAxisPOL});
      histos.add("hSparseBASASameEvent_costheta_SA_MCRec", "hSparseBASASameEvent_costheta_SA_MCRec", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisCentrality, thnAxisPOL});
      histos.add("hSparseRndASASameEvent_costheta_SA_MCRec", "hSparseRndASASameEvent_costheta_SA_MCRec", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisCentrality, thnAxisPOL});
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

  ROOT::Math::PxPyPzMVector daughter1, daughter2;

  void processSameEvent(EventCandidates::iterator const& collision, TrackCandidates const& tracks, aod::BCs const&)
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

    histos.fill(HIST("hFTOCvsTPC"), centrality, multTPC);
    if (additionalEvsel && !eventSelected(collision, centrality)) {
      return;
    }
    histos.fill(HIST("hFTOCvsTPCSelected"), centrality, multTPC);
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
      histos.fill(HIST("hNsigmaKaonTPC"), track1.tpcNSigmaKa(), track1.pt());
      histos.fill(HIST("hNsigmaKaonTOF"), track1.tofNSigmaKa(), track1.pt());
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

        daughter1 = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), massKa); // Kplus
        daughter2 = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), massKa); // Kminus
        PhiMesonMother = daughter1 + daughter2;

        histos.fill(HIST("hpTvsRapidity"), PhiMesonMother.Pt(), PhiMesonMother.Rapidity());

        if (TMath::Abs(PhiMesonMother.Rapidity()) > confRapidity) {
          continue;
        }
        auto phiRandom = gRandom->Uniform(0.f, constants::math::TwoPI);
        auto thetaRandom = gRandom->Uniform(0.f, constants::math::PI);
        ROOT::Math::PxPyPzMVector fourVecDau = ROOT::Math::PxPyPzMVector(daughter1.px(), daughter1.py(), daughter1.pz(), massKa); // Kplus
        ROOT::Math::PxPyPzMVector fourVecMother = ROOT::Math::PxPyPzMVector(PhiMesonMother.px(), PhiMesonMother.py(), PhiMesonMother.pz(), PhiMesonMother.M());
        ROOT::Math::Boost boost{fourVecMother.BoostToCM()};
        ROOT::Math::PxPyPzMVector fourVecDauCM = boost(fourVecDau);
        ROOT::Math::XYZVector threeVecDauCM = fourVecDauCM.Vect();

        if (activateTHnSparseCosThStarHelicity) {
          ROOT::Math::XYZVector helicityVec = fourVecMother.Vect();
          auto cosThetaStarHelicity = helicityVec.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2()) / std::sqrt(helicityVec.Mag2());
          histos.fill(HIST("hSparseHESASameEvent"), PhiMesonMother.M(), PhiMesonMother.Pt(), centrality, cosThetaStarHelicity);
          // ROOT::Math::XYZVectorF zaxis_HE{helicityVec.Unit()};
          // ROOT::Math::XYZVectorF v1_CM{threeVecDauCM.Unit()};
          // float cosThetaStarHelicity = zaxis_HE.Dot(v1_CM);
        }
        if (activateTHnSparseCosThStarProduction) {
          ROOT::Math::XYZVector normalVec = ROOT::Math::XYZVector(PhiMesonMother.py(), -PhiMesonMother.px(), 0.f);
          auto cosThetaStarProduction = normalVec.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2()) / std::sqrt(normalVec.Mag2());
          histos.fill(HIST("hSparsePPSASameEvent"), PhiMesonMother.M(), PhiMesonMother.Pt(), centrality, cosThetaStarProduction);
        }
        if (activateTHnSparseCosThStarBeam) {
          ROOT::Math::XYZVector beamVec = ROOT::Math::XYZVector(0.f, 0.f, 1.f);
          auto cosThetaStarBeam = beamVec.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2());
          histos.fill(HIST("hSparseBASASameEvent"), PhiMesonMother.M(), PhiMesonMother.Pt(), centrality, cosThetaStarBeam);
        }
        if (activateTHnSparseCosThStarRandom) {
          ROOT::Math::XYZVector randomVec = ROOT::Math::XYZVector(std::sin(thetaRandom) * std::cos(phiRandom), std::sin(thetaRandom) * std::sin(phiRandom), std::cos(thetaRandom));
          auto cosThetaStarRandom = randomVec.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2());
          histos.fill(HIST("hSparseRndASASameEvent"), PhiMesonMother.M(), PhiMesonMother.Pt(), centrality, cosThetaStarRandom);
        }
      }
    }
  }
  PROCESS_SWITCH(phipbpbpol, processSameEvent, "Process Same event", true);
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
        daughter1 = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), massKa);
        daughter2 = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), massKa);
        PhiMesonMother = daughter1 + daughter2;
        if (TMath::Abs(PhiMesonMother.Rapidity()) > confRapidity) {
          continue;
        }

        auto phiRandom = gRandom->Uniform(0.f, constants::math::TwoPI);
        auto thetaRandom = gRandom->Uniform(0.f, constants::math::PI);
        ROOT::Math::PxPyPzMVector fourVecDau = ROOT::Math::PxPyPzMVector(daughter1.px(), daughter1.py(), daughter1.pz(), massKa); // Kplus
        ROOT::Math::PxPyPzMVector fourVecMother = ROOT::Math::PxPyPzMVector(PhiMesonMother.px(), PhiMesonMother.py(), PhiMesonMother.pz(), PhiMesonMother.M());

        ROOT::Math::Boost boost{fourVecMother.BoostToCM()};
        ROOT::Math::PxPyPzMVector fourVecDauCM = boost(fourVecDau);
        ROOT::Math::XYZVector threeVecDauCM = fourVecDauCM.Vect();

        if (activateTHnSparseCosThStarHelicity) {
          ROOT::Math::XYZVector helicityVec = fourVecMother.Vect();
          auto cosThetaStarHelicity = helicityVec.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2()) / std::sqrt(helicityVec.Mag2());
          histos.fill(HIST("hSparseHESAMixedEvent"), PhiMesonMother.M(), PhiMesonMother.Pt(), centrality, cosThetaStarHelicity);
        }

        if (activateTHnSparseCosThStarProduction) {
          ROOT::Math::XYZVector normalVec = ROOT::Math::XYZVector(PhiMesonMother.py(), -PhiMesonMother.px(), 0.f);
          auto cosThetaStarProduction = normalVec.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2()) / std::sqrt(normalVec.Mag2());
          histos.fill(HIST("hSparsePPSAMixedEvent"), PhiMesonMother.M(), PhiMesonMother.Pt(), centrality, cosThetaStarProduction);
        }
        if (activateTHnSparseCosThStarBeam) {
          ROOT::Math::XYZVector beamVec = ROOT::Math::XYZVector(0.f, 0.f, 1.f);
          auto cosThetaStarBeam = beamVec.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2());
          histos.fill(HIST("hSparseBASAMixedEvent"), PhiMesonMother.M(), PhiMesonMother.Pt(), centrality, cosThetaStarBeam);
        }
        if (activateTHnSparseCosThStarRandom) {
          ROOT::Math::XYZVector randomVec = ROOT::Math::XYZVector(std::sin(thetaRandom) * std::cos(phiRandom), std::sin(thetaRandom) * std::sin(phiRandom), std::cos(thetaRandom));
          auto cosThetaStarRandom = randomVec.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2());
          histos.fill(HIST("hSparseRndASAMixedEvent"), PhiMesonMother.M(), PhiMesonMother.Pt(), centrality, cosThetaStarRandom);
        }
      }
    }
  }
  PROCESS_SWITCH(phipbpbpol, processMixedEvent, "Process Mixed event", true);
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
          daughter1 = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), massKa);
          daughter2 = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), massKa);
        } else if (track1.sign() < 0 && track2.sign() > 0) {
          daughter2 = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), massKa);
          daughter1 = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), massKa);
        }

        PhiMesonMother = daughter1 + daughter2;
        if (TMath::Abs(PhiMesonMother.Rapidity()) > confRapidity) {
          continue;
        }
        auto phiRandom = gRandom->Uniform(0.f, constants::math::TwoPI);
        auto thetaRandom = gRandom->Uniform(0.f, constants::math::PI);
        ROOT::Math::PxPyPzMVector fourVecDau = ROOT::Math::PxPyPzMVector(daughter1.px(), daughter1.py(), daughter1.pz(), massKa);
        ROOT::Math::PxPyPzMVector fourVecMother = ROOT::Math::PxPyPzMVector(PhiMesonMother.px(), PhiMesonMother.py(), PhiMesonMother.pz(), PhiMesonMother.M());
        ROOT::Math::Boost boost{fourVecMother.BoostToCM()};
        ROOT::Math::PxPyPzMVector fourVecDauCM = boost(fourVecDau);
        ROOT::Math::XYZVector threeVecDauCM = fourVecDauCM.Vect();
        if (activateTHnSparseCosThStarHelicity) {
          ROOT::Math::XYZVector helicityVec = fourVecMother.Vect();
          auto cosThetaStarHelicity = helicityVec.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2()) / std::sqrt(helicityVec.Mag2());
          histos.fill(HIST("hSparseHESALikeSignPair"), PhiMesonMother.M(), PhiMesonMother.Pt(), centrality, cosThetaStarHelicity);
        }

        if (activateTHnSparseCosThStarProduction) {
          ROOT::Math::XYZVector normalVec = ROOT::Math::XYZVector(PhiMesonMother.py(), -PhiMesonMother.px(), 0.f);
          auto cosThetaStarProduction = normalVec.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2()) / std::sqrt(normalVec.Mag2());
          histos.fill(HIST("hSparsePPSALikeSignPair"), PhiMesonMother.M(), PhiMesonMother.Pt(), centrality, cosThetaStarProduction);
        }
        if (activateTHnSparseCosThStarBeam) {
          ROOT::Math::XYZVector beamVec = ROOT::Math::XYZVector(0.f, 0.f, 1.f);
          auto cosThetaStarBeam = beamVec.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2());
          histos.fill(HIST("hSparseBASALikeSignPair"), PhiMesonMother.M(), PhiMesonMother.Pt(), centrality, cosThetaStarBeam);
        }
        if (activateTHnSparseCosThStarRandom) {
          ROOT::Math::XYZVector randomVec = ROOT::Math::XYZVector(std::sin(thetaRandom) * std::cos(phiRandom), std::sin(thetaRandom) * std::sin(phiRandom), std::cos(thetaRandom));
          auto cosThetaStarRandom = randomVec.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2());
          histos.fill(HIST("hSparseRndASALikeSignPair"), PhiMesonMother.M(), PhiMesonMother.Pt(), centrality, cosThetaStarRandom);
        }
      }
    }
  }
  PROCESS_SWITCH(phipbpbpol, processMixedEventOpti, "Process Mixed event new", true);
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
      // auto psiFT0C = 0.0;
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
                daughter1 = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), massKa);
                daughter2 = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), massKa);
              }
              if (track1.sign() < 0 && track2.sign() > 0) {
                daughter2 = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), massKa);
                daughter1 = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), massKa);
              }
              PhiMesonMother = daughter1 + daughter2;

              if (TMath::Abs(PhiMesonMother.Rapidity()) > confRapidity) {
                continue;
              }

              auto phiRandom = gRandom->Uniform(0.f, constants::math::TwoPI);
              auto thetaRandom = gRandom->Uniform(0.f, constants::math::PI);

              ROOT::Math::PxPyPzMVector fourVecDau = ROOT::Math::PxPyPzMVector(daughter1.px(), daughter1.py(), daughter1.pz(), massKa);
              ROOT::Math::PxPyPzMVector fourVecMother = ROOT::Math::PxPyPzMVector(PhiMesonMother.px(), PhiMesonMother.py(), PhiMesonMother.pz(), PhiMesonMother.M());
              ROOT::Math::Boost boost{fourVecMother.BoostToCM()};
              ROOT::Math::PxPyPzMVector fourVecDauCM = boost(fourVecDau);
              ROOT::Math::XYZVector threeVecDauCM = fourVecDauCM.Vect();
              if (activateTHnSparseCosThStarHelicity) {
                ROOT::Math::XYZVector helicityVec = fourVecMother.Vect();
                auto cosThetaStarHelicity = helicityVec.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2()) / std::sqrt(helicityVec.Mag2());
                histos.fill(HIST("hSparseHESASameEvent_costheta_SA_MCRec"), PhiMesonMother.M(), PhiMesonMother.Pt(), centrality, cosThetaStarHelicity);
              }

              if (activateTHnSparseCosThStarProduction) {
                ROOT::Math::XYZVector normalVec = ROOT::Math::XYZVector(PhiMesonMother.py(), -PhiMesonMother.px(), 0.f);
                auto cosThetaStarProduction = normalVec.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2()) / std::sqrt(normalVec.Mag2());
                histos.fill(HIST("hSparsePPSASameEvent_costheta_SA_MCRec"), PhiMesonMother.M(), PhiMesonMother.Pt(), centrality, cosThetaStarProduction);
              }
              if (activateTHnSparseCosThStarBeam) {
                ROOT::Math::XYZVector beamVec = ROOT::Math::XYZVector(0.f, 0.f, 1.f);
                auto cosThetaStarBeam = beamVec.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2());
                histos.fill(HIST("hSparseBASASameEvent_costheta_SA_MCRec"), PhiMesonMother.M(), PhiMesonMother.Pt(), centrality, cosThetaStarBeam);
              }
              if (activateTHnSparseCosThStarRandom) {
                ROOT::Math::XYZVector randomVec = ROOT::Math::XYZVector(std::sin(thetaRandom) * std::cos(phiRandom), std::sin(thetaRandom) * std::sin(phiRandom), std::cos(thetaRandom));
                auto cosThetaStarRandom = randomVec.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2());
                histos.fill(HIST("hSparseRndASASameEvent_costheta_SA_MCRec"), PhiMesonMother.M(), PhiMesonMother.Pt(), centrality, cosThetaStarRandom);
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
            daughter1 = ROOT::Math::PxPyPzMVector(kCurrentDaughter.px(), kCurrentDaughter.py(), kCurrentDaughter.pz(), massKa);
          } else if (kCurrentDaughter.pdgCode() == -321) {
            if (genacceptancecut && kCurrentDaughter.pt() > cfgCutPT && TMath::Abs(kCurrentDaughter.eta()) < cfgCutEta) {
              daughtm = true;
            }
            if (!genacceptancecut) {
              daughtm = true;
            }
            daughter2 = ROOT::Math::PxPyPzMVector(kCurrentDaughter.px(), kCurrentDaughter.py(), kCurrentDaughter.pz(), massKa);
          }
        }
        if (daughtp && daughtm) {
          PhiMesonMother = daughter1 + daughter2;
          auto phiRandom = gRandom->Uniform(0.f, constants::math::TwoPI);
          auto thetaRandom = gRandom->Uniform(0.f, constants::math::PI);

          ROOT::Math::PxPyPzMVector fourVecDau = ROOT::Math::PxPyPzMVector(daughter1.px(), daughter1.py(), daughter1.pz(), massKa);
          ROOT::Math::PxPyPzMVector fourVecMother = ROOT::Math::PxPyPzMVector(PhiMesonMother.px(), PhiMesonMother.py(), PhiMesonMother.pz(), PhiMesonMother.M());
          ROOT::Math::Boost boost{fourVecMother.BoostToCM()};
          ROOT::Math::PxPyPzMVector fourVecDauCM = boost(fourVecDau);
          ROOT::Math::XYZVector threeVecDauCM = fourVecDauCM.Vect();

          if (activateTHnSparseCosThStarHelicity) {
            ROOT::Math::XYZVector helicityVec = fourVecMother.Vect();
            auto cosThetaStarHelicity = helicityVec.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2()) / std::sqrt(helicityVec.Mag2());
            histos.fill(HIST("hSparseHESASameEvent_costheta_SA_MCGen"), PhiMesonMother.M(), PhiMesonMother.Pt(), centrality, cosThetaStarHelicity);
          }

          if (activateTHnSparseCosThStarProduction) {
            ROOT::Math::XYZVector normalVec = ROOT::Math::XYZVector(PhiMesonMother.py(), -PhiMesonMother.px(), 0.f);
            auto cosThetaStarProduction = normalVec.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2()) / std::sqrt(normalVec.Mag2());
            histos.fill(HIST("hSparsePPSASameEvent_costheta_SA_MCGen"), PhiMesonMother.M(), PhiMesonMother.Pt(), centrality, cosThetaStarProduction);
          }
          if (activateTHnSparseCosThStarBeam) {
            ROOT::Math::XYZVector beamVec = ROOT::Math::XYZVector(0.f, 0.f, 1.f);
            auto cosThetaStarBeam = beamVec.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2());
            histos.fill(HIST("hSparseBASASameEvent_costheta_SA_MCGen"), PhiMesonMother.M(), PhiMesonMother.Pt(), centrality, cosThetaStarBeam);
          }
          if (activateTHnSparseCosThStarRandom) {
            ROOT::Math::XYZVector randomVec = ROOT::Math::XYZVector(std::sin(thetaRandom) * std::cos(phiRandom), std::sin(thetaRandom) * std::sin(phiRandom), std::cos(thetaRandom));
            auto cosThetaStarRandom = randomVec.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2());
            histos.fill(HIST("hSparseRndASASameEvent_costheta_SA_MCGen"), PhiMesonMother.M(), PhiMesonMother.Pt(), centrality, cosThetaStarRandom);
          }
        }
      }
    } // rec collision loop

  } // process MC
  PROCESS_SWITCH(phipbpbpol, processMC, "Process MC", false);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<phipbpbpol>(cfgc, TaskName{"phipbpb"})};
}
