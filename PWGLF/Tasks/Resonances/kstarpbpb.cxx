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
struct kstarpbpb {

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
  // track
  Configurable<float> cfgCutCharge{"cfgCutCharge", 0.0, "cut on Charge"};
  Configurable<float> cfgCutPT{"cfgCutPT", 0.2, "PT cut on daughter track"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8, "Eta cut on daughter track"};
  Configurable<float> cfgCutDCAxy{"cfgCutDCAxy", 2.0f, "DCAxy range for tracks"};
  Configurable<float> cfgCutDCAz{"cfgCutDCAz", 2.0f, "DCAz range for tracks"};
  Configurable<bool> useGlobalTrack{"useGlobalTrack", true, "use Global track"};
  Configurable<float> nsigmaCutTPC{"nsigmacutTPC", 3.0, "Value of the TPC Nsigma cut"};
  Configurable<float> nsigmaCutCombined{"nsigmaCutCombined", 3.0, "Value of the TOF Nsigma cut"};
  Configurable<int> cfgNoMixedEvents{"cfgNoMixedEvents", 1, "Number of mixed events per event"};
  Configurable<int> cfgITScluster{"cfgITScluster", 0, "Number of ITS cluster"};
  Configurable<int> cfgTPCcluster{"cfgTPCcluster", 70, "Number of TPC cluster"};
  Configurable<double> confRapidity{"confRapidity", 0.5, "Rapidity cut"};
  ConfigurableAxis configThnAxisInvMass{"configThnAxisInvMass", {180, 0.6, 1.5}, "#it{M} (GeV/#it{c}^{2})"};
  ConfigurableAxis configThnAxisPt{"configThnAxisPt", {100, 0.0, 10.}, "#it{p}_{T} (GeV/#it{c})"};
  ConfigurableAxis configThnAxisCentrality{"configThnAxisCentrality", {8, 0., 80}, "Centrality"};
  ConfigurableAxis configThnAxisPhiminusPsi{"configThnAxisPhiminusPsi", {6, 0.0, TMath::Pi()}, "#phi - #psi"};
  ConfigurableAxis configThnAxisV2{"configThnAxisV2", {200, -1, 1}, "V2"};
  Configurable<bool> additionalEvsel{"additionalEvsel", false, "Additional event selcection"};
  Configurable<bool> timFrameEvsel{"timFrameEvsel", false, "TPC Time frame boundary cut"};
  Configurable<bool> ispTdepPID{"ispTdepPID", true, "pT dependent PID"};
  Configurable<bool> additionalQAplots{"additionalQAplots", true, "Additional QA plots"};
  Configurable<float> confMinRot{"confMinRot", 5.0 * TMath::Pi() / 6.0, "Minimum of rotation"};
  Configurable<float> confMaxRot{"confMaxRot", 7.0 * TMath::Pi() / 6.0, "Maximum of rotation"};
  Configurable<int> nBkgRotations{"nBkgRotations", 9, "Number of rotated copies (background) per each original candidate"};
  Configurable<bool> fillRotation{"fillRotation", true, "fill rotation"};

  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter acceptanceFilter = (nabs(aod::track::eta) < cfgCutEta && nabs(aod::track::pt) > cfgCutPT);
  Filter DCAcutFilter = (nabs(aod::track::dcaXY) < cfgCutDCAxy) && (nabs(aod::track::dcaZ) < cfgCutDCAz);

  using EventCandidates = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::FV0Mults, aod::TPCMults, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0Cs, aod::CentFT0As, aod::EPCalibrationTables, aod::Mults>>;
  using TrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullKa, aod::pidTOFFullKa, aod::pidTPCFullPi, aod::pidTOFFullPi>>;

  using EventCandidatesMC = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>;
  using TrackCandidatesMC = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullKa, aod::pidTOFFullKa, aod::McTrackLabels>>;

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
    const AxisSpec thnAxisPhiminusPsi{configThnAxisPhiminusPsi, "#phi - #psi"};
    const AxisSpec thnAxisCentrality{configThnAxisCentrality, "Centrality (%)"};
    const AxisSpec thnAxisV2{configThnAxisV2, "V2"};
    AxisSpec phiAxis = {500, -6.28, 6.28, "phi"};
    AxisSpec resAxis = {400, -2, 2, "Res"};
    AxisSpec centAxis = {8, 0, 80, "V0M (%)"};
    AxisSpec occupancyAxis = {1500, 0, 1500, "Occupancy"};

    histos.add("hpTvsRapidity", "pT vs Rapidity", kTH2F, {{100, 0.0f, 10.0f}, {300, -1.5f, 1.5f}});
    histos.add("hFTOCvsTPC", "Mult correlation FT0C vs. TPC", kTH2F, {{80, 0.0f, 80.0f}, {100, -0.5f, 5999.5f}});
    histos.add("hFTOCvsTPCSelected", "Mult correlation FT0C vs. TPC after selection", kTH2F, {{80, 0.0f, 80.0f}, {100, -0.5f, 5999.5f}});
    histos.add("hCentrality", "Centrality distribution", kTH1F, {{200, 0.0, 200.0}});
    histos.add("hOccupancy", "Occupancy distribution", kTH1F, {occupancyAxis});
    histos.add("hVtxZ", "Vertex distribution in Z;Z (cm)", kTH1F, {{400, -20.0, 20.0}});
    histos.add("hPsiFT0C", "PsiFT0C", kTH3F, {centAxis, occupancyAxis, phiAxis});
    histos.add("hPsiFT0A", "PsiFT0A", kTH3F, {centAxis, occupancyAxis, phiAxis});
    histos.add("hPsiTPC", "PsiTPC", kTH3F, {centAxis, occupancyAxis, phiAxis});
    histos.add("hSparseV2SASameEvent_V2", "hSparseV2SASameEvent_V2", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisV2, thnAxisCentrality});
    histos.add("hSparseV2SAMixedEvent_V2", "hSparseV2SAMixedEvent_V2", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisV2, thnAxisCentrality});
    histos.add("hSparseV2SASameEventRotational_V2", "hSparseV2SASameEventRotational_V2", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisV2, thnAxisCentrality});
    histos.add("hMC", "MC Event statistics", kTH1F, {{6, 0.0f, 6.0f}});

    // histogram for resolution
    histos.add("ResFT0CTPC", "ResFT0CTPC", kTH2F, {centAxis, resAxis});
    histos.add("ResFT0CFT0A", "ResFT0CFT0A", kTH2F, {centAxis, resAxis});
    histos.add("ResFT0ATPC", "ResFT0ATPC", kTH2F, {centAxis, resAxis});
    if (additionalQAplots) {
      // DCA QA
      histos.add("QAbefore/trkDCAxyka", "DCAxy distribution of kaon track candidates", HistType::kTH1F, {{150, 0.0f, 1.0f}});
      histos.add("QAbefore/trkDCAzka", "DCAz distribution of kaon track candidates", HistType::kTH1F, {{150, 0.0f, 1.0f}});
      histos.add("QAafter/trkDCAxyka", "DCAxy distribution of kaon track candidates", HistType::kTH1F, {{150, 0.0f, 1.0f}});
      histos.add("QAafter/trkDCAzka", "DCAz distribution of kaon track candidates", HistType::kTH1F, {{150, 0.0f, 1.0f}});
      // PID QA before cuts
      histos.add("QAbefore/TOF_TPC_Mapka_allka", "TOF + TPC Combined PID for Kaon;#sigma_{TOF}^{Kaon};#sigma_{TPC}^{Kaon}", {HistType::kTH2D, {{100, -6, 6}, {100, -6, 6}}});
      histos.add("QAbefore/TOF_Nsigma_allka", "TOF NSigma for Kaon;#it{p}_{T} (GeV/#it{c});#sigma_{TOF}^{Kaon};", {HistType::kTH2D, {{200, 0.0, 20.0}, {100, -6, 6}}});
      histos.add("QAbefore/TPC_Nsigma_allka", "TPC NSigma for Kaon;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Kaon};", {HistType::kTH2D, {{200, 0.0, 20.0}, {100, -6, 6}}});
      // PID QA after cuts
      histos.add("QAafter/TOF_TPC_Mapka_allka", "TOF + TPC Combined PID for Kaon;#sigma_{TOF}^{Kaon};#sigma_{TPC}^{Kaon}", {HistType::kTH2D, {{100, -6, 6}, {100, -6, 6}}});
      histos.add("QAafter/TOF_Nsigma_allka", "TOF NSigma for Kaon;#it{p}_{T} (GeV/#it{c});#sigma_{TOF}^{Kaon};", {HistType::kTH2D, {{200, 0.0, 20.0}, {100, -6, 6}}});
      histos.add("QAafter/TPC_Nsigma_allka", "TPC NSigma for Kaon;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Kaon};", {HistType::kTH2D, {{200, 0.0, 20.0}, {100, -6, 6}}});

      // DCA QA
      histos.add("QAbefore/trkDCAxypi", "DCAxy distribution of pion track candidates", HistType::kTH1F, {{150, 0.0f, 1.0f}});
      histos.add("QAbefore/trkDCAzpi", "DCAz distribution of pion track candidates", HistType::kTH1F, {{150, 0.0f, 1.0f}});
      histos.add("QAafter/trkDCAxypi", "DCAxy distribution of pion track candidates", HistType::kTH1F, {{150, 0.0f, 1.0f}});
      histos.add("QAafter/trkDCAzpi", "DCAz distribution of pion track candidates", HistType::kTH1F, {{150, 0.0f, 1.0f}});
      // PID QA before cuts
      histos.add("QAbefore/TOF_TPC_Mapka_allpi", "TOF + TPC Combined PID for pion;#sigma_{TOF}^{pion};#sigma_{TPC}^{pion}", {HistType::kTH2D, {{100, -6, 6}, {100, -6, 6}}});
      histos.add("QAbefore/TOF_Nsigma_allpi", "TOF NSigma for pion;#it{p}_{T} (GeV/#it{c});#sigma_{TOF}^{pion};", {HistType::kTH2D, {{200, 0.0, 20.0}, {100, -6, 6}}});
      histos.add("QAbefore/TPC_Nsigma_allpi", "TPC NSigma for pion;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{pion};", {HistType::kTH2D, {{200, 0.0, 20.0}, {100, -6, 6}}});
      // PID QA after cuts
      histos.add("QAafter/TOF_TPC_Mapka_allpi", "TOF + TPC Combined PID for pion;#sigma_{TOF}^{pion};#sigma_{TPC}^{pion}", {HistType::kTH2D, {{100, -6, 6}, {100, -6, 6}}});
      histos.add("QAafter/TOF_Nsigma_allpi", "TOF NSigma for pion;#it{p}_{T} (GeV/#it{c});#sigma_{TOF}^{pion};", {HistType::kTH2D, {{200, 0.0, 20.0}, {100, -6, 6}}});
      histos.add("QAafter/TPC_Nsigma_allpi", "TPC NSigma for pion;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{pion};", {HistType::kTH2D, {{200, 0.0, 20.0}, {100, -6, 6}}});
    }
    if (fillRotation) {
      histos.add("hRotation", "hRotation", kTH1F, {{360, 0.0, 2.0 * TMath::Pi()}});
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
  double massPi = o2::constants::physics::MassPiMinus;

  template <typename TCollision>
  bool eventSelected(TCollision collision, const int& /*multTrk*/, const float& centrality)
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
  bool selectionPIDpTdependent(const T& candidate, int PID)
  {
    if (PID == 0) {
      if (candidate.p() < 0.5 && TMath::Abs(candidate.tpcNSigmaKa()) < nsigmaCutTPC) {
        return true;
      }
      if (candidate.p() >= 0.5 && candidate.hasTOF() && ((candidate.tofNSigmaKa() * candidate.tofNSigmaKa()) + (candidate.tpcNSigmaKa() * candidate.tpcNSigmaKa())) < (nsigmaCutCombined * nsigmaCutCombined)) {
        return true;
      }
    } else if (PID == 1) {
      if (candidate.p() < 0.5 && TMath::Abs(candidate.tpcNSigmaPi()) < nsigmaCutTPC) {
        return true;
      }
      if (candidate.p() >= 0.5 && candidate.hasTOF() && ((candidate.tofNSigmaPi() * candidate.tofNSigmaPi()) + (candidate.tpcNSigmaPi() * candidate.tpcNSigmaPi())) < (nsigmaCutCombined * nsigmaCutCombined)) {
        return true;
      }
    }
    return false;
  }

  template <typename T>
  bool selectionPID(const T& candidate, int PID)
  {
    if (PID == 0) {
      if (candidate.hasTOF() && ((candidate.tofNSigmaKa() * candidate.tofNSigmaKa()) + (candidate.tpcNSigmaKa() * candidate.tpcNSigmaKa())) < (nsigmaCutCombined * nsigmaCutCombined)) {
        return true;
      }
    } else if (PID == 1) {
      if (candidate.hasTOF() && ((candidate.tofNSigmaPi() * candidate.tofNSigmaPi()) + (candidate.tpcNSigmaPi() * candidate.tpcNSigmaPi())) < (nsigmaCutCombined * nsigmaCutCombined)) {
        return true;
      }
    }
    return false;
  }

  template <typename T>
  bool selectionPIDNew(const T& candidate, int PID)
  {
    if (PID == 0) {
      if (candidate.p() < 0.6 && TMath::Abs(candidate.tpcNSigmaKa()) < 2.0) {
        return true;
      }
      if (candidate.p() >= 0.6 && candidate.p() < 3.0 && candidate.hasTOF() && candidate.tpcNSigmaKa() > -2.0 && candidate.tpcNSigmaKa() < 3.0 && TMath::Abs(candidate.tofNSigmaKa()) < 2.0) {
        return true;
      }
    } else if (PID == 1) {
      if (candidate.p() < 1.0 && TMath::Abs(candidate.tpcNSigmaPi()) < 2.0) {
        return true;
      }
      if (candidate.p() >= 1.0 && candidate.p() < 3.0 && candidate.hasTOF() && candidate.tpcNSigmaPi() > -2.0 && candidate.tpcNSigmaPi() < 3.0 && TMath::Abs(candidate.tofNSigmaPi()) < 2.0) {
        return true;
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

  ConfigurableAxis axisVertex{"axisVertex", {20, -10, 10}, "vertex axis for bin"};
  ConfigurableAxis axisMultiplicityClass{"axisMultiplicityClass", {20, 0, 100}, "multiplicity percentile for bin"};
  ConfigurableAxis axisEPAngle{"axisEPAngle", {6, -TMath::Pi() / 2, TMath::Pi() / 2}, "event plane angle"};

  using BinningTypeVertexContributor = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0C, aod::epcalibrationtable::PsiFT0C>;
  ROOT::Math::PxPyPzMVector KstarMother, daughter1, daughter2, kaonrot, kstarrot;

  void processSameEvent(EventCandidates::iterator const& collision, TrackCandidates const& tracks, aod::BCs const&)
  {
    histos.fill(HIST("hMC"), 0.5);
    if (!collision.sel8()) {
      return;
    }
    histos.fill(HIST("hMC"), 1.5);
    if (!collision.triggereventep()) {
      return;
    }
    histos.fill(HIST("hMC"), 2.5);
    if (timFrameEvsel && (!collision.selection_bit(aod::evsel::kNoTimeFrameBorder) || !collision.selection_bit(aod::evsel::kNoITSROFrameBorder))) {
      return;
    }
    histos.fill(HIST("hMC"), 3.5);
    auto posThisColl = posTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto negThisColl = negTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto centrality = collision.centFT0C();
    auto multTPC = collision.multNTracksPV();
    auto psiFT0C = collision.psiFT0C();
    auto psiFT0A = collision.psiFT0A();
    auto psiTPC = collision.psiTPC();
    int occupancy = collision.trackOccupancyInTimeRange();
    if (occupancy >= 1500) // occupancy info is available for this collision (*)
    {
      return;
    }
    histos.fill(HIST("hFTOCvsTPC"), centrality, multTPC);
    if (additionalEvsel && !eventSelected(collision, tracks.size(), centrality)) {
      return;
    }
    histos.fill(HIST("hFTOCvsTPCSelected"), centrality, multTPC);
    histos.fill(HIST("hPsiFT0C"), centrality, occupancy, psiFT0C);
    histos.fill(HIST("hPsiFT0A"), centrality, occupancy, psiFT0A);
    histos.fill(HIST("hPsiTPC"), centrality, occupancy, psiTPC);
    histos.fill(HIST("ResFT0CTPC"), centrality, TMath::Cos(2.0 * (psiFT0C - psiTPC)));
    histos.fill(HIST("ResFT0CFT0A"), centrality, TMath::Cos(2.0 * (psiFT0C - psiFT0A)));
    histos.fill(HIST("ResFT0ATPC"), centrality, TMath::Cos(2.0 * (psiTPC - psiFT0A)));
    histos.fill(HIST("hCentrality"), centrality);
    histos.fill(HIST("hOccupancy"), occupancy);
    histos.fill(HIST("hVtxZ"), collision.posZ());

    for (auto track1 : posThisColl) {
      if (!selectionTrack(track1)) {
        continue;
      }
      if (additionalQAplots) {
        histos.fill(HIST("QAbefore/TPC_Nsigma_allka"), track1.pt(), track1.tpcNSigmaKa());
        histos.fill(HIST("QAbefore/TOF_Nsigma_allka"), track1.pt(), track1.tofNSigmaKa());
        histos.fill(HIST("QAbefore/trkDCAxyka"), track1.dcaXY());
        histos.fill(HIST("QAbefore/trkDCAzka"), track1.dcaZ());
        histos.fill(HIST("QAbefore/TOF_TPC_Mapka_allka"), track1.tofNSigmaKa(), track1.tpcNSigmaKa());
      }

      bool track1pion = false;
      bool track1kaon = false;
      if (ispTdepPID && !(selectionPIDNew(track1, 0) || selectionPIDNew(track1, 1))) {
        continue;
      }
      if (!ispTdepPID && !(selectionPID(track1, 0) || selectionPID(track1, 1))) {
        continue;
      }
      for (auto track2 : negThisColl) {
        bool track2pion = false;
        bool track2kaon = false;
        if (!selectionTrack(track2)) {
          continue;
        }
        if (additionalQAplots) {
          histos.fill(HIST("QAbefore/TOF_TPC_Mapka_allpi"), track2.tofNSigmaPi(), track2.tpcNSigmaPi());
          histos.fill(HIST("QAbefore/TPC_Nsigma_allpi"), track2.pt(), track2.tpcNSigmaPi());
          histos.fill(HIST("QAbefore/TOF_Nsigma_allpi"), track2.pt(), track2.tofNSigmaPi());
          histos.fill(HIST("QAbefore/trkDCAxypi"), track2.dcaXY());
          histos.fill(HIST("QAbefore/trkDCAzpi"), track2.dcaZ());
        }
        if (ispTdepPID && !(selectionPIDNew(track2, 0) || selectionPIDNew(track2, 1))) {
          continue;
        }
        if (!ispTdepPID && !(selectionPID(track2, 0) || selectionPID(track2, 1))) {
          continue;
        }
        if (track1.sign() * track2.sign() > 0) {
          continue;
        }

        if (ispTdepPID) {
          if (selectionPIDNew(track1, 1) && selectionPIDNew(track2, 0)) {
            track1pion = true;
            track2kaon = true;
          }
          if (selectionPIDNew(track2, 1) && selectionPIDNew(track1, 0)) {
            track2pion = true;
            track1kaon = true;
          }
        }
        if (!ispTdepPID) {
          if (selectionPID(track1, 1) && selectionPID(track2, 0)) {
            track1pion = true;
            track2kaon = true;
          }
          if (selectionPID(track2, 1) && selectionPID(track1, 0)) {
            track2pion = true;
            track1kaon = true;
          }
        }

        if (track1kaon && track2pion) {
          if (additionalQAplots) {
            histos.fill(HIST("QAafter/TPC_Nsigma_allka"), track1.pt(), track1.tpcNSigmaKa());
            histos.fill(HIST("QAafter/TOF_Nsigma_allka"), track1.pt(), track1.tofNSigmaKa());
            histos.fill(HIST("QAafter/trkDCAxyka"), track1.dcaXY());
            histos.fill(HIST("QAafter/trkDCAzka"), track1.dcaZ());
            histos.fill(HIST("QAafter/TOF_TPC_Mapka_allka"), track1.tofNSigmaKa(), track1.tpcNSigmaKa());
            histos.fill(HIST("QAafter/TOF_TPC_Mapka_allpi"), track2.tofNSigmaPi(), track2.tpcNSigmaPi());
            histos.fill(HIST("QAafter/TPC_Nsigma_allpi"), track2.pt(), track2.tpcNSigmaPi());
            histos.fill(HIST("QAafter/TOF_Nsigma_allpi"), track2.pt(), track2.tofNSigmaPi());
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
        histos.fill(HIST("hpTvsRapidity"), KstarMother.Pt(), KstarMother.Rapidity());
        if (TMath::Abs(KstarMother.Rapidity()) > confRapidity) {
          continue;
        }
        auto phiminuspsi = GetPhiInRange(KstarMother.Phi() - psiFT0C);
        auto v2 = TMath::Cos(2.0 * phiminuspsi);
        histos.fill(HIST("hSparseV2SASameEvent_V2"), KstarMother.M(), KstarMother.Pt(), v2, centrality);

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
            auto phiminuspsiRot = GetPhiInRange(kstarrot.Phi() - psiFT0C);
            auto v2Rot = TMath::Cos(2.0 * phiminuspsiRot);
            histos.fill(HIST("hSparseV2SASameEventRotational_V2"), kstarrot.M(), kstarrot.Pt(), v2Rot, centrality);
          }
        }
      }
    }
  }
  PROCESS_SWITCH(kstarpbpb, processSameEvent, "Process Same event", true);
  void processMixedEvent(EventCandidates const& collisions, TrackCandidates const& tracks)
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
      int occupancy = collision1.trackOccupancyInTimeRange();
      if (occupancy >= 1500) {
        return;
      }
      auto centrality = collision1.centFT0C();
      auto centrality2 = collision2.centFT0C();
      auto psiFT0C = collision1.psiFT0C();
      bool track1pion = false;
      bool track1kaon = false;
      bool track2pion = false;
      bool track2kaon = false;

      if (additionalEvsel && !eventSelected(collision1, tracks.size(), centrality)) {
        continue;
      }
      if (additionalEvsel && !eventSelected(collision2, tracks.size(), centrality2)) {
        continue;
      }

      for (auto& [track1, track2] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {
        if (track1.sign() * track2.sign() > 0) {
          continue;
        }
        if (!selectionTrack(track1) || !selectionTrack(track2)) {
          continue;
        }
        if (ispTdepPID && !(selectionPIDNew(track1, 0) || selectionPIDNew(track1, 1))) {
          continue;
        }
        if (ispTdepPID && !(selectionPIDNew(track2, 1) || selectionPIDNew(track2, 0))) {
          continue;
        }
        if (!ispTdepPID && !(selectionPID(track1, 0) || selectionPID(track1, 1))) {
          continue;
        }
        if (!ispTdepPID && !(selectionPID(track2, 1) || selectionPID(track2, 0))) {
          continue;
        }

        if (ispTdepPID) {
          if (selectionPIDNew(track1, 1) && selectionPIDNew(track2, 0)) {
            track1pion = true;
            track2kaon = true;
          }
          if (selectionPIDNew(track2, 1) && selectionPIDNew(track1, 0)) {
            track2pion = true;
            track1kaon = true;
          }
        }
        if (!ispTdepPID) {
          if (selectionPID(track1, 1) && selectionPID(track2, 0)) {
            track1pion = true;
            track2kaon = true;
          }
          if (selectionPID(track2, 1) && selectionPID(track1, 0)) {
            track2pion = true;
            track1kaon = true;
          }
        }
        if (track1kaon && track2pion) {
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
        auto v2 = TMath::Cos(2.0 * phiminuspsi);
        histos.fill(HIST("hSparseV2SAMixedEvent_V2"), KstarMother.M(), KstarMother.Pt(), v2, centrality);
      }
    }
  }
  PROCESS_SWITCH(kstarpbpb, processMixedEvent, "Process Mixed event", true);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<kstarpbpb>(cfgc, TaskName{"kstarpbpb"})};
}
