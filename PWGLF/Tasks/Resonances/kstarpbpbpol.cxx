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
// polarization task for resonances, dukhishyam.mallick@cern.ch

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
struct kstarpbpbpol {

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
  Configurable<float> nsigmaCutTPC{"nsigmacutTPC", 3.0, "Value of the TPC Nsigma cut"};
  Configurable<float> nsigmaCutCombined{"nsigmaCutCombined", 3.0, "Value of the TOF Nsigma cut"};
  Configurable<int> cfgNoMixedEvents{"cfgNoMixedEvents", 1, "Number of mixed events per event"};
  Configurable<int> cfgITScluster{"cfgITScluster", 0, "Number of ITS cluster"};
  Configurable<bool> ispTdepPID{"ispTdepPID", true, "pT dependent PID"};
  Configurable<double> confRapidity{"confRapidity", 0.5, "Rapidity cut"};

  // AxisSpec axisMult{cCentBins, cCentLow, cCentHigh};
  // AxisSpec axisPt{cpTbins, cpTlow, cpThigh};
  // AxisSpec axisMass{cMassbins, cMasslow, cMasshigh};
  // AxisSpec axisMass{cCosbins, cCoslow, cCoshigh};

  ConfigurableAxis configThnAxisInvMass{"configThnAxisInvMass", {120, 0.66, 1.2}, "#it{M} (GeV/#it{c}^{2})"};
  ConfigurableAxis configThnAxisPt{"configThnAxisPt", {100, 0.0, 10.}, "#it{p}_{T} (GeV/#it{c})"};
  ConfigurableAxis configThnAxisCentrality{"configThnAxisCentrality", {8, 0., 80}, "Centrality"};
  ConfigurableAxis configThnAxisPOL{"configThnAxisPOL", {20, -1.0, 1.0}, "Costheta"};

  // Define axes with variable binning
  // ConfigurableAxis configThnAxisInvMass{"configThnAxisInvMass", {cMassbins,cMasslow,cMasshigh}, "#it{M} (GeV/#it{c}^{2})"}; // need to define variable binning
  // ConfigurableAxis configThnAxisPt{"configThnAxisPt", {cpTbins,cpTlow, cpThigh}, "#it{p}_{T} (GeV/#it{c})"};
  // ConfigurableAxis configThnAxisCentrality{"configThnAxisCentrality", {cCentBins,cCentLow,cCentHigh}, "Centrality"};
  // ConfigurableAxis configThnAxisPOL{"configThnAxisPOL", {cCosbins,cCoslow,cCoshigh}, "Costheta"};

  Configurable<bool> additionalEvsel{"additionalEvsel", false, "Additional event selcection"};
  Configurable<bool> timFrameEvsel{"timFrameEvsel", false, "TPC Time frame boundary cut"};

  /// output THnSparses
  Configurable<bool> activateTHnSparseCosThStarHelicity{"activateTHnSparseCosThStarHelicity", true, "Activate the THnSparse with cosThStar w.r.t. helicity axis"};
  Configurable<bool> activateTHnSparseCosThStarProduction{"activateTHnSparseCosThStarProduction", false, "Activate the THnSparse with cosThStar w.r.t. production axis"};
  Configurable<bool> activateTHnSparseCosThStarBeam{"activateTHnSparseCosThStarBeam", false, "Activate the THnSparse with cosThStar w.r.t. beam axis"};
  Configurable<bool> activateTHnSparseCosThStarRandom{"activateTHnSparseCosThStarRandom", false, "Activate the THnSparse with cosThStar w.r.t. random axis"};

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

    /// check output THnSparses
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

    const AxisSpec thnAxisInvMass{configThnAxisInvMass, "#it{M} (GeV/#it{c}^{2})"};
    const AxisSpec thnAxisPt{configThnAxisPt, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec thnAxisCentrality{configThnAxisCentrality, "Centrality (%)"};
    const AxisSpec thnAxisPOL{configThnAxisPOL, "POL"};
    // const AxisSpec thnAxisV2{configThnAxisV2, "V2"};
    //  AxisSpec phiAxis = {500, -6.28, 6.28, "phi"};
    //  AxisSpec resAxis = {400, -2, 2, "Res"};
    //  AxisSpec centAxis = {8, 0, 80, "V0M (%)"};

    histos.add("hpTvsRapidity", "pT vs Rapidity", kTH2F, {{100, 0.0f, 10.0f}, {300, -1.5f, 1.5f}});
    histos.add("hFTOCvsTPC", "Mult correlation FT0C vs. TPC", kTH2F, {{80, 0.0f, 80.0f}, {100, -0.5f, 5999.5f}});
    histos.add("hFTOCvsTPCSelected", "Mult correlation FT0C vs. TPC after selection", kTH2F, {{80, 0.0f, 80.0f}, {100, -0.5f, 5999.5f}});
    histos.add("hCentrality", "Centrality distribution", kTH1F, {{200, 0.0, 200.0}});
    histos.add("hVtxZ", "Vertex distribution in Z;Z (cm)", kTH1F, {{400, -20.0, 20.0}});

    histos.add("hSparseHESASameEvent", "hSparseHESASameEvent", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisCentrality, thnAxisPOL});
    histos.add("hSparsePPSASameEvent", "hSparsePPSASameEvent", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisCentrality, thnAxisPOL});
    histos.add("hSparseBASASameEvent", "hSparseBASASameEvent", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisCentrality, thnAxisPOL});
    histos.add("hSparseRndASASameEvent", "hSparseRndASASameEvent", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisCentrality, thnAxisPOL});

    histos.add("hSparseHESAMixedEvent", "hSparseHESAMixedEvent", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisCentrality, thnAxisPOL});
    histos.add("hSparsePPSAMixedEvent", "hSparsePPSAMixedEvent", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisCentrality, thnAxisPOL});
    histos.add("hSparseBASAMixedEvent", "hSparseBASAMixedEvent", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisCentrality, thnAxisPOL});
    histos.add("hSparseRndASAMixedEvent", "hSparseRndASAMixedEvent", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisCentrality, thnAxisPOL});

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
  // double massKstar = o2::constants::physics::MassKstar;
  //  double massKstar = 0.892;

  template <typename TCollision>
  bool eventSelected(TCollision collision, const int& multTrk, const float& centrality)
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
    if (!(candidate.isGlobalTrack() && candidate.isPVContributor() && candidate.itsNCls() > cfgITScluster)) {
      return false;
    }
    return true;
  }

  template <typename T>
  bool selectionPIDpTdependent(const T& candidate, int PID)
  {
    if (PID == 0) {
      if (candidate.pt() < 0.5 && TMath::Abs(candidate.tpcNSigmaKa()) < nsigmaCutTPC) {
        return true;
      }
      if (candidate.pt() >= 0.5 && candidate.hasTOF() && ((candidate.tofNSigmaKa() * candidate.tofNSigmaKa()) + (candidate.tpcNSigmaKa() * candidate.tpcNSigmaKa())) < (nsigmaCutCombined * nsigmaCutCombined)) {
        return true;
      }
    } else if (PID == 1) {
      if (candidate.pt() < 0.5 && TMath::Abs(candidate.tpcNSigmaPi()) < nsigmaCutTPC) {
        return true;
      }
      if (candidate.pt() >= 0.5 && candidate.hasTOF() && ((candidate.tofNSigmaPi() * candidate.tofNSigmaPi()) + (candidate.tpcNSigmaPi() * candidate.tpcNSigmaPi())) < (nsigmaCutCombined * nsigmaCutCombined)) {
        return true;
      }
    }
    return false;
  }

  template <typename T>
  bool selectionPID(const T& candidate, int PID)
  {
    if (PID == 0) {
      if (!candidate.hasTOF() && TMath::Abs(candidate.tpcNSigmaKa()) < nsigmaCutTPC) {
        return true;
      }
      if (candidate.hasTOF() && ((candidate.tofNSigmaKa() * candidate.tofNSigmaKa()) + (candidate.tpcNSigmaKa() * candidate.tpcNSigmaKa())) < (nsigmaCutCombined * nsigmaCutCombined)) {
        return true;
      }
    } else if (PID == 1) {
      if (!candidate.hasTOF() && TMath::Abs(candidate.tpcNSigmaPi()) < nsigmaCutTPC) {
        return true;
      }
      if (candidate.hasTOF() && ((candidate.tofNSigmaPi() * candidate.tofNSigmaPi()) + (candidate.tpcNSigmaPi() * candidate.tpcNSigmaPi())) < (nsigmaCutCombined * nsigmaCutCombined)) {
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
  // Here also need to define different mixing criteria
  ConfigurableAxis axisVertex{"axisVertex", {20, -10, 10}, "vertex axis for bin"};
  ConfigurableAxis axisMultiplicityClass{"axisMultiplicityClass", {20, 0, 100}, "multiplicity percentile for bin"};
  ConfigurableAxis axisEPAngle{"axisEPAngle", {6, -TMath::Pi() / 2, TMath::Pi() / 2}, "event plane angle"};

  using BinningTypeVertexContributor = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0C, aod::epcalibrationtable::PsiFT0C>;
  ROOT::Math::PxPyPzMVector KstarMother, daughter1, daughter2;

  void processSameEvent(EventCandidates::iterator const& collision, TrackCandidates const& tracks, aod::BCs const&)
  {
    if (!collision.sel8()) {
      return;
    }

    if (!collision.triggereventep()) {
      return;
    }
    if (timFrameEvsel && (!collision.selection_bit(aod::evsel::kNoTimeFrameBorder) || !collision.selection_bit(aod::evsel::kNoITSROFrameBorder))) {
      return;
    }
    auto posThisColl = posTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto negThisColl = negTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto centrality = collision.centFT0C();
    auto multTPC = collision.multNTracksPV();
    histos.fill(HIST("hFTOCvsTPC"), centrality, multTPC);
    if (additionalEvsel && !eventSelected(collision, tracks.size(), centrality)) {
      return;
    }
    histos.fill(HIST("hFTOCvsTPCSelected"), centrality, multTPC);
    histos.fill(HIST("hCentrality"), centrality);
    histos.fill(HIST("hVtxZ"), collision.posZ());
    for (auto track1 : posThisColl) {
      if (!selectionTrack(track1)) {
        continue;
      }
      bool track1pion = false;
      bool track1kaon = false;
      if (ispTdepPID && !(selectionPIDpTdependent(track1, 0) || selectionPIDpTdependent(track1, 1))) {
        continue;
      }
      for (auto track2 : negThisColl) {
        bool track2pion = false;
        bool track2kaon = false;
        if (ispTdepPID && !(selectionPIDpTdependent(track2, 0) || selectionPIDpTdependent(track2, 1))) {
          continue;
        }
        // track selection
        if (!selectionTrack(track2)) {
          continue;
        }
        if (track1.sign() * track2.sign() > 0) {
          continue;
        }
        if (selectionPIDpTdependent(track1, 1)) {
          track1pion = true;
        } else if (selectionPIDpTdependent(track1, 0)) {
          track1kaon = true;
        }
        if (selectionPIDpTdependent(track2, 1)) {
          track2pion = true;
        } else if (selectionPIDpTdependent(track2, 0)) {
          track2kaon = true;
        }
        if (track1pion && track2kaon) {
          daughter1 = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), massPi);
          daughter2 = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), massKa);
        } else if (track1kaon && track2pion) {
          daughter1 = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), massKa);
          daughter2 = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), massPi);
        } else {
          continue;
        }
        KstarMother = daughter1 + daughter2;
        histos.fill(HIST("hpTvsRapidity"), KstarMother.Pt(), KstarMother.Rapidity());

        if (TMath::Abs(KstarMother.Rapidity()) > confRapidity) {
          continue;
        }

        float phiRandom = gRandom->Uniform(0.f, constants::math::TwoPI);
        float thetaRandom = gRandom->Uniform(0.f, constants::math::PI);

        ROOT::Math::PxPyPzMVector fourVecDau = ROOT::Math::PxPyPzMVector(daughter1.px(), daughter1.py(), daughter1.pz(), massKa);
        ROOT::Math::PxPyPzMVector fourVecMother = ROOT::Math::PxPyPzMVector(KstarMother.px(), KstarMother.py(), KstarMother.pz(), KstarMother.M());
        ROOT::Math::Boost boost{fourVecMother.BoostToCM()};
        ROOT::Math::PxPyPzMVector fourVecDauCM = boost(fourVecDau);
        ROOT::Math::XYZVector threeVecDauCM = fourVecDauCM.Vect();

        if (activateTHnSparseCosThStarHelicity) {
          ROOT::Math::XYZVector helicityVec = fourVecMother.Vect();
          // ROOT::Math::XYZVectorF zaxis_HE{helicityVec.Unit()};
          // ROOT::Math::XYZVectorF v1_CM{threeVecDauCM.Unit()};
          // float cosThetaStarHelicity = zaxis_HE.Dot(v1_CM);

          float cosThetaStarHelicity = helicityVec.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2()) / std::sqrt(helicityVec.Mag2());
          histos.fill(HIST("hSparseHESASameEvent"), KstarMother.M(), KstarMother.Pt(), centrality, cosThetaStarHelicity);
        }

        if (activateTHnSparseCosThStarProduction) {
          ROOT::Math::XYZVector normalVec = ROOT::Math::XYZVector(KstarMother.py(), -KstarMother.px(), 0.f);
          float cosThetaStarProduction = normalVec.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2()) / std::sqrt(normalVec.Mag2());
          histos.fill(HIST("hSparsePPSASameEvent"), KstarMother.M(), KstarMother.Pt(), centrality, cosThetaStarProduction);
        }
        if (activateTHnSparseCosThStarBeam) {
          ROOT::Math::XYZVector beamVec = ROOT::Math::XYZVector(0.f, 0.f, 1.f);
          float cosThetaStarBeam = beamVec.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2());
          histos.fill(HIST("hSparseBASASameEvent"), KstarMother.M(), KstarMother.Pt(), centrality, cosThetaStarBeam);
        }
        if (activateTHnSparseCosThStarRandom) {
          ROOT::Math::XYZVector randomVec = ROOT::Math::XYZVector(std::sin(thetaRandom) * std::cos(phiRandom), std::sin(thetaRandom) * std::sin(phiRandom), std::cos(thetaRandom));
          float cosThetaStarRandom = randomVec.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2());
          histos.fill(HIST("hSparseRndASASameEvent"), KstarMother.M(), KstarMother.Pt(), centrality, cosThetaStarRandom);
        }
        // auto phiminuspsi = GetPhiInRange(KstarMother.Phi() - psiFT0C);
        // auto v2 = TMath::Cos(2.0 * phiminuspsi);
        // histos.fill(HIST("hSparseV2SASameEvent_V2"), KstarMother.M(), KstarMother.Pt(), v2, centrality);
      }
    }
  }
  PROCESS_SWITCH(kstarpbpbpol, processSameEvent, "Process Same event", true);
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
      auto centrality = collision1.centFT0C();
      auto centrality2 = collision2.centFT0C();
      // auto psiFT0C = collision1.psiFT0C();
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
        if (ispTdepPID && !(selectionPIDpTdependent(track1, 0) || selectionPIDpTdependent(track1, 1))) {
          continue;
        }
        if (ispTdepPID && !(selectionPIDpTdependent(track2, 1) || selectionPIDpTdependent(track2, 0))) {
          continue;
        }

        if (selectionPIDpTdependent(track1, 1)) {
          track1pion = true;
        } else if (selectionPIDpTdependent(track1, 0)) {
          track1kaon = true;
        }
        if (selectionPIDpTdependent(track2, 1)) {
          track2pion = true;
        } else if (selectionPIDpTdependent(track2, 0)) {
          track2kaon = true;
        }

        if (track1pion && track2kaon) {
          daughter1 = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), massPi);
          daughter2 = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), massKa);
        } else if (track1kaon && track2pion) {
          daughter1 = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), massKa);
          daughter2 = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), massPi);
        } else {
          continue;
        }
        KstarMother = daughter1 + daughter2;
        if (TMath::Abs(KstarMother.Rapidity()) > confRapidity) {
          continue;
        }

        float phiRandom = gRandom->Uniform(0.f, constants::math::TwoPI);
        float thetaRandom = gRandom->Uniform(0.f, constants::math::PI);

        ROOT::Math::PxPyPzMVector fourVecDau = ROOT::Math::PxPyPzMVector(daughter1.px(), daughter1.py(), daughter1.pz(), massKa);
        ROOT::Math::PxPyPzMVector fourVecMother = ROOT::Math::PxPyPzMVector(KstarMother.px(), KstarMother.py(), KstarMother.pz(), KstarMother.M());
        ROOT::Math::Boost boost{fourVecMother.BoostToCM()};
        ROOT::Math::PxPyPzMVector fourVecDauCM = boost(fourVecDau);
        ROOT::Math::XYZVector threeVecDauCM = fourVecDauCM.Vect();

        if (activateTHnSparseCosThStarHelicity) {
          ROOT::Math::XYZVector helicityVec = fourVecMother.Vect();
          float cosThetaStarHelicity = helicityVec.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2()) / std::sqrt(helicityVec.Mag2());
          histos.fill(HIST("hSparseHESAMixedEvent"), KstarMother.M(), KstarMother.Pt(), centrality, cosThetaStarHelicity);
        }

        if (activateTHnSparseCosThStarProduction) {
          ROOT::Math::XYZVector normalVec = ROOT::Math::XYZVector(KstarMother.py(), -KstarMother.px(), 0.f);
          float cosThetaStarProduction = normalVec.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2()) / std::sqrt(normalVec.Mag2());
          histos.fill(HIST("hSparsePPSAMixedEvent"), KstarMother.M(), KstarMother.Pt(), centrality, cosThetaStarProduction);
        }
        if (activateTHnSparseCosThStarBeam) {
          ROOT::Math::XYZVector beamVec = ROOT::Math::XYZVector(0.f, 0.f, 1.f);
          float cosThetaStarBeam = beamVec.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2());
          histos.fill(HIST("hSparseBASAMixedEvent"), KstarMother.M(), KstarMother.Pt(), centrality, cosThetaStarBeam);
        }
        if (activateTHnSparseCosThStarRandom) {
          ROOT::Math::XYZVector randomVec = ROOT::Math::XYZVector(std::sin(thetaRandom) * std::cos(phiRandom), std::sin(thetaRandom) * std::sin(phiRandom), std::cos(thetaRandom));
          float cosThetaStarRandom = randomVec.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2());
          histos.fill(HIST("hSparseRndASAMixedEvent"), KstarMother.M(), KstarMother.Pt(), centrality, cosThetaStarRandom);
        }
      }
    }
  }
  PROCESS_SWITCH(kstarpbpbpol, processMixedEvent, "Process Mixed event", true);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<kstarpbpbpol>(cfgc, TaskName{"kstarpbpbpol"})};
}
