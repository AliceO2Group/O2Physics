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
// Particle flow task
// prottay.das@cern.ch

#include "Math/GenVector/Boost.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "TF1.h"
#include "TRandom3.h"
#include <TDatabasePDG.h>
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

// #include "Common/DataModel/Qvectors.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/SPCalibrationTables.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/FT0Corrected.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "CommonConstants/PhysicsConstants.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/StepTHn.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;

struct lambdav2 {

  int mRunNumber;
  int multEstimator;
  float d_bz;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Service<o2::framework::O2DatabasePDG> pdg;

  // fill output
  Configurable<bool> additionalEvSel{"additionalEvSel", false, "additionalEvSel"};
  Configurable<bool> additionalEvSel2{"additionalEvSel2", false, "additionalEvSel2"};
  Configurable<bool> additionalEvSel3{"additionalEvSel3", false, "additionalEvSel3"};
  Configurable<bool> correction1{"correction1", false, "fill histograms including corrections 1"};
  Configurable<bool> correction2{"correction2", false, "fill histograms including corrections 2"};
  Configurable<bool> QA{"QA", false, "flag for QA"};
  Configurable<bool> mycut{"mycut", false, "select tracks based on my cuts"};
  Configurable<bool> tofhit{"tofhit", true, "select tracks based on tof hit"};
  Configurable<bool> globalpt{"globalpt", true, "select tracks based on pt global vs tpc"};
  Configurable<int> useprofile{"useprofile", 3, "flag to select profile vs Sparse"};
  Configurable<int> QxyNbins{"QxyNbins", 100, "Number of bins in QxQy histograms"};
  Configurable<float> lbinQxy{"lbinQxy", -5.0, "lower bin value in QxQy histograms"};
  Configurable<float> hbinQxy{"hbinQxy", 5.0, "higher bin value in QxQy histograms"};
  Configurable<int> cfgMaxOccupancy{"cfgMaxOccupancy", 1000, "maximum occupancy of tracks in neighbouring collisions in a given time range"};
  Configurable<int> cfgMinOccupancy{"cfgMinOccupancy", 0, "maximum occupancy of tracks in neighbouring collisions in a given time range"};

  // events
  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};
  Configurable<float> cfgCutCentralityMax{"cfgCutCentralityMax", 50.0f, "Accepted maximum Centrality"};
  Configurable<float> cfgCutCentralityMin{"cfgCutCentralityMin", 30.0f, "Accepted minimum Centrality"};
  // proton track cut
  Configurable<float> confRapidity{"confRapidity", 0.8, "cut on Rapidity"};
  Configurable<float> cfgCutPT{"cfgCutPT", 0.15, "PT cut on daughter track"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8, "Eta cut on daughter track"};
  Configurable<float> cfgCutDCAxy{"cfgCutDCAxy", 0.1f, "DCAxy range for tracks"};
  Configurable<float> cfgCutDCAz{"cfgCutDCAz", 0.1f, "DCAz range for tracks"};
  Configurable<int> cfgITScluster{"cfgITScluster", 5, "Number of ITS cluster"};
  Configurable<int> cfgTPCcluster{"cfgTPCcluster", 70, "Number of TPC cluster"};
  Configurable<bool> isPVContributor{"isPVContributor", true, "is PV contributor"};
  Configurable<bool> checkwithpub{"checkwithpub", true, "checking results with published"};
  Configurable<float> nsigmaCutTPCPi{"nsigmaCutTPCPi", 3, "PID selections for Pions"};
  Configurable<float> nsigmaCutTPCKa{"nsigmaCutTPCKa", 3, "PID selections for Kaons"};
  Configurable<float> nsigmaCutTPCPr{"nsigmaCutTPCPr", 3, "PID selections for Protons"};
  Configurable<float> nsigmaCutTOFPi{"nsigmaCutTOFPi", 3, "PID selections for TOF Pions"};
  Configurable<float> nsigmaCutTOFKa{"nsigmaCutTOFKa", 3, "PID selections for TOF Kaons"};
  Configurable<float> nsigmaCutTOFPr{"nsigmaCutTOFPr", 3, "PID selections for TOF Protons"};
  Configurable<float> cfgCutTOFBeta{"cfgCutTOFBeta", 0, "Beta selections for Particles"};

  Configurable<int> CentNbins{"CentNbins", 16, "Number of bins in cent histograms"};
  Configurable<float> lbinCent{"lbinCent", 0.0, "lower bin value in cent histograms"};
  Configurable<float> hbinCent{"hbinCent", 80.0, "higher bin value in cent histograms"};
  Configurable<int> SANbins{"SANbins", 20, "Number of bins in costhetastar"};
  Configurable<float> lbinSA{"lbinSA", -1.0, "lower bin value in costhetastar histograms"};
  Configurable<float> hbinSA{"hbinSA", 1.0, "higher bin value in costhetastar histograms"};
  Configurable<int> PolNbins{"PolNbins", 20, "Number of bins in polarisation"};
  Configurable<float> lbinPol{"lbinPol", -1.0, "lower bin value in #phi-#psi histograms"};
  Configurable<float> hbinPol{"hbinPol", 1.0, "higher bin value in #phi-#psi histograms"};
  Configurable<int> IMNbins{"IMNbins", 100, "Number of bins in invariant mass"};
  Configurable<float> lbinIM{"lbinIM", 1.0, "lower bin value in IM histograms"};
  Configurable<float> hbinIM{"hbinIM", 1.2, "higher bin value in IM histograms"};
  Configurable<int> ptNbins{"ptNbins", 50, "Number of bins in pt"};
  Configurable<float> lbinpt{"lbinpt", 0.0, "lower bin value in pt histograms"};
  Configurable<float> hbinpt{"hbinpt", 10.0, "higher bin value in pt histograms"};
  Configurable<int> resNbins{"resNbins", 50, "Number of bins in reso"};
  Configurable<float> lbinres{"lbinres", 0.0, "lower bin value in reso histograms"};
  Configurable<float> hbinres{"hbinres", 10.0, "higher bin value in reso histograms"};
  Configurable<int> etaNbins{"etaNbins", 20, "Number of bins in eta"};
  Configurable<float> lbineta{"lbineta", -1.0, "lower bin value in eta histograms"};
  Configurable<float> hbineta{"hbineta", 1.0, "higher bin value in eta histograms"};
  Configurable<int> spNbins{"spNbins", 2000, "Number of bins in sp"};
  Configurable<float> lbinsp{"lbinsp", -1.0, "lower bin value in sp histograms"};
  Configurable<float> hbinsp{"hbinsp", 1.0, "higher bin value in sp histograms"};
  Configurable<int> phiNbins{"phiNbins", 30, "Number of bins in phi"};

  ConfigurableAxis configcentAxis{"configcentAxis", {VARIABLE_WIDTH, 0.0, 10.0, 40.0, 80.0}, "Cent V0M"};
  ConfigurableAxis configthnAxispT{"configthnAxisPt", {VARIABLE_WIDTH, 0.2, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.5, 8.0, 10.0, 100.0}, "#it{p}_{T} (GeV/#it{c})"};
  ConfigurableAxis configetaAxis{"configetaAxis", {VARIABLE_WIDTH, -0.8, -0.4, -0.2, 0, 0.2, 0.4, 0.8}, "Eta"};
  ConfigurableAxis configthnAxisPol{"configthnAxisPol", {VARIABLE_WIDTH, -1.0, -0.6, -0.2, 0, 0.2, 0.4, 0.8}, "Pol"};
  ConfigurableAxis configphiAxis{"configphiAxis", {VARIABLE_WIDTH, 0.0, 0.2, 0.4, 0.8, 1.0, 2.0, 2.5, 3.0, 4.0, 5.0, 5.5, 6.28}, "PhiAxis"};

  SliceCache cache;
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(o2::framework::InitContext&)
  {
    AxisSpec thnAxispT{ptNbins, lbinpt, hbinpt, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec thnAxisRap{100, -0.5, 0.5, "Rapidity"};
    AxisSpec thnAxisres{resNbins, lbinres, hbinres, "Reso"};
    AxisSpec thnAxisInvMass{IMNbins, lbinIM, hbinIM, "#it{M} (GeV/#it{c}^{2})"};
    AxisSpec thnAxisPol{PolNbins, lbinPol, hbinPol, "Sin(#phi - #psi)"};
    AxisSpec thnAxisCosThetaStar{SANbins, lbinSA, hbinSA, "SA"};
    AxisSpec centAxis = {CentNbins, lbinCent, hbinCent, "V0M (%)"};
    AxisSpec etaAxis = {etaNbins, lbineta, hbineta, "Eta"};
    AxisSpec spAxis = {spNbins, lbinsp, hbinsp, "Sp"};
    AxisSpec qxZDCAxis = {QxyNbins, lbinQxy, hbinQxy, "Qx"};
    AxisSpec phiAxis = {phiNbins, 0.0, 6.28, "phi-phiStar"};
    /*
      histos.add("hpuxQxpvscentpteta", "hpuxQxpvscentpteta", HistType::kTHnSparseF, {centAxis, thnAxispT, etaAxis, spAxis}, true);
      histos.add("hpuyQypvscentpteta", "hpuyQypvscentpteta", HistType::kTHnSparseF, {centAxis, thnAxispT, etaAxis, spAxis}, true);
      histos.add("hpuxQxtvscentpteta", "hpuxQxtvscentpteta", HistType::kTHnSparseF, {centAxis, thnAxispT, etaAxis, spAxis}, true);
      histos.add("hpuyQytvscentpteta", "hpuyQytvscentpteta", HistType::kTHnSparseF, {centAxis, thnAxispT, etaAxis, spAxis}, true);
      histos.add("hpuxyQxytvscentpteta", "hpuxyQxytvscentpteta", HistType::kTHnSparseF, {centAxis, thnAxispT, etaAxis, spAxis}, true);
      histos.add("hpuxyQxypvscentpteta", "hpuxyQxypvscentpteta", HistType::kTHnSparseF, {centAxis, thnAxispT, etaAxis, spAxis}, true);*/
    histos.add("hpoddv1vscentpteta", "hpoddv1vscentpteta", HistType::kTHnSparseF, {centAxis, thnAxispT, etaAxis, spAxis}, true);
    histos.add("hpevenv1vscentpteta", "hpevenv1vscentpteta", HistType::kTHnSparseF, {centAxis, thnAxispT, etaAxis, spAxis}, true);
    histos.add("hpoddv1vscentptetakaon", "hpoddv1vscentptetakaon", HistType::kTHnSparseF, {centAxis, thnAxispT, etaAxis, spAxis}, true);
    histos.add("hpevenv1vscentptetakaon", "hpevenv1vscentptetakaon", HistType::kTHnSparseF, {centAxis, thnAxispT, etaAxis, spAxis}, true);
    histos.add("hpoddv1vscentptetaproton", "hpoddv1vscentptetaproton", HistType::kTHnSparseF, {centAxis, thnAxispT, etaAxis, spAxis}, true);
    histos.add("hpevenv1vscentptetaproton", "hpevenv1vscentptetaproton", HistType::kTHnSparseF, {centAxis, thnAxispT, etaAxis, spAxis}, true);

    /* histos.add("hpuxQxpvscentptetaneg", "hpuxQxpvscentptetaneg", HistType::kTHnSparseF, {centAxis, thnAxispT, etaAxis, spAxis}, true);
       histos.add("hpuyQypvscentptetaneg", "hpuyQypvscentptetaneg", HistType::kTHnSparseF, {centAxis, thnAxispT, etaAxis, spAxis}, true);
       histos.add("hpuxQxtvscentptetaneg", "hpuxQxtvscentptetaneg", HistType::kTHnSparseF, {centAxis, thnAxispT, etaAxis, spAxis}, true);
       histos.add("hpuyQytvscentptetaneg", "hpuyQytvscentptetaneg", HistType::kTHnSparseF, {centAxis, thnAxispT, etaAxis, spAxis}, true);
       histos.add("hpuxyQxytvscentptetaneg", "hpuxyQxytvscentptetaneg", HistType::kTHnSparseF, {centAxis, thnAxispT, etaAxis, spAxis}, true);
       histos.add("hpuxyQxypvscentptetaneg", "hpuxyQxypvscentptetaneg", HistType::kTHnSparseF, {centAxis, thnAxispT, etaAxis, spAxis}, true);*/
    histos.add("hpoddv1vscentptetaneg", "hpoddv1vscentptetaneg", HistType::kTHnSparseF, {centAxis, thnAxispT, etaAxis, spAxis}, true);
    histos.add("hpevenv1vscentptetaneg", "hpevenv1vscentptetaneg", HistType::kTHnSparseF, {centAxis, thnAxispT, etaAxis, spAxis}, true);
    histos.add("hpoddv1vscentptetakaonneg", "hpoddv1vscentptetakaonneg", HistType::kTHnSparseF, {centAxis, thnAxispT, etaAxis, spAxis}, true);
    histos.add("hpevenv1vscentptetakaonneg", "hpevenv1vscentptetakaonneg", HistType::kTHnSparseF, {centAxis, thnAxispT, etaAxis, spAxis}, true);
    histos.add("hpoddv1vscentptetaprotonneg", "hpoddv1vscentptetaprotonneg", HistType::kTHnSparseF, {centAxis, thnAxispT, etaAxis, spAxis}, true);
    histos.add("hpevenv1vscentptetaprotonneg", "hpevenv1vscentptetaprotonneg", HistType::kTHnSparseF, {centAxis, thnAxispT, etaAxis, spAxis}, true);

    histos.add("hpQxtQxpvscent", "hpQxtQxpvscent", HistType::kTHnSparseF, {centAxis, spAxis}, true);
    histos.add("hpQytQypvscent", "hpQytQypvscent", HistType::kTHnSparseF, {centAxis, spAxis}, true);
    histos.add("hpQxytpvscent", "hpQxytpvscent", HistType::kTHnSparseF, {centAxis, spAxis}, true);
    histos.add("hpQxtQypvscent", "hpQxtQypvscent", HistType::kTHnSparseF, {centAxis, spAxis}, true);
    histos.add("hpQxpQytvscent", "hpQxpQytvscent", HistType::kTHnSparseF, {centAxis, spAxis}, true);

    histos.add("hCentrality", "Centrality distribution", kTH1F, {{centAxis}});
  }

  template <typename T>
  bool selectionTrack(const T& candidate)
  {

    if (mycut) {
      if (!candidate.isGlobalTrack() || !candidate.isPVContributor() || !(candidate.itsNCls() > cfgITScluster) || !(candidate.tpcNClsFound() > cfgTPCcluster) || !(candidate.itsNClsInnerBarrel() >= 1)) {
        return false;
      }
    } else {
      if (!(candidate.isGlobalTrack() && candidate.isPVContributor() && candidate.itsNCls() > cfgITScluster && candidate.tpcNClsFound() > cfgTPCcluster && candidate.itsNClsInnerBarrel() >= 1)) {
        return false;
      }
    }
    return true;
  }

  double GetPhiInRange(double phi)
  {
    double result = phi;
    while (result < 0) {
      result = result + 2. * TMath::Pi();
    }
    while (result > 2. * TMath::Pi()) {
      result = result - 2. * TMath::Pi();
    }
    return result;
  }

  template <typename T>
  bool SelectionPID(const T& candidate, int PID)
  {
    if (PID == 0) // pion
    {
      auto combPIDPi = TMath::Sqrt(TMath::Abs(candidate.tofNSigmaPi() * candidate.tofNSigmaPi() + candidate.tpcNSigmaPi() * candidate.tpcNSigmaPi()));
      if (!candidate.hasTOF() && candidate.tpcInnerParam() < 0.6 && TMath::Abs(candidate.tpcNSigmaPi()) < nsigmaCutTPCPi) {
        return true;
      }
      if (candidate.hasTOF() && candidate.beta() > cfgCutTOFBeta && combPIDPi < nsigmaCutTOFPi) {
        return true;
      }
    } else if (PID == 1) // kaon
    {
      auto combPIDKa = TMath::Sqrt(TMath::Abs(candidate.tofNSigmaKa() * candidate.tofNSigmaKa() + candidate.tpcNSigmaKa() * candidate.tpcNSigmaKa()));
      if (!candidate.hasTOF() && candidate.tpcInnerParam() < 0.45 && TMath::Abs(candidate.tpcNSigmaKa()) < nsigmaCutTPCKa) {
        return true;
      }
      if (candidate.hasTOF() && candidate.beta() > cfgCutTOFBeta && combPIDKa < nsigmaCutTOFKa) {
        return true;
      }
    } else // proton
    {
      auto combPIDPr = TMath::Sqrt(TMath::Abs(candidate.tofNSigmaPr() * candidate.tofNSigmaPr() + candidate.tpcNSigmaPr() * candidate.tpcNSigmaPr()));
      if (!candidate.hasTOF() && candidate.tpcInnerParam() < 0.6 && TMath::Abs(candidate.tpcNSigmaPr()) < nsigmaCutTPCPr) {
        return true;
      }
      if (candidate.hasTOF() && candidate.beta() > cfgCutTOFBeta && combPIDPr < nsigmaCutTOFPr) {
        return true;
      }
    }
    return false;
  }

  ROOT::Math::PxPyPzMVector Lambda, Proton, Pion, fourVecDauCM;
  // ROOT::Math::XYZVector threeVecDauCM, threeVecDauCMXY, eventplaneVec, eventplaneVecNorm, beamvector;
  ROOT::Math::XYZVector threeVecDauCM, threeVecDauCMXY;
  double phiangle = 0.0;
  double massPi = o2::constants::physics::MassPionCharged;
  double massKa = o2::constants::physics::MassKaonCharged;
  double massPr = o2::constants::physics::MassProton;

  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter centralityFilter = (nabs(aod::cent::centFT0C) < cfgCutCentralityMax && nabs(aod::cent::centFT0C) > cfgCutCentralityMin);
  Filter acceptanceFilter = (nabs(aod::track::eta) < cfgCutEta && nabs(aod::track::pt) > cfgCutPT);
  Filter dcaCutFilter = (nabs(aod::track::dcaXY) < cfgCutDCAxy) && (nabs(aod::track::dcaZ) < cfgCutDCAz);

  using EventCandidates = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::FV0Mults, aod::TPCMults, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0Cs, aod::CentFT0As, aod::SPCalibrationTables, aod::Mults>>;
  // using AllTrackCandidates = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullPi, aod::pidTPCFullPr>;
  using AllTrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullPi, aod::pidTPCFullPr, aod::pidTPCFullKa, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFbeta>>;
  using ResoV0s = aod::V0Datas;

  // void processData(EventCandidates::iterator const& collision, AllTrackCandidates const&, ResoV0s const& V0s, aod::BCs const&)
  void processData(EventCandidates::iterator const& collision, AllTrackCandidates const& tracks, ResoV0s const& /*V0s*/, aod::BCs const&)
  {

    if (!collision.sel8()) {
      return;
    }
    auto centrality = collision.centFT0C();
    // histos.fill(HIST("hCentrality0"), centrality);
    if (!collision.triggereventsp()) {
      return;
    }
    // histos.fill(HIST("hCentrality1"), centrality);

    if (additionalEvSel && (!collision.selection_bit(aod::evsel::kNoSameBunchPileup) || !collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV))) {
      return;
    }
    // histos.fill(HIST("hCentrality2"), centrality);
    //  if (additionalEvSel2 && (!collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard))) {
    if (additionalEvSel2 && (collision.trackOccupancyInTimeRange() > cfgMaxOccupancy || collision.trackOccupancyInTimeRange() < cfgMinOccupancy)) {
      return;
    }
    // histos.fill(HIST("hCentrality3"), centrality);
    if (additionalEvSel3 && (!collision.selection_bit(aod::evsel::kNoTimeFrameBorder) || !collision.selection_bit(aod::evsel::kNoITSROFrameBorder))) {
      return;
    }

    auto qxZDCA = collision.qxZDCA();
    auto qxZDCC = collision.qxZDCC();
    auto qyZDCA = collision.qyZDCA();
    auto qyZDCC = collision.qyZDCC();
    // auto psiZDCC = collision.psiZDCC();
    // auto psiZDCA = collision.psiZDCA();

    histos.fill(HIST("hCentrality"), centrality);

    auto QxtQxp = qxZDCA * qxZDCC;
    auto QytQyp = qyZDCA * qyZDCC;
    auto Qxytp = QxtQxp + QytQyp;
    auto QxpQyt = qxZDCA * qyZDCC;
    auto QxtQyp = qxZDCC * qyZDCA;

    histos.fill(HIST("hpQxtQxpvscent"), centrality, QxtQxp);
    histos.fill(HIST("hpQytQypvscent"), centrality, QytQyp);
    histos.fill(HIST("hpQxytpvscent"), centrality, Qxytp);
    histos.fill(HIST("hpQxpQytvscent"), centrality, QxpQyt);
    histos.fill(HIST("hpQxtQypvscent"), centrality, QxtQyp);

    for (auto track : tracks) {
      if (!selectionTrack(track)) {
        continue;
      }

      bool ispion = 0;
      bool iskaon = 0;
      bool isproton = 0;

      if (SelectionPID(track, 0))
        ispion = 1;
      if (SelectionPID(track, 1))
        iskaon = 1;
      if (SelectionPID(track, 2))
        isproton = 1;

      if (ispion && iskaon)
        continue;
      if (ispion && isproton)
        continue;
      if (iskaon && isproton)
        continue;

      float sign = track.sign();
      if (sign == 0.0) // removing neutral particles
        continue;

      auto ux = TMath::Cos(GetPhiInRange(track.phi()));
      auto uy = TMath::Sin(GetPhiInRange(track.phi()));

      // auto uxQxp = ux * qxZDCA;
      // auto uyQyp = uy * qyZDCA;
      // auto uxyQxyp = uxQxp + uyQyp;
      // auto uxQxt = ux * qxZDCC;
      // auto uyQyt = uy * qyZDCC;
      // auto uxyQxyt = uxQxt + uyQyt;
      auto oddv1 = ux * (qxZDCA - qxZDCC) + uy * (qyZDCA - qyZDCC);
      auto evenv1 = ux * (qxZDCA + qxZDCC) + uy * (qyZDCA + qyZDCC);

      if (sign > 0) {
        if (ispion) {
          histos.fill(HIST("hpoddv1vscentpteta"), centrality, track.pt(), track.rapidity(massPi), oddv1);
          histos.fill(HIST("hpevenv1vscentpteta"), centrality, track.pt(), track.rapidity(massPi), evenv1);
        } else if (iskaon) {
          histos.fill(HIST("hpoddv1vscentptetakaon"), centrality, track.pt(), track.rapidity(massKa), oddv1);
          histos.fill(HIST("hpevenv1vscentptetakaon"), centrality, track.pt(), track.rapidity(massKa), evenv1);
        } else if (isproton) {
          histos.fill(HIST("hpoddv1vscentptetaproton"), centrality, track.pt(), track.rapidity(massPr), oddv1);
          histos.fill(HIST("hpevenv1vscentptetaproton"), centrality, track.pt(), track.rapidity(massPr), evenv1);
        }

      } else {
        if (ispion) {
          histos.fill(HIST("hpoddv1vscentptetaneg"), centrality, track.pt(), track.rapidity(massPi), oddv1);
          histos.fill(HIST("hpevenv1vscentptetaneg"), centrality, track.pt(), track.rapidity(massPi), evenv1);
        } else if (iskaon) {
          histos.fill(HIST("hpoddv1vscentptetakaonneg"), centrality, track.pt(), track.rapidity(massKa), oddv1);
          histos.fill(HIST("hpevenv1vscentptetakaonneg"), centrality, track.pt(), track.rapidity(massKa), evenv1);
        } else if (isproton) {
          histos.fill(HIST("hpoddv1vscentptetaprotonneg"), centrality, track.pt(), track.rapidity(massPr), oddv1);
          histos.fill(HIST("hpevenv1vscentptetaprotonneg"), centrality, track.pt(), track.rapidity(massPr), evenv1);
        }
      }
    }
  }
  PROCESS_SWITCH(lambdav2, processData, "Process data", true);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<lambdav2>(cfgc, TaskName{"lambdav2"})};
}
