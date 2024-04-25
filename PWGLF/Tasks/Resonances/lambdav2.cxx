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
/// \brief this is a code for flow of lambda baryons for exotic pheno
/// \author prottay das
/// \since 25/04/2024

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
#include "TF1.h"

#include <array>
#include <cmath>
#include <cstdlib>

#include "CCDB/BasicCCDBManager.h"
#include "CCDB/CcdbApi.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "CommonConstants/PhysicsConstants.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/StepTHn.h"
#include "Framework/runDataProcessing.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "ReconstructionDataFormats/Track.h"

#include "PWGLF/DataModel/EPCalibrationTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;

struct lambdav2 {

  // Connect to ccdb
  Service<ccdb::BasicCCDBManager> ccdb;
  Configurable<int64_t> nolaterthan{
    "ccdb-no-later-than",
    std::chrono::duration_cast<std::chrono::milliseconds>(
      std::chrono::system_clock::now().time_since_epoch())
      .count(),
    "latest acceptable timestamp of creation for the object"};
  Configurable<std::string> url{"ccdb-url", "http://ccdb-test.cern.ch:8080",
                                "url of the ccdb repository"};

  SliceCache cache;

  // Histograms are defined with HistogramRegistry
  HistogramRegistry rEventSelection{"eventSelection",
                                    {},
                                    OutputObjHandlingPolicy::AnalysisObject,
                                    true,
                                    true};
  HistogramRegistry histos{
    "histos",
    {},
    OutputObjHandlingPolicy::AnalysisObject,
    true,
    true};

  // Configurable for histograms
  Configurable<int> nBins{"nBins", 100, "N bins in all histos"};
  // Confugrable for QA histograms
  Configurable<bool> QA{"QA", true, "QA"};
  Configurable<bool> QAv0{"QAv0", false, "QAv0"};
  // Configurable for event selection
  Configurable<float> cutzvertex{"cutzvertex", 10.0f,
                                 "Accepted z-vertex range (cm)"};
  // Configurable parameters for V0 selection
  Configurable<float> ConfV0PtMin{"ConfV0PtMin", 0.f,
                                  "Minimum transverse momentum of V0"};
  Configurable<float> ConfV0DCADaughMax{"ConfV0DCADaughMax", 1.0f,
                                        "Maximum DCA between the V0 daughters"};
  Configurable<float> ConfV0CPAMin{"ConfV0CPAMin", 0.985f, "Minimum CPA of V0"};
  Configurable<float> ConfV0TranRadV0Min{"ConfV0TranRadV0Min", 0.5f,
                                         "Minimum transverse radius"};
  Configurable<float> ConfV0TranRadV0Max{"ConfV0TranRadV0Max", 200.f,
                                         "Maximum transverse radius"};
  Configurable<double> cMaxV0LifeTime{"cMaxV0LifeTime", 4,
                                      "Maximum V0 life time"};
  Configurable<double> cMaxV0DCA{"cMaxV0DCA", 0.3, "DCA V0 to PV"};
  Configurable<double> cSigmaMassLambda{"cSigmaMassLambda", 4,
                                        "n Sigma cut on Lambda mass"};
  Configurable<double> cWidthLambda{"cWidthLambda", 0.005, "Width of Lambda"};
  Configurable<float> ConfDaughEta{"ConfDaughEta", 0.8f,
                                   "V0 Daugh sel: max eta"};
  Configurable<float> ConfDaughTPCnclsMin{"ConfDaughTPCnclsMin", 70.f,
                                          "V0 Daugh sel: Min. nCls TPC"};
  Configurable<float> ConfDaughDCAMin{
    "ConfDaughDCAMin", 0.06f, "V0 Daugh sel:  Max. DCA Daugh to PV (cm)"};
  Configurable<float> ConfDaughPIDCuts{"ConfDaughPIDCuts", 4,
                                       "PID selections for KS0 daughters"};

  // Configurables for track selections
  Configurable<float> cfgCutPT{"cfgCutPT", 0.2f, "PT cut on daughter track"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8f, "Eta cut on daughter track"};
  Configurable<float> cfgCutDCAxy{"cfgCutDCAxy", 2.0f,
                                  "DCAxy range for tracks"};
  Configurable<float> cfgCutDCAz{"cfgCutDCAz", 2.0f, "DCAz range for tracks"};
  Configurable<float> nsigmaCutTPC{"nsigmacutTPC", 3.0,
                                   "Value of the TPC Nsigma cut"};
  Configurable<float> nsigmaCutTOF{"nsigmacutTOF", 3.0,
                                   "Value of the TOF Nsigma cut"};
  Configurable<float> nsigmaCutCombined{"nsigmaCutCombined", 3.0,
                                        "Value of the Combined Nsigma cut"};
  Configurable<bool> iscustomDCAcut{"iscustomDCAcut", false, "iscustomDCAcut"};
  Configurable<bool> ismanualDCAcut{"ismanualDCAcut", true, "ismanualDCAcut"};
  Configurable<int> cfgITScluster{"cfgITScluster", 0, "Number of ITS cluster"};
  Configurable<int> cfgTPCcluster{"cfgTPCcluster", 70, "Number of TPC cluster"};
  Configurable<bool> isNoTOF{"isNoTOF", false, "isNoTOF"};
  Configurable<bool> timFrameEvsel{"timFrameEvsel", false, "TPC Time frame boundary cut"};
  Configurable<bool> additionalEvsel{"additionalEvsel", false, "Additional event selcection"};

  // Event selection cuts - Alex
  TF1* fMultPVCutLow = nullptr;
  TF1* fMultPVCutHigh = nullptr;
  TF1* fMultCutLow = nullptr;
  TF1* fMultCutHigh = nullptr;
  TF1* fMultMultPVCut = nullptr;

  void init(InitContext const&)
  {
    // Axes
    AxisSpec vertexZAxis = {nBins, -10., 10., "vrtx_{Z} [cm]"};
    AxisSpec ptAxis = {200, 0.0f, 20.0f, "#it{p}_{T} (GeV/#it{c})"};

    // Event selection
    rEventSelection.add("hVertexZRec", "hVertexZRec",
                        {HistType::kTH1F, {vertexZAxis}});
    rEventSelection.add("hmult", "Centrality distribution", kTH1F,
                        {{100, 0.0f, 100.0f}});

    // from sourav da
    AxisSpec phiAxis = {500, -6.28, 6.28, "phi"};
    AxisSpec resAxis = {400, -2, 2, "Res"};
    AxisSpec centAxis = {100, 0, 100, "V0M (%)"};

    if (QA) {
      histos.add("hFTOCvsTPCNoCut", "Mult correlation FT0C vs. TPC without any cut", kTH2F, {{80, 0.0f, 80.0f}, {100, -0.5f, 5999.5f}});
      histos.add("hFTOCvsTPCSelected", "Mult correlation FT0C vs. TPC after selection", kTH2F, {{80, 0.0f, 80.0f}, {100, -0.5f, 5999.5f}});
      histos.add("hPsiFT0C", "PsiFT0C", kTH2F, {centAxis, phiAxis});
      histos.add("hPsiFT0A", "PsiFT0A", kTH2F, {centAxis, phiAxis});
      histos.add("hPsiTPC", "PsiTPC", kTH2F, {centAxis, phiAxis});
      histos.add("hPsiTPCR", "PsiTPCR", kTH2F, {centAxis, phiAxis});
      histos.add("hPsiTPCL", "PsiTPCL", kTH2F, {centAxis, phiAxis});
      // histogram for resolution
      histos.add("ResFT0CTPC", "ResFT0CTPC", kTH2F, {centAxis, resAxis});
      histos.add("ResFT0CTPCR", "ResFT0CTPCR", kTH2F, {centAxis, resAxis});
      histos.add("ResFT0CTPCL", "ResFT0CTPCL", kTH2F, {centAxis, resAxis});
      histos.add("ResTPCRTPCL", "ResTPCRTPCL", kTH2F, {centAxis, resAxis});
      histos.add("ResFT0CFT0A", "ResFT0CFT0A", kTH2F, {centAxis, resAxis});
      histos.add("ResFT0ATPC", "ResFT0ATPC", kTH2F, {centAxis, resAxis});
    }

    if (QAv0) {
      // Lambda topological cuts
      histos.add("hDCAV0Daughters", "hDCAV0Daughters",
                 {HistType::kTH1F, {{50, 0.0f, 5.0f}}});
      histos.add("hLT", "hLT", {HistType::kTH1F, {{100, 0.0f, 50.0f}}});
      histos.add("hV0CosPA", "hV0CosPA",
                 {HistType::kTH1F, {{100, 0.95f, 1.f}}});
    }

    // CKStar histograms
    histos.add("hLamInvMass",
               "Invariant mass of Lambda baryon", kTHnSparseF,
               {{100, 0.0, 100.0}, {200, 0.0f, 20.0f}, {200, 1.0, 1.2}, {100, -1, 1}}, true);

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

  double massLambda = TDatabasePDG::Instance()
                        ->GetParticle(kLambda0)
                        ->Mass(); // FIXME: Get from the common header

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
    if (iscustomDCAcut &&
        (!candidate.isGlobalTrack() || !candidate.isPVContributor() ||
         candidate.itsNCls() < cfgITScluster)) {
      return false;
    }

    if (ismanualDCAcut && !(candidate.isGlobalTrackWoDCA() && candidate.isPVContributor() && std::abs(candidate.dcaXY()) < cfgCutDCAxy && std::abs(candidate.dcaZ()) < cfgCutDCAz && candidate.itsNCls() > cfgITScluster && candidate.tpcNClsFound() > cfgTPCcluster)) {
      return false;
    }

    return true;
  }

  template <typename T>
  bool selectionPID(const T& candidate)
  {

    if (!candidate.hasTOF() &&
        std::abs(candidate.tpcNSigmaPi()) < nsigmaCutTPC) {
      return true;
    } else if (candidate.hasTOF() &&
               std::abs(candidate.tofNSigmaPi()) < nsigmaCutTOF) {
      return true;
    }

    return false;
  }

  template <typename Collision, typename V0>
  bool SelectionV0(Collision const& collision, V0 const& candidate,
                   float multiplicity)
  {
    if (fabs(candidate.dcav0topv()) > cMaxV0DCA) {
      return false;
    }

    const float pT = candidate.pt();
    const float tranRad = candidate.v0radius();
    const float dcaDaughv0 = candidate.dcaV0daughters();
    const float cpav0 = candidate.v0cosPA();
    float CtauK0s = candidate.distovertotmom(collision.posX(), collision.posY(),
                                             collision.posZ()) *
                    TDatabasePDG::Instance()
                      ->GetParticle(kLambda0)
                      ->Mass(); // FIXME: Get from the common header
    // float lowmasscutks0 = 0.497 - cWidthKs0 * cSigmaMassKs0;
    // float highmasscutks0 = 0.497 + cWidthKs0 * cSigmaMassKs0;
    float lowmasscutlambda = 1.0;
    float highmasscutlambda = 1.2;
    // float decayLength = candidate.distovertotmom(collision.posX(),
    // collision.posY(), collision.posZ()) *
    // RecoDecay::sqrtSumOfSquares(candidate.px(), candidate.py(),
    // candidate.pz());

    if (pT < ConfV0PtMin) {
      return false;
    }
    if (dcaDaughv0 > ConfV0DCADaughMax) {
      return false;
    }
    if (cpav0 < ConfV0CPAMin) {
      return false;
    }
    if (tranRad < ConfV0TranRadV0Min) {
      return false;
    }
    if (tranRad > ConfV0TranRadV0Max) {
      return false;
    }
    if (fabs(CtauK0s) > cMaxV0LifeTime ||
        candidate.mLambda() < lowmasscutlambda ||
        candidate.mLambda() > highmasscutlambda) {
      return false;
    }

    if (QAv0) {
      histos.fill(HIST("hLT"), CtauK0s);
      histos.fill(HIST("hDCAV0Daughters"), candidate.dcaV0daughters());
      histos.fill(HIST("hV0CosPA"), candidate.v0cosPA());
    }
    return true;
  }

  template <typename T>
  bool isSelectedV0Daughter(T const& track, float charge,
                            double nsigmaV0Daughter)
  {
    const auto eta = track.eta();
    const auto tpcNClsF = track.tpcNClsFound();
    const auto dcaXY = track.dcaXY();
    const auto sign = track.sign();

    if (!track.hasTPC())
      return false;
    if (track.tpcNClsCrossedRows() < 70)
      return false;
    if (track.tpcCrossedRowsOverFindableCls() < 0.8)
      return false;

    if (charge < 0 && sign > 0) {
      return false;
    }
    if (charge > 0 && sign < 0) {
      return false;
    }
    if (std::abs(eta) > ConfDaughEta) {
      return false;
    }
    if (tpcNClsF < ConfDaughTPCnclsMin) {
      return false;
    }
    if (std::abs(dcaXY) < ConfDaughDCAMin) {
      return false;
    }
    if (std::abs(nsigmaV0Daughter) > ConfDaughPIDCuts) {
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

  // Defining filters for events (event selection)
  // Processed events will be already fulfilling the event selection
  // requirements
  // Filter eventFilter = (o2::aod::evsel::sel8 == true);
  Filter posZFilter = (nabs(o2::aod::collision::posZ) < cutzvertex);

  Filter acceptanceFilter =
    (nabs(aod::track::eta) < cfgCutEta && nabs(aod::track::pt) > cfgCutPT);

  using EventCandidates = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::FV0Mults, aod::TPCMults, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0Cs, aod::CentFT0As, aod::Mults, aod::EPCalibrationTables>>;

  using TrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA,
                                                  aod::TrackSelection, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullPr, aod::pidTOFFullPr>>;
  using V0TrackCandidate = aod::V0Datas;

  ROOT::Math::PxPyPzMVector Lambda;

  void processSE(EventCandidates::iterator const& collision,
                 TrackCandidates const& tracks, aod::V0Datas const& V0s,
                 aod::BCs const&)

  {

    if (!collision.sel8()) {
      return;
    }

    float centrality = 0.0f;
    centrality = collision.centFT0C();
    auto multTPC = collision.multNTracksPV();

    if (timFrameEvsel && (!collision.selection_bit(aod::evsel::kNoTimeFrameBorder) || !collision.selection_bit(aod::evsel::kNoITSROFrameBorder))) {
      return;
    }

    if (additionalEvsel && !eventSelected(collision, centrality)) {
      return;
    }

    if (QA) {
      histos.fill(HIST("hFTOCvsTPCNoCut"), centrality, multTPC);
    }
    if (!collision.triggereventep()) {
      return;
    }
    auto psiFT0C = collision.psiFT0C();
    auto psiFT0A = collision.psiFT0A();
    auto psiTPC = collision.psiTPC();
    auto psiTPCR = collision.psiTPCR();
    auto psiTPCL = collision.psiTPCL();
    if (QA) {
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
    }

    // Fill the event counter
    rEventSelection.fill(HIST("hVertexZRec"), collision.posZ());
    rEventSelection.fill(HIST("hmult"), centrality);

    for (auto& v0 : V0s) {

      auto postrack = v0.template posTrack_as<TrackCandidates>();
      auto negtrack = v0.template negTrack_as<TrackCandidates>();
      double nTPCSigmaPos[1]{postrack.tpcNSigmaPr()};
      double nTPCSigmaNeg[1]{negtrack.tpcNSigmaPi()};

      if (!isSelectedV0Daughter(postrack, 1, nTPCSigmaPos[0])) {
        continue;
      }
      if (!isSelectedV0Daughter(negtrack, -1, nTPCSigmaNeg[0])) {
        continue;
      }

      if (!SelectionV0(collision, v0, centrality)) {
        continue;
      }

      Lambda = ROOT::Math::PxPyPzMVector(v0.px(), v0.py(), v0.pz(), massLambda);

      auto phiminuspsi = GetPhiInRange(Lambda.Phi() - psiFT0C);
      auto v2 = TMath::Cos(2.0 * phiminuspsi);

      if (TMath::Abs(Lambda.Rapidity()) < 0.5) {
        histos.fill(HIST("hLamInvMass"), centrality,
                    Lambda.Pt(), v0.mLambda(), v2);
      }
    }
  }

  PROCESS_SWITCH(lambdav2, processSE, "Process Same event", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<lambdav2>(cfgc)};
}
