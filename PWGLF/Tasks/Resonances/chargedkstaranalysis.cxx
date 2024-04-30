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
/// \brief this is a code for the CKS resonance
/// \author prottay das
/// \since 13/01/2024

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

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;

struct chargedkstaranalysis {

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

  HistogramRegistry rGenParticles{"genParticles", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rRecParticles{"recParticles", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  // Configurable for histograms
  Configurable<int> nBins{"nBins", 100, "N bins in all histos"};

  // Confugrable for QA histograms
  Configurable<bool> QAbefore{"QAbefore", false, "QAbefore"};
  Configurable<bool> QAafter{"QAafter", false, "QAafter"};
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
  Configurable<double> cMaxV0LifeTime{"cMaxV0LifeTime", 15,
                                      "Maximum V0 life time"};
  Configurable<double> cMaxV0DCA{"cMaxV0DCA", 0.3, "DCA V0 to PV"};
  Configurable<double> cSigmaMassKs0{"cSigmaMassKs0", 4,
                                     "n Sigma cut on KS0 mass"};
  Configurable<double> cWidthKs0{"cWidthKs0", 0.005, "Width of KS0"};

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
  Configurable<float> nsigmaCutCombined{"nsigmaCutCombined", 3.0,
                                        "Value of the Combined Nsigma cut"};
  Configurable<int> cfgNoMixedEvents{"cfgNoMixedEvents", 5,
                                     "Number of mixed events per event"};
  Configurable<bool> cfgMultFT0{"cfgMultFT0", false, "cfgMultFT0"};
  Configurable<bool> cfgCentFT0C{"cfgCentFT0C", true, "cfgCentFT0C"};
  Configurable<bool> iscustomDCAcut{"iscustomDCAcut", false, "iscustomDCAcut"};
  Configurable<bool> ismanualDCAcut{"ismanualDCAcut", true, "ismanualDCAcut"};
  Configurable<int> cfgITScluster{"cfgITScluster", 0, "Number of ITS cluster"};
  ConfigurableAxis cMixMultBins{"cMixMultBins", {VARIABLE_WIDTH, 0.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f, 100.0f}, "Mixing bins - multiplicity"};
  Configurable<bool> isMC{"isMC", true, "Run MC"};
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
    AxisSpec K0ShortMassAxis = {200, 0.45f, 0.55f,
                                "#it{M}_{inv} [GeV/#it{c}^{2}]"};
    AxisSpec vertexZAxis = {nBins, -10., 10., "vrtx_{Z} [cm]"};
    AxisSpec ptAxis = {200, 0.0f, 20.0f, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec multAxis = {100, 0.0f, 100.0f, "Multiplicity"};

    // Histograms
    // Event selection
    rEventSelection.add("hVertexZRec", "hVertexZRec",
                        {HistType::kTH1F, {vertexZAxis}});
    rEventSelection.add("hmult", "Centrality distribution", kTH1F,
                        {{200, 0.0f, 200.0f}});

    // for primary tracks
    if (QAbefore && QAafter) {
      histos.add("hNsigmaPionTPC_before", "NsigmaPion TPC distribution before",
                 kTH1F, {{200, -10.0f, 10.0f}});
      histos.add("hNsigmaPionTOF_before", "NsigmaPion TOF distribution before",
                 kTH1F, {{200, -10.0f, 10.0f}});

      histos.add("hEta_after", "Eta distribution", kTH1F, {{200, -1.0f, 1.0f}});
      histos.add("hDcaxy_after", "Dcaxy distribution", kTH1F,
                 {{200, -10.0f, 10.0f}});
      histos.add("hDcaz_after", "Dcaz distribution", kTH1F,
                 {{200, -10.0f, 10.0f}});
      histos.add("hNsigmaPionTPC_after", "NsigmaPion TPC distribution", kTH1F,
                 {{200, -10.0f, 10.0f}});
      histos.add("hNsigmaPionTOF_after", "NsigmaPion TOF distribution", kTH1F,
                 {{200, -10.0f, 10.0f}});
    }

    if (QAv0) {
      // K0s reconstruction
      histos.add(
        "hMassvsptvsmult", "hMassvsptvsmult",
        {HistType::kTHnSparseF, {{K0ShortMassAxis}, {ptAxis}, {multAxis}}},
        true);
      // K0s topological/PID cuts
      histos.add("hDCAV0Daughters", "hDCAV0Daughters",
                 {HistType::kTH1F, {{50, 0.0f, 5.0f}}});
      histos.add("hLT", "hLT", {HistType::kTH1F, {{100, 0.0f, 50.0f}}});
      histos.add("hV0CosPA", "hV0CosPA",
                 {HistType::kTH1F, {{100, 0.95f, 1.f}}});
    }

    // CKStar histograms
    histos.add("h3CKSInvMassUnlikeSign",
               "Invariant mass of CKS meson Unlike Sign", kTHnSparseF,
               {{200, 0.0, 200.0}, {200, 0.0f, 20.0f}, {90, 0.6, 1.5}}, true);
    histos.add("h3CKSInvMassMixed", "Invariant mass of CKS meson Mixed",
               kTHnSparseF,
               {{200, 0.0, 200.0}, {200, 0.0f, 20.0f}, {90, 0.6, 1.5}}, true);

    if (isMC) {
      rGenParticles.add("hMC", "Gen MC Event statistics", kTH1F, {{10, 0.0f, 10.0f}});
      rRecParticles.add("hMCRec", "Rec MC Event statistics", kTH1F, {{10, 0.0f, 10.0f}});
      rGenParticles.add("hPtK0ShortGen", "hPtK0ShortGen", {HistType::kTH1F, {{ptAxis}}});
      rGenParticles.add("hCKSGen", "hCKSGen", {HistType::kTH1F, {{ptAxis}}});
      rRecParticles.add("hCKSRec", "hCKSRec", {HistType::kTH1F, {{ptAxis}}});
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

  double massPi = TDatabasePDG::Instance()
                    ->GetParticle(kPiPlus)
                    ->Mass(); // FIXME: Get from the common header
  double massK0s = TDatabasePDG::Instance()
                     ->GetParticle(kK0Short)
                     ->Mass(); // FIXME: Get from the common header
  double massKa = o2::constants::physics::MassKPlus;
  ROOT::Math::PtEtaPhiMVector CKSVector;

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
    if (ismanualDCAcut &&
        (!candidate.isGlobalTrackWoDCA() || !candidate.isPVContributor() ||
         std::abs(candidate.dcaXY()) > cfgCutDCAxy ||
         std::abs(candidate.dcaZ()) > cfgCutDCAz ||
         candidate.itsNCls() < cfgITScluster)) {
      return false;
    }
    return true;
  }

  template <typename T>
  bool selectionPID(const T& candidate)
  {

    if (candidate.hasTOF() &&
        (candidate.tofNSigmaPi() * candidate.tofNSigmaPi() +
         candidate.tpcNSigmaPi() * candidate.tpcNSigmaPi()) <
          (nsigmaCutCombined * nsigmaCutCombined)) {
      return true;
    }
    if (!candidate.hasTOF() &&
        std::abs(candidate.tpcNSigmaPi()) < nsigmaCutTPC) {
      return true;
    }

    /*
    if (candidate.hasTOF() &&
        (candidate.tofNSigmaKa() * candidate.tofNSigmaKa() +
         candidate.tpcNSigmaKa() * candidate.tpcNSigmaKa()) <
          (nsigmaCutCombined * nsigmaCutCombined)) {
      return true;
    }
    if (!candidate.hasTOF() &&
        std::abs(candidate.tpcNSigmaKa()) < nsigmaCutTPC) {
      return true;
    }
    */
    return false;
  }

  template <typename Collision, typename V0>
  bool SelectionV0(Collision const& collision, V0 const& candidate,
                   float multiplicity)
  {
    if (fabs(candidate.dcav0topv()) > cMaxV0DCA) {
      return false;
    }

    if (TMath::Abs(candidate.yK0Short()) > 0.5) {
      return false;
    }

    const float qtarm = candidate.qtarm();
    const float alph = candidate.alpha();
    float arm = qtarm / alph;
    const float pT = candidate.pt();
    const float tranRad = candidate.v0radius();
    const float dcaDaughv0 = candidate.dcaV0daughters();
    const float cpav0 = candidate.v0cosPA();
    float CtauK0s = candidate.distovertotmom(collision.posX(), collision.posY(),
                                             collision.posZ()) *
                    TDatabasePDG::Instance()
                      ->GetParticle(kK0Short)
                      ->Mass(); // FIXME: Get from the common header
    float lowmasscutks0 = 0.497 - cWidthKs0 * cSigmaMassKs0;
    float highmasscutks0 = 0.497 + cWidthKs0 * cSigmaMassKs0;
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
        candidate.mK0Short() < lowmasscutks0 ||
        candidate.mK0Short() > highmasscutks0) {
      return false;
    }
    if (arm < 0.2) {
      return false;
    }

    if (QAv0) {
      histos.fill(HIST("hLT"), CtauK0s);
      histos.fill(HIST("hMassvsptvsmult"), candidate.mK0Short(), candidate.pt(),
                  multiplicity);
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

  // Defining filters for events (event selection)
  // Processed events will be already fulfilling the event selection
  // requirements
  Filter eventFilter = (o2::aod::evsel::sel8 == true);
  Filter posZFilter = (nabs(o2::aod::collision::posZ) < cutzvertex);
  Filter posZFilterMC = (nabs(o2::aod::mccollision::posZ) < cutzvertex);

  Filter acceptanceFilter =
    (nabs(aod::track::eta) < cfgCutEta && nabs(aod::track::pt) > cfgCutPT);
  Filter DCAcutFilter = (nabs(aod::track::dcaXY) < cfgCutDCAxy) &&
                        (nabs(aod::track::dcaZ) < cfgCutDCAz);

  using EventCandidatesMC = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>;
  /*using EventCandidates = soa::Filtered<
    soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::MultZeqs,
    aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs>>;*/
  using EventCandidates = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::FV0Mults, aod::TPCMults, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0Cs, aod::CentFT0As, aod::Mults>>;

  using TrackCandidates = soa::Filtered<
    soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA,
              aod::TrackSelection, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullKa, aod::pidTOFFullKa>>;
  using TrackCandidatesMC = soa::Filtered<
    soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA,
              aod::TrackSelection, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullKa, aod::pidTOFFullKa, aod::McTrackLabels>>;

  using V0TrackCandidatesMC = soa::Join<aod::V0Datas, aod::McV0Labels>;
  using V0TrackCandidate = aod::V0Datas;

  ConfigurableAxis axisVertex{
    "axisVertex",
    {20, -10, 10},
    "vertex axis for bin"};
  ConfigurableAxis axisMultiplicityClass{
    "axisMultiplicityClass",
    {1, 0, 100},
    "multiplicity percentile for bin"};
  ConfigurableAxis axisMultiplicity{
    "axisMultiplicity",
    {2, 0, 100},
    "TPC multiplicity  for bin"};

  using BinningTypeTPCMultiplicity =
    ColumnBinningPolicy<aod::collision::PosZ, aod::mult::MultTPC>;
  // using BinningTypeVertexContributor =
  // ColumnBinningPolicy<aod::collision::PosZ, aod::collision::NumContrib>;
  using BinningTypeCentralityM =
    ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0M>;
  using BinningTypeVertexContributor =
    ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0C>;

  BinningTypeVertexContributor binningOnPositions{
    {axisVertex, axisMultiplicity},
    true};

  Pair<EventCandidates, TrackCandidates, V0TrackCandidate,
       BinningTypeVertexContributor>
    pair{binningOnPositions, cfgNoMixedEvents, -1, &cache};

  /*
  SameKindPair<EventCandidates, TrackCandidates,
       BinningTypeVertexContributor>
    pair{binningOnPositions, cfgNoMixedEvents, -1, &cache};
  */

  void processSE(EventCandidates::iterator const& collision,
                 TrackCandidates const& tracks, aod::V0Datas const& V0s,
                 aod::BCs const&)

  {

    if (!collision.sel8()) {
      return;
    }

    std::vector<ROOT::Math::PtEtaPhiMVector> pions, kshorts;
    std::vector<ROOT::Math::PtEtaPhiMVector> pions2;
    std::vector<int64_t> PionIndex = {};
    std::vector<int64_t> PionSign = {};
    std::vector<int64_t> PioncollIndex = {};
    std::vector<int64_t> PionIndex2 = {};
    std::vector<int64_t> PionSign2 = {};
    std::vector<int64_t> PioncollIndex2 = {};
    std::vector<int64_t> V0collIndex = {};
    std::vector<int64_t> KshortPosDaughIndex = {};
    std::vector<int64_t> KshortNegDaughIndex = {};
    /*
    float multiplicity = 0.0f;
    if (cfgMultFT0)
      multiplicity = collision.multZeqFT0A() + collision.multZeqFT0C();
    if (cfgMultFT0 == 0 && cfgCentFT0C == 1)
      multiplicity = collision.centFT0C();
    if (cfgMultFT0 == 0 && cfgCentFT0C == 0)
      multiplicity = collision.centFT0M();
    */
    float centrality = 0.0f;
    centrality = collision.centFT0C();

    if (timFrameEvsel && (!collision.selection_bit(aod::evsel::kNoTimeFrameBorder) || !collision.selection_bit(aod::evsel::kNoITSROFrameBorder))) {
      return;
    }

    if (additionalEvsel && !eventSelected(collision, centrality)) {
      return;
    }

    // Fill the event counter
    rEventSelection.fill(HIST("hVertexZRec"), collision.posZ());
    rEventSelection.fill(HIST("hmult"), centrality);

    for (auto track1 : tracks) {

      if (QAbefore) {
        histos.fill(HIST("hNsigmaPionTPC_before"), track1.tpcNSigmaPi());
        histos.fill(HIST("hNsigmaPionTOF_before"), track1.tofNSigmaPi());
      }

      if (!selectionPID(track1))
        continue; // for primary particle PID

      if (!selectionTrack(track1)) {
        continue;
      }

      if (QAafter) {
        histos.fill(HIST("hEta_after"), track1.eta());
        histos.fill(HIST("hDcaxy_after"), track1.dcaXY());
        histos.fill(HIST("hDcaz_after"), track1.dcaZ());
        histos.fill(HIST("hNsigmaPionTPC_after"), track1.tpcNSigmaPi());
        histos.fill(HIST("hNsigmaPionTOF_after"), track1.tofNSigmaPi());
      }

      ROOT::Math::PtEtaPhiMVector temp1(track1.pt(), track1.eta(), track1.phi(),
                                        massPi);
      pions.push_back(temp1);
      PionIndex.push_back(track1.globalIndex());
      PioncollIndex.push_back(track1.collisionId());
      PionSign.push_back(track1.sign());
    } // track loop ends

    for (auto& v0 : V0s) {

      auto postrack = v0.template posTrack_as<TrackCandidates>();
      auto negtrack = v0.template negTrack_as<TrackCandidates>();
      double nTPCSigmaPos[1]{postrack.tpcNSigmaPi()};
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

      ROOT::Math::PtEtaPhiMVector temp2(v0.pt(), v0.eta(), v0.phi(), massK0s);
      kshorts.push_back(temp2);
      V0collIndex.push_back(v0.collisionId());
      KshortPosDaughIndex.push_back(postrack.globalIndex());
      KshortNegDaughIndex.push_back(negtrack.globalIndex());
    }

    if (pions.size() != 0 && kshorts.size() != 0) {
      // if (pions.size() != 0 && pions2.size() != 0) {
      for (auto ipion = pions.begin(); ipion != pions.end(); ++ipion) {
        auto i1 = std::distance(pions.begin(), ipion);
        if (PionSign.at(i1) == 0)
          continue;
        for (auto ikshort = kshorts.begin(); ikshort != kshorts.end();
             ++ikshort) {
          // for (auto ikshort = pions2.begin(); ikshort != pions2.end();
          //    ++ikshort) {
          auto i3 = std::distance(kshorts.begin(), ikshort);
          if (PionIndex.at(i1) == KshortPosDaughIndex.at(i3))
            continue;
          if (PionIndex.at(i1) == KshortNegDaughIndex.at(i3))
            continue;
          // if (PioncollIndex.at(i1) != V0collIndex.at(i3))
          // continue;

          CKSVector = pions.at(i1) + kshorts.at(i3);
          if (TMath::Abs(CKSVector.Rapidity()) < 0.5) {
            histos.fill(HIST("h3CKSInvMassUnlikeSign"), centrality,
                        CKSVector.Pt(), CKSVector.M());
          }
        }
      }
    }
  }

  PROCESS_SWITCH(chargedkstaranalysis, processSE, "Process Same event", false);

  void processME(EventCandidates const& /*collisions*/,
                 TrackCandidates const& /*tracks*/, V0TrackCandidate const& /*V0s*/)

  {

    for (auto& [c1, tracks1, c2, tracks2] : pair) {

      if (!c1.sel8()) {
        continue;
      }
      if (!c2.sel8()) {
        continue;
      }
      /*
      float multiplicity = 0.0f;
      if (cfgMultFT0)
        multiplicity = c1.multZeqFT0A() + c1.multZeqFT0C();
      if (cfgMultFT0 == 0 && cfgCentFT0C == 1)
        multiplicity = c1.centFT0C();
      if (cfgMultFT0 == 0 && cfgCentFT0C == 0)
        multiplicity = c1.centFT0M();
      */
      auto centrality = c1.centFT0C();
      auto centrality2 = c2.centFT0C();

      if (timFrameEvsel && (!c1.selection_bit(aod::evsel::kNoTimeFrameBorder) || !c2.selection_bit(aod::evsel::kNoTimeFrameBorder) || !c1.selection_bit(aod::evsel::kNoITSROFrameBorder) || !c2.selection_bit(aod::evsel::kNoITSROFrameBorder))) {
        continue;
      }

      if (additionalEvsel && !eventSelected(c1, centrality)) {
        continue;
      }
      if (additionalEvsel && !eventSelected(c2, centrality2)) {
        continue;
      }

      for (auto& [t1, t2] : o2::soa::combinations(
             o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {

        if (!selectionTrack(t1))
          continue;
        if (!selectionPID(t1))
          continue;
        if (t1.sign() == 0)
          continue;

        if (!SelectionV0(c2, t2, centrality2))
          continue;

        auto postrack = t2.template posTrack_as<TrackCandidates>();
        auto negtrack = t2.template negTrack_as<TrackCandidates>();
        double nTPCSigmaPos[1]{postrack.tpcNSigmaPi()};
        double nTPCSigmaNeg[1]{negtrack.tpcNSigmaPi()};

        if (!isSelectedV0Daughter(postrack, 1, nTPCSigmaPos[0])) {
          continue;
        }
        if (!isSelectedV0Daughter(negtrack, -1, nTPCSigmaNeg[0])) {
          continue;
        }

        TLorentzVector pi;
        pi.SetPtEtaPhiM(t1.pt(), t1.eta(), t1.phi(), massPi);
        TLorentzVector KSh;
        KSh.SetPtEtaPhiM(t2.pt(), t2.eta(), t2.phi(), massK0s);

        TLorentzVector CKSmix = pi + KSh;

        if (TMath::Abs(CKSmix.Rapidity()) < 0.5) {
          histos.fill(HIST("h3CKSInvMassMixed"), centrality, CKSmix.Pt(),
                      CKSmix.M());
        }
      }
    }
  }

  PROCESS_SWITCH(chargedkstaranalysis, processME, "Process Mixed event", false);

  // taken from Sourav da
  void processGenMC(aod::McCollision const& mcCollision, aod::McParticles& mcParticles, const soa::SmallGroups<EventCandidatesMC>& collisions)
  {

    if (std::abs(mcCollision.posZ()) < cutzvertex)
      rGenParticles.fill(HIST("hMC"), 0.5);
    std::vector<int64_t> SelectedEvents(collisions.size());
    int nevts = 0;
    for (const auto& collision : collisions) {
      if (!collision.sel8() || std::abs(collision.mcCollision().posZ()) > cutzvertex) {
        continue;
      }
      rGenParticles.fill(HIST("hMC"), 1.5);
      if (timFrameEvsel && (!collision.selection_bit(aod::evsel::kNoTimeFrameBorder) || !collision.selection_bit(aod::evsel::kNoITSROFrameBorder))) {
        continue;
      }
      rGenParticles.fill(HIST("hMC"), 2.5);
      SelectedEvents[nevts++] = collision.mcCollision_as<aod::McCollisions>().globalIndex();
    }
    SelectedEvents.resize(nevts);
    const auto evtReconstructedAndSelected = std::find(SelectedEvents.begin(), SelectedEvents.end(), mcCollision.globalIndex()) != SelectedEvents.end();

    if (!evtReconstructedAndSelected) { // Check that the event is reconstructed and that the reconstructed events pass the selection
      return;
    }

    rGenParticles.fill(HIST("hMC"), 3.5);
    for (auto& mcParticle : mcParticles) {
      if (std::abs(mcParticle.y()) >= 0.5) {
        continue;
      }
      rGenParticles.fill(HIST("hMC"), 4.5);
      if (std::abs(mcParticle.pdgCode()) != 323) {
        continue;
      }
      rGenParticles.fill(HIST("hMC"), 5.5);
      auto kDaughters = mcParticle.daughters_as<aod::McParticles>();
      if (kDaughters.size() != 2) {
        continue;
      }

      rGenParticles.fill(HIST("hMC"), 6.5);
      auto daughts = false;
      auto daughtp = false;
      // int count = 0;
      for (auto kCurrentDaughter : kDaughters) {
        // LOG(info) << "Daughters PDG:\t" << count<<" "<<kCurrentDaughter.pdgCode();
        if (kCurrentDaughter.pdgCode() == 311) {
          auto kDaughter2 = kCurrentDaughter.daughters_as<aod::McParticles>();
          for (auto kCurrentDaughter2 : kDaughter2) {
            if (kCurrentDaughter2.pdgCode() == 310)
              daughts = true;
          }
        } else if (std::abs(kCurrentDaughter.pdgCode()) == 211) {
          if (kCurrentDaughter.isPhysicalPrimary() == 1)
            daughtp = true;
        }
        // count += 1;
      }
      rGenParticles.fill(HIST("hMC"), 7.5);
      if (daughtp && daughts) {
        rGenParticles.fill(HIST("hCKSGen"), mcParticle.pt());
      }
    }
  }

  void processRecMC(EventCandidatesMC::iterator const& collision,
                    TrackCandidatesMC const& tracks, V0TrackCandidatesMC const& V0s,
                    aod::McParticles const& /*mcParticles*/, aod::McCollisions const& /*mcCollisions*/)

  {

    if (!collision.has_mcCollision()) {
      return;
    }
    if (std::abs(collision.mcCollision().posZ()) > cutzvertex || !collision.sel8()) {
      return;
    }

    rRecParticles.fill(HIST("hMCRec"), 0.5);

    if (timFrameEvsel && (!collision.selection_bit(aod::evsel::kNoTimeFrameBorder) || !collision.selection_bit(aod::evsel::kNoITSROFrameBorder))) {
      return;
    }

    rRecParticles.fill(HIST("hMCRec"), 1.5);

    float centrality = 0.0f;

    for (auto track1 : tracks) {

      if (!selectionPID(track1))
        continue; // for primary particle PID

      if (!track1.has_mcParticle()) {
        continue;
      }

      if (!selectionTrack(track1)) {
        continue;
      }

      auto mctrack1 = track1.mcParticle();

      if (!mctrack1.isPhysicalPrimary()) {
        continue;
      }

      for (auto& v0 : V0s) {

        if (!v0.has_mcParticle()) {
          continue;
        }

        auto postrack = v0.template posTrack_as<TrackCandidatesMC>();
        auto negtrack = v0.template negTrack_as<TrackCandidatesMC>();

        if (!postrack.has_mcParticle())
          continue;                     // Checking that the daughter tracks come from particles and are not fake
        if (!negtrack.has_mcParticle()) // Checking that the daughter tracks come from particles and are not fake
          continue;

        // auto posParticle = postrack.mcParticle();
        // auto negParticle = negtrack.mcParticle();

        double nTPCSigmaPos[1]{postrack.tpcNSigmaPi()};
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

        auto mctrackv0 = v0.mcParticle();

        /*
            for (auto track2 : tracks) {

            if (!selectionPID(track2))
        continue; // for primary particle PID

            if (!track2.has_mcParticle()) {
              continue;
            }

            if (!selectionTrack(track2)) {
              continue;
            }

            auto mctrack2 = track2.mcParticle();


            if (!mctrack2.isPhysicalPrimary()) {
        continue;
            }
        */

        int track1PDG = std::abs(mctrack1.pdgCode());
        // int track2PDG = std::abs(mctrack2.pdgCode());
        int trackv0PDG = std::abs(mctrackv0.pdgCode());

        if (postrack.globalIndex() == track1.globalIndex())
          continue;
        if (negtrack.globalIndex() == track1.globalIndex())
          continue;

        rRecParticles.fill(HIST("hMCRec"), 2.5);

        if (track1PDG != 211) {
          continue;
        }
        // if (track2PDG != 321) {
        // continue;
        // }
        if (trackv0PDG != 310) {
          continue;
        }

        rRecParticles.fill(HIST("hMCRec"), 3.5);

        for (auto& mothertrack1 : mctrack1.mothers_as<aod::McParticles>()) {
          // for (auto& mothertrack2 : mctrack2.mothers_as<aod::McParticles>()) {
          for (auto& mothertrack2 : mctrackv0.mothers_as<aod::McParticles>()) {

            rRecParticles.fill(HIST("hMCRec"), 4.5);
            // LOG(info) << "Initial Mothers PDG:\t" <<mothertrack1.pdgCode()<<" "<<mothertrack2.pdgCode();

            if (mothertrack2.pdgCode() != 311) // K0
              continue;

            for (auto& mothertrack3 : mothertrack2.mothers_as<aod::McParticles>()) {

              // LOG(info) << "final Mothers PDG:\t" <<mothertrack3.pdgCode();

              if (std::abs(mothertrack3.pdgCode()) != 323) {
                continue;
              }

              if (mothertrack3.pdgCode() != mothertrack1.pdgCode()) {
                continue;
              }

              // LOG(info) << "final Mothers PDG:\t" <<mothertrack3.pdgCode()<<" "<<mothertrack3.globalIndex()<<" "<<mothertrack1.globalIndex();

              rRecParticles.fill(HIST("hMCRec"), 5.5);

              if (mothertrack3.globalIndex() != mothertrack1.globalIndex()) {
                continue;
              }

              rRecParticles.fill(HIST("hMCRec"), 6.5);
              /*
                    if (!mothertrack1.producedByGenerator()) {
                      continue;
                    }
                */
              // rRecParticles.fill(HIST("hMCRec"), 7.5);

              // LOG(info) << "Mothers PDG:\t" <<mothertrack1.pdgCode()<<" "<<mothertrack2.pdgCode();

              if (std::abs(mothertrack1.y()) >= 0.5) {
                continue;
              }

              rRecParticles.fill(HIST("hMCRec"), 7.5);

              rRecParticles.fill(HIST("hCKSRec"), mothertrack1.pt());
            }
          }
        }
      }
    } // track loop ends
  }

  PROCESS_SWITCH(chargedkstaranalysis, processGenMC, "Process Gen event", true);
  PROCESS_SWITCH(chargedkstaranalysis, processRecMC, "Process Rec event", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<chargedkstaranalysis>(cfgc)};
}
