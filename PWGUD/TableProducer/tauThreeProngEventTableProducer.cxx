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
//
/// \file tauThreeProngEventTableProducer.cxx
/// \brief Produces derived table from UD tables for tau pair production (3 prong)
///
/// \author Adam Matyja <adam.tomasz.matyja@cern.ch>, IFJ PAN, Cracow
/// \since  2025-09-06
//
// to run it execute:
// copts="--configuration json://tautauMC_modified_new.json -b"
// for MC
// oopts="--aod-writer-json saveDerivedConfig.json"
//  for data
// oopts="--aod-writer-json saveDerivedConfigData.json"
// o2-analysis-ud-tau-three-prong-event-table-producer $copts $oopts > output.log

//// C++ headers
#include "Math/Vector4D.h"

#include <algorithm>
#include <random>
#include <set>
#include <utility>
#include <vector>
//
//// O2 headers
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
// #include "Framework/HistogramRegistry.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"
//
//// O2Physics headers
#include "Common/CCDB/EventSelectionParams.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"
// #include "PWGUD/Core/UPCTauCentralBarrelHelperRL.h"
#include "PWGUD/Core/DGPIDSelector.h"
#include "PWGUD/Core/SGSelector.h"
#include "PWGUD/Core/UDHelpers.h"
#include "PWGUD/DataModel/TauThreeProngEventTables.h"
#include "PWGUD/DataModel/UDTables.h"

#include "Common/Core/RecoDecay.h"

#include "TPDGCode.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;

enum MyRecoProblem {
  NO_PROBLEM = 0,         // no problem
  MANY_RECO = 1,          // more than 1 reconstructed collision
  TOO_MANY_DAUGHTERS = 2, // more than 6 daughters from 2 taus
  TWO_TRACKS = 3          // more than 1 associated track to MC particle (tau daughter)
};

enum MyParticle {
  MyOtherParticle = -1,
  MyElectron = 1,
  MyMuon = 2,
  MyPion = 3,
  MyKaon = 4
};

struct TauThreeProngEventTableProducer {
  Produces<o2::aod::TrueTauFourTracks> trueTauFourTracks;
  Produces<o2::aod::DataTauFourTracks> dataTauFourTracks;

  // Global varialbes
  // Service<o2::framework::O2DatabasePDG> pdg;
  SGSelector sgSelector;

  // initialize histogram registry
  HistogramRegistry registrySkim{
    "registrySkim",
    {}};

  //// declare configurables
  // Configurable<bool> verboseInfo{"verboseInfo", false, {"Print general info to terminal; default it false."}};
  //
  // struct : ConfigurableGroup {
  //   Configurable<int> whichGapSide{"whichGapSide", 2, {"0 for side A, 1 for side C, 2 for both sides"}};
  //   Configurable<bool> useTrueGap{"useTrueGap", true, {"Calculate gapSide for a given FV0/FT0/ZDC thresholds"}};
  //   Configurable<int> cutNumContribs{"cutNumContribs", 2, {"How many contributors event has"}};
  //   Configurable<bool> useNumContribs{"useNumContribs", false, {"Use coll.numContribs as event cut"}};
  //   Configurable<int> cutRecoFlag{"cutRecoFlag", 1, {"0 = std mode, 1 = upc mode"}};
  //   Configurable<bool> useRecoFlag{"useRecoFlag", false, {"Use coll.flags as event cut"}};
  //   Configurable<int> cutRCTflag{"cutRCTflag", 0, {"0 = off, 1 = CBT, 2 = CBT+ZDC, 3 = CBThadron, 4 = CBThadron+ZDC"}};
  //   Configurable<float> cutTrueGapSideFV0{"cutTrueGapSideFV0", 180000, "FV0A threshold for SG selector"};
  //   Configurable<float> cutTrueGapSideFT0A{"cutTrueGapSideFT0A", 150., "FT0A threshold for SG selector"};
  //   Configurable<float> cutTrueGapSideFT0C{"cutTrueGapSideFT0C", 50., "FT0C threshold for SG selector"};
  //   Configurable<float> cutTrueGapSideZDC{"cutTrueGapSideZDC", 10000., "ZDC threshold for SG selector. 0 is <1n, 4.2 is <2n, 6.7 is <3n, 9.5 is <4n, 12.5 is <5n"};
  //   Configurable<float> cutFITtime{"cutFITtime", 40., "Maximum FIT time allowed. Default is 40ns"};
  //   Configurable<bool> cutEvTFb{"cutEvTFb", true, {"Event selection bit kNoTimeFrameBorder"}};
  //   Configurable<bool> cutEvITSROFb{"cutEvITSROFb", true, {"Event selection bit kNoITSROFrameBorder"}};
  //   Configurable<bool> cutEvSbp{"cutEvSbp", true, {"Event selection bit kNoSameBunchPileup"}};
  //   Configurable<bool> cutEvZvtxFT0vPV{"cutEvZvtxFT0vPV", false, {"Event selection bit kIsGoodZvtxFT0vsPV"}};
  //   Configurable<bool> cutEvVtxITSTPC{"cutEvVtxITSTPC", true, {"Event selection bit kIsVertexITSTPC"}};
  //   Configurable<float> cutEvOccupancy{"cutEvOccupancy", 100000., "Maximum allowed occupancy"};
  //   Configurable<bool> cutEvTrs{"cutEvTrs", false, {"Event selection bit kNoCollInTimeRangeStandard"}};
  //   Configurable<bool> cutEvTrofs{"cutEvTrofs", false, {"Event selection bit kNoCollInRofStandard"}};
  //   Configurable<bool> cutEvHmpr{"cutEvHmpr", false, {"Event selection bit kNoHighMultCollInPrevRof"}};
  // } cutSample;
  //
  // struct : ConfigurableGroup {
  //   Configurable<bool> applyGlobalTrackSelection{"applyGlobalTrackSelection", false, {"Applies cut on here defined global tracks"}};
  //   Configurable<float> cutMinPt{"cutMinPt", 0.1f, {"Global track cut"}};
  //   Configurable<float> cutMaxPt{"cutMaxPt", 1e10f, {"Global track cut"}};
  //   Configurable<float> cutMinEta{"cutMinEta", -0.8f, {"Global track cut"}};
  //   Configurable<float> cutMaxEta{"cutMaxEta", 0.8f, {"Global track cut"}};
  //   Configurable<float> cutMaxDCAz{"cutMaxDCAz", 2.f, {"Global track cut"}};
  //   Configurable<float> cutMaxDCAxy{"cutMaxDCAxy", 1e10f, {"Global track cut"}};
  //   Configurable<bool> applyPtDependentDCAxy{"applyPtDependentDCAxy", false, {"Global track cut"}};
  //   Configurable<bool> cutHasITS{"cutHasITS", true, {"Global track cut"}};
  //   Configurable<int> cutMinITSnCls{"cutMinITSnCls", 1, {"Global track cut"}};
  //   Configurable<float> cutMaxITSchi2{"cutMaxITSchi2", 36.f, {"Global track cut"}};
  //   Configurable<int> cutITShitsRule{"cutITShitsRule", 0, {"Global track cut"}};
  //   Configurable<bool> cutHasTPC{"cutHasTPC", true, {"Global track cut"}};
  //   Configurable<int> cutMinTPCnCls{"cutMinTPCnCls", 1, {"Global track cut"}};
  //   Configurable<int> cutMinTPCnClsXrows{"cutMinTPCnClsXrows", 70, {"Global track cut"}};
  //   Configurable<float> cutMinTPCnClsXrowsOverNcls{"cutMinTPCnClsXrowsOverNcls", 0.8f, {"Global track cut"}};
  //   Configurable<float> cutMaxTPCchi2{"cutMaxTPCchi2", 4.f, {"Global track cut"}};
  //   Configurable<bool> cutGoodITSTPCmatching{"cutGoodITSTPCmatching", true, {"Global track cut"}};
  //   Configurable<float> cutMaxTOFchi2{"cutMaxTOFchi2", 3.f, {"Global track cut"}};
  // } cutGlobalTrack;
  //
  // struct : ConfigurableGroup {
  //   Configurable<bool> preselUseTrackPID{"preselUseTrackPID", true, {"Apply weak PID check on tracks."}};
  //   Configurable<int> preselNgoodPVtracs{"preselNgoodPVtracs", 2, {"How many good PV tracks to select."}};
  //   Configurable<float> preselMinElectronNsigmaEl{"preselMinElectronNsigmaEl", 4.0, {"Good el candidate hypo in. Upper n sigma cut on el hypo of selected electron. What is more goes away."}};
  //   Configurable<float> preselMaxElectronNsigmaEl{"preselMaxElectronNsigmaEl", -2.0, {"Good el candidate hypo in. Lower n sigma cut on el hypo of selected electron. What is less goes away."}};
  //   Configurable<bool> preselElectronHasTOF{"preselElectronHasTOF", true, {"Electron candidated is required to hit TOF."}};
  //   Configurable<float> preselMinPionNsigmaEl{"preselMinPionNsigmaEl", 5.0, {"Good pi candidate hypo in. Upper n sigma cut on pi hypo of selected electron. What is more goes away."}};
  //   Configurable<float> preselMaxPionNsigmaEl{"preselMaxPionNsigmaEl", -5.0, {"Good pi candidate hypo in. Lower n sigma cut on pi hypo of selected electron. What is less goes away."}};
  //   Configurable<float> preselMinMuonNsigmaEl{"preselMinMuonNsigmaEl", 5.0, {"Good pi candidate hypo in. Upper n sigma cut on pi hypo of selected electron. What is more goes away."}};
  //   Configurable<float> preselMaxMuonNsigmaEl{"preselMaxMuonNsigmaEl", -5.0, {"Good pi candidate hypo in. Lower n sigma cut on pi hypo of selected electron. What is less goes away."}};
  //   Configurable<bool> preselMupionHasTOF{"preselMupionHasTOF", true, {"Mupion candidate is required to hit TOF."}};
  // } cutPreselect;

  // configurables
  Configurable<float> cutFV0{"cutFV0", 10000., "FV0A threshold"};
  Configurable<float> cutFT0A{"cutFT0A", 150., "FT0A threshold"};
  Configurable<float> cutFT0C{"cutFT0C", 50., "FT0C threshold"};
  Configurable<float> cutZDC{"cutZDC", 10000., "ZDC threshold"};
  Configurable<float> mGapSide{"mGapSide", 2, "gap selection"};

  //  ConfigurableAxis ptAxis{"ptAxis", {120, 0., 4.}, "#it{p} (GeV/#it{c})"};
  //  ConfigurableAxis dedxAxis{"dedxAxis", {100, 20., 160.}, "dE/dx"};
  //  ConfigurableAxis minvAxis{"minvAxis", {100, 0.5, 5.0}, "M_{inv} (GeV/#it{c}^{2})"};
  //  ConfigurableAxis phiAxis{"phiAxis", {120, 0., 3.2}, "#phi"};
  Configurable<bool> verbose{"verbose", {}, "Additional print outs"};

  // cut selection configurables
  //  Configurable<float> zvertexcut{"zvertexcut", 10., "Z vertex cut"};
  Configurable<float> trkEtacut{"trkEtacut", 0.9, "max track eta cut"};
  //   Configurable<bool> sameSign{"sameSign", {}, "Switch: same (true) - BG or opposite (false) - SIGNAL sign"};
  //  Configurable<float> ptTotcut{"ptTotcut", 0.15, "min pt of all 4 tracks cut"};
  //  Configurable<float> minAnglecut{"minAnglecut", 0.05, "min angle between tracks cut"};
  //  Configurable<float> minNsigmaElcut{"minNsigmaElcut", -2., "min Nsigma for Electrons cut"};
  //  Configurable<float> maxNsigmaElcut{"maxNsigmaElcut", 3., "max Nsigma for Electrons cut"};
  //  Configurable<float> maxNsigmaPiVetocut{"maxNsigmaPiVetocut", 4., "max Nsigma for Pion veto cut"};
  //  Configurable<float> maxNsigmaPrVetocut{"maxNsigmaPrVetocut", 3., "max Nsigma for Proton veto cut"};
  //  Configurable<float> maxNsigmaKaVetocut{"maxNsigmaKaVetocut", 3., "max Nsigma for Kaon veto cut"};
  //  Configurable<float> minPtEtrkcut{"minPtEtrkcut", 0.25, "min Pt for El track cut"};
  Configurable<float> minTrkPtcut{"minTrkPtcut", 0.1, "min Pt for each charged track cut, default=100MeV/c"};
  Configurable<int> nPVtrackscut{"nPVtrackscut", 4, "number of PV contributors, default=4"};
  Configurable<bool> mFITvetoFlag{"mFITvetoFlag", true, "To apply FIT veto"};
  Configurable<int> mFITvetoWindow{"mFITvetoWindow", 2, "FIT veto window"};
  Configurable<bool> useFV0ForVeto{"useFV0ForVeto", 0, "use FV0 for veto"};
  Configurable<bool> useFDDAForVeto{"useFDDAForVeto", 0, "use FDDA for veto"};
  Configurable<bool> useFDDCForVeto{"useFDDCForVeto", 0, "use FDDC for veto"};
  Configurable<int> nTofTrkMinCut{"nTofTrkMinCut", 2, "min TOF tracks"};

  //  Configurable<bool> invMass3piSignalRegion{"invMass3piSignalRegion", 1, "1-use inv mass 3pi in signal region, 0-in background region"};
  //  Configurable<float> invMass3piMaxcut{"invMass3piMaxcut", 1.8, "Z invariant mass of 3 pi cut"};
  //  Configurable<float> deltaPhiMincut{"deltaPhiMincut", 0., "delta phi electron - 3 pi direction cut"};
  //  Configurable<int> nTPCcrossedRowsMinCut{"nTPCcrossedRowsMinCut", 50, "min N_crossed TPC rows for electron candidate"};
  //  Configurable<float> nSigma3piMaxCut{"nSigma3piMaxCut", 5., "n sigma 3 pi max cut"};
  //  Configurable<float> whichPIDCut{"whichPIDCut", 1., "type of PID selection: 1-TPC,2-sigma(TPC+TOF),3-hardcoded ptCut,default=1"};

  //  Configurable<int> generatorIDMC{"generatorIDMC", -1, "MC generator ID"};
  //  Configurable<bool> removeNoTOFrunsInData{"removeNoTOFrunsInData", 1, "1-remove or 0-keep no TOF runs"};
  //  Configurable<float> occupancyCut{"occupancyCut", 10000., "occupancy cut"};

  // adam
  using UDTracksFull = soa::Join<aod::UDTracks, aod::UDTracksPID, aod::UDTracksExtra, aod::UDTracksFlags, aod::UDTracksDCA>;
  // using UDCollisionsFull2 = soa::Join<aod::UDCollisions, aod::UDCollisionsSels, aod::SGCollisions, aod::UDZdcsReduced>; // without occupancy cut
  using UDCollisionsFull2 = soa::Join<aod::UDCollisions, aod::UDCollisionsSels, aod::SGCollisions, aod::UDZdcsReduced, aod::UDCollisionSelExtras>;
  using UDCollisionFull2 = UDCollisionsFull2::iterator;
  // PVContributors
  Filter pVContributorFilter = aod::udtrack::isPVContributor == true;
  using PVTracks = soa::Filtered<UDTracksFull>;

  // roman
  ////using FullUDTracks = soa::Join<aod::UDTracks, aod::UDTracksExtra, aod::UDTracksDCA, aod::UDTracksPID, aod::UDTracksFlags>;
  ////using FullSGUDCollisions = soa::Join<aod::UDCollisions, aod::UDCollisionsSels, aod::UDCollisionSelExtras, aod::SGCollisions, aod::UDZdcsReduced>;
  ////using FullSGUDCollision = FullSGUDCollisions::iterator;
  using FullMCUDTracks = soa::Join<aod::UDTracks, aod::UDTracksExtra, aod::UDTracksDCA, aod::UDTracksPID, aod::UDTracksFlags, aod::UDMcTrackLabels>;
  using FullMCSGUDCollisions = soa::Join<aod::UDCollisions, aod::UDCollisionsSels, aod::UDCollisionSelExtras, aod::SGCollisions, aod::UDMcCollsLabels>;
  // using FullMCSGUDCollision = FullMCSGUDCollisions::iterator;

  // init
  void init(InitContext&)
  {

    // dgcandidates histograms
    // const AxisSpec axisp{100, 0., 5., "#it{p} (GeV/#it{c})"};
    // const AxisSpec axispt{ptAxis, "p_{T} axis"};
    // const AxisSpec axiseta{etaAxis, "#eta - pseudo rapidity axis"};
    // const AxisSpec axiseta{100, -2., 2., "#eta"};
    // const AxisSpec axisdedx{dedxAxis, "dEdx axis"};
    // const AxisSpec axisminv{minvAxis, "invariant mass axis"};
    //    const AxisSpec axisphi{phiAxis, "phi axis"};
    //    const AxisSpec axisav{vectorAxis, "AV axis"};
    //    const AxisSpec axisas{scalarAxis, "AS axis"};
    // const AxisSpec vectorAxis{100, 0., 2., "A_{V}"};
    // const AxisSpec scalarAxis{100, -1., 1., "A_{S}"};
    // const AxisSpec axisZDC{50, -1., 14., "#it{E} (TeV)"};
    // const AxisSpec axisInvMass4trk{160, 0.5, 8.5, "#it{M}^{4trk}_{inv} (GeV/#it{c}^{2})"};
    // const AxisSpec acoAxis{100, 0., 1., "A^{1+3}"};

    // mySetITShitsRule(cutGlobalTrack.cutITShitsRule);

    if (doprocessDoSkim) {
      registrySkim.add("skim/efficiency", ";efficeincy;events", {HistType::kTH1F, {{10, 0., 10.}}});
      registrySkim.get<TH1>(HIST("skim/efficiency"))->GetXaxis()->SetBinLabel(1, "1: All");
      registrySkim.get<TH1>(HIST("skim/efficiency"))->GetXaxis()->SetBinLabel(2, "2: Gap=012");
      registrySkim.get<TH1>(HIST("skim/efficiency"))->GetXaxis()->SetBinLabel(3, "3: Gap=2");
      registrySkim.get<TH1>(HIST("skim/efficiency"))->GetXaxis()->SetBinLabel(4, "4: PVcont=4");
      registrySkim.get<TH1>(HIST("skim/efficiency"))->GetXaxis()->SetBinLabel(5, "5: |#eta^{tr}|<0.9");
      registrySkim.get<TH1>(HIST("skim/efficiency"))->GetXaxis()->SetBinLabel(6, "6: p_{T}^{tr}>0.1");
      registrySkim.get<TH1>(HIST("skim/efficiency"))->GetXaxis()->SetBinLabel(7, "7: N_{TOF}^{tr}>1");
      registrySkim.get<TH1>(HIST("skim/efficiency"))->GetXaxis()->SetBinLabel(8, "8: FIT veto");

      registrySkim.add("skim/gapSide", ";Gap;events", {HistType::kTH1F, {{10, -1., 9.}}});
      registrySkim.add("skim/trueGapSide", ";TrueGap;events", {HistType::kTH1F, {{10, -1., 9.}}});
      registrySkim.add("skim/etaTrk", ";#eta^{trk};events", {HistType::kTH1F, {{100, -1.5, 1.5}}});
      registrySkim.add("skim/ptTrk", ";p_{T}^{trk};events", {HistType::kTH1F, {{100, 0., 5.}}});
      registrySkim.add("skim/phiTrk", ";#phi^{trk};events", {HistType::kTH1F, {{128, -3.2, 3.2}}});
      registrySkim.add("skim/nTof", ";N_{TOFtrk};events", {HistType::kTH1F, {{10, -1., 9.}}});
    }
    if (doprocessMonteCarlo) {
      registrySkim.add("skim/efficiencyMC", ";efficeincy;events", {HistType::kTH1F, {{10, 0., 10.}}});
      registrySkim.get<TH1>(HIST("skim/efficiencyMC"))->GetXaxis()->SetBinLabel(1, "1: All");
      registrySkim.get<TH1>(HIST("skim/efficiencyMC"))->GetXaxis()->SetBinLabel(2, "2: N^{#tau}=2");
      registrySkim.get<TH1>(HIST("skim/efficiencyMC"))->GetXaxis()->SetBinLabel(3, "3: |y^{#tau}| <= 0.9");
      registrySkim.get<TH1>(HIST("skim/efficiencyMC"))->GetXaxis()->SetBinLabel(4, "4: 4 or 6 trk");
      registrySkim.get<TH1>(HIST("skim/efficiencyMC"))->GetXaxis()->SetBinLabel(5, "5: 4 trk");
      registrySkim.get<TH1>(HIST("skim/efficiencyMC"))->GetXaxis()->SetBinLabel(6, "6: 6 trk");
      registrySkim.get<TH1>(HIST("skim/efficiencyMC"))->GetXaxis()->SetBinLabel(7, "7: |#eta^{ch}|<0.9");
      registrySkim.get<TH1>(HIST("skim/efficiencyMC"))->GetXaxis()->SetBinLabel(8, "8: 7+4 trk");
      registrySkim.get<TH1>(HIST("skim/efficiencyMC"))->GetXaxis()->SetBinLabel(9, "9: 7+6 trk");

      registrySkim.add("skim/problemMC", ";problem;events", {HistType::kTH1F, {{10, 0., 10.}}});

      registrySkim.add("skim/nTauMC", ";N_{#tau};events", {HistType::kTH1F, {{10, 0., 10.}}});
      registrySkim.add("skim/tauRapidityMC", ";y_{#tau};events", {HistType::kTH1F, {{100, -2.5, 2.5}}});
      registrySkim.add("skim/tauPhiMC", ";#phi^{#tau};events", {HistType::kTH1F, {{100, 0, 6.4}}});
      registrySkim.add("skim/tauEtaMC", ";#eta^{#tau};events", {HistType::kTH1F, {{100, -2.5, 2.5}}});
      registrySkim.add("skim/tauPtMC", ";p_{T}^{#tau};events", {HistType::kTH1F, {{100, 0, 5.}}});
      registrySkim.add("skim/tauDeltaEtaMC", ";#Delta#eta^{#tau};events ", {HistType::kTH1F, {{100, -5., 5.}}});
      registrySkim.add("skim/tauDeltaPhiMC", ";#Delta#phi^{#tau}(deg.);events", {HistType::kTH1F, {{100, 131., 181}}});
      registrySkim.add("skim/nChPartMC", ";N^{ch. part};events", {HistType::kTH1F, {{10, 0, 10.}}});
      registrySkim.add("skim/daughterPhiMC", ";#phi^{daughter};events", {HistType::kTH1F, {{100, 0, 6.4}}});
      registrySkim.add("skim/daughterEtaMC", ";#eta^{daughter};events", {HistType::kTH1F, {{100, -4., 4.}}});
      registrySkim.add("skim/daughterPtMC", ";p_{T}^{daughter};events", {HistType::kTH1F, {{100, 0, 5.0}}});
    }

    // histos.add("Truth/hTroubles", "Counter of unwanted issues;;Number of  troubles (-)", HistType::kTH1D, {{15, 0.5, 15.5}});

  } // end init

  float rapidity(float energy, float pz)
  // Just a simple function to return track rapidity
  {
    return 0.5 * std::log((energy + pz) / (energy - pz));
  }

  // helper function to calculate delta alpha
  float deltaAlpha(auto particle1, auto particle2)
  {

    TVector3 vtmp(particle1.px(), particle1.py(), particle1.pz());
    TVector3 v1(particle2.px(), particle2.py(), particle2.pz());
    auto angle = v1.Angle(vtmp);

    return angle;
  }

  float calculateDeltaPhi(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> p, ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> p1)
  {
    //    float delta = p.Phi();
    float delta = RecoDecay::constrainAngle(p.Phi());
    // if (delta < 0)
    //   delta += o2::constants::math::TwoPI;
    // if (p1.Phi() < 0)
    //   delta -= (p1.Phi() + o2::constants::math::TwoPI);
    // else
    delta -= RecoDecay::constrainAngle(p1.Phi());
    delta = RecoDecay::constrainAngle(delta);
    // if (delta < 0)
    //   delta += o2::constants::math::TwoPI;
    if (delta > o2::constants::math::PI)
      delta = o2::constants::math::TwoPI - delta;
    return delta;
  }

  float calculateDeltaPhi(float p, float p1)
  {
    float delta = RecoDecay::constrainAngle(p);
    // if (delta < 0)
    //   delta += o2::constants::math::TwoPI;
    // if (p1 < 0)
    //   delta -= (p1 + o2::constants::math::TwoPI);
    // else
    delta -= RecoDecay::constrainAngle(p1);
    delta = RecoDecay::constrainAngle(delta);
    // if (delta < 0)
    //   delta += o2::constants::math::TwoPI;
    if (delta > o2::constants::math::PI)
      delta = o2::constants::math::TwoPI - delta;
    return delta;
  }

  // helper function to calculate scalar asymmetry
  float scalarAsymMC(auto particle1, auto particle2)
  {
    // auto pt1 = pt(particle1.px(), particle1.py());
    auto pt1 = RecoDecay::pt(particle1.px(), particle1.py());
    // auto pt2 = pt(particle2.px(), particle2.py());
    auto pt2 = RecoDecay::pt(particle2.px(), particle2.py());
    auto delta = pt1 - pt2;
    return std::abs(delta) / (pt1 + pt2);
  }

  // helper function to calculate vector asymmetry
  float vectorAsym(auto particle1, auto particle2)
  {
    auto delta = std::sqrt((particle1.px() - particle2.px()) * (particle1.px() - particle2.px()) +
                           (particle1.py() - particle2.py()) * (particle1.py() - particle2.py()));
    auto sum = std::sqrt((particle1.px() + particle2.px()) * (particle1.px() + particle2.px()) +
                         (particle1.py() + particle2.py()) * (particle1.py() + particle2.py()));
    return sum / delta;
  }

  template <typename T>
  int trackCheck(T track)
  {
    // 1
    if (track.hasITS() && !track.hasTPC() && !track.hasTRD() && !track.hasTOF())
      return 0;
    else if (!track.hasITS() && track.hasTPC() && !track.hasTRD() && !track.hasTOF())
      return 1;
    else if (!track.hasITS() && !track.hasTPC() && track.hasTRD() && !track.hasTOF())
      return 2;
    else if (!track.hasITS() && !track.hasTPC() && !track.hasTRD() && track.hasTOF())
      return 3;
    // 2
    else if (track.hasITS() && track.hasTPC() && !track.hasTRD() && !track.hasTOF())
      return 4;
    else if (track.hasITS() && !track.hasTPC() && track.hasTRD() && !track.hasTOF())
      return 5;
    else if (track.hasITS() && !track.hasTPC() && !track.hasTRD() && track.hasTOF())
      return 6;
    else if (!track.hasITS() && track.hasTPC() && track.hasTRD() && !track.hasTOF())
      return 7;
    else if (!track.hasITS() && track.hasTPC() && !track.hasTRD() && track.hasTOF())
      return 8;
    else if (!track.hasITS() && !track.hasTPC() && track.hasTRD() && track.hasTOF())
      return 9;
    // 3
    else if (track.hasITS() && track.hasTPC() && track.hasTRD() && !track.hasTOF())
      return 10;
    else if (track.hasITS() && track.hasTPC() && !track.hasTRD() && track.hasTOF())
      return 11;
    else if (track.hasITS() && !track.hasTPC() && track.hasTRD() && track.hasTOF())
      return 12;
    else if (!track.hasITS() && track.hasTPC() && track.hasTRD() && track.hasTOF())
      return 13;
    // 4
    else if (track.hasITS() && track.hasTPC() && track.hasTRD() && track.hasTOF())
      return 14;
    return -1;
  }

  // to be appied later
  //  // analysis track quality check
  //  template <typename T>
  //  bool isGoodTrackCheck(T track)
  //  {
  //    if (!track.hasTPC())
  //      return false;
  //    if (track.tpcChi2NCl() >= 4.)
  //      return false;
  //    if (track.itsChi2NCl() >= 36.)
  //      return false;
  //    // if (track.dcaZ() >= 2.) return false;
  //    // if (track.dcaXY() >= 0.0105 * 0.035 / std::pow(track.pt(), 1.1)) return false;
  //    if (track.tpcNClsCrossedRows() <= 50)
  //      return false;
  //    if (track.tpcNClsFindable() == 0)
  //      return false;
  //    if (track.tpcNClsCrossedRows() / track.tpcNClsFindable() <= 0.8)
  //      return false;
  //    return true;
  //  }

  //  // analysis track quality check
  //  template <typename T>
  //  bool isGoodTOFTrackCheckHisto(T track)
  //  {
  //    bool isGoodTrack = true;
  //    if (track.hasTOF()) {
  //      registry.get<TH1>(HIST("global/hTrackPVGood"))->Fill(8., 1.);
  //    } else {
  //      isGoodTrack = false;
  //    }
  //    if (track.hasTOF() && track.tofChi2() < 3) {
  //      registry.get<TH1>(HIST("global/hTrackPVGood"))->Fill(9., 1.);
  //    } else {
  //      isGoodTrack = false;
  //    }
  //    return isGoodTrack;
  //  }

  //  // analysis track quality check
  //  template <typename T>
  //  bool isGoodTOFTrackCheck(T track)
  //  {
  //    if (!track.hasTOF())
  //      return false;
  //    if (track.tofChi2() >= 3)
  //      return false;
  //    return true;
  //  }

  // check ITS clusters, how many -1,0,1,7 + 10 if 0,1,2 layers were fired
  // analysis track quality check
  template <typename T>
  int numberOfItsClustersCheck(T track)
  {
    if (!track.hasITS())
      return -1;
    int nITSbits = 0;
    int firstThreeLayers = 0;
    const int threeLayers = 3;
    const int maxITSlayers = 7;
    uint32_t clusterSizes = track.itsClusterSizes();
    for (int layer = 0; layer < maxITSlayers; layer++) {
      if ((clusterSizes >> (layer * 4)) & 0xf) {
        nITSbits++;
        if (layer < threeLayers) // 3
          firstThreeLayers++;
      }
    } // end of loop over ITS bits
    if (firstThreeLayers == threeLayers) // 3
      nITSbits += 10;

    return nITSbits;
  }

  // RCT check
  template <typename C>
  int isGoodRCTflag(C const& coll)
  {
    if (sgSelector.isCBTHadronZdcOk(coll))
      return 4;
    else if (sgSelector.isCBTHadronOk(coll))
      return 3;
    else if (sgSelector.isCBTZdcOk(coll))
      return 2;
    else if (sgSelector.isCBTOk(coll))
      return 1;
    else
      return 0;
  }

  //////////////////////////////////////////

  //  template <typename C>
  //  bool isGoodFITtime(C const& coll, float maxFITtime)
  //  {
  //
  //    // FTOA
  //    if ((std::abs(coll.timeFT0A()) > maxFITtime) && coll.timeFT0A() > -998.)
  //      return false;
  //
  //    // FTOC
  //    if ((std::abs(coll.timeFT0C()) > maxFITtime) && coll.timeFT0C() > -998.)
  //      return false;
  //
  //    return true;
  //  }

  //   template <typename C>
  //   bool isGoodROFtime(C const& coll)
  //   {
  //
  //     // kNoTimeFrameBorder
  //     if (cutSample.cutEvTFb && !coll.tfb())
  //       return false;
  //
  //     // kNoITSROFrameBorder
  //     if (cutSample.cutEvITSROFb && !coll.itsROFb())
  //       return false;
  //
  //     // kNoSameBunchPileup
  //     if (cutSample.cutEvSbp && !coll.sbp())
  //       return false;
  //
  //     // kIsGoodZvtxFT0vsPV
  //     if (cutSample.cutEvZvtxFT0vPV && !coll.zVtxFT0vPV())
  //       return false;
  //
  //     // kIsVertexITSTPC
  //     if (cutSample.cutEvVtxITSTPC && !coll.vtxITSTPC())
  //       return false;
  //
  //     // Occupancy
  //     if (coll.occupancyInTime() > cutSample.cutEvOccupancy)
  //       return false;
  //
  //     // kNoCollInTimeRangeStandard
  //     if (cutSample.cutEvTrs && !coll.trs())
  //       return false;
  //
  //     // kNoCollInRofStandard
  //     if (cutSample.cutEvTrofs && !coll.trofs())
  //       return false;
  //
  //     // kNoHighMultCollInPrevRof
  //     if (cutSample.cutEvHmpr && !coll.hmpr())
  //       return false;
  //
  //     return true;
  //   }

  std::vector<std::pair<int8_t, std::set<uint8_t>>> cutMyRequiredITSHits{};

  void mySetRequireHitsInITSLayers(int8_t minNRequiredHits, std::set<uint8_t> requiredLayers)
  {
    // layer 0 corresponds to the the innermost ITS layer
    cutMyRequiredITSHits.push_back(std::make_pair(minNRequiredHits, requiredLayers));
  }

  void mySetITShitsRule(int matching)
  {
    switch (matching) {
      case 0: // Run3ITSibAny
        mySetRequireHitsInITSLayers(1, {0, 1, 2});
        break;
      case 1: // Run3ITSibTwo
        mySetRequireHitsInITSLayers(2, {0, 1, 2});
        break;
      case 2: // Run3ITSallAny
        mySetRequireHitsInITSLayers(1, {0, 1, 2, 3, 4, 5, 6});
        break;
      case 3: // Run3ITSall7Layers
        mySetRequireHitsInITSLayers(7, {0, 1, 2, 3, 4, 5, 6});
        break;
      default:
        LOG(fatal) << "You chose wrong ITS matching";
        break;
    }
  }

  bool isFulfillsITSHitRequirementsReinstatement(uint8_t itsClusterMap) const
  {
    constexpr uint8_t KBit = 1;
    for (const auto& kITSrequirement : cutMyRequiredITSHits) {
      auto hits = std::count_if(kITSrequirement.second.begin(), kITSrequirement.second.end(), [&](auto&& requiredLayer) { return itsClusterMap & (KBit << requiredLayer); });
      if ((kITSrequirement.first == -1) && (hits > 0)) {
        return false; // no hits were required in specified layers
      } else if (hits < kITSrequirement.first) {
        return false; // not enough hits found in specified layers
      }
    }
    return true;
  }

  //  template <typename T>
  //  bool isGlobalTrackReinstatement(T const& track)
  //  {
  //    // kInAcceptance copy
  //    if (track.pt() < cutGlobalTrack.cutMinPt || track.pt() > cutGlobalTrack.cutMaxPt)
  //      return false;
  //    if (eta(track.px(), track.py(), track.pz()) < cutGlobalTrack.cutMinEta || eta(track.px(), track.py(), track.pz()) > cutGlobalTrack.cutMaxEta)
  //      return false;
  //    // kPrimaryTracks
  //    // GoldenChi2 cut is only for Run 2
  //    if (std::abs(track.dcaZ()) > cutGlobalTrack.cutMaxDCAz)
  //      return false;
  //    if (cutGlobalTrack.applyPtDependentDCAxy) {
  //      float maxDCA = 0.0182f + 0.0350f / std::pow(track.pt(), 1.01f);
  //      if (std::abs(track.dcaXY()) > maxDCA)
  //        return false;
  //    } else {
  //      if (std::abs(track.dcaXY()) > cutGlobalTrack.cutMaxDCAxy)
  //        return false;
  //    }
  //    // kQualityTrack
  //    // TrackType is always 1 as per definition of processed Run3 AO2Ds
  //    // ITS
  //    if (cutGlobalTrack.cutHasITS && !track.hasITS())
  //      return false; // ITS refit
  //    if (track.itsNCls() < cutGlobalTrack.cutMinITSnCls)
  //      return false;
  //    if (track.itsChi2NCl() > cutGlobalTrack.cutMaxITSchi2)
  //      return false;
  //    if (!isFulfillsITSHitRequirementsReinstatement(track.itsClusterMap()))
  //      return false;
  //    //  TPC
  //    if (cutGlobalTrack.cutHasTPC && !track.hasTPC())
  //      return false; // TPC refit
  //    if ((track.tpcNClsFindable() - track.tpcNClsFindableMinusFound()) < cutGlobalTrack.cutMinTPCnCls)
  //      return false; // tpcNClsFound()
  //    if (track.tpcNClsCrossedRows() < cutGlobalTrack.cutMinTPCnClsXrows)
  //      return false;
  //    if ((static_cast<float>(track.tpcNClsCrossedRows()) / static_cast<float>(track.tpcNClsFindable())) < cutGlobalTrack.cutMinTPCnClsXrowsOverNcls)
  //      return false;
  //    if (track.tpcChi2NCl() > cutGlobalTrack.cutMaxTPCchi2)
  //      return false; // TPC chi2
  //    if (cutGlobalTrack.cutGoodITSTPCmatching) {
  //      if (track.itsChi2NCl() < 0.)
  //        return false; // good ITS-TPC matching means ITS ch2 is not negative
  //    }
  //    //  TOF
  //    if (track.hasTOF()) {
  //      if (track.tpcChi2NCl() > cutGlobalTrack.cutMaxTOFchi2)
  //        return false; // TOF chi2
  //    }
  //
  //    return true;
  //  }

  //  template <typename T>
  //  bool isElectronCandidate(T const& electronCandidate)
  //  // Loose criterium to find electron-like particle
  //  // Requiring TOF to avoid double-counting pions/electrons and for better timing
  //  {
  //    if (electronCandidate.tpcNSigmaEl() < cutPreselect.preselMaxElectronNsigmaEl || electronCandidate.tpcNSigmaEl() > cutPreselect.preselMinElectronNsigmaEl)
  //      return false;
  //    if (cutPreselect.preselElectronHasTOF && !electronCandidate.hasTOF())
  //      return false;
  //    return true;
  //  }
  //
  //  template <typename T>
  //  bool isMuPionCandidate(T const& muPionCandidate)
  //  // Loose criterium to find muon/pion-like particle
  //  // Requiring TOF for better timing
  //  {
  //    if (muPionCandidate.tpcNSigmaMu() < cutPreselect.preselMaxMuonNsigmaEl || muPionCandidate.tpcNSigmaMu() > cutPreselect.preselMinMuonNsigmaEl)
  //      return false;
  //    if (muPionCandidate.tpcNSigmaPi() < cutPreselect.preselMaxPionNsigmaEl || muPionCandidate.tpcNSigmaPi() > cutPreselect.preselMinPionNsigmaEl)
  //      return false;
  //    if (cutPreselect.preselMupionHasTOF && !muPionCandidate.hasTOF())
  //      return false;
  //    return true;
  //  }

  int enumMyParticle(int valuePDG)
  // reads pdg value and returns particle number as in enumMyParticle
  {
    if (std::abs(valuePDG) == kElectron) {         // 11 e+ or e-
      return MyElectron;                           // 1
    } else if (std::abs(valuePDG) == kMuonMinus) { // 13 mu+ or mu -
      return MyMuon;                               // 2
    } else if (std::abs(valuePDG) == kPiPlus) {    // 211 pi+ (or pi-)
      return MyPion;                               // 3
    } else if (std::abs(valuePDG) == kKPlus) {     // 321 K+ (or K-)
      return MyKaon;                               // 4
    } else {
      if (verbose)
        LOGF(info, "PDG value not found in enumMyParticle. Returning -1.");
      return MyOtherParticle; // -1
    }
  }

  // skimming: only 4 tracks selection in data
  void processDoSkim(UDCollisionFull2 const& dgcand, UDTracksFull const& dgtracks, PVTracks const& PVContributors)
  {
    registrySkim.get<TH1>(HIST("skim/efficiency"))->Fill(0., 1.);

    int gapSide = dgcand.gapSide();
    registrySkim.get<TH1>(HIST("skim/gapSide"))->Fill(gapSide);
    //    if (gapSide < 0 || gapSide > 2)
    if (gapSide < o2::aod::sgselector::SingleGapA || gapSide > o2::aod::sgselector::DoubleGap)
      return;
    registrySkim.get<TH1>(HIST("skim/efficiency"))->Fill(1., 1.);

    int truegapSide = sgSelector.trueGap(dgcand, cutFV0, cutFT0A, cutFT0C, cutZDC);
    gapSide = truegapSide;
    registrySkim.get<TH1>(HIST("skim/trueGapSide"))->Fill(gapSide);
    if (gapSide != mGapSide)
      return;
    registrySkim.get<TH1>(HIST("skim/efficiency"))->Fill(2., 1.);

    // zdc information
    float energyZNA = dgcand.energyCommonZNA();
    float energyZNC = dgcand.energyCommonZNC();
    if (energyZNA < 0)
      energyZNA = -1.;
    if (energyZNC < 0)
      energyZNC = -1.;

    int nTofTrk = 0;
    int nEtaIn15 = 0;
    int npT100 = 0;
    // //   int qtot = 0;
    // int8_t qtot = 0;
    // TLorentzVector p;
    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> p;
    for (const auto& trk : PVContributors) {
      // qtot += trk.sign();
      // p.SetXYZM(trk.px(), trk.py(), trk.pz(), MassPiPlus);
      p.SetXYZT(trk.px(), trk.py(), trk.pz(), RecoDecay::e(trk.px(), trk.py(), trk.pz(), MassPiPlus));
      registrySkim.get<TH1>(HIST("skim/etaTrk"))->Fill(p.Eta());
      registrySkim.get<TH1>(HIST("skim/ptTrk"))->Fill(trk.pt());
      registrySkim.get<TH1>(HIST("skim/phiTrk"))->Fill(p.Phi());

      if (std::abs(p.Eta()) < trkEtacut)
        nEtaIn15++;               // 0.9 is a default
      if (trk.pt() > minTrkPtcut) // 0.1 GeV/c
        npT100++;
      if (trk.hasTOF())
        nTofTrk++;
    } // end of loop over PV tracks

    if (PVContributors.size() != nPVtrackscut) // 4
      return;
    registrySkim.get<TH1>(HIST("skim/efficiency"))->Fill(3., 1.);

    if (nEtaIn15 != nPVtrackscut) // 4
      return;
    registrySkim.get<TH1>(HIST("skim/efficiency"))->Fill(4., 1.);

    if (npT100 != nPVtrackscut) // 4
      return;
    registrySkim.get<TH1>(HIST("skim/efficiency"))->Fill(5., 1.);

    registrySkim.get<TH1>(HIST("skim/nTof"))->Fill(nTofTrk);
    if (nTofTrk < nTofTrkMinCut)
      return;
    registrySkim.get<TH1>(HIST("skim/efficiency"))->Fill(6., 1.);

    //
    // FIT informaton
    //
    auto bitMin = 16 - mFITvetoWindow; // default is +- 2 bc (2 bit)
    auto bitMax = 16 + mFITvetoWindow;
    bool flagFITveto = false;
    // check FIT information
    for (auto bit = bitMin; bit <= bitMax; bit++) {
      if (TESTBIT(dgcand.bbFT0Apf(), bit))
        flagFITveto = true;
      if (TESTBIT(dgcand.bbFT0Cpf(), bit))
        flagFITveto = true;
      if (useFV0ForVeto && TESTBIT(dgcand.bbFV0Apf(), bit))
        flagFITveto = true;
      if (useFDDAForVeto && TESTBIT(dgcand.bbFDDApf(), bit))
        flagFITveto = true;
      if (useFDDCForVeto && TESTBIT(dgcand.bbFDDCpf(), bit))
        flagFITveto = true;
    } // end of loop over FIT bits

    if (mFITvetoFlag && flagFITveto)
      return;
    registrySkim.get<TH1>(HIST("skim/efficiency"))->Fill(7., 1.);

    // RCT variable
    int rct = 0;
    rct = isGoodRCTflag(dgcand);

    //
    // variables per track
    //
    int counterTmp = 0;
    float px[6] = {-999., -999., -999., -999., -999., -999.};
    float py[6] = {-999., -999., -999., -999., -999., -999.};
    float pz[6] = {-999., -999., -999., -999., -999., -999.};
    int sign[6] = {-99, -99, -99, -99, -99, -99};
    float dcaZ[6] = {-999., -999., -999., -999., -999., -999.};
    float dcaXY[6] = {-999., -999., -999., -999., -999., -999.};

    float tmpDedx[6] = {-999., -999., -999., -999., -999., -999.};
    float tmpTofNsigmaEl[6] = {-999., -999., -999., -999., -999., -999.};
    float tmpTofNsigmaPi[6] = {-999., -999., -999., -999., -999., -999.};
    float tmpTofNsigmaKa[6] = {-999., -999., -999., -999., -999., -999.};
    float tmpTofNsigmaPr[6] = {-999., -999., -999., -999., -999., -999.};
    float tmpTofNsigmaMu[6] = {-999., -999., -999., -999., -999., -999.};

    float nSigmaEl[6] = {-999., -999., -999., -999., -999., -999.};
    float nSigmaPi[6] = {-999., -999., -999., -999., -999., -999.};
    float nSigmaPr[6] = {-999., -999., -999., -999., -999., -999.};
    float nSigmaKa[6] = {-999., -999., -999., -999., -999., -999.};
    float nSigmaMu[6] = {-999., -999., -999., -999., -999., -999.};
    float chi2TOF[6] = {-999., -999., -999., -999., -999., -999.};
    int nclTPCcrossedRows[6] = {-999, -999, -999, -999, -999, -999};
    int nclTPCfind[6] = {-999, -999, -999, -999, -999, -999};
    float nclTPCchi2[6] = {-999., -999., -999., -999., -999., -999.};
    float trkITSchi2[6] = {-999., -999., -999., -999., -999., -999.};
    int trkITScl[6] = {-999, -999, -999, -999, -999, -999};

    // double trkTime[4];
    // float trkTimeRes[4];
    float trkTofSignal[6] = {-999., -999., -999., -999., -999., -999.};

    for (const auto& trk : PVContributors) {
      if (counterTmp > nPVtrackscut - 1)
        continue; // >5 -> continue default
      px[counterTmp] = trk.px();
      py[counterTmp] = trk.py();
      pz[counterTmp] = trk.pz();
      sign[counterTmp] = trk.sign();
      dcaZ[counterTmp] = trk.dcaZ();
      dcaXY[counterTmp] = trk.dcaXY();

      tmpDedx[counterTmp] = trk.tpcSignal();
      nSigmaEl[counterTmp] = trk.tpcNSigmaEl();
      nSigmaPi[counterTmp] = trk.tpcNSigmaPi();
      nSigmaPr[counterTmp] = trk.tpcNSigmaPr();
      nSigmaKa[counterTmp] = trk.tpcNSigmaKa();
      nSigmaMu[counterTmp] = trk.tpcNSigmaMu();

      trkTofSignal[counterTmp] = trk.beta();
      tmpTofNsigmaEl[counterTmp] = trk.tofNSigmaEl();
      tmpTofNsigmaPi[counterTmp] = trk.tofNSigmaPi();
      tmpTofNsigmaKa[counterTmp] = trk.tofNSigmaKa();
      tmpTofNsigmaPr[counterTmp] = trk.tofNSigmaPr();
      tmpTofNsigmaMu[counterTmp] = trk.tofNSigmaMu();

      if (trk.hasTOF())
        chi2TOF[counterTmp] = trk.tofChi2();

      nclTPCcrossedRows[counterTmp] = trk.tpcNClsCrossedRows();
      nclTPCfind[counterTmp] = trk.tpcNClsFindable();
      nclTPCchi2[counterTmp] = trk.tpcChi2NCl();
      trkITSchi2[counterTmp] = trk.itsChi2NCl();
      trkITScl[counterTmp] = numberOfItsClustersCheck(trk);

      // trkTime[counterTmp] = trk.trackTime();
      // trkTimeRes[counterTmp] = trk.trackTimeRes();
      counterTmp++;
    }

    dataTauFourTracks(dgcand.runNumber(),
                      dgcand.globalBC(), // is it necessary
                      dgtracks.size(),
                      dgcand.numContrib(),
                      rct,
                      // dgcand.posX(), dgcand.posY(),
                      dgcand.posZ(),
                      dgcand.flags(),
                      dgcand.occupancyInTime(),
                      dgcand.hadronicRate(),                       // is it necessary
                      dgcand.trs(), dgcand.trofs(), dgcand.hmpr(), // to test it
                      dgcand.tfb(), dgcand.itsROFb(), dgcand.sbp(), dgcand.zVtxFT0vPV(), dgcand.vtxITSTPC(),
                      energyZNA, energyZNC,
                      // qtot, <<-------- comment out
                      dgcand.totalFT0AmplitudeA(), dgcand.totalFT0AmplitudeC(), dgcand.totalFV0AmplitudeA(),
                      // dgcand.timeFT0A(), dgcand.timeFT0C(), dgcand.timeFV0A(),
                      px, py, pz, sign,
                      dcaXY, dcaZ,
                      nclTPCcrossedRows, nclTPCfind, nclTPCchi2, trkITSchi2, trkITScl,
                      tmpDedx, nSigmaEl, nSigmaPi, nSigmaKa, nSigmaPr, nSigmaMu,
                      trkTofSignal, tmpTofNsigmaEl, tmpTofNsigmaPi, tmpTofNsigmaKa, tmpTofNsigmaPr, tmpTofNsigmaMu,
                      chi2TOF);
  } // end of skim process processDoSkim
  PROCESS_SWITCH(TauThreeProngEventTableProducer, processDoSkim, "Run over SG Producer tables to produce skimmed data", false);

  //  void processDataSG(FullSGUDCollision const& collision,
  //                     FullUDTracks const& tracks)
  //  {
  //
  //    if (!isGoodRCTflag(collision))
  //      return;
  //
  //    if (!isGoodROFtime(collision))
  //      return;
  //
  //    int gapSide = collision.gapSide();
  //    int trueGapSide = sgSelector.trueGap(collision, cutSample.cutTrueGapSideFV0, cutSample.cutTrueGapSideFT0A, cutSample.cutTrueGapSideFT0C, cutSample.cutTrueGapSideZDC);
  //    if (cutSample.useTrueGap)
  //      gapSide = trueGapSide;
  //    if (gapSide != cutSample.whichGapSide)
  //      return;
  //
  //    if (!isGoodFITtime(collision, cutSample.cutFITtime))
  //      return;
  //
  //    if (cutSample.useNumContribs && (collision.numContrib() != cutSample.cutNumContribs))
  //      return;
  //
  //    if (cutSample.useRecoFlag && (collision.flags() != cutSample.cutRecoFlag))
  //      return;
  //
  //    int countTracksPerCollision = 0;
  //    int countGoodNonPVtracks = 0;
  //    int countGoodPVtracks = 0;
  //    std::vector<int> vecTrkIdx;
  //    // Loop over tracks with selections
  //    for (const auto& track : tracks) {
  //      countTracksPerCollision++;
  //      if (!isGlobalTrackReinstatement(track))
  //        continue;
  //      if (!track.isPVContributor()) {
  //        countGoodNonPVtracks++;
  //        continue;
  //      }
  //      countGoodPVtracks++;
  //      vecTrkIdx.push_back(track.index());
  //    } // Loop over tracks with selections
  //
  //    // Apply weak condition on track PID
  //    int countPVGTel = 0;
  //    int countPVGTmupi = 0;
  //    if (countGoodPVtracks == 2) {
  //      for (const auto& vecMember : vecTrkIdx) {
  //        const auto& thisTrk = tracks.iteratorAt(vecMember);
  //        if (isElectronCandidate(thisTrk)) {
  //          countPVGTel++;
  //          continue;
  //        }
  //        if (isMuPionCandidate(thisTrk)) {
  //          countPVGTmupi++;
  //        }
  //      }
  //    }
  //
  //    if (cutPreselect.preselUseTrackPID ? ((countPVGTel == 2 && countPVGTmupi == 0) || (countPVGTel == 1 && countPVGTmupi == 1)) : countGoodPVtracks == cutPreselect.preselNgoodPVtracs) {
  //      const auto& trk1 = tracks.iteratorAt(vecTrkIdx[0]);
  //      const auto& trk2 = tracks.iteratorAt(vecTrkIdx[1]);
  //
  //      float px[2] = {trk1.px(), trk2.px()};
  //      float py[2] = {trk1.py(), trk2.py()};
  //      float pz[2] = {trk1.pz(), trk2.pz()};
  //      int sign[2] = {trk1.sign(), trk2.sign()};
  //      float dcaxy[2] = {trk1.dcaXY(), trk2.dcaXY()};
  //      float dcaz[2] = {trk1.dcaZ(), trk2.dcaZ()};
  //      float trkTimeRes[2] = {trk1.trackTimeRes(), trk2.trackTimeRes()};
  //      uint32_t itsClusterSizesTrk1 = trk1.itsClusterSizes();
  //      uint32_t itsClusterSizesTrk2 = trk2.itsClusterSizes();
  //      float tpcSignal[2] = {trk1.tpcSignal(), trk2.tpcSignal()};
  //      float tpcEl[2] = {trk1.tpcNSigmaEl(), trk2.tpcNSigmaEl()};
  //      float tpcMu[2] = {trk1.tpcNSigmaMu(), trk2.tpcNSigmaMu()};
  //      float tpcPi[2] = {trk1.tpcNSigmaPi(), trk2.tpcNSigmaPi()};
  //      float tpcKa[2] = {trk1.tpcNSigmaKa(), trk2.tpcNSigmaKa()};
  //      float tpcPr[2] = {trk1.tpcNSigmaPr(), trk2.tpcNSigmaPr()};
  //      float tpcIP[2] = {trk1.tpcInnerParam(), trk2.tpcInnerParam()};
  //      float tofSignal[2] = {trk1.tofSignal(), trk2.tofSignal()};
  //      float tofEl[2] = {trk1.tofNSigmaEl(), trk2.tofNSigmaEl()};
  //      float tofMu[2] = {trk1.tofNSigmaMu(), trk2.tofNSigmaMu()};
  //      float tofPi[2] = {trk1.tofNSigmaPi(), trk2.tofNSigmaPi()};
  //      float tofKa[2] = {trk1.tofNSigmaKa(), trk2.tofNSigmaKa()};
  //      float tofPr[2] = {trk1.tofNSigmaPr(), trk2.tofNSigmaPr()};
  //      float tofEP[2] = {trk1.tofExpMom(), trk2.tofExpMom()};
  //      //      float infoZDC[4] = {-999., -999., -999., -999.};
  //      //      if constexpr (requires { collision.udZdcsReduced(); }) {
  //      //        infoZDC[0] = collision.energyCommonZNA();
  //      //        infoZDC[1] = collision.energyCommonZNC();
  //      //        infoZDC[2] = collision.timeZNA();
  //      //        infoZDC[3] = collision.timeZNC();
  //      //      }
  //      float infoZDC[4] = {collision.energyCommonZNA(), collision.energyCommonZNC(), collision.timeZNA(), collision.timeZNC()};
  //
  //      tauTwoTracks(collision.runNumber(), collision.globalBC(), countTracksPerCollision, collision.numContrib(), countGoodNonPVtracks, collision.posX(), collision.posY(), collision.posZ(),
  //                   collision.flags(), collision.occupancyInTime(), collision.hadronicRate(), collision.trs(), collision.trofs(), collision.hmpr(),
  //                   collision.tfb(), collision.itsROFb(), collision.sbp(), collision.zVtxFT0vPV(), collision.vtxITSTPC(),
  //                   collision.totalFT0AmplitudeA(), collision.totalFT0AmplitudeC(), collision.totalFV0AmplitudeA(), infoZDC[0], infoZDC[1],
  //                   collision.timeFT0A(), collision.timeFT0C(), collision.timeFV0A(), infoZDC[2], infoZDC[3],
  //                   px, py, pz, sign, dcaxy, dcaz, trkTimeRes,
  //                   itsClusterSizesTrk1, itsClusterSizesTrk2,
  //                   tpcSignal, tpcEl, tpcMu, tpcPi, tpcKa, tpcPr, tpcIP,
  //                   tofSignal, tofEl, tofMu, tofPi, tofKa, tofPr, tofEP);
  //    }
  //  }
  //  PROCESS_SWITCH(TauEventTableProducer, processDataSG, "Iterate UD tables with measured data created by SG-Candidate-Producer.", false);

  //  PresliceUnsorted<aod::UDMcParticles> partPerMcCollision = aod::udmcparticle::udMcCollisionId;
  //  PresliceUnsorted<FullMCSGUDCollisions> colPerMcCollision = aod::udcollision::udMcCollisionId;
  //  PresliceUnsorted<FullMCUDTracks> trackPerMcParticle = aod::udmctracklabel::udMcParticleId;
  //  Preslice<FullMCUDTracks> trackPerCollision = aod::udtrack::udCollisionId; // sorted preslice used because the pair track-collision is already sorted in processDataSG function
  //
  //  void processMonteCarlo(aod::UDMcCollisions const& mccollisions,
  //                         aod::UDMcParticles const& parts,
  //                         FullMCSGUDCollisions const& recolls,
  //                         FullMCUDTracks const& trks)

  PresliceUnsorted<aod::UDMcParticles> partPerMcCollision = aod::udmcparticle::udMcCollisionId;
  PresliceUnsorted<FullMCSGUDCollisions> colPerMcCollision = aod::udcollision::udMcCollisionId;
  PresliceUnsorted<FullMCUDTracks> trackPerMcParticle = aod::udmctracklabel::udMcParticleId;
  Preslice<FullMCUDTracks> trackPerCollision = aod::udtrack::udCollisionId; // sorted preslice used because the pair track-collision is already sorted in processDataSG function

  void processMonteCarlo(aod::UDMcCollisions const& mcCollisions,
                         aod::UDMcParticles const& mcParticles,
                         FullMCSGUDCollisions const& collisions,
                         FullMCUDTracks const& tracks)
  {
    // registrySkim.get<TH1>(HIST("skim/efficiencyMC"))->Fill(0., 1.);

    const int fourTracks = 4;
    const int sixTracks = 6;
    const int oneProng = 1;
    const int threeProng = 3;

    if (verbose)
      LOGF(info, "0. <MC> UDMcCollision size %d, Collisions size %d, UDtracks %d, UDMcParticles %d", mcCollisions.size(), collisions.size(), tracks.size(), mcParticles.size());

    // temporary variables
    float tmpRapidity = -999.;
    float trueTauEta[2] = {-999., -999.};
    float trueTauPhi[2] = {-999., -999.};

    // init variables for tree
    float trueTauX[2] = {-999., -999.};
    float trueTauY[2] = {-999., -999.};
    float trueTauZ[2] = {-999., -999.};

    bool tauInRapidity = true;
    bool partFromTauInEta = true;

    // start loop over generated collisions
    for (const auto& mccoll : mcCollisions) {
      registrySkim.get<TH1>(HIST("skim/efficiencyMC"))->Fill(0., 1.); // all MC collisions

      // set up default values per colission
      trueTauX[0] = -999.;
      trueTauY[0] = -999.;
      trueTauZ[0] = -999.;
      trueTauX[1] = -999.;
      trueTauY[1] = -999.;
      trueTauZ[1] = -999.;

      tauInRapidity = true;
      partFromTauInEta = true;

      // get particles associated to generated collision
      auto const& tmpPartsFromMcColl = mcParticles.sliceBy(partPerMcCollision, mccoll.globalIndex());
      if (verbose)
        LOGF(info, "1. part from MC coll %d", tmpPartsFromMcColl.size());
      int countMothers = 0;
      const int desiredNMothers = 2;
      for (const auto& particle : tmpPartsFromMcColl) {
        if (verbose)
          LOGF(info, "2. MC part pdg %d", particle.pdgCode());
        if (std::abs(particle.pdgCode()) != kTauMinus)
          continue; // 15 = tau_minus
        // if (std::abs(particle.pdgCode()) != 15) continue; // 15 = tau_minus
        if (countMothers < desiredNMothers) { // < 2
          // fill info for each tau
          trueTauX[countMothers] = particle.px();
          trueTauY[countMothers] = particle.py();
          trueTauZ[countMothers] = particle.pz();
          tmpRapidity = rapidity(particle.e(), trueTauZ[countMothers]);
          trueTauEta[countMothers] = RecoDecay::eta(std::array<double, 3>{particle.px(), particle.py(), particle.pz()});
          trueTauPhi[countMothers] = RecoDecay::phi(particle.px(), particle.py());

          if (verbose)
            LOGF(info, "tau P(%f,%f,%f), e %f, y %f", particle.px(), particle.py(), particle.pz(), particle.e(), tmpRapidity);
          registrySkim.get<TH1>(HIST("skim/tauRapidityMC"))->Fill(tmpRapidity);
          registrySkim.get<TH1>(HIST("skim/tauPhiMC"))->Fill(trueTauPhi[countMothers]);
          registrySkim.get<TH1>(HIST("skim/tauEtaMC"))->Fill(trueTauEta[countMothers]);
          registrySkim.get<TH1>(HIST("skim/tauPtMC"))->Fill(RecoDecay::pt(particle.px(), particle.py()));
          if (std::abs(tmpRapidity) > trkEtacut) { // 0.9
            tauInRapidity = false;
            if (verbose)
              LOGF(info, "tau y %f", tmpRapidity);
          } // rapidity check
        } // number of taus
        countMothers++;
      } // end of loop over MC paricles
      registrySkim.get<TH1>(HIST("skim/nTauMC"))->Fill(countMothers);
      if (countMothers != desiredNMothers) { // 2
        if (verbose)
          LOGF(info, "Truth collision has number of mother particles (taus) %d different than 2. Jump to the next MC event.", countMothers);
        continue;
      }

      registrySkim.get<TH1>(HIST("skim/efficiencyMC"))->Fill(1., 1.); // exactly 2 taus

      if (!tauInRapidity) { // tau NOT in rapidity -> continue
        if (verbose)
          LOGF(info, "At least one mother particle (taus) out of rapidity (|y|<0.9). Jump to the next MC event.");
        continue;
      }

      // delta eta and delta phi between taus
      registrySkim.get<TH1>(HIST("skim/tauDeltaEtaMC"))->Fill(trueTauEta[0] - trueTauEta[1]);
      registrySkim.get<TH1>(HIST("skim/tauDeltaPhiMC"))->Fill(calculateDeltaPhi(trueTauPhi[0], trueTauPhi[1]) * 180. / o2::constants::math::PI);

      registrySkim.get<TH1>(HIST("skim/efficiencyMC"))->Fill(2., 1.); // |y_tau| <= 0.9
      countMothers = 0;
      int nChargedDaughtersTau[2] = {0, 0};
      int nElec = 0;
      int nMuon = 0;
      int nPi = 0;
      int particleType = -1;
      int zerothTau = -10;
      int trueChannel = -1;
      int countPi0 = -1;

      for (const auto& particle : tmpPartsFromMcColl) {
        if (std::abs(particle.pdgCode()) != kTauMinus)
          continue; // 15 = tau_minus
        const auto& daughters = particle.daughters_as<aod::UDMcParticles>();
        for (const auto& daughter : daughters) {
          particleType = enumMyParticle(daughter.pdgCode());
          if (particleType == MyOtherParticle) { // -1
            continue;
          } else {
            nChargedDaughtersTau[countMothers]++;
            if (particleType == MyElectron) // 1
              nElec++;
            else if (particleType == MyMuon) // 2
              nMuon++;
            else if (particleType == MyPion) // 3
              nPi++;
          }

          if (std::abs(RecoDecay::eta(std::array<double, 3>{daughter.px(), daughter.py(), daughter.pz()})) > trkEtacut) // 0.9
            partFromTauInEta = false;
          registrySkim.get<TH1>(HIST("skim/daughterPhiMC"))->Fill(RecoDecay::phi(daughter.px(), daughter.py()));
          registrySkim.get<TH1>(HIST("skim/daughterEtaMC"))->Fill(RecoDecay::eta(std::array<double, 3>{daughter.px(), daughter.py(), daughter.pz()}));
          registrySkim.get<TH1>(HIST("skim/daughterPtMC"))->Fill(RecoDecay::pt(daughter.px(), daughter.py()));
        }
        countMothers++;
        if (countMothers >= desiredNMothers) // 2
          break;
      } // end of loop over MC particles

      registrySkim.get<TH1>(HIST("skim/nChPartMC"))->Fill(nChargedDaughtersTau[0] + nChargedDaughtersTau[1]); // N charged particles from taus
      // check number of charged particles in MC event
      if ((nChargedDaughtersTau[0] + nChargedDaughtersTau[1] != fourTracks) && (nChargedDaughtersTau[0] + nChargedDaughtersTau[1] != sixTracks)) {
        if (verbose)
          LOGF(info, "Different from 4/6 charged particles (%d) from both taus. Jump to the next MC event.", nChargedDaughtersTau[0] + nChargedDaughtersTau[1]);
        continue;
      }
      registrySkim.get<TH1>(HIST("skim/efficiencyMC"))->Fill(3., 1.);        // 1+3 (3+3) topology
      if (nChargedDaughtersTau[0] + nChargedDaughtersTau[1] == fourTracks) { // 4
        registrySkim.get<TH1>(HIST("skim/efficiencyMC"))->Fill(4., 1.);
      } else if (nChargedDaughtersTau[0] + nChargedDaughtersTau[1] == sixTracks) { // 6
        registrySkim.get<TH1>(HIST("skim/efficiencyMC"))->Fill(5., 1.);
      }

      if (!partFromTauInEta) {
        if (verbose)
          LOGF(info, "At least one daughter particle from taus out of pseudo-rapidity (|eta|<0.9). Jump to the next MC event.");
        continue;
      }
      registrySkim.get<TH1>(HIST("skim/efficiencyMC"))->Fill(6., 1.);        // particles from tau in |eta|<0.9
      if (nChargedDaughtersTau[0] + nChargedDaughtersTau[1] == fourTracks) { // 4
        registrySkim.get<TH1>(HIST("skim/efficiencyMC"))->Fill(7., 1.);
      } else if (nChargedDaughtersTau[0] + nChargedDaughtersTau[1] == sixTracks) { // 6
        registrySkim.get<TH1>(HIST("skim/efficiencyMC"))->Fill(8., 1.);
      }

      if (nChargedDaughtersTau[0] == oneProng) // 1
        zerothTau = 0;
      else if (nChargedDaughtersTau[1] == oneProng) // 1
        zerothTau = 1;
      else if (nChargedDaughtersTau[0] == threeProng && nChargedDaughtersTau[1] == threeProng) // 3 and 3
        zerothTau = 0;

      // prepare local variables for output table
      int32_t runNumber = -999;
      int bc = -999;
      // int nTrks[3] = {-999, -999, -999}; // totalTracks, numContrib, globalNonPVtracks
      int totalTracks = -999;
      int8_t nPVcontrib = -99;
      int rct = -999;
      // float vtxPos[3] = {-999., -999., -999.};
      float zVertex = -999;

      int8_t recoMode = -99;
      int occupancy = -999;
      double hadronicRate = -999.;
      int8_t bcSels[8] = {-99, -99, -99, -99, -99, -99, -99, -99};
      // zdc information - there i sno information in MC
      float energyZNA = -999.;
      float energyZNC = -999.;

      float amplitudesFIT[3] = {-999., -999., -999.}; // FT0A, FT0C, FV0
      // float timesFIT[3] = {-999., -999., -999.};      // FT0A, FT0C, FV0
      // track momentum and sign
      float px[6] = {-999., -999., -999., -999., -999., -999.};
      float py[6] = {-999., -999., -999., -999., -999., -999.};
      float pz[6] = {-999., -999., -999., -999., -999., -999.};
      int sign[6] = {-99, -99, -99, -99, -99, -99};

      float dcaXY[6] = {-999., -999., -999., -999., -999., -999.};
      float dcaZ[6] = {-999., -999., -999., -999., -999., -999.};

      int nclTPCcrossedRows[6] = {-999, -999, -999, -999, -999, -999};
      int nclTPCfind[6] = {-999, -999, -999, -999, -999, -999};
      float nclTPCchi2[6] = {-999., -999., -999., -999., -999., -999.};
      float trkITSchi2[6] = {-999., -999., -999., -999., -999., -999.};
      int trkITScl[6] = {-999, -999, -999, -999, -999, -999};

      //      float trkTimeRes[2] = {-999., -999.};
      //      uint32_t itsClusterSizesTrk1 = 4294967295;
      //      uint32_t itsClusterSizesTrk2 = 4294967295;
      float tpcSignal[6] = {-999, -999, -999, -999, -999, -999};
      float tpcEl[6] = {-999, -999, -999, -999, -999, -999};
      float tpcMu[6] = {-999, -999, -999, -999, -999, -999};
      float tpcPi[6] = {-999, -999, -999, -999, -999, -999};
      float tpcKa[6] = {-999, -999, -999, -999, -999, -999};
      float tpcPr[6] = {-999, -999, -999, -999, -999, -999};
      // float tpcIP[2] = {-999, -999};
      float tofSignal[6] = {-999, -999, -999, -999, -999, -999};
      float tofEl[6] = {-999, -999, -999, -999, -999, -999};
      float tofMu[6] = {-999, -999, -999, -999, -999, -999};
      float tofPi[6] = {-999, -999, -999, -999, -999, -999};
      float tofKa[6] = {-999, -999, -999, -999, -999, -999};
      float tofPr[6] = {-999, -999, -999, -999, -999, -999};
      // float tofEP[2] = {-999, -999};
      float chi2TOF[6] = {-999., -999., -999., -999., -999., -999.};

      bool trueHasRecoColl = false;
      float trueDaugX[6] = {-999., -999., -999., -999., -999., -999.};
      float trueDaugY[6] = {-999., -999., -999., -999., -999., -999.};
      float trueDaugZ[6] = {-999., -999., -999., -999., -999., -999.};
      int trueDaugPdgCode[6] = {-999, -999, -999, -999, -999, -999};
      //      bool problem = false;
      MyRecoProblem problem = NO_PROBLEM;
      registrySkim.get<TH1>(HIST("skim/problemMC"))->Fill(NO_PROBLEM);

      // tau tau event type
      // 1 = e+3pi
      // 2 = mu+3pi
      // 3 = pi+3pi
      // 4 = 3pi+3pi

      if (nElec == oneProng && nPi == threeProng) // 1 + 3
        trueChannel = 1;
      else if (nMuon == oneProng && nPi == threeProng) // 1 + 3
        trueChannel = 2;
      else if (nPi == fourTracks) // 4
        trueChannel = 3;
      else if (nPi == sixTracks) // 6
        trueChannel = 4;

      // find reconstructed collisions associated to the generated collision
      auto const& collFromMcColls = collisions.sliceBy(colPerMcCollision, mccoll.globalIndex());
      if (verbose)
        LOGF(info, "coll from MC Coll %d", collFromMcColls.size());
      // check the generated collision was reconstructed
      if (collFromMcColls.size() > 0) { // get the truth and reco-level info
        if (verbose)
          LOGF(info, "MC Collision has reconstructed collision!");
        trueHasRecoColl = true;
        // check there is exactly one reco-level collision associated to generated collision
        if (collFromMcColls.size() > 1) {
          if (verbose)
            LOGF(info, "Truth collision has more than 1 reco collision. Skipping this event.");
          // histos.get<TH1>(HIST("Truth/hTroubles"))->Fill(1);
          //  problem = true;
          problem = MANY_RECO;
          registrySkim.get<TH1>(HIST("skim/problemMC"))->Fill(MANY_RECO);
          continue;
        }
        // grap reco-level collision
        auto const& collFromMcColl = collFromMcColls.iteratorAt(0);
        // grab tracks from the reco-level collision to get info to match measured data tables (processDataSG function)
        auto const& trksFromColl = tracks.sliceBy(trackPerCollision, collFromMcColl.globalIndex());
        // int countTracksPerCollision = 0;
        // int countGoodNonPVtracks = 0;
        // for (auto const& trkFromColl : trksFromColl) {
        //   // countTracksPerCollision++;
        //   if (!trkFromColl.isPVContributor()) {
        //     countGoodNonPVtracks++;
        //     continue;
        //   }
        // }

        // fill info for reconstructed collision
        runNumber = collFromMcColl.runNumber();
        bc = collFromMcColl.globalBC();
        totalTracks = trksFromColl.size();
        //        nTrks[0] = countTracksPerCollision;
        nPVcontrib = collFromMcColl.numContrib();
        //        nTrks[1] = collFromMcColl.numContrib();
        //        nTrks[2] = countGoodNonPVtracks;
        rct = isGoodRCTflag(collFromMcColl);
        zVertex = collFromMcColl.posZ();
        //        vtxPos[0] = collFromMcColl.posX();
        //        vtxPos[1] = collFromMcColl.posY();
        //        vtxPos[2] = collFromMcColl.posZ();

        recoMode = collFromMcColl.flags();
        occupancy = collFromMcColl.occupancyInTime();
        hadronicRate = collFromMcColl.hadronicRate();
        bcSels[0] = collFromMcColl.trs();
        bcSels[1] = collFromMcColl.trofs();
        bcSels[2] = collFromMcColl.hmpr();
        bcSels[3] = collFromMcColl.tfb();
        bcSels[4] = collFromMcColl.itsROFb();
        bcSels[5] = collFromMcColl.sbp();
        bcSels[6] = collFromMcColl.zVtxFT0vPV();
        bcSels[7] = collFromMcColl.vtxITSTPC();
        // energyZNA = collFromMcColl.energyCommonZNA();
        // energyZNC = collFromMcColl.energyCommonZNC();
        // if (energyZNA < 0)
        //   energyZNA = -1.;
        // if (energyZNC < 0)
        //   energyZNC = -1.;

        amplitudesFIT[0] = collFromMcColl.totalFT0AmplitudeA();
        amplitudesFIT[1] = collFromMcColl.totalFT0AmplitudeC();
        amplitudesFIT[2] = collFromMcColl.totalFV0AmplitudeA();
        // timesFIT[0] = collFromMcColl.timeFT0A();
        // timesFIT[1] = collFromMcColl.timeFT0C();
        // timesFIT[2] = collFromMcColl.timeFV0A();

        // get particles associated to generated collision
        auto const& partsFromMcColl = mcParticles.sliceBy(partPerMcCollision, mccoll.globalIndex());
        if (verbose)
          LOGF(info, "part from MC coll %d", partsFromMcColl.size());
        // int countMothers = 0;
        int countDaughters = 0;
        countPi0 = 0;
        for (const auto& particle : partsFromMcColl) {
          if (verbose)
            LOGF(info, "Reco coll; part pdg %d", particle.pdgCode());
          // select only tauons with checking if particle has no mother
          // in UPC MC taus have mothers
          // if (particle.has_mothers())
          if (std::abs(particle.pdgCode()) != kTauMinus)
            continue; // 15 = tau_minus

          // get daughters of the tau
          const auto& daughters = particle.daughters_as<aod::UDMcParticles>();
          // int countDaughters = 0;
          for (const auto& daughter : daughters) {
            if (verbose)
              LOGF(info, "With Coll; daug pdg %d", daughter.pdgCode());
            // check if it is the charged particle (= no pi0 or neutrino)
            if (enumMyParticle(daughter.pdgCode()) == MyOtherParticle) // -1
              continue;
            countDaughters++;
            if (daughter.pdgCode() == kPi0)
              countPi0++;

            // check whether 1+3 or 3+3 topology is present
            if (countDaughters > sixTracks) { // 6
              if (verbose)
                LOGF(info, "Truth collision has more than 6 charged daughters from 2 taus. Breaking the daughter loop.");
              //               histos.get<TH1>(HIST("Truth/hTroubles"))->Fill(3);
              //               problem = true;
              problem = TOO_MANY_DAUGHTERS;
              registrySkim.get<TH1>(HIST("skim/problemMC"))->Fill(TOO_MANY_DAUGHTERS);

              break;
            }

            // fill info for each daughter
            trueDaugX[countDaughters - 1] = daughter.px();
            trueDaugY[countDaughters - 1] = daughter.py();
            trueDaugZ[countDaughters - 1] = daughter.pz();
            trueDaugPdgCode[countDaughters - 1] = daughter.pdgCode();

            // get tracks associated to MC daughter (how well the daughter was reconstructed)
            auto const& tracksFromDaughter = tracks.sliceBy(trackPerMcParticle, daughter.globalIndex());
            // check there is exactly 1 track per 1 particle
            if (tracksFromDaughter.size() > 1) {
              if (verbose)
                LOGF(info, "Daughter has more than 1 associated track. Skipping this daughter.");
              //              histos.get<TH1>(HIST("Truth/hTroubles"))->Fill(4);
              //              problem = true;
              problem = TWO_TRACKS;
              registrySkim.get<TH1>(HIST("skim/problemMC"))->Fill(TWO_TRACKS);
              continue;
            }
            // grab the track and fill info for reconstructed track (should be done 4 or 6 times)
            const auto& trk = tracksFromDaughter.iteratorAt(0);
            if (verbose)
              LOGF(info, "p(%f,%f,%f)", trk.px(), trk.py(), trk.pz());
            px[countDaughters - 1] = trk.px();
            py[countDaughters - 1] = trk.py();
            pz[countDaughters - 1] = trk.pz();
            sign[countDaughters - 1] = trk.sign();
            dcaXY[countDaughters - 1] = trk.dcaXY();
            dcaZ[countDaughters - 1] = trk.dcaZ();
            //            trkTimeRes[countMothers - 1] = trk.trackTimeRes();
            //            if (countMothers == 1) {
            //              itsClusterSizesTrk1 = trk.itsClusterSizes();
            //            } else {
            //              itsClusterSizesTrk2 = trk.itsClusterSizes();
            //            }

            nclTPCcrossedRows[countDaughters - 1] = trk.tpcNClsCrossedRows();
            nclTPCfind[countDaughters - 1] = trk.tpcNClsFindable();
            nclTPCchi2[countDaughters - 1] = trk.tpcChi2NCl();
            trkITSchi2[countDaughters - 1] = trk.itsChi2NCl();
            trkITScl[countDaughters - 1] = numberOfItsClustersCheck(trk);

            tpcSignal[countDaughters - 1] = trk.tpcSignal();
            tpcEl[countDaughters - 1] = trk.tpcNSigmaEl();
            tpcMu[countDaughters - 1] = trk.tpcNSigmaMu();
            tpcPi[countDaughters - 1] = trk.tpcNSigmaPi();
            tpcKa[countDaughters - 1] = trk.tpcNSigmaKa();
            tpcPr[countDaughters - 1] = trk.tpcNSigmaPr();
            // tpcIP[countDaughters - 1] = trk.tpcInnerParam();

            tofSignal[countDaughters - 1] = trk.beta();
            tofEl[countDaughters - 1] = trk.tofNSigmaEl();
            tofMu[countDaughters - 1] = trk.tofNSigmaMu();
            tofPi[countDaughters - 1] = trk.tofNSigmaPi();
            tofKa[countDaughters - 1] = trk.tofNSigmaKa();
            tofPr[countDaughters - 1] = trk.tofNSigmaPr();
            // tofEP[countMothers - 1] = trk.tofExpMom();
            if (trk.hasTOF())
              chi2TOF[countDaughters - 1] = trk.tofChi2();

          } // daughters
        } // particles
      } else { // get only the truth information. The reco-level info is left on default
        if (verbose)
          LOGF(info, "MC Collision has NO reconstructed collision!");
        // get particles associated to generated collision
        auto const& partsFromMcColl = mcParticles.sliceBy(partPerMcCollision, mccoll.globalIndex());
        if (verbose)
          LOGF(info, "NO Coll; partsFromMcColl in MC %d", partsFromMcColl.size());
        // int countMothers = 0;
        int countDaughters = 0;
        countPi0 = 0;
        for (const auto& particle : partsFromMcColl) {
          if (verbose)
            LOGF(info, "No Coll; part Gid %d, Id %d, pdg %d, hasM %d, hasD %d", particle.globalIndex(), particle.index(), particle.pdgCode(), particle.has_mothers(), particle.has_daughters());
          // select only tauons with checking if particle has no mother
          // in UPC MC taus have mothers
          // if (particle.has_mothers())
          if (std::abs(particle.pdgCode()) != kTauMinus) // 15
            continue;
          // countMothers++;
          // check the generated collision does not have more than 2 tauons
          // if (countMothers > 2) {
          //   if (verbose)
          //     LOGF(info,"Truth collision has more than 2 no mother particles. Breaking the particle loop.");
          //        //     histos.get<TH1>(HIST("Truth/hTroubles"))->Fill(12);
          //        //     problem = true;
          //   break;
          // }
          // // fill info for each tau
          // trueTauX[countMothers - 1] = particle.px();
          // trueTauY[countMothers - 1] = particle.py();
          // trueTauZ[countMothers - 1] = particle.pz();

          // get daughters of the tau
          const auto& daughters = particle.daughters_as<aod::UDMcParticles>();
          if (verbose)
            LOGF(info, "NO coll; N_daughters %d", daughters.size());
          // int countDaughters = 0;
          for (const auto& daughter : daughters) {
            if (verbose)
              LOGF(info, "NO Coll; daug id %d, pdg %d", daughter.globalIndex(), daughter.pdgCode());
            // select only the charged particle (= no pi0 or neutrino)
            if (enumMyParticle(daughter.pdgCode()) == -1)
              continue;
            countDaughters++;
            if (daughter.pdgCode() == kPi0)
              countPi0++;

            // check whether 1+3 or 3+3 topology is present
            if (countDaughters > sixTracks) { // 6
              if (verbose)
                LOGF(info, "Truth collision has more than 6 charged daughters from taus. Breaking the daughter loop.");
              //              histos.get<TH1>(HIST("Truth/hTroubles"))->Fill(13);
              //              problem = true;
              registrySkim.get<TH1>(HIST("skim/problemMC"))->Fill(TOO_MANY_DAUGHTERS);
              problem = TOO_MANY_DAUGHTERS;
              break;
            }
            // fill info for each daughter
            trueDaugX[countDaughters - 1] = daughter.px();
            trueDaugY[countDaughters - 1] = daughter.py();
            trueDaugZ[countDaughters - 1] = daughter.pz();
            trueDaugPdgCode[countDaughters - 1] = daughter.pdgCode();
          } // daughters
          if (verbose)
            LOGF(info, "End of daughters");
        } // particles
      } // collisions

      // decide the channel and set the variable.
      trueChannel = trueChannel + countPi0 * 10 + zerothTau * 100;

      // LOGF(info, "Should be written!");
      trueTauFourTracks(runNumber,
                        bc, // is it necessary
                        totalTracks,
                        nPVcontrib,
                        rct,
                        // dgcand.posX(), dgcand.posY(),
                        zVertex,
                        recoMode,
                        occupancy,
                        hadronicRate,                    // is it necessary
                        bcSels[0], bcSels[1], bcSels[2], // to test it
                        bcSels[3], bcSels[4], bcSels[5], bcSels[6], bcSels[7],
                        energyZNA, energyZNC,
                        // qtot, <<-------- comment out
                        amplitudesFIT[0], amplitudesFIT[1], amplitudesFIT[2],
                        // timesFIT[0], timesFIT[1], timesFIT[2],
                        px, py, pz, sign,
                        dcaXY, dcaZ,
                        nclTPCcrossedRows, nclTPCfind, nclTPCchi2, trkITSchi2, trkITScl,
                        tpcSignal, tpcEl, tpcPi, tpcKa, tpcPr, tpcMu,
                        tofSignal, tofEl, tofPi, tofKa, tofPr, tofMu,
                        chi2TOF,
                        //
                        trueChannel,
                        trueHasRecoColl,
                        mccoll.posZ(),
                        trueTauX, trueTauY, trueTauZ,
                        trueDaugX, trueDaugY, trueDaugZ,
                        trueDaugPdgCode, problem);
    } // mccollisions
  } // end of  processMonteCarlo
  PROCESS_SWITCH(TauThreeProngEventTableProducer, processMonteCarlo, "Iterate UD tables with simulated data created by SG-Candidate-Producer.", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<TauThreeProngEventTableProducer>(cfgc)};
}
