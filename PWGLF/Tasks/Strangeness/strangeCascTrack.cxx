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
/// \file strangecasctrack.cxx
/// \brief Analysis of strangeness tracking efficiency via primary production of Omega and Xi in Run 3
/// \author Yakiv Paroviak (yakiv.paroviak@cern.ch)

#include "PWGLF/DataModel/LFStrangenessPIDTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/StaticFor.h"
#include "Framework/runDataProcessing.h"
#include <Framework/Configurable.h>

#include "TF1.h"
#include "TF2.h"
#include <Math/Vector4D.h>
#include <TPDGCode.h>

#include <string>
#include <vector>

using namespace o2;
using namespace o2::constants::math;
using namespace o2::framework;
using namespace o2::framework::expressions;

// tables for derived data
using DerCollisionWMults = soa::Join<aod::StraCollisions, aod::StraStamps, aod::StraEvSels, aod::StraCents>;
using DerCascDatas = soa::Join<aod::CascCores, aod::CascCollRefs, aod::CascExtras, aod::CascBBs, aod::CascTOFNSigmas>;
using DerTraCascDatas = soa::Join<aod::TraCascCores, aod::TraCascCollRefs, aod::StraTrackExtras, aod::TraToCascRefs>;

// tables for derived MC
using DerMCGenCascades = soa::Join<aod::CascMCCores, aod::CascMCCollRefs>;
using DerMCRecCollisions = soa::Join<aod::StraCollisions, aod::StraStamps, aod::StraCollLabels, aod::StraEvSels, aod::StraCents>;
using DerMCRecCascDatas = soa::Join<aod::CascCollRefs, aod::CascCores, aod::CascCoreMCLabels, aod::CascExtras, aod::CascBBs, aod::CascTOFNSigmas>;
using DerMCRecTraCascDatas = soa::Join<aod::TraCascCores, aod::TraCascCollRefs, aod::StraTrackExtras, aod::TraToCascRefs>;

// tables for PID selection
using DauTracks = soa::Join<aod::DauTrackExtras, aod::DauTrackTPCPIDs>;

struct StrangeCascTrack {

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Service<o2::framework::O2DatabasePDG> pdgDB;

  PresliceUnsorted<soa::Join<aod::StraCollisions, aod::StraCollLabels, aod::StraCents>> perMcCollision = aod::v0data::straMCCollisionId;

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // subprocess switches:
  Configurable<bool> doProcesspp{"doProcesspp", true, "true for pp"};
  Configurable<bool> doProcessPbPb{"doProcessPbPb", false, "true for PbPb"};
  Configurable<bool> doProcessOO{"doProcessOO", false, "true for OO"};
  Configurable<bool> doProcesspO{"doProcesspO", false, "true for pO"};

  Configurable<bool> doApplyEventCuts{"doApplyEventCuts", true, "apply general event cuts"}; // general cascade cuts - dca, cosPA etc.
  // Xi selections
  Configurable<bool> doApplyGenCutsXi{"doApplyGenCutsXi", true, "apply general cuts (Omega)"}; // general cascade cuts - dca, cosPA etc.
  Configurable<bool> doApplyPtCutsXi{"doApplyPtCutsXi", true, "apply pt cuts (Xi)"};           // ignore particles with extremely low efficiencies
  Configurable<bool> doApplyTPCPIDXi{"doApplyTPCPIDXi", true, "apply tpc pid to dau tracks (Xi)"};
  Configurable<bool> doApplyTOFPIDXi{"doApplyTOFPIDXi", true, "apply tof pid to dau tracks (Xi)"};
  // Omega selections
  Configurable<bool> doApplyGenCutsOmega{"doApplyGenCutsOmega", true, "apply general cuts (Omega)"}; // general cascade cuts - dca, cosPA etc.
  Configurable<bool> doApplyPtCutsOmega{"doApplyPtCutsOmega", true, "apply pt cuts (Omega)"};        // ignore particles with extremely low efficiencies
  Configurable<bool> doApplyTPCPIDOmega{"doApplyTPCPIDOmega", true, "apply tpc pid to dau tracks (Omega)"};
  Configurable<bool> doApplyTOFPIDOmega{"doApplyTOFPIDOmega", true, "apply tof pid to dau tracks (Omega)"};
  Configurable<bool> doCompetingMassRej{"doCompetingMassRej", true, "competing mass rejection (Omega)"};
  // efficiency and purity corrections:
  // only correct by pt
  Configurable<bool> doApplyEfficiency1D{"doApplyEfficiency1D", false, "apply efficiency correction"};
  Configurable<bool> doPropagateEfficiency1D{"doPropagateEfficiency1D", false, "apply efficiency propagation"};
  Configurable<bool> doApplyPurity1D{"doApplyPurity1D", false, "apply purity correction"};
  Configurable<bool> doPropagatePurity1D{"doPropagatePurity1D", false, "apply purity propagation"};
  // correct by both pt and mult
  Configurable<bool> doApplyEfficiency2D{"doApplyEfficiency2D", false, "apply efficiency correction"};
  Configurable<bool> doPropagateEfficiency2D{"doPropagateEfficiency2D", false, "apply efficiency propagation"};
  Configurable<bool> doApplyPurity2D{"doApplyPurity2D", false, "apply purity correction"};
  Configurable<bool> doPropagatePurity2D{"doPropagatePurity2D", false, "apply purity propagation"};
  Configurable<std::string> ccdburl{"ccdburl", "http://alice-ccdb.cern.ch", "url of the ccdb repository to use"};
  Configurable<std::string> efficiencyCCDBPathpp{"efficiencyCCDBPathpp", "Users/y/yparovia/LHC24f4d", "Path of the efficiency and purity corrections (pp)"};
  Configurable<std::string> efficiencyCCDBPathPbPb{"efficiencyCCDBPathPbPb", "Users/y/yparovia/LHC25f3", "Path of the efficiency and purity corrections (PbPb)"};
  Configurable<std::string> efficiencyCCDBPathOO{"efficiencyCCDBPathOO", "Users/y/yparovia/LHC25h3", "Path of the efficiency and purity corrections (OO)"};
  Configurable<std::string> efficiencyCCDBPathpO{"efficiencyCCDBPathpO", "Users/y/yparovia/LHC25h2", "Path of the efficiency and purity corrections (pO)"};

  // event and dau track selection
  struct : ConfigurableGroup {
    Configurable<double> cutZVertex{"cutZVertex", 10.0f, "max Z-vertex position"};
    Configurable<bool> cutSel8{"cutSel8", true, "choose events with sel8"};
    Configurable<double> cutDCAtoPVxy{"cutDCAtoPVxy", 0.02f, "max cascade dca to PV in xy"};
    Configurable<double> cutDCAtoPVz{"cutDCAtoPVz", 0.02f, "max cascade dca to PV in z"};
    Configurable<double> cutV0CosPA{"cutV0CosPA", 0.97f, "max V0 cosPA"};
    Configurable<double> cutBachCosPA{"cutBachCosPA", 0.97f, "max Bachelor cosPA"};
    Configurable<double> cutRapidity{"cutRapidity", 0.5f, "max rapidity"};
    Configurable<double> cutCompMassRej{"cutCompMassRej", 0.008f, "Competing mass rejection"};
    // minimum and maximum desired pt
    Configurable<double> cutMinPtXiStd{"cutMinPtXiStd", 0.0f, "min pt for standard Xi"};
    Configurable<double> cutMaxPtXiStd{"cutMaxPtXiStd", 15.0f, "min pt for standard Xi"};
    Configurable<double> cutMinPtXiTra{"cutMinPtXiTra", 0.5f, "min pt for tracked Xi"};
    Configurable<double> cutMaxPtXiTra{"cutMaxPtXiTra", 15.0f, "min pt for standard Xi"};
    Configurable<double> cutMinPtOmegaStd{"cutMinPtOmegaStd", 0.5f, "min pt for standard Omega"};
    Configurable<double> cutMaxPtOmegaStd{"cutMaxPtOmegaStd", 15.0f, "min pt for standard Omega"};
    Configurable<double> cutMinPtOmegaTra{"cutMinPtOmegaTra", 1.0f, "min pt for tracked Omega"};
    Configurable<double> cutMaxPtOmegaTra{"cutMaxPtOmegaTra", 15.0f, "min pt for standard Omega"};
    // TPC PID selection
    Configurable<float> nSigmaTPCPion{"nSigmaTPCPion", 4, "NSigmaTPCPion"};
    Configurable<float> nSigmaTPCKaon{"nSigmaTPCKaon", 4, "NSigmaTPCKaon"};
    Configurable<float> nSigmaTPCProton{"nSigmaTPCProton", 4, "NSigmaTPCProton"};
    // TOF PID selection
    Configurable<float> nSigmaTOFXi{"nSigmaTOFXi", 3, "nSigmaTOFXi"};
    Configurable<float> nSigmaTOFOmega{"nSigmaTOFOmega", 3, "nSigmaTOFOmega"};
  } selCuts;

  // axes
  struct : ConfigurableGroup {
    ConfigurableAxis axisPhi{"axisPhi", {72, 0, TwoPI}, "#phi"};
    ConfigurableAxis axisEta{"axisEta", {102, -2.01, 2.01}, "#eta"};
    ConfigurableAxis axisDCAxy{"axisDCAxy", {500, 0., 0.5}, "cm"};
    ConfigurableAxis axisDCAz{"axisDCAz", {500, 0., 0.5}, "cm"};
    ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 10.0}, "p_{T} (GeV/c)"};
    ConfigurableAxis axisMult{"axisMult", {VARIABLE_WIDTH, 0.0f, 5.0, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 70.0f, 100.0f}, "FT0 mult %"};
    ConfigurableAxis axisOmegaMass{"axisOmegaMass", {2000, 1.6, 1.8}, "#Omega M_{inv} (GeV/c^{2})"};
    ConfigurableAxis axisXiMass{"axisXiMass", {2000, 1.2, 1.4}, "#Xi M_{inv} (GeV/c^{2})"};
  } axesConfig;

  // cascade reconstruction types
  static constexpr std::string_view kTypeNames[] = {"Standard", "Tracked"};

  // for efficiency and purity corrections
  TH1F* hEfficiencyOmegaStd1D;
  TH1F* hEfficiencyOmegaTra1D;
  TH1F* hEfficiencyXiStd1D;
  TH1F* hEfficiencyXiTra1D;

  TH1F* hEfficiencyErrOmegaStd1D;
  TH1F* hEfficiencyErrOmegaTra1D;
  TH1F* hEfficiencyErrXiStd1D;
  TH1F* hEfficiencyErrXiTra1D;

  TH1F* hPurityOmegaStd1D;
  TH1F* hPurityOmegaTra1D;
  TH1F* hPurityXiStd1D;
  TH1F* hPurityXiTra1D;

  TH1F* hPurityErrOmegaStd1D;
  TH1F* hPurityErrOmegaTra1D;
  TH1F* hPurityErrXiStd1D;
  TH1F* hPurityErrXiTra1D;

  TH2F* hEfficiencyOmegaStd2D;
  TH2F* hEfficiencyOmegaTra2D;
  TH2F* hEfficiencyXiStd2D;
  TH2F* hEfficiencyXiTra2D;

  TH2F* hEfficiencyErrOmegaStd2D;
  TH2F* hEfficiencyErrOmegaTra2D;
  TH2F* hEfficiencyErrXiStd2D;
  TH2F* hEfficiencyErrXiTra2D;

  TH2F* hPurityOmegaStd2D;
  TH2F* hPurityOmegaTra2D;
  TH2F* hPurityXiStd2D;
  TH2F* hPurityXiTra2D;

  TH2F* hPurityErrOmegaStd2D;
  TH2F* hPurityErrOmegaTra2D;
  TH2F* hPurityErrXiStd2D;
  TH2F* hPurityErrXiTra2D;

  int mRunNumber;
  // loads efficiencies and purities
  void initEfficiencyFromCCDB(int64_t runNumber, int64_t timestamp)
  {
    if (mRunNumber == runNumber) {
      return;
    }
    mRunNumber = runNumber;
    LOG(info) << "Loading efficiencies and purities from CCDB for run " << mRunNumber << " now...";
    auto timeStamp = timestamp;

    std::string efficiencyCCDBPath = [&]() {
      if (doProcesspp) {
        return efficiencyCCDBPathpp;
      } else if (doProcesspO) {
        return efficiencyCCDBPathpO;
      } else if (doProcessPbPb) {
        return efficiencyCCDBPathPbPb;
      }
      return efficiencyCCDBPathOO;
    }();

    TList* listEfficiencies = ccdb->getForTimeStamp<TList>(efficiencyCCDBPath, timeStamp);

    if (!listEfficiencies) {
      LOG(fatal) << "Problem getting TList object with efficiencies and purities!";
    }

    hEfficiencyOmegaStd1D = static_cast<TH1F*>(listEfficiencies->FindObject("Eff_Omega_Standard_byPt"));
    hEfficiencyOmegaTra1D = static_cast<TH1F*>(listEfficiencies->FindObject("Eff_Omega_Tracked_byPt"));
    hEfficiencyXiStd1D = static_cast<TH1F*>(listEfficiencies->FindObject("Eff_Xi_Standard_byPt"));
    hEfficiencyXiTra1D = static_cast<TH1F*>(listEfficiencies->FindObject("Eff_Xi_Tracked_byPt"));
    hEfficiencyErrOmegaStd1D = static_cast<TH1F*>(listEfficiencies->FindObject("EffErr_Omega_Standard_byPt"));
    hEfficiencyErrOmegaTra1D = static_cast<TH1F*>(listEfficiencies->FindObject("EffErr_Omega_Tracked_byPt"));
    hEfficiencyErrXiStd1D = static_cast<TH1F*>(listEfficiencies->FindObject("EffErr_Xi_Standard_byPt"));
    hEfficiencyErrXiTra1D = static_cast<TH1F*>(listEfficiencies->FindObject("EffErr_Xi_Tracked_byPt"));
    hPurityOmegaStd1D = static_cast<TH1F*>(listEfficiencies->FindObject("Pur_Omega_Standard_byPt"));
    hPurityOmegaTra1D = static_cast<TH1F*>(listEfficiencies->FindObject("Pur_Omega_Tracked_byPt"));
    hPurityXiStd1D = static_cast<TH1F*>(listEfficiencies->FindObject("Pur_Xi_Standard_byPt"));
    hPurityXiTra1D = static_cast<TH1F*>(listEfficiencies->FindObject("Pur_Xi_Tracked_byPt"));
    hPurityErrOmegaStd1D = static_cast<TH1F*>(listEfficiencies->FindObject("PurErr_Omega_Standard_byPt"));
    hPurityErrOmegaTra1D = static_cast<TH1F*>(listEfficiencies->FindObject("PurErr_Omega_Tracked_byPt"));
    hPurityErrXiStd1D = static_cast<TH1F*>(listEfficiencies->FindObject("PurErr_Xi_Standard_byPt"));
    hPurityErrXiTra1D = static_cast<TH1F*>(listEfficiencies->FindObject("PurErr_Xi_Tracked_byPt"));

    hEfficiencyOmegaStd2D = static_cast<TH2F*>(listEfficiencies->FindObject("Eff_Omega_Standard_byPtMult"));
    hEfficiencyOmegaTra2D = static_cast<TH2F*>(listEfficiencies->FindObject("Eff_Omega_Tracked_byPtMult"));
    hEfficiencyXiStd2D = static_cast<TH2F*>(listEfficiencies->FindObject("Eff_Xi_Standard_byPtMult"));
    hEfficiencyXiTra2D = static_cast<TH2F*>(listEfficiencies->FindObject("Eff_Xi_Tracked_byPtMult"));
    hEfficiencyErrOmegaStd2D = static_cast<TH2F*>(listEfficiencies->FindObject("EffErr_Omega_Standard_byPtMult"));
    hEfficiencyErrOmegaTra2D = static_cast<TH2F*>(listEfficiencies->FindObject("EffErr_Omega_Tracked_byPtMult"));
    hEfficiencyErrXiStd2D = static_cast<TH2F*>(listEfficiencies->FindObject("EffErr_Xi_Standard_byPtMult"));
    hEfficiencyErrXiTra2D = static_cast<TH2F*>(listEfficiencies->FindObject("EffErr_Xi_Tracked_byPtMult"));
    hPurityOmegaStd2D = static_cast<TH2F*>(listEfficiencies->FindObject("Pur_Omega_Standard_byPtMult"));
    hPurityOmegaTra2D = static_cast<TH2F*>(listEfficiencies->FindObject("Pur_Omega_Tracked_byPtMult"));
    hPurityXiStd2D = static_cast<TH2F*>(listEfficiencies->FindObject("Pur_Xi_Standard_byPtMult"));
    hPurityXiTra2D = static_cast<TH2F*>(listEfficiencies->FindObject("Pur_Xi_Tracked_byPtMult"));
    hPurityErrOmegaStd2D = static_cast<TH2F*>(listEfficiencies->FindObject("PurErr_Omega_Standard_byPtMult"));
    hPurityErrOmegaTra2D = static_cast<TH2F*>(listEfficiencies->FindObject("PurErr_Omega_Tracked_byPtMult"));
    hPurityErrXiStd2D = static_cast<TH2F*>(listEfficiencies->FindObject("PurErr_Xi_Standard_byPtMult"));
    hPurityErrXiTra2D = static_cast<TH2F*>(listEfficiencies->FindObject("PurErr_Xi_Tracked_byPtMult"));

    if (doPropagateEfficiency1D && (!hEfficiencyErrOmegaStd1D || !hEfficiencyErrOmegaTra1D || !hEfficiencyErrXiStd1D || !hEfficiencyErrXiTra1D))
      LOG(fatal) << "Problem getting hEfficiencyUncertainty!";
    if (doPropagatePurity1D && (!hPurityErrOmegaStd1D || !hPurityErrOmegaTra1D || !hPurityErrXiStd1D || !hPurityErrXiTra1D))
      LOG(fatal) << "Problem getting hPurityUncertainty!";
    LOG(info) << "Efficiencies and purities now loaded for " << mRunNumber;

    if (doPropagateEfficiency2D && (!hEfficiencyErrOmegaStd2D || !hEfficiencyErrOmegaTra2D || !hEfficiencyErrXiStd2D || !hEfficiencyErrXiTra2D))
      LOG(fatal) << "Problem getting hEfficiencyUncertainty!";
    if (doPropagatePurity2D && (!hPurityErrOmegaStd2D || !hPurityErrOmegaTra2D || !hPurityErrXiStd2D || !hPurityErrXiTra2D))
      LOG(fatal) << "Problem getting hPurityUncertainty!";
    LOG(info) << "Efficiencies and purities now loaded for " << mRunNumber;
  }
  // general info about processed events
  template <typename TEvent>
  void fillEvents(TEvent const& collision)
  {
    histos.fill(HIST("Events/EvCounter"), 0.5);
    double mult = (doProcesspp || doProcesspO) ? collision.centFT0M() : collision.centFT0C();
    histos.fill(HIST("Events/Mult"), mult);
    double pvx = collision.posX();
    double pvy = collision.posY();
    double pvz = collision.posZ();
    histos.fill(HIST("Events/PVx"), pvx);
    histos.fill(HIST("Events/PVy"), pvy);
    histos.fill(HIST("Events/PVz"), pvz);
  }
  // checks general selection criteria for collisions
  template <typename TEvent>
  bool isValidEvent(TEvent collision)
  {
    if (std::abs(collision.posZ()) > selCuts.cutZVertex)
      return false;
    if (selCuts.cutSel8 && !collision.sel8())
      return false;
    return true;
  }
  // checks general selection criteria for cascades
  template <typename TEvent, typename TCascade>
  bool isValidCasc(TEvent collision, TCascade cascade, TString particle)
  {
    if (cascade.dcaXYCascToPV() > selCuts.cutDCAtoPVxy)
      return false;
    if (cascade.dcaZCascToPV() > selCuts.cutDCAtoPVz)
      return false;
    if (cascade.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) < selCuts.cutV0CosPA)
      return false;
    if (cascade.bachBaryonCosPA() < selCuts.cutBachCosPA)
      return false;
    ROOT::Math::PxPyPzMVector momentum;
    if (particle == "xi")
      momentum.SetCoordinates(cascade.px(), cascade.py(), cascade.pz(), o2::constants::physics::MassXiMinus);
    else
      momentum.SetCoordinates(cascade.px(), cascade.py(), cascade.pz(), o2::constants::physics::MassOmegaMinus);
    if (std::abs(momentum.Rapidity()) > selCuts.cutRapidity)
      return false;
    return true;
  }
  // checks TPC PID of dau tracks
  template <typename TCascade>
  bool passesTPC(TCascade cascade)
  {
    const auto& posDaughterTrackCasc = cascade.template posTrackExtra_as<DauTracks>();
    const auto& negDaughterTrackCasc = cascade.template negTrackExtra_as<DauTracks>();
    if (cascade.sign() < 0) {
      if (std::abs(posDaughterTrackCasc.tpcNSigmaPr()) > selCuts.nSigmaTPCProton) {
        return false;
      }
      if (std::abs(negDaughterTrackCasc.tpcNSigmaPi()) > selCuts.nSigmaTPCPion) {
        return false;
      }
    } else {
      if (std::abs(negDaughterTrackCasc.tpcNSigmaPr()) > selCuts.nSigmaTPCProton) {
        return false;
      }
      if (std::abs(posDaughterTrackCasc.tpcNSigmaPi()) > selCuts.nSigmaTPCPion) {
        return false;
      }
    }
    return true;
  }
  // checks TOF PID of dau tracks
  template <typename TCascade>
  bool passesTOF(TCascade cascade, TString particle)
  {
    if (particle == "xi")
      return cascade.tofXiCompatibility(selCuts.nSigmaTOFXi);
    if (particle == "omega")
      return cascade.tofOmegaCompatibility(selCuts.nSigmaTOFOmega);
    return true;
  }
  // checks whether gen cascade corresponds to PDG code
  template <typename TCascade>
  bool isValidPDG(TCascade cascade, TString particle)
  {
    if (particle == "xi" && std::abs(cascade.pdgCode()) == PDG_t::kXiMinus)
      return true;
    if (particle == "omega" && std::abs(cascade.pdgCode()) == PDG_t::kOmegaMinus)
      return true;
    return false;
  }
  // checks whether rec cascade is a truth primary xi or omega
  template <typename TCascade>
  bool isMCTruth(const TCascade& cascade, TString particle)
  {
    if constexpr (requires { cascade.has_cascMCCore(); }) { // safety check: discard rec cascade without gen reference
      auto cascmccore = cascade.template cascMCCore_as<DerMCGenCascades>();
      if (!cascmccore.isPhysicalPrimary())
        return false;
      int pdg = std::abs(cascmccore.pdgCode());
      if (particle == "xi")
        return (pdg == PDG_t::kXiMinus);
      if (particle == "omega")
        return (pdg == PDG_t::kOmegaMinus);
    }
    return false;
  }
  // applies purities and efficiencies
  void fillHist(std::shared_ptr<THn> hist, double binFillThn[], float efficiency, float effUncert, float purity, float purityUncert)
  {
    float previousContent, previousError2, currentContent, currentError2;
    int bin = hist->GetBin(binFillThn);
    previousContent = hist->GetBinContent(bin);
    previousError2 = hist->GetBinError2(bin);
    currentContent = previousContent + purity / (efficiency);
    currentError2 = previousError2 + std::pow(purity / (efficiency), 2) + std::pow(purityUncert / (efficiency), 2) + std::pow(effUncert * purity, 2) / std::pow(efficiency, 4);
    hist->SetBinContent(bin, currentContent);
    hist->SetBinError2(bin, currentError2);
  }
  // applies selections for and fills histograms
  template <typename TEvent, typename TCascs>
  void analyseCascs(TEvent collision, TCascs cascades)
  {
    int64_t casccollid = 0;
    for (auto const& cascade : cascades) {

      if constexpr (requires { cascade.topologyChi2(); }) {
        if (!cascade.has_standardCascade())
          continue; // safety check: dismisses tracked cascades without proper reference
      }

      // for tracked cascades, make a reference to standard table
      auto stdCasc = [&]() {
        if constexpr (requires { cascade.topologyChi2(); }) {
          if constexpr (requires { collision.straMCCollisionId(); }) {
            return cascade.template standardCascade_as<DerMCRecCascDatas>();
          } else {
            return cascade.template standardCascade_as<DerCascDatas>();
          }
        } else {
          return cascade;
        }
      }();

      // fill cascade statistics without any selections
      static constexpr int type = [&]() {
        if constexpr (requires { cascade.topologyChi2(); }) {
          return 1;
        } else {
          return 0;
        }
      }();

      double mult = (doProcesspp || doProcesspO) ? collision.centFT0M() : collision.centFT0C(); // ion collisions use FT0C for multiplicity, pp uses both

      float efficiencyOmega = 1.0f;
      float efficiencyXi = 1.0f;
      float efficiencyOmegaErr = 0.0f;
      float efficiencyXiErr = 0.0f;
      float purityOmega = 1.0f;
      float purityXi = 1.0f;
      float purityOmegaErr = 0.0f;
      float purityXiErr = 0.0f;

      if (doApplyEfficiency1D) {
        if constexpr (requires { cascade.topologyChi2(); }) {
          efficiencyOmega = hEfficiencyOmegaTra1D->Interpolate(cascade.pt());
          efficiencyXi = hEfficiencyXiTra1D->Interpolate(cascade.pt());
          if (doPropagateEfficiency1D) {
            efficiencyOmegaErr = hEfficiencyErrOmegaTra1D->Interpolate(cascade.pt());
            efficiencyXiErr = hEfficiencyErrXiTra1D->Interpolate(cascade.pt());
          }
          if (efficiencyOmega == 0) { // check for zero efficiency, do not apply if the case
            efficiencyOmega = 1.;
            efficiencyOmegaErr = 0.;
          }
          if (efficiencyXi == 0) { // check for zero efficiency, do not apply if the case
            efficiencyXi = 1.;
            efficiencyXiErr = 0.;
          }
        } else {
          efficiencyOmega = hEfficiencyOmegaStd1D->Interpolate(cascade.pt());
          efficiencyXi = hEfficiencyXiStd1D->Interpolate(cascade.pt());
          if (doPropagateEfficiency1D) {
            efficiencyOmegaErr = hEfficiencyErrOmegaStd1D->Interpolate(cascade.pt());
            efficiencyXiErr = hEfficiencyErrXiStd1D->Interpolate(cascade.pt());
          }
          if (efficiencyOmega == 0) { // check for zero efficiency, do not apply if the case
            efficiencyOmega = 1.;
            efficiencyOmegaErr = 0.;
          }
          if (efficiencyXi == 0) { // check for zero efficiency, do not apply if the case
            efficiencyXi = 1.;
            efficiencyXiErr = 0.;
          }
        }
      }

      if (doApplyEfficiency2D) {
        if constexpr (requires { cascade.topologyChi2(); }) {
          efficiencyOmega = hEfficiencyOmegaTra2D->Interpolate(cascade.pt(), mult);
          efficiencyXi = hEfficiencyXiTra2D->Interpolate(cascade.pt(), mult);
          if (doPropagateEfficiency2D) {
            efficiencyOmegaErr = hEfficiencyErrOmegaTra2D->Interpolate(cascade.pt(), mult);
            efficiencyXiErr = hEfficiencyErrXiTra2D->Interpolate(cascade.pt(), mult);
          }
          if (efficiencyOmega == 0) { // check for zero efficiency, do not apply if the case
            efficiencyOmega = 1.;
            efficiencyOmegaErr = 0.;
          }
          if (efficiencyXi == 0) { // check for zero efficiency, do not apply if the case
            efficiencyXi = 1.;
            efficiencyXiErr = 0.;
          }
        } else {
          efficiencyOmega = hEfficiencyOmegaStd2D->Interpolate(cascade.pt(), mult);
          efficiencyXi = hEfficiencyXiStd2D->Interpolate(cascade.pt(), mult);
          if (doPropagateEfficiency2D) {
            efficiencyOmegaErr = hEfficiencyErrOmegaStd2D->Interpolate(cascade.pt(), mult);
            efficiencyXiErr = hEfficiencyErrXiStd2D->Interpolate(cascade.pt(), mult);
          }
          if (efficiencyOmega == 0) { // check for zero efficiency, do not apply if the case
            efficiencyOmega = 1.;
            efficiencyOmegaErr = 0.;
          }
          if (efficiencyXi == 0) { // check for zero efficiency, do not apply if the case
            efficiencyXi = 1.;
            efficiencyXiErr = 0.;
          }
        }
      }

      if (doApplyPurity1D) {
        if constexpr (requires { cascade.topologyChi2(); }) {
          purityOmega = hPurityOmegaTra1D->Interpolate(cascade.pt(), mult);
          purityXi = hPurityXiTra1D->Interpolate(cascade.pt(), mult);
          if (doPropagatePurity1D) {
            purityOmegaErr = hPurityErrOmegaTra1D->Interpolate(cascade.pt(), mult);
            purityXiErr = hPurityErrXiTra1D->Interpolate(cascade.pt(), mult);
          }
          if (purityOmega == 0) { // check for zero purity, do not apply if the case
            purityOmega = 1.;
            purityOmegaErr = 0.;
          }
          if (purityXi == 0) { // check for zero purity, do not apply if the case
            purityXi = 1.;
            purityXiErr = 0.;
          }
        } else {
          purityOmega = hPurityOmegaStd1D->Interpolate(cascade.pt(), mult);
          purityXi = hPurityXiStd1D->Interpolate(cascade.pt(), mult);
          if (doPropagatePurity1D) {
            purityOmegaErr = hPurityErrOmegaStd1D->Interpolate(cascade.pt(), mult);
            purityXiErr = hPurityErrXiStd1D->Interpolate(cascade.pt(), mult);
          }
          if (purityOmega == 0) { // check for zero purity, do not apply if the case
            purityOmega = 1.;
            purityOmegaErr = 0.;
          }
          if (purityXi == 0) { // check for zero purity, do not apply if the case
            purityXi = 1.;
            purityXiErr = 0.;
          }
        }
      }

      if (doApplyPurity2D) {
        if constexpr (requires { cascade.topologyChi2(); }) {
          purityOmega = hPurityOmegaTra2D->Interpolate(cascade.pt(), mult);
          purityXi = hPurityXiTra2D->Interpolate(cascade.pt(), mult);
          if (doPropagatePurity2D) {
            purityOmegaErr = hPurityErrOmegaTra2D->Interpolate(cascade.pt(), mult);
            purityXiErr = hPurityErrXiTra2D->Interpolate(cascade.pt(), mult);
          }
          if (purityOmega == 0) { // check for zero purity, do not apply if the case
            purityOmega = 1.;
            purityOmegaErr = 0.;
          }
          if (purityXi == 0) { // check for zero purity, do not apply if the case
            purityXi = 1.;
            purityXiErr = 0.;
          }
        } else {
          purityOmega = hPurityOmegaStd2D->Interpolate(cascade.pt(), mult);
          purityXi = hPurityXiStd2D->Interpolate(cascade.pt(), mult);
          if (doPropagatePurity2D) {
            purityOmegaErr = hPurityErrOmegaStd2D->Interpolate(cascade.pt(), mult);
            purityXiErr = hPurityErrXiStd2D->Interpolate(cascade.pt(), mult);
          }
          if (purityOmega == 0) { // check for zero purity, do not apply if the case
            purityOmega = 1.;
            purityOmegaErr = 0.;
          }
          if (purityXi == 0) { // check for zero purity, do not apply if the case
            purityXi = 1.;
            purityXiErr = 0.;
          }
        }
      }

      if (collision.index() != casccollid) {
        histos.fill(HIST(kTypeNames[type]) + HIST("/NoSel/EvMult"), mult);
        casccollid = collision.index();
      }

      double massXi = cascade.mXi();
      double massOmega = cascade.mOmega();
      double pt = cascade.pt();
      double v0cosPA = stdCasc.v0cosPA(collision.posX(), collision.posY(), collision.posZ());

      histos.fill(HIST(kTypeNames[type]) + HIST("/NoSel/DCAxy"), cascade.dcaXYCascToPV());
      histos.fill(HIST(kTypeNames[type]) + HIST("/NoSel/DCAz"), cascade.dcaZCascToPV());
      histos.fill(HIST(kTypeNames[type]) + HIST("/NoSel/Phi"), cascade.phi());
      histos.fill(HIST(kTypeNames[type]) + HIST("/NoSel/Eta"), cascade.eta());
      histos.fill(HIST(kTypeNames[type]) + HIST("/NoSel/BachCosPA"), stdCasc.bachBaryonCosPA());
      histos.fill(HIST(kTypeNames[type]) + HIST("/NoSel/V0CosPA"), v0cosPA);
      histos.fill(HIST(kTypeNames[type]) + HIST("/NoSel/MassXi"), massXi);
      histos.fill(HIST(kTypeNames[type]) + HIST("/NoSel/MassOmega"), massOmega);

      if constexpr (requires { collision.straMCCollisionId(); }) {
        if (isMCTruth(stdCasc, "xi") || isMCTruth(stdCasc, "omega")) {
          histos.fill(HIST(kTypeNames[type]) + HIST("/NoSel-Truth/DCAxy"), cascade.dcaXYCascToPV());
          histos.fill(HIST(kTypeNames[type]) + HIST("/NoSel-Truth/DCAz"), cascade.dcaZCascToPV());
          histos.fill(HIST(kTypeNames[type]) + HIST("/NoSel-Truth/DCAzVSpt"), pt, cascade.dcaZCascToPV());
          histos.fill(HIST(kTypeNames[type]) + HIST("/NoSel-Truth/Phi"), cascade.phi());
          histos.fill(HIST(kTypeNames[type]) + HIST("/NoSel-Truth/Eta"), cascade.eta());
          histos.fill(HIST(kTypeNames[type]) + HIST("/NoSel-Truth/EvMult"), mult);
          histos.fill(HIST(kTypeNames[type]) + HIST("/NoSel-Truth/BachCosPA"), stdCasc.bachBaryonCosPA());
          histos.fill(HIST(kTypeNames[type]) + HIST("/NoSel-Truth/V0CosPA"), v0cosPA);
          if (isMCTruth(stdCasc, "xi"))
            histos.fill(HIST(kTypeNames[type]) + HIST("/NoSel-Truth/MassXi"), massXi);
          if (isMCTruth(stdCasc, "omega"))
            histos.fill(HIST(kTypeNames[type]) + HIST("/NoSel-Truth/MassOmega"), massOmega);
        }
      }

      // start checking selections
      bool passedAllSelsXi = true;
      bool passedAllSelsOmega = true;
      // apply general selection criteria
      if (doApplyEventCuts) {
        if (!isValidEvent(collision)) {
          passedAllSelsXi = false;
          passedAllSelsOmega = false;
        }
      }
      if (doApplyGenCutsXi) {
        if (!isValidCasc(collision, stdCasc, "xi"))
          passedAllSelsXi = false;
      }
      if (doApplyGenCutsOmega) {
        if (!isValidCasc(collision, stdCasc, "omega"))
          passedAllSelsOmega = false;
      }
      // apply pt cuts
      if constexpr (requires { cascade.topologyChi2(); }) {
        if (doApplyPtCutsXi) {
          if (pt < selCuts.cutMinPtXiTra || pt > selCuts.cutMaxPtXiTra)
            passedAllSelsXi = false;
        }
        if (doApplyPtCutsOmega) {
          if (pt < selCuts.cutMinPtOmegaTra || pt > selCuts.cutMaxPtOmegaTra)
            passedAllSelsOmega = false;
        }
      } else {
        if (doApplyPtCutsXi) {
          if (pt < selCuts.cutMinPtXiStd || pt > selCuts.cutMaxPtXiStd)
            passedAllSelsXi = false;
        }
        if (doApplyPtCutsOmega) {
          if (pt < selCuts.cutMinPtOmegaStd || pt > selCuts.cutMaxPtOmegaStd)
            passedAllSelsOmega = false;
        }
      }
      // apply tpc pid
      if (doApplyTPCPIDXi) {
        if (!passesTPC(stdCasc))
          passedAllSelsXi = false;
      }
      if (doApplyTPCPIDOmega) {
        if (!passesTPC(stdCasc))
          passedAllSelsOmega = false;
      }
      // apply tof pid
      if (doApplyTOFPIDXi) {
        if (!passesTOF(stdCasc, "xi"))
          passedAllSelsXi = false;
      }
      if (doApplyTOFPIDOmega) {
        if (!passesTOF(stdCasc, "omega"))
          passedAllSelsOmega = false;
      }
      // apply competing mass rej
      if (doCompetingMassRej) {
        if (!(std::abs(massXi - o2::constants::physics::MassXiMinus) > selCuts.cutCompMassRej))
          passedAllSelsOmega = false;
      }

      // fill truth w/ cascs that passed all applied sels
      double binFillXi[3] = {massXi, pt, mult};

      if constexpr (requires { collision.straMCCollisionId(); }) {
        if (passedAllSelsXi || passedAllSelsOmega) { // fill once for every desired cascade
          if (isMCTruth(stdCasc, "xi") || isMCTruth(stdCasc, "omega")) {
            histos.fill(HIST(kTypeNames[type]) + HIST("/Rec-Truth/DCAxy"), cascade.dcaXYCascToPV());
            histos.fill(HIST(kTypeNames[type]) + HIST("/Rec-Truth/DCAz"), cascade.dcaZCascToPV());
            histos.fill(HIST(kTypeNames[type]) + HIST("/Rec-Truth/DCAzVSpt"), pt, cascade.dcaZCascToPV());
            histos.fill(HIST(kTypeNames[type]) + HIST("/Rec-Truth/Phi"), cascade.phi());
            histos.fill(HIST(kTypeNames[type]) + HIST("/Rec-Truth/Eta"), cascade.eta());
            histos.fill(HIST(kTypeNames[type]) + HIST("/Rec-Truth/EvMult"), mult);
            histos.fill(HIST(kTypeNames[type]) + HIST("/Rec-Truth/BachCosPA"), stdCasc.bachBaryonCosPA());
            histos.fill(HIST(kTypeNames[type]) + HIST("/Rec-Truth/V0CosPA"), v0cosPA);
          }
        }
      }

      // fill rec
      if (passedAllSelsXi) {
        histos.fill(HIST(kTypeNames[type]) + HIST("/Rec/MassXi"), massXi);
        fillHist(histos.get<THn>(HIST(kTypeNames[type]) + HIST("/Rec/Xi")), binFillXi, efficiencyXi, efficiencyXiErr, purityXi, purityXiErr);
        if constexpr (requires { collision.straMCCollisionId(); }) {
          if (isMCTruth(stdCasc, "xi")) {
            histos.fill(HIST(kTypeNames[type]) + HIST("/Rec-Truth/MassXi"), massXi);
            histos.fill(HIST(kTypeNames[type]) + HIST("/Rec-Truth/Xi"), massXi, pt, mult);
          }
        }
      }
      double binFillOmega[3] = {massOmega, pt, mult};
      if (passedAllSelsOmega) {
        histos.fill(HIST(kTypeNames[type]) + HIST("/Rec/MassOmega"), massOmega);
        fillHist(histos.get<THn>(HIST(kTypeNames[type]) + HIST("/Rec/Omega")), binFillOmega, efficiencyOmega, efficiencyOmegaErr, purityOmega, purityOmegaErr);
        if constexpr (requires { collision.straMCCollisionId(); }) {
          if (isMCTruth(stdCasc, "omega")) {
            histos.fill(HIST(kTypeNames[type]) + HIST("/Rec-Truth/MassOmega"), massOmega);
            histos.fill(HIST(kTypeNames[type]) + HIST("/Rec-Truth/Omega"), massOmega, pt, mult);
          }
        }
      }
    }
  }

  void init(InitContext const&)
  {
    // for all events processing
    histos.add("Events/EvCounter", "Event Counter", kTH1F, {{1, 0, 1}});
    histos.add("Events/PVx", "PV x position", kTH1F, {{200, -0.1, 0.1}});
    histos.add("Events/PVy", "PV y position", kTH1F, {{200, -0.1, 0.1}});
    histos.add("Events/PVz", "PV z position", kTH1F, {{100, -20, 20}});
    histos.add("Events/Mult", "Multiplicity", kTH1F, {axesConfig.axisMult});
    // for cascade processing
    static_for<0, 1>([&](auto type) {
      // no selections applied
      histos.add(Form("%s/NoSel/Phi", kTypeNames[type].data()), "Phi", kTH1F, {axesConfig.axisPhi});
      histos.add(Form("%s/NoSel/Eta", kTypeNames[type].data()), "Eta", kTH1F, {axesConfig.axisEta});
      histos.add(Form("%s/NoSel/DCAxy", kTypeNames[type].data()), "DCA to xy", kTH1F, {axesConfig.axisDCAxy});
      histos.add(Form("%s/NoSel/DCAz", kTypeNames[type].data()), "DCA to z", kTH1F, {axesConfig.axisDCAz});
      histos.add(Form("%s/NoSel/BachCosPA", kTypeNames[type].data()), "Bachelor cosPA", kTH1F, {{202, -1.1, 1.1}});
      histos.add(Form("%s/NoSel/V0CosPA", kTypeNames[type].data()), "V0 cosPA", kTH1F, {{202, -1.1, 1.1}});
      histos.add(Form("%s/NoSel/EvMult", kTypeNames[type].data()), "Multiplicity of events with >=1 cascade", kTH1F, {axesConfig.axisMult});
      histos.add(Form("%s/NoSel/MassXi", kTypeNames[type].data()), "Invariant mass hypothesis", kTH1F, {axesConfig.axisXiMass});
      histos.add(Form("%s/NoSel/MassOmega", kTypeNames[type].data()), "Invariant mass hypothesis", kTH1F, {axesConfig.axisOmegaMass});
      // mc truth for no selectrion
      histos.add(Form("%s/NoSel-Truth/Phi", kTypeNames[type].data()), "Phi", kTH1F, {axesConfig.axisPhi});
      histos.add(Form("%s/NoSel-Truth/Eta", kTypeNames[type].data()), "Eta", kTH1F, {axesConfig.axisEta});
      histos.add(Form("%s/NoSel-Truth/DCAxy", kTypeNames[type].data()), "DCA to xy", kTH1F, {axesConfig.axisDCAxy});
      histos.add(Form("%s/NoSel-Truth/DCAz", kTypeNames[type].data()), "DCA to z", kTH1F, {axesConfig.axisDCAz});
      histos.add(Form("%s/NoSel-Truth/DCAzVSpt", kTypeNames[type].data()), "DCA to z vs pT", kTH2F, {axesConfig.axisPt, axesConfig.axisDCAz});
      histos.add(Form("%s/NoSel-Truth/BachCosPA", kTypeNames[type].data()), "Bachelor cosPA", kTH1F, {{202, -1.1, 1.1}});
      histos.add(Form("%s/NoSel-Truth/V0CosPA", kTypeNames[type].data()), "V0 cosPA", kTH1F, {{202, -1.1, 1.1}});
      histos.add(Form("%s/NoSel-Truth/EvMult", kTypeNames[type].data()), "Multiplicity of events with >=1 cascade", kTH1F, {axesConfig.axisMult});
      histos.add(Form("%s/NoSel-Truth/MassXi", kTypeNames[type].data()), "Invariant mass hypothesis", kTH1F, {axesConfig.axisXiMass});
      histos.add(Form("%s/NoSel-Truth/MassOmega", kTypeNames[type].data()), "Invariant mass hypothesis", kTH1F, {axesConfig.axisOmegaMass});
      // passed all applied sels
      histos.add(Form("%s/Rec/MassOmega", kTypeNames[type].data()), "Invariant mass hypothesis", kTH1F, {axesConfig.axisOmegaMass});
      histos.add(Form("%s/Rec/MassXi", kTypeNames[type].data()), "Invariant mass hypothesis", kTH1F, {axesConfig.axisXiMass});
      histos.add(Form("%s/Rec/Omega", kTypeNames[type].data()), "", kTHnD, {axesConfig.axisOmegaMass, axesConfig.axisPt, axesConfig.axisMult});
      histos.add(Form("%s/Rec/Xi", kTypeNames[type].data()), "", kTHnD, {axesConfig.axisXiMass, axesConfig.axisPt, axesConfig.axisMult});
      // mc truth for all passed selections
      histos.add(Form("%s/Rec-Truth/Phi", kTypeNames[type].data()), "Phi", kTH1F, {axesConfig.axisPhi});
      histos.add(Form("%s/Rec-Truth/Eta", kTypeNames[type].data()), "Eta", kTH1F, {axesConfig.axisEta});
      histos.add(Form("%s/Rec-Truth/DCAxy", kTypeNames[type].data()), "DCA to xy", kTH1F, {axesConfig.axisDCAxy});
      histos.add(Form("%s/Rec-Truth/DCAz", kTypeNames[type].data()), "DCA to z", kTH1F, {axesConfig.axisDCAz});
      histos.add(Form("%s/Rec-Truth/DCAzVSpt", kTypeNames[type].data()), "DCA to z vs pT", kTH2F, {axesConfig.axisPt, axesConfig.axisDCAz});
      histos.add(Form("%s/Rec-Truth/BachCosPA", kTypeNames[type].data()), "Bachelor cosPA", kTH1F, {{202, -1.1, 1.1}});
      histos.add(Form("%s/Rec-Truth/V0CosPA", kTypeNames[type].data()), "V0 cosPA", kTH1F, {{202, -1.1, 1.1}});
      histos.add(Form("%s/Rec-Truth/EvMult", kTypeNames[type].data()), "Multiplicity of events with >=1 cascade", kTH1F, {axesConfig.axisMult});
      histos.add(Form("%s/Rec-Truth/MassXi", kTypeNames[type].data()), "Invariant mass hypothesis", kTH1F, {axesConfig.axisXiMass});
      histos.add(Form("%s/Rec-Truth/MassOmega", kTypeNames[type].data()), "Invariant mass hypothesis", kTH1F, {axesConfig.axisOmegaMass});
      histos.add(Form("%s/Rec-Truth/Omega", kTypeNames[type].data()), "", kTHnD, {axesConfig.axisOmegaMass, axesConfig.axisPt, axesConfig.axisMult});
      histos.add(Form("%s/Rec-Truth/Xi", kTypeNames[type].data()), "", kTHnD, {axesConfig.axisXiMass, axesConfig.axisPt, axesConfig.axisMult});
    });
    // for MC-specific processing
    histos.add("MC/Gen/EvCounter", "Event Counter", kTH1F, {{1, 0, 1}});
    histos.add("MC/Gen/Xi", "Xi", kTH2F, {axesConfig.axisPt, axesConfig.axisMult});                        // generated Xis
    histos.add("MC/Gen/Omega", "Omega", kTH2F, {axesConfig.axisPt, axesConfig.axisMult});                  // generated Omegas
    histos.add("MC/Gen/PrimaryXi", "Xi primaries", kTH2F, {axesConfig.axisPt, axesConfig.axisMult});       // generated primary Xis
    histos.add("MC/Gen/PrimaryOmega", "Omega primaries", kTH2F, {axesConfig.axisPt, axesConfig.axisMult}); // generated primary Omegas
    histos.add("MC/Rec/EvCounter", "Event Counter", kTH1F, {{1, 0, 1}});                                   // counter of all recreated events
    histos.add("MC/Rec/EvMult", "Multiplicity", kTH1F, {axesConfig.axisMult});                             // multiplicity of all recreated events
  }

  void processDerivedData(DerCollisionWMults::iterator const& collision, DerCascDatas const& allCascs, DerTraCascDatas const& traCascs, DauTracks const&)
  {
    fillEvents(collision); // save info about all processed events
    if (doApplyEfficiency1D || doApplyPurity1D || doApplyEfficiency2D || doApplyPurity2D) {
      initEfficiencyFromCCDB(collision.runNumber(), collision.timestamp());
    }
    analyseCascs(collision, allCascs); // process all cascades
    analyseCascs(collision, traCascs); // process tracked cascades
  }

  void processDerivedMCGen(aod::StraMCCollisions const& genColls, DerMCGenCascades const& genCascs, soa::Join<aod::StraCollisions, aod::StraCollLabels, aod::StraCents> const& recColls)
  {
    for (auto const& genColl : genColls) {
      histos.fill(HIST("MC/Gen/EvCounter"), 0.5); // generated events statistics
      // (for efficiency calculation) only take reconstructed events
      auto slicedRecColls = recColls.sliceBy(perMcCollision, genColl.globalIndex());
      for (auto const& recColl : slicedRecColls) {
        histos.fill(HIST("MC/Rec/EvCounter"), 0.5);
        double casc_mult = (doProcesspp || doProcesspO) ? recColl.centFT0M() : recColl.centFT0C();
        histos.fill(HIST("MC/Rec/EvMult"), casc_mult);
        int64_t genCollId = recColl.straMCCollisionId();
        for (auto const& casc : genCascs) {
          if (casc.straMCCollisionId() != genCollId)
            continue; // safety check
          double casc_pt = std::sqrt(std::pow(casc.pxMC(), 2) + std::pow(casc.pyMC(), 2));
          if (isValidPDG(casc, "xi"))
            histos.fill(HIST("MC/Gen/Xi"), casc_pt, casc_mult);
          if (isValidPDG(casc, "omega"))
            histos.fill(HIST("MC/Gen/Omega"), casc_pt, casc_mult);
          if (casc.isPhysicalPrimary()) {
            if (isValidPDG(casc, "xi"))
              histos.fill(HIST("MC/Gen/PrimaryXi"), casc_pt, casc_mult);
            if (isValidPDG(casc, "omega"))
              histos.fill(HIST("MC/Gen/PrimaryOmega"), casc_pt, casc_mult);
          }
        }
      }
    }
  }

  void processDerivedMCRec(DerMCRecCollisions::iterator const& collision, DerMCRecCascDatas const& allCascs, DerMCRecTraCascDatas const& traCascs, DauTracks const&, DerMCGenCascades const&)
  {
    fillEvents(collision); // save info about all processed events
    if (doApplyEfficiency1D || doApplyPurity1D || doApplyEfficiency2D || doApplyPurity2D) {
      initEfficiencyFromCCDB(collision.runNumber(), collision.timestamp());
    }
    analyseCascs(collision, allCascs); // process all cascades
    analyseCascs(collision, traCascs); // process tracked cascades
  }

  PROCESS_SWITCH(StrangeCascTrack, processDerivedData, "process derived data", true);
  PROCESS_SWITCH(StrangeCascTrack, processDerivedMCGen, "process derived generated mc data", false);
  PROCESS_SWITCH(StrangeCascTrack, processDerivedMCRec, "process derived reconstructed mc data", false); // mc and data are mutually exclusive!
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<StrangeCascTrack>(cfgc),
  };
}
