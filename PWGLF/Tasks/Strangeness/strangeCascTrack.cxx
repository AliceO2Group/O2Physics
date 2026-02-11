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
/// \file strangeCascTrack.cxx
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

#include <array>
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

  //* Service<o2::ccdb::BasicCCDBManager> ccdb;
  Service<o2::framework::O2DatabasePDG> pdgDB;

  PresliceUnsorted<soa::Join<aod::StraCollisions, aod::StraCollLabels, aod::StraCents>> perMcCollision = aod::v0data::straMCCollisionId;
  PresliceUnsorted<DerMCGenCascades> cascsPerMcCollision = aod::cascdata::straMCCollisionId;

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // subprocess switches:

  Configurable<bool> doProcessIons{"doProcessIons", false, "true for PbPb and OO, false for pp and pO"};

  Configurable<bool> doApplyEventCuts{"doApplyEventCuts", true, "apply general event cuts"}; // event filter - PVz, sel8, INEL>0
  // Xi selections
  Configurable<bool> doApplyPtCutsXi{"doApplyPtCutsXi", true, "apply pt cuts (Xi)"};           // ignore particles with extremely low efficiencies
  Configurable<bool> doApplyGenCutsXi{"doApplyGenCutsXi", true, "apply general cuts (Xi)"};    // general cascade cuts - cosPA, TPC hits etc.
  Configurable<bool> doApplyTPCPIDXi{"doApplyTPCPIDXi", true, "apply tpc pid to dau tracks (Xi)"};
  Configurable<bool> doApplyTOFPIDXi{"doApplyTOFPIDXi", true, "apply tof pid to dau tracks (Xi)"};
  // Omega selections
  Configurable<bool> doApplyPtCutsOmega{"doApplyPtCutsOmega", true, "apply pt cuts (Omega)"};
  Configurable<bool> doApplyGenCutsOmega{"doApplyGenCutsOmega", true, "apply general cuts (Omega)"};
  Configurable<bool> doApplyTPCPIDOmega{"doApplyTPCPIDOmega", true, "apply tpc pid to dau tracks (Omega)"};
  Configurable<bool> doApplyTOFPIDOmega{"doApplyTOFPIDOmega", true, "apply tof pid to dau tracks (Omega)"};
  Configurable<bool> doCompetingMassRej{"doCompetingMassRej", true, "competing mass rejection (Omega)"};
  //* // efficiency and purity corrections on the fly (warning: to be avoided because interpolation causes errors):
  //* // only correct by pt
  //* Configurable<bool> doApplyEfficiency1D{"doApplyEfficiency1D", false, "apply efficiency correction"};
  //* Configurable<bool> doPropagateEfficiency1D{"doPropagateEfficiency1D", false, "apply efficiency propagation"};
  //* Configurable<bool> doApplyPurity1D{"doApplyPurity1D", false, "apply purity correction"};
  //* Configurable<bool> doPropagatePurity1D{"doPropagatePurity1D", false, "apply purity propagation"};
  //* // correct by both pt and mult
  //* Configurable<bool> doApplyEfficiency2D{"doApplyEfficiency2D", false, "apply efficiency correction"};
  //* Configurable<bool> doPropagateEfficiency2D{"doPropagateEfficiency2D", false, "apply efficiency propagation"};
  //* Configurable<bool> doApplyPurity2D{"doApplyPurity2D", false, "apply purity correction"};
  //* Configurable<bool> doPropagatePurity2D{"doPropagatePurity2D", false, "apply purity propagation"};
  //* Configurable<std::string> ccdburl{"ccdburl", "http://alice-ccdb.cern.ch", "url of the ccdb repository to use"};
  //* Configurable<std::string> efficiencyCCDBPathpp{"efficiencyCCDBPathpp", "Users/y/yparovia/LHC24f4d", "Path of the efficiency and purity corrections (pp)"};
  //* Configurable<std::string> efficiencyCCDBPathPbPb{"efficiencyCCDBPathPbPb", "Users/y/yparovia/LHC25f3", "Path of the efficiency and purity corrections (PbPb)"};
  //* Configurable<std::string> efficiencyCCDBPathOO{"efficiencyCCDBPathOO", "Users/y/yparovia/LHC25h3", "Path of the efficiency and purity corrections (OO)"};
  //* Configurable<std::string> efficiencyCCDBPathpO{"efficiencyCCDBPathpO", "Users/y/yparovia/LHC25h2", "Path of the efficiency and purity corrections (pO)"};

  // axes
  struct : ConfigurableGroup {
    ConfigurableAxis axisEta{"axisEta", {102, -2.01, 2.01}, "#eta"};
    ConfigurableAxis axisDCAxy{"axisDCAxy", {500, 0., 0.5}, "cm"};
    ConfigurableAxis axisDCAz{"axisDCAz", {500, 0., 0.5}, "cm"};
    ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 10.0}, "p_{T} (GeV/c)"};
    ConfigurableAxis axisMult{"axisMult", {VARIABLE_WIDTH, 0.0f, 5.0, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 70.0f, 100.0f}, "FT0 mult %"};
    ConfigurableAxis axisOmegaMass{"axisOmegaMass", {2000, 1.6, 1.8}, "#Omega M_{inv} (GeV/c^{2})"};
    ConfigurableAxis axisXiMass{"axisXiMass", {2000, 1.2, 1.4}, "#Xi M_{inv} (GeV/c^{2})"};
  } axesConfig;

  // event and dau track selection
  struct : ConfigurableGroup {
    // event cuts
    Configurable<bool> cutDoINEL{"cutDoINEL", true, "choose events with INEL>0"};
    Configurable<float> cutZVertex{"cutZVertex", 10.0f, "max Z-vertex position"};
    Configurable<bool> cutDoSel8{"cutDoSel8", true, "choose events with sel8"};
    Configurable<bool> cutDoNoPileup{"cutDoNoPileup", true, "choose events with no bunch pileup"};
    Configurable<bool> cutDoGoodFT0 {"cutDoGoodFT0", true, "choose events with IsGoodZvtxFT0VsPV"};
    // cascade cuts
    Configurable<bool> cutDoPropagateDCA{"cutDoPropagateDCA", false, "choose events with sel8"};
    Configurable<float> cutPropDCAtoPVxy{"cutPropDCAtoPVxy", 0.02f, "max cascade dca to PV in xy - propagated"};
    Configurable<float> cutPropDCAtoPVz{"cutPropDCAtoPVz", 0.02f, "max cascade dca to PV in z - propagated"};
    Configurable<int> cutNClsTPC{"cutNClsTPC", 70, "min number of found TPC clusters for dau tracks"};
    Configurable<float> cutMinV0CosPA{"cutMinV0CosPA", -1.1f, "min V0 cosPA"};
    Configurable<float> cutMaxV0CosPA{"cutMaxV0CosPA", 1.1f, "max V0 cosPA"};
    Configurable<float> cutMinBachCosPA{"cutMinBachCosPA", -1.1f, "min Bachelor cosPA"};
    Configurable<float> cutMaxBachCosPA{"cutMaxBachCosPA", 1.1f, "max Bachelor cosPA"};
    Configurable<std::vector<float>> cutMinCascCosPaVsPt{
      "cutMinCascCosPaVsPt",
      {0.993, 0.994, 0.995, 0.996, 0.997, 0.997, 0.998, 0.998, 0.999, 0.999, 0.999},
      "Min Casc CosPA per pT bin (same binning as axisPt)"
    };
    Configurable<float> cutRapidity{"cutRapidity", 0.5f, "max rapidity"};
    Configurable<float> cutDauEta{"cutDauEta", 1.0f, "max eta of dau tracks"};
    Configurable<float> cutCompMassRej{"cutCompMassRej", 0.008f, "Competing mass rejection"};
    // minimum and maximum desired pt
    Configurable<float> cutMinPtXiStd{"cutMinPtXiStd", 0.0f, "min pt for standard Xi"};
    Configurable<float> cutMaxPtXiStd{"cutMaxPtXiStd", 15.0f, "min pt for standard Xi"};
    Configurable<float> cutMinPtXiTra{"cutMinPtXiTra", 0.5f, "min pt for tracked Xi"};
    Configurable<float> cutMaxPtXiTra{"cutMaxPtXiTra", 15.0f, "min pt for standard Xi"};
    Configurable<float> cutMinPtOmegaStd{"cutMinPtOmegaStd", 0.5f, "min pt for standard Omega"};
    Configurable<float> cutMaxPtOmegaStd{"cutMaxPtOmegaStd", 15.0f, "min pt for standard Omega"};
    Configurable<float> cutMinPtOmegaTra{"cutMinPtOmegaTra", 1.0f, "min pt for tracked Omega"};
    Configurable<float> cutMaxPtOmegaTra{"cutMaxPtOmegaTra", 15.0f, "min pt for standard Omega"};
    // TPC PID selection
    Configurable<float> cutNSigmaTPCPion{"cutNSigmaTPCPion", 4, "cutNSigmaTPCPion"};
    Configurable<float> cutNSigmaTPCKaon{"cutNSigmaTPCKaon", 4, "cutNSigmaTPCKaon"};
    Configurable<float> cutNSigmaTPCProton{"cutNSigmaTPCProton", 4, "cutNSigmaTPCProton"};
    // TOF PID selection
    Configurable<float> cutNSigmaTOFXi{"cutNSigmaTOFXi", 3, "cutNSigmaTOFXi"};
    Configurable<float> cutNSigmaTOFOmega{"cutNSigmaTOFOmega", 3, "cutNSigmaTOFOmega"};
  } selCuts;

  // cascade reconstruction Types
  static constexpr std::string_view TypeNames[] = {"Standard", "Tracked"};

  //* // for on-the-fly efficiency and purity corrections - working but deprecated
  //* TH1F* hEfficiencyOmegaStd1D;
  //* TH1F* hEfficiencyOmegaTra1D;
  //* TH1F* hEfficiencyXiStd1D;
  //* TH1F* hEfficiencyXiTra1D;
  //*
  //* TH1F* hEfficiencyErrOmegaStd1D;
  //* TH1F* hEfficiencyErrOmegaTra1D;
  //* TH1F* hEfficiencyErrXiStd1D;
  //* TH1F* hEfficiencyErrXiTra1D;
  //*
  //* TH1F* hPurityOmegaStd1D;
  //* TH1F* hPurityOmegaTra1D;
  //* TH1F* hPurityXiStd1D;
  //* TH1F* hPurityXiTra1D;
  //*
  //* TH1F* hPurityErrOmegaStd1D;
  //* TH1F* hPurityErrOmegaTra1D;
  //* TH1F* hPurityErrXiStd1D;
  //* TH1F* hPurityErrXiTra1D;
  //*
  //* TH2F* hEfficiencyOmegaStd2D;
  //* TH2F* hEfficiencyOmegaTra2D;
  //* TH2F* hEfficiencyXiStd2D;
  //* TH2F* hEfficiencyXiTra2D;
  //*
  //* TH2F* hEfficiencyErrOmegaStd2D;
  //* TH2F* hEfficiencyErrOmegaTra2D;
  //* TH2F* hEfficiencyErrXiStd2D;
  //* TH2F* hEfficiencyErrXiTra2D;
  //*
  //* TH2F* hPurityOmegaStd2D;
  //* TH2F* hPurityOmegaTra2D;
  //* TH2F* hPurityXiStd2D;
  //* TH2F* hPurityXiTra2D;
  //*
  //* TH2F* hPurityErrOmegaStd2D;
  //* TH2F* hPurityErrOmegaTra2D;
  //* TH2F* hPurityErrXiStd2D;
  //* TH2F* hPurityErrXiTra2D;
  //*
  //* int mRunNumber;
  //* // loads efficiencies and purities
  //* void initEfficiencyFromCCDB(int64_t runNumber, int64_t timestamp)
  //* {
  //*   if (mRunNumber == runNumber) {
  //*     return;
  //*   }
  //*   mRunNumber = runNumber;
  //*   LOG(info) << "Loading efficiencies and purities from CCDB for run " << mRunNumber << " now...";
  //*   auto timeStamp = timestamp;
  //*
  //*   std::string efficiencyCCDBPath = [&]() {
  //*     if (doProcesspp) {
  //*       return efficiencyCCDBPathpp;
  //*     } else if (doProcesspO) {
  //*       return efficiencyCCDBPathpO;
  //*     } else if (doProcessPbPb) {
  //*       return efficiencyCCDBPathPbPb;
  //*     }
  //*     return efficiencyCCDBPathOO;
  //*   }();
  //*
  //*   TList* listEfficiencies = ccdb->getForTimeStamp<TList>(efficiencyCCDBPath, timeStamp);
  //*
  //*   if (!listEfficiencies) {
  //*     LOG(fatal) << "Problem getting TList object with efficiencies and purities!";
  //*   }
  //*
  //*   hEfficiencyOmegaStd1D = static_cast<TH1F*>(listEfficiencies->FindObject("Eff_Omega_Standard_byPt"));
  //*   hEfficiencyOmegaTra1D = static_cast<TH1F*>(listEfficiencies->FindObject("Eff_Omega_Tracked_byPt"));
  //*   hEfficiencyXiStd1D = static_cast<TH1F*>(listEfficiencies->FindObject("Eff_Xi_Standard_byPt"));
  //*   hEfficiencyXiTra1D = static_cast<TH1F*>(listEfficiencies->FindObject("Eff_Xi_Tracked_byPt"));
  //*   hEfficiencyErrOmegaStd1D = static_cast<TH1F*>(listEfficiencies->FindObject("EffErr_Omega_Standard_byPt"));
  //*   hEfficiencyErrOmegaTra1D = static_cast<TH1F*>(listEfficiencies->FindObject("EffErr_Omega_Tracked_byPt"));
  //*   hEfficiencyErrXiStd1D = static_cast<TH1F*>(listEfficiencies->FindObject("EffErr_Xi_Standard_byPt"));
  //*   hEfficiencyErrXiTra1D = static_cast<TH1F*>(listEfficiencies->FindObject("EffErr_Xi_Tracked_byPt"));
  //*   hPurityOmegaStd1D = static_cast<TH1F*>(listEfficiencies->FindObject("Pur_Omega_Standard_byPt"));
  //*   hPurityOmegaTra1D = static_cast<TH1F*>(listEfficiencies->FindObject("Pur_Omega_Tracked_byPt"));
  //*   hPurityXiStd1D = static_cast<TH1F*>(listEfficiencies->FindObject("Pur_Xi_Standard_byPt"));
  //*   hPurityXiTra1D = static_cast<TH1F*>(listEfficiencies->FindObject("Pur_Xi_Tracked_byPt"));
  //*   hPurityErrOmegaStd1D = static_cast<TH1F*>(listEfficiencies->FindObject("PurErr_Omega_Standard_byPt"));
  //*   hPurityErrOmegaTra1D = static_cast<TH1F*>(listEfficiencies->FindObject("PurErr_Omega_Tracked_byPt"));
  //*   hPurityErrXiStd1D = static_cast<TH1F*>(listEfficiencies->FindObject("PurErr_Xi_Standard_byPt"));
  //*   hPurityErrXiTra1D = static_cast<TH1F*>(listEfficiencies->FindObject("PurErr_Xi_Tracked_byPt"));
  //*
  //*   hEfficiencyOmegaStd2D = static_cast<TH2F*>(listEfficiencies->FindObject("Eff_Omega_Standard_byPtMult"));
  //*   hEfficiencyOmegaTra2D = static_cast<TH2F*>(listEfficiencies->FindObject("Eff_Omega_Tracked_byPtMult"));
  //*   hEfficiencyXiStd2D = static_cast<TH2F*>(listEfficiencies->FindObject("Eff_Xi_Standard_byPtMult"));
  //*   hEfficiencyXiTra2D = static_cast<TH2F*>(listEfficiencies->FindObject("Eff_Xi_Tracked_byPtMult"));
  //*   hEfficiencyErrOmegaStd2D = static_cast<TH2F*>(listEfficiencies->FindObject("EffErr_Omega_Standard_byPtMult"));
  //*   hEfficiencyErrOmegaTra2D = static_cast<TH2F*>(listEfficiencies->FindObject("EffErr_Omega_Tracked_byPtMult"));
  //*   hEfficiencyErrXiStd2D = static_cast<TH2F*>(listEfficiencies->FindObject("EffErr_Xi_Standard_byPtMult"));
  //*   hEfficiencyErrXiTra2D = static_cast<TH2F*>(listEfficiencies->FindObject("EffErr_Xi_Tracked_byPtMult"));
  //*   hPurityOmegaStd2D = static_cast<TH2F*>(listEfficiencies->FindObject("Pur_Omega_Standard_byPtMult"));
  //*   hPurityOmegaTra2D = static_cast<TH2F*>(listEfficiencies->FindObject("Pur_Omega_Tracked_byPtMult"));
  //*   hPurityXiStd2D = static_cast<TH2F*>(listEfficiencies->FindObject("Pur_Xi_Standard_byPtMult"));
  //*   hPurityXiTra2D = static_cast<TH2F*>(listEfficiencies->FindObject("Pur_Xi_Tracked_byPtMult"));
  //*   hPurityErrOmegaStd2D = static_cast<TH2F*>(listEfficiencies->FindObject("PurErr_Omega_Standard_byPtMult"));
  //*   hPurityErrOmegaTra2D = static_cast<TH2F*>(listEfficiencies->FindObject("PurErr_Omega_Tracked_byPtMult"));
  //*   hPurityErrXiStd2D = static_cast<TH2F*>(listEfficiencies->FindObject("PurErr_Xi_Standard_byPtMult"));
  //*   hPurityErrXiTra2D = static_cast<TH2F*>(listEfficiencies->FindObject("PurErr_Xi_Tracked_byPtMult"));
  //*
  //*   if (doPropagateEfficiency1D && (!hEfficiencyErrOmegaStd1D || !hEfficiencyErrOmegaTra1D || !hEfficiencyErrXiStd1D || !hEfficiencyErrXiTra1D))
  //*     LOG(fatal) << "Problem getting hEfficiencyUncertainty!";
  //*   if (doPropagatePurity1D && (!hPurityErrOmegaStd1D || !hPurityErrOmegaTra1D || !hPurityErrXiStd1D || !hPurityErrXiTra1D))
  //*     LOG(fatal) << "Problem getting hPurityUncertainty!";
  //*   LOG(info) << "Efficiencies and purities now loaded for " << mRunNumber;
  //*
  //*   if (doPropagateEfficiency2D && (!hEfficiencyErrOmegaStd2D || !hEfficiencyErrOmegaTra2D || !hEfficiencyErrXiStd2D || !hEfficiencyErrXiTra2D))
  //*     LOG(fatal) << "Problem getting hEfficiencyUncertainty!";
  //*   if (doPropagatePurity2D && (!hPurityErrOmegaStd2D || !hPurityErrOmegaTra2D || !hPurityErrXiStd2D || !hPurityErrXiTra2D))
  //*     LOG(fatal) << "Problem getting hPurityUncertainty!";
  //*   LOG(info) << "Efficiencies and purities now loaded for " << mRunNumber;
  //* }

  // general info about processed events
  template <typename TEvent>
  void fillEvents(TEvent const& collision)
  {
    histos.fill(HIST("NoSel-Events/EvCounter"), 0.5);
    double mult = (doProcessIons) ? collision.centFT0C() : collision.centFT0M();
    histos.fill(HIST("NoSel-Events/Mult"), mult);
    double pvx = collision.posX();
    double pvy = collision.posY();
    double pvz = collision.posZ();
    histos.fill(HIST("NoSel-Events/PVxy"), pvx, pvy);
    histos.fill(HIST("NoSel-Events/PVz"), pvz);
  }
  // checks general selection criteria for collisions
  template <typename TEvent>
  bool isValidEvent(TEvent collision)
  {
    bool passedAllSels = true;
    //* inel>0 cut
    if (!selCuts.cutDoINEL || collision.multNTracksPVeta1() > 0)
      histos.fill(HIST("Rec-Events/EvFilter"), 0.5);
    else
      passedAllSels = false;
    //* pvz cut
    if (std::abs(collision.posZ()) < selCuts.cutZVertex)
      histos.fill(HIST("Rec-Events/EvFilter"), 1.5);
    else
      passedAllSels = false;
    //* sel8 cut
    if (!selCuts.cutDoSel8 || collision.sel8())
      histos.fill(HIST("Rec-Events/EvFilter"), 2.5);
    else
      passedAllSels = false;
    //* pileup cut
    if (!selCuts.cutDoNoPileup || collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup))
      histos.fill(HIST("Rec-Events/EvFilter"), 3.5);
    else
      passedAllSels = false;
    //* good ft0 z-vertex vs pv cut
    if (!selCuts.cutDoGoodFT0 || collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV))
      histos.fill(HIST("Rec-Events/EvFilter"), 4.5);
    else
      passedAllSels = false;
    //* all cuts
    if (passedAllSels)
      histos.fill(HIST("Rec-Events/EvFilter"), 5.5);
    return passedAllSels;
  }
  // checks cascade pt
  template <typename TCascade>
  bool isValidPt(TCascade cascade, TString particle, int Type)
  {
    bool passedSel = true;
    double ptMin = 0.0;
    double ptMax = 0.0;
    if (Type == 1 && particle == "Xi") {
      ptMin = selCuts.cutMinPtXiTra;
      ptMax = selCuts.cutMaxPtXiTra;
    }
    if (Type == 1 && particle == "Omega") {
      ptMin = selCuts.cutMinPtOmegaTra;
      ptMax = selCuts.cutMaxPtOmegaTra;
    }
    if (Type == 0 && particle == "Xi") {
      ptMin = selCuts.cutMinPtXiStd;
      ptMax = selCuts.cutMaxPtXiStd;
    }
    if (Type == 0 && particle == "Omega") {
      ptMin = selCuts.cutMinPtOmegaStd;
      ptMax = selCuts.cutMaxPtOmegaStd;
    }
    if (cascade.pt() < ptMin || cascade.pt() > ptMax)
      passedSel = false;
    //*  histos.fill(HIST(Form("%s/Rec/Filters%s", TypeNames[Type].data(), particle)), 0.5);
    return passedSel;
  }
  // checks general selection criteria for cascades
  template <typename TEvent, typename TCascade, typename TStdCascade>
  std::array<bool, 9> isValidCasc(TEvent collision, TCascade cascade, TStdCascade stdcasc, TString particle)
  {
    bool passedAllSels = true;
    // cascade rapidity
    bool passedRapidity = true;
    double y;
    if (particle == "Xi")
      y = std::abs(cascade.yXi());
    else
      y = std::abs(cascade.yOmega());
    if (y > selCuts.cutRapidity) {
      passedRapidity = false;
      passedAllSels = false;
    }
    // daughter track pseudorapidity
    bool passedDauEta = true;
    double bachEta = std::abs(cascade.bacheloreta());
    double negEta = std::abs(cascade.negativeeta());
    double posEta = std::abs(cascade.positiveeta());
    if (bachEta > selCuts.cutDauEta || negEta > selCuts.cutDauEta || posEta > selCuts.cutDauEta) {
      passedDauEta = false;
      passedAllSels = false;
    }
    // daughter found TPC clusters
    bool passedTPCCls = true;
    const auto& posTrack = stdcasc.template posTrackExtra_as<DauTracks>();
    const auto& negTrack = stdcasc.template negTrackExtra_as<DauTracks>();
    const auto& bachTrack = stdcasc.template bachTrackExtra_as<DauTracks>();
    double posCls = posTrack.tpcClusters();
    double negCls = negTrack.tpcClusters();
    double bachCls = bachTrack.tpcClusters();
    if (posCls < selCuts.cutNClsTPC || negCls < selCuts.cutNClsTPC || bachCls < selCuts.cutNClsTPC) {
      passedTPCCls = false;
      passedAllSels = false;
    }
    // V0 cosPA
    bool passedV0CosPA = true;
    double v0cospa = cascade.v0cosPA(collision.posX(), collision.posY(), collision.posZ());
    if (v0cospa < selCuts.cutMinV0CosPA || v0cospa > selCuts.cutMaxV0CosPA) {
      passedV0CosPA = false;
      passedAllSels = false;
    }
    // Bachelor cosPA
    bool passedBachCosPA = true;
    double bachcospa = stdcasc.bachBaryonCosPA();
    if (bachcospa < selCuts.cutMinBachCosPA || bachcospa > selCuts.cutMaxBachCosPA) {
      passedBachCosPA = false;
      passedAllSels = false;
    }
    // Cascade cosPA
    bool passedCascCosPA = true;
    double casccospa = cascade.casccosPA(collision.posX(), collision.posY(), collision.posZ());
    const auto& edges = axesConfig.axisPt.value;
    int ptBin = std::upper_bound(edges.begin(), edges.end(), cascade.pt())
                - edges.begin() - 1;
    if (ptBin < 0 || ptBin >= (int)selCuts.cutMinCascCosPaVsPt->size()) {
      ptBin = 1;
    } // safety check - if pt bin not determined, default to loosest cut 
    if (casccospa < selCuts.cutMinCascCosPaVsPt->at(ptBin)) {
      passedCascCosPA = false;
      passedAllSels = false;
    }
    // Propagated cacade DCAxy to PV
    bool passedPropDCAxy = true;
    if (selCuts.cutDoPropagateDCA && cascade.dcaXYCascToPV() > selCuts.cutPropDCAtoPVxy) {
      passedPropDCAxy = false;
      passedAllSels = false;
    }
    // Propagated cacade DCAz to PV
    bool passedPropDCAz = true;
    if (selCuts.cutDoPropagateDCA && cascade.dcaZCascToPV() > selCuts.cutPropDCAtoPVz) {
      passedPropDCAz = false;
      passedAllSels = false;
    }
    return {passedRapidity, passedDauEta, passedTPCCls, passedV0CosPA, passedBachCosPA, passedCascCosPA, passedPropDCAxy, passedPropDCAz, passedAllSels};
  }
  // checks TPC PID of dau tracks
  template <typename TCascade>
  bool passesTPC(TCascade cascade, TString particle)
  {
    bool passedSel = true;
    const auto& posTrack = cascade.template posTrackExtra_as<DauTracks>();
    const auto& negTrack = cascade.template negTrackExtra_as<DauTracks>();
    const auto& bachTrack = cascade.template bachTrackExtra_as<DauTracks>();
    if (cascade.sign() < 0) {
      if (std::abs(posTrack.tpcNSigmaPr()) > selCuts.cutNSigmaTPCProton) {
        passedSel = false;
      }
      if (std::abs(negTrack.tpcNSigmaPi()) > selCuts.cutNSigmaTPCPion) {
        passedSel = false;
      }
      if ((particle == "Xi" && std::abs(bachTrack.tpcNSigmaPi()) > selCuts.cutNSigmaTPCPion) ||
          (particle == "Omega" && std::abs(bachTrack.tpcNSigmaKa()) > selCuts.cutNSigmaTPCKaon)) {
        passedSel = false;
      }
    } else {
      if (std::abs(negTrack.tpcNSigmaPr()) > selCuts.cutNSigmaTPCProton) {
        passedSel = false;
      }
      if (std::abs(posTrack.tpcNSigmaPi()) > selCuts.cutNSigmaTPCPion) {
        passedSel = false;
      }
      if ((particle == "Xi" && std::abs(bachTrack.tpcNSigmaPi()) > selCuts.cutNSigmaTPCPion) ||
          (particle == "Omega" && std::abs(bachTrack.tpcNSigmaKa()) > selCuts.cutNSigmaTPCKaon)) {
        passedSel = false;
      }
    }
    return passedSel;
  }
  // checks TOF PID of dau tracks
  template <typename TCascade>
  bool passesTOF(TCascade cascade, TString particle)
  {
    bool passedSel = true;
    if (particle == "Xi")
      passedSel = cascade.tofXiCompatibility(selCuts.cutNSigmaTOFXi);
    if (particle == "Omega")
      passedSel = cascade.tofOmegaCompatibility(selCuts.cutNSigmaTOFOmega);
    return passedSel;
  }
  // checks whether gen cascade corresponds to PDG code
  template <typename TCascade>
  bool isValidPDG(TCascade cascade, TString particle)
  {
    if (particle == "Xi" && std::abs(cascade.pdgCode()) == PDG_t::kXiMinus)
      return true;
    if (particle == "Omega" && std::abs(cascade.pdgCode()) == PDG_t::kOmegaMinus)
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
      if (particle == "Xi")
        return (pdg == PDG_t::kXiMinus);
      if (particle == "Omega")
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

  ////* calculates DCA from cosPA - deprecated due to low precision
  ////* template <typename TEvent>
  ////* double calculateDCA(TEvent collision, double cosPA, double decX, double decY, double decZ)
  ////* {
  ////*   double pvX = collision.posX();
  ////*   double pvY = collision.posX();
  ////*   double pvZ = collision.posX();
  ////*   double sinPA = std::sqrt(1 - cosPA * cosPA);
  ////*   double dca = sinPA * std::sqrt(std::pow(decX - pvX, 2) + std::pow(decY - pvY, 2) + std::pow(decZ - pvZ, 2));
  ////*   return dca;
  ////* }

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

      // Type 1 for tracked cascades, Type 0 for standard
      static constexpr int Type = [&]() {
        if constexpr (requires { cascade.topologyChi2(); }) {
          return 1;
        } else {
          return 0;
        }
      }();

      float efficiencyOmega = 1.0f;
      float efficiencyXi = 1.0f;
      float efficiencyOmegaErr = 0.0f;
      float efficiencyXiErr = 0.0f;
      float purityOmega = 1.0f;
      float purityXi = 1.0f;
      float purityOmegaErr = 0.0f;
      float purityXiErr = 0.0f;

      double mult = (doProcessIons) ? collision.centFT0C() : collision.centFT0M(); // ion collisions use FT0C for multiplicity, pp uses both

      //* if (doApplyEfficiency1D) {
      //*   if constexpr (requires { cascade.topologyChi2(); }) {
      //*     efficiencyOmega = hEfficiencyOmegaTra1D->Interpolate(cascade.pt());
      //*     efficiencyXi = hEfficiencyXiTra1D->Interpolate(cascade.pt());
      //*     if (doPropagateEfficiency1D) {
      //*       efficiencyOmegaErr = hEfficiencyErrOmegaTra1D->Interpolate(cascade.pt());
      //*       efficiencyXiErr = hEfficiencyErrXiTra1D->Interpolate(cascade.pt());
      //*     }
      //*     if (efficiencyOmega == 0) { // check for zero efficiency, do not apply if the case
      //*       efficiencyOmega = 1.;
      //*       efficiencyOmegaErr = 0.;
      //*     }
      //*     if (efficiencyXi == 0) { // check for zero efficiency, do not apply if the case
      //*       efficiencyXi = 1.;
      //*       efficiencyXiErr = 0.;
      //*     }
      //*   } else {
      //*     efficiencyOmega = hEfficiencyOmegaStd1D->Interpolate(cascade.pt());
      //*     efficiencyXi = hEfficiencyXiStd1D->Interpolate(cascade.pt());
      //*     if (doPropagateEfficiency1D) {
      //*       efficiencyOmegaErr = hEfficiencyErrOmegaStd1D->Interpolate(cascade.pt());
      //*       efficiencyXiErr = hEfficiencyErrXiStd1D->Interpolate(cascade.pt());
      //*     }
      //*     if (efficiencyOmega == 0) { // check for zero efficiency, do not apply if the case
      //*       efficiencyOmega = 1.;
      //*       efficiencyOmegaErr = 0.;
      //*     }
      //*     if (efficiencyXi == 0) { // check for zero efficiency, do not apply if the case
      //*       efficiencyXi = 1.;
      //*       efficiencyXiErr = 0.;
      //*     }
      //*   }
      //* }
      //*
      //* if (doApplyEfficiency2D) {
      //*   if constexpr (requires { cascade.topologyChi2(); }) {
      //*     efficiencyOmega = hEfficiencyOmegaTra2D->Interpolate(cascade.pt(), mult);
      //*     efficiencyXi = hEfficiencyXiTra2D->Interpolate(cascade.pt(), mult);
      //*     if (doPropagateEfficiency2D) {
      //*       efficiencyOmegaErr = hEfficiencyErrOmegaTra2D->Interpolate(cascade.pt(), mult);
      //*       efficiencyXiErr = hEfficiencyErrXiTra2D->Interpolate(cascade.pt(), mult);
      //*     }
      //*     if (efficiencyOmega == 0) { // check for zero efficiency, do not apply if the case
      //*       efficiencyOmega = 1.;
      //*       efficiencyOmegaErr = 0.;
      //*     }
      //*     if (efficiencyXi == 0) { // check for zero efficiency, do not apply if the case
      //*       efficiencyXi = 1.;
      //*       efficiencyXiErr = 0.;
      //*     }
      //*   } else {
      //*     efficiencyOmega = hEfficiencyOmegaStd2D->Interpolate(cascade.pt(), mult);
      //*     efficiencyXi = hEfficiencyXiStd2D->Interpolate(cascade.pt(), mult);
      //*     if (doPropagateEfficiency2D) {
      //*       efficiencyOmegaErr = hEfficiencyErrOmegaStd2D->Interpolate(cascade.pt(), mult);
      //*       efficiencyXiErr = hEfficiencyErrXiStd2D->Interpolate(cascade.pt(), mult);
      //*     }
      //*     if (efficiencyOmega == 0) { // check for zero efficiency, do not apply if the case
      //*       efficiencyOmega = 1.;
      //*       efficiencyOmegaErr = 0.;
      //*     }
      //*     if (efficiencyXi == 0) { // check for zero efficiency, do not apply if the case
      //*       efficiencyXi = 1.;
      //*       efficiencyXiErr = 0.;
      //*     }
      //*   }
      //* }
      //*
      //* if (doApplyPurity1D) {
      //*   if constexpr (requires { cascade.topologyChi2(); }) {
      //*     purityOmega = hPurityOmegaTra1D->Interpolate(cascade.pt());
      //*     purityXi = hPurityXiTra1D->Interpolate(cascade.pt());
      //*     if (doPropagatePurity1D) {
      //*       purityOmegaErr = hPurityErrOmegaTra1D->Interpolate(cascade.pt());
      //*       purityXiErr = hPurityErrXiTra1D->Interpolate(cascade.pt());
      //*     }
      //*     if (purityOmega == 0) { // check for zero purity, do not apply if the case
      //*       purityOmega = 1.;
      //*       purityOmegaErr = 0.;
      //*     }
      //*     if (purityXi == 0) { // check for zero purity, do not apply if the case
      //*       purityXi = 1.;
      //*       purityXiErr = 0.;
      //*     }
      //*   } else {
      //*     purityOmega = hPurityOmegaStd1D->Interpolate(cascade.pt());
      //*     purityXi = hPurityXiStd1D->Interpolate(cascade.pt());
      //*     if (doPropagatePurity1D) {
      //*       purityOmegaErr = hPurityErrOmegaStd1D->Interpolate(cascade.pt());
      //*       purityXiErr = hPurityErrXiStd1D->Interpolate(cascade.pt());
      //*     }
      //*     if (purityOmega == 0) { // check for zero purity, do not apply if the case
      //*       purityOmega = 1.;
      //*       purityOmegaErr = 0.;
      //*     }
      //*     if (purityXi == 0) { // check for zero purity, do not apply if the case
      //*       purityXi = 1.;
      //*       purityXiErr = 0.;
      //*     }
      //*   }
      //* }
      //*
      //* if (doApplyPurity2D) {
      //*   if constexpr (requires { cascade.topologyChi2(); }) {
      //*     purityOmega = hPurityOmegaTra2D->Interpolate(cascade.pt(), mult);
      //*     purityXi = hPurityXiTra2D->Interpolate(cascade.pt(), mult);
      //*     if (doPropagatePurity2D) {
      //*       purityOmegaErr = hPurityErrOmegaTra2D->Interpolate(cascade.pt(), mult);
      //*       purityXiErr = hPurityErrXiTra2D->Interpolate(cascade.pt(), mult);
      //*     }
      //*     if (purityOmega == 0) { // check for zero purity, do not apply if the case
      //*       purityOmega = 1.;
      //*       purityOmegaErr = 0.;
      //*     }
      //*     if (purityXi == 0) { // check for zero purity, do not apply if the case
      //*       purityXi = 1.;
      //*       purityXiErr = 0.;
      //*     }
      //*   } else {
      //*     purityOmega = hPurityOmegaStd2D->Interpolate(cascade.pt(), mult);
      //*     purityXi = hPurityXiStd2D->Interpolate(cascade.pt(), mult);
      //*     if (doPropagatePurity2D) {
      //*       purityOmegaErr = hPurityErrOmegaStd2D->Interpolate(cascade.pt(), mult);
      //*       purityXiErr = hPurityErrXiStd2D->Interpolate(cascade.pt(), mult);
      //*     }
      //*     if (purityOmega == 0) { // check for zero purity, do not apply if the case
      //*       purityOmega = 1.;
      //*       purityOmegaErr = 0.;
      //*     }
      //*     if (purityXi == 0) { // check for zero purity, do not apply if the case
      //*       purityXi = 1.;
      //*       purityXiErr = 0.;
      //*     }
      //*   }
      //* }

      // fill multiplicity for events with >=1 cascade
      if (collision.index() != casccollid) {
        histos.fill(HIST(TypeNames[Type]) + HIST("/NoSel/EvMult"), mult);
        if constexpr (requires { collision.straMCCollisionId(); }) {
          if (isMCTruth(stdCasc, "Xi") || isMCTruth(stdCasc, "Omega")) {
            histos.fill(HIST(TypeNames[Type]) + HIST("/NoSel-Truth/EvMult"), mult);
          }
        }
        casccollid = collision.index();
      }

      double massXi = cascade.mXi();
      double massOmega = cascade.mOmega();
      double pt = cascade.pt();
      double v0cosPA = cascade.v0cosPA(collision.posX(), collision.posY(), collision.posZ());
      double casccosPA = cascade.casccosPA(collision.posX(), collision.posY(), collision.posZ());
      ////* double calcDCA = calculateDCA(collision, casccosPA, cascade.x(), cascade.y(), cascade.z());
      double bachEta = cascade.bacheloreta();
      double negEta = cascade.negativeeta();
      double posEta = cascade.positiveeta();

      // fill filters for no cascade selections
      histos.fill(HIST(TypeNames[Type]) + HIST("/NoSel/Filters/PropDCAxy"), cascade.dcaXYCascToPV());
      histos.fill(HIST(TypeNames[Type]) + HIST("/NoSel/Filters/PropDCAz"), cascade.dcaZCascToPV());
      ////* histos.fill(HIST(TypeNames[Type]) + HIST("/NoSel/Filters/CalcDCA"), calcDCA);
      histos.fill(HIST(TypeNames[Type]) + HIST("/NoSel/Filters/BachCosPA"), stdCasc.bachBaryonCosPA());
      histos.fill(HIST(TypeNames[Type]) + HIST("/NoSel/Filters/V0CosPA"), v0cosPA);
      histos.fill(HIST(TypeNames[Type]) + HIST("/NoSel/Filters/CascCosPA"), casccosPA);
      histos.fill(HIST(TypeNames[Type]) + HIST("/NoSel/Filters/RapidityXi"), cascade.yXi());
      histos.fill(HIST(TypeNames[Type]) + HIST("/NoSel/Filters/PtVsRapidityXi"), pt, cascade.yXi());
      histos.fill(HIST(TypeNames[Type]) + HIST("/NoSel/Filters/RapidityOmega"), cascade.yOmega());
      histos.fill(HIST(TypeNames[Type]) + HIST("/NoSel/Filters/PtVsRapidityOmega"), pt, cascade.yOmega());
      histos.fill(HIST(TypeNames[Type]) + HIST("/NoSel/Filters/EtaDau"), bachEta);
      histos.fill(HIST(TypeNames[Type]) + HIST("/NoSel/Filters/EtaDau"), negEta);
      histos.fill(HIST(TypeNames[Type]) + HIST("/NoSel/Filters/EtaDau"), posEta);
      // fill inv mass for no cascade selections
      histos.fill(HIST(TypeNames[Type]) + HIST("/NoSel/MassXi"), massXi);
      histos.fill(HIST(TypeNames[Type]) + HIST("/NoSel/MassOmega"), massOmega);
      // fill filters and inv mass for no cascade selections (MC truth)
      if constexpr (requires { collision.straMCCollisionId(); }) {
        if (isMCTruth(stdCasc, "Xi") || isMCTruth(stdCasc, "Omega")) {
          histos.fill(HIST(TypeNames[Type]) + HIST("/NoSel-Truth/Filters/PropDCAxy"), cascade.dcaXYCascToPV());
          histos.fill(HIST(TypeNames[Type]) + HIST("/NoSel-Truth/Filters/PropDCAz"), cascade.dcaZCascToPV());
          ////* histos.fill(HIST(TypeNames[Type]) + HIST("/NoSel-Truth/Filters/CalcDCA"), calcDCA);
          histos.fill(HIST(TypeNames[Type]) + HIST("/NoSel-Truth/Filters/BachCosPA"), stdCasc.bachBaryonCosPA());
          histos.fill(HIST(TypeNames[Type]) + HIST("/NoSel-Truth/Filters/V0CosPA"), v0cosPA);
          histos.fill(HIST(TypeNames[Type]) + HIST("/NoSel-Truth/Filters/CascCosPA"), casccosPA);
          histos.fill(HIST(TypeNames[Type]) + HIST("/NoSel-Truth/Filters/EtaDau"), bachEta);
          histos.fill(HIST(TypeNames[Type]) + HIST("/NoSel-Truth/Filters/EtaDau"), negEta);
          histos.fill(HIST(TypeNames[Type]) + HIST("/NoSel-Truth/Filters/EtaDau"), posEta);
          if (isMCTruth(stdCasc, "Xi")) {
            histos.fill(HIST(TypeNames[Type]) + HIST("/NoSel-Truth/Filters/RapidityXi"), cascade.yXi());
            histos.fill(HIST(TypeNames[Type]) + HIST("/NoSel-Truth/MassXi"), massXi);
          }
          if (isMCTruth(stdCasc, "Omega")) {
            histos.fill(HIST(TypeNames[Type]) + HIST("/NoSel-Truth/Filters/RapidityOmega"), cascade.yOmega());
            histos.fill(HIST(TypeNames[Type]) + HIST("/NoSel-Truth/MassOmega"), massOmega);
          }
        }
      }
      // start checking selections
      bool passedAllSelsXi = true;
      bool passedAllSelsOmega = true;
      bool fillTruthXi = false;
      bool fillTruthOmega = false;
      if constexpr (requires { collision.straMCCollisionId(); }) {
        if (isMCTruth(stdCasc, "Xi")) {
          fillTruthXi = true;
        }
      }
      if constexpr (requires { collision.straMCCollisionId(); }) {
        if (isMCTruth(stdCasc, "Omega")) {
          fillTruthOmega = true;
        }
      }
      // apply pt cuts
      if (doApplyPtCutsXi) {
        if (isValidPt(cascade, "Xi", Type)) {
          histos.fill(HIST(TypeNames[Type]) + HIST("/Rec/FiltersXi"), 0.5);
          if (fillTruthXi)
            histos.fill(HIST(TypeNames[Type]) + HIST("/Rec-Truth/FiltersXi"), 0.5);
        } else {
          passedAllSelsXi = false;
        }
      }
      if (doApplyPtCutsOmega) {
        if (isValidPt(cascade, "Omega", Type)) {
          histos.fill(HIST(TypeNames[Type]) + HIST("/Rec/FiltersOmega"), 0.5);
          if (fillTruthOmega)
            histos.fill(HIST(TypeNames[Type]) + HIST("/Rec-Truth/FiltersOmega"), 0.5);
        } else {
          passedAllSelsOmega = false;
        }
      }
      // apply general cascade cuts
      if (doApplyGenCutsXi) {
        auto genSels = isValidCasc(collision, cascade, stdCasc, "Xi");
        for (size_t i = 0; i < std::size(genSels); ++i) {
          if (genSels[i]) {
            histos.fill(HIST(TypeNames[Type]) + HIST("/Rec/GenFiltersXi"), (i + 0.5));
            if (fillTruthXi)
              histos.fill(HIST(TypeNames[Type]) + HIST("/Rec-Truth/GenFiltersXi"), (i + 0.5));
          }
        }
        if (genSels[std::size(genSels) - 1]) {
          histos.fill(HIST(TypeNames[Type]) + HIST("/Rec/FiltersXi"), 1.5);
          if (fillTruthXi)
            histos.fill(HIST(TypeNames[Type]) + HIST("/Rec-Truth/FiltersXi"), 1.5);
        } else {
          passedAllSelsXi = false;
        }
      }
      if (doApplyGenCutsOmega) {
        auto genSels = isValidCasc(collision, cascade, stdCasc, "Omega");
        for (size_t i = 0; i < std::size(genSels); ++i) {
          if (genSels[i]) {
            histos.fill(HIST(TypeNames[Type]) + HIST("/Rec/GenFiltersOmega"), (i + 0.5));
            if (fillTruthOmega)
              histos.fill(HIST(TypeNames[Type]) + HIST("/Rec-Truth/GenFiltersOmega"), (i + 0.5));
          }
        }
        if (genSels[std::size(genSels) - 1]) {
          histos.fill(HIST(TypeNames[Type]) + HIST("/Rec/FiltersOmega"), 1.5);
          if (fillTruthOmega)
            histos.fill(HIST(TypeNames[Type]) + HIST("/Rec-Truth/FiltersOmega"), 1.5);
        } else {
          passedAllSelsOmega = false;
        }
      }
      // apply tpc pid
      if (doApplyTPCPIDXi) {
        if (passesTPC(stdCasc, "Xi")) {
          histos.fill(HIST(TypeNames[Type]) + HIST("/Rec/FiltersXi"), 2.5);
          if (fillTruthXi)
            histos.fill(HIST(TypeNames[Type]) + HIST("/Rec-Truth/FiltersXi"), 2.5);
        } else {
          passedAllSelsXi = false;
        }
      }
      if (doApplyTPCPIDOmega) {
        if (passesTPC(stdCasc, "Omega")) {
          histos.fill(HIST(TypeNames[Type]) + HIST("/Rec/FiltersOmega"), 2.5);
          if (fillTruthOmega)
            histos.fill(HIST(TypeNames[Type]) + HIST("/Rec-Truth/FiltersOmega"), 2.5);
        } else {
          passedAllSelsOmega = false;
        }
      }
      // apply tof pid
      if (doApplyTOFPIDXi) {
        if (passesTOF(stdCasc, "Xi")) {
          histos.fill(HIST(TypeNames[Type]) + HIST("/Rec/FiltersXi"), 3.5);
          if (fillTruthXi)
            histos.fill(HIST(TypeNames[Type]) + HIST("/Rec-Truth/FiltersXi"), 3.5);
        } else {
          passedAllSelsXi = false;
        }
      }
      if (doApplyTOFPIDOmega) {
        if (passesTOF(stdCasc, "Omega")) {
          histos.fill(HIST(TypeNames[Type]) + HIST("/Rec/FiltersOmega"), 3.5);
          if (fillTruthOmega)
            histos.fill(HIST(TypeNames[Type]) + HIST("/Rec-Truth/FiltersOmega"), 3.5);
        } else {
          passedAllSelsOmega = false;
        }
      }
      // apply competing mass rej
      if (doCompetingMassRej) {
        if ((std::abs(massXi - o2::constants::physics::MassXiMinus) > selCuts.cutCompMassRej)) {
          histos.fill(HIST(TypeNames[Type]) + HIST("/Rec/FiltersOmega"), 4.5);
          if (fillTruthOmega)
            histos.fill(HIST(TypeNames[Type]) + HIST("/Rec-Truth/FiltersOmega"), 4.5);
        } else {
          passedAllSelsOmega = false;
        }
      }

      // fil rec histograms
      if (passedAllSelsXi || passedAllSelsOmega) {
        histos.fill(HIST(TypeNames[Type]) + HIST("/Rec/EvMult"), mult);
        histos.fill(HIST(TypeNames[Type]) + HIST("/Rec/Filters/PropDCAxy"), cascade.dcaXYCascToPV());
        histos.fill(HIST(TypeNames[Type]) + HIST("/Rec/Filters/PropDCAz"), cascade.dcaZCascToPV());
        ////* histos.fill(HIST(TypeNames[Type]) + HIST("/Rec/Filters/CalcDCA"), calcDCA);
        histos.fill(HIST(TypeNames[Type]) + HIST("/Rec/Filters/BachCosPA"), stdCasc.bachBaryonCosPA());
        histos.fill(HIST(TypeNames[Type]) + HIST("/Rec/Filters/V0CosPA"), v0cosPA);
        histos.fill(HIST(TypeNames[Type]) + HIST("/Rec/Filters/CascCosPA"), casccosPA);
        histos.fill(HIST(TypeNames[Type]) + HIST("/Rec/Filters/EtaDau"), bachEta);
        histos.fill(HIST(TypeNames[Type]) + HIST("/Rec/Filters/EtaDau"), negEta);
        histos.fill(HIST(TypeNames[Type]) + HIST("/Rec/Filters/EtaDau"), posEta);
        if (passedAllSelsXi) {
          histos.fill(HIST(TypeNames[Type]) + HIST("/Rec/Filters/RapidityXi"), cascade.yXi());
          histos.fill(HIST(TypeNames[Type]) + HIST("/Rec/Filters/PtVsRapidityXi"), pt, cascade.yXi());
        }
        if (passedAllSelsOmega) {
          histos.fill(HIST(TypeNames[Type]) + HIST("/Rec/Filters/RapidityOmega"), cascade.yOmega());
          histos.fill(HIST(TypeNames[Type]) + HIST("/Rec/Filters/PtVsRapidityOmega"), pt, cascade.yOmega());
        }
        if (fillTruthXi || fillTruthOmega) {
          histos.fill(HIST(TypeNames[Type]) + HIST("/Rec-Truth/EvMult"), mult);
          histos.fill(HIST(TypeNames[Type]) + HIST("/Rec-Truth/Filters/PropDCAxy"), cascade.dcaXYCascToPV());
          histos.fill(HIST(TypeNames[Type]) + HIST("/Rec-Truth/Filters/PropDCAz"), cascade.dcaZCascToPV());
          ////* histos.fill(HIST(TypeNames[Type]) + HIST("/Rec-Truth/Filters/CalcDCA"), calcDCA);
          histos.fill(HIST(TypeNames[Type]) + HIST("/Rec-Truth/Filters/BachCosPA"), stdCasc.bachBaryonCosPA());
          histos.fill(HIST(TypeNames[Type]) + HIST("/Rec-Truth/Filters/V0CosPA"), v0cosPA);
          histos.fill(HIST(TypeNames[Type]) + HIST("/Rec-Truth/Filters/CascCosPA"), casccosPA);
          histos.fill(HIST(TypeNames[Type]) + HIST("/Rec-Truth/Filters/EtaDau"), bachEta);
          histos.fill(HIST(TypeNames[Type]) + HIST("/Rec-Truth/Filters/EtaDau"), negEta);
          histos.fill(HIST(TypeNames[Type]) + HIST("/Rec-Truth/Filters/EtaDau"), posEta);
          if (passedAllSelsXi)
            histos.fill(HIST(TypeNames[Type]) + HIST("/Rec-Truth/Filters/RapidityXi"), cascade.yXi());
          if (passedAllSelsOmega)
            histos.fill(HIST(TypeNames[Type]) + HIST("/Rec-Truth/Filters/RapidityOmega"), cascade.yOmega());
        }
      }
      double binFillXi[3] = {massXi, pt, mult};
      if (passedAllSelsXi) {
        histos.fill(HIST(TypeNames[Type]) + HIST("/Rec/FiltersXi"), 4.5);
        histos.fill(HIST(TypeNames[Type]) + HIST("/Rec/MassXi"), massXi);
        fillHist(histos.get<THn>(HIST(TypeNames[Type]) + HIST("/Rec/Xi")), binFillXi, efficiencyXi, efficiencyXiErr, purityXi, purityXiErr);
        // fill for particle-antiparticle type
        if (cascade.sign() < 0) {
          histos.fill(HIST(TypeNames[Type]) + HIST("/Rec/MassXiMinus"), massXi);
          fillHist(histos.get<THn>(HIST(TypeNames[Type]) + HIST("/Rec/XiMinus")), binFillXi, efficiencyXi, efficiencyXiErr, purityXi, purityXiErr);
        }
        if (cascade.sign() > 0) {
          histos.fill(HIST(TypeNames[Type]) + HIST("/Rec/MassXiPlus"), massXi);
          fillHist(histos.get<THn>(HIST(TypeNames[Type]) + HIST("/Rec/XiPlus")), binFillXi, efficiencyXi, efficiencyXiErr, purityXi, purityXiErr);
        }
        // fill truth
        if (fillTruthXi) {
          histos.fill(HIST(TypeNames[Type]) + HIST("/Rec-Truth/FiltersXi"), 4.5);
          histos.fill(HIST(TypeNames[Type]) + HIST("/Rec-Truth/MassXi"), massXi);
          histos.fill(HIST(TypeNames[Type]) + HIST("/Rec-Truth/Xi"), massXi, pt, mult);
        }
      }
      double binFillOmega[3] = {massOmega, pt, mult};
      if (passedAllSelsOmega) {
        histos.fill(HIST(TypeNames[Type]) + HIST("/Rec/FiltersOmega"), 5.5);
        histos.fill(HIST(TypeNames[Type]) + HIST("/Rec/MassOmega"), massOmega);
        fillHist(histos.get<THn>(HIST(TypeNames[Type]) + HIST("/Rec/Omega")), binFillOmega, efficiencyOmega, efficiencyOmegaErr, purityOmega, purityOmegaErr);
        // fill for particle-antiparticle type
        if (cascade.sign() < 0) {
          histos.fill(HIST(TypeNames[Type]) + HIST("/Rec/MassOmegaMinus"), massOmega);
          fillHist(histos.get<THn>(HIST(TypeNames[Type]) + HIST("/Rec/OmegaMinus")), binFillOmega, efficiencyOmega, efficiencyOmegaErr, purityOmega, purityOmegaErr);
        }
        if (cascade.sign() > 0) {
          histos.fill(HIST(TypeNames[Type]) + HIST("/Rec/MassOmegaPlus"), massOmega);
          fillHist(histos.get<THn>(HIST(TypeNames[Type]) + HIST("/Rec/OmegaPlus")), binFillOmega, efficiencyOmega, efficiencyOmegaErr, purityOmega, purityOmegaErr);
        }
        // fill truth
        if (fillTruthOmega) {
          histos.fill(HIST(TypeNames[Type]) + HIST("/Rec-Truth/FiltersOmega"), 5.5);
          histos.fill(HIST(TypeNames[Type]) + HIST("/Rec-Truth/MassOmega"), massOmega);
          histos.fill(HIST(TypeNames[Type]) + HIST("/Rec-Truth/Omega"), massOmega, pt, mult);
        }
      }

      // statistics - compare gen and reco pt and rapidity
      if constexpr (requires { collision.straMCCollisionId(); }) {
        if constexpr (requires { stdCasc.has_cascMCCore(); }) {
          auto cascmccore = stdCasc.template cascMCCore_as<DerMCGenCascades>();
          double genPt = cascmccore.ptMC();
          double genYXi = cascmccore.rapidityMC(0);
          double genYOmega = cascmccore.rapidityMC(2);
          histos.fill(HIST(TypeNames[Type]) + HIST("/NoSel/GenRecPt"), genPt, pt);
          histos.fill(HIST(TypeNames[Type]) + HIST("/NoSel/GenRecRapidityXi"), genYXi, cascade.yXi());
          histos.fill(HIST(TypeNames[Type]) + HIST("/NoSel/GenRecRapidityOmega"), genYOmega, cascade.yOmega());
          if (fillTruthXi || fillTruthOmega)
            histos.fill(HIST(TypeNames[Type]) + HIST("/NoSel-Truth/GenRecPt"), genPt, pt);
          if (passedAllSelsOmega || passedAllSelsXi) {
            histos.fill(HIST(TypeNames[Type]) + HIST("/Rec/GenRecPt"), genPt, pt);
              if (fillTruthXi || fillTruthOmega)
                histos.fill(HIST(TypeNames[Type]) + HIST("/Rec-Truth/GenRecPt"), genPt, pt);
          }
          if (fillTruthXi)
            histos.fill(HIST(TypeNames[Type]) + HIST("/NoSel-Truth/GenRecRapidityXi"), genYXi, cascade.yXi());
          if (passedAllSelsXi) {
            histos.fill(HIST(TypeNames[Type]) + HIST("/Rec/GenRecRapidityXi"), genYXi, cascade.yXi());
              if (fillTruthXi)
                histos.fill(HIST(TypeNames[Type]) + HIST("/Rec-Truth/GenRecRapidityXi"), genYXi, cascade.yXi());
          }
          if (fillTruthOmega)
            histos.fill(HIST(TypeNames[Type]) + HIST("/NoSel-Truth/GenRecRapidityOmega"), genYOmega, cascade.yOmega());
          if (passedAllSelsOmega) {
            histos.fill(HIST(TypeNames[Type]) + HIST("/Rec/GenRecRapidityOmega"), genYOmega, cascade.yOmega());
              if (fillTruthOmega)
                histos.fill(HIST(TypeNames[Type]) + HIST("/Rec-Truth/GenRecRapidityOmega"), genYOmega, cascade.yOmega());
          }
        }
      }
    }
  }

  void init(InitContext const&)
  {
    // check if cut configuration is valid
    if (selCuts.cutMinCascCosPaVsPt->size() != axesConfig.axisPt.value.size() - 1) {
    LOGF(fatal, "cutMinCascCosPaVsPt size does not match axisPt binning");
    }
    // for all events processing
    histos.add("NoSel-Events/EvCounter", "Event Counter", kTH1D, {{1, 0, 1}});
    histos.add("NoSel-Events/PVxy", "PV xy position", kTH2D, {{200, -0.1, 0.1}, {200, -0.1, 0.1}});
    histos.add("NoSel-Events/PVz", "PV z position", kTH1D, {{100, -20, 20}});
    histos.add("NoSel-Events/Mult", "Multiplicity", kTH1D, {axesConfig.axisMult});
    // for all events processing
    histos.add("Rec-Events/EvCounter", "Event Counter", kTH1D, {{1, 0, 1}});
    histos.add("Rec-Events/PVxy", "PV xy position", kTH2D, {{200, -0.1, 0.1}, {200, -0.1, 0.1}});
    histos.add("Rec-Events/PVz", "PV z position", kTH1D, {{100, -20, 20}});
    histos.add("Rec-Events/Mult", "Multiplicity", kTH1D, {axesConfig.axisMult});
    histos.add("Rec-Events/EvFilter", "Event Filter", kTH1D, {{6, 0, 6}});
    histos.get<TH1>(HIST("Rec-Events/EvFilter"))->GetXaxis()->SetBinLabel(1, "INEL>0");
    histos.get<TH1>(HIST("Rec-Events/EvFilter"))->GetXaxis()->SetBinLabel(2, "PVz cut");
    histos.get<TH1>(HIST("Rec-Events/EvFilter"))->GetXaxis()->SetBinLabel(3, "sel8");
    histos.get<TH1>(HIST("Rec-Events/EvFilter"))->GetXaxis()->SetBinLabel(4, "NoSameBunchPileup");
    histos.get<TH1>(HIST("Rec-Events/EvFilter"))->GetXaxis()->SetBinLabel(5, "IsGoodZvtxFT0VsPV");
    histos.get<TH1>(HIST("Rec-Events/EvFilter"))->GetXaxis()->SetBinLabel(6, "all");
    // for cascade processing
    static_for<0, 1>([&](auto Type) {
      // no selections applied
      histos.add(Form("%s/NoSel/Filters/PropDCAxy", TypeNames[Type].data()), "DCA to xy (propagated)", kTH1D, {axesConfig.axisDCAxy});
      histos.add(Form("%s/NoSel/Filters/PropDCAz", TypeNames[Type].data()), "DCA to z (propagated)", kTH1D, {axesConfig.axisDCAz});
      ////* histos.add(Form("%s/NoSel/Filters/CalcDCA", TypeNames[Type].data()), "DCA (calculated)", kTH1D, {axesConfig.axisDCAxy});
      histos.add(Form("%s/NoSel/Filters/BachCosPA", TypeNames[Type].data()), "Bachelor cosPA", kTH1D, {{200, -1.0, 1.0}});
      histos.add(Form("%s/NoSel/Filters/V0CosPA", TypeNames[Type].data()), "V0 cosPA", kTH1D, {{200, -1.0, 1.0}});
      histos.add(Form("%s/NoSel/Filters/CascCosPA", TypeNames[Type].data()), "V0 cosPA", kTH1D, {{200, -1.0, 1.0}});
      histos.add(Form("%s/NoSel/Filters/RapidityXi", TypeNames[Type].data()), "y under Xi hypothesis", kTH1D, {{200, -1.0, 1.0}});
      histos.add(Form("%s/NoSel/Filters/PtVsRapidityXi", TypeNames[Type].data()), "pt vs y under Xi hypothesis", kTH2D, {axesConfig.axisPt, {200, -1.0, 1.0}});
      histos.add(Form("%s/NoSel/Filters/RapidityOmega", TypeNames[Type].data()), "y under Omega hypothesis", kTH1D, {{200, -1.0, 1.0}});
      histos.add(Form("%s/NoSel/Filters/PtVsRapidityOmega", TypeNames[Type].data()), "pt vs y under Xi hypothesis", kTH2D, {axesConfig.axisPt, {200, -1.0, 1.0}});
      histos.add(Form("%s/NoSel/Filters/EtaDau", TypeNames[Type].data()), "|#eta| of dau tracks", kTH1D, {axesConfig.axisEta});
      histos.add(Form("%s/NoSel/EvMult", TypeNames[Type].data()), "Multiplicity of events with >=1 cascade", kTH1D, {axesConfig.axisMult});
      histos.add(Form("%s/NoSel/GenRecPt", TypeNames[Type].data()), "Generated vs reconstructed pt", kTH2D, {axesConfig.axisPt, axesConfig.axisPt});
      histos.add(Form("%s/NoSel/GenRecRapidityXi", TypeNames[Type].data()), "Generated vs reconstructed y (Xi)", kTH2D, {{200, -1.0, 1.0}, {200, -1.0, 1.0}});
      histos.add(Form("%s/NoSel/GenRecRapidityOmega", TypeNames[Type].data()), "Generated vs reconstructed y (Omega)", kTH2D, {{200, -1.0, 1.0}, {200, -1.0, 1.0}});
      histos.add(Form("%s/NoSel/MassXi", TypeNames[Type].data()), "Invariant mass hypothesis", kTH1D, {axesConfig.axisXiMass});
      histos.add(Form("%s/NoSel/MassOmega", TypeNames[Type].data()), "Invariant mass hypothesis", kTH1D, {axesConfig.axisOmegaMass});
      // mc truth for no selectrion
      histos.add(Form("%s/NoSel-Truth/Filters/PropDCAxy", TypeNames[Type].data()), "DCA to xy (propagated)", kTH1D, {axesConfig.axisDCAxy});
      histos.add(Form("%s/NoSel-Truth/Filters/PropDCAz", TypeNames[Type].data()), "DCA to z (propagated)", kTH1D, {axesConfig.axisDCAz});
      ////* histos.add(Form("%s/NoSel-Truth/Filters/CalcDCA", TypeNames[Type].data()), "DCA (calculated)", kTH1D, {axesConfig.axisDCAxy});
      histos.add(Form("%s/NoSel-Truth/Filters/BachCosPA", TypeNames[Type].data()), "Bachelor cosPA", kTH1D, {{200, -1.0, 1.0}});
      histos.add(Form("%s/NoSel-Truth/Filters/V0CosPA", TypeNames[Type].data()), "V0 cosPA", kTH1D, {{200, -1.0, 1.0}});
      histos.add(Form("%s/NoSel-Truth/Filters/CascCosPA", TypeNames[Type].data()), "V0 cosPA", kTH1D, {{200, -1.0, 1.0}});
      histos.add(Form("%s/NoSel-Truth/Filters/RapidityXi", TypeNames[Type].data()), "y under Xi hypothesis", kTH1D, {{200, -1.0, 1.0}});
      histos.add(Form("%s/NoSel-Truth/Filters/RapidityOmega", TypeNames[Type].data()), "y under Omega hypothesis", kTH1D, {{200, -1.0, 1.0}});
      histos.add(Form("%s/NoSel-Truth/Filters/EtaDau", TypeNames[Type].data()), "|#eta| of dau tracks", kTH1D, {axesConfig.axisEta});
      histos.add(Form("%s/NoSel-Truth/EvMult", TypeNames[Type].data()), "Multiplicity of events with >=1 cascade", kTH1D, {axesConfig.axisMult});
      histos.add(Form("%s/NoSel-Truth/GenRecPt", TypeNames[Type].data()), "Generated vs reconstructed pt", kTH2D, {axesConfig.axisPt, axesConfig.axisPt});
      histos.add(Form("%s/NoSel-Truth/GenRecRapidityXi", TypeNames[Type].data()), "Generated vs reconstructed y (Xi)", kTH2D, {{200, -1.0, 1.0}, {200, -1.0, 1.0}});
      histos.add(Form("%s/NoSel-Truth/GenRecRapidityOmega", TypeNames[Type].data()), "Generated vs reconstructed y (Omega)", kTH2D, {{200, -1.0, 1.0}, {200, -1.0, 1.0}});
      histos.add(Form("%s/NoSel-Truth/MassXi", TypeNames[Type].data()), "Invariant mass hypothesis", kTH1D, {axesConfig.axisXiMass});
      histos.add(Form("%s/NoSel-Truth/MassOmega", TypeNames[Type].data()), "Invariant mass hypothesis", kTH1D, {axesConfig.axisOmegaMass});
      // xi and omega selection statistics
      histos.add(Form("%s/Rec/Filters/PropDCAxy", TypeNames[Type].data()), "DCA to xy (propagated)", kTH1D, {axesConfig.axisDCAxy});
      histos.add(Form("%s/Rec/Filters/PropDCAz", TypeNames[Type].data()), "DCA to z (propagated)", kTH1D, {axesConfig.axisDCAz});
      ////* histos.add(Form("%s/Rec/Filters/CalcDCA", TypeNames[Type].data()), "DCA (calculated)", kTH1D, {axesConfig.axisDCAxy});
      histos.add(Form("%s/Rec/Filters/BachCosPA", TypeNames[Type].data()), "Bachelor cosPA", kTH1D, {{200, -1.0, 1.0}});
      histos.add(Form("%s/Rec/Filters/V0CosPA", TypeNames[Type].data()), "V0 cosPA", kTH1D, {{200, -1.0, 1.0}});
      histos.add(Form("%s/Rec/Filters/CascCosPA", TypeNames[Type].data()), "V0 cosPA", kTH1D, {{200, -1.0, 1.0}});
      histos.add(Form("%s/Rec/Filters/RapidityXi", TypeNames[Type].data()), "y under Xi hypothesis", kTH1D, {{200, -1.0, 1.0}});
      histos.add(Form("%s/Rec/Filters/PtVsRapidityXi", TypeNames[Type].data()), "pt vs y under Xi hypothesis", kTH2D, {axesConfig.axisPt, {200, -1.0, 1.0}});
      histos.add(Form("%s/Rec/Filters/RapidityOmega", TypeNames[Type].data()), "y under Omega hypothesis", kTH1D, {{200, -1.0, 1.0}});
      histos.add(Form("%s/Rec/Filters/PtVsRapidityOmega", TypeNames[Type].data()), "pt vs y under Omega hypothesis", kTH2D, {axesConfig.axisPt, {200, -1.0, 1.0}});
      histos.add(Form("%s/Rec/Filters/EtaDau", TypeNames[Type].data()), "|#eta| of dau tracks", kTH1D, {axesConfig.axisEta});
      histos.add(Form("%s/Rec/EvMult", TypeNames[Type].data()), "Multiplicity of events with >=1 selected cascade", kTH1D, {axesConfig.axisMult});
      histos.add(Form("%s/Rec/GenRecPt", TypeNames[Type].data()), "Generated vs reconstructed pt", kTH2D, {axesConfig.axisPt, axesConfig.axisPt});
      histos.add(Form("%s/Rec/GenRecRapidityXi", TypeNames[Type].data()), "Generated vs reconstructed y (Xi)", kTH2D, {{200, -1.0, 1.0}, {200, -1.0, 1.0}});
      histos.add(Form("%s/Rec/GenRecRapidityOmega", TypeNames[Type].data()), "Generated vs reconstructed y (Omega)", kTH2D, {{200, -1.0, 1.0}, {200, -1.0, 1.0}});
      // selected xi
      histos.add(Form("%s/Rec/FiltersXi", TypeNames[Type].data()), "main cascade filters for Xi", kTH1D, {{5, 0, 5}});
      histos.add(Form("%s/Rec/GenFiltersXi", TypeNames[Type].data()), "general cascade filters for Xi", kTH1D, {{9, 0, 9}});
      histos.add(Form("%s/Rec/MassXi", TypeNames[Type].data()), "Invariant mass hypothesis", kTH1D, {axesConfig.axisXiMass});
      histos.add(Form("%s/Rec/MassXiMinus", TypeNames[Type].data()), "Invariant mass hypothesis", kTH1D, {axesConfig.axisXiMass});
      histos.add(Form("%s/Rec/MassXiPlus", TypeNames[Type].data()), "Invariant mass hypothesis", kTH1D, {axesConfig.axisXiMass});
      histos.add(Form("%s/Rec/Xi", TypeNames[Type].data()), "", kTHnD, {axesConfig.axisXiMass, axesConfig.axisPt, axesConfig.axisMult});
      histos.add(Form("%s/Rec/XiMinus", TypeNames[Type].data()), "", kTHnD, {axesConfig.axisXiMass, axesConfig.axisPt, axesConfig.axisMult});
      histos.add(Form("%s/Rec/XiPlus", TypeNames[Type].data()), "", kTHnD, {axesConfig.axisXiMass, axesConfig.axisPt, axesConfig.axisMult});
      // selected omega
      histos.add(Form("%s/Rec/FiltersOmega", TypeNames[Type].data()), "main cascade filters for Omega", kTH1D, {{6, 0, 6}});
      histos.add(Form("%s/Rec/GenFiltersOmega", TypeNames[Type].data()), "general cascade filters for Omega", kTH1D, {{9, 0, 9}});
      histos.add(Form("%s/Rec/MassOmega", TypeNames[Type].data()), "Invariant mass hypothesis", kTH1D, {axesConfig.axisOmegaMass});
      histos.add(Form("%s/Rec/MassOmegaMinus", TypeNames[Type].data()), "Invariant mass hypothesis", kTH1D, {axesConfig.axisOmegaMass});
      histos.add(Form("%s/Rec/MassOmegaPlus", TypeNames[Type].data()), "Invariant mass hypothesis", kTH1D, {axesConfig.axisOmegaMass});
      histos.add(Form("%s/Rec/Omega", TypeNames[Type].data()), "", kTHnD, {axesConfig.axisOmegaMass, axesConfig.axisPt, axesConfig.axisMult});
      histos.add(Form("%s/Rec/OmegaMinus", TypeNames[Type].data()), "", kTHnD, {axesConfig.axisOmegaMass, axesConfig.axisPt, axesConfig.axisMult});
      histos.add(Form("%s/Rec/OmegaPlus", TypeNames[Type].data()), "", kTHnD, {axesConfig.axisOmegaMass, axesConfig.axisPt, axesConfig.axisMult});
      // mc truth for all passed selections
      // xi and omega truth selection statistics
      histos.add(Form("%s/Rec-Truth/Filters/PropDCAxy", TypeNames[Type].data()), "DCA to xy (propagated)", kTH1D, {axesConfig.axisDCAxy});
      histos.add(Form("%s/Rec-Truth/Filters/PropDCAz", TypeNames[Type].data()), "DCA to z (propagated)", kTH1D, {axesConfig.axisDCAz});
      ////* histos.add(Form("%s/Rec-Truth/Filters/CalcDCA", TypeNames[Type].data()), "DCA (calculated)", kTH1D, {axesConfig.axisDCAxy});
      histos.add(Form("%s/Rec-Truth/Filters/BachCosPA", TypeNames[Type].data()), "Bachelor cosPA", kTH1D, {{200, -1.0, 1.0}});
      histos.add(Form("%s/Rec-Truth/Filters/V0CosPA", TypeNames[Type].data()), "V0 cosPA", kTH1D, {{200, -1.0, 1.0}});
      histos.add(Form("%s/Rec-Truth/Filters/CascCosPA", TypeNames[Type].data()), "V0 cosPA", kTH1D, {{200, -1.0, 1.0}});
      histos.add(Form("%s/Rec-Truth/Filters/RapidityXi", TypeNames[Type].data()), "y under Xi hypothesis", kTH1D, {{200, -1.0, 1.0}});
      histos.add(Form("%s/Rec-Truth/Filters/RapidityOmega", TypeNames[Type].data()), "y under Omega hypothesis", kTH1D, {{200, -1.0, 1.0}});
      histos.add(Form("%s/Rec-Truth/Filters/EtaDau", TypeNames[Type].data()), "|#eta| of dau tracks", kTH1D, {axesConfig.axisEta});
      histos.add(Form("%s/Rec-Truth/EvMult", TypeNames[Type].data()), "Multiplicity of events with >=1 cascade", kTH1D, {axesConfig.axisMult});
      histos.add(Form("%s/Rec-Truth/GenRecPt", TypeNames[Type].data()), "Generated vs reconstructed pt", kTH2D, {axesConfig.axisPt, axesConfig.axisPt});
      histos.add(Form("%s/Rec-Truth/GenRecRapidityXi", TypeNames[Type].data()), "Generated vs reconstructed y (Xi)", kTH2D, {{200, -1.0, 1.0}, {200, -1.0, 1.0}});
      histos.add(Form("%s/Rec-Truth/GenRecRapidityOmega", TypeNames[Type].data()), "Generated vs reconstructed y (Omega)", kTH2D, {{200, -1.0, 1.0}, {200, -1.0, 1.0}});
      histos.add(Form("%s/Rec-Truth/FiltersXi", TypeNames[Type].data()), "main cascade filters for Xi", kTH1D, {{5, 0, 5}});
      histos.add(Form("%s/Rec-Truth/GenFiltersXi", TypeNames[Type].data()), "general cascade filters for Xi", kTH1D, {{9, 0, 9}});
      histos.add(Form("%s/Rec-Truth/MassXi", TypeNames[Type].data()), "Invariant mass hypothesis", kTH1D, {axesConfig.axisXiMass});
      histos.add(Form("%s/Rec-Truth/Xi", TypeNames[Type].data()), "", kTHnD, {axesConfig.axisXiMass, axesConfig.axisPt, axesConfig.axisMult});
      histos.add(Form("%s/Rec-Truth/FiltersOmega", TypeNames[Type].data()), "main cascade filters for Omega", kTH1D, {{6, 0, 6}});
      histos.add(Form("%s/Rec-Truth/GenFiltersOmega", TypeNames[Type].data()), "general cascade filters for Omega", kTH1D, {{9, 0, 9}});
      histos.add(Form("%s/Rec-Truth/MassOmega", TypeNames[Type].data()), "Invariant mass hypothesis", kTH1D, {axesConfig.axisOmegaMass});
      histos.add(Form("%s/Rec-Truth/Omega", TypeNames[Type].data()), "", kTHnD, {axesConfig.axisOmegaMass, axesConfig.axisPt, axesConfig.axisMult});
      
    });
    // for MC-specific processing
    // all generated events:
    histos.add("MC/Gen/EvCounter", "Event Counter", kTH1D, {{1, 0, 1}});
    histos.add("MC/Gen/Mult", "Event charged multiplicty (all gen events)", kTH1D, {{3200, 0, 3200}});
    histos.add("MC/Gen/RapidityXi", "Generated y (Xi)", kTH1D, {{200, -1.0, 1.0}});
    histos.add("MC/Gen/Xi", "Xi", kTH2D, {{3200, 0, 3200}, {100, 0.0, 10.0}});                                             // generated Xis
    histos.add("MC/Gen/PrimaryXi", "Xi primaries", kTH2D, {{3200, 0, 3200}, {100, 0.0, 10.0}});                            // generated primary Xis
    histos.add("MC/Gen/PrimaryXiRapidity", "Xi primaries in |y|", kTH2D, {{3200, 0, 3200}, {100, 0.0, 10.0}});             // generated primary Xis in selected rapidity range
    histos.add("MC/Gen/PrimaryXiMinusRapidity", "Xi- primaries in |y| ", kTH2D, {{3200, 0, 3200}, {100, 0.0, 10.0}});      // generated primary Xi-minus in selected rapidity range
    histos.add("MC/Gen/PrimaryXiPlusRapidity", "Xi+ primaries in |y|", kTH2D, {{3200, 0, 3200}, {100, 0.0, 10.0}});        // generated primary Xi-plus in selected rapidity range
    histos.add("MC/Gen/RapidityOmega", "Generated y (Omega)", kTH1D, {{200, -1.0, 1.0}});
    histos.add("MC/Gen/Omega", "Omega", kTH2D, {{3200, 0, 3200}, {100, 0.0, 10.0}});                                             // generated Omegas
    histos.add("MC/Gen/PrimaryOmega", "Omega primaries", kTH2D, {{3200, 0, 3200}, {100, 0.0, 10.0}});                            // generated primary Omegas
    histos.add("MC/Gen/PrimaryOmegaRapidity", "Omega primaries in |y|", kTH2D, {{3200, 0, 3200}, {100, 0.0, 10.0}});             // generated primary Omegas in selected rapidity range
    histos.add("MC/Gen/PrimaryOmegaMinusRapidity", "Omega- primaries in |y| ", kTH2D, {{3200, 0, 3200}, {100, 0.0, 10.0}});      // generated primary Omega-minus in selected rapidity range
    histos.add("MC/Gen/PrimaryOmegaPlusRapidity", "Omega+ primaries in |y|", kTH2D, {{3200, 0, 3200}, {100, 0.0, 10.0}});        // generated primary Omega-plus in selected rapidity range
    // generated events with reco >= 1:
    histos.add("MC/EvRec/EvCounter", "Event Counter", kTH1D, {{1, 0, 1}});
    histos.add("MC/EvRec/Mult", "Event charged multiplicty (gen events with reco>=1)", kTH1D, {{3200, 0, 3200}});
    histos.add("MC/EvRec/MultCent", "Gen multiplicity vs reco centrality", kTH2D, {{100, 0, 100}, {3200, 0, 3200}});
    histos.add("MC/EvRec/Xi", "Xi", kTH2D, {axesConfig.axisPt, axesConfig.axisMult});                                             // generated Xis in reco>=1
    histos.add("MC/EvRec/PrimaryXi", "Xi primaries", kTH2D, {axesConfig.axisPt, axesConfig.axisMult});                            // generated primary Xis in reco>=1
    histos.add("MC/EvRec/PrimaryXiRapidity", "Xi primaries in |y|", kTH2D, {axesConfig.axisPt, axesConfig.axisMult});             // generated primary Xis in selected rapidity range in reco>=1
    histos.add("MC/EvRec/PrimaryXiMinusRapidity", "Xi- primaries in |y| ", kTH2D, {axesConfig.axisPt, axesConfig.axisMult});      // generated primary Xi-minus in selected rapidity range in reco>=1
    histos.add("MC/EvRec/PrimaryXiPlusRapidity", "Xi+ primaries in |y|", kTH2D, {axesConfig.axisPt, axesConfig.axisMult});        // generated primary Xi-plus in selected rapidity range in reco>=1
    histos.add("MC/EvRec/Omega", "Omega", kTH2D, {axesConfig.axisPt, axesConfig.axisMult});                                             // generated Omegas in reco>=1
    histos.add("MC/EvRec/PrimaryOmega", "Omega primaries", kTH2D, {axesConfig.axisPt, axesConfig.axisMult});                            // generated primary Omegas in reco>=1
    histos.add("MC/EvRec/PrimaryOmegaRapidity", "Omega primaries in |y|", kTH2D, {axesConfig.axisPt, axesConfig.axisMult});             // generated primary Omegas in selected rapidity range in reco>=1
    histos.add("MC/EvRec/PrimaryOmegaMinusRapidity", "Omega- primaries in |y|", kTH2D, {axesConfig.axisPt, axesConfig.axisMult});      // generated primary Omega-minus in selected rapidity range in reco>=1
    histos.add("MC/EvRec/PrimaryOmegaPlusRapidity", "Omega+ primaries in |y|", kTH2D, {axesConfig.axisPt, axesConfig.axisMult});        // generated primary Omega-plus in selected rapidity range in reco>=1
    // label filter statistic bins for standard cascs
    histos.get<TH1>(HIST("Standard/Rec/FiltersXi"))->GetXaxis()->SetBinLabel(1, "p_{T}");
    histos.get<TH1>(HIST("Standard/Rec/FiltersXi"))->GetXaxis()->SetBinLabel(2, "gen");
    histos.get<TH1>(HIST("Standard/Rec/FiltersXi"))->GetXaxis()->SetBinLabel(3, "TPC");
    histos.get<TH1>(HIST("Standard/Rec/FiltersXi"))->GetXaxis()->SetBinLabel(4, "TOF");
    histos.get<TH1>(HIST("Standard/Rec/FiltersXi"))->GetXaxis()->SetBinLabel(5, "all");
    histos.get<TH1>(HIST("Standard/Rec/GenFiltersXi"))->GetXaxis()->SetBinLabel(1, "casc |y|");
    histos.get<TH1>(HIST("Standard/Rec/GenFiltersXi"))->GetXaxis()->SetBinLabel(2, "dau |#eta|");
    histos.get<TH1>(HIST("Standard/Rec/GenFiltersXi"))->GetXaxis()->SetBinLabel(3, "dau TPC cls");
    histos.get<TH1>(HIST("Standard/Rec/GenFiltersXi"))->GetXaxis()->SetBinLabel(4, "V0CosPA");
    histos.get<TH1>(HIST("Standard/Rec/GenFiltersXi"))->GetXaxis()->SetBinLabel(5, "BachCosPA");
    histos.get<TH1>(HIST("Standard/Rec/GenFiltersXi"))->GetXaxis()->SetBinLabel(6, "CascCosPA");
    histos.get<TH1>(HIST("Standard/Rec/GenFiltersXi"))->GetXaxis()->SetBinLabel(7, "PropDCAxy");
    histos.get<TH1>(HIST("Standard/Rec/GenFiltersXi"))->GetXaxis()->SetBinLabel(8, "PropDCAz");
    histos.get<TH1>(HIST("Standard/Rec/GenFiltersXi"))->GetXaxis()->SetBinLabel(9, "all");
    histos.get<TH1>(HIST("Standard/Rec/FiltersOmega"))->GetXaxis()->SetBinLabel(1, "p_{T}");
    histos.get<TH1>(HIST("Standard/Rec/FiltersOmega"))->GetXaxis()->SetBinLabel(2, "gen");
    histos.get<TH1>(HIST("Standard/Rec/FiltersOmega"))->GetXaxis()->SetBinLabel(3, "TPC");
    histos.get<TH1>(HIST("Standard/Rec/FiltersOmega"))->GetXaxis()->SetBinLabel(4, "TOF");
    histos.get<TH1>(HIST("Standard/Rec/FiltersOmega"))->GetXaxis()->SetBinLabel(5, "CMR");
    histos.get<TH1>(HIST("Standard/Rec/FiltersOmega"))->GetXaxis()->SetBinLabel(6, "all");
    histos.get<TH1>(HIST("Standard/Rec/GenFiltersOmega"))->GetXaxis()->SetBinLabel(1, "casc |y|");
    histos.get<TH1>(HIST("Standard/Rec/GenFiltersOmega"))->GetXaxis()->SetBinLabel(2, "dau |#eta|");
    histos.get<TH1>(HIST("Standard/Rec/GenFiltersOmega"))->GetXaxis()->SetBinLabel(3, "dau TPC cls");
    histos.get<TH1>(HIST("Standard/Rec/GenFiltersOmega"))->GetXaxis()->SetBinLabel(4, "V0CosPA");
    histos.get<TH1>(HIST("Standard/Rec/GenFiltersOmega"))->GetXaxis()->SetBinLabel(5, "BachCosPA");
    histos.get<TH1>(HIST("Standard/Rec/GenFiltersOmega"))->GetXaxis()->SetBinLabel(6, "CascCosPA");
    histos.get<TH1>(HIST("Standard/Rec/GenFiltersOmega"))->GetXaxis()->SetBinLabel(7, "PropDCAxy");
    histos.get<TH1>(HIST("Standard/Rec/GenFiltersOmega"))->GetXaxis()->SetBinLabel(8, "PropDCAz");
    histos.get<TH1>(HIST("Standard/Rec/GenFiltersOmega"))->GetXaxis()->SetBinLabel(9, "all");
    histos.get<TH1>(HIST("Standard/Rec-Truth/FiltersXi"))->GetXaxis()->SetBinLabel(1, "p_{T}");
    histos.get<TH1>(HIST("Standard/Rec-Truth/FiltersXi"))->GetXaxis()->SetBinLabel(2, "gen");
    histos.get<TH1>(HIST("Standard/Rec-Truth/FiltersXi"))->GetXaxis()->SetBinLabel(3, "TPC");
    histos.get<TH1>(HIST("Standard/Rec-Truth/FiltersXi"))->GetXaxis()->SetBinLabel(4, "TOF");
    histos.get<TH1>(HIST("Standard/Rec-Truth/FiltersXi"))->GetXaxis()->SetBinLabel(5, "all");
    histos.get<TH1>(HIST("Standard/Rec-Truth/GenFiltersXi"))->GetXaxis()->SetBinLabel(1, "casc |y|");
    histos.get<TH1>(HIST("Standard/Rec-Truth/GenFiltersXi"))->GetXaxis()->SetBinLabel(2, "dau |#eta|");
    histos.get<TH1>(HIST("Standard/Rec-Truth/GenFiltersXi"))->GetXaxis()->SetBinLabel(3, "dau TPC cls");
    histos.get<TH1>(HIST("Standard/Rec-Truth/GenFiltersXi"))->GetXaxis()->SetBinLabel(4, "V0CosPA");
    histos.get<TH1>(HIST("Standard/Rec-Truth/GenFiltersXi"))->GetXaxis()->SetBinLabel(5, "BachCosPA");
    histos.get<TH1>(HIST("Standard/Rec-Truth/GenFiltersXi"))->GetXaxis()->SetBinLabel(6, "CascCosPA");
    histos.get<TH1>(HIST("Standard/Rec-Truth/GenFiltersXi"))->GetXaxis()->SetBinLabel(7, "PropDCAxy");
    histos.get<TH1>(HIST("Standard/Rec-Truth/GenFiltersXi"))->GetXaxis()->SetBinLabel(8, "PropDCAz");
    histos.get<TH1>(HIST("Standard/Rec-Truth/GenFiltersXi"))->GetXaxis()->SetBinLabel(9, "all");
    histos.get<TH1>(HIST("Standard/Rec-Truth/FiltersOmega"))->GetXaxis()->SetBinLabel(1, "p_{T}");
    histos.get<TH1>(HIST("Standard/Rec-Truth/FiltersOmega"))->GetXaxis()->SetBinLabel(2, "gen");
    histos.get<TH1>(HIST("Standard/Rec-Truth/FiltersOmega"))->GetXaxis()->SetBinLabel(3, "TPC");
    histos.get<TH1>(HIST("Standard/Rec-Truth/FiltersOmega"))->GetXaxis()->SetBinLabel(4, "TOF");
    histos.get<TH1>(HIST("Standard/Rec-Truth/FiltersOmega"))->GetXaxis()->SetBinLabel(5, "CMR");
    histos.get<TH1>(HIST("Standard/Rec-Truth/FiltersOmega"))->GetXaxis()->SetBinLabel(6, "all");
    histos.get<TH1>(HIST("Standard/Rec-Truth/GenFiltersOmega"))->GetXaxis()->SetBinLabel(1, "casc |y|");
    histos.get<TH1>(HIST("Standard/Rec-Truth/GenFiltersOmega"))->GetXaxis()->SetBinLabel(2, "dau |#eta|");
    histos.get<TH1>(HIST("Standard/Rec-Truth/GenFiltersOmega"))->GetXaxis()->SetBinLabel(3, "dau TPC cls");
    histos.get<TH1>(HIST("Standard/Rec-Truth/GenFiltersOmega"))->GetXaxis()->SetBinLabel(4, "V0CosPA");
    histos.get<TH1>(HIST("Standard/Rec-Truth/GenFiltersOmega"))->GetXaxis()->SetBinLabel(5, "BachCosPA");
    histos.get<TH1>(HIST("Standard/Rec-Truth/GenFiltersOmega"))->GetXaxis()->SetBinLabel(6, "CascCosPA");
    histos.get<TH1>(HIST("Standard/Rec-Truth/GenFiltersOmega"))->GetXaxis()->SetBinLabel(7, "PropDCAxy");
    histos.get<TH1>(HIST("Standard/Rec-Truth/GenFiltersOmega"))->GetXaxis()->SetBinLabel(8, "PropDCAz");
    histos.get<TH1>(HIST("Standard/Rec-Truth/GenFiltersOmega"))->GetXaxis()->SetBinLabel(9, "all");
    // label filter statistic bins for tracked cascs
    histos.get<TH1>(HIST("Tracked/Rec/FiltersXi"))->GetXaxis()->SetBinLabel(1, "p_{T}");
    histos.get<TH1>(HIST("Tracked/Rec/FiltersXi"))->GetXaxis()->SetBinLabel(2, "gen");
    histos.get<TH1>(HIST("Tracked/Rec/FiltersXi"))->GetXaxis()->SetBinLabel(3, "TPC");
    histos.get<TH1>(HIST("Tracked/Rec/FiltersXi"))->GetXaxis()->SetBinLabel(4, "TOF");
    histos.get<TH1>(HIST("Tracked/Rec/FiltersXi"))->GetXaxis()->SetBinLabel(5, "all");
    histos.get<TH1>(HIST("Tracked/Rec/GenFiltersXi"))->GetXaxis()->SetBinLabel(1, "casc |y|");
    histos.get<TH1>(HIST("Tracked/Rec/GenFiltersXi"))->GetXaxis()->SetBinLabel(2, "dau |#eta|");
    histos.get<TH1>(HIST("Tracked/Rec/GenFiltersXi"))->GetXaxis()->SetBinLabel(3, "dau TPC cls");
    histos.get<TH1>(HIST("Tracked/Rec/GenFiltersXi"))->GetXaxis()->SetBinLabel(4, "V0CosPA");
    histos.get<TH1>(HIST("Tracked/Rec/GenFiltersXi"))->GetXaxis()->SetBinLabel(5, "BachCosPA");
    histos.get<TH1>(HIST("Tracked/Rec/GenFiltersXi"))->GetXaxis()->SetBinLabel(6, "CascCosPA");
    histos.get<TH1>(HIST("Tracked/Rec/GenFiltersXi"))->GetXaxis()->SetBinLabel(7, "PropDCAxy");
    histos.get<TH1>(HIST("Tracked/Rec/GenFiltersXi"))->GetXaxis()->SetBinLabel(8, "PropDCAz");
    histos.get<TH1>(HIST("Tracked/Rec/GenFiltersXi"))->GetXaxis()->SetBinLabel(9, "all");
    histos.get<TH1>(HIST("Tracked/Rec/FiltersOmega"))->GetXaxis()->SetBinLabel(1, "p_{T}");
    histos.get<TH1>(HIST("Tracked/Rec/FiltersOmega"))->GetXaxis()->SetBinLabel(2, "gen");
    histos.get<TH1>(HIST("Tracked/Rec/FiltersOmega"))->GetXaxis()->SetBinLabel(3, "TPC");
    histos.get<TH1>(HIST("Tracked/Rec/FiltersOmega"))->GetXaxis()->SetBinLabel(4, "TOF");
    histos.get<TH1>(HIST("Tracked/Rec/FiltersOmega"))->GetXaxis()->SetBinLabel(5, "CMR");
    histos.get<TH1>(HIST("Tracked/Rec/FiltersOmega"))->GetXaxis()->SetBinLabel(6, "all");
    histos.get<TH1>(HIST("Tracked/Rec/GenFiltersOmega"))->GetXaxis()->SetBinLabel(1, "casc |y|");
    histos.get<TH1>(HIST("Tracked/Rec/GenFiltersOmega"))->GetXaxis()->SetBinLabel(2, "dau |#eta|");
    histos.get<TH1>(HIST("Tracked/Rec/GenFiltersOmega"))->GetXaxis()->SetBinLabel(3, "dau TPC cls");
    histos.get<TH1>(HIST("Tracked/Rec/GenFiltersOmega"))->GetXaxis()->SetBinLabel(4, "V0CosPA");
    histos.get<TH1>(HIST("Tracked/Rec/GenFiltersOmega"))->GetXaxis()->SetBinLabel(5, "BachCosPA");
    histos.get<TH1>(HIST("Tracked/Rec/GenFiltersOmega"))->GetXaxis()->SetBinLabel(6, "CascCosPA");
    histos.get<TH1>(HIST("Tracked/Rec/GenFiltersOmega"))->GetXaxis()->SetBinLabel(7, "PropDCAxy");
    histos.get<TH1>(HIST("Tracked/Rec/GenFiltersOmega"))->GetXaxis()->SetBinLabel(8, "PropDCAz");
    histos.get<TH1>(HIST("Tracked/Rec/GenFiltersOmega"))->GetXaxis()->SetBinLabel(9, "all");
    histos.get<TH1>(HIST("Tracked/Rec-Truth/FiltersXi"))->GetXaxis()->SetBinLabel(1, "p_{T}");
    histos.get<TH1>(HIST("Tracked/Rec-Truth/FiltersXi"))->GetXaxis()->SetBinLabel(2, "gen");
    histos.get<TH1>(HIST("Tracked/Rec-Truth/FiltersXi"))->GetXaxis()->SetBinLabel(3, "TPC");
    histos.get<TH1>(HIST("Tracked/Rec-Truth/FiltersXi"))->GetXaxis()->SetBinLabel(4, "TOF");
    histos.get<TH1>(HIST("Tracked/Rec-Truth/FiltersXi"))->GetXaxis()->SetBinLabel(5, "all");
    histos.get<TH1>(HIST("Tracked/Rec-Truth/GenFiltersXi"))->GetXaxis()->SetBinLabel(1, "casc |y|");
    histos.get<TH1>(HIST("Tracked/Rec-Truth/GenFiltersXi"))->GetXaxis()->SetBinLabel(2, "dau |#eta|");
    histos.get<TH1>(HIST("Tracked/Rec-Truth/GenFiltersXi"))->GetXaxis()->SetBinLabel(3, "dau TPC cls");
    histos.get<TH1>(HIST("Tracked/Rec-Truth/GenFiltersXi"))->GetXaxis()->SetBinLabel(4, "V0CosPA");
    histos.get<TH1>(HIST("Tracked/Rec-Truth/GenFiltersXi"))->GetXaxis()->SetBinLabel(5, "BachCosPA");
    histos.get<TH1>(HIST("Tracked/Rec-Truth/GenFiltersXi"))->GetXaxis()->SetBinLabel(6, "CascCosPA");
    histos.get<TH1>(HIST("Tracked/Rec-Truth/GenFiltersXi"))->GetXaxis()->SetBinLabel(7, "PropDCAxy");
    histos.get<TH1>(HIST("Tracked/Rec-Truth/GenFiltersXi"))->GetXaxis()->SetBinLabel(8, "PropDCAz");
    histos.get<TH1>(HIST("Tracked/Rec-Truth/GenFiltersXi"))->GetXaxis()->SetBinLabel(9, "all");
    histos.get<TH1>(HIST("Tracked/Rec-Truth/FiltersOmega"))->GetXaxis()->SetBinLabel(1, "p_{T}");
    histos.get<TH1>(HIST("Tracked/Rec-Truth/FiltersOmega"))->GetXaxis()->SetBinLabel(2, "gen");
    histos.get<TH1>(HIST("Tracked/Rec-Truth/FiltersOmega"))->GetXaxis()->SetBinLabel(3, "TPC");
    histos.get<TH1>(HIST("Tracked/Rec-Truth/FiltersOmega"))->GetXaxis()->SetBinLabel(4, "TOF");
    histos.get<TH1>(HIST("Tracked/Rec-Truth/FiltersOmega"))->GetXaxis()->SetBinLabel(5, "CMR");
    histos.get<TH1>(HIST("Tracked/Rec-Truth/FiltersOmega"))->GetXaxis()->SetBinLabel(6, "all");
    histos.get<TH1>(HIST("Tracked/Rec-Truth/GenFiltersOmega"))->GetXaxis()->SetBinLabel(1, "casc |y|");
    histos.get<TH1>(HIST("Tracked/Rec-Truth/GenFiltersOmega"))->GetXaxis()->SetBinLabel(2, "dau |#eta|");
    histos.get<TH1>(HIST("Tracked/Rec-Truth/GenFiltersOmega"))->GetXaxis()->SetBinLabel(3, "dau TPC cls");
    histos.get<TH1>(HIST("Tracked/Rec-Truth/GenFiltersOmega"))->GetXaxis()->SetBinLabel(4, "V0CosPA");
    histos.get<TH1>(HIST("Tracked/Rec-Truth/GenFiltersOmega"))->GetXaxis()->SetBinLabel(5, "BachCosPA");
    histos.get<TH1>(HIST("Tracked/Rec-Truth/GenFiltersOmega"))->GetXaxis()->SetBinLabel(6, "CascCosPA");
    histos.get<TH1>(HIST("Tracked/Rec-Truth/GenFiltersOmega"))->GetXaxis()->SetBinLabel(7, "PropDCAxy");
    histos.get<TH1>(HIST("Tracked/Rec-Truth/GenFiltersOmega"))->GetXaxis()->SetBinLabel(8, "PropDCAz");
    histos.get<TH1>(HIST("Tracked/Rec-Truth/GenFiltersOmega"))->GetXaxis()->SetBinLabel(9, "all");
  }

  void processDerivedData(DerCollisionWMults::iterator const& collision, DerCascDatas const& allCascs, DerTraCascDatas const& traCascs, DauTracks const&)
  {
    fillEvents(collision); // save info about all processed events
    //* if (doApplyEfficiency1D || doApplyPurity1D || doApplyEfficiency2D || doApplyPurity2D) {
    //*   initEfficiencyFromCCDB(collision.runNumber(), collision.timestamp());
    //* }
    if (isValidEvent(collision)) {
      histos.fill(HIST("Rec-Events/EvCounter"), 0.5);
      histos.fill(HIST("Rec-Events/PVxy"), collision.posX(), collision.posY());
      histos.fill(HIST("Rec-Events/PVz"), collision.posZ());
      double mult = (doProcessIons) ? collision.centFT0C() : collision.centFT0M();
      histos.fill(HIST("Rec-Events/Mult"), mult);
      analyseCascs(collision, allCascs); // process all cascades
      analyseCascs(collision, traCascs); // process tracked cascades
    }
  }

  void processDerivedMCGen(soa::Join<aod::StraMCCollisions, aod::StraMCCollMults> const& genColls, DerMCGenCascades const& genCascs, soa::Join<aod::StraCollisions, aod::StraCollLabels, aod::StraCents> const& recColls)
  {
    for (auto const& genColl : genColls) {
      histos.fill(HIST("MC/Gen/EvCounter"), 0.5); // generated events statistics
      double genMult = genColl.multMCNParticlesEta05();
      histos.fill(HIST("MC/Gen/Mult"), genMult);
      // loop over generated cascades to obtain pt-genMult distribution
      auto slicedGenCascs = genCascs.sliceBy(cascsPerMcCollision, genColl.globalIndex());
      for (auto const& casc : slicedGenCascs) {
        if (casc.straMCCollisionId() != genColl.index())
             continue; // safety check for correct slicing
        // fill generated Xi and Omega (all) for statistics
        double cascPt = casc.ptMC();
        if (isValidPDG(casc, "Xi"))
            histos.fill(HIST("MC/Gen/Xi"), genMult, cascPt);
        if (isValidPDG(casc, "Omega"))
            histos.fill(HIST("MC/Gen/Omega"), genMult, cascPt);
        if (casc.isPhysicalPrimary()) {
          // fill generated primary Xi and Omega for statistics
            if (isValidPDG(casc, "Xi"))
              histos.fill(HIST("MC/Gen/PrimaryXi"), genMult, cascPt);
            if (isValidPDG(casc, "Omega"))
              histos.fill(HIST("MC/Gen/PrimaryOmega"), genMult, cascPt);
          // fill generated primary Xi within rapidity for corrections
            if (isValidPDG(casc, "Xi") && std::abs(casc.rapidityMC(0)) < selCuts.cutRapidity) {
              histos.fill(HIST("MC/Gen/PrimaryXiRapidity"), genMult, cascPt);
              if (casc.pdgCode() == PDG_t::kXiMinus)
                histos.fill(HIST("MC/Gen/PrimaryXiMinusRapidity"), genMult, cascPt);
              if (casc.pdgCode() == PDG_t::kXiPlusBar)
                histos.fill(HIST("MC/Gen/PrimaryXiPlusRapidity"), genMult, cascPt);
            }
          // fill generated primary Omega within rapidity for corrections
            if (isValidPDG(casc, "Omega") && std::abs(casc.rapidityMC(2)) < selCuts.cutRapidity) {
              histos.fill(HIST("MC/Gen/PrimaryOmegaRapidity"), genMult, cascPt);
              if (casc.pdgCode() == PDG_t::kOmegaMinus)
                histos.fill(HIST("MC/Gen/PrimaryOmegaMinusRapidity"), genMult, cascPt);
              if (casc.pdgCode() == PDG_t::kOmegaPlusBar)
                histos.fill(HIST("MC/Gen/PrimaryOmegaPlusRapidity"), genMult, cascPt);
            }
        }
        // fill gen rapidity for statistics
        if (isValidPDG(casc, "Xi"))
          histos.fill(HIST("MC/Gen/RapidityXi"), casc.rapidityMC(0));
        if (isValidPDG(casc, "Omega"))
          histos.fill(HIST("MC/Gen/RapidityOmega"), casc.rapidityMC(2));
      }
      // look at rec. ev. >= 1 for corrections
      auto slicedRecColls = recColls.sliceBy(perMcCollision, genColl.globalIndex());
      if (slicedRecColls.size() > 0) {
      histos.fill(HIST("MC/EvRec/EvCounter"), 0.5);
      histos.fill(HIST("MC/EvRec/Mult"), genMult);
      }
      for (auto const& recColl : slicedRecColls) {
        int64_t genCollId = recColl.straMCCollisionId();
        if (genColl.index() != genCollId)
            continue; // safety check for correct slicing
        // fill multiplicity-centrality distribution
        double recoMult = (doProcessIons) ? recColl.centFT0C() : recColl.centFT0M();
        histos.fill(HIST("MC/EvRec/MultCent"), recoMult, genMult);
        // look at generated cascades within reconstructed events
        for (auto const& casc : slicedGenCascs) {
          if (casc.straMCCollisionId() != genCollId)
            continue; // check for cascades belonging to the correct reco collision
          double cascPt = casc.ptMC();
          if (isValidPDG(casc, "Xi"))
            histos.fill(HIST("MC/EvRec/Xi"), cascPt, recoMult);
          if (isValidPDG(casc, "Omega"))
            histos.fill(HIST("MC/EvRec/Omega"), cascPt, recoMult);
          if (casc.isPhysicalPrimary()) {
            if (isValidPDG(casc, "Xi")) {
              histos.fill(HIST("MC/EvRec/PrimaryXi"), cascPt, recoMult);
              if (std::abs(casc.rapidityMC(0)) < selCuts.cutRapidity) {
                histos.fill(HIST("MC/EvRec/PrimaryXiRapidity"), cascPt, recoMult);
                if (casc.pdgCode() == PDG_t::kXiMinus)
                  histos.fill(HIST("MC/EvRec/PrimaryXiMinusRapidity"), cascPt, recoMult);
                if (casc.pdgCode() == PDG_t::kXiPlusBar)
                  histos.fill(HIST("MC/EvRec/PrimaryXiPlusRapidity"), cascPt, recoMult);
              }
            }
            if (isValidPDG(casc, "Omega")) {
              histos.fill(HIST("MC/EvRec/PrimaryOmega"), cascPt, recoMult);
              if (std::abs(casc.rapidityMC(2)) < selCuts.cutRapidity) {
                histos.fill(HIST("MC/EvRec/PrimaryOmegaRapidity"), cascPt, recoMult);
                if (casc.pdgCode() == PDG_t::kOmegaMinus)
                  histos.fill(HIST("MC/EvRec/PrimaryOmegaMinusRapidity"), cascPt, recoMult);
                if (casc.pdgCode() == PDG_t::kOmegaPlusBar)
                  histos.fill(HIST("MC/EvRec/PrimaryOmegaPlusRapidity"), cascPt, recoMult);
              }
            }
          }
        }
      }
    }
  }

  void processDerivedMCRec(DerMCRecCollisions::iterator const& collision, DerMCRecCascDatas const& allCascs, DerMCRecTraCascDatas const& traCascs, DauTracks const&, DerMCGenCascades const&)
  {
    fillEvents(collision); // save info about all processed events
    //* if (doApplyEfficiency1D || doApplyPurity1D || doApplyEfficiency2D || doApplyPurity2D) {
    //*   initEfficiencyFromCCDB(collision.runNumber(), collision.timestamp());
    //* }
    if (isValidEvent(collision)) {
      histos.fill(HIST("Rec-Events/EvCounter"), 0.5);
      histos.fill(HIST("Rec-Events/PVxy"), collision.posX(), collision.posY());
      histos.fill(HIST("Rec-Events/PVz"), collision.posZ());
      double mult = (doProcessIons) ? collision.centFT0C() : collision.centFT0M();
      histos.fill(HIST("Rec-Events/Mult"), mult);
      analyseCascs(collision, allCascs); // process all cascades
      analyseCascs(collision, traCascs); // process tracked cascades
    }
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