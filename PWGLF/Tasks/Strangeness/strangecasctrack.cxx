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

struct strangecasctrack {

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Service<o2::framework::O2DatabasePDG> pdgDB;

  PresliceUnsorted<soa::Join<aod::StraCollisions, aod::StraCollLabels, aod::StraCents>> perMcCollision = aod::v0data::straMCCollisionId;

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // subprocess switches:
  Configurable<bool> doProcesspp{"doProcesspp", true, "true for pp"};
  Configurable<bool> doProcessPbPb{"doProcessPbPb", false, "true for PbPb"};
  Configurable<bool> doProcessOO{"doProcessOO", false, "true for OO"};
  Configurable<bool> doProcesspO{"doProcesspO", false, "true for pO"};
  // selections
  Configurable<bool> doApplyCuts{"doApplyCuts", true, "apply cuts"}; // general cascade cuts - dca, cosPA etc.
  Configurable<bool> doApplyTPCPID{"doApplyTPCPID", true, "apply tpc pid to dau tracks"};
  Configurable<bool> doApplyTOFPID{"doApplyTOFPID", true, "apply tof pid to dau tracks"};
  Configurable<bool> doCompetingMassRej{"doCompetingMassRej", true, "competing mass rejection for omegas"};
  // corrections
  Configurable<bool> doApplyEfficiency{"doApplyEfficiency", false, "apply efficiency correction"};
  Configurable<bool> doPropagateEfficiency{"doPropagateEfficiency", false, "apply efficiency propagation"};
  Configurable<bool> doApplyPurity{"doApplyPurity", false, "apply purity correction"};
  Configurable<bool> doPropagatePurity{"doPropagatePurity", false, "apply purity propagation"};
  Configurable<std::string> ccdburl{"ccdburl", "http://alice-ccdb.cern.ch", "url of the ccdb repository to use"};
  Configurable<std::string> efficiencyCCDBPath_pp{"efficiencyCCDBPath_pp", "Users/y/yparovia/LHC24f4d/Efficiency", "Path of the efficiency corrections"};
  Configurable<std::string> efficiencyCCDBPath_PbPb{"efficiencyCCDBPath_PbPb", "Users/y/yparovia/LHC25f3/Efficiency", "Path of the efficiency corrections"};
  Configurable<std::string> efficiencyCCDBPath_OO{"efficiencyCCDBPath_OO", "Users/y/yparovia/LHC25h3/Efficiency", "Path of the efficiency corrections"};
  Configurable<std::string> efficiencyCCDBPath_pO{"efficiencyCCDBPath_pO", "Users/y/yparovia/LHC25h2/Efficiency", "Path of the efficiency corrections"};

  // event and dau track selection
  struct : ConfigurableGroup {
    Configurable<double> cutZVertex{"cutZVertex", 10.0f, "max Z-vertex position"};
    Configurable<double> cutDCAtoPVxy{"cutDCAtoPVxy", 0.02f, "max cascade dca to PV in xy"};
    Configurable<double> cutDCAtoPVz{"cutDCAtoPVz", 0.02f, "max cascade dca to PV in z"};
    Configurable<double> compMassRej{"compMassRej", 0.008f, "Competing mass rejection"};
    Configurable<double> cutV0CosPA{"cutV0CosPA", 0.97f, "max V0 cosPA"};
    Configurable<double> cutBachCosPA{"cutBachCosPA", 0.97f, "max Bachelor cosPA"};
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

  // Filters events
  // Filter eventFilter = (o2::aod::evsel::sel8 == true);
  // Filter posZFilter = (nabs(o2::aod::collision::posZ) < selCuts.cutZVertex);
  // Filter posZFilterMC = (nabs(o2::aod::mccollision::posZ) < selCuts.cutZVertex);

  // cascade reconstruction types
  static constexpr std::string_view kTypeNames[] = {"Standard", "Tracked"};

  // for efficiency and purity corrections
  TH2F* hEfficiencyOmegaStd;
  TH2F* hEfficiencyOmegaTra;
  TH2F* hEfficiencyXiStd;
  TH2F* hEfficiencyXiTra;

  TH2F* hEfficiencyErrOmegaStd;
  TH2F* hEfficiencyErrOmegaTra;
  TH2F* hEfficiencyErrXiStd;
  TH2F* hEfficiencyErrXiTra;

  TH2F* hPurityOmegaStd;
  TH2F* hPurityOmegaTra;
  TH2F* hPurityXiStd;
  TH2F* hPurityXiTra;

  TH2F* hPurityErrOmegaStd;
  TH2F* hPurityErrOmegaTra;
  TH2F* hPurityErrXiStd;
  TH2F* hPurityErrXiTra;

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
        return efficiencyCCDBPath_pp;
      } else if (doProcesspO) {
        return efficiencyCCDBPath_pO;
      } else if (doProcessPbPb) {
        return efficiencyCCDBPath_PbPb;
      }
      return efficiencyCCDBPath_OO;
    }();

    TList* listEfficiencies = ccdb->getForTimeStamp<TList>(efficiencyCCDBPath, timeStamp);

    if (!listEfficiencies) {
      LOG(fatal) << "Problem getting TList object with efficiencies and purities!";
    }

    hEfficiencyOmegaStd = static_cast<TH2F*>(listEfficiencies->FindObject("Eff_Omega_Standard"));
    hEfficiencyOmegaTra = static_cast<TH2F*>(listEfficiencies->FindObject("Eff_Omega_Tracked"));
    hEfficiencyXiStd = static_cast<TH2F*>(listEfficiencies->FindObject("Eff_Xi_Standard"));
    hEfficiencyXiTra = static_cast<TH2F*>(listEfficiencies->FindObject("Eff_Xi_Tracked"));

    hEfficiencyErrOmegaStd = static_cast<TH2F*>(listEfficiencies->FindObject("EffErr_Omega_Standard"));
    hEfficiencyErrOmegaTra = static_cast<TH2F*>(listEfficiencies->FindObject("EffErr_Omega_Tracked"));
    hEfficiencyErrXiStd = static_cast<TH2F*>(listEfficiencies->FindObject("EffErr_Xi_Standard"));
    hEfficiencyErrXiTra = static_cast<TH2F*>(listEfficiencies->FindObject("EffErr_Xi_Tracked"));

    hPurityOmegaStd = static_cast<TH2F*>(listEfficiencies->FindObject("Pur_Omega_Standard"));
    hPurityOmegaTra = static_cast<TH2F*>(listEfficiencies->FindObject("Pur_Omega_Tracked"));
    hPurityXiStd = static_cast<TH2F*>(listEfficiencies->FindObject("Pur_Xi_Standard"));
    hPurityXiTra = static_cast<TH2F*>(listEfficiencies->FindObject("Pur_Xi_Tracked"));

    hPurityErrOmegaStd = static_cast<TH2F*>(listEfficiencies->FindObject("PurErr_Omega_Standard"));
    hPurityErrOmegaTra = static_cast<TH2F*>(listEfficiencies->FindObject("PurErr_Omega_Tracked"));
    hPurityErrXiStd = static_cast<TH2F*>(listEfficiencies->FindObject("PurErr_Xi_Standard"));
    hPurityErrXiTra = static_cast<TH2F*>(listEfficiencies->FindObject("PurErr_Xi_Tracked"));

    if (doPropagateEfficiency && (!hEfficiencyErrOmegaStd || !hEfficiencyErrOmegaTra || !hEfficiencyErrXiStd || !hEfficiencyErrXiTra))
      LOG(fatal) << "Problem getting hEfficiencyUncertainty!";
    if (doPropagatePurity && (!hPurityErrOmegaStd || !hPurityErrOmegaTra || !hPurityErrXiStd || !hPurityErrXiTra))
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
  // checks general selection criteria
  template <typename TEvent, typename TCascade>
  bool isValidCasc(TEvent collision, TCascade cascade)
  {
    if (std::abs(collision.posZ()) > selCuts.cutZVertex)
      return false;
    if (cascade.dcaXYCascToPV() > selCuts.cutDCAtoPVxy)
      return false;
    if (cascade.dcaZCascToPV() > selCuts.cutDCAtoPVz)
      return false;
    if (cascade.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) < selCuts.cutV0CosPA)
      return false;
    if (cascade.bachBaryonCosPA() < selCuts.cutBachCosPA)
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
  // applies selections and fills histograms
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

      if (doApplyEfficiency) {
        if constexpr (requires { cascade.topologyChi2(); }) {
          efficiencyOmega = hEfficiencyOmegaTra->Interpolate(cascade.pt(), mult);
          efficiencyXi = hEfficiencyXiTra->Interpolate(cascade.pt(), mult);
          if (doPropagateEfficiency) {
            efficiencyOmegaErr = hEfficiencyErrOmegaTra->Interpolate(cascade.pt(), mult);
            efficiencyXiErr = hEfficiencyErrXiTra->Interpolate(cascade.pt(), mult);
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
          efficiencyOmega = hEfficiencyOmegaStd->Interpolate(cascade.pt(), mult);
          efficiencyXi = hEfficiencyXiStd->Interpolate(cascade.pt(), mult);
          if (doPropagateEfficiency) {
            efficiencyOmegaErr = hEfficiencyErrOmegaStd->Interpolate(cascade.pt(), mult);
            efficiencyXiErr = hEfficiencyErrXiStd->Interpolate(cascade.pt(), mult);
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

      if (doApplyPurity) {
        if constexpr (requires { cascade.topologyChi2(); }) {
          purityOmega = hPurityOmegaTra->Interpolate(cascade.pt(), mult);
          purityXi = hPurityXiTra->Interpolate(cascade.pt(), mult);
          if (doPropagatePurity) {
            purityOmegaErr = hPurityErrOmegaTra->Interpolate(cascade.pt(), mult);
            purityXiErr = hPurityErrXiTra->Interpolate(cascade.pt(), mult);
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
          purityOmega = hPurityOmegaStd->Interpolate(cascade.pt(), mult);
          purityXi = hPurityXiStd->Interpolate(cascade.pt(), mult);
          if (doPropagatePurity) {
            purityOmegaErr = hPurityErrOmegaStd->Interpolate(cascade.pt(), mult);
            purityXiErr = hPurityErrXiStd->Interpolate(cascade.pt(), mult);
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
      bool passedAllSels = true;
      // apply general selection criteria
      if (doApplyCuts) {
        if (!isValidCasc(collision, stdCasc))
          passedAllSels = false;
      }
      // apply tpc pid
      if (doApplyTPCPID) {
        if (!passesTPC(stdCasc))
          passedAllSels = false;
      }
      // apply tof pid
      bool passedAllSelsXi = passedAllSels;
      bool passedAllSelsOmega = passedAllSels;
      if (doApplyTOFPID) {
        if (!passesTOF(stdCasc, "xi"))
          passedAllSelsXi = false;
        if (!passesTOF(stdCasc, "omega"))
          passedAllSelsOmega = false;
      }
      // apply competing mass rej
      if (doCompetingMassRej) {
        if (!(std::abs(massXi - o2::constants::physics::MassXiMinus) > selCuts.compMassRej))
          passedAllSelsOmega = false;
      }

      // fill truth w/ cascs that passed all applied sels
      double binFillXi[3] = {massXi, pt, mult};

      if constexpr (requires { collision.straMCCollisionId(); }) {
        if (passedAllSels && (passedAllSelsXi || passedAllSelsOmega)) { // fill once for every desired cascade
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
    if (doApplyEfficiency) {
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
    if (doApplyEfficiency) {
      initEfficiencyFromCCDB(collision.runNumber(), collision.timestamp());
    }
    analyseCascs(collision, allCascs); // process all cascades
    analyseCascs(collision, traCascs); // process tracked cascades
  }

  PROCESS_SWITCH(strangecasctrack, processDerivedData, "process derived data", true);
  PROCESS_SWITCH(strangecasctrack, processDerivedMCGen, "process derived generated mc data", false);
  PROCESS_SWITCH(strangecasctrack, processDerivedMCRec, "process derived reconstructed mc data", false); // mc and data are mutually exclusive!
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<strangecasctrack>(cfgc),
  };
}
