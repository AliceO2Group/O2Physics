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
using DerCollisionWMult = soa::Join<aod::StraCollisions, aod::StraStamps, aod::StraEvSels, aod::StraCents>::iterator;
using DerCascDatas = soa::Join<aod::CascCores, aod::CascCollRefs, aod::CascExtras, aod::CascBBs>;
using DerTraCascDatas = soa::Join<aod::TraCascCores, aod::TraCascCollRefs, aod::StraTrackExtras, aod::TraToCascRefs>;

// tables for derived MC
using DerMCGenCascades = soa::Join<aod::CascMCCores, aod::CascMCCollRefs>;
using DerMCRecCollision = soa::Join<aod::StraCollisions, aod::StraStamps, aod::StraCollLabels, aod::StraEvSels, aod::StraCents>::iterator;
using DerMCRecCascDatas = soa::Join<aod::CascCollRefs, aod::CascCores, aod::CascCoreMCLabels, aod::CascExtras, aod::CascBBs>;
using DerMCRecTraCascDatas = soa::Join<aod::TraCascCores, aod::TraCascCollRefs, aod::StraTrackExtras, aod::TraToCascRefs>;

// tables for PID selection
using DauTracks = soa::Join<aod::DauTrackExtras, aod::DauTrackTPCPIDs>;

struct StrangeCascTrack {

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Service<o2::framework::O2DatabasePDG> pdgDB;

  PresliceUnsorted<soa::Join<aod::StraCollisions, aod::StraCollLabels, aod::StraCents>> perMcCollision = aod::v0data::straMCCollisionId;

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // subprocess switches:
  Configurable<bool> doProcessPP{"doProcessPP", true, "true for pp, false for PbPb and OO"};
  // selections
  Configurable<bool> doApplyCuts{"doApplyCuts", true, "apply cuts"}; // dca for filtering data primaries
  Configurable<bool> doApplyTPCPID{"doApplyTPCPID", true, "apply tpc pid to dau tracks"};
  Configurable<bool> doApplyTOFPID{"doApplyTOFPID", true, "apply tof pid to dau tracks"};
  Configurable<bool> doCompetingMassRej{"doCompetingMassRej", true, "competing mass rejection for omegas"};
  // corrections
  Configurable<bool> doApplyEfficiency{"doApplyEfficiency", false, "apply efficiency correction"};
  Configurable<bool> doPropagateEfficiency{"doPropagateEfficiency", false, "apply efficiency propagation"};
  Configurable<bool> doApplyPurity{"doApplyPurity", false, "apply purity correction"};
  Configurable<bool> doPropagatePurity{"doPropagatePurity", false, "apply purity propagation"};
  Configurable<std::string> ccdburl{"ccdburl", "http://alice-ccdb.cern.ch", "url of the ccdb repository to use"};
  Configurable<std::string> efficiencyCCDBPath{"efficiencyCCDBPath", "GLO/Config/GeometryAligned", "Path of the efficiency corrections"};
  // mc settings
  Configurable<bool> doFillTruth{"doFillTruth", false, "require MC truth for reco"};

  // event and dau track selection
  struct : ConfigurableGroup {
    Configurable<double> cutDCAtoPVxy{"cutDCAtoPVxy", 0.02f, "max cascade dca to PV in xy"};
    Configurable<double> cutDCAtoPVz{"cutDCAtoPVz", 0.02f, "max cascade dca to PV in z"};
    Configurable<float> compMassRej{"compMassRej", 0.008, "Competing mass rejection"};
    // TPC PID selection
    Configurable<float> NSigmaTPCPion{"NSigmaTPCPion", 4, "NSigmaTPCPion"};
    Configurable<float> NSigmaTPCKaon{"NSigmaTPCKaon", 4, "NSigmaTPCKaon"};
    Configurable<float> NSigmaTPCProton{"NSigmaTPCProton", 4, "NSigmaTPCProton"};
    // TOF PID selection
    Configurable<float> NSigmaTOFPion{"NSigmaTOFPion", 3, "NSigmaTOFPion"};
    Configurable<float> NSigmaTOFKaon{"NSigmaTOFKaon", 3, "NSigmaTOFKaon"};
    Configurable<float> NSigmaTOFProton{"NSigmaTOFProton", 3, "NSigmaTOFProton"};
  } selCuts;

  // axes
  struct : ConfigurableGroup {
    ConfigurableAxis axisPhi{"Phi", {72, 0, TwoPI}, "#phi"};
    ConfigurableAxis axisEta{"Eta", {102, -2.01, 2.01}, "#eta"};
    ConfigurableAxis axisDCAxy{"DCA to xy plane", {500, 0., 0.5}, "cm"};
    ConfigurableAxis axisDCAz{"DCA to z plane", {500, 0., 0.5}, "cm"};
    ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 10.0}, "p_{T} (GeV/c)"};
    ConfigurableAxis axisMult{"axisMult", {VARIABLE_WIDTH, 0.0f, 0.01f, 1.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 70.0f, 100.0f}, "FT0 mult %"};
    ConfigurableAxis axisOmegaMass{"axisOmegaMass", {2000, 1.6, 1.8}, "#Omega M_{inv} (GeV/c^{2})"};
    ConfigurableAxis axisXiMass{"axisXiMass", {2000, 1.2, 1.4}, "#Xi M_{inv} (GeV/c^{2})"};
  } axesConfig;

  // // Filters events
  // if (doFilterEvents) {
  //   Filter eventFilter = (o2::aod::evsel::sel8 == true);
  //   Filter posZFilter = (nabs(o2::aod::collision::posZ) < cutzvertex);
  //   Filter posZFilterMC = (nabs(o2::aod::mccollision::posZ) < cutzvertex);
  // }

  // cascade reconstruction types
  static constexpr std::string_view kTypeNames[] = {"Standard", "Tracked"};

  // for efficiency and purity corrections
  TH2F* hEfficiency;
  TH2F* hEfficiencyUncertainty;
  TH2F* hPurity;
  TH2F* hPurityUncertainty;
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

    TList* listEfficiencies = ccdb->getForTimeStamp<TList>(efficiencyCCDBPath, timeStamp);

    if (!listEfficiencies) {
      LOG(fatal) << "Problem getting TList object with efficiencies and purities!";
    }

    hEfficiency = static_cast<TH2F*>(listEfficiencies->FindObject("hEfficiency"));
    hPurity = static_cast<TH2F*>(listEfficiencies->FindObject("hPurity"));
    hEfficiencyUncertainty = static_cast<TH2F*>(listEfficiencies->FindObject("hEfficiencyUncertainty"));
    hPurityUncertainty = static_cast<TH2F*>(listEfficiencies->FindObject("hPurityUncertainty"));
    if (doPropagateEfficiency && !hEfficiencyUncertainty)
      LOG(fatal) << "Problem getting hEfficiencyUncertainty!";
    if (doPropagatePurity && !hPurityUncertainty)
      LOG(fatal) << "Problem getting hPurityUncertainty!";
    LOG(info) << "Efficiencies and purities now loaded for " << mRunNumber;
  }
  // general info about processed events
  template <typename TEvent>
  void fillEvents(TEvent const& collision)
  {
    histos.fill(HIST("Events/EvCounter"), 0.5);
    double mult = doProcessPP ? collision.centFT0M() : collision.centFT0C();
    histos.fill(HIST("Events/Mult"), mult);
    double pvx = collision.posX();
    double pvy = collision.posY();
    double pvz = collision.posZ();
    histos.fill(HIST("Events/PVx"), pvx);
    histos.fill(HIST("Events/PVy"), pvy);
    histos.fill(HIST("Events/PVz"), pvz);
  }
  // checks general selection criteria
  template <typename TCascade>
  bool isValidCasc(TCascade cascade)
  {
    if (cascade.dcaXYCascToPV() > selCuts.cutDCAtoPVxy)
      return false;
    if (cascade.dcaZCascToPV() > selCuts.cutDCAtoPVz)
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
      if (std::abs(posDaughterTrackCasc.tpcNSigmaPr()) > selCuts.NSigmaTPCProton) {
        return false;
      }
      if (std::abs(negDaughterTrackCasc.tpcNSigmaPi()) > selCuts.NSigmaTPCPion) {
        return false;
      }
    } else {
      if (std::abs(negDaughterTrackCasc.tpcNSigmaPr()) > selCuts.NSigmaTPCProton) {
        return false;
      }
      if (std::abs(posDaughterTrackCasc.tpcNSigmaPi()) > selCuts.NSigmaTPCPion) {
        return false;
      }
    }
    return true;
  }
  // checks TOF PID of dau tracks
  // template <typename TCascade>
  // bool passesTOF(TCascade cascade, TString particle)
  // {
  //   return true;
  //   // const auto& bachDaughterTrackCasc = cascade.bachTrackExtra_as<DauTracks>();
  //   // const auto& posDaughterTrackCasc = cascade.posTrackExtra_as<DauTracks>();
  //   // const auto& negDaughterTrackCasc = cascade.negTrackExtra_as<DauTracks>();
  //   // bool xiPassTOFSelection = true;
  //   // bool omegaPassTOFSelection = true;
  //   // if (cascade.sign() < 0) {
  //   //   if (posDaughterTrackCasc.hasTOF()) {
  //   //     if (std::abs(cascade.tofNSigmaXiLaPr()) > selCuts.NSigmaTOFProton) {
  //   //       xiPassTOFSelection &= false;
  //   //     }
  //   //     if (std::abs(cascade.tofNSigmaOmLaPr()) > selCuts.NSigmaTOFProton) {
  //   //       omegaPassTOFSelection &= false;
  //   //     }
  //   //   }
  //   //   if (negDaughterTrackCasc.hasTOF()) {
  //   //     if (std::abs(cascade.tofNSigmaXiLaPi()) > selCuts.NSigmaTOFPion) {
  //   //       xiPassTOFSelection &= false;
  //   //     }
  //   //     if (std::abs(cascade.tofNSigmaOmLaPi()) > selCuts.NSigmaTOFPion) {
  //   //       omegaPassTOFSelection &= false;
  //   //     }
  //   //   }
  //   // } else {
  //   //   if (posDaughterTrackCasc.hasTOF()) {
  //   //     if (std::abs(cascade.tofNSigmaXiLaPi()) > selCuts.NSigmaTOFPion) {
  //   //       xiPassTOFSelection &= false;
  //   //     }
  //   //     if (std::abs(cascade.tofNSigmaOmLaPi()) > selCuts.NSigmaTOFPion) {
  //   //       omegaPassTOFSelection &= false;
  //   //     }
  //   //   }
  //   //   if (negDaughterTrackCasc.hasTOF()) {
  //   //     if (std::abs(cascade.tofNSigmaXiLaPr()) > selCuts.NSigmaTOFProton) {
  //   //       xiPassTOFSelection &= false;
  //   //     }
  //   //     if (std::abs(cascade.tofNSigmaOmLaPr()) > selCuts.NSigmaTOFProton) {
  //   //       omegaPassTOFSelection &= false;
  //   //     }
  //   //   }
  //   // }
  //
  //   // if (bachDaughterTrackCasc.hasTOF()) {
  //   //   if (std::abs(cascade.tofNSigmaXiPi()) > selCuts.NSigmaTOFPion) {
  //   //     xiPassTOFSelection &= false;
  //   //   }
  //   //   if (std::abs(cascade.tofNSigmaOmKa()) > selCuts.NSigmaTOFKaon) {
  //   //     omegaPassTOFSelection &= false;
  //   //   }
  //   // }
  //
  //   // if (bachDaughterTrackCasc.hasTOF()) {
  //   //   if (std::abs(cascade.tofNSigmaXiPi()) > selCuts.NSigmaTOFPion) {
  //   //     xiPassTOFSelection &= false;
  //   //   }
  //   //   if (std::abs(cascade.tofNSigmaOmKa()) > selCuts.NSigmaTOFKaon) {
  //   //     omegaPassTOFSelection &= false;
  //   //   }
  //   // }
  //
  //   // if (particle == "xi") {return xiPassTOFSelection;} else {return omegaPassTOFSelection;}
  // }
  // checks whether gen cascade corresponds to PDG code
  template <typename TCascade>
  bool isValidPDG(TCascade cascade, TString particle)
  {
    static constexpr int kPdgCodes[] = {3312, 3334}; // "XiMinus", "OmegaMinus"
    if (particle == "xi" && std::abs(cascade.pdgCode()) == kPdgCodes[0])
      return true;
    if (particle == "omega" && std::abs(cascade.pdgCode()) == kPdgCodes[1])
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
        return (pdg == 3312);
      if (particle == "omega")
        return (pdg == 3334);
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

      double mult = doProcessPP ? collision.centFT0M() : collision.centFT0C(); // ion collisions use FT0C for multiplicity, pp uses both

      float efficiency = 1.0f;
      float efficiencyError = 0.0f;
      float purity = 1.0f;
      float purityError = 0.0f;
      if (doApplyEfficiency) {
        efficiency = hEfficiency->Interpolate(cascade.pt(), mult);
        if (doPropagateEfficiency) {
          efficiencyError = hEfficiencyUncertainty->Interpolate(cascade.pt(), mult);
        }
        if (efficiency == 0) { // check for zero efficiency, do not apply if the case
          efficiency = 1.;
          efficiencyError = 0.;
        }
      }
      if (doApplyPurity) {
        purity = hPurity->Interpolate(cascade.pt(), mult);
        if (doPropagatePurity) {
          purityError = hPurityUncertainty->Interpolate(cascade.pt(), mult);
        }
        if (purity == 0) { // check for zero purity, do not apply if the case
          purity = 1.;
          purityError = 0.;
        }
      }

      if (collision.index() != casccollid) {
        histos.fill(HIST(kTypeNames[type]) + HIST("/EvMult"), mult);
        casccollid = collision.index();
      }

      double massXi = cascade.mXi();
      double massOmega = cascade.mOmega();
      double pt = cascade.pt();

      histos.fill(HIST(kTypeNames[type]) + HIST("/DCAxy"), cascade.dcaXYCascToPV());
      histos.fill(HIST(kTypeNames[type]) + HIST("/DCAz"), cascade.dcaZCascToPV());
      histos.fill(HIST(kTypeNames[type]) + HIST("/Phi"), cascade.phi());
      histos.fill(HIST(kTypeNames[type]) + HIST("/Eta"), cascade.eta());
      histos.fill(HIST(kTypeNames[type]) + HIST("/MassXiNoSel"), massXi);
      histos.fill(HIST(kTypeNames[type]) + HIST("/MassOmegaNoSel"), massOmega);

      // start checking selections
      bool passedAllSels = true;
      // apply general selection criteria
      if (doApplyCuts) {
        if (isValidCasc(cascade)) {
          histos.fill(HIST(kTypeNames[type]) + HIST("/MassXiGenSel"), massXi);
          histos.fill(HIST(kTypeNames[type]) + HIST("/MassOmegaGenSel"), massOmega);
        } else {
          passedAllSels = false;
        }
      }
      // apply tpc pid
      if (doApplyTPCPID) {
        if (passesTPC(stdCasc)) {
          histos.fill(HIST(kTypeNames[type]) + HIST("/MassXiTPCPID"), massXi);
          histos.fill(HIST(kTypeNames[type]) + HIST("/MassOmegaTPCPID"), massOmega);
        } else {
          passedAllSels = false;
        }
      }
      // apply tof pid
      bool passedAllSelsXi = passedAllSels;
      bool passedAllSelsOmega = passedAllSels;
      // if (doApplyTOFPID) {
      //   if (passesTOF(cascade, "xi")) {
      //     histos.fill(HIST(kTypeNames[type]) + HIST("/MassXiTOFPID"), massXi);
      //   } else {
      //     passedAllSelsXi = false;
      //   }
      //   if (passesTOF(cascade, "omega")) {
      //     histos.fill(HIST(kTypeNames[type]) + HIST("/MassOmegaTOFPID"), massOmega);
      //   } else {
      //     passedAllSelsOmega = false;
      //   }
      // }
      // apply competing mass rej
      if (doCompetingMassRej) {
        if (std::abs(massXi - pdgDB->Mass(3312)) > selCuts.compMassRej) {
          histos.fill(HIST(kTypeNames[type]) + HIST("/MassOmegaMassSel"), massOmega);
        } else {
          passedAllSelsOmega = false;
        }
      }

      // fill w/ cascs that passed all applied sels
      double binFillXi[3] = {massXi, pt, mult};
      if (passedAllSelsXi) {
        histos.fill(HIST(kTypeNames[type]) + HIST("/MassXi"), massXi);
        fillHist(histos.get<THn>(HIST(kTypeNames[type]) + HIST("/Xi")), binFillXi, efficiency, efficiencyError, purity, purityError);
        if constexpr (requires { collision.straMCCollisionId(); }) {
          if (doFillTruth && isMCTruth(stdCasc, "xi"))
            histos.fill(HIST("MC/RecTruth/") + HIST(kTypeNames[type]) + HIST("/Xi"), massXi, pt, mult);
        }
      }
      double binFillOmega[3] = {massOmega, pt, mult};
      if (passedAllSelsOmega) {
        histos.fill(HIST(kTypeNames[type]) + HIST("/MassOmega"), massOmega);
        fillHist(histos.get<THn>(HIST(kTypeNames[type]) + HIST("/Omega")), binFillOmega, efficiency, efficiencyError, purity, purityError);
        if constexpr (requires { collision.straMCCollisionId(); }) {
          if (doFillTruth && isMCTruth(stdCasc, "omega"))
            histos.fill(HIST("MC/RecTruth/") + HIST(kTypeNames[type]) + HIST("/Omega"), massOmega, pt, mult);
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
      histos.add(Form("%s/Phi", kTypeNames[type].data()), "Phi", kTH1F, {axesConfig.axisPhi});
      histos.add(Form("%s/Eta", kTypeNames[type].data()), "Eta", kTH1F, {axesConfig.axisEta});
      histos.add(Form("%s/DCAxy", kTypeNames[type].data()), "DCA to xy", kTH1F, {axesConfig.axisDCAxy});
      histos.add(Form("%s/DCAz", kTypeNames[type].data()), "DCA to z", kTH1F, {axesConfig.axisDCAz});
      histos.add(Form("%s/EvMult", kTypeNames[type].data()), "Multiplicity of events with >=1 cascade", kTH1F, {axesConfig.axisMult});
      // no selection applied
      histos.add(Form("%s/MassOmegaNoSel", kTypeNames[type].data()), "Invariant mass hypothesis", kTH1F, {axesConfig.axisOmegaMass});
      histos.add(Form("%s/MassXiNoSel", kTypeNames[type].data()), "Invariant mass hypothesis", kTH1F, {axesConfig.axisXiMass});
      // only gen selection applied
      histos.add(Form("%s/MassOmegaGenSel", kTypeNames[type].data()), "Invariant mass hypothesis", kTH1F, {axesConfig.axisOmegaMass});
      histos.add(Form("%s/MassXiGenSel", kTypeNames[type].data()), "Invariant mass hypothesis", kTH1F, {axesConfig.axisXiMass});
      // only tpc pid selection applied
      histos.add(Form("%s/MassOmegaTPCPID", kTypeNames[type].data()), "Invariant mass hypothesis", kTH1F, {axesConfig.axisOmegaMass});
      histos.add(Form("%s/MassXiTPCPID", kTypeNames[type].data()), "Invariant mass hypothesis", kTH1F, {axesConfig.axisXiMass});
      // only tof pid selection applied
      histos.add(Form("%s/MassOmegaTOFPID", kTypeNames[type].data()), "Invariant mass hypothesis", kTH1F, {axesConfig.axisOmegaMass});
      histos.add(Form("%s/MassXiTOFPID", kTypeNames[type].data()), "Invariant mass hypothesis", kTH1F, {axesConfig.axisXiMass});
      // only competing mass rejection selection applied
      histos.add(Form("%s/MassOmegaMassSel", kTypeNames[type].data()), "Invariant mass hypothesis", kTH1F, {axesConfig.axisOmegaMass});
      // passed all applied sels
      histos.add(Form("%s/MassOmega", kTypeNames[type].data()), "Invariant mass hypothesis", kTH1F, {axesConfig.axisOmegaMass});
      histos.add(Form("%s/MassXi", kTypeNames[type].data()), "Invariant mass hypothesis", kTH1F, {axesConfig.axisXiMass});
      histos.add(Form("%s/Omega", kTypeNames[type].data()), "", kTHnD, {axesConfig.axisOmegaMass, axesConfig.axisPt, axesConfig.axisMult});
      histos.add(Form("%s/Xi", kTypeNames[type].data()), "", kTHnD, {axesConfig.axisXiMass, axesConfig.axisPt, axesConfig.axisMult});
    });
    // for MC-specific processing
    histos.add("MC/Gen/EvCounter", "Event Counter", kTH1F, {{1, 0, 1}});
    histos.add("MC/Gen/Xi", "Xi", kTH2F, {axesConfig.axisPt, axesConfig.axisMult});       // generated primary Xis
    histos.add("MC/Gen/Omega", "Omega", kTH2F, {axesConfig.axisPt, axesConfig.axisMult}); // generated primary Omegas
    histos.add("MC/Rec/EvCounter", "Event Counter", kTH1F, {{1, 0, 1}});
    histos.add("MC/Rec/EvMult", "Multiplicity", kTH1F, {axesConfig.axisMult});
    histos.add("MC/RecTruth/Standard/Omega", "", kTHnD, {axesConfig.axisOmegaMass, axesConfig.axisPt, axesConfig.axisMult});
    histos.add("MC/RecTruth/Standard/Xi", "", kTHnD, {axesConfig.axisXiMass, axesConfig.axisPt, axesConfig.axisMult});
    histos.add("MC/RecTruth/Tracked/Omega", "", kTHnD, {axesConfig.axisOmegaMass, axesConfig.axisPt, axesConfig.axisMult});
    histos.add("MC/RecTruth/Tracked/Xi", "", kTHnD, {axesConfig.axisXiMass, axesConfig.axisPt, axesConfig.axisMult});
  }

  void processDerivedData(DerCollisionWMult const& collision, DerCascDatas const& allCascs, DerTraCascDatas const& traCascs, DauTracks const&)
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
        double casc_mult = doProcessPP ? recColl.centFT0M() : recColl.centFT0C();
        histos.fill(HIST("MC/Rec/EvMult"), casc_mult);
        int64_t genCollId = recColl.straMCCollisionId();
        for (auto const& casc : genCascs) {
          if (casc.straMCCollisionId() != genCollId)
            continue; // safety check
          if (!casc.isPhysicalPrimary())
            continue; // skip non-primary particles
          double casc_pt = std::sqrt(std::pow(casc.pxMC(), 2) + std::pow(casc.pyMC(), 2));
          if (isValidPDG(casc, "xi"))
            histos.fill(HIST("MC/Gen/Xi"), casc_pt, casc_mult);
          if (isValidPDG(casc, "omega"))
            histos.fill(HIST("MC/Gen/Omega"), casc_pt, casc_mult);
        }
      }
    }
  }

  void processDerivedMCRec(DerMCRecCollision const& collision, DerMCRecCascDatas const& allCascs, DerMCRecTraCascDatas const& traCascs, DauTracks const&, DerMCGenCascades const&)
  {
    fillEvents(collision); // save info about all processed events
    if (doApplyEfficiency) {
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
