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
//
/// \file lmeeHFCocktail.cxx
/// \analysis task for lmee heavy flavour cocktail
/// \author Daniel Samitz, <daniel.samitz@cern.ch>, SMI Vienna
///         Elisa Meninno, <elisa.meninno@cern.ch>, SMI Vienna

#include <vector>

#include "Math/Vector4D.h"
#include "MathUtils/Utils.h"
#include "Framework/Task.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "PWGEM/Dilepton/Utils/MCUtilities.h"
#include "PWGEM/Dilepton/DataModel/dileptonTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::pwgem::dilepton::utils::mcutil;

using McParticlesSmeared = soa::Join<aod::McParticles, aod::SmearedElectrons>;

enum EFromHFType {
  kNoE = -1,
  kNoHFE = 0,
  kCE = 1,
  kBE = 2,
  kBCE = 3
};

namespace o2::aod
{
namespace hftable
{
DECLARE_SOA_COLUMN(IsHF, isHF, int);
DECLARE_SOA_COLUMN(BHadronId, bHadronId, int);
DECLARE_SOA_COLUMN(CHadronId, cHadronId, int);
DECLARE_SOA_COLUMN(BQuarkId, bQuarkId, int);
DECLARE_SOA_COLUMN(CQuarkId, cQuarkId, int);
DECLARE_SOA_COLUMN(BQuarkOriginId, bQuarkOriginId, int);
DECLARE_SOA_COLUMN(CQuarkOriginId, cQuarkOriginId, int);
DECLARE_SOA_COLUMN(BQuarkRap, bQuarkRap, float);
DECLARE_SOA_COLUMN(CQuarkRap, cQuarkRap, float);
} // namespace hftable
DECLARE_SOA_TABLE(HfTable, "AOD", "HFTABLE",
                  hftable::IsHF,
                  hftable::BHadronId,
                  hftable::CHadronId,
                  hftable::BQuarkId,
                  hftable::CQuarkId,
                  hftable::BQuarkOriginId,
                  hftable::CQuarkOriginId,
                  hftable::BQuarkRap,
                  hftable::CQuarkRap);
} // namespace o2::aod

const int Nstages = 3;
const char* stageNames[Nstages] = {"gen", "meas", "meas_and_acc"};
// gen: normal weights applied. By default set to 1.
// meas: normal and efficiency weights applied. By default set to 1.
// meas_acc: normal, efficiency weights and acceptance cuts applied

template <typename T>
void doQuark(T& p, std::vector<std::shared_ptr<TH1>> hRapQuark, float ptMin, float etaMax, int pdg)
{
  float pt[Nstages] = {p.pt(), p.ptSmeared(), p.ptSmeared()};
  float eta[Nstages] = {p.eta(), p.etaSmeared(), p.etaSmeared()};
  float cut_pt[Nstages] = {0., 0., ptMin};
  float cut_eta[Nstages] = {9999., 99999., etaMax};
  for (int i = 0; i < Nstages; i++) {
    if (pt[i] > cut_pt[i] && fabs(eta[i]) < cut_eta[i]) {
      if (pdg == 4)
        hRapQuark[i]->Fill(p.cQuarkRap());
      else if (pdg == 5)
        hRapQuark[i]->Fill(p.bQuarkRap());
      else
        hRapQuark[i]->Fill(999.);
    }
  }
}

template <typename T>
void doSingle(T& p, std::vector<std::shared_ptr<TH1>> hEta, std::vector<std::shared_ptr<TH1>> hPt, std::vector<std::shared_ptr<TH2>> hPtEta, float ptMin, float etaMax)
{
  float weight[Nstages] = {p.weight(), p.efficiency() * p.weight(), p.efficiency() * p.weight()};
  float pt[Nstages] = {p.pt(), p.ptSmeared(), p.ptSmeared()};
  float eta[Nstages] = {p.eta(), p.etaSmeared(), p.etaSmeared()};
  float cut_pt[Nstages] = {0., 0., ptMin};
  float cut_eta[Nstages] = {9999., 99999., etaMax};
  for (int i = 0; i < Nstages; i++) {
    if (pt[i] > cut_pt[i] && fabs(eta[i]) < cut_eta[i]) {
      hEta[i]->Fill(eta[i], weight[i]);
      hPt[i]->Fill(pt[i], weight[i]);
      hPtEta[i]->Fill(pt[i], eta[i], weight[i]);
    }
  }
}

template <typename T>
void doPair(T& p1, T& p2, std::vector<std::shared_ptr<TH1>> hMee, std::vector<std::shared_ptr<TH2>> hMeePtee, float ptMin, float etaMax, bool apply_detadphi, float min_deta, float min_dphi)
{

  ROOT::Math::PtEtaPhiMVector v1(p1.ptSmeared(), p1.etaSmeared(), p1.phiSmeared(), o2::constants::physics::MassElectron);
  ROOT::Math::PtEtaPhiMVector v2(p2.ptSmeared(), p2.etaSmeared(), p2.phiSmeared(), o2::constants::physics::MassElectron);
  ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
  ROOT::Math::PtEtaPhiMVector v1_gen(p1.pt(), p1.eta(), p1.phi(), o2::constants::physics::MassElectron);
  ROOT::Math::PtEtaPhiMVector v2_gen(p2.pt(), p2.eta(), p2.phi(), o2::constants::physics::MassElectron);
  ROOT::Math::PtEtaPhiMVector v12_gen = v1_gen + v2_gen;

  double mass[Nstages] = {v12_gen.M(), v12.M(), v12.M()};
  double pt[Nstages] = {v12_gen.Pt(), v12.Pt(), v12.Pt()};
  float pt1[Nstages] = {p1.pt(), p1.ptSmeared(), p1.ptSmeared()};
  float pt2[Nstages] = {p2.pt(), p2.ptSmeared(), p2.ptSmeared()};
  float eta1[Nstages] = {p1.eta(), p1.etaSmeared(), p1.etaSmeared()};
  float eta2[Nstages] = {p2.eta(), p2.etaSmeared(), p2.etaSmeared()};
  float weight[Nstages] = {p1.weight() * p2.weight(), p1.efficiency() * p2.efficiency() * p1.weight() * p2.weight(), p1.efficiency() * p2.efficiency() * p1.weight() * p2.weight()};
  float cut_pt[Nstages] = {0., 0., ptMin};
  float cut_eta[Nstages] = {9999., 99999., etaMax};

  float deta = v1.Eta() - v2.Eta();
  float dphi = v1.Phi() - v2.Phi();
  o2::math_utils::bringToPMPi(dphi);
  if (apply_detadphi && std::pow(deta / min_deta, 2) + std::pow(dphi / min_dphi, 2) < 1.f) {
    return;
  }

  for (int i = 0; i < Nstages; i++) {
    if (pt1[i] > cut_pt[i] && pt2[i] > cut_pt[i] && fabs(eta1[i]) < cut_eta[i] && fabs(eta2[i]) < cut_eta[i]) {
      hMee[i]->Fill(mass[i], weight[i]);
      hMeePtee[i]->Fill(mass[i], pt[i], weight[i]);
    }
  }
}

struct MyConfigs : ConfigurableGroup {
  Configurable<float> fConfigPtMin{"cfgPtMin", 0.2, "min. pT of single electrons"};
  Configurable<float> fConfigEtaMax{"cfgEtaMax", 0.8, "max. |eta| of single electrons"};
  ConfigurableAxis fConfigPtBins{"cfgPtBins", {200, 0.f, 10.f}, "single electron pT binning"};
  ConfigurableAxis fConfigEtaBins{"cfgEtaBins", {200, -10.f, 10.f}, "single electron eta binning"};
  ConfigurableAxis fConfigRapBins{"cfgRapBins", {200, -10.f, 10.f}, "Quark rapidity binning"};
  ConfigurableAxis fConfigMeeBins{"cfgMeeBins", {800, 0.f, 8.f}, "Mee binning"};
  ConfigurableAxis fConfigPteeBins{"cfgPteeBins", {400, 0.f, 10.f}, "pTee binning"};
  Configurable<bool> fConfigApplyDEtaDPhi{"cfgApplyDEtaDPhi", false, "flag to apply deta-phi cut"};
  Configurable<float> fConfigMinDEta{"cfgMinDEta", 0.02, "minimum deta"};
  Configurable<float> fConfigMinDPhi{"cfgMinDPhi", 0.2, "minimum dphi"};
};

struct lmeehfcocktailprefilter {

  Produces<o2::aod::HfTable> hfTable;
  void process(aod::McParticles const& mcParticles)
  {
    for (auto const& p : mcParticles) {

      if (abs(p.pdgCode()) != 11 || o2::mcgenstatus::getHepMCStatusCode(p.statusCode()) != 1 || !p.has_mothers()) {
        hfTable(EFromHFType::kNoE, -1, -1, -1, -1, -1, -1, -999., -999.);
        continue;
      }

      int mother_pdg = mcParticles.iteratorAt(p.mothersIds()[0]).pdgCode();
      bool direct_charm_mother = abs(mother_pdg) < 1e+9 && (std::to_string(mother_pdg)[std::to_string(mother_pdg).length() - 3] == '4' || std::to_string(mother_pdg)[std::to_string(mother_pdg).length() - 4] == '4');
      if (abs(mother_pdg) == 443) {
        direct_charm_mother = false; // we don't want JPsi here
      }
      int cHadronId = -1;
      if (direct_charm_mother) {
        cHadronId = p.mothersIds()[0];
      }
      bool direct_beauty_mother = abs(mother_pdg) < 1e+9 && (std::to_string(mother_pdg)[std::to_string(mother_pdg).length() - 3] == '5' || std::to_string(mother_pdg)[std::to_string(mother_pdg).length() - 4] == '5');
      int bHadronId = IsFromBeauty(p, mcParticles);

      int bQuarkId = -1;
      int cQuarkId = -1;
      int cQuarkOriginId = -1;
      int bQuarkOriginId = -1;

      float bQuarkRap = -999.;
      float cQuarkRap = -999.;

      int isHF = EFromHFType::kNoHFE;

      if (direct_beauty_mother) { // b->e
        isHF = EFromHFType::kBE;
        auto bHadron = mcParticles.iteratorAt(bHadronId);
        bQuarkId = searchMothers(bHadron, mcParticles, 5, true);
        if (bQuarkId > -1) {
          auto bQuark = mcParticles.iteratorAt(bQuarkId);
          bQuarkOriginId = searchMothers(bQuark, mcParticles, 5, false);
          bQuarkRap = bQuark.y();
        }
      } else if (bHadronId > -1 && direct_charm_mother) { // b->c->e
        isHF = EFromHFType::kBCE;
        auto bHadron = mcParticles.iteratorAt(bHadronId);
        bQuarkId = searchMothers(bHadron, mcParticles, 5, true);
        if (bQuarkId > -1) {
          auto bQuark = mcParticles.iteratorAt(bQuarkId);
          bQuarkOriginId = searchMothers(bQuark, mcParticles, 5, false);
          bQuarkRap = bQuark.y();
        }
        auto cHadron = mcParticles.iteratorAt(cHadronId);
        cQuarkId = searchMothers(cHadron, mcParticles, 4, true);
        if (cQuarkId > -1) {
          auto cQuark = mcParticles.iteratorAt(cQuarkId);
          cQuarkOriginId = searchMothers(cQuark, mcParticles, 4, false);
          cQuarkRap = cQuark.y();
        }
      } else if (bHadronId < 0 && direct_charm_mother) { // c->e
        isHF = EFromHFType::kCE;
        auto cHadron = mcParticles.iteratorAt(cHadronId);
        cQuarkId = searchMothers(cHadron, mcParticles, 4, true);
        if (cQuarkId > -1) {
          auto cQuark = mcParticles.iteratorAt(cQuarkId);
          cQuarkOriginId = searchMothers(cQuark, mcParticles, 4, false);
          cQuarkRap = cQuark.y();
        }
      }

      hfTable(isHF, bHadronId, cHadronId, bQuarkId, cQuarkId, bQuarkOriginId, cQuarkOriginId, bQuarkRap, cQuarkRap);
    }
  }
};

struct lmeehfcocktailbeauty {

  HistogramRegistry registry{"registry", {}};

  std::vector<std::vector<std::shared_ptr<TH1>>> hRapQuark;
  std::vector<std::vector<std::shared_ptr<TH1>>> hEta, hPt;
  std::vector<std::vector<std::shared_ptr<TH2>>> hPtEta;
  std::vector<std::shared_ptr<TH1>> hLS_Mee, hULS_Mee;
  std::vector<std::shared_ptr<TH2>> hLS_MeePtee, hULS_MeePtee;

  std::vector<std::shared_ptr<TH1>> hLS_Mee_wPartonicCheck, hULS_Mee_wPartonicCheck;
  std::vector<std::shared_ptr<TH2>> hLS_MeePtee_wPartonicCheck, hULS_MeePtee_wPartonicCheck;

  MyConfigs myConfigs;

  Filter hfFilter = o2::aod::hftable::isHF == static_cast<int>(EFromHFType::kBE) || o2::aod::hftable::isHF == static_cast<int>(EFromHFType::kBCE);
  using MyFilteredMcParticlesSmeared = soa::Filtered<soa::Join<aod::McParticles, aod::SmearedElectrons, aod::HfTable>>;

  Preslice<MyFilteredMcParticlesSmeared> perCollision = aod::mcparticle::mcCollisionId;

  Partition<MyFilteredMcParticlesSmeared> Electrons = (aod::mcparticle::pdgCode == 11);
  Partition<MyFilteredMcParticlesSmeared> Positrons = (aod::mcparticle::pdgCode == -11);

  void init(o2::framework::InitContext&)
  {
    registry.add<TH1>("NEvents", "NEvents", HistType::kTH1F, {{1, 0, 1}}, false);

    const int Nchannels = 2;
    const char* typeNamesSingle[Nchannels] = {"be", "bce"};
    const char* typeTitlesSingle[Nchannels] = {"b->e", "b->c->e"};

    AxisSpec rap_axis = {myConfigs.fConfigRapBins, "y_{b}"};
    AxisSpec eta_axis = {myConfigs.fConfigEtaBins, "#eta_{e}"};
    AxisSpec pt_axis = {myConfigs.fConfigPtBins, "#it{p}_{T,e} (GeV/c)"};
    AxisSpec mass_axis = {myConfigs.fConfigMeeBins, "m_{ee} (GeV/c^{2})"};
    AxisSpec ptee_axis = {myConfigs.fConfigPteeBins, "#it{p}_{T,ee} (GeV/c)"};

    // quark histograms
    for (int i = 0; i < Nchannels; i++) {
      std::vector<std::shared_ptr<TH1>> hRap_temp;
      for (int j = 0; j < Nstages; j++) {
        hRap_temp.push_back(registry.add<TH1>(Form("BeautyQuark_Rap_%s_%s", typeNamesSingle[i], stageNames[j]), Form("Rap Beauty Quark %s %s", typeTitlesSingle[i], stageNames[j]), HistType::kTH1F, {rap_axis}, true));
      }
      hRapQuark.push_back(hRap_temp);
    }

    // single histograms
    for (int i = 0; i < Nchannels; i++) {
      std::vector<std::shared_ptr<TH1>> hEta_temp, hPt_temp;
      std::vector<std::shared_ptr<TH2>> hPtEta_temp;
      for (int j = 0; j < Nstages; j++) {
        hEta_temp.push_back(registry.add<TH1>(Form("Electron_Eta_%s_%s", typeNamesSingle[i], stageNames[j]), Form("Single Electron Eta %s %s", typeTitlesSingle[i], stageNames[j]), HistType::kTH1F, {eta_axis}, true));
        hPt_temp.push_back(registry.add<TH1>(Form("Electron_Pt_%s_%s", typeNamesSingle[i], stageNames[j]), Form("Single Electron Pt %s %s", typeTitlesSingle[i], stageNames[j]), HistType::kTH1F, {pt_axis}, true));
        hPtEta_temp.push_back(registry.add<TH2>(Form("Electron_PtEta_%s_%s", typeNamesSingle[i], stageNames[j]), Form("Single Electron Pt vs. Eta %s %s", typeTitlesSingle[i], stageNames[j]), HistType::kTH2F, {pt_axis, eta_axis}, true));
      }
      hEta.push_back(hEta_temp);
      hPt.push_back(hPt_temp);
      hPtEta.push_back(hPtEta_temp);
    }

    // pair histograms
    // ULS
    for (int j = 0; j < Nstages; j++) {
      hULS_Mee.push_back(registry.add<TH1>(Form("ULS_Mee_%s", stageNames[j]), Form("ULS Mee %s", stageNames[j]), HistType::kTH1F, {mass_axis}, true));
      hULS_MeePtee.push_back(registry.add<TH2>(Form("ULS_MeePtee_%s", stageNames[j]), Form("ULS Mee vs. Ptee %s", stageNames[j]), HistType::kTH2F, {mass_axis, ptee_axis}, true));

      hULS_Mee_wPartonicCheck.push_back(registry.add<TH1>(Form("ULS_Mee_wPartonicCheck_%s", stageNames[j]), Form("ULS Mee wPartonicCheck %s", stageNames[j]), HistType::kTH1F, {mass_axis}, true));
      hULS_MeePtee_wPartonicCheck.push_back(registry.add<TH2>(Form("ULS_MeePtee_wPartonicCheck_%s", stageNames[j]), Form("ULS Mee vs. Ptee wPartonicCheck %s", stageNames[j]), HistType::kTH2F, {mass_axis, ptee_axis}, true));
    }
    // LS
    for (int j = 0; j < Nstages; j++) {
      hLS_Mee.push_back(registry.add<TH1>(Form("LS_Mee_%s", stageNames[j]), Form("LS Mee %s", stageNames[j]), HistType::kTH1F, {mass_axis}, true));
      hLS_MeePtee.push_back(registry.add<TH2>(Form("LS_MeePtee_%s", stageNames[j]), Form("LS Mee vs. Ptee %s", stageNames[j]), HistType::kTH2F, {mass_axis, ptee_axis}, true));

      hLS_Mee_wPartonicCheck.push_back(registry.add<TH1>(Form("LS_Mee_wPartonicCheck_%s", stageNames[j]), Form("LS Mee wPartonicCheck %s", stageNames[j]), HistType::kTH1F, {mass_axis}, true));
      hLS_MeePtee_wPartonicCheck.push_back(registry.add<TH2>(Form("LS_MeePtee_wPartonicCheck_%s", stageNames[j]), Form("LS Mee vs. Ptee wPartonicCheck %s", stageNames[j]), HistType::kTH2F, {mass_axis, ptee_axis}, true));
    }
  }

  void processBeauty(aod::McCollisions const& collisions, MyFilteredMcParticlesSmeared const& mcParticles)
  {
    for (auto const& p : mcParticles) {
      int from_quark = p.isHF() - 2;
      doSingle(p, hEta[from_quark], hPt[from_quark], hPtEta[from_quark], myConfigs.fConfigPtMin, myConfigs.fConfigEtaMax);
      doQuark(p, hRapQuark[from_quark], myConfigs.fConfigPtMin, myConfigs.fConfigEtaMax, 5);
    }

    for (auto const& collision : collisions) {

      registry.fill(HIST("NEvents"), 0.5);

      auto const electronsGrouped = Electrons->sliceBy(perCollision, collision.globalIndex());
      auto const positronsGrouped = Positrons->sliceBy(perCollision, collision.globalIndex());
      // ULS spectrum
      for (auto const& [particle1, particle2] : combinations(o2::soa::CombinationsFullIndexPolicy(electronsGrouped, positronsGrouped))) {

        if (particle1.mothersIds()[0] == particle2.mothersIds()[0]) {
          LOG(error) << "Something is wrong here. There should not be dielectrons with same mother.";
        }

        doPair(particle1, particle2, hULS_Mee, hULS_MeePtee, myConfigs.fConfigPtMin, myConfigs.fConfigEtaMax, myConfigs.fConfigApplyDEtaDPhi, myConfigs.fConfigMinDEta, myConfigs.fConfigMinDPhi);

        if (particle1.bQuarkOriginId() < 0 || particle2.bQuarkOriginId() < 0 || particle1.bQuarkOriginId() != particle2.bQuarkOriginId())
          continue;
        doPair(particle1, particle2, hULS_Mee_wPartonicCheck, hULS_MeePtee_wPartonicCheck, myConfigs.fConfigPtMin, myConfigs.fConfigEtaMax, myConfigs.fConfigApplyDEtaDPhi, myConfigs.fConfigMinDEta, myConfigs.fConfigMinDPhi);
      }
      // LS spectrum
      for (auto const& [particle1, particle2] : combinations(o2::soa::CombinationsStrictlyUpperIndexPolicy(electronsGrouped, electronsGrouped))) {

        if (particle1.mothersIds()[0] == particle2.mothersIds()[0]) {
          LOG(error) << "Something is wrong here. There should not be dielectrons with same mother.";
        }

        doPair(particle1, particle2, hLS_Mee, hLS_MeePtee, myConfigs.fConfigPtMin, myConfigs.fConfigEtaMax, myConfigs.fConfigApplyDEtaDPhi, myConfigs.fConfigMinDEta, myConfigs.fConfigMinDPhi);

        if (particle1.bQuarkOriginId() < 0 || particle2.bQuarkOriginId() < 0 || particle1.bQuarkOriginId() != particle2.bQuarkOriginId())
          continue;
        doPair(particle1, particle2, hLS_Mee_wPartonicCheck, hLS_MeePtee_wPartonicCheck, myConfigs.fConfigPtMin, myConfigs.fConfigEtaMax, myConfigs.fConfigApplyDEtaDPhi, myConfigs.fConfigMinDEta, myConfigs.fConfigMinDPhi);
      }
      for (auto const& [particle1, particle2] : combinations(o2::soa::CombinationsStrictlyUpperIndexPolicy(positronsGrouped, positronsGrouped))) {

        if (particle1.mothersIds()[0] == particle2.mothersIds()[0]) {
          LOG(error) << "Something is wrong here. There should not be dielectrons with same mother.";
        }

        doPair(particle1, particle2, hLS_Mee, hLS_MeePtee, myConfigs.fConfigPtMin, myConfigs.fConfigEtaMax, myConfigs.fConfigApplyDEtaDPhi, myConfigs.fConfigMinDEta, myConfigs.fConfigMinDPhi);

        if (particle1.bQuarkOriginId() < 0 || particle2.bQuarkOriginId() < 0 || particle1.bQuarkOriginId() != particle2.bQuarkOriginId())
          continue;
        doPair(particle1, particle2, hLS_Mee_wPartonicCheck, hLS_MeePtee_wPartonicCheck, myConfigs.fConfigPtMin, myConfigs.fConfigEtaMax, myConfigs.fConfigApplyDEtaDPhi, myConfigs.fConfigMinDEta, myConfigs.fConfigMinDPhi);
      }
    }
  }

  void processDummy(aod::McCollisions const&)
  {
    // dummy
  }

  PROCESS_SWITCH(lmeehfcocktailbeauty, processBeauty, "process beauty cocktail", true);
  PROCESS_SWITCH(lmeehfcocktailbeauty, processDummy, "dummy", false);
};

struct lmeehfcocktailcharm {

  HistogramRegistry registry{"registry", {}};

  std::vector<std::shared_ptr<TH1>> hRapQuark;
  std::vector<std::shared_ptr<TH1>> hEta, hPt, hULS_Mee, hLS_Mee;
  std::vector<std::shared_ptr<TH2>> hPtEta, hULS_MeePtee, hLS_MeePtee;
  std::vector<std::shared_ptr<TH1>> hULS_Mee_wPartonicCheck, hLS_Mee_wPartonicCheck;
  std::vector<std::shared_ptr<TH2>> hULS_MeePtee_wPartonicCheck, hLS_MeePtee_wPartonicCheck;

  MyConfigs myConfigs;

  Filter hfFilter = o2::aod::hftable::isHF == static_cast<int>(EFromHFType::kCE);
  using MyFilteredMcParticlesSmeared = soa::Filtered<soa::Join<aod::McParticles, aod::SmearedElectrons, aod::HfTable>>;

  Preslice<MyFilteredMcParticlesSmeared> perCollision = aod::mcparticle::mcCollisionId;

  Partition<MyFilteredMcParticlesSmeared> Electrons = (aod::mcparticle::pdgCode == 11);
  Partition<MyFilteredMcParticlesSmeared> Positrons = (aod::mcparticle::pdgCode == -11);

  void init(o2::framework::InitContext&)
  {
    registry.add<TH1>("NEvents", "NEvents", HistType::kTH1F, {{1, 0, 1}}, false);

    const char* typeNamesSingle = "ce";
    const char* typeTitlesSingle = "c->e";

    AxisSpec rap_axis = {myConfigs.fConfigRapBins, "y_{c}"};
    AxisSpec eta_axis = {myConfigs.fConfigEtaBins, "#eta_{e}"};
    AxisSpec pt_axis = {myConfigs.fConfigPtBins, "#it{p}_{T,e} (GeV/c)"};
    AxisSpec mass_axis = {myConfigs.fConfigMeeBins, "m_{ee} (GeV/c^{2})"};
    AxisSpec ptee_axis = {myConfigs.fConfigPteeBins, "#it{p}_{T,ee} (GeV/c)"};

    // quark histograms
    for (int j = 0; j < Nstages; j++) {
      hRapQuark.push_back(registry.add<TH1>(Form("CharmQuark_Rap_%s_%s", typeNamesSingle, stageNames[j]), Form("Rapidity Charm Quark %s %s", typeTitlesSingle, stageNames[j]), HistType::kTH1F, {rap_axis}, true));
    }

    // single histograms
    for (int j = 0; j < Nstages; j++) {
      hEta.push_back(registry.add<TH1>(Form("Electron_Eta_%s_%s", typeNamesSingle, stageNames[j]), Form("Single Electron Eta %s %s", typeTitlesSingle, stageNames[j]), HistType::kTH1F, {eta_axis}, true));
      hPt.push_back(registry.add<TH1>(Form("Electron_Pt_%s_%s", typeNamesSingle, stageNames[j]), Form("Single Electron Pt %s %s", typeTitlesSingle, stageNames[j]), HistType::kTH1F, {pt_axis}, true));
      hPtEta.push_back(registry.add<TH2>(Form("Electron_PtEta_%s_%s", typeNamesSingle, stageNames[j]), Form("Single Electron Pt vs. Eta %s %s", typeTitlesSingle, stageNames[j]), HistType::kTH2F, {pt_axis, eta_axis}, true));
    }

    // pair histograms
    // ULS
    for (int j = 0; j < Nstages; j++) {
      hULS_Mee.push_back(registry.add<TH1>(Form("ULS_Mee_%s", stageNames[j]), Form("ULS Mee %s", stageNames[j]), HistType::kTH1F, {mass_axis}, true));
      hULS_MeePtee.push_back(registry.add<TH2>(Form("ULS_MeePtee_%s", stageNames[j]), Form("ULS Mee vs Ptee %s", stageNames[j]), HistType::kTH2F, {mass_axis, ptee_axis}, true));

      hULS_Mee_wPartonicCheck.push_back(registry.add<TH1>(Form("ULS_Mee_wPartonicCheck_%s", stageNames[j]), Form("ULS Mee wPartonicCheck %s", stageNames[j]), HistType::kTH1F, {mass_axis}, true));
      hULS_MeePtee_wPartonicCheck.push_back(registry.add<TH2>(Form("ULS_MeePtee_wPartonicCheck_%s", stageNames[j]), Form("ULS Mee vs Ptee wPartonicCheck %s", stageNames[j]), HistType::kTH2F, {mass_axis, ptee_axis}, true));
    }
    // LS
    for (int j = 0; j < Nstages; j++) {
      hLS_Mee.push_back(registry.add<TH1>(Form("LS_Mee_%s", stageNames[j]), Form("LS Mee %s", stageNames[j]), HistType::kTH1F, {mass_axis}, true));
      hLS_MeePtee.push_back(registry.add<TH2>(Form("LS_MeePtee_%s", stageNames[j]), Form("LS Mee vs Ptee %s", stageNames[j]), HistType::kTH2F, {mass_axis, ptee_axis}, true));

      hLS_Mee_wPartonicCheck.push_back(registry.add<TH1>(Form("LS_Mee_wPartonicCheck_%s", stageNames[j]), Form("LS Mee wPartonicCheck %s", stageNames[j]), HistType::kTH1F, {mass_axis}, true));
      hLS_MeePtee_wPartonicCheck.push_back(registry.add<TH2>(Form("LS_MeePtee_wPartonicCheck_%s", stageNames[j]), Form("LS Mee vs Ptee wPartonicCheck %s", stageNames[j]), HistType::kTH2F, {mass_axis, ptee_axis}, true));
    }
  }

  void processCharm(aod::McCollisions const& collisions, MyFilteredMcParticlesSmeared const& mcParticles)
  {
    for (auto const& p : mcParticles) {
      doSingle(p, hEta, hPt, hPtEta, myConfigs.fConfigPtMin, myConfigs.fConfigEtaMax);
      doQuark(p, hRapQuark, myConfigs.fConfigPtMin, myConfigs.fConfigEtaMax, 4);
    }

    for (auto const& collision : collisions) {

      registry.fill(HIST("NEvents"), 0.5);

      auto const electronsGrouped = Electrons->sliceBy(perCollision, collision.globalIndex());
      auto const positronsGrouped = Positrons->sliceBy(perCollision, collision.globalIndex());
      // ULS spectrum
      for (auto const& [particle1, particle2] : combinations(o2::soa::CombinationsFullIndexPolicy(electronsGrouped, positronsGrouped))) {

        if (particle1.mothersIds()[0] == particle2.mothersIds()[0]) {
          LOG(error) << "Something is wrong here. There should not be dielectrons with same mother.";
        }

        doPair(particle1, particle2, hULS_Mee, hULS_MeePtee, myConfigs.fConfigPtMin, myConfigs.fConfigEtaMax, myConfigs.fConfigApplyDEtaDPhi, myConfigs.fConfigMinDEta, myConfigs.fConfigMinDPhi);

        if (particle1.cQuarkOriginId() < 0 || particle2.cQuarkOriginId() < 0 || particle1.cQuarkOriginId() != particle2.cQuarkOriginId())
          continue;
        doPair(particle1, particle2, hULS_Mee_wPartonicCheck, hULS_MeePtee_wPartonicCheck, myConfigs.fConfigPtMin, myConfigs.fConfigEtaMax, myConfigs.fConfigApplyDEtaDPhi, myConfigs.fConfigMinDEta, myConfigs.fConfigMinDPhi);
      }
      // LS
      for (auto const& [particle1, particle2] : combinations(o2::soa::CombinationsStrictlyUpperIndexPolicy(electronsGrouped, electronsGrouped))) {

        if (particle1.mothersIds()[0] == particle2.mothersIds()[0]) {
          LOG(error) << "Something is wrong here. There should not be dielectrons with same mother.";
        }

        doPair(particle1, particle2, hLS_Mee, hLS_MeePtee, myConfigs.fConfigPtMin, myConfigs.fConfigEtaMax, myConfigs.fConfigApplyDEtaDPhi, myConfigs.fConfigMinDEta, myConfigs.fConfigMinDPhi);

        if (particle1.cQuarkOriginId() < 0 || particle2.cQuarkOriginId() < 0 || particle1.cQuarkOriginId() != particle2.cQuarkOriginId())
          continue;
        doPair(particle1, particle2, hLS_Mee_wPartonicCheck, hLS_MeePtee_wPartonicCheck, myConfigs.fConfigPtMin, myConfigs.fConfigEtaMax, myConfigs.fConfigApplyDEtaDPhi, myConfigs.fConfigMinDEta, myConfigs.fConfigMinDPhi);
      }
      for (auto const& [particle1, particle2] : combinations(o2::soa::CombinationsStrictlyUpperIndexPolicy(positronsGrouped, positronsGrouped))) {

        if (particle1.mothersIds()[0] == particle2.mothersIds()[0]) {
          LOG(error) << "Something is wrong here. There should not be dielectrons with same mother.";
        }

        doPair(particle1, particle2, hLS_Mee, hLS_MeePtee, myConfigs.fConfigPtMin, myConfigs.fConfigEtaMax, myConfigs.fConfigApplyDEtaDPhi, myConfigs.fConfigMinDEta, myConfigs.fConfigMinDPhi);

        if (particle1.cQuarkOriginId() < 0 || particle2.cQuarkOriginId() < 0 || particle1.cQuarkOriginId() != particle2.cQuarkOriginId())
          continue;
        doPair(particle1, particle2, hLS_Mee_wPartonicCheck, hLS_MeePtee_wPartonicCheck, myConfigs.fConfigPtMin, myConfigs.fConfigEtaMax, myConfigs.fConfigApplyDEtaDPhi, myConfigs.fConfigMinDEta, myConfigs.fConfigMinDPhi);
      }
    }
  }

  void processDummy(aod::McCollisions const&)
  {
    // dummy
  }

  PROCESS_SWITCH(lmeehfcocktailcharm, processCharm, "process charm cocktail", true);
  PROCESS_SWITCH(lmeehfcocktailcharm, processDummy, "dummy", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<lmeehfcocktailcharm>(cfgc, TaskName("em-lmee-hf-cocktail-charm")),
    adaptAnalysisTask<lmeehfcocktailbeauty>(cfgc, TaskName("em-lmee-hf-cocktail-beauty")),
    adaptAnalysisTask<lmeehfcocktailprefilter>(cfgc, TaskName("em-lmee-hf-cocktail-prefilter"))};
}
