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

#include "Math/Vector4D.h"
#include "Framework/Task.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "PWGDQ/DataModel/ReducedInfoTables.h"
#include "PWGEM/Dilepton/Utils/MCUtilities.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::pwgem::dilepton::mcutil;

using McParticlesSmeared = soa::Join<aod::McParticles, aod::SmearedTracks>;

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
}
DECLARE_SOA_TABLE(HfTable, "AOD", "HFTABLE",
                  hftable::IsHF);
} // namespace o2::aod

const char* stageNames[3] = {"gen", "eff", "eff_and_acc"};

template <typename T>
void doSingle(T& p, std::vector<std::shared_ptr<TH1>> hEta, std::vector<std::shared_ptr<TH1>> hPt, std::vector<std::shared_ptr<TH2>> hPtEta, float ptMin, float etaMax)
{
  float weight[3] = {p.weight(), p.efficiency() * p.weight(), p.efficiency() * p.weight()};
  float pt[3] = {p.pt(), p.ptSmeared(), p.ptSmeared()};
  float eta[3] = {p.eta(), p.etaSmeared(), p.etaSmeared()};
  float cut_pt[3] = {0., 0., ptMin};
  float cut_eta[3] = {9999., 99999., etaMax};
  for (int i = 0; i < 3; i++) {
    if (pt[i] > cut_pt[i] && abs(eta[i]) < cut_eta[i]) {
      hEta[i]->Fill(eta[i], weight[i]);
      hPt[i]->Fill(pt[i], weight[i]);
      hPtEta[i]->Fill(pt[i], eta[i], weight[i]);
    }
  }
}

template <typename T>
void doPair(T& p1, T& p2, std::vector<std::shared_ptr<TH1>> hMee, std::vector<std::shared_ptr<TH2>> hMeePtee, float ptMin, float etaMax)
{

  ROOT::Math::PtEtaPhiMVector v1(p1.ptSmeared(), p1.etaSmeared(), p1.phiSmeared(), o2::constants::physics::MassElectron);
  ROOT::Math::PtEtaPhiMVector v2(p2.ptSmeared(), p2.etaSmeared(), p2.phiSmeared(), o2::constants::physics::MassElectron);
  ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
  ROOT::Math::PtEtaPhiMVector v1_gen(p1.pt(), p1.eta(), p1.phi(), o2::constants::physics::MassElectron);
  ROOT::Math::PtEtaPhiMVector v2_gen(p2.pt(), p2.eta(), p2.phi(), o2::constants::physics::MassElectron);
  ROOT::Math::PtEtaPhiMVector v12_gen = v1_gen + v2_gen;

  double mass[3] = {v12_gen.M(), v12.M(), v12.M()};
  double pt[3] = {v12_gen.Pt(), v12.Pt(), v12.Pt()};
  float pt1[3] = {p1.pt(), p1.ptSmeared(), p1.ptSmeared()};
  float pt2[3] = {p2.pt(), p2.ptSmeared(), p2.ptSmeared()};
  float eta1[3] = {p1.eta(), p1.etaSmeared(), p1.etaSmeared()};
  float eta2[3] = {p2.eta(), p2.etaSmeared(), p2.etaSmeared()};
  float weight[3] = {p1.weight() * p2.weight(), p1.efficiency() * p2.efficiency() * p1.weight() * p2.weight(), p1.efficiency() * p2.efficiency() * p1.weight() * p2.weight()};
  float cut_pt[3] = {0., 0., ptMin};
  float cut_eta[3] = {9999., 99999., etaMax};

  for (int i = 0; i < 3; i++) {
    if (pt1[i] > cut_pt[i] && pt2[i] > cut_pt[i] && abs(eta1[i]) < cut_eta[i] && abs(eta2[i]) < cut_eta[i]) {
      hMee[i]->Fill(mass[i], weight[i]);
      hMeePtee[i]->Fill(mass[i], pt[i], weight[i]);
    }
  }
}

struct lmeehfcocktailprefilter {

  Produces<o2::aod::HfTable> hfTable;
  void process(aod::McParticles const& mcParticles)
  {
    for (auto const& p : mcParticles) {

      if (abs(p.pdgCode()) != 11 || o2::mcgenstatus::getHepMCStatusCode(p.statusCode()) != 1) {
        hfTable(EFromHFType::kNoE);
        continue;
      }

      // no HF = 0; c->e = 1; b->e = 2; b->c->e = 3;
      int isHF = 0;
      if (IsFromCharm(p, mcParticles) > -1) {
        isHF = isHF + 1;
      }
      if (IsFromBeauty(p, mcParticles) > -1) {
        isHF = isHF + 2;
      }
      hfTable(isHF);
    }
  }
};

struct lmeehfcocktailbeauty {

  HistogramRegistry registry{"registry", {}};

  std::vector<std::vector<std::shared_ptr<TH1>>> hEta, hPt, hULS_Mee;
  std::vector<std::vector<std::shared_ptr<TH2>>> hPtEta, hULS_MeePtee;
  std::vector<std::shared_ptr<TH1>> hLSpp_Mee, hLSmm_Mee;
  std::vector<std::shared_ptr<TH2>> hLSpp_MeePtee, hLSmm_MeePtee;

  Configurable<float> fConfigPtMin{"cfgPtMin", 0.2, "min. pT"};
  Configurable<float> fConfigEtaMax{"cfgEtaMax", 0.8, "max. |eta|"};
  ConfigurableAxis fConfigPtBins{"cfgPtBins", {200, 0.f, 10.f}, "pT binning"};
  ConfigurableAxis fConfigEtaBins{"cfgEtaBins", {200, -10.f, 10.f}, "eta binning"};
  ConfigurableAxis fConfigMeeBins{"cfgMeeBins", {800, 0.f, 8.f}, "Mee binning"};
  ConfigurableAxis fConfigPteeBins{"cfgPteeBins", {400, 0.f, 10.f}, "pTee binning"};

  Filter hfFilter = o2::aod::hftable::isHF == static_cast<int>(EFromHFType::kBE) || o2::aod::hftable::isHF == static_cast<int>(EFromHFType::kBCE);
  using MyFilteredMcParticlesSmeared = soa::Filtered<soa::Join<aod::McParticles, aod::SmearedTracks, aod::HfTable>>;

  Preslice<MyFilteredMcParticlesSmeared> perCollision = aod::mcparticle::mcCollisionId;

  Partition<MyFilteredMcParticlesSmeared> Electrons = (aod::mcparticle::pdgCode == 11);
  Partition<MyFilteredMcParticlesSmeared> Positrons = (aod::mcparticle::pdgCode == -11);

  void init(o2::framework::InitContext&)
  {
    registry.add<TH1>("NEvents", "NEvents", HistType::kTH1F, {{1, 0, 1}}, false);

    const char* typeNamesPairULS[4] = {"Be_Be", "BCe_BCe", "BCe_Be_SameB", "allB"};
    const char* typeNamesPairLS = "BCe_Be_DiffB";
    const char* typeNamesSingle[2] = {"be", "bce"};
    const char* typeTitlesSingle[2] = {"b->e", "b->c->e"};

    AxisSpec eta_axis = {fConfigEtaBins, "#eta"};
    AxisSpec pt_axis = {fConfigPtBins, "#it{p}_{T} (GeV/c)"};
    AxisSpec mass_axis = {fConfigMeeBins, "m_{ee} (GeV/c^{2})"};
    AxisSpec ptee_axis = {fConfigPteeBins, "#it{p}_{T,ee} (GeV/c)"};

    // single histograms
    for (int i = 0; i < 2; i++) {
      std::vector<std::shared_ptr<TH1>> hEta_temp, hPt_temp;
      std::vector<std::shared_ptr<TH2>> hPtEta_temp;
      for (int j = 0; j < 3; j++) {
        hEta_temp.push_back(registry.add<TH1>(Form("Eta_%s_%s", typeNamesSingle[i], stageNames[j]), Form("Eta %s %s", typeTitlesSingle[i], stageNames[j]), HistType::kTH1F, {eta_axis}, true));
        hPt_temp.push_back(registry.add<TH1>(Form("Pt_%s_%s", typeNamesSingle[i], stageNames[j]), Form("Pt %s %s", typeTitlesSingle[i], stageNames[j]), HistType::kTH1F, {pt_axis}, true));
        hPtEta_temp.push_back(registry.add<TH2>(Form("PtEta_%s_%s", typeNamesSingle[i], stageNames[j]), Form("Pt vs. Eta %s %s", typeTitlesSingle[i], stageNames[j]), HistType::kTH2F, {pt_axis, eta_axis}, true));
      }
      hEta.push_back(hEta_temp);
      hPt.push_back(hPt_temp);
      hPtEta.push_back(hPtEta_temp);
    }

    // pair histograms
    // ULS
    for (int i = 0; i < 4; i++) {
      std::vector<std::shared_ptr<TH1>> hMee_temp;
      std::vector<std::shared_ptr<TH2>> hMeePtee_temp;
      for (int j = 0; j < 3; j++) {
        hMee_temp.push_back(registry.add<TH1>(Form("ULS_Mee_%s_%s", typeNamesPairULS[i], stageNames[j]), Form("ULS Mee %s %s", typeNamesPairULS[i], stageNames[j]), HistType::kTH1F, {mass_axis}, true));
        hMeePtee_temp.push_back(registry.add<TH2>(Form("ULS_MeePtee_%s_%s", typeNamesPairULS[i], stageNames[j]), Form("ULS Mee vs. Ptee %s %s", typeNamesPairULS[i], stageNames[j]), HistType::kTH2F, {mass_axis, ptee_axis}, true));
      }
      hULS_Mee.push_back(hMee_temp);
      hULS_MeePtee.push_back(hMeePtee_temp);
    }
    // LS
    for (int j = 0; j < 3; j++) {
      hLSpp_Mee.push_back(registry.add<TH1>(Form("LSpp_Mee_%s_%s", typeNamesPairLS, stageNames[j]), Form("LS++ Mee %s %s", typeNamesPairLS, stageNames[j]), HistType::kTH1F, {mass_axis}, true));
      hLSmm_Mee.push_back(registry.add<TH1>(Form("LSmm_Mee_%s_%s", typeNamesPairLS, stageNames[j]), Form("LS-- Mee %s %s", typeNamesPairLS, stageNames[j]), HistType::kTH1F, {mass_axis}, true));
      hLSpp_MeePtee.push_back(registry.add<TH2>(Form("LSpp_MeePtee_%s_%s", typeNamesPairLS, stageNames[j]), Form("LS++ Mee vs. Ptee %s %s", typeNamesPairLS, stageNames[j]), HistType::kTH2F, {mass_axis, ptee_axis}, true));
      hLSmm_MeePtee.push_back(registry.add<TH2>(Form("LSmm_MeePtee_%s_%s", typeNamesPairLS, stageNames[j]), Form("LS-- Mee vs. Ptee %s %s", typeNamesPairLS, stageNames[j]), HistType::kTH2F, {mass_axis, ptee_axis}, true));
    }
  }

  void processBeauty(aod::McCollisions const& collisions, MyFilteredMcParticlesSmeared const& mcParticles, aod::McParticles const& mcParticlesAll)
  {
    for (auto const& collision : collisions) {

      registry.fill(HIST("NEvents"), 0.5);

      for (auto const& p : mcParticles) {
        int from_quark = p.isHF() - 2;
        doSingle(p, hEta[from_quark], hPt[from_quark], hPtEta[from_quark], fConfigPtMin, fConfigEtaMax);
      }

      auto const electronsGrouped = Electrons->sliceBy(perCollision, collision.globalIndex());
      auto const positronsGrouped = Positrons->sliceBy(perCollision, collision.globalIndex());
      // ULS spectrum
      for (auto const& [particle1, particle2] : combinations(o2::soa::CombinationsFullIndexPolicy(electronsGrouped, positronsGrouped))) {
        int type = IsHF(particle1, particle2, mcParticlesAll);
        if (type == static_cast<int>(EM_HFeeType::kUndef))
          continue;
        if ((type < static_cast<int>(EM_HFeeType::kBe_Be)) || (type > static_cast<int>(EM_HFeeType::kBCe_Be_SameB))) {
          LOG(error) << "Something is wrong here. There should only be pairs of type kBe_Be = 1, kBCe_BCe = 2 and kBCe_Be_SameB = 3 left at this point.";
        }
        doPair(particle1, particle2, hULS_Mee[type - 1], hULS_MeePtee[type - 1], fConfigPtMin, fConfigEtaMax);
        doPair(particle1, particle2, hULS_Mee[3], hULS_MeePtee[3], fConfigPtMin, fConfigEtaMax); // fill the 'allB' histograms that holds the sum of the others
      }
      // LS spectrum
      for (auto const& [particle1, particle2] : combinations(o2::soa::CombinationsStrictlyUpperIndexPolicy(electronsGrouped, electronsGrouped))) {
        int type = IsHF(particle1, particle2, mcParticlesAll);
        if (type == static_cast<int>(EM_HFeeType::kUndef))
          continue;
        if (type != static_cast<int>(EM_HFeeType::kBCe_Be_DiffB)) {
          LOG(error) << "Something is wrong here. There should only be pairs of type kBCe_Be_DiffB = 4 left at this point.";
        }
        doPair(particle1, particle2, hLSmm_Mee, hLSmm_MeePtee, fConfigPtMin, fConfigEtaMax);
      }
      for (auto const& [particle1, particle2] : combinations(o2::soa::CombinationsStrictlyUpperIndexPolicy(positronsGrouped, positronsGrouped))) {
        int type = IsHF(particle1, particle2, mcParticlesAll);
        if (type == static_cast<int>(EM_HFeeType::kUndef))
          continue;
        if (type != static_cast<int>(EM_HFeeType::kBCe_Be_DiffB)) {
          LOG(error) << "Something is wrong here. There should only be pairs of type kBCe_Be_DiffB = 4 left at this point.";
        }
        doPair(particle1, particle2, hLSpp_Mee, hLSpp_MeePtee, fConfigPtMin, fConfigEtaMax);
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

  std::vector<std::shared_ptr<TH1>> hEta, hPt, hULS_Mee;
  std::vector<std::shared_ptr<TH2>> hPtEta, hULS_MeePtee;

  Configurable<float> fConfigPtMin{"cfgPtMin", 0.2f, "min. pT"};
  Configurable<float> fConfigEtaMax{"cfgEtaMax", 0.8f, "max. |eta|"};
  ConfigurableAxis fConfigPtBins{"cfgPtBins", {200, 0.f, 10.f}, "pT binning"};
  ConfigurableAxis fConfigEtaBins{"cfgEtaBins", {200, -10.f, 10.f}, "eta binning"};
  ConfigurableAxis fConfigMeeBins{"cfgMeeBins", {800, 0.f, 8.f}, "Mee binning"};
  ConfigurableAxis fConfigPteeBins{"cfgPteeBins", {400, 0.f, 10.f}, "pTee binning"};

  Filter hfFilter = o2::aod::hftable::isHF == static_cast<int>(EFromHFType::kCE);
  using MyFilteredMcParticlesSmeared = soa::Filtered<soa::Join<aod::McParticles, aod::SmearedTracks, aod::HfTable>>;

  Preslice<MyFilteredMcParticlesSmeared> perCollision = aod::mcparticle::mcCollisionId;

  Partition<MyFilteredMcParticlesSmeared> Electrons = (aod::mcparticle::pdgCode == 11);
  Partition<MyFilteredMcParticlesSmeared> Positrons = (aod::mcparticle::pdgCode == -11);

  void init(o2::framework::InitContext&)
  {
    registry.add<TH1>("NEvents", "NEvents", HistType::kTH1F, {{1, 0, 1}}, false);

    const char* typeNamesPairULS = "Ce_Ce";
    const char* typeNamesSingle = "ce";
    const char* typeTitlesSingle = "c->e";

    AxisSpec eta_axis = {fConfigEtaBins, "#eta"};
    AxisSpec pt_axis = {fConfigPtBins, "#it{p}_{T} (GeV/c)"};
    AxisSpec mass_axis = {fConfigMeeBins, "m_{ee} (GeV/c^{2})"};
    AxisSpec ptee_axis = {fConfigPteeBins, "#it{p}_{T,ee} (GeV/c)"};

    // single histograms
    for (int j = 0; j < 3; j++) {
      hEta.push_back(registry.add<TH1>(Form("Eta_%s_%s", typeNamesSingle, stageNames[j]), Form("Eta %s %s", typeTitlesSingle, stageNames[j]), HistType::kTH1F, {eta_axis}, true));
      hPt.push_back(registry.add<TH1>(Form("Pt_%s_%s", typeNamesSingle, stageNames[j]), Form("Pt %s %s", typeTitlesSingle, stageNames[j]), HistType::kTH1F, {pt_axis}, true));
      hPtEta.push_back(registry.add<TH2>(Form("PtEta_%s_%s", typeNamesSingle, stageNames[j]), Form("Pt vs. Eta %s %s", typeTitlesSingle, stageNames[j]), HistType::kTH2F, {pt_axis, eta_axis}, true));
    }

    // pair histograms
    // ULS
    for (int j = 0; j < 3; j++) {
      hULS_Mee.push_back(registry.add<TH1>(Form("ULS_Mee_%s_%s", typeNamesPairULS, stageNames[j]), Form("ULS Mee %s %s", typeNamesPairULS, stageNames[j]), HistType::kTH1F, {mass_axis}, true));
      hULS_MeePtee.push_back(registry.add<TH2>(Form("ULS_MeePtee_%s_%s", typeNamesPairULS, stageNames[j]), Form("ULS MeePtee %s %s", typeNamesPairULS, stageNames[j]), HistType::kTH2F, {mass_axis, ptee_axis}, true));
    }
  }

  void processCharm(aod::McCollisions const& collisions, MyFilteredMcParticlesSmeared const& mcParticles, aod::McParticles const& mcParticlesAll)
  {
    for (auto const& collision : collisions) {

      registry.fill(HIST("NEvents"), 0.5);

      for (auto const& p : mcParticles) {
        doSingle(p, hEta, hPt, hPtEta, fConfigPtMin, fConfigEtaMax);
      }

      auto const electronsGrouped = Electrons->sliceBy(perCollision, collision.globalIndex());
      auto const positronsGrouped = Positrons->sliceBy(perCollision, collision.globalIndex());
      // ULS spectrum
      for (auto const& [particle1, particle2] : combinations(o2::soa::CombinationsFullIndexPolicy(electronsGrouped, positronsGrouped))) {
        int type = IsHF(particle1, particle2, mcParticlesAll);
        if (type == static_cast<int>(EM_HFeeType::kUndef))
          continue;
        if (type != static_cast<int>(EM_HFeeType::kCe_Ce)) {
          LOG(error) << "Something is wrong here. There should only be pairs of type kCe_Ce = 0 left at this point.";
        }
        doPair(particle1, particle2, hULS_Mee, hULS_MeePtee, fConfigPtMin, fConfigEtaMax);
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
