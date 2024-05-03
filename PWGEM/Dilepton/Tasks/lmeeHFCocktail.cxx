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
// Analysis task for lmee heavy flavour cocktail

#include "Math/Vector4D.h"
#include "Framework/Task.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "PWGDQ/DataModel/ReducedInfoTables.h"
#include  "PWGEM/Dilepton/Utils/MCUtilities.h"


using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using McParticlesSmeared = soa::Join<aod::McParticles, aod::SmearedTracks>;

namespace o2::aod
{
namespace hftable // reconstructed collision information
{
DECLARE_SOA_COLUMN(IsHF, isHF, int);
}
DECLARE_SOA_TABLE(HfTable, "AOD", "HFTABLE" ,
  hftable::IsHF);
}

struct lmeehfcocktailprefilter {

  Produces<o2::aod::HfTable> hfTable;
  void process(aod::McParticles const& mcParticles){
    for (auto const& p : mcParticles){

      if (abs(p.pdgCode()) != 11 || o2::mcgenstatus::getHepMCStatusCode(p.statusCode()) != 1){
        hfTable(-1);
        continue;
      }

      // no HF = 0; c->e = 1; b->e = 2; b->c->e = 3;
      int isHF = 0;
      if (o2::aod::pwgem::dilepton::mcutil::IsFromCharm(p, mcParticles) > -1){
        isHF = isHF + 1;
      }
      if (o2::aod::pwgem::dilepton::mcutil::IsFromBeauty(p, mcParticles) > -1){
        isHF = isHF + 2;
      }
      hfTable(isHF);
    }
  }
};

struct lmeehfcocktail {

  HistogramRegistry registry{"registry", {}};
  std::vector<std::shared_ptr<TH1>> hType, hMee_allB, hEta, hPt;
  std::vector<std::vector<std::shared_ptr<TH1>>> hMee;

  Filter hfFilter = o2::aod::hftable::isHF > 0;
  using MyFilteredMcParticlesSmeared = soa::Filtered<soa::Join<aod::McParticles, aod::SmearedTracks, aod::HfTable>>;

  Preslice<MyFilteredMcParticlesSmeared> perCollision = aod::mcparticle::mcCollisionId;

  Partition<MyFilteredMcParticlesSmeared> Electrons = (aod::mcparticle::pdgCode == 11);
  Partition<MyFilteredMcParticlesSmeared> Positrons = (aod::mcparticle::pdgCode == -11);

  template <typename T, typename U>
  void doPair(T &p1, T &p2, U &allParticles , int spectrum_type){
    int type = o2::aod::pwgem::dilepton::mcutil::IsHF(p1, p2, allParticles);
    hType[spectrum_type]->Fill(type);

    ROOT::Math::PtEtaPhiMVector v1(p1.ptSmeared(), p1.etaSmeared(), p1.phiSmeared(), o2::constants::physics::MassElectron);
    ROOT::Math::PtEtaPhiMVector v2(p2.ptSmeared(), p2.etaSmeared(), p2.phiSmeared(), o2::constants::physics::MassElectron);
    ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
    double mass = v12.M();
    hMee[spectrum_type][type+1]->Fill(mass);
    if (type>0){
      hMee_allB[spectrum_type]->Fill(mass);
    }
  }

  template <typename T>
  void doSingle(T &p){
    hEta[p.isHF()-1]->Fill(p.etaSmeared());
    hPt[p.isHF()-1]->Fill(p.ptSmeared());
    registry.fill(HIST("TypeQ"), p.isHF());
  }

  void init(o2::framework::InitContext&)
  {
    registry.add<TH1>("NEvents", "NEvents", HistType::kTH1F, {{1, 0, 1}}, false);
    AxisSpec type_axis = {6, -1.5, 4.5, "type"};
    AxisSpec mass_axis = {1000, 0., 10., "inv. mass (GeV/c^{2})"};
    AxisSpec eta_axis = {100,-6.,6., "#{eta}"};
    AxisSpec pt_axis = {100,0.,10., "p_T (GeV/c)"};
    AxisSpec typeQ_axis = {3,0.5,3.4,"type"};

    const char* typeNamesQ[3] = {"c", "b", "bc"};
    registry.add<TH1>(Form("TypeQ"),Form("TypeQ"), HistType::kTH1F, {eta_axis}, true);
    for(int i=0; i<3; i++){
      hEta.push_back(registry.add<TH1>(Form("Eta_%se",typeNamesQ[i]),Form("Eta_%se",typeNamesQ[i]), HistType::kTH1F, {eta_axis}, true));
      hPt.push_back(registry.add<TH1>(Form("Pt_%se",typeNamesQ[i]),Form("Pt_%se",typeNamesQ[i]), HistType::kTH1F, {pt_axis}, true));
      registry.get<TH1>(HIST("TypeQ"))->GetXaxis()->SetBinLabel(i + 1, typeNamesQ[i]);
    }
  

    const char* typeNames[6] = {"Undef", "Ce_Ce", "Be_Be", "BCe_BCe", "BCe_Be_SameB", "BCe_Be_DiffB"};
    for (int i=0; i<2; i++){
      std::string name = "ULS";
      if (i==1){
        name="LS";
      }
      hType.push_back(registry.add<TH1>(Form("%s_Type", name.c_str()), Form("%s_Type", name.c_str()), HistType::kTH1F, {type_axis}, true));
      for (int j=0; j<6; j++){
        hType[i]->GetXaxis()->SetBinLabel(j+1, typeNames[j]);
      }
      std::vector<std::shared_ptr<TH1>> hMee_temp;
      std::vector<std::shared_ptr<TH1>> hMee_MC_temp;
      for (int j=0; j<6; j++){
        hMee_temp.push_back(registry.add<TH1>(Form("%s_%s_Mee", name.c_str(), typeNames[j]), Form("%s_%s_Mee", name.c_str(), typeNames[j]), HistType::kTH1F, {mass_axis}, true));
      }
      hMee.push_back(hMee_temp);
      hMee_allB.push_back(registry.add<TH1>(Form("%s_allB_Mee", name.c_str()), Form("%s_allB_Mee", name.c_str()), HistType::kTH1F, {mass_axis}, true));
    }
  }

  void process(aod::McCollisions const& collisions, MyFilteredMcParticlesSmeared const& mcParticles, aod::McParticles const& mcParticlesAll){
    for (auto const& collision : collisions){
      registry.fill(HIST("NEvents"), 0.5);

      for (auto const& p : mcParticles){
        doSingle(p);
      }
      auto const electronsGrouped = Electrons->sliceBy(perCollision, collision.globalIndex());
      auto const positronsGrouped = Positrons->sliceBy(perCollision, collision.globalIndex());
      // ULS spectrum
      for (auto const& [particle1, particle2] : combinations(o2::soa::CombinationsFullIndexPolicy(electronsGrouped, positronsGrouped))) {
        doPair(particle1, particle2, mcParticlesAll, 0);
      }
      // LS spectrum
      for (auto const& [particle1, particle2] : combinations(o2::soa::CombinationsStrictlyUpperIndexPolicy(electronsGrouped, electronsGrouped))) {
        doPair(particle1, particle2, mcParticlesAll, 1);
      }
      for (auto const& [particle1, particle2] : combinations(o2::soa::CombinationsStrictlyUpperIndexPolicy(positronsGrouped, positronsGrouped))) {
        doPair(particle1, particle2, mcParticlesAll, 1);
      }
    }
  }

};








WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<lmeehfcocktail>(cfgc, TaskName("em-lmee-hf-cocktail")),
    adaptAnalysisTask<lmeehfcocktailprefilter>(cfgc, TaskName("em-lmee-hf-cocktail-prefilter"))
  };
}

