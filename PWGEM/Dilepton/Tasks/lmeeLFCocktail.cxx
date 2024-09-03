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
/// \file lmeeLFCocktail.cxx
/// \analysis task for lmee light flavour cocktail
/// \author Daniel Samitz, <daniel.samitz@cern.ch>, SMI Vienna

#include "Math/Vector4D.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "PWGDQ/DataModel/ReducedInfoTables.h"

using namespace o2;
using namespace o2::framework;

using McParticlesSmeared = soa::Join<aod::McParticles, aod::SmearedTracks>;

struct lmeelfcocktail {

  struct mesonInfo {
    TString name;
    std::vector<int> decayModes;
  };

  std::map<int, TString> decays = {
    {-1, "e_e/"},
    {-2, "e_e_e_e/"},
    {22, "gamma_e_e/"},
    {223, "omega_e_e/"},
    {211 * 211, "pi_pi_e_e/"},
    {111, "pi0_e_e/"},
    {221, "eta_e_e/"}};

  std::map<int, mesonInfo> mesons = {
    {111, {"pi0/", {22, -2}}},
    {221, {"eta/", {22, 211 * 211, -2}}},
    {331, {"etaP/", {22, 223}}},
    {113, {"rho/", {-1}}},
    {223, {"omega/", {-1, 111}}},
    {333, {"phi/", {-1, 111, 221}}}};

  std::map<TString, int> histogramId;

  std::vector<TString> stage = {"gen/", "rec/"};

  HistogramRegistry registry{"registry", {}};

  Configurable<float> fConfigMaxEta{"cfgMaxEta", 0.8, "maxium |eta|"};
  Configurable<float> fConfigMinPt{"cfgMinPt", 0.2, "minium pT"};
  Configurable<float> fConfigMaxPt{"cfgMaxPt", 8.0, "maximum pT"};
  Configurable<float> fConfigMinOpAng{"cfgMinOpAng", 0.050, "minimum opening angle"};
  Configurable<float> fConfigMinPtee{"cfgMinPtee", 0.0, "minimum pair pT"};
  ConfigurableAxis fConfigMeeBins{"cfgMeeBins", {800, 0.f, 8.f}, "Mee binning"};
  ConfigurableAxis fConfigPteeBins{"cfgPteeBins", {400, 0.f, 10.f}, "pTee binning"};
  ConfigurableAxis fConfigPtBins{"cfgPtBins", {200, 0.f, 10.f}, "pT binning"};
  ConfigurableAxis fConfigEtaBins{"cfgEtaBins", {200, -10.f, 10.f}, "eta binning"};
  ConfigurableAxis fConfigPhiBins{"cfgPhiBins", {200, -TMath::Pi(), TMath::Pi()}, "phi binning"};
  ConfigurableAxis fConfigPhiVBins{"cfgPhiVBins", {200, -TMath::Pi(), TMath::Pi()}, "phiV binning"};
  ConfigurableAxis fConfigOpAngBins{"cfgOpAngBins", {200, -TMath::Pi(), TMath::Pi()}, "opening angle binning"};

  std::vector<std::shared_ptr<TH1>> histograms1D;
  std::vector<std::shared_ptr<TH2>> histograms2D;


  template <typename T, typename U>
  bool from_primary(T& p1, U& mcParticles)
  {
    if (!p1.has_mothers()){
      return false;
    }
    auto mother = mcParticles.iteratorAt(p1.mothersIds()[0]);
    if (mother.has_mothers()){
      return false;
    }
    return true;
  }

  bool isAcceptedSingle(ROOT::Math::PtEtaPhiMVector p1){
    if (p1.Pt() < fConfigMinPt)
      return false;
    if (p1.Pt() > fConfigMaxPt)
      return false;
    if (abs(p1.Eta()) > fConfigMaxEta)
      return false;
    return true;
  }


  bool isAcceptedPair(ROOT::Math::PtEtaPhiMVector p1, ROOT::Math::PtEtaPhiMVector p2){
    if (!isAcceptedSingle(p1)){
      return false;
    }
    if (!isAcceptedSingle(p2)){
      return false;
    }
    ROOT::Math::PtEtaPhiMVector p12 = p1+p2;
    if (p12.Pt() < fConfigMinPtee)
      return false;
    if (TMath::ACos(p1.Vect().Unit().Dot(p2.Vect().Unit())) < fConfigMinOpAng)
      return false;
    return true;
  }

  template <typename T>
  bool isAcceptedPair(T& p1, T& p2)
  {
    ROOT::Math::PtEtaPhiMVector v1(p1.ptSmeared(), p1.etaSmeared(), p1.phiSmeared(), o2::constants::physics::MassElectron);
    ROOT::Math::PtEtaPhiMVector v2(p2.ptSmeared(), p2.etaSmeared(), p2.phiSmeared(), o2::constants::physics::MassElectron);
    return isAcceptedPair(v1,v2);
  }

  template <typename T>
  bool isAcceptedSingle(T& p1)
  {
    ROOT::Math::PtEtaPhiMVector v1(p1.ptSmeared(), p1.etaSmeared(), p1.phiSmeared(), o2::constants::physics::MassElectron);
    return isAcceptedSingle(v1);
  }


  void addHistogram1D(TString histname, AxisSpec axis, int& i)
  {
    for (auto s : stage) {
      i++;
      TString name = s + histname;
      histogramId[name] = i;
      histograms1D.push_back(registry.add<TH1>(name, histname, HistType::kTH1F, {axis}, true));
      for (auto const& [pdg, meson] : mesons) {
        i++;
        name = s + meson.name + histname;
        histogramId[name] = i;
        histograms1D.push_back(registry.add<TH1>(name, histname, HistType::kTH1F, {axis}, true));
        for (auto const& mode : meson.decayModes) {
          i++;
          name = s + meson.name + decays[mode] + histname;
          histogramId[name] = i;
          histograms1D.push_back(registry.add<TH1>(name, histname, HistType::kTH1F, {axis}, true));
        }
      }
    }
  }

  void addHistogram1D_mother(TString histname, AxisSpec axis, int& i)
  {
    for (auto const& [pdg, meson] : mesons) {
      i++;
      TString name = stage[0] + meson.name + histname;
      histogramId[name] = i;
      histograms1D.push_back(registry.add<TH1>(name, histname, HistType::kTH1F, {axis}, true));
    }
  }

  void addHistogram2D(TString histname, AxisSpec axis1, AxisSpec axis2, int& i)
  {
    for (auto s : stage) {
      i++;
      TString name = s + histname;
      histogramId[name] = i;
      histograms2D.push_back(registry.add<TH2>(name, histname, HistType::kTH2F, {axis1, axis2}, true));
      for (auto const& [pdg, meson] : mesons) {
        i++;
        name = s + meson.name + histname;
        histogramId[name] = i;
        histograms2D.push_back(registry.add<TH2>(name, histname, HistType::kTH2F, {axis1, axis2}, true));
        for (auto const& mode : meson.decayModes) {
          i++;
          name = s + meson.name + decays[mode] + histname;
          histogramId[name] = i;
          histograms2D.push_back(registry.add<TH2>(name, histname, HistType::kTH2F, {axis1, axis2}, true));
        }
      }
    }
  }

  void fillHistogram1D(TString histname, int s, int pdg, int other_daughter_pdg, float value, float weight){
    histograms1D[histogramId[stage[s]+histname]]->Fill(value,weight);
    histograms1D[histogramId[stage[s]+mesons[pdg].name+histname]]->Fill(value,weight);
    histograms1D[histogramId[stage[s]+mesons[pdg].name+decays[other_daughter_pdg]+histname]]->Fill(value,weight);
  }

  void fillHistogram2D(TString histname, int s, int pdg, int other_daughter_pdg, float value1, float value2, float weight){
    histograms1D[histogramId[stage[s]+histname]]->Fill(value1, value2 ,weight);
    histograms1D[histogramId[stage[s]+mesons[pdg].name+histname]]->Fill(value1, value2, ,weight);
    histograms1D[histogramId[stage[s]+mesons[pdg].name+decays[other_daughter_pdg]+histname]]->Fill(value1 ,value2, weight);
  }

  void init(InitContext& context)
  {
    AxisSpec mass_axis = {fConfigMeeBins, "m_{ee} (GeV/#it{c}^{2})"};
    AxisSpec ptee_axis = {fConfigPteeBins, "#it{p}_{T,ee} (GeV/#it{c})"};
    AxisSpec eta_axis = {fConfigEtaBins, "#it{#eta}_{e}"};
    AxisSpec pt_axis = {fConfigPtBins, "#it{p}_{T,e} (GeV/c)"};
    AxisSpec phi_axis = {fConfigPhiBins, "#it{#varphi}_{e}"};
    AxisSpec phiV_axis = {fConfigPhiVBins, "#it{#varphi}_{V,ee}"};
    AxisSpec opAng_axis = {fConfigOpAngBins, "#it{#omega}_{ee}"};
    AxisSpec eta_axis_mother = {fConfigEtaBins, "#it{#eta}_{mother}"};
    AxisSpec pt_axis_mother = {fConfigPtBins, "#it{p}_{T,mother} (GeV/#it{c})"};
    AxisSpec phi_axis_mother = {fConfigPhiBins, "#it{varphi}_{mother}"};

    if (context.mOptions.get<bool>("processPairing")) {
      registry.add<TH2>("gen/ULS", "ULS gen.", HistType::kTH2F, {mass_axis, ptee_axis}, true);
      registry.add<TH2>("gen/LSpp", "LS++ gen.", HistType::kTH2F, {mass_axis, ptee_axis}, true);
      registry.add<TH2>("gen/LSmm", "LS-- gen.", HistType::kTH2F, {mass_axis, ptee_axis}, true);
      registry.add<TH2>("rec/ULS", "ULS rec.", HistType::kTH2F, {mass_axis, ptee_axis}, true);
      registry.add<TH2>("rec/LSpp", "LS++ rec.", HistType::kTH2F, {mass_axis, ptee_axis}, true);
      registry.add<TH2>("rec/LSmm", "LS-- rec.", HistType::kTH2F, {mass_axis, ptee_axis}, true);
    }
    if (context.mOptions.get<bool>("processCocktail")) {
      int i = -1;
      addHistogram1D("Pt", pt_axis, i);
      addHistogram1D("Eta", eta_axis, i);
      addHistogram1D("Phi", phi_axis, i);
      addHistogram1D_mother("Mohter_Pt", pt_axis_mother, i);
      addHistogram1D_mother("Mother_Eta", eta_axis_mother, i);
      addHistogram1D_mother("Mother_Phi", phi_axis_mother, i);
      addHistogram1D("PhiV", phiV_axis, i);
      addHistogram1D("OpAng", opAng_axis, i);
      addHistogram1D("Mee", mass_axis, i);
      addHistogram1D("Ptee", ptee_axis, i);
      i = -1;
      addHistogram2D("MeeVsPtee", mass_axis, ptee_axis, i);
      // missing DCA
    }
  }

  void processCocktail(McParticlesSmeared const& mcParticles)
  {
    for (auto const& particle : mcParticles) {
      if (particle.has_mothers()) {
        continue;
      }
      int pdg = abs(particle.pdgCode());
      if (mesons.find(pdg) == mesons.end()) {
        LOG(error) << "Found mother particle with pdg = " << pdg << " that is not in list of mesons";
      }
      if (!particle.has_daughters()) {
        LOG(error) << "Found meson that has no daughters";
      }

      int other_daughter_pdg = -1;
      int nEle = 0;
      int nPos = 0;

      ROOT::Math::PtEtaPhiMVector pEle[2], pPos[2];
      float wEle[2], wPos[2];
      for (const auto& daughter : particle.daughters_as<McParticlesSmeared>()) {
        int temp_pdg = daughter.pdgCode();
        if (temp_pdg == 11) {
          ROOT::Math::PtEtaPhiMVector temp_p_gen(daughter.pt(), daughter.eta(), daughter.phi(), o2::constants::physics::MassElectron);
          ROOT::Math::PtEtaPhiMVector temp_p(daughter.ptSmeared(), daughter.etaSmeared(), daughter.phiSmeared(), o2::constants::physics::MassElectron);
          pEle[0] = temp_p_gen;
          pEle[1] = temp_p;
          wEle[0] = daughter.weight();
          wEle[1] = daughter.weight()*daughter.efficiency();
          nEle++;
          continue;
        }
        if (temp_pdg == -11) {
          ROOT::Math::PtEtaPhiMVector temp_p_gen(daughter.pt(), daughter.eta(), daughter.phi(), o2::constants::physics::MassElectron);
          ROOT::Math::PtEtaPhiMVector temp_p(daughter.ptSmeared(), daughter.etaSmeared(), daughter.phiSmeared(), o2::constants::physics::MassElectron);
          pPos[0] = temp_p_gen;
          pPos[1] = temp_p;
          wPos[0] = daughter.weight();
          wPos[1] = daughter.weight()*daughter.efficiency();
          nPos++;
          continue;
        }
        other_daughter_pdg = abs(other_daughter_pdg * temp_pdg);
      }
      if (!(((nEle == 1) && (nPos == 1)) || ((nEle == 2) && (nPos == 2)))) {
        LOG(error) << "Found decay with wrong number of electrons: nElectrons = " << nEle << ", nPositrons = " << nPos;
        continue;
      }
      if ((nEle == 2) && (nPos == 2) && (other_daughter_pdg == -1)) {
        other_daughter_pdg = -2;
        // ToDO
        // LOG(warning) << "Found a decay to four electrons of meson with pdg = " << pdg <<". Not included yet";
        continue;
      }
      if (decays.find(other_daughter_pdg) == decays.end()) {
        LOG(error) << "Found decay with code = " << other_daughter_pdg << " that is not in list of decays";
        continue;
      }

      for (int s=0; s<2;s++){
        TString histname;

       /* if ((s==1) && (!isAccepted(pEle[s],pPos[s]))){
          continue; 
        }*/
        // missing acceptance
        // not OK now. need separate for Ele and Pos and only later for pair

        histname="Pt";
        fillHistogram1D(histname, s, pdg, other_daughter_pdg, pEle[s].Pt(),wEle[s]);
        fillHistogram1D(histname, s, pdg, other_daughter_pdg, pPos[s].Pt(),wPos[s]);

        histname="Eta";
        fillHistogram1D(histname, s, pdg, other_daughter_pdg, pEle[s].Eta(),wEle[s]);
        fillHistogram1D(histname, s, pdg, other_daughter_pdg, pPos[s].Eta(),wPos[s]);

        histname="Phi";
        fillHistogram1D(histname, s, pdg, other_daughter_pdg, pEle[s].Phi(),wEle[s]);
        fillHistogram1D(histname, s, pdg, other_daughter_pdg, pPos[s].Phi(),wPos[s]);


        ROOT::Math::PtEtaPhiMVector p12=pEle[s]+pPos[s];
        float mee = p12.M();
        float ptee = p12.Pt();


        histname="MeeVsPtee";
        fillHistogram1D(histname, s, pdg, other_daughter_pdg, mee, ptee,wEle[s]*wPos[s]);
      }

    } // end particle loop
  }
  PROCESS_SWITCH(lmeelfcocktail, processCocktail, "Process cocktail", true);

  // ULS and LS spectra
  Preslice<McParticlesSmeared> perCollision = aod::mcparticle::mcCollisionId;
  Partition<McParticlesSmeared> Electrons = (aod::mcparticle::pdgCode == 11);
  Partition<McParticlesSmeared> Positrons = (aod::mcparticle::pdgCode == -11);

  void processPairing(aod::McCollisions const& mcCollisions, McParticlesSmeared const& mcParticles)
  {

    for (auto const& mcCollision : mcCollisions) {
      auto const electronsGrouped = Electrons->sliceBy(perCollision, mcCollision.globalIndex());
      auto const positronsGrouped = Positrons->sliceBy(perCollision, mcCollision.globalIndex());

      // ULS spectrum
      for (auto const& [p1, p2] : combinations(o2::soa::CombinationsFullIndexPolicy(electronsGrouped, positronsGrouped))) {
        if (!( from_primary(p1, mcParticles) && from_primary(p2, mcParticles))){
          continue;
        }
        ROOT::Math::PtEtaPhiMVector v1_gen(p1.pt(), p1.eta(), p1.phi(), o2::constants::physics::MassElectron);
        ROOT::Math::PtEtaPhiMVector v2_gen(p2.pt(), p2.eta(), p2.phi(), o2::constants::physics::MassElectron);
        ROOT::Math::PtEtaPhiMVector v12_gen = v1_gen + v2_gen;
        registry.fill(HIST("gen/ULS"), v12_gen.M(), v12_gen.Pt(), p1.weight() * p2.weight());
        if (isAcceptedPair(p1, p2)) {
          ROOT::Math::PtEtaPhiMVector v1(p1.ptSmeared(), p1.etaSmeared(), p1.phiSmeared(), o2::constants::physics::MassElectron);
          ROOT::Math::PtEtaPhiMVector v2(p2.ptSmeared(), p2.etaSmeared(), p2.phiSmeared(), o2::constants::physics::MassElectron);
          ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
          registry.fill(HIST("rec/ULS"), v12.M(), v12.Pt(), p1.weight() * p2.weight() * p1.efficiency() * p2.efficiency());
        }
      }
      // LS spectra
      for (auto& [p1, p2] : combinations(o2::soa::CombinationsStrictlyUpperIndexPolicy(electronsGrouped, electronsGrouped))) {
        if (!( from_primary(p1, mcParticles) && from_primary(p2, mcParticles))){
          continue;
        }
        ROOT::Math::PtEtaPhiMVector v1_gen(p1.pt(), p1.eta(), p1.phi(), o2::constants::physics::MassElectron);
        ROOT::Math::PtEtaPhiMVector v2_gen(p2.pt(), p2.eta(), p2.phi(), o2::constants::physics::MassElectron);
        ROOT::Math::PtEtaPhiMVector v12_gen = v1_gen + v2_gen;
        registry.fill(HIST("gen/LSmm"), v12_gen.M(), v12_gen.Pt(), p1.weight() * p2.weight());
        if (isAcceptedPair(p1, p2)) {
          ROOT::Math::PtEtaPhiMVector v1(p1.ptSmeared(), p1.etaSmeared(), p1.phiSmeared(), o2::constants::physics::MassElectron);
          ROOT::Math::PtEtaPhiMVector v2(p2.ptSmeared(), p2.etaSmeared(), p2.phiSmeared(), o2::constants::physics::MassElectron);
          ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
          registry.fill(HIST("rec/LSmm"), v12.M(), v12.Pt(), p1.weight() * p2.weight() * p1.efficiency() * p2.efficiency());
        }
      }
      for (auto& [p1, p2] : combinations(o2::soa::CombinationsStrictlyUpperIndexPolicy(positronsGrouped, positronsGrouped))) {
        if (!( from_primary(p1, mcParticles) && from_primary(p2, mcParticles))){
          continue;
        }
        ROOT::Math::PtEtaPhiMVector v1_gen(p1.pt(), p1.eta(), p1.phi(), o2::constants::physics::MassElectron);
        ROOT::Math::PtEtaPhiMVector v2_gen(p2.pt(), p2.eta(), p2.phi(), o2::constants::physics::MassElectron);
        ROOT::Math::PtEtaPhiMVector v12_gen = v1_gen + v2_gen;
        registry.fill(HIST("gen/LSpp"), v12_gen.M(), v12_gen.Pt(), p1.weight() * p2.weight());
        if (isAcceptedPair(p1, p2)) {
          ROOT::Math::PtEtaPhiMVector v1(p1.ptSmeared(), p1.etaSmeared(), p1.phiSmeared(), o2::constants::physics::MassElectron);
          ROOT::Math::PtEtaPhiMVector v2(p2.ptSmeared(), p2.etaSmeared(), p2.phiSmeared(), o2::constants::physics::MassElectron);
          ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
          registry.fill(HIST("rec/LSpp"), v12.M(), v12.Pt(), p1.weight() * p2.weight() * p1.efficiency() * p2.efficiency());
        }
      }
    }
  }
  PROCESS_SWITCH(lmeelfcocktail, processPairing, "Process ULS and LS pairing", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<lmeelfcocktail>(cfgc, TaskName("em-lmee-lf-cocktail"))};
}
