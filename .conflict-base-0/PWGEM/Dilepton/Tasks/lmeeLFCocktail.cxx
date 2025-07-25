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

#include <map>
#include <vector>
#include <string>

#include "Math/Vector4D.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "PWGEM/Dilepton/DataModel/dileptonTables.h"
#include "PWGEM/Dilepton/Utils/PairUtilities.h"

using namespace o2;
using namespace o2::framework;

using McParticlesSmeared = soa::Join<aod::McParticles, aod::SmearedElectrons>;

struct lmeelfcocktail {

  enum recLevel {
    kGen = 0,
    kRec,
    kNRecLevels
  };

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
    {331, {"etaP/", {22, 223, 211 * 211}}},
    {113, {"rho/", {-1}}},
    {223, {"omega/", {-1, 111}}},
    {333, {"phi/", {-1, 111, 221}}},
    {443, {"Jpsi/", {-1}}},
    {100443, {"psi2S/", {-1}}},
    {553, {"Upsilon/", {-1}}}};

  std::map<TString, int> histogramId;

  std::vector<TString> stage = {"gen/", "rec/"};

  HistogramRegistry registry{"registry", {}};

  Configurable<bool> fConfigApplyDEtaDPhi{"cfgApplyDEtaDPhi", false, "flag to apply deta-phi cut"};
  Configurable<float> fConfigMinDEta{"cfgMinDEta", 0.02, "minimum deta"};
  Configurable<float> fConfigMinDPhi{"cfgMinDPhi", 0.2, "minimum dphi"};

  Configurable<float> fConfigMaxEta{"cfgMaxEta", 0.8, "maxium |eta|"};
  Configurable<float> fConfigMinPt{"cfgMinPt", 0.2, "minium pT"};
  Configurable<float> fConfigMaxPt{"cfgMaxPt", 8.0, "maximum pT"};
  Configurable<float> fConfigMinOpAng{"cfgMinOpAng", 0, "minimum opening angle"};
  Configurable<float> fConfigMinPtee{"cfgMinPtee", 0, "minimum pair pT"};
  Configurable<float> fConfigMaxRapee{"cfgMaxRapee", 999., "maximum pair rapidity"};
  ConfigurableAxis fConfigMeeBins{"cfgMeeBins", {600, 0.f, 6.f}, "Mee binning"};
  ConfigurableAxis fConfigPteeBins{"cfgPteeBins", {400, 0.f, 10.f}, "pTee binning"};
  ConfigurableAxis fConfigCos2DPhi{"cfgCos2DPhi", {200, -1.f, +1.f}, "cos(2x(phiee-PsiRP)) binning"}; // for dilepton v2.
  ConfigurableAxis fConfigPtBins{"cfgPtBins", {200, 0.f, 10.f}, "pT binning"};
  ConfigurableAxis fConfigEtaBins{"cfgEtaBins", {200, -5.f, 5.f}, "eta binning"};
  ConfigurableAxis fConfigPhiBins{"cfgPhiBins", {200, -TMath::Pi(), TMath::Pi()}, "phi binning"};
  ConfigurableAxis fConfigPhiVBins{"cfgPhiVBins", {200, 0, TMath::Pi()}, "phiV binning"};
  ConfigurableAxis fConfigOpAngBins{"cfgOpAngBins", {200, 0, TMath::Pi()}, "opening angle binning"};
  ConfigurableAxis fConfigDcaBins{"cfgDcaBins", {VARIABLE_WIDTH, 0., 0.4, 0.8, 1.2, 1.6, 2.0, 2.4, 3., 4., 5., 7., 10.}, "dca binning"};

  std::vector<std::shared_ptr<TH1>> histograms1D;
  std::vector<std::shared_ptr<TH2>> histograms2D;
  std::vector<std::shared_ptr<THnSparse>> histogramsND;

  template <typename T, typename U>
  bool from_primary(T& p1, U& mcParticles)
  {
    if (!p1.has_mothers()) {
      return false;
    }
    auto mother = mcParticles.iteratorAt(p1.mothersIds()[0]);
    if (mother.has_mothers()) {
      return false;
    }
    return true;
  }

  bool isAcceptedSingle(ROOT::Math::PtEtaPhiMVector p1)
  {
    if (p1.Pt() < fConfigMinPt)
      return false;
    if (p1.Pt() > fConfigMaxPt)
      return false;
    if (fabs(p1.Eta()) > fConfigMaxEta)
      return false;
    return true;
  }

  bool isAcceptedPair(ROOT::Math::PtEtaPhiMVector p1, ROOT::Math::PtEtaPhiMVector p2)
  {
    if (!isAcceptedSingle(p1)) {
      return false;
    }
    if (!isAcceptedSingle(p2)) {
      return false;
    }
    ROOT::Math::PtEtaPhiMVector p12 = p1 + p2;
    if (p12.Pt() < fConfigMinPtee) {
      return false;
    }
    if (o2::aod::pwgem::dilepton::utils::pairutil::getOpeningAngle(p1.Px(), p1.Py(), p1.Pz(), p2.Px(), p2.Py(), p2.Pz()) < fConfigMinOpAng) {
      return false;
    }
    if (fabs(p12.Rapidity()) > fConfigMaxRapee) {
      return false;
    }

    float deta = p1.Eta() - p2.Eta();
    float dphi = p1.Phi() - p2.Phi();
    o2::math_utils::bringToPMPi(dphi);
    if (fConfigApplyDEtaDPhi && std::pow(deta / fConfigMinDEta, 2) + std::pow(dphi / fConfigMinDPhi, 2) < 1.f) {
      return false;
    }
    return true;
  }

  template <typename T>
  bool isAcceptedPair(T& p1, T& p2)
  {
    ROOT::Math::PtEtaPhiMVector v1(p1.ptSmeared(), p1.etaSmeared(), p1.phiSmeared(), o2::constants::physics::MassElectron);
    ROOT::Math::PtEtaPhiMVector v2(p2.ptSmeared(), p2.etaSmeared(), p2.phiSmeared(), o2::constants::physics::MassElectron);
    return isAcceptedPair(v1, v2);
  }

  template <typename T>
  bool isAcceptedSingle(T& p1)
  {
    ROOT::Math::PtEtaPhiMVector v1(p1.ptSmeared(), p1.etaSmeared(), p1.phiSmeared(), o2::constants::physics::MassElectron);
    return isAcceptedSingle(v1);
  }

  void addHistogram1D_stage(TString histname, AxisSpec axis, int& i, TString s)
  {
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

  void addHistogram1D(TString histname, AxisSpec axis, int& i)
  {
    for (auto s : stage) {
      addHistogram1D_stage(histname, axis, i, s);
    }
  }

  void addHistogram1D_mother(TString histname, AxisSpec axis, int& i) // mother histograms only for gen. level, no decay channels
  {
    i++;
    TString name = stage[0] + histname;
    histogramId[name] = i;
    histograms1D.push_back(registry.add<TH1>(name, histname, HistType::kTH1F, {axis}, true));
    for (auto const& [pdg, meson] : mesons) {
      i++;
      name = stage[0] + meson.name + histname;
      histogramId[name] = i;
      histograms1D.push_back(registry.add<TH1>(name, histname, HistType::kTH1F, {axis}, true));
    }
  }

  void addHistogram2D_stage(TString histname, AxisSpec axis1, AxisSpec axis2, int& i, TString s)
  {
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

  void addHistogram2D(TString histname, AxisSpec axis1, AxisSpec axis2, int& i)
  {
    for (auto s : stage) {
      addHistogram2D_stage(histname, axis1, axis2, i, s);
    }
  }

  template <typename TAxes>
  void addHistogramND_stage(TString histname, TAxes const& axes, int& i, TString s)
  {
    i++;
    TString name = s + histname;
    histogramId[name] = i;
    histogramsND.push_back(registry.add<THnSparse>(name, histname, HistType::kTHnSparseF, axes, true));
    for (auto const& [pdg, meson] : mesons) {
      i++;
      name = s + meson.name + histname;
      histogramId[name] = i;
      histogramsND.push_back(registry.add<THnSparse>(name, histname, HistType::kTHnSparseF, axes, true));
      for (auto const& mode : meson.decayModes) {
        i++;
        name = s + meson.name + decays[mode] + histname;
        histogramId[name] = i;
        histogramsND.push_back(registry.add<THnSparse>(name, histname, HistType::kTHnSparseF, axes, true));
      }
    }
  }

  template <typename TAxes>
  void addHistogramND(TString histname, TAxes const& axes, int& i)
  {
    for (auto s : stage) {
      addHistogramND_stage(histname, axes, i, s);
    }
  }

  void fillHistogram1D(TString histname, int s, int pdg, int other_daughter_pdg, float value, float weight)
  {
    histograms1D[histogramId[stage[s] + histname]]->Fill(value, weight);
    histograms1D[histogramId[stage[s] + mesons[pdg].name + histname]]->Fill(value, weight);
    histograms1D[histogramId[stage[s] + mesons[pdg].name + decays[other_daughter_pdg] + histname]]->Fill(value, weight);
  }

  void fillHistogram1D_mother(TString histname, int pdg, float value, float weight)
  {
    histograms1D[histogramId[stage[0] + histname]]->Fill(value, weight);
    histograms1D[histogramId[stage[0] + mesons[pdg].name + histname]]->Fill(value, weight);
  }

  void fillHistogram2D(TString histname, int s, int pdg, int other_daughter_pdg, float value1, float value2, float weight)
  {
    histograms2D[histogramId[stage[s] + histname]]->Fill(value1, value2, weight);
    histograms2D[histogramId[stage[s] + mesons[pdg].name + histname]]->Fill(value1, value2, weight);
    histograms2D[histogramId[stage[s] + mesons[pdg].name + decays[other_daughter_pdg] + histname]]->Fill(value1, value2, weight);
  }

  void fillHistogramND(TString histname, int s, int pdg, int other_daughter_pdg, double* values, double weight)
  {
    histogramsND[histogramId[stage[s] + histname]]->Fill(values, weight);
    histogramsND[histogramId[stage[s] + mesons[pdg].name + histname]]->Fill(values, weight);
    histogramsND[histogramId[stage[s] + mesons[pdg].name + decays[other_daughter_pdg] + histname]]->Fill(values, weight);
  }

  void init(InitContext& context)
  {
    registry.add<TH1>("NEvents", "NEvents", HistType::kTH1F, {{1, 0., 1.}}, false);

    AxisSpec mass_axis = {fConfigMeeBins, "m_{ee} (GeV/#it{c}^{2})"};
    AxisSpec ptee_axis = {fConfigPteeBins, "#it{p}_{T,ee} (GeV/#it{c})"};
    AxisSpec cos2dphi_axis = {fConfigCos2DPhi, "cos(2(#varphi_{ee} - #Psi_{RP}))"}; // PsiRP = 0 rad. in generator.
    AxisSpec eta_axis = {fConfigEtaBins, "#it{#eta}_{e}"};
    AxisSpec pt_axis = {fConfigPtBins, "#it{p}_{T,e} (GeV/c)"};
    AxisSpec phi_axis = {fConfigPhiBins, "#it{#varphi}_{e}"};
    AxisSpec phiV_axis = {fConfigPhiVBins, "#it{#varphi}_{V,ee}"};
    AxisSpec opAng_axis = {fConfigOpAngBins, "#it{#omega}_{ee}"};
    AxisSpec eta_axis_mother = {fConfigEtaBins, "#it{#eta}_{mother}"};
    AxisSpec pt_axis_mother = {fConfigPtBins, "#it{p}_{T,mother} (GeV/#it{c})"};
    AxisSpec phi_axis_mother = {fConfigPhiBins, "#it{#varphi}_{mother}"};
    AxisSpec dca_axis = {fConfigDcaBins, "DCA_{e}"};
    AxisSpec dcaee_axis = {fConfigDcaBins, "DCA_{ee}"};

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
      addHistogram1D_mother("Mother_Pt", pt_axis_mother, i);
      addHistogram1D_mother("Mother_Eta", eta_axis_mother, i);
      addHistogram1D_mother("Mother_Phi", phi_axis_mother, i);
      addHistogram1D("PhiV", phiV_axis, i);
      addHistogram1D("OpAng", opAng_axis, i);
      addHistogram1D("Mee", mass_axis, i);
      addHistogram1D("Ptee", ptee_axis, i);
      // addHistogram1D_stage("Dca", dca_axis, i, "rec/");
      // addHistogram1D_stage("Dcaee", dcaee_axis, i, "rec/");
      i = -1;
      // addHistogram2D_stage("DcaVsPt", dca_axis, pt_axis, i, "rec/");
      // addHistogram2D_stage("DcaeeVsPtee", dcaee_axis, ptee_axis, i, "rec/");
      // addHistogram2D_stage("DcaeeVsMee", dcaee_axis, mass_axis, i, "rec/");
      i = -1;
      addHistogramND("MeeVsPteeVsCos2DPhiRP", std::vector<AxisSpec>{mass_axis, ptee_axis, cos2dphi_axis}, i);
    }
  }

  void processCocktail(aod::McCollision const&, McParticlesSmeared const& mcParticles)
  {
    registry.fill(HIST("NEvents"), 0.5);

    for (auto const& particle : mcParticles) {
      if (particle.has_mothers()) {
        continue;
      }
      int pdg = abs(particle.pdgCode());
      if (mesons.find(pdg) == mesons.end()) {
        LOG(error) << "Found mother particle with pdg = " << pdg << " that is not in list of mesons";
      }
      if (!particle.has_daughters()) {
        LOG(error) << "Found meson with pdg = " << pdg << "that has no daughters";
      }

      int other_daughter_pdg = -1;
      int nEle = 0;
      int nPos = 0;

      ROOT::Math::PtEtaPhiMVector pEleGen, pPosGen, pEleRec, pPosRec;
      float weight(1.), effEle(1.), effPos(1.); //, dcaEle(0.), dcaPos(0.);
      for (const auto& daughter : particle.daughters_as<McParticlesSmeared>()) {
        int temp_pdg = daughter.pdgCode();
        if (temp_pdg == 11) {
          ROOT::Math::PtEtaPhiMVector temp_p_gen(daughter.pt(), daughter.eta(), daughter.phi(), o2::constants::physics::MassElectron);
          ROOT::Math::PtEtaPhiMVector temp_p(daughter.ptSmeared(), daughter.etaSmeared(), daughter.phiSmeared(), o2::constants::physics::MassElectron);
          pEleGen = temp_p_gen;
          pEleRec = temp_p;
          weight = daughter.weight();
          effEle = daughter.efficiency();
          // dcaEle = daughter.dca();
          nEle++;
          continue;
        }
        if (temp_pdg == -11) {
          ROOT::Math::PtEtaPhiMVector temp_p_gen(daughter.pt(), daughter.eta(), daughter.phi(), o2::constants::physics::MassElectron);
          ROOT::Math::PtEtaPhiMVector temp_p(daughter.ptSmeared(), daughter.etaSmeared(), daughter.phiSmeared(), o2::constants::physics::MassElectron);
          pPosGen = temp_p_gen;
          pPosRec = temp_p;
          effPos = daughter.efficiency();
          // dcaPos = daughter.dca();
          nPos++;
          continue;
        }
        other_daughter_pdg = abs(other_daughter_pdg * temp_pdg);
      }
      if (!(((nEle == 1) && (nPos == 1)) || ((nEle == 2) && (nPos == 2)))) {
        LOG(error) << "Found decay with wrong number of electrons in decay of meson with pdg " << pdg << ": nElectrons = " << nEle << ", nPositrons = " << nPos;
        continue;
      }
      if ((nEle == 2) && (nPos == 2) && (other_daughter_pdg == -1)) {
        other_daughter_pdg = -2;
        weight = 2 * weight;
      }
      auto this_meson_decays = mesons[pdg].decayModes;
      if (std::find(this_meson_decays.begin(), this_meson_decays.end(), other_daughter_pdg) == this_meson_decays.end()) {
        LOG(error) << "Found decay with code = " << other_daughter_pdg << " that is not in list of decays of meson with pdg " << pdg;
        continue;
      }

      for (int s = 0; s < kNRecLevels; s++) { // s=0: gen, s=1: rec

        ROOT::Math::PtEtaPhiMVector pEle, pPos;
        float pairWeight(1.), weightEle(1.), weightPos(1.);
        bool acceptedEle(true), acceptedPos(true), acceptedPair(true);

        if (s == kGen) {
          pEle = pEleGen;
          pPos = pPosGen;
          pairWeight = weight;
          weightEle = weight;
          weightPos = weight;
          acceptedEle = true;
          acceptedPos = true;
          acceptedPair = true;
        } else if (s == kRec) {
          pEle = pEleRec;
          pPos = pPosRec;
          pairWeight = weight * effEle * effPos;
          weightEle = weight * effEle;
          weightPos = weight * effPos;
          acceptedEle = isAcceptedSingle(pEle);
          acceptedPos = isAcceptedSingle(pPos);
          acceptedPair = isAcceptedPair(pEle, pPos);
        }

        // single track histograms
        if (acceptedEle)
          fillHistogram1D("Pt", s, pdg, other_daughter_pdg, pEle.Pt(), weightEle);
        if (acceptedPos)
          fillHistogram1D("Pt", s, pdg, other_daughter_pdg, pPos.Pt(), weightPos);

        if (acceptedEle)
          fillHistogram1D("Eta", s, pdg, other_daughter_pdg, pEle.Eta(), weightEle);
        if (acceptedPos)
          fillHistogram1D("Eta", s, pdg, other_daughter_pdg, pPos.Eta(), weightPos);

        if (acceptedEle)
          fillHistogram1D("Phi", s, pdg, other_daughter_pdg, pEle.Phi(), weightEle);
        if (acceptedPos)
          fillHistogram1D("Phi", s, pdg, other_daughter_pdg, pPos.Phi(), weightPos);

        if (s == kRec) { // dca only at rec. level
          if (acceptedEle) {
            // fillHistogram1D("Dca", s, pdg, other_daughter_pdg, dcaEle, weightEle);
            // fillHistogram2D("DcaVsPt", s, pdg, other_daughter_pdg, dcaEle, pEle.Pt(), weightEle);
          }
          if (acceptedPos) {
            // fillHistogram1D("Dca", s, pdg, other_daughter_pdg, dcaPos, weightPos);
            // fillHistogram2D("DcaVsPt", s, pdg, other_daughter_pdg, dcaPos, pPos.Pt(), weightPos);
          }
        }

        // mother histograms
        if (s == kGen) { // only at gen. level
          fillHistogram1D_mother("Mother_Pt", pdg, particle.pt(), weight);
          fillHistogram1D_mother("Mother_Eta", pdg, particle.eta(), weight);
          fillHistogram1D_mother("Mother_Phi", pdg, particle.phi(), weight);
        }

        // pair historams
        if (acceptedPair) {
          ROOT::Math::PtEtaPhiMVector p12 = pEle + pPos;
          float mee = p12.M();
          float ptee = p12.Pt();
          float opAng = o2::aod::pwgem::dilepton::utils::pairutil::getOpeningAngle(pPos.Px(), pPos.Py(), pPos.Pz(), pEle.Px(), pEle.Py(), pEle.Pz());
          float phiV = o2::aod::pwgem::dilepton::utils::pairutil::getPhivPair(pPos.Px(), pPos.Py(), pPos.Pz(), pEle.Px(), pEle.Py(), pEle.Pz(), 1, -1, 1);
          // float dcaee = sqrt((pow(dcaEle, 2) + pow(dcaPos, 2)) / 2);
          float cos2dphi = std::cos(2.f * p12.Phi()); // PsiRP = 0 rad.
          double values[3] = {mee, ptee, cos2dphi};
          fillHistogramND("MeeVsPteeVsCos2DPhiRP", s, pdg, other_daughter_pdg, values, pairWeight);
          fillHistogram1D("Mee", s, pdg, other_daughter_pdg, mee, pairWeight);
          fillHistogram1D("Ptee", s, pdg, other_daughter_pdg, ptee, pairWeight);
          fillHistogram1D("PhiV", s, pdg, other_daughter_pdg, phiV, pairWeight);
          fillHistogram1D("OpAng", s, pdg, other_daughter_pdg, opAng, pairWeight);
          if (s == kRec) { // dca only at rec. level
            // fillHistogram1D("Dcaee", s, pdg, other_daughter_pdg, dcaee, pairWeight);
            // fillHistogram2D("DcaeeVsPtee", s, pdg, other_daughter_pdg, dcaee, ptee, pairWeight);
            // fillHistogram2D("DcaeeVsMee", s, pdg, other_daughter_pdg, dcaee, mee, pairWeight);
          }
        }
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
        if (!(from_primary(p1, mcParticles) && from_primary(p2, mcParticles))) {
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
        if (!(from_primary(p1, mcParticles) && from_primary(p2, mcParticles))) {
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
        if (!(from_primary(p1, mcParticles) && from_primary(p2, mcParticles))) {
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
