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
// Analysis task for lmee light flavour cocktail

#include <vector>
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/Logger.h"
#include "Math/Vector4D.h"
#include "Math/Vector3D.h"
#include "TFile.h"
#include "TF1.h"
#include "TRandom.h"
#include "TDatabasePDG.h"
#include "TGenPhaseSpace.h"

using namespace o2;
using namespace o2::framework;
using namespace ROOT::Math;

namespace o2::aod
{
namespace eeTable
{
DECLARE_SOA_COLUMN(fd1DCA, d1DCA, float);
DECLARE_SOA_COLUMN(fd2DCA, d2DCA, float);
DECLARE_SOA_COLUMN(fpairDCA, pairDCA, float);
DECLARE_SOA_COLUMN(fd1origpt, d1origpt, float);
DECLARE_SOA_COLUMN(fd1origp, d1origp, float);
DECLARE_SOA_COLUMN(fd1origeta, d1origega, float);
DECLARE_SOA_COLUMN(fd1origphi, d1origphi, float);
DECLARE_SOA_COLUMN(fd2origpt, d2origpt, float);
DECLARE_SOA_COLUMN(fd2origp, d2origp, float);
DECLARE_SOA_COLUMN(fd2origeta, d2origeta, float);
DECLARE_SOA_COLUMN(fd2origphi, d2origphi, float);
DECLARE_SOA_COLUMN(fd1pt, d1pt, float);
DECLARE_SOA_COLUMN(fd1p, d1p, float);
DECLARE_SOA_COLUMN(fd1eta, d1eta, float);
DECLARE_SOA_COLUMN(fd1phi, d1phi, float);
DECLARE_SOA_COLUMN(fd2pt, d2pt, float);
DECLARE_SOA_COLUMN(fd2p, d2p, float);
DECLARE_SOA_COLUMN(fd2eta, d2eta, float);
DECLARE_SOA_COLUMN(fd2phi, d2phi, float);
DECLARE_SOA_COLUMN(feeorigpt, eeorigpt, float);
DECLARE_SOA_COLUMN(feeorigp, eeorigp, float);
DECLARE_SOA_COLUMN(feeorigm, eeorigm, float);
DECLARE_SOA_COLUMN(feeorigeta, eeorigeta, float);
DECLARE_SOA_COLUMN(feeorigphi, eeorigphi, float);
DECLARE_SOA_COLUMN(feeorigphiv, eeorigphiv, float);
DECLARE_SOA_COLUMN(feept, eept, float);
DECLARE_SOA_COLUMN(feemt, eemt, float);
DECLARE_SOA_COLUMN(feep, eep, float);
DECLARE_SOA_COLUMN(feem, eem, float);
DECLARE_SOA_COLUMN(feeeta, eeeta, float);
DECLARE_SOA_COLUMN(feephi, eephi, float);
DECLARE_SOA_COLUMN(feephiv, eephiv, float);
DECLARE_SOA_COLUMN(fmotherpt, motherpt, float);
DECLARE_SOA_COLUMN(fmothermt, mothermt, float);
DECLARE_SOA_COLUMN(fmotherp, motherp, float);
DECLARE_SOA_COLUMN(fmotherm, motherm, float);
DECLARE_SOA_COLUMN(fmothereta, mothereta, float);
DECLARE_SOA_COLUMN(fmotherphi, motherphi, float);
DECLARE_SOA_COLUMN(fID, ID, int);
DECLARE_SOA_COLUMN(fdectyp, dectyp, uint);
DECLARE_SOA_COLUMN(fdau3pdg, dau3pdg, int);
DECLARE_SOA_COLUMN(fweight, weight, double);
DECLARE_SOA_COLUMN(fwEffpT, wEffpT, double);
DECLARE_SOA_COLUMN(fwMultpT, wMultpT, double);
DECLARE_SOA_COLUMN(fwMultmT, wMultmT, double);
DECLARE_SOA_COLUMN(fwMultpT2, wMultpT2, double);
DECLARE_SOA_COLUMN(fwMultmT2, wMultmT2, double);
DECLARE_SOA_COLUMN(fpass, pass, bool);

} // namespace eeTable
DECLARE_SOA_TABLE(eeTTree, "AOD", "EETTREE",
                  eeTable::fd1DCA, eeTable::fd2DCA, eeTable::fpairDCA, eeTable::fd1origpt, eeTable::fd1origp, eeTable::fd1origeta, eeTable::fd1origphi, eeTable::fd2origpt, eeTable::fd2origp, eeTable::fd2origeta, eeTable::fd2origphi, eeTable::fd1pt, eeTable::fd1p, eeTable::fd1eta, eeTable::fd1phi, eeTable::fd2pt, eeTable::fd2p, eeTable::fd2eta, eeTable::fd2phi, eeTable::feeorigpt, eeTable::feeorigp, eeTable::feeorigm, eeTable::feeorigeta, eeTable::feeorigphi, eeTable::feeorigphiv, eeTable::feept, eeTable::feemt, eeTable::feep, eeTable::feem, eeTable::feeeta, eeTable::feephi, eeTable::feephiv, eeTable::fmotherpt, eeTable::fmothermt, eeTable::fmotherp, eeTable::fmotherm, eeTable::fmothereta, eeTable::fmotherphi, eeTable::fID, eeTable::fdectyp, eeTable::fdau3pdg, eeTable::fweight, eeTable::fwEffpT, eeTable::fwMultpT, eeTable::fwMultmT, eeTable::fwMultpT2, eeTable::fwMultmT2, eeTable::fpass);

} // namespace o2::aod

struct lmeelfcocktail {

  Produces<aod::eeTTree> tree;
  HistogramRegistry registry{"registry", {}};
  Int_t nInputParticles = 17;
  std::vector<TString> fParticleListNames = {"Pi0", "Eta", "EtaP", "EtaP_dalitz_photon", "EtaP_dalitz_omega", "Rho", "Omega", "Omega_2body", "Omega_dalitz", "Phi", "Phi_2body", "Phi_dalitz_eta", "Phi_dalitz_pi0", "Jpsi", "Jpsi_2body", "Jpsi_radiative", "Virtual_Photon"};
  // std::vector<Int_t> fParticleList = {111, 221, 331, 223331, 2233331, 113, 223, 2223, 3223, 333, 2333, 2213333, 1113333, 443,2443,223443,000};
  TH1F* fhwEffpT;
  TH1F* fhwMultpT;
  TH1F* fhwMultmT;
  TH1F* fhwMultpT2;
  TH1F* fhwMultmT2;
  TH1F* fhKW;
  TF1* ffVPHpT;
  TObjArray* fArr;
  TObjArray* fArrResoPt;
  TObjArray* fArrResoEta;
  TObjArray* fArrResoPhi_Pos;
  TObjArray* fArrResoPhi_Neg;

  std::vector<std::shared_ptr<TH1>> fmee_orig, fmotherpT_orig, fphi_orig, frap_orig, fmee_orig_wALT, fmotherpT_orig_wALT, fmee, fphi, frap, fmee_wALT;
  std::vector<std::shared_ptr<TH2>> fpteevsmee_wALT, fpteevsmee_orig_wALT, fpteevsmee_orig, fpteevsmee;

  std::vector<double> DCATemplateEdges;
  int nbDCAtemplate;
  TH1F** fh_DCAtemplates;

  Configurable<int> fCollisionSystem{"cfgCollisionSystem", 200, "set the collision system"};
  Configurable<bool> fConfigWriteTTree{"cfgWriteTTree", false, "whether tree output should be written"};
  Configurable<float> fConfigMaxEta{"cfgMaxEta", 0.8, "maxium |eta|"};
  Configurable<float> fConfigMinPt{"cfgMinPt", 0.2, "minium pT"};
  Configurable<float> fConfigMaxPt{"cfgMaxPt", 8.0, "maximum pT"};
  Configurable<int> fConfigNBinsMee{"cfgNBinsMee", 1200, "number of bins in invariant mass"};
  Configurable<float> fConfigMinMee{"cfgMinMee", 0.0, "lowest bin in invariant mass"};
  Configurable<float> fConfigMaxMee{"cfgMaxMee", 6.0, "highest bin in invariant mass"};
  Configurable<int> fConfigNBinsPtee{"cfgNBinsPtee", 400, "number of bins in pT"};
  Configurable<float> fConfigMinPtee{"cfgMinPtee", 0.0, "lowest bin in pT"};
  Configurable<float> fConfigMaxPtee{"cfgMaxPtee", 10.0, "hightest bin in pT"};
  Configurable<int> fConfigResolType{"cfgResolType", 2, "set resolution type"};
  Configurable<int> fConfigALTweight{"cfgALTweight", 1, "set alternative weighting type"};
  Configurable<std::string> fConfigResFileName{"cfgResFileName", "", "name of resolution file"};
  // Configurable<bool> fConfigResFileLocal{"cfgResFileLocal", false, "..."};
  Configurable<std::string> fConfigEffFileName{"cfgEffFileName", "", "name of efficiency file"};
  Configurable<float> fConfigMinOpAng{"cfgMinOpAng", 0.050, "minimum opening angle"};
  Configurable<bool> fConfigDoRapidityCut{"cfgDoRapidityCut", false, "apply rapidity cut"};
  Configurable<float> fConfigRapidityCut{"cfgRapidityCut", 1.0, "rapdity cut"};
  Configurable<int> fConfigNBinsPhi{"cfgNBinsPhi", 240, "number of bins in phi"};
  Configurable<float> fConfigMinPhi{"cfgMinPhi", 0.0, "lowerst bin in phi"};
  Configurable<float> fConfigMaxPhi{"cfgMaxPhi", TMath::TwoPi(), "hightest bin in phi"};
  Configurable<int> fConfigNBinsRap{"cfgNBinsRap", 240, "number of bins in rap"};
  Configurable<float> fConfigMinRap{"cfgMinRap", -1.2, "lowest bin in rap"};
  Configurable<float> fConfigMaxRap{"cfgMaxRap", 1.2, "hightest bin in rap"};
  Configurable<std::string> fConfigEffHistName{"cfgEffHistName", "fhwEffpT", "hisogram name in efficiency file"};
  Configurable<std::string> fConfigResPHistName{"cfgResPHistName", "ptSlices", "histogram name for p in resolution file"};
  Configurable<std::string> fConfigResPtHistName{"cfgResPtHistName", "RelPtResArrCocktail", "histogram name for pt in resolution file"};
  Configurable<std::string> fConfigResEtaHistName{"cfgResEtaHistName", "EtaResArr", "histogram name for eta in resolution file"};
  Configurable<std::string> fConfigResPhiPosHistName{"cfgResPhiPosHistName", "PhiPosResArr", "histogram name for phi pos in resolution file"};
  Configurable<std::string> fConfigResPhiNegHistName{"cfgResPhiNegHistName", "PhiEleResArr", "hisogram for phi neg in resolution file"};
  Configurable<std::string> fConfigDCAFileName{"cfgDCAFileName", "", "DCA file name"};
  Configurable<std::string> fConfigDCAHistName{"cfgDCAHistName", "fh_DCAtemplate", "histogram name in DCA file"};
  Configurable<std::string> fConfigMultFileName{"cfgMultFileName", "", "multiplicity file name"};
  Configurable<std::string> fConfigMultHistPtName{"cfgMultHistPtName", "fhwMultpT", "hisogram name for pt in multiplicity file"};
  Configurable<std::string> fConfigMultHistPt2Name{"cfgMultHistPt2Name", "fhwMultpT_upperlimit", "histogram name for pt 2 in multiplicity file"};
  Configurable<std::string> fConfigMultHistMtName{"cfgMultHistMtName", "fhwMultmT", "histogram name for mt in multiplicity file"};
  Configurable<std::string> fConfigMultHistMt2Name{"cfgMultHistMt2Name", "fhwMultmT_upperlimit", "histogram name for mt 2 in multiplicity file"};
  Configurable<int> fConfigKWNBins{"cfgKWNBins", 10000, "number of bins for Kroll-Wada"};
  Configurable<float> fConfigKWMax{"cfgKWMax", 1.1, "upper bound of Kroll-Wada"};
  Configurable<std::string> fConfigPhotonPtFileName{"cfgPhotonPtFileName", "", "file name for photon pT parametrization"};
  Configurable<std::string> fConfigPhotonPtDirName{"cfgPhotonPtDirName", "7TeV_Comb", "directory name for photon pT parametrization"};
  Configurable<std::string> fConfigPhotonPtFuncName{"cfgPhotonPtFuncName", "111_pt", "function name for photon pT parametrization"};

  ConfigurableAxis fConfigPtBins{"cfgPtBins", {0., 0.5, 1, 1.5, 2., 2.5, 3., 3.5, 4., 4.5, 5., 5.5, 6., 6.5, 7., 7.5, 8.}, "pT bins"};
  ConfigurableAxis fConfigMBins{"cfgMBins", {0., 0.08, 0.14, 0.2, 1.1, 2.7, 2.8, 3.2, 5.0}, "mee bins"};
  ConfigurableAxis fConfigDCABins{"cfgDCABins", {0., 0.4, 0.8, 1.2, 1.6, 2.0, 2.4, 3., 4., 5., 7., 10.}, "DCA bins"};

  Configurable<std::vector<double>> fConfigDCATemplateEdges{"cfgDCATemplateEdges", {0., .3, .4, .6, 1., 2.}, "DCA template edges"};

  Preslice<aod::McParticles> perCollision = aod::mcparticle::mcCollisionId;

  void init(o2::framework::InitContext&)
  {

    DCATemplateEdges = fConfigDCATemplateEdges;
    nbDCAtemplate = DCATemplateEdges.size();
    DCATemplateEdges.push_back(10000000.);

    SetHistograms();
    GetEffHisto(TString(fConfigEffFileName), TString(fConfigEffHistName));
    if (fConfigResolType == 1) {
      GetResHisto(TString(fConfigResFileName), TString(fConfigResPHistName));
    } else if (fConfigResolType == 2) {
      GetResHisto(TString(fConfigResFileName), TString(fConfigResPtHistName), TString(fConfigResEtaHistName), TString(fConfigResPhiPosHistName), TString(fConfigResPhiNegHistName));
    }
    GetDCATemplates(TString(fConfigDCAFileName), TString(fConfigDCAHistName));
    GetMultHisto(TString(fConfigMultFileName), TString(fConfigMultHistPtName), TString(fConfigMultHistPt2Name), TString(fConfigMultHistMtName), TString(fConfigMultHistMt2Name));
    GetPhotonPtParametrization(TString(fConfigPhotonPtFileName), TString(fConfigPhotonPtDirName), TString(fConfigPhotonPtFuncName));
    fillKrollWada();
  }

  void processCocktail(aod::McCollision const&, aod::McParticles const& mcParticles)
  {

    double fwEffpT, fd1origpt, fd1origp, fd1origeta, fd1origphi, fd2origpt, fd2origp, fd2origeta, fd2origphi, feeorigpt, feeorigp, feeorigm, feeorigeta, feeorigphi, feeorigphiv, fpairDCA, fd1DCA = 0., fd2DCA = 0., fd1pt, fd1p, fd1eta, fd1phi, fd2pt, fd2p, fd2eta, fd2phi, feept, feemt, feep, feem, feeeta, feephi, feephiv, fmotherpt, fmothermt, fmotherp, fmotherm, fmothereta, fmotherphi, fID, fweight, fwMultpT, fwMultpT2, fwMultmT, fwMultmT2;
    bool fpass;

    Partition<aod::McParticles> Mothers = ((aod::mcparticle::pdgCode == 111) || (aod::mcparticle::pdgCode == 221) || (aod::mcparticle::pdgCode == 331) || (aod::mcparticle::pdgCode == 113) || (aod::mcparticle::pdgCode == 223) || (aod::mcparticle::pdgCode == 333) || (aod::mcparticle::pdgCode == 443));
    Mothers.bindTable(mcParticles);
    for (auto& mother : Mothers) {

      if (!mother.has_daughters())
        continue;
      if (mother.has_mothers())
        continue;

      int fdectyp = 0;
      int fdau3pdg = 0;
      bool has_e = false;
      bool has_p = false;
      PxPyPzEVector dau1, dau2, ee, ee_orig;

      for (auto& d : mother.daughters_as<aod::McParticles>()) {
        fdectyp++;
        if (d.pdgCode() == 11) {
          has_e = true;
          dau1.SetPxPyPzE(d.px(), d.py(), d.pz(), d.e());
          fweight = d.weight(); // get particle weight from generator
        } else if (d.pdgCode() == -11) {
          has_p = true;
          dau2.SetPxPyPzE(d.px(), d.py(), d.pz(), d.e());
        } else {
          fdau3pdg = d.pdgCode();
        }
      }
      if ((!has_e) || (!has_p))
        continue;

      if (fdectyp > 4)
        continue; // here dectype==4 is included, but when filling histograms it is excluded?

      // Not sure about this cut. From GammaConv group. Harmless a priori.
      if (!(fabs(dau1.E() - dau1.Pz()) > 0.))
        continue; // taken from Run2 version. Does this make sense?

      // taken from Run2 version. why rap cut is applied only to dau1 and not dau2?
      Double_t yPre = (dau1.E() + dau1.Pz()) / (dau1.E() - dau1.Pz());
      Double_t y = 0.5 * TMath::Log(yPre);
      if (fConfigDoRapidityCut) { // Apply rapidity cut on mother consistent with GammaConv group.
        if (yPre <= 0.)
          continue;
        if (TMath::Abs(y) > fConfigRapidityCut)
          continue;
      } else {
        if (yPre == 0.)
          continue;
      }

      // create dielectron before resolution effects:
      ee = dau1 + dau2;
      ee_orig = ee;

      // get index for histograms
      Int_t hindex[3];
      for (Int_t jj = 0; jj < 3; jj++) {
        hindex[jj] = -1;
      }
      switch (mother.pdgCode()) {
        case 111:
          hindex[0] = 0;
          break;
        case 221:
          hindex[0] = 1;
          break;
        case 331:
          hindex[0] = 2;
          if (fdectyp == 3 && fdau3pdg == 22)
            hindex[1] = 3;
          if (fdectyp == 3 && fdau3pdg == 223)
            hindex[1] = 4;
          break;
        case 113:
          hindex[0] = 5;
          break;
        case 223:
          hindex[0] = 6;
          if (fdectyp == 2)
            hindex[1] = 7;
          if (fdectyp == 3 && fdau3pdg == 111)
            hindex[1] = 8;
          break;
        case 333:
          hindex[0] = 9;
          if (fdectyp == 2)
            hindex[1] = 10;
          if (fdectyp == 3 && fdau3pdg == 221)
            hindex[1] = 11;
          if (fdectyp == 3 && fdau3pdg == 111)
            hindex[1] = 12;
          break;
        case 443:
          hindex[0] = 13;
          if (fdectyp == 2)
            hindex[1] = 14;
          if (fdectyp == 3 && fdau3pdg == 22)
            hindex[1] = 15;
          break;
      }

      hindex[2] = nInputParticles;

      if (hindex[0] < 0) {
        LOGP(info, "Error LMeeCocktail hindex[0]<0");
        continue;
      }

      // Fill tree words before resolution/acceptance
      fd1origpt = dau1.Pt();
      fd1origp = dau1.P();
      fd1origeta = dau1.Eta();
      fd1origphi = dau1.Phi();
      fd2origpt = dau2.Pt();
      fd2origp = dau2.P();
      fd2origeta = dau2.Eta();
      fd2origphi = dau2.Phi();
      feeorigpt = ee.Pt();
      feeorigp = ee.P();
      feeorigm = ee.M();
      feeorigeta = ee.Eta();
      feeorigphi = ee.Phi();
      feeorigphiv = PhiV(dau1, dau2);

      // get the efficiency weight
      Int_t effbin = fhwEffpT->FindBin(dau1.Pt());
      fwEffpT = fhwEffpT->GetBinContent(effbin);
      effbin = fhwEffpT->FindBin(dau2.Pt());
      fwEffpT = fwEffpT * fhwEffpT->GetBinContent(effbin);

      // Resolution and acceptance
      dau1 = ApplyResolution(dau1, -1, fConfigResolType);
      dau2 = ApplyResolution(dau2, 1, fConfigResolType);
      fpass = true;
      if (dau1.Pt() < fConfigMinPt || dau2.Pt() < fConfigMinPt)
        fpass = false; // leg pT cut
      if (dau1.Pt() > fConfigMaxPt || dau2.Pt() > fConfigMaxPt)
        fpass = false; // leg pT cut
      if (dau1.Vect().Unit().Dot(dau2.Vect().Unit()) > TMath::Cos(fConfigMinOpAng))
        fpass = false; // opening angle cut
      if (TMath::Abs(dau1.Eta()) > fConfigMaxEta || TMath::Abs(dau2.Eta()) > fConfigMaxEta)
        fpass = false;

      // get the pair DCA (based in smeared pT)
      for (int jj = 0; jj < nbDCAtemplate; jj++) { // loop over DCA templates
        if (dau1.Pt() >= DCATemplateEdges[jj] && dau1.Pt() < DCATemplateEdges[jj + 1]) {
          fd1DCA = fh_DCAtemplates[jj]->GetRandom();
        }
        if (dau2.Pt() >= DCATemplateEdges[jj] && dau2.Pt() < DCATemplateEdges[jj + 1]) {
          fd2DCA = fh_DCAtemplates[jj]->GetRandom();
        }
      }
      fpairDCA = sqrt((pow(fd1DCA, 2) + pow(fd2DCA, 2)) / 2);

      // Fill tree words after resolution/acceptance
      ee = dau1 + dau2;
      fd1pt = dau1.Pt();
      fd1p = dau1.P();
      fd1eta = dau1.Eta();
      fd1phi = dau1.Phi();
      fd2pt = dau2.Pt();
      fd2p = dau2.P();
      fd2eta = dau2.Eta();
      fd2phi = dau2.Phi();
      feept = ee.Pt();
      feemt = ee.Mt();
      feep = ee.P();
      feem = ee.M();
      feeeta = ee.Eta();
      feephi = ee.Phi();
      feephiv = PhiV(dau1, dau2);
      fmotherpt = mother.pt();
      fmotherm = sqrt(pow(mother.e(), 2) + pow(mother.p(), 2)); // run2: GetCalcMass() ??
      fmothermt = sqrt(pow(fmotherm, 2) + pow(fmotherpt, 2));
      fmotherp = mother.p();
      fmothereta = mother.eta();
      fmotherphi = mother.phi();
      fID = mother.pdgCode();

      // get multiplicity based weight:
      int iwbin = fhwMultpT->FindBin(fmotherpt);
      fwMultpT = fhwMultpT->GetBinContent(iwbin);   // pT weight
      fwMultpT2 = fhwMultpT2->GetBinContent(iwbin); // pT weight
      double min_mT = fhwMultmT->GetBinLowEdge(1);  // consider as minimum valid mT value the edge of the weight histo.
      if (fmothermt > min_mT) {
        iwbin = fhwMultmT->FindBin(fmothermt);
        fwMultmT = fhwMultmT->GetBinContent(iwbin);   // mT weight
        fwMultmT2 = fhwMultmT2->GetBinContent(iwbin); // mT weight
      } else {
        LOGP(error, "Generated particle with mT < Pion mass cannot be weighted");
        fwMultmT = 0.;
        fwMultmT2 = 0.;
      }

      // Which ALT weight to use?:
      Double_t fwALT = fwEffpT; // by default use pt efficiency weight
      if (fConfigALTweight == 1)
        fwALT = fwMultmT; // mT multiplicity weight
      if (fConfigALTweight == 11)
        fwALT = fwMultmT2; // mT multiplicity weight, higher mult
      if (fConfigALTweight == 2)
        fwALT = fwMultpT; // pT multiplicity weight
      if (fConfigALTweight == 22)
        fwALT = fwMultpT2; // pT multiplicity weight, higher mult

      // fill the tree
      if (fConfigWriteTTree) {
        tree(fd1DCA, fd2DCA, fpairDCA, fd1origpt, fd1origp, fd1origeta, fd1origphi, fd2origpt, fd2origp, fd2origeta, fd2origphi, fd1pt, fd1p, fd1eta, fd1phi, fd2pt, fd2p, fd2eta, fd2phi, feeorigpt, feeorigp, feeorigm, feeorigeta, feeorigphi, feeorigphiv, feept, feemt, feep, feem, feeeta, feephi, feephiv, fmotherpt, fmothermt, fmotherp, fmotherm, fmothereta, fmotherphi, fID, fdectyp, fdau3pdg, fweight, fwEffpT, fwMultpT, fwMultmT, fwMultpT2, fwMultmT2, fpass);
      }

      // fill the histograms
      if (fdectyp < 4) {                   // skip for the moment 4-particle decays
        for (Int_t jj = 0; jj < 3; jj++) { // fill the different hindex -> particles
          if (hindex[jj] > -1) {
            fmee_orig[hindex[jj]]->Fill(ee_orig.M(), fweight);
            if (fConfigALTweight == 1 || fConfigALTweight == 11) {
              fmotherpT_orig[hindex[jj]]->Fill(fmothermt, fweight);
            } else if (fConfigALTweight == 2 || fConfigALTweight == 22 || fConfigALTweight == 0) {
              fmotherpT_orig[hindex[jj]]->Fill(fmotherpt, fweight);
            }
            fpteevsmee_orig[hindex[jj]]->Fill(ee_orig.M(), ee.Pt(), fweight);
            fphi_orig[hindex[jj]]->Fill(ee_orig.Phi(), fweight);
            frap_orig[hindex[jj]]->Fill(ee_orig.Rapidity(), fweight);
            fmee_orig_wALT[hindex[jj]]->Fill(ee_orig.M(), fweight * fwALT);
            fpteevsmee_orig_wALT[hindex[jj]]->Fill(ee_orig.M(), ee.Pt(), fweight * fwALT);
            if (fConfigALTweight == 1 || fConfigALTweight == 11) {
              fmotherpT_orig_wALT[hindex[jj]]->Fill(fmothermt, fweight * fwALT);
            } else if (fConfigALTweight == 2 || fConfigALTweight == 22 || fConfigALTweight == 0) {
              fmotherpT_orig_wALT[hindex[jj]]->Fill(fmotherpt, fweight * fwALT);
            }
            if (fpass) {
              fmee[hindex[jj]]->Fill(ee.M(), fweight);
              fpteevsmee[hindex[jj]]->Fill(ee.M(), ee.Pt(), fweight);
              fphi[hindex[jj]]->Fill(ee.Phi(), fweight);
              frap[hindex[jj]]->Fill(ee.Rapidity(), fweight);
              registry.fill(HIST("DCAeevsmee"), ee.M(), fpairDCA, fweight);
              registry.fill(HIST("DCAeevsptee"), ee.Pt(), fpairDCA, fweight);
              fmee_wALT[hindex[jj]]->Fill(ee.M(), fweight * fwALT);
              fpteevsmee_wALT[hindex[jj]]->Fill(ee.M(), ee.Pt(), fweight * fwALT);
            }
          }
        }
      }

      // Virtual photon generation
      //-------------------------
      // We will generate one virtual photon per histogrammed pion
      if (mother.pdgCode() == 111) {
        // get mass and pt from histos and flat eta and phi
        Double_t VPHpT = ffVPHpT->GetRandom();
        Double_t VPHmass = fhKW->GetRandom();
        Double_t VPHeta = -1. + gRandom->Rndm() * 2.;
        Double_t VPHphi = 2.0 * TMath::ACos(-1.) * gRandom->Rndm();
        TLorentzVector beam;
        beam.SetPtEtaPhiM(VPHpT, VPHeta, VPHphi, VPHmass);
        Double_t decaymasses[2] = {(TDatabasePDG::Instance()->GetParticle(11))->Mass(), (TDatabasePDG::Instance()->GetParticle(11))->Mass()};
        TGenPhaseSpace VPHgen;
        Bool_t SetDecay;
        SetDecay = VPHgen.SetDecay(beam, 2, decaymasses);
        if (SetDecay == 0) {
          LOGP(error, "decay not permitted by kinematics");
        }
        Double_t VPHweight = VPHgen.Generate();
        // get electrons from the decay
        TLorentzVector *decay1, *decay2;
        decay1 = VPHgen.GetDecay(0);
        decay2 = VPHgen.GetDecay(1);
        dau1.SetPxPyPzE(decay1->Px(), decay1->Py(), decay1->Pz(), decay1->E());
        dau2.SetPxPyPzE(decay2->Px(), decay2->Py(), decay2->Pz(), decay2->E());

        // create dielectron before resolution effects:
        ee = dau1 + dau2;
        ee_orig = ee;

        // get index for histograms
        hindex[0] = nInputParticles - 1;
        hindex[1] = -1;
        hindex[2] = -1;

        // Fill tree words before resolution/acceptance
        fd1origpt = dau1.Pt();
        fd1origp = dau1.P();
        fd1origeta = dau1.Eta();
        fd1origphi = dau1.Phi();
        fd2origpt = dau2.Pt();
        fd2origp = dau2.P();
        fd2origeta = dau2.Eta();
        fd2origphi = dau2.Phi();
        feeorigpt = ee.Pt();
        feeorigp = ee.P();
        feeorigm = ee.M();
        feeorigeta = ee.Eta();
        feeorigphi = ee.Phi();
        feeorigphiv = PhiV(dau1, dau2);

        // get the efficiency weight
        Int_t effbin = fhwEffpT->FindBin(dau1.Pt());
        fwEffpT = fhwEffpT->GetBinContent(effbin);
        effbin = fhwEffpT->FindBin(dau2.Pt());
        fwEffpT = fwEffpT * fhwEffpT->GetBinContent(effbin);

        // Resolution and acceptance
        //-------------------------
        dau1 = ApplyResolution(dau1, 1, fConfigResolType);
        dau2 = ApplyResolution(dau2, -1, fConfigResolType);
        fpass = true;
        if (dau1.Pt() < fConfigMinPt || dau2.Pt() < fConfigMinPt)
          fpass = false; // leg pT cut
        if (dau1.Pt() > fConfigMaxPt || dau2.Pt() > fConfigMaxPt)
          fpass = false; // leg pT cut
        if (dau1.Vect().Unit().Dot(dau2.Vect().Unit()) > TMath::Cos(fConfigMinOpAng))
          fpass = false; // opening angle cut
        if (TMath::Abs(dau1.Eta()) > fConfigMaxEta || TMath::Abs(dau2.Eta()) > fConfigMaxEta)
          fpass = false;

        fpairDCA = 10000.; // ??

        // Fill tree words after resolution/acceptance
        ee = dau1 + dau2;
        fd1pt = dau1.Pt();
        fd1p = dau1.P();
        fd1eta = dau1.Eta();
        fd1phi = dau1.Phi();
        fd2pt = dau2.Pt();
        fd2p = dau2.P();
        fd2eta = dau2.Eta();
        fd2phi = dau2.Phi();
        feept = ee.Pt();
        feemt = ee.Mt();
        feep = ee.P();
        feem = ee.M();
        feeeta = ee.Eta();
        feephi = ee.Phi();
        feephiv = PhiV(dau1, dau2);
        fmotherpt = beam.Pt();
        fmothermt = sqrt(pow(beam.M(), 2) + pow(beam.Pt(), 2));
        fmotherp = beam.P();
        fmotherm = beam.M();
        fmothereta = beam.Eta();
        fmotherphi = beam.Phi();
        fID = 0; // set ID to Zero for VPH
        fweight = VPHweight;
        // get multiplicity based weight:
        fwMultmT = 1; // no weight for photons so far

        // Fill the tree
        if (fConfigWriteTTree) { // many parameters not set for photons: d1DCA,fd2DCA, fdectyp,fdau3pdg,fwMultpT,fwMultpT2,fwMultmT2
          tree(fd1DCA, fd2DCA, fpairDCA, fd1origpt, fd1origp, fd1origeta, fd1origphi, fd2origpt, fd2origp, fd2origeta, fd2origphi, fd1pt, fd1p, fd1eta, fd1phi, fd2pt, fd2p, fd2eta, fd2phi, feeorigpt, feeorigp, feeorigm, feeorigeta, feeorigphi, feeorigphiv, feept, feemt, feep, feem, feeeta, feephi, feephiv, fmotherpt, fmothermt, fmotherp, fmotherm, fmothereta, fmotherphi, fID, fdectyp, fdau3pdg, fweight, fwEffpT, fwMultpT, fwMultmT, fwMultpT2, fwMultmT2, fpass);
        }

        // Fill the histograms
        for (Int_t jj = 0; jj < 3; jj++) { // fill the different hindex -> particles
          if (hindex[jj] > -1) {
            fmee_orig[hindex[jj]]->Fill(ee_orig.M(), VPHweight);
            fpteevsmee_orig[hindex[jj]]->Fill(ee_orig.M(), ee.Pt(), VPHweight);
            fphi_orig[hindex[jj]]->Fill(ee_orig.Phi(), VPHweight);
            frap_orig[hindex[jj]]->Fill(ee_orig.Rapidity(), VPHweight);
            if (fpass) {
              fmee[hindex[jj]]->Fill(ee.M(), VPHweight);
              fpteevsmee[hindex[jj]]->Fill(ee.M(), ee.Pt(), VPHweight);
              fphi[hindex[jj]]->Fill(ee.Phi(), VPHweight);
              frap[hindex[jj]]->Fill(ee.Rapidity(), VPHweight);
            }
          }
        }
      }
    }
  }

  Double_t PhiV(PxPyPzEVector e1, PxPyPzEVector e2)
  {
    Double_t outPhiV;
    XYZVector p1 = e1.Vect();
    XYZVector p2 = e2.Vect();
    XYZVector p12 = p1 + p2;
    XYZVector u = p12.Unit();
    XYZVector p1u = p1.Unit();
    XYZVector p2u = p2.Unit();
    XYZVector v = p1u.Cross(p2u);
    XYZVector w = u.Cross(v);
    XYZVector zu(0, 0, 1);
    XYZVector wc = u.Cross(zu);
    outPhiV = TMath::ACos(wc.Unit().Dot(w.Unit()));
    return outPhiV;
  }

  void SetHistograms()
  {

    AxisSpec ptAxis = {fConfigNBinsPtee, fConfigMinPtee, fConfigMaxPtee, "#it{p}_{T,ee} (GeV/c)"};
    AxisSpec mAxis = {fConfigNBinsMee, fConfigMinMee, fConfigMaxMee, "#it{m}_{ee} (GeV/c^{2})"};
    AxisSpec phiAxis = {fConfigNBinsPhi, fConfigMinPhi, fConfigMaxPhi, "#it{phi}_{ee}"};
    AxisSpec rapAxis = {fConfigNBinsRap, fConfigMinRap, fConfigMaxRap, "#it{y}_{ee}"};

    registry.add<TH1>("NEvents", "NEvents", HistType::kTH1F, {{1, 0, 1}}, false);

    registry.add<TH2>("ULS", "ULS", HistType::kTH2F, {mAxis, ptAxis}, true);
    registry.add<TH2>("LSpp", "LSpp", HistType::kTH2F, {mAxis, ptAxis}, true);
    registry.add<TH2>("LSmm", "LSmm", HistType::kTH2F, {mAxis, ptAxis}, true);

    registry.add<TH2>("ULS_orig", "ULS_orig", HistType::kTH2F, {mAxis, ptAxis}, true);
    registry.add<TH2>("LSpp_orig", "LSpp_orig", HistType::kTH2F, {mAxis, ptAxis}, true);
    registry.add<TH2>("LSmm_orig", "LSmm_orig", HistType::kTH2F, {mAxis, ptAxis}, true);

    registry.add<TH2>("DCAeevsmee", "DCAeevsmee", HistType::kTH2F, {{fConfigMBins, "#it{m}_{ee} (GeV/c^{2})"}, {fConfigDCABins, "DCA_{xy}^{ee} (cm)"}}, true);
    registry.add<TH2>("DCAeevsptee", "DCAeevsptee", HistType::kTH2F, {{fConfigPtBins, "#it{p}_{T,ee} (GeV/c)"}, {fConfigDCABins, "DCA_{xy}^{ee} (cm)"}}, true);

    for (auto& particle : fParticleListNames) {
      fmee.push_back(registry.add<TH1>(Form("mee_%s", particle.Data()), Form("mee_%s", particle.Data()), HistType::kTH1F, {mAxis}, true));
      fmee_orig.push_back(registry.add<TH1>(Form("mee_orig_%s", particle.Data()), Form("mee_orig_%s", particle.Data()), HistType::kTH1F, {mAxis}, true));
      fmotherpT_orig.push_back(registry.add<TH1>(Form("mohterpT_orig_%s", particle.Data()), Form("motherpT_orig_%s", particle.Data()), HistType::kTH1F, {ptAxis}, true));
      fphi.push_back(registry.add<TH1>(Form("phi_%s", particle.Data()), Form("phi_%s", particle.Data()), HistType::kTH1F, {phiAxis}, true));
      fphi_orig.push_back(registry.add<TH1>(Form("phi_orig_%s", particle.Data()), Form("phi_orig_%s", particle.Data()), HistType::kTH1F, {phiAxis}, true));
      frap.push_back(registry.add<TH1>(Form("rap_%s", particle.Data()), Form("rap_%s", particle.Data()), HistType::kTH1F, {rapAxis}, true));
      frap_orig.push_back(registry.add<TH1>(Form("rap_orig_%s", particle.Data()), Form("rap_orig_%s", particle.Data()), HistType::kTH1F, {rapAxis}, true));
      fpteevsmee.push_back(registry.add<TH2>(Form("pteevsmee_%s", particle.Data()), Form("pteevsmee_%s", particle.Data()), HistType::kTH2F, {mAxis, ptAxis}, true));
      fpteevsmee_orig.push_back(registry.add<TH2>(Form("pteevsmee_orig_%s", particle.Data()), Form("pteevsmee_orig_%s", particle.Data()), HistType::kTH2F, {mAxis, ptAxis}, true));
      fmee_wALT.push_back(registry.add<TH1>(Form("mee_wALT_%s", particle.Data()), Form("mee_wALT_%s", particle.Data()), HistType::kTH1F, {mAxis}, true));
      fmee_orig_wALT.push_back(registry.add<TH1>(Form("mee_orig_wALT_%s", particle.Data()), Form("mee_orig_wALT_%s", particle.Data()), HistType::kTH1F, {mAxis}, true));
      fmotherpT_orig_wALT.push_back(registry.add<TH1>(Form("mohterpT_orig_wALT_%s", particle.Data()), Form("motherpT_orig_wALT%s", particle.Data()), HistType::kTH1F, {ptAxis}, true));
      fpteevsmee_wALT.push_back(registry.add<TH2>(Form("pteevsmee_wALT%s", particle.Data()), Form("pteevsmee_wALT_%s", particle.Data()), HistType::kTH2F, {mAxis, ptAxis}, true));
      fpteevsmee_orig_wALT.push_back(registry.add<TH2>(Form("pteevsmee_orig_wALT%s", particle.Data()), Form("pteevsmee_orig_wALT_%s", particle.Data()), HistType::kTH2F, {mAxis, ptAxis}, true));
    }

    fmee.push_back(registry.add<TH1>("mee", "mee", HistType::kTH1F, {mAxis}, true));
    fmee_orig.push_back(registry.add<TH1>("mee_orig", "mee_orig", HistType::kTH1F, {mAxis}, true));
    fmotherpT_orig.push_back(registry.add<TH1>("motherpT_orig", "motherpT_orig", HistType::kTH1F, {ptAxis}, true));
    fphi.push_back(registry.add<TH1>("phi", "phi", HistType::kTH1F, {phiAxis}, true));
    fphi_orig.push_back(registry.add<TH1>("phi_orig", "phi_orig", HistType::kTH1F, {phiAxis}, true));
    frap.push_back(registry.add<TH1>("rap", "rap", HistType::kTH1F, {rapAxis}, true));
    frap_orig.push_back(registry.add<TH1>("rap_orig", "rap_orig", HistType::kTH1F, {rapAxis}, true));
    fpteevsmee.push_back(registry.add<TH2>("pteevsmee", "pteevsmee", HistType::kTH2F, {mAxis, ptAxis}, true));
    fpteevsmee_orig.push_back(registry.add<TH2>("pteevsmee_orig", "pteevsmee_orig", HistType::kTH2F, {mAxis, ptAxis}, true));
    fmee_wALT.push_back(registry.add<TH1>("mee_wALT", "mee_wALT", HistType::kTH1F, {mAxis}, true));
    fmee_orig_wALT.push_back(registry.add<TH1>("mee_orig_wALT", "mee_orig_wALT", HistType::kTH1F, {mAxis}, true));
    fmotherpT_orig_wALT.push_back(registry.add<TH1>("motherpT_orig_wALT", "motherpT_orig_wALT", HistType::kTH1F, {ptAxis}, true));
    fpteevsmee_wALT.push_back(registry.add<TH2>("pteevsmee_wALT", "pteevsmee_wALT", HistType::kTH2F, {mAxis, ptAxis}, true));
    fpteevsmee_orig_wALT.push_back(registry.add<TH2>("pteevsmee_orig_wALT", "pteevsmee_orig_wALT", HistType::kTH2F, {mAxis, ptAxis}, true));
  }

  PxPyPzEVector ApplyResolution(PxPyPzEVector vec, Char_t ch = 0, Int_t Run = 2)
  {

    Double_t theta, phi, pt, p, px, py, pz, E, mass, eta;
    PxPyPzEVector resvec;

    mass = 0.51099906e-3;
    pt = vec.Pt();
    p = vec.P();
    theta = vec.Theta();
    phi = vec.Phi();
    eta = vec.Eta();

    if (Run == 0) {
      resvec = vec;
    } else if (Run == 1) {
      TH1D* hisSlice(nullptr);
      if (fArr) {
        TH2D* hDeltaPtvsPt = static_cast<TH2D*>(fArr->At(0));
        // Get the momentum slice histogram for sampling of the smearing
        // (in some input file versions this also contains a landau fit )
        // since histogram bins start at 1, we can use the bin number directly for the array index (first slice stored in position 1).
        Int_t histIndex = hDeltaPtvsPt->GetXaxis()->FindBin(pt);
        if (histIndex < 1)
          histIndex = 1; // in case some track is below the first p-bin (which currently starts at 100 MeV).
        if (histIndex > fArr->GetLast()) {
          histIndex = fArr->GetLast();
        }
        hisSlice = static_cast<TH1D*>(fArr->At(histIndex));
      }
      // get smear parameter via random selection from the p slices retreived from the deltaP/P plot
      Double_t SmearingPPercent(0.);
      if (hisSlice) {
        SmearingPPercent = hisSlice->GetRandom();
      } else {
        SmearingPPercent = gRandom->Gaus(0, p * sqrt(0.004 * 0.004 + (0.012 * p) * (0.012 * p))) / p; // for B=0.5 Tesla
      }

      Double_t SmearingP = p * SmearingPPercent;
      p -= SmearingP;
      px = p * sin(theta) * cos(phi);
      py = p * sin(theta) * sin(phi);
      pz = p * cos(theta);
      E = sqrt(p * p + mass * mass);

      resvec.SetPxPyPzE(px, py, pz, E);

    } else if (Run == 2) {
      // smear pt
      Int_t ptbin = reinterpret_cast<TH2D*>(fArrResoPt->At(0))->GetXaxis()->FindBin(pt);
      if (ptbin < 1)
        ptbin = 1;
      if (ptbin > fArrResoPt->GetLast())
        ptbin = fArrResoPt->GetLast();
      Double_t smearing = reinterpret_cast<TH1D*>(fArrResoPt->At(ptbin))->GetRandom() * pt;
      Double_t sPt = pt - smearing;

      // smear eta
      ptbin = reinterpret_cast<TH2D*>(fArrResoEta->At(0))->GetXaxis()->FindBin(pt);
      if (ptbin < 1)
        ptbin = 1;
      if (ptbin > fArrResoEta->GetLast())
        ptbin = fArrResoEta->GetLast();
      smearing = reinterpret_cast<TH1D*>(fArrResoEta->At(ptbin))->GetRandom();
      Double_t sEta = eta - smearing;

      // smear phi
      ptbin = reinterpret_cast<TH2D*>(fArrResoPhi_Pos->At(0))->GetXaxis()->FindBin(pt);
      if (ptbin < 1)
        ptbin = 1;
      if (ptbin > fArrResoPhi_Pos->GetLast())
        ptbin = fArrResoPhi_Pos->GetLast();
      if (ch > 0) {
        smearing = reinterpret_cast<TH1D*>(fArrResoPhi_Pos->At(ptbin))->GetRandom();
      } else if (ch < 0) {
        smearing = reinterpret_cast<TH1D*>(fArrResoPhi_Neg->At(ptbin))->GetRandom();
      }
      Double_t sPhi = phi - smearing;

      Double_t sPx = sPt * cos(sPhi);
      Double_t sPy = sPt * sin(sPhi);
      Double_t sPz = sPt * sinh(sEta);
      Double_t sP = sPt * cosh(sEta);
      Double_t sE = sqrt(sP * sP + mass * mass);

      resvec.SetPxPyPzE(sPx, sPy, sPz, sE);
    }

    return resvec;
  }

  void GetEffHisto(TString filename, TString histname)
  {
    // get efficiency histo
    LOGP(detail, "Set Efficiency histo");
    TString fFileName = filename;
    TString fFileNameLocal = filename;
    // Get Efficiency
    if (fFileName.Contains("alien")) { // TODO
      LOGP(error, "Copying from alien not yet implemented");
    }
    // Get Efficiency weight file:
    TFile* fFile = TFile::Open(fFileNameLocal.Data());
    if (!fFile) {
      LOGP(error, "Could not open Efficiency file {}", fFileNameLocal.Data());
      return;
    }
    if (fFile->GetListOfKeys()->Contains(histname.Data())) {
      fhwEffpT = reinterpret_cast<TH1F*>(fFile->Get(histname.Data())); // histo: eff weight in function of pT.
      fhwEffpT->SetDirectory(nullptr);
    } else {
      LOGP(error, "Could not open histogram {} from file {}", histname.Data(), fFileNameLocal.Data());
    }

    fFile->Close();
  }

  void GetResHisto(TString filename, TString ptHistName, TString etaHistName, TString phiPosHistName, TString phiNegHistName)
  {
    // get resolutoin histo
    LOGP(detail, "Set Resolution histo");
    TString fFileName = filename;
    TString fFileNameLocal = filename;
    // Get Resolution map
    if (fFileName.Contains("alien")) { // TODO
      LOGP(error, "Copying from alien not yet implemented");
    }
    // Get file:
    TFile* fFile = TFile::Open(fFileNameLocal.Data());
    if (!fFile) {
      LOGP(error, "Could not open Resolution file {}", fFileNameLocal.Data());
      return;
    }
    TObjArray* ArrResoPt = nullptr;
    if (fFile->GetListOfKeys()->Contains(ptHistName.Data())) {
      ArrResoPt = reinterpret_cast<TObjArray*>(fFile->Get(ptHistName.Data()));
    } else {
      LOGP(error, "Could not open {} from file {}", ptHistName.Data(), fFileNameLocal.Data());
    }

    TObjArray* ArrResoEta = nullptr;
    if (fFile->GetListOfKeys()->Contains(etaHistName.Data())) {
      ArrResoEta = reinterpret_cast<TObjArray*>(fFile->Get(etaHistName.Data()));
    } else {
      LOGP(error, "Could not open {} from file {}", etaHistName.Data(), fFileNameLocal.Data());
    }

    TObjArray* ArrResoPhi_Pos = nullptr;
    if (fFile->GetListOfKeys()->Contains(phiPosHistName.Data())) {
      ArrResoPhi_Pos = reinterpret_cast<TObjArray*>(fFile->Get(phiPosHistName.Data()));
    } else {
      LOGP(error, "Could not open {} from file {}", phiPosHistName.Data(), fFileNameLocal.Data());
    }

    TObjArray* ArrResoPhi_Neg = nullptr;
    if (fFile->GetListOfKeys()->Contains(phiNegHistName.Data())) {
      ArrResoPhi_Neg = reinterpret_cast<TObjArray*>(fFile->Get(phiNegHistName.Data()));
    } else {
      LOGP(error, "Could not open {} from file {}", phiNegHistName.Data(), fFileNameLocal.Data());
    }
    fArrResoPt = ArrResoPt;
    fArrResoEta = ArrResoEta;
    fArrResoPhi_Pos = ArrResoPhi_Pos;
    fArrResoPhi_Neg = ArrResoPhi_Neg;
    fFile->Close();
  }

  void GetResHisto(TString filename, TString pHistName)
  {
    // get resolutoin histo
    LOGP(detail, "Set Resolution histo\n");
    TString fFileName = filename;
    TString fFileNameLocal = filename;
    // Get Resolution map
    if (fFileName.Contains("alien")) { // TODO
      LOGP(error, "copying from alien not yet implemented");
    }
    // Get file:
    TFile* fFile = TFile::Open(fFileNameLocal.Data());
    if (!fFile) {
      LOGP(error, "Could not open Resolution file {}", fFileNameLocal.Data());
      LOGP(error, "No resolution array set! Using internal parametrization");
      return;
    }
    TObjArray* ArrResoP = 0x0;
    if (fFile->GetListOfKeys()->Contains(pHistName.Data())) {
      ArrResoP = reinterpret_cast<TObjArray*>(fFile->Get(pHistName.Data()));
    } else {
      LOGP(error, "Could not open {} from file {}", pHistName.Data(), fFileNameLocal.Data());
      LOGP(error, "No resolution array set! Using internal parametrization");
    }

    fArr = ArrResoP;
    fFile->Close();
  }

  void GetDCATemplates(TString filename, TString histname)
  {
    // get dca tamplates
    LOGP(detail, "Set DCA templates");
    TString fFileName = filename;
    TString fFileNameLocal = filename;
    if (fFileName.Contains("alien")) { // TODO
      LOGP(error, "copying from alien not yet implemented");
    }
    // Get  file:
    TFile* fFile = TFile::Open(fFileNameLocal.Data());
    if (!fFile) {
      LOGP(error, "Could not open DCATemplate file {}", fFileNameLocal.Data());
      return;
    }
    fh_DCAtemplates = new TH1F*[nbDCAtemplate];
    for (int jj = 0; jj < nbDCAtemplate; jj++) {
      if (fFile->GetListOfKeys()->Contains(Form("%s%d", histname.Data(), jj + 1))) {
        fh_DCAtemplates[jj] = reinterpret_cast<TH1F*>(fFile->Get(Form("%s%d", histname.Data(), jj + 1)));
      } else {
        LOGP(error, "Could not open {}{} from file {}", histname.Data(), jj + 1, fFileNameLocal.Data());
      }
    }
    fFile->Close();
  }

  void GetMultHisto(TString filename, TString histnamept, TString histnamept2, TString histnamemt, TString histnamemt2)
  {
    // get multiplicity weights
    LOGP(detail, "Set Multiplicity weight files");
    TString fFileName = filename;
    TString fFileNameLocal = filename;
    if (fFileName.Contains("alien")) { // TODO
      LOGP(error, "copying from alien not yet implemented");
    }
    // Get  file:
    TFile* fFile = TFile::Open(fFileNameLocal.Data());
    if (!fFile) {
      LOGP(error, "Could not open Multiplicity weight file {}", fFileNameLocal.Data());
      return;
    }

    if (fFile->GetListOfKeys()->Contains(histnamept.Data())) {
      fhwMultpT = reinterpret_cast<TH1F*>(fFile->Get(histnamept.Data())); // histo: multiplicity weight in function of pT.
    } else {
      LOGP(error, "Could not open {} from file {}", histnamept.Data(), fFileNameLocal.Data());
    }

    if (fFile->GetListOfKeys()->Contains(histnamemt.Data())) {
      fhwMultmT = reinterpret_cast<TH1F*>(fFile->Get(histnamemt.Data())); // histo: multiplicity weight in function of mT.
    } else {
      LOGP(error, "Could not open {} from file {}", histnamemt.Data(), fFileNameLocal.Data());
    }

    if (fFile->GetListOfKeys()->Contains(histnamept2.Data())) {
      fhwMultpT2 = reinterpret_cast<TH1F*>(fFile->Get(histnamept2.Data())); // histo: multiplicity weight in function of pT.
    } else {
      LOGP(error, "Could not open {} from file {}", histnamept2.Data(), fFileNameLocal.Data());
    }

    if (fFile->GetListOfKeys()->Contains(histnamemt2.Data())) {
      fhwMultmT2 = reinterpret_cast<TH1F*>(fFile->Get(histnamemt2)); // histo: multiplicity weight in function of mT.
    } else {
      LOGP(error, "Could not open {} from file {}", histnamemt2.Data(), fFileNameLocal.Data());
    }
    fFile->Close();
  }

  void GetPhotonPtParametrization(TString filename, TString dirname, TString funcname)
  {
    LOGP(detail, "Set photon parametrization");
    TString fFileName = filename;
    TString fFileNameLocal = filename;
    if (fFileName.Contains("alien")) { // TODO
      LOGP(error, "copying from alien not yet implemented");
    }
    // Get file:
    TFile* fFile = TFile::Open(fFileNameLocal.Data());
    if (!fFile) {
      LOGP(error, "Could not open photon parametrization from file {}", fFileNameLocal.Data());
      return;
    }
    bool good = false;
    if (fFile->GetListOfKeys()->Contains(dirname.Data())) {
      TDirectory* dir = fFile->GetDirectory(dirname.Data());
      if (dir->GetListOfKeys()->Contains(funcname.Data())) {
        ffVPHpT = reinterpret_cast<TF1*>(dir->Get(funcname.Data()));
        good = true;
      }
    }
    if (!good) {
      LOGP(error, "Could not open photon parametrization {}/{} from file {}", dirname.Data(), funcname.Data(), fFileNameLocal.Data());
    }
    fFile->Close();
  }

  void fillKrollWada()
  {
    // Build Kroll-wada for virtual photon mass parametrization:
    Double_t KWmass = 0.;
    Double_t emass = (TDatabasePDG::Instance()->GetParticle(11))->Mass();
    // Int_t KWnbins = 10000;
    Float_t KWmin = 2. * emass;
    // Float_t KWmax         = 1.1;
    Double_t KWbinwidth = (fConfigKWMax - KWmin) / (Double_t)fConfigKWNBins;
    fhKW = new TH1F("fhKW", "fhKW", fConfigKWNBins, KWmin, fConfigKWMax);
    for (Int_t ibin = 1; ibin <= fConfigKWNBins; ibin++) {
      KWmass = KWmin + (Double_t)(ibin - 1) * KWbinwidth + KWbinwidth / 2.0;
      fhKW->AddBinContent(ibin, 2. * (1. / 137.03599911) / 3. / 3.14159265359 / KWmass * sqrt(1. - 4. * emass * emass / KWmass / KWmass) * (1. + 2. * emass * emass / KWmass / KWmass));
    }
  }

  // ULS and LS spectra
  Partition<aod::McParticles> Electrons = (aod::mcparticle::pdgCode == 11) && (aod::mcparticle::statusCode == 1);
  Partition<aod::McParticles> Positrons = (aod::mcparticle::pdgCode == -11) && (aod::mcparticle::statusCode == 1);

  void processPairing(aod::McCollision const& collision, aod::McParticles const& mcParticles)
  {

    registry.fill(HIST("NEvents"), 0.5);

    auto const electronsGrouped = Electrons->sliceBy(perCollision, collision.globalIndex());
    auto const positronsGrouped = Positrons->sliceBy(perCollision, collision.globalIndex());

    // ULS spectrum
    for (auto const& [particle1, particle2] : combinations(o2::soa::CombinationsFullIndexPolicy(electronsGrouped, positronsGrouped))) {
      PxPyPzEVector p1, p2, dielectron;
      p1.SetPxPyPzE(particle1.px(), particle1.py(), particle1.pz(), particle1.e());
      p2.SetPxPyPzE(particle2.px(), particle2.py(), particle2.pz(), particle2.e());
      dielectron = p1 + p2;
      registry.fill(HIST("ULS_orig"), dielectron.M(), dielectron.Pt(), particle1.weight() * particle2.weight());
      if (p1.Pt() > fConfigMinPt && p2.Pt() > fConfigMinPt && p1.Pt() < fConfigMaxPt && p2.Pt() < fConfigMaxPt && TMath::Abs(p1.Eta()) < fConfigMaxEta && TMath::Abs(p2.Eta()) < fConfigMaxEta && p1.Vect().Unit().Dot(p2.Vect().Unit()) < TMath::Cos(fConfigMinOpAng)) {
        registry.fill(HIST("ULS"), dielectron.M(), dielectron.Pt(), particle1.weight() * particle2.weight());
      }
    }

    // LS spectra
    for (auto& [particle1, particle2] : combinations(o2::soa::CombinationsStrictlyUpperIndexPolicy(electronsGrouped, electronsGrouped))) {
      PxPyPzEVector p1, p2, dielectron;
      p1.SetPxPyPzE(particle1.px(), particle1.py(), particle1.pz(), particle1.e());
      p2.SetPxPyPzE(particle2.px(), particle2.py(), particle2.pz(), particle2.e());
      dielectron = p1 + p2;
      registry.fill(HIST("LSmm_orig"), dielectron.M(), dielectron.Pt(), particle1.weight() * particle2.weight());
      if (p1.Pt() > fConfigMinPt && p2.Pt() > fConfigMinPt && p1.Pt() < fConfigMaxPt && p2.Pt() < fConfigMaxPt && TMath::Abs(p1.Eta()) < fConfigMaxEta && TMath::Abs(p2.Eta()) < fConfigMaxEta && p1.Vect().Unit().Dot(p2.Vect().Unit()) < TMath::Cos(fConfigMinOpAng)) {
        registry.fill(HIST("LSmm"), dielectron.M(), dielectron.Pt(), particle1.weight() * particle2.weight());
      }
    }
    for (auto& [particle1, particle2] : combinations(o2::soa::CombinationsStrictlyUpperIndexPolicy(positronsGrouped, positronsGrouped))) {
      PxPyPzEVector p1, p2, dielectron;
      p1.SetPxPyPzE(particle1.px(), particle1.py(), particle1.pz(), particle1.e());
      p2.SetPxPyPzE(particle2.px(), particle2.py(), particle2.pz(), particle2.e());
      dielectron = p1 + p2;
      registry.fill(HIST("LSpp_orig"), dielectron.M(), dielectron.Pt(), particle1.weight() * particle2.weight());
      if (p1.Pt() > fConfigMinPt && p2.Pt() > fConfigMinPt && p1.Pt() < fConfigMaxPt && p2.Pt() < fConfigMaxPt && TMath::Abs(p1.Eta()) < fConfigMaxEta && TMath::Abs(p2.Eta()) < fConfigMaxEta && p1.Vect().Unit().Dot(p2.Vect().Unit()) < TMath::Cos(fConfigMinOpAng)) {
        registry.fill(HIST("LSpp"), dielectron.M(), dielectron.Pt(), particle1.weight() * particle2.weight());
      }
    }
  }

  PROCESS_SWITCH(lmeelfcocktail, processPairing, "Process pairing", true);
  PROCESS_SWITCH(lmeelfcocktail, processCocktail, "Process cocktail", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<lmeelfcocktail>(cfgc, TaskName{"em-lmee-lf-cocktail"})};
}
