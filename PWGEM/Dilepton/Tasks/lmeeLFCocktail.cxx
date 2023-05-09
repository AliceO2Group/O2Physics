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
#include "Framework/Task.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/Logger.h"
#include "SimulationDataFormat/MCTrack.h"
#include "Math/Vector4D.h"
#include "Math/Vector3D.h"
#include "TFile.h"
#include "TF1.h"
#include "TRandom.h"
#include "TDatabasePDG.h"
#include "TGenPhaseSpace.h"
#include "TGrid.h"
#include "TTree.h"
#include <nlohmann/json.hpp>

using namespace o2::framework;
using namespace ROOT::Math;

struct eeTTree {
  float fd1DCA;
  float fd2DCA;
  float fpairDCA;
  float fd1origpt;
  float fd1origp;
  float fd1origeta;
  float fd1origphi;
  float fd2origpt;
  float fd2origp;
  float fd2origeta;
  float fd2origphi;
  float fd1pt;
  float fd1p;
  float fd1eta;
  float fd1phi;
  float fd2pt;
  float fd2p;
  float fd2eta;
  float fd2phi;
  float feeorigpt;
  float feeorigp;
  float feeorigm;
  float feeorigeta;
  float feeorigphi;
  float feeorigphiv;
  float feept;
  float feemt;
  float feep;
  float feem;
  float feeeta;
  float feephi;
  float feephiv;
  float fmotherpt;
  float fmothermt;
  float fmotherp;
  float fmotherm;
  float fmothereta;
  float fmotherphi;
  int fID;
  int fdectyp;
  int fdau3pdg;
  float fweight;
  float fwEffpT;
  float fwMultpT;
  float fwMultmT;
  float fwMultpT2;
  float fwMultmT2;
  bool fpass;
  float feeorigrap; // only in histogram, not in tree?
  float feerap;     // only in histogram, not in tree?
};

struct lmeelfcocktail {
  OutputObj<TTree> tree{"eeTTree"};

  HistogramRegistry registry{"registry", {}};
  Int_t nInputParticles = 17;
  std::vector<TString> fParticleListNames = {"Pi0", "Eta", "EtaP", "EtaP_dalitz_photon", "EtaP_dalitz_omega", "Rho", "Omega", "Omega_2body", "Omega_dalitz", "Phi", "Phi_2body", "Phi_dalitz_eta", "Phi_dalitz_pi0", "Jpsi", "Jpsi_2body", "Jpsi_radiative", "Virtual_Photon"};
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

  eeTTree treeWords;

  std::vector<double> DCATemplateEdges;
  int nbDCAtemplate;
  TH1F** fh_DCAtemplates;

  Configurable<int> fCollisionSystem{"cfgCollisionSystem", 200, "set the collision system"};
  Configurable<bool> fConfigWriteTTree{"cfgWriteTTree", false, "write tree output"};
  Configurable<bool> fConfigDoPairing{"cfgDoPairing", true, "do like and unlike sign pairing"};
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
  Configurable<bool> fConfigDoVirtPh{"cfgDoVirtPh", false, "generate one virt. photon for each pion"};
  Configurable<std::string> fConfigPhotonPtFileName{"cfgPhotonPtFileName", "", "file name for photon pT parametrization"};
  Configurable<std::string> fConfigPhotonPtDirName{"cfgPhotonPtDirName", "", "directory name for photon pT parametrization"};
  Configurable<std::string> fConfigPhotonPtFuncName{"cfgPhotonPtFuncName", "111_pt", "function name for photon pT parametrization"};

  // Configurable axes crashed the task. Take them out for the moment
  // ConfigurableAxis fConfigPtBins{"cfgPtBins", {VARIABLE_WIDTH, 0., 0.5, 1, 1.5, 2., 2.5, 3., 3.5, 4., 4.5, 5., 5.5, 6., 6.5, 7., 7.5, 8.}, "pT bins"};
  // ConfigurableAxis fConfigMBins{"cfgMBins", {VARIABLE_WIDTH, 0., 0.08, 0.14, 0.2, 1.1, 2.7, 2.8, 3.2, 5.0}, "mee bins"};
  // ConfigurableAxis fConfigDCABins{"cfgDCABins", {VARIABLE_WIDTH, 0., 0.4, 0.8, 1.2, 1.6, 2.0, 2.4, 3., 4., 5., 7., 10.}, "DCA bins"};

  Configurable<std::vector<double>> fConfigDCATemplateEdges{"cfgDCATemplateEdges", {0., .3, .4, .6, 1., 2.}, "DCA template edges"};

  void init(o2::framework::InitContext& ic)
  {
    if (fConfigWriteTTree) {
      SetTree();
    }
    SetHistograms();
    DCATemplateEdges = fConfigDCATemplateEdges;
    nbDCAtemplate = DCATemplateEdges.size();
    DCATemplateEdges.push_back(10000000.);

    if ((TString(fConfigEffFileName).BeginsWith("alien://") && TString(fConfigEffFileName).EndsWith(".root")) || (TString(fConfigResFileName).BeginsWith("alien://") && TString(fConfigResFileName).EndsWith(".root")) || (TString(fConfigDCAFileName).BeginsWith("alien://") && TString(fConfigDCAFileName).EndsWith(".root")) || (TString(fConfigMultFileName).BeginsWith("alien://") && TString(fConfigMultFileName).EndsWith(".root")) || (TString(fConfigPhotonPtFileName).BeginsWith("alien://") && TString(fConfigPhotonPtFileName).EndsWith(".root"))) {
      LOGP(info, "Connecting to grid via TGrid");
      TGrid::Connect("alien://");
    }

    GetEffHisto(TString(fConfigEffFileName), TString(fConfigEffHistName));
    if (fConfigResolType == 1) {
      GetResHisto(TString(fConfigResFileName), TString(fConfigResPHistName));
    } else if (fConfigResolType == 2) {
      GetResHisto(TString(fConfigResFileName), TString(fConfigResPtHistName), TString(fConfigResEtaHistName), TString(fConfigResPhiPosHistName), TString(fConfigResPhiNegHistName));
    }
    GetDCATemplates(TString(fConfigDCAFileName), TString(fConfigDCAHistName));
    GetMultHisto(TString(fConfigMultFileName), TString(fConfigMultHistPtName), TString(fConfigMultHistPt2Name), TString(fConfigMultHistMtName), TString(fConfigMultHistMt2Name));
    if (fConfigDoVirtPh) {
      GetPhotonPtParametrization(TString(fConfigPhotonPtFileName), TString(fConfigPhotonPtDirName), TString(fConfigPhotonPtFuncName));
    }
    fillKrollWada();
  }

  void run(o2::framework::ProcessingContext& pc)
  {
    // get the tracks
    auto mctracks = pc.inputs().get<std::vector<o2::MCTrack>>("mctracks");
    registry.fill(HIST("NEvents"), 0.5);

    std::vector<PxPyPzEVector> eBuff;
    std::vector<Char_t> echBuff;
    std::vector<Double_t> eweightBuff;

    bool skipNext = false;

    int trackID = -1;
    //  Loop over all MC particle
    for (auto& mctrack : mctracks) {
      trackID++;
      if (o2::mcgenstatus::getHepMCStatusCode(mctrack.getStatusCode()) != 1)
        continue;
      if (abs(mctrack.GetPdgCode()) == 11) {
        // get the electron
        //---------------
        if (fConfigDoPairing) {
          // LS and ULS spectra
          PxPyPzEVector e, dielectron;
          Char_t ech, dielectron_ch;
          Double_t eweight, dielectron_weight;
          e.SetPxPyPzE(mctrack.Px(), mctrack.Py(), mctrack.Pz(),
                       mctrack.GetEnergy());
          if (mctrack.GetPdgCode() > 0) {
            ech = 1.;
          } else {
            ech = -1.;
          }
          eweight = mctrack.getWeight();
          // put in the buffer
          //-----------------
          eBuff.push_back(e);
          echBuff.push_back(ech);
          eweightBuff.push_back(eweight);
          // loop the buffer and pair
          //------------------------
          for (Int_t jj = eBuff.size() - 2; jj >= 0; jj--) {
            dielectron = eBuff.at(jj) + e;
            dielectron_ch = (echBuff.at(jj) + ech) / 2;
            dielectron_weight = eweightBuff.at(jj) * eweight;

            if (dielectron_ch == 0)
              registry.fill(HIST("ULS_orig"), dielectron.M(), dielectron.Pt(), dielectron_weight);
            if (dielectron_ch > 0)
              registry.fill(HIST("LSpp_orig"), dielectron.M(), dielectron.Pt(), dielectron_weight);
            if (dielectron_ch < 0)
              registry.fill(HIST("LSmm_orig"), dielectron.M(), dielectron.Pt(), dielectron_weight);
            if (e.Pt() > fConfigMinPt && eBuff.at(jj).Pt() > fConfigMinPt && e.Pt() < fConfigMaxPt && eBuff.at(jj).Pt() < fConfigMaxPt && TMath::Abs(e.Eta()) < fConfigMaxEta && TMath::Abs(eBuff.at(jj).Eta()) < fConfigMaxEta && e.Vect().Unit().Dot(eBuff.at(jj).Vect().Unit()) < TMath::Cos(fConfigMinOpAng)) {
              if (dielectron_ch == 0)
                registry.fill(HIST("ULS"), dielectron.M(), dielectron.Pt(), dielectron_weight);
              if (dielectron_ch > 0)
                registry.fill(HIST("LSpp"), dielectron.M(), dielectron.Pt(), dielectron_weight);
              if (dielectron_ch < 0)
                registry.fill(HIST("LSmm"), dielectron.M(), dielectron.Pt(), dielectron_weight);
            }
          }
        }

        if (skipNext) {
          skipNext = false;
          continue; // skip if marked as second electron
        }

        if (!(mctrack.getMotherTrackId() > -1))
          continue; // has no mother

        auto const& mother = mctracks[mctrack.getMotherTrackId()];

        if (mother.getMotherTrackId() > -1)
          continue; // mother is not primary

        if (mctrack.getSecondMotherTrackId() - mctrack.getMotherTrackId() > 0)
          continue; // more than one mother

        // skip for the moment other particles rather than pi0, eta, etaprime,
        // omega, rho, phi.
        switch (mother.GetPdgCode()) {
          case 111:
            break;
          case 221:
            break;
          case 331:
            break;
          case 113:
            break;
          case 223:
            break;
          case 333:
            break;
          case 443:
            break;
          default:
            continue;
        }

        // Not sure about this cut. From GammaConv group. Harmless a priori.
        if (!(fabs(mctrack.GetEnergy() - mctrack.Pz()) > 0.))
          continue;

        // ???? this applied only to first daughter!
        Double_t yPre = (mctrack.GetEnergy() + mctrack.Pz()) / (mctrack.GetEnergy() - mctrack.Pz());
        Double_t y = 0.5 * TMath::Log(yPre);
        if (fConfigDoRapidityCut) { // Apply rapidity cut on mother consistent with GammaConv group. (??? but it is not applied on mother?)
          if (yPre <= 0.)
            continue;
          if (TMath::Abs(y) > fConfigRapidityCut)
            continue;
        } else {
          if (yPre == 0.)
            continue;
        }

        treeWords.fdectyp = mother.getLastDaughterTrackId() - mother.getFirstDaughterTrackId() + 1; // fdectyp: decay type (based on number of daughters).
        if (treeWords.fdectyp > 4)
          continue; // exclude five or more particles decay

        if (trackID == mctracks.size())
          continue; // no particle left in the list
        auto mctrack2 = mctracks[trackID + 1];
        if (!(mctrack2.getMotherTrackId() == mctrack.getMotherTrackId()))
          continue; // no matching second electron
        if (!(mctrack.getSecondMotherTrackId() == -1))
          continue; // second daughter has more than one mother
        if (!(abs(mctrack2.GetPdgCode()) == 11))
          continue; // not an electron

        skipNext = true; // is matching electron --> next particle in list will be skipped

        PxPyPzEVector dau1, dau2, ee;
        dau1.SetPxPyPzE(mctrack.Px(), mctrack.Py(), mctrack.Pz(), mctrack.GetEnergy());
        dau2.SetPxPyPzE(mctrack2.Px(), mctrack2.Py(), mctrack2.Pz(), mctrack2.GetEnergy());

        // create dielectron before resolution effects:
        ee = dau1 + dau2;

        // get info of the other particles in the decay:
        treeWords.fdau3pdg = 0;
        for (Int_t jj = mother.getFirstDaughterTrackId(); jj <= mother.getLastDaughterTrackId(); jj++) {
          if (jj == trackID || jj == trackID + 1) {
            continue; // first or second electron
          }
          auto mctrack3 = mctracks[jj];
          treeWords.fdau3pdg = abs(mctrack3.GetPdgCode());
        }

        // get index for histograms
        Int_t hindex[3];
        for (Int_t jj = 0; jj < 3; jj++) {
          hindex[jj] = -1;
        }
        switch (mother.GetPdgCode()) {
          case 111:
            hindex[0] = 0;
            break;
          case 221:
            hindex[0] = 1;
            break;
          case 331:
            hindex[0] = 2;
            if (treeWords.fdectyp == 3 && treeWords.fdau3pdg == 22)
              hindex[1] = 3;
            if (treeWords.fdectyp == 3 && treeWords.fdau3pdg == 223)
              hindex[1] = 4;
            break;
          case 113:
            hindex[0] = 5;
            break;
          case 223:
            hindex[0] = 6;
            if (treeWords.fdectyp == 2)
              hindex[1] = 7;
            if (treeWords.fdectyp == 3 && treeWords.fdau3pdg == 111)
              hindex[1] = 8;
            break;
          case 333:
            hindex[0] = 9;
            if (treeWords.fdectyp == 2)
              hindex[1] = 10;
            if (treeWords.fdectyp == 3 && treeWords.fdau3pdg == 221)
              hindex[1] = 11;
            if (treeWords.fdectyp == 3 && treeWords.fdau3pdg == 111)
              hindex[1] = 12;
            break;
          case 443:
            hindex[0] = 13;
            if (treeWords.fdectyp == 2)
              hindex[1] = 14;
            if (treeWords.fdectyp == 3 && treeWords.fdau3pdg == 22)
              hindex[1] = 15;
            break;
        }

        hindex[2] = nInputParticles;

        if (hindex[0] < 0) {
          LOGP(error, "hindex[0]<0");
          continue;
        }

        // Fill tree words before resolution/acceptance
        treeWords.fd1origpt = dau1.Pt();
        treeWords.fd1origp = dau1.P();
        treeWords.fd1origeta = dau1.Eta();
        treeWords.fd1origphi = dau1.Phi();
        treeWords.fd2origpt = dau2.Pt();
        treeWords.fd2origp = dau2.P();
        treeWords.fd2origeta = dau2.Eta();
        treeWords.fd2origphi = dau2.Phi();
        treeWords.feeorigpt = ee.Pt();
        treeWords.feeorigp = ee.P();
        treeWords.feeorigm = ee.M();
        treeWords.feeorigeta = ee.Eta();
        treeWords.feeorigrap = ee.Rapidity();
        treeWords.feeorigphi = ee.Phi();
        if (mctrack.GetPdgCode() > 0) {
          treeWords.feeorigphiv = PhiV(dau1, dau2);
        } else {
          treeWords.feeorigphiv = PhiV(dau2, dau1);
        }

        // get the efficiency weight
        Int_t effbin = fhwEffpT->FindBin(treeWords.fd1origpt);
        treeWords.fwEffpT = fhwEffpT->GetBinContent(effbin);
        effbin = fhwEffpT->FindBin(treeWords.fd2origpt);
        treeWords.fwEffpT = treeWords.fwEffpT * fhwEffpT->GetBinContent(effbin);

        // Resolution and acceptance
        //-------------------------
        if (mctrack.GetPdgCode() > 0) {
          dau1 = ApplyResolution(dau1, -1, fConfigResolType);
        } else {
          dau1 = ApplyResolution(dau1, 1, fConfigResolType);
        }
        if (mctrack2.GetPdgCode() > 0) {
          dau2 = ApplyResolution(dau2, -1, fConfigResolType);
        } else {
          dau2 = ApplyResolution(dau2, 1, fConfigResolType);
        }
        treeWords.fpass = true;
        if (dau1.Pt() < fConfigMinPt || dau2.Pt() < fConfigMinPt)
          treeWords.fpass = false; // leg pT cut
        if (dau1.Pt() > fConfigMaxPt || dau2.Pt() > fConfigMaxPt)
          treeWords.fpass = false; // leg pT cut
        if (dau1.Vect().Unit().Dot(dau2.Vect().Unit()) > TMath::Cos(fConfigMinOpAng))
          treeWords.fpass = false; // opening angle cut
        if (TMath::Abs(dau1.Eta()) > fConfigMaxEta || TMath::Abs(dau2.Eta()) > fConfigMaxEta)
          treeWords.fpass = false;

        // get the pair DCA (based in smeared pT)
        for (int jj = 0; jj < nbDCAtemplate; jj++) { // loop over DCA templates
          if (dau1.Pt() >= DCATemplateEdges[jj] && dau1.Pt() < DCATemplateEdges[jj + 1]) {
            treeWords.fd1DCA = fh_DCAtemplates[jj]->GetRandom();
          }
          if (dau2.Pt() >= DCATemplateEdges[jj] && dau2.Pt() < DCATemplateEdges[jj + 1]) {
            treeWords.fd2DCA = fh_DCAtemplates[jj]->GetRandom();
          }
        }
        treeWords.fpairDCA = sqrt((pow(treeWords.fd1DCA, 2) + pow(treeWords.fd2DCA, 2)) / 2);

        // Fill tree words after resolution/acceptance
        ee = dau1 + dau2;
        treeWords.fd1pt = dau1.Pt();
        treeWords.fd1p = dau1.P();
        treeWords.fd1eta = dau1.Eta();
        treeWords.fd1phi = dau1.Phi();
        treeWords.fd2pt = dau2.Pt();
        treeWords.fd2p = dau2.P();
        treeWords.fd2eta = dau2.Eta();
        treeWords.fd2phi = dau2.Phi();
        treeWords.feept = ee.Pt();
        treeWords.feemt = ee.Mt();
        treeWords.feep = ee.P();
        treeWords.feem = ee.M();
        treeWords.feeeta = ee.Eta();
        treeWords.feerap = ee.Rapidity();
        treeWords.feephi = ee.Phi();
        if (mctrack.GetPdgCode() > 0) {
          treeWords.feephiv = PhiV(dau1, dau2);
        } else {
          treeWords.feephiv = PhiV(dau2, dau1);
        }
        treeWords.fmotherpt = mother.GetPt();
        treeWords.fmotherm = sqrt(pow(mother.GetEnergy(), 2) + pow(mother.GetP(), 2));
        treeWords.fmothermt = sqrt(pow(treeWords.fmotherm, 2) + pow(treeWords.fmotherpt, 2));
        treeWords.fmotherp = mother.GetP();
        treeWords.fmothereta = mother.GetEta();
        treeWords.fmotherphi = mother.GetPhi();
        treeWords.fID = mother.GetPdgCode();
        treeWords.fweight = mctrack.getWeight(); // get particle weight from generator

        // get multiplicity based weight:
        int iwbin = fhwMultpT->FindBin(treeWords.fmotherpt);
        treeWords.fwMultpT = fhwMultpT->GetBinContent(iwbin);   // pT weight
        treeWords.fwMultpT2 = fhwMultpT2->GetBinContent(iwbin); // pT weight
        double min_mT = fhwMultmT->GetBinLowEdge(1);            // consider as minimum valid mT value the edge of the weight histo.
        if (treeWords.fmothermt > min_mT) {
          iwbin = fhwMultmT->FindBin(treeWords.fmothermt);
          treeWords.fwMultmT = fhwMultmT->GetBinContent(iwbin);   // mT weight
          treeWords.fwMultmT2 = fhwMultmT2->GetBinContent(iwbin); // mT weight
        } else {
          LOGP(error, "Generated particle with mT < Pion mass cannot be weighted");
          treeWords.fwMultmT = 0.;
          treeWords.fwMultmT2 = 0.;
        }

        // Which ALT weight to use?:
        Double_t fwALT = treeWords.fwEffpT; // by default use pt efficiency weight
        if (fConfigALTweight == 1)
          fwALT = treeWords.fwMultmT; // mT multiplicity weight
        if (fConfigALTweight == 11)
          fwALT = treeWords.fwMultmT2; // mT multiplicity weight, higher mult
        if (fConfigALTweight == 2)
          fwALT = treeWords.fwMultpT; // pT multiplicity weight
        if (fConfigALTweight == 22)
          fwALT = treeWords.fwMultpT2; // pT multiplicity weight, higher mult

        // fill the tree
        if (fConfigWriteTTree) {
          tree->Fill();
        }

        // fill the histograms
        if (treeWords.fdectyp < 4) {         // why here <4 and before <5 ???
          for (Int_t jj = 0; jj < 3; jj++) { // fill the different hindex -> particles
            if (hindex[jj] > -1) {
              fmee_orig[hindex[jj]]->Fill(treeWords.feeorigm, treeWords.fweight);
              if (fConfigALTweight == 1 || fConfigALTweight == 11) {
                fmotherpT_orig[hindex[jj]]->Fill(treeWords.fmothermt, treeWords.fweight);
              } else if (fConfigALTweight == 2 || fConfigALTweight == 22 || fConfigALTweight == 0) {
                fmotherpT_orig[hindex[jj]]->Fill(treeWords.fmotherpt, treeWords.fweight);
              }
              fpteevsmee_orig[hindex[jj]]->Fill(treeWords.feeorigm, treeWords.feept, treeWords.fweight);
              fphi_orig[hindex[jj]]->Fill(treeWords.feeorigphi, treeWords.fweight);
              frap_orig[hindex[jj]]->Fill(treeWords.feeorigrap, treeWords.fweight);
              fmee_orig_wALT[hindex[jj]]->Fill(treeWords.feeorigm, treeWords.fweight * fwALT);
              fpteevsmee_orig_wALT[hindex[jj]]->Fill(treeWords.feeorigm, treeWords.feept, treeWords.fweight * fwALT);
              if (fConfigALTweight == 1 || fConfigALTweight == 11) {
                fmotherpT_orig_wALT[hindex[jj]]->Fill(treeWords.fmothermt, treeWords.fweight * fwALT);
              } else if (fConfigALTweight == 2 || fConfigALTweight == 22 || fConfigALTweight == 0) {
                fmotherpT_orig_wALT[hindex[jj]]->Fill(treeWords.fmotherpt, treeWords.fweight * fwALT);
              }
              if (treeWords.fpass) {
                fmee[hindex[jj]]->Fill(treeWords.feem, treeWords.fweight);
                fpteevsmee[hindex[jj]]->Fill(treeWords.feem, treeWords.feept, treeWords.fweight);
                fphi[hindex[jj]]->Fill(treeWords.feephi, treeWords.fweight);
                frap[hindex[jj]]->Fill(treeWords.feerap, treeWords.fweight);
                registry.fill(HIST("DCAeevsmee"), treeWords.feem, treeWords.fpairDCA, treeWords.fweight);
                registry.fill(HIST("DCAeevsptee"), treeWords.feept, treeWords.fpairDCA, treeWords.fweight);
                fmee_wALT[hindex[jj]]->Fill(treeWords.feem, treeWords.fweight * fwALT);
                fpteevsmee_wALT[hindex[jj]]->Fill(treeWords.feem, treeWords.feept, treeWords.fweight * fwALT);
              }
            }
          }
        }

        if (fConfigDoVirtPh) {
          // Virtual photon generation
          //-------------------------
          // We will generate one virtual photon per histogrammed pion
          if (mother.GetPdgCode() == 111) {
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
            if (SetDecay == 0)
              LOGP(error, "Decay not permitted by kinematics");
            Double_t VPHweight = VPHgen.Generate();
            // get electrons from the decay
            TLorentzVector *decay1, *decay2;
            decay1 = VPHgen.GetDecay(0);
            decay2 = VPHgen.GetDecay(1);
            dau1.SetPxPyPzE(decay1->Px(), decay1->Py(), decay1->Pz(), decay1->E());
            dau2.SetPxPyPzE(decay2->Px(), decay2->Py(), decay2->Pz(), decay2->E());

            // create dielectron before resolution effects:
            ee = dau1 + dau2;

            // get index for histograms
            hindex[0] = nInputParticles - 1;
            hindex[1] = -1;
            hindex[2] = -1;

            // Fill tree words before resolution/acceptance
            treeWords.fd1origpt = dau1.Pt();
            treeWords.fd1origp = dau1.P();
            treeWords.fd1origeta = dau1.Eta();
            treeWords.fd1origphi = dau1.Phi();
            treeWords.fd2origpt = dau2.Pt();
            treeWords.fd2origp = dau2.P();
            treeWords.fd2origeta = dau2.Eta();
            treeWords.fd2origphi = dau2.Phi();
            treeWords.feeorigpt = ee.Pt();
            treeWords.feeorigp = ee.P();
            treeWords.feeorigm = ee.M();
            treeWords.feeorigeta = ee.Eta();
            treeWords.feeorigrap = ee.Rapidity();
            treeWords.feeorigphi = ee.Phi();
            treeWords.feeorigphiv = PhiV(dau1, dau2);

            // get the efficiency weight
            Int_t effbin = fhwEffpT->FindBin(treeWords.fd1origpt);
            treeWords.fwEffpT = fhwEffpT->GetBinContent(effbin);
            effbin = fhwEffpT->FindBin(treeWords.fd2origpt);
            treeWords.fwEffpT = treeWords.fwEffpT * fhwEffpT->GetBinContent(effbin);

            // Resolution and acceptance
            //-------------------------
            dau1 = ApplyResolution(dau1, 1, fConfigResolType);
            dau2 = ApplyResolution(dau2, -1, fConfigResolType);
            treeWords.fpass = true;
            if (dau1.Pt() < fConfigMinPt || dau2.Pt() < fConfigMinPt)
              treeWords.fpass = false; // leg pT cut
            if (dau1.Pt() > fConfigMaxPt || dau2.Pt() > fConfigMaxPt)
              treeWords.fpass = false; // leg pT cut
            if (dau1.Vect().Unit().Dot(dau2.Vect().Unit()) > TMath::Cos(fConfigMinOpAng))
              treeWords.fpass = false; // opening angle cut
            if (TMath::Abs(dau1.Eta()) > fConfigMaxEta || TMath::Abs(dau2.Eta()) > fConfigMaxEta)
              treeWords.fpass = false;

            treeWords.fpairDCA = 10000.; // ??

            // Fill tree words after resolution/acceptance
            ee = dau1 + dau2;
            treeWords.fd1pt = dau1.Pt();
            treeWords.fd1p = dau1.P();
            treeWords.fd1eta = dau1.Eta();
            treeWords.fd1phi = dau1.Phi();
            treeWords.fd2pt = dau2.Pt();
            treeWords.fd2p = dau2.P();
            treeWords.fd2eta = dau2.Eta();
            treeWords.fd2phi = dau2.Phi();
            treeWords.feept = ee.Pt();
            treeWords.feemt = ee.Mt();
            treeWords.feep = ee.P();
            treeWords.feem = ee.M();
            treeWords.feeeta = ee.Eta();
            treeWords.feerap = ee.Rapidity();
            treeWords.feephi = ee.Phi();
            treeWords.feephiv = PhiV(dau1, dau2);
            treeWords.fmotherpt = beam.Pt();
            treeWords.fmothermt = sqrt(pow(beam.M(), 2) + pow(beam.Pt(), 2));
            treeWords.fmotherp = beam.P();
            treeWords.fmotherm = beam.M();
            treeWords.fmothereta = beam.Eta();
            treeWords.fmotherphi = beam.Phi();
            treeWords.fID = 0; // set ID to Zero for VPH
            treeWords.fweight = VPHweight;
            // get multiplicity based weight:
            treeWords.fwMultmT = 1; // no weight for photons so far

            // Fill the tree
            if (fConfigWriteTTree) { // many parameters not set for photons: d1DCA,fd2DCA, fdectyp,fdau3pdg,fwMultpT,fwMultpT2,fwMultmT2
              tree->Fill();
            }

            // Fill the histograms
            for (Int_t jj = 0; jj < 3; jj++) { // fill the different hindex -> particles
              if (hindex[jj] > -1) {
                fmee_orig[hindex[jj]]->Fill(treeWords.feeorigm, VPHweight);
                fpteevsmee_orig[hindex[jj]]->Fill(treeWords.feeorigm, treeWords.feept, VPHweight);
                fphi_orig[hindex[jj]]->Fill(treeWords.feeorigphi, VPHweight);
                frap_orig[hindex[jj]]->Fill(treeWords.feeorigrap, VPHweight);
                fmotherpT_orig[hindex[jj]]->Fill(treeWords.fmotherpt, treeWords.fweight);
                if (treeWords.fpass) {
                  fmee[hindex[jj]]->Fill(treeWords.feem, VPHweight);
                  fpteevsmee[hindex[jj]]->Fill(treeWords.feem, treeWords.feept, VPHweight);
                  fphi[hindex[jj]]->Fill(treeWords.feephi, VPHweight);
                  frap[hindex[jj]]->Fill(treeWords.feerap, VPHweight);
                }
              }
            }

          } // mother.pdgCode()==111
        }   // fConfigDoVirtPh

      } // abs(pdgCode())==11

    } // loop over mctracks

    // Clear buffers
    eBuff.clear();
    echBuff.clear();
    eweightBuff.clear();
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

    if (fConfigDoPairing) {
      registry.add<TH2>("ULS", "ULS", HistType::kTH2F, {mAxis, ptAxis}, true);
      registry.add<TH2>("LSpp", "LSpp", HistType::kTH2F, {mAxis, ptAxis}, true);
      registry.add<TH2>("LSmm", "LSmm", HistType::kTH2F, {mAxis, ptAxis}, true);

      registry.add<TH2>("ULS_orig", "ULS_orig", HistType::kTH2F, {mAxis, ptAxis}, true);
      registry.add<TH2>("LSpp_orig", "LSpp_orig", HistType::kTH2F, {mAxis, ptAxis}, true);
      registry.add<TH2>("LSmm_orig", "LSmm_orig", HistType::kTH2F, {mAxis, ptAxis}, true);
    }

    // configurable axes crashed the task. Take them out for the moment
    // registry.add<TH2>("DCAeevsmee", "DCAeevsmee", HistType::kTH2F, {{fConfigMBins, "#it{m}_{ee} (GeV/c^{2})"}, {fConfigDCABins, "DCA_{xy}^{ee} (cm)"}}, true);
    // registry.add<TH2>("DCAeevsptee", "DCAeevsptee", HistType::kTH2F, {{fConfigPtBins, "#it{p}_{T,ee} (GeV/c)"}, {fConfigDCABins, "DCA_{xy}^{ee} (cm)"}}, true);
    // replace them with hard coded axes
    AxisSpec configPtBins = {{0., 0.5, 1, 1.5, 2., 2.5, 3., 3.5, 4., 4.5, 5., 5.5, 6., 6.5, 7., 7.5, 8.}, "#it{p}_{T,ee} (GeV/c)"};
    AxisSpec configMBins = {{0., 0.08, 0.14, 0.2, 1.1, 2.7, 2.8, 3.2, 5.0}, "#it{m}_{ee} (GeV/c^{2})"};
    AxisSpec configDcaBins = {{0., 0.4, 0.8, 1.2, 1.6, 2.0, 2.4, 3., 4., 5., 7., 10.}, "DCA_{xy}^{ee} (cm)"};
    registry.add<TH2>("DCAeevsmee", "DCAeevsmee", HistType::kTH2F, {configMBins, configDcaBins}, true);
    registry.add<TH2>("DCAeevsptee", "DCAeevsptee", HistType::kTH2F, {configPtBins, configDcaBins}, true);

    for (auto& particle : fParticleListNames) {
      fmee.push_back(registry.add<TH1>(Form("mee_%s", particle.Data()), Form("mee_%s", particle.Data()), HistType::kTH1F, {mAxis}, true));
      fmee_orig.push_back(registry.add<TH1>(Form("mee_orig_%s", particle.Data()), Form("mee_orig_%s", particle.Data()), HistType::kTH1F, {mAxis}, true));
      fmotherpT_orig.push_back(registry.add<TH1>(Form("motherpT_orig_%s", particle.Data()), Form("motherpT_orig_%s", particle.Data()), HistType::kTH1F, {ptAxis}, true));
      fphi.push_back(registry.add<TH1>(Form("phi_%s", particle.Data()), Form("phi_%s", particle.Data()), HistType::kTH1F, {phiAxis}, true));
      fphi_orig.push_back(registry.add<TH1>(Form("phi_orig_%s", particle.Data()), Form("phi_orig_%s", particle.Data()), HistType::kTH1F, {phiAxis}, true));
      frap.push_back(registry.add<TH1>(Form("rap_%s", particle.Data()), Form("rap_%s", particle.Data()), HistType::kTH1F, {rapAxis}, true));
      frap_orig.push_back(registry.add<TH1>(Form("rap_orig_%s", particle.Data()), Form("rap_orig_%s", particle.Data()), HistType::kTH1F, {rapAxis}, true));
      fpteevsmee.push_back(registry.add<TH2>(Form("pteevsmee_%s", particle.Data()), Form("pteevsmee_%s", particle.Data()), HistType::kTH2F, {mAxis, ptAxis}, true));
      fpteevsmee_orig.push_back(registry.add<TH2>(Form("pteevsmee_orig_%s", particle.Data()), Form("pteevsmee_orig_%s", particle.Data()), HistType::kTH2F, {mAxis, ptAxis}, true));
      fmee_wALT.push_back(registry.add<TH1>(Form("mee_wALT_%s", particle.Data()), Form("mee_wALT_%s", particle.Data()), HistType::kTH1F, {mAxis}, true));
      fmee_orig_wALT.push_back(registry.add<TH1>(Form("mee_orig_wALT_%s", particle.Data()), Form("mee_orig_wALT_%s", particle.Data()), HistType::kTH1F, {mAxis}, true));
      fmotherpT_orig_wALT.push_back(registry.add<TH1>(Form("motherpT_orig_wALT_%s", particle.Data()), Form("motherpT_orig_wALT%s", particle.Data()), HistType::kTH1F, {ptAxis}, true));
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

  void SetTree()
  {
    tree.setObject(new TTree("eeTTree", "eeTTree"));

    tree->Branch("fd1DCA", &treeWords.fd1DCA, "fd1DCA/F");
    tree->Branch("fd2DCA", &treeWords.fd2DCA, "fd2DCA/F");
    tree->Branch("fpairDCA", &treeWords.fpairDCA, "fpairDCA/F");
    tree->Branch("fd1origpt", &treeWords.fd1origpt, "fd1origpt/F");
    tree->Branch("fd1origp", &treeWords.fd1origp, "fd1origp/F");
    tree->Branch("fd1origeta", &treeWords.fd1origeta, "fd1origeta/F");
    tree->Branch("fd1origphi", &treeWords.fd1origphi, "fd1origphi/F");
    tree->Branch("fd2origpt", &treeWords.fd2origpt, "fd2origpt/F");
    tree->Branch("fd2origp", &treeWords.fd2origp, "fd2origp/F");
    tree->Branch("fd2origeta", &treeWords.fd2origeta, "fd2origeta/F");
    tree->Branch("fd2origphi", &treeWords.fd2origphi, "fd2origphi/F");
    tree->Branch("fd1pt", &treeWords.fd1pt, "fd1pt/F");
    tree->Branch("fd1p", &treeWords.fd1p, "fd1p/F");
    tree->Branch("fd1eta", &treeWords.fd1eta, "fd1eta/F");
    tree->Branch("fd1phi", &treeWords.fd1phi, "fd1phi/F");
    tree->Branch("fd2pt", &treeWords.fd2pt, "fd2pt/F");
    tree->Branch("fd2p", &treeWords.fd2p, "fd2p/F");
    tree->Branch("fd2eta", &treeWords.fd2eta, "fd2eta/F");
    tree->Branch("fd2phi", &treeWords.fd2phi, "fd2phi/F");
    tree->Branch("feeorigpt", &treeWords.feeorigpt, "feeorigpt/F");
    tree->Branch("feeorigp", &treeWords.feeorigp, "feeorigp/F");
    tree->Branch("feeorigm", &treeWords.feeorigm, "feeorigm/F");
    tree->Branch("feeorigeta", &treeWords.feeorigeta, "feeorigeta/F");
    tree->Branch("feeorigphi", &treeWords.feeorigphi, "feeorigphi/F");
    tree->Branch("feeorigphiv", &treeWords.feeorigphiv, "feeorigphiv/F");
    tree->Branch("feept", &treeWords.feept, "feept/F");
    tree->Branch("feemt", &treeWords.feemt, "feemt/F");
    tree->Branch("feep", &treeWords.feep, "feep/F");
    tree->Branch("feem", &treeWords.feem, "feem/F");
    tree->Branch("feeeta", &treeWords.feeeta, "feeeta/F");
    tree->Branch("feephi", &treeWords.feephi, "feephi/F");
    tree->Branch("feephiv", &treeWords.feephiv, "feephiv/F");
    tree->Branch("fmotherpt", &treeWords.fmotherpt, "fmotherpt/F");
    tree->Branch("fmothermt", &treeWords.fmothermt, "fmothermt/F");
    tree->Branch("fmotherp", &treeWords.fmotherp, "fmotherp/F");
    tree->Branch("fmotherm", &treeWords.fmotherm, "fmotherm/F");
    tree->Branch("fmothereta", &treeWords.fmothereta, "fmothereta/F");
    tree->Branch("fmotherphi", &treeWords.fmotherphi, "fmotherphi/F");
    tree->Branch("fID", &treeWords.fID, "fID/I");
    tree->Branch("fdectyp", &treeWords.fdectyp, "fdectyp/I");
    tree->Branch("fdau3pdg", &treeWords.fdau3pdg, "fdau3pdg/I");
    tree->Branch("fweight", &treeWords.fweight, "fweight/F");
    tree->Branch("fwEffpT", &treeWords.fwEffpT, "fwEffpT/F");
    tree->Branch("fwMultpT", &treeWords.fwMultpT, "fwMultpT/F");
    tree->Branch("fwMultmT", &treeWords.fwMultmT, "fwMultmT/F");
    tree->Branch("fwMultpT2", &treeWords.fwMultpT2, "fwMultpT2/F");
    tree->Branch("fwMultmT2", &treeWords.fwMultmT2, "fwMultmT2/F");
    tree->Branch("fpass", &treeWords.fpass, "fpass/B");
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
      Double_t smearing = 0.;
      TH1D* thisHist = reinterpret_cast<TH1D*>(fArrResoPt->At(ptbin));
      if (thisHist->GetEntries() > 0) {
        smearing = thisHist->GetRandom() * pt;
      }
      Double_t sPt = pt - smearing;

      // smear eta
      ptbin = reinterpret_cast<TH2D*>(fArrResoEta->At(0))->GetXaxis()->FindBin(pt);
      if (ptbin < 1)
        ptbin = 1;
      if (ptbin > fArrResoEta->GetLast())
        ptbin = fArrResoEta->GetLast();
      smearing = 0.;
      thisHist = reinterpret_cast<TH1D*>(fArrResoEta->At(ptbin));
      if (thisHist->GetEntries() > 0) {
        smearing = thisHist->GetRandom();
      }
      Double_t sEta = eta - smearing;

      // smear phi
      ptbin = reinterpret_cast<TH2D*>(fArrResoPhi_Pos->At(0))->GetXaxis()->FindBin(pt);
      if (ptbin < 1)
        ptbin = 1;
      if (ptbin > fArrResoPhi_Pos->GetLast())
        ptbin = fArrResoPhi_Pos->GetLast();
      smearing = 0.;
      if (ch > 0) {
        thisHist = reinterpret_cast<TH1D*>(fArrResoPhi_Pos->At(ptbin));
      } else if (ch < 0) {
        thisHist = reinterpret_cast<TH1D*>(fArrResoPhi_Neg->At(ptbin));
      }
      if (thisHist->GetEntries() > 0) {
        smearing = thisHist->GetRandom();
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
    LOGP(info, "Set Efficiency histo");
    // Get Efficiency weight file:
    TFile* fFile = TFile::Open(filename.Data());
    if (!fFile) {
      LOGP(error, "Could not open Efficiency file {}", filename.Data());
      return;
    }
    if (fFile->GetListOfKeys()->Contains(histname.Data())) {
      fhwEffpT = reinterpret_cast<TH1F*>(fFile->Get(histname.Data())); // histo: eff weight in function of pT.
      fhwEffpT->SetDirectory(nullptr);
    } else {
      LOGP(error, "Could not open histogram {} from file {}", histname.Data(), filename.Data());
    }

    fFile->Close();
  }

  void GetResHisto(TString filename, TString ptHistName, TString etaHistName, TString phiPosHistName, TString phiNegHistName)
  {
    // get resolutoin histo
    LOGP(info, "Set Resolution histo");
    // Get Resolution map
    TFile* fFile = TFile::Open(filename.Data());
    if (!fFile) {
      LOGP(error, "Could not open Resolution file {}", filename.Data());
      return;
    }
    TObjArray* ArrResoPt = nullptr;
    if (fFile->GetListOfKeys()->Contains(ptHistName.Data())) {
      ArrResoPt = reinterpret_cast<TObjArray*>(fFile->Get(ptHistName.Data()));
    } else {
      LOGP(error, "Could not open {} from file {}", ptHistName.Data(), filename.Data());
    }

    TObjArray* ArrResoEta = nullptr;
    if (fFile->GetListOfKeys()->Contains(etaHistName.Data())) {
      ArrResoEta = reinterpret_cast<TObjArray*>(fFile->Get(etaHistName.Data()));
    } else {
      LOGP(error, "Could not open {} from file {}", etaHistName.Data(), filename.Data());
    }

    TObjArray* ArrResoPhi_Pos = nullptr;
    if (fFile->GetListOfKeys()->Contains(phiPosHistName.Data())) {
      ArrResoPhi_Pos = reinterpret_cast<TObjArray*>(fFile->Get(phiPosHistName.Data()));
    } else {
      LOGP(error, "Could not open {} from file {}", phiPosHistName.Data(), filename.Data());
    }

    TObjArray* ArrResoPhi_Neg = nullptr;
    if (fFile->GetListOfKeys()->Contains(phiNegHistName.Data())) {
      ArrResoPhi_Neg = reinterpret_cast<TObjArray*>(fFile->Get(phiNegHistName.Data()));
    } else {
      LOGP(error, "Could not open {} from file {}", phiNegHistName.Data(), filename.Data());
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
    LOGP(info, "Set Resolution histo");
    // Get Resolution map
    TFile* fFile = TFile::Open(filename.Data());
    if (!fFile) {
      LOGP(error, "Could not open Resolution file {}", filename.Data());
      LOGP(error, "No resolution array set! Using internal parametrization");
      return;
    }
    TObjArray* ArrResoP = 0x0;
    if (fFile->GetListOfKeys()->Contains(pHistName.Data())) {
      ArrResoP = reinterpret_cast<TObjArray*>(fFile->Get(pHistName.Data()));
    } else {
      LOGP(error, "Could not open {} from file {}", pHistName.Data(), filename.Data());
      LOGP(error, "No resolution array set! Using internal parametrization");
    }

    fArr = ArrResoP;
    fFile->Close();
  }

  void GetDCATemplates(TString filename, TString histname)
  {
    // get dca tamplates
    LOGP(info, "Set DCA templates");
    // Get  file:
    TFile* fFile = TFile::Open(filename.Data());
    if (!fFile) {
      LOGP(error, "Could not open DCATemplate file {}", filename.Data());
      return;
    }
    fh_DCAtemplates = new TH1F*[nbDCAtemplate];
    for (int jj = 0; jj < nbDCAtemplate; jj++) {
      if (fFile->GetListOfKeys()->Contains(Form("%s%d", histname.Data(), jj + 1))) {
        fh_DCAtemplates[jj] = reinterpret_cast<TH1F*>(fFile->Get(Form("%s%d", histname.Data(), jj + 1)));
      } else {
        LOGP(error, "Could not open {}{} from file {}", histname.Data(), jj + 1, filename.Data());
      }
    }
    fFile->Close();
  }

  void GetMultHisto(TString filename, TString histnamept, TString histnamept2, TString histnamemt, TString histnamemt2)
  {
    // get multiplicity weights
    LOGP(info, "Set Multiplicity weight files");
    TFile* fFile = TFile::Open(filename.Data());
    if (!fFile) {
      LOGP(error, "Could not open Multiplicity weight file {}", filename.Data());
      return;
    }

    if (fFile->GetListOfKeys()->Contains(histnamept.Data())) {
      fhwMultpT = reinterpret_cast<TH1F*>(fFile->Get(histnamept.Data())); // histo: multiplicity weight in function of pT.
    } else {
      LOGP(error, "Could not open {} from file {}", histnamept.Data(), filename.Data());
    }

    if (fFile->GetListOfKeys()->Contains(histnamemt.Data())) {
      fhwMultmT = reinterpret_cast<TH1F*>(fFile->Get(histnamemt.Data())); // histo: multiplicity weight in function of mT.
    } else {
      LOGP(error, "Could not open {} from file {}", histnamemt.Data(), filename.Data());
    }

    if (fFile->GetListOfKeys()->Contains(histnamept2.Data())) {
      fhwMultpT2 = reinterpret_cast<TH1F*>(fFile->Get(histnamept2.Data())); // histo: multiplicity weight in function of pT.
    } else {
      LOGP(error, "Could not open {} from file {}", histnamept2.Data(), filename.Data());
    }

    if (fFile->GetListOfKeys()->Contains(histnamemt2.Data())) {
      fhwMultmT2 = reinterpret_cast<TH1F*>(fFile->Get(histnamemt2)); // histo: multiplicity weight in function of mT.
    } else {
      LOGP(error, "Could not open {} from file {}", histnamemt2.Data(), filename.Data());
    }
    fFile->Close();
  }

  void GetPhotonPtParametrization(TString filename, TString dirname, TString funcname)
  {
    LOGP(info, "Set photon parametrization");

    if (filename.EndsWith(".root")) { // read from ROOT file
      TFile* fFile = TFile::Open(filename.Data());
      if (!fFile) {
        LOGP(error, "Could not open photon parametrization from file {}", filename.Data());
        return;
      }
      bool good = false;
      if (fFile->GetListOfKeys()->Contains(dirname.Data())) {
        TDirectory* dir = fFile->GetDirectory(dirname.Data());
        if (dir->GetListOfKeys()->Contains(funcname.Data())) {
          ffVPHpT = reinterpret_cast<TF1*>(dir->Get(funcname.Data()));
          ffVPHpT->SetNpx(10000);
          good = true;
        }
      }
      if (!good) {
        LOGP(error, "Could not open photon parametrization {}/{} from file {}", dirname.Data(), funcname.Data(), filename.Data());
      }
      fFile->Close();
    } else if (filename.EndsWith(".json")) { // read from JSON file
      std::ifstream fFile(filename.Data());
      if (!fFile) {
        LOGP(error, "Could not open photon parametrization from file {}", filename.Data());
        return;
      }
      nlohmann::json paramfile = nlohmann::json::parse(fFile);
      if (paramfile.contains(dirname.Data())) {
        nlohmann::json dir = paramfile[dirname.Data()];
        if (dir.contains(funcname.Data())) {
          std::string formula = dir[funcname.Data()];
          ffVPHpT = new TF1(TString(funcname.Data()), TString(formula), 0, 100);
          if (ffVPHpT) {
            ffVPHpT->SetNpx(10000);
            return;
          }
        }
      }
      LOGP(error, "Could not open photon parametrization {}/{} from file {}", dirname.Data(), funcname.Data(), filename.Data());
      return;
    } else { // neither ROOT nor JSON
      LOGP(error, "Not compatible file format for {}", filename.Data());
    }
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
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec specs;
  std::vector<InputSpec> inputs;
  inputs.emplace_back("mctracks", "MC", "MCTRACKS", 0., Lifetime::Timeframe);
  DataProcessorSpec dSpec = adaptAnalysisTask<lmeelfcocktail>(cfgc, TaskName{"em-lmee-lf-cocktail"});
  dSpec.inputs = inputs;
  specs.emplace_back(dSpec);
  return specs;
}
