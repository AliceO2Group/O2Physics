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

/// \file
/// \brief macro for the calculation of FONLL+PYTHIA predictions for non-prompt charm hadrons
///
/// \author Fabrizio Grosa <fgrosa@cern.ch>, CERN

#if !defined(__CINT__) || defined(__CLING__)

#include <CommonConstants/MathConstants.h>
#include <Framework/Logger.h>

#include <TDatabasePDG.h>
#include <TDirectoryFile.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1.h>
#include <TROOT.h>
#include <TRandom3.h>
#include <TString.h>

#include <Pythia8/Pythia.h>

#include <array>
#include <cmath>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#endif

enum BHadrons {
  Bplus = 0,
  Bzero,
  Bs,
  Lb,
  NBeautyHadrons
};

enum CHadrons {
  Dplus = 0,
  Dzero,
  Ds,
  Lc,
  DstarPlus,
  NCharmHadrons
};

enum BeautyFragFracOptions {
  EpEm = 0,
  PPbar,
  LHCb,
  LHCbMin,
  LHCbMax,
  NBeutyFragFracOptions
};

const std::array<float, NBeautyHadrons> beautyFragFracEpEm = {0.408, 0.408, 0.100, 0.084};
const std::array<float, NBeautyHadrons> beautyFragFracPPbar = {0.344, 0.344, 0.115, 0.198};
const std::array<float, NCharmHadrons> charmFragFracEpEm = {0.542, 0.225, 0.092, 0.057, 0.236}; // Values from e+e- ARXIV:1404.3888 (D0, D+, Ds, Lc, D*+)
const std::array<int, NBeautyHadrons> beautyHadPdgs = {511, 521, 531, 5122};
const std::array<int, NCharmHadrons> charmHadPdgs = {411, 421, 431, 4122, 413};
const std::array<std::string, NBeautyHadrons> beautyHadNames = {"Bzero", "Bplus", "Bs", "Lb"};
const std::array<std::string, NCharmHadrons> charmHadNames = {"Dplus", "Dzero", "Ds", "Lc", "DstarPlus"};
std::array<std::string, 3> namesFonll = {"Central", "Min", "Max"};

// FUNCTION PROTOTYPES
//__________________________________________________________________________________________________
void computeFonllPlusPythiaPredictions(int nDecays = 10000000,
                                       int seed = 42,
                                       std::string inFileFonllBhad = "fonll_bhadron_5dot5teV_y1.txt",
                                       int fragFracOpt = EpEm,
                                       bool addPromptCharmHadrons = true,
                                       std::string inFileFonllPromptDzero = "fonll_prompt_dzero_5dot5teV_y05.txt",
                                       std::string inFileFonllPromptDplus = "fonll_prompt_dplus_5dot5teV_y05.txt",
                                       std::string inFileFonllPromptDstarPlus = "fonll_prompt_dstar_5dot5teV_y05.txt",
                                       std::string outFileName = "fonll_pythia_beautyFFee_charmhadrons_5dot5tev_y0dot5.root");
std::vector<std::string> splitString(const std::string& str, char delimiter);
std::array<TH1D*, 3> readFonll(std::string inFile, std::string histName = "hFonllBhadron");

// FUNCTION IMPLEMENTATIONS
//__________________________________________________________________________________________________
std::vector<std::string> splitString(const std::string& str, char delimiter)
{
  std::vector<std::string> tokens;
  std::stringstream ss(str);
  std::string token;

  while (std::getline(ss, token, delimiter)) {
    tokens.push_back(token);
  }
  tokens.erase(std::remove_if(tokens.begin(), tokens.end(), [](const std::string& str) {
                 return str.find_first_not_of(' ') == std::string::npos; // Check if the string contains only spaces
               }),
               tokens.end());

  return tokens;
}

//__________________________________________________________________________________________________
std::array<TH1D*, 3> readFonll(const std::string& inFile, const std::string& histName)
{

  std::array<TH1D*, 3> hFonll{nullptr, nullptr, nullptr};

  std::ifstream inputFile(inFile);

  if (!inputFile) {
    LOGP(fatal, "Error opening file {}", inFile);
    return hFonll;
  }

  std::string line;

  std::vector<float> ptCent{};
  std::vector<float> crossSecCent{};
  std::vector<float> crossSecMin{};
  std::vector<float> crossSecMax{};
  while (std::getline(inputFile, line)) {
    if (line.find("#") != std::string::npos) {
      continue;
    }
    std::vector<std::string> elements = splitString(line, ' ');
    ptCent.push_back(std::stof(elements[0]));
    crossSecCent.push_back(std::stof(elements[1]));
    crossSecMin.push_back(std::stof(elements[2]));
    crossSecMax.push_back(std::stof(elements[3]));
  }
  inputFile.close();

  if (ptCent.size() < 2) {
    LOGP(fatal, "Only one pT value in FONLL file {}, cannot deduce binning.", inFile);
  }

  float ptWidth = ptCent[1] - ptCent[0];
  float ptMin = ptCent.front() - ptWidth / 2;
  float ptMax = ptCent.back() + ptWidth / 2;

  hFonll[0] = new TH1D(Form("%sCentral", histName.data()), ";#it{p}_{T} (GeV/#it{c});d#sigma/d#it{p}_{T} (#it{c}/GeV)", ptCent.size(), ptMin, ptMax);
  hFonll[1] = new TH1D(Form("%sMin", histName.data()), ";#it{p}_{T} (GeV/#it{c});d#sigma/d#it{p}_{T} (#it{c}/GeV)", ptCent.size(), ptMin, ptMax);
  hFonll[2] = new TH1D(Form("%sMax", histName.data()), ";#it{p}_{T} (GeV/#it{c});d#sigma/d#it{p}_{T} (#it{c}/GeV)", ptCent.size(), ptMin, ptMax);
  for (auto iPt{0u}; iPt < ptCent.size(); ++iPt) {
    hFonll[0]->SetBinContent(iPt + 1, crossSecCent[iPt]);
    hFonll[1]->SetBinContent(iPt + 1, crossSecMin[iPt]);
    hFonll[2]->SetBinContent(iPt + 1, crossSecMax[iPt]);
  }

  return hFonll;
}

//__________________________________________________________________________________________________
void computeFonllPlusPythiaPredictions(int nDecays, int seed, std::string inFileFonllBhad, int fragFracOpt, bool addPromptCharmHadrons, std::string inFileFonllPromptDzero, std::string inFileFonllPromptDplus, std::string inFileFonllPromptDstarPlus, std::string outFileName)
{

  gROOT->SetBatch(true);
  gRandom->SetSeed(seed);

  if (fragFracOpt >= NBeutyFragFracOptions) {
    LOGP(fatal, "Fragmentation fraction option not supported! Exit");
  }

  // init pythia object for the decayer
  Pythia8::Pythia pythia;
  pythia.readString("SoftQCD:inelastic = on");
  pythia.readString("Tune:pp = 14");
  pythia.readString("Random:setSeed = on");
  pythia.readString(Form("Random:seed %d", seed));
  pythia.init();

  // get histograms from FONLL
  auto hFonllBhad = readFonll(inFileFonllBhad);
  if (hFonllBhad[0] == nullptr) {
    return;
  }

  std::map<int, std::array<TH1D*, 3>> hFonllPromptChad{};
  if (addPromptCharmHadrons) {
    hFonllPromptChad[411] = readFonll(inFileFonllPromptDplus, "hFonllPromptDplus");
    hFonllPromptChad[421] = readFonll(inFileFonllPromptDzero, "hFonllPromptDzero");
    hFonllPromptChad[413] = readFonll(inFileFonllPromptDstarPlus, "hFonllPromptDstarPlus");
    // TODO: cook something for Ds and Lc
    hFonllPromptChad[431] = {nullptr, nullptr, nullptr};
    hFonllPromptChad[4122] = {nullptr, nullptr, nullptr};
    for (auto iChad{0}; iChad < NCharmHadrons; ++iChad) {
      if (charmHadPdgs[iChad] == 431 || charmHadPdgs[iChad] == 4122) {
        continue;
      }
      for (auto iFonll{0}; iFonll < 3; ++iFonll) {
        hFonllPromptChad[charmHadPdgs[iChad]][iFonll]->Scale(charmFragFracEpEm[iChad]);
      }
    }
  }

  // initialise histograms for non-prompt charm hadrons
  std::map<int, std::array<std::array<TH1D*, 3>, NBeautyHadrons + 1>> hFonllPythiaNonPromptChad{};
  for (auto iChad{0}; iChad < NCharmHadrons; ++iChad) {
    for (auto iBHad{0}; iBHad < NBeautyHadrons; ++iBHad) {
      for (auto iFonll{0}; iFonll < 3; ++iFonll) {
        hFonllPythiaNonPromptChad[charmHadPdgs[iChad]][iBHad][iFonll] = new TH1D(
          Form("hFonll%sFrom%s%s", charmHadNames[iChad].data(), beautyHadNames[iBHad].data(), namesFonll[iFonll].data()),
          ";#it{p}_{T} (GeV/#it{c});d#sigma/d#it{p}_{T} (#it{c}/GeV)", 1000, 0., 100.);
        hFonllPythiaNonPromptChad[charmHadPdgs[iChad]][iBHad][iFonll]->Sumw2();
      }
    }
  }

  // compute fractions for normalisation
  std::array<float, NBeautyHadrons> fragFracs{};
  std::array<std::array<TF1*, 15>, NBeautyHadrons> fragFracFuncs{}; // used in case of LHCb due to the pT dependence
  if (fragFracOpt == EpEm) {
    fragFracs = beautyFragFracEpEm;
  } else if (fragFracOpt == PPbar) {
    fragFracs = beautyFragFracPPbar;
  } else { // LHCb
    for (int iPar{0}; iPar < 15; ++iPar) {
      fragFracFuncs[0][iPar] = new TF1(Form("fracBz_%d", iPar), "1 /  (2 * (([0] * ([1] + [2] * (x - [3])))  + ([4] * ([5] + exp([6] + [7] * x))) + 1) )", 0.f, 300.f);                           // B0
      fragFracFuncs[1][iPar] = new TF1(Form("fracBp_%d", iPar), "1 /  (2 * (([0] * ([1] + [2] * (x - [3])))  + ([4] * ([5] + exp([6] + [7] * x))) + 1) )", 0.f, 300.f);                           // B+
      fragFracFuncs[2][iPar] = new TF1(Form("fracBs_%d", iPar), "([0] * ([1] + [2] * (x - [3]))) /  (([0] * [1] + [2] * (x - [3]))  + ([4] * ([5] + exp([6] + [7] * x))) + 1)", 0.f, 300.f);      // Bs0
      fragFracFuncs[3][iPar] = new TF1(Form("fracLb_%d", iPar), "([4] * ([5] + exp([6] + [7] * x))) /  (([0] * ([1] + [2] * (x - [3])))  + ([4] * ([5] + exp([6] + [7] * x))) + 1)", 0.f, 300.f); // Lb

      // parameters from https://arxiv.org/pdf/1902.06794.pdf
      float sign = (iPar < 8) ? 1.f : -1.f;
      float parLbA = 1.f + ((iPar == 1 || iPar == 8) ? sign * 0.061f : 0.f);
      float parLbp1 = 0.0793f + ((iPar == 2 || iPar == 9) ? sign * 0.0141f : 0.f);
      float parLbp2 = -1.022f + ((iPar == 3 || iPar == 10) ? sign * 0.0047f : 0.f);
      float parLbp3 = -0.107f + ((iPar == 4 || iPar == 11) ? sign * 0.002f : 0.f);
      float parBsA = 1.f + ((iPar == 5 || iPar == 12) ? sign * 0.043f : 0.f);
      float parBsp1 = 0.119f + ((iPar == 6 || iPar == 13) ? sign * 0.001f : 0.f);
      float parBsp2 = -0.00091f + ((iPar == 7 || iPar == 14) ? sign * 0.00025f : 0.f);
      float parBsAvePt = 10.1f;

      for (int iBHad{0}; iBHad < NBeautyHadrons; ++iBHad) {
        fragFracFuncs[iBHad][iPar]->SetParameters(parBsA, parBsp1, parBsp2, parBsAvePt, parLbA, parLbp1, parLbp2, parLbp3);
      }
    }
  }

  std::array<float, NBeautyHadrons> beautyHadMasses{};
  for (auto iFonll{0}; iFonll < 3; ++iFonll) {
    for (auto iBHad{0}; iBHad < NBeautyHadrons; ++iBHad) {
      beautyHadMasses[iBHad] = TDatabasePDG::Instance()->GetParticle(beautyHadPdgs[iBHad])->Mass();
      for (auto iDecay{0}; iDecay < nDecays; ++iDecay) {
        auto ptB = hFonllBhad[iFonll]->GetRandom();
        auto yB = gRandom->Uniform(-1., 1.); // we might consider to use more realistic shape from FONLL in the future
        auto phiB = gRandom->Rndm() * o2::constants::math::TwoPI;
        auto pxB = ptB * std::cos(phiB);
        auto pyB = ptB * std::sin(phiB);
        auto mtB = std::sqrt(beautyHadMasses[iBHad] * beautyHadMasses[iBHad] + ptB * ptB);
        auto pzB = mtB * std::sinh(yB);
        auto pB = std::sqrt(ptB * ptB + pzB * pzB);
        auto eB = std::sqrt(beautyHadMasses[iBHad] * beautyHadMasses[iBHad] + pB * pB);

        Pythia8::Particle bHad;
        bHad.id(beautyHadPdgs[iBHad]);
        bHad.status(81);
        bHad.m(beautyHadMasses[iBHad]);
        bHad.xProd(0.);
        bHad.yProd(0.);
        bHad.zProd(0.);
        bHad.tProd(0.);
        bHad.e(eB);
        bHad.px(pxB);
        bHad.py(pyB);
        bHad.pz(pzB);

        pythia.event.reset();
        pythia.event.append(bHad);
        auto idPart = pythia.event[1].id();
        pythia.particleData.mayDecay(idPart, true);
        pythia.moreDecays();

        auto fracB = fragFracs[iBHad];
        if (fragFracOpt == LHCb) {
          fracB = fragFracFuncs[iBHad][0]->Eval(ptB > 5.f ? ptB : 5);
        } else if (fragFracOpt == LHCbMin) {
          fracB = 2.f;
          for (int iPar{0}; iPar < 15; ++iPar) {
            auto tmpFrac = fragFracFuncs[iBHad][iPar]->Eval(ptB > 5.f ? ptB : 5);
            if (tmpFrac < fracB) {
              fracB = tmpFrac;
            }
          }
        } else if (fragFracOpt == LHCbMax) {
          fracB = -1.f;
          for (int iPar{0}; iPar < 15; ++iPar) {
            auto tmpFrac = fragFracFuncs[iBHad][iPar]->Eval(ptB > 5.f ? ptB : 5);
            if (tmpFrac > fracB) {
              fracB = tmpFrac;
            }
          }
        }

        for (int iPart{1}; iPart < pythia.event.size(); ++iPart) {
          if (std::abs(pythia.event[iPart].y()) > 0.5) {
            continue;
          }
          auto absPdg = std::abs(pythia.event[iPart].id());
          if (std::find(charmHadPdgs.begin(), charmHadPdgs.end(), absPdg) != charmHadPdgs.end()) { // we found a charm hadron, let's fill the corresponding histogram
            hFonllPythiaNonPromptChad[absPdg][iBHad][iFonll]->Fill(pythia.event[iPart].pT(), fracB);
          }
        }
      }
    }
  }

  std::array<float, 3> normCrossSec{};
  for (auto iFonll{0}; iFonll < 3; ++iFonll) {
    normCrossSec[iFonll] = hFonllBhad[iFonll]->Integral();
    for (auto iChad{0}; iChad < NCharmHadrons; ++iChad) {
      hFonllPythiaNonPromptChad[charmHadPdgs[iChad]][NBeautyHadrons][iFonll] = reinterpret_cast<TH1D*>(hFonllPythiaNonPromptChad[charmHadPdgs[iChad]][0][iFonll]->Clone(Form("hFonllNonPrompt%s%s", charmHadNames[iChad].data(), namesFonll[iFonll].data())));
      hFonllPythiaNonPromptChad[charmHadPdgs[iChad]][NBeautyHadrons][iFonll]->Reset();
      for (auto iBHad{0}; iBHad < NBeautyHadrons; ++iBHad) {
        hFonllPythiaNonPromptChad[charmHadPdgs[iChad]][iBHad][iFonll]->Scale(normCrossSec[iFonll] / nDecays);
        hFonllPythiaNonPromptChad[charmHadPdgs[iChad]][NBeautyHadrons][iFonll]->Add(hFonllPythiaNonPromptChad[charmHadPdgs[iChad]][iBHad][iFonll]);
      }
    }
  }

  TFile outFile(outFileName.data(), "recreate");
  for (auto iFonll{0}; iFonll < 3; ++iFonll) {
    hFonllBhad[iFonll]->Write();
  }
  auto* dirNonPrompt = new TDirectoryFile("NonPrompt", "NonPrompt");
  dirNonPrompt->Write();
  for (auto iChad{0}; iChad < NCharmHadrons; ++iChad) {
    dirNonPrompt->cd();
    auto* dirCharmHad = new TDirectoryFile(charmHadNames[iChad].data(), charmHadNames[iChad].data());
    dirCharmHad->Write();
    dirCharmHad->cd();
    for (auto iBHad{0}; iBHad < NBeautyHadrons + 1; ++iBHad) {
      for (auto iFonll{0}; iFonll < 3; ++iFonll) {
        hFonllPythiaNonPromptChad[charmHadPdgs[iChad]][iBHad][iFonll]->Write();
      }
    }
  }
  if (addPromptCharmHadrons) {
    outFile.cd();
    auto* dirPrompt = new TDirectoryFile("Prompt", "Prompt");
    dirPrompt->Write();
    for (auto iChad{0}; iChad < NCharmHadrons; ++iChad) {
      dirPrompt->cd();
      auto* dirCharmHad = new TDirectoryFile(charmHadNames[iChad].data(), charmHadNames[iChad].data());
      dirCharmHad->Write();
      dirCharmHad->cd();
      for (auto iFonll{0}; iFonll < 3; ++iFonll) {
        if (hFonllPromptChad[charmHadPdgs[iChad]][iFonll] != nullptr) {
          hFonllPromptChad[charmHadPdgs[iChad]][iFonll]->Write();
        }
      }
    }
  }
  outFile.Close();
}
