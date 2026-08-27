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

// NOLINTBEGIN(readability/braces)
#include "OutputHandler.h"

#include <algorithm>
#include <cstdio>
#include <iostream>
#include <map>
#include <string>
#include <unordered_set>
#include <vector>

ClassImp(OutputHandler);

namespace
{
void deleteSubsampleHistograms(const std::vector<TH1D*>& sampleHistograms)
{
  // StatHandler::getNominal() returns sampleHistograms[0]. Keep that histogram
  // for the caller and release only the bootstrap/subsample histograms.
  for (std::size_t i = 1; i < sampleHistograms.size(); ++i) {
    delete sampleHistograms[i];
  }
}
} // namespace

OutputHandler::OutputHandler(std::string system) : TNamed("", ""),
                                                   flowcont(nullptr),
                                                   flowptcont(nullptr),
                                                   stathandler(nullptr),
                                                   sysHandler(nullptr),
                                                   fileName(""),
                                                   isInitialised(0),
                                                   isGen(0),
                                                   useFSRebin(true),
                                                   useLateCentralityRebin(false),
                                                   doSystemRatio(false),
                                                   radialFlowAnalysis(0),
                                                   fNrb(0),
                                                   fBinsrb(0),
                                                   xAxisTitle("Centrality (%)")
{
  calculator.setSystem(system);
}
OutputHandler::OutputHandler(std::string filename, std::string system) : TNamed("", ""),
                                                                         flowcont(nullptr),
                                                                         flowptcont(nullptr),
                                                                         stathandler(nullptr),
                                                                         sysHandler(nullptr),
                                                                         fileName(filename),
                                                                         isInitialised(0),
                                                                         isGen(0),
                                                                         useFSRebin(true),
                                                                         useLateCentralityRebin(false),
                                                                         doSystemRatio(false),
                                                                         radialFlowAnalysis(0),
                                                                         fNrb(0),
                                                                         fBinsrb(0),
                                                                         xAxisTitle("Centrality (%)")
{
  calculator.setSystem(system);
  initialise(filename);
}
OutputHandler::~OutputHandler()
{
  for (auto& [index, variation] : additionalSystematicVariations) {
    delete variation;
  }
  delete stathandler;
  delete flowcont;
  delete flowptcont;
  delete sysHandler;
  for (auto& [name, file] : inputFiles) {
    if (file) {
      file->Close();
      delete file;
    }
  }
};

void OutputHandler::setInvMassFitFunctions(const std::string& species, const std::string& signalFunction, const std::string& backgroundFunction)
{
  calculator.setInvMassFitFunctions(species, signalFunction, backgroundFunction);
}

void OutputHandler::setInvMassFitFunctionsForSystematic(int systematicIndex,
                                                        const std::string& species,
                                                        const std::string& signalFunction,
                                                        const std::string& backgroundFunction)
{
  calculator.setInvMassFitFunctionsForSystematic(systematicIndex, species, signalFunction, backgroundFunction);
}

void OutputHandler::setInvMassFitCentralityBins(const std::vector<double>& centralityBins)
{
  calculator.setInvMassFitCentralityBins(centralityBins);
}

void OutputHandler::setAdditionalSystematicVariation(int index, TH1* variation)
{
  if (!variation) {
    return;
  }
  auto existing = additionalSystematicVariations.find(index);
  if (existing != additionalSystematicVariations.end()) {
    delete existing->second;
  }
  additionalSystematicVariations[index] = variation;
}

void OutputHandler::setSystematicVariations(const std::vector<SysHandler::SystematicVariation>& variations)
{
  systematicVariations = variations;
  hasCustomSystematicVariations = true;
  if (sysHandler) {
    sysHandler->setSystematicVariations(systematicVariations);
  }
}

void OutputHandler::initialise(std::string filename, std::string dirname, bool hasSystematics, SysHandler::SystSettings settings)
{
  delete stathandler;
  stathandler = nullptr;
  delete sysHandler;
  sysHandler = nullptr;
  isInitialised = false;
  if (hasSystematics) {
    sys_indices = getKeyIndices(filename);
    sysHandler = new SysHandler();
    sysHandler->setSystSettings(settings);
    if (hasCustomSystematicVariations) {
      sysHandler->setSystematicVariations(systematicVariations);
    }
    if (!sysHandler) {
      printf("Outputhandler failed to initialise SysHandler object\n");
      return;
    }
  }
  fileName = filename;
  dirName = dirname;
  TFile* file = readFile(fileName);
  if (!file)
    return;
  if (!calculator.readData(file, dirname, isGen, radialFlowAnalysis)) {
    printf("Outputhandler failed to initialise StatHandler object\n");
    return;
  }
  stathandler = new StatHandler();
  isInitialised = 1;
  printf("Outputhandler initialised succesfully\n");
  return;
}
TH1D* OutputHandler::getHist(Calculator::ObsFunc func, StatHandler::StatMethod method)
{
  if (!isInitialised) {
    printf("Initialise OutputHandler before getting histograms\n");
    return nullptr;
  }
  Calculator::ObsFunc tmp_func = func;
  if (isVN(func) && fNrb)
    func = getCorrespondingCumulant(func);
  std::vector<TH1D*> samplehists;
  for (int i = -1; i < 10; ++i) {
    TH1D* samplehist = calculator.calculate(func, i);
    if (doSystemRatio && (!isVN(tmp_func) || !fNrb)) {
      calculator.readData(readFile(fileName2), dirName, isGen, radialFlowAnalysis);
      TH1D* denhist = calculator.calculate(func, i);
      samplehist->Divide(denhist);
      delete denhist;
      calculator.readData(readFile(fileName), dirName, isGen, radialFlowAnalysis);
    }
    samplehists.push_back(samplehist);
  }
  auto reth = stathandler->getNominal(samplehists, method);
  deleteSubsampleHistograms(samplehists);
  if (isVN(tmp_func) && fNrb) {
    if (useFSRebin && !doSystemRatio) {
      printf("Rebinning on cumulant level\n");
      FSRebin(&reth);
    }
    calculateVN(tmp_func, reth);
    if (doSystemRatio) {
      std::vector<TH1D*> samplehists2;
      calculator.readData(readFile(fileName2), dirName, isGen, radialFlowAnalysis);
      for (int i = -1; i < 10; ++i) {
        TH1D* samplehist2 = calculator.calculate(func, i);
        samplehists2.push_back(samplehist2);
      }
      auto reth2 = stathandler->getNominal(samplehists2, method);
      deleteSubsampleHistograms(samplehists2);
      calculator.readData(readFile(fileName), dirName, isGen, radialFlowAnalysis);
      calculateVN(tmp_func, reth2);
      reth->Divide(reth2);
      delete reth2;
      FSRebin(&reth);
    }
  }
  func = tmp_func;
  // sysHandler->applySyst(reth);
  reth->SetName(getName(func).c_str());
  if (reth->GetNbinsX() > 100)
    xAxisTitle = "N_{ch}";
  reth->SetTitle(Form("; %s;%s", xAxisTitle.c_str(), getLatexName(func).c_str()));
  if (useFSRebin && !isVN(func))
    FSRebin(&reth);
  return reth;
}
TH1D* OutputHandler::getHist(Calculator::PtDiffFunc func, double min, double max, StatHandler::StatMethod method)
{
  if (!isInitialised) {
    printf("Initialise OutputHandler before getting histograms\n");
    return nullptr;
  }
  Calculator::PtDiffFunc tmp_func = func;
  if (isVN(func))
    func = getCorrespondingCumulant(func);
  std::vector<TH1D*> samplehists;
  for (int i = -1; i < 10; ++i) {
    TH1D* samplehist = useLateCentralityRebin
                         ? calculator.calculateWithLateCentralityRebin(func, i, min, max)
                         : calculator.calculate(func, i, min, max);
    if (!samplehist) {
      printf("Failed to calculate sample %d for centrality %.1f-%.1f\n", i, min, max);
      for (auto* hist : samplehists)
        delete hist;
      return nullptr;
    }
    if (doSystemRatio) {
      calculator.readData(readFile(fileName2), dirName, isGen, radialFlowAnalysis);
      TH1D* denhist = useLateCentralityRebin
                        ? calculator.calculateWithLateCentralityRebin(func, i, min, max)
                        : calculator.calculate(func, i, min, max);
      samplehist->Divide(denhist);
      delete denhist;
      calculator.readData(readFile(fileName), dirName, isGen, radialFlowAnalysis);
    }
    samplehists.push_back(samplehist);
  }
  auto reth = stathandler->getNominal(samplehists, method);
  deleteSubsampleHistograms(samplehists);
  if (isVN(tmp_func)) {
    if (useFSRebin) {
      if (fNrb)
        printf("Rebinning on cumulant level\n");
      FSRebin(&reth);
    }
    calculateVN(tmp_func, reth, min, max);
  }
  func = tmp_func;
  reth->SetName(getName(func, min, max).c_str());
  if (reth->GetNbinsX() > 100)
    xAxisTitle = "N_{ch}";
  reth->SetTitle(Form("; %s;%s", xAxisTitle.c_str(), getLatexName(func).c_str()));
  if (useFSRebin && !isVN(func))
    FSRebin(&reth);
  return reth;
}
TH1D* OutputHandler::getSysHist(Calculator::ObsFunc func, bool use_bootstrap_for_systs, StatHandler::StatMethod method, bool relative, const TH1* nominal)
{
  if (!sysHandler) {
    printf("Warning: SysHandler is not initialised\n");
    return nullptr;
  }
  current_obs = getName(func);
  std::cout << current_obs.c_str() << std::endl;
  std::map<int, TH1*> sysMap;
  TH1* nom = nominal ? dynamic_cast<TH1*>(nominal->Clone(Form("%s_sysNominal", nominal->GetName()))) : getHist(func);
  if (nominal)
    nom->SetDirectory(nullptr);
  // if(useFSRebin) FSRebin(&nom);

  sysMap[0] = nom;
  const std::size_t systematicCount = std::count_if(sys_indices.begin(), sys_indices.end(), [this](int index) {
    return sysHandler->getTarget(index) != 0 && !additionalSystematicVariations.count(index);
  });
  std::size_t processedSystematics = 0;
  for (auto& i : sys_indices) {
    if (sysHandler->getTarget(i) == 0 || additionalSystematicVariations.count(i))
      continue;
    if (systematicProgressCallback)
      systematicProgressCallback(processedSystematics, systematicCount, i);
    std::string index = std::to_string(i);
    calculator.readData(readFile(fileName), index, isGen, radialFlowAnalysis);
    if (use_bootstrap_for_systs) {
      std::vector<TH1D*> sys_sample_hists;
      for (int j = -1; j < 10; ++j) {
        TH1D* sys_samplehist = calculator.calculate(func, j);
        if (doSystemRatio) {
          calculator.readData(readFile(fileName2), dirName, isGen, radialFlowAnalysis);
          TH1D* denhist = calculator.calculate(func, j);
          sys_samplehist->Divide(denhist);
          delete denhist;
          calculator.readData(readFile(fileName), index, isGen, radialFlowAnalysis);
        }
        sys_sample_hists.push_back(sys_samplehist);
      }
      auto reth = stathandler->getNominal(sys_sample_hists, method);
      deleteSubsampleHistograms(sys_sample_hists);
      if (useFSRebin)
        FSRebin(&reth);
      sysMap[i] = reth;
    } else {
      TH1* reth = calculator.calculate(func, -1);
      if (doSystemRatio) {
        calculator.readData(readFile(fileName2), dirName, isGen, radialFlowAnalysis);
        TH1D* denhist = calculator.calculate(func, -1);
        reth->Divide(denhist);
        delete denhist;
        calculator.readData(readFile(fileName), index, isGen, radialFlowAnalysis);
      }
      if (useFSRebin)
        FSRebin(&reth);
      sysMap[i] = reth;
    }
    if (systematicProgressCallback)
      systematicProgressCallback(++processedSystematics, systematicCount, i);
  }
  for (auto& [index, variation] : additionalSystematicVariations) {
    sysMap[index] = variation;
  }
  additionalSystematicVariations.clear();
  calculator.readData(readFile(fileName), dirName, isGen, radialFlowAnalysis);
  sysHandler->setSysMap(sysMap);
  TH1D* sysh = sysHandler->getNominal(relative);
  if (sysh->GetNbinsX() > 100)
    xAxisTitle = "N_{ch}";
  return sysh;
}
TH1D* OutputHandler::getSysHist(Calculator::PtDiffFunc func, double min, double max, bool use_bootstrap_for_systs, StatHandler::StatMethod method, bool relative, const TH1* nominal)
{
  if (!sysHandler) {
    printf("Warning: SysHandler is not initialised\n");
    return nullptr;
  }
  current_obs = getName(func, min, max);
  std::cout << current_obs.c_str() << std::endl;
  std::map<int, TH1*> sysMap;
  TH1* nom = nominal ? dynamic_cast<TH1*>(nominal->Clone(Form("%s_sysNominal", nominal->GetName()))) : getHist(func, min, max);
  if (nominal)
    nom->SetDirectory(nullptr);
  // if(useFSRebin) FSRebin(&nom);
  sysMap[0] = nom;
  const std::size_t systematicCount = std::count_if(sys_indices.begin(), sys_indices.end(), [this](int index) {
    return sysHandler->getTarget(index) != 0 && !additionalSystematicVariations.count(index);
  });
  std::size_t processedSystematics = 0;
  for (auto& i : sys_indices) {
    if (sysHandler->getTarget(i) == 0 || additionalSystematicVariations.count(i))
      continue;
    if (systematicProgressCallback)
      systematicProgressCallback(processedSystematics, systematicCount, i);
    std::string index = std::to_string(i);
    calculator.readData(readFile(fileName), index, isGen, radialFlowAnalysis);
    if (use_bootstrap_for_systs) {
      std::vector<TH1D*> sys_sample_hists;
      for (int j = -1; j < 10; ++j) {
        TH1D* sys_samplehist = useLateCentralityRebin
                                 ? calculator.calculateWithLateCentralityRebin(func, j, min, max)
                                 : calculator.calculate(func, j, min, max);
        if (doSystemRatio) {
          calculator.readData(readFile(fileName2), dirName, isGen, radialFlowAnalysis);
          TH1D* denhist = useLateCentralityRebin
                            ? calculator.calculateWithLateCentralityRebin(func, j, min, max)
                            : calculator.calculate(func, j, min, max);
          sys_samplehist->Divide(denhist);
          delete denhist;
          calculator.readData(readFile(fileName), index, isGen, radialFlowAnalysis);
        }
        sys_sample_hists.push_back(sys_samplehist);
      }
      auto reth = stathandler->getNominal(sys_sample_hists, method);
      deleteSubsampleHistograms(sys_sample_hists);
      if (useFSRebin)
        FSRebin(&reth);
      sysMap[i] = reth;
    } else {
      TH1* reth = useLateCentralityRebin
                    ? calculator.calculateWithLateCentralityRebin(func, -1, min, max)
                    : calculator.calculate(func, -1, min, max);
      if (doSystemRatio) {
        calculator.readData(readFile(fileName2), dirName, isGen, radialFlowAnalysis);
        TH1D* denhist = useLateCentralityRebin
                          ? calculator.calculateWithLateCentralityRebin(func, -1, min, max)
                          : calculator.calculate(func, -1, min, max);
        reth->Divide(denhist);
        delete denhist;
        calculator.readData(readFile(fileName), index, isGen, radialFlowAnalysis);
      }
      if (useFSRebin)
        FSRebin(&reth);
      sysMap[i] = reth;
    }
    if (systematicProgressCallback)
      systematicProgressCallback(++processedSystematics, systematicCount, i);
  }
  for (auto& [index, variation] : additionalSystematicVariations) {
    sysMap[index] = variation;
  }
  additionalSystematicVariations.clear();
  calculator.readData(readFile(fileName), dirName, isGen, radialFlowAnalysis);
  sysHandler->setSysMap(sysMap);
  TH1D* sysh = sysHandler->getNominal(relative);
  if (sysh->GetNbinsX() > 100)
    xAxisTitle = "N_{ch}";
  return sysh;
}

void OutputHandler::RebinMulti(Int_t nrb, Double_t* chbins)
{
  // New approach:
  fNrb = nrb;
  fBinsrb = chbins;
  if (useFSRebin)
    return;
  printf("Rebinning on correlator level\n");
  calculator.Rebin(nrb, chbins);
}
template <typename T>
void OutputHandler::FSRebin(T** inh)
{
  if (!fNrb)
    return; // Rebin only when needed
  T* sH = *inh;
  TString hName(sH->GetName());
  sH->SetName(Form("%s_BU", hName.Data()));
  T* tH = dynamic_cast<TH1D*>(sH->Rebin(fNrb, hName.Data(), fBinsrb));
  tH->Reset();
  for (Int_t i = 1; i <= sH->GetNbinsX(); i++) {
    Double_t bc = sH->GetBinContent(i);
    Double_t be = sH->GetBinError(i);
    if (be == 0)
      continue;
    be = (1. / (be * be));
    Double_t bind = tH->FindBin(sH->GetBinCenter(i));
    tH->SetBinContent(bind, tH->GetBinContent(bind) + bc * be);
    tH->SetBinError(bind, tH->GetBinError(bind) + be);
  }
  for (Int_t i = 1; i <= tH->GetNbinsX(); i++) {
    Double_t bc = tH->GetBinContent(i);
    Double_t be = tH->GetBinError(i);
    if (be == 0)
      continue;
    tH->SetBinContent(i, bc / be);
    tH->SetBinError(i, TMath::Sqrt(1. / be));
  }
  delete sH;
  (*inh) = tH;
};
TH1* OutputHandler::getRelSysHist(int sys_index, bool sys_from_fit, bool relative)
{
  if (!sysHandler) {
    printf("Warning: SysHandler is not initialised\n");
    return nullptr;
  }
  if (sysHandler->getSystMap().empty()) {
    printf("Warning: SysMap not filled! Call getSysHist() before this call\n");
    return nullptr;
  }
  TH1* relsys = (sys_from_fit) ? sysHandler->getSystematicsFromFit(sys_index, true, false, relative) : sysHandler->getSystematics(sys_index, relative, true, true);
  if (!relsys) {
    printf("Warning: Could not calculate systematic index %i\n", sys_index);
    return nullptr;
  }
  // FSRebin(&relsys);
  relsys->SetTitle(Form(";%s;%s", xAxisTitle.c_str(), relative ? "Rel. sys." : "Abs. sys."));
  return relsys;
}
TH1* OutputHandler::getBarlowHist(int sys_index)
{
  if (!sysHandler) {
    printf("Warning: SysHandler is not initialised\n");
    return nullptr;
  }
  if (sysHandler->getSystMap().empty()) {
    printf("Warning: SysMap not filled! Call getSysHist() before this call\n");
    return nullptr;
  }
  TH1* barlow = sysHandler->getBarlowTest(sys_index, sysHandler->getCorrelationError(sys_index));
  // FSRebin(&barlow);
  barlow->SetTitle(Form(";%s;Barlow", xAxisTitle.c_str()));
  return barlow;
}
TH1* OutputHandler::getVarNomRatioHist(int sys_index)
{
  if (!sysHandler) {
    printf("Warning: SysHandler is not initialised\n");
    return nullptr;
  }
  if (sysHandler->getSystMap().empty()) {
    printf("Warning: SysMap not filled! Call getSysHist() before this call\n");
    return nullptr;
  }
  TH1* ratio = sysHandler->getSystematicsFromFit(sys_index, false, true);
  if (useFSRebin) {
    FSRebin(&ratio);
  }
  SysHandler::SystSettings sys_settings = sysHandler->getSystSettings();
  ratio->Fit("pol0", "Q0+R", "", sys_settings.sys_avg_range.first, sys_settings.sys_avg_range.second);
  ratio->SetTitle(Form(";%s; nominal / variation", xAxisTitle.c_str()));
  return ratio;
}
TH1* OutputHandler::getVarHist(int sys_index)
{
  if (!sysHandler) {
    printf("Warning: SysHandler is not initialised\n");
    return nullptr;
  }
  if (sysHandler->getSystMap().empty()) {
    printf("Warning: SysMap not filled! Call getSysHist() before this call\n");
    return nullptr;
  }
  TH1* storedVariation = sysHandler->getVariation(sys_index);
  if (!storedVariation) {
    printf("Warning: Systematic variation %d is missing\n", sys_index);
    return nullptr;
  }
  TH1* var = dynamic_cast<TH1*>(storedVariation->Clone(Form("VariationPlot_Ind%d_%s", sys_index, storedVariation->GetName())));
  var->SetDirectory(nullptr);
  if (useFSRebin)
    FSRebin(&var);
  var->SetTitle(Form(";%s; variation", xAxisTitle.c_str()));
  return var;
}
TFile* OutputHandler::readFile(std::string name)
{
  auto cachedFile = inputFiles.find(name);
  if (cachedFile != inputFiles.end()) {
    return cachedFile->second;
  }
  TFile* infile = new TFile(name.c_str(), "READ");
  if (!infile) {
    printf("File %s not found!\n", name.c_str());
    return nullptr;
  };
  if (infile->IsZombie()) {
    printf("File %s is a zombie!\n", name.c_str());
    delete infile;
    return nullptr;
  };
  inputFiles[name] = infile;
  return infile;
}
std::vector<int> OutputHandler::getKeyIndices(std::string name)
{
  TFile* f = readFile(name);
  std::vector<int> indices;
  std::string dirname = dirName;
  for (auto&& key : *f->GetListOfKeys()) {
    std::string keyname = key->GetName();
    if (keyname.find(dirname) == std::string::npos)
      continue;
    auto tokens = split(keyname, '_');
    if (tokens.size() > 1) {
      indices.push_back(std::atoi(tokens[1].c_str()));
    }
  }
  std::unordered_set<int> seen;
  indices.erase(std::remove_if(indices.begin(), indices.end(), [&seen](int num) {
                  return !seen.insert(num).second; // Insert returns false if num is already in the set
                }),
                indices.end());
  std::sort(indices.begin(), indices.end());
  return indices;
}
std::vector<std::string> OutputHandler::split(const std::string& s, char delimiter)
{
  std::vector<std::string> tokens;
  std::stringstream ss(s);
  std::string item;
  while (std::getline(ss, item, delimiter)) {
    tokens.push_back(item);
  }
  return tokens;
}
std::string OutputHandler::getLatexName(Calculator::ObsFunc func)
{
  std::string str;
  if (func == &Calculator::getV12)
    str = "v_{1}{2}";
  if (func == &Calculator::getV22)
    str = "v_{2}{2}";
  if (func == &Calculator::getV32)
    str = "v_{3}{2}";
  if (func == &Calculator::getV42)
    str = "v_{4}{2}";
  if (func == &Calculator::getV24)
    str = "v_{2}{4}";
  if (func == &Calculator::getV26)
    str = "v_{2}{6}";
  if (func == &Calculator::getV28)
    str = "v_{2}{8}";
  if (func == &Calculator::getV24overV22)
    str = "v_{2}{4}/v_{2}{2}";
  if (func == &Calculator::getV32overV22)
    str = "v_{3}{2}/v_{2}{2}";
  if (func == &Calculator::getV24_3sub)
    str = "v_{2}{4}^{3sub}";
  if (func == &Calculator::getV2LLMR)
    str = "v_{2}{LLMR}";
  if (func == &Calculator::getV2LMMR)
    str = "v_{2}{LMMR}";
  if (func == &Calculator::getV2LMRR)
    str = "v_{2}{LMRR}";
  if (func == &Calculator::getVarV22)
    str = "Var(v_{2}^{2})";
  if (func == &Calculator::getC12)
    str = "c_{1}{2}";
  if (func == &Calculator::getC12Full)
    str = "c_{1}{2}";
  if (func == &Calculator::getC22)
    str = "c_{2}{2}";
  if (func == &Calculator::getC22Full)
    str = "c_{2}{2}";
  if (func == &Calculator::getC32)
    str = "c_{3}{2}";
  if (func == &Calculator::getC42)
    str = "c_{4}{2}";
  if (func == &Calculator::getC24)
    str = "c_{2}{4}";
  if (func == &Calculator::getC14)
    str = "c_{1}{4}";
  if (func == &Calculator::getC26)
    str = "c_{2}{6}";
  if (func == &Calculator::getC28)
    str = "c_{2}{8}";
  if (func == &Calculator::get22)
    str = "#LT#LT2#GT#GT_{2,-2}";
  if (func == &Calculator::get32)
    str = "#LT#LT2#GT#GT_{3,-3}";
  if (func == &Calculator::get42)
    str = "#LT#LT2#GT#GT_{4,-4}";
  if (func == &Calculator::get24)
    str = "#LT#LT4#GT#GT_{2,2,-2,-2}";
  if (func == &Calculator::getC24_3sub)
    str = "c_{2}{4}^{3sub}";
  if (func == &Calculator::getC2LLMR)
    str = "c_{2}{LLMR}";
  if (func == &Calculator::getC2LMMR)
    str = "c_{2}{LMMR}";
  if (func == &Calculator::getC2LMRR)
    str = "c_{2}{LMRR}";
  if (func == &Calculator::get2LR)
    str = "#LT LR#GT_{2,-2}";
  if (func == &Calculator::get2LM)
    str = "#LT LM#GT_{2,-2}";
  if (func == &Calculator::get2MR)
    str = "#LT MR#GT_{2,-2}";
  if (func == &Calculator::getMpt)
    str = "#LT[#it{p}]#GT";
  if (func == &Calculator::getMptSubAvg)
    str = "#LT[#it{p}]#GT^{avg}";
  if (func == &Calculator::getK2)
    str = "c_{2}";
  if (func == &Calculator::getK2Norm)
    str = "k_{2}";
  if (func == &Calculator::getK3)
    str = "c_{3}";
  if (func == &Calculator::getK3Norm)
    str = "k_{3}";
  if (func == &Calculator::getK4)
    str = "c_{4}";
  if (func == &Calculator::getMu2)
    str = "#mu_{2}";
  if (func == &Calculator::getDeltaPtOverPt)
    str = "#sqrt{#delta#it{p}_{T}^{2}}/#LT[#it{p}]#GT";
  if (func == &Calculator::getMu3)
    str = "#mu_{3}";
  if (func == &Calculator::getMu4)
    str = "#mu_{4}";
  if (func == &Calculator::getStdSkew)
    str = "#gamma_{#LT[#it{p}_{T}]#GT}";
  if (func == &Calculator::getIntSkew)
    str = "#Gamma_{#LT[#it{p}_{T}]#GT}";
  if (func == &Calculator::getExKur)
    str = "#kappa_{#LT[#it{p}_{T}]#GT}";
  if (func == &Calculator::getK2sub)
    str = "c_{2}^{sub}";
  if (func == &Calculator::getK2NormSub)
    str = "k_{2}^{sub}";
  if (func == &Calculator::getK3sub)
    str = "c_{3}^{sub}";
  if (func == &Calculator::getK3_3sub)
    str = "c_{3}^{3sub}";
  if (func == &Calculator::getK3Norm_3sub)
    str = "k_{3}^{3sub}";
  if (func == &Calculator::getK3subA)
    str = "c_{3}^{sub,A}";
  if (func == &Calculator::getK3subB)
    str = "c_{3}^{sub,B}";
  if (func == &Calculator::getK4sub)
    str = "c_{4}^{sub}";
  if (func == &Calculator::getK4_4sub)
    str = "c_{4}^{4sub}";
  if (func == &Calculator::getMu2sub)
    str = "#mu_{2}^{sub}";
  if (func == &Calculator::getMu3sub)
    str = "#mu_{3}^{sub}";
  if (func == &Calculator::getMu4sub)
    str = "#mu_{4}^{sub}";
  if (func == &Calculator::getStdSkewSub)
    str = "#gamma_{#LT[#it{p}_{T}]#GT}^{sub}";
  if (func == &Calculator::getIntSkewSub)
    str = "#Gamma_{#LT[#it{p}_{T}]#GT}^{sub}";
  if (func == &Calculator::getExKurSub)
    str = "#kappa_{#LT[#it{p}_{T}]#GT}^{sub}";
  if (func == &Calculator::getStdSkew_3sub)
    str = "#gamma_{#LT[#it{p}_{T}]#GT}^{3sub}";
  if (func == &Calculator::getIntSkew_3sub)
    str = "#Gamma_{#LT[#it{p}_{T}]#GT}^{3sub}";
  if (func == &Calculator::getExKur_4sub)
    str = "#kappa_{#LT[#it{p}_{T}]#GT}^{4sub}";
  if (func == &Calculator::getCV12pt)
    str = "C(v_{1}^{2},#delta#it{p}_{T})";
  if (func == &Calculator::getCV22pt)
    str = "C(v_{2}^{2},#delta#it{p}_{T})";
  if (func == &Calculator::getRhoV22pt)
    str = "#rho(v_{2}^{2},#delta#it{p}_{T})";

  if (func == &Calculator::getIntegratedV0)
    str = "v_{0}";
  if (func == &Calculator::getIntegratedV02)
    str = "v_{02}";
  return str;
}
std::string OutputHandler::getLatexName(Calculator::PtDiffFunc func)
{
  std::string str;
  if (func == &Calculator::getVPt22)
    str = "v_{2}{2}";
  if (func == &Calculator::getVPt32)
    str = "v_{3}{2}";
  if (func == &Calculator::getVPt42)
    str = "v_{4}{2}";
  if (func == &Calculator::getVPt24)
    str = "v_{2}{4}";
  if (func == &Calculator::getD22)
    str = "d_{2}{2}";
  if (func == &Calculator::getD32)
    str = "d_{3}{2}";
  if (func == &Calculator::getD42)
    str = "d_{4}{2}";
  if (func == &Calculator::getD24)
    str = "d_{2}{4}";

  if (func == &Calculator::getVPt24_3sub)
    str = "v_{2}{4}^{3sub}";
  if (func == &Calculator::getD24_3sub)
    str = "d_{2}{4}^{3sub}";
  if (func == &Calculator::getD2LpLMR)
    str = "d_{2}{L'LMR}";
  if (func == &Calculator::getD2LLMpR)
    str = "d_{2}{LLM'R}";
  if (func == &Calculator::getD2LLMRp)
    str = "d_{2}{LLMR'}";
  if (func == &Calculator::getD2LpMMR)
    str = "d_{2}{L'MMR}";
  if (func == &Calculator::getD2LMpMR)
    str = "d_{2}{LM'MR}";
  if (func == &Calculator::getD2LMMRp)
    str = "d_{2}{LMMR'}";
  if (func == &Calculator::getD2LpMRR)
    str = "d_{2}{L'MRR}";
  if (func == &Calculator::getD2LMpRR)
    str = "d_{2}{LM'RR}";
  if (func == &Calculator::getD2LMRpR)
    str = "d_{2}{LMR'R}";

  if (func == &Calculator::get2LpLMR)
    str = "#LT L'LMR#GT_{2,2,-2,-2}";
  if (func == &Calculator::get2LLMpR)
    str = "#LT LLM'R#GT_{2,2,-2,-2}";
  if (func == &Calculator::get2LLMRp)
    str = "#LT LLMR'#GT_{2,2,-2,-2}";
  if (func == &Calculator::get2LpMMR)
    str = "#LT L'MMR#GT_{2,2,-2,-2}";
  if (func == &Calculator::get2LMpMR)
    str = "#LT LM'MR#GT_{2,2,-2,-2}";
  if (func == &Calculator::get2LMMRp)
    str = "#LT LMMR'#GT_{2,2,-2,-2}";
  if (func == &Calculator::get2LpMRR)
    str = "#LT L'MRR#GT_{2,2,-2,-2}";
  if (func == &Calculator::get2LMpRR)
    str = "#LT LM'RR#GT_{2,2,-2,-2}";
  if (func == &Calculator::get2LMRpR)
    str = "#LT LMR'R#GT_{2,2,-2,-2}";

  if (func == &Calculator::get2LRp)
    str = "#LT LR'#GT_{2,-2}";
  if (func == &Calculator::get2LpR)
    str = "#LT L'R#GT_{2,-2}";
  if (func == &Calculator::get2LpM)
    str = "#LT L'M#GT_{2,-2}";
  if (func == &Calculator::get2LMp)
    str = "#LT LM'#GT_{2,-2}";
  if (func == &Calculator::get2MpR)
    str = "#LT M'R#GT_{2,-2}";
  if (func == &Calculator::get2MRp)
    str = "#LT MR'#GT_{2,-2}";

  if (func == &Calculator::getnptCh)
    str = "n(#it{p}_{T})";
  if (func == &Calculator::getnptPi)
    str = "n(#it{p}_{T})";
  if (func == &Calculator::getnptKa)
    str = "n(#it{p}_{T})";
  if (func == &Calculator::getnptPr)
    str = "n(#it{p}_{T})";
  if (func == static_cast<Calculator::PtDiffFunc>(&Calculator::getnptK0))
    str = "n(#it{p}_{T})";
  if (func == static_cast<Calculator::PtDiffFunc>(&Calculator::getnptLambda))
    str = "n(#it{p}_{T})";

  if (func == &Calculator::getV0Ch)
    str = "v_{0}";
  if (func == &Calculator::getV0Pi)
    str = "v_{0}";
  if (func == &Calculator::getV0Ka)
    str = "v_{0}";
  if (func == &Calculator::getV0Pr)
    str = "v_{0}";
  if (func == &Calculator::getV0K0)
    str = "v_{0}";
  if (func == &Calculator::getV0Lambda)
    str = "v_{0}";
  if (func == &Calculator::getV0WeightedCh || func == &Calculator::getV0WeightedPi ||
      func == &Calculator::getV0WeightedKa || func == &Calculator::getV0WeightedPr ||
      func == &Calculator::getV0WeightedK0 || func == &Calculator::getV0WeightedLambda)
    str = "v_{0}";
  if (func == &Calculator::getV0K0Alt)
    str = "v_{0}";
  if (func == &Calculator::getV0LambdaAlt)
    str = "v_{0}";

  if (func == &Calculator::getV02Ch)
    str = "v_{02}";
  if (func == &Calculator::getV02Pi)
    str = "v_{02}";
  if (func == &Calculator::getV02Ka)
    str = "v_{02}";
  if (func == &Calculator::getV02Pr)
    str = "v_{02}";
  if (func == &Calculator::getV02K0)
    str = "v_{02}";
  if (func == &Calculator::getV02Lambda)
    str = "v_{02}";
  if (func == &Calculator::getV02WeightedCh || func == &Calculator::getV02WeightedPi ||
      func == &Calculator::getV02WeightedKa || func == &Calculator::getV02WeightedPr ||
      func == &Calculator::getV02WeightedK0 || func == &Calculator::getV02WeightedLambda)
    str = "v_{02}";
  if (func == &Calculator::getV02K0Alt)
    str = "v_{02}";
  if (func == &Calculator::getV02LambdaAlt)
    str = "v_{02}";

  if (func == &Calculator::getV0CorrCh)
    str = "#LT [#it{p}_{T}]n(#it{p}_{T})#GT";
  if (func == &Calculator::getV0CorrPi)
    str = "#LT [#it{p}_{T}]n(#it{p}_{T})#GT";
  if (func == &Calculator::getV0CorrKa)
    str = "#LT [#it{p}_{T}]n(#it{p}_{T})#GT";
  if (func == &Calculator::getV0CorrPr)
    str = "#LT [#it{p}_{T}]n(#it{p}_{T})#GT";
  if (func == &Calculator::getV0CorrK0)
    str = "#LT [#it{p}_{T}]n(#it{p}_{T})#GT";
  if (func == &Calculator::getV0CorrLambda)
    str = "#LT [#it{p}_{T}]n(#it{p}_{T})#GT";

  if (func == &Calculator::getV0CovCh)
    str = "Cov([#it{p}_{T}],n(#it{p}_{T}))";
  if (func == &Calculator::getV0CovPi)
    str = "Cov([#it{p}_{T}],n(#it{p}_{T}))";
  if (func == &Calculator::getV0CovKa)
    str = "Cov([#it{p}_{T}],n(#it{p}_{T}))";
  if (func == &Calculator::getV0CovPr)
    str = "Cov([#it{p}_{T}],n(#it{p}_{T}))";
  if (func == &Calculator::getV0CovK0)
    str = "Cov([#it{p}_{T}],n(#it{p}_{T}))";
  if (func == &Calculator::getV0CovLambda)
    str = "Cov([#it{p}_{T}],n(#it{p}_{T}))";

  if (func == &Calculator::getV0MeansCh)
    str = "#LT [#it{p}_{T}]#GT#LT n(#it{p}_{T})#GT";
  if (func == &Calculator::getV0MeansPi)
    str = "#LT [#it{p}_{T}]#GT#LT n(#it{p}_{T})#GT";
  if (func == &Calculator::getV0MeansKa)
    str = "#LT [#it{p}_{T}]#GT#LT n(#it{p}_{T})#GT";
  if (func == &Calculator::getV0MeansPr)
    str = "#LT [#it{p}_{T}]#GT#LT n(#it{p}_{T})#GT";
  if (func == &Calculator::getV0MeansK0)
    str = "#LT [#it{p}_{T}]#GT#LT n(#it{p}_{T})#GT";
  if (func == &Calculator::getV0MeansLambda)
    str = "#LT [#it{p}_{T}]#GT#LT n(#it{p}_{T})#GT";

  if (func == &Calculator::getV02CorrCh)
    str = "#LT v_{2}^{2}n(#it{p}_{T})#GT";
  if (func == &Calculator::getV02CorrPi)
    str = "#LT v_{2}^{2}n(#it{p}_{T})#GT";
  if (func == &Calculator::getV02CorrKa)
    str = "#LT v_{2}^{2}n(#it{p}_{T})#GT";
  if (func == &Calculator::getV02CorrPr)
    str = "#LT v_{2}^{2}n(#it{p}_{T})#GT";
  if (func == &Calculator::getV02CorrK0)
    str = "#LT v_{2}^{2}n(#it{p}_{T})#GT";
  if (func == &Calculator::getV02CorrLambda)
    str = "#LT v_{2}^{2}n(#it{p}_{T})#GT";

  if (func == &Calculator::getV02CovCh)
    str = "Cov(v_{2}^{2},n(#it{p}_{T}))";
  if (func == &Calculator::getV02CovPi)
    str = "Cov(v_{2}^{2},n(#it{p}_{T}))";
  if (func == &Calculator::getV02CovKa)
    str = "Cov(v_{2}^{2},n(#it{p}_{T}))";
  if (func == &Calculator::getV02CovPr)
    str = "Cov(v_{2}^{2},n(#it{p}_{T}))";
  if (func == &Calculator::getV02CovK0)
    str = "Cov(v_{2}^{2},n(#it{p}_{T}))";
  if (func == &Calculator::getV02CovLambda)
    str = "Cov(v_{2}^{2},n(#it{p}_{T}))";

  if (func == &Calculator::getV02MeansCh)
    str = "#LT v_{2}^{2}#GT#LT n(#it{p}_{T})#GT";
  if (func == &Calculator::getV02MeansPi)
    str = "#LT v_{2}^{2}#GT#LT n(#it{p}_{T})#GT";
  if (func == &Calculator::getV02MeansKa)
    str = "#LT v_{2}^{2}#GT#LT n(#it{p}_{T})#GT";
  if (func == &Calculator::getV02MeansPr)
    str = "#LT v_{2}^{2}#GT#LT n(#it{p}_{T})#GT";
  if (func == &Calculator::getV02MeansK0)
    str = "#LT v_{2}^{2}#GT#LT n(#it{p}_{T})#GT";
  if (func == &Calculator::getV02MeansLambda)
    str = "#LT v_{2}^{2}#GT#LT n(#it{p}_{T})#GT";

  return str;
}
std::string OutputHandler::getName(Calculator::ObsFunc func)
{
  std::string str;
  if (func == &Calculator::getV12)
    str = "v12";
  if (func == &Calculator::getV22)
    str = "v22";
  if (func == &Calculator::getV32)
    str = "v32";
  if (func == &Calculator::getV42)
    str = "v42";
  if (func == &Calculator::getV24)
    str = "v24";
  if (func == &Calculator::getV26)
    str = "v26";
  if (func == &Calculator::getV28)
    str = "v28";
  if (func == &Calculator::getV24overV22)
    str = "v24overv22";
  if (func == &Calculator::getV32overV22)
    str = "v32overv22";
  if (func == &Calculator::getV24_3sub)
    str = "v24_3sub";
  if (func == &Calculator::getV2LLMR)
    str = "v2LLMR";
  if (func == &Calculator::getV2LMMR)
    str = "v2LMMR";
  if (func == &Calculator::getV2LMRR)
    str = "v2LMRR";
  if (func == &Calculator::getVarV22)
    str = "varv22";
  if (func == &Calculator::getC12)
    str = "c12";
  if (func == &Calculator::getC12Full)
    str = "c12full";
  if (func == &Calculator::getC22)
    str = "c22";
  if (func == &Calculator::getC22Full)
    str = "c22full";
  if (func == &Calculator::getC32)
    str = "c32";
  if (func == &Calculator::getC42)
    str = "c42";
  if (func == &Calculator::getC24)
    str = "c24";
  if (func == &Calculator::getC14)
    str = "c14";
  if (func == &Calculator::getC26)
    str = "c26";
  if (func == &Calculator::getC28)
    str = "c28";
  if (func == &Calculator::get22)
    str = "22";
  if (func == &Calculator::get32)
    str = "32";
  if (func == &Calculator::get42)
    str = "42";
  if (func == &Calculator::get24)
    str = "24";
  if (func == &Calculator::getC24_3sub)
    str = "c24_3sub";
  if (func == &Calculator::getC2LLMR)
    str = "c2LLMR";
  if (func == &Calculator::getC2LMMR)
    str = "c2LMMR";
  if (func == &Calculator::getC2LMRR)
    str = "c2LMRR";
  if (func == &Calculator::get2LR)
    str = "2LR";
  if (func == &Calculator::get2LM)
    str = "2LM";
  if (func == &Calculator::get2MR)
    str = "2MR";
  if (func == &Calculator::getMpt)
    str = "mpt";
  if (func == &Calculator::getMptSubAvg)
    str = "mptSubAvg";
  if (func == &Calculator::getK2)
    str = "k2";
  if (func == &Calculator::getK2Norm)
    str = "k2norm";
  if (func == &Calculator::getK3)
    str = "k3";
  if (func == &Calculator::getK3Norm)
    str = "k3norm";
  if (func == &Calculator::getK4)
    str = "k4";
  if (func == &Calculator::getMu2)
    str = "mu2";
  if (func == &Calculator::getDeltaPtOverPt)
    str = "deltaptoverpt";
  if (func == &Calculator::getMu3)
    str = "mu3";
  if (func == &Calculator::getMu4)
    str = "mu4";
  if (func == &Calculator::getStdSkew)
    str = "stdskew";
  if (func == &Calculator::getIntSkew)
    str = "intskew";
  if (func == &Calculator::getExKur)
    str = "exkur";
  if (func == &Calculator::getK2sub)
    str = "k2sub";
  if (func == &Calculator::getK2NormSub)
    str = "k2normsub";
  if (func == &Calculator::getK3sub)
    str = "k3sub";
  if (func == &Calculator::getK3_3sub)
    str = "k3_3sub";
  if (func == &Calculator::getK3Norm_3sub)
    str = "k3norm_3sub";
  if (func == &Calculator::getK3subA)
    str = "k3subA";
  if (func == &Calculator::getK3subB)
    str = "k3subB";
  if (func == &Calculator::getK4sub)
    str = "k4sub";
  if (func == &Calculator::getK4_4sub)
    str = "k4_4sub";
  if (func == &Calculator::getMu2sub)
    str = "mu2sub";
  if (func == &Calculator::getMu3sub)
    str = "mu3sub";
  if (func == &Calculator::getMu4sub)
    str = "mu4sub";
  if (func == &Calculator::getStdSkewSub)
    str = "stdskewsub";
  if (func == &Calculator::getIntSkewSub)
    str = "intskewsub";
  if (func == &Calculator::getExKurSub)
    str = "kursub";
  if (func == &Calculator::getStdSkew_3sub)
    str = "stdskew_3sub";
  if (func == &Calculator::getIntSkew_3sub)
    str = "intskew_3sub";
  if (func == &Calculator::getExKur_4sub)
    str = "kur_4sub";
  if (func == &Calculator::getCV12pt)
    str = "cv12dpt";
  if (func == &Calculator::getCV22pt)
    str = "cv22dpt";
  if (func == &Calculator::getRhoV22pt)
    str = "rhov22dpt";
  if (func == &Calculator::getIntegratedV0)
    str = "v0";
  if (func == &Calculator::getIntegratedV02)
    str = "v02";
  return str;
}
std::string OutputHandler::getName(Calculator::PtDiffFunc func, double min, double max, bool toplevelname)
{
  std::string str;
  if (func == &Calculator::getVPt22)
    str = toplevelname ? "vpt22" : Form("vpt22_%g%g", min, max);
  if (func == &Calculator::getVPt32)
    str = toplevelname ? "vpt32" : Form("vpt32_%g%g", min, max);
  if (func == &Calculator::getVPt42)
    str = toplevelname ? "vpt42" : Form("vpt42_%g%g", min, max);
  if (func == &Calculator::getVPt24)
    str = toplevelname ? "vpt24" : Form("vpt24_%g%g", min, max);
  if (func == &Calculator::getD22)
    str = toplevelname ? "d22" : Form("d22_%g%g", min, max);
  if (func == &Calculator::getD32)
    str = toplevelname ? "d32" : Form("d32_%g%g", min, max);
  if (func == &Calculator::getD42)
    str = toplevelname ? "d42" : Form("d42_%g%g", min, max);
  if (func == &Calculator::getD24)
    str = toplevelname ? "d24" : Form("d24_%g%g", min, max);

  if (func == &Calculator::getVPt24_3sub)
    str = toplevelname ? "vpt24_3sub" : Form("vpt24_3sub_%g%g", min, max);
  if (func == &Calculator::getD24_3sub)
    str = toplevelname ? "d24_3sub" : Form("d24_3sub_%g%g", min, max);
  if (func == &Calculator::getD2LpLMR)
    str = toplevelname ? "d2LpLMR" : Form("d2LpLMR_%g%g", min, max);
  if (func == &Calculator::getD2LLMpR)
    str = toplevelname ? "d2LLMpR" : Form("d2LLMpR_%g%g", min, max);
  if (func == &Calculator::getD2LLMRp)
    str = toplevelname ? "d2LLMRp" : Form("d2LLMRp_%g%g", min, max);
  if (func == &Calculator::getD2LpMMR)
    str = toplevelname ? "d2LpMMR" : Form("d2LpMMR_%g%g", min, max);
  if (func == &Calculator::getD2LMpMR)
    str = toplevelname ? "d2LMpMR" : Form("d2LMpMR_%g%g", min, max);
  if (func == &Calculator::getD2LMMRp)
    str = toplevelname ? "d2LMMRp" : Form("d2LMMRp_%g%g", min, max);
  if (func == &Calculator::getD2LpMRR)
    str = toplevelname ? "d2LpMRR" : Form("d2LpMRR_%g%g", min, max);
  if (func == &Calculator::getD2LMpRR)
    str = toplevelname ? "d2LMpRR" : Form("d2LMpRR_%g%g", min, max);
  if (func == &Calculator::getD2LMRpR)
    str = toplevelname ? "d2LMRpR" : Form("d2LMRpR_%g%g", min, max);

  if (func == &Calculator::get2LpLMR)
    str = toplevelname ? "2LpLMR" : Form("2LpLMR_%g%g", min, max);
  if (func == &Calculator::get2LLMpR)
    str = toplevelname ? "2LLMpR" : Form("2LLMpR_%g%g", min, max);
  if (func == &Calculator::get2LLMRp)
    str = toplevelname ? "2LLMRp" : Form("2LLMRp_%g%g", min, max);
  if (func == &Calculator::get2LpMMR)
    str = toplevelname ? "2LpMMR" : Form("2LpMMR_%g%g", min, max);
  if (func == &Calculator::get2LMpMR)
    str = toplevelname ? "2LMpMR" : Form("2LMpMR_%g%g", min, max);
  if (func == &Calculator::get2LMMRp)
    str = toplevelname ? "2LMMRp" : Form("2LMMRp_%g%g", min, max);
  if (func == &Calculator::get2LpMRR)
    str = toplevelname ? "2LpMRR" : Form("2LpMRR_%g%g", min, max);
  if (func == &Calculator::get2LMpRR)
    str = toplevelname ? "2LMpRR" : Form("2LMpRR_%g%g", min, max);
  if (func == &Calculator::get2LMRpR)
    str = toplevelname ? "2LMRpR" : Form("2LMRpR_%g%g", min, max);

  if (func == &Calculator::get2LRp)
    str = toplevelname ? "2LRp" : Form("2LRp_%g%g", min, max);
  if (func == &Calculator::get2LpR)
    str = toplevelname ? "2LpR" : Form("2LpR_%g%g", min, max);
  if (func == &Calculator::get2LpM)
    str = toplevelname ? "2LpM" : Form("2LpM_%g%g", min, max);
  if (func == &Calculator::get2LMp)
    str = toplevelname ? "2LMp" : Form("2LMp_%g%g", min, max);
  if (func == &Calculator::get2MpR)
    str = toplevelname ? "2MpR" : Form("2MpR_%g%g", min, max);
  if (func == &Calculator::get2MRp)
    str = toplevelname ? "2MRp" : Form("2MRp_%g%g", min, max);

  if (func == &Calculator::getnptCh)
    str = toplevelname ? "nptch" : Form("nptch_%g%g", min, max);
  if (func == &Calculator::getnptPi)
    str = toplevelname ? "nptpi" : Form("nptpi_%g%g", min, max);
  if (func == &Calculator::getnptKa)
    str = toplevelname ? "nptka" : Form("nptka_%g%g", min, max);
  if (func == &Calculator::getnptPr)
    str = toplevelname ? "nptpr" : Form("nptpr_%g%g", min, max);
  if (func == static_cast<Calculator::PtDiffFunc>(&Calculator::getnptK0))
    str = toplevelname ? "nptk0" : Form("nptk0_%g%g", min, max);
  if (func == static_cast<Calculator::PtDiffFunc>(&Calculator::getnptLambda))
    str = toplevelname ? "nptlambda" : Form("nptlambda_%g%g", min, max);

  if (func == &Calculator::getV0Ch)
    str = toplevelname ? "v0ch" : Form("v0ch_%g%g", min, max);
  if (func == &Calculator::getV0Pi)
    str = toplevelname ? "v0pi" : Form("v0pi_%g%g", min, max);
  if (func == &Calculator::getV0Ka)
    str = toplevelname ? "v0ka" : Form("v0ka_%g%g", min, max);
  if (func == &Calculator::getV0Pr)
    str = toplevelname ? "v0pr" : Form("v0pr_%g%g", min, max);
  if (func == &Calculator::getV0K0)
    str = toplevelname ? "v0k0" : Form("v0k0_%g%g", min, max);
  if (func == &Calculator::getV0Lambda)
    str = toplevelname ? "v0lambda" : Form("v0lambda_%g%g", min, max);
  if (func == &Calculator::getV0WeightedCh)
    str = toplevelname ? "v0weightedch" : Form("v0weightedch_%g%g", min, max);
  if (func == &Calculator::getV0WeightedPi)
    str = toplevelname ? "v0weightedpi" : Form("v0weightedpi_%g%g", min, max);
  if (func == &Calculator::getV0WeightedKa)
    str = toplevelname ? "v0weightedka" : Form("v0weightedka_%g%g", min, max);
  if (func == &Calculator::getV0WeightedPr)
    str = toplevelname ? "v0weightedpr" : Form("v0weightedpr_%g%g", min, max);
  if (func == &Calculator::getV0WeightedK0)
    str = toplevelname ? "v0weightedk0" : Form("v0weightedk0_%g%g", min, max);
  if (func == &Calculator::getV0WeightedLambda)
    str = toplevelname ? "v0weightedlambda" : Form("v0weightedlambda_%g%g", min, max);
  if (func == &Calculator::getV0K0Alt)
    str = toplevelname ? "v0k0alt" : Form("v0k0alt_%g%g", min, max);
  if (func == &Calculator::getV0LambdaAlt)
    str = toplevelname ? "v0lambdaalt" : Form("v0lambdaalt_%g%g", min, max);

  if (func == &Calculator::getV02Ch)
    str = toplevelname ? "v02ch" : Form("v02ch_%g%g", min, max);
  if (func == &Calculator::getV02Pi)
    str = toplevelname ? "v02pi" : Form("v02pi_%g%g", min, max);
  if (func == &Calculator::getV02Ka)
    str = toplevelname ? "v02ka" : Form("v02ka_%g%g", min, max);
  if (func == &Calculator::getV02Pr)
    str = toplevelname ? "v02pr" : Form("v02pr_%g%g", min, max);
  if (func == &Calculator::getV02K0)
    str = toplevelname ? "v02k0" : Form("v02k0_%g%g", min, max);
  if (func == &Calculator::getV02Lambda)
    str = toplevelname ? "v02lambda" : Form("v02lambda_%g%g", min, max);
  if (func == &Calculator::getV02WeightedCh)
    str = toplevelname ? "v02weightedch" : Form("v02weightedch_%g%g", min, max);
  if (func == &Calculator::getV02WeightedPi)
    str = toplevelname ? "v02weightedpi" : Form("v02weightedpi_%g%g", min, max);
  if (func == &Calculator::getV02WeightedKa)
    str = toplevelname ? "v02weightedka" : Form("v02weightedka_%g%g", min, max);
  if (func == &Calculator::getV02WeightedPr)
    str = toplevelname ? "v02weightedpr" : Form("v02weightedpr_%g%g", min, max);
  if (func == &Calculator::getV02WeightedK0)
    str = toplevelname ? "v02weightedk0" : Form("v02weightedk0_%g%g", min, max);
  if (func == &Calculator::getV02WeightedLambda)
    str = toplevelname ? "v02weightedlambda" : Form("v02weightedlambda_%g%g", min, max);
  if (func == &Calculator::getV02K0Alt)
    str = toplevelname ? "v02k0alt" : Form("v02k0alt_%g%g", min, max);
  if (func == &Calculator::getV02LambdaAlt)
    str = toplevelname ? "v02lambdaalt" : Form("v02lambdaalt_%g%g", min, max);

  if (func == &Calculator::getV0CorrCh)
    str = toplevelname ? "v0corrch" : Form("v0corrch_%g%g", min, max);
  if (func == &Calculator::getV0CorrPi)
    str = toplevelname ? "v0corrpi" : Form("v0corrpi_%g%g", min, max);
  if (func == &Calculator::getV0CorrKa)
    str = toplevelname ? "v0corrka" : Form("v0corrka_%g%g", min, max);
  if (func == &Calculator::getV0CorrPr)
    str = toplevelname ? "v0corrpr" : Form("v0corrpr_%g%g", min, max);
  if (func == &Calculator::getV0CorrK0)
    str = toplevelname ? "v0corrk0" : Form("v0corrk0_%g%g", min, max);
  if (func == &Calculator::getV0CorrLambda)
    str = toplevelname ? "v0corrlambda" : Form("v0corrlambda_%g%g", min, max);

  if (func == &Calculator::getV0CovCh)
    str = toplevelname ? "v0covch" : Form("v0covch_%g%g", min, max);
  if (func == &Calculator::getV0CovPi)
    str = toplevelname ? "v0covpi" : Form("v0covpi_%g%g", min, max);
  if (func == &Calculator::getV0CovKa)
    str = toplevelname ? "v0covka" : Form("v0covka_%g%g", min, max);
  if (func == &Calculator::getV0CovPr)
    str = toplevelname ? "v0covpr" : Form("v0covpr_%g%g", min, max);
  if (func == &Calculator::getV0CovK0)
    str = toplevelname ? "v0covk0" : Form("v0covk0_%g%g", min, max);
  if (func == &Calculator::getV0CovLambda)
    str = toplevelname ? "v0covlambda" : Form("v0covlambda_%g%g", min, max);

  if (func == &Calculator::getV02CorrCh)
    str = toplevelname ? "v02corrch" : Form("v02corrch_%g%g", min, max);
  if (func == &Calculator::getV02CorrPi)
    str = toplevelname ? "v02corrpi" : Form("v02corrpi_%g%g", min, max);
  if (func == &Calculator::getV02CorrKa)
    str = toplevelname ? "v02corrka" : Form("v02corrka_%g%g", min, max);
  if (func == &Calculator::getV02CorrPr)
    str = toplevelname ? "v02corrpr" : Form("v02corrpr_%g%g", min, max);
  if (func == &Calculator::getV02CorrK0)
    str = toplevelname ? "v02corrk0" : Form("v02corrk0_%g%g", min, max);
  if (func == &Calculator::getV02CorrLambda)
    str = toplevelname ? "v02corrlambda" : Form("v02corrlambda_%g%g", min, max);

  if (func == &Calculator::getV02CovCh)
    str = toplevelname ? "v02covch" : Form("v02covch_%g%g", min, max);
  if (func == &Calculator::getV02CovPi)
    str = toplevelname ? "v02covpi" : Form("v02covpi_%g%g", min, max);
  if (func == &Calculator::getV02CovKa)
    str = toplevelname ? "v02covka" : Form("v02covka_%g%g", min, max);
  if (func == &Calculator::getV02CovPr)
    str = toplevelname ? "v02covpr" : Form("v02covpr_%g%g", min, max);
  if (func == &Calculator::getV02CovK0)
    str = toplevelname ? "v02covk0" : Form("v02covk0_%g%g", min, max);
  if (func == &Calculator::getV02CovLambda)
    str = toplevelname ? "v02covlambda" : Form("v02covlambda_%g%g", min, max);

  if (func == &Calculator::getV0MeansCh)
    str = toplevelname ? "v0meansch" : Form("v0meansch_%g%g", min, max);
  if (func == &Calculator::getV0MeansPi)
    str = toplevelname ? "v0meanspi" : Form("v0meanspi_%g%g", min, max);
  if (func == &Calculator::getV0MeansKa)
    str = toplevelname ? "v0meanska" : Form("v0meanska_%g%g", min, max);
  if (func == &Calculator::getV0MeansPr)
    str = toplevelname ? "v0meanspr" : Form("v0meanspr_%g%g", min, max);
  if (func == &Calculator::getV0MeansK0)
    str = toplevelname ? "v0meansk0" : Form("v0meansk0_%g%g", min, max);
  if (func == &Calculator::getV0MeansLambda)
    str = toplevelname ? "v0meanslambda" : Form("v0meanslambda_%g%g", min, max);

  if (func == &Calculator::getV02MeansCh)
    str = toplevelname ? "v02meansch" : Form("v02meansch_%g%g", min, max);
  if (func == &Calculator::getV02MeansPi)
    str = toplevelname ? "v02meanspi" : Form("v02meanspi_%g%g", min, max);
  if (func == &Calculator::getV02MeansKa)
    str = toplevelname ? "v02meanska" : Form("v02meanska_%g%g", min, max);
  if (func == &Calculator::getV02MeansPr)
    str = toplevelname ? "v02meanspr" : Form("v02meanspr_%g%g", min, max);
  if (func == &Calculator::getV02MeansK0)
    str = toplevelname ? "v02meansk0" : Form("v02meansk0_%g%g", min, max);
  if (func == &Calculator::getV02MeansLambda)
    str = toplevelname ? "v02meanslambda" : Form("v02meanslambda_%g%g", min, max);

  return str;
}
void OutputHandler::PrintSystematics(bool relative)
{
  std::map<int, double> printmap = sysHandler->getPrintMap();
  if (!relative) {
    const std::map<int, std::vector<double>> printBinMap = sysHandler->getPrintBinMap();
    const std::map<int, TH1*> systMap = sysHandler->getSystMap();
    const auto nominalIt = systMap.find(0);
    const TH1* nominal = nominalIt != systMap.end() ? nominalIt->second : nullptr;
    for (const auto& [index, values] : printBinMap) {
      std::cout << "index " << index << ":";
      for (std::size_t bin = 0; bin < values.size(); ++bin) {
        if (nominal && static_cast<int>(bin) < nominal->GetNbinsX()) {
          std::cout << " [" << nominal->GetXaxis()->GetBinLowEdge(bin + 1) << ", "
                    << nominal->GetXaxis()->GetBinUpEdge(bin + 1) << "]: ";
        } else {
          std::cout << " bin " << bin + 1 << ": ";
        }
        std::cout << values[bin];
      }
      std::cout << " (absolute)" << std::endl;
    }
    return;
  }
  for (auto& sys : printmap) {
    std::cout << "index " << sys.first << ": " << sys.second * 100. << "%" << std::endl;
  }
}
bool OutputHandler::isVN(Calculator::ObsFunc func)
{
  // func == &Calculator::getV24 ||
  if (func == &Calculator::getV24 || func == &Calculator::getV26 || func == &Calculator::getV28 || func == &Calculator::getV24_3sub || func == &Calculator::getV2LLMR || func == &Calculator::getV2LMMR || func == &Calculator::getV2LMRR)
    return true;
  else
    return false;
}
Calculator::ObsFunc OutputHandler::getCorrespondingCumulant(Calculator::ObsFunc func)
{
  if (func == &Calculator::getV24)
    return &Calculator::getC24;
  if (func == &Calculator::getV26)
    return &Calculator::getC26;
  if (func == &Calculator::getV28)
    return &Calculator::getC28;
  if (func == &Calculator::getV24_3sub)
    return &Calculator::getC24_3sub;
  if (func == &Calculator::getV2LLMR)
    return &Calculator::getC2LLMR;
  if (func == &Calculator::getV2LMMR)
    return &Calculator::getC2LMMR;
  if (func == &Calculator::getV2LMRR)
    return &Calculator::getC2LMRR;
  return 0;
}
void OutputHandler::calculateVN(Calculator::ObsFunc func, TH1D*& inh)
{
  if (func == &Calculator::getV24 || func == &Calculator::getV24_3sub || func == &Calculator::getV2LLMR || func == &Calculator::getV2LMMR || func == &Calculator::getV2LMRR) {
    inh->Scale(-1.);
    inh = calculator.powerHistogram(inh, 0.25);
  } else if (func == &Calculator::getV26) {
    inh = calculator.powerHistogram(inh, 1. / 6);
  } else if (func == &Calculator::getV28) {
    inh->Scale(-1.);
    inh = calculator.powerHistogram(inh, 1. / 8);
  } else
    return;
}
bool OutputHandler::isVN(Calculator::PtDiffFunc func)
{
  // func == &Calculator::getVPt24 ||
  if (func == &Calculator::getVPt24_3sub)
    return true;
  else
    return false;
}
Calculator::PtDiffFunc OutputHandler::getCorrespondingCumulant(Calculator::PtDiffFunc func)
{
  if (func == &Calculator::getVPt24)
    return &Calculator::getD24;
  if (func == &Calculator::getVPt24_3sub)
    return &Calculator::getD24_3sub;
  return 0;
}
// TODO fix this rebinning
void OutputHandler::calculateVN(Calculator::PtDiffFunc func, TH1D*& inh, double min, double max)
{
  if (func == &Calculator::getVPt24) {
    inh->Scale(-1.);
    inh = calculator.powerHistogram(inh, 0.25);
  } else if (func == &Calculator::getVPt24_3sub) {
    inh->Scale(-1.);
    calculator.setFlowID("ChFull");
    TH1D* c24 = getHist(&Calculator::getC24);
    double c24_val = calculator.getIntegratedValue(c24, min, max);
    inh->Scale(std::pow(-c24_val, -3. / 4));
  } else
    return;
}
// NOLINTEND(readability/braces)
