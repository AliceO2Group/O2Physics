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

// NOLINTBEGIN(readability/casting, readability/braces)
#include "Calculator.h"

#include "TDirectoryFile.h"
#include "TF1.h"
#include "TFile.h"
#include "THashList.h"
#include "TObjString.h"
#include "TROOT.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <functional>
#include <iostream>
#include <iterator>
#include <memory>
#include <regex>
#include <string>
#include <utility>
#include <vector>

namespace
{
constexpr int ChargedSpeciesCount = 4;
constexpr int K0SpeciesIndex = ChargedSpeciesCount;
constexpr int LambdaSpeciesIndex = ChargedSpeciesCount + 1;

int fitCentEdge(double value)
{
  return static_cast<int>(std::lround(value));
}

double fitPtEdge(double value)
{
  return std::round(value * 10.) / 10.;
}

void applySignalFractionCorrection(Calculator& calculator, TH1D* totalHist, TH1D* backgroundHist, const std::string& species, double centMin, double centMax, const char* caller)
{
  if (!totalHist || !backgroundHist) {
    return;
  }

  const int centMinInt = fitCentEdge(centMin);
  const int centMaxInt = fitCentEdge(centMax);
  for (int bin = 1; bin <= totalHist->GetNbinsX(); ++bin) {
    const double ptMin = totalHist->GetXaxis()->GetBinLowEdge(bin); // fitPtEdge(totalHist->GetXaxis()->GetBinLowEdge(bin));
    const double ptMax = totalHist->GetXaxis()->GetBinUpEdge(bin);  // fitPtEdge(totalHist->GetXaxis()->GetBinUpEdge(bin));
    const double signalFraction = calculator.readInvMassFit(species, -1., -1., ptMin, ptMax, centMinInt, centMaxInt);
    if (signalFraction <= 0. || signalFraction >= 1.) {
      std::cout << "[" << caller << "] Invalid S/(S+B) value for " << species
                << " cent " << centMinInt << "-" << centMaxInt
                << " pt " << ptMin << "-" << ptMax
                << ": " << signalFraction << std::endl;
      totalHist->SetBinContent(bin, 0.);
      totalHist->SetBinError(bin, 0.);
      continue;
    }

    const double backgroundFraction = 1. - signalFraction;
    const double totalValue = totalHist->GetBinContent(bin);
    const double totalError = totalHist->GetBinError(bin);
    const double backgroundValue = backgroundHist->GetBinContent(bin);
    const double backgroundError = backgroundHist->GetBinError(bin);

    const double correctedValue = (totalValue - backgroundFraction * backgroundValue) / signalFraction;
    const double correctedError = std::sqrt(totalError * totalError + backgroundFraction * backgroundFraction * backgroundError * backgroundError) / signalFraction;
    totalHist->SetBinContent(bin, correctedValue);
    totalHist->SetBinError(bin, correctedError);
  }
}

TH1D* projectFraction(TH2D* histogram, const std::string& species, double min, double max)
{
  if (!histogram) {
    return nullptr;
  }
  std::unique_ptr<TProfile> profile(histogram->ProfileX(Form("pf_%s_mean_fraction_pt_%.1f_%.1f", species.c_str(), min, max), min, max));
  return profile ? profile->ProjectionX(Form("%s_mean_fraction_pt_%.1f_%.1f", species.c_str(), min, max)) : nullptr;
}
} // namespace

Calculator::Calculator() : flowcont(nullptr), flowptcont(nullptr), flowptcontfull(nullptr), corr_list(nullptr), flowid(""), collision_system(""), invMassFitPath(""), useTotal(false), useInverseVarianceWeight(true)
{
}
Calculator::~Calculator()
{
  clearInputData();
};

std::unique_ptr<FlowContainer> Calculator::prepareFlowContainerSample(int sampleIndex)
{
  if (sampleIndex < 0) {
    return nullptr;
  }
  auto prepared = std::unique_ptr<FlowContainer>(dynamic_cast<FlowContainer*>(flowcont->Clone("TempClone")));
  if (!prepared || !prepared->OverrideMainWithSub(sampleIndex, kFALSE)) {
    return nullptr;
  }
  Int_t nBins = 0;
  Double_t* binEdges = flowcont->GetMultiRebin(nBins);
  prepared->SetMultiRebin(nBins, binEdges);
  delete[] binEdges;
  return prepared;
}

FlowContainer* Calculator::getFlowContainerSample(int sampleIndex, std::unique_ptr<FlowContainer>& ownedContainer)
{
  if (activeFlowContainer && activeFlowContainerSample == sampleIndex) {
    return activeFlowContainer;
  }
  if (sampleIndex < 0) {
    return flowcont;
  }
  ownedContainer = prepareFlowContainerSample(sampleIndex);
  return ownedContainer.get();
}

void Calculator::setInvMassFitFunctions(const std::string& species, const std::string& signalFunction, const std::string& backgroundFunction)
{
  const std::string key = (species == "AntiLambda" || species == "LambdaAntiLambda") ? "Lambda" : species;
  invMassFitFunctions[key] = {signalFunction, backgroundFunction};
  invMassFitCache.clear();
  std::cout << "[Calculator::setInvMassFitFunctions] " << key << ": "
            << signalFunction << "+" << backgroundFunction << std::endl;
}

void Calculator::setInvMassFitFunctionsForSystematic(int index,
                                                     const std::string& species,
                                                     const std::string& signalFunction,
                                                     const std::string& backgroundFunction)
{
  const std::string key = (species == "AntiLambda" || species == "LambdaAntiLambda") ? "Lambda" : species;
  invMassFitFunctionsBySystematic[index][key] = {signalFunction, backgroundFunction};
  invMassFitCache.clear();
}

void Calculator::setInvMassFitCentralityBins(const std::vector<double>& centralityBins)
{
  if (centralityBins.empty()) {
    invMassFitCentralityBins.clear();
    return;
  }
  if (centralityBins.size() < 2 ||
      !std::is_sorted(centralityBins.begin(), centralityBins.end()) ||
      std::adjacent_find(centralityBins.begin(), centralityBins.end(), std::greater_equal<double>()) != centralityBins.end()) {
    std::cerr << "[Calculator::setInvMassFitCentralityBins] Centrality bin edges must be strictly increasing" << std::endl;
    return;
  }
  invMassFitCentralityBins = centralityBins;
}

TH2D* Calculator::getNptFile(TFile* infile,
                             const std::string& dirname,
                             const std::string& suffix,
                             const std::string& species,
                             const std::string& fallbackSuffix,
                             const std::string& profileSuffix)
{
  auto* hist = dynamic_cast<TH2D*>(
    infile->Get(Form("%s/npt%s_%s%s", dirname.c_str(), suffix.c_str(), species.c_str(), profileSuffix.c_str())));

  if (!hist && suffix != fallbackSuffix) {
    hist = dynamic_cast<TH2D*>(
      infile->Get(Form("%s/npt%s_%s%s", dirname.c_str(), fallbackSuffix.c_str(), species.c_str(), profileSuffix.c_str())));
  }

  if (hist) {
    hist->SetDirectory(0);
  }

  return hist;
}

void Calculator::loadFractionHists(TFile* infile, std::string dirname, std::string suffix, RadialFlowAnalysis& radialFlowStruct, const std::string& profileSuffix)
{
  const std::string fallbackSuffix = profileSuffix.empty() ? "" : suffix;
  std::vector<std::string> species = {"ch", "pi", "ka", "pr"};
  int i = 0;
  for (const auto& sp : species) {
    radialFlowStruct.fraction_cent[i] = getNptFile(infile, dirname, suffix, sp, fallbackSuffix, profileSuffix);
    ++i;
  }

  radialFlowStruct.fraction_cent_K0[0] = getNptFile(infile, dirname, suffix, "K0_sb1", fallbackSuffix, profileSuffix);
  radialFlowStruct.fraction_cent_K0[1] = getNptFile(infile, dirname, suffix, "K0_sig", fallbackSuffix, profileSuffix);
  radialFlowStruct.fraction_cent_K0[2] = getNptFile(infile, dirname, suffix, "K0_sb2", fallbackSuffix, profileSuffix);
  radialFlowStruct.fraction_cent_Lambda[0] = getNptFile(infile, dirname, suffix, "Lambda_sb1", fallbackSuffix, profileSuffix);
  radialFlowStruct.fraction_cent_Lambda[1] = getNptFile(infile, dirname, suffix, "Lambda_sig", fallbackSuffix, profileSuffix);
  radialFlowStruct.fraction_cent_Lambda[2] = getNptFile(infile, dirname, suffix, "Lambda_sb2", fallbackSuffix, profileSuffix);
}

void Calculator::clearFractionHists(RadialFlowAnalysis& radialFlowStruct)
{
  for (auto*& hist : radialFlowStruct.fraction_cent) {
    delete hist;
    hist = nullptr;
  }
  for (auto*& hist : radialFlowStruct.fraction_cent_K0) {
    delete hist;
    hist = nullptr;
  }
  for (auto*& hist : radialFlowStruct.fraction_cent_Lambda) {
    delete hist;
    hist = nullptr;
  }
}

void Calculator::clearInputData()
{
  activeFlowContainer = nullptr;
  activeFlowContainerSample = -2;
  hasChGap3pcW22 = false;
  corr_list = nullptr; // Owned by flowptcont.
  delete flowcont;
  flowcont = nullptr;
  delete flowptcont;
  flowptcont = nullptr;
  delete flowptcontfull;
  flowptcontfull = nullptr;
  clearFractionHists(radialFlowV02);
  clearFractionHists(radialFlowV0);
  clearFractionHists(radialFlowV02Weighted);
  clearFractionHists(radialFlowV0Weighted);
}

bool Calculator::readData(TFile* infile, std::string dirname, bool isGen, bool isRadialFlow)
{
  if (!infile) {
    return false;
  }
  clearInputData();
  systematicIndex = is_digits(dirname) ? std::atoi(dirname.c_str()) : 0;
  if (dirname.empty()) {
    // std::cout << "Accessing base directory ";
    dirname = (collision_system == "OO" || collision_system == "NeNe" || collision_system == "pO") ? "flow-gfw-light-ions" : "flow-generic-framework";
    // std::cout << dirname.c_str() << std::endl;
  } else {
    if (is_digits(dirname)) {
      std::cout << "Reading sys index " << dirname.c_str();
      dirname = padWithZeros(std::atoi(dirname.c_str()), 3);
      //(collision_system == "OO" || collision_system == "NeNe" || collision_system == "pO") ? dirname.insert(0, "flow-gfw-light-ions_") : dirname.insert(0, "flow-generic-framework_");
      (collision_system == "LI") ? dirname.insert(0, "flow-gfw-light-ions_") : dirname.insert(0, "flow-generic-framework_");
      std::cout << " in " << dirname.c_str() << std::endl;
    }
    // else {
    //   (collision_system == "OO" || collision_system == "NeNe" || collision_system == "pO") ? dirname.insert(0, "flow-gfw-light-ions_") : dirname.insert(0, "flow-generic-framework_");
    //   std::cout << "Reading from directory " << dirname.c_str() << std::endl;
    // }
  }
  TDirectory* inputDirectory = infile->GetDirectory(dirname.c_str());
  if (!inputDirectory) {
    std::cout << "No input directory " << dirname << std::endl;
    return false;
  }
  auto takeFromInputDirectory = [inputDirectory](const char* name) -> TObject* {
    TObject* object = inputDirectory->Get(name);
    if (object) {
      // The Calculator owns the active variation. Removing it from the ROOT
      // directory lets us release it when the next variation is loaded.
      inputDirectory->Remove(object);
    }
    return object;
  };

  // Get flow containers and transfer ownership from the input directory.
  flowcont = dynamic_cast<FlowContainer*>(takeFromInputDirectory(Form("FlowContainer%s", isGen ? "_gen" : "")));
  if (flowcont && flowcont->fProf) {
    const std::string label = "ChGap3pcW22";
    for (TObject* object : TRangeDynCast<TObjString>(flowcont->fProf->GetYaxis()->GetLabels())) {
      if (label == object->GetName()) {
        hasChGap3pcW22 = true;
        break;
      }
    }
  }
  // Combined processMC wagons keep reco and generated containers together.
  // Older generated-only wagons used the unsuffixed FlowPtContainer name.
  auto* generatedFlowPt = isGen ? dynamic_cast<FlowPtContainer*>(inputDirectory->Get("FlowPtContainer_gen")) : nullptr;
  const bool hasGeneratedFlowPt = generatedFlowPt && generatedFlowPt->getCovList() && generatedFlowPt->getCovList()->GetEntries() > 0;
  const char* flowPtName = hasGeneratedFlowPt ? "FlowPtContainer_gen" : "FlowPtContainer";
  flowptcont = dynamic_cast<FlowPtContainer*>(takeFromInputDirectory(flowPtName));
  if (!flowptcont)
    std::cout << "No FlowPtContainer with name " << Form("%s/%s", dirname.c_str(), flowPtName) << std::endl;
  const char* flowPtFullName = isGen && inputDirectory->Get("FlowPtContainerFull_gen") ? "FlowPtContainerFull_gen" : "FlowPtContainerFull";
  flowptcontfull = dynamic_cast<FlowPtContainer*>(takeFromInputDirectory(flowPtFullName));
  if (!flowptcont && !flowptcontfull) {
    std::cout << "No FlowPtContainer! Take care" << std::endl;
    return 1;
  }
  if (flowptcontfull)
    subeventInfo = checkSubeventInput(flowptcontfull);
  if (!flowptcontfull)
    subeventInfo = checkSubeventInput(flowptcont);

  // Get the vn-pt correlator profiles
  corr_list = flowptcont->getCovList();

  if (isRadialFlow) {
    // Generated fractional spectra from processMC are isolated under MCGen.
    // Keep the wagon root as a fallback for legacy generated-only files.
    std::string fractionDirectory = dirname;
    if (isGen && infile->Get((dirname + "/MCGen/npt_v02_K0_sig").c_str())) {
      fractionDirectory += "/MCGen";
    }
    loadFractionHists(infile, fractionDirectory, "_v02", radialFlowV02);
    loadFractionHists(infile, fractionDirectory, "_v0", radialFlowV0);
    loadFractionHists(infile, fractionDirectory, "_v02", radialFlowV02Weighted, "_w3pc");
    loadFractionHists(infile, fractionDirectory, "_v0", radialFlowV0Weighted, "_w2pc");

    if (collision_system == "MCGen")
      setUseTotal(true);
  }
  // if(fNrb) this->Rebin(fNrb,fBinEdges);
  return 1;
}
TH1D* Calculator::calculate(ObsFunc func, int i)
{
  std::unique_ptr<FlowContainer> prepared = prepareFlowContainerSample(i);
  if (i >= 0 && !prepared)
    return nullptr;
  FlowContainer* previousContainer = activeFlowContainer;
  const int previousSample = activeFlowContainerSample;
  struct ActiveContainerRestore {
    FlowContainer*& container;
    int& sample;
    FlowContainer* previousContainer;
    int previousSample;
    ~ActiveContainerRestore()
    {
      container = previousContainer;
      sample = previousSample;
    }
  } restore{activeFlowContainer, activeFlowContainerSample, previousContainer, previousSample};
  activeFlowContainer = i < 0 ? flowcont : prepared.get();
  activeFlowContainerSample = i;
  return (this->*func)(i);
}
TH1D* Calculator::calculate(PtDiffFunc func, int i, double min, double max)
{
  std::unique_ptr<FlowContainer> prepared = prepareFlowContainerSample(i);
  if (i >= 0 && !prepared)
    return nullptr;
  FlowContainer* previousContainer = activeFlowContainer;
  const int previousSample = activeFlowContainerSample;
  struct ActiveContainerRestore {
    FlowContainer*& container;
    int& sample;
    FlowContainer* previousContainer;
    int previousSample;
    ~ActiveContainerRestore()
    {
      container = previousContainer;
      sample = previousSample;
    }
  } restore{activeFlowContainer, activeFlowContainerSample, previousContainer, previousSample};
  activeFlowContainer = i < 0 ? flowcont : prepared.get();
  activeFlowContainerSample = i;
  return (this->*func)(i, min, max);
}
TH1D* Calculator::calculateWithLateCentralityRebin(PtDiffFunc func, int i, double min, double max)
{
  if (!flowcont || !flowcont->GetProfile() || min >= max) {
    std::cerr << "[Calculator::calculateWithLateCentralityRebin] Invalid FlowContainer or centrality range "
              << min << "-" << max << std::endl;
    return nullptr;
  }

  std::unique_ptr<FlowContainer> prepared = prepareFlowContainerSample(i);
  if (i >= 0 && !prepared)
    return nullptr;
  FlowContainer* previousContainer = activeFlowContainer;
  const int previousSample = activeFlowContainerSample;
  struct ActiveContainerRestore {
    FlowContainer*& container;
    int& sample;
    FlowContainer* previousContainer;
    int previousSample;
    ~ActiveContainerRestore()
    {
      container = previousContainer;
      sample = previousSample;
    }
  } restore{activeFlowContainer, activeFlowContainerSample, previousContainer, previousSample};
  activeFlowContainer = i < 0 ? flowcont : prepared.get();
  activeFlowContainerSample = i;

  TAxis* centralityAxis = flowcont->GetProfile()->GetXaxis();
  const int firstBin = centralityAxis->FindFixBin(min + 1.e-6);
  const int lastBin = centralityAxis->FindFixBin(max - 1.e-6);
  const int previousFitLookupCentMin = invMassFitLookupCentMin;
  const int previousFitLookupCentMax = invMassFitLookupCentMax;
  invMassFitLookupCentMin = fitCentEdge(min);
  invMassFitLookupCentMax = fitCentEdge(max);
  TH1D* result = nullptr;
  std::vector<double> weightedSums;
  std::vector<double> weightSums;
  std::vector<double> errorSquareSums;

  for (int centralityBin = firstBin; centralityBin <= lastBin; ++centralityBin) {
    if (centralityBin < 1 || centralityBin > centralityAxis->GetNbins())
      continue;
    const double binMin = centralityAxis->GetBinLowEdge(centralityBin);
    const double binMax = centralityAxis->GetBinUpEdge(centralityBin);
    TH1D* fineBinObservable = nullptr;
    {
      // ROOT normally registers every projection and clone in gDirectory.
      // A late-rebin calculation creates many such intermediates per native
      // centrality bin. Keep them in a short-lived directory and retain only
      // the observable returned by the requested Calculator function. Attach
      // it to gROOT so an open output file does not persist the scratch keys.
      TDirectoryFile scratchDirectory(
        Form("lateCentralityScratch_%d_%d", i, centralityBin), "", "", gROOT);
      TDirectory::TContext directoryContext(&scratchDirectory);
      fineBinObservable = (this->*func)(i, binMin, binMax);
      if (fineBinObservable) {
        fineBinObservable->SetDirectory(nullptr);
      }
    }
    if (!fineBinObservable)
      continue;

    if (!result) {
      result = dynamic_cast<TH1D*>(fineBinObservable->Clone(
        Form("%s_lateCentralityRebin_%.1f_%.1f", fineBinObservable->GetName(), min, max)));
      result->SetDirectory(nullptr);
      result->Reset();
      weightedSums.assign(result->GetNbinsX() + 2, 0.);
      weightSums.assign(result->GetNbinsX() + 2, 0.);
      errorSquareSums.assign(result->GetNbinsX() + 2, 0.);
    }
    for (int ptBin = 0; ptBin <= fineBinObservable->GetNbinsX() + 1; ++ptBin) {
      const double value = fineBinObservable->GetBinContent(ptBin);
      const double error = fineBinObservable->GetBinError(ptBin);
      if (!std::isfinite(value) || !std::isfinite(error) || error <= 0.)
        continue;
      const double weight = (useInverseVarianceWeight) ? 1. / (error * error) : 1.;
      weightedSums[ptBin] += value * weight;
      weightSums[ptBin] += weight;
      if (!useInverseVarianceWeight)
        errorSquareSums[ptBin] += error * error;
    }
    delete fineBinObservable;
  }

  if (!result) {
    invMassFitLookupCentMin = previousFitLookupCentMin;
    invMassFitLookupCentMax = previousFitLookupCentMax;
    std::cerr << "[Calculator::calculateWithLateCentralityRebin] No native centrality bins in range "
              << min << "-" << max << std::endl;
    return nullptr;
  }
  for (int ptBin = 0; ptBin <= result->GetNbinsX() + 1; ++ptBin) {
    if (weightSums[ptBin] <= 0.)
      continue;
    result->SetBinContent(ptBin, weightedSums[ptBin] / weightSums[ptBin]);
    result->SetBinError(ptBin, useInverseVarianceWeight ? std::sqrt(1. / weightSums[ptBin])
                                                        : std::sqrt(errorSquareSums[ptBin]) / weightSums[ptBin]);
  }
  invMassFitLookupCentMin = previousFitLookupCentMin;
  invMassFitLookupCentMax = previousFitLookupCentMax;
  return result;
}
void Calculator::Rebin(int nrb, double* chbins)
{
  if (!flowcont || !flowptcont || !corr_list)
    return;
  flowcont->SetMultiRebin(nrb, chbins);
  flowptcont->rebinMulti(nrb, chbins);
  TIter next(corr_list);
  TObject* obj;
  while ((obj = next())) {
    BootstrapProfile* bspf = dynamic_cast<BootstrapProfile*>(obj);
    bspf->RebinMulti(nrb, chbins);
  }
  return;
}

TH1D* Calculator::getK2Norm(int i)
{
  TH1* reth = getK2(i);
  TH1* mpt = getMpt(i);
  mpt->Multiply(mpt);
  reth->Divide(mpt);
  return dynamic_cast<TH1D*>(reth);
}
TH1D* Calculator::getK3Norm(int i)
{
  TH1* reth = getK3(i);
  TH1* mpt = getMpt(i);
  mpt->Multiply(mpt);
  mpt->Multiply(mpt);
  reth->Divide(mpt);
  return dynamic_cast<TH1D*>(reth);
}
TH1D* Calculator::getDeltaPtOverPt(int i)
{
  TH1D* mu2 = getMu2(i);
  TH1D* mpt = getMpt(i);
  TH1D* reth = powerHistogram(mu2, 0.5);
  reth->Divide(mpt);
  return dynamic_cast<TH1D*>(reth);
}

TH1D* Calculator::getStdSkew(int i)
{
  TH1D* reth = getK3(i);
  TH1D* k2 = getK2(i);
  reth->Divide(powerHistogram(k2, 3. / 2));
  reth->SetName(Form("stdskew%i", i));
  return reth;
}
TH1D* Calculator::getIntSkew(int i)
{
  TH1D* reth = getK3(i);
  TH1D* k2 = getK2(i);
  TH1D* k1 = getMpt(i);
  reth->Multiply(k1);
  reth->Divide(powerHistogram(k2, 2));
  reth->SetName(Form("intskew%i", i));
  return reth;
}
TH1D* Calculator::getExKur(int i)
{
  TH1D* reth = getK4(i);
  TH1D* k2 = getK2(i);
  reth->Divide(powerHistogram(k2, 2));
  reth->SetName(Form("exkur%i", i));
  return reth;
}
TH1D* Calculator::getMptSubAvg(int i)
{
  if (subeventInfo.first < 2) {
    std::cout << "Not enough subevents (" << subeventInfo.first << ")" << std::endl;
    return nullptr;
  }
  TList* l = (TList*)((!flowptcontfull) ? flowptcont : flowptcontfull)->fSubList;
  if (subeventInfo.second == 2) {

    TH1* reth = ((BootstrapProfile*)l->FindObject("ptptsub_1sub1_0sub2"))->getHist(i);
    TH1* pt1_sub2 = ((BootstrapProfile*)l->FindObject("ptptsub_0sub1_1sub2"))->getHist(i);
    reth->Add(pt1_sub2, 1);
    reth->Scale(0.5);
    reth->SetName(Form("MptSubAvg%i", i));
    return dynamic_cast<TH1D*>(reth);
  }
  if (subeventInfo.second == 3) {
    TH1* reth = ((BootstrapProfile*)l->FindObject("ptptsub_1sub1_0sub2_0sub3"))->getHist(i);
    TH1* pt1_sub2 = ((BootstrapProfile*)l->FindObject("ptptsub_0sub1_0sub2_1sub3"))->getHist(i);
    reth->Add(pt1_sub2, 1);
    reth->Scale(0.5);
    reth->SetName(Form("MptSubAvg%i", i));
    return dynamic_cast<TH1D*>(reth);
  }
  if (subeventInfo.second == 4) {
    TH1* reth = ((BootstrapProfile*)l->FindObject("ptptsub_1sub1_0sub2_0sub3_0sub4"))->getHist(i);
    TH1* pt1_sub2 = ((BootstrapProfile*)l->FindObject("ptptsub_0sub1_0sub2_0sub3_1sub4"))->getHist(i);
    reth->Add(pt1_sub2, 1);
    reth->Scale(0.5);
    reth->SetName(Form("MptSubAvg%i", i));
    return dynamic_cast<TH1D*>(reth);
  }
  TH1* reth = ((BootstrapProfile*)l->FindObject("mpt_1sub1_1sub2_2par"))->getHist(i);
  TH1* pt1_sub1 = ((BootstrapProfile*)l->FindObject("mpt_sub1_1par"))->getHist(i);
  TH1* pt1_sub2 = ((BootstrapProfile*)l->FindObject("mpt_sub2_1par"))->getHist(i);
  pt1_sub1->Multiply(pt1_sub2);
  reth->Add(pt1_sub1, -1);
  reth->SetName(Form("k2sub%i", i));
  return dynamic_cast<TH1D*>(reth);
}
// Subevent pt-pt correlations
TH1D* Calculator::getK2sub(int i)
{
  if (subeventInfo.first < 2) {
    std::cout << "Not enough particles (" << subeventInfo.first << ") correlated" << std::endl;
    return nullptr;
  }
  TList* l = (TList*)((!flowptcontfull) ? flowptcont : flowptcontfull)->fSubList;
  if (subeventInfo.second == 2) {

    TH1* reth = ((BootstrapProfile*)l->FindObject("ptptsub_1sub1_1sub2"))->getHist(i);
    TH1* pt1_sub1 = ((BootstrapProfile*)l->FindObject("ptptsub_1sub1_0sub2"))->getHist(i);
    TH1* pt1_sub2 = ((BootstrapProfile*)l->FindObject("ptptsub_0sub1_1sub2"))->getHist(i);
    pt1_sub1->Multiply(pt1_sub2);
    reth->Add(pt1_sub1, -1);
    reth->SetName(Form("k2sub%i", i));
    return dynamic_cast<TH1D*>(reth);
  }
  if (subeventInfo.second == 3) {
    TH1* reth = ((BootstrapProfile*)l->FindObject("ptptsub_1sub1_0sub2_1sub3"))->getHist(i);
    TH1* pt1_sub1 = ((BootstrapProfile*)l->FindObject("ptptsub_1sub1_0sub2_0sub3"))->getHist(i);
    TH1* pt1_sub2 = ((BootstrapProfile*)l->FindObject("ptptsub_0sub1_0sub2_1sub3"))->getHist(i);
    pt1_sub1->Multiply(pt1_sub2);
    reth->Add(pt1_sub1, -1);
    reth->SetName(Form("k2sub%i", i));
    return dynamic_cast<TH1D*>(reth);
  }
  if (subeventInfo.second == 4) {
    TH1* reth = ((BootstrapProfile*)l->FindObject("ptptsub_1sub1_0sub2_0sub3_1sub4"))->getHist(i);
    TH1* pt1_sub1 = ((BootstrapProfile*)l->FindObject("ptptsub_1sub1_0sub2_0sub3_0sub4"))->getHist(i);
    TH1* pt1_sub2 = ((BootstrapProfile*)l->FindObject("ptptsub_0sub1_0sub2_0sub3_1sub4"))->getHist(i);
    pt1_sub1->Multiply(pt1_sub2);
    reth->Add(pt1_sub1, -1);
    reth->SetName(Form("k2sub%i", i));
    return dynamic_cast<TH1D*>(reth);
  }
  TH1* reth = ((BootstrapProfile*)l->FindObject("mpt_1sub1_1sub2_2par"))->getHist(i);
  TH1* pt1_sub1 = ((BootstrapProfile*)l->FindObject("mpt_sub1_1par"))->getHist(i);
  TH1* pt1_sub2 = ((BootstrapProfile*)l->FindObject("mpt_sub2_1par"))->getHist(i);
  pt1_sub1->Multiply(pt1_sub2);
  reth->Add(pt1_sub1, -1);
  reth->SetName(Form("k2sub%i", i));
  return dynamic_cast<TH1D*>(reth);
}
TH1D* Calculator::getMu2sub(int i)
{
  TList* l = (TList*)((!flowptcontfull) ? flowptcont : flowptcontfull)->fSubCMList;
  TH1* reth = ((BootstrapProfile*)l->FindObject("cm2_11sub1_11sub2"))->getHist(i);
  TH1* pt1_sub1 = ((BootstrapProfile*)l->FindObject("cm1_sub1_Mpt0"))->getHist(i);
  TH1* pt1_sub2 = ((BootstrapProfile*)l->FindObject("cm1_sub2_Mpt0"))->getHist(i);
  pt1_sub1->Multiply(pt1_sub2);
  reth->Add(pt1_sub1, -1);
  return dynamic_cast<TH1D*>(reth);
}
TH1D* Calculator::getK2NormSub(int i)
{
  TList* l = (TList*)flowptcontfull->fSubList;
  TH1* reth = getK2sub(i);
  if (subeventInfo.second == 2) {
    TH1* pt1_sub1 = ((BootstrapProfile*)l->FindObject("ptptsub_1sub1_0sub2"))->getHist(i);
    TH1* pt1_sub2 = ((BootstrapProfile*)l->FindObject("ptptsub_0sub1_1sub2"))->getHist(i);
    pt1_sub1->Multiply(pt1_sub2);
    reth->Divide(pt1_sub1);
    reth->SetName(Form("k2sub%i", i));
    return dynamic_cast<TH1D*>(reth);
  }
  if (subeventInfo.second == 3) {
    TH1* pt1_sub1 = ((BootstrapProfile*)l->FindObject("ptptsub_1sub1_0sub2_0sub3"))->getHist(i);
    TH1* pt1_sub2 = ((BootstrapProfile*)l->FindObject("ptptsub_0sub1_0sub2_1sub3"))->getHist(i);
    pt1_sub1->Multiply(pt1_sub2);
    reth->Divide(pt1_sub1);
    reth->SetName(Form("k2sub%i", i));
    return dynamic_cast<TH1D*>(reth);
  }
  if (subeventInfo.second == 4) {
    TH1* pt1_sub1 = ((BootstrapProfile*)l->FindObject("ptptsub_1sub1_0sub2_0sub3_0sub4"))->getHist(i);
    TH1* pt1_sub2 = ((BootstrapProfile*)l->FindObject("ptptsub_0sub1_0sub2_0sub3_1sub4"))->getHist(i);
    pt1_sub1->Multiply(pt1_sub2);
    reth->Divide(pt1_sub1);
    reth->SetName(Form("k2sub%i", i));
    return dynamic_cast<TH1D*>(reth);
  }
  TH1* pt1_sub1 = ((BootstrapProfile*)l->FindObject("mpt_sub1_1par"))->getHist(i);
  TH1* pt1_sub2 = ((BootstrapProfile*)l->FindObject("mpt_sub2_1par"))->getHist(i);
  pt1_sub1->Multiply(pt1_sub2);
  reth->Divide(pt1_sub1);
  reth->SetName(Form("k2normsub%i", i));
  return dynamic_cast<TH1D*>(reth);
}
TH1D* Calculator::getK3NormSub(int i)
{
  if (subeventInfo.first < 3) {
    std::cout << "Not enough particles (" << subeventInfo.first << ") correlated" << std::endl;
    return nullptr;
  }
  if (subeventInfo.second < 3) {
    std::cout << "Not enough subevents (" << subeventInfo.second << ") to calculate" << std::endl;
    return nullptr;
  }
  TList* l = (TList*)flowptcontfull->fSubList;
  TH1* reth = getK3sub(i);
  if (subeventInfo.second == 3) {
    TH1* pt1_sub1 = ((BootstrapProfile*)l->FindObject("ptptsub_1sub1_0sub2_0sub3"))->getHist(i);
    TH1* pt1_sub2 = ((BootstrapProfile*)l->FindObject("ptptsub_0sub1_1sub2_0sub3"))->getHist(i);
    TH1* pt1_sub3 = ((BootstrapProfile*)l->FindObject("ptptsub_0sub1_0sub2_1sub3"))->getHist(i);
    pt1_sub1->Multiply(pt1_sub2);
    pt1_sub1->Multiply(pt1_sub3);
    reth->Divide(pt1_sub1);
    reth->SetName(Form("k2sub%i", i));
    return dynamic_cast<TH1D*>(reth);
  }
  if (subeventInfo.second == 4) {
    TH1* pt1_sub1 = ((BootstrapProfile*)l->FindObject("ptptsub_1sub1_0sub2_0sub3_0sub4"))->getHist(i);
    TH1* pt1_sub2 = ((BootstrapProfile*)l->FindObject("ptptsub_0sub1_1sub2_0sub3_0sub4"))->getHist(i);
    TH1* pt1_sub3 = ((BootstrapProfile*)l->FindObject("ptptsub_0sub1_0sub2_1sub3_0sub4"))->getHist(i);
    TH1* pt1_sub4 = ((BootstrapProfile*)l->FindObject("ptptsub_0sub1_0sub2_0sub3_1sub4"))->getHist(i);
    TH1D* pt_sub124 = (TH1D*)pt1_sub1->Clone(Form("sub124_%s", pt1_sub1->GetName()));
    TH1D* pt_sub134 = (TH1D*)pt1_sub1->Clone(Form("sub134_%s", pt1_sub1->GetName()));
    pt_sub124->Multiply(pt1_sub2);
    pt_sub124->Multiply(pt1_sub4);
    pt_sub134->Multiply(pt1_sub3);
    pt_sub134->Multiply(pt1_sub4);
    pt_sub134->Add(pt_sub124);
    pt_sub134->Scale(0.5);
    reth->Divide(pt_sub134);
    reth->SetName(Form("k3normsub%i", i));
    return dynamic_cast<TH1D*>(reth);
  }
  TH1* pt1_sub1 = ((BootstrapProfile*)l->FindObject("mpt_sub1_1par"))->getHist(i);
  TH1* pt1_sub2 = ((BootstrapProfile*)l->FindObject("mpt_sub2_1par"))->getHist(i);
  pt1_sub1->Multiply(pt1_sub2);
  reth->Divide(pt1_sub1);
  reth->SetName(Form("k2normsub%i", i));
  return dynamic_cast<TH1D*>(reth);
}
TH1D* Calculator::getK3sub(int i)
{

  TList* l = (TList*)((!flowptcontfull) ? flowptcont : flowptcontfull)->fSubList;
  // Consider 2 particles from subevent 1
  TH1* term1_1 = ((BootstrapProfile*)l->FindObject("mpt_2sub1_1sub2_3par"))->getHist(i);

  TH1* term2_1 = ((BootstrapProfile*)l->FindObject("mpt_1sub1_1sub2_2par"))->getHist(i);
  term2_1->Multiply(((BootstrapProfile*)l->FindObject("mpt_sub1_1par"))->getHist(i));

  TH1* term3_1 = ((BootstrapProfile*)l->FindObject("mpt_sub1_2par"))->getHist(i);
  term3_1->Multiply(((BootstrapProfile*)l->FindObject("mpt_sub2_1par"))->getHist(i));

  TH1* term4_1 = ((BootstrapProfile*)l->FindObject("mpt_sub1_1par"))->getHist(i);
  term4_1->Multiply(((BootstrapProfile*)l->FindObject("mpt_sub1_1par"))->getHist(i));
  term4_1->Multiply(((BootstrapProfile*)l->FindObject("mpt_sub2_1par"))->getHist(i));

  term1_1->Add(term2_1, -2);
  term1_1->Add(term3_1, -1);
  term1_1->Add(term4_1, 2);

  // Now swap subevents to consider two particles from subevent 2
  TH1* term1_2 = ((BootstrapProfile*)l->FindObject("mpt_1sub1_2sub2_3par"))->getHist(i);

  TH1* term2_2 = ((BootstrapProfile*)l->FindObject("mpt_1sub1_1sub2_2par"))->getHist(i);
  term2_2->Multiply(((BootstrapProfile*)l->FindObject("mpt_sub2_1par"))->getHist(i));

  TH1* term3_2 = ((BootstrapProfile*)l->FindObject("mpt_sub2_2par"))->getHist(i);
  term3_2->Multiply(((BootstrapProfile*)l->FindObject("mpt_sub1_1par"))->getHist(i));

  TH1* term4_2 = ((BootstrapProfile*)l->FindObject("mpt_sub2_1par"))->getHist(i);
  term4_2->Multiply(((BootstrapProfile*)l->FindObject("mpt_sub2_1par"))->getHist(i));
  term4_2->Multiply(((BootstrapProfile*)l->FindObject("mpt_sub1_1par"))->getHist(i));

  term1_2->Add(term2_2, -2);
  term1_2->Add(term3_2, -1);
  term1_2->Add(term4_2, 2);

  // Take average of two cases
  term1_2->Add(term1_1);
  term1_2->Scale(0.5);
  term1_2->SetName(Form("k3sub_%i", i));
  return dynamic_cast<TH1D*>(term1_2);
}
TH1D* Calculator::getK3subA(int i)
{
  TList* l = (TList*)((!flowptcontfull) ? flowptcont : flowptcontfull)->fSubList;
  // Consider 2 particles from subevent 1
  TH1* term1 = ((BootstrapProfile*)l->FindObject("mpt_2sub1_1sub2_3par"))->getHist(i);

  TH1* term2 = ((BootstrapProfile*)l->FindObject("mpt_1sub1_1sub2_2par"))->getHist(i);
  term2->Multiply(((BootstrapProfile*)l->FindObject("mpt_sub1_1par"))->getHist(i));

  TH1* term3 = ((BootstrapProfile*)l->FindObject("mpt_sub1_2par"))->getHist(i);
  term3->Multiply(((BootstrapProfile*)l->FindObject("mpt_sub2_1par"))->getHist(i));

  TH1* term4 = ((BootstrapProfile*)l->FindObject("mpt_sub1_1par"))->getHist(i);
  term4->Multiply(((BootstrapProfile*)l->FindObject("mpt_sub1_1par"))->getHist(i));
  term4->Multiply(((BootstrapProfile*)l->FindObject("mpt_sub2_1par"))->getHist(i));

  term1->Add(term2, -2);
  term1->Add(term3, -1);
  term1->Add(term4, 2);
  TH1* reth = (TH1*)term1->Clone(Form("k3subA_%i", i));
  return dynamic_cast<TH1D*>(reth);
}
TH1D* Calculator::getK3subB(int i)
{
  TList* l = (TList*)((!flowptcontfull) ? flowptcont : flowptcontfull)->fSubList;

  TH1* term1 = ((BootstrapProfile*)l->FindObject("mpt_1sub1_2sub2_3par"))->getHist(i);

  TH1* term2 = ((BootstrapProfile*)l->FindObject("mpt_1sub1_1sub2_2par"))->getHist(i);
  term2->Multiply(((BootstrapProfile*)l->FindObject("mpt_sub2_1par"))->getHist(i));

  TH1* term3 = ((BootstrapProfile*)l->FindObject("mpt_sub2_2par"))->getHist(i);
  term3->Multiply(((BootstrapProfile*)l->FindObject("mpt_sub1_1par"))->getHist(i));

  TH1* term4 = ((BootstrapProfile*)l->FindObject("mpt_sub2_1par"))->getHist(i);
  term4->Multiply(((BootstrapProfile*)l->FindObject("mpt_sub2_1par"))->getHist(i));
  term4->Multiply(((BootstrapProfile*)l->FindObject("mpt_sub1_1par"))->getHist(i));

  term1->Add(term2, -2);
  term1->Add(term3, -1);
  term1->Add(term4, 2);
  TH1* reth = (TH1*)term1->Clone(Form("k3subB_%i", i));
  return dynamic_cast<TH1D*>(reth);
}
TH1D* Calculator::getK3_3sub(int i)
{
  if (subeventInfo.first < 3) {
    std::cout << "Not enough particles (" << subeventInfo.first << ") correlated" << std::endl;
    return nullptr;
  }
  if (subeventInfo.second < 3) {
    std::cout << "Not enough subevents (" << subeventInfo.second << ") to calculate" << std::endl;
    return nullptr;
  }
  TList* l = (TList*)((!flowptcontfull) ? flowptcont : flowptcontfull)->fSubList;

  if (subeventInfo.second == 3) {
    TH1* term1 = ((BootstrapProfile*)l->FindObject("ptptsub_1sub1_1sub2_1sub3"))->getHist(i);

    TH1* term2 = ((BootstrapProfile*)l->FindObject("ptptsub_1sub1_1sub2_0sub3"))->getHist(i);
    term2->Multiply(((BootstrapProfile*)l->FindObject("ptptsub_0sub1_0sub2_1sub3"))->getHist(i));

    TH1* term3 = ((BootstrapProfile*)l->FindObject("ptptsub_1sub1_0sub2_1sub3"))->getHist(i);
    term3->Multiply(((BootstrapProfile*)l->FindObject("ptptsub_0sub1_1sub2_0sub3"))->getHist(i));

    TH1* term4 = ((BootstrapProfile*)l->FindObject("ptptsub_0sub1_1sub2_1sub3"))->getHist(i);
    term4->Multiply(((BootstrapProfile*)l->FindObject("ptptsub_1sub1_0sub2_0sub3"))->getHist(i));

    TH1* term5 = ((BootstrapProfile*)l->FindObject("ptptsub_1sub1_0sub2_0sub3"))->getHist(i);
    term5->Multiply(((BootstrapProfile*)l->FindObject("ptptsub_0sub1_1sub2_0sub3"))->getHist(i));
    term5->Multiply(((BootstrapProfile*)l->FindObject("ptptsub_0sub1_0sub2_1sub3"))->getHist(i));

    term1->Add(term2, -1);
    term1->Add(term3, -1);
    term1->Add(term4, -1);
    term1->Add(term5, 2);
    TH1* reth = (TH1*)term1->Clone(Form("k3sub_3sub_%i", i));
    return dynamic_cast<TH1D*>(reth);
  }
  if (subeventInfo.second == 4) {
    TH1* term1_A = ((BootstrapProfile*)l->FindObject("ptptsub_1sub1_1sub2_0sub3_1sub4"))->getHist(i);

    TH1* term2_A = ((BootstrapProfile*)l->FindObject("ptptsub_1sub1_1sub2_0sub3_0sub4"))->getHist(i);
    term2_A->Multiply(((BootstrapProfile*)l->FindObject("ptptsub_0sub1_0sub2_0sub3_1sub4"))->getHist(i));

    TH1* term3_A = ((BootstrapProfile*)l->FindObject("ptptsub_1sub1_0sub2_0sub3_1sub4"))->getHist(i);
    term3_A->Multiply(((BootstrapProfile*)l->FindObject("ptptsub_0sub1_1sub2_0sub3_0sub4"))->getHist(i));

    TH1* term4_A = ((BootstrapProfile*)l->FindObject("ptptsub_0sub1_1sub2_0sub3_1sub4"))->getHist(i);
    term4_A->Multiply(((BootstrapProfile*)l->FindObject("ptptsub_1sub1_0sub2_0sub3_0sub4"))->getHist(i));

    TH1* term5_A = ((BootstrapProfile*)l->FindObject("ptptsub_1sub1_0sub2_0sub3_0sub4"))->getHist(i);
    term5_A->Multiply(((BootstrapProfile*)l->FindObject("ptptsub_0sub1_1sub2_0sub3_0sub4"))->getHist(i));
    term5_A->Multiply(((BootstrapProfile*)l->FindObject("ptptsub_0sub1_0sub2_0sub3_1sub4"))->getHist(i));

    term1_A->Add(term2_A, -1);
    term1_A->Add(term3_A, -1);
    term1_A->Add(term4_A, -1);
    term1_A->Add(term5_A, 2);

    TH1* term1_B = ((BootstrapProfile*)l->FindObject("ptptsub_1sub1_0sub2_1sub3_1sub4"))->getHist(i);

    TH1* term2_B = ((BootstrapProfile*)l->FindObject("ptptsub_1sub1_0sub2_1sub3_0sub4"))->getHist(i);
    term2_B->Multiply(((BootstrapProfile*)l->FindObject("ptptsub_0sub1_0sub2_0sub3_1sub4"))->getHist(i));

    TH1* term3_B = ((BootstrapProfile*)l->FindObject("ptptsub_1sub1_0sub2_0sub3_1sub4"))->getHist(i);
    term3_B->Multiply(((BootstrapProfile*)l->FindObject("ptptsub_0sub1_0sub2_1sub3_0sub4"))->getHist(i));

    TH1* term4_B = ((BootstrapProfile*)l->FindObject("ptptsub_0sub1_0sub2_1sub3_1sub4"))->getHist(i);
    term4_B->Multiply(((BootstrapProfile*)l->FindObject("ptptsub_1sub1_0sub2_0sub3_0sub4"))->getHist(i));

    TH1* term5_B = ((BootstrapProfile*)l->FindObject("ptptsub_1sub1_0sub2_0sub3_0sub4"))->getHist(i);
    term5_B->Multiply(((BootstrapProfile*)l->FindObject("ptptsub_0sub1_0sub2_1sub3_0sub4"))->getHist(i));
    term5_B->Multiply(((BootstrapProfile*)l->FindObject("ptptsub_0sub1_0sub2_0sub3_1sub4"))->getHist(i));

    term1_B->Add(term2_B, -1);
    term1_B->Add(term3_B, -1);
    term1_B->Add(term4_B, -1);
    term1_B->Add(term5_B, 2);

    term1_A->Add(term1_B);
    term1_A->Scale(0.5);
    TH1* reth = (TH1*)term1_A->Clone(Form("k3sub_3sub_%i", i));
    return dynamic_cast<TH1D*>(reth);
  }
}
TH1D* Calculator::getK3Norm_3sub(int i)
{
  if (subeventInfo.first < 3) {
    std::cout << "Not enough particles (" << subeventInfo.first << ") correlated" << std::endl;
    return nullptr;
  }
  if (subeventInfo.second < 3) {
    std::cout << "Not enough subevents (" << subeventInfo.second << ") to calculate" << std::endl;
    return nullptr;
  }
  TList* l = (TList*)((!flowptcontfull) ? flowptcont : flowptcontfull)->fSubList;
  TH1* reth = getK3_3sub(i);
  if (subeventInfo.second == 3) {
    TH1* pt1_sub1 = ((BootstrapProfile*)l->FindObject("ptptsub_1sub1_0sub2_0sub3"))->getHist(i);
    TH1* pt1_sub2 = ((BootstrapProfile*)l->FindObject("ptptsub_0sub1_1sub2_0sub3"))->getHist(i);
    TH1* pt1_sub3 = ((BootstrapProfile*)l->FindObject("ptptsub_0sub1_0sub2_1sub3"))->getHist(i);
    pt1_sub1->Multiply(pt1_sub2);
    pt1_sub1->Multiply(pt1_sub3);
    reth->Divide(pt1_sub1);
    reth->SetName(Form("k2sub%i", i));
    return dynamic_cast<TH1D*>(reth);
  }
  if (subeventInfo.second == 4) {
    TH1* pt1_sub1 = ((BootstrapProfile*)l->FindObject("ptptsub_1sub1_0sub2_0sub3_0sub4"))->getHist(i);
    TH1* pt1_sub2 = ((BootstrapProfile*)l->FindObject("ptptsub_0sub1_1sub2_0sub3_0sub4"))->getHist(i);
    TH1* pt1_sub3 = ((BootstrapProfile*)l->FindObject("ptptsub_0sub1_0sub2_1sub3_0sub4"))->getHist(i);
    TH1* pt1_sub4 = ((BootstrapProfile*)l->FindObject("ptptsub_0sub1_0sub2_0sub3_1sub4"))->getHist(i);
    TH1D* pt_sub124 = (TH1D*)pt1_sub1->Clone(Form("sub124_%s", pt1_sub1->GetName()));
    TH1D* pt_sub134 = (TH1D*)pt1_sub1->Clone(Form("sub134_%s", pt1_sub1->GetName()));
    pt_sub124->Multiply(pt1_sub2);
    pt_sub124->Multiply(pt1_sub4);
    pt_sub134->Multiply(pt1_sub3);
    pt_sub134->Multiply(pt1_sub4);
    pt_sub134->Add(pt_sub124);
    pt_sub134->Scale(0.5);
    reth->Divide(pt_sub134);
    reth->SetName(Form("k3normsub%i", i));
    return dynamic_cast<TH1D*>(reth);
  }
  TH1* pt1_sub1 = ((BootstrapProfile*)l->FindObject("mpt_sub1_1par"))->getHist(i);
  TH1* pt1_sub2 = ((BootstrapProfile*)l->FindObject("mpt_sub2_1par"))->getHist(i);
  pt1_sub1->Multiply(pt1_sub2);
  reth->Divide(pt1_sub1);
  reth->SetName(Form("k2normsub%i", i));
  return dynamic_cast<TH1D*>(reth);
}

TH1D* Calculator::getMu3sub(int i)
{
  TList* l = (TList*)((!flowptcontfull) ? flowptcont : flowptcontfull)->fSubCMList;
  // Consider 2 particles from subevent 1
  TH1* term1_1 = ((BootstrapProfile*)l->FindObject("cm3_22sub1_11sub2"))->getHist(i);

  TH1* term2_1 = ((BootstrapProfile*)l->FindObject("cm3_12sub1_11sub2"))->getHist(i);
  term2_1->Multiply(((BootstrapProfile*)l->FindObject("cm1_sub1_Mpt0"))->getHist(i));

  TH1* term3_1 = ((BootstrapProfile*)l->FindObject("cm2_sub1_Mpt0"))->getHist(i);
  term3_1->Multiply(((BootstrapProfile*)l->FindObject("cm1_sub2_Mpt0"))->getHist(i));

  TH1* term4_1 = ((BootstrapProfile*)l->FindObject("cm2_sub1_Mpt1"))->getHist(i);
  term4_1->Multiply(((BootstrapProfile*)l->FindObject("cm1_sub1_Mpt0"))->getHist(i));
  term4_1->Multiply(((BootstrapProfile*)l->FindObject("cm1_sub2_Mpt0"))->getHist(i));

  term1_1->Add(term2_1, -2);
  term1_1->Add(term3_1, -1);
  term1_1->Add(term4_1, 2);

  // Now swap subevents to consider two particles from subevent 2
  TH1* term1_2 = ((BootstrapProfile*)l->FindObject("cm3_11sub1_22sub2"))->getHist(i);

  TH1* term2_2 = ((BootstrapProfile*)l->FindObject("cm3_11sub1_12sub2"))->getHist(i);
  term2_2->Multiply(((BootstrapProfile*)l->FindObject("cm1_sub2_Mpt0"))->getHist(i));

  TH1* term3_2 = ((BootstrapProfile*)l->FindObject("cm2_sub2_Mpt0"))->getHist(i);
  term3_2->Multiply(((BootstrapProfile*)l->FindObject("cm1_sub1_Mpt0"))->getHist(i));

  TH1* term4_2 = ((BootstrapProfile*)l->FindObject("cm2_sub2_Mpt1"))->getHist(i);
  term4_2->Multiply(((BootstrapProfile*)l->FindObject("cm1_sub2_Mpt0"))->getHist(i));
  term4_2->Multiply(((BootstrapProfile*)l->FindObject("cm1_sub1_Mpt0"))->getHist(i));

  term1_2->Add(term2_2, -2);
  term1_2->Add(term3_2, -1);
  term1_2->Add(term4_2, 2);

  // Take average of two cases
  term1_2->Add(term1_1);
  term1_2->Scale(0.5);
  return dynamic_cast<TH1D*>(term1_2);
}
TH1D* Calculator::getK4sub(int i)
{
  TList* l = (TList*)((!flowptcontfull) ? flowptcont : flowptcontfull)->fSubList;
  /*   TH1* term1 = ((BootstrapProfile*)l->FindObject("corr_2sub1_2sub2_4par"))->getHist(i);

    TH1* term2 = ((BootstrapProfile*)l->FindObject("corr_2sub1_1sub2_3par"))->getHist(i);
    term2->Multiply(((BootstrapProfile*)l->FindObject("corr_sub2_1par"))->getHist(i));

    TH1* term3 = ((BootstrapProfile*)l->FindObject("corr_1sub1_2sub2_3par"))->getHist(i);
    term3->Multiply(((BootstrapProfile*)l->FindObject("corr_sub1_1par"))->getHist(i));

    TH1* term4 = ((BootstrapProfile*)l->FindObject("corr_sub1_2par"))->getHist(i);
    term4->Multiply(((BootstrapProfile*)l->FindObject("corr_sub2_1par"))->getHist(i));
    term4->Multiply(((BootstrapProfile*)l->FindObject("corr_sub2_1par"))->getHist(i));

    TH1* term5 = ((BootstrapProfile*)l->FindObject("corr_sub2_2par"))->getHist(i);
    term5->Multiply(((BootstrapProfile*)l->FindObject("corr_sub1_1par"))->getHist(i));
    term5->Multiply(((BootstrapProfile*)l->FindObject("corr_sub1_1par"))->getHist(i));

    TH1* term6 = ((BootstrapProfile*)l->FindObject("corr_1sub1_1sub2_2par"))->getHist(i);
    term6->Multiply(((BootstrapProfile*)l->FindObject("corr_sub2_1par"))->getHist(i));
    term6->Multiply(((BootstrapProfile*)l->FindObject("corr_sub1_1par"))->getHist(i));

    TH1* term7 = ((BootstrapProfile*)l->FindObject("corr_sub1_1par"))->getHist(i);
    term7->Multiply(((BootstrapProfile*)l->FindObject("corr_sub1_1par"))->getHist(i));
    term7->Multiply(((BootstrapProfile*)l->FindObject("corr_sub2_1par"))->getHist(i));
    term7->Multiply(((BootstrapProfile*)l->FindObject("corr_sub2_1par"))->getHist(i));

    term1->Add(term2,-2);
    term1->Add(term3,-2);
    term1->Add(term4,1);
    term1->Add(term5,1);
    term1->Add(term6,4);
    term1->Add(term7,-3); */

  TH1* term1 = ((BootstrapProfile*)l->FindObject("mpt_2sub1_2sub2_4par"))->getHist(i);

  TH1* term2 = ((BootstrapProfile*)l->FindObject("mpt_2sub1_1sub2_3par"))->getHist(i);
  term2->Multiply(((BootstrapProfile*)l->FindObject("mpt_sub2_1par"))->getHist(i));

  TH1* term3 = ((BootstrapProfile*)l->FindObject("mpt_1sub1_2sub2_3par"))->getHist(i);
  term3->Multiply(((BootstrapProfile*)l->FindObject("mpt_sub1_1par"))->getHist(i));

  TH1* term4 = ((BootstrapProfile*)l->FindObject("mpt_sub1_2par"))->getHist(i);
  term4->Multiply(((BootstrapProfile*)l->FindObject("mpt_sub2_2par"))->getHist(i));

  TH1* term5 = ((BootstrapProfile*)l->FindObject("mpt_1sub1_1sub2_2par"))->getHist(i);
  term5->Multiply(term5);

  TH1* term6 = ((BootstrapProfile*)l->FindObject("mpt_sub1_2par"))->getHist(i);
  term6->Multiply(((BootstrapProfile*)l->FindObject("mpt_sub2_1par"))->getHist(i));
  term6->Multiply(((BootstrapProfile*)l->FindObject("mpt_sub2_1par"))->getHist(i));

  TH1* term7 = ((BootstrapProfile*)l->FindObject("mpt_sub2_2par"))->getHist(i);
  term7->Multiply(((BootstrapProfile*)l->FindObject("mpt_sub1_1par"))->getHist(i));
  term7->Multiply(((BootstrapProfile*)l->FindObject("mpt_sub1_1par"))->getHist(i));

  TH1* term8 = ((BootstrapProfile*)l->FindObject("mpt_1sub1_1sub2_2par"))->getHist(i);
  term8->Multiply(((BootstrapProfile*)l->FindObject("mpt_sub1_1par"))->getHist(i));
  term8->Multiply(((BootstrapProfile*)l->FindObject("mpt_sub2_1par"))->getHist(i));

  TH1* term9 = ((BootstrapProfile*)l->FindObject("mpt_sub1_1par"))->getHist(i);
  term9->Multiply(((BootstrapProfile*)l->FindObject("mpt_sub1_1par"))->getHist(i));
  term9->Multiply(((BootstrapProfile*)l->FindObject("mpt_sub2_1par"))->getHist(i));
  term9->Multiply(((BootstrapProfile*)l->FindObject("mpt_sub2_1par"))->getHist(i));

  term1->Add(term2, -2);
  term1->Add(term3, -2);
  term1->Add(term4, -1);
  term1->Add(term5, -2);
  term1->Add(term6, 2);
  term1->Add(term7, 2);
  term1->Add(term8, 8);
  term1->Add(term9, -6);

  return dynamic_cast<TH1D*>(term1);
}
TH1D* Calculator::getK4_4sub(int i)
{
  if (subeventInfo.first < 4) {
    std::cout << "Not enough particles (" << subeventInfo.first << ") correlated" << std::endl;
    return nullptr;
  }
  if (subeventInfo.second < 4) {
    std::cout << "Not enough subevents (" << subeventInfo.second << ") to calculate" << std::endl;
    return nullptr;
  }
  TList* l = (TList*)((!flowptcontfull) ? flowptcont : flowptcontfull)->fSubList;
  // Consider 4 particles from subevent 1
  TH1* term1 = ((BootstrapProfile*)l->FindObject("ptptsub_1sub1_1sub2_1sub3_1sub4"))->getHist(i);

  TH1* term2 = ((BootstrapProfile*)l->FindObject("ptptsub_1sub1_1sub2_1sub3_0sub4"))->getHist(i);
  term2->Multiply(((BootstrapProfile*)l->FindObject("ptptsub_0sub1_0sub2_0sub3_1sub4"))->getHist(i));

  TH1* term3 = ((BootstrapProfile*)l->FindObject("ptptsub_1sub1_1sub2_0sub3_1sub4"))->getHist(i);
  term3->Multiply(((BootstrapProfile*)l->FindObject("ptptsub_0sub1_0sub2_1sub3_0sub4"))->getHist(i));

  TH1* term4 = ((BootstrapProfile*)l->FindObject("ptptsub_1sub1_0sub2_1sub3_1sub4"))->getHist(i);
  term4->Multiply(((BootstrapProfile*)l->FindObject("ptptsub_0sub1_1sub2_0sub3_0sub4"))->getHist(i));

  TH1* term5 = ((BootstrapProfile*)l->FindObject("ptptsub_0sub1_1sub2_1sub3_1sub4"))->getHist(i);
  term5->Multiply(((BootstrapProfile*)l->FindObject("ptptsub_1sub1_0sub2_0sub3_0sub4"))->getHist(i));

  TH1* term6 = ((BootstrapProfile*)l->FindObject("ptptsub_1sub1_1sub2_0sub3_0sub4"))->getHist(i);
  term6->Multiply(((BootstrapProfile*)l->FindObject("ptptsub_0sub1_0sub2_1sub3_1sub4"))->getHist(i));

  TH1* term7 = ((BootstrapProfile*)l->FindObject("ptptsub_1sub1_0sub2_1sub3_0sub4"))->getHist(i);
  term7->Multiply(((BootstrapProfile*)l->FindObject("ptptsub_0sub1_1sub2_0sub3_1sub4"))->getHist(i));

  TH1* term8 = ((BootstrapProfile*)l->FindObject("ptptsub_1sub1_0sub2_0sub3_1sub4"))->getHist(i);
  term8->Multiply(((BootstrapProfile*)l->FindObject("ptptsub_0sub1_1sub2_1sub3_0sub4"))->getHist(i));

  TH1* term9 = ((BootstrapProfile*)l->FindObject("ptptsub_1sub1_1sub2_0sub3_0sub4"))->getHist(i);
  term9->Multiply(((BootstrapProfile*)l->FindObject("ptptsub_0sub1_0sub2_1sub3_0sub4"))->getHist(i));
  term9->Multiply(((BootstrapProfile*)l->FindObject("ptptsub_0sub1_0sub2_0sub3_1sub4"))->getHist(i));

  TH1* term10 = ((BootstrapProfile*)l->FindObject("ptptsub_1sub1_0sub2_1sub3_0sub4"))->getHist(i);
  term10->Multiply(((BootstrapProfile*)l->FindObject("ptptsub_0sub1_1sub2_0sub3_0sub4"))->getHist(i));
  term10->Multiply(((BootstrapProfile*)l->FindObject("ptptsub_0sub1_0sub2_0sub3_1sub4"))->getHist(i));

  TH1* term11 = ((BootstrapProfile*)l->FindObject("ptptsub_1sub1_0sub2_0sub3_1sub4"))->getHist(i);
  term11->Multiply(((BootstrapProfile*)l->FindObject("ptptsub_0sub1_1sub2_0sub3_0sub4"))->getHist(i));
  term11->Multiply(((BootstrapProfile*)l->FindObject("ptptsub_0sub1_0sub2_1sub3_0sub4"))->getHist(i));

  TH1* term12 = ((BootstrapProfile*)l->FindObject("ptptsub_0sub1_1sub2_1sub3_0sub4"))->getHist(i);
  term12->Multiply(((BootstrapProfile*)l->FindObject("ptptsub_1sub1_0sub2_0sub3_0sub4"))->getHist(i));
  term12->Multiply(((BootstrapProfile*)l->FindObject("ptptsub_0sub1_0sub2_0sub3_1sub4"))->getHist(i));

  TH1* term13 = ((BootstrapProfile*)l->FindObject("ptptsub_0sub1_1sub2_0sub3_1sub4"))->getHist(i);
  term13->Multiply(((BootstrapProfile*)l->FindObject("ptptsub_1sub1_0sub2_0sub3_0sub4"))->getHist(i));
  term13->Multiply(((BootstrapProfile*)l->FindObject("ptptsub_0sub1_0sub2_1sub3_0sub4"))->getHist(i));

  TH1* term14 = ((BootstrapProfile*)l->FindObject("ptptsub_0sub1_0sub2_1sub3_1sub4"))->getHist(i);
  term14->Multiply(((BootstrapProfile*)l->FindObject("ptptsub_1sub1_0sub2_0sub3_0sub4"))->getHist(i));
  term14->Multiply(((BootstrapProfile*)l->FindObject("ptptsub_0sub1_1sub2_0sub3_0sub4"))->getHist(i));

  TH1* term15 = ((BootstrapProfile*)l->FindObject("ptptsub_1sub1_0sub2_0sub3_0sub4"))->getHist(i);
  term15->Multiply(((BootstrapProfile*)l->FindObject("ptptsub_0sub1_1sub2_0sub3_0sub4"))->getHist(i));
  term15->Multiply(((BootstrapProfile*)l->FindObject("ptptsub_0sub1_0sub2_1sub3_0sub4"))->getHist(i));
  term15->Multiply(((BootstrapProfile*)l->FindObject("ptptsub_0sub1_0sub2_0sub3_1sub4"))->getHist(i));

  term1->Add(term2, -1);
  term1->Add(term3, -1);
  term1->Add(term4, -1);
  term1->Add(term5, -1);

  term1->Add(term6, -1);
  term1->Add(term7, -1);
  term1->Add(term8, -1);

  term1->Add(term9, 2);
  term1->Add(term10, 2);
  term1->Add(term11, 2);
  term1->Add(term12, 2);
  term1->Add(term13, 2);
  term1->Add(term14, 2);
  term1->Add(term15, -6);

  TH1* reth = (TH1*)term1->Clone(Form("k4sub_4sub_%i", i));
  return dynamic_cast<TH1D*>(reth);
}

TH1D* Calculator::getMu4sub(int i)
{
  TList* l = (TList*)((!flowptcontfull) ? flowptcont : flowptcontfull)->fSubCMList;
  // Consider 2 particles from subevent 1
  TH1* term1 = ((BootstrapProfile*)l->FindObject("cm4_22sub1_11sub2"))->getHist(i);

  return dynamic_cast<TH1D*>(term1);
}
TH1D* Calculator::getStdSkewSub(int i)
{
  TH1D* reth = getK3sub(i);
  reth->SetName(Form("stdskew_k3sub_%i", i));
  TH1D* k2sub = getK2sub(i);
  k2sub->SetName(Form("stdskew_k2sub_%i", i));
  reth->Divide(powerHistogram(k2sub, 3. / 2));
  reth->SetName(Form("stdskewsub%i", i));
  return reth;
}
TH1D* Calculator::getStdSkew_3sub(int i)
{
  TH1D* reth = getK3_3sub(i);
  reth->SetName(Form("stdskew_k3_3sub_%i", i));
  // should take the appropriate subevent config for K2 based on the input file
  TH1D* k2sub = getK2sub(i);
  k2sub->SetName(Form("stdskew_k2_3sub_%i", i));
  reth->Divide(powerHistogram(k2sub, 3. / 2));
  reth->SetName(Form("stdskew_3sub%i", i));
  return reth;
}
TH1D* Calculator::getIntSkewSub(int i)
{
  TH1D* reth = getK3sub(i);
  TH1D* k2sub = getK2sub(i);
  TList* l = (TList*)flowptcontfull->fSubList;
  TH1* pt1_sub1 = ((BootstrapProfile*)l->FindObject("mpt_sub1_1par"))->getHist(i);
  TH1* pt1_sub2 = ((BootstrapProfile*)l->FindObject("mpt_sub2_1par"))->getHist(i);
  pt1_sub1->Add(pt1_sub2, 1);
  pt1_sub1->Scale(0.5);
  reth->Multiply(pt1_sub1);
  reth->Divide(powerHistogram(k2sub, 2));
  reth->SetName(Form("intskewsub%i", i));
  return reth;
}
TH1D* Calculator::getIntSkew_3sub(int i)
{
  TH1D* reth = getK3_3sub(i);
  TH1D* k2sub = getK2sub(i);
  TList* l = (TList*)flowptcontfull->fSubList;
  TH1* pt1_subMid;
  if (subeventInfo.second == 3)
    pt1_subMid = ((BootstrapProfile*)l->FindObject("ptptsub_0sub1_1sub2_0sub3"))->getHist(i);
  if (subeventInfo.second == 4) {
    pt1_subMid = ((BootstrapProfile*)l->FindObject("ptptsub_0sub1_1sub2_0sub3_0sub4"))->getHist(i);
    TH1* pt1_subMidB = ((BootstrapProfile*)l->FindObject("ptptsub_0sub1_0sub2_1sub3_0sub4"))->getHist(i);
    pt1_subMid->Add(pt1_subMidB);
    pt1_subMid->Scale(0.5);
  }
  reth->Multiply(pt1_subMid);
  reth->Divide(powerHistogram(k2sub, 2));
  reth->SetName(Form("intskew_3sub%i", i));
  return reth;
}
TH1D* Calculator::getIntSkewAvgSub(int i)
{
  TH1D* k3sub1 = getK3subA(i);
  TH1D* k2sub = getK2sub(i);
  TList* l = (TList*)flowptcontfull->fSubList;
  TH1* pt1_sub1 = ((BootstrapProfile*)l->FindObject("mpt_sub1_1par"))->getHist(i);

  k3sub1->Multiply(pt1_sub1);
  k3sub1->Divide(powerHistogram(k2sub, 2));

  TH1D* k3sub2 = getK3subB(i);
  TH1* pt1_sub2 = ((BootstrapProfile*)l->FindObject("mpt_sub2_1par"))->getHist(i);
  k3sub2->Multiply(pt1_sub2);
  k3sub2->Divide(powerHistogram(k2sub, 2));

  TH1* reth = (TH1*)k3sub1->Clone(Form("intskew_avg_sub_%i", i));
  reth->Add(k3sub2);
  reth->Scale(0.5);
  reth->SetName(Form("intskewavgsub%i", i));
  return dynamic_cast<TH1D*>(reth);
}
TH1D* Calculator::getExKurSub(int i)
{
  TH1D* reth = getK4sub(i);
  TH1D* k2sub = getK2sub(i);
  k2sub->SetName(Form("k2subExKur_%i", i));
  reth->Divide(powerHistogram(k2sub, 2));
  reth->SetName(Form("exkursub%i", i));
  return reth;
}
TH1D* Calculator::getExKur_4sub(int i)
{
  TH1D* reth = getK4_4sub(i);
  TH1D* k2sub = getK2sub(i);
  k2sub->SetName(Form("k2_4sub_ExKur_%i", i));
  reth->Divide(powerHistogram(k2sub, 2));
  reth->SetName(Form("exkur_4sub%i", i));
  return reth;
}
TH1D* Calculator::getV32overV22(int i)
{
  TH1D* reth = getVNN(3, 2, "ChGap", i);
  TH1D* v2 = getVNN(2, 2, "ChGap", i);
  reth->Divide(v2);
  reth->GetYaxis()->SetTitle("v_{3}{2}/v_{2}{2}");
  return reth;
}
TH1D* Calculator::getV24overV22(int i)
{
  TH1D* reth = getVNN(2, 4, "ChFull", i);
  TH1D* v2 = getVNN(2, 2, "ChGap", i);
  reth->Divide(v2);
  reth->GetYaxis()->SetTitle("v_{2}{4}/v_{2}{2}");
  return reth;
}
TH1D* Calculator::getVarV22(int i)
{
  auto reth = getCNN(2, 4, "ChFull", i);
  auto c22gap = getCNN(2, 2, "ChGap", i);
  c22gap->Multiply(c22gap);
  reth->Add(c22gap);
  return reth;
}
TH1D* Calculator::getV24_3sub(int i)
{
  TH1D* reth = getC24_3sub(i);
  reth->Scale(-1.);
  reth = powerHistogram(reth, 0.25);
  return reth;
}
TH1D* Calculator::getVPt24_3sub(int i, double min, double max)
{
  TH1D* reth = getD24_3sub(i, min, max);
  reth->Scale(-1.);
  TH1D* c4 = getC24_3sub(i);
  double c4_val = getIntegratedValue(c4, min, max);
  reth->Scale(std::pow(-c4_val, -3. / 4));
  return reth;
}
TH1D* Calculator::getVPt24_3sub_avgOnVn(int i, double min, double max)
{
  TH1D* hD2LpLMR = getD2LpLMR(i, min, max);
  hD2LpLMR->Scale(-1.);
  TH1D* c2LLMR = getC2LLMR(i);
  double c2LLMR_val = getIntegratedValue(c2LLMR, min, max);
  hD2LpLMR->Scale(std::pow(-c2LLMR_val, -3. / 4));

  TH1D* hD2LLMpR = getD2LLMpR(i, min, max);
  hD2LLMpR->Scale(-1.);
  hD2LLMpR->Scale(std::pow(-c2LLMR_val, -3. / 4));

  TH1D* hD2LLMRp = getD2LLMRp(i, min, max);
  hD2LLMRp->Scale(-1.);
  hD2LLMRp->Scale(std::pow(-c2LLMR_val, -3. / 4));

  TH1D* hD2LpMMR = getD2LpMMR(i, min, max);
  hD2LpMMR->Scale(-1.);
  TH1D* c2LMMR = getC2LMMR(i);
  double c2LMMR_val = getIntegratedValue(c2LMMR, min, max);
  hD2LpMMR->Scale(std::pow(-c2LMMR_val, -3. / 4));

  TH1D* hD2LMpMR = getD2LMpMR(i, min, max);
  hD2LMpMR->Scale(-1.);
  hD2LMpMR->Scale(std::pow(-c2LMMR_val, -3. / 4));

  TH1D* hD2LMMRp = getD2LMMRp(i, min, max);
  hD2LMMRp->Scale(-1.);
  hD2LMMRp->Scale(std::pow(-c2LMMR_val, -3. / 4));

  TH1D* hD2LpMRR = getD2LpMRR(i, min, max);
  hD2LpMRR->Scale(-1.);
  TH1D* c2LMRR = getC2LMRR(i);
  double c2LMRR_val = getIntegratedValue(c2LMRR, min, max);
  hD2LpMRR->Scale(std::pow(-c2LMRR_val, -3. / 4));

  TH1D* hD2LMpRR = getD2LMpRR(i, min, max);
  hD2LMpRR->Scale(-1.);
  hD2LMpRR->Scale(std::pow(-c2LMRR_val, -3. / 4));

  TH1D* hD2LMRpR = getD2LMRpR(i, min, max);
  hD2LMRpR->Scale(-1.);
  hD2LMRpR->Scale(std::pow(-c2LMRR_val, -3. / 4));

  /*   for(int bin = 1; bin<=reth->GetNbinsX(); ++bin){
      double be_LpLMR = d2LpLMR->GetBinError(bin);
      double be_LLMpR = d2LLMpR->GetBinError(bin);
      double be_LLMRp = d2LLMRp->GetBinError(bin);

      double be_LpMMR = d2LpMMR->GetBinError(bin);
      double be_LMpMR = d2LMpMR->GetBinError(bin);
      double be_LMMRp = d2LMMRp->GetBinError(bin);

      double be_LpMRR = d2LpMRR->GetBinError(bin);
      double be_LMpRR = d2LMpRR->GetBinError(bin);
      double be_LMRRp = d2LMRRp->GetBinError(bin);

      if( be_LpLMR == 0.|| be_LLMpR == 0. || be_LLMRp == 0. || be_LpMMR == 0. || be_LMpMR == 0. || be_LMMRp == 0. || be_LpMRR == 0. || be_LMpRR == 0. || be_LMRRp == 0.) { reth->SetBinContent(bin,-2.); reth->SetBinError(bin,0.); }

      double bc_LpLMR = d2LpLMR->GetBinContent(bin);
      double bc_LLMpR = d2LLMpR->GetBinContent(bin);
      double bc_LLMRp = d2LLMRp->GetBinContent(bin);

      double bc_LpMMR = d2LpMMR->GetBinContent(bin);
      double bc_LMpMR = d2LMpMR->GetBinContent(bin);
      double bc_LMMRp = d2LMMRp->GetBinContent(bin);

      double bc_LpMRR = d2LpMRR->GetBinContent(bin);
      double bc_LMpRR = d2LMpRR->GetBinContent(bin);
      double bc_LMRRp = d2LMRRp->GetBinContent(bin);

      double sum = bc_LpLMR/(be_LpLMR*be_LpLMR) + bc_LLMpR/(be_LLMpR*be_LLMpR) + bc_LLMRp/(be_LLMRp*be_LLMRp) +
      bc_LpMMR/(be_LpMMR*be_LpMMR) + bc_LMpMR/(be_LMpMR*be_LMpMR) + bc_LMMRp/(be_LMMRp*be_LMMRp) +
      bc_LpMRR/(be_LpMRR*be_LpMRR) + bc_LMpRR/(be_LMpRR*be_LMpRR) + bc_LMRRp/(be_LMRRp*be_LMRRp);
      double sumw = 1/(be_LpLMR*be_LpLMR) + 1/(be_LLMpR*be_LLMpR) + 1/(be_LLMRp*be_LLMRp) +
      1/(be_LpMMR*be_LpMMR) + 1/(be_LMpMR*be_LMpMR) + 1/(be_LMMRp*be_LMMRp) +
      1/(be_LpMRR*be_LpMRR) + 1/(be_LMpRR*be_LMpRR) + 1/(be_LMRRp*be_LMRRp);
      double avg = sum/sumw;
      reth->SetBinContent(bin,avg);
      reth->SetBinError(bin,std::sqrt(1/sumw));
  } */

  TH1D* reth = (TH1D*)hD2LpLMR->Clone("vpt24_3sub_avgOnVn");
  reth->Sumw2();
  reth->Add(hD2LLMpR);
  reth->Add(hD2LLMRp);
  reth->Add(hD2LpMMR);
  reth->Add(hD2LMpMR);
  reth->Add(hD2LMMRp);
  reth->Add(hD2LpMRR);
  reth->Add(hD2LMpRR);
  reth->Add(hD2LMRpR);
  reth->Scale(1.0 / 9);
  return reth;
}
TH1D* Calculator::getC24_3sub(int i)
{
  TH1D* c2llmr = getC2LLMR(i);
  TH1D* c2lmmr = getC2LMMR(i);
  TH1D* c2lmrr = getC2LMRR(i);

  TH1D* reth = (TH1D*)c2llmr->Clone("c24_3sub");
  for (int bin = 1; bin <= reth->GetNbinsX(); ++bin) {
    double be_llmr = c2llmr->GetBinError(bin);
    double be_lmmr = c2lmmr->GetBinError(bin);
    double be_lmrr = c2lmrr->GetBinError(bin);
    if (be_llmr == 0. || be_lmmr == 0. || be_lmrr == 0.) {
      reth->SetBinContent(bin, -2.);
      reth->SetBinError(bin, 0.);
    }

    double bc_llmr = c2llmr->GetBinContent(bin);
    double bc_lmmr = c2lmmr->GetBinContent(bin);
    double bc_lmrr = c2lmrr->GetBinContent(bin);

    double sum = bc_llmr / (be_llmr * be_llmr) + bc_lmmr / (be_lmmr * be_lmmr) + bc_lmrr / (be_lmrr * be_lmrr);
    double sumw = 1 / (be_llmr * be_llmr) + 1 / (be_lmmr * be_lmmr) + 1 / (be_lmrr * be_lmrr);
    double avg = sum / sumw;
    reth->SetBinContent(bin, avg);
    reth->SetBinError(bin, std::sqrt(1 / sumw));
  }
  return reth;
}
TH1D* Calculator::getD24_3sub(int i, double min, double max)
{
  TH1D* d2LpLMR = getD2LpLMR(i, min, max);
  TH1D* d2LLMpR = getD2LLMpR(i, min, max);
  TH1D* d2LLMRp = getD2LLMRp(i, min, max);

  TH1D* d2LpMMR = getD2LpMMR(i, min, max);
  TH1D* d2LMpMR = getD2LMpMR(i, min, max);
  TH1D* d2LMMRp = getD2LMMRp(i, min, max);

  TH1D* d2LpMRR = getD2LpMRR(i, min, max);
  TH1D* d2LMpRR = getD2LMpRR(i, min, max);
  TH1D* d2LMRRp = getD2LMRpR(i, min, max);

  TH1D* reth = (TH1D*)d2LpLMR->Clone("dpt24_3sub");
  for (int bin = 1; bin <= reth->GetNbinsX(); ++bin) {
    double be_LpLMR = d2LpLMR->GetBinError(bin);
    double be_LLMpR = d2LLMpR->GetBinError(bin);
    double be_LLMRp = d2LLMRp->GetBinError(bin);

    double be_LpMMR = d2LpMMR->GetBinError(bin);
    double be_LMpMR = d2LMpMR->GetBinError(bin);
    double be_LMMRp = d2LMMRp->GetBinError(bin);

    double be_LpMRR = d2LpMRR->GetBinError(bin);
    double be_LMpRR = d2LMpRR->GetBinError(bin);
    double be_LMRRp = d2LMRRp->GetBinError(bin);

    if (be_LpLMR == 0. || be_LLMpR == 0. || be_LLMRp == 0. || be_LpMMR == 0. || be_LMpMR == 0. || be_LMMRp == 0. || be_LpMRR == 0. || be_LMpRR == 0. || be_LMRRp == 0.) {
      reth->SetBinContent(bin, -2.);
      reth->SetBinError(bin, 0.);
    }

    double bc_LpLMR = d2LpLMR->GetBinContent(bin);
    double bc_LLMpR = d2LLMpR->GetBinContent(bin);
    double bc_LLMRp = d2LLMRp->GetBinContent(bin);

    double bc_LpMMR = d2LpMMR->GetBinContent(bin);
    double bc_LMpMR = d2LMpMR->GetBinContent(bin);
    double bc_LMMRp = d2LMMRp->GetBinContent(bin);

    double bc_LpMRR = d2LpMRR->GetBinContent(bin);
    double bc_LMpRR = d2LMpRR->GetBinContent(bin);
    double bc_LMRRp = d2LMRRp->GetBinContent(bin);

    double sum = bc_LpLMR / (be_LpLMR * be_LpLMR) + bc_LLMpR / (be_LLMpR * be_LLMpR) + bc_LLMRp / (be_LLMRp * be_LLMRp) +
                 bc_LpMMR / (be_LpMMR * be_LpMMR) + bc_LMpMR / (be_LMpMR * be_LMpMR) + bc_LMMRp / (be_LMMRp * be_LMMRp) +
                 bc_LpMRR / (be_LpMRR * be_LpMRR) + bc_LMpRR / (be_LMpRR * be_LMpRR) + bc_LMRRp / (be_LMRRp * be_LMRRp);
    double sumw = 1 / (be_LpLMR * be_LpLMR) + 1 / (be_LLMpR * be_LLMpR) + 1 / (be_LLMRp * be_LLMRp) +
                  1 / (be_LpMMR * be_LpMMR) + 1 / (be_LMpMR * be_LMpMR) + 1 / (be_LMMRp * be_LMMRp) +
                  1 / (be_LpMRR * be_LpMRR) + 1 / (be_LMpRR * be_LMpRR) + 1 / (be_LMRRp * be_LMRRp);
    double avg = sum / sumw;
    reth->SetBinContent(bin, avg);
    reth->SetBinError(bin, std::sqrt(1 / sumw));
  }
  return reth;
}

TH1D* Calculator::getV2LLMR(int i)
{
  TH1D* reth = getC2LLMR(i);
  reth->Scale(-1.);
  reth = powerHistogram(reth, 0.25);
  return reth;
}
TH1D* Calculator::getV2LMMR(int i)
{
  TH1D* reth = getC2LMMR(i);
  reth->Scale(-1.);
  reth = powerHistogram(reth, 0.25);
  return reth;
}
TH1D* Calculator::getV2LMRR(int i)
{
  TH1D* reth = getC2LMRR(i);
  reth->Scale(-1.);
  reth = powerHistogram(reth, 0.25);
  return reth;
}
TH1D* Calculator::getC2LLMR(int i)
{
  TH1D* reth = nullptr;
  auto clone = (i < 0) ? flowcont : dynamic_cast<FlowContainer*>(flowcont->Clone("TempClone"));
  if (i >= 0) {
    clone->OverrideMainWithSub(i, kFALSE);
    Int_t nxb;
    Double_t* xbs = flowcont->GetMultiRebin(nxb);
    clone->SetMultiRebin(nxb, xbs);
    delete[] xbs;
  };
  clone->SetIDName("LLMR");
  reth = dynamic_cast<TH1D*>(clone->GetHistCorrXXVsMulti("24", 0));
  clone->SetIDName("LM");
  TH1D* c2LM = dynamic_cast<TH1D*>(clone->GetHistCorrXXVsMulti("22", 0));
  clone->SetIDName("ChGap");
  TH1D* c2LR = dynamic_cast<TH1D*>(clone->GetHistCorrXXVsMulti("22", 0));
  c2LM->Multiply(c2LR);
  reth->Add(c2LM, -2);
  return reth;
}
TH1D* Calculator::getC2LMMR(int i)
{
  TH1D* reth = nullptr;
  auto clone = (i < 0) ? flowcont : dynamic_cast<FlowContainer*>(flowcont->Clone("TempClone"));
  if (i >= 0) {
    clone->OverrideMainWithSub(i, kFALSE);
    Int_t nxb;
    Double_t* xbs = flowcont->GetMultiRebin(nxb);
    clone->SetMultiRebin(nxb, xbs);
    delete[] xbs;
  };
  clone->SetIDName("LMMR");
  reth = dynamic_cast<TH1D*>(clone->GetHistCorrXXVsMulti("24", 0));
  clone->SetIDName("LM");
  TH1D* c2LM = dynamic_cast<TH1D*>(clone->GetHistCorrXXVsMulti("22", 0));
  clone->SetIDName("MR");
  TH1D* c2MR = dynamic_cast<TH1D*>(clone->GetHistCorrXXVsMulti("22", 0));
  c2LM->Multiply(c2MR);
  reth->Add(c2LM, -2);
  return reth;
}
TH1D* Calculator::getC2LMRR(int i)
{
  TH1D* reth = nullptr;
  auto clone = (i < 0) ? flowcont : dynamic_cast<FlowContainer*>(flowcont->Clone("TempClone"));
  if (i >= 0) {
    clone->OverrideMainWithSub(i, kFALSE);
    Int_t nxb;
    Double_t* xbs = flowcont->GetMultiRebin(nxb);
    clone->SetMultiRebin(nxb, xbs);
    delete[] xbs;
  };
  clone->SetIDName("LMRR");
  reth = dynamic_cast<TH1D*>(clone->GetHistCorrXXVsMulti("24", 0));
  clone->SetIDName("ChGap");
  TH1D* c2LR = dynamic_cast<TH1D*>(clone->GetHistCorrXXVsMulti("22", 0));
  clone->SetIDName("MR");
  TH1D* c2MR = dynamic_cast<TH1D*>(clone->GetHistCorrXXVsMulti("22", 0));
  c2LR->Multiply(c2MR);
  reth->Add(c2LR, -2);
  return reth;
}

TH1D* Calculator::getD2LpLMR(int i, double min, double max)
{
  TH1D* reth = nullptr;
  auto clone = (i < 0) ? flowcont : dynamic_cast<FlowContainer*>(flowcont->Clone("TempClone"));
  if (i >= 0) {
    clone->OverrideMainWithSub(i, kFALSE);
    Int_t nxb;
    Double_t* xbs = flowcont->GetMultiRebin(nxb);
    clone->SetMultiRebin(nxb, xbs);
    delete[] xbs;
  };
  clone->SetIDName("L'LMR");
  reth = dynamic_cast<TH1D*>(clone->GetHistCorrXXVsPt("24", min, max));
  clone->SetIDName("L'M");
  TH1D* c2LpM = dynamic_cast<TH1D*>(clone->GetHistCorrXXVsPt("22", min, max));
  clone->SetIDName("L'R");
  TH1D* c2LR = dynamic_cast<TH1D*>(clone->GetHistCorrXXVsMulti("22", 0));
  double c2LR_val = getIntegratedValue(c2LR, min, max);
  c2LpM->Scale(c2LR_val);
  reth->Add(c2LpM, -1);

  clone->SetIDName("L'M");
  TH1D* c2LM = dynamic_cast<TH1D*>(clone->GetHistCorrXXVsMulti("22", 0));
  clone->SetIDName("L'R");
  TH1D* c2LpR = dynamic_cast<TH1D*>(clone->GetHistCorrXXVsPt("22", min, max));
  double c2LM_val = getIntegratedValue(c2LM, min, max);
  c2LpR->Scale(c2LM_val);
  reth->Add(c2LpR, -1);

  return reth;
}
TH1D* Calculator::getD2LLMpR(int i, double min, double max)
{
  TH1D* reth = nullptr;
  auto clone = (i < 0) ? flowcont : dynamic_cast<FlowContainer*>(flowcont->Clone("TempClone"));
  if (i >= 0) {
    clone->OverrideMainWithSub(i, kFALSE);
    Int_t nxb;
    Double_t* xbs = flowcont->GetMultiRebin(nxb);
    clone->SetMultiRebin(nxb, xbs);
    delete[] xbs;
  };
  clone->SetIDName("LLM'R");
  reth = dynamic_cast<TH1D*>(clone->GetHistCorrXXVsPt("24", min, max));
  clone->SetIDName("LM'");
  TH1D* c2LMp = dynamic_cast<TH1D*>(clone->GetHistCorrXXVsPt("22", min, max));
  clone->SetIDName("L'R");
  TH1D* c2LR = dynamic_cast<TH1D*>(clone->GetHistCorrXXVsMulti("22", 0));
  double c2LR_val = getIntegratedValue(c2LR, min, max);
  c2LMp->Scale(c2LR_val);
  reth->Add(c2LMp, -2);

  return reth;
}
TH1D* Calculator::getD2LLMRp(int i, double min, double max)
{
  TH1D* reth = nullptr;
  auto clone = (i < 0) ? flowcont : dynamic_cast<FlowContainer*>(flowcont->Clone("TempClone"));
  if (i >= 0) {
    clone->OverrideMainWithSub(i, kFALSE);
    Int_t nxb;
    Double_t* xbs = flowcont->GetMultiRebin(nxb);
    clone->SetMultiRebin(nxb, xbs);
    delete[] xbs;
  };
  clone->SetIDName("LLMR'");
  reth = dynamic_cast<TH1D*>(clone->GetHistCorrXXVsPt("24", min, max));
  clone->SetIDName("L'M");
  TH1D* c2LM = dynamic_cast<TH1D*>(clone->GetHistCorrXXVsMulti("22", 0));
  double c2LM_val = getIntegratedValue(c2LM, min, max);
  clone->SetIDName("LR'");
  TH1D* c2LRp = dynamic_cast<TH1D*>(clone->GetHistCorrXXVsPt("22", min, max));
  c2LRp->Scale(c2LM_val);
  reth->Add(c2LRp, -2);
  return reth;
}

TH1D* Calculator::getD2LpMMR(int i, double min, double max)
{
  TH1D* reth = nullptr;
  auto clone = (i < 0) ? flowcont : dynamic_cast<FlowContainer*>(flowcont->Clone("TempClone"));
  if (i >= 0) {
    clone->OverrideMainWithSub(i, kFALSE);
    Int_t nxb;
    Double_t* xbs = flowcont->GetMultiRebin(nxb);
    clone->SetMultiRebin(nxb, xbs);
    delete[] xbs;
  };
  clone->SetIDName("L'MMR");
  reth = dynamic_cast<TH1D*>(clone->GetHistCorrXXVsPt("24", min, max));
  clone->SetIDName("L'M");
  TH1D* c2LM = dynamic_cast<TH1D*>(clone->GetHistCorrXXVsPt("22", min, max));
  clone->SetIDName("M'R");
  TH1D* c2LR = dynamic_cast<TH1D*>(clone->GetHistCorrXXVsMulti("22", 0));
  double c2LR_val = getIntegratedValue(c2LR, min, max);
  c2LM->Scale(c2LR_val);
  reth->Add(c2LM, -2);
  return reth;
}
TH1D* Calculator::getD2LMpMR(int i, double min, double max)
{
  TH1D* reth = nullptr;
  auto clone = (i < 0) ? flowcont : dynamic_cast<FlowContainer*>(flowcont->Clone("TempClone"));
  if (i >= 0) {
    clone->OverrideMainWithSub(i, kFALSE);
    Int_t nxb;
    Double_t* xbs = flowcont->GetMultiRebin(nxb);
    clone->SetMultiRebin(nxb, xbs);
    delete[] xbs;
  };
  clone->SetIDName("LM'MR");
  reth = dynamic_cast<TH1D*>(clone->GetHistCorrXXVsPt("24", min, max));
  clone->SetIDName("LM'");
  TH1D* c2LMp = dynamic_cast<TH1D*>(clone->GetHistCorrXXVsPt("22", min, max));
  clone->SetIDName("M'R");
  TH1D* c2MR = dynamic_cast<TH1D*>(clone->GetHistCorrXXVsMulti("22", 0));
  double c2MR_val = getIntegratedValue(c2MR, min, max);
  c2LMp->Scale(c2MR_val);
  reth->Add(c2LMp, -1);

  clone->SetIDName("L'M");
  TH1D* c2LM = dynamic_cast<TH1D*>(clone->GetHistCorrXXVsMulti("22", 0));
  double c2LM_val = getIntegratedValue(c2LM, min, max);
  clone->SetIDName("M'R");
  TH1D* c2MpR = dynamic_cast<TH1D*>(clone->GetHistCorrXXVsPt("22", min, max));
  c2MpR->Scale(c2LM_val);
  reth->Add(c2MpR, -1);

  return reth;
}
TH1D* Calculator::getD2LMMRp(int i, double min, double max)
{
  TH1D* reth = nullptr;
  auto clone = (i < 0) ? flowcont : dynamic_cast<FlowContainer*>(flowcont->Clone("TempClone"));
  if (i >= 0) {
    clone->OverrideMainWithSub(i, kFALSE);
    Int_t nxb;
    Double_t* xbs = flowcont->GetMultiRebin(nxb);
    clone->SetMultiRebin(nxb, xbs);
    delete[] xbs;
  };
  clone->SetIDName("LMMR'");
  reth = dynamic_cast<TH1D*>(clone->GetHistCorrXXVsPt("24", min, max));
  clone->SetIDName("L'M");
  TH1D* c2LM = dynamic_cast<TH1D*>(clone->GetHistCorrXXVsMulti("22", 0));
  double c2LM_val = getIntegratedValue(c2LM, min, max);
  clone->SetIDName("MR'");
  TH1D* c2MRp = dynamic_cast<TH1D*>(clone->GetHistCorrXXVsPt("22", min, max));
  c2MRp->Scale(c2LM_val);
  reth->Add(c2MRp, -2);
  return reth;
}

TH1D* Calculator::getD2LpMRR(int i, double min, double max)
{
  TH1D* reth = nullptr;
  auto clone = (i < 0) ? flowcont : dynamic_cast<FlowContainer*>(flowcont->Clone("TempClone"));
  if (i >= 0) {
    clone->OverrideMainWithSub(i, kFALSE);
    Int_t nxb;
    Double_t* xbs = flowcont->GetMultiRebin(nxb);
    clone->SetMultiRebin(nxb, xbs);
    delete[] xbs;
  };
  clone->SetIDName("L'MRR");
  reth = dynamic_cast<TH1D*>(clone->GetHistCorrXXVsPt("24", min, max));
  clone->SetIDName("L'R");
  TH1D* c2LpR = dynamic_cast<TH1D*>(clone->GetHistCorrXXVsPt("22", min, max));
  clone->SetIDName("M'R");
  TH1D* c2MR = dynamic_cast<TH1D*>(clone->GetHistCorrXXVsMulti("22", 0));
  double c2MR_val = getIntegratedValue(c2MR, min, max);
  c2LpR->Scale(c2MR_val);
  reth->Add(c2LpR, -2);
  return reth;
}
TH1D* Calculator::getD2LMpRR(int i, double min, double max)
{
  TH1D* reth = nullptr;
  auto clone = (i < 0) ? flowcont : dynamic_cast<FlowContainer*>(flowcont->Clone("TempClone"));
  if (i >= 0) {
    clone->OverrideMainWithSub(i, kFALSE);
    Int_t nxb;
    Double_t* xbs = flowcont->GetMultiRebin(nxb);
    clone->SetMultiRebin(nxb, xbs);
    delete[] xbs;
  };
  clone->SetIDName("LM'RR");
  reth = dynamic_cast<TH1D*>(clone->GetHistCorrXXVsPt("24", min, max));
  clone->SetIDName("L'R");
  TH1D* c2LM = dynamic_cast<TH1D*>(clone->GetHistCorrXXVsMulti("22", 0));
  double c2LM_val = getIntegratedValue(c2LM, min, max);
  clone->SetIDName("M'R");
  TH1D* c2MpR = dynamic_cast<TH1D*>(clone->GetHistCorrXXVsPt("22", min, max));
  c2MpR->Scale(c2LM_val);
  reth->Add(c2MpR, -2);
  return reth;
}
TH1D* Calculator::getD2LMRpR(int i, double min, double max)
{
  TH1D* reth = nullptr;
  auto clone = (i < 0) ? flowcont : dynamic_cast<FlowContainer*>(flowcont->Clone("TempClone"));
  if (i >= 0) {
    clone->OverrideMainWithSub(i, kFALSE);
    Int_t nxb;
    Double_t* xbs = flowcont->GetMultiRebin(nxb);
    clone->SetMultiRebin(nxb, xbs);
    delete[] xbs;
  };
  clone->SetIDName("LMR'R");
  reth = dynamic_cast<TH1D*>(clone->GetHistCorrXXVsPt("24", min, max));
  clone->SetIDName("LR'");
  TH1D* c2LRp = dynamic_cast<TH1D*>(clone->GetHistCorrXXVsPt("22", min, max));
  clone->SetIDName("M'R");
  TH1D* c2MR = dynamic_cast<TH1D*>(clone->GetHistCorrXXVsMulti("22", 0));
  double c2MR_val = getIntegratedValue(c2MR, min, max);
  c2LRp->Scale(c2MR_val);
  reth->Add(c2LRp, -1);

  clone->SetIDName("L'R");
  TH1D* c2LR = dynamic_cast<TH1D*>(clone->GetHistCorrXXVsMulti("22", 0));
  double c2LR_val = getIntegratedValue(c2LR, min, max);
  clone->SetIDName("MR'");
  TH1D* c2MRp = dynamic_cast<TH1D*>(clone->GetHistCorrXXVsPt("22", min, max));
  c2MRp->Scale(c2LR_val);
  reth->Add(c2MRp, -1);
  return reth;
}

TH1D* Calculator::get2LpLMR(int i, double min, double max)
{
  TH1D* reth = nullptr;
  auto clone = (i < 0) ? flowcont : dynamic_cast<FlowContainer*>(flowcont->Clone("TempClone"));
  if (i >= 0) {
    clone->OverrideMainWithSub(i, kFALSE);
    Int_t nxb;
    Double_t* xbs = flowcont->GetMultiRebin(nxb);
    clone->SetMultiRebin(nxb, xbs);
    delete[] xbs;
  };
  clone->SetIDName("L'LMR");
  reth = dynamic_cast<TH1D*>(clone->GetHistCorrXXVsPt("24", min, max));
  return reth;
}
TH1D* Calculator::get2LLMpR(int i, double min, double max)
{
  TH1D* reth = nullptr;
  auto clone = (i < 0) ? flowcont : dynamic_cast<FlowContainer*>(flowcont->Clone("TempClone"));
  if (i >= 0) {
    clone->OverrideMainWithSub(i, kFALSE);
    Int_t nxb;
    Double_t* xbs = flowcont->GetMultiRebin(nxb);
    clone->SetMultiRebin(nxb, xbs);
    delete[] xbs;
  };
  clone->SetIDName("LLM'R");
  reth = dynamic_cast<TH1D*>(clone->GetHistCorrXXVsPt("24", min, max));
  return reth;
}
TH1D* Calculator::get2LLMRp(int i, double min, double max)
{
  TH1D* reth = nullptr;
  auto clone = (i < 0) ? flowcont : dynamic_cast<FlowContainer*>(flowcont->Clone("TempClone"));
  if (i >= 0) {
    clone->OverrideMainWithSub(i, kFALSE);
    Int_t nxb;
    Double_t* xbs = flowcont->GetMultiRebin(nxb);
    clone->SetMultiRebin(nxb, xbs);
    delete[] xbs;
  };
  clone->SetIDName("LLMR'");
  reth = dynamic_cast<TH1D*>(clone->GetHistCorrXXVsPt("24", min, max));
  return reth;
}

TH1D* Calculator::get2LpMMR(int i, double min, double max)
{
  TH1D* reth = nullptr;
  auto clone = (i < 0) ? flowcont : dynamic_cast<FlowContainer*>(flowcont->Clone("TempClone"));
  if (i >= 0) {
    clone->OverrideMainWithSub(i, kFALSE);
    Int_t nxb;
    Double_t* xbs = flowcont->GetMultiRebin(nxb);
    clone->SetMultiRebin(nxb, xbs);
    delete[] xbs;
  };
  clone->SetIDName("L'MMR");
  reth = dynamic_cast<TH1D*>(clone->GetHistCorrXXVsPt("24", min, max));
  return reth;
}
TH1D* Calculator::get2LMpMR(int i, double min, double max)
{
  TH1D* reth = nullptr;
  auto clone = (i < 0) ? flowcont : dynamic_cast<FlowContainer*>(flowcont->Clone("TempClone"));
  if (i >= 0) {
    clone->OverrideMainWithSub(i, kFALSE);
    Int_t nxb;
    Double_t* xbs = flowcont->GetMultiRebin(nxb);
    clone->SetMultiRebin(nxb, xbs);
    delete[] xbs;
  };
  clone->SetIDName("LM'MR");
  reth = dynamic_cast<TH1D*>(clone->GetHistCorrXXVsPt("24", min, max));
  return reth;
}
TH1D* Calculator::get2LMMRp(int i, double min, double max)
{
  TH1D* reth = nullptr;
  auto clone = (i < 0) ? flowcont : dynamic_cast<FlowContainer*>(flowcont->Clone("TempClone"));
  if (i >= 0) {
    clone->OverrideMainWithSub(i, kFALSE);
    Int_t nxb;
    Double_t* xbs = flowcont->GetMultiRebin(nxb);
    clone->SetMultiRebin(nxb, xbs);
    delete[] xbs;
  };
  clone->SetIDName("LMMR'");
  reth = dynamic_cast<TH1D*>(clone->GetHistCorrXXVsPt("24", min, max));
  return reth;
}

TH1D* Calculator::get2LpMRR(int i, double min, double max)
{
  TH1D* reth = nullptr;
  auto clone = (i < 0) ? flowcont : dynamic_cast<FlowContainer*>(flowcont->Clone("TempClone"));
  if (i >= 0) {
    clone->OverrideMainWithSub(i, kFALSE);
    Int_t nxb;
    Double_t* xbs = flowcont->GetMultiRebin(nxb);
    clone->SetMultiRebin(nxb, xbs);
    delete[] xbs;
  };
  clone->SetIDName("L'MRR");
  reth = dynamic_cast<TH1D*>(clone->GetHistCorrXXVsPt("24", min, max));
  return reth;
}
TH1D* Calculator::get2LMpRR(int i, double min, double max)
{
  TH1D* reth = nullptr;
  auto clone = (i < 0) ? flowcont : dynamic_cast<FlowContainer*>(flowcont->Clone("TempClone"));
  if (i >= 0) {
    clone->OverrideMainWithSub(i, kFALSE);
    Int_t nxb;
    Double_t* xbs = flowcont->GetMultiRebin(nxb);
    clone->SetMultiRebin(nxb, xbs);
    delete[] xbs;
  };
  clone->SetIDName("LM'RR");
  reth = dynamic_cast<TH1D*>(clone->GetHistCorrXXVsPt("24", min, max));
  return reth;
}
TH1D* Calculator::get2LMRpR(int i, double min, double max)
{
  TH1D* reth = nullptr;
  auto clone = (i < 0) ? flowcont : dynamic_cast<FlowContainer*>(flowcont->Clone("TempClone"));
  if (i >= 0) {
    clone->OverrideMainWithSub(i, kFALSE);
    Int_t nxb;
    Double_t* xbs = flowcont->GetMultiRebin(nxb);
    clone->SetMultiRebin(nxb, xbs);
    delete[] xbs;
  };
  clone->SetIDName("LMR'R");
  reth = dynamic_cast<TH1D*>(clone->GetHistCorrXXVsPt("24", min, max));
  return reth;
}

TH1D* Calculator::get2LRp(int i, double min, double max)
{
  TH1D* reth = nullptr;
  auto clone = (i < 0) ? flowcont : dynamic_cast<FlowContainer*>(flowcont->Clone("TempClone"));
  if (i >= 0) {
    clone->OverrideMainWithSub(i, kFALSE);
    Int_t nxb;
    Double_t* xbs = flowcont->GetMultiRebin(nxb);
    clone->SetMultiRebin(nxb, xbs);
    delete[] xbs;
  };
  clone->SetIDName("ChPGap");
  reth = dynamic_cast<TH1D*>(clone->GetHistCorrXXVsPt("22", min, max));
  return reth;
}
TH1D* Calculator::get2LpR(int i, double min, double max)
{
  TH1D* reth = nullptr;
  auto clone = (i < 0) ? flowcont : dynamic_cast<FlowContainer*>(flowcont->Clone("TempClone"));
  if (i >= 0) {
    clone->OverrideMainWithSub(i, kFALSE);
    Int_t nxb;
    Double_t* xbs = flowcont->GetMultiRebin(nxb);
    clone->SetMultiRebin(nxb, xbs);
    delete[] xbs;
  };
  clone->SetIDName("ChNGap");
  reth = dynamic_cast<TH1D*>(clone->GetHistCorrXXVsPt("22", min, max));
  return reth;
}
TH1D* Calculator::get2LpM(int i, double min, double max)
{
  TH1D* reth = nullptr;
  auto clone = (i < 0) ? flowcont : dynamic_cast<FlowContainer*>(flowcont->Clone("TempClone"));
  if (i >= 0) {
    clone->OverrideMainWithSub(i, kFALSE);
    Int_t nxb;
    Double_t* xbs = flowcont->GetMultiRebin(nxb);
    clone->SetMultiRebin(nxb, xbs);
    delete[] xbs;
  };
  clone->SetIDName("L'M");
  reth = dynamic_cast<TH1D*>(clone->GetHistCorrXXVsPt("22", min, max));
  return reth;
}
TH1D* Calculator::get2LMp(int i, double min, double max)
{
  TH1D* reth = nullptr;
  auto clone = (i < 0) ? flowcont : dynamic_cast<FlowContainer*>(flowcont->Clone("TempClone"));
  if (i >= 0) {
    clone->OverrideMainWithSub(i, kFALSE);
    Int_t nxb;
    Double_t* xbs = flowcont->GetMultiRebin(nxb);
    clone->SetMultiRebin(nxb, xbs);
    delete[] xbs;
  };
  clone->SetIDName("LM'");
  reth = dynamic_cast<TH1D*>(clone->GetHistCorrXXVsPt("22", min, max));
  return reth;
}
TH1D* Calculator::get2MpR(int i, double min, double max)
{
  TH1D* reth = nullptr;
  auto clone = (i < 0) ? flowcont : dynamic_cast<FlowContainer*>(flowcont->Clone("TempClone"));
  if (i >= 0) {
    clone->OverrideMainWithSub(i, kFALSE);
    Int_t nxb;
    Double_t* xbs = flowcont->GetMultiRebin(nxb);
    clone->SetMultiRebin(nxb, xbs);
    delete[] xbs;
  };
  clone->SetIDName("M'R");
  reth = dynamic_cast<TH1D*>(clone->GetHistCorrXXVsPt("22", min, max));
  return reth;
}
TH1D* Calculator::get2MRp(int i, double min, double max)
{
  TH1D* reth = nullptr;
  auto clone = (i < 0) ? flowcont : dynamic_cast<FlowContainer*>(flowcont->Clone("TempClone"));
  if (i >= 0) {
    clone->OverrideMainWithSub(i, kFALSE);
    Int_t nxb;
    Double_t* xbs = flowcont->GetMultiRebin(nxb);
    clone->SetMultiRebin(nxb, xbs);
    delete[] xbs;
  };
  clone->SetIDName("MR'");
  reth = dynamic_cast<TH1D*>(clone->GetHistCorrXXVsPt("22", min, max));
  return reth;
}

TH1D* Calculator::get2LR(int i)
{
  TH1D* reth = nullptr;
  auto clone = (i < 0) ? flowcont : dynamic_cast<FlowContainer*>(flowcont->Clone("TempClone"));
  if (i >= 0) {
    clone->OverrideMainWithSub(i, kFALSE);
    Int_t nxb;
    Double_t* xbs = flowcont->GetMultiRebin(nxb);
    clone->SetMultiRebin(nxb, xbs);
    delete[] xbs;
  };
  clone->SetIDName("ChGap");
  reth = dynamic_cast<TH1D*>(clone->GetHistCorrXXVsMulti("22", 0));
  return reth;
}
TH1D* Calculator::get2LM(int i)
{
  TH1D* reth = nullptr;
  auto clone = (i < 0) ? flowcont : dynamic_cast<FlowContainer*>(flowcont->Clone("TempClone"));
  if (i >= 0) {
    clone->OverrideMainWithSub(i, kFALSE);
    Int_t nxb;
    Double_t* xbs = flowcont->GetMultiRebin(nxb);
    clone->SetMultiRebin(nxb, xbs);
    delete[] xbs;
  };
  clone->SetIDName("LM");
  reth = dynamic_cast<TH1D*>(clone->GetHistCorrXXVsMulti("22", 0));
  return reth;
}
TH1D* Calculator::get2MR(int i)
{
  TH1D* reth = nullptr;
  auto clone = (i < 0) ? flowcont : dynamic_cast<FlowContainer*>(flowcont->Clone("TempClone"));
  if (i >= 0) {
    clone->OverrideMainWithSub(i, kFALSE);
    Int_t nxb;
    Double_t* xbs = flowcont->GetMultiRebin(nxb);
    clone->SetMultiRebin(nxb, xbs);
    delete[] xbs;
  };
  clone->SetIDName("MR");
  reth = dynamic_cast<TH1D*>(clone->GetHistCorrXXVsMulti("22", 0));
  return reth;
}

TH1D* Calculator::getCVN2pt(int nhar, int i)
{
  if (flowptcont->usesCentralMoments()) {
    auto vn2mpt_profile = dynamic_cast<BootstrapProfile*>(corr_list->FindObject(Form("%s%i2pt1_Mpt0", flowid.c_str(), nhar)));
    auto vn2_mpt_profile = dynamic_cast<BootstrapProfile*>(corr_list->FindObject(Form("%s%i2pt1_Mpt1", flowid.c_str(), nhar)));
    auto mpt = getMpt(i);
    TH1D* reth = dynamic_cast<TH1D*>(vn2mpt_profile->getHist(i));
    TH1D* vn2_mpt = dynamic_cast<TH1D*>(vn2_mpt_profile->getHist(i));
    vn2_mpt->Multiply(mpt);
    reth->Add(vn2_mpt, -1);
    return reth;
  }
  auto vn2mpt_profile = dynamic_cast<BootstrapProfile*>(corr_list->FindObject(Form("%s%i2_mpt1", flowid.c_str(), nhar)));
  TH1D* reth = dynamic_cast<TH1D*>(vn2mpt_profile->getHist(i));
  auto vn2 = getNN(nhar, 2, flowid.c_str(), i);
  auto mpt = getMpt(i);
  vn2->Multiply(mpt);
  reth->Add(vn2, -1);
  return reth;
}
TH1D* Calculator::getRhoVN2pt(int nhar, int i)
{
  TH1D* reth = getCVN2pt(nhar, i);
  TH1D* varv22 = getVarV22(i);
  TH1D* ck = getMu2(i);
  TH1D* sqrtck = powerHistogram(ck, 0.5);
  TH1D* sqrtvarv22 = powerHistogram(varv22, 0.5);
  sqrtck->Multiply(sqrtvarv22);
  reth->Divide(sqrtck);
  return reth;
}
TH1D* Calculator::getCV22pt2(int i)
{
  if (flowptcont->usesCentralMoments()) {
    auto v22pt2_Mpt0_profile = dynamic_cast<BootstrapProfile*>(corr_list->FindObject("ChGap22pt2_Mpt0"));
    // First get the <v2^2*deltapt^2> term
    TH1D* reth = dynamic_cast<TH1D*>(v22pt2_Mpt0_profile->getHist(i));
    auto mpt_hist = getMpt(i);
    auto v22pt2_Mpt1_profile = dynamic_cast<BootstrapProfile*>(corr_list->FindObject("ChGap22pt2_Mpt1"));
    auto v22pt_pt = dynamic_cast<TH1D*>(v22pt2_Mpt1_profile->getHist(i));
    v22pt_pt->Multiply(mpt_hist);
    reth->Add(v22pt_pt);
    auto v22_hist = getNN(2, 2, "ChGap", i);
    auto mpt_sq_hist = dynamic_cast<TH1D*>(mpt_hist->Clone("mpt_sq_hist"));
    mpt_sq_hist->Multiply(mpt_hist);
    auto v22_mptsq_hist = dynamic_cast<TH1D*>(v22_hist->Clone("v22_mptsq_hist"));
    v22_mptsq_hist->Multiply(mpt_sq_hist);
    reth->Add(v22_mptsq_hist);
    // Then subtract <v2^2><deltapt^2>
    auto v22_cm2_hist = dynamic_cast<TH1D*>(v22_hist->Clone("v22_cm2_hist"));
    v22_cm2_hist->Multiply(flowptcont->getCentralMomentHist(i, 2));
    reth->Add(v22_cm2_hist, -1);
    return reth;
  }
  return nullptr;
}
TH1D* Calculator::getCV22pt3(int i)
{
  if (flowptcont->usesCentralMoments()) {
    // First get the <v2^2*deltapt^3> term
    auto v22pt3_Mpt0_profile = dynamic_cast<BootstrapProfile*>(corr_list->FindObject("ChGap22pt3_Mpt0"));
    auto v22pt3_Mpt1_profile = dynamic_cast<BootstrapProfile*>(corr_list->FindObject("ChGap22pt3_Mpt1"));
    auto v22pt3_Mpt2_profile = dynamic_cast<BootstrapProfile*>(corr_list->FindObject("ChGap22pt3_Mpt2"));
    TH1D* reth = dynamic_cast<TH1D*>(v22pt3_Mpt0_profile->getHist(i));
    auto mpt_hist = getMpt(i);
    auto v22pt2_pt = dynamic_cast<TH1D*>(v22pt3_Mpt1_profile->getHist(i));
    v22pt2_pt->Multiply(mpt_hist);
    mpt_hist->Multiply(mpt_hist);
    auto v22pt_pt2 = dynamic_cast<TH1D*>(v22pt3_Mpt2_profile->getHist(i));
    v22pt_pt2->Multiply(mpt_hist);
    mpt_hist->Multiply(mpt_hist);
    auto v22_pt3 = getNN(2, 2, "ChGap", i);
    v22_pt3->Multiply(mpt_hist);
    reth->Add(v22pt2_pt, -1);
    reth->Add(v22pt_pt2, -1);
    reth->Add(v22_pt3, -1);
    // Then subtract the 3 x <v2^2*deltapt>*<deltapt^2> term
    auto v22deltapt_deltapt2 = getCV22pt(i);
    v22deltapt_deltapt2->Multiply(flowptcont->getCentralMomentHist(i, 2));
    reth->Add(v22deltapt_deltapt2, -3);
    // And then subtract the <v2^2>*<deltapt^3> term
    auto v22_deltapt3 = getNN(2, 2, "ChGap", i);
    v22_deltapt3->Multiply(flowptcont->getCentralMomentHist(i, 3));
    reth->Add(v22_deltapt3, -1);
    return reth;
  }
  return nullptr;
}
TH1D* Calculator::getCV22pt4(int i)
{
  if (flowptcont->usesCentralMoments()) {
    // First get the <v2^2*deltapt^4> term
    auto v22pt4_Mpt0_profile = dynamic_cast<BootstrapProfile*>(corr_list->FindObject("ChGap22pt4_Mpt0"));
    auto v22pt4_Mpt1_profile = dynamic_cast<BootstrapProfile*>(corr_list->FindObject("ChGap22pt4_Mpt1"));
    auto v22pt4_Mpt2_profile = dynamic_cast<BootstrapProfile*>(corr_list->FindObject("ChGap22pt4_Mpt2"));
    auto v22pt4_Mpt3_profile = dynamic_cast<BootstrapProfile*>(corr_list->FindObject("ChGap22pt4_Mpt3"));
    TH1D* reth = dynamic_cast<TH1D*>(v22pt4_Mpt0_profile->getHist(i));
    auto mpt_hist = getMpt(i);
    auto v22pt3_pt = dynamic_cast<TH1D*>(v22pt4_Mpt1_profile->getHist(i));
    v22pt3_pt->Multiply(mpt_hist);
    mpt_hist->Multiply(mpt_hist);
    auto v22pt2_pt2 = dynamic_cast<TH1D*>(v22pt4_Mpt2_profile->getHist(i));
    v22pt2_pt2->Multiply(mpt_hist);
    auto v22pt_pt3 = dynamic_cast<TH1D*>(v22pt4_Mpt3_profile->getHist(i));
    v22pt_pt3->Multiply(mpt_hist);
    mpt_hist->Multiply(mpt_hist);
    auto v22_pt4 = getNN(2, 2, "ChGap", i);
    v22_pt4->Multiply(mpt_hist);
    reth->Add(v22pt3_pt, -1);
    reth->Add(v22pt2_pt2, -1);
    reth->Add(v22pt_pt3, -1);
    reth->Add(v22_pt4, -1);
    // Then subtract the 6 x <v2^2*deltapt2>*<deltapt^2> term
    auto v22deltapt2_deltapt2 = getCV22pt2(i);
    v22deltapt2_deltapt2->Multiply(flowptcont->getCentralMomentHist(i, 2));
    reth->Add(v22deltapt2_deltapt2, -6);
    // Then subtract the 4 x <v2^2*deltapt>*<deltapt^3> term
    auto v22deltapt_deltapt3 = getCV22pt(i);
    v22deltapt_deltapt3->Multiply(flowptcont->getCentralMomentHist(i, 3));
    reth->Add(v22deltapt_deltapt3, -4);
    // Then add the 6 x <v2^2>*<deltapt^2>*<deltapt^2> term
    auto v22_deltapt2_deltapt2 = getNN(2, 2, "ChGap", i);
    v22_deltapt2_deltapt2->Multiply(flowptcont->getCentralMomentHist(i, 2));
    v22_deltapt2_deltapt2->Multiply(flowptcont->getCentralMomentHist(i, 2));
    reth->Add(v22_deltapt2_deltapt2, -4);
    // Then subtract the <v2^2>*<deltapt^4> term
    auto v22_deltapt4 = getNN(2, 2, "ChGap", i);
    v22_deltapt4->Multiply(flowptcont->getCentralMomentHist(i, 4));
    reth->Add(v22_deltapt4, 1);
    return reth;
  }
  return nullptr;
}
TH1D* Calculator::getCV24pt(int i)
{
  if (flowptcont->usesCentralMoments()) {
    // First get the <v2^4*deltapt> term
    BootstrapProfile* v24pt_Mpt0_profile = dynamic_cast<BootstrapProfile*>(corr_list->FindObject("ChGap24pt_Mpt0"));
    TH1D* reth = dynamic_cast<TH1D*>(v24pt_Mpt0_profile->getHist(i));
    auto v24_pt = getNN(2, 4, "ChGap", i);
    auto mpt_hist = getMpt(i);
    v24_pt->Multiply(mpt_hist);
    reth->Add(v24_pt, -1);
    // Then subtract the 4 x <v2^2*deltapt>*<v2^2> term
    auto v22pt_v22 = getCV22pt(i);
    v22pt_v22->Multiply(getNN(2, 2, "ChGap", i));
    reth->Add(v22pt_v22, -4);
    return reth;
  }
  return nullptr;
}
TH1D* Calculator::getCV24pt2(int i)
{
  if (flowptcont->usesCentralMoments()) {
    // First get the <v2^4*deltapt^2> term
    auto v24pt2_profile = dynamic_cast<BootstrapProfile*>(corr_list->FindObject("ChFull24pt2_Mpt0"));
    auto v24pt_6w_profile = dynamic_cast<BootstrapProfile*>(corr_list->FindObject("ChFull24pt2_Mpt1"));
    auto v24_6w_profile = dynamic_cast<BootstrapProfile*>(corr_list->FindObject("ChFull24pt2_Mpt2"));
    auto reth = dynamic_cast<TH1D*>(v24pt2_profile->getHist(i));
    auto v24pt_mpt = dynamic_cast<TH1D*>(v24pt_6w_profile->getHist(i));
    auto mpt = getMpt(i);
    v24pt_mpt->Multiply(mpt);
    auto v24_mpt_mpt = dynamic_cast<TH1D*>(v24_6w_profile->getHist(i));
    v24_mpt_mpt->Multiply(mpt);
    v24_mpt_mpt->Multiply(mpt);
    reth->Add(v24pt_mpt, -2);
    reth->Add(v24_mpt_mpt);

    // Then subtract the 4 x <v2^2>*<vn^2*deltapt^2> term
    auto v22pt2_profile = dynamic_cast<BootstrapProfile*>(corr_list->FindObject("ChFull22pt2_Mpt0"));
    auto v22pt_4w_profile = dynamic_cast<BootstrapProfile*>(corr_list->FindObject("ChFull22pt2_Mpt1"));
    auto v22_4w_profile = dynamic_cast<BootstrapProfile*>(corr_list->FindObject("ChFull22pt2_Mpt2"));
    auto v22deltapt2_v22 = dynamic_cast<TH1D*>(v22pt2_profile->getHist(i));
    auto v22pt_mpt = dynamic_cast<TH1D*>(v22pt_4w_profile->getHist(i));
    v22pt_mpt->Multiply(mpt);
    auto v22_mpt_mpt = dynamic_cast<TH1D*>(v22_4w_profile->getHist(i));
    v22_mpt_mpt->Multiply(mpt);
    v22_mpt_mpt->Multiply(mpt);
    v22deltapt2_v22->Add(v22pt_mpt, -2);
    v22deltapt2_v22->Add(v22_mpt_mpt);
    v22deltapt2_v22->Multiply(getNN(2, 2, "ChFull", i));
    reth->Add(v22deltapt2_v22, -4);

    // Then subtract the 4 x <v2^2*deltapt>*<v2^2*deltapt> term
    auto v22pt_profile = dynamic_cast<BootstrapProfile*>(corr_list->FindObject("ChFull22pt1_Mpt0"));
    auto v22_3w_profile = dynamic_cast<BootstrapProfile*>(corr_list->FindObject("ChFull22pt1_Mpt1"));
    auto v22deltapt_v22deltapt = dynamic_cast<TH1D*>(v22pt_profile->getHist(i));
    auto v22_mpt = dynamic_cast<TH1D*>(v22_3w_profile->getHist(i));
    v22_mpt->Multiply(mpt);
    v22deltapt_v22deltapt->Add(v22_mpt, -1);
    v22deltapt_v22deltapt->Multiply(v22deltapt_v22deltapt);
    reth->Add(v22deltapt_v22deltapt, -4);

    // Then add the 4 x <v2^2>*<v2^2>*<deltapt^2> term
    auto v22_v22_deltapt2 = getNN(2, 2, "ChFull", i);
    v22_v22_deltapt2->Multiply(v22_v22_deltapt2);
    v22_v22_deltapt2->Multiply(flowptcont->getCentralMomentHist(i, 2));
    reth->Add(v22_v22_deltapt2, 4);

    // Then subtract the <v2^4>*<deltapt^2> term
    auto v24_deltapt2 = getNN(2, 4, "ChFull", i);
    v24_deltapt2->Multiply(flowptcont->getCentralMomentHist(i, 2));
    reth->Add(v24_deltapt2, -1);

    return reth;
  }
  return nullptr;
}
// Radial flow analyis//
TH1D* Calculator::getProjectedFractionalSpectra(std::string sp, int pidIndex, double min, double max, FractionScheme scheme, bool combinedWeight)
{
  auto& analysis = scheme == FractionScheme::V02 ? (combinedWeight ? radialFlowV02Weighted : radialFlowV02)
                                                 : (combinedWeight ? radialFlowV0Weighted : radialFlowV0);
  return pidIndex >= 0 && pidIndex < static_cast<int>(analysis.fraction_cent.size())
           ? projectFraction(analysis.fraction_cent[pidIndex], sp, min, max)
           : nullptr;
}
TH1D* Calculator::getProjectedFractionalK0Spectra(std::string sp, int resoIndex, double min, double max, FractionScheme scheme, bool combinedWeight)
{
  auto& analysis = scheme == FractionScheme::V02 ? (combinedWeight ? radialFlowV02Weighted : radialFlowV02)
                                                 : (combinedWeight ? radialFlowV0Weighted : radialFlowV0);
  return resoIndex >= 0 && resoIndex < static_cast<int>(analysis.fraction_cent_K0.size())
           ? projectFraction(analysis.fraction_cent_K0[resoIndex], sp, min, max)
           : nullptr;
}
TH1D* Calculator::getProjectedFractionalLambdaSpectra(std::string sp, int resoIndex, double min, double max, FractionScheme scheme, bool combinedWeight)
{
  auto& analysis = scheme == FractionScheme::V02 ? (combinedWeight ? radialFlowV02Weighted : radialFlowV02)
                                                 : (combinedWeight ? radialFlowV0Weighted : radialFlowV0);
  return resoIndex >= 0 && resoIndex < static_cast<int>(analysis.fraction_cent_Lambda.size())
           ? projectFraction(analysis.fraction_cent_Lambda[resoIndex], sp, min, max)
           : nullptr;
}

TH1D* Calculator::getnpt(int i, int pidIndex, double min, double max)
{
  return getnpt(i, pidIndex, min, max, FractionScheme::V02);
}
TH1D* Calculator::getnpt(int i, int pidIndex, double min, double max, FractionScheme scheme, bool combinedWeight)
{
  std::vector<std::string> species = {"Ch", "Pi", "Ka", "Pr"};
  if (pidIndex > 3) {
    std::cout << "[Calculator::getV0]: Invalid pid index - must be 0 (ch), 1 (pi), 2 (ka), or 3 (pr)" << std::endl;
    return nullptr;
  }
  std::string sp = species[pidIndex];
  TH1D* mean_fraction_pt = getProjectedFractionalSpectra(sp, pidIndex, min, max, scheme, combinedWeight);
  if (!mean_fraction_pt)
    return nullptr;
  TH1D* reth = (TH1D*)mean_fraction_pt->Clone(Form("mpt_%snpt_%.1f_%.1f", sp.c_str(), min, max));
  reth->SetName(Form("npt%s_%i_%.1f_%.1f", sp.c_str(), i, min, max));
  return reth;
}
TH1D* Calculator::getnptK0(int i, double min, double max)
{
  return getnptK0(i, min, max, FractionScheme::V02);
}
TH1D* Calculator::getnptK0(int i, double min, double max, FractionScheme scheme, bool combinedWeight)
{
  TH1D* mean_fraction_pt_tot = getProjectedFractionalK0Spectra("K0_sig", 1, min, max, scheme, combinedWeight);
  if (!mean_fraction_pt_tot)
    return nullptr;

  TH1D* npt_tot = (TH1D*)mean_fraction_pt_tot->Clone(Form("npt_K0_tot_%.1f_%.1f", min, max));
  if (useTotal) {
    npt_tot->SetName(Form("nptK0_%i", i));
    return npt_tot;
  }
  TH1D* mean_fraction_pt_sb1 = getProjectedFractionalK0Spectra("K0_sb1", 0, min, max, scheme, combinedWeight);
  TH1D* mean_fraction_pt_sb2 = getProjectedFractionalK0Spectra("K0_sb2", 2, min, max, scheme, combinedWeight);
  if (!mean_fraction_pt_sb1 || !mean_fraction_pt_sb2) {
    delete npt_tot;
    return nullptr;
  }
  TH1D* npt_sb1 = (TH1D*)mean_fraction_pt_sb1->Clone(Form("npt_K0_sb1_%.1f_%.1f", min, max));
  TH1D* npt_sb2 = (TH1D*)mean_fraction_pt_sb2->Clone(Form("npt_K0_sb2_%.1f_%.1f", min, max));

  TH1D* npt_bkg = (TH1D*)npt_sb1->Clone(Form("v0_K0_bkg_%.1f_%.1f", min, max));
  npt_bkg->Add(npt_sb2, 1);
  npt_bkg->Scale(0.5);

  TH1D* npt_sig = (TH1D*)npt_tot->Clone(Form("npt_K0_sig_%.1f_%.1f", min, max));
  if (!useTotal)
    applySignalFractionCorrection(*this, npt_sig, npt_bkg, "K0", min, max, "Calculator::getnptK0");
  npt_sig->SetName(Form("nptK0_%i", i));
  return npt_sig;
}
TH1D* Calculator::getnptLambda(int i, double min, double max)
{
  return getnptLambda(i, min, max, FractionScheme::V02);
}
TH1D* Calculator::getnptLambda(int i, double min, double max, FractionScheme scheme, bool combinedWeight)
{
  TH1D* mean_fraction_pt_tot = getProjectedFractionalLambdaSpectra("Lambda_sig", 1, min, max, scheme, combinedWeight);
  if (!mean_fraction_pt_tot)
    return nullptr;

  TH1D* npt_tot = (TH1D*)mean_fraction_pt_tot->Clone(Form("npt_Lambda_tot_%.1f_%.1f", min, max));
  if (useTotal) {
    npt_tot->SetName(Form("nptLambda_%i", i));
    return npt_tot;
  }
  TH1D* mean_fraction_pt_sb1 = getProjectedFractionalLambdaSpectra("Lambda_sb1", 0, min, max, scheme, combinedWeight);
  TH1D* mean_fraction_pt_sb2 = getProjectedFractionalLambdaSpectra("Lambda_sb2", 2, min, max, scheme, combinedWeight);
  if (!mean_fraction_pt_sb1 || !mean_fraction_pt_sb2) {
    delete npt_tot;
    return nullptr;
  }
  TH1D* npt_sb1 = (TH1D*)mean_fraction_pt_sb1->Clone(Form("npt_Lambda_sb1_%.1f_%.1f", min, max));
  TH1D* npt_sb2 = (TH1D*)mean_fraction_pt_sb2->Clone(Form("npt_Lambda_sb2_%.1f_%.1f", min, max));

  TH1D* npt_bkg = (TH1D*)npt_sb1->Clone(Form("npt_Lambda_bkg_%.1f_%.1f", min, max));
  npt_bkg->Add(npt_sb2, 1);
  npt_bkg->Scale(0.5);

  TH1D* npt_sig = (TH1D*)npt_tot->Clone(Form("npt_Lambda_sig_%.1f_%.1f", min, max));
  if (!useTotal)
    applySignalFractionCorrection(*this, npt_sig, npt_bkg, "Lambda", min, max, "Calculator::getnptLambda");
  npt_sig->SetName(Form("nptLambda_%i", i));
  return npt_sig;
}

TH1D* Calculator::getV0Corr(int i, int pidIndex, double min, double max)
{
  std::vector<std::string> species = {"Ch", "Pi", "Ka", "Pr"};
  if (pidIndex > 3) {
    std::cout << "[Calculator::getV0Corr]: Invalid pid index - must be 0 (ch), 1 (pi), 2 (ka), or 3 (pr)" << std::endl;
    return nullptr;
  }
  std::string sp = species[pidIndex];
  TH1D* reth = getCorrVsPt(2, 2, Form("ChMpt%s", sp.c_str()), i, min, max);
  reth->SetName(Form("V0Corr%s_%i_%.1f_%.1f", sp.c_str(), i, min, max));
  return reth;
}
TH1D* Calculator::getV0CorrK0(int i, double min, double max)
{
  TH1D* v0_tot = getCorrVsPt(2, 2, Form("ChMptK0s"), i, min, max);
  if (useTotal) {
    v0_tot->SetName(Form("V0CorrK0_%i", i));
    return v0_tot;
  }
  TH1D* v0_sb1 = getCorrVsPt(2, 2, Form("ChMptK0sSBL"), i, min, max);
  TH1D* v0_sb2 = getCorrVsPt(2, 2, Form("ChMptK0sSBR"), i, min, max);

  TH1D* v0_bkg = (TH1D*)v0_sb1->Clone(Form("v0corr_K0_bkg_%.1f_%.1f", min, max));
  v0_bkg->Add(v0_sb2, 1);
  v0_bkg->Scale(0.5);

  TH1D* v0_sig = (TH1D*)v0_tot->Clone(Form("v0corr_K0_sig_%.1f_%.1f", min, max));
  if (!useTotal)
    applySignalFractionCorrection(*this, v0_sig, v0_bkg, "K0", min, max, "Calculator::getV0CorrK0");
  v0_sig->SetName(Form("V0CorrK0_%i", i));
  return v0_sig;
}
TH1D* Calculator::getV0CorrLambda(int i, double min, double max)
{
  TH1D* v0_tot = getCorrVsPt(2, 2, Form("ChMptLambda"), i, min, max);
  if (useTotal) {
    v0_tot->SetName(Form("V0CorrLambda_%i", i));
    return v0_tot;
  }
  TH1D* v0_sb1 = getCorrVsPt(2, 2, Form("ChMptLambdaSBL"), i, min, max);
  TH1D* v0_sb2 = getCorrVsPt(2, 2, Form("ChMptLambdaSBR"), i, min, max);

  TH1D* v0_bkg = (TH1D*)v0_sb1->Clone(Form("v0corr_Lambda_bkg_%.1f_%.1f", min, max));
  v0_bkg->Add(v0_sb2, 1);
  v0_bkg->Scale(0.5);

  TH1D* v0_sig = (TH1D*)v0_tot->Clone(Form("v0corr_Lambda_sig_%.1f_%.1f", min, max));
  if (!useTotal)
    applySignalFractionCorrection(*this, v0_sig, v0_bkg, "Lambda", min, max, "Calculator::getV0CorrLambda");
  v0_sig->SetName(Form("V0CorrLambda_%i", i));
  return v0_sig;
}
TH1D* Calculator::getV0Cov(int i, int pidIndex, double min, double max)
{
  std::vector<std::string> species = {"Ch", "Pi", "Ka", "Pr"};
  if (pidIndex > 3) {
    std::cout << "[Calculator::getV0Cov]: Invalid pid index - must be 0 (ch), 1 (pi), 2 (ka), or 3 (pr)" << std::endl;
    return nullptr;
  }
  std::string sp = species[pidIndex];
  TH1D* reth = getCorrVsPt(2, 2, Form("ChMpt%s", sp.c_str()), i, min, max);
  TH1D* sub = getV0Means(i, pidIndex, min, max);
  reth->Add(sub, -1);
  reth->SetName(Form("V0Cov%s_%i_%.1f_%.1f", sp.c_str(), i, min, max));
  return reth;
}
TH1D* Calculator::getV0CovK0(int i, double min, double max)
{
  TH1D* reth = getV0CorrK0(i, min, max);
  TH1D* sub = getV0MeansK0(i, min, max);
  reth->Add(sub, -1);
  reth->SetName(Form("V0CovK0_%i", i));
  return reth;
}
TH1D* Calculator::getV0CovLambda(int i, double min, double max)
{
  TH1D* reth = getV0CorrLambda(i, min, max);
  TH1D* sub = getV0MeansLambda(i, min, max);
  reth->Add(sub, -1);
  reth->SetName(Form("V0CovLambda_%i", i));
  return reth;
}

TH1D* Calculator::getV0Means(int i, int pidIndex, double min, double max)
{
  std::vector<std::string> species = {"Ch", "Pi", "Ka", "Pr"};
  if (pidIndex > 3) {
    std::cout << "[Calculator::getV0Means]: Invalid pid index - must be 0 (ch), 1 (pi), 2 (ka), or 3 (pr)" << std::endl;
    return nullptr;
  }
  std::string sp = species[pidIndex];
  TH1D* mean_fraction_pt = getProjectedFractionalSpectra(sp, pidIndex, min, max, FractionScheme::V0);
  TH1D* mpt = getMpt(i);
  TH1D* reth = (TH1D*)mean_fraction_pt->Clone(Form("mpt_%snpt_%.1f_%.1f", sp.c_str(), min, max));
  reth->Scale(getIntegratedValue(mpt, min, max));
  reth->SetName(Form("V0Means%s_%i_%.1f_%.1f", sp.c_str(), i, min, max));
  return reth;
}
TH1D* Calculator::getV0MeansK0(int i, double min, double max)
{
  TH1D* mean_fraction_pt_tot = getProjectedFractionalK0Spectra("K0_sig", 1, min, max, FractionScheme::V0);
  TH1D* mpt = getMpt(i);
  TH1D* mpt_npt_tot = (TH1D*)mean_fraction_pt_tot->Clone(Form("mpt_npt_K0_tot_%.1f_%.1f", min, max));
  mpt_npt_tot->Scale(getIntegratedValue(mpt, min, max));
  if (useTotal) {
    mpt_npt_tot->SetName(Form("V0MeansK0_%i", i));
    return mpt_npt_tot;
  }
  TH1D* mean_fraction_pt_sb1 = getProjectedFractionalK0Spectra("K0_sb1", 0, min, max, FractionScheme::V0);
  TH1D* mean_fraction_pt_sb2 = getProjectedFractionalK0Spectra("K0_sb2", 2, min, max, FractionScheme::V0);
  TH1D* mpt_npt_sb1 = (TH1D*)mean_fraction_pt_sb1->Clone(Form("mpt_npt_K0_sb1_%.1f_%.1f", min, max));
  mpt_npt_sb1->Scale(getIntegratedValue(mpt, min, max));
  TH1D* mpt_npt_sb2 = (TH1D*)mean_fraction_pt_sb2->Clone(Form("mpt_npt_K0_sb2_%.1f_%.1f", min, max));
  mpt_npt_sb2->Scale(getIntegratedValue(mpt, min, max));

  TH1D* mpt_npt_bkg = (TH1D*)mpt_npt_sb1->Clone(Form("mpt_npt_K0_bkg_%.1f_%.1f", min, max));
  mpt_npt_bkg->Add(mpt_npt_sb2, 1);
  mpt_npt_bkg->Scale(0.5);

  TH1D* mpt_npt_sig = (TH1D*)mpt_npt_tot->Clone(Form("mpt_npt_K0_sig_%.1f_%.1f", min, max));
  if (!useTotal)
    applySignalFractionCorrection(*this, mpt_npt_sig, mpt_npt_bkg, "K0", min, max, "Calculator::getV0MeansK0");
  mpt_npt_sig->SetName(Form("V0MeansK0_%i", i));
  return mpt_npt_sig;
}
TH1D* Calculator::getV0MeansLambda(int i, double min, double max)
{
  TH1D* mean_fraction_pt_tot = getProjectedFractionalLambdaSpectra("Lambda_sig", 1, min, max, FractionScheme::V0);
  TH1D* mpt = getMpt(i);
  TH1D* mpt_npt_tot = (TH1D*)mean_fraction_pt_tot->Clone(Form("mpt_npt_Lambda_tot_%.1f_%.1f", min, max));
  mpt_npt_tot->Scale(getIntegratedValue(mpt, min, max));
  if (useTotal) {
    mpt_npt_tot->SetName(Form("V0MeansLambda_%i", i));
    return mpt_npt_tot;
  }
  TH1D* mean_fraction_pt_sb1 = getProjectedFractionalLambdaSpectra("Lambda_sb1", 0, min, max, FractionScheme::V0);
  TH1D* mean_fraction_pt_sb2 = getProjectedFractionalLambdaSpectra("Lambda_sb2", 2, min, max, FractionScheme::V0);
  TH1D* mpt_npt_sb1 = (TH1D*)mean_fraction_pt_sb1->Clone(Form("mpt_npt_Lambda_sb1_%.1f_%.1f", min, max));
  mpt_npt_sb1->Scale(getIntegratedValue(mpt, min, max));
  TH1D* mpt_npt_sb2 = (TH1D*)mean_fraction_pt_sb2->Clone(Form("mpt_npt_Lambda_sb2_%.1f_%.1f", min, max));
  mpt_npt_sb2->Scale(getIntegratedValue(mpt, min, max));

  TH1D* mpt_npt_bkg = (TH1D*)mpt_npt_sb1->Clone(Form("mpt_npt_Lambda_bkg_%.1f_%.1f", min, max));
  mpt_npt_bkg->Add(mpt_npt_sb2, 1);
  mpt_npt_bkg->Scale(0.5);

  TH1D* mpt_npt_sig = (TH1D*)mpt_npt_tot->Clone(Form("mpt_npt_Lambda_sig_%.1f_%.1f", min, max));
  if (!useTotal)
    applySignalFractionCorrection(*this, mpt_npt_sig, mpt_npt_bkg, "Lambda", min, max, "Calculator::getV0MeansLambda");
  mpt_npt_sig->SetName(Form("V0MeansLambda_%i", i));
  return mpt_npt_sig;
}
TH1D* Calculator::getV0(int i, int pidIndex, double min, double max)
{
  std::vector<std::string> species = {"Ch", "Pi", "Ka", "Pr"};
  if (pidIndex > 3) {
    std::cout << "[Calculator::getV0]: Invalid pid index - must be 0 (ch), 1 (pi), 2 (ka), or 3 (pr)" << std::endl;
    return nullptr;
  }
  std::string sp = species[pidIndex];
  TH1D* reth = getV0Corr(i, pidIndex, min, max);
  TH1D* npt = getnpt(i, pidIndex, min, max, FractionScheme::V0);
  TH1D* mpt = getMpt(i);

  TH1D* mpt_npt = (TH1D*)npt->Clone(Form("mpt_%snpt_%.1f_%.1f", sp.c_str(), min, max));
  mpt_npt->Scale(getIntegratedValue(mpt, min, max));

  TH1D* mu2 = getMu2sub(i);
  mu2 = powerHistogram(mu2, 0.5);
  TH1D* mu_npt = (TH1D*)npt->Clone(Form("mu_npt_%.1f_%.1f", min, max));
  mu_npt->Scale(getIntegratedValue(mu2, min, max));

  reth->Add(mpt_npt, -1);
  reth->Divide(mu_npt);
  reth->SetName(Form("V0%s_%i_%.1f_%.1f", sp.c_str(), i, min, max));
  return reth;
}
TH1D* Calculator::getFractionForSpecies(int i, int pidIndex, double min, double max, FractionScheme scheme, bool combinedWeight)
{
  // Fraction projections are temporary; keep them out of the output TFile.
  TDirectoryFile scratchDirectory(Form("weightedFractionScratch_%d_%d", i, pidIndex), "", "", gROOT);
  TDirectory::TContext directoryContext(&scratchDirectory);
  TH1D* result = nullptr;
  if (pidIndex >= 0 && pidIndex < ChargedSpeciesCount)
    result = getnpt(i, pidIndex, min, max, scheme, combinedWeight);
  else if (pidIndex == K0SpeciesIndex)
    result = getnptK0(i, min, max, scheme, combinedWeight);
  else if (pidIndex == LambdaSpeciesIndex)
    result = getnptLambda(i, min, max, scheme, combinedWeight);
  else
    std::cerr << "[Calculator::getFractionForSpecies] Invalid species index " << pidIndex << std::endl;
  if (result)
    result->SetDirectory(nullptr);
  return result;
}

TH1D* Calculator::getCorrelationForSpecies(int i, int pidIndex, double min, double max, FractionScheme scheme)
{
  if (scheme == FractionScheme::V02) {
    if (pidIndex >= 0 && pidIndex < ChargedSpeciesCount)
      return getV02Corr(i, pidIndex, min, max);
    if (pidIndex == K0SpeciesIndex)
      return getV02CorrK0(i, min, max);
    if (pidIndex == LambdaSpeciesIndex)
      return getV02CorrLambda(i, min, max);
  } else {
    if (pidIndex >= 0 && pidIndex < ChargedSpeciesCount)
      return getV0Corr(i, pidIndex, min, max);
    if (pidIndex == K0SpeciesIndex)
      return getV0CorrK0(i, min, max);
    if (pidIndex == LambdaSpeciesIndex)
      return getV0CorrLambda(i, min, max);
  }
  std::cerr << "[Calculator::getCorrelationForSpecies] Invalid species index " << pidIndex << std::endl;
  return nullptr;
}

TH1D* Calculator::getV0NumeratorWeighted(int i, int pidIndex, double min, double max)
{
  // <mpt*n>_AB - <mpt>_AB<n>_B - <mpt>_A<n>_AB + <mpt>_A<n>_B.
  static constexpr const char* ProfileNames[] = {"mpt1_w2pc_ch", "mpt1_w2pc_pi", "mpt1_w2pc_ka", "mpt1_w2pc_pr"};
  if (pidIndex < 0 || pidIndex > LambdaSpeciesIndex || !flowptcont || !flowptcont->getCorrList())
    return nullptr;
  // K0 and Lambda use the charged-track normalization when cfgUsePIDTotal is false.
  const char* profileName = ProfileNames[pidIndex < ChargedSpeciesCount ? pidIndex : 0];
  auto* profile = dynamic_cast<BootstrapProfile*>(flowptcont->getCorrList()->FindObject(profileName));
  if (!profile) {
    std::cerr << "[Calculator::getV0NumeratorWeighted] Missing " << profileName << std::endl;
    return nullptr;
  }
  std::unique_ptr<TH1D> mA(dynamic_cast<TH1D*>(flowptcont->getCorrHist(i, 1)));
  std::unique_ptr<TH1D> mAB(dynamic_cast<TH1D*>(profile->getHist(i)));
  std::unique_ptr<TH1D> nB(getFractionForSpecies(i, pidIndex, min, max, FractionScheme::V0, false));
  std::unique_ptr<TH1D> nAB(getFractionForSpecies(i, pidIndex, min, max, FractionScheme::V0, true));
  if (!mA || !mAB || !nB || !nAB) {
    std::cerr << "[Calculator::getV0NumeratorWeighted] Missing V0 fraction or mean-pT input for species " << pidIndex << std::endl;
    return nullptr;
  }
  std::unique_ptr<TH1D> result(getCorrelationForSpecies(i, pidIndex, min, max, FractionScheme::V0));
  if (!result)
    return nullptr;
  const double meanA = getIntegratedValue(mA.get(), min, max);
  const double meanAB = getIntegratedValue(mAB.get(), min, max);
  if (!std::isfinite(meanA) || !std::isfinite(meanAB))
    return nullptr;
  result->Add(nB.get(), meanA - meanAB);
  result->Add(nAB.get(), -meanA);
  result->SetName(Form("V0NumeratorWeighted_%d_%d_%.1f_%.1f", pidIndex, i, min, max));
  return result.release();
}

TH1D* Calculator::getV02NumeratorWeighted(int i, int pidIndex, double min, double max)
{
  // <v2^2*n>_ABC - <v2^2>_ABC<n>_C - <v2^2>_AB<n>_ABC + <v2^2>_AB<n>_C.
  if (pidIndex < 0 || pidIndex > LambdaSpeciesIndex || !hasChGap3pcW22) {
    std::cerr << "[Calculator::getV02NumeratorWeighted] Missing ChGap3pcW22 or invalid species index" << std::endl;
    return nullptr;
  }
  std::unique_ptr<TH1D> vAB(getCNN(2, 2, "ChGap", i));
  std::unique_ptr<TH1D> vABC(getCNN(2, 2, "ChGap3pcW", i));
  std::unique_ptr<TH1D> nC(getFractionForSpecies(i, pidIndex, min, max, FractionScheme::V02, false));
  std::unique_ptr<TH1D> nABC(getFractionForSpecies(i, pidIndex, min, max, FractionScheme::V02, true));
  if (!vAB || !vABC || !nC || !nABC) {
    std::cerr << "[Calculator::getV02NumeratorWeighted] Missing V02 flow or fraction input for species " << pidIndex << std::endl;
    return nullptr;
  }
  std::unique_ptr<TH1D> result(getCorrelationForSpecies(i, pidIndex, min, max, FractionScheme::V02));
  if (!result)
    return nullptr;
  const double meanAB = getIntegratedValue(vAB.get(), min, max);
  const double meanABC = getIntegratedValue(vABC.get(), min, max);
  if (!std::isfinite(meanAB) || !std::isfinite(meanABC))
    return nullptr;
  result->Add(nC.get(), meanAB - meanABC);
  result->Add(nABC.get(), -meanAB);
  result->SetName(Form("V02NumeratorWeighted_%d_%d_%.1f_%.1f", pidIndex, i, min, max));
  return result.release();
}

TH1D* Calculator::getV0Weighted(int i, int pidIndex, double min, double max)
{
  std::unique_ptr<TH1D> result(getV0NumeratorWeighted(i, pidIndex, min, max));
  std::unique_ptr<TH1D> nB(getFractionForSpecies(i, pidIndex, min, max, FractionScheme::V0, false));
  std::unique_ptr<TH1D> mu2(getMu2sub(i));
  if (!result || !nB || !mu2)
    return nullptr;
  std::unique_ptr<TH1D> sigmaPt(powerHistogram(mu2.get(), 0.5));
  if (!sigmaPt)
    return nullptr;
  const double denominatorScale = getIntegratedValue(sigmaPt.get(), min, max);
  if (!std::isfinite(denominatorScale))
    return nullptr;
  nB->Scale(denominatorScale);
  result->Divide(nB.get());
  result->SetName(Form("V0Weighted_%d_%d_%.1f_%.1f", pidIndex, i, min, max));
  result->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  return result.release();
}

TH1D* Calculator::getV02Weighted(int i, int pidIndex, double min, double max)
{
  std::unique_ptr<TH1D> result(getV02NumeratorWeighted(i, pidIndex, min, max));
  std::unique_ptr<TH1D> nC(getFractionForSpecies(i, pidIndex, min, max, FractionScheme::V02, false));
  std::unique_ptr<TH1D> vAB(getCNN(2, 2, "ChGap", i));
  if (!result || !nC || !vAB)
    return nullptr;
  const double denominatorScale = getIntegratedValue(vAB.get(), min, max);
  if (!std::isfinite(denominatorScale))
    return nullptr;
  nC->Scale(denominatorScale);
  result->Divide(nC.get());
  result->SetName(Form("V02Weighted_%d_%d_%.1f_%.1f", pidIndex, i, min, max));
  result->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  return result.release();
}

// Alternatives using the final V0
TH1D* Calculator::getV0K0Alt(int i, double min, double max)
{
  TH1D* mean_fraction_pt_sb1 = getProjectedFractionalK0Spectra("K0_sb1", 0, min, max, FractionScheme::V0);
  TH1D* mean_fraction_pt_tot = getProjectedFractionalK0Spectra("K0_sig", 1, min, max, FractionScheme::V0);
  TH1D* mean_fraction_pt_sb2 = getProjectedFractionalK0Spectra("K0_sb2", 2, min, max, FractionScheme::V0);
  TH1D* v0_tot = getCorrVsPt(2, 2, Form("ChMptK0s"), i, min, max);
  TH1D* v0_sb1 = getCorrVsPt(2, 2, Form("ChMptK0sSBL"), i, min, max);
  TH1D* v0_sb2 = getCorrVsPt(2, 2, Form("ChMptK0sSBR"), i, min, max);

  TList* l = (TList*)flowptcont->fSubList;
  TH1* pt1_sub1 = ((BootstrapProfile*)l->FindObject("ptptsub_1sub1_0sub2"))->getHist(i);
  TH1* pt1_sub2 = ((BootstrapProfile*)l->FindObject("ptptsub_0sub1_1sub2"))->getHist(i);
  TH1D* mpt = (TH1D*)pt1_sub1->Clone(Form("mpt_avg_%i", i));
  mpt->Add(pt1_sub2);
  mpt->Scale(0.5);
  TH1D* mpt_npt_tot = (TH1D*)mean_fraction_pt_tot->Clone(Form("v0_K0_tot_%.1f_%.1f", min, max));
  mpt_npt_tot->Scale(getIntegratedValue(mpt, min, max));
  TH1D* mpt_npt_sb1 = (TH1D*)mean_fraction_pt_sb1->Clone(Form("v0_K0_sb1_%.1f_%.1f", min, max));
  mpt_npt_sb1->Scale(getIntegratedValue(mpt, min, max));
  TH1D* mpt_npt_sb2 = (TH1D*)mean_fraction_pt_sb2->Clone(Form("v0_K0_sb2_%.1f_%.1f", min, max));
  mpt_npt_sb2->Scale(getIntegratedValue(mpt, min, max));

  TH1D* mu2 = getMu2sub(i);
  mu2 = powerHistogram(mu2, 0.5);
  TH1D* mu_npt_tot = (TH1D*)mean_fraction_pt_tot->Clone(Form("v0_K0_den_tot_%.1f_%.1f", min, max));
  mu_npt_tot->Scale(getIntegratedValue(mu2, min, max));
  TH1D* mu_npt_sb1 = (TH1D*)mean_fraction_pt_sb1->Clone(Form("v0_K0_den_sb1_%.1f_%.1f", min, max));
  mu_npt_sb1->Scale(getIntegratedValue(mu2, min, max));
  TH1D* mu_npt_sb2 = (TH1D*)mean_fraction_pt_sb2->Clone(Form("v0_K0_den_sb2_%.1f_%.1f", min, max));
  mu_npt_sb2->Scale(getIntegratedValue(mu2, min, max));

  v0_tot->Add(mpt_npt_tot, -1);
  v0_tot->Divide(mu_npt_tot);

  v0_sb1->Add(mpt_npt_sb1, -1);
  v0_sb1->Divide(mu_npt_sb1);
  v0_sb2->Add(mpt_npt_sb2, -1);
  v0_sb2->Divide(mu_npt_sb2);

  TH1D* v0_bkg = (TH1D*)v0_sb1->Clone(Form("v0_K0_bkg_%.1f_%.1f", min, max));
  v0_bkg->Add(v0_sb2, 1);
  v0_bkg->Scale(0.5);

  TH1D* v0_sig = (TH1D*)v0_tot->Clone(Form("v0_K0_sig_%.1f_%.1f", min, max));
  if (!useTotal)
    applySignalFractionCorrection(*this, v0_sig, v0_bkg, "K0", min, max, "Calculator::getV0K0");
  v0_sig->SetName(Form("V0K0_%i", i));
  return v0_sig;
}
TH1D* Calculator::getV0LambdaAlt(int i, double min, double max)
{
  TH1D* mean_fraction_pt_sb1 = getProjectedFractionalLambdaSpectra("Lambda_sb1", 0, min, max, FractionScheme::V0);
  TH1D* mean_fraction_pt_tot = getProjectedFractionalLambdaSpectra("Lambda_sig", 1, min, max, FractionScheme::V0);
  TH1D* mean_fraction_pt_sb2 = getProjectedFractionalLambdaSpectra("Lambda_sb2", 2, min, max, FractionScheme::V0);
  TH1D* v0_tot = getCorrVsPt(2, 2, Form("ChMptLambda"), i, min, max);
  TH1D* v0_sb1 = getCorrVsPt(2, 2, Form("ChMptLambdaSBL"), i, min, max);
  TH1D* v0_sb2 = getCorrVsPt(2, 2, Form("ChMptLambdaSBR"), i, min, max);

  TList* l = (TList*)flowptcont->fSubList;
  TH1* pt1_sub1 = ((BootstrapProfile*)l->FindObject("ptptsub_1sub1_0sub2"))->getHist(i);
  TH1* pt1_sub2 = ((BootstrapProfile*)l->FindObject("ptptsub_0sub1_1sub2"))->getHist(i);
  TH1D* mpt = (TH1D*)pt1_sub1->Clone(Form("mpt_avg_%i", i));
  mpt->Add(pt1_sub2);
  mpt->Scale(0.5);
  TH1D* mpt_npt_tot = (TH1D*)mean_fraction_pt_tot->Clone(Form("v0_Lambda_tot_%.1f_%.1f", min, max));
  mpt_npt_tot->Scale(getIntegratedValue(mpt, min, max));
  TH1D* mpt_npt_sb1 = (TH1D*)mean_fraction_pt_sb1->Clone(Form("v0_Lambda_sb1_%.1f_%.1f", min, max));
  mpt_npt_sb1->Scale(getIntegratedValue(mpt, min, max));
  TH1D* mpt_npt_sb2 = (TH1D*)mean_fraction_pt_sb2->Clone(Form("v0_Lambda_sb2_%.1f_%.1f", min, max));
  mpt_npt_sb2->Scale(getIntegratedValue(mpt, min, max));

  TH1D* mu2 = getMu2sub(i);
  mu2 = powerHistogram(mu2, 0.5);
  TH1D* mu_npt_tot = (TH1D*)mean_fraction_pt_tot->Clone(Form("v0_Lambda_den_tot_%.1f_%.1f", min, max));
  mu_npt_tot->Scale(getIntegratedValue(mu2, min, max));
  TH1D* mu_npt_sb1 = (TH1D*)mean_fraction_pt_sb1->Clone(Form("v0_Lambda_den_sb1_%.1f_%.1f", min, max));
  mu_npt_sb1->Scale(getIntegratedValue(mu2, min, max));
  TH1D* mu_npt_sb2 = (TH1D*)mean_fraction_pt_sb2->Clone(Form("v0_Lambda_den_sb2_%.1f_%.1f", min, max));
  mu_npt_sb2->Scale(getIntegratedValue(mu2, min, max));

  v0_tot->Add(mpt_npt_tot, -1);
  v0_tot->Divide(mu_npt_tot);

  v0_sb1->Add(mpt_npt_sb1, -1);
  v0_sb1->Divide(mu_npt_sb1);
  v0_sb2->Add(mpt_npt_sb2, -1);
  v0_sb2->Divide(mu_npt_sb2);

  TH1D* v0_bkg = (TH1D*)v0_sb1->Clone(Form("v0_Lambda_bkg_%.1f_%.1f", min, max));
  v0_bkg->Add(v0_sb2, 1);
  v0_bkg->Scale(0.5);

  TH1D* v0_sig = (TH1D*)v0_tot->Clone(Form("v0_Lambda_sig_%.1f_%.1f", min, max));
  if (!useTotal)
    applySignalFractionCorrection(*this, v0_sig, v0_bkg, "Lambda", min, max, "Calculator::getV0Lambda");
  v0_sig->SetName(Form("V0Lambda_%i", i));
  return v0_sig;
}

TH1D* Calculator::getV02Corr(int i, int pidIndex, double min, double max)
{
  std::vector<std::string> species = {"Ch", "Pi", "Ka", "Pr"};
  if (pidIndex > 3) {
    std::cout << "[Calculator::getV02Corr]: Invalid pid index - must be 0 (ch), 1 (pi), 2 (ka), or 3 (pr)" << std::endl;
    return nullptr;
  }
  std::string sp = species[pidIndex];
  TH1D* reth = getCorrVsPt(2, 2, Form("ChGap%s", sp.c_str()), i, min, max);
  reth->SetName(Form("V02Corr%s_%i_%.1f_%.1f", sp.c_str(), i, min, max));
  return reth;
}
TH1D* Calculator::getV02CorrK0(int i, double min, double max)
{
  TH1D* v02_tot = getCorrVsPt(2, 2, Form("ChGapK0s"), i, min, max);
  if (useTotal) {
    v02_tot->SetName(Form("V02CorrK0_%i", i));
    return v02_tot;
  }
  TH1D* v02_sb1 = getCorrVsPt(2, 2, Form("ChGapK0sSBL"), i, min, max);
  TH1D* v02_sb2 = getCorrVsPt(2, 2, Form("ChGapK0sSBR"), i, min, max);

  TH1D* v02_bkg = (TH1D*)v02_sb1->Clone(Form("v02corr_K0_bkg_%.1f_%.1f", min, max));
  v02_bkg->Add(v02_sb2, 1);
  v02_bkg->Scale(0.5);

  TH1D* v02_sig = (TH1D*)v02_tot->Clone(Form("v02_K0_sig_%.1f_%.1f", min, max));
  if (!useTotal)
    applySignalFractionCorrection(*this, v02_sig, v02_bkg, "K0", min, max, "Calculator::getV02CorrK0");
  v02_sig->SetName(Form("V02CorrK0_%i", i));
  return v02_sig;
}
TH1D* Calculator::getV02CorrLambda(int i, double min, double max)
{
  TH1D* v02_tot = getCorrVsPt(2, 2, Form("ChGapLambda"), i, min, max);
  if (useTotal) {
    v02_tot->SetName(Form("V02CorrLambda_%i", i));
    return v02_tot;
  }
  TH1D* v02_sb1 = getCorrVsPt(2, 2, Form("ChGapLambdaSBL"), i, min, max);
  TH1D* v02_sb2 = getCorrVsPt(2, 2, Form("ChGapLambdaSBR"), i, min, max);
  TH1D* v02_bkg = (TH1D*)v02_sb1->Clone(Form("v02corr_Lambda_bkg_%.1f_%.1f", min, max));
  v02_bkg->Add(v02_sb2, 1);
  v02_bkg->Scale(0.5);

  TH1D* v02_sig = (TH1D*)v02_tot->Clone(Form("v02corr_Lambda_sig_%.1f_%.1f", min, max));
  if (!useTotal)
    applySignalFractionCorrection(*this, v02_sig, v02_bkg, "Lambda", min, max, "Calculator::getV02Lambda");
  v02_sig->SetName(Form("V02CorrLambda_%i", i));
  return v02_sig;
}
TH1D* Calculator::getV02Cov(int i, int pidIndex, double min, double max)
{
  std::vector<std::string> species = {"Ch", "Pi", "Ka", "Pr"};
  if (pidIndex > 3) {
    std::cout << "[Calculator::getV02Cov]: Invalid pid index - must be 0 (ch), 1 (pi), 2 (ka), or 3 (pr)" << std::endl;
    return nullptr;
  }
  std::string sp = species[pidIndex];
  TH1D* reth = getCorrVsPt(2, 2, Form("ChGap%s", sp.c_str()), i, min, max);
  TH1D* npt = getnpt(i, pidIndex, min, max, FractionScheme::V02);
  TH1D* v22_3pc = hasChGap3pcW22 ? getCNN(2, 2, "ChGap3pcW", i) : getCNN(2, 2, "ChGap", i);
  TH1D* sub = (TH1D*)npt->Clone(Form("v22_3pc_%snpt_%.1f_%.1f", sp.c_str(), min, max));
  sub->Scale(getIntegratedValue(v22_3pc, min, max));
  reth->Add(sub, -1);
  reth->SetName(Form("V02Cov%s_%i_%.1f_%.1f", sp.c_str(), i, min, max));
  return reth;
}
TH1D* Calculator::getV02CovK0(int i, double min, double max)
{
  TH1D* reth = getV02CorrK0(i, min, max);
  TH1D* npt = getnptK0(i, min, max, FractionScheme::V02);
  TH1D* v22_3pc = hasChGap3pcW22 ? getCNN(2, 2, "ChGap3pcW", i) : getCNN(2, 2, "ChGap", i);
  TH1D* sub = (TH1D*)npt->Clone(Form("v22_3pc_K0npt_%.1f_%.1f", min, max));
  sub->Scale(getIntegratedValue(v22_3pc, min, max));
  reth->Add(sub, -1);
  reth->SetName(Form("V02CovK0_%i", i));
  return reth;
}
TH1D* Calculator::getV02CovLambda(int i, double min, double max)
{
  TH1D* reth = getV02CorrLambda(i, min, max);
  TH1D* npt = getnptLambda(i, min, max, FractionScheme::V02);
  TH1D* v22_3pc = hasChGap3pcW22 ? getCNN(2, 2, "ChGap3pcW", i) : getCNN(2, 2, "ChGap", i);
  TH1D* sub = (TH1D*)npt->Clone(Form("v22_3pc_Lambdanpt_%.1f_%.1f", min, max));
  sub->Scale(getIntegratedValue(v22_3pc, min, max));
  reth->Add(sub, -1);
  reth->SetName(Form("V02CovLambda_%i", i));
  return reth;
}
TH1D* Calculator::getV02Means(int i, int pidIndex, double min, double max)
{
  std::vector<std::string> species = {"Ch", "Pi", "Ka", "Pr"};
  if (pidIndex > 3) {
    std::cout << "[Calculator::getV02]: Invalid pid index - must be 0 (ch), 1 (pi), 2 (ka), or 3 (pr)" << std::endl;
    return nullptr;
  }
  std::string sp = species[pidIndex];

  TH1D* mean_fraction_pt = getProjectedFractionalSpectra(sp, pidIndex, min, max, FractionScheme::V02);
  TH1D* v22 = getCNN(2, 2, "ChGap", i);
  TH1D* reth = (TH1D*)mean_fraction_pt->Clone(Form("v22_%snpt_%.1f_%.1f", sp.c_str(), min, max));
  reth->Scale(getIntegratedValue(v22, min, max));
  reth->SetName(Form("V02Means%s_%i_%.1f_%.1f", sp.c_str(), i, min, max));
  return reth;
}
TH1D* Calculator::getV02MeansK0(int i, double min, double max)
{
  TH1D* mean_fraction_pt_tot = getProjectedFractionalK0Spectra("K0_sig", 1, min, max, FractionScheme::V02);
  TH1D* v22 = getCNN(2, 2, "ChGap", i);
  TH1D* v22_npt_tot = (TH1D*)mean_fraction_pt_tot->Clone(Form("v22_npt_K0_tot_%.1f_%.1f", min, max));
  v22_npt_tot->Scale(getIntegratedValue(v22, min, max));
  if (useTotal) {
    v22_npt_tot->SetName(Form("V02MeansK0_%i", i));
    return v22_npt_tot;
  }
  TH1D* mean_fraction_pt_sb1 = getProjectedFractionalK0Spectra("K0_sb1", 0, min, max, FractionScheme::V02);
  TH1D* mean_fraction_pt_sb2 = getProjectedFractionalK0Spectra("K0_sb2", 2, min, max, FractionScheme::V02);
  TH1D* v22_npt_sb1 = (TH1D*)mean_fraction_pt_sb1->Clone(Form("v22_npt_K0_sb1_%.1f_%.1f", min, max));
  v22_npt_sb1->Scale(getIntegratedValue(v22, min, max));
  TH1D* v22_npt_sb2 = (TH1D*)mean_fraction_pt_sb2->Clone(Form("v22_npt_K0_sb2_%.1f_%.1f", min, max));
  v22_npt_sb2->Scale(getIntegratedValue(v22, min, max));

  TH1D* v22_npt_bkg = (TH1D*)v22_npt_sb1->Clone(Form("v22_npt_K0_bkg_%.1f_%.1f", min, max));
  v22_npt_bkg->Add(v22_npt_sb2, 1);
  v22_npt_bkg->Scale(0.5);

  TH1D* v22_npt_sig = (TH1D*)v22_npt_tot->Clone(Form("v22_npt_K0_sig_%.1f_%.1f", min, max));
  if (!useTotal)
    applySignalFractionCorrection(*this, v22_npt_sig, v22_npt_bkg, "K0", min, max, "Calculator::getV02MeansK0");
  v22_npt_sig->SetName(Form("V02MeansK0_%i", i));
  return v22_npt_sig;
}
TH1D* Calculator::getV02MeansLambda(int i, double min, double max)
{
  TH1D* mean_fraction_pt_tot = getProjectedFractionalLambdaSpectra("Lambda_sig", 1, min, max, FractionScheme::V02);
  TH1D* v22 = getCNN(2, 2, "ChGap", i);
  TH1D* v22_npt_tot = (TH1D*)mean_fraction_pt_tot->Clone(Form("v22_npt_Lambda_tot_%.1f_%.1f", min, max));
  v22_npt_tot->Scale(getIntegratedValue(v22, min, max));
  if (useTotal) {
    v22_npt_tot->SetName(Form("V02MeansLambda_%i", i));
    return v22_npt_tot;
  }
  TH1D* mean_fraction_pt_sb1 = getProjectedFractionalLambdaSpectra("Lambda_sb1", 0, min, max, FractionScheme::V02);
  TH1D* mean_fraction_pt_sb2 = getProjectedFractionalLambdaSpectra("Lambda_sb2", 2, min, max, FractionScheme::V02);
  TH1D* v22_npt_sb1 = (TH1D*)mean_fraction_pt_sb1->Clone(Form("v22_npt_Lambda_sb1_%.1f_%.1f", min, max));
  v22_npt_sb1->Scale(getIntegratedValue(v22, min, max));
  TH1D* v22_npt_sb2 = (TH1D*)mean_fraction_pt_sb2->Clone(Form("v22_npt_Lambda_sb2_%.1f_%.1f", min, max));
  v22_npt_sb2->Scale(getIntegratedValue(v22, min, max));

  TH1D* v22_npt_bkg = (TH1D*)v22_npt_sb1->Clone(Form("v22_npt_Lambda_bkg_%.1f_%.1f", min, max));
  v22_npt_bkg->Add(v22_npt_sb2, 1);
  v22_npt_bkg->Scale(0.5);

  TH1D* v22_npt_sig = (TH1D*)v22_npt_tot->Clone(Form("v22_npt_Lambda_sig_%.1f_%.1f", min, max));
  if (!useTotal)
    applySignalFractionCorrection(*this, v22_npt_sig, v22_npt_bkg, "Lambda", min, max, "Calculator::getV02MeansLambda");
  v22_npt_sig->SetName(Form("V02MeansLambda_%i", i));
  return v22_npt_sig;
}
TH1D* Calculator::getV02(int i, int pidIndex, double min, double max)
{
  std::vector<std::string> species = {"Ch", "Pi", "Ka", "Pr"};
  if (pidIndex > 3) {
    std::cout << "[Calculator::getV02]: Invalid pid index - must be 0 (ch), 1 (pi), 2 (ka), or 3 (pr)" << std::endl;
    return nullptr;
  }
  std::string sp = species[pidIndex];

  TH1D* mean_fraction_pt = getnpt(i, pidIndex, min, max, FractionScheme::V02);
  TH1D* reth = getCorrVsPt(2, 2, Form("ChGap%s", sp.c_str()), i, min, max);
  TH1D* v22_3pc = hasChGap3pcW22 ? getCNN(2, 2, "ChGap3pcW", i) : getCNN(2, 2, "ChGap", i);
  TH1D* sub = (TH1D*)mean_fraction_pt->Clone(Form("v22_%snpt_%.1f_%.1f", sp.c_str(), min, max));
  sub->Scale(getIntegratedValue(v22_3pc, min, max));
  TH1D* den = getV02Means(i, pidIndex, min, max);
  reth->Add(sub, -1);
  reth->Divide(den);
  reth->SetName(Form("V02%s_%i_%.1f_%.1f", sp.c_str(), i, min, max));
  return reth;
}
// Alternatives using the final V02
TH1D* Calculator::getV02K0Alt(int i, double min, double max)
{
  TH1D* mean_fraction_pt_sb1 = getProjectedFractionalK0Spectra("K0_sb1", 0, min, max, FractionScheme::V02);
  TH1D* mean_fraction_pt_tot = getProjectedFractionalK0Spectra("K0_sig", 1, min, max, FractionScheme::V02);
  TH1D* mean_fraction_pt_sb2 = getProjectedFractionalK0Spectra("K0_sb2", 2, min, max, FractionScheme::V02);
  TH1D* v02_tot = getCorrVsPt(2, 2, Form("ChGapK0s"), i, min, max);
  TH1D* v02_sb1 = getCorrVsPt(2, 2, Form("ChGapK0sSBL"), i, min, max);
  TH1D* v02_sb2 = getCorrVsPt(2, 2, Form("ChGapK0sSBR"), i, min, max);

  TH1D* v22 = getCNN(2, 2, "ChGap", i);
  TH1D* v22_npt_tot = (TH1D*)mean_fraction_pt_tot->Clone(Form("v02_K0_tot_%.1f_%.1f", min, max));
  v22_npt_tot->Scale(getIntegratedValue(v22, min, max));
  TH1D* v22_npt_sb1 = (TH1D*)mean_fraction_pt_sb1->Clone(Form("v02_K0_sb1_%.1f_%.1f", min, max));
  v22_npt_sb1->Scale(getIntegratedValue(v22, min, max));
  TH1D* v22_npt_sb2 = (TH1D*)mean_fraction_pt_sb2->Clone(Form("v02_K0_sb2_%.1f_%.1f", min, max));
  v22_npt_sb2->Scale(getIntegratedValue(v22, min, max));

  v02_tot->Add(v22_npt_tot, -1);
  v02_tot->Divide(v22_npt_tot);

  v02_sb1->Add(v22_npt_sb1, -1);
  v02_sb1->Divide(v22_npt_sb1);
  v02_sb2->Add(v22_npt_sb2, -1);
  v02_sb2->Divide(v22_npt_sb2);

  TH1D* v02_bkg = (TH1D*)v02_sb1->Clone(Form("v02_K0_bkg_%.1f_%.1f", min, max));
  v02_bkg->Add(v02_sb2, 1);
  v02_bkg->Scale(0.5);

  TH1D* v02_sig = (TH1D*)v02_tot->Clone(Form("v02_K0_sig_%.1f_%.1f", min, max));
  if (!useTotal)
    applySignalFractionCorrection(*this, v02_sig, v02_bkg, "K0", min, max, "Calculator::getV02K0");
  v02_sig->SetName(Form("V02K0_%i", i));
  return v02_sig;
}
TH1D* Calculator::getV02LambdaAlt(int i, double min, double max)
{
  TH1D* mean_fraction_pt_sb1 = getProjectedFractionalLambdaSpectra("Lambda_sb1", 0, min, max, FractionScheme::V02);
  TH1D* mean_fraction_pt_tot = getProjectedFractionalLambdaSpectra("Lambda_sig", 1, min, max, FractionScheme::V02);
  TH1D* mean_fraction_pt_sb2 = getProjectedFractionalLambdaSpectra("Lambda_sb2", 2, min, max, FractionScheme::V02);
  TH1D* v02_tot = getCorrVsPt(2, 2, Form("ChGapLambda"), i, min, max);
  TH1D* v02_sb1 = getCorrVsPt(2, 2, Form("ChGapLambdaSBL"), i, min, max);
  TH1D* v02_sb2 = getCorrVsPt(2, 2, Form("ChGapLambdaSBR"), i, min, max);

  TH1D* v22 = getCNN(2, 2, "ChGap", i);
  TH1D* v22_npt_tot = (TH1D*)mean_fraction_pt_tot->Clone(Form("v02_Lambda_tot_%.1f_%.1f", min, max));
  v22_npt_tot->Scale(getIntegratedValue(v22, min, max));
  TH1D* v22_npt_sb1 = (TH1D*)mean_fraction_pt_sb1->Clone(Form("v02_Lambda_sb1_%.1f_%.1f", min, max));
  v22_npt_sb1->Scale(getIntegratedValue(v22, min, max));
  TH1D* v22_npt_sb2 = (TH1D*)mean_fraction_pt_sb2->Clone(Form("v02_Lambda_sb2_%.1f_%.1f", min, max));
  v22_npt_sb2->Scale(getIntegratedValue(v22, min, max));

  v02_tot->Add(v22_npt_tot, -1);
  v02_tot->Divide(v22_npt_tot);

  v02_sb1->Add(v22_npt_sb1, -1);
  v02_sb1->Divide(v22_npt_sb1);
  v02_sb2->Add(v22_npt_sb2, -1);
  v02_sb2->Divide(v22_npt_sb2);

  TH1D* v02_bkg = (TH1D*)v02_sb1->Clone(Form("v02_Lambda_bkg_%.1f_%.1f", min, max));
  v02_bkg->Add(v02_sb2, 1);
  v02_bkg->Scale(0.5);

  TH1D* v02_sig = (TH1D*)v02_tot->Clone(Form("v02_Lambda_sig_%.1f_%.1f", min, max));
  if (!useTotal)
    applySignalFractionCorrection(*this, v02_sig, v02_bkg, "Lambda", min, max, "Calculator::getV02Lambda");
  v02_sig->SetName(Form("V02Lambda_%i", i));
  return v02_sig;
}

// Using signal fraction and fraction correlations
TH1D* Calculator::getV02K0(int i, double min, double max)
{
  TH1D* reth = getV02CorrK0(i, min, max);
  TH1D* v22 = getCNN(2, 2, "ChGap", i);
  TH1D* v22_3pc = hasChGap3pcW22 ? getCNN(2, 2, "ChGap3pcW", i) : getCNN(2, 2, "ChGap", i);
  TH1D* npt = getnptK0(i, min, max, FractionScheme::V02);
  TH1D* v22_npt = (TH1D*)npt->Clone(Form("v22_K0npt_%.1f_%.1f", min, max));
  v22_npt->Scale(getIntegratedValue(v22_3pc, min, max));
  TH1D* den = (TH1D*)npt->Clone(Form("den_K0_%.1f_%.1f", min, max));
  den->Scale(getIntegratedValue(v22, min, max));
  reth->Add(v22_npt, -1);
  reth->Divide(den);
  reth->SetName(Form("V02K0_%i", i));
  return reth;
}
TH1D* Calculator::getV02Lambda(int i, double min, double max)
{
  TH1D* reth = getV02CorrLambda(i, min, max);
  TH1D* v22 = getCNN(2, 2, "ChGap", i);
  TH1D* v22_3pc = hasChGap3pcW22 ? getCNN(2, 2, "ChGap3pcW", i) : getCNN(2, 2, "ChGap", i);
  TH1D* npt = getnptLambda(i, min, max, FractionScheme::V02);
  TH1D* v22_npt = (TH1D*)npt->Clone(Form("v22_Lambdanpt_%.1f_%.1f", min, max));
  v22_npt->Scale(getIntegratedValue(v22_3pc, min, max));
  TH1D* den = (TH1D*)npt->Clone(Form("den_Lambda_%.1f_%.1f", min, max));
  den->Scale(getIntegratedValue(v22, min, max));
  reth->Add(v22_npt, -1);
  reth->Divide(den);
  reth->SetName(Form("V02Lambda_%i", i));
  return reth;
}
TH1D* Calculator::getV0K0(int i, double min, double max)
{
  TH1D* reth = getV0CorrK0(i, min, max);
  TH1D* npt = getnptK0(i, min, max, FractionScheme::V0);
  TH1D* mpt = getMpt(i);
  /*     TList* l = (TList*)flowptcont->fSubList;
      TH1* pt1_sub1 = ((BootstrapProfile*)l->FindObject("ptptsub_1sub1_0sub2"))->getHist(i);
      TH1* pt1_sub2 = ((BootstrapProfile*)l->FindObject("ptptsub_0sub1_1sub2"))->getHist(i);
      TH1D* mpt = (TH1D*)pt1_sub1->Clone(Form("mpt_avg_%i_%.1f_%.1f",i));
      mpt->Add(pt1_sub2);
      mpt->Scale(0.5); */

  TH1D* mpt_npt = (TH1D*)npt->Clone(Form("mpt_K0npt_%.1f_%.1f", min, max));
  mpt_npt->Scale(getIntegratedValue(mpt, min, max));

  TH1D* mu2 = getMu2sub(i);
  mu2 = powerHistogram(mu2, 0.5);
  TH1D* mu_npt = (TH1D*)npt->Clone(Form("mu_K0npt_%.1f_%.1f", min, max));
  mu_npt->Scale(getIntegratedValue(mu2, min, max));

  reth->Add(mpt_npt, -1);
  reth->Divide(mu_npt);
  reth->SetName(Form("V0K0_%i", i));
  return reth;
}
TH1D* Calculator::getV0Lambda(int i, double min, double max)
{
  TH1D* reth = getV0CorrLambda(i, min, max);
  TH1D* npt = getnptLambda(i, min, max, FractionScheme::V0);
  TH1D* mpt = getMpt(i);
  /*     TList* l = (TList*)flowptcont->fSubList;
      TH1* pt1_sub1 = ((BootstrapProfile*)l->FindObject("ptptsub_1sub1_0sub2"))->getHist(i);
      TH1* pt1_sub2 = ((BootstrapProfile*)l->FindObject("ptptsub_0sub1_1sub2"))->getHist(i);
      TH1D* mpt = (TH1D*)pt1_sub1->Clone(Form("mpt_avg_%i_%.1f_%.1f",i));
      mpt->Add(pt1_sub2);
      mpt->Scale(0.5); */

  TH1D* mpt_npt = (TH1D*)npt->Clone(Form("mpt_Lambdanpt_%.1f_%.1f", min, max));
  mpt_npt->Scale(getIntegratedValue(mpt, min, max));

  TH1D* mu2 = getMu2sub(i);
  mu2 = powerHistogram(mu2, 0.5);
  TH1D* mu_npt = (TH1D*)npt->Clone(Form("mu_Lambdanpt_%.1f_%.1f", min, max));
  mu_npt->Scale(getIntegratedValue(mu2, min, max));

  reth->Add(mpt_npt, -1);
  reth->Divide(mu_npt);
  reth->SetName(Form("V0Lambda_%i", i));
  return reth;
}

// Get integrated and normalized V0(2)
TH1D* Calculator::getIntegratedV0(int i)
{
  TH1D* mu2 = getMu2sub(i);
  TH1D* mpt = getMpt(i);
  TH1D* reth = powerHistogram(mu2, 0.5);
  reth->Divide(mpt);
  reth->SetTitle("; Centrality (%); v_{0}");
  return reth;
}
TH1D* Calculator::getIntegratedV02(int i)
{
  TH1D* reth = getCVN2pt(2, i);
  TH1D* v22 = getV22(i);
  TH1D* mpt = getMpt(i);
  v22->Multiply(mpt);
  reth->Divide(v22);
  reth->SetTitle("; Centrality (%); v_{02}");
  return reth;
}

TH1D* Calculator::getAlpha(int i, double min, double max)
{
  // centrality-dependent observables
  TH1D* v0 = getIntegratedV0(i);
  TH1D* v02 = getIntegratedV02(i);
  TH1D* cv22pt = getCV22pt(i);
  TH1D* varV22 = getVarV22(i);
  TH1D* c22 = getC22(i);
  TH1D* getVarPt = getMu2sub(i);
  TH1D* sigmav22 = powerHistogram(varV22, 0.5);
  TH1D* sigmapt = powerHistogram(getVarPt, 0.5);

  // Get constant value in given centrality
  double v0_const = getIntegratedValue(v0, min, max);
  double v02_const = getIntegratedValue(v02, min, max);
  double cv22pt_const = getIntegratedValue(cv22pt, min, max);
  double varv22_const = getIntegratedValue(varV22, min, max);
  double c22_const = getIntegratedValue(c22, min, max);
  double sigmapt_const = getIntegratedValue(sigmapt, min, max);

  // pt-dependent quantities (Charged-particle only)
  TH1D* v0pt = getV0Weighted(i, 0, min, max);
  TH1D* v02pt = getV02Weighted(i, 0, min, max);

  double denominator = sigmapt_const * varv22_const * v0_const - cv22pt_const * v02_const * c22_const;
  v0pt->Scale(varv22_const * sigmapt_const / denominator);
  v02pt->Scale(cv22pt_const * c22_const / denominator);
  TH1D* alpha = (TH1D*)v0pt->Clone(Form("alphapt_%i_%g%g", i, min, max));
  alpha->Add(v02pt, -1);
  alpha->SetTitle("; #it{p}_{T} (GeV/#it{c}); #alpha");
  return alpha;
}

TH1D* Calculator::getBeta(int i, double min, double max)
{
  // centrality-dependent observables
  TH1D* v0 = getIntegratedV0(i);
  TH1D* v02 = getIntegratedV02(i);
  TH1D* cv22pt = getCV22pt(i);
  TH1D* varV22 = getVarV22(i);
  TH1D* c22 = getC22(i);
  TH1D* getVarPt = getMu2sub(i);
  TH1D* sigmav22 = powerHistogram(varV22, 0.5);
  TH1D* sigmapt = powerHistogram(getVarPt, 0.5);

  // Get constant value in given centrality
  double v0_const = getIntegratedValue(v0, min, max);
  double v02_const = getIntegratedValue(v02, min, max);
  double cv22pt_const = getIntegratedValue(cv22pt, min, max);
  double varv22_const = getIntegratedValue(varV22, min, max);
  double c22_const = getIntegratedValue(c22, min, max);
  double sigmapt_const = getIntegratedValue(sigmapt, min, max);

  // pt-dependent quantities (Charged-particle only)
  TH1D* v0pt = getV0Weighted(i, 0, min, max);
  TH1D* v02pt = getV02Weighted(i, 0, min, max);

  double denominator = sigmapt_const * varv22_const * v0_const - cv22pt_const * v02_const * c22_const;
  v0pt->Scale(sigmapt_const * v02_const * c22_const / denominator);
  v02pt->Scale(sigmapt_const * v0_const * c22_const / denominator);
  TH1D* beta = (TH1D*)v02pt->Clone(Form("betapt_%i_%g%g", i, min, max));
  beta->Add(v0pt, -1);
  beta->SetTitle("; #it{p}_{T} (GeV/#it{c}); #beta");
  return beta;
}

TH1D* Calculator::getVNN(int nhar, int npar, std::string id, int i)
{
  TH1D* reth = nullptr;
  auto clone = (i < 0) ? flowcont : dynamic_cast<FlowContainer*>(flowcont->Clone("TempClone"));
  // Should be set via input config (json?)
  // if(fPropagateFCerrors) clone->SetPropagateErrors(kTRUE);
  if (i >= 0) {
    clone->OverrideMainWithSub(i, kFALSE);
    Int_t nxb;
    Double_t* xbs = flowcont->GetMultiRebin(nxb);
    clone->SetMultiRebin(nxb, xbs);
    delete[] xbs;
  };
  clone->SetIDName(id.c_str());
  switch (npar) {
    case 2:
      reth = dynamic_cast<TH1D*>(clone->GetVN2VsX(nhar, kFALSE, 0));
      break;
    case 4:
      reth = dynamic_cast<TH1D*>(clone->GetVN4VsX(nhar, kFALSE, 0));
      break;
    case 6:
      reth = dynamic_cast<TH1D*>(clone->GetVN6VsX(nhar, kFALSE, 0));
      break;
    case 8:
      reth = dynamic_cast<TH1D*>(clone->GetVN8VsX(nhar, kFALSE, 0));
      break;
    default:
      printf("Number of particles in correlation should be multiple of 2 and only up to 8 (for now)\n");
      reth = nullptr;
      break;
  }
  return reth;
}
TH1D* Calculator::getVNN(int nhar, int npar, std::string id, int i, double min, double max)
{
  printf("here\n");
  TH1D* reth = nullptr;
  TH1D* rethP = nullptr;
  auto clone = (i < 0) ? flowcont : dynamic_cast<FlowContainer*>(flowcont->Clone("TempClone"));
  // Should be set via input config (json?)
  // if(fPropagateFCerrors) clone->SetPropagateErrors(kTRUE);
  if (i >= 0) {
    clone->OverrideMainWithSub(i, kFALSE);
    Int_t nxb;
    Double_t* xbs = flowcont->GetMultiRebin(nxb);
    clone->SetMultiRebin(nxb, xbs);
    delete[] xbs;
  };
  clone->SetIDName(id.c_str());
  switch (npar) {
    case 2:
      // clone->SetIDName("ChNGap");
      // reth = dynamic_cast<TH1D*>(clone->GetVN2VsX(nhar, kTRUE, min, max));
      // clone->SetIDName("ChPGap");
      reth = dynamic_cast<TH1D*>(clone->GetVN2VsX(nhar, kTRUE, min, max));
      // reth->Sumw2();
      // rethP->Sumw2();
      // reth->Add(rethP);
      // reth->Scale(0.5);
      break;
    case 4:
      reth = dynamic_cast<TH1D*>(clone->GetVN4VsX(nhar, kTRUE, min, max));
      break;
    case 6:
      reth = dynamic_cast<TH1D*>(clone->GetVN6VsX(nhar, kTRUE, min, max));
      break;
    case 8:
      reth = dynamic_cast<TH1D*>(clone->GetVN8VsX(nhar, kTRUE, min, max));
      break;
    default:
      printf("Number of particles in correlation should be multiple of 2 and only up to 8 (for now)\n");
      reth = nullptr;
      break;
  }
  return reth;
}
TH1D* Calculator::getCNN(int nhar, int npar, std::string id, int i)
{
  TH1D* reth = nullptr;
  std::unique_ptr<FlowContainer> ownedClone;
  FlowContainer* clone = getFlowContainerSample(i, ownedClone);
  if (!clone)
    return nullptr;
  // Should be set via input config (json?)
  // if(fPropagateFCerrors) clone->SetPropagateErrors(kTRUE);
  clone->SetIDName(id.c_str());
  switch (npar) {
    case 2:
      reth = dynamic_cast<TH1D*>(clone->GetCN2VsX(nhar, kFALSE, 0));
      break;
    case 4:
      reth = dynamic_cast<TH1D*>(clone->GetCN4VsX(nhar, kFALSE, 0));
      break;
    case 6:
      reth = dynamic_cast<TH1D*>(clone->GetCN6VsX(nhar, kFALSE, 0));
      break;
    case 8:
      reth = dynamic_cast<TH1D*>(clone->GetCN8VsX(nhar, kFALSE, 0));
      break;
    default:
      printf("Number of particles in correlation should be multiple of 2 and only up to 8 (for now)\n");
      reth = nullptr;
      break;
  }
  return reth;
}
TH1D* Calculator::getCNN(int nhar, int npar, std::string id, int i, double min, double max)
{
  TH1D* reth = nullptr;
  TH1D* rethP = nullptr;
  std::unique_ptr<FlowContainer> ownedClone;
  FlowContainer* clone = getFlowContainerSample(i, ownedClone);
  if (!clone)
    return nullptr;
  // Should be set via input config (json?)
  // if(fPropagateFCerrors) clone->SetPropagateErrors(kTRUE);
  clone->SetIDName(id.c_str());
  switch (npar) {
    case 2:
      clone->SetIDName("ChNGap");
      reth = dynamic_cast<TH1D*>(clone->GetCN2VsX(nhar, kTRUE, min, max));
      clone->SetIDName("ChPGap");
      rethP = dynamic_cast<TH1D*>(clone->GetCN2VsX(nhar, kTRUE, min, max));
      reth->Sumw2();
      rethP->Sumw2();
      reth->Add(rethP);
      reth->Scale(0.5);
      delete rethP;
      break;
    case 4:
      reth = dynamic_cast<TH1D*>(clone->GetCN4VsX(nhar, kTRUE, min, max));
      break;
    case 6:
      reth = dynamic_cast<TH1D*>(clone->GetCN6VsX(nhar, kTRUE, min, max));
      break;
    case 8:
      reth = dynamic_cast<TH1D*>(clone->GetCN8VsX(nhar, kTRUE, min, max));
      break;
    default:
      printf("Number of particles in correlation should be multiple of 2 and only up to 8 (for now)\n");
      reth = nullptr;
      break;
  }
  return reth;
}
TH1D* Calculator::getCorrVsPt(int nhar, int npar, std::string id, int i, double min, double max)
{
  TH1D* reth = nullptr;
  std::unique_ptr<FlowContainer> ownedClone;
  FlowContainer* clone = getFlowContainerSample(i, ownedClone);
  if (!clone)
    return nullptr;
  // Should be set via input config (json?)
  // if(fPropagateFCerrors) clone->SetPropagateErrors(kTRUE);
  clone->SetIDName(id.c_str());
  std::string order = std::to_string(nhar) + std::to_string(npar);
  reth = dynamic_cast<TH1D*>(clone->GetHistCorrXXVsPt(order.c_str(), min, max));
  return reth;
}
TH1D* Calculator::getNN(int nhar, int npar, std::string id, int i)
{
  auto clone = (i < 0) ? flowcont : dynamic_cast<FlowContainer*>(flowcont->Clone("TempClone"));
  // Should be set via input config (json?)
  // if(fPropagateFCerrors) clone->SetPropagateErrors(kTRUE);
  if (i >= 0) {
    clone->OverrideMainWithSub(i, kFALSE);
    Int_t nxb;
    Double_t* xbs = flowcont->GetMultiRebin(nxb);
    clone->SetMultiRebin(nxb, xbs);
    delete[] xbs;
  };
  clone->SetIDName(id.c_str());
  return dynamic_cast<TH1D*>(clone->GetHistCorrXXVsMulti(Form("%i%u", nhar, npar), 0));
}
TH1D* Calculator::powerHistogram(TH1D* inhist, double power)
{
  TH1D* outhist = dynamic_cast<TH1D*>(inhist->Clone(Form("power%.1f%s", power, inhist->GetName())));
  for (int bin = 1; bin <= inhist->GetNbinsX(); ++bin) {
    if (inhist->GetBinContent(bin) >= 0) {
      outhist->SetBinContent(bin, std::pow(inhist->GetBinContent(bin), power));
      outhist->SetBinError(bin, power * std::pow(inhist->GetBinContent(bin), power - 1) * inhist->GetBinError(bin));
    } else {
      outhist->SetBinContent(bin, -2);
      outhist->SetBinError(bin, 0);
    }
  }
  return outhist;
}
double Calculator::getIntegratedValue(TH1D* inhist, double min, double max)
{
  double sum = 0.;
  double sumw = 0.;
  int minbin = inhist->GetXaxis()->FindBin(min + 0.00001);
  int maxbin = inhist->GetXaxis()->FindBin(max - 0.00001);
  for (int bin = minbin; bin <= maxbin; ++bin) {
    double bc = inhist->GetBinContent(bin);
    double be = inhist->GetBinError(bin);
    if (be == 0.0)
      continue;
    sum += bc / (be * be);
    sumw += 1. / (be * be);
  }
  return sum / sumw;
}
double Calculator::readInvMassFit(std::string species, double massMin, double massMax, double pTMin, double pTMax, int centMin, int centMax)
{
  if (massMin < 0. || massMax < 0.) {
    if (species == "K0") {
      massMin = 0.488;
      massMax = 0.507;
    } else if (species == "Lambda" || species == "AntiLambda") {
      massMin = 1.112;
      massMax = 1.120;
    } else {
      std::cout << "[Calculator::readInvMassFit] Unknown species " << species << " and no mass window provided" << std::endl;
      return -1.;
    }
  }

  // Late-centrality calculations are performed in the native (usually 1%)
  // FlowContainer bins. Reuse the invariant-mass fit from the configured bin
  // that fully contains the native interval, rather than requiring unstable
  // invariant-mass fits for every native bin.
  bool foundConfiguredFitBin = false;
  for (std::size_t i = 0; i + 1 < invMassFitCentralityBins.size(); ++i) {
    const int fitCentMin = fitCentEdge(invMassFitCentralityBins[i]);
    const int fitCentMax = fitCentEdge(invMassFitCentralityBins[i + 1]);
    if (centMin >= fitCentMin && centMax <= fitCentMax) {
      centMin = fitCentMin;
      centMax = fitCentMax;
      foundConfiguredFitBin = true;
      break;
    }
  }
  if (!foundConfiguredFitBin && invMassFitLookupCentMin >= 0 && invMassFitLookupCentMax > invMassFitLookupCentMin &&
      centMin >= invMassFitLookupCentMin && centMax <= invMassFitLookupCentMax) {
    centMin = invMassFitLookupCentMin;
    centMax = invMassFitLookupCentMax;
  }

  const std::string fitFunctionKey = (species == "AntiLambda" || species == "LambdaAntiLambda") ? "Lambda" : species;
  const std::pair<std::string, std::string>* defaultFitFunctions = nullptr;
  const auto defaultFunctionsIt = invMassFitFunctions.find(fitFunctionKey);
  if (defaultFunctionsIt != invMassFitFunctions.end()) {
    defaultFitFunctions = &defaultFunctionsIt->second;
  }
  const std::pair<std::string, std::string>* fitFunctions = defaultFitFunctions;
  const auto systematicFunctionsIt = invMassFitFunctionsBySystematic.find(systematicIndex);
  if (systematicFunctionsIt != invMassFitFunctionsBySystematic.end()) {
    const auto functionsIt = systematicFunctionsIt->second.find(fitFunctionKey);
    if (functionsIt != systematicFunctionsIt->second.end()) {
      fitFunctions = &functionsIt->second;
    }
  }
  const bool hasFitFunctions = fitFunctions != nullptr;
  const std::string signalFunction = hasFitFunctions ? fitFunctions->first : "";
  const std::string backgroundFunction = hasFitFunctions ? fitFunctions->second : "";
  std::string modelTag = hasFitFunctions ? signalFunction + "_" + backgroundFunction : "";
  if (!hasFitFunctions) {
    std::cout << "[Calculator::readInvMassFit] No fit-function configuration for "
              << fitFunctionKey << "; using legacy untagged lookup" << std::endl;
  }

  TString cacheKey = Form("%s|%.6f|%.6f|%g|%g|%i|%i|%i|%s|%s",
                          species.c_str(), massMin, massMax, pTMin, pTMax, centMin, centMax,
                          systematicIndex, invMassFitPath.c_str(), modelTag.c_str());
  auto cacheIt = invMassFitCache.find(cacheKey.Data());
  if (cacheIt != invMassFitCache.end()) {
    return cacheIt->second;
  }

  // Get normalised system name
  std::string system_str = collision_system;
  for (char& c : system_str) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  }
  if (system_str == "oo" || system_str == "o-o")
    system_str = "OO";
  else if (system_str == "ne-ne" || system_str == "nene")
    system_str = "Ne-Ne";
  else if (system_str == "mc" || system_str == "mcreco" || system_str == "reco")
    system_str = "MC";
  else if (system_str == "mcgen" || system_str == "gen")
    system_str = "MCGen";
  else if (system_str == "LHC25h3c" || system_str == "lhc25h3c")
    system_str = "LHC25h3c";
  else if (system_str == "LHC25h4c" || system_str == "lhc25h4c")
    system_str = "LHC25h4c";
  else if (system_str == "LHC26c4" || system_str == "lhc26c4")
    system_str = "LHC26c4";
  else
    system_str = "Pb-Pb";

  TString fileName;
  const TString systematicDirectory = systematicIndex > 0 ? Form("sys_%03d/", systematicIndex) : "";
  if (hasFitFunctions) {
    fileName = Form("../pidQA/invmass_outputs/%s/%sinvmass_%s_%s_%i_%i_Pt_%g_%g.root",
                    system_str.c_str(), systematicDirectory.Data(), species.c_str(), modelTag.c_str(),
                    centMin, centMax, pTMin, pTMax);
  } else {
    fileName = Form("../pidQA/invmass_outputs/%s/%sinvmass_%s_%i_%i_Pt_%g_%g.root",
                    system_str.c_str(), systematicDirectory.Data(), species.c_str(), centMin, centMax, pTMin, pTMax);
  }
  if (!invMassFitPath.empty()) {
    TString basePath = invMassFitPath.c_str();
    if (!basePath.EndsWith("/")) {
      basePath += "/";
    }
    fileName = basePath + fileName;
  }
  TFile* infile = TFile::Open(fileName, "READ");
  if (!infile || infile->IsZombie()) {
    delete infile;
    infile = nullptr;

    TString sourceLog = Form("../pidQA/invmass_outputs/%s/sys_%03d/mass_fit_source.txt",
                             system_str.c_str(), systematicIndex);
    if (!invMassFitPath.empty()) {
      TString basePath = invMassFitPath.c_str();
      if (!basePath.EndsWith("/")) {
        basePath += "/";
      }
      sourceLog = basePath + sourceLog;
    }
    std::ifstream sourceLogFile(sourceLog.Data());
    std::string sourceDescription((std::istreambuf_iterator<char>(sourceLogFile)),
                                  std::istreambuf_iterator<char>());
    const std::string fallbackSpecies = species == "LambdaAntiLambda" ? "Lambda" : species;
    const bool useNominalFallback = systematicIndex > 0 && defaultFitFunctions &&
                                    sourceDescription.find("fallback") != std::string::npos &&
                                    sourceDescription.find(" " + fallbackSpecies) != std::string::npos;
    if (useNominalFallback) {
      modelTag = defaultFitFunctions->first + "_" + defaultFitFunctions->second;
      fileName = Form("../pidQA/invmass_outputs/%s/invmass_%s_%s_%i_%i_Pt_%g_%g.root",
                      system_str.c_str(), species.c_str(), modelTag.c_str(),
                      centMin, centMax, pTMin, pTMax);
      if (!invMassFitPath.empty()) {
        TString basePath = invMassFitPath.c_str();
        if (!basePath.EndsWith("/")) {
          basePath += "/";
        }
        fileName = basePath + fileName;
      }
      infile = TFile::Open(fileName, "READ");
      const std::string warningKey = std::to_string(systematicIndex) + "|" + fitFunctionKey;
      if (infile && !infile->IsZombie() && invMassFitFallbackWarnings.insert(warningKey).second) {
        std::cout << "[Calculator::readInvMassFit] Systematic index " << systematicIndex
                  << " uses nominal invariant-mass fits for " << fitFunctionKey << std::endl;
      }
    }
    if (!infile || infile->IsZombie()) {
      std::cout << "[Calculator::readInvMassFit] Could not open " << fileName
                << " in invmass_outputs/" << system_str << std::endl;
      delete infile;
      invMassFitCache[cacheKey.Data()] = -1.;
      return -1.;
    }
  }

  TString totalName;
  TString bkgName;
  if (hasFitFunctions) {
    totalName = Form("InvMassRecoFits/f_%s_%s_total_pt_%g_%g",
                     species.c_str(), modelTag.c_str(), pTMin, pTMax);
    bkgName = Form("InvMassRecoFits/f_%s_%s_background_pt_%g_%g",
                   species.c_str(), modelTag.c_str(), pTMin, pTMax);
  } else {
    totalName = Form("InvMassRecoFits/f_%s_sig_pt_%g_%g", species.c_str(), pTMin, pTMax);
    bkgName = Form("InvMassRecoFits/f_%s_bkg_final_pt_%g_%g", species.c_str(), pTMin, pTMax);
  }
  TF1* totalFit = dynamic_cast<TF1*>(infile->Get(totalName));
  TF1* bkgFit = dynamic_cast<TF1*>(infile->Get(bkgName));
  if (!totalFit || !bkgFit) {
    std::cout << "[Calculator::readInvMassFit] Missing fit functions in " << fileName << std::endl;
    infile->Close();
    delete infile;
    invMassFitCache[cacheKey.Data()] = -1.;
    return -1.;
  }

  const double totalIntegral = totalFit->Integral(massMin, massMax);
  const double bkgIntegral = bkgFit->Integral(massMin, massMax);
  infile->Close();
  delete infile;

  if (totalIntegral <= 0.) {
    std::cout << "[Calculator::readInvMassFit] Non-positive total integral for " << species << std::endl;
    invMassFitCache[cacheKey.Data()] = -1.;
    return -1.;
  }

  const double signalFraction = (totalIntegral - bkgIntegral) / totalIntegral;
  invMassFitCache[cacheKey.Data()] = signalFraction;
  return signalFraction;
}
bool Calculator::is_digits(std::string& str)
{
  for (char ch : str) {
    if (!isdigit(ch)) {
      return false;
    }
  }
  return true;
}
std::string Calculator::padWithZeros(int number, int totalDigits)
{
  std::ostringstream oss;
  oss << std::setw(totalDigits) << std::setfill('0') << number;
  return oss.str();
}
std::pair<int, int> Calculator::checkSubeventInput(FlowPtContainer* cont)
{
  TList* sublist = cont->fSubList;
  std::vector<std::string> names;
  TIter next(sublist);
  TObject* obj = 0;
  while ((obj = next()))
    names.push_back(obj->GetName());
  std::regex re(R"((\d+)sub(\d+))");
  int mpar = 0;
  int nsub = 0;
  for (const auto& s : names) {
    auto begin = std::sregex_iterator(s.begin(), s.end(), re);
    auto end = std::sregex_iterator();

    for (auto it = begin; it != end; ++it) {
      int m = std::stoi((*it)[1].str());
      int sub = std::stoi((*it)[2].str());
      mpar = std::max(mpar, m);
      nsub = std::max(nsub, sub);
    }
  }
  std::cout << "Maximum number of correlated particles: " << mpar << std::endl;
  std::cout << "Maximum number of subevents: " << nsub << std::endl;
  return std::make_pair(mpar, nsub);
}
// NOLINTEND(readability/casting, readability/braces)
