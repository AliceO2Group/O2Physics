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
#include "SysHandler.h"

#include "TF1.h"
#include "TList.h"
#include "TMath.h"

#include <cstdio>
#include <string>
#include <utility>
#include <vector>

namespace
{
bool hasValidSystematicRange(const SysHandler::SystSettings& settings)
{
  return settings.sys_avg_range.second > settings.sys_avg_range.first;
}

bool isInSystematicRange(const TH1* hist, int bin, const SysHandler::SystSettings& settings)
{
  if (!hasValidSystematicRange(settings)) {
    return true;
  }
  const double center = hist->GetXaxis()->GetBinCenter(bin);
  return center >= settings.sys_avg_range.first && center <= settings.sys_avg_range.second;
}

double getReportedGroupUncertainty(const TH1* hist, const SysHandler::SystSettings& settings, bool relative)
{
  if (!hist) {
    return 0.;
  }
  if (settings.use_fit_for_syst) {
    for (int bin = 1; bin <= hist->GetNbinsX(); ++bin) {
      if (isInSystematicRange(hist, bin, settings)) {
        return hist->GetBinContent(bin);
      }
    }
    return 0.;
  }
  if (!hasValidSystematicRange(settings)) {
    return hist->GetMaximum();
  }

  double sum = 0.;
  int acceptedBins = 0;
  for (int bin = 1; bin <= hist->GetNbinsX(); ++bin) {
    if (!isInSystematicRange(hist, bin, settings)) {
      continue;
    }
    const double value = hist->GetBinContent(bin);
    if (value == 0. || (relative && value > settings.max_allowed_syst)) {
      continue;
    }
    sum += value;
    ++acceptedBins;
  }
  return acceptedBins > 0 ? sum / acceptedBins : 0.;
}
} // namespace

SysHandler::SysHandler() : current_settings{{0., 90.}, true, 1, 1, true}
{
  setSystematicVariations({
    {1, 1, -1, "VtxZ < 9cm", "VtxZ"},
    {2, 1, -1, "VtxZ < 8cm", "VtxZ"},
    {3, 1, -1, "VtxZ < 7cm", "VtxZ"},
    {4, 2, -1, "DCAxy < 5#sigma", "DCAxy"},
    {5, 2, -1, "DCAxy < 10#sigma", "DCAxy"},
    {6, 3, -1, "DCAz < 0.1cm", "DCAz"},
    {7, 3, -1, "DCAz < 2cm", "DCAz"},
    {8, 4, -1, "TPC X-rows > 60", "TPCXrows"},
    {9, 4, -1, "TPC X-rows > 90", "TPCXrows"},
    {10, 5, -1, "TPCncls > 40", "TPC cluster"},
    {11, 5, -1, "TPCncls > 70", "TPC cluster"},
    {12, 6, -1, "#chi^2 pr. TPCcls < 4", "TPCChi2"},
    {13, 10, -1, "c#tau", "Lifetime"},
    {14, 0, -1, "FT0M", "centEst"},
    {15, 11, -1, "Competing mass cut", "Competing mass cut"},
    {16, 12, -1, "#cos(#theta_{p})", "Pointing Angle"},
    {17, 0, -1, "Daughter PID rejection", "Daughter PID rejection"},
    {18, 8, -1, "FitCentrality", "fitCentrality"},
    {19, 9, -1, "MassFitFunction", "fitFunction"},
  });
}
SysHandler::~SysHandler()
{
  clearSysMap();
};

void SysHandler::clearSysMap()
{
  for (auto& [index, histogram] : fSystMap) {
    delete histogram;
  }
  fSystMap.clear();
}

void SysHandler::setSystematicVariations(const std::vector<SystematicVariation>& variations)
{
  fSystematicVariations.clear();
  fIndMap.clear();
  for (const auto& variation : variations) {
    fSystematicVariations[variation.index] = variation;
  }
}
TH1D* SysHandler::getNominal(bool relative, bool apply_barlow)
{
  if (fSystMap.empty()) {
    printf("Warning: SystMap not set!");
    return nullptr;
  }
  TH1* systs = getSystematicsSummed(relative, apply_barlow);
  systs->SetName("Systematics_Summed");
  TH1* nominal = dynamic_cast<TH1*>(fSystMap[0]->Clone(Form("Nominal_fullSyst_%s", fSystMap[0]->GetName())));
  // FSRebin(&cvals);
  applyErrors(nominal, systs, relative);
  delete systs;
  return dynamic_cast<TH1D*>(nominal);
}
TH1* SysHandler::getSystematics(int sys_index, bool rel, bool apply_barlow, bool abs)
{
  TH1* var = dynamic_cast<TH1*>(fSystMap[sys_index]->Clone(Form("RelSyst_Ind%i_%s", sys_index, fSystMap[0]->GetName())));
  TH1* nom = dynamic_cast<TH1*>(fSystMap[0]->Clone(Form("Nominal_Ind%i_%s", sys_index, fSystMap[0]->GetName())));
  var->Add(nom, -1);
  if (rel)
    var->Divide(nom);
  histOper(var, 0, 1, abs, kTRUE);
  delete nom;
  if (!apply_barlow)
    return var;
  TH1* barlow = getBarlowTest(sys_index, getCorrelationError(sys_index));
  if (!barlow)
    return var;
  if (!passesBarlowSelection(barlow)) {
    for (int i = 1; i <= var->GetNbinsX(); ++i) {
      var->SetBinContent(i, 0.);
    }
  }
  delete barlow;
  return var;
}
void SysHandler::applyBarlow(TH1* target, TH1* barlow)
{
  if (!target || !barlow) {
    printf("ApplyBarlow: missing a histogram...\n");
    return;
  };
  if (target->GetNbinsX() != barlow->GetNbinsX()) {
    printf("applyBarlow: different number of bins\n");
    return;
  };
  for (Int_t i = 1; i <= target->GetNbinsX(); i++) {
    if (barlow->GetBinContent(i) > current_settings.fNbarlow)
      continue;
    target->SetBinContent(i, 0);
  };
}
bool SysHandler::passesBarlowSelection(const TH1* barlow) const
{
  if (!barlow) {
    return true;
  }
  double sum = 0.;
  int binCount = 0;
  int barlowCount = 0;
  for (int i = 1; i <= barlow->GetNbinsX(); ++i) {
    if (!isInSystematicRange(barlow, i, current_settings)) {
      continue;
    }
    sum += barlow->GetBinContent(i);
    ++binCount;
    if (barlow->GetBinContent(i) > current_settings.fNbarlow) {
      ++barlowCount;
    }
  }
  if (binCount == 0) {
    return true;
  }
  if (current_settings.fUseBarlowCount) {
    const int barlowCutoff = ceil(0.5 * binCount);
    return barlowCount >= barlowCutoff;
  }
  return sum / binCount >= current_settings.fNbarlow;
}
void SysHandler::histOper(TH1* reth, Double_t cvl, Double_t lpower, Bool_t Abs, Bool_t RemoveErrors)
{
  for (Int_t i = 1; i <= reth->GetNbinsX(); i++) {
    Double_t bc = reth->GetBinContent(i);
    if (!bc)
      continue;
    bc += cvl;
    Double_t bcn = TMath::Power(bc, lpower);
    if (Abs)
      bcn = TMath::Abs(bcn);
    reth->SetBinContent(i, bcn);
    if (RemoveErrors) {
      reth->SetBinError(i, 0);
    } else {
      Double_t err = reth->GetBinError(i);
      err = err * lpower * TMath::Power(bc, lpower - 1);
      reth->SetBinError(i, err);
    }
  }
}
TH1* SysHandler::getSystematicsFromFit(Int_t sys_index, Bool_t apply_barlow, bool ratio, bool relative)
{
  if (fSystMap.empty())
    printf("Sys map not set\n");
  const auto nominalIt = fSystMap.find(0);
  const auto variationIt = fSystMap.find(sys_index);
  if (nominalIt == fSystMap.end() || !nominalIt->second || variationIt == fSystMap.end() || !variationIt->second) {
    printf("Could not fetch nominal or systematic index %i\n", sys_index);
    return nullptr;
  }
  TH1* var = dynamic_cast<TH1*>(variationIt->second->Clone(Form("RelSystFit_sys%i_%s", sys_index, nominalIt->second->GetName())));
  if (!var) {
    printf("Could not fetch systematic index %i\n", sys_index);
    return 0;
  };
  // If we have final state rebinning, then need bootstrapping on the nominal here. Otherwise, final state rebinning will be bogus
  TH1* nom = dynamic_cast<TH1*>(nominalIt->second->Clone(Form("NominalFit_sys%i_%s", sys_index, nominalIt->second->GetName())));
  // FSRebin(&nom);
  var->SetName(Form("RelSyst_%s_Ind%i", nom->GetName(), sys_index));
  if (relative || ratio) {
    var->Divide(nom);
  } else {
    var->Add(nom, -1.);
  }
  histOper(var, 0, 1, kTRUE, kFALSE);
  if (ratio) {
    delete nom;
    return var;
  }
  var->Fit("pol0", "Q0+R", "", current_settings.sys_avg_range.first, current_settings.sys_avg_range.second);
  TF1* func = (TF1*)var->GetListOfFunctions()->FindObject("pol0");
  const double fittedUncertainty = relative ? TMath::Abs(1. - func->GetParameter(0)) : TMath::Abs(func->GetParameter(0));
  for (int i = 1; i <= var->GetNbinsX(); ++i) {
    var->SetBinContent(i, isInSystematicRange(var, i, current_settings) ? fittedUncertainty : 0.);
  }
  histOper(var, 0, 1, kFALSE, kTRUE);
  delete nom;
  if (!apply_barlow)
    return var;
  TH1* barlow = getBarlowTest(sys_index, getCorrelationError(sys_index));
  if (!barlow)
    return var;
  if (!passesBarlowSelection(barlow)) {
    for (int i = 1; i <= var->GetNbinsX(); ++i) {
      var->SetBinContent(i, 0.);
    }
  }
  delete barlow;
  return var;
}
TH1* SysHandler::getBarlowTest(Int_t sys_index, Double_t corrFactor)
{
  TH1* var = dynamic_cast<TH1*>(fSystMap[sys_index]->Clone(Form("VarBarlow_Ind%i_%s", sys_index, fSystMap[0]->GetName())));
  if (!var)
    return 0;
  TH1* nom = dynamic_cast<TH1*>(fSystMap[0]->Clone(Form("NominalBarlow_Ind%i_%s", sys_index, fSystMap[0]->GetName())));
  TH1* absDif = getSystematics(sys_index, kFALSE, kFALSE, kTRUE);
  if (!absDif) {
    delete nom;
    delete var;
    return 0;
  }
  for (Int_t i = 1; i <= absDif->GetNbinsX(); i++) {
    Double_t bc = absDif->GetBinContent(i);
    if (!bc)
      continue;
    Double_t s1 = var->GetBinError(i);
    Double_t s2 = nom->GetBinError(i);
    Double_t denom = TMath::Sqrt(TMath::Abs((s1 * s1 + corrFactor * s2 * s2)));
    if (!denom) {
      absDif->SetBinError(i, 0);
      continue;
    };
    absDif->SetBinContent(i, bc / denom);
    absDif->SetBinError(i, 0);
  }
  delete nom;
  delete var;
  return absDif;
}
Bool_t SysHandler::BuildIndexMap(Bool_t force)
{
  // if(!fSysts) {printf("Systematics not read!\n"); return 0; };
  // if(!fSysts->GetEntries()) {printf("No systematics available!\n"); return 0; };
  if (force && !((int)fIndMap.size()))
    fIndMap.clear();
  for (auto& el : fSystMap) {
    Int_t tval = getTarget(el.first);
    if (!fIndMap.count(tval))
      fIndMap.insert(make_pair(tval, std::vector<int>{}));
    fIndMap.at(tval).push_back(el.first);
  };
  return kTRUE;
};

std::vector<TH1*> SysHandler::getSystSubset(Int_t KeyInd, Bool_t relative, Bool_t apply_barlow)
{
  if (!(int)fIndMap.size())
    if (!BuildIndexMap(kTRUE)) {
      printf("getSystSubset: could not build map!\n");
      return nullvec;
    };
  if (!fIndMap.count(KeyInd)) {
    printf("getSystSubset: could not find key %i in the map!\n", KeyInd);
    return nullvec;
  };
  std::vector<int> inds = fIndMap.at(KeyInd);
  std::vector<TH1*> retvec;
  for (Int_t j = 0; j < (int)inds.size(); j++)
    retvec.push_back(current_settings.use_fit_for_syst ? getSystematicsFromFit(inds.at(j), apply_barlow, false, relative) : getSystematics(inds.at(j), relative, apply_barlow));
  return retvec;
}
std::vector<TH1*> SysHandler::getMergedErrors(Bool_t relative, Bool_t apply_barlow)
{
  if (!(int)fIndMap.size())
    if (!BuildIndexMap(kTRUE)) {
      printf("getMergedErrors: could not build map!\n");
      return nullvec;
    };
  if (!fPrintMap.empty())
    fPrintMap.erase(fPrintMap.begin(), fPrintMap.end());
  fPrintBinMap.clear();
  const auto storePrintValues = [this, relative](int index, const TH1* hist) {
    fPrintMap[index] = getReportedGroupUncertainty(hist, current_settings, relative);
    std::vector<double>& binValues = fPrintBinMap[index];
    binValues.reserve(hist->GetNbinsX());
    for (int bin = 1; bin <= hist->GetNbinsX(); ++bin) {
      binValues.push_back(hist->GetBinContent(bin));
    }
  };
  std::vector<TH1*> reth;
  for (auto errcind = fIndMap.begin(); errcind != fIndMap.end(); errcind++) {
    Int_t targetInd = errcind->first;
    if (targetInd == 0)
      continue; // 0 contains all those that are not relevant
    std::vector<TH1*> lSubsets = getSystSubset(targetInd, relative, apply_barlow);
    if ((int)lSubsets.size() == 0)
      continue; // if no subsets, then continue
    TH1* rh = (TH1*)lSubsets.at(0)->Clone(Form("MergedErrors_%i", targetInd));
    if (!rh) {
      for (auto* subset : lSubsets)
        delete subset;
      continue;
    }
    RemoveErrors(rh); // To remove errors
    reth.push_back(rh);
    if (fSystOverride.count(targetInd)) {
      Double_t newerr = fSystOverride.at(targetInd);
      printf("Found a bins with key index %i, value %f\n", targetInd, newerr);
      for (Int_t i = 1; i <= rh->GetNbinsX(); i++)
        rh->SetBinContent(i, newerr);
      storePrintValues(targetInd, rh);
      for (auto* subset : lSubsets)
        delete subset;
      continue;
    }
    if ((int)lSubsets.size() == 1) { // if only one, then no need to continue, but need to add to print map
      storePrintValues(targetInd, rh);
      delete lSubsets.front();
      continue;
    }
    // loop over subsets and pick max
    for (auto sh = lSubsets.begin() + 1; sh != lSubsets.end(); sh++) {
      // Looping over bins
      for (Int_t i = 1; i <= rh->GetNbinsX(); i++) {
        Double_t vmax = TMath::Max(rh->GetBinContent(i), (*sh)->GetBinContent(i));
        rh->SetBinContent(i, vmax);
      }
    }
    // Ideally only use this for systematics from fit as it has a single value in all bins, otherwise we get the max value but might not be as informative
    storePrintValues(targetInd, rh);
    for (auto* subset : lSubsets)
      delete subset;
  }
  return reth;
}
TH1* SysHandler::getSystematicsSummed(Bool_t relative, Bool_t apply_barlow)
{
  std::vector<TH1*> mergedErrors = getMergedErrors(relative, apply_barlow);
  if ((int)mergedErrors.size() == 0) {
    printf("getSystematicsSummed: could not find any errors to sum\n");
    return 0;
  };
  TH1* reth = (TH1*)mergedErrors.at(0)->Clone("MergedErrors_Summed");
  reth->Reset();
  for (auto hist : mergedErrors) {
    addInQuad(reth, hist);
    hist->Delete();
  };
  if (relative && hasValidSystematicRange(current_settings) && !current_settings.use_fit_for_syst) {
    int acceptedBins = 0;
    double sum = 0.;
    for (int bin = 1; bin <= reth->GetNbinsX(); ++bin) {
      if (!isInSystematicRange(reth, bin, current_settings)) {
        continue;
      }
      const double value = reth->GetBinContent(bin);
      if (value == 0. || value > current_settings.max_allowed_syst) {
        continue;
      }
      sum += value;
      ++acceptedBins;
    }
    if (acceptedBins > 0) {
      const double average = sum / acceptedBins;
      for (int bin = 1; bin <= reth->GetNbinsX(); ++bin) {
        reth->SetBinContent(bin, isInSystematicRange(reth, bin, current_settings) ? average : 0.);
      }
    }
  }
  return reth;
};
void SysHandler::addInQuad(TH1* h1, TH1* h2)
{
  if (!h1 || !h2) {
    printf("addInQuad: One of the histograms is missing!\n");
    return;
  };
  if (h1->GetNbinsX() != h2->GetNbinsX()) {
    printf("addInQuad: number of bins in provided histograms are not the same!\n");
    return;
  };
  for (Int_t i = 1; i <= h1->GetNbinsX(); i++) {
    Double_t bc1 = h1->GetBinContent(i);
    Double_t bc2 = h2->GetBinContent(i);
    if (bc1 == 0 && bc2 == 0)
      continue;
    h1->SetBinContent(i, TMath::Sqrt(bc1 * bc1 + bc2 * bc2));
  }
}
void SysHandler::applyErrors(TH1* inh, TH1* errors, Bool_t relative)
{
  if (!inh || !errors) {
    printf("ApplyErrors: One of the histograms is missing!\n");
    return;
  };
  // if(inh->GetNbinsX()!=errors->GetNbinsX() && fNRebinForSyst<2) { printf("ApplyErrors: number of bins in provided histograms are not the same!\n"); return; };
  if (inh->GetNbinsX() != errors->GetNbinsX()) {
    printf("ApplyErrors: number of bins in provided histograms are not the same!\n");
    return;
  };
  for (Int_t i = 1; i <= inh->GetNbinsX(); i++) {
    Double_t bc1 = inh->GetBinContent(i);
    Double_t bce = inh->GetBinError(i);
    if (bc1 == 0 && bce == 0)
      continue;
    double bc2 = errors->GetBinContent(i);
    /*     Double_t bc2;
        if(fNRebinForSyst>1) {
          Double_t bcent = inh->GetBinCenter(i);
          Int_t tind = errors->FindBin(bcent);
          bc2 = errors->GetBinContent(tind);
        } else bc2 = errors->GetBinContent(i); */
    bce = relative ? (bc1 * bc2) : bc2;
    inh->SetBinError(i, bce);
  }
}
int SysHandler::getTarget(int fSystInd) const
{
  const auto variation = fSystematicVariations.find(fSystInd);
  return variation == fSystematicVariations.end() ? 0 : variation->second.group;
}
int SysHandler::getCorrelationError(int fSystInd) const
{
  const auto variation = fSystematicVariations.find(fSystInd);
  return variation == fSystematicVariations.end() ? 0 : variation->second.correlation;
}
std::string SysHandler::getSystDescription(int fSystInd) const
{
  const auto variation = fSystematicVariations.find(fSystInd);
  return variation == fSystematicVariations.end() ? "Not a valid sys ID" : variation->second.description;
}
// NOLINTEND(readability/casting, readability/braces)
std::string SysHandler::getSystGroupDescription(int targetInd) const
{
  if (targetInd == 0) {
    return "Not a valid sys ID";
  }
  for (const auto& [index, variation] : fSystematicVariations) {
    if (variation.group == targetInd) {
      return variation.groupDescription;
    }
  }
  return "Not a valid sys ID";
}
