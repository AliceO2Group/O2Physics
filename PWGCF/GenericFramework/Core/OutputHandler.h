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

#ifndef PWGCF_GENERICFRAMEWORK_CORE_OUTPUTHANDLER_H_
#define PWGCF_GENERICFRAMEWORK_CORE_OUTPUTHANDLER_H_

#include "Calculator.h"
#include "FlowContainer.h"
#include "FlowPtContainer.h"
#include "StatHandler.h"
#include "SysHandler.h"

#include "TH1.h"

#include <cstdio>
#include <functional>
#include <map>
#include <string>
#include <utility>
#include <vector>

class OutputHandler : public TNamed
{
 public:
  explicit OutputHandler(std::string system = "");
  OutputHandler(std::string filename, std::string system);
  ~OutputHandler();

  void initialise(std::string filename, std::string dirname = "", bool hasSystematics = false, SysHandler::SystSettings settings = {{0., 90.}, true, 1, 1, true});
  void initialiseGen(std::string filename, std::string dirname = "")
  {
    isGen = true;
    initialise(filename, dirname);
  }
  void addFileForRatio(std::string name)
  {
    doSystemRatio = true;
    fileName2 = name;
  };

  TH1D* getHist(Calculator::ObsFunc func, StatHandler::StatMethod method = StatHandler::kVarianceOfSamples);
  TH1D* getHist(Calculator::PtDiffFunc func, double min, double max, StatHandler::StatMethod method = StatHandler::kVarianceOfSamples);
  TH1D* getSysHist(Calculator::ObsFunc func, bool use_bootstrap_for_systs, StatHandler::StatMethod method = StatHandler::kVarianceOfSamples, bool relative = true, const TH1* nominal = nullptr);
  TH1D* getSysHist(Calculator::PtDiffFunc func, double min, double max, bool use_bootstrap_for_systs, StatHandler::StatMethod method = StatHandler::kVarianceOfSamples, bool relative = true, const TH1* nominal = nullptr);
  void RebinMulti(int nrb, double* chbins);
  template <typename T>
  void FSRebin(T** inh);
  TH1* getRelSysHist(int sys_index, bool sys_from_fit, bool relative = true);
  TH1* getBarlowHist(int sys_index);
  TH1* getVarNomRatioHist(int sys_index);
  TH1* getVarHist(int sys_index);
  std::vector<std::string> split(const std::string& s, char delimiter);
  std::string getLatexName(Calculator::ObsFunc func);
  std::string getLatexName(Calculator::PtDiffFunc func);
  std::string getName(Calculator::ObsFunc func);
  std::string getName(Calculator::PtDiffFunc func, double min, double max, bool toplevelname = 0);
  void setFlowID(std::string id)
  {
    if (isInitialised)
      calculator.setFlowID(id);
    else
      printf("Please initialise OutputHandler before setting Calculator data members\n");
  }
  void setInvMassFitFunctions(const std::string& species, const std::string& signalFunction, const std::string& backgroundFunction);
  void setInvMassFitFunctionsForSystematic(int systematicIndex,
                                           const std::string& species,
                                           const std::string& signalFunction,
                                           const std::string& backgroundFunction);
  void setInvMassFitCentralityBins(const std::vector<double>& centralityBins);
  void setAdditionalSystematicVariation(int index, TH1* variation);
  void setSystematicVariations(const std::vector<SysHandler::SystematicVariation>& variations);
  void setSystematicProgressCallback(std::function<void(std::size_t, std::size_t, int)> callback) { systematicProgressCallback = std::move(callback); }
  void setUseFSRebin(bool newval) { useFSRebin = newval; }
  void setUseLateCentralityRebin(bool newval) { useLateCentralityRebin = newval; }
  void setUseInverseVarianceWeighting(bool newval) { calculator.setUseInverseVarianceWeight(newval); }
  void setRadialFlowAnalysis(bool newval, bool useTotal = false)
  {
    radialFlowAnalysis = newval;
    calculator.setUseTotal(useTotal);
  }
  std::vector<int> getSysIndices() { return sys_indices; }
  int getSysColor(int index)
  {
    if (index < 0 || index >= static_cast<int>(syscolors.size())) {
      throw std::out_of_range("Index out of range");
    }
    return syscolors[index];
  }
  std::string getSysName(int index) { return sysHandler->getSystDescription(index); }
  void PrintSystematics(bool relative = true);

 private:
  FlowContainer* flowcont;     //!
  FlowPtContainer* flowptcont; //!
  StatHandler* stathandler;    //!
  SysHandler* sysHandler;      //!
  std::string fileName;
  std::string dirName;
  Calculator calculator; //!
  bool isInitialised;
  bool isGen;
  bool useFSRebin;
  bool useLateCentralityRebin;
  bool doSystemRatio;
  bool radialFlowAnalysis;
  bool isVN(Calculator::ObsFunc func);
  bool isVN(Calculator::PtDiffFunc func);
  void calculateVN(Calculator::ObsFunc func, TH1D*& inh);
  Calculator::PtDiffFunc getCorrespondingCumulant(Calculator::PtDiffFunc func);
  Calculator::ObsFunc getCorrespondingCumulant(Calculator::ObsFunc func);
  void calculateVN(Calculator::PtDiffFunc func, TH1D*& inh, double min, double max);
  std::string current_obs;
  std::vector<int> sys_indices;
  int fNrb;
  double* fBinsrb;
  std::string xAxisTitle;
  std::string fileName2;
  std::map<std::string, TFile*> inputFiles;                                      //!
  std::map<int, TH1*> additionalSystematicVariations;                            //!
  std::vector<SysHandler::SystematicVariation> systematicVariations;             //!
  bool hasCustomSystematicVariations = false;                                    //!
  std::function<void(std::size_t, std::size_t, int)> systematicProgressCallback; //!
  // predefined systematic colors
  std::vector<int> syscolors = {kBlack, kRed + 1, kBlue + 1, kGreen + 3, kMagenta + 1, kOrange - 1, kCyan + 2, kYellow + 2, kGray + 1, kRed - 10, kBlue - 9, kGreen - 8, kMagenta - 9, kOrange - 9, kCyan - 8, kYellow - 7, kAzure + 1, kPink, kP6Blue, kP6Yellow, kP6Red, kP6Grape, kP6Violet};
  // std::vector<int> syscolors = {126,125,124,123,122,121,120,119,118,117};

  TFile* readFile(std::string name);
  std::vector<int> getKeyIndices(std::string name);

  ClassDef(OutputHandler, 3);
};

#endif // PWGCF_GENERICFRAMEWORK_CORE_OUTPUTHANDLER_H_
