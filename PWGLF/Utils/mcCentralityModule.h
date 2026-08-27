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

///
/// \file   mcCentralityModule.cxx
/// \author Romain Schotter romain.schotter@cern.ch
/// \brief  Module to produce MC centrality table based on ANY table containing MC multiplicity
///

#ifndef PWGLF_UTILS_MCCENTRALITYMODULE_H_
#define PWGLF_UTILS_MCCENTRALITYMODULE_H_

#include "PWGLF/DataModel/mcCentrality.h"

#include <CCDB/BasicCCDBManager.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/Array2D.h>
#include <Framework/Configurable.h>
#include <Framework/DataSpecUtils.h>
#include <Framework/DeviceSpec.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/RunningWorkflowInfo.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TRandom3.h>
#include <TString.h>

#include <RtypesCore.h>

#include <concepts>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <memory>
#include <string>
#include <type_traits>
#include <vector>

//__________________________________________
// strangeness builder module

namespace o2
{
namespace pwglf
{
namespace mccentrality // avoid polluting other namespaces
{

// statics necessary for the configurables in this namespace
static constexpr int nParameters = 1;
static const std::vector<std::string> tableNames{
  "McCentFV0As",
  "McCentFT0Ms",
  "McCentFT0As",
  "McCentFT0Cs",
  "McCentFT0CVariant1s",
  "McCentFT0CVariant2s",
  "McCentFDDMs",
  "McCentNTPVs",
  "McCentNGlobals",
  "McCentMFTs"};

static constexpr int nTablesConst = 10;
static const std::vector<std::string> parameterNames{"enable"};
static const int defaultParameters[nTablesConst][nParameters]{
  {-1},
  {-1},
  {-1},
  {-1},
  {-1},
  {-1},
  {-1},
  {-1},
  {-1},
  {-1}};

// table index : match order above
enum tableIndex { kFV0A = 0,
                  kFT0M,
                  kFT0A,
                  kFT0C,
                  kFT0CVariant1,
                  kFT0CVariant2,
                  kFDDM,
                  kNTPV,
                  kNGlobal,
                  kMFT,
                  kNestimators };

static constexpr const char* DirList[] = {
  "FV0A",
  "FT0M",
  "FT0A",
  "FT0C",
  "FT0CVariant1",
  "FT0CVariant2",
  "FDDM",
  "NTPV",
  "NGlobal",
  "MFT"};

// mcCentralityModule: 1st-order configurables
struct coreConfigurables : o2::framework::ConfigurableGroup {
  o2::framework::Configurable<o2::framework::LabeledArray<int>> enabledTables{"enabledTables",
                                                                              {defaultParameters[0], nTablesConst, nParameters, tableNames, parameterNames},
                                                                              "Produce this table: -1 for autodetect; otherwise, 0/1 is false/true"};
  std::vector<int> mEnabledTables; // Vector of enabled tables

  o2::framework::Configurable<bool> recalibrateCentrality{"recalibrateCentrality", false, "If true, re-calibrate the MC centrality for the binning in binPercentile."};
  o2::framework::Configurable<int> recalibrateMode{"recalibrateMode", 1, "Strategy to calibrate MC centrality? 0: from low to high mult.; 1: from high to low mult."};
  o2::framework::Configurable<int> minEntries{"minEntries", 100, "Minimum number of entries for estimating the mean value for the recalibration"};
  o2::framework::Configurable<bool> doNotCrashOnNull{"doNotCrashOnNull", false, "If ccdb object does not exist, fill with dummy values"};

  o2::framework::Configurable<bool> assignCentralityPerCandidate{"assignCentralityPerCandidate", false, "If true, assign centrality by sampling P(centrality|generated multiplicity) per MC collision instead of the class-averaged (statistical) value"};
  o2::framework::Configurable<ULong64_t> centralitySamplingSeed{"centralitySamplingSeed", 137, "Base seed of the RNG used to sample centrality in candidate-by-candidate mode (assignCentralityPerCandidate)"};
  o2::framework::ConfigurableAxis binsPercentile{"binsPercentile", {o2::framework::VARIABLE_WIDTH, 0, 1.0, 5.0, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0}, "Binning of the percentile axis"};
  o2::framework::ConfigurableAxis binsPercentileFine{"binsPercentileFine", {o2::framework::VARIABLE_WIDTH, 0, 0.001, 0.01, 1.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0}, "Binning of the percentile axis"};
  o2::framework::ConfigurableAxis binsMultiplicity{"binsMultiplicity", {1000, 0, 5000}, "Binning of the multiplicity axis"};

  // ccdb configurables
  o2::framework::Configurable<std::string> path{"path", "/tmp/InputCalibMC.root", "path to calib file or ccdb path if begins with ccdb://"};

  // debug option
  o2::framework::Configurable<bool> verbose{"verbose", false, "If true, display more messages"};
};

struct products : o2::framework::ProducesGroup {
  // Tables to produce
  o2::framework::Produces<aod::McCentFV0As> centFV0A;
  o2::framework::Produces<aod::McCentFT0Ms> centFT0M;
  o2::framework::Produces<aod::McCentFT0As> centFT0A;
  o2::framework::Produces<aod::McCentFT0Cs> centFT0C;
  o2::framework::Produces<aod::McCentFT0CVariant1s> centFT0CVariant1;
  o2::framework::Produces<aod::McCentFT0CVariant2s> centFT0CVariant2;
  o2::framework::Produces<aod::McCentFDDMs> centFDDM;
  o2::framework::Produces<aod::McCentNTPVs> centNTPV;
  o2::framework::Produces<aod::McCentNGlobals> centNGlobal;
  o2::framework::Produces<aod::McCentMFTs> centMFT;
};

template <typename T>
concept HasMcMults = requires(typename T::iterator a) {
  { a.multMCFT0A() } -> std::convertible_to<int>;
  { a.multMCFT0C() } -> std::convertible_to<int>;
  { a.multMCFV0A() } -> std::convertible_to<int>;
  { a.multMCFDDA() } -> std::convertible_to<int>;
  { a.multMCFDDC() } -> std::convertible_to<int>;
  { a.multMCNParticlesEta08() } -> std::convertible_to<int>;
  { a.multMCNParticlesEta05() } -> std::convertible_to<int>;
};

/// Task to produce the response table
struct BuilderModule {
  // Input parameters
  o2::framework::Service<o2::ccdb::BasicCCDBManager> ccdb;

  // declaration of structs here
  // (N.B.: will be invisible to the outside, create your own copies)
  coreConfigurables baseOpts;

  TList* MCCentralityCalibObjects = nullptr;
  TH1D* h1dFV0A = nullptr;
  TH1D* h1dFT0M = nullptr;
  TH1D* h1dFT0A = nullptr;
  TH1D* h1dFT0C = nullptr;
  TH1D* h1dFT0CVariant1 = nullptr;
  TH1D* h1dFT0CVariant2 = nullptr;
  TH1D* h1dFDDM = nullptr;
  TH1D* h1dNTPV = nullptr;
  TH1D* h1dNGlobal = nullptr;
  TH1D* h1dMFT = nullptr;

  // Calibration objects for candidate-by-candidate sampling (assignCentralityPerCandidate):
  // x = centrality (%) of the matched reconstructed MC collision, y = generated multiplicity (|eta| < 0.5)
  std::vector<TH2D*> h2dCentVsGenMult;
  TRandom3 fRandomSampler;

  // QA of the on-the-fly calibration (recalibrateCentrality): mean <N_PV>_{|eta|<0.5} in data and MC
  // (reco), one entry per centrality class -- filled in extractCentralityCalibration().
  // Must be TProfile, not TH1D+SetBinContent/SetBinError: every grid job redundantly derives the exact
  // same calibration from the same CCDB input, and grid-merged output histograms are combined via a
  // plain bin-content sum (TH1::Merge), which would turn a mean into N_jobs*mean. TProfile stores the
  // underlying sum(w)/sum(wy)/sum(wy^2) per bin and merges those correctly, so the mean/error recomputed
  // from the merged output stay correct however many jobs contributed. For the same reason there is no
  // ratio histogram here: a ratio of two already-merged means can't be reconstructed after the fact, so
  // compute PVMC/PVData downstream from the merged output if you need it, not inside the task.
  std::vector<std::shared_ptr<TProfile>> hCalibPVData;
  std::vector<std::shared_ptr<TProfile>> hCalibPVMC;

  int nEnabledTables = 0;
  int mRunNumber;

  // TAxis
  int nCentBins;
  std::vector<double> centralityBins;

  // Registers the on-the-fly-calibration QA histograms for one estimator (no-op unless recalibrateCentrality
  // is on): mean <N_PV>_{|eta|<0.5} in data and in MC (reco), one entry per centrality class, as TProfile
  // (see the comment on hCalibPVData/hCalibPVMC for why TH1D+SetBinContent doesn't survive grid merging).
  // Uses the runtime add() overload (returning the type-erased HistPtr) since the estimator name/directory
  // is only known at runtime here, unlike fillHistograms() which relies on a compile-time HIST() lookup.
  template <typename THistoRegistry>
  void registerCalibQAHistos(THistoRegistry& histos, int idx, const char* name)
  {
    if (!baseOpts.recalibrateCentrality) {
      return;
    }
    o2::framework::AxisSpec axisCentClass{baseOpts.binsPercentile, Form("%s percentile (%%)", name)};
    hCalibPVData[idx] = std::get<std::shared_ptr<TProfile>>(histos.add(Form("%s/calibPVData", name), Form("#LT #it{N}_{PV}^{Data} #GT vs %s percentile;%s percentile (%%);#LT #it{N}_{PV}^{Data} #GT_{|#it{#eta}|<0.5}", name, name), o2::framework::HistType::kTProfile, {axisCentClass}));
    hCalibPVMC[idx] = std::get<std::shared_ptr<TProfile>>(histos.add(Form("%s/calibPVMC", name), Form("#LT #it{N}_{PV}^{MC} #GT (reco) vs %s percentile;%s percentile (%%);#LT #it{N}_{PV}^{MC} #GT_{|#it{#eta}|<0.5} (reco)", name, name), o2::framework::HistType::kTProfile, {axisCentClass}));
  }

  // Registers the candidate-by-candidate sampling QA (no-op unless assignCentralityPerCandidate is on):
  // sampled percentile vs the generated multiplicity (|eta| < 0.5) actually used to key the sampling, so
  // it can be visually compared against the input calibration (hGenMultEta05VsCentrality<estimator>) to
  // check that the sampled distribution reproduces the expected P(centrality|genMult). This is a plain
  // TH2D (not a TProfile): unlike the calibration QA above, it stores a distribution to look at, not a
  // mean to merge, so ordinary bin-content-sum merging across grid jobs is exactly what's wanted here.
  template <typename THistoRegistry>
  void registerCandidateQAHistos(THistoRegistry& histos, const char* name)
  {
    if (!baseOpts.assignCentralityPerCandidate) {
      return;
    }
    histos.add(Form("%s/percentileVsGenMultEta05", name), Form("Sampled %s percentile vs generated mult (|#it{#eta}|<0.5);%s percentile (%%);generated mult (|#it{#eta}|<0.5)", name, name), o2::framework::HistType::kTH2D, {{baseOpts.binsPercentileFine, Form("%s percentile", name)}, {baseOpts.binsMultiplicity, "generated mult (|#eta|<0.5)"}});
  }

  template <typename TBaseConfigurables, typename THistoRegistry, typename TInitContext>
  void init(TBaseConfigurables const& inputBaseOpts, THistoRegistry& histos, TInitContext& context)
  {
    // read in configurations from the task where it's used
    // could be grouped even further, but should work
    baseOpts = inputBaseOpts;

    baseOpts.mEnabledTables.resize(nTablesConst, 0);
    h2dCentVsGenMult.resize(nTablesConst, nullptr);
    hCalibPVData.resize(nTablesConst);
    hCalibPVMC.resize(nTablesConst);
    fRandomSampler.SetSeed(baseOpts.centralitySamplingSeed);

    LOGF(info, "Checking if MC centrality is required");
    auto& workflows = context.services().template get<o2::framework::RunningWorkflowInfo const>();

    nEnabledTables = 0;

    TString listOfRequestors[nTablesConst];
    for (int i = 0; i < nTablesConst; i++) {
      int f = baseOpts.enabledTables->get(tableNames[i].c_str(), "enable");
      if (f == 1) {
        baseOpts.mEnabledTables[i] = 1;
        listOfRequestors[i] = "manual enabling";
        nEnabledTables++;
      }
      if (f == -1) {
        // autodetect this table in other devices
        for (o2::framework::DeviceSpec const& device : workflows.devices) {
          // Step 1: check if this device subscribed to the V0data table
          for (auto const& input : device.inputs) {
            if (o2::framework::DataSpecUtils::partialMatch(input.matcher, o2::header::DataOrigin("AOD"))) {
              auto&& [origin, description, version] = o2::framework::DataSpecUtils::asConcreteDataMatcher(input.matcher);
              std::string tableNameWithVersion = tableNames[i];
              if (version > 0) {
                tableNameWithVersion += Form("_%03d", version);
              }
              if (input.matcher.binding == tableNameWithVersion) {
                LOGF(info, "Device %s has subscribed to %s (version %i)", device.name, tableNames[i], version);
                listOfRequestors[i].Append(Form("%s ", device.name.c_str()));
                baseOpts.mEnabledTables[i] = 1;
                nEnabledTables++;
              }
            }
          }
        }
      }
    }

    if (nEnabledTables == 0) {
      LOGF(info, "MC centrality not required. Will suppress all functionality, including logs, from this point forward.");
      return;
    }
    mRunNumber = 0;

    // TAxis
    o2::framework::AxisSpec axisBinsPercentile{baseOpts.binsPercentile, "axisBinsPercentile"};
    nCentBins = axisBinsPercentile.binEdges.size() - 1;
    for (std::size_t iCent = 0; iCent < axisBinsPercentile.binEdges.size(); iCent++) {
      centralityBins.push_back(axisBinsPercentile.binEdges[iCent]);
    }

    if (baseOpts.mEnabledTables[kFV0A]) {
      histos.add("FV0A/percentile", "FV0A percentile.", o2::framework::HistType::kTH1D, {{baseOpts.binsPercentileFine, "FV0A percentile"}});
      histos.add("FV0A/percentilevsMult", "FV0A percentile.", o2::framework::HistType::kTH2D, {{baseOpts.binsPercentileFine, "FV0A percentile"}, {baseOpts.binsMultiplicity, "FV0A mult."}});
      registerCalibQAHistos(histos, kFV0A, "FV0A");
      registerCandidateQAHistos(histos, "FV0A");
    }
    if (baseOpts.mEnabledTables[kFT0M]) {
      histos.add("FT0M/percentile", "FT0M percentile.", o2::framework::HistType::kTH1D, {{baseOpts.binsPercentileFine, "FT0M percentile"}});
      histos.add("FT0M/percentilevsMult", "FT0M percentile.", o2::framework::HistType::kTH2D, {{baseOpts.binsPercentileFine, "FT0M percentile"}, {baseOpts.binsMultiplicity, "FT0M mult."}});
      registerCalibQAHistos(histos, kFT0M, "FT0M");
      registerCandidateQAHistos(histos, "FT0M");
    }
    if (baseOpts.mEnabledTables[kFT0A]) {
      histos.add("FT0A/percentile", "FT0A percentile.", o2::framework::HistType::kTH1D, {{baseOpts.binsPercentileFine, "FT0A percentile"}});
      histos.add("FT0A/percentilevsMult", "FT0A percentile.", o2::framework::HistType::kTH2D, {{baseOpts.binsPercentileFine, "FT0A percentile"}, {baseOpts.binsMultiplicity, "FT0A mult."}});
      registerCalibQAHistos(histos, kFT0A, "FT0A");
      registerCandidateQAHistos(histos, "FT0A");
    }
    if (baseOpts.mEnabledTables[kFT0C]) {
      histos.add("FT0C/percentile", "FT0C percentile.", o2::framework::HistType::kTH1D, {{baseOpts.binsPercentileFine, "FT0C percentile"}});
      histos.add("FT0C/percentilevsMult", "FT0C percentile.", o2::framework::HistType::kTH2D, {{baseOpts.binsPercentileFine, "FT0C percentile"}, {baseOpts.binsMultiplicity, "FT0C mult."}});
      registerCalibQAHistos(histos, kFT0C, "FT0C");
      registerCandidateQAHistos(histos, "FT0C");
    }
    if (baseOpts.mEnabledTables[kFT0CVariant1]) {
      histos.add("FT0CVariant1/percentile", "FT0CVariant1 percentile.", o2::framework::HistType::kTH1D, {{baseOpts.binsPercentileFine, "FT0CVariant1 percentile"}});
      histos.add("FT0CVariant1/percentilevsMult", "FT0CVariant1 percentile.", o2::framework::HistType::kTH2D, {{baseOpts.binsPercentileFine, "FT0CVariant1 percentile"}, {baseOpts.binsMultiplicity, "FT0C mult."}});
      registerCalibQAHistos(histos, kFT0CVariant1, "FT0CVariant1");
      registerCandidateQAHistos(histos, "FT0CVariant1");
    }
    if (baseOpts.mEnabledTables[kFT0CVariant2]) {
      histos.add("FT0CVariant2/percentile", "FT0CVariant2 percentile.", o2::framework::HistType::kTH1D, {{baseOpts.binsPercentileFine, "FT0CVariant2 percentile"}});
      histos.add("FT0CVariant2/percentilevsMult", "FT0CVariant2 percentile.", o2::framework::HistType::kTH2D, {{baseOpts.binsPercentileFine, "FT0CVariant2 percentile"}, {baseOpts.binsMultiplicity, "FT0C mult."}});
      registerCalibQAHistos(histos, kFT0CVariant2, "FT0CVariant2");
      registerCandidateQAHistos(histos, "FT0CVariant2");
    }
    if (baseOpts.mEnabledTables[kFDDM]) {
      histos.add("FDDM/percentile", "FDDM percentile.", o2::framework::HistType::kTH1D, {{baseOpts.binsPercentileFine, "FDDM percentile"}});
      histos.add("FDDM/percentilevsMult", "FDDM percentile.", o2::framework::HistType::kTH2D, {{baseOpts.binsPercentileFine, "FDDM percentile"}, {baseOpts.binsMultiplicity, "FDDM mult."}});
      registerCalibQAHistos(histos, kFDDM, "FDDM");
      registerCandidateQAHistos(histos, "FDDM");
    }
    if (baseOpts.mEnabledTables[kNTPV]) {
      histos.add("NTPV/percentile", "NTPV percentile.", o2::framework::HistType::kTH1D, {{baseOpts.binsPercentileFine, "NTPV percentile"}});
      histos.add("NTPV/percentilevsMult", "NTPV percentile.", o2::framework::HistType::kTH2D, {{baseOpts.binsPercentileFine, "NTPV percentile"}, {baseOpts.binsMultiplicity, "NTPV mult."}});
      registerCalibQAHistos(histos, kNTPV, "NTPV");
      registerCandidateQAHistos(histos, "NTPV");
    }
    if (baseOpts.mEnabledTables[kNGlobal]) {
      histos.add("NGlobal/percentile", "NGlobal percentile.", o2::framework::HistType::kTH1D, {{baseOpts.binsPercentileFine, "NGlobal percentile"}});
      histos.add("NGlobal/percentilevsMult", "NGlobal percentile.", o2::framework::HistType::kTH2D, {{baseOpts.binsPercentileFine, "NGlobal percentile"}, {baseOpts.binsMultiplicity, "NGlobal mult."}});
      registerCalibQAHistos(histos, kNGlobal, "NGlobal");
      registerCandidateQAHistos(histos, "NGlobal");
    }
    if (baseOpts.mEnabledTables[kMFT]) {
      histos.add("MFT/percentile", "MFT percentile.", o2::framework::HistType::kTH1D, {{baseOpts.binsPercentileFine, "MFT percentile"}});
      histos.add("MFT/percentilevsMult", "MFT percentile.", o2::framework::HistType::kTH2D, {{baseOpts.binsPercentileFine, "MFT percentile"}, {baseOpts.binsMultiplicity, "MFT mult."}});
    }

    LOGF(info, "*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*");
    LOGF(info, " MC centrality: basic configuration listing");
    LOGF(info, "*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*");
    // list enabled tables
    for (int i = 0; i < nTablesConst; i++) {
      // printout to be improved in the future
      if (baseOpts.mEnabledTables[i]) {
        LOGF(info, " -~> Table enabled: %s, requested by %s", tableNames[i], listOfRequestors[i].Data());
      }
    }
  }

  template <typename THist>
  THist* getHist(const char* name)
  {
    if (!this->MCCentralityCalibObjects) {
      return (THist*)0x0;
    }
    auto hist = reinterpret_cast<THist*>(this->MCCentralityCalibObjects->FindObject(name));
    if (!hist) {
      if (this->baseOpts.verbose) {
        this->MCCentralityCalibObjects->ls();
      }
      if (this->baseOpts.doNotCrashOnNull) {
        LOG(info) << "Could not open histogram " << name << " from TList, will fill tables with dummy values";
      } else {
        LOG(fatal) << "Could not open histogram " << name << " from TList";
      }
    }
    return hist;
  }

  TH1D* extractCentralityCalibration(int idx, const char* name, bool reverse = false)
  {

    auto CalibMC = [&](TString const& estimator) {
      LOGF(info, "Calibration for %s estimator -> Starting...\n", estimator.Data());

      std::vector<double> percentile_center(nCentBins);
      std::vector<double> epercentile_center(nCentBins);

      // Histograms
      TH2D* h2dMultVsCent_Data = getHist<TH2D>(Form("hMultEta05VsCent%s_Data", estimator.Data()));
      TH2D* h2dMultRecoVsMultGen_MC = getHist<TH2D>(Form("hMultEta05VsGenMult%s", estimator.Data()));
      if (!h2dMultVsCent_Data || !h2dMultRecoVsMultGen_MC) {
        return (TH1D*)0x0;
      }

      // QA of this on-the-fly calibration: mean <N_PV>_{|eta|<0.5} in data/MC per centrality class,
      // registered by registerCalibQAHistos() (only non-null when recalibrateCentrality is on). Filled
      // bin-by-bin from projData/projMC below (TProfile::Fill(x, y, weight)), not via SetBinContent, so
      // that the mean/error reconstructed from the merged grid output stay statistically correct.
      TProfile* hPVData = hCalibPVData[idx].get();
      TProfile* hPVMC = hCalibPVMC[idx].get();

      TH1D* h1dCalib = h2dMultRecoVsMultGen_MC->ProjectionX(Form("h1d%s", estimator.Data()), 1, h2dMultRecoVsMultGen_MC->GetNbinsX());
      h1dCalib->Reset();
      h1dCalib->SetTitle(Form("%s calibration object", estimator.Data()));
      h1dCalib->GetXaxis()->SetTitle(Form("#it{N}_{%s, gen.}", estimator.Data()));
      h1dCalib->GetYaxis()->SetTitle(Form("%s percentile (%%)", estimator.Data()));

      // NOTE: candidate-by-candidate assignment (assignCentralityPerCandidate) does not go through this
      // mean-matching path at all. It samples directly from the "hGenMultEta05VsCentrality<estimator>"
      // joint histogram (x = reco centrality of the matched MC collision, y = generated mult |eta|<0.5)
      // retrieved in initCCDB() below, see h2dCentVsGenMult and dataProcess().

      if (reverse) {
        for (int i = 0; i < nCentBins; i++) {
          int irev = i;
          percentile_center[i] = (centralityBins[irev] + centralityBins[irev + 1]) / 2;
          epercentile_center[i] = (centralityBins[irev] - centralityBins[irev + 1]) / 2;
        }

        int startBinMc = h1dCalib->GetNbinsX();
        for (int i = 0; i < nCentBins; i++) { // Loop over centrality bins
          // start from the end (from the high multiplicity collisions)
          int irev = i;
          TH1D* projData = h2dMultVsCent_Data->ProjectionY(Form("projData_%d", i), h2dMultVsCent_Data->GetXaxis()->FindBin(centralityBins[irev] + 1e-5),
                                                           h2dMultVsCent_Data->GetXaxis()->FindBin(centralityBins[irev + 1] - 1e-5));

          double meanMult_Data = projData->GetMean();

          double meanMult_MC = -1;
          double diffMultDataVsMC = 1e+09;
          int endBinMc = startBinMc;
          for (int j = startBinMc; j >= 1; j--) {
            // Loop over MC bins
            TH1D* projMC = h2dMultRecoVsMultGen_MC->ProjectionY("", j, startBinMc);
            int nEntries = projMC->Integral();
            double meanMC = projMC->GetMean();
            double ldiff = std::abs(meanMC - meanMult_Data);
            // std::cout << ldiff << std::endl;

            // if less than 100 entries for estimating the mean, do not consider it
            if (nEntries < baseOpts.minEntries)
              continue;

            // Find the minimum difference and corresponding MC bin
            if (ldiff < diffMultDataVsMC) {
              diffMultDataVsMC = ldiff;
              meanMult_MC = meanMC;
              endBinMc = j;
            }
          }
          if (i == nCentBins - 1) {
            endBinMc = 1;
          }

          TH1D* projMC = h2dMultRecoVsMultGen_MC->ProjectionY(Form("projMC_%d", i), endBinMc, startBinMc);

          // Refill the per-class mean from the underlying (data/MC) multiplicity distributions, bin by
          // bin, so the TProfile accumulates real sum(w)/sum(wy)/sum(wy^2) rather than a precomputed mean.
          if (hPVData) {
            for (int ybin = 1; ybin <= projData->GetNbinsX(); ybin++) {
              double w = projData->GetBinContent(ybin);
              if (w > 0) {
                hPVData->Fill(percentile_center[i], projData->GetBinCenter(ybin), w);
              }
            }
          }
          if (hPVMC) {
            for (int ybin = 1; ybin <= projMC->GetNbinsX(); ybin++) {
              double w = projMC->GetBinContent(ybin);
              if (w > 0) {
                hPVMC->Fill(percentile_center[i], projMC->GetBinCenter(ybin), w);
              }
            }
          }

          LOGF(info, "Calibration for %s estimator -> Data centrality bin %g-%g%%\n", estimator.Data(), centralityBins[irev], centralityBins[irev + 1]);
          LOGF(info, "Calibration for %s estimator -> MC multiplicity range %g-%g\n", estimator.Data(), projMC->GetBinLowEdge(endBinMc), projMC->GetBinLowEdge(startBinMc));
          LOGF(info, "Calibration for %s estimator -> <PV> data = %.4f Vs <PV> MC = %.4f (MC/Data = %.4f%%)\n", estimator.Data(), meanMult_Data, meanMult_MC, (meanMult_MC - meanMult_Data) * 100 / meanMult_Data);
          LOGF(info, "Calibration for %s estimator -> N entries Data = %g Vs N entries MC = %g\n", estimator.Data(), projData->Integral(), projMC->Integral());
          for (int ibin = 1; ibin <= h1dCalib->GetNbinsX(); ibin++) {
            if (ibin <= startBinMc && ibin >= endBinMc) {
              h1dCalib->SetBinContent(ibin, percentile_center[i]);
            }
          }
          startBinMc = endBinMc;
        } // End loop over centrality bins
      } else {
        for (int i = 0; i < nCentBins; i++) {
          int irev = nCentBins - i;
          percentile_center[i] = (centralityBins[irev - 1] + centralityBins[irev]) / 2;
          epercentile_center[i] = (centralityBins[irev] - centralityBins[irev - 1]) / 2;
        }
        int startBinMc = 1;
        for (int i = 0; i < nCentBins; i++) { // Loop over centrality bins
          // start from the end (from the low multiplicity collisions)
          int irev = nCentBins - i;
          TH1D* projData = h2dMultVsCent_Data->ProjectionY(Form("projData_%d", i), h2dMultVsCent_Data->GetXaxis()->FindBin(centralityBins[irev - 1] + 1e-5),
                                                           h2dMultVsCent_Data->GetXaxis()->FindBin(centralityBins[irev] - 1e-5));

          double meanMult_Data = projData->GetMean();

          double meanMult_MC = -1;
          double diffMultDataVsMC = 1e+09;
          int endBinMc = h1dCalib->GetNbinsX();
          for (int j = startBinMc; j <= h1dCalib->GetNbinsX(); j++) {
            // Loop over MC bins
            TH1D* projMC = h2dMultRecoVsMultGen_MC->ProjectionY("", startBinMc, endBinMc);
            int nEntries = projMC->Integral();
            double meanMC = projMC->GetMean();
            double ldiff = std::abs(meanMC - meanMult_Data);
            // std::cout << ldiff << std::endl;

            // if less than 100 entries for estimating the mean, do not consider it
            if (nEntries < baseOpts.minEntries)
              continue;

            // Find the minimum difference and corresponding MC bin
            if (ldiff < diffMultDataVsMC) {
              diffMultDataVsMC = ldiff;
              meanMult_MC = meanMC;
              endBinMc = j;
            }
          }
          if (i == nCentBins - 1) {
            endBinMc = h2dMultRecoVsMultGen_MC->GetNbinsX();
          }

          TH1D* projMC = h2dMultRecoVsMultGen_MC->ProjectionY(Form("projMC_%d", i), startBinMc, endBinMc);

          // Refill the per-class mean from the underlying (data/MC) multiplicity distributions, bin by
          // bin, so the TProfile accumulates real sum(w)/sum(wy)/sum(wy^2) rather than a precomputed mean.
          if (hPVData) {
            for (int ybin = 1; ybin <= projData->GetNbinsX(); ybin++) {
              double w = projData->GetBinContent(ybin);
              if (w > 0) {
                hPVData->Fill(percentile_center[i], projData->GetBinCenter(ybin), w);
              }
            }
          }
          if (hPVMC) {
            for (int ybin = 1; ybin <= projMC->GetNbinsX(); ybin++) {
              double w = projMC->GetBinContent(ybin);
              if (w > 0) {
                hPVMC->Fill(percentile_center[i], projMC->GetBinCenter(ybin), w);
              }
            }
          }

          LOGF(info, "Calibration for %s estimator -> Data centrality bin %g-%g%%\n", estimator.Data(), centralityBins[irev - 1], centralityBins[irev]);
          LOGF(info, "Calibration for %s estimator -> MC multiplicity range %g-%g\n", estimator.Data(), projMC->GetBinLowEdge(startBinMc), projMC->GetBinLowEdge(endBinMc + 1));
          LOGF(info, "Calibration for %s estimator -> <PV> data = %.4f Vs <PV> MC = %.4f\n", estimator.Data(), meanMult_Data, meanMult_MC);
          LOGF(info, "Calibration for %s estimator -> N entries Data = %g Vs N entries MC = %g\n", estimator.Data(), projData->Integral(), projMC->Integral());
          for (int ibin = 1; ibin <= h1dCalib->GetNbinsX(); ibin++) {
            if (ibin <= endBinMc && ibin >= startBinMc) {
              h1dCalib->SetBinContent(ibin, percentile_center[i]);
            }
          }

          startBinMc = endBinMc + 1;
        } // End loop over centrality bins
      }

      LOGF(info, "Calibration for %s estimator -> Done!\n", estimator.Data());

      return h1dCalib;
    };
    return CalibMC(name);
  }

  // Loads the calibration object(s) for a single estimator: either the class-averaged TH1D (statistical
  // assignment) or, in candidate-by-candidate mode, the raw joint TH2D used to sample P(centrality|genMult).
  void loadEstimatorCalibration(int idx, const char* name, TH1D*& h1dOut)
  {
    if (!baseOpts.mEnabledTables[idx]) {
      return;
    }
    if (baseOpts.assignCentralityPerCandidate) {
      // x = reco centrality of the matched MC collision, y = generated multiplicity (|eta| < 0.5)
      // filled directly by the calibration-producing task (e.g. centralityQa.cxx); used as-is, no
      // mean-matching required since it is already a genuine joint distribution.
      h2dCentVsGenMult[idx] = getHist<TH2D>(Form("hGenMultEta05VsCentrality%s", name));
      return;
    }
    h1dOut = baseOpts.recalibrateCentrality ? extractCentralityCalibration(idx, name, baseOpts.recalibrateMode) : getHist<TH1D>(Form("h1d%s", name));
  }

  template <typename TCCDB, typename TBCs>
  bool initCCDB(TCCDB& ccdb, TBCs const& bcs)
  {
    if (!bcs.size()) {
      LOGF(warn, "No BC found, skipping this DF.");
      return false; // signal to skip this DF
    }

    if (mRunNumber == bcs.iteratorAt(0).runNumber()) {
      return true;
    }
    // mark this run as configured
    mRunNumber = bcs.iteratorAt(0).runNumber();

    MCCentralityCalibObjects = ccdb->template getForRun<TList>(baseOpts.path, mRunNumber);
    if (MCCentralityCalibObjects) {
      LOGF(info, "loaded TList with this many objects: %i", MCCentralityCalibObjects->GetEntries());
    }

    loadEstimatorCalibration(kFV0A, "FV0A", h1dFV0A);
    loadEstimatorCalibration(kFT0M, "FT0M", h1dFT0M);
    loadEstimatorCalibration(kFT0A, "FT0A", h1dFT0A);
    loadEstimatorCalibration(kFT0C, "FT0C", h1dFT0C);
    loadEstimatorCalibration(kFT0CVariant1, "FT0CVariant1", h1dFT0CVariant1);
    loadEstimatorCalibration(kFT0CVariant2, "FT0CVariant2", h1dFT0CVariant2);
    loadEstimatorCalibration(kFDDM, "FDDM", h1dFDDM);
    loadEstimatorCalibration(kNTPV, "NTPV", h1dNTPV);
    loadEstimatorCalibration(kNGlobal, "NGlobal", h1dNGlobal);
    if (baseOpts.mEnabledTables[kMFT]) {
      // to be added later (multMCMFT not yet available), kept as a dedicated block on purpose
    }
    LOG(info) << "Fully configured for run: " << mRunNumber;

    return true;
  }

  template <int tableIndex>
  void fillHistograms(o2::framework::HistogramRegistry& histos, double percentile, double multiplicity, double genMultEta05)
  {
    histos.fill(HIST(DirList[tableIndex]) + HIST("/percentile"), percentile);
    histos.fill(HIST(DirList[tableIndex]) + HIST("/percentilevsMult"), percentile, multiplicity);
    if (baseOpts.assignCentralityPerCandidate) {
      // QA of the sampling itself: sampled percentile vs the generated multiplicity it was actually keyed
      // on (registered by registerCandidateQAHistos()), so it can be compared to the input calibration.
      histos.fill(HIST(DirList[tableIndex]) + HIST("/percentileVsGenMultEta05"), percentile, genMultEta05);
    }
  }

  // Combines the run number, a per-DF collision sequence number and the estimator index into a single RNG
  // seed, so that a given generated collision always draws the same sampled centrality for a given
  // estimator, regardless of processing order (note: the sequence number is only unique within one DF,
  // it is not a persistent table index).
  // This is a simple hash-combine, not a cryptographic hash: the strides just keep the three fields from
  // landing on the same seed value for the collision counts typically seen in one run, they don't
  // guarantee it. A rare accidental collision only means two unrelated collisions share a random draw,
  // which has no systematic effect on the sampled distribution.
  ULong64_t computeSamplingSeed(int64_t collisionIndex, int tableIdx) const
  {
    static constexpr ULong64_t kCollisionStride = 131ull; // > kNestimators, so tableIdx can't alias into the collision term
    static constexpr ULong64_t kRunStride = 1000003ull;   // prime, well above kCollisionStride * (typical collisions per DF)
    return baseOpts.centralitySamplingSeed.value + kRunStride * static_cast<ULong64_t>(mRunNumber) +
           kCollisionStride * static_cast<ULong64_t>(collisionIndex) + static_cast<ULong64_t>(tableIdx);
  }

  //__________________________________________________
  template <typename TCCDB, typename THistoRegistry, typename TBCs, typename TProducts>
  void dataProcess(TCCDB& ccdb, THistoRegistry& histos, TBCs const& bcs, const HasMcMults auto& mcCollisions, TProducts& products)
  {
    if (nEnabledTables == 0) {
      return; // fully suppressed
    }

    if (!initCCDB(ccdb, bcs))
      return;

    // Local per-DF sequence number used only to seed the candidate-by-candidate sampling below: not every
    // joined mcCollisions type carries an index column (globalIndex()), so we cannot rely on that instead.
    int64_t mcCollisionCounter = 0;
    for (const auto& mcCollision : mcCollisions) {
      const double nFV0A = mcCollision.multMCFV0A();
      const double nFT0A = mcCollision.multMCFT0A();
      const double nFT0C = mcCollision.multMCFT0C();
      const double nFDDA = mcCollision.multMCFDDA();
      const double nFDDC = mcCollision.multMCFDDC();
      const double nGlobal = mcCollision.multMCNParticlesEta08();
      const double nGenMultEta05 = mcCollision.multMCNParticlesEta05();
      // const double nMFT = mcCollision.multMCMFT(); // to be added later
      const double nFT0M = nFT0A + nFT0C;
      const double nFDDM = nFDDA + nFDDC;

      auto populateTable = [&](auto& table, TH1D* h1dMultCalib, TH2D* h2dCentCalib, double multiplicity, auto tableIndex) {
        double percentile = 105.0f;
        if (baseOpts.mEnabledTables[tableIndex]) {
          if (baseOpts.assignCentralityPerCandidate && h2dCentCalib) {
            // Sample P(centrality | genMult) for this generated collision: slice h2dCentCalib at the
            // y-bin (generated mult |eta|<0.5) matching this collision and draw from that 1-D distribution.
            fRandomSampler.SetSeed(computeSamplingSeed(mcCollisionCounter, tableIndex.value));
            int yBin = h2dCentCalib->GetYaxis()->FindBin(nGenMultEta05);
            std::unique_ptr<TH1D> pdf{h2dCentCalib->ProjectionX(Form("pdfCent_%d", tableIndex.value), yBin, yBin)};
            if (pdf->Integral() > 0) {
              percentile = pdf->GetRandom(&fRandomSampler);
            }
            fillHistograms<tableIndex.value>(histos, percentile, multiplicity, nGenMultEta05);
          } else if (h1dMultCalib) {
            int bin = h1dMultCalib->FindBin(multiplicity);
            percentile = h1dMultCalib->GetBinContent(bin);
            fillHistograms<tableIndex.value>(histos, percentile, multiplicity, nGenMultEta05);
          }
          table(percentile);
        }
        return percentile;
      };

      populateTable(products.centFV0A, h1dFV0A, h2dCentVsGenMult[kFV0A], nFV0A, std::integral_constant<int, kFV0A>{});
      populateTable(products.centFT0M, h1dFT0M, h2dCentVsGenMult[kFT0M], nFT0M, std::integral_constant<int, kFT0M>{});
      populateTable(products.centFT0A, h1dFT0A, h2dCentVsGenMult[kFT0A], nFT0A, std::integral_constant<int, kFT0A>{});
      populateTable(products.centFT0C, h1dFT0C, h2dCentVsGenMult[kFT0C], nFT0C, std::integral_constant<int, kFT0C>{});
      populateTable(products.centFT0CVariant1, h1dFT0CVariant1, h2dCentVsGenMult[kFT0CVariant1], nFT0C, std::integral_constant<int, kFT0CVariant1>{});
      populateTable(products.centFT0CVariant2, h1dFT0CVariant2, h2dCentVsGenMult[kFT0CVariant2], nFT0C, std::integral_constant<int, kFT0CVariant2>{});
      populateTable(products.centFDDM, h1dFDDM, h2dCentVsGenMult[kFDDM], nFDDM, std::integral_constant<int, kFDDM>{});
      populateTable(products.centNTPV, h1dNTPV, h2dCentVsGenMult[kNTPV], nGlobal, std::integral_constant<int, kNTPV>{});
      populateTable(products.centNGlobal, h1dNGlobal, h2dCentVsGenMult[kNGlobal], nGlobal, std::integral_constant<int, kNGlobal>{});
      // populateTable(products.centMFT, h1dMFT, h2dCentVsGenMult[kMFT], nMFT, std::integral_constant<int, kMFT>{}); // to be added later

      mcCollisionCounter++;
    }
  }
}; // end mcCentralityModule

} // namespace mccentrality
} // namespace pwglf
} // namespace o2

#endif // PWGLF_UTILS_MCCENTRALITYMODULE_
