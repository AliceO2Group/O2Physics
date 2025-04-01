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

#include <CCDB/BasicCCDBManager.h>
#include <DataFormatsParameters/GRPObject.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <algorithm>
#include <numeric>
#include <vector>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/HistogramRegistry.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"

#include "GFWPowerArray.h"
#include "GFW.h"
#include "GFWCumulant.h"
#include "FlowContainer.h"
#include "FlowPtContainer.h"
#include "GFWConfig.h"
#include "GFWWeights.h"
#include <TProfile.h>
#include <TRandom3.h>
#include <TF1.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::analysis;

#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};
namespace o2::analysis::flowanalysis
{
std::vector<double> ptbinning = {
  0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55,
  0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95,
  1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9,
  2, 2.2, 2.4, 2.6, 2.8, 3, 3.5, 4, 5, 6, 8, 10};
float ptpoilow = 0.2, ptpoiup = 10.0;
float ptreflow = 0.2, ptrefup = 3.0;
float ptlow = 0.2, ptup = 10.0;
int etabins = 16;
float etalow = -0.8, etaup = 0.8;
int vtxZbins = 20;
float vtxZlow = -10.0, vtxZup = 10.0;
int phibins = 72;
float philow = 0.0;
float phiup = constants::math::TwoPI;
int nchbins = 300;
float nchlow = 0;
float nchup = 3000;
std::vector<double> centbinning(90);
int nBootstrap = 10;
GFWRegions regions;
GFWCorrConfigs configs;
int eventBin = 0;
} // namespace o2::analysis::flowanalysis

using namespace o2::analysis::flowanalysis;

struct flowAnalysisGF {

  O2_DEFINE_CONFIGURABLE(cfgNbootstrap, int, 10, "Number of subsamples")
  O2_DEFINE_CONFIGURABLE(cfgMpar, int, 8, "Highest order of pt-pt correlations")
  O2_DEFINE_CONFIGURABLE(cfgUseNch, bool, false, "Do correlations as function of Nch")
  O2_DEFINE_CONFIGURABLE(cfgFillWeights, bool, false, "Fill NUA weights")
  O2_DEFINE_CONFIGURABLE(cfgFillQA, bool, false, "Fill QA histograms")
  O2_DEFINE_CONFIGURABLE(cfgUseAdditionalEventCut, bool, true, "Use additional event cut on mult correlations")
  O2_DEFINE_CONFIGURABLE(cfgUseAdditionalTrackCut, bool, true, "Use additional track cut on phi")
  O2_DEFINE_CONFIGURABLE(cfgEfficiency, std::string, "", "CCDB path to efficiency object")
  O2_DEFINE_CONFIGURABLE(cfgAcceptance, std::string, "", "CCDB path to acceptance object")
  O2_DEFINE_CONFIGURABLE(cfgDCAxy, float, 0.2, "Cut on DCA in the transverse direction (cm)");
  O2_DEFINE_CONFIGURABLE(cfgDCAz, float, 2, "Cut on DCA in the longitudinal direction (cm)");
  O2_DEFINE_CONFIGURABLE(cfgNcls, float, 70, "Cut on number of TPC clusters found");
  O2_DEFINE_CONFIGURABLE(cfgPtmin, float, 0.2, "minimum pt (GeV/c)");
  O2_DEFINE_CONFIGURABLE(cfgPtmax, float, 10, "maximum pt (GeV/c)");
  O2_DEFINE_CONFIGURABLE(cfgEta, float, 0.8, "eta cut");
  O2_DEFINE_CONFIGURABLE(cfgVtxZ, float, 10, "vertex cut (cm)");
  O2_DEFINE_CONFIGURABLE(cfgMagField, float, 99999, "Configurable magnetic field; default CCDB will be queried");
  // added
  O2_DEFINE_CONFIGURABLE(cfgDoubleTrackFunction, bool, true, "Include track cut at low pt");
  O2_DEFINE_CONFIGURABLE(cfgTrackCutSize, float, 0.06, "Spread of track cut");
  O2_DEFINE_CONFIGURABLE(cfgMaxOccupancy, int, 500, "Maximum occupancy of selected events");
  O2_DEFINE_CONFIGURABLE(cfgNoSameBunchPileupCut, bool, true, "kNoSameBunchPileupCut");
  O2_DEFINE_CONFIGURABLE(cfgIsGoodZvtxFT0vsPV, bool, true, "kIsGoodZvtxFT0vsPV");
  O2_DEFINE_CONFIGURABLE(cfgNoCollInTimeRangeStandard, bool, true, "kNoCollInTimeRangeStandard");
  O2_DEFINE_CONFIGURABLE(cfgDoOccupancySel, bool, true, "Bool for event selection on detector occupancy");
  O2_DEFINE_CONFIGURABLE(cfgMultCut, bool, true, "Use additional evenr cut on mult correlations");
  O2_DEFINE_CONFIGURABLE(cfgTVXinTRD, bool, true, "Use kTVXinTRD (reject TRD triggered events)");
  O2_DEFINE_CONFIGURABLE(cfgIsVertexITSTPC, bool, true, "Selects collisions with at least one ITS-TPC track");

  Configurable<GFWBinningCuts> cfgGFWBinning{"cfgGFWBinning", {40, 16, 72, 300, 0, 3000, 0.2, 10.0, 0.2, 3.0, {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3, 3.25, 3.5, 3.75, 4, 4.5, 5, 5.5, 6, 7, 8, 9, 10}, {0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90}}, "Configuration for binning"};
  Configurable<GFWRegions> cfgRegions{"cfgRegions", {{"refN", "refP", "refFull", "poiN", "poiP", "poiFull", "olN", "olP", "olFull"}, {-0.8, 0.5, -0.8, -0.8, -0.8, 0.5, -0.8, -0.8}, {-0.5, 0.8, 0.8, -0.5, 0.8, 0.8, -0.5}, {0, 0, 0, 1, 1, 1, 1, 1, 1}, {1, 1, 1, 2, 2, 2, 4, 4, 4}}, "Configurations for GFW regions"};

  Configurable<GFWCorrConfigs> cfgCorrConfig{"cfgCorrConfig", {{"refP {2} refN {-2}", "poiP refP | olP {2} refN {-2}", "refP {3} refN {-3}", "refP {4} refN {-4}", "refFull {2 -2}", "poiFull refFull | olFull {2 -2}", "refFull {2 2 -2 -2}", "poiFull refFull | olFull {2 2 -2 -2}"}, {"ChGap22", "ChGap22", "ChGap32", "ChGap42", "ChFull22", "ChFull22", "ChFull24", "ChFull24"}, {0, 0, 0, 0, 0, 0, 0, 0}, {0, 1, 0, 0, 0, 1, 0, 1}}, "Configurations for each correlation to calculate"};

  //  Connect to ccdb
  Service<ccdb::BasicCCDBManager> ccdb;

  struct Config {
    TH1D* mEfficiency = nullptr;
    GFWWeights* mAcceptance = nullptr;
    bool correctionsLoaded = false;
  } cfg;

  // Define output
  OutputObj<FlowContainer> fFC{FlowContainer("FlowContainer")};
  OutputObj<FlowPtContainer> fFCpt{FlowPtContainer("FlowPtContainer")};
  OutputObj<FlowContainer> fFC_gen{FlowContainer("FlowContainer_gen")};
  OutputObj<GFWWeights> fWeights{GFWWeights("weights")};
  HistogramRegistry registry{"registry"};

  // define global variables
  GFW* fGFW = new GFW();
  std::vector<GFW::CorrConfig> corrconfigs;
  TRandom3* fRndm = new TRandom3(0);
  TAxis* fPtAxis;

  // Event selection cuts - Alex
  TF1* fPhiCutLow = nullptr;
  TF1* fPhiCutHigh = nullptr;
  // for double track cut
  TF1* fPhiCutLow2 = nullptr;
  TF1* fPhiCutHigh2 = nullptr;

  TF1* fMultPVCutLow = nullptr;
  TF1* fMultPVCutHigh = nullptr;
  TF1* fMultCutLow = nullptr;
  TF1* fMultCutHigh = nullptr;
  TF1* fMultMultPVCut = nullptr;

  void init(InitContext const&)
  {
    LOGF(info, "flowAnalysisGF::init()");
    regions.SetNames(cfgRegions->GetNames());
    regions.SetEtaMin(cfgRegions->GetEtaMin());
    regions.SetEtaMax(cfgRegions->GetEtaMax());
    regions.SetpTDifs(cfgRegions->GetpTDifs());
    regions.SetBitmasks(cfgRegions->GetBitmasks());
    configs.SetCorrs(cfgCorrConfig->GetCorrs());
    configs.SetHeads(cfgCorrConfig->GetHeads());
    configs.SetpTDifs(cfgCorrConfig->GetpTDifs());
    configs.SetpTCorrMasks(cfgCorrConfig->GetpTCorrMasks());

    regions.Print();
    configs.Print();

    ptbinning = cfgGFWBinning->GetPtBinning();
    ptpoilow = cfgGFWBinning->GetPtPOImin();
    ptpoiup = cfgGFWBinning->GetPtPOImax();
    ptreflow = cfgGFWBinning->GetPtRefMin();
    ptrefup = cfgGFWBinning->GetPtRefMax();
    ptlow = cfgPtmin;
    ptup = cfgPtmax;
    etabins = cfgGFWBinning->GetEtaBins();
    vtxZbins = cfgGFWBinning->GetVtxZbins();
    phibins = cfgGFWBinning->GetPhiBins();
    philow = 0.0f;
    phiup = constants::math::TwoPI;
    nchbins = cfgGFWBinning->GetNchBins();
    nchlow = cfgGFWBinning->GetNchMin();
    nchup = cfgGFWBinning->GetNchMax();
    centbinning = cfgGFWBinning->GetCentBinning();
    cfgGFWBinning->Print();

    AxisSpec phiAxis = {phibins, philow, phiup, "#phi"};
    AxisSpec phiModAxis = {100, 0, constants::math::PI / 9, "fmod(#varphi,#pi/9)"};
    AxisSpec etaAxis = {etabins, -cfgEta, cfgEta, "#eta"};
    AxisSpec vtxAxis = {vtxZbins, -cfgVtxZ, cfgVtxZ, "Vtx_{z} (cm)"};
    AxisSpec ptAxis = {ptbinning, "#it{p}_{T} GeV/#it{c}"};
    AxisSpec centAxis = {centbinning, "Centrality (%)"};
    std::vector<double> nchbinning;
    int nchskip = (nchup - nchlow) / nchbins;
    for (int i = 0; i <= nchbins; ++i) {
      nchbinning.push_back(nchskip * i + nchlow + 0.5);
    }

    AxisSpec nchAxis = {nchbinning, "N_{ch}"};
    AxisSpec t0cAxis = {70, 0, 70000, "N_{ch} (T0C)"};
    AxisSpec t0aAxis = {200, 0, 200, "N_{ch}"};
    AxisSpec multpvAxis = {4000, 0, 4000, "N_{ch} (PV)"};
    AxisSpec multAxis = (cfgUseNch) ? nchAxis : centAxis;
    AxisSpec dcaZAxis = {200, -2, 2, "DCA_{z} (cm)"};
    AxisSpec dcaXYAxis = {200, -1, 1, "DCA_{xy} (cm)"};
    AxisSpec tpcAxis = {200, 0, 200, "TPC Clusters"};
    AxisSpec tpcAxisCl = {10, 0, 10, "TPC Crossed rows / findable Clusters"};
    AxisSpec dEdxAxis = {400, 0, 200, "TPC: dE/dx"};
    AxisSpec itsClAxis = {10, 0, 10, "Number of ITS clusters"};

    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    int64_t now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    ccdb->setCreatedNotAfter(now);

    int ptbins = ptbinning.size() - 1;
    fPtAxis = new TAxis(ptbins, &ptbinning[0]);

    if (cfgFillWeights) {
      fWeights->setPtBins(ptbins, &ptbinning[0]);
      fWeights->init(true, false);
    }

    if (doprocessMCGen) {
      registry.add("trackMCGen/pt_gen", "", {HistType::kTH1D, {ptAxis}});
      registry.add("trackMCGen/phi_eta_vtxZ_gen", "", {HistType::kTH3D, {phiAxis, etaAxis, vtxAxis}});
    }
    if (doprocessMCReco || doprocessData || doprocessRun2) {
      registry.add("trackQA/phi_eta_vtxZ", "", {HistType::kTH3D, {phiAxis, etaAxis, vtxAxis}});
      registry.add("trackQA/pt_dcaXY_dcaZ", "", {HistType::kTH3D, {ptAxis, dcaXYAxis, dcaZAxis}});
      registry.add("trackQA/pt_phi_bef", "", {HistType::kTH2D, {ptAxis, phiModAxis}});
      registry.add("trackQA/pt_phi_aft", "", {HistType::kTH2D, {ptAxis, phiModAxis}});
      registry.add("trackQA/phi_eta_vtxZ_corrected", "", {HistType::kTH3D, {phiAxis, etaAxis, vtxAxis}});
      registry.add("trackQA/pt_reco", "", {HistType::kTH1D, {ptAxis}});

      registry.add("hEventCount", "Number of Event;; Count", {HistType::kTH1D, {{10, 0, 10}}});
      registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(1, "Filtered event");
      registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(2, "sel8");
      registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(3, "occupancy");
      registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(4, "kTVXinTRD");
      registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(5, "kNoSameBunchPileup");
      registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(6, "kIsGoodZvtxFT0vsPV");
      registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(7, "kNoCollInTimeRangeStandard");
      registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(8, "kIsVertexITSTPC");
      registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(9, "after Mult cuts");
      registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(10, "corr + cent");

      registry.add("eventQA/nITS_nTPC_Tracks", "", {HistType::kTH2D, {nchAxis, nchAxis}});

      // track QA before and after selections
      registry.add("trackQA/before/CrossedTPCRows_pt", "", {HistType::kTH2D, {ptAxis, tpcAxis}});
      registry.add("trackQA/before/NumberSharedClustersTPC_pt", "", {HistType::kTH2D, {ptAxis, tpcAxis}});
      registry.add("trackQA/before/ITSClusters_pt", "", {HistType::kTH2D, {ptAxis, itsClAxis}});
      registry.add("trackQA/before/DCAxy_pt", "", {HistType::kTH2D, {ptAxis, dcaXYAxis}});
      registry.add("trackQA/before/DCAz_pt", "", {HistType::kTH2D, {ptAxis, dcaZAxis}});
      registry.add("trackQA/before/tpcSignal_pt", "", {HistType::kTH2D, {ptAxis, dEdxAxis}});
      registry.add("trackQA/before/phi_pt", "", {HistType::kTH2D, {ptAxis, phiAxis}});
      registry.add("trackQA/before/phi_eta", "", {HistType::kTH2D, {etaAxis, phiAxis}});
      registry.add("trackQA/before/tpcCrossedRowsOverFindableCls_pt", "", {HistType::kTH2D, {ptAxis, tpcAxisCl}});
      registry.addClone("trackQA/before/", "trackQA/after/");

      // track QA after selections
      registry.add("eventQA/before/globalTracks_centT0C", "", {HistType::kTH2D, {centAxis, nchAxis}});
      registry.add("eventQA/before/PVTracks_centT0C", "", {HistType::kTH2D, {centAxis, multpvAxis}});
      registry.add("eventQA/before/globalTracks_PVTracks", "", {HistType::kTH2D, {multpvAxis, nchAxis}});
      registry.add("eventQA/before/globalTracks_multT0A", "", {HistType::kTH2D, {t0aAxis, nchAxis}});
      registry.add("eventQA/before/globalTracks_multV0A", "", {HistType::kTH2D, {t0aAxis, nchAxis}});
      registry.add("eventQA/before/multV0A_multT0A", "", {HistType::kTH2D, {t0aAxis, t0aAxis}});
      registry.add("eventQA/before/multT0C_centT0C", "", {HistType::kTH2D, {centAxis, t0cAxis}});
      registry.addClone("eventQA/before/", "eventQA/after/");
    }

    if (regions.GetSize() < 0)
      LOGF(error, "Configuration contains vectors of different size - check the GFWRegions configurable");
    for (auto i(0); i < regions.GetSize(); ++i) {
      fGFW->AddRegion(regions.GetNames()[i], regions.GetEtaMin()[i], regions.GetEtaMax()[i], (regions.GetpTDifs()[i]) ? ptbins + 1 : 1, regions.GetBitmasks()[i]);
    }
    for (auto i = 0; i < configs.GetSize(); ++i) {
      corrconfigs.push_back(fGFW->GetCorrelatorConfig(configs.GetCorrs()[i], configs.GetHeads()[i], configs.GetpTDifs()[i]));
    }
    if (corrconfigs.empty())
      LOGF(error, "Configuration contains vectors of different size - check the GFWCorrConfig configurable");
    fGFW->CreateRegions();
    TObjArray* oba = new TObjArray();
    AddConfigObjectsToObjArray(oba, corrconfigs);

    fFC->SetName("FlowContainer");
    fFC->SetXAxis(fPtAxis);
    fFC->Initialize(oba, multAxis, cfgNbootstrap);

    delete oba;
    fFCpt->initialise(multAxis, cfgMpar, configs, cfgNbootstrap);

    // Event selection - Alex
    if (cfgUseAdditionalEventCut) {
      // zzh pas3
      fMultPVCutLow = new TF1("fMultPVCutLow", "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x - 3.5*([5]+[6]*x+[7]*x*x+[8]*x*x*x+[9]*x*x*x*x)", 0, 100);
      fMultPVCutLow->SetParameters(3307.55, -122.1, 1.97653, -0.0172405, 6.57892e-05, 147.955, -2.40658, 0.00626944, 7.19406e-05, -3.92605e-07);

      fMultPVCutHigh = new TF1("fMultPVCutHigh", "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x + 3.5*([5]+[6]*x+[7]*x*x+[8]*x*x*x+[9]*x*x*x*x)", 0, 100);
      fMultPVCutHigh->SetParameters(3307.55, -122.1, 1.97653, -0.0172405, 6.57892e-05, 147.955, -2.40658, 0.00626944, 7.19406e-05, -3.92605e-07);

      fMultCutLow = new TF1("fMultCutLow", "[0]+[1]*x+[2]*x*x+[3]*x*x*x - 2.*([4]+[5]*x+[6]*x*x+[7]*x*x*x+[8]*x*x*x*x)", 0, 100);
      fMultCutLow->SetParameters(1770.02, -49.5537, 0.460941, -0.00140622, 105.477, -1.58301, 0.011655, -0.000190804, 1.36003e-06);

      fMultCutHigh = new TF1("fMultCutHigh", "[0]+[1]*x+[2]*x*x+[3]*x*x*x + 3.*([4]+[5]*x+[6]*x*x+[7]*x*x*x+[8]*x*x*x*x)", 0, 100);
      fMultCutHigh->SetParameters(1770.02, -49.5537, 0.460941, -0.00140622, 105.477, -1.58301, 0.011655, -0.000190804, 1.36003e-06);
    }

    if (cfgUseAdditionalTrackCut) {
      // 0.06
      fPhiCutLow = new TF1("fPhiCutLow", TString::Format("0.06/x+pi/18.0-%f", static_cast<float>(cfgTrackCutSize)), 0, 100);
      fPhiCutHigh = new TF1("fPhiCutHigh", TString::Format("0.1/x+pi/18.0+%f", static_cast<float>(cfgTrackCutSize)), 0, 100);

      if (cfgDoubleTrackFunction) {
        fPhiCutLow2 = new TF1("fPhiCutLow2", TString::Format("0.06/x+pi/18.0-%f-pi/9", static_cast<float>(cfgTrackCutSize)), 0, 100); // for the lower part of the gap (pt<1)
        fPhiCutHigh2 = new TF1("fPhiCutHigh2", TString::Format("0.1/x+pi/18.0+%f-pi/9", static_cast<float>(cfgTrackCutSize)), 0, 100);
      }
    }
  }

  static constexpr std::string_view moment[] = {
    "before/",
    "after/"};

  enum QAtime {
    kBefore,
    kAfter
  };

  void AddConfigObjectsToObjArray(TObjArray* oba, const std::vector<GFW::CorrConfig>& configs)
  {
    for (auto it = configs.begin(); it != configs.end(); ++it) {
      if (it->pTDif) {
        std::string suffix = "_ptDiff";
        for (auto i = 0; i < fPtAxis->GetNbins(); ++i) {
          std::string index = Form("_pt_%i", i + 1);
          oba->Add(new TNamed(it->Head.c_str() + index, it->Head.c_str() + suffix));
        }
      } else {
        oba->Add(new TNamed(it->Head.c_str(), it->Head.c_str()));
      }
    }
  }

  int getMagneticField(uint64_t timestamp)
  {
    // TODO done only once (and not per run). Will be replaced by CCDBConfigurable
    // static o2::parameters::GRPObject* grpo = nullptr;
    static o2::parameters::GRPMagField* grpo = nullptr;
    if (grpo == nullptr) {
      grpo = ccdb->getForTimeStamp<o2::parameters::GRPMagField>("GLO/Config/GRPMagField", timestamp);
      if (grpo == nullptr) {
        LOGF(fatal, "GRP object not found for timestamp %llu", timestamp);
        return 0;
      }
      LOGF(info, "Retrieved GRP for timestamp %llu with magnetic field of %d kG", timestamp, grpo->getNominalL3Field());
    }
    return grpo->getNominalL3Field();
  }

  void loadCorrections(uint64_t timestamp)
  {
    if (cfg.correctionsLoaded)
      return;
    if (cfgAcceptance.value.empty() == false) {
      cfg.mAcceptance = ccdb->getForTimeStamp<GFWWeights>(cfgAcceptance, timestamp);
      if (cfg.mAcceptance)
        LOGF(info, "Loaded acceptance weights from %s (%p)", cfgAcceptance.value.c_str(), (void*)cfg.mAcceptance);
      else
        LOGF(warning, "Could not load acceptance weights from %s (%p)", cfgAcceptance.value.c_str(), (void*)cfg.mAcceptance);
    }
    if (cfgEfficiency.value.empty() == false) {
      cfg.mEfficiency = ccdb->getForTimeStamp<TH1D>(cfgEfficiency, timestamp);
      if (cfg.mEfficiency == nullptr) {
        LOGF(fatal, "Could not load efficiency histogram for trigger particles from %s", cfgEfficiency.value.c_str());
      }
      LOGF(info, "Loaded efficiency histogram from %s (%p)", cfgEfficiency.value.c_str(), (void*)cfg.mEfficiency);
    }
    cfg.correctionsLoaded = true;
  }

  bool setCurrentParticleWeights(float& weight_nue, float& weight_nua, const float& phi, const float& eta, const float& pt, const float& vtxz)
  {
    float eff = 1.;
    if (cfg.mEfficiency)
      eff = cfg.mEfficiency->GetBinContent(cfg.mEfficiency->FindBin(pt));
    else
      eff = 1.0;
    if (eff == 0)
      return false;
    weight_nue = 1. / eff;
    if (cfg.mAcceptance)
      weight_nua = cfg.mAcceptance->getNUA(phi, eta, vtxz);
    else
      weight_nua = 1;
    return true;
  }

  template <typename TCollision>
  bool eventSelected(TCollision collision, const int& multTrk, const float& centrality)
  {
    if (cfgTVXinTRD) {
      if (collision.alias_bit(kTVXinTRD)) {
        // TRD triggered
        // "CMTVX-B-NOPF-TRD,minbias_TVX"
        return 0;
      }
      registry.fill(HIST("hEventCount"), 3.5);
    }

    if (cfgNoSameBunchPileupCut) {
      if (!collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
        // rejects collisions which are associated with the same "found-by-T0" bunch crossing
        // https://indico.cern.ch/event/1396220/#1-event-selection-with-its-rof
        return 0;
      }
      registry.fill(HIST("hEventCount"), 4.5);
    }
    if (cfgIsGoodZvtxFT0vsPV) {
      if (!collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
        // removes collisions with large differences between z of PV by tracks and z of PV from FT0 A-C time difference
        // use this cut at low multiplicities with caution
        return 0;
      }
      registry.fill(HIST("hEventCount"), 5.5);
    }
    if (cfgNoCollInTimeRangeStandard) {
      if (!collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
        //  Rejection of the collisions which have other events nearby
        return 0;
      }
      registry.fill(HIST("hEventCount"), 6.5);
    }

    if (cfgIsVertexITSTPC) {
      if (!collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) {
        // selects collisions with at least one ITS-TPC track, and thus rejects vertices built from ITS-only tracks
        return 0;
      }
      registry.fill(HIST("hEventCount"), 7.5);
    }

    float vtxz = -999;
    if (collision.numContrib() > 1) {
      vtxz = collision.posZ();
      float zRes = TMath::Sqrt(collision.covZZ());
      if (zRes > 0.25 && collision.numContrib() < 20)
        vtxz = -999;
    }

    auto multNTracksPV = collision.multNTracksPV();

    if (vtxz > vtxZup || vtxz < vtxZlow)
      return 0;

    if (cfgMultCut) {
      if (multNTracksPV < fMultPVCutLow->Eval(centrality))
        return 0;
      if (multNTracksPV > fMultPVCutHigh->Eval(centrality))
        return 0;
      if (multTrk < fMultCutLow->Eval(centrality))
        return 0;
      if (multTrk > fMultCutHigh->Eval(centrality))
        return 0;
      registry.fill(HIST("hEventCount"), 8.5);
    }

    return 1;
  }

  template <typename TTrack>
  bool trackSelected(TTrack track, const int& field)
  {
    double phimodn = track.phi();
    if (field < 0) // for negative polarity field
      phimodn = TMath::TwoPi() - phimodn;
    if (track.sign() < 0) // for negative charge
      phimodn = TMath::TwoPi() - phimodn;
    if (phimodn < 0)
      LOGF(debug, "phi < 0: %g", phimodn);

    phimodn += TMath::Pi() / 18.0; // to center gap in the middle
    phimodn = fmod(phimodn, TMath::Pi() / 9.0);
    registry.fill(HIST("trackQA/pt_phi_bef"), track.pt(), phimodn);
    if (phimodn < fPhiCutHigh->Eval(track.pt()) && phimodn > fPhiCutLow->Eval(track.pt()))
      return false; // reject track
    if (cfgDoubleTrackFunction) {
      if (phimodn < fPhiCutHigh2->Eval(track.pt()) && phimodn > fPhiCutLow2->Eval(track.pt()))
        return false; // reject track
    }
    registry.fill(HIST("trackQA/pt_phi_aft"), track.pt(), phimodn);
    return true;
  }

  enum datatype {
    kReco,
    kGen
  };

  void FillOutputContainers(datatype dt, const float& centmult, const double& rndm)
  {
    fFCpt->calculateCorrelations();
    fFCpt->fillPtProfiles(centmult, rndm);
    fFCpt->fillCMProfiles(centmult, rndm);
    for (uint l_ind = 0; l_ind < corrconfigs.size(); ++l_ind) {
      auto dnx = fGFW->Calculate(corrconfigs.at(l_ind), 0, kTRUE).real();
      if (dnx == 0)
        continue;
      if (!corrconfigs.at(l_ind).pTDif) {
        auto val = fGFW->Calculate(corrconfigs.at(l_ind), 0, kFALSE).real() / dnx;
        if (TMath::Abs(val) < 1) {
          (dt == kGen) ? fFC_gen->FillProfile(corrconfigs.at(l_ind).Head.c_str(), centmult, val, dnx, rndm) : fFC->FillProfile(corrconfigs.at(l_ind).Head.c_str(), centmult, val, dnx, rndm);
          fFCpt->fillVnPtProfiles(centmult, val, dnx, rndm, configs.GetpTCorrMasks()[l_ind]);
        }
        continue;
      }
      for (Int_t i = 1; i <= fPtAxis->GetNbins(); i++) {
        dnx = fGFW->Calculate(corrconfigs.at(l_ind), i - 1, kTRUE).real();
        if (dnx == 0)
          continue;
        auto val = fGFW->Calculate(corrconfigs.at(l_ind), i - 1, kFALSE).real() / dnx;
        if (TMath::Abs(val) < 1)
          (dt == kGen) ? fFC_gen->FillProfile(Form("%s_pt_%i", corrconfigs.at(l_ind).Head.c_str(), i), centmult, val, dnx, rndm) : fFC->FillProfile(Form("%s_pt_%i", corrconfigs.at(l_ind).Head.c_str(), i), centmult, val, dnx, rndm);
      }
    }
    return;
  }

  template <typename TCollision, typename TTracks>
  void processCollision(datatype dt, TCollision const& collision, TTracks tracks, const float& centrality, const int& field)
  {

    if (tracks.size() < 1)
      return;
    if (centrality < centbinning.front() || centrality > centbinning.back())
      return;
    registry.fill(HIST("hEventCount"), 9.5);
    float vtxz = collision.posZ();
    fGFW->Clear();
    fFCpt->clearVector();
    float l_Random = fRndm->Rndm();

    std::vector<int> ITS_TPC_tracks(6, 0);

    for (auto& track : tracks) {
      ProcessTrack(track, centrality, vtxz, field);
      FillITSTPCQA(track, ITS_TPC_tracks);
    }
    registry.fill(HIST("eventQA/nITS_nTPC_Tracks"), ITS_TPC_tracks[1], ITS_TPC_tracks[2]);

    FillOutputContainers(dt, (cfgUseNch) ? tracks.size() : centrality, l_Random);
  }

  template <typename TrackObject>
  inline void FillITSTPCQA(TrackObject const& track, std::vector<int>& ITS_TPC_tracks)
  {
    if constexpr (!framework::has_type_v<aod::mcparticle::Vy, typename TrackObject::all_columns>) {
      // if track its and tpc required then these three are the same
      if (track.hasITS() || track.hasTPC())
        ITS_TPC_tracks[0] += 1;
      if (track.hasITS())
        ITS_TPC_tracks[1] += 1;
      if (track.hasTPC())
        ITS_TPC_tracks[2] += 1;
      if (track.hasITS() && track.hasTPC())
        ITS_TPC_tracks[3] += 1;
      if (track.hasITS() && !track.hasTPC())
        ITS_TPC_tracks[4] += 1;
      if (track.hasTPC() && !track.hasITS())
        ITS_TPC_tracks[5] += 1;
    }
  }

  template <typename TrackObject>
  inline void ProcessTrack(TrackObject const& track, const float& centrality, const float& vtxz, const int& field)
  {
    float weff = 1, wacc = 1;

    auto handleReco = [&](auto const& particle) -> bool {
      if (cfgFillQA)
        FillTrackQA<kBefore>(track, vtxz);

      if (cfgUseAdditionalTrackCut && !trackSelected(track, field))
        return false;

      if (cfgFillWeights)
        fWeights->fill(particle.phi(), particle.eta(), vtxz, particle.pt(), centrality, 0);

      if (!setCurrentParticleWeights(weff, wacc, particle.phi(), particle.eta(), particle.pt(), vtxz))
        return false;

      if (cfgFillQA)
        FillTrackQA<kAfter>(track, vtxz);

      FillGFW(particle, weff, wacc);

      return true;
    };

    if constexpr (framework::has_type_v<aod::mctracklabel::McParticleId, typename TrackObject::all_columns>) {
      if (track.mcParticleId() < 0 || !(track.has_mcParticle()))
        return;

      auto mcParticle = track.mcParticle();
      if (!mcParticle.isPhysicalPrimary() || mcParticle.eta() < etalow || mcParticle.eta() > etaup || mcParticle.pt() < ptlow || mcParticle.pt() > ptup || track.tpcNClsFound() < cfgNcls)
        return;

      if (!handleReco(mcParticle))
        return;

      registry.fill(HIST("trackQA/phi_eta_vtxZ_corrected"), mcParticle.phi(), mcParticle.eta(), vtxz, wacc);

    } else if constexpr (framework::has_type_v<aod::mcparticle::McCollisionId, typename TrackObject::all_columns>) {
      if (!track.isPhysicalPrimary() || track.eta() < etalow || track.eta() > etaup || track.pt() < ptlow || track.pt() > ptup)
        return;

      if (cfgFillQA)
        FillTrackQA<kBefore>(track, vtxz);

      FillGFW(track, 1., 1.);
    } else {
      if (track.tpcNClsFound() < cfgNcls)
        return;

      if (!handleReco(track))
        return;

      registry.fill(HIST("trackQA/phi_eta_vtxZ_corrected"), track.phi(), track.eta(), vtxz, wacc);
    }
  }

  template <typename TrackObject>
  inline void FillGFW(TrackObject track, float weff, float wacc)
  {
    fFCpt->fill(weff, track.pt());
    bool WithinPtPOI = (ptpoilow < track.pt()) && (track.pt() < ptpoiup); // within POI pT range
    bool WithinPtRef = (ptreflow < track.pt()) && (track.pt() < ptrefup); // within RF pT range
    if (WithinPtRef)
      fGFW->Fill(track.eta(), fPtAxis->FindBin(track.pt()) - 1, track.phi(), weff * wacc, 1);
    if (WithinPtPOI)
      fGFW->Fill(track.eta(), fPtAxis->FindBin(track.pt()) - 1, track.phi(), weff * wacc, 2);
    if (WithinPtPOI && WithinPtRef)
      fGFW->Fill(track.eta(), fPtAxis->FindBin(track.pt()) - 1, track.phi(), weff * wacc, 4);
    return;
  }

  template <QAtime qt, typename TrackObject>
  inline void FillTrackQA(TrackObject track, const float vtxz)
  {
    if constexpr (framework::has_type_v<aod::mcparticle::McCollisionId, typename TrackObject::all_columns>) {
      registry.fill(HIST("trackMCGen/phi_eta_vtxZ_gen"), track.phi(), track.eta(), vtxz);
      registry.fill(HIST("trackMCGen/pt_gen"), track.pt());
    } else {
      registry.fill(HIST("trackQA/phi_eta_vtxZ"), track.phi(), track.eta(), vtxz);
      registry.fill(HIST("trackQA/pt_dcaXY_dcaZ"), track.pt(), track.dcaXY(), track.dcaZ());
      registry.fill(HIST("trackQA/pt_reco"), track.pt());
      registry.fill(HIST("trackQA/") + HIST(moment[qt]) + HIST("CrossedTPCRows_pt"), track.pt(), track.tpcNClsCrossedRows());
      registry.fill(HIST("trackQA/") + HIST(moment[qt]) + HIST("NumberSharedClustersTPC_pt"), track.pt(), track.tpcNClsShared());
      registry.fill(HIST("trackQA/") + HIST(moment[qt]) + HIST("ITSClusters_pt"), track.pt(), track.itsNCls());
      registry.fill(HIST("trackQA/") + HIST(moment[qt]) + HIST("DCAxy_pt"), track.pt(), track.dcaXY());
      registry.fill(HIST("trackQA/") + HIST(moment[qt]) + HIST("DCAz_pt"), track.pt(), track.dcaZ());
      registry.fill(HIST("trackQA/") + HIST(moment[qt]) + HIST("tpcSignal_pt"), track.pt(), track.tpcSignal());
      registry.fill(HIST("trackQA/") + HIST(moment[qt]) + HIST("phi_pt"), track.pt(), track.phi());
      registry.fill(HIST("trackQA/") + HIST(moment[qt]) + HIST("phi_eta"), track.eta(), track.phi());
      registry.fill(HIST("trackQA/") + HIST(moment[qt]) + HIST("tpcCrossedRowsOverFindableCls_pt"), track.pt(), track.tpcCrossedRowsOverFindableCls());
    }
  }

  template <QAtime qt, typename CollisionObject, typename TracksObject>
  inline void FillEventQA(CollisionObject collision, TracksObject tracks)
  {

    if constexpr (framework::has_type_v<aod::cent::CentFT0C, typename CollisionObject::all_columns>) {
      registry.fill(HIST("eventQA/") + HIST(moment[qt]) + HIST("globalTracks_centT0C"), collision.centFT0C(), tracks.size());
      registry.fill(HIST("eventQA/") + HIST(moment[qt]) + HIST("PVTracks_centT0C"), collision.centFT0C(), collision.multNTracksPV());
      registry.fill(HIST("eventQA/") + HIST(moment[qt]) + HIST("globalTracks_PVTracks"), collision.multNTracksPV(), tracks.size());
      registry.fill(HIST("eventQA/") + HIST(moment[qt]) + HIST("globalTracks_multT0A"), collision.multFT0A(), tracks.size());
      registry.fill(HIST("eventQA/") + HIST(moment[qt]) + HIST("globalTracks_multV0A"), collision.multFV0A(), tracks.size());
      registry.fill(HIST("eventQA/") + HIST(moment[qt]) + HIST("multV0A_multT0A"), collision.multFT0A(), collision.multFV0A());
      registry.fill(HIST("eventQA/") + HIST(moment[qt]) + HIST("multT0C_centT0C"), collision.centFT0C(), collision.multFT0C());
    }
  }

  template <typename CollisionObject, typename TracksObject>
  inline void RunProcess(CollisionObject collision, TracksObject tracks)
  {

    float centrality;

    if constexpr (framework::has_type_v<aod::cent::CentRun2V0M, typename CollisionObject::all_columns>) {
      if (!collision.sel7())
        return;
      centrality = collision.centRun2V0M();
    } else if constexpr (framework::has_type_v<aod::cent::CentFT0C, typename CollisionObject::all_columns>) {
      registry.fill(HIST("hEventCount"), .5);
      if (!collision.sel8())
        return;
      registry.fill(HIST("hEventCount"), 1.5);
      centrality = collision.centFT0C();

      if (cfgFillQA)
        FillEventQA<kBefore>(collision, tracks);
    }

    if (cfgDoOccupancySel) {
      int occupancy = collision.trackOccupancyInTimeRange();
      if (occupancy < 0 || occupancy > cfgMaxOccupancy)
        return;
      registry.fill(HIST("hEventCount"), 2.5);
    }

    auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
    if (cfgUseAdditionalEventCut && !eventSelected(collision, tracks.size(), centrality))
      return;
    if (cfgFillQA)
      FillEventQA<kAfter>(collision, tracks);

    loadCorrections(bc.timestamp());
    auto field = (cfgMagField == 99999) ? getMagneticField(bc.timestamp()) : cfgMagField;
    processCollision(kReco, collision, tracks, centrality, field);
  }

  Filter collisionFilter = nabs(aod::collision::posZ) < cfgVtxZ;
  Filter trackFilter = nabs(aod::track::eta) < cfgEta && aod::track::pt > cfgPtmin&& aod::track::pt < cfgPtmax && ((requireGlobalTrackInFilter()) || (aod::track::isGlobalTrackSDD == (uint8_t) true)) && nabs(aod::track::dcaXY) < cfgDCAxy&& nabs(aod::track::dcaZ) < cfgDCAz;

  using myTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA>>;

  void processData(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Cs>>::iterator const& collision, aod::BCsWithTimestamps const&, myTracks const& tracks)
  {
    RunProcess(collision, tracks);
  }
  PROCESS_SWITCH(flowAnalysisGF, processData, "Process analysis for non-derived data", true);

  void processMCReco(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Cs>>::iterator const& collision, aod::BCsWithTimestamps const&, soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA, aod::McTrackLabels>> const& tracks, aod::McParticles const&)
  {
    RunProcess(collision, tracks);
  }
  PROCESS_SWITCH(flowAnalysisGF, processMCReco, "Process analysis for MC reconstructed events", false);

  Filter mcCollFilter = nabs(aod::mccollision::posZ) < cfgVtxZ;
  void processMCGen(soa::Filtered<aod::McCollisions>::iterator const& mcCollision, soa::SmallGroups<soa::Join<aod::McCollisionLabels, aod::Collisions, aod::CentFT0Cs>> const& collisions, aod::McParticles const& particles)
  {
    if (collisions.size() != 1)
      return;
    float centrality = -1;
    for (auto& collision : collisions) {
      centrality = collision.centFT0C();
    }
    processCollision(kGen, mcCollision, particles, centrality, -999);
  }
  PROCESS_SWITCH(flowAnalysisGF, processMCGen, "Process analysis for MC generated events", false);

  void processRun2(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentRun2V0Ms>>::iterator const& collision, aod::BCsWithTimestamps const&, myTracks const& tracks)
  {
    RunProcess(collision, tracks);
  }
  PROCESS_SWITCH(flowAnalysisGF, processRun2, "Process analysis for Run 2 converted data", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<flowAnalysisGF>(cfgc),
  };
}
