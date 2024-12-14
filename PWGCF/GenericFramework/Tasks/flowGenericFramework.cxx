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
#include <complex>
#include <string>

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
#include "Common/DataModel/PIDResponse.h"

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

namespace o2::analysis::genericframework
{
std::vector<double> ptbinning = {0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.5, 4, 5, 6, 8, 10};
float ptpoilow = 0.2, ptpoiup = 10.0;
float ptreflow = 0.2, ptrefup = 3.0;
float ptlow = 0.2, ptup = 10.0;
int etabins = 16;
float etalow = -0.8, etaup = 0.8;
int vtxZbins = 40;
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
} // namespace o2::analysis::genericframework

using namespace o2::analysis::genericframework;

struct GenericFramework {

  O2_DEFINE_CONFIGURABLE(cfgNbootstrap, int, 10, "Number of subsamples")
  O2_DEFINE_CONFIGURABLE(cfgMpar, int, 8, "Highest order of pt-pt correlations")
  O2_DEFINE_CONFIGURABLE(cfgUseNch, bool, false, "Do correlations as function of Nch")
  O2_DEFINE_CONFIGURABLE(cfgFillWeights, bool, false, "Fill NUA weights")
  O2_DEFINE_CONFIGURABLE(cfgRunByRunWeights, bool, false, "Use run by run NUA corrections")
  O2_DEFINE_CONFIGURABLE(cfgFillQA, bool, false, "Fill QA histograms")
  O2_DEFINE_CONFIGURABLE(cfgUseAdditionalEventCut, bool, false, "Use additional event cut on mult correlations")
  O2_DEFINE_CONFIGURABLE(cfgUseAdditionalTrackCut, bool, false, "Use additional track cut on phi")
  O2_DEFINE_CONFIGURABLE(cfgUseCentralMoments, bool, true, "Use central moments in vn-pt calculations")
  O2_DEFINE_CONFIGURABLE(cfgUsePID, bool, true, "Enable PID information")
  O2_DEFINE_CONFIGURABLE(cfgUseGapMethod, bool, false, "Use gap method in vn-pt calculations")
  O2_DEFINE_CONFIGURABLE(cfgEfficiency, std::string, "", "CCDB path to efficiency object")
  O2_DEFINE_CONFIGURABLE(cfgAcceptance, std::string, "", "CCDB path to acceptance object")
  O2_DEFINE_CONFIGURABLE(cfgDCAxy, float, 0.2, "Cut on DCA in the transverse direction (cm)");
  O2_DEFINE_CONFIGURABLE(cfgDCAz, float, 2, "Cut on DCA in the longitudinal direction (cm)");
  O2_DEFINE_CONFIGURABLE(cfgNcls, float, 70, "Cut on number of TPC clusters found");
  O2_DEFINE_CONFIGURABLE(cfgPtmin, float, 0.2, "minimum pt (GeV/c)");
  O2_DEFINE_CONFIGURABLE(cfgPtmax, float, 10, "maximum pt (GeV/c)");
  O2_DEFINE_CONFIGURABLE(cfgEta, float, 0.8, "eta cut");
  O2_DEFINE_CONFIGURABLE(cfgEtaPtPt, float, 0.4, "eta cut for pt-pt correlations");
  O2_DEFINE_CONFIGURABLE(cfgVtxZ, float, 10, "vertex cut (cm)");
  O2_DEFINE_CONFIGURABLE(cfgOccupancySelection, int, -999, "Max occupancy selection, -999 to disable");
  O2_DEFINE_CONFIGURABLE(cfgNoSameBunchPileupCut, bool, true, "kNoSameBunchPileupCut");
  O2_DEFINE_CONFIGURABLE(cfgIsGoodZvtxFT0vsPV, bool, true, "kIsGoodZvtxFT0vsPV");
  O2_DEFINE_CONFIGURABLE(cfgNoCollInTimeRangeStandard, bool, true, "kNoCollInTimeRangeStandard");
  O2_DEFINE_CONFIGURABLE(cfgDoOccupancySel, bool, true, "Bool for event selection on detector occupancy");
  O2_DEFINE_CONFIGURABLE(cfgMultCut, bool, true, "Use additional evenr cut on mult correlations");
  O2_DEFINE_CONFIGURABLE(cfgTVXinTRD, bool, true, "Use kTVXinTRD (reject TRD triggered events)");
  O2_DEFINE_CONFIGURABLE(cfgIsVertexITSTPC, bool, true, "Selects collisions with at least one ITS-TPC track");
  O2_DEFINE_CONFIGURABLE(cfgMagField, float, 99999, "Configurable magnetic field; default CCDB will be queried");

  Configurable<GFWBinningCuts> cfgGFWBinning{"cfgGFWBinning", {40, 16, 72, 300, 0, 3000, 0.2, 10.0, 0.2, 3.0, {0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.5, 4, 5, 6, 8, 10}, {0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90}}, "Configuration for binning"};
  Configurable<GFWRegions> cfgRegions{"cfgRegions", {{"refN", "refP", "refFull"}, {-0.8, 0.4, -0.8}, {-0.4, 0.8, 0.8}, {0, 0, 0}, {1, 1, 1}}, "Configurations for GFW regions"};

  Configurable<GFWCorrConfigs> cfgCorrConfig{"cfgCorrConfig", {{"refP {2} refN {-2}", "refP {3} refN {-3}", "refP {4} refN {-4}", "refFull {2 -2}", "refFull {2 2 -2 -2}"}, {"ChGap22", "ChGap32", "ChGap42", "ChFull22", "ChFull24"}, {0, 0, 0, 0, 0}, {15, 1, 1, 0, 0}}, "Configurations for each correlation to calculate"};

  // #include "PWGCF/TwoParticleCorrelations/TableProducer/Productions/skimmingconf_20221115.cxx" // NOLINT
  //  Connect to ccdb
  Service<ccdb::BasicCCDBManager> ccdb;

  struct Config {
    TH1D* mEfficiency = nullptr;
    std::vector<GFWWeights*> mAcceptance;
    bool correctionsLoaded = false;
  } cfg;

  // Define output
  OutputObj<FlowContainer> fFC{FlowContainer("FlowContainer")};
  OutputObj<FlowPtContainer> fFCpt{FlowPtContainer("FlowPtContainer")};
  OutputObj<FlowContainer> fFC_gen{FlowContainer("FlowContainer_gen")};
  OutputObj<TList> fWeightList{"WeightList", OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry registry{"registry"};

  // define global variables
  GFW* fGFW = new GFW();
  std::vector<GFW::CorrConfig> corrconfigs;
  TRandom3* fRndm = new TRandom3(0);
  TAxis* fPtAxis;
  int lastRun = 0;

  // Event selection cuts - Alex
  TF1* fPhiCutLow = nullptr;
  TF1* fPhiCutHigh = nullptr;
  TF1* fMultPVCutLow = nullptr;
  TF1* fMultPVCutHigh = nullptr;
  TF1* fMultCutLow = nullptr;
  TF1* fMultCutHigh = nullptr;
  TF1* fMultMultPVCut = nullptr;

  void init(InitContext const&)
  {
    LOGF(info, "flowGenericFramework::init()");
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
    AxisSpec dcaZAXis = {200, -2, 2, "DCA_{z} (cm)"};
    AxisSpec dcaXYAXis = {200, -1, 1, "DCA_{xy} (cm)"};
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    int64_t now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    ccdb->setCreatedNotAfter(now);

    int ptbins = ptbinning.size() - 1;
    fPtAxis = new TAxis(ptbins, &ptbinning[0]);

    TList* weightlist = new TList();
    weightlist->SetOwner(true);
    fWeightList.setObject(weightlist);

    if (!cfgRunByRunWeights && cfgFillWeights) {
      if (cfgUsePID) {
        std::vector<GFWWeights*> weights;
        std::vector<std::string> species = {"ref", "ch", "pi", "ka", "pr"};
        for (size_t i = 0; i < species.size(); ++i) {
          weights.push_back(new GFWWeights(Form("w_%s", species[i].c_str())));
          if (i == 0) {
            auto it = std::find(ptbinning.begin(), ptbinning.end(), ptrefup);
            std::vector<double> refpt(ptbinning.begin(), it + 1);
            weights[i]->SetPtBins(refpt.size() - 1, &refpt[0]);
          } else {
            weights[i]->SetPtBins(fPtAxis->GetNbins(), &ptbinning[0]);
          }
          weights[i]->Init(true, false);
          fWeightList->Add(weights[i]);
        }
      } else {
        GFWWeights* weight = new GFWWeights("w_ch");
        weight->SetPtBins(fPtAxis->GetNbins(), &ptbinning[0]);
        weight->Init(true, false);
        fWeightList->Add(weight);
      }
    }

    if (doprocessMCGen) {
      registry.add("MCGen/before/pt_gen", "", {HistType::kTH1D, {ptAxis}});
      registry.add("MCGen/before/phi_eta_vtxZ_gen", "", {HistType::kTH3D, {phiAxis, etaAxis, vtxAxis}});
      registry.addClone("MCGen/before/", "MCGen/after/");
    }
    if (doprocessMCReco || doprocessData || doprocessRun2) {
      registry.add("trackQA/before/phi_eta_vtxZ", "", {HistType::kTH3D, {phiAxis, etaAxis, vtxAxis}});
      registry.add("trackQA/before/pt_dcaXY_dcaZ", "", {HistType::kTH3D, {ptAxis, dcaXYAXis, dcaZAXis}});
      registry.add("trackQA/before/pt_phi", "", {HistType::kTH2D, {ptAxis, phiModAxis}});
      registry.addClone("trackQA/before/", "trackQA/after/");
      registry.add("trackQA/after/pt_ref", "", {HistType::kTH1D, {{100, ptreflow, ptrefup}}});
      registry.add("trackQA/after/pt_poi", "", {HistType::kTH1D, {{100, ptpoilow, ptpoiup}}});

      registry.add("eventQA/before/globalTracks_centT0C", "", {HistType::kTH2D, {centAxis, nchAxis}});
      registry.add("eventQA/before/PVTracks_centT0C", "", {HistType::kTH2D, {centAxis, multpvAxis}});
      registry.add("eventQA/before/globalTracks_PVTracks", "", {HistType::kTH2D, {multpvAxis, nchAxis}});
      registry.add("eventQA/before/globalTracks_multT0A", "", {HistType::kTH2D, {t0aAxis, nchAxis}});
      registry.add("eventQA/before/globalTracks_multV0A", "", {HistType::kTH2D, {t0aAxis, nchAxis}});
      registry.add("eventQA/before/multV0A_multT0A", "", {HistType::kTH2D, {t0aAxis, t0aAxis}});
      registry.add("eventQA/before/multT0C_centT0C", "", {HistType::kTH2D, {centAxis, t0cAxis}});
      registry.addClone("eventQA/before/", "eventQA/after/");
      registry.add("eventQA/eventSel", "Number of Events;; Counts", {HistType::kTH1D, {{10, 0, 10}}});
      registry.get<TH1>(HIST("eventQA/eventSel"))->GetXaxis()->SetBinLabel(1, "Filtered event");
      registry.get<TH1>(HIST("eventQA/eventSel"))->GetXaxis()->SetBinLabel(2, "sel8");
      registry.get<TH1>(HIST("eventQA/eventSel"))->GetXaxis()->SetBinLabel(3, "occupancy");
      registry.get<TH1>(HIST("eventQA/eventSel"))->GetXaxis()->SetBinLabel(4, "kTVXinTRD");
      registry.get<TH1>(HIST("eventQA/eventSel"))->GetXaxis()->SetBinLabel(5, "kNoSameBunchPileup");
      registry.get<TH1>(HIST("eventQA/eventSel"))->GetXaxis()->SetBinLabel(6, "kIsGoodZvtxFT0vsPV");
      registry.get<TH1>(HIST("eventQA/eventSel"))->GetXaxis()->SetBinLabel(7, "kNoCollInTimeRangeStandard");
      registry.get<TH1>(HIST("eventQA/eventSel"))->GetXaxis()->SetBinLabel(8, "kIsVertexITSTPC");
      registry.get<TH1>(HIST("eventQA/eventSel"))->GetXaxis()->SetBinLabel(9, "after Mult cuts");
      registry.get<TH1>(HIST("eventQA/eventSel"))->GetXaxis()->SetBinLabel(10, "has track + within cent");
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
    if (doprocessData || doprocessRun2 || doprocessMCReco) {
      fFC->SetName("FlowContainer");
      fFC->SetXAxis(fPtAxis);
      fFC->Initialize(oba, multAxis, cfgNbootstrap);
    }
    if (doprocessMCGen) {
      fFC_gen->SetName("FlowContainer_gen");
      fFC_gen->SetXAxis(fPtAxis);
      fFC_gen->Initialize(oba, multAxis, cfgNbootstrap);
    }
    delete oba;
    fFCpt->SetUseCentralMoments(cfgUseCentralMoments);
    fFCpt->SetUseGapMethod(cfgUseGapMethod);
    fFCpt->Initialise(multAxis, cfgMpar, configs, cfgNbootstrap);
    // Event selection - Alex
    if (cfgUseAdditionalEventCut) {
      /*
      //22s cuts
      fMultPVCutLow = new TF1("fMultPVCutLow", "[0]+[1]*x+[2]*x*x+[3]*x*x*x - 2.5*([4]+[5]*x+[6]*x*x+[7]*x*x*x+[8]*x*x*x*x)", 0, 100);
      fMultPVCutLow->SetParameters(2834.66, -87.0127, 0.915126, -0.00330136, 332.513, -12.3476, 0.251663, -0.00272819, 1.12242e-05);
      fMultPVCutHigh = new TF1("fMultPVCutHigh", "[0]+[1]*x+[2]*x*x+[3]*x*x*x + 2.5*([4]+[5]*x+[6]*x*x+[7]*x*x*x+[8]*x*x*x*x)", 0, 100);
      fMultPVCutHigh->SetParameters(2834.66, -87.0127, 0.915126, -0.00330136, 332.513, -12.3476, 0.251663, -0.00272819, 1.12242e-05);

      fMultCutLow = new TF1("fMultCutLow", "[0]+[1]*x+[2]*x*x+[3]*x*x*x - 2.5*([4]+[5]*x)", 0, 100);
      fMultCutLow->SetParameters(1893.94, -53.86, 0.502913, -0.0015122, 109.625, -1.19253);
      fMultCutHigh = new TF1("fMultCutHigh", "[0]+[1]*x+[2]*x*x+[3]*x*x*x + 3.*([4]+[5]*x)", 0, 100);
      fMultCutHigh->SetParameters(1893.94, -53.86, 0.502913, -0.0015122, 109.625, -1.19253);
      fMultMultPVCut = new TF1("fMultMultPVCut", "[0]+[1]*x+[2]*x*x", 0, 5000);
      fMultMultPVCut->SetParameters(-0.1, 0.785, -4.7e-05);
      */
      fMultPVCutLow = new TF1("fMultPVCutLow", "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x - 3.5*([5]+[6]*x+[7]*x*x+[8]*x*x*x+[9]*x*x*x*x)", 0, 100);
      fMultPVCutLow->SetParameters(3257.29, -121.848, 1.98492, -0.0172128, 6.47528e-05, 154.756, -1.86072, -0.0274713, 0.000633499, -3.37757e-06);
      fMultPVCutHigh = new TF1("fMultPVCutHigh", "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x + 3.5*([5]+[6]*x+[7]*x*x+[8]*x*x*x+[9]*x*x*x*x)", 0, 100);
      fMultPVCutHigh->SetParameters(3257.29, -121.848, 1.98492, -0.0172128, 6.47528e-05, 154.756, -1.86072, -0.0274713, 0.000633499, -3.37757e-06);

      fMultCutLow = new TF1("fMultCutLow", "[0]+[1]*x+[2]*x*x+[3]*x*x*x - 2.*([4]+[5]*x+[6]*x*x+[7]*x*x*x+[8]*x*x*x*x)", 0, 100);
      fMultCutLow->SetParameters(1654.46, -47.2379, 0.449833, -0.0014125, 150.773, -3.67334, 0.0530503, -0.000614061, 3.15956e-06);
      fMultCutHigh = new TF1("fMultCutHigh", "[0]+[1]*x+[2]*x*x+[3]*x*x*x + 3.*([4]+[5]*x+[6]*x*x+[7]*x*x*x+[8]*x*x*x*x)", 0, 100);
      fMultCutHigh->SetParameters(1654.46, -47.2379, 0.449833, -0.0014125, 150.773, -3.67334, 0.0530503, -0.000614061, 3.15956e-06);
    }

    if (cfgUseAdditionalTrackCut) {
      fPhiCutLow = new TF1("fPhiCutLow", "0.06/x+pi/18.0-0.06", 0, 100);
      fPhiCutHigh = new TF1("fPhiCutHigh", "0.1/x+pi/18.0+0.06", 0, 100);
    }
  }

  static constexpr std::string_view moment[] = {"before/", "after/"};

  enum QAFillTime {
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
      // grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>("GLO/GRP/GRP", timestamp);
      grpo = ccdb->getForTimeStamp<o2::parameters::GRPMagField>("GLO/Config/GRPMagField", timestamp);
      if (grpo == nullptr) {
        LOGF(fatal, "GRP object not found for timestamp %llu", timestamp);
        return 0;
      }
      LOGF(info, "Retrieved GRP for timestamp %llu with magnetic field of %d kG", timestamp, grpo->getNominalL3Field());
    }
    return grpo->getNominalL3Field();
  }

  void loadCorrections(aod::BCsWithTimestamps::iterator const& bc)
  {
    uint64_t timestamp = bc.timestamp();
    int run = bc.runNumber();
    if (cfg.correctionsLoaded) {
      if (!cfgRunByRunWeights)
        return;
      if (run == lastRun)
        return;
    }
    if (cfgUsePID) {
      if (cfgAcceptance.value.empty() == false) {
        if (cfgRunByRunWeights) { // run-by-run NUA weights from ccdb, stored in TList to hold PID weights
          TList* weightlist = ccdb->getForTimeStamp<TList>(cfgAcceptance, timestamp);
          cfg.mAcceptance.push_back(dynamic_cast<GFWWeights*>(weightlist->FindObject(Form("w%i_ref", run))));
          cfg.mAcceptance.push_back(dynamic_cast<GFWWeights*>(weightlist->FindObject(Form("w%i_ch", run))));
          cfg.mAcceptance.push_back(dynamic_cast<GFWWeights*>(weightlist->FindObject(Form("w%i_pi", run))));
          cfg.mAcceptance.push_back(dynamic_cast<GFWWeights*>(weightlist->FindObject(Form("w%i_ka", run))));
          cfg.mAcceptance.push_back(dynamic_cast<GFWWeights*>(weightlist->FindObject(Form("w%i_pr", run))));
        } else { // run-averaged weights, stored in TList to hold PID weights
          TList* weightlist = ccdb->getForTimeStamp<TList>(cfgAcceptance, timestamp);
          weightlist->ls();
          cfg.mAcceptance.push_back(dynamic_cast<GFWWeights*>(weightlist->FindObject("weights_ref")));
          cfg.mAcceptance.push_back(dynamic_cast<GFWWeights*>(weightlist->FindObject("weights_ch")));
          cfg.mAcceptance.push_back(dynamic_cast<GFWWeights*>(weightlist->FindObject("weights_pi")));
          cfg.mAcceptance.push_back(dynamic_cast<GFWWeights*>(weightlist->FindObject("weights_ka")));
          cfg.mAcceptance.push_back(dynamic_cast<GFWWeights*>(weightlist->FindObject("weights_pr")));
        }
        if (!cfg.mAcceptance.empty())
          LOGF(info, "Loaded acceptance weights from %s", cfgAcceptance.value.c_str());
        else
          LOGF(warning, "Could not load acceptance weights from %s", cfgAcceptance.value.c_str());
      }
    } else {
      if (cfgAcceptance.value.empty() == false) {
        if (cfgRunByRunWeights) { // run-by-run NUA weights from ccdb, stored in TList
          TList* weightlist = ccdb->getForTimeStamp<TList>(cfgAcceptance, timestamp);
          cfg.mAcceptance.push_back(dynamic_cast<GFWWeights*>(weightlist->FindObject(Form("w%i_ch", run))));
        } else { // run-averaged weights, stored in TList
          TList* weightlist = ccdb->getForTimeStamp<TList>(cfgAcceptance, timestamp);
          cfg.mAcceptance.push_back(dynamic_cast<GFWWeights*>(weightlist->FindObject("w_ch")));
        }
        if (!cfg.mAcceptance.empty())
          LOGF(info, "Loaded acceptance weights from %s", cfgAcceptance.value.c_str());
        else
          LOGF(warning, "Could not load acceptance weights from %s", cfgAcceptance.value.c_str());
      }
    }
    if (cfgEfficiency.value.empty() == false) {
      cfg.mEfficiency = ccdb->getForTimeStamp<TH1D>(cfgEfficiency, timestamp);
      if (cfg.mEfficiency == nullptr) {
        LOGF(fatal, "Could not load efficiency histogram from %s", cfgEfficiency.value.c_str());
      }
      LOGF(info, "Loaded efficiency histogram from %s (%p)", cfgEfficiency.value.c_str(), (void*)cfg.mEfficiency);
    }
    cfg.correctionsLoaded = true;
  }

  template <typename TTrack>
  double getAcceptance(TTrack track, const double& vtxz, int index)
  { //-1 ref, 0 ch, 1 pi, 2 ka, 3 pr
    double wacc = 1;
    index += 1;
    if (!cfg.mAcceptance.empty()) {
      if (cfgUsePID) {
        wacc = cfg.mAcceptance[index]->GetNUA(track.phi(), track.eta(), vtxz);
      } else {
        wacc = cfg.mAcceptance[0]->GetNUA(track.phi(), track.eta(), vtxz);
      }
    }
    return wacc;
  }

  template <typename TTrack>
  double getEfficiency(TTrack track)
  { //-1 ref, 0 ch, 1 pi, 2 ka, 3 pr
    double eff = 1.;
    if (cfg.mEfficiency)
      eff = cfg.mEfficiency->GetBinContent(cfg.mEfficiency->FindBin(track.pt()));
    if (eff == 0)
      return -1.;
    else
      return 1. / eff;
  }
  // Obsolete for now untill service wagons get added
  /*   template<typename TTrack>
    int getBayesPIDIndex(TTrack track) {
      float maxProb[3] = {0.95,0.85,0.85};
      int pidID = 0;
      if(track.bayesID()==o2::track::PID::Pion || track.bayesID()==o2::track::PID::Kaon || track.bayesID()==o2::track::PID::Proton){
        pidID = track.bayesID()-1;  //Realign
        float nsigmaTPC[3] = {track.tpcNSigmaPi(),track.tpcNSigmaKa(),track.tpcNSigmaPr()};
        float nsigmaTOF[3] = {track.tofNSigmaPi(),track.tofNSigmaKa(),track.tofNSigmaPr()};
        if(track.bayesProb() > maxProb[pidID-1]) {
          if(abs(nsigmaTPC[pidID-1]) > 3) return 0;
          if(abs(nsigmaTOF[pidID-1]) > 3) return 0;
          return pidID;
        }
        else return 0;
      }
      return 0;
    } */

  template <typename TTrack>
  int GetNsigmaPID(TTrack track)
  {
    // Computing Nsigma arrays for pion, kaon, and protons
    std::array<float, 3> nSigmaTPC = {track.tpcNSigmaPi(), track.tpcNSigmaKa(), track.tpcNSigmaPr()};
    std::array<float, 3> nSigmaCombined = {std::hypot(track.tpcNSigmaPi(), track.tofNSigmaPi()), std::hypot(track.tpcNSigmaKa(), track.tofNSigmaKa()), std::hypot(track.tpcNSigmaPr(), track.tofNSigmaPr())};
    int pid = -1;
    float nsigma = 3.0;

    // Choose which nSigma to use
    std::array<float, 3> nSigmaToUse = (track.pt() > 0.4 && track.hasTOF()) ? nSigmaCombined : nSigmaTPC;

    // Select particle with the lowest nsigma
    for (int i = 0; i < 3; ++i) {
      if (std::abs(nSigmaToUse[i]) < nsigma) {
        pid = i;
        nsigma = std::abs(nSigmaToUse[i]);
      }
    }
    return pid + 1; // shift the pid by 1, 1 = pion, 2 = kaon, 3 = proton
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
      registry.fill(HIST("eventQA/eventSel"), 3.5);
    }

    if (cfgNoSameBunchPileupCut) {
      if (!collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
        // rejects collisions which are associated with the same "found-by-T0" bunch crossing
        // https://indico.cern.ch/event/1396220/#1-event-selection-with-its-rof
        return 0;
      }
      registry.fill(HIST("eventQA/eventSel"), 4.5);
    }
    if (cfgIsGoodZvtxFT0vsPV) {
      if (!collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
        // removes collisions with large differences between z of PV by tracks and z of PV from FT0 A-C time difference
        // use this cut at low multiplicities with caution
        return 0;
      }
      registry.fill(HIST("eventQA/eventSel"), 5.5);
    }
    if (cfgNoCollInTimeRangeStandard) {
      if (!collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
        //  Rejection of the collisions which have other events nearby
        return 0;
      }
      registry.fill(HIST("eventQA/eventSel"), 6.5);
    }

    if (cfgIsVertexITSTPC) {
      if (!collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) {
        // selects collisions with at least one ITS-TPC track, and thus rejects vertices built from ITS-only tracks
        return 0;
      }
      registry.fill(HIST("eventQA/eventSel"), 7.5);
    }
    float vtxz = -999;
    if (collision.numContrib() > 1) {
      vtxz = collision.posZ();
      float zRes = TMath::Sqrt(collision.covZZ());
      if (zRes > 0.25 && collision.numContrib() < 20)
        vtxz = -999;
    }
    // auto multV0A = collision.multFV0A();
    // auto multT0A = collision.multFT0A();
    // auto multT0C = collision.multFT0C();
    auto multNTracksPV = collision.multNTracksPV();

    if (vtxz > vtxZup || vtxz < vtxZlow)
      return 0;
    if (multNTracksPV < fMultPVCutLow->Eval(centrality))
      return 0;
    if (multNTracksPV > fMultPVCutHigh->Eval(centrality))
      return 0;
    if (multTrk < fMultCutLow->Eval(centrality))
      return 0;
    if (multTrk > fMultCutHigh->Eval(centrality))
      return 0;
    registry.fill(HIST("eventQA/eventSel"), 8.5);
    /* 22s
    if (multNTracksPV < fMultPVCutLow->Eval(centrality))
      return 0;
    if (multNTracksPV > fMultPVCutHigh->Eval(centrality))
      return 0;
    if (multTrk < fMultCutLow->Eval(centrality))
      return 0;
    if (multTrk > fMultCutHigh->Eval(centrality))
      return 0;
    if (multTrk > fMultMultPVCut->Eval(multNTracksPV))
      return 0;
    */
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
      LOGF(warning, "phi < 0: %g", phimodn);

    phimodn += TMath::Pi() / 18.0; // to center gap in the middle
    phimodn = fmod(phimodn, TMath::Pi() / 9.0);
    if (cfgFillQA)
      registry.fill(HIST("trackQA/before/pt_phi"), track.pt(), phimodn);
    if (phimodn < fPhiCutHigh->Eval(track.pt()) && phimodn > fPhiCutLow->Eval(track.pt()))
      return false; // reject track
    if (cfgFillQA)
      registry.fill(HIST("trackQA/after/pt_phi"), track.pt(), phimodn);
    return true;
  }

  enum datatype {
    kReco,
    kGen
  };

  template <typename TTrack>
  void FillWeights(const TTrack track, const double vtxz, const double multcent, int pid_index)
  {
    if (cfgUsePID) {
      std::vector<std::string> species = {"ref", "ch", "pi", "ka", "pr"};
      double ptpidmins[] = {ptpoilow, ptpoilow, 0.3, 0.5};                                           // min pt for ch, pi, ka, pr
      double ptpidmaxs[] = {ptpoiup, ptpoiup, 6.0, 6.0};                                             // max pt for ch, pi, ka, pr
      bool withinPtPOI = (ptpidmins[pid_index] < track.pt()) && (track.pt() < ptpidmaxs[pid_index]); // within POI pT range
      bool withinPtRef = (ptreflow < track.pt()) && (track.pt() < ptrefup);                          // within RF pT range
      if (cfgRunByRunWeights) {
        if (withinPtRef && !pid_index)
          dynamic_cast<GFWWeights*>(fWeightList->FindObject(Form("w%i_%s", lastRun, species[pid_index].c_str())))->Fill(track.phi(), track.eta(), vtxz, track.pt(), multcent, 0); // pt-subset of charged particles for ref flow
        if (withinPtPOI)
          dynamic_cast<GFWWeights*>(fWeightList->FindObject(Form("w%i_%s", lastRun, species[pid_index + 1].c_str())))->Fill(track.phi(), track.eta(), vtxz, track.pt(), multcent, 0); // charged and id'ed particle weights
      } else {
        if (withinPtRef && !pid_index)
          dynamic_cast<GFWWeights*>(fWeightList->FindObject(Form("w_%s", species[pid_index].c_str())))->Fill(track.phi(), track.eta(), vtxz, track.pt(), multcent, 0); // pt-subset of charged particles for ref flow
        if (withinPtPOI)
          dynamic_cast<GFWWeights*>(fWeightList->FindObject(Form("w_%s", species[pid_index + 1].c_str())))->Fill(track.phi(), track.eta(), vtxz, track.pt(), multcent, 0); // charged and id'ed particle weights
      }
    } else {
      if (cfgRunByRunWeights)
        dynamic_cast<GFWWeights*>(fWeightList->FindObject(Form("w%i_ch", lastRun)))->Fill(track.phi(), track.eta(), vtxz, track.pt(), multcent, 0);
      else
        dynamic_cast<GFWWeights*>(fWeightList->FindObject("w_ch"))->Fill(track.phi(), track.eta(), vtxz, track.pt(), multcent, 0);
    }
    return;
  }

  void CreateRunByRunWeights()
  {
    if (cfgUsePID) {
      std::vector<GFWWeights*> weights;
      std::vector<std::string> species = {"ref", "ch", "pi", "ka", "pr"};
      for (size_t i = 0; i < species.size(); ++i) {
        weights.push_back(new GFWWeights(Form("w%i_%s", lastRun, species[i].c_str())));
        if (i == 0) {
          auto it = std::find(ptbinning.begin(), ptbinning.end(), ptrefup);
          std::vector<double> refpt(ptbinning.begin(), it + 1);
          weights[i]->SetPtBins(refpt.size() - 1, &refpt[0]);
        } else {
          weights[i]->SetPtBins(fPtAxis->GetNbins(), &ptbinning[0]);
        }
        weights[i]->Init(true, false);
        fWeightList->Add(weights[i]);
      }
    } else {
      GFWWeights* weight = new GFWWeights(Form("w%i_ch", lastRun));
      weight->SetPtBins(fPtAxis->GetNbins(), &ptbinning[0]);
      weight->Init(true, false);
      fWeightList->Add(weight);
    }

    return;
  }

  template <datatype dt>
  void FillOutputContainers(const float& centmult, const double& rndm)
  {
    fFCpt->CalculateCorrelations();
    fFCpt->FillPtProfiles(centmult, rndm);
    fFCpt->FillCMProfiles(centmult, rndm);
    if (!cfgUseGapMethod)
      fFCpt->FillVnPtStdProfiles(centmult, rndm);
    for (uint l_ind = 0; l_ind < corrconfigs.size(); ++l_ind) {
      if (!corrconfigs.at(l_ind).pTDif) {
        auto dnx = fGFW->Calculate(corrconfigs.at(l_ind), 0, kTRUE).real();
        if (dnx == 0)
          continue;
        auto val = fGFW->Calculate(corrconfigs.at(l_ind), 0, kFALSE).real() / dnx;
        if (TMath::Abs(val) < 1) {
          (dt == kGen) ? fFC_gen->FillProfile(corrconfigs.at(l_ind).Head.c_str(), centmult, val, dnx, rndm) : fFC->FillProfile(corrconfigs.at(l_ind).Head.c_str(), centmult, val, dnx, rndm);
          if (cfgUseGapMethod)
            fFCpt->FillVnPtProfiles(centmult, val, dnx, rndm, configs.GetpTCorrMasks()[l_ind]);
        }
        continue;
      }
      for (Int_t i = 1; i <= fPtAxis->GetNbins(); i++) {
        auto dnx = fGFW->Calculate(corrconfigs.at(l_ind), i - 1, kTRUE).real();
        if (dnx == 0)
          continue;
        auto val = fGFW->Calculate(corrconfigs.at(l_ind), i - 1, kFALSE).real() / dnx;
        if (TMath::Abs(val) < 1)
          (dt == kGen) ? fFC_gen->FillProfile(Form("%s_pt_%i", corrconfigs.at(l_ind).Head.c_str(), i), centmult, val, dnx, rndm) : fFC->FillProfile(Form("%s_pt_%i", corrconfigs.at(l_ind).Head.c_str(), i), centmult, val, dnx, rndm);
      }
    }
    return;
  }

  template <datatype dt, typename TCollision, typename TTracks>
  void processCollision(TCollision collision, TTracks tracks, const float& centrality, const int& field)
  {
    if (tracks.size() < 1)
      return;
    if (centrality < centbinning.front() || centrality > centbinning.back())
      return;
    registry.fill(HIST("eventQA/eventSel"), 9.5);
    float vtxz = collision.posZ();
    fGFW->Clear();
    fFCpt->ClearVector();
    fFCpt->ClearArray();
    float l_Random = fRndm->Rndm();
    for (auto& track : tracks) {
      ProcessTrack(track, centrality, vtxz, field);
    }
    if (!cfgFillWeights)
      FillOutputContainers<dt>((cfgUseNch) ? tracks.size() : centrality, l_Random);
  }

  template <typename TTrack>
  inline void ProcessTrack(TTrack const& track, const float& centrality, const float& vtxz, const int& field)
  {
    if constexpr (framework::has_type_v<aod::mctracklabel::McParticleId, typename TTrack::all_columns>) {
      if (track.mcParticleId() < 0 || !(track.has_mcParticle()))
        return;

      auto mcParticle = track.mcParticle();
      if (!mcParticle.isPhysicalPrimary())
        return;
      if (cfgFillQA)
        FillTrackQA<kReco, kBefore>(track, vtxz);

      if (mcParticle.eta() < etalow || mcParticle.eta() > etaup || mcParticle.pt() < ptlow || mcParticle.pt() > ptup || track.tpcNClsFound() < cfgNcls)
        return;

      if (cfgUseAdditionalTrackCut && !trackSelected(track, field))
        return;

      int pid_index = 0;
      if (cfgUsePID) {
        if (mcParticle.pdgCode() == 211)
          pid_index = 1;
        if (mcParticle.pdgCode() == 321)
          pid_index = 2;
        if (mcParticle.pdgCode() == 2212)
          pid_index = 3;
      }

      if (cfgFillWeights) {
        FillWeights(mcParticle, vtxz, centrality, 0);
      } else {
        FillPtSums<kReco>(track, vtxz);
        FillGFW<kReco>(mcParticle, vtxz, pid_index);
      }

      if (cfgFillQA)
        FillTrackQA<kReco, kAfter>(track, vtxz);

    } else if constexpr (framework::has_type_v<aod::mcparticle::McCollisionId, typename TTrack::all_columns>) {
      if (!track.isPhysicalPrimary())
        return;
      if (cfgFillQA)
        FillTrackQA<kGen, kBefore>(track, vtxz);

      if (track.eta() < etalow || track.eta() > etaup || track.pt() < ptlow || track.pt() > ptup)
        return;

      int pid_index = 0;
      if (cfgUsePID) {
        if (track.pdgCode() == 211)
          pid_index = 1;
        if (track.pdgCode() == 321)
          pid_index = 2;
        if (track.pdgCode() == 2212)
          pid_index = 3;
      }

      FillPtSums<kGen>(track, vtxz);
      FillGFW<kGen>(track, vtxz, pid_index);

      if (cfgFillQA)
        FillTrackQA<kGen, kAfter>(track, vtxz);
    } else {
      if (cfgFillQA)
        FillTrackQA<kReco, kBefore>(track, vtxz);
      if (track.tpcNClsFound() < cfgNcls)
        return;

      if (cfgUseAdditionalTrackCut && !trackSelected(track, field))
        return;

      int pid_index = 0;
      if (cfgUsePID) {
        // pid_index = getBayesPIDIndex(track);
        pid_index = GetNsigmaPID(track);
      }
      if (cfgFillWeights) {
        FillWeights(track, vtxz, centrality, pid_index);
      } else {
        FillPtSums<kReco>(track, vtxz);
        FillGFW<kReco>(track, vtxz, pid_index);
      }
      if (cfgFillQA)
        FillTrackQA<kReco, kAfter>(track, vtxz);
    }
  }

  template <datatype dt, typename TTrack>
  inline void FillPtSums(TTrack track, const double& vtxz)
  {
    double wacc = (dt == kGen) ? 1. : getAcceptance(track, vtxz, -1);
    double weff = (dt == kGen) ? 1. : getEfficiency(track);
    if (weff < 0)
      return;
    if (fabs(track.eta()) < cfgEtaPtPt) {
      fFCpt->Fill(weff, track.pt());
    }
    if (!cfgUseGapMethod) {
      std::complex<double> q2p = {weff * wacc * cos(2 * track.phi()), weff * wacc * sin(2 * track.phi())};
      std::complex<double> q2n = {weff * wacc * cos(-2 * track.phi()), weff * wacc * sin(-2 * track.phi())};
      fFCpt->FillArray(q2p, q2n, weff * track.pt(), weff);
      fFCpt->FillArray(weff * wacc, weff * wacc, weff, weff);
    }
  }

  template <datatype dt, typename TTrack>
  inline void FillGFW(TTrack track, const double& vtxz, int pid_index)
  {
    if (cfgUsePID) { // Analysing POI flow with id'ed particles
      double ptmins[] = {ptpoilow, ptpoilow, 0.3, 0.5};
      double ptmaxs[] = {ptpoiup, ptpoiup, 6.0, 6.0};
      bool WithinPtRef = (track.pt() > ptreflow && track.pt() < ptrefup);
      bool WithinPtPOI = (track.pt() > ptmins[pid_index] && track.pt() < ptmaxs[pid_index]);
      bool WithinPtNch = (track.pt() > ptmins[0] && track.pt() < ptmaxs[0]);
      if (!WithinPtPOI && !WithinPtRef)
        return;
      double wacc_ref = (dt == kGen) ? 1. : getAcceptance(track, vtxz, -1);
      double wacc_poi = (dt == kGen) ? 1. : WithinPtPOI ? getAcceptance(track, vtxz, pid_index)
                                                        : getAcceptance(track, vtxz, 0); //
      if (WithinPtRef && WithinPtPOI && pid_index)
        wacc_ref = wacc_poi; // if particle is both (then it's overlap), override ref with POI
      if (WithinPtRef)
        fGFW->Fill(track.eta(), fPtAxis->FindBin(track.pt()) - 1, track.phi(), wacc_ref, 1);
      if (WithinPtPOI && pid_index)
        fGFW->Fill(track.eta(), fPtAxis->FindBin(track.pt()) - 1, track.phi(), wacc_poi, (1 << (pid_index + 1)));
      if (WithinPtNch)
        fGFW->Fill(track.eta(), fPtAxis->FindBin(track.pt()) - 1, track.phi(), wacc_poi, 2);
      if (WithinPtPOI && WithinPtRef && pid_index)
        fGFW->Fill(track.eta(), fPtAxis->FindBin(track.pt()) - 1, track.phi(), wacc_poi, (1 << (pid_index + 5)));
      if (WithinPtNch && WithinPtRef)
        fGFW->Fill(track.eta(), fPtAxis->FindBin(track.pt()) - 1, track.phi(), wacc_poi, 32);
    } else { // Analysing only integrated flow
      double weff = (dt == kGen) ? 1. : getEfficiency(track);
      if (weff < 0)
        return;
      double wacc = (dt == kGen) ? 1. : getAcceptance(track, vtxz, -1);
      fGFW->Fill(track.eta(), fPtAxis->FindBin(track.pt()) - 1, track.phi(), weff * wacc, 1);
    }
    return;
  }

  template <datatype dt, QAFillTime ft, typename TTrack>
  inline void FillTrackQA(TTrack track, const float vtxz)
  {
    if constexpr (dt == kGen) {
      registry.fill(HIST("MCGen/") + HIST(moment[ft]) + HIST("phi_eta_vtxZ_gen"), track.phi(), track.eta(), vtxz);
      registry.fill(HIST("MCGen/") + HIST(moment[ft]) + HIST("pt_gen"), track.pt());
    } else {
      double wacc = getAcceptance(track, vtxz, -1);
      registry.fill(HIST("trackQA/") + HIST(moment[ft]) + HIST("phi_eta_vtxZ"), track.phi(), track.eta(), vtxz, wacc);
      registry.fill(HIST("trackQA/") + HIST(moment[ft]) + HIST("pt_dcaXY_dcaZ"), track.pt(), track.dcaXY(), track.dcaZ());
      if (ft == kAfter) {
        registry.fill(HIST("trackQA/") + HIST(moment[ft]) + HIST("pt_ref"), track.pt());
        registry.fill(HIST("trackQA/") + HIST(moment[ft]) + HIST("pt_poi"), track.pt());
      }
    }
  }

  template <QAFillTime ft, typename CollisionObject, typename TracksObject>
  inline void FillEventQA(CollisionObject collision, TracksObject tracks)
  {
    registry.fill(HIST("eventQA/") + HIST(moment[ft]) + HIST("globalTracks_centT0C"), collision.centFT0C(), tracks.size());
    registry.fill(HIST("eventQA/") + HIST(moment[ft]) + HIST("PVTracks_centT0C"), collision.centFT0C(), collision.multNTracksPV());
    registry.fill(HIST("eventQA/") + HIST(moment[ft]) + HIST("globalTracks_PVTracks"), collision.multNTracksPV(), tracks.size());
    registry.fill(HIST("eventQA/") + HIST(moment[ft]) + HIST("globalTracks_multT0A"), collision.multFT0A(), tracks.size());
    registry.fill(HIST("eventQA/") + HIST(moment[ft]) + HIST("globalTracks_multV0A"), collision.multFV0A(), tracks.size());
    registry.fill(HIST("eventQA/") + HIST(moment[ft]) + HIST("multV0A_multT0A"), collision.multFT0A(), collision.multFV0A());
    registry.fill(HIST("eventQA/") + HIST(moment[ft]) + HIST("multT0C_centT0C"), collision.centFT0C(), collision.multFT0C());
    return;
  }

  Filter collisionFilter = nabs(aod::collision::posZ) < cfgVtxZ;
  Filter trackFilter = nabs(aod::track::eta) < cfgEta && aod::track::pt > cfgPtmin&& aod::track::pt < cfgPtmax && ((requireGlobalTrackInFilter()) || (aod::track::isGlobalTrackSDD == (uint8_t) true)) && nabs(aod::track::dcaXY) < cfgDCAxy&& nabs(aod::track::dcaZ) < cfgDCAz;
  using myTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA, aod::pidTOFPi, aod::pidTPCPi, aod::pidTOFKa, aod::pidTPCKa, aod::pidTOFPr, aod::pidTPCPr>>;

  void processData(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Cs>>::iterator const& collision, aod::BCsWithTimestamps const&, myTracks const& tracks)
  {
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    int run = bc.runNumber();
    if (run != lastRun) {
      lastRun = run;
      if (cfgFillWeights && cfgRunByRunWeights)
        CreateRunByRunWeights();
    }
    registry.fill(HIST("eventQA/eventSel"), 0.5);
    if (!collision.sel8())
      return;
    registry.fill(HIST("eventQA/eventSel"), 1.5);

    if (cfgOccupancySelection != -999) {
      int occupancy = collision.trackOccupancyInTimeRange();
      if (occupancy < 0 || occupancy > cfgOccupancySelection)
        return;
    }
    registry.fill(HIST("eventQA/eventSel"), 2.5);
    const auto centrality = collision.centFT0C();

    if (cfgFillQA)
      FillEventQA<kBefore>(collision, tracks);
    if (cfgUseAdditionalEventCut && !eventSelected(collision, tracks.size(), centrality))
      return;
    if (cfgFillQA)
      FillEventQA<kAfter>(collision, tracks);
    if (!cfgFillWeights)
      loadCorrections(bc);
    auto field = (cfgMagField == 99999) ? getMagneticField(bc.timestamp()) : cfgMagField;
    processCollision<kReco>(collision, tracks, centrality, field);
  }
  PROCESS_SWITCH(GenericFramework, processData, "Process analysis for non-derived data", true);

  void processMCReco(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Cs>>::iterator const& collision, aod::BCsWithTimestamps const&, soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA, aod::McTrackLabels>> const& tracks, aod::McParticles const&)
  {
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    int run = bc.runNumber();
    if (run != lastRun) {
      lastRun = run;
      if (cfgFillWeights && cfgRunByRunWeights)
        CreateRunByRunWeights();
    }
    if (!collision.sel8())
      return;
    const auto centrality = collision.centFT0C();
    if (cfgFillQA)
      FillEventQA<kBefore>(collision, tracks);
    if (cfgUseAdditionalEventCut && !eventSelected(collision, tracks.size(), centrality))
      return;
    if (cfgFillQA)
      FillEventQA<kAfter>(collision, tracks);

    if (!cfgFillWeights)
      loadCorrections(bc);
    auto field = (cfgMagField == 99999) ? getMagneticField(bc.timestamp()) : cfgMagField;
    processCollision<kReco>(collision, tracks, centrality, field);
  }
  PROCESS_SWITCH(GenericFramework, processMCReco, "Process analysis for MC reconstructed events", false);

  Filter mcCollFilter = nabs(aod::mccollision::posZ) < cfgVtxZ;
  void processMCGen(soa::Filtered<aod::McCollisions>::iterator const& mcCollision, soa::SmallGroups<soa::Join<aod::McCollisionLabels, aod::Collisions, aod::CentFT0Cs>> const& collisions, aod::McParticles const& particles)
  {
    if (collisions.size() != 1)
      return;
    float centrality = -1;
    for (auto& collision : collisions) {
      centrality = collision.centFT0C();
    }
    processCollision<kGen>(mcCollision, particles, centrality, -999);
  }
  PROCESS_SWITCH(GenericFramework, processMCGen, "Process analysis for MC generated events", false);

  void processRun2(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentRun2V0Ms>>::iterator const& collision, aod::BCsWithTimestamps const&, myTracks const& tracks)
  {
    if (!collision.sel7())
      return;
    const auto centrality = collision.centRun2V0M();
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    loadCorrections(bc);
    auto field = (cfgMagField == 99999) ? getMagneticField(bc.timestamp()) : cfgMagField;
    processCollision<kReco>(collision, tracks, centrality, field);
  }
  PROCESS_SWITCH(GenericFramework, processRun2, "Process analysis for Run 2 converted data", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<GenericFramework>(cfgc),
  };
}
