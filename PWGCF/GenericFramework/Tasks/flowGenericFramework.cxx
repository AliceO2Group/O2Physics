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

namespace o2::analysis::genericframework
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
  O2_DEFINE_CONFIGURABLE(cfgFillQA, bool, false, "Fill QA histograms")
  O2_DEFINE_CONFIGURABLE(cfgUse22sEventCut, bool, false, "Use 22s event cut on mult correlations")
  O2_DEFINE_CONFIGURABLE(cfgEfficiency, std::string, "", "CCDB path to efficiency object")
  O2_DEFINE_CONFIGURABLE(cfgAcceptance, std::string, "", "CCDB path to acceptance object")
  O2_DEFINE_CONFIGURABLE(cfgDCAxy, float, 0.2, "Cut on DCA in the transverse direction");
  O2_DEFINE_CONFIGURABLE(cfgPtmin, float, 0.2, "minimum pt");
  O2_DEFINE_CONFIGURABLE(cfgPtmax, float, 10, "maximum pt");

  Configurable<GFWBinningCuts> cfgGFWBinning{"cfgGFWBinning", {40, -10.0, 10.0, 16, -0.8, 0.8, 72, 300, 0.5, 3000.5, 0.2, 10.0, 0.2, 3.0, {0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.5, 4, 5, 6, 8, 10}, {0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90}}, "Configuration for binning"};
  Configurable<GFWRegions> cfgRegions{"cfgRegions", {{"refN", "refP", "refFull"}, {-0.8, 0.4, -0.8}, {-0.4, 0.8, 0.8}, {0, 0, 0}, {1, 1, 1}}, "Configurations for GFW regions"};

  Configurable<GFWCorrConfigs> cfgCorrConfig{"cfgCorrConfig", {{"refP {2} refN {-2}", "refP {3} refN {-3}", "refP {4} refN {-4}", "refFull {2 -2}", "refFull {2 2 -2 -2}"}, {"ChGap22", "ChGap32", "ChGap42", "ChFull22", "ChFull24"}, {0, 0, 0, 0, 0}, {15, 1, 1, 0, 0}}, "Configurations for each correlation to calculate"};

  // #include "PWGCF/TwoParticleCorrelations/TableProducer/Productions/skimmingconf_20221115.cxx" // NOLINT
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
    etalow = cfgGFWBinning->GetEtaMin();
    etaup = cfgGFWBinning->GetEtaMax();
    vtxZbins = cfgGFWBinning->GetVtxZbins();
    vtxZlow = cfgGFWBinning->GetVtxZmin();
    vtxZup = cfgGFWBinning->GetVtxZmax();
    phibins = cfgGFWBinning->GetPhiBins();
    philow = 0.0f;
    phiup = constants::math::TwoPI;
    nchbins = cfgGFWBinning->GetNchBins();
    nchlow = cfgGFWBinning->GetNchMin();
    nchup = cfgGFWBinning->GetNchMax();
    centbinning = cfgGFWBinning->GetCentBinning();
    cfgGFWBinning->Print();

    AxisSpec phiAxis = {phibins, philow, phiup, "#phi"};
    AxisSpec etaAxis = {etabins, etalow, etaup, "#eta"};
    AxisSpec vtxAxis = {vtxZbins, vtxZlow, vtxZup, "Vtx_{z} (cm)"};
    AxisSpec ptAxis = {ptbinning, "#it{p}_{T} GeV/#it{c}"};
    AxisSpec centAxis = {centbinning, "Centrality (%)"};
    std::vector<double> nchbinning;
    int nchskip = (nchup - nchlow) / nchbins;
    for (int i = 0; i <= nchbins; ++i) {
      nchbinning.push_back(nchskip * i + nchlow + 0.5);
    }
    AxisSpec nchAxis = {nchbinning, "N_{ch}"};
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

    if (cfgFillWeights) {
      fWeights->SetPtBins(ptbins, &ptbinning[0]);
      fWeights->Init(true, false);
    }

    if (doprocessMCGen) {
      registry.add("pt_gen", "", {HistType::kTH1D, {ptAxis}});
      registry.add("phi_gen", "", {HistType::kTH1D, {phiAxis}});
      registry.add("eta_gen", "", {HistType::kTH1D, {etaAxis}});
      registry.add("vtxZ_gen", "", {HistType::kTH1D, {vtxAxis}});
    }
    if (doprocessMCReco || doprocessData || doprocessRun2) {
      registry.add("phi", "", {HistType::kTH1D, {phiAxis}});
      registry.add("eta", "", {HistType::kTH1D, {etaAxis}});
      registry.add("vtxZ", "", {HistType::kTH1D, {vtxAxis}});
      registry.add("pt_dcaXY_dcaZ", "", {HistType::kTH3D, {ptAxis, dcaXYAXis, dcaZAXis}});
      registry.add("phi_eta_vtxZ_corrected", "", {HistType::kTH3D, {phiAxis, etaAxis, vtxAxis}});
      registry.add("cent_nch", "", {HistType::kTH2D, {nchAxis, centAxis}});
      registry.add("globalTracks_PVTracks", "", {HistType::kTH2D, {nchAxis, nchAxis}});
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
    fFCpt->Initialise(multAxis, cfgMpar, configs, cfgNbootstrap);
    // Event selection - Alex
    if (cfgUse22sEventCut) {
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
    }
  }

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
      weight_nua = cfg.mAcceptance->GetNUA(phi, eta, vtxz);
    else
      weight_nua = 1;
    return true;
  }

  template <typename TCollision>
  bool eventSelected(TCollision collision, const int& multTrk, const float& centrality)
  {
    if (collision.alias_bit(kTVXinTRD)) {
      // TRD triggered
      return 0;
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
    if (multTrk > fMultMultPVCut->Eval(multNTracksPV))
      return 0;

    return 1;
  }

  enum datatype {
    kReco,
    kGen
  };

  template <datatype dt>
  void FillOutputContainers(const GFW::CorrConfig& corrconf, const float& centmult, const double& rndm, uint8_t mask)
  {
    fFCpt->CalculateCorrelations();
    fFCpt->FillPtProfiles(centmult, rndm);
    fFCpt->FillCMProfiles(centmult, rndm);
    auto dnx = fGFW->Calculate(corrconf, 0, kTRUE).real();
    if (dnx == 0)
      return;
    if (!corrconf.pTDif) {
      auto val = fGFW->Calculate(corrconf, 0, kFALSE).real() / dnx;
      if (TMath::Abs(val) < 1) {
        (dt == kGen) ? fFC_gen->FillProfile(corrconf.Head.c_str(), centmult, val, dnx, rndm) : fFC->FillProfile(corrconf.Head.c_str(), centmult, val, dnx, rndm);
        fFCpt->FillVnPtProfiles(centmult, val, dnx, rndm, mask);
      }
      return;
    }
    for (Int_t i = 1; i <= fPtAxis->GetNbins(); i++) {
      dnx = fGFW->Calculate(corrconf, i - 1, kTRUE).real();
      if (dnx == 0)
        continue;
      auto val = fGFW->Calculate(corrconf, i - 1, kFALSE).real() / dnx;
      if (TMath::Abs(val) < 1)
        (dt == kGen) ? fFC_gen->FillProfile(Form("%s_pt_%i", corrconf.Head.c_str(), i), centmult, val, dnx, rndm) : fFC->FillProfile(Form("%s_pt_%i", corrconf.Head.c_str(), i), centmult, val, dnx, rndm);
    }
    return;
  }

  template <datatype dt, typename TCollision, typename TTracks>
  void processCollision(TCollision collision, TTracks tracks, const float& centrality)
  {
    if (tracks.size() < 1)
      return;
    if (centrality < centbinning.front() || centrality > centbinning.back())
      return;
    float vtxz = collision.posZ();
    if (cfgFillQA)
      (dt == kGen) ? registry.fill(HIST("vtxZ_gen"), vtxz) : registry.fill(HIST("vtxZ"), vtxz);
    fGFW->Clear();
    fFCpt->ClearVector();
    float l_Random = fRndm->Rndm();
    for (auto& track : tracks) {
      ProcessTrack(track, centrality, vtxz);
    }
    for (uint l_ind = 0; l_ind < corrconfigs.size(); ++l_ind) {
      FillOutputContainers<dt>(corrconfigs.at(l_ind), (cfgUseNch) ? tracks.size() : centrality, l_Random, configs.GetpTCorrMasks()[l_ind]);
    }
  }

  template <typename TrackObject>
  inline void ProcessTrack(TrackObject const& track, const float& centrality, const float& vtxz)
  {
    float weff = 1, wacc = 1;
    if constexpr (framework::has_type_v<aod::mctracklabel::McParticleId, typename TrackObject::all_columns>) {
      if (track.mcParticleId() < 0 || !(track.has_mcParticle()))
        return;
      auto mcParticle = track.mcParticle();
      if (!mcParticle.isPhysicalPrimary() || mcParticle.eta() < etalow || mcParticle.eta() > etaup || mcParticle.pt() < ptlow || mcParticle.pt() > ptup)
        return;
      if (cfgFillWeights)
        fWeights->Fill(mcParticle.phi(), mcParticle.eta(), vtxz, mcParticle.pt(), centrality, 0);
      if (!setCurrentParticleWeights(weff, wacc, mcParticle.phi(), mcParticle.eta(), mcParticle.pt(), vtxz))
        return;
      registry.fill(HIST("phi_eta_vtxZ_corrected"), mcParticle.phi(), mcParticle.eta(), vtxz, wacc);
      if (cfgFillQA)
        FillQA<kReco>(track);
      FillGFW(mcParticle, weff, wacc);
    } else if constexpr (framework::has_type_v<aod::mcparticle::McCollisionId, typename TrackObject::all_columns>) {
      if (!track.isPhysicalPrimary() || track.eta() < etalow || track.eta() > etaup || track.pt() < ptlow || track.pt() > ptup)
        return;
      if (cfgFillQA)
        FillQA<kGen>(track);
      FillGFW(track, 1., 1.);
    } else {
      if (cfgFillWeights)
        fWeights->Fill(track.phi(), track.eta(), vtxz, track.pt(), centrality, 0);
      if (!setCurrentParticleWeights(weff, wacc, track.phi(), track.eta(), track.pt(), vtxz))
        return;
      registry.fill(HIST("phi_eta_vtxZ_corrected"), track.phi(), track.eta(), vtxz, wacc);
      if (cfgFillQA)
        FillQA<kReco>(track);
      FillGFW(track, weff, wacc);
    }
  }

  template <typename TrackObject>
  inline void FillGFW(TrackObject track, float weff, float wacc)
  {
    fFCpt->Fill(weff, track.pt());
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

  template <datatype dt, typename TrackObject>
  inline void FillQA(TrackObject track)
  {
    if constexpr (dt == kGen) {
      registry.fill(HIST("phi_gen"), track.phi());
      registry.fill(HIST("eta_gen"), track.eta());
      registry.fill(HIST("pt_gen"), track.pt());
    } else {
      registry.fill(HIST("phi"), track.phi());
      registry.fill(HIST("eta"), track.eta());
      registry.fill(HIST("pt_dcaXY_dcaZ"), track.pt(), track.dcaXY(), track.dcaZ());
    }
  }

  Filter collisionFilter = aod::collision::posZ < cfgGFWBinning->GetVtxZmax() && aod::collision::posZ > cfgGFWBinning->GetVtxZmin();
  Filter trackFilter = aod::track::eta < cfgGFWBinning->GetEtaMax() && aod::track::eta > cfgGFWBinning->GetEtaMin() && aod::track::pt > cfgPtmin&& aod::track::pt < cfgPtmax && ((requireGlobalTrackInFilter()) || (aod::track::isGlobalTrackSDD == (uint8_t) true)) && nabs(aod::track::dcaXY) < cfgDCAxy;
  using myTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA>>;

  void processData(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Cs>>::iterator const& collision, aod::BCsWithTimestamps const&, myTracks const& tracks)
  {
    if (!collision.sel8())
      return;
    const auto centrality = collision.centFT0C();
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    if (cfgUse22sEventCut && !eventSelected(collision, tracks.size(), centrality))
      return;
    registry.fill(HIST("globalTracks_PVTracks"), collision.multNTracksPV(), tracks.size());
    registry.fill(HIST("cent_nch"), tracks.size(), centrality);
    loadCorrections(bc.timestamp());
    processCollision<kReco>(collision, tracks, centrality);
  }
  PROCESS_SWITCH(GenericFramework, processData, "Process analysis for non-derived data", true);

  void processMCReco(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Cs>>::iterator const& collision, aod::BCsWithTimestamps const&, soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA, aod::McTrackLabels>> const& tracks, aod::McParticles const&)
  {
    if (!collision.sel8())
      return;
    const auto centrality = collision.centFT0C();
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    if (cfgUse22sEventCut && !eventSelected(collision, tracks.size(), centrality))
      return;
    registry.fill(HIST("globalTracks_PVTracks"), collision.multNTracksPV(), tracks.size());
    registry.fill(HIST("cent_nch"), tracks.size(), centrality);
    loadCorrections(bc.timestamp());
    processCollision<kReco>(collision, tracks, centrality);
  }
  PROCESS_SWITCH(GenericFramework, processMCReco, "Process analysis for MC reconstructed events", false);

  Filter mcCollFilter = aod::mccollision::posZ < cfgGFWBinning->GetVtxZmax() && aod::mccollision::posZ > cfgGFWBinning->GetVtxZmin();
  void processMCGen(soa::Filtered<aod::McCollisions>::iterator const& mcCollision, soa::SmallGroups<soa::Join<aod::McCollisionLabels, aod::Collisions, aod::CentFT0Cs>> const& collisions, aod::McParticles const& particles)
  {
    if (collisions.size() != 1)
      return;
    float centrality = -1;
    for (auto& collision : collisions) {
      centrality = collision.centFT0C();
    }
    processCollision<kGen>(mcCollision, particles, centrality);
  }
  PROCESS_SWITCH(GenericFramework, processMCGen, "Process analysis for MC generated events", false);

  void processRun2(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentRun2V0Ms>>::iterator const& collision, aod::BCsWithTimestamps const&, myTracks const& tracks)
  {
    if (!collision.sel7())
      return;
    const auto centrality = collision.centRun2V0M();
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    if (cfgUse22sEventCut && !eventSelected(collision, tracks.size(), centrality))
      return;
    registry.fill(HIST("globalTracks_PVTracks"), collision.multNTracksPV(), tracks.size());
    registry.fill(HIST("cent_nch"), tracks.size(), centrality);
    loadCorrections(bc.timestamp());
    processCollision<kReco>(collision, tracks, centrality);
  }
  PROCESS_SWITCH(GenericFramework, processRun2, "Process analysis for Run 2 converted data", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<GenericFramework>(cfgc),
  };
}
