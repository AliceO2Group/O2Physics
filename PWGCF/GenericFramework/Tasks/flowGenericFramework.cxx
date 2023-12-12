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
#include "GFWConfig.h"
#include "GFWWeights.h"
#include <TProfile.h>
#include <TRandom3.h>

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
float ptlow = 0.2, ptup = 3.0;
int etabins = 16;
float etalow = -0.8, etaup = 0.8;
int vtxZbins = 40;
float vtxZlow = -10.0, vtxZup = 10.0;
int phibins = 72;
float philow = 0.0;
float phiup = constants::math::TwoPI;
int nchbins = 300;
float nchlow = 0.5;
float nchup = 3000.5;
std::vector<double> centbinning(90);
int nBootstrap = 10;
GFWRegions regions;
GFWCorrConfigs configs;
} // namespace o2::analysis::genericframework

using namespace o2::analysis::genericframework;

struct GenericFramework {

  O2_DEFINE_CONFIGURABLE(cfgNbootstrap, int, 10, "Number of subsamples")
  O2_DEFINE_CONFIGURABLE(cfgUseNch, bool, false, "Do correlations as function of Nch")
  O2_DEFINE_CONFIGURABLE(cfgFillWeights, bool, false, "Fill NUA weights")
  O2_DEFINE_CONFIGURABLE(cfgFillQA, bool, false, "Fill QA histograms")
  O2_DEFINE_CONFIGURABLE(cfgEfficiency, std::string, "", "CCDB path to efficiency object")
  O2_DEFINE_CONFIGURABLE(cfgAcceptance, std::string, "", "CCDB path to acceptance object")

  Configurable<GFWBinningCuts> cfgBinning{"cfgBinning",
                                          {40, -10.0, 10.0, 0.2, 2.0, {0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.5, 4, 5, 6, 8, 10}, 16, -0.8, 0.8, 72, 0.2, 5.0, 300, 0.5, 3000.5, {0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90}},
                                          "triplets - nbins, min, max - for z_vtx, PtPOI, eta - nbins for phi - ptRef min and max"};

  Configurable<GFWRegions> cfgRegions{"cfgRegions", {{"refN", "refP", "refFull"}, {-0.8, 0.4, -0.8}, {-0.4, 0.8, 0.8}, {0, 0, 0}, {1, 1, 1}}, "Configurations for GFW regions"};

  Configurable<GFWCorrConfigs> cfgCorrConfig{"cfgCorrConfig", {{"refP {2} refN {-2}", "refP {3} refN {-3}", "refP {4} refN {-4}", "refFull {2 -2}", "refFull {2 2 -2 -2}"}, {"ChGap22", "ChGap32", "ChGap42", "ChFull22", "ChFull24"}, {0, 0, 0, 0, 0}}, "Configurations for each correlation to calculate"};

#include "PWGCF/TwoParticleCorrelations/TableProducer/Productions/skimmingconf_20221115.cxx" // NOLINT
  // Connect to ccdb
  Service<ccdb::BasicCCDBManager> ccdb;

  struct Config {
    TH1D* mEfficiency = nullptr;
    GFWWeights* mAcceptance = nullptr;
    bool correctionsLoaded = false;
  } cfg;

  // Define output
  OutputObj<FlowContainer> fFC{FlowContainer("FlowContainer")};
  OutputObj<FlowContainer> fFC_gen{FlowContainer("FlowContainer_gen")};
  OutputObj<FlowContainer> fFC_reco{FlowContainer("FlowContainer_reco")};
  OutputObj<GFWWeights> fWeights{GFWWeights("weights")};
  HistogramRegistry registry{"registry"};
  HistogramRegistry QAregistry{"QAregistry", {}, OutputObjHandlingPolicy::QAObject};

  // define global variables
  GFW* fGFW = new GFW();
  std::vector<GFW::CorrConfig> corrconfigs;
  TRandom3* fRndm = new TRandom3(0);
  TAxis* fPtAxis;

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
    LOGF(info, "Regions:\n%s", regions.Print());
    LOGF(info, "CorrConfigs:\n%s", configs.Print());
    ptbinning = cfgBinning->GetPtBinning();
    ptpoilow = cfgBinning->GetPtPOImin();
    ptpoiup = cfgBinning->GetPtPOImax();
    ptreflow = cfgBinning->GetPtRefMin();
    ptrefup = cfgBinning->GetPtRefMax();
    ptlow = static_cast<float>(std::min(ptpoilow, ptreflow));
    ptup = static_cast<float>(std::max(ptpoiup, ptrefup));
    etabins = cfgBinning->GetEtaBins();
    etalow = cfgBinning->GetEtaMin();
    etaup = cfgBinning->GetEtaMax();
    vtxZbins = cfgBinning->GetVtxZbins();
    vtxZlow = cfgBinning->GetVtxZmin();
    vtxZup = cfgBinning->GetVtxZmax();
    phibins = cfgBinning->GetPhiBins();
    philow = 0.0f;
    phiup = constants::math::TwoPI;
    nchbins = cfgBinning->GetNchBins();
    nchlow = cfgBinning->GetNchMin();
    nchup = cfgBinning->GetNchMax();
    centbinning = cfgBinning->GetCentBinning();

    AxisSpec phiAxis = {phibins, philow, phiup, "#phi"};
    AxisSpec etaAxis = {etabins, etalow, etaup, "#eta"};
    AxisSpec vtxAxis = {vtxZbins, vtxZlow, vtxZup, "Vtx_{z} (cm)"};
    AxisSpec ptAxis = {ptbinning, "#it{p}_{T} GeV/#it{c}"};
    AxisSpec centAxis = {centbinning, "Centrality (%)"};
    AxisSpec nchAxis = {nchbins, nchlow, nchup};
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

    if (doprocessGen || doprocessReco) {
      if (doprocessGen) {
        QAregistry.add("pt_gen", "", {HistType::kTH1D, {ptAxis}});
        QAregistry.add("phi_gen", "", {HistType::kTH1D, {phiAxis}});
        QAregistry.add("eta_gen", "", {HistType::kTH1D, {etaAxis}});
        QAregistry.add("vtxZ_gen", "", {HistType::kTH1D, {vtxAxis}});
      }
      if (doprocessReco) {
        QAregistry.add("pt_reco", "", {HistType::kTH1D, {ptAxis}});
        QAregistry.add("phi_reco", "", {HistType::kTH1D, {phiAxis}});
        QAregistry.add("eta_reco", "", {HistType::kTH1D, {etaAxis}});
        QAregistry.add("vtxZ_reco", "", {HistType::kTH1D, {vtxAxis}});
      }
    } else {
      QAregistry.add("phi", "", {HistType::kTH1D, {phiAxis}});
      QAregistry.add("eta", "", {HistType::kTH1D, {etaAxis}});
      QAregistry.add("vtxZ", "", {HistType::kTH1D, {vtxAxis}});

      QAregistry.add("pt_dcaXY_dcaZ", "", {HistType::kTH3D, {ptAxis, dcaXYAXis, dcaZAXis}});
      registry.add("phi_eta_vtxZ_corrected", "", {HistType::kTH3D, {phiAxis, etaAxis, vtxAxis}});
      QAregistry.add("cent_nch", "", {HistType::kTH2D, {nchAxis, centAxis}});
      QAregistry.add("globalTracks_PVTracks", "", {HistType::kTH2D, {nchAxis, nchAxis}});
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
    if (!doprocessGen && !doprocessReco) {
      fFC->SetName("FlowContainer");
      fFC->SetXAxis(fPtAxis);
      fFC->Initialize(oba, multAxis, cfgNbootstrap);
    }
    if (doprocessGen) {
      fFC_gen->SetName("FlowContainer_gen");
      fFC_gen->SetXAxis(fPtAxis);
      fFC_gen->Initialize(oba, multAxis, cfgNbootstrap);
    }
    if (doprocessReco) {
      fFC_reco->SetName("FlowContainer_reco");
      fFC_reco->SetXAxis(fPtAxis);
      fFC_reco->Initialize(oba, multAxis, cfgNbootstrap);
    }
    delete oba;
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

  void FillFC(OutputObj<FlowContainer> fc, const GFW::CorrConfig& corrconf, const double& cent, const double& rndm)
  {
    double dnx, val;
    dnx = fGFW->Calculate(corrconf, 0, kTRUE).real();
    if (dnx == 0)
      return;
    if (!corrconf.pTDif) {
      val = fGFW->Calculate(corrconf, 0, kFALSE).real() / dnx;
      if (TMath::Abs(val) < 1)
        fc->FillProfile(corrconf.Head.c_str(), cent, val, dnx, rndm);
      return;
    }
    for (Int_t i = 1; i <= fPtAxis->GetNbins(); i++) {
      dnx = fGFW->Calculate(corrconf, i - 1, kTRUE).real();
      if (dnx == 0)
        continue;
      val = fGFW->Calculate(corrconf, i - 1, kFALSE).real() / dnx;
      if (TMath::Abs(val) < 1)
        fc->FillProfile(Form("%s_pt_%i", corrconf.Head.c_str(), i), cent, val, dnx, rndm);
    }
    return;
  }

  Filter collisionFilter = aod::collision::posZ < vtxZup && aod::collision::posZ > vtxZlow;
  Filter trackFilter = aod::track::eta < etaup && aod::track::eta > etalow&& aod::track::pt > ptlow&& aod::track::pt < ptup && ((requireGlobalTrackInFilter()) || (aod::track::isGlobalTrackSDD == (uint8_t) true)) && nabs(aod::track::dcaXY) < 0.2f;
  using myTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA>>;

  enum datatype {
    kRaw,
    kDerived
  };

  template <datatype dt, typename TCollision, typename TTracks>
  void processData(TCollision collision, TTracks tracks, float centrality)
  {
    if (tracks.size() < 1)
      return;
    if (centrality < centbinning.front() || centrality > centbinning.back())
      return;
    float vtxz = collision.posZ();

    if (cfgFillQA) {
      QAregistry.fill(HIST("cent_nch"), tracks.size(), centrality);
      QAregistry.fill(HIST("vtxZ"), vtxz);
    }
    fGFW->Clear();
    float l_Random = fRndm->Rndm();
    float weff = 1, wacc = 1;

    for (auto& track : tracks) {
      if (cfgFillWeights)
        fWeights->Fill(track.phi(), track.eta(), vtxz, track.pt(), centrality, 0);
      if (cfg.mEfficiency)
        weff = cfg.mEfficiency->GetBinContent(cfg.mEfficiency->FindBin(track.pt()));
      else
        weff = 1.0;
      if (weff == 0)
        continue;
      weff = 1. / weff;
      if (cfg.mAcceptance)
        wacc = cfg.mAcceptance->GetNUA(track.phi(), track.eta(), vtxz);
      else
        wacc = 1;
      registry.fill(HIST("phi_eta_vtxZ_corrected"), track.phi(), track.eta(), vtxz, wacc);
      bool WithinPtPOI = (ptpoilow < track.pt()) && (track.pt() < ptpoiup); // within POI pT range
      bool WithinPtRef = (ptreflow < track.pt()) && (track.pt() < ptrefup); // within RF pT range
      if (WithinPtRef)
        fGFW->Fill(track.eta(), fPtAxis->FindBin(track.pt() - 1), track.phi(), wacc * weff, 1);
      if (WithinPtPOI)
        fGFW->Fill(track.eta(), fPtAxis->FindBin(track.pt()) - 1, track.phi(), wacc * weff, 2);
      if (WithinPtPOI && WithinPtRef)
        fGFW->Fill(track.eta(), fPtAxis->FindBin(track.pt()) - 1, track.phi(), wacc * weff, 4);

      if (cfgFillQA) {
        QAregistry.fill(HIST("phi"), track.phi());
        QAregistry.fill(HIST("eta"), track.eta());
        QAregistry.fill(HIST("pt_dcaXY_dcaZ"), track.pt(), track.dcaXY(), track.dcaZ());
      }
    }
    for (uint l_ind = 0; l_ind < corrconfigs.size(); l_ind++) {
      FillFC(fFC, corrconfigs.at(l_ind), centrality, l_Random);
    }
  }

  void processRawData(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Cs>>::iterator const& collision, aod::BCsWithTimestamps const&, myTracks const& tracks)
  {
    const auto centrality = collision.centFT0C();
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    if (!collision.sel8())
      return;
    QAregistry.fill(HIST("globalTracks_PVTracks"), collision.multNTracksPV(), tracks.size());
    loadCorrections(bc.timestamp());
    processData<kRaw>(collision, tracks, centrality);
  }
  PROCESS_SWITCH(GenericFramework, processRawData, "Process analysis for derived data", true);

  void processGen(soa::Filtered<soa::Join<aod::Collisions, aod::CentFT0Cs, aod::McCollisionLabels>>::iterator const& collision, aod::McParticles const& particles)
  {
    if (collision.has_mcCollision()) {
      float centrality = collision.centFT0C();
      if (particles.size() < 1)
        return;
      if (centrality < centbinning.front() || centrality > centbinning.back())
        return;
      float vtxz = collision.posZ();
      if (cfgFillQA) {
        QAregistry.fill(HIST("vtxZ_gen"), vtxz);
      }
      fGFW->Clear();
      float l_Random = fRndm->Rndm();

      for (auto& particle : particles) {
        if (!particle.isPhysicalPrimary() || particle.eta() < etalow || particle.eta() > etaup || particle.pt() < ptlow || particle.pt() > ptup)
          continue;
        bool WithinPtPOI = (ptpoilow < particle.pt()) && (particle.pt() < ptpoiup); // within POI pT range
        bool WithinPtRef = (ptreflow < particle.pt()) && (particle.pt() < ptrefup); // within RF pT range
        if (WithinPtRef)
          fGFW->Fill(particle.eta(), fPtAxis->FindBin(particle.pt()) - 1, particle.phi(), 1, 1);
        if (WithinPtPOI)
          fGFW->Fill(particle.eta(), fPtAxis->FindBin(particle.pt()) - 1, particle.phi(), 1, 2);
        if (WithinPtPOI && WithinPtRef)
          fGFW->Fill(particle.eta(), fPtAxis->FindBin(particle.pt()) - 1, particle.phi(), 1, 4);
        if (cfgFillQA) {
          QAregistry.fill(HIST("phi_gen"), particle.phi());
          QAregistry.fill(HIST("eta_gen"), particle.eta());
          QAregistry.fill(HIST("pt_gen"), particle.pt());
        }
      }
      for (uint l_ind = 0; l_ind < corrconfigs.size(); l_ind++) {
        FillFC(fFC_gen, corrconfigs.at(l_ind), centrality, l_Random);
      }
    }
  }
  PROCESS_SWITCH(GenericFramework, processGen, "Process analysis for MC generated events", false);

  Filter dcaFilter = nabs(aod::track::dcaXY) < 0.2f;
  void processReco(soa::Filtered<soa::Join<aod::Collisions, aod::CentFT0Cs>>::iterator const& collision, soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::McTrackLabels>> const& tracks, aod::McParticles const&)
  {
    float centrality = collision.centFT0C();
    if (tracks.size() < 1)
      return;
    if (centrality < centbinning.front() || centrality > centbinning.back())
      return;
    float vtxz = collision.posZ();
    if (cfgFillQA) {
      QAregistry.fill(HIST("vtxZ_reco"), vtxz);
    }
    fGFW->Clear();
    float l_Random = fRndm->Rndm();

    for (auto& track : tracks) {
      if (track.has_mcParticle()) {
        auto particle = track.mcParticle();
        if (!particle.isPhysicalPrimary() || particle.eta() < etalow || particle.eta() > etaup || particle.pt() < ptlow || particle.pt() > ptup)
          continue;
        bool WithinPtPOI = (ptpoilow < particle.pt()) && (particle.pt() < ptpoiup); // within POI pT range
        bool WithinPtRef = (ptreflow < particle.pt()) && (particle.pt() < ptrefup); // within RF pT range
        if (WithinPtRef)
          fGFW->Fill(particle.eta(), fPtAxis->FindBin(particle.pt()) - 1, particle.phi(), 1, 1);
        if (WithinPtPOI)
          fGFW->Fill(particle.eta(), fPtAxis->FindBin(particle.pt()) - 1, particle.phi(), 1, 2);
        if (WithinPtPOI && WithinPtRef)
          fGFW->Fill(particle.eta(), fPtAxis->FindBin(particle.pt()) - 1, particle.phi(), 1, 4);
        if (cfgFillQA) {
          QAregistry.fill(HIST("phi_reco"), particle.phi());
          QAregistry.fill(HIST("eta_reco"), particle.eta());
          QAregistry.fill(HIST("pt_reco"), particle.pt());
        }
      }
    }
    for (uint l_ind = 0; l_ind < corrconfigs.size(); l_ind++) {
      FillFC(fFC_reco, corrconfigs.at(l_ind), centrality, l_Random);
    }
  }
  PROCESS_SWITCH(GenericFramework, processReco, "Process analysis for MC reconstructed events", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<GenericFramework>(cfgc),
  };
}
