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
/// \author Antonio Ortiz (antonio.ortiz@nucleares.unam.mx)
/// \since November 2021
/// \last update: July 2024

#include <cmath>
#include <vector>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/StaticFor.h"
#include "Framework/O2DatabasePDGPlugin.h"

#include "ReconstructionDataFormats/Track.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"

#include "TDatabasePDG.h"
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TRandom.h>

// TODO: implement 50% stat for MC closure vs 50% for testing, add flag for weak decays

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct ueCharged {

  TrackSelection myTrackSelectionPrim()
  {
    TrackSelection selectedTracks;
    selectedTracks.SetPtRange(0.1f, 1e10f);
    selectedTracks.SetEtaRange(-0.8f, 0.8f);
    selectedTracks.SetRequireITSRefit(true);
    selectedTracks.SetRequireTPCRefit(true);
    // selectedTracks.SetRequireGoldenChi2(true);
    selectedTracks.SetMinNCrossedRowsTPC(70);
    selectedTracks.SetMinNCrossedRowsOverFindableClustersTPC(0.8f);
    selectedTracks.SetMaxChi2PerClusterTPC(4.f);
    selectedTracks.SetRequireHitsInITSLayers(1, {0, 1}); // one hit in any SPD layer
    selectedTracks.SetMaxChi2PerClusterITS(36.f);
    selectedTracks.SetMaxDcaXYPtDep([](float pt) { return 0.0105f + 0.0350f / pow(pt, 1.1f); });
    selectedTracks.SetMaxDcaZ(2.f);
    return selectedTracks;
  }
  TrackSelection myTrackSelectionOpenDCA()
  {
    TrackSelection selectedTracks;
    selectedTracks.SetPtRange(0.1f, 1e10f);
    selectedTracks.SetEtaRange(-0.8f, 0.8f);
    selectedTracks.SetRequireITSRefit(true);
    selectedTracks.SetRequireTPCRefit(true);
    // selectedTracks.SetRequireGoldenChi2(true);
    selectedTracks.SetMinNCrossedRowsTPC(70);
    selectedTracks.SetMinNCrossedRowsOverFindableClustersTPC(0.8f);
    selectedTracks.SetMaxChi2PerClusterTPC(4.f);
    selectedTracks.SetRequireHitsInITSLayers(1, {0, 1}); // one hit in any SPD layer
    selectedTracks.SetMaxChi2PerClusterITS(36.f);
    // selectedTracks.SetMaxDcaXYPtDep([](float pt) { return 0.0105f + 0.0350f / pow(pt, 1.1f); });
    selectedTracks.SetMaxDcaZ(2.f);
    return selectedTracks;
  }

  TrackSelection mySelectionPrim;
  TrackSelection mySelectionOpenDCA;

  Service<o2::framework::O2DatabasePDG> pdg;
  float DeltaPhi(float phia, float phib, float rangeMin, float rangeMax);
  // Configurable for event selection
  Configurable<bool> isRun3{"isRun3", true, "is Run3 dataset"};
  Configurable<bool> piluprejection{"piluprejection", true, "Pileup rejection"};
  Configurable<bool> goodzvertex{"goodzvertex", true, "removes collisions with large differences between z of PV by tracks and z of PV from FT0 A-C time difference."};
  Configurable<bool> sel8{"sel8", true, "Apply the sel8 event selection"};
  Configurable<bool> removeITSROFBorder{"removeITSROFBorder", false, "Remove ITS Read-Out Frame border and only apply kIsTriggerTVX & kNoTimeFrameBorder (recommended for MC)"};
  Configurable<bool> manuallyApplysel8{"manuallyApplysel8", false, "Apply manually the event selection criteria considered in sel8, ie. kIsTriggerTVX & kNoTimeFrameBorder & kNoITSROFrameBorder"};

  // acceptance cuts
  Configurable<float> cfgTrkEtaCut{"cfgTrkEtaCut", 0.8f, "Eta range for tracks"};
  Configurable<float> cfgTrkLowPtCut{"cfgTrkLowPtCut", 0.15f, "Minimum constituent pT"};

  HistogramRegistry ue;
  static constexpr std::string_view pNumDenMeasuredPS[3] = {"pNumDenMeasuredPS_NS", "pNumDenMeasuredPS_AS", "pNumDenMeasuredPS_TS"};
  static constexpr std::string_view pSumPtMeasuredPS[3] = {"pSumPtMeasuredPS_NS", "pSumPtMeasuredPS_AS", "pSumPtMeasuredPS_TS"};
  static constexpr std::string_view hPhi[3] = {"hPhi_NS", "hPhi_AS", "hPhi_TS"};
  // data driven correction
  static constexpr std::string_view hNumDenMCDd[3] = {"hNumDenMCDd_NS", "hNumDenMCDd_AS", "hNumDenMCDd_TS"};
  static constexpr std::string_view hSumPtMCDd[3] = {"hSumPtMCDd_NS", "hSumPtMCDd_AS", "hSumPtMCDd_TS"};
  static constexpr std::string_view hNumDenMCMatchDd[3] = {"hNumDenMCMatchDd_NS", "hNumDenMCMatchDd_AS", "hNumDenMCMatchDd_TS"};
  static constexpr std::string_view hSumPtMCMatchDd[3] = {"hSumPtMCMatchDd_NS", "hSumPtMCMatchDd_AS", "hSumPtMCMatchDd_TS"};
  // hist data for corrections
  static constexpr std::string_view hPtVsPtLeadingData[3] = {"hPtVsPtLeadingData_NS", "hPtVsPtLeadingData_AS", "hPtVsPtLeadingData_TS"};
  static constexpr std::string_view pNumDenData[3] = {"pNumDenData_NS", "pNumDenData_AS", "pNumDenData_TS"};
  static constexpr std::string_view pSumPtData[3] = {"pSumPtData_NS", "pSumPtData_AS", "pSumPtData_TS"};
  // hist data true
  static constexpr std::string_view hPtVsPtLeadingTrue[3] = {"hPtVsPtLeadingTrue_NS", "hPtVsPtLeadingTrue_AS", "hPtVsPtLeadingTrue_TS"};
  // all wo detector effects
  static constexpr std::string_view pNumDenTrueAll[3] = {"pNumDenTrueAll_NS", "pNumDenTrueAll_AS", "pNumDenTrueAll_TS"};
  static constexpr std::string_view pSumPtTrueAll[3] = {"pSumPtTrueAll_NS", "pSumPtTrueAll_AS", "pSumPtTrueAll_TS"};
  // true, 50%
  static constexpr std::string_view pNumDenTrue[3] = {"pNumDenTrue_NS", "pNumDenTrue_AS", "pNumDenTrue_TS"};
  static constexpr std::string_view pSumPtTrue[3] = {"pSumPtTrue_NS", "pSumPtTrue_AS", "pSumPtTrue_TS"};

  // this must have all event selection effects, but it has not been implemented 50%
  static constexpr std::string_view pNumDenTruePS[3] = {"pNumDenTruePS_NS", "pNumDenTruePS_AS", "pNumDenTruePS_TS"};
  static constexpr std::string_view pSumPtTruePS[3] = {"pSumPtTruePS_NS", "pSumPtTruePS_AS", "pSumPtTruePS_TS"};
  static constexpr std::string_view hPhiTrue[3] = {"hPhiTrue_NS", "hPhiTrue_AS", "hPhiTrue_TS"};

  OutputObj<TF1> f_Eff{"fpara"};
  void init(InitContext const&);

  template <bool IS_MC, typename C, typename T, typename P>
  void processMeasMC(const C& collision, const T& tracks, const P& particles);

  template <typename C, typename T>
  void processMeas(const C& collision, const T& tracks);

  template <typename C, typename P>
  void processTrue(const C& mcCollision, const P& particles);

  Filter trackFilter = (nabs(aod::track::eta) < cfgTrkEtaCut) && (aod::track::pt > cfgTrkLowPtCut);

  using CollisionTableMCTrue = aod::McCollisions;
  using CollisionTableMC = soa::SmallGroups<soa::Join<aod::McCollisionLabels, aod::Collisions, aod::EvSels>>;
  using TrackTableMC = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::McTrackLabels, aod::TrackSelection>>;
  using ParticleTableMC = aod::McParticles;
  Preslice<TrackTableMC> perCollision = aod::track::collisionId;
  void processMC(CollisionTableMCTrue::iterator const& mcCollision, CollisionTableMC const& collisions, TrackTableMC const& tracks, ParticleTableMC const& particles);
  PROCESS_SWITCH(ueCharged, processMC, "process MC", false);

  using CollisionTableMCData = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels>;
  using TrackTableMCData = soa::Filtered<soa::Join<aod::Tracks, aod::McTrackLabels, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection>>;
  void processDataMC(CollisionTableMCData::iterator const& collision, TrackTableMCData const& tracks, ParticleTableMC const& particles, aod::McCollisions const& mcCollisions);
  PROCESS_SWITCH(ueCharged, processDataMC, "process data MC", false);

  using CollisionTableData = soa::Join<aod::Collisions, aod::EvSels>;
  using TrackTableData = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection>>;
  void processData(CollisionTableData::iterator const& collision, TrackTableData const& tracks);
  PROCESS_SWITCH(ueCharged, processData, "process data", false);

  // add new method
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{};
  workflow.push_back(adaptAnalysisTask<ueCharged>(cfgc));
  return workflow;
}
// implementation
float ueCharged::DeltaPhi(float phia, float phib,
                          float rangeMin = -M_PI / 2.0, float rangeMax = 3.0 * M_PI / 2.0)
{
  float dphi = -999;

  if (phia < 0) {
    phia += 2 * M_PI;
  } else if (phia > 2 * M_PI) {
    phia -= 2 * M_PI;
  }
  if (phib < 0) {
    phib += 2 * M_PI;
  } else if (phib > 2 * M_PI) {
    phib -= 2 * M_PI;
  }
  dphi = phib - phia;
  if (dphi < rangeMin) {
    dphi += 2 * M_PI;
  } else if (dphi > rangeMax) {
    dphi -= 2 * M_PI;
  }

  return dphi;
}

void ueCharged::init(InitContext const&)
{

  mySelectionPrim = myTrackSelectionPrim();
  mySelectionOpenDCA = myTrackSelectionOpenDCA();

  ConfigurableAxis ptBinningt{"ptBinningt", {0, 0.15, 0.50, 1.00, 1.50, 2.00, 2.50, 3.00, 3.50, 4.00, 4.50, 5.00, 6.00, 7.00, 8.00, 9.00, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 25.0, 30.0, 40.0, 50.0}, "pTtrig bin limits"};
  AxisSpec ptAxist = {ptBinningt, "#it{p}_{T}^{trig} (GeV/#it{c})"};

  ConfigurableAxis ptBinning{"ptBinning", {0, 0.0, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 25.0, 30.0, 40.0, 50.0}, "pTassoc bin limits"};
  AxisSpec ptAxis = {ptBinning, "#it{p}_{T}^{assoc} (GeV/#it{c})"};

  f_Eff.setObject(new TF1("fpara", "(x<0.3)*((0.283781)+x*(3.0492)+x*x*(-6.17018)) + (x>=0.3&&x<1.8)*((0.597121)+(0.200737)*x+(-0.11255)*x*x+(0.0242807)*x*x*x) + (x>=1.8&&x<14.)*((0.729892)+(0.0018516)*x+(0.000257896)*x*x+(-2.05202e-05)*x*x*x) + (x>=14)*(0.749313)", 0., 1e5));

  if (doprocessMC) {
    ue.add("hPtOut", "pT all rec; pT; Nch", HistType::kTH1D, {ptAxis});
    ue.add("hPtInPrim", "pT mc prim; pT; Nch", HistType::kTH1D, {ptAxis});
    ue.add("hPtInPrimGen", "pT mc prim all gen; pT; Nch", HistType::kTH1D, {ptAxis});
    ue.add("hPtOutPrim", "pT rec prim; pT; Nch", HistType::kTH1D, {ptAxis});
    ue.add("hPtOutSec", "pT rec sec; pT; Nch", HistType::kTH1D, {ptAxis});
    ue.add("hPtDCAall", "all MC; DCA_xy; Nch", HistType::kTH2D, {{ptAxis}, {121, -3.025, 3.025, "#it{DCA}_{xy} (cm)"}});
    ue.add("hPtDCAPrimary", "primary; DCA_xy; Nch", HistType::kTH2D, {{ptAxis}, {121, -3.025, 3.025, "#it{DCA}_{xy} (cm)"}});
    ue.add("hPtDCAWeak", "Weak decays; DCA_xy; Nch", HistType::kTH2D, {{ptAxis}, {121, -3.025, 3.025, "#it{DCA}_{xy} (cm)"}});
    ue.add("hPtDCAMat", "Material; DCA_xy; Nch", HistType::kTH2D, {{ptAxis}, {121, -3.025, 3.025, "#it{DCA}_{xy} (cm)"}});

    ue.add("hmultTrue", "mult true", HistType::kTH1F, {{200, -0.5, 199.5, " "}});
    ue.add("hmultTrueGen", "mult true all Gen", HistType::kTH1F, {{200, -0.5, 199.5, " "}});
    ue.add("hvtxZmc", "vtxZ mctrue", HistType::kTH1F, {{40, -20.0, 20.0, " "}});
    ue.add("hPtLeadingTrue", "true pTleading after physics selection", HistType::kTH1D, {ptAxist});

    for (int i = 0; i < 3; ++i) {
      ue.add(hPtVsPtLeadingTrue[i].data(), " ", HistType::kTH2D, {{ptAxist}, {ptAxis}});
      ue.add(pNumDenTrueAll[i].data(), "", HistType::kTProfile, {ptAxist});
      ue.add(pSumPtTrueAll[i].data(), "", HistType::kTProfile, {ptAxist});
      ue.add(pNumDenTrue[i].data(), "", HistType::kTProfile, {ptAxist});
      ue.add(pSumPtTrue[i].data(), "", HistType::kTProfile, {ptAxist});
      ue.add(pNumDenTruePS[i].data(), "", HistType::kTProfile, {ptAxist});
      ue.add(pSumPtTruePS[i].data(), "", HistType::kTProfile, {ptAxist});
    }
    for (int i = 0; i < 3; ++i) {
      ue.add(hPhiTrue[i].data(), "all charged true; #Delta#phi; Counts", HistType::kTH1D, {{64, -M_PI / 2.0, 3.0 * M_PI / 2.0, ""}});
    }
  }
  ue.add("hStat", "TotalEvents", HistType::kTH1F, {{1, 0.5, 1.5, " "}});
  ue.add("hmultRec", "mult rec", HistType::kTH1F, {{200, -0.5, 199.5, " "}});
  ue.add("hdNdeta", "dNdeta", HistType::kTH1F, {{50, -2.5, 2.5, " "}});
  ue.add("vtxZEta", ";#eta;vtxZ", HistType::kTH2F, {{50, -2.5, 2.5, " "}, {60, -30, 30, " "}});
  ue.add("phiEta", ";#eta;#varphi", HistType::kTH2F, {{50, -2.5, 2.5}, {200, 0., 2 * M_PI, " "}});
  ue.add("hvtxZ", "vtxZ", HistType::kTH1F, {{40, -20.0, 20.0, " "}});

  ue.add("hCounter", "Counter; sel; Nev", HistType::kTH1D, {{7, 0, 7, " "}});
  ue.add("hPtLeadingRecPS", "rec pTleading after physics selection", HistType::kTH1D, {ptAxist});
  ue.add("hPtLeadingMeasured", "measured pTleading after physics selection", HistType::kTH1D, {ptAxist});

  for (int i = 0; i < 3; ++i) {
    ue.add(pNumDenMeasuredPS[i].data(), "Number Density; ; #LT #it{N}_{trk} #GT", HistType::kTProfile, {ptAxist});
    ue.add(pSumPtMeasuredPS[i].data(), "Total #it{p}_{T}; ; #LT#sum#it{p}_{T}#GT", HistType::kTProfile, {ptAxist});
    ue.add(hPhi[i].data(), "all charged; #Delta#phi; Counts", HistType::kTH1D, {{64, -M_PI / 2.0, 3.0 * M_PI / 2.0, ""}});
  }

  // Data driven
  for (int i = 0; i < 3; ++i) {
    ue.add(hNumDenMCDd[i].data(), " ", HistType::kTH2D, {{ptAxist}, {100, -0.5, 99.5, "#it{N}_{trk}"}});
    ue.add(hSumPtMCDd[i].data(), " ", HistType::kTH2D, {{ptAxist}, {ptAxis}});
    ue.add(hNumDenMCMatchDd[i].data(), " ", HistType::kTH2D, {{ptAxist}, {100, -0.5, 99.5, "#it{N}_{trk}"}});
    ue.add(hSumPtMCMatchDd[i].data(), " ", HistType::kTH2D, {{ptAxist}, {ptAxis}});
  }

  for (int i = 0; i < 3; ++i) {
    ue.add(hPtVsPtLeadingData[i].data(), " ", HistType::kTH2D, {{ptAxist}, {ptAxis}});
    ue.add(pNumDenData[i].data(), "", HistType::kTProfile, {ptAxist});
    ue.add(pSumPtData[i].data(), "", HistType::kTProfile, {ptAxist});
  }

  ue.add("hPtLeadingData", " ", HistType::kTH1D, {{ptAxist}});
  ue.add("hPTVsDCAData", " ", HistType::kTH2D, {{ptAxis}, {121, -3.025, 3.025, "#it{DCA}_{xy} (cm)"}});
  ue.add("hEtaLeadingVsPtLeading", " ", HistType::kTH2D, {{ptAxist}, {50, -2.5, 2.5, "#eta"}});
}

void ueCharged::processMC(CollisionTableMCTrue::iterator const& mcCollision, CollisionTableMC const& collisions, TrackTableMC const& tracks, ParticleTableMC const& particles)
{
  if (collisions.size() != 0) {
    for (auto& collision : collisions) {
      auto curTracks = tracks.sliceBy(perCollision, collision.globalIndex());
      processMeasMC<true>(collision, curTracks, particles);
      break; // for now look only at first collision...
    }
  }
  processTrue(mcCollision, particles);
}

void ueCharged::processDataMC(CollisionTableMCData::iterator const& collision, TrackTableMCData const& tracks, ParticleTableMC const& particles, aod::McCollisions const& /*mcCollisions*/)
{
  processMeasMC<false>(collision, tracks, particles);
}
void ueCharged::processData(CollisionTableData::iterator const& collision, TrackTableData const& tracks)
{
  processMeas(collision, tracks);
}

template <typename C, typename P>
void ueCharged::processTrue(const C& mcCollision, const P& particles)
{
  int multTrue = 0;
  int multTrueINEL = 0;
  for (auto& particle : particles) {
    auto pdgParticle = pdg->GetParticle(particle.pdgCode());
    if (!pdgParticle || pdgParticle->Charge() == 0.) {
      continue;
    }
    if (!particle.isPhysicalPrimary()) {
      continue;
    }
    if (std::abs(particle.eta()) <= 1.0) {
      multTrueINEL++;
    }
    if (std::abs(particle.eta()) >= cfgTrkEtaCut) {
      continue;
    }
    multTrue++;
    if (particle.pt() < cfgTrkLowPtCut) {
      continue;
    }
    ue.fill(HIST("hPtInPrimGen"), particle.pt());
  }
  ue.fill(HIST("hmultTrueGen"), multTrue);
  if (std::abs(mcCollision.posZ()) > 10.f && multTrueINEL <= 0) {
    return;
  }

  ue.fill(HIST("hvtxZmc"), mcCollision.posZ());

  double flPtTrue = 0; // leading pT
  double flPhiTrue = 0;
  int flIndexTrue = 0;

  for (auto& particle : particles) {
    auto pdgParticle = pdg->GetParticle(particle.pdgCode());
    if (!pdgParticle || pdgParticle->Charge() == 0.) {
      continue;
    }
    if (!particle.isPhysicalPrimary()) {
      continue;
    }
    if (std::abs(particle.eta()) >= cfgTrkEtaCut) {
      continue;
    }

    if (particle.pt() < cfgTrkLowPtCut) {
      continue;
    }

    // ue.fill(HIST("hPtInPrim"), particle.pt());

    if (flPtTrue < particle.pt()) {
      flPtTrue = particle.pt();
      flPhiTrue = particle.phi();
      flIndexTrue = particle.globalIndex();
    }
  }
  ue.fill(HIST("hPtLeadingTrue"), flPtTrue);
  std::vector<double> ue_true;
  ue_true.clear();
  int nchm_toptrue[3];
  double sumptm_toptrue[3];
  for (int i = 0; i < 3; ++i) {
    nchm_toptrue[i] = 0;
    sumptm_toptrue[i] = 0;
  }
  std::vector<Float_t> ptArrayTrue;
  std::vector<Float_t> phiArrayTrue;
  std::vector<int> indexArrayTrue;

  for (auto& particle : particles) {

    auto pdgParticle = pdg->GetParticle(particle.pdgCode());
    if (!pdgParticle || pdgParticle->Charge() == 0.) {
      continue;
    }
    if (!particle.isPhysicalPrimary()) {
      continue;
    }
    if (std::abs(particle.eta()) >= cfgTrkEtaCut) {
      continue;
    }
    if (particle.pt() < cfgTrkLowPtCut) {
      continue;
    }
    // remove the autocorrelation
    if (flIndexTrue == particle.globalIndex()) {
      continue;
    }
    double DPhi = DeltaPhi(particle.phi(), flPhiTrue);

    // definition of the topological regions
    if (TMath::Abs(DPhi) < M_PI / 3.0) { // near side
      ue.fill(HIST(hPhiTrue[0]), DPhi);
      ue.fill(HIST(hPtVsPtLeadingTrue[0]), flPtTrue, particle.pt());
      nchm_toptrue[0]++;
      sumptm_toptrue[0] += particle.pt();
    } else if (TMath::Abs(DPhi - M_PI) < M_PI / 3.0) { // away side
      ue.fill(HIST(hPhiTrue[1]), DPhi);
      ue.fill(HIST(hPtVsPtLeadingTrue[1]), flPtTrue, particle.pt());
      nchm_toptrue[1]++;
      sumptm_toptrue[1] += particle.pt();
    } else { // transverse side
      ue.fill(HIST(hPhiTrue[2]), DPhi);
      ue.fill(HIST(hPtVsPtLeadingTrue[2]), flPtTrue, particle.pt());
      nchm_toptrue[2]++;
      sumptm_toptrue[2] += particle.pt();
    }
  }

  for (int i_reg = 0; i_reg < 3; ++i_reg) {
    ue_true.push_back(1.0 * nchm_toptrue[i_reg]);
  }
  for (int i_reg = 0; i_reg < 3; ++i_reg) {
    ue_true.push_back(sumptm_toptrue[i_reg]);
  }
  ue.fill(HIST(pNumDenTrueAll[0]), flPtTrue, ue_true[0]);
  ue.fill(HIST(pSumPtTrueAll[0]), flPtTrue, ue_true[3]);

  ue.fill(HIST(pNumDenTrueAll[1]), flPtTrue, ue_true[1]);
  ue.fill(HIST(pSumPtTrueAll[1]), flPtTrue, ue_true[4]);

  ue.fill(HIST(pNumDenTrueAll[2]), flPtTrue, ue_true[2]);
  ue.fill(HIST(pSumPtTrueAll[2]), flPtTrue, ue_true[5]);

  ptArrayTrue.clear();
  phiArrayTrue.clear();
  indexArrayTrue.clear();
}
template <typename C, typename T>
void ueCharged::processMeas(const C& collision, const T& tracks)
{

  ue.fill(HIST("hCounter"), 0);

  if (sel8 && !collision.sel8()) {
    return;
  }

  if (removeITSROFBorder && (!collision.selection_bit(o2::aod::evsel::kIsTriggerTVX) || !collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder))) {
    return;
  }

  if (manuallyApplysel8 && (!collision.selection_bit(o2::aod::evsel::kIsTriggerTVX) || !collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder) || !collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder))) {
    return;
  }

  ue.fill(HIST("hCounter"), 1);

  if (piluprejection && !collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
    return;
  }

  ue.fill(HIST("hCounter"), 2);
  if (goodzvertex && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
    return;
  }

  ue.fill(HIST("hCounter"), 3);

  ue.fill(HIST("hStat"), collision.size());
  auto vtxZ = collision.posZ();
  // + vtx
  if ((std::abs(collision.posZ()) >= 10.f)) {
    return;
  }

  ue.fill(HIST("hCounter"), 4);

  ue.fill(HIST("hvtxZ"), vtxZ);

  // loop over selected tracks
  double flPt = 0; // leading pT
  double flPhi = 0;
  double flEta = 0;
  int flIndex = 0;
  int multRec = 0;
  for (auto& track : tracks) {
    if (!mySelectionPrim.IsSelected(track)) {
      continue;
    }

    ue.fill(HIST("hdNdeta"), track.eta());
    ue.fill(HIST("vtxZEta"), track.eta(), vtxZ);
    ue.fill(HIST("phiEta"), track.eta(), track.phi());

    if (flPt < track.pt()) {
      multRec++;
      flPt = track.pt();
      flPhi = track.phi();
      flIndex = track.globalIndex();
      flEta = track.eta();
    }
  }
  ue.fill(HIST("hmultRec"), multRec);
  ue.fill(HIST("hPtLeadingMeasured"), flPt);
  ue.fill(HIST("hPtLeadingRecPS"), flPt);
  ue.fill(HIST("hEtaLeadingVsPtLeading"), flPt, flEta);

  std::vector<double> ue_rec;
  ue_rec.clear();
  int nchm_top[3];
  double sumptm_top[3];
  for (int i = 0; i < 3; ++i) {
    nchm_top[i] = 0;
    sumptm_top[i] = 0;
  }
  std::vector<Float_t> ptArray;
  std::vector<Float_t> phiArray;
  std::vector<int> indexArray;

  for (auto& track : tracks) {

    if (mySelectionOpenDCA.IsSelected(track)) {
      ue.fill(HIST("hPTVsDCAData"), track.pt(), track.dcaXY());
    }

    if (mySelectionPrim.IsSelected(track)) {
      // applying the efficiency twice for the misrec of leading particle
      if (f_Eff->Eval(track.pt()) > gRandom->Uniform(0, 1)) {
        ptArray.push_back(track.pt());
        phiArray.push_back(track.phi());
        indexArray.push_back(track.globalIndex());
      }

      // remove the autocorrelation
      if (flIndex == track.globalIndex()) {
        continue;
      }

      double DPhi = DeltaPhi(track.phi(), flPhi);

      // definition of the topological regions
      if (TMath::Abs(DPhi) < M_PI / 3.0) { // near side
        ue.fill(HIST(hPhi[0]), DPhi);
        ue.fill(HIST(hPtVsPtLeadingData[0]), flPt, track.pt());
        nchm_top[0]++;
        sumptm_top[0] += track.pt();
      } else if (TMath::Abs(DPhi - M_PI) < M_PI / 3.0) { // away side
        ue.fill(HIST(hPhi[1]), DPhi);
        ue.fill(HIST(hPtVsPtLeadingData[1]), flPt, track.pt());
        nchm_top[1]++;
        sumptm_top[1] += track.pt();
      } else { // transverse side
        ue.fill(HIST(hPhi[2]), DPhi);
        ue.fill(HIST(hPtVsPtLeadingData[2]), flPt, track.pt());
        nchm_top[2]++;
        sumptm_top[2] += track.pt();
      }
    }
  }
  for (int i_reg = 0; i_reg < 3; ++i_reg) {
    ue_rec.push_back(1.0 * nchm_top[i_reg]);
  }
  for (int i_reg = 0; i_reg < 3; ++i_reg) {
    ue_rec.push_back(sumptm_top[i_reg]);
  }

  // add flags for Vtx, PS, ev sel
  ue.fill(HIST(pNumDenMeasuredPS[0]), flPt, ue_rec[0]);
  ue.fill(HIST(pNumDenData[0]), flPt, ue_rec[0]);
  ue.fill(HIST(pSumPtMeasuredPS[0]), flPt, ue_rec[3]);
  ue.fill(HIST(pSumPtData[0]), flPt, ue_rec[3]);

  ue.fill(HIST(pNumDenMeasuredPS[1]), flPt, ue_rec[1]);
  ue.fill(HIST(pNumDenData[1]), flPt, ue_rec[1]);
  ue.fill(HIST(pSumPtMeasuredPS[1]), flPt, ue_rec[4]);
  ue.fill(HIST(pSumPtData[1]), flPt, ue_rec[4]);

  ue.fill(HIST(pNumDenMeasuredPS[2]), flPt, ue_rec[2]);
  ue.fill(HIST(pNumDenData[2]), flPt, ue_rec[2]);
  ue.fill(HIST(pSumPtMeasuredPS[2]), flPt, ue_rec[5]);
  ue.fill(HIST(pSumPtData[2]), flPt, ue_rec[5]);

  ue.fill(HIST("hPtLeadingData"), flPt);

  // Compute data driven (DD) missidentification correction
  Float_t flPtdd = 0; // leading pT
  Float_t flPhidd = 0;
  int flIndexdd = 0;
  int ntrkdd = ptArray.size();

  for (int i = 0; i < ntrkdd; ++i) {
    if (flPtdd < ptArray[i]) {
      flPtdd = ptArray[i];
      flPhidd = phiArray[i];
      flIndexdd = indexArray[i];
    }
  }
  int nchm_topdd[3];
  double sumptm_topdd[3];
  for (int i = 0; i < 3; ++i) {
    nchm_topdd[i] = 0;
    sumptm_topdd[i] = 0;
  }
  for (int i = 0; i < ntrkdd; ++i) {
    if (indexArray[i] == flIndexdd) {
      continue;
    }
    double DPhi = DeltaPhi(phiArray[i], flPhidd);
    if (TMath::Abs(DPhi) < M_PI / 3.0) { // near side
      nchm_topdd[0]++;
      sumptm_topdd[0] += ptArray[i];
    } else if (TMath::Abs(DPhi - M_PI) < M_PI / 3.0) { // away side
      nchm_topdd[1]++;
      sumptm_topdd[1] += ptArray[i];
    } else { // transverse side
      nchm_topdd[2]++;
      sumptm_topdd[2] += ptArray[i];
    }
  }

  ue.fill(HIST(hNumDenMCDd[0]), flPtdd, nchm_topdd[0]);
  ue.fill(HIST(hSumPtMCDd[0]), flPtdd, sumptm_topdd[0]);

  ue.fill(HIST(hNumDenMCDd[1]), flPtdd, nchm_topdd[1]);
  ue.fill(HIST(hSumPtMCDd[1]), flPtdd, sumptm_topdd[1]);

  ue.fill(HIST(hNumDenMCDd[2]), flPtdd, nchm_topdd[2]);
  ue.fill(HIST(hSumPtMCDd[2]), flPtdd, sumptm_topdd[2]);

  if (flIndexdd == flIndex) {
    ue.fill(HIST(hNumDenMCMatchDd[0]), flPtdd, nchm_topdd[0]);
    ue.fill(HIST(hSumPtMCMatchDd[0]), flPtdd, sumptm_topdd[0]);

    ue.fill(HIST(hNumDenMCMatchDd[1]), flPtdd, nchm_topdd[1]);
    ue.fill(HIST(hSumPtMCMatchDd[1]), flPtdd, sumptm_topdd[1]);

    ue.fill(HIST(hNumDenMCMatchDd[2]), flPtdd, nchm_topdd[2]);
    ue.fill(HIST(hSumPtMCMatchDd[2]), flPtdd, sumptm_topdd[2]);
  }
  ptArray.clear();
  phiArray.clear();
  indexArray.clear();
}

template <bool IS_MC, typename C, typename T, typename P>
void ueCharged::processMeasMC(const C& collision, const T& tracks, const P& particles)
{

  ue.fill(HIST("hCounter"), 0);

  ue.fill(HIST("hStat"), collision.size());
  auto vtxZ = collision.posZ();

  // int multTrue = 0;
  double flPtTrue = 0; // leading pT
  double flPhiTrue = 0;
  int flIndexTrue = 0;

  for (auto& particle : particles) {
    auto pdgParticle = pdg->GetParticle(particle.pdgCode());
    if (!pdgParticle || pdgParticle->Charge() == 0.) {
      continue;
    }
    if (!particle.isPhysicalPrimary()) {
      continue;
    }
    if (std::abs(particle.eta()) >= cfgTrkEtaCut) {
      continue;
    }
    // multTrue++;
    if (particle.pt() < cfgTrkLowPtCut) {
      continue;
    }

    if (flPtTrue < particle.pt()) {
      flPtTrue = particle.pt();
      flPhiTrue = particle.phi();
      flIndexTrue = particle.globalIndex();
    }
  }
  // ue.fill(HIST("hmultTrue"), multTrue);
  std::vector<double> ue_true;
  ue_true.clear();
  int nchm_toptrue[3];
  double sumptm_toptrue[3];
  for (int i = 0; i < 3; ++i) {
    nchm_toptrue[i] = 0;
    sumptm_toptrue[i] = 0;
  }
  std::vector<Float_t> ptArrayTrue;
  std::vector<Float_t> phiArrayTrue;
  std::vector<int> indexArrayTrue;

  for (auto& particle : particles) {

    auto pdgParticle = pdg->GetParticle(particle.pdgCode());
    if (!pdgParticle || pdgParticle->Charge() == 0.) {
      continue;
    }
    if (!particle.isPhysicalPrimary()) {
      continue;
    }
    if (std::abs(particle.eta()) >= cfgTrkEtaCut) {
      continue;
    }
    if (particle.pt() < cfgTrkLowPtCut) {
      continue;
    }
    // remove the autocorrelation
    if (flIndexTrue == particle.globalIndex()) {
      continue;
    }
    double DPhi = DeltaPhi(particle.phi(), flPhiTrue);

    // definition of the topological regions
    if (TMath::Abs(DPhi) < M_PI / 3.0) { // near side
      nchm_toptrue[0]++;
      sumptm_toptrue[0] += particle.pt();
    } else if (TMath::Abs(DPhi - M_PI) < M_PI / 3.0) { // away side
      nchm_toptrue[1]++;
      sumptm_toptrue[1] += particle.pt();
    } else { // transverse side
      nchm_toptrue[2]++;
      sumptm_toptrue[2] += particle.pt();
    }
  }

  for (int i_reg = 0; i_reg < 3; ++i_reg) {
    ue_true.push_back(1.0 * nchm_toptrue[i_reg]);
  }
  for (int i_reg = 0; i_reg < 3; ++i_reg) {
    ue_true.push_back(sumptm_toptrue[i_reg]);
  }

  ue.fill(HIST(pNumDenTrue[0]), flPtTrue, ue_true[0]);
  ue.fill(HIST(pSumPtTrue[0]), flPtTrue, ue_true[3]);

  ue.fill(HIST(pNumDenTrue[1]), flPtTrue, ue_true[1]);
  ue.fill(HIST(pSumPtTrue[1]), flPtTrue, ue_true[4]);

  ue.fill(HIST(pNumDenTrue[2]), flPtTrue, ue_true[2]);
  ue.fill(HIST(pSumPtTrue[2]), flPtTrue, ue_true[5]);

  ptArrayTrue.clear();
  phiArrayTrue.clear();
  indexArrayTrue.clear();

  if (sel8 && !collision.sel8()) {
    return;
  }

  if (removeITSROFBorder && (!collision.selection_bit(o2::aod::evsel::kIsTriggerTVX) || !collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder))) {
    return;
  }

  if (manuallyApplysel8 && (!collision.selection_bit(o2::aod::evsel::kIsTriggerTVX) || !collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder) || !collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder))) {
    return;
  }

  ue.fill(HIST("hCounter"), 1);

  if (piluprejection && !collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
    return;
  }
  ue.fill(HIST("hCounter"), 2);

  if (goodzvertex && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
    return;
  }
  ue.fill(HIST("hCounter"), 3);

  // only PS
  if ((std::abs(collision.posZ()) >= 10.f)) {
    return;
  }
  ue.fill(HIST("hCounter"), 4);

  ue.fill(HIST(pNumDenTruePS[0]), flPtTrue, ue_true[0]);
  ue.fill(HIST(pSumPtTruePS[0]), flPtTrue, ue_true[3]);

  ue.fill(HIST(pNumDenTruePS[1]), flPtTrue, ue_true[1]);
  ue.fill(HIST(pSumPtTruePS[1]), flPtTrue, ue_true[4]);

  ue.fill(HIST(pNumDenTruePS[2]), flPtTrue, ue_true[2]);
  ue.fill(HIST(pSumPtTruePS[2]), flPtTrue, ue_true[5]);

  // ue.fill(HIST("hCounter"), 2);

  ue.fill(HIST("hvtxZ"), vtxZ);
  // loop over MC true particles
  int multTrue = 0;
  for (auto& particle : particles) {
    auto pdgParticle = pdg->GetParticle(particle.pdgCode());
    if (!pdgParticle || pdgParticle->Charge() == 0.) {
      continue;
    }
    if (!particle.isPhysicalPrimary()) {
      continue;
    }
    if (std::abs(particle.eta()) >= cfgTrkEtaCut) {
      continue;
    }
    multTrue++;
    if (particle.pt() < cfgTrkLowPtCut) {
      continue;
    }
    ue.fill(HIST("hPtInPrim"), particle.pt());
  }
  ue.fill(HIST("hmultTrue"), multTrue);

  // loop over selected tracks
  double flPt = 0; // leading pT
  double flPhi = 0;
  int flIndex = 0;
  int multRec = 0;
  for (auto& track : tracks) {
    if (!mySelectionPrim.IsSelected(track)) {
      continue;
    }

    ue.fill(HIST("hdNdeta"), track.eta());
    ue.fill(HIST("vtxZEta"), track.eta(), vtxZ);
    ue.fill(HIST("phiEta"), track.eta(), track.phi());

    if (flPt < track.pt()) {
      multRec++;
      flPt = track.pt();
      flPhi = track.phi();
      flIndex = track.globalIndex();
    }
  }
  ue.fill(HIST("hmultRec"), multRec);
  ue.fill(HIST("hPtLeadingMeasured"), flPt);
  ue.fill(HIST("hPtLeadingRecPS"), flPt);
  std::vector<double> ue_rec;
  ue_rec.clear();
  int nchm_top[3];
  double sumptm_top[3];
  for (int i = 0; i < 3; ++i) {
    nchm_top[i] = 0;
    sumptm_top[i] = 0;
  }
  std::vector<Float_t> ptArray;
  std::vector<Float_t> phiArray;
  std::vector<int> indexArray;

  for (auto& track : tracks) {

    if (mySelectionOpenDCA.IsSelected(track)) { // TODO: set cuts w/o DCA cut
      ue.fill(HIST("hPTVsDCAData"), track.pt(), track.dcaXY());
    }
    if constexpr (IS_MC) {
      if (track.has_mcParticle()) {
        if (mySelectionPrim.IsSelected(track)) {
          ue.fill(HIST("hPtOut"), track.pt());
        }
        if (mySelectionOpenDCA.IsSelected(track)) {
          ue.fill(HIST("hPtDCAall"), track.pt(), track.dcaXY());
        }
        const auto& particle = track.template mcParticle_as<aod::McParticles>();
        if (particle.isPhysicalPrimary()) {

          if (mySelectionPrim.IsSelected(track)) {
            // ue.fill(HIST("hPtOutPrim"), track.pt());
            ue.fill(HIST("hPtOutPrim"), particle.pt());
          }
          if (mySelectionOpenDCA.IsSelected(track)) {
            ue.fill(HIST("hPtDCAPrimary"), track.pt(), track.dcaXY());
          }
          // LOGP(info, "this track has MC particle {}", 1);
        } else {
          if (mySelectionPrim.IsSelected(track)) {
            ue.fill(HIST("hPtOutSec"), track.pt());
          }
          if (mySelectionOpenDCA.IsSelected(track)) {
            if (particle.producedByGenerator()) { // i guess these are from decays
              ue.fill(HIST("hPtDCAWeak"), track.pt(), track.dcaXY());
            } else { //// i guess these are from material
              ue.fill(HIST("hPtDCAMat"), track.pt(), track.dcaXY());
            }
          }
        }
      }
    }

    if (mySelectionPrim.IsSelected(track)) {
      // applying the efficiency twice for the misrec of leading particle
      if (f_Eff->Eval(track.pt()) > gRandom->Uniform(0, 1)) {
        ptArray.push_back(track.pt());
        phiArray.push_back(track.phi());
        indexArray.push_back(track.globalIndex());
      }

      // remove the autocorrelation
      if (flIndex == track.globalIndex()) {
        continue;
      }

      double DPhi = DeltaPhi(track.phi(), flPhi);

      // definition of the topological regions
      if (TMath::Abs(DPhi) < M_PI / 3.0) { // near side
        ue.fill(HIST(hPhi[0]), DPhi);
        ue.fill(HIST(hPtVsPtLeadingData[0]), flPt, track.pt());
        nchm_top[0]++;
        sumptm_top[0] += track.pt();
      } else if (TMath::Abs(DPhi - M_PI) < M_PI / 3.0) { // away side
        ue.fill(HIST(hPhi[1]), DPhi);
        ue.fill(HIST(hPtVsPtLeadingData[1]), flPt, track.pt());
        nchm_top[1]++;
        sumptm_top[1] += track.pt();
      } else { // transverse side
        ue.fill(HIST(hPhi[2]), DPhi);
        ue.fill(HIST(hPtVsPtLeadingData[2]), flPt, track.pt());
        nchm_top[2]++;
        sumptm_top[2] += track.pt();
      }
    }
  }
  for (int i_reg = 0; i_reg < 3; ++i_reg) {
    ue_rec.push_back(1.0 * nchm_top[i_reg]);
  }
  for (int i_reg = 0; i_reg < 3; ++i_reg) {
    ue_rec.push_back(sumptm_top[i_reg]);
  }

  // add flags for Vtx, PS, ev sel
  ue.fill(HIST(pNumDenMeasuredPS[0]), flPt, ue_rec[0]);
  ue.fill(HIST(pNumDenData[0]), flPt, ue_rec[0]);
  ue.fill(HIST(pSumPtMeasuredPS[0]), flPt, ue_rec[3]);
  ue.fill(HIST(pSumPtData[0]), flPt, ue_rec[3]);

  ue.fill(HIST(pNumDenMeasuredPS[1]), flPt, ue_rec[1]);
  ue.fill(HIST(pNumDenData[1]), flPt, ue_rec[1]);
  ue.fill(HIST(pSumPtMeasuredPS[1]), flPt, ue_rec[4]);
  ue.fill(HIST(pSumPtData[1]), flPt, ue_rec[4]);

  ue.fill(HIST(pNumDenMeasuredPS[2]), flPt, ue_rec[2]);
  ue.fill(HIST(pNumDenData[2]), flPt, ue_rec[2]);
  ue.fill(HIST(pSumPtMeasuredPS[2]), flPt, ue_rec[5]);
  ue.fill(HIST(pSumPtData[2]), flPt, ue_rec[5]);

  ue.fill(HIST("hPtLeadingData"), flPt);

  // Compute data driven (DD) missidentification correction
  Float_t flPtdd = 0; // leading pT
  Float_t flPhidd = 0;
  int flIndexdd = 0;
  int ntrkdd = ptArray.size();

  for (int i = 0; i < ntrkdd; ++i) {
    if (flPtdd < ptArray[i]) {
      flPtdd = ptArray[i];
      flPhidd = phiArray[i];
      flIndexdd = indexArray[i];
    }
  }
  int nchm_topdd[3];
  double sumptm_topdd[3];
  for (int i = 0; i < 3; ++i) {
    nchm_topdd[i] = 0;
    sumptm_topdd[i] = 0;
  }
  for (int i = 0; i < ntrkdd; ++i) {
    if (indexArray[i] == flIndexdd) {
      continue;
    }
    double DPhi = DeltaPhi(phiArray[i], flPhidd);
    if (TMath::Abs(DPhi) < M_PI / 3.0) { // near side
      nchm_topdd[0]++;
      sumptm_topdd[0] += ptArray[i];
    } else if (TMath::Abs(DPhi - M_PI) < M_PI / 3.0) { // away side
      nchm_topdd[1]++;
      sumptm_topdd[1] += ptArray[i];
    } else { // transverse side
      nchm_topdd[2]++;
      sumptm_topdd[2] += ptArray[i];
    }
  }

  ue.fill(HIST(hNumDenMCDd[0]), flPtdd, nchm_topdd[0]);
  ue.fill(HIST(hSumPtMCDd[0]), flPtdd, sumptm_topdd[0]);

  ue.fill(HIST(hNumDenMCDd[1]), flPtdd, nchm_topdd[1]);
  ue.fill(HIST(hSumPtMCDd[1]), flPtdd, sumptm_topdd[1]);

  ue.fill(HIST(hNumDenMCDd[2]), flPtdd, nchm_topdd[2]);
  ue.fill(HIST(hSumPtMCDd[2]), flPtdd, sumptm_topdd[2]);

  if (flIndexdd == flIndex) {
    ue.fill(HIST(hNumDenMCMatchDd[0]), flPtdd, nchm_topdd[0]);
    ue.fill(HIST(hSumPtMCMatchDd[0]), flPtdd, sumptm_topdd[0]);

    ue.fill(HIST(hNumDenMCMatchDd[1]), flPtdd, nchm_topdd[1]);
    ue.fill(HIST(hSumPtMCMatchDd[1]), flPtdd, sumptm_topdd[1]);

    ue.fill(HIST(hNumDenMCMatchDd[2]), flPtdd, nchm_topdd[2]);
    ue.fill(HIST(hSumPtMCMatchDd[2]), flPtdd, sumptm_topdd[2]);
  }
  ptArray.clear();
  phiArray.clear();
  indexArray.clear();
}
