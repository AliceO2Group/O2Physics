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
//
/// \author Antonio Ortiz (antonio.ortiz@nucleares.unam.mx)
/// \file uecharged.cxx
/// \brief Underlying event analysis task
/// \since November 2021
/// \last update: September 2025

#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/StaticFor.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

#include "TDatabasePDG.h"
#include <TF1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TRandom.h>

#include <cmath>
#include <vector>

// TODO: implement 50% stat for MC closure vs 50% for testing, add flag for weak
// decays

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using BCsRun3 = soa::Join<aod::BCs, aod::Timestamps, aod::BcSels,
                          aod::Run3MatchedToBCSparse>;

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
    selectedTracks.SetMaxDcaXYPtDep([](float pt) { return 0.0105f + 0.0350f / std::pow(pt, 1.1f); });
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
    // selectedTracks.SetMaxDcaXYPtDep([](float pt) { return 0.0105f + 0.0350f /
    // pow(pt, 1.1f); });
    selectedTracks.SetMaxDcaZ(2.f);
    return selectedTracks;
  }

  TrackSelection mySelectionPrim;
  TrackSelection mySelectionOpenDCA;

  Service<o2::framework::O2DatabasePDG> pdg;
  float deltaPhi(float phia, float phib, float rangeMin, float rangeMax);
  // Configurable for event selection
  Configurable<bool> isRun3{"isRun3", true, "is Run3 dataset"};
  Configurable<bool> piluprejection{"piluprejection", true, "Pileup rejection"};
  Configurable<bool> goodzvertex{"goodzvertex", true, "removes collisions with large differences between z of PV by tracks and z of PV from FT0 A-C time difference"};
  Configurable<bool> sel8{"sel8", true, "Apply the sel8 event selection"};
  Configurable<bool> removeITSROFBorder{"removeITSROFBorder", false, "Remove ITS Read-Out Frame border and only apply kIsTriggerTVX & kNoTimeFrameBorder (recommended for MC)"};
  Configurable<bool> analyzeEvandTracksel{"analyzeEvandTracksel", true, "Analyze the event and track selection"};

  // acceptance cuts
  Configurable<float> cfgTrkEtaCut{"cfgTrkEtaCut", 0.8f, "Eta range for tracks"};
  Configurable<float> cfgTrkLowPtCut{"cfgTrkLowPtCut", 0.15f, "Minimum constituent pT"};

  HistogramRegistry ue;
  static constexpr std::string_view pNumDenMeasuredPS[3] = {
    "pNumDenMeasuredPS_NS", "pNumDenMeasuredPS_AS", "pNumDenMeasuredPS_TS"};
  static constexpr std::string_view pSumPtMeasuredPS[3] = {
    "pSumPtMeasuredPS_NS", "pSumPtMeasuredPS_AS", "pSumPtMeasuredPS_TS"};
  static constexpr std::string_view hPhi[3] = {"hPhi_NS", "hPhi_AS", "hPhi_TS"};
  // data driven correction
  static constexpr std::string_view hNumDenMCDd[3] = {
    "hNumDenMCDd_NS", "hNumDenMCDd_AS", "hNumDenMCDd_TS"};
  static constexpr std::string_view hSumPtMCDd[3] = {
    "hSumPtMCDd_NS", "hSumPtMCDd_AS", "hSumPtMCDd_TS"};
  static constexpr std::string_view hNumDenMCMatchDd[3] = {
    "hNumDenMCMatchDd_NS", "hNumDenMCMatchDd_AS", "hNumDenMCMatchDd_TS"};
  static constexpr std::string_view hSumPtMCMatchDd[3] = {
    "hSumPtMCMatchDd_NS", "hSumPtMCMatchDd_AS", "hSumPtMCMatchDd_TS"};
  // hist data for corrections
  static constexpr std::string_view hPtVsPtLeadingData[3] = {
    "hPtVsPtLeadingData_NS", "hPtVsPtLeadingData_AS",
    "hPtVsPtLeadingData_TS"};
  static constexpr std::string_view pNumDenData[3] = {
    "pNumDenData_NS", "pNumDenData_AS", "pNumDenData_TS"};
  static constexpr std::string_view pSumPtData[3] = {
    "pSumPtData_NS", "pSumPtData_AS", "pSumPtData_TS"};
  // hist data true
  static constexpr std::string_view hPtVsPtLeadingTrue[3] = {
    "hPtVsPtLeadingTrue_NS", "hPtVsPtLeadingTrue_AS",
    "hPtVsPtLeadingTrue_TS"};
  // all wo detector effects
  static constexpr std::string_view pNumDenTrueAll[3] = {
    "pNumDenTrueAll_NS", "pNumDenTrueAll_AS", "pNumDenTrueAll_TS"};
  static constexpr std::string_view pSumPtTrueAll[3] = {
    "pSumPtTrueAll_NS", "pSumPtTrueAll_AS", "pSumPtTrueAll_TS"};
  // true, 50%
  static constexpr std::string_view pNumDenTrue[3] = {
    "pNumDenTrue_NS", "pNumDenTrue_AS", "pNumDenTrue_TS"};
  static constexpr std::string_view pSumPtTrue[3] = {
    "pSumPtTrue_NS", "pSumPtTrue_AS", "pSumPtTrue_TS"};

  // this must have all event selection effects, but it has not been implemented
  // 50%
  static constexpr std::string_view pNumDenTruePS[3] = {
    "pNumDenTruePS_NS", "pNumDenTruePS_AS", "pNumDenTruePS_TS"};
  static constexpr std::string_view pSumPtTruePS[3] = {
    "pSumPtTruePS_NS", "pSumPtTruePS_AS", "pSumPtTruePS_TS"};
  static constexpr std::string_view hPhiTrue[3] = {"hPhiTrue_NS", "hPhiTrue_AS",
                                                   "hPhiTrue_TS"};

  OutputObj<TF1> fEff{"fpara"};
  void init(InitContext const&);

  template <bool IS_MC, typename C, typename T, typename P>
  void processMeasMC(const C& collision, const T& tracks, const P& particles);

  template <typename C, typename T>
  void processMeas(const C& collision, const T& tracks);

  template <typename C, typename P>
  void processTrue(const C& mcCollision, const P& particles);

  template <typename C, typename T>
  void analyzeEventAndTrackSelection(const C& collision, const T& tracks);

  Filter trackFilter = (nabs(aod::track::eta) < cfgTrkEtaCut) &&
                       (aod::track::pt > cfgTrkLowPtCut);

  using CollisionTableMCTrue = aod::McCollisions;
  using CollisionTableMC =
    soa::SmallGroups<soa::Join<aod::McCollisionLabels, aod::Collisions,
                               aod::EvSels, aod::PVMults>>;
  using TrackTableMC =
    soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA,
                            aod::McTrackLabels, aod::TrackSelection>>;
  using ParticleTableMC = aod::McParticles;
  Preslice<TrackTableMC> perCollision = aod::track::collisionId;
  void processMC(CollisionTableMCTrue::iterator const& mcCollision,
                 CollisionTableMC const& collisions, TrackTableMC const& tracks,
                 ParticleTableMC const& particles, aod::FT0s const&,
                 BCsRun3 const&);
  PROCESS_SWITCH(ueCharged, processMC, "process MC", false);

  using CollisionTableMCData =
    soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels,
              aod::PVMults>;
  using TrackTableMCData =
    soa::Filtered<soa::Join<aod::Tracks, aod::McTrackLabels, aod::TracksExtra,
                            aod::TracksDCA, aod::TrackSelection>>;
  void processDataMC(CollisionTableMCData::iterator const& collision,
                     TrackTableMCData const& tracks,
                     ParticleTableMC const& particles,
                     aod::McCollisions const& mcCollisions);
  PROCESS_SWITCH(ueCharged, processDataMC, "process data MC", false);

  using CollisionTableData =
    soa::Join<aod::Collisions, aod::EvSels, aod::PVMults>;
  using TrackTableData =
    soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA,
                            aod::TrackSelection>>;
  void processData(CollisionTableData::iterator const& collision,
                   TrackTableData const& tracks, aod::FT0s const&,
                   BCsRun3 const&);
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
float ueCharged::deltaPhi(float phia, float phib,
                          float rangeMin = -o2::constants::math::PI / 2.0,
                          float rangeMax = 3.0 * o2::constants::math::PI /
                                           2.0)
{
  float dphi = -999;

  if (phia < 0) {
    phia += 2 * o2::constants::math::PI;
  } else if (phia > 2 * o2::constants::math::PI) {
    phia -= 2 * o2::constants::math::PI;
  }
  if (phib < 0) {
    phib += 2 * o2::constants::math::PI;
  } else if (phib > 2 * o2::constants::math::PI) {
    phib -= 2 * o2::constants::math::PI;
  }
  dphi = phib - phia;
  if (dphi < rangeMin) {
    dphi += 2 * o2::constants::math::PI;
  } else if (dphi > rangeMax) {
    dphi -= 2 * o2::constants::math::PI;
  }

  return dphi;
}

void ueCharged::init(InitContext const&)
{

  mySelectionPrim = myTrackSelectionPrim();
  mySelectionOpenDCA = myTrackSelectionOpenDCA();

  ConfigurableAxis ptBinningt{"ptBinningt",
                              {0, 0.15, 0.50, 1.00, 1.50, 2.00, 2.50,
                               3.00, 3.50, 4.00, 4.50, 5.00, 6.00, 7.00,
                               8.00, 9.00, 10.0, 12.0, 14.0, 16.0, 18.0,
                               20.0, 25.0, 30.0, 40.0, 50.0},
                              "pTtrig bin limits"};
  AxisSpec ptAxist = {ptBinningt, "#it{p}_{T}^{trig} (GeV/#it{c})"};

  ConfigurableAxis ptBinning{
    "ptBinning",
    {0, 0.0, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45,
     0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4,
     1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8,
     3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0,
     12.0, 14.0, 16.0, 18.0, 20.0, 25.0, 30.0, 40.0, 50.0},
    "pTassoc bin limits"};
  AxisSpec ptAxis = {ptBinning, "#it{p}_{T}^{assoc} (GeV/#it{c})"};

  fEff.setObject(new TF1("fpara",
                         "(x<0.3)*((0.315318)+x*(2.38596)+x*x*(-4.388)) +"
                         "(x>=0.3&&x<1.8)*((0.604051)+(0.154763)*x+(-0.103004)*"
                         "x*x+(0.0266487)*x*x*x) +"
                         "(x>=1.8&&x<14.)*((0.700444)+(-0.00115506)*x+(0."
                         "000667608)*x*x+(-3.82915e-05)*x*x*x) +"
                         "(x>=14)*((0.731778)+(-0.000994634)*x)",
                         0., 1e5));

  if (doprocessMC) {
    ue.add("hPtOut", "pT all rec; pT; Nch", HistType::kTH1D, {ptAxis});
    ue.add("hPtInPrim", "pT mc prim; pT; Nch", HistType::kTH1D, {ptAxis});
    ue.add("hPtInPrimGen", "pT mc prim all gen; pT; Nch", HistType::kTH1D,
           {ptAxis});
    ue.add("hPtOutPrim", "pT rec prim; pT; Nch", HistType::kTH1D, {ptAxis});
    ue.add("hPtOutSec", "pT rec sec; pT; Nch", HistType::kTH1D, {ptAxis});
    ue.add("hPtDCAall", "all MC; DCA_xy; Nch", HistType::kTH2D,
           {{ptAxis}, {121, -3.025, 3.025, "#it{DCA}_{xy} (cm)"}});
    ue.add("hPtDCAPrimary", "primary; DCA_xy; Nch", HistType::kTH2D,
           {{ptAxis}, {121, -3.025, 3.025, "#it{DCA}_{xy} (cm)"}});
    ue.add("hPtDCAWeak", "Weak decays; DCA_xy; Nch", HistType::kTH2D,
           {{ptAxis}, {121, -3.025, 3.025, "#it{DCA}_{xy} (cm)"}});
    ue.add("hPtDCAMat", "Material; DCA_xy; Nch", HistType::kTH2D,
           {{ptAxis}, {121, -3.025, 3.025, "#it{DCA}_{xy} (cm)"}});

    ue.add("hmultTrue", "mult true", HistType::kTH1F,
           {{200, -0.5, 199.5, " "}});
    ue.add("hmultTrueGen", "mult true all Gen", HistType::kTH1F,
           {{200, -0.5, 199.5, " "}});
    ue.add("hvtxZmc", "vtxZ mctrue", HistType::kTH1F, {{40, -20.0, 20.0, " "}});
    ue.add("hPtLeadingTrue", "true pTleading after physics selection",
           HistType::kTH1D, {ptAxist});

    for (int i = 0; i < 3; ++i) {
      ue.add(hPtVsPtLeadingTrue[i].data(), " ", HistType::kTH2D,
             {{ptAxist}, {ptAxis}});
      ue.add(pNumDenTrueAll[i].data(), "", HistType::kTProfile, {ptAxist});
      ue.add(pSumPtTrueAll[i].data(), "", HistType::kTProfile, {ptAxist});
      ue.add(pNumDenTrue[i].data(), "", HistType::kTProfile, {ptAxist});
      ue.add(pSumPtTrue[i].data(), "", HistType::kTProfile, {ptAxist});
      ue.add(pNumDenTruePS[i].data(), "", HistType::kTProfile, {ptAxist});
      ue.add(pSumPtTruePS[i].data(), "", HistType::kTProfile, {ptAxist});
    }
    for (int i = 0; i < 3; ++i) {
      ue.add(hPhiTrue[i].data(), "all charged true; #Delta#phi; Counts",
             HistType::kTH1D,
             {{64, -o2::constants::math::PI / 2.0,
               3.0 * o2::constants::math::PI / 2.0, ""}});
    }
  }

  ue.add("hStat", "TotalEvents", HistType::kTH1F, {{1, 0.5, 1.5, " "}});
  ue.add("hmultRec", "mult rec", HistType::kTH1F, {{200, -0.5, 199.5, " "}});
  ue.add("hdNdeta", "dNdeta", HistType::kTH1F, {{50, -2.5, 2.5, " "}});
  ue.add("vtxZEta", ";#eta;vtxZ", HistType::kTH2F,
         {{50, -2.5, 2.5, " "}, {60, -30, 30, " "}});
  ue.add("phiEta", ";#eta;#varphi", HistType::kTH2F,
         {{50, -2.5, 2.5}, {200, 0., 2 * o2::constants::math::PI, " "}});
  ue.add("hvtxZ", "vtxZ", HistType::kTH1F, {{40, -20.0, 20.0, " "}});

  ue.add("hCounter", "Counter; sel; Nev", HistType::kTH1D, {{7, 0, 7, " "}});
  ue.add("hPtLeadingRecPS", "rec pTleading after physics selection",
         HistType::kTH1D, {ptAxist});
  ue.add("hPtLeadingMeasured", "measured pTleading after physics selection",
         HistType::kTH1D, {ptAxist});
  ue.add("hPtLeadingVsTracks", "", HistType::kTProfile, {{ptAxist}});

  for (int i = 0; i < 3; ++i) {
    ue.add(pNumDenMeasuredPS[i].data(),
           "Number Density; ; #LT #it{N}_{trk} #GT", HistType::kTProfile,
           {ptAxist});
    ue.add(pSumPtMeasuredPS[i].data(),
           "Total #it{p}_{T}; ; #LT#sum#it{p}_{T}#GT", HistType::kTProfile,
           {ptAxist});
    ue.add(hPhi[i].data(), "all charged; #Delta#phi; Counts", HistType::kTH1D,
           {{64, -o2::constants::math::PI / 2.0,
             3.0 * o2::constants::math::PI / 2.0, ""}});
  }

  // Data driven
  for (int i = 0; i < 3; ++i) {
    ue.add(hNumDenMCDd[i].data(), " ", HistType::kTH2D,
           {{ptAxist}, {100, -0.5, 99.5, "#it{N}_{trk}"}});
    ue.add(hSumPtMCDd[i].data(), " ", HistType::kTH2D, {{ptAxist}, {ptAxis}});
    ue.add(hNumDenMCMatchDd[i].data(), " ", HistType::kTH2D,
           {{ptAxist}, {100, -0.5, 99.5, "#it{N}_{trk}"}});
    ue.add(hSumPtMCMatchDd[i].data(), " ", HistType::kTH2D,
           {{ptAxist}, {ptAxis}});
  }

  for (int i = 0; i < 3; ++i) {
    ue.add(hPtVsPtLeadingData[i].data(), " ", HistType::kTH2D,
           {{ptAxist}, {ptAxis}});
    ue.add(pNumDenData[i].data(), "", HistType::kTProfile, {ptAxist});
    ue.add(pSumPtData[i].data(), "", HistType::kTProfile, {ptAxist});
  }

  ue.add("hPtLeadingData", " ", HistType::kTH1D, {{ptAxist}});
  ue.add("hPTVsDCAData", " ", HistType::kTH2D,
         {{ptAxis}, {121, -3.025, 3.025, "#it{DCA}_{xy} (cm)"}});
  ue.add("hEtaLeadingVsPtLeading", " ", HistType::kTH2D,
         {{ptAxist}, {50, -2.5, 2.5, "#eta"}});

  if (analyzeEvandTracksel) {

    const AxisSpec axisVtxZ{500, -25., 25., ""};
    ue.add("hVtxFT0VsVtxCol", " ", HistType::kTH2D,
           {{axisVtxZ}, {axisVtxZ}}); // FT0-vertex vs z-vertex from collisions
    ue.add("hVtxFT0VsVtxCol_afterSel8", " ", HistType::kTH2D,
           {{axisVtxZ}, {axisVtxZ}});
    ue.add("hVtxFT0VsVtxCol_afterPile", " ", HistType::kTH2D,
           {{axisVtxZ}, {axisVtxZ}});
    ue.add("hVtxFT0VsVtxCol_afterGoodZvtx", " ", HistType::kTH2D,
           {{axisVtxZ}, {axisVtxZ}});
    ue.add("hvtxZ_before", "vtxZ befer ev selection", HistType::kTH1F,
           {{40, -20.0, 20.0, " "}});
    ue.add("hvtxZ_after", "vtxZ befer ev after", HistType::kTH1F,
           {{40, -20.0, 20.0, " "}});

    const AxisSpec axisMultT0M{1000, 0., 8000., "T0M multiplicity"};
    ue.add("hVtxFT0MinusVtxColVsMultT0M", "", kTH2F,
           {{axisVtxZ}, {axisMultT0M}}); // FT0-vertex minus z-vertex from
                                         // collisions vs multiplicity
    ue.add("postselection_track/hT0MVsTracks", "", HistType::kTH2D,
           {{axisMultT0M}, {200, 0., 200.}});
    ue.add("postselection_track/hVtxFT0VsTracks", "", HistType::kTH2D,
           {{axisVtxZ}, {200, 0., 200.}});
    ue.add("postselection_track/hVtxVsTracks", "", HistType::kTH2D,
           {{axisVtxZ}, {200, 0., 200.}});

    // its histograms
    ue.add("preselection_track/ITS/itsNCls",
           "number of found ITS clusters;# clusters ITS", kTH1D,
           {{8, -0.5, 7.5}});
    ue.add("preselection_track/ITS/itsChi2NCl",
           "chi2 per ITS cluster;chi2 / cluster ITS", kTH1D, {{100, 0, 40}});
    ue.add("preselection_track/ITS/itsClusterMap", "ITS cluster map", kTH1D,
           {{128, -0.5, 127.5}});
    ue.add("postselection_track/ITS/itsNCls",
           "number of found ITS clusters;# clusters ITS", kTH1D,
           {{8, -0.5, 7.5}});
    ue.add("postselection_track/ITS/itsChi2NCl",
           "chi2 per ITS cluster;chi2 / cluster ITS", kTH1D, {{100, 0, 40}});
    ue.add("postselection_track/ITS/itsClusterMap", "ITS cluster map", kTH1D,
           {{128, -0.5, 127.5}});

    // tpc histograms
    ue.add("preselection_track/TPC/tpcNClsFindable",
           "number of findable TPC clusters;# findable clusters TPC", kTH1D,
           {{165, -0.5, 164.5}});
    ue.add("preselection_track/TPC/tpcNClsFound",
           "number of found TPC clusters;# clusters TPC", kTH1D,
           {{165, -0.5, 164.5}});
    ue.add("preselection_track/TPC/tpcNClsShared",
           "number of shared TPC clusters;# shared clusters TPC", kTH1D,
           {{165, -0.5, 164.5}});
    ue.add("preselection_track/TPC/tpcCrossedRows",
           "number of crossed TPC rows;# crossed rows TPC", kTH1D,
           {{165, -0.5, 164.5}});
    ue.add("preselection_track/TPC/tpcFractionSharedCls",
           "fraction of shared TPC clusters;fraction shared clusters TPC",
           kTH1D, {{100, 0., 1.}});
    ue.add("preselection_track/TPC/tpcCrossedRowsOverFindableCls",
           "crossed TPC rows over findable clusters;crossed rows / findable "
           "clusters TPC",
           kTH1D, {{60, 0.7, 1.3}});
    ue.add("preselection_track/TPC/tpcChi2NCl",
           "chi2 per cluster in TPC;chi2 / cluster TPC", kTH1D, {{100, 0, 10}});
    ue.add("postselection_track/TPC/tpcNClsFindable",
           "number of findable TPC clusters;# findable clusters TPC", kTH1D,
           {{165, -0.5, 164.5}});
    ue.add("postselection_track/TPC/tpcNClsFound",
           "number of found TPC clusters;# clusters TPC", kTH1D,
           {{165, -0.5, 164.5}});
    ue.add("postselection_track/TPC/tpcNClsShared",
           "number of shared TPC clusters;# shared clusters TPC", kTH1D,
           {{165, -0.5, 164.5}});
    ue.add("postselection_track/TPC/tpcCrossedRows",
           "number of crossed TPC rows;# crossed rows TPC", kTH1D,
           {{165, -0.5, 164.5}});
    ue.add("postselection_track/TPC/tpcFractionSharedCls",
           "fraction of shared TPC clusters;fraction shared clusters TPC",
           kTH1D, {{100, 0., 1.}});
    ue.add("postselection_track/TPC/tpcCrossedRowsOverFindableCls",
           "crossed TPC rows over findable clusters;crossed rows / findable "
           "clusters TPC",
           kTH1D, {{60, 0.7, 1.3}});
    ue.add("postselection_track/TPC/tpcChi2NCl",
           "chi2 per cluster in TPC;chi2 / cluster TPC", kTH1D, {{100, 0, 10}});

    // general
    ue.add("preselection_track/hvtxZ", "vtxZ before track selection",
           HistType::kTH1D, {{40, -20.0, 20.0, " "}});
    ue.add("preselection_track/hvtxXY", "vtxXY before track selection",
           HistType::kTH1D, {{121, -3.025, 3.025, "#it{DCA}_{xy} (cm)"}});
    ue.add("preselection_track/htracks", "tracks before track selection",
           HistType::kTH1D, {{100, 0., 100., "N_{tracks}"}});
    ue.add("postselection_track/hvtxZ", "vtxZ after track selection",
           HistType::kTH1D, {{40, -20.0, 20.0, "#it{DCA}_{z} (cm) "}});
    ue.add("postselection_track/hvtxXY", "vtxXY after track selection",
           HistType::kTH1D, {{121, -3.025, 3.025, "#it{DCA}_{xy} (cm)"}});
    ue.add("postselection_track/htracks", "tracks after track selection",
           HistType::kTH1D, {{100, 0., 100., "N_{tracks}"}});
  }
}

void ueCharged::processMC(CollisionTableMCTrue::iterator const& mcCollision,
                          CollisionTableMC const& collisions,
                          TrackTableMC const& tracks,
                          ParticleTableMC const& particles, aod::FT0s const&,
                          BCsRun3 const&)
{
  if (collisions.size() != 0) {
    for (const auto& collision : collisions) {
      auto curTracks = tracks.sliceBy(perCollision, collision.globalIndex());
      processMeasMC<true>(collision, curTracks, particles);
      if (analyzeEvandTracksel) {
        analyzeEventAndTrackSelection(collision, curTracks);
      }
      break; // for now look only at first collision...
    }
  }
  processTrue(mcCollision, particles);
}

void ueCharged::processDataMC(CollisionTableMCData::iterator const& collision,
                              TrackTableMCData const& tracks,
                              ParticleTableMC const& particles,
                              aod::McCollisions const& /*mcCollisions*/)
{
  processMeasMC<false>(collision, tracks, particles);
}

void ueCharged::processData(CollisionTableData::iterator const& collision,
                            TrackTableData const& tracks, aod::FT0s const&,
                            BCsRun3 const&)
{
  processMeas(collision, tracks);
  if (analyzeEvandTracksel) {
    analyzeEventAndTrackSelection(collision, tracks);
  }
}

template <typename C, typename P>
void ueCharged::processTrue(const C& mcCollision, const P& particles)
{
  int multTrue = 0;
  int multTrueINEL = 0;
  for (const auto& particle : particles) {
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

  for (const auto& particle : particles) {
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
  std::vector<double> ueTrue;
  ueTrue.clear();
  int nchmTopTrue[3];
  double sumptmTopTrue[3];
  for (int i = 0; i < 3; ++i) {
    nchmTopTrue[i] = 0;
    sumptmTopTrue[i] = 0;
  }
  std::vector<float> ptArrayTrue;
  std::vector<float> phiArrayTrue;
  std::vector<int> indexArrayTrue;

  for (const auto& particle : particles) {

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
    double dPhi = deltaPhi(particle.phi(), flPhiTrue);

    // definition of the topological regions
    if (std::abs(dPhi) < o2::constants::math::PI / 3.0) { // near side
      ue.fill(HIST(hPhiTrue[0]), dPhi);
      ue.fill(HIST(hPtVsPtLeadingTrue[0]), flPtTrue, particle.pt());
      nchmTopTrue[0]++;
      sumptmTopTrue[0] += particle.pt();
    } else if (std::abs(dPhi - o2::constants::math::PI) <
               o2::constants::math::PI / 3.0) { // away side
      ue.fill(HIST(hPhiTrue[1]), dPhi);
      ue.fill(HIST(hPtVsPtLeadingTrue[1]), flPtTrue, particle.pt());
      nchmTopTrue[1]++;
      sumptmTopTrue[1] += particle.pt();
    } else { // transverse side
      ue.fill(HIST(hPhiTrue[2]), dPhi);
      ue.fill(HIST(hPtVsPtLeadingTrue[2]), flPtTrue, particle.pt());
      nchmTopTrue[2]++;
      sumptmTopTrue[2] += particle.pt();
    }
  }

  for (int i_reg = 0; i_reg < 3; ++i_reg) {
    ueTrue.push_back(1.0 * nchmTopTrue[i_reg]);
  }
  for (int i_reg = 0; i_reg < 3; ++i_reg) {
    ueTrue.push_back(sumptmTopTrue[i_reg]);
  }
  ue.fill(HIST(pNumDenTrueAll[0]), flPtTrue, ueTrue[0]);
  ue.fill(HIST(pSumPtTrueAll[0]), flPtTrue, ueTrue[3]);

  ue.fill(HIST(pNumDenTrueAll[1]), flPtTrue, ueTrue[1]);
  ue.fill(HIST(pSumPtTrueAll[1]), flPtTrue, ueTrue[4]);

  ue.fill(HIST(pNumDenTrueAll[2]), flPtTrue, ueTrue[2]);
  ue.fill(HIST(pSumPtTrueAll[2]), flPtTrue, ueTrue[5]);

  ptArrayTrue.clear();
  phiArrayTrue.clear();
  indexArrayTrue.clear();
}

template <typename C, typename T>
void ueCharged::processMeas(const C& collision, const T& tracks)
{

  ue.fill(HIST("hCounter"), 0);

  ue.fill(HIST("hCounter"), 1);

  if (sel8 && !collision.sel8()) {
    return;
  }

  if (removeITSROFBorder &&
      (!collision.selection_bit(o2::aod::evsel::kIsTriggerTVX) ||
       !collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder))) {
    return;
  }

  ue.fill(HIST("hCounter"), 2);

  if (piluprejection &&
      !collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
    return;
  }

  ue.fill(HIST("hCounter"), 3);
  if (goodzvertex &&
      !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
    return;
  }

  ue.fill(HIST("hCounter"), 4);

  ue.fill(HIST("hStat"), collision.size());
  auto vtxZ = collision.posZ();
  // + vtx
  if ((std::abs(collision.posZ()) >= 10.f)) {
    return;
  }

  ue.fill(HIST("hCounter"), 5);

  ue.fill(HIST("hvtxZ"), vtxZ);

  // loop over selected tracks
  double flPt = 0; // leading pT
  double flPhi = 0;
  double flEta = 0;
  int flIndex = 0;
  int multRec = 0;
  int track_multiplicity = 0;
  for (const auto& track : tracks) {
    if (!mySelectionPrim.IsSelected(track)) {
      continue;
    }
    track_multiplicity++;
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
  ue.fill(HIST("hPtLeadingVsTracks"), flPt, track_multiplicity);

  std::vector<double> ueRec;
  ueRec.clear();
  int nchmTop[3];
  double sumptmTop[3];
  for (int i = 0; i < 3; ++i) {
    nchmTop[i] = 0;
    sumptmTop[i] = 0;
  }
  std::vector<float> ptArray;
  std::vector<float> phiArray;
  std::vector<int> indexArray;

  for (const auto& track : tracks) {

    if (mySelectionOpenDCA.IsSelected(track)) {
      ue.fill(HIST("hPTVsDCAData"), track.pt(), track.dcaXY());
    }

    if (mySelectionPrim.IsSelected(track)) {
      // applying the efficiency twice for the misrec of leading particle
      if (fEff->Eval(track.pt()) > gRandom->Uniform(0, 1)) {
        ptArray.push_back(track.pt());
        phiArray.push_back(track.phi());
        indexArray.push_back(track.globalIndex());
      }

      // remove the autocorrelation
      if (flIndex == track.globalIndex()) {
        continue;
      }

      double dPhi = deltaPhi(track.phi(), flPhi);

      // definition of the topological regions
      if (std::abs(dPhi) < o2::constants::math::PI / 3.0) { // near side
        ue.fill(HIST(hPhi[0]), dPhi);
        ue.fill(HIST(hPtVsPtLeadingData[0]), flPt, track.pt());
        nchmTop[0]++;
        sumptmTop[0] += track.pt();
      } else if (std::abs(dPhi - o2::constants::math::PI) <
                 o2::constants::math::PI / 3.0) { // away side
        ue.fill(HIST(hPhi[1]), dPhi);
        ue.fill(HIST(hPtVsPtLeadingData[1]), flPt, track.pt());
        nchmTop[1]++;
        sumptmTop[1] += track.pt();
      } else { // transverse side
        ue.fill(HIST(hPhi[2]), dPhi);
        ue.fill(HIST(hPtVsPtLeadingData[2]), flPt, track.pt());
        nchmTop[2]++;
        sumptmTop[2] += track.pt();
      }
    }
  }
  for (int i_reg = 0; i_reg < 3; ++i_reg) {
    ueRec.push_back(1.0 * nchmTop[i_reg]);
  }
  for (int i_reg = 0; i_reg < 3; ++i_reg) {
    ueRec.push_back(sumptmTop[i_reg]);
  }

  // add flags for Vtx, PS, ev sel
  ue.fill(HIST(pNumDenMeasuredPS[0]), flPt, ueRec[0]);
  ue.fill(HIST(pNumDenData[0]), flPt, ueRec[0]);
  ue.fill(HIST(pSumPtMeasuredPS[0]), flPt, ueRec[3]);
  ue.fill(HIST(pSumPtData[0]), flPt, ueRec[3]);

  ue.fill(HIST(pNumDenMeasuredPS[1]), flPt, ueRec[1]);
  ue.fill(HIST(pNumDenData[1]), flPt, ueRec[1]);
  ue.fill(HIST(pSumPtMeasuredPS[1]), flPt, ueRec[4]);
  ue.fill(HIST(pSumPtData[1]), flPt, ueRec[4]);

  ue.fill(HIST(pNumDenMeasuredPS[2]), flPt, ueRec[2]);
  ue.fill(HIST(pNumDenData[2]), flPt, ueRec[2]);
  ue.fill(HIST(pSumPtMeasuredPS[2]), flPt, ueRec[5]);
  ue.fill(HIST(pSumPtData[2]), flPt, ueRec[5]);

  ue.fill(HIST("hPtLeadingData"), flPt);

  // Compute data driven (DD) missidentification correction
  float flPtdd = 0; // leading pT
  float flPhidd = 0;
  int flIndexdd = 0;
  int ntrkdd = ptArray.size();

  for (int i = 0; i < ntrkdd; ++i) {
    if (flPtdd < ptArray[i]) {
      flPtdd = ptArray[i];
      flPhidd = phiArray[i];
      flIndexdd = indexArray[i];
    }
  }
  int nchmTopdd[3];
  double sumptmTopdd[3];
  for (int i = 0; i < 3; ++i) {
    nchmTopdd[i] = 0;
    sumptmTopdd[i] = 0;
  }
  for (int i = 0; i < ntrkdd; ++i) {
    if (indexArray[i] == flIndexdd) {
      continue;
    }
    double dPhi = deltaPhi(phiArray[i], flPhidd);
    if (std::abs(dPhi) < o2::constants::math::PI / 3.0) { // near side
      nchmTopdd[0]++;
      sumptmTopdd[0] += ptArray[i];
    } else if (std::abs(dPhi - o2::constants::math::PI) <
               o2::constants::math::PI / 3.0) { // away side
      nchmTopdd[1]++;
      sumptmTopdd[1] += ptArray[i];
    } else { // transverse side
      nchmTopdd[2]++;
      sumptmTopdd[2] += ptArray[i];
    }
  }

  ue.fill(HIST(hNumDenMCDd[0]), flPtdd, nchmTopdd[0]);
  ue.fill(HIST(hSumPtMCDd[0]), flPtdd, sumptmTopdd[0]);

  ue.fill(HIST(hNumDenMCDd[1]), flPtdd, nchmTopdd[1]);
  ue.fill(HIST(hSumPtMCDd[1]), flPtdd, sumptmTopdd[1]);

  ue.fill(HIST(hNumDenMCDd[2]), flPtdd, nchmTopdd[2]);
  ue.fill(HIST(hSumPtMCDd[2]), flPtdd, sumptmTopdd[2]);

  if (flIndexdd == flIndex) {
    ue.fill(HIST(hNumDenMCMatchDd[0]), flPtdd, nchmTopdd[0]);
    ue.fill(HIST(hSumPtMCMatchDd[0]), flPtdd, sumptmTopdd[0]);

    ue.fill(HIST(hNumDenMCMatchDd[1]), flPtdd, nchmTopdd[1]);
    ue.fill(HIST(hSumPtMCMatchDd[1]), flPtdd, sumptmTopdd[1]);

    ue.fill(HIST(hNumDenMCMatchDd[2]), flPtdd, nchmTopdd[2]);
    ue.fill(HIST(hSumPtMCMatchDd[2]), flPtdd, sumptmTopdd[2]);
  }
  ptArray.clear();
  phiArray.clear();
  indexArray.clear();
}

template <bool IS_MC, typename C, typename T, typename P>
void ueCharged::processMeasMC(const C& collision, const T& tracks,
                              const P& particles)
{

  ue.fill(HIST("hStat"), collision.size());
  auto vtxZ = collision.posZ();

  // int multTrue = 0;
  double flPtTrue = 0; // leading pT
  double flPhiTrue = 0;
  int flIndexTrue = 0;

  for (const auto& particle : particles) {
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
  std::vector<double> ueTrue;
  ueTrue.clear();
  int nchmTopTrue[3];
  double sumptmTopTrue[3];
  for (int i = 0; i < 3; ++i) {
    nchmTopTrue[i] = 0;
    sumptmTopTrue[i] = 0;
  }
  std::vector<float> ptArrayTrue;
  std::vector<float> phiArrayTrue;
  std::vector<int> indexArrayTrue;

  for (const auto& particle : particles) {

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
    double dPhi = deltaPhi(particle.phi(), flPhiTrue);

    // definition of the topological regions
    if (std::abs(dPhi) < o2::constants::math::PI / 3.0) { // near side
      nchmTopTrue[0]++;
      sumptmTopTrue[0] += particle.pt();
    } else if (std::abs(dPhi - o2::constants::math::PI) <
               o2::constants::math::PI / 3.0) { // away side
      nchmTopTrue[1]++;
      sumptmTopTrue[1] += particle.pt();
    } else { // transverse side
      nchmTopTrue[2]++;
      sumptmTopTrue[2] += particle.pt();
    }
  }

  for (int i_reg = 0; i_reg < 3; ++i_reg) {
    ueTrue.push_back(1.0 * nchmTopTrue[i_reg]);
  }
  for (int i_reg = 0; i_reg < 3; ++i_reg) {
    ueTrue.push_back(sumptmTopTrue[i_reg]);
  }

  ue.fill(HIST(pNumDenTrue[0]), flPtTrue, ueTrue[0]);
  ue.fill(HIST(pSumPtTrue[0]), flPtTrue, ueTrue[3]);

  ue.fill(HIST(pNumDenTrue[1]), flPtTrue, ueTrue[1]);
  ue.fill(HIST(pSumPtTrue[1]), flPtTrue, ueTrue[4]);

  ue.fill(HIST(pNumDenTrue[2]), flPtTrue, ueTrue[2]);
  ue.fill(HIST(pSumPtTrue[2]), flPtTrue, ueTrue[5]);

  ptArrayTrue.clear();
  phiArrayTrue.clear();
  indexArrayTrue.clear();

  ue.fill(HIST("hCounter"), 0);

  ue.fill(HIST("hCounter"), 1);

  if (sel8 && !collision.sel8()) {
    return;
  }

  if (removeITSROFBorder &&
      (!collision.selection_bit(o2::aod::evsel::kIsTriggerTVX) ||
       !collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder))) {
    return;
  }

  ue.fill(HIST("hCounter"), 2);

  if (piluprejection &&
      !collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
    return;
  }
  ue.fill(HIST("hCounter"), 3);

  if (goodzvertex &&
      !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
    return;
  }
  ue.fill(HIST("hCounter"), 4);

  // only PS
  if ((std::abs(collision.posZ()) >= 10.f)) {
    return;
  }
  ue.fill(HIST("hCounter"), 5);

  ue.fill(HIST(pNumDenTruePS[0]), flPtTrue, ueTrue[0]);
  ue.fill(HIST(pSumPtTruePS[0]), flPtTrue, ueTrue[3]);

  ue.fill(HIST(pNumDenTruePS[1]), flPtTrue, ueTrue[1]);
  ue.fill(HIST(pSumPtTruePS[1]), flPtTrue, ueTrue[4]);

  ue.fill(HIST(pNumDenTruePS[2]), flPtTrue, ueTrue[2]);
  ue.fill(HIST(pSumPtTruePS[2]), flPtTrue, ueTrue[5]);

  ue.fill(HIST("hvtxZ"), vtxZ);
  // loop over MC true particles
  int multTrue = 0;
  for (const auto& particle : particles) {
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
  int track_multiplicity = 0;

  for (const auto& track : tracks) {
    if (!mySelectionPrim.IsSelected(track)) {
      continue;
    }
    track_multiplicity++;
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

  ue.fill(HIST("hPtLeadingVsTracks"), flPt, track_multiplicity);
  ue.fill(HIST("hmultRec"), multRec);
  ue.fill(HIST("hPtLeadingMeasured"), flPt);
  ue.fill(HIST("hPtLeadingRecPS"), flPt);
  std::vector<double> ueRec;
  ueRec.clear();
  int nchmTop[3];
  double sumptmTop[3];
  for (int i = 0; i < 3; ++i) {
    nchmTop[i] = 0;
    sumptmTop[i] = 0;
  }
  std::vector<float> ptArray;
  std::vector<float> phiArray;
  std::vector<int> indexArray;

  for (const auto& track : tracks) {

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
        const auto& particle = track.mcParticle();
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
            if (particle
                  .producedByGenerator()) { // i guess these are from decays
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
      if (fEff->Eval(track.pt()) > gRandom->Uniform(0, 1)) {
        ptArray.push_back(track.pt());
        phiArray.push_back(track.phi());
        indexArray.push_back(track.globalIndex());
      }

      // remove the autocorrelation
      if (flIndex == track.globalIndex()) {
        continue;
      }

      double dPhi = deltaPhi(track.phi(), flPhi);

      // definition of the topological regions
      if (std::abs(dPhi) < o2::constants::math::PI / 3.0) { // near side
        ue.fill(HIST(hPhi[0]), dPhi);
        ue.fill(HIST(hPtVsPtLeadingData[0]), flPt, track.pt());
        nchmTop[0]++;
        sumptmTop[0] += track.pt();
      } else if (std::abs(dPhi - o2::constants::math::PI) <
                 o2::constants::math::PI / 3.0) { // away side
        ue.fill(HIST(hPhi[1]), dPhi);
        ue.fill(HIST(hPtVsPtLeadingData[1]), flPt, track.pt());
        nchmTop[1]++;
        sumptmTop[1] += track.pt();
      } else { // transverse side
        ue.fill(HIST(hPhi[2]), dPhi);
        ue.fill(HIST(hPtVsPtLeadingData[2]), flPt, track.pt());
        nchmTop[2]++;
        sumptmTop[2] += track.pt();
      }
    }
  }
  for (int i_reg = 0; i_reg < 3; ++i_reg) {
    ueRec.push_back(1.0 * nchmTop[i_reg]);
  }
  for (int i_reg = 0; i_reg < 3; ++i_reg) {
    ueRec.push_back(sumptmTop[i_reg]);
  }

  // add flags for Vtx, PS, ev sel
  ue.fill(HIST(pNumDenMeasuredPS[0]), flPt, ueRec[0]);
  ue.fill(HIST(pNumDenData[0]), flPt, ueRec[0]);
  ue.fill(HIST(pSumPtMeasuredPS[0]), flPt, ueRec[3]);
  ue.fill(HIST(pSumPtData[0]), flPt, ueRec[3]);

  ue.fill(HIST(pNumDenMeasuredPS[1]), flPt, ueRec[1]);
  ue.fill(HIST(pNumDenData[1]), flPt, ueRec[1]);
  ue.fill(HIST(pSumPtMeasuredPS[1]), flPt, ueRec[4]);
  ue.fill(HIST(pSumPtData[1]), flPt, ueRec[4]);

  ue.fill(HIST(pNumDenMeasuredPS[2]), flPt, ueRec[2]);
  ue.fill(HIST(pNumDenData[2]), flPt, ueRec[2]);
  ue.fill(HIST(pSumPtMeasuredPS[2]), flPt, ueRec[5]);
  ue.fill(HIST(pSumPtData[2]), flPt, ueRec[5]);

  ue.fill(HIST("hPtLeadingData"), flPt);

  // Compute data driven (DD) missidentification correction
  float flPtdd = 0; // leading pT
  float flPhidd = 0;
  int flIndexdd = 0;
  int ntrkdd = ptArray.size();

  for (int i = 0; i < ntrkdd; ++i) {
    if (flPtdd < ptArray[i]) {
      flPtdd = ptArray[i];
      flPhidd = phiArray[i];
      flIndexdd = indexArray[i];
    }
  }
  int nchmTopdd[3];
  double sumptmTopdd[3];
  for (int i = 0; i < 3; ++i) {
    nchmTopdd[i] = 0;
    sumptmTopdd[i] = 0;
  }
  for (int i = 0; i < ntrkdd; ++i) {
    if (indexArray[i] == flIndexdd) {
      continue;
    }
    double dPhi = deltaPhi(phiArray[i], flPhidd);
    if (std::abs(dPhi) < o2::constants::math::PI / 3.0) { // near side
      nchmTopdd[0]++;
      sumptmTopdd[0] += ptArray[i];
    } else if (std::abs(dPhi - o2::constants::math::PI) <
               o2::constants::math::PI / 3.0) { // away side
      nchmTopdd[1]++;
      sumptmTopdd[1] += ptArray[i];
    } else { // transverse side
      nchmTopdd[2]++;
      sumptmTopdd[2] += ptArray[i];
    }
  }

  ue.fill(HIST(hNumDenMCDd[0]), flPtdd, nchmTopdd[0]);
  ue.fill(HIST(hSumPtMCDd[0]), flPtdd, sumptmTopdd[0]);

  ue.fill(HIST(hNumDenMCDd[1]), flPtdd, nchmTopdd[1]);
  ue.fill(HIST(hSumPtMCDd[1]), flPtdd, sumptmTopdd[1]);

  ue.fill(HIST(hNumDenMCDd[2]), flPtdd, nchmTopdd[2]);
  ue.fill(HIST(hSumPtMCDd[2]), flPtdd, sumptmTopdd[2]);

  if (flIndexdd == flIndex) {
    ue.fill(HIST(hNumDenMCMatchDd[0]), flPtdd, nchmTopdd[0]);
    ue.fill(HIST(hSumPtMCMatchDd[0]), flPtdd, sumptmTopdd[0]);

    ue.fill(HIST(hNumDenMCMatchDd[1]), flPtdd, nchmTopdd[1]);
    ue.fill(HIST(hSumPtMCMatchDd[1]), flPtdd, sumptmTopdd[1]);

    ue.fill(HIST(hNumDenMCMatchDd[2]), flPtdd, nchmTopdd[2]);
    ue.fill(HIST(hSumPtMCMatchDd[2]), flPtdd, sumptmTopdd[2]);
  }
  ptArray.clear();
  phiArray.clear();
  indexArray.clear();
}

template <typename C, typename T>
void ueCharged::analyzeEventAndTrackSelection(const C& collision,
                                              const T& tracks)
{

  // z-vertex from FT0 vs PV analysis

  const auto& foundBC = collision.template foundBC_as<BCsRun3>();

  if (foundBC.has_ft0()) {
    ue.fill(HIST("hVtxFT0VsVtxCol"), foundBC.ft0().posZ(), collision.posZ());
  }

  ue.fill(HIST("hvtxZ_before"), collision.posZ());

  if (sel8 && !collision.sel8()) {
    return;
  }

  if (removeITSROFBorder &&
      (!collision.selection_bit(o2::aod::evsel::kIsTriggerTVX) ||
       !collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder))) {
    return;
  }

  if (foundBC.has_ft0()) {
    ue.fill(HIST("hVtxFT0VsVtxCol_afterSel8"), foundBC.ft0().posZ(),
            collision.posZ());
  }

  if (piluprejection &&
      !collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
    return;
  }

  if (foundBC.has_ft0()) {
    ue.fill(HIST("hVtxFT0VsVtxCol_afterPile"), foundBC.ft0().posZ(),
            collision.posZ());
  }

  if (goodzvertex &&
      !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
    return;
  }

  if (foundBC.has_ft0()) {
    ue.fill(HIST("hVtxFT0VsVtxCol_afterGoodZvtx"), foundBC.ft0().posZ(),
            collision.posZ());
  }

  ue.fill(HIST("hvtxZ_after"), collision.posZ());

  // analysis of the track selection

  int tracks_before = 0;
  int tracks_after = 0;

  for (auto& track : tracks) {

    if (track.hasITS() && track.hasTPC()) {
      ue.fill(HIST("preselection_track/ITS/itsNCls"), track.itsNCls());
      ue.fill(HIST("preselection_track/ITS/itsChi2NCl"), track.itsChi2NCl());
      ue.fill(HIST("preselection_track/ITS/itsClusterMap"),
              track.itsClusterMap());
      ue.fill(HIST("preselection_track/TPC/tpcNClsFindable"),
              track.tpcNClsFindable());
      ue.fill(HIST("preselection_track/TPC/tpcNClsFound"),
              track.tpcNClsFound());
      ue.fill(HIST("preselection_track/TPC/tpcNClsShared"),
              track.tpcNClsShared());
      ue.fill(HIST("preselection_track/TPC/tpcCrossedRows"),
              track.tpcNClsCrossedRows());
      ue.fill(HIST("preselection_track/TPC/tpcCrossedRowsOverFindableCls"),
              track.tpcCrossedRowsOverFindableCls());
      ue.fill(HIST("preselection_track/TPC/tpcFractionSharedCls"),
              track.tpcFractionSharedCls());
      ue.fill(HIST("preselection_track/TPC/tpcChi2NCl"), track.tpcChi2NCl());
      ue.fill(HIST("preselection_track/hvtxZ"), track.dcaZ());
      ue.fill(HIST("preselection_track/hvtxXY"), track.dcaXY());
      tracks_before++;
    }

    if (mySelectionPrim.IsSelected(track)) {
      if (track.hasITS() && track.hasTPC()) {
        ue.fill(HIST("postselection_track/ITS/itsNCls"), track.itsNCls());
        ue.fill(HIST("postselection_track/ITS/itsChi2NCl"), track.itsChi2NCl());
        ue.fill(HIST("postselection_track/ITS/itsClusterMap"),
                track.itsClusterMap());
        ue.fill(HIST("postselection_track/TPC/tpcNClsFindable"),
                track.tpcNClsFindable());
        ue.fill(HIST("postselection_track/TPC/tpcNClsFound"),
                track.tpcNClsFound());
        ue.fill(HIST("postselection_track/TPC/tpcNClsShared"),
                track.tpcNClsShared());
        ue.fill(HIST("postselection_track/TPC/tpcCrossedRows"),
                track.tpcNClsCrossedRows());
        ue.fill(HIST("postselection_track/TPC/tpcCrossedRowsOverFindableCls"),
                track.tpcCrossedRowsOverFindableCls());
        ue.fill(HIST("postselection_track/TPC/tpcFractionSharedCls"),
                track.tpcFractionSharedCls());
        ue.fill(HIST("postselection_track/TPC/tpcChi2NCl"), track.tpcChi2NCl());
        ue.fill(HIST("postselection_track/hvtxZ"), track.dcaZ());
        ue.fill(HIST("postselection_track/hvtxXY"), track.dcaXY());
        tracks_after++;
      }
    }
  }

  ue.fill(HIST("postselection_track/htracks"), tracks_after);
  ue.fill(HIST("preselection_track/htracks"), tracks_before);

  // FT0
  float multT0A = foundBC.has_ft0() ? foundBC.ft0().sumAmpA() : -999.f;
  float multT0C = foundBC.has_ft0() ? foundBC.ft0().sumAmpC() : -999.f;

  // z-vertex from FT0 vs PV
  if (foundBC.has_ft0()) {
    ue.fill(HIST("hVtxFT0MinusVtxColVsMultT0M"),
            foundBC.ft0().posZ() - collision.posZ(), multT0A + multT0C);
  }

  ue.fill(HIST("postselection_track/hT0MVsTracks"), multT0A + multT0C,
          tracks_after);
  ue.fill(HIST("postselection_track/hVtxFT0VsTracks"), foundBC.ft0().posZ(),
          tracks_after);
  ue.fill(HIST("postselection_track/hVtxVsTracks"), collision.posZ(),
          tracks_after);
}
