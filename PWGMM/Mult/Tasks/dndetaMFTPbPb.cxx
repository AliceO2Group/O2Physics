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
///
/// \file dndetaMFTPbPb.cxx
/// \brief  Task for calculating dNdeta in Pb-Pb collisions using MFT detector
/// \author Gyula Bencedi, gyula.bencedi@cern.ch
/// \since  Nov 2024

#include <algorithm>
#include <chrono>
#include <cmath>
#include <numeric>
#include <vector>
#include <string>

#include "CCDB/BasicCCDBManager.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/Configurable.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/RuntimeError.h"
#include "Framework/runDataProcessing.h"

#include "Common/CCDB/ctpRateFetcher.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "CommonConstants/MathConstants.h"

#include "MathUtils/Utils.h"
#include "ReconstructionDataFormats/GlobalTrackID.h"
#include "TPDGCode.h"

#include "Functions.h"
#include "Index.h"
#include "bestCollisionTable.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::track;
using namespace o2::aod::fwdtrack;
using namespace o2::constants::physics;
using namespace o2::constants::math;
using namespace pwgmm::mult;

AxisSpec ptAxis = {1001, -0.005, 10.005};
AxisSpec multAxis = {701, -0.5, 700.5, "N_{trk}"};
AxisSpec zAxis = {60, -30., 30.};
AxisSpec deltaZAxis = {61, -6.1, 6.1};
AxisSpec dcaxyAxis = {500, -1, 50};
AxisSpec phiAxis = {629, 0, TwoPI, "Rad", "#phi"};
AxisSpec etaAxis = {20, -4., -2.};
AxisSpec centAxis{100, 0, 100, "centrality"};

struct DndetaMFTPbPb {
  SliceCache cache;

  enum OccupancyEst { TrkITS = 1,
                      Ft0C };

  HistogramRegistry registry{
    "registry",
    {},
    OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry qaregistry{
    "qaregistry",
    {},
    OutputObjHandlingPolicy::AnalysisObject,
    false,
    true};

  struct : ConfigurableGroup {
    Configurable<bool> usephiCut{"usephiCut", false, "use azimuthal angle cut"};
    Configurable<float> phiCut{"phiCut", 0.1f,
                               "Cut on azimuthal angle of MFT tracks"};
    Configurable<float> minPhi{"minPhi", 0.f, ""};
    Configurable<float> maxPhi{"maxPhi", 6.2832, ""};
    Configurable<float> minEta{"minEta", -3.6f, ""};
    Configurable<float> maxEta{"maxEta", -2.5f, ""};
    Configurable<int> minNclusterMft{"minNclusterMft", 5,
                                     "minimum number of MFT clusters"};
    Configurable<double> minPt{"minPt", 0., "minimum pT of the MFT tracks"};
    Configurable<bool> requireCA{
      "requireCA", false, "Use Cellular Automaton track-finding algorithm"};
    Configurable<float> maxDCAxy{"maxDCAxy", 2.0f, "Cut on dcaXY"};
  } trackCuts;

  struct : ConfigurableGroup {
    Configurable<float> maxZvtx{"maxZvtx", 10.0f, "Cut on z-vtx"};
    Configurable<bool> useZDiffCut{"useZDiffCut", false,
                                   "use Zvtx reco-mc diff. cut"};
    Configurable<float> maxZvtxDiff{
      "maxZvtxDiff", 1.0f,
      "max allowed Z vtx difference for reconstruced collisions (cm)"};
    Configurable<bool> requireNoCollInTimeRangeStd{
      "requireNoCollInTimeRangeStd", false,
      "reject collisions corrupted by the cannibalism, with other collisions "
      "within +/- 10 microseconds"};
    Configurable<bool> requireNoCollInTimeRangeNarrow{
      "requireNoCollInTimeRangeNarrow", false,
      "reject collisions corrupted by the cannibalism, with other collisions "
      "within +/- 10 microseconds"};
    Configurable<uint> occupancyEstimator{
      "occupancyEstimator", 1,
      "Occupancy estimator: 1 = trackOccupancyInTimeRange, 2 = "
      "ft0cOccupancyInTimeRange"};
    Configurable<float> minOccupancy{
      "minOccupancy", -1, "minimum occupancy from neighbouring collisions"};
    Configurable<float> maxOccupancy{
      "maxOccupancy", -1, "maximum occupancy from neighbouring collisions"};
    Configurable<float> minIR{"minIR", -1, "minimum IR (kHz) collisions"};
    Configurable<float> maxIR{"maxIR", -1, "maximum IR (kHz) collisions"};
  } eventCuts;

  ConfigurableAxis occupancyBins{"occupancyBins",
                                 {VARIABLE_WIDTH, 0.0f, 250.0f, 500.0f, 750.0f,
                                  1000.0f, 1500.0f, 2000.0f, 3000.0f, 4500.0f,
                                  6000.0f, 8000.0f, 10000.0f, 50000.0f},
                                 "Occupancy"};
  ConfigurableAxis centralityBins{
    "centralityBins",
    {VARIABLE_WIDTH, 0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100},
    "Centrality"};
  ConfigurableAxis irBins{"irBins", {500, 0, 50}, "Interaction rate (kHz)"};

  Service<o2::framework::O2DatabasePDG> pdg;

  Service<ccdb::BasicCCDBManager> ccdb;
  Configurable<int64_t> ccdbNoLaterThan{
    "ccdbNoLaterThan",
    std::chrono::duration_cast<std::chrono::milliseconds>(
      std::chrono::system_clock::now().time_since_epoch())
      .count(),
    "latest acceptable timestamp of creation for the object"};
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch",
                                    "url of the ccdb repository"};
  ctpRateFetcher rateFetcher;

  /// @brief init function, definition of histograms
  void init(InitContext&)
  {
    ccdb->setURL(ccdbUrl.value);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setCreatedNotAfter(ccdbNoLaterThan.value);
    ccdb->setFatalWhenNull(false);

    if (static_cast<int>(doprocessDataInclusive) +
          static_cast<int>(doprocessDatawBestTracksInclusive) >
        1) {
      LOGP(fatal,
           "Either processDataInclusive OR "
           "processDatawBestTracksInclusive should be enabled!");
    }
    if ((static_cast<int>(doprocessDataCentFT0C) +
           static_cast<int>(doprocessDatawBestTracksCentFT0C) >
         1) ||
        (static_cast<int>(doprocessDataCentFT0CVariant1) +
           static_cast<int>(doprocessDatawBestTracksCentFT0CVariant1) >
         1) ||
        (static_cast<int>(doprocessDataCentFT0M) +
           static_cast<int>(doprocessDatawBestTracksCentFT0M) >
         1) ||
        (static_cast<int>(doprocessDataCentNGlobal) +
           static_cast<int>(doprocessDatawBestTracksCentNGlobal) >
         1) ||
        (static_cast<int>(doprocessDataCentMFT) +
           static_cast<int>(doprocessDatawBestTracksCentMFT) >
         1)) {
      LOGP(fatal,
           "Either processDataCent[ESTIMATOR] OR "
           "processDatawBestTracksCent[ESTIMATOR] should "
           "be enabled!");
    }
    if (static_cast<int>(doprocessMCInclusive) +
          static_cast<int>(doprocessMCwBestTracksInclusive) >
        1) {
      LOGP(fatal,
           "Either processMCInclusive OR processMCwBestTracksInclusive "
           "should be enabled!");
    }
    if ((static_cast<int>(doprocessMCCentFT0C) +
           static_cast<int>(doprocessMCwBestTracksCentFT0C) >
         1) ||
        (static_cast<int>(doprocessMCCentFT0CVariant1) +
           static_cast<int>(doprocessMCwBestTracksCentFT0CVariant1) >
         1) ||
        (static_cast<int>(doprocessMCCentFT0M) +
           static_cast<int>(doprocessMCwBestTracksCentFT0M) >
         1) ||
        (static_cast<int>(doprocessMCCentNGlobal) +
           static_cast<int>(doprocessMCwBestTracksCentNGlobal) >
         1) ||
        (static_cast<int>(doprocessMCCentMFT) +
           static_cast<int>(doprocessMCwBestTracksCentMFT) >
         1)) {
      LOGP(fatal,
           "Either processMCCent[ESTIMATOR] OR "
           "processMCwBestTracksCent[ESTIMATOR] should "
           "be enabled!");
    }

    auto hev = registry.add<TH1>("hEvtSel", "hEvtSel", HistType::kTH1F,
                                 {{12, -0.5f, +11.5f}});
    hev->GetXaxis()->SetBinLabel(1, "All collisions");
    hev->GetXaxis()->SetBinLabel(2, "Ev. sel.");
    hev->GetXaxis()->SetBinLabel(3, "kIsGoodZvtxFT0vsPV");
    hev->GetXaxis()->SetBinLabel(4, "NoSameBunchPileup");
    hev->GetXaxis()->SetBinLabel(5, "Z-vtx cut");
    hev->GetXaxis()->SetBinLabel(6, "kNoCollInTimeRangeStd");
    hev->GetXaxis()->SetBinLabel(7, "kNoCollInTimeRangeNarrow");
    hev->GetXaxis()->SetBinLabel(8, "Below min occup.");
    hev->GetXaxis()->SetBinLabel(9, "Above max occup.");
    hev->GetXaxis()->SetBinLabel(10, "Below min IR (kHz)");
    hev->GetXaxis()->SetBinLabel(11, "Above max IR (kHz)");

    auto hBcSel = registry.add<TH1>("hBcSel", "hBcSel", HistType::kTH1F,
                                    {{3, -0.5f, +2.5f}});
    hBcSel->GetXaxis()->SetBinLabel(1, "Good BCs");
    hBcSel->GetXaxis()->SetBinLabel(2, "BCs with collisions");
    hBcSel->GetXaxis()->SetBinLabel(3, "BCs with pile-up/splitting");

    AxisSpec centralityAxis = {centralityBins, "Centrality", "centralityAxis"};
    AxisSpec occupancyAxis = {occupancyBins, "Occupancy", "occupancyAxis"};

    if (doprocessDataInclusive || doprocessDatawBestTracksInclusive ||
        doprocessMCInclusive || doprocessMCwBestTracksInclusive) {
      registry.add({"Events/Selection",
                    ";status;occupancy",
                    {HistType::kTH2F, {{2, 0.5, 2.5}, occupancyAxis}}});
      auto hstat = registry.get<TH2>(HIST("Events/Selection"));
      auto* x = hstat->GetXaxis();
      x->SetBinLabel(1, "All");
      x->SetBinLabel(2, "Selected");

      qaregistry.add("hOccIRate", "hOccIRate", HistType::kTH2F,
                     {occupancyAxis, irBins});

      registry.add({"Events/NtrkZvtx",
                    "; N_{trk}; Z_{vtx} (cm); occupancy",
                    {HistType::kTHnSparseF, {multAxis, zAxis, occupancyAxis}}});
      registry.add({"Tracks/EtaZvtx",
                    "; #eta; Z_{vtx} (cm); occupancy",
                    {HistType::kTHnSparseF, {etaAxis, zAxis, occupancyAxis}}});
      registry.add(
        {"Tracks/PhiEta",
         "; #varphi; #eta; occupancy",
         {HistType::kTHnSparseF, {phiAxis, etaAxis, occupancyAxis}}});

      qaregistry.add(
        {"Tracks/Chi2Eta",
         "; #chi^{2}; #it{#eta}; occupancy",
         {HistType::kTHnSparseF, {{600, 0, 20}, etaAxis, occupancyAxis}}});
      qaregistry.add({"Tracks/Chi2",
                      "; #chi^{2};",
                      {HistType::kTH2F, {{600, 0, 20}, occupancyAxis}}});
      qaregistry.add(
        {"Tracks/NclustersEta",
         "; nClusters; #eta; occupancy",
         {HistType::kTHnSparseF, {{7, 4, 10}, etaAxis, occupancyAxis}}});
      qaregistry.add({"Tracks/NchSel",
                      "; N_{ch}; occupancy",
                      {HistType::kTH2F, {multAxis, occupancyAxis}}});

      if (doprocessDatawBestTracksInclusive) {
        registry.add(
          {"Events/NtrkZvtxBest",
           "; N_{trk}; Z_{vtx} (cm); occupancy",
           {HistType::kTHnSparseF, {multAxis, zAxis, occupancyAxis}}});
        registry.add(
          {"Tracks/EtaZvtxBest",
           "; #eta; Z_{vtx} (cm); occupancy",
           {HistType::kTHnSparseF, {etaAxis, zAxis, occupancyAxis}}});
        registry.add(
          {"Tracks/PhiEtaBest",
           "; #varphi; #eta; occupancy",
           {HistType::kTHnSparseF, {phiAxis, etaAxis, occupancyAxis}}});
        qaregistry.add(
          {"Tracks/NclustersEtaBest",
           "; nClusters; #eta; occupancy",
           {HistType::kTHnSparseF, {{7, 4, 10}, etaAxis, occupancyAxis}}});
        qaregistry.add(
          {"Tracks/DCAXYPt",
           "; p_{T} (GeV/c) ; DCA_{XY} (cm); occupancy",
           {HistType::kTHnSparseF, {ptAxis, dcaxyAxis, occupancyAxis}}});
        qaregistry.add({"Tracks/DCAXY",
                        "; DCA_{XY} (cm); occupancy",
                        {HistType::kTH2F, {dcaxyAxis, occupancyAxis}}});
        qaregistry.add(
          {"Tracks/ReTracksEtaZvtx",
           "; #eta; #it{z}_{vtx} (cm); occupancy",
           {HistType::kTHnSparseF, {etaAxis, zAxis, occupancyAxis}}});
        qaregistry.add(
          {"Tracks/ReTracksPhiEta",
           "; #varphi; #eta; occupancy",
           {HistType::kTHnSparseF, {phiAxis, etaAxis, occupancyAxis}}});
        qaregistry.add({"Tracks/TrackAmbDegree",
                        "; N_{coll}^{comp}; occupancy",
                        {HistType::kTH2F, {{51, -0.5, 50.5}, occupancyAxis}}});
      }
    }

    if (doprocessDataCentFT0C || doprocessDatawBestTracksCentFT0C ||
        doprocessMCCentFT0C || doprocessMCwBestTracksCentFT0C ||
        doprocessDataCentFT0CVariant1 ||
        doprocessDatawBestTracksCentFT0CVariant1 ||
        doprocessMCCentFT0CVariant1 || doprocessMCwBestTracksCentFT0CVariant1 ||
        doprocessDataCentFT0M || doprocessDatawBestTracksCentFT0M ||
        doprocessMCCentFT0M || doprocessMCwBestTracksCentFT0M ||
        doprocessDataCentNGlobal || doprocessDatawBestTracksCentNGlobal ||
        doprocessMCCentNGlobal || doprocessMCwBestTracksCentNGlobal ||
        doprocessDataCentMFT || doprocessDatawBestTracksCentMFT ||
        doprocessMCCentMFT || doprocessMCwBestTracksCentMFT) {
      registry.add({"Events/Centrality/Selection",
                    ";status;centrality;occupancy",
                    {HistType::kTHnSparseF,
                     {{2, 0.5, 2.5}, centralityAxis, occupancyAxis}}});
      auto hstat = registry.get<THnSparse>(HIST("Events/Centrality/Selection"));
      hstat->GetAxis(0)->SetBinLabel(1, "All");
      hstat->GetAxis(0)->SetBinLabel(2, "Selected");

      qaregistry.add("hCentOccIRate", "hCentOccIRate", HistType::kTHnSparseF,
                     {centralityAxis, occupancyAxis, irBins});

      qaregistry.add({"Events/Centrality/hCent",
                      "; centrality; occupancy",
                      {HistType::kTH2F, {centAxis, occupancyAxis}},
                      true});
      qaregistry.add(
        {"Events/Centrality/hZvtxCent",
         "; Z_{vtx} (cm); centrality; occupancy",
         {HistType::kTHnSparseF, {zAxis, centralityAxis, occupancyAxis}}});
      registry.add({"Events/Centrality/NtrkZvtx",
                    "; N_{trk}; Z_{vtx} (cm); centrality; occupancy",
                    {HistType::kTHnSparseF,
                     {multAxis, zAxis, centralityAxis, occupancyAxis}}});
      registry.add({"Tracks/Centrality/EtaZvtx",
                    "; #eta; Z_{vtx} (cm); centrality; occupancy",
                    {HistType::kTHnSparseF,
                     {etaAxis, zAxis, centralityAxis, occupancyAxis}}});
      registry.add({"Tracks/Centrality/PhiEta",
                    "; #varphi; #eta; centrality; occupancy",
                    {HistType::kTHnSparseF,
                     {phiAxis, etaAxis, centralityAxis, occupancyAxis}}});

      qaregistry.add(
        {"Tracks/Centrality/NchSel",
         "; N_{ch}; centrality; occupancy",
         {HistType::kTHnSparseF, {multAxis, centralityAxis, occupancyAxis}}});
      qaregistry.add(
        {"Tracks/Centrality/Chi2Eta",
         "; #chi^{2}; #it{#eta}; centrality; occupancy",
         {HistType::kTHnSparseF,
          {{600, 0, 20}, etaAxis, centralityAxis, occupancyAxis}}});
      qaregistry.add({"Tracks/Centrality/Chi2",
                      "; #chi^{2}; centrality; occupancy",
                      {HistType::kTHnSparseF,
                       {{600, 0, 20}, centralityAxis, occupancyAxis}}});
      qaregistry.add({"Tracks/Centrality/NclustersEta",
                      "; nClusters; #eta; centrality; occupancy",
                      {HistType::kTHnSparseF,
                       {{7, 4, 10}, etaAxis, centralityAxis, occupancyAxis}}});

      if (doprocessDatawBestTracksCentFT0C ||
          doprocessDatawBestTracksCentFT0CVariant1 ||
          doprocessDatawBestTracksCentFT0M ||
          doprocessDatawBestTracksCentNGlobal ||
          doprocessDatawBestTracksCentMFT) {
        registry.add({"Events/Centrality/NtrkZvtxBest",
                      "; N_{trk}; Z_{vtx} (cm); centrality; occupancy",
                      {HistType::kTHnSparseF,
                       {multAxis, zAxis, centralityAxis, occupancyAxis}}});
        registry.add({"Tracks/Centrality/EtaZvtxBest",
                      "; #eta; Z_{vtx} (cm); centrality; occupancy",
                      {HistType::kTHnSparseF,
                       {etaAxis, zAxis, centralityAxis, occupancyAxis}}});
        registry.add({"Tracks/Centrality/PhiEtaBest",
                      "; #varphi; #eta; centrality; occupancy",
                      {HistType::kTHnSparseF,
                       {phiAxis, etaAxis, centralityAxis, occupancyAxis}}});
        qaregistry.add(
          {"Tracks/Centrality/NclustersEtaBest",
           "; nClusters; #eta; centrality; occupancy",
           {HistType::kTHnSparseF,
            {{7, 4, 10}, etaAxis, centralityAxis, occupancyAxis}}});
        qaregistry.add({"Tracks/Centrality/TrackAmbDegree",
                        "; N_{coll}^{comp}; centrality; occupancy",
                        {HistType::kTHnSparseF,
                         {{51, -0.5, 50.5}, centralityAxis, occupancyAxis}}});
        qaregistry.add({"Tracks/Centrality/DCAXY",
                        "; DCA_{XY} (cm); centrality; occupancy",
                        {HistType::kTHnSparseF,
                         {dcaxyAxis, centralityAxis, occupancyAxis}}});
        qaregistry.add(
          {"Tracks/Centrality/DCAXYPt",
           "; p_{T} (GeV/c) ; DCA_{XY} (cm); centrality; occupancy",
           {HistType::kTHnSparseF,
            {ptAxis, dcaxyAxis, centralityAxis, occupancyAxis}}});
        qaregistry.add({"Tracks/Centrality/ReTracksEtaZvtx",
                        "; #eta; #it{z}_{vtx} (cm); occupancy",
                        {HistType::kTHnSparseF,
                         {etaAxis, zAxis, centralityAxis, occupancyAxis}}});
        qaregistry.add({"Tracks/Centrality/ReTracksPhiEta",
                        "; #varphi; #eta; occupancy",
                        {HistType::kTHnSparseF,
                         {phiAxis, etaAxis, centralityAxis, occupancyAxis}}});
      }
    }

    if (doprocessMCInclusive || doprocessMCwBestTracksInclusive) {
      registry.add({"Events/EvtEffGen",
                    ";status;occupancy",
                    {HistType::kTH2F, {{3, 0.5, 3.5}, occupancyAxis}}});
      auto heff = registry.get<TH2>(HIST("Events/EvtEffGen"));
      auto* h = heff->GetXaxis();
      h->SetBinLabel(1, "All reconstructed");
      h->SetBinLabel(2, "Selected reconstructed");
      h->SetBinLabel(3, "All generated");

      registry.add({"Events/NtrkZvtxGen_t",
                    "; N_{trk}; Z_{vtx} (cm);",
                    {HistType::kTH2F, {multAxis, zAxis}}});
      registry.add({"Events/NtrkZvtxGen",
                    "; N_{trk}; Z_{vtx} (cm); occupancy",
                    {HistType::kTHnSparseF, {multAxis, zAxis, occupancyAxis}}});
      registry.add({"Tracks/EtaZvtxGen",
                    "; #eta; Z_{vtx} (cm); occupancy",
                    {HistType::kTHnSparseF, {etaAxis, zAxis, occupancyAxis}}});
      registry.add(
        {"Tracks/PhiEtaGen",
         "; #varphi; #eta;",
         {HistType::kTHnSparseF, {phiAxis, etaAxis, occupancyAxis}}});
      registry.add({"Tracks/EtaZvtxGen_t",
                    "; #eta; Z_{vtx} (cm);",
                    {HistType::kTH2F, {etaAxis, zAxis}}});
      registry.add({"Tracks/PhiEtaGen_t",
                    "; #varphi; #eta;",
                    {HistType::kTH2F, {phiAxis, etaAxis}}});
      qaregistry.add({"Events/NotFoundEventZvtx",
                      "; #it{z}_{vtx} (cm)",
                      {HistType::kTH1F, {zAxis}}});
      qaregistry.add({"Events/ZvtxDiff",
                      "; Z_{rec} - Z_{gen} (cm)",
                      {HistType::kTH1F, {deltaZAxis}}});
      qaregistry.add({"Events/SplitMult",
                      "; N_{gen}; #it{z}_{vtx} (cm)",
                      {HistType::kTH2F, {multAxis, zAxis}}});
    }

    if (doprocessMCCentFT0C || doprocessMCwBestTracksCentFT0C ||
        doprocessMCCentFT0CVariant1 || doprocessMCwBestTracksCentFT0CVariant1 ||
        doprocessMCCentFT0M || doprocessMCwBestTracksCentFT0M ||
        doprocessMCCentNGlobal || doprocessMCwBestTracksCentNGlobal ||
        doprocessMCCentMFT || doprocessMCwBestTracksCentMFT) {
      registry.add({"Events/Centrality/EvtEffGen",
                    ";status;centrality;occupancy",
                    {HistType::kTHnSparseF,
                     {{3, 0.5, 3.5}, centralityAxis, occupancyAxis}}});
      auto heff = registry.get<THnSparse>(HIST("Events/Centrality/EvtEffGen"));
      heff->GetAxis(0)->SetBinLabel(1, "All reconstructed");
      heff->GetAxis(0)->SetBinLabel(2, "Selected reconstructed");
      heff->GetAxis(0)->SetBinLabel(3, "All generated");

      registry.add(
        {"Events/Centrality/NtrkZvtxGen_t",
         "; N_{trk}; Z_{vtx} (cm); centrality",
         {HistType::kTHnSparseF, {multAxis, zAxis, centralityAxis}}});
      registry.add({"Events/Centrality/NtrkZvtxGen",
                    "; N_{trk}; Z_{vtx} (cm); centrality; occupancy",
                    {HistType::kTHnSparseF,
                     {multAxis, zAxis, centralityAxis, occupancyAxis}}});
      registry.add({"Events/Centrality/hRecCent",
                    "; centrality; occupancy",
                    {HistType::kTH2F, {centralityAxis, occupancyAxis}}});
      registry.add(
        {"Events/Centrality/hRecZvtxCent",
         "; Z_{vtx} (cm); centrality; occupancy",
         {HistType::kTHnSparseF, {zAxis, centralityAxis, occupancyAxis}}});
      registry.add({"Tracks/Centrality/EtaZvtxGen",
                    "; #eta; Z_{vtx} (cm); centrality; occupancy",
                    {HistType::kTHnSparseF,
                     {etaAxis, zAxis, centralityAxis, occupancyAxis}}});
      registry.add({"Tracks/Centrality/PhiEtaGen",
                    "; #varphi; #eta; centrality; occupancy",
                    {HistType::kTHnSparseF,
                     {phiAxis, etaAxis, centralityAxis, occupancyAxis}}});
      registry.add({"Tracks/Centrality/EtaZvtxGen_t",
                    "; #eta; Z_{vtx} (cm); centrality",
                    {HistType::kTHnSparseF, {etaAxis, zAxis, centralityAxis}}});
      registry.add(
        {"Tracks/Centrality/PhiEtaGen_t",
         "; #varphi; #eta; centrality",
         {HistType::kTHnSparseF, {phiAxis, etaAxis, centralityAxis}}});
      qaregistry.add({"Events/Centrality/NotFoundEventZvtx",
                      "; #it{z}_{vtx} (cm); centrality",
                      {HistType::kTH2F, {zAxis, centralityAxis}}});
      qaregistry.add({"Events/Centrality/ZvtxDiff",
                      "; Z_{rec} - Z_{gen} (cm); centrality",
                      {HistType::kTH2F, {deltaZAxis, centralityAxis}}});
      qaregistry.add(
        {"Events/Centrality/SplitMult",
         "; N_{gen}; #it{z}_{vtx} (cm); centrality",
         {HistType::kTHnSparseF, {multAxis, zAxis, centralityAxis}}});
    }

    if (doprocessTrkEffIdxInlusive) {
      qaregistry.add({"Tracks/hPtPhiEtaZvtxEffGen",
                      "; p_{T} (GeV/c); #varphi; #eta; Z_{vtx} (cm); occupancy",
                      {HistType::kTHnSparseF,
                       {ptAxis, phiAxis, etaAxis, zAxis, occupancyAxis}}});
      qaregistry.add({"Tracks/hPtPhiEtaZvtxEffRec",
                      "; p_{T} (GeV/c); #varphi; #eta; Z_{vtx} (cm); occupancy",
                      {HistType::kTHnSparseF,
                       {ptAxis, phiAxis, etaAxis, zAxis, occupancyAxis}}});
      qaregistry.add(
        {"Tracks/hPhiEtaDuplicates",
         "; #varphi; #eta; occupancy",
         {HistType::kTHnSparseF, {phiAxis, etaAxis, occupancyAxis}}});
      qaregistry.add({"Tracks/hPtPhiEtaZvtxEffDuplicates",
                      "; p_{T} (GeV/c); #varphi; #eta; Z_{vtx} (cm); occupancy",
                      {HistType::kTHnSparseF,
                       {ptAxis, phiAxis, etaAxis, zAxis, occupancyAxis}}});
      qaregistry.add({"Tracks/hPtPhiEtaZvtxEffGenDuplicates",
                      "; p_{T} (GeV/c); #varphi; #eta; Z_{vtx} (cm); occupancy",
                      {HistType::kTHnSparseF,
                       {ptAxis, phiAxis, etaAxis, zAxis, occupancyAxis}}});
      qaregistry.add({"Tracks/NmftTrkPerPart",
                      "; #it{N}_{mft tracks per particle}; occupancy",
                      {HistType::kTH2F, {multAxis, occupancyAxis}}});
    }

    if (doprocessTrkEffIdxCentFT0C) {
      qaregistry.add(
        {"Tracks/Centrality/hPtPhiEtaZvtxEffGen",
         "; p_{T} (GeV/c); #varphi; #eta; Z_{vtx} (cm); centrality; "
         "occupancy",
         {HistType::kTHnSparseF,
          {ptAxis, phiAxis, etaAxis, zAxis, centralityAxis, occupancyAxis}}});
      qaregistry.add(
        {"Tracks/Centrality/hPtPhiEtaZvtxEffRec",
         "; p_{T} (GeV/c); #varphi; #eta; Z_{vtx} (cm); centrality; "
         "occupancy",
         {HistType::kTHnSparseF,
          {ptAxis, phiAxis, etaAxis, zAxis, centralityAxis, occupancyAxis}}});
      qaregistry.add({"Tracks/Centrality/hPhiEtaDuplicates",
                      "; #varphi; #eta; centrality; occupancy",
                      {HistType::kTHnSparseF,
                       {phiAxis, etaAxis, centralityAxis, occupancyAxis}}});
      qaregistry.add(
        {"Tracks/Centrality/hPtPhiEtaZvtxEffDuplicates",
         "; p_{T} (GeV/c); #varphi; #eta; Z_{vtx} (cm); centrality; "
         "occupancy",
         {HistType::kTHnSparseF,
          {ptAxis, phiAxis, etaAxis, zAxis, centralityAxis, occupancyAxis}}});
      qaregistry.add(
        {"Tracks/Centrality/hPtPhiEtaZvtxEffGenDuplicates",
         "; p_{T} (GeV/c); #varphi; #eta; Z_{vtx} (cm); centrality; "
         "occupancy",
         {HistType::kTHnSparseF,
          {ptAxis, phiAxis, etaAxis, zAxis, centralityAxis, occupancyAxis}}});
      qaregistry.add(
        {"Tracks/Centrality/NmftTrkPerPart",
         "; #it{N}_{mft tracks per particle}; centrality; occupancy",
         {HistType::kTHnSparseF, {multAxis, centralityAxis, occupancyAxis}}});
    }

    if (doprocessTrkEffBestInclusive) {
      qaregistry.add({"Tracks/hPtPhiEtaZvtxEffBestGen",
                      "; p_{T} (GeV/c); #varphi; #eta; Z_{vtx} (cm); occupancy",
                      {HistType::kTHnSparseF,
                       {ptAxis, phiAxis, etaAxis, zAxis, occupancyAxis}}});
      qaregistry.add({"Tracks/hPtPhiEtaZvtxEffBestRec",
                      "; p_{T} (GeV/c); #varphi; #eta; Z_{vtx} (cm); occupancy",
                      {HistType::kTHnSparseF,
                       {ptAxis, phiAxis, etaAxis, zAxis, occupancyAxis}}});
      qaregistry.add({"Tracks/hPtEffBestFakeRec",
                      " ; p_{T} (GeV/c); occupancy",
                      {HistType::kTH2F, {ptAxis, occupancyAxis}}});
    }

    if (doprocessTrkEffBestCentFT0C) {
      qaregistry.add(
        {"Tracks/Centrality/hPtPhiEtaZvtxEffBestGen",
         "; p_{T} (GeV/c); #varphi; #eta; Z_{vtx} (cm); centrality; "
         "occupancy",
         {HistType::kTHnSparseF,
          {ptAxis, phiAxis, etaAxis, zAxis, centralityAxis, occupancyAxis}}});
      qaregistry.add(
        {"Tracks/Centrality/hPtPhiEtaZvtxEffBestRec",
         "; p_{T} (GeV/c); #varphi; #eta; Z_{vtx} (cm); centrality; "
         "occupancy",
         {HistType::kTHnSparseF,
          {ptAxis, phiAxis, etaAxis, zAxis, centralityAxis, occupancyAxis}}});
      qaregistry.add(
        {"Tracks/Centrality/hPtEffBestFakeRec",
         "; p_{T} (GeV/c); centrality; occupancy",
         {HistType::kTHnSparseF, {ptAxis, centralityAxis, occupancyAxis}}});
    }

    if (doprocessMcQAInclusive) {
      qaregistry.add({"Events/hRecPerGenColls",
                      "; #it{N}_{reco collisions} / #it{N}_{gen collisions};",
                      {HistType::kTH2F, {{200, 0., 2.}, occupancyAxis}}});
      qaregistry.add({"Tracks/hNmftTrks",
                      "; #it{N}_{mft tracks};",
                      {HistType::kTH2F, {{200, -0.5, 200.}, occupancyAxis}}});
      qaregistry.add({"Tracks/hFracAmbiguousMftTrks",
                      "; #it{N}_{ambiguous tracks} / #it{N}_{tracks};",
                      {HistType::kTH2F, {{100, 0., 1.}, occupancyAxis}}});
    }

    if (doprocessMcQACentFT0C) {
      qaregistry.add(
        {"Events/Centrality/hRecPerGenColls",
         "; #it{N}_{reco collisions} / #it{N}_{gen collisions}; centrality",
         {HistType::kTHnSparseF,
          {{200, 0., 2.}, centralityAxis, occupancyAxis}}});
      qaregistry.add({"Tracks/Centrality/hNmftTrks",
                      "; #it{N}_{mft tracks}; centrality",
                      {HistType::kTHnSparseF,
                       {{200, -0.5, 200.}, centralityAxis, occupancyAxis}}});
      qaregistry.add(
        {"Tracks/Centrality/hFracAmbiguousMftTrks",
         "; #it{N}_{ambiguous tracks} / #it{N}_{tracks}; centrality",
         {HistType::kTHnSparseF,
          {{100, 0., 1.}, centralityAxis, occupancyAxis}}});
    }

    if (doprocessCheckAmbiguousMftTracksInclusive) {
      qaregistry.add({"Tracks/hMftTracksAmbDegree",
                      " ; N_{coll}^{comp}; occupancy",
                      {HistType::kTH2F, {{41, -0.5, 40.5}, occupancyAxis}}});
      qaregistry.add({"Tracks/hAmbTrackType",
                      " ; Ambiguous track type; occupancy",
                      {HistType::kTH2F, {{5, -0.5, 4.5}, occupancyAxis}}});
      qaregistry.add({"Tracks/histAmbZvtx",
                      "#it{z}_{vtx} of collisions associated to a "
                      "track;#it{z}_{vtx} (cm);",
                      {HistType::kTH1F, {zAxis}}});
    }

    if (doprocessCheckAmbiguousMftTracksCentFT0C) {
      qaregistry.add({"Tracks/Centrality/hMftTracksAmbDegree",
                      " ; N_{coll}^{comp}; occupancy",
                      {HistType::kTHnSparseF,
                       {{41, -0.5, 40.5}, centralityAxis, occupancyAxis}}});
      qaregistry.add({"Tracks/Centrality/hAmbTrackType",
                      " ; Ambiguous track type; occupancy",
                      {HistType::kTHnSparseF,
                       {{5, -0.5, 4.5}, centralityAxis, occupancyAxis}}});
    }

    if (doprocessEfficiencyInclusive) {
      qaregistry.add({"Tracks/hEffRec",
                      "; p_{T} (GeV/c); #varphi; #eta; Z_{vtx} (cm); occupancy",
                      {HistType::kTHnSparseF,
                       {ptAxis, phiAxis, etaAxis, zAxis, occupancyAxis}}});
      qaregistry.add({"Tracks/hEffFake",
                      "; p_{T} (GeV/c); #varphi; #eta; Z_{vtx} (cm); occupancy",
                      {HistType::kTHnSparseF,
                       {ptAxis, phiAxis, etaAxis, zAxis, occupancyAxis}}});
    }

    if (doprocessEfficiencyCentFT0C) {
      qaregistry.add({"Tracks/Centrality/hEffRec",
                      "; p_{T} (GeV/c); #varphi; #eta; Z_{vtx} (cm); occupancy",
                      {HistType::kTHnSparseF,
                       {ptAxis, phiAxis, etaAxis, zAxis, occupancyAxis}}});
      qaregistry.add({"Tracks/Centrality/hEffFake",
                      "; p_{T} (GeV/c); #varphi; #eta; Z_{vtx} (cm); occupancy",
                      {HistType::kTHnSparseF,
                       {ptAxis, phiAxis, etaAxis, zAxis, occupancyAxis}}});
    }
  }

  /// Filters - tracks
  Filter filtTrkEta = (aod::fwdtrack::eta < trackCuts.maxEta) &&
                      (aod::fwdtrack::eta > trackCuts.minEta);
  Filter filtATrackID = (aod::fwdtrack::bestCollisionId >= 0);
  Filter filtATrackDCA = (nabs(aod::fwdtrack::bestDCAXY) < trackCuts.maxDCAxy);

  /// Filters - mc particles
  Filter primaries = (aod::mcparticle::flags &
                      (uint8_t)o2::aod::mcparticle::enums::PhysicalPrimary) ==
                     (uint8_t)o2::aod::mcparticle::enums::PhysicalPrimary;

  /// Joined tables
  using FullBCs = soa::Join<aod::BCsWithTimestamps, aod::BcSels>;
  using CollBCs = soa::Join<aod::BCsWithTimestamps, aod::Run3MatchedToBCSparse>;
  using Colls = soa::Join<aod::Collisions, aod::EvSels>;
  using Coll = Colls::iterator;
  using CollsCentFT0C = soa::Join<aod::Collisions, aod::CentFT0Cs, aod::EvSels>;
  using CollsCentFT0CVariant1 =
    soa::Join<aod::Collisions, aod::CentFT0CVariant1s, aod::EvSels>;
  using CollsCentFT0M = soa::Join<aod::Collisions, aod::CentFT0Ms, aod::EvSels>;
  using CollsCentNGlobal =
    soa::Join<aod::Collisions, aod::CentNGlobals, aod::EvSels>;
  using CollsCentMFT = soa::Join<aod::Collisions, aod::CentMFTs, aod::EvSels>;
  using CollCentFT0C = CollsCentFT0C::iterator;
  using CollsGenCentFT0C = soa::Join<aod::McCollisionLabels, aod::Collisions,
                                     aod::CentFT0Cs, aod::EvSels>;
  using CollGenCent = CollsGenCentFT0C::iterator;
  using MFTTracksLabeled = soa::Join<aod::MFTTracks, aod::McMFTTrackLabels>;
  using MftTracksWColls = soa::Join<aod::MFTTracks, aod::MFTTrkCompColls>;

  /// Filtered tables
  using FiltMftTracks = soa::Filtered<aod::MFTTracks>;
  using FiltMcMftTracks = soa::Filtered<MFTTracksLabeled>;
  using FiltBestTracks = soa::Filtered<aod::BestCollisionsFwd>;
  using FiltParticles = soa::Filtered<aod::McParticles>;

  bool isIRSelected(CollBCs::iterator const& bc, bool fillHis = false)
  {
    double ir = (eventCuts.minIR >= 0 || eventCuts.maxIR >= 0)
                  ? rateFetcher.fetch(ccdb.service, bc.timestamp(),
                                      bc.runNumber(), "ZNC hadronic") *
                      1.e-3
                  : -1;
    if (eventCuts.minIR >= 0 && ir < eventCuts.minIR) {
      return false;
    }
    if (fillHis) {
      registry.fill(HIST("hEvtSel"), 9);
    }
    if (eventCuts.maxIR >= 0 && ir > eventCuts.maxIR) {
      return false;
    }
    if (fillHis) {
      registry.fill(HIST("hEvtSel"), 10);
    }
    return true;
  }

  template <typename T>
  bool isTrackSelected(const T& track)
  {
    if (track.eta() < trackCuts.minEta || track.eta() > trackCuts.maxEta)
      return false;
    if (trackCuts.requireCA && !track.isCA())
      return false;
    if (track.nClusters() < trackCuts.minNclusterMft)
      return false;
    if (track.pt() < trackCuts.minPt)
      return false;
    if (trackCuts.usephiCut) {
      float phi = track.phi();
      o2::math_utils::bringTo02Pi(phi);
      if (phi < trackCuts.minPhi || trackCuts.maxPhi < phi) {
        return false;
      }
      if ((phi < trackCuts.phiCut) ||
          ((phi > PI - trackCuts.phiCut) && (phi < PI + trackCuts.phiCut)) ||
          (phi > TwoPI - trackCuts.phiCut) ||
          ((phi > ((PIHalf - 0.1) * PI) - trackCuts.phiCut) &&
           (phi < ((PIHalf - 0.1) * PI) + trackCuts.phiCut)))
        return false;
    }
    return true;
  }

  template <typename C, bool fillHis = false, typename T>
  int countTracks(T const& tracks, float z, float c, float occ)
  {
    auto nTrk = 0;
    if (tracks.size() > 0) {
      for (auto const& track : tracks) {
        if (fillHis) {
          if constexpr (has_reco_cent<C>) {
            qaregistry.fill(HIST("Tracks/Centrality/Chi2Eta"), track.chi2(),
                            track.eta(), c, occ);
            qaregistry.fill(HIST("Tracks/Centrality/Chi2"), track.chi2(), c,
                            occ);
            qaregistry.fill(HIST("Tracks/Centrality/NclustersEta"),
                            track.nClusters(), track.eta(), c, occ);
          } else {
            qaregistry.fill(HIST("Tracks/Chi2Eta"), track.chi2(), track.eta(),
                            occ);
            qaregistry.fill(HIST("Tracks/Chi2"), track.chi2(), occ);
            qaregistry.fill(HIST("Tracks/NclustersEta"), track.nClusters(),
                            track.eta(), occ);
          }
        }
        if (!isTrackSelected(track)) {
          continue;
        }
        if (fillHis) {
          float phi = track.phi();
          o2::math_utils::bringTo02Pi(phi);
          if (phi < 0.f || TwoPI < phi) {
            continue;
          }
          if constexpr (has_reco_cent<C>) {
            registry.fill(HIST("Tracks/Centrality/EtaZvtx"), track.eta(), z, c,
                          occ);
            registry.fill(HIST("Tracks/Centrality/PhiEta"), phi, track.eta(), c,
                          occ);
          } else {
            registry.fill(HIST("Tracks/EtaZvtx"), track.eta(), z, occ);
            registry.fill(HIST("Tracks/PhiEta"), phi, track.eta(), occ);
          }
        }
        ++nTrk;
      }
    }
    if (fillHis) {
      if constexpr (has_reco_cent<C>) {
        qaregistry.fill(HIST("Tracks/Centrality/NchSel"), nTrk, c, occ);
      } else {
        qaregistry.fill(HIST("Tracks/NchSel"), nTrk, occ);
      }
    }
    return nTrk;
  }

  template <typename C, bool fillHis = false, typename T, typename B>
  int countBestTracks(T const& /*tracks*/, B const& besttracks, float z,
                      float c, float occ)
  {
    auto nATrk = 0;
    if (besttracks.size() > 0) {
      for (auto const& atrack : besttracks) {
        auto itrack = atrack.template mfttrack_as<T>();
        if (!isTrackSelected(itrack)) {
          continue;
        }
        if (fillHis) {
          float phi = itrack.phi();
          o2::math_utils::bringTo02Pi(phi);
          if (phi < 0.f || TwoPI < phi) {
            continue;
          }
          if constexpr (has_reco_cent<C>) {
            registry.fill(HIST("Tracks/Centrality/EtaZvtxBest"), itrack.eta(),
                          z, c, occ);
            registry.fill(HIST("Tracks/Centrality/PhiEtaBest"), phi,
                          itrack.eta(), c, occ);
            qaregistry.fill(HIST("Tracks/Centrality/DCAXYPt"), itrack.pt(),
                            atrack.bestDCAXY(), c, occ);
            qaregistry.fill(HIST("Tracks/Centrality/DCAXY"), atrack.bestDCAXY(),
                            c, occ);
            qaregistry.fill(HIST("Tracks/Centrality/NclustersEtaBest"),
                            itrack.nClusters(), itrack.eta(), c, occ);
            if (itrack.collisionId() != atrack.bestCollisionId()) {
              qaregistry.fill(HIST("Tracks/Centrality/ReTracksEtaZvtx"),
                              itrack.eta(), z, c, occ);
              qaregistry.fill(HIST("Tracks/Centrality/ReTracksPhiEta"), phi,
                              itrack.eta(), c, occ);
            }
            qaregistry.fill(HIST("Tracks/Centrality/TrackAmbDegree"),
                            atrack.ambDegree(), c, occ);
          } else {
            registry.fill(HIST("Tracks/EtaZvtxBest"), itrack.eta(), z, occ);
            registry.fill(HIST("Tracks/PhiEtaBest"), phi, itrack.eta(), occ);
            qaregistry.fill(HIST("Tracks/DCAXYPt"), itrack.pt(),
                            atrack.bestDCAXY(), occ);
            qaregistry.fill(HIST("Tracks/DCAXY"), atrack.bestDCAXY(), occ);
            qaregistry.fill(HIST("Tracks/NclustersEtaBest"), itrack.nClusters(),
                            itrack.eta(), occ);
            if (itrack.collisionId() != atrack.bestCollisionId()) {
              qaregistry.fill(HIST("Tracks/ReTracksEtaZvtx"), itrack.eta(), z,
                              occ);
              qaregistry.fill(HIST("Tracks/ReTracksPhiEta"), phi, itrack.eta(),
                              occ);
            }
            qaregistry.fill(HIST("Tracks/TrackAmbDegree"), atrack.ambDegree(),
                            occ);
          }
        }
        ++nATrk;
      }
    }
    return nATrk;
  }

  template <typename P>
  int countPart(P const& particles)
  {
    auto nCharged = 0;
    for (auto const& particle : particles) {
      if (!isChrgParticle(particle.pdgCode())) {
        continue;
      }
      nCharged++;
    }
    return nCharged;
  }

  template <typename C>
  float getOccupancy(C const& collision, uint occEstimator)
  {
    switch (occEstimator) {
      case OccupancyEst::TrkITS:
        return collision.trackOccupancyInTimeRange();
      case OccupancyEst::Ft0C:
        return collision.ft0cOccupancyInTimeRange();
      default:
        LOG(fatal) << "No valid occupancy estimator ";
        break;
    }
    return -1.f;
  }

  template <bool fillHis = false, typename C>
  bool isGoodEvent(C const& collision)
  {
    if constexpr (fillHis) {
      registry.fill(HIST("hEvtSel"), 0);
    }
    if (!collision.sel8()) {
      return false;
    }
    if constexpr (fillHis) {
      registry.fill(HIST("hEvtSel"), 1);
    }
    if (!collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return false;
    }
    if constexpr (fillHis) {
      registry.fill(HIST("hEvtSel"), 2);
    }
    if (!collision.selection_bit(aod::evsel::kNoSameBunchPileup)) {
      return false;
    }
    if constexpr (fillHis) {
      registry.fill(HIST("hEvtSel"), 3);
    }
    if (std::abs(collision.posZ()) >= eventCuts.maxZvtx) {
      return false;
    }
    if constexpr (fillHis) {
      registry.fill(HIST("hEvtSel"), 4);
    }
    if (eventCuts.requireNoCollInTimeRangeStd &&
        !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
      return false;
    }
    if constexpr (fillHis) {
      registry.fill(HIST("hEvtSel"), 5);
    }
    if (eventCuts.requireNoCollInTimeRangeNarrow &&
        !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeNarrow)) {
      return false;
    }
    if constexpr (fillHis) {
      registry.fill(HIST("hEvtSel"), 6);
    }
    if (eventCuts.minOccupancy >= 0 &&
        getOccupancy(collision, eventCuts.occupancyEstimator) <
          eventCuts.minOccupancy) {
      return false;
    }
    if constexpr (fillHis) {
      registry.fill(HIST("hEvtSel"), 7);
    }
    if (eventCuts.maxOccupancy >= 0 &&
        getOccupancy(collision, eventCuts.occupancyEstimator) >
          eventCuts.maxOccupancy) {
      return false;
    }
    if constexpr (fillHis) {
      registry.fill(HIST("hEvtSel"), 8);
    }
    return true;
  }

  /// @brief Selection of charged particles
  /// @return true: charged; false: not charged
  bool isChrgParticle(int code)
  {
    auto p = pdg->GetParticle(code);
    auto charge = 0.;
    if (p != nullptr) {
      charge = p->Charge();
    }
    return std::abs(charge) >= 3.;
  }

  template <bool isCent, typename P>
  void fillHistMC(P const& particles, float c, float occ, float zvtx,
                  bool const gtZeroColl)
  {
    for (auto const& particle : particles) {
      if (!isChrgParticle(particle.pdgCode())) {
        continue;
      }

      float phi = particle.phi();
      o2::math_utils::bringTo02Pi(phi);
      if (phi < 0.f || TwoPI < phi) {
        continue;
      }
      if constexpr (isCent) {
        registry.fill(HIST("Tracks/Centrality/EtaZvtxGen_t"), particle.eta(),
                      zvtx, c);
        registry.fill(HIST("Tracks/Centrality/PhiEtaGen_t"), phi,
                      particle.eta(), c);
      } else {
        registry.fill(HIST("Tracks/EtaZvtxGen_t"), particle.eta(), zvtx);
        registry.fill(HIST("Tracks/PhiEtaGen_t"), phi, particle.eta());
      }

      if (gtZeroColl) {
        float phi = particle.phi();
        o2::math_utils::bringTo02Pi(phi);
        if (phi < 0.f || TwoPI < phi) {
          continue;
        }
        if constexpr (isCent) {
          registry.fill(HIST("Tracks/Centrality/EtaZvtxGen"), particle.eta(),
                        zvtx, c, occ);
          registry.fill(HIST("Tracks/Centrality/PhiEtaGen"), phi,
                        particle.eta(), c, occ);
        } else {
          registry.fill(HIST("Tracks/EtaZvtxGen"), particle.eta(), zvtx, occ);
          registry.fill(HIST("Tracks/PhiEtaGen"), phi, particle.eta(), occ);
        }
      }
    }
  }

  /// @brief process function for general event statistics
  void processTagging(FullBCs const& bcs, CollsCentFT0C const& collisions)
  {
    std::vector<typename std::decay_t<decltype(collisions)>::iterator> cols;
    for (auto const& bc : bcs) {
      if ((bc.selection_bit(aod::evsel::kIsBBT0A) &&
           bc.selection_bit(aod::evsel::kIsBBT0C)) != 0) {
        registry.fill(HIST("hBcSel"), 0);
        cols.clear();
        for (auto const& collision : collisions) {
          if (collision.has_foundBC()) {
            if (collision.foundBCId() == bc.globalIndex()) {
              cols.emplace_back(collision);
            }
          } else if (collision.bcId() == bc.globalIndex()) {
            cols.emplace_back(collision);
          }
        }
        LOGP(debug, "BC {} has {} collisions", bc.globalBC(), cols.size());
        if (!cols.empty()) {
          registry.fill(HIST("hBcSel"), 1);
          if (cols.size() > 1) {
            registry.fill(HIST("hBcSel"), 2);
          }
        }
      }
    }
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processTagging, "Collect event sample stats",
                 true);

  /// @brief process function for counting tracks
  template <typename C>
  void processData(typename C::iterator const& collision,
                   FiltMftTracks const& tracks, CollBCs const& /*bcs*/)
  {
    auto occ = getOccupancy(collision, eventCuts.occupancyEstimator);
    float c = getRecoCent(collision);
    auto bc = collision.template foundBC_as<CollBCs>();
    double ir = rateFetcher.fetch(ccdb.service, bc.timestamp(), bc.runNumber(),
                                  "ZNC hadronic") *
                1.e-3;

    if constexpr (has_reco_cent<C>) {
      registry.fill(HIST("Events/Centrality/Selection"), 1., c, occ);
    } else {
      registry.fill(HIST("Events/Selection"), 1., occ);
    }

    if (!isGoodEvent<true>(collision)) {
      return;
    }
    if (!isIRSelected(bc, true)) {
      return;
    }

    auto z = collision.posZ();
    if constexpr (has_reco_cent<C>) {
      registry.fill(HIST("Events/Centrality/Selection"), 2., c, occ);
      qaregistry.fill(HIST("Events/Centrality/hZvtxCent"), z, c, occ);
      qaregistry.fill(HIST("Events/Centrality/hCent"), c, occ);
      qaregistry.fill(HIST("hCentOccIRate"), c, occ, ir);

    } else {
      qaregistry.fill(HIST("hOccIRate"), occ, ir);
      registry.fill(HIST("Events/Selection"), 2., occ);
    }

    auto nTrk = countTracks<C, true>(tracks, z, c, occ);

    if constexpr (has_reco_cent<C>) {
      registry.fill(HIST("Events/Centrality/NtrkZvtx"), nTrk, z, c, occ);
    } else {
      registry.fill(HIST("Events/NtrkZvtx"), nTrk, z, occ);
    }
  }

  /// @brief process function for counting tracks (based on BestCollisionsFwd
  /// table)
  template <typename C>
  void processDatawBestTracks(
    typename C::iterator const& collision, FiltMftTracks const& tracks,
    soa::SmallGroups<aod::BestCollisionsFwd> const& besttracks,
    CollBCs const& /*bcs*/)
  {
    auto occ = getOccupancy(collision, eventCuts.occupancyEstimator);
    float c = getRecoCent(collision);
    auto bc = collision.template foundBC_as<CollBCs>();
    double ir = rateFetcher.fetch(ccdb.service, bc.timestamp(), bc.runNumber(),
                                  "ZNC hadronic") *
                1.e-3;

    if constexpr (has_reco_cent<C>) {
      registry.fill(HIST("Events/Centrality/Selection"), 1., c, occ);
    } else {
      registry.fill(HIST("Events/Selection"), 1., occ);
    }

    if (!isGoodEvent<false>(collision)) {
      return;
    }
    if (!isIRSelected(bc, true)) {
      return;
    }

    auto z = collision.posZ();
    if constexpr (has_reco_cent<C>) {
      registry.fill(HIST("Events/Centrality/Selection"), 2., c, occ);
      qaregistry.fill(HIST("hCentOccIRate"), c, occ, ir);
    } else {
      registry.fill(HIST("Events/Selection"), 2., occ);
      qaregistry.fill(HIST("hOccIRate"), occ, ir);
    }

    auto nBestTrks = countBestTracks<C, true>(tracks, besttracks, z, c, occ);

    if constexpr (has_reco_cent<C>) {
      registry.fill(HIST("Events/Centrality/NtrkZvtxBest"), nBestTrks, z, c,
                    occ);
    } else {
      registry.fill(HIST("Events/NtrkZvtxBest"), nBestTrks, z, occ);
    }
  }

  void processDataInclusive(Colls::iterator const& collision,
                            FiltMftTracks const& tracks, CollBCs const& bcs)
  {
    processData<Colls>(collision, tracks, bcs);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processDataInclusive,
                 "Count tracks (inclusive)", false);

  void processDataCentFT0C(CollsCentFT0C::iterator const& collision,
                           FiltMftTracks const& tracks, CollBCs const& bcs)
  {
    processData<CollsCentFT0C>(collision, tracks, bcs);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processDataCentFT0C,
                 "Count tracks in FT0C centrality bins", false);

  void
    processDataCentFT0CVariant1(CollsCentFT0CVariant1::iterator const& collision,
                                FiltMftTracks const& tracks, CollBCs const& bcs)
  {
    processData<CollsCentFT0CVariant1>(collision, tracks, bcs);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processDataCentFT0CVariant1,
                 "Count tracks in FT0CVariant1 centrality bins", false);

  void processDataCentFT0M(CollsCentFT0M::iterator const& collision,
                           FiltMftTracks const& tracks, CollBCs const& bcs)
  {
    processData<CollsCentFT0M>(collision, tracks, bcs);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processDataCentFT0M,
                 "Count tracks in FT0M centrality bins", false);

  void processDataCentNGlobal(CollsCentNGlobal::iterator const& collision,
                              FiltMftTracks const& tracks, CollBCs const& bcs)
  {
    processData<CollsCentNGlobal>(collision, tracks, bcs);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processDataCentNGlobal,
                 "Count tracks in NGlobal centrality bins", false);

  void processDataCentMFT(CollsCentMFT::iterator const& collision,
                          FiltMftTracks const& tracks, CollBCs const& bcs)
  {
    processData<CollsCentMFT>(collision, tracks, bcs);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processDataCentMFT,
                 "Count tracks in MFT centrality bins", false);

  void processDatawBestTracksInclusive(
    Colls::iterator const& collision, FiltMftTracks const& tracks,
    soa::SmallGroups<aod::BestCollisionsFwd> const& besttracks,
    CollBCs const& bcs)
  {
    processDatawBestTracks<Colls>(collision, tracks, besttracks, bcs);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processDatawBestTracksInclusive,
                 "Count tracks based on BestCollisionsFwd table (inclusive)",
                 false);

  void processDatawBestTracksCentFT0C(
    CollsCentFT0C::iterator const& collision, FiltMftTracks const& tracks,
    soa::SmallGroups<aod::BestCollisionsFwd> const& besttracks,
    CollBCs const& bcs)
  {
    processDatawBestTracks<CollsCentFT0C>(collision, tracks, besttracks, bcs);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processDatawBestTracksCentFT0C,
                 "Count tracks in FT0C centrality bins based on BestCollisionsFwd table",
                 false);

  void processDatawBestTracksCentFT0CVariant1(
    CollsCentFT0CVariant1::iterator const& collision,
    FiltMftTracks const& tracks,
    soa::SmallGroups<aod::BestCollisionsFwd> const& besttracks,
    CollBCs const& bcs)
  {
    processDatawBestTracks<CollsCentFT0CVariant1>(collision, tracks, besttracks,
                                                  bcs);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processDatawBestTracksCentFT0CVariant1,
                 "Count tracks in FT0CVariant1 centrality bins based on "
                 "BestCollisionsFwd table",
                 false);

  void processDatawBestTracksCentFT0M(
    CollsCentFT0M::iterator const& collision, FiltMftTracks const& tracks,
    soa::SmallGroups<aod::BestCollisionsFwd> const& besttracks,
    CollBCs const& bcs)
  {
    processDatawBestTracks<CollsCentFT0M>(collision, tracks, besttracks, bcs);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processDatawBestTracksCentFT0M,
                 "Count tracks in FT0M centrality bins based on BestCollisionsFwd table",
                 false);

  void processDatawBestTracksCentNGlobal(
    CollsCentNGlobal::iterator const& collision, FiltMftTracks const& tracks,
    soa::SmallGroups<aod::BestCollisionsFwd> const& besttracks,
    CollBCs const& bcs)
  {
    processDatawBestTracks<CollsCentNGlobal>(collision, tracks, besttracks,
                                             bcs);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processDatawBestTracksCentNGlobal,
                 "Count tracks in NGlobal centrality bins based on "
                 "BestCollisionsFwd table",
                 false);

  void processDatawBestTracksCentMFT(
    CollsCentMFT::iterator const& collision, FiltMftTracks const& tracks,
    soa::SmallGroups<aod::BestCollisionsFwd> const& besttracks,
    CollBCs const& bcs)
  {
    processDatawBestTracks<CollsCentMFT>(collision, tracks, besttracks, bcs);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processDatawBestTracksCentMFT,
                 "Count tracks in MFT centrality bins based on BestCollisionsFwd table",
                 false);

  Preslice<FiltMcMftTracks> perCol = o2::aod::fwdtrack::collisionId;
  PresliceUnsorted<CollsGenCentFT0C> recColPerMcCol =
    aod::mccollisionlabel::mcCollisionId;
  Partition<FiltParticles> mcSample = nabs(aod::mcparticle::eta) < 1.0f;

  /// @brief process template function to run on MC gen
  template <typename MC, typename C>
  void processMC(
    typename MC::iterator const& mcCollision,
    soa::SmallGroups<soa::Join<C, aod::McCollisionLabels>> const& collisions,
    FiltParticles const& particles, FiltMcMftTracks const& tracks)
  {
    bool gtZeroColl = false;
    int gtOneColl = 0;

    float cgen = -1;
    if constexpr (has_reco_cent<C>) {
      float crec_min = 105.f;
      for (const auto& collision : collisions) {
        if (isGoodEvent<false>(collision)) {
          float c = getRecoCent(collision);
          if (c < crec_min) {
            crec_min = c;
          }
        }
      }
      if (cgen < 0)
        cgen = crec_min;
    }

    float occgen = -1.;
    for (const auto& collision : collisions) {
      if (isGoodEvent<false>(collision)) {
        float o = getOccupancy(collision, eventCuts.occupancyEstimator);
        if (o > occgen) {
          occgen = o;
        }
      }
    }

    for (auto const& collision : collisions) {
      float occrec = getOccupancy(collision, eventCuts.occupancyEstimator);
      float crec = getRecoCent(collision);

      if constexpr (has_reco_cent<C>) {
        registry.fill(HIST("Events/Centrality/EvtEffGen"), 1., crec, occrec);
      } else {
        registry.fill(HIST("Events/EvtEffGen"), 1., occrec);
      }

      if (isGoodEvent<false>(collision)) {
        gtZeroColl = true;
        ++gtOneColl;
        auto z = collision.posZ();

        if constexpr (has_reco_cent<C>) {
          registry.fill(HIST("Events/Centrality/EvtEffGen"), 2., crec, occrec);
          registry.fill(HIST("Events/Centrality/hRecCent"), crec, occrec);
          registry.fill(HIST("Events/Centrality/hRecZvtxCent"), z, crec,
                        occrec);
        } else {
          registry.fill(HIST("Events/EvtEffGen"), 2., occrec);
        }

        auto perColSample = tracks.sliceBy(perCol, collision.globalIndex());
        auto nTrkRec = countTracks<C, true>(perColSample, z, crec, occrec);

        if constexpr (has_reco_cent<C>) {
          qaregistry.fill(HIST("Events/Centrality/ZvtxDiff"),
                          collision.posZ() - mcCollision.posZ(), crec);
        } else {
          qaregistry.fill(HIST("Events/ZvtxDiff"),
                          collision.posZ() - mcCollision.posZ());
        }

        if (eventCuts.useZDiffCut) {
          if (std::abs(collision.posZ() - mcCollision.posZ()) >
              eventCuts.maxZvtxDiff) {
            continue;
          }
        }

        if constexpr (has_reco_cent<C>) {
          registry.fill(HIST("Events/Centrality/NtrkZvtxGen"), nTrkRec,
                        collision.posZ(), crec, occrec);
        } else {
          registry.fill(HIST("Events/NtrkZvtxGen"), nTrkRec, collision.posZ(),
                        occrec);
        }
      }
    }

    if constexpr (has_reco_cent<C>) {
      registry.fill(HIST("Events/Centrality/EvtEffGen"), 3., cgen, occgen);
    } else {
      registry.fill(HIST("Events/EvtEffGen"), 3., occgen);
    }

    auto perCollMCsample = mcSample->sliceByCached(
      aod::mcparticle::mcCollisionId, mcCollision.globalIndex(), cache);
    auto nchrg = countPart(perCollMCsample);
    auto zvtxMC = mcCollision.posZ();

    if (gtOneColl > 1) {
      if constexpr (has_reco_cent<C>) {
        qaregistry.fill(HIST("Events/Centrality/SplitMult"), nchrg, zvtxMC, cgen);
      } else {
        qaregistry.fill(HIST("Events/SplitMult"), nchrg, zvtxMC);
      }
    }

    auto nCharged = countPart(particles);
    if constexpr (has_reco_cent<C>) {
      registry.fill(HIST("Events/Centrality/NtrkZvtxGen_t"), nCharged, zvtxMC,
                    cgen);
    } else {
      registry.fill(HIST("Events/NtrkZvtxGen_t"), nCharged, zvtxMC);
    }

    fillHistMC<has_reco_cent<C>>(particles, cgen, occgen, zvtxMC, gtZeroColl);

    if (collisions.size() == 0) {
      if constexpr (has_reco_cent<C>) {
        qaregistry.fill(HIST("Events/Centrality/NotFoundEventZvtx"),
                        mcCollision.posZ(), cgen);
      } else {
        qaregistry.fill(HIST("Events/NotFoundEventZvtx"), mcCollision.posZ());
      }
    }
  }

  void processMCInclusive(
    aod::McCollisions::iterator const& mccollision,
    soa::SmallGroups<soa::Join<Colls, aod::McCollisionLabels>> const& collisions,
    FiltParticles const& particles, FiltMcMftTracks const& tracks)
  {
    processMC<aod::McCollisions, Colls>(mccollision, collisions, particles,
                                        tracks);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processMCInclusive,
                 "Count MC particles (inclusive)", false);

  void processMCCentFT0C(
    aod::McCollisions::iterator const& mccollision,
    soa::SmallGroups<soa::Join<CollsCentFT0C, aod::McCollisionLabels>> const& collisions,
    FiltParticles const& particles, FiltMcMftTracks const& tracks)
  {
    processMC<aod::McCollisions, CollsCentFT0C>(mccollision, collisions,
                                                particles, tracks);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processMCCentFT0C,
                 "Count MC particles in FT0C centrality bins", false);

  void processMCCentFT0CVariant1(
    aod::McCollisions::iterator const& mccollision,
    soa::SmallGroups<soa::Join<CollsCentFT0CVariant1,
                               aod::McCollisionLabels>> const& collisions,
    FiltParticles const& particles, FiltMcMftTracks const& tracks)
  {
    processMC<aod::McCollisions, CollsCentFT0CVariant1>(mccollision, collisions,
                                                        particles, tracks);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processMCCentFT0CVariant1,
                 "Count MC particles in FT0CVariant1 centrality bins", false);

  void processMCCentFT0M(
    aod::McCollisions::iterator const& mccollision,
    soa::SmallGroups<soa::Join<CollsCentFT0M, aod::McCollisionLabels>> const& collisions,
    FiltParticles const& particles, FiltMcMftTracks const& tracks)
  {
    processMC<aod::McCollisions, CollsCentFT0M>(mccollision, collisions,
                                                particles, tracks);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processMCCentFT0M,
                 "Count MC particles in FT0M centrality bins", false);

  void processMCCentNGlobal(
    aod::McCollisions::iterator const& mccollision,
    soa::SmallGroups<soa::Join<CollsCentNGlobal,
                               aod::McCollisionLabels>> const& collisions,
    FiltParticles const& particles, FiltMcMftTracks const& tracks)
  {
    processMC<aod::McCollisions, CollsCentNGlobal>(mccollision, collisions,
                                                   particles, tracks);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processMCCentNGlobal,
                 "Count MC particles in NGlobal centrality bins", false);

  void processMCCentMFT(
    aod::McCollisions::iterator const& mccollision,
    soa::SmallGroups<soa::Join<CollsCentMFT, aod::McCollisionLabels>> const& collisions,
    FiltParticles const& particles, FiltMcMftTracks const& tracks)
  {
    processMC<aod::McCollisions, CollsCentMFT>(mccollision, collisions,
                                               particles, tracks);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processMCCentMFT,
                 "Count MC particles in MFT centrality bins", false);

  PresliceUnsorted<aod::BestCollisionsFwd> perColU =
    aod::fwdtrack::bestCollisionId;

  /// @brief process template function to run on MC truth using
  /// aod::BestCollisionsFwd tracks
  template <typename MC, typename C>
  void processMCwBestTracks(
    typename MC::iterator const& mcCollision,
    soa::SmallGroups<soa::Join<C, aod::McCollisionLabels>> const& collisions,
    FiltParticles const& particles, FiltMcMftTracks const& tracks,
    FiltBestTracks const& besttracks)
  {
    bool gtZeroColl = false;
    float cgen = -1;
    if constexpr (has_reco_cent<C>) {
      float crec_min = 105.f;
      for (const auto& collision : collisions) {
        if (isGoodEvent<false>(collision)) {
          float c = getRecoCent(collision);
          if (c < crec_min) {
            crec_min = c;
          }
        }
      }
      if (cgen < 0)
        cgen = crec_min;
    }

    float occgen = -1.;
    for (const auto& collision : collisions) {
      if (isGoodEvent<false>(collision)) {
        float o = getOccupancy(collision, eventCuts.occupancyEstimator);
        if (o > occgen) {
          occgen = o;
        }
      }
    }

    for (auto const& collision : collisions) {
      auto occrec = getOccupancy(collision, eventCuts.occupancyEstimator);
      float crec = getRecoCent(collision);

      if constexpr (has_reco_cent<C>) {
        registry.fill(HIST("Events/Centrality/EvtEffGen"), 1., crec, occrec);
      } else {
        registry.fill(HIST("Events/EvtEffGen"), 1., occrec);
      }

      if (isGoodEvent<false>(collision)) {
        gtZeroColl = true;
        auto z = collision.posZ();

        if constexpr (has_reco_cent<C>) {
          registry.fill(HIST("Events/Centrality/EvtEffGen"), 2., crec, occrec);
        } else {
          registry.fill(HIST("Events/EvtEffGen"), 2., occrec);
        }

        auto perCollisionSample =
          tracks.sliceBy(perCol, collision.globalIndex());
        auto perCollisionASample =
          besttracks.sliceBy(perColU, collision.globalIndex());
        auto nTrkRec = countBestTracks<C, false>(
          perCollisionSample, perCollisionASample, z, crec,
          collision.trackOccupancyInTimeRange());

        if constexpr (has_reco_cent<C>) {
          registry.fill(HIST("Events/Centrality/NtrkZvtxGen"), nTrkRec, z,
                        crec, occrec);
        } else {
          registry.fill(HIST("Events/NtrkZvtxGen"), nTrkRec, z, occrec);
        }
      }
    }

    if constexpr (has_reco_cent<C>) {
      registry.fill(HIST("Events/Centrality/EvtEffGen"), 3., cgen, occgen);
    } else {
      registry.fill(HIST("Events/EvtEffGen"), 3., occgen);
    }

    auto zvtxMC = mcCollision.posZ();
    auto nCharged = countPart(particles);
    if constexpr (has_reco_cent<C>) {
      registry.fill(HIST("Events/Centrality/NtrkZvtxGen_t"), nCharged, zvtxMC,
                    cgen);
    } else {
      registry.fill(HIST("Events/NtrkZvtxGen_t"), nCharged, zvtxMC);
    }

    fillHistMC<has_reco_cent<C>>(particles, cgen, occgen, zvtxMC, gtZeroColl);

    if (collisions.size() == 0) {
      if constexpr (has_reco_cent<C>) {
        qaregistry.fill(HIST("Events/Centrality/NotFoundEventZvtx"),
                        mcCollision.posZ(), cgen);
      } else {
        qaregistry.fill(HIST("Events/NotFoundEventZvtx"), mcCollision.posZ());
      }
    }
  }

  void processMCwBestTracksInclusive(
    aod::McCollisions::iterator const& mccollision,
    soa::SmallGroups<soa::Join<Colls, aod::McCollisionLabels>> const& collisions,
    FiltParticles const& particles, FiltMcMftTracks const& tracks,
    FiltBestTracks const& besttracks)
  {
    processMCwBestTracks<aod::McCollisions, Colls>(
      mccollision, collisions, particles, tracks, besttracks);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processMCwBestTracksInclusive,
                 "Count MC particles using aod::BestCollisionsFwd (inclusive)",
                 false);

  void processMCwBestTracksCentFT0C(
    aod::McCollisions::iterator const& mccollision,
    soa::SmallGroups<soa::Join<CollsCentFT0C, aod::McCollisionLabels>> const& collisions,
    FiltParticles const& particles, FiltMcMftTracks const& tracks,
    FiltBestTracks const& besttracks)
  {
    processMCwBestTracks<aod::McCollisions, CollsCentFT0C>(
      mccollision, collisions, particles, tracks, besttracks);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processMCwBestTracksCentFT0C,
                 "Count MC particles in FT0C centrality bins using aod::BestCollisionsFwd",
                 false);

  void processMCwBestTracksCentFT0CVariant1(
    aod::McCollisions::iterator const& mccollision,
    soa::SmallGroups<soa::Join<CollsCentFT0CVariant1,
                               aod::McCollisionLabels>> const& collisions,
    FiltParticles const& particles, FiltMcMftTracks const& tracks,
    FiltBestTracks const& besttracks)
  {
    processMCwBestTracks<aod::McCollisions, CollsCentFT0CVariant1>(
      mccollision, collisions, particles, tracks, besttracks);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processMCwBestTracksCentFT0CVariant1,
                 "Count MC particles in FT0CVariant1 centrality bins using "
                 "aod::BestCollisionsFwd",
                 false);

  void processMCwBestTracksCentFT0M(
    aod::McCollisions::iterator const& mccollision,
    soa::SmallGroups<soa::Join<CollsCentFT0M, aod::McCollisionLabels>> const& collisions,
    FiltParticles const& particles, FiltMcMftTracks const& tracks,
    FiltBestTracks const& besttracks)
  {
    processMCwBestTracks<aod::McCollisions, CollsCentFT0M>(
      mccollision, collisions, particles, tracks, besttracks);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processMCwBestTracksCentFT0M,
                 "Count MC particles in FT0M centrality bins using aod::BestCollisionsFwd",
                 false);

  void processMCwBestTracksCentNGlobal(
    aod::McCollisions::iterator const& mccollision,
    soa::SmallGroups<soa::Join<CollsCentNGlobal,
                               aod::McCollisionLabels>> const& collisions,
    FiltParticles const& particles, FiltMcMftTracks const& tracks,
    FiltBestTracks const& besttracks)
  {
    processMCwBestTracks<aod::McCollisions, CollsCentNGlobal>(
      mccollision, collisions, particles, tracks, besttracks);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processMCwBestTracksCentNGlobal,
                 "Count MC particles in NGlobal centrality bins using "
                 "aod::BestCollisionsFwd",
                 false);

  void processMCwBestTracksCentMFT(
    aod::McCollisions::iterator const& mccollision,
    soa::SmallGroups<soa::Join<CollsCentMFT, aod::McCollisionLabels>> const& collisions,
    FiltParticles const& particles, FiltMcMftTracks const& tracks,
    FiltBestTracks const& besttracks)
  {
    processMCwBestTracks<aod::McCollisions, CollsCentMFT>(
      mccollision, collisions, particles, tracks, besttracks);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processMCwBestTracksCentMFT,
                 "Count MC particles in MFT centrality bins using aod::BestCollisionsFwd",
                 false);

  using ParticlesI = soa::Join<aod::McParticles, aod::ParticlesToMftTracks>;
  Partition<ParticlesI> primariesI =
    ((aod::mcparticle::flags &
      (uint8_t)o2::aod::mcparticle::enums::PhysicalPrimary) ==
     (uint8_t)o2::aod::mcparticle::enums::PhysicalPrimary);

  /// @brief process template function to calculate tracking efficiency (indexed
  /// as particle-to-MFT-tracks)
  template <typename C, typename MC>
  void processTrkEffIdx(
    typename soa::Join<C, aod::McCollisionLabels> const& collisions,
    MC const& /*mccollisions*/, ParticlesI const& /*particles*/,
    MFTTracksLabeled const& tracks)
  {
    for (auto const& collision : collisions) {
      if (!isGoodEvent<false>(collision)) {
        continue;
      }
      if (!collision.has_mcCollision()) {
        continue;
      }

      float crec = getRecoCent(collision);
      auto occrec = getOccupancy(collision, eventCuts.occupancyEstimator);
      auto mcCollision = collision.mcCollision();

      float cgen = -1;
      if constexpr (has_reco_cent<C>) {
        float crec_min = 105.f;
        for (const auto& collision : collisions) {
          if (isGoodEvent<false>(collision)) {
            float c = getRecoCent(collision);
            if (c < crec_min) {
              crec_min = c;
            }
          }
        }
        if (cgen < 0)
          cgen = crec_min;
      }

      float occgen = -1.;
      for (const auto& collision : collisions) {
        if (isGoodEvent<false>(collision)) {
          float o = getOccupancy(collision, eventCuts.occupancyEstimator);
          if (o > occgen) {
            occgen = o;
          }
        }
      }

      auto partsPerCol = primariesI->sliceByCached(
        aod::mcparticle::mcCollisionId, mcCollision.globalIndex(), cache);
      partsPerCol.bindExternalIndices(&tracks);

      for (auto const& particle : partsPerCol) {
        if (!isChrgParticle(particle.pdgCode())) {
          continue;
        }

        // MC gen
        if constexpr (has_reco_cent<C>) {
          qaregistry.fill(HIST("Tracks/Centrality/hPtPhiEtaZvtxEffGen"),
                          particle.pt(), particle.phi(), particle.eta(),
                          mcCollision.posZ(), cgen, occgen);
        } else {
          qaregistry.fill(HIST("Tracks/hPtPhiEtaZvtxEffGen"), particle.pt(),
                          particle.phi(), particle.eta(), mcCollision.posZ(),
                          occgen);
        }
        // MC rec
        if (particle.has_mfttracks()) {
          auto iscounted = false;
          auto ncnt = 0;
          auto relatedTracks =
            particle.template mfttracks_as<MFTTracksLabeled>();
          for (auto const& track : relatedTracks) {
            if (!isTrackSelected(track)) {
              continue;
            }
            ++ncnt;
            if constexpr (has_reco_cent<C>) {
              if (!iscounted) { // primaries
                qaregistry.fill(HIST("Tracks/Centrality/hPtPhiEtaZvtxEffRec"),
                                particle.pt(), particle.phi(), particle.eta(),
                                mcCollision.posZ(), crec, occrec);
                iscounted = true;
              }
              if (ncnt > 1) { // duplicates
                qaregistry.fill(HIST("Tracks/Centrality/hPhiEtaDuplicates"),
                                track.phi(), track.eta(), crec, occrec);
                qaregistry.fill(
                  HIST("Tracks/Centrality/hPtPhiEtaZvtxEffDuplicates"),
                  particle.pt(), particle.phi(), particle.eta(),
                  mcCollision.posZ(), crec, occrec);
              }
            } else {
              if (!iscounted) { // primaries
                qaregistry.fill(HIST("Tracks/hPtPhiEtaZvtxEffRec"),
                                particle.pt(), particle.phi(), particle.eta(),
                                mcCollision.posZ(), occrec);
                iscounted = true;
              }
              if (ncnt > 1) { // duplicates
                qaregistry.fill(HIST("Tracks/hPhiEtaDuplicates"), track.phi(),
                                track.eta(), occrec);
                qaregistry.fill(HIST("Tracks/hPtPhiEtaZvtxEffDuplicates"),
                                particle.pt(), particle.phi(), particle.eta(),
                                mcCollision.posZ(), occrec);
              }
            }
          }
          if constexpr (has_reco_cent<C>) {
            qaregistry.fill(HIST("Tracks/Centrality/NmftTrkPerPart"), ncnt,
                            crec, occrec);
          } else {
            qaregistry.fill(HIST("Tracks/NmftTrkPerPart"), ncnt, occrec);
          }
          if (relatedTracks.size() > 1) {
            if constexpr (has_reco_cent<C>) {
              qaregistry.fill(
                HIST("Tracks/Centrality/hPtPhiEtaZvtxEffGenDuplicates"),
                particle.pt(), particle.phi(), particle.eta(),
                mcCollision.posZ(), crec, occrec);
            } else {
              qaregistry.fill(HIST("Tracks/hPtPhiEtaZvtxEffGenDuplicates"),
                              particle.pt(), particle.phi(), particle.eta(),
                              mcCollision.posZ(), occrec);
            }
          }
        }
      }
    }
  }

  void processTrkEffIdxInlusive(
    soa::Join<Colls, aod::McCollisionLabels> const& collisions,
    aod::McCollisions const& mccollisions, ParticlesI const& particles,
    MFTTracksLabeled const& tracks)
  {
    processTrkEffIdx<Colls, aod::McCollisions>(collisions, mccollisions,
                                               particles, tracks);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processTrkEffIdxInlusive,
                 "Process tracking efficiency (inclusive, indexed)", false);

  void processTrkEffIdxCentFT0C(
    soa::Join<CollsCentFT0C, aod::McCollisionLabels> const& collisions,
    aod::McCollisions const& mccollisions, ParticlesI const& particles,
    MFTTracksLabeled const& tracks)
  {
    processTrkEffIdx<CollsCentFT0C, aod::McCollisions>(collisions, mccollisions,
                                                       particles, tracks);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processTrkEffIdxCentFT0C,
                 "Process tracking efficiency (in FT0C centrality bins, indexed)", false);

  /// @brief process function to calculate tracking efficiency (indexed) based
  /// on BestCollisionsFwd in FT0C bins
  template <typename C, typename MC>
  void processTrkEffBest(
    typename soa::Join<C, aod::McCollisionLabels>::iterator const& collision,
    MC const& /*mccollisions*/, FiltParticles const& particles,
    FiltMcMftTracks const& /*tracks*/,
    soa::SmallGroups<aod::BestCollisionsFwd> const& besttracks)
  {
    if (!isGoodEvent<false>(collision)) {
      return;
    }
    if (!collision.has_mcCollision()) {
      return;
    }

    float crec = getRecoCent(collision);
    auto occrec = getOccupancy(collision, eventCuts.occupancyEstimator);
    auto mcCollision = collision.mcCollision();
    auto partsPerCol = particles.sliceByCached(
      aod::mcparticle::mcCollisionId, mcCollision.globalIndex(), cache);

    for (auto const& particle : partsPerCol) {
      if (!isChrgParticle(particle.pdgCode())) {
        continue;
      }
      if constexpr (has_reco_cent<C>) {
        qaregistry.fill(HIST("Tracks/Centrality/hPtPhiEtaZvtxEffBestGen"),
                        particle.pt(), particle.phi(), particle.eta(),
                        mcCollision.posZ(), crec, occrec);
      } else {
        qaregistry.fill(HIST("Tracks/hPtPhiEtaZvtxEffBestGen"), particle.pt(),
                        particle.phi(), particle.eta(), mcCollision.posZ(),
                        occrec);
      }
    }

    for (auto const& track : besttracks) {
      auto itrack = track.mfttrack_as<FiltMcMftTracks>();
      if (!isTrackSelected(itrack)) {
        continue;
      }
      if (itrack.has_mcParticle()) {
        auto particle = itrack.mcParticle_as<FiltParticles>();
        if constexpr (has_reco_cent<C>) {
          qaregistry.fill(HIST("Tracks/Centrality/hPtPhiEtaZvtxEffBestRec"),
                          particle.pt(), itrack.phi(), itrack.eta(),
                          mcCollision.posZ(), crec, occrec);
        } else {
          qaregistry.fill(HIST("Tracks/hPtPhiEtaZvtxEffBestRec"), particle.pt(),
                          itrack.phi(), itrack.eta(), mcCollision.posZ(),
                          occrec);
        }
      } else {
        if constexpr (has_reco_cent<C>) {
          qaregistry.fill(HIST("Tracks/Centrality/hPtEffBestFakeRec"),
                          itrack.pt(), crec, occrec);
        } else {
          qaregistry.fill(HIST("Tracks/hPtEffBestFakeRec"), itrack.pt(),
                          occrec);
        }
      }
    }
  }

  void processTrkEffBestInclusive(
    soa::Join<Colls, aod::McCollisionLabels>::iterator const& collision,
    aod::McCollisions const& mccollisions, FiltParticles const& particles,
    FiltMcMftTracks const& tracks,
    soa::SmallGroups<aod::BestCollisionsFwd> const& besttracks)
  {
    processTrkEffBest<Colls, aod::McCollisions>(collision, mccollisions,
                                                particles, tracks, besttracks);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processTrkEffBestInclusive,
                 "Process tracking efficiency (inclusive, based on BestCollisionsFwd)",
                 false);

  void processTrkEffBestCentFT0C(
    soa::Join<CollsCentFT0C, aod::McCollisionLabels>::iterator const& collision,
    aod::McCollisions const& mccollisions, FiltParticles const& particles,
    FiltMcMftTracks const& tracks,
    soa::SmallGroups<aod::BestCollisionsFwd> const& besttracks)
  {
    processTrkEffBest<CollsCentFT0C, aod::McCollisions>(
      collision, mccollisions, particles, tracks, besttracks);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processTrkEffBestCentFT0C,
                 "Process tracking efficiency (in FT0 centrality bins, based "
                 "on BestCollisionsFwd)",
                 false);

  Preslice<FiltMcMftTracks> filtMcTrkperCol = o2::aod::fwdtrack::collisionId;

  /// @brief process function to calculate MC efficiency and fraction of fake
  /// tracks
  template <typename C, typename MC>
  void processEfficiency(
    typename soa::Join<C, aod::McCollisionLabels> const& collisions,
    MC const& /*mccollisions*/, FiltParticles const& /*particles*/,
    FiltMcMftTracks const& tracks)
  {
    for (auto const& collision : collisions) {
      if (!isGoodEvent<false>(collision)) {
        continue;
      }

      float crec = getRecoCent(collision);
      auto occrec = getOccupancy(collision, eventCuts.occupancyEstimator);
      auto mcCollision = collision.mcCollision();
      auto perColTrks =
        tracks.sliceBy(filtMcTrkperCol, collision.globalIndex());

      for (auto const& track : perColTrks) {
        if (!isTrackSelected(track)) {
          continue;
        }
        if (track.has_mcParticle()) {
          auto particle = track.template mcParticle_as<FiltParticles>();
          if constexpr (has_reco_cent<C>) {
            qaregistry.fill(HIST("Tracks/Centrality/hEffRec"), particle.pt(),
                            particle.phi(), particle.eta(), mcCollision.posZ(),
                            crec, occrec);
          } else {
            qaregistry.fill(HIST("Tracks/hEffRec"), particle.pt(),
                            particle.phi(), particle.eta(), mcCollision.posZ(),
                            crec, occrec);
          }
        } else {
          if constexpr (has_reco_cent<C>) {
            qaregistry.fill(HIST("Tracks/Centrality/hEffFake"), track.pt(),
                            track.phi(), track.eta(), mcCollision.posZ(), crec,
                            occrec);
          } else {
            qaregistry.fill(HIST("Tracks/hEffFake"), track.pt(), track.phi(),
                            track.eta(), mcCollision.posZ(), crec, occrec);
          }
        }
      }
    }
  }

  void processEfficiencyInclusive(
    soa::Join<Colls, aod::McCollisionLabels> const& collisions,
    aod::McCollisions const& mccollisions, FiltParticles const& particles,
    FiltMcMftTracks const& tracks)
  {
    processEfficiency<Colls, aod::McCollisions>(collisions, mccollisions,
                                                particles, tracks);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processEfficiencyInclusive,
                 "Process efficiencies (inclusive)", false);

  void processEfficiencyCentFT0C(
    soa::Join<CollsCentFT0C, aod::McCollisionLabels> const& collisions,
    aod::McCollisions const& mccollisions, FiltParticles const& particles,
    FiltMcMftTracks const& tracks)
  {
    processEfficiency<CollsCentFT0C, aod::McCollisions>(
      collisions, mccollisions, particles, tracks);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processEfficiencyCentFT0C,
                 "Process efficiencies in FT0C centrality bins", false);

  /// @brief process function to check ambiguous tracks
  template <typename C, typename allC>
  void processCheckAmbiguousMftTracks(typename C::iterator const& collision,
                                      allC const& allcollisions,
                                      MftTracksWColls const& tracks)
  {
    auto occ = getOccupancy(collision, eventCuts.occupancyEstimator);
    float c = getRecoCent(collision);

    bool ambTrk = false;
    int typeAmbTrk = 0;
    for (auto const& track : tracks) {
      if constexpr (has_reco_cent<C>) {
        qaregistry.fill(HIST("Tracks/Centrality/hMftTracksAmbDegree"),
                        track.compatibleCollIds().size(), c, occ);
      } else {
        qaregistry.fill(HIST("Tracks/hMftTracksAmbDegree"),
                        track.compatibleCollIds().size(), occ);
      }
      if (track.compatibleCollIds().size() > 0) {
        if (track.compatibleCollIds().size() == 1) {
          if (track.collisionId() != track.compatibleCollIds()[0]) {
            ambTrk = true;
            typeAmbTrk = 2;
          } else {
            typeAmbTrk = 1;
          }
        } else {
          ambTrk = true;
          typeAmbTrk = 3;

          for (const auto& collIdx : track.compatibleCollIds()) {
            auto ambColl = allcollisions.rawIteratorAt(collIdx);
            qaregistry.fill(HIST("Tracks/histAmbZvtx"), ambColl.posZ());
          }
        }
      }
    }

    if (ambTrk) {
      if constexpr (has_reco_cent<C>) {
        qaregistry.fill(HIST("Tracks/Centrality/hAmbTrackType"), typeAmbTrk, c,
                        occ);
      } else {
        qaregistry.fill(HIST("Tracks/hAmbTrackType"), typeAmbTrk, occ);
      }
    } else {
      if constexpr (has_reco_cent<C>) {
        qaregistry.fill(HIST("Tracks/Centrality/hAmbTrackType"), typeAmbTrk, c,
                        occ);
      } else {
        qaregistry.fill(HIST("Tracks/hAmbTrackType"), typeAmbTrk, occ);
      }
    }
  }

  void processCheckAmbiguousMftTracksInclusive(Colls::iterator const& collision,
                                               Colls const& allcollisions,
                                               MftTracksWColls const& track)
  {
    processCheckAmbiguousMftTracks<Colls, Colls>(collision, allcollisions,
                                                 track);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processCheckAmbiguousMftTracksInclusive,
                 "Process checks for Ambiguous MFT tracks (inclusive)", false);

  void processCheckAmbiguousMftTracksCentFT0C(
    CollsCentFT0C::iterator const& collision,
    CollsCentFT0C const& allcollisions, MftTracksWColls const& track)
  {
    processCheckAmbiguousMftTracks<CollsCentFT0C, CollsCentFT0C>(
      collision, allcollisions, track);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processCheckAmbiguousMftTracksCentFT0C,
                 "Process checks for Ambiguous MFT tracks (in FT0C centrality bins)",
                 false);

  Preslice<FiltMftTracks> filtTrkperCol = o2::aod::fwdtrack::collisionId;

  /// @brief process template function for MC QA checks
  template <typename C>
  void
    processMcQA(typename soa::Join<C, aod::McCollisionLabels> const& collisions,
                MFTTracksLabeled const& tracks,
                aod::AmbiguousMFTTracks const& atracks,
                aod::McCollisions const& mcCollisions,
                FiltParticles const& /*particles*/)
  {
    for (const auto& collision : collisions) {
      float crec = getRecoCent(collision);
      auto occrec = getOccupancy(collision, eventCuts.occupancyEstimator);

      if constexpr (has_reco_cent<C>) {
        qaregistry.fill(HIST("Events/Centrality/hRecPerGenColls"),
                        static_cast<float>(collisions.size()) /
                          mcCollisions.size(),
                        crec, occrec);
      } else {
        qaregistry.fill(HIST("Events/hRecPerGenColls"),
                        static_cast<float>(collisions.size()) /
                          mcCollisions.size(),
                        occrec);
      }

      if (!isGoodEvent<false>(collision)) {
        return;
      }

      auto trkPerColl = tracks.sliceBy(filtTrkperCol, collision.globalIndex());
      uint ntracks{0u}, nAtracks{0u};
      for (const auto& track : trkPerColl) {
        ntracks++;
        for (const auto& atrack : atracks) {
          if (atrack.mfttrackId() == track.globalIndex()) {
            nAtracks++;
            break;
          }
        }
      }
      if constexpr (has_reco_cent<C>) {
        qaregistry.fill(HIST("Tracks/Centrality/hNmftTrks"), ntracks, crec,
                        occrec);
        qaregistry.fill(HIST("Tracks/Centrality/hFracAmbiguousMftTrks"),
                        static_cast<float>(nAtracks) / ntracks, crec, occrec);
      } else {
        qaregistry.fill(HIST("Tracks/hNmftTrks"), ntracks, occrec);
        qaregistry.fill(HIST("Tracks/hFracAmbiguousMftTrks"),
                        static_cast<float>(nAtracks) / ntracks, occrec);
      }
    }
  }

  void processMcQAInclusive(
    soa::Join<Colls, aod::McCollisionLabels> const& collisions,
    MFTTracksLabeled const& tracks, aod::AmbiguousMFTTracks const& atracks,
    aod::McCollisions const& mcCollisions, FiltParticles const& particles)
  {
    processMcQA<Colls>(collisions, tracks, atracks, mcCollisions, particles);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processMcQAInclusive,
                 "Process MC QA checks (inclusive)", false);

  void processMcQACentFT0C(
    soa::Join<CollsCentFT0C, aod::McCollisionLabels> const& collisions,
    MFTTracksLabeled const& tracks, aod::AmbiguousMFTTracks const& atracks,
    aod::McCollisions const& mcCollisions, FiltParticles const& particles)
  {
    processMcQA<CollsCentFT0C>(collisions, tracks, atracks, mcCollisions,
                               particles);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processMcQACentFT0C,
                 "Process MC QA checks (in FT0 centrality bins)", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<DndetaMFTPbPb>(cfgc)};
}
