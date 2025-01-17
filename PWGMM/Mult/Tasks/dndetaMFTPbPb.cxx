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
/// \brief  Task for calculating dNdeta in Pb-Pb collisions using MFT detector
/// \author Gyula Bencedi, gyula.bencedi@cern.ch
/// \since  Nov 2024

#include <chrono>
#include <cmath>
#include <vector>

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/Configurable.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/RuntimeError.h"
#include "Framework/runDataProcessing.h"

#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "CommonConstants/MathConstants.h"

#include "MathUtils/Utils.h"
#include "ReconstructionDataFormats/GlobalTrackID.h"

#include <TPDGCode.h>

#include "Index.h"
#include "bestCollisionTable.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::track;
using namespace o2::aod::fwdtrack;

AxisSpec ptAxis = {1001, -0.005, 10.005};
AxisSpec multAxis = {701, -0.5, 700.5, "N_{trk}"};
AxisSpec zAxis = {60, -30., 30.};
AxisSpec deltaZAxis = {61, -6.1, 6.1};
AxisSpec dcaxyAxis = {500, -1, 50};
AxisSpec phiAxis = {629, 0, o2::constants::math::TwoPI, "Rad", "#phi"};
AxisSpec etaAxis = {20, -4., -2.};

struct DndetaMFTPbPb {
  SliceCache cache;

  // Histogram registry
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

  // analysis specific conf.
  Configurable<bool> usePhiCut{"usePhiCut", false, "use azimuthal angle cut"};
  Configurable<float> cfgPhiCut{"cfgPhiCut", 0.1f,
                                "Cut on azimuthal angle of MFT tracks"};

  // track selection conf.
  struct : ConfigurableGroup {
    Configurable<float> cfgPhiMin{"cfgPhiMin", 0.f, ""};
    Configurable<float> cfgPhiMax{"cfgPhiMax", 6.2832, ""};
    Configurable<float> cfgEtaMin{"cfgEtaMin", -3.6f, ""};
    Configurable<float> cfgEtaMax{"cfgEtaMax", -2.5f, ""};
    Configurable<int> cfgMinNclusterMft{"cfgMinNclusterMft", 5,
                                        "minimum number of MFT clusters"};
    Configurable<double> cfgPtMin{"cfgPtMin", 0.,
                                  "minimum pT of the MFT tracks"};
    Configurable<bool> cfgRequireCA{"cfgRequireCA", false,
                                    "Use Cellular Automaton track-finding algorithm"};
    Configurable<float> cfgDCAxyMax{"cfgDCAxyMax", 2.0f, "Cut on dcaXY"};
  } trkcuts;

  // event selection conf.
  Configurable<float> cfgCutZvtx{"cfgCutZvtx", 10.0f, "Cut on z-vtx"};
  Configurable<float> cfgCutCent{"cfgCutCent", 80.0f,
                                 "Cut on maximum centrality"};
  Configurable<bool> useZDiffCut{"useZDiffCut", false,
                                 "use Zvtx reco-mc diff. cut"};
  Configurable<float> maxZvtxDiff{"maxZvtxDiff", 1.0f,
                                  "max allowed Z vtx difference for reconstruced collisions (cm)"};
  Configurable<bool> requireNoCollInTimeRangeStd{"requireNoCollInTimeRangeStd", false,
                                                 "reject collisions corrupted by the cannibalism, with other collisions "
                                                 "within +/- 10 microseconds"};
  Configurable<bool> requireNoCollInTimeRangeNarrow{"requireNoCollInTimeRangeNarrow", false,
                                                    "reject collisions corrupted by the cannibalism, with other collisions "
                                                    "within +/- 10 microseconds"};
  ConfigurableAxis occupancyBins{"occupancyBins",
                                 {VARIABLE_WIDTH, 0.0f, 250.0f, 500.0f, 750.0f,
                                  1000.0f, 1500.0f, 2000.0f, 3000.0f, 4500.0f,
                                  6000.0f, 8000.0f, 10000.0f, 50000.0f},
                                 "Occupancy"};
  Configurable<float> minOccupancy{"minOccupancy", -1, "minimum occupancy from neighbouring collisions"};
  Configurable<float> maxOccupancy{"maxOccupancy", -1, "maximum occupancy from neighbouring collisions"};
  ConfigurableAxis centBins{"centBins",
                            {VARIABLE_WIDTH, 0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100},
                            ""};

  Service<o2::framework::O2DatabasePDG> pdg;

  /// @brief init function, definition of histograms
  void init(InitContext&)
  {
    if (static_cast<int>(doprocessDataInclusive) +
          static_cast<int>(doprocessDatawBestTracksInclusive) >
        1) {
      LOGP(fatal,
           "Either processDataInclusive OR "
           "processDatawBestTracksInclusive should be enabled!");
    }
    if (static_cast<int>(doprocessDataCent) +
          static_cast<int>(doprocessDatawBestTracksCent) >
        1) {
      LOGP(fatal,
           "Either processDataCent OR processDatawBestTracksCent should "
           "be enabled!");
    }
    if (static_cast<int>(doprocessMCInclusive) +
          static_cast<int>(doprocessMCwBestTracksInclusive) >
        1) {
      LOGP(fatal,
           "Either processMCInclusive OR processMCwBestTracksInclusive "
           "should be enabled!");
    }
    if (static_cast<int>(doprocessMCCent) +
          static_cast<int>(doprocessMCwBestTracksCent) >
        1) {
      LOGP(fatal,
           "Either processMCCent OR processMCwBestTracksCent should be "
           "enabled!");
    }

    auto hev = registry.add<TH1>("hEvtSel", "hEvtSel", HistType::kTH1F,
                                 {{10, -0.5f, +9.5f}});
    hev->GetXaxis()->SetBinLabel(1, "All collisions");
    hev->GetXaxis()->SetBinLabel(2, "Ev. sel.");
    hev->GetXaxis()->SetBinLabel(3, "kIsGoodZvtxFT0vsPV");
    hev->GetXaxis()->SetBinLabel(4, "NoSameBunchPileup");
    hev->GetXaxis()->SetBinLabel(5, "Z-vtx cut");
    hev->GetXaxis()->SetBinLabel(6, "kNoCollInTimeRangeStd");
    hev->GetXaxis()->SetBinLabel(7, "kNoCollInTimeRangeNarrow");
    hev->GetXaxis()->SetBinLabel(8, "Below min occup.");
    hev->GetXaxis()->SetBinLabel(9, "Above max occup.");

    auto hBcSel = registry.add<TH1>("hBcSel", "hBcSel", HistType::kTH1F,
                                    {{3, -0.5f, +2.5f}});
    hBcSel->GetXaxis()->SetBinLabel(1, "Good BCs");
    hBcSel->GetXaxis()->SetBinLabel(2, "BCs with collisions");
    hBcSel->GetXaxis()->SetBinLabel(3, "BCs with pile-up/splitting");

    AxisSpec centAxis = {centBins, "Centrality", "CentralityAxis"};
    AxisSpec occupancyAxis = {occupancyBins, "Occupancy", "occupancyAxis"};

    if (doprocessDataInclusive || doprocessDatawBestTracksInclusive) {
      registry.add({"Events/Selection",
                    ";status;events",
                    {HistType::kTH1F, {{2, 0.5, 2.5}}}});
      auto hstat = registry.get<TH1>(HIST("Events/Selection"));
      auto* x = hstat->GetXaxis();
      x->SetBinLabel(1, "All");
      x->SetBinLabel(2, "Selected");

      registry.add({"Events/NtrkZvtx",
                    "; N_{trk}; Z_{vtx} (cm);",
                    {HistType::kTH2F, {multAxis, zAxis}}});
      registry.add({"Tracks/EtaZvtx",
                    "; #eta; Z_{vtx} (cm);",
                    {HistType::kTH2F, {etaAxis, zAxis}}});
      registry.add({"Tracks/PhiEta",
                    "; #varphi; #eta;",
                    {HistType::kTH2F, {phiAxis, etaAxis}}});
      qaregistry.add({"Tracks/Chi2Eta",
                      "; #chi^{2}; #it{#eta};",
                      {HistType::kTH2F, {{600, 0, 20}, {100, -8, 8}}}});
      qaregistry.add(
        {"Tracks/Chi2", "; #chi^{2};", {HistType::kTH1F, {{600, 0, 20}}}});
      qaregistry.add({"Tracks/NclustersEta",
                      "; nClusters; #eta;",
                      {HistType::kTH2F, {{7, 4, 10}, {100, -8, 8}}}});
      qaregistry.add("Events/Occupancy", "; Z_{vtx} (cm); Occupancy",
                     HistType::kTH2F, {zAxis, occupancyAxis}, false);

      if (doprocessDatawBestTracksInclusive) {
        registry.add({"Events/NtrkZvtxBest",
                      "; N_{trk}; Z_{vtx} (cm);",
                      {HistType::kTH2F, {multAxis, zAxis}}});
        registry.add({"Tracks/EtaZvtxBest",
                      "; #eta; Z_{vtx} (cm);",
                      {HistType::kTH2F, {etaAxis, zAxis}}});
        registry.add({"Tracks/PhiEtaBest",
                      "; #varphi; #eta;",
                      {HistType::kTH2F, {phiAxis, etaAxis}}});
        qaregistry.add({"Tracks/NclustersEtaBest",
                        "; nClusters; #eta;",
                        {HistType::kTH2F, {{7, 4, 10}, {100, -8, 8}}}});
        qaregistry.add({"Tracks/DCAXYPt",
                        " ; p_{T} (GeV/c) ; DCA_{XY} (cm); centrality",
                        {HistType::kTH2F, {ptAxis, dcaxyAxis}}});
        qaregistry.add({"Tracks/DCAXY",
                        " ; DCA_{XY} (cm)",
                        {HistType::kTH1F, {dcaxyAxis}}});
        qaregistry.add({"Tracks/ReTracksEtaZvtx",
                        "; #eta; #it{z}_{vtx} (cm); tracks",
                        {HistType::kTH2F, {etaAxis, zAxis}}});
        qaregistry.add({"Tracks/ReTracksPhiEta",
                        "; #varphi; #eta; tracks",
                        {HistType::kTH2F, {phiAxis, etaAxis}}});
        qaregistry.add({"Tracks/TrackAmbDegree",
                        " ; N_{coll}^{comp}",
                        {HistType::kTH1F, {{51, -0.5, 50.5}}}});
      }
    }

    if (doprocessDataCent || doprocessDatawBestTracksCent) {
      registry.add({"Events/Centrality/Selection",
                    ";status;centrality;events",
                    {HistType::kTH2F, {{2, 0.5, 2.5}, centAxis}}});
      auto hstat = registry.get<TH2>(HIST("Events/Centrality/Selection"));
      auto* x = hstat->GetXaxis();
      x->SetBinLabel(1, "All");
      x->SetBinLabel(2, "Selected");

      registry.add({"Events/Centrality/NtrkZvtx",
                    "; N_{trk}; Z_{vtx} (cm); centrality",
                    {HistType::kTH3F, {multAxis, zAxis, centAxis}}});
      registry.add({"Tracks/Centrality/EtaZvtx",
                    "; #eta; Z_{vtx} (cm); centrality",
                    {HistType::kTH3F, {etaAxis, zAxis, centAxis}}});
      registry.add({"Tracks/Centrality/PhiEta",
                    "; #varphi; #eta; centrality",
                    {HistType::kTH3F, {phiAxis, etaAxis, centAxis}}});
      qaregistry.add({"Events/Centrality/hcentFT0C",
                      " ; cent FT0C",
                      {HistType::kTH1F, {centAxis}},
                      true});
      qaregistry.add(
        {"Tracks/Centrality/Chi2Eta",
         "; #chi^{2}; #it{#eta}; centrality",
         {HistType::kTH3F, {{600, 0, 20}, {100, -8, 8}, centAxis}}});
      qaregistry.add({"Tracks/Centrality/Chi2",
                      "; #chi^{2}; centrality",
                      {HistType::kTH2F, {{600, 0, 20}, centAxis}}});
      qaregistry.add({"Tracks/Centrality/NclustersEta",
                      "; nClusters; #eta; centrality",
                      {HistType::kTH3F, {{7, 4, 10}, {100, -8, 8}, centAxis}}});
      qaregistry.add("Events/Centrality/Occupancy",
                     "; Z_{vtx} (cm); centrality; Occupancy", HistType::kTH3F,
                     {zAxis, centAxis, occupancyAxis}, false);
      qaregistry.add("Tracks/Centrality/Occupancy", "dndeta occupancy",
                     HistType::kTHnSparseF,
                     {zAxis, centAxis, etaAxis, phiAxis, occupancyAxis}, false);

      if (doprocessDatawBestTracksCent) {
        registry.add({"Events/Centrality/NtrkZvtxBest",
                      "; N_{trk}; Z_{vtx} (cm); centrality",
                      {HistType::kTH3F, {multAxis, zAxis, centAxis}}});
        registry.add({"Tracks/Centrality/EtaZvtxBest",
                      "; #eta; Z_{vtx} (cm); centrality",
                      {HistType::kTH3F, {etaAxis, zAxis, centAxis}}});
        registry.add({"Tracks/Centrality/PhiEtaBest",
                      "; #varphi; #eta; centrality",
                      {HistType::kTH3F, {phiAxis, etaAxis, centAxis}}});
        qaregistry.add(
          {"Tracks/Centrality/NclustersEtaBest",
           "; nClusters; #eta; centrality",
           {HistType::kTH3F, {{7, 4, 10}, {100, -8, 8}, centAxis}}});
        qaregistry.add({"Tracks/Centrality/TrackAmbDegree",
                        " ; N_{coll}^{comp}",
                        {HistType::kTH2F, {{51, -0.5, 50.5}, centAxis}}});
        qaregistry.add({"Tracks/Centrality/DCAXY",
                        " ; DCA_{XY} (cm)",
                        {HistType::kTH2F, {dcaxyAxis, centAxis}}});
        qaregistry.add({"Tracks/Centrality/DCAXYPt",
                        " ; p_{T} (GeV/c) ; DCA_{XY} (cm); centrality",
                        {HistType::kTH3F, {ptAxis, dcaxyAxis, centAxis}}});
        qaregistry.add({"Tracks/Centrality/ReTracksEtaZvtx",
                        "; #eta; #it{z}_{vtx} (cm); tracks",
                        {HistType::kTH3F, {etaAxis, zAxis, centAxis}}});
        qaregistry.add({"Tracks/Centrality/ReTracksPhiEta",
                        "; #varphi; #eta; tracks",
                        {HistType::kTH3F, {phiAxis, etaAxis, centAxis}}});
        qaregistry.add("Events/Centrality/OccupancyBest",
                       "; Z_{vtx} (cm); centrality; Occupancy", HistType::kTH3F,
                       {zAxis, centAxis, occupancyAxis}, false);
        qaregistry.add("Tracks/Centrality/OccupancyBest", "dndeta occupancy",
                       HistType::kTHnSparseF,
                       {zAxis, centAxis, etaAxis, phiAxis, occupancyAxis},
                       false);
      }
    }

    if (doprocessMCInclusive || doprocessMCwBestTracksInclusive) {
      registry.add({"Events/EvtEffGen",
                    ";status;events",
                    {HistType::kTH1F, {{3, 0.5, 3.5}}}});
      auto heff = registry.get<TH1>(HIST("Events/EvtEffGen"));
      auto* x = heff->GetXaxis();
      x->SetBinLabel(1, "All reconstructed");
      x->SetBinLabel(2, "Selected reconstructed");
      x->SetBinLabel(3, "All generated");

      registry.add({"Events/NtrkZvtxGen_t",
                    "; N_{trk}; Z_{vtx} (cm);",
                    {HistType::kTH2F, {multAxis, zAxis}}});
      registry.add({"Events/NtrkZvtxGen",
                    "; N_{trk}; Z_{vtx} (cm);",
                    {HistType::kTH2F, {multAxis, zAxis}}});
      registry.add({"Tracks/EtaZvtxGen",
                    "; #eta; Z_{vtx} (cm);",
                    {HistType::kTH2F, {etaAxis, zAxis}}});
      registry.add({"Tracks/PhiEtaGen",
                    "; #varphi; #eta;",
                    {HistType::kTH2F, {phiAxis, etaAxis}}});
      registry.add({"Tracks/EtaZvtxGen_t",
                    "; #eta; Z_{vtx} (cm);",
                    {HistType::kTH2F, {etaAxis, zAxis}}});
      registry.add({"Tracks/PhiEtaGen_t",
                    "; #varphi; #eta;",
                    {HistType::kTH2F, {phiAxis, etaAxis}}});
      qaregistry.add({"Events/NotFoundEventZvtx",
                      " ; #it{z}_{vtx} (cm)",
                      {HistType::kTH1F, {zAxis}}});
      qaregistry.add({"Events/ZvtxDiff",
                      " ; Z_{rec} - Z_{gen} (cm)",
                      {HistType::kTH1F, {deltaZAxis}}});
      qaregistry.add(
        {"Events/SplitMult", " ; N_{gen}", {HistType::kTH1F, {multAxis}}});
    }

    if (doprocessMCCent || doprocessMCwBestTracksCent) {
      registry.add({"Events/Centrality/EvtEffGen",
                    ";status;events",
                    {HistType::kTH2F, {{3, 0.5, 3.5}, centAxis}}});
      auto heff = registry.get<TH2>(HIST("Events/Centrality/EvtEffGen"));
      auto* x = heff->GetXaxis();
      x->SetBinLabel(1, "All reconstructed");
      x->SetBinLabel(2, "Selected reconstructed");
      x->SetBinLabel(3, "All generated");

      registry.add({"Events/Centrality/NtrkZvtxGen_t",
                    "; N_{trk}; Z_{vtx} (cm);",
                    {HistType::kTH3F, {multAxis, zAxis, centAxis}}});
      registry.add({"Events/Centrality/NtrkZvtxGen",
                    "; N_{trk}; Z_{vtx} (cm);",
                    {HistType::kTH3F, {multAxis, zAxis, centAxis}}});
      registry.add({"Events/Centrality/hRecCent",
                    "Events/Centrality/hRecCent",
                    {HistType::kTH1F, {centAxis}}});
      registry.add({"Events/Centrality/hRecZvtxCent",
                    "Events/Centrality/hRecZvtxCent",
                    {HistType::kTH2F, {zAxis, centAxis}}});
      registry.add({"Tracks/Centrality/EtaZvtxGen",
                    "; #eta; Z_{vtx} (cm);",
                    {HistType::kTH3F, {etaAxis, zAxis, centAxis}}});
      registry.add({"Tracks/Centrality/PhiEtaGen",
                    "; #varphi; #eta;",
                    {HistType::kTH3F, {phiAxis, etaAxis, centAxis}}});
      registry.add({"Tracks/Centrality/EtaZvtxGen_t",
                    "; #eta; Z_{vtx} (cm);",
                    {HistType::kTH3F, {etaAxis, zAxis, centAxis}}});
      registry.add({"Tracks/Centrality/PhiEtaGen_t",
                    "; #varphi; #eta;",
                    {HistType::kTH3F, {phiAxis, etaAxis, centAxis}}});
      qaregistry.add({"Events/Centrality/NotFoundEventZvtx",
                      " ; #it{z}_{vtx} (cm)",
                      {HistType::kTH2F, {zAxis, centAxis}}});
      qaregistry.add({"Events/Centrality/ZvtxDiff",
                      " ; Z_{rec} - Z_{gen} (cm)",
                      {HistType::kTH2F, {deltaZAxis, centAxis}}});
      qaregistry.add({"Events/Centrality/SplitMult",
                      " ; N_{gen}",
                      {HistType::kTH2F, {multAxis, centAxis}}});
    }

    if (doprocessTrkEffIdxInlusive) {
      qaregistry.add(
        {"Tracks/hPtPhiEtaZvtxEffGen",
         "hPtPhiEtaZvtxEffGen",
         {HistType::kTHnSparseF, {ptAxis, phiAxis, etaAxis, zAxis}}});
      qaregistry.add(
        {"Tracks/hPtPhiEtaZvtxEffRec",
         "hPtPhiEtaZvtxEffRec",
         {HistType::kTHnSparseF, {ptAxis, phiAxis, etaAxis, zAxis}}});
      qaregistry.add({"Tracks/hPhiEtaDuplicates",
                      " ; p_{T} (GeV/c);",
                      {HistType::kTH2F, {phiAxis, etaAxis}}});
      qaregistry.add(
        {"Tracks/hPtPhiEtaZvtxEffDuplicates",
         "hPtPhiEtaZvtxEffDuplicates",
         {HistType::kTHnSparseF, {ptAxis, phiAxis, etaAxis, zAxis}}});
      qaregistry.add(
        {"Tracks/hPtPhiEtaZvtxEffGenDuplicates",
         "hPtPhiEtaZvtxEffGenDuplicates",
         {HistType::kTHnSparseF, {ptAxis, phiAxis, etaAxis, zAxis}}});
      qaregistry.add({"Tracks/NmftTrkPerPart",
                      "; #it{N}_{mft tracks per particle};",
                      {HistType::kTH1F, {{200, -0.5, 200.}}}});
    }

    if (doprocessTrkEffIdxCent) {
      qaregistry.add({"Tracks/Centrality/hPtPhiEtaZvtxEffGen",
                      "hPtPhiEtaZvtxEffGen",
                      {HistType::kTHnSparseF,
                       {ptAxis, phiAxis, etaAxis, zAxis, centAxis}}});
      qaregistry.add({"Tracks/Centrality/hPtPhiEtaZvtxEffRec",
                      "hPtPhiEtaZvtxEffRec",
                      {HistType::kTHnSparseF,
                       {ptAxis, phiAxis, etaAxis, zAxis, centAxis}}});
      qaregistry.add({"Tracks/Centrality/hPhiEtaDuplicates",
                      " ; p_{T} (GeV/c);",
                      {HistType::kTH3F, {phiAxis, etaAxis, centAxis}}});
      qaregistry.add({"Tracks/Centrality/hPtPhiEtaZvtxEffDuplicates",
                      "hPtPhiEtaZvtxEffDuplicates",
                      {HistType::kTHnSparseF,
                       {ptAxis, phiAxis, etaAxis, zAxis, centAxis}}});
      qaregistry.add({"Tracks/Centrality/hPtPhiEtaZvtxEffGenDuplicates",
                      "hPtPhiEtaZvtxEffGenDuplicates",
                      {HistType::kTHnSparseF,
                       {ptAxis, phiAxis, etaAxis, zAxis, centAxis}}});
      qaregistry.add({"Tracks/Centrality/NmftTrkPerPart",
                      "; #it{N}_{mft tracks per particle};",
                      {HistType::kTH2F, {{200, -0.5, 200.}, centAxis}}});
    }

    if (doprocessTrkEffBestInclusive) {
      qaregistry.add(
        {"Tracks/hPtPhiEtaZvtxEffBestGen",
         "hPtPhiEtaZvtxEffGen",
         {HistType::kTHnSparseF, {ptAxis, phiAxis, etaAxis, zAxis}}});
      qaregistry.add(
        {"Tracks/hPtPhiEtaZvtxEffBestRec",
         "hPtPhiEtaZvtxEffRec",
         {HistType::kTHnSparseF, {ptAxis, phiAxis, etaAxis, zAxis}}});
      qaregistry.add({"Tracks/hPtEffBestFakeRec",
                      " ; p_{T} (GeV/c);",
                      {HistType::kTH1F, {ptAxis}}});
    }

    if (doprocessTrkEffBestCent) {
      qaregistry.add({"Tracks/Centrality/hPtPhiEtaZvtxEffBestGen",
                      "hPtPhiEtaZvtxEffGen",
                      {HistType::kTHnSparseF,
                       {ptAxis, phiAxis, etaAxis, zAxis, centAxis}}});
      qaregistry.add({"Tracks/Centrality/hPtPhiEtaZvtxEffBestRec",
                      "hPtPhiEtaZvtxEffRec",
                      {HistType::kTHnSparseF,
                       {ptAxis, phiAxis, etaAxis, zAxis, centAxis}}});
      qaregistry.add({"Tracks/Centrality/hPtEffBestFakeRec",
                      " ; p_{T} (GeV/c);",
                      {HistType::kTH2F, {ptAxis, centAxis}}});
    }

    if (doprocessMcQAInclusive) {
      qaregistry.add({"Events/hRecPerGenColls",
                      "; #it{N}_{reco collisions} / #it{N}_{gen collisions};",
                      {HistType::kTH1F, {{200, 0., 2.}}}});
      qaregistry.add({"Tracks/hNmftTrks",
                      "; #it{N}_{mft tracks};",
                      {HistType::kTH1F, {{200, -0.5, 200.}}}});
      qaregistry.add({"Tracks/hFracAmbiguousMftTrks",
                      "; #it{N}_{ambiguous tracks} / #it{N}_{tracks};",
                      {HistType::kTH1F, {{100, 0., 1.}}}});
    }

    if (doprocessMcQACent) {
      qaregistry.add(
        {"Events/Centrality/hRecPerGenColls",
         "; #it{N}_{reco collisions} / #it{N}_{gen collisions}; centrality",
         {HistType::kTH2F, {{200, 0., 2.}, centAxis}}});
      qaregistry.add({"Tracks/Centrality/hNmftTrks",
                      "; #it{N}_{mft tracks}; centrality",
                      {HistType::kTH2F, {{200, -0.5, 200.}, centAxis}}});
      qaregistry.add(
        {"Tracks/Centrality/hFracAmbiguousMftTrks",
         "; #it{N}_{ambiguous tracks} / #it{N}_{tracks}; centrality",
         {HistType::kTH2F, {{100, 0., 1.}, centAxis}}});
    }
  }

  /// Filters - tracks
  Filter filtTrkEta = (aod::fwdtrack::eta < trkcuts.cfgEtaMax) &&
                      (aod::fwdtrack::eta > trkcuts.cfgEtaMin);
  Filter filtATrackID = (aod::fwdtrack::bestCollisionId >= 0);
  Filter filtATrackDCA = (nabs(aod::fwdtrack::bestDCAXY) < trkcuts.cfgDCAxyMax);

  /// Filters - mc particles
  Filter primaries = (aod::mcparticle::flags &
                      (uint8_t)o2::aod::mcparticle::enums::PhysicalPrimary) ==
                     (uint8_t)o2::aod::mcparticle::enums::PhysicalPrimary;

  /// Joined tables
  using FullBCs = soa::Join<aod::BCsWithTimestamps, aod::BcSels>;
  using CollBCs =
    soa::Join<aod::BCsWithTimestamps, aod::MatchedBCCollisionsSparseMulti>;
  using Colls = soa::Join<aod::Collisions, aod::EvSels>;
  using Coll = Colls::iterator;
  using CollsCent = soa::Join<aod::Collisions, aod::CentFT0Cs, aod::EvSels>;
  using CollCent = CollsCent::iterator;
  using CollsGenCent = soa::Join<aod::McCollisionLabels, aod::Collisions,
                                 aod::CentFT0Cs, aod::EvSels>;
  using CollGenCent = CollsGenCent::iterator;
  using MFTTracksLabeled = soa::Join<o2::aod::MFTTracks, aod::McMFTTrackLabels>;

  /// Filtered tables
  using FiltMftTracks = soa::Filtered<aod::MFTTracks>;
  using FiltMcMftTracks = soa::Filtered<MFTTracksLabeled>;
  using FiltBestTracks = soa::Filtered<aod::BestCollisionsFwd>;
  using FiltParticles = soa::Filtered<aod::McParticles>;

  template <typename T>
  bool isTrackSelected(const T& track)
  {
    if (track.eta() < trkcuts.cfgEtaMin || track.eta() > trkcuts.cfgEtaMax)
      return false;
    if (trkcuts.cfgRequireCA && !track.isCA())
      return false;
    if (track.nClusters() < trkcuts.cfgMinNclusterMft)
      return false;
    if (track.pt() < trkcuts.cfgPtMin)
      return false;
    if (usePhiCut) {
      float phi = track.phi();
      o2::math_utils::bringTo02Pi(phi);
      if (phi < trkcuts.cfgPhiMin || trkcuts.cfgPhiMax < phi) {
        return false;
      }
      if ((phi < cfgPhiCut) ||
          ((phi > o2::constants::math::PI - cfgPhiCut) &&
           (phi < o2::constants::math::PI + cfgPhiCut)) ||
          (phi > o2::constants::math::TwoPI - cfgPhiCut) ||
          ((phi >
            ((o2::constants::math::PIHalf - 0.1) * o2::constants::math::PI) -
              cfgPhiCut) &&
           (phi <
            ((o2::constants::math::PIHalf - 0.1) * o2::constants::math::PI) +
              cfgPhiCut)))
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
          if constexpr (C::template contains<aod::CentFT0Cs>()) {
            qaregistry.fill(HIST("Tracks/Centrality/Chi2Eta"), track.chi2(),
                            track.eta(), c);
            qaregistry.fill(HIST("Tracks/Centrality/Chi2"), track.chi2(), c);
            qaregistry.fill(HIST("Tracks/Centrality/NclustersEta"),
                            track.nClusters(), track.eta(), c);
          } else {
            qaregistry.fill(HIST("Tracks/Chi2Eta"), track.chi2(), track.eta());
            qaregistry.fill(HIST("Tracks/Chi2"), track.chi2());
            qaregistry.fill(HIST("Tracks/NclustersEta"), track.nClusters(),
                            track.eta());
          }
        }
        if (!isTrackSelected(track)) {
          continue;
        }
        if (fillHis) {
          float phi = track.phi();
          o2::math_utils::bringTo02Pi(phi);
          if (phi < 0.f || o2::constants::math::TwoPI < phi) {
            continue;
          }
          if constexpr (C::template contains<aod::CentFT0Cs>()) {
            registry.fill(HIST("Tracks/Centrality/EtaZvtx"), track.eta(), z, c);
            registry.fill(HIST("Tracks/Centrality/PhiEta"), phi, track.eta(),
                          c);
            qaregistry.fill(HIST("Tracks/Centrality/Occupancy"), z, c,
                            track.eta(), track.phi(), occ);
          } else {
            registry.fill(HIST("Tracks/EtaZvtx"), track.eta(), z);
            registry.fill(HIST("Tracks/PhiEta"), phi, track.eta());
          }
        }
        ++nTrk;
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
          if (phi < 0.f || o2::constants::math::TwoPI < phi) {
            continue;
          }
          if constexpr (C::template contains<aod::CentFT0Cs>()) {
            registry.fill(HIST("Tracks/Centrality/EtaZvtxBest"), itrack.eta(),
                          z, c);
            registry.fill(HIST("Tracks/Centrality/PhiEtaBest"), phi,
                          itrack.eta(), c);
            qaregistry.fill(HIST("Tracks/Centrality/OccupancyBest"), z, c,
                            itrack.eta(), itrack.phi(), occ);
            qaregistry.fill(HIST("Tracks/Centrality/DCAXYPt"), itrack.pt(),
                            atrack.bestDCAXY(), c);
            qaregistry.fill(HIST("Tracks/Centrality/DCAXY"), atrack.bestDCAXY(),
                            c);
            qaregistry.fill(HIST("Tracks/Centrality/NclustersEtaBest"),
                            itrack.nClusters(), itrack.eta(), c);
            if (itrack.collisionId() != atrack.bestCollisionId()) {
              qaregistry.fill(HIST("Tracks/Centrality/ReTracksEtaZvtx"),
                              itrack.eta(), z, c);
              qaregistry.fill(HIST("Tracks/Centrality/ReTracksPhiEta"), phi,
                              itrack.eta(), c);
            }
            qaregistry.fill(HIST("Tracks/Centrality/TrackAmbDegree"),
                            atrack.ambDegree(), c);
          } else {
            registry.fill(HIST("Tracks/EtaZvtxBest"), itrack.eta(), z);
            registry.fill(HIST("Tracks/PhiEtaBest"), phi, itrack.eta());
            qaregistry.fill(HIST("Tracks/DCAXYPt"), itrack.pt(),
                            atrack.bestDCAXY());
            qaregistry.fill(HIST("Tracks/DCAXY"), atrack.bestDCAXY());
            qaregistry.fill(HIST("Tracks/NclustersEtaBest"), itrack.nClusters(),
                            itrack.eta());
            if (itrack.collisionId() != atrack.bestCollisionId()) {
              qaregistry.fill(HIST("Tracks/ReTracksEtaZvtx"), itrack.eta(), z);
              qaregistry.fill(HIST("Tracks/ReTracksPhiEta"), phi, itrack.eta());
            }
            qaregistry.fill(HIST("Tracks/TrackAmbDegree"), atrack.ambDegree());
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
    if (std::abs(collision.posZ()) >= cfgCutZvtx) {
      return false;
    }
    if constexpr (fillHis) {
      registry.fill(HIST("hEvtSel"), 4);
    }
    if (requireNoCollInTimeRangeStd &&
        !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
      return false;
    }
    if constexpr (fillHis) {
      registry.fill(HIST("hEvtSel"), 5);
    }
    if (requireNoCollInTimeRangeNarrow &&
        !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeNarrow)) {
      return false;
    }
    if constexpr (fillHis) {
      registry.fill(HIST("hEvtSel"), 6);
    }
    if (minOccupancy > 0 &&
        collision.trackOccupancyInTimeRange() < minOccupancy) {
      return false;
    }
    if constexpr (fillHis) {
      registry.fill(HIST("hEvtSel"), 7);
    }
    if (maxOccupancy > 0 &&
        collision.trackOccupancyInTimeRange() > maxOccupancy) {
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
  void fillHistMC(P const& particles, float cent, float zvtx,
                  bool const atLeastOne)
  {
    for (auto const& particle : particles) {
      if (!isChrgParticle(particle.pdgCode())) {
        continue;
      }

      float phi = particle.phi();
      o2::math_utils::bringTo02Pi(phi);
      if (phi < 0.f || o2::constants::math::TwoPI < phi) {
        continue;
      }
      if constexpr (isCent) {
        registry.fill(HIST("Tracks/Centrality/EtaZvtxGen_t"), particle.eta(),
                      zvtx, cent);
        registry.fill(HIST("Tracks/Centrality/PhiEtaGen_t"), phi,
                      particle.eta(), cent);
      } else {
        registry.fill(HIST("Tracks/EtaZvtxGen_t"), particle.eta(), zvtx);
        registry.fill(HIST("Tracks/PhiEtaGen_t"), phi, particle.eta());
      }

      if (atLeastOne) {
        float phi = particle.phi();
        o2::math_utils::bringTo02Pi(phi);
        if (phi < 0.f || o2::constants::math::TwoPI < phi) {
          continue;
        }
        if constexpr (isCent) {
          registry.fill(HIST("Tracks/Centrality/EtaZvtxGen"), particle.eta(),
                        zvtx, cent);
          registry.fill(HIST("Tracks/Centrality/PhiEtaGen"), phi,
                        particle.eta(), cent);
        } else {
          registry.fill(HIST("Tracks/EtaZvtxGen"), particle.eta(), zvtx);
          registry.fill(HIST("Tracks/PhiEtaGen"), phi, particle.eta());
        }
      }
    }
  }

  /// @brief process fnc. for general event statistics
  void processTagging(FullBCs const& bcs, CollsCent const& collisions)
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

  PROCESS_SWITCH(DndetaMFTPbPb, processTagging,
                 "Collect event sample stats", true);

  template <typename C>
  void processData(typename C::iterator const& collision,
                   FiltMftTracks const& tracks)
  {
    float c = -1;
    if constexpr (C::template contains<aod::CentFT0Cs>()) {
      c = collision.centFT0C();
      registry.fill(HIST("Events/Centrality/Selection"), 1., c);
    } else {
      registry.fill(HIST("Events/Selection"), 1.);
    }
    if (!isGoodEvent<true>(collision)) {
      return;
    }
    auto z = collision.posZ();
    auto occ = collision.trackOccupancyInTimeRange();
    if constexpr (C::template contains<aod::CentFT0Cs>()) {
      registry.fill(HIST("Events/Centrality/Selection"), 2., c);
      qaregistry.fill(HIST("Events/Centrality/Occupancy"), z, c, occ);
      qaregistry.fill(HIST("Events/Centrality/hcentFT0C"), c);
    } else {
      registry.fill(HIST("Events/Selection"), 2.);
    }

    auto nTrk = countTracks<C, true>(
      tracks, z, c, occ); //!@note here we obtain eta-z and phi-eta
    if constexpr (C::template contains<aod::CentFT0Cs>()) {
      registry.fill(HIST("Events/Centrality/NtrkZvtx"), nTrk, z, c);
    } else {
      registry.fill(HIST("Events/NtrkZvtx"), nTrk, z);
    }
  }

  template <typename C>
  void processDatawBestTracks(
    typename C::iterator const& collision, FiltMftTracks const& tracks,
    soa::SmallGroups<aod::BestCollisionsFwd> const& besttracks)
  {
    float c = -1;
    if constexpr (C::template contains<aod::CentFT0Cs>()) {
      c = collision.centFT0C();
      registry.fill(HIST("Events/Centrality/Selection"), 1., c);
    } else {
      registry.fill(HIST("Events/Selection"), 1.);
    }
    if (!isGoodEvent<false>(collision)) {
      return;
    }
    auto z = collision.posZ();
    auto occ = collision.trackOccupancyInTimeRange();
    if constexpr (C::template contains<aod::CentFT0Cs>()) {
      registry.fill(HIST("Events/Centrality/Selection"), 2., c);
      qaregistry.fill(HIST("Events/Centrality/OccupancyBest"), z, c, occ);
    } else {
      registry.fill(HIST("Events/Selection"), 2.);
    }

    auto nBestTrks =
      countBestTracks<C, true>(tracks, besttracks, z, c,
                               occ); //!@note here we obtain eta-z and phi-eta
    if constexpr (C::template contains<aod::CentFT0Cs>()) {
      registry.fill(HIST("Events/Centrality/NtrkZvtxBest"), nBestTrks, z, c);
    } else {
      registry.fill(HIST("Events/NtrkZvtxBest"), nBestTrks, z);
    }
  }

  /// @brief process fnc. to run on DATA and REC MC w/o centrality selection
  void processDataInclusive(Colls::iterator const& collision,
                            FiltMftTracks const& tracks)
  {
    processData<Colls>(collision, tracks);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processDataInclusive, "Count tracks",
                 false);

  /// @brief process fnc. to run on DATA and REC MC w/ FT0C centrality selection
  void processDataCent(CollsCent::iterator const& collision,
                       FiltMftTracks const& tracks)
  {
    processData<CollsCent>(collision, tracks);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processDataCent,
                 "Count tracks in FT0C bins", false);

  /// @brief process fnc. to run on DATA and REC MC based on BestCollisionsFwd
  /// table w/o centrality selection
  void processDatawBestTracksInclusive(
    Colls::iterator const& collision, FiltMftTracks const& tracks,
    soa::SmallGroups<aod::BestCollisionsFwd> const& besttracks)
  {
    processDatawBestTracks<Colls>(collision, tracks, besttracks);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processDatawBestTracksInclusive,
                 "Count tracks based on BestCollisionsFwd table", false);

  /// @brief process fnc. to run on DATA and REC MC based on BestCollisionsFwd
  /// table w/ FT0C centrality selection
  void processDatawBestTracksCent(
    CollsCent::iterator const& collision, FiltMftTracks const& tracks,
    soa::SmallGroups<aod::BestCollisionsFwd> const& besttracks)
  {
    processDatawBestTracks<CollsCent>(collision, tracks, besttracks);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processDatawBestTracksCent,
                 "Count tracks in FT0C bins based on BestCollisionsFwd table",
                 false);

  Preslice<FiltMcMftTracks> perCol = o2::aod::fwdtrack::collisionId;
  Partition<FiltParticles> mcSample = nabs(aod::mcparticle::eta) < 1.0f;

  /// @brief process template function to run on MC truth
  /// @param cols subscribe to the collisions
  /// @param parts subscribe to filtered MC particle table
  template <typename MC, typename C>
  void processMC(
    typename MC::iterator const& mcCollision,
    soa::SmallGroups<soa::Join<C, aod::McCollisionLabels>> const& collisions,
    FiltParticles const& particles, FiltMcMftTracks const& tracks)
  {
    float cgen = -1;
    bool atLeastOne = false;
    int moreThanOne = 0;
    for (auto const& collision : collisions) {
      float crec = -1;
      if constexpr (C::template contains<aod::CentFT0Cs>()) {
        crec = collision.centFT0C();
        registry.fill(HIST("Events/Centrality/EvtEffGen"), 1., crec);
      } else {
        registry.fill(HIST("Events/EvtEffGen"), 1.);
      }

      if (isGoodEvent<false>(collision)) {
        if constexpr (C::template contains<aod::CentFT0Cs>()) {
          if (!atLeastOne) {
            cgen = crec;
          }
        }
        atLeastOne = true;
        ++moreThanOne;
        auto z = collision.posZ();

        if constexpr (C::template contains<aod::CentFT0Cs>()) {
          registry.fill(HIST("Events/Centrality/EvtEffGen"), 2., crec);
          registry.fill(HIST("Events/Centrality/hRecCent"), crec);
          registry.fill(HIST("Events/Centrality/hRecZvtxCent"), z, crec);
        } else {
          registry.fill(HIST("Events/EvtEffGen"), 2.);
        }

        auto perCollisionSample =
          tracks.sliceBy(perCol, collision.globalIndex());
        auto nTrkRec =
          countTracks<C, true>(perCollisionSample, z, crec,
                               collision.trackOccupancyInTimeRange());

        if constexpr (C::template contains<aod::CentFT0Cs>()) {
          qaregistry.fill(HIST("Events/Centrality/ZvtxDiff"),
                          collision.posZ() - mcCollision.posZ(), crec);
        } else {
          qaregistry.fill(HIST("Events/ZvtxDiff"),
                          collision.posZ() - mcCollision.posZ());
        }

        if (useZDiffCut) {
          if (std::abs(collision.posZ() - mcCollision.posZ()) > maxZvtxDiff) {
            continue;
          }
        }

        if constexpr (C::template contains<aod::CentFT0Cs>()) {
          registry.fill(HIST("Events/Centrality/NtrkZvtxGen"), nTrkRec,
                        collision.posZ(), crec);
        } else {
          registry.fill(HIST("Events/NtrkZvtxGen"), nTrkRec, collision.posZ());
        }
      }
    }

    if constexpr (C::template contains<aod::CentFT0Cs>()) {
      registry.fill(HIST("Events/Centrality/EvtEffGen"), 3., cgen);
    } else {
      registry.fill(HIST("Events/EvtEffGen"), 3.);
    }

    auto perCollMCsample = mcSample->sliceByCached(
      aod::mcparticle::mcCollisionId, mcCollision.globalIndex(), cache);
    auto nchrg = countPart(perCollMCsample);
    if (moreThanOne > 1) {
      if constexpr (C::template contains<aod::CentFT0Cs>()) {
        qaregistry.fill(HIST("Events/Centrality/SplitMult"), nchrg, cgen);
      } else {
        qaregistry.fill(HIST("Events/SplitMult"), nchrg);
      }
    }

    auto zvtxMC = mcCollision.posZ();
    auto nCharged = countPart(particles);
    if constexpr (C::template contains<aod::CentFT0Cs>()) {
      registry.fill(HIST("Events/Centrality/NtrkZvtxGen_t"), nCharged, zvtxMC,
                    cgen);
    } else {
      registry.fill(HIST("Events/NtrkZvtxGen_t"), nCharged, zvtxMC);
    }

    fillHistMC<C::template contains<aod::CentFT0Cs>()>(particles, cgen,
                                                       zvtxMC, atLeastOne);

    if (collisions.size() == 0) {
      if constexpr (C::template contains<aod::CentFT0Cs>()) {
        qaregistry.fill(HIST("Events/Centrality/NotFoundEventZvtx"),
                        mcCollision.posZ(), cgen);
      } else {
        qaregistry.fill(HIST("Events/NotFoundEventZvtx"), mcCollision.posZ());
      }
    }
  }

  /// @brief process fnc. to run on MC w/o centrality selection
  void processMCInclusive(
    aod::McCollisions::iterator const& mccollision,
    soa::SmallGroups<soa::Join<Colls, aod::McCollisionLabels>> const& collisions,
    FiltParticles const& particles, FiltMcMftTracks const& tracks)
  {
    processMC<aod::McCollisions, Colls>(mccollision, collisions, particles,
                                        tracks);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processMCInclusive,
                 "Count MC particles", false);

  /// @brief process fnc. to run on MC w FT0C centrality selection
  void processMCCent(
    aod::McCollisions::iterator const& mccollision,
    soa::SmallGroups<soa::Join<CollsCent, aod::McCollisionLabels>> const& collisions,
    FiltParticles const& particles, FiltMcMftTracks const& tracks)
  {
    processMC<aod::McCollisions, CollsCent>(mccollision, collisions, particles,
                                            tracks);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processMCCent,
                 "Count MC particles in FT0C bins", false);

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
    float cgen = -1;
    bool atLeastOne = false;
    // int moreThanOne = 0;
    for (auto const& collision : collisions) {
      float crec = -1;
      if constexpr (C::template contains<aod::CentFT0Cs>()) {
        crec = collision.centFT0C();
        registry.fill(HIST("Events/Centrality/EvtEffGen"), 1., crec);
      } else {
        registry.fill(HIST("Events/EvtEffGen"), 1.);
      }

      if (isGoodEvent<false>(collision)) {
        if constexpr (C::template contains<aod::CentFT0Cs>()) {
          if (!atLeastOne) {
            cgen = crec;
          }
        }
        atLeastOne = true;
        // ++moreThanOne;
        auto z = collision.posZ();

        if constexpr (C::template contains<aod::CentFT0Cs>()) {
          registry.fill(HIST("Events/Centrality/EvtEffGen"), 2., crec);
        } else {
          registry.fill(HIST("Events/EvtEffGen"), 2.);
        }

        auto perCollisionSample =
          tracks.sliceBy(perCol, collision.globalIndex());
        auto perCollisionASample =
          besttracks.sliceBy(perColU, collision.globalIndex());
        auto nTrkRec = countBestTracks<C, false>(
          perCollisionSample, perCollisionASample, z, crec,
          collision.trackOccupancyInTimeRange());

        if constexpr (C::template contains<aod::CentFT0Cs>()) {
          registry.fill(HIST("Events/Centrality/NtrkZvtxGen"), nTrkRec, z,
                        crec);
        } else {
          registry.fill(HIST("Events/NtrkZvtxGen"), nTrkRec, z);
        }
      }
    }

    if constexpr (C::template contains<aod::CentFT0Cs>()) {
      registry.fill(HIST("Events/Centrality/EvtEffGen"), 3., cgen);
    } else {
      registry.fill(HIST("Events/EvtEffGen"), 3.);
    }

    auto zvtxMC = mcCollision.posZ();
    auto nCharged = countPart(particles);
    if constexpr (C::template contains<aod::CentFT0Cs>()) {
      registry.fill(HIST("Events/Centrality/NtrkZvtxGen_t"), nCharged, zvtxMC,
                    cgen);
    } else {
      registry.fill(HIST("Events/NtrkZvtxGen_t"), nCharged, zvtxMC);
    }

    fillHistMC<C::template contains<aod::CentFT0Cs>()>(particles, cgen,
                                                       zvtxMC, atLeastOne);

    if (collisions.size() == 0) {
      if constexpr (C::template contains<aod::CentFT0Cs>()) {
        qaregistry.fill(HIST("Events/Centrality/NotFoundEventZvtx"),
                        mcCollision.posZ(), cgen);
      } else {
        qaregistry.fill(HIST("Events/NotFoundEventZvtx"), mcCollision.posZ());
      }
    }
  }

  /// @brief process fnc. to run on MC (inclusive, using aod::BestCollisionsFwd
  /// tracks)
  void processMCwBestTracksInclusive(
    aod::McCollisions::iterator const& mccollision,
    soa::SmallGroups<soa::Join<Colls, aod::McCollisionLabels>> const& collisions,
    FiltParticles const& particles, FiltMcMftTracks const& tracks,
    //                                      aod::BestCollisionsFwd const
    //                                      &besttracks
    FiltBestTracks const& besttracks)
  {
    processMCwBestTracks<aod::McCollisions, Colls>(
      mccollision, collisions, particles, tracks, besttracks);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processMCwBestTracksInclusive,
                 "Count MC particles using aod::BestCollisionsFwd", false);

  /// @brief process fnc. to run on MC (FT0C centrality, using
  /// aod::BestCollisionsFwd tracks)
  void processMCwBestTracksCent(
    aod::McCollisions::iterator const& mccollision,
    soa::SmallGroups<soa::Join<CollsCent, aod::McCollisionLabels>> const& collisions,
    FiltParticles const& particles, FiltMcMftTracks const& tracks,
    FiltBestTracks const& besttracks)
  {
    processMCwBestTracks<aod::McCollisions, CollsCent>(
      mccollision, collisions, particles, tracks, besttracks);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processMCwBestTracksCent,
                 "Count MC particles in FT0C bins using aod::BestCollisionsFwd",
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

      float crec = -1;
      if constexpr (C::template contains<aod::CentFT0Cs>()) {
        crec = collision.centFT0C();
      }

      auto mcCollision = collision.mcCollision();
      auto particlesPerCol = primariesI->sliceByCached(
        aod::mcparticle::mcCollisionId, mcCollision.globalIndex(), cache);
      particlesPerCol.bindExternalIndices(&tracks);

      for (auto const& particle : particlesPerCol) {
        if (!isChrgParticle(particle.pdgCode())) {
          continue;
        }
        // MC gen
        if constexpr (C::template contains<aod::CentFT0Cs>()) {
          qaregistry.fill(HIST("Tracks/Centrality/hPtPhiEtaZvtxEffGen"),
                          particle.pt(), particle.phi(), particle.eta(),
                          mcCollision.posZ(), crec);
        } else {
          qaregistry.fill(HIST("Tracks/hPtPhiEtaZvtxEffGen"), particle.pt(),
                          particle.phi(), particle.eta(), mcCollision.posZ());
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
            if constexpr (C::template contains<aod::CentFT0Cs>()) {
              if (!iscounted) { // primaries
                qaregistry.fill(HIST("Tracks/Centrality/hPtPhiEtaZvtxEffRec"),
                                particle.pt(), particle.phi(), particle.eta(),
                                mcCollision.posZ(), crec);
                iscounted = true;
              }
              if (ncnt > 1) { // duplicates
                qaregistry.fill(HIST("Tracks/Centrality/hPhiEtaDuplicates"),
                                track.phi(), track.eta(), crec);
                qaregistry.fill(
                  HIST("Tracks/Centrality/hPtPhiEtaZvtxEffDuplicates"),
                  particle.pt(), particle.phi(), particle.eta(),
                  mcCollision.posZ(), crec);
              }
            } else {
              if (!iscounted) { // primaries
                qaregistry.fill(HIST("Tracks/hPtPhiEtaZvtxEffRec"),
                                particle.pt(), particle.phi(), particle.eta(),
                                mcCollision.posZ());
                iscounted = true;
              }
              if (ncnt > 1) { // duplicates
                qaregistry.fill(HIST("Tracks/hPhiEtaDuplicates"), track.phi(),
                                track.eta());
                qaregistry.fill(HIST("Tracks/hPtPhiEtaZvtxEffDuplicates"),
                                particle.pt(), particle.phi(), particle.eta(),
                                mcCollision.posZ());
              }
            }
          }
          if constexpr (C::template contains<aod::CentFT0Cs>()) {
            qaregistry.fill(HIST("Tracks/Centrality/NmftTrkPerPart"), ncnt,
                            crec);
          } else {
            qaregistry.fill(HIST("Tracks/NmftTrkPerPart"), ncnt);
          }
          if (relatedTracks.size() > 1) {
            if constexpr (C::template contains<aod::CentFT0Cs>()) {
              qaregistry.fill(
                HIST("Tracks/Centrality/hPtPhiEtaZvtxEffGenDuplicates"),
                particle.pt(), particle.phi(), particle.eta(),
                mcCollision.posZ(), crec);
            } else {
              qaregistry.fill(HIST("Tracks/hPtPhiEtaZvtxEffGenDuplicates"),
                              particle.pt(), particle.phi(), particle.eta(),
                              mcCollision.posZ());
            }
          }
        }
      }
    }
  }

  /// @brief process function to calculate tracking efficiency (inclusive,
  /// indexed)
  void processTrkEffIdxInlusive(
    soa::Join<Colls, aod::McCollisionLabels> const& collisions,
    aod::McCollisions const& mccollisions, ParticlesI const& particles,
    MFTTracksLabeled const& tracks)
  {
    processTrkEffIdx<Colls, aod::McCollisions>(collisions, mccollisions,
                                               particles, tracks);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processTrkEffIdxInlusive,
                 "Process tracking efficiency (inclusive)", false);

  /// @brief process function to calculate tracking efficiency (FT0 bins,
  /// indexed)
  void processTrkEffIdxCent(
    soa::Join<CollsCent, aod::McCollisionLabels> const& collisions,
    aod::McCollisions const& mccollisions, ParticlesI const& particles,
    MFTTracksLabeled const& tracks)
  {
    processTrkEffIdx<CollsCent, aod::McCollisions>(collisions, mccollisions,
                                                   particles, tracks);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processTrkEffIdxCent,
                 "Process tracking efficiency in FT0 bins", false);

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

    float crec = -1;
    if constexpr (C::template contains<aod::CentFT0Cs>()) {
      crec = collision.centFT0C();
    }

    auto mcCollision = collision.mcCollision();
    auto particlesPerCol = particles.sliceByCached(
      aod::mcparticle::mcCollisionId, mcCollision.globalIndex(), cache);
    for (auto const& particle : particlesPerCol) {
      if (!isChrgParticle(particle.pdgCode())) {
        continue;
      }
      if constexpr (C::template contains<aod::CentFT0Cs>()) {
        qaregistry.fill(HIST("Tracks/Centrality/hPtPhiEtaZvtxEffBestGen"),
                        particle.pt(), particle.phi(), particle.eta(),
                        mcCollision.posZ(), crec);
      } else {
        qaregistry.fill(HIST("Tracks/hPtPhiEtaZvtxEffBestGen"), particle.pt(),
                        particle.phi(), particle.eta(), mcCollision.posZ());
      }
    }

    for (auto const& track : besttracks) {
      auto itrack = track.mfttrack_as<FiltMcMftTracks>();
      if (!isTrackSelected(itrack)) {
        continue;
      }
      if (itrack.has_mcParticle()) {
        auto particle = itrack.mcParticle_as<FiltParticles>();
        if constexpr (C::template contains<aod::CentFT0Cs>()) {
          qaregistry.fill(HIST("Tracks/Centrality/hPtPhiEtaZvtxEffBestRec"),
                          particle.pt(), itrack.phi(), itrack.eta(),
                          mcCollision.posZ(), crec);
        } else {
          qaregistry.fill(HIST("Tracks/hPtPhiEtaZvtxEffBestRec"), particle.pt(),
                          itrack.phi(), itrack.eta(), mcCollision.posZ());
        }
      } else {
        if constexpr (C::template contains<aod::CentFT0Cs>()) {
          qaregistry.fill(HIST("Tracks/Centrality/hPtEffBestFakeRec"),
                          itrack.pt(), crec);
        } else {
          qaregistry.fill(HIST("Tracks/hPtEffBestFakeRec"), itrack.pt());
        }
      }
    }
  }

  /// @brief process function to calculate tracking efficiency (inclusive, based
  /// on BestCollisionsFwd)
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

  /// @brief process function to calculate tracking efficiency (in FT0 bins,
  /// based on BestCollisionsFwd)
  void processTrkEffBestCent(
    soa::Join<CollsCent, aod::McCollisionLabels>::iterator const& collision,
    aod::McCollisions const& mccollisions, FiltParticles const& particles,
    FiltMcMftTracks const& tracks,
    soa::SmallGroups<aod::BestCollisionsFwd> const& besttracks)
  {
    processTrkEffBest<CollsCent, aod::McCollisions>(
      collision, mccollisions, particles, tracks, besttracks);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processTrkEffBestCent,
                 "Process tracking efficiency (in FT0 bins, based on BestCollisionsFwd)",
                 false);

  Preslice<FiltMftTracks> filtTrkperCol = o2::aod::fwdtrack::collisionId;

  /// @brief process template function for MC QA checks
  template <typename C>
  void processMcQA(
    typename soa::SmallGroups<soa::Join<C, aod::McCollisionLabels>> const& collisions,
    aod::McCollisions const& mcCollisions,
    FiltParticles const& /*particles*/, MFTTracksLabeled const& tracks,
    aod::AmbiguousMFTTracks const& atracks)
  {
    for (const auto& collision : collisions) {
      float crec = -1;
      if constexpr (C::template contains<aod::CentFT0Cs>()) {
        crec = collision.centFT0C();
        qaregistry.fill(
          HIST("Events/Centrality/hRecPerGenColls"),
          static_cast<float>(collisions.size()) / mcCollisions.size(), crec);
      } else {
        qaregistry.fill(HIST("Events/hRecPerGenColls"),
                        static_cast<float>(collisions.size()) /
                          mcCollisions.size());
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
      if constexpr (C::template contains<aod::CentFT0Cs>()) {
        qaregistry.fill(HIST("Tracks/Centrality/hNmftTrks"), ntracks, crec);
        qaregistry.fill(HIST("Tracks/Centrality/hFracAmbiguousMftTrks"),
                        static_cast<float>(nAtracks) / ntracks, crec);
      } else {
        qaregistry.fill(HIST("Tracks/hNmftTrks"), ntracks);
        qaregistry.fill(HIST("Tracks/hFracAmbiguousMftTrks"),
                        static_cast<float>(nAtracks) / ntracks);
      }
    }
  }

  /// @brief process function for QA checks (inclusive)
  void processMcQAInclusive(
    soa::SmallGroups<soa::Join<Colls, aod::McCollisionLabels>> const& collisions,
    aod::McCollisions const& mcCollisions, FiltParticles const& particles,
    MFTTracksLabeled const& tracks, aod::AmbiguousMFTTracks const& atracks)
  {
    processMcQA<Colls>(collisions, mcCollisions, particles, tracks, atracks);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processMcQAInclusive,
                 "Process MC QA checks (inclusive)", false);

  /// @brief process function for QA checks (in FT0 bins)
  void processMcQACent(
    soa::SmallGroups<soa::Join<CollsCent, aod::McCollisionLabels>> const& collisions,
    aod::McCollisions const& mcCollisions, FiltParticles const& particles,
    MFTTracksLabeled const& tracks, aod::AmbiguousMFTTracks const& atracks)
  {
    processMcQA<CollsCent>(collisions, mcCollisions, particles, tracks,
                           atracks);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processMcQACent,
                 "Process MC QA checks (in FT0 bins)", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<DndetaMFTPbPb>(cfgc)};
}
