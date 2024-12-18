// Copyright 2020-2022 CERN and copyright holders of ALICE O2.
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
/// \file   dndeta-mft-pbpb.cxx
/// \struct dndeta analysis at forward pseudorapidity
/// \brief Task for calculating dNdeta in Pb-Pb collisions using MFT detector
/// \author Gyula Bencedi <gyula.bencedi@cern.ch>
/// \since Nov 2024
/// @note based on dndeta-mft.cxx
///

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

AxisSpec PtAxis = {1001, -0.005, 10.005};
AxisSpec MultAxis = {701, -0.5, 700.5, "N_{trk}"};
AxisSpec ZAxis = {60, -30., 30.};
AxisSpec DeltaZAxis = {61, -6.1, 6.1};
AxisSpec DCAxyAxis = {500, -1, 50};
AxisSpec PhiAxis = {629, 0, o2::constants::math::TwoPI, "Rad", "#phi"};
AxisSpec EtaAxis = {20, -4., -2.};

struct PseudorapidityDensityMFT {
  SliceCache cache;

  // Histogram registry
  HistogramRegistry registry{
    "registry",
    {},
    OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry QAregistry{
    "QAregistry",
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
    Configurable<float> cfg_eta_min{"cfg_eta_min", -3.6f, ""};
    Configurable<float> cfg_eta_max{"cfg_eta_max", -2.5f, ""};
    Configurable<int> cfg_min_ncluster_mft{"cfg_min_ncluster_mft", 5,
                                           "minimum number of MFT clusters"};
    Configurable<double> cfg_min_pt{"cfg_min_pt", 0.,
                                    "minimum pT of the MFT tracks"};
    Configurable<bool> cfg_require_ca{
      "cfg_require_ca", false,
      "Use Cellular Automaton track-finding algorithm"};
    Configurable<float> cfg_max_dcaxy{"cfg_max_dcaxy", 2.0f, "Cut on dcaXY"};
  } trkcuts;

  // event selection conf.
  Configurable<float> cfgCutZvtx{"cfgCutZvtx", 10.0f, "Cut on z-vtx"};
  Configurable<float> cfgCutCent{"cfgCutCent", 80.0f,
                                 "Cut on maximum centrality"};
  Configurable<bool> useZDiffCut{"useZvtxDiffCut", false,
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
  ConfigurableAxis OccupancyBins{"OccupancyBins",
                                 {VARIABLE_WIDTH, 0.0f, 250.0f, 500.0f, 750.0f,
                                  1000.0f, 1500.0f, 2000.0f, 3000.0f, 4500.0f,
                                  6000.0f, 8000.0f, 10000.0f, 50000.0f},
                                 "Occupancy"};
  Configurable<float> minOccupancy{
    "minOccupancy", -1, "minimum occupancy from neighbouring collisions"};
  Configurable<float> maxOccupancy{
    "maxOccupancy", -1, "maximum occupancy from neighbouring collisions"};
  ConfigurableAxis CentBins{
    "CentBins",
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

    AxisSpec CentAxis = {CentBins, "Centrality", "CentralityAxis"};
    AxisSpec OccupancyAxis = {OccupancyBins, "Occupancy", "OccupancyAxis"};

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
                    {HistType::kTH2F, {MultAxis, ZAxis}}});
      registry.add({"Tracks/EtaZvtx",
                    "; #eta; Z_{vtx} (cm);",
                    {HistType::kTH2F, {EtaAxis, ZAxis}}});
      registry.add({"Tracks/PhiEta",
                    "; #varphi; #eta;",
                    {HistType::kTH2F, {PhiAxis, EtaAxis}}});
      QAregistry.add({"Tracks/Chi2Eta",
                      "; #chi^{2}; #it{#eta};",
                      {HistType::kTH2F, {{600, 0, 20}, {100, -8, 8}}}});
      QAregistry.add(
        {"Tracks/Chi2", "; #chi^{2};", {HistType::kTH1F, {{600, 0, 20}}}});
      QAregistry.add({"Tracks/NclustersEta",
                      "; nClusters; #eta;",
                      {HistType::kTH2F, {{7, 4, 10}, {100, -8, 8}}}});
      QAregistry.add("Events/Occupancy", "; Z_{vtx} (cm); Occupancy",
                     HistType::kTH2F, {ZAxis, OccupancyAxis}, false);

      if (doprocessDatawBestTracksInclusive) {
        registry.add({"Events/NtrkZvtxBest",
                      "; N_{trk}; Z_{vtx} (cm);",
                      {HistType::kTH2F, {MultAxis, ZAxis}}});
        registry.add({"Tracks/EtaZvtxBest",
                      "; #eta; Z_{vtx} (cm);",
                      {HistType::kTH2F, {EtaAxis, ZAxis}}});
        registry.add({"Tracks/PhiEtaBest",
                      "; #varphi; #eta;",
                      {HistType::kTH2F, {PhiAxis, EtaAxis}}});
        QAregistry.add({"Tracks/NclustersEtaBest",
                        "; nClusters; #eta;",
                        {HistType::kTH2F, {{7, 4, 10}, {100, -8, 8}}}});
        QAregistry.add({"Tracks/DCAXYPt",
                        " ; p_{T} (GeV/c) ; DCA_{XY} (cm); centrality",
                        {HistType::kTH2F, {PtAxis, DCAxyAxis}}});
        QAregistry.add({"Tracks/DCAXY",
                        " ; DCA_{XY} (cm)",
                        {HistType::kTH1F, {DCAxyAxis}}});
        QAregistry.add({"Tracks/ReTracksEtaZvtx",
                        "; #eta; #it{z}_{vtx} (cm); tracks",
                        {HistType::kTH2F, {EtaAxis, ZAxis}}});
        QAregistry.add({"Tracks/ReTracksPhiEta",
                        "; #varphi; #eta; tracks",
                        {HistType::kTH2F, {PhiAxis, EtaAxis}}});
        QAregistry.add({"Tracks/TrackAmbDegree",
                        " ; N_{coll}^{comp}",
                        {HistType::kTH1F, {{51, -0.5, 50.5}}}});
      }
    }

    if (doprocessDataCent || doprocessDatawBestTracksCent) {
      registry.add({"Events/Centrality/Selection",
                    ";status;centrality;events",
                    {HistType::kTH2F, {{2, 0.5, 2.5}, CentAxis}}});
      auto hstat = registry.get<TH2>(HIST("Events/Centrality/Selection"));
      auto* x = hstat->GetXaxis();
      x->SetBinLabel(1, "All");
      x->SetBinLabel(2, "Selected");

      registry.add({"Events/Centrality/NtrkZvtx",
                    "; N_{trk}; Z_{vtx} (cm); centrality",
                    {HistType::kTH3F, {MultAxis, ZAxis, CentAxis}}});
      registry.add({"Tracks/Centrality/EtaZvtx",
                    "; #eta; Z_{vtx} (cm); centrality",
                    {HistType::kTH3F, {EtaAxis, ZAxis, CentAxis}}});
      registry.add({"Tracks/Centrality/PhiEta",
                    "; #varphi; #eta; centrality",
                    {HistType::kTH3F, {PhiAxis, EtaAxis, CentAxis}}});
      QAregistry.add({"Events/Centrality/hcentFT0C",
                      " ; cent FT0C",
                      {HistType::kTH1F, {CentAxis}},
                      true});
      QAregistry.add(
        {"Tracks/Centrality/Chi2Eta",
         "; #chi^{2}; #it{#eta}; centrality",
         {HistType::kTH3F, {{600, 0, 20}, {100, -8, 8}, CentAxis}}});
      QAregistry.add({"Tracks/Centrality/Chi2",
                      "; #chi^{2}; centrality",
                      {HistType::kTH2F, {{600, 0, 20}, CentAxis}}});
      QAregistry.add({"Tracks/Centrality/NclustersEta",
                      "; nClusters; #eta; centrality",
                      {HistType::kTH3F, {{7, 4, 10}, {100, -8, 8}, CentAxis}}});
      QAregistry.add("Events/Centrality/Occupancy",
                     "; Z_{vtx} (cm); centrality; Occupancy", HistType::kTH3F,
                     {ZAxis, CentAxis, OccupancyAxis}, false);
      QAregistry.add("Tracks/Centrality/Occupancy", "dndeta occupancy",
                     HistType::kTHnSparseF,
                     {ZAxis, CentAxis, EtaAxis, PhiAxis, OccupancyAxis}, false);

      if (doprocessDatawBestTracksCent) {
        registry.add({"Events/Centrality/NtrkZvtxBest",
                      "; N_{trk}; Z_{vtx} (cm); centrality",
                      {HistType::kTH3F, {MultAxis, ZAxis, CentAxis}}});
        registry.add({"Tracks/Centrality/EtaZvtxBest",
                      "; #eta; Z_{vtx} (cm); centrality",
                      {HistType::kTH3F, {EtaAxis, ZAxis, CentAxis}}});
        registry.add({"Tracks/Centrality/PhiEtaBest",
                      "; #varphi; #eta; centrality",
                      {HistType::kTH3F, {PhiAxis, EtaAxis, CentAxis}}});
        QAregistry.add(
          {"Tracks/Centrality/NclustersEtaBest",
           "; nClusters; #eta; centrality",
           {HistType::kTH3F, {{7, 4, 10}, {100, -8, 8}, CentAxis}}});
        QAregistry.add({"Tracks/Centrality/TrackAmbDegree",
                        " ; N_{coll}^{comp}",
                        {HistType::kTH2F, {{51, -0.5, 50.5}, CentAxis}}});
        QAregistry.add({"Tracks/Centrality/DCAXY",
                        " ; DCA_{XY} (cm)",
                        {HistType::kTH2F, {DCAxyAxis, CentAxis}}});
        QAregistry.add({"Tracks/Centrality/DCAXYPt",
                        " ; p_{T} (GeV/c) ; DCA_{XY} (cm); centrality",
                        {HistType::kTH3F, {PtAxis, DCAxyAxis, CentAxis}}});
        QAregistry.add({"Tracks/Centrality/ReTracksEtaZvtx",
                        "; #eta; #it{z}_{vtx} (cm); tracks",
                        {HistType::kTH3F, {EtaAxis, ZAxis, CentAxis}}});
        QAregistry.add({"Tracks/Centrality/ReTracksPhiEta",
                        "; #varphi; #eta; tracks",
                        {HistType::kTH3F, {PhiAxis, EtaAxis, CentAxis}}});
        QAregistry.add("Events/Centrality/OccupancyBest",
                       "; Z_{vtx} (cm); centrality; Occupancy", HistType::kTH3F,
                       {ZAxis, CentAxis, OccupancyAxis}, false);
        QAregistry.add("Tracks/Centrality/OccupancyBest", "dndeta occupancy",
                       HistType::kTHnSparseF,
                       {ZAxis, CentAxis, EtaAxis, PhiAxis, OccupancyAxis},
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
                    {HistType::kTH2F, {MultAxis, ZAxis}}});
      registry.add({"Events/NtrkZvtxGen",
                    "; N_{trk}; Z_{vtx} (cm);",
                    {HistType::kTH2F, {MultAxis, ZAxis}}});
      registry.add({"Tracks/EtaZvtxGen",
                    "; #eta; Z_{vtx} (cm);",
                    {HistType::kTH2F, {EtaAxis, ZAxis}}});
      registry.add({"Tracks/PhiEtaGen",
                    "; #varphi; #eta;",
                    {HistType::kTH2F, {PhiAxis, EtaAxis}}});
      registry.add({"Tracks/EtaZvtxGen_t",
                    "; #eta; Z_{vtx} (cm);",
                    {HistType::kTH2F, {EtaAxis, ZAxis}}});
      registry.add({"Tracks/PhiEtaGen_t",
                    "; #varphi; #eta;",
                    {HistType::kTH2F, {PhiAxis, EtaAxis}}});
      QAregistry.add({"Events/NotFoundEventZvtx",
                      " ; #it{z}_{vtx} (cm)",
                      {HistType::kTH1F, {ZAxis}}});
      QAregistry.add({"Events/ZvtxDiff",
                      " ; Z_{rec} - Z_{gen} (cm)",
                      {HistType::kTH1F, {DeltaZAxis}}});
      QAregistry.add(
        {"Events/SplitMult", " ; N_{gen}", {HistType::kTH1F, {MultAxis}}});
    }

    if (doprocessMCCent || doprocessMCwBestTracksCent) {
      registry.add({"Events/Centrality/EvtEffGen",
                    ";status;events",
                    {HistType::kTH2F, {{3, 0.5, 3.5}, CentAxis}}});
      auto heff = registry.get<TH2>(HIST("Events/Centrality/EvtEffGen"));
      auto* x = heff->GetXaxis();
      x->SetBinLabel(1, "All reconstructed");
      x->SetBinLabel(2, "Selected reconstructed");
      x->SetBinLabel(3, "All generated");

      registry.add({"Events/Centrality/NtrkZvtxGen_t",
                    "; N_{trk}; Z_{vtx} (cm);",
                    {HistType::kTH3F, {MultAxis, ZAxis, CentAxis}}});
      registry.add({"Events/Centrality/NtrkZvtxGen",
                    "; N_{trk}; Z_{vtx} (cm);",
                    {HistType::kTH3F, {MultAxis, ZAxis, CentAxis}}});
      registry.add({"Events/Centrality/hRecCent",
                    "Events/Centrality/hRecCent",
                    {HistType::kTH1F, {CentAxis}}});
      registry.add({"Events/Centrality/hRecZvtxCent",
                    "Events/Centrality/hRecZvtxCent",
                    {HistType::kTH2F, {ZAxis, CentAxis}}});
      registry.add({"Tracks/Centrality/EtaZvtxGen",
                    "; #eta; Z_{vtx} (cm);",
                    {HistType::kTH3F, {EtaAxis, ZAxis, CentAxis}}});
      registry.add({"Tracks/Centrality/PhiEtaGen",
                    "; #varphi; #eta;",
                    {HistType::kTH3F, {PhiAxis, EtaAxis, CentAxis}}});
      registry.add({"Tracks/Centrality/EtaZvtxGen_t",
                    "; #eta; Z_{vtx} (cm);",
                    {HistType::kTH3F, {EtaAxis, ZAxis, CentAxis}}});
      registry.add({"Tracks/Centrality/PhiEtaGen_t",
                    "; #varphi; #eta;",
                    {HistType::kTH3F, {PhiAxis, EtaAxis, CentAxis}}});
      QAregistry.add({"Events/Centrality/NotFoundEventZvtx",
                      " ; #it{z}_{vtx} (cm)",
                      {HistType::kTH2F, {ZAxis, CentAxis}}});
      QAregistry.add({"Events/Centrality/ZvtxDiff",
                      " ; Z_{rec} - Z_{gen} (cm)",
                      {HistType::kTH2F, {DeltaZAxis, CentAxis}}});
      QAregistry.add({"Events/Centrality/SplitMult",
                      " ; N_{gen}",
                      {HistType::kTH2F, {MultAxis, CentAxis}}});
    }

    if (doprocessTrkEffIdxInlusive) {
      QAregistry.add(
        {"Tracks/hPtPhiEtaZvtxEffGen",
         "hPtPhiEtaZvtxEffGen",
         {HistType::kTHnSparseF, {PtAxis, PhiAxis, EtaAxis, ZAxis}}});
      QAregistry.add(
        {"Tracks/hPtPhiEtaZvtxEffRec",
         "hPtPhiEtaZvtxEffRec",
         {HistType::kTHnSparseF, {PtAxis, PhiAxis, EtaAxis, ZAxis}}});
      QAregistry.add({"Tracks/hPhiEtaDuplicates",
                      " ; p_{T} (GeV/c);",
                      {HistType::kTH2F, {PhiAxis, EtaAxis}}});
      QAregistry.add(
        {"Tracks/hPtPhiEtaZvtxEffDuplicates",
         "hPtPhiEtaZvtxEffDuplicates",
         {HistType::kTHnSparseF, {PtAxis, PhiAxis, EtaAxis, ZAxis}}});
      QAregistry.add(
        {"Tracks/hPtPhiEtaZvtxEffGenDuplicates",
         "hPtPhiEtaZvtxEffGenDuplicates",
         {HistType::kTHnSparseF, {PtAxis, PhiAxis, EtaAxis, ZAxis}}});
      QAregistry.add({"Tracks/NmftTrkPerPart",
                      "; #it{N}_{mft tracks per particle};",
                      {HistType::kTH1F, {{200, -0.5, 200.}}}});
    }

    if (doprocessTrkEffIdxCent) {
      QAregistry.add({"Tracks/Centrality/hPtPhiEtaZvtxEffGen",
                      "hPtPhiEtaZvtxEffGen",
                      {HistType::kTHnSparseF,
                       {PtAxis, PhiAxis, EtaAxis, ZAxis, CentAxis}}});
      QAregistry.add({"Tracks/Centrality/hPtPhiEtaZvtxEffRec",
                      "hPtPhiEtaZvtxEffRec",
                      {HistType::kTHnSparseF,
                       {PtAxis, PhiAxis, EtaAxis, ZAxis, CentAxis}}});
      QAregistry.add({"Tracks/Centrality/hPhiEtaDuplicates",
                      " ; p_{T} (GeV/c);",
                      {HistType::kTH3F, {PhiAxis, EtaAxis, CentAxis}}});
      QAregistry.add({"Tracks/Centrality/hPtPhiEtaZvtxEffDuplicates",
                      "hPtPhiEtaZvtxEffDuplicates",
                      {HistType::kTHnSparseF,
                       {PtAxis, PhiAxis, EtaAxis, ZAxis, CentAxis}}});
      QAregistry.add({"Tracks/Centrality/hPtPhiEtaZvtxEffGenDuplicates",
                      "hPtPhiEtaZvtxEffGenDuplicates",
                      {HistType::kTHnSparseF,
                       {PtAxis, PhiAxis, EtaAxis, ZAxis, CentAxis}}});
      QAregistry.add({"Tracks/Centrality/NmftTrkPerPart",
                      "; #it{N}_{mft tracks per particle};",
                      {HistType::kTH2F, {{200, -0.5, 200.}, CentAxis}}});
    }

    if (doprocessTrkEffBestInclusive) {
      QAregistry.add(
        {"Tracks/hPtPhiEtaZvtxEffBestGen",
         "hPtPhiEtaZvtxEffGen",
         {HistType::kTHnSparseF, {PtAxis, PhiAxis, EtaAxis, ZAxis}}});
      QAregistry.add(
        {"Tracks/hPtPhiEtaZvtxEffBestRec",
         "hPtPhiEtaZvtxEffRec",
         {HistType::kTHnSparseF, {PtAxis, PhiAxis, EtaAxis, ZAxis}}});
      QAregistry.add({"Tracks/hPtEffBestFakeRec",
                      " ; p_{T} (GeV/c);",
                      {HistType::kTH1F, {PtAxis}}});
    }

    if (doprocessTrkEffBestCent) {
      QAregistry.add({"Tracks/Centrality/hPtPhiEtaZvtxEffBestGen",
                      "hPtPhiEtaZvtxEffGen",
                      {HistType::kTHnSparseF,
                       {PtAxis, PhiAxis, EtaAxis, ZAxis, CentAxis}}});
      QAregistry.add({"Tracks/Centrality/hPtPhiEtaZvtxEffBestRec",
                      "hPtPhiEtaZvtxEffRec",
                      {HistType::kTHnSparseF,
                       {PtAxis, PhiAxis, EtaAxis, ZAxis, CentAxis}}});
      QAregistry.add({"Tracks/Centrality/hPtEffBestFakeRec",
                      " ; p_{T} (GeV/c);",
                      {HistType::kTH2F, {PtAxis, CentAxis}}});
    }

    if (doprocessMcQAInclusive) {
      QAregistry.add({"Events/hRecPerGenColls",
                      "; #it{N}_{reco collisions} / #it{N}_{gen collisions};",
                      {HistType::kTH1F, {{200, 0., 2.}}}});
      QAregistry.add({"Tracks/hNmftTrks",
                      "; #it{N}_{mft tracks};",
                      {HistType::kTH1F, {{200, -0.5, 200.}}}});
      QAregistry.add({"Tracks/hFracAmbiguousMftTrks",
                      "; #it{N}_{ambiguous tracks} / #it{N}_{tracks};",
                      {HistType::kTH1F, {{100, 0., 1.}}}});
    }

    if (doprocessMcQACent) {
      QAregistry.add(
        {"Events/Centrality/hRecPerGenColls",
         "; #it{N}_{reco collisions} / #it{N}_{gen collisions}; centrality",
         {HistType::kTH2F, {{200, 0., 2.}, CentAxis}}});
      QAregistry.add({"Tracks/Centrality/hNmftTrks",
                      "; #it{N}_{mft tracks}; centrality",
                      {HistType::kTH2F, {{200, -0.5, 200.}, CentAxis}}});
      QAregistry.add(
        {"Tracks/Centrality/hFracAmbiguousMftTrks",
         "; #it{N}_{ambiguous tracks} / #it{N}_{tracks}; centrality",
         {HistType::kTH2F, {{100, 0., 1.}, CentAxis}}});
    }
  }

  /// Filters - tracks
  Filter filtTrkEta = (aod::fwdtrack::eta < trkcuts.cfg_eta_max) &&
                      (aod::fwdtrack::eta > trkcuts.cfg_eta_min);
  Filter filtATrackID = (aod::fwdtrack::bestCollisionId >= 0);
  Filter filtATrackDCA =
    (nabs(aod::fwdtrack::bestDCAXY) < trkcuts.cfg_max_dcaxy);

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
  using filtMftTracks = soa::Filtered<aod::MFTTracks>;
  using filtMcMftTracks = soa::Filtered<MFTTracksLabeled>;
  using filtBestTracks = soa::Filtered<aod::BestCollisionsFwd>;
  using filtParticles = soa::Filtered<aod::McParticles>;

  template <typename T>
  bool isTrackSelected(const T& track)
  {
    if (track.eta() < trkcuts.cfg_eta_min || track.eta() > trkcuts.cfg_eta_max)
      return false;
    if (trkcuts.cfg_require_ca && !track.isCA())
      return false;
    if (track.nClusters() < trkcuts.cfg_min_ncluster_mft)
      return false;
    if (track.pt() < trkcuts.cfg_min_pt)
      return false;
    if (usePhiCut) {
      float phi = track.phi();
      o2::math_utils::bringTo02Pi(phi);
      if (phi < 0.f || 2.f * M_PI < phi) {
        return false;
      }
      if ((phi < cfgPhiCut) ||
          ((phi > M_PI - cfgPhiCut) && (phi < M_PI + cfgPhiCut)) ||
          (phi > 2. * M_PI - cfgPhiCut) ||
          ((phi > ((M_PI / 2. - 0.1) * M_PI) - cfgPhiCut) &&
           (phi < ((M_PI / 2. - 0.1) * M_PI) + cfgPhiCut)))
        return false;
    }
    return true;
  }

  template <typename C, bool fillHis = false, typename T>
  int countTracks(T const& tracks, float z, float c, float occ)
  {
    auto nTrk = 0;
    if (tracks.size() > 0) {
      for (auto& track : tracks) {
        if (fillHis) {
          if constexpr (C::template contains<aod::CentFT0Cs>()) {
            QAregistry.fill(HIST("Tracks/Centrality/Chi2Eta"), track.chi2(),
                            track.eta(), c);
            QAregistry.fill(HIST("Tracks/Centrality/Chi2"), track.chi2(), c);
            QAregistry.fill(HIST("Tracks/Centrality/NclustersEta"),
                            track.nClusters(), track.eta(), c);
          } else {
            QAregistry.fill(HIST("Tracks/Chi2Eta"), track.chi2(), track.eta());
            QAregistry.fill(HIST("Tracks/Chi2"), track.chi2());
            QAregistry.fill(HIST("Tracks/NclustersEta"), track.nClusters(),
                            track.eta());
          }
        }
        if (!isTrackSelected(track)) {
          continue;
        }
        if (fillHis) {
          float phi = track.phi();
          o2::math_utils::bringTo02Pi(phi);
          if (phi < 0.f || 2.f * M_PI < phi) {
            continue;
          }
          if constexpr (C::template contains<aod::CentFT0Cs>()) {
            registry.fill(HIST("Tracks/Centrality/EtaZvtx"), track.eta(), z, c);
            registry.fill(HIST("Tracks/Centrality/PhiEta"), phi, track.eta(),
                          c);
            QAregistry.fill(HIST("Tracks/Centrality/Occupancy"), z, c,
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
      for (auto& atrack : besttracks) {
        auto itrack = atrack.template mfttrack_as<T>();
        if (!isTrackSelected(itrack)) {
          continue;
        }
        if (fillHis) {
          float phi = itrack.phi();
          o2::math_utils::bringTo02Pi(phi);
          if (phi < 0.f || 2.f * M_PI < phi) {
            continue;
          }
          if constexpr (C::template contains<aod::CentFT0Cs>()) {
            registry.fill(HIST("Tracks/Centrality/EtaZvtxBest"), itrack.eta(),
                          z, c);
            registry.fill(HIST("Tracks/Centrality/PhiEtaBest"), phi,
                          itrack.eta(), c);
            QAregistry.fill(HIST("Tracks/Centrality/OccupancyBest"), z, c,
                            itrack.eta(), itrack.phi(), occ);
            QAregistry.fill(HIST("Tracks/Centrality/DCAXYPt"), itrack.pt(),
                            atrack.bestDCAXY(), c);
            QAregistry.fill(HIST("Tracks/Centrality/DCAXY"), atrack.bestDCAXY(),
                            c);
            QAregistry.fill(HIST("Tracks/Centrality/NclustersEtaBest"),
                            itrack.nClusters(), itrack.eta(), c);
            if (itrack.collisionId() != atrack.bestCollisionId()) {
              QAregistry.fill(HIST("Tracks/Centrality/ReTracksEtaZvtx"),
                              itrack.eta(), z, c);
              QAregistry.fill(HIST("Tracks/Centrality/ReTracksPhiEta"), phi,
                              itrack.eta(), c);
            }
            QAregistry.fill(HIST("Tracks/Centrality/TrackAmbDegree"),
                            atrack.ambDegree(), c);
          } else {
            registry.fill(HIST("Tracks/EtaZvtxBest"), itrack.eta(), z);
            registry.fill(HIST("Tracks/PhiEtaBest"), phi, itrack.eta());
            QAregistry.fill(HIST("Tracks/DCAXYPt"), itrack.pt(),
                            atrack.bestDCAXY());
            QAregistry.fill(HIST("Tracks/DCAXY"), atrack.bestDCAXY());
            QAregistry.fill(HIST("Tracks/NclustersEtaBest"), itrack.nClusters(),
                            itrack.eta());
            if (itrack.collisionId() != atrack.bestCollisionId()) {
              QAregistry.fill(HIST("Tracks/ReTracksEtaZvtx"), itrack.eta(), z);
              QAregistry.fill(HIST("Tracks/ReTracksPhiEta"), phi, itrack.eta());
            }
            QAregistry.fill(HIST("Tracks/TrackAmbDegree"), atrack.ambDegree());
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
    for (auto& particle : particles) {
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
  void fillHist_MC(P const& particles, float cent, float zvtx,
                   bool const atLeastOne)
  {
    for (auto& particle : particles) {
      if (!isChrgParticle(particle.pdgCode())) {
        continue;
      }

      float phi = particle.phi();
      o2::math_utils::bringTo02Pi(phi);
      if (phi < 0.f || 2.f * M_PI < phi) {
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
        if (phi < 0.f || 2.f * M_PI < phi) {
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
    for (auto& bc : bcs) {
      if ((bc.selection_bit(aod::evsel::kIsBBT0A) &&
           bc.selection_bit(aod::evsel::kIsBBT0C)) != 0) {
        registry.fill(HIST("hBcSel"), 0);
        cols.clear();
        for (auto& collision : collisions) {
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

  PROCESS_SWITCH(PseudorapidityDensityMFT, processTagging,
                 "Collect event sample stats", true);

  template <typename C>
  void processData(typename C::iterator const& collision,
                   filtMftTracks const& tracks)
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
      QAregistry.fill(HIST("Events/Centrality/Occupancy"), z, c, occ);
      QAregistry.fill(HIST("Events/Centrality/hcentFT0C"), c);
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
    typename C::iterator const& collision, filtMftTracks const& tracks,
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
      QAregistry.fill(HIST("Events/Centrality/OccupancyBest"), z, c, occ);
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
                            filtMftTracks const& tracks)
  {
    processData<Colls>(collision, tracks);
  }

  PROCESS_SWITCH(PseudorapidityDensityMFT, processDataInclusive, "Count tracks",
                 false);

  /// @brief process fnc. to run on DATA and REC MC w/ FT0C centrality selection
  void processDataCent(CollsCent::iterator const& collision,
                       filtMftTracks const& tracks)
  {
    processData<CollsCent>(collision, tracks);
  }

  PROCESS_SWITCH(PseudorapidityDensityMFT, processDataCent,
                 "Count tracks in FT0C bins", false);

  /// @brief process fnc. to run on DATA and REC MC based on BestCollisionsFwd
  /// table w/o centrality selection
  void processDatawBestTracksInclusive(
    Colls::iterator const& collision, filtMftTracks const& tracks,
    soa::SmallGroups<aod::BestCollisionsFwd> const& besttracks)
  {
    processDatawBestTracks<Colls>(collision, tracks, besttracks);
  }

  PROCESS_SWITCH(PseudorapidityDensityMFT, processDatawBestTracksInclusive,
                 "Count tracks based on BestCollisionsFwd table", false);

  /// @brief process fnc. to run on DATA and REC MC based on BestCollisionsFwd
  /// table w/ FT0C centrality selection
  void processDatawBestTracksCent(
    CollsCent::iterator const& collision, filtMftTracks const& tracks,
    soa::SmallGroups<aod::BestCollisionsFwd> const& besttracks)
  {
    processDatawBestTracks<CollsCent>(collision, tracks, besttracks);
  }

  PROCESS_SWITCH(PseudorapidityDensityMFT, processDatawBestTracksCent,
                 "Count tracks in FT0C bins based on BestCollisionsFwd table",
                 false);

  Preslice<filtMcMftTracks> perCol = o2::aod::fwdtrack::collisionId;
  Partition<filtParticles> mcSample = nabs(aod::mcparticle::eta) < 1.0f;

  /// @brief process template function to run on MC truth
  /// @param cols subscribe to the collisions
  /// @param parts subscribe to filtered MC particle table
  template <typename MC, typename C>
  void processMC(
    typename MC::iterator const& mcCollision,
    soa::SmallGroups<soa::Join<C, aod::McCollisionLabels>> const& collisions,
    filtParticles const& particles, filtMcMftTracks const& tracks)
  {
    float c_gen = -1;
    bool atLeastOne = false;
    int moreThanOne = 0;
    for (auto& collision : collisions) {
      float c_rec = -1;
      if constexpr (C::template contains<aod::CentFT0Cs>()) {
        c_rec = collision.centFT0C();
        registry.fill(HIST("Events/Centrality/EvtEffGen"), 1., c_rec);
      } else {
        registry.fill(HIST("Events/EvtEffGen"), 1.);
      }

      if (isGoodEvent<false>(collision)) {
        if constexpr (C::template contains<aod::CentFT0Cs>()) {
          if (!atLeastOne) {
            c_gen = c_rec;
          }
        }
        atLeastOne = true;
        ++moreThanOne;
        auto z = collision.posZ();

        if constexpr (C::template contains<aod::CentFT0Cs>()) {
          registry.fill(HIST("Events/Centrality/EvtEffGen"), 2., c_rec);
          registry.fill(HIST("Events/Centrality/hRecCent"), c_rec);
          registry.fill(HIST("Events/Centrality/hRecZvtxCent"), z, c_rec);
        } else {
          registry.fill(HIST("Events/EvtEffGen"), 2.);
        }

        auto perCollisionSample =
          tracks.sliceBy(perCol, collision.globalIndex());
        auto nTrkRec =
          countTracks<C, true>(perCollisionSample, z, c_rec,
                               collision.trackOccupancyInTimeRange());

        if constexpr (C::template contains<aod::CentFT0Cs>()) {
          QAregistry.fill(HIST("Events/Centrality/ZvtxDiff"),
                          collision.posZ() - mcCollision.posZ(), c_rec);
        } else {
          QAregistry.fill(HIST("Events/ZvtxDiff"),
                          collision.posZ() - mcCollision.posZ());
        }

        if (useZDiffCut) {
          if (std::abs(collision.posZ() - mcCollision.posZ()) > maxZvtxDiff) {
            continue;
          }
        }

        if constexpr (C::template contains<aod::CentFT0Cs>()) {
          registry.fill(HIST("Events/Centrality/NtrkZvtxGen"), nTrkRec,
                        collision.posZ(), c_rec);
        } else {
          registry.fill(HIST("Events/NtrkZvtxGen"), nTrkRec, collision.posZ());
        }
      }
    }

    if constexpr (C::template contains<aod::CentFT0Cs>()) {
      registry.fill(HIST("Events/Centrality/EvtEffGen"), 3., c_gen);
    } else {
      registry.fill(HIST("Events/EvtEffGen"), 3.);
    }

    auto perCollMCsample = mcSample->sliceByCached(
      aod::mcparticle::mcCollisionId, mcCollision.globalIndex(), cache);
    auto Nchrg = countPart(perCollMCsample);
    if (moreThanOne > 1) {
      if constexpr (C::template contains<aod::CentFT0Cs>()) {
        QAregistry.fill(HIST("Events/Centrality/SplitMult"), Nchrg, c_gen);
      } else {
        QAregistry.fill(HIST("Events/SplitMult"), Nchrg);
      }
    }

    auto zvtxMC = mcCollision.posZ();
    auto nCharged = countPart(particles);
    if constexpr (C::template contains<aod::CentFT0Cs>()) {
      registry.fill(HIST("Events/Centrality/NtrkZvtxGen_t"), nCharged, zvtxMC,
                    c_gen);
    } else {
      registry.fill(HIST("Events/NtrkZvtxGen_t"), nCharged, zvtxMC);
    }

    fillHist_MC<C::template contains<aod::CentFT0Cs>()>(particles, c_gen,
                                                        zvtxMC, atLeastOne);

    if (collisions.size() == 0) {
      if constexpr (C::template contains<aod::CentFT0Cs>()) {
        QAregistry.fill(HIST("Events/Centrality/NotFoundEventZvtx"),
                        mcCollision.posZ(), c_gen);
      } else {
        QAregistry.fill(HIST("Events/NotFoundEventZvtx"), mcCollision.posZ());
      }
    }
  }

  /// @brief process fnc. to run on MC w/o centrality selection
  void processMCInclusive(
    aod::McCollisions::iterator const& mccollision,
    soa::SmallGroups<soa::Join<Colls, aod::McCollisionLabels>> const& collisions,
    filtParticles const& particles, filtMcMftTracks const& tracks)
  {
    processMC<aod::McCollisions, Colls>(mccollision, collisions, particles,
                                        tracks);
  }

  PROCESS_SWITCH(PseudorapidityDensityMFT, processMCInclusive,
                 "Count MC particles", false);

  /// @brief process fnc. to run on MC w FT0C centrality selection
  void processMCCent(
    aod::McCollisions::iterator const& mccollision,
    soa::SmallGroups<soa::Join<CollsCent, aod::McCollisionLabels>> const& collisions,
    filtParticles const& particles, filtMcMftTracks const& tracks)
  {
    processMC<aod::McCollisions, CollsCent>(mccollision, collisions, particles,
                                            tracks);
  }

  PROCESS_SWITCH(PseudorapidityDensityMFT, processMCCent,
                 "Count MC particles in FT0C bins", false);

  PresliceUnsorted<aod::BestCollisionsFwd> perColU =
    aod::fwdtrack::bestCollisionId;

  /// @brief process template function to run on MC truth using
  /// aod::BestCollisionsFwd tracks
  template <typename MC, typename C>
  void processMCwBestTracks(
    typename MC::iterator const& mcCollision,
    soa::SmallGroups<soa::Join<C, aod::McCollisionLabels>> const& collisions,
    filtParticles const& particles, filtMcMftTracks const& tracks,
    filtBestTracks const& besttracks)
  {
    float c_gen = -1;
    bool atLeastOne = false;
    // int moreThanOne = 0;
    for (auto& collision : collisions) {
      float c_rec = -1;
      if constexpr (C::template contains<aod::CentFT0Cs>()) {
        c_rec = collision.centFT0C();
        registry.fill(HIST("Events/Centrality/EvtEffGen"), 1., c_rec);
      } else {
        registry.fill(HIST("Events/EvtEffGen"), 1.);
      }

      if (isGoodEvent<false>(collision)) {
        if constexpr (C::template contains<aod::CentFT0Cs>()) {
          if (!atLeastOne) {
            c_gen = c_rec;
          }
        }
        atLeastOne = true;
        // ++moreThanOne;
        auto z = collision.posZ();

        if constexpr (C::template contains<aod::CentFT0Cs>()) {
          registry.fill(HIST("Events/Centrality/EvtEffGen"), 2., c_rec);
        } else {
          registry.fill(HIST("Events/EvtEffGen"), 2.);
        }

        auto perCollisionSample =
          tracks.sliceBy(perCol, collision.globalIndex());
        auto perCollisionASample =
          besttracks.sliceBy(perColU, collision.globalIndex());
        auto nTrkRec = countBestTracks<C, false>(
          perCollisionSample, perCollisionASample, z, c_rec,
          collision.trackOccupancyInTimeRange());

        if constexpr (C::template contains<aod::CentFT0Cs>()) {
          registry.fill(HIST("Events/Centrality/NtrkZvtxGen"), nTrkRec, z,
                        c_rec);
        } else {
          registry.fill(HIST("Events/NtrkZvtxGen"), nTrkRec, z);
        }
      }
    }

    if constexpr (C::template contains<aod::CentFT0Cs>()) {
      registry.fill(HIST("Events/Centrality/EvtEffGen"), 3., c_gen);
    } else {
      registry.fill(HIST("Events/EvtEffGen"), 3.);
    }

    auto zvtxMC = mcCollision.posZ();
    auto nCharged = countPart(particles);
    if constexpr (C::template contains<aod::CentFT0Cs>()) {
      registry.fill(HIST("Events/Centrality/NtrkZvtxGen_t"), nCharged, zvtxMC,
                    c_gen);
    } else {
      registry.fill(HIST("Events/NtrkZvtxGen_t"), nCharged, zvtxMC);
    }

    fillHist_MC<C::template contains<aod::CentFT0Cs>()>(particles, c_gen,
                                                        zvtxMC, atLeastOne);

    if (collisions.size() == 0) {
      if constexpr (C::template contains<aod::CentFT0Cs>()) {
        QAregistry.fill(HIST("Events/Centrality/NotFoundEventZvtx"),
                        mcCollision.posZ(), c_gen);
      } else {
        QAregistry.fill(HIST("Events/NotFoundEventZvtx"), mcCollision.posZ());
      }
    }
  }

  /// @brief process fnc. to run on MC (inclusive, using aod::BestCollisionsFwd
  /// tracks)
  void processMCwBestTracksInclusive(
    aod::McCollisions::iterator const& mccollision,
    soa::SmallGroups<soa::Join<Colls, aod::McCollisionLabels>> const& collisions,
    filtParticles const& particles, filtMcMftTracks const& tracks,
    //                                      aod::BestCollisionsFwd const
    //                                      &besttracks
    filtBestTracks const& besttracks)
  {
    processMCwBestTracks<aod::McCollisions, Colls>(
      mccollision, collisions, particles, tracks, besttracks);
  }

  PROCESS_SWITCH(PseudorapidityDensityMFT, processMCwBestTracksInclusive,
                 "Count MC particles using aod::BestCollisionsFwd", false);

  /// @brief process fnc. to run on MC (FT0C centrality, using
  /// aod::BestCollisionsFwd tracks)
  void processMCwBestTracksCent(
    aod::McCollisions::iterator const& mccollision,
    soa::SmallGroups<soa::Join<CollsCent, aod::McCollisionLabels>> const& collisions,
    filtParticles const& particles, filtMcMftTracks const& tracks,
    filtBestTracks const& besttracks)
  {
    processMCwBestTracks<aod::McCollisions, CollsCent>(
      mccollision, collisions, particles, tracks, besttracks);
  }

  PROCESS_SWITCH(PseudorapidityDensityMFT, processMCwBestTracksCent,
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
    typename soa::Filtered<soa::Join<C, aod::McCollisionLabels>> const& collisions,
    MC const& /*mccollisions*/, ParticlesI const& /*particles*/,
    MFTTracksLabeled const& tracks)
  {
    for (auto& collision : collisions) {
      if (!isGoodEvent<false>(collision)) {
        continue;
      }
      if (!collision.has_mcCollision()) {
        continue;
      }

      float c_rec = -1;
      if constexpr (C::template contains<aod::CentFT0Cs>()) {
        c_rec = collision.centFT0C();
      }

      auto mcCollision = collision.mcCollision();
      auto particlesPerCol = primariesI->sliceByCached(
        aod::mcparticle::mcCollisionId, mcCollision.globalIndex(), cache);
      particlesPerCol.bindExternalIndices(&tracks);

      for (auto& particle : particlesPerCol) {
        if (!isChrgParticle(particle.pdgCode())) {
          continue;
        }
        // MC gen
        if constexpr (C::template contains<aod::CentFT0Cs>()) {
          QAregistry.fill(HIST("Tracks/Centrality/hPtPhiEtaZvtxEffGen"),
                          particle.pt(), particle.phi(), particle.eta(),
                          mcCollision.posZ(), c_rec);
        } else {
          QAregistry.fill(HIST("Tracks/hPtPhiEtaZvtxEffGen"), particle.pt(),
                          particle.phi(), particle.eta(), mcCollision.posZ());
        }
        // MC rec
        if (particle.has_mfttracks()) {
          auto iscounted = false;
          auto ncnt = 0;
          auto relatedTracks =
            particle.template mfttracks_as<MFTTracksLabeled>();
          for (auto& track : relatedTracks) {
            if (!isTrackSelected(track)) {
              continue;
            }
            ++ncnt;
            if constexpr (C::template contains<aod::CentFT0Cs>()) {
              if (!iscounted) { // primaries
                QAregistry.fill(HIST("Tracks/Centrality/hPtPhiEtaZvtxEffRec"),
                                particle.pt(), particle.phi(), particle.eta(),
                                mcCollision.posZ(), c_rec);
                iscounted = true;
              }
              if (ncnt > 1) { // duplicates
                QAregistry.fill(HIST("Tracks/Centrality/hPhiEtaDuplicates"),
                                track.phi(), track.eta(), c_rec);
                QAregistry.fill(
                  HIST("Tracks/Centrality/hPtPhiEtaZvtxEffDuplicates"),
                  particle.pt(), particle.phi(), particle.eta(),
                  mcCollision.posZ(), c_rec);
              }
            } else {
              if (!iscounted) { // primaries
                QAregistry.fill(HIST("Tracks/hPtPhiEtaZvtxEffRec"),
                                particle.pt(), particle.phi(), particle.eta(),
                                mcCollision.posZ());
                iscounted = true;
              }
              if (ncnt > 1) { // duplicates
                QAregistry.fill(HIST("Tracks/hPhiEtaDuplicates"), track.phi(),
                                track.eta());
                QAregistry.fill(HIST("Tracks/hPtPhiEtaZvtxEffDuplicates"),
                                particle.pt(), particle.phi(), particle.eta(),
                                mcCollision.posZ());
              }
            }
          }
          if constexpr (C::template contains<aod::CentFT0Cs>()) {
            QAregistry.fill(HIST("Tracks/Centrality/NmftTrkPerPart"), ncnt,
                            c_rec);
          } else {
            QAregistry.fill(HIST("Tracks/NmftTrkPerPart"), ncnt);
          }
          if (relatedTracks.size() > 1) {
            if constexpr (C::template contains<aod::CentFT0Cs>()) {
              QAregistry.fill(
                HIST("Tracks/Centrality/hPtPhiEtaZvtxEffGenDuplicates"),
                particle.pt(), particle.phi(), particle.eta(),
                mcCollision.posZ(), c_rec);
            } else {
              QAregistry.fill(HIST("Tracks/hPtPhiEtaZvtxEffGenDuplicates"),
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
    soa::Filtered<soa::Join<Colls, aod::McCollisionLabels>> const& collisions,
    aod::McCollisions const& mccollisions, ParticlesI const& particles,
    MFTTracksLabeled const& tracks)
  {
    processTrkEffIdx<Colls, aod::McCollisions>(collisions, mccollisions,
                                               particles, tracks);
  }

  PROCESS_SWITCH(PseudorapidityDensityMFT, processTrkEffIdxInlusive,
                 "Process tracking efficiency (inclusive)", false);

  /// @brief process function to calculate tracking efficiency (FT0 bins,
  /// indexed)
  void processTrkEffIdxCent(
    soa::Filtered<soa::Join<CollsCent, aod::McCollisionLabels>> const& collisions,
    aod::McCollisions const& mccollisions, ParticlesI const& particles,
    MFTTracksLabeled const& tracks)
  {
    processTrkEffIdx<CollsCent, aod::McCollisions>(collisions, mccollisions,
                                                   particles, tracks);
  }

  PROCESS_SWITCH(PseudorapidityDensityMFT, processTrkEffIdxCent,
                 "Process tracking efficiency in FT0 bins", false);

  /// @brief process function to calculate tracking efficiency (indexed) based
  /// on BestCollisionsFwd in FT0C bins
  template <typename C, typename MC>
  void processTrkEffBest(
    typename soa::Filtered<
      soa::Join<C, aod::McCollisionLabels>>::iterator const& collision,
    MC const& /*mccollisions*/, filtParticles const& particles,
    filtMcMftTracks const& /*tracks*/,
    soa::SmallGroups<aod::BestCollisionsFwd> const& besttracks)
  {
    if (!isGoodEvent<false>(collision)) {
      return;
    }
    if (!collision.has_mcCollision()) {
      return;
    }

    float c_rec = -1;
    if constexpr (C::template contains<aod::CentFT0Cs>()) {
      c_rec = collision.centFT0C();
    }

    auto mcCollision = collision.mcCollision();
    auto particlesPerCol = particles.sliceByCached(
      aod::mcparticle::mcCollisionId, mcCollision.globalIndex(), cache);
    for (auto& particle : particlesPerCol) {
      if (!isChrgParticle(particle.pdgCode())) {
        continue;
      }
      if constexpr (C::template contains<aod::CentFT0Cs>()) {
        QAregistry.fill(HIST("Tracks/Centrality/hPtPhiEtaZvtxEffBestGen"),
                        particle.pt(), particle.phi(), particle.eta(),
                        mcCollision.posZ(), c_rec);
      } else {
        QAregistry.fill(HIST("Tracks/hPtPhiEtaZvtxEffBestGen"), particle.pt(),
                        particle.phi(), particle.eta(), mcCollision.posZ());
      }
    }

    for (auto const& track : besttracks) {
      auto itrack = track.mfttrack_as<filtMcMftTracks>();
      if (!isTrackSelected(itrack)) {
        continue;
      }
      if (itrack.has_mcParticle()) {
        auto particle = itrack.mcParticle_as<filtParticles>();
        if constexpr (C::template contains<aod::CentFT0Cs>()) {
          QAregistry.fill(HIST("Tracks/Centrality/hPtPhiEtaZvtxEffBestRec"),
                          particle.pt(), itrack.phi(), itrack.eta(),
                          mcCollision.posZ(), c_rec);
        } else {
          QAregistry.fill(HIST("Tracks/hPtPhiEtaZvtxEffBestRec"), particle.pt(),
                          itrack.phi(), itrack.eta(), mcCollision.posZ());
        }
      } else {
        if constexpr (C::template contains<aod::CentFT0Cs>()) {
          QAregistry.fill(HIST("Tracks/Centrality/hPtEffBestFakeRec"),
                          itrack.pt(), c_rec);
        } else {
          QAregistry.fill(HIST("Tracks/hPtEffBestFakeRec"), itrack.pt());
        }
      }
    }
  }

  /// @brief process function to calculate tracking efficiency (inclusive, based
  /// on BestCollisionsFwd)
  void processTrkEffBestInclusive(
    soa::Filtered<soa::Join<Colls, aod::McCollisionLabels>>::iterator const& collision,
    aod::McCollisions const& mccollisions, filtParticles const& particles,
    filtMcMftTracks const& tracks,
    soa::SmallGroups<aod::BestCollisionsFwd> const& besttracks)
  {
    processTrkEffBest<Colls, aod::McCollisions>(collision, mccollisions,
                                                particles, tracks, besttracks);
  }

  PROCESS_SWITCH(
    PseudorapidityDensityMFT, processTrkEffBestInclusive,
    "Process tracking efficiency (inclusive, based on BestCollisionsFwd)",
    false);

  /// @brief process function to calculate tracking efficiency (in FT0 bins,
  /// based on BestCollisionsFwd)
  void processTrkEffBestCent(
    soa::Filtered<soa::Join<CollsCent, aod::McCollisionLabels>>::
      iterator const& collision,
    aod::McCollisions const& mccollisions, filtParticles const& particles,
    filtMcMftTracks const& tracks,
    soa::SmallGroups<aod::BestCollisionsFwd> const& besttracks)
  {
    processTrkEffBest<CollsCent, aod::McCollisions>(
      collision, mccollisions, particles, tracks, besttracks);
  }

  PROCESS_SWITCH(
    PseudorapidityDensityMFT, processTrkEffBestCent,
    "Process tracking efficiency (in FT0 bins, based on BestCollisionsFwd)",
    false);

  Preslice<filtMftTracks> filtTrkperCol = o2::aod::fwdtrack::collisionId;

  /// @brief process template function for MC QA checks
  template <typename C>
  void processMcQA(
    typename soa::SmallGroups<soa::Join<C, aod::McCollisionLabels>> const& collisions,
    aod::McCollisions const& mcCollisions,
    filtParticles const& /*particles*/, MFTTracksLabeled const& tracks,
    aod::AmbiguousMFTTracks const& atracks)
  {
    for (const auto& collision : collisions) {
      float c_rec = -1;
      if constexpr (C::template contains<aod::CentFT0Cs>()) {
        c_rec = collision.centFT0C();
        QAregistry.fill(
          HIST("Events/Centrality/hRecPerGenColls"),
          static_cast<float>(collisions.size()) / mcCollisions.size(), c_rec);
      } else {
        QAregistry.fill(HIST("Events/hRecPerGenColls"),
                        static_cast<float>(collisions.size()) /
                          mcCollisions.size());
      }

      if (!isGoodEvent<false>(collision)) {
        return;
      }

      auto trkPerColl = tracks.sliceBy(filtTrkperCol, collision.globalIndex());
      uint Ntracks{0u}, Natracks{0u};
      for (const auto& track : trkPerColl) {
        Ntracks++;
        for (const auto& atrack : atracks) {
          if (atrack.mfttrackId() == track.globalIndex()) {
            Natracks++;
            break;
          }
        }
      }
      if constexpr (C::template contains<aod::CentFT0Cs>()) {
        QAregistry.fill(HIST("Tracks/Centrality/hNmftTrks"), Ntracks, c_rec);
        QAregistry.fill(HIST("Tracks/Centrality/hFracAmbiguousMftTrks"),
                        static_cast<float>(Natracks) / Ntracks, c_rec);
      } else {
        QAregistry.fill(HIST("Tracks/hNmftTrks"), Ntracks);
        QAregistry.fill(HIST("Tracks/hFracAmbiguousMftTrks"),
                        static_cast<float>(Natracks) / Ntracks);
      }
    }
  }

  /// @brief process function for QA checks (inclusive)
  void processMcQAInclusive(
    soa::SmallGroups<soa::Join<Colls, aod::McCollisionLabels>> const& collisions,
    aod::McCollisions const& mcCollisions, filtParticles const& particles,
    MFTTracksLabeled const& tracks, aod::AmbiguousMFTTracks const& atracks)
  {
    processMcQA<Colls>(collisions, mcCollisions, particles, tracks, atracks);
  }

  PROCESS_SWITCH(PseudorapidityDensityMFT, processMcQAInclusive,
                 "Process MC QA checks (inclusive)", false);

  /// @brief process function for QA checks (in FT0 bins)
  void processMcQACent(
    soa::SmallGroups<soa::Join<CollsCent, aod::McCollisionLabels>> const& collisions,
    aod::McCollisions const& mcCollisions, filtParticles const& particles,
    MFTTracksLabeled const& tracks, aod::AmbiguousMFTTracks const& atracks)
  {
    processMcQA<CollsCent>(collisions, mcCollisions, particles, tracks,
                           atracks);
  }

  PROCESS_SWITCH(PseudorapidityDensityMFT, processMcQACent,
                 "Process MC QA checks (in FT0 bins)", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<PseudorapidityDensityMFT>(cfgc)};
}
