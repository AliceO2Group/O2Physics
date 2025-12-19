// Copyright 2020-2025 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

// \file   dndeta-mft.cxx
// \author Sarah Herrmann <sarah.herrmann@cern.ch>
//
// \brief This code loops over MFT tracks and collisions and fills histograms
//        useful to compute dNdeta

#include "PWGMM/Mult/DataModel/bestCollisionTable.h"

#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CommonConstants/MathConstants.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/Configurable.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/RuntimeError.h"
#include "Framework/runDataProcessing.h"
#include "MathUtils/Utils.h"
#include "ReconstructionDataFormats/GlobalTrackID.h"

#include "TFile.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <unordered_set>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::track;

AxisSpec PtAxis = {1001, -0.005, 10.005};
AxisSpec DeltaZAxis = {61, -6.1, 6.1};
AxisSpec ZAxis = {301, -30.1, 30.1};
AxisSpec PhiAxis = {629, 0, o2::constants::math::TwoPI, "Rad", "phi axis"};
// AxisSpec EtaAxis = {18, -4.6, -1.};
AxisSpec DCAxyAxis = {5000, -1, 500};
AxisSpec DCAzAxis = {5000, -251, 250};
AxisSpec CentAxis = {{0, 10, 20, 30, 40, 50, 60, 70, 80, 100}};

static constexpr TrackSelectionFlags::flagtype trackSelectionITS =
  TrackSelectionFlags::kITSNCls | TrackSelectionFlags::kITSChi2NDF |
  TrackSelectionFlags::kITSHits;

static constexpr TrackSelectionFlags::flagtype trackSelectionTPC =
  TrackSelectionFlags::kTPCNCls |
  TrackSelectionFlags::kTPCCrossedRowsOverNCls |
  TrackSelectionFlags::kTPCChi2NDF;

static constexpr TrackSelectionFlags::flagtype trackSelectionDCA =
  TrackSelectionFlags::kDCAz | TrackSelectionFlags::kDCAxy;

using MFTTracksLabeled = soa::Join<o2::aod::MFTTracks, aod::McMFTTrackLabels>;

struct PseudorapidityDensityMFT {
  SliceCache cache;
  Preslice<aod::MFTTracks> perCol = o2::aod::fwdtrack::collisionId;
  Preslice<aod::McParticles> perMcCol = aod::mcparticle::mcCollisionId;
  Preslice<aod::Tracks> perColCentral = aod::track::collisionId;

  Service<o2::framework::O2DatabasePDG> pdg;

  Configurable<float> estimatorEta{"estimatorEta", 1.0,
                                   "eta range for INEL>0 sample definition"};

  Configurable<bool> useEvSel{"useEvSel", true, "use event selection"};
  Configurable<bool> disableITSROFCut{"disableITSROFCut", false, "Disable ITS ROF cut for event selection"};
  ConfigurableAxis multBinning{"multBinning", {701, -0.5, 700.5}, ""};
  ConfigurableAxis EtaAxis = {"etaBinning", {36, -4.6, -1.}, ""};

  Configurable<bool> useZDiffCut{"useZDiffCut", true, "use Z difference cut"};
  Configurable<float> maxZDiff{
    "maxZDiff", 1.0f,
    "max allowed Z difference for reconstructed collisions (cm)"};

  Configurable<bool> usePhiCut{"usePhiCut", true, "use azimuthal angle cut"};
  Configurable<bool> useDCAxyCut{"useDCAxyCut", false, "use DCAxy cut"};
  Configurable<bool> useDCAzCut{"useDCAzCut", false, "use DCAz cut"};

  Configurable<float> cfgPhiCut{"cfgPhiCut", 0.1f,
                                "Cut on azimuthal angle of MFT tracks"};
  Configurable<float> cfgPhiCut1{"cfgPhiCut1", 0.0f,
                                 "low Cut on azimuthal angle of MFT tracks"};
  Configurable<float> cfgPhiCut2{"cfgPhiCut2", 6.3f,
                                 "high Cut on azimuthal angle of MFT tracks"};
  Configurable<float> cfgVzCut1{"cfgVzCut1", -30.0f,
                                "Cut1 on vertex position of MFT tracks"};
  Configurable<float> cfgVzCut2{"cfgVzCut2", 30.0f,
                                "Cut2 on vertex position of MFT tracks"};
  Configurable<float> cfgnCluster{"cfgnCluster", 5.0f,
                                  "Cut on no of clusters per MFT track"};
  Configurable<float> cfgnEta1{"cfgnEta1", -4.5f,
                               "Cut on eta1"};
  Configurable<float> cfgnEta2{"cfgnEta2", -1.0f,
                               "Cut on eta1"};
  Configurable<float> cfgChi2NDFMax{"cfgChi2NDFMax", 2000.0f, "Max allowed chi2/NDF for MFT tracks"};
  Configurable<float> maxDCAxy{"maxDCAxy", 2.0f, "Cut on dcaXY"};
  Configurable<float> maxDCAz{"maxDCAz", 2.0f, "Cut on dcaZ"};

  HistogramRegistry registry{
    "registry",
    {{"TracksEtaZvtx",
      "; #eta; #it{z}_{vtx} (cm); tracks",
      {HistType::kTH2F, {EtaAxis, ZAxis}}}, //
     {"Tracks/EtaZvtx_gt0",
      "; #eta; #it{z}_{vtx} (cm); tracks",
      {HistType::kTH2F, {EtaAxis, ZAxis}}}, //
     {"TracksPhiEta",
      "; #varphi; #eta; tracks",
      {HistType::kTH2F, {PhiAxis, EtaAxis}}}, //
     {"TracksPhiZvtx",
      "; #varphi; #it{z}_{vtx} (cm); tracks",
      {HistType::kTH2F, {PhiAxis, ZAxis}}}, //
     {"TracksPtEta",
      " ; p_{T} (GeV/c); #eta",
      {HistType::kTH2F, {PtAxis, EtaAxis}}}, //
     {"EventSelection",
      ";status;events",
      {HistType::kTH1F, {{15, 0.5, 15.5}}}},
     {"EventCounts",
      ";status;events",
      {HistType::kTH1F, {{2, 0.5, 2.5}}}},
     {"Tracks/Control/TrackCount", ";status;Track counts", {HistType::kTH1F, {{15, 0.5, 15.5}}}}, // added
     // Purity-related histograms
     {"Purity/SelectedAfterDCAxy/All",
      ";bin;counts",
      {HistType::kTH1F, {{1, 0.5, 1.5}}}},
     {"Purity/SelectedAfterDCAxy/AllEta",
      ";#eta;counts",
      {HistType::kTH1F, {EtaAxis}}},
     {"Purity/Gen/PrimaryEta",
      ";#eta;primaries",
      {HistType::kTH1F, {EtaAxis}}},
     {"Purity/Gen/All",
      ";bin;counts",
      {HistType::kTH1F, {{1, 0.5, 1.5}}}},
     {"Purity/Gen/AllEta",
      ";#eta;counts",
      {HistType::kTH1F, {EtaAxis}}}}};

  void init(InitContext&)
  {
    if (static_cast<int>(doprocessMult) +
          static_cast<int>(doprocessMultReassoc) +
          static_cast<int>(doprocessMultReassoc3d) +
          static_cast<int>(doprocessCountingCentrality) >
        1) {
      LOGP(fatal,
           "Exactly one process function between processMult, "
           "processMultReassoc, processMultReassoc3d and processCountingCentrality should be "
           "enabled!");
    }
    AxisSpec MultAxis = {multBinning, "N_{trk}"};
    auto hstat = registry.get<TH1>(HIST("EventSelection"));
    auto* x = hstat->GetXaxis();
    x->SetBinLabel(1, "All");
    x->SetBinLabel(2, "Vz");
    x->SetBinLabel(3, "Vz+ITSRof");
    x->SetBinLabel(4, "Vz+Selected");
    x->SetBinLabel(5, "Sel8+Vz+INEL>0");
    x->SetBinLabel(6, "Sel INEL,INEL_fwd>0");
    x->SetBinLabel(7, "Rejected");
    x->SetBinLabel(8, "Good BCs");
    x->SetBinLabel(9, "BCs with collisions");
    x->SetBinLabel(10, "BCs with pile-up/splitting");
    x->SetBinLabel(11, "percollisionSample>0");
    x->SetBinLabel(12, "midtracks+percollisionSample>0");
    registry.add({"EventsNtrkZvtx",
                  "; N_{trk}; #it{z}_{vtx} (cm); events",
                  {HistType::kTH2F, {MultAxis, ZAxis}}});
    registry.add({"EventsNtrkZvtx_gt0",
                  "; N_{trk}; #it{z}_{vtx} (cm); events",
                  {HistType::kTH2F, {MultAxis, ZAxis}}});
    registry.add({"Tracks/2Danalysis/EventsNtrkZvtx_all",
                  "; N_{trk}; #it{z}_{vtx} (cm); events",
                  {HistType::kTH2F, {MultAxis, ZAxis}}});
    registry.add({"Tracks/2Danalysis/EventsNtrkZvtx_sel8",
                  "; N_{trk}; #it{z}_{vtx} (cm); events",
                  {HistType::kTH2F, {MultAxis, ZAxis}}});
    registry.add({"Tracks/2Danalysis/EventsNtrkZvtx_sel8_inelgt0",
                  "; N_{trk}; #it{z}_{vtx} (cm); events",
                  {HistType::kTH2F, {MultAxis, ZAxis}}});
    registry.add({"Tracks/2Danalysis/EventsNtrkZvtx_sel8_inelfwdgt0",
                  "; N_{trk}; #it{z}_{vtx} (cm); events",
                  {HistType::kTH2F, {MultAxis, ZAxis}}});
    registry.add({"Tracks/Control/DCAXY",
                  " ; DCA_{XY} (cm)",
                  {HistType::kTH1F, {DCAxyAxis}}});
    if (doprocessGen) {
      registry.add({"EventsNtrkZvtxGen",
                    "; N_{trk}; #it{z}_{vtx} (cm); events",
                    {HistType::kTH2F, {MultAxis, ZAxis}}});
      registry.add({"EventsNtrkZvtxGen_t",
                    "; N_{trk}; #it{z}_{vtx} (cm); events",
                    {HistType::kTH2F, {MultAxis, ZAxis}}});
      registry.add({"EventsNtrkZvtxGen_gt0",
                    "; N_{trk}; #it{z}_{vtx} (cm); events",
                    {HistType::kTH2F, {MultAxis, ZAxis}}});
      registry.add({"EventsNtrkZvtxGen_gt0t",
                    "; N_{trk}; #it{z}_{vtx} (cm); events",
                    {HistType::kTH2F, {MultAxis, ZAxis}}});
      registry.add({"TracksEtaZvtxGen",
                    "; #eta; #it{z}_{vtx} (cm); tracks",
                    {HistType::kTH2F, {EtaAxis, ZAxis}}});
      registry.add({"TracksEtaZvtxGen_t",
                    "; #eta; #it{z}_{vtx} (cm); tracks",
                    {HistType::kTH2F, {EtaAxis, ZAxis}}});
      registry.add({"TracksEtaZvtxGen_gt0",
                    "; #eta; #it{z}_{vtx} (cm); tracks",
                    {HistType::kTH2F, {EtaAxis, ZAxis}}});
      registry.add({"TracksEtaZvtxGen_gt0t",
                    "; #eta; #it{z}_{vtx} (cm); tracks",
                    {HistType::kTH2F, {EtaAxis, ZAxis}}});
      registry.add({"TracksPhiEtaGen",
                    "; #varphi; #eta; tracks",
                    {HistType::kTH2F, {PhiAxis, EtaAxis}}});
      registry.add({"TracksPhiEtaGen_gt0",
                    "; #varphi; #eta; tracks",
                    {HistType::kTH2F, {PhiAxis, EtaAxis}}});
      registry.add({"TracksPhiEtaGen_gt0t",
                    "; #varphi; #eta; tracks",
                    {HistType::kTH2F, {PhiAxis, EtaAxis}}});
      registry.add({"TracksPhiZvtxGen",
                    "; #varphi; #it{z}_{vtx} (cm); tracks",
                    {HistType::kTH2F, {PhiAxis, ZAxis}}}); //
      registry.add({"TracksToPartPtEta",
                    " ; p_{T} (GeV/c); #eta",
                    {HistType::kTH2F, {PtAxis, EtaAxis}}}); //
      registry.add({"TracksPtEtaGen",
                    " ; p_{T} (GeV/c); #eta",
                    {HistType::kTH2F, {PtAxis, EtaAxis}}});
      registry.add({"TracksPtEtaGen_t",
                    " ; p_{T} (GeV/c); #eta",
                    {HistType::kTH2F, {PtAxis, EtaAxis}}});
      registry.add({"EventEfficiency",
                    "; status; events",
                    {HistType::kTH1F, {{5, 0.5, 5.5}}}});
      registry.add({"NotFoundEventZvtx",
                    " ; #it{z}_{vtx} (cm)",
                    {HistType::kTH1F, {ZAxis}}});
      registry.add({"EventsZposDiff",
                    " ; Z_{rec} - Z_{gen} (cm)",
                    {HistType::kTH1F, {DeltaZAxis}}});
      registry.add({"EventsSplitMult", " ; N_{gen}", {HistType::kTH1F, {MultAxis}}});
      auto heff = registry.get<TH1>(HIST("EventEfficiency"));
      x = heff->GetXaxis();
      x->SetBinLabel(1, "Generated");
      x->SetBinLabel(2, "Generated INEL>0");
      x->SetBinLabel(3, "Reconstructed");
      x->SetBinLabel(4, "Selected");
      x->SetBinLabel(5, "Selected INEL>0");
    }

    if (doprocessMultReassoc || doprocessMultReassoc3d) {
      registry.add({"Tracks/Control/DeltaZ",
                    " ; #it{z_{orig}}-#it{z_{reass}}",
                    {HistType::kTH1F, {ZAxis}}});

      registry.add({"Tracks/Control/TrackAmbDegree",
                    " ; N_{coll}^{comp}",
                    {HistType::kTH1F, {{51, -0.5, 50.5}}}});
      registry.add({"Tracks/Control/TrackIsAmb",
                    " ; isAmbiguous",
                    {HistType::kTH1I, {{2, -0.5, 1.5}}}});

      auto htrk = registry.get<TH1>(HIST("Tracks/Control/TrackCount"));
      auto* x = htrk->GetXaxis();
      x->SetBinLabel(0, "All");
      x->SetBinLabel(1, "Reass");
      x->SetBinLabel(2, "Not Reass");
      x->SetBinLabel(3, "Amb");
      x->SetBinLabel(4, "Amb+Not-reass");
      x->SetBinLabel(5, "Non-Amb");
      x->SetBinLabel(6, "Not-Reass+Non-Amb");
      x->SetBinLabel(7, "Amb+Non-Amb");
      x->SetBinLabel(8, "colid<0");
      x->SetBinLabel(9, "wo orphan");

      registry.add({"Tracks/Control/ReassignedTracksEtaZvtx",
                    "; #eta; #it{z}_{vtx} (cm); tracks",
                    {HistType::kTH2F, {EtaAxis, ZAxis}}});
      registry.add({"Tracks/Control/ReassignedTracksPhiEta",
                    "; #varphi; #eta; tracks",
                    {HistType::kTH2F, {PhiAxis, EtaAxis}}});
      registry.add({"Tracks/Control/ReassignedVertexCorr",
                    "; #it{z}_{vtx}^{orig} (cm); #it{z}_{vtx}^{re} (cm)",
                    {HistType::kTH2F, {ZAxis, ZAxis}}});

      registry.add({"Tracks/Control/notReassignedTracksEtaZvtx",
                    "; #eta; #it{z}_{vtx} (cm); tracks",
                    {HistType::kTH2F, {EtaAxis, ZAxis}}});
      registry.add({"Tracks/Control/notReassignedTracksPhiEta",
                    "; #varphi; #eta; tracks",
                    {HistType::kTH2F, {PhiAxis, EtaAxis}}});
      registry.add({"Tracks/Control/notReassignedVertexCorr",
                    "; #it{z}_{vtx}^{orig} (cm); #it{z}_{vtx}^{re} (cm)",
                    {HistType::kTH2F, {ZAxis, ZAxis}}});
      registry.add({"Tracks/Control/Chi2NDF",
                    " ; #chi^{2}/ndf",
                    {HistType::kTH1F, {{5000, 0.0, 5000.0}}}});
      registry.add({"Tracks/Control/amb/AmbTracksEtaZvtx",
                    "; #eta; #it{z}_{vtx} (cm); tracks",
                    {HistType::kTH2F, {EtaAxis, ZAxis}}}); //

      registry.add({"Tracks/Control/woOrp/nTrk",
                    " ; N_{Trk}^{all}",
                    {HistType::kTH1F, {{701, -0.5, 700.5}}}}); //
      registry.add({"Tracks/Control/amb/nTrkAmb",
                    " ; N_{Trk}^{amb}",
                    {HistType::kTH1F, {{701, -0.5, 700.5}}}}); //
      registry.add({"Tracks/Control/nonamb/nTrkNonAmb",
                    " ; N_{Trk}^{nonamb}",
                    {HistType::kTH1F, {{701, -0.5, 700.5}}}}); //

      registry.add({"Tracks/Control/amb/AmbTracksPhiEta",
                    "; #varphi; #eta; tracks",
                    {HistType::kTH2F, {PhiAxis, EtaAxis}}}); //
      registry.add({"Tracks/Control/amb/AmbVertexCorr",
                    "; #it{z}_{vtx}^{orig} (cm); #it{z}_{vtx}^{re} (cm)",
                    {HistType::kTH2F, {ZAxis, ZAxis}}}); //
      registry.add({"Tracks/Control/amb/EtaZvtxAmb_gt0",
                    "; #eta; #it{z}_{vtx} (cm); tracks",
                    {HistType::kTH2F, {EtaAxis, ZAxis}}}); //
      registry.add({"Tracks/Control/amb/DCAxy_amb", " ; DCA_{xy} (cm) ambiguous",
                    //  {HistType::kTH1F,{{100000, 0.5, 100000.0}}}}); //
                    {HistType::kTH1F, {DCAxyAxis}}}); //

      registry.add({"Tracks/Control/nonamb/nonAmbTracksEtaZvtx",
                    "; #eta; #it{z}_{vtx} (cm); tracks",
                    {HistType::kTH2F, {EtaAxis, ZAxis}}}); //

      registry.add({"Tracks/Control/nonamb/nonAmbTracksPhiEta",
                    "; #varphi; #eta; tracks",
                    {HistType::kTH2F, {PhiAxis, EtaAxis}}}); //
      registry.add({"Tracks/Control/nonamb/nonAmbVertexCorr",
                    "; #it{z}_{vtx}^{orig} (cm); #it{z}_{vtx}^{re} (cm)",
                    {HistType::kTH2F, {ZAxis, ZAxis}}}); //
      registry.add({"Tracks/Control/nonamb/EtaZvtxNonAmb_gt0",
                    "; #eta; #it{z}_{vtx} (cm); tracks",
                    {HistType::kTH2F, {EtaAxis, ZAxis}}}); //
      registry.add({"Tracks/Control/nonamb/DCAxy_nonamb", " ; DCA_{xy}(cm) non-ambiguous",
                    //  {HistType::kTH1F,{{100000, 0.5, 100000.0}}}}); //
                    {HistType::kTH1F, {{DCAxyAxis}}}}); //

      registry.add({"Tracks/Control/woOrp/woOrpTracksEtaZvtx",
                    "; #eta; #it{z}_{vtx} (cm); tracks",
                    {HistType::kTH2F, {EtaAxis, ZAxis}}}); //
      registry.add({"Tracks/Control/woOrp/woOrpEtaZvtx_gt0",
                    "; #eta; #it{z}_{vtx} (cm); tracks",
                    {HistType::kTH2F, {EtaAxis, ZAxis}}}); //
      registry.add({"Tracks/2Danalysis/EtaZvtx",
                    "; #eta; #it{z}_{vtx} (cm); tracks",
                    {HistType::kTH2F, {EtaAxis, ZAxis}}}); //
      registry.add({"Tracks/2Danalysis/EtaZvtx_sel8",
                    "; #eta; #it{z}_{vtx} (cm); tracks",
                    {HistType::kTH2F, {EtaAxis, ZAxis}}}); //
      registry.add({"Tracks/2Danalysis/EtaZvtx_sel8_inelgt0",
                    "; #eta; #it{z}_{vtx} (cm); tracks",
                    {HistType::kTH2F, {EtaAxis, ZAxis}}}); //
      registry.add({"Tracks/2Danalysis/EtaZvtx_sel8_inelfwdgt0",
                    "; #eta; #it{z}_{vtx} (cm); tracks",
                    {HistType::kTH2F, {EtaAxis, ZAxis}}}); //
      registry.add({"Tracks/Control/woOrp/woOrpTracksPhiEta",
                    "; #varphi; #eta; tracks",
                    {HistType::kTH2F, {PhiAxis, EtaAxis}}}); //
      registry.add({"Tracks/Control/woOrp/woOrpVertexCorr",
                    "; #it{z}_{vtx}^{orig} (cm); #it{z}_{vtx}^{re} (cm)",
                    {HistType::kTH2F, {ZAxis, ZAxis}}}); //
      registry.add({"Tracks/Control/woOrp/DCAxy_woOrp", " ; DCA_{xy}(cm) w/o orphan",
                    //  {HistType::kTH1F,{{100000, 0.5, 100000.0}}}}); //
                    {HistType::kTH1F, {{DCAxyAxis}}}}); //

      if (doprocessMultReassoc3d) {
        // DCAz histograms analogous to DCAxy, only for 3D reassociation
        registry.add({"Tracks/Control/DCAZ",
                      " ; DCA_{Z} (cm)",
                      {HistType::kTH1F, {DCAzAxis}}});
        registry.add({"Tracks/Control/amb/DCAz_amb",
                      " ; DCA_{z} (cm) ambiguous",
                      {HistType::kTH1F, {DCAzAxis}}});
        registry.add({"Tracks/Control/nonamb/DCAz_nonamb",
                      " ; DCA_{z}(cm) non-ambiguous",
                      {HistType::kTH1F, {DCAzAxis}}});
        registry.add({"Tracks/Control/woOrp/DCAz_woOrp",
                      " ; DCA_{z}(cm) w/o orphan",
                      {HistType::kTH1F, {DCAzAxis}}});
      }

      registry.add({"collisionID", " ; Collision ID",
                    //  {HistType::kTH1F,{{100000, 0.5, 100000.0}}}}); //
                    {HistType::kTH1F, {{100000, -50000.0, 50000.0}}}}); //
      registry.add({"collisionIDamb", " ; Collision ID amb",
                    //  {HistType::kTH1F,{{100000, 0.5, 100000.0}}}}); //
                    {HistType::kTH1F, {{100000, -50000.0, 50000.0}}}});                                                                                        //
      registry.add({"NonambEventCounts", " ; EventCounts Nonamb", {HistType::kTH1F, {{1, 0.5, 1.5}}}});                                                        //
      registry.add({"hNumCollisionsNonAmb_InelMFT", " ; Number of Collisions with Non-Ambiguous Tracks;Count;Frequency", {HistType::kTH1F, {{1, 0.5, 1.5}}}}); //
      registry.add({"hNumCollisionsAmb_InelMFT", " ; Number of Collisions with Non-Ambiguous Tracks;Count;Frequency", {HistType::kTH1F, {{1, 0.5, 1.5}}}});    //
      registry.add({"hNumCollisions_InelMFT", " ; Number of selected events with Inel>0 and MFT>0;Count;Frequency", {HistType::kTH1F, {{1, 0.5, 1.5}}}});      //
      registry.add({"hNumCollisions_Inel", " ; Number of selected events with Inel>0;Count;Frequency", {HistType::kTH1F, {{1, 0.5, 1.5}}}});                   //
      registry.add({"ambEventCounts", " ; EventCounts Nonamb", {HistType::kTH1F, {{1, 0.5, 1.5}}}});                                                           //
    }

    if (doprocessCountingCentrality) {
      registry.add({"Events/Centrality/Selection",
                    ";status;centrality;events",
                    {HistType::kTH2F, {{3, 0.5, 3.5}, CentAxis}}});
      auto hstat = registry.get<TH2>(HIST("Events/Centrality/Selection"));
      auto* x = hstat->GetXaxis();
      x->SetBinLabel(1, "All");
      x->SetBinLabel(2, "Selected");
      x->SetBinLabel(3, "Rejected");

      registry.add({"Events/Centrality/NtrkZvtx",
                    "; N_{trk}; Z_{vtx} (cm); centrality",
                    {HistType::kTH3F, {MultAxis, ZAxis, CentAxis}}});
      registry.add({"Tracks/Centrality/EtaZvtx",
                    "; #eta; Z_{vtx} (cm); centrality",
                    {HistType::kTH3F, {EtaAxis, ZAxis, CentAxis}}});
      registry.add({"Tracks/Centrality/PhiEta",
                    "; #varphi; #eta; centrality",
                    {HistType::kTH3F, {PhiAxis, EtaAxis, CentAxis}}});
      registry.add({"Tracks/Centrality/Control/PtEta",
                    " ; p_{T} (GeV/c); #eta; centrality",
                    {HistType::kTH3F, {PtAxis, EtaAxis, CentAxis}}});
      registry.add({"Tracks/Centrality/Control/DCAXYPt",
                    " ; p_{T} (GeV/c) ; DCA_{XY} (cm); centrality",
                    {HistType::kTH3F, {PtAxis, DCAxyAxis, CentAxis}}});
      registry.add({"Tracks/Centrality/Control/ReassignedDCAXYPt",
                    " ; p_{T} (GeV/c) ; DCA_{XY} (cm); centrality",
                    {HistType::kTH3F, {PtAxis, DCAxyAxis, CentAxis}}});
      registry.add({"Tracks/Centrality/Control/ExtraDCAXYPt",
                    " ; p_{T} (GeV/c) ; DCA_{XY} (cm); centrality",
                    {HistType::kTH3F, {PtAxis, DCAxyAxis, CentAxis}}});
      registry.add({"Tracks/Centrality/Control/ExtraTracksEtaZvtx",
                    "; #eta; Z_{vtx} (cm); centrality",
                    {HistType::kTH3F, {EtaAxis, ZAxis, CentAxis}}});
      registry.add({"Tracks/Centrality/Control/ExtraTracksPhiEta",
                    "; #varphi; #eta; centrality",
                    {HistType::kTH3F, {PhiAxis, EtaAxis, CentAxis}}});
      registry.add({"Tracks/Centrality/Control/ReassignedTracksEtaZvtx",
                    "; #eta; Z_{vtx} (cm); centrality",
                    {HistType::kTH3F, {EtaAxis, ZAxis, CentAxis}}});
      registry.add({"Tracks/Centrality/Control/ReassignedTracksPhiEta",
                    "; #varphi; #eta; centrality",
                    {HistType::kTH3F, {PhiAxis, EtaAxis, CentAxis}}});
      registry.add({"Tracks/Centrality/Control/ReassignedVertexCorr",
                    "; Z_{vtx}^{orig} (cm); Z_{vtx}^{re} (cm); centrality",
                    {HistType::kTH3F, {ZAxis, ZAxis, CentAxis}}});
    }

    if (doprocessGenCent) {
      registry.add({"Events/Centrality/EventEfficiency",
                    ";status;centrality;events",
                    {HistType::kTH2F, {{2, 0.5, 2.5}, CentAxis}}});
      auto heff = registry.get<TH2>(HIST("Events/Centrality/EventEfficiency"));
      auto* x = heff->GetXaxis();
      x->SetBinLabel(1, "Generated");
      x->SetBinLabel(2, "Selected");

      registry.add("Events/Centrality/CentPercentileMCGen",
                   "CentPercentileMCGen", kTH1D, {CentAxis}, false);
      registry.add({"Events/Centrality/NtrkZvtxGen",
                    "; N_{trk}; Z_{vtx} (cm); centrality",
                    {HistType::kTH3F, {MultAxis, ZAxis, CentAxis}}});
      registry.add({"Events/Centrality/NtrkZvtxGen_t",
                    "; N_{trk}; Z_{vtx} (cm); centrality",
                    {HistType::kTH3F, {MultAxis, ZAxis, CentAxis}}});
      registry.add({"Tracks/Centrality/EtaZvtxGen_t",
                    "; #eta; Z_{vtx} (cm); centrality",
                    {HistType::kTH3F, {EtaAxis, ZAxis, CentAxis}}});
      registry.add({"Tracks/Centrality/EtaZvtxGen",
                    "; #eta; Z_{vtx} (cm); centrality",
                    {HistType::kTH3F, {EtaAxis, ZAxis, CentAxis}}});
      registry.add({"Tracks/Centrality/PhiEtaGen",
                    "; #varphi; #eta; centrality",
                    {HistType::kTH3F, {PhiAxis, EtaAxis, CentAxis}}});
    }
  }

  using FullBCs = soa::Join<aod::BCsWithTimestamps, aod::BcSels>;
  void processTagging(FullBCs const& bcs,
                      soa::Join<aod::Collisions, aod::EvSels> const& collisions)
  {

    std::vector<typename std::decay_t<decltype(collisions)>::iterator> cols;
    for (const auto& bc : bcs) {
      if (!useEvSel ||
          (useEvSel && ((bc.selection_bit(aod::evsel::kIsBBT0A) &&
                         bc.selection_bit(aod::evsel::kIsBBT0C)) != 0))) {
        registry.fill(HIST("EventSelection"), 8); // added 5->12
        cols.clear();
        for (const auto& collision : collisions) {
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
          registry.fill(HIST("EventSelection"), 9); // added 6->13
          if (cols.size() > 1) {
            registry.fill(HIST("EventSelection"), 10); // added 7->14
          }
        }
      }
    }
  }

  PROCESS_SWITCH(PseudorapidityDensityMFT, processTagging,
                 "Collect event sample stats", true);

  Partition<aod::MFTTracks> sample =
    (aod::fwdtrack::eta < -2.8f) && (aod::fwdtrack::eta > -3.2f);

  Partition<aod::Tracks> sampleCentral = (nabs(aod::track::eta) < 1.f);

  expressions::Filter atrackFilter =
    (aod::fwdtrack::bestCollisionId >= 0) && (aod::fwdtrack::eta < -2.0f) &&
    (aod::fwdtrack::eta > -3.9f) && (nabs(aod::fwdtrack::bestDCAXY) <= 2.f);

  using CollwEv = soa::Join<aod::Collisions, aod::EvSels>;

  expressions::Filter trackSelectionCentral =
    ((aod::track::trackCutFlag & trackSelectionITS) == trackSelectionITS) &&
    ifnode((aod::track::v001::detectorMap & (uint8_t)o2::aod::track::TPC) ==
             (uint8_t)o2::aod::track::TPC,
           (aod::track::trackCutFlag & trackSelectionTPC) ==
             trackSelectionTPC,
           true) &&
    ((aod::track::trackCutFlag & trackSelectionDCA) == trackSelectionDCA) &&
    (nabs(aod::track::eta) < estimatorEta);

  using FiCentralTracks = soa::Filtered<
    soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection,
              aod::TracksDCA>>; // central tracks for INEL>0

  void processMult(CollwEv::iterator const& collision,
                   aod::MFTTracks const& tracks,
                   FiCentralTracks const& midtracks, aod::Tracks const&)
  {

    registry.fill(HIST("EventSelection"), 1.);
    if (!useEvSel || (useEvSel && collision.sel8())) {
      registry.fill(HIST("EventSelection"), 2.);
      auto z = collision.posZ();
      auto perCollisionSample = sampleCentral->sliceByCached(
        o2::aod::track::collisionId, collision.globalIndex(), cache);
      auto Ntrk = perCollisionSample.size();

      registry.fill(HIST("EventsNtrkZvtx"), Ntrk, z);

      if (midtracks.size() > 0) // INEL>0
      {
        registry.fill(HIST("EventSelection"), 3.);
        registry.fill(HIST("EventsNtrkZvtx_gt0"), Ntrk, z);
      }

      if (tracks.size() > 0) {
        for (const auto& track : tracks) {

          float phi = track.phi();
          o2::math_utils::bringTo02Pi(phi);

          if (usePhiCut) {
            if ((phi < cfgPhiCut) ||
                ((phi > o2::constants::math::PI - cfgPhiCut) && (phi < o2::constants::math::PI + cfgPhiCut)) ||
                (phi > o2::constants::math::TwoPI - cfgPhiCut) ||
                ((phi > ((o2::constants::math::PIHalf - 0.1) * o2::constants::math::PI) - cfgPhiCut) &&
                 (phi < ((o2::constants::math::PIHalf - 0.1) * o2::constants::math::PI) + cfgPhiCut)))
              continue;
          }

          registry.fill(HIST("TracksEtaZvtx"), track.eta(), z);
          if (midtracks.size() > 0) // INEL>0
          {
            registry.fill(HIST("Tracks/EtaZvtx_gt0"), track.eta(), z);
          }
          registry.fill(HIST("TracksPhiEta"), phi, track.eta());
          registry.fill(HIST("TracksPtEta"), track.pt(), track.eta());
          if ((track.eta() < -2.0f) && (track.eta() > -3.9f)) {
            registry.fill(HIST("TracksPhiZvtx"), phi, z);
          }
        }
      }

    } else {
      registry.fill(HIST("EventSelection"), 4.);
    }
  }

  PROCESS_SWITCH(PseudorapidityDensityMFT, processMult,
                 "Process reco or data info", true);
  // Common implementation for both BestCollisionsFwd and BestCollisionsFwd3d
  template <typename RetracksT>
  void processMultReassocCommon(CollwEv::iterator const& collision,
                                o2::aod::MFTTracks const&,
                                RetracksT const& retracks,
                                FiCentralTracks const& midtracks, aod::Tracks const&)
  {
    registry.fill(HIST("EventSelection"), 1.);
    auto perCollisionSample = sampleCentral->sliceByCached(
      o2::aod::track::collisionId, collision.globalIndex(), cache);
    auto Ntrk = perCollisionSample.size();
    auto z = collision.posZ();
    registry.fill(HIST("EventsNtrkZvtx"), Ntrk, z);
    if ((z >= cfgVzCut1) && (z <= cfgVzCut2)) {
      registry.fill(HIST("Tracks/2Danalysis/EventsNtrkZvtx_all"), Ntrk, z);
      registry.fill(HIST("EventSelection"), 2.);
      for (const auto& retrack : retracks) {
        auto track = retrack.mfttrack();
        float ndf = std::max(2.0f * track.nClusters() - 5.0f, 1.0f);
        float chi2ndf = track.chi2() / ndf;
        float phi = track.phi();
        o2::math_utils::bringTo02Pi(phi);
        if (usePhiCut) {
          if ((phi <= 0.02) || ((phi >= 3.10) && (phi <= 3.23)) || (phi >= 6.21))
            continue;
        }
        float dcaxy_cut = retrack.bestDCAXY();
        if (useDCAxyCut) {
          if (dcaxy_cut > maxDCAxy)
            continue;
        }
        if constexpr (std::is_same_v<RetracksT, soa::SmallGroups<aod::BestCollisionsFwd3d>>) {
          float dcaz_cut = retrack.bestDCAZ();
          if (useDCAzCut) {
            if (dcaz_cut > maxDCAz)
              continue;
          }
        }
        if ((cfgnEta1 < track.eta()) && (track.eta() < cfgnEta2) && track.nClusters() >= cfgnCluster && retrack.ambDegree() > 0 && chi2ndf < cfgChi2NDFMax && (phi > cfgPhiCut1 && phi < cfgPhiCut2)) {
          registry.fill(HIST("Tracks/2Danalysis/EtaZvtx"), track.eta(), z);
        }
      }
      if (!disableITSROFCut && !collision.selection_bit(aod::evsel::kNoITSROFrameBorder)) {
        return;
      }
      registry.fill(HIST("EventSelection"), 3.);
      if (!useEvSel || (useEvSel && collision.selection_bit(aod::evsel::kIsTriggerTVX) && collision.selection_bit(aod::evsel::kNoTimeFrameBorder) && collision.selection_bit(aod::evsel::kNoSameBunchPileup))) {
        registry.fill(HIST("EventSelection"), 4.);
        registry.fill(HIST("Tracks/2Danalysis/EventsNtrkZvtx_sel8"), Ntrk, z);
        std::unordered_set<int> uniqueEvents;
        std::unordered_set<int> uniqueEventsAmb;
        std::unordered_set<int> uniqueCollisions;
        std::unordered_set<int> uniqueCollisionsAmb;
        std::unordered_set<int> eventsInelMFT;
        std::unordered_set<int> eventsInel;
        if (midtracks.size() > 0) {
          registry.fill(HIST("EventSelection"), 5.);
          registry.fill(HIST("EventsNtrkZvtx_gt0"), Ntrk, z);
          registry.fill(HIST("Tracks/2Danalysis/EventsNtrkZvtx_sel8_inelgt0"), Ntrk, z);
          eventsInel.insert(collision.globalIndex());
        }
        if (perCollisionSample.size() > 0) {
          registry.fill(HIST("EventSelection"), 11.);
        }
        if (midtracks.size() > 0 && perCollisionSample.size() > 0) {
          registry.fill(HIST("EventSelection"), 12.);
        }
        int64_t i = 0.0, j = 0.0, k = 0.0;
        for (const auto& retrack : retracks) {
          auto track = retrack.mfttrack();
          float ndf = std::max(2.0f * track.nClusters() - 5.0f, 1.0f);
          float chi2ndf = track.chi2() / ndf;
          float phi = track.phi();
          o2::math_utils::bringTo02Pi(phi);
          if (usePhiCut) {
            if ((phi <= 0.02) || ((phi >= 3.10) && (phi <= 3.23)) || (phi >= 6.21))
              continue;
          }
          float dcaxy_cut = retrack.bestDCAXY();
          if (useDCAxyCut) {
            if (dcaxy_cut > maxDCAxy)
              continue;
          }
          if constexpr (std::is_same_v<RetracksT, soa::SmallGroups<aod::BestCollisionsFwd3d>>) {
            float dcaz_cut = retrack.bestDCAZ();
            if (useDCAzCut) {
              if (dcaz_cut > maxDCAz)
                continue;
            }
          }
          if ((cfgnEta1 < track.eta()) && (track.eta() < cfgnEta2) && track.nClusters() >= cfgnCluster && retrack.ambDegree() > 0 && chi2ndf < cfgChi2NDFMax && (phi > cfgPhiCut1 && phi < cfgPhiCut2)) {
            registry.fill(HIST("Tracks/Control/Chi2NDF"), chi2ndf);
            registry.fill(HIST("Tracks/2Danalysis/EtaZvtx_sel8"), track.eta(), z);
            if (midtracks.size() > 0 && retrack.ambDegree() > 0) {
              registry.fill(HIST("Tracks/2Danalysis/EtaZvtx_sel8_inelgt0"), track.eta(), z);
            }
          }
        }
        if (retracks.size() > 0) {
          registry.fill(HIST("EventSelection"), 6.);
          if (midtracks.size() > 0) {
            registry.fill(HIST("Tracks/2Danalysis/EventsNtrkZvtx_sel8_inelfwdgt0"), Ntrk, z);
          }
          for (const auto& retrack : retracks) {
            auto track = retrack.mfttrack();
            float ndf = std::max(2.0f * track.nClusters() - 5.0f, 1.0f);
            float chi2ndf = track.chi2() / ndf;
            float phi = track.phi();
            float dcaxy_cut = retrack.bestDCAXY();
            o2::math_utils::bringTo02Pi(phi);
            // Declare dcaz_cut only if needed below.
            if ((cfgnEta1 < track.eta()) && (track.eta() < cfgnEta2) && track.nClusters() >= cfgnCluster && chi2ndf < cfgChi2NDFMax && (phi > cfgPhiCut1 && phi < cfgPhiCut2)) {
              if (usePhiCut) {
                if ((phi <= 0.02) || ((phi >= 3.10) && (phi <= 3.23)) || (phi >= 6.21))
                  continue;
              }
              if (useDCAxyCut) {
                if (dcaxy_cut > maxDCAxy)
                  continue;
              }
              if constexpr (std::is_same_v<RetracksT, soa::SmallGroups<aod::BestCollisionsFwd3d>>) {
                float dcaz_cut = retrack.bestDCAZ();
                if (useDCAzCut) {
                  if (dcaz_cut > maxDCAz)
                    continue;
                }
              }
              // Purity denominator: all tracks that pass the DCA selection and other quality cuts
              registry.fill(HIST("Purity/SelectedAfterDCAxy/All"), 1.);
              registry.fill(HIST("Purity/SelectedAfterDCAxy/AllEta"), track.eta());
              registry.fill(HIST("TracksEtaZvtx"), track.eta(), z);
              if (midtracks.size() > 0 && retrack.ambDegree() > 0) {
                registry.fill(HIST("Tracks/EtaZvtx_gt0"), track.eta(), z);
                registry.fill(HIST("Tracks/2Danalysis/EtaZvtx_sel8_inelfwdgt0"), track.eta(), z);
                eventsInelMFT.insert(retrack.bestCollisionId());
              }
              if (retrack.ambDegree() != 0) {
                registry.fill(HIST("Tracks/Control/woOrp/woOrpEtaZvtx_gt0"), track.eta(), z);
                ++k;
              }
              float phi = track.phi();
              o2::math_utils::bringTo02Pi(phi);
              registry.fill(HIST("Tracks/Control/TrackCount"), 0);
              registry.fill(HIST("TracksPhiEta"), phi, track.eta());
              registry.fill(HIST("TracksPtEta"), track.pt(), track.eta());
              if ((track.eta() < -2.0f) && (track.eta() > -3.9f)) {
                registry.fill(HIST("TracksPhiZvtx"), phi, z);
              }
              if (track.collisionId() > -1 && retrack.ambDegree() == 1) {
                registry.fill(HIST("Tracks/Control/TrackCount"), 8);
                registry.fill(HIST("collisionID"), track.collisionId());
              }
              if (track.collisionId() > -1 && retrack.ambDegree() > 1) {
                registry.fill(HIST("collisionIDamb"), track.collisionId());
              }
              if (track.collisionId() != retrack.bestCollisionId()) {
                registry.fill(HIST("Tracks/Control/ReassignedTracksEtaZvtx"),
                              track.eta(), z);
                registry.fill(HIST("Tracks/Control/ReassignedTracksPhiEta"), phi,
                              track.eta());
                registry.fill(HIST("Tracks/Control/ReassignedVertexCorr"),
                              track.template collision_as<CollwEv>().posZ(), z);

                registry.fill(HIST("Tracks/Control/DeltaZ"),
                              track.template collision_as<CollwEv>().posZ() -
                                collision.posZ());
                registry.fill(HIST("Tracks/Control/TrackCount"), 1);
              }
              if (track.collisionId() == retrack.bestCollisionId()) {
                registry.fill(HIST("Tracks/Control/notReassignedTracksEtaZvtx"),
                              track.eta(), z);
                registry.fill(HIST("Tracks/Control/notReassignedTracksPhiEta"), phi,
                              track.eta());
                registry.fill(HIST("Tracks/Control/notReassignedVertexCorr"),
                              track.template collision_as<CollwEv>().posZ(), z);
                registry.fill(HIST("Tracks/Control/TrackCount"), 2);
              }

              registry.fill(HIST("Tracks/Control/TrackAmbDegree"),
                            retrack.ambDegree());
              registry.fill(HIST("Tracks/Control/DCAXY"), retrack.bestDCAXY());
              if constexpr (std::is_same_v<RetracksT, soa::SmallGroups<aod::BestCollisionsFwd3d>>) {
                registry.fill(HIST("Tracks/Control/DCAZ"), retrack.bestDCAZ());
              }
              int isAmbiguous = 0;

              if (retrack.ambDegree() > 1 && retrack.ambDegree() != 0) {
                isAmbiguous = 1;
                ++i;

                registry.fill(HIST("Tracks/Control/amb/EtaZvtxAmb_gt0"), track.eta(), z);

                registry.fill(HIST("Tracks/Control/amb/AmbTracksEtaZvtx"),
                              track.eta(), z);
                registry.fill(HIST("Tracks/Control/amb/AmbTracksPhiEta"), phi,
                              track.eta());
                registry.fill(HIST("Tracks/Control/amb/AmbVertexCorr"),
                              track.template collision_as<CollwEv>().posZ(), z);
                registry.fill(HIST("Tracks/Control/amb/DCAxy_amb"), retrack.bestDCAXY());
                if constexpr (std::is_same_v<RetracksT, soa::SmallGroups<aod::BestCollisionsFwd3d>>) {
                  registry.fill(HIST("Tracks/Control/amb/DCAz_amb"), retrack.bestDCAZ());
                }
                registry.fill(HIST("Tracks/Control/TrackCount"), 3);
                if (track.collisionId() == retrack.bestCollisionId()) {
                  registry.fill(HIST("Tracks/Control/TrackCount"), 5);
                }
                uniqueEventsAmb.insert(retrack.bestCollisionId());
              }
              if (midtracks.size() > 0 && retrack.ambDegree() > 1 && retrack.ambDegree() != 0) {
                uniqueCollisionsAmb.insert(collision.globalIndex());
              }

              registry.fill(HIST("Tracks/Control/TrackIsAmb"), isAmbiguous);
              if (retrack.ambDegree() == 1 && retrack.ambDegree() != 0) {
                ++j;
                registry.fill(HIST("Tracks/Control/nonamb/EtaZvtxNonAmb_gt0"), track.eta(), z);
                registry.fill(HIST("Tracks/Control/nonamb/nonAmbTracksEtaZvtx"),
                              track.eta(), z);
                registry.fill(HIST("Tracks/Control/nonamb/nonAmbTracksPhiEta"), phi,
                              track.eta());
                registry.fill(HIST("Tracks/Control/nonamb/nonAmbVertexCorr"),
                              track.template collision_as<CollwEv>().posZ(), z);
                registry.fill(HIST("Tracks/Control/nonamb/DCAxy_nonamb"), retrack.bestDCAXY());
                if constexpr (std::is_same_v<RetracksT, soa::SmallGroups<aod::BestCollisionsFwd3d>>) {
                  registry.fill(HIST("Tracks/Control/nonamb/DCAz_nonamb"), retrack.bestDCAZ());
                }
                registry.fill(HIST("Tracks/Control/TrackCount"), 4);
                if (track.collisionId() == retrack.bestCollisionId()) {
                  registry.fill(HIST("Tracks/Control/TrackCount"), 6);
                }
                uniqueEvents.insert(retrack.bestCollisionId());
              }
              if (midtracks.size() > 0 && retrack.ambDegree() == 1 && retrack.ambDegree() != 0) {
                uniqueCollisions.insert(collision.globalIndex());
              }
              if ((retrack.ambDegree() > 1) || (retrack.ambDegree() <= 1))
                registry.fill(HIST("Tracks/Control/TrackCount"), 7);
              if (retrack.ambDegree() != 0) {
                registry.fill(HIST("Tracks/Control/woOrp/woOrpTracksEtaZvtx"),
                              track.eta(), z);
                registry.fill(HIST("Tracks/Control/woOrp/woOrpTracksPhiEta"), phi,
                              track.eta());
                registry.fill(HIST("Tracks/Control/woOrp/woOrpVertexCorr"),
                              track.template collision_as<CollwEv>().posZ(), z);
                registry.fill(HIST("Tracks/Control/TrackCount"), 9); // without orphan
                registry.fill(HIST("Tracks/Control/woOrp/DCAxy_woOrp"), retrack.bestDCAXY());
                if constexpr (std::is_same_v<RetracksT, soa::SmallGroups<aod::BestCollisionsFwd3d>>) {
                  registry.fill(HIST("Tracks/Control/woOrp/DCAz_woOrp"), retrack.bestDCAZ());
                }
              }
            }
          }
          registry.fill(HIST("ambEventCounts"), 1, uniqueEventsAmb.size());
          registry.fill(HIST("NonambEventCounts"), 1, uniqueEvents.size());
          registry.fill(HIST("hNumCollisionsNonAmb_InelMFT"), 1, uniqueCollisions.size());
          registry.fill(HIST("hNumCollisionsAmb_InelMFT"), 1, uniqueCollisionsAmb.size());
          registry.fill(HIST("hNumCollisions_InelMFT"), 1, eventsInelMFT.size());
        }
        registry.fill(HIST("Tracks/Control/amb/nTrkAmb"), i);
        registry.fill(HIST("Tracks/Control/nonamb/nTrkNonAmb"), j);
        registry.fill(HIST("Tracks/Control/woOrp/nTrk"), k);
        registry.fill(HIST("hNumCollisions_Inel"), 1, eventsInel.size());
      }
    } else {
      registry.fill(HIST("EventSelection"), 7);
    }
  }

  void processMultReassoc(CollwEv::iterator const& collision,
                          o2::aod::MFTTracks const& mft,
                          soa::SmallGroups<aod::BestCollisionsFwd> const& retracks,
                          FiCentralTracks const& midtracks, aod::Tracks const& trk)
  {
    processMultReassocCommon(collision, mft, retracks, midtracks, trk);
  }

  void processMultReassoc3d(CollwEv::iterator const& collision,
                            o2::aod::MFTTracks const& mft,
                            soa::SmallGroups<aod::BestCollisionsFwd3d> const& retracks,
                            FiCentralTracks const& midtracks, aod::Tracks const& trk)
  {
    processMultReassocCommon(collision, mft, retracks, midtracks, trk);
  }
  PROCESS_SWITCH(PseudorapidityDensityMFT, processMultReassoc,
                 "Process reco or data info", false);

  PROCESS_SWITCH(PseudorapidityDensityMFT, processMultReassoc3d,
                 "Process reco or data info (3d)", false);

  using ExColsCent = soa::Join<aod::Collisions, aod::CentFT0Cs, aod::EvSels>;

  void processCountingCentrality(ExColsCent::iterator const& collision,
                                 aod::MFTTracks const& tracks)
  {
    auto c = collision.centFT0C();
    registry.fill(HIST("Events/Centrality/Selection"), 1., c);

    if (!useEvSel || collision.sel8()) {
      auto z = collision.posZ();
      registry.fill(HIST("Events/Centrality/Selection"), 2., c);
      auto perCollisionSample = sample->sliceByCached(
        o2::aod::fwdtrack::collisionId, collision.globalIndex(), cache);
      auto Ntrk = perCollisionSample.size();

      registry.fill(HIST("Events/Centrality/NtrkZvtx"), Ntrk, z, c);

      for (const auto& track : tracks) {

        float phi = track.phi();
        o2::math_utils::bringTo02Pi(phi);

        if (usePhiCut) {
          if ((phi < cfgPhiCut) ||
              ((phi > o2::constants::math::PI - cfgPhiCut) && (phi < o2::constants::math::PI + cfgPhiCut)) ||
              (phi > o2::constants::math::TwoPI - cfgPhiCut) ||
              ((phi > ((o2::constants::math::PIHalf - 0.1) * o2::constants::math::PI) - cfgPhiCut) &&
               (phi < ((o2::constants::math::PIHalf - 0.1) * o2::constants::math::PI) + cfgPhiCut)))
            continue;
        }

        registry.fill(HIST("Tracks/Centrality/EtaZvtx"), track.eta(), z, c);
        registry.fill(HIST("Tracks/Centrality/PhiEta"), phi, track.eta(), c);
      }

    } else {
      registry.fill(HIST("Events/Centrality/Selection"), 3.,
                    c); // rejected events
    }
  }

  PROCESS_SWITCH(PseudorapidityDensityMFT, processCountingCentrality,
                 "Count tracks in centrality bins", false);

  using Particles = soa::Filtered<aod::McParticles>;
  expressions::Filter primaries =
    (aod::mcparticle::flags &
     (uint8_t)o2::aod::mcparticle::enums::PhysicalPrimary) ==
    (uint8_t)o2::aod::mcparticle::enums::PhysicalPrimary;
  Partition<Particles> mcSample = nabs(aod::mcparticle::eta) < 1.1f;
  Partition<Particles> mcSampleCentral =
    nabs(aod::mcparticle::eta) < estimatorEta;

  void processGen(
    aod::McCollisions::iterator const& mcCollision,
    o2::soa::SmallGroups<soa::Join<aod::Collisions, aod::EvSels,
                                   aod::McCollisionLabels>> const& collisions,
    Particles const& particles, aod::MFTTracks const& /*tracks*/,
    FiCentralTracks const& midtracks)
  {
    registry.fill(HIST("EventEfficiency"), 1.);

    auto perCollisionMCSample = mcSample->sliceByCached(
      aod::mcparticle::mcCollisionId, mcCollision.globalIndex(), cache);
    auto nCharged = 0;
    for (const auto& particle : perCollisionMCSample) {
      auto charge = 0.;
      auto p = pdg->GetParticle(particle.pdgCode());
      if (p != nullptr) {
        charge = p->Charge();
      }
      if (std::abs(charge) < 3.) {
        continue;
      }
      nCharged++;
    }
    registry.fill(HIST("EventsNtrkZvtxGen_t"), nCharged, mcCollision.posZ());

    //--------for INEL>0
    auto perCollisionMCSampleCentral = mcSampleCentral->sliceByCached(
      aod::mcparticle::mcCollisionId, mcCollision.globalIndex(), cache);
    auto nChargedCentral = 0;
    for (const auto& particle : perCollisionMCSample) {
      auto charge = 0.;
      auto p = pdg->GetParticle(particle.pdgCode());
      if (p != nullptr) {
        charge = p->Charge();
      }
      if (std::abs(charge) < 3.) {
        continue;
      }
      nChargedCentral++;
    }
    if ((mcCollision.posZ() >= cfgVzCut1) && (mcCollision.posZ() <= cfgVzCut2)) {
      if (nChargedCentral > 0) {
        registry.fill(HIST("EventEfficiency"), 2.);
        registry.fill(HIST("EventsNtrkZvtxGen_gt0t"), nCharged,
                      mcCollision.posZ());
      }
    }
    //-----------
    bool atLeastOne = false;
    bool atLeastOne_gt0 = false;
    int moreThanOne = 0;

    LOGP(debug, "MC col {} has {} reco cols", mcCollision.globalIndex(),
         collisions.size());
    for (const auto& collision : collisions) {
      registry.fill(HIST("EventEfficiency"), 3.);
      if (!disableITSROFCut && !collision.selection_bit(aod::evsel::kNoITSROFrameBorder)) {
        return;
      }
      if (!useEvSel || (useEvSel && collision.selection_bit(aod::evsel::kIsTriggerTVX) && collision.selection_bit(aod::evsel::kNoTimeFrameBorder) && collision.selection_bit(aod::evsel::kNoSameBunchPileup))) {
        atLeastOne = true;
        auto perCollisionSample = sample->sliceByCached(
          o2::aod::fwdtrack::collisionId, collision.globalIndex(), cache);

        registry.fill(HIST("EventEfficiency"), 4.);

        auto perCollisionSampleCentral =
          midtracks.sliceBy(perColCentral, collision.globalIndex());
        if ((collision.posZ() >= cfgVzCut1) && (collision.posZ() <= cfgVzCut2) && (mcCollision.posZ() >= cfgVzCut1) && (mcCollision.posZ() <= cfgVzCut2)) {
          if (perCollisionSampleCentral.size() > 0) {
            registry.fill(HIST("EventEfficiency"), 5.);
            atLeastOne_gt0 = true;
            registry.fill(HIST("EventsNtrkZvtxGen_gt0"),
                          perCollisionSample.size(), collision.posZ());
          }

          registry.fill(HIST("EventsZposDiff"),
                        collision.posZ() - mcCollision.posZ());
          if (useZDiffCut) {
            if (std::abs(collision.posZ() - mcCollision.posZ()) > maxZDiff) {
              continue;
            }
          }
          registry.fill(HIST("EventsNtrkZvtxGen"), perCollisionSample.size(),
                        collision.posZ());
          ++moreThanOne;
        }
      }
    }
    if (collisions.size() == 0) {
      registry.fill(HIST("NotFoundEventZvtx"), mcCollision.posZ());
    }
    if (moreThanOne > 1) {
      registry.fill(HIST("EventsSplitMult"), nCharged);
    }
    if ((mcCollision.posZ() >= cfgVzCut1) && (mcCollision.posZ() <= cfgVzCut2)) {
      for (const auto& particle : particles) {
        auto p = pdg->GetParticle(particle.pdgCode());
        auto charge = 0;
        if (p != nullptr) {
          charge = static_cast<int>(p->Charge());
        }
        if (std::abs(charge) < 3.) {
          continue;
        }
        float phi = particle.phi();
        o2::math_utils::bringTo02Pi(phi);
        if (usePhiCut) {
          if ((phi <= 0.02) || ((phi >= 3.10) && (phi <= 3.23)) || (phi >= 6.21))
            continue;
        }
        if (cfgnEta1 < particle.eta() && particle.eta() < cfgnEta2 && (phi > cfgPhiCut1 && phi < cfgPhiCut2)) {
          // Purity numerator reference at generator level: physical primaries in the same eta window
          if (particle.isPhysicalPrimary()) {
            registry.fill(HIST("Purity/Gen/PrimaryEta"), particle.eta());
            // Truth-side total counters for primaries in acceptance (for purity calculations)
            registry.fill(HIST("Purity/Gen/All"), 1.);
            registry.fill(HIST("Purity/Gen/AllEta"), particle.eta());
          }
          registry.fill(HIST("TracksEtaZvtxGen_t"), particle.eta(),
                        mcCollision.posZ());
          if (perCollisionMCSampleCentral.size() > 0) {
            registry.fill(HIST("TracksEtaZvtxGen_gt0t"), particle.eta(),
                          mcCollision.posZ());
            registry.fill(HIST("TracksPhiEtaGen_gt0t"), particle.phi(), particle.eta());
          }
          if (atLeastOne) {
            registry.fill(HIST("TracksEtaZvtxGen"), particle.eta(),
                          mcCollision.posZ());
            registry.fill(HIST("TracksPtEtaGen"), particle.pt(), particle.eta());
            if (atLeastOne_gt0) {
              registry.fill(HIST("TracksEtaZvtxGen_gt0"), particle.eta(),
                            mcCollision.posZ());
              registry.fill(HIST("TracksPhiEtaGen_gt0"), particle.phi(), particle.eta());
            }
          }

          registry.fill(HIST("TracksPhiEtaGen"), particle.phi(), particle.eta());
          registry.fill(HIST("TracksPhiZvtxGen"), particle.phi(),
                        mcCollision.posZ());
          registry.fill(HIST("TracksPtEtaGen_t"), particle.pt(), particle.eta());
        }
      }
    }
  }

  PROCESS_SWITCH(PseudorapidityDensityMFT, processGen,
                 "Process generator-level info", false);

  using ExColsGenCent =
    soa::SmallGroups<soa::Join<aod::McCollisionLabels, aod::Collisions,
                               aod::CentFT0Cs, aod::EvSels>>;

  void processGenCent(aod::McCollisions::iterator const& mcCollision,
                      ExColsGenCent const& collisions,
                      Particles const& particles,
                      MFTTracksLabeled const& /*tracks*/)
  {

    LOGP(debug, "MC col {} has {} reco cols", mcCollision.globalIndex(),
         collisions.size());

    float c_gen = -1;
    bool atLeastOne = false;
    for (const auto& collision : collisions) {
      float c_rec = -1;
      if constexpr (ExColsGenCent::template contains<aod::CentFT0Cs>()) {
        c_rec = collision.centFT0C();
      }
      if (!useEvSel || (useEvSel && collision.sel8())) {
        if constexpr (ExColsGenCent::template contains<aod::CentFT0Cs>()) {
          if (!atLeastOne) {
            c_gen = c_rec;
          }
        }
        atLeastOne = true;

        registry.fill(HIST("Events/Centrality/EventEfficiency"), 2., c_gen);
        registry.fill(HIST("Events/Centrality/CentPercentileMCGen"), c_gen);

        auto perCollisionSample = sample->sliceByCached(
          o2::aod::fwdtrack::collisionId, collision.globalIndex(), cache);
        registry.fill(HIST("Events/Centrality/NtrkZvtxGen"),
                      perCollisionSample.size(), collision.posZ(), c_gen);
      }
    }

    registry.fill(HIST("Events/Centrality/EventEfficiency"), 1., c_gen);

    auto perCollisionMCSample = mcSample->sliceByCached(
      aod::mcparticle::mcCollisionId, mcCollision.globalIndex(), cache);
    auto nCharged = 0;

    for (const auto& particle : perCollisionMCSample) {
      auto p = pdg->GetParticle(particle.pdgCode());
      auto charge = 0;
      if (p != nullptr) {
        charge = static_cast<int>(p->Charge());
      }
      if (std::abs(charge) < 3.) {
        continue;
      }
      nCharged++;
    }

    if constexpr (ExColsGenCent::template contains<aod::CentFT0Cs>()) {
      registry.fill(HIST("Events/Centrality/NtrkZvtxGen_t"), nCharged,
                    mcCollision.posZ(), c_gen);
    }

    for (const auto& particle : particles) {
      auto p = pdg->GetParticle(particle.pdgCode());
      auto charge = 0;
      if (p != nullptr) {
        charge = static_cast<int>(p->Charge());
      }
      if (std::abs(charge) < 3.) {
        continue;
      }

      if constexpr (ExColsGenCent::template contains<aod::CentFT0Cs>()) {
        registry.fill(HIST("Tracks/Centrality/EtaZvtxGen_t"), particle.eta(),
                      mcCollision.posZ(), c_gen);
      }

      if (atLeastOne) {
        if constexpr (ExColsGenCent::template contains<aod::CentFT0Cs>()) {
          registry.fill(HIST("Tracks/Centrality/EtaZvtxGen"), particle.eta(),
                        mcCollision.posZ(), c_gen);
          float phi = particle.phi();
          o2::math_utils::bringTo02Pi(phi);
          registry.fill(HIST("Tracks/Centrality/PhiEtaGen"), phi,
                        particle.eta(), c_gen);
        }
      }
    }
  }

  PROCESS_SWITCH(PseudorapidityDensityMFT, processGenCent,
                 "Process generator-level info in centrality bins", false);

  void processGenPt(
    soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision,
    MFTTracksLabeled const& tracks, aod::McParticles const&)
  {
    if (!useEvSel || (useEvSel && collision.sel8())) {
      for (const auto& track : tracks) {
        if (!track.has_mcParticle()) {
          continue;
        }
        auto particle = track.mcParticle();
        if (!particle.isPhysicalPrimary()) {
          continue;
        }
        registry.fill(HIST("TracksToPartPtEta"), particle.pt(), particle.eta());
      }
    }
  }

  PROCESS_SWITCH(PseudorapidityDensityMFT, processGenPt,
                 "Process particle-level info of pt", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<PseudorapidityDensityMFT>(cfgc)};
}
