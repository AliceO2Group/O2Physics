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
// 2-particle correlations for cascades
// =============================
//
// Author: Rik Spijkers (rik.spijkers@cern.ch)
//

#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/Utils/inelGt.h"

#include "Common/Core/Zorro.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

#include <TFile.h>
#include <TH2F.h>
#include <TList.h>
#include <TPDGCode.h>
#include <TProfile.h>

#include <array>
#include <cmath>
#include <cstdlib>
#include <string>
#include <utility>

using namespace o2;
using namespace o2::soa;
using namespace o2::constants::math;
using namespace o2::constants::physics;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;

// use parameters + cov mat non-propagated, aux info + (extension propagated)
using FullTracksExt = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCov>;
using FullTracksExtIU = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU>;
using FullTracksExtWithPID = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCov, aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr>;
using FullTracksExtIUWithPID = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr>;

Zorro zorro;

// Add a column to the cascdataext table: IsSelected.
// 0 = not selected, 1 = Xi, 2 = both, 3 = Omega
namespace o2::aod
{
namespace cascadeflags
{
DECLARE_SOA_COLUMN(IsSelected, isSelected, int); //~!
} // namespace cascadeflags
DECLARE_SOA_TABLE(CascadeFlags, "AOD", "CASCADEFLAGS", //!
                  cascadeflags::IsSelected);
using CascDataExtSelected = soa::Join<CascDataExt, CascadeFlags>;
} // namespace o2::aod

using MyCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::PVMults>;
using MyCollisionsMult = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms>;
using MyCascades = soa::Filtered<aod::CascDataExtSelected>;
using LabeledCascades = soa::Join<aod::CascDataExt, aod::McCascLabels>;

struct CascadeSelector {
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Service<o2::framework::O2DatabasePDG> pdgDB;

  Produces<aod::CascadeFlags> cascflags;

  // Configurables
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "CCDB url"};
  Configurable<bool> useTrigger{"useTrigger", false, "Use trigger selection on skimmed data"};
  Configurable<std::string> triggerList{"triggerList", "fDoubleXi, fDoubleOmega, fOmegaXi", "List of triggers used to select events"};
  Configurable<bool> doTFBorderCut{"doTFBorderCut", true, "Switch to apply TimeframeBorderCut event selection"};
  Configurable<bool> doSel8{"doSel8", true, "Switch to apply sel8 event selection"};
  Configurable<bool> doNoSameBunchPileUp{"doNoSameBunchPileUp", true, "Switch to apply NoSameBunchPileUp event selection"};
  Configurable<int> INEL{"INEL", 0, "Number of charged tracks within |eta| < 1 has to be greater than value"};
  Configurable<double> maxVertexZ{"maxVertexZ", 10., "Maximum value of z coordinate of PV"};
  Configurable<float> etaCascades{"etaCascades", 0.8, "min/max of eta for cascades"};
  Configurable<bool> doCompetingMassCut{"doCompetingMassCut", true, "Switch to apply a competing mass cut for the Omega's"};
  Configurable<float> competingMassWindow{"competingMassWindow", 0.01, "Mass window for the competing mass cut"};

  // Tracklevel
  Configurable<float> tpcNsigmaBachelor{"tpcNsigmaBachelor", 3, "TPC NSigma bachelor"};
  Configurable<float> tpcNsigmaProton{"tpcNsigmaProton", 3, "TPC NSigma proton <- lambda"};
  Configurable<float> tpcNsigmaPion{"tpcNsigmaPion", 3, "TPC NSigma pion <- lambda"};
  Configurable<int> minTPCCrossedRows{"minTPCCrossedRows", 80, "min N TPC crossed rows"}; // TODO: finetune! 80 > 159/2, so no split tracks?
  Configurable<int> minITSClusters{"minITSClusters", 4, "minimum number of ITS clusters"};
  Configurable<float> etaTracks{"etaTracks", 1.0, "min/max of eta for tracks"};
  Configurable<float> tpcChi2{"tpcChi2", 4, "TPC Chi2"};
  Configurable<float> itsChi2{"itsChi2", 36, "ITS Chi2"};

  // Selection criteria - compatible with core wagon autodetect - copied from cascadeanalysis.cxx
  //*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
  Configurable<double> v0setting_cospa{"v0setting_cospa", 0.995, "v0setting_cospa"};
  Configurable<float> v0setting_dcav0dau{"v0setting_dcav0dau", 1.0, "v0setting_dcav0dau"};
  Configurable<float> v0setting_dcapostopv{"v0setting_dcapostopv", 0.1, "v0setting_dcapostopv"};
  Configurable<float> v0setting_dcanegtopv{"v0setting_dcanegtopv", 0.1, "v0setting_dcanegtopv"};
  Configurable<float> v0setting_radius{"v0setting_radius", 0.9, "v0setting_radius"};
  Configurable<double> cascadesetting_cospa{"cascadesetting_cospa", 0.95, "cascadesetting_cospa"};
  Configurable<float> cascadesetting_dcacascdau{"cascadesetting_dcacascdau", 1.0, "cascadesetting_dcacascdau"};
  Configurable<float> cascadesetting_dcabachtopv{"cascadesetting_dcabachtopv", 0.05, "cascadesetting_dcabachtopv"};
  Configurable<float> cascadesetting_cascradius{"cascadesetting_cascradius", 0.9, "cascadesetting_cascradius"};
  Configurable<float> cascadesetting_v0masswindow{"cascadesetting_v0masswindow", 0.01, "cascadesetting_v0masswindow"};
  Configurable<float> cascadesetting_mindcav0topv{"cascadesetting_mindcav0topv", 0.01, "cascadesetting_mindcav0topv"};
  //*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*

  // TODO: variables as function of Omega mass, only do Xi for now
  ConfigurableAxis radiusAxis = {"radiusAxis", {100, 0.0f, 50.0f}, "cm"};
  ConfigurableAxis cpaAxis = {"cpaAxis", {100, 0.95f, 1.0f}, "CPA"};
  ConfigurableAxis vertexAxis = {"vertexAxis", {100, -10.0f, 10.0f}, "cm"};
  ConfigurableAxis dcaAxis = {"dcaAxis", {100, 0.0f, 2.0f}, "cm"};
  ConfigurableAxis invXiMassAxis = {"invXiMassAxis", {100, 1.28f, 1.38f}, "Inv. Mass (GeV/c^{2})"};
  ConfigurableAxis invOmegaMassAxis = {"invOmegaMassAxis", {100, 1.62f, 1.72f}, "Inv. Mass (GeV/c^{2})"};
  ConfigurableAxis ptAxis = {"ptAxis", {150, 0, 15}, "#it{p}_{T}"};
  ConfigurableAxis etaAxis{"etaAxis", {100, -1.f, 1.f}, "#eta"};
  ConfigurableAxis invLambdaMassAxis{"invLambdaMassAxis", {100, 1.07f, 1.17f}, "Inv. Mass (GeV/c^{2})"};
  ConfigurableAxis chi2Axis{"chi2Axis", {100, 0.f, 10.f}, "Chi2"};
  AxisSpec itsClustersAxis{8, -0.5, 7.5, "number of ITS clusters"};
  AxisSpec tpcRowsAxis{160, -0.5, 159.5, "TPC crossed rows"};
  HistogramRegistry registry{
    "registry",
    {
      // basic selection variables
      {"hV0Radius", "hV0Radius", {HistType::kTH3F, {radiusAxis, invXiMassAxis, ptAxis}}},
      {"hCascRadius", "hCascRadius", {HistType::kTH3F, {radiusAxis, invXiMassAxis, ptAxis}}},
      {"hV0CosPA", "hV0CosPA", {HistType::kTH3F, {cpaAxis, invXiMassAxis, ptAxis}}},
      {"hCascCosPA", "hCascCosPA", {HistType::kTH3F, {cpaAxis, invXiMassAxis, ptAxis}}},
      {"hDCAPosToPV", "hDCAPosToPV", {HistType::kTH3F, {dcaAxis, invXiMassAxis, ptAxis}}},
      {"hDCANegToPV", "hDCANegToPV", {HistType::kTH3F, {dcaAxis, invXiMassAxis, ptAxis}}},
      {"hDCABachToPV", "hDCABachToPV", {HistType::kTH3F, {dcaAxis, invXiMassAxis, ptAxis}}},
      {"hDCAV0ToPV", "hDCAV0ToPV", {HistType::kTH3F, {dcaAxis, invXiMassAxis, ptAxis}}},
      {"hDCAV0Dau", "hDCAV0Dau", {HistType::kTH3F, {dcaAxis, invXiMassAxis, ptAxis}}},
      {"hDCACascDau", "hDCACascDau", {HistType::kTH3F, {dcaAxis, invXiMassAxis, ptAxis}}},
      {"hLambdaMass", "hLambdaMass", {HistType::kTH3F, {invLambdaMassAxis, invXiMassAxis, ptAxis}}},

      {"hMassXiMinus", "hMassXiMinus", {HistType::kTH3F, {invXiMassAxis, ptAxis, etaAxis}}},
      {"hMassXiPlus", "hMassXiPlus", {HistType::kTH3F, {invXiMassAxis, ptAxis, etaAxis}}},
      {"hMassOmegaMinus", "hMassOmegaMinus", {HistType::kTH3F, {invOmegaMassAxis, ptAxis, etaAxis}}},
      {"hMassOmegaPlus", "hMassOmegaPlus", {HistType::kTH3F, {invOmegaMassAxis, ptAxis, etaAxis}}},

      // ITS & TPC clusters, with Xi inv mass
      {"hTPCnCrossedRowsPos", "hTPCnCrossedRowsPos", {HistType::kTH3F, {tpcRowsAxis, invXiMassAxis, ptAxis}}},
      {"hTPCnCrossedRowsNeg", "hTPCnCrossedRowsNeg", {HistType::kTH3F, {tpcRowsAxis, invXiMassAxis, ptAxis}}},
      {"hTPCnCrossedRowsBach", "hTPCnCrossedRowsBach", {HistType::kTH3F, {tpcRowsAxis, invXiMassAxis, ptAxis}}},
      {"hITSnClustersPos", "hITSnClustersPos", {HistType::kTH3F, {itsClustersAxis, invXiMassAxis, ptAxis}}},
      {"hITSnClustersNeg", "hITSnClustersNeg", {HistType::kTH3F, {itsClustersAxis, invXiMassAxis, ptAxis}}},
      {"hITSnClustersBach", "hITSnClustersBach", {HistType::kTH3F, {itsClustersAxis, invXiMassAxis, ptAxis}}},
      {"hTPCChi2Pos", "hTPCChi2Pos", {HistType::kTH1F, {chi2Axis}}},
      {"hTPCChi2Neg", "hTPCChi2Neg", {HistType::kTH1F, {chi2Axis}}},
      {"hTPCChi2Bach", "hTPCChi2Bach", {HistType::kTH1F, {chi2Axis}}},
      {"hITSChi2Pos", "hITSChi2Pos", {HistType::kTH1F, {chi2Axis}}},
      {"hITSChi2Neg", "hITSChi2Neg", {HistType::kTH1F, {chi2Axis}}},
      {"hITSChi2Bach", "hITSChi2Bach", {HistType::kTH1F, {chi2Axis}}},

      {"hTriggerQA", "hTriggerQA", {HistType::kTH1F, {{2, -0.5, 1.5, "Trigger y/n"}}}},
    },
  };

  // Keep track of which selections the candidates pass
  void init(InitContext const&)
  {
    ccdb->setURL(ccdbUrl);
    ccdb->setCaching(true);

    auto h = registry.add<TH1>("hSelectionStatus", "hSelectionStatus", HistType::kTH1F, {{10, 0, 10, "status"}});
    h->GetXaxis()->SetBinLabel(1, "All");
    h->GetXaxis()->SetBinLabel(2, "nTPC OK");
    h->GetXaxis()->SetBinLabel(3, "nITS OK");
    h->GetXaxis()->SetBinLabel(4, "track Chi2 OK");
    h->GetXaxis()->SetBinLabel(5, "Topo OK");
    h->GetXaxis()->SetBinLabel(6, "Track eta OK");
    h->GetXaxis()->SetBinLabel(7, "Cascade eta OK");
    h->GetXaxis()->SetBinLabel(8, "V0 PID OK");
    h->GetXaxis()->SetBinLabel(9, "Bach PID OK");

    auto hEventSel = registry.add<TH1>("hEventSel", "hEventSel", HistType::kTH1F, {{10, 0, 10, "selection criteria"}});
    hEventSel->GetXaxis()->SetBinLabel(1, "All");
    hEventSel->GetXaxis()->SetBinLabel(2, "sel8");
    hEventSel->GetXaxis()->SetBinLabel(3, "INEL0");
    hEventSel->GetXaxis()->SetBinLabel(4, "V_z");
    hEventSel->GetXaxis()->SetBinLabel(5, "NoSameBunchPileUp");
    hEventSel->GetXaxis()->SetBinLabel(6, "Selected events");

    if (doprocessRecMC) {
      // only create the rec matched to gen histograms if relevant
      registry.add("truerec/hV0Radius", "hV0Radius", HistType::kTH1F, {radiusAxis});
      registry.add("truerec/hCascRadius", "hCascRadius", HistType::kTH1F, {radiusAxis});
      registry.add("truerec/hV0CosPA", "hV0CosPA", HistType::kTH1F, {cpaAxis});
      registry.add("truerec/hCascCosPA", "hCascCosPA", HistType::kTH1F, {cpaAxis});
      registry.add("truerec/hDCAPosToPV", "hDCAPosToPV", HistType::kTH1F, {dcaAxis});
      registry.add("truerec/hDCANegToPV", "hDCANegToPV", HistType::kTH1F, {dcaAxis});
      registry.add("truerec/hDCABachToPV", "hDCABachToPV", HistType::kTH1F, {dcaAxis});
      registry.add("truerec/hDCAV0ToPV", "hDCAV0ToPV", HistType::kTH1F, {dcaAxis});
      registry.add("truerec/hDCAV0Dau", "hDCAV0Dau", HistType::kTH1F, {dcaAxis});
      registry.add("truerec/hDCACascDau", "hDCACascDau", HistType::kTH1F, {dcaAxis});
      registry.add("truerec/hLambdaMass", "hLambdaMass", HistType::kTH1F, {invLambdaMassAxis});
      registry.add("truerec/hTPCnCrossedRowsPos", "hTPCnCrossedRowsPos", HistType::kTH1F, {tpcRowsAxis});
      registry.add("truerec/hTPCnCrossedRowsNeg", "hTPCnCrossedRowsNeg", HistType::kTH1F, {tpcRowsAxis});
      registry.add("truerec/hTPCnCrossedRowsBach", "hTPCnCrossedRowsBach", HistType::kTH1F, {tpcRowsAxis});
      registry.add("truerec/hITSnClustersPos", "hITSnClustersPos", HistType::kTH1F, {itsClustersAxis});
      registry.add("truerec/hITSnClustersNeg", "hITSnClustersNeg", HistType::kTH1F, {itsClustersAxis});
      registry.add("truerec/hITSnClustersBach", "hITSnClustersBach", HistType::kTH1F, {itsClustersAxis});
      registry.add("truerec/hTPCChi2Pos", "hTPCChi2Pos", HistType::kTH1F, {chi2Axis});
      registry.add("truerec/hTPCChi2Neg", "hTPCChi2Neg", HistType::kTH1F, {chi2Axis});
      registry.add("truerec/hTPCChi2Bach", "hTPCChi2Bach", HistType::kTH1F, {chi2Axis});
      registry.add("truerec/hITSChi2Pos", "hITSChi2Pos", HistType::kTH1F, {chi2Axis});
      registry.add("truerec/hITSChi2Neg", "hITSChi2Neg", HistType::kTH1F, {chi2Axis});
      registry.add("truerec/hITSChi2Bach", "hITSChi2Bach", HistType::kTH1F, {chi2Axis});
      registry.add("truerec/hXiMinus", "hXiMinus", HistType::kTH2F, {ptAxis, etaAxis});
      registry.add("truerec/hXiPlus", "hXiPlus", HistType::kTH2F, {ptAxis, etaAxis});
      registry.add("truerec/hOmegaMinus", "hOmegaMinus", HistType::kTH2F, {ptAxis, etaAxis});
      registry.add("truerec/hOmegaPlus", "hOmegaPlus", HistType::kTH2F, {ptAxis, etaAxis});
    }

    if (doprocessGenMC) {
      // only create the MC gen histograms if relevant
      registry.add("gen/hXiMinus", "hXiMinus", HistType::kTH2F, {ptAxis, etaAxis});
      registry.add("gen/hXiPlus", "hXiPlus", HistType::kTH2F, {ptAxis, etaAxis});
      registry.add("gen/hOmegaMinus", "hOmegaMinus", HistType::kTH2F, {ptAxis, etaAxis});
      registry.add("gen/hOmegaPlus", "hOmegaPlus", HistType::kTH2F, {ptAxis, etaAxis});

      registry.add("genwithrec/hXiMinus", "hXiMinus", HistType::kTH2F, {ptAxis, etaAxis});
      registry.add("genwithrec/hXiPlus", "hXiPlus", HistType::kTH2F, {ptAxis, etaAxis});
      registry.add("genwithrec/hOmegaMinus", "hOmegaMinus", HistType::kTH2F, {ptAxis, etaAxis});
      registry.add("genwithrec/hOmegaPlus", "hOmegaPlus", HistType::kTH2F, {ptAxis, etaAxis});

      registry.add("genwithrec/hNevents", "hNevents", HistType::kTH1F, {{1, 0, 1, "N generated events with reconstructed event"}});
      registry.add("gen/hNevents", "hNevents", HistType::kTH1F, {{1, 0, 1, "N generated events"}});
    }
  }

  template <typename TCollision>
  bool eventSelection(TCollision const& collision, bool fillHistos)
  {
    if (useTrigger) {
      auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
      zorro.initCCDB(ccdb.service, bc.runNumber(), bc.timestamp(), triggerList);
      bool eventTrigger = zorro.isSelected(bc.globalBC());
      if (eventTrigger) {
        if (fillHistos)
          registry.fill(HIST("hTriggerQA"), 1);
      } else {
        if (fillHistos)
          registry.fill(HIST("hTriggerQA"), 0);
        return false;
      }
    }
    // fill event selection based on which selection criteria are applied and passed
    if (fillHistos)
      registry.fill(HIST("hEventSel"), 0);
    if (doSel8 && !collision.sel8()) {
      if (fillHistos)
        registry.fill(HIST("hEventSel"), 1);
      return false;
    } else if (collision.multNTracksPVeta1() <= INEL) {
      if (fillHistos)
        registry.fill(HIST("hEventSel"), 2);
      return false;
    } else if (std::abs(collision.posZ()) > maxVertexZ) {
      if (fillHistos)
        registry.fill(HIST("hEventSel"), 3);
      return false;
    } else if (doNoSameBunchPileUp && !collision.selection_bit(aod::evsel::kNoSameBunchPileup)) {
      if (fillHistos)
        registry.fill(HIST("hEventSel"), 4);
      return false;
    }
    // passes all selections
    if (fillHistos)
      registry.fill(HIST("hEventSel"), 5);
    return true;
  }

  template <typename TCollision>
  void fillMatchedHistos(LabeledCascades::iterator rec, int flag, TCollision collision)
  {
    if (flag == 0)
      return;
    if (!rec.has_mcParticle())
      return;
    auto gen = rec.mcParticle();
    if (!gen.isPhysicalPrimary())
      return;
    int genpdg = gen.pdgCode();
    if ((flag < 3 && std::abs(genpdg) == kXiMinus) || (flag > 1 && std::abs(genpdg) == kOmegaMinus)) {
      // if casc is consistent with Xi and has matched gen Xi OR cand is consistent with Omega and has matched gen omega
      // have to do this in case we reco true Xi with only Omega hypothesis (or vice versa) (very unlikely)
      registry.fill(HIST("truerec/hV0Radius"), rec.v0radius());
      registry.fill(HIST("truerec/hCascRadius"), rec.cascradius());
      registry.fill(HIST("truerec/hV0CosPA"), rec.v0cosPA(collision.posX(), collision.posY(), collision.posZ()));
      registry.fill(HIST("truerec/hCascCosPA"), rec.casccosPA(collision.posX(), collision.posY(), collision.posZ()));
      registry.fill(HIST("truerec/hDCAPosToPV"), rec.dcapostopv());
      registry.fill(HIST("truerec/hDCANegToPV"), rec.dcanegtopv());
      registry.fill(HIST("truerec/hDCABachToPV"), rec.dcabachtopv());
      registry.fill(HIST("truerec/hDCAV0ToPV"), rec.dcav0topv(collision.posX(), collision.posY(), collision.posZ()));
      registry.fill(HIST("truerec/hDCAV0Dau"), rec.dcaV0daughters());
      registry.fill(HIST("truerec/hDCACascDau"), rec.dcacascdaughters());
      registry.fill(HIST("truerec/hLambdaMass"), rec.mLambda());
      registry.fill(HIST("truerec/hITSnClustersPos"), rec.posTrack_as<FullTracksExtIUWithPID>().itsNCls());
      registry.fill(HIST("truerec/hITSnClustersNeg"), rec.negTrack_as<FullTracksExtIUWithPID>().itsNCls());
      registry.fill(HIST("truerec/hITSnClustersBach"), rec.bachelor_as<FullTracksExtIUWithPID>().itsNCls());
      registry.fill(HIST("truerec/hTPCnCrossedRowsPos"), rec.posTrack_as<FullTracksExtIUWithPID>().tpcNClsCrossedRows());
      registry.fill(HIST("truerec/hTPCnCrossedRowsNeg"), rec.negTrack_as<FullTracksExtIUWithPID>().tpcNClsCrossedRows());
      registry.fill(HIST("truerec/hTPCnCrossedRowsBach"), rec.bachelor_as<FullTracksExtIUWithPID>().tpcNClsCrossedRows());
      registry.fill(HIST("truerec/hITSChi2Pos"), rec.posTrack_as<FullTracksExtIUWithPID>().itsChi2NCl());
      registry.fill(HIST("truerec/hITSChi2Neg"), rec.negTrack_as<FullTracksExtIUWithPID>().itsChi2NCl());
      registry.fill(HIST("truerec/hITSChi2Bach"), rec.bachelor_as<FullTracksExtIUWithPID>().itsChi2NCl());
      registry.fill(HIST("truerec/hTPCChi2Pos"), rec.posTrack_as<FullTracksExtIUWithPID>().tpcChi2NCl());
      registry.fill(HIST("truerec/hTPCChi2Neg"), rec.negTrack_as<FullTracksExtIUWithPID>().tpcChi2NCl());
      registry.fill(HIST("truerec/hTPCChi2Bach"), rec.bachelor_as<FullTracksExtIUWithPID>().tpcChi2NCl());
      switch (genpdg) { // is matched so we can use genpdg
        case kXiMinus:
          registry.fill(HIST("truerec/hXiMinus"), rec.pt(), rec.eta());
          break;
        case kXiPlusBar:
          registry.fill(HIST("truerec/hXiPlus"), rec.pt(), rec.eta());
          break;
        case kOmegaMinus:
          registry.fill(HIST("truerec/hOmegaMinus"), rec.pt(), rec.eta());
          break;
        case kOmegaPlusBar:
          registry.fill(HIST("truerec/hOmegaPlus"), rec.pt(), rec.eta());
          break;
      }
    }
  }

  template <typename TCascade, typename TCollision>
  int processCandidate(TCascade const& casc, TCollision const& collision)
  {
    // these are the tracks:
    auto bachTrack = casc.template bachelor_as<FullTracksExtIUWithPID>();
    auto posTrack = casc.template posTrack_as<FullTracksExtIUWithPID>();
    auto negTrack = casc.template negTrack_as<FullTracksExtIUWithPID>();

    // topo variables before cuts:
    registry.fill(HIST("hV0Radius"), casc.v0radius(), casc.mXi(), casc.pt());
    registry.fill(HIST("hCascRadius"), casc.cascradius(), casc.mXi(), casc.pt());
    registry.fill(HIST("hV0CosPA"), casc.v0cosPA(collision.posX(), collision.posY(), collision.posZ()), casc.mXi(), casc.pt());
    registry.fill(HIST("hCascCosPA"), casc.casccosPA(collision.posX(), collision.posY(), collision.posZ()), casc.mXi(), casc.pt());
    registry.fill(HIST("hDCAPosToPV"), casc.dcapostopv(), casc.mXi(), casc.pt());
    registry.fill(HIST("hDCANegToPV"), casc.dcanegtopv(), casc.mXi(), casc.pt());
    registry.fill(HIST("hDCABachToPV"), casc.dcabachtopv(), casc.mXi(), casc.pt());
    registry.fill(HIST("hDCAV0ToPV"), casc.dcav0topv(collision.posX(), collision.posY(), collision.posZ()), casc.mXi(), casc.pt());
    registry.fill(HIST("hDCAV0Dau"), casc.dcaV0daughters(), casc.mXi(), casc.pt());
    registry.fill(HIST("hDCACascDau"), casc.dcacascdaughters(), casc.mXi(), casc.pt());
    registry.fill(HIST("hLambdaMass"), casc.mLambda(), casc.mXi(), casc.pt());

    registry.fill(HIST("hITSnClustersPos"), posTrack.itsNCls(), casc.mXi(), casc.pt());
    registry.fill(HIST("hITSnClustersNeg"), negTrack.itsNCls(), casc.mXi(), casc.pt());
    registry.fill(HIST("hITSnClustersBach"), bachTrack.itsNCls(), casc.mXi(), casc.pt());
    registry.fill(HIST("hTPCnCrossedRowsPos"), posTrack.tpcNClsCrossedRows(), casc.mXi(), casc.pt());
    registry.fill(HIST("hTPCnCrossedRowsNeg"), negTrack.tpcNClsCrossedRows(), casc.mXi(), casc.pt());
    registry.fill(HIST("hTPCnCrossedRowsBach"), bachTrack.tpcNClsCrossedRows(), casc.mXi(), casc.pt());
    registry.fill(HIST("hITSChi2Pos"), posTrack.itsChi2NCl());
    registry.fill(HIST("hITSChi2Neg"), negTrack.itsChi2NCl());
    registry.fill(HIST("hITSChi2Bach"), bachTrack.itsChi2NCl());
    registry.fill(HIST("hTPCChi2Pos"), posTrack.tpcChi2NCl());
    registry.fill(HIST("hTPCChi2Neg"), negTrack.tpcChi2NCl());
    registry.fill(HIST("hTPCChi2Bach"), bachTrack.tpcChi2NCl());

    registry.fill(HIST("hSelectionStatus"), 0); // all the cascade before selections
    // registry.fill(HIST("hMassXi0"), casc.mXi(), casc.pt());

    // TPC N crossed rows todo: check if minTPCCrossedRows > 50
    if (posTrack.tpcNClsCrossedRows() < minTPCCrossedRows || negTrack.tpcNClsCrossedRows() < minTPCCrossedRows || bachTrack.tpcNClsCrossedRows() < minTPCCrossedRows)
      return 0;

    registry.fill(HIST("hSelectionStatus"), 1); // passes nTPC crossed rows
    // registry.fill(HIST("hMassXi1"), casc.mXi(), casc.pt());

    // ITS N clusters todo: check if minITSClusters > 0
    if (posTrack.itsNCls() < minITSClusters || negTrack.itsNCls() < minITSClusters || bachTrack.itsNCls() < minITSClusters)
      return 0;

    registry.fill(HIST("hSelectionStatus"), 2); // passes nITS clusters
    // registry.fill(HIST("hMassXi2"), casc.mXi(), casc.pt());

    // Chi2 cuts
    if (posTrack.itsChi2NCl() > itsChi2 || negTrack.itsChi2NCl() > itsChi2 || bachTrack.itsChi2NCl() > itsChi2)
      return 0;
    if (posTrack.tpcChi2NCl() > tpcChi2 || negTrack.tpcChi2NCl() > tpcChi2 || bachTrack.tpcChi2NCl() > tpcChi2)
      return 0;

    registry.fill(HIST("hSelectionStatus"), 3); // passes Chi2 cuts

    //// TOPO CUTS //// TODO: improve!
    double pvx = collision.posX();
    double pvy = collision.posY();
    double pvz = collision.posZ();
    if (casc.v0radius() < v0setting_radius ||
        casc.cascradius() < cascadesetting_cascradius ||
        casc.v0cosPA(pvx, pvy, pvz) < v0setting_cospa ||
        casc.casccosPA(pvx, pvy, pvz) < cascadesetting_cospa ||
        std::abs(casc.dcav0topv(pvx, pvy, pvz)) < cascadesetting_mindcav0topv ||
        std::abs(casc.mLambda() - o2::constants::physics::MassLambda) > cascadesetting_v0masswindow ||
        std::abs(casc.dcapostopv()) < v0setting_dcapostopv ||
        std::abs(casc.dcanegtopv()) < v0setting_dcanegtopv ||
        casc.dcaV0daughters() > v0setting_dcav0dau ||
        std::abs(casc.dcabachtopv()) < cascadesetting_dcabachtopv ||
        casc.dcacascdaughters() > cascadesetting_dcacascdau)
      return 0; // It failed at least one topo selection

    registry.fill(HIST("hSelectionStatus"), 4); // passes topo
    // registry.fill(HIST("hMassXi3"), casc.mXi(), casc.pt());

    if (std::abs(posTrack.eta()) > etaTracks || std::abs(negTrack.eta()) > etaTracks || std::abs(bachTrack.eta()) > etaTracks)
      return 0;

    registry.fill(HIST("hSelectionStatus"), 5); // passes track eta

    if (std::abs(casc.eta()) > etaCascades)
      return 0;

    registry.fill(HIST("hSelectionStatus"), 6); // passes candidate eta

    // TODO: TOF (for pT > 2 GeV per track?)

    //// TPC PID ////
    // Lambda check
    if (casc.sign() < 0) {
      // Proton check:
      if (std::abs(posTrack.tpcNSigmaPr()) > tpcNsigmaProton)
        return 0;
      // Pion check:
      if (std::abs(negTrack.tpcNSigmaPi()) > tpcNsigmaPion)
        return 0;
    } else {
      // Proton check:
      if (std::abs(negTrack.tpcNSigmaPr()) > tpcNsigmaProton)
        return 0;
      // Pion check:
      if (std::abs(posTrack.tpcNSigmaPi()) > tpcNsigmaPion)
        return 0;
    }
    registry.fill(HIST("hSelectionStatus"), 7); // passes V0 daughters PID
    // registry.fill(HIST("hMassXi4"), casc.mXi(), casc.pt());

    // setting selection flag based on bachelor PID (and competing mass cut for omega's)
    int flag = 0;
    if (std::abs(bachTrack.tpcNSigmaPi()) < tpcNsigmaBachelor)
      flag = 1;
    if (std::abs(bachTrack.tpcNSigmaKa()) < tpcNsigmaBachelor && (!doCompetingMassCut || std::abs(o2::constants::physics::MassXiMinus - casc.mXi()) > competingMassWindow))
      flag = 3 - flag; // 3 if only consistent with omega, 2 if consistent with both

    switch (flag) {
      case 1:                                       // only Xi
        registry.fill(HIST("hSelectionStatus"), 8); // passes bach PID
        if (casc.sign() < 0) {
          registry.fill(HIST("hMassXiMinus"), casc.mXi(), casc.pt(), casc.eta());
        } else {
          registry.fill(HIST("hMassXiPlus"), casc.mXi(), casc.pt(), casc.eta());
        }
        break;
      case 2:                                       // Xi or Omega
        registry.fill(HIST("hSelectionStatus"), 8); // passes bach PID
        if (casc.sign() < 0) {
          registry.fill(HIST("hMassXiMinus"), casc.mXi(), casc.pt(), casc.eta());
          registry.fill(HIST("hMassOmegaMinus"), casc.mOmega(), casc.pt(), casc.eta());
        } else {
          registry.fill(HIST("hMassXiPlus"), casc.mXi(), casc.pt(), casc.eta());
          registry.fill(HIST("hMassOmegaPlus"), casc.mOmega(), casc.pt(), casc.eta());
        }
        break;
      case 3:                                       // only Omega
        registry.fill(HIST("hSelectionStatus"), 8); // passes bach PID
        if (casc.sign() < 0) {
          registry.fill(HIST("hMassOmegaMinus"), casc.mOmega(), casc.pt(), casc.eta());
        } else {
          registry.fill(HIST("hMassOmegaPlus"), casc.mOmega(), casc.pt(), casc.eta());
        }
        break;
    }

    return flag;

  } // processCandidate

  void processGenMC(aod::McCollision const& mcCollision, soa::SmallGroups<soa::Join<aod::McCollisionLabels, MyCollisions>> const& collisions, aod::McParticles const& mcParticles)
  {
    // evsel
    if (INEL >= 0 && !pwglf::isINELgtNmc(mcParticles, INEL, pdgDB))
      return;
    if (std::abs(mcCollision.posZ()) > maxVertexZ)
      return;

    registry.fill(HIST("gen/hNevents"), 0);

    for (auto const& mcPart : mcParticles) {
      if (!mcPart.isPhysicalPrimary())
        continue;
      if (std::abs(mcPart.eta()) > etaCascades)
        continue;

      switch (mcPart.pdgCode()) {
        case kXiMinus:
          registry.fill(HIST("gen/hXiMinus"), mcPart.pt(), mcPart.eta());
          break;
        case kXiPlusBar:
          registry.fill(HIST("gen/hXiPlus"), mcPart.pt(), mcPart.eta());
          break;
        case kOmegaMinus:
          registry.fill(HIST("gen/hOmegaMinus"), mcPart.pt(), mcPart.eta());
          break;
        case kOmegaPlusBar:
          registry.fill(HIST("gen/hOmegaPlus"), mcPart.pt(), mcPart.eta());
          break;
      }
    }

    // Do the same thing, but now making sure there is at least one matched reconstructed event:
    if (collisions.size() < 1) {
      return;
    } else {
      bool evSel = false; // will be true if at least one rec. collision passes evsel
      for (auto const& collision : collisions) {
        // can be more than 1 rec. collisions due to event splitting
        evSel = eventSelection(collision, false);
        if (evSel) // exit loop if we find 1 rec. event that passes evsel
          break;
      }
      if (evSel) {
        // N gen events with a reconstructed event
        registry.fill(HIST("genwithrec/hNevents"), 0);

        for (auto const& mcPart : mcParticles) {
          if (!mcPart.isPhysicalPrimary())
            continue;
          if (std::abs(mcPart.eta()) > etaCascades)
            continue;

          switch (mcPart.pdgCode()) {
            case kXiMinus:
              registry.fill(HIST("genwithrec/hXiMinus"), mcPart.pt(), mcPart.eta());
              break;
            case kXiPlusBar:
              registry.fill(HIST("genwithrec/hXiPlus"), mcPart.pt(), mcPart.eta());
              break;
            case kOmegaMinus:
              registry.fill(HIST("genwithrec/hOmegaMinus"), mcPart.pt(), mcPart.eta());
              break;
            case kOmegaPlusBar:
              registry.fill(HIST("genwithrec/hOmegaPlus"), mcPart.pt(), mcPart.eta());
              break;
          }
        }
      }
    }
  } // processGen

  // wrappers for data/MC processes on reco level
  void processRecData(MyCollisions::iterator const& collision, aod::CascDataExt const& Cascades, FullTracksExtIUWithPID const&, aod::BCsWithTimestamps const&)
  {
    bool evSel = eventSelection(collision, true);
    // do not skip the collision if event selection fails - this will lead to the cascadeFlag table having less entries than the Cascade table, and therefor not joinable.
    for (auto const& casc : Cascades) {
      if (!evSel) {
        cascflags(0);
        continue;
      }
      int flag = processCandidate(casc, collision);
      cascflags(flag);
    }
  }

  void processRecMC(MyCollisions::iterator const& collision, LabeledCascades const& Cascades, FullTracksExtIUWithPID const&, aod::BCsWithTimestamps const&, aod::McParticles const&)
  {
    bool evSel = eventSelection(collision, true);
    // do not skip the collision if event selection fails - this will lead to the cascadeFlag table having less entries than the Cascade table, and therefor not joinable.
    for (auto const& casc : Cascades) {
      if (!evSel) {
        cascflags(0);
        continue;
      }
      int flag = processCandidate(casc, collision);
      cascflags(flag);
      // do mc matching here
      fillMatchedHistos(casc, flag, collision); // if sign < 0 then pdg > 0
    }
  }

  PROCESS_SWITCH(CascadeSelector, processRecData, "Process rec data", true);
  PROCESS_SWITCH(CascadeSelector, processRecMC, "Process rec MC", false);
  PROCESS_SWITCH(CascadeSelector, processGenMC, "Process gen MC", false);
}; // struct

struct CascadeCorrelations {
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  OutputObj<ZorroSummary> zorroSummary{"zorroSummary"};

  // Configurables
  Configurable<float> maxRapidity{"maxRapidity", 0.5, "|y| < maxRapidity"};
  Configurable<float> zVertexCut{"zVertexCut", 10, "Cut on PV position"};
  Configurable<int> nMixedEvents{"nMixedEvents", 10, "Number of events to be mixed"};
  Configurable<bool> doEfficiencyCorrection{"doEfficiencyCorrection", true, "flag to do efficiency corrections"};
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "CCDB url"};
  Configurable<bool> useTrigger{"useTrigger", false, "Use trigger selection on skimmed data"};
  Configurable<std::string> triggerList{"triggerList", "fDoubleXi, fDoubleOmega, fOmegaXi", "List of triggers used to select events"};
  Configurable<std::string> efficiencyCCDBPath{"efficiencyCCDBPath", "Users/r/rspijker/test/EffTest", "Path of the efficiency corrections"};
  Configurable<bool> doTFBorderCut{"doTFBorderCut", true, "Switch to apply TimeframeBorderCut event selection"};
  Configurable<bool> doSel8{"doSel8", true, "Switch to apply sel8 event selection"};
  Configurable<int> INEL{"INEL", 0, "Number of charged tracks within |eta| < 1 has to be greater than value"}; // used in MC closure
  Configurable<float> etaGenCascades{"etaGenCascades", 0.8, "min/max of eta for generated cascades"};

  ConfigurableAxis radiusAxis = {"radiusAxis", {100, 0.0f, 50.0f}, "cm"};
  ConfigurableAxis cpaAxis = {"cpaAxis", {100, 0.95f, 1.0f}, "CPA"};
  ConfigurableAxis invMassAxis = {"invMassAxis", {1000, 1.0f, 2.0f}, "Inv. Mass (GeV/c^{2})"};
  ConfigurableAxis deltaPhiAxis = {"deltaPhiAxis", {180, -PIHalf, 3 * PIHalf}, "#Delta#varphi"}; // 180 is divisible by 18 (tpc sectors) and 20 (run 2 binning)
  ConfigurableAxis ptAxis = {"ptAxis", {150, 0, 15}, "#it{p}_{T}"};
  ConfigurableAxis vertexAxis = {"vertexAxis", {200, -10.0f, 10.0f}, "cm"};
  ConfigurableAxis dcaAxis = {"dcaAxis", {100, 0.0f, 2.0f}, "cm"};
  ConfigurableAxis multiplicityAxis{"multiplicityAxis", {100, 0, 100}, "Multiplicity (centFT0M?)"};
  ConfigurableAxis invLambdaMassAxis{"invLambdaMassAxis", {100, 1.07f, 1.17f}, "Inv. Mass (GeV/c^{2})"};
  AxisSpec signAxis{3, -1.5, 1.5, "sign of cascade"};
  AxisSpec deltaYAxis{40, -2.f, 2.f, "#Delta y"};
  AxisSpec etaAxis{100, -1.f, 1.f, "#eta"};
  AxisSpec rapidityAxis{100, -1.f, 1.f, "y"};
  AxisSpec selectionFlagAxis{4, -0.5f, 3.5f, "Selection flag of casc candidate"};
  AxisSpec itsClustersAxis{8, -0.5, 7.5, "number of ITS clusters"};
  AxisSpec tpcRowsAxis{160, -0.5, 159.5, "TPC crossed rows"};

  // initialize efficiency maps
  TH1D* hEffXiMin;
  TH1D* hEffXiPlus;
  TH1D* hEffOmegaMin;
  TH1D* hEffOmegaPlus;

  // used in MC closure test
  Service<o2::framework::O2DatabasePDG> pdgDB;
  o2::pwglf::ParticleCounter<o2::framework::O2DatabasePDG> mCounter;

  void init(InitContext const&)
  {
    ccdb->setURL(ccdbUrl);
    ccdb->setCaching(true);
    if (doEfficiencyCorrection) {
      TList* effList = ccdb->getForTimeStamp<TList>(efficiencyCCDBPath, 1);
      if (!effList) {
        LOGF(fatal, "null ptr in efficiency list!");
      }
      hEffXiMin = static_cast<TH1D*>(effList->FindObject("hXiMinEff"));
      hEffXiPlus = static_cast<TH1D*>(effList->FindObject("hXiPlusEff"));
      hEffOmegaMin = static_cast<TH1D*>(effList->FindObject("hOmegaMinEff"));
      hEffOmegaPlus = static_cast<TH1D*>(effList->FindObject("hOmegaPlusEff"));
    }

    zorroSummary.setObject(zorro.getZorroSummary());

    mCounter.mPdgDatabase = pdgDB.service;
    mCounter.mSelectPrimaries = true;
  }

  double getEfficiency(TH1* h, double pT, double eta = 0)
  {
    // This function returns 1 / eff
    double eff = h->GetBinContent(h->FindFixBin(pT, eta));
    if (eff == 0)
      return 0;
    else
      return 1. / eff;
  }

  bool autoCorrelation(std::array<int, 3> triggerTracks, std::array<int, 3> assocTracks)
  {
    // function that loops over 2 arrays of track indices, checking for common elements
    for (const int& triggerTrack : triggerTracks) {
      for (const int& assocTrack : assocTracks) {
        if (triggerTrack == assocTrack)
          return true;
      }
    }
    return false;
  }

  template <typename TCascade, typename TCollision>
  void doSameEventCorrelation(const TCascade& trigger, const TCascade& assoc, const TCollision& collision){
    // autocorrelation check
    std::array<int, 3> triggerTracks = {trigger.posTrackId(), trigger.negTrackId(), trigger.bachelorId()};
    std::array<int, 3> assocTracks = {assoc.posTrackId(), assoc.negTrackId(), assoc.bachelorId()};
    if (autoCorrelation(triggerTracks, assocTracks))
      return;

    // calculate angular correlations
    double dphi = RecoDecay::constrainAngle(trigger.phi() - assoc.phi(), -PIHalf);

    double invMassXiTrigg = trigger.mXi();
    double invMassOmTrigg = trigger.mOmega();
    double invMassXiAssoc = assoc.mXi();
    double invMassOmAssoc = assoc.mOmega();

    double weightTrigg = 1.;
    double weightAssoc = 1.;

    if (trigger.isSelected() <= 2 && std::abs(trigger.yXi()) < maxRapidity) { // trigger Xi
      if (doEfficiencyCorrection)
        weightTrigg = trigger.sign() < 0 ? getEfficiency(hEffXiMin, trigger.pt(), trigger.eta()) : getEfficiency(hEffXiPlus, trigger.pt(), trigger.eta());
      if (assoc.isSelected() <= 2 && std::abs(assoc.yXi()) < maxRapidity) { // assoc Xi
        if (doEfficiencyCorrection)
          weightAssoc = assoc.sign() < 0 ? getEfficiency(hEffXiMin, assoc.pt(), assoc.eta()) : getEfficiency(hEffXiPlus, assoc.pt(), assoc.eta());
        registry.fill(HIST("hXiXi"), dphi, trigger.yXi() - assoc.yXi(), trigger.sign(), assoc.sign(), trigger.pt(), assoc.pt(), invMassXiTrigg, invMassXiAssoc, collision.posZ(), collision.centFT0M(), weightTrigg * weightAssoc);
      }
      if (assoc.isSelected() >= 2 && std::abs(assoc.yOmega()) < maxRapidity) { // assoc Omega
        if (doEfficiencyCorrection)
          weightAssoc = assoc.sign() < 0 ? getEfficiency(hEffOmegaMin, assoc.pt(), assoc.eta()) : getEfficiency(hEffOmegaPlus, assoc.pt(), assoc.eta());
        registry.fill(HIST("hXiOm"), dphi, trigger.yXi() - assoc.yOmega(), trigger.sign(), assoc.sign(), trigger.pt(), assoc.pt(), invMassXiTrigg, invMassOmAssoc, collision.posZ(), collision.centFT0M(), weightTrigg * weightAssoc);
      }
    }
    if (trigger.isSelected() >= 2 && std::abs(trigger.yOmega()) < maxRapidity) { // trigger Omega
      if (doEfficiencyCorrection)
        weightTrigg = trigger.sign() < 0 ? getEfficiency(hEffOmegaMin, trigger.pt(), trigger.eta()) : getEfficiency(hEffOmegaPlus, trigger.pt(), trigger.eta());
      if (assoc.isSelected() <= 2 && std::abs(assoc.yXi()) < maxRapidity) { // assoc Xi
        if (doEfficiencyCorrection)
          weightAssoc = assoc.sign() < 0 ? getEfficiency(hEffXiMin, assoc.pt(), assoc.eta()) : getEfficiency(hEffXiPlus, assoc.pt(), assoc.eta());
        // if Omega-Xi, fill the Xi-Omega histogram (flip the trigger/assoc and dphy,dy signs)
        registry.fill(HIST("hXiOm"), RecoDecay::constrainAngle(assoc.phi() - trigger.phi(), -PIHalf), -(trigger.yOmega() - assoc.yXi()), assoc.sign(), trigger.sign(), assoc.pt(), trigger.pt(), invMassXiAssoc, invMassOmTrigg, collision.posZ(), collision.centFT0M(), weightTrigg * weightAssoc);
      }
      if (assoc.isSelected() >= 2 && std::abs(assoc.yOmega()) < maxRapidity) { // assoc Omega
        if (doEfficiencyCorrection)
          weightAssoc = assoc.sign() < 0 ? getEfficiency(hEffOmegaMin, assoc.pt(), assoc.eta()) : getEfficiency(hEffOmegaPlus, assoc.pt(), assoc.eta());
        registry.fill(HIST("hOmOm"), dphi, trigger.yOmega() - assoc.yOmega(), trigger.sign(), assoc.sign(), trigger.pt(), assoc.pt(), invMassOmTrigg, invMassOmAssoc, collision.posZ(), collision.centFT0M(), weightTrigg * weightAssoc);
      }
    }

    // QA plots
    if (trigger.sign() * assoc.sign() < 0) {
      registry.fill(HIST("hDeltaPhiOS"), dphi);
    } else {
      registry.fill(HIST("hDeltaPhiSS"), dphi);
    }
  }

  template <typename TCascade, typename TCollision>
  void doMixedEventCorrelation(const TCascade& trigger, const TCascade& assoc, const TCollision& col1){
    if (trigger.collisionId() == assoc.collisionId()) {
      registry.fill(HIST("hMEQA"), 1.5);
      return;
    }

    std::array<int, 3> triggerTracks = {trigger.posTrackId(), trigger.negTrackId(), trigger.bachelorId()};
    std::array<int, 3> assocTracks = {assoc.posTrackId(), assoc.negTrackId(), assoc.bachelorId()};
    if (autoCorrelation(triggerTracks, assocTracks))
      return;

    double dphi = RecoDecay::constrainAngle(trigger.phi() - assoc.phi(), -PIHalf);

    double invMassXiTrigg = trigger.mXi();
    double invMassOmTrigg = trigger.mOmega();
    double invMassXiAssoc = assoc.mXi();
    double invMassOmAssoc = assoc.mOmega();

    double weightTrigg = 1.;
    double weightAssoc = 1.;

    if (trigger.isSelected() <= 2 && std::abs(trigger.yXi()) < maxRapidity) { // trigger Xi
      if (doEfficiencyCorrection)
        weightTrigg = trigger.sign() < 0 ? getEfficiency(hEffXiMin, trigger.pt(), trigger.eta()) : getEfficiency(hEffXiPlus, trigger.pt(), trigger.eta());
      if (assoc.isSelected() <= 2 && std::abs(assoc.yXi()) < maxRapidity) { // assoc Xi
        if (doEfficiencyCorrection)
          weightAssoc = assoc.sign() < 0 ? getEfficiency(hEffXiMin, assoc.pt(), assoc.eta()) : getEfficiency(hEffXiPlus, assoc.pt(), assoc.eta());
        registry.fill(HIST("MixedEvents/hMEXiXi"), dphi, trigger.yXi() - assoc.yXi(), trigger.sign(), assoc.sign(), trigger.pt(), assoc.pt(), invMassXiTrigg, invMassXiAssoc, col1.posZ(), col1.centFT0M(), weightTrigg * weightAssoc);
      }
      if (assoc.isSelected() >= 2 && std::abs(assoc.yOmega()) < maxRapidity) { // assoc Omega
        if (doEfficiencyCorrection)
          weightAssoc = assoc.sign() < 0 ? getEfficiency(hEffOmegaMin, assoc.pt(), assoc.eta()) : getEfficiency(hEffOmegaPlus, assoc.pt(), assoc.eta());
        registry.fill(HIST("MixedEvents/hMEXiOm"), dphi, trigger.yXi() - assoc.yOmega(), trigger.sign(), assoc.sign(), trigger.pt(), assoc.pt(), invMassXiTrigg, invMassOmAssoc, col1.posZ(), col1.centFT0M(), weightTrigg * weightAssoc);
      }
    }
    if (trigger.isSelected() >= 2 && std::abs(trigger.yOmega()) < maxRapidity) { // trigger Omega
      if (doEfficiencyCorrection)
        weightTrigg = trigger.sign() < 0 ? getEfficiency(hEffOmegaMin, trigger.pt(), trigger.eta()) : getEfficiency(hEffOmegaPlus, trigger.pt(), trigger.eta());
      if (assoc.isSelected() <= 2 && std::abs(assoc.yXi()) < maxRapidity) { // assoc Xi
        if (doEfficiencyCorrection)
          weightAssoc = assoc.sign() < 0 ? getEfficiency(hEffXiMin, assoc.pt(), assoc.eta()) : getEfficiency(hEffXiPlus, assoc.pt(), assoc.eta());
        // if Omega-Xi, fill the Xi-Omega histogram (flip the trigger/assoc and dphy,dy signs)
        registry.fill(HIST("MixedEvents/hMEXiOm"), RecoDecay::constrainAngle(assoc.phi() - trigger.phi(), -PIHalf), -(trigger.yOmega() - assoc.yXi()), assoc.sign(), trigger.sign(), assoc.pt(), trigger.pt(), invMassXiAssoc, invMassOmTrigg, col1.posZ(), col1.centFT0M(), weightTrigg * weightAssoc);
      }
      if (assoc.isSelected() >= 2 && std::abs(assoc.yOmega()) < maxRapidity) { // assoc Omega
        if (doEfficiencyCorrection)
          weightAssoc = assoc.sign() < 0 ? getEfficiency(hEffOmegaMin, assoc.pt(), assoc.eta()) : getEfficiency(hEffOmegaPlus, assoc.pt(), assoc.eta());
        registry.fill(HIST("MixedEvents/hMEOmOm"), dphi, trigger.yOmega() - assoc.yOmega(), trigger.sign(), assoc.sign(), trigger.pt(), assoc.pt(), invMassOmTrigg, invMassOmAssoc, col1.posZ(), col1.centFT0M(), weightTrigg * weightAssoc);
      }
    }

    // QA plots
    if (trigger.sign() * assoc.sign() < 0) {
      registry.fill(HIST("MixedEvents/hMEDeltaPhiOS"), dphi);
    } else {
      registry.fill(HIST("MixedEvents/hMEDeltaPhiSS"), dphi);
    }
  }

  template <typename TCascade>
  void doMCCorrelation(const TCascade& trigger, const TCascade& assoc, double vtxz, int FT0mult){
    if (!trigger.isPhysicalPrimary() || !assoc.isPhysicalPrimary())
      return; // require the cascades to be primaries
    if (std::abs(trigger.eta()) > etaGenCascades)
      return; // only apply eta cut to trigger - trigger normalization still valid without introducing 2-particle-acceptance effects

    double dphi = RecoDecay::constrainAngle(trigger.phi() - assoc.phi(), -PIHalf);

    if (trigger.pdgCode() < 0) { // anti-trigg --> Plus
      if (assoc.pdgCode() < 0) { // anti-assoc --> Plus
        registry.fill(HIST("MC/hMCPlusPlus"), dphi, trigger.y() - assoc.y(), trigger.pt(), assoc.pt(), vtxz, FT0mult);
      } else { // assoc --> Minus
        registry.fill(HIST("MC/hMCPlusMinus"), dphi, trigger.y() - assoc.y(), trigger.pt(), assoc.pt(), vtxz, FT0mult);
      }
    } else {                     // trig --> Minus
      if (assoc.pdgCode() < 0) { // anti-assoc --> Plus
        registry.fill(HIST("MC/hMCMinusPlus"), dphi, trigger.y() - assoc.y(), trigger.pt(), assoc.pt(), vtxz, FT0mult);
      } else {
        registry.fill(HIST("MC/hMCMinusMinus"), dphi, trigger.y() - assoc.y(), trigger.pt(), assoc.pt(), vtxz, FT0mult);
      }
    }
  }

  HistogramRegistry registry{
    "registry",
    {
      // inv mass
      {"hMassXiMinus", "hMassXiMinus", {HistType::kTH3F, {{200, 1.24, 1.44, "Inv. Mass (GeV/c^{2})"}, ptAxis, etaAxis}}},
      {"hMassXiPlus", "hMassXiPlus", {HistType::kTH3F, {{200, 1.24, 1.44, "Inv. Mass (GeV/c^{2})"}, ptAxis, etaAxis}}},
      {"hMassOmegaMinus", "hMassOmegaMinus", {HistType::kTH3F, {{200, 1.6, 1.8, "Inv. Mass (GeV/c^{2})"}, ptAxis, etaAxis}}},
      {"hMassOmegaPlus", "hMassOmegaPlus", {HistType::kTH3F, {{200, 1.6, 1.8, "Inv. Mass (GeV/c^{2})"}, ptAxis, etaAxis}}},
      // efficiency corrected inv mass
      {"hMassXiEffCorrected", "hMassXiEffCorrected", {HistType::kTHnSparseF, {invMassAxis, signAxis, ptAxis, etaAxis, vertexAxis, multiplicityAxis}}, true},
      {"hMassOmegaEffCorrected", "hMassOmegaEffCorrected", {HistType::kTHnSparseF, {invMassAxis, signAxis, ptAxis, etaAxis, vertexAxis, multiplicityAxis}}, true},

      // trigger QA
      {"hTriggerQA", "hTriggerQA", {HistType::kTH1F, {{2, -0.5, 1.5, "Trigger y/n"}}}},

      // basic selection variables (after cuts)
      {"hV0Radius", "hV0Radius", {HistType::kTH1F, {radiusAxis}}},
      {"hCascRadius", "hCascRadius", {HistType::kTH1F, {radiusAxis}}},
      {"hV0CosPA", "hV0CosPA", {HistType::kTH1F, {cpaAxis}}},
      {"hCascCosPA", "hCascCosPA", {HistType::kTH1F, {cpaAxis}}},
      {"hDCAPosToPV", "hDCAPosToPV", {HistType::kTH1F, {vertexAxis}}},
      {"hDCANegToPV", "hDCANegToPV", {HistType::kTH1F, {vertexAxis}}},
      {"hDCABachToPV", "hDCABachToPV", {HistType::kTH1F, {vertexAxis}}},
      {"hDCAV0ToPV", "hDCAV0ToPV", {HistType::kTH1F, {vertexAxis}}},
      {"hDCAV0Dau", "hDCAV0Dau", {HistType::kTH1F, {dcaAxis}}},
      {"hDCACascDau", "hDCACascDau", {HistType::kTH1F, {dcaAxis}}},
      {"hLambdaMass", "hLambdaMass", {HistType::kTH1F, {invLambdaMassAxis}}},
      {"hTPCnCrossedRowsPos", "hTPCnCrossedRowsPos", {HistType::kTH1F, {tpcRowsAxis}}},
      {"hTPCnCrossedRowsNeg", "hTPCnCrossedRowsNeg", {HistType::kTH1F, {tpcRowsAxis}}},
      {"hTPCnCrossedRowsBach", "hTPCnCrossedRowsBach", {HistType::kTH1F, {tpcRowsAxis}}},
      {"hITSnClustersPos", "hITSnClustersPos", {HistType::kTH1F, {itsClustersAxis}}},
      {"hITSnClustersNeg", "hITSnClustersNeg", {HistType::kTH1F, {itsClustersAxis}}},
      {"hITSnClustersBach", "hITSnClustersBach", {HistType::kTH1F, {itsClustersAxis}}},

      {"hSelectionFlag", "hSelectionFlag", {HistType::kTH1I, {selectionFlagAxis}}},
      {"hAutoCorrelation", "hAutoCorrelation", {HistType::kTH1I, {{4, -0.5f, 3.5f, "Types of SS autocorrelation"}}}},
      {"hAutoCorrelationOS", "hAutoCorrelationOS", {HistType::kTH1I, {{2, -1.f, 1.f, "Charge of OS autocorrelated track"}}}},
      {"hPhi", "hPhi", {HistType::kTH1F, {{180, 0, TwoPI, "#varphi"}}}},
      {"hEta", "hEta", {HistType::kTH1F, {{100, -2, 2, "#eta"}}}},
      {"hRapidityXi", "hRapidityXi", {HistType::kTH1F, {rapidityAxis}}},
      {"hRapidityOmega", "hRapidityOmega", {HistType::kTH1F, {rapidityAxis}}},

      // correlation histos
      {"hDeltaPhiSS", "hDeltaPhiSS", {HistType::kTH1F, {deltaPhiAxis}}},
      {"hDeltaPhiOS", "hDeltaPhiOS", {HistType::kTH1F, {deltaPhiAxis}}},

      {"hXiXi", "hXiXi", {HistType::kTHnSparseF, {deltaPhiAxis, deltaYAxis, signAxis, signAxis, ptAxis, ptAxis, invMassAxis, invMassAxis, vertexAxis, multiplicityAxis}}, true},
      {"hXiOm", "hXiOm", {HistType::kTHnSparseF, {deltaPhiAxis, deltaYAxis, signAxis, signAxis, ptAxis, ptAxis, invMassAxis, invMassAxis, vertexAxis, multiplicityAxis}}, true},
      {"hOmOm", "hOmOm", {HistType::kTHnSparseF, {deltaPhiAxis, deltaYAxis, signAxis, signAxis, ptAxis, ptAxis, invMassAxis, invMassAxis, vertexAxis, multiplicityAxis}}, true},

      // Mixed events
      {"MixedEvents/hMEVz1", "hMEVz1", {HistType::kTH1F, {vertexAxis}}},
      {"MixedEvents/hMEVz2", "hMEVz2", {HistType::kTH1F, {vertexAxis}}},
      {"MixedEvents/hMEDeltaPhiSS", "hMEDeltaPhiSS", {HistType::kTH1F, {deltaPhiAxis}}},
      {"MixedEvents/hMEDeltaPhiOS", "hMEDeltaPhiOS", {HistType::kTH1F, {deltaPhiAxis}}},
      {"MixedEvents/hMEQA", "hMEQA", {HistType::kTH1I, {{2, 0, 2, "QA for exceptions in ME (this histogram should have 0 entries!)"}}}},
      {"MixedEvents/hMEAutoCorrelation", "hMEAutoCorrelation", {HistType::kTH1I, {{4, -0.5f, 3.5f, "Types of SS autocorrelation"}}}},
      {"MixedEvents/hMEAutoCorrelationOS", "hMEAutoCorrelationOS", {HistType::kTH1I, {{2, -1.f, 1.f, "Charge of OS autocorrelated track"}}}},

      {"MixedEvents/hMEXiXi", "hMEXiXi", {HistType::kTHnSparseF, {deltaPhiAxis, deltaYAxis, signAxis, signAxis, ptAxis, ptAxis, invMassAxis, invMassAxis, vertexAxis, multiplicityAxis}}, true},
      {"MixedEvents/hMEXiOm", "hMEXiOm", {HistType::kTHnSparseF, {deltaPhiAxis, deltaYAxis, signAxis, signAxis, ptAxis, ptAxis, invMassAxis, invMassAxis, vertexAxis, multiplicityAxis}}, true},
      {"MixedEvents/hMEOmOm", "hMEOmOm", {HistType::kTHnSparseF, {deltaPhiAxis, deltaYAxis, signAxis, signAxis, ptAxis, ptAxis, invMassAxis, invMassAxis, vertexAxis, multiplicityAxis}}, true},

      // MC closure
      {"MC/hMCPlusMinus", "hMCPlusMinus", {HistType::kTHnSparseF, {deltaPhiAxis, deltaYAxis, ptAxis, ptAxis, vertexAxis, multiplicityAxis}}, true},
      {"MC/hMCPlusPlus", "hMCPlusPlus", {HistType::kTHnSparseF, {deltaPhiAxis, deltaYAxis, ptAxis, ptAxis, vertexAxis, multiplicityAxis}}, true},
      {"MC/hMCMinusPlus", "hMCMinusPlus", {HistType::kTHnSparseF, {deltaPhiAxis, deltaYAxis, ptAxis, ptAxis, vertexAxis, multiplicityAxis}}, true},
      {"MC/hMCMinusMinus", "hMCMinusMinus", {HistType::kTHnSparseF, {deltaPhiAxis, deltaYAxis, ptAxis, ptAxis, vertexAxis, multiplicityAxis}}, true},

      {"MC/hGenMultNoReco", "hGenMultNoReco", {HistType::kTH1I, {{100, 0, 100, "Number of generated charged primaries"}}}},
      {"MC/hGenMultOneReco", "hGenMultOneReco", {HistType::kTH1I, {{100, 0, 100, "Number of generated charged primaries"}}}},
      {"MC/hSplitEvents", "hSplitEvents", {HistType::kTH1I, {{10, 0, 10, "Number of rec. events per gen event"}}}},

      // debug
      {"MC/hPhi", "hPhi", {HistType::kTH1F, {{180, 0, TwoPI}}}},
      {"MC/hEta", "hEta", {HistType::kTH1F, {{100, -2, 2}}}},
      {"MC/hRapidity", "hRapidity", {HistType::kTH1F, {{100, -2, 2}}}},
    },
  };

  // cascade filter
  Filter cascadeSelector = aod::cascadeflags::isSelected > 0;

  // Warning: it is not possible to use this axis as configurable due to a bug - however, default values are sensible.
  SliceCache cache;
  ConfigurableAxis axisVtxZ{"axisVtxZ", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
  // ConfigurableAxis axisMult{"axisMult", {VARIABLE_WIDTH, 0, 5, 10, 20, 30, 40, 50, 100, 1000}, "Mixing bins - multiplicity"};
  // using BinningType = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::centFT0M>;
  // BinningType colBinning{{axisVtxZ, axisMult}, true}; // true is for 'ignore overflows' (true by default). Underflows and overflows will have bin -1.
  using BinningType = ColumnBinningPolicy<aod::collision::PosZ>;
  BinningType colBinning{{axisVtxZ}, true}; // true is for 'ignore overflows' (true by default). Underflows and overflows will have bin -1.

  void processSameEvent(MyCollisionsMult::iterator const& collision, MyCascades const& Cascades, FullTracksExtIU const&, aod::BCsWithTimestamps const&)
  {
    if (useTrigger) {
      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      zorro.initCCDB(ccdb.service, bc.runNumber(), bc.timestamp(), triggerList);
      bool eventTrigger = zorro.isSelected(bc.globalBC());
      if (eventTrigger) {
        registry.fill(HIST("hTriggerQA"), 1);
      } else {
        registry.fill(HIST("hTriggerQA"), 0);
        return;
      }
    }

    double weight;
    // Some QA on the cascades
    for (auto const& casc : Cascades) {
      if (casc.isSelected() <= 2) { // not exclusively an Omega --> consistent with Xi or both
        if (casc.sign() < 0) {
          registry.fill(HIST("hMassXiMinus"), casc.mXi(), casc.pt(), casc.eta());
          weight = getEfficiency(hEffXiMin, casc.pt(), casc.eta());
        } else {
          registry.fill(HIST("hMassXiPlus"), casc.mXi(), casc.pt(), casc.eta());
          weight = getEfficiency(hEffXiPlus, casc.pt(), casc.eta());
        }
        // LOGF(info, "casc pt %f, weight %f", casc.pt(), weight);
        registry.fill(HIST("hMassXiEffCorrected"), casc.mXi(), casc.sign(), casc.pt(), casc.eta(), collision.posZ(), collision.centFT0M(), weight);
        registry.fill(HIST("hRapidityXi"), casc.yXi());
      }
      if (casc.isSelected() >= 2) { // consistent with Omega or both
        if (casc.sign() < 0) {
          registry.fill(HIST("hMassOmegaMinus"), casc.mOmega(), casc.pt(), casc.eta());
          weight = getEfficiency(hEffOmegaMin, casc.pt(), casc.eta());
        } else {
          registry.fill(HIST("hMassOmegaPlus"), casc.mOmega(), casc.pt(), casc.eta());
          weight = getEfficiency(hEffOmegaPlus, casc.pt(), casc.eta());
        }
        registry.fill(HIST("hMassOmegaEffCorrected"), casc.mOmega(), casc.sign(), casc.pt(), casc.eta(), collision.posZ(), collision.centFT0M(), weight);
        registry.fill(HIST("hRapidityOmega"), casc.yOmega());
      }
      registry.fill(HIST("hV0Radius"), casc.v0radius());
      registry.fill(HIST("hCascRadius"), casc.cascradius());
      registry.fill(HIST("hV0CosPA"), casc.v0cosPA(collision.posX(), collision.posY(), collision.posZ()));
      registry.fill(HIST("hCascCosPA"), casc.casccosPA(collision.posX(), collision.posY(), collision.posZ()));
      registry.fill(HIST("hDCAPosToPV"), casc.dcapostopv());
      registry.fill(HIST("hDCANegToPV"), casc.dcanegtopv());
      registry.fill(HIST("hDCABachToPV"), casc.dcabachtopv());
      registry.fill(HIST("hDCAV0ToPV"), casc.dcav0topv(collision.posX(), collision.posY(), collision.posZ()));
      registry.fill(HIST("hDCAV0Dau"), casc.dcaV0daughters());
      registry.fill(HIST("hDCACascDau"), casc.dcacascdaughters());
      registry.fill(HIST("hLambdaMass"), casc.mLambda());
      registry.fill(HIST("hITSnClustersPos"), casc.posTrack_as<FullTracksExtIU>().itsNCls());
      registry.fill(HIST("hITSnClustersNeg"), casc.negTrack_as<FullTracksExtIU>().itsNCls());
      registry.fill(HIST("hITSnClustersBach"), casc.bachelor_as<FullTracksExtIU>().itsNCls());
      registry.fill(HIST("hTPCnCrossedRowsPos"), casc.posTrack_as<FullTracksExtIU>().tpcNClsCrossedRows());
      registry.fill(HIST("hTPCnCrossedRowsNeg"), casc.negTrack_as<FullTracksExtIU>().tpcNClsCrossedRows());
      registry.fill(HIST("hTPCnCrossedRowsBach"), casc.bachelor_as<FullTracksExtIU>().tpcNClsCrossedRows());

      registry.fill(HIST("hSelectionFlag"), casc.isSelected());
      registry.fill(HIST("hPhi"), casc.phi());
      registry.fill(HIST("hEta"), casc.eta());
    } // casc loop

    for (const auto& [c0, c1] : combinations(Cascades, Cascades)) { // combinations automatically applies strictly upper in case of 2 identical tables
      // Define the trigger as the particle with the highest pT. 
      // As we can't swap the cascade tables themselves, we have created a function that we can call with the correct order.
      if (c0.pt() >= c1.pt()) {
        doSameEventCorrelation(c0, c1, collision);
      } else {
        doSameEventCorrelation(c1, c0, collision);
      }
    } // correlations
  } // process same event

  void processMixedEvent(MyCollisionsMult const& collisions, MyCascades const& Cascades, FullTracksExtIU const&)
  {
    auto cascadesTuple = std::make_tuple(Cascades);
    SameKindPair<MyCollisionsMult, MyCascades, BinningType> pair{colBinning, nMixedEvents, -1, collisions, cascadesTuple, &cache};

    for (auto const& [col1, cascades1, col2, cascades2] : pair) {
      if (!col1.sel8() || !col2.sel8())
        continue;
      if (std::abs(col1.posZ()) > zVertexCut || std::abs(col2.posZ()) > zVertexCut)
        continue;
      if (col1.globalIndex() == col2.globalIndex()) {
        registry.fill(HIST("hMEQA"), 0.5);
        continue;
      }
      registry.fill(HIST("MixedEvents/hMEVz1"), col1.posZ());
      registry.fill(HIST("MixedEvents/hMEVz2"), col2.posZ());

      for (const auto& [casc1, casc2] : combinations(CombinationsFullIndexPolicy(cascades1, cascades2))) {
        // specify FullIndexPolicy since the cascades are from different collisions

        // Define the trigger as the particle with the highest pT.
        // As we can't swap the cascade tables themselves, we have created a function that we can call with the correct order.
        if(casc1.pt() >= casc2.pt()) {
          doMixedEventCorrelation(casc1, casc2, col1);
        } else {
          doMixedEventCorrelation(casc2, casc1, col2);
        }
      } // correlations
    } // collisions
  } // process mixed events

  Filter genCascadesFilter = nabs(aod::mcparticle::pdgCode) == static_cast<int>(kXiMinus);

  void processMC(aod::McCollision const& mcCollision, soa::SmallGroups<soa::Join<aod::McCollisionLabels, MyCollisionsMult>> const& collisions, soa::Filtered<aod::McParticles> const& genCascades, aod::McParticles const& mcParticles)
  {
    // apply evsel
    if (INEL >= 0 && !pwglf::isINELgtNmc(mcParticles, INEL, pdgDB))
      return;
    if (std::abs(mcCollision.posZ()) > zVertexCut)
      return;

    // Let's do some logic on matched reconstructed collisions - if there less or more than one, fill some QA and skip the rest
    double FT0mult = -1; // non-sensible default value just in case
    double vtxz = mcCollision.posZ();
    if (collisions.size() < 1) {
      registry.fill(HIST("MC/hSplitEvents"), 0);
      registry.fill(HIST("MC/hGenMultNoReco"), mCounter.countFT0A(mcParticles) + mCounter.countFT0C(mcParticles));
    } else if (collisions.size() == 1) {
      registry.fill(HIST("MC/hSplitEvents"), 1);
      registry.fill(HIST("MC/hGenMultOneReco"), mCounter.countFT0A(mcParticles) + mCounter.countFT0C(mcParticles));
      for (auto const& collision : collisions) { // not really a loop, as there is only one collision
        FT0mult = collision.centFT0M();
      }
    } else if (collisions.size() > 1) {
      registry.fill(HIST("MC/hSplitEvents"), collisions.size());
    }

    // QA
    for (const auto& casc : genCascades) {
      if (!casc.isPhysicalPrimary())
        continue;
      registry.fill(HIST("MC/hPhi"), casc.phi());
      registry.fill(HIST("MC/hEta"), casc.eta());
      registry.fill(HIST("MC/hRapidity"), casc.y());
    }

    for (const auto& [c0, c1] : combinations(genCascades, genCascades)) { // combinations automatically applies strictly upper in case of 2 identical tables
      // Define the trigger as the particle with the highest pT.
      // As we can't swap the cascade tables themselves, we have created a function that we can call with the correct order.
      if (c0.pt() >= c1.pt()) {
        doMCCorrelation(c0, c1, vtxz, FT0mult);
      } else {
        doMCCorrelation(c1, c0, vtxz, FT0mult);
      }
    }
  }

  PROCESS_SWITCH(CascadeCorrelations, processSameEvent, "Process same events", true);
  PROCESS_SWITCH(CascadeCorrelations, processMixedEvent, "Process mixed events", true);
  PROCESS_SWITCH(CascadeCorrelations, processMC, "Process MC", false);

}; // struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<CascadeSelector>(cfgc),
    adaptAnalysisTask<CascadeCorrelations>(cfgc)};
}
