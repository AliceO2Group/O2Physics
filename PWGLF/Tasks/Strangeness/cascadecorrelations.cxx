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

#include <Math/Vector4D.h>
#include <cmath>
#include <array>
#include <cstdlib>
#include <string>
#include <utility>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponse.h"
#include "CCDB/BasicCCDBManager.h"
#include "EventFiltering/Zorro.h"

#include <TFile.h>
#include <TList.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TLorentzVector.h>
#include <TPDGCode.h>
// #include <TDatabasePDG.h>

using namespace o2;
using namespace o2::soa;
using namespace o2::constants::math;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;

// use parameters + cov mat non-propagated, aux info + (extension propagated)
using FullTracksExt = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCov, aod::TracksDCA>;
using FullTracksExtIU = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::TracksDCA>;
using FullTracksExtWithPID = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCov, aod::TracksDCA, aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr>;
using FullTracksExtIUWithPID = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::TracksDCA, aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr>;

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
using MyCollisionsMult = soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults>;
using MyCascades = soa::Filtered<aod::CascDataExtSelected>;

struct CascadeSelector {
  Service<o2::ccdb::BasicCCDBManager> ccdb;

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

  // Tracklevel
  Configurable<float> tpcNsigmaBachelor{"tpcNsigmaBachelor", 3, "TPC NSigma bachelor"};
  Configurable<float> tpcNsigmaProton{"tpcNsigmaProton", 3, "TPC NSigma proton <- lambda"};
  Configurable<float> tpcNsigmaPion{"tpcNsigmaPion", 3, "TPC NSigma pion <- lambda"};
  Configurable<int> minTPCCrossedRows{"minTPCCrossedRows", 80, "min N TPC crossed rows"}; // TODO: finetune! 80 > 159/2, so no split tracks?
  Configurable<int> minITSClusters{"minITSClusters", 4, "minimum number of ITS clusters"};
  Configurable<float> etaTracks{"etaTracks", 1.0, "min/max of eta for tracks"};

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
  AxisSpec vertexAxis = {100, -10.0f, 10.0f, "cm"};
  AxisSpec dcaAxis = {50, 0.0f, 5.0f, "cm"};
  // AxisSpec invMassAxis = {1000, 1.0f, 2.0f, "Inv. Mass (GeV/c^{2})"};
  AxisSpec invXiMassAxis = {100, 1.28f, 1.38f, "Inv. Mass (GeV/c^{2})"};
  AxisSpec invOmegaMassAxis = {100, 1.62f, 1.72f, "Inv. Mass (GeV/c^{2})"};
  AxisSpec ptAxis = {150, 0, 15, "#it{p}_{T}"};
  AxisSpec rapidityAxis{100, -1.f, 1.f, "y"};
  HistogramRegistry registry{
    "registry",
    {
      // basic selection variables
      {"hV0Radius", "hV0Radius", {HistType::kTH3F, {{100, 0.0f, 100.0f, "cm"}, invXiMassAxis, ptAxis}}},
      {"hCascRadius", "hCascRadius", {HistType::kTH3F, {{100, 0.0f, 100.0f, "cm"}, invXiMassAxis, ptAxis}}},
      {"hV0CosPA", "hV0CosPA", {HistType::kTH3F, {{100, 0.95f, 1.0f}, invXiMassAxis, ptAxis}}},
      {"hCascCosPA", "hCascCosPA", {HistType::kTH3F, {{100, 0.95f, 1.0f}, invXiMassAxis, ptAxis}}},
      {"hDCAPosToPV", "hDCAPosToPV", {HistType::kTH3F, {vertexAxis, invXiMassAxis, ptAxis}}},
      {"hDCANegToPV", "hDCANegToPV", {HistType::kTH3F, {vertexAxis, invXiMassAxis, ptAxis}}},
      {"hDCABachToPV", "hDCABachToPV", {HistType::kTH3F, {vertexAxis, invXiMassAxis, ptAxis}}},
      {"hDCAV0ToPV", "hDCAV0ToPV", {HistType::kTH3F, {vertexAxis, invXiMassAxis, ptAxis}}},
      {"hDCAV0Dau", "hDCAV0Dau", {HistType::kTH3F, {dcaAxis, invXiMassAxis, ptAxis}}},
      {"hDCACascDau", "hDCACascDau", {HistType::kTH3F, {dcaAxis, invXiMassAxis, ptAxis}}},
      {"hLambdaMass", "hLambdaMass", {HistType::kTH3F, {{100, 1.0f, 1.2f, "Inv. Mass (GeV/c^{2})"}, invXiMassAxis, ptAxis}}},

      {"hMassXiMinus", "hMassXiMinus", {HistType::kTH3F, {invXiMassAxis, ptAxis, rapidityAxis}}},
      {"hMassXiPlus", "hMassXiPlus", {HistType::kTH3F, {invXiMassAxis, ptAxis, rapidityAxis}}},
      {"hMassOmegaMinus", "hMassOmegaMinus", {HistType::kTH3F, {invOmegaMassAxis, ptAxis, rapidityAxis}}},
      {"hMassOmegaPlus", "hMassOmegaPlus", {HistType::kTH3F, {invOmegaMassAxis, ptAxis, rapidityAxis}}},

      // // invariant mass per cut, start with Xi
      // {"hMassXi0", "Xi inv mass before selections", {HistType::kTH2F, {invMassAxis, ptAxis}}},
      // {"hMassXi1", "Xi inv mass after TPCnCrossedRows cut", {HistType::kTH2F, {invMassAxis, ptAxis}}},
      // {"hMassXi2", "Xi inv mass after ITSnClusters cut", {HistType::kTH2F, {invMassAxis, ptAxis}}},
      // {"hMassXi3", "Xi inv mass after topo cuts", {HistType::kTH2F, {invMassAxis, ptAxis}}},
      // {"hMassXi4", "Xi inv mass after V0 daughters PID cut", {HistType::kTH2F, {invMassAxis, ptAxis}}},
      // {"hMassXi5", "Xi inv mass after bachelor PID cut", {HistType::kTH2F, {invMassAxis, ptAxis}}},

      // ITS & TPC clusters, with Xi inv mass
      {"hTPCnCrossedRowsPos", "hTPCnCrossedRowsPos", {HistType::kTH3F, {{160, -0.5, 159.5, "TPC crossed rows"}, invXiMassAxis, ptAxis}}},
      {"hTPCnCrossedRowsNeg", "hTPCnCrossedRowsNeg", {HistType::kTH3F, {{160, -0.5, 159.5, "TPC crossed rows"}, invXiMassAxis, ptAxis}}},
      {"hTPCnCrossedRowsBach", "hTPCnCrossedRowsBach", {HistType::kTH3F, {{160, -0.5, 159.5, "TPC crossed rows"}, invXiMassAxis, ptAxis}}},
      {"hITSnClustersPos", "hITSnClustersPos", {HistType::kTH3F, {{8, -0.5, 7.5, "number of ITS clusters"}, invXiMassAxis, ptAxis}}},
      {"hITSnClustersNeg", "hITSnClustersNeg", {HistType::kTH3F, {{8, -0.5, 7.5, "number of ITS clusters"}, invXiMassAxis, ptAxis}}},
      {"hITSnClustersBach", "hITSnClustersBach", {HistType::kTH3F, {{8, -0.5, 7.5, "number of ITS clusters"}, invXiMassAxis, ptAxis}}},

      {"hTriggerQA", "hTriggerQA", {HistType::kTH1F, {{2, -0.5, 1.5, "Trigger y/n"}}}},
    },
  };

  // Keep track of which selections the candidates pass
  void init(InitContext const&)
  {
    ccdb->setURL(ccdbUrl);
    ccdb->setCaching(true);

    auto h = registry.add<TH1>("hSelectionStatus", "hSelectionStatus", HistType::kTH1I, {{10, 0, 10, "status"}});
    h->GetXaxis()->SetBinLabel(1, "All");
    h->GetXaxis()->SetBinLabel(2, "nTPC OK");
    h->GetXaxis()->SetBinLabel(3, "nITS OK");
    h->GetXaxis()->SetBinLabel(4, "Topo OK");
    h->GetXaxis()->SetBinLabel(5, "Track eta OK");
    h->GetXaxis()->SetBinLabel(6, "V0 PID OK");
    h->GetXaxis()->SetBinLabel(7, "Bach PID OK");

    auto hEventSel = registry.add<TH1>("hEventSel", "hEventSel", HistType::kTH1I, {{10, 0, 10, "selection criteria"}});
    hEventSel->GetXaxis()->SetBinLabel(1, "All");
    hEventSel->GetXaxis()->SetBinLabel(2, "sel8");
    hEventSel->GetXaxis()->SetBinLabel(3, "INEL0");
    hEventSel->GetXaxis()->SetBinLabel(4, "V_z");
    hEventSel->GetXaxis()->SetBinLabel(5, "NoSameBunchPileUp");
    hEventSel->GetXaxis()->SetBinLabel(6, "Selected events");
  }

  void process(MyCollisions::iterator const& collision, aod::CascDataExt const& Cascades, FullTracksExtIUWithPID const&, aod::BCsWithTimestamps const&)
  {
    bool evSel = true;
    if (useTrigger) {
      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      zorro.initCCDB(ccdb.service, bc.runNumber(), bc.timestamp(), triggerList);
      bool eventTrigger = zorro.isSelected(bc.globalBC());
      if (eventTrigger) {
        registry.fill(HIST("hTriggerQA"), 1);
      } else {
        registry.fill(HIST("hTriggerQA"), 0);
        evSel = false;
      }
    }

    // fill event selection based on which selection criteria are applied and passed
    // do not skip the collision - this will lead to the cascadeFlag table having less entries than the Cascade table, and therefor not joinable.
    registry.fill(HIST("hEventSel"), 0);
    if (doSel8 && !collision.sel8()) {
      evSel = false;
      registry.fill(HIST("hEventSel"), 1);
    } else if (collision.multNTracksPVeta1() <= INEL) {
      evSel = false;
      registry.fill(HIST("hEventSel"), 2);
    } else if (std::abs(collision.posZ()) > maxVertexZ) {
      evSel = false;
      registry.fill(HIST("hEventSel"), 3);
    } else if (doNoSameBunchPileUp && !collision.selection_bit(aod::evsel::kNoSameBunchPileup)) {
      evSel = false;
      registry.fill(HIST("hEventSel"), 4);
    }
    if (evSel) // passes all selections
      registry.fill(HIST("hEventSel"), 5);

    for (auto const& casc : Cascades) {
      if (!evSel) {
        cascflags(0);
        continue;
      }

      // these are the tracks:
      auto bachTrack = casc.bachelor_as<FullTracksExtIUWithPID>();
      auto posTrack = casc.posTrack_as<FullTracksExtIUWithPID>();
      auto negTrack = casc.negTrack_as<FullTracksExtIUWithPID>();

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

      registry.fill(HIST("hSelectionStatus"), 0); // all the cascade before selections
      // registry.fill(HIST("hMassXi0"), casc.mXi(), casc.pt());

      // TPC N crossed rows todo: check if minTPCCrossedRows > 50
      if (posTrack.tpcNClsCrossedRows() < minTPCCrossedRows || negTrack.tpcNClsCrossedRows() < minTPCCrossedRows || bachTrack.tpcNClsCrossedRows() < minTPCCrossedRows) {
        cascflags(0);
        continue;
      }
      registry.fill(HIST("hSelectionStatus"), 1); // passes nTPC crossed rows
      // registry.fill(HIST("hMassXi1"), casc.mXi(), casc.pt());

      // ITS N clusters todo: check if minITSClusters > 0
      if (posTrack.itsNCls() < minITSClusters || negTrack.itsNCls() < minITSClusters || bachTrack.itsNCls() < minITSClusters) {
        cascflags(0);
        continue;
      }
      registry.fill(HIST("hSelectionStatus"), 2); // passes nITS clusters
      // registry.fill(HIST("hMassXi2"), casc.mXi(), casc.pt());

      //// TOPO CUTS //// TODO: improve!
      double pvx = collision.posX();
      double pvy = collision.posY();
      double pvz = collision.posZ();
      if (casc.v0radius() < v0setting_radius ||
          casc.cascradius() < cascadesetting_cascradius ||
          casc.v0cosPA(pvx, pvy, pvz) < v0setting_cospa ||
          casc.casccosPA(pvx, pvy, pvz) < cascadesetting_cospa ||
          casc.dcav0topv(pvx, pvy, pvz) < cascadesetting_mindcav0topv ||
          TMath::Abs(casc.mLambda() - 1.115683) > cascadesetting_v0masswindow) {
        // It failed at least one topo selection
        cascflags(0);
        continue;
      }
      registry.fill(HIST("hSelectionStatus"), 3); // passes topo
      // registry.fill(HIST("hMassXi3"), casc.mXi(), casc.pt());

      if (TMath::Abs(posTrack.eta()) > etaTracks || TMath::Abs(negTrack.eta()) > etaTracks || TMath::Abs(bachTrack.eta()) > etaTracks) {
        cascflags(0);
        continue;
      }
      registry.fill(HIST("hSelectionStatus"), 4); // passes track eta

      // TODO: TOF (for pT > 2 GeV per track?)

      //// TPC PID ////
      // Lambda check
      if (casc.sign() < 0) {
        // Proton check:
        if (TMath::Abs(posTrack.tpcNSigmaPr()) > tpcNsigmaProton) {
          cascflags(0);
          continue;
        }
        // Pion check:
        if (TMath::Abs(negTrack.tpcNSigmaPi()) > tpcNsigmaPion) {
          cascflags(0);
          continue;
        }
      } else {
        // Proton check:
        if (TMath::Abs(negTrack.tpcNSigmaPr()) > tpcNsigmaProton) {
          cascflags(0);
          continue;
        }
        // Pion check:
        if (TMath::Abs(posTrack.tpcNSigmaPi()) > tpcNsigmaPion) {
          cascflags(0);
          continue;
        }
      }
      registry.fill(HIST("hSelectionStatus"), 5); // passes V0 daughters PID
      // registry.fill(HIST("hMassXi4"), casc.mXi(), casc.pt());

      // Bachelor check
      if (TMath::Abs(bachTrack.tpcNSigmaPi()) < tpcNsigmaBachelor) {
        if (TMath::Abs(bachTrack.tpcNSigmaKa()) < tpcNsigmaBachelor) {
          // consistent with both!
          cascflags(2);
          registry.fill(HIST("hSelectionStatus"), 6); // passes bach PID
          // registry.fill(HIST("hMassXi5"), casc.mXi(), casc.pt());
          if (casc.sign() < 0) {
            registry.fill(HIST("hMassXiMinus"), casc.mXi(), casc.pt(), casc.yXi());
            registry.fill(HIST("hMassOmegaMinus"), casc.mOmega(), casc.pt(), casc.yOmega());
          } else {
            registry.fill(HIST("hMassXiPlus"), casc.mXi(), casc.pt(), casc.yXi());
            registry.fill(HIST("hMassOmegaPlus"), casc.mOmega(), casc.pt(), casc.yOmega());
          }
          continue;
        }
        cascflags(1);
        registry.fill(HIST("hSelectionStatus"), 6); // passes bach PID
        // registry.fill(HIST("hMassXi5"), casc.mXi(), casc.pt());
        if (casc.sign() < 0) {
          registry.fill(HIST("hMassXiMinus"), casc.mXi(), casc.pt(), casc.yXi());
        } else {
          registry.fill(HIST("hMassXiPlus"), casc.mXi(), casc.pt(), casc.yXi());
        }
        continue;
      } else if (TMath::Abs(bachTrack.tpcNSigmaKa()) < tpcNsigmaBachelor) {
        cascflags(3);
        registry.fill(HIST("hSelectionStatus"), 6); // passes bach PID
        if (casc.sign() < 0) {
          registry.fill(HIST("hMassOmegaMinus"), casc.mOmega(), casc.pt(), casc.yOmega());
        } else {
          registry.fill(HIST("hMassOmegaPlus"), casc.mOmega(), casc.pt(), casc.yOmega());
        }
        continue;
      }
      // if we reach here, the bachelor was neither pion nor kaon
      cascflags(0);
    } // cascade loop
  } // process
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

  AxisSpec invMassAxis = {1000, 1.0f, 2.0f, "Inv. Mass (GeV/c^{2})"};
  AxisSpec deltaPhiAxis = {180, -PIHalf, 3 * PIHalf, "#Delta#varphi"};       // 180 is divisible by 18 (tpc sectors) and 20 (run 2 binning)
  AxisSpec deltaYAxis = {40, -2 * maxRapidity, 2 * maxRapidity, "#Delta y"}; // TODO: narrower range?
  AxisSpec ptAxis = {150, 0, 15, "#it{p}_{T}"};
  AxisSpec selectionFlagAxis = {4, -0.5f, 3.5f, "Selection flag of casc candidate"};
  AxisSpec vertexAxis = {200, -10.0f, 10.0f, "cm"};
  AxisSpec multiplicityAxis{100, 0, 100, "Multiplicity (MultFT0M?)"};
  AxisSpec rapidityAxis{100, -maxRapidity, maxRapidity, "y"};

  // initialize efficiency maps
  TH1D* hEffXiMin;
  TH1D* hEffXiPlus;
  TH1D* hEffOmegaMin;
  TH1D* hEffOmegaPlus;

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
  }

  double getEfficiency(TH1D* h, double pT)
  { // TODO: make 2D (rapidity)
    // This function returns the value of histogram h corresponding to the x-coordinate pT
    return h->GetBinContent(h->GetXaxis()->FindFixBin(pT));
  }

  HistogramRegistry registry{
    "registry",
    {
      // inv mass
      {"hMassXiMinus", "hMassXiMinus", {HistType::kTH2F, {invMassAxis, ptAxis}}},
      {"hMassXiPlus", "hMassXiPlus", {HistType::kTH2F, {invMassAxis, ptAxis}}},
      {"hMassOmegaMinus", "hMassOmegaMinus", {HistType::kTH2F, {invMassAxis, ptAxis}}},
      {"hMassOmegaPlus", "hMassOmegaPlus", {HistType::kTH2F, {invMassAxis, ptAxis}}},
      // efficiency corrected inv mass
      {"hMassXiEffCorrected", "hMassXiEffCorrected", {HistType::kTHnSparseF, {invMassAxis, ptAxis, rapidityAxis, vertexAxis, multiplicityAxis}}, true},
      {"hMassOmegaEffCorrected", "hMassOmegaEffCorrected", {HistType::kTHnSparseF, {invMassAxis, ptAxis, rapidityAxis, vertexAxis, multiplicityAxis}}, true},

      // trigger QA
      {"hTriggerQA", "hTriggerQA", {HistType::kTH1F, {{2, -0.5, 1.5, "Trigger y/n"}}}},

      // basic selection variables
      {"hV0Radius", "hV0Radius", {HistType::kTH1F, {{1000, 0.0f, 100.0f, "cm"}}}},
      {"hCascRadius", "hCascRadius", {HistType::kTH1F, {{1000, 0.0f, 100.0f, "cm"}}}},
      {"hV0CosPA", "hV0CosPA", {HistType::kTH1F, {{100, 0.95f, 1.0f}}}},
      {"hCascCosPA", "hCascCosPA", {HistType::kTH1F, {{100, 0.95f, 1.0f}}}},
      {"hDCAPosToPV", "hDCAPosToPV", {HistType::kTH1F, {vertexAxis}}},
      {"hDCANegToPV", "hDCANegToPV", {HistType::kTH1F, {vertexAxis}}},
      {"hDCABachToPV", "hDCABachToPV", {HistType::kTH1F, {vertexAxis}}},
      {"hDCAV0ToPV", "hDCAV0ToPV", {HistType::kTH1F, {vertexAxis}}},
      {"hDCAV0Dau", "hDCAV0Dau", {HistType::kTH1F, {{100, 0.0f, 10.0f, "cm^{2}"}}}},
      {"hDCACascDau", "hDCACascDau", {HistType::kTH1F, {{100, 0.0f, 10.0f, "cm^{2}"}}}},
      {"hLambdaMass", "hLambdaMass", {HistType::kTH1F, {{500, 1.0f, 1.5f, "Inv. Mass (GeV/c^{2})"}}}},

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

      {"hXiXiOS", "hXiXiOS", {HistType::kTHnSparseF, {deltaPhiAxis, deltaYAxis, ptAxis, ptAxis, invMassAxis, invMassAxis, vertexAxis, multiplicityAxis}}, true},
      {"hXiXiSS", "hXiXiSS", {HistType::kTHnSparseF, {deltaPhiAxis, deltaYAxis, ptAxis, ptAxis, invMassAxis, invMassAxis, vertexAxis, multiplicityAxis}}, true},
      {"hXiOmOS", "hXiOmOS", {HistType::kTHnSparseF, {deltaPhiAxis, deltaYAxis, ptAxis, ptAxis, invMassAxis, invMassAxis, vertexAxis, multiplicityAxis}}, true},
      {"hXiOmSS", "hXiOmSS", {HistType::kTHnSparseF, {deltaPhiAxis, deltaYAxis, ptAxis, ptAxis, invMassAxis, invMassAxis, vertexAxis, multiplicityAxis}}, true},
      {"hOmXiOS", "hOmXiOS", {HistType::kTHnSparseF, {deltaPhiAxis, deltaYAxis, ptAxis, ptAxis, invMassAxis, invMassAxis, vertexAxis, multiplicityAxis}}, true},
      {"hOmXiSS", "hOmXiSS", {HistType::kTHnSparseF, {deltaPhiAxis, deltaYAxis, ptAxis, ptAxis, invMassAxis, invMassAxis, vertexAxis, multiplicityAxis}}, true},
      {"hOmOmOS", "hOmOmOS", {HistType::kTHnSparseF, {deltaPhiAxis, deltaYAxis, ptAxis, ptAxis, invMassAxis, invMassAxis, vertexAxis, multiplicityAxis}}, true},
      {"hOmOmSS", "hOmOmSS", {HistType::kTHnSparseF, {deltaPhiAxis, deltaYAxis, ptAxis, ptAxis, invMassAxis, invMassAxis, vertexAxis, multiplicityAxis}}, true},

      // Mixed events
      {"MixedEvents/hMEVz1", "hMEVz1", {HistType::kTH1F, {vertexAxis}}},
      {"MixedEvents/hMEVz2", "hMEVz2", {HistType::kTH1F, {vertexAxis}}},
      {"MixedEvents/hMEDeltaPhiSS", "hMEDeltaPhiSS", {HistType::kTH1F, {deltaPhiAxis}}},
      {"MixedEvents/hMEDeltaPhiOS", "hMEDeltaPhiOS", {HistType::kTH1F, {deltaPhiAxis}}},
      {"MixedEvents/hMEQA", "hMEQA", {HistType::kTH1I, {{2, 0, 2, "QA for exceptions in ME (this histogram should have 0 entries!)"}}}},
      {"MixedEvents/hMEAutoCorrelation", "hMEAutoCorrelation", {HistType::kTH1I, {{4, -0.5f, 3.5f, "Types of SS autocorrelation"}}}},
      {"MixedEvents/hMEAutoCorrelationOS", "hMEAutoCorrelationOS", {HistType::kTH1I, {{2, -1.f, 1.f, "Charge of OS autocorrelated track"}}}},

      {"MixedEvents/hMEXiXiOS", "hMEXiXiOS", {HistType::kTHnSparseF, {deltaPhiAxis, deltaYAxis, ptAxis, ptAxis, invMassAxis, invMassAxis, vertexAxis, multiplicityAxis}}, true},
      {"MixedEvents/hMEXiXiSS", "hMEXiXiSS", {HistType::kTHnSparseF, {deltaPhiAxis, deltaYAxis, ptAxis, ptAxis, invMassAxis, invMassAxis, vertexAxis, multiplicityAxis}}, true},
      {"MixedEvents/hMEXiOmOS", "hMEXiOmOS", {HistType::kTHnSparseF, {deltaPhiAxis, deltaYAxis, ptAxis, ptAxis, invMassAxis, invMassAxis, vertexAxis, multiplicityAxis}}, true},
      {"MixedEvents/hMEXiOmSS", "hMEXiOmSS", {HistType::kTHnSparseF, {deltaPhiAxis, deltaYAxis, ptAxis, ptAxis, invMassAxis, invMassAxis, vertexAxis, multiplicityAxis}}, true},
      {"MixedEvents/hMEOmXiOS", "hMEOmXiOS", {HistType::kTHnSparseF, {deltaPhiAxis, deltaYAxis, ptAxis, ptAxis, invMassAxis, invMassAxis, vertexAxis, multiplicityAxis}}, true},
      {"MixedEvents/hMEOmXiSS", "hMEOmXiSS", {HistType::kTHnSparseF, {deltaPhiAxis, deltaYAxis, ptAxis, ptAxis, invMassAxis, invMassAxis, vertexAxis, multiplicityAxis}}, true},
      {"MixedEvents/hMEOmOmOS", "hMEOmOmOS", {HistType::kTHnSparseF, {deltaPhiAxis, deltaYAxis, ptAxis, ptAxis, invMassAxis, invMassAxis, vertexAxis, multiplicityAxis}}, true},
      {"MixedEvents/hMEOmOmSS", "hMEOmOmSS", {HistType::kTHnSparseF, {deltaPhiAxis, deltaYAxis, ptAxis, ptAxis, invMassAxis, invMassAxis, vertexAxis, multiplicityAxis}}, true},
    },
  };

  // cascade filter
  Filter cascadeSelector = aod::cascadeflags::isSelected > 0;

  // Warning: it is not possible to use this axis as configurable due to a bug - however, default values are sensible.
  SliceCache cache;
  ConfigurableAxis axisVtxZ{"axisVtxZ", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
  // ConfigurableAxis axisMult{"axisMult", {VARIABLE_WIDTH, 0, 5, 10, 20, 30, 40, 50, 100, 1000}, "Mixing bins - multiplicity"};
  // using BinningType = ColumnBinningPolicy<aod::collision::PosZ, aod::mult::MultFT0M<aod::mult::MultFT0A, aod::mult::MultFT0C>>;
  // BinningType colBinning{{axisVtxZ, axisMult}, true}; // true is for 'ignore overflows' (true by default). Underflows and overflows will have bin -1.
  using BinningType = ColumnBinningPolicy<aod::collision::PosZ>;
  BinningType colBinning{{axisVtxZ}, true}; // true is for 'ignore overflows' (true by default). Underflows and overflows will have bin -1.
  SameKindPair<MyCollisionsMult, MyCascades, BinningType> pair{colBinning, nMixedEvents, -1, &cache};

  void processSameEvent(MyCollisionsMult::iterator const& collision, MyCascades const& Cascades, aod::V0sLinked const&, aod::V0Datas const&, FullTracksExtIU const&, aod::BCsWithTimestamps const&)
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
          registry.fill(HIST("hMassXiMinus"), casc.mXi(), casc.pt());
          weight = 1. / getEfficiency(hEffXiMin, casc.pt());
        } else {
          registry.fill(HIST("hMassXiPlus"), casc.mXi(), casc.pt());
          weight = 1. / getEfficiency(hEffXiPlus, casc.pt());
        }
        // LOGF(info, "casc pt %f, weight %f", casc.pt(), weight);
        registry.fill(HIST("hMassXiEffCorrected"), casc.mXi(), casc.pt(), casc.yXi(), collision.posZ(), collision.multFT0M(), weight);
        registry.fill(HIST("hRapidityXi"), casc.yXi());
      }
      if (casc.isSelected() >= 2) { // consistent with Omega or both
        if (casc.sign() < 0) {
          registry.fill(HIST("hMassOmegaMinus"), casc.mOmega(), casc.pt());
          weight = 1. / getEfficiency(hEffOmegaMin, casc.pt());
        } else {
          registry.fill(HIST("hMassOmegaPlus"), casc.mOmega(), casc.pt());
          weight = 1. / getEfficiency(hEffOmegaPlus, casc.pt());
        }
        registry.fill(HIST("hMassOmegaEffCorrected"), casc.mOmega(), casc.pt(), casc.yOmega(), collision.posZ(), collision.multFT0M(), weight);
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

      registry.fill(HIST("hSelectionFlag"), casc.isSelected());
      registry.fill(HIST("hPhi"), casc.phi());
      registry.fill(HIST("hEta"), casc.eta());
    } // casc loop

    for (auto& [c0, c1] : combinations(Cascades, Cascades)) { // combinations automatically applies strictly upper in case of 2 identical tables
      // Define the trigger as the particle with the highest pT. As we can't swap the cascade tables themselves, we swap the addresses and later dereference them
      auto* triggerAddress = &c0;
      auto* assocAddress = &c1;
      if (assocAddress->pt() > triggerAddress->pt()) {
        std::swap(triggerAddress, assocAddress);
      }
      auto trigger = *triggerAddress;
      auto assoc = *assocAddress;

      // track indices for posterior checks
      // retains logic of V0 index while being safe wrt data model
      int posIdTrigg = trigger.posTrackId();
      int negIdTrigg = trigger.negTrackId();
      int posIdAssoc = assoc.posTrackId();
      int negIdAssoc = assoc.negTrackId();

      // calculate angular correlations
      double dphi = RecoDecay::constrainAngle(trigger.phi() - assoc.phi(), -PIHalf);

      double invMassXiTrigg = trigger.mXi();
      double invMassOmTrigg = trigger.mOmega();
      double invMassXiAssoc = assoc.mXi();
      double invMassOmAssoc = assoc.mOmega();

      double weightTrigg = 1.;
      double weightAssoc = 1.;

      // split into opposite-sign or same-sign
      if (trigger.sign() * assoc.sign() < 0) { // opposite-sign
        // check for autocorrelations between mis-identified kaons (omega bach) and protons (lambda daughter) TODO: improve logic?
        if (trigger.isSelected() >= 2) {
          if (trigger.sign() > 0 && trigger.bachelorId() == posIdAssoc) {
            // K+ from trigger Omega is the same as proton from assoc lambda
            registry.fill(HIST("hAutoCorrelationOS"), 1);
            continue;
          }
          if (trigger.sign() < 0 && trigger.bachelorId() == negIdAssoc) {
            // K- from trigger Omega is the same as antiproton from assoc antilambda
            registry.fill(HIST("hAutoCorrelationOS"), -1);
            continue;
          }
        }
        if (assoc.isSelected() >= 2) {
          if (assoc.sign() > 0 && assoc.bachelorId() == posIdTrigg) {
            // K+ from assoc Omega is the same as proton from trigger lambda
            registry.fill(HIST("hAutoCorrelationOS"), 1);
            continue;
          }
          if (assoc.sign() < 0 && assoc.bachelorId() == negIdTrigg) {
            // K- from assoc Omega is the same as antiproton from trigger antilambda
            registry.fill(HIST("hAutoCorrelationOS"), -1);
            continue;
          }
        }
        registry.fill(HIST("hDeltaPhiOS"), dphi);
        // Fill the different THnSparses depending on PID logic (important for rapidity & inv mass information)
        if (trigger.isSelected() <= 2 && TMath::Abs(trigger.yXi()) < maxRapidity) { // trigger Xi
          if (doEfficiencyCorrection)
            weightTrigg = trigger.sign() < 0 ? 1. / getEfficiency(hEffXiMin, trigger.pt()) : 1. / getEfficiency(hEffXiPlus, trigger.pt());
          if (assoc.isSelected() <= 2 && TMath::Abs(assoc.yXi()) < maxRapidity) { // assoc Xi
            if (doEfficiencyCorrection)
              weightAssoc = assoc.sign() < 0 ? 1. / getEfficiency(hEffXiMin, assoc.pt()) : 1. / getEfficiency(hEffXiPlus, assoc.pt());
            registry.fill(HIST("hXiXiOS"), dphi, trigger.yXi() - assoc.yXi(), trigger.pt(), assoc.pt(), invMassXiTrigg, invMassXiAssoc, collision.posZ(), collision.multFT0M(), weightTrigg * weightAssoc);
          }
          if (assoc.isSelected() >= 2 && TMath::Abs(assoc.yOmega()) < maxRapidity) { // assoc Omega
            if (doEfficiencyCorrection)
              weightAssoc = assoc.sign() < 0 ? 1. / getEfficiency(hEffOmegaMin, assoc.pt()) : 1. / getEfficiency(hEffOmegaPlus, assoc.pt());
            registry.fill(HIST("hXiOmOS"), dphi, trigger.yXi() - assoc.yOmega(), trigger.pt(), assoc.pt(), invMassXiTrigg, invMassOmAssoc, collision.posZ(), collision.multFT0M(), weightTrigg * weightAssoc);
          }
        }
        if (trigger.isSelected() >= 2 && TMath::Abs(trigger.yOmega()) < maxRapidity) { // trigger Omega
          if (doEfficiencyCorrection)
            weightTrigg = trigger.sign() < 0 ? 1. / getEfficiency(hEffOmegaMin, trigger.pt()) : 1. / getEfficiency(hEffOmegaPlus, trigger.pt());
          if (assoc.isSelected() <= 2 && TMath::Abs(assoc.yXi()) < maxRapidity) { // assoc Xi
            if (doEfficiencyCorrection)
              weightAssoc = assoc.sign() < 0 ? 1. / getEfficiency(hEffXiMin, assoc.pt()) : 1. / getEfficiency(hEffXiPlus, assoc.pt());
            registry.fill(HIST("hOmXiOS"), dphi, trigger.yOmega() - assoc.yXi(), trigger.pt(), assoc.pt(), invMassOmTrigg, invMassXiAssoc, collision.posZ(), collision.multFT0M(), weightTrigg * weightAssoc);
          }
          if (assoc.isSelected() >= 2 && TMath::Abs(assoc.yOmega()) < maxRapidity) { // assoc Omega
            if (doEfficiencyCorrection)
              weightAssoc = assoc.sign() < 0 ? 1. / getEfficiency(hEffOmegaMin, assoc.pt()) : 1. / getEfficiency(hEffOmegaPlus, assoc.pt());
            registry.fill(HIST("hOmOmOS"), dphi, trigger.yOmega() - assoc.yOmega(), trigger.pt(), assoc.pt(), invMassOmTrigg, invMassOmAssoc, collision.posZ(), collision.multFT0M(), weightTrigg * weightAssoc);
          }
        }
      } else { // same-sign
        // make sure to check for autocorrelations - only possible in same-sign correlations (if PID is correct)
        if (posIdTrigg == posIdAssoc && negIdTrigg == negIdAssoc) {
          // LOGF(info, "same v0 in SS correlation! %d %d", v0dataTrigg.v0Id(), v0dataAssoc.v0Id());
          registry.fill(HIST("hAutoCorrelation"), 0);
          continue;
        }
        int bachIdTrigg = trigger.bachelorId();
        int bachIdAssoc = assoc.bachelorId();

        if (bachIdTrigg == bachIdAssoc) {
          // LOGF(info, "same bachelor in SS correlation! %d %d", bachIdTrigg, bachIdAssoc);
          registry.fill(HIST("hAutoCorrelation"), 1);
          continue;
        }
        // check for same tracks in v0's of cascades
        if (negIdTrigg == negIdAssoc || posIdTrigg == posIdAssoc) {
          // LOGF(info, "cascades have a v0-track in common in SS correlation!");
          registry.fill(HIST("hAutoCorrelation"), 2);
          continue;
        }
        if (trigger.sign() < 0) { // neg cascade
          if (negIdTrigg == bachIdAssoc || negIdAssoc == bachIdTrigg) {
            // LOGF(info, "bach of casc == v0-pion of other casc in neg SS correlation!");
            registry.fill(HIST("hAutoCorrelation"), 3);
            continue;
          }
        } else { // pos cascade
          if (posIdTrigg == bachIdAssoc || posIdAssoc == bachIdTrigg) {
            // LOGF(info, "bach of casc == v0-pion of other casc in pos SS correlation!");
            registry.fill(HIST("hAutoCorrelation"), 3);
            continue;
          }
        }
        registry.fill(HIST("hDeltaPhiSS"), dphi);
        // Fill the different THnSparses depending on PID logic (important for rapidity & inv mass information)
        if (trigger.isSelected() <= 2 && TMath::Abs(trigger.yXi()) < maxRapidity) { // trigger Xi
          if (doEfficiencyCorrection)
            weightTrigg = trigger.sign() < 0 ? 1. / getEfficiency(hEffXiMin, trigger.pt()) : 1. / getEfficiency(hEffXiPlus, trigger.pt());
          if (assoc.isSelected() <= 2 && TMath::Abs(assoc.yXi()) < maxRapidity) { // assoc Xi
            if (doEfficiencyCorrection)
              weightAssoc = assoc.sign() < 0 ? 1. / getEfficiency(hEffXiMin, assoc.pt()) : 1. / getEfficiency(hEffXiPlus, assoc.pt());
            registry.fill(HIST("hXiXiSS"), dphi, trigger.yXi() - assoc.yXi(), trigger.pt(), assoc.pt(), invMassXiTrigg, invMassXiAssoc, collision.posZ(), collision.multFT0M(), weightTrigg * weightAssoc);
          }
          if (assoc.isSelected() >= 2 && TMath::Abs(assoc.yOmega()) < maxRapidity) { // assoc Omega
            if (doEfficiencyCorrection)
              weightAssoc = assoc.sign() < 0 ? 1. / getEfficiency(hEffOmegaMin, assoc.pt()) : 1. / getEfficiency(hEffOmegaPlus, assoc.pt());
            registry.fill(HIST("hXiOmSS"), dphi, trigger.yXi() - assoc.yOmega(), trigger.pt(), assoc.pt(), invMassXiTrigg, invMassOmAssoc, collision.posZ(), collision.multFT0M(), weightTrigg * weightAssoc);
          }
        }
        if (trigger.isSelected() >= 2 && TMath::Abs(trigger.yOmega()) < maxRapidity) { // trigger Omega
          if (doEfficiencyCorrection)
            weightTrigg = trigger.sign() < 0 ? 1. / getEfficiency(hEffOmegaMin, trigger.pt()) : 1. / getEfficiency(hEffOmegaPlus, trigger.pt());
          if (assoc.isSelected() <= 2 && TMath::Abs(assoc.yXi()) < maxRapidity) { // assoc Xi
            if (doEfficiencyCorrection)
              weightAssoc = assoc.sign() < 0 ? 1. / getEfficiency(hEffXiMin, assoc.pt()) : 1. / getEfficiency(hEffXiPlus, assoc.pt());
            registry.fill(HIST("hOmXiSS"), dphi, trigger.yOmega() - assoc.yXi(), trigger.pt(), assoc.pt(), invMassOmTrigg, invMassXiAssoc, collision.posZ(), collision.multFT0M(), weightTrigg * weightAssoc);
          }
          if (assoc.isSelected() >= 2 && TMath::Abs(assoc.yOmega()) < maxRapidity) { // assoc Omega
            if (doEfficiencyCorrection)
              weightAssoc = assoc.sign() < 0 ? 1. / getEfficiency(hEffOmegaMin, assoc.pt()) : 1. / getEfficiency(hEffOmegaPlus, assoc.pt());
            registry.fill(HIST("hOmOmSS"), dphi, trigger.yOmega() - assoc.yOmega(), trigger.pt(), assoc.pt(), invMassOmTrigg, invMassOmAssoc, collision.posZ(), collision.multFT0M(), weightTrigg * weightAssoc);
          }
        }
      }
    } // correlations
  }   // process same event

  void processMixedEvent(MyCollisionsMult const& /*collisions*/, MyCascades const& /*Cascades*/,
                         aod::V0sLinked const&, aod::V0Datas const&, FullTracksExtIU const&)
  {
    for (auto const& [col1, cascades1, col2, cascades2] : pair) {
      if (!col1.sel8() || !col2.sel8())
        continue;
      if (TMath::Abs(col1.posZ()) > zVertexCut || TMath::Abs(col2.posZ()) > zVertexCut)
        continue;
      if (col1.globalIndex() == col2.globalIndex()) {
        registry.fill(HIST("hMEQA"), 0.5);
        continue;
      }

      registry.fill(HIST("MixedEvents/hMEVz1"), col1.posZ());
      registry.fill(HIST("MixedEvents/hMEVz2"), col2.posZ());

      for (auto& [casc1, casc2] : combinations(CombinationsFullIndexPolicy(cascades1, cascades2))) {
        // specify FullIndexPolicy since the cascades are from different collisions
        auto* triggerAddress = &casc1;
        auto* assocAddress = &casc2;
        if (assocAddress->pt() > triggerAddress->pt()) {
          std::swap(triggerAddress, assocAddress);
        }
        auto trigger = *triggerAddress;
        auto assoc = *assocAddress;

        if (trigger.collisionId() == assoc.collisionId()) {
          registry.fill(HIST("hMEQA"), 1.5);
          continue;
        }

        double dphi = RecoDecay::constrainAngle(trigger.phi() - assoc.phi(), -PIHalf);

        double invMassXiTrigg = trigger.mXi();
        double invMassOmTrigg = trigger.mOmega();
        double invMassXiAssoc = assoc.mXi();
        double invMassOmAssoc = assoc.mOmega();

        // V0 daughter track ID's used for autocorrelation check
        int posIdTrigg = trigger.posTrackId();
        int negIdTrigg = trigger.negTrackId();
        int posIdAssoc = assoc.posTrackId();
        int negIdAssoc = assoc.negTrackId();

        double weightTrigg = 1.;
        double weightAssoc = 1.;

        if (trigger.sign() * assoc.sign() < 0) { // opposite-sign

          // check for autocorrelations between mis-identified kaons (omega bach) and protons (lambda daughter) TODO: improve logic?
          if (trigger.isSelected() >= 2) {
            if (trigger.sign() > 0 && trigger.bachelorId() == posIdAssoc) {
              // K+ from trigger Omega is the same as proton from assoc lambda
              registry.fill(HIST("MixedEvents/hMEAutoCorrelationOS"), 1);
              continue;
            }
            if (trigger.sign() < 0 && trigger.bachelorId() == negIdAssoc) {
              // K- from trigger Omega is the same as antiproton from assoc antilambda
              registry.fill(HIST("MixedEvents/hMEAutoCorrelationOS"), -1);
              continue;
            }
          }
          if (assoc.isSelected() >= 2) {
            if (assoc.sign() > 0 && assoc.bachelorId() == posIdTrigg) {
              // K+ from assoc Omega is the same as proton from trigger lambda
              registry.fill(HIST("MixedEvents/hMEAutoCorrelationOS"), 1);
              continue;
            }
            if (assoc.sign() < 0 && assoc.bachelorId() == negIdTrigg) {
              // K- from assoc Omega is the same as antiproton from trigger antilambda
              registry.fill(HIST("MixedEvents/hMEAutoCorrelationOS"), -1);
              continue;
            }
          }

          registry.fill(HIST("MixedEvents/hMEDeltaPhiOS"), dphi);

          // Fill the different THnSparses depending on PID logic (important for rapidity & inv mass information)
          if (trigger.isSelected() <= 2 && TMath::Abs(trigger.yXi()) < maxRapidity) { // trigger Xi
            if (doEfficiencyCorrection)
              weightTrigg = trigger.sign() < 0 ? 1. / getEfficiency(hEffXiMin, trigger.pt()) : 1. / getEfficiency(hEffXiPlus, trigger.pt());
            if (assoc.isSelected() <= 2 && TMath::Abs(assoc.yXi()) < maxRapidity) { // assoc Xi
              if (doEfficiencyCorrection)
                weightAssoc = assoc.sign() < 0 ? 1. / getEfficiency(hEffXiMin, assoc.pt()) : 1. / getEfficiency(hEffXiPlus, assoc.pt());
              registry.fill(HIST("MixedEvents/hMEXiXiOS"), dphi, trigger.yXi() - assoc.yXi(), trigger.pt(), assoc.pt(), invMassXiTrigg, invMassXiAssoc, col1.posZ(), col1.multFT0M(), weightTrigg * weightAssoc);
            }
            if (assoc.isSelected() >= 2 && TMath::Abs(assoc.yOmega()) < maxRapidity) { // assoc Omega
              if (doEfficiencyCorrection)
                weightAssoc = assoc.sign() < 0 ? 1. / getEfficiency(hEffOmegaMin, assoc.pt()) : 1. / getEfficiency(hEffOmegaPlus, assoc.pt());
              registry.fill(HIST("MixedEvents/hMEXiOmOS"), dphi, trigger.yXi() - assoc.yOmega(), trigger.pt(), assoc.pt(), invMassXiTrigg, invMassOmAssoc, col1.posZ(), col1.multFT0M(), weightTrigg * weightAssoc);
            }
          }
          if (trigger.isSelected() >= 2 && TMath::Abs(trigger.yOmega()) < maxRapidity) { // trigger Omega
            if (doEfficiencyCorrection)
              weightTrigg = trigger.sign() < 0 ? 1. / getEfficiency(hEffOmegaMin, trigger.pt()) : 1. / getEfficiency(hEffOmegaPlus, trigger.pt());
            if (assoc.isSelected() <= 2 && TMath::Abs(assoc.yXi()) < maxRapidity) { // assoc Xi
              if (doEfficiencyCorrection)
                weightAssoc = assoc.sign() < 0 ? 1. / getEfficiency(hEffXiMin, assoc.pt()) : 1. / getEfficiency(hEffXiPlus, assoc.pt());
              registry.fill(HIST("MixedEvents/hMEOmXiOS"), dphi, trigger.yOmega() - assoc.yXi(), trigger.pt(), assoc.pt(), invMassOmTrigg, invMassXiAssoc, col1.posZ(), col1.multFT0M(), weightTrigg * weightAssoc);
            }
            if (assoc.isSelected() >= 2 && TMath::Abs(assoc.yOmega()) < maxRapidity) { // assoc Omega
              if (doEfficiencyCorrection)
                weightAssoc = assoc.sign() < 0 ? 1. / getEfficiency(hEffOmegaMin, assoc.pt()) : 1. / getEfficiency(hEffOmegaPlus, assoc.pt());
              registry.fill(HIST("MixedEvents/hMEOmOmOS"), dphi, trigger.yOmega() - assoc.yOmega(), trigger.pt(), assoc.pt(), invMassOmTrigg, invMassOmAssoc, col1.posZ(), col1.multFT0M(), weightTrigg * weightAssoc);
            }
          }
        } else { // same sign
          // make sure to check for autocorrelations - only possible in same-sign correlations (if PID is correct)
          if (posIdTrigg == posIdAssoc && negIdTrigg == negIdAssoc) {
            // LOGF(info, "same v0 in SS correlation! %d %d", v0dataTrigg.v0Id(), v0dataAssoc.v0Id());
            registry.fill(HIST("MixedEvents/hMEAutoCorrelation"), 0);
            continue;
          }
          int bachIdTrigg = trigger.bachelorId();
          int bachIdAssoc = assoc.bachelorId();

          if (bachIdTrigg == bachIdAssoc) {
            // LOGF(info, "same bachelor in SS correlation! %d %d", bachIdTrigg, bachIdAssoc);
            registry.fill(HIST("MixedEvents/hMEAutoCorrelation"), 1);
            continue;
          }
          // check for same tracks in v0's of cascades
          if (negIdTrigg == negIdAssoc || posIdTrigg == posIdAssoc) {
            // LOGF(info, "cascades have a v0-track in common in SS correlation!");
            registry.fill(HIST("MixedEvents/hMEAutoCorrelation"), 2);
            continue;
          }
          if (trigger.sign() < 0) { // neg cascade
            if (negIdTrigg == bachIdAssoc || negIdAssoc == bachIdTrigg) {
              // LOGF(info, "bach of casc == v0-pion of other casc in neg SS correlation!");
              registry.fill(HIST("MixedEvents/hMEAutoCorrelation"), 3);
              continue;
            }
          } else { // pos cascade
            if (posIdTrigg == bachIdAssoc || posIdAssoc == bachIdTrigg) {
              // LOGF(info, "bach of casc == v0-pion of other casc in pos SS correlation!");
              registry.fill(HIST("MixedEvents/hMEAutoCorrelation"), 3);
              continue;
            }
          }

          registry.fill(HIST("MixedEvents/hMEDeltaPhiSS"), dphi);

          if (trigger.isSelected() <= 2 && TMath::Abs(trigger.yXi()) < maxRapidity) { // trigger Xi
            if (doEfficiencyCorrection)
              weightTrigg = trigger.sign() < 0 ? 1. / getEfficiency(hEffXiMin, trigger.pt()) : 1. / getEfficiency(hEffXiPlus, trigger.pt());
            if (assoc.isSelected() <= 2 && TMath::Abs(assoc.yXi()) < maxRapidity) { // assoc Xi
              if (doEfficiencyCorrection)
                weightAssoc = assoc.sign() < 0 ? 1. / getEfficiency(hEffXiMin, assoc.pt()) : 1. / getEfficiency(hEffXiPlus, assoc.pt());
              registry.fill(HIST("MixedEvents/hMEXiXiSS"), dphi, trigger.yXi() - assoc.yXi(), trigger.pt(), assoc.pt(), invMassXiTrigg, invMassXiAssoc, col1.posZ(), col1.multFT0M(), weightTrigg * weightAssoc);
            }
            if (assoc.isSelected() >= 2 && TMath::Abs(assoc.yOmega()) < maxRapidity) { // assoc Omega
              if (doEfficiencyCorrection)
                weightAssoc = assoc.sign() < 0 ? 1. / getEfficiency(hEffOmegaMin, assoc.pt()) : 1. / getEfficiency(hEffOmegaPlus, assoc.pt());
              registry.fill(HIST("MixedEvents/hMEXiOmSS"), dphi, trigger.yXi() - assoc.yOmega(), trigger.pt(), assoc.pt(), invMassXiTrigg, invMassOmAssoc, col1.posZ(), col1.multFT0M(), weightTrigg * weightAssoc);
            }
          }
          if (trigger.isSelected() >= 2 && TMath::Abs(trigger.yOmega()) < maxRapidity) { // trigger Omega
            if (doEfficiencyCorrection)
              weightTrigg = trigger.sign() < 0 ? 1. / getEfficiency(hEffOmegaMin, trigger.pt()) : 1. / getEfficiency(hEffOmegaPlus, trigger.pt());
            if (assoc.isSelected() <= 2 && TMath::Abs(assoc.yXi()) < maxRapidity) { // assoc Xi
              if (doEfficiencyCorrection)
                weightAssoc = assoc.sign() < 0 ? 1. / getEfficiency(hEffXiMin, assoc.pt()) : 1. / getEfficiency(hEffXiPlus, assoc.pt());
              registry.fill(HIST("MixedEvents/hMEOmXiSS"), dphi, trigger.yOmega() - assoc.yXi(), trigger.pt(), assoc.pt(), invMassOmTrigg, invMassXiAssoc, col1.posZ(), col1.multFT0M(), weightTrigg * weightAssoc);
            }
            if (assoc.isSelected() >= 2 && TMath::Abs(assoc.yOmega()) < maxRapidity) { // assoc Omega
              if (doEfficiencyCorrection)
                weightAssoc = assoc.sign() < 0 ? 1. / getEfficiency(hEffOmegaMin, assoc.pt()) : 1. / getEfficiency(hEffOmegaPlus, assoc.pt());
              registry.fill(HIST("MixedEvents/hMEOmOmSS"), dphi, trigger.yOmega() - assoc.yOmega(), trigger.pt(), assoc.pt(), invMassOmTrigg, invMassOmAssoc, col1.posZ(), col1.multFT0M(), weightTrigg * weightAssoc);
            }
          }
        } // same sign
      } // correlations
    } // collisions
  } // process mixed events

  PROCESS_SWITCH(CascadeCorrelations, processSameEvent, "Process same events", true);
  PROCESS_SWITCH(CascadeCorrelations, processMixedEvent, "Process mixed events", true);

}; // struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<CascadeSelector>(cfgc),
    adaptAnalysisTask<CascadeCorrelations>(cfgc)};
}
