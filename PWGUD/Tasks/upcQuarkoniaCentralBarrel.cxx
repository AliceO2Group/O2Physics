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
/// \file upcQuarkoniaCentralBarrel.cxx
/// \brief quarkonia --> ppbar task
///
/// \author David Dobrigkeit Chinellato <david.dobrigkeit.chinellato@cern.ch>, Austrian Academy of Sciences & SMI
/// \author Roman Lavicka <roman.lavicka@cern.ch>, Austrian Academy of Sciences & SMI
/// \author Romain Schotter <romain.schotter@cern.ch>, Austrian Academy of Sciences & SMI
//
// V0 analysis task
// ================
//
// This code loops over a V0Cores table and produces some
// standard analysis output. It is meant to be run over
// derived data.
//
//
//    Comments, questions, complaints, suggestions?
//    Please write to:
//    romain.schotter@cern.ch
//    david.dobrigkeit.chinellato@cern.ch
//

#include "PWGUD/Core/SGSelector.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "EventFiltering/Zorro.h"
#include "EventFiltering/ZorroSummary.h"

#include "CCDB/BasicCCDBManager.h"
#include "CommonConstants/PhysicsConstants.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

#include <Math/Vector4D.h>
#include <TFile.h>
#include <TH2F.h>
#include <TLorentzVector.h>
#include <TPDGCode.h>
#include <TProfile.h>

#include <array>
#include <cmath>
#include <cstdlib>
#include <map>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;

using UDCollisions = soa::Join<aod::UDCollisions, aod::UDCollisionsSels, aod::UDZdcsReduced>;
using UDCollision = UDCollisions::iterator;
using UDTracks = soa::Join<aod::UDTracks, aod::UDTracksPID, aod::UDTracksExtra, aod::UDTracksFlags, aod::UDTracksDCA>;
using UDTrack = UDTracks::iterator;

// simple checkers, but ensure 64 bit integers
#define BITSET(var, nbit) ((var) |= (static_cast<uint64_t>(1) << static_cast<uint64_t>(nbit)))
#define BITCHECK(var, nbit) ((var) & (static_cast<uint64_t>(1) << static_cast<uint64_t>(nbit)))

struct upcQuarkoniaCentralBarrel {
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<bool> buildSameSignPairs{"buildSameSignPairs", false, "If true: build same-sign pairs, otherwise consider only opposite-sign pairs"};

  // rapidity cut on the hyperon-antiHyperon pair
  Configurable<float> rapidityCut{"rapidityCut", 0.5, "rapidity cut on the ppbar pair"};

  struct : ConfigurableGroup {
    // Selection criteria: acceptance
    Configurable<float> etaCut{"trackSelections.etaCut", 0.8, "max eta for daughters"};

    // Track quality
    Configurable<float> dcaxytopv{"trackSelections.dcaxytopv", .05, "max transverse DCA to PV (cm)"};
    Configurable<float> dcaztopv{"trackSelections.dcaztopv", .05, "max longitudinal DCA to PV (cm)"};
    Configurable<int> minTPCrows{"trackSelections.minTPCrows", 70, "minimum TPC crossed rows"};
    Configurable<int> minTPCclusters{"trackSelections.minTPCclusters", 3, "minimum TPC clusters"};
    Configurable<float> minTPCchi2clusters{"trackSelections.minTPCchi2clusters", 4.0, "minimum TPC chi2/clusters"};
    Configurable<float> minTPCrowsoverfindable{"trackSelections.minTPCrowsoverfindable", 0.8, "minimum TPC rows/findable clusters"};
    Configurable<int> minITSclusters{"trackSelections.minITSclusters", -1, "minimum ITS clusters"};
    Configurable<float> minITSchi2clusters{"trackSelections.minITSchi2clusters", -1.0, "minimum ITS chi2/clusters"};
    Configurable<bool> requirePVcontributor{"trackSelections.requirePVcontributor", false, "require that track is a PV contributor"};
    Configurable<bool> applyDCAptdepsel{"trackSelections.applyDCAptdepsel", false, "apply DCA pt dep. cut à la Run2"};
    Configurable<bool> skipTPConly{"trackSelections.skipTPConly", false, "skip V0s comprised of at least one TPC only prong"};
    Configurable<bool> requireITSonly{"trackSelections.requirePosITSonly", false, "require that track is ITSonly (overrides TPC quality)"};
    Configurable<bool> rejectITSafterburner{"trackSelections.rejectNegITSafterburner", false, "reject track formed out of afterburner ITS tracks"};

    // PID (TPC/TOF)
    Configurable<float> tpcPidNsigmaCut{"trackSelections.tpcPidNsigmaCut", 5, "tpcPidNsigmaCut"};
    Configurable<float> tofPidNsigmaCut{"trackSelections.tofPidNsigmaCut", 1e+6, "tofPidNsigmaCut"};
  } trackSelections;

  // for MC
  Configurable<bool> doMCAssociation{"doMCAssociation", true, "if MC, do MC association"};

  // UPC selections
  SGSelector sgSelector;
  struct : ConfigurableGroup {
    Configurable<float> fv0Cut{"upcCuts.fv0Cut", 100., "FV0A threshold"};
    Configurable<float> ft0aCut{"upcCuts.ft0aCut", 200., "FT0A threshold"};
    Configurable<float> ft0cCut{"upcCuts.ft0cCut", 100., "FT0C threshold"};
    Configurable<float> zdcCut{"upcCuts.zdcCut", 10., "ZDC threshold"};
    // Configurable<float> gapSel{"upcCuts.gapSel", 2, "Gap selection"};
  } upcCuts;

  static constexpr float defaultLifetimeCuts[1][2] = {{30., 20.}};
  Configurable<LabeledArray<float>> lifetimecut{"lifetimecut", {defaultLifetimeCuts[0], 2, {"lifetimecutLambda", "lifetimecutK0S"}}, "lifetimecut"};

  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.0f, 0.2f, 0.4f, 0.6f, 0.8f, 1.0f, 1.2f, 1.4f, 1.6f, 1.8f, 2.0f, 2.4f, 2.8f, 3.2f, 3.6f, 4.0f, 4.8f, 5.6f, 6.5f, 7.5f, 9.0f, 11.0f, 13.0f, 15.0f, 19.0f, 23.0f, 30.0f, 40.0f, 50.0f}, "pt axis for analysis"};
  ConfigurableAxis axisQuarkoniumMass{"axisQuarkoniumMass", {500, 2.600f, 4.000f}, "M (hyp. #bar{hyp.} ) (GeV/#it{c}^{2})"};
  ConfigurableAxis axisOccupancy{"axisOccupancy", {VARIABLE_WIDTH, 0.0f, 250.0f, 500.0f, 750.0f, 1000.0f, 1500.0f, 2000.0f, 3000.0f, 4500.0f, 6000.0f, 8000.0f, 10000.0f, 50000.0f}, "Occupancy"};

  ConfigurableAxis axisDCAXYtoPV{"axisDCAXYtoPV", {20, 0.0f, 1.0f}, "DCAxy (cm)"};
  ConfigurableAxis axisDCAZtoPV{"axisDCAZtoPV", {20, 0.0f, 1.0f}, "DCAz (cm)"};
  ConfigurableAxis axisTPCrows{"axisTPCrows", {160, 0.0f, 160.0f}, "N TPC rows"};
  ConfigurableAxis axisTPCclus{"axisTPCclus", {160, 0.0f, 160.0f}, "N TPC Clusters"};
  ConfigurableAxis axisTPCChi2clus{"axisTPCChi2clus", {100, 0.0f, 50.0f}, "TPC Chi2/Clusters"};
  ConfigurableAxis axisTPCrowsOverFindable{"axisTPCrowsOverFindable", {100, 0.0f, 1.0f}, "TPC Rows/Findable"};
  ConfigurableAxis axisITSclus{"axisITSclus", {7, 0.0f, 7.0f}, "N ITS Clusters"};
  ConfigurableAxis axisITSChi2clus{"axisITSChi2clus", {100, 0.0f, 50.0f}, "ITS Chi2/Clusters"};
  ConfigurableAxis axisNsigmaTPC{"axisNsigmaTPC", {200, -10.0f, 10.0f}, "N sigma TPC"};
  ConfigurableAxis axisNsigmaTOF{"axisNsigmaTOF", {200, -10.0f, 10.0f}, "N sigma TOF"};

  // UPC axes
  ConfigurableAxis axisSelGap{"axisSelGap", {4, -1.5, 2.5}, "Gap side"};

  // PDG database
  Service<o2::framework::O2DatabasePDG> pdgDB;

  void init(InitContext const&)
  {
    // Event Counters
    histos.add("hEventPVz", "hEventPVz", kTH1F, {{100, -20.0f, +20.0f}});
    histos.add("hSelGapSideVsPVz", "hSelGapSideVsPVz", kTH2F, {axisSelGap, {100, -20.0f, +20.0f}});

    histos.add("hEventOccupancy", "hEventOccupancy", kTH1F, {axisOccupancy});
    histos.add("hSelGapSideVsOccupancy", "hSelGapSideVsOccupancy", kTH2F, {axisSelGap, axisOccupancy});

    histos.add("hGapSide", "Gap side; Entries", kTH1F, {{5, -0.5, 4.5}});
    histos.add("hSelGapSide", "Selected gap side; Entries", kTH1F, {axisSelGap});

    // histograms versus mass
    histos.add("PPbar/h2dMassPPbar", "h2dMassPPbar", kTH2F, {axisPt, axisQuarkoniumMass});
    // Non-UPC info
    histos.add("PPbar/h2dMassPPbarHadronic", "h2dMassPPbarHadronic", kTH2F, {axisPt, axisQuarkoniumMass});
    // UPC info
    histos.add("PPbar/h2dMassPPbarSGA", "h2dMassPPbarSGA", kTH2F, {axisPt, axisQuarkoniumMass});
    histos.add("PPbar/h2dMassPPbarSGC", "h2dMassPPbarSGC", kTH2F, {axisPt, axisQuarkoniumMass});
    histos.add("PPbar/h2dMassPPbarDG", "h2dMassPPbarDG", kTH2F, {axisPt, axisQuarkoniumMass});

    histos.add("PPbar/h2dNbrOfProtonsVsSelGapSide", "h2dNbrOfProtonsVsSelGapSide", kTH2F, {axisSelGap, {100, -0.5f, 99.5f}});
    histos.add("PPbar/h2dNbrOfAntiProtonsVsSelGapSide", "h2dNbrOfAntiProtonsVsSelGapSide", kTH2F, {axisSelGap, {100, -0.5f, 99.5f}});
    // QA plot
    // Proton Candidates before selections
    histos.add("PPbar/Proton/hDCAxyToPV", "hDCAxyToPV", kTH1F, {axisDCAXYtoPV});
    histos.add("PPbar/Proton/hDCAzToPV", "hDCAzToPV", kTH1F, {axisDCAZtoPV});
    histos.add("PPbar/Proton/hTPCCrossedRows", "hTPCCrossedRows", kTH1F, {axisTPCrows});
    histos.add("PPbar/Proton/hTPCNClusters", "hTPCNClusters", kTH1F, {axisTPCclus});
    histos.add("PPbar/Proton/hTPCChi2Clusters", "hTPCChi2Clusters", kTH1F, {axisTPCChi2clus});
    histos.add("PPbar/Proton/hTPCCrossedRowsOverFindable", "hTPCCrossedRowsOverFindable", kTH1F, {axisTPCrowsOverFindable});
    histos.add("PPbar/Proton/hITSNClusters", "hITSNClusters", kTH1F, {axisITSclus});
    histos.add("PPbar/Proton/hITSChi2Clusters", "hITSChi2Clusters", kTH1F, {axisITSChi2clus});
    histos.add("PPbar/Proton/hTPCNsigma", "hTPCNsigma", kTH1F, {axisNsigmaTPC});
    histos.add("PPbar/Proton/hTOFNsigma", "hTOFNsigma", kTH1F, {axisNsigmaTOF});
    histos.add("PPbar/Proton/h2dITSvsTPCpts", "h2dITSvsTPCpts", kTH2F, {axisTPCrows, axisITSclus});
    // Anti Proton before selections
    histos.add("PPbar/AntiProton/hDCAxyToPV", "hDCAxyToPV", kTH1F, {axisDCAXYtoPV});
    histos.add("PPbar/AntiProton/hDCAzToPV", "hDCAzToPV", kTH1F, {axisDCAZtoPV});
    histos.add("PPbar/AntiProton/hTPCCrossedRows", "hTPCCrossedRows", kTH1F, {axisTPCrows});
    histos.add("PPbar/AntiProton/hTPCNClusters", "hTPCNClusters", kTH1F, {axisTPCclus});
    histos.add("PPbar/AntiProton/hTPCChi2Clusters", "hTPCChi2Clusters", kTH1F, {axisTPCChi2clus});
    histos.add("PPbar/AntiProton/hTPCCrossedRowsOverFindable", "hTPCCrossedRowsOverFindable", kTH1F, {axisTPCrowsOverFindable});
    histos.add("PPbar/AntiProton/hITSNClusters", "hITSNClusters", kTH1F, {axisITSclus});
    histos.add("PPbar/AntiProton/hITSChi2Clusters", "hITSChi2Clusters", kTH1F, {axisITSChi2clus});
    histos.add("PPbar/AntiProton/hTPCNsigma", "hTPCNsigma", kTH1F, {axisNsigmaTPC});
    histos.add("PPbar/AntiProton/hTOFNsigma", "hTOFNsigma", kTH1F, {axisNsigmaTOF});
    histos.add("PPbar/AntiProton/h2dITSvsTPCpts", "h2dITSvsTPCpts", kTH2F, {axisTPCrows, axisITSclus});
    // Proton Candidates after selections
    histos.add("PPbar/Proton/hDCAxyToPV_aftersel", "hDCAxyToPV", kTH1F, {axisDCAXYtoPV});
    histos.add("PPbar/Proton/hDCAzToPV_aftersel", "hDCAzToPV", kTH1F, {axisDCAZtoPV});
    histos.add("PPbar/Proton/hTPCCrossedRows_aftersel", "hTPCCrossedRows", kTH1F, {axisTPCrows});
    histos.add("PPbar/Proton/hTPCNClusters_aftersel", "hTPCNClusters", kTH1F, {axisTPCclus});
    histos.add("PPbar/Proton/hTPCChi2Clusters_aftersel", "hTPCChi2Clusters", kTH1F, {axisTPCChi2clus});
    histos.add("PPbar/Proton/hTPCCrossedRowsOverFindable_aftersel", "hTPCCrossedRowsOverFindable", kTH1F, {axisTPCrowsOverFindable});
    histos.add("PPbar/Proton/hITSNClusters_aftersel", "hITSNClusters", kTH1F, {axisITSclus});
    histos.add("PPbar/Proton/hITSChi2Clusters_aftersel", "hITSChi2Clusters", kTH1F, {axisITSChi2clus});
    histos.add("PPbar/Proton/hTPCNsigma_aftersel", "hTPCNsigma", kTH1F, {axisNsigmaTPC});
    histos.add("PPbar/Proton/hTOFNsigma_aftersel", "hTOFNsigma", kTH1F, {axisNsigmaTOF});
    histos.add("PPbar/Proton/h2dITSvsTPCpts_aftersel", "h2dITSvsTPCpts", kTH2F, {axisTPCrows, axisITSclus});
    // Anti Proton after selections
    histos.add("PPbar/AntiProton/hDCAxyToPV_aftersel", "hDCAxyToPV", kTH1F, {axisDCAXYtoPV});
    histos.add("PPbar/AntiProton/hDCAzToPV_aftersel", "hDCAzToPV", kTH1F, {axisDCAZtoPV});
    histos.add("PPbar/AntiProton/hTPCCrossedRows_aftersel", "hTPCCrossedRows", kTH1F, {axisTPCrows});
    histos.add("PPbar/AntiProton/hTPCNClusters_aftersel", "hTPCNClusters", kTH1F, {axisTPCclus});
    histos.add("PPbar/AntiProton/hTPCChi2Clusters_aftersel", "hTPCChi2Clusters", kTH1F, {axisTPCChi2clus});
    histos.add("PPbar/AntiProton/hTPCCrossedRowsOverFindable_aftersel", "hTPCCrossedRowsOverFindable", kTH1F, {axisTPCrowsOverFindable});
    histos.add("PPbar/AntiProton/hITSNClusters_aftersel", "hITSNClusters", kTH1F, {axisITSclus});
    histos.add("PPbar/AntiProton/hITSChi2Clusters_aftersel", "hITSChi2Clusters", kTH1F, {axisITSChi2clus});
    histos.add("PPbar/AntiProton/hTPCNsigma_aftersel", "hTPCNsigma", kTH1F, {axisNsigmaTPC});
    histos.add("PPbar/AntiProton/hTOFNsigma_aftersel", "hTOFNsigma", kTH1F, {axisNsigmaTOF});
    histos.add("PPbar/AntiProton/h2dITSvsTPCpts_aftersel", "h2dITSvsTPCpts", kTH2F, {axisTPCrows, axisITSclus});
    if (doMCAssociation) {
      histos.add("PPbar/h2dInvMassTrueEtaC1S", "h2dInvMassTrueEtaC1S", kTH2F, {axisPt, axisQuarkoniumMass});
      histos.add("PPbar/h2dInvMassTrueJPsi", "h2dInvMassTrueJPsi", kTH2F, {axisPt, axisQuarkoniumMass});
      histos.add("PPbar/h2dInvMassTrueChiC0", "h2dInvMassTrueChiC0", kTH2F, {axisPt, axisQuarkoniumMass});
      histos.add("PPbar/h2dInvMassTrueChiC1", "h2dInvMassTrueChiC1", kTH2F, {axisPt, axisQuarkoniumMass});
      histos.add("PPbar/h2dInvMassTrueHC", "h2dInvMassTrueHC", kTH2F, {axisPt, axisQuarkoniumMass});
      histos.add("PPbar/h2dInvMassTrueChiC2", "h2dInvMassTrueChiC2", kTH2F, {axisPt, axisQuarkoniumMass});
      histos.add("PPbar/h2dInvMassTrueEtaC2S", "h2dInvMassTrueEtaC2S", kTH2F, {axisPt, axisQuarkoniumMass});
      histos.add("PPbar/h2dInvMassTruePsi2S", "h2dInvMassTruePsi2S", kTH2F, {axisPt, axisQuarkoniumMass});
    }

    // inspect histogram sizes, please
    histos.print();
  }

  template <typename TCollision>
  void fillEventHistograms(TCollision collision, int& selGapSide)
  {
    // in case we want to push the analysis to Pb-Pb UPC
    int gapSide = collision.gapSide();
    // -1 --> Hadronic
    // 0 --> Single Gap - A side
    // 1 --> Single Gap - C side
    // 2 --> Double Gap - both A & C sides
    selGapSide = sgSelector.trueGap(collision, upcCuts.fv0Cut, upcCuts.ft0aCut, upcCuts.ft0cCut, upcCuts.zdcCut);
    histos.fill(HIST("hGapSide"), gapSide);
    histos.fill(HIST("hSelGapSide"), selGapSide);

    histos.fill(HIST("hSelGapSideVsPVz"), selGapSide, collision.posZ());
    histos.fill(HIST("hEventPVz"), collision.posZ());

    histos.fill(HIST("hEventOccupancy"), collision.occupancyInTime());
    histos.fill(HIST("hSelGapSideVsOccupancy"), selGapSide, collision.occupancyInTime());

    return;
  }

  template <typename TTrack>
  bool isTrackSelected(TTrack track)
  {
    //
    // acceptance cut
    //
    if (std::fabs(RecoDecay::eta(std::array{track.px(), track.py(), track.pz()})) > trackSelections.etaCut)
      return false;

    // PV contributor selection
    if (trackSelections.requirePVcontributor && !track.isPVContributor())
      return false;

    // dca XY to PV
    if (trackSelections.applyDCAptdepsel) { // apply pt dep. selection on DCAxy
      float dcaXYPtCut = 0.0105f + 0.0350f / pow(track.pt(), 1.1f);
      if (std::fabs(track.dcaXY()) > dcaXYPtCut)
        return false;
    } else {
      if (std::fabs(track.dcaXY()) > trackSelections.dcaxytopv)
        return false;
    }
    // dca Z to PV
    if (std::fabs(track.dcaZ()) > trackSelections.dcaztopv)
      return false;

    //
    // ITS quality flags
    //
    if (track.itsNCls() < trackSelections.minITSclusters)
      return false;
    if (track.itsChi2NCl() < trackSelections.minITSchi2clusters)
      return false;
    if (trackSelections.rejectITSafterburner && track.itsChi2NCl() < 0)
      return false;
    if (trackSelections.requireITSonly && track.tpcNClsCrossedRows() > 0)
      return false;

    //
    // TPC quality flags
    //
    if (track.tpcNClsCrossedRows() < trackSelections.minTPCrows)
      return false;
    if (track.tpcChi2NCl() < trackSelections.minTPCchi2clusters)
      return false;
    if (track.tpcNClsFindable() - track.tpcNClsFindableMinusFound() < trackSelections.minTPCclusters)
      return false;
    if (static_cast<float>(track.tpcNClsCrossedRows()) / static_cast<float>(track.tpcNClsFindable()) < trackSelections.minTPCrowsoverfindable)
      return false;
    if (trackSelections.skipTPConly && track.detectorMap() == o2::aod::track::TPC)
      return false;

    //
    // TPC PID
    //
    if (std::fabs(track.tpcNSigmaPr()) > trackSelections.tpcPidNsigmaCut)
      return false;

    //
    // TOF PID in NSigma
    // Bachelor track
    if (track.hasTOF()) {
      if (std::fabs(track.tofNSigmaPr()) > trackSelections.tofPidNsigmaCut)
        return false;
    }

    return true;
  }

  template <typename TTrack, typename TTrackMC>
  bool checkMCAssociation(TTrack track, TTrackMC trackMC)
  // MC association (if asked)
  {
    if (track.sign() * trackMC.pdgCode() != 2212)
      return false;
    if (!trackMC.isPhysicalPrimary())
      return false;
    return true;
  }

  template <typename TTrack>
  void fillQAplot(TTrack track, bool afterSel = false)
  { // fill QA information about proton/antiproton track
    if (afterSel) {
      if (track.sign() > 0) { // Proton Candidates after selections
        histos.fill(HIST("PPbar/Proton/hDCAxyToPV_aftersel"), track.dcaXY());
        histos.fill(HIST("PPbar/Proton/hDCAzToPV_aftersel"), track.dcaZ());
        histos.fill(HIST("PPbar/Proton/hTPCCrossedRows_aftersel"), track.tpcNClsCrossedRows());
        histos.fill(HIST("PPbar/Proton/hTPCNClusters_aftersel"), track.tpcNClsFindable() - track.tpcNClsFindableMinusFound());
        histos.fill(HIST("PPbar/Proton/hTPCChi2Clusters_aftersel"), track.tpcChi2NCl());
        histos.fill(HIST("PPbar/Proton/hTPCCrossedRowsOverFindable_aftersel"), static_cast<float>(track.tpcNClsCrossedRows()) / static_cast<float>(track.tpcNClsFindable()));
        histos.fill(HIST("PPbar/Proton/hITSNClusters_aftersel"), track.itsNCls());
        histos.fill(HIST("PPbar/Proton/hITSChi2Clusters_aftersel"), track.itsChi2NCl());
        histos.fill(HIST("PPbar/Proton/hTPCNsigma_aftersel"), track.tpcNSigmaPr());
        if (track.hasTOF())
          histos.fill(HIST("PPbar/Proton/hTOFNsigma_aftersel"), track.tofNSigmaPr());
        histos.fill(HIST("PPbar/Proton/h2dITSvsTPCpts_aftersel"), track.tpcNClsCrossedRows(), track.itsNCls());
      } else { // Anti Proton after selections
        histos.fill(HIST("PPbar/AntiProton/hDCAxyToPV_aftersel"), track.dcaXY());
        histos.fill(HIST("PPbar/AntiProton/hDCAzToPV_aftersel"), track.dcaZ());
        histos.fill(HIST("PPbar/AntiProton/hTPCCrossedRows_aftersel"), track.tpcNClsCrossedRows());
        histos.fill(HIST("PPbar/AntiProton/hTPCNClusters_aftersel"), track.tpcNClsFindable() - track.tpcNClsFindableMinusFound());
        histos.fill(HIST("PPbar/AntiProton/hTPCChi2Clusters_aftersel"), track.tpcChi2NCl());
        histos.fill(HIST("PPbar/AntiProton/hTPCCrossedRowsOverFindable_aftersel"), static_cast<float>(track.tpcNClsCrossedRows()) / static_cast<float>(track.tpcNClsFindable()));
        histos.fill(HIST("PPbar/AntiProton/hITSNClusters_aftersel"), track.itsNCls());
        histos.fill(HIST("PPbar/AntiProton/hITSChi2Clusters_aftersel"), track.itsChi2NCl());
        histos.fill(HIST("PPbar/AntiProton/hTPCNsigma_aftersel"), track.tpcNSigmaPr());
        if (track.hasTOF())
          histos.fill(HIST("PPbar/AntiProton/hTOFNsigma_aftersel"), track.tofNSigmaPr());
        histos.fill(HIST("PPbar/AntiProton/h2dITSvsTPCpts_aftersel"), track.tpcNClsCrossedRows(), track.itsNCls());
      }
    } else {
      if (track.sign() > 0) { // Proton Candidates before selections
        histos.fill(HIST("PPbar/Proton/hDCAxyToPV"), track.dcaXY());
        histos.fill(HIST("PPbar/Proton/hDCAzToPV"), track.dcaZ());
        histos.fill(HIST("PPbar/Proton/hTPCCrossedRows"), track.tpcNClsCrossedRows());
        histos.fill(HIST("PPbar/Proton/hTPCNClusters"), track.tpcNClsFindable() - track.tpcNClsFindableMinusFound());
        histos.fill(HIST("PPbar/Proton/hTPCChi2Clusters"), track.tpcChi2NCl());
        histos.fill(HIST("PPbar/Proton/hTPCCrossedRowsOverFindable"), static_cast<float>(track.tpcNClsCrossedRows()) / static_cast<float>(track.tpcNClsFindable()));
        histos.fill(HIST("PPbar/Proton/hITSNClusters"), track.itsNCls());
        histos.fill(HIST("PPbar/Proton/hITSChi2Clusters"), track.itsChi2NCl());
        histos.fill(HIST("PPbar/Proton/hTPCNsigma"), track.tpcNSigmaPr());
        if (track.hasTOF())
          histos.fill(HIST("PPbar/Proton/hTOFNsigma"), track.tofNSigmaPr());
        histos.fill(HIST("PPbar/Proton/h2dITSvsTPCpts"), track.tpcNClsCrossedRows(), track.itsNCls());
      } else {
        // Anti Proton before selections
        histos.fill(HIST("PPbar/AntiProton/hDCAxyToPV"), track.dcaXY());
        histos.fill(HIST("PPbar/AntiProton/hDCAzToPV"), track.dcaZ());
        histos.fill(HIST("PPbar/AntiProton/hTPCCrossedRows"), track.tpcNClsCrossedRows());
        histos.fill(HIST("PPbar/AntiProton/hTPCNClusters"), track.tpcNClsFindable() - track.tpcNClsFindableMinusFound());
        histos.fill(HIST("PPbar/AntiProton/hTPCChi2Clusters"), track.tpcChi2NCl());
        histos.fill(HIST("PPbar/AntiProton/hTPCCrossedRowsOverFindable"), static_cast<float>(track.tpcNClsCrossedRows()) / static_cast<float>(track.tpcNClsFindable()));
        histos.fill(HIST("PPbar/AntiProton/hITSNClusters"), track.itsNCls());
        histos.fill(HIST("PPbar/AntiProton/hITSChi2Clusters"), track.itsChi2NCl());
        histos.fill(HIST("PPbar/AntiProton/hTPCNsigma"), track.tpcNSigmaPr());
        if (track.hasTOF())
          histos.fill(HIST("PPbar/AntiProton/hTOFNsigma"), track.tofNSigmaPr());
        histos.fill(HIST("PPbar/AntiProton/h2dITSvsTPCpts"), track.tpcNClsCrossedRows(), track.itsNCls());
      }
    }
  }

  template <typename TTrack, typename TTrackMCs>
  void analyseTrackPairCandidate(TTrack proton, TTrack antiProton, TTrackMCs const& fullTrackMCs, uint8_t gapSide)
  // fill information related to the quarkonium mother
  {
    float pt = RecoDecay::pt(proton.px() + antiProton.px(), proton.py() + antiProton.py());
    float invmass = RecoDecay::m(std::array{std::array{proton.px(), proton.py(), proton.pz()}, std::array{antiProton.px(), antiProton.py(), antiProton.pz()}}, std::array{o2::constants::physics::MassProton, o2::constants::physics::MassProtonBar});
    float rapidity = RecoDecay::y(std::array{proton.px() + antiProton.px(), proton.py() + antiProton.py(), proton.pz() + antiProton.pz()}, invmass);

    // rapidity cut on the quarkonium mother
    if (!doMCAssociation && std::fabs(rapidity) > rapidityCut)
      return;

    // __________________________________________
    // main analysis
    if (doMCAssociation) {
      if constexpr (requires { proton.udMcParticle(); }) { // check if MC information is available
        auto protonMC = fullTrackMCs.iteratorAt(proton.udMcParticle().globalIndex());
        auto antiProtonMC = fullTrackMCs.iteratorAt(antiProton.udMcParticle().globalIndex());

        if (!protonMC.has_mothers())
          return;
        if (!antiProtonMC.has_mothers())
          return;

        float ptmc = RecoDecay::pt(protonMC.px() + antiProtonMC.px(), protonMC.py() + antiProtonMC.py());

        auto protonMothers = protonMC.template mothers_as<aod::UDMcParticles>();
        auto antiProtonMothers = antiProtonMC.template mothers_as<aod::UDMcParticles>();
        for (const auto& protonMother : protonMothers) {
          for (const auto& antiProtonMother : antiProtonMothers) {
            if (protonMother.globalIndex() != antiProtonMother.globalIndex()) {
              continue;
            }

            float rapiditymc = RecoDecay::y(std::array{protonMC.px() + antiProtonMC.px(), protonMC.py() + antiProtonMC.py(), protonMC.pz() + antiProtonMC.pz()}, pdgDB->Mass(protonMother.pdgCode()));

            if (std::fabs(rapiditymc) > rapidityCut)
              continue;

            if (protonMother.pdgCode() == 441 && protonMother.pdgCode() == antiProtonMother.pdgCode()) { // EtaC(1S)
              histos.fill(HIST("PPbar/h2dInvMassTrueEtaC1S"), ptmc, invmass);
            }
            if (protonMother.pdgCode() == 443 && protonMother.pdgCode() == antiProtonMother.pdgCode()) { // J/psi
              histos.fill(HIST("PPbar/h2dInvMassTrueJPsi"), ptmc, invmass);
            }
            if (protonMother.pdgCode() == 10441 && protonMother.pdgCode() == antiProtonMother.pdgCode()) { // ChiC0
              histos.fill(HIST("PPbar/h2dInvMassTrueChiC0"), ptmc, invmass);
            }
            if (protonMother.pdgCode() == 20443 && protonMother.pdgCode() == antiProtonMother.pdgCode()) { // ChiC1
              histos.fill(HIST("PPbar/h2dInvMassTrueChiC1"), ptmc, invmass);
            }
            if (protonMother.pdgCode() == 10443 && protonMother.pdgCode() == antiProtonMother.pdgCode()) { // hC
              histos.fill(HIST("PPbar/h2dInvMassTrueHC"), ptmc, invmass);
            }
            if (protonMother.pdgCode() == 445 && protonMother.pdgCode() == antiProtonMother.pdgCode()) { // ChiC2
              histos.fill(HIST("PPbar/h2dInvMassTrueChiC2"), ptmc, invmass);
            }
            if (protonMother.pdgCode() == 100441 && protonMother.pdgCode() == antiProtonMother.pdgCode()) { // EtaC(2S)
              histos.fill(HIST("PPbar/h2dInvMassTrueEtaC2S"), ptmc, invmass);
            }
            if (protonMother.pdgCode() == 100443 && protonMother.pdgCode() == antiProtonMother.pdgCode()) { // Psi(2S)
              histos.fill(HIST("PPbar/h2dInvMassTruePsi2S"), ptmc, invmass);
            }
          }
        }
      }
    }

    histos.fill(HIST("PPbar/h2dMassPPbar"), pt, invmass);
    if (gapSide == 0)
      histos.fill(HIST("PPbar/h2dMassPPbarSGA"), pt, invmass);
    else if (gapSide == 1)
      histos.fill(HIST("PPbar/h2dMassPPbarSGC"), pt, invmass);
    else if (gapSide == 2)
      histos.fill(HIST("PPbar/h2dMassPPbarDG"), pt, invmass);
    else
      histos.fill(HIST("PPbar/h2dMassPPbarHadronic"), pt, invmass);

    fillQAplot(proton, true);
    fillQAplot(antiProton, true);
  }

  template <typename TTracks, typename TTrackMCs>
  void buildProtonAntiProtonPairs(TTracks const& fullTracks, TTrackMCs const& fullMCTracks, std::vector<bool> selProtonIndices, std::vector<bool> selAntiProtonIndices, uint8_t gapSide)
  {
    // 1st loop over all protons
    for (const auto& proton : fullTracks) {
      // select only protons
      if (!selProtonIndices[proton.globalIndex() - fullTracks.offset()]) { // local index needed due to collisions grouping
        continue;
      }

      // 2nd loop over all protons
      for (const auto& antiProton : fullTracks) {
        // select only anti-protons
        if (!selAntiProtonIndices[antiProton.globalIndex() - fullTracks.offset()]) { // local index needed due to collisions grouping
          continue;
        }

        // check we don't look at the same protons
        if (proton.globalIndex() == antiProton.globalIndex()) {
          continue;
        }

        // form proton-antiproton pairs and fill histograms
        analyseTrackPairCandidate(proton, antiProton, fullMCTracks, gapSide);
      } // end antiProton loop
    } // end proton loop

    return;
  }

  // ______________________________________________________
  // Real data processing - no MC subscription
  void processRealData(soa::Join<aod::UDCollisions, aod::UDCollisionsSels, aod::UDZdcsReduced, aod::SGCollisions, aod::UDCollisionSelExtras>::iterator const& collision, soa::Join<aod::UDTracks, aod::UDTracksPID, aod::UDTracksExtra, aod::UDTracksFlags, aod::UDTracksDCA> const& fullTracks)
  {
    int selGapSide = -1; // only useful in case one wants to use this task in Pb-Pb UPC
    fillEventHistograms(collision, selGapSide);
    // __________________________________________
    // perform main analysis
    //
    // Look at tracks and tag those passing the selections
    std::vector<bool> selProtonIndices(fullTracks.size(), false);
    std::vector<bool> selAntiProtonIndices(fullTracks.size(), false);
    for (const auto& track : fullTracks) {
      if (track.sign() > 0) {
        selProtonIndices[track.globalIndex() - fullTracks.offset()] = isTrackSelected(track);
      } else {
        selAntiProtonIndices[track.globalIndex() - fullTracks.offset()] = isTrackSelected(track);
      }
      fillQAplot(track, false);
    } // end track loop

    // count the number of Proton and antiProton passsing the selections
    int nProtons = std::count(selProtonIndices.begin(), selProtonIndices.end(), true);
    int nAntiProtons = std::count(selAntiProtonIndices.begin(), selAntiProtonIndices.end(), true);

    // fill the histograms with the number of reconstructed protons/antiprotons per collision
    histos.fill(HIST("PPbar/h2dNbrOfProtonsVsSelGapSide"), selGapSide, nProtons);
    histos.fill(HIST("PPbar/h2dNbrOfAntiProtonsVsSelGapSide"), selGapSide, nAntiProtons);

    // Check the number of Protons and antiProtons
    // needs at least 1 of each
    if (!buildSameSignPairs && nProtons >= 1 && nAntiProtons >= 1) {
      buildProtonAntiProtonPairs(fullTracks, (TObject*)nullptr, selProtonIndices, selAntiProtonIndices, selGapSide);
    }
    if (buildSameSignPairs && nProtons > 1) {
      buildProtonAntiProtonPairs(fullTracks, (TObject*)nullptr, selProtonIndices, selProtonIndices, selGapSide);
    }
    if (buildSameSignPairs && nAntiProtons > 1) {
      buildProtonAntiProtonPairs(fullTracks, (TObject*)nullptr, selAntiProtonIndices, selAntiProtonIndices, selGapSide);
    }
  }

  // ______________________________________________________
  // Simulated processing (subscribes to MC information too)
  void processMonteCarlo(soa::Join<aod::UDCollisions, aod::UDCollisionsSels, aod::UDZdcsReduced, aod::SGCollisions, aod::UDCollisionSelExtras, aod::UDMcCollsLabels>::iterator const& collision, soa::Join<aod::UDTracks, aod::UDTracksPID, aod::UDTracksExtra, aod::UDTracksFlags, aod::UDTracksDCA, aod::UDMcTrackLabels> const& fullTracks, aod::UDMcCollisions const& /*mccollisions*/, aod::UDMcParticles const& fullMCTracks)
  {
    int selGapSide = -1; // only useful in case one wants to use this task in Pb-Pb UPC
    fillEventHistograms(collision, selGapSide);
    // __________________________________________
    // perform main analysis
    //
    // Look at tracks and tag those passing the selections
    std::vector<bool> selProtonIndices(fullTracks.size(), false);
    std::vector<bool> selAntiProtonIndices(fullTracks.size(), false);
    for (const auto& track : fullTracks) {
      if (!track.has_udMcParticle())
        continue;

      auto trackMC = fullMCTracks.iteratorAt(track.udMcParticle().globalIndex());

      if (track.sign() > 0) {
        selProtonIndices[track.globalIndex() - fullTracks.offset()] = isTrackSelected(track) && (!doMCAssociation || checkMCAssociation(track, trackMC));
      } else {
        selAntiProtonIndices[track.globalIndex() - fullTracks.offset()] = isTrackSelected(track) && (!doMCAssociation || checkMCAssociation(track, trackMC));
      }
      fillQAplot(track, false);
    } // end track loop

    // count the number of Proton and antiProton passsing the selections
    int nProtons = std::count(selProtonIndices.begin(), selProtonIndices.end(), true);
    int nAntiProtons = std::count(selAntiProtonIndices.begin(), selAntiProtonIndices.end(), true);

    // fill the histograms with the number of reconstructed protons/antiprotons per collision
    histos.fill(HIST("PPbar/h2dNbrOfProtonsVsSelGapSide"), selGapSide, nProtons);
    histos.fill(HIST("PPbar/h2dNbrOfAntiProtonsVsSelGapSide"), selGapSide, nAntiProtons);

    // Check the number of Protons and antiProtons
    // needs at least 1 of each
    if (!buildSameSignPairs && nProtons >= 1 && nAntiProtons >= 1) {
      buildProtonAntiProtonPairs(fullTracks, fullMCTracks, selProtonIndices, selAntiProtonIndices, selGapSide);
    }
    if (buildSameSignPairs && nProtons > 1) {
      buildProtonAntiProtonPairs(fullTracks, fullMCTracks, selProtonIndices, selProtonIndices, selGapSide);
    }
    if (buildSameSignPairs && nAntiProtons > 1) {
      buildProtonAntiProtonPairs(fullTracks, fullMCTracks, selAntiProtonIndices, selAntiProtonIndices, selGapSide);
    }
  }

  PROCESS_SWITCH(upcQuarkoniaCentralBarrel, processRealData, "process as if real data", true);
  PROCESS_SWITCH(upcQuarkoniaCentralBarrel, processMonteCarlo, "process as if MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<upcQuarkoniaCentralBarrel>(cfgc)};
}
