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

/// \author Alice Caluisi <alice.caluisi@cern.ch>
/// \since July 2023
/// \author Johanna Lömker <johanna.lomker@cern.ch>
/// \since  2023-10-02
///  \brief Task producing jet tracking qa histograms
///
#include "PWGJE/DataModel/TrackJetQa.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include <Framework/Configurable.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <math.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct TrackJetQa {
  Configurable<double> ValVtx{"ValVtx", 10, "Value of the vertex position"};
  Configurable<float> ValCutEta{"ValCutEta", 0.8f, "Eta range for tracks"};
  Configurable<float> minPt{"minPt", 0.15f, "minimum pT for tracks"};
  Configurable<float> maxPt{"maxPt", 10e10, "maximum pT for tracks"};
  Configurable<bool> fillMultiplicity{"fillMultiplicity", true, "To fill multiplicity and centrality histograms"};

  Configurable<bool> globalTrack{"globalTrack", false, "to enable the isGlobalTrack() selection"};
  Configurable<bool> globalTrackWoPtEta{"globalTrackWoPtEta", false, "to enable the isGlobalTrackWoPtEta() selection"};
  Configurable<bool> globalTrackWoDCA{"globalTrackWoDCA", false, "to enable the isGlobalTrackWoDCA() selection"};
  Configurable<bool> customTrack{"customTrack", true, "to enable the trackselection based on the customTrack cuts (selected via configurables)"};

  // Custom track cuts for the cut variation study
  TrackSelection customTrackCuts;
  Configurable<int> itsPattern{"itsPattern", 2, "0 = Run3ITSibAny, 1 = Run3ITSibTwo, 2 = Run3ITSallAny, 3 = Run3ITSall7Layers"};
  Configurable<bool> requireITS{"requireITS", true, "Additional cut on the ITS requirement"};
  Configurable<bool> requireTPC{"requireTPC", true, "Additional cut on the TPC requirement"};
  Configurable<bool> requireGoldenChi2{"requireGoldenChi2", true, "Additional cut on the GoldenChi2"};
  Configurable<float> minNCrossedRowsTPC{"minNCrossedRowsTPC", 60.f, "Additional cut on the minimum number of crossed rows in the TPC"};
  Configurable<float> minNCrossedRowsOverFindableClustersTPC{"minNCrossedRowsOverFindableClustersTPC", 0.7f, "Additional cut on the minimum value of the ratio between crossed rows and findable clusters in the TPC"};
  Configurable<float> maxChi2PerClusterTPC{"maxChi2PerClusterTPC", 7.f, "Additional cut on the maximum value of the chi2 per cluster in the TPC"};
  Configurable<float> maxChi2PerClusterITS{"maxChi2PerClusterITS", 36.f, "Additional cut on the maximum value of the chi2 per cluster in the ITS"};
  Configurable<float> maxDcaXY{"maxDcaXY", 0.25f, "Cut on the maximum value of the DCA xy "};
  Configurable<float> maxDcaZ{"maxDcaZ", 3.f, "Additional cut on the maximum value of the DCA z"};
  Configurable<float> minTPCNClsFound{"minTPCNClsFound", 0.f, "Additional cut on the minimum value of the number of found clusters in the TPC"};

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  ConfigurableAxis binsMultiplicity{"binsMultiplicity", {2000, 0, 2000}, "Binning for multiplicity"};
  ConfigurableAxis binsMultPV{"binsMultNTracksPV", {10000, 0, 50000}, "Binning for the multNTracksPV axis"};
  ConfigurableAxis binsPercentile{"binsPercentile", {100, 0, 100}, "Binning for percentiles"};
  ConfigurableAxis binsVtx{"binsVtx", {200, -20, 20}, "Binning for the vertex position z axis"};
  ConfigurableAxis binsPt{"binsPt", {200, 0, 200}, "Binning for the pT axis"};
  ConfigurableAxis binsSigma1OverPt{"binsSigma1OverPt", {100, 0, 1}, "Binning for the sigma 1 over pT * pT"};
  ConfigurableAxis binsPhi{"binsPhi", {180, 0, 2 * M_PI}, "Binning for the phi axis"};
  ConfigurableAxis binsEta{"binsEta", {100, -1, 1}, "Binning for the eta axis"};
  ConfigurableAxis binsTrackXY{"binsTrackXY", {100, -0.5, 0.5}, "Binning for the x and y track position at dca in local coordinate system axis"};
  ConfigurableAxis binsTrackZ{"binsTrackZ", {100, -11, 11}, "Binning for the z track position at dca in local coordinate system axis"};
  ConfigurableAxis binsRot{"binsRot", {36, -M_PI, M_PI}, "Binning for the rotation angle axis"};
  ConfigurableAxis binsSignedPt{"binsSignedPt", {200, -8, 8}, "Binning for the q over pt axis"};
  ConfigurableAxis binsDcaXY{"binsDcaXY", {100, -0.5, 0.5}, "Binning for the dcaXY axis"};
  ConfigurableAxis binsDcaZ{"binsDcaZ", {100, -5, 5}, "Binning for the dcaXY axis"};
  ConfigurableAxis binsLength{"binsLength", {200, 0, 1000}, "Binning for the track length axis"};

  Filter ptFilter = (aod::track::pt >= minPt) && (aod::track::pt <= maxPt);
  Filter etaFilter = (aod::track::eta <= ValCutEta) && (aod::track::eta >= -ValCutEta);

  void init(o2::framework::InitContext&)
  {
    if (customTrack) {
      // Custom track cuts
      LOG(info) << "Using custom track cuts from values:";
      LOG(info) << "\trequireITS=" << requireITS.value;
      LOG(info) << "\trequireTPC=" << requireTPC.value;
      LOG(info) << "\trequireGoldenChi2=" << requireGoldenChi2.value;
      LOG(info) << "\tmaxChi2PerClusterTPC=" << maxChi2PerClusterTPC.value;
      LOG(info) << "\tminNCrossedRowsTPC=" << minNCrossedRowsTPC.value;
      LOG(info) << "\tminTPCNClsFound=" << minTPCNClsFound.value;
      LOG(info) << "\tmaxChi2PerClusterITS=" << maxChi2PerClusterITS.value;
      LOG(info) << "\tRequireHitsInITSLayers=" << maxChi2PerClusterITS.value;
      LOG(info) << "\tmaxDcaXY=" << maxDcaXY.value;
      LOG(info) << "\tmaxDcaZ=" << maxDcaZ.value;
      LOG(info) << "\tminPt=" << minPt.value;
      LOG(info) << "\tmaxPt=" << maxPt.value;
      LOG(info) << "\tmaxEta=" << ValCutEta.value;

      customTrackCuts = getGlobalTrackSelectionRun3ITSMatch(itsPattern.value);
      LOG(info) << "Customizing track cuts:";
      customTrackCuts.SetEtaRange(-ValCutEta.value, ValCutEta.value);
      customTrackCuts.SetPtRange(minPt.value, maxPt.value);
      customTrackCuts.SetRequireITSRefit(requireITS.value);
      customTrackCuts.SetRequireTPCRefit(requireTPC.value);
      customTrackCuts.SetRequireGoldenChi2(requireGoldenChi2.value);
      customTrackCuts.SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC.value);
      customTrackCuts.SetMaxChi2PerClusterITS(maxChi2PerClusterITS.value);
      customTrackCuts.SetMinNCrossedRowsTPC(minNCrossedRowsTPC.value);
      customTrackCuts.SetMinNClustersTPC(minTPCNClsFound.value);
      // customTrackCuts.SetRequireHitsInITSLayers(nHits.value, {0, 1}); // one hit in any SPD layer (#hits, {layer0, layer1,...})
      customTrackCuts.SetMinNCrossedRowsOverFindableClustersTPC(minNCrossedRowsOverFindableClustersTPC.value);
      customTrackCuts.SetMaxDcaXY(maxDcaXY.value);
      customTrackCuts.SetMaxDcaZ(maxDcaZ.value);
      customTrackCuts.print();
    }
    if (globalTrack) {
      LOG(info) << "Using globalTracks";
    }
    if (globalTrackWoPtEta) {
      LOG(info) << "Using globalTracksWoPtEta";
    }
    if (globalTrackWoDCA) {
      LOG(info) << "Using globalTracksWoDCA";
    } else {
      LOG(info) << "No trackselection enabled !";
    }

    // Common axes
    const AxisSpec axisPercentileFT0M{binsPercentile, "Centrality FT0M"};
    const AxisSpec axisPercentileFT0A{binsPercentile, "Centrality FT0A"};
    const AxisSpec axisPercentileFT0C{binsPercentile, "Centrality FT0C"};
    const AxisSpec axisMultiplicityPV{binsMultPV, "Multiplicity N TracksPV"};
    const AxisSpec axisMultiplicityTracks{binsMultiplicity, "Multiplicity tracks.size()"};
    const AxisSpec axisMultiplicityFT0M{binsMultiplicity, "Multiplicity FT0M"};
    const AxisSpec axisMultiplicityFT0A{binsMultiplicity, "Multiplicity FT0A"};
    const AxisSpec axisMultiplicityFT0C{binsMultiplicity, "Multiplicity FT0C"};
    const AxisSpec axisVtx{binsVtx, "#it{Vtx}_{z} [cm]"};
    const AxisSpec axisPt{binsPt, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec axisPhi{binsPhi, "#phi [rad]"};
    const AxisSpec axisEta{binsEta, " #eta"};
    const AxisSpec axisSigma1OverPt{binsSigma1OverPt, "#sigma(1/#it{p}_{T})*#it{p}_{T}"};
    const AxisSpec axisTrackX{binsTrackXY, "track #it{x} [cm]"};
    const AxisSpec axisTrackY{binsTrackXY, "track #it{y} [cm]"};
    const AxisSpec axisTrackZ{binsTrackZ, "track #it{z} [cm]"};
    const AxisSpec axisRotation{binsRot, "#alpha [rad]"};
    const AxisSpec axisSignedPt{binsSignedPt, "#it{q}/#it{p}_{T}"};
    const AxisSpec axisDcaXY{binsDcaXY, "#it{dcaXY} [cm]"};
    const AxisSpec axisDcaZ{binsDcaZ, "#it{dcaZ} [cm]"};
    const AxisSpec axisTrackLength{binsLength, "#it{Length} [cm]"};

    // event property histograms
    histos.add("EventProp/collisionVtxZ", "Collsion Vertex Z position", HistType::kTHnSparseD, {axisVtx, axisPercentileFT0A, axisPercentileFT0C});
    histos.add("EventProp/collisionVtxZnoSel", "Collsion Vertex Z position without event selection", HistType::kTHnSparseD, {axisVtx, axisPercentileFT0A, axisPercentileFT0C});
    histos.add("EventProp/collisionVtxZSel8", "Collsion Vertex Z position with event selection", HistType::kTHnSparseD, {axisVtx, axisPercentileFT0A, axisPercentileFT0C});
    histos.add("EventProp/rejectedCollId", "CollisionId of collisions that did not pass the event selection; collisionId; number of entries", HistType::kTH1F, {{10, 0, 5}});
    histos.add("EventProp/MultCorrelations", "Multiplicity and Centrality Correlations", HistType::kTHnSparseD, {axisPercentileFT0A, axisPercentileFT0C, axisPercentileFT0M, axisMultiplicityFT0A, axisMultiplicityFT0C, axisMultiplicityFT0M, axisMultiplicityPV, axisMultiplicityTracks});

    histos.add("TrackEventPar/MultCorrelations", "Sigma1Pt*pT vs Multiplicity and Centrality Correlations", HistType::kTHnSparseD, {axisPt, axisSigma1OverPt, axisPercentileFT0A, axisPercentileFT0C, axisPercentileFT0M, axisMultiplicityFT0A, axisMultiplicityFT0C, axisMultiplicityFT0M, axisMultiplicityPV, axisMultiplicityTracks});

    // kinetic histograms
    histos.add("Kine/pt", "#it{p}_{T}", HistType::kTHnSparseD, {axisPt, axisSigma1OverPt, axisPercentileFT0A, axisPercentileFT0C});
    histos.add("Kine/pt_TRD", "#it{p}_{T} if track has a TRD match", HistType::kTHnSparseD, {axisPt, axisSigma1OverPt, axisPercentileFT0A, axisPercentileFT0C});
    histos.add("Kine/EtaPhiPt", "Correlation of #eta #phi and #it{p}_{T}", HistType::kTHnSparseD, {axisPt, axisSigma1OverPt, axisEta, axisPhi, axisPercentileFT0A, axisPercentileFT0C});

    // track parameter histograms - add sigma 1pt to all of them ! then cp this part to below
    histos.add("TrackPar/xyz", "track #it{x}, #it{y}, #it{z} position at dca in local coordinate system", HistType::kTHnSparseD, {axisPt, axisSigma1OverPt, axisTrackX, axisTrackY, axisTrackZ, axisPercentileFT0A, axisPercentileFT0C});
    histos.add("TrackPar/alpha", "rotation angle of local wrt. global coordinate system", HistType::kTHnSparseD, {axisPt, axisSigma1OverPt, axisRotation, axisPercentileFT0A, axisPercentileFT0C});
    histos.add("TrackPar/signed1Pt", "track signed 1/#it{p}_{T}", HistType::kTHnSparseD, {axisPt, axisSigma1OverPt, axisSignedPt, axisPercentileFT0A, axisPercentileFT0C});
    histos.add("TrackPar/snp", "sinus of track momentum azimuthal angle (snp)", HistType::kTHnSparseD, {axisPt, axisSigma1OverPt, {50, -0.5, 0.5, "snp"}, axisPercentileFT0A, axisPercentileFT0C});
    histos.add("TrackPar/tgl", "tangent of the track momentum dip angle (tgl)", HistType::kTHnSparseD, {axisPt, axisSigma1OverPt, {200, -1., 1., "tgl"}, axisPercentileFT0A, axisPercentileFT0C});
    histos.add("TrackPar/flags", "track flag;#it{p}_{T} [GeV/c];flag bit", {HistType::kTH2F, {{200, 0, 200}, {64, -0.5, 63.5}}});
    histos.add("TrackPar/dcaXY", "distance of closest approach in #it{xy} plane", HistType::kTHnSparseD, {axisPt, axisSigma1OverPt, axisDcaXY, axisPercentileFT0A, axisPercentileFT0C});
    histos.add("TrackPar/dcaZ", "distance of closest approach in #it{z}", HistType::kTHnSparseD, {axisPt, axisSigma1OverPt, axisDcaZ, axisPercentileFT0A, axisPercentileFT0C});
    histos.add("TrackPar/length", "track length in cm", HistType::kTHnSparseD, {axisPt, axisSigma1OverPt, axisTrackLength, axisPercentileFT0A, axisPercentileFT0C});
    histos.add("TrackPar/Sigma1Pt", "uncertainty over #it{p}_{T}", HistType::kTHnSparseD, {axisPt, axisSigma1OverPt, axisPercentileFT0A, axisPercentileFT0C});
    histos.add("TrackPar/Sigma1Pt_hasTRD", "uncertainty over #it{p}_{T} for tracks with TRD", HistType::kTHnSparseD, {axisPt, axisSigma1OverPt, axisPercentileFT0A, axisPercentileFT0C});
    histos.add("TrackPar/Sigma1Pt_hasNoTRD", "uncertainty over #it{p}_{T} for tracks without TRD", HistType::kTHnSparseD, {axisPt, axisSigma1OverPt, axisPercentileFT0A, axisPercentileFT0C});
    histos.add("TrackPar/Sigma1Pt_Layer1", "uncertainty over #it{p}_{T} with 1st ITS layer active", HistType::kTHnSparseD, {axisPt, axisSigma1OverPt, axisPercentileFT0A, axisPercentileFT0C});
    histos.add("TrackPar/Sigma1Pt_Layer2", "uncertainty over #it{p}_{T} with 2nd ITS layer active", HistType::kTHnSparseD, {axisPt, axisSigma1OverPt, axisPercentileFT0A, axisPercentileFT0C});
    histos.add("TrackPar/Sigma1Pt_Layer3", "uncertainty over #it{p}_{T} with 3rd ITS layer active", HistType::kTHnSparseD, {axisPt, axisSigma1OverPt, axisPercentileFT0A, axisPercentileFT0C});
    histos.add("TrackPar/Sigma1Pt_Layer4", "uncertainty over #it{p}_{T} with 4th ITS layer active", HistType::kTHnSparseD, {axisPt, axisSigma1OverPt, axisPercentileFT0A, axisPercentileFT0C});
    histos.add("TrackPar/Sigma1Pt_Layer5", "uncertainty over #it{p}_{T} with 5th ITS layer active", HistType::kTHnSparseD, {axisPt, axisSigma1OverPt, axisPercentileFT0A, axisPercentileFT0C});
    histos.add("TrackPar/Sigma1Pt_Layer6", "uncertainty over #it{p}_{T} with 6th ITS layer active", HistType::kTHnSparseD, {axisPt, axisSigma1OverPt, axisPercentileFT0A, axisPercentileFT0C});
    histos.add("TrackPar/Sigma1Pt_Layer7", "uncertainty over #it{p}_{T} with 7th ITS layer active", HistType::kTHnSparseD, {axisPt, axisSigma1OverPt, axisPercentileFT0A, axisPercentileFT0C});
    histos.add("TrackPar/Sigma1Pt_Layers12", "uncertainty over #it{p}_{T} with 1st and 2nd ITS layers active", HistType::kTHnSparseD, {axisPt, axisSigma1OverPt, axisPercentileFT0A, axisPercentileFT0C});
    histos.add("TrackPar/Sigma1Pt_Layers34", "uncertainty over #it{p}_{T} with 3rd and 4th ITS layers active", HistType::kTHnSparseD, {axisPt, axisSigma1OverPt, axisPercentileFT0A, axisPercentileFT0C});
    histos.add("TrackPar/Sigma1Pt_Layers45", "uncertainty over #it{p}_{T} with 4th and 5th ITS layers active", HistType::kTHnSparseD, {axisPt, axisSigma1OverPt, axisPercentileFT0A, axisPercentileFT0C});
    histos.add("TrackPar/Sigma1Pt_Layers56", "uncertainty over #it{p}_{T} with 5th and 6th ITS layers active", HistType::kTHnSparseD, {axisPt, axisSigma1OverPt, axisPercentileFT0A, axisPercentileFT0C});
    histos.add("TrackPar/Sigma1Pt_Layers67", "uncertainty over #it{p}_{T} with 6th and 7th ITS layers active", HistType::kTHnSparseD, {axisPt, axisSigma1OverPt, axisPercentileFT0A, axisPercentileFT0C});
    histos.add("TrackPar/Sigma1Pt_Layers1or2and7", "uncertainty over #it{p}_{T} with 1st or 2nd and 7th ITS layers active", HistType::kTHnSparseD, {axisPt, axisSigma1OverPt, axisPercentileFT0A, axisPercentileFT0C});
    histos.add("TrackPar/Sigma1Pt_Layers1or2and6or7", "uncertainty over #it{p}_{T} with 1st or 2nd and 6th or 7th ITS layers active", HistType::kTHnSparseD, {axisPt, axisSigma1OverPt, axisPercentileFT0A, axisPercentileFT0C});
    histos.add("TrackPar/Sigma1Pt_Layers456", "uncertainty over #it{p}_{T} with 4th, 5th and 6th ITS layers active", HistType::kTHnSparseD, {axisPt, axisSigma1OverPt, axisPercentileFT0A, axisPercentileFT0C});

    // ITS histograms
    histos.add("ITS/itsNCls", "number of found ITS clusters", HistType::kTHnSparseD, {axisPt, axisSigma1OverPt, {8, -0.5, 7.5, "# clusters ITS"}, axisPercentileFT0A, axisPercentileFT0C});
    histos.add("ITS/itsChi2NCl", "chi2 per ITS cluster", HistType::kTHnSparseD, {axisPt, axisSigma1OverPt, {100, 0, 40, "chi2 / cluster ITS"}, axisPercentileFT0A, axisPercentileFT0C});
    histos.add("ITS/itsHits", "hitmap ITS", HistType::kTHnSparseD, {axisPt, axisSigma1OverPt, {7, -0.5, 6.5, "layer ITS"}, axisPercentileFT0A, axisPercentileFT0C});

    // TPC histograms
    histos.add("TPC/tpcNClsFindable", "number of findable TPC clusters", HistType::kTHnSparseD, {axisPt, axisSigma1OverPt, {165, -0.5, 164.5, "# findable clusters TPC"}, axisPercentileFT0A, axisPercentileFT0C});
    histos.add("TPC/tpcNClsFound", "number of found TPC clusters", HistType::kTHnSparseD, {axisPt, axisSigma1OverPt, {165, -0.5, 164.5, "# clusters TPC"}, axisPercentileFT0A, axisPercentileFT0C});
    histos.add("TPC/tpcNClsShared", "number of shared TPC clusters", HistType::kTHnSparseD, {axisPt, axisSigma1OverPt, {165, -0.5, 164.5, "# shared clusters TPC"}, axisPercentileFT0A, axisPercentileFT0C});
    histos.add("TPC/tpcNClsCrossedRows", "number of crossed TPC rows", HistType::kTHnSparseD, {axisPt, axisSigma1OverPt, {165, -0.5, 164.5, "# crossed rows TPC"}, axisPercentileFT0A, axisPercentileFT0C});
    histos.add("TPC/tpcFractionSharedCls", "fraction of shared TPC clusters", HistType::kTHnSparseD, {axisPt, axisSigma1OverPt, {100, 0., 1., "fraction shared clusters TPC"}, axisPercentileFT0A, axisPercentileFT0C});
    histos.add("TPC/tpcCrossedRowsOverFindableCls", "crossed TPC rows over findable clusters", HistType::kTHnSparseD, {axisPt, axisSigma1OverPt, {120, 0.0, 1.2, "crossed rows / findable clusters TPC"}, axisPercentileFT0A, axisPercentileFT0C});
    histos.add("TPC/tpcChi2NCl", "chi2 per cluster in TPC", HistType::kTHnSparseD, {axisPt, axisSigma1OverPt, {100, 0, 10, "chi2 / cluster TPC"}, axisPercentileFT0A, axisPercentileFT0C});

    histos.print();
  }

  template <typename eventInfo>
  bool checkEventSelection(eventInfo const& collision)
  {
    // fill event property variables
    histos.fill(HIST("EventProp/collisionVtxZnoSel"), collision.posZ(), collision.centFT0A(), collision.centFT0C());
    if (!collision.sel8()) {
      histos.fill(HIST("EventProp/rejectedCollId"), 2);
      return false;
    }
    histos.fill(HIST("EventProp/collisionVtxZSel8"), collision.posZ(), collision.centFT0A(), collision.centFT0C());
    if (fabs(collision.posZ()) > ValVtx) {
      histos.fill(HIST("EventProp/rejectedCollId"), 3);
      return false;
    }
    histos.fill(HIST("EventProp/collisionVtxZ"), collision.posZ(), collision.centFT0A(), collision.centFT0C());
    return true;
  }

  template <typename Tracks>
  bool checkTrackSelection(Tracks const& track)
  {
    // check track selection
    if ((globalTrack == true) && (!track.isGlobalTrack())) {
      return false;
    }
    if ((globalTrackWoDCA == true) && (!track.isGlobalTrackWoDCA())) {
      return false;
    }
    if ((globalTrackWoPtEta == true) && (!track.isGlobalTrackWoPtEta())) {
      return false;
    }
    if ((customTrack == true) && (!customTrackCuts.IsSelected(track))) {
      return false;
    }
    return true;
  }

  template <typename Tracks, typename eventInfo>
  void fillTrackQa(Tracks const& track, eventInfo const& collision)
  {
    // fill kinematic variables
    histos.fill(HIST("Kine/pt"), track.pt(), track.sigma1Pt() * track.pt(), collision.centFT0A(), collision.centFT0C());
    if (track.hasTRD()) {
      histos.fill(HIST("Kine/pt_TRD"), track.pt(), track.sigma1Pt() * track.pt(), collision.centFT0A(), collision.centFT0C());
      histos.fill(HIST("TrackPar/Sigma1Pt_hasTRD"), track.pt(), track.sigma1Pt() * track.pt(), track.sigma1Pt() * track.pt(), collision.centFT0A(), collision.centFT0C());
    }
    if (!track.hasTRD()) {
      histos.fill(HIST("TrackPar/Sigma1Pt_hasNoTRD"), track.pt(), track.sigma1Pt() * track.pt(), track.sigma1Pt() * track.pt(), collision.centFT0A(), collision.centFT0C());
    }
    histos.fill(HIST("Kine/EtaPhiPt"), track.pt(), track.sigma1Pt() * track.pt(), track.eta(), track.phi(), collision.centFT0A(), collision.centFT0C());
    // fill track parameter variables
    histos.fill(HIST("TrackPar/alpha"), track.pt(), track.sigma1Pt() * track.pt(), track.alpha(), collision.centFT0A(), collision.centFT0C());
    histos.fill(HIST("TrackPar/xyz"), track.pt(), track.sigma1Pt() * track.pt(), track.x(), track.y(), track.z(), collision.centFT0A(), collision.centFT0C());
    histos.fill(HIST("TrackPar/signed1Pt"), track.pt(), track.sigma1Pt() * track.pt(), track.signed1Pt(), collision.centFT0A(), collision.centFT0C());
    histos.fill(HIST("TrackPar/snp"), track.pt(), track.sigma1Pt() * track.pt(), track.snp(), collision.centFT0A(), collision.centFT0C());
    histos.fill(HIST("TrackPar/tgl"), track.pt(), track.sigma1Pt() * track.pt(), track.tgl(), collision.centFT0A(), collision.centFT0C());
    for (unsigned uint32_t int i = 0; i < 32; i++) {
      if (track.flags() & (1 << i)) {
        histos.fill(HIST("TrackPar/flags"), track.pt(), track.sigma1Pt() * track.pt(), i);
      }
    }
    histos.fill(HIST("TrackPar/dcaXY"), track.pt(), track.sigma1Pt() * track.pt(), track.dcaXY(), collision.centFT0A(), collision.centFT0C());
    histos.fill(HIST("TrackPar/dcaZ"), track.pt(), track.sigma1Pt() * track.pt(), track.dcaZ(), collision.centFT0A(), collision.centFT0C());
    histos.fill(HIST("TrackPar/length"), track.pt(), track.sigma1Pt() * track.pt(), track.length(), collision.centFT0A(), collision.centFT0C());
    histos.fill(HIST("TrackPar/Sigma1Pt"), track.pt(), track.sigma1Pt() * track.pt(), collision.centFT0A(), collision.centFT0C());
    //// check the uncertainty over pT activating several ITS layers
    bool firstLayerActive = track.itsClusterMap() & (1 << 0);
    bool secondLayerActive = track.itsClusterMap() & (1 << 1);
    bool thirdLayerActive = track.itsClusterMap() & (1 << 2);
    bool fourthLayerActive = track.itsClusterMap() & (1 << 3);
    bool fifthLayerActive = track.itsClusterMap() & (1 << 4);
    bool sixthLayerActive = track.itsClusterMap() & (1 << 5);
    bool seventhLayerActive = track.itsClusterMap() & (1 << 6);

    if (firstLayerActive) {
      histos.fill(HIST("TrackPar/Sigma1Pt_Layer1"), track.pt(), track.sigma1Pt() * track.pt(), collision.centFT0A(), collision.centFT0C());
    }
    if (secondLayerActive) {
      histos.fill(HIST("TrackPar/Sigma1Pt_Layer2"), track.pt(), track.sigma1Pt() * track.pt(), collision.centFT0A(), collision.centFT0C());
    }
    if (thirdLayerActive) {
      histos.fill(HIST("TrackPar/Sigma1Pt_Layer3"), track.pt(), track.sigma1Pt() * track.pt(), collision.centFT0A(), collision.centFT0C());
    }
    if (fourthLayerActive) {
      histos.fill(HIST("TrackPar/Sigma1Pt_Layer4"), track.pt(), track.sigma1Pt() * track.pt(), collision.centFT0A(), collision.centFT0C());
    }
    if (fifthLayerActive) {
      histos.fill(HIST("TrackPar/Sigma1Pt_Layer5"), track.pt(), track.sigma1Pt() * track.pt(), collision.centFT0A(), collision.centFT0C());
    }
    if (sixthLayerActive) {
      histos.fill(HIST("TrackPar/Sigma1Pt_Layer6"), track.pt(), track.sigma1Pt() * track.pt(), collision.centFT0A(), collision.centFT0C());
    }
    if (seventhLayerActive) {
      histos.fill(HIST("TrackPar/Sigma1Pt_Layer7"), track.pt(), track.sigma1Pt() * track.pt(), collision.centFT0A(), collision.centFT0C());
    }
    if (firstLayerActive && secondLayerActive) {
      histos.fill(HIST("TrackPar/Sigma1Pt_Layers12"), track.pt(), track.sigma1Pt() * track.pt(), collision.centFT0A(), collision.centFT0C());
    }
    if (thirdLayerActive && fourthLayerActive) {
      histos.fill(HIST("TrackPar/Sigma1Pt_Layers34"), track.pt(), track.sigma1Pt() * track.pt(), collision.centFT0A(), collision.centFT0C());
    }
    if (fourthLayerActive && fifthLayerActive) {
      histos.fill(HIST("TrackPar/Sigma1Pt_Layers45"), track.pt(), track.sigma1Pt() * track.pt(), collision.centFT0A(), collision.centFT0C());
    }
    if (fifthLayerActive && sixthLayerActive) {
      histos.fill(HIST("TrackPar/Sigma1Pt_Layers56"), track.pt(), track.sigma1Pt() * track.pt(), collision.centFT0A(), collision.centFT0C());
    }
    if (sixthLayerActive && seventhLayerActive) {
      histos.fill(HIST("TrackPar/Sigma1Pt_Layers67"), track.pt(), track.sigma1Pt() * track.pt(), collision.centFT0A(), collision.centFT0C());
    }
    if ((firstLayerActive || secondLayerActive) && seventhLayerActive) {
      histos.fill(HIST("TrackPar/Sigma1Pt_Layers1or2and7"), track.pt(), track.sigma1Pt() * track.pt(), collision.centFT0A(), collision.centFT0C());
    }
    if ((firstLayerActive || secondLayerActive) && (sixthLayerActive || seventhLayerActive)) {
      histos.fill(HIST("TrackPar/Sigma1Pt_Layers1or2and6or7"), track.pt(), track.sigma1Pt() * track.pt(), collision.centFT0A(), collision.centFT0C());
    }
    if (fourthLayerActive && fifthLayerActive && sixthLayerActive) {
      histos.fill(HIST("TrackPar/Sigma1Pt_Layers456"), track.pt(), track.sigma1Pt() * track.pt(), collision.centFT0A(), collision.centFT0C());
    }
    // fill ITS variables
    histos.fill(HIST("ITS/itsNCls"), track.pt(), track.sigma1Pt() * track.pt(), track.itsNCls(), collision.centFT0A(), collision.centFT0C());
    histos.fill(HIST("ITS/itsChi2NCl"), track.pt(), track.sigma1Pt() * track.pt(), track.itsChi2NCl(), collision.centFT0A(), collision.centFT0C());
    for (unsigned int i = 0; i < 7; i++) {
      if (track.itsClusterMap() & (1 << i)) {
        histos.fill(HIST("ITS/itsHits"), track.pt(), track.sigma1Pt() * track.pt(), i, collision.centFT0A(), collision.centFT0C());
      }
    }
    // fill TPC variables
    histos.fill(HIST("TPC/tpcNClsFindable"), track.pt(), track.sigma1Pt() * track.pt(), track.tpcNClsFindable(), collision.centFT0A(), collision.centFT0C());
    histos.fill(HIST("TPC/tpcNClsFound"), track.pt(), track.sigma1Pt() * track.pt(), track.tpcNClsFound(), collision.centFT0A(), collision.centFT0C());
    histos.fill(HIST("TPC/tpcNClsShared"), track.pt(), track.sigma1Pt() * track.pt(), track.tpcNClsShared(), collision.centFT0A(), collision.centFT0C());
    histos.fill(HIST("TPC/tpcNClsCrossedRows"), track.pt(), track.sigma1Pt() * track.pt(), track.tpcNClsCrossedRows(), collision.centFT0A(), collision.centFT0C());
    histos.fill(HIST("TPC/tpcCrossedRowsOverFindableCls"), track.pt(), track.sigma1Pt() * track.pt(), track.tpcCrossedRowsOverFindableCls(), collision.centFT0A(), collision.centFT0C());
    histos.fill(HIST("TPC/tpcFractionSharedCls"), track.pt(), track.sigma1Pt() * track.pt(), track.tpcFractionSharedCls(), collision.centFT0A(), collision.centFT0C());
    histos.fill(HIST("TPC/tpcChi2NCl"), track.pt(), track.sigma1Pt() * track.pt(), track.tpcChi2NCl(), collision.centFT0A(), collision.centFT0C());
  }

  Preslice<aod::Track> trackPerColl = aod::track::collisionId;
  using CollisionCandidate = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs>;
  using TrackCandidates = soa::Join<aod::FullTracks, aod::TracksDCA, aod::TrackSelection, aod::TracksCov>;

  void processFull(CollisionCandidate const& collisions,
                   soa::Filtered<TrackCandidates> const& tracks)
  {
    for (const auto& collision : collisions) {
      auto tracksInCollision = tracks.sliceBy(trackPerColl, collision.globalIndex());
      if (checkEventSelection(collision)) {
        histos.fill(HIST("EventProp/MultCorrelations"), collision.centFT0A(), collision.centFT0C(), collision.centFT0M(), collision.multFT0A(), collision.multFT0C(), collision.multFT0M(), collision.multNTracksPV(), tracksInCollision.size());
        for (const auto& track : tracksInCollision) {
          if (track.has_collision() && (collision.globalIndex() == track.collisionId())) {
            if (checkTrackSelection(track)) {
              histos.fill(HIST("TrackEventPar/MultCorrelations"), track.pt(), track.sigma1Pt() * track.pt(), collision.centFT0A(), collision.centFT0C(), collision.centFT0M(), collision.multFT0A(), collision.multFT0C(), collision.multFT0M(), collision.multNTracksPV(), tracksInCollision.size());
              fillTrackQa(track, collision);
            }
          }
        }
      } else {
        histos.fill(HIST("EventProp/rejectedCollId"), 1);
      }
    }
  }
  PROCESS_SWITCH(TrackJetQa, processFull, "Standard data processor", true);

  void processDerived(aod::JeColls const& collisions, aod::JeTracks const& tracks)
  {
    for (const auto& collision : collisions) {
      if (checkEventSelection(collision)) {
        histos.fill(HIST("EventProp/MultCorrelations"), collision.centFT0A(), collision.centFT0C(), collision.centFT0A() + collision.centFT0C(), collision.multFT0A(), collision.multFT0C(), collision.multFT0A() + collision.multFT0C(), collision.multNTracksPV(), collision.multTracks());
        for (const auto& track : tracks) {
          if (!(track.collisionId() == collision.globalIdx())) {
            continue;
          }
          // LOGF(info, "Compatible Id's: %d (tracks) and %d (collisions)", track.collisionId(), collision.globalIdx());
          if (checkTrackSelection(track)) {
            histos.fill(HIST("TrackEventPar/MultCorrelations"), track.pt(), track.sigma1Pt() * track.pt(), collision.centFT0A(), collision.centFT0C(), collision.centFT0A() + collision.centFT0C(), collision.multFT0A(), collision.multFT0C(), collision.multFT0A() + collision.multFT0C(), collision.multNTracksPV(), collision.multTracks());
            fillTrackQa(track, collision);
          }
        }
      } else {
        histos.fill(HIST("EventProp/rejectedCollId"), 1);
      }
    }
  }
  PROCESS_SWITCH(TrackJetQa, processDerived, "Derived data processor", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow;
  workflow.push_back(adaptAnalysisTask<TrackJetQa>(cfgc));
  return workflow;
}
