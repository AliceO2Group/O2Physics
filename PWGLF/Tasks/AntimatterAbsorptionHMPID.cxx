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
///
/// \author Alberto Caliva (alberto.caliva@cern.ch)
/// \since June 27, 2023

#include <cmath>
#include <TMath.h>
#include <TPDGCode.h>
#include <TObjArray.h>
#include <TLorentzVector.h>
#include <TDatabasePDG.h>

#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsBase/Propagator.h"

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/DataTypes.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Core/PID/PIDTOF.h"
#include "Common/TableProducer/PID/pidTOFBase.h"

#include "ReconstructionDataFormats/Track.h"
#include "ReconstructionDataFormats/PID.h"
#include "ReconstructionDataFormats/TrackParametrization.h"
#include "ReconstructionDataFormats/TrackParametrizationWithError.h"
#include "ReconstructionDataFormats/DCA.h"

#include "PWGLF/DataModel/LFParticleIdentification.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;

struct AntimatterAbsorptionHMPID {

  // Registry Data
  HistogramRegistry pos_reg{"positive_tracks", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry neg_reg{"negative_tracks", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry pion_pos_reg{"pion_pos", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry kaon_pos_reg{"kaon_pos", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry proton_reg{"proton", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry deuteron_reg{"deuteron", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry pion_neg_reg{"pion_neg", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry kaon_neg_reg{"kaon_neg", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry antiproton_reg{"antiproton", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry antideuteron_reg{"antideuteron", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  // Registry MC
  HistogramRegistry pion_plus_MC_reg{"pion_plus_MC_reg", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry pion_minus_MC_reg{"pion_minus_MC_reg", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry kaon_plus_MC_reg{"kaon_plus_MC_reg", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry kaon_minus_MC_reg{"kaon_minus_MC_reg", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry proton_MC_reg{"proton_MC_reg", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry antiproton_MC_reg{"antiproton_MC_reg", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry deuteron_MC_reg{"deuteron_MC_reg", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry antideuteron_MC_reg{"antideuteron_MC_reg", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  void init(o2::framework::InitContext&)
  {

    std::vector<double> momentum_bin_edges = {0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 1.65, 1.7, 1.75, 1.8, 1.85, 1.9, 1.95, 2.0, 2.05, 2.1, 2.15, 2.2, 2.25, 2.3, 2.35, 2.4, 2.45, 2.5, 2.55, 2.6, 2.65, 2.7, 2.75, 2.8, 2.85, 2.9, 2.95, 3.0};
    AxisSpec pAxis = {momentum_bin_edges, "#it{p} (GeV/#it{c})"};

    // General Histogram (Positive Tracks)
    pos_reg.add("histTpcSignalData", "dE/dx", HistType::kTH2F,
                {{500, 0.0, 5.0, "#it{p} (GeV/#it{c})"},
                 {1400, 0, 1400, "d#it{E}/d#it{x} (a. u.)"}});
    pos_reg.add("histRecVtxZData", "collision z position", HistType::kTH1F, {{200, -20.0, +20.0, "z_{vtx} (cm)"}});
    pos_reg.add("histTofSignalData", "TOF signal", HistType::kTH2F,
                {{500, 0.0, 5.0, "#it{p} (GeV/#it{c})"},
                 {550, 0.0, 1.1, "#beta (TOF)"}});
    pos_reg.add("histDcaxyVsPData", "dca_xy vs p", HistType::kTH2F,
                {pAxis, {250, -0.5, 0.5, "DCA_{xy} (cm)"}});
    pos_reg.add("histDcaZVsPtData", "dca_z vs p", HistType::kTH2F,
                {pAxis, {1000, -2.0, 2.0, "DCA_{z} (cm)"}});
    pos_reg.add("histNClusterTPC", "Number of TPC Clusters vs p",
                HistType::kTH2F, {pAxis, {160, 0.0, 160.0, "nClustersTPC"}});
    pos_reg.add("histNCrossedRowTPC", "Number of TPC crossed row vs p",
                HistType::kTH2F, {pAxis, {160, 0.0, 160.0, "nCrossedRowsTPC"}});
    pos_reg.add("histNClusterITS", "Number of ITS Clusters vs p",
                HistType::kTH2F, {pAxis, {10, 0.0, 10.0, "nClustersITS"}});
    pos_reg.add("histChi2TPC", "chi^2 TPC vs p", HistType::kTH2F,
                {pAxis, {100, 0.0, 5.0, "chi^2_{TPC}"}});
    pos_reg.add("histChi2ITS", "chi^2 ITS vs p", HistType::kTH2F,
                {pAxis, {500, 0.0, 50.0, "chi^2_{ITS}"}});
    pos_reg.add("histEta", "eta", HistType::kTH2F,
                {pAxis, {100, -1.0, 1.0, "eta"}});
    pos_reg.add("histPhi", "phi", HistType::kTH2F,
                {pAxis, {200, 0.0, TMath::TwoPi(), "phi"}});

    // General Histogram (Negative Tracks)
    neg_reg.add("histTpcSignalData", "dE/dx", HistType::kTH2F,
                {{500, 0.0, 5.0, "#it{p} (GeV/#it{c})"},
                 {1400, 0, 1400, "d#it{E}/d#it{x} (a. u.)"}});
    neg_reg.add("histTofSignalData", "TOF signal", HistType::kTH2F,
                {{500, 0.0, 5.0, "#it{p} (GeV/#it{c})"},
                 {550, 0.0, 1.1, "#beta (TOF)"}});
    neg_reg.add("histDcaxyVsPData", "dca_xy vs p", HistType::kTH2F,
                {pAxis, {250, -0.5, 0.5, "DCA_{xy} (cm)"}});
    neg_reg.add("histDcaZVsPtData", "dca_z vs p", HistType::kTH2F,
                {pAxis, {1000, -2.0, 2.0, "DCA_{z} (cm)"}});
    neg_reg.add("histNClusterTPC", "Number of TPC Clusters vs p",
                HistType::kTH2F, {pAxis, {160, 0.0, 160.0, "nClustersTPC"}});
    neg_reg.add("histNCrossedRowTPC", "Number of TPC crossed row vs p",
                HistType::kTH2F, {pAxis, {160, 0.0, 160.0, "nCrossedRowsTPC"}});
    neg_reg.add("histNClusterITS", "Number of ITS Clusters vs p",
                HistType::kTH2F, {pAxis, {10, 0.0, 10.0, "nClustersITS"}});
    neg_reg.add("histChi2TPC", "chi^2 TPC vs p", HistType::kTH2F,
                {pAxis, {100, 0.0, 5.0, "chi^2_{TPC}"}});
    neg_reg.add("histChi2ITS", "chi^2 ITS vs p", HistType::kTH2F,
                {pAxis, {500, 0.0, 50.0, "chi^2_{ITS}"}});
    neg_reg.add("histEta", "eta", HistType::kTH2F,
                {pAxis, {100, -1.0, 1.0, "eta"}});
    neg_reg.add("histPhi", "phi", HistType::kTH2F,
                {pAxis, {200, 0.0, TMath::TwoPi(), "phi"}});

    // Pion Pos
    pion_pos_reg.add("histTpcNsigmaData", "nsigmaTPC (#pi)", HistType::kTH2F,
                     {pAxis, {160, -20.0, +20.0, "n#sigma_{TPC} (#pi)"}});
    pion_pos_reg.add("histTofNsigmaData", "nsigmaTOF (#pi)", HistType::kTH2F,
                     {pAxis, {160, -20.0, +20.0, "n#sigma_{TOF} (#pi)"}});

    // Kaon Pos
    kaon_pos_reg.add("histTpcNsigmaData", "nsigmaTPC (K)", HistType::kTH2F,
                     {pAxis, {160, -20.0, +20.0, "n#sigma_{TPC} (K)"}});
    kaon_pos_reg.add("histTofNsigmaData", "nsigmaTOF (K)", HistType::kTH2F,
                     {pAxis, {160, -20.0, +20.0, "n#sigma_{TOF} (K)"}});

    // Protons
    proton_reg.add("histTpcNsigmaData", "nsigmaTPC (p)", HistType::kTH2F,
                   {pAxis, {160, -20.0, +20.0, "n#sigma_{TPC} (p)"}});
    proton_reg.add("histTofNsigmaData", "nsigmaTOF (p)", HistType::kTH2F,
                   {pAxis, {160, -20.0, +20.0, "n#sigma_{TOF} (p)"}});

    // Deuterons
    deuteron_reg.add("histTpcNsigmaData", "nsigmaTPC (d)", HistType::kTH2F,
                     {pAxis, {160, -20.0, +20.0, "n#sigma_{TPC} (d)"}});
    deuteron_reg.add("histTofNsigmaData", "nsigmaTOF (d)", HistType::kTH2F,
                     {pAxis, {160, -20.0, +20.0, "n#sigma_{TOF} (d)"}});

    // Pion Neg
    pion_neg_reg.add("histTpcNsigmaData", "nsigmaTPC (#pi)", HistType::kTH2F,
                     {pAxis, {160, -20.0, +20.0, "n#sigma_{TPC} (#pi)"}});
    pion_neg_reg.add("histTofNsigmaData", "nsigmaTOF (#pi)", HistType::kTH2F,
                     {pAxis, {160, -20.0, +20.0, "n#sigma_{TOF} (#pi)"}});

    // Kaon Neg
    kaon_neg_reg.add("histTpcNsigmaData", "nsigmaTPC (K)", HistType::kTH2F,
                     {pAxis, {160, -20.0, +20.0, "n#sigma_{TPC} (K)"}});
    kaon_neg_reg.add("histTofNsigmaData", "nsigmaTOF (K)", HistType::kTH2F,
                     {pAxis, {160, -20.0, +20.0, "n#sigma_{TOF} (K)"}});

    // Antiprotons
    antiproton_reg.add("histTpcNsigmaData", "nsigmaTPC (antip)",
                       HistType::kTH2F,
                       {pAxis, {160, -20.0, +20.0, "n#sigma_{TPC} (p)"}});
    antiproton_reg.add("histTofNsigmaData", "nsigmaTOF (antip)",
                       HistType::kTH2F,
                       {pAxis, {160, -20.0, +20.0, "n#sigma_{TOF} (p)"}});

    // Antideuterons
    antideuteron_reg.add("histTpcNsigmaData", "nsigmaTPC (d)", HistType::kTH2F,
                         {pAxis, {160, -20.0, +20.0, "n#sigma_{TPC} (d)"}});
    antideuteron_reg.add("histTofNsigmaData", "nsigmaTOF (d)", HistType::kTH2F,
                         {pAxis, {160, -20.0, +20.0, "n#sigma_{TOF} (d)"}});

    // MC Energy loss correction maps
    pion_plus_MC_reg.add("energy_loss_corr", "energy_loss_corr_piplus", HistType::kTH2F, {{300, 0.0, +3.0, "#it{p}_{vtx} (GeV/#it{c})"}, {300, 0.0, +3.0, "#it{p}_{absorber} (GeV/#it{c})"}});
    pion_minus_MC_reg.add("energy_loss_corr", "energy_loss_corr_piminus", HistType::kTH2F, {{300, 0.0, +3.0, "#it{p}_{vtx} (GeV/#it{c})"}, {300, 0.0, +3.0, "#it{p}_{absorber} (GeV/#it{c})"}});
    kaon_plus_MC_reg.add("energy_loss_corr", "energy_loss_corr_kaonplus", HistType::kTH2F, {{300, 0.0, +3.0, "#it{p}_{vtx} (GeV/#it{c})"}, {300, 0.0, +3.0, "#it{p}_{absorber} (GeV/#it{c})"}});
    kaon_minus_MC_reg.add("energy_loss_corr", "energy_loss_corr_kaonminus", HistType::kTH2F, {{300, 0.0, +3.0, "#it{p}_{vtx} (GeV/#it{c})"}, {300, 0.0, +3.0, "#it{p}_{absorber} (GeV/#it{c})"}});
    proton_MC_reg.add("energy_loss_corr", "energy_loss_corr_proton", HistType::kTH2F, {{300, 0.0, +3.0, "#it{p}_{vtx} (GeV/#it{c})"}, {300, 0.0, +3.0, "#it{p}_{absorber} (GeV/#it{c})"}});
    antiproton_MC_reg.add("energy_loss_corr", "energy_loss_antiproton", HistType::kTH2F, {{300, 0.0, +3.0, "#it{p}_{vtx} (GeV/#it{c})"}, {300, 0.0, +3.0, "#it{p}_{absorber} (GeV/#it{c})"}});
    deuteron_MC_reg.add("energy_loss_corr", "energy_loss_corr_deuteron", HistType::kTH2F, {{300, 0.0, +3.0, "#it{p}_{vtx} (GeV/#it{c})"}, {300, 0.0, +3.0, "#it{p}_{absorber} (GeV/#it{c})"}});
    antideuteron_MC_reg.add("energy_loss_corr", "energy_loss_antideuteron", HistType::kTH2F, {{300, 0.0, +3.0, "#it{p}_{vtx} (GeV/#it{c})"}, {300, 0.0, +3.0, "#it{p}_{absorber} (GeV/#it{c})"}});
  }

  // Configurables
  Configurable<float> zVertexRange{"zVertexRange", 10.0f, "zVertexRange"};
  Configurable<float> pmin{"pmin", 0.1, "pmin"};
  Configurable<float> pmax{"pmax", 2.0, "pmax"};
  Configurable<float> etaMin{"etaMin", -0.8, "etaMin"};
  Configurable<float> etaMax{"etaMax", +0.8, "etaMax"};
  Configurable<float> phiMin{"phiMin", 0.0, "phiMin"};
  Configurable<float> phiMax{"phiMax", TMath::TwoPi(), "phiMax"};
  Configurable<float> nsigmaTPCMin{"nsigmaTPCMin", -3.0, "nsigmaTPCMin"};
  Configurable<float> nsigmaTPCMax{"nsigmaTPCMax", +3.0, "nsigmaTPCMax"};
  Configurable<float> nsigmaTOFMin{"nsigmaTOFMin", -3.0, "nsigmaTOFMin"};
  Configurable<float> nsigmaTOFMax{"nsigmaTOFMax", +3.5, "nsigmaTOFMax"};
  Configurable<float> minReqClusterITS{"minReqClusterITS", 1.0, "min number of clusters required in ITS"};
  Configurable<float> minTPCnClsFound{"minTPCnClsFound", 0.0f, "minTPCnClsFound"};
  Configurable<float> minNCrossedRowsTPC{"minNCrossedRowsTPC", 70.0f, "min number of crossed rows TPC"};
  Configurable<float> maxChi2ITS{"maxChi2ITS", 36.0f, "max chi2 per cluster ITS"};
  Configurable<float> maxChi2TPC{"maxChi2TPC", 4.0f, "max chi2 per cluster TPC"};
  Configurable<float> maxDCA_xy{"maxDCA_xy", 0.5f, "maxDCA_xy"};
  Configurable<float> maxDCA_z{"maxDCA_z", 2.0f, "maxDCA_z"};
  Configurable<int> lastRequiredTrdCluster{"lastRequiredTrdCluster", 5, "Last TRD cluster (-1 = no requirement)"};
  Configurable<bool> enable_PVcontributor_global{"enable_PVcontributor_global", true, "is PV contributor (global)"};
  Configurable<bool> enable_PVcontributor_proton{"enable_PVcontributor_proton", true, "is PV contributor (proton)"};
  Configurable<bool> enable_PVcontributor_antiproton{"enable_PVcontributor_antiproton", true, "is PV contributor (antiproton)"};
  Configurable<bool> enable_PVcontributor_deuteron{"enable_PVcontributor_deuteron", true, "is PV contributor (deuteron)"};
  Configurable<bool> enable_PVcontributor_antideuteron{"enable_PVcontributor_antideuteron", true, "is PV contributor (antideuteron)"};

  template <typename CollisionType, typename TracksType>
  void fillHistograms(const CollisionType& event, const TracksType& tracks)
  {
    // Event Selection
    if (!event.sel8())
      return;

    // Event Counter
    pos_reg.fill(HIST("histRecVtxZData"), event.posZ());

    // Loop over Reconstructed Tracks
    for (auto track : tracks) {

      // Loose Track Selection
      if (!track.isGlobalTrackWoDCA())
        continue;
      if (!track.passedITSRefit())
        continue;
      if (!track.passedTPCRefit())
        continue;
      if (track.itsNCls() < 1)
        continue;
      if (track.tpcNClsFound() < 0)
        continue;
      if (track.tpcNClsCrossedRows() < 60)
        continue;
      if (TMath::Abs(track.dcaXY()) > 1.0)
        continue;
      if (TMath::Abs(track.dcaZ()) > 1.0)
        continue;
      if (enable_PVcontributor_global && !(track.isPVContributor()))
        continue;

      // Fill QA Histograms (Positive Tracks)
      if (track.sign() > 0) {

        pos_reg.fill(HIST("histTpcSignalData"), track.tpcInnerParam(), track.tpcSignal());
        pos_reg.fill(HIST("histTofSignalData"), track.p(), track.beta());
        pos_reg.fill(HIST("histDcaxyVsPData"), track.p(), track.dcaXY());
        pos_reg.fill(HIST("histDcaZVsPtData"), track.p(), track.dcaZ());
        pos_reg.fill(HIST("histNClusterTPC"), track.p(), track.tpcNClsFound());
        pos_reg.fill(HIST("histNCrossedRowTPC"), track.p(), track.tpcNClsCrossedRows());
        pos_reg.fill(HIST("histNClusterITS"), track.p(), track.itsNCls());
        pos_reg.fill(HIST("histChi2TPC"), track.p(), track.tpcChi2NCl());
        pos_reg.fill(HIST("histChi2ITS"), track.p(), track.itsChi2NCl());
        pos_reg.fill(HIST("histEta"), track.p(), track.eta());
        pos_reg.fill(HIST("histPhi"), track.p(), track.phi());
      }

      // Fill QA Histograms (Negative Tracks)
      if (track.sign() < 0) {

        neg_reg.fill(HIST("histTpcSignalData"), track.tpcInnerParam(), track.tpcSignal());
        neg_reg.fill(HIST("histTofSignalData"), track.p(), track.beta());
        neg_reg.fill(HIST("histDcaxyVsPData"), track.p(), track.dcaXY());
        neg_reg.fill(HIST("histDcaZVsPtData"), track.p(), track.dcaZ());
        neg_reg.fill(HIST("histNClusterTPC"), track.p(), track.tpcNClsFound());
        neg_reg.fill(HIST("histNCrossedRowTPC"), track.p(), track.tpcNClsCrossedRows());
        neg_reg.fill(HIST("histNClusterITS"), track.p(), track.itsNCls());
        neg_reg.fill(HIST("histChi2TPC"), track.p(), track.tpcChi2NCl());
        neg_reg.fill(HIST("histChi2ITS"), track.p(), track.itsChi2NCl());
        neg_reg.fill(HIST("histEta"), track.p(), track.eta());
        neg_reg.fill(HIST("histPhi"), track.p(), track.phi());
      }

      // Track Selection
      if (!track.passedITSRefit())
        continue;
      if (!track.passedTPCRefit())
        continue;
      if (track.itsNCls() < minReqClusterITS)
        continue;
      if (track.tpcNClsFound() < minTPCnClsFound)
        continue;
      if (track.tpcNClsCrossedRows() < minNCrossedRowsTPC)
        continue;
      if (track.tpcChi2NCl() > maxChi2TPC)
        continue;
      if (track.itsChi2NCl() > maxChi2ITS)
        continue;
      if (TMath::Abs(track.dcaXY()) > maxDCA_xy)
        continue;
      if (TMath::Abs(track.dcaZ()) > maxDCA_z)
        continue;
      if (track.eta() < etaMin)
        continue;
      if (track.eta() > etaMax)
        continue;
      if (track.phi() < phiMin)
        continue;
      if (track.phi() > phiMax)
        continue;

      if (track.sign() > 0) {
        pion_pos_reg.fill(HIST("histTpcNsigmaData"), track.p(), track.tpcNSigmaPi());
        kaon_pos_reg.fill(HIST("histTpcNsigmaData"), track.p(), track.tpcNSigmaKa());
        proton_reg.fill(HIST("histTpcNsigmaData"), track.p(), track.tpcNSigmaPr());
        deuteron_reg.fill(HIST("histTpcNsigmaData"), track.p(), track.tpcNSigmaDe());
      }

      if (track.sign() < 0) {
        pion_neg_reg.fill(HIST("histTpcNsigmaData"), track.p(), track.tpcNSigmaPi());
        kaon_neg_reg.fill(HIST("histTpcNsigmaData"), track.p(), track.tpcNSigmaKa());
        antiproton_reg.fill(HIST("histTpcNsigmaData"), track.p(), track.tpcNSigmaPr());
        antideuteron_reg.fill(HIST("histTpcNsigmaData"), track.p(), track.tpcNSigmaDe());
      }

      bool passedPionTPCsel = false;
      bool passedKaonTPCsel = false;
      bool passedProtTPCsel = false;
      bool passedDeutTPCsel = false;

      if (track.tpcNSigmaPi() > nsigmaTPCMin &&
          track.tpcNSigmaPi() < nsigmaTPCMax)
        passedPionTPCsel = true;

      if (track.tpcNSigmaKa() > nsigmaTPCMin &&
          track.tpcNSigmaKa() < nsigmaTPCMax)
        passedKaonTPCsel = true;

      if (track.tpcNSigmaPr() > nsigmaTPCMin &&
          track.tpcNSigmaPr() < nsigmaTPCMax)
        passedProtTPCsel = true;

      if (track.tpcNSigmaDe() > nsigmaTPCMin &&
          track.tpcNSigmaDe() < nsigmaTPCMax)
        passedDeutTPCsel = true;

      if (!track.hasTOF())
        continue;

      if (track.sign() > 0) {
        if (passedPionTPCsel)
          pion_pos_reg.fill(HIST("histTofNsigmaData"), track.p(), track.tofNSigmaPi());
        if (passedKaonTPCsel)
          kaon_pos_reg.fill(HIST("histTofNsigmaData"), track.p(), track.tofNSigmaKa());
        if (passedProtTPCsel)
          proton_reg.fill(HIST("histTofNsigmaData"), track.p(), track.tofNSigmaPr());
        if (passedDeutTPCsel)
          deuteron_reg.fill(HIST("histTofNsigmaData"), track.p(), track.tofNSigmaDe());
      }

      if (track.sign() < 0) {
        if (passedPionTPCsel)
          pion_neg_reg.fill(HIST("histTofNsigmaData"), track.p(), track.tofNSigmaPi());
        if (passedKaonTPCsel)
          kaon_neg_reg.fill(HIST("histTofNsigmaData"), track.p(), track.tofNSigmaKa());
        if (passedProtTPCsel)
          antiproton_reg.fill(HIST("histTofNsigmaData"), track.p(), track.tofNSigmaPr());
        if (passedDeutTPCsel)
          antideuteron_reg.fill(HIST("histTofNsigmaData"), track.p(), track.tofNSigmaDe());
      }
    }
  }

  // Filter collisionFilter = (nabs(aod::collision::posZ) < zVertexRange);
  // Filter trackFilter = (nabs(aod::track::eta) < 0.8f && requireGlobalTrackWoDCAInFilter());

  using EventCandidates = soa::Join<aod::Collisions, aod::EvSels>;
  // using EventCandidates = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>;

  /*
  using TrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TracksCov, aod::TracksIU, aod::TracksCovIU, aod::TracksExtra, aod::TracksDCA, aod::McTrackLabels, aod::TrackSelection, aod::TrackSelectionExtension, aod::pidTPCLfFullPi, aod::pidTOFFullPi, aod::pidTPCLfFullKa, aod::pidTOFFullKa, aod::pidTPCLfFullPr, aod::pidTOFFullPr, aod::pidTPCLfFullDe, aod::pidTOFFullDe, aod::pidTPCLfFullTr, aod::pidTOFFullTr, aod::pidTPCLfFullHe, aod::pidTOFFullHe, aod::pidTPCLfFullAl, aod::pidTOFFullAl, aod::TrackSelection, aod::TrackSelectionExtension, aod::TOFSignal, aod::pidTOFmass, aod::pidTOFbeta>>;
  */
  // Info for TPC PID
  using PidInfoTPC = soa::Join<aod::pidTPCLfFullPi, aod::pidTPCLfFullKa,
                               aod::pidTPCLfFullPr, aod::pidTPCLfFullDe,
                               aod::pidTPCLfFullTr, aod::pidTPCLfFullHe,
                               aod::pidTPCLfFullAl>;

  // Info for TOF PID
  using PidInfoTOF = soa::Join<aod::pidTOFFullPi, aod::pidTOFFullKa,
                               aod::pidTOFFullPr, aod::pidTOFFullDe,
                               aod::pidTOFFullTr, aod::pidTOFFullHe,
                               aod::pidTOFFullAl,
                               aod::TOFSignal, aod::pidTOFmass, aod::pidTOFbeta>;

  // Propagated tracks
  using TrackCandidatesIU = soa::Join<aod::TracksIU, aod::TracksCovIU, aod::TracksExtra, aod::TracksDCA,
                                      PidInfoTPC, PidInfoTOF,
                                      aod::TrackSelection, aod::TrackSelectionExtension>;

  // Propagated to PV tracks
  using TrackCandidates = soa::Join<aod::TracksIU, aod::TracksCovIU, aod::TracksExtra, aod::TracksDCA,
                                    PidInfoTPC, PidInfoTOF,
                                    aod::TrackSelection, aod::TrackSelectionExtension>;

  // using TrackCandidates = soa::Join<aod::Tracks, aod::TracksCov, aod::TracksIU, aod::TracksCovIU, aod::TracksExtra, aod::TracksDCA, aod::pidTPCLfFullPi, aod::pidTOFFullPi, aod::pidTPCLfFullKa, aod::pidTOFFullKa, aod::pidTPCLfFullPr, aod::pidTOFFullPr, aod::pidTPCLfFullDe, aod::pidTOFFullDe, aod::pidTPCLfFullTr, aod::pidTOFFullTr, aod::pidTPCLfFullHe, aod::pidTOFFullHe, aod::pidTPCLfFullAl, aod::pidTOFFullAl, aod::TrackSelection, aod::TrackSelectionExtension, aod::TOFSignal, aod::pidTOFmass, aod::pidTOFbeta>;

  void processData(EventCandidates::iterator const& event,
                   TrackCandidatesIU const& tracksIU)
  //  , TrackCandidates const& tracks)
  {
    fillHistograms(event, tracksIU);
  }
  PROCESS_SWITCH(AntimatterAbsorptionHMPID, processData, "process data", true);

  // Process MC
  void processMC(soa::Join<EventCandidates, aod::McCollisionLabels>::iterator const& event,
                 soa::Join<TrackCandidatesIU, aod::McTrackLabels> const& tracksIU,
                 //  soa::Join<TrackCandidates, aod::McTrackLabels> const& tracks,
                 aod::McParticles& mcParticles,
                 aod::McCollisions const& mcCollisions)
  {

    float radius_hmpid = 500;

    for (auto track : tracksIU) {

      // Get MC Particle
      if (!track.has_mcParticle())
        continue;
      const auto particle = track.mcParticle();

      // Track Selection
      if (!track.isGlobalTrackWoDCA())
        continue;
      if (!track.passedITSRefit())
        continue;
      if (!track.passedTPCRefit())
        continue;
      if (track.itsNCls() < minReqClusterITS)
        continue;
      if (track.tpcNClsFound() < minTPCnClsFound)
        continue;
      if (track.tpcNClsCrossedRows() < minNCrossedRowsTPC)
        continue;
      if (track.tpcChi2NCl() > maxChi2TPC)
        continue;
      if (track.itsChi2NCl() > maxChi2ITS)
        continue;
      if (TMath::Abs(track.dcaXY()) > maxDCA_xy)
        continue;
      if (TMath::Abs(track.dcaZ()) > maxDCA_z)
        continue;
      if (enable_PVcontributor_global && !(track.isPVContributor()))
        continue;
      if (track.eta() < etaMin)
        continue;
      if (track.eta() > etaMax)
        continue;
      if (track.phi() < phiMin)
        continue;
      if (track.phi() > phiMax)
        continue;

      // TOF Hit Requirement
      if (!track.hasTOF())
        continue;

      // Pi Plus
      if (particle.pdgCode() == 211) {
        if (track.tpcNSigmaPi() < nsigmaTPCMin)
          continue;
        if (track.tpcNSigmaPi() > nsigmaTPCMax)
          continue;
        if (track.tofNSigmaPi() < nsigmaTOFMin)
          continue;
        if (track.tofNSigmaPi() > nsigmaTOFMax)
          continue;

        double p_vtx = track.p();
        auto par_cov = getTrackParCov(track);
        par_cov.setPID(2); // Pion
        auto prop = o2::base::Propagator::Instance();
        float xtogo = 0;
        if (!par_cov.getXatLabR(radius_hmpid, xtogo, prop->getNominalBz(), o2::track::DirType::DirOutward) || !prop->PropagateToXBxByBz(par_cov, xtogo, 0.95, 10, o2::base::Propagator::MatCorrType::USEMatCorrLUT)) {
          continue;
        }
        double p_hmpid = track.p();
        pion_plus_MC_reg.fill(HIST("energy_loss_corr"), p_vtx, p_hmpid);
      }

      // Pi Minus
      if (particle.pdgCode() == -211) {
        if (track.tpcNSigmaPi() < nsigmaTPCMin)
          continue;
        if (track.tpcNSigmaPi() > nsigmaTPCMax)
          continue;
        if (track.tofNSigmaPi() < nsigmaTOFMin)
          continue;
        if (track.tofNSigmaPi() > nsigmaTOFMax)
          continue;

        double p_vtx = track.p();
        auto par_cov = getTrackParCov(track);
        par_cov.setPID(2); // Pion
        auto prop = o2::base::Propagator::Instance();
        float xtogo = 0;
        if (!par_cov.getXatLabR(radius_hmpid, xtogo, prop->getNominalBz(), o2::track::DirType::DirOutward) || !prop->PropagateToXBxByBz(par_cov, xtogo, 0.95, 10, o2::base::Propagator::MatCorrType::USEMatCorrLUT)) {
          continue;
        }
        double p_hmpid = track.p();
        pion_minus_MC_reg.fill(HIST("energy_loss_corr"), p_vtx, p_hmpid);
      }

      // Kaon Plus
      if (particle.pdgCode() == 321) {
        if (track.tpcNSigmaKa() < nsigmaTPCMin)
          continue;
        if (track.tpcNSigmaKa() > nsigmaTPCMax)
          continue;
        if (track.tofNSigmaKa() < nsigmaTOFMin)
          continue;
        if (track.tofNSigmaKa() > nsigmaTOFMax)
          continue;

        double p_vtx = track.p();
        auto par_cov = getTrackParCov(track);
        par_cov.setPID(3); // Kaon
        auto prop = o2::base::Propagator::Instance();
        float xtogo = 0;
        if (!par_cov.getXatLabR(radius_hmpid, xtogo, prop->getNominalBz(), o2::track::DirType::DirOutward) || !prop->PropagateToXBxByBz(par_cov, xtogo, 0.95, 10, o2::base::Propagator::MatCorrType::USEMatCorrLUT)) {
          continue;
        }
        double p_hmpid = track.p();
        kaon_plus_MC_reg.fill(HIST("energy_loss_corr"), p_vtx, p_hmpid);
      }

      // Kaon Minus
      if (particle.pdgCode() == -321) {
        if (track.tpcNSigmaKa() < nsigmaTPCMin)
          continue;
        if (track.tpcNSigmaKa() > nsigmaTPCMax)
          continue;
        if (track.tofNSigmaKa() < nsigmaTOFMin)
          continue;
        if (track.tofNSigmaKa() > nsigmaTOFMax)
          continue;

        double p_vtx = track.p();
        auto par_cov = getTrackParCov(track);
        par_cov.setPID(3); // Kaon
        auto prop = o2::base::Propagator::Instance();
        float xtogo = 0;
        if (!par_cov.getXatLabR(radius_hmpid, xtogo, prop->getNominalBz(), o2::track::DirType::DirOutward) || !prop->PropagateToXBxByBz(par_cov, xtogo, 0.95, 10, o2::base::Propagator::MatCorrType::USEMatCorrLUT)) {
          continue;
        }
        double p_hmpid = track.p();
        kaon_minus_MC_reg.fill(HIST("energy_loss_corr"), p_vtx, p_hmpid);
      }

      // Proton
      if (particle.pdgCode() == 2212) {
        if (track.tpcNSigmaPr() < nsigmaTPCMin)
          continue;
        if (track.tpcNSigmaPr() > nsigmaTPCMax)
          continue;
        if (track.tofNSigmaPr() < nsigmaTOFMin)
          continue;
        if (track.tofNSigmaPr() > nsigmaTOFMax)
          continue;

        double p_vtx = track.p();
        auto par_cov = getTrackParCov(track);
        par_cov.setPID(4); // Proton
        auto prop = o2::base::Propagator::Instance();
        float xtogo = 0;
        if (!par_cov.getXatLabR(radius_hmpid, xtogo, prop->getNominalBz(), o2::track::DirType::DirOutward) || !prop->PropagateToXBxByBz(par_cov, xtogo, 0.95, 10, o2::base::Propagator::MatCorrType::USEMatCorrLUT)) {
          continue;
        }
        double p_hmpid = track.p();
        proton_MC_reg.fill(HIST("energy_loss_corr"), p_vtx, p_hmpid);
      }

      // AntiProton
      if (particle.pdgCode() == -2212) {
        if (track.tpcNSigmaPr() < nsigmaTPCMin)
          continue;
        if (track.tpcNSigmaPr() > nsigmaTPCMax)
          continue;
        if (track.tofNSigmaPr() < nsigmaTOFMin)
          continue;
        if (track.tofNSigmaPr() > nsigmaTOFMax)
          continue;

        double p_vtx = track.p();
        auto par_cov = getTrackParCov(track);
        par_cov.setPID(4); // Proton
        auto prop = o2::base::Propagator::Instance();
        float xtogo = 0;
        if (!par_cov.getXatLabR(radius_hmpid, xtogo, prop->getNominalBz(), o2::track::DirType::DirOutward) || !prop->PropagateToXBxByBz(par_cov, xtogo, 0.95, 10, o2::base::Propagator::MatCorrType::USEMatCorrLUT)) {
          continue;
        }
        double p_hmpid = track.p();
        antiproton_MC_reg.fill(HIST("energy_loss_corr"), p_vtx, p_hmpid);
      }

      // Deuteron
      if (particle.pdgCode() == 1000010020) {
        if (track.tpcNSigmaDe() < nsigmaTPCMin)
          continue;
        if (track.tpcNSigmaDe() > nsigmaTPCMax)
          continue;
        if (track.tofNSigmaDe() < nsigmaTOFMin)
          continue;
        if (track.tofNSigmaDe() > nsigmaTOFMax)
          continue;

        double p_vtx = track.p();
        auto par_cov = getTrackParCov(track);
        par_cov.setPID(5); // Deuteron
        auto prop = o2::base::Propagator::Instance();
        float xtogo = 0;
        if (!par_cov.getXatLabR(radius_hmpid, xtogo, prop->getNominalBz(), o2::track::DirType::DirOutward) || !prop->PropagateToXBxByBz(par_cov, xtogo, 0.95, 10, o2::base::Propagator::MatCorrType::USEMatCorrLUT)) {
          continue;
        }
        double p_hmpid = track.p();
        deuteron_MC_reg.fill(HIST("energy_loss_corr"), p_vtx, p_hmpid);
      }

      // AntiDeuteron
      if (particle.pdgCode() == -1000010020) {
        if (track.tpcNSigmaDe() < nsigmaTPCMin)
          continue;
        if (track.tpcNSigmaDe() > nsigmaTPCMax)
          continue;
        if (track.tofNSigmaDe() < nsigmaTOFMin)
          continue;
        if (track.tofNSigmaDe() > nsigmaTOFMax)
          continue;

        double p_vtx = track.p();
        auto par_cov = getTrackParCov(track);
        par_cov.setPID(5); // Deuteron
        auto prop = o2::base::Propagator::Instance();
        float xtogo = 0;
        if (!par_cov.getXatLabR(radius_hmpid, xtogo, prop->getNominalBz(), o2::track::DirType::DirOutward) || !prop->PropagateToXBxByBz(par_cov, xtogo, 0.95, 10, o2::base::Propagator::MatCorrType::USEMatCorrLUT)) {
          continue;
        }
        double p_hmpid = track.p();
        antideuteron_MC_reg.fill(HIST("energy_loss_corr"), p_vtx, p_hmpid);
      }
    }
  }
  PROCESS_SWITCH(AntimatterAbsorptionHMPID, processMC, "process MC", false);
};

//**********************************************************************************************************************************************

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<AntimatterAbsorptionHMPID>(
    cfgc)};
}
