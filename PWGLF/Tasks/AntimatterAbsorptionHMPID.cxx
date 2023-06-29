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
#include <TLorentzVector.h>
#include <TMath.h>
#include <TObjArray.h>

#include "ReconstructionDataFormats/Track.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"

#include "Framework/HistogramRegistry.h"

#include "PWGLF/DataModel/LFParticleIdentification.h"
#include "PWGDQ/DataModel/ReducedInfoTables.h"
#include "TPDGCode.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct AntimatterAbsorptionHMPID {

  // Registry (Data)
  HistogramRegistry pos_reg{
    "positive_tracks",
    {},
    OutputObjHandlingPolicy::AnalysisObject,
    true,
    true};
  HistogramRegistry neg_reg{
    "negative_tracks",
    {},
    OutputObjHandlingPolicy::AnalysisObject,
    true,
    true};
  HistogramRegistry proton_reg{
    "proton",
    {},
    OutputObjHandlingPolicy::AnalysisObject,
    true,
    true};
  HistogramRegistry deuteron_reg{
    "deuteron",
    {},
    OutputObjHandlingPolicy::AnalysisObject,
    true,
    true};
  HistogramRegistry antiproton_reg{
    "antiproton",
    {},
    OutputObjHandlingPolicy::AnalysisObject,
    true,
    true};
  HistogramRegistry antideuteron_reg{
    "antideuteron",
    {},
    OutputObjHandlingPolicy::AnalysisObject,
    true,
    true};

  void init(o2::framework::InitContext&)
  {

    std::vector<double> momentum_bin_edges = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7,
                                              0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4,
                                              1.5, 1.6, 1.7, 1.8, 1.9, 2.0};
    AxisSpec pAxis = {momentum_bin_edges, "#it{p} (GeV/#it{c})"};

    // General Histogram (Positive Tracks)
    pos_reg.add("histTpcSignalData", "dE/dx", HistType::kTH2F,
                {{500, 0.0, 5.0, "#it{p} (GeV/#it{c})"},
                 {1400, 0, 1400, "d#it{E}/d#it{x} (a. u.)"}});
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
  Configurable<float> minReqClusterITS{
    "minReqClusterITS", 1.0, "min number of clusters required in ITS"};
  Configurable<float> minTPCnClsFound{"minTPCnClsFound", 0.0f,
                                      "minTPCnClsFound"};
  Configurable<float> minNCrossedRowsTPC{"minNCrossedRowsTPC", 70.0f,
                                         "min number of crossed rows TPC"};
  Configurable<float> maxChi2ITS{"maxChi2ITS", 36.0f,
                                 "max chi2 per cluster ITS"};
  Configurable<float> maxChi2TPC{"maxChi2TPC", 4.0f,
                                 "max chi2 per cluster TPC"};
  Configurable<float> maxDCA_xy{"maxDCA_xy", 0.5f, "maxDCA_xy"};
  Configurable<float> maxDCA_z{"maxDCA_z", 2.0f, "maxDCA_z"};
  Configurable<int> lastRequiredTrdCluster{
    "lastRequiredTrdCluster", 5, "Last TRD cluster (-1 = no requirement)"};
  Configurable<bool> enable_PVcontributor_global{
    "enable_PVcontributor_global", true, "is PV contributor (global)"};
  Configurable<bool> enable_PVcontributor_proton{
    "enable_PVcontributor_proton", true, "is PV contributor (proton)"};
  Configurable<bool> enable_PVcontributor_antiproton{
    "enable_PVcontributor_antiproton", true, "is PV contributor (antiproton)"};
  Configurable<bool> enable_PVcontributor_deuteron{
    "enable_PVcontributor_deuteron", true, "is PV contributor (deuteron)"};
  Configurable<bool> enable_PVcontributor_antideuteron{
    "enable_PVcontributor_antideuteron", true, "is PV contributor (antideuteron)"};

  template <typename CollisionType, typename TracksType>
  void fillHistograms(const CollisionType& event, const TracksType& tracks)
  {

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

        pos_reg.fill(HIST("histTpcSignalData"), track.tpcInnerParam(),
                     track.tpcSignal());
        pos_reg.fill(HIST("histTofSignalData"), track.p(), track.beta());
        pos_reg.fill(HIST("histDcaxyVsPData"), track.p(), track.dcaXY());
        pos_reg.fill(HIST("histDcaZVsPtData"), track.p(), track.dcaZ());
        pos_reg.fill(HIST("histNClusterTPC"), track.p(), track.tpcNClsFound());
        pos_reg.fill(HIST("histNCrossedRowTPC"), track.p(),
                     track.tpcNClsCrossedRows());
        pos_reg.fill(HIST("histNClusterITS"), track.p(), track.itsNCls());
        pos_reg.fill(HIST("histChi2TPC"), track.p(), track.tpcChi2NCl());
        pos_reg.fill(HIST("histChi2ITS"), track.p(), track.itsChi2NCl());
        pos_reg.fill(HIST("histEta"), track.p(), track.eta());
        pos_reg.fill(HIST("histPhi"), track.p(), track.phi());
      }

      // Fill QA Histograms (Negative Tracks)
      if (track.sign() < 0) {

        neg_reg.fill(HIST("histTpcSignalData"), track.tpcInnerParam(),
                     track.tpcSignal());
        neg_reg.fill(HIST("histTofSignalData"), track.p(), track.beta());
        neg_reg.fill(HIST("histDcaxyVsPData"), track.p(), track.dcaXY());
        neg_reg.fill(HIST("histDcaZVsPtData"), track.p(), track.dcaZ());
        neg_reg.fill(HIST("histNClusterTPC"), track.p(), track.tpcNClsFound());
        neg_reg.fill(HIST("histNCrossedRowTPC"), track.p(),
                     track.tpcNClsCrossedRows());
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
        proton_reg.fill(HIST("histTpcNsigmaData"), track.p(),
                        track.tpcNSigmaPr());
        deuteron_reg.fill(HIST("histTpcNsigmaData"), track.p(),
                          track.tpcNSigmaDe());
      }

      if (track.sign() < 0) {
        antiproton_reg.fill(HIST("histTpcNsigmaData"), track.p(),
                            track.tpcNSigmaPr());
        antideuteron_reg.fill(HIST("histTpcNsigmaData"), track.p(),
                              track.tpcNSigmaDe());
      }

      bool passedProtTPCsel = false;
      bool passedDeutTPCsel = false;

      if (track.tpcNSigmaPr() > nsigmaTPCMin &&
          track.tpcNSigmaPr() < nsigmaTPCMax)
        passedProtTPCsel = true;
      if (track.tpcNSigmaDe() > nsigmaTPCMin &&
          track.tpcNSigmaDe() < nsigmaTPCMax)
        passedDeutTPCsel = true;

      if (track.hasTOF()) {
        if (track.sign() > 0) {
          if (passedProtTPCsel)
            proton_reg.fill(HIST("histTofNsigmaData"), track.p(),
                            track.tofNSigmaPr());
          if (passedDeutTPCsel)
            deuteron_reg.fill(HIST("histTofNsigmaData"), track.p(),
                              track.tofNSigmaDe());
        }

        if (track.sign() < 0) {
          if (passedProtTPCsel)
            antiproton_reg.fill(HIST("histTofNsigmaData"), track.p(),
                                track.tofNSigmaPr());
          if (passedDeutTPCsel)
            antideuteron_reg.fill(HIST("histTofNsigmaData"), track.p(),
                                  track.tofNSigmaDe());
        }
      }
    }
  }

  Filter collisionFilter = (nabs(aod::collision::posZ) < zVertexRange);
  Filter trackFilter =
    (nabs(aod::track::eta) < 0.8f && requireGlobalTrackWoDCAInFilter());

  using EventCandidates =
    soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>;

  using TrackCandidates = soa::Filtered<soa::Join<
    aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::pidTPCLfFullPr,
    aod::pidTOFFullPr, aod::pidTPCLfFullDe, aod::pidTOFFullDe,
    aod::pidTPCLfFullTr, aod::pidTOFFullTr, aod::pidTPCLfFullHe,
    aod::pidTOFFullHe, aod::pidTPCLfFullAl, aod::pidTOFFullAl,
    aod::TrackSelection, aod::TrackSelectionExtension, aod::TOFSignal,
    aod::pidTOFmass, aod::pidTOFbeta>>;

  void processData(EventCandidates::iterator const& event,
                   TrackCandidates const& tracks)
  {
    fillHistograms(event, tracks);
  }
  PROCESS_SWITCH(AntimatterAbsorptionHMPID, processData, "process data", true);
};

//**********************************************************************************************************************************************

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<AntimatterAbsorptionHMPID>(
    cfgc)};
}
