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
// O2 includes

#include <cmath>
#include <string>

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "DataFormatsTPC/BetheBlochAleph.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

#include "PWGLF/DataModel/Vtx3BodyTables.h"
#include "../filterTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace
{

static constexpr int nNuclei{3};
static constexpr int nHyperNuclei{1};
static constexpr int nCutsPID{5};
static constexpr std::array<float, nNuclei> masses{
  constants::physics::MassDeuteron, constants::physics::MassTriton,
  constants::physics::MassHelium3};
static constexpr std::array<int, nNuclei> charges{1, 1, 2};
static const std::vector<std::string> matterOrNot{"Matter", "Antimatter"};
static const std::vector<std::string> nucleiNames{"H2", "H3", "Helium"};
static const std::vector<std::string> hypernucleiNames{"H3L"}; // 3-body decay case
static const std::vector<std::string> cutsNames{
  "TPCnSigmaMin", "TPCnSigmaMax", "TOFnSigmaMin", "TOFnSigmaMax", "TOFpidStartPt"};
constexpr double betheBlochDefault[nNuclei][6]{
  {-1.e32, -1.e32, -1.e32, -1.e32, -1.e32, -1.e32},
  {-1.e32, -1.e32, -1.e32, -1.e32, -1.e32, -1.e32},
  {-1.e32, -1.e32, -1.e32, -1.e32, -1.e32, -1.e32}};
static constexpr float cutsPID[nNuclei][nCutsPID]{
  {-3.f, +3.f, -4.f, +4.f, 1.0f},    /*H2*/
  {-3.f, +3.f, -4.f, +4.f, 1.6f},    /*H3*/
  {-5.f, +5.f, -4.f, +4.f, 14000.f}, /*He3*/
};
constexpr double bbMomScalingDefault[nNuclei][2]{
  {1., 1.},
  {1., 1.},
  {1., 1.}};
constexpr double minTPCmom[nNuclei][2]{
  {0., 0.},
  {0., 0.},
  {0.8, 0.}};
static const std::vector<std::string> betheBlochParNames{"p0", "p1", "p2", "p3", "p4", "resolution"};

std::shared_ptr<TH2> h2TPCsignal[nNuclei];
std::shared_ptr<TH2> h2TPCnSigma[nNuclei];

} // namespace

struct nucleiFilter {

  Produces<aod::NucleiFilters> tags;

  // configurable for nuclei
  Configurable<float> cfgCutVertex{"cfgCutVertex", 12.0f, "Accepted z-vertex range"};
  Configurable<float> cfgCutEta{"cfgCutEta", 1.f, "Eta range for tracks"};

  Configurable<float> cfgCutNclusITS{"cfgCutNclusITS", 2, "Minimum number of ITS clusters"};
  Configurable<float> cfgCutNclusTPC{"cfgCutNclusTPC", 80, "Minimum number of TPC clusters"};
  Configurable<float> cfgCutDCAxy{"cfgCutDCAxy", 3, "Max DCAxy"};
  Configurable<float> cfgCutDCAz{"cfgCutDCAz", 10, "Max DCAz"};

  Configurable<LabeledArray<double>> cfgBetheBlochParams{"cfgBetheBlochParams", {betheBlochDefault[0], nNuclei, 6, nucleiNames, betheBlochParNames}, "TPC Bethe-Bloch parameterisation for light nuclei"};
  Configurable<LabeledArray<double>> cfgMomentumScalingBetheBloch{"cfgMomentumScalingBetheBloch", {bbMomScalingDefault[0], nNuclei, 2, nucleiNames, matterOrNot}, "TPC Bethe-Bloch momentum scaling for light nuclei"};
  Configurable<LabeledArray<double>> cfgMinTPCmom{"cfgMinTPCmom", {minTPCmom[0], nNuclei, 2, nucleiNames, matterOrNot}, "Minimum TPC p/Z for nuclei PID"};

  Configurable<LabeledArray<float>> cfgCutsPID{"nucleiCutsPID", {cutsPID[0], nNuclei, nCutsPID, nucleiNames, cutsNames}, "Nuclei PID selections"};

  // configurable for hypertriton 3body decay
  Configurable<float> minCosPA3body{"minCosPA3body", 0.99, "minCosPA3body"};
  Configurable<float> dcavtxdau{"dcavtxdau", 1.0, "DCA Vtx Daughters"};
  Configurable<float> dcapiontopv{"dcapiontopv", 0.05, "DCA Pion To PV"};
  Configurable<float> TofPidNsigmaMin{"TofPidNsigmaMin", -5, "TofPidNsigmaMin"};
  Configurable<float> TofPidNsigmaMax{"TofPidNsigmaMax", 5, "TofPidNsigmaMax"};
  Configurable<float> TpcPidNsigmaCut{"TpcPidNsigmaCut", 5, "TpcPidNsigmaCut"};
  Configurable<float> lifetimecut{"lifetimecut", 40., "lifetimecut"}; // ct
  Configurable<float> minProtonPt{"minProtonPt", 0.3, "minProtonPt"};
  Configurable<float> maxProtonPt{"maxProtonPt", 5, "maxProtonPt"};
  Configurable<float> minPionPt{"minPionPt", 0.1, "minPionPt"};
  Configurable<float> maxPionPt{"maxPionPt", 1.2, "maxPionPt"};
  Configurable<float> minDeuteronPt{"minDeuteronPt", 0.6, "minDeuteronPt"};
  Configurable<float> maxDeuteronPt{"maxDeuteronPt", 10, "maxDeuteronPt"};
  Configurable<float> minDeuteronPUseTOF{"minDeuteronPUseTOF", 1, "minDeuteronPt Enable TOF PID"};
  Configurable<float> h3LMassLowerlimit{"h3LMassLowerlimit", 2.96, "Hypertriton mass lower limit"};
  Configurable<float> h3LMassUpperlimit{"h3LMassUpperlimit", 3.04, "Hypertriton mass upper limit"};
  Configurable<int> mincrossedrowsproton{"mincrossedrowsproton", 90, "min tpc crossed rows for pion"};
  Configurable<int> mincrossedrowspion{"mincrossedrowspion", 70, "min tpc crossed rows"};
  Configurable<int> mincrossedrowsdeuteron{"mincrossedrowsdeuteron", 100, "min tpc crossed rows for deuteron"};

  HistogramRegistry qaHists{"qaHists", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  void init(o2::framework::InitContext&)
  {
    std::vector<double> ptBinning = {0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.8, 2.0, 2.2, 2.4, 2.8, 3.2, 3.6, 4., 5.};

    AxisSpec ptAxis = {ptBinning, "#it{p}_{T} (GeV/#it{c})"};

    qaHists.add("fCollZpos", "collision z position", HistType::kTH1F, {{600, -20., +20., "z position (cm)"}});
    qaHists.add("fTPCsignal", "Specific energy loss", HistType::kTH2F, {{1200, -6, 6, "#it{p} (GeV/#it{c})"}, {1400, 0, 1400, "d#it{E} / d#it{X} (a. u.)"}});
    qaHists.add("fDeuTOFNsigma", "Deuteron TOF Nsigma distribution", HistType::kTH2F, {{1200, -6, 6, "#it{p} (GeV/#it{c})"}, {2000, -100, 100, "TOF n#sigma"}});
    qaHists.add("fH3LMassVsPt", "Hypertrion mass Vs pT", HistType::kTH2F, {{100, 0, 10, "#it{p}_{T} (GeV/#it{c})"}, {80, 2.96, 3.04, "Inv. Mass (GeV/c^{2})"}});

    for (int iN{0}; iN < nNuclei; ++iN) {
      h2TPCsignal[iN] = qaHists.add<TH2>(Form("fTPCsignal_%s", nucleiNames[iN].data()), "Specific energy loss", HistType::kTH2F, {{1200, -6, 6., "#it{p}/Z (GeV/#it{c})"}, {1400, 0, 1400, "d#it{E} / d#it{X} (a. u.)"}});
      h2TPCnSigma[iN] = qaHists.add<TH2>(Form("fTPCcounts_%s", nucleiNames[iN].data()), "n-sigma TPC", HistType::kTH2F, {{100, -5, 5, "#it{p} /Z (GeV/#it{c})"}, {200, -10., +10., "n#sigma_{He} (a. u.)"}});
    }

    auto scalers{std::get<std::shared_ptr<TH1>>(qaHists.add("fProcessedEvents", ";;Number of filtered events", HistType::kTH1F, {{nNuclei + nHyperNuclei + 1, -0.5, nNuclei + nHyperNuclei + 0.5}}))};
    scalers->GetXaxis()->SetBinLabel(1, "Processed events");
    for (uint32_t iS{0}; iS < nucleiNames.size(); ++iS) {
      scalers->GetXaxis()->SetBinLabel(iS + 2, nucleiNames[iS].data());
    }
    for (uint32_t iS{0}; iS < hypernucleiNames.size(); ++iS) {
      scalers->GetXaxis()->SetBinLabel(iS + nucleiNames.size() + 2, hypernucleiNames[iS].data());
    }
  }

  // Filter trackFilter = (nabs(aod::track::eta) < cfgCutEta);
  // using TrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA, aod::pidTPCFullPi, aod::pidTPCFullPr, aod::pidTPCFullDe, aod::pidTPCFullTr, aod::pidTPCFullHe, aod::pidTPCFullAl, aod::pidTOFFullDe, aod::pidTOFFullTr, aod::pidTOFFullHe, aod::pidTOFFullAl>>;
  using TrackCandidates = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA, aod::pidTPCFullPi, aod::pidTPCFullPr, aod::pidTPCFullDe, aod::pidTPCFullTr, aod::pidTPCFullHe, aod::pidTPCFullAl, aod::pidTOFFullDe, aod::pidTOFFullTr, aod::pidTOFFullHe, aod::pidTOFFullAl>;
  void process(aod::Collisions::iterator const& collision, aod::Vtx3BodyDatas const& vtx3bodydatas, TrackCandidates const& tracks)
  {
    // collision process loop
    bool keepEvent[nNuclei + nHyperNuclei]{false};
    //
    qaHists.fill(HIST("fCollZpos"), collision.posZ());
    qaHists.fill(HIST("fProcessedEvents"), 0);
    //
    const double bgScalings[nNuclei][2]{
      {charges[0] * cfgMomentumScalingBetheBloch->get(0u, 0u) / masses[0], charges[0] * cfgMomentumScalingBetheBloch->get(0u, 1u) / masses[0]},
      {charges[1] * cfgMomentumScalingBetheBloch->get(1u, 0u) / masses[1], charges[1] * cfgMomentumScalingBetheBloch->get(1u, 1u) / masses[1]},
      {charges[2] * cfgMomentumScalingBetheBloch->get(2u, 0u) / masses[2], charges[2] * cfgMomentumScalingBetheBloch->get(2u, 1u) / masses[2]}};

    for (auto& track : tracks) { // start loop over tracks
      if (track.itsNCls() < cfgCutNclusITS ||
          track.tpcNClsFound() < cfgCutNclusTPC) {
        continue;
      }

      if (std::abs(track.tpcNSigmaDe()) < 5) {
        qaHists.fill(HIST("fDeuTOFNsigma"), track.p() * track.sign(), track.tofNSigmaDe());
      }

      if (track.sign() > 0 && (std::abs(track.dcaXY()) > cfgCutDCAxy ||
                               std::abs(track.dcaZ()) > cfgCutDCAz)) {
        continue;
      }

      float nSigmaTPC[nNuclei]{
        track.tpcNSigmaDe(), track.tpcNSigmaTr(), track.tpcNSigmaHe()};
      const float nSigmaTOF[nNuclei]{
        track.tofNSigmaDe(), track.tofNSigmaTr(), track.tofNSigmaHe()};
      const int iC{track.sign() < 0};

      for (int iN{0}; iN < nNuclei; ++iN) {
        /// Cheap checks first
        if (track.tpcInnerParam() < cfgMinTPCmom->get(iN, iC)) {
          continue;
        }

        if (cfgBetheBlochParams->get(iN, 5u) > 0.f) {
          double expBethe{tpc::BetheBlochAleph(static_cast<double>(track.tpcInnerParam() * bgScalings[iN][iC]), cfgBetheBlochParams->get(iN, 0u), cfgBetheBlochParams->get(iN, 1u), cfgBetheBlochParams->get(iN, 2u), cfgBetheBlochParams->get(iN, 3u), cfgBetheBlochParams->get(iN, 4u))};
          double expSigma{expBethe * cfgBetheBlochParams->get(iN, 5u)};
          nSigmaTPC[iN] = static_cast<float>((track.tpcSignal() - expBethe) / expSigma);
        }
        h2TPCnSigma[iN]->Fill(track.sign() * track.tpcInnerParam(), nSigmaTPC[iN]);
        if (nSigmaTPC[iN] < cfgCutsPID->get(iN, 0u) || nSigmaTPC[iN] > cfgCutsPID->get(iN, 1u)) {
          continue;
        }
        if (track.pt() > cfgCutsPID->get(iN, 4u) && (nSigmaTOF[iN] < cfgCutsPID->get(iN, 2u) || nSigmaTOF[iN] > cfgCutsPID->get(iN, 3u))) {
          continue;
        }
        keepEvent[iN] = true;
        if (keepEvent[iN]) {
          h2TPCsignal[iN]->Fill(track.sign() * track.tpcInnerParam(), track.tpcSignal());
        }
      }

      //
      // fill QA histograms
      //
      qaHists.fill(HIST("fTPCsignal"), track.sign() * track.tpcInnerParam(), track.tpcSignal());

    } // end loop over tracks

    // hypertriton 3body loop
    for (auto& vtx : vtx3bodydatas) {

      auto track0 = vtx.track0_as<TrackCandidates>();
      auto track1 = vtx.track1_as<TrackCandidates>();
      auto track2 = vtx.track2_as<TrackCandidates>();

      if (vtx.vtxcosPA(collision.posX(), collision.posY(), collision.posZ()) < minCosPA3body) {
        continue;
      }
      float ct = vtx.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassHyperTriton;
      if (ct > lifetimecut) {
        continue;
      }
      if (vtx.dcaVtxdaughters() > dcavtxdau) {
        continue;
      }
      if ((track2.tofNSigmaDe() < TofPidNsigmaMin || track2.tofNSigmaDe() > TofPidNsigmaMax) && track2.p() > minDeuteronPUseTOF) {
        continue;
      }
      if (TMath::Abs(track0.tpcNSigmaPr()) < TpcPidNsigmaCut && TMath::Abs(track1.tpcNSigmaPi()) < TpcPidNsigmaCut && TMath::Abs(track2.tpcNSigmaDe()) < TpcPidNsigmaCut && vtx.mHypertriton() > h3LMassLowerlimit && vtx.mHypertriton() < h3LMassUpperlimit) {
        if (track0.tpcNClsCrossedRows() > mincrossedrowsproton && track1.tpcNClsCrossedRows() > mincrossedrowspion && track2.tpcNClsCrossedRows() > mincrossedrowsdeuteron) {
          if (TMath::Abs(vtx.dcatrack1topv()) > dcapiontopv) {
            keepEvent[3] = true;
            qaHists.fill(HIST("fH3LMassVsPt"), vtx.pt(), vtx.mHypertriton());
          }
        }
      }
      if (TMath::Abs(track0.tpcNSigmaPi()) < TpcPidNsigmaCut && TMath::Abs(track1.tpcNSigmaPr()) < TpcPidNsigmaCut && TMath::Abs(track2.tpcNSigmaDe()) < TpcPidNsigmaCut && vtx.mAntiHypertriton() > h3LMassLowerlimit && vtx.mAntiHypertriton() < h3LMassUpperlimit) {
        if (track0.tpcNClsCrossedRows() > mincrossedrowspion && track1.tpcNClsCrossedRows() > mincrossedrowsproton && track2.tpcNClsCrossedRows() > mincrossedrowsdeuteron) {
          if (TMath::Abs(vtx.dcatrack0topv()) > dcapiontopv) {
            keepEvent[3] = true;
            qaHists.fill(HIST("fH3LMassVsPt"), vtx.pt(), vtx.mAntiHypertriton());
          }
        }
      }
    } // end loop over hypertriton 3body decay candidates

    for (int iDecision{0}; iDecision < nNuclei + nHyperNuclei; ++iDecision) {
      if (keepEvent[iDecision]) {
        qaHists.fill(HIST("fProcessedEvents"), iDecision + 1);
      }
    }
    tags(keepEvent[2], keepEvent[3]);
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfg)
{
  return WorkflowSpec{
    adaptAnalysisTask<nucleiFilter>(cfg)};
}
