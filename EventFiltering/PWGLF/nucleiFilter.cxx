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

#include "DataFormatsTPC/BetheBlochAleph.h"
#include "ReconstructionDataFormats/Track.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Common/DataModel/EventSelection.h"
#include "../filterTables.h"

#include "Framework/HistogramRegistry.h"

#include <cmath>
#include <string>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace
{

static constexpr int nNuclei{3};
static constexpr int nCutsPID{5};
static constexpr std::array<float, nNuclei> masses{
  constants::physics::MassDeuteron, constants::physics::MassTriton,
  constants::physics::MassHelium3};
static constexpr std::array<int, nNuclei> charges{1, 1, 2};
static const std::vector<std::string> matterOrNot{"Matter", "Antimatter"};
static const std::vector<std::string> nucleiNames{"H2", "H3", "Helium"};
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
static const std::vector<std::string> betheBlochParNames{"p0", "p1", "p2", "p3", "p4", "resolution"};
} // namespace

struct nucleiFilter {

  Produces<aod::NucleiFilters> tags;

  Configurable<float> yBeam{"yBeam", 0., "Beam rapidity"};

  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};
  Configurable<float> cfgCutEta{"cfgCutEta", 1.f, "Eta range for tracks"};

  Configurable<LabeledArray<double>> cfgBetheBlochParams{"cfgBetheBlochParams", {betheBlochDefault[0], nNuclei, 6, nucleiNames, betheBlochParNames}, "TPC Bethe-Bloch parameterisation for light nuclei"};
  Configurable<LabeledArray<double>> cfgMomentumScalingBetheBloch{"cfgMomentumScalingBetheBloch", {bbMomScalingDefault[0], nNuclei, 2, nucleiNames, matterOrNot}, "TPC Bethe-Bloch momentum scaling for light nuclei"};

  Configurable<LabeledArray<float>> cfgCutsPID{"nucleiCutsPID", {cutsPID[0], nNuclei, nCutsPID, nucleiNames, cutsNames}, "Nuclei PID selections"};

  HistogramRegistry spectra{"spectra", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  void init(o2::framework::InitContext&)
  {
    std::vector<double> ptBinning = {0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.8, 2.0, 2.2, 2.4, 2.8, 3.2, 3.6, 4., 5.};
    std::vector<double> centBinning = {0., 1., 5., 10., 20., 30., 40., 50., 70., 100.};

    AxisSpec ptAxis = {ptBinning, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec centAxis = {centBinning, "V0M (%)"};

    spectra.add("fCollZpos", "collision z position", HistType::kTH1F, {{600, -20., +20., "z position (cm)"}});
    spectra.add("fTPCsignal", "Specific energy loss", HistType::kTH2F, {{600, 0., 3, "#it{p} (GeV/#it{c})"}, {1400, 0, 1400, "d#it{E} / d#it{X} (a. u.)"}});
    spectra.add("fTPCcounts", "n-sigma TPC", HistType::kTH2F, {ptAxis, {200, -100., +100., "n#sigma_{He} (a. u.)"}});

    auto scalers{std::get<std::shared_ptr<TH1>>(spectra.add("fProcessedEvents", ";;Number of filtered events", HistType::kTH1F, {{4, -0.5, 3.5}}))};
    for (uint32_t iS{1}; iS <= nucleiNames.size(); ++iS) {
      scalers->GetXaxis()->SetBinLabel(iS, nucleiNames[iS - 1].data());
    }
  }

  Filter trackFilter = (nabs(aod::track::eta) < cfgCutEta) && (requireGlobalTrackInFilter());

  using TrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA, aod::pidTPCFullDe, aod::pidTPCFullTr, aod::pidTPCFullHe, aod::pidTPCFullAl, aod::pidTOFFullDe, aod::pidTOFFullTr, aod::pidTOFFullHe, aod::pidTOFFullAl>>;
  void process(aod::Collisions::iterator const& collision, TrackCandidates const& tracks)
  {
    // collision process loop
    bool keepEvent[nNuclei]{false};
    //
    spectra.fill(HIST("fCollZpos"), collision.posZ());
    //
    const double bgScalings[nNuclei][2]{
      {charges[0] * cfgMomentumScalingBetheBloch->get(0u, 0u) / masses[0], charges[0] * cfgMomentumScalingBetheBloch->get(0u, 1u) / masses[0]},
      {charges[1] * cfgMomentumScalingBetheBloch->get(1u, 0u) / masses[1], charges[1] * cfgMomentumScalingBetheBloch->get(1u, 1u) / masses[1]},
      {charges[2] * cfgMomentumScalingBetheBloch->get(2u, 0u) / masses[2], charges[2] * cfgMomentumScalingBetheBloch->get(2u, 1u) / masses[2]}};

    for (auto& track : tracks) { // start loop over tracks

      float nSigmaTPC[nNuclei]{
        track.tpcNSigmaDe(), track.tpcNSigmaTr(), track.tpcNSigmaHe()};
      const float nSigmaTOF[nNuclei]{
        track.tofNSigmaDe(), track.tofNSigmaTr(), track.tofNSigmaHe()};
      const int iC{track.sign() < 0};

      for (int iN{0}; iN < nNuclei; ++iN) {
        if (cfgBetheBlochParams->get(iN, 5u) > 0.f) {
          double expBethe{tpc::BetheBlochAleph(static_cast<double>(track.tpcInnerParam() * bgScalings[iN][iC]), cfgBetheBlochParams->get(iN, 0u), cfgBetheBlochParams->get(iN, 1u), cfgBetheBlochParams->get(iN, 2u), cfgBetheBlochParams->get(iN, 3u), cfgBetheBlochParams->get(iN, 4u))};
          double expSigma{expBethe * cfgBetheBlochParams->get(iN, 5u)};
          nSigmaTPC[iN] = static_cast<float>((track.tpcSignal() - expBethe) / expSigma);
        }
        if (nSigmaTPC[iN] < cfgCutsPID->get(iN, 0u) || nSigmaTPC[iN] > cfgCutsPID->get(iN, 1u)) {
          continue;
        }
        if (track.pt() > cfgCutsPID->get(iN, 4u) && (nSigmaTOF[iN] < cfgCutsPID->get(iN, 2u) || nSigmaTOF[iN] > cfgCutsPID->get(iN, 3u))) {
          continue;
        }
        keepEvent[iN] = true;
      }

      //
      // fill QA histograms
      //
      spectra.fill(HIST("fTPCsignal"), track.tpcInnerParam(), track.tpcSignal());
      spectra.fill(HIST("fTPCcounts"), track.tpcInnerParam(), nSigmaTPC[2]);

    } // end loop over tracks
    //
    for (int iDecision{0}; iDecision < 4; ++iDecision) {
      if (keepEvent[iDecision]) {
        spectra.fill(HIST("fProcessedEvents"), iDecision);
      }
    }
    tags(keepEvent[0], keepEvent[1], keepEvent[2]);
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfg)
{
  return WorkflowSpec{
    adaptAnalysisTask<nucleiFilter>(cfg)};
}
